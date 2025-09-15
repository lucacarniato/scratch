# OpenVDB Python bindings (name can be openvdb or pyopenvdb depending on your build)
pip install openvdb || pip install pyopenvdb

# Core deps
pip install numpy torch trimesh

# Kaolin (includes FlexiCubes)
pip install kaolin


python vdb_to_flexicubes.py path/to/your.vdb --grid <grid_name> --res 64 64 64 --iso 0.0 --out mesh.obj



#!/usr/bin/env python3
import argparse, sys, numpy as np, pathlib, math
import torch
import trimesh

# ---- OpenVDB import (module name varies) ----
try:
    import openvdb as vdb
except ImportError:
    import pyopenvdb as vdb  # older builds

# ---- Kaolin FlexiCubes ----
from kaolin.ops.conversions import FlexiCubes  # documented API
# See: kaolin.readthedocs.io "kaolin.ops.conversions.FlexiCubes"  # :contentReference[oaicite:1]{index=1}

def list_grids(vdb_path):
    try:
        allg = vdb.readAll(vdb_path)
        return list(allg.keys())
    except Exception as e:
        print(f"Failed to read {vdb_path}: {e}")
        return []

def get_grid(vdb_path, grid_name=None):
    if grid_name is None:
        # pick the first FloatGrid if no name provided
        allg = vdb.readAll(vdb_path)
        for k, g in allg.items():
            if g.__class__.__name__ == 'FloatGrid':
                return g, k
        # fallback to first grid
        if allg:
            k = next(iter(allg))
            return allg[k], k
        raise RuntimeError("No grids found in VDB.")
    else:
        return vdb.read(vdb_path, grid_name), grid_name

def tri_linear_sample_floatgrid(grid, pts_world):
    """
    Tri-linearly sample a FloatGrid at given world-space positions.
    OpenVDB values live at voxel centers; we convert to index space and interpolate.
    """
    xfm = grid.transform
    # accessor (const if available)
    try:
        acc = grid.getConstAccessor()
    except Exception:
        acc = grid.getAccessor()

    def get(i, j, k):
        # Returns background automatically for inactive voxels
        return acc.getValue((int(i), int(j), int(k)))

    out = np.empty((pts_world.shape[0],), dtype=np.float32)
    for n in range(pts_world.shape[0]):
        xw, yw, zw = float(pts_world[n,0]), float(pts_world[n,1]), float(pts_world[n,2])
        # world -> fractional index coordinates
        xi, yi, zi = xfm.worldToIndex((xw, yw, zw))
        i0, j0, k0 = math.floor(xi), math.floor(yi), math.floor(zi)
        fx, fy, fz = xi - i0, yi - j0, zi - k0
        i1, j1, k1 = i0 + 1, j0 + 1, k0 + 1
        # fetch 8 corners
        c000 = get(i0, j0, k0); c100 = get(i1, j0, k0)
        c010 = get(i0, j1, k0); c110 = get(i1, j1, k0)
        c001 = get(i0, j0, k1); c101 = get(i1, j0, k1)
        c011 = get(i0, j1, k1); c111 = get(i1, j1, k1)
        # trilinear blend
        c00 = c000 * (1 - fx) + c100 * fx
        c01 = c001 * (1 - fx) + c101 * fx
        c10 = c010 * (1 - fx) + c110 * fx
        c11 = c011 * (1 - fx) + c111 * fx
        c0  = c00  * (1 - fy) + c10  * fy
        c1  = c01  * (1 - fy) + c11  * fy
        out[n] = c0 * (1 - fz) + c1 * fz
    return out

def main():
    ap = argparse.ArgumentParser(description="VDB -> FlexiCubes mesh")
    ap.add_argument("vdb", help="Path to .vdb file")
    ap.add_argument("--grid", help="Grid name (FloatGrid level set). If omitted, first FloatGrid is used.")
    ap.add_argument("--res", nargs="+", type=int, default=[64, 64, 64],
                    help="Voxel grid resolution as 1 or 3 ints (e.g. 128 128 128 or 128).")
    ap.add_argument("--iso", type=float, default=0.0, help="Isovalue for extraction (0 for SDF/level set).")
    ap.add_argument("--out", default="flexicubes_mesh.obj", help="Output OBJ/PLY path")
    args = ap.parse_args()

    # Normalize resolution
    if len(args.res) == 1:
        res = [args.res[0]] * 3
    elif len(args.res) == 3:
        res = args.res
    else:
        print("--res must be 1 or 3 integers"); sys.exit(2)
    res = [int(max(2, r)) for r in res]  # min 2

    print(f"Loading grid from {args.vdb}...")
    grid, used_name = get_grid(args.vdb, args.grid)
    print(f"  Using grid: {used_name} ({grid.__class__.__name__})")

    # Active voxel bbox in INDEX space (inclusive mins/max)
    (imin, jmin, kmin), (imax, jmax, kmax) = grid.evalActiveVoxelBoundingBox()
    # Map bbox corners to WORLD space
    min_world = np.array(grid.transform.indexToWorld((imin, jmin, kmin)), dtype=np.float64)
    # include top corner: +1 in index to reach the far face corner
    max_world = np.array(grid.transform.indexToWorld((imax + 1, jmax + 1, kmax + 1)), dtype=np.float64)
    extent    = max_world - min_world

    # ---- Build a voxel grid (vertices & cube indices) with Kaolin ----
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    fc = FlexiCubes(device=device)
    # Returns normalized vertices in [-0.5, 0.5]^3 and cube topology
    vox_vertices, cube_idx = fc.construct_voxel_grid(resolution=res)  # tensors on device
    # Map to WORLD space for sampling SDF:
    # verts_norm -> in [-0.5,0.5] ; map to [min_world, max_world]
    verts_world = (vox_vertices.detach().cpu().numpy() + 0.5) * extent + min_world

    print(f"Sampling VDB at {verts_world.shape[0]:,} vertices (this may take a minute at higher resolutions)...")
    sdf_vals = tri_linear_sample_floatgrid(grid, verts_world).astype(np.float32)

    # Move inputs to device for FlexiCubes
    voxelgrid_vertices = vox_vertices.to(device=device)
    scalar_field = torch.from_numpy(sdf_vals).to(device=device)
    cube_idx = cube_idx.to(device=device)

    print("Running FlexiCubes...")
    # FlexiCubes call; we pass resolution and leave advanced params as defaults.
    V_fc, F_fc, *_ = fc(voxelgrid_vertices=voxelgrid_vertices,
                        scalar_field=scalar_field,
                        cube_idx=cube_idx,
                        resolution=res,
                        training=False)

    # FlexiCubes vertices are in the same normalized cube as construct_voxel_grid returned:
    V_world = (V_fc.detach().cpu().numpy() + 0.5) * extent + min_world
    F = F_fc.detach().cpu().numpy()

    print(f"Writing mesh: {args.out}")
    pathlib.Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    mesh = trimesh.Trimesh(V_world, F, process=False)
    mesh.export(args.out)
    print("Done.")

if __name__ == "__main__":
    main()
