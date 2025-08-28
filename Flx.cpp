import sys, pathlib, numpy as np, torch, trimesh

# --- OpenVDB import ---
try:
    import openvdb as vdb
except ImportError:
    import pyopenvdb as vdb  # older wheels

# --- FlexiCubes via Kaolin ---
from kaolin.non_commercial.flexicubes import ops as fc_ops

def main(vdb_path, grid_name, iso=0.0, target_res=256, out_dir="out"):
    out = pathlib.Path(out_dir); out.mkdir(exist_ok=True)

    # 1) Load VDB grid
    grid = vdb.read(vdb_path, grid_name)
    imin, imax = grid.evalActiveVoxelBoundingBox()
    pad = 2
    imin = (imin[0]-pad, imin[1]-pad, imin[2]-pad)
    imax = (imax[0]+pad, imax[1]+pad, imax[2]+pad)

    # 2) Decide dense sampling shape
    dims = (imax[0]-imin[0]+1, imax[1]-imin[1]+1, imax[2]-imin[2]+1)
    scale = max(dims) / float(target_res)
    dense_shape = tuple(int(np.ceil(d/scale)) for d in dims)

    # 3) Pull dense SDF from VDB into NumPy
    sdf_np = np.empty(dense_shape, dtype=np.float32)
    grid.copyToArray(sdf_np, ijk=imin)  # index-space sampling

    # 4) -> torch tensor [1, D, H, W]
    sdf_t = torch.from_numpy(sdf_np).float()
    if torch.cuda.is_available():
        sdf_t = sdf_t.cuda()
    sdf_t = sdf_t.unsqueeze(0)

    # 5) FlexiCubes extraction
    verts_t, faces_t = fc_ops.extract_mesh_from_sdf(sdf_t, iso)
    V = verts_t.detach().cpu().numpy()
    F = faces_t.detach().cpu().numpy()

    # 6) Map index-space -> world-space using VDB transform
    xfm = grid.transform
    def ijk_to_world(p):
        return np.asarray(xfm.indexToWorld((float(p[0]), float(p[1]), float(p[2]))), dtype=np.float32)
    V_world = np.vstack([ijk_to_world(p) for p in V])

    # 7) Save
    mesh = trimesh.Trimesh(V_world, F, process=False)
    obj_path = out / "flexicubes_mesh.obj"
    mesh.export(obj_path)
    print(f"✓ Wrote {obj_path}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python vdb_to_flexicubes.py <path.vdb> <grid_name> [iso] [target_res]")
        sys.exit(1)
    vdb_path = sys.argv[1]
    grid_name = sys.argv[2]
    iso = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
    target_res = int(sys.argv[4]) if len(sys.argv) > 4 else 256
    main(vdb_path, grid_name, iso, target_rep1x+
