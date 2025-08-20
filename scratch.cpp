#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>   // GridSampler
#include <manifold/manifold.h>              // Manifold::LevelSet

// Assume: grid is a FloatGrid::Ptr that stores a *level set* (SDF).
manifold::Manifold meshFromVDB(openvdb::FloatGrid::ConstPtr grid)
{
    // 1) Continuous sampler over the VDB in *world* space.
    openvdb::tools::GridSampler<
        openvdb::FloatGrid, openvdb::tools::BoxSampler> sampler(*grid);

    // 2) SDF callback for Manifold (expects signed distance; zero isosurface).
    // OpenVDB convention: negative inside, positive outside — also what Manifold expects.
    auto sdf = [&](manifold::vec3 p) -> double {
        return sampler.wsSample(openvdb::Vec3d(p.x, p.y, p.z));
    };

    // 3) Bounds (world-space) from the grid’s active voxel bbox, padded a bit.
    openvdb::CoordBBox ibox;
    grid->tree().evalActiveVoxelBoundingBox(ibox);
    const auto& xform = grid->transform();
    openvdb::Vec3d wmin = xform.indexToWorld(ibox.min().asVec3d());
    openvdb::Vec3d wmax = xform.indexToWorld(ibox.max().asVec3d());

    // pad by a few voxels to ensure the narrow band is fully inside
    const openvdb::Vec3d vox = xform.voxelSize();
    const openvdb::Vec3d pad = 3.0 * vox; // adjust if your band is wider
    wmin -= pad; wmax += pad;

    manifold::Box bounds{
        manifold::vec3(wmin.x(), wmin.y(), wmin.z()),
        manifold::vec3(wmax.x(), wmax.y(), wmax.z())
    };

    // 4) Grid step for marching tetrahedra. Start near voxel size.
    const double edgeLength = (vox.x() + vox.y() + vox.z()) / 3.0;

    // 5) Extract zero level (level=0). Leave tolerance default; enable parallel.
    return manifold::Manifold::LevelSet(sdf, bounds, edgeLength, /*level=*/0.0,
                                        /*tolerance=*/-1.0, /*canParallel=*/true);
}

