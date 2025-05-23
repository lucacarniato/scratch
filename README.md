bool isValidSDF(const openvdb::FloatGrid& grid,
                double  gradTol    = 0.25,
                int     bandWidth  = 3,
                size_t  maxSamples = 500'000)
{
    // 0. quick metadata sanity
    if (grid.getGridClass() != openvdb::GRID_LEVEL_SET) return false;

    const double voxelSize = grid.voxelSize()[0];              // assume isotropic
    const double band      = bandWidth * voxelSize;

    // 1. make sure we do have both negative and positive samples
    bool hasNeg = false, hasPos = false;
    for (auto it = grid.cbeginValueOn(); it; ++it) {
        const float v = *it;
        hasNeg |= v < 0.0f;
        hasPos |= v > 0.0f;
        if (hasNeg && hasPos) break;
    }
    if (!(hasNeg && hasPos)) return false;

    // 2. check the Eikonal property near the surface
    const auto acc = grid.getConstAccessor();
    size_t tested = 0;
    for (auto it = grid.cbeginValueOn(); it; ++it) {
        const float v = *it;
        if (std::abs(v) > band) continue;                 // only near φ = 0
        if (++tested > maxSamples) break;

        const openvdb::Coord ijk = it.getCoord();
        // central differences, world-space
        const double dx = (acc.getValue(ijk.offsetBy(1,0,0)) -
                           acc.getValue(ijk.offsetBy(-1,0,0))) / (2.0 * voxelSize);
        const double dy = (acc.getValue(ijk.offsetBy(0,1,0)) -
                           acc.getValue(ijk.offsetBy(0,-1,0))) / (2.0 * voxelSize);
        const double dz = (acc.getValue(ijk.offsetBy(0,0,1)) -
                           acc.getValue(ijk.offsetBy(0,0,-1))) / (2.0 * voxelSize);

        const double mag = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (std::abs(mag - 1.0) > gradTol) return false; // fails Eikonal
    }
    return true;
}

} // anonymous namespace

//-----------------------------------------------------------------------

bool checkSDFInFile(const std::string& vdbPath,
                    const std::string& gridName = "distance")
{
    openvdb::initialize();

    openvdb::io::File file(vdbPath);
    file.open();
    openvdb::GridBase::Ptr base = file.readGrid(gridName);
    file.close();

    if (!base) {
        std::cerr << "Grid \"" << gridName << "\" not found.\n";
        return false;
    }
    if (base->valueType() != openvdb::FloatGrid::valueType()) {
        std::cerr << "Grid is not a FloatGrid.\n";
        return false;
    }
    const auto sdfGrid = openvdb::gridConstPtrCast<openvdb::FloatGrid>(base);
    return isValidSDF(*sdfGrid);
}
