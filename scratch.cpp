#include <iostream>
#include <cmath>

// Normalize angle to [0, 2π)
double normalizeAngle(double angle) {
    while (angle < 0) angle += 2.0 * M_PI;
    while (angle >= 2.0 * M_PI) angle -= 2.0 * M_PI;
    return angle;
}

// Linear interpolation helper
double lerp(double a, double b, double t) {
    return a * (1.0 - t) + b * t;
}

// Spherical coordinate interpolation with 4 coefficients
double sphericalInterpolate(
    double coeffNorth, double coeffSouth, double coeffEast, double coeffWest,
    double nx, double ny, double nz)
{
    // Normalize vector
    double len = std::sqrt(nx * nx + ny * ny + nz * nz);
    if (len < 1e-6) return (coeffNorth + coeffSouth + coeffEast + coeffWest) * 0.25;
    nx /= len;
    ny /= len;
    nz /= len;

    // Spherical coordinates
    double theta = atan2(ny, nx); // Azimuth: [-π, π]
    theta = normalizeAngle(theta); // [0, 2π)
    double phi = std::asin(nz);   // Elevation: [-π/2, π/2]

    // Map elevation phi to a vertical blend factor:
    // phi = π/2 (up) -> 1 (North)
    // phi = -π/2 (down) -> 1 (South)
    // phi = 0 (flat) -> 0
    double verticalBlendNorth = (nz > 0) ? nz : 0.0;  // Only up contributes to North
    double verticalBlendSouth = (nz < 0) ? -nz : 0.0; // Only down contributes to South
    double horizontalWeight = 1.0 - (verticalBlendNorth + verticalBlendSouth);

    // Determine azimuthal sector (East/North/West/South)
    const double sectorSize = M_PI / 2.0;
    double sectorPos = theta / sectorSize;
    int sector = static_cast<int>(sectorPos);
    double t = sectorPos - sector;

    double coeffA, coeffB;

    switch (sector) {
        case 0: coeffA = coeffEast;  coeffB = coeffNorth; break; // East to North
        case 1: coeffA = coeffNorth; coeffB = coeffWest;  break; // North to West
        case 2: coeffA = coeffWest;  coeffB = coeffSouth; break; // West to South
        case 3: coeffA = coeffSouth; coeffB = coeffEast;  break; // South to East
        default: coeffA = coeffB = 0.0; break;
    }

    // Interpolate in azimuth (longitude)
    double azimuthInterp = lerp(coeffA, coeffB, t);

    // Final blend: vertical North/South influence + horizontal azimuth interpolation
    double finalCoeff =
        verticalBlendNorth * coeffNorth +
        verticalBlendSouth * coeffSouth +
        horizontalWeight * azimuthInterp;

    return finalCoeff;
}

int main() {
    // Example coefficients
    double coeffNorth = 1.0;
    double coeffSouth = 0.5;
    double coeffEast  = 0.8;
    double coeffWest  = 0.3;

    // Example normal vectors:
    std::cout << "Up (0,0,1): " << sphericalInterpolate(coeffNorth, coeffSouth, coeffEast, coeffWest, 0, 0, 1) << std::endl;
    std::cout << "Down (0,0,-1): " << sphericalInterpolate(coeffNorth, coeffSouth, coeffEast, coeffWest, 0, 0, -1) << std::endl;
    std::cout << "East (1,0,0): " << sphericalInterpolate(coeffNorth, coeffSouth, coeffEast, coeffWest, 1, 0, 0) << std::endl;
    std::cout << "North-East (1,1,0): " << sphericalInterpolate(coeffNorth, coeffSouth, coeffEast, coeffWest, 1, 1, 0) << std::endl;
    std::cout << "North-East-Up (1,1,1): " << sphericalInterpolate(coeffNorth, coeffSouth, coeffEast, coeffWest, 1, 1, 1) << std::endl;

    return 0;
}
