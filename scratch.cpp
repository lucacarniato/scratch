#include <iostream>
#include <cmath>
#include <algorithm>

// Clamp negative dot products to 0
double saturate(double x) {
    return std::max(0.0, x);
}

double sphericalInterpolate6(
    double coeffNorth, double coeffSouth,
    double coeffEast, double coeffWest,
    double nx, double ny, double nz)
{
    // Normalize the direction vector
    double len = std::sqrt(nx * nx + ny * ny + nz * nz);
    if (len < 1e-6) {
        // Degenerate vector: return average
        return (coeffNorth + coeffSouth + coeffEast + coeffWest) * 0.25;
    }

    nx /= len;
    ny /= len;
    nz /= len;

    // Dot products with axis directions
    double dotEast  = saturate( nx);     // +X
    double dotWest  = saturate(-nx);     // -X
    double dotUp    = saturate( ny);     // +Y (mapped to East)
    double dotDown  = saturate(-ny);     // -Y (mapped to West)
    double dotNorth = saturate( nz);     // +Z
    double dotSouth = saturate(-nz);     // -Z

    // Use assigned mappings for Up and Down
    dotEast += dotUp;   // Up contributes to East
    dotWest += dotDown; // Down contributes to West

    // Total weight
    double totalWeight = dotEast + dotWest + dotNorth + dotSouth;
    if (totalWeight < 1e-6) {
        return (coeffNorth + coeffSouth + coeffEast + coeffWest) * 0.25;
    }

    // Weighted sum
    double interpolated =
        (coeffEast  * dotEast +
         coeffWest  * dotWest +
         coeffNorth * dotNorth +
         coeffSouth * dotSouth) / totalWeight;

    return interpolated;
}

int main() {
    // Example coefficients
    double coeffNorth = 1.0;
    double coeffSouth = 0.3;
    double coeffEast  = 0.8;
    double coeffWest  = 0.4;

    // Example directions
    std::cout << "North (0, 0, 1):  " << sphericalInterpolate6(coeffNorth, coeffSouth, coeffEast, coeffWest, 0, 0, 1) << std::endl;
    std::cout << "South (0, 0, -1): " << sphericalInterpolate6(coeffNorth, coeffSouth, coeffEast, coeffWest, 0, 0, -1) << std::endl;
    std::cout << "East (1, 0, 0):   " << sphericalInterpolate6(coeffNorth, coeffSouth, coeffEast, coeffWest, 1, 0, 0) << std::endl;
    std::cout << "West (-1, 0, 0):  " << sphericalInterpolate6(coeffNorth, coeffSouth, coeffEast, coeffWest, -1, 0, 0) << std::endl;
    std::cout << "Up (0, 1, 0):     " << sphericalInterpolate6(coeffNorth, coeffSouth, coeffEast, coeffWest, 0, 1, 0) << std::endl;
    std::cout << "Down (0, -1, 0):  " << sphericalInterpolate6(coeffNorth, coeffSouth, coeffEast, coeffWest, 0, -1, 0) << std::endl;
    std::cout << "Diagonal (1, 1, 1): " << sphericalInterpolate6(coeffNorth, coeffSouth, coeffEast, coeffWest, 1, 1, 1) << std::endl;

    return 0;
}
