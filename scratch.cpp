
#include <iostream>
#include <cmath>

// Helper struct for 3D vectors
struct Vec3 {
    double x, y, z;

    // Normalize the vector
    void normalize() {
        double len = std::sqrt(x * x + y * y + z * z);
        if (len > 1e-6) {
            x /= len;
            y /= len;
            z /= len;
        }
    }

    // Length in XY-plane
    double lengthXY() const {
        return std::sqrt(x * x + y * y);
    }
};

// Interpolation function in 3D with only 4 coefficients (XY projection)
double interpolateCoefficientXY(const Vec3& normal, double coeffNorth, double coeffSouth, double coeffEast, double coeffWest) {
    Vec3 n = normal;
    n.normalize();

    // Project onto XY-plane
    Vec3 projected = { n.x, n.y, 0 };
    double lenXY = projected.lengthXY();

    // If XY projection is too small, default to average
    if (lenXY < 1e-6) {
        return (coeffNorth + coeffSouth + coeffEast + coeffWest) * 0.25;
    }

    // Normalize projected vector
    projected.x /= lenXY;
    projected.y /= lenXY;

    // Define axis directions
    Vec3 north = { 0, 1, 0 };
    Vec3 south = { 0, -1, 0 };
    Vec3 east  = { 1, 0, 0 };
    Vec3 west  = { -1, 0, 0 };

    // Compute positive dot products (weights)
    double weightN = std::max(0.0, projected.y);
    double weightS = std::max(0.0, -projected.y);
    double weightE = std::max(0.0, projected.x);
    double weightW = std::max(0.0, -projected.x);

    double weightSum = weightN + weightS + weightE + weightW;

    if (weightSum < 1e-6)
        return 0.0; // Degenerate case

    // Interpolate
    double interpolated = (coeffNorth * weightN + coeffSouth * weightS + coeffEast * weightE + coeffWest * weightW) / weightSum;

    return interpolated;
}

int main() {
    // Example coefficients
    double coeffN = 1.0;
    double coeffS = 0.5;
    double coeffE = 0.8;
    double coeffW = 0.3;

    // Example 3D direction vector
    Vec3 direction = { 1.0, 1.0, 1.0 }; // 45 degrees NE, tilted up

    double result = interpolateCoefficientXY(direction, coeffN, coeffS, coeffE, coeffW);
    std::cout << "Interpolated Coefficient: " << result << std::endl;

    return 0;
}
