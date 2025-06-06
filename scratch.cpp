#include "manifold.h"
#include <glm/glm.hpp>
#include <vector>
#include <iostream>

using namespace manifold;

// Your mesh data: fill these with your actual geometry
std::vector<glm::vec3> verticesA;
std::vector<glm::ivec3> trianglesA;
std::vector<glm::vec3> verticesB;
std::vector<glm::ivec3> trianglesB;

Manifold BooleanDifferenceAndClean(
    const std::vector<glm::vec3>& vertsA,
    const std::vector<glm::ivec3>& trisA,
    const std::vector<glm::vec3>& vertsB,
    const std::vector<glm::ivec3>& trisB,
    float relativeThreshold = 0.01f  // 1% of largest volume
) {
    // Step 1: Create mesh structs
    Mesh meshA, meshB;
    meshA.vertPos = vertsA;
    meshA.triVerts = trisA;

    meshB.vertPos = vertsB;
    meshB.triVerts = trisB;

    // Step 2: Construct Manifolds
    Manifold A(meshA);
    Manifold B(meshB);

    // Step 3: Boolean Difference
    Manifold result = A - B;

    // Step 4: Split into components
    std::vector<Manifold> components = result.Split();

    // Step 5: Find largest volume
    float maxVolume = 0.0f;
    for (const auto& c : components) {
        float vol = c.GetProperties().volume;
        if (vol > maxVolume) {
            maxVolume = vol;
        }
    }

    // Step 6: Set threshold based on largest volume
    float volumeThreshold = relativeThreshold * maxVolume;

    // Step 7: Filter components
    std::vector<Manifold> filtered;
    for (const auto& c : components) {
        float vol = c.GetProperties().volume;
        if (vol >= volumeThreshold) {
            filtered.push_back(c);
        }
    }

    // Step 8: Recombine components
    if (filtered.empty()) {
        std::cerr << "Warning: All components filtered out. Returning empty manifold.\n";
        return Manifold();
    }

    Manifold cleaned = filtered[0];
    for (size_t i = 1; i < filtered.size(); ++i) {
        cleaned = cleaned + filtered[i];
    }

    return cleaned;
}

int main() {
    // TODO: Fill verticesA, trianglesA, verticesB, trianglesB with your mesh data.

    Manifold cleaned = BooleanDifferenceAndClean(verticesA, trianglesA, verticesB, trianglesB, 0.01f);

    // Export result
    cleaned.ExportMesh("cleaned_difference.obj");

    std::cout << "Cleaned mesh exported.\n";
    return 0;
}

