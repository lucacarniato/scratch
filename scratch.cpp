#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/OBJ_reader.h>

#include "manifold.h"
#include <unordered_map>
#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Surface_mesh<Point> SurfaceMesh;

// Convert CGAL Surface_mesh to manifold::Mesh
manifold::Mesh ConvertToManifold(const SurfaceMesh& mesh) {
    using namespace manifold;
    std::vector<glm::vec3> vertices;
    std::vector<glm::ivec3> triangles;
    std::unordered_map<SurfaceMesh::Vertex_index, int> vmap;

    int i = 0;
    for (auto v : mesh.vertices()) {
        Point p = mesh.point(v);
        vertices.emplace_back(p.x(), p.y(), p.z());
        vmap[v] = i++;
    }

    for (auto f : mesh.faces()) {
        std::vector<SurfaceMesh::Vertex_index> vs;
        for (auto v : CGAL::vertices_around_face(mesh.halfedge(f), mesh))
            vs.push_back(v);

        if (vs.size() != 3) continue; // skip non-triangular
        triangles.emplace_back(vmap[vs[0]], vmap[vs[1]], vmap[vs[2]]);
    }

    return {vertices, triangles};
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " mesh1.obj mesh2.obj result.obj\n";
        return 1;
    }

    SurfaceMesh mesh1, mesh2;
    std::ifstream in1(argv[1]), in2(argv[2]);
    if (!CGAL::IO::read_OBJ(in1, mesh1) || !CGAL::IO::read_OBJ(in2, mesh2)) {
        std::cerr << "Error reading OBJ files.\n";
        return 1;
    }

    manifold::Mesh m1 = ConvertToManifold(mesh1);
    manifold::Mesh m2 = ConvertToManifold(mesh2);

    manifold::Manifold solidA(m1);
    manifold::Manifold solidB(m2);

    manifold::Manifold result = solidA - solidB;
    manifold::Mesh output = result.GetMesh();

    manifold::ExportMesh(argv[3], output);

    std::cout << "Boolean difference complete. Output written to " << argv[3] << "\n";
    return 0;
}
