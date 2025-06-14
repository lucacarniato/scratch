#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/minkowski_sum_3.h>
#include <CGAL/boost/graph/convert_surface_mesh.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/make_icosahedron.h>

typedef CGAL::Simple_cartesian<double>         Kernel;
typedef Kernel::Point_3                        Point;
typedef CGAL::Surface_mesh<Point>              SurfaceMesh;
typedef CGAL::Polyhedron_3<Kernel>             Polyhedron;
typedef CGAL::Nef_polyhedron_3<Kernel>         Nef_polyhedron;

SurfaceMesh offset_surface(const SurfaceMesh& input_sm, double offset_radius) {
    // Step 1: Convert Surface_mesh → Polyhedron_3
    Polyhedron input_poly;
    CGAL::copy_face_graph(input_sm, input_poly);

    // Step 2: Create small sphere polyhedron for Minkowski sum
    Polyhedron sphere_poly;
    CGAL::make_icosahedron(sphere_poly); // Low-res sphere (~20 triangles)

    // Scale the icosahedron to desired offset radius
    for (auto it = sphere_poly.points_begin(); it != sphere_poly.points_end(); ++it) {
        *it = Point(it->x() * offset_radius, it->y() * offset_radius, it->z() * offset_radius);
    }

    // Step 3: Convert to Nef polyhedra
    Nef_polyhedron nef_input(input_poly);
    Nef_polyhedron nef_sphere(sphere_poly);

    // Step 4: Perform Minkowski sum
    Nef_polyhedron nef_result = CGAL::minkowski_sum_3(nef_input, nef_sphere);

    // Step 5: Convert result back to Polyhedron_3
    Polyhedron result_poly;
    nef_result.convert_to_polyhedron(result_poly);

    // Step 6: Convert back to Surface_mesh
    SurfaceMesh result_sm;
    CGAL::copy_face_graph(result_poly, result_sm);

    return result_sm;
}



