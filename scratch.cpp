#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/boost/graph/copy_face_graph.h>          // copy_face_graph()
#include <CGAL/Polygon_mesh_processing/repair.h>        // optional clean-up
namespace PMP = CGAL::Polygon_mesh_processing;

using Kernel       = CGAL::Exact_predicates_exact_constructions_kernel;
using Point        = Kernel::Point_3;
using Surface_mesh = CGAL::Surface_mesh<Point>;
using Polyhedron   = CGAL::Polyhedron_3<Kernel>;
using Nef          = CGAL::Nef_polyhedron_3<Kernel>;

Polyhedron P, Q;
CGAL::copy_face_graph(smA, P);      // Surface_mesh  -> Polyhedron
CGAL::copy_face_graph(smB, Q);


Nef nA(P);            // ctor from Polyhedron_3
Nef nB(Q);

Nef nDiff = nA - nB;  // Boolean: A \ B


Polyhedron Pdiff;
if ( !nDiff.convert_to_polyhedron(Pdiff) )
{
  std::cerr << "Boolean result is not a 2-manifold – cannot export\n";
  return EXIT_FAILURE;
}

Surface_mesh diff;
CGAL::copy_face_graph(Pdiff, diff);    // Polyhedron -> Surface_mesh


