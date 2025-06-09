##include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Polygon_mesh_processing/corefine_and_compute_difference.h>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_exact_constructions_kernel EPECK;
typedef CGAL::Surface_mesh<EPECK::Point_3> ExactMesh;

// Your original mesh:
typedef CGAL::Filtered_kernel<CartesianKernel> FilteredKernel;
typedef CGAL::Surface_mesh<FilteredKernel::Point_3> FilteredMesh;

FilteredMesh inputA, inputB; // your input meshes
ExactMesh exactA, exactB;

CGAL::copy_face_graph(inputA, exactA);
CGAL::copy_face_graph(inputB, exactB);

// Boolean
ExactMesh result;
bool ok = PMP::corefine_and_compute_difference(exactA, exactB, result);

// (Optional) Copy result back
FilteredMesh finalResult;
CGAL::copy_face_graph(result, finalResult);


