// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>

// Manifold
#include "manifold/meshGL.h"
#include "manifold/manifold.h"

#include <unordered_map>
#include <vector>
#include <array>
#include <cmath>

namespace PMP = CGAL::Polygon_mesh_processing;
using K   = CGAL::Simple_cartesian<double>;
using SM  = CGAL::Surface_mesh<K::Point_3>;
using VI  = SM::Vertex_index;
using FI  = SM::Face_index;
using EI  = SM::Edge_index;
using HI  = SM::Halfedge_index;

// --- Step 1: read/prepare your DC mesh into `SM mesh` (omitted) ---

// Optional: basic repairs / checks
void repair_mesh(SM& mesh) {
  PMP::remove_isolated_vertices(mesh);
  PMP::stitch_borders(mesh);
  // You can also call PMP::merge_duplicate_vertices, etc., as needed.
}

// Step 2: detect sharp feature edges in CGAL
std::unordered_set<EI> detect_sharp_edges(const SM& mesh, double sharp_angle_deg) {
  std::unordered_set<EI> sharp;
  PMP::detect_sharp_edges(mesh, CGAL::parameters::angle_in_degrees(sharp_angle_deg)
                                      .face_patch_map(nullptr) /* default */,
                          std::inserter(sharp, sharp.end()));
  // (Border edges aren’t expected on a closed DC mesh; if present, treat as sharp.)
  for (EI e : mesh.edges()) {
    if (PMP::is_border_edge(e, mesh)) sharp.insert(e);
  }
  return sharp;
}

// Utility: build triangle arrays and map (undirected edge) -> incident triangles
struct TriRef { int tri; int localEdge; }; // localEdge: 0:(v0,v1), 1:(v1,v2), 2:(v2,v0)

struct MeshIndexing {
  std::vector<K::Point_3> V;
  std::vector<std::array<int,3>> F;
  // map undirected edge -> up to two {tri, localEdge}
  std::unordered_map<uint64_t, std::array<TriRef,2>> edge2halves;
  std::unordered_map<uint64_t,int> eCount; // 1 or 2
};

static inline uint64_t key_undirected(int a, int b) {
  if (a>b) std::swap(a,b);
  return (uint64_t(a) << 32) | uint32_t(b);
}

MeshIndexing build_indexing(const SM& mesh) {
  MeshIndexing M;
  M.V.reserve(mesh.number_of_vertices());
  std::vector<int> vmap(mesh.num_vertices(), -1);

  // Map CGAL vertex ids to packed indices
  int vid = 0;
  for (VI v : mesh.vertices()) {
    M.V.push_back(mesh.point(v));
    vmap[v] = vid++;
  }

  M.F.reserve(mesh.number_of_faces());
  int tid = 0;
  for (FI f : mesh.faces()) {
    std::array<int,3> tri{-1,-1,-1};
    int i = 0;
    for (VI v : vertices_around_face(mesh.halfedge(f), mesh)) {
      tri[i++] = vmap[v];
    }
    // ensure triangles (DC outputs triangles; if polygons appear, triangulate earlier)
    if (i == 3) {
      int t = (int)M.F.size();
      M.F.push_back(tri);

      auto addEdge = [&](int le, int a, int b){
        uint64_t k = key_undirected(a,b);
        int cnt = M.eCount[k];
        if (cnt == 0) { M.edge2halves[k][0] = {t, le}; M.eCount[k] = 1; }
        else if (cnt == 1) { M.edge2halves[k][1] = {t, le}; M.eCount[k] = 2; }
      };

      addEdge(0, tri[0], tri[1]);
      addEdge(1, tri[1], tri[2]);
      addEdge(2, tri[2], tri[0]);
      ++tid;
    }
  }
  return M;
}

// Step 3+4: Build Manifold mesh + sharpened edge list
manifold::MeshGL64 to_meshGL64(const MeshIndexing& M) {
  manifold::MeshGL64 m;
  m.vertPos.reserve(M.V.size());
  for (auto const& p : M.V) {
    m.vertPos.push_back(float(p.x()));
    m.vertPos.push_back(float(p.y()));
    m.vertPos.push_back(float(p.z()));
  }
  m.triVerts.reserve(M.F.size()*3);
  for (auto const& t : M.F) {
    m.triVerts.push_back(t[0]);
    m.triVerts.push_back(t[1]);
    m.triVerts.push_back(t[2]);
  }
  return m;
}

// Optional: map dihedral angle -> partial smoothness (0 sharp … 1 smooth)
inline double angle_to_smoothness(double ang_deg,
                                  double sharp_deg=60.0,
                                  double smooth_deg=25.0) {
  if (ang_deg <= smooth_deg) return 1.0;
  if (ang_deg >= sharp_deg)  return 0.0;
  double t = (ang_deg - smooth_deg) / (sharp_deg - smooth_deg);
  double s = 1.0 - (t*t*(3.0 - 2.0*t)); // smoothstep
  return s;
}

// Compute `sharpenedEdges` from CGAL sharp-edge set
std::vector<manifold::Smoothness> build_sharpened_edges(
    const SM& mesh,
    const MeshIndexing& M,
    const std::unordered_set<EI>& sharpCGAL,
    double sharpAngleDeg = 60.0,
    double fullySmoothBelowDeg = 25.0)
{
  // Precompute per-edge dihedral angles in degrees (optional; for partial smoothness)
  std::unordered_map<uint64_t, double> edgeAngleDeg;

  // CGAL gives us EI; convert to our (a,b) indexing
  for (EI e : mesh.edges()) {
    auto h = mesh.halfedge(e);
    VI va = source(h, mesh), vb = target(h, mesh);
    int a = (int)va; int b = (int)vb;  // Surface_mesh indices are small ints
    uint64_t k = key_undirected(a,b);

    // dihedral; if border, force 180 to make it sharp
    double ang = 180.0;
    if (!PMP::is_border_edge(e, mesh)) {
      // CGAL has dihedral angle utilities in PMP; fallback: compute via face normals.
      // Here we use the built-in:
      ang = CGAL::to_double(PMP::dihedral_angle(e, mesh)) * 180.0 / CGAL_PI;
      // dihedral_angle returns signed angle; we take abs
      ang = std::fabs(ang);
    }
    edgeAngleDeg[k] = ang;
  }

  // For each *sharp* CGAL edge, emit both halfedges with smoothness in [0,1]
  std::vector<manifold::Smoothness> out;
  out.reserve(sharpCGAL.size()*2);

  for (EI e : sharpCGAL) {
    auto h = mesh.halfedge(e);
    int a = (int)source(h, mesh);
    int b = (int)target(h, mesh);
    uint64_t k = key_undirected(a,b);

    auto itCnt = M.eCount.find(k);
    if (itCnt == M.eCount.end()) continue; // shouldn't happen
    int cnt = itCnt->second;

    double ang = edgeAngleDeg[k];
    double s   = angle_to_smoothness(ang, sharpAngleDeg, fullySmoothBelowDeg); // 0..1

    const auto& halves = M.edge2halves.at(k);
    // halfedge index in Manifold: 3*tri + localEdge
    auto push_half = [&](const TriRef& r){
      uint32_t he = 3u * (uint32_t)r.tri + (uint32_t)r.localEdge;
      out.push_back({ s, he });
    };
    push_half(halves[0]);
    if (cnt == 2) push_half(halves[1]);
  }

  // Optional: promote corners — if a vertex has >=3 incident sharp edges, force s=0 on all its incident halfedges
  // (Left to implement if you need stronger corners.)

  return out;
}

int main() {
  SM mesh;
  // TODO: load from your DC output (OFF/OBJ/PLY)

  repair_mesh(mesh);

  // Detect sharp edges in CGAL
  const double sharp_angle_deg = 60.0;
  auto sharpCGAL = detect_sharp_edges(mesh, sharp_angle_deg);

  // Build indexing for Manifold
  auto M = build_indexing(mesh);

  // Convert to Manifold MeshGL64
  auto m64 = to_meshGL64(M);

  // Build sharpenedEdges for Manifold::Smooth
  auto sharpened = build_sharpened_edges(mesh, M, sharpCGAL, /*sharp*/60.0, /*smooth*/25.0);

  // Manifold: Smooth → Refine → export
  auto smoothed   = manifold::Manifold::Smooth(m64, sharpened);
  auto refined    = smoothed.Refine(2);             // try 2 first
  auto tri_out    = refined.GetMeshGL64();          // or ToTriangleMesh(), depending on API
  // TODO: write tri_out as OBJ/PLY
  return 0;
}


