// ================================================================
// Crease-aware SDF reprojection for Manifold MeshGL64 + OpenVDB SDF
// ================================================================

#include <manifold/meshGL.h>     // MeshGL64
#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cstdint>

// OpenVDB
#include <openvdb/openvdb.h>
#include <openvdb/tools/GridSampler.h>

// Optional CGAL global self-intersection pass
#ifdef USE_CGAL_SELF_INTERSECTION
  #include <CGAL/Simple_cartesian.h>
  #include <CGAL/Surface_mesh.h>
  #include <CGAL/Polygon_mesh_processing/self_intersections.h>
  #include <set>
  namespace PMP = CGAL::Polygon_mesh_processing;
  using Kernel = CGAL::Simple_cartesian<double>;
  using SurfMesh = CGAL::Surface_mesh<Kernel::Point_3>;
#endif

// ---------- Parameters ----------
struct ProjectionParams {
  int    newtonIters        = 5;     // Newton steps per vertex
  double tolVox             = 0.25;  // stop when |phi| <= tolVox * voxelMin
  double clampFrac          = 0.3;   // clamp step to clampFrac * minAdjEdgeLen
  double dihedralSharpDeg   = 60.0;  // >= → sharp edge
  double dihedralSmoothDeg  = 25.0;  // <= → fully smooth (for partial creases)
  double maxNormalFlipDeg   = 80.0;  // guard: new normal must not flip beyond this
  double minAreaFactor      = 0.2;   // guard: area must not drop below this factor
  int    lineSearchMax      = 5;     // backtracking halving steps
};

// ---------- Small math helpers ----------
static inline void normalizeVec(double v[3]) {
  const double n = std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  if (n > 0) { v[0]/=n; v[1]/=n; v[2]/=n; }
}
static inline double dot3(const double a[3], const double b[3]) {
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
static inline void sub3(const double a[3], const double b[3], double out[3]) {
  out[0]=a[0]-b[0]; out[1]=a[1]-b[1]; out[2]=a[2]-b[2];
}
static inline void add3(const double a[3], const double b[3], double out[3]) {
  out[0]=a[0]+b[0]; out[1]=a[1]+b[1]; out[2]=a[2]+b[2];
}
static inline void mul3(double s, const double a[3], double out[3]) {
  out[0]=s*a[0]; out[1]=s*a[1]; out[2]=s*a[2];
}
static inline void cross3(const double a[3], const double b[3], double out[3]) {
  out[0] = a[1]*b[2]-a[2]*b[1];
  out[1] = a[2]*b[0]-a[0]*b[2];
  out[2] = a[0]*b[1]-a[1]*b[0];
}
static inline double norm3(const double a[3]) {
  return std::sqrt(dot3(a,a));
}

// ---------- Mesh helpers ----------
struct Tri { uint32_t v[3]; };
static inline uint64_t edgeKey(uint32_t a, uint32_t b) {
  if (a>b) std::swap(a,b);
  return (uint64_t(a)<<32) | b;
}

struct LocalTopology {
  std::vector<Tri>    F;
  std::vector<std::vector<uint32_t>> v2Tris;
  std::vector<std::vector<uint32_t>> v2Nbrs;
  std::unordered_map<uint64_t, std::array<uint32_t,2>> edge2Tris; // up to 2 tris
  std::unordered_map<uint64_t, int> edgeCount;
  std::vector<std::array<double,3>> faceN;   // unit normals
  std::vector<double>               faceA;   // area
  std::vector<double>               vMinEdge;// per-vertex min adjacent edge length
};

static LocalTopology buildTopology(const manifold::MeshGL64& mesh) {
  LocalTopology topo;
  const size_t nV = mesh.vertPos.size()/3;
  const size_t nT = mesh.triVerts.size()/3;

  topo.F.resize(nT);
  for (size_t t=0; t<nT; ++t) {
    topo.F[t].v[0] = mesh.triVerts[3*t+0];
    topo.F[t].v[1] = mesh.triVerts[3*t+1];
    topo.F[t].v[2] = mesh.triVerts[3*t+2];
  }

  topo.v2Tris.assign(nV, {});
  topo.v2Nbrs.assign(nV, {});
  topo.faceN.resize(nT);
  topo.faceA.resize(nT, 0.0);
  topo.vMinEdge.assign(nV, std::numeric_limits<double>::infinity());

  auto V = [&](uint32_t i){ return &mesh.vertPos[3*i]; };

  // Fill adjacency and face normals/areas
  for (uint32_t t=0; t<nT; ++t) {
    const auto& tri = topo.F[t];
    for (int k=0;k<3;++k) topo.v2Tris[tri.v[k]].push_back(t);

    // neighbors + edge lengths
    for (int k=0;k<3;++k) {
      uint32_t a = tri.v[k], b = tri.v[(k+1)%3];
      topo.v2Nbrs[a].push_back(b);
      topo.v2Nbrs[b].push_back(a);

      double e[3]; sub3(V(a), V(b), e);
      double l = norm3(e);
      topo.vMinEdge[a] = std::min(topo.vMinEdge[a], l);
      topo.vMinEdge[b] = std::min(topo.vMinEdge[b], l);

      auto kkey = edgeKey(a,b);
      int cnt = topo.edgeCount[kkey];
      if (cnt==0) { topo.edge2Tris[kkey][0]=t; topo.edgeCount[kkey]=1; }
      else if (cnt==1) { topo.edge2Tris[kkey][1]=t; topo.edgeCount[kkey]=2; }
    }

    // face normal + area
    double e0[3], e1[3], n[3];
    sub3(V(tri.v[1]), V(tri.v[0]), e0);
    sub3(V(tri.v[2]), V(tri.v[0]), e1);
    cross3(e0,e1,n);
    double A2 = norm3(n);
    topo.faceA[t] = 0.5 * A2;
    if (A2>0) { n[0]/=A2; n[1]/=A2; n[2]/=A2; }
    topo.faceN[t] = {n[0],n[1],n[2]};
  }

  // Dedup neighbors
  for (auto& nbrs : topo.v2Nbrs) {
    std::sort(nbrs.begin(), nbrs.end());
    nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
  }
  return topo;
}

// Dihedral in degrees between two triangle normals
static inline double dihedralDeg(const std::array<double,3>& n0,
                                 const std::array<double,3>& n1) {
  double c = std::clamp(n0[0]*n1[0]+n0[1]*n1[1]+n0[2]*n1[2], -1.0, 1.0);
  return std::acos(c) * 180.0 / M_PI;
}

// Tag edges/vertices as features
struct FeatureTags {
  std::unordered_set<uint64_t> sharpEdges; // undirected keys
  std::vector<uint8_t> isCorner;           // >=3 incident sharp edges
  std::vector<uint8_t> isCrease;           // exactly 2 incident sharp edges
  std::vector<std::array<double,3>> creaseTangent; // for crease vertices
};

static FeatureTags detectFeatures(const manifold::MeshGL64& mesh,
                                  const LocalTopology& topo,
                                  double sharpDeg)
{
  FeatureTags tags;
  const size_t nV = mesh.vertPos.size()/3;
  tags.isCorner.assign(nV, 0);
  tags.isCrease.assign(nV, 0);
  tags.creaseTangent.assign(nV, {0,0,0});
  auto V = [&](uint32_t i){ return &mesh.vertPos[3*i]; };

  // Edge sharpness by dihedral
  for (const auto& kv : topo.edge2Tris) {
    const uint64_t k = kv.first;
    int cnt = topo.edgeCount.at(k);
    if (cnt<2) { // border -> treat as sharp
      tags.sharpEdges.insert(k);
      continue;
    }
    uint32_t t0 = kv.second[0], t1 = kv.second[1];
    double ang = dihedralDeg(topo.faceN[t0], topo.faceN[t1]);
    if (ang >= sharpDeg) tags.sharpEdges.insert(k);
  }

  // Per-vertex classification & crease tangent
  for (uint32_t v=0; v<nV; ++v) {
    // count sharp incident edges & record their directions
    int sharpCount = 0;
    std::array<double,3> tavg{0,0,0};
    for (uint32_t nb : topo.v2Nbrs[v]) {
      if (tags.sharpEdges.count(edgeKey(v, nb))) {
        ++sharpCount;
        double dir[3]; sub3(V(nb), V(v), dir);
        double n = norm3(dir); if (n>0) { dir[0]/=n; dir[1]/=n; dir[2]/=n; }
        tavg[0]+=dir[0]; tavg[1]+=dir[1]; tavg[2]+=dir[2];
      }
    }
    if (sharpCount>=3) tags.isCorner[v] = 1;
    else if (sharpCount==2) {
      tags.isCrease[v] = 1;
      normalizeVec(tavg.data());
      tags.creaseTangent[v] = tavg;
    }
  }
  return tags;
}

// ---------- OpenVDB sampling ----------
struct SDFSampler {
  using GridT = openvdb::FloatGrid;
  using TreeT = GridT::TreeType;
  openvdb::tools::GridSampler<TreeT, openvdb::tools::BoxSampler> sampler;
  openvdb::math::Vec3d voxel;
  double voxelMin;

  explicit SDFSampler(openvdb::FloatGrid::ConstPtr sdf)
    : sampler(sdf->tree(), sdf->transform()),
      voxel(sdf->voxelSize())
  {
    voxelMin = std::min({voxel.x(), voxel.y(), voxel.z()});
  }

  inline double phiWS(const openvdb::math::Vec3d& p) const {
    return sampler.wsSample(p);
  }

  // Central differences in WORLD space with h = 0.5 * voxelMin
  inline openvdb::math::Vec3d gradWS(const openvdb::math::Vec3d& p) const {
    const double h = 0.5 * voxelMin;
    const openvdb::math::Vec3d ex(h,0,0), ey(0,h,0), ez(0,0,h);
    const double px = phiWS(p+ex), mx = phiWS(p-ex);
    const double py = phiWS(p+ey), my = phiWS(p-ey);
    const double pz = phiWS(p+ez), mz = phiWS(p-ez);
    return openvdb::math::Vec3d((px-mx)/(2*h), (py-my)/(2*h), (pz-mz)/(2*h));
  }
};

// ---------- Guards ----------
static bool preservesIncidentFaces(const manifold::MeshGL64& mesh,
                                   const LocalTopology& topo,
                                   uint32_t vIdx,
                                   const double newPos[3],
                                   double cosMaxFlip, double minAreaFactor)
{
  const auto& tris = topo.v2Tris[vIdx];
  auto V = [&](uint32_t i){ return &mesh.vertPos[3*i]; };

  for (uint32_t t : tris) {
    const auto& tri = topo.F[t];
    // copy verts, swap vIdx with newPos
    const double* p0 = V(tri.v[0]);
    const double* p1 = V(tri.v[1]);
    const double* p2 = V(tri.v[2]);
    double a[3]={p0[0],p0[1],p0[2]};
    double b[3]={p1[0],p1[1],p1[2]};
    double c[3]={p2[0],p2[1],p2[2]};
    if (tri.v[0]==vIdx) { a[0]=newPos[0]; a[1]=newPos[1]; a[2]=newPos[2]; }
    if (tri.v[1]==vIdx) { b[0]=newPos[0]; b[1]=newPos[1]; b[2]=newPos[2]; }
    if (tri.v[2]==vIdx) { c[0]=newPos[0]; c[1]=newPos[1]; c[2]=newPos[2]; }

    double e0[3], e1[3], nNew[3];
    sub3(b,a,e0); sub3(c,a,e1); cross3(e0,e1,nNew);
    const double A2 = norm3(nNew);
    if (A2==0.0) return false;
    double nNewUnit[3] = {nNew[0]/A2, nNew[1]/A2, nNew[2]/A2};

    const auto& nOld = topo.faceN[t];
    double cosang = nOld[0]*nNewUnit[0]+nOld[1]*nNewUnit[1]+nOld[2]*nNewUnit[2];
    if (cosang < cosMaxFlip) return false; // flipped too much

    const double Aold = topo.faceA[t];
    const double Anew = 0.5 * A2;
    if (Anew < minAreaFactor * Aold) return false;
  }
  return true;
}

// ---------- Main projection ----------
void ProjectMeshToSDF_ManifoldOpenVDB(
    manifold::MeshGL64& mesh,                       // in/out
    openvdb::FloatGrid::ConstPtr sdf,               // in
    const ProjectionParams& P = ProjectionParams{}) // params
{
  openvdb::initialize(); // safe if already initialized

  const size_t nV = mesh.vertPos.size()/3;
  if (nV==0 || mesh.triVerts.empty()) return;

  // Build topology and features
  LocalTopology topo = buildTopology(mesh);
  FeatureTags   tags = detectFeatures(mesh, topo, P.dihedralSharpDeg);

  SDFSampler s(sdf);
  const double tol = P.tolVox * s.voxelMin;
  const double cosMaxFlip = std::cos(P.maxNormalFlipDeg * M_PI/180.0);

  // Cache previous positions for potential global backoff
  std::vector<std::array<double,3>> prevPos(nV);
  for (size_t i=0;i<nV;++i) {
    prevPos[i] = {mesh.vertPos[3*i+0], mesh.vertPos[3*i+1], mesh.vertPos[3*i+2]};
  }

  // Per-vertex projection
  for (uint32_t vi=0; vi<nV; ++vi) {
    if (tags.isCorner[vi]) continue; // keep corners fixed

    double* p = &mesh.vertPos[3*vi];
    openvdb::math::Vec3d pw(p[0],p[1],p[2]);

    for (int it=0; it<P.newtonIters; ++it) {
      const double phi = s.phiWS(pw);
      if (std::abs(phi) <= tol) break;

      openvdb::math::Vec3d g = s.gradWS(pw);
      const double g2 = g.lengthSqr();
      if (!(g2>1e-14)) break; // avoid division by tiny gradient

      // Newton step along normal
      openvdb::math::Vec3d dn = - (phi / g2) * g;

      // Crease constraint: remove along-crease component
      if (tags.isCrease[vi]) {
        const auto& t = tags.creaseTangent[vi]; // unit (approx)
        openvdb::math::Vec3d T(t[0],t[1],t[2]);
        const double along = dn.dot(T);
        dn -= along * T;
      }

      // Clamp to local scale
      const double maxStep = P.clampFrac * topo.vMinEdge[vi];
      const double dlen = dn.length();
      if (dlen > maxStep && dlen>0) dn *= (maxStep / dlen);

      // Backtracking line search with local guards
      double alpha = 1.0;
      bool accepted = false;
      for (int ls=0; ls<P.lineSearchMax; ++ls) {
        openvdb::math::Vec3d trial = pw + alpha * dn;

        double np[3] = {trial.x(), trial.y(), trial.z()};
        // local orientation/area guard
        if (!preservesIncidentFaces(mesh, topo, vi, np, cosMaxFlip, P.minAreaFactor)) {
          alpha *= 0.5;
          continue;
        }

        const double phiNew = s.phiWS(trial);
        if (std::abs(phiNew) < std::abs(phi)) {
          // accept
          pw = trial;
          accepted = true;
          break;
        } else {
          alpha *= 0.5;
        }
      }
      if (!accepted) break; // give up on this vertex this iteration

      // Write back position
      p[0] = pw.x(); p[1] = pw.y(); p[2] = pw.z();
    }
  }

#ifdef USE_CGAL_SELF_INTERSECTION
  // -------- Global self-intersection pass (optional) --------
  auto buildSM = [&](){
    SurfMesh sm;
    std::vector<SurfMesh::Vertex_index> smV(mesh.vertPos.size()/3);
    for (size_t i=0;i<mesh.vertPos.size()/3;++i) {
      smV[i] = sm.add_vertex(Kernel::Point_3(
        mesh.vertPos[3*i+0], mesh.vertPos[3*i+1], mesh.vertPos[3*i+2]));
    }
    for (size_t t=0;t<mesh.triVerts.size()/3;++t) {
      const uint32_t a = mesh.triVerts[3*t+0];
      const uint32_t b = mesh.triVerts[3*t+1];
      const uint32_t c = mesh.triVerts[3*t+2];
      if (a!=b && b!=c && a!=c) sm.add_face(smV[a], smV[b], smV[c]);
    }
    return sm;
  };

  // try a few backoff rounds if intersections exist
  for (int round=0; round<3; ++round) {
    SurfMesh sm = buildSM();
    std::set<std::pair<SurfMesh::Face_index, SurfMesh::Face_index>> hits;
    PMP::self_intersections(sm, std::inserter(hits, hits.end()));

    if (hits.empty()) break;

    // Collect involved vertices and back off their displacement by 50%
    std::vector<uint8_t> mark(nV, 0);
    for (auto const& pr : hits) {
      for (int k=0;k<3;++k) {
        auto v = vertices_around_face(halfedge(pr.first, sm), sm)[k];
        // Map CGAL indices back: we created vertices in order
        const size_t idx = (size_t)v;
        if (idx < nV) mark[idx] = 1;
      }
      for (int k=0;k<3;++k) {
        auto v = vertices_around_face(halfedge(pr.second, sm), sm)[k];
        const size_t idx = (size_t)v;
        if (idx < nV) mark[idx] = 1;
      }
    }
    for (size_t i=0;i<nV;++i) if (mark[i]) {
      double* p = &mesh.vertPos[3*i];
      // back off towards previous position
      p[0] = 0.5*p[0] + 0.5*prevPos[i][0];
      p[1] = 0.5*p[1] + 0.5*prevPos[i][1];
      p[2] = 0.5*p[2] + 0.5*prevPos[i][2];
    }
    // update prevPos for next round
    for (size_t i=0;i<nV;++i) prevPos[i] = {mesh.vertPos[3*i], mesh.vertPos[3*i+1], mesh.vertPos[3*i+2]};
  }
#endif
}
