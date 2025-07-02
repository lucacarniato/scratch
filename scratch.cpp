
#include <BOPDS_Pair.hxx>
#include <functional> // for std::hash

struct BOPDS_PairHasher {
  static Standard_Integer HashCode(const BOPDS_Pair& theKey, const Standard_Integer theUpper) {
    std::size_t h1 = std::hash<int>{}(theKey.Index1());
    std::size_t h2 = std::hash<int>{}(theKey.Index2());

    // Combine hashes using a simple mixing formula (e.g. boost-style)
    std::size_t combined = h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));

    // Map to OpenCascade hash range
    return static_cast<Standard_Integer>(combined % theUpper) + 1;
  }

  static Standard_Boolean IsEqual(const BOPDS_Pair& k1, const BOPDS_Pair& k2) {
    return k1.Index1() == k2.Index1() && k1.Index2() == k2.Index2();
  }
};


