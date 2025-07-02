
#include <BOPDS_Pair.hxx>

struct BOPDS_PairHasher {
  static Standard_Integer HashCode(const BOPDS_Pair& theKey, const Standard_Integer theUpper)
  {
    return ::HashCode(theKey.Index1() * 31 + theKey.Index2(), theUpper);
  }

  static Standard_Boolean IsEqual(const BOPDS_Pair& k1, const BOPDS_Pair& k2)
  {
    return k1.Index1() == k2.Index1() && k1.Index2() == k2.Index2();
  }
};

