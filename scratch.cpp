// 0) (Optional) pre-simplify if your B-Rep came from triangles:
ShapeUpgrade_UnifySameDomain uL(left,  true, true, false);  uL.Build();
ShapeUpgrade_UnifySameDomain uR(right, true, true, false);  uR.Build();
TopoDS_Shape Ls = uL.Shape(), Rs = uR.Shape();

// 1) Explode into solids/components
TopTools_ListOfShape Lparts = ExplodeSolidsOrComponents(Ls);
TopTools_ListOfShape Rparts = ExplodeSolidsOrComponents(Rs);

// 2) Precompute AABBs/BVH for Rparts (your helper)
auto Rbins = BuildSpatialIndex(Rparts);

BRep_Builder bb;
TopoDS_Compound out; bb.MakeCompound(out);

for (TopTools_ListIteratorOfListOfShape it(Lparts); it.More(); it.Next()) {
  const TopoDS_Shape& SL = it.Value();

  // pick only overlapping right components (+optional small halo)
  TopTools_ListOfShape Rchunk = PickOverlapping(SL, Rbins);

  if (Rchunk.IsEmpty()) { bb.Add(out, SL); continue; }

  // --- build per-batch argument list
  TopTools_ListOfShape batchArgs; batchArgs.Append(SL);
  for (TopTools_ListIteratorOfListOfShape ir(Rchunk); ir.More(); ir.Next())
    batchArgs.Append(ir.Value());

  // --- pave only this small batch
  BOPAlgo_PaveFiller pf;
  pf.SetArguments(batchArgs);
  pf.SetNonDestructive(false);
  pf.SetRunParallel(false);
  pf.SetUseOBB(true);
  pf.SetFuzzyValue(1e-4);
  pf.Perform();

  // --- build cells from the same batch arguments
  BOPAlgo_CellsBuilder cb;
  cb.SetArguments(batchArgs);
  cb.PerformWithFiller(pf);

  bb.Add(out, cb.Shape());

  // free temps from this batch before moving on
  cb.Clear();
  pf.Clear();
}

// Optional: clean seams from chunking
ShapeUpgrade_UnifySameDomain uOut(out, true, true, false); uOut.Build();
TopoDS_Shape finalShape = uOut.Shape();





TopTools_ListOfShape ExplodeSolids(const TopoDS_Shape& shape) {
  TopTools_ListOfShape solids;
  for (TopExp_Explorer ex(shape, TopAbs_SOLID); ex.More(); ex.Next()) {
    solids.Append(ex.Current());
  }
  if (solids.IsEmpty())
    solids.Append(shape);  // fallback if it's a single shell
  return solids;
}

struct ShapeBox {
  TopoDS_Shape shape;
  Bnd_Box box;
};

std::vector<ShapeBox> BuildSpatialIndex(const TopTools_ListOfShape& shapes) {
  std::vector<ShapeBox> bins;
  for (TopTools_ListIteratorOfListOfShape it(shapes); it.More(); it.Next()) {
    ShapeBox sb;
    sb.shape = it.Value();
    BRepBndLib::Add(sb.shape, sb.box);
    sb.box.SetGap(1.0e-4); // optional margin
    bins.push_back(sb);
  }
  return bins;
}


TopTools_ListOfShape PickOverlapping(const TopoDS_Shape& solid,
                                     const std::vector<ShapeBox>& bins,
                                     double halo = 0.0)
{
  Bnd_Box boxL;
  BRepBndLib::Add(solid, boxL);
  if (halo > 0.0) boxL.Enlarge(halo);

  TopTools_ListOfShape result;
  for (const auto& sb : bins) {
    if (!boxL.IsOut(sb.box)) {
      result.Append(sb.shape);
    }
  }
  return result;
}


auto Lsolids = ExplodeSolids(left);
auto Rsolids = ExplodeSolids(right);
auto Rbins   = BuildSpatialIndex(Rsolids);

for (TopTools_ListIteratorOfListOfShape it(Lsolids); it.More(); it.Next()) {
  const TopoDS_Shape& SL = it.Value();
  TopTools_ListOfShape Rchunk = PickOverlapping(SL, Rbins, /*halo*/ 1e-3);
  if (Rchunk.IsEmpty()) continue;

  // pave filler + cells builder on (SL + Rchunk)
  ...
}





auto Lsolids = ExplodeSolids(left);
auto Rsolids = ExplodeSolids(right);
auto Rbins   = BuildSpatialIndex(Rsolids);

for (TopTools_ListIteratorOfListOfShape it(Lsolids); it.More(); it.Next()) {
  const TopoDS_Shape& SL = it.Value();
  TopTools_ListOfShape Rchunk = PickOverlapping(SL, Rbins, /*halo*/ 1e-3);
  if (Rchunk.IsEmpty()) continue;

  // pave filler + cells builder on (SL + Rchunk)
  ...
}



// accumulators
  BRep_Builder bb;
  TopoDS_Compound cutOut;    bb.MakeCompound(cutOut);
  TopoDS_Compound commonOut; bb.MakeCompound(commonOut);

  for (TopTools_ListIteratorOfListOfShape it(Lparts); it.More(); it.Next()) {
    const TopoDS_Shape& SL = it.Value();
    TopTools_ListOfShape Rchunk = PickOverlapping(SL, rItems, rTree, /*halo*/ halo);

    if (Rchunk.IsEmpty()) {
      // nothing overlaps: whole left solid goes to the CUT result, no COMMON
      bb.Add(cutOut, SL);
      continue;
    }

    // --- build batch arguments
    TopTools_ListOfShape args;  args.Append(SL);          // left
    TopTools_ListOfShape tools;                           // right chunk
    for (TopTools_ListIteratorOfListOfShape jr(Rchunk); jr.More(); jr.Next())
      tools.Append(jr.Value());

    // --- pave once
    BOPAlgo_PaveFiller pf;
    pf.SetNonDestructive(nonDestr);
    pf.SetRunParallel(runParallel);
    pf.SetUseOBB(useOBB);
    pf.SetFuzzyValue(fuzzy);

    TopTools_ListOfShape batchForFiller = args;
    for (TopTools_ListIteratorOfListOfShape jt(tools); jt.More(); jt.Next())
      batchForFiller.Append(jt.Value());
    pf.SetArguments(batchForFiller);
    pf.Perform();

    // --- CUT: L \ R (reuse filler)
    BOPAlgo_BOP opCut;
    opCut.SetOperation(BOPAlgo_CUT);
    opCut.SetNonDestructive(nonDestr);
    opCut.SetRunParallel(runParallel);
    opCut.SetUseOBB(useOBB);
    opCut.SetFuzzyValue(fuzzy);
    opCut.SetArguments(args);
    opCut.SetTools(tools);
    opCut.PerformWithFiller(pf);
    bb.Add(cutOut, opCut.Shape());
    opCut.Clear();

    // --- COMMON: L ∩ R (reuse the same filler again)
    BOPAlgo_BOP opCommon;
    opCommon.SetOperation(BOPAlgo_COMMON);
    opCommon.SetNonDestructive(nonDestr);
    opCommon.SetRunParallel(runParallel);
    opCommon.SetUseOBB(useOBB);
    opCommon.SetFuzzyValue(fuzzy);
    opCommon.SetArguments(args);
    opCommon.SetTools(tools);
    opCommon.PerformWithFiller(pf);
    bb.Add(commonOut, opCommon.Shape());
    opCommon.Clear();

    // free pave DS before next batch
    pf.Clear();
  }

  // optional cleanup to remove batch seams / duplicates
  TopoDS_Shape cutRes    = cutOut;
  TopoDS_Shape commonRes = commonOut;
  if (unifyOut) {
    ShapeUpgrade_UnifySameDomain uC(cutRes,    Standard_True, Standard_True, Standard_False);    uC.Build();    cutRes    = uC.Shape();
    ShapeUpgrade_UnifySameDomain uI(commonRes, Standard_True, Standard_True, Standard_False);    uI.Build();    commonRes = uI.Shape();
  }

  return { cutRes, commonRes };
}




