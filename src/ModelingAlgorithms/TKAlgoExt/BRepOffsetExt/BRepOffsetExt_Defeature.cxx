// Copyright (c) 2026 Timon Gentzsch
//
// This file is part of Open CASCADE Technology software library.
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
// distribution for complete text of the license and disclaimer of any warranty.
//
// Alternatively, this file may be used under the terms of Open CASCADE
// commercial license or contractual agreement.

#include <GeomAlgoExt_Cancel.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
// defeature: delete faces and heal.
//
// BRepAlgoAPI_Defeaturing removes the given faces and extends the adjacent ones
// until they close the gap. A feature it cannot close is left in place with
// only a warning, so a call can "succeed" and return the input; we compare
// against the request and fail when nothing went.
//
// It rebuilds the whole solid with a General Fuse for every feature, one after
// the other, which grows with the square of their number. A feature that meets
// the rest of the body only in holes of its neighbours (a hole, a pocket, a
// boss) needs none of that: the neighbour's own surface spans the gap, which is
// what extending it would give, so dropping those loops removes the feature
// (fillHoles). The others go to BRepAlgoAPI_Defeaturing.

#include <BRepToolsExt_Checks.hxx>
#include <BRepAlgoExt_FillHoles.hxx>
#include <BRepToolsExt_History.hxx>
#include <BRepToolsExt_Operations.hxx>

#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Defeaturing.hxx>
#include <BRepBndLib.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepLib.hxx>
#include <BRepTools.hxx>
#include <BRepTools_ReShape.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <BndLib_AddSurface.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <NCollection_IndexedMap.hxx>
#include <ShapeAnalysis.hxx>
#include <ShapeAnalysis_FreeBounds.hxx>
#include <ShapeExtend_WireData.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Wire.hxx>

#include <numeric>

namespace BRepExt
{

namespace
{

using Ancestors =
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>;
using ShapeSet = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;

} // namespace

bool holeOf(const TopoDS_Shape& wire, const TopoDS_Face& owner)
{
  for (TopoDS_Iterator wires(owner); wires.More(); wires.Next())
    if (wires.Value().IsSame(wire))
      return ShapeAnalysis::TotCross2D(new ShapeExtend_WireData(TopoDS::Wire(wires.Value()), false),
                                       owner)
             < 0;
  return false;
}

TopoDS_Face holeCap(const TopoDS_Face& owner, const TopoDS_Shape& wire)
{
  TopoDS_Face cap = TopoDS::Face(owner.EmptyCopied());
  for (TopoDS_Iterator wires(owner); wires.More(); wires.Next())
    if (wires.Value().IsSame(wire))
      BRep_Builder().Add(cap, wires.Value().Reversed());
  return cap;
}

std::vector<TopoDS_Shape> featuresOf(const TopoDS_Shape&              shape,
                                     const std::vector<TopoDS_Shape>& faces)
{
  Ancestors facesOf, wiresOf, ownerOf;
  TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, facesOf);
  TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_WIRE, wiresOf);
  TopExp::MapShapesAndAncestors(shape, TopAbs_WIRE, TopAbs_FACE, ownerOf);
  ShapeSet feature;
  for (const TopoDS_Shape& face : faces)
    feature.Add(face);
  for (int i = 1; i <= feature.Extent(); i++)
    for (TopExp_Explorer edges(feature(i), TopAbs_EDGE); edges.More(); edges.Next())
    {
      if (BRep_Tool::Degenerated(TopoDS::Edge(edges.Current())))
        continue;
      for (const TopoDS_Shape& other : facesOf.FindFromKey(edges.Current()))
      {
        if (feature.Contains(other))
          continue;
        const TopoDS_Face owner = TopoDS::Face(other.Oriented(TopAbs_FORWARD));
        bool              hole  = false;
        for (const TopoDS_Shape& wire : wiresOf.FindFromKey(edges.Current()))
          for (const TopoDS_Shape& holder : ownerOf.FindFromKey(wire))
            hole = hole || (holder.IsSame(owner) && holeOf(wire, owner));
        if (!hole)
          feature.Add(other);
      }
    }
  std::vector<TopoDS_Shape> out;
  for (int i = 1; i <= feature.Extent(); i++)
    out.push_back(feature(i));
  return out;
}

Filled fillHoles(const TopoDS_Shape& shape, const std::vector<TopoDS_Shape>& faces)
{
  Ancestors facesOf, wiresOf, ownerOf;
  TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, facesOf);
  TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_WIRE, wiresOf);
  TopExp::MapShapesAndAncestors(shape, TopAbs_WIRE, TopAbs_FACE, ownerOf);
  ShapeSet chosen, loops;
  for (const TopoDS_Shape& face : faces)
    chosen.Add(face);

  // Union-find over the chosen faces, then the loops as they are met; and what rules out
  // the feature a member falls into.
  std::vector<int> root(chosen.Extent());
  std::iota(root.begin(), root.end(), 0);
  std::vector<bool> refused(root.size(), false);
  auto              find = [&](int i) {
    while (root[i] != i)
      i = root[i] = root[root[i]];
    return i;
  };
  auto join = [&](int a, int b) { root[find(a)] = find(b); };
  // `owner` forward, so that its wires wind as it has them; their edges taken in the order
  // they join (not chained: ShapeExtend_WireData then follows BRepTools_WireExplorer).
  auto loopOf = [&](const TopoDS_Shape& wire, const TopoDS_Face& owner) {
    if (const int known = loops.FindIndex(wire))
      return chosen.Extent() + known - 1;
    loops.Add(wire);
    root.push_back(static_cast<int>(root.size()));
    bool hole = holeOf(wire, owner);
    for (TopoDS_Iterator edges(wire); hole && edges.More(); edges.Next())
      for (const TopoDS_Shape& face : facesOf.FindFromKey(edges.Value()))
        hole = hole && (face.IsSame(owner) || chosen.Contains(face));
    refused.push_back(!hole);
    return static_cast<int>(root.size()) - 1;
  };

  for (int i = 1; i <= chosen.Extent(); i++)
  {
    const TopoDS_Face& face = TopoDS::Face(chosen(i));
    for (TopExp_Explorer edges(face, TopAbs_EDGE); edges.More(); edges.Next())
    {
      const TopoDS_Edge& edge = TopoDS::Edge(edges.Current());
      if (BRep_Tool::Degenerated(edge))
        continue;
      const NCollection_List<TopoDS_Shape>& owners = facesOf.FindFromKey(edge);
      if (owners.Extent() == 1 && !BRep_Tool::IsClosed(edge, face))
        refused[i - 1] = true; // a free edge: taking the face off opens the body
      for (const TopoDS_Shape& other : owners)
      {
        if (const int j = chosen.FindIndex(other))
        {
          join(i - 1, j - 1);
          continue;
        }
        const TopoDS_Face owner = TopoDS::Face(other.Oriented(TopAbs_FORWARD));
        for (const TopoDS_Shape& wire : wiresOf.FindFromKey(edge))
          for (const TopoDS_Shape& holder : ownerOf.FindFromKey(wire))
            if (holder.IsSame(owner))
              join(i - 1, loopOf(wire, owner));
      }
    }
  }
  for (size_t i = 0; i < root.size(); i++)
    if (refused[i])
      refused[find(static_cast<int>(i))] = true;

  // A feature with no loop is not one that this takes off.
  const int         first = chosen.Extent();
  std::vector<bool> filled(root.size(), false);
  for (int l = 0; l < loops.Extent(); l++)
    filled[find(first + l)] = !refused[find(first + l)];
  // The loops to fill span at most these boxes of their faces' surfaces, measured on the
  // surface itself: the boxes of a spline's poles would take in much of the body.
  std::vector<Bnd_Box>      spans(loops.Extent());
  std::vector<TopoDS_Shape> owners(loops.Extent());
  bool                      any = false;
  for (int l = 0; l < loops.Extent(); l++)
  {
    owners[l] = ownerOf.FindFromKey(loops(l + 1)).First();
    if (!filled[find(first + l)])
      continue;
    any                     = true;
    const TopoDS_Face owner = TopoDS::Face(owners[l].Oriented(TopAbs_FORWARD));
    double            u0, u1, v0, v1;
    BRepTools::UVBounds(owner, TopoDS::Wire(loops(l + 1)), u0, u1, v0, v1);
    BndLib_AddSurface::AddOptimal(BRepAdaptor_Surface(owner),
                                  u0,
                                  u1,
                                  v0,
                                  v1,
                                  BRep_Tool::Tolerance(owner),
                                  spans[l]);
  }
  ShapeSet             all;
  std::vector<Bnd_Box> boxes;
  if (any)
  {
    TopExp::MapShapes(shape, TopAbs_FACE, all);
    boxes.resize(all.Extent());
    for (int f = 1; f <= all.Extent(); f++)
      BRepBndLib::Add(all(f), boxes[f - 1], false);
  }
  // Where the boxes meet, whether the filled loop and the face (or another filled loop) do:
  // the loop is filled with its owner's surface, so what touches that cap within their
  // tolerances would be cut through.
  std::vector<TopoDS_Face> caps(loops.Extent());
  const auto               capOf = [&](int l) -> const TopoDS_Face& {
    if (caps[l].IsNull())
      caps[l] = holeCap(TopoDS::Face(owners[l].Oriented(TopAbs_FORWARD)), loops(l + 1));
    return caps[l];
  };
  const auto meet = [&](const TopoDS_Face& cap, const TopoDS_Shape& face) {
    BRepExtrema_DistShapeShape distance(cap, face);
    return !distance.IsDone()
           || distance.Value()
                <= BRep_Tool::Tolerance(cap) + BRep_Tool::Tolerance(TopoDS::Face(face));
  };
  std::vector<bool> clash(root.size(), false);
  for (int l = 0; l < loops.Extent(); l++)
  {
    const int feature = find(first + l);
    if (!filled[feature])
      continue;
    for (int f = 1; f <= all.Extent() && !clash[feature]; f++)
    {
      const int j = chosen.FindIndex(all(f));
      if (!all(f).IsSame(owners[l]) && (j == 0 || find(j - 1) != feature)
          && !boxes[f - 1].IsOut(spans[l]) && meet(capOf(l), all(f)))
        clash[feature] = true;
    }
    for (int m = 0; m < loops.Extent() && !clash[feature]; m++)
    {
      const int other = find(first + m);
      if (other != feature && filled[other] && !spans[m].IsOut(spans[l])
          && meet(capOf(l), capOf(m)))
        clash[feature] = true;
    }
  }
  for (size_t i = 0; i < root.size(); i++)
    filled[i] = filled[i] && !clash[i];
  BRepTools_ReShape reshape;
  Filled            out;
  for (int i = 1; i <= chosen.Extent(); i++)
  {
    if (filled[find(i - 1)])
      reshape.Remove(chosen(i));
    else
      out.rest.push_back(chosen(i));
  }
  ShapeSet done;
  for (int l = 0; l < loops.Extent(); l++)
  {
    if (!filled[find(first + l)] || done.Contains(owners[l]))
      continue;
    done.Add(owners[l]);
    const TopoDS_Face owner = TopoDS::Face(owners[l].Oriented(TopAbs_FORWARD));
    TopoDS_Shape      whole = owner.EmptyCopied();
    for (TopoDS_Iterator wires(owner); wires.More(); wires.Next())
    {
      const int m = loops.FindIndex(wires.Value());
      if (m == 0 || !filled[find(first + m - 1)])
        BRep_Builder().Add(whole, wires.Value());
    }
    reshape.Replace(owner, whole);
  }
  out.shape   = reshape.Apply(shape);
  out.history = reshape.History();
  return out;
}

EvolutionData defeature(OcctKernel&                  kernel,
                        uint32_t                     shapeId,
                        const std::vector<uint32_t>& faceIds,
                        const std::vector<int>&      inputFaceHashes,
                        int                          hashUpperBound)
{
  const char* const               op    = "defeature";
  const TopoDS_Shape&             shape = get(kernel, shapeId);
  FaceHistory                     history(op, inputFaceHashes, hashUpperBound);
  const std::vector<TopoDS_Shape> faces  = subShapes(kernel, op, shape, faceIds, TopAbs_FACE);
  const Filled                    filled = fillHoles(shape, faces);
  if (!faces.empty() && filled.rest.empty())
  {
    history.collect(*filled.history, shape);
    return history.finish(store(kernel, filled.shape));
  }
  const bool sheet = !TopExp_Explorer(shape, TopAbs_SOLID).More();
  // A sheet is closed below as it came, and what fillHoles took off goes with the rest.
  const bool                       staged = !sheet && filled.rest.size() < faces.size();
  const std::vector<TopoDS_Shape>& rest   = staged ? filled.rest : faces;

  // Defeaturing requires a volume. Temporarily close planar openings of an
  // open shell, rebuild the feature, then remove the evolved temporary caps.
  // General nonplanar openings cannot be closed by this bounded path.
  TopoDS_Shape              input = staged ? filled.shape : shape;
  std::vector<TopoDS_Shape> caps;
  BRepBuilderAPI_Sewing     sewing(BRepExt::tolerance::sewing);
  if (sheet)
  {
    ShapeAnalysis_FreeBounds free(shape, BRepExt::tolerance::sewing, false, true);
    if (TopExp_Explorer(free.GetOpenWires(), TopAbs_EDGE).More())
      fail(op, unsupported_geometry, "shell has open boundary chains");
    sewing.Add(shape);
    for (TopExp_Explorer it(free.GetClosedWires(), TopAbs_WIRE); it.More(); it.Next())
    {
      BRepBuilderAPI_MakeFace cap(TopoDS::Wire(it.Current()), true);
      if (!cap.IsDone())
        fail(op, unsupported_geometry, "shell openings must be planar to remove blends");
      caps.push_back(cap.Face());
      sewing.Add(cap.Face());
    }
    sewing.Perform();
    if (sewing.SewedShape().ShapeType() != TopAbs_SHELL)
      fail(op, unsupported_geometry, "shell cannot be closed as one volume");
    auto solid = BRepBuilderAPI_MakeSolid(TopoDS::Shell(sewing.SewedShape())).Solid();
    if (!BRepLib::OrientClosedSolid(solid) || !BRepCheck_Analyzer(solid).IsValid())
      fail(op, unsupported_geometry, "temporary shell closure is invalid");
    input = solid;
  }
  auto sewn = [&](const TopoDS_Shape& face) {
    return sheet && sewing.IsModifiedSubShape(face) ? sewing.ModifiedSubShape(face) : face;
  };
  BRepAlgoAPI_Defeaturing maker;
  maker.SetShape(input);
  for (const TopoDS_Shape& face : rest)
    maker.AddFaceToRemove(sewn(face));
  maker.SetToFillHistory(true);
  // Cancellable: defeaturing is where a Delete Face can take long.
  maker.Build(progress(true)->Start());
  checkCancelled(op);
  if (!maker.IsDone() || maker.HasErrors())
  {
    fail(op, no_solution, "the neighbouring faces cannot be extended to close the gap");
  }

  size_t removed = faces.size() - rest.size();
  for (const TopoDS_Shape& face : rest)
    removed += maker.IsDeleted(sewn(face)) ? 1 : 0;
  if (removed == 0)
  {
    fail(op,
         nothing_removed,
         "none of the faces could be removed: their neighbours cannot be "
         "extended to close "
         "the gap. Select the whole feature (every face of the hole, pocket or "
         "blend)");
  }
  if (!sheet)
  {
    if (staged)
      history.collectComposed(*filled.history, maker, shape);
    else
      history.collect(maker, shape);
    return history.finish(store(kernel, maker.Shape()));
  }
  NCollection_List<TopoDS_Shape> capFaces;
  for (const auto& cap : caps)
  {
    const auto original = sewn(cap);
    if (maker.Modified(original).IsEmpty() && !maker.IsDeleted(original))
      capFaces.Append(original);
    for (const auto& f : maker.Modified(original))
      capFaces.Append(f);
  }
  BRepBuilderAPI_Sewing open(BRepExt::tolerance::sewing);
  for (TopExp_Explorer it(maker.Shape(), TopAbs_FACE); it.More(); it.Next())
  {
    bool temporary = false;
    for (const auto& cap : capFaces)
      if (it.Current().IsSame(cap))
        temporary = true;
    if (!temporary)
      open.Add(it.Current());
  }
  open.Perform();
  if (open.SewedShape().IsNull() || !BRepCheck_Analyzer(open.SewedShape()).IsValid())
    fail(op, no_solution, "rebuilt shell is invalid");
  for (TopExp_Explorer it(shape, TopAbs_FACE); it.More(); it.Next())
  {
    const auto old = it.Current();
    if (!history.tracks(old))
      continue;
    const auto                     middle = sewn(old);
    NCollection_List<TopoDS_Shape> next   = maker.Modified(middle), final;
    if (next.IsEmpty() && !maker.IsDeleted(middle))
      next.Append(middle);
    for (const auto& f : next)
      final.Append(open.IsModifiedSubShape(f) ? open.ModifiedSubShape(f) : f);
    if (final.IsEmpty())
      history.deleted(old);
    else
      history.modified(old, final);
  }
  return history.finish(store(kernel, open.SewedShape()));
}

} // namespace BRepExt
