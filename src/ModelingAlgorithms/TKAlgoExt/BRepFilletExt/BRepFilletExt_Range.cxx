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

#include <BRepToolsExt_Checks.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_History.hxx>
#include <BRepFilletExt_Untouched.hxx>
#include <BRepFilletExt_Blend.hxx>
#include <Standard_Failure.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepTools_ReShape.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRep_Builder.hxx>
#include <BOPAlgo_PaveFiller.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeHalfSpace.hxx>
#include <GProp_GProps.hxx>
#include <TopoDS.hxx>
#include <cmath>
#include <NCollection_Array1.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Pln.hxx>
#include <BRep_Tool.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopTools_ShapeMapHasher.hxx>

namespace BRepExt
{
//! Carry intent through Boolean edge modifications, retaining only final edges.
template<class Maker>
static std::vector<blending::Join> rangeJoins(const std::vector<blending::Join>& before,
                                            Maker& maker)
{
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> retained;
  TopExp::MapShapes(maker.Shape(), TopAbs_EDGE, retained);
  std::vector<blending::Join> after;
  for (const auto& join : before)
  {
    if (retained.Contains(join.edge))
      after.push_back(join);
    for (const auto& edge : maker.Modified(join.edge))
      if (edge.ShapeType() == TopAbs_EDGE && retained.Contains(edge))
        after.push_back({edge, join.role, join.order});
  }
  return after;
}

//! True when `shape` is one solid bounded by one shell.
static bool oneShell(const TopoDS_Shape& shape)
{
  int solids = 0, shells = 0;
  for (TopExp_Explorer solid(shape, TopAbs_SOLID); solid.More(); solid.Next())
    ++solids;
  for (TopExp_Explorer shell(shape, TopAbs_SHELL); shell.More(); shell.Next())
    ++shells;
  return solids == 1 && shells == 1;
}

//! Assemble the signed boundary difference from support-face provenance.
//! Unchanged faces cancel. The remaining skins share one intersection inventory;
//! new sheets enter with reversed orientation in the removed-material boundary.
//! Each connected component is then a whole outer shell of the removed or the added
//! material, told apart by its signed volume. That holds when both bodies are single
//! shells sharing an unchanged face: neither difference then has a cavity shell
//! (derived; see the operations documentation).
struct RangeDelta
{
  TopoDS_Compound             removed, added;
  std::vector<blending::Join> joins;
};

static RangeDelta localDelta(const char*              op,
                             const TopoDS_Shape&      shape,
                             const blending::Blended& made)
{
  RangeDelta   delta;
  BRep_Builder builder;
  builder.MakeCompound(delta.removed);
  builder.MakeCompound(delta.added);
  TopoDS_Compound oldSkin, newSkin;
  builder.MakeCompound(oldSkin);
  builder.MakeCompound(newSkin);
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> cancelled;
  for (const auto& [old, now] : made.unchanged)
  {
    cancelled.Add(old);
    cancelled.Add(now);
  }
  for (TopExp_Explorer face(shape, TopAbs_FACE); face.More(); face.Next())
    if (!cancelled.Contains(face.Current()))
      builder.Add(oldSkin, face.Current());
  for (TopExp_Explorer face(made.shape, TopAbs_FACE); face.More(); face.Next())
    if (!cancelled.Contains(face.Current()))
      builder.Add(newSkin, face.Current());

  // Both oriented boundary differences must share the split edge inventory.
  NCollection_List<TopoDS_Shape> arguments;
  arguments.Append(oldSkin);
  arguments.Append(newSkin);
  BOPAlgo_PaveFiller intersections;
  intersections.SetArguments(arguments);
  intersections.SetNonDestructive(true);
  intersections.Perform();
  BRepAlgoAPI_Cut oldOnly(oldSkin, newSkin, intersections);
  BRepAlgoAPI_Cut newOnly(newSkin, oldSkin, intersections);
  if (!oldOnly.IsDone() || !newOnly.IsDone())
    fail(op, no_solution, "support boundary difference failed");

  // Substitute only retained source-face pieces. No sewing is allowed to
  // rewrite the surfaces or the edge representation lists during assembly.
  BRepTools_ReShape reshape;
  for (TopExp_Explorer face(shape, TopAbs_FACE); face.More(); face.Next())
  {
    if (cancelled.Contains(face.Current()) || oldOnly.IsDeleted(face.Current()))
      reshape.Remove(face.Current());
    else if (!oldOnly.Modified(face.Current()).IsEmpty())
    {
      TopoDS_Compound images;
      builder.MakeCompound(images);
      for (const auto& image : oldOnly.Modified(face.Current()))
        if (image.ShapeType() == TopAbs_FACE)
          builder.Add(images, image);
      reshape.Replace(face.Current(), images);
    }
  }
  std::vector<TopoDS_Shape> patches;
  const auto                retained = reshape.Apply(shape);
  for (TopExp_Explorer face(retained, TopAbs_FACE); face.More(); face.Next())
    patches.push_back(face.Current());
  for (TopExp_Explorer face(newOnly.Shape(), TopAbs_FACE); face.More(); face.Next())
    patches.push_back(face.Current().Reversed());

  // Connected components are determined by topological edges, never by a
  // tolerance merge. Each component carries the signed boundary orientation.
  std::vector<bool> used(patches.size(), false);
  for (size_t seed = 0; seed < patches.size(); ++seed)
  {
    if (used[seed])
      continue;
    TopoDS_Shell shell;
    builder.MakeShell(shell);
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> edges;
    std::vector<size_t>                                        pending{seed};
    used[seed] = true;
    for (size_t cursor = 0; cursor < pending.size(); ++cursor)
    {
      const auto& face = patches[pending[cursor]];
      builder.Add(shell, face);
      TopExp::MapShapes(face, TopAbs_EDGE, edges);
      for (size_t i = 0; i < patches.size(); ++i)
      {
        if (used[i])
          continue;
        for (TopExp_Explorer edge(patches[i], TopAbs_EDGE); edge.More(); edge.Next())
          if (edges.Contains(edge.Current()))
          {
            used[i] = true;
            pending.push_back(i);
            break;
          }
      }
    }
    auto solid = BRepBuilderAPI_MakeSolid(shell).Solid();
    if (!BRepCheck_Analyzer(solid).IsValid())
      fail(op, no_solution, "material boundary is not closed");
    GProp_GProps mass;
    BRepGProp::VolumeProperties(solid, mass);
    const bool removed = mass.Mass() > 0;
    if (!removed)
      solid.Reverse();
    builder.Add(removed ? delta.removed : delta.added, solid);
  }
  delta.joins         = rangeJoins(made.joins, oldOnly);
  const auto newJoins = rangeJoins(made.joins, newOnly);
  delta.joins.insert(delta.joins.end(), newJoins.begin(), newJoins.end());
  return delta;
}

//! The part of the change `made` makes to `shape` that lies between the planes square to `edge`
//! at the fractions `start` and `end` of it (the contour's other edges, `others`, must not come
//! back between them): the material removed and added there only, the rest of `shape` as it was.
//! The faces' history goes to `history`; the result is stored.
static EvolutionData withinRange(OcctKernel&                     kernel,
                                 const char*                     op,
                                 const TopoDS_Shape&             shape,
                                 const TopoDS_Shape&             made,
                                 const TopoDS_Edge&              edge,
                                 double                          start,
                                 double                          end,
                                 const std::vector<TopoDS_Edge>& others,
                                 FaceHistory&                    history,
                                 std::vector<blending::Join>*     joins = nullptr,
                                 const blending::Blended*         provenance = nullptr)
{
  BRepAdaptor_Curve curve(edge);
  // Isolate both sides of the change. Concave blends add material; convex
  // blends remove it. Mixed contours can do both in one operation.
  // Use the explicit boundary provenance where the blend owns it and the
  // signed-shell ownership is proven. Otherwise, and for the stock builder,
  // the two differences share one safe intersection data structure.
  TopoDS_Shape removedShape, addedShape;
  std::vector<blending::Join> removedIntent, addedIntent;
  if (provenance && !provenance->unchanged.empty() && oneShell(shape) && oneShell(made))
  {
    const auto delta = localDelta(op, shape, *provenance);
    removedShape = delta.removed;
    addedShape = delta.added;
    removedIntent = delta.joins;
    addedIntent = delta.joins;
  }
  else
  {
    NCollection_List<TopoDS_Shape> arguments;
    arguments.Append(shape);
    arguments.Append(made);
    BOPAlgo_PaveFiller intersections;
    intersections.SetArguments(arguments);
    intersections.SetNonDestructive(true);
    intersections.Perform();
    BRepAlgoAPI_Cut removed(shape, made, intersections);
    BRepAlgoAPI_Cut added(made, shape, intersections);
    if (!removed.IsDone() || !added.IsDone())
      fail(op, no_solution, "could not isolate fillet volume");
    removedShape = removed.Shape();
    addedShape = added.Shape();
    if (joins)
    {
      removedIntent = rangeJoins(*joins, removed);
      addedIntent = rangeJoins(*joins, added);
    }
  }
  double a = curve.FirstParameter(), b = curve.LastParameter();
  if (edge.Orientation() == TopAbs_REVERSED)
    std::swap(a, b);
  gp_Pnt p0, p1;
  gp_Vec t0, t1;
  curve.D1(a + (b - a) * start, p0, t0);
  curve.D1(a + (b - a) * end, p1, t1);
  if (t0.Magnitude() < 1e-12 || t1.Magnitude() < 1e-12)
    fail(op, unsupported_geometry, "a limit has a singular tangent");
  if (b < a)
  {
    t0.Reverse();
    t1.Reverse();
  }
  const gp_Dir n0(t0), n1(t1);
  auto         inside = [&](const gp_Pnt& p, double tol) {
    return gp_Vec(p0, p).Dot(gp_Vec(n0)) > tol && gp_Vec(p1, p).Dot(gp_Vec(n1)) < -tol;
  };
  // A pair of normal stopping planes cannot describe a winding interval.
  // Refuse those inputs instead of accidentally rounding a second arc.
  for (int i = 0; i <= 256; ++i)
  {
    const double f = double(i) / 256;
    const auto   p = curve.Value(a + (b - a) * f);
    if (f > start + 1e-8 && f < end - 1e-8)
    {
      if (!inside(p, -BRepExt::tolerance::length))
        fail(op, unsupported_geometry, "range crosses its normal stopping planes");
    }
    else if ((f < start - 1e-8 || f > end + 1e-8) && inside(p, BRepExt::tolerance::length))
      fail(op, unsupported_geometry, "range returns inside its normal stopping planes");
  }
  for (const TopoDS_Edge& onward : others)
  {
    BRepAdaptor_Curve other(onward);
    for (int i = 1; i < 64; ++i)
      if (inside(other.Value(other.FirstParameter()
                             + (other.LastParameter() - other.FirstParameter()) * i / 64),
                 BRepExt::tolerance::length))
        fail(op, unsupported_geometry, "tangent continuation re-enters the limited interval");
  }
  auto clip = [&](TopoDS_Shape tool, std::vector<blending::Join>& intent) {
    for (int i = 0; i < 2; ++i)
    {
      const auto         point   = i ? p1 : p0;
      const auto         normal  = i ? n1 : n0;
      const auto         plane   = BRepBuilderAPI_MakeFace(gp_Pln(point, normal)).Face();
      const auto         witness = point.Translated(gp_Vec(normal) * (i ? -1 : 1));
      const auto         half    = BRepPrimAPI_MakeHalfSpace(plane, witness).Solid();
      BRepAlgoAPI_Common common(tool, half);
      if (!common.IsDone())
        fail(op, no_solution, "fillet limit intersection failed");
      if (joins)
      {
        intent = rangeJoins(intent, common);
        // The exact plane supplied to Common owns these new faces. They close
        // the limited material change sharply; they are not support contacts.
        const auto support = BRep_Tool::Surface(plane);
        for (TopExp_Explorer f(common.Shape(), TopAbs_FACE); f.More(); f.Next())
          if (BRep_Tool::Surface(TopoDS::Face(f.Current())) == support)
            for (TopExp_Explorer e(f.Current(), TopAbs_EDGE); e.More(); e.Next())
              intent.push_back({e.Current(), blending::JoinRole::Stop, 0});
      }
      tool = common.Shape();
    }
    return tool;
  };
  Untouched<BRepAlgoAPI_Cut> cut(shape, clip(removedShape, removedIntent));
  if (!cut.IsDone())
    fail(op, no_solution, "limited fillet subtraction failed");
  BRepAlgoAPI_Fuse fuse(cut.Shape(), clip(addedShape, addedIntent));
  if (!fuse.IsDone() || !BRepExt::invalidity(fuse.Shape(), shape).empty())
    fail(op, no_solution, "limited fillet is invalid");
  history.collectComposed(cut, fuse, shape);
  if (joins)
  {
    *joins = rangeJoins(removedIntent, cut);
    joins->insert(joins->end(), addedIntent.begin(), addedIntent.end());
    *joins = rangeJoins(*joins, fuse);
  }
  return history.finish(store(kernel, fuse.Shape()));
}

EvolutionData filletRange(OcctKernel&                kernel,
                          uint32_t                   shapeId,
                          uint32_t                   edgeId,
                          double                     radius,
                          double                     start,
                          double                     end,
                          bool                       propagate,
                          const std::vector<int>&    hashes,
                          int                        bound,
                          const std::vector<double>& law)
{
  const char* op       = "filletRange";
  const auto& shape    = get(kernel, shapeId);
  const auto  selected = subShapes(kernel, op, shape, {edgeId}, TopAbs_EDGE);
  const auto  edge     = TopoDS::Edge(selected[0]);
  if (!std::isfinite(radius) || radius <= 0 || !std::isfinite(start) || !std::isfinite(end)
      || start < 0 || end > 1 || start >= end)
    fail(op, invalid_input, "positive radius and ordered range within [0,1] required");
  BRepAdaptor_Curve rangeCurve(edge);
  const double      rangeLength = GCPnts_AbscissaPoint::Length(
    rangeCurve,
    rangeCurve.FirstParameter()
      + (rangeCurve.LastParameter() - rangeCurve.FirstParameter()) * start,
    rangeCurve.FirstParameter() + (rangeCurve.LastParameter() - rangeCurve.FirstParameter()) * end);
  if (rangeLength <= BRepExt::tolerance::length)
    fail(op, invalid_input, "empty range: limits are within kernel resolution");
  BRepFilletAPI_MakeFillet maker(shape);
  maker.Add(radius, edge);
  const int contour = maker.Contour(edge);
  if (!law.empty())
  {
    if (law.size() < 4 || law.size() % 2 || law.front() != 0 || law[law.size() - 2] != 1)
      fail(op, invalid_input, "radius points must run from 0 to 1");
    const int size  = int(law.size() / 2);
    int       index = 1;
    while (index <= maker.NbEdges(contour) && !maker.Edge(contour, index).IsSame(edge))
      ++index;
    if (index > maker.NbEdges(contour))
      fail(op, no_solution, "selected edge is absent from contour");
    const bool reverse = maker.Edge(contour, index).Orientation() != edge.Orientation();
    NCollection_Array1<gp_Pnt2d> points(1, size);
    for (int i = 0; i < size; ++i)
    {
      if (!std::isfinite(law[2 * i]) || !std::isfinite(law[2 * i + 1]) || law[2 * i + 1] <= 0
          || (i && law[2 * i] <= law[2 * i - 2]))
        fail(op, invalid_input, "ordered finite parameters and positive radii required");
      points(reverse ? size - i : i + 1) =
        gp_Pnt2d(reverse ? 1 - law[2 * i] : law[2 * i], law[2 * i + 1]);
    }
    maker.SetRadius(points, contour, index);
  }
  const bool limited = start > 0 || end < 1 || (!propagate && maker.NbEdges(contour) > 1);
  maker.Build();
  if (!maker.IsDone() || !BRepExt::invalidity(maker.Shape(), shape).empty())
    fail(op, no_solution, "fillet construction failed");
  FaceHistory history(op, hashes, bound);
  if (!limited)
  {
    history.collect(maker, shape);
    return history.finish(store(kernel, maker.Shape()));
  }
  std::vector<TopoDS_Edge> others;
  for (int e = 1; e <= maker.NbEdges(contour); ++e)
    if (!maker.Edge(contour, e).IsSame(edge))
      others.push_back(maker.Edge(contour, e));
  return withinRange(kernel, op, shape, maker.Shape(), edge, start, end, others, history);
}

EvolutionData blendRange(OcctKernel&                kernel,
                         uint32_t                   shapeId,
                         uint32_t                   edgeId,
                         double                     first,
                         double                     second,
                         int                        law,
                         double                     rho,
                         uint32_t                   reference,
                         double                     angle,
                         double                     start,
                         double                     end,
                         const std::vector<int>&    hashes,
                         int                        bound,
                         const std::vector<double>& sizes,
                         std::vector<blending::Join>* joins)
{
  const char* op    = "filletRange";
  const auto& shape = get(kernel, shapeId);
  const auto  edge  = TopoDS::Edge(subShapes(kernel, op, shape, {edgeId}, TopAbs_EDGE)[0]);
  if (!std::isfinite(start) || !std::isfinite(end) || start < 0 || end > 1 || start >= end)
    fail(op, invalid_input, "an ordered range within [0,1] required");
  if (law < 0 || law > 2)
    fail(op, invalid_input, "a range takes a conic, G2 or chamfer section");
  if (sizes.size() % 2)
    fail(op, invalid_input, "size places come in (place, size) pairs");
  // The edge alone, no tangent edges added: the range lies on it.
  blending::Request request;
  request.operation = op;
  request.body      = shape;
  request.edges     = {edge};
  request.mayAdd    = {edge};
  request.first     = first;
  request.second    = second;
  request.law       = static_cast<blending::Law>(law);
  request.rho       = rho;
  request.angle     = angle;
  if (reference)
    request.reference = subShapes(kernel, op, shape, {reference}, TopAbs_FACE)[0];
  if (!sizes.empty())
  {
    auto& points = request.radii.emplace_back();
    for (size_t k = 0; k + 1 < sizes.size(); k += 2)
      points.emplace_back(sizes[k], sizes[k + 1]);
  }
  blending::Blended blended;
  try
  {
    blended = blending::blend(request);
  }
  catch (const Standard_Failure& error)
  {
    fail(op,
         no_solution,
         std::string("OpenCascade could not build the blend (") + error.ExceptionType() + ")");
  }
  double size = std::max(first, second);
  for (const auto& points : request.radii)
    for (const auto& [place, at] : points)
      size = std::max(size, std::max(first, second) * at / first);
  blending::accepted(op, shape, blended.shape, size, &blended);
  if (joins)
    *joins = blended.joins;
  FaceHistory history(op, hashes, bound);
  if (start == 0 && end == 1)
  {
    history.collect(*blended.history, shape);
    return history.finish(store(kernel, blended.shape));
  }
  return withinRange(kernel, op, shape, blended.shape, edge, start, end, {}, history, joins, &blended);
}
} // namespace BRepExt
