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

// filletVariable: fillets whose radius follows a law along the edge.
//
// BRepFilletAPI_MakeFillet::SetRadius(UandR, contour, edge) interpolates the radius through
// (parameter, radius) points along one edge of a contour, the chain of tangent edges OCCT
// fillets as one. One maker takes all edges, so fillets that meet at a vertex get their corner
// built together.

#include <BRepToolsExt_History.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepFilletExt_BlendChecks.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRep_Tool.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <BRepTools.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <GC_MakeCircle.hxx>
#include <Geom_Circle.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Torus.hxx>
#include <Law_Function.hxx>
#include <Standard_Failure.hxx>
#include <TopExp.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <BRepToolsExt_Operations.hxx>

#include <BRepCheck_Analyzer.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Map.hxx>
#include <TopoDS.hxx>
#include <gp_Pnt2d.hxx>

#include <cmath>

namespace BRepExt
{
namespace
{
// OCCT's variable fillet uses circular sections in the guide plane. Its
// approximated guide need not be the source edge's exact normal plane: cutting
// the two spring curves there samples different stations (not one radius).
// Measure the actual section, at the law's spine parameter, after Build.
void checkRadiusLaw(const char*                                   op,
                    const TopoDS_Shape&                           before,
                    const TopoDS_Shape&                           result,
                    const TopoDS_Edge&                            source,
                    BRepFilletAPI_MakeFillet&                     maker,
                    double                                        first,
                    double                                        last,
                    const std::vector<std::pair<double, double>>& expected)
{
  ShapeSet generated;
  for (TopExp_Explorer e(before, TopAbs_EDGE); e.More(); e.Next())
    for (const auto& f : maker.Generated(e.Current()))
      if (f.ShapeType() == TopAbs_FACE)
        generated.Add(f);
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    old, adjacent;
  TopExp::MapShapesAndUniqueAncestors(before, TopAbs_EDGE, TopAbs_FACE, old);
  TopExp::MapShapesAndUniqueAncestors(result, TopAbs_EDGE, TopAbs_FACE, adjacent);
  if (!old.Contains(source) || old.FindFromKey(source).Extent() != 2)
    fail(op,
         no_solution,
         "radius law has no pair of support faces",
         {source},
         FailureKind::SupportContact);
  const auto&       supports     = old.FindFromKey(source);
  const TopoDS_Face originals[2] = {TopoDS::Face(supports.First()), TopoDS::Face(supports.Last())};
  std::vector<TopoDS_Face> stripes;
  for (const auto& face : generated)
  {
    bool touches[2] = {false, false};
    for (int i = 1; i <= adjacent.Extent(); ++i)
    {
      const auto& pair = adjacent.FindFromIndex(i);
      if (pair.Extent() != 2)
        continue;
      TopoDS_Face support;
      if (pair.First().IsSame(face))
        support = TopoDS::Face(pair.Last());
      else if (pair.Last().IsSame(face))
        support = TopoDS::Face(pair.First());
      else
        continue;
      for (int side = 0; side < 2; ++side)
      {
        touches[side] |= sameBlendSurface(support, originals[side]);
        for (const auto& changed : maker.Modified(originals[side]))
          touches[side] |= support.IsSame(changed);
      }
    }
    if (touches[0] && touches[1])
      stripes.push_back(TopoDS::Face(face));
  }
  const bool constant = std::all_of(expected.begin(), expected.end(), [&](const auto& sample) {
    return std::abs(sample.second - expected.front().second) < tolerance::length;
  });
  int        checked  = 0;
  for (const auto& [fraction, wanted] : expected)
  {
    bool   located = false;
    double actual = 0, residual = 0;
    for (const auto& face : stripes)
    {
      BRepAdaptor_Surface adaptor(face, false);
      if (constant && adaptor.GetType() == GeomAbs_Cylinder)
      {
        actual  = adaptor.Cylinder().Radius();
        located = true;
        break;
      }
      if (constant && adaptor.GetType() == GeomAbs_Torus)
      {
        actual  = adaptor.Torus().MinorRadius();
        located = true;
        break;
      }
      double u0, u1, v0, v1;
      BRepTools::UVBounds(face, u0, u1, v0, v1);
      const double v = constant ? v0 + (v1 - v0) * fraction : first + (last - first) * fraction;
      if (v <= v0 + tolerance::parameter || v >= v1 - tolerance::parameter)
        continue;
      const double middle = (u0 + u1) / 2;
      if (BRepClass_FaceClassifier(face, gp_Pnt2d(middle, v), tolerance::length).State()
          == TopAbs_OUT)
        continue;
      const auto a = adaptor.Value(u0, v), b = adaptor.Value(middle, v), c = adaptor.Value(u1, v);
      GC_MakeCircle circle(a, b, c);
      if (!circle.IsDone())
        continue;
      const auto measured = circle.Value()->Circ();
      actual              = measured.Radius();
      residual            = 0;
      for (int j = 0; j <= 8; ++j)
      {
        const auto p = adaptor.Value(u0 + (u1 - u0) * j / 8, v);
        residual     = std::max(residual, std::abs(p.Distance(measured.Location()) - actual));
        residual     = std::max(
          residual,
          std::abs(gp_Vec(measured.Location(), p).Dot(gp_Vec(measured.Axis().Direction()))));
      }
      located = true;
      break;
    }
    // Corner patches may trim the ends; certify ample interior sections.
    if (!located)
      continue;
    const double allowance = std::max(tolerance::sewing, 1e-3 * wanted);
    if (residual > allowance || std::abs(actual - wanted) > allowance)
      fail(op,
           no_solution,
           "radius law mismatch at " + std::to_string(fraction) + ": requested "
             + std::to_string(wanted) + ", measured " + std::to_string(actual)
             + " mm; section residual " + std::to_string(residual),
           {source},
           FailureKind::RadiusLaw);
    ++checked;
  }
  if (checked < 8)
    fail(op,
         no_solution,
         "radius law has only " + std::to_string(checked)
           + " certified interior sections of the edge",
         {source},
         FailureKind::RadiusLaw);
}

} // namespace

EvolutionData filletVariable(OcctKernel&                  kernel,
                             uint32_t                     shapeId,
                             const std::vector<uint32_t>& edgeIds,
                             const std::vector<int>&      lawSizes,
                             const std::vector<double>&   laws,
                             const std::vector<int>&      inputFaceHashes,
                             int                          hashUpperBound)
{
  const char* const   op    = "filletVariable";
  const TopoDS_Shape& shape = get(kernel, shapeId);
  if (shape.ShapeType() != TopAbs_SOLID)
    fail(op, invalid_input, "the shape is not a solid");
  FaceHistory                     history(op, inputFaceHashes, hashUpperBound);
  const std::vector<TopoDS_Shape> edges = subShapes(kernel, op, shape, edgeIds, TopAbs_EDGE);
  if (lawSizes.size() != edges.size())
    fail(op, invalid_input, "one law per edge is required");

  BRepFilletAPI_MakeFillet                               maker(shape);
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> seen;
  size_t                                                 at = 0;
  std::vector<bool>                                      reversed;
  for (size_t i = 0; i < edges.size(); i++)
  {
    const TopoDS_Edge& edge = TopoDS::Edge(edges[i]);
    const std::string  which =
      "edge " + std::to_string(i + 1) + " of " + std::to_string(edges.size());
    if (!seen.Add(edge))
      fail(op, invalid_input, which + " is given twice");
    const int size = lawSizes[i];
    if (size < 2 || at + 2 * static_cast<size_t>(size) > laws.size())
    {
      fail(op,
           invalid_input,
           "the law of " + which + " needs at least two (parameter, radius) points");
    }
    NCollection_Array1<gp_Pnt2d> law(1, size);
    for (int j = 0; j < size; j++, at += 2)
    {
      const double parameter = laws[at], radius = laws[at + 1];
      const bool   ordered = j == 0 ? parameter == 0.0 : parameter > laws[at - 2];
      if (!ordered || parameter > 1.0 || (j == size - 1 && parameter != 1.0))
      {
        fail(op,
             invalid_input,
             "the law of " + which + " must run from parameter 0 to 1 in increasing order");
      }
      if (!std::isfinite(radius) || !(radius > 0.0))
        fail(op, invalid_input, "the law of " + which + " has a radius that is not positive");
      law.SetValue(j + 1, gp_Pnt2d(parameter, radius));
    }

    BRepAdaptor_Curve source(edge);
    const bool        closed =
      BRep_Tool::IsClosed(edge)
      || source.Value(source.FirstParameter()).Distance(source.Value(source.LastParameter()))
           < tolerance::length;
    if (closed && std::abs(law.First().Y() - law.Last().Y()) > tolerance::length)
      fail(op,
           invalid_input,
           "closed edge radius-law endpoints conflict: " + std::to_string(law.First().Y()) + " and "
             + std::to_string(law.Last().Y()) + " mm refer to the same point",
           {edge},
           FailureKind::RadiusEndpoints);

    // A contour is the chain of tangent edges through an edge; an earlier edge may have
    // created the one this edge is on.
    if (maker.Contour(edge) == 0)
      maker.Add(edge);
    const int contour = maker.Contour(edge);
    if (contour == 0)
    {
      fail(op,
           unsupported_geometry,
           which + " cannot carry a fillet (it must be a sharp edge between two faces)");
    }
    int indexInContour = 1;
    while (!maker.Edge(contour, indexInContour).IsSame(edge))
      indexInContour++;
    // The law is given from the start of the edge's curve; the contour may run the other way.
    reversed.push_back(maker.Edge(contour, indexInContour).Orientation() == TopAbs_REVERSED);
    if (reversed.back())
    {
      NCollection_Array1<gp_Pnt2d> flipped(1, size);
      for (int j = 1; j <= size; j++)
      {
        const gp_Pnt2d& point = law(size + 1 - j);
        flipped.SetValue(j, gp_Pnt2d(1.0 - point.X(), point.Y()));
      }
      law = flipped;
    }
    maker.SetRadius(law, contour, indexInContour);
  }
  try
  {
    maker.Build();
  }
  catch (const Standard_Failure& error)
  {
    fail(op,
         no_solution,
         std::string("OCCT variable fillet: ") + error.what(),
         edges,
         FailureKind::NativeBuilder);
  }
  if (!maker.IsDone())
  {
    fail(op,
         no_solution,
         "the fillet could not be built; a radius is probably too large for its neighbourhood");
  }
  std::vector<std::vector<std::pair<double, double>>> expected(edges.size());
  std::vector<std::pair<double, double>>              spineBounds;
  size_t                                              inputAt = 0;
  for (size_t i = 0; i < edges.size(); ++i)
  {
    const auto edge    = TopoDS::Edge(edges[i]);
    const int  contour = maker.Contour(edge);
    double     first = 0, last = 0;
    for (int j = 1; j <= maker.NbEdges(contour); ++j)
    {
      BRepAdaptor_Curve curve(maker.Edge(contour, j));
      last = first + GCPnts_AbscissaPoint::Length(curve);
      if (maker.Edge(contour, j).IsSame(edge))
        break;
      first = last;
    }
    spineBounds.emplace_back(first, last);
    const bool                 constant = maker.IsConstant(contour, edge);
    const Handle(Law_Function) law =
      constant ? Handle(Law_Function)() : maker.GetLaw(contour, edge);
    const auto value = [&](double along) {
      return constant ? maker.Radius(contour, edge) : law->Value(first + (last - first) * along);
    };
    for (int j = 0; j < lawSizes[i]; ++j, inputAt += 2)
    {
      const double fraction = laws[inputAt], wanted = laws[inputAt + 1];
      const double along = reversed[i] ? 1 - fraction : fraction;
      if (std::abs(value(along) - wanted) > tolerance::length)
        fail(op,
             no_solution,
             "the configured radius law does not retain its requested interpolation points",
             {edge},
             FailureKind::RadiusLaw);
    }
    for (int j = 1; j < 18; ++j)
    {
      const double fraction = j / 18.0, along = reversed[i] ? 1 - fraction : fraction;
      expected[i].emplace_back(fraction, value(along));
    }
  }
  // MakeFillet can report success with an open shell where a blend runs off its support.
  // Topology-only check, like upstream's validateFilletResult.
  // The maker wraps its result in a compound; hand back the solid, like upstream's fillet.
  const TopoDS_Shape result = solidIfSingle(maker.Shape());
  if (result.IsNull() || !BRepCheck_Analyzer(result, false).IsValid())
  {
    fail(op, no_solution, "the fillet left the solid open");
  }
  for (size_t i = 0; i < edges.size(); ++i)
  {
    // Surface V follows the contour, while the requested fraction follows the
    // source curve. Reverse the bounds along with that edge's law.
    const auto [first, last] = spineBounds[i];
    checkRadiusLaw(op,
                   shape,
                   result,
                   TopoDS::Edge(edges[i]),
                   maker,
                   reversed[i] ? last : first,
                   reversed[i] ? first : last,
                   expected[i]);
  }
  BooleanResult descendants;
  for (TopExp_Explorer f(shape, TopAbs_FACE); f.More(); f.Next())
  {
    const auto& changed = maker.Modified(f.Current());
    if (!changed.IsEmpty())
      descendants.faces.Bind(f.Current(), changed);
  }
  std::vector<TopoDS_Shape> generated;
  for (const auto& edge : edges)
    for (const auto& made : maker.Generated(edge))
      if (made.ShapeType() == TopAbs_FACE)
        generated.push_back(made);
  auto roles                  = blendBoundaryRoles(shape, edges, result, descendants, generated);
  roles.constructionTolerance = 1e-4; // Same native ChFi3d approximation contract.
  certifyBlendBoundaries(op, result, roles, 1);
  history.collect(maker, shape);
  return history.finish(store(kernel, result));
}

} // namespace BRepExt
