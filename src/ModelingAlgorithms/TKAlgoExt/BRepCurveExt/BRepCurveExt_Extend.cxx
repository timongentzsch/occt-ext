// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>

#include <BRepAlgoAPI_Section.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Tool.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomConvert.hxx>
#include <GeomConvert_BSplineCurveToBezierCurve.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Line.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <BRepBndLib.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <Bnd_Box.hxx>
#include <NCollection_Array1.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Wire.hxx>
#include <algorithm>
#include <cmath>
#include <gp_Ax2.hxx>

namespace BRepExt
{
namespace
{
// A copy in walking direction, bounded to the edge rather than its underlying
// spline.
Handle(Geom_TrimmedCurve) directed(const TopoDS_Edge& edge)
{
  double     a, b;
  const auto basis = BRep_Tool::Curve(edge, a, b);
  if (basis.IsNull())
    fail("extendCurve", invalid_input, "edge has no 3D curve");
  Handle(Geom_TrimmedCurve) curve = new Geom_TrimmedCurve(basis, a, b);
  if (edge.Orientation() == TopAbs_REVERSED)
    curve->Reverse();
  return curve;
}

double length(const Handle(Geom_Curve)& curve, double a, double b, double tolerance)
{
  GeomAdaptor_Curve adaptor(curve);
  return GCPnts_AbscissaPoint::Length(adaptor, a, b, tolerance);
}

// One extension by a distance; with `reach` a shape that cannot be that long is as long as it
// can be.
EvolutionData extendBy(OcctKernel& kernel,
                       uint32_t    wireId,
                       bool        atStart,
                       int         mode,
                       double      distance,
                       uint32_t    targetId,
                       double      tolerance,
                       bool        reach)
{
  const char* const op    = "extendCurve";
  const auto&       shape = get(kernel, wireId);
  if (shape.ShapeType() != TopAbs_WIRE || mode < 0 || mode > 4 || !std::isfinite(distance)
      || distance <= tolerance || !std::isfinite(tolerance) || tolerance <= 0)
    fail(op, invalid_input, "an open wire, positive distance and tolerance are required");
  std::vector<TopoDS_Edge> edges;
  for (BRepTools_WireExplorer it(TopoDS::Wire(shape)); it.More(); it.Next())
    edges.push_back(it.Current());
  if (edges.empty())
    fail(op, invalid_input, "wire has no edges");
  auto begin = directed(edges.front()), finish = directed(edges.back());
  if (begin->StartPoint().Distance(finish->EndPoint()) <= tolerance)
    fail(op, invalid_input, "a closed curve has no endpoint to extend");
  auto source = atStart ? begin : finish;
  if (atStart)
    source->Reverse();
  const double end = source->LastParameter();
  gp_Pnt       p;
  gp_Vec       d1, d2;
  source->D2(end, p, d1, d2);
  if (d1.Magnitude() <= 1e-15)
    fail(op, invalid_input, "endpoint tangent is singular");
  const gp_Vec       tangent   = d1.Normalized();
  const gp_Vec       curvature = (d2 - tangent * d2.Dot(tangent)) / d1.SquareMagnitude();
  Handle(Geom_Curve) extension;
  double             first = 0, last = 1;
  if (mode == 1 || ((mode == 2 || mode == 4) && curvature.Magnitude() < 1e-14))
  {
    extension = new Geom_Line(p, gp_Dir(tangent));
    last      = distance;
  }
  else if (mode == 4)
  {
    const double radius = 1 / curvature.Magnitude();
    // The circle of an arc's own curvature is the arc's circle: what the arc already
    // takes of it is not free for the extension.
    const double taken = source->BasisCurve()->IsKind(STANDARD_TYPE(Geom_Circle))
                           ? length(source, source->FirstParameter(), end, tolerance)
                           : 0;
    if (reach)
      distance = std::min(distance, (2 * M_PI * radius - taken) * 0.999);
    if (distance + taken >= 2 * M_PI * radius - tolerance)
      fail(op,
           invalid_input,
           taken > 0 ? "an arc extension would close or overlap the arc's own circle"
                     : "an arc extension must be shorter than a full circle");
    const gp_Pnt center = p.Translated(curvature * (radius * radius));
    extension =
      new Geom_Circle(gp_Ax2(center, gp_Dir(tangent.Crossed(curvature)), gp_Dir(gp_Vec(center, p))),
                      radius);
    last = distance / radius;
  }
  else if (mode == 3)
  {
    const double available = length(source, source->FirstParameter(), end, tolerance);
    if (reach)
      distance = std::min(distance, available);
    if (distance > available + tolerance)
      fail(op, invalid_input, "reflective distance exceeds the endpoint segment's length");
    GeomAdaptor_Curve    adaptor(source);
    GCPnts_AbscissaPoint cut(tolerance, adaptor, -std::min(distance, available), end);
    if (!cut.IsDone())
      fail(op, no_solution, "could not locate reflective length");
    Handle(Geom_TrimmedCurve) reflected = new Geom_TrimmedCurve(source, cut.Parameter(), end);
    reflected->Reverse();
    reflected->Mirror(gp_Ax2(p, gp_Dir(tangent)));
    extension = reflected;
    first     = reflected->FirstParameter();
    last      = reflected->LastParameter();
  }
  else if (mode == 0 && source->BasisCurve()->IsKind(STANDARD_TYPE(Geom_Line)))
  {
    extension = new Geom_Line(p, gp_Dir(tangent));
    last      = distance;
  }
  else if (mode == 0 && source->BasisCurve()->IsKind(STANDARD_TYPE(Geom_Circle)))
  {
    auto circle = Handle(Geom_Circle)::DownCast(source->BasisCurve());
    if (reach)
      distance = std::min(distance,
                          0.999
                            * (2 * M_PI * circle->Radius()
                               - length(source, source->FirstParameter(), end, tolerance)));
    if (distance + length(source, source->FirstParameter(), end, tolerance)
        >= 2 * M_PI * circle->Radius() - tolerance)
      fail(op, invalid_input, "natural extension would close or overlap the circle");
    extension = circle;
    first     = end;
    last      = end + distance / circle->Radius();
  }
  else if (mode == 0 && GeomAdaptor_Curve(source->BasisCurve()).GetType() != GeomAbs_BSplineCurve
           && GeomAdaptor_Curve(source->BasisCurve()).GetType() != GeomAbs_BezierCurve)
  {
    // Keep analytic conics analytic; a rational Bezier continuation has a
    // finite pole.
    extension = source->BasisCurve();
    GeomAdaptor_Curve adaptor(extension);
    if (extension->IsPeriodic()
        && distance + length(source, source->FirstParameter(), end, tolerance)
             >= length(extension, 0, extension->Period(), tolerance) - tolerance)
      fail(op, invalid_input, "natural extension would close or overlap the conic");
    GCPnts_AbscissaPoint cut(adaptor, distance, end, end + distance / d1.Magnitude(), tolerance);
    if (!cut.IsDone())
      fail(op, no_solution, "could not locate analytic extension length");
    first = end;
    last  = cut.Parameter();
  }
  else
  {
    Handle(Geom_BezierCurve) terminal;
    if (mode == 0)
    {
      GeomConvert_BSplineCurveToBezierCurve arcs(GeomConvert::CurveToBSplineCurve(source));
      terminal = arcs.Arc(arcs.NbArcs());
    }
    // Natural continues the last rational polynomial span exactly. Soft
    // preserves G2 at the join and reduces curvature to zero at its far end.
    // Solve the actual arc length.
    auto candidate = [&](double extent) -> Handle(Geom_BezierCurve) {
      if (mode == 0)
      {
        auto copy = Handle(Geom_BezierCurve)::DownCast(terminal->Copy());
        copy->Segment(1, 1 + extent);
        return copy;
      }
      NCollection_Array1<gp_Pnt> poles(1, 4);
      poles(1) = p;
      poles(2) = p.Translated(tangent * extent);
      poles(3) = p.Translated(tangent * (2 * extent) + curvature * (1.5 * extent * extent));
      poles(4) = p.Translated(tangent * (3 * extent) + curvature * (3 * extent * extent));
      return new Geom_BezierCurve(poles);
    };
    double lo = 0, hi = mode == 0 ? 1 : distance / 3;
    for (int i = 0; length(candidate(hi), 0, 1, tolerance / 10) < distance; ++i)
    {
      if (i == 40)
        fail(op, no_solution, "extension length did not converge");
      hi *= 2;
    }
    for (int i = 0; i < 60; ++i)
    {
      const double mid = (lo + hi) / 2;
      if (length(candidate(mid), 0, 1, tolerance / 10) < distance)
        lo = mid;
      else
        hi = mid;
    }
    extension = candidate((lo + hi) / 2);
  }
  TopoDS_Edge added = BRepBuilderAPI_MakeEdge(extension, first, last);
  if (targetId != 0)
  {
    BRepAlgoAPI_Section section(added, get(kernel, targetId), false);
    section.Approximation(true);
    section.Build();
    if (!section.IsDone())
      fail(op, no_solution, "target intersection failed");
    double nearest = last + 1;
    for (TopExp_Explorer it(section.Shape(), TopAbs_VERTEX); it.More(); it.Next())
    {
      GeomAPI_ProjectPointOnCurve projection(BRep_Tool::Pnt(TopoDS::Vertex(it.Current())),
                                             extension,
                                             first,
                                             last);
      for (int i = 1; i <= projection.NbPoints(); ++i)
      {
        const double u = projection.Parameter(i);
        if (projection.Distance(i) <= tolerance
            && length(extension, first, u, tolerance) > tolerance)
          nearest = std::min(nearest, u);
      }
    }
    if (nearest > last)
      fail(op, no_solution, "no target intersection within the extension distance");
    added = BRepBuilderAPI_MakeEdge(extension, first, nearest);
  }
  // Reuse the wire endpoint itself: MakeWire otherwise replaces the first
  // original edge when merging its vertex with the extension added before it.
  const auto joint =
    atStart ? TopExp::FirstVertex(edges.front(), true) : TopExp::LastVertex(edges.back(), true);
  double     ea, eb;
  const auto ec = BRep_Tool::Curve(added, ea, eb);
  if (ec.IsNull())
    fail(op, no_solution, "extension has no 3D curve");
  added = BRepBuilderAPI_MakeEdge(ec, joint, TopExp::LastVertex(added), ea, eb);
  BRepBuilderAPI_MakeWire result;
  if (atStart)
    result.Add(TopoDS::Edge(added.Reversed()));
  for (const auto& edge : edges)
    result.Add(edge);
  if (!atStart)
    result.Add(added);
  if (!result.IsDone() || !BRepCheck_Analyzer(result.Wire()).IsValid())
    fail(op, no_solution, "extension is not a valid wire");
  return stored(kernel, result.Wire());
}
} // namespace

EvolutionData extendCurve(OcctKernel& kernel,
                          uint32_t    wireId,
                          bool        atStart,
                          int         mode,
                          double      distance,
                          uint32_t    targetId,
                          double      tolerance)
{
  // Up to a target no distance need be given: the extension is tried at an eighth of the
  // curve's own size and doubled until it meets the target, up to twice the size of curve
  // and target together. (One long extension at once is no good: a spline carried on
  // naturally for many times its length runs wild, and a soft one bends by its length.)
  if (targetId == 0 || !std::isfinite(distance) || distance > tolerance || !(tolerance > 0))
    return extendBy(kernel, wireId, atStart, mode, distance, targetId, tolerance, false);
  Bnd_Box own, box;
  BRepBndLib::Add(get(kernel, wireId), own);
  BRepBndLib::Add(get(kernel, wireId), box);
  BRepBndLib::Add(get(kernel, targetId), box);
  const double size = std::sqrt(box.SquareExtent()) + 1;
  // No point of an extension is further from the curve than the extension is long: one
  // shorter than the gap between curve and target cannot meet the target, and is not tried.
  BRepExtrema_DistShapeShape gap(get(kernel, wireId), get(kernel, targetId));
  const double               least = gap.IsDone() && gap.NbSolution() > 0 ? gap.Value() : 0;
  for (double reach = std::max(std::sqrt(own.SquareExtent()) / 8, 2 * tolerance);; reach *= 2)
  {
    if (reach < least && reach < 2 * size)
      continue;
    try
    {
      return extendBy(kernel, wireId, atStart, mode, reach, targetId, tolerance, true);
    }
    catch (const Failure&)
    {
      if (reach >= 2 * size)
        throw;
    }
  }
}
} // namespace BRepExt
