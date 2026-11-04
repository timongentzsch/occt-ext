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

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepPatchExt_Inside.hxx>
#include <BRepSurfaceExt_Tools.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <Geom2d_Curve.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <algorithm>
#include <cmath>
#include <cstdio>

namespace BRepExt::surfacing
{
TopoDS_Face faceOf(OcctKernel& k, const char* op, uint32_t id)
{
  if (get(k, id).ShapeType() != TopAbs_FACE)
    fail(op, invalid_input, "a face is required");
  return TopoDS::Face(get(k, id));
}

Boundary boundaryOf(OcctKernel& k, const char* op, uint32_t edgeId, uint32_t faceId)
{
  if (get(k, edgeId).ShapeType() != TopAbs_EDGE)
    fail(op, invalid_input, "an edge is required");
  auto face = faceOf(k, op, faceId);
  for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next())
    if (it.Current().IsSame(get(k, edgeId)))
    {
      auto   edge = TopoDS::Edge(it.Current());
      double a, b;
      auto   curve = BRep_Tool::Curve(edge, a, b);
      if (curve.IsNull() || BRep_Tool::Degenerated(edge))
        fail(op, unsupported_geometry, "a regular edge with a curve in space is required");
      return {edge, face, curve, a, b};
    }
  fail(op, invalid_input, "the edge must belong to its support face");
}

gp_Pnt Boundary::point(double f) const
{
  return curve->Value(first + f * (last - first));
}

SurfaceJet Boundary::jet(const char* op, double f) const
{
  return surfaceJet(op, edge, face, first + f * (last - first));
}

gp_XY boundaryUV(const char* op, const TopoDS_Face& face, const TopoDS_Edge& edge, bool end)
{
  double     a, b;
  const auto pc = BRep_Tool::CurveOnSurface(edge, face, a, b);
  if (pc.IsNull())
    fail(op, invalid_input, "the edge has no parameter curve on the face");
  return pc->Value(end ? b : a).XY();
}

int boundarySide(const char* op, const TopoDS_Face& face, const TopoDS_Edge& edge)
{
  double u0, u1, v0, v1;
  BRepTools::UVBounds(face, u0, u1, v0, v1);
  double a, b;
  auto   pc = BRep_Tool::CurveOnSurface(edge, face, a, b);
  if (pc.IsNull())
    fail(op, invalid_input, "the edge must belong to the face");
  int          side  = -1;
  const double slack = BRepExt::tolerance::parameterMerge * std::max(u1 - u0, v1 - v0);
  for (int s = 0; s < 4; ++s)
  {
    bool matches = true;
    for (int i = 0; i <= 16; ++i)
    {
      auto         uv  = pc->Value(a + (b - a) * i / 16);
      const double gap = s == 0   ? uv.X() - u0
                         : s == 1 ? uv.X() - u1
                         : s == 2 ? uv.Y() - v0
                                  : uv.Y() - v1;
      matches          = matches && std::abs(gap) <= slack;
    }
    if (matches)
      side = s;
  }
  auto         start = pc->Value(a), end = pc->Value(b);
  const double run = side < 2 ? std::abs(end.Y() - start.Y()) : std::abs(end.X() - start.X());
  if (side < 0 || std::abs(run - (side < 2 ? v1 - v0 : u1 - u0)) > slack)
    fail(op, unsupported_geometry, "requires a complete isoparametric boundary edge");
  return side;
}

occ::handle<Geom_BSplineSurface> editable(const char* op, const TopoDS_Face& face, bool fitting)
{
  auto s = occ::down_cast<Geom_BSplineSurface>(BRep_Tool::Surface(face));
  if (s.IsNull())
    fail(op, unsupported_geometry, "requires a B-spline face; rebuild the surface first");
  if (fitting && (s->IsURational() || s->IsVRational() || s->UDegree() > 5 || s->VDegree() > 5))
    fail(op,
         unsupported_geometry,
         "the Patch fit requires non-rational surfaces of degree at most five");
  if (s->IsUPeriodic() || s->IsVPeriodic())
    fail(op, unsupported_geometry, "periodic surfaces have no editable boundary rows");
  double u0, u1, v0, v1, a, b, c, d;
  BRepTools::UVBounds(face, u0, u1, v0, v1);
  s->Bounds(a, b, c, d);
  const double slack = BRepExt::tolerance::parameterMerge * std::max(b - a, d - c);
  if (std::abs(a - u0) + std::abs(b - u1) + std::abs(c - v0) + std::abs(d - v1) > slack)
    fail(op, unsupported_geometry, "requires an untrimmed rectangular B-spline face");
  int count = 0, mask = 0;
  for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next())
  {
    ++count;
    mask |= 1 << boundarySide(op, face, TopoDS::Edge(it.Current()));
  }
  if (count != 4 || mask != 15)
    fail(op, unsupported_geometry, "trimmed faces and holes are unsupported");
  return occ::down_cast<Geom_BSplineSurface>(s->Copy());
}

SplineNet netOf(const occ::handle<Geom_BSplineSurface>& s)
{
  SplineNet net;
  net.degreeU = s->UDegree();
  net.degreeV = s->VDegree();
  auto ku = s->UKnotSequence(), kv = s->VKnotSequence();
  net.knotsU.assign(ku.begin(), ku.end());
  net.knotsV.assign(kv.begin(), kv.end());
  for (int i = 1; i <= s->NbUPoles(); ++i)
    for (int j = 1; j <= s->NbVPoles(); ++j)
      net.poles.push_back(s->Pole(i, j).XYZ());
  return net;
}

Constraint constraint(const char*        op,
                      const TopoDS_Face& face,
                      const TopoDS_Edge& edge,
                      Boundary           ref,
                      int                order)
{
  if (order < 0 || order > 2)
    fail(op, invalid_input, "continuity must be G0, G1 or G2");
  boundarySide(op, face, edge);
  auto       a = boundaryUV(op, face, edge, false), b = boundaryUV(op, face, edge, true);
  auto       s = BRep_Tool::Surface(face);
  const auto p = s->Value(a.X(), a.Y()), q = s->Value(b.X(), b.Y());
  if (p.Distance(ref.point(1)) + q.Distance(ref.point(0))
      < p.Distance(ref.point(0)) + q.Distance(ref.point(1)))
    std::swap(ref.first, ref.last);
  return {ref, a, b, order};
}

namespace
{
// Derivative constraints only. Values are equation tolerances, scaled by the current metric.
constexpr double normalWeight = 1e-3, normalWithin = 1e-7;
constexpr double curvatureWeight = 1e-6, curvatureWithin = 1e-9, metricFloor = 1;

void appendJetConditions(std::vector<FitCondition>& out,
                         const gp_XY&               uv,
                         const SurfaceJet&          j,
                         const SurfaceJet&          current,
                         int                        order)
{
  if (order == 0)
    return;
  const auto n = j.n.XYZ();
  for (const auto [du, dv] : {std::pair{1, 0}, std::pair{0, 1}})
    out.push_back({uv.X(),
                   uv.Y(),
                   {{du, dv, n}},
                   1,
                   0,
                   normalWeight,
                   false,
                   normalWithin * (du ? current.u.Magnitude() : current.v.Magnitude())});
  if (order < 2)
    return;
  // The same linearisation as Patch, applied to all entries of the second form.
  const auto e1 = j.u.Normalized(), e2 = j.n.Crossed(e1);
  const auto w = [&](gp_Vec t) {
    return (e1 * j.curvature(t, e1) + e2 * j.curvature(t, e2)).XYZ();
  };
  for (const auto [du, dv] : {std::pair{2, 0}, std::pair{1, 1}, std::pair{0, 2}})
  {
    auto         x = du == 0 ? current.v : current.u, y = dv == 0 ? current.u : current.v;
    const double value = j.curvature(x, y);
    out.push_back({uv.X(),
                   uv.Y(),
                   {{du, dv, n},
                    {1, 0, du == 0 ? gp_XYZ(0, 0, 0) : w(y) * -1},
                    {0, 1, dv == 0 ? gp_XYZ(0, 0, 0) : w(x) * -1},
                    {du == 2 ? 1 : 0,
                     du == 2 ? 0 : 1,
                     du == 2   ? w(x) * -1
                     : dv == 2 ? w(y) * -1
                               : gp_XYZ(0, 0, 0)}},
                   4,
                   -value,
                   curvatureWeight,
                   false,
                   curvatureWithin * std::max(metricFloor, x.Magnitude() * y.Magnitude())});
  }
}
} // namespace

std::vector<FitCondition> conditions(const char*                             op,
                                     const std::vector<Constraint>&          constraints,
                                     const occ::handle<Geom_BSplineSurface>& s,
                                     double                                  tolerance,
                                     bool                                    normalFlow)
{
  std::vector<FitCondition> out;
  double                    a, b, c, d;
  s->Bounds(a, b, c, d);
  for (const auto& target : constraints)
  {
    const int count = std::max(129, 8 * std::max(s->NbUPoles(), s->NbVPoles()));
    for (int i = 0; i <= count; ++i)
    {
      const double f  = double(i) / count;
      const auto   uv = target.start + (target.end - target.start) * f;
      const auto   p  = target.reference.point(f);
      for (int axis = 1; axis <= 3; ++axis)
      {
        gp_XYZ direction(0, 0, 0);
        direction.SetCoord(axis, 1);
        out.push_back(
          {uv.X(), uv.Y(), {{0, 0, direction}}, 1, p.Coord(axis), 1, false, tolerance * 0.05});
      }
      if (normalFlow || target.skewed)
      {
        gp_Pnt p;
        gp_Vec along;
        target.reference.curve->D1(target.reference.first
                                     + f * (target.reference.last - target.reference.first),
                                   p,
                                   along);
        const bool acrossU =
          std::abs(target.end.X() - target.start.X()) < std::abs(target.end.Y() - target.start.Y());
        const auto   jet = surfaceJet(op, s, uv.X(), uv.Y());
        const gp_Vec towards =
          (target.end - target.start).X() * jet.u + (target.end - target.start).Y() * jet.v;
        gp_Vec running = along.Normalized();
        if (running.Dot(towards) < 0)
          running.Reverse();
        // The way off the boundary is held square to this: `along` itself where it goes
        // straight across, turned by the skew where it leans.
        gp_Vec square = running;
        if (target.skewed)
        {
          const gp_Vec across = acrossU ? jet.u : jet.v;
          gp_Vec       inward = jet.n.Crossed(running);
          if (inward.Dot(across) < 0)
            inward.Reverse();
          square = running - inward * target.skewAt(f);
          square.Normalize();
        }
        out.push_back({uv.X(),
                       uv.Y(),
                       {{acrossU ? 1 : 0, acrossU ? 0 : 1, square.XYZ()}},
                       1,
                       0,
                       1e-3,
                       false,
                       1e-7 * (acrossU ? jet.u.Magnitude() : jet.v.Magnitude())});
      }
      if (target.order == 0)
        continue;
      const auto j = target.reference.jet(op, f), current = surfaceJet(op, s, uv.X(), uv.Y());
      appendJetConditions(out, uv, j, current, target.order);
    }
  }
  return out;
}

TopoDS_Face makeFace(const char*                             op,
                     const occ::handle<Geom_BSplineSurface>& s,
                     double                                  tolerance,
                     TopAbs_Orientation                      orientation)
{
  if (!std::isfinite(tolerance) || tolerance <= 0)
    fail(op, invalid_input, "tolerance must be positive and finite");
  BRepBuilderAPI_MakeFace builder(s, Precision::Confusion());
  if (!builder.IsDone() || !BRepCheck_Analyzer(builder.Face()).IsValid())
    fail(op, no_solution, "the result is not a valid face");
  auto face = builder.Face();
  face.Orientation(orientation);
  if (auto fault = insideFault(face, 0))
    fail(op, no_solution, *fault);
  return face;
}

void measure(const char*                    op,
             const TopoDS_Face&             face,
             const std::vector<Constraint>& constraints,
             double                         tolerance,
             bool                           normalFlow)
{
  const auto s = BRep_Tool::Surface(face);
  for (const auto& target : constraints)
  {
    double gap = 0, angle = 0, bend = 0, relative = 0, skew = 0;
    for (int i = 0; i <= 512; ++i)
    {
      const double f  = double(i) / 512;
      auto         uv = target.start + (target.end - target.start) * f;
      auto         at = s->Value(uv.X(), uv.Y());
      gap             = std::max(gap, at.Distance(target.reference.point(f)));
      if (!target.reference.face.IsNull())
      {
        double a, b;
        auto   pc = BRep_Tool::CurveOnSurface(target.reference.edge, target.reference.face, a, b);
        if (pc.IsNull())
          fail(op, invalid_input, "the reference edge has no parameter curve");
        auto on =
          pc->Value(target.reference.first + f * (target.reference.last - target.reference.first));
        gap =
          std::max(gap,
                   at.Distance(BRep_Tool::Surface(target.reference.face)->Value(on.X(), on.Y())));
      }
      if (target.skewed)
      {
        gp_Pnt on;
        gp_Vec along;
        target.reference.curve->D1(target.reference.first
                                     + f * (target.reference.last - target.reference.first),
                                   on,
                                   along);
        const auto   jet = surfaceJet(op, s, uv.X(), uv.Y());
        const gp_Vec towards =
          (target.end - target.start).X() * jet.u + (target.end - target.start).Y() * jet.v;
        if (along.Dot(towards) < 0)
          along.Reverse();
        along.Normalize();
        const bool acrossU =
          std::abs(target.end.X() - target.start.X()) < std::abs(target.end.Y() - target.start.Y());
        const gp_Vec across = acrossU ? jet.u : jet.v;
        gp_Vec       inward = jet.n.Crossed(along);
        if (inward.Dot(across) < 0)
          inward.Reverse();
        skew = std::max(
          skew,
          std::abs(std::atan2(across.Dot(along), across.Dot(inward)) - std::atan(target.skewAt(f)))
            * 180 / M_PI);
      }
      else if (normalFlow)
      {
        const auto jet = surfaceJet(op, s, uv.X(), uv.Y());
        skew           = std::max(
          skew,
          std::asin(std::clamp(std::abs(jet.u.Normalized().Dot(jet.v.Normalized())), 0.0, 1.0))
            * 180 / M_PI);
      }
      if (target.order < 1)
        continue;
      const auto   j = surfaceJet(op, s, uv.X(), uv.Y()), r = target.reference.jet(op, f);
      const double cosine = j.n.Dot(r.n);
      angle = std::max(angle, std::acos(std::clamp(std::abs(cosine), 0.0, 1.0)) * 180 / M_PI);
      gp_Pnt p;
      gp_Vec tangent;
      target.reference.curve->D1(target.reference.first
                                   + f * (target.reference.last - target.reference.first),
                                 p,
                                 tangent);
      const auto [absolute, share] =
        bendDifference(j, r, tangent.Normalized(), cosine < 0 ? -1 : 1);
      bend     = std::max(bend, absolute);
      relative = std::max(relative, share);
    }
    if ((normalFlow || target.skewed) && skew > .1)
    {
      char message[160];
      std::snprintf(message,
                    sizeof message,
                    "boundary UV flow misses %s by %.9g degrees (limit 0.1)",
                    target.skewed ? "its skew" : "perpendicular",
                    skew);
      fail(op, no_solution, message);
    }
    if (gap > tolerance || (target.order > 0 && angle > 0.1)
        || (target.order > 1 && relative > 0.05 && bend > 1e-7))
    {
      char message[256];
      std::snprintf(
        message,
        sizeof message,
        "boundary misses: gap %.9g (limit %.9g), angle %.9g degrees (0.1), curvature %.9g (share "
        "%.9g, limit 0.05); the available poles cannot meet the requested continuity",
        gap,
        tolerance,
        angle,
        bend,
        relative);
      fail(op, no_solution, message);
    }
  }
}
} // namespace BRepExt::surfacing
