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

// Rounding the corners of a planar wire.
//
// A circle of radius r tangent to the two pieces that meet at a corner has its centre M on the
// offset of each piece by r towards the inside of the turn, and touches the pieces at the feet of M
// on them. With C_i(u) the pieces, T_i their unit tangents and n the normal of the plane, the
// offset is O_i(u) = C_i(u) + sigma r (n x T_i(u)), sigma = +1 for a left turn and -1 for a right
// one, and the centre solves O_1(s) = O_2(t): two equations in the two parameters. Differentiating
// the unit normal gives dO_i/du = C_i' + sigma r (n x dT_i/du), with dT_i/du = (C_i'' - (C_i'' .
// T_i) T_i) / |C_i'|, so the system is solved by Newton's method with exact derivatives, from the
// closed-form line-line guess and from a grid of starts over the two pieces; of the roots whose
// feet are on the pieces, the one nearest the corner is taken (it is the circle that touches first)
// if the circle does not cut what is kept of either piece. For two lines the closed form w = r /
// tan(theta / 2) of the distance from the corner to each foot is the first start and already the
// root. A negative radius cuts the corner along the chord between the same two feet. Corners are
// made in the order the wire runs; a piece that two selected corners share is trimmed at both ends
// in turn. Curved pieces whose curvature jumps (a spline of degree 2) are covered: only two
// derivatives are needed and none is differentiated again.

#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_Failure.hxx>

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepLib_FindSurface.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Tool.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Plane.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Circ.hxx>
#include <gp_Pln.hxx>
#include <gp_Vec.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <string>
#include <vector>

namespace BRepExt
{
namespace
{

// A piece of the wire as a curve whose parameter increases the way the wire runs.
struct Piece
{
  Handle(Geom_Curve) curve;
  double             first = 0, last = 0;
};

Piece pieceOf(const TopoDS_Edge& edge)
{
  TopLoc_Location    location;
  double             first, last;
  Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, location, first, last);
  if (curve.IsNull())
    fail("roundWireCorners", unsupported_geometry, "an edge of the wire has no curve");
  curve = Handle(Geom_Curve)::DownCast(curve->Transformed(location.Transformation()));
  if (edge.Orientation() == TopAbs_REVERSED)
  {
    const double a = curve->ReversedParameter(last), b = curve->ReversedParameter(first);
    curve = curve->Reversed();
    first = a;
    last  = b;
  }
  return {curve, first, last};
}

// The offset of a piece towards the side `sigma` n x T, with its derivative in the parameter.
struct Offset
{
  gp_Pnt point;
  gp_Vec rate;
};

Offset offsetAt(const Piece& piece, double u, const gp_Dir& n, double sigma, double radius)
{
  gp_Pnt p;
  gp_Vec d1, d2;
  piece.curve->D2(u, p, d1, d2);
  const double speed = d1.Magnitude();
  const gp_Vec T     = d1 / speed;
  const gp_Vec dT    = (d2 - T * d2.Dot(T)) / speed;
  const gp_Vec N     = gp_Vec(n).Crossed(T) * sigma;
  return {p.Translated(N * radius), d1 + gp_Vec(n).Crossed(dT) * (sigma * radius)};
}

struct Fillet
{
  double s = 0, t = 0; // parameters of the feet on the first and the second piece
  gp_Pnt centre;
};

// The circle of `radius` tangent to `a` (ending at the corner) and `b` (starting there), or none.
std::optional<Fillet> circleAt(const Piece&  a,
                               const Piece&  b,
                               const gp_Pln& plane,
                               double        radius,
                               double        tolerance)
{
  const gp_Dir n = plane.Axis().Direction();
  const gp_Dir X = plane.XAxis().Direction(), Y = plane.YAxis().Direction();
  gp_Vec       ta, tb;
  {
    gp_Pnt p;
    a.curve->D1(a.last, p, ta);
    b.curve->D1(b.first, p, tb);
  }
  const double speedA = ta.Magnitude(), speedB = tb.Magnitude();
  if (speedA <= 0 || speedB <= 0)
    return std::nullopt;
  ta /= speedA;
  tb /= speedB;
  const double cross = gp_Vec(n).Dot(ta.Crossed(tb));
  const double dot   = ta.Dot(tb);
  const double turn  = std::atan2(std::abs(cross), dot); // the angle the wire turns through
  if (turn < 1e-9)
    return std::nullopt; // already tangent: nothing to round
  if (turn > M_PI - 1e-9)
    return std::nullopt; // a reversal: the circle would be infinitely far
  const double sigma = cross > 0 ? 1.0 : -1.0;

  const auto residual = [&](double s, double t, Offset& oa, Offset& ob) {
    oa = offsetAt(a, s, n, sigma, radius);
    ob = offsetAt(b, t, n, sigma, radius);
    const gp_Vec f(ob.point, oa.point);
    return std::array<double, 2>{f.Dot(gp_Vec(X)), f.Dot(gp_Vec(Y))};
  };
  const auto newton = [&](double s, double t) -> std::optional<Fillet> {
    Offset oa, ob;
    for (int step = 0; step < 60; ++step)
    {
      auto   f    = residual(s, t, oa, ob);
      double size = std::hypot(f[0], f[1]);
      if (size < 1e-13 * std::max(1.0, radius))
        return Fillet{s, t, oa.point};
      // d(oa - ob) = oa.rate ds - ob.rate dt
      const double a11 = oa.rate.Dot(gp_Vec(X)), a21 = oa.rate.Dot(gp_Vec(Y));
      const double a12 = -ob.rate.Dot(gp_Vec(X)), a22 = -ob.rate.Dot(gp_Vec(Y));
      const double det = a11 * a22 - a12 * a21;
      if (std::abs(det) < 1e-14 * (std::hypot(a11, a21) * std::hypot(a12, a22)))
        return std::nullopt;
      double ds = -(a22 * f[0] - a12 * f[1]) / det, dt = -(-a21 * f[0] + a11 * f[1]) / det;
      double lambda = 1;
      for (int halve = 0; halve < 12; ++halve, lambda /= 2)
      {
        Offset     pa, pb;
        const auto g = residual(s + lambda * ds, t + lambda * dt, pa, pb);
        if (std::hypot(g[0], g[1]) < size)
          break;
      }
      s += lambda * ds;
      t += lambda * dt;
      // Out of the curves' own parameters there is no root to find on these pieces.
      const double ra = a.last - a.first, rb = b.last - b.first;
      if (s < a.first - ra || s > a.last + ra || t < b.first - rb || t > b.last + rb)
        return std::nullopt;
    }
    return std::nullopt;
  };

  // Starts: the line-line guess, then a grid over both pieces from the corner outwards.
  const double                           w = radius * std::tan(turn / 2);
  std::vector<std::pair<double, double>> starts;
  starts.emplace_back(a.last - w / speedA, b.first + w / speedB);
  constexpr int grid = 12;
  for (int i = 1; i <= grid; ++i)
    for (int j = 1; j <= grid; ++j)
      starts.emplace_back(a.last - (a.last - a.first) * i / (grid + 1),
                          b.first + (b.last - b.first) * j / (grid + 1));

  const GeomAdaptor_Curve adaptA(a.curve, a.first, a.last), adaptB(b.curve, b.first, b.last);
  const double            slack = 1e-9 * std::max(1.0, a.last - a.first);
  std::optional<Fillet>   best;
  double                  bestReach = HUGE_VAL;
  for (const auto& [s0, t0] : starts)
  {
    const auto root = newton(s0, t0);
    if (!root)
      continue;
    if (root->s < a.first - slack || root->s > a.last + slack || root->t < b.first - slack
        || root->t > b.last + slack)
      continue;
    // The circle must touch from the inside of the turn and not cut what is kept of either piece.
    const double reach = GCPnts_AbscissaPoint::Length(adaptA, root->s, a.last)
                         + GCPnts_AbscissaPoint::Length(adaptB, b.first, root->t);
    if (reach >= bestReach)
      continue;
    bool cuts = false;
    for (int k = 0; k <= 64 && !cuts; ++k)
    {
      const gp_Pnt pa = a.curve->Value(a.first + (root->s - a.first) * k / 64.0);
      const gp_Pnt pb = b.curve->Value(root->t + (b.last - root->t) * k / 64.0);
      cuts            = pa.Distance(root->centre) < radius * (1 - 1e-6) - tolerance
             || pb.Distance(root->centre) < radius * (1 - 1e-6) - tolerance;
    }
    if (cuts)
      continue;
    best      = root;
    bestReach = reach;
  }
  return best;
}

} // namespace

EvolutionData roundWireCorners(OcctKernel&                  kernel,
                               uint32_t                     wireId,
                               const std::vector<uint32_t>& vertexIds,
                               double                       radius,
                               double                       tolerance)
{
  const char* const   op    = "roundWireCorners";
  const TopoDS_Shape& shape = get(kernel, wireId);
  if (shape.ShapeType() != TopAbs_WIRE)
    fail(op, invalid_input, "a wire is required");
  if (!std::isfinite(radius) || radius == 0.0 || !std::isfinite(tolerance) || tolerance <= 0.0)
    fail(op, invalid_input, "a nonzero finite radius and a positive tolerance are required");
  if (vertexIds.empty())
    fail(op, invalid_input, "no vertex was chosen");
  const TopoDS_Wire wire = TopoDS::Wire(shape);

  const std::vector<TopoDS_Shape> chosen = subShapes(kernel, op, wire, vertexIds, TopAbs_VERTEX);
  for (size_t i = 0; i < chosen.size(); ++i)
    for (size_t k = i + 1; k < chosen.size(); ++k)
      if (chosen[i].IsSame(chosen[k]))
        fail(op, invalid_input, "a vertex is given twice");

  // The plane every piece lies in: the circles are in it.
  const BRepLib_FindSurface surface(wire, tolerance, true);
  if (!surface.Found())
    fail(op, unsupported_geometry, "the wire is not planar");
  const gp_Pln plane = Handle(Geom_Plane)::DownCast(surface.Surface())->Pln();

  // The pieces in the order the wire runs, and the vertex each starts at.
  std::vector<Piece>         pieces;
  std::vector<TopoDS_Vertex> starts;
  for (BRepTools_WireExplorer it(wire); it.More(); it.Next())
  {
    pieces.push_back(pieceOf(it.Current()));
    starts.push_back(it.CurrentVertex());
  }
  if (pieces.empty())
    fail(op, invalid_input, "the wire has no edges");
  const size_t n      = pieces.size();
  const bool   closed = BRep_Tool::IsClosed(wire);

  std::vector<std::optional<TopoDS_Edge>> arcs(n); // what follows piece i
  std::vector<bool>                       used(chosen.size(), false);
  for (size_t i = 0; i < n; ++i)
  {
    if (!closed && i + 1 == n)
      break; // the end of an open wire is not a corner
    const size_t        next   = (i + 1) % n;
    const TopoDS_Vertex corner = starts[next];
    int                 which  = -1;
    for (size_t k = 0; k < chosen.size(); ++k)
      if (chosen[k].IsSame(corner))
        which = int(k);
    if (which < 0)
      continue;
    used[size_t(which)] = true;
    const auto found    = circleAt(pieces[i], pieces[next], plane, std::abs(radius), tolerance);
    if (!found)
      fail(op,
           no_solution,
           "a corner is too tight for the radius, or its pieces are tangent: no circle of that "
           "radius touches both and stays clear of what is kept");
    const gp_Pnt a = pieces[i].curve->Value(found->s), b = pieces[next].curve->Value(found->t);
    pieces[i].last     = found->s;
    pieces[next].first = found->t;
    if (radius < 0.0)
      arcs[i] = BRepBuilderAPI_MakeEdge(a, b).Edge();
    else
    {
      // The circle is made from its centre and radius, not through three points: a short arc would
      // lose its radius to the rounding of the points. Its axis is turned so that the arc runs
      // counterclockwise from a to b, by the shorter way.
      const gp_Vec toA(found->centre, a), toB(found->centre, b);
      const bool   left = gp_Vec(plane.Axis().Direction()).Dot(toA.Crossed(toB)) > 0;
      gp_Dir       axis = plane.Axis().Direction();
      if (!left)
        axis.Reverse();
      const gp_Circ      circle(gp_Ax2(found->centre, axis), std::abs(radius));
      GC_MakeArcOfCircle arc(circle, a, b, true);
      if (!arc.IsDone())
        fail(op, no_solution, "the arc of a corner could not be made");
      arcs[i] = BRepBuilderAPI_MakeEdge(arc.Value()).Edge();
    }
  }
  for (size_t k = 0; k < chosen.size(); ++k)
    if (!used[k])
      fail(op, invalid_input, "a vertex is not a corner of the wire");

  BRepBuilderAPI_MakeWire made;
  for (size_t i = 0; i < n; ++i)
  {
    const GeomAdaptor_Curve along(pieces[i].curve, pieces[i].first, pieces[i].last);
    if (GCPnts_AbscissaPoint::Length(along) <= tolerance)
      fail(op, no_solution, "a corner takes the whole of a piece: the radius is too large");
    made.Add(BRepBuilderAPI_MakeEdge(pieces[i].curve, pieces[i].first, pieces[i].last).Edge());
    if (arcs[i])
      made.Add(*arcs[i]);
  }
  if (!made.IsDone())
    fail(op, no_solution, "the trimmed pieces and the arcs do not join up");
  const TopoDS_Wire result = made.Wire();
  if (!invalidity(result).empty())
    fail(op, no_solution, "the rounded wire is not valid");
  return stored(kernel, result);
}

} // namespace BRepExt
