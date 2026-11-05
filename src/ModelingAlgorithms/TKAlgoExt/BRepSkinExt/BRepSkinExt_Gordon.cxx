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
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_Cancel.hxx>
#include <GeomAlgoExt_CubicThrough.hxx>
#include <NCollection_Array2.hxx>
#include <NCollection_HArray1.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <GCPnts_UniformAbscissa.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepFill_CompatibleWires.hxx>
#include <BSplCLib.hxx>
#include <BRep_Builder.hxx>
#include <NCollection_Sequence.hxx>
#include <functional>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_Shell.hxx>
#include <BRepOffsetExt_Focal.hxx>
#include <BRepLib.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <Bnd_Box.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Tool.hxx>
#include <GeomAPI_ExtremaCurveCurve.hxx>
#include <GeomConvert.hxx>
#include <GeomConvert_CompCurveToBSplineCurve.hxx>
#include <GeomFill_Gordon.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <NCollection_Array1.hxx>
#include <Precision.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Wire.hxx>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <optional>
#include <string>
#include <vector>

// A skin through profiles along any number of guides (two or more), as a Gordon surface: the
// surface through a network of curves that is the skin of the profiles plus the skin of the
// guides minus the surface through their crossings. OpenCascade 8 has it as GeomFill_Gordon
// (after TiGL's method), which is used here; nothing was ported.
//
// The net is put in order here first, since GeomFill_Gordon takes closed profiles only in one
// of two forms and returns a surface through its own axis for others: where every guide meets
// every profile is found (and said where one does not); closed profiles are turned to run the
// same way round and given their seams in the same place between the same two guides, away
// from every guide; what an open profile or guide has beyond the outer curves of the net is
// left out, since the surface ends at them. Closed profiles that do not run alike, whose weights
// differ at the crossings, are beyond GeomFill_Gordon: ringsNet builds their surface instead.
//
// What comes back is checked apart from GeomFill_Gordon's own report: every profile and guide,
// at eight places to each of its spans, has to lie on the surface to the tolerance; the
// normals are sampled in every actual knot span for rejection witnesses. Interval
// derivatives establish local regularity separately. At each crossing the normal must
// agree with the crossing curves. Closed flat profiles are capped for a solid.

namespace BRepExt
{
namespace
{
const char* op = "skinGordon";

// One curve for an edge or for a wire of edges that run on from each other.
Handle(Geom_BSplineCurve) curveOf(const TopoDS_Shape& shape, double tolerance)
{
  auto of = [&](const TopoDS_Edge& edge) -> Handle(Geom_BSplineCurve) {
    double     a, b;
    const auto curve = BRep_Tool::Curve(edge, a, b);
    if (curve.IsNull())
      fail(op, invalid_input, "a profile or guide has an edge without a curve");
    Handle(Geom_TrimmedCurve) trimmed = new Geom_TrimmedCurve(curve, a, b);
    if (edge.Orientation() == TopAbs_REVERSED)
      trimmed->Reverse();
    return GeomConvert::CurveToBSplineCurve(trimmed);
  };
  Handle(Geom_BSplineCurve) made;
  if (shape.ShapeType() == TopAbs_EDGE)
    made = of(TopoDS::Edge(shape));
  else if (shape.ShapeType() != TopAbs_WIRE)
    fail(op, invalid_input, "profiles and guides must be curves or edges");
  else
  {
    GeomConvert_CompCurveToBSplineCurve joined;
    bool                                any = false;
    for (BRepTools_WireExplorer it(TopoDS::Wire(shape)); it.More(); it.Next())
    {
      const auto piece = of(it.Current());
      if (!any)
        joined = GeomConvert_CompCurveToBSplineCurve(piece);
      else if (!joined.Add(piece, std::max(tolerance, Precision::Confusion()), true))
        fail(op, invalid_input, "the edges of a profile or guide do not run on from each other");
      any = true;
    }
    if (!any)
      fail(op, invalid_input, "a profile or guide is empty");
    made = joined.BSplineCurve();
  }
  if (made->IsPeriodic())
    made->SetNotPeriodic();
  return made;
}

bool closedCurve(const Handle(Geom_BSplineCurve)& c, double tolerance)
{
  return c->StartPoint().Distance(c->EndPoint()) <= tolerance;
}

//! Where two curves meet: the parameter on each, and how far apart they are there.
struct Crossing
{
  double onFirst = 0, onSecond = 0, gap = 0;
};

Crossing crossingOf(const Handle(Geom_BSplineCurve)& a, const Handle(Geom_BSplineCurve)& b)
{
  Crossing                  found;
  GeomAPI_ExtremaCurveCurve nearest(a, b);
  nearest.TotalLowerDistanceParameters(found.onFirst, found.onSecond);
  found.gap = nearest.TotalLowerDistance();
  return found;
}

gp_Vec tangentAt(const Handle(Geom_BSplineCurve)& c, double u)
{
  gp_Pnt p;
  gp_Vec d;
  c->D1(u, p, d);
  return d;
}

//! A closed curve begun at `u` instead: the same points, the same way round.
Handle(Geom_BSplineCurve) reseamed(const Handle(Geom_BSplineCurve)& c, double u, double tolerance)
{
  Handle(Geom_BSplineCurve) head = Handle(Geom_BSplineCurve)::DownCast(c->Copy());
  Handle(Geom_BSplineCurve) tail = Handle(Geom_BSplineCurve)::DownCast(c->Copy());
  head->Segment(u, c->LastParameter());
  tail->Segment(c->FirstParameter(), u);
  GeomConvert_CompCurveToBSplineCurve joined(head);
  if (!joined.Add(tail, std::max(tolerance, Precision::Confusion()), true))
    fail(op, no_solution, "a closed profile could not be begun between two guides");
  Handle(Geom_BSplineCurve) whole = joined.BSplineCurve();
  // The old seam is now a knot of full multiplicity inside the curve. Where the curve is
  // smooth there the knot is taken down as far as it goes without moving the curve.
  const double size =
    c->StartPoint().Distance(c->Value((c->FirstParameter() + c->LastParameter()) / 2));
  for (int k = whole->FirstUKnotIndex() + 1; k < whole->LastUKnotIndex(); ++k)
    if (whole->Multiplicity(k) == whole->Degree())
      for (int m = whole->Degree() - 1; m >= 1; --m)
        if (!whole->RemoveKnot(k, m, 1e-12 * std::max(size, 1.0)))
          break;
  return whole;
}

const char* said(GeomFill_Gordon::ResultStatus status)
{
  using S = GeomFill_Gordon::ResultStatus;
  switch (status)
  {
    case S::InvalidInput:
      return "too few profiles or guides";
    case S::IntersectionFailed:
      return "not every guide meets every profile";
    case S::OrderingFailed:
      return "the profiles and guides cannot be put in one order";
    case S::CompatibilityFailed:
    case S::CurveCompatibilityFailed:
      return "the profiles and guides do not form a consistent net";
    case S::RationalReparametrizationFailed:
    case S::RationalDegreeOverflow:
    case S::RationalConstructionFailed:
      return "arcs and other rational curves in the net are beyond the exact construction";
    case S::PeriodicityFailed:
      return "the closed direction could not be closed smoothly";
    default:
      return "the surface through the net could not be built";
  }
}

//! What comes back is checked apart from how it was built: every profile and guide, at eight
//! places to each of its spans, lies on the surface to the tolerance; the surface's normal nowhere
//! vanishes or turns by a right angle between neighbours of a grid of four to each span; and at
//! every crossing it is on the same side of the two curves there.
void certify(const Handle(Geom_BSplineSurface)&            surface,
             const std::vector<Handle(Geom_BSplineCurve)>& profiles,
             const std::vector<Handle(Geom_BSplineCurve)>& guides,
             const std::vector<std::vector<Crossing>>&     at,
             const std::vector<double>&                    profileWay,
             const std::vector<double>&                    guideWay,
             double                                        tolerance)
{
  const size_t np = profiles.size(), ng = guides.size();
  // Every curve of the net lies on the surface: eight places to each of its spans. A place is
  // within the tolerance of the surface as soon as some point of the surface is, so the point that
  // Gauss-Newton steps reach from the place before settles it; only where that is farther is the
  // whole surface searched, and that distance is the one held to the tolerance.
  ShapeAnalysis_Surface on(surface);
  double                u0, u1, v0, v1;
  surface->Bounds(u0, u1, v0, v1);
  const auto nearFrom = [&](gp_Pnt2d uv, const gp_Pnt& p) {
    for (int step = 0; step < 16; ++step)
    {
      gp_Pnt q;
      gp_Vec su, sv;
      surface->D1(uv.X(), uv.Y(), q, su, sv);
      const gp_Vec r(q, p);
      const double a = su.SquareMagnitude(), b = su.Dot(sv), c = sv.SquareMagnitude();
      const double det = a * c - b * b;
      if (!(det > 0))
        break;
      const double du = (c * r.Dot(su) - b * r.Dot(sv)) / det,
                   dv = (a * r.Dot(sv) - b * r.Dot(su)) / det;
      uv.SetCoord(std::clamp(uv.X() + du, u0, u1), std::clamp(uv.Y() + dv, v0, v1));
      if (std::abs(du) <= 1e-12 * (u1 - u0) && std::abs(dv) <= 1e-12 * (v1 - v0))
        break;
    }
    return uv;
  };
  double off = 0;
  for (const auto* curves : {&profiles, &guides})
    for (const auto& c : *curves)
    {
      std::optional<gp_Pnt2d> before;
      for (int k = c->FirstUKnotIndex(); k < c->LastUKnotIndex(); ++k)
        for (int s = 0; s <= 8; ++s)
        {
          const gp_Pnt p    = c->Value(c->Knot(k) + (c->Knot(k + 1) - c->Knot(k)) * s / 8);
          double       away = Precision::Infinite();
          if (before)
          {
            const gp_Pnt2d uv = nearFrom(*before, p);
            away              = surface->Value(uv.X(), uv.Y()).Distance(p);
            if (away <= tolerance)
              before = uv;
          }
          if (!(away <= tolerance))
          {
            before = on.ValueOfUV(p, Precision::Confusion());
            away   = on.Value(*before).Distance(p);
          }
          off = std::max(off, away);
        }
    }
  if (off > tolerance)
    fail(op,
         no_solution,
         "the surface built misses the profiles or guides by " + std::to_string(off));

  // Sample each actual knot rectangle, including its interior. Knot counts
  // cannot distribute samples over a nonuniform parameter domain.
  // These samples are rejection witnesses; interval regularity is checked below.
  {
    double u0, u1, v0, v1;
    surface->Bounds(u0, u1, v0, v1);
    const auto samples = [](const auto& knots) {
      std::vector<double> result;
      const double full = knots.Last() - knots.First();
      for (int k = knots.Lower(); k < knots.Upper(); ++k)
      {
        const int count = std::max(4, int(std::ceil(32 * (knots(k + 1) - knots(k)) / full)));
        for (int j = 0; j < count; ++j)
          result.push_back(knots(k) + (knots(k + 1) - knots(k)) * j / count);
      }
      result.push_back(knots.Last());
      return result;
    };
    const auto us = samples(surface->UKnots()), vs = samples(surface->VKnots());
    const int nu = int(us.size()) - 1, nv = int(vs.size()) - 1;
    std::vector<gp_Vec> normals(size_t(nu + 1) * size_t(nv + 1));
    for (int i = 0; i <= nu; ++i)
      for (int j = 0; j <= nv; ++j)
      {
        gp_Pnt p;
        gp_Vec du, dv;
        surface->D1(us[size_t(i)], vs[size_t(j)], p, du, dv);
        const gp_Vec n = du.Crossed(dv);
        if (n.Magnitude() <= 1e-9 * du.Magnitude() * dv.Magnitude())
          fail(op,
               no_solution,
               "the surface through the profiles and guides folds: it has no direction somewhere");
        normals[size_t(i) * size_t(nv + 1) + size_t(j)] = n.Normalized();
      }
    for (int i = 0; i <= nu; ++i)
      for (int j = 0; j <= nv; ++j)
      {
        const gp_Vec& n = normals[size_t(i) * size_t(nv + 1) + size_t(j)];
        if ((i < nu && n.Dot(normals[size_t(i + 1) * size_t(nv + 1) + size_t(j)]) <= 0)
            || (j < nv && n.Dot(normals[size_t(i) * size_t(nv + 1) + size_t(j) + 1]) <= 0))
          fail(op, no_solution, "the surface through the profiles and guides folds over itself");
      }
    // At every crossing the surface is on the same side of the profile and the guide.
    int side = 0;
    for (size_t i = 0; i < np; ++i)
      for (size_t j = 0; j < ng; ++j)
      {
        const gp_Vec across = (tangentAt(profiles[i], at[i][j].onFirst) * profileWay[i])
                                .Crossed(tangentAt(guides[j], at[i][j].onSecond) * guideWay[j]);
        if (across.Magnitude() <= 0)
          continue;
        const gp_Pnt2d uv =
          on.ValueOfUV(profiles[i]->Value(at[i][j].onFirst), Precision::Confusion());
        gp_Pnt p;
        gp_Vec du, dv;
        surface->D1(uv.X(), uv.Y(), p, du, dv);
        const int here = across.Dot(du.Crossed(dv)) > 0 ? 1 : -1;
        if (side != 0 && here != side)
          fail(op,
               no_solution,
               "the surface through the profiles and guides turns over between them: it would pass "
               "through itself");
        side = here;
      }
  }
}

//! A rising map through rising points: the cubic spline through them, smooth to the second
//! derivative and as steep at its ends as the first and last stretches, where that rises
//! everywhere; otherwise the monotone cubic of Fritsch and Carlson, smooth to the first only.
struct Monotone
{
  std::vector<double> x, y, m;

  Monotone(std::vector<double> xs, std::vector<double> ys)
      : x(std::move(xs)),
        y(std::move(ys)),
        m(x.size())
  {
    const size_t        n = x.size();
    std::vector<double> d(n - 1), h(n - 1);
    for (size_t k = 0; k + 1 < n; ++k)
    {
      h[k] = x[k + 1] - x[k];
      d[k] = (y[k + 1] - y[k]) / h[k];
    }
    // The spline's slopes: m0 = d0, h_k m_{k-1} + 2 (h_{k-1} + h_k) m_k + h_{k-1} m_{k+1} =
    // 3 (h_k d_{k-1} + h_{k-1} d_k), m_{n-1} = d_{n-2}.
    std::vector<double> lower(n, 0), main(n, 0), upper(n, 0), right(n, 0);
    main[0]  = 1;
    right[0] = d[0];
    for (size_t k = 1; k + 1 < n; ++k)
    {
      lower[k] = h[k];
      main[k]  = 2 * (h[k - 1] + h[k]);
      upper[k] = h[k - 1];
      right[k] = 3 * (h[k] * d[k - 1] + h[k - 1] * d[k]);
    }
    main[n - 1]  = 1;
    right[n - 1] = d[n - 2];
    for (size_t k = 1; k < n; ++k)
    {
      const double f = lower[k] / main[k - 1];
      main[k] -= f * upper[k - 1];
      right[k] -= f * right[k - 1];
    }
    m[n - 1] = right[n - 1] / main[n - 1];
    for (size_t k = n - 1; k-- > 0;)
      m[k] = (right[k] - upper[k] * m[k + 1]) / main[k];
    bool rising = true;
    for (int k = 0; k <= 256 && rising; ++k)
      rising = slope(x.front() + (x.back() - x.front()) * k / 256) > 0;
    if (rising)
      return;
    m[0]     = d[0];
    m[n - 1] = d[n - 2];
    for (size_t k = 1; k + 1 < n; ++k)
      m[k] = d[k - 1] * d[k] <= 0
               ? 0
               : 3 * (d[k - 1] + d[k]) / (2 * (d[k] / d[k - 1] + 1) + (d[k - 1] / d[k] + 1));
  }

  size_t span(double t) const
  {
    return size_t(std::clamp(std::upper_bound(x.begin(), x.end(), t) - x.begin() - 1,
                             std::ptrdiff_t(0),
                             std::ptrdiff_t(x.size()) - 2));
  }

  double slope(double t) const
  {
    const size_t k = span(t);
    const double h = x[k + 1] - x[k], s = (t - x[k]) / h;
    return (y[k] * (6 * s * s - 6 * s) + y[k + 1] * (-6 * s * s + 6 * s)) / h
           + m[k] * (3 * s * s - 4 * s + 1) + m[k + 1] * (3 * s * s - 2 * s);
  }

  double operator()(double t) const
  {
    const size_t k = span(t);
    const double h = x[k + 1] - x[k], s = (t - x[k]) / h;
    return y[k] * (2 * s * s * s - 3 * s * s + 1) + m[k] * h * (s * s * s - 2 * s * s + s)
           + y[k + 1] * (-2 * s * s * s + 3 * s * s) + m[k + 1] * h * (s * s * s - s * s);
  }
};

//! A curve between two of its parameters made again as a spline whose parameter is `to(s)`, s the
//! share of its length up to the point, or of its own parameter where `byLength` is false: the
//! spline through points of the curve, the given parameters of the curve among them, more of them
//! until the spline is within a quarter of the tolerance of the curve between them too.
Handle(Geom_BSplineCurve) refit(const Handle(Geom_BSplineCurve)& c,
                                double                           a,
                                double                           b,
                                const Monotone&                  to,
                                const std::vector<double>&       through,
                                bool                             byLength,
                                double                           tolerance)
{
  GeomAdaptor_Curve along(c, a, b);
  auto              lengthOf = [&](double from, double to) {
    return byLength ? GCPnts_AbscissaPoint::Length(along, from, to, 1e-12) : to - from;
  };
  if (!(lengthOf(a, b) > tolerance))
    fail(op, invalid_input, "a profile or guide has no length between its crossings");
  for (int n = 16 * std::max(1, c->NbKnots() - 1);; n *= 2)
  {
    if (n > 4096)
      fail(op,
           no_solution,
           "a profile or guide could not be made again with the guides at one parameter");
    checkCancelled(op);
    // Even places, and the given ones among them in place of any even one within a quarter
    // of a step: points closer than that would only unsettle the spline.
    std::vector<double> ts;
    for (const double t : through)
      if (t > a && t < b)
        ts.push_back(t);
    const size_t given = ts.size();
    for (int k = 0; k <= n; ++k)
    {
      const double t = a + (b - a) * k / n;
      if (std::none_of(ts.begin(), ts.begin() + std::ptrdiff_t(given), [&](double g) {
            return std::abs(g - t) < (b - a) / (4.0 * n);
          }))
        ts.push_back(t);
    }
    std::sort(ts.begin(), ts.end());
    // The share of the length at each, and halfway to the next, of one sum of the pieces.
    std::vector<double> shares{0.0}, halves;
    for (size_t k = 1; k < ts.size(); ++k)
    {
      const double middle = (ts[k - 1] + ts[k]) / 2;
      halves.push_back(shares.back() + lengthOf(ts[k - 1], middle));
      shares.push_back(halves.back() + lengthOf(middle, ts[k]));
    }
    const double length = shares.back();
    for (double& x : shares)
      x /= length;
    for (double& x : halves)
      x /= length;
    Handle(NCollection_HArray1<gp_Pnt>) points = new NCollection_HArray1<gp_Pnt>(1, int(ts.size()));
    Handle(NCollection_HArray1<double>) params = new NCollection_HArray1<double>(1, int(ts.size()));
    for (size_t k = 0; k < ts.size(); ++k)
    {
      points->SetValue(int(k) + 1, c->Value(ts[k]));
      params->SetValue(int(k) + 1, to(shares[k]));
    }
    GeomAPI_Interpolate fit(points, params, false, Precision::Confusion());
    fit.Perform();
    if (!fit.IsDone())
      continue;
    const Handle(Geom_BSplineCurve) made  = fit.Curve();
    double                          worst = 0;
    for (size_t k = 0; k + 1 < ts.size(); ++k)
      worst =
        std::max(worst, made->Value(to(halves[k])).Distance(c->Value((ts[k] + ts[k + 1]) / 2)));
    if (worst <= tolerance / 4)
      return made;
  }
}

//! Whether two curves run alike: one degree, and the same knots and weights at the same shares of
//! their parameters.
bool runAlike(const Handle(Geom_BSplineCurve)& a, const Handle(Geom_BSplineCurve)& b)
{
  if (a->Degree() != b->Degree() || a->NbKnots() != b->NbKnots() || a->NbPoles() != b->NbPoles())
    return false;
  const auto share = [](const Handle(Geom_BSplineCurve)& c, int k) {
    return (c->Knot(k) - c->Knot(1)) / (c->Knot(c->NbKnots()) - c->Knot(1));
  };
  for (int k = 1; k <= a->NbKnots(); ++k)
    if (a->Multiplicity(k) != b->Multiplicity(k)
        || std::abs(share(a, k) - share(b, k)) > tolerance::parameterMerge)
      return false;
  for (int k = 1; k <= a->NbPoles(); ++k)
    if (std::abs(a->Weight(k) - b->Weight(k)) > tolerance::parameterMerge)
      return false;
  return true;
}

//! Curves made again as cubics on one grid of parameters from 0 to 1, each parameter at the share
//! of its curve's length that `shares[i]` gives, closed ones periodic: same knots, no weights, the
//! parameters `given` among the grid's. Each is within 1e-8 of its length of its curve; the grid
//! is refined for all together. Curves that are not made periodic keep their directions at their
//! ends. Empty where 4096 points do not bring one there (a corner along it). The accuracy is set by
//! the closed forms of rings on a cone or a cylinder. A grid found too coarse is doubled and every
//! placement made again, so the cost is about twice that of the grid that is accepted.
std::vector<Handle(Geom_BSplineCurve)> onOneGrid(
  const std::vector<Handle(Geom_BSplineCurve)>& curves,
  const std::vector<Monotone>&                  shares,
  const std::vector<double>&                    given,
  bool                                          periodic)
{
  std::vector<GeomAdaptor_Curve> along;
  std::vector<double>            whole, within;
  for (const auto& c : curves)
  {
    along.emplace_back(c);
    whole.push_back(
      GCPnts_AbscissaPoint::Length(along.back(), c->FirstParameter(), c->LastParameter()));
    within.push_back(tolerance::resolved(1e-8, whole.back(), c->StartPoint().XYZ().Modulus()));
  }
  // The points of a curve at rising parameters, each found on from the one before.
  const auto pointsAt = [&](size_t i, const std::vector<double>& rising) {
    std::vector<gp_Pnt> points;
    double              from = curves[i]->FirstParameter(), reached = 0;
    for (const double t : rising)
    {
      const double length = std::clamp(shares[i](t), 0.0, 1.0) * whole[i];
      // Each step to its share of the tolerance: the errors of the steps add up along the curve.
      GCPnts_AbscissaPoint place(along[i],
                                 length - reached,
                                 from,
                                 within[i] / (16 * double(rising.size())));
      if (place.IsDone())
      {
        from    = place.Parameter();
        reached = length;
      }
      points.push_back(curves[i]->Value(from));
    }
    return points;
  };
  const auto endTangent = [&](size_t i, bool last) {
    const double u = last ? curves[i]->LastParameter() : curves[i]->FirstParameter();
    return tangentAt(curves[i], u).Normalized() * (whole[i] * shares[i].slope(last ? 1.0 : 0.0));
  };
  for (int n = periodic ? 4 : 1; n <= 4096; n *= 2)
  {
    checkCancelled(op);
    std::vector<double> grid;
    for (const double g : given)
      if (g > 0 && g < 1)
        grid.push_back(g);
    const size_t kept = grid.size();
    for (int k = 0; k <= n; ++k)
    {
      const double t = double(k) / n;
      if (std::none_of(grid.begin(), grid.begin() + std::ptrdiff_t(kept), [&](double g) {
            return std::abs(g - t) < 0.25 / n;
          }))
        grid.push_back(t);
    }
    std::sort(grid.begin(), grid.end());
    const int                           count = int(grid.size()) - (periodic ? 1 : 0);
    Handle(NCollection_HArray1<double>) params =
      new NCollection_HArray1<double>(1, int(grid.size()));
    std::vector<double> rising;
    for (size_t k = 0; k < grid.size(); ++k)
    {
      params->SetValue(int(k) + 1, grid[k]);
      rising.push_back(grid[k]);
      if (k + 1 < grid.size())
        rising.push_back((grid[k] + grid[k + 1]) / 2);
    }
    std::vector<Handle(Geom_BSplineCurve)> made;
    bool                                   close = true;
    for (size_t i = 0; i < curves.size() && close; ++i)
    {
      const std::vector<gp_Pnt>           points  = pointsAt(i, rising);
      Handle(NCollection_HArray1<gp_Pnt>) through = new NCollection_HArray1<gp_Pnt>(1, count);
      for (int k = 0; k < count; ++k)
        through->SetValue(k + 1, points[size_t(2 * k)]);
      GeomAPI_Interpolate fit(through, params, periodic, within[i]);
      if (!periodic)
        fit.Load(endTangent(i, false), endTangent(i, true), false);
      fit.Perform();
      if (!fit.IsDone())
        return {};
      const Handle(Geom_BSplineCurve) curve = fit.Curve();
      for (size_t k = 1; k < rising.size() && close; k += 2)
        close = curve->Value(rising[k]).Distance(points[k]) <= within[i];
      made.push_back(curve);
    }
    if (close)
      return made;
  }
  return {};
}

//! The surface through closed rings and guides, as the sum of the skin of the rings and the skin of
//! the guides less the skin of their crossings. All of them are made again on one grid each, every
//! guide at one parameter of the rings, the mean of the shares of their lengths at which it crosses
//! them, and every ring at one of the guides likewise. Null where a curve cannot be made again or a
//! guide does not run from ring to ring in their order, for GeomFill_Gordon to say what it finds.
Handle(Geom_BSplineSurface) ringsNet(std::vector<Handle(Geom_BSplineCurve)>&       rings,
                                     const std::vector<Handle(Geom_BSplineCurve)>& guidesGiven,
                                     std::vector<std::vector<Crossing>>&           at,
                                     double                                        tolerance)
{
  const size_t np = rings.size(), ng = guidesGiven.size();
  const auto   shareTo = [](const Handle(Geom_BSplineCurve)& c, double u) {
    GeomAdaptor_Curve along(c);
    return GCPnts_AbscissaPoint::Length(along, c->FirstParameter(), u)
           / GCPnts_AbscissaPoint::Length(along, c->FirstParameter(), c->LastParameter());
  };
  // Every seam halfway along its ring's stretch between the two guides around it, so that rings
  // of one shape have every guide at one share of their length, and nothing to warp.
  for (size_t i = 0; i < np; ++i)
  {
    double low = 1, high = 0;
    for (size_t j = 0; j < ng; ++j)
    {
      const double share = shareTo(rings[i], at[i][j].onFirst);
      low                = std::min(low, share);
      high               = std::max(high, share);
    }
    const double      middle = (low + 1 + high) / 2;
    const double      share  = middle >= 1 ? middle - 1 : middle;
    GeomAdaptor_Curve along(rings[i]);
    const double      whole =
      GCPnts_AbscissaPoint::Length(along, rings[i]->FirstParameter(), rings[i]->LastParameter());
    GCPnts_AbscissaPoint place(along, share * whole, rings[i]->FirstParameter(), tolerance / 16);
    if (!place.IsDone())
      return nullptr;
    if (rings[i]->Value(place.Parameter()).Distance(rings[i]->StartPoint()) > tolerance)
    {
      rings[i] = reseamed(rings[i], place.Parameter(), tolerance);
      for (size_t j = 0; j < ng; ++j)
        at[i][j] = crossingOf(rings[i], guidesGiven[j]);
    }
  }
  // The guides in the order they cross the first ring from its seam.
  std::vector<size_t> order(ng);
  std::iota(order.begin(), order.end(), size_t(0));
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    return at[0][a].onFirst < at[0][b].onFirst;
  });
  std::vector<double>              nodes{0.0};
  std::vector<std::vector<double>> onRing(np, std::vector<double>{0.0});
  for (size_t k = 0; k < ng; ++k)
  {
    double mean = 0;
    for (size_t i = 0; i < np; ++i)
    {
      onRing[i].push_back(shareTo(rings[i], at[i][order[k]].onFirst));
      mean += onRing[i].back() / double(np);
    }
    nodes.push_back(mean);
  }
  nodes.push_back(1.0);
  std::vector<Monotone> ringShares;
  for (size_t i = 0; i < np; ++i)
  {
    onRing[i].push_back(1.0);
    for (size_t k = 1; k < onRing[i].size(); ++k)
      if (onRing[i][k] <= onRing[i][k - 1] || nodes[k] <= nodes[k - 1])
        return nullptr;
    ringShares.emplace_back(nodes, onRing[i]);
  }
  // Each guide from the first ring to the last, and the rings' places along them.
  std::vector<Handle(Geom_BSplineCurve)> runs;
  std::vector<std::vector<double>>       onGuide(ng);
  std::vector<double>                    places(np, 0.0);
  for (size_t k = 0; k < ng; ++k)
  {
    const size_t j = order[k];
    if (closedCurve(guidesGiven[j], tolerance))
      return nullptr;
    const Handle(Geom_BSplineCurve) run =
      Handle(Geom_BSplineCurve)::DownCast(guidesGiven[j]->Copy());
    const bool   back  = at[0][j].onSecond > at[np - 1][j].onSecond;
    const double first = run->FirstParameter(), last = run->LastParameter();
    if (back)
      run->Reverse();
    for (size_t i = 0; i < np; ++i)
    {
      onGuide[k].push_back(
        shareTo(run, back ? first + last - at[i][j].onSecond : at[i][j].onSecond));
      places[i] += onGuide[k].back() / double(ng);
    }
    runs.push_back(run);
  }
  std::vector<Monotone> guideShares;
  for (size_t k = 0; k < ng; ++k)
  {
    for (size_t i = 1; i < np; ++i)
      if (onGuide[k][i] <= onGuide[k][i - 1] || places[i] <= places[i - 1])
        return nullptr;
    guideShares.emplace_back(places, onGuide[k]);
  }
  const auto newRings  = onOneGrid(rings, ringShares, nodes, false);
  const auto newGuides = onOneGrid(runs, guideShares, places, false);
  if (newRings.empty() || newGuides.empty())
    return nullptr;

  // The skin of the rings: each column of their poles through the rings' places.
  const GeomAlgoExt::CubicThrough across(places);
  const int                       nu = newRings[0]->NbPoles();
  NCollection_Array2<gp_Pnt>      ringPoles(1, nu, 1, across.poles);
  for (int a = 1; a <= nu; ++a)
  {
    std::vector<gp_Pnt> values;
    for (const auto& c : newRings)
      values.push_back(c->Pole(a));
    const auto column = across.through(values);
    for (int b = 1; b <= across.poles; ++b)
      ringPoles(a, b) = column[size_t(b - 1)];
  }
  Handle(Geom_BSplineSurface) ringSkin = new Geom_BSplineSurface(ringPoles,
                                                                 newRings[0]->Knots(),
                                                                 across.knots,
                                                                 newRings[0]->Multiplicities(),
                                                                 across.mults,
                                                                 3,
                                                                 3);

  // The closed spline round through values at the seam and at the guides' parameters.
  const auto round = [&](const std::vector<gp_Pnt>& values) {
    Handle(NCollection_HArray1<gp_Pnt>) points =
      new NCollection_HArray1<gp_Pnt>(1, int(values.size()));
    Handle(NCollection_HArray1<double>) params =
      new NCollection_HArray1<double>(1, int(nodes.size()));
    for (size_t k = 0; k < values.size(); ++k)
      points->SetValue(int(k) + 1, values[k]);
    for (size_t k = 0; k < nodes.size(); ++k)
      params->SetValue(int(k) + 1, nodes[k]);
    GeomAPI_Interpolate fit(points, params, true, 0.0);
    fit.Perform();
    Handle(Geom_BSplineCurve) curve;
    if (fit.IsDone())
    {
      curve = fit.Curve();
      // Not periodic it still runs on past its period: clamped by a segment of one period.
      curve->SetNotPeriodic();
      curve->Segment(nodes.front(), nodes.back());
    }
    return curve;
  };
  // The skin of the guides, and with the seam as one of them, the iso of the rings' skin there.
  const Handle(Geom_BSplineCurve) seam = Handle(Geom_BSplineCurve)::DownCast(ringSkin->UIso(0.0));
  for (int k = 2; k < newGuides[0]->NbKnots(); ++k)
    seam->InsertKnot(newGuides[0]->Knot(k), newGuides[0]->Multiplicity(k), 1e-12, false);
  const int nv = newGuides[0]->NbPoles();
  if (seam->NbPoles() != nv)
    return nullptr;
  Handle(Geom_BSplineSurface) guideSkin, crossings;
  {
    std::vector<Handle(Geom_BSplineCurve)> rows;
    for (int b = 1; b <= nv; ++b)
    {
      std::vector<gp_Pnt> values{seam->Pole(b)};
      for (const auto& g : newGuides)
        values.push_back(g->Pole(b));
      rows.push_back(round(values));
      if (rows.back().IsNull())
        return nullptr;
    }
    NCollection_Array2<gp_Pnt> poles(1, rows[0]->NbPoles(), 1, nv);
    for (int a = 1; a <= rows[0]->NbPoles(); ++a)
      for (int b = 1; b <= nv; ++b)
        poles(a, b) = rows[size_t(b - 1)]->Pole(a);
    guideSkin = new Geom_BSplineSurface(poles,
                                        rows[0]->Knots(),
                                        newGuides[0]->Knots(),
                                        rows[0]->Multiplicities(),
                                        newGuides[0]->Multiplicities(),
                                        3,
                                        3);
  }
  {
    // The crossings: closed round each ring, then each column across.
    std::vector<Handle(Geom_BSplineCurve)> rows;
    for (size_t i = 0; i < np; ++i)
    {
      std::vector<gp_Pnt> values{newRings[i]->StartPoint()};
      for (size_t k = 1; k + 1 < nodes.size(); ++k)
        values.push_back(newRings[i]->Value(nodes[k]));
      rows.push_back(round(values));
      if (rows.back().IsNull())
        return nullptr;
    }
    const int                  count = rows[0]->NbPoles();
    NCollection_Array2<gp_Pnt> poles(1, count, 1, across.poles);
    for (int a = 1; a <= count; ++a)
    {
      std::vector<gp_Pnt> values;
      for (const auto& c : rows)
        values.push_back(c->Pole(a));
      const auto column = across.through(values);
      for (int b = 1; b <= across.poles; ++b)
        poles(a, b) = column[size_t(b - 1)];
    }
    crossings = new Geom_BSplineSurface(poles,
                                        rows[0]->Knots(),
                                        across.knots,
                                        rows[0]->Multiplicities(),
                                        across.mults,
                                        3,
                                        3);
  }
  // On one set of knots, and summed.
  const std::vector<Handle(Geom_BSplineSurface)> skins{ringSkin, guideSkin, crossings};
  for (const auto& from : skins)
    for (const auto& to : skins)
    {
      for (int k = 2; k < from->NbUKnots(); ++k)
        to->InsertUKnot(from->UKnot(k), from->UMultiplicity(k), 1e-12, false);
      for (int k = 2; k < from->NbVKnots(); ++k)
        to->InsertVKnot(from->VKnot(k), from->VMultiplicity(k), 1e-12, false);
    }
  if (ringSkin->NbUPoles() != guideSkin->NbUPoles() || ringSkin->NbUPoles() != crossings->NbUPoles()
      || ringSkin->NbVPoles() != guideSkin->NbVPoles()
      || ringSkin->NbVPoles() != crossings->NbVPoles())
    return nullptr;
  for (int a = 1; a <= ringSkin->NbUPoles(); ++a)
    for (int b = 1; b <= ringSkin->NbVPoles(); ++b)
      ringSkin->SetPole(a,
                        b,
                        gp_Pnt(ringSkin->Pole(a, b).XYZ() + guideSkin->Pole(a, b).XYZ()
                               - crossings->Pole(a, b).XYZ()));
  return ringSkin;
}

} // namespace

EvolutionData skinGordon(OcctKernel&                  kernel,
                         const std::vector<uint32_t>& profileIds,
                         const std::vector<uint32_t>& guideIds,
                         bool                         solid,
                         double                       tolerance)
{
  if (profileIds.size() < 2 || guideIds.size() < 2 || !std::isfinite(tolerance) || tolerance <= 0)
    fail(op, invalid_input, "two or more profiles and two or more guides are required");
  const size_t                           np = profileIds.size(), ng = guideIds.size();
  std::vector<Handle(Geom_BSplineCurve)> profiles(np), guides(ng);
  for (size_t i = 0; i < np; ++i)
    profiles[i] = curveOf(get(kernel, profileIds[i]), tolerance);
  for (size_t j = 0; j < ng; ++j)
    guides[j] = curveOf(get(kernel, guideIds[j]), tolerance);

  // Where every guide meets every profile.
  std::vector<std::vector<Crossing>> at(np, std::vector<Crossing>(ng));
  auto                               meet = [&](size_t i) {
    for (size_t j = 0; j < ng; ++j)
    {
      at[i][j] = crossingOf(profiles[i], guides[j]);
      if (at[i][j].gap > tolerance)
        fail(
          op,
          no_solution,
          "no surface through the profiles and guides: not every guide meets every profile (guide "
            + std::to_string(j + 1) + " passes profile " + std::to_string(i + 1)
            + " at a distance of " + std::to_string(at[i][j].gap) + ")");
    }
  };
  for (size_t i = 0; i < np; ++i)
    meet(i);

  size_t closedProfiles = 0;
  for (const auto& c : profiles)
    closedProfiles += closedCurve(c, tolerance) ? 1 : 0;
  if (closedProfiles != 0 && closedProfiles != np)
    fail(op, invalid_input, "the profiles have to be all closed or all open");
  const bool tube = closedProfiles == np;

  // Which way each curve runs against the first of its kind: by its direction at the
  // crossings they share.
  std::vector<double> profileWay(np, 1.0), guideWay(ng, 1.0);
  for (size_t i = 1; i < np; ++i)
  {
    double agree = 0;
    for (size_t j = 0; j < ng; ++j)
      agree += tangentAt(profiles[i], at[i][j].onFirst)
                 .Normalized()
                 .Dot(tangentAt(profiles[0], at[0][j].onFirst).Normalized());
    profileWay[i] = agree < 0 ? -1.0 : 1.0;
  }
  for (size_t j = 1; j < ng; ++j)
  {
    double agree = 0;
    for (size_t i = 0; i < np; ++i)
      agree += tangentAt(guides[j], at[i][j].onSecond)
                 .Normalized()
                 .Dot(tangentAt(guides[0], at[i][0].onSecond).Normalized());
    guideWay[j] = agree < 0 ? -1.0 : 1.0;
  }

  if (tube)
  {
    // All the same way round.
    for (size_t i = 1; i < np; ++i)
      if (profileWay[i] < 0)
      {
        profiles[i]->Reverse();
        profileWay[i] = 1.0;
        meet(i);
      }
    // The guides in their order round the first profile, and the widest stretch of it
    // between two neighbours: the seams go into that stretch.
    std::vector<size_t> round(ng);
    std::iota(round.begin(), round.end(), size_t(0));
    std::sort(round.begin(), round.end(), [&](size_t a, size_t b) {
      return at[0][a].onFirst < at[0][b].onFirst;
    });
    const double period0 = profiles[0]->LastParameter() - profiles[0]->FirstParameter();
    size_t       widest  = 0;
    double       width   = -1;
    for (size_t k = 0; k < ng; ++k)
    {
      double w = at[0][round[(k + 1) % ng]].onFirst - at[0][round[k]].onFirst;
      if (w <= 0)
        w += period0;
      if (w > width)
      {
        width  = w;
        widest = k;
      }
    }
    const size_t before = round[widest], after = round[(widest + 1) % ng];
    for (size_t i = 0; i < np; ++i)
    {
      const double first = profiles[i]->FirstParameter(), last = profiles[i]->LastParameter();
      const double a = at[i][before].onFirst, b = at[i][after].onFirst;
      // The stretch from `before` on to `after`; it holds the seam already where it wraps
      // round the curve's ends and neither guide is at them.
      const gp_Pnt seam  = profiles[i]->StartPoint();
      const bool   holds = a > b && profiles[i]->Value(a).Distance(seam) > tolerance
                         && profiles[i]->Value(b).Distance(seam) > tolerance;
      if (!holds)
      {
        double middle = b > a ? (a + b) / 2 : (a + b + (last - first)) / 2;
        if (b <= a && a - b < BRepExt::tolerance::parameterMerge * (last - first))
          middle = a + (last - first) / 2; // one guide given twice, or both at the seam
        if (middle > last)
          middle -= last - first;
        if (middle <= first || middle >= last)
          fail(op, no_solution, "a closed profile could not be begun between two guides");
        profiles[i] = reseamed(profiles[i], middle, tolerance);
        meet(i);
      }
      // Round every profile the guides come in the same order.
      std::vector<size_t> here(ng);
      std::iota(here.begin(), here.end(), size_t(0));
      std::sort(here.begin(), here.end(), [&](size_t x, size_t y) {
        return at[i][x].onFirst < at[i][y].onFirst;
      });
      const size_t shift = size_t(std::find(round.begin(), round.end(), here[0]) - round.begin());
      for (size_t k = 0; k < ng; ++k)
        if (here[k] != round[(shift + k) % ng])
          fail(
            op,
            no_solution,
            "no surface through the profiles and guides: the guides do not come in the same order "
            "round every profile");
    }
  }
  else
  {
    // What a profile has beyond the outer guides is left out.
    for (size_t i = 0; i < np; ++i)
    {
      double low = at[i][0].onFirst, high = low;
      for (size_t j = 1; j < ng; ++j)
      {
        low  = std::min(low, at[i][j].onFirst);
        high = std::max(high, at[i][j].onFirst);
      }
      const auto& c = profiles[i];
      if (c->Value(low).Distance(c->StartPoint()) > tolerance
          || c->Value(high).Distance(c->EndPoint()) > tolerance)
      {
        if (high - low
            <= BRepExt::tolerance::parameter * (c->LastParameter() - c->FirstParameter()))
          fail(op, no_solution, "the guides all meet a profile in one point");
        profiles[i]->Segment(low, high);
        meet(i);
      }
    }
  }
  // And what a guide has beyond the outer profiles.
  for (size_t j = 0; j < ng; ++j)
  {
    if (closedCurve(guides[j], tolerance))
      continue;
    double low = at[0][j].onSecond, high = low;
    for (size_t i = 1; i < np; ++i)
    {
      low  = std::min(low, at[i][j].onSecond);
      high = std::max(high, at[i][j].onSecond);
    }
    const auto& c = guides[j];
    if (c->Value(low).Distance(c->StartPoint()) > tolerance
        || c->Value(high).Distance(c->EndPoint()) > tolerance)
    {
      if (high - low <= BRepExt::tolerance::parameter * (c->LastParameter() - c->FirstParameter()))
        fail(op, no_solution, "the profiles all meet a guide in one point");
      guides[j]->Segment(low, high);
      for (size_t i = 0; i < np; ++i)
        at[i][j] = crossingOf(profiles[i], guides[j]);
    }
  }

  // Rings that do not run alike, of unlike spans or with arcs at unlike angles, have weights that
  // differ at the crossings, which GeomFill_Gordon takes as no net.
  Handle(Geom_BSplineSurface) network;
  if (tube
      && std::any_of(profiles.begin(),
                     profiles.end(),
                     [](const auto& c) { return c->IsRational(); })
      && !std::all_of(profiles.begin(), profiles.end(), [&](const auto& c) {
           return runAlike(profiles[0], c);
         }))
    network = ringsNet(profiles, guides, at, tolerance);

  NCollection_Array1<Handle(Geom_Curve)> profileList(1, int(np)), guideList(1, int(ng));
  for (size_t i = 0; i < np; ++i)
    profileList(int(i) + 1) = profiles[i];
  for (size_t j = 0; j < ng; ++j)
    guideList(int(j) + 1) = guides[j];
  Handle(Geom_BSplineSurface) surface = network;
  if (surface.IsNull())
  {
    GeomFill_Gordon gordon;
    NCollection_Array2<double> profileParams(1, int(np), 1, int(ng));
    NCollection_Array2<double> guideParams(1, int(np), 1, int(ng));
    for (size_t i = 0; i < np; ++i)
      for (size_t j = 0; j < ng; ++j)
      {
        profileParams.ChangeAt(i, j) = at[i][j].onFirst;
        guideParams.ChangeAt(i, j) = at[i][j].onSecond;
      }
    gordon.Init(profileList, guideList, profileParams, guideParams, tolerance);
    gordon.Perform();
    if (!gordon.IsDone() || gordon.Surface().IsNull())
      fail(op,
           no_solution,
           std::string("no surface through the profiles and guides: ") + said(gordon.Status()));
    if (gordon.IsApproximate())
      fail(op, no_solution, "the surface through the profiles and guides is only approximate");
    surface = gordon.Surface();
  }

  certify(surface, profiles, guides, at, profileWay, guideWay, tolerance);

  BRepBuilderAPI_MakeFace face(surface, Precision::Confusion());
  if (!face.IsDone())
    fail(op, no_solution, "the surface through the net could not be made a face");
  TopoDS_Shape result = face.Face();
  if (solid)
  {
    // The first and last profile, closed and flat, close the tube.
    BRepBuilderAPI_Sewing sewing(tolerance);
    sewing.Add(result);
    for (const uint32_t id : {profileIds.front(), profileIds.back()})
    {
      const auto& shape = get(kernel, id);
      if (shape.ShapeType() != TopAbs_WIRE || !BRep_Tool::IsClosed(TopoDS::Wire(shape)))
        fail(op, invalid_input, "a solid needs closed first and last profiles");
      // Of a copy: a face made of the caller's wire marks that wire as part of something.
      BRepBuilderAPI_MakeFace cap(TopoDS::Wire(BRepBuilderAPI_Copy(shape).Shape()), true);
      if (!cap.IsDone())
        fail(op, invalid_input, "a solid needs flat first and last profiles");
      sewing.Add(cap.Face());
    }
    sewing.Perform();
    TopExp_Explorer shell(sewing.SewedShape(), TopAbs_SHELL);
    if (!shell.More() || !BRep_Tool::IsClosed(shell.Current()))
      fail(op, no_solution, "the skin and its two ends do not close into a solid");
    BRepBuilderAPI_MakeSolid made(TopoDS::Shell(shell.Current()));
    if (!made.IsDone())
      fail(op, no_solution, "the skin and its two ends do not close into a solid");
    TopoDS_Solid body = made.Solid();
    BRepLib::OrientClosedSolid(body);
    result = body;
  }
  // The zero-distance normal-sweep certificate proves a nonsingular first
  // fundamental form on every knot rectangle using outward-rounded intervals.
  // Unknown is not a regular verdict. This is a local immersion certificate.
  if (focalSweep(face.Face(), 0) != FocalSweep::Regular)
    fail(op, no_solution, "regularity of the network surface could not be established");
  if (!BRepCheck_Analyzer(result).IsValid())
    fail(op, no_solution, "the skin through the net is not a valid body");
  return stored(kernel, result);
}

// A skin through open profiles along guides that cross them anywhere. Each guide is put at one
// parameter across the profiles, the mean of the shares of their lengths at which it crosses them;
// the profiles are made again with that parameter (a monotone smooth map from their length to it),
// and skinned. Each guide then moves the skin by how far it is off it, blended across to the next
// guide or to the profiles' ends by the cubic that is one at the guide and nothing at the next,
// with no slope at either (a convention of this code, not derived).
EvolutionData skinGuided(OcctKernel&                  kernel,
                         const std::vector<uint32_t>& profileIds,
                         const std::vector<uint32_t>& guideIds,
                         double                       tolerance)
{
  if (profileIds.size() < 2 || guideIds.empty() || !std::isfinite(tolerance) || tolerance <= 0)
    fail(op, invalid_input, "two or more profiles and one or more guides are required");
  checkCancelled(op);
  const size_t                           np = profileIds.size(), ng = guideIds.size();
  std::vector<Handle(Geom_BSplineCurve)> profiles(np), guides(ng);
  for (size_t i = 0; i < np; ++i)
    profiles[i] = curveOf(get(kernel, profileIds[i]), tolerance);
  for (size_t j = 0; j < ng; ++j)
    guides[j] = curveOf(get(kernel, guideIds[j]), tolerance);
  for (const auto* curves : {&profiles, &guides})
    for (const auto& c : *curves)
      if (closedCurve(c, tolerance))
        fail(op, invalid_input, "this skin is of open profiles along open guides");
  std::vector<std::vector<Crossing>> at(np, std::vector<Crossing>(ng));
  auto                               meet = [&]() {
    for (size_t i = 0; i < np; ++i)
      for (size_t j = 0; j < ng; ++j)
      {
        at[i][j] = crossingOf(profiles[i], guides[j]);
        if (at[i][j].gap > tolerance)
          fail(op,
               no_solution,
               "no surface through the profiles and guides: not every guide meets every profile "
                                             "(guide "
                 + std::to_string(j + 1) + " passes profile " + std::to_string(i + 1)
                 + " at a distance of " + std::to_string(at[i][j].gap) + ")");
      }
  };
  meet();
  // The profiles run alike: by their directions where the first guide crosses them.
  for (size_t i = 1; i < np; ++i)
    if (tangentAt(profiles[i], at[i][0].onFirst).Dot(tangentAt(profiles[0], at[0][0].onFirst)) < 0)
      profiles[i]->Reverse();
  // Each guide from the first profile to the last.
  meet();
  for (size_t j = 0; j < ng; ++j)
    if (at[0][j].onSecond > at[np - 1][j].onSecond)
      guides[j]->Reverse();
  meet();

  // Where each guide crosses each profile, as a share of the profile's length; the guides in
  // the order they cross the first, which has to be the order they cross every other.
  auto lengthTo = [](const Handle(Geom_BSplineCurve)& c, double u) {
    GeomAdaptor_Curve along(c);
    return GCPnts_AbscissaPoint::Length(along, c->FirstParameter(), u);
  };
  std::vector<std::vector<double>> share(np, std::vector<double>(ng));
  for (size_t i = 0; i < np; ++i)
  {
    const double whole = lengthTo(profiles[i], profiles[i]->LastParameter());
    for (size_t j = 0; j < ng; ++j)
      share[i][j] = lengthTo(profiles[i], at[i][j].onFirst) / whole;
  }
  std::vector<size_t> order(ng);
  std::iota(order.begin(), order.end(), size_t(0));
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    return share[0][a] < share[0][b];
  });
  // A guide that meets the profiles' ends is their boundary there; otherwise the profiles' ends
  // are, and the surface is free to follow them.
  auto atEnd = [&](size_t j, bool last) {
    for (size_t i = 0; i < np; ++i)
    {
      const gp_Pnt end = last ? profiles[i]->EndPoint() : profiles[i]->StartPoint();
      if (profiles[i]->Value(at[i][j].onFirst).Distance(end) > tolerance)
        return false;
    }
    return true;
  };
  std::vector<double> nodes;
  std::vector<size_t> nodeOf(ng);
  if (!atEnd(order[0], false))
    nodes.push_back(0.0);
  for (size_t k = 0; k < ng; ++k)
  {
    const size_t j    = order[k];
    double       mean = 0;
    for (size_t i = 0; i < np; ++i)
      mean += share[i][j] / double(np);
    nodeOf[k] = nodes.size();
    nodes.push_back(atEnd(j, false) ? 0.0 : atEnd(j, true) ? 1.0 : mean);
  }
  if (nodes.back() < 1.0)
    nodes.push_back(1.0);
  // Where each profile has the nodes.
  std::vector<std::vector<double>> on(np);
  for (size_t i = 0; i < np; ++i)
  {
    if (nodes.front() == 0.0 && (ng == 0 || nodeOf[0] != 0))
      on[i].push_back(0.0);
    for (size_t k = 0; k < ng; ++k)
    {
      const size_t j = order[k];
      on[i].push_back(nodes[nodeOf[k]] == 0.0 || nodes[nodeOf[k]] == 1.0 ? nodes[nodeOf[k]]
                                                                         : share[i][j]);
    }
    if (on[i].size() < nodes.size())
      on[i].push_back(1.0);
    for (size_t k = 1; k < on[i].size(); ++k)
      if (on[i][k] - on[i][k - 1] <= BRepExt::tolerance::parameter)
        fail(op,
             no_solution,
             "no surface through the profiles and guides: the guides meet a profile at one point, "
             "or in "
             "another order than they meet the first");
  }

  // The profiles again, each guide at its own parameter.
  std::vector<Handle(Geom_BSplineCurve)> sections(np);
  for (size_t i = 0; i < np; ++i)
  {
    std::vector<double> crossings;
    for (size_t j = 0; j < ng; ++j)
      crossings.push_back(at[i][j].onFirst);
    sections[i] = refit(profiles[i],
                        profiles[i]->FirstParameter(),
                        profiles[i]->LastParameter(),
                        Monotone(on[i], nodes),
                        crossings,
                        true,
                        tolerance);
  }
  // The profiles' places across: by how far apart they are on the average.
  std::vector<double> across{0.0};
  for (size_t i = 1; i < np; ++i)
  {
    double apart = 0;
    for (int k = 0; k <= 8; ++k)
      apart += sections[i]->Value(k / 8.0).Distance(sections[i - 1]->Value(k / 8.0)) / 9;
    if (apart <= tolerance)
      fail(op, invalid_input, "two profiles lie on each other");
    across.push_back(across.back() + apart);
  }
  for (double& v : across)
    v /= across.back();
  GeomAlgoExt::unify(sections);

  // The profiles' own surface: each column of poles through the profiles, at their places.
  const Handle(Geom_BSplineCurve)&       base = sections[0];
  const int                              nu   = base->NbPoles();
  std::vector<Handle(Geom_BSplineCurve)> columns(static_cast<size_t>(nu));
  for (int a = 1; a <= nu; ++a)
  {
    if (np == 2)
    {
      NCollection_Array1<gp_Pnt> poles(1, 2);
      NCollection_Array1<double> knots(1, 2);
      NCollection_Array1<int>    mults(1, 2);
      poles(1) = sections[0]->Pole(a);
      poles(2) = sections[1]->Pole(a);
      knots(1) = 0;
      knots(2) = 1;
      mults(1) = mults(2)    = 2;
      columns[size_t(a - 1)] = new Geom_BSplineCurve(poles, knots, mults, 1);
      continue;
    }
    Handle(NCollection_HArray1<gp_Pnt>) points = new NCollection_HArray1<gp_Pnt>(1, int(np));
    Handle(NCollection_HArray1<double>) params = new NCollection_HArray1<double>(1, int(np));
    for (size_t i = 0; i < np; ++i)
    {
      points->SetValue(int(i) + 1, sections[i]->Pole(a));
      params->SetValue(int(i) + 1, across[i]);
    }
    GeomAPI_Interpolate through(points, params, false, Precision::Confusion());
    through.Perform();
    if (!through.IsDone())
      fail(op, no_solution, "the profiles could not be skinned");
    columns[size_t(a - 1)] = through.Curve();
  }
  GeomAlgoExt::unify(columns);
  const Handle(Geom_BSplineCurve)& column = columns[0];
  const int                        nv     = column->NbPoles();
  NCollection_Array2<gp_Pnt>       skinPoles(1, nu, 1, nv);
  for (int a = 1; a <= nu; ++a)
    for (int b = 1; b <= nv; ++b)
      skinPoles(a, b) = columns[size_t(a - 1)]->Pole(b);
  Handle(Geom_BSplineSurface) skin = new Geom_BSplineSurface(skinPoles,
                                                             base->Knots(),
                                                             column->Knots(),
                                                             base->Multiplicities(),
                                                             column->Multiplicities(),
                                                             base->Degree(),
                                                             column->Degree());

  // How far each guide is off the skin along it, with the profiles' places as its parameter.
  std::vector<Handle(Geom_BSplineCurve)> offsets;
  for (size_t k = 0; k < ng; ++k)
  {
    const size_t        j     = order[k];
    const auto&         g     = guides[j];
    const double        first = at[0][j].onSecond, last = at[np - 1][j].onSecond;
    std::vector<double> from;
    for (size_t i = 0; i < np; ++i)
      from.push_back((at[i][j].onSecond - first) / (last - first));
    for (size_t i = 1; i < np; ++i)
      if (from[i] - from[i - 1] <= BRepExt::tolerance::parameter)
        fail(op, no_solution, "a guide does not run from profile to profile in their order");
    std::vector<double> crossings;
    for (size_t i = 0; i < np; ++i)
      crossings.push_back(at[i][j].onSecond);
    std::vector<Handle(Geom_BSplineCurve)> pair{
      refit(g, first, last, Monotone(from, across), crossings, false, tolerance),
      Handle(Geom_BSplineCurve)::DownCast(skin->UIso(nodes[nodeOf[k]]))};
    GeomAlgoExt::unify(pair);
    Handle(Geom_BSplineCurve) off = Handle(Geom_BSplineCurve)::DownCast(pair[0]->Copy());
    for (int b = 1; b <= off->NbPoles(); ++b)
      off->SetPole(b, gp_Pnt(pair[0]->Pole(b).XYZ() - pair[1]->Pole(b).XYZ()));
    offsets.push_back(off);
  }
  GeomAlgoExt::unify(offsets);
  const Handle(Geom_BSplineCurve)& offset = offsets[0];

  // The blend across: on each stretch between two guides, or a guide and an end, the cubic
  // with no slope at either end, as the spline with the stretches' ends as triple knots.
  const int                  stretches = int(nodes.size()) - 1;
  NCollection_Array1<double> hatKnots(1, stretches + 1);
  NCollection_Array1<int>    hatMults(1, stretches + 1);
  for (int k = 0; k <= stretches; ++k)
  {
    hatKnots(k + 1) = nodes[size_t(k)];
    hatMults(k + 1) = k == 0 || k == stretches ? 4 : 3;
  }
  const int nh = 3 * stretches + 1;
  // Poles of the k-th guide's blend: one at its node, the same next to it, nothing elsewhere.
  auto hat = [&](size_t k, int pole) {
    const int node =
      (pole + 1) / 3; // the node a pole belongs to: poles 3s and 3s+1 to node s, 3s-1 to s
    const int owner = pole % 3 == 2 ? node : pole / 3;
    return owner == int(nodeOf[k]) ? 1.0 : 0.0;
  };
  NCollection_Array2<gp_Pnt> movePoles(1, nh, 1, offset->NbPoles());
  for (int a = 1; a <= nh; ++a)
    for (int b = 1; b <= offset->NbPoles(); ++b)
    {
      gp_XYZ sum(0, 0, 0);
      for (size_t k = 0; k < ng; ++k)
        sum += offsets[k]->Pole(b).XYZ() * hat(k, a - 1);
      movePoles(a, b) = gp_Pnt(sum);
    }
  Handle(Geom_BSplineSurface) move = new Geom_BSplineSurface(movePoles,
                                                             hatKnots,
                                                             offset->Knots(),
                                                             hatMults,
                                                             offset->Multiplicities(),
                                                             3,
                                                             offset->Degree());
  // The two of one degree with the same knots, and added.
  const int du = std::max(skin->UDegree(), move->UDegree()),
            dv = std::max(skin->VDegree(), move->VDegree());
  skin->IncreaseDegree(du, dv);
  move->IncreaseDegree(du, dv);
  for (const auto& [from, to] : {std::pair{skin, move}, std::pair{move, skin}})
  {
    for (int k = 1; k <= from->NbUKnots(); ++k)
      to->InsertUKnot(from->UKnot(k), std::min(from->UMultiplicity(k), du), 1e-12, false);
    for (int k = 1; k <= from->NbVKnots(); ++k)
      to->InsertVKnot(from->VKnot(k), std::min(from->VMultiplicity(k), dv), 1e-12, false);
  }
  if (skin->NbUPoles() != move->NbUPoles() || skin->NbVPoles() != move->NbVPoles())
    fail(op, no_solution, "the skin and the guides' move could not be put together");
  for (int a = 1; a <= skin->NbUPoles(); ++a)
    for (int b = 1; b <= skin->NbVPoles(); ++b)
      skin->SetPole(a, b, gp_Pnt(skin->Pole(a, b).XYZ() + move->Pole(a, b).XYZ()));

  std::vector<double> profileWay(np, 1.0), guideWay(ng, 1.0);
  certify(skin, profiles, guides, at, profileWay, guideWay, tolerance);
  BRepBuilderAPI_MakeFace face(skin, Precision::Confusion());
  if (!face.IsDone())
    fail(op, no_solution, "the surface through the net could not be made a face");
  if (!BRepCheck_Analyzer(face.Face()).IsValid())
    fail(op, no_solution, "the skin through the net is not a valid body");
  return stored(kernel, face.Face());
}

namespace
{
//! Curves ready to be brought to one basis exactly: those with weights that differ have no such
//! basis, so they are made again by length, without weights, on one grid; those without weights
//! keep their own parameter. False where one cannot be made again (a corner along it).
bool withoutUnlikeWeights(std::vector<Handle(Geom_BSplineCurve)>& curves, bool periodic)
{
  std::vector<size_t> weighted;
  for (size_t i = 0; i < curves.size(); ++i)
    if (curves[i]->IsRational())
      weighted.push_back(i);
  if (weighted.empty() || std::all_of(curves.begin(), curves.end(), [&](const auto& c) {
        return runAlike(curves[0], c);
      }))
    return true;
  std::vector<Handle(Geom_BSplineCurve)> those;
  for (const size_t i : weighted)
    those.push_back(curves[i]);
  const std::vector<Monotone> byLength(those.size(), Monotone({0.0, 1.0}, {0.0, 1.0}));
  const auto                  made = onOneGrid(those, byLength, {}, periodic);
  if (made.empty())
    return false;
  for (size_t k = 0; k < weighted.size(); ++k)
  {
    // The curve made again is near the one it replaces, not equal to it, and neighbouring edges of
    // a section meet at its vertices: left alone, their ends miss each other by the fit's error
    // (1e-6 on circles beside lines, the size of the sewing tolerance) and the faces do not close.
    // The end poles of a clamped curve are its end points, which are exact on the one it replaces.
    const Handle(Geom_BSplineCurve)& old = curves[weighted[k]];
    if (!made[k]->IsPeriodic() && !old->IsPeriodic())
    {
      made[k]->SetPole(1, old->Pole(1));
      made[k]->SetPole(made[k]->NbPoles(), old->Pole(old->NbPoles()));
    }
    curves[weighted[k]] = made[k];
  }
  return true;
}

//! Closed rings of one edge each, made fit to be joined point to point. Their seams go in one line,
//! as a convention of this code: the first keeps its seam, each next has its seam at its point
//! nearest to the seam before, carried over by the difference of the centres of the two rings'
//! world boxes. This convention is not a theorem: it is pose dependent and greedy from the first
//! ring. The alternative with a proof is cyclic alignment, minimising the squared mismatch over
//! shifts as a cyclic cross-correlation (Dogan, Bernal and Hagwood 2025), which is globally optimal
//! pairwise but moves the seam elsewhere than this convention does. Rational rings that then do not
//! run alike are made again with the share of their length as their parameter, so that the same
//! share of each is joined and rings of one circle make a cylinder. Other sections are left to the
//! skinning, which matches their vertices.
std::vector<TopoDS_Wire> ringsInLine(std::vector<TopoDS_Wire> wires)
{
  const auto soleEdge = [](const TopoDS_Wire& wire) {
    TopoDS_Edge edge;
    int         count = 0;
    for (TopExp_Explorer it(wire, TopAbs_EDGE); it.More(); it.Next(), ++count)
      edge = TopoDS::Edge(it.Current());
    return count == 1 && BRep_Tool::IsClosed(edge) && !BRep_Tool::Degenerated(edge) ? edge
                                                                                    : TopoDS_Edge();
  };
  if (wires.size() < 2 || std::any_of(wires.begin(), wires.end(), [&](const auto& w) {
        return soleEdge(w).IsNull();
      }))
    return wires;
  const auto centre = [](const TopoDS_Wire& wire) {
    Bnd_Box box;
    BRepBndLib::AddOptimal(wire, box, false, false);
    double x0, y0, z0, x1, y1, z1;
    box.Get(x0, y0, z0, x1, y1, z1);
    return gp_XYZ((x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2);
  };
  const auto startOf = [](const TopoDS_Wire& wire) {
    return curveOf(wire, Precision::Confusion())->StartPoint();
  };
  // A skin that comes back to its first ring ends on that ring as it is.
  const bool   loops = wires.size() > 2 && wires.front().IsSame(wires.back());
  const size_t rings = wires.size() - (loops ? 1 : 0);
  gp_Pnt       seam  = startOf(wires[0]);
  gp_XYZ       from  = centre(wires[0]);
  for (size_t i = 1; i < rings; ++i)
  {
    const TopoDS_Edge           edge = soleEdge(wires[i]);
    const gp_XYZ                to   = centre(wires[i]);
    const gp_Pnt                target(seam.XYZ() - from + to);
    double                      first, last;
    const Handle(Geom_Curve)    curve = BRep_Tool::Curve(edge, first, last);
    GeomAPI_ProjectPointOnCurve nearest(target, curve, first, last);
    // The extrema within the range; the ring's own seam, at its ends, is not always among them.
    double at = first;
    if (nearest.NbPoints() > 0 && nearest.LowerDistance() < curve->Value(first).Distance(target))
      at = nearest.LowerDistanceParameter();
    if (curve->Value(at).Distance(curve->Value(first)) > Precision::Confusion())
    {
      TopoDS_Edge moved;
      if (curve->DynamicType() != STANDARD_TYPE(Geom_BSplineCurve) && curve->IsPeriodic())
      {
        // A conic stays a conic: the same curve, from there once round.
        moved = BRepBuilderAPI_MakeEdge(curve, at, at + curve->Period()).Edge();
      }
      else
      {
        Handle(Geom_BSplineCurve) whole =
          GeomConvert::CurveToBSplineCurve(new Geom_TrimmedCurve(curve, first, last));
        if (whole->IsPeriodic())
          whole->SetNotPeriodic();
        moved = BRepBuilderAPI_MakeEdge(reseamed(whole, at, Precision::Confusion())).Edge();
      }
      moved.Orientation(edge.Orientation());
      wires[i] = BRepBuilderAPI_MakeWire(moved).Wire();
    }
    seam = startOf(wires[i]);
    from = to;
  }
  std::vector<Handle(Geom_BSplineCurve)> curves;
  for (size_t i = 0; i < rings; ++i)
    curves.push_back(curveOf(wires[i], Precision::Confusion()));
  const std::vector<Handle(Geom_BSplineCurve)> before = curves;
  if (!withoutUnlikeWeights(curves, true))
    return wires;
  for (size_t i = 0; i < rings; ++i)
    if (curves[i] != before[i])
      wires[i] = BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(curves[i]).Edge()).Wire();
  if (loops)
    wires.back() = wires.front();
  return wires;
}
} // namespace

namespace
{
//! A smooth skin through sections: a tensor product whose every line
//! across the sections is the natural C2 cubic through them, with knots at chord-length places. The
//! law is measured (its poles fit it to four decimals over 2 to 5 rings), not a
//! theorem of the skin's; the natural spline has its own variational theorem. Sections are matched
//! by OpenCascade's compatible wires, a face to each edge; rational sections that do not run alike
//! are made again by length. A section may be a point, which is a row of equal poles. A skin that
//! comes back to its first section closes its line across: `Closed` (periodic form). How far apart
//! two sections are is empirical: for sections of several edges the mean of the distance at the
//! edges' middles (fits squares to 4e-4 of a share), for one edge the mean over it at 64 even
//! shares, an estimate of an integral, not tuned. How the points of rational sections of unlike
//! size or shape are paired is a choice of this code, not derived. Null where a section with
//! weights that differ has a corner along an edge: not made again by length, those stay with
//! OpenCascade's skinning.
TopoDS_Shape skin(const std::vector<TopoDS_Wire>& given,
                  const gp_Pnt*                   start,
                  const gp_Pnt*                   end,
                  bool                            solid)
{
  const char* const skinOp = "skin";
  const double resolution  = 1e-6; // BRepOffsetAPI_ThruSections' default, what it sews and caps to
  NCollection_Sequence<TopoDS_Shape> sections;
  for (const auto& wire : given)
    sections.Append(wire);
  BRepFill_CompatibleWires compatible(sections);
  compatible.Perform();
  if (!compatible.IsDone())
    fail(skinOp, no_solution, "the sections could not be matched edge to edge");
  const auto&  matched = compatible.Shape();
  const bool   loops   = matched.Length() > 2 && matched.First().IsSame(matched.Last());
  const size_t rings   = size_t(matched.Length()) - (loops ? 1 : 0);
  // The curves: `curves[e][s]` of edge e in section s.
  std::vector<std::vector<Handle(Geom_BSplineCurve)>> curves;
  for (size_t s = 0; s < rings; ++s)
  {
    size_t e = 0;
    for (BRepTools_WireExplorer it(TopoDS::Wire(matched(int(s) + 1))); it.More(); it.Next(), ++e)
    {
      if (e == curves.size())
        curves.emplace_back();
      Handle(Geom_BSplineCurve) c = curveOf(it.Current(), Precision::Confusion());
      c->Segment(c->FirstParameter(), c->LastParameter());
      curves[e].push_back(c);
    }
  }
  if (curves.empty())
    fail(skinOp, invalid_input, "a section has no edges");
  for (auto& edge : curves)
  {
    if (edge.size() != rings)
      fail(skinOp, no_solution, "the sections could not be matched edge to edge");
    if (!withoutUnlikeWeights(edge, false))
      return {};
    for (auto& c : edge)
    {
      NCollection_Array1<double> knots(1, c->NbKnots());
      c->Knots(knots);
      BSplCLib::Reparametrize(0, 1, knots);
      c->SetKnots(knots);
    }
    GeomAlgoExt::unify(edge);
  }
  // The places across.
  const auto sampled = [&](size_t s, const std::vector<Handle(Geom_BSplineCurve)>& edge, int k) {
    return edge[s]->Value(k / 64.0);
  };
  std::vector<double> places{0.0};
  const auto          stretch =
    [&](const std::function<gp_Pnt(const std::vector<Handle(Geom_BSplineCurve)>&, int)>& from,
        const std::function<gp_Pnt(const std::vector<Handle(Geom_BSplineCurve)>&, int)>& to) {
      double sum   = 0;
      int    count = 0;
      for (const auto& edge : curves)
        if (curves.size() > 1)
        {
          sum += from(edge, 32).Distance(to(edge, 32));
          ++count;
        }
        else
          for (int k = 0; k < 64; ++k, ++count)
            sum += from(edge, k).Distance(to(edge, k));
      places.push_back(places.back() + sum / count);
    };
  const auto onSection = [&](size_t s) {
    return [&, s](const std::vector<Handle(Geom_BSplineCurve)>& edge, int k) {
      return sampled(s, edge, k);
    };
  };
  const auto atPoint = [](const gp_Pnt& p) {
    return [p](const std::vector<Handle(Geom_BSplineCurve)>&, int) { return p; };
  };
  if (start)
    stretch(atPoint(*start), onSection(0));
  for (size_t s = 1; s < rings; ++s)
    stretch(onSection(s - 1), onSection(s));
  if (end)
    stretch(onSection(rings - 1), atPoint(*end));
  if (loops)
    stretch(onSection(rings - 1), onSection(0));
  for (size_t k = 1; k < places.size(); ++k)
    if (places[k] - places[k - 1] <= Precision::Confusion())
      fail(skinOp, invalid_input, "two sections lie on each other");
  for (double& v : places)
    v /= places.back();
  const GeomAlgoExt::CubicThrough across(places,
                                         loops ? GeomAlgoExt::CubicThrough::Ends::Closed
                                               : GeomAlgoExt::CubicThrough::Ends::Free);
  // A face to each edge.
  BRepBuilderAPI_Sewing    sewing(resolution);
  std::vector<TopoDS_Face> faces;
  for (const auto& edge : curves)
  {
    const auto&                base = edge[0];
    const int                  nu   = base->NbPoles();
    NCollection_Array2<gp_Pnt> poles(1, nu, 1, across.poles);
    NCollection_Array2<double> weights(1, nu, 1, across.poles);
    for (int a = 1; a <= nu; ++a)
    {
      std::vector<gp_Pnt> values;
      if (start)
        values.push_back(*start);
      for (const auto& c : edge)
        values.push_back(c->Pole(a));
      if (end)
        values.push_back(*end);
      if (loops)
        values.push_back(edge[0]->Pole(a));
      const auto column = across.through(values);
      for (int b = 1; b <= across.poles; ++b)
      {
        poles(a, b)   = column[size_t(b - 1)];
        weights(a, b) = base->Weight(a);
      }
    }
    Handle(Geom_BSplineSurface) surface = new Geom_BSplineSurface(poles,
                                                                  weights,
                                                                  base->Knots(),
                                                                  across.knots,
                                                                  base->Multiplicities(),
                                                                  across.mults,
                                                                  base->Degree(),
                                                                  3);
    BRepBuilderAPI_MakeFace     face(surface, Precision::Confusion());
    if (!face.IsDone())
      fail(skinOp, no_solution, "the surface through the sections could not be made a face");
    faces.push_back(face.Face());
    sewing.Add(face.Face());
  }
  // The ends, flat and closed, of what is asked to be solid.
  const auto capOf = [&](size_t s) {
    BRepBuilderAPI_MakeWire wire;
    for (const auto& edge : curves)
      wire.Add(BRepBuilderAPI_MakeEdge(edge[s]).Edge());
    BRepBuilderAPI_MakeFace cap(wire.Wire(), true);
    return cap.IsDone() ? cap.Face() : TopoDS_Face();
  };
  const bool capped = solid && !loops;
  if (capped)
  {
    if (!start)
    {
      const auto cap = capOf(0);
      if (cap.IsNull())
        fail(skinOp, no_solution, "a nonplanar end cannot be closed by a planar cap");
      sewing.Add(cap);
    }
    if (!end)
    {
      const auto cap = capOf(rings - 1);
      if (cap.IsNull())
        fail(skinOp, no_solution, "a nonplanar end cannot be closed by a planar cap");
      sewing.Add(cap);
    }
  }
  sewing.Perform();
  TopoDS_Shape    sewn = sewing.SewedShape();
  TopoDS_Shell    shell;
  TopExp_Explorer found(sewn, TopAbs_SHELL);
  if (found.More())
  {
    shell = TopoDS::Shell(found.Current());
    found.Next();
    if (found.More())
      fail(skinOp, no_solution, "the skin contains disconnected shells");
  }
  else
  {
    BRep_Builder builder;
    builder.MakeShell(shell);
    for (TopExp_Explorer it(sewn, TopAbs_FACE); it.More(); it.Next())
      builder.Add(shell, it.Current());
  }
  if (!solid)
    return shell;
  if (BRepCheck_Shell(shell).Closed() != BRepCheck_NoError)
    fail(skinOp, no_solution, "the skin and its ends leave an open shell");
  BRepBuilderAPI_MakeSolid made(shell);
  if (!made.IsDone())
    fail(skinOp, no_solution, "the skin and its ends do not close into a solid");
  TopoDS_Solid body = made.Solid();
  if (!BRepLib::OrientClosedSolid(body) || !BRepCheck_Analyzer(body).IsValid())
    fail(skinOp, no_solution, "the skin does not bound a valid oriented solid");
  return body;
}
} // namespace

// A skin through sections: smooth by `skin`; ruled, or of two sections, by OpenCascade's.
EvolutionData skinThrough(OcctKernel&                  kernel,
                          const std::vector<uint32_t>& ids,
                          uint32_t                     startVertex,
                          uint32_t                     endVertex,
                          bool                         solid,
                          bool                         ruled)
{
  std::vector<TopoDS_Wire> wires;
  for (const uint32_t id : ids)
    wires.push_back(TopoDS::Wire(get(kernel, id)));
  if (!ruled && wires.size() + (startVertex ? 1 : 0) + (endVertex ? 1 : 0) >= 3)
  {
    gp_Pnt start, end;
    if (startVertex)
      start = BRep_Tool::Pnt(TopoDS::Vertex(get(kernel, startVertex)));
    if (endVertex)
      end = BRep_Tool::Pnt(TopoDS::Vertex(get(kernel, endVertex)));
    const TopoDS_Shape skinned =
      skin(ringsInLine(wires), startVertex ? &start : nullptr, endVertex ? &end : nullptr, solid);
    if (!skinned.IsNull())
    {
      if (!BRepCheck_Analyzer(skinned).IsValid())
        fail("skin", no_solution, "the skin is invalid");
      return stored(kernel, skinned);
    }
  }
  BRepOffsetAPI_ThruSections maker(solid, ruled);
  if (startVertex != 0)
    maker.AddVertex(TopoDS::Vertex(get(kernel, startVertex)));
  for (const auto& wire : ringsInLine(wires))
    maker.AddWire(wire);
  if (endVertex != 0)
    maker.AddVertex(TopoDS::Vertex(get(kernel, endVertex)));
  maker.Build();
  if (!maker.IsDone())
    fail("skin", no_solution, "operation failed");
  const auto result = maker.Shape();
  if (result.IsNull() || !BRepCheck_Analyzer(result).IsValid())
    fail("skin", no_solution, "the skin is invalid");
  if (solid)
  {
    for (TopExp_Explorer shells(result, TopAbs_SHELL); shells.More(); shells.Next())
      if (BRepCheck_Shell(TopoDS::Shell(shells.Current())).Closed() != BRepCheck_NoError)
        fail("skin", no_solution, "the skin leaves an open shell");
  }
  return stored(kernel, result);
}
} // namespace BRepExt
