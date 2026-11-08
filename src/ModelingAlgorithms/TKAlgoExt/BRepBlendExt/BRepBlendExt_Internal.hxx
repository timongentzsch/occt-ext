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

#pragma once

// The internals of the rolling-ball blend, shared by its modules:
//   BRepBlendExt_Ball.cxx      the ball on the supports and its section (Station, Section)
//   BRepBlendExt_Fit.cxx       the blend approximated from the exact sections (Fitted)
//   BRepBlendExt_Spine.cxx     independent closed offset-intersection station carriers
//   BRepBlendExt_Sheet.cxx     exact sheets and frames that read back bit for bit (Sheet)
//   BRepBlendExt_Split.cxx     inherited support breaks, section frames and joined sheets
//   BRepBlendExt_Topology.cxx  meetings, edges, trims and the stop faces' ends
//   BRepBlendExt_Attach.cxx    the blends attached to the body (Attach)
//   BRepBlendExt_Solver.cxx           the pieces of a request, their joints and splits (blend)

#include <Standard_Macro.hxx>
#include <BRepBlendExt_Request.hxx>
#include <GeomAlgoExt_Timing.hxx>
#include <GeomAlgoExt_Cancel.hxx>
#include <GeomAlgoExt_Interval.hxx>
#include <BRepAlgoExt_Geometry.hxx>
#include <Geom_ElementarySurface.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <BRepAlgoExt_ReadBack.hxx>
#include <GeomAlgoExt_SurfaceJet.hxx>
#include <math_Gauss.hxx>
#include <math_Matrix.hxx>
#include <math_Vector.hxx>
#include <GeomAlgoExt_Failure.hxx>
#ifndef BREPBLENDEXT_POLICY
  #include <BRepBlendExt_Policy.hxx>
  #define BREPBLENDEXT_POLICY BRepExt::blending::NativePolicy
#endif
#include <GeomAlgoExt_ModelTolerance.hxx>
#include <deque>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepTools.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>
#include <BRepBndLib.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <Bnd_Box.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <BRep_Builder.hxx>
#include <BRep_CurveRepresentation.hxx>
#include <BRep_TEdge.hxx>
#include <BRep_Tool.hxx>
#include <BSplCLib.hxx>
#include <Geom2dAPI_InterCurveCurve.hxx>
#include <Geom2dConvert.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_Circle.hxx>
#include <Geom2d_Conic.hxx>
#include <Geom2d_Line.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <GeomAPI.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <GeomAPI_IntCS.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <Geom2dAPI_Interpolate.hxx>
#include <NCollection_HArray1.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomInt_IntSS.hxx>
#include <GeomProjLib.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Conic.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Line.hxx>
#include <Geom_Plane.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_SphericalSurface.hxx>
#include <ElSLib.hxx>
#include <gp_Sphere.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <Geom_Surface.hxx>
#include <Geom_SurfaceOfLinearExtrusion.hxx>
#include <Geom_SurfaceOfRevolution.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <IntRes2d_IntersectionPoint.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <NCollection_DataMap.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <NCollection_IndexedMap.hxx>
#include <NCollection_List.hxx>
#include <NCollection_Map.hxx>
#include <ShapeAnalysis_Edge.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopAbs.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <gce_MakeCirc.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <math_Gauss.hxx>
#include <math_Matrix.hxx>
#include <math_Vector.hxx>
#include <algorithm>
#include <functional>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <type_traits>
#include <optional>

namespace BRepExt::blending::detail
{
using SectionPolicy = BREPBLENDEXT_POLICY;
using rho           = SectionPolicy;

namespace tol = BRepExt::tolerance;

//! The blend's own thresholds, each named for what it bounds (the kernel's shared ones are in
//! GeomAlgoExt_ModelTolerance.hxx). Those marked "x size" are shares of the blend's size
//! (Rail::scale), so that a blend is held to the same relative accuracy at any scale, and on
//! lengths between points no finer than their coordinates' rounding (Rail::resolved,
//! tolerance::resolved).
namespace limit
{
constexpr double together = 1e-6; // x size: places two constructions must share (corners, joins)
constexpr double patch    = 1e-6; // x size: Patch's distance along the loop where a blend folds
constexpr double patchDegrees = 0.05; // degrees: Patch's G1 along that loop
constexpr double agree      = 1e-9;  // relative: computed values that must agree (weights, reaches)
constexpr double exactSheet = 1e-10; // x size: feet on an analytic sheet's lines
constexpr double solved =
  1e-12; // x size: Newton's residual for the ball and its corners (Rail::resolved)
constexpr double degenerate = 1e-14; // a vector too short to have a direction; a search's end
constexpr double equal      = 1e-15; // relative: distances and laws taken as given equal
constexpr double knotMerge  = 1e-6;  // the fit's u: a break this near a uniform knot replaces it
constexpr double marched =
  1e-10; // where two sheets meet, each place solved to this length (or its rounding)
// x size: past an edge's end, how far a ball may be from its rail where the faces meet at a right
// angle or wider; where they meet more sharply, that many times the ball's own distance from the
// rail at the end over sqrt 2 (Stations::walk). A judgement.
constexpr double strayed = 4;
constexpr double bend =
  4e-6; // x 1/size: a G2 fit's normal curvature along a contact against its support's
constexpr double seamBend =
  1e-4; // x 1/size: the second forms across a delivered G2 seam (the G2 matrix's own limit)
constexpr double seamDegrees =
  0.01; // degrees: the normals across a delivered G2 seam (the G2 matrix's own limit)
constexpr double lawShare =
  SectionPolicy::lawShare; // a sheet's break at an inner join against its law's own (sheetJoins;
                           // measured, not derived)
constexpr double section =
  5e-6; // x size: a G2 sheet off its section law after coupling (the G2 matrix's own limit)
constexpr int bendRefinements = 6; // refinements a G2 fit takes for that bend alone
constexpr int replays =
  8; // replays (replayOverrun) one inside another: a pocket's four floor corners take four

//! How near two balls, or two places they make, must be to be one: within the model's length
//! resolution, as a vertex lies on its edges, and within what two constructions share at the
//! blend's size (`together`); a ball off a face by more does not rest on it. Never finer than the
//! rounding of points within `magnitude` of the origin (tolerance::rounding).
inline double one(double scale, double magnitude)
{
  return std::max(tol::rounding(magnitude + scale), std::min(tol::length, together * scale));
}
} // namespace limit

//! Whether a difference of curvature over a span of length `span` is no shape: it moves the span's
//! middle by bend span^2 / 8, here no more than two constructions of a blend of `size` share.
inline bool unresolved(double bend, double span, double size)
{
  return bend * span * span / 8 <= limit::together * size;
}

//! A support of the blend: its face, its surface and the sign that makes the surface's normal
//! point out of the body, and its offset distance.
struct Side
{
  TopoDS_Face          face;
  Handle(Geom_Surface) surface;
  double               outward  = 1;
  double               distance = 0;
  //! Explicit polynomial spans of a spline support, for a continued section branch.
  //! Zero retains the support's ordinary span ownership and boundary extension.
  std::array<int, 2> polynomialSpan = {0, 0};
  //! A cliff: the ball rests on this edge of `face`, its far edge, instead of on the face. The
  //! parameters are then the edge's and an angle about it: the foot is on the edge and the
  //! normal points from it at that angle (the ball's centre is on the edge's tube).
  TopoDS_Edge        cliff;
  Handle(Geom_Curve) cliffCurve;
  gp_Vec cliffReference; // the fixed direction the angle is measured from, made normal to the edge
  //! Whether the face is only the strip between the edge and its far edge, which the blend then
  //! replaces whole; else the face goes on past the blend's ends and keeps the rest (for a cliff
  //! and for a notch).
  bool cliffWhole = true;
  //! A notch: the face beyond the face's far edge, which trims the blend where it runs past that
  //! edge (a law's section on the edge, or the round past both ends), and that far edge.
  TopoDS_Face trim;
  TopoDS_Edge trimEdge;
};

//! A size along an edge (Request::radii): the cubic spline through the places, natural (no
//! curvature at its ends, as the reference modeller's variable blends are, measured) or round an edge that
//! closes periodic. Past the ends the natural one goes on along its end tangents, which keeps it
//! twice differentiable there (its second derivative is zero at the ends).
struct SizeLaw
{
  std::vector<double> places, values, bends; // t_i, f_i and the second derivatives there
  bool                periodic = false;

  //! The spline through (places[i], values[i]): its second derivatives from the continuity of the
  //! first, h_{i-1} m_{i-1} + 2 (h_{i-1} + h_i) m_i + h_i m_{i+1} = 6 (s_i - s_{i-1}), s_i the
  //! chords' slopes, with m = 0 at the ends or, periodic, the conditions wrapped round (values
  //! ending as they begin).
  void solve()
  {
    const int n = int(places.size()), unknowns = periodic ? n - 1 : n;
    bends.assign(n, 0);
    if (n < 3)
      return;
    math_Matrix A(1, unknowns, 1, unknowns, 0);
    math_Vector b(1, unknowns, 0);
    const auto  h     = [&](int i) { return places[i + 1] - places[i]; };
    const auto  slope = [&](int i) { return (values[i + 1] - values[i]) / h(i); };
    for (int i = 0; i < unknowns; ++i)
    {
      if (!periodic && (i == 0 || i == n - 1))
      {
        A(i + 1, i + 1) = 1;
        continue;
      }
      const int before = i == 0 ? n - 2 : i - 1, after = (i + 1) % unknowns;
      A(i + 1, before + 1) += h(before);
      A(i + 1, i + 1) += 2 * (h(before) + h(i));
      A(i + 1, after + 1) += h(i);
      b(i + 1) = 6 * (slope(i) - slope(before));
    }
    math_Gauss  gauss(A);
    math_Vector m(1, unknowns);
    gauss.Solve(b, m);
    for (int i = 0; i < unknowns; ++i)
      bends[i] = m(i + 1);
    if (periodic)
      bends[n - 1] = bends[0];
  }

  //! Enclose every span's extrema. Bernstein derivative signs isolate all stationary
  //! points by subdivision, without an absolute parameter sampling scale. Outward
  //! interval arithmetic includes coefficient and subdivision roundoff.
  GeomAlgoExt::interval::Range bounds() const
  {
    using namespace GeomAlgoExt::interval;
    Range result{values.front(), values.front()};
    using Cubic = std::array<Range, 4>;
    struct Span { Cubic b; double lo, hi; };
    std::vector<Span> pending;
    const auto half = [](Range a, Range b) { return multiply(add(a, b), point(0.5)); };
    const auto include = [&](Range range) {
      result.low = std::min(result.low, range.low);
      result.high = std::max(result.high, range.high);
    };
    for (size_t i = 0; i + 1 < places.size(); ++i)
    {
      const Range h = subtract(point(places[i + 1]), point(places[i]));
      const Range h2 = multiply(h, h), y0 = point(values[i]), y1 = point(values[i + 1]);
      const Range m0 = point(bends[i]), m1 = point(bends[i + 1]);
      const Range delta = subtract(y1, y0);
      const Range left = subtract(delta, dividePositive(
        multiply(h2, add(multiply(point(2), m0), m1)), point(6)));
      const Range right = add(delta, dividePositive(
        multiply(h2, add(m0, multiply(point(2), m1))), point(6)));
      pending.push_back({{y0, add(y0, dividePositive(left, point(3))),
                             subtract(y1, dividePositive(right, point(3))), y1}, 0, 1});
      while (!pending.empty())
      {
        const Span span = pending.back();
        pending.pop_back();
        const auto& b = span.b;
        bool increasing = true, decreasing = true;
        double low = b[0].low, high = b[0].high;
        double commonLow = low, commonHigh = high;
        for (int j = 0; j < 4; ++j)
        {
          if (!valid(b[j])) return invalid();
          low = std::min(low, b[j].low); high = std::max(high, b[j].high);
          commonLow = std::max(commonLow, b[j].low);
          commonHigh = std::min(commonHigh, b[j].high);
          if (j < 3)
          {
            const Range derivative = subtract(b[j + 1], b[j]);
            increasing &= derivative.low >= 0;
            decreasing &= derivative.high <= 0;
          }
        }
        if (increasing || decreasing)
        {
          include({std::min(b[0].low, b[3].low), std::max(b[0].high, b[3].high)});
          continue;
        }
        const double middle = (span.lo + span.hi) / 2;
        if (commonLow <= commonHigh || middle == span.lo || middle == span.hi)
        {
          // The remaining stationary region is enclosed at arithmetic resolution.
          // Keeping its hull remains conservative, including a multiple derivative root.
          include({low, high});
          continue;
        }
        const Range a0 = half(b[0], b[1]), a1 = half(b[1], b[2]), a2 = half(b[2], b[3]);
        const Range c0 = half(a0, a1), c1 = half(a1, a2), mid = half(c0, c1);
        pending.push_back({{mid, c1, a2, b[3]}, middle, span.hi});
        pending.push_back({{b[0], a0, c0, mid}, span.lo, middle});
      }
    }
    return result;
  }

  double value(double t) const
  {
    const int    n     = int(places.size());
    const double first = places.front(), last = places.back();
    if (periodic)
      t = first + std::fmod(std::fmod(t - first, last - first) + (last - first), last - first);
    else if (t < first || t > last)
    {
      // Along the end's tangent: the spline's slope there, its second derivative being zero.
      const int    i  = t < first ? 0 : n - 2;
      const double at = t < first ? first : last, hi = places[i + 1] - places[i];
      const double slope =
        (values[i + 1] - values[i]) / hi
        + (t < first ? -(2 * bends[i] + bends[i + 1]) : bends[i] + 2 * bends[i + 1]) * hi / 6;
      return (t < first ? values.front() : values.back()) + slope * (t - at);
    }
    int i           = int(std::upper_bound(places.begin(), places.end(), t) - places.begin()) - 1;
    i               = std::clamp(i, 0, n - 2);
    const double hi = places[i + 1] - places[i], a = (places[i + 1] - t) / hi, b = 1 - a;
    return a * values[i] + b * values[i + 1]
           + ((a * a * a - a) * bends[i] + (b * b * b - b) * bends[i + 1]) * hi * hi / 6;
  }
};

//! Numerical centre-path carrier; its parameter remains the original edge coordinate.
struct CentreParameter
{
  Handle(Geom_BSplineCurve) curve;
  std::vector<double>       knots, lengths;
  double                    tolerance = 0;
  double                    length(double t) const;
  double                    parameter(double length) const;
};

//! A full round (Law::Full) is a notch of its face between the two edges, `s[1]`, whose face beyond
//! is the second support: that side's surface is the face beyond's, where the ball's foot is, and
//! the ball touches the face between too (`Rail::middle`).

struct Rail
{
  //! Finite supports entirely between the source rail and a crossed smooth boundary.
  std::vector<TopoDS_Face> consumedSupports;
  //! Experimental fixed-contact shoulder homotopy on a closed offset circuit.
  struct ConicWindow
  {
    int cylinder = 0;
    gp_Vec axis, acute, planeNormal;
    gp_Pnt planeOrigin;
    double halfWidth = 0;
  };
  std::optional<ConicWindow> conicWindow;
  //! Where a G2 section law crosses a knot of a spline support that is only C2 there, the law is
  //! C0 in its poles; within `half` of the knot's `value`, in the support's parameter, the
  //! section is that of both continued spans in transition (supportKnotWindows).
  struct KnotWindow
  {
    int    side = 0, parameter = 0, knot = 0;
    double value = 0, half = 0;
  };
  std::shared_ptr<const std::vector<KnotWindow>> knotWindows;
  bool boundedFit = false;
  const char*       op = "blendEdges";
  TopoDS_Edge       edge;
  BRepAdaptor_Curve curve;
  //! Offset-intersection carrier; the source edge remains authoritative for source coordinates.
  Handle(Geom_Curve) spine;
  double spineFirst = 0, spineLast = 0;
  double spineParameterFirst = 0, spineParameterLast = 0;

  Handle(Geom_Curve) terminalCarrier[2];

  void stationD1(double t, gp_Pnt& point, gp_Vec& tangent) const
  {
    for (int end = 0; end < 2; ++end)
      if (!terminalCarrier[end].IsNull()
          && (end == 0 ? t <= terminalCarrier[end]->LastParameter()
                       : t >= terminalCarrier[end]->FirstParameter()))
        return terminalCarrier[end]->D1(t, point, tangent);
    if (spine.IsNull())
      return curve.D1(t, point, tangent);
    const double rate = (spineLast - spineFirst) / (spineParameterLast - spineParameterFirst);
    double parameter = spineFirst + (t - spineParameterFirst) * rate;
    // A closed SSI spline can be clamped rather than periodic. Evaluate its one circuit,
    // never the polynomial continuation of its end span.
    const double begin = spine->FirstParameter(), period = spine->LastParameter() - begin;
    parameter -= period * std::floor((parameter - begin) / period);
    spine->D1(parameter, point, tangent);
    tangent *= rate;
  }

  Side              s[2];
  double sigma = -1; // the centre is sigma times the distance along the outward normal: -1 convex
  double first = 0, last = 0, from = 0, to = 0;
  bool   closed = false; // an edge that closes on itself: no ends, and the blend closes too
  Law    law    = Law::Conic;
  double rho    = 0.5;
  double angle  = 0; // a chamfer's bevel off the reference face, side `reference` (Request::angle)
  int    reference = 0;
  double scale     = 1; // the larger distance, for tolerances that scale with the blend
  double magnitude =
    0; // how far the edge reaches from the origin, which its points' rounding scales with

  //! A tolerance on lengths between points of the blend: `share` of its size, and never finer
  //! than the rounding of points where the blend sits (tolerance::resolved).
  double resolved(double share) const { return tol::resolved(share, scale, magnitude); }

  Side middle; // a full round's face between its supports, which the ball touches too
  //! The ball at `t` is the one whose section's plane passes through the edge's point there, or,
  //! where that point lies past the centres' path's centre of curvature (an acute wedge on a
  //! tight bend: the planes cross before they reach it, and no ball's plane passes through the
  //! points there), through the centre a ball would have between the faces' tangent planes at
  //! it, which lies within a fraction of the ball of its own centre (ballAt).
  bool byCentre = false;
  //! Centre interpolation resolves the source path to the station-solve tolerance.
  //! Its extrapolated parameter curves use the endpoint polynomial evaluator.
  bool centreInterpolation = false;
  //! A split changes the rolling supports, not the original edge's parameter frame.
  std::optional<std::array<Side, 2>> parameterSupports;
  //! A contact transferred over its entire interval has a new rolling support.
  //! Its old parameter guide cannot excuse curvature jumps in the delivered sheet.
  bool wholeSupportContinuation = false;
  //! Whether a conic section is made of two spans, split at its middle: where the faces meet so
  //! sharply somewhere along the edge that it turns by more than 120 degrees (Blend::rails), past
  //! which one rational quadratic's middle weight cos(turn / 2) falls below a half and its apex
  //! runs far off the section, ill-conditioned for the fit. Shared along a chain.
  bool halves = false;

  //! The factor on the distances along the edge, by its parameter (Request::radii); empty where
  //! they keep one size.
  SizeLaw sizes;

  //! The factor on the distances at `t`: 1 where they keep one size.
  double factor(double t) const { return sizes.places.empty() ? 1 : sizes.value(t); }

  //! An end where the faces are tangent to each other (two crossing cylinders of one radius): the
  //! ball rests on both at the edge's vertex, and the blend closes to that point [at first, last].
  bool pinch[2] = {false, false};

  //! Whether the faces are tangent at some place inside the edge for a blend of this size
  //! (tolerance::tangentFor) and turn elsewhere. Where they are tangent to the angular resolution
  //! (touchesAt) the ball rests on both at the edge's point and the blend narrows to it, as at a
  //! pinch: the plane of the ball's section, across n0 x n1, is not resolved there.
  bool touching = false;

  //! Whether the faces' normals at the edge's point at `t` agree to the angular resolution
  //! (tolerance::angle), or for `forSize` are tangent for a blend of this size
  //! (tolerance::tangentFor).
  Standard_EXPORT bool tangentAt(double t, bool forSize = false) const;

  //! The edge's parameter at the fit's u, which runs from 0 to 1 over the carried rail, and back.
  std::shared_ptr<const CentreParameter> centreParameter;

  double t(double u) const
  {
    if (!centreParameter)
      return from + (to - from) * u;
    const double a = centreParameter->length(from), b = centreParameter->length(to);
    return centreParameter->parameter(a + (b - a) * u);
  }

  double u(double t) const
  {
    if (!centreParameter)
      return (t - from) / (to - from);
    const double a = centreParameter->length(from), b = centreParameter->length(to);
    return (centreParameter->length(t) - a) / (b - a);
  }

  //! Whether `t` is an end where the blend closes to a point (pinch).
  bool pinchedAt(double t) const
  {
    const double near = tol::parameterMerge * (last - first);
    return (pinch[0] && std::abs(t - first) <= near) || (pinch[1] && std::abs(t - last) <= near)
           || (touching && tangentAt(t));
  }
};

//! A copy of the rail adaptor with a parameter range that includes [lo, hi]. Keep the
//! underlying spline's knots even when they were outside the original edge trim. Trim() on
//! BRepAdaptor_Curve temporarily changes the receiver, so operate on a local copy.
inline Handle(Adaptor3d_Curve) railRange(const Rail& rail, double lo, double hi)
{
  BRepAdaptor_Curve copy(rail.curve);
  return copy.Trim(std::min(lo, copy.FirstParameter()),
                   std::max(hi, copy.LastParameter()),
                   tol::parameter);
}

//! Arclength along the continued rail, including continuity intervals hidden by its edge trim.
inline double railLength(const Rail& rail, double t0, double t1)
{
  const double lo = std::min(t0, t1), hi = std::max(t0, t1);
  if (rail.curve.IsPeriodic()
      || (lo >= rail.curve.FirstParameter() && hi <= rail.curve.LastParameter()))
    return GCPnts_AbscissaPoint::Length(rail.curve, lo, hi);
  return GCPnts_AbscissaPoint::Length(*railRange(rail, lo, hi), lo, hi, rail.resolved(limit::solved));
}

//! Physical separation of two normalized rail stations used at a corner.
inline double railStationSeparation(const Rail& rail, double a, double b)
{
  return railLength(rail, rail.t(a), rail.t(b));
}

//! The parameter `distance` of arclength from t. A regular curve's length is monotone in
//! each direction, so enlarging the range until it contains that length brackets the inverse.
//! This changes only the adaptor's range: point and derivative evaluation keep the same curve.
inline std::optional<double> railParameterAt(const Rail& rail, double t, double distance)
{
  if (!std::isfinite(t) || !std::isfinite(distance))
    return {};
  if (distance == 0)
    return t;
  const double first = rail.curve.FirstParameter(), last = rail.curve.LastParameter();
  if (rail.curve.IsPeriodic()
      || (t >= first && t <= last
          && std::abs(distance) <= railLength(rail, t, distance < 0 ? first : last)))
  {
    GCPnts_AbscissaPoint along(rail.curve, distance, t);
    return along.IsDone() ? std::optional<double>(along.Parameter()) : std::nullopt;
  }
  gp_Pnt p;
  gp_Vec d;
  rail.curve.D1(t, p, d);
  if (!(d.Magnitude() > 0))
    return {};
  double step = distance / d.Magnitude();
  for (;; step *= 2)
  {
    const double end = t + step;
    if (!std::isfinite(end) || end == t)
      return {};
    const auto curve = railRange(rail, std::min(t, end), std::max(t, end));
    const double length = GCPnts_AbscissaPoint::Length(*curve, t, end, rail.resolved(limit::solved));
    if (!std::isfinite(length))
      return {};
    if (length < std::abs(distance))
      continue;
    GCPnts_AbscissaPoint along(rail.resolved(limit::solved), *curve, distance, t);
    return along.IsDone() ? std::optional<double>(along.Parameter()) : std::nullopt;
  }
}

//! Raise the multiplicity required by located continuity events.
inline bool raiseFitBreaks(std::map<double, int>& knots, const std::map<double, int>& breaks)
{
  bool more = false;
  for (const auto& [u, m] : breaks)
  {
    // An optional simple fit knot can yield to an event. Existing continuity
    // events remain at their own stations; proximity is not geometric equality.
    for (auto near = knots.lower_bound(u - limit::knotMerge);
         near != knots.end() && near->first <= u + limit::knotMerge;)
    {
      if (near->first != u && near->second == 1)
        near = knots.erase(near);
      else
        ++near;
    }
    int& slot = knots[u];
    if (slot < m)
    {
      slot = m;
      more = true;
    }
  }
  return more;
}

//! A physical arclength step on a rail, expressed in the fit's normalized parameter.
inline std::optional<double> railUAtDistance(const Rail& rail, double t, double distance)
{
  const auto at = railParameterAt(rail, t, distance);
  if (!at)
    return {};
  const double u = rail.u(*at);
  if (!std::isfinite(u))
    return {};
  return u;
}

//! The ball at one place along the edge: its centre and its feet on the supports.
struct Station
{
  double t        = 0;
  double uv[2][2] = {{0, 0}, {0, 0}};
  gp_Pnt foot[2], centre;
  gp_Vec normal[2], along;
  double uvMiddle[2] = {0, 0}; // a full round's foot on the face between (Rail::middle)
  gp_Pnt middle;
};

//! A surface at one place, with its outward unit normal and that normal's derivatives.
struct Jet
{
  gp_Pnt p;
  gp_Vec du, dv, duu, dvv, duv, n, ndu, ndv;
};

Standard_EXPORT bool jetAt(const Side& side, double u, double v, Jet& j);

//! Where a support's finite spline parameters end and run on along their boundary tangents:
//! nonperiodic spline/Bezier surfaces, and spline profiles of extrusions and revolutions.
//! This keeps carried contacts on a G1 extension instead of an unbounded last-span polynomial.
//! Periodic parameters and other surface types continue on their original surfaces.
Standard_EXPORT std::array<std::optional<std::pair<double, double>>, 2> runsOn(
  const Handle(Geom_Surface)& surface);

//! The parameter v where a surface closes to a point (a sphere's pole, a cone's apex, a
//! revolution's profile ending on its axis), where v lies past it: there (u, v) is the point (u +
//! pi, 2 point - v), the surface's other side, its parameters running the other way round.
Standard_EXPORT std::optional<double> pastPoint(const Handle(Geom_Surface)& surface, double v);

//! A support's surface and its derivatives at (u, v), run on past its ends (runsOn).
Standard_EXPORT void   supportD2(const Handle(Geom_Surface)& surface,
                                 double                      u,
                                 double                      v,
                                 gp_Pnt&                     p,
                                 gp_Vec&                     du,
                                 gp_Vec&                     dv,
                                 gp_Vec&                     duu,
                                 gp_Vec&                     dvv,
                                 gp_Vec&                     duv);
Standard_EXPORT gp_Pnt supportPoint(const Handle(Geom_Surface)& surface, double u, double v);

//! A side's foot at its parameters: on its face's surface, or on a cliff's edge.
inline gp_Pnt footAt(const Side& side, double u, double v)
{
  return side.cliffCurve ? side.cliffCurve->Value(u) : supportPoint(side.surface, u, v);
}

//! The angle about a cliff's edge at `w` of the direction from it to `towards`.
Standard_EXPORT double cliffAngle(const Side& side, double w, const gp_Pnt& towards);

//! An N x N system, its right-hand side in the last column, by elimination with partial pivoting;
//! false where it is singular.
template <size_t N>
bool solve(std::array<std::array<double, N + 1>, N> m, std::array<double, N>& x)
{
  for (size_t c = 0; c < N; ++c)
  {
    size_t pivot = c;
    for (size_t r = c + 1; r < N; ++r)
      if (std::abs(m[r][c]) > std::abs(m[pivot][c]))
        pivot = r;
    if (!(std::abs(m[pivot][c]) > std::numeric_limits<double>::min()))
      return false;
    std::swap(m[c], m[pivot]);
    for (size_t r = c + 1; r < N; ++r)
    {
      const double f = m[r][c] / m[c][c];
      for (size_t k = c; k <= N; ++k)
        m[r][k] -= f * m[c][k];
    }
  }
  for (size_t r = N; r-- > 0;)
  {
    double sum = m[r][N];
    for (size_t k = r + 1; k < N; ++k)
      sum -= m[r][k] * x[k];
    x[r] = sum / m[r][r];
  }
  return std::all_of(x.begin(), x.end(), [](double v) { return std::isfinite(v); });
}

//! Newton's method on N unknowns, with the step halved while it does not reduce the residual.
//! `residual` fills the N equations and their Jacobian at `x`. Within the tolerance it goes on
//! with whole steps while they still reduce the residual (at most two, quadratic convergence
//! takes it to rounding), so that the solution does not depend on the seed.
template <size_t N, typename Residual>
bool newton(std::array<double, N>& x, double tolerance, Residual residual)
{
  std::array<std::array<double, N + 1>, N> m{};
  double                                   norm = 0;
  if (!residual(x, m, norm))
    return false;
  int polished = 0;
  for (int iteration = 0; iteration < 50; ++iteration)
  {
    const bool within = norm < tolerance;
    if (within && polished++ == 2)
      return true;
    std::array<double, N> dx{};
    if (!solve<N>(m, dx))
      return within;
    double step     = 1;
    bool   improved = false;
    for (int halving = 0; halving < (within ? 1 : 12) && !improved; ++halving, step /= 2)
    {
      std::array<double, N> trial = x;
      for (size_t i = 0; i < N; ++i)
        trial[i] += step * dx[i];
      std::array<std::array<double, N + 1>, N> mt{};
      double                                   nt = 0;
      if (residual(trial, mt, nt) && nt < norm)
      {
        x        = trial;
        m        = mt;
        norm     = nt;
        improved = true;
      }
    }
    if (!improved)
    {
      // A rejected trial may have changed the residual's cached surface jets.
      // Restore them at the accepted point before reporting convergence.
      return norm < tolerance && residual(x, m, norm) && norm < tolerance;
    }
  }
  return norm < tolerance;
}

//! Where the march of a rail's stations starts: the edge's middle, or where the faces turn most
//! when they are tangent there (Rail::touching), since no ball has a section of its own there.
Standard_EXPORT double seedPlace(const Rail& rail);

//! The ball touching both supports whose section's plane passes through the edge's point at `t`.
Standard_EXPORT bool contactAt(const Rail& rail, double t, Station& s, double scale = 1);

//! The section at a station, from the foot on the first support to the foot on the second.
struct Section
{
  std::vector<gp_Pnt> poles;
  std::vector<double> weights;
  std::vector<double> knots;
  std::vector<int>    multiplicities;
  int                 degree = 2;
  int                 regime =
    0; // which limits a curved G2 section's handles and depth derivative are at (curvedG2)
};

//! A G2 section leaving a split contact along the inherited support seam.
struct SectionFrame
{
  double t, reach; // rail parameter and compact neighborhood's arclength
  int    side;
  gp_Vec tangent;
};

Standard_EXPORT Handle(Geom_BSplineCurve) sectionCurve(const Section& section);

Standard_EXPORT bool sectionAt(const Rail& rail, const Station& s, Section& out);
Standard_EXPORT bool sectionAt(const Rail&                      rail,
                               const Station&                   s,
                               Section&                         out,
                               const std::vector<SectionFrame>& frames);

//! The sections of both support spans beside a knot, each continued over the other's domain
//! with its own station solve and frame.
Standard_EXPORT bool sectionsBesideSupportKnot(const Rail&                      rail,
                                               const Station&                   seed,
                                               int                              side,
                                               int                              parameter,
                                               int                              knot,
                                               Section                          out[2],
                                               const std::vector<SectionFrame>& frames = {});

//! Continue both adjacent support spans, including their station solves and frames,
//! then blend the complete homogeneous section data with a quintic activation.
//! `position` is the transition coordinate in [0,1]; sectionAt takes it from the rail's
//! knot windows (supportKnotWindows), which bound the departure from the exact law.
Standard_EXPORT bool sectionAcrossSupportKnot(const Rail&                      rail,
                                              const Station&                   seed,
                                              int                              side,
                                              int                              parameter,
                                              int                              knot,
                                              double                           position,
                                              Section&                         out,
                                              const std::vector<SectionFrame>& frames = {});

//! A ball resting on several faces at once: its centre and its foot and outward normal on each.
struct Ball
{
  gp_Pnt              centre;
  std::vector<gp_Pnt> feet;
  std::vector<gp_Vec> normals;
};

//! The ball of radius `distance` resting on every one of `faces`, its centre `sigma` times the
//! distance along each one's outward normal from its foot (-1 convex), by Newton's method on the
//! signed distances from `seed`, whose gradients are the faces' normals at the feet; where more
//! than three faces meet, on their least squares, so that a ball is found only where one rests
//! on all of them (within the length resolution, tol::length). Empty where none does.
Standard_EXPORT std::optional<Ball> restingBall(const std::vector<Side>& faces,
                                                const gp_Pnt&            seed,
                                                double                   distance,
                                                double                   sigma,
                                                double                   scale);

//! Select a regular closed branch of the two constant-distance offset supports through the
//! existing middle station. Does not change the source edge or its parameter range.
Standard_EXPORT bool intersectSpine(Rail& rail);
Standard_EXPORT bool prepareConicWindow(Rail& rail);
Standard_EXPORT double conicPhase(const Rail& rail, const Station& station);

//! The stations of a rail, marched from the middle of the edge to both carried ends, and
//! solved again at any parameter from the nearest one found.
class Stations
{
public:
  explicit Stations(const Rail& rail)
      : rail_(rail)
  {
  }

  //! From the middle of the edge, or from `seed` (a station of the piece before, where a
  //! contact crossed onto this piece's face) to the end.
  void march(const std::optional<Station>& seed = {})
  {
    if (seed)
    {
      Station start = *seed;
      if (!contactAt(rail_, start.t, start))
        fail(rail_.op,
             no_solution,
             "the ball does not rest on the next face where its contact crosses onto it",
             {rail_.edge});
      known_[start.t] = start;
      walk(start, rail_.to);
      walk(start, rail_.from);
      return;
    }
    const double middle = seedPlace(rail_);
    Station      start;
    start.t = middle;
    for (int i = 0; i < 2; ++i)
    {
      double     first, last;
      const auto pcurve = BRep_Tool::CurveOnSurface(rail_.edge, rail_.s[i].face, first, last);
      if (pcurve.IsNull())
        fail(rail_.op,
             invalid_input,
             "an edge has no parameter curve on a face beside it",
             {rail_.edge});
      const gp_Pnt2d uv = pcurve->Value(middle);
      start.uv[i][0]    = uv.X();
      start.uv[i][1]    = uv.Y();
    }
    // The ball grows from the edge to its size, so that the first contact is the one beside it;
    // a ball on a cliff, which a small one does not reach, is found at its size from where it
    // rests on the other face alone.
    const int cliff = rail_.s[0].cliffCurve ? 0 : rail_.s[1].cliffCurve ? 1 : -1;
    if (rail_.law == Law::Full)
    {
      // A full round's ball from the face between's foot at the edge, the second support's
      // at the far edge's nearest place, and half the way across.
      double     first, last;
      const auto onMiddle = BRep_Tool::CurveOnSurface(rail_.edge, rail_.middle.face, first, last);
      gp_Pnt     edge;
      rail_.curve.D0(middle, edge);
      double     f, l;
      const auto far = BRep_Tool::Curve(rail_.s[1].trimEdge, f, l);
      const auto onFar =
        BRep_Tool::CurveOnSurface(rail_.s[1].trimEdge, rail_.s[1].trim, first, last);
      if (onMiddle.IsNull() || far.IsNull() || onFar.IsNull())
        fail(rail_.op,
             invalid_input,
             "an edge has no parameter curve on a face beside it",
             {rail_.edge});
      GeomAPI_ProjectPointOnCurve across(edge, far, f, l);
      if (across.NbPoints() == 0)
        fail(rail_.op,
             no_solution,
             "no ball rests on the face between the edges and both faces beyond",
             {rail_.edge});
      const gp_Pnt2d m = onMiddle->Value(middle), b = onFar->Value(across.LowerDistanceParameter());
      start.uvMiddle[0] = m.X();
      start.uvMiddle[1] = m.Y();
      start.uv[1][0]    = b.X();
      start.uv[1][1]    = b.Y();
      start.middle      = edge;
      start.centre      = gp_Pnt((edge.XYZ() + across.NearestPoint().XYZ()) / 2);
      if (!contactAt(rail_, middle, start))
        fail(rail_.op,
             no_solution,
             "no ball rests on the face between the edges and both faces beyond",
             {rail_.edge});
    }
    else if (cliff >= 0)
    {
      // From the other face's foot a ball's size into it, away from the narrow face (for
      // a convex blend; towards it for a concave one), where the ball nearly rests.
      const int other  = 1 - cliff;
      Side      narrow = rail_.s[cliff];
      narrow.cliffCurve.Nullify();
      Jet j, k;
      if (!jetAt(rail_.s[other], start.uv[other][0], start.uv[other][1], j)
          || !jetAt(narrow, start.uv[cliff][0], start.uv[cliff][1], k))
        fail(rail_.op,
             no_solution,
             "no ball of this size rests on both faces beside the edge",
             {rail_.edge});
      gp_Vec into = k.n - j.n * j.n.Dot(k.n);
      if (!(into.Magnitude() > tol::angle))
        fail(rail_.op,
             no_solution,
             "no ball of this size rests on the face beside the edge and on the face's far edge",
             {rail_.edge});
      into *= (rail_.sigma < 0 ? -1.0 : 1.0) * rail_.s[other].distance / into.Magnitude();
      GeomAPI_ProjectPointOnSurf onto(j.p.Translated(into), rail_.s[other].surface);
      if (onto.NbPoints() == 0
          || !jetAt(rail_.s[other],
                    (onto.LowerDistanceParameters(start.uv[other][0], start.uv[other][1]),
                     start.uv[other][0]),
                    start.uv[other][1],
                    j))
        fail(rail_.op,
             no_solution,
             "no ball of this size rests on the face beside the edge and on the face's far edge",
             {rail_.edge});
      const gp_Pnt guess = j.p.Translated(j.n * (rail_.sigma * rail_.s[other].distance));
      start.uv[cliff][0] =
        GeomAPI_ProjectPointOnCurve(guess, rail_.s[cliff].cliffCurve).LowerDistanceParameter();
      // The centre is the foot plus sigma times the distance along the normal: the normal
      // points from the edge towards the centre, or away from it for a convex blend.
      start.uv[cliff][1] =
        cliffAngle(rail_.s[cliff], start.uv[cliff][0], guess) + (rail_.sigma < 0 ? M_PI : 0);
      if (!contactAt(rail_, middle, start))
        fail(rail_.op,
             no_solution,
             "no ball of this size rests on the face beside the edge and on the face's far edge",
             {rail_.edge});
    }
    else
      for (int step = 1; step <= 8; ++step)
        if (!contactAt(rail_, middle, start, step / 8.0))
          fail(rail_.op,
               no_solution,
               "no ball of this size rests on both faces beside the edge",
               {rail_.edge});
    // Each contact on its own face's side of the other face: offsets that meet beyond the

    gp_Pnt edge;
    rail_.curve.D0(middle, edge);
    for (int i = 0; i < 2 && rail_.law != Law::Full; ++i)
      if (!(gp_Vec(edge, start.foot[i]).Dot(start.normal[1 - i]) * rail_.sigma > 0))
        fail(rail_.op,
             no_solution,
             "these offset distances meet one support before the edge",
             {rail_.edge},
             FailureKind::Overflow);
    known_[middle] = start;
    walk(start, rail_.to);
    walk(start, rail_.from);
  }

  //! The stations found so far, by the edge's parameter.
  const std::map<double, Station>& known() const { return known_; }

  //! Where the rail's carry past an end (towards `from` for 0, `to` for 1) stopped short: where
  //! the ball strays from the rail or no longer rolls on, past which the carry is not made.
  std::optional<double> reached(int end) const { return reached_[end]; }

  //! Whether the march failed where the ball stopped rolling on along the edge (no step found
  //! it again, however short): where no ball's plane passes through the edge's points
  //! (Rail::byCentre), unlike where no ball rests at all.
  bool stalled() const { return stalled_; }

  //! The exact station at `t`, from the nearest station the march found: never from one an
  //! earlier call found, so that it does not depend on which stations were asked for first.
  bool at(double t, Station& s)
  {
    for (const auto* map : {&known_, &asked_})
      if (const auto found = map->find(t); found != map->end())
      {
        s = found->second;
        return true;
      }
    auto           after = known_.lower_bound(t);
    const Station& near =
      after == known_.end() ? std::prev(after)->second
      : after == known_.begin()
        ? after->second
        : (t - std::prev(after)->first < after->first - t ? std::prev(after)->second
                                                          : after->second);
    // Not from a blend's end that closes to a point (Rail::pinch), which seeds no other: from
    // the station on the other side.
    const bool     inside = after != known_.end() && after != known_.begin();
    const Station* other  = !inside                   ? nullptr
                            : &near == &after->second ? &std::prev(after)->second
                                                      : &after->second;
    const Station& from   = rail_.pinchedAt(near.t) && other ? *other : near;
    s                     = from;
    // Between two of them, neither a pinch, the feet seeded where they are on the way from one
    // to the other.
    if (inside && !rail_.pinchedAt(std::prev(after)->first) && !rail_.pinchedAt(after->first))
    {
      const Station &a = std::prev(after)->second, &b = after->second;
      const double   f = (t - a.t) / (b.t - a.t);
      for (int i = 0; i < 2; ++i)
      {
        for (int k = 0; k < 2; ++k)
          s.uv[i][k] = a.uv[i][k] + (b.uv[i][k] - a.uv[i][k]) * f;
        s.uvMiddle[i] = a.uvMiddle[i] + (b.uvMiddle[i] - a.uvMiddle[i]) * f;
      }
    }
    // Between two stations the ball found is theirs only where it lies near the line between
    // their centres (as `onward` asks of a march's step: off it by no more than a quarter of
    // how far they are apart); past a bend tighter than the ball another ball's reach is near,
    // and a jump there is walked instead.
    const auto between = [&](const Station& found) {
      if (!inside)
        return true;
      const Station& a = std::prev(after)->second;
      const Station& b = after->second;
      const double   f = (t - a.t) / (b.t - a.t);
      const gp_Pnt   line(a.centre.XYZ() + (b.centre.XYZ() - a.centre.XYZ()) * f);
      return found.centre.Distance(line) <= a.centre.Distance(b.centre) / 4;
    };
    if (!contactAt(rail_, t, s) || !between(s))
    {
      // Walk there in steps from that station instead of jumping, else from the station on
      // its other side (the nearest may be an end the ball barely reaches).
      const auto walkFrom = [&](const Station& start) {
        s = start;
        for (int k = 1; k <= 16; ++k)
          if (!contactAt(rail_, start.t + (t - start.t) * k / 16, s))
            return false;
        return true;
      };
      const Station* second = &from == &near ? other : &near;
      if (!walkFrom(from) && !(second && !rail_.pinchedAt(second->t) && walkFrom(*second)))
        return false;
    }
    asked_[t] = s;
    return true;
  }

private:
  void walk(Station current, double end)
  {
    const double span     = rail_.last - rail_.first;
    double       step     = span / 16 * (end > current.t ? 1 : -1);
    const int    towards  = end > current.t ? 1 : 0;
    Station      previous = current;
    // How far the ball may stray from the rail past the end: at a right angle between the faces
    // its centre is sqrt 2 times its size from the edge, and limit::strayed times its size away
    // it has left the rail. In a wedge of angle psi the centre is r / sin(psi / 2) from the edge,
    // beyond that bound already below 29 degrees, so where the faces meet more sharply than at a
    // right angle the bound keeps its ratio to the ball's own distance at the end (from the last
    // station on the edge).
    double strayed = -1;
    while ((end - current.t) * step > 0)
    {
      checkCancelled(rail_.op);
      const double t = std::abs(end - current.t) < std::abs(step) ? end : current.t + step;
      // Past the edge's end the rail is carried only as far as the ball goes on with it.
      const bool carried = t < rail_.first || t > rail_.last;
      // Where its size changes along the edge, it goes on past the end as the law's last span
      // (Rail::factor), and the ball only while its size stays within half and twice its size at
      // the end (a judgement): not into the nothing or the bulge a cubic's continuation reaches.
      if (carried && !rail_.sizes.places.empty())
      {
        const double atEnd = rail_.factor(t < rail_.first ? rail_.first : rail_.last),
                     here  = rail_.factor(t);
        if (!(here >= atEnd / 2 && here <= 2 * atEnd))
        {
          reached_[towards] = current.t;
          return;
        }
      }
      Station next = current;
      // Predict the feet from the last two stations.
      if (previous.t != current.t)
      {
        const double f = (t - current.t) / (current.t - previous.t);
        for (int i = 0; i < 2; ++i)
          for (int k = 0; k < 2; ++k)
            next.uv[i][k] += (current.uv[i][k] - previous.uv[i][k]) * f;
      }
      const Station guess = next;
      if (contactAt(rail_, t, next) && onward(previous, current, guess, next))
      {
        // In the turn of the station before: a prediction near a pole can be half a turn off.
        for (int i = 0; i < 2; ++i)
          if (!rail_.s[i].cliffCurve && rail_.s[i].surface->IsUPeriodic())
            next.uv[i][0] =
              BRepExt::nearPeriod(next.uv[i][0], current.uv[i][0], rail_.s[i].surface->UPeriod());
        if (carried && strayed < 0)
          strayed = limit::strayed
                    * std::max(rail_.scale,
                               current.centre.Distance(rail_.curve.Value(current.t)) * M_SQRT1_2);
        if (carried && next.centre.Distance(rail_.curve.Value(t)) > strayed)
        {
          reached_[towards] = current.t;
          return;
        }
        known_[t] = next;
        previous  = current;
        current   = next;
        step      = std::clamp(step * 1.5, -span / 8, span / 8);
      }
      else
      {
        step /= 2;
        if (std::abs(step) < tol::parameterMerge * span)
        {
          if (carried && (current.t <= rail_.first || current.t >= rail_.last))
          {
            reached_[towards] = current.t;
            return;
          }
          stalled_ = true;
          fail(rail_.op,
               no_solution,
               "the ball cannot roll the length of the edge at this size",
               {rail_.edge});
        }
      }
    }
  }

  //! Whether `next` is where the ball rolls on to from `current`: found again from halfway. A
  //! long step's prediction can fall into another ball's reach (where the dihedral turns along
  //! the edge, the feet then leap across the faces); the step is then halved. Asked only where
  //! the solve moved the feet far from the prediction `guess`, against how far they went: on a
  //! smooth run the prediction is good to second order.
  bool onward(const Station& previous,
              const Station& current,
              const Station& guess,
              const Station& next) const
  {
    double off = 0, went = 0;
    for (int i = 0; i < 2; ++i)
    {
      off =
        std::max(off, next.foot[i].Distance(footAt(rail_.s[i], guess.uv[i][0], guess.uv[i][1])));
      went = std::max(went, next.foot[i].Distance(current.foot[i]));
    }
    if (previous.t != current.t && off <= went / 4)
      return true;
    const auto predicted = [](const Station& a, const Station& b, double t) {
      Station s = b;
      s.t       = t;
      if (a.t != b.t)
      {
        const double f = (t - b.t) / (b.t - a.t);
        for (int i = 0; i < 2; ++i)
          for (int k = 0; k < 2; ++k)
            s.uv[i][k] += (b.uv[i][k] - a.uv[i][k]) * f;
      }
      return s;
    };
    Station half = predicted(previous, current, (current.t + next.t) / 2);
    if (!contactAt(rail_, half.t, half))
      return false;
    Station again = predicted(current, half, next.t);
    if (!contactAt(rail_, next.t, again))
      return false;
    return again.centre.Distance(next.centre) <= rail_.resolved(limit::together);
  }

  const Rail&               rail_;
  std::map<double, Station> known_;      // the march's stations
  std::map<double, Station> asked_;      // the stations at() found from them
  std::optional<double>     reached_[2]; // [towards from, towards to]: where a carry stopped short
  bool                      stalled_ = false;
};

//! The blend approximated from the exact sections: degree 5 along the rail, the section's own
//! basis across it, interpolated at the Greville abscissae and refined until the exact sections
//! and feet are met within the length tolerance (de Boor; Piegl & Tiller ch. 9).
struct Fitted
{
  Handle(Geom_BSplineSurface) blend;
  Handle(Geom2d_BSplineCurve) contact[2];
};

constexpr int kDegree = 5;

Standard_EXPORT std::vector<double> flatten(const std::map<double, int>& knots);

//! The breaks of the data along the rail (de Boor): where a foot crosses a knot line of a spline
//! support, and where the edge's own curve has a knot. A support C^c there leaves the section's
//! poles C^(c - k), k being the order of the supports' derivatives the section uses (1 for the
//! conic and the bevel, which follow the normals; 2 for G2, which follows the curvature too), so
//! the fit takes a knot of multiplicity 5 - (c - k) at that place. Each as (u, multiplicity).
Standard_EXPORT std::map<double, int> breaksOf(const Rail& rail,
                                               Stations&   stations,
                                               bool        sought = false,
                                               bool        centre = false);

//! How closely a blend is fitted to its exact sections (and the sections of its parts meet).
Standard_EXPORT double fitTolerance(const Rail& rail);

//! The transition windows of a G2 rail across the C2 knots of its spline supports
//! (Rail::KnotWindow), each as wide as the section and contact bounds allow.
Standard_EXPORT std::shared_ptr<const std::vector<Rail::KnotWindow>> supportKnotWindows(
  const Rail& rail,
  Stations&   stations);

// Contact discovery keeps the requested tolerance but uses circular pilot
// sections; the caller must replace these with the requested law before
// attachment.
Standard_EXPORT Fitted fit(Rail& rail, Stations& stations, bool contactPilot = false);
Standard_EXPORT Fitted fit(const Rail&                      rail,
                           Stations&                        stations,
                           bool                             contactPilot,
                           const std::vector<SectionFrame>& frames,
                           const Fitted*                    terminalReference = nullptr,
                           bool                             reparameterizedTerminal = false);

//! With clamped transverse ends, preserveContactJets retains both contact jets
//! through order two. The remaining split residual requires an independent check.
//! `rates` gives each part's centre length against its u at the join, (dl/du, d2l/du2) for
//! the first part's end then the second's start; without it both carry a centre parameter,
//! along which u is proportional to centre length.
Standard_EXPORT void coupleSplitJets(Handle(Geom_BSplineSurface)& a,
                                      Handle(Geom_BSplineSurface)& b,
                                      const Rail& railA, const Rail& railB,
                                      bool preserveContactJets = false,
                                      const std::array<double, 4>* rates = nullptr);
Standard_EXPORT bool defectiveRailKnot(const Rail& rail, const Fitted& fitted);
Standard_EXPORT std::shared_ptr<CentreParameter> centreParameter(
  const Rail& rail, const std::function<gp_Pnt(double)>& pointAt,
  const std::vector<double>& joins = {}, const std::map<double, int>& events = {});

//! Where each contact's speed the way the ball rolls is least, (speed, u) at
//! each place: worked out once for a fitted piece, for both checks below.
struct Dips
{
  const Rail&                                           rail;
  const Fitted&                                         fitted;
  std::optional<std::vector<std::pair<double, double>>> found[2];
  const std::vector<std::pair<double, double>>&         of(int side);
};

//! A contact that runs back against the way the ball rolls (across both supports' normals at the
//! contacts) folds the blend over itself: the ball is larger than a face beside the edge is round,
//! or a flat cap has no room left. Its speed that way, from its own derivatives, is sampled, and
//! each least speed found is refined by golden-section search to the true least, so that a fold
//! shorter than the samples' spacing (a ball just past the bend) is found.
Standard_EXPORT void checkFold(Dips& dips);

//! Whether a contact that checkFold reads (not a cliff's, nor a notch's) runs back anywhere.
Standard_EXPORT bool runsBack(Dips& dips);

//! Where a contact on a flat cap folds back over itself while the other rolls on (a swallowtail of
//! the blend surface). For each fold, in order along the edge, the fit's u at either end of the
//! loop where the contact crosses itself, where the two balls are one; the blend between is cut
//! out and patched. Empty where a fold lies in no such loop. With `open`, a place where the
//! contact runs back but crosses itself nowhere within the piece is put there (from = to = its u)
//! instead: its loop may run on across a joint into the next piece (Blend::sheets).
struct Fold
{
  int    side = 0;
  double from = 0, to = 0;
};

Standard_EXPORT std::vector<Fold> capFolds(Dips& dips, std::vector<Fold>* open = nullptr);

//! The blend face's surface, and where the rail and the contacts lie in its parameters: along
//! the rail the parameter is `alpha + beta u`, u being the fit's, and each contact is the line
//! where the other parameter is `across[side]`. A cylinder runs the rail along its v.
struct Sheet
{
  Handle(Geom_Surface) surface;
  bool                 swapped = false;
  double               alpha = 0, beta = 1;
  double               across[2] = {0, 1};

  double rail(double u) const { return alpha + beta * u; }

  gp_Pnt2d uv(double rail, double across) const
  {
    return swapped ? gp_Pnt2d(across, rail) : gp_Pnt2d(rail, across);
  }

  //! A contact, parameterised along the rail.
  Handle(Geom_Curve) contact(int side) const
  {
    return swapped ? surface->UIso(across[side]) : surface->VIso(across[side]);
  }

  Handle(Geom2d_Curve) contactUV(int side) const
  {
    return new Geom2d_Line(uv(0, across[side]), swapped ? gp_Dir2d(0, 1) : gp_Dir2d(1, 0));
  }

  //! The section at a rail parameter, parameterised across.
  Handle(Geom_Curve) section(double rail) const
  {
    return swapped ? surface->VIso(rail) : surface->UIso(rail);
  }

  Handle(Geom2d_Curve) sectionUV(double rail) const
  {
    return new Geom2d_Line(uv(rail, 0), swapped ? gp_Dir2d(1, 0) : gp_Dir2d(0, 1));
  }
};

//! The blend's sheet as an exact cylinder or torus, where its spine is a line or a circle. The
//! blend is then exact, not approximated. Empty where the spine is neither, or the feet stray from
//! the surface's lines (limit::exactSheet).
Standard_EXPORT std::optional<Sheet> exactSheet(const Rail& rail, Stations& stations);

//! A cylinder or a torus, whose meetings with another OpenCascade finds as exact conics; other
//! sheets (fitted, extruded) that touch at a mitre's inner end are marched (marchedMeeting).
inline bool conical(const Sheet& sheet)
{
  const GeomAbs_SurfaceType type = GeomAdaptor_Surface(sheet.surface).GetType();
  return type == GeomAbs_Cylinder || type == GeomAbs_Torus;
}

//! Along a straight edge between planes every section is the middle one moved along the edge:
//! the blend is that section's extrusion, whatever its law, its contacts exact lines. Empty
//! elsewhere, and where the feet stray from it (limit::exactSheet).
Standard_EXPORT std::optional<Sheet> extrudedSheet(const Rail& rail, Stations& stations);

//! On an exact sheet, the contacts on the supports from the sheet's own (no fit is needed): each
//! laid on its support and carried to the fit's u, checked against the ball's feet. Empty where
//! they cannot be laid exactly; the blend is fitted then.
Standard_EXPORT std::optional<Fitted> exactContacts(const Rail&  rail,
                                                    Stations&    stations,
                                                    const Sheet& sheet);

//! A parameter curve through `at` over [first, last]: quintic, interpolated at its Schoenberg
//! sites, its spans halved until it is within `tolerance` of `at` between them (empty past 1024).
Standard_EXPORT Handle(Geom2d_BSplineCurve) interpolated(const std::function<gp_Pnt2d(double)>& at,
                                                         double first,
                                                         double last,
                                                         double tolerance);

//! A parameter curve of the fit's u carried to the sheet's rail parameter.
Standard_EXPORT Handle(Geom2d_BSplineCurve) alongRail(const Handle(Geom2d_BSplineCurve)& curve,
                                                      const Sheet&                       sheet);

//! Curves of a split edge's parts, each over its own span of the sheet's rail parameter, end to
//! end as one: their poles and knots as they are, so the parameter is kept (joining them to match
//! speeds at the joints would move it off the edge's). Empty where they do not meet.
Standard_EXPORT Handle(Geom2d_BSplineCurve) endToEnd(std::vector<Handle(Geom2d_BSplineCurve)> parts,
                                                     double tolerance);
Standard_EXPORT Handle(Geom_BSplineCurve)   endToEnd(std::vector<Handle(Geom_BSplineCurve)> parts,
                                                     double tolerance);

//! An edge of the body cut back to `at` on its curve, `to` taking the place of its vertex `from`.
Standard_EXPORT TopoDS_Edge trimmed(const TopoDS_Edge&   edge,
                                    const TopoDS_Vertex& from,
                                    const TopoDS_Vertex& to,
                                    double               at);

using EdgeMap = NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>;

//! A face with some of its edges replaced, each running the way the one it replaces ran (or gone,
//! replaced by nothing), and with `bridges` closing the gaps that leaves where vertices were
//! taken out.
//! Several faces where the bridges leave islands of it, each with the holes inside it.
Standard_EXPORT std::vector<TopoDS_Face> rebuiltFaces(const TopoDS_Face&              face,
                                                      const EdgeMap&                  replace,
                                                      const std::vector<TopoDS_Edge>& bridges);

//! The face's surface without its placement, and the placement, as parameter curves are kept.
Standard_EXPORT std::pair<Handle(Geom_Surface), TopLoc_Location> surfaceOf(const TopoDS_Face& face);

//! `p` moved by whole periods of `surface`'s parameters to lie nearest `target`.
Standard_EXPORT gp_Pnt2d nearPeriod(gp_Pnt2d                    p,
                                    const gp_Pnt2d&             target,
                                    const Handle(Geom_Surface)& surface);

//! A parameter curve moved by whole periods of its surface so that it passes `expected` at `t`.
Standard_EXPORT Handle(Geom2d_Curve) aligned(const Handle(Geom2d_Curve)& curve,
                                             double                      t,
                                             const gp_Pnt2d&             expected,
                                             const Handle(Geom_Surface)& surface);

//! A curve's piece between two of its parameters.
Standard_EXPORT Handle(Geom2d_Curve) between(const Handle(Geom2d_Curve)& curve,
                                             double                      from,
                                             double                      to);

//! Where two curves in one surface's parameters meet: the parameter on each, and the whole
//! periods the second was moved by to meet the first there.
struct Meeting
{
  double   first = 0, second = 0;
  gp_Vec2d shift;
};

//! Where `one` meets `two` in `surface`'s parameters, `two` taken also a whole period either way
//! where the surface is periodic: a curve laid a period from another still meets it.
Standard_EXPORT std::vector<Meeting> meetings(const Handle(Geom2d_Curve)& one,
                                              const Handle(Geom2d_Curve)& two,
                                              const Handle(Geom_Surface)& surface);

//! Of `found`, the one nearest `near` on the first curve.
Standard_EXPORT std::optional<Meeting> nearest(const std::vector<Meeting>& found, double near);

//! The place of a curve nearest a point, and how far it is.
struct OnCurve
{
  double parameter = 0, distance = 0;
};

//! Where `curve` comes nearest `p`, its ends included: projection finds only where the distance
//! is stationary, which a meeting line need not be at its end (two surfaces touching there).
Standard_EXPORT std::optional<OnCurve> onCurve(const gp_Pnt& p, const Handle(Geom_Curve)& curve);

//! A curve of an edge on a surface, as the edge keeps it.
struct OnSurface
{
  Handle(Geom2d_Curve) curve;
  Handle(Geom_Surface) surface;
  TopLoc_Location      location;
};

//! An edge along `curve` from `first` at `t0` to `second` at `t1`, either way along it, with its
//! curves on surfaces over the same parameter: it runs the way the parameter rises.
Standard_EXPORT TopoDS_Edge edgeBetween(const Handle(Geom_Curve)&     curve,
                                        const TopoDS_Vertex&          first,
                                        double                        t0,
                                        const TopoDS_Vertex&          second,
                                        double                        t1,
                                        const std::vector<OnSurface>& on,
                                        double                        tolerance);

//! A curve that lies on an elementary `surface` (a corner's ball), laid on it over [lo, hi] at its
//! own parameter: its projection where that meets it exactly (a line of the surface's parameters),
//! else through its points' own parameters, interpolated (the projection approximates it only to
//! about the length tolerance), where that is nearer.
Standard_EXPORT Handle(Geom2d_Curve) laidOn(const Handle(Geom_Curve)&   curve,
                                            double                      lo,
                                            double                      hi,
                                            const Handle(Geom_Surface)& surface);

//! The straight parameter curve from `a` at `t0` to `b` at `t1`.
Standard_EXPORT Handle(Geom2d_Curve) segment2d(const gp_Pnt2d& a,
                                               const gp_Pnt2d& b,
                                               double          t0 = 0,
                                               double          t1 = 1);

//! Where two surfaces meet (GeomInt_IntSS), each branch with its curve on each surface: where the
//! meeting is analytic IntSS gives none, and the branch is laid on the surface over [lo, hi].
class SurfaceMeeting
{
public:
  //! `lay` is the tolerance an analytic branch is laid on the surfaces to (by default `tolerance`).
  SurfaceMeeting(const Handle(Geom_Surface)& one,
                 const Handle(Geom_Surface)& two,
                 double                      tolerance,
                 double                      lay = 0)
      : one_(one),
        two_(two),
        lay_(lay > 0 ? lay : tolerance),
        meet_(untrimmed(one), untrimmed(two), tolerance, true, true, true)
  {
  }

  int count() const { return meet_.IsDone() ? meet_.NbLines() : 0; }

  const Handle(Geom_Curve)& curve(int i) const { return meet_.Line(i); }

  //! On the first surface (1) or the second (2).
  Handle(Geom2d_Curve) on(int i, int surface, double lo, double hi) const
  {
    // An analytic branch (a conic) is laid on the surface to `lay`: GeomInt's own curves on a
    // round's cylinder miss its ellipse by 3e-2 of the radius at 1e-3 mm, as they do not scale.
    const bool analytic = !meet_.Line(i)->IsKind(STANDARD_TYPE(Geom_BSplineCurve));
    if (!analytic && (surface == 1 ? meet_.HasLineOnS1(i) : meet_.HasLineOnS2(i)))
      return surface == 1 ? meet_.LineOnS1(i) : meet_.LineOnS2(i);
    double tolerance =
      std::min(lay_, Precision::PConfusion()); // never looser than GeomProjLib's own
    return GeomProjLib::Curve2d(meet_.Line(i), lo, hi, surface == 1 ? one_ : two_, tolerance);
  }

private:
  //! A trimmed elementary surface as itself: the intersector finds no line of a trimmed cone.
  static Handle(Geom_Surface) untrimmed(const Handle(Geom_Surface)& surface)
  {
    const auto trimmed = Handle(Geom_RectangularTrimmedSurface)::DownCast(surface);
    return !trimmed.IsNull()
               && trimmed->BasisSurface()->IsKind(STANDARD_TYPE(Geom_ElementarySurface))
             ? trimmed->BasisSurface()
             : surface;
  }

  Handle(Geom_Surface) one_, two_;
  double               lay_;
  GeomInt_IntSS        meet_;
};

//! Where a contact crosses an edge of its support: (u on the contact, parameter on the edge).
Standard_EXPORT std::optional<Meeting> crossing(const Handle(Geom2d_Curve)& contact,
                                                double                      from,
                                                double                      to,
                                                double                      near,
                                                const TopoDS_Edge&          edge,
                                                const TopoDS_Face&          support);

//! Whether a point of a face's surface lies off the face, the point first taken by whole periods
//! of the surface to the face's own parameters: a contact that crosses a sphere's or a
//! cylinder's seam stays on its face.
class Inside
{
public:
  explicit Inside(const TopoDS_Face& face)
      : classifier_(face, uvTolerance(face)),
        surface_(BRep_Tool::Surface(face))
  {
    double u0, u1, v0, v1;
    BRepTools::UVBounds(face, u0, u1, v0, v1);
    mid_ = gp_Pnt2d((u0 + u1) / 2, (v0 + v1) / 2);
  }

  TopAbs_State classify(const gp_Pnt2d& p)
  {
    return classifier_.Perform(nearPeriod(p, mid_, surface_));
  }

  bool out(const gp_Pnt2d& p) { return classify(p) == TopAbs_OUT; }

private:
  //! The length tolerance as a distance in the face's parameters: over its fastest surface
  //! speed, sampled at its parameter box's corners and middle (a cylinder of radius 1000 turns
  //! 1000 mm a radian, and the same face of radius 1e-3 turns 1e-3), so that the classifier takes
  //! as on the face the places within the length tolerance of its boundary at any size.
  static double uvTolerance(const TopoDS_Face& face)
  {
    const Handle(Geom_Surface) surface = BRep_Tool::Surface(face);
    double                     u0, u1, v0, v1;
    BRepTools::UVBounds(face, u0, u1, v0, v1);
    double speed = 0;
    for (const double u : {u0, (u0 + u1) / 2, u1})
      for (const double v : {v0, (v0 + v1) / 2, v1})
      {
        gp_Pnt p;
        gp_Vec du, dv;
        surface->D1(u, v, p, du, dv);
        speed = std::max({speed, du.Magnitude(), dv.Magnitude()});
      }
    return speed > 0 ? tol::length / speed : tol::length;
  }

  BRepTopAdaptor_FClass2d classifier_;
  Handle(Geom_Surface)    surface_;
  gp_Pnt2d                mid_;
};

//! Whether the retained UV contact lies in the finite support. Boundary intersections,
//! including the ends of tangential overlaps, partition its constant-membership intervals.
Standard_EXPORT bool contactWithinFace(const TopoDS_Face& face,
                                      const Handle(Geom2d_Curve)& contact,
                                      double from,
                                      double to);

//! Where a contact meets an edge of its support carried on past its vertex `from`, for a stop
//! face that the blend runs past (the stop face is extended to cap the blend). Only a
//! line or a circle on a line or circle of the support, which runs on exactly; the point found
//! must lie on the stop face's surface too. (u on the contact, parameter on the edge beyond its
//! range). It is looked for as far on as the contact itself lies from the vertex.
Standard_EXPORT std::optional<Meeting> beyond(const Handle(Geom2d_Curve)& contact,
                                              double                      from,
                                              double                      to,
                                              double                      near,
                                              const TopoDS_Edge&          edge,
                                              const TopoDS_Vertex&        vertex,
                                              const TopoDS_Face&          support,
                                              const TopoDS_Face&          stop,
                                              double                      tolerance);

//! The same where the edge between the support and the stop face is not a line or circle: its
//! continuation past `vertex` is where the two surfaces meet (GeomInt_IntSS's branch through the
//! vertex), and the contact meets that. The piece from there back to the vertex becomes an edge
//! of both faces.
struct Continuation
{
  double               u    = 0;         // on the contact
  double               from = 0, to = 0; // on `curve`, from the contact to the vertex
  Handle(Geom_Curve)   curve;
  Handle(Geom2d_Curve) onSupport, onStop;
};

Standard_EXPORT std::optional<Continuation> continuation(const Handle(Geom2d_Curve)& contact,
                                                         double                      from,
                                                         double                      to,
                                                         double                      near,
                                                         const TopoDS_Edge&          edge,
                                                         const TopoDS_Vertex&        vertex,
                                                         const TopoDS_Face&          support,
                                                         const TopoDS_Face&          stop,
                                                         double                      tolerance);

//! A body rebuilt around a solid that one blend changed.
Standard_EXPORT TopoDS_Shape withSolid(const TopoDS_Shape& body,
                                       const TopoDS_Shape& old,
                                       const TopoDS_Shape& solid);

//! How an end of a blended edge meets the rest: a face stops it, it runs on into another piece
//! tangent to it, it is mitred against another piece at an angle, it ends in a corner patch with
//! two others, it goes on as the next part of a split edge on the next face, or it ends where its
//! contact on a flat cap folds over itself and a patch fills the loop to the next part (Fold).
enum class Meets
{
  Stop,
  Joint,
  Miter,
  Corner,
  Split,
  Fold
};

//! One edge of a blend, or the part of it on one pair of faces: its rail, fit and sheet, and at
//! each end what it meets: for a joint or a miter the other piece and its end there (and for a
//! joint whether their sides swap), for a corner the corner's index, and for a split the other
//! part of the same edge, where one contact crosses the edge `crossed` onto the next face.
struct Piece
{
  Rail   rail;
  Fitted fitted;
  Sheet  sheet;

  struct Joint
  {
    Meets       meets = Meets::Stop;
    int         piece = -1, end = 0;
    bool        swapped = false;
    int         side    = 0;
    TopoDS_Edge crossed;
  } joint[2];

  struct TerminalContact
  {
    TopoDS_Face support;
    TopoDS_Edge arc;
    double parameter = 0, station = 0;
  };
  std::optional<TerminalContact> terminalContact[2][2]; // [end][side]


  //! A blend between two faces that meet only at this vertex of the body, where more than three
  //! blended edges meet and the corner parts into balls on three faces or more
  //! (Blend::vertexBlend): its rail is the line where the faces' planes meet, not an edge of the
  //! body.
  TopoDS_Vertex apex;
  //! The vertex a face stops the blend at, where it has run on over the next edges of its chain
  //! (Blend::absorbEnds), and those edges.
  TopoDS_Vertex            stopVertex[2];
  std::vector<TopoDS_Edge> swallowed;

  //! The fit's u at an end of the edge itself.
  double uAt(int end) const { return rail.u(end == 0 ? rail.first : rail.last); }

  //! Whether the sheet is exact (analytic or extruded): not the fit, nor the fits of a split
  //! edge's parts joined, which a split chain always is (Blend::chains).
  bool exact() const
  {
    return sheet.surface.get() != fitted.blend.get() && joint[0].meets != Meets::Split
           && joint[1].meets != Meets::Split;
  }
};

//! Whether the seam runs to a pole in the strip the blend replaces, on the edge's side of the
//! contact: it goes with the strip (Attach::sweepToPole), and the contact, which passes the seam's
//! place and comes back, stays in one run of the face's parameters.
Standard_EXPORT bool poleInStrip(const Piece& piece, int side, const TopoDS_Edge& seam);

//! Split parts expressed as one sheet, with inherited support breaks continued internally.
Standard_EXPORT void joinSplitPieces(const char*                                op,
                                     std::vector<Piece>&                        pieces,
                                     const std::vector<std::optional<Station>>& seeds);

//! A face the blends change: its edges replaced only in it, and the edges that close its gaps.
struct Edit
{
  TopoDS_Face              face;
  EdgeMap                  replace;
  std::vector<TopoDS_Edge> bridges;
};

//! The faces to rebuild. A deque: `editOf` hands out references that must stay valid while it
//! adds the next face's edit, which a vector's growth would leave dangling.
using Edits = std::deque<Edit>;

Standard_EXPORT Edit& editOf(Edits& edits, const TopoDS_Shape& face);

//! A surface as an edge's curve on it is kept: as placed in space for computing, and its basis
//! and placement for storing.
struct Placed
{
  Handle(Geom_Surface) surface, basis;
  TopLoc_Location      location;
};

Standard_EXPORT Placed placedSheet(const Sheet& sheet);

Standard_EXPORT Placed placedFace(const TopoDS_Face& face);

//! The edge where two surfaces meet between two corners, with its curve on each: the shorter way
//! round on a closed curve, and each parameter curve moved by whole periods to pass where the
//! first corner is on that surface. Empty where they do not meet through both corners.
Standard_EXPORT std::optional<TopoDS_Edge> meetingEdge(const Placed&        one,
                                                       const gp_Pnt2d&      firstOnOne,
                                                       const Placed&        two,
                                                       const gp_Pnt2d&      firstOnTwo,
                                                       const TopoDS_Vertex& first,
                                                       const TopoDS_Vertex& second,
                                                       double               scale);

//! A sheet/plane meeting in a transverse graph chart, split at the sheet's spline knots.
Standard_EXPORT std::optional<TopoDS_Edge> planeMeeting(const Sheet&         sheet,
                                                        const TopoDS_Face&   face,
                                                        const gp_Pnt2d&      firstUV,
                                                        const gp_Pnt2d&      lastUV,
                                                        const TopoDS_Vertex& first,
                                                        const TopoDS_Vertex& last);

//! Where two sheets meet between two corners where both are known in each sheet's parameters,
//! marched along its length from the first corner, in even steps, each place on the plane across
//! the meeting's direction there, all four parameters by Newton's method from the last place. The
//! two sheets may touch at a corner (a mitre's inner end, where both are tangent to the face they
//! share): that corner is taken as given. The curve and both parameter curves interpolate the
//! same places at their chord lengths, the steps halved until, between the places, the curve
//! opens no gap between the sheets beyond the length tolerance (its distance to the meeting times
//! the sine of their angle) and its parameter curves lie on it.
Standard_EXPORT std::optional<TopoDS_Edge> marchedMeeting(const Sheet&         one,
                                                          const gp_Pnt2d&      oneFirst,
                                                          const gp_Pnt2d&      oneSecond,
                                                          const Sheet&         two,
                                                          const gp_Pnt2d&      twoFirst,
                                                          const gp_Pnt2d&      twoSecond,
                                                          const TopoDS_Vertex& first,
                                                          const TopoDS_Vertex& second);

//! Where two sheets cross between two places where both meet (a crease, where the blends of two
//! parts of an outline run into each other over a narrow neck of their cap), each place known in
//! both sheets' parameters: marched along the chord between them, each place where the plane
//! across the chord meets both sheets, interpolated as marchedMeeting does.
Standard_EXPORT std::optional<TopoDS_Edge> crossedMeeting(const Sheet&         one,
                                                          const gp_Pnt2d&      oneFirst,
                                                          const gp_Pnt2d&      oneSecond,
                                                          const Sheet&         two,
                                                          const gp_Pnt2d&      twoFirst,
                                                          const gp_Pnt2d&      twoSecond,
                                                          const TopoDS_Vertex& first,
                                                          const TopoDS_Vertex& second);

//! The section of a sheet at a rail parameter between two corners, running the way its own
//! across parameter does.
Standard_EXPORT TopoDS_Edge isoSection(const Sheet&         sheet,
                                       double               rail,
                                       const TopoDS_Vertex& first,
                                       const TopoDS_Vertex& second);

//! An edge's unit direction away from a vertex of it, at its first end (0), its last (1) or
//! whichever is there (-1); false where it does not end there or has no direction.
Standard_EXPORT bool leaving(const TopoDS_Edge&   edge,
                             const TopoDS_Vertex& vertex,
                             gp_Vec&              out,
                             int                  end = -1);

//! Two edges that go on from each other at a vertex, tangent within the blend tangency angle.
Standard_EXPORT bool tangent(const gp_Vec& a, const gp_Vec& b);

//! A face's outward unit normal at a vertex of an edge of it; false where it is singular there.
Standard_EXPORT bool normalAt(const TopoDS_Face&   face,
                              const TopoDS_Edge&   edge,
                              const TopoDS_Vertex& vertex,
                              gp_Vec&              out);

//! Whether two faces of an edge go on smoothly into each other at a vertex of it: their outward
//! normals there the same within the blend tangency angle (two faces of one wall).
Standard_EXPORT bool smoothAt(const TopoDS_Face&   a,
                              const TopoDS_Face&   b,
                              const TopoDS_Edge&   edge,
                              const TopoDS_Vertex& vertex);

//! Where three or more blended edges meet: the vertex, and each piece's end there.
struct Corner
{
  TopoDS_Vertex                    vertex;
  std::vector<std::pair<int, int>> ends;
  //! The faces there are tangent to each other: every blend closes to the vertex (Rail::pinch),
  //! and the corner has no patch.
  bool pinch = false;
};

//! What attach refuses where a blend runs into a face of the body beside none of its edges (a
//! hole's wall): the faces it runs into, which blend() may lift off and put back (liftFeatures).
struct Collision : Failure
{
  std::vector<TopoDS_Shape> faces;

  Collision(const char*               op,
            const std::string&        message,
            const TopoDS_Shape&       at,
            std::vector<TopoDS_Shape> faces_)
      : Failure(op, unsupported_geometry, message, {failureHash(at)}),
        faces(std::move(faces_))
  {
    shapes = {at};
  }
};

//! A split contact reaches an unselected tangent edge between its current supports:
//! those faces provide no stop there. The caller may continue the rail on that edge.
struct ContinuationRequired : Failure
{
  ContinuationRequired(const char* op, const TopoDS_Shape& rail, const TopoDS_Shape& next)
      : Failure(op,
                unsupported_geometry,
                "the blend reaches a continuation between its supports",
                {failureHash(rail), failureHash(next)})
  {
    shapes = {rail, next};
  }
};

//! The blends of the pieces attached to the body (BRepBlendExt_Attach.cxx).
Standard_EXPORT TopoDS_Shape attach(TopoDS_Shape                     body,
                                    const std::vector<Piece>&        pieces,
                                    const std::vector<Corner>&       cornersAt,
                                    const Handle(BRepTools_History)& history,
                                    std::vector<Join>* joins = nullptr,
                                    bool refineContacts = false);

} // namespace BRepExt::blending::detail
