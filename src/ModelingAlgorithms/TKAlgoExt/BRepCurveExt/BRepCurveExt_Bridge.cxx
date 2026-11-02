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

#include <BRepCurveExt_Bridge.hxx>

#include <GeomAlgoExt_Failure.hxx>

#include <NCollection_Array1.hxx>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace BRepExt
{
namespace
{

// The jet of a curve in arc length from its derivatives c', c'', c''' in any regular parameter
// (v = |c'|): T = c' / v, r2 = q / v^2 with q the part of c'' normal to T, and
// r3 = (c''' - (c''' . T) T) / v^3 - 3 (c'' . T) q / v^4 - (q . q) T / v^4. This has no division by
// the curvature, so a straight piece gives r2 = 0 and r3 = 0 without a special case.
void jetFromDerivatives(const gp_Vec& d1,
                        const gp_Vec& d2,
                        const gp_Vec& d3,
                        gp_Vec&       T,
                        gp_Vec&       r2,
                        gp_Vec&       r3)
{
  const double v  = d1.Magnitude();
  T               = d1 / v;
  const gp_Vec q  = d2 - T * d2.Dot(T);
  r2              = q / (v * v);
  const double v4 = v * v * v * v;
  r3 = (d3 - T * d3.Dot(T)) / (v * v * v) - q * (3.0 * d2.Dot(T) / v4) - T * (q.Dot(q) / v4);
}

// The control points past an end A (degree n) for a given jet and tensions x1, x2, x3, numbered
// from the end: result[0] = A. With k the order, only the first k + 1 are set.
std::vector<gp_Pnt> endPoles(const BridgeEnd& end, int k, int n, const double x[3])
{
  std::vector<gp_Pnt> poles{end.point};
  if (k < 1)
    return poles;
  const gp_Vec T = end.tangent;
  // alpha = |B'(0)| = n x1; the rest of the derivatives follow from the three scalars.
  const double alpha = n * x[0];
  poles.push_back(end.point.Translated(T * x[0]));
  if (k < 2)
    return poles;
  const double beta = n * (n - 1) * (x[1] - 2 * x[0]);
  const gp_Vec D2   = end.curvature * (alpha * alpha) + T * beta;
  const gp_Pnt P2   = gp_Pnt(poles[1].XYZ() * 2.0 - end.point.XYZ() + D2.XYZ() / (n * (n - 1.0)));
  poles.push_back(P2);
  if (k < 3)
    return poles;
  // x3 is the T-coordinate of P3, which fixes gamma; the normal part is fixed by the jet.
  const double gamma = n * (n - 1.0) * (n - 2.0) * (x[2] - 3 * x[1] + 3 * x[0])
                       - alpha * alpha * alpha * end.third.Dot(T);
  const gp_Vec D3 =
    end.third * (alpha * alpha * alpha) + end.curvature * (3 * alpha * beta) + T * gamma;
  poles.push_back(gp_Pnt(poles[2].XYZ() * 3.0 - poles[1].XYZ() * 3.0 + end.point.XYZ()
                         + D3.XYZ() / (n * (n - 1.0) * (n - 2.0))));
  return poles;
}

} // namespace

BridgeEnd bridgeEndOf(const Handle(Geom_Curve)& curve, double u, bool against)
{
  if (curve.IsNull())
    fail("bridgeEndOf", invalid_input, "no curve to continue");
  gp_Pnt p;
  gp_Vec d1, d2, d3;
  curve->D3(u, p, d1, d2, d3);
  const double sign = against ? -1.0 : 1.0;
  d1 *= sign;
  d3 *= sign;
  if (d1.Magnitude() < 1e-14)
    fail("bridgeEndOf", unsupported_geometry, "the curve is singular where the bridge begins");
  return bridgeEndFromDerivatives(p, d1, d2, d3);
}

BridgeEnd bridgeEndFromDerivatives(const gp_Pnt& point,
                                   const gp_Vec& d1,
                                   const gp_Vec& d2,
                                   const gp_Vec& d3)
{
  if (d1.Magnitude() < 1e-14)
    fail("bridgeEndFromDerivatives",
         unsupported_geometry,
         "the curve is singular where the bridge begins");
  BridgeEnd end;
  end.point = point;
  jetFromDerivatives(d1, d2, d3, end.tangent, end.curvature, end.third);
  end.order = 3;
  return end;
}

Handle(Geom_BezierCurve) bridgeCurve(const BridgeEnd&             first,
                                     int                          firstOrder,
                                     const BridgeEnd&             second,
                                     int                          secondOrder,
                                     const std::array<double, 6>& tensions,
                                     std::array<double, 6>&       used)
{
  const char* const op    = "bridgeCurve";
  const int         ka    = std::clamp(std::min(firstOrder, first.order), 0, 3);
  const int         kb    = std::clamp(std::min(secondOrder, second.order), 0, 3);
  const int         n     = ka + kb + 1;
  const double      chord = first.point.Distance(second.point);
  if (chord < 1e-12)
    fail(op, invalid_input, "the ends of the bridge coincide");
  double x[6];
  for (int i = 0; i < 6; ++i)
  {
    x[i] = std::isnan(tensions[size_t(i)]) ? (i % 3 + 1) * chord / n : tensions[size_t(i)];
    if (!std::isfinite(x[i]))
      fail(op, invalid_input, "a tension is not a finite number");
    used[size_t(i)] = x[i];
  }
  if ((ka >= 1 && x[0] <= 0) || (kb >= 1 && x[3] <= 0))
    fail(op,
         invalid_input,
         "the first tension of an end must be positive: the bridge would leave backwards");

  const std::vector<gp_Pnt>  a = endPoles(first, ka, n, x);
  const std::vector<gp_Pnt>  b = endPoles(second, kb, n, x + 3);
  NCollection_Array1<gp_Pnt> poles(1, n + 1);
  for (int i = 0; i < n + 1; ++i)
  {
    if (i <= ka)
      poles(i + 1) = a[size_t(i)];
    else if (i >= n - kb)
      poles(i + 1) = b[size_t(n - i)];
    else
      fail(op, no_solution, "the ends share a control point");
  }
  // The ends' own control points cover all n + 1 by n = ka + kb + 1, none shared.
  Handle(Geom_BezierCurve) curve = new Geom_BezierCurve(poles);

  // A bridge that doubles back has a zero of its hodograph, a cusp. The speed is sampled and every
  // sampled local minimum is refined by golden-section search, since a zero lies between samples.
  const auto speed = [&](double t) {
    gp_Pnt p;
    gp_Vec d1;
    curve->D1(t, p, d1);
    return d1.Magnitude();
  };
  constexpr int       samples = 400;
  double              least   = std::numeric_limits<double>::infinity();
  std::vector<double> v(samples + 1);
  for (int i = 0; i <= samples; ++i)
    v[size_t(i)] = speed(double(i) / samples);
  for (int i = 0; i <= samples; ++i)
  {
    least                  = std::min(least, v[size_t(i)]);
    const bool lowerBefore = i == 0 || v[size_t(i)] <= v[size_t(i - 1)];
    const bool lowerAfter  = i == samples || v[size_t(i)] <= v[size_t(i + 1)];
    if (!lowerBefore || !lowerAfter)
      continue;
    double lo = std::max(0.0, double(i - 1) / samples), hi = std::min(1.0, double(i + 1) / samples);
    constexpr double golden = 0.6180339887498949;
    for (int step = 0; step < 60; ++step)
    {
      const double m1 = hi - golden * (hi - lo), m2 = lo + golden * (hi - lo);
      if (speed(m1) < speed(m2))
        hi = m2;
      else
        lo = m1;
    }
    least = std::min(least, speed((lo + hi) / 2));
  }
  if (n > 1 && least < 1e-6 * chord)
    fail(op,
         no_solution,
         "the bridge has a cusp: the ends face each other, or the tensions are too large");
  return curve;
}

} // namespace BRepExt
