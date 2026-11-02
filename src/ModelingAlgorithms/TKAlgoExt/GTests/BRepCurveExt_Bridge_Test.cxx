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

// Bridge curves of geometric continuity G0 to G3 (BRepCurveExt_Bridge.hxx), against what the
// continuity conditions state: the curvature, its change and the torsion of the bridge at an end
// are those of the curve it continues, whatever the tensions, and the tensions move only the
// control points they are the distances of.

#include <BRepCurveExt_Bridge.hxx>

#include <Geom_Circle.hxx>
#include <gp_Ax2.hxx>

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <limits>

namespace
{

constexpr double notGiven = std::numeric_limits<double>::quiet_NaN();

// Curvature, its change with arc length, and torsion of a curve at parameter t from its first three
// derivatives: kappa = |c' x c''| / v^3, tau = (c' x c'') . c''' / |c' x c''|^2 and
// kappa' = [(w . (c' x c''')) / (|w| v^3) - 3 |w| (c' . c'') / v^5] / v with w = c' x c''.
struct Shape
{
  double kappa, change, torsion;
};

Shape shapeOf(const Handle(Geom_Curve)& curve, double t)
{
  gp_Pnt p;
  gp_Vec d1, d2, d3;
  curve->D3(t, p, d1, d2, d3);
  const gp_Vec w = d1.Crossed(d2);
  const double v = d1.Magnitude(), l = w.Magnitude();
  return {l / (v * v * v),
          (w.Dot(d1.Crossed(d3)) / (l * v * v * v) - 3 * l * d1.Dot(d2) / std::pow(v, 5)) / v,
          w.Dot(d3) / (l * l)};
}

// An end that is a point with a direction: order 1.
BRepExt::BridgeEnd withDirection(const gp_Pnt& at, const gp_Vec& leaving)
{
  BRepExt::BridgeEnd end;
  end.point   = at;
  end.tangent = leaving.Normalized();
  end.order   = 1;
  return end;
}

} // namespace

// A cubic G1 bridge from (0, 0) along +x to (3, 1) arriving along +x has the control points
// (0,0), (l,0), (3-l,1), (3,1) for the tensions l, which is point symmetric about (1.5, 0.5)
// whatever l is, and the curvature at the end is (n - 1) h / (n x1^2) = 2 / (3 l^2), h = 1 the
// distance of the third control point from the tangent line.
TEST(BRepCurveExtBridge, CubicG1HasTheClosedFormPolesMidpointAndCurvature)
{
  for (const double l : {0.4, 0.9, 1.3})
  {
    std::array<double, 6> used;
    const auto bridge = BRepExt::bridgeCurve(withDirection(gp_Pnt(0, 0, 0), gp_Vec(1, 0, 0)),
                                             1,
                                             withDirection(gp_Pnt(3, 1, 0), gp_Vec(-1, 0, 0)),
                                             1,
                                             {l, notGiven, notGiven, l, notGiven, notGiven},
                                             used);
    ASSERT_EQ(bridge->Degree(), 3);
    EXPECT_LT(bridge->Pole(2).Distance(gp_Pnt(l, 0, 0)), 1e-14);
    EXPECT_LT(bridge->Pole(3).Distance(gp_Pnt(3 - l, 1, 0)), 1e-14);
    EXPECT_LT(bridge->Value(0.5).Distance(gp_Pnt(1.5, 0.5, 0)), 1e-14);
    EXPECT_NEAR(shapeOf(bridge, 0).kappa, 2 / (3 * l * l), 1e-12) << "l=" << l;
    EXPECT_NEAR(used[0], l, 1e-15);
  }
}

// Continuing a unit circle at both ends, G2 on both: the bridge leaves (1, 0) heading up with the
// circle's curvature 1 and arrives at (0, 1) with it, whatever the tensions. The curvature is the
// end condition, so it is exact to the rounding of the derivatives (1e-12), and so is the tangent.
TEST(BRepCurveExtBridge, G2BridgeOfACircleHasTheCirclesCurvatureAtBothEnds)
{
  Handle(Geom_Circle) circle = new Geom_Circle(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 1.0);
  const auto          first  = BRepExt::bridgeEndOf(circle, 0.0, false);
  const auto          second = BRepExt::bridgeEndOf(circle, M_PI / 2, true);
  for (const std::array<double, 6>& tensions :
       {std::array<double, 6>{0.3, notGiven, notGiven, 0.3, notGiven, notGiven},
        std::array<double, 6>{0.2, 0.5, notGiven, 0.4, 0.9, notGiven}})
  {
    std::array<double, 6> used;
    const auto            bridge = BRepExt::bridgeCurve(first, 2, second, 2, tensions, used);
    ASSERT_EQ(bridge->Degree(), 5);
    EXPECT_NEAR(shapeOf(bridge, 0).kappa, 1.0, 1e-12);
    EXPECT_NEAR(shapeOf(bridge, 1).kappa, 1.0, 1e-12);
    gp_Pnt p;
    gp_Vec d1;
    bridge->D1(0, p, d1);
    EXPECT_LT(gp_Vec(d1.Normalized()).Subtracted(gp_Vec(0, 1, 0)).Magnitude(), 1e-12);
    bridge->D1(1, p, d1);
    EXPECT_LT(gp_Vec(d1.Normalized()).Subtracted(gp_Vec(-1, 0, 0)).Magnitude(), 1e-12);
  }
}

// G3 at the start of a bridge off a helix c(u) = (R cos u, R sin u, b u): curvature R / (R^2 +
// b^2), torsion b / (R^2 + b^2) and no change of curvature, at the end of a degree-4 bridge to a
// point. Changing the three tensions moves the poles but not these three numbers: the condition
// fixes the normal parts and leaves only the distances along the tangent free.
TEST(BRepCurveExtBridge, G3BridgeOfAHelixKeepsCurvatureChangeAndTorsionForAnyTensions)
{
  const double R = 2, b = 0.7, u = 0.9;
  const gp_Pnt at(R * std::cos(u), R * std::sin(u), b * u);
  const auto   first =
    BRepExt::bridgeEndFromDerivatives(at,
                                      gp_Vec(-R * std::sin(u), R * std::cos(u), b),
                                      gp_Vec(-R * std::cos(u), -R * std::sin(u), 0),
                                      gp_Vec(R * std::sin(u), -R * std::cos(u), 0));
  BRepExt::BridgeEnd second;
  second.point                = gp_Pnt(1.0, -2.0, 5.0);
  second.tangent              = gp_Vec(0, 0, 1);
  second.order                = 0;
  const double          kappa = R / (R * R + b * b), tau = b / (R * R + b * b);
  std::array<double, 6> poles0;
  for (const std::array<double, 6>& tensions :
       {std::array<double, 6>{0.3, 0.7, 1.2, notGiven, notGiven, notGiven},
        std::array<double, 6>{0.5, 1.4, 1.9, notGiven, notGiven, notGiven}})
  {
    std::array<double, 6> used;
    const auto            bridge = BRepExt::bridgeCurve(first, 3, second, 0, tensions, used);
    ASSERT_EQ(bridge->Degree(), 4);
    const Shape s = shapeOf(bridge, 0);
    EXPECT_NEAR(s.kappa, kappa, 1e-12);
    EXPECT_NEAR(s.change, 0.0, 1e-11);
    EXPECT_NEAR(s.torsion, tau, 1e-12);
    // The poles that carry the tensions are exactly as far from the end along T as asked.
    const gp_Vec T = first.tangent;
    EXPECT_NEAR(gp_Vec(at, bridge->Pole(2)).Dot(T), tensions[0], 1e-13);
    EXPECT_NEAR(gp_Vec(at, bridge->Pole(3)).Dot(T), tensions[1], 1e-13);
    EXPECT_NEAR(gp_Vec(at, bridge->Pole(4)).Dot(T), tensions[2], 1e-13);
  }
}

// Two collinear ends with a gap, G3 on both: every control point is on the line (the curvatures and
// their change are zero, so the normal parts of the conditions vanish), whatever increasing
// tensions are given, and the curve runs along it without turning back.
TEST(BRepCurveExtBridge, G3BridgeBetweenCollinearEndsIsTheSegment)
{
  BRepExt::BridgeEnd a, b;
  a.point   = gp_Pnt(0, 0, 0);
  a.tangent = gp_Vec(1, 0, 0);
  a.order   = 3;
  b.point   = gp_Pnt(10, 0, 0);
  b.tangent = gp_Vec(-1, 0, 0);
  b.order   = 3;
  std::array<double, 6> used;
  const auto            bridge = BRepExt::bridgeCurve(a, 3, b, 3, {1, 2.5, 4, 1.5, 3, 4.5}, used);
  ASSERT_EQ(bridge->Degree(), 7);
  for (int i = 1; i <= 8; ++i)
  {
    EXPECT_NEAR(bridge->Pole(i).Y(), 0, 1e-13);
    EXPECT_NEAR(bridge->Pole(i).Z(), 0, 1e-13);
  }
  for (int i = 1; i < 8; ++i)
    EXPECT_GT(bridge->Pole(i + 1).X(), bridge->Pole(i).X());
}

// A bridge cannot arrive at a point heading the way it left when the point is where it left from
// along that line: the control polygon (0), (l), (1 + l), (1) on one line goes forward and back, so
// its hodograph has a zero (a cusp). That, ends that coincide, and a tension that would send the
// bridge backwards are refused with an error, never a curve.
TEST(BRepCurveExtBridge, CuspsCoincidentEndsAndBackwardTensionsAreRefused)
{
  std::array<double, 6> used;
  const auto            start = withDirection(gp_Pnt(0, 0, 0), gp_Vec(1, 0, 0));
  EXPECT_THROW(BRepExt::bridgeCurve(start,
                                    1,
                                    withDirection(gp_Pnt(1, 0, 0), gp_Vec(1, 0, 0)),
                                    1,
                                    {0.3, notGiven, notGiven, 0.3, notGiven, notGiven},
                                    used),
               std::exception);
  EXPECT_THROW(BRepExt::bridgeCurve(start,
                                    1,
                                    withDirection(gp_Pnt(0, 0, 0), gp_Vec(0, 1, 0)),
                                    1,
                                    {0.3, notGiven, notGiven, 0.3, notGiven, notGiven},
                                    used),
               std::exception);
  EXPECT_THROW(BRepExt::bridgeCurve(start,
                                    1,
                                    withDirection(gp_Pnt(3, 1, 0), gp_Vec(-1, 0, 0)),
                                    1,
                                    {-0.3, notGiven, notGiven, 0.3, notGiven, notGiven},
                                    used),
               std::exception);
}
