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

// Curve operations of BRepCurveExt against closed forms: location by arc length on a polyline and a
// circle, and shape-preserving degree elevation.

#include "BRepExtTest_Support.hxx"

#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_Failure.hxx>

#include <BRepAdaptor_CompCurve.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_Circle.hxx>
#include <NCollection_Array1.hxx>

#include <gtest/gtest.h>

using namespace BRepExtTest;

namespace
{

TopoDS_Wire polyline(const std::vector<gp_Pnt>& points, bool closed)
{
  BRepBuilderAPI_MakeWire maker;
  const size_t            n = points.size();
  for (size_t i = 0; i + (closed ? 0 : 1) < n; ++i)
    maker.Add(BRepBuilderAPI_MakeEdge(points[i], points[(i + 1) % n]).Edge());
  return maker.Wire();
}

} // namespace

// Rounding the four corners of an a x b rectangle by arcs of radius r replaces, at each right-angle
// corner, two stretches of length r (corner to tangent point) by a quarter arc of length pi r / 2:
// P - 8 r + 2 pi r. A negative radius cuts along the chord of the same arc, r sqrt(2) long: P - 8 r
// + 4 sqrt(2) r. Wrong tangent points change the length at first order in r. Eight edges remain:
// four lines and four arcs or chords.
TEST(BRepCurveExt, FilletWireOfARectangleHasTheClosedFormLength)
{
  const double      a = 60, b = 40, r = 5;
  const TopoDS_Wire wire =
    polyline({gp_Pnt(0, 0, 0), gp_Pnt(a, 0, 0), gp_Pnt(a, b, 0), gp_Pnt(0, b, 0)}, true);
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> corner;
  TopExp::MapShapes(wire, TopAbs_VERTEX, corner);
  ASSERT_EQ(corner.Extent(), 4);
  std::vector<uint32_t> vertices;
  for (int i = 1; i <= corner.Extent(); ++i)
    vertices.push_back(BRepExt::store(session(), corner(i)));
  const uint32_t      id      = BRepExt::store(session(), wire);
  const double        P       = 2 * (a + b);
  const auto          rounded = BRepExt::filletWire(session(), id, vertices, r, 1e-7);
  const TopoDS_Shape& arcs    = BRepExt::get(session(), rounded.resultId);
  EXPECT_NEAR(length(arcs), P - 8 * r + 2 * M_PI * r, 1e-6);
  int edges = 0;
  for (TopExp_Explorer it(arcs, TopAbs_EDGE); it.More(); it.Next())
    ++edges;
  EXPECT_EQ(edges, 8);
  const auto cut = BRepExt::filletWire(session(), id, vertices, -r, 1e-7);
  EXPECT_NEAR(length(BRepExt::get(session(), cut.resultId)),
              P - 8 * r + 4 * std::sqrt(2.) * r,
              1e-6);
}

// A corner that is too tight for the radius is refused, naming no result, and a vertex given twice
// is an input error: neither may return a body that rounds some other corner.
TEST(BRepCurveExt, FilletWireRefusesATightCornerAndARepeatedVertex)
{
  const double      a = 10;
  const TopoDS_Wire wire =
    polyline({gp_Pnt(0, 0, 0), gp_Pnt(a, 0, 0), gp_Pnt(a, a, 0), gp_Pnt(0, a, 0)}, true);
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> corner;
  TopExp::MapShapes(wire, TopAbs_VERTEX, corner);
  std::vector<uint32_t> vertices;
  for (int i = 1; i <= corner.Extent(); ++i)
    vertices.push_back(BRepExt::store(session(), corner(i)));
  const uint32_t id = BRepExt::store(session(), wire);
  // r > a / 2: the two arcs at the ends of a side would overlap.
  EXPECT_THROW(BRepExt::filletWire(session(), id, vertices, 6, 1e-7), std::exception);
  std::vector<uint32_t> twice = vertices;
  twice.insert(twice.end(), vertices.begin(), vertices.end());
  EXPECT_THROW(BRepExt::filletWire(session(), id, twice, -1, 1e-7), std::exception);
}

// Locations are by arc length, not by parameter. On a two-segment polyline of lengths 1 and 3 the
// point at a fraction f of the length is at distance 4 f along it, so f = 1/2 is one unit into the
// second segment (parameter-uniform stations would put it at the join), and the tangent is that
// segment's direction.
TEST(BRepCurveExt, LocationsAreByArcLength)
{
  const uint32_t id =
    BRepExt::store(session(), polyline({gp_Pnt(0, 0, 0), gp_Pnt(1, 0, 0), gp_Pnt(1, 3, 0)}, false));
  const auto v = BRepExt::curveLocations(session(), id, {0.0, 0.5, 1.0}, false);
  ASSERT_EQ(v.size(), 36u);
  EXPECT_NEAR(v[0], 0, 1e-9);
  EXPECT_NEAR(v[12], 1, 1e-7);
  EXPECT_NEAR(v[13], 1, 1e-7);
  EXPECT_NEAR(v[15], 0, 1e-7);
  EXPECT_NEAR(v[16], 1, 1e-7);
  EXPECT_NEAR(v[24], 1, 1e-9);
  EXPECT_NEAR(v[25], 3, 1e-9);
}

// On a circle of radius R in the xy-plane at fraction f the point is R (cos 2 pi f, sin 2 pi f),
// the tangent (-sin, cos), the Frenet normal points to the centre (-cos, -sin) and the binormal is
// the plane's normal. The rotation-minimising frame of a planar curve keeps its normal on the
// plane's normal: no twist is added where the curve has none.
TEST(BRepCurveExt, FramesOnACircleAreFrenetAndTwistFree)
{
  const double      R = 7;
  const TopoDS_Wire circle =
    BRepBuilderAPI_MakeWire(
      BRepBuilderAPI_MakeEdge(new Geom_Circle(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), R)).Edge())
      .Wire();
  const uint32_t            id          = BRepExt::store(session(), circle);
  const std::vector<double> fractions   = {0.0, 0.125, 0.3, 0.7};
  const auto                frenet      = BRepExt::curveLocations(session(), id, fractions, false);
  const auto                transported = BRepExt::curveLocations(session(), id, fractions, true);
  for (size_t i = 0; i < fractions.size(); ++i)
  {
    const double  th = 2 * M_PI * fractions[i];
    const double* f  = &frenet[12 * i];
    EXPECT_NEAR(f[0], R * std::cos(th), 1e-6);
    EXPECT_NEAR(f[1], R * std::sin(th), 1e-6);
    EXPECT_NEAR(f[3], -std::sin(th), 1e-6);
    EXPECT_NEAR(f[4], std::cos(th), 1e-6);
    EXPECT_NEAR(f[6], -std::cos(th), 1e-6);
    EXPECT_NEAR(f[7], -std::sin(th), 1e-6);
    EXPECT_NEAR(std::abs(f[11]), 1, 1e-6);
    const double* g = &transported[12 * i];
    EXPECT_NEAR(std::abs(g[8]), 1, 1e-6) << "fraction " << fractions[i];
    EXPECT_NEAR(g[6], 0, 1e-6);
    EXPECT_NEAR(g[7], 0, 1e-6);
  }
}

// Raising the degree of a cubic Bezier by k re-expresses the same curve: the degree grows by k and
// every point is where it was (degree elevation is exact, so the distance is rounding only).
TEST(BRepCurveExt, DegreeElevationKeepsTheCurve)
{
  NCollection_Array1<gp_Pnt> poles(1, 4);
  poles(1) = gp_Pnt(0, 0, 0);
  poles(2) = gp_Pnt(2, 5, 1);
  poles(3) = gp_Pnt(6, 4, -1);
  poles(4) = gp_Pnt(8, 0, 2);
  const TopoDS_Wire wire =
    BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(new Geom_BezierCurve(poles)).Edge()).Wire();
  const uint32_t              id     = BRepExt::store(session(), wire);
  const auto                  result = BRepExt::refineCurve(session(), id, 0, 2, 1e-9);
  const TopoDS_Shape&         raised = BRepExt::get(session(), result.resultId);
  const BRepAdaptor_CompCurve before(wire), after(TopoDS::Wire(raised));
  for (int k = 0; k <= 40; ++k)
  {
    const gp_Pnt p = before.Value(before.FirstParameter()
                                  + (before.LastParameter() - before.FirstParameter()) * k / 40);
    const gp_Pnt q = after.Value(after.FirstParameter()
                                 + (after.LastParameter() - after.FirstParameter()) * k / 40);
    EXPECT_LT(p.Distance(q), 1e-7) << "k=" << k;
  }
  // Lengths by adaptive quadrature: the fixed rule of BRepGProp differs by 1e-5 between degrees 3
  // and 5.
  const auto arc = [](const BRepAdaptor_CompCurve& curve) {
    return GCPnts_AbscissaPoint::Length(curve,
                                        curve.FirstParameter(),
                                        curve.LastParameter(),
                                        1e-10);
  };
  EXPECT_NEAR(arc(after), arc(before), 1e-8);
}
