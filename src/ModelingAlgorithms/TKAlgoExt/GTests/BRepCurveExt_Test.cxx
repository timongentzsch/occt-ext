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
// circle, shape-preserving degree elevation, the cells of overlapping and nested loops, and the
// outline of a sphere and a cylinder, and a curve carried from one plane to another.

#include "BRepExtTest_Support.hxx"

#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_Failure.hxx>

#include <BRepAdaptor_CompCurve.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <Geom_Plane.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <BRep_Tool.hxx>
#include <TopoDS.hxx>
#include <GeomAbs_CurveType.hxx>
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
TEST(BRepCurveExt, RoundWireCornersOfARectangleHasTheClosedFormLength)
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
  const auto          rounded = BRepExt::roundWireCorners(session(), id, vertices, r, 1e-7);
  const TopoDS_Shape& arcs    = BRepExt::get(session(), rounded.resultId);
  EXPECT_NEAR(length(arcs), P - 8 * r + 2 * M_PI * r, 1e-6);
  int edges = 0;
  for (TopExp_Explorer it(arcs, TopAbs_EDGE); it.More(); it.Next())
    ++edges;
  EXPECT_EQ(edges, 8);
  const auto cut = BRepExt::roundWireCorners(session(), id, vertices, -r, 1e-7);
  EXPECT_NEAR(length(BRepExt::get(session(), cut.resultId)),
              P - 8 * r + 4 * std::sqrt(2.) * r,
              1e-6);
}

// A line along the x axis meets, at the origin, an arc of the circle of radius R about (R, 0),
// which leaves the origin upwards and bends to the right. A circle of radius r on the inside of the
// turn touches the line from above and the arc from outside: its centre is (m, r) with (m - R)^2 +
// r^2 = (R + r)^2, that is m = R - sqrt(R^2 + 2 R r), the feet are (m, 0) and (R, 0) + R (M - (R,
// 0)) / (R + r). The wire is trimmed exactly there: the corner is solved for the pieces' own
// geometry, not a polygon through them.
TEST(BRepCurveExt, RoundWireCornersBetweenALineAndAnArcTouchesAtTheClosedFormFeet)
{
  const double R = 10, r = 2;
  const gp_Pnt origin(0, 0, 0), top(R, R, 0),
    middle(R + R * std::cos(3 * M_PI / 4), R * std::sin(3 * M_PI / 4), 0);
  BRepBuilderAPI_MakeWire wire;
  wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(-R, 0, 0), origin).Edge());
  wire.Add(BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(origin, middle, top).Value()).Edge());
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> vertices;
  TopExp::MapShapes(wire.Wire(), TopAbs_VERTEX, vertices);
  uint32_t corner = 0;
  for (int i = 1; i <= vertices.Extent(); ++i)
    if (BRep_Tool::Pnt(TopoDS::Vertex(vertices(i))).Distance(origin) < 1e-9)
      corner = BRepExt::store(session(), vertices(i));
  ASSERT_NE(corner, 0u);
  const uint32_t      id     = BRepExt::store(session(), wire.Wire());
  const auto          result = BRepExt::roundWireCorners(session(), id, {corner}, r, 1e-7);
  const double        m      = R - std::sqrt(R * R + 2 * R * r);
  const gp_Pnt        centre(m, r, 0), bigCentre(R, 0, 0);
  const gp_Pnt        footOnLine(m, 0, 0);
  const gp_Pnt        footOnArc = bigCentre.Translated(gp_Vec(bigCentre, centre) * (R / (R + r)));
  const TopoDS_Shape& made      = BRepExt::get(session(), result.resultId);
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> feet;
  TopExp::MapShapes(made, TopAbs_VERTEX, feet);
  double nearLine = 1e9, nearArc = 1e9;
  for (int i = 1; i <= feet.Extent(); ++i)
  {
    const gp_Pnt p = BRep_Tool::Pnt(TopoDS::Vertex(feet(i)));
    nearLine       = std::min(nearLine, p.Distance(footOnLine));
    nearArc        = std::min(nearArc, p.Distance(footOnArc));
  }
  EXPECT_LT(nearLine, 1e-9);
  EXPECT_LT(nearArc, 1e-9);
  // The rounding arc has the radius asked, about the centre found.
  double arcLength = 0;
  for (TopExp_Explorer it(made, TopAbs_EDGE); it.More(); it.Next())
  {
    BRepAdaptor_Curve curve(TopoDS::Edge(it.Current()));
    if (curve.GetType() == GeomAbs_Circle && std::abs(curve.Circle().Radius() - r) < 1e-9)
    {
      EXPECT_LT(curve.Circle().Location().Distance(centre), 1e-9);
      arcLength += length(it.Current());
    }
  }
  // It sweeps the angle between the radii to the two feet: (0, -1) and (C - M) / (R + r), whose dot
  // product is r / (R + r).
  EXPECT_NEAR(arcLength, r * std::acos(r / (R + r)), 1e-6);
}

// A corner that is too tight for the radius is refused, naming no result, and a vertex given twice
// is an input error: neither may return a body that rounds some other corner.
TEST(BRepCurveExt, RoundWireCornersRefusesATightCornerAndARepeatedVertex)
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
  EXPECT_THROW(BRepExt::roundWireCorners(session(), id, vertices, 6, 1e-7), std::exception);
  std::vector<uint32_t> twice = vertices;
  twice.insert(twice.end(), vertices.begin(), vertices.end());
  EXPECT_THROW(BRepExt::roundWireCorners(session(), id, twice, -1, 1e-7), std::exception);
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

namespace
{

TopoDS_Wire rectangle(double x0, double y0, double x1, double y1)
{
  return polyline({gp_Pnt(x0, y0, 0), gp_Pnt(x1, y0, 0), gp_Pnt(x1, y1, 0), gp_Pnt(x0, y1, 0)},
                  true);
}

std::vector<double> areasOfFaces(const TopoDS_Shape& shape)
{
  std::vector<double> areas;
  for (TopExp_Explorer it(shape, TopAbs_FACE); it.More(); it.Next())
    areas.push_back(area(it.Current()));
  std::sort(areas.begin(), areas.end());
  return areas;
}

} // namespace

// Two overlapping rectangles, [0,6] x [0,4] and [4,10] x [2,6], cut the plane into three bounded
// cells: A alone (24 - 4 = 20), B alone (20) and the overlap [4,6] x [2,4] (4). A cell lost at an
// intersection, or one counted twice, changes the number or the areas at first order.
TEST(BRepCurveExt, RegionsOfTwoOverlappingRectanglesAreTheThreeCells)
{
  const std::vector<uint32_t> curves = {BRepExt::store(session(), rectangle(0, 0, 6, 4)),
                                        BRepExt::store(session(), rectangle(4, 2, 10, 6))};
  const auto result = BRepExt::regions(session(), curves, {0, 0, 0}, {0, 0, 1}, 1e-7);
  const auto areas  = areasOfFaces(BRepExt::get(session(), result.resultId));
  ASSERT_EQ(areas.size(), 3u);
  EXPECT_NEAR(areas[0], 4, 1e-9);
  EXPECT_NEAR(areas[1], 20, 1e-9);
  EXPECT_NEAR(areas[2], 20, 1e-9);
}

// A rectangle with a loop inside it that touches nothing gives two cells, the inner one and the
// ring between: 3 x 2 and 10 x 8 - 3 x 2. Treating the inner loop as a hole of nothing, or
// dropping the ring's hole, gives areas 6 and 80 or 6 and 6.
TEST(BRepCurveExt, RegionsOfANestedLoopAreTheInnerCellAndTheRing)
{
  const std::vector<uint32_t> curves = {BRepExt::store(session(), rectangle(0, 0, 10, 8)),
                                        BRepExt::store(session(), rectangle(3.5, 3, 6.5, 5))};
  const auto result = BRepExt::regions(session(), curves, {0, 0, 0}, {0, 0, 1}, 1e-7);
  const auto areas  = areasOfFaces(BRepExt::get(session(), result.resultId));
  ASSERT_EQ(areas.size(), 2u);
  EXPECT_NEAR(areas[0], 6, 1e-9);
  EXPECT_NEAR(areas[1], 74, 1e-9);
}

// Curves that are not in the plane asked for are refused rather than projected.
TEST(BRepCurveExt, RegionsRefuseCurvesOutsideThePlane)
{
  const std::vector<uint32_t> curves = {BRepExt::store(session(), rectangle(0, 0, 6, 4))};
  EXPECT_THROW(BRepExt::regions(session(), curves, {0, 0, 1}, {0, 0, 1}, 1e-7), std::exception);
}

// The outline of a sphere of radius R seen from any direction is its great circle of radius R
// about the centre, in the plane through the centre normal to the view: length 2 pi R, every
// point at distance R from the centre.
TEST(BRepCurveExt, OutlineOfASphereIsItsGreatCircle)
{
  const double       R = 5;
  const gp_Pnt       centre(1, 2, 3);
  const TopoDS_Shape sphere = BRepPrimAPI_MakeSphere(centre, R).Shape();
  const uint32_t     id     = BRepExt::store(session(), sphere);
  const auto result        = BRepExt::outline(session(), id, {0, 0, 3}, {0, 0, 1}, {1, 0, 0}, true);
  const TopoDS_Shape& made = BRepExt::get(session(), result.resultId);
  EXPECT_NEAR(length(made), 2 * M_PI * R, 1e-6 * 2 * M_PI * R);
  for (TopExp_Explorer it(made, TopAbs_EDGE); it.More(); it.Next())
  {
    BRepAdaptor_Curve edge(TopoDS::Edge(it.Current()));
    for (int k = 0; k <= 8; ++k)
    {
      const gp_Pnt p =
        edge.Value(edge.FirstParameter() + (edge.LastParameter() - edge.FirstParameter()) * k / 8);
      EXPECT_NEAR(p.Distance(centre), R, 1e-7);
      EXPECT_NEAR(p.Z(), centre.Z(), 1e-7);
    }
  }
}

// A cylinder of radius 3 and length 10 along x, seen along z, has the outline of a rectangle
// 10 x 6: the two generators at y = +-3 and the end circles seen edge on. Its box is
// [0,10] x [-3,3] and the two long sides are 10 each.
TEST(BRepCurveExt, OutlineOfACylinderSeenSideOnIsItsTwoGenerators)
{
  const TopoDS_Shape cylinder =
    BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)), 3, 10).Shape();
  const auto          result = BRepExt::outline(session(),
                                       BRepExt::store(session(), cylinder),
                                                {0, 0, 0},
                                                {0, 0, 1},
                                                {1, 0, 0},
                                       true);
  const TopoDS_Shape& made   = BRepExt::get(session(), result.resultId);
  Bnd_Box             box;
  BRepBndLib::Add(made, box);
  double b[6];
  box.Get(b[0], b[1], b[2], b[3], b[4], b[5]);
  EXPECT_NEAR(b[0], 0, 1e-6);
  EXPECT_NEAR(b[3], 10, 1e-6);
  EXPECT_NEAR(b[1], -3, 1e-6);
  EXPECT_NEAR(b[4], 3, 1e-6);
  int generators = 0;
  for (TopExp_Explorer it(made, TopAbs_EDGE); it.More(); it.Next())
    if (std::abs(length(it.Current()) - 10) < 1e-6)
      ++generators;
  EXPECT_EQ(generators, 2);
}

namespace
{

TopoDS_Face squarePlane(double z, double side)
{
  return BRepBuilderAPI_MakeFace(
           new Geom_Plane(gp_Ax3(gp_Pnt(0, 0, z), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0))),
           0,
           side,
           0,
           side,
           1e-7)
    .Face();
}

//! The largest distance from the points of `shape`'s edges to the segment from `a` to `b`, and the
//! farthest the edges' ends reach from them.
double offSegment(const TopoDS_Shape& shape, const gp_Pnt& a, const gp_Pnt& b)
{
  const gp_Lin line(a, gp_Dir(gp_Vec(a, b)));
  double       worst = 0;
  for (TopExp_Explorer it(shape, TopAbs_EDGE); it.More(); it.Next())
  {
    BRepAdaptor_Curve curve(TopoDS::Edge(it.Current()));
    for (int k = 0; k <= 16; ++k)
    {
      const gp_Pnt p     = curve.Value(curve.FirstParameter()
                                   + (curve.LastParameter() - curve.FirstParameter()) * k / 16);
      worst              = std::max(worst, line.Distance(p));
      const double along = gp_Vec(a, p).Dot(gp_Vec(a, b)) / gp_Vec(a, b).SquareMagnitude();
      worst              = std::max(worst, std::max(-along, along - 1.0));
    }
  }
  return worst;
}

} // namespace

// A line from (1, 1, 0) to (9, 1, 3) over a 10 x 10 reference square, carried to a 20 x 20 square
// at height 5 by shares of the parameters: a point (x, y) h above the reference goes to (2x, 2y)
// on the target and h along its normal. The map is linear, so the image is the straight line
// from (2, 2, 5) to (18, 2, 8). Mirroring turns u into 1 - u, so x goes to 20 - 2x; swapping
// exchanges u and v, so the line runs along y; flipping the normal puts the heights below the
// target. Each of these is a change of the shares or of the sign that a wrong order of the
// operations would not give.
TEST(BRepCurveExt, DeformCurveUvnCarriesALineToTheLineOfTheLinearMap)
{
  const uint32_t line =
    BRepExt::store(session(), polyline({gp_Pnt(1, 1, 0), gp_Pnt(9, 1, 3)}, false));
  const uint32_t source = BRepExt::store(session(), squarePlane(0, 10));
  const uint32_t target = BRepExt::store(session(), squarePlane(5, 20));
  const auto     carry  = [&](bool swap, bool mirror, bool flip) {
    return BRepExt::get(session(),
                        BRepExt::deformCurveUVN(session(),
                                                line,
                                                source,
                                                target,
                                                     {1, 1, 1},
                                                     {0, 0, 0},
                                                swap,
                                                mirror,
                                                flip,
                                                1e-7)
                          .resultId);
  };
  EXPECT_LT(offSegment(carry(false, false, false), gp_Pnt(2, 2, 5), gp_Pnt(18, 2, 8)), 1e-6);
  EXPECT_LT(offSegment(carry(false, true, false), gp_Pnt(18, 2, 5), gp_Pnt(2, 2, 8)), 1e-6);
  EXPECT_LT(offSegment(carry(true, false, false), gp_Pnt(2, 2, 5), gp_Pnt(2, 18, 8)), 1e-6);
  EXPECT_LT(offSegment(carry(false, false, true), gp_Pnt(2, 2, 5), gp_Pnt(18, 2, 2)), 1e-6);
}
