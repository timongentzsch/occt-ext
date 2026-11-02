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

// Skins and sweeps of BRepSkinExt against closed forms: stacks of coaxial circles (frusta) and a
// disc swept along a circular arc (Pappus).

#include "BRepExtTest_Support.hxx"

#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_Failure.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>

#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <NCollection_HArray1.hxx>
#include <Geom_Circle.hxx>
#include <Geom_TrimmedCurve.hxx>

#include <gtest/gtest.h>

using namespace BRepExtTest;

namespace
{

TopoDS_Wire circleAt(double z, double radius)
{
  return BRepBuilderAPI_MakeWire(
           BRepBuilderAPI_MakeEdge(
             new Geom_Circle(gp_Ax2(gp_Pnt(0, 0, z), gp_Dir(0, 0, 1)), radius))
             .Edge())
    .Wire();
}

TopoDS_Wire squareAt(double z, double half)
{
  const gp_Pnt            p[4] = {gp_Pnt(-half, -half, z),
                                  gp_Pnt(half, -half, z),
                                  gp_Pnt(half, half, z),
                                  gp_Pnt(-half, half, z)};
  BRepBuilderAPI_MakeWire maker;
  for (int i = 0; i < 4; ++i)
    maker.Add(BRepBuilderAPI_MakeEdge(p[i], p[(i + 1) % 4]).Edge());
  return maker.Wire();
}

double frustum(double height, double r0, double r1)
{
  return M_PI * height / 3 * (r0 * r0 + r0 * r1 + r1 * r1);
}

} // namespace

// A ruled skin through coaxial circles at unequal heights is a stack of frusta, each of volume
// pi h (r0^2 + r0 r1 + r1^2) / 3. A section placed at the wrong height, a smooth blend where ruled
// was asked, or caps missing change the volume at first order.
TEST(BRepSkinExt, RuledSkinThroughCoaxialCirclesIsAStackOfFrusta)
{
  const double          z[3] = {0, 4, 11}, r[3] = {6, 3, 5};
  std::vector<uint32_t> sections;
  for (int i = 0; i < 3; ++i)
    sections.push_back(BRepExt::store(session(), circleAt(z[i], r[i])));
  const auto          result = BRepExt::skinThrough(session(), sections, 0, 0, true, true);
  const TopoDS_Shape& made   = BRepExt::get(session(), result.resultId);
  EXPECT_TRUE(BRepExt::invalidity(made).empty());
  const double expected = frustum(z[1] - z[0], r[0], r[1]) + frustum(z[2] - z[1], r[1], r[2]);
  EXPECT_NEAR(volume(made), expected, 1e-6 * expected);
}

// Through two sections a smooth skin has nothing to bend about and is the ruled one.
TEST(BRepSkinExt, SmoothSkinThroughTwoCirclesIsTheFrustum)
{
  std::vector<uint32_t> sections = {BRepExt::store(session(), circleAt(0, 6)),
                                    BRepExt::store(session(), circleAt(9, 2))};
  const auto            result   = BRepExt::skinThrough(session(), sections, 0, 0, true, false);
  const TopoDS_Shape&   made     = BRepExt::get(session(), result.resultId);
  EXPECT_TRUE(BRepExt::invalidity(made).empty());
  EXPECT_NEAR(volume(made), frustum(9, 6, 2), 1e-6 * frustum(9, 6, 2));
}

// A smooth skin through a square, a circle and a square. The circle's arcs are rational and the
// square's edges are not, so the arcs are made again without weights, and neighbouring arcs must
// still meet at the section's vertices: ends that miss each other by the size of the fit's error
// (1.4e-6 here) leave three pairs of coincident, unsewn edges and a solid that is not connected.
// The sections are exact on the result: every sampled point of them lies on its faces.
TEST(BRepSkinExt, SmoothSkinThroughSquareCircleSquareIsAClosedSolidThroughTheSections)
{
  const std::vector<TopoDS_Wire> wires = {squareAt(0, 10), circleAt(20, 6), squareAt(40, 9)};
  std::vector<uint32_t>          sections;
  for (const TopoDS_Wire& wire : wires)
    sections.push_back(BRepExt::store(session(), wire));
  const auto          result = BRepExt::skinThrough(session(), sections, 0, 0, true, false);
  const TopoDS_Shape& made   = BRepExt::get(session(), result.resultId);
  EXPECT_TRUE(BRepExt::invalidity(made).empty());
  EXPECT_GT(volume(made), 0);
  TopoDS_Compound skin;
  BRep_Builder    builder;
  builder.MakeCompound(skin);
  for (TopExp_Explorer f(made, TopAbs_FACE); f.More(); f.Next())
    builder.Add(skin, f.Current());
  for (const TopoDS_Wire& wire : wires)
    for (TopExp_Explorer e(wire, TopAbs_EDGE); e.More(); e.Next())
    {
      BRepAdaptor_Curve curve(TopoDS::Edge(e.Current()));
      for (int k = 0; k <= 8; ++k)
      {
        const gp_Pnt                     p = curve.Value(curve.FirstParameter()
                                     + (curve.LastParameter() - curve.FirstParameter()) * k / 8);
        const BRepExtrema_DistShapeShape d(BRepBuilderAPI_MakeVertex(p).Vertex(), skin);
        ASSERT_TRUE(d.IsDone());
        EXPECT_LE(d.Value(), BRepExt::tolerance::sewing);
      }
    }
}

// A disc of radius r swept along a circular arc of radius rho through angle theta, square to the
// path all along, has the volume pi r^2 rho theta (Pappus: area times the path of its centroid). A
// twist turns the disc about its own axis and changes nothing; a scale s at the far end, linear
// along the path, makes the area pi r^2 (1 + (s - 1) t)^2 at fraction t, so the volume is pi r^2
// rho theta (1 + s + s^2) / 3.
TEST(BRepSkinExt, DiscSweptAlongAnArcFollowsPappus)
{
  const double      rho = 30, r = 4, theta = M_PI / 2;
  const TopoDS_Wire path =
    BRepBuilderAPI_MakeWire(
      BRepBuilderAPI_MakeEdge(
        new Geom_TrimmedCurve(new Geom_Circle(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), rho),
                              0,
                              theta))
        .Edge())
      .Wire();
  const TopoDS_Wire disc =
    BRepBuilderAPI_MakeWire(
      BRepBuilderAPI_MakeEdge(new Geom_Circle(gp_Ax2(gp_Pnt(rho, 0, 0), gp_Dir(0, 1, 0)), r))
        .Edge())
      .Wire();
  const uint32_t profile = BRepExt::store(session(), disc), along = BRepExt::store(session(), path);
  const double   base = M_PI * r * r * rho * theta;
  for (const double twist : {0.0, 0.7})
  {
    const auto result =
      BRepExt::sweepFramed(session(), profile, along, twist, 1, {}, 0, 1e-6, true);
    const TopoDS_Shape& made = BRepExt::get(session(), result.resultId);
    EXPECT_TRUE(BRepExt::invalidity(made).empty());
    EXPECT_NEAR(volume(made), base, 1e-5 * base) << "twist " << twist;
  }
  const double s      = 1.5;
  const auto   scaled = BRepExt::sweepFramed(session(), profile, along, 0, s, {}, 0, 1e-6, true);
  EXPECT_NEAR(volume(BRepExt::get(session(), scaled.resultId)),
              base * (1 + s + s * s) / 3,
              1e-4 * base);
}

// A Gordon skin of three quarter circles of radius 4 at heights 0, 5 and 10 and the three
// generators at 0, 45 and 90 degrees through their crossings is a patch of the cylinder x^2 + y^2
// = 16: the profiles interpolated along the height are exact (equal arcs at equally spaced
// heights), and the guides, being straight, add nothing the grid points do not already carry. Its
// area is r (pi / 2) h = 20 pi and every sampled point is at radius 4 from the axis. A
// profiles and generators are given out of order, which the skin puts in order. A
// network that is crossed wrongly, a guide taken for a profile, or the correction term (the
// tensor interpolation of the grid) left out, gives a surface off the cylinder at first order.
TEST(BRepSkinExt, GordonSkinOfArcsAndGeneratorsIsAPatchOfTheCylinder)
{
  const double          r = 4;
  std::vector<uint32_t> profiles, guides;
  for (const double z : {5.0, 10.0, 0.0})
  {
    Handle(Geom_Circle) circle = new Geom_Circle(gp_Ax2(gp_Pnt(0, 0, z), gp_Dir(0, 0, 1)), r);
    profiles.push_back(BRepExt::store(
      session(),
      BRepBuilderAPI_MakeWire(
        BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(circle->Circ(), 0, M_PI / 2, true).Value())
          .Edge())
        .Wire()));
  }
  for (const double angle : {M_PI / 2, 0.0, M_PI / 4})
  {
    const gp_Pnt a(r * std::cos(angle), r * std::sin(angle), 0), b(a.X(), a.Y(), 10);
    guides.push_back(
      BRepExt::store(session(),
                     BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(a, b).Edge()).Wire()));
  }
  const auto          result = BRepExt::skinGordon(session(), profiles, guides, false, 1e-7);
  const TopoDS_Shape& made   = BRepExt::get(session(), result.resultId);
  EXPECT_NEAR(area(made), 20 * M_PI, 1e-6 * 20 * M_PI);
  double worst = 0;
  for (TopExp_Explorer it(made, TopAbs_FACE); it.More(); it.Next())
  {
    BRepAdaptor_Surface surface(TopoDS::Face(it.Current()));
    for (int i = 0; i <= 8; ++i)
      for (int j = 0; j <= 8; ++j)
      {
        const gp_Pnt p =
          surface.Value(surface.FirstUParameter()
                          + (surface.LastUParameter() - surface.FirstUParameter()) * i / 8,
                        surface.FirstVParameter()
                          + (surface.LastVParameter() - surface.FirstVParameter()) * j / 8);
        worst = std::max(worst, std::abs(std::hypot(p.X(), p.Y()) - r));
      }
  }
  EXPECT_LT(worst, 1e-7);
}

// A guided skin through three straight profiles (the segment from (0,0,z) to (10,0,z) at z = 0, 5
// and 10) along a guide y = 2 sin(pi z / 5) at x = 5 that crosses each profile where y = 0 and
// bulges between them. The profiles alone make the plane y = 0; the guide's move off that plane is
// blended to the next crossing, so the result contains the profiles and the guide: points of each
// lie on it to the tolerance. A guide that is ignored, or a result that passes the crossings
// but not the bulge, leaves the guide's middle 2 away from the surface.
TEST(BRepSkinExt, GuidedSkinContainsItsProfilesAndABulgingGuide)
{
  std::vector<uint32_t> profiles, guides;
  for (const double z : {0.0, 5.0, 10.0})
    profiles.push_back(BRepExt::store(
      session(),
      BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(gp_Pnt(0, 0, z), gp_Pnt(10, 0, z)).Edge())
        .Wire()));
  NCollection_Array1<gp_Pnt> through(1, 21);
  for (int i = 0; i <= 20; ++i)
  {
    const double z = 0.5 * i;
    through(i + 1) = gp_Pnt(5, 2 * std::sin(M_PI * z / 5), z);
  }
  Handle(NCollection_HArray1<gp_Pnt>) points = new NCollection_HArray1<gp_Pnt>(through);
  GeomAPI_Interpolate                 fit(points, false, 1e-9);
  fit.Perform();
  ASSERT_TRUE(fit.IsDone());
  const TopoDS_Edge guide = BRepBuilderAPI_MakeEdge(fit.Curve()).Edge();
  guides.push_back(BRepExt::store(session(), BRepBuilderAPI_MakeWire(guide).Wire()));
  const auto          result = BRepExt::skinGuided(session(), profiles, guides, 1e-6);
  const TopoDS_Shape& made   = BRepExt::get(session(), result.resultId);
  double              worst  = 0;
  for (int i = 0; i <= 20; ++i)
  {
    BRepExtrema_DistShapeShape distance(
      BRepBuilderAPI_MakeVertex(gp_Pnt(5, 2 * std::sin(M_PI * 0.5 * i / 5), 0.5 * i)).Vertex(),
      made);
    ASSERT_TRUE(distance.IsDone());
    worst = std::max(worst, distance.Value());
  }
  for (const double z : {0.0, 5.0, 10.0})
    for (const double x : {0.0, 2.5, 5.0, 10.0})
    {
      BRepExtrema_DistShapeShape distance(BRepBuilderAPI_MakeVertex(gp_Pnt(x, 0, z)).Vertex(),
                                          made);
      ASSERT_TRUE(distance.IsDone());
      worst = std::max(worst, distance.Value());
    }
  EXPECT_LT(worst, 1e-5);
}
