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
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRep_Tool.hxx>
#include <BRepOffsetAPI_MakeFilling.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <gp_Trsf.hxx>

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

namespace
{
//! A seeded draw in [0, 1): the minimal standard generator, the same on every platform.
struct Draw
{
  uint64_t state;

  double operator()()
  {
    state = state * 48271 % 2147483647;
    return double(state) / 2147483647;
  }
};

//! The Bezier curve of these poles, as a spline.
Handle(Geom_BSplineCurve) bezier(const std::vector<gp_Pnt>& poles)
{
  NCollection_Array1<gp_Pnt> array(1, int(poles.size()));
  for (size_t i = 0; i < poles.size(); ++i)
    array(int(i) + 1) = poles[i];
  NCollection_Array1<double> knots(1, 2);
  NCollection_Array1<int>    mults(1, 2);
  knots(1) = 0;
  knots(2) = 1;
  mults(1) = mults(2) = int(poles.size());
  return new Geom_BSplineCurve(array, knots, mults, int(poles.size()) - 1);
}

TopoDS_Edge bezierEdge(const std::vector<gp_Pnt>& poles)
{
  return BRepBuilderAPI_MakeEdge(bezier(poles)).Edge();
}
} // namespace

// A guided skin through three to five profiles between two planar supports at any tilt, G1 or G2
// at each end, in any pose and at any size: guides of degree 5 leave each support in its plane
// with their first three poles in it (no bend out of it), interior profiles bulge between the
// guides. For newly supported nets every profile and guide lies on the skin, and at each support
// the skin's normal is the
// plane's, and at G2 its bend across is the plane's, nothing. The bug: through more than two
// profiles G2 was a plate filling that missed the bend limit, or the boundary, on nearly all of
// these.
TEST(BRepSkinExt, GuidedSkinThroughSeveralProfilesMeetsPlanarSupports)
{
  for (uint64_t seed = 1; seed <= 12; ++seed)
  {
    SCOPED_TRACE(seed);
    Draw         draw{seed * 7919 + 13};
    const double size = std::pow(10.0, 2 * draw() - 1), width = size * (2 + 2 * draw()),
                 length = size * (4 + 3 * draw()), drop = size * (draw() - 0.5) * 2;
    const int    inner     = 1 + int(draw() * 3);
    const int    orders[2] = {1 + int(draw() * 2), seed % 2 == 0 ? 2 : 1 + int(draw() * 2)};
    const double tilt[2]   = {(draw() - 0.5) * 1.2, (draw() - 0.5) * 1.2};
    const gp_Vec way[2]    = {gp_Vec(0, std::cos(tilt[0]), std::sin(tilt[0])),
                              gp_Vec(0, std::cos(tilt[1]), std::sin(tilt[1]))};
    gp_Trsf      pose;
    pose.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(draw() - 0.5, draw() - 0.5, draw() + 0.5)),
                     2 * M_PI * draw());
    pose.SetTranslationPart(gp_Vec(size * draw(), size * draw(), size * draw()));
    const gp_Pnt rim[2] = {gp_Pnt(0, 0, 0), gp_Pnt(0, length, drop)};
    // The supports: strips of width `size` beyond each end profile, in the planes the guides leave.
    const auto strip = [&](int end) {
      const gp_Vec            out  = way[end] * (end == 0 ? -size : size);
      const gp_Pnt            c[4] = {rim[end],
                                      rim[end].Translated(gp_Vec(width, 0, 0)),
                                      rim[end].Translated(gp_Vec(width, 0, 0) + out),
                                      rim[end].Translated(out)};
      BRepBuilderAPI_MakeWire wire;
      for (int i = 0; i < 4; ++i)
        wire.Add(BRepBuilderAPI_MakeEdge(c[i], c[(i + 1) % 4]).Edge());
      return TopoDS::Face(
        BRepBuilderAPI_MakeFace(wire.Wire(), true).Face().Moved(TopLoc_Location(pose)));
    };
    const TopoDS_Face                supports[2] = {strip(0), strip(1)};
    std::vector<std::vector<gp_Pnt>> guidePoles;
    for (const double x : {0.0, width})
    {
      const gp_Vec along(x, 0, 0);
      const double s0 = length * (0.1 + 0.1 * draw()), s1 = length * (0.1 + 0.1 * draw());
      const double lift[2] = {size * (draw() - 0.5), size * (draw() - 0.5)};
      const gp_Pnt a = rim[0].Translated(along), b = rim[1].Translated(along);
      guidePoles.push_back({a,
                            a.Translated(way[0] * s0),
                            a.Translated(way[0] * (2 * s0 + lift[0])),
                            b.Translated(way[1] * -(2 * s1 + lift[1])),
                            b.Translated(way[1] * -s1),
                            b});
    }
    std::vector<TopoDS_Edge> guides, inners;
    for (const auto& poles : guidePoles)
      guides.push_back(TopoDS::Edge(bezierEdge(poles).Moved(TopLoc_Location(pose))));
    for (int k = 1; k <= inner; ++k)
    {
      const double t = double(k) / (inner + 1) + (draw() - 0.5) * 0.1;
      const gp_Pnt p = bezier(guidePoles[0])->Value(t), q = bezier(guidePoles[1])->Value(t);
      const gp_Pnt middle((p.XYZ() + q.XYZ()) / 2
                          + gp_XYZ(0, size * (draw() - 0.5) * 0.6, size * (draw() - 0.5) * 1.2));
      inners.push_back(TopoDS::Edge(bezierEdge({p, middle, q}).Moved(TopLoc_Location(pose))));
    }
    const auto rimOf = [&](int end) {
      for (TopExp_Explorer it(supports[end], TopAbs_EDGE); it.More(); it.Next())
      {
        BRepAdaptor_Curve c(TopoDS::Edge(it.Current()));
        const gp_Pnt      m = c.Value((c.FirstParameter() + c.LastParameter()) / 2);
        if (m.Distance(rim[end].Translated(gp_Vec(width / 2, 0, 0)).Transformed(pose))
            < 1e-9 * size)
          return TopoDS::Edge(it.Current());
      }
      return TopoDS_Edge();
    };
    std::vector<uint32_t> edges{BRepExt::store(session(), rimOf(0))};
    for (const auto& e : inners)
      edges.push_back(BRepExt::store(session(), e));
    edges.push_back(BRepExt::store(session(), rimOf(1)));
    for (const auto& g : guides)
      edges.push_back(BRepExt::store(session(), g));
    const uint32_t     faceIds[2] = {BRepExt::store(session(), supports[0]),
                                     BRepExt::store(session(), supports[1])};
    // Keep accepted plate results to the original acceptance contract. The
    // stronger exact-net/jet checks below apply to newly supported nets.
    const double tol = 1e-6*size;
    BRepOffsetAPI_MakeFilling reference(3,50,5,false,tol/10,tol/10,1e-5,1e-4,8,128);
    const int lastProfile = int(edges.size())-3;
    for (int i=0; i<int(edges.size()); ++i)
    {
      const auto edge = TopoDS::Edge(BRepExt::get(session(),edges[i]));
      if (i==0 || i==lastProfile)
        reference.Add(edge, supports[i==0 ? 0 : 1], static_cast<GeomAbs_Shape>(orders[i==0 ? 0 : 1]));
      else
        reference.Add(edge, GeomAbs_C0, i>=lastProfile);
    }
    reference.Build();
    const bool accepted = reference.IsDone() && reference.G0Error()<=tol
      && BRepCheck_Analyzer(reference.Shape()).IsValid() && reference.G1Error()<=1e-4
      && ((orders[0]<2 && orders[1]<2) || reference.G2Error()<=1e-3);
    const TopoDS_Shape made       = BRepExt::get(session(),
                                           BRepExt::skinGuides(session(),
                                                               edges,
                                                               faceIds[0],
                                                               faceIds[1],
                                                               orders[0],
                                                               orders[1],
                                                               1e-6 * size)
                                             .resultId);
    ASSERT_TRUE(BRepCheck_Analyzer(made).IsValid());
    if (accepted)
    {
      const auto expected = occ::handle<Geom_BSplineSurface>::DownCast(
        BRep_Tool::Surface(TopoDS::Face(reference.Shape())));
      const auto actual = occ::handle<Geom_BSplineSurface>::DownCast(
        BRep_Tool::Surface(TopoDS::Face(made)));
      ASSERT_FALSE(expected.IsNull()); ASSERT_FALSE(actual.IsNull());
      ASSERT_EQ(actual->UDegree(),expected->UDegree());
      ASSERT_EQ(actual->VDegree(),expected->VDegree());
      ASSERT_EQ(actual->NbUPoles(),expected->NbUPoles());
      ASSERT_EQ(actual->NbVPoles(),expected->NbVPoles());
      for (int i=1; i<=actual->NbUPoles(); ++i)
        for (int j=1; j<=actual->NbVPoles(); ++j)
          EXPECT_LE(actual->Pole(i,j).Distance(expected->Pole(i,j)),1e-9*size);
      const auto oldBend = BRepExt::faceFairness(session(), BRepExt::store(session(),reference.Shape()),128);
      const auto newBend = BRepExt::faceFairness(session(), BRepExt::store(session(),made),128);
      EXPECT_NEAR(newBend[1]/oldBend[1],1.0,1e-3);
      continue;
    }
    // The net on the skin.
    double off = 0;
    for (const auto* list : {&inners, &guides})
      for (const auto& e : *list)
      {
        BRepAdaptor_Curve c(e);
        for (int i = 0; i <= 16; ++i)
        {
          BRepExtrema_DistShapeShape d(
            BRepBuilderAPI_MakeVertex(
              c.Value(c.FirstParameter() + (c.LastParameter() - c.FirstParameter()) * i / 16))
              .Vertex(),
            made);
          off = std::max(off, d.Value());
        }
      }
    EXPECT_LT(off, 1e-7 * size);
    // At each support: the normal of the plane, and at G2 no bend across.
    TopoDS_Face skin;
    for (TopExp_Explorer it(made, TopAbs_FACE); it.More(); it.Next())
      skin = TopoDS::Face(it.Current());
    const Handle(Geom_Surface) surface = BRep_Tool::Surface(skin);
    for (int end = 0; end < 2; ++end)
    {
      const gp_Dir      n    = gp_Dir(way[end].Crossed(gp_Vec(1, 0, 0))).Transformed(pose);
      double            turn = 0, bend = 0;
      BRepAdaptor_Curve c(rimOf(end));
      for (int i = 0; i <= 16; ++i)
      {
        const gp_Pnt p =
          c.Value(c.FirstParameter() + (c.LastParameter() - c.FirstParameter()) * i / 16);
        ShapeAnalysis_Surface on(surface);
        const gp_Pnt2d        uv = on.ValueOfUV(p, 1e-9 * size);
        gp_Pnt                q;
        gp_Vec                su, sv, suu, svv, suv;
        surface->D2(uv.X(), uv.Y(), q, su, sv, suu, svv, suv);
        const gp_Vec normal = su.Crossed(sv).Normalized();
        turn                = std::max(turn, normal.Crossed(gp_Vec(n)).Magnitude());
        // The bend across: the normal curvature in the direction square to the rim.
        const gp_Vec along =
          c.DN(c.FirstParameter() + (c.LastParameter() - c.FirstParameter()) * i / 16, 1);
        const gp_Vec across = normal.Crossed(along).Normalized();
        const double E = su.Dot(su), F = su.Dot(sv), G = sv.Dot(sv), det = E * G - F * F;
        const double a = (G * across.Dot(su) - F * across.Dot(sv)) / det,
                     b = (E * across.Dot(sv) - F * across.Dot(su)) / det;
        bend =
          std::max(bend, std::abs(normal.Dot(suu * (a * a) + suv * (2 * a * b) + svv * (b * b))));
      }
      EXPECT_LT(turn, 1e-6);
      if (orders[end] == 2)
        EXPECT_LT(bend * size, 1e-5);
    }
  }
}

// Exact interpolation alone does not minimize bending energy. Preserve a plate
// that already satisfies the requested constraints, including its spline space.
TEST(BRepSkinExt, GuidedSkinRetainsAnAcceptedPlate)
{
  for (int np : {3, 4, 5})
  {
    const auto support = [](double y0, double y1) {
      const gp_Pnt points[] = {gp_Pnt(0,y0,0), gp_Pnt(4,y0,0), gp_Pnt(4,y1,0), gp_Pnt(0,y1,0)};
      BRepBuilderAPI_MakeWire w;
      for (int i=0; i<4; ++i) w.Add(BRepBuilderAPI_MakeEdge(points[i],points[(i+1)%4]));
      return BRepBuilderAPI_MakeFace(w.Wire()).Face();
    };
    const TopoDS_Face first = support(0,-2), last = support(8,10);
    std::vector<TopoDS_Edge> edgeShapes;
    for (int i=0; i<np; ++i)
      edgeShapes.push_back(BRepBuilderAPI_MakeEdge(gp_Pnt(0,8.0*i/(np-1),0),
                                                  gp_Pnt(4,8.0*i/(np-1),0)));
    edgeShapes.push_back(BRepBuilderAPI_MakeEdge(gp_Pnt(0,0,0),gp_Pnt(0,8,0)));
    edgeShapes.push_back(BRepBuilderAPI_MakeEdge(gp_Pnt(4,0,0),gp_Pnt(4,8,0)));
    const auto edgeAt = [](const TopoDS_Face& f, double y) {
      for (TopExp_Explorer it(f,TopAbs_EDGE); it.More(); it.Next())
      {
        const auto e = TopoDS::Edge(it.Current());
        BRepAdaptor_Curve c(e);
        if (std::abs(c.Value((c.FirstParameter()+c.LastParameter())/2).Y()-y)<1e-12) return e;
      }
      return TopoDS_Edge();
    };
    edgeShapes.front() = edgeAt(first,0);
    edgeShapes[np-1] = edgeAt(last,8);
    BRepOffsetAPI_MakeFilling plate(3,50,5,false,1e-7,1e-7,1e-5,1e-4,8,128);
    std::vector<uint32_t> ids;
    for (int i=0; i<int(edgeShapes.size()); ++i)
    {
      ids.push_back(BRepExt::store(session(), edgeShapes[i]));
      if (i==0 || i==np-1)
        plate.Add(edgeShapes[i], i==0 ? first : last, GeomAbs_C1);
      else
        plate.Add(edgeShapes[i], GeomAbs_C0, i>=np-1);
    }
    plate.Build();
    ASSERT_TRUE(plate.IsDone());
    ASSERT_LE(plate.G0Error(), 1e-6);
    ASSERT_LE(plate.G1Error(), 1e-4);
    const auto result = BRepExt::skinGuides(session(), ids, BRepExt::store(session(),first),
                                           BRepExt::store(session(),last), 1,1,1e-6);
    const auto expected = occ::handle<Geom_BSplineSurface>::DownCast(
      BRep_Tool::Surface(TopoDS::Face(plate.Shape())));
    const auto actual = occ::handle<Geom_BSplineSurface>::DownCast(
      BRep_Tool::Surface(TopoDS::Face(BRepExt::get(session(),result.resultId))));
    ASSERT_FALSE(expected.IsNull()); ASSERT_FALSE(actual.IsNull());
    ASSERT_EQ(actual->UDegree(),expected->UDegree());
    ASSERT_EQ(actual->VDegree(),expected->VDegree());
    ASSERT_EQ(actual->NbUPoles(),expected->NbUPoles());
    ASSERT_EQ(actual->NbVPoles(),expected->NbVPoles());
    for (int i=1; i<=actual->NbUPoles(); ++i)
      for (int j=1; j<=actual->NbVPoles(); ++j)
        EXPECT_LE(actual->Pole(i,j).Distance(expected->Pole(i,j)), 1e-9);

  }
}

// A nonplanar closed end cannot be capped by the planar-cap construction.
TEST(BRepSkinExt, NonplanarEndsCannotPublishAnInvalidSolid)
{
  for (const double scale : {0.1, 1.0, 10.0})
  {
    std::vector<uint32_t> sections;
    for (const double h : {0.0, 5.0, 10.0})
    {
      const gp_Pnt p[] = {gp_Pnt(0, 0, h * scale), gp_Pnt(10 * scale, 0, h * scale),
                          gp_Pnt(10 * scale, 10 * scale, (h + 1) * scale),
                          gp_Pnt(0, 10 * scale, h * scale)};
      BRepBuilderAPI_MakeWire wire;
      for (int i = 0; i < 4; ++i)
        wire.Add(BRepBuilderAPI_MakeEdge(p[i], p[(i + 1) % 4]).Edge());
      sections.push_back(BRepExt::store(session(), wire.Wire()));
    }
    try
    {
      const auto result = BRepExt::skinThrough(session(), sections, 0, 0, true, false);
      EXPECT_TRUE(BRepCheck_Analyzer(BRepExt::get(session(), result.resultId)).IsValid());
    }
    catch (const BRepExt::Failure& failure)
    {
      EXPECT_EQ(failure.reason, BRepExt::no_solution);
    }
    const auto sheet = BRepExt::skinThrough(session(), sections, 0, 0, false, false);
    EXPECT_TRUE(BRepCheck_Analyzer(BRepExt::get(session(), sheet.resultId)).IsValid());
  }
}
