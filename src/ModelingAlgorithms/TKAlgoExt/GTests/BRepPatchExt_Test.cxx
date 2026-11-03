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

// Patch against what its contract states, measured without its help: a planar boundary gives the
// disc, and a boundary asked to be smooth with a support is smooth with it, by normals found from
// the two surfaces alone.

#include "BRepExtTest_Support.hxx"

#include <BRepToolsExt_Operations.hxx>
#include <BRepPatchExt_Inside.hxx>
#include <BRepPatchExt_Side.hxx>
#include <GeomAlgoExt_Failure.hxx>

#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRep_Tool.hxx>
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomLProp_SLProps.hxx>
#include <GeomConvert.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Circle.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Plane.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>

#include <gtest/gtest.h>

#include <regex>

using namespace BRepExtTest;

namespace
{

TopoDS_Face onlyFace(const TopoDS_Shape& shape)
{
  TopoDS_Face face;
  for (TopExp_Explorer it(shape, TopAbs_FACE); it.More(); it.Next())
    face = TopoDS::Face(it.Current());
  return face;
}

} // namespace

// The face inside a circle is the disc, of area pi r^2, whether Patch takes its shortcut for planar
// boundaries (method 0) or is made to build a surface (method 1): a fitted surface that bulges
// has more area, and one that leaves the boundary has a different area.
TEST(BRepPatchExt, FillOfACircleHasTheAreaOfTheDisc)
{
  const double      r = 12;
  const TopoDS_Edge circle =
    BRepBuilderAPI_MakeEdge(new Geom_Circle(gp_Ax2(gp_Pnt(3, -2, 5), gp_Dir(0, 0, 1)), r)).Edge();
  const uint32_t edge = BRepExt::store(session(), circle);
  for (const int method : {0, 1})
  {
    const auto result = BRepExt::patch(session(), {edge}, {0}, {0}, {}, {}, 1e-6, 1, 0.1, method);
    const TopoDS_Face face = onlyFace(BRepExt::get(session(), result.evolution.resultId));
    ASSERT_FALSE(face.IsNull()) << "method " << method;
    EXPECT_NEAR(area(face), M_PI * r * r, 1e-4 * M_PI * r * r) << "method " << method;
  }
}

// A sphere with the cap above latitude 40 degrees missing, and the hole filled asking for G1 with
// the sphere along its rim: the angle between the normals of the sphere and of the fill, found by
// projecting points of the rim onto both surfaces, stays within the g1Degrees asked.
TEST(BRepPatchExt, FillOfACapHoleIsTangentToTheSphere)
{
  const double       R = 10, latitude = 40 * M_PI / 180, g1Degrees = 1.0;
  const TopoDS_Shape sphere =
    BRepPrimAPI_MakeSphere(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), R, -M_PI / 2, latitude)
      .Shape();
  const TopoDS_Face support = onlyFace(sphere);
  const auto        rim     = edgesWhere(sphere, [&](const double* b) {
    return near(b[2], R * std::sin(latitude)) && near(b[5], R * std::sin(latitude));
  });
  ASSERT_FALSE(rim.empty());
  std::vector<uint32_t> edges, supports;
  for (const TopoDS_Edge& e : rim)
  {
    edges.push_back(BRepExt::store(session(), e));
    supports.push_back(BRepExt::store(session(), support));
  }
  const auto        result = BRepExt::patch(session(),
                                     edges,
                                     supports,
                                     std::vector<int>(edges.size(), 1),
                                            {},
                                            {},
                                     1e-4,
                                     g1Degrees,
                                     0.1,
                                     0);
  const TopoDS_Face fill   = onlyFace(BRepExt::get(session(), result.evolution.resultId));
  ASSERT_FALSE(fill.IsNull());
  const Handle(Geom_Surface) sphereSurface = BRep_Tool::Surface(support),
                             fillSurface   = BRep_Tool::Surface(fill);
  double worst                             = 0;
  for (const TopoDS_Edge& e : rim)
  {
    BRepAdaptor_Curve curve(e);
    for (int k = 1; k < 32; ++k)
    {
      const gp_Pnt               p = curve.Value(curve.FirstParameter()
                                   + (curve.LastParameter() - curve.FirstParameter()) * k / 32);
      gp_Dir                     normal[2];
      const Handle(Geom_Surface) surfaces[2] = {sphereSurface, fillSurface};
      for (int s = 0; s < 2; ++s)
      {
        GeomAPI_ProjectPointOnSurf project(p, surfaces[s]);
        ASSERT_GT(project.NbPoints(), 0);
        ASSERT_LT(project.LowerDistance(), 1e-4);
        double u, v;
        project.LowerDistanceParameters(u, v);
        GeomLProp_SLProps props(surfaces[s], u, v, 1, 1e-9);
        ASSERT_TRUE(props.IsNormalDefined());
        normal[s] = props.Normal();
      }
      worst =
        std::max(worst, std::acos(std::min(1., std::abs(normal[0].Dot(normal[1])))) * 180 / M_PI);
    }
  }
  EXPECT_LE(worst, g1Degrees);
}

// A boundary that does not close: three quarter arcs of a circle of radius r, given in any order
// and any direction, with the fourth missing. They are chained into one run and the gap between its
// ends is closed by a chord, so the face is the circular segment cut off by the chord of a 270
// degree arc: area r^2 (theta - sin theta) / 2 with theta = 3 pi / 2. A wrong chain (an arc left
// out, a gap closed to another end) or a gap closed off its chord changes the area at first order.
TEST(BRepPatchExt, OpenBoundaryOfThreeArcsIsClosedByTheChordOfTheGap)
{
  const double             r = 9;
  const gp_Ax2             axis(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
  std::vector<TopoDS_Edge> arcs;
  for (const int k : {0, 1, 2})
  {
    Handle(Geom_Circle) circle = new Geom_Circle(axis, r);
    arcs.push_back(BRepBuilderAPI_MakeEdge(circle, k * M_PI / 2, (k + 1) * M_PI / 2).Edge());
  }
  // Out of order, and the middle one against the others.
  arcs[1].Reverse();
  const std::vector<TopoDS_Edge> given = {arcs[2], arcs[1], arcs[0]};
  std::vector<uint32_t>          edges;
  for (const TopoDS_Edge& edge : given)
    edges.push_back(BRepExt::store(session(), edge));
  BRepExt::PatchControls controls;
  controls.open            = true;
  const auto        result = BRepExt::patch(session(),
                                     edges,
                                     std::vector<uint32_t>(3, 0),
                                     std::vector<int>(3, 0),
                                            {},
                                            {},
                                     1e-6,
                                     1,
                                     0.1,
                                     0,
                                     -1,
                                     controls);
  const TopoDS_Face face   = onlyFace(BRepExt::get(session(), result.evolution.resultId));
  ASSERT_FALSE(face.IsNull());
  const double theta = 3 * M_PI / 2;
  EXPECT_NEAR(area(face), r * r * (theta - std::sin(theta)) / 2, 1e-3);
}

namespace
{

//! The six plates of a dish round a regular hexagon of circumradius 10: each rises from its edge
//! outwards at 30 degrees. Returns the edges and the plates; `normals` gets each plate's normal.
void dish(std::vector<uint32_t>& edges,
          std::vector<uint32_t>& supports,
          std::vector<gp_Dir>&   normals)
{
  const double rho = 10, rise = 30 * M_PI / 180, reach = 6;
  gp_Pnt       vertex[6];
  for (int k = 0; k < 6; ++k)
    vertex[k] = gp_Pnt(rho * std::cos(k * M_PI / 3), rho * std::sin(k * M_PI / 3), 0);
  for (int k = 0; k < 6; ++k)
  {
    const gp_Pnt &    a = vertex[k], &b = vertex[(k + 1) % 6];
    const gp_Vec      along(a, b);
    const gp_Vec      outward(gp_Pnt(0, 0, 0), gp_Pnt((a.X() + b.X()) / 2, (a.Y() + b.Y()) / 2, 0));
    const gp_Vec      up = outward.Normalized() * std::cos(rise) + gp_Vec(0, 0, std::sin(rise));
    const gp_Pnt      c = b.Translated(up * reach), d = a.Translated(up * reach);
    const TopoDS_Face plate =
      BRepBuilderAPI_MakeFace(BRepBuilderAPI_MakePolygon(a, b, c, d, true).Wire()).Face();
    const auto edge = edgesWhere(plate, [&](const double* l) {
      return near(l[0], std::min(a.X(), b.X())) && near(l[3], std::max(a.X(), b.X()))
             && near(l[1], std::min(a.Y(), b.Y())) && near(l[4], std::max(a.Y(), b.Y()))
             && near(l[2], 0) && near(l[5], 0);
    });
    edges.push_back(BRepExt::store(session(), edge.at(0)));
    supports.push_back(BRepExt::store(session(), plate));
    normals.push_back(gp_Dir(along.Crossed(up)));
  }
}

} // namespace

// Two neighbouring plates of the dish are 29 degrees apart at the vertex they share, so no surface
// can be tangent to both there (it has one tangent plane at the vertex). Asking G1 along all six is
// refused with that reason; filling the loop on a wrong guess, or silently as G0, would break the
// contract the caller asked for.
TEST(BRepPatchExt, G1FillBetweenPlatesThatDisagreeAtACornerIsRefused)
{
  std::vector<uint32_t> edges, supports;
  std::vector<gp_Dir>   normals;
  dish(edges, supports, normals);
  try
  {
    BRepExt::patch(session(), edges, supports, std::vector<int>(6, 1), {}, {}, 1e-4, 1.0, 0.1, 0);
    FAIL() << "a G1 fill was returned";
  }
  catch (const BRepExt::Failure& failure)
  {
    EXPECT_NE(std::string(failure.what()).find("no surface is tangent to both"), std::string::npos)
      << failure.what();
  }
}

// The same loop asked only for G0 is filled: it passes through every hexagon vertex and, the six
// plates being equal under a rotation of 60 degrees about the axis, has its middle on that axis.
// A fill that loses a side, or is pulled to one plate, misses a vertex or the axis.
TEST(BRepPatchExt, G0FillOfASixSidedLoopHoldsEveryVertexAndTheSymmetry)
{
  std::vector<uint32_t> edges, supports;
  std::vector<gp_Dir>   normals;
  dish(edges, supports, normals);
  const auto result =
    BRepExt::patch(session(), edges, supports, std::vector<int>(6, 0), {}, {}, 1e-4, 1.0, 0.1, 0);
  const TopoDS_Face fill = onlyFace(BRepExt::get(session(), result.evolution.resultId));
  ASSERT_FALSE(fill.IsNull());
  const Handle(Geom_Surface) surface = BRep_Tool::Surface(fill);
  for (int k = 0; k < 6; ++k)
  {
    const gp_Pnt               corner(10 * std::cos(k * M_PI / 3), 10 * std::sin(k * M_PI / 3), 0);
    GeomAPI_ProjectPointOnSurf project(corner, surface);
    ASSERT_GT(project.NbPoints(), 0);
    EXPECT_LT(project.LowerDistance(), 1e-4) << "vertex " << k;
  }
  GeomAPI_ProjectPointOnSurf middle(gp_Pnt(0, 0, 3), surface);
  ASSERT_GT(middle.NbPoints(), 0);
  EXPECT_NEAR(middle.NearestPoint().X(), 0, 1e-3);
  EXPECT_NEAR(middle.NearestPoint().Y(), 0, 1e-3);
}

// A six-sided hole in a sphere, its rim six arcs of the circle of latitude 40 degrees of a sphere
// of radius 10 given in scrambled order and direction, filled for G1 with the sphere: the angle
// between the sphere's and the fill's normals along every arc stays within g1Degrees. Chaining
// six edges into the loop, and carrying the support of each, is what is tested beyond the one-edge
// rim above.
TEST(BRepPatchExt, SixArcRimOnASphereIsFilledTangentToItThroughEveryArc)
{
  const double       R = 10, latitude = 40 * M_PI / 180, g1Degrees = 1.0;
  const TopoDS_Shape sphere =
    BRepPrimAPI_MakeSphere(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), R, -M_PI / 2, latitude)
      .Shape();
  const TopoDS_Face   support = onlyFace(sphere);
  const double        ringR = R * std::cos(latitude), ringZ = R * std::sin(latitude);
  Handle(Geom_Circle) ring = new Geom_Circle(gp_Ax2(gp_Pnt(0, 0, ringZ), gp_Dir(0, 0, 1)), ringR);
  std::vector<TopoDS_Edge> arcs;
  for (const int k : {3, 0, 5, 1, 4, 2})
  {
    TopoDS_Edge arc = BRepBuilderAPI_MakeEdge(ring, k * M_PI / 3, (k + 1) * M_PI / 3).Edge();
    arcs.push_back(k % 2 ? TopoDS::Edge(arc.Reversed()) : arc);
  }
  std::vector<uint32_t> edges, supports;
  for (const TopoDS_Edge& arc : arcs)
  {
    edges.push_back(BRepExt::store(session(), arc));
    supports.push_back(BRepExt::store(session(), support));
  }
  const auto        result = BRepExt::patch(session(),
                                     edges,
                                     supports,
                                     std::vector<int>(6, 1),
                                            {},
                                            {},
                                     1e-4,
                                     g1Degrees,
                                     0.1,
                                     0);
  const TopoDS_Face fill   = onlyFace(BRepExt::get(session(), result.evolution.resultId));
  ASSERT_FALSE(fill.IsNull());
  const Handle(Geom_Surface) sphereSurface = BRep_Tool::Surface(support),
                             fillSurface   = BRep_Tool::Surface(fill);
  double worst                             = 0;
  for (const TopoDS_Edge& arc : arcs)
  {
    BRepAdaptor_Curve curve(arc);
    for (int i = 1; i < 16; ++i)
    {
      const gp_Pnt               p = curve.Value(curve.FirstParameter()
                                   + (curve.LastParameter() - curve.FirstParameter()) * i / 16);
      gp_Dir                     normal[2];
      const Handle(Geom_Surface) surfaces[2] = {sphereSurface, fillSurface};
      for (int s = 0; s < 2; ++s)
      {
        GeomAPI_ProjectPointOnSurf project(p, surfaces[s]);
        ASSERT_GT(project.NbPoints(), 0);
        ASSERT_LT(project.LowerDistance(), 1e-4);
        double u, v;
        project.LowerDistanceParameters(u, v);
        GeomLProp_SLProps props(surfaces[s], u, v, 1, 1e-9);
        ASSERT_TRUE(props.IsNormalDefined());
        normal[s] = props.Normal();
      }
      worst =
        std::max(worst, std::acos(std::min(1., std::abs(normal[0].Dot(normal[1])))) * 180 / M_PI);
    }
  }
  EXPECT_LE(worst, g1Degrees);
}

namespace
{

//! The end of a half tube of radius 20 below z = 0, on the side y <= 0, and the end of a plate 26
//! wide at height 50 that leaves the edge y = 0 along +y: two edges, a gap of about 50 between
//! them on either side. A fill from the one to the other is an elbow, turning from upright to
//! level. The tube is a cylinder, or with `sheet` the cubic spline through 25 by 5 of its points,
//! whose rim has more poles at degree 5 than a fitted quad may have.
struct Elbow
{
  TopoDS_Face           tube, plate;
  std::vector<uint32_t> edges, supports;

  explicit Elbow(bool sheet = false)
  {
    const Handle(Geom_CylindricalSurface) cylinder =
      new Geom_CylindricalSurface(gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, -1), gp_Dir(1, 0, 0)), 20);
    tube = BRepBuilderAPI_MakeFace(cylinder, 0, M_PI, 0, 20, 1e-7).Face();
    if (sheet)
    {
      NCollection_Array2<gp_Pnt> points(1, 25, 1, 5);
      for (int i = 1; i <= 25; ++i)
        for (int j = 1; j <= 5; ++j)
          points(i, j) = cylinder->Value(M_PI * (i - 1) / 24, 5.0 * (j - 1));
      GeomAPI_PointsToBSplineSurface through;
      through.Interpolate(points, Approx_IsoParametric);
      tube = BRepBuilderAPI_MakeFace(through.Surface(), 1e-7).Face();
    }
    plate = BRepBuilderAPI_MakeFace(BRepBuilderAPI_MakePolygon(gp_Pnt(-13, 0, 50),
                                                               gp_Pnt(13, 0, 50),
                                                               gp_Pnt(13, 22, 50),
                                                               gp_Pnt(-13, 22, 50),
                                                               true)
                                      .Wire())
              .Face();
    for (const TopoDS_Face& face : {tube, plate})
    {
      const double height = face.IsSame(tube) ? 0 : 50;
      for (const TopoDS_Edge& edge : edgesWhere(face, [&](const double* b) {
             return near(b[2], height) && near(b[5], height) && b[4] < 1e-6;
           }))
      {
        edges.push_back(BRepExt::store(session(), edge));
        supports.push_back(BRepExt::store(session(), face));
      }
    }
  }

  TopoDS_Face fill(int flow) const
  {
    BRepExt::PatchControls controls;
    controls.quad = 1;
    controls.flow = flow;
    controls.open = true;
    const auto result =
      BRepExt::patch(session(), edges, supports, {2, 2}, {}, {}, 1e-3, 0.1, 0.05, 0, -1, controls);
    return onlyFace(BRepExt::get(session(), result.evolution.resultId));
  }
};

//! The principal curvatures of `surface` at the point nearest `p`, least first, signed for its
//! normal turned to the side of `normal`, which then becomes that normal.
std::pair<double, double> bends(const Handle(Geom_Surface)& surface,
                                const gp_Pnt&               p,
                                gp_Dir&                     normal)
{
  GeomAPI_ProjectPointOnSurf project(p, surface);
  double                     u, v;
  project.LowerDistanceParameters(u, v);
  GeomLProp_SLProps props(surface, u, v, 2, 1e-9);
  EXPECT_TRUE(props.IsCurvatureDefined());
  const gp_Dir own   = props.Normal();
  const double sign  = own.Dot(normal) < 0 ? -1 : 1;
  const double least = sign * props.MinCurvature(), most = sign * props.MaxCurvature();
  normal = sign < 0 ? own.Reversed() : own;
  return {std::min(least, most), std::max(least, most)};
}

} // namespace

// The elbow asked G2 with both faces and Quad sided, the open boundary bridged across its gaps
// untrimmed: for Natural and Adjacent flow (with no side next to either rim, Adjacent leaves both
// rims square), and for the cylinder and the spline through its points, the fill keeps each face's
// tangent plane within the 0.1 degrees asked and its two principal curvatures within the 5 % asked
// (the plate's, which are none, within the flat share of a thousandth of that over 50, less than
// the loop's size), all found by projecting points of the rims onto both surfaces; and it is its
// own mirror image in x = 0, as the elbow is. A fill that missed a face, leant to one side or was
// refused fails.
TEST(BRepPatchExt, OpenElbowBetweenATubeAndAPlateIsFilledG2AcrossItsGaps)
{
  for (const bool sheet : {false, true})
    for (const int flow : {0, 3})
    {
      const Elbow elbow(sheet);
      ASSERT_EQ(elbow.edges.size(), 2u);
      const TopoDS_Face fill = elbow.fill(flow);
      ASSERT_FALSE(fill.IsNull());
      const Handle(Geom_Surface) surface = BRep_Tool::Surface(fill);
      for (const TopoDS_Face& face : {elbow.tube, elbow.plate})
      {
        const Handle(Geom_Surface) support = BRep_Tool::Surface(face);
        const bool                 tube    = face.IsSame(elbow.tube);
        const double               scale   = tube ? 1.0 / 20 : 1e-3 / 50;
        for (int k = 1; k < 32; ++k)
        {
          const double angle  = M_PI * k / 32;
          const gp_Pnt p      = tube ? gp_Pnt(20 * std::cos(angle), -20 * std::sin(angle), 0)
                                     : gp_Pnt(-13 + 26.0 * k / 32, 0, 50);
          gp_Dir       normal = tube ? gp_Dir(p.X(), p.Y(), 0) : gp_Dir(0, 0, 1);
          const auto   theirs = bends(support, p, normal);
          gp_Dir       own    = normal;
          const auto   mine   = bends(surface, p, own);
          const double tilt   = std::acos(std::min(1.0, own.Dot(normal))) * 180 / M_PI;
          EXPECT_LE(tilt, 0.1) << "flow " << flow << " at " << k;
          EXPECT_LE(std::abs(mine.first - theirs.first), 0.05 * scale)
            << "flow " << flow << " " << k;
          EXPECT_LE(std::abs(mine.second - theirs.second), 0.05 * scale)
            << "flow " << flow << " " << k;
        }
      }
      double u0, u1, v0, v1, skew = 0;
      BRepTools::UVBounds(fill, u0, u1, v0, v1);
      for (int i = 0; i <= 12; ++i)
        for (int j = 0; j <= 12; ++j)
        {
          const gp_Pnt p = surface->Value(u0 + (u1 - u0) * i / 12, v0 + (v1 - v0) * j / 12);
          skew           = std::max(
            skew,
            GeomAPI_ProjectPointOnSurf(gp_Pnt(-p.X(), p.Y(), p.Z()), surface).LowerDistance());
        }
      EXPECT_LE(skew, 0.2) << "flow " << flow << " sheet " << sheet;
    }
}

// What Patch's inside check makes of the elbow, asked directly. The elbow bends about five times
// as tightly where it turns as the tube does, the tightest thing on its rims, so against the rims
// alone it is refused; across its crossing from rim to rim it is held to the bridge between them
// and kept. The same face with a ripple of a hundredth of its size laid over its inside, six times
// across, by moving only poles three or more rows from either rim (the rims are checked to be
// unchanged to their second derivatives), is refused across the crossing all the same.
TEST(BRepPatchExt, InsideCheckKeepsTheElbowAcrossItsCrossingAndRefusesARippledOne)
{
  const Elbow                       elbow;
  const TopoDS_Face                 fill = elbow.fill(0);
  const Handle(Geom_BSplineSurface) surface =
    Handle(Geom_BSplineSurface)::DownCast(BRep_Tool::Surface(fill));
  ASSERT_FALSE(surface.IsNull());
  double u0, u1, v0, v1;
  surface->Bounds(u0, u1, v0, v1);
  // Which way the face crosses from the tube (z = 0) to the plate (z = 50).
  const bool acrossV = std::abs(surface->Value((u0 + u1) / 2, v0).Z()) < 1e-3;
  ASSERT_TRUE(acrossV || std::abs(surface->Value(u0, (v0 + v1) / 2).Z()) < 1e-3);
  const BRepExt::Crossing crossing{acrossV, 2, 2};
  const double            tube = 1.0 / 20;
  EXPECT_TRUE(BRepExt::insideFault(fill, tube).has_value());
  EXPECT_FALSE(BRepExt::insideFault(fill, tube, crossing).has_value());

  auto rippled = Handle(Geom_BSplineSurface)::DownCast(surface->Copy());
  for (int k = 1; k < 24; ++k)
  {
    rippled->InsertUKnot(u0 + (u1 - u0) * k / 24, 1, 1e-9, false);
    rippled->InsertVKnot(v0 + (v1 - v0) * k / 24, 1, 1e-9, false);
  }
  const int nu = rippled->NbUPoles(), nv = rippled->NbVPoles();
  for (int i = 1; i <= nu; ++i)
    for (int j = 1; j <= nv; ++j)
    {
      if ((acrossV ? std::min(j - 1, nv - j) : std::min(i - 1, nu - i)) < 3)
        continue;
      const double a = acrossV ? (i - 1.0) / (nu - 1) : (j - 1.0) / (nv - 1),
                   b = acrossV ? (j - 1.0) / (nv - 1) : (i - 1.0) / (nu - 1);
      const gp_Vec lift(0, 0, 0.64 * std::sin(2 * M_PI * 6 * a) * std::sin(M_PI * b));
      rippled->SetPole(i, j, rippled->Pole(i, j).Translated(lift));
    }
  for (int k = 0; k <= 16; ++k)
    for (const double rim : {acrossV ? v0 : u0, acrossV ? v1 : u1})
    {
      const double along = (acrossV ? u0 : v0) + (acrossV ? u1 - u0 : v1 - v0) * k / 16;
      const double u = acrossV ? along : rim, v = acrossV ? rim : along;
      gp_Pnt       p, q;
      gp_Vec       a[5], b[5];
      surface->D2(u, v, p, a[0], a[1], a[2], a[3], a[4]);
      rippled->D2(u, v, q, b[0], b[1], b[2], b[3], b[4]);
      ASSERT_LT(p.Distance(q), 1e-9);
      for (int d = 0; d < 5; ++d)
        ASSERT_LT((a[d] - b[d]).Magnitude(), 1e-6 * (1 + a[d].Magnitude()));
    }
  EXPECT_TRUE(BRepExt::insideFault(BRepBuilderAPI_MakeFace(rippled, 1e-7).Face(), tube, crossing)
                .has_value());
}

namespace
{

//! A loop of edges and the faces they are held G2 to.
struct Rounds
{
  std::vector<TopoDS_Face> faces;
  std::vector<TopoDS_Edge> edges;
  //! Per edge, the share of it either side of its middle on which the faces' bends are met.
  std::vector<double> middles;

  //! The edge of `face` whose middle is nearest `middle`.
  void hold(const TopoDS_Face& face, const gp_Pnt& middle, double share = 0.15)
  {
    TopoDS_Edge best;
    double      least = HUGE_VAL;
    for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next())
    {
      const BRepAdaptor_Curve curve(TopoDS::Edge(it.Current()));
      const double            distance =
        curve.Value((curve.FirstParameter() + curve.LastParameter()) / 2).Distance(middle);
      if (distance < least)
      {
        least = distance;
        best  = TopoDS::Edge(it.Current());
      }
    }
    faces.push_back(face);
    edges.push_back(best);
    middles.push_back(share);
  }

  BRepExt::PatchResult fill() const
  {
    std::vector<uint32_t> ids, supports;
    for (size_t i = 0; i < edges.size(); ++i)
    {
      ids.push_back(BRepExt::store(session(), edges[i]));
      supports.push_back(BRepExt::store(session(), faces[i]));
    }
    return BRepExt::patch(session(),
                          ids,
                          supports,
                          std::vector<int>(ids.size(), 2),
                          {},
                          {},
                          1e-3,
                          0.1,
                          0.05,
                          0);
  }

  //! Along every edge, 33 places: the fill passes through it within `tolerance` and keeps its
  //! face's tangent plane within what was asked; on its middle, away from the corners (30 % of it
  //! unless held with another share), its two principal curvatures are within 5 % of the tighter
  //! of the face's. All from projections onto both.
  void expectMet(const TopoDS_Face& fill, double tolerance = 1e-3) const
  {
    const Handle(Geom_Surface) surface = BRep_Tool::Surface(fill);
    for (size_t i = 0; i < edges.size(); ++i)
    {
      const BRepAdaptor_Curve    curve(edges[i]);
      const Handle(Geom_Surface) face = BRep_Tool::Surface(faces[i]);
      for (int k = 0; k <= 32; ++k)
      {
        const double               f = k / 32.0;
        const gp_Pnt               p = curve.Value(curve.FirstParameter()
                                     + (curve.LastParameter() - curve.FirstParameter()) * f);
        GeomAPI_ProjectPointOnSurf onto(p, surface);
        ASSERT_GT(onto.NbPoints(), 0);
        EXPECT_LE(onto.LowerDistance(), tolerance) << "edge " << i << " at " << f;
        if (k == 0 || k == 32)
          continue;
        gp_Dir     normal = gp_Dir(0, 0, 1), own;
        const auto theirs = bends(face, p, normal);
        own               = normal;
        const auto mine   = bends(surface, p, own);
        EXPECT_LE(std::acos(std::min(1.0, own.Dot(normal))) * 180 / M_PI, 0.1)
          << "edge " << i << " at " << f;
        if (std::abs(f - 0.5) > middles[i])
          continue;
        const double tightest = std::max(std::abs(theirs.first), std::abs(theirs.second));
        EXPECT_LE(std::abs(mine.first - theirs.first), 0.05 * tightest) << "edge " << i << " " << f;
        EXPECT_LE(std::abs(mine.second - theirs.second), 0.05 * tightest)
          << "edge " << i << " " << f;
      }
    }
  }
};

//! Three arms of quarter rounds of radius r about the axes, each from height r + q, with concave
//! quarter rounds of radius q between them (the teddy example), placed by `pose`. At each
//! of the six corners an arm, which does not bend along its axis, meets a strip whose edge bends
//! across that axis, and the other way round: the faces bend differently there.
Rounds teddy(double r, double q, const gp_Trsf& pose)
{
  Rounds       out;
  const double start = r + q, length = 60;
  for (int arm = 0; arm < 3; ++arm)
  {
    const auto to = [arm](double x, double y, double z) {
      return arm == 0 ? gp_Pnt(x, y, z) : arm == 1 ? gp_Pnt(y, z, x) : gp_Pnt(z, x, y);
    };
    const gp_Pnt o = to(0, 0, 0);
    const gp_Dir x(gp_Vec(o, to(1, 0, 0))), y(gp_Vec(o, to(0, 1, 0))), z(gp_Vec(o, to(0, 0, 1)));
    const auto   placed = [&](const TopoDS_Face& face) {
      return TopoDS::Face(BRepBuilderAPI_Transform(face, pose, true).Shape());
    };
    out.hold(placed(BRepBuilderAPI_MakeFace(new Geom_CylindricalSurface(gp_Ax3(o, z, x), r),
                                            0,
                                            M_PI / 2,
                                            start,
                                            start + length,
                                            1e-7)
                      .Face()),
             to(r * std::sqrt(0.5), r * std::sqrt(0.5), start).Transformed(pose));
    out.hold(
      placed(
        BRepBuilderAPI_MakeFace(
          new Geom_CylindricalSurface(gp_Ax3(to(start, 0, start), y.Reversed(), x.Reversed()), q),
          0,
          M_PI / 2,
          0,
          length,
          1e-7)
          .Face()),
      to(start - q * std::sqrt(0.5), 0, start - q * std::sqrt(0.5)).Transformed(pose));
  }
  return out;
}

//! A regular n-gon of circumradius 20 in z = 0, placed by `pose`, each side held to a round of
//! radius radii[i] tangent to the plane along it, above the plane with signs[i] 1 (the fill
//! rises inward) and below it with -1. The sides are straight and the rounds bend across them:
//! at every corner the faces bend differently. With `split` the first side is two edges, its
//! second half on a round of twice the radius: where they meet the edge runs straight on and the
//! two rounds bend differently across it.
Rounds pillow(int                        n,
              const std::vector<double>& radii,
              const std::vector<int>&    signs,
              const gp_Trsf&             pose,
              bool                       split = false)
{
  Rounds out;
  for (int i = split ? -1 : 0; i < n; ++i)
  {
    const int    k = std::max(i, 0);
    const double a = 2 * M_PI * k / n, b = 2 * M_PI * (k + 1) / n;
    gp_Pnt       p(20 * std::cos(a), 20 * std::sin(a), 0), q(20 * std::cos(b), 20 * std::sin(b), 0);
    if (split && i < 1)
      (i < 0 ? q : p) = gp_Pnt((p.XYZ() + q.XYZ()) / 2);
    const double radius = radii[k % radii.size()] * (split && i == 0 ? 2 : 1);
    const gp_Dir along(gp_Vec(p, q)), up(0, 0, signs[k % signs.size()]);
    // The round's parameter turns from its line on the plane toward along x -up: outward when
    // that is along x z, which points out of a loop that runs counterclockwise.
    const bool outward =
      gp_Vec(along.Crossed(up.Reversed())).Dot(gp_Vec(along.Crossed(gp_Dir(0, 0, 1)))) > 0;
    const TopoDS_Face face =
      BRepBuilderAPI_MakeFace(
        new Geom_CylindricalSurface(gp_Ax3(p.Translated(gp_Vec(up) * radius), along, up.Reversed()),
                                    radius),
        outward ? 0 : -0.6,
        outward ? 0.6 : 0,
        0,
        p.Distance(q),
        1e-7)
        .Face();
    out.hold(TopoDS::Face(BRepBuilderAPI_Transform(face, pose, true).Shape()),
             gp_Pnt((p.XYZ() + q.XYZ()) / 2).Transformed(pose),
             split && i < 1 ? 0 : 0.15);
  }
  return out;
}

//! The normal curvature of `surface` at the point nearest `p` in the tangent direction nearest
//! `direction`, for its normal turned to the side of `normal`.
double bendAlong(const Handle(Geom_Surface)& surface,
                 const gp_Pnt&               p,
                 const gp_Vec&               direction,
                 const gp_Dir&               normal)
{
  double u, v;
  GeomAPI_ProjectPointOnSurf(p, surface).LowerDistanceParameters(u, v);
  gp_Pnt at;
  gp_Vec su, sv, suu, svv, suv;
  surface->D2(u, v, at, su, sv, suu, svv, suv);
  gp_Vec n = su.Crossed(sv).Normalized();
  if (n.Dot(gp_Vec(normal)) < 0)
    n.Reverse();
  const double E = su.Dot(su), F = su.Dot(sv), G = sv.Dot(sv), det = E * G - F * F;
  const double du = (G * direction.Dot(su) - F * direction.Dot(sv)) / det,
               dv = (E * direction.Dot(sv) - F * direction.Dot(su)) / det;
  const gp_Vec t  = su * du + sv * dv;
  return n.Dot(suu * (du * du) + suv * (2 * du * dv) + svv * (dv * dv)) / t.SquareMagnitude();
}

gp_Trsf turned(double scale)
{
  gp_Trsf rotation, size, shift;
  rotation.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, -2, 0.5)), 0.9);
  size.SetScale(gp_Pnt(0, 0, 0), scale);
  shift.SetTranslation(gp_Vec(30, -12, 7));
  return shift * rotation * size;
}

} // namespace

// The teddy, turned in space, and with narrower strips at rest: at its six corners an arm and a
// strip bend differently, which no surface twice differentiable there can have both of. The fill
// is delivered and not refused: it meets every edge and both tangent planes along it, and the
// bend of each face on the middle of each edge, away from the corners (expectMet). What it misses
// near the corners is reported in the per-edge numbers, a share of the curvature above the 5 %
// asked.
TEST(BRepPatchExt, TeddyCornersThatBendDifferentlyAreEasedAndTheFacesMetBetween)
{
  for (const auto& [q, pose] : {std::pair<double, gp_Trsf>(50, turned(1)), {30, gp_Trsf()}})
  {
    SCOPED_TRACE(q);
    const Rounds      loop  = teddy(50, q, pose);
    const auto        built = loop.fill();
    const TopoDS_Face fill  = onlyFace(BRepExt::get(session(), built.evolution.resultId));
    ASSERT_FALSE(fill.IsNull());
    loop.expectMet(fill);
    for (size_t i = 0; i < loop.edges.size(); ++i)
      EXPECT_GT(built.values[5 + 4 * i + 2], 0.05) << "edge " << i;
  }
}

// Straight-sided loops of three to six sides on rounds tangent to their plane, of one radius or
// several, rising or falling, at rest or turned and scaled: every corner joins two straight edges
// whose rounds bend across them, so the one bend a surface can have there is none, and the
// rounds' is not. Each is delivered meeting every edge and tangent plane, and the rounds' bends
// on the middle of every edge; the corners' misses are reported. Four equal sides at rest are the
// untrimmed quad, the others fitted over their plane.
TEST(BRepPatchExt, StraightSidesOnRoundsAreEasedAtEveryCornerAndMeetTheRoundsBetween)
{
  struct Case
  {
    int                 n;
    std::vector<double> radii;
    std::vector<int>    signs;
    double              scale;
  };

  const Case cases[] = {{3, {30}, {1}, 0},
                        {4, {30}, {1}, 0},
                        {5, {30}, {1}, 0},
                        {6, {30}, {-1}, 3},
                        {4, {20, 60}, {1}, 1},
                        {5, {25, 40, 80}, {1, -1}, 1},
                        {3, {15, 45}, {1, 1, -1}, 0.5}};
  for (const Case& c : cases)
  {
    const gp_Trsf     pose  = c.scale > 0 ? turned(c.scale) : gp_Trsf();
    const Rounds      loop  = pillow(c.n, c.radii, c.signs, pose);
    const auto        built = loop.fill();
    const TopoDS_Face fill  = onlyFace(BRepExt::get(session(), built.evolution.resultId));
    ASSERT_FALSE(fill.IsNull()) << c.n;
    SCOPED_TRACE(c.n);
    loop.expectMet(fill);
    for (size_t i = 0; i < loop.edges.size(); ++i)
      EXPECT_GT(built.values[5 + 4 * i + 2], 0.05) << "edge " << i;
  }
}

// A square and a pentagon, at rest and turned, on rounds as above, with one side split at its
// middle between a round and one of twice the radius: at the split the edge runs straight on, its
// direction gives nothing of how the surface bends across it, and the bend asked across is the
// mean of the two rounds', as near to each as one bend is. The fill meets every edge, tangent
// plane and, at the middle of every edge, the round's bends; at the split its bend across the side
// is that mean, within 5 % of the tighter round's.
TEST(BRepPatchExt, ASideSplitBetweenTwoRoundsTakesTheMeanOfTheirBendsAcrossItWhereTheyMeet)
{
  for (const auto& [n, pose] : {std::pair<int, gp_Trsf>(4, gp_Trsf()), {5, turned(1)}})
  {
    SCOPED_TRACE(n);
    const Rounds      loop  = pillow(n, {30}, {1}, pose, true);
    const auto        built = loop.fill();
    const TopoDS_Face fill  = onlyFace(BRepExt::get(session(), built.evolution.resultId));
    ASSERT_FALSE(fill.IsNull());
    loop.expectMet(fill);
    const gp_Pnt split =
      gp_Pnt(10 * (1 + std::cos(2 * M_PI / n)), 10 * std::sin(2 * M_PI / n), 0).Transformed(pose);
    const gp_Vec across = gp_Vec(std::cos(M_PI / n), std::sin(M_PI / n), 0).Transformed(pose);
    const gp_Dir up     = gp_Dir(0, 0, 1).Transformed(pose);
    const double mean   = (bendAlong(BRep_Tool::Surface(loop.faces[0]), split, across, up)
                         + bendAlong(BRep_Tool::Surface(loop.faces[1]), split, across, up))
                        / 2;
    EXPECT_NEAR(bendAlong(BRep_Tool::Surface(fill), split, across, up), mean, 0.05 / 30);
  }
}

namespace
{

//! The rim of the corner of a box of 20 with its three edges at that corner rounded by 4: the
//! edges of the ball's eighth, each with the round it shares.
void cornerOfRounds(std::vector<uint32_t>& edges, std::vector<uint32_t>& supports)
{
  const TopoDS_Shape       box = BRepPrimAPI_MakeBox(20, 20, 20).Shape();
  BRepFilletAPI_MakeFillet rounds(box);
  for (const TopoDS_Edge& edge : edgesWhere(box, [](const double* b) {
         return b[3] > 19.9 && b[4] > 19.9 && b[5] > 19.9
                && (b[0] > 19.9) + (b[1] > 19.9) + (b[2] > 19.9) == 2;
       }))
    rounds.Add(4, edge);
  const TopoDS_Shape body = rounds.Shape();
  const auto         sphere =
    facesWhere(body, [](const double* b) { return b[0] > 15.9 && b[1] > 15.9 && b[2] > 15.9; });
  ASSERT_EQ(sphere.size(), 1u);
  for (TopExp_Explorer edge(sphere.front(), TopAbs_EDGE); edge.More(); edge.Next())
    for (const TopoDS_Face& face : facesWhere(body, [](const double*) { return true; }))
      if (!face.IsSame(sphere.front()))
        for (TopExp_Explorer other(face, TopAbs_EDGE); other.More(); other.Next())
          if (other.Current().IsSame(edge.Current()))
          {
            edges.push_back(BRepExt::store(session(), edge.Current()));
            supports.push_back(BRepExt::store(session(), face));
          }
}

} // namespace

// The corner of a box with its three edges rounded by 4, filled G2 with the three rounds at a
// tolerance of a millionth: no fill that eases at the corners meets the rounds beyond them that
// near, and the refusal names the corner where the rounds bend differently, not a miss of the
// net, and says that G1 can be had.
TEST(BRepPatchExt, ThreeRoundsAtAMillionthAreRefusedWithTheCornerWhereTheyBendDifferently)
{
  std::vector<uint32_t> edges, supports;
  cornerOfRounds(edges, supports);
  ASSERT_EQ(edges.size(), 3u);
  try
  {
    BRepExt::patch(session(), edges, supports, {2, 2, 2}, {}, {}, 1e-6, 0.1, 0.05, 0);
    FAIL() << "a G2 fill was returned";
  }
  catch (const BRepExt::Failure& failure)
  {
    const std::string what = failure.what();
    EXPECT_NE(what.find("bend differently"), std::string::npos) << what;
    EXPECT_NE(what.find("G1 can be had"), std::string::npos) << what;
  }
}

namespace
{

//! Held out from the design of the ease: an irregular polygon of the xy plane, its sides lifted
//! onto z = a x^2 + b y^2, each side held to a strip that leaves it outward tangent to that
//! surface but bends k across it, placed by `pose`. Every corner is then G1-compatible, and
//! where the strips' bends across differ from the surface's the faces bend differently there.
Rounds strips(const std::vector<gp_XY>&  polygon,
              double                     a,
              double                     b,
              const std::vector<double>& bends,
              const gp_Trsf&             pose)
{
  Rounds    out;
  const int n = static_cast<int>(polygon.size());
  for (int i = 0; i < n; ++i)
  {
    const gp_XY  p = polygon[i], q = polygon[(i + 1) % n];
    const double width = 0.3 * (q - p).Modulus(), k = bends[i % bends.size()];
    const auto   at = [&](double s, double t) {
      const gp_XY  xy = p + (q - p) * s;
      const gp_Pnt c(xy.X(), xy.Y(), a * xy.X() * xy.X() + b * xy.Y() * xy.Y());
      const gp_Vec normal = gp_Vec(-2 * a * xy.X(), -2 * b * xy.Y(), 1).Normalized();
      const gp_Vec along(q.X() - p.X(),
                         q.Y() - p.Y(),
                         2 * a * xy.X() * (q.X() - p.X()) + 2 * b * xy.Y() * (q.Y() - p.Y()));
      return c.Translated(along.Crossed(normal).Normalized() * t + normal * (k * t * t / 2))
        .Transformed(pose);
    };
    NCollection_Array2<gp_Pnt> grid(1, 17, 1, 5);
    for (int r = 1; r <= 17; ++r)
      for (int s = 1; s <= 5; ++s)
        grid(r, s) = at((r - 1) / 16.0, width * (s - 1) / 4);
    GeomAPI_PointsToBSplineSurface through;
    through.Interpolate(grid, Approx_IsoParametric);
    out.hold(BRepBuilderAPI_MakeFace(through.Surface(), 1e-7).Face(), at(0.5, 0));
  }
  return out;
}

} // namespace

// The family the ease was not designed on, on a dome and on a saddle, triangle and pentagon,
// at rest, mirrored, and turned at a thousandth and a thousand times the size (the tolerance in
// step): every one is delivered meeting its edges, tangent planes and, on the middle of every
// edge, the strips' bends, with the corners' misses reported. On a plane, where a strip that
// does not bend meets one that does, no eased fill meets the flat strip beyond its ease within
// the flat G2 limit: that is refused, naming the corner.
TEST(BRepPatchExt, HeldOutStripsOnDomesAndSaddlesAreEasedAtAnyScaleAndPose)
{
  const std::vector<std::vector<gp_XY>> polygons = {
    {{20, 0}, {-8, 17}, {-12, -15}},
    {{20, 0}, {8, 18}, {-14, 13}, {-17, -9}, {5, -19}}};
  const std::vector<double> bends = {0.05, -0.03, 0.0};
  gp_Trsf                   mirror, small, large;
  mirror.SetMirror(gp_Ax2(gp_Pnt(3, -2, 1), gp_Dir(1, 2, -1)));
  small.SetScale(gp_Pnt(0, 0, 0), 1e-3);
  large.SetScale(gp_Pnt(0, 0, 0), 1e3);
  const std::pair<gp_Trsf, double> poses[] = {{gp_Trsf(), 1},
                                              {mirror, 1},
                                              {turned(1) * small, 1e-3},
                                              {turned(1) * large, 1e3}};
  for (const auto& polygon : polygons)
    for (const auto& [a, b] : {std::pair<double, double>(0.012, 0.012), {0.012, -0.008}})
      for (const auto& [pose, scale] : poses)
      {
        SCOPED_TRACE(testing::Message() << polygon.size() << " sides, " << b << ", " << scale);
        Rounds                loop = strips(polygon, a, b, bends, pose);
        std::vector<uint32_t> edges, supports;
        for (size_t i = 0; i < loop.edges.size(); ++i)
        {
          edges.push_back(BRepExt::store(session(), loop.edges[i]));
          supports.push_back(BRepExt::store(session(), loop.faces[i]));
        }
        const auto        built = BRepExt::patch(session(),
                                          edges,
                                          supports,
                                          std::vector<int>(edges.size(), 2),
                                                 {},
                                                 {},
                                          1e-3 * scale,
                                          0.1,
                                          0.05,
                                          0);
        const TopoDS_Face fill  = onlyFace(BRepExt::get(session(), built.evolution.resultId));
        ASSERT_FALSE(fill.IsNull());
        loop.expectMet(fill, 1e-3 * scale);
        double reported = 0;
        for (size_t i = 0; i < edges.size(); ++i)
          reported = std::max(reported, built.values[5 + 4 * i + 2]);
        EXPECT_GT(reported, 0.05);
      }
  const Rounds          flat = strips(polygons[0], 0, 0, bends, gp_Trsf());
  std::vector<uint32_t> edges, supports;
  for (size_t i = 0; i < flat.edges.size(); ++i)
  {
    edges.push_back(BRepExt::store(session(), flat.edges[i]));
    supports.push_back(BRepExt::store(session(), flat.faces[i]));
  }
  try
  {
    BRepExt::patch(session(), edges, supports, {2, 2, 2}, {}, {}, 1e-3, 0.1, 0.05, 0);
    FAIL() << "the flat strips were filled";
  }
  catch (const BRepExt::Failure& failure)
  {
    EXPECT_NE(std::string(failure.what()).find("bend differently"), std::string::npos)
      << failure.what();
  }
}

namespace
{

//! A fill as it is asked of Patch: edges and the faces they are smooth with by their ids in the
//! session (0 for none), the order of each, the limits and how it is to be made.
struct Asking
{
  std::vector<uint32_t>  edges, supports;
  std::vector<int>       orders;
  double                 tolerance, g1Degrees, g2Share;
  int                    method   = 0;
  BRepExt::PatchControls controls = {};
};

//! `loop` asked as `Rounds::fill` asks it.
Asking asking(const Rounds& loop)
{
  Asking out{{}, {}, std::vector<int>(loop.edges.size(), 2), 1e-3, 0.1, 0.05};
  for (size_t i = 0; i < loop.edges.size(); ++i)
  {
    out.edges.push_back(BRepExt::store(session(), loop.edges[i]));
    out.supports.push_back(BRepExt::store(session(), loop.faces[i]));
  }
  return out;
}

//! What Patch makes of `asked` placed by `pose`, its tolerance scaled with it: whether it is
//! made, and why not with every number left out; what it made, its degrees and its poles.
struct Outcome
{
  bool                made = false;
  std::string         why;
  std::vector<double> built;
};

Outcome placed(const Asking& asked, const gp_Trsf& pose)
{
  BRep_Builder    builder;
  TopoDS_Compound all;
  builder.MakeCompound(all);
  for (const auto* ids : {&asked.edges, &asked.supports})
    for (const uint32_t id : *ids)
      if (id)
        builder.Add(all, BRepExt::get(session(), id));
  // Placed together, an edge stays the edge of its face.
  BRepBuilderAPI_Transform moved(all, pose, true);
  const auto               ids = [&](const std::vector<uint32_t>& given) {
    std::vector<uint32_t> out;
    for (const uint32_t id : given)
      out.push_back(id ? BRepExt::store(session(), moved.ModifiedShape(BRepExt::get(session(), id)))
                       : 0);
    return out;
  };
  Outcome out;
  try
  {
    const auto made = BRepExt::patch(session(),
                                     ids(asked.edges),
                                     ids(asked.supports),
                                     asked.orders,
                                     {},
                                     {},
                                     asked.tolerance * pose.ScaleFactor(),
                                     asked.g1Degrees,
                                     asked.g2Share,
                                     asked.method,
                                     -1,
                                     asked.controls);
    out.made        = true;
    out.built       = {made.values.begin(), made.values.begin() + 5};
  }
  catch (const BRepExt::Failure& failure)
  {
    out.why = std::regex_replace(failure.what(), std::regex("[0-9][0-9.]*(e[-+]?[0-9]+)?"), "#");
  }
  return out;
}

//! `at` is what `base` is: made or refused alike and for the same reason, by the same
//! construction of the same degrees with as many poles.
void expectSame(const Outcome& at, const Outcome& base)
{
  EXPECT_EQ(at.made, base.made) << (at.made ? base.why : at.why);
  EXPECT_EQ(at.why, base.why);
  EXPECT_EQ(at.built, base.built);
}

} // namespace

// The fills above, with every length and the tolerance a thousandth and a thousand times as
// large: Patch measures what it asks in the loop's size, so each is made or refused alike, by the
// same construction with as many poles. A length that is not a share of the size or of the
// tolerance shows as another net or another answer at one of the sizes. Not yet so, and asked at
// a thousand times only: at a thousandth the dish at G0 (its tolerance then OpenCascade's least
// distance) is refused as folding over, the corner of rounds (its tolerance below that) is made,
// the square on rounds misses G2 by its quad (a share of 0.109) and is trimmed instead, and the
// pentagon and the split side take other nets. The teddy is left out: it takes another net at a
// thousandth and is refused at a thousand times.
TEST(BRepPatchExt, FillsAreMadeTheSameAtAThousandthAndAThousandTimesTheirSize)
{
  struct Case
  {
    std::string         name;
    Asking              asked;
    std::vector<double> scales = {1e-3, 1e3};
  };

  std::vector<Case> cases;
  const uint32_t    circle = BRepExt::store(
    session(),
    BRepBuilderAPI_MakeEdge(new Geom_Circle(gp_Ax2(gp_Pnt(3, -2, 5), gp_Dir(0, 0, 1)), 12)).Edge());
  for (const int method : {0, 1})
    cases.push_back({"circle", {{circle}, {0}, {0}, 1e-6, 1, 0.1, method}});
  {
    const double       R = 10, latitude = 40 * M_PI / 180;
    const TopoDS_Shape sphere =
      BRepPrimAPI_MakeSphere(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), R, -M_PI / 2, latitude)
        .Shape();
    Asking cap{{}, {}, {}, 1e-4, 1, 0.1};
    for (const TopoDS_Edge& e : edgesWhere(sphere, [&](const double* b) {
           return near(b[2], R * std::sin(latitude)) && near(b[5], R * std::sin(latitude));
         }))
    {
      cap.edges.push_back(BRepExt::store(session(), e));
      cap.supports.push_back(BRepExt::store(session(), onlyFace(sphere)));
      cap.orders.push_back(1);
    }
    cases.push_back({"cap of a sphere", cap});
  }
  {
    Asking arcs{{}, {0, 0, 0}, {0, 0, 0}, 1e-6, 1, 0.1};
    for (const int k : {2, 1, 0})
    {
      const TopoDS_Edge arc =
        BRepBuilderAPI_MakeEdge(new Geom_Circle(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 9),
                                k * M_PI / 2,
                                (k + 1) * M_PI / 2)
          .Edge();
      arcs.edges.push_back(BRepExt::store(session(), k == 1 ? arc.Reversed() : TopoDS_Shape(arc)));
    }
    arcs.controls.open = true;
    cases.push_back({"three arcs", arcs});
  }
  for (const int order : {1, 0})
  {
    Asking              plates{{}, {}, std::vector<int>(6, order), 1e-4, 1, 0.1};
    std::vector<gp_Dir> normals;
    dish(plates.edges, plates.supports, normals);
    if (order)
      cases.push_back({"dish at G1", plates});
    else
      cases.push_back({"dish at G0", plates, {1e3}});
  }
  {
    const Elbow elbow;
    Asking      bend{elbow.edges, elbow.supports, {2, 2}, 1e-3, 0.1, 0.05};
    bend.controls.quad = 1;
    bend.controls.open = true;
    cases.push_back({"elbow", bend});
  }
  cases.push_back({"square on rounds", asking(pillow(4, {30}, {1}, gp_Trsf())), {1e3}});
  cases.push_back(
    {"pentagon on rounds", asking(pillow(5, {25, 40, 80}, {1, -1}, turned(1))), {1e3}});
  cases.push_back({"split side", asking(pillow(4, {30}, {1}, gp_Trsf(), true)), {1e3}});
  {
    Asking corner{{}, {}, {2, 2, 2}, 1e-6, 0.1, 0.05};
    cornerOfRounds(corner.edges, corner.supports);
    cases.push_back({"corner of rounds", corner, {1e3}});
  }
  for (const Case& c : cases)
  {
    SCOPED_TRACE(c.name);
    const Outcome base = placed(c.asked, gp_Trsf());
    for (const double scale : c.scales)
    {
      SCOPED_TRACE(scale);
      gp_Trsf pose;
      pose.SetScale(gp_Pnt(0, 0, 0), scale);
      expectSame(placed(c.asked, pose), base);
    }
  }
}

namespace
{

//! A circle of radius 10 in z = 0 on a cone that stands at `tilt` degrees against that plane
//! below it, narrowing upward, given as the rational spline it is.
Asking coneRim(double tilt)
{
  const occ::handle<Geom_ConicalSurface> cone =
    new Geom_ConicalSurface(gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, -1), gp_Dir(1, 0, 0)),
                            (90 - tilt) * M_PI / 180,
                            10);
  const TopoDS_Face face =
    BRepBuilderAPI_MakeFace(GeomConvert::SurfaceToBSplineSurface(
                              new Geom_RectangularTrimmedSurface(cone, 0.0, 2 * M_PI, 0.0, 4.0)),
                            1e-7)
      .Face();
  const auto rim = edgesWhere(face, [](const double* b) { return near(b[2], 0) && near(b[5], 0); });
  Asking     out{{BRepExt::store(session(), rim.at(0))},
                 {BRepExt::store(session(), face)},
                 {1},
             1e-4,
             0.1,
             0.05,
             1};
  out.controls.quad = 0;
  return out;
}

} // namespace

// The cap over a circle on a cone, G1 with it, the cone standing from 60 to 85 degrees against
// the circle's plane, at rest and turned at a thousandth, one and a thousand times the size,
// fitted over the plane. Above some tilt a face to be smooth with stands too steeply and the cap
// is refused at once: at every pose from the same tilt on and at none below it. Below it each
// tilt is made or refused alike at every pose, by the same net.
TEST(BRepPatchExt, TheCapOfASteepConeIsRefusedFromOneTiltOnAtEveryPoseAndSize)
{
  const double         tilts[] = {60, 66, 72, 76, 78, 80, 80.4, 82, 85};
  const gp_Trsf        poses[] = {gp_Trsf(), turned(1), turned(1e-3), turned(1e3)};
  std::vector<Asking>  rims;
  std::vector<Outcome> base;
  for (const double tilt : tilts)
  {
    rims.push_back(coneRim(tilt));
    base.push_back(placed(rims.back(), poses[0]));
  }
  const auto steep = [](const Outcome& outcome) {
    return outcome.why.find("stands too steeply") != std::string::npos;
  };
  ASSERT_FALSE(steep(base.front()));
  ASSERT_TRUE(steep(base.back()));
  for (size_t i = 1; i < base.size(); ++i)
    EXPECT_TRUE(!steep(base[i - 1]) || steep(base[i])) << tilts[i];
  for (size_t p = 1; p < std::size(poses); ++p)
    for (size_t i = 0; i < std::size(tilts); ++i)
    {
      SCOPED_TRACE(testing::Message() << "pose " << p << ", tilt " << tilts[i]);
      expectSame(placed(rims[i], poses[p]), base[i]);
    }
}
