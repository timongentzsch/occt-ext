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
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRep_Tool.hxx>
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomLProp_SLProps.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Circle.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Plane.hxx>

#include <gtest/gtest.h>

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
