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

// Rolling-ball blends of BRepBlendExt against definitions that do not use the blend: a closed form
// for the volume where the body has one, and the ball defect (BRepExtTest_Support.hxx) where it has
// none. The bound on the defect is the default sewing tolerance, the tolerance a blend face is
// built to.

#include "BRepExtTest_Support.hxx"

#include <BRepBlendExt_Request.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepToolsExt_Checks.hxx>

#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>

#include <gtest/gtest.h>

using namespace BRepExtTest;

namespace
{

BRepExt::blend2::Request roundOf(const TopoDS_Shape&             body,
                                 const std::vector<TopoDS_Edge>& edges,
                                 double                          radius)
{
  BRepExt::blend2::Request request;
  request.operation = "round";
  request.body      = body;
  for (const TopoDS_Edge& edge : edges)
    request.edges.push_back(edge);
  request.first = radius;
  request.law   = BRepExt::blend2::Law::Conic;
  request.rho   = 0.5;
  return request;
}

//! A prism `height` high over the parallelogram with a corner at the origin and sides `side` long
//! at angles 0 and `angle` (below a half turn).
TopoDS_Shape wedge(double angle, double side, double height)
{
  const gp_Pnt a(side, 0, 0), c(side * std::cos(angle), side * std::sin(angle), 0);
  const gp_Pnt b(a.X() + c.X(), a.Y() + c.Y(), 0);
  return BRepPrimAPI_MakePrism(
           BRepBuilderAPI_MakeFace(BRepBuilderAPI_MakePolygon(gp::Origin(), a, b, c, true).Wire()),
           gp_Vec(0, 0, height))
    .Shape();
}

//! Faces of one surface merged, as a modeller leaves a Boolean's result.
TopoDS_Shape unified(const TopoDS_Shape& shape)
{
  ShapeUpgrade_UnifySameDomain unify(shape, true, true, true);
  unify.Build();
  return unify.Shape();
}

//! The edge of `shape` through `at`.
TopoDS_Edge edgeThrough(const TopoDS_Shape& shape, const gp_Pnt& at, double size)
{
  for (const TopoDS_Edge& edge : edgesWhere(shape, [](const double*) { return true; }))
  {
    const BRepExtrema_DistShapeShape d(BRepBuilderAPI_MakeVertex(at).Vertex(), edge);
    if (d.IsDone() && d.Value() < 1e-9 * size)
      return edge;
  }
  return {};
}

//! A rigid motion and a uniform scale `s`: turned about a skew axis and moved off the origin.
gp_Trsf posed(double s)
{
  gp_Trsf turn, shift, scale;
  turn.SetRotation(gp_Ax1(gp::Origin(), gp_Dir(1, 2, 3)), 0.7);
  shift.SetTranslation(gp_Vec(3, -5, 7));
  scale.SetScale(gp::Origin(), s);
  return scale * shift * turn;
}

} // namespace

// Every edge of a box rounded at once: the result is the parallel body of the inner box, whose
// volume is Steiner's polynomial. Eight corners where three blends meet are eighths of a sphere,
// and a wrong corner patch, a wrong trim or a blend that stops short of the corner changes the
// volume by far more than the tolerance (the corner patch alone is 4/3 pi r^3 of it).
TEST(BRepBlendExt, RoundingEveryEdgeOfABoxHasSteinersVolume)
{
  const double       a = 40, b = 30, c = 20;
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(a, b, c).Shape();
  for (const double r : {6.0, 9.0})
  {
    const TopoDS_Shape made =
      BRepExt::blend2::blend(roundOf(box, edgesWhere(box, [](const double*) { return true; }), r))
        .shape;
    EXPECT_TRUE(BRepExt::invalidity(made).empty());
    EXPECT_NEAR(volume(made), roundedBoxVolume(a, b, c, r), 1e-6 * volume(box)) << "r=" << r;
  }
}

// The rolling-ball definition on a body with no closed form: a sphere standing on a box, the ring
// where they meet rounded. The ring is a space curve between a plane and a sphere; the ball's
// centre must stay at distance r from both along all of it.
TEST(BRepBlendExt, RoundOfASphereOnABoxRollsAtDistanceR)
{
  const TopoDS_Shape standing =
    BRepAlgoAPI_Fuse(BRepPrimAPI_MakeBox(gp_Pnt(-25, -25, -10), 50, 50, 10).Shape(),
                     BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), 14).Shape())
      .Shape();
  const auto ring = [](const double* b) {
    return near(b[2], b[5]) && near(b[2], 0) && b[3] - b[0] < 29 && b[3] - b[0] > 20;
  };
  for (const double r : {2.0, 6.0})
  {
    const auto edges = edgesWhere(standing, ring);
    ASSERT_FALSE(edges.empty());
    const TopoDS_Shape made = BRepExt::blend2::blend(roundOf(standing, edges, r)).shape;
    EXPECT_TRUE(BRepExt::invalidity(made).empty());
    const Defect defect = ballDefect(made, standing, r);
    ASSERT_GT(defect.samples, 0);
    EXPECT_LE(std::abs(defect.worst), BRepExt::tolerance::sewing) << "r=" << r;
  }
}

// Two perpendicular cylinders fused: the edges are the intersection curves of two curved surfaces,
// where the contact curves of the blend change with the position along the edge.
TEST(BRepBlendExt, RoundOfCrossedCylindersRollsAtDistanceR)
{
  const TopoDS_Shape along =
    BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -20), gp_Dir(0, 0, 1)), 10, 40).Shape();
  const TopoDS_Shape across =
    BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-20, 0, 0), gp_Dir(1, 0, 0)), 6, 40).Shape();
  const TopoDS_Shape crossed = BRepAlgoAPI_Fuse(along, across).Shape();
  const auto         inside  = [](const double* b) {
    return b[3] - b[0] < 39 && b[5] - b[2] < 39 && b[3] - b[0] > 1;
  };
  for (const double r : {1.0, 2.5})
  {
    const auto edges = edgesWhere(crossed, inside);
    ASSERT_FALSE(edges.empty());
    const TopoDS_Shape made = BRepExt::blend2::blend(roundOf(crossed, edges, r)).shape;
    EXPECT_TRUE(BRepExt::invalidity(made).empty());
    const Defect defect = ballDefect(made, crossed, r);
    ASSERT_GT(defect.samples, 0);
    EXPECT_LE(std::abs(defect.worst), BRepExt::tolerance::sewing) << "r=" << r;
  }
}

// A blend on a body that already carries blends: the top edges rounded by 8, then the uprights by
// 3, so the second blend runs out onto the first one's surfaces. The ball rolls where it does not
// reach the bottom face (z = 0).
TEST(BRepBlendExt, ReblendOfRoundedEdgesRollsAtDistanceR)
{
  const TopoDS_Shape box      = BRepPrimAPI_MakeBox(40, 30, 20).Shape();
  const auto         top      = [](const double* b) { return near(b[2], 20) && near(b[5], 20); };
  const auto         uprights = [](const double* b) {
    return near(b[0], b[3]) && near(b[1], b[4]) && b[5] - b[2] > 5;
  };
  const TopoDS_Shape first = BRepExt::blend2::blend(roundOf(box, edgesWhere(box, top), 8)).shape;
  const auto         edges = edgesWhere(first, uprights);
  ASSERT_FALSE(edges.empty());
  const TopoDS_Shape second = BRepExt::blend2::blend(roundOf(first, edges, 3)).shape;
  EXPECT_TRUE(BRepExt::invalidity(second).empty());
  const Defect defect =
    ballDefect(second, first, 3, [](const gp_Pnt& centre) { return centre.Z() >= 3; });
  ASSERT_GT(defect.samples, 0);
  EXPECT_LE(std::abs(defect.worst), BRepExt::tolerance::sewing);
}

// A straight edge between two planes at any angle: the ball of radius r in a wedge of opening psi
// removes (convex) or adds (concave) the section r^2 (cot(psi / 2) - (pi - psi) / 2) along the
// whole edge, and the faces at its ends, square to it, stop it. In a sharp wedge the ball's
// centre is r / sin(psi / 2) from the edge, more than four times r below 29 degrees, where the
// blend's carry past the ends once stopped at the ends and the blend did not reach the faces that
// stop it. Held to the closed form at rest and turned, moved and scaled by a thousandth and a
// thousand, with the tolerance scaled.
TEST(BRepBlendExt, RoundInAWedgeOfAnyAngleHasTheClosedFormSection)
{
  const double side = 10, height = 10;
  for (const double degrees : {5., 10., 20., 28., 30., 60., 90., 120., 150., 170., 175.})
    for (const bool convex : {true, false})
    {
      const double psi = degrees * M_PI / 180;
      // Concave: a block with a wedge of opening psi cut out to its middle.
      gp_Trsf below;
      below.SetTranslation(gp_Vec(0, 0, -1));
      const TopoDS_Shape body =
        convex
          ? wedge(psi, side, height)
          : unified(
              BRepAlgoAPI_Cut(
                BRepPrimAPI_MakeBox(gp_Pnt(-side, -side, 0), 2 * side, 2 * side, height).Shape(),
                BRepBuilderAPI_Transform(wedge(psi, 4 * side, height + 2), below, true).Shape())
                .Shape());
      // As large as fits the faces: the contacts are r cot(psi / 2) from the edge.
      const double r       = std::min(1.0, 0.3 * std::tan(psi / 2) * side);
      const double section = r * r * (1 / std::tan(psi / 2) - (M_PI - psi) / 2);
      const double before  = volume(body);
      for (const double s : {1.0, 1e-3, 1e3})
        for (const bool turned : {false, true})
        {
          const gp_Trsf      move  = turned ? posed(s) : gp_Trsf();
          const TopoDS_Shape moved = BRepBuilderAPI_Transform(body, move, true).Shape();
          const TopoDS_Edge  edge =
            edgeThrough(moved, gp_Pnt(0, 0, height / 2).Transformed(move), s * side);
          ASSERT_FALSE(edge.IsNull());
          if (!turned && s != 1)
            continue;
          const TopoDS_Shape made = BRepExt::blend2::blend(roundOf(moved, {edge}, r * s)).shape;
          EXPECT_TRUE(BRepExt::invalidity(made).empty()) << degrees << (convex ? " convex" : "");
          EXPECT_NEAR(volume(made) / (s * s * s),
                      before + (convex ? -1 : 1) * section * height,
                      1e-9 * before)
            << degrees << (convex ? " convex" : " concave") << " scale " << s;
        }
    }
}

// The same sharp and flat wedges between a plane and a cylinder and between two cylinders (lines
// along the cylinders' axes, so the angle stays the same along the edge), not used to find the
// cause: the ball rolls at distance r from both faces.
TEST(BRepBlendExt, RoundInAWedgeOfACylinderRollsAtDistanceR)
{
  const double R = 10, height = 10;
  for (const double degrees : {10., 20., 160.})
    for (const bool convex : {true, false})
      for (const bool cylinders : {false, true})
      {
        const double       psi = degrees * M_PI / 180;
        const double       phi = convex ? psi : 2 * M_PI - psi; // the angle inside the body
        TopoDS_Shape       body;
        gp_Pnt             on;
        const TopoDS_Shape round = BRepPrimAPI_MakeCylinder(R, height).Shape();
        if (cylinders)
        {
          // Axes D apart: the faces meet at the angle whose normals are D / R apart.
          const double D =
            convex ? 2 * R * std::cos(phi / 2) : 2 * R * std::sin(phi / 2 - M_PI / 2);
          const TopoDS_Shape other =
            BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(D, 0, convex ? -1 : 0), gp::DZ()),
                                     R,
                                     convex ? height + 2 : height)
              .Shape();
          body = convex ? BRepAlgoAPI_Common(round, other).Shape()
                        : BRepAlgoAPI_Fuse(round, other).Shape();
          on   = gp_Pnt(D / 2, std::sqrt(R * R - D * D / 4), height / 2);
        }
        else
        {
          // The plane y = h through the cylinder, at the point where the circle's angle is t.
          const double       t = convex ? phi - M_PI / 2 : 1.5 * M_PI - phi, h = R * std::sin(t);
          const TopoDS_Shape slab = BRepPrimAPI_MakeBox(gp_Pnt(-2 * R, -2 * R, convex ? -1 : 0),
                                                        4 * R,
                                                        h + 2 * R,
                                                        convex ? height + 2 : height)
                                      .Shape();
          body = convex ? BRepAlgoAPI_Common(round, slab).Shape()
                        : BRepAlgoAPI_Fuse(round, slab).Shape();
          on   = gp_Pnt(std::sqrt(R * R - h * h), h, height / 2);
        }
        body                   = unified(body);
        const double      r    = 0.3 * std::min(1.0, 4 * std::pow(std::sin(psi / 2), 2));
        const TopoDS_Edge edge = edgeThrough(body, on, R);
        ASSERT_FALSE(edge.IsNull()) << degrees;
        const TopoDS_Shape made = BRepExt::blend2::blend(roundOf(body, {edge}, r)).shape;
        EXPECT_TRUE(BRepExt::invalidity(made).empty()) << degrees;
        // Where the ball rolls: between the faces that stop the edge.
        const Defect defect = ballDefect(made, body, r, [&](const gp_Pnt& centre) {
          return centre.Z() > r && centre.Z() < height - r;
        });
        ASSERT_GT(defect.samples, 0) << degrees;
        EXPECT_LE(std::abs(defect.worst), BRepExt::tolerance::sewing)
          << degrees << (convex ? " convex" : " concave") << (cylinders ? " cylinders" : " plane");
      }
}
