// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Rolling-ball blends of BRepBlendExt against definitions that do not use the blend: a closed form
// for the volume where the body has one, and the ball defect (BRepExtTest_Support.hxx) where it has
// none. The bound on the defect is the default sewing tolerance, the tolerance a blend face is
// built to.

#include "BRepExtTest_Support.hxx"

#include <BRepBlendExt_Request.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepToolsExt_Checks.hxx>

#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>

#include <gtest/gtest.h>

using namespace BRepExtTest;

namespace {

BRepExt::blend2::Request roundOf(const TopoDS_Shape& body, const std::vector<TopoDS_Edge>& edges, double radius) {
  BRepExt::blend2::Request request;
  request.operation = "round";
  request.body = body;
  for (const TopoDS_Edge& edge : edges) request.edges.push_back(edge);
  request.first = radius;
  request.law = BRepExt::blend2::Law::Conic;
  request.rho = 0.5;
  return request;
}

} // namespace

// Every edge of a box rounded at once: the result is the parallel body of the inner box, whose
// volume is Steiner's polynomial. Eight corners where three blends meet are eighths of a sphere,
// and a wrong corner patch, a wrong trim or a blend that stops short of the corner changes the
// volume by far more than the tolerance (the corner patch alone is 4/3 pi r^3 of it).
TEST(BRepBlendExt, RoundingEveryEdgeOfABoxHasSteinersVolume) {
  const double a = 40, b = 30, c = 20;
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(a, b, c).Shape();
  for (const double r : {6.0, 9.0}) {
    const TopoDS_Shape made =
        BRepExt::blend2::blend(roundOf(box, edgesWhere(box, [](const double*) { return true; }), r)).shape;
    EXPECT_TRUE(BRepExt::invalidity(made).empty());
    EXPECT_NEAR(volume(made), roundedBoxVolume(a, b, c, r), 1e-6 * volume(box)) << "r=" << r;
  }
}

// The rolling-ball definition on a body with no closed form: a sphere standing on a box, the ring
// where they meet rounded. The ring is a space curve between a plane and a sphere; the ball's
// centre must stay at distance r from both along all of it.
TEST(BRepBlendExt, RoundOfASphereOnABoxRollsAtDistanceR) {
  const TopoDS_Shape standing = BRepAlgoAPI_Fuse(BRepPrimAPI_MakeBox(gp_Pnt(-25, -25, -10), 50, 50, 10).Shape(),
                                                 BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), 14).Shape()).Shape();
  const auto ring = [](const double* b) { return near(b[2], b[5]) && near(b[2], 0) && b[3] - b[0] < 29 && b[3] - b[0] > 20; };
  for (const double r : {2.0, 6.0}) {
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
TEST(BRepBlendExt, RoundOfCrossedCylindersRollsAtDistanceR) {
  const TopoDS_Shape along = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -20), gp_Dir(0, 0, 1)), 10, 40).Shape();
  const TopoDS_Shape across = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-20, 0, 0), gp_Dir(1, 0, 0)), 6, 40).Shape();
  const TopoDS_Shape crossed = BRepAlgoAPI_Fuse(along, across).Shape();
  const auto inside = [](const double* b) { return b[3] - b[0] < 39 && b[5] - b[2] < 39 && b[3] - b[0] > 1; };
  for (const double r : {1.0, 2.5}) {
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
TEST(BRepBlendExt, ReblendOfRoundedEdgesRollsAtDistanceR) {
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(40, 30, 20).Shape();
  const auto top = [](const double* b) { return near(b[2], 20) && near(b[5], 20); };
  const auto uprights = [](const double* b) { return near(b[0], b[3]) && near(b[1], b[4]) && b[5] - b[2] > 5; };
  const TopoDS_Shape first = BRepExt::blend2::blend(roundOf(box, edgesWhere(box, top), 8)).shape;
  const auto edges = edgesWhere(first, uprights);
  ASSERT_FALSE(edges.empty());
  const TopoDS_Shape second = BRepExt::blend2::blend(roundOf(first, edges, 3)).shape;
  EXPECT_TRUE(BRepExt::invalidity(second).empty());
  const Defect defect = ballDefect(second, first, 3, [](const gp_Pnt& centre) { return centre.Z() >= 3; });
  ASSERT_GT(defect.samples, 0);
  EXPECT_LE(std::abs(defect.worst), BRepExt::tolerance::sewing);
}
