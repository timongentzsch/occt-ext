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

// Whole-operation checks of TKAlgoExt against closed forms, through the host of
// BRepExtTest_Support.hxx. Every expected value below is derived from the input alone; none is read
// from the code under test.

#include "BRepExtTest_Support.hxx"

#include <BRepToolsExt_Operations.hxx>
#include <BRepPatchExt_Side.hxx>
#include <GeomAlgoExt_Failure.hxx>

#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Geom_BezierCurve.hxx>
#include <NCollection_Array1.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>

#include <gtest/gtest.h>

#include <cmath>
#include <ctime>
#include <vector>

using BRepExtTest::length;
using BRepExtTest::session;
using BRepExtTest::volume;

// A fillet of radius r on a straight convex box edge of length L removes the prism whose
// cross-section is the square r x r minus the quarter disc: area r^2 (1 - pi/4), so the volume
// is V - r^2 (1 - pi/4) L. A wrong contact, section or trim changes this by far more than the
// tolerance.
TEST(BRepToolsExtOperations, FilletOfBoxEdgeRemovesTheQuarterDiscComplement)
{
  for (const double scale : {0.1, 1.0, 10.0})
  {
    const double       a = 6 * scale, b = 4 * scale, c = 3 * scale, r = 0.5 * scale;
    const TopoDS_Shape box = BRepPrimAPI_MakeBox(a, b, c).Shape();
    TopoDS_Edge        edge;
    for (TopExp_Explorer it(box, TopAbs_EDGE); it.More(); it.Next())
    {
      const TopoDS_Edge candidate = TopoDS::Edge(it.Current());
      BRepAdaptor_Curve curve(candidate);
      const gp_Pnt      p0 = curve.Value(curve.FirstParameter()),
                   p1      = curve.Value(curve.LastParameter());
      if (std::abs(p0.X() - a) < 1e-9 && std::abs(p1.X() - a) < 1e-9 && std::abs(p0.Y() - b) < 1e-9
          && std::abs(p1.Y() - b) < 1e-9)
        edge = candidate;
    }
    ASSERT_FALSE(edge.IsNull());
    const uint32_t body = BRepExt::store(session(), box), picked = BRepExt::store(session(), edge);
    const auto     result =
      BRepExt::filletRange(session(), body, picked, r, 0, 1, true, {}, BRepExt::failureHashBound);
    const TopoDS_Shape& made = BRepExt::get(session(), result.resultId);
    EXPECT_TRUE(BRepCheck_Analyzer(made).IsValid());
    EXPECT_NEAR(volume(made), a * b * c - r * r * (1 - M_PI / 4) * c, 1e-4 * r * r * c);
  }
}

// Hollowing a box to a wall t leaves the volume abc - (a - 2t)(b - 2t)(c - 2t).
TEST(BRepToolsExtOperations, HollowOfBoxLeavesTheClosedFormWall)
{
  const double   a = 60, b = 40, c = 30, t = 2;
  const uint32_t body   = BRepExt::store(session(), BRepPrimAPI_MakeBox(a, b, c).Shape());
  const auto     result = BRepExt::hollow(session(), body, t, 1e-4, {}, BRepExt::failureHashBound);
  const TopoDS_Shape& made = BRepExt::get(session(), result.resultId);
  EXPECT_TRUE(BRepCheck_Analyzer(made).IsValid());
  EXPECT_NEAR(volume(made), a * b * c - (a - 2 * t) * (b - 2 * t) * (c - 2 * t), 1e-6 * a * b * c);
}

// The outward offset of a convex polygon by d with round corners has length P + 2 pi d: each
// side keeps its length and the corner arcs add one full turn of radius d.
TEST(BRepToolsExtOperations, RoundOffsetOfRectangleAddsOneTurnOfArc)
{
  const double a = 60, b = 40;
  const gp_Pnt p[4] = {gp_Pnt(0, 0, 0), gp_Pnt(a, 0, 0), gp_Pnt(a, b, 0), gp_Pnt(0, b, 0)};
  BRepBuilderAPI_MakeWire maker;
  for (int i = 0; i < 4; ++i)
    maker.Add(BRepBuilderAPI_MakeEdge(p[i], p[(i + 1) % 4]).Edge());
  const uint32_t wire = BRepExt::store(session(), maker.Wire());
  for (const double d : {2.0, 10.0})
  {
    const auto result = BRepExt::offsetCurve(session(), wire, d, 0);
    EXPECT_NEAR(length(BRepExt::get(session(), result.resultId)),
                2 * (a + b) + 2 * M_PI * d,
                1e-6 * (a + b));
  }
}

// A helix of radius R, `turns` turns and rise h has length turns * sqrt((2 pi R)^2 + (h /
// turns)^2).
TEST(BRepToolsExtOperations, HelixHasTheClosedFormLength)
{
  const double radius = 3, height = 8, turns = 4;
  const auto   result = BRepExt::spiral(session(),
                                        {0, 0, 0},
                                        {0, 0, 1},
                                        {1, 0, 0},
                                      radius,
                                      height,
                                      turns,
                                      0.0,
                                      false,
                                      1e-7);
  EXPECT_NEAR(length(BRepExt::get(session(), result.resultId)),
              turns * std::sqrt(std::pow(2 * M_PI * radius, 2) + std::pow(height / turns, 2)),
              1e-5);
}

// Four cubic Bezier sides of one analytic surface S(u,v) = (u, v, a (u^3 - 3 u v^2)). Patch
// returns a face whose sides lie on it: the distance of each side's points to the face is below
// the tolerance asked, measured here from the face's surface alone.
TEST(BRepToolsExtOperations, PatchOfFourCubicSidesContainsItsSides)
{
  const double a     = 0.5;
  const auto   point = [&](double u, double v) {
    return gp_Pnt(u, v, a * (u * u * u - 3 * u * v * v));
  };
  const auto side = [&](bool alongU, double fixed, bool reverse) {
    const auto   at = [&](double t) { return alongU ? point(t, fixed) : point(fixed, t); };
    const gp_Pnt p0 = at(0), p1 = at(1. / 3), p2 = at(2. / 3), p3 = at(1);
    NCollection_Array1<gp_Pnt> poles(1, 4);
    poles(1) = p0;
    poles(2) = gp_Pnt((-5 * p0.XYZ() + 18 * p1.XYZ() - 9 * p2.XYZ() + 2 * p3.XYZ()) / 6);
    poles(3) = gp_Pnt((2 * p0.XYZ() - 9 * p1.XYZ() + 18 * p2.XYZ() - 5 * p3.XYZ()) / 6);
    poles(4) = p3;
    Handle(Geom_BezierCurve) curve = new Geom_BezierCurve(poles);
    if (reverse)
      curve->Reverse();
    return BRepBuilderAPI_MakeEdge(curve).Edge();
  };
  const TopoDS_Edge     edges[4] = {side(true, 0, false),
                                    side(false, 1, false),
                                    side(true, 1, true),
                                    side(false, 0, true)};
  std::vector<uint32_t> ids;
  for (const auto& edge : edges)
    ids.push_back(BRepExt::store(session(), edge));
  const auto  result = BRepExt::patch(session(),
                                     ids,
                                     std::vector<uint32_t>(4, 0),
                                     std::vector<int>(4, 0),
                                      {},
                                      {},
                                     1e-6,
                                     1,
                                     0.1,
                                     0);
  TopoDS_Face face;
  for (TopExp_Explorer it(BRepExt::get(session(), result.evolution.resultId), TopAbs_FACE);
       it.More();
       it.Next())
    face = TopoDS::Face(it.Current());
  ASSERT_FALSE(face.IsNull());
  double worst = 0;
  for (const auto& edge : edges)
  {
    BRepAdaptor_Curve curve(edge);
    for (int k = 0; k <= 64; ++k)
    {
      const gp_Pnt               p = curve.Value(curve.FirstParameter()
                                   + (curve.LastParameter() - curve.FirstParameter()) * k / 64.);
      GeomAPI_ProjectPointOnSurf project(p, BRep_Tool::Surface(face));
      ASSERT_GT(project.NbPoints(), 0);
      worst = std::max(worst, project.LowerDistance());
    }
  }
  EXPECT_LE(worst, 1e-5);
}
