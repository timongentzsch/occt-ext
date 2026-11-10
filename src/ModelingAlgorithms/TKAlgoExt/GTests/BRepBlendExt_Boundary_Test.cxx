// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepBlendExt_BoundaryCurve.hxx>
#include <BRepBlendExt_BoundaryGraph.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom2d_BezierCurve.hxx>
#include <BRepTools.hxx>
#include <gtest/gtest.h>
#include <sstream>

using namespace BRepExt::blending::detail;

namespace
{
Handle(Geom_BSplineSurface) paraboloid()
{
  NCollection_Array2<gp_Pnt> poles(1, 3, 1, 2);
  for (int i = 1; i <= 3; ++i)
    for (int j = 1; j <= 2; ++j)
      poles(i, j) = gp_Pnt((i - 1) / 2., j - 1, (i == 2 ? -.25 : .25) + j - 1);
  NCollection_Array1<double> knots(1, 2);
  knots(1) = 0;
  knots(2) = 1;
  NCollection_Array1<int> um(1, 2), vm(1, 2);
  um.Init(3);
  vm.Init(2);
  return new Geom_BSplineSurface(poles, knots, knots, um, vm, 2, 1);
}

std::string serialized(const TopoDS_Shape& shape)
{
  std::ostringstream stream;
  BRepTools::Write(shape, stream);
  return stream.str();
}
} // namespace

TEST(BRepBlendExt, BoundaryCurveBoundContainsIndependentEvaluations)
{
  for (bool rational : {false, true})
    for (bool split : {false, true})
    {
      auto surface = paraboloid();
      if (rational)
        surface->SetWeight(2, 2, 0.7);
      if (split)
        surface->InsertUKnot(.43, 2, 0);
      const std::array<gp_Pnt2d, 4> pc = {gp_Pnt2d(.2, .1),
                                          gp_Pnt2d(.35, .3),
                                          gp_Pnt2d(.55, .25),
                                          gp_Pnt2d(.7, .2)};
      std::array<gp_Pnt, 4>         curve;
      NCollection_Array1<gp_Pnt>    cp(1, 4);
      NCollection_Array1<gp_Pnt2d>  uv(1, 4);
      for (int i = 0; i < 4; ++i)
      {
        curve[i]  = surface->Value(pc[i].X(), pc[i].Y());
        cp(i + 1) = curve[i];
        uv(i + 1) = pc[i];
      }
      const double bound = boundaryCurveGap(surface, curve, pc);
      ASSERT_TRUE(std::isfinite(bound));
      Geom_BezierCurve   c(cp);
      Geom2d_BezierCurve p(uv);
      double             sampled = 0;
      for (int i = 0; i <= 1024; ++i)
      {
        const double t = i / 1024.;
        const auto   q = p.Value(t);
        sampled        = std::max(sampled, c.Value(t).Distance(surface->Value(q.X(), q.Y())));
      }
      EXPECT_GT(sampled, 1e-4);
      EXPECT_LE(sampled, bound);
    }
}

TEST(BRepBlendExt, BoundaryCurveBoundDetectsDefectBetweenZeros)
{
  const auto surface = paraboloid();
  // A constant pcurve and a cubic displacement zero at 0, 1/2 and 1.
  std::array<gp_Pnt2d, 4> pc;
  pc.fill(gp_Pnt2d(.5, .5));
  std::array<gp_Pnt, 4> curve;
  curve.fill(surface->Value(.5, .5));
  curve[1].SetZ(curve[1].Z() + 1e-5);
  curve[2].SetZ(curve[2].Z() - 1e-5);
  const double bound = boundaryCurveGap(surface, curve, pc);
  EXPECT_GT(bound, 2e-6);
  EXPECT_LT(bound, 1.01e-5);
  pc[1].SetX(-.1);
  EXPECT_FALSE(std::isfinite(boundaryCurveGap(surface, curve, pc)));
}

TEST(BRepBlendExt, BoundaryCurveFitsSameContactAcrossPosesAndScales)
{
  for (double scale : {.001, 1., 1000.})
    for (bool rotated : {false, true})
    {
      auto    surface = paraboloid();
      gp_Trsf tr;
      tr.SetScale(gp::Origin(), scale);
      if (rotated)
      {
        gp_Trsf rotation;
        rotation.SetRotation(gp_Ax1(gp::Origin(), gp_Dir(1, 2, 3)), .73);
        gp_Trsf shift;
        shift.SetTranslation(gp_Vec(17 * scale, -23 * scale, 31 * scale));
        tr = shift * rotation * tr;
      }
      surface->Transform(tr);
      const gp_Pln plane = gp_Pln(gp_Pnt(0, 0, .25), gp::DZ()).Transformed(tr);
      const auto   face  = BRepBuilderAPI_MakeFace(plane, -1e5, 1e5, -1e5, 1e5).Face();
      Sheet        sheet;
      sheet.surface = surface;
      std::vector<gp_Pnt2d> guide;
      for (int i = 0; i <= 8; ++i)
      {
        const double u = i / 8.;
        guide.emplace_back(u, u * (1 - u));
      }
      const auto first  = BRepBuilderAPI_MakeVertex(surface->Value(0, 0)).Vertex(),
                 last   = BRepBuilderAPI_MakeVertex(surface->Value(1, 0)).Vertex();
      const auto result = planeBoundaryCurve(sheet, face, guide, first, last);
      ASSERT_TRUE(result);
      EXPECT_LE(result->sheetGap, tol::length);
      EXPECT_LE(result->planeGap, tol::length);
      double     f, l;
      const auto c  = BRep_Tool::Curve(result->edge, f, l);
      const auto s  = BRep_Tool::Surface(face);
      const auto pc = BRep_Tool::CurveOnSurface(result->edge, surface, TopLoc_Location(), f, l),
                 pp = BRep_Tool::CurveOnSurface(result->edge, face, f, l);
      for (int i = 0; i <= 2048; ++i)
      {
        const double t  = f + (l - f) * i / 2048;
        const auto   uv = pc->Value(t), planeUV = pp->Value(t);
        EXPECT_LE(c->Value(t).Distance(surface->Value(uv.X(), uv.Y())),
                  result->sheetGap + tol::rounding(100 * scale));
        EXPECT_LE(c->Value(t).Distance(s->Value(planeUV.X(), planeUV.Y())),
                  result->planeGap + tol::rounding(100 * scale));
      }
    }
}

TEST(BRepBlendExt, BoundaryGraphPreservesOccurrencesAndSeparateSeamPorts)
{
  const auto    a       = BRepBuilderAPI_MakeEdge(gp_Pnt(0, 0, 0), gp_Pnt(1, 0, 0)).Edge(),
                b       = BRepBuilderAPI_MakeEdge(gp_Pnt(1, 0, 0), gp_Pnt(0, 1, 0)).Edge(),
                c       = BRepBuilderAPI_MakeEdge(gp_Pnt(0, 1, 0), gp_Pnt(0, 0, 0)).Edge();
  const auto    reverse = TopoDS::Edge(a.Reversed());
  BoundaryGraph graph{{{a, 4, 5}, {c, 6, 4}, {b, 5, 6}, {a, 9, 10}, {reverse, 10, 9}}};
  const auto    loops = graph.loops();
  ASSERT_TRUE(loops);
  ASSERT_EQ(loops->size(), 2u);
  ASSERT_EQ((*loops)[0].size(), 3u);
  ASSERT_EQ((*loops)[1].size(), 2u);
  EXPECT_TRUE((*loops)[0][0].IsEqual(a));
  EXPECT_TRUE((*loops)[0][1].IsEqual(b));
  EXPECT_TRUE((*loops)[0][2].IsEqual(c));
  EXPECT_TRUE((*loops)[1][1].IsEqual(reverse));
  BRep_Builder builder;
  TopoDS_Wire  before, after;
  builder.MakeWire(before);
  builder.MakeWire(after);
  for (const auto& edge : {a, b, c})
    builder.Add(before, edge);
  for (const auto& edge : (*loops)[0])
    builder.Add(after, edge);
  before.Closed(true);
  after.Closed(true);
  EXPECT_EQ(serialized(before), serialized(after));
  graph.arcs.push_back({b, 4, 12});
  EXPECT_FALSE(graph.loops());
  graph.arcs.pop_back();
  graph.arcs.pop_back();
  EXPECT_FALSE(graph.loops());
}
