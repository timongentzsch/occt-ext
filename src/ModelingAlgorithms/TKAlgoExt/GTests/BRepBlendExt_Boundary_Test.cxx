// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepBlendExt_BoundaryCurve.hxx>
#include <BRepBlendExt_BoundaryGraph.hxx>
#include <BRepBlendExt_Attach.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_Line.hxx>
#include <Geom2d_BezierCurve.hxx>
#include <BRepTools.hxx>
#include <gtest/gtest.h>
#include <sstream>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>

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

TEST(BRepBlendExt, SharedSectionEnclosesBothParameterCharts)
{
  for (double scale : {.1, 1., 10.})
    for (bool moved : {false, true})
    {
      gp_Trsf pose;
      pose.SetScale(gp::Origin(), scale);
      if (moved)
      {
        gp_Trsf rotation, translation;
        rotation.SetRotation(gp_Ax1(gp::Origin(), gp_Dir(1, 2, 3)), .73);
        translation.SetTranslation(gp_Vec(17 * scale, -23 * scale, 31 * scale));
        pose = translation * rotation * pose;
      }
      auto surface = paraboloid();
      surface->Transform(pose);
      auto swapped = Handle(Geom_BSplineSurface)::DownCast(surface->Copy());
      swapped->ExchangeUV();
      std::array<Sheet, 2> sheets;
      sheets[0].surface = surface;
      sheets[1].surface = swapped;
      sheets[1].swapped = true;
      for (auto& sheet : sheets)
      {
        sheet.across[0] = 0;
        sheet.across[1] = 1;
      }
      const double u = .37, first = .13, last = .89;
      const auto   a      = BRepBuilderAPI_MakeVertex(surface->Value(u, first)).Vertex();
      const auto   b      = BRepBuilderAPI_MakeVertex(surface->Value(u, last)).Vertex();
      const auto   shared = sharedSectionCurve(sheets, {u, u}, first, last, a, b);
      ASSERT_TRUE(shared);
      EXPECT_LE(shared->gaps[0], tol::length);
      EXPECT_LE(shared->gaps[1], tol::length);
      double     f, l;
      const auto curve = BRep_Tool::Curve(shared->edge, f, l);
      for (int side = 0; side < 2; ++side)
      {
        const auto pc = BRep_Tool::CurveOnSurface(shared->edge, sheets[side].surface, {}, f, l);
        ASSERT_FALSE(pc.IsNull());
        for (int i = 0; i <= 512; ++i)
        {
          const double t  = first + (last - first) * i / 512.;
          const auto   uv = pc->Value(t);
          EXPECT_LE(curve->Value(t).Distance(sheets[side].surface->Value(uv.X(), uv.Y())),
                    side ? shared->gaps[1] : shared->gaps[0]);
        }
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

TEST(BRepBlendExt, CompleteCarrierRetainsContactChartsAtFiniteRestrictions)
{
  for (double scale : {.001, 1., 1000.})
    for (bool moved : {false, true})
    {
      SCOPED_TRACE(testing::Message() << scale << " " << moved);
      gp_Trsf pose;
      if (moved)
      {
        pose.SetRotation(gp_Ax1(gp::Origin(), gp_Dir(1, 2, 3)), .73);
        pose.SetTranslationPart(gp_Vec(17 * scale, -23 * scale, 31 * scale));
      }
      Piece piece;
      Rail& rail = piece.rail;
      rail.first = rail.from = -scale;
      rail.last = rail.to = scale;
      rail.scale = .2 * scale;
      rail.magnitude = 100 * scale;
      Handle(Geom_Curve) source = new Geom_Line(gp::Origin(), gp::DX());
      source->Transform(pose);
      rail.edge = BRepBuilderAPI_MakeEdge(source, -scale, scale);
      rail.curve.Initialize(rail.edge);
      for (int side = 0; side < 2; ++side)
      {
        NCollection_Array2<gp_Pnt> poles(1, 2, 1, 2);
        NCollection_Array1<double> uk(1, 2), vk(1, 2);
        NCollection_Array1<int> mult(1, 2);
        uk(1) = -2 * scale; uk(2) = 2 * scale;
        vk(1) = -scale; vk(2) = scale;
        mult.Init(2);
        for (int i = 1; i <= 2; ++i)
          for (int j = 1; j <= 2; ++j)
            poles(i, j) = (side == 0 ? gp_Pnt(uk(i), vk(j), 0)
                                      : gp_Pnt(uk(i), 0, vk(j))).Transformed(pose);
        rail.s[side].surface = new Geom_BSplineSurface(poles, uk, vk, mult, mult, 1, 1);
        rail.s[side].face = BRepBuilderAPI_MakeFace(rail.s[side].surface, tol::length);
        rail.s[side].distance = rail.scale;
        NCollection_Array1<gp_Pnt2d> feet(1, 2);
        NCollection_Array1<double> knots(1, 2);
        knots(1) = 0; knots(2) = 1;
        feet(1) = gp_Pnt2d(-scale, side == 0 ? rail.scale : -rail.scale);
        feet(2) = gp_Pnt2d(scale, feet(1).Y());
        piece.fitted.contact[side] = new Geom2d_BSplineCurve(feet, knots, mult, 1);
      }
      const auto completed = completeCarrier(piece);
      ASSERT_TRUE(completed);
      EXPECT_NEAR(completed->rail.spineLast, 4 * scale, 1e-8 * scale);
      for (int side = 0; side < 2; ++side)
        for (int end = 0; end < 2; ++end)
        {
          const auto uv = completed->fitted.contact[side]->Value(end);
          EXPECT_NEAR(uv.X(), (end ? 2 : -2) * scale, 1e-8 * scale);
          EXPECT_NEAR(uv.Y(), side == 0 ? rail.scale : -rail.scale, 1e-8 * scale);
          EXPECT_LE(completed->sheet.contact(side)->Value(completed->sheet.rail(end)).Distance(
                      rail.s[side].surface->Value(uv.X(), uv.Y())), 1e-7);
        }
      // The contact builder must preserve a genuine gap between two retained runs.
      std::vector<Piece> pieces = {*completed};
      std::vector<Corner> corners;
      Handle(BRepTools_History) history = new BRepTools_History;
      Attach attachment(rail.s[0].face, pieces, corners, history, false, nullptr);
      auto& intervals = attachment.contactIntervals[0][0];
      intervals.resize(2);
      const double ends[2][2] = {{.1, .3}, {.6, .9}};
      for (int run = 0; run < 2; ++run)
        for (int end = 0; end < 2; ++end)
        {
          auto& event = intervals[run].ends[end];
          event.parameter = ends[run][end];
          event.vertex = BRepBuilderAPI_MakeVertex(completed->sheet.contact(0)->Value(
            completed->sheet.rail(event.parameter))).Vertex();
        }
      attachment.plainContact(0, 0, rail.s[0].surface, TopLoc_Location());
      ASSERT_EQ(attachment.sideEdges(0, 0).size(), 2u);
      for (const auto& interval : intervals)
      {
        EXPECT_TRUE(TopExp::FirstVertex(interval.edge).IsSame(interval.ends[0].vertex));
        EXPECT_TRUE(TopExp::LastVertex(interval.edge).IsSame(interval.ends[1].vertex));
        double first, last;
        const auto curve = BRep_Tool::Curve(interval.edge, first, last);
        const auto pc = BRep_Tool::CurveOnSurface(interval.edge, rail.s[0].face, first, last);
        ASSERT_FALSE(pc.IsNull());
        for (int i = 0; i <= 16; ++i)
        {
          const double parameter = first + (last - first) * i / 16.;
          const auto uv = pc->Value(parameter);
          EXPECT_LE(curve->Value(parameter).Distance(
                      rail.s[0].surface->Value(uv.X(), uv.Y())), tol::length);
        }
      }
    }
}

TEST(BRepBlendExt, CarrierContinuationRetainsSharedBallAcrossSupportCharts)
{
  for (double scale : {.1, 1., 10.})
    for (bool moved : {false, true})
      for (auto law : {BRepExt::blending::Law::Conic, BRepExt::blending::Law::G2})
      {
        if (law == BRepExt::blending::Law::G2)
        {
          // A build without a supplied section policy has no G2 section to trace.
          try { (void)BREPBLENDEXT_POLICY::g2Poles(.5); }
          catch (const BRepExt::Failure&) { continue; }
        }
        SCOPED_TRACE(testing::Message() << scale << " " << moved << " " << int(law));
        gp_Trsf pose;
        if (moved)
        {
          pose.SetRotation(gp_Ax1(gp::Origin(), gp_Dir(1, 2, 3)), .73);
          pose.SetTranslationPart(gp_Vec(17 * scale, -23 * scale, 31 * scale));
        }
        Piece piece;
        Rail& rail = piece.rail;
        rail.first = rail.from = -scale;
        rail.last = rail.to = scale;
        rail.scale = .2 * scale;
        rail.magnitude = 100 * scale;
        rail.law = law;
        Handle(Geom_Curve) source = new Geom_Line(gp::Origin(), gp::DX());
        source->Transform(pose);
        rail.edge = BRepBuilderAPI_MakeEdge(source, -scale, scale);
        rail.curve.Initialize(rail.edge);
        const auto wall = [&](double first, double last) {
          NCollection_Array2<gp_Pnt> poles(1, 2, 1, 2);
          NCollection_Array1<double> uk(1, 2), vk(1, 2);
          NCollection_Array1<int> mult(1, 2);
          uk(1) = first; uk(2) = last;
          vk(1) = -scale; vk(2) = scale;
          mult.Init(2);
          for (int i = 1; i <= 2; ++i)
            for (int j = 1; j <= 2; ++j)
              poles(i, j) = gp_Pnt(uk(i), 0, vk(j)).Transformed(pose);
          return Handle(Geom_Surface)(new Geom_BSplineSurface(poles, uk, vk, mult, mult, 1, 1));
        };
        rail.s[0].surface = new Geom_Plane(gp_Pln(gp::Origin(), gp::DZ()).Transformed(pose));
        rail.s[0].face = BRepBuilderAPI_MakeFace(rail.s[0].surface, -6 * scale, 6 * scale,
                                               -2 * scale, 2 * scale, tol::length);
        rail.s[1].surface = wall(-2 * scale, 2 * scale);
        rail.s[1].face = BRepBuilderAPI_MakeFace(rail.s[1].surface, tol::length);
        for (int side = 0; side < 2; ++side)
        {
          rail.s[side].distance = rail.scale;
          NCollection_Array1<gp_Pnt2d> feet(1, 2);
          NCollection_Array1<double> knots(1, 2);
          NCollection_Array1<int> mult(1, 2);
          knots(1) = 0; knots(2) = 1;
          mult.Init(2);
          feet(1) = gp_Pnt2d(-scale, side == 0 ? rail.scale : -rail.scale);
          feet(2) = gp_Pnt2d(scale, feet(1).Y());
          piece.fitted.contact[side] = new Geom2d_BSplineCurve(feet, knots, mult, 1);
        }
        const auto completed = completeCarrier(piece);
        ASSERT_TRUE(completed);
        for (int end = 0; end < 2; ++end)
        {
          const auto next = BRepBuilderAPI_MakeFace(
            end ? wall(2 * scale, 4 * scale) : wall(-4 * scale, -2 * scale), tol::length).Face();
          const auto continued = continueCarrier(*completed, end, 1, next);
          ASSERT_TRUE(continued);
          EXPECT_NEAR(continued->rail.spineLast, 2 * scale, 1e-8 * scale);
          // Exact sheets can swap their axes and use an angular section parameter.
          const auto a = completed->sheet.section(completed->sheet.rail(end));
          const auto b = continued->sheet.section(continued->sheet.rail(1 - end));
          for (int i = 0; i <= 32; ++i)
          {
            const double w = i / 32.;
            EXPECT_LE(a->Value((1 - w) * completed->sheet.across[0] + w * completed->sheet.across[1]).Distance(
              b->Value((1 - w) * continued->sheet.across[0] + w * continued->sheet.across[1])), tol::length);
          }
          Station oldEvent, newEvent;
          for (int side = 0; side < 2; ++side)
          {
            auto uv = completed->fitted.contact[side]->Value(end);
            oldEvent.uv[side][0] = uv.X(); oldEvent.uv[side][1] = uv.Y();
            uv = continued->fitted.contact[side]->Value(1 - end);
            newEvent.uv[side][0] = uv.X(); newEvent.uv[side][1] = uv.Y();
          }
          ASSERT_TRUE(contactAt(completed->rail, completed->rail.t(end), oldEvent));
          ASSERT_TRUE(contactAt(continued->rail, continued->rail.t(1 - end), newEvent));
          EXPECT_LE(oldEvent.centre.Distance(newEvent.centre), tol::length);
          for (int side = 0; side < 2; ++side)
            EXPECT_LE(oldEvent.foot[side].Distance(newEvent.foot[side]), tol::length);
          for (int side = 0; side < 2; ++side)
            for (int i = 0; i <= 32; ++i)
            {
              const auto uv = continued->fitted.contact[side]->Value(i / 32.);
              EXPECT_LE(continued->sheet.contact(side)->Value(continued->sheet.rail(i / 32.)).Distance(
                continued->rail.s[side].surface->Value(uv.X(), uv.Y())), tol::length);
            }
        }
      }
}

// A common straight edge is represented with one or many spline spans. The
// surfaces have an exact, constant second-form jump along it; knot insertion
// must not turn a failed continuity claim into an accepted one.
TEST(BRepBlendExt, CurvatureAcceptanceIsIndependentOfSeamPartition)
{
  for (double scale : {.1, 1., 10.})
    for (double bend : {5e-5, 2.5e-4})
      for (int spans : {1, 64})
      {
        BRepBuilderAPI_Sewing sew;
        for (int side : {-1, 1})
        {
          NCollection_Array2<gp_Pnt> poles(1, 3, 1, 2);
          for (int i = 1; i <= 3; ++i)
            for (int j = 1; j <= 2; ++j)
              poles(i, j) = gp_Pnt(side * (i - 1) * scale / 2, (j - 1) * scale,
                                   side > 0 && i == 3 ? bend * scale / 2 : 0);
          NCollection_Array1<double> knots(1, 2);
          knots(1) = 0; knots(2) = 1;
          NCollection_Array1<int> um(1, 2), vm(1, 2);
          um.Init(3); vm.Init(2);
          Handle(Geom_BSplineSurface) surface =
            new Geom_BSplineSurface(poles, knots, knots, um, vm, 2, 1);
          sew.Add(BRepBuilderAPI_MakeFace(surface, 1e-7).Face());
        }
        sew.Perform();
        const auto body = sew.SewedShape();
        ASSERT_TRUE(BRepCheck_Analyzer(body).IsValid());
        NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
                                   TopTools_ShapeMapHasher> faces;
        TopExp::MapShapesAndAncestors(body, TopAbs_EDGE, TopAbs_FACE, faces);
        int shared = 0;
        BRep_Builder builder;
        for (int i = 1; i <= faces.Extent(); ++i)
        {
          const auto& on = faces.FindFromIndex(i);
          if (on.Extent() != 2) continue;
          ++shared;
          const auto edge = TopoDS::Edge(faces.FindKey(i));
          double first, last;
          const auto old = BRep_Tool::Curve(edge, first, last);
          NCollection_Array1<gp_Pnt> points(1, spans + 1);
          NCollection_Array1<double> knots(1, spans + 1);
          NCollection_Array1<int> mult(1, spans + 1);
          for (int k = 1; k <= spans + 1; ++k)
          {
            knots(k) = first + (last - first) * (k - 1) / spans;
            points(k) = old->Value(knots(k));
            mult(k) = k == 1 || k == spans + 1 ? 2 : 1;
          }
          Handle(Geom_BSplineCurve) curve = new Geom_BSplineCurve(points, knots, mult, 1);
          builder.UpdateEdge(edge, curve, BRep_Tool::Tolerance(edge));
          builder.Continuity(edge, TopoDS::Face(on.First()), TopoDS::Face(on.Last()), GeomAbs_G2);
        }
        ASSERT_EQ(shared, 1);
        ASSERT_TRUE(BRepCheck_Analyzer(body).IsValid());
        if (bend <= limit::seamBend)
          EXPECT_NO_THROW(BRepExt::blending::accepted("blend", body, body, scale));
        else
          EXPECT_THROW(BRepExt::blending::accepted("blend", body, body, scale), BRepExt::Failure);
      }
}
