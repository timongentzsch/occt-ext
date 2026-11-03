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
#include <BRepAdaptor_Curve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomAlgoExt_Failure.hxx>
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

//! The blend, or a null shape and the test failed with the refusal (so that the case is named).
TopoDS_Shape blended(const BRepExt::blend2::Request& request)
{
  try
  {
    return BRepExt::blend2::blend(request).shape;
  }
  catch (const std::exception& refused)
  {
    ADD_FAILURE() << refused.what();
    return {};
  }
}

//! A size law (Request::radii) as the reference modeller's variable blends take it, measured: the cubic spline
//! through the places with no curvature at its ends (periodic round a closed edge), by the
//! place's fraction of [first, last]. Built here by eliminating the slopes' continuity conditions
//! in turn (the Thomas algorithm), apart from the kernel's own solve; periodic by repeating the
//! period five times on each side, natural at the far ends, whose effect falls by 2 - sqrt 3 a
//! place: the periodic spline to rounding over the middle period.
std::function<double(double)> sizeLaw(const std::vector<std::pair<double, double>>& points,
                                      double                                        first,
                                      double                                        last,
                                      bool periodic = false)
{
  std::vector<double> x, y;
  const int           copies = periodic ? 11 : 1;
  for (int c = 0; c < copies; ++c)
    for (size_t i = 0; i < points.size(); ++i)
    {
      if (c > 0 && i == 0)
        continue;
      x.push_back(points[i].first + c - (periodic ? 5 : 0));
      y.push_back(points[i].second);
    }
  const size_t        n = x.size();
  std::vector<double> m(n, 0), c(n, 0), d(n, 0);
  // h_{i-1} m_{i-1} + 2 (h_{i-1} + h_i) m_i + h_i m_{i+1} = 6 (slope_i - slope_{i-1}), m = 0 at
  // ends.
  for (size_t i = 1; i + 1 < n; ++i)
  {
    const double h0 = x[i] - x[i - 1], h1 = x[i + 1] - x[i];
    const double rhs  = 6 * ((y[i + 1] - y[i]) / h1 - (y[i] - y[i - 1]) / h0);
    const double diag = 2 * (h0 + h1) - h0 * c[i - 1];
    c[i]              = h1 / diag;
    d[i]              = (rhs - h0 * d[i - 1]) / diag;
  }
  for (size_t i = n - 2; i >= 1; --i)
    m[i] = d[i] - c[i] * m[i + 1];
  return [=](double t) {
    double s = (t - first) / (last - first);
    if (periodic)
      s -= std::floor(s);
    size_t i = 0;
    while (i + 2 < n && s > x[i + 1])
      ++i;
    const double h = x[i + 1] - x[i], a = (x[i + 1] - s) / h, b = 1 - a;
    return a * y[i] + b * y[i + 1]
           + ((a * a * a - a) * m[i] + (b * b * b - b) * m[i + 1]) * h * h / 6;
  };
}

//! How far the faces of `made` that are not on `input` are from their discs: each sample P lies
//! in the section of the station station(P), on the circle about centre(t) of radius radius(t).
//! The largest |P - centre(t)| - radius(t), and the number of samples.
std::pair<double, int> discDefect(const TopoDS_Shape&                  made,
                                  const TopoDS_Shape&                  input,
                                  const std::function<double(gp_Pnt)>& station,
                                  const std::function<gp_Pnt(double)>& centre,
                                  const std::function<double(double)>& radius)
{
  TopoDS_Compound skin;
  BRep_Builder    builder;
  builder.MakeCompound(skin);
  for (TopExp_Explorer f(input, TopAbs_FACE); f.More(); f.Next())
    builder.Add(skin, f.Current());
  double worst = 0;
  int    count = 0;
  for (TopExp_Explorer f(made, TopAbs_FACE); f.More(); f.Next())
  {
    const TopoDS_Face face = TopoDS::Face(f.Current());
    double            u0, u1, v0, v1;
    BRepTools::UVBounds(face, u0, u1, v0, v1);
    const BRepAdaptor_Surface surface(face);
    // A face of the input, kept or trimmed, is not the blend: a point inside it lies on the input.
    bool onInput = false;
    for (int i = 1; i < 8 && !onInput; ++i)
    {
      const double u = u0 + (u1 - u0) * i / 8, v = v0 + (v1 - v0) * (8 - i) / 8;
      if (BRepClass_FaceClassifier(face, gp_Pnt2d(u, v), 1e-7).State() == TopAbs_IN)
        onInput =
          BRepExtrema_DistShapeShape(BRepBuilderAPI_MakeVertex(surface.Value(u, v)).Vertex(), skin)
            .Value()
          < 1e-6;
    }
    if (onInput)
      continue;
    constexpr int N = 16;
    for (int i = 1; i < N; ++i)
      for (int j = 1; j < N; ++j)
      {
        const double u = u0 + (u1 - u0) * i / N, v = v0 + (v1 - v0) * j / N;
        if (BRepClass_FaceClassifier(face, gp_Pnt2d(u, v), 1e-7).State() != TopAbs_IN)
          continue;
        const gp_Pnt p      = surface.Value(u, v);
        const double t      = station(p);
        const double defect = p.Distance(centre(t)) - radius(t);
        if (std::abs(defect) > std::abs(worst))
          worst = defect;
        ++count;
      }
  }
  return {worst, count};
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

// A size that changes along a box's edge (planes x = a, y = b, the edge along z), as the reference modeller
// makes it (its exported bodies, sectioned): each section across the edge is the circle of radius r(z)
// tangent to both faces, about (a - r, b - r, z), and r(z) the natural cubic spline through the
// places. The section removed is r^2 (1 - pi / 4), so the volume is (1 - pi / 4) times the integral
// of r^2, here in closed form: linear, and the natural splines through a peak and a valley, which
// are r = 1 + 6 s - 8 s^3 and 2 - 3.6 s + 4.8 s^3 on the first half and mirrored. At rest and
// turned, moved and scaled by a thousandth and a thousand, with the tolerances scaled.
TEST(BRepBlendExt, RoundWhoseSizeChangesHasCircularSectionsAcrossTheEdge)
{
  const double       a = 20, b = 20, c = 20;
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(a, b, c).Shape();
  const auto         up  = [&](const double* l) {
    return near(l[0], a) && near(l[1], b) && near(l[3], a) && near(l[4], b);
  };
  const auto half = [](double s, double r0, double k1, double k3) {
    const double x = s <= 0.5 ? s : 1 - s;
    return r0 + k1 * x + k3 * x * x * x;
  };
  const std::vector<
    std::pair<std::vector<std::pair<double, double>>, std::function<double(double)>>>
    laws = {{{{0, 1}, {1, 3}}, [](double s) { return 1 + 2 * s; }},
            {{{0, 0.5}, {1, 2.5}}, [](double s) { return 0.5 + 2 * s; }},
            {{{0, 1}, {0.5, 3}, {1, 1}}, [&](double s) { return half(s, 1, 6, -8); }},
            {{{0, 2}, {0.5, 0.8}, {1, 2}}, [&](double s) { return half(s, 2, -3.6, 4.8); }}};
  for (const auto& [points, law] : laws)
  {
    SCOPED_TRACE(points.size());
    const auto edges = edgesWhere(box, up);
    ASSERT_EQ(edges.size(), 1u);
    const BRepAdaptor_Curve curve(edges[0]);
    // The edge runs up z with its parameter.
    ASSERT_NEAR(curve.Value(curve.FirstParameter()).Z(), 0, 1e-12);
    const auto               r       = [&](double z) { return law(z / c); };
    BRepExt::blend2::Request request = roundOf(box, edges, 1);
    request.radii                    = {points};
    const TopoDS_Shape made          = blended(request);
    if (made.IsNull())
      continue;
    EXPECT_TRUE(BRepExt::invalidity(made).empty());
    const auto [worst, count] = discDefect(
      made,
      box,
      [](const gp_Pnt& p) { return p.Z(); },
      [&](double z) { return gp_Pnt(a - r(z), b - r(z), z); },
      r);
    EXPECT_GT(count, 50);
    EXPECT_LE(std::abs(worst), BRepExt::tolerance::sewing);
    // (1 - pi / 4) times the integral of r^2: Simpson's rule on 2000 spans, exact to rounding for
    // these polynomials of degree 6 at most on each half.
    double        integral = 0;
    constexpr int spans    = 2000;
    for (int k = 0; k <= spans; ++k)
      integral += (k == 0 || k == spans ? 1 : k % 2 ? 4 : 2) * std::pow(r(c * k / spans), 2);
    integral *= c / spans / 3;
    const double before = a * b * c - (1 - M_PI / 4) * integral;
    EXPECT_NEAR(volume(made), before, 1e-9 * a * b * c);
    for (const double s : {1.0, 1e-3, 1e3})
    {
      SCOPED_TRACE(s);
      const gp_Trsf      move   = posed(s);
      const TopoDS_Shape moved  = BRepBuilderAPI_Transform(box, move, true).Shape();
      const auto         turned = edgeThrough(moved, gp_Pnt(a, b, c / 2).Transformed(move), s * a);
      ASSERT_FALSE(turned.IsNull());
      // The turned edge may run the other way: its places then from its other end.
      const BRepAdaptor_Curve along(turned);
      const bool              flipped =
        along.Value(along.FirstParameter()).Distance(gp_Pnt(a, b, 0).Transformed(move))
        > 1e-9 * s * a;
      std::vector<std::pair<double, double>> mine;
      for (size_t i = 0; i < points.size(); ++i)
      {
        const auto& [place, size] = flipped ? points[points.size() - 1 - i] : points[i];
        mine.emplace_back(flipped ? 1 - place : place, size * s);
      }
      BRepExt::blend2::Request again = roundOf(moved, {turned}, s);
      again.radii                    = {mine};
      const TopoDS_Shape other       = blended(again);
      if (other.IsNull())
        continue;
      EXPECT_NEAR(volume(other) / (s * s * s), before, 1e-9 * a * b * c);
    }
  }
}

// The chamfer of the peak law keeps its feet r(z) from the edge on both faces, every section the
// line between them, (a - x) + (b - y) = r(z); a conic of rho 0.7 of the same law has its feet
// there too and meets the faces in their normals.
TEST(BRepBlendExt, ChamferAndConicWhoseSizeChangesKeepTheirFeet)
{
  const double       a = 20, b = 20, c = 20;
  const TopoDS_Shape box   = BRepPrimAPI_MakeBox(a, b, c).Shape();
  const auto         edges = edgesWhere(box, [&](const double* l) {
    return near(l[0], a) && near(l[1], b) && near(l[3], a) && near(l[4], b);
  });
  ASSERT_EQ(edges.size(), 1u);
  const std::vector<std::pair<double, double>> points = {{0, 1}, {0.5, 3}, {1, 1}};
  const auto                                   r      = [&](double z) {
    const double x = std::min(z / c, 1 - z / c);
    return 1 + 6 * x - 8 * x * x * x;
  };
  for (const auto law : {BRepExt::blend2::Law::Chamfer, BRepExt::blend2::Law::Conic})
  {
    BRepExt::blend2::Request request = roundOf(box, edges, 1);
    request.law                      = law;
    request.rho                      = 0.7;
    request.radii                    = {points};
    const TopoDS_Shape made          = blended(request);
    ASSERT_FALSE(made.IsNull());
    ASSERT_TRUE(BRepExt::invalidity(made).empty());
    double worst = 0, turn = 0;
    int    count = 0;
    for (const TopoDS_Face& face : facesWhere(made, [&](const double* l) {
           return l[3] - l[0] > 1e-6 && l[4] - l[1] > 1e-6 && l[3] - l[0] < 4;
         }))
    {
      double u0, u1, v0, v1;
      BRepTools::UVBounds(face, u0, u1, v0, v1);
      const BRepAdaptor_Surface surface(face);
      if (law == BRepExt::blend2::Law::Chamfer)
        for (int i = 1; i < 12; ++i)
          for (int j = 1; j < 12; ++j)
          {
            const gp_Pnt p = surface.Value(u0 + (u1 - u0) * i / 12, v0 + (v1 - v0) * j / 12);
            worst          = std::max(worst, std::abs((a - p.X()) + (b - p.Y()) - r(p.Z())));
            ++count;
          }
      // Along the contacts: on x = a the foot is (a, b - r, z), and the blend's normal x there.
      for (TopExp_Explorer e(face, TopAbs_EDGE); e.More(); e.Next())
      {
        const BRepAdaptor_Curve on(TopoDS::Edge(e.Current()));
        for (int k = 1; k < 16; ++k)
        {
          const gp_Pnt p =
            on.Value(on.FirstParameter() + (on.LastParameter() - on.FirstParameter()) * k / 16);
          const bool x = near(p.X(), a), y = near(p.Y(), b);
          if (x == y || p.Z() < 1e-6 || p.Z() > c - 1e-6)
            continue;
          worst = std::max(worst, std::abs((x ? b - p.Y() : a - p.X()) - r(p.Z())));
          ++count;
          if (law != BRepExt::blend2::Law::Conic)
            continue;
          GeomAPI_ProjectPointOnSurf onto(p, BRep_Tool::Surface(face));
          double                     u, v;
          onto.LowerDistanceParameters(u, v);
          BRepLProp_SLProps props(surface, u, v, 1, 1e-9);
          turn =
            std::max(turn, std::abs(std::abs(props.Normal().Dot(x ? gp::DX() : gp::DY())) - 1));
        }
      }
    }
    EXPECT_GT(count, 20);
    EXPECT_LE(worst, BRepExt::tolerance::sewing) << int(law);
    EXPECT_LE(turn, 1e-10) << "tangent to the faces";
  }
}

// Round a cylinder's top rim with a size that changes round it and ends as it begins (periodic):
// the section at the rim's angle t, in the half-plane through the axis, is the circle of radius
// r(t) about ((R - r) cos t, (R - r) sin t, H - r). A rim whose size does not end as it begins is
// refused by name.
TEST(BRepBlendExt, RoundOfARimWhoseSizeChangesHasCircularSections)
{
  const double       R = 10, H = 20;
  const TopoDS_Shape cylinder = BRepPrimAPI_MakeCylinder(R, H).Shape();
  const auto         rim =
    edgesWhere(cylinder, [&](const double* l) { return near(l[2], H) && near(l[5], H); });
  ASSERT_EQ(rim.size(), 1u);
  const std::vector<std::pair<double, double>> points = {{0, 1}, {0.25, 2}, {0.6, 0.7}, {1, 1}};
  const BRepAdaptor_Curve                      curve(rim[0]);
  const auto               r = sizeLaw(points, curve.FirstParameter(), curve.LastParameter(), true);
  BRepExt::blend2::Request request = roundOf(cylinder, rim, 1);
  request.radii                    = {points};
  const TopoDS_Shape made          = blended(request);
  ASSERT_FALSE(made.IsNull());
  EXPECT_TRUE(BRepExt::invalidity(made).empty());
  const auto angle = [&](const gp_Pnt& p) {
    const double t = std::atan2(p.Y(), p.X());
    return t < curve.FirstParameter() ? t + 2 * M_PI : t;
  };
  const auto [worst, count] = discDefect(
    made,
    cylinder,
    angle,
    [&](double t) { return gp_Pnt((R - r(t)) * std::cos(t), (R - r(t)) * std::sin(t), H - r(t)); },
    r);
  EXPECT_GT(count, 50);
  EXPECT_LE(std::abs(worst), BRepExt::tolerance::sewing);
  request.radii = {{{0, 1}, {1, 3}}};
  EXPECT_THROW(BRepExt::blend2::blend(request), BRepExt::Failure);
}

// Held out: a size that changes round the ring where a sphere of radius R stands on a plane, a
// concave edge between a plane and a sphere. The section at the ring's point E(t) lies in the
// half-plane through the axis and E(t): the circle of radius r about the point r above the plane
// and R + r from the sphere's centre, sqrt((R + r)^2 - r^2) from the axis.
TEST(BRepBlendExt, RoundOfASphereOnAPlaneWhoseSizeChangesHasCircularSections)
{
  const double       R = 14;
  const TopoDS_Shape standing =
    BRepAlgoAPI_Fuse(BRepPrimAPI_MakeBox(gp_Pnt(-25, -25, -10), 50, 50, 10).Shape(),
                     BRepPrimAPI_MakeSphere(gp_Pnt(0, 0, 0), R).Shape())
      .Shape();
  const auto ring = edgesWhere(standing, [](const double* b) {
    return near(b[2], b[5]) && near(b[2], 0) && b[3] - b[0] < 29 && b[3] - b[0] > 20;
  });
  ASSERT_EQ(ring.size(), 1u);
  const std::vector<std::pair<double, double>> points = {{0, 2}, {0.3, 3}, {0.7, 1.2}, {1, 2}};
  const BRepAdaptor_Curve                      curve(ring[0]);
  const auto               r = sizeLaw(points, curve.FirstParameter(), curve.LastParameter(), true);
  BRepExt::blend2::Request request = roundOf(standing, ring, 2);
  request.radii                    = {points};
  const TopoDS_Shape made          = blended(request);
  ASSERT_FALSE(made.IsNull());
  EXPECT_TRUE(BRepExt::invalidity(made).empty());
  // The ring's parameter at the angle of a point: the ring's point at that angle, found along it.
  const auto at = [&](const gp_Pnt& p) {
    const double want = std::atan2(p.Y(), p.X());
    double       best = curve.FirstParameter(), off = 1e9;
    for (int k = 0; k <= 3600; ++k)
    {
      const double t =
        curve.FirstParameter() + (curve.LastParameter() - curve.FirstParameter()) * k / 3600;
      const gp_Pnt E = curve.Value(t);
      const double d = std::abs(std::remainder(std::atan2(E.Y(), E.X()) - want, 2 * M_PI));
      if (d < off)
      {
        off  = d;
        best = t;
      }
    }
    // Refined by halving about it.
    double step = (curve.LastParameter() - curve.FirstParameter()) / 3600;
    for (int k = 0; k < 60; ++k, step /= 2)
      for (const double t : {best - step, best + step})
      {
        const gp_Pnt E = curve.Value(t);
        const double d = std::abs(std::remainder(std::atan2(E.Y(), E.X()) - want, 2 * M_PI));
        if (d < off)
        {
          off  = d;
          best = t;
        }
      }
    return best;
  };
  const auto [worst, count] = discDefect(
    made,
    standing,
    at,
    [&](double t) {
      const gp_Pnt E = curve.Value(t);
      const double a = std::atan2(E.Y(), E.X()), size = r(t);
      const double d = std::sqrt((R + size) * (R + size) - size * size);
      return gp_Pnt(d * std::cos(a), d * std::sin(a), size);
    },
    r);
  EXPECT_GT(count, 50);
  EXPECT_LE(std::abs(worst), BRepExt::tolerance::sewing);
}
