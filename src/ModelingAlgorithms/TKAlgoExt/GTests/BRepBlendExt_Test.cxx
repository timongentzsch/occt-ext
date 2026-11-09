// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Rolling-ball blends of BRepBlendExt against definitions that do not use the blend: a closed form
// for the volume where the body has one, and the ball defect (BRepExtTest_Support.hxx) where it has
// none. The bound on the defect is the default sewing tolerance, the tolerance a blend face is
// built to.

#include "BRepExtTest_Support.hxx"

#include <BRepBlendExt_Request.hxx>
#include <BRepBlendExt_Smooth.hxx>
#include <BRepBlendExt_Internal.hxx>
#include <GeomConvert.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Ellipse.hxx>
#include <IntAna_QuadQuadGeo.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_SurfaceOfLinearExtrusion.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>

#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomAlgoExt_Failure.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>

#include <gtest/gtest.h>
#include <BRepTools.hxx>
#include <sstream>

using namespace BRepExtTest;

namespace
{

BRepExt::blending::Request roundOf(const TopoDS_Shape&             body,
                                   const std::vector<TopoDS_Edge>& edges,
                                   double                          radius)
{
  BRepExt::blending::Request request;
  request.operation = "round";
  request.body      = body;
  for (const TopoDS_Edge& edge : edges)
    request.edges.push_back(edge);
  request.first = radius;
  request.law   = BRepExt::blending::Law::Conic;
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
TopoDS_Shape blended(const BRepExt::blending::Request& request)
{
  try
  {
    return BRepExt::blending::blend(request).shape;
  }
  catch (const std::exception& refused)
  {
    ADD_FAILURE() << refused.what();
    return {};
  }
}

//! A size law (Request::radii): the cubic spline
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
      BRepExt::blending::blend(roundOf(box, edgesWhere(box, [](const double*) { return true; }), r))
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
    const TopoDS_Shape made = BRepExt::blending::blend(roundOf(standing, edges, r)).shape;
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
    const TopoDS_Shape made = BRepExt::blending::blend(roundOf(crossed, edges, r)).shape;
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
  const TopoDS_Shape first = BRepExt::blending::blend(roundOf(box, edgesWhere(box, top), 8)).shape;
  const auto         edges = edgesWhere(first, uprights);
  ASSERT_FALSE(edges.empty());
  const TopoDS_Shape second = BRepExt::blending::blend(roundOf(first, edges, 3)).shape;
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
          const TopoDS_Shape made = BRepExt::blending::blend(roundOf(moved, {edge}, r * s)).shape;
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
        const TopoDS_Shape made = BRepExt::blending::blend(roundOf(body, {edge}, r)).shape;
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

// A size that changes along a box's edge (planes x = a, y = b, the edge along z): each section across the edge is the circle of radius r(z)
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
    BRepExt::blending::Request request = roundOf(box, edges, 1);
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
      BRepExt::blending::Request again = roundOf(moved, {turned}, s);
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
  for (const auto law : {BRepExt::blending::Law::Chamfer, BRepExt::blending::Law::Conic})
  {
    BRepExt::blending::Request request = roundOf(box, edges, 1);
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
      if (law == BRepExt::blending::Law::Chamfer)
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
          if (law != BRepExt::blending::Law::Conic)
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
  BRepExt::blending::Request request = roundOf(cylinder, rim, 1);
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
  EXPECT_THROW(BRepExt::blending::blend(request), BRepExt::Failure);
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
  BRepExt::blending::Request request = roundOf(standing, ring, 2);
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

// A hole of radius 3 beside the top edge x = 0 of a block, a gap g from the face x = 0. Tangent
// (g = 0) the body is pinched: the hole meets the face in a line of four faces, which no blend
// passes; it is refused by name. With a gap the round is made and rolls
// at distance r from both faces where its ball stays clear of the hole and of the faces that stop
// the edge.
TEST(BRepBlendExt, RoundBesideAHoleTangentToItsFaceIsRefusedByName)
{
  for (const double gap : {0.0, 0.25, 1.0})
  {
    const TopoDS_Shape body =
      BRepAlgoAPI_Cut(
        BRepPrimAPI_MakeBox(20, 20, 20).Shape(),
        BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(3 + gap, 10, -5), gp::DZ()), 3, 30).Shape())
        .Shape();
    const auto edges = edgesWhere(body, [](const double* b) {
      return near(b[2], 20) && near(b[5], 20) && b[3] < 1e-6;
    });
    ASSERT_FALSE(edges.empty());
    const double r = 2;
    if (gap == 0)
    {
      try
      {
        BRepExt::blending::blend(roundOf(body, edges, r));
        ADD_FAILURE() << "a pinched body was blended";
      }
      catch (const BRepExt::Failure& refused)
      {
        EXPECT_NE(std::string(refused.what()).find("pinched"), std::string::npos) << refused.what();
      }
      continue;
    }
    const TopoDS_Shape made = BRepExt::blending::blend(roundOf(body, edges, r)).shape;
    EXPECT_TRUE(BRepExt::invalidity(made).empty()) << gap;
    const Defect defect = ballDefect(made, body, r, [&](const gp_Pnt& c) {
      return std::hypot(c.X() - 3 - gap, c.Y() - 10) > 3 + r && c.Y() > r && c.Y() < 20 - r;
    });
    ASSERT_GT(defect.samples, 0);
    EXPECT_LE(std::abs(defect.worst), BRepExt::tolerance::sewing) << gap;
  }
}

// The top rim of a cylinder of radius R whose side is h high, rounded or chamfered by r >= h: the
// side is used up (the blend is one face between the two caps). The round's
// ball rests on the top and on the bottom rim, its centre (R - sqrt(r^2 - (r - h)^2), h - r) in the
// plane through the axis, and the chamfer's bevel from the top's foot at 45 degrees runs on to the
// bottom, which cuts it. Each removes pi times the integral over z of R^2 - x(z)^2, x(z) the
// section's distance from the axis: in closed form for the bevel, and for the arc by its angle (z =
// c + r sin t), where the integrand is smooth. h = r is the side exactly as high as the ball is
// wide.
TEST(BRepBlendExt, BlendOfARimOverASideNarrowerThanItUsesTheSideUp)
{
  const double R = 10;
  for (const auto& [h, r] :
       std::vector<std::pair<double, double>>{{1, 1}, {1, 3}, {2, 3}, {0.5, 1}})
    for (const bool chamfer : {false, true})
    {
      const TopoDS_Shape cylinder = BRepPrimAPI_MakeCylinder(R, h).Shape();
      const auto         rim =
        edgesWhere(cylinder, [&](const double* b) { return near(b[2], h) && near(b[5], h); });
      ASSERT_EQ(rim.size(), 1u);
      BRepExt::blending::Request request = roundOf(cylinder, rim, r);
      if (chamfer)
        request.law = BRepExt::blending::Law::Chamfer;
      const TopoDS_Shape made = blended(request);
      if (made.IsNull())
        continue;
      EXPECT_TRUE(BRepExt::invalidity(made).empty()) << h << " " << r << " " << chamfer;
      int faces = 0;
      for (TopExp_Explorer f(made, TopAbs_FACE); f.More(); f.Next())
        ++faces;
      EXPECT_EQ(faces, 3) << h << " " << r << " " << chamfer;
      double removed = 0;
      if (chamfer)
      {
        // x(z) = R - r + (h - z) for z down to 0: pi times the integral of R^2 - x^2.
        const auto F = [&](double z) {
          const double x = R - r + h - z;
          return R * R * z + x * x * x / 3;
        };
        removed = M_PI * (F(h) - F(0));
      }
      else
      {
        const double  cz = h - r, cx = R - std::sqrt(r * r - cz * cz);
        const double  t0    = std::asin(-cz / r);
        constexpr int spans = 2000;
        for (int i = 0; i <= spans; ++i)
        {
          const double t = t0 + (M_PI / 2 - t0) * i / spans, x = cx + r * std::cos(t);
          removed += (i == 0 || i == spans ? 1 : i % 2 ? 4 : 2) * (R * R - x * x) * r * std::cos(t);
        }
        removed *= M_PI * (M_PI / 2 - t0) / spans / 3;
      }
      // The round's sheet is exact (a torus); the bevel's edge on the bottom is interpolated
      // through its crossings with the cone to the length tolerance, so the volume is within that
      // over the bevel's area (at most the cone's from the top's foot to the side's).
      EXPECT_NEAR(volume(made),
                  M_PI * R * R * h - removed,
                  chamfer ? BRepExt::tolerance::length * M_PI * (2 * R - r) * r * M_SQRT2
                          : 1e-9 * M_PI * R * R * h)
        << h << " " << r << " " << chamfer;
    }
}

// The rim of a hole of radius R leaning a through a plate t thick (the slanted-rim family), rounded by r. The balls' centres lie r below the top and R + r from the
// hole's axis: on an ellipse of semi-axes R + r and (R + r) / cos a, whose radius of curvature at
// its acute end is (R + r) cos a. Below that the round is the balls' envelope, every point r from
// the ellipse (at 75 degrees its sections turn by up to 165 degrees and are made of two spans);
// past r/(R+r) >= cos(a), a compact fixed-contact conic replaces its acute sector.
// The retained sector must still satisfy the original rolling-ball bound.
TEST(BRepBlendExt, RoundOfASteepRimIsTheBallsWhereItDoesNotCrossItself)
{
  const auto text = [](const TopoDS_Shape& shape) {
    std::ostringstream stream;
    BRepTools::Write(shape, stream);
    return stream.str();
  };
  struct Member
  {
    double t, R, r;
  };

  for (const double degrees : {70.0, 75.0, 80.0})
    for (const Member& m :
         {Member{4, 1, 0.5}, Member{4, 2, 0.5}, Member{8, 1, 0.3}, Member{8, 2, 1}})
    {
      const double a = degrees * M_PI / 180;
      const double size =
        std::max(20.0, 2 * std::ceil(m.t / 2 * std::tan(a) + m.R / std::cos(a) + 3));
      const double height = std::max(4 * m.t + 8 * m.R, 2 * m.t / std::cos(a) + 8 * m.R);
      gp_Trsf      tilt;
      tilt.SetRotation(gp_Ax1(gp_Pnt(size / 2, size / 2, m.t / 2), gp::DX()), a);
      const TopoDS_Shape hole =
        BRepBuilderAPI_Transform(
          BRepPrimAPI_MakeCylinder(
            gp_Ax2(gp_Pnt(size / 2, size / 2, m.t / 2 - height / 2), gp::DZ()),
            m.R,
            height)
            .Shape(),
          tilt,
          true)
          .Shape();
      const TopoDS_Shape plate =
        BRepAlgoAPI_Cut(BRepPrimAPI_MakeBox(size, size, m.t).Shape(), hole).Shape();
      std::vector<TopoDS_Edge> rim;
      for (const TopoDS_Edge& e :
           edgesWhere(plate, [&](const double* b) { return b[5] > m.t - 1e-6; }))
      {
        const BRepAdaptor_Curve c(e);
        bool                    top = c.GetType() != GeomAbs_Line;
        for (const double s : {0.1, 0.5, 0.9})
          top = top
                && c.Value(c.FirstParameter() + (c.LastParameter() - c.FirstParameter()) * s).Z()
                     > m.t - 1e-6;
        if (top)
          rim.push_back(e);
      }
      ASSERT_FALSE(rim.empty());
      const bool crosses = m.r / (m.R + m.r) >= std::cos(a);
      SCOPED_TRACE(std::to_string(degrees) + " t" + std::to_string(m.t) + " R"
                   + std::to_string(m.R));
      const auto inputBefore = text(plate);
      try
      {
        const auto result = BRepExt::blending::blend(roundOf(plate, rim, m.r));
        const TopoDS_Shape made = result.shape;
        EXPECT_EQ(inputBefore, text(plate));
        EXPECT_NO_THROW(BRepExt::blending::accepted("round", plate, made, m.r, &result));
        EXPECT_TRUE(BRepExt::invalidity(made).empty());
        // Every point of the blend r from the centres' ellipse: about the axis through (size/2,
        // size/2, t/2) along (0, -sin a, cos a), the circle of radius R + r in the plane z = t - r.
        const gp_Pnt O(size / 2, size / 2, m.t / 2);
        const gp_Dir d(0, -std::sin(a), std::cos(a));
        const gp_Dir e1 = gp::DX(), e2 = d.Crossed(e1);
        const auto   centre = [&](double q) {
          const gp_Vec radial =
            gp_Vec(e1) * std::cos(q) * (m.R + m.r) + gp_Vec(e2) * std::sin(q) * (m.R + m.r);
          const double along = (m.t - m.r - O.Z() - radial.Z()) / d.Z();
          return O.Translated(gp_Vec(d) * along + radial);
        };
        double worst = 0;
        int    count = 0;
        for (const TopoDS_Face& face : facesWhere(made, [](const double*) { return true; }))
        {
          const BRepAdaptor_Surface surface(face);
          if (surface.GetType() == GeomAbs_Plane || surface.GetType() == GeomAbs_Cylinder)
            continue;
          double u0, u1, v0, v1;
          BRepTools::UVBounds(face, u0, u1, v0, v1);
          for (int i = 1; i < 16; ++i)
            for (int j = 1; j < 16; ++j)
            {
              const double u = u0 + (u1 - u0) * i / 16, v = v0 + (v1 - v0) * j / 16;
              if (BRepClass_FaceClassifier(face, gp_Pnt2d(u, v), 1e-7).State() != TopAbs_IN)
                continue;
              if (crosses)
              {
                // Read the station phase from the cylindrical contact, independently of
                // the spline's longitudinal coordinate and of closest-centre projection.
                const gp_Lin axis(O, d);
                const gp_Pnt p0 = surface.Value(u, v0), p1 = surface.Value(u, v1);
                const gp_Pnt q = std::abs(axis.Distance(p0)-m.R) < std::abs(axis.Distance(p1)-m.R)
                                   ? p0 : p1;
                const gp_Vec arm(O, q);
                const double phase = std::atan2(arm.Dot(gp_Vec(e1)), arm.Dot(gp_Vec(e2)));
                if (std::abs(phase) <= 2 * std::asin(m.r / (m.R + m.r)))
                  continue;
              }
              const gp_Pnt p    = surface.Value(u, v);
              const auto   gap  = [&](double q) { return p.Distance(centre(q)); };
              int          best = 0;
              for (int k = 1; k < 720; ++k)
                if (gap(k * M_PI / 360) < gap(best * M_PI / 360))
                  best = k;
              double lo = (best - 1) * M_PI / 360, hi = (best + 1) * M_PI / 360;
              for (int k = 0; k < 100; ++k)
              {
                const double x = lo + (hi - lo) / 3, y = hi - (hi - lo) / 3;
                if (gap(x) < gap(y))
                  hi = y;
                else
                  lo = x;
              }
              worst = std::max(worst, std::abs(gap((lo + hi) / 2) - m.r));
              ++count;
            }
        }
        EXPECT_GT(count, 50);
        EXPECT_LE(worst, 1e-6 * std::max(1.0, m.r));
      }
      catch (const BRepExt::Failure& refused)
      {
        ADD_FAILURE() << refused.what();
      }
    }
}

// The rail can continue beyond a trimmed, multi-span edge. Compare its integral and inverse
// with the closed form for (t,t^2): F(t) = t sqrt(1+4t^2)/2 + asinh(2t)/4.
TEST(BRepBlendExt, ExtendedRailArcLength)
{
  using namespace BRepExt::blending::detail;
  const auto primitive = [](double t) {
    return t * std::sqrt(1 + 4 * t * t) / 2 + std::asinh(2 * t) / 4;
  };
  for (double scale : {0.001, 1.0, 1000.0})
    for (bool turned : {false, true})
    {
      NCollection_Array1<gp_Pnt> poles(1, 3);
      poles(1) = gp_Pnt(0, 0, 0);
      poles(2) = gp_Pnt(0.5, 0, 0);
      poles(3) = gp_Pnt(1, 1, 0);
      auto curve = GeomConvert::CurveToBSplineCurve(new Geom_BezierCurve(poles));
      curve->InsertKnot(0.5, 1);
      gp_Trsf transform;
      if (turned)
        transform.SetRotation(gp_Ax1(gp::Origin(), gp_Dir(1, 2, 3)), 0.7);
      gp_Trsf scaling;
      scaling.SetScale(gp::Origin(), scale);
      curve->Transform(scaling * transform);
      Rail rail;
      rail.scale = scale;
      rail.magnitude = 4 * scale;
      rail.edge = BRepBuilderAPI_MakeEdge(curve, 0.25, 0.75);
      rail.curve.Initialize(rail.edge);
      for (double start : {-0.5, 0.0, 0.25, 0.5, 0.75, 1.0, 1.5})
        for (double end : {-0.5, 0.0, 0.25, 0.5, 0.75, 1.0, 1.5})
        {
          SCOPED_TRACE(::testing::Message() << scale << " " << turned << " " << start << " " << end);
          const double distance = scale * (primitive(end) - primitive(start));
          EXPECT_NEAR(railLength(rail, start, end), std::abs(distance), scale * 1e-7);
          const auto parameter = railParameterAt(rail, start, distance);
          ASSERT_TRUE(parameter.has_value());
          EXPECT_NEAR(*parameter, end, 1e-6);
        }
    }
}

// Rational held-out family: a circular arc written as a quadratic rational spline. Length is
// radius times the difference of atan2 angles of the exact evaluated endpoints, independent of
// the spline parameter and the integrator. Its denominator stays positive on the tested extension.
TEST(BRepBlendExt, ExtendedRationalRailArcLength)
{
  using namespace BRepExt::blending::detail;
  for (double radius : {0.001, 1.0, 1000.0})
  {
    NCollection_Array1<gp_Pnt> poles(1, 3);
    poles(1) = gp_Pnt(radius, 0, 0);
    poles(2) = gp_Pnt(radius, radius, 0);
    poles(3) = gp_Pnt(0, radius, 0);
    NCollection_Array1<double> weights(1, 3), knots(1, 2);
    weights(1) = weights(3) = 1;
    weights(2) = std::sqrt(0.5);
    knots(1) = 0;
    knots(2) = 1;
    NCollection_Array1<int> mults(1, 2);
    mults(1) = mults(2) = 3;
    Handle(Geom_BSplineCurve) curve = new Geom_BSplineCurve(poles, weights, knots, mults, 2);
    curve->InsertKnot(0.5, 1);
    Rail rail;
    rail.scale = radius;
    rail.magnitude = 2 * radius;
    rail.edge = BRepBuilderAPI_MakeEdge(curve, 0.25, 0.75);
    rail.curve.Initialize(rail.edge);
    for (double start : {-0.5, 0.0, 0.25, 0.5, 0.75, 1.0, 1.5})
      for (double end : {-0.5, 0.0, 0.25, 0.5, 0.75, 1.0, 1.5})
      {
        SCOPED_TRACE(::testing::Message() << radius << " " << start << " " << end);
        const gp_Pnt a = curve->Value(start), b = curve->Value(end);
        const double distance = radius * (std::atan2(b.Y(), b.X()) - std::atan2(a.Y(), a.X()));
        EXPECT_NEAR(railLength(rail, start, end), std::abs(distance), radius * 1e-7);
        const auto parameter = railParameterAt(rail, start, distance);
        ASSERT_TRUE(parameter.has_value());
        EXPECT_NEAR(*parameter, end, 1e-6);
      }
  }
}

// A trim can hide a continuity break of the underlying curve. This C1 curve is a line until
// t=1/4, then (t,(t-1/4)^2); its length on either side has a closed form.
TEST(BRepBlendExt, ExtendedRailKeepsKnotsOutsideTheEdgeTrim)
{
  using namespace BRepExt::blending::detail;
  NCollection_Array1<gp_Pnt> poles(1, 4);
  poles(1) = gp_Pnt(0, 0, 0);
  poles(2) = gp_Pnt(0.125, 0, 0);
  poles(3) = gp_Pnt(0.625, 0, 0);
  poles(4) = gp_Pnt(1, 0.5625, 0);
  NCollection_Array1<double> knots(1, 3);
  knots(1) = 0;
  knots(2) = 0.25;
  knots(3) = 1;
  NCollection_Array1<int> mults(1, 3);
  mults(1) = mults(3) = 3;
  mults(2) = 1;
  Handle(Geom_BSplineCurve) curve = new Geom_BSplineCurve(poles, knots, mults, 2);
  const auto primitive = [](double t) {
    const double x = t - 0.25;
    return x <= 0 ? t : 0.25 + x * std::sqrt(1 + 4 * x * x) / 2 + std::asinh(2 * x) / 4;
  };
  for (const auto ends : {std::pair{0.05, 0.2}, std::pair{0.4, 0.6}})
  {
    Rail rail;
    rail.edge = BRepBuilderAPI_MakeEdge(curve, ends.first, ends.second);
    rail.curve.Initialize(rail.edge);
    for (double start : {-0.5, 0.0, 0.2, 0.4, 0.6, 1.0, 1.5})
      for (double end : {-0.5, 0.0, 0.2, 0.4, 0.6, 1.0, 1.5})
      {
        SCOPED_TRACE(::testing::Message() << ends.first << " " << start << " " << end);
        const double distance = primitive(end) - primitive(start);
        EXPECT_NEAR(railLength(rail, start, end), std::abs(distance), 1e-7);
        const auto parameter = railParameterAt(rail, start, distance);
        ASSERT_TRUE(parameter.has_value());
        EXPECT_NEAR(*parameter, end, 1e-6);
      }
  }
}

TEST(BRepBlendExt, RadiusLawCannotHideANegativeSpanBetweenSamples)
{
  for (double scale : {0.001, 1.0, 1000.0})
  {
    const TopoDS_Shape body = BRepPrimAPI_MakeBox(20 * scale, 20 * scale, 20 * scale).Shape();
    const auto edge = edgesWhere(body, [](const double*) { return true; }).front();
    auto request = roundOf(body, {edge}, scale);
    request.radii = {{{0, scale}, {.5001, .01 * scale}, {.5002, .001 * scale},
                       {.5003, scale}, {1, scale}}};
    try
    {
      BRepExt::blending::blend(request);
      ADD_FAILURE() << "a negative radius law reached construction";
    }
    catch (const BRepExt::Failure& failure)
    {
      EXPECT_EQ(failure.kind, BRepExt::FailureKind::RadiusLaw) << failure.what();
    }
  }
}

TEST(BRepBlendExt, RadiusSpanBoundsEncloseStationaryValues)
{
  using BRepExt::blending::detail::SizeLaw;
  // Quadratics with an interior minimum, their negatives, and a monotone cubic
  // with an interior stationary point. Vary span parameter widths independently.
  for (double width : {1e-6, 1.0, 1e6})
    for (double at : {0.1, 0.5, 0.9})
      for (double sign : {-1.0, 1.0})
      {
        SizeLaw law;
        law.places = {2 * width, 3 * width};
        law.values = {sign * (at * at + .01), sign * ((1 - at) * (1 - at) + .01)};
        law.bends = {2 * sign / (width * width), 2 * sign / (width * width)};
        const auto range = law.bounds();
        const double end = std::max(at * at, (1 - at) * (1 - at)) + .01;
        const double low = sign > 0 ? .01 : -end, high = sign > 0 ? end : -.01;
        ASSERT_TRUE(GeomAlgoExt::interval::valid(range));
        EXPECT_LE(range.low, low + 1e-14);
        EXPECT_GE(range.high, high - 1e-14);
        EXPECT_NEAR(range.low, low, 1e-12);
        EXPECT_NEAR(range.high, high, 1e-12);
      }
  for (double at : {0.1, 0.5, 0.9})
  {
    SizeLaw cubic;
    cubic.places = {0, 1};
    cubic.values = {1 - at * at * at, 1 + std::pow(1 - at, 3)};
    cubic.bends = {-6 * at, 6 * (1 - at)};
    const auto range = cubic.bounds();
    EXPECT_NEAR(range.low, cubic.values[0], 1e-12);
    EXPECT_NEAR(range.high, cubic.values[1], 1e-12);
  }
  SizeLaw flat;
  flat.places = {0, 1}; flat.values = {1, 1}; flat.bends = {0, 0};
  const auto constant = flat.bounds();
  EXPECT_GT(constant.low, 0);
  EXPECT_NEAR(constant.low, 1, 1e-14);
  EXPECT_NEAR(constant.high, 1, 1e-14);
}

TEST(BRepBlendExt, CornerStationDistanceUsesPhysicalUnits)
{
  using namespace BRepExt::blending::detail;
  for (double speed : {0.001, 1.0, 1000.0})
    for (double scale : {0.001, 1.0, 1000.0})
    {
      NCollection_Array1<gp_Pnt> poles(1, 2);
      poles(1) = gp_Pnt(0, 0, 0);
      poles(2) = gp_Pnt(10 * scale, 0, 0);
      NCollection_Array1<double> knots(1, 2);
      knots(1) = 0;
      knots(2) = 10 / speed;
      NCollection_Array1<int> multiplicities(1, 2);
      multiplicities.Init(2);
      Handle(Geom_BSplineCurve) curve = new Geom_BSplineCurve(poles, knots, multiplicities, 1);
      Rail rail;
      rail.scale = scale;
      rail.edge = BRepBuilderAPI_MakeEdge(curve);
      rail.curve.Initialize(rail.edge);
      rail.from = 0;
      rail.to = 10 / speed;
      for (double delta : {1e-7, 1e-5, 0.1})
        EXPECT_NEAR(railStationSeparation(rail, 0.5, 0.5 + delta), 10 * scale * delta,
                    1e-10 * scale);
    }
}

TEST(BRepBlendExt, LocatedFitEventsKeepTheirPhysicalStation)
{
  using namespace BRepExt::blending::detail;
  for (double offset : {-9e-7, -4e-7, 4e-7, 9e-7})
    for (int existing : {1, 3, 5})
    {
      std::map<double, int> knots = {{0, 6}, {0.5, existing}, {1, 6}};
      const double event = 0.5 + offset;
      EXPECT_TRUE(raiseFitBreaks(knots, {{event, 5}}));
      ASSERT_TRUE(knots.count(event));
      EXPECT_EQ(knots.at(event), 5);
      if (existing > 1)
        EXPECT_EQ(knots.at(0.5), existing);
      EXPECT_FALSE(raiseFitBreaks(knots, {{event, 5}}));
    }
}

// Intent survives independently of continuity metadata, for single stopped rails
// and sharp mitred loops, across scale and section shape.
TEST(BRepBlendExt, ConstructionRolesDistinguishContactsFromSharpEnds)
{
  using namespace BRepExt::blending;
  for (double scale : {0.01, 1.0, 100.0})
    for (double rho : {0.3, 0.7})
      for (bool loop : {false, true})
      {
        const TopoDS_Shape box = BRepPrimAPI_MakeBox(20*scale, 20*scale, 20*scale).Shape();
        const auto picked = edgesWhere(box, [&](const double* b) {
          return near(b[2], 20*scale) && near(b[5], 20*scale)
                 && (loop || (near(b[0], 0) && near(b[3], 0)));
        });
        Request request = roundOf(box, picked, scale);
        request.rho = rho;
        const Blended made = blend(request);
        ASSERT_FALSE(made.joins.empty());
        int contacts = 0, sharp = 0;
        for (const auto& join : made.joins)
        {
          EXPECT_EQ(join.edge.ShapeType(), TopAbs_EDGE);
          bool retained = false;
          for (TopExp_Explorer e(made.shape, TopAbs_EDGE); e.More(); e.Next())
            retained = retained || e.Current().IsSame(join.edge);
          EXPECT_TRUE(retained);
          if (join.role == JoinRole::Contact)
          {
            EXPECT_EQ(join.order, 1);
            ++contacts;
          }
          if (join.role == (loop ? JoinRole::Miter : JoinRole::Stop))
          {
            EXPECT_EQ(join.order, 0);
            ++sharp;
          }
        }
        EXPECT_GE(contacts, loop ? 8 : 2);
        EXPECT_GE(sharp, loop ? 4 : 2);
      }
}

// Small intervals of a regular offset intersection can need many local refinements without
// needing a large global net. Generate the supports, rather than retaining a serialized body;
// repeat at nearby crossing angles and at different scales and placements.
TEST(BRepBlendExt, LocalFitRefinementOfATorusAndASplineSkin)
{
  const auto oval = [](double z, double scale) {
    constexpr double q[4][2] = {{1, 0}, {1, 0.5522847498307936},
                                {0.5522847498307936, 1}, {0, 1}};
    BRepBuilderAPI_MakeWire wire;
    for (int quadrant = 0; quadrant < 4; ++quadrant)
    {
      NCollection_Array1<gp_Pnt> poles(1, 4);
      const double a = quadrant * M_PI / 2;
      for (int i = 0; i < 4; ++i)
        poles(i + 1) = gp_Pnt(scale * 10 * (q[i][0] * std::cos(a) - q[i][1] * std::sin(a)),
                              scale * 7 * (q[i][0] * std::sin(a) + q[i][1] * std::cos(a)), z);
      wire.Add(BRepBuilderAPI_MakeEdge(new Geom_BezierCurve(poles)).Edge());
    }
    return wire.Wire();
  };
  std::vector<uint32_t> sections;
  for (const auto zr : {std::pair{-10., .8}, std::pair{0., 1.1}, std::pair{10., .9}})
    sections.push_back(BRepExt::store(session(), oval(zr.first, zr.second)));
  const auto skin = BRepExt::get(session(),
    BRepExt::skinThrough(session(), sections, 0, 0, true, false).resultId);
  for (const auto angleRadius : {std::pair{.69, .55}, std::pair{.7, .5}})
  {
    gp_Trsf turn, move;
    turn.SetRotation(gp_Ax1(gp::Origin(), gp_Dir(1, .3, 0)), angleRadius.first);
    move.SetTranslation(gp_Vec(6, 3, 2));
    const auto support = BRepBuilderAPI_Transform(skin, move * turn, true).Shape();
    const auto united = BRepExt::booleanChecked(session(), 0,
      {BRepExt::store(session(), BRepPrimAPI_MakeTorus(9, 4).Shape())},
      {BRepExt::store(session(), support)}, 0, 0, 0, {}, BRepExt::failureHashBound);
    const auto whole = BRepExt::get(session(), united.resultId);
    TopExp_Explorer solids(whole, TopAbs_SOLID);
    ASSERT_TRUE(solids.More());
    const auto body = solids.Current();
    TopoDS_Edge seed;
    double longest = 0;
    for (const auto& edge : edgesWhere(body, [](const double*) { return true; }))
    {
      bool torus = false, spline = false;
      for (TopExp_Explorer faces(body, TopAbs_FACE); faces.More(); faces.Next())
        for (TopExp_Explorer boundary(faces.Current(), TopAbs_EDGE); boundary.More(); boundary.Next())
          if (boundary.Current().IsSame(edge))
          {
            const auto type = BRepAdaptor_Surface(TopoDS::Face(faces.Current())).GetType();
            torus = torus || type == GeomAbs_Torus;
            spline = spline || type == GeomAbs_BSplineSurface;
          }
      const double size = length(edge);
      if (torus && spline && size > longest)
      {
        seed = edge;
        longest = size;
      }
    }
    ASSERT_FALSE(seed.IsNull());
    for (double scale : {.5, 1., 2.})
      for (bool turned : {false, true})
      {
        SCOPED_TRACE(::testing::Message() << angleRadius.first << " " << scale << " " << turned);
        gp_Trsf pose;
        if (turned)
          pose = posed(scale);
        else
          pose.SetScale(gp::Origin(), scale);
        BRepBuilderAPI_Transform transform(body, pose, true);
        const auto input = transform.Shape();
        const auto edge = TopoDS::Edge(transform.ModifiedShape(seed));
        const auto request = roundOf(input, {edge}, angleRadius.second * scale);
        BRepExt::blending::Blended result;
        try
        {
          result = BRepExt::blending::blend(request);
        }
        catch (const std::exception& error)
        {
          ADD_FAILURE() << error.what();
          continue;
        }
        EXPECT_NO_THROW(BRepExt::blending::accepted("round", input, result.shape,
                                                   request.first, &result));
        EXPECT_TRUE(BRepExt::invalidity(result.shape).empty());
        if (!turned && scale == 1)
        {
          const auto defect = ballDefect(result.shape, input, request.first);
          ASSERT_GT(defect.samples, 0);
          EXPECT_LE(std::abs(defect.worst), 1e-6);
        }
      }
  }
}

// A deferred construction must preserve a successful bounded construction byte for byte.
TEST(BRepBlendExt, DeferredFitPreservesBoundedDeliveries)
{
  using namespace BRepExt::blending;
  const auto text = [](const TopoDS_Shape& shape) {
    std::ostringstream out;
    BRepTools::Write(shape, out);
    return out.str();
  };
  for (double scale : {.5, 1., 2.})
    for (bool turned : {false, true})
      for (Law law : {Law::Conic, Law::Chamfer})
      {
        SCOPED_TRACE(::testing::Message() << scale << " " << turned << " " << int(law));
        gp_Trsf pose;
        if (turned)
          pose = posed(scale);
        else
          pose.SetScale(gp::Origin(), scale);
        const auto body = BRepBuilderAPI_Transform(BRepPrimAPI_MakeBox(10, 12, 14).Shape(),
                                                   pose, true).Shape();
        const auto edges = edgesWhere(body, [](const double*) { return true; });
        auto request = roundOf(body, {edges.front()}, scale);
        request.law = law;
        request.fitWork = FitWork::Bounded;
        const auto bounded = blend(request);
        request.fitWork = FitWork::Deferred;
        const auto deferred = blend(request);
        EXPECT_EQ(text(bounded.shape), text(deferred.shape));
        EXPECT_NO_THROW(accepted("round", body, deferred.shape, scale, &deferred));
      }
}

TEST(BRepBlendExt, ContactOwnershipHasACommonPartitionOnANarrowSupport)
{
  using namespace BRepExt::blending;
  const auto text = [](const TopoDS_Shape& shape) {
    std::ostringstream out;
    BRepTools::Write(shape, out);
    return out.str();
  };
  for (double scale : {.5, 1., 2.})
    for (bool turned : {false, true})
      for (double factor : {.98, 1., 1.02})
      {
        SCOPED_TRACE(::testing::Message() << scale << " " << turned << " " << factor);
        const auto box = BRepPrimAPI_MakeBox(gp_Pnt(-10, -10, -10), 20, 20, 20).Shape();
        const auto cone = BRepPrimAPI_MakeCone(gp_Ax2(gp_Pnt(6, 3, -8), gp::DZ()), 12, 6, 20).Shape();
        BRepAlgoAPI_Fuse fused(box, cone);
        const auto body = fused.Shape();
        NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher> edges;
        TopExp::MapShapesAndAncestors(body, TopAbs_EDGE, TopAbs_FACE, edges);
        TopoDS_Edge picked;
        double longest = -1;
        for (int i = 1; i <= edges.Extent(); ++i)
        {
          const auto& faces = edges.FindFromIndex(i);
          if (faces.Extent() != 2)
            continue;
          const auto a = BRepAdaptor_Surface(TopoDS::Face(faces.First())).GetType();
          const auto b = BRepAdaptor_Surface(TopoDS::Face(faces.Last())).GetType();
          if (!((a == GeomAbs_Plane && b == GeomAbs_Cone) || (a == GeomAbs_Cone && b == GeomAbs_Plane)))
            continue;
          const auto plane = TopoDS::Face(a == GeomAbs_Plane ? faces.First() : faces.Last());
          if (BRepAdaptor_Surface(plane).Plane().Distance(gp_Pnt(0, 0, 10)) > Precision::Confusion())
            continue;
          GProp_GProps length;
          BRepGProp::LinearProperties(edges.FindKey(i), length);
          if (length.Mass() > longest)
          {
            longest = length.Mass();
            picked = TopoDS::Edge(edges.FindKey(i));
          }
        }
        ASSERT_FALSE(picked.IsNull());
        gp_Trsf scaling;
        scaling.SetScale(gp::Origin(), scale);
        gp_Trsf pose;
        if (turned)
        {
          pose.SetRotation(gp_Ax1(gp::Origin(), gp_Dir(1, 2, 3)), .43);
          pose.SetTranslationPart(gp_Vec(3*scale, -2*scale, 5*scale));
        }
        BRepBuilderAPI_Transform moved(body, pose*scaling, true);
        const auto input = moved.Shape();
        const auto before = text(input);
        auto request = roundOf(input, {TopoDS::Edge(moved.ModifiedShape(picked))}, scale*factor);
        const auto result = blend(request);
        EXPECT_NO_THROW(accepted("round", input, result.shape, request.first, &result));
        EXPECT_EQ(before, text(input));
      }
}

// The oracle is the independent analytic intersection, not the carrier returned by the SSI
// implementation. A full circuit of contacts must not retrace at a tight bend.
TEST(BRepBlendExt, ClosedOffsetIntersectionHasOneCircuitOfContacts)
{
  using namespace BRepExt::blending::detail;
  for (double scale : {.01, 1., 100.})
    for (bool turned : {false, true})
      for (double degrees : {75., 80.})
        for (double ratio : {.25, .5})
        {
          SCOPED_TRACE(::testing::Message() << scale << " " << turned << " " << degrees
                                            << " " << ratio);
          const double radius = 2 * scale, r = radius * ratio;
          const double angle = degrees * M_PI / 180;
          gp_Cylinder cylinder(gp_Ax3(gp::Origin(), gp::DZ()), radius);
          gp_Pln plane(gp::Origin(), gp_Dir(0, std::sin(angle), std::cos(angle)));
          if (turned)
          {
            const gp_Trsf pose = posed(1);
            cylinder.Transform(pose);
            plane.Transform(pose);
          }
          IntAna_QuadQuadGeo source(plane, cylinder, 1e-12, 1e-10 * scale);
          ASSERT_TRUE(source.IsDone());
          ASSERT_EQ(source.TypeInter(), IntAna_Ellipse);
          const Handle(Geom_Curve) curve = new Geom_Ellipse(source.Ellipse(1));
          Rail rail;
          rail.edge = BRepBuilderAPI_MakeEdge(curve);
          rail.curve.Initialize(rail.edge);
          rail.first = rail.from = rail.curve.FirstParameter();
          rail.last = rail.to = rail.curve.LastParameter();
          rail.closed = true;
          rail.scale = r;
          rail.magnitude = radius / std::cos(angle) + 10;
          rail.s[0].surface = new Geom_CylindricalSurface(cylinder);
          rail.s[1].surface = new Geom_Plane(plane);
          rail.s[0].outward = -1;
          for (int side = 0; side < 2; ++side)
          {
            rail.s[side].distance = r;
            rail.s[side].face = BRepBuilderAPI_MakeFace(rail.s[side].surface, 1e-10 * scale);
            const auto pc = GeomProjLib::Curve2d(curve, rail.first, rail.last, rail.s[side].surface);
            ASSERT_FALSE(pc.IsNull());
            BRep_Builder().UpdateEdge(rail.edge, pc, rail.s[side].face, 1e-9 * scale);
          }
          ASSERT_TRUE(intersectSpine(rail));
          // Restricting a piece at a support seam must not stretch another complete
          // carrier circuit across that piece. Compare its physical station frame.
          Rail restricted = rail;
          restricted.first = rail.t(.2);
          restricted.last = rail.t(.7);
          restricted.from = restricted.first;
          restricted.to = restricted.last;
          for (double u : {.21, .4, .69})
          {
            gp_Pnt wholePoint, partPoint;
            gp_Vec wholeTangent, partTangent;
            rail.stationD1(rail.t(u), wholePoint, wholeTangent);
            restricted.stationD1(rail.t(u), partPoint, partTangent);
            EXPECT_LE(wholePoint.Distance(partPoint), 1e-12 * scale);
            EXPECT_LE((wholeTangent-partTangent).Magnitude(), 1e-12 * scale);
          }
          Stations stations(rail);
          ASSERT_NO_THROW(stations.march());
          double first = 0, previous = 0, minimum = 1e100;
          Station start, end;
          for (int i = 0; i <= 256; ++i)
          {
            Station station;
            ASSERT_TRUE(stations.at(rail.t(i / 256.), station));
            // Analytic offset equations: radial distance R+r and signed plane distance -r.
            EXPECT_NEAR(gp_Lin(cylinder.Axis()).Distance(station.centre), radius + r, 1e-9 * scale);
            const gp_Vec arm(plane.Location(), station.centre);
            EXPECT_NEAR(arm.Dot(gp_Vec(plane.Axis().Direction())), -r, 1e-9 * scale);
            double phase = station.uv[0][0];
            if (i == 0)
            {
              first = phase;
              start = station;
            }
            else
            {
              phase = BRepExt::nearPeriod(phase, previous, 2 * M_PI);
              minimum = std::min(minimum, phase - previous);
            }
            previous = phase;
            end = station;
          }
          // Source orientation can be either way; all increments must have that same sign.
          EXPECT_NEAR(std::abs(previous - first), 2 * M_PI, 1e-9);
          EXPECT_GT(minimum, 0);
          EXPECT_LE(start.centre.Distance(end.centre), 1e-9 * scale);
        }
}

TEST(BRepBlendExt, CentreLengthIgnoresTangentialParameterBreaks)
{
  using namespace BRepExt::blending::detail;
  // x(t)=t+max(0,t-1/2)^2/5 is C1, with a nonzero tangential
  // acceleration jump. Its geometric trace is a straight segment.
  for (double scale : {0.001, 1.0, 1000.0})
    for (bool turned : {false, true})
    {
      gp_Trsf pose;
      if (turned)
        pose.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 2, 3)), 0.7);
      const gp_Pnt               origin(2 * scale, -3 * scale, scale);
      const gp_Vec               direction = gp_Vec(scale, 0, 0).Transformed(pose);
      NCollection_Array1<gp_Pnt> poles(1, 4);
      const double               x[] = {0, 0.25, 0.75, 1.05};
      for (int i = 1; i <= 4; ++i)
        poles(i) = origin.Translated(direction * x[i - 1]);
      NCollection_Array1<double> knots(1, 3);
      NCollection_Array1<int>    mult(1, 3);
      knots(1)       = 0;
      knots(2)       = 0.5;
      knots(3)       = 1;
      mult(1)        = 3;
      mult(2)        = 1;
      mult(3)        = 3;
      auto arc       = std::make_shared<CentreParameter>();
      arc->curve     = new Geom_BSplineCurve(poles, knots, mult, 2);
      arc->knots     = {0, 0.5, 1};
      arc->lengths   = {0, 0.5 * scale, 1.05 * scale};
      arc->tolerance = scale * 1e-12;
      Rail rail;
      rail.from            = 0.1;
      rail.to              = 0.9;
      rail.centreParameter = arc;
      const auto xAt       = [](double t) { return t + std::pow(std::max(0.0, t - 0.5), 2) / 5; };
      for (int i = 0; i <= 32; ++i)
      {
        const double u = i / 32.0, t = rail.t(u);
        const double expected = xAt(rail.from) + u * (xAt(rail.to) - xAt(rail.from));
        EXPECT_NEAR(xAt(t), expected, 1e-10);
        EXPECT_NEAR(rail.u(t), u, 1e-10);
        EXPECT_NEAR(arc->length(t), scale * xAt(t), scale * 1e-10);
      }
    }
}

TEST(BRepBlendExt, CompactRoundingPreservesBoundsAndEndpointJets)
{
  using namespace BRepExt::blending;
  for (double scale : {0.001, 1.0, 1000.0})
  {
    for (int i = -200; i <= 200; ++i)
    {
      const double a = scale * i / 100;
      const double upper = smoothMaximum(a, 0, scale);
      const double lower = smoothMinimum(a, 0, scale);
      EXPECT_GE(upper + 1e-14 * scale, std::max(a, 0.0));
      EXPECT_LE(upper - std::max(a, 0.0), 3 * scale / 16 + 1e-14 * scale);
      EXPECT_LE(lower - 1e-14 * scale, std::min(a, 0.0));
      EXPECT_DOUBLE_EQ(upper, smoothMaximum(0, a, scale));
      if (std::abs(a) >= scale)
        EXPECT_DOUBLE_EQ(upper, std::max(a, 0.0));
    }
    const auto f = [&](double x) { return smoothMaximum(scale * x, 0, scale) / scale; };
    for (double x : {-1.0, 1.0})
    {
      const double h = 1e-4;
      const double left = (3 * f(x) - 4 * f(x-h) + f(x-2*h)) / (2*h);
      const double right = (-3 * f(x) + 4 * f(x+h) - f(x+2*h)) / (2*h);
      EXPECT_NEAR(left, right, 3e-8);
      const double ll = (2*f(x) - 5*f(x-h) + 4*f(x-2*h) - f(x-3*h)) / (h*h);
      const double rr = (2*f(x) - 5*f(x+h) + 4*f(x+2*h) - f(x+3*h)) / (h*h);
      EXPECT_NEAR(ll, rr, 3e-7);
    }
  }
  for (double x : {0.0, 1.0})
  {
    const double h = 1e-4;
    const auto f = smoothActivation;
    const double left = (3*f(x) - 4*f(x-h) + f(x-2*h)) / (2*h);
    const double right = (-3*f(x) + 4*f(x+h) - f(x+2*h)) / (2*h);
    EXPECT_NEAR(left, right, 3e-7);
    const double ll = (2*f(x) - 5*f(x-h) + 4*f(x-2*h) - f(x-3*h)) / (h*h);
    const double rr = (2*f(x) - 5*f(x+h) + 4*f(x+2*h) - f(x+3*h)) / (h*h);
    EXPECT_NEAR(ll, rr, 5e-6);
  }
}

TEST(BRepBlendExt, ReenteredSupportRetainsTheIntervalBetweenCrossings)
{
  using namespace BRepExt::blending;
  BRepBuilderAPI_MakeWire boundary;
  const double q[4][2] = {{1,0},{1,.5522847498307936},{.5522847498307936,1},{0,1}};
  for (int quadrant = 0; quadrant < 4; ++quadrant)
  {
    NCollection_Array1<gp_Pnt> poles(1,4);
    const double a = quadrant * M_PI / 2;
    for (int j = 0; j < 4; ++j)
      poles(j+1) = gp_Pnt(10*(q[j][0]*std::cos(a)-q[j][1]*std::sin(a)),
                           7*(q[j][0]*std::sin(a)+q[j][1]*std::cos(a)), -10);
    boundary.Add(BRepBuilderAPI_MakeEdge(new Geom_BezierCurve(poles)).Edge());
  }
  const auto prism = BRepPrimAPI_MakePrism(BRepBuilderAPI_MakeFace(boundary.Wire()),gp_Vec(0,0,20)).Shape();
  gp_Trsf turn;
  turn.SetRotation(gp_Ax1(gp::Origin(),gp_Dir(1,.3,0)),.7);
  turn.SetTranslationPart(gp_Vec(6,3,2));
  const auto tool = BRepBuilderAPI_Transform(prism,turn,true).Shape();
  BRepAlgoAPI_Common common(BRepPrimAPI_MakeTorus(9,4).Shape(),tool);
  const auto body = unified(common.Shape());
  NCollection_IndexedDataMap<TopoDS_Shape,NCollection_List<TopoDS_Shape>,TopTools_ShapeMapHasher> adjacent;
  TopExp::MapShapesAndUniqueAncestors(body,TopAbs_EDGE,TopAbs_FACE,adjacent);
  std::vector<TopoDS_Edge> junction;
  TopoDS_Edge upper;
  double highest = -HUGE_VAL;
  for (int i=1;i<=adjacent.Extent();++i)
  {
    const auto& faces = adjacent.FindFromIndex(i);
    if (faces.Extent()!=2) continue;
    const auto a=BRepAdaptor_Surface(TopoDS::Face(faces.First())).GetType();
    const auto b=BRepAdaptor_Surface(TopoDS::Face(faces.Last())).GetType();
    if (!((a==GeomAbs_Torus && b==GeomAbs_SurfaceOfExtrusion)
          ||(b==GeomAbs_Torus && a==GeomAbs_SurfaceOfExtrusion))) continue;
    const auto edge=TopoDS::Edge(adjacent.FindKey(i));
    junction.push_back(edge);
    BRepAdaptor_Curve curve(edge);
    const double z=curve.Value((curve.FirstParameter()+curve.LastParameter())/2).Z();
    if (z>highest) { highest=z;upper=edge; }
  }
  ASSERT_FALSE(upper.IsNull());
  std::vector<TopoDS_Edge> chain={upper};
  for(bool grew=true;grew;)
  {
    grew=false;
    for(const auto& candidate:junction)
    {
      if(std::any_of(chain.begin(),chain.end(),[&](const auto& e){return e.IsSame(candidate);}))continue;
      if(std::any_of(chain.begin(),chain.end(),[&](const auto& e){
        return TopExp::FirstVertex(e).IsSame(TopExp::FirstVertex(candidate))
          ||TopExp::FirstVertex(e).IsSame(TopExp::LastVertex(candidate))
          ||TopExp::LastVertex(e).IsSame(TopExp::FirstVertex(candidate))
          ||TopExp::LastVertex(e).IsSame(TopExp::LastVertex(candidate));}))
      {chain.push_back(candidate);grew=true;}
    }
  }
  for(double radius:{.2,.5}) for(double scale:{.5,2.}) for(bool turned:{false,true}) for(bool all:{false,true})
  {
    SCOPED_TRACE(::testing::Message()<<radius<<" "<<scale<<" "<<turned<<" "<<all);
    gp_Trsf transform;
    transform.SetScale(gp::Origin(),scale);
    if(turned)transform=posed(scale);
    BRepBuilderAPI_Transform moved(body,transform,true);
    std::vector<TopoDS_Edge> selected;
    for(const auto& edge:all?chain:std::vector<TopoDS_Edge>{upper})
      selected.push_back(TopoDS::Edge(moved.ModifiedShape(edge)));
    auto request=roundOf(moved.Shape(),selected,radius*scale);
    request.fitWork=FitWork::Deferred;
    const auto made=blend(request);
    EXPECT_NO_THROW(accepted("round",request.body,made.shape,request.first,&made));
    // The ball is cut by a stopping plane near the ends. Away from those planes,
    // compare its envelope to distances from the input skin, independently of the solver.
    const auto rolls = [&](const gp_Pnt& centre) {
      for (TopExp_Explorer f(request.body, TopAbs_FACE); f.More(); f.Next())
      {
        const BRepAdaptor_Surface support(TopoDS::Face(f.Current()));
        if (support.GetType() == GeomAbs_Plane
            && support.Plane().Distance(centre) <= request.first + 1e-6 * scale)
          return false;
      }
      return true;
    };
    const auto defect = ballDefect(made.shape, request.body, request.first, rolls);
    EXPECT_GT(defect.samples, 0);
    EXPECT_LT(std::abs(defect.worst), 5e-6 * request.first);
  }
}

// The affine corner is a rational quadratic tensor, not a fitted many-pole fill.
// Its representation matters as well as continuity: a general fill can refuse
// curvature data which this exact tensor reproduces.
TEST(BRepBlendExt, PlanarConicCornerRetainsItsExactTensor)
{
  using namespace BRepExt::blending;
  for (double scale : {.2, 1., 5.})
    for (bool turned : {false, true})
    {
      SCOPED_TRACE(::testing::Message() << scale << " " << turned);
      gp_Trsf pose;
      if (turned)
        pose = posed(scale);
      else
        pose.SetScale(gp::Origin(), scale);
      const auto body = BRepBuilderAPI_Transform(BRepPrimAPI_MakeBox(10, 10, 10).Shape(),
                                                 pose, true).Shape();
      const gp_Pnt at = gp::Origin().Transformed(pose);
      std::vector<TopoDS_Edge> selected;
      for (const auto& edge : edgesWhere(body, [](const double*) { return true; }))
      {
        TopoDS_Vertex a, b;
        TopExp::Vertices(edge, a, b);
        if ((!a.IsNull() && BRep_Tool::Pnt(a).Distance(at) < 1e-9 * scale)
            || (!b.IsNull() && BRep_Tool::Pnt(b).Distance(at) < 1e-9 * scale))
          selected.push_back(edge);
      }
      ASSERT_EQ(selected.size(), 3u);
      auto request = roundOf(body, selected, .2 * scale);
      request.rho = .7;
      const auto made = blend(request);
      EXPECT_NO_THROW(accepted("round", body, made.shape, request.first, &made));
      EXPECT_TRUE(BRepExt::invalidity(made.shape).empty());
      int tensors = 0;
      for (TopExp_Explorer it(made.shape, TopAbs_FACE); it.More(); it.Next())
      {
        const auto surface = Handle(Geom_BSplineSurface)::DownCast(
          BRep_Tool::Surface(TopoDS::Face(it.Current())));
        if (!surface.IsNull() && surface->UDegree() == 2 && surface->VDegree() == 2
            && surface->NbUPoles() == 3 && surface->NbVPoles() == 3)
          ++tensors;
      }
      EXPECT_EQ(tensors, 1);
    }
}

// Where more than three conic blends meet between planes the corner is split about a centre into
// one rational patch per crossing of contacts, made directly and G1 everywhere (accepted reads
// every seam), on regular and uneven apexes of four to six faces, posed and scaled, and on a steep
// pentagon whose faces span 12 to 38 degrees at the apex (its middle normals once misread the
// split as folded). The G2 law keeps its refusal there.
TEST(BRepBlendExt, ConicCornerOfMoreThanThreeBlendsIsSplitIntoExactPatches)
{
  using namespace BRepExt::blending;
  const auto pyramid = [](int n, double off, double uneven, double height) {
    BRepBuilderAPI_Sewing      sew(1e-6);
    std::vector<gp_Pnt>        base;
    BRepBuilderAPI_MakePolygon bottom;
    for (int i = 0; i < n; ++i)
    {
      const double angle = 2 * M_PI * i / n + uneven * std::sin(3.1 * i + n);
      base.emplace_back(6 * std::cos(angle), 6 * std::sin(angle), 0);
    }
    const gp_Pnt apex(off, off / 2, height);
    for (int i = n - 1; i >= 0; --i)
      bottom.Add(base[i]);
    bottom.Close();
    sew.Add(BRepBuilderAPI_MakeFace(bottom.Wire()).Face());
    for (int i = 0; i < n; ++i)
      sew.Add(BRepBuilderAPI_MakeFace(
                BRepBuilderAPI_MakePolygon(base[i], base[(i + 1) % n], apex, true).Wire())
                .Face());
    sew.Perform();
    BRepBuilderAPI_MakeSolid solid;
    for (TopExp_Explorer x(sew.SewedShape(), TopAbs_SHELL); x.More(); x.Next())
      solid.Add(TopoDS::Shell(x.Current()));
    return std::make_pair(TopoDS_Shape(solid.Solid()), apex);
  };
  struct Member
  {
    int    n;
    double off, uneven, height;
  };
  std::vector<Member> members;
  for (int n = 4; n <= 6; ++n)
    for (double off : {0., 1.5})
      members.push_back({n, off, 0.3, 9});
  members.push_back({5, 0., 0.37, 14});
  for (const Member& m : members)
    for (double scale : {1., 5.})
      {
        const int  n             = m.n;
        const auto [shape, apex] = pyramid(n, m.off, m.uneven, m.height);
        const gp_Trsf pose       = scale == 1. ? gp_Trsf() : posed(scale);
        const auto    body       = BRepBuilderAPI_Transform(shape, pose, true).Shape();
        const gp_Pnt  at         = apex.Transformed(pose);
        std::vector<TopoDS_Edge> selected;
        for (const auto& edge : edgesWhere(body, [](const double*) { return true; }))
        {
          TopoDS_Vertex a, b;
          TopExp::Vertices(edge, a, b);
          if (BRep_Tool::Pnt(a).Distance(at) < 1e-9 * scale
              || BRep_Tool::Pnt(b).Distance(at) < 1e-9 * scale)
            selected.push_back(edge);
        }
        ASSERT_EQ(int(selected.size()), n);
        for (double rho : {.3, .7})
        {
          SCOPED_TRACE(::testing::Message() << n << " " << m.off << " " << m.uneven << " "
                                            << m.height << " " << scale << " " << rho);
          auto request = roundOf(body, selected, .8 * scale);
          request.rho  = rho;
          const auto made = blend(request);
          EXPECT_NO_THROW(accepted("round", body, made.shape, request.first, &made));
          EXPECT_TRUE(BRepExt::invalidity(made.shape).empty());
          int patches = 0;
          for (TopExp_Explorer it(made.shape, TopAbs_FACE); it.More(); it.Next())
          {
            const auto surface = Handle(Geom_BSplineSurface)::DownCast(
              BRep_Tool::Surface(TopoDS::Face(it.Current())));
            patches += !surface.IsNull() && surface->UDegree() == 6 && surface->VDegree() == 6
                       && surface->NbUPoles() == 7 && surface->NbVPoles() == 7;
          }
          EXPECT_EQ(patches, n);
        }
        auto g2 = roundOf(body, selected, .8 * scale);
        g2.law  = Law::G2;
        EXPECT_THROW(blend(g2), BRepExt::Failure);
      }
}

TEST(BRepBlendExt, ConnectedRationalSheetsShareTwoCentreLengthJets)
{
  using namespace BRepExt::blending::detail;
  for (double scale : {0.001, 1., 1000.})
    for (bool rational : {false, true})
      for (bool turned : {false, true})
      for (bool preserveContact : {false, true})
      {
        SCOPED_TRACE(::testing::Message() << scale << " " << rational << " " << turned << " " << preserveContact);
        gp_Trsf transform;
        transform.SetScale(gp::Origin(), scale);
        if (turned) transform = posed(scale);
        NCollection_Array1<gp_Pnt> centres(1, 2);
        centres(1) = gp_Pnt(0,0,0).Transformed(transform);
        centres(2) = gp_Pnt(3,0,0).Transformed(transform);
        NCollection_Array1<double> knots(1, 2); knots(1)=0; knots(2)=3;
        NCollection_Array1<int> mults(1, 2); mults.Init(2);
        auto parameter = std::make_shared<CentreParameter>();
        parameter->curve = new Geom_BSplineCurve(centres, knots, mults, 1);
        parameter->knots = {0,3}; parameter->lengths = {0,3*scale};
        parameter->tolerance = 1e-12*scale;
        const TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(parameter->curve);
        Rail a, b;
        a.curve.Initialize(edge); b.curve.Initialize(edge);
        a.from=0; a.to=1; b.from=1; b.to=3;
        a.scale=b.scale=scale;
        for (auto* rail : {&a,&b})
        {
          rail->centreParameter=parameter;
          rail->s[0].distance=rail->s[1].distance=.2*scale;
        }
        const auto sheet = [&](bool last) {
          NCollection_Array2<gp_Pnt> poles(1,6,1,6);
          NCollection_Array2<double> weights(1,6,1,6);
          for (int i=1;i<=6;++i) for (int j=1;j<=6;++j)
          {
            const double u=(i-1)/5., v=(j-1)/5.;
            const double z=(!preserveContact && (last ? i==3 : i==4)) ? (last ? -1 : 1)*1e-6*v*(1-v) : 0;
            poles(i,j)=gp_Pnt(last ? 1+2*u : u,v,z).Transformed(transform);
            weights(i,j)=rational ? 1+.1*v*(1-v) : 1;
            if (rational && !preserveContact && (last ? i==2 : i==5)) weights(i,j)*=1.000001;
          }
          NCollection_Array1<double> k(1,2);k(1)=0;k(2)=1;
          NCollection_Array1<int> m(1,2);m.Init(6);
          auto result = Handle(Geom_BSplineSurface)(new Geom_BSplineSurface(poles,weights,k,k,m,m,5,5));
          if (preserveContact)
          {
            result->InsertVKnot(.25, 3, 1e-12, false);
            result->InsertVKnot(.75, 3, 1e-12, false);
            // The held contact ribbons are compatible; the interior second
            // jets differ. Both polynomial and rational interior fields vary.
            for (int j = 4; j <= result->NbVPoles()-3; ++j)
            {
              const int i = last ? 3 : 4;
              const gp_Vec delta = gp_Vec(0, 0, (last ? -1 : 1)*1e-6).Transformed(transform);
              result->SetPole(i, j, result->Pole(i, j).Translated(delta),
                              result->Weight(i, j) * (rational ? 1.000001 : 1));
            }
          }
          return result;
        };
        auto sa=sheet(false), sb=sheet(true);
        const auto oldA=Handle(Geom_BSplineSurface)::DownCast(sa->Copy());
        const auto oldB=Handle(Geom_BSplineSurface)::DownCast(sb->Copy());
        ASSERT_NO_THROW(coupleSplitJets(sa,sb,a,b,preserveContact));
        if (preserveContact)
          for (double u : {0., .03, .1, .5, .9, .97, 1.})
            for (double v : {0., 1.})
              for (const auto& pair : {std::make_pair(sa, oldA), std::make_pair(sb, oldB)})
              {
                const auto now = BRepExt::surfaceJet("test", pair.first, u, v);
                const auto old = BRepExt::surfaceJet("test", pair.second, u, v);
                EXPECT_LT(now.p.Distance(old.p), 1e-11*scale);
                EXPECT_LT((now.u-old.u).Magnitude(), 1e-10*scale);
                EXPECT_LT((now.v-old.v).Magnitude(), 1e-10*scale);
                EXPECT_LT((now.uu-old.uu).Magnitude(), 1e-8*scale);
                EXPECT_LT((now.vv-old.vv).Magnitude(), 1e-8*scale);
                EXPECT_LT((now.uv-old.uv).Magnitude(), 1e-8*scale);
              }
        for (double v : {0.,.13,.37,.5,.83,1.})
        {
          gp_Pnt pa,pb;gp_Vec au,av,auu,avv,auv,bu,bv,buu,bvv,buv;
          sa->D2(1,v,pa,au,av,auu,avv,auv);
          sb->D2(0,v,pb,bu,bv,buu,bvv,buv);
          EXPECT_LT(pa.Distance(pb),1e-10*scale);
          EXPECT_LT((au-bu*.5).Magnitude(),1e-9*scale);
          EXPECT_LT((auu-buu*.25).Magnitude(),1e-8*scale);
          EXPECT_LT((av-bv).Magnitude(),1e-9*scale);
          EXPECT_LT((auv-buv*.5).Magnitude(),1e-8*scale);
          for(double u:{0.,.2,.5})
            EXPECT_LT(sa->Value(u,v).Distance(oldA->Value(u,v)),1e-11*scale);
          for(double u:{.5,.8,1.})
            EXPECT_LT(sb->Value(u,v).Distance(oldB->Value(u,v)),1e-11*scale);
        }
      }
}

TEST(BRepBlendExt, CentreLengthCarrierRetainsPositionalSupportJoins)
{
  using namespace BRepExt::blending::detail;
  // Two straight centre branches share a position, with independent tangents.
  // The exact length is piecewise linear; a smooth interpolant cannot represent
  // the corner. This checks the carrier independently of section fitting.
  constexpr double join = .37;
  for (double scale : {.001, 1., 1000.})
    for (bool turned : {false, true})
    {
      SCOPED_TRACE(::testing::Message() << scale << " " << turned);
      gp_Trsf transform;
      transform.SetScale(gp::Origin(), scale);
      if (turned)
        transform = posed(scale);
      Rail rail;
      rail.curve.Initialize(BRepBuilderAPI_MakeEdge(gp_Pnt(0, 0, 0), gp_Pnt(1, 0, 0)));
      rail.first = rail.from = 0;
      rail.last = rail.to = 1;
      rail.scale = scale;
      const auto point = [&](double t) {
        return gp_Pnt(t, t <= join ? 0 : 2 * (t - join), 0).Transformed(transform);
      };
      const auto carrier = centreParameter(rail, point, {join});
      ASSERT_FALSE(carrier->curve.IsNull());
      for (double t : {0., .1, .3, join, .4, .6, .9, 1.})
      {
        const double length = scale * (t <= join ? t : join + std::sqrt(5.) * (t - join));
        EXPECT_LT(carrier->curve->Value(t).Distance(point(t)), 1e-10 * scale);
        EXPECT_NEAR(carrier->length(t), length, 1e-10 * scale);
        EXPECT_NEAR(carrier->parameter(length), t, 1e-9);
      }
    }
}


TEST(BRepBlendExt, CentreLengthRetainsTheDeliveringCircuit)
{
  using namespace BRepExt::blending::detail;
  for (double scale : {.001, 1., 1000.})
    for (bool reversed : {false, true})
      for (bool turned : {false, true})
      {
        SCOPED_TRACE(::testing::Message() << scale << " " << reversed << " " << turned);
        Rail rail;
        rail.curve.Initialize(BRepBuilderAPI_MakeEdge(gp_Pnt(0,0,0),gp_Pnt(1,0,0)));
        rail.first=rail.from=0;
        rail.last=rail.to=1;
        rail.scale=scale;
        const gp_Circ circle(gp_Ax2(gp_Pnt(0,0,0),gp_Dir(0,0,1)),scale);
        rail.spine=GeomConvert::CurveToBSplineCurve(new Geom_Circle(circle));
        if (turned) rail.spine->Transform(posed(1.));
        const double first=rail.spine->FirstParameter(), span=rail.spine->LastParameter()-first;
        rail.spineFirst=first+.37*span;
        rail.spineLast=rail.spineFirst+(reversed ? -span : span);
        rail.spineParameterFirst=0;
        rail.spineParameterLast=1;
        bool remarched=false;
        const auto parameter=centreParameter(rail,[&](double) { remarched=true; return gp_Pnt(); });
        EXPECT_FALSE(remarched);
        EXPECT_NEAR(parameter->length(1),2*M_PI*scale,1e-10*scale);
        for (int i=0;i<=32;++i)
        {
          const double t=i/32.;
          gp_Pnt point; gp_Vec tangent;
          rail.stationD1(t,point,tangent);
          EXPECT_LT(parameter->curve->Value(t).Distance(point),1e-10*scale);
          EXPECT_NEAR(parameter->parameter(parameter->length(t)),t,1e-9);
        }
      }
}


TEST(BRepBlendExt, CentreCarrierRetainsSourceNormalEventsBeforeContactEvents)
{
  using namespace BRepExt::blending::detail;
  // Two perpendicular planes, one represented as a finite spline extrusion.
  // The offset contact crosses its profile end at t=sqrt(2)*scale-radius;
  // the source pcurve crosses at t=sqrt(2)*scale. These are distinct data paths.
  for (double scale : {.001, 1., 1000.})
    for (bool turned : {false, true})
    {
      SCOPED_TRACE(::testing::Message() << scale << " " << turned);
      const double radius = .1*scale;
      const gp_Trsf pose = turned ? posed(1) : gp_Trsf();
      NCollection_Array1<gp_Pnt> poles(1,2);
      poles(1)=gp_Pnt(0,0,0); poles(2)=gp_Pnt(0,0,scale);
      Rail rail;
      rail.s[0].surface = new Geom_SurfaceOfLinearExtrusion(new Geom_BezierCurve(poles),gp::DY());
      rail.s[1].surface = new Geom_Plane(gp_Pln(gp::Origin(),gp_Dir(0,1,-1)));
      const Handle(Geom_Curve) curve = new Geom_Line(gp::Origin(),gp_Dir(0,1,1));
      curve->Transform(pose);
      rail.edge = BRepBuilderAPI_MakeEdge(curve,0,2*scale);
      rail.curve.Initialize(rail.edge);
      rail.first=rail.from=0; rail.last=rail.to=2*scale;
      rail.scale=radius; rail.magnitude=10+3*scale;
      rail.byCentre=true;
      for (int side=0;side<2;++side)
      {
        rail.s[side].surface->Transform(pose);
        rail.s[side].distance=radius;
        rail.s[side].face=BRepBuilderAPI_MakeFace(rail.s[side].surface,1e-10*scale);
        const auto pc=GeomProjLib::Curve2d(curve,rail.first,rail.last,rail.s[side].surface);
        ASSERT_FALSE(pc.IsNull());
        BRep_Builder().UpdateEdge(rail.edge,pc,rail.s[side].face,1e-9*scale);
      }
      Stations stations(rail);
      ASSERT_NO_THROW(stations.march());
      const auto sections=breaksOf(rail,stations,true);
      const auto centres=breaksOf(rail,stations,true,true);
      const double source=std::sqrt(.5),contact=(std::sqrt(2.)*scale-radius)/(2*scale);
      const auto has = [](const auto& events,double u) {
        return std::any_of(events.begin(),events.end(),[&](const auto& event) {
          return std::abs(event.first-u)<1e-10 && event.second==kDegree;
        });
      };
      EXPECT_TRUE(has(sections,contact));
      EXPECT_FALSE(has(sections,source));
      EXPECT_TRUE(has(centres,contact));
      EXPECT_TRUE(has(centres,source));
    }
}

TEST(BRepBlendExt, RoundContinuesPastAnEntireSmoothSupport)
{
  using namespace BRepExt::blending;
  const double          q[4][2] = {{1, 0}, {1, .5522847498307936}, {.5522847498307936, 1}, {0, 1}};
  std::vector<uint32_t> sections;
  for (const auto [z, size] : {std::pair{-10., .8}, std::pair{0., 1.1}, std::pair{10., .9}})
  {
    BRepBuilderAPI_MakeWire boundary;
    for (int quadrant = 0; quadrant < 4; ++quadrant)
    {
      NCollection_Array1<gp_Pnt> poles(1, 4);
      const double               a = quadrant * M_PI / 2;
      for (int j = 0; j < 4; ++j)
        poles(j + 1) = gp_Pnt(size * 10 * (q[j][0] * std::cos(a) - q[j][1] * std::sin(a)),
                              size * 7 * (q[j][0] * std::sin(a) + q[j][1] * std::cos(a)),
                              z);
      boundary.Add(BRepBuilderAPI_MakeEdge(new Geom_BezierCurve(poles)).Edge());
    }
    sections.push_back(BRepExt::store(session(), boundary.Wire()));
  }
  const auto skin =
    BRepExt::get(session(), BRepExt::skinThrough(session(), sections, 0, 0, true, false).resultId);
  gp_Trsf translate;
  translate.SetTranslation(gp_Vec(6, 3, 2));
  const auto tool = BRepBuilderAPI_Transform(skin, translate, true).Shape();
  const auto body = unified(BRepAlgoAPI_Common(BRepPrimAPI_MakeTorus(9, 4).Shape(), tool).Shape());
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    adjacent;
  TopExp::MapShapesAndUniqueAncestors(body, TopAbs_EDGE, TopAbs_FACE, adjacent);
  TopoDS_Edge seed;
  for (int i = 1; i <= adjacent.Extent(); ++i)
  {
    const auto& faces = adjacent.FindFromIndex(i);
    if (faces.Extent() != 2)
      continue;
    const auto a = BRepAdaptor_Surface(TopoDS::Face(faces.First())).GetType();
    const auto b = BRepAdaptor_Surface(TopoDS::Face(faces.Last())).GetType();
    if (!((a == GeomAbs_Torus && b == GeomAbs_BSplineSurface)
          || (b == GeomAbs_Torus && a == GeomAbs_BSplineSurface)))
      continue;
    seed = TopoDS::Edge(adjacent.FindKey(i));
    break; // Selecting one edge propagates around its entire tangent chain.
  }
  ASSERT_FALSE(seed.IsNull());
  for (double radius : {.4, .5, .6})
    for (double scale : {.5, 2.})
      for (bool turned : {false, true})
      {
        SCOPED_TRACE(::testing::Message() << radius << " " << scale << " " << turned);
        gp_Trsf transform;
        transform.SetScale(gp::Origin(), scale);
        if (turned)
          transform = posed(scale);
        BRepBuilderAPI_Transform moved(body, transform, true);
        auto                     request =
          roundOf(moved.Shape(), {TopoDS::Edge(moved.ModifiedShape(seed))}, radius * scale);
        request.fitWork = FitWork::Deferred;
        const auto made = blend(request);
        EXPECT_NO_THROW(accepted("round", request.body, made.shape, request.first, &made));
        const auto defect = ballDefect(made.shape, request.body, request.first);
        EXPECT_GT(defect.samples, 0);
        EXPECT_LT(std::abs(defect.worst), 5e-6 * request.first);
      }
}

TEST(BRepBlendExt, RetainedContactsRespectTrimmedDomains)
{
  using namespace BRepExt::blending::detail;
  // The hole is away from the endpoints and midpoint of the contact. Sampling
  // those three places alone would incorrectly accept the whole interval.
  for (double scale : {0.001, 1.0, 1000.0})
  {
    BRepBuilderAPI_MakeFace face(gp_Pln(gp::XOY()), 0, scale, 0, scale);
    BRepBuilderAPI_MakePolygon hole;
    for (const auto [x, y] : {std::pair{.23,.2}, std::pair{.23,.8},
                             std::pair{.24,.8}, std::pair{.24,.2}})
      hole.Add(gp_Pnt(x * scale, y * scale, 0));
    hole.Close();
    face.Add(hole.Wire());
    const Handle(Geom2d_Curve) contact =
      new Geom2d_Line(gp_Pnt2d(0, .5 * scale), gp_Dir2d(1, 0));
    EXPECT_FALSE(contactWithinFace(face.Face(), contact, .1 * scale, .9 * scale));
    EXPECT_FALSE(contactWithinFace(face.Face(), contact, .9 * scale, .1 * scale));
    EXPECT_TRUE(contactWithinFace(face.Face(), contact, .1 * scale, .2 * scale));
    EXPECT_TRUE(contactWithinFace(face.Face(), contact, .3 * scale, .9 * scale));
    const Handle(Geom2d_Curve) onBoundary =
      new Geom2d_Line(gp_Pnt2d(0, .2 * scale), gp_Dir2d(1, 0));
    EXPECT_TRUE(contactWithinFace(face.Face(), onBoundary, .1 * scale, .9 * scale));
    EXPECT_FALSE(contactWithinFace(face.Face(), contact, -.1 * scale, .2 * scale));
  }
  // A periodic chart can be displaced by arbitrarily many whole periods.
  const Handle(Geom_Surface) cylinder = new Geom_CylindricalSurface(gp::XOY(), 2);
  const auto face = BRepBuilderAPI_MakeFace(cylinder, .2, 1.2, -1, 1, 1e-7).Face();
  for (int turns : {-7, 0, 5})
  {
    const Handle(Geom2d_Curve) contact =
      new Geom2d_Line(gp_Pnt2d(turns * 2 * M_PI, 0), gp_Dir2d(1, 0));
    EXPECT_TRUE(contactWithinFace(face, contact, .3, 1.1));
    EXPECT_FALSE(contactWithinFace(face, contact, .3, 1.3));
  }
}

TEST(BRepBlendExt, FiniteContactToleranceIsMeasuredOnTheFace)
{
  using namespace BRepExt::blending::detail;
  // Along a sphere's meridian the U speed is R*cos(v), not the maximum R
  // over the face. Near a pole an angular tolerance based on R is too small.
  for (double scale : {.001, 1., 1000.})
    for (double latitude : {-1.495, .3, 1.495})
    {
      const double               radius = 11 * scale;
      const Handle(Geom_Surface) sphere = new Geom_SphericalSurface(gp::XOY(), radius);
      const auto   face = BRepBuilderAPI_MakeFace(sphere, 0, 1, -M_PI / 2, M_PI / 2, 1e-7).Face();
      const double tolerance = BRep_Tool::Tolerance(face);
      for (int turns : {-7, 0, 5})
      {
        const Handle(Geom2d_Curve) contact =
          new Geom2d_Line(gp_Pnt2d(turns * 2 * M_PI, latitude), gp_Dir2d(1, 0));
        for (double share : {.25, .75, 2.})
        {
          const double angle = share * tolerance / (radius * std::cos(latitude));
          // Independent closed form: distance from a point on the sphere to
          // the meridian great circle. The nearest arc is inside the V trim.
          const double x        = std::pow(std::cos(latitude) * std::sin(angle), 2);
          const double distance = radius * std::sqrt(2 * x / (1 + std::sqrt(1 - x)));
          SCOPED_TRACE(::testing::Message() << scale << " " << latitude << " " << turns << " "
                                            << share << " " << distance);
          EXPECT_EQ(contactWithinFace(face, contact, -angle, .2), distance <= tolerance);
          EXPECT_EQ(contactWithinFace(face, contact, .2, -angle), distance <= tolerance);
        }
      }
    }
}

TEST(BRepBlendExt, FiniteContactsUseRetainedCurves)
{
  // Exact sheets replace the fitted endpoint maps. The finite-face check must
  // use the same pcurve and retained interval that attachment stores.
  for (bool plane : {false, true})
    for (bool crossed : {false, true})
      for (double radius : {.5, 1.})
      {
        SCOPED_TRACE(::testing::Message() << plane << " " << crossed << " " << radius);
        const auto   left  = plane ? BRepPrimAPI_MakeBox(gp_Pnt(-10, -10, -10), 20, 20, 20).Shape()
                                   : BRepPrimAPI_MakeSphere(11).Shape();
        TopoDS_Shape right = BRepPrimAPI_MakeSphere(11).Shape();
        if (crossed)
        {
          gp_Trsf turn;
          turn.SetRotation(gp_Ax1(gp::Origin(), gp_Dir(1, .3, 0)), .7);
          right = BRepBuilderAPI_Transform(right, turn, true).Shape();
        }
        gp_Trsf move;
        move.SetTranslation(gp_Vec(6, 3, 2));
        right = BRepBuilderAPI_Transform(right, move, true).Shape();
        const TopoDS_Shape shape =
          plane ? BRepAlgoAPI_Fuse(left, right).Shape() : BRepAlgoAPI_Common(left, right).Shape();
        NCollection_IndexedDataMap<TopoDS_Shape,
                                   NCollection_List<TopoDS_Shape>,
                                   TopTools_ShapeMapHasher>
          adjacent;
        TopExp::MapShapesAndUniqueAncestors(shape, TopAbs_EDGE, TopAbs_FACE, adjacent);
        TopoDS_Edge seed;
        double      longest = 0;
        // Choose the intersection edge for this regression test only.
        for (int i = 1; i <= adjacent.Extent(); ++i)
        {
          const auto& faces = adjacent.FindFromIndex(i);
          if (faces.Size() != 2)
            continue;
          const auto a = BRepAdaptor_Surface(TopoDS::Face(faces.First())).GetType();
          const auto b = BRepAdaptor_Surface(TopoDS::Face(faces.Last())).GetType();
          if (plane ? !((a == GeomAbs_Plane && b == GeomAbs_Sphere)
                        || (b == GeomAbs_Plane && a == GeomAbs_Sphere))
                    : !(a == GeomAbs_Sphere && b == GeomAbs_Sphere))
            continue;
          const auto   edge   = TopoDS::Edge(adjacent.FindKey(i));
          const double length = GCPnts_AbscissaPoint::Length(BRepAdaptor_Curve(edge));
          if (length > longest)
          {
            longest = length;
            seed    = edge;
          }
        }
        ASSERT_FALSE(seed.IsNull());
        auto       request = roundOf(shape, {seed}, radius);
        const auto made    = blend(request);
        EXPECT_NO_THROW(accepted("round", request.body, made.shape, request.first, &made));
        const auto defect = ballDefect(made.shape, request.body, request.first);
        EXPECT_GT(defect.samples, 0);
        EXPECT_LT(std::abs(defect.worst), 5e-6 * radius);
      }
}

TEST(BRepBlendExt, FiniteContactUsesLiftedTrims)
{
  using namespace BRepExt::blending::detail;
  for (double scale : {.1, 1., 10.})
  {
    const auto   face = BRepBuilderAPI_MakeFace(gp_Pln(gp::XOY()), 0, scale, 0, scale).Face();
    BRep_Builder builder;
    // The spatial representation is independent of the face's pcurves and
    // may differ within the edge tolerance. It must not move the face trim.
    for (TopExp_Explorer e(face, TopAbs_EDGE); e.More(); e.Next())
    {
      const auto edge = TopoDS::Edge(e.Current());
      double     first, last;
      const auto curve = Handle(Geom_Curve)::DownCast(BRep_Tool::Curve(edge, first, last)->Copy());
      curve->Translate(gp_Vec(0, 0, .25 * scale));
      builder.UpdateEdge(edge, curve, .3 * scale);
    }
    for (double share : {0., .75, 2.})
    {
      const Handle(Geom2d_Curve) contact =
        new Geom2d_Line(gp_Pnt2d(0, -share * BRep_Tool::Tolerance(face)), gp_Dir2d(1, 0));
      EXPECT_EQ(contactWithinFace(face, contact, .2 * scale, .8 * scale), share <= 1);
    }
    // Classification does not repair the caller's independent edge curves.
    for (TopExp_Explorer e(face, TopAbs_EDGE); e.More(); e.Next())
    {
      double     first, last;
      const auto curve = BRep_Tool::Curve(TopoDS::Edge(e.Current()), first, last);
      EXPECT_DOUBLE_EQ(curve->Value((first + last) / 2).Z(), .25 * scale);
    }
  }
}

// A block along x whose edge on the x axis lies between the plane z = 0 (y < 0) and the ruled
// face through the points (x, W, W s(x)) (y > 0): convex where s < 0, concave where s > 0, the
// faces tangent where s = 0.
static TopoDS_Shape bentBlock(const std::function<double(double)>& s,
                              double                               L,
                              double                               W,
                              double                               D,
                              int                                  sections = 12)
{
  BRepOffsetAPI_ThruSections through(true, false);
  for (int k = 0; k <= sections; ++k)
  {
    const double               x = L * k / sections;
    BRepBuilderAPI_MakePolygon p;
    p.Add(gp_Pnt(x, 0, 0));
    p.Add(gp_Pnt(x, -W, 0));
    p.Add(gp_Pnt(x, -W, -D));
    p.Add(gp_Pnt(x, W, -D));
    p.Add(gp_Pnt(x, W, W * s(x)));
    p.Close();
    through.AddWire(p.Wire());
  }
  through.Build();
  return through.Shape();
}

// Convexity is the faces' turn all along the edge, not at its middle alone. Where it changes sign
// the ball would change sides where the faces are tangent: refused by that name, wherever along
// the edge the change lies (the middle one was taken for an edge between tangent faces and left
// out without a word). Where the faces are tangent at the middle only and turn one way elsewhere,
// the blend narrows to the edge's point there and is the rolling ball's all along.
TEST(BRepBlendExt, ConvexityIsDecidedAlongTheWholeEdge)
{
  const double L = 40, W = 10, D = 10;
  for (const double s : {0.1, 1.0, 10.0})
    for (const double a : {0.5, -0.5})
      for (const double x0 : {0.3, 0.5, 0.7})
      {
        const gp_Trsf      pose = s == 1 ? gp_Trsf() : posed(s);
        const TopoDS_Shape body =
          BRepBuilderAPI_Transform(
            bentBlock([=](double x) { return 2 * a * (x / L - x0); }, L, W, D),
            pose,
            true)
            .Shape();
        const TopoDS_Edge edge = edgeThrough(body, gp_Pnt(L * 0.37, 0, 0).Transformed(pose), s * L);
        ASSERT_FALSE(edge.IsNull());
        try
        {
          BRepExt::blending::blend(roundOf(body, {edge}, s));
          ADD_FAILURE() << "delivered: s=" << s << " a=" << a << " x0=" << x0;
        }
        catch (const std::exception& refused)
        {
          EXPECT_NE(std::string(refused.what()).find("from convex to concave"), std::string::npos)
            << refused.what();
        }
      }
  // A concave dip between two of any 16 places along the edge, convex elsewhere: refused all the
  // same, the sign being proven on every part of the edge rather than read at places.
  {
    const TopoDS_Shape dipped = bentBlock(
      [=](double x) { return -0.05 + 0.2 * std::exp(-std::pow((x / L - 0.53) / 0.015, 2)); },
      L,
      W,
      D,
      48);
    const TopoDS_Edge on = edgeThrough(dipped, gp_Pnt(L * 0.37, 0, 0), L);
    ASSERT_FALSE(on.IsNull());
    try
    {
      BRepExt::blending::blend(roundOf(dipped, {on}, 1));
      ADD_FAILURE() << "the dipped edge delivered";
    }
    catch (const std::exception& refused)
    {
      EXPECT_NE(std::string(refused.what()).find("from convex to concave"), std::string::npos)
        << refused.what();
    }
  }
  // Faces tangent at the middle only (to 4.7e-8 rad, the swept surface's own), turned and
  // scaled: the same blend in every pose.
  for (const double s : {1.0, 0.1, 10.0})
  {
    const gp_Trsf      pose = s == 1 ? gp_Trsf() : posed(s);
    const TopoDS_Shape body =
      BRepBuilderAPI_Transform(
        bentBlock([=](double x) { return -2 * (x / L - 0.5) * (x / L - 0.5); }, L, W, D),
        pose,
        true)
        .Shape();
    const TopoDS_Edge  edge = edgeThrough(body, gp_Pnt(L * 0.37, 0, 0).Transformed(pose), s * L);
    const TopoDS_Shape made = blended(roundOf(body, {edge}, s));
    ASSERT_FALSE(made.IsNull()) << "s=" << s;
    EXPECT_TRUE(BRepExt::invalidity(made).empty());
    EXPECT_NEAR((volume(body) - volume(made)) / (s * s * s), 0.0260952, 1e-6) << "s=" << s;
    EXPECT_LT(std::abs(ballDefect(made, body, s).worst), 1e-6 * s) << "s=" << s;
  }
}

// The graph z=a*(x_+)^3 is C2 at zero. Its two polynomial extensions give
// independently known contact positions/normals on either side of that knot.
TEST(BRepBlendExt, ContinuedSectionsCarryTheirCompleteSupportBranch)
{
  using namespace BRepExt::blending;
  using namespace BRepExt::blending::detail;
  for (double scale : {0.01, 1.0, 100.0})
    for (double amplitude : {-0.1, 0.1, 0.2})
      for (bool turned : {false, true})
      {
        SCOPED_TRACE(::testing::Message() << scale << " " << amplitude << " " << turned);
        NCollection_Array2<gp_Pnt> poles(1, 5, 1, 2);
        const double x[5] = {-1, -2.0 / 3, 0, 2.0 / 3, 1};
        for (int i = 1; i <= 5; ++i)
          for (int j = 1; j <= 2; ++j)
            poles(i, j) = gp_Pnt(scale * x[i - 1], scale * (j == 1 ? -1 : 1),
                                 i == 5 ? amplitude * scale : 0);
        NCollection_Array1<double> uk(1, 3), vk(1, 2);
        NCollection_Array1<int> um(1, 3), vm(1, 2);
        uk(1) = -scale; uk(2) = 0; uk(3) = scale;
        um(1) = 4; um(2) = 1; um(3) = 4;
        vk(1) = -scale; vk(2) = scale; vm(1) = vm(2) = 2;
        Handle(Geom_BSplineSurface) support =
          new Geom_BSplineSurface(poles, uk, vk, um, vm, 3, 1);
        Handle(Geom_Surface) plane = new Geom_Plane(
          gp_Ax3(gp::Origin(), gp_Dir(0, -1, 0), gp_Dir(1, 0, 0)));
        auto source = support->VIso(0);
        gp_Trsf pose;
        if (turned)
          pose = posed(1);
        support->Transform(pose); plane->Transform(pose); source->Transform(pose);
        Rail rail;
        rail.edge = BRepBuilderAPI_MakeEdge(source, -scale / 2, scale / 2);
        rail.curve = BRepAdaptor_Curve(rail.edge);
        rail.first = rail.from = -scale / 2;
        rail.last = rail.to = scale / 2;
        rail.scale = 0.2 * scale;
        rail.s[0].surface = support; rail.s[1].surface = plane;
        rail.s[0].distance = rail.s[1].distance = rail.scale;
        rail.law = Law::Conic;
        for (double fraction : {-0.05, 0.05})
        {
          Station seed;
          seed.t = fraction * scale;
          seed.uv[0][0] = seed.uv[1][0] = seed.t;
          seed.uv[0][1] = rail.scale; seed.uv[1][1] = -rail.scale;
          ASSERT_TRUE(contactAt(rail, seed.t, seed));
          for (int end = 0; end < 2; ++end)
          {
            Section section;
            ASSERT_TRUE(sectionAcrossSupportKnot(rail, seed, 0, 0, 2, end, section));
            const gp_Pnt foot = section.poles.front().Transformed(pose.Inverted());
            const gp_Pnt other = section.poles.back().Transformed(pose.Inverted());
            const double a = end ? amplitude / (scale * scale) : 0;
            EXPECT_NEAR(foot.Z(), a * foot.X() * foot.X() * foot.X(), 1e-8 * scale);
            EXPECT_NEAR(other.Y(), 0, 1e-8 * scale);
            const gp_Pnt centre(other.X(), rail.scale, other.Z());
            EXPECT_NEAR(centre.Distance(foot), rail.scale, 1e-8 * scale);
            const gp_Vec normal = gp_Vec(-3 * a * foot.X() * foot.X(), 0, 1).Normalized();
            EXPECT_LT((gp_Vec(foot, centre) + normal * rail.scale).Magnitude(), 1e-8 * scale);
            gp_Pnt p; gp_Vec tangent;
            sectionCurve(section)->D1(0, p, tangent);
            tangent.Transform(pose.Inverted());
            EXPECT_NEAR(tangent.Normalized().Dot(normal), 0, 1e-8);
          }
          Section mixed;
          ASSERT_TRUE(sectionAcrossSupportKnot(rail, seed, 0, 0, 2, 0.5, mixed));
          for (double weight : mixed.weights)
            EXPECT_GT(weight, 0);
        }
      }
}
