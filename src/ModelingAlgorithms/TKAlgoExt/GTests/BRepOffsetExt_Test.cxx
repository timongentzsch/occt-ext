// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Push-pull, draft, defeature and hollow of BRepOffsetExt on bodies whose result has a closed form
// of its input: a box, a cylinder, a box with a boss.

#include "BRepExtTest_Support.hxx"

#include <BRepOffsetExt_Kept.hxx>
#include <BRepOffset_MakeOffset.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBlendExt_Request.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <gp_Ax2.hxx>
#include <BRepOffsetExt_Parts.hxx>
#include <BRepOffsetExt_Trim.hxx>
#include <Geom_Plane.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom2d_Line.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <Geom_BezierCurve.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <BRepOffsetExt_Focal.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <Geom_BezierSurface.hxx>
#include <GeomConvert.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_Failure.hxx>

#include <BRepAdaptor_Surface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <Geom2dConvert.hxx>
#include <NCollection_Array2.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepFilletAPI_MakeChamfer.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>

#include <gtest/gtest.h>

using namespace BRepExtTest;

namespace
{
const double a = 60, b = 40, c = 30;

uint32_t keep(const TopoDS_Shape& shape)
{
  return BRepExt::store(session(), shape);
}

//! Every face of `body` moved by `distance` with its topology kept (BRepOffsetExt_Kept.hxx), or a
//! null shape and the test failed with the reason.
TopoDS_Shape kept(const TopoDS_Shape& body, double distance, double tolerance)
{
  BRepExt::FaceSet all;
  for (TopExp_Explorer f(body, TopAbs_FACE); f.More(); f.Next())
    all.Add(f.Current());
  BRepExt::Lineage   faces(body);
  std::string        why;
  const TopoDS_Shape made = BRepExt::keptOffset(body, all, distance, tolerance, faces, why);
  if (made.IsNull())
    ADD_FAILURE() << why;
  else
    for (const auto type : {TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX})
    {
      BRepExt::ShapeSet before, after, covered;
      TopExp::MapShapes(body, type, before);
      TopExp::MapShapes(made, type, after);
      for (const auto& source : before)
      {
        const auto& products = faces.of(source);
        EXPECT_FALSE(products.empty());
        for (const auto& product : products)
        {
          EXPECT_TRUE(after.Contains(product));
          covered.Add(product);
        }
      }
      // Separated contacts add endpoints at one original vertex. Every new endpoint
      // must be traced to that vertex, including the original common endpoint.
      if (type == TopAbs_VERTEX)
        EXPECT_EQ(covered.Extent(), after.Extent());
    }
  return made;
}

int faceCount(const TopoDS_Shape& shape)
{
  int count = 0;
  for (TopExp_Explorer f(shape, TopAbs_FACE); f.More(); f.Next())
    ++count;
  return count;
}

std::vector<uint32_t> keep(const std::vector<TopoDS_Face>& faces)
{
  std::vector<uint32_t> ids;
  for (const TopoDS_Face& face : faces)
    ids.push_back(keep(face));
  return ids;
}
} // namespace

// Pulling the top face of a box out by d lengthens the box to c + d and leaves every side face
// planar and square to it: a x b x (c + d). Pushing it in by d shortens it likewise. Neighbours
// extended or trimmed wrongly change the volume by a multiple of d.
TEST(BRepOffsetExt, PushPullOfTheTopFaceMovesItAlongItsNormal)
{
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(a, b, c).Shape();
  const auto top = facesWhere(box, [](const double* l) { return near(l[2], c) && near(l[5], c); });
  ASSERT_EQ(top.size(), 1u);
  for (const double d : {7.0, -7.0})
  {
    const auto result =
      BRepExt::offsetFaces(session(), keep(box), keep(top), d, 1e-6, {}, BRepExt::failureHashBound);
    const TopoDS_Shape& made = BRepExt::get(session(), result.resultId);
    EXPECT_TRUE(BRepExt::invalidity(made).empty());
    EXPECT_NEAR(volume(made), a * b * (c + d), 1e-9 * a * b * c) << "d=" << d;
  }
}

// Drafting the four sides of a box by alpha about the plane z = 0, pulling along +z, tilts each
// side inward by alpha as z grows, so a section at height z is (a - 2 z tan alpha)(b - 2 z tan
// alpha) and the volume is its integral over [0, c]: abc - (a + b) c^2 tan alpha + 4/3 c^3 tan^2
// alpha. A draft about the wrong plane, in the wrong direction or by the wrong angle changes it at
// first order in alpha.
TEST(BRepOffsetExt, DraftOfTheSidesOfABoxIsAFrustumOfClosedFormVolume)
{
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(a, b, c).Shape();
  const auto sides       = facesWhere(box, [](const double* l) { return l[5] - l[2] > c - 1e-6; });
  ASSERT_EQ(sides.size(), 4u);
  for (const double degrees : {2.0, 5.0})
  {
    const double        alpha = degrees * M_PI / 180, t = std::tan(alpha);
    const auto          result = BRepExt::draftFaces(session(),
                                            keep(box),
                                            keep(sides),
                                            alpha,
                                                     {0, 0, 1},
                                                     {0, 0, 0},
                                                     {0, 0, 1},
                                                     {},
                                            BRepExt::failureHashBound);
    const TopoDS_Shape& made   = BRepExt::get(session(), result.resultId);
    EXPECT_TRUE(BRepExt::invalidity(made).empty());
    EXPECT_NEAR(volume(made),
                a * b * c - (a + b) * c * c * t + 4. / 3 * c * c * c * t * t,
                1e-9 * a * b * c)
      << degrees << " degrees";
  }
}

// A boss fused onto a box and then defeatured leaves the box: the boss's wall and top are removed
// and the neighbours closed over it. A defeature that leaves a stub or cuts the box changes the
// volume by that of the boss.
TEST(BRepOffsetExt, DefeatureOfABossLeavesTheBox)
{
  const double       radius = 8, height = 6;
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(a, b, c).Shape();
  const TopoDS_Shape boss =
    BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(a / 2, b / 2, c), gp_Dir(0, 0, 1)), radius, height)
      .Shape();
  const TopoDS_Shape fused = BRepAlgoAPI_Fuse(box, boss).Shape();
  const auto         features =
    facesWhere(fused, [](const double* l) { return l[5] > c + 1e-6 && l[2] >= c - 1e-6; });
  ASSERT_EQ(features.size(), 2u);
  const auto result =
    BRepExt::defeature(session(), keep(fused), keep(features), {}, BRepExt::failureHashBound);
  const TopoDS_Shape& made = BRepExt::get(session(), result.resultId);
  EXPECT_TRUE(BRepExt::invalidity(made).empty());
  EXPECT_NEAR(volume(made), a * b * c, 1e-9 * a * b * c);
}

// Hollowing a cylinder of radius R and height H to a wall t leaves the cylinder of radius R - t
// and height H - 2t missing: pi R^2 H - pi (R - t)^2 (H - 2t). The curved wall must be offset
// along its normals and the cavity closed by planar caps at t from the ends.
TEST(BRepOffsetExt, HollowOfACylinderLeavesTheClosedFormWall)
{
  const double       radius = 15, height = 40, t = 2.5;
  const TopoDS_Shape cylinder = BRepPrimAPI_MakeCylinder(radius, height).Shape();
  const auto         result =
    BRepExt::hollow(session(), keep(cylinder), t, 1e-4, {}, BRepExt::failureHashBound);
  const TopoDS_Shape& made = BRepExt::get(session(), result.resultId);
  EXPECT_TRUE(BRepExt::invalidity(made).empty());
  const double expected =
    M_PI * radius * radius * height - M_PI * (radius - t) * (radius - t) * (height - 2 * t);
  EXPECT_NEAR(volume(made), expected, 1e-6 * expected);
}

// Every face of a solid moved out by d: with round joins the result is the parallel body, whose
// volume is Steiner's polynomial abc + 2d(ab + bc + ca) + pi d^2 (a + b + c) + 4/3 pi d^3; with
// sharp joins the box simply grows to (a + 2d)(b + 2d)(c + 2d); moved in by d it shrinks to (a -
// 2d)(b - 2d)(c - 2d). The three differ by far more than the tolerance in the corners.
TEST(BRepOffsetExt, OffsetOfABoxIsSteinersParallelBodyOrTheGrownBox)
{
  const TopoDS_Shape box    = BRepPrimAPI_MakeBox(a, b, c).Shape();
  const double       d      = 4;
  const auto         offset = [&](double distance, BRepExt::OffsetJoin join) {
    const auto result = BRepExt::offsetSolid(session(),
                                             keep(box),
                                             distance,
                                             1e-6,
                                             join,
                                                     {},
                                             BRepExt::failureHashBound);
    return BRepExt::get(session(), result.resultId);
  };
  const TopoDS_Shape rounded = offset(d, BRepExt::OffsetJoin::Arc),
                     sharp   = offset(d, BRepExt::OffsetJoin::Intersection),
                     inward  = offset(-d, BRepExt::OffsetJoin::Intersection);
  EXPECT_TRUE(BRepExt::invalidity(rounded).empty());
  EXPECT_TRUE(BRepExt::invalidity(sharp).empty());
  EXPECT_TRUE(BRepExt::invalidity(inward).empty());
  EXPECT_NEAR(volume(rounded),
              a * b * c + 2 * d * (a * b + b * c + c * a) + M_PI * d * d * (a + b + c)
                + 4. / 3 * M_PI * d * d * d,
              1e-6 * a * b * c);
  EXPECT_NEAR(volume(sharp), (a + 2 * d) * (b + 2 * d) * (c + 2 * d), 1e-9 * a * b * c);
  EXPECT_NEAR(volume(inward), (a - 2 * d) * (b - 2 * d) * (c - 2 * d), 1e-9 * a * b * c);
}

// Every face of a box with one edge rounded by r moved by d, with every face kept:
// the same faces, the sharp edges where the moved planes cross and the round of radius r + d
// meeting them along the offsets of its edges. The volume is that of the box grown by d each way
// less the corner the round leaves out along its edge, now a + 2d long:
// (a + 2d)(b + 2d)(c + 2d) - (1 - pi/4)(r + d)^2 (a + 2d). A round left at r, an edge joined by an
// arc or a vertex left where it was changes it by d^2 times a length or more. The same at a
// thousandth and a thousand times the size, turned and moved off the origin.
TEST(BRepOffsetExt, EveryFaceOfARoundedBoxMovedKeepsItsFacesAndHasTheClosedFormVolume)
{
  const double       r = 6;
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(a, b, c).Shape();
  BRepFilletAPI_MakeFillet round(box);
  for (const TopoDS_Edge& edge : edgesWhere(box, [](const double* l) {
         return near(l[1], 0) && near(l[4], 0) && near(l[2], c) && near(l[5], c);
       }))
    round.Add(r, edge);
  const TopoDS_Shape rounded = round.Shape();
  for (const double s : {1e-3, 1.0, 1e3})
  {
    gp_Trsf turn, shift, scale;
    turn.SetRotation(gp_Ax1(gp::Origin(), gp_Dir(1, 2, 3)), 0.7);
    shift.SetTranslation(gp_Vec(3, -5, 7));
    scale.SetScale(gp::Origin(), s);
    const TopoDS_Shape posed = BRepBuilderAPI_Transform(rounded, scale * shift * turn, true).Shape();
    for (const double d : {2.0, -2.0})
    {
      const TopoDS_Shape made = kept(posed, d * s, 1e-6 * s);
      if (made.IsNull())
        continue;
      EXPECT_TRUE(BRepExt::invalidity(made).empty());
      EXPECT_EQ(faceCount(made), faceCount(posed)) << "s=" << s << " d=" << d;
      const double grown = (a + 2 * d) * (b + 2 * d) * (c + 2 * d)
                           - (1 - M_PI / 4) * (r + d) * (r + d) * (a + 2 * d);
      EXPECT_NEAR(volume(made) / (s * s * s), grown, 1e-9 * a * b * c) << "s=" << s << " d=" << d;
    }
  }
}

// Held out from the box above: prisms over parallelograms of angles alpha and pi - alpha, h high,
// one upright edge rounded by r, every face moved by d. Moving the sides of a convex polygon out
// by d and letting them cross adds P d + d^2 sum cot(a_i / 2) to its area (a_i its angles, P its
// perimeter); a corner of angle a rounded by rho leaves out rho^2 (cot(a/2) - (pi - a)/2), and the
// round's radius becomes r + d. So the volume is
// (h + 2d) (A + P d + d^2 sum cot(a_i / 2) - (r + d)^2 (cot(alpha/2) - (pi - alpha)/2)).
TEST(BRepOffsetExt, EveryFaceOfARoundedPrismMovedHasTheClosedFormVolumeAtAnyAngle)
{
  const double side = 30, h = 20, r = 5;
  for (const double degrees : {50.0, 70.0, 110.0, 130.0})
  {
    const double alpha = degrees * M_PI / 180;
    const gp_Pnt p0(0, 0, 0), p1(side, 0, 0),
      p3(side * std::cos(alpha), side * std::sin(alpha), 0),
      p2(p1.X() + p3.X(), p1.Y() + p3.Y(), 0);
    const TopoDS_Shape prism =
      BRepPrimAPI_MakePrism(
        BRepBuilderAPI_MakeFace(BRepBuilderAPI_MakePolygon(p0, p1, p2, p3, true).Wire()),
        gp_Vec(0, 0, h))
        .Shape();
    BRepFilletAPI_MakeFillet round(prism);
    for (const TopoDS_Edge& edge : edgesWhere(prism, [](const double* l) {
           return near(l[0], 0) && near(l[3], 0) && near(l[1], 0) && near(l[4], 0);
         }))
      round.Add(r, edge);
    const TopoDS_Shape rounded = round.Shape();
    const double area = side * side * std::sin(alpha), perimeter = 4 * side;
    const double corners = 2 / std::tan(alpha / 2) + 2 / std::tan((M_PI - alpha) / 2);
    for (const double d : {1.5, -1.5})
    {
      const TopoDS_Shape made = kept(rounded, d, 1e-6);
      if (made.IsNull())
        continue;
      EXPECT_TRUE(BRepExt::invalidity(made).empty());
      EXPECT_EQ(faceCount(made), faceCount(rounded)) << degrees << " degrees, d=" << d;
      const double cut = (r + d) * (r + d) * (1 / std::tan(alpha / 2) - (M_PI - alpha) / 2);
      EXPECT_NEAR(volume(made),
                  (h + 2 * d) * (area + perimeter * d + d * d * corners - cut),
                  1e-9 * area * h)
        << degrees << " degrees, d=" << d;
    }
  }
}

// Held out: a square prism a x a x h whose corner at the origin is rounded by an arc that meets
// its two sides at a known angle theta instead of tangentially: the circle through (r, 0) and
// (0, r) centred at (c, c), c = r / (1 + tan theta). Its edges along the round claim tangency to
// a tolerance tau (as an input that is G1 only to a few 1e-5 rad does). Moved by d, the faces
// there come apart by |d| 2 sin(theta / 2), and the result's edge carries tau + |d| theta, the
// input's tolerance and what the move adds, where |d| theta is within tau; past that the edge is
// a crossing and keeps the tolerance. At three scales, both ways.
TEST(BRepOffsetExt, EveryFaceOfAPrismWithANearlyTangentRoundCarriesTheInputsTolerance)
{
  for (const double scale : {1e-3, 1.0, 1e3})
    for (const double theta : {1e-5, 3e-5, 1e-4, 1e-3})
      for (const double d : {0.5, -0.5})
      {
        const double a = 30 * scale, h = 20 * scale, r = 5 * scale, tau = 1e-4 * scale;
        const double c = r / (1 + std::tan(theta)), rho = std::hypot(r - c, c);
        const gp_Pnt mid(c - rho / std::sqrt(2.0), c - rho / std::sqrt(2.0), 0);
        BRepBuilderAPI_MakeWire wire;
        wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(0, r, 0), gp_Pnt(0, a, 0)));
        wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(0, a, 0), gp_Pnt(a, a, 0)));
        wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(a, a, 0), gp_Pnt(a, 0, 0)));
        wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(a, 0, 0), gp_Pnt(r, 0, 0)));
        wire.Add(BRepBuilderAPI_MakeEdge(
          GC_MakeArcOfCircle(gp_Pnt(r, 0, 0), mid, gp_Pnt(0, r, 0)).Value()));
        const TopoDS_Shape prism =
          BRepPrimAPI_MakePrism(BRepBuilderAPI_MakeFace(wire.Wire()), gp_Vec(0, 0, h)).Shape();
        // The two edges along the round: lines between the cylinder and a plane.
        BRep_Builder builder;
        const auto   alongRound = [](const TopoDS_Shape& shape) {
          std::vector<TopoDS_Edge> found;
          for (TopExp_Explorer f(shape, TopAbs_FACE); f.More(); f.Next())
          {
            if (!BRep_Tool::Surface(TopoDS::Face(f.Current()))
                   ->IsKind(STANDARD_TYPE(Geom_CylindricalSurface)))
              continue;
            for (TopExp_Explorer e(f.Current(), TopAbs_EDGE); e.More(); e.Next())
            {
              double first, last;
              const Handle(Geom_Curve) curve = BRep_Tool::Curve(TopoDS::Edge(e.Current()), first, last);
              if (curve->Value(first).Z() != curve->Value(last).Z())
                found.push_back(TopoDS::Edge(e.Current()));
            }
          }
          return found;
        };
        const std::vector<TopoDS_Edge> claimed = alongRound(prism);
        ASSERT_EQ(claimed.size(), 2u);
        for (const TopoDS_Edge& edge : claimed)
        {
          builder.UpdateEdge(edge, tau);
          for (TopExp_Explorer v(edge, TopAbs_VERTEX); v.More(); v.Next())
            builder.UpdateVertex(TopoDS::Vertex(v.Current()), tau);
        }
        const double       spread = std::abs(d) * scale * theta;
        const TopoDS_Shape made   = kept(prism, d * scale, 1e-6 * scale);
        ASSERT_FALSE(made.IsNull()) << scale << " " << theta << " " << d;
        EXPECT_TRUE(BRepExt::invalidity(made).empty());
        EXPECT_EQ(faceCount(made), faceCount(prism));
        const std::vector<TopoDS_Edge> moved = alongRound(made);
        ASSERT_EQ(moved.size(), 2u);
        for (const TopoDS_Edge& edge : moved)
        {
          if (spread > tau)
          {
            EXPECT_LE(BRep_Tool::Tolerance(edge), tau * (1 + 1e-9)) << scale << " " << theta << " " << d;
            continue;
          }
          EXPECT_NEAR(BRep_Tool::Tolerance(edge), tau + spread, 1e-9 * tau)
            << scale << " " << theta << " " << d;
          // The two faces' points of the edge, apart by what the move adds.
          std::vector<TopoDS_Face> faces;
          for (TopExp_Explorer f(made, TopAbs_FACE); f.More(); f.Next())
            for (TopExp_Explorer e(f.Current(), TopAbs_EDGE); e.More(); e.Next())
              if (e.Current().IsSame(edge))
              {
                faces.push_back(TopoDS::Face(f.Current()));
                break;
              }
          ASSERT_EQ(faces.size(), 2u);
          double first, last, apart = 0;
          BRep_Tool::Range(edge, first, last);
          for (int k = 0; k <= 8; ++k)
          {
            const double t = first + (last - first) * k / 8;
            gp_Pnt       on[2];
            for (int j = 0; j < 2; ++j)
            {
              double                   c0, c1;
              const Handle(Geom2d_Curve) uv = BRep_Tool::CurveOnSurface(edge, faces[j], c0, c1);
              on[j] = BRep_Tool::Surface(faces[j])->Value(uv->Value(t).X(), uv->Value(t).Y());
            }
            apart = std::max(apart, on[0].Distance(on[1]));
          }
          EXPECT_NEAR(apart, 2 * std::abs(d) * scale * std::sin(theta / 2), 1e-3 * spread)
            << scale << " " << theta << " " << d;
        }
      }
}

// A sheet of two rectangles L x w meeting at a right angle along an edge of length w, thickened
// by t. Offset to the inside of the bend, the two slabs overlap in a cube t x t x w, so the solid
// is their union: 2 L w t - t^2 w. Offset to the outside, the slabs leave a gap at the crease that
// an arc join fills with a quarter cylinder of radius t: 2 L w t + pi t^2 w / 4. The two signs of
// the thickness give these two volumes. A crease joined by a chamfer or left open changes the
// outside volume by the difference between that quarter cylinder and a triangle or nothing.
TEST(BRepOffsetExt, ThickenOfABentSheetFillsOrTrimsTheCreaseByTheArcJoin)
{
  const double          L = 12, w = 7, t = 2;
  BRepBuilderAPI_Sewing sewing(1e-7);
  sewing.Add(BRepBuilderAPI_MakeFace(BRepBuilderAPI_MakePolygon(gp_Pnt(0, 0, 0),
                                                                gp_Pnt(L, 0, 0),
                                                                gp_Pnt(L, w, 0),
                                                                gp_Pnt(0, w, 0),
                                                                true)
                                       .Wire())
               .Face());
  sewing.Add(BRepBuilderAPI_MakeFace(BRepBuilderAPI_MakePolygon(gp_Pnt(0, 0, 0),
                                                                gp_Pnt(0, 0, L),
                                                                gp_Pnt(0, w, L),
                                                                gp_Pnt(0, w, 0),
                                                                true)
                                       .Wire())
               .Face());
  sewing.Perform();
  const uint32_t id      = BRepExt::store(session(), sewing.SewedShape());
  const double   inside  = 2 * L * w * t - t * t * w;
  const double   outside = 2 * L * w * t + M_PI * t * t * w / 4;
  double         got[2];
  int            k = 0;
  for (const double sign : {1.0, -1.0})
  {
    const auto result =
      BRepExt::thicken(session(), id, sign * t, 1e-7, {}, BRepExt::failureHashBound);
    got[k++] = std::abs(volume(BRepExt::get(session(), result.resultId)));
  }
  const auto close = [](double x, double y) { return std::abs(x - y) < 1e-6 * y; };
  EXPECT_TRUE((close(got[0], inside) && close(got[1], outside))
              || (close(got[0], outside) && close(got[1], inside)))
    << got[0] << " " << got[1] << " expected " << inside << " and " << outside;
}

// The same bend as one face: a degree 1 surface with an interior knot, so C0 along the fold. It is
// divided along the fold before it is offset, which gives the two faces of the sheet above, and
// then has the same two volumes; an offset of the undivided surface is not defined at the knot
// and is refused or wrong.
TEST(BRepOffsetExt, ThickenOfOneFoldedFaceIsDividedAtTheFold)
{
  const double               L = 12, w = 7, t = 2;
  NCollection_Array2<gp_Pnt> poles(1, 3, 1, 2);
  for (int j = 1; j <= 2; ++j)
  {
    const double y = (j - 1) * w;
    poles(1, j)    = gp_Pnt(0, y, L);
    poles(2, j)    = gp_Pnt(0, y, 0);
    poles(3, j)    = gp_Pnt(L, y, 0);
  }
  NCollection_Array1<double> uKnots(1, 3), vKnots(1, 2);
  NCollection_Array1<int>    uMult(1, 3), vMult(1, 2);
  uKnots(1) = 0;
  uKnots(2) = 1;
  uKnots(3) = 2;
  uMult(1)  = 2;
  uMult(2)  = 1;
  uMult(3)  = 2;
  vKnots(1) = 0;
  vKnots(2) = 1;
  vMult(1)  = 2;
  vMult(2)  = 2;
  Handle(Geom_BSplineSurface) fold =
    new Geom_BSplineSurface(poles, uKnots, vKnots, uMult, vMult, 1, 1);
  const uint32_t id      = BRepExt::store(session(), BRepBuilderAPI_MakeFace(fold, 1e-7).Face());
  const double   inside  = 2 * L * w * t - t * t * w;
  const double   outside = 2 * L * w * t + M_PI * t * t * w / 4;
  double         got[2];
  int            k = 0;
  for (const double sign : {1.0, -1.0})
  {
    const auto result =
      BRepExt::thicken(session(), id, sign * t, 1e-7, {}, BRepExt::failureHashBound);
    got[k++] = std::abs(volume(BRepExt::get(session(), result.resultId)));
  }
  const auto close = [](double x, double y) { return std::abs(x - y) < 1e-6 * y; };
  EXPECT_TRUE((close(got[0], inside) && close(got[1], outside))
              || (close(got[0], outside) && close(got[1], inside)))
    << got[0] << " " << got[1] << " expected " << inside << " and " << outside;
}

// Two bosses on one box, removed in one call: the volume is that of the box, and removing only
// the faces of one leaves the other boss standing (box plus its cylinder, pi r^2 h). A group that
// stops after the first feature, or extends the neighbours of one over the other, shows in the
// volume at first order.
TEST(BRepOffsetExt, DefeatureOfTwoBossesAtOnceLeavesTheBoxAndOneAtATimeLeavesTheOther)
{
  const double       radius = 6, height = 5;
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(a, b, c).Shape();
  const TopoDS_Shape bossA =
    BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(a / 4, b / 2, c), gp_Dir(0, 0, 1)), radius, height)
      .Shape();
  const TopoDS_Shape bossB =
    BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(3 * a / 4, b / 2, c), gp_Dir(0, 0, 1)), radius, height)
      .Shape();
  const TopoDS_Shape fused = BRepAlgoAPI_Fuse(BRepAlgoAPI_Fuse(box, bossA).Shape(), bossB).Shape();
  const auto         left  = facesWhere(fused, [](const double* l) {
    return l[5] > c + 1e-6 && l[2] >= c - 1e-6 && l[3] < a / 2;
  });
  const auto         right = facesWhere(fused, [](const double* l) {
    return l[5] > c + 1e-6 && l[2] >= c - 1e-6 && l[0] > a / 2;
  });
  ASSERT_EQ(left.size(), 2u);
  ASSERT_EQ(right.size(), 2u);
  std::vector<TopoDS_Face> both = left;
  both.insert(both.end(), right.begin(), right.end());
  const double bossVolume = M_PI * radius * radius * height;
  {
    const auto result =
      BRepExt::defeature(session(), keep(fused), keep(both), {}, BRepExt::failureHashBound);
    const TopoDS_Shape& made = BRepExt::get(session(), result.resultId);
    EXPECT_TRUE(BRepExt::invalidity(made).empty());
    EXPECT_NEAR(volume(made), a * b * c, 1e-9 * a * b * c);
  }
  {
    const auto result =
      BRepExt::defeature(session(), keep(fused), keep(left), {}, BRepExt::failureHashBound);
    const TopoDS_Shape& made = BRepExt::get(session(), result.resultId);
    EXPECT_TRUE(BRepExt::invalidity(made).empty());
    EXPECT_NEAR(volume(made), a * b * c + bossVolume, 1e-9 * a * b * c);
  }
}

// Quadratic graphs have principal curvatures 2a/s and 2b/s at their centre. Convex,
// saddle and cylindrical graphs, positive rational weights, sizes, poses and orientations
// exercise the signed normal sweep, independently of the solid builder's validity test.
TEST(BRepOffsetExt, SplineFocalWitnessIsSignedAndIndependentOfSizeAndPose)
{
  for (double scale : {0.01, 1.0, 100.0})
    for (double b : {-0.3, 0.0, 0.3})
      for (bool rational : {false, true})
        for (bool posed : {false, true})
        {
          NCollection_Array2<gp_Pnt> poles(1, 3, 1, 3);
          NCollection_Array2<double> weights(1, 3, 1, 3);
          const double               x[]{-1, 0, 1}, z[]{1, -1, 1};
          for (int i = 1; i <= 3; ++i)
            for (int j = 1; j <= 3; ++j)
            {
              poles(i, j) =
                gp_Pnt(scale * x[i - 1], scale * x[j - 1], scale * (0.5 * z[i - 1] + b * z[j - 1]));
              weights(i, j) = rational && i == 2 && j == 2 ? 1.125 : 1.0;
            }
          Handle(Geom_BezierSurface)  patch  = new Geom_BezierSurface(poles, weights);
          Handle(Geom_BSplineSurface) spline = GeomConvert::SurfaceToBSplineSurface(patch);
          // Extra spans exercise repeated knots without converting before the proof.
          spline->InsertUKnot(0.37, 1, 1e-12);
          spline->InsertVKnot(0.61, 1, 1e-12);
          if (posed)
          {
            gp_Trsf turn;
            turn.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 2, 3)), 0.71);
            turn.SetTranslationPart(gp_Vec(7 * scale, -3 * scale, 11 * scale));
            spline->Transform(turn);
          }
          const TopoDS_Face face = BRepBuilderAPI_MakeFace(spline, 1e-9 * scale).Face();
          EXPECT_EQ(BRepExt::focalSweep(face, 1.5 * scale), BRepExt::FocalSweep::Crosses);
          try
          {
            BRepExt::thicken(session(),
                             keep(face),
                             1.5 * scale,
                             1e-8 * scale,
                             {},
                             BRepExt::failureHashBound);
            ADD_FAILURE() << "a focal crossing was accepted";
          }
          catch (const BRepExt::Failure& failure)
          {
            EXPECT_EQ(failure.operation, "thicken");
            EXPECT_EQ(failure.reason, BRepExt::no_solution);
            EXPECT_EQ(failure.kind, BRepExt::FailureKind::Curvature);
            ASSERT_EQ(failure.shapes.size(), 1u);
            EXPECT_TRUE(failure.shapes[0].IsSame(face));
          }
          EXPECT_EQ(BRepExt::focalSweep(TopoDS::Face(face.Reversed()), -1.5 * scale),
                    BRepExt::FocalSweep::Crosses);
          EXPECT_EQ(BRepExt::focalSweep(face, 0.05 * scale), BRepExt::FocalSweep::Regular);
          gp_Trsf placement;
          placement.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(3, 1, 2)), 0.43);
          placement.SetTranslationPart(gp_Vec(-13 * scale, 4 * scale, 7 * scale));
          const auto located = TopoDS::Face(face.Moved(TopLoc_Location(placement)));
          EXPECT_EQ(BRepExt::focalSweep(located, 1.5 * scale), BRepExt::FocalSweep::Crosses);
          EXPECT_EQ(BRepExt::focalSweep(located, 0.05 * scale), BRepExt::FocalSweep::Regular);
          // The outward side of a convex graph has no positive focal distance.
          if (b >= 0 && !rational)
            EXPECT_EQ(BRepExt::focalSweep(face, -1.5 * scale), BRepExt::FocalSweep::Regular);
        }
}

// Independent meridian integral: 2 pi times the first moment of the disk of radius
// a centred at (R,0), clipped at x=0. This is a parallel BODY after the hole closes.
TEST(BRepOffsetExt, RingParallelBodyClosesItsHoleWithoutKeepingTheCrossingSheet)
{
  for (double scale : {0.1, 1.0, 10.0})
    for (double ratio : {1.001, 1.1, 2.0})
      for (bool posed : {false, true})
      {
        const double R = 10 * scale, r = 4 * scale, a = ratio * R, d = a - r;
        gp_Ax2       axes(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
        if (posed)
          axes = gp_Ax2(gp_Pnt(7 * scale, -3 * scale, 11 * scale), gp_Dir(1, 2, 3));
        const TopoDS_Shape source = BRepPrimAPI_MakeTorus(axes, R, r).Shape();
        const double       before = volume(source);
        const auto         result = BRepExt::offsetSolid(session(),
                                                 keep(source),
                                                 d,
                                                 1e-7 * scale,
                                                 BRepExt::OffsetJoin::Arc,
                                                         {},
                                                 BRepExt::failureHashBound);
        const TopoDS_Shape made   = BRepExt::get(session(), result.resultId);
        EXPECT_TRUE(BRepCheck_Analyzer(made).IsValid());
        EXPECT_EQ(faceCount(made), 1);
        const double z = std::sqrt((a - R) * (a + R));
        const double expected =
          4 * M_PI
          * (z * z * z / 3 + R * (M_PI * a * a / 4 + R * z / 2 + a * a * std::asin(R / a) / 2));
        EXPECT_NEAR(volume(made), expected, expected * 1e-8);
        EXPECT_EQ(volume(source), before);
      }
}


// A valid imported C0 support may already agree with its ideal plane within the input
// tolerance, while its normal changes slightly at a knot. The offset must extend that support
// to meet the adjacent moved planes. A trimming wrapper rejects the extended UV box before
// approximation starts; an adaptor evaluates the same end-span continuation used by the edges.
TEST(BRepOffsetExt, NormalOffsetApproximationCoversTheExtendedFaceDomain)
{
  for (const double scale : {0.1, 1.0, 10.0})
  {
    const double side = 20 * scale, ridge = 1e-8 * scale, tolerance = 1e-6 * scale;
    TopoDS_Shape body = BRepPrimAPI_MakeBox(side, side, side).Shape();
    BRep_Builder builder;
    for (TopExp_Explorer f(body, TopAbs_FACE); f.More(); f.Next())
    {
      const TopoDS_Face face = TopoDS::Face(f.Current());
      const Handle(Geom_Surface) plane = BRep_Tool::Surface(face);
      if (std::abs(plane->Value(side / 2, side / 2).Z() - side) > tolerance)
        continue;
      NCollection_Array2<gp_Pnt> poles(1, 3, 1, 2);
      NCollection_Array1<double> uk(1, 3), vk(1, 2);
      NCollection_Array1<int> um(1, 3), vm(1, 2);
      uk(1) = vk(1) = 0.0;
      uk(2) = side / 2;
      uk(3) = vk(2) = side;
      um(1) = um(3) = vm(1) = vm(2) = 2;
      um(2) = 1;
      for (int i = 1; i <= 3; ++i)
        for (int j = 1; j <= 2; ++j)
        {
          gp_Pnt p = plane->Value(uk(i), vk(j));
          if (i == 2)
            p.SetZ(p.Z() + ridge);
          poles(i, j) = p;
        }
      const Handle(Geom_Surface) support = new Geom_BSplineSurface(poles, uk, vk, um, vm, 1, 1);
      struct OnFace
      {
        TopoDS_Edge edge;
        Handle(Geom2d_Curve) curve;
        double first, last;
      };
      std::vector<OnFace> curves;
      for (TopExp_Explorer e(face, TopAbs_EDGE); e.More(); e.Next())
      {
        OnFace on;
        on.edge = TopoDS::Edge(e.Current());
        on.curve = BRep_Tool::CurveOnSurface(on.edge, face, on.first, on.last);
        curves.push_back(on);
      }
      builder.UpdateFace(face, support, TopLoc_Location(), 1e-7 * scale);
      for (const OnFace& on : curves)
      {
        builder.UpdateEdge(on.edge, on.curve, face, 1e-7 * scale);
        builder.Range(on.edge, face, on.first, on.last);
      }
    }
    gp_Trsf pose;
    pose.SetRotation(gp_Ax1(gp::Origin(), gp_Dir(1, 2, 3)), 0.7);
    body = BRepBuilderAPI_Transform(body, pose, true).Shape();
    ASSERT_TRUE(BRepCheck_Analyzer(body).IsValid());
    for (const double d : {-0.5 * scale, 0.05 * scale, 0.5 * scale})
    {
      BRepExt::FaceSet selected;
      for (TopExp_Explorer f(body, TopAbs_FACE); f.More(); f.Next())
        selected.Add(f.Current());
      BRepExt::Lineage lineage(body);
      std::string why;
      const TopoDS_Shape result = BRepExt::keptOffset(body, selected, d, tolerance, lineage, why);
      ASSERT_FALSE(result.IsNull()) << why << " scale=" << scale << " d=" << d;
      EXPECT_TRUE(BRepCheck_Analyzer(result).IsValid());
      EXPECT_EQ(faceCount(result), faceCount(body));
      // The input differs from the exact cube by at most ridge; its parallel boundary differs
      // by ridge plus |d| times the normal's change, below tolerance. The resulting volume is
      // therefore within its area times tolerance of the expanded/contracted cube.
      EXPECT_NEAR(volume(result), std::pow(side + 2 * d, 3),
                  6 * std::pow(side + 2 * d, 2) * tolerance);
      for (TopExp_Explorer f(body, TopAbs_FACE); f.More(); f.Next())
      {
        const TopoDS_Face source = TopoDS::Face(f.Current());
        const auto successors = lineage.of(source);
        ASSERT_EQ(successors.size(), 1u);
        const TopoDS_Face made = TopoDS::Face(successors[0]);
        double u0, u1, v0, v1;
        BRepTools::UVBounds(made, u0, u1, v0, v1);
        const Handle(Geom_Surface) before = BRep_Tool::Surface(source), after = BRep_Tool::Surface(made);
        for (const double x : {0.1, 0.3, 0.7, 0.9})
          for (const double y : {0.1, 0.3, 0.7, 0.9})
          {
            const double u = u0 + x * (u1 - u0), v = v0 + y * (v1 - v0);
            gp_Pnt p;
            gp_Vec su, sv;
            before->D1(u, v, p, su, sv);
            const double along = source.Orientation() == TopAbs_REVERSED ? -d : d;
            const gp_Pnt exact = p.Translated(su.Crossed(sv).Normalized() * along);
            EXPECT_LE(after->Value(u, v).Distance(exact), tolerance);
          }
      }
    }
  }
}

// A finite spline pcurve on a periodic seam has two distinct parameter branches. Moving
// the end planes extends both branches, without moving either to the other periodic bound.
// This changes only the representation of the cylinder's seam, not its geometry.
TEST(BRepOffsetExt, FiniteSeamBranchesExtendWithTheirEndPlanes)
{
  const double radius = 3.0, height = 7.0;
  TopoDS_Shape body = BRepPrimAPI_MakeCylinder(radius, height).Shape();
  BRep_Builder builder;
  for (TopExp_Explorer f(body, TopAbs_FACE); f.More(); f.Next())
  {
    const TopoDS_Face face = TopoDS::Face(f.Current());
    for (TopExp_Explorer e(face, TopAbs_EDGE); e.More(); e.Next())
    {
      const TopoDS_Edge edge = TopoDS::Edge(e.Current());
      if (!BRep_Tool::IsClosed(edge, face))
        continue;
      double first, last;
      BRep_Tool::Range(edge, first, last);
      const auto branch = [&](TopAbs_Orientation orientation) {
        return Geom2dConvert::CurveToBSplineCurve(new Geom2d_TrimmedCurve(
          BRep_Tool::CurveOnSurface(TopoDS::Edge(edge.Oriented(orientation)), face, first, last),
          first, last));
      };
      const Handle(Geom2d_Curve) one = branch(TopAbs_FORWARD), two = branch(TopAbs_REVERSED);
      builder.UpdateEdge(edge, one, two, face, BRep_Tool::Tolerance(edge));
    }
  }
  ASSERT_TRUE(BRepCheck_Analyzer(body).IsValid());
  for (const double scale : {0.01, 1.0, 100.0})
  {
    gp_Trsf turn, shift, size;
    turn.SetRotation(gp_Ax1(gp::Origin(), gp_Dir(1, 2, 3)), 0.47);
    shift.SetTranslation(gp_Vec(2, -3, 5));
    size.SetScale(gp::Origin(), scale);
    const TopoDS_Shape input = BRepBuilderAPI_Transform(body, size * shift * turn, true).Shape();
    for (const double d : {-0.2, 0.2})
    {
      const TopoDS_Shape result = kept(input, d * scale, 1e-6 * scale);
      ASSERT_FALSE(result.IsNull());
      EXPECT_TRUE(BRepCheck_Analyzer(result).IsValid());
      EXPECT_EQ(faceCount(result), 3);
      EXPECT_NEAR(volume(result) / (scale * scale * scale),
                  M_PI * (radius + d) * (radius + d) * (height + 2 * d), 1e-7);
    }
  }
}

TEST(BRepOffsetExt, ExactSphereKeepsEveryPositiveRadiusBeforeItsFocalDistance)
{
  for (const double radius : {100.0, 10000.0, 1000000.0})
  {
    const auto sphere = BRepPrimAPI_MakeSphere(radius).Shape();
    const auto fate = BRepExt::classifyForOffset(sphere, -radius * (1 - 5e-5));
    EXPECT_TRUE(fate.vanish.empty());
    EXPECT_TRUE(fate.partial.empty());
    const auto made  = BRepExt::offsetSolid(session(),
                                           keep(sphere),
                                           -radius * (1 - 5e-5),
                                           1e-6,
                                           BRepExt::OffsetJoin::Intersection,
                                            {},
                                           BRepExt::failureHashBound);
    const auto result=BRepExt::get(session(),made.resultId);
    EXPECT_TRUE(BRepCheck_Analyzer(result).IsValid());
    const double expectedRadius=radius*5e-5;
    EXPECT_NEAR(volume(result),4*std::acos(-1.0)*std::pow(expectedRadius,3)/3,
                std::max(1e-10,1e-7*std::pow(expectedRadius,3)));
    const auto collapsed = BRepExt::classifyForOffset(sphere, -radius);
    EXPECT_EQ(collapsed.vanish.size(), 1u);
  }
}

TEST(BRepOffsetExt, UnmeasuredClearanceIsNeverAnInfiniteSuccess)
{
  const auto box = BRepPrimAPI_MakeBox(10, 10, 10).Shape();
  EXPECT_FALSE(BRepExt::clearance({}, box).satisfies(1));
  EXPECT_FALSE(BRepExt::clearance({}, box, {}, 0, 1).satisfies(1));
  TopoDS_Face face = BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)),
                                               0, 1, 0, 1).Face();
  EXPECT_FALSE(BRepExt::clearance({face}, box, {face}, 100, 1).satisfies(1));
}

TEST(BRepOffsetExt, ClearanceIncludesTheTrimmedBoundaryBetweenSampleSites)
{
  const auto body=BRepPrimAPI_MakeBox(gp_Pnt(0,0,-10),10,10,10).Shape();
  BRepBuilderAPI_MakePolygon boundary;
  for (const auto& point : {gp_Pnt(0,0,0),gp_Pnt(10,0,10),gp_Pnt(10,10,10),gp_Pnt(0,10,0)})
    boundary.Add(point);
  boundary.Close();
  const auto face=BRepBuilderAPI_MakeFace(boundary.Wire()).Face();
  const auto measured=BRepExt::clearance({face},body,{},0,1);
  EXPECT_TRUE(measured.resolved);
  EXPECT_EQ(measured.covered,1u);
  EXPECT_NEAR(measured.minimum,0,1e-9);
  EXPECT_FALSE(measured.satisfies(1));
}

TEST(BRepOffsetExt, SweptFaceUsesSlantedNeighbourSupports)
{
  for (const double scale : {0.1, 1.0, 10.0})
  {
    BRepBuilderAPI_MakePolygon outline;
    for (const auto& p : {gp_Pnt(0,0,0), gp_Pnt(20,0,0), gp_Pnt(20,0,19),
                          gp_Pnt(19,0,20), gp_Pnt(0,0,20)})
      outline.Add(gp_Pnt(p.X()*scale, p.Y()*scale, p.Z()*scale));
    outline.Close();
    const auto body = BRepPrimAPI_MakePrism(BRepBuilderAPI_MakeFace(outline.Wire()).Face(),
                                            gp_Vec(0, 20*scale, 0)).Shape();
    const auto moved = facesWhere(body, [&](const double* b) {
      return std::abs(b[0]-20*scale) < 1e-5 && std::abs(b[3]-20*scale) < 1e-5;
    });
    ASSERT_EQ(moved.size(), 1u);
    for (const double distance : {-2.0, -1.0, 1.0, 3.0})
    {
      BRepExt::Lineage lineage(body);
      std::string why;
      std::vector<TopoDS_Shape> at;
      const auto result = BRepExt::sweepFaces(body, {moved[0]}, distance*scale,
                                             1e-6*scale, lineage, why, at);
      ASSERT_FALSE(result.IsNull()) << why;
      EXPECT_TRUE(BRepCheck_Analyzer(result).IsValid());
      const double width=20+distance;
      const double area=width<=19 ? 20*width : 380+39*(width-19)-(width*width-361)/2;
      EXPECT_NEAR(volume(result),20*area*std::pow(scale,3),1e-7*std::pow(scale,3));
      EXPECT_EQ(lineage.of(moved[0]).size(),1u);
      for (const auto type : {TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX})
      {
        BRepExt::ShapeSet sources, products;
        TopExp::MapShapes(body, type, sources);
        TopExp::MapShapes(result, type, products);
        int retained = 0;
        for (const auto& source : sources)
        {
          const auto& images = lineage.of(source);
          retained += !images.empty();
          for (const auto& image : images) EXPECT_TRUE(products.Contains(image));
          if (type == TopAbs_VERTEX && width >= 19)
          {
            ASSERT_EQ(images.size(), 1u);
            auto expected = BRep_Tool::Pnt(TopoDS::Vertex(source));
            if (near(expected.X(), 20*scale))
            {
              expected.SetX(width*scale);
              if (near(expected.Z(), 19*scale)) expected.SetZ((39-width)*scale);
            }
            EXPECT_LE(expected.Distance(BRep_Tool::Pnt(TopoDS::Vertex(images[0]))),1e-6*scale);
          }
        }
        const int expected = type == TopAbs_FACE ? (width > 19 ? 7 : 6)
                           : type == TopAbs_EDGE ? (width > 19 ? 15 : width == 19 ? 13 : 11)
                           : (width >= 19 ? 10 : 6);
        EXPECT_EQ(retained, expected) << "type=" << type << " distance=" << distance;
      }
    }
  }
}

TEST(BRepOffsetExt, ClearanceUsesSurfaceTrimsInsteadOfTolerancedSpatialEdges)
{
  for (const double scale : {0.1, 1.0, 10.0})
  {
    auto source = BRepBuilderAPI_MakeFace(gp_Pln(gp::Origin(), gp::DZ()),
                                         0, 10 * scale, 0, 10 * scale).Face();
    const auto target = BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 0, scale), gp::DZ()),
                                               0, 10 * scale, 0, 10 * scale).Face();
    BRep_Builder builder;
    for (TopExp_Explorer e(source, TopAbs_EDGE); e.More(); e.Next())
    {
      const auto edge = TopoDS::Edge(e.Current());
      double first, last;
      const auto curve = Handle(Geom_Curve)::DownCast(BRep_Tool::Curve(edge, first, last)->Copy());
      curve->Translate(gp_Vec(0, 0, 0.25 * scale));
      builder.UpdateEdge(edge, curve, 0.3 * scale);
    }
    const auto measured = BRepExt::clearance({target}, source, {}, 0, 0.9 * scale);
    EXPECT_TRUE(measured.satisfies(0.9 * scale)) << measured.issue;
    EXPECT_NEAR(measured.minimum, scale, Precision::Confusion());
    // Measurement copies must not repair or otherwise mutate caller geometry.
    for (TopExp_Explorer e(source, TopAbs_EDGE); e.More(); e.Next())
    {
      double first, last;
      const auto curve = BRep_Tool::Curve(TopoDS::Edge(e.Current()), first, last);
      EXPECT_NEAR(curve->Value((first + last) / 2).Z(), 0.25 * scale, Precision::Confusion());
    }
  }
}

TEST(BRepOffsetExt, ElementaryFocalClassificationIncludesTrimExtrema)
{
  for (const double scale : {0.1, 1.0, 10.0})
  {
    const double alpha = 0.4;
    const Handle(Geom_Surface) cone = new Geom_ConicalSurface(gp::XOY(), alpha, 10 * scale);
    const auto coneFace = BRepBuilderAPI_MakeFace(cone, 0, 2 * M_PI, 0, 100 * scale,
                                                Precision::Confusion()).Face();
    const auto coneFate = BRepExt::classifyForOffset(coneFace, -10.1 * scale / std::cos(alpha));
    EXPECT_EQ(coneFate.partial.size(), 1u);
    EXPECT_TRUE(coneFate.vanish.empty());
    const Handle(Geom_Surface) torus = new Geom_ToroidalSurface(gp::XOY(), 10 * scale, 2 * scale);
    const auto torusFace = BRepBuilderAPI_MakeFace(torus, 0, 2 * M_PI, 0, M_PI,
                                                 Precision::Confusion()).Face();
    const auto torusFate = BRepExt::classifyForOffset(torusFace, 8.001 * scale);
    EXPECT_EQ(torusFate.partial.size(), 1u);
    EXPECT_TRUE(torusFate.vanish.empty());
  }
}

TEST(BRepOffsetExt, PocketedBodyOpeningClearance)
{
  auto body =
    BRepAlgoAPI_Cut(BRepPrimAPI_MakeBox(10, 10, 10).Shape(),
                    BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(5, 5, 0), gp::DZ()), 2, 10).Shape())
      .Shape();
  auto edges = edgesWhere(body, [](const double* b) {
    return near(b[0], 0) && near(b[3], 0) && near(b[1], 0) && near(b[4], 0);
  });
  ASSERT_EQ(edges.size(), 1u);
  BRepFilletAPI_MakeFillet fillet(body);
  fillet.Add(1, edges[0]);
  body  = fillet.Shape();
  body  = BRepAlgoAPI_Cut(body, BRepPrimAPI_MakeBox(gp_Pnt(9, 3, -5), 2, 4, 20).Shape()).Shape();
  edges = edgesWhere(body, [](const double* b) {
    return near(b[1], 10) && near(b[4], 10) && near(b[2], 10) && near(b[5], 10);
  });
  ASSERT_EQ(edges.size(), 1u);
  BRepFilletAPI_MakeChamfer chamfer(body);
  chamfer.Add(1, edges[0]);
  body = chamfer.Shape();
  body = BRepAlgoAPI_Cut(body, BRepPrimAPI_MakeBox(gp_Pnt(1, 7, -1), 2, 2, 2).Shape()).Shape();
  for (const double scale : {0.1, 1.0, 10.0})
    for (const double angle : {0.0, 0.47})
      for (const double wall : {0.1, 0.2, 0.3})
      {
        SCOPED_TRACE(::testing::Message() << scale << ", " << angle << ", " << wall);
        gp_Trsf turn, shift, size;
        turn.SetRotation(gp_Ax1(gp::Origin(), gp_Dir(1, 2, 3)), angle);
        shift.SetTranslation(gp_Vec(2, -3, 5));
        size.SetScale(gp::Origin(), scale);
        const auto input = BRepBuilderAPI_Transform(body, size * shift * turn, true).Shape();
        const auto hash  = [](const TopoDS_Shape& shape) {
          return static_cast<int>(TopTools_ShapeMapHasher{}(shape) % BRepExt::failureHashBound);
        };
        std::vector<int> hashes;
        TopoDS_Shape     opened;
        TopoDS_Compound  skin;
        BRep_Builder     builder;
        builder.MakeCompound(skin);
        for (TopExp_Explorer f(input, TopAbs_FACE); f.More(); f.Next())
        {
          opened = f.Current();
          hashes.push_back(hash(f.Current()));
        }
        for (TopExp_Explorer f(input, TopAbs_FACE); f.More(); f.Next())
          if (!f.Current().IsSame(opened))
            builder.Add(skin, f.Current());
        const auto made   = BRepExt::shell(session(),
                                         keep(input),
                                           {keep(opened)},
                                         wall * scale,
                                         1e-6 * scale,
                                         BRepExt::OffsetJoin::Intersection,
                                         hashes,
                                         BRepExt::failureHashBound,
                                         false);
        const auto result = BRepExt::get(session(), made.resultId);
        EXPECT_TRUE(BRepCheck_Analyzer(result).IsValid());
        EXPECT_LT(volume(result), 0.45 * volume(input));
        size_t walls = 0, rims = 0;
        double minimum = 1e100;
        for (size_t i = 0; i < made.generated.size();)
        {
          const bool rim   = made.generated[i++] == hash(opened);
          const int  count = made.generated[i++];
          for (int j = 0; j < count; ++j)
          {
            const int image = made.generated[i++];
            bool      found = false;
            for (TopExp_Explorer f(result, TopAbs_FACE); f.More(); f.Next())
              if (hash(f.Current()) == image)
              {
                found = true;
                BRepExtrema_DistShapeShape distance(f.Current(), skin);
                ASSERT_TRUE(distance.IsDone());
                if (rim)
                {
                  ++rims;
                  EXPECT_NEAR(distance.Value(), 0, 1e-6 * scale);
                }
                else
                {
                  ++walls;
                  minimum = std::min(minimum, distance.Value());
                  EXPECT_GE(distance.Value(), wall * scale - 1e-6 * scale);
                }
              }
            EXPECT_TRUE(found);
          }
        }
        EXPECT_EQ(walls, 17u);
        EXPECT_EQ(rims, 1u);
        EXPECT_NEAR(minimum, wall * scale, 1e-6 * scale);
      }
}

TEST(BRepOffsetExt, OpeningExclusionKeepsSharedWallDescendants)
{
  const auto      box = BRepPrimAPI_MakeBox(10, 10, 10).Shape();
  TopExp_Explorer f(box, TopAbs_FACE);
  const auto      opened = f.Current();
  f.Next();
  const auto       retained = f.Current();
  BRepExt::Lineage lineage(box);
  lineage.set(opened, {opened, retained});
  const auto walls = lineage.all({opened});
  EXPECT_EQ(walls.size(), 5u);
  EXPECT_TRUE(
    std::any_of(walls.begin(), walls.end(), [&](const auto& s) { return s.IsSame(retained); }));
  // A known thin retained wall cannot be hidden by an opening's provenance.
  const auto thin = BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 0, 0.1), gp::DZ()), 1, 9, 1, 9).Face();
  lineage.set(retained, {thin});
  const auto measured = BRepExt::clearance(lineage.of(retained), box, {}, 0, 0.2);
  EXPECT_TRUE(measured.resolved);
  EXPECT_FALSE(measured.satisfies(0.2));
  EXPECT_NEAR(measured.minimum, 0.1, Precision::Confusion());
}

namespace
{
// Independent of output tolerance metadata. Traverse each oriented use so a closed edge's
// final occurrence and the second pcurve of a periodic seam cannot disappear in a shape map.
void checkEndpointOccurrences(const TopoDS_Shape& shape, double bound)
{
  for (TopExp_Explorer f(shape, TopAbs_FACE); f.More(); f.Next())
  {
    const auto face    = TopoDS::Face(f.Current());
    const auto surface = BRep_Tool::Surface(face);
    for (TopExp_Explorer e(face, TopAbs_EDGE); e.More(); e.Next())
    {
      const auto edge = TopoDS::Edge(e.Current());
      double     first, last;
      const auto curve  = BRep_Tool::Curve(edge, first, last);
      const auto pcurve = BRep_Tool::CurveOnSurface(edge, face, first, last);
      ASSERT_FALSE(pcurve.IsNull());
      for (TopExp_Explorer v(edge, TopAbs_VERTEX); v.More(); v.Next())
      {
        const auto   vertex = TopoDS::Vertex(v.Current());
        const auto   point  = BRep_Tool::Pnt(vertex);
        const double t      = BRep_Tool::Parameter(vertex, edge);
        if (!curve.IsNull())
        {
          EXPECT_LE(point.Distance(curve->Value(t)), bound);
        }
        const auto uv = pcurve->Value(t);
        EXPECT_LE(point.Distance(surface->Value(uv.X(), uv.Y())), bound);
      }
    }
  }
}

void checkJunctionFamily(bool spherical)
{
  // The closed-contact test varies distance, size and pose on one supported blended source.
  // Wider carrier-placement diagnostics include source blends that are not constructible.
  const std::vector<double> centers =
    spherical ? std::vector<double>{5.7} : std::vector<double>{5.7, 6.0, 6.3};
  for (const double dx : centers)
    for (const double scale : {0.1, 1.0, 10.0})
      for (const bool turned : {false, true})
      {
        SCOPED_TRACE(::testing::Message()
                     << spherical << " " << dx << " " << scale << " " << turned);
        const TopoDS_Shape left =
          spherical ? BRepPrimAPI_MakeSphere(11).Shape()
                    : BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -10), gp::DZ()), 10, 20).Shape();
        const TopoDS_Shape right =
          BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(dx, 3, 2), gp::DZ()), 9, 4).Shape();
        const TopoDS_Shape joined = spherical ? BRepAlgoAPI_Common(left, right).Shape()
                                              : BRepAlgoAPI_Fuse(left, right).Shape();
        TopExp_Explorer    solids(joined, TopAbs_SOLID);
        ASSERT_TRUE(solids.More());
        const TopoDS_Shape body = solids.Current();
        NCollection_IndexedDataMap<TopoDS_Shape,
                                   NCollection_List<TopoDS_Shape>,
                                   TopTools_ShapeMapHasher>
          adjacent;
        TopExp::MapShapesAndUniqueAncestors(body, TopAbs_EDGE, TopAbs_FACE, adjacent);
        TopoDS_Shape chosen;
        for (int j = 1; j <= adjacent.Extent(); ++j)
        {
          if (adjacent(j).Size() != 2)
            continue;
          const auto a    = BRepAdaptor_Surface(TopoDS::Face(adjacent(j).First())).GetType();
          const auto b    = BRepAdaptor_Surface(TopoDS::Face(adjacent(j).Last())).GetType();
          const auto kind = spherical ? GeomAbs_Sphere : GeomAbs_Cylinder;
          if ((a == kind && b == GeomAbs_Torus) || (b == kind && a == GeomAbs_Torus))
          {
            chosen = adjacent.FindKey(j);
            break;
          }
        }
        ASSERT_FALSE(chosen.IsNull());
        BRepExt::blending::Request request;
        request.body        = body;
        request.edges       = {chosen};
        request.first       = 0.5;
        TopoDS_Shape source = BRepExt::blending::blend(request).shape;
        gp_Trsf      scaling;
        scaling.SetScale(gp::Origin(), scale);
        source = BRepBuilderAPI_Transform(source, scaling, true).Shape();
        if (turned)
        {
          gp_Trsf pose;
          pose.SetRotation(gp_Ax1(gp::Origin(), gp_Dir(1, 2, 3)), 0.7);
          source = BRepBuilderAPI_Transform(source, pose, true).Shape();
        }
        ASSERT_TRUE(BRepCheck_Analyzer(source).IsValid());
        const double tolerance = 1e-6 * scale;
        // Keep planar end caps fixed. This isolates junction/contact reconstruction from
        // the separate global all-face clearance enclosure on large trimmed spline faces.
        BRepExt::FaceSet moved;
        for (TopExp_Explorer f(source, TopAbs_FACE); f.More(); f.Next())
          if (spherical
              || BRepAdaptor_Surface(TopoDS::Face(f.Current())).GetType() != GeomAbs_Plane)
            moved.Add(f.Current());
        const std::vector<double> distances =
          spherical ? std::vector<double>{-0.025, -0.05, -0.1} : std::vector<double>{-0.05};
        for (const double distance : distances)
        {
          BRepExt::Lineage   lineage(source);
          std::string        why;
          const TopoDS_Shape result =
            BRepExt::keptOffset(source, moved, distance * scale, tolerance, lineage, why);
          ASSERT_FALSE(result.IsNull()) << why;
          EXPECT_TRUE(BRepCheck_Analyzer(result).IsValid());
          EXPECT_EQ(faceCount(result), faceCount(source));
          checkEndpointOccurrences(result, tolerance);
        }
      }
}
} // namespace

TEST(BRepOffsetExt, PeriodicChartCutsFollowTheirIncidentJunctions)
{
  checkJunctionFamily(false);
}

TEST(BRepOffsetExt, ClosedIntersectionOccurrencesShareTheirJunction)
{
  checkJunctionFamily(true);
}

namespace
{
TopoDS_Shape separatedRevolutionContact(double bulge,
                                        double top,
                                        bool   concave = false,
                                        bool   conical = false)
{
  NCollection_Array1<gp_Pnt> pts(1, 4);
  pts(1)                           = gp_Pnt(6, 0, -10);
  pts(2)                           = gp_Pnt(bulge, 0, -4);
  pts(3)                           = gp_Pnt(11, 0, 5);
  pts(4)                           = gp_Pnt(top, 0, 10);
  Handle(Geom_BezierCurve) profile = new Geom_BezierCurve(pts);
  BRepBuilderAPI_MakeWire  wire;
  wire.Add(BRepBuilderAPI_MakeEdge(profile));
  wire.Add(BRepBuilderAPI_MakeEdge(pts(4), gp_Pnt(0, 0, 10)));
  wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(0, 0, 10), gp_Pnt(0, 0, -10)));
  wire.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(0, 0, -10), pts(1)));
  auto right =
    BRepPrimAPI_MakeRevol(BRepBuilderAPI_MakeFace(wire.Wire()).Face(), gp::OZ(), 2 * M_PI).Shape();
  gp_Trsf shift;
  shift.SetTranslation(gp_Vec(6, 3, 2));
  right = BRepBuilderAPI_Transform(right, shift, true).Shape();
  const auto left =
    conical ? BRepPrimAPI_MakeCone(gp_Ax2(gp_Pnt(0, 0, -10), gp::DZ()), 12, 6, 20).Shape()
            : BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -10), gp::DZ()), 10, 20).Shape();
  auto common =
    concave ? BRepAlgoAPI_Fuse(left, right).Shape() : BRepAlgoAPI_Common(left, right).Shape();
  TopExp_Explorer solids(common, TopAbs_SOLID);
  if (!solids.More())
    throw std::runtime_error("no solid");
  auto body = solids.Current();
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    ef;
  TopExp::MapShapesAndUniqueAncestors(body, TopAbs_EDGE, TopAbs_FACE, ef);
  TopoDS_Shape chosen;
  double       longest = 0;
  const auto   kind    = conical ? GeomAbs_Cone : GeomAbs_Cylinder;
  for (int i = 1; i <= ef.Extent(); ++i)
  {
    if (ef(i).Size() != 2)
      continue;
    auto a = BRepAdaptor_Surface(TopoDS::Face(ef(i).First())).GetType(),
         b = BRepAdaptor_Surface(TopoDS::Face(ef(i).Last())).GetType();
    if (!((a == kind && b == GeomAbs_SurfaceOfRevolution)
          || (b == kind && a == GeomAbs_SurfaceOfRevolution)))
      continue;
    GProp_GProps props;
    BRepGProp::LinearProperties(ef.FindKey(i), props);
    if (props.Mass() > longest)
    {
      longest = props.Mass();
      chosen  = ef.FindKey(i);
    }
  }
  if (chosen.IsNull())
    throw std::runtime_error("no contact");
  BRepExt::blending::Request request;
  request.body  = body;
  request.edges = {chosen};
  request.first = .5;
  auto source   = BRepExt::blending::blend(request).shape;

  return source;
}

void checkContactBoundaryLifts(const TopoDS_Shape& result, double bound)
{
  for (TopExp_Explorer f(result, TopAbs_FACE); f.More(); f.Next())
  {
    const auto face = TopoDS::Face(f.Current());
    for (TopExp_Explorer e(face, TopAbs_EDGE); e.More(); e.Next())
    {
      const auto edge = TopoDS::Edge(e.Current());
      double     first, last, lo, hi;
      const auto curve = BRep_Tool::Curve(edge, first, last);
      const auto pc    = BRep_Tool::CurveOnSurface(edge, face, lo, hi);
      ASSERT_FALSE(pc.IsNull());
      if (curve.IsNull())
        continue;
      for (int k = 0; k <= 128; ++k)
      {
        const double t  = first + (last - first) * k / 128.;
        const auto   uv = pc->Value(t);
        EXPECT_LE(curve->Value(t).Distance(BRep_Tool::Surface(face)->Value(uv.X(), uv.Y())), bound);
      }
    }
  }
}

} // namespace

TEST(BRepOffsetExt, SeparatedTangentContactsAddBoundedStrips)
{
  for (const double bulge : {12.0, 14.0})
    for (const double top : {7.0, 9.0})
    {
      const auto source = separatedRevolutionContact(bulge, top);
      for (const double scale : (bulge == 12.0 && top == 7.0 ? std::vector<double>{0.1, 1.0, 10.0}
                                                             : std::vector<double>{1.0}))
        for (const bool turned : {false, true})
        {
          SCOPED_TRACE(::testing::Message()
                       << bulge << " " << top << " " << scale << " " << turned);
          gp_Trsf size;
          size.SetScale(gp::Origin(), scale);
          auto body = BRepBuilderAPI_Transform(source, size, true).Shape();
          if (turned)
          {
            gp_Trsf pose;
            pose.SetRotation(gp_Ax1(gp::Origin(), gp_Dir(1, 2, 3)), 0.6);
            pose.SetTranslationPart(gp_Vec(13, -7, 19) * scale);
            body = BRepBuilderAPI_Transform(body, pose, true).Shape();
          }
          const double bound  = 1e-6 * scale;
          const auto   result = kept(body, 0.05 * scale, bound);
          ASSERT_FALSE(result.IsNull());
          EXPECT_GT(faceCount(result), faceCount(body));
          EXPECT_TRUE(BRepCheck_Analyzer(result).IsValid());
          EXPECT_GT(volume(result), volume(body));
          checkEndpointOccurrences(result, bound);
          checkContactBoundaryLifts(result, bound);
        }
    }
}

TEST(BRepOffsetExt, ContactShellRetainsUnchangedFaceOrientations)
{
  for (const double bulge : {12.0, 14.0})
    for (const double top : {7.0, 9.0})
      for (const bool conical : {false, true})
      {
        const auto source = separatedRevolutionContact(bulge, top, true, conical);
        ASSERT_TRUE(BRepCheck_Analyzer(source).IsValid());
        for (const double scale : (bulge == 12.0 && top == 7.0 ? std::vector<double>{0.1, 1.0, 10.0}
                                                               : std::vector<double>{1.0}))
          for (const bool turned : {false, true})
          {
            SCOPED_TRACE(::testing::Message()
                         << bulge << " " << top << " " << conical << " " << scale << " " << turned);
            gp_Trsf size;
            size.SetScale(gp::Origin(), scale);
            auto body = BRepBuilderAPI_Transform(source, size, true).Shape();
            if (turned)
            {
              gp_Trsf pose;
              pose.SetRotation(gp_Ax1(gp::Origin(), gp_Dir(1, 2, 3)), 0.6);
              pose.SetTranslationPart(gp_Vec(13, -7, 19) * scale);
              body = BRepBuilderAPI_Transform(body, pose, true).Shape();
            }
            const double bound  = 1e-6 * scale;
            const auto   result = kept(body, 0.05 * scale, bound);
            ASSERT_FALSE(result.IsNull());
            EXPECT_TRUE(BRepCheck_Analyzer(result).IsValid());
            EXPECT_GT(faceCount(result), faceCount(body));
            EXPECT_GT(volume(result), volume(body));
            checkEndpointOccurrences(result, bound);
            checkContactBoundaryLifts(result, bound);
          }
      }
}

// The image of an affine iso trim depends on its coordinate endpoints, not on
// the inherited edge parameter. Compare with direct support evaluation, at
// scales and poses that exercise coefficient rounding as well as both senses.
TEST(BRepOffsetExt, AnalyticTrimImagesEncloseTheFullAffineParameterImage)
{
  for (double scale : {0.01, 1.0, 100.0})
    for (bool turned : {false, true})
    {
      const gp_Ax3                            axes(gp_Pnt(13 * scale, -7 * scale, 19 * scale),
                                                   turned ? gp_Dir(1, 2, 3) : gp::DZ());
      const std::vector<Handle(Geom_Surface)> supports{
        new Geom_Plane(axes),
        new Geom_CylindricalSurface(axes, 2 * scale),
        new Geom_ConicalSurface(axes, 0.4, 2 * scale),
        new Geom_SphericalSurface(axes, 2 * scale),
        new Geom_ToroidalSurface(axes, 4 * scale, scale)};
      for (const auto& support : supports)
        for (bool uIso : {false, true})
          for (double sense : {-1.0, 1.0})
          {
            const Handle(Geom2d_Curve) pcurve =
              new Geom2d_Line(uIso ? gp_Pnt2d(0.7, -0.2) : gp_Pnt2d(-0.2, 0.7),
                              uIso ? gp_Dir2d(0, sense) : gp_Dir2d(sense, 0));
            const double first = -0.5, last = 0.8;
            const auto   image = BRepExt::isoTrimImage(support, pcurve, first, last, 1e-8);
            ASSERT_TRUE(image.has_value());
            EXPECT_TRUE(GeomAdaptor_Curve(image->curve).GetType() == GeomAbs_Line
                        || GeomAdaptor_Curve(image->curve).GetType() == GeomAbs_Circle);
            EXPECT_LE(image->error, 1e-8);
            for (int i = 0; i <= 32; ++i)
            {
              const double t  = first + (last - first) * i / 32;
              const auto   uv = pcurve->Value(t);
              const double at = uIso ? uv.Y() : uv.X();
              EXPECT_LE(image->first, at);
              EXPECT_GE(image->last, at);
              EXPECT_LE(image->curve->Value(at).Distance(support->Value(uv.X(), uv.Y())),
                        image->error);
            }
          }
    }
}

// Every point of the cylindrical patch is below z=2s, and its top arc
// projects inside the plane z=3s: the distance is exactly s. The spatial
// edge curves deliberately disagree with that trimmed support.
TEST(BRepOffsetExt, ClearanceOfAnAnalyticTrimIgnoresItsDisplacedSpatialCurve)
{
  for (double scale : {0.1, 1.0, 10.0})
  {
    const Handle(Geom_Surface) support =
      new Geom_CylindricalSurface(gp_Ax3(gp::Origin(), gp::DZ()), 2 * scale);
    const auto source =
      BRepBuilderAPI_MakeFace(support, 0.7, 1.4, -0.5 * scale, 2 * scale, Precision::Confusion())
        .Face();
    const auto   target = BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 0, 3 * scale), gp::DZ()),
                                                  -3 * scale,
                                                  3 * scale,
                                                  -3 * scale,
                                                  3 * scale)
                            .Face();
    BRep_Builder builder;
    for (TopExp_Explorer e(source, TopAbs_EDGE); e.More(); e.Next())
    {
      const auto edge = TopoDS::Edge(e.Current());
      double     first, last;
      const auto curve = Handle(Geom_Curve)::DownCast(BRep_Tool::Curve(edge, first, last)->Copy());
      curve->Translate(gp_Vec(0, 0, 0.25 * scale));
      builder.UpdateEdge(edge, curve, 0.3 * scale);
    }
    const auto measured = BRepExt::clearance({target}, source, {}, 0, 0.9 * scale);
    EXPECT_TRUE(measured.satisfies(0.9 * scale)) << measured.issue;
    EXPECT_NEAR(measured.minimum, scale, Precision::Confusion());
  }
}

// Each vertex of an axis-aligned box moves by d on each of its three support planes.
// Check this closed form independently of the modifier's correspondence, also after a
// rigid transformation. Face/edge/vertex images must be unique and belong to the output.
TEST(BRepOffsetExt, KeptOffsetCarriesExactHistoryInEveryDimension)
{
  for (const double scale : {0.1, 1.0, 10.0})
    for (const double sign : {-1.0, 1.0})
      for (const bool turned : {false, true})
      {
        const double distance = sign * scale;
        const auto box = BRepPrimAPI_MakeBox(6 * scale, 4 * scale, 3 * scale).Shape();
        gp_Trsf transform;
        if (turned)
          transform.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 2, 3)), 0.7);
        const auto body = BRepBuilderAPI_Transform(box, transform, true).Shape();
        BRepExt::FaceSet moved;
        for (TopExp_Explorer face(body, TopAbs_FACE); face.More(); face.Next())
          moved.Add(face.Current());
        BRepExt::Lineage lineage(body);
        std::string why;
        const auto result = BRepExt::keptOffset(body, moved, distance, 1e-6 * scale, lineage, why);
        ASSERT_FALSE(result.IsNull()) << why;
        const auto history = lineage.history();
        for (const auto type : {TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX})
        {
          BRepExt::ShapeSet before, after, images;
          TopExp::MapShapes(body, type, before);
          TopExp::MapShapes(result, type, after);
          EXPECT_EQ(before.Extent(), after.Extent());
          for (const auto& source : before)
          {
            const auto& products = lineage.of(source);
            ASSERT_EQ(products.size(), 1u);
            const auto& image = products.front();
            EXPECT_EQ(image.ShapeType(), type);
            EXPECT_TRUE(after.Contains(image));
            EXPECT_TRUE(images.Add(image));
            EXPECT_FALSE(history->IsRemoved(source));
            ASSERT_EQ(history->Modified(source).Size(), 1);
            EXPECT_TRUE(history->Modified(source).First().IsSame(image));
            if (type == TopAbs_VERTEX)
            {
              const auto local = BRep_Tool::Pnt(TopoDS::Vertex(source)).Transformed(transform.Inverted());
              gp_Pnt expected(local.X() + (local.X() < 3 * scale ? -distance : distance),
                              local.Y() + (local.Y() < 2 * scale ? -distance : distance),
                              local.Z() + (local.Z() < 1.5 * scale ? -distance : distance));
              expected.Transform(transform);
              EXPECT_LE(expected.Distance(BRep_Tool::Pnt(TopoDS::Vertex(image))), 1e-6 * scale);
            }
          }
          EXPECT_EQ(images.Extent(), after.Extent());
        }
      }
}

// A top-face move changes only z=c to z=c+d. Compose its offset images with a deep copy:
// looking at Modified alone would lose every vertex and edge reported by the offset maker.
TEST(BRepOffsetExt, OffsetImagesComposeThroughADeepCopyInEveryDimension)
{
  for (const double scale : {0.1, 1.0, 10.0})
    for (const double sign : {-1.0, 1.0})
    {
      const double height = 3 * scale, distance = sign * scale;
      const auto body = BRepPrimAPI_MakeBox(6 * scale, 4 * scale, height).Shape();
      const auto top = facesWhere(body, [height](const double* l) {
        return near(l[2], height) && near(l[5], height);
      });
      ASSERT_EQ(top.size(), 1u);
      BRepOffset_MakeOffset maker;
      maker.Initialize(body, 0.0, 1e-6 * scale, BRepOffset_Skin, false, false, GeomAbs_Intersection);
      maker.SetOffsetOnFace(top[0], distance);
      maker.MakeOffsetShape();
      ASSERT_TRUE(maker.IsDone());
      BRepExt::Lineage lineage(body);
      lineage.images(maker);
      lineage.within(maker.Shape());
      BRepBuilderAPI_Copy copy(maker.Shape(), true, true);
      BRepExt::Lineage copied(maker.Shape());
      copied.through(copy);
      lineage.then(copied);
      lineage.within(copy.Shape());
      const auto history = lineage.history();
      for (const auto type : {TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX})
      {
        BRepExt::ShapeSet before, after, images;
        TopExp::MapShapes(body, type, before);
        TopExp::MapShapes(copy.Shape(), type, after);
        for (const auto& source : before)
        {
          const auto& products = lineage.of(source);
          ASSERT_EQ(products.size(), 1u);
          const auto& image = products.front();
          EXPECT_TRUE(after.Contains(image));
          EXPECT_TRUE(images.Add(image));
          ASSERT_EQ(history->Modified(source).Size(), 1);
          EXPECT_TRUE(history->Modified(source).First().IsSame(image));
          if (type == TopAbs_VERTEX)
          {
            auto expected = BRep_Tool::Pnt(TopoDS::Vertex(source));
            if (expected.Z() > height / 2) expected.SetZ(expected.Z() + distance);
            EXPECT_LE(expected.Distance(BRep_Tool::Pnt(TopoDS::Vertex(image))), 1e-6 * scale);
          }
        }
        EXPECT_EQ(images.Extent(), after.Extent());
      }
    }
}

// The endpoint/corner sets below are computed from the input, independently of
// the history reader, including rotations and both directions of the move.
TEST(BRepOffsetExt, RebuiltAndSweptCapsRetainConstructorBoundaries)
{
  for (const double scale : {0.1, 1.0, 10.0})
    for (const double sign : {-1.0, 1.0})
      for (const bool turned : {false, true})
        for (const bool rebuilt : {false, true})
        {
          const auto box = BRepPrimAPI_MakeBox(6*scale, 4*scale, 3*scale).Shape();
          const auto top = facesWhere(box, [scale](const double* l) {
            return near(l[2], 3*scale) && near(l[5], 3*scale);
          });
          ASSERT_EQ(top.size(), 1u);
          gp_Trsf transform;
          if (turned) transform.SetRotation(gp_Ax1(gp_Pnt(0,0,0), gp_Dir(1,2,3)), 0.7);
          BRepBuilderAPI_Transform move(box, transform, true);
          const auto body = move.Shape();
          const auto moved = move.ModifiedShape(top[0]);
          BRepExt::Lineage lineage(body);
          std::string why;
          std::vector<TopoDS_Shape> at;
          const auto result = rebuilt
            ? BRepExt::reconstructOffset(body, {moved}, sign*scale, 1e-6*scale, lineage)
            : BRepExt::sweepFaces(body, {moved}, sign*scale, 1e-6*scale, lineage, why, at);
          ASSERT_FALSE(result.IsNull()) << why;
          for (const auto type : {TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX})
          {
            BRepExt::ShapeSet before, after, images;
            TopExp::MapShapes(body, type, before);
            TopExp::MapShapes(result, type, after);
            ASSERT_EQ(before.Extent(), after.Extent());
            for (const auto& source : before)
            {
              const auto& products = lineage.of(source);
              ASSERT_EQ(products.size(), 1u) << "type=" << type << " rebuilt=" << rebuilt;
              const auto& image = products.front();
              ASSERT_TRUE(after.Contains(image));
              images.Add(image);
              BRepExt::ShapeSet oldVertices, newVertices;
              TopExp::MapShapes(source, TopAbs_VERTEX, oldVertices);
              TopExp::MapShapes(image, TopAbs_VERTEX, newVertices);
              ASSERT_EQ(oldVertices.Extent(), newVertices.Extent());
              for (const auto& vertex : oldVertices)
              {
                auto expected = BRep_Tool::Pnt(TopoDS::Vertex(vertex)).Transformed(transform.Inverted());
                if (expected.Z() > 1.5*scale) expected.SetZ(expected.Z()+sign*scale);
                expected.Transform(transform);
                int matches = 0;
                for (const auto& next : newVertices)
                  matches += expected.Distance(BRep_Tool::Pnt(TopoDS::Vertex(next))) <= 1e-6*scale;
                EXPECT_EQ(matches, 1);
              }
            }
            EXPECT_EQ(images.Extent(), after.Extent());
          }
        }
}

TEST(BRepOffsetExt, RecutCylinderRetainsPrimitiveBoundaryRoles)
{
  for (const double scale : {0.1, 1.0, 10.0})
    for (const bool turned : {false, true})
    {
      const auto holed = BRepAlgoAPI_Cut(
        BRepPrimAPI_MakeBox(10*scale, 10*scale, 10*scale).Shape(),
        BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(5*scale, 5*scale, 0), gp::DZ()),
                                2*scale, 10*scale).Shape()).Shape();
      TopoDS_Shape wall;
      for (TopExp_Explorer f(holed, TopAbs_FACE); f.More(); f.Next())
        if (BRepAdaptor_Surface(TopoDS::Face(f.Current())).GetType() == GeomAbs_Cylinder)
          wall = f.Current();
      ASSERT_FALSE(wall.IsNull());
      gp_Trsf transform;
      if (turned) transform.SetRotation(gp_Ax1(gp::Origin(), gp_Dir(1,2,3)), 0.7);
      BRepBuilderAPI_Transform move(holed, transform, true);
      const auto body = move.Shape();
      BRepExt::Lineage lineage(body);
      const auto result = BRepExt::reconstructOffset(body, {move.ModifiedShape(wall).Oriented(wall.Orientation())},
                                                    -scale, 1e-6*scale, lineage);
      ASSERT_FALSE(result.IsNull());
      EXPECT_NEAR(volume(result), (1000-90*M_PI)*scale*scale*scale, 1e-6*scale*scale*scale);
      for (const auto type : {TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX})
      {
        BRepExt::ShapeSet sources, products, images;
        TopExp::MapShapes(body, type, sources);
        TopExp::MapShapes(result, type, products);
        EXPECT_EQ(sources.Extent(), products.Extent());
        for (const auto& source : sources)
        {
          const auto& next = lineage.of(source);
          ASSERT_EQ(next.size(), 1u) << "type=" << type;
          ASSERT_TRUE(products.Contains(next[0]));
          images.Add(next[0]);
          if (type == TopAbs_VERTEX)
          {
            const auto before = BRep_Tool::Pnt(TopoDS::Vertex(source)).Transformed(transform.Inverted());
            const auto after = BRep_Tool::Pnt(TopoDS::Vertex(next[0])).Transformed(transform.Inverted());
            EXPECT_NEAR(before.Z(), after.Z(), 1e-6*scale);
            const auto radius = [scale](const gp_Pnt& p) { return std::hypot(p.X()-5*scale, p.Y()-5*scale); };
            // Recut chooses its own angular frame. Its seam retains its primitive
            // role and altitude, rather than promising a radial point trajectory.
            if (radius(before) < 3*scale) EXPECT_NEAR(radius(after), 3*scale, 1e-6*scale);
            else EXPECT_LE(before.Distance(after), 1e-6*scale);
          }
        }
        EXPECT_EQ(images.Extent(), products.Extent());
      }
    }
}
