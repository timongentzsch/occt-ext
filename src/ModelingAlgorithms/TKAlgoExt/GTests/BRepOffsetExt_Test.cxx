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

// Push-pull, draft, defeature and hollow of BRepOffsetExt on bodies whose result has a closed form
// of its input: a box, a cylinder, a box with a boss.

#include "BRepExtTest_Support.hxx"

#include <BRepOffsetExt_Kept.hxx>
#include <BRepOffsetExt_Focal.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <Geom_BezierSurface.hxx>
#include <GeomConvert.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_Failure.hxx>

#include <BRepAdaptor_Surface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <NCollection_Array2.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
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
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
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
  const auto         offset = [&](double distance, int join) {
    const auto result = BRepExt::offsetSolid(session(),
                                             keep(box),
                                             distance,
                                             1e-6,
                                             join,
                                                     {},
                                             BRepExt::failureHashBound);
    return BRepExt::get(session(), result.resultId);
  };
  const TopoDS_Shape rounded = offset(d, 0), sharp = offset(d, 2), inward = offset(-d, 2);
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

// Every face of a box with one edge rounded by r moved by d, as the reference modeller's face offset keeps it:
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
                                                 0,
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
