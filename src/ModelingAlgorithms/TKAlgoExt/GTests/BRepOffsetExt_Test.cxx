// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Push-pull, draft, defeature and hollow of BRepOffsetExt on bodies whose result has a closed form
// of its input: a box, a cylinder, a box with a boss.

#include "BRepExtTest_Support.hxx"

#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_Failure.hxx>

#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>

#include <gtest/gtest.h>

using namespace BRepExtTest;

namespace {
const double a = 60, b = 40, c = 30;

uint32_t keep(const TopoDS_Shape& shape) { return BRepExt::store(session(), shape); }

std::vector<uint32_t> keep(const std::vector<TopoDS_Face>& faces) {
  std::vector<uint32_t> ids;
  for (const TopoDS_Face& face : faces) ids.push_back(keep(face));
  return ids;
}
} // namespace

// Pulling the top face of a box out by d lengthens the box to c + d and leaves every side face
// planar and square to it: a x b x (c + d). Pushing it in by d shortens it likewise. Neighbours
// extended or trimmed wrongly change the volume by a multiple of d.
TEST(BRepOffsetExt, PushPullOfTheTopFaceMovesItAlongItsNormal) {
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(a, b, c).Shape();
  const auto top = facesWhere(box, [](const double* l) { return near(l[2], c) && near(l[5], c); });
  ASSERT_EQ(top.size(), 1u);
  for (const double d : {7.0, -7.0}) {
    const auto result = BRepExt::offsetFaces(session(), keep(box), keep(top), d, 1e-6, {}, BRepExt::failureHashBound);
    const TopoDS_Shape& made = BRepExt::get(session(), result.resultId);
    EXPECT_TRUE(BRepExt::invalidity(made).empty());
    EXPECT_NEAR(volume(made), a * b * (c + d), 1e-9 * a * b * c) << "d=" << d;
  }
}

// Drafting the four sides of a box by alpha about the plane z = 0, pulling along +z, tilts each
// side inward by alpha as z grows, so a section at height z is (a - 2 z tan alpha)(b - 2 z tan alpha)
// and the volume is its integral over [0, c]: abc - (a + b) c^2 tan alpha + 4/3 c^3 tan^2 alpha.
// A draft about the wrong plane, in the wrong direction or by the wrong angle changes it at first
// order in alpha.
TEST(BRepOffsetExt, DraftOfTheSidesOfABoxIsAFrustumOfClosedFormVolume) {
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(a, b, c).Shape();
  const auto sides = facesWhere(box, [](const double* l) { return l[5] - l[2] > c - 1e-6; });
  ASSERT_EQ(sides.size(), 4u);
  for (const double degrees : {2.0, 5.0}) {
    const double alpha = degrees * M_PI / 180, t = std::tan(alpha);
    const auto result = BRepExt::draftFaces(session(), keep(box), keep(sides), alpha, {0, 0, 1}, {0, 0, 0}, {0, 0, 1},
                                            {}, BRepExt::failureHashBound);
    const TopoDS_Shape& made = BRepExt::get(session(), result.resultId);
    EXPECT_TRUE(BRepExt::invalidity(made).empty());
    EXPECT_NEAR(volume(made), a * b * c - (a + b) * c * c * t + 4. / 3 * c * c * c * t * t, 1e-9 * a * b * c)
        << degrees << " degrees";
  }
}

// A boss fused onto a box and then defeatured leaves the box: the boss's wall and top are removed
// and the neighbours closed over it. A defeature that leaves a stub or cuts the box changes the
// volume by that of the boss.
TEST(BRepOffsetExt, DefeatureOfABossLeavesTheBox) {
  const double radius = 8, height = 6;
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(a, b, c).Shape();
  const TopoDS_Shape boss =
      BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(a / 2, b / 2, c), gp_Dir(0, 0, 1)), radius, height).Shape();
  const TopoDS_Shape fused = BRepAlgoAPI_Fuse(box, boss).Shape();
  const auto features = facesWhere(fused, [](const double* l) { return l[5] > c + 1e-6 && l[2] >= c - 1e-6; });
  ASSERT_EQ(features.size(), 2u);
  const auto result = BRepExt::defeature(session(), keep(fused), keep(features), {}, BRepExt::failureHashBound);
  const TopoDS_Shape& made = BRepExt::get(session(), result.resultId);
  EXPECT_TRUE(BRepExt::invalidity(made).empty());
  EXPECT_NEAR(volume(made), a * b * c, 1e-9 * a * b * c);
}

// Hollowing a cylinder of radius R and height H to a wall t leaves the cylinder of radius R - t
// and height H - 2t missing: pi R^2 H - pi (R - t)^2 (H - 2t). The curved wall must be offset
// along its normals and the cavity closed by planar caps at t from the ends.
TEST(BRepOffsetExt, HollowOfACylinderLeavesTheClosedFormWall) {
  const double radius = 15, height = 40, t = 2.5;
  const TopoDS_Shape cylinder = BRepPrimAPI_MakeCylinder(radius, height).Shape();
  const auto result = BRepExt::hollow(session(), keep(cylinder), t, 1e-4, {}, BRepExt::failureHashBound);
  const TopoDS_Shape& made = BRepExt::get(session(), result.resultId);
  EXPECT_TRUE(BRepExt::invalidity(made).empty());
  const double expected = M_PI * radius * radius * height - M_PI * (radius - t) * (radius - t) * (height - 2 * t);
  EXPECT_NEAR(volume(made), expected, 1e-6 * expected);
}

// Every face of a solid moved out by d: with round joins the result is the parallel body, whose
// volume is Steiner's polynomial abc + 2d(ab + bc + ca) + pi d^2 (a + b + c) + 4/3 pi d^3; with sharp
// joins the box simply grows to (a + 2d)(b + 2d)(c + 2d); moved in by d it shrinks to
// (a - 2d)(b - 2d)(c - 2d). The three differ by far more than the tolerance in the corners.
TEST(BRepOffsetExt, OffsetOfABoxIsSteinersParallelBodyOrTheGrownBox) {
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(a, b, c).Shape();
  const double d = 4;
  const auto offset = [&](double distance, int join) {
    const auto result = BRepExt::offsetSolid(session(), keep(box), distance, 1e-6, join, {}, BRepExt::failureHashBound);
    return BRepExt::get(session(), result.resultId);
  };
  const TopoDS_Shape rounded = offset(d, 0), sharp = offset(d, 2), inward = offset(-d, 2);
  EXPECT_TRUE(BRepExt::invalidity(rounded).empty());
  EXPECT_TRUE(BRepExt::invalidity(sharp).empty());
  EXPECT_TRUE(BRepExt::invalidity(inward).empty());
  EXPECT_NEAR(volume(rounded),
              a * b * c + 2 * d * (a * b + b * c + c * a) + M_PI * d * d * (a + b + c) + 4. / 3 * M_PI * d * d * d,
              1e-6 * a * b * c);
  EXPECT_NEAR(volume(sharp), (a + 2 * d) * (b + 2 * d) * (c + 2 * d), 1e-9 * a * b * c);
  EXPECT_NEAR(volume(inward), (a - 2 * d) * (b - 2 * d) * (c - 2 * d), 1e-9 * a * b * c);
}
