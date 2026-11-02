// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Chamfers and variable-radius fillets of BRepFilletExt on a straight box edge, where the removed
// material is a prism and has a closed form: its cross-section area times the edge length L.

#include "BRepExtTest_Support.hxx"

#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_Failure.hxx>

#include <BRepAdaptor_Curve.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <TopExp.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <NCollection_List.hxx>
#include <TopTools_ShapeMapHasher.hxx>

#include <gtest/gtest.h>

using namespace BRepExtTest;

namespace
{

// The edge of the box a x b x c that runs along z at (a, b).
TopoDS_Edge uprightAtFarCorner(const TopoDS_Shape& box, double a, double b)
{
  const auto at = [&](const double* l) {
    return near(l[0], a) && near(l[3], a) && near(l[1], b) && near(l[4], b);
  };
  const auto edges = edgesWhere(box, at);
  return edges.empty() ? TopoDS_Edge() : edges.front();
}

} // namespace

// An equal-distance chamfer of distance d removes the triangular prism of area d^2 / 2: a wrong
// setback on either face, or a bevel that is not planar, changes the volume by far more than the
// tolerance.
TEST(BRepFilletExt, EqualChamferRemovesATriangularPrism)
{
  const double       a = 6, b = 4, c = 3;
  const TopoDS_Shape box  = BRepPrimAPI_MakeBox(a, b, c).Shape();
  const TopoDS_Edge  edge = uprightAtFarCorner(box, a, b);
  ASSERT_FALSE(edge.IsNull());
  for (const double d : {0.3, 1.0})
  {
    const uint32_t body = BRepExt::store(session(), box), picked = BRepExt::store(session(), edge);
    const auto     result =
      BRepExt::chamfer(session(), body, {picked}, d, d, 0, {}, BRepExt::failureHashBound, 0, false);
    const TopoDS_Shape& made = BRepExt::get(session(), result.resultId);
    EXPECT_TRUE(BRepExt::invalidity(made).empty());
    EXPECT_NEAR(volume(made), a * b * c - 0.5 * d * d * c, 1e-9 * a * b * c) << "d=" << d;
  }
}

// An asymmetric chamfer takes d1 along the reference face and d2 along the other: the prism
// removed has area d1 d2 / 2, whichever of the two faces is the reference.
TEST(BRepFilletExt, AsymmetricChamferRemovesTheRightTriangle)
{
  const double       a = 6, b = 4, c = 3, d1 = 1.5, d2 = 0.5;
  const TopoDS_Shape box  = BRepPrimAPI_MakeBox(a, b, c).Shape();
  const TopoDS_Edge  edge = uprightAtFarCorner(box, a, b);
  ASSERT_FALSE(edge.IsNull());
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    faces;
  TopExp::MapShapesAndAncestors(box, TopAbs_EDGE, TopAbs_FACE, faces);
  const auto& adjacent = faces.FindFromKey(edge);
  ASSERT_EQ(adjacent.Extent(), 2);
  for (const TopoDS_Shape& reference : adjacent)
  {
    const uint32_t body = BRepExt::store(session(), box), picked = BRepExt::store(session(), edge),
                   ref         = BRepExt::store(session(), reference);
    const auto          result = BRepExt::chamfer(session(),
                                         body,
                                                  {picked},
                                         d1,
                                         d2,
                                         ref,
                                                  {},
                                         BRepExt::failureHashBound,
                                         0,
                                         false);
    const TopoDS_Shape& made   = BRepExt::get(session(), result.resultId);
    EXPECT_TRUE(BRepExt::invalidity(made).empty());
    EXPECT_NEAR(volume(made), a * b * c - 0.5 * d1 * d2 * c, 1e-9 * a * b * c);
  }
}

// A variable fillet whose law is constant is the circular fillet: the quarter-disc complement of
// area r^2 (1 - pi/4). It exercises the variable-radius path (law points, sections along the
// spine) on a case whose answer does not depend on it.
TEST(BRepFilletExt, ConstantLawVariableFilletIsTheCircularFillet)
{
  const double       a = 6, b = 4, c = 3, r = 0.5;
  const TopoDS_Shape box  = BRepPrimAPI_MakeBox(a, b, c).Shape();
  const TopoDS_Edge  edge = uprightAtFarCorner(box, a, b);
  ASSERT_FALSE(edge.IsNull());
  const uint32_t body = BRepExt::store(session(), box), picked = BRepExt::store(session(), edge);
  const auto     result    = BRepExt::filletVariable(session(),
                                              body,
                                                     {picked},
                                                     {3},
                                                     {0, r, 0.5, r, 1, r},
                                                     {},
                                              BRepExt::failureHashBound);
  const TopoDS_Shape& made = BRepExt::get(session(), result.resultId);
  EXPECT_TRUE(BRepExt::invalidity(made).empty());
  EXPECT_NEAR(volume(made), a * b * c - r * r * (1 - M_PI / 4) * c, 1e-6 * r * r * c);
}

// A law that grows along the edge removes more material than its smaller end alone and less than
// its larger end alone, since every cross-section of the blend lies between the two circular
// ones: the volume lies between the two constant-radius closed forms. A law ignored, or a radius
// taken from the wrong parameter range, leaves it outside.
TEST(BRepFilletExt, GrowingLawRemovesBetweenTheTwoConstantRadiusVolumes)
{
  const double       a = 6, b = 4, c = 3, r0 = 0.3, r1 = 0.9;
  const TopoDS_Shape box  = BRepPrimAPI_MakeBox(a, b, c).Shape();
  const TopoDS_Edge  edge = uprightAtFarCorner(box, a, b);
  ASSERT_FALSE(edge.IsNull());
  const uint32_t body = BRepExt::store(session(), box), picked = BRepExt::store(session(), edge);
  const auto     result    = BRepExt::filletVariable(session(),
                                              body,
                                                     {picked},
                                                     {2},
                                                     {0, r0, 1, r1},
                                                     {},
                                              BRepExt::failureHashBound);
  const TopoDS_Shape& made = BRepExt::get(session(), result.resultId);
  EXPECT_TRUE(BRepExt::invalidity(made).empty());
  const double removed = a * b * c - volume(made), k = (1 - M_PI / 4) * c;
  EXPECT_GT(removed, k * r0 * r0);
  EXPECT_LT(removed, k * r1 * r1);
}

// A chamfer, like every operation, refuses what it cannot do with a typed failure instead of
// returning a wrong body: an asymmetric chamfer without a reference face has no side to take d1 on.
TEST(BRepFilletExt, AsymmetricChamferWithoutReferenceIsRefused)
{
  const TopoDS_Shape box  = BRepPrimAPI_MakeBox(6, 4, 3).Shape();
  const TopoDS_Edge  edge = uprightAtFarCorner(box, 6, 4);
  ASSERT_FALSE(edge.IsNull());
  const uint32_t body = BRepExt::store(session(), box), picked = BRepExt::store(session(), edge);
  EXPECT_THROW(BRepExt::chamfer(session(),
                                body,
                                {picked},
                                1.5,
                                0.5,
                                0,
                                {},
                                BRepExt::failureHashBound,
                                0,
                                false),
               std::exception);
}
