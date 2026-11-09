// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// The checked Boolean of BRepBooleanExt on bodies whose result has a closed form, including curved
// intersection curves: crossed cylinders (the Steinmetz solid) and a box pierced by a cylinder.

#include "BRepExtTest_Support.hxx"

#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_Failure.hxx>

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <GC_MakeCircle.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>

#include <gtest/gtest.h>

using namespace BRepExtTest;

namespace
{

TopoDS_Shape combine(BRepExt::Select     operation,
                     const TopoDS_Shape& source,
                     const TopoDS_Shape& tool)
{
  const uint32_t a = BRepExt::store(session(), source), b = BRepExt::store(session(), tool);
  const auto     result = BRepExt::booleanChecked(session(),
                                              operation,
                                                  {a},
                                                  {b},
                                              BRepExt::BooleanGlue::Off,
                                              0,
                                              0,
                                                  {},
                                              BRepExt::failureHashBound);
  return BRepExt::get(session(), result.resultId);
}

} // namespace

// Two cylinders of radius r with perpendicular, intersecting axes meet in the Steinmetz solid, of
// volume 16 r^3 / 3 (Archimedes). Its edges are two ellipses where the surfaces cross; a Boolean
// that misplaces them or drops a lens changes the volume by far more than the tolerance.
TEST(BRepBooleanExt, IntersectionOfCrossedCylindersIsTheSteinmetzSolid)
{
  const double       r = 5, length = 4 * r;
  const TopoDS_Shape along =
    BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -length / 2), gp_Dir(0, 0, 1)), r, length).Shape();
  const TopoDS_Shape across =
    BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-length / 2, 0, 0), gp_Dir(1, 0, 0)), r, length).Shape();
  const TopoDS_Shape made = combine(BRepExt::Select::Intersection, along, across);
  EXPECT_TRUE(BRepExt::invalidity(made).empty());
  EXPECT_NEAR(volume(made), 16 * r * r * r / 3, 1e-6 * r * r * r);
}

// The union of the same two cylinders has the volume of both minus the Steinmetz solid counted
// twice: 2 pi r^2 L - 16 r^3 / 3, by inclusion-exclusion.
TEST(BRepBooleanExt, UnionOfCrossedCylindersFollowsInclusionExclusion)
{
  const double       r = 5, length = 4 * r;
  const TopoDS_Shape along =
    BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, -length / 2), gp_Dir(0, 0, 1)), r, length).Shape();
  const TopoDS_Shape across =
    BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-length / 2, 0, 0), gp_Dir(1, 0, 0)), r, length).Shape();
  const TopoDS_Shape made = combine(BRepExt::Select::Union, along, across);
  EXPECT_TRUE(BRepExt::invalidity(made).empty());
  EXPECT_NEAR(volume(made), 2 * M_PI * r * r * length - 16 * r * r * r / 3, 1e-6 * r * r * length);
}

// A cylinder through the whole height of a box removes pi r^2 c: the hole's wall must run exactly
// between the two faces it pierces.
TEST(BRepBooleanExt, DifferenceOfABoxAndAThroughCylinderRemovesTheBore)
{
  const double       a = 60, b = 40, c = 30, r = 7;
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(a, b, c).Shape();
  const TopoDS_Shape bore =
    BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(a / 2, b / 2, -5), gp_Dir(0, 0, 1)), r, c + 10).Shape();
  const TopoDS_Shape made = combine(BRepExt::Select::Difference, box, bore);
  EXPECT_TRUE(BRepExt::invalidity(made).empty());
  EXPECT_NEAR(volume(made), a * b * c - M_PI * r * r * c, 1e-9 * a * b * c);
}

// Two boxes that overlap in a box: inclusion-exclusion again, with the overlap found from the
// coordinates alone.
TEST(BRepBooleanExt, UnionOfOverlappingBoxesRemovesTheOverlapOnce)
{
  const TopoDS_Shape first  = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 10, 10, 10).Shape();
  const TopoDS_Shape second = BRepPrimAPI_MakeBox(gp_Pnt(6, 4, 2), 10, 10, 10).Shape();
  const TopoDS_Shape made   = combine(BRepExt::Select::Union, first, second);
  EXPECT_TRUE(BRepExt::invalidity(made).empty());
  EXPECT_NEAR(volume(made), 2000 - 4 * 6 * 8, 1e-9 * 2000);
}

// A circle of radius 4 above the middle of the top face of a 20 x 20 x 10 box, projected down onto
// it and imprinted, splits the face into the disc of area 16 pi and the ring around it; the solid
// is the same body (volume 4000, area unchanged) with one more face. A projection along the
// wrong direction, or a body re-made rather than imprinted, changes the volume or the areas of the
// two faces at first order.
TEST(BRepBooleanExt, ImprintOfACircleSplitsTheFaceIntoDiscAndRing)
{
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(gp_Pnt(-10, -10, 0), 20, 20, 10).Shape();
  const TopoDS_Wire  circle =
    BRepBuilderAPI_MakeWire(
      BRepBuilderAPI_MakeEdge(GC_MakeCircle(gp_Ax2(gp_Pnt(0, 0, 15), gp_Dir(0, 0, 1)), 4.0).Value())
        .Edge())
      .Wire();
  const auto          result = BRepExt::imprintCurves(session(),
                                             BRepExt::store(session(), box),
                                                      {BRepExt::store(session(), circle)},
                                             BRepExt::ImprintMethod::Direction,
                                                      {0, 0, -1},
                                             false,
                                             false,
                                             BRepExt::ImprintCompletion::None,
                                                      {},
                                             BRepExt::failureHashBound);
  const TopoDS_Shape& made   = BRepExt::get(session(), result.resultId);
  EXPECT_NEAR(volume(made), 4000, 1e-6 * 4000);
  EXPECT_NEAR(area(made), area(box), 1e-6 * area(box));
  const auto top =
    facesWhere(made, [](const double* l) { return near(l[2], 10) && near(l[5], 10); });
  ASSERT_EQ(top.size(), 2u);
  const double small = std::min(area(top[0]), area(top[1])),
               large = std::max(area(top[0]), area(top[1]));
  EXPECT_NEAR(small, 16 * M_PI, 1e-6 * 16 * M_PI);
  EXPECT_NEAR(large, 400 - 16 * M_PI, 1e-6 * 400);
}
