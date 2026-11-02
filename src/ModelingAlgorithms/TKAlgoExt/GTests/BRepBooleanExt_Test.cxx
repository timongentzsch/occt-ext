// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// The checked Boolean of BRepBooleanExt on bodies whose result has a closed form, including curved
// intersection curves: crossed cylinders (the Steinmetz solid) and a box pierced by a cylinder.

#include "BRepExtTest_Support.hxx"

#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_Failure.hxx>

#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>

#include <gtest/gtest.h>

using namespace BRepExtTest;

namespace
{

enum Operation
{
  Union        = 0,
  Difference   = 1,
  Intersection = 2
};

TopoDS_Shape combine(Operation operation, const TopoDS_Shape& source, const TopoDS_Shape& tool)
{
  const uint32_t a = BRepExt::store(session(), source), b = BRepExt::store(session(), tool);
  const auto     result =
    BRepExt::booleanChecked(session(), operation, {a}, {b}, 0, 0, 0, {}, BRepExt::failureHashBound);
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
  const TopoDS_Shape made = combine(Intersection, along, across);
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
  const TopoDS_Shape made = combine(Union, along, across);
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
  const TopoDS_Shape made = combine(Difference, box, bore);
  EXPECT_TRUE(BRepExt::invalidity(made).empty());
  EXPECT_NEAR(volume(made), a * b * c - M_PI * r * r * c, 1e-9 * a * b * c);
}

// Two boxes that overlap in a box: inclusion-exclusion again, with the overlap found from the
// coordinates alone.
TEST(BRepBooleanExt, UnionOfOverlappingBoxesRemovesTheOverlapOnce)
{
  const TopoDS_Shape first  = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), 10, 10, 10).Shape();
  const TopoDS_Shape second = BRepPrimAPI_MakeBox(gp_Pnt(6, 4, 2), 10, 10, 10).Shape();
  const TopoDS_Shape made   = combine(Union, first, second);
  EXPECT_TRUE(BRepExt::invalidity(made).empty());
  EXPECT_NEAR(volume(made), 2000 - 4 * 6 * 8, 1e-9 * 2000);
}
