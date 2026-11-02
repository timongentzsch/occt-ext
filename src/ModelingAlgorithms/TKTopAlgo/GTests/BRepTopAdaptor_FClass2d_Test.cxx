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

#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepTools.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Pnt2d.hxx>

#include <gtest/gtest.h>

#include <cmath>

// The wall of a half cylinder is the face u in [0, pi], v in [0, 4] of a surface that is periodic
// in u. With the periodic reframing on, a point is classified by its place within the period: u =
// pi/2 plus any whole number of periods is on the face, u = 3 pi/2 plus any is not, and v decides
// on its own, since the surface is not periodic in v.
TEST(BRepTopAdaptor_FClass2dTest, PeriodicReframingClassifiesByThePlaceWithinThePeriod)
{
  const TopoDS_Shape aHalf = BRepPrimAPI_MakeCylinder(3.0, 4.0, M_PI).Shape();
  TopoDS_Face        aWall;
  for (TopExp_Explorer anIt(aHalf, TopAbs_FACE); anIt.More(); anIt.Next())
  {
    const TopoDS_Face aFace = TopoDS::Face(anIt.Current());
    double            aU0, aU1, aV0, aV1;
    BRepTools::UVBounds(aFace, aU0, aU1, aV0, aV1);
    if (std::abs(aU1 - aU0 - M_PI) < 1e-9 && std::abs(aV1 - aV0 - 4.0) < 1e-9)
    {
      aWall = aFace;
    }
  }
  ASSERT_FALSE(aWall.IsNull());
  const BRepTopAdaptor_FClass2d aClass(aWall, 1e-7);
  for (int aPeriods = -3; aPeriods <= 3; ++aPeriods)
  {
    const double aShift = 2.0 * M_PI * aPeriods;
    EXPECT_EQ(aClass.Perform(gp_Pnt2d(M_PI / 2 + aShift, 2.0), true), TopAbs_IN) << aPeriods;
    EXPECT_EQ(aClass.Perform(gp_Pnt2d(3 * M_PI / 2 + aShift, 2.0), true), TopAbs_OUT) << aPeriods;
    EXPECT_EQ(aClass.Perform(gp_Pnt2d(M_PI / 2 + aShift, 5.0), true), TopAbs_OUT) << aPeriods;
  }
}
