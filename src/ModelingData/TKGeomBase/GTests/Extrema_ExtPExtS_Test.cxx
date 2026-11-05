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

#include <Extrema_ExtPExtS.hxx>
#include <Extrema_ExtPC.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomAdaptor_SurfaceOfLinearExtrusion.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_SurfaceOfLinearExtrusion.hxx>
#include <gtest/gtest.h>

TEST(Extrema_ExtPExtS, TransverseSplineReturnsAllCurveExtrema)
{
  for (const int spans : {6, 8, 12})
  {
    NCollection_Array1<gp_Pnt> poles(1, 2 * spans + 1);
    NCollection_Array1<double> knots(1, spans + 1);
    NCollection_Array1<int> mults(1, spans + 1);
    for (int i = 0; i <= 2 * spans; ++i)
      poles(i + 1) = gp_Pnt(i, i % 2 ? 2 : 0, 0);
    for (int i = 0; i <= spans; ++i)
    {
      knots(i + 1) = i;
      mults(i + 1) = i == 0 || i == spans ? 3 : 2;
    }
    Handle(Geom_BSplineCurve) curve = new Geom_BSplineCurve(poles, knots, mults, 2);
    const gp_Pnt point(spans, 100, 2);
    GeomAdaptor_Curve adapted(curve);
    Extrema_ExtPC curveExt(gp_Pnt(spans, 100, 0), adapted, 0, spans, 1e-9);
    ASSERT_TRUE(curveExt.IsDone());
    ASSERT_GT(curveExt.NbExt(), 4);
    Handle(Geom_SurfaceOfLinearExtrusion) surface =
      new Geom_SurfaceOfLinearExtrusion(curve, gp_Dir(0, 0, 1));
    Handle(GeomAdaptor_SurfaceOfLinearExtrusion) adaptor =
      new GeomAdaptor_SurfaceOfLinearExtrusion(new GeomAdaptor_Curve(curve), gp_Dir(0, 0, 1));
    Extrema_ExtPExtS extrema(point, adaptor, 0, spans, -10, 10, 1e-9, 1e-9);
    ASSERT_TRUE(extrema.IsDone());
    EXPECT_EQ(extrema.NbExt(), curveExt.NbExt());
    for (int i = 1; i <= extrema.NbExt(); ++i)
      EXPECT_NEAR(extrema.Point(i).Value().Z(), 2, 1e-12);
  }
}
