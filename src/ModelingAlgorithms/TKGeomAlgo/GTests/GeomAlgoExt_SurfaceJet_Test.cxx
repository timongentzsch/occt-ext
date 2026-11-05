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
#include <GeomAlgoExt_SurfaceJet.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <gtest/gtest.h>
#include <cmath>

// Piecewise quadratic graph with C1 joins: f''=2 on the left, 4 on the
// right, h''=0 below and 2 above. A broad V span would extrapolate h''=0.
TEST(GeomAlgoExtSurfaceJet, IndependentOneSidedPolynomialDerivatives)
{
  for (double scale : {.1, 1., 10.})
  {
    NCollection_Array1<double> knots(1, 3);
    knots(1) = 0;
    knots(2) = .5;
    knots(3) = 1;
    NCollection_Array1<int> mult(1, 3);
    mult(1) = 3;
    mult(2) = 1;
    mult(3) = 3;
    NCollection_Array2<gp_Pnt> poles(1, 4, 1, 4);
    const double               coordinate[] = {0, .25, .75, 1};
    const double               f[] = {0, 0, .5, 1.25}, h[] = {0, 0, 0, .25};
    for (int i = 1; i <= 4; ++i)
      for (int j = 1; j <= 4; ++j)
        poles(i, j) = gp_Pnt(scale * coordinate[i - 1],
                             scale * coordinate[j - 1],
                             scale * (f[i - 1] + h[j - 1]));
    const Handle(Geom_BSplineSurface) surface =
      new Geom_BSplineSurface(poles, knots, knots, mult, mult, 2, 2);
    for (int left : {1, 2})
    {
      BRepExt::SurfaceJet jet{"test", {}, {}, {}, {}, {}, {}, {}};
      BRepExt::surfaceD2OnUSpan(surface, .5, .75, left, left + 1, jet);
      jet.orient(false);
      EXPECT_NEAR(jet.u.X(), scale, 1e-12 * scale);
      EXPECT_NEAR(jet.u.Z(), scale, 1e-12 * scale);
      EXPECT_NEAR(jet.v.Y(), scale, 1e-12 * scale);
      EXPECT_NEAR(jet.v.Z(), .5 * scale, 1e-12 * scale);
      EXPECT_NEAR(jet.uu.Z(), (left == 1 ? 2 : 4) * scale, 1e-12 * scale);
      EXPECT_NEAR(jet.vv.Z(), 2 * scale, 1e-12 * scale);
      EXPECT_NEAR(jet.uv.Magnitude(), 0, 1e-12 * scale);
      // II(Su,Su) = n dot Suu for the independent graph normal (-1,-.5,1)/1.5.
      EXPECT_NEAR(jet.curvature(jet.u, jet.u), (left == 1 ? 2 : 4) * scale / 1.5, 1e-11 * scale);
      jet.orient(true);
      EXPECT_NEAR(jet.curvature(jet.v, jet.v), -2 * scale / 1.5, 1e-11 * scale);
    }
  }
}

#include <GeomAlgoExt_PlanarCurve.hxx>
#include <Geom_BezierCurve.hxx>
#include <gp_Pln.hxx>

// z(t)=70 t^4(1-t)^4 has maximum 70/256. Its initial tenth is below .2.
// This catches checking only endpoints, using control poles as exact extrema,
// or bounding the untrimmed curve when a trimmed interval was requested.
TEST(GeomAlgoExtSurfaceJet, TrimmedPolynomialPlaneSlab)
{
  NCollection_Array1<gp_Pnt> poles(1, 9);
  for (int i = 1; i <= 9; ++i)
    poles(i) = gp_Pnt((i - 1) / 8., 0, i == 5 ? 1 : 0);
  const Handle(Geom_BezierCurve) curve = new Geom_BezierCurve(poles);
  const gp_Pln                   plane(gp::XOY());
  EXPECT_NO_THROW(BRepExt::requirePlanarCurve("test", curve, 0, 1, plane, .3));
  EXPECT_THROW(BRepExt::requirePlanarCurve("test", curve, 0, 1, plane, .2), BRepExt::Failure);
  EXPECT_NO_THROW(BRepExt::requirePlanarCurve("test", curve, 0, .1, plane, .2));
}

// At the origin z=k*x*y has principal curvatures +/-k. Rotating the
// observation frame cannot change either the absolute or relative mismatch.
TEST(GeomAlgoExtSurfaceJet, MixedCurvatureRelativeDeviation)
{
  for (double k : {1e-18, .01, 1., 100.})
    for (double angle : {0., .13, .7, 1.2})
    {
      BRepExt::SurfaceJet saddle{"test", {}, {1, 0, 0}, {0, 1, 0}, {}, {}, {0, 0, k}, {}};
      saddle.orient(false);
      auto doubled = saddle;
      doubled.uv *= 2;
      auto plane = saddle;
      plane.uv.SetCoord(0, 0, 0);
      const gp_Vec t(std::cos(angle), std::sin(angle), 0);
      const auto mismatch = BRepExt::bendDifference(saddle, doubled, t);
      EXPECT_NEAR(mismatch.first, k, k * 1e-12);
      EXPECT_NEAR(mismatch.second, .5, 1e-12);
      EXPECT_NEAR(BRepExt::bendDifference(saddle, plane, t).second, 1., 1e-12);
      EXPECT_EQ(BRepExt::bendDifference(plane, plane, t).second, 0.);
    }
}
