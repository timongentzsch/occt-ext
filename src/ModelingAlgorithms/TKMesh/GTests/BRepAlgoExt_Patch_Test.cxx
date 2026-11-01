// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepAlgoExt_PatchFlow.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <gp_Pln.hxx>
#include <gp_Sphere.hxx>
#include <Precision.hxx>
#include <gtest/gtest.h>
#include <cmath>
// A flattening of a Euclidean rectangle must preserve all pair distances.
// Wrong cotangent signs, boundary turning or UV lookup strides fail this oracle.
TEST(BRepAlgoExtPatch, PlanarChartIsometry) {
 for(double scale:{.1,1.,10.}) {
  const auto face=BRepBuilderAPI_MakeFace(gp_Pln(gp::XOY()),0,2*scale,0,scale).Face();
  const auto chart=BRepExt::conformalChart(face,.2*scale,.01*scale);
  ASSERT_FALSE(chart.triangles.empty());
  for(int i=0;i<=8;++i)for(int j=0;j<=8;++j) {
   const gp_XY a(2*scale*i/8.,scale*j/8.);
   for(const gp_XY b:{gp_XY(0,0),gp_XY(2*scale,scale),gp_XY(scale,.5*scale)})
    EXPECT_NEAR((chart.map(a)-chart.map(b)).Modulus(),(a-b).Modulus(),Precision::Confusion()*scale);
  }
 }
}
// Sphere principal curvature has magnitude 1/r, independent of chart or sample grid.
// Sign-insensitive magnitude detects wrong surface/normal order; ratio gate tests policy.
