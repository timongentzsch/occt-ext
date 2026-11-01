// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <BRepAlgoExt_TrimDomain.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepAlgoExt_ReadBack.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_Circ.hxx>
#include <gp_Pln.hxx>
#include <gtest/gtest.h>
#include <cmath>

// Exact moments of a rectangle with an eccentric circular hole. The independent
// signed subtraction catches ignoring inner wires, wrong winding and UV Jacobian.
TEST(BRepAlgoExtShared, EccentricHoleMoments) {
  for(double s : {.1,1.,10.}) {
    BRepBuilderAPI_MakeFace maker(gp_Pln(gp_Pnt(0,0,2*s),gp::DZ()),0,4*s,0,3*s);
    const auto hole=BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(gp_Circ(gp_Ax2(gp_Pnt(s,s,2*s),gp::DZ()),.4*s)).Edge()).Wire();
    maker.Add(TopoDS::Wire(hole.Reversed()));
    const gp_Pnt origin(0,0,0);
    const auto result=BRepExt::trimmedProperties(maker.Face(),&origin);
    const double disk=M_PI*.16*s*s,area=12*s*s-disk;
    ASSERT_TRUE(result.converged);
    EXPECT_NEAR(result.area,area,1e-9*s*s);
    EXPECT_NEAR(result.center.X(),(24*s*s*s-disk*s)/area,1e-9*s);
    EXPECT_NEAR(result.center.Y(),(18*s*s*s-disk*s)/area,1e-9*s);
    EXPECT_NEAR(result.center.Z(),2*s,1e-9*s);
    EXPECT_NEAR(result.volume,area*2*s/3,1e-9*s*s*s);
  }
}
// Stored tolerance is a maximum over topology, without a document floor.
// A nested vertex, not the root container or its first child, holds the maximum.
TEST(BRepAlgoExtShared, NestedToleranceMaximum) {
  BRep_Builder b;TopoDS_Compound outer,inner;b.MakeCompound(outer);b.MakeCompound(inner);
  for(double t : {1e-5,3e-4,2e-5}) {
    TopoDS_Vertex v;b.MakeVertex(v,gp_Pnt(t,0,0),t);b.Add(inner,v);
  }
  b.Add(outer,inner);
  EXPECT_DOUBLE_EQ(BRepExt::tolerance::own(outer),0);
  EXPECT_DOUBLE_EQ(BRepExt::tolerance::largest(outer),3e-4);
}
// The reader normalizes coordinates again. These non-axis-aligned direction
// families catch returning an unstable input or perturbing the direction grossly.
TEST(BRepAlgoExtShared, DirectionReaderFixedPoints) {
  for(int i=1;i<=128;++i) {
    const gp_Dir input(std::sin(i*.37),std::cos(i*.19),std::sin(i*.11)+.4);
    const auto settled=BRepExt::readBack(input);gp_Dir read;
    read.SetCoord(settled.X(),settled.Y(),settled.Z());
    EXPECT_TRUE(BRepExt::same(settled,read));
    EXPECT_LT(gp_Vec(input).Crossed(gp_Vec(settled)).Magnitude(),1e-12);
  }
}

#include <BRepAlgoExt_Properties.hxx>
#include <gp_Sphere.hxx>
// Integral r^2 cos(v) du dv on a sphere rectangle; probes both closed-form
// and adaptive paths without deriving expected values from another OCCT measure.
TEST(BRepAlgoExtShared, SphereRectangleArea) {
  for(double r : {.1,1.,10.}) {
    const auto face=BRepBuilderAPI_MakeFace(gp_Sphere(gp_Ax3(),r),.2,1.3,-.4,.7).Face();
    const double expected=r*r*1.1*(std::sin(.7)-std::sin(-.4));
    EXPECT_NEAR(BRepExt::faceArea(face),expected,1e-10*r*r);
    bool converged=true;
    EXPECT_NEAR(BRepExt::faceMeasure(face,true,gp_Pnt(),&converged),expected,1e-8*r*r);
    EXPECT_TRUE(converged);
  }
}
