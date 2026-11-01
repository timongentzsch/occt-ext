// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_PatchFit.hxx>
#include <gtest/gtest.h>
#include <GeomAlgoExt_PatchRim.hxx>
#include <GeomAlgoExt_PatchCompact.hxx>
#include <Message_ProgressIndicator.hxx>
#include <GeomAlgoExt_Progress.hxx>
#include <GeomAlgoExt_Failure.hxx>
#include <Precision.hxx>
#include <cmath>

// A wrong basis offset or derivative row fails this closed-form Bernstein oracle.
TEST(GeomAlgoExtSplineFit, CubicBasisJets) {
 const std::vector<double> knots{0,0,0,0,1,1,1,1};
 for (double u : {0.0, .125, .5, .875, 1.0}) {
  const auto b = loft::splineBasis("basis", knots, 3, u);
  const double v = 1-u;
  const double value[4]{v*v*v,3*u*v*v,3*u*u*v,u*u*u};
  const double first[4]{-3*v*v,3*v*v-6*u*v,6*u*v-3*u*u,3*u*u};
  const double second[4]{6*v,18*u-12,6-18*u,6*u};
  const double third[4]{-6,18,-18,6};
  EXPECT_EQ(b.first,0);
  for (int i=0;i<4;++i) {
   EXPECT_NEAR(b.value[0][i],value[i],Precision::PConfusion());
   EXPECT_NEAR(b.value[1][i],first[i],Precision::PConfusion());
   EXPECT_NEAR(b.value[2][i],second[i],Precision::PConfusion());
   EXPECT_NEAR(b.value[3][i],third[i],Precision::PConfusion());
  }
 }
}
// Tridiagonal SPD matrices: the expected x is independent of the factor/solve.
TEST(GeomAlgoExtSplineFit, BandedCholeskyResidual) {
 for (bool paired : {false,true}) for (int n : {4,17,65}) {
  GeomAlgoExt::Band a(n,1,paired);
  std::vector<double> x(n), expected(n);
  for (int i=0;i<n;++i) expected[i]=std::sin(double(i));
  for (int i=0;i<n;++i) {
   a.add(i,i,2);
   if(i) a.add(i,i-1,-.25);
   x[i]=2*expected[i]-(i?.25*expected[i-1]:0)-(i+1<n?.25*expected[i+1]:0);
  }
  ASSERT_TRUE(a.factor());a.solve(x);
  for(int i=0;i<n;++i) EXPECT_NEAR(x[i],expected[i],Precision::PConfusion());
 }
}
// A cubic tensor polynomial is uniquely determined by its 4x4 values.
// Crossed U/V strides, incorrect fixed-pole elimination or non-finite results fail.
TEST(GeomAlgoExtSplineFit, TensorPolynomialInterpolation) {
 for (double scale : {.01,1.,100.}) {
  loft::SplineNet net;
  net.degreeU=net.degreeV=3;
  net.knotsU=net.knotsV={0,0,0,0,1,1,1,1};net.poles.resize(16,gp_XYZ(0,0,0));
  auto exact=[&](double u,double v) {return gp_Pnt(scale*(u+v*v),scale*(v-u*u),scale*(u*u*u-2*u*v+v*v*v));};
  std::vector<loft::FitCondition> conditions;
  for(int i=0;i<4;++i) for(int j=0;j<4;++j) {
   const double u=i/3.,v=j/3.;const auto p=exact(u,v);
   for(int c=0;c<3;++c) {
    loft::FitCondition row{}; row.u=u;row.v=v;row.count=1;
    row.terms[0]={0,0,gp_XYZ(c==0,c==1,c==2)};
    row.value=p.Coord(c+1);row.weight=1;row.within=Precision::Confusion()*scale;
    conditions.push_back(row);
   }
  }
  loft::fitNet("interpolation",net,std::vector<char>(16,0),conditions,{0,0},{});
  auto surface=loft::surfaceOf(net);
  for(int i=0;i<=16;++i) for(int j=0;j<=16;++j)
   EXPECT_LE(surface->Value(i/16.,j/16.).Distance(exact(i/16.,j/16.)),Precision::Confusion()*scale);
 }
}

namespace {
loft::SplineNet planeNet(double height) {
 loft::SplineNet net;net.degreeU=net.degreeV=3;
 net.knotsU=net.knotsV={0,0,0,0,1,1,1,1};
 for(int i=0;i<4;++i)for(int j=0;j<4;++j)net.poles.emplace_back(i/3.,j/3.,height);
 return net;
}
class StopIndicator : public Message_ProgressIndicator {
 bool UserBreak() override { return true; }
 void Show(const Message_ProgressScope&,const bool) override {}
};
}
// A wrong normal-system RHS or pole update fails the defining plane distances.
TEST(GeomAlgoExtSplineFit, RimCorrectionToParallelPlane) {
 auto net=planeNet(0);
 const auto places=[](double) {
  std::vector<loft::RimPlace> out;
  for(int i=0;i<=8;++i)for(int j=0;j<=8;++j)
   out.push_back({gp_Pnt(i/8.,j/8.,.125),gp_XY(i/8.,j/8.),gp_XY(0,0),false});
  return out;
 };
 const double error=loft::correctRim("rim",net,places,gp_XYZ(0,0,1),Precision::Confusion(),0,16);
 EXPECT_LE(error,Precision::Confusion());
 auto surface=loft::surfaceOf(net);
 for(int i=0;i<=16;++i)for(int j=0;j<=16;++j)
  EXPECT_NEAR(surface->Value(i/16.,j/16.).Z(),.125,Precision::Confusion());
}
// A false interval certificate would accept displacement greater than its bound.
TEST(GeomAlgoExtSplineFit, KnotRemovalBoundOnAffineSurface) {
 auto surface=loft::surfaceOf(planeNet(.125));
 surface->InsertUKnot(.5,1,Precision::PConfusion());
 surface->InsertVKnot(.5,1,Precision::PConfusion());
 auto reduced=loft::compactPatchSurface("compact",surface,Precision::Confusion());
 ASSERT_TRUE(reduced.has_value());
 for(int i=0;i<=32;++i)for(int j=0;j<=32;++j) {
  const auto point=reduced->surface->Value(i/32.,j/32.);
  EXPECT_LE(point.Distance(gp_Pnt(i/32.,j/32.,.125)),Precision::Confusion());
  EXPECT_LE(point.Distance(surface->Value(i/32.,j/32.)),reduced->displacement);
 }
}
// Optional optimization catches OCCT builder failures, but must never swallow a stop.
TEST(GeomAlgoExtSplineFit, CompactCancellationEscapesOptionalFailureCatch) {
 auto surface=loft::surfaceOf(planeNet(0));
 surface->InsertUKnot(.5,1,Precision::PConfusion());
 Handle(Message_ProgressIndicator) stop=new StopIndicator;
 const auto range=stop->Start();
 loft::ProgressScope progress(range);
 try { loft::compactPatchSurface("compact",surface,Precision::Confusion()); FAIL()<<"cancel did not escape"; }
 catch(const loft::Failure& failure) { EXPECT_EQ(failure.reason,"cancelled"); EXPECT_EQ(failure.message,"Operation cancelled"); }
 for(int i=0;i<=8;++i)for(int j=0;j<=8;++j)
  EXPECT_LE(surface->Value(i/8.,j/8.).Distance(gp_Pnt(i/8.,j/8.,0)),Precision::Confusion());
}
