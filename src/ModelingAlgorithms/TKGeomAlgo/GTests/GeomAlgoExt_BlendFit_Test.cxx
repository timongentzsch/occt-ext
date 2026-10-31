// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <GeomAlgoExt_BlendFit.hxx>
#include <GeomAlgoExt_Execution.hxx>
#include <Message_ProgressIndicator.hxx>
#include <Precision.hxx>
#include <gtest/gtest.h>
#include <cmath>
namespace {
GeomAlgoExt::BlendFitInput cylinder(double radius,double length) {
 GeomAlgoExt::BlendFitInput in{};
 in.degree=5;in.sectionDegree=2;in.sectionCount=3;
 in.sectionKnots={0,1};in.sectionMultiplicities={3,3};
 in.knots={{0,6},{1,6}};for(int k=1;k<8;++k)in.knots[k/8.]=1;
 in.scale=radius;in.tolerance=Precision::Confusion()*radius;
 in.knotMerge=1e-12;in.degenerate=1e-12;in.bendTolerance=1e-3;in.bendRefinements=6;
 in.row=[=](double u) {
  const double w=std::sqrt(.5);
  return std::vector<double>{radius,0,length*u,1,0,0,length*u*w,w,0,radius,length*u,1,u,0,u,0};
 };
 in.discoveredBreaks=[] {return std::map<double,int>{};};
 in.regime=[](double){return 0;};in.addBreak=[](std::map<double,int>& m,double u,int n){m[u]=n;};
 in.contactJet=[](int,double,double,GeomAlgoExt::BlendContactJet&){return false;};
 in.contactPoint=[=](int side,double u,double) {return gp_Pnt(side?0:radius,side?radius:0,length*u);};
 return in;
}
class Stop : public Message_ProgressIndicator {
 bool UserBreak() override {return true;}
 void Show(const Message_ProgressScope&,const bool) override {}
};
}
// Rational homogeneous interpolation must retain the defining circle and axial line.
// Swapped section/rail indices, lost weights and wrong contact coefficients fail.
TEST(GeomAlgoExtBlendFit, ExactQuarterCylinder) {
 for(double scale:{.01,1.,100.}) {
  const auto in=cylinder(scale,20*scale);
  const auto fit=GeomAlgoExt::fitBlendSections(in,{});
  for(int i=0;i<=24;++i)for(int j=0;j<=24;++j) {
   const double u=i/24.,v=j/24.;const auto p=fit.blend->Value(u,v);
   EXPECT_NEAR(std::hypot(p.X()-scale,p.Y()-scale),scale,Precision::Confusion()*scale);
   EXPECT_NEAR(p.Z(),20*scale*u,Precision::Confusion()*scale);
   if(j==0) {EXPECT_NEAR(p.X(),scale,Precision::Confusion()*scale);EXPECT_NEAR(p.Y(),0,Precision::Confusion()*scale);}
   if(j==24) {EXPECT_NEAR(p.X(),0,Precision::Confusion()*scale);EXPECT_NEAR(p.Y(),scale,Precision::Confusion()*scale);}
  }
  for(int side=0;side<2;++side)for(int k=0;k<=24;++k) {
   const auto uv=fit.contact[side]->Value(k/24.);
   EXPECT_NEAR(uv.X(),k/24.,Precision::PConfusion());EXPECT_NEAR(uv.Y(),0,Precision::PConfusion());
  }
 }
}
TEST(GeomAlgoExtBlendFit, CancellationPrecedesSampling) {
 auto in=cylinder(1,20);int calls=0;
 in.row=[&](double){++calls;return std::vector<double>{};};
 Handle(Message_ProgressIndicator) stop=new Stop;
 try {GeomAlgoExt::fitBlendSections(in,stop->Start());FAIL()<<"stop ignored";}
 catch(const GeomAlgoExt_Failure& e){EXPECT_STREQ(e.Reason(),"cancelled");EXPECT_STREQ(e.what(),"Operation cancelled");}
 EXPECT_EQ(calls,0);
}
