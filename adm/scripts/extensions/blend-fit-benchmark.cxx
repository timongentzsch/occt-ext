// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef OCCT_EXTENSION_STOCK_ONLY
#include <GeomAlgoExt_BlendFit.hxx>
#endif
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepTools.hxx>
#include <TopoDS.hxx>
#include <TopExp_Explorer.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <vector>
#ifndef OCCT_EXTENSION_STOCK_ONLY
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

#endif
int main() {
 for(double scale:{.1,1.,10.}) {
#ifdef OCCT_EXTENSION_STOCK_ONLY
  const bool stock=true;
#else
  const bool stock=false;
#endif
  std::vector<double> times;double largest=0;int refused=0;
  for(int sample=-1;sample<10;++sample) {
   const auto box=BRepPrimAPI_MakeBox(10*scale,10*scale,20*scale).Shape();
   TopoDS_Edge edge;
   for(TopExp_Explorer it(box,TopAbs_EDGE);it.More();it.Next()) {
    const auto candidate=TopoDS::Edge(it.Current());BRepAdaptor_Curve c(candidate);
    const auto a=c.Value(c.FirstParameter()),b=c.Value(c.LastParameter());
    if(std::abs(a.X())<Precision::Confusion()*scale&&std::abs(a.Y())<Precision::Confusion()*scale&&
       std::abs(b.X())<Precision::Confusion()*scale&&std::abs(b.Y())<Precision::Confusion()*scale&&
       std::abs(a.Z()-b.Z())>19*scale){edge=candidate;break;}
   }
   BRepFilletAPI_MakeFillet fillet(box);fillet.Add(scale,edge);
#ifndef OCCT_EXTENSION_STOCK_ONLY
   const auto input=cylinder(scale,20*scale);
   GeomAlgoExt::BlendFitResult result;
#endif
   try {
    const auto start=std::clock();
#ifdef OCCT_EXTENSION_STOCK_ONLY
    fillet.Build();
#else
    result=GeomAlgoExt::fitBlendSections(input,{});
#endif
    const double cpu=1000.*(std::clock()-start)/CLOCKS_PER_SEC;
    if(sample<0)continue;
    times.push_back(cpu);
#ifdef OCCT_EXTENSION_STOCK_ONLY
    if(!fillet.IsDone()){++refused;continue;}
    int cylinders=0;
    for(TopExp_Explorer it(fillet.Shape(),TopAbs_FACE);it.More();it.Next()) {
     const auto face=TopoDS::Face(it.Current());BRepAdaptor_Surface surface(face);
     if(surface.GetType()!=GeomAbs_Cylinder)continue;++cylinders;
     double u0,u1,v0,v1;BRepTools::UVBounds(face,u0,u1,v0,v1);
     for(int i=0;i<=32;++i)for(int j=0;j<=32;++j) {
      const auto p=surface.Value(u0+(u1-u0)*i/32.,v0+(v1-v0)*j/32.);
      largest=std::max(largest,std::abs(std::hypot(p.X()-scale,p.Y()-scale)-scale));
     }
    }
    if(cylinders!=1)++refused;
#else
    for(int i=0;i<=32;++i)for(int j=0;j<=32;++j) {
     const auto p=result.blend->Value(i/32.,j/32.);
     largest=std::max({largest,std::abs(std::hypot(p.X()-scale,p.Y()-scale)-scale),std::abs(p.Z()-20*scale*i/32.)});
    }
#endif
   }catch(const Standard_Failure&){if(sample>=0)++refused;}
  }
  std::sort(times.begin(),times.end());
  std::printf("{\"case\":\"quarter-cylinder-box-edge\",\"scale\":%.17g,\"method\":\"%s\",\"samples\":%zu,\"refused\":%d,\"cpuMedianMs\":%.17g,\"cpuMinMs\":%.17g,\"cpuMaxMs\":%.17g,\"maxDefiningDistance\":%.17g}\n",scale,stock?"BRepFilletAPI_MakeFillet":"BlendSectionFit",times.size(),refused,times.empty()?0:(times[4]+times[5])/2,times.empty()?0:times.front(),times.empty()?0:times.back(),largest);
 }
}
