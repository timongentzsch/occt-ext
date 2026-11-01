// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef OCCT_EXTENSION_STOCK_ONLY
#include <GeomAlgoExt_PatchFit.hxx>
#endif
#include <GeomPlate_BuildPlateSurface.hxx>
#include <GeomPlate_CurveConstraint.hxx>
#include <GeomPlate_Surface.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_Plane.hxx>
#include <gp_Pln.hxx>
#include <gp_XY.hxx>
#include <Precision.hxx>
#include <algorithm>
#include <cstdio>
#include <ctime>
#include <vector>
#include <exception>

static gp_Pnt exact(double u,double v,double scale) { return gp_Pnt(scale*u,scale*v,scale*u*v); }
int main() {
 for(double scale:{.1,1.,10.}) for(bool stock:
#ifdef OCCT_EXTENSION_STOCK_ONLY
 {true}
#else
 {false}
#endif
 ) {
  std::vector<double> times;double largest=0, boundary=0;int refused=0;
  for(int sample=-1;sample<10;++sample) {
#ifndef OCCT_EXTENSION_STOCK_ONLY
   loft::SplineNet net; net.degreeU=net.degreeV=3;
   net.knotsU=net.knotsV={0,0,0,0,1,1,1,1};
   for(int i=0;i<4;++i)for(int j=0;j<4;++j)net.poles.emplace_back(scale*i/3.,scale*j/3.,0);
   std::vector<loft::FitCondition> rows;
#endif
   const double tolerance=Precision::Confusion()*scale;
   GeomPlate_BuildPlateSurface plate(new Geom_Plane(gp_Pln(gp::XOY())),3,17,3,Precision::PConfusion(),tolerance,.01,.1,false);
   for(int edge=0;edge<4;++edge) {
    const auto uv=[&](double t) {return gp_XY(edge==0?t:edge==1?1:edge==2?1-t:0,edge==0?0:edge==1?t:edge==2?1:1-t);};
    NCollection_Array1<gp_Pnt> ends(1,2);const auto a=uv(0),b=uv(1);
    ends(1)=exact(a.X(),a.Y(),scale);ends(2)=exact(b.X(),b.Y(),scale);
    Handle(Geom_Curve) line=new Geom_BezierCurve(ends);
    plate.Add(new GeomPlate_CurveConstraint(new GeomAdaptor_Curve(line),0,17,tolerance));
#ifndef OCCT_EXTENSION_STOCK_ONLY
    for(int k=0;k<=16;++k) {
     const auto p=uv(k/16.);const auto target=exact(p.X(),p.Y(),scale);
     for(int c=0;c<3;++c) {
      loft::FitCondition row{};row.u=p.X();row.v=p.Y();row.count=1;
      row.terms[0]={0,0,gp_XYZ(c==0,c==1,c==2)};row.value=target.Coord(c+1);
      row.weight=1;row.within=tolerance;rows.push_back(row);
     }
    }
 #endif
   }
   try {
    const auto begin=std::clock();
    if(stock)plate.Perform();
#ifndef OCCT_EXTENSION_STOCK_ONLY
    else loft::fitNet("benchmark",net,std::vector<char>(16,0),rows,{1e-6,1e-8},{});
#endif
    const double cpu=1000.*(std::clock()-begin)/CLOCKS_PER_SEC;
    if(sample<0)continue;
    times.push_back(cpu);
    Handle(Geom_Surface) surface;
    if(stock) {if(!plate.IsDone()||plate.Surface().IsNull()){++refused;continue;}surface=plate.Surface();}
#ifndef OCCT_EXTENSION_STOCK_ONLY
    else surface=loft::surfaceOf(net);
#endif
    for(int i=0;i<=32;++i)for(int j=0;j<=32;++j) {
     const double u=i/32.,v=j/32.;
     const auto got=surface->Value(stock?u*scale:u,stock?v*scale:v);
     const double error=got.Distance(exact(u,v,scale));largest=std::max(largest,error);
     if(i==0||j==0||i==32||j==32)boundary=std::max(boundary,error);
    }
   }catch(const Standard_Failure&){if(sample>=0)++refused;}
   catch(const std::exception&){if(sample>=0)++refused;}
  }
  std::sort(times.begin(),times.end());
  const double median=times.empty()?0:(times[(times.size()-1)/2]+times[times.size()/2])/2;
  std::printf("{\"case\":\"bilinear-boundary\",\"scale\":%.17g,\"method\":\"%s\",\"samples\":%zu,\"refused\":%d,\"cpuMedianMs\":%.17g,\"cpuMinMs\":%.17g,\"cpuMaxMs\":%.17g,\"maxGridError\":%.17g,\"maxBoundaryError\":%.17g}\n",scale,stock?"GeomPlate":"SplineFit",times.size(),refused,median,times.empty()?0:times.front(),times.empty()?0:times.back(),largest,boundary);
 }
}
