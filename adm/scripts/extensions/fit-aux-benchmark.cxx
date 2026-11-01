// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef OCCT_EXTENSION_STOCK_ONLY
#include <GeomAlgoExt_PatchRim.hxx>
#include <GeomAlgoExt_PatchCompact.hxx>
#endif
#include <Geom_BSplineSurface.hxx>
#include <Geom_Plane.hxx>
#include <GeomPlate_BuildPlateSurface.hxx>
#include <GeomPlate_PointConstraint.hxx>
#include <GeomPlate_Surface.hxx>
#include <gp_Pln.hxx>
#include <Precision.hxx>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <vector>
#include <exception>
static Handle(Geom_BSplineSurface) plane(double scale) {
 NCollection_Array2<gp_Pnt> poles(1,4,1,4);NCollection_Array1<double> knots(1,2);NCollection_Array1<int> mults(1,2);
 knots(1)=0;knots(2)=1;mults.Init(4);
 for(int i=1;i<=4;++i)for(int j=1;j<=4;++j)poles(i,j)=gp_Pnt(scale*(i-1)/3.,scale*(j-1)/3.,0);
 return new Geom_BSplineSurface(poles,knots,knots,mults,mults,3,3);
}
int main() {
 for(bool compact:{false,true})for(double scale:{.1,1.,10.}) {
#ifdef OCCT_EXTENSION_STOCK_ONLY
  const bool stock=true;
#else
  const bool stock=false;
#endif
  std::vector<double> times;double error=0;int refused=0;
  for(int sample=-1;sample<10;++sample) {
   auto input=plane(scale);const double tolerance=Precision::Confusion()*scale;
   if(compact)for(double u:{.25,.5,.75}){input->InsertUKnot(u,1,Precision::PConfusion());input->InsertVKnot(u,1,Precision::PConfusion());}
   GeomPlate_BuildPlateSurface plate(new Geom_Plane(gp_Pln(gp::XOY())),3,17,3,Precision::PConfusion(),tolerance,.01,.1,false);
   for(int i=0;i<=8;++i)for(int j=0;j<=8;++j)plate.Add(new GeomPlate_PointConstraint(gp_Pnt(scale*i/8.,scale*j/8.,scale*.125),0,tolerance));
#ifndef OCCT_EXTENSION_STOCK_ONLY
   auto net=BRepExt::netOf(input);
   const auto places=[=](double) {std::vector<BRepExt::RimPlace> out;
    for(int i=0;i<=8;++i)for(int j=0;j<=8;++j)out.push_back({gp_Pnt(scale*i/8.,scale*j/8.,scale*.125),gp_XY(i/8.,j/8.),gp_XY(0,0),false});return out;};
#endif
   try {
    Handle(Geom_Surface) surface;const auto begin=std::clock();
#ifdef OCCT_EXTENSION_STOCK_ONLY
    if(compact) {
     for(int i=input->NbUKnots()-1;i>1;--i)if(!input->RemoveUKnot(i,0,tolerance))++refused;
     for(int i=input->NbVKnots()-1;i>1;--i)if(!input->RemoveVKnot(i,0,tolerance))++refused;
    }else plate.Perform();
#else
    if(compact){const auto reduced=BRepExt::compactPatchSurface("compact",input,tolerance);if(reduced)input=reduced->surface;else if(sample>=0)++refused;}
    else BRepExt::correctRim("benchmark",net,places,gp_XYZ(0,0,1),tolerance,0,16);
#endif
    const double cpu=1000.*(std::clock()-begin)/CLOCKS_PER_SEC;
    if(sample<0)continue;times.push_back(cpu);
    if(compact)surface=input;
#ifdef OCCT_EXTENSION_STOCK_ONLY
    else {if(!plate.IsDone()){++refused;continue;}surface=plate.Surface();}
#else
    else surface=BRepExt::surfaceOf(net);
#endif
    for(int i=0;i<=32;++i)for(int j=0;j<=32;++j) {
     const double u=i/32.,v=j/32.;const auto p=surface->Value(stock&&!compact?scale*u:u,stock&&!compact?scale*v:v);
     error=std::max(error,p.Distance(gp_Pnt(scale*u,scale*v,compact?0:scale*.125)));
    }
   }catch(const Standard_Failure&){if(sample>=0)++refused;}
   catch(const std::exception&){if(sample>=0)++refused;}
  }
  std::sort(times.begin(),times.end());
  std::printf("{\"case\":\"%s\",\"scale\":%.17g,\"method\":\"%s\",\"samples\":%zu,\"refused\":%d,\"cpuMedianMs\":%.17g,\"cpuMinMs\":%.17g,\"cpuMaxMs\":%.17g,\"maxGridDistance\":%.17g}\n",compact?"affine-knot-compaction":"parallel-plane-rim",scale,stock?(compact?"RemoveUKnot-RemoveVKnot":"GeomPlate-points"):(compact?"CertifiedCompaction":"RimFit"),times.size(),refused,times.empty()?0:(times[(times.size()-1)/2]+times[times.size()/2])/2,times.empty()?0:times.front(),times.empty()?0:times.back(),error);
 }
}
