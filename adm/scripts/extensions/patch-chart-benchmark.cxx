// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef OCCT_EXTENSION_STOCK_ONLY
#include <BRepAlgoExt_PatchFlow.hxx>
#include <BRepAlgoExt_Interior.hxx>
#endif
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <GeomLProp_SLProps.hxx>
#include <Poly_Triangulation.hxx>
#include <gp_Pln.hxx>
#include <gp_Cylinder.hxx>
#include <TopoDS_Face.hxx>
#include <Precision.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <vector>
struct Chart {std::vector<gp_Pnt> point;std::vector<gp_XY> flat;std::vector<std::array<int,3>> triangles;};
int main() {
 for(bool cylinder:{false,true})for(double scale:{.1,1.,10.})for(bool survey:{false,true}) {
#ifdef OCCT_EXTENSION_STOCK_ONLY
  const bool stock=true;
#else
  const bool stock=false;
#endif
  std::vector<double> times;double angle=0,edgeError=0,curvatureError=0;size_t nodes=0,triangles=0;int refused=0;
  for(int sample=-1;sample<10;++sample) {
   const auto face=cylinder?BRepBuilderAPI_MakeFace(gp_Cylinder(gp_Ax3(),scale),0,M_PI/2,0,2*scale).Face():BRepBuilderAPI_MakeFace(gp_Pln(gp::XOY()),0,2*scale,0,scale).Face();
   Chart out;double bend=0;
   try {
    const auto start=std::clock();
#ifdef OCCT_EXTENSION_STOCK_ONLY
    if(survey) {
     const auto surface=BRep_Tool::Surface(face);double u0,u1,v0,v1;BRepTools::UVBounds(face,u0,u1,v0,v1);
     for(int i=0;i<40;++i)for(int j=0;j<40;++j) {
      GeomLProp_SLProps p(surface,u0+(u1-u0)*(i+.5)/40,v0+(v1-v0)*(j+.5)/40,2,1e-12);
      if(p.IsCurvatureDefined())bend=std::max({bend,std::abs(p.MaxCurvature()),std::abs(p.MinCurvature())});
     }
    }else { BRepMesh_IncrementalMesh mesh(face,.01*scale,false,.5,false); }
#else
    if(survey)bend=BRepAlgoExt::maximumInteriorCurvature(face,Precision::PConfusion());
    else {auto c=loft::conformalChart(face,.2*scale,.01*scale);out.point=std::move(c.point);out.flat=std::move(c.flat);out.triangles=std::move(c.triangles);}
#endif
    const double cpu=1000.*(std::clock()-start)/CLOCKS_PER_SEC;
    if(sample<0)continue;times.push_back(cpu);
    if(survey) {curvatureError=std::max(curvatureError,std::abs(bend-(cylinder?1/scale:0)));continue;}
#ifdef OCCT_EXTENSION_STOCK_ONLY
    TopLoc_Location location;const auto mesh=BRep_Tool::Triangulation(face,location);
    if(mesh.IsNull()||!mesh->HasUVNodes()){++refused;continue;}
    for(int k=1;k<=mesh->NbNodes();++k) {
     out.point.push_back(mesh->Node(k).Transformed(location.Transformation()));const auto uv=mesh->UVNode(k);
     out.flat.emplace_back(cylinder?scale*uv.X():uv.X(),uv.Y());
    }
    for(int k=1;k<=mesh->NbTriangles();++k) {int a,b,c;mesh->Triangle(k).Get(a,b,c);out.triangles.push_back({a-1,b-1,c-1});}
#endif
    nodes=out.point.size();triangles=out.triangles.size();
    for(const auto& t:out.triangles)for(int k=0;k<3;++k) {
     const int a=t[k],b=t[(k+1)%3],c=t[(k+2)%3];const gp_Vec p(out.point[a],out.point[b]),q(out.point[a],out.point[c]);
     const gp_XY x=out.flat[b]-out.flat[a],y=out.flat[c]-out.flat[a];
     angle=std::max(angle,std::abs(p.Dot(q)/(p.Magnitude()*q.Magnitude())-x.Dot(y)/(x.Modulus()*y.Modulus())));
     edgeError=std::max(edgeError,std::abs(x.Modulus()/p.Magnitude()-1));
    }
   }catch(const Standard_Failure&){if(sample>=0)++refused;}
  }
  std::sort(times.begin(),times.end());
  std::printf("{\"case\":\"%s-%s\",\"scale\":%.17g,\"method\":\"%s\",\"samples\":%zu,\"refused\":%d,\"cpuMedianMs\":%.17g,\"cpuMinMs\":%.17g,\"cpuMaxMs\":%.17g,\"nodes\":%zu,\"triangles\":%zu,\"maxCosineDifference\":%.17g,\"maxRelativeEdgeError\":%.17g,\"maxCurvatureError\":%.17g}\n",cylinder?"cylinder":"plane",survey?"curvature-survey":"chart",scale,stock?(survey?"GeomLProp":"BRepMesh-and-carrier-UV"):(survey?"InteriorSurvey":"BoundaryFirstChart"),times.size(),refused,times.empty()?0:(times[(times.size()-1)/2]+times[times.size()/2])/2,times.empty()?0:times.front(),times.empty()?0:times.back(),nodes,triangles,angle,edgeError,curvatureError);
 }
}
