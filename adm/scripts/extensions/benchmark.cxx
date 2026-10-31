// Copyright (c) 2026 Timon Gentzsch
// This file is part of Open CASCADE Technology software library.
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult LICENSE_LGPL_21.txt for the full license.
#include <BRepAdaptor_Curve.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepGProp.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>
#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>
#include <ElSLib.hxx>
#include <Geom_Plane.hxx>
#include <GProp_GProps.hxx>
#include <Poly_Triangulation.hxx>
#include <Precision.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

// Workload dimensions and sample counts are arbitrary benchmark inputs, not solver constants.
// All measured intervals use process CPU time; setup is excluded except sphere construction; parallel meshing is off.
template <class F> void Measure(const char* theName, F theOperation, double theExpected = std::numeric_limits<double>::quiet_NaN())
{
  std::vector<double> aTimes;
  double aResult = 0;
  for (int i = 0; i < 11; ++i)
  {
    const std::clock_t aStart = std::clock();
    aResult = theOperation();
    const double aMs = 1000.0 * (std::clock() - aStart) / CLOCKS_PER_SEC;
    if (std::isfinite(theExpected) && std::abs(aResult-theExpected)>Precision::Confusion())
      throw std::runtime_error("result violates analytic volume oracle");
    if (i != 0) aTimes.push_back(aMs); // one warmup, ten observations
  }
  std::sort(aTimes.begin(), aTimes.end());
  std::cout << std::setprecision(17) << "{\"case\":\"" << theName
            << "\",\"samples\":10,\"cpuMedianMs\":" << (aTimes[4]+aTimes[5])/2
            << ",\"cpuMinMs\":" << aTimes.front() << ",\"cpuMaxMs\":" << aTimes.back()
            << ",\"value\":" << aResult << "}\n";
}
int main()
{
  TopoDS_Compound aTools;
  BRep_Builder aBuilder;
  aBuilder.MakeCompound(aTools);
  for (int j = 0; j < 4; ++j)
    for (int i = 0; i < 8; ++i)
      aBuilder.Add(aTools, BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(8+12*i,8+18*j,-1),gp_Dir(0,0,1)),1,4).Shape());
  const TopoDS_Shape aPlate = BRepAlgoAPI_Cut(BRepPrimAPI_MakeBox(100,80,2).Shape(),aTools).Shape();
  if (!BRepCheck_Analyzer(aPlate).IsValid()) throw std::runtime_error("invalid input plate");
  Measure("plate32-properties", [&] {
    GProp_GProps aSurface, aVolume;
    BRepGProp::SurfaceProperties(aPlate,aSurface);
    BRepGProp::VolumeProperties(aPlate,aVolume);
    const double aPi = std::acos(-1.0);
    if (std::abs(aVolume.Mass()-(16000-64*aPi)) > Precision::Confusion()
        || std::abs(aSurface.Mass()-(16720+64*aPi)) > Precision::Confusion())
      throw std::runtime_error("plate mass violates analytic cylindrical-hole formula");
    return aVolume.Mass();
  });
  TopoDS_Face aTop;
  for (TopExp_Explorer anIt(aPlate,TopAbs_FACE); anIt.More(); anIt.Next())
  {
    const auto aPlane = occ::down_cast<Geom_Plane>(BRep_Tool::Surface(TopoDS::Face(anIt.Current())));
    if (!aPlane.IsNull() && std::abs(aPlane->Location().Z()-2) < Precision::Confusion())
      aTop = TopoDS::Face(anIt.Current());
  }
  if (aTop.IsNull()) throw std::runtime_error("missing top face");
  const gp_Pln aPlane = occ::down_cast<Geom_Plane>(BRep_Tool::Surface(aTop))->Pln();
  std::vector<gp_Pnt2d> aQueries;
  int aExpected = 0;
  for (int j = 0; j < 100; ++j)
    for (int i = 0; i < 200; ++i)
    {
      const double x=(i+0.31)*0.5, y=(j+0.27)*0.8;
      double u,v;
      ElSLib::Parameters(aPlane,gp_Pnt(x,y,2),u,v);
      aQueries.emplace_back(u,v);
      bool isInside=true;
      for (int b=0;b<4;++b) for(int a=0;a<8;++a)
        if (std::hypot(x-(8+12*a),y-(8+18*b))<1) isInside=false;
      aExpected+=isInside;
    }
  Measure("plate32-classification", [&] {
    BRepTopAdaptor_FClass2d aClassifier(aTop,Precision::PConfusion());
    int aCount=0;
    for(const auto& p:aQueries) aCount+=aClassifier.Perform(p)==TopAbs_IN;
    if(aCount!=aExpected) throw std::runtime_error("classification violates analytic disk complement");
    return double(aCount);
  });
  Measure("plate32-fillet", [&] {
    BRepFilletAPI_MakeFillet aFillet(aPlate);
    for(TopExp_Explorer anIt(aTop,TopAbs_EDGE);anIt.More();anIt.Next())
    {
      const TopoDS_Edge e=TopoDS::Edge(anIt.Current());
      if(BRepAdaptor_Curve(e).GetType()==GeomAbs_Circle) aFillet.Add(0.2,e);
    }
    aFillet.Build();
    if(!aFillet.IsDone()) throw std::runtime_error("fillet failed");
    GProp_GProps aMass;
    BRepGProp::VolumeProperties(aFillet.Shape(),aMass);
    return aMass.Mass();
  }, 16000-64*std::acos(-1.0)
       -32*std::acos(-1.0)*(2*0.2*0.2+(5.0/3)*0.2*0.2*0.2
                           -(std::acos(-1.0)/2)*(1+0.2)*0.2*0.2));
  Measure("sphere-mesh", [&] {
    const TopoDS_Shape aSphere=BRepPrimAPI_MakeSphere(10).Shape();
    BRepMesh_IncrementalMesh aMesh(aSphere,0.01,false,0.1,false);
    int aTriangles=0;
    for(TopExp_Explorer anIt(aSphere,TopAbs_FACE);anIt.More();anIt.Next())
    {
      TopLoc_Location aLocation;
      const auto t=BRep_Tool::Triangulation(TopoDS::Face(anIt.Current()),aLocation);
      if(t.IsNull()) throw std::runtime_error("missing mesh");
      aTriangles+=t->NbTriangles();
      for(int i=1;i<=t->NbNodes();++i)
        if(std::abs(t->Node(i).Transformed(aLocation.Transformation()).Distance(gp_Pnt())-10)>Precision::Confusion())
          throw std::runtime_error("mesh vertex is off analytic sphere");
    }
    return double(aTriangles);
  });
}
