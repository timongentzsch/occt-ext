// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#include <GeomAlgoExt_Bernstein.hxx>
#include <GeomAlgoExt_Interval.hxx>
#include <GeomAlgoExt_UniformKnots.hxx>
#include <GeomAlgoExt_Hermite.hxx>
#include <GeomAlgoExt_CubicThrough.hxx>
#include <GeomAlgoExt_ParameterRectangle.hxx>
#include <GeomAlgoExt_SurfaceBounds.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <Geom_Plane.hxx>
#include <Geom_BezierSurface.hxx>
#include <Precision.hxx>
#include <gp_Pln.hxx>
#include <iostream>
#include <stdexcept>
static void require(bool condition) { if (!condition) throw std::runtime_error("mathematical oracle failed"); }
static bool near(double a, double b) { return std::abs(a-b) <= Precision::Confusion(); }
static double eval(const std::vector<double>& c, int n, double t) {
  double value=0;
  for(int i=0;i<=n;++i) value+=c[i]*GeomAlgoExt::choose(n,i)*std::pow(t,i)*std::pow(1-t,n-i);
  return value;
}
int main() {
  using namespace GeomAlgoExt;
  // f(t)=1+2t; g(t)=2-t. Compare product/elevation against the polynomial itself.
  const auto p=product(std::vector<double>{1,3},1,0,std::vector<double>{2,1},1,0);
  const auto e=elevate(std::vector<double>{1,3},1,0,5,0);
  for(int i=0;i<=32;++i) {
    const double t=double(i)/32;
    require(near(eval(p,2,t),(1+2*t)*(2-t)));
    require(near(eval(e,5,t),1+2*t));
  }
  // Long-double exact reference for sums/products of these binary64 inputs.
  for(double a: {-1.25,0.1,3.0}) for(double b: {-2.0,0.2,7.0}) {
    const auto sum=interval::add(interval::point(a),interval::point(b));
    const auto mul=interval::multiply(interval::point(a),interval::point(b));
    const long double s=(long double)a+b,m=(long double)a*b;
    require(sum.low<=s && s<=sum.high && mul.low<=m && m<=mul.high);
  }
  require(!interval::valid(interval::multiply(interval::invalid(),interval::point(0))));
  require(uniformKnots(2,2,-1,3)==std::vector<double>({-1,-1,-1,1,3,3,3}));
  // Quintic Hermite recovery of (t,t^3,t^5) from both endpoint jets.
  const auto knots=uniformKnots(5,1);
  const auto sites=hermiteSites(5,knots,2,2);
  std::vector<gp_Vec> jets;
  for(const auto& site:sites) {
    const double t=site.parameter;
    if(site.order==0) jets.emplace_back(t,t*t*t,std::pow(t,5));
    if(site.order==1) jets.emplace_back(1,3*t*t,5*std::pow(t,4));
    if(site.order==2) jets.emplace_back(0,6*t,20*std::pow(t,3));
  }
  require(interpolateHermite(5,knots,sites,jets));
  for(int axis=1;axis<=3;++axis) {
    std::vector<double> c;
    for(const auto& pole:jets) c.push_back(pole.Coord(axis));
    for(int i=0;i<=32;++i) require(near(eval(c,5,double(i)/32),std::pow(double(i)/32,2*axis-1)));
  }
  // Natural interpolation of an affine function must reproduce its whole graph.
  const std::vector<double> places{0,0.25,0.75,1};
  CubicThrough natural(places);
  std::vector<gp_Pnt> values;
  for(double t:places) values.emplace_back(t,2*t+1,-t);
  const auto poles=natural.through(values);
  NCollection_Array1<gp_Pnt> points(1,int(poles.size()));
  for(int i=1;i<=points.Length();++i) points(i)=poles[i-1];
  Handle(Geom_BSplineCurve) curve=new Geom_BSplineCurve(points,natural.knots,natural.mults,3);
  for(int i=0;i<=32;++i) {
    const double t=double(i)/32;const auto point=curve->Value(t);
    require(point.Distance(gp_Pnt(t,2*t+1,-t))<=Precision::Confusion());
  }
  const TopoDS_Face face=BRepBuilderAPI_MakeFace(gp_Pln(gp::XOY()),-2,3,-4,5).Face();
  require(parameterRectangle(face));
  NCollection_Array2<gp_Pnt> net(1,2,1,2);
  net(1,1)=gp_Pnt(-2,-4,0);net(2,1)=gp_Pnt(3,-4,0);
  net(1,2)=gp_Pnt(-2,5,0);net(2,2)=gp_Pnt(3,5,0);
  Handle(Geom_BezierSurface) patch=new Geom_BezierSurface(net);
  const auto box=bezierBounds(patch);
  for(int i=0;i<=8;++i) for(int j=0;j<=8;++j) require(!box.IsOut(patch->Value(double(i)/8,double(j)/8)));
  std::cout<<"7 utility headers; polynomial, interval, Hermite, affine spline and hull oracles passed\n";
}
