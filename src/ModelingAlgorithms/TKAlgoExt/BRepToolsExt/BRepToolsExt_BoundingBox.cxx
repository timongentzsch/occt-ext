// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_SurfaceBounds.hxx>
#include <GeomAlgoExt_Bernstein.hxx>
#include <BRepAlgoExt_TrimDomain.hxx>
#include <GeomAlgoExt_ParameterRectangle.hxx>
#include <BRepToolsExt_Session.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <TopExp_Explorer.hxx>
#include <BRepBndLib.hxx>
#include <Standard_Failure.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Curve2d.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <BndLib_Add3dCurve.hxx>
#include <BndLib_AddSurface.hxx>
#include <GeomConvert_BSplineSurfaceToBezierSurface.hxx>
#include <Geom_BezierSurface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <GeomConvert.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <gp_Torus.hxx>
#include <limits>
#include <queue>

namespace loft {
namespace {
std::vector<double> extremaAngles(double first,double last,double phase) {
    std::vector<double> out{first,last};
    const double pi=std::acos(-1.0);
    // Both signs of the trigonometric extremum, translated into the trim range.
    for(int sign=0;sign<2;++sign) {
        const double a=phase+sign*pi;
        const double t=a+2*pi*std::ceil((first-a)/(2*pi));
        if(t>=first && t<=last) out.push_back(t);
    }
    return out;
}
void torusBox(const BRepAdaptor_Surface& surface,double u0,double u1,double v0,double v1,Bnd_Box& box) {
    const auto torus=surface.Torus();
    const auto x=torus.Position().XDirection(),y=torus.Position().YDirection(),z=torus.Axis().Direction();
    for(int axis=1;axis<=3;++axis) {
        for(double u:extremaAngles(u0,u1,std::atan2(y.Coord(axis),x.Coord(axis)))) {
            const double radial=x.Coord(axis)*std::cos(u)+y.Coord(axis)*std::sin(u);
            for(double v:extremaAngles(v0,v1,std::atan2(z.Coord(axis),radial))) box.Add(surface.Value(u,v));
        }
    }
    box.Enlarge(1e-10);
}
using Limits=std::array<double,6>;
Limits limits(const Bnd_Box& b) {Limits out;b.Get(out[0],out[1],out[2],out[3],out[4],out[5]);return out;}
// Positive-weight rational Bezier patches lie inside their pole hull, and their corner poles
// are on them. Each of the six extremes is its own best-first search: the patch whose hull
// reaches furthest is halved, until no hull reaches beyond the furthest point found by more
// than the requested error. Returned bounds are hulls, never just the points found: a search
// that runs out of halvings keeps the hull of what is left, wider than the error but enclosing.
bool splineBox(const Handle(Geom_BSplineSurface)& surface,double u0,double u1,double v0,double v1,Bnd_Box& result) {
    const auto patches=GeomAlgoExt::surfacePatches(surface,u0,u1,v0,v1);
    Bnd_Box sampled;
    for(const auto& p:patches) for(double u:{0.0,.5,1.0}) for(double v:{0.0,.5,1.0}) sampled.Add(p->Value(u,v));
    if(sampled.IsVoid()) return false;
    const auto initial=limits(sampled);
    struct Cell {Handle(Geom_BezierSurface) patch;double reach;int depth;};
    const auto nearer=[](const Cell& a,const Cell& b){return a.reach<b.reach;};
    Limits bounds;
    for(int axis=1;axis<=3;++axis) for(double sign:{-1.0,1.0}) {
        // Signed so that every search is for a maximum.
        const int slot=sign>0?axis+2:axis-1;
        const double tolerance=std::max(1e-9,(initial[axis+2]-initial[axis-1])*1e-8);
        double found=sign*initial[slot];
        const auto reach=[&](const Handle(Geom_BezierSurface)& p) {
            double out=-std::numeric_limits<double>::infinity();
            for(int i=1;i<=p->NbUPoles();++i) for(int j=1;j<=p->NbVPoles();++j) out=std::max(out,sign*p->Pole(i,j).Coord(axis));
            return out;
        };
        std::priority_queue<Cell,std::vector<Cell>,decltype(nearer)> open(nearer);
        for(const auto& p:patches) open.push({p,reach(p),0});
        for(size_t visits=0;!open.empty() && open.top().reach>found+tolerance && visits<20000;++visits) {
            const auto cell=open.top();open.pop();
            if(cell.depth>=48) {found=std::max(found,cell.reach);continue;}
            // Halve across the way this coordinate of the poles bends most: along a way it is
            // linear in, halving brings the hull no nearer.
            const auto& p=cell.patch;double uBend=0,vBend=0;
            for(int i=1;i<=p->NbUPoles();++i) for(int j=1;j<=p->NbVPoles();++j) {
                const double c=2*p->Pole(i,j).Coord(axis);
                if(i>1 && i<p->NbUPoles()) uBend=std::max(uBend,std::abs(p->Pole(i-1,j).Coord(axis)+p->Pole(i+1,j).Coord(axis)-c));
                if(j>1 && j<p->NbVPoles()) vBend=std::max(vBend,std::abs(p->Pole(i,j-1).Coord(axis)+p->Pole(i,j+1).Coord(axis)-c));
            }
            for(int half=0;half<2;++half) {
                auto child=Handle(Geom_BezierSurface)::DownCast(p->Copy());
                if(uBend>=vBend) child->Segment(half*.5,(half+1)*.5,0,1);
                else child->Segment(0,1,half*.5,(half+1)*.5);
                for(int i:{1,child->NbUPoles()}) for(int j:{1,child->NbVPoles()}) found=std::max(found,sign*child->Pole(i,j).Coord(axis));
                open.push({child,reach(child),cell.depth+1});
            }
        }
        bounds[slot]=sign*(open.empty()?found:std::max(found,open.top().reach));
    }
    result.Update(bounds[0],bounds[1],bounds[2],bounds[3],bounds[4],bounds[5]);
    result.Enlarge(1e-10);
    return true;
}
// Whether a patch's coordinate is monotone along u or along v across it: its extremes over
// any region of the patch are then on that region's boundary. With positive weights the
// derivative has the sign of X'W-XW', itself a tensor Bernstein polynomial.
unsigned monotoneDirections(const Handle(Geom_BezierSurface)& patch,int axis) {
    const int n=patch->NbUPoles()-1,m=patch->NbVPoles()-1;
    const double origin=patch->Pole(1,1).Coord(axis);
    std::vector<double> x,w;
    for(int i=1;i<=n+1;++i) for(int j=1;j<=m+1;++j) {w.push_back(patch->Weight(i,j));x.push_back((patch->Pole(i,j).Coord(axis)-origin)*w.back());}
    const bool rational=patch->IsURational() || patch->IsVRational();
    auto difference=[&](const std::vector<double>& c,bool along) {
        std::vector<double> out;
        for(int i=0;i<=n-(along?1:0);++i) for(int j=0;j<=m-(along?0:1);++j)
            out.push_back(along ? c[(i+1)*(m+1)+j]-c[i*(m+1)+j] : c[i*(m+1)+j+1]-c[i*(m+1)+j]);
        return out;
    };
    unsigned directions=0;
    for(bool along:{true,false}) {
        if((along?n:m)==0) {directions|=along?3u:12u;continue;}
        auto sign=difference(x,along);
        if(rational) {
            const int dn=along?n-1:n,dm=along?m:m-1;
            sign=GeomAlgoExt::product(sign,dn,dm,w,n,m);const auto other=GeomAlgoExt::product(x,n,m,difference(w,along),dn,dm);
            for(size_t i=0;i<sign.size();++i) sign[i]-=other[i];
        }
        const auto [lo,hi]=std::minmax_element(sign.begin(),sign.end());
        if(*lo>=0) directions|=along?1u:4u;
        if(*hi<=0) directions|=along?2u:8u;
    }
    return directions;
}
// A spline face trimmed to less than its parameter rectangle: its edges and the surface along
// their pcurves bound it wherever its extremes are on its boundary, and pole hulls of Bezier
// cells inside the trim find the extremes inside it. As for an untrimmed face, each of the six
// extremes is its own best-first search: the cell whose hull reaches furthest is halved across
// the way its poles bend most in that coordinate, its corners inside the trim raising what is
// found, and a cell in which the coordinate is monotone is settled by the exact extremes of its
// side curves inside the trim, its others being on the rim. Cells outside the trim are dropped,
// so nothing of the support beyond the trim is enclosed.
bool trimmedSplineBox(const TopoDS_Face& face,const Handle(Geom_BSplineSurface)& surface,double u0,double u1,double v0,double v1,Bnd_Box& result) try {
    Bnd_Box edges;
    for(TopExp_Explorer it(face,TopAbs_EDGE);it.More();it.Next()) {
        const auto edge=TopoDS::Edge(it.Current());
        if(BRep_Tool::Degenerated(edge)) continue;
        if(!BRep_Tool::IsGeometric(edge)) return false;
        // The edge, and the surface along its pcurve, which is within the edge's tolerance of it.
        BRepAdaptor_Curve curve(edge),along(edge,face);
        BndLib_Add3dCurve::AddOptimal(curve,0,edges);BndLib_Add3dCurve::AddOptimal(along,0,edges);
    }
    if(edges.IsVoid()) return false;
    const TrimDomain domain(face);const auto rows=domain.bands({});
    const auto rim=limits(edges);
    struct Cell {Handle(Geom_BezierSurface) patch;double a,b,c,d;int depth;bool crossed;double reach;};
    std::vector<Cell> cells;
    auto opened=Handle(Geom_BSplineSurface)::DownCast(surface->Copy());
    if(opened->IsUPeriodic()) opened->SetUNotPeriodic();
    if(opened->IsVPeriodic()) opened->SetVNotPeriodic();
    GeomConvert_BSplineSurfaceToBezierSurface convert(opened,u0,u1,v0,v1,1e-10);
    NCollection_Array1<double> uk(1,convert.NbUPatches()+1),vk(1,convert.NbVPatches()+1);convert.UKnots(uk);convert.VKnots(vk);
    // A cell of the trim, or none: one the boundary does not cross is all inside or all out.
    const auto keep=[&](const Handle(Geom_BezierSurface)& patch,double a,double b,double c,double d,int depth,std::vector<Cell>& into) {
        const bool crossed=domain.boundaryIntersects(a,b,c,d);
        if(!crossed && !domain.inside(rows,(a+b)/2,(c+d)/2)) return;
        into.push_back({patch,a,b,c,d,depth,crossed,0});
    };
    for(int i=1;i<=convert.NbUPatches();++i) for(int j=1;j<=convert.NbVPatches();++j)
        keep(convert.Patch(i,j),uk(i),uk(i+1),vk(j),vk(j+1),0,cells);
    Limits bounds;
    for(int axis=1;axis<=3;++axis) for(double sign:{-1.0,1.0}) {
        // Signed so that every search is for a maximum; the rim's box bounds the boundary.
        const int slot=sign>0?axis+2:axis-1;
        const double tolerance=std::max(1e-9,(rim[axis+2]-rim[axis-1])*1e-8);
        double found=sign*rim[slot];
        const auto reach=[&](const Handle(Geom_BezierSurface)& p) {
            double out=-std::numeric_limits<double>::infinity();
            for(int i=1;i<=p->NbUPoles();++i) for(int j=1;j<=p->NbVPoles();++j) out=std::max(out,sign*p->Pole(i,j).Coord(axis));
            return out;
        };
        // Its corner poles are on the surface: they count where they are in the trim.
        const auto corners=[&](const Cell& cell) {
            const auto& p=cell.patch;
            for(int i:{0,1}) for(int j:{0,1}) {
                if(cell.crossed && !domain.inside(rows,i?cell.b:cell.a,j?cell.d:cell.c)) continue;
                found=std::max(found,sign*p->Pole(i?p->NbUPoles():1,j?p->NbVPoles():1).Coord(axis));
            }
        };
        const auto nearer=[](const Cell& a,const Cell& b){return a.reach<b.reach;};
        std::priority_queue<Cell,std::vector<Cell>,decltype(nearer)> open(nearer);
        for(auto cell:cells) {corners(cell);cell.reach=reach(cell.patch);open.push(cell);}
        for(size_t visits=0;!open.empty() && open.top().reach>found+tolerance && visits<20000;++visits) {
            const auto cell=open.top();open.pop();
            if(cell.depth>=48) {found=std::max(found,cell.reach);continue;}
            if(monotoneDirections(cell.patch,axis)) {
                // Its extremes are on its sides: those inside the trim, read exactly from the
                // side curves; those on the rim, in the rim's box.
                bool read=true;
                for(bool alongU:{true,false}) for(double end:{0.0,1.0}) {
                    auto curve=Handle(Geom_BezierCurve)::DownCast(alongU ? cell.patch->UIso(end) : cell.patch->VIso(end));
                    if(curve.IsNull()) {read=false;continue;}
                    for(double t:coordinateExtrema(curve,axis)) {
                        const double u=cell.a+(alongU?end:t)*(cell.b-cell.a);
                        const double v=cell.c+(alongU?t:end)*(cell.d-cell.c);
                        if(domain.inside(rows,u,v)) found=std::max(found,sign*curve->Value(t).Coord(axis));
                    }
                }
                if(read) continue;
            }
            // Halve across the way this coordinate of the poles bends most.
            const auto& p=cell.patch;double uBend=0,vBend=0;
            for(int i=1;i<=p->NbUPoles();++i) for(int j=1;j<=p->NbVPoles();++j) {
                const double c=2*p->Pole(i,j).Coord(axis);
                if(i>1 && i<p->NbUPoles()) uBend=std::max(uBend,std::abs(p->Pole(i-1,j).Coord(axis)+p->Pole(i+1,j).Coord(axis)-c));
                if(j>1 && j<p->NbVPoles()) vBend=std::max(vBend,std::abs(p->Pole(i,j-1).Coord(axis)+p->Pole(i,j+1).Coord(axis)-c));
            }
            std::vector<Cell> children;
            for(int half=0;half<2;++half) {
                auto child=Handle(Geom_BezierSurface)::DownCast(p->Copy());
                if(uBend>=vBend) {
                    child->Segment(half*.5,(half+1)*.5,0,1);const double h=(cell.b-cell.a)/2;
                    keep(child,cell.a+half*h,cell.a+(half+1)*h,cell.c,cell.d,cell.depth+1,children);
                } else {
                    child->Segment(0,1,half*.5,(half+1)*.5);const double h=(cell.d-cell.c)/2;
                    keep(child,cell.a,cell.b,cell.c+half*h,cell.c+(half+1)*h,cell.depth+1,children);
                }
            }
            for(auto& child:children) {corners(child);child.reach=reach(child.patch);open.push(child);}
        }
        bounds[slot]=sign*(open.empty()?found:std::max(found,open.top().reach));
    }
    result.Update(bounds[0],bounds[1],bounds[2],bounds[3],bounds[4],bounds[5]);
    result.Enlarge(1e-10);
    return true;
} catch(const Failure&) {
    return false;
} catch(const Standard_Failure&) {
    return false;
}
}
std::vector<double> boundingBox(OcctKernel& kernel,uint32_t id) {
    const auto shape=get(kernel,id);Bnd_Box box;
    NCollection_IndexedMap<TopoDS_Shape,TopTools_ShapeMapHasher> boundedEdges;
    const auto addEdges=[&](const TopoDS_Shape& part) {
        for(TopExp_Explorer it(part,TopAbs_EDGE);it.More();it.Next()) {
            const auto edge=TopoDS::Edge(it.Current());
            if(boundedEdges.Contains(edge) || BRep_Tool::Degenerated(edge) || !BRep_Tool::IsGeometric(edge)) continue;
            boundedEdges.Add(edge);BRepAdaptor_Curve c(edge);BndLib_Add3dCurve::AddOptimal(c,0,box);
        }
    };
    for(TopExp_Explorer it(shape,TopAbs_FACE);it.More();it.Next()) {
        const auto face=TopoDS::Face(it.Current());BRepAdaptor_Surface surface(face);
        const auto type=surface.GetType();
        if(type==GeomAbs_Plane || type==GeomAbs_SurfaceOfExtrusion) {addEdges(face);continue;}
        double u0,u1,v0,v1;BRepTools::UVBounds(face,u0,u1,v0,v1);
        if(GeomAlgoExt::parameterRectangle(face,u0,u1,v0,v1)) {
            Bnd_Box local;
            if(type==GeomAbs_BSplineSurface && splineBox(surface.BSpline(),u0,u1,v0,v1,local)) {box.Add(local);continue;}
            if(type==GeomAbs_BezierSurface && splineBox(GeomConvert::SurfaceToBSplineSurface(surface.Bezier()),u0,u1,v0,v1,local)) {box.Add(local);continue;}
            if(type==GeomAbs_SurfaceOfRevolution) {
                // Convert in world coordinates. Rotating a box computed in the
                // revolution's local frame encloses far more than the surface.
                Handle(Geom_RectangularTrimmedSurface) trim=new Geom_RectangularTrimmedSurface(
                    BRep_Tool::Surface(face),u0,u1,v0,v1);
                const auto spline=GeomConvert::SurfaceToBSplineSurface(trim);
                if(splineBox(spline,u0,u1,v0,v1,local)) {box.Add(local);continue;}
            }
            if(type==GeomAbs_Torus) {torusBox(surface,u0,u1,v0,v1,box);continue;}
            if(type<=GeomAbs_Torus) {BndLib_AddSurface::AddOptimal(surface,u0,u1,v0,v1,0,box);continue;}
        } else if(type==GeomAbs_BSplineSurface || type==GeomAbs_BezierSurface) {
            Bnd_Box local;
            const auto spline=type==GeomAbs_BSplineSurface ? surface.BSpline() : GeomConvert::SurfaceToBSplineSurface(surface.Bezier());
            if(trimmedSplineBox(face,spline,u0,u1,v0,v1,local)) {box.Add(local);continue;}
        }
        BRepBndLib::AddOptimal(face,box,false,false);
    }
    for(TopExp_Explorer it(shape,TopAbs_EDGE,TopAbs_FACE);it.More();it.Next()) {
        BRepAdaptor_Curve curve(TopoDS::Edge(it.Current()));BndLib_Add3dCurve::AddOptimal(curve,0,box);
    }
    for(TopExp_Explorer it(shape,TopAbs_VERTEX,TopAbs_EDGE);it.More();it.Next()) box.Add(BRep_Tool::Pnt(TopoDS::Vertex(it.Current())));
    if(box.IsVoid()) fail("boundingBox",invalid_input,"shape has no geometry");
    const auto out=limits(box);return {out.begin(),out.end()};
}
}
