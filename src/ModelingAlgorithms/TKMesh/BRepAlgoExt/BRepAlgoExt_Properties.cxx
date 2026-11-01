// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_Failure.hxx>
#include <BRepAlgoExt_TrimDomain.hxx>
#include <GeomAlgoExt_ParameterRectangle.hxx>
#include <BRepAlgoExt_Properties.hxx>

#include <BRepGProp.hxx>
#include <vector>
#include <GProp_GProps.hxx>
#include <BRepBndLib.hxx>
#include <Standard_Failure.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <Bnd_Box.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRep_Tool.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <BRepAdaptor_Curve2d.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <algorithm>

#include <cmath>

namespace BRepExt {

namespace {

// Boundary integration of the surface Jacobian. Split at spline spans before
// adaptive quadrature; the legacy adaptive SurfaceProperties path can report a
// small error while missing entire span contributions on trimmed extrusions.
class AreaIntegral {
    bool fixed;
    int evaluations=0;
    double noise=1e-12;
    double absoluteNoise=1e-12;
    bool converged=true;
    static constexpr double nodes[4]={0.1834346424956498,0.5255324099163290,0.7966664774136267,0.9602898564975363};
    static constexpr double weights[4]={0.3626837833783620,0.3137066458778873,0.2223810344533745,0.1012285362903763};
    template<class F> double gauss(F& f,double a,double b) {
        const double mid=(a+b)/2,half=(b-a)/2;double sum=0;
        for(int i=0;i<4;++i) sum+=weights[i]*(f(mid-half*nodes[i])+f(mid+half*nodes[i]));
        return sum*half;
    }
    template<class F> double integrate(F& f,double a,double b,double tol,int depth=0) {
        if(std::abs(b-a)<1e-14) return 0;
        if(fixed) return gauss(f,a,b);
        const double whole=gauss(f,a,b),mid=(a+b)/2,parts=gauss(f,a,mid)+gauss(f,mid,b);
        if(std::abs(parts-whole)<=std::max(tol,absoluteNoise+noise*std::abs(parts))) return parts;
        if(depth>=18 || evaluations>2000000) {converged=false;return parts;}
        return integrate(f,a,mid,tol/2,depth+1)+integrate(f,mid,b,tol/2,depth+1);
    }
public:
    explicit AreaIntegral(bool fixedRule=false) : fixed(fixedRule) {}
    bool didConverge() const {return converged;}
    double face(const TopoDS_Face& original,bool volume=false,const gp_Pnt& origin=gp_Pnt(),const double* areaEstimate=nullptr) {
        evaluations=0;
        const auto face=TopoDS::Face(original.Oriented(TopAbs_FORWARD));
        BRepAdaptor_Surface surface(face);
        double u0,u1,v0,v1;BRepTools::UVBounds(face,u0,u1,v0,v1);
        // The Boolean balance already measured this forward face with the same GProp call.
        // Reuse that exact estimate, retaining every tolerance and quadrature evaluation.
        GProp_GProps estimate;
        if(!areaEstimate && !fixed) BRepGProp::SurfaceProperties(face,estimate);
        const double area=areaEstimate ? *areaEstimate : estimate.Mass();
        double tolerance=std::max(1e-12,std::abs(area)*1e-8);
        noise=std::max(1e-14,std::abs(area)*1e-12);
        absoluteNoise=noise;
        if(volume) {
            Bnd_Box box;BRepBndLib::Add(face,box,false);
            double x0,y0,z0,x1,y1,z1;box.Get(x0,y0,z0,x1,y1,z1);
            const double radius=std::hypot(std::max(std::abs(x0-origin.X()),std::abs(x1-origin.X())),
                std::max(std::abs(y0-origin.Y()),std::abs(y1-origin.Y())),
                std::max(std::abs(z0-origin.Z()),std::abs(z1-origin.Z())));
            tolerance=std::max(1e-24,std::abs(area)*radius*1e-11);
            // This is a relative noise floor. Multiplying it by face area
            // made the accepted error grow as scale^5 for a volume integral.
            noise=1e-13;
            absoluteNoise=tolerance*1e-3;
        }
        const int count=surface.NbUIntervals(GeomAbs_CN);
        NCollection_Array1<double> knots(1,count+1);surface.UIntervals(knots,GeomAbs_CN);
        auto density=[&](double u,double v) {
            ++evaluations;
            gp_Pnt p;gp_Vec du,dv;surface.D1(u,v,p,du,dv);
            const auto normal=du.Crossed(dv);
            return volume ? gp_Vec(origin,p).Dot(normal)/3 : normal.Magnitude();
        };
        std::vector<double> breaks{u0};
        for(int i=1;i<=count+1;++i) if(knots(i)>u0 && knots(i)<u1) breaks.push_back(knots(i));
        breaks.push_back(u1);
        const bool rectangle=GeomAlgoExt::parameterRectangle(face);
        const bool planeArea=!volume && surface.GetType()==GeomAbs_Plane;
        if(rectangle && planeArea) return (u1-u0)*(v1-v0);
        // On an extrusion, a cylinder among them, the density is independent of V: the normal
        // is the same along a ruling, and so is the part of the position along it.
        const bool extrusion=surface.GetType()==GeomAbs_SurfaceOfExtrusion || surface.GetType()==GeomAbs_Cylinder;
        if(rectangle) {
            const int nv=surface.NbVIntervals(GeomAbs_CN);
            NCollection_Array1<double> vk(1,nv+1);surface.VIntervals(vk,GeomAbs_CN);
            std::vector<double> vb{v0};
            for(int i=1;i<=nv+1;++i) if(vk(i)>v0 && vk(i)<v1) vb.push_back(vk(i));
            vb.push_back(v1);
            double total=0;
            for(size_t i=1;i<breaks.size();++i) for(size_t j=1;j<vb.size();++j) {
                if(extrusion) {
                    // The integral over U, the same at every V, times the height.
                    const double height=vb[j]-vb[j-1],v=(vb[j-1]+vb[j])/2;
                    auto inner=[&](double u){return density(u,v);};
                    total+=height*integrate(inner,breaks[i-1],breaks[i],tolerance/(10*breaks.size()*vb.size()*std::max(1.0,height)));
                    continue;
                }
                auto outer=[&](double v) {
                    auto inner=[&](double u){return density(u,v);};
                    return integrate(inner,breaks[i-1],breaks[i],tolerance/(100*breaks.size()*vb.size()));
                };
                total+=integrate(outer,vb[j-1],vb[j],tolerance/(10*breaks.size()*vb.size()));
            }
            return volume && original.Orientation()==TopAbs_REVERSED ? -total : total;
        }
        // So on an extrusion the integral of complete U spans is cached once; each boundary
        // evaluation integrates only its last partial span, avoiding the former spans-squared
        // cost.
        // On a non-rational ruled B-spline, p=P(u)+v*D(u). Its volume density
        // p.(du x dv)/3 is affine in V (the terms with two D factors vanish).
        // Cache the two endpoint primitives instead of reintegrating every U
        // span at each contour sample. Area density is not affine in V.
        const bool ruledVolume=volume && surface.GetType()==GeomAbs_BSplineSurface &&
            surface.VDegree()==1 && !surface.IsURational() && !surface.IsVRational() &&
            surface.NbVIntervals(GeomAbs_CN)==1 && v1-v0>1e-14;
        std::vector<double> prefix(breaks.size(),0);
        std::vector<double> other(breaks.size(),0);
        if(extrusion || ruledVolume) {
            auto f=[&](double u){return density(u,ruledVolume?v0:(v0+v1)/2);};
            auto g=[&](double u){return density(u,v1);};
            for(size_t i=1;i<breaks.size();++i)
                prefix[i]=prefix[i-1]+integrate(f,breaks[i-1],breaks[i],tolerance/(100*breaks.size()));
            if(ruledVolume) for(size_t i=1;i<breaks.size();++i)
                other[i]=other[i-1]+integrate(g,breaks[i-1],breaks[i],tolerance/(100*breaks.size()));
        }
        double total=0;
        for(TopExp_Explorer wires(face,TopAbs_WIRE);wires.More();wires.Next())
        for(BRepTools_WireExplorer walk(TopoDS::Wire(wires.Current()),face);walk.More();walk.Next()) {
            BRepAdaptor_Curve2d boundary(walk.Current(),face);
            auto contour=[&](double t) {
                gp_Pnt2d uv;gp_Vec2d derivative;boundary.D1(t,uv,derivative);
                if(std::abs(derivative.Y())<1e-18) return 0.0;
                // A plane's orthonormal parameter frame has unit Jacobian.
                // Its primitive is exact; only the boundary needs quadrature.
                if(planeArea) return (uv.X()-u0)*derivative.Y();
                auto inner=[&](double u) {return density(u,uv.Y());};
                double sum=0,left=u0;
                const double right=std::clamp(uv.X(),u0,u1);
                const double tol=tolerance/(100*std::max(1.0,v1-v0)*(count+1));
                if(extrusion || ruledVolume) {
                    const auto at=std::min(breaks.size()-2,size_t(std::upper_bound(breaks.begin(),breaks.end(),right)-breaks.begin()-1));
                    sum=prefix[at];left=breaks[at];
                    if(ruledVolume) sum+=(other[at]-prefix[at])*(uv.Y()-v0)/(v1-v0);
                } else for(int i=1;i<=count+1;++i) if(knots(i)>left && knots(i)<right) {
                    sum+=integrate(inner,left,knots(i),tol);left=knots(i);
                }
                sum+=integrate(inner,left,right,tol);
                return sum*derivative.Y();
            };
            const int n=boundary.NbIntervals(GeomAbs_CN);
            NCollection_Array1<double> breaks(1,n+1);boundary.Intervals(breaks,GeomAbs_CN);
            const double sign=walk.Current().Orientation()==TopAbs_REVERSED?-1:1;
            for(int i=1;i<=n;++i) total+=sign*integrate(contour,breaks(i),breaks(i+1),tolerance/(10*n));
        }
        return volume ? (original.Orientation()==TopAbs_REVERSED ? -total:total) : std::abs(total);
    }
};

// Faces of a support other than a plane, cylinder, cone, sphere or torus. The boundary
// integral above sweeps the support between the trim and the side of its parameter
// rectangle, which on an extrapolated or folded spline or offset is far larger than the face
// and cancels out of it only approximately, and crosses every span of the support at each of
// its samples; the facade's fixed rule misplaces the centroid of a rational one.
bool trimmedSupport(const TopoDS_Face& face) {
    TopLoc_Location location;
    return GeomAdaptor_Surface(BRep_Tool::Surface(face,location)).GetType()>GeomAbs_Torus;
}
// Integrated over the trim only; where the trim cannot be read, the boundary integral and the
// facade's centroid stand in, as not converged.
TrimProperties trimmedFace(const TopoDS_Face& face, const gp_Pnt* volumeOrigin=nullptr) {
    try {
        return trimmedProperties(face,volumeOrigin);
    } catch(const Failure&) {
    } catch(const Standard_Failure&) {
    }
    TrimProperties out;AreaIntegral integral;out.area=integral.face(face);
    if(volumeOrigin) out.volume=integral.face(TopoDS::Face(face.Oriented(TopAbs_FORWARD)),true,*volumeOrigin);
    GProp_GProps props;BRepGProp::SurfaceProperties(face,props);out.center=props.CentreOfMass();
    return out;
}

} // namespace

/// What one face adds to its body's area, or to its volume about `origin` (signed by the
/// face's orientation). The Boolean check adds these up over the cells of a partition.
double faceMeasure(const TopoDS_Face& face, bool surface, const gp_Pnt& origin, bool* converged, const double* areaEstimate) {
    double total=0;
    {
        BRepAdaptor_Surface adaptor(face);
        if(adaptor.GetType()==GeomAbs_Plane) {
            // On a plane, (p-origin).normal is constant. Integrate only its
            // trimmed area rather than a nested adaptive volume integral.
            // Besides avoiding pathological near-tangent trim subdivisions,
            // this keeps small removed volumes accurate on large flat bodies.
            // Integrate each pcurve span. OCCT's nested adaptive surface rule can
            // report a small error while missing oscillations in a spline trim.
            AreaIntegral integral;
            const double area=integral.face(face);
            const auto plane=adaptor.Plane();const auto position=plane.Position();
            const auto normal=position.XDirection().Crossed(position.YDirection());
            const double sign=face.Orientation()==TopAbs_REVERSED ? -1 : 1;
            total+=surface ? area : sign*area*gp_Vec(origin,plane.Location()).Dot(gp_Vec(normal))/3;
            if(converged) *converged &= integral.didConverge();
        } else if(surface && adaptor.GetType()<=GeomAbs_Torus) {
            GProp_GProps props;BRepGProp::SurfaceProperties(face,props,1e-9);
            total+=props.Mass();
        } else if(surface && trimmedSupport(face)) {
            const auto measured=trimmedFace(face);
            total+=measured.area;
            if(converged) *converged &= measured.converged;
        } else if(!surface && trimmedSupport(face) && !GeomAlgoExt::parameterRectangle(face)) {
            // Integrate only material inside the trim. A boundary primitive over
            // the whole support reintegrates every U span at every trim sample,
            // and can spend minutes cancelling extrapolated spline regions.
            const auto measured=trimmedFace(face,&origin);
            total+=(face.Orientation()==TopAbs_REVERSED ? -1 : 1)*measured.volume;
            if(converged) *converged &= measured.converged;
        } else {
            AreaIntegral integral;
            total+=integral.face(face,!surface,origin,areaEstimate);
            if(converged) *converged &= integral.didConverge();
        }
    }
    return total;
}

/// The area of a face. On a rectangle of the parameters of a plane, cylinder, cone, sphere or
/// torus it has a closed form; of any other face it is OpenCascade's Gauss rule, as before.
double faceArea(const TopoDS_Face& original) {
    const auto face=TopoDS::Face(original.Oriented(TopAbs_FORWARD));
    const BRepAdaptor_Surface surface(face);
    if(surface.GetType()<=GeomAbs_Torus && GeomAlgoExt::parameterRectangle(face)) {
        double u0,u1,v0,v1;BRepTools::UVBounds(face,u0,u1,v0,v1);
        const double du=u1-u0,dv=v1-v0;
        switch(surface.GetType()) {
        case GeomAbs_Plane: return du*dv;
        case GeomAbs_Cylinder: return surface.Cylinder().Radius()*du*dv;
        case GeomAbs_Cone: {
            // The distance from the axis is R + v sin(a), along a generatrix of unit speed.
            const auto cone=surface.Cone();const double s=std::sin(cone.SemiAngle()),r=cone.RefRadius();
            const auto ring=[&](double v){return r*v+s*v*v/2;};
            if((r+s*v0)*(r+s*v1)<0) break; // across the apex
            return std::abs(du*(ring(v1)-ring(v0)));
        }
        case GeomAbs_Sphere: {
            const double r=surface.Sphere().Radius();
            return r*r*du*(std::sin(v1)-std::sin(v0));
        }
        case GeomAbs_Torus: {
            const auto torus=surface.Torus();const double big=torus.MajorRadius(),r=torus.MinorRadius();
            if(big<r) break; // a spindle torus folds over its axis
            return r*du*(big*dv+r*(std::sin(v1)-std::sin(v0)));
        }
        default: break;
        }
    }
    GProp_GProps props;
    BRepGProp::SurfaceProperties(face,props);
    return props.Mass();
}

// Coarse material bracket only: eight Gauss nodes in each direction per span/cell.
// Exact property operations and Boolean balance keep their adaptive convergence checks.
double faceMeasureBracket(const TopoDS_Face& face, const gp_Pnt& origin) {
    if (BRepAdaptor_Surface(face).GetType() == GeomAbs_Plane)
        return faceMeasure(face, false, origin, nullptr, nullptr);
    AreaIntegral integral(true);
    return integral.face(face, true, origin);
}

double measure(const TopoDS_Shape& shape, bool surface, bool* converged) {
    if(converged) *converged=true;
    Bnd_Box box; BRepBndLib::Add(shape,box,false);
    if(box.IsVoid()) return 0;
    if(box.IsOpen()) fail(surface ? "surfaceArea":"volume",invalid_input,"finite body required");
    double x0,y0,z0,x1,y1,z1;box.Get(x0,y0,z0,x1,y1,z1);
    const gp_Pnt origin((x0+x1)/2,(y0+y1)/2,(z0+z1)/2);
    double total=0;
    for(TopExp_Explorer it(shape,TopAbs_FACE);it.More();it.Next())
        total+=faceMeasure(TopoDS::Face(it.Current()),surface,origin,converged);
    return total;
}

std::vector<double> surfaceCentre(const TopoDS_Shape& shape) {
    double area=0;gp_XYZ moment(0,0,0);
    for(TopExp_Explorer it(shape,TopAbs_FACE);it.More();it.Next()) {
        const auto face=TopoDS::Face(it.Current());TrimProperties measured;
        if(trimmedSupport(face)) measured=trimmedFace(face);
        else {GProp_GProps props;BRepGProp::SurfaceProperties(face,props);measured.area=props.Mass();measured.center=props.CentreOfMass();}
        area+=measured.area;moment+=measured.center.XYZ()*measured.area;
    }
    if(!(area>0) || !std::isfinite(area)) {
        // With no area to weigh, the centre of mass of the whole shape's surface.
        GProp_GProps props;BRepGProp::SurfaceProperties(shape,props);const auto center=props.CentreOfMass();
        return {center.X(),center.Y(),center.Z()};
    }
    const auto center=moment/area;
    return {center.X(),center.Y(),center.Z()};
}

} // namespace BRepExt
