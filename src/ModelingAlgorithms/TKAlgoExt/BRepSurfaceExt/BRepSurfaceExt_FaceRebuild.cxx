// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_History.hxx>
#include <BRepSurfaceExt_Tools.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <NCollection_IndexedMap.hxx>
#include <Standard_Failure.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>

#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Builder.hxx>
#include <BRepLib.hxx>
#include <Geom2dConvert.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <TopoDS_Wire.hxx>
#include <TopExp_Explorer.hxx>
#include <BRep_Tool.hxx>
#include <Geom_BSplineCurve.hxx>
#include <NCollection_Array2.hxx>
#include <math_Gauss.hxx>
#include <math_Matrix.hxx>
#include <math_Vector.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_BezierSurface.hxx>
#include <GeomConvert.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <TopoDS.hxx>
#include <gp_Pln.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Cone.hxx>
#include <gp_Circ.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <TopExp.hxx>
#include <NCollection_DataMap.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <cmath>
#include <algorithm>
namespace loft {
namespace {
TopoDS_Face rectangular(OcctKernel &k,uint32_t id,double tolerance) {
    if(get(k,id).ShapeType()!=TopAbs_FACE || !std::isfinite(tolerance) || tolerance<=0)
        fail("faceSupport",invalid_input,"face and positive tolerance required");
    const auto face=TopoDS::Face(get(k,id));
    double a,b,c,d;BRepTools::UVBounds(face,a,b,c,d);
    const auto rectangle=BRepBuilderAPI_MakeFace(BRep_Tool::Surface(face),a,b,c,d,tolerance).Face();
    const double area=surfaceArea(k,id),envelope=surfaceArea(k,store(k,rectangle));
    if(std::abs(area-envelope)>tolerance*std::max(1.0,area))
        fail("faceSupport",unsupported_geometry,"requires a complete rectangular UV face; trimmed holes and nonrectangular boundaries are unsupported");
    return face;
}
EvolutionData result(OcctKernel &k,const TopoDS_Shape &shape) {
    if(!BRepCheck_Analyzer(shape).IsValid()) fail("faceSupport",no_solution,"invalid surface result");
    return stored(k, shape);
}
}
namespace {
/// The face's support interpolated over the face's UV box, as a B-spline on [0,1]x[0,1], and
/// checked against the support on an independent grid. `seams` lets a closed face through: its
/// two seam rows are interpolated from the same points, so the fit closes as the support does.
struct Fit { Handle(Geom_BSplineSurface) surface; double a,b,c,d; };
Fit fitSupport(const char* op,const TopoDS_Face& face,int uDegree,int vDegree,int spans,double tolerance,bool seams) {
    if(uDegree<3 || vDegree<3 || uDegree>8 || vDegree>8 || spans<1 || spans>100)
        fail(op,invalid_input,"degree limits 3..8 and segment limit 1..100 required");
    double a,b,c,d;BRepTools::UVBounds(face,a,b,c,d);
    if(!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c) || !std::isfinite(d) || b<=a || d<=c)
        fail(op,unsupported_geometry,"finite nonempty trim bounds required");
    for(TopExp_Explorer it(face,TopAbs_EDGE);it.More();it.Next()) {
        const auto edge=TopoDS::Edge(it.Current());
        if(BRep_Tool::Degenerated(edge) || (!seams && BRep_Tool::IsClosed(edge,face)))
            fail(op,unsupported_geometry,seams?"degenerate trim edges are unsupported":"seam and degenerate trim edges are unsupported");
    }
    BRepAdaptor_Surface source(face);
    // Retain the source's breakpoints and continuity when the requested span
    // count can contain them. Raising a cubic to degree five requires repeated
    // knots at its C2 joins; uniform simple knots cannot represent that surface.
    auto basis=[&](int degree,bool alongU) {
        const double low=alongU?a:c, high=alongU?b:d;
        const int count=alongU?source.NbUIntervals(GeomAbs_CN):source.NbVIntervals(GeomAbs_CN);
        NCollection_Array1<double> intervals(1,count+1);
        if(alongU) source.UIntervals(intervals,GeomAbs_CN); else source.VIntervals(intervals,GeomAbs_CN);
        std::vector<std::pair<double,int>> values{{0,degree+1},{1,degree+1}};
        if(count<=spans) for(int i=2;i<=count;++i) {
            const double t=(intervals(i)-low)/(high-low);
            if(t<=loft::tolerance::parameter || t>=1-loft::tolerance::parameter) continue;
            int continuity=2;
            if(source.GetType()==GeomAbs_BSplineSurface) {
                const auto bs=source.BSpline();
                const int n=alongU?bs->NbUKnots():bs->NbVKnots();
                for(int j=1;j<=n;++j) if(std::abs((alongU?bs->UKnot(j):bs->VKnot(j))-intervals(i))<1e-10)
                    continuity=(alongU?bs->UDegree()-bs->UMultiplicity(j):bs->VDegree()-bs->VMultiplicity(j));
            } else if(alongU && source.GetType()==GeomAbs_SurfaceOfExtrusion && source.BasisCurve()->GetType()==GeomAbs_BSplineCurve) {
                const auto bs=source.BasisCurve()->BSpline();
                for(int j=1;j<=bs->NbKnots();++j) if(std::abs(bs->Knot(j)-intervals(i))<1e-10)
                    continuity=bs->Degree()-bs->Multiplicity(j);
            }
            values.push_back({t,std::max(1,degree-continuity)});
        }
        std::sort(values.begin(),values.end());
        while(int(values.size())<spans+1) {
            size_t at=1;
            for(size_t i=2;i<values.size();++i)
                if(values[i].first-values[i-1].first>values[at].first-values[at-1].first) at=i;
            values.insert(values.begin()+at,{(values[at-1].first+values[at].first)/2,1});
        }
        int n=-degree-1;for(const auto& entry:values) n+=entry.second;
        NCollection_Array1<gp_Pnt> poles(1,n);poles.Init(gp_Pnt());
        NCollection_Array1<double> knots(1,values.size());NCollection_Array1<int> mult(1,values.size());
        for(size_t i=0;i<values.size();++i) {knots(i+1)=values[i].first;mult(i+1)=values[i].second;}
        return Handle(Geom_BSplineCurve)(new Geom_BSplineCurve(poles,knots,mult,degree));
    };
    auto u=basis(uDegree,true),v=basis(vDegree,false);
    auto sites=[](const Handle(Geom_BSplineCurve)& curve) {
        std::vector<double> expanded,sites;
        for(int i=1;i<=curve->NbKnots();++i)
            for(int j=0;j<curve->Multiplicity(i);++j) expanded.push_back(curve->Knot(i));
        for(int i=0;i<curve->NbPoles();++i) {
            double t=0;for(int j=1;j<=curve->Degree();++j) t+=expanded[i+j]/curve->Degree();
            sites.push_back(t);
        }
        return sites;
    };
    const auto us=sites(u),vs=sites(v);
    auto matrix=[](const Handle(Geom_BSplineCurve)& curve,const std::vector<double>& sites) {
        const int n=curve->NbPoles(); math_Matrix m(1,n,1,n);
        for(int col=1;col<=n;++col) {
            curve->SetPole(col,gp_Pnt(1,0,0));
            for(int row=1;row<=n;++row) m(row,col)=curve->Value(sites[row-1]).X();
            curve->SetPole(col,gp_Pnt());
        }
        return m;
    };
    math_Gauss solveU(matrix(u,us)),solveV(matrix(v,vs));
    if(!solveU.IsDone() || !solveV.IsDone()) fail(op,no_solution,"singular surface fit");
    const int nu=u->NbPoles(),nv=v->NbPoles();
    NCollection_Array2<gp_Pnt> poles(1,nu,1,nv),middle(1,nu,1,nv);
    for(int axis=1;axis<=3;++axis) {
        for(int j=1;j<=nv;++j) {
            math_Vector rhs(1,nu),solution(1,nu);
            for(int i=1;i<=nu;++i) rhs(i)=source.Value(a+(b-a)*us[i-1],c+(d-c)*vs[j-1]).Coord(axis);
            solveU.Solve(rhs,solution);
            for(int i=1;i<=nu;++i) middle(i,j).SetCoord(axis,solution(i));
        }
        for(int i=1;i<=nu;++i) {
            math_Vector rhs(1,nv),solution(1,nv);
            for(int j=1;j<=nv;++j) rhs(j)=middle(i,j).Coord(axis);
            solveV.Solve(rhs,solution);
            for(int j=1;j<=nv;++j) poles(i,j).SetCoord(axis,solution(j));
        }
    }
    Handle(Geom_BSplineSurface) fitted=new Geom_BSplineSurface(poles,u->Knots(),v->Knots(),u->Multiplicities(),v->Multiplicities(),uDegree,vDegree);
    // Independent sites between the interpolation sites. This is a sampled
    // geometric tolerance check, not a certified bound for arbitrary surfaces.
    const int checks=std::max(64,spans*8);
    for(int i=0;i<=checks;++i) for(int j=0;j<=checks;++j) {
        const double x=double(i)/checks,y=double(j)/checks;
        if(fitted->Value(x,y).Distance(source.Value(a+(b-a)*x,c+(d-c)*y))>tolerance)
            fail(op,no_solution,"surface fit did not meet tolerance at the requested degree and segment count");
    }
    return {fitted,a,b,c,d};
}
}
EvolutionData rebuildFace(OcctKernel &k,uint32_t id,int uDegree,int vDegree,int spans,double tolerance) {
    if(get(k,id).ShapeType()!=TopAbs_FACE || !std::isfinite(tolerance) || tolerance<=0)
        fail("rebuildFace",invalid_input,"face and positive tolerance required");
    const auto face=TopoDS::Face(get(k,id));
    const auto fit=fitSupport("rebuildFace",face,uDegree,vDegree,spans,tolerance,false);
    const auto& fitted=fit.surface;
    const double a=fit.a,b=fit.b,c=fit.c,d=fit.d;
    BRepAdaptor_Surface source(face);
    // Transfer the original trim domain, including inner wires, by its exact
    // affine UV map. Reprojecting 3D edges would change the requested boundary.
    BRep_Builder builder;TopoDS_Face rebuilt;
    builder.MakeFace(rebuilt,fitted,loft::tolerance::length);
    builder.NaturalRestriction(rebuilt,false);
    NCollection_DataMap<TopoDS_Shape,TopoDS_Shape,TopTools_ShapeMapHasher> vertices;
    const auto forward=TopoDS::Face(face.Oriented(TopAbs_FORWARD));
    for(TopExp_Explorer wi(forward,TopAbs_WIRE);wi.More();wi.Next()) {
        BRepBuilderAPI_MakeWire wire;
        for(BRepTools_WireExplorer ei(TopoDS::Wire(wi.Current()),forward);ei.More();ei.Next()) {
            const auto edge=ei.Current();
            double first,last;const auto pc=BRep_Tool::CurveOnSurface(edge,forward,first,last);
            if(pc.IsNull()) fail("rebuildFace",no_solution,"trim edge has no parameter curve");
            // Check world-space boundary sites separately from the support grid:
            // a narrow trim can lie entirely between its independent grid sites.
            for(int i=0;i<=128;++i) {
                const auto uv=pc->Value(first+(last-first)*i/128);
                if(fitted->Value((uv.X()-a)/(b-a),(uv.Y()-c)/(d-c)).Distance(source.Value(uv.X(),uv.Y()))>tolerance)
                    fail("rebuildFace",no_solution,"trim boundary fit did not meet tolerance");
            }
            auto curve=Handle(Geom2d_BSplineCurve)::DownCast(
                Geom2dConvert::CurveToBSplineCurve(new Geom2d_TrimmedCurve(pc,first,last))->Copy());
            for(int i=1;i<=curve->NbPoles();++i) {
                const auto p=curve->Pole(i);
                curve->SetPole(i,gp_Pnt2d((p.X()-a)/(b-a),(p.Y()-c)/(d-c)));
            }
            // Preserve the source's vertex sharing, especially a closed trim's
            // single vertex. Its two p-curve ends may differ within the input's
            // existing tolerance; guessing connectivity from proximity loses it.
            TopoDS_Vertex v0,v1;
            TopExp::Vertices(TopoDS::Edge(edge.Oriented(TopAbs_FORWARD)),v0,v1);
            auto vertex=[&](const TopoDS_Vertex& original,double t) {
                const auto uv=curve->Value(t);const auto point=fitted->Value(uv.X(),uv.Y());
                const double inherited=std::max(loft::tolerance::length,BRep_Tool::Tolerance(original));
                if(!vertices.IsBound(original)) {
                    auto made=BRepBuilderAPI_MakeVertex(point).Vertex();
                    builder.UpdateVertex(made,inherited);vertices.Bind(original,made);
                }
                auto made=TopoDS::Vertex(vertices.Find(original));
                const double gap=BRep_Tool::Pnt(made).Distance(point);
                if(gap>tolerance+inherited)
                    fail("rebuildFace",no_solution,"trim endpoints exceed the fit tolerance");
                builder.UpdateVertex(made,std::max(inherited,gap+loft::tolerance::length));
                return made;
            };
            const auto firstVertex=vertex(v0,curve->FirstParameter());
            const auto lastVertex=vertex(v1,curve->LastParameter());
            BRepBuilderAPI_MakeEdge make(curve,fitted,firstVertex,lastVertex,curve->FirstParameter(),curve->LastParameter());
            if(!make.IsDone()) fail("rebuildFace",no_solution,"could not rebuild trim edge");
            auto made=make.Edge();made.Orientation(edge.Orientation());wire.Add(made);
        }
        if(!wire.IsDone()) fail("rebuildFace",no_solution,"could not rebuild trim wire");
        builder.Add(rebuilt,wire.Wire());
    }
    if(!BRepLib::BuildCurves3d(rebuilt,loft::tolerance::length))
        fail("rebuildFace",no_solution,"could not construct trim curves");
    rebuilt.Orientation(face.Orientation());
    return result(k,rebuilt);
}
EvolutionData rebuildBodyFace(OcctKernel &k,uint32_t bodyId,uint32_t faceId,int uDegree,int vDegree,int spans,
                              double tolerance,const std::vector<int>& hashes,int bound) {
    const char* op="rebuildFace";
    if(!std::isfinite(tolerance) || tolerance<=0) fail(op,invalid_input,"positive tolerance required");
    const auto& input=get(k,bodyId);
    NCollection_IndexedMap<TopoDS_Shape,TopTools_ShapeMapHasher> members;
    TopExp::MapShapes(input,TopAbs_FACE,members);
    const auto& original=get(k,faceId);
    if(original.ShapeType()!=TopAbs_FACE || !members.Contains(original))
        fail(op,invalid_input,"select a face of the body to rebuild");
    FaceHistory history(op,hashes,bound);
    BRepBuilderAPI_Copy copy(input,true,false);
    const auto output=copy.Shape();
    auto face=TopoDS::Face(copy.ModifiedShape(original));
    const auto fit=fitSupport(op,face,uDegree,vDegree,spans,tolerance,true);
    // The fit is laid on the face's own UV box, so every trim keeps its p-curve, seams included,
    // and the edges stay shared with the neighbours: only their tolerance grows by the fit's
    // distance from them.
    auto domain=[](const NCollection_Array1<double>& knots,double low,double high) {
        NCollection_Array1<double> out(knots.Lower(),knots.Upper());
        for(int i=knots.Lower();i<=knots.Upper();++i) out(i)=low+(high-low)*knots(i);
        return out;
    };
    const auto& fitted=fit.surface;
    Handle(Geom_BSplineSurface) support=new Geom_BSplineSurface(fitted->Poles(),
        domain(fitted->UKnots(),fit.a,fit.b),domain(fitted->VKnots(),fit.c,fit.d),
        fitted->UMultiplicities(),fitted->VMultiplicities(),fitted->UDegree(),fitted->VDegree());
    struct Curve { TopoDS_Edge edge; Handle(Geom2d_Curve) curve,second; double first,last; };
    std::vector<Curve> curves;
    const auto forward=TopoDS::Face(face.Oriented(TopAbs_FORWARD));
    NCollection_IndexedMap<TopoDS_Shape,TopTools_ShapeMapHasher> faceEdges;
    TopExp::MapShapes(forward,TopAbs_EDGE,faceEdges);
    for(const auto& shape:faceEdges) {
        const auto edge=TopoDS::Edge(shape.Oriented(TopAbs_FORWARD));double first,last;
        const auto curve=BRep_Tool::CurveOnSurface(edge,forward,first,last);
        if(curve.IsNull()) fail(op,no_solution,"face edge has no parameter curve");
        Handle(Geom2d_Curve) second;
        if(BRep_Tool::IsClosed(edge,forward)) {
            double s,t;second=BRep_Tool::CurveOnSurface(TopoDS::Edge(edge.Reversed()),forward,s,t);
            if(second.IsNull()) fail(op,no_solution,"seam has no second parameter curve");
        }
        curves.push_back({edge,curve,second,first,last});
    }
    BRep_Builder builder;
    builder.UpdateFace(face,support,TopLoc_Location(),BRep_Tool::Tolerance(face));
    builder.NaturalRestriction(face,false);
    for(const auto& curve:curves) {
        if(curve.second.IsNull()) builder.UpdateEdge(curve.edge,curve.curve,face,BRep_Tool::Tolerance(curve.edge));
        else builder.UpdateEdge(curve.edge,curve.curve,curve.second,face,BRep_Tool::Tolerance(curve.edge));
        builder.Range(curve.edge,face,curve.first,curve.last);
        // How far the edge is now from the rebuilt support, along its whole length.
        BRepAdaptor_Curve onSpace(curve.edge),onFace(curve.edge,face);
        double gap=0;
        for(int i=0;i<=128;++i) {
            const double t=curve.first+(curve.last-curve.first)*i/128;
            gap=std::max(gap,onSpace.Value(t).Distance(onFace.Value(t)));
        }
        const double edgeTolerance=BRep_Tool::Tolerance(curve.edge);
        if(gap>tolerance+edgeTolerance) fail(op,no_solution,"trim boundary fit did not meet tolerance");
        if(gap>edgeTolerance) builder.UpdateEdge(curve.edge,gap+loft::tolerance::length);
    }
    // A vertex covers its edges' ends, which may now sit up to their new tolerance away.
    for(TopExp_Explorer it(forward,TopAbs_VERTEX);it.More();it.Next()) {
        const auto vertex=TopoDS::Vertex(it.Current());
        double need=BRep_Tool::Tolerance(vertex);
        for(const auto& curve:curves) {
            TopoDS_Vertex v0,v1;TopExp::Vertices(curve.edge,v0,v1);
            if(v0.IsSame(vertex) || v1.IsSame(vertex)) need=std::max(need,BRep_Tool::Tolerance(curve.edge));
        }
        builder.UpdateVertex(vertex,need);
    }
    if(!BRepCheck_Analyzer(output).IsValid()) fail(op,no_solution,"rebuilding the face produced invalid topology");
    for(const auto& each:members) if(history.tracks(each)) {
        NCollection_List<TopoDS_Shape> successor;successor.Append(copy.ModifiedShape(each));
        history.modified(each,successor);
    }
    return history.finish(store(k,output));
}
std::vector<double> faceSpline(OcctKernel& k,uint32_t id) {
    const char* op="faceSpline";
    const auto& shape=get(k,id);
    if(shape.ShapeType()!=TopAbs_FACE) fail(op,invalid_input,"face required");
    const auto face=TopoDS::Face(shape);
    TopLoc_Location location;auto basis=BRep_Tool::Surface(face,location);
    while(auto trim=Handle(Geom_RectangularTrimmedSurface)::DownCast(basis)) basis=trim->BasisSurface();
    auto spline=Handle(Geom_BSplineSurface)::DownCast(basis);
    if(spline.IsNull() && !Handle(Geom_BezierSurface)::DownCast(basis).IsNull())
        spline=GeomConvert::SurfaceToBSplineSurface(basis);
    if(spline.IsNull()) fail(op,unsupported_geometry,"B-spline or Bezier surface required");
    spline=Handle(Geom_BSplineSurface)::DownCast(spline->Transformed(location.Transformation()));
    const int nu=spline->NbUPoles(),nv=spline->NbVPoles(),ku=spline->NbUKnots(),kv=spline->NbVKnots();
    std::vector<double> out{double(spline->UDegree()),double(spline->VDegree()),double(ku),double(kv),double(nu),double(nv)};
    for(int i=1;i<=ku;++i) out.push_back(spline->UKnot(i));
    for(int i=1;i<=ku;++i) out.push_back(spline->UMultiplicity(i));
    for(int i=1;i<=kv;++i) out.push_back(spline->VKnot(i));
    for(int i=1;i<=kv;++i) out.push_back(spline->VMultiplicity(i));
    for(int i=1;i<=nu;++i) for(int j=1;j<=nv;++j) {
        const auto p=spline->Pole(i,j);out.insert(out.end(),{p.X(),p.Y(),p.Z(),spline->Weight(i,j)});
    }
    return out;
}
occ::handle<Geom_BSplineSurface> splineFromData(const char* op,const std::vector<double>& data) {
    auto at=[&](size_t i) {
        if(i>=data.size() || !std::isfinite(data[i])) fail(op,invalid_input,"incomplete or non-finite spline data");
        return data[i];
    };
    auto whole=[&](size_t i,int low,int high) {
        const double value=at(i);
        if(value!=std::floor(value) || value<low || value>high) fail(op,invalid_input,"spline counts out of range");
        return int(value);
    };
    const int du=whole(0,1,25),dv=whole(1,1,25),ku=whole(2,2,10000),kv=whole(3,2,10000),nu=whole(4,2,10000),nv=whole(5,2,10000);
    size_t next=6;
    auto knots=[&](int count,NCollection_Array1<double>& values,NCollection_Array1<int>& mults) {
        for(int i=1;i<=count;++i) values(i)=at(next++);
        for(int i=1;i<=count;++i) mults(i)=whole(next++,1,26);
    };
    NCollection_Array1<double> uk(1,ku),vk(1,kv);NCollection_Array1<int> um(1,ku),vm(1,kv);
    knots(ku,uk,um);knots(kv,vk,vm);
    NCollection_Array2<gp_Pnt> poles(1,nu,1,nv);NCollection_Array2<double> weights(1,nu,1,nv);
    bool rational=false;
    for(int i=1;i<=nu;++i) for(int j=1;j<=nv;++j) {
        poles(i,j)=gp_Pnt(at(next),at(next+1),at(next+2));weights(i,j)=at(next+3);next+=4;
        if(weights(i,j)<=0) fail(op,invalid_input,"weights must be positive");
        rational=rational || std::abs(weights(i,j)-1)>1e-12;
    }
    if(next!=data.size()) fail(op,invalid_input,"spline data has the wrong length");
    try {
        return rational?new Geom_BSplineSurface(poles,weights,uk,vk,um,vm,du,dv)
                       :new Geom_BSplineSurface(poles,uk,vk,um,vm,du,dv);
    } catch(const Standard_Failure&) { fail(op,invalid_input,"inconsistent knots, multiplicities and poles"); }
}
EvolutionData splineFace(OcctKernel& k,const std::vector<double>& data) {
    BRepBuilderAPI_MakeFace make(splineFromData("splineFace",data),loft::tolerance::length);
    if(!make.IsDone()) fail("splineFace",no_solution,"could not make a face on the spline");
    return stored(k, make.Face());
}
std::vector<int> surfaceSplineInfo(OcctKernel& k,uint32_t id) {
    const auto& shape=get(k,id);
    if(shape.ShapeType()!=TopAbs_FACE) fail("surfaceSplineInfo",invalid_input,"face required");
    BRepAdaptor_Surface surface(TopoDS::Face(shape));
    if(surface.GetType()!=GeomAbs_BSplineSurface) fail("surfaceSplineInfo",unsupported_geometry,"B-spline surface required");
    const auto spline=surface.BSpline();
    return {spline->UDegree(),spline->VDegree(),spline->NbUKnots()-1,spline->NbVKnots()-1};
}
EvolutionData unwrapFace(OcctKernel &k,uint32_t id,double tolerance) {
    const auto face=TopoDS::Face(get(k,id));
    BRepAdaptor_Surface surface(face);
    if(surface.GetType()==GeomAbs_Plane) {
        gp_Trsf move;
        auto plane=surface.Plane();
        auto frame=plane.Position();
        if(face.Orientation()==TopAbs_REVERSED) frame.ZReverse();
        // Keep a direct frame: flipping only Z would reflect the developed face.
        if(!frame.Direct()) frame.YReverse();
        move.SetTransformation(frame);
        auto flat=BRepBuilderAPI_Transform(face,move,true).Shape();
        GProp_GProps properties;BRepGProp::SurfaceProperties(flat,properties);
        const auto centre=properties.CentreOfMass();gp_Trsf centreMove;
        centreMove.SetTranslation(gp_Vec(-centre.X(),-centre.Y(),-centre.Z()));
        return result(k,BRepBuilderAPI_Transform(flat,centreMove,true).Shape());
    }
    rectangular(k,id,tolerance);
    if(surface.GetType()==GeomAbs_Cone) {
        const auto cone=surface.Cone();const double sine=std::sin(cone.SemiAngle());
        const double angle=(surface.LastUParameter()-surface.FirstUParameter())*std::abs(sine);
        const double r0=cone.RefRadius()+surface.FirstVParameter()*sine;
        const double r1=cone.RefRadius()+surface.LastVParameter()*sine;
        if(r0<-tolerance || r1<-tolerance || angle<=loft::tolerance::angle)
            fail("unwrapFace",unsupported_geometry,"cone development cannot cross the apex");
        const double inner=std::max(0.0,std::min(r0,r1))/std::abs(sine);
        const double outer=std::max(r0,r1)/std::abs(sine);
        auto point=[](double radius,double t){return gp_Pnt(radius*std::cos(t),radius*std::sin(t),0);};
        BRepBuilderAPI_MakeWire wire;
        wire.Add(BRepBuilderAPI_MakeEdge(gp_Circ(gp::XOY(),outer),0,angle).Edge());
        wire.Add(BRepBuilderAPI_MakeEdge(point(outer,angle),point(inner,angle)).Edge());
        if(inner>tolerance) wire.Add(TopoDS::Edge(BRepBuilderAPI_MakeEdge(gp_Circ(gp::XOY(),inner),0,angle).Edge().Reversed()));
        wire.Add(BRepBuilderAPI_MakeEdge(point(inner,0),point(outer,0)).Edge());
        if(!wire.IsDone()) fail("unwrapFace",no_solution,"could not construct conical development");
        const auto flat=BRepBuilderAPI_MakeFace(gp_Pln(gp::XOY()),wire.Wire()).Face();
        GProp_GProps properties;BRepGProp::SurfaceProperties(flat,properties);
        auto centre=properties.CentreOfMass();gp_Trsf move;
        move.SetTranslation(gp_Vec(-centre.X(),-centre.Y(),0));
        return result(k,BRepBuilderAPI_Transform(flat,move,true).Shape());
    }
    if(surface.GetType()!=GeomAbs_Cylinder)
        fail("unwrapFace",unsupported_geometry,"exact unwrap supports planes and rectangular cylindrical or conical faces only");
    const double width=surface.Cylinder().Radius()*(surface.LastUParameter()-surface.FirstUParameter());
    const double height=surface.LastVParameter()-surface.FirstVParameter();
    return result(k,BRepBuilderAPI_MakeFace(gp_Pln(gp::XOY()),-width/2,width/2,-height/2,height/2).Face());
}
}
