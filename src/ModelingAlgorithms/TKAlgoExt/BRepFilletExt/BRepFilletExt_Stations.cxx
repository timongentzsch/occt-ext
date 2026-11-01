// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_History.hxx>
#include <BRepAlgoExt_Cells.hxx>
#include <BRepFilletExt_Untouched.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepPrimAPI_MakeHalfSpace.hxx>
#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <NCollection_Array1.hxx>
#include <TopoDS.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Pln.hxx>
#include <algorithm>
#include <cmath>

namespace BRepExt {
// Stations are spatial points. OCCT's radius laws use contour abscissae, which we derive
// by projecting onto the named edge and measuring the curve, never from its tessellation.
EvolutionData filletStations(OcctKernel &kernel, uint32_t shapeId,
    const std::vector<uint32_t> &edgeIds, double radius,
    const std::vector<uint32_t> &stationEdges, const std::vector<double> &stations,
    const std::vector<uint32_t> &limitEdges, const std::vector<double> &limits,
    const std::vector<int> &hashes, int bound) {
    const char *op="filletStations";
    const auto &shape=get(kernel,shapeId);
    const auto edges=subShapes(kernel,op,shape,edgeIds,TopAbs_EDGE);
    auto references=[&](const std::vector<uint32_t> &ids) {
        std::vector<TopoDS_Shape> result;
        for(auto id:ids) result.push_back(subShapes(kernel,op,shape,{id},TopAbs_EDGE)[0]);
        return result;
    };
    const auto pointEdges=references(stationEdges);
    const auto stopEdges=references(limitEdges);
    if (!(radius>0) || !std::isfinite(radius) || stations.size()!=4*pointEdges.size() || limits.size()!=6*stopEdges.size())
        fail(op,invalid_input,"positive radius and complete spatial stations required");
    for(double value:stations) if(!std::isfinite(value)) fail(op,invalid_input,"nonfinite station");
    for(double value:limits) if(!std::isfinite(value)) fail(op,invalid_input,"nonfinite limit");
    const auto configure=[](BRepFilletAPI_MakeFillet& maker) {
        maker.SetParams(1e-2, BRepExt::tolerance::length, BRepExt::tolerance::uv,
                        BRepExt::tolerance::length, BRepExt::tolerance::uv, 1e-3);
    };
    BRepFilletAPI_MakeFillet all(shape);configure(all);
    for(const auto &edge:edges) if(!all.Contour(TopoDS::Edge(edge))) all.Add(radius,TopoDS::Edge(edge));
    for(const auto &edge:pointEdges) if(!all.Contour(TopoDS::Edge(edge))) fail(op,invalid_input,"station is outside the selected chains");
    for(const auto &edge:stopEdges) if(!all.Contour(TopoDS::Edge(edge))) fail(op,invalid_input,"limit is outside the selected chains");
    auto projected=[&](const TopoDS_Edge &edge,const gp_Pnt &p) {
        double first,last;
        const auto curve=BRep_Tool::Curve(edge,first,last);
        GeomAPI_ProjectPointOnCurve projection(p,curve,first,last);
        if(!projection.NbPoints() || projection.LowerDistance()>std::max(1e-5,1e-3*GCPnts_AbscissaPoint::Length(BRepAdaptor_Curve(edge))))
            fail(op,invalid_input,"station must lie on its named edge");
        return projection.LowerDistanceParameter();
    };
    auto setLaw=[&](BRepFilletAPI_MakeFillet &maker,int contour) {
        struct Piece { TopoDS_Edge edge; double from,length; };
        std::vector<Piece> pieces;
        double total=0;
        for(int j=1;j<=maker.NbEdges(contour);++j) {
            auto edge=maker.Edge(contour,j);
            BRepAdaptor_Curve curve(edge);
            const double length=GCPnts_AbscissaPoint::Length(curve);
            pieces.push_back({edge,total,length}); total+=length;
        }
        if(total<=BRepExt::tolerance::length) fail(op,invalid_input,"empty chain");
        std::vector<std::pair<double,double>> law;
        for(size_t i=0;i<pointEdges.size();++i) {
            const auto piece=std::find_if(pieces.begin(),pieces.end(),[&](const Piece &p){return p.edge.IsSame(pointEdges[i]);});
            if(piece==pieces.end()) continue;
            if(!(stations[4*i+3]>0)) fail(op,invalid_input,"station radius must be positive");
            BRepAdaptor_Curve curve(piece->edge);
            const auto t=projected(piece->edge,gp_Pnt(stations[4*i],stations[4*i+1],stations[4*i+2]));
            double length=GCPnts_AbscissaPoint::Length(curve,curve.FirstParameter(),t);
            if(piece->edge.Orientation()==TopAbs_REVERSED) length=piece->length-length;
            double at=piece->from+length;
            if(maker.Closed(contour) && total-at<BRepExt::tolerance::length) at=0;
            law.emplace_back(at,stations[4*i+3]);
        }
        if(law.empty()) return;
        const auto range=std::minmax_element(law.begin(),law.end(),[](const auto& a,const auto& b){return a.second<b.second;});
        if(range.second->second-range.first->second>BRepExt::tolerance::length &&
           std::any_of(pieces.begin(),pieces.end(),[](const Piece& p){
               const auto type=BRepAdaptor_Curve(p.edge).GetType();
               return type==GeomAbs_BSplineCurve || type==GeomAbs_BezierCurve;
           }))
            fail(op,unsupported_geometry,"variable radii on spline chains do not satisfy the rolling-ball envelope; construction is unavailable");

        std::sort(law.begin(),law.end());
        for(size_t i=1;i<law.size();++i) if(law[i].first-law[i-1].first<BRepExt::tolerance::length) {
            if(std::abs(law[i].second-law[i-1].second)>BRepExt::tolerance::length) fail(op,invalid_input,"conflicting radii at one station");
            law.erase(law.begin()+i); --i;
        }
        const bool closed=maker.Closed(contour);
        if(closed) {
            const auto first=law.front(),last=law.back();
            law.insert(law.begin(),{last.first-total,last.second});
            law.emplace_back(first.first+total,first.second);
        } else {
            if(law.front().first>BRepExt::tolerance::length) law.insert(law.begin(),{0,law.front().second});
            if(law.back().first<total-BRepExt::tolerance::length) law.emplace_back(total,law.back().second);
        }
        auto radiusAt=[&](double at) {
            for(size_t i=1;i<law.size();++i) if(at<=law[i].first+1e-8) {
                const auto a=law[i-1],b=law[i];
                return a.second+(b.second-a.second)*(at-a.first)/(b.first-a.first);
            }
            return law.back().second;
        };
        for(size_t i=0;i<pieces.size();++i) {
            const auto &p=pieces[i];
            std::vector<gp_Pnt2d> points{{0,radiusAt(p.from)}};
            for(const auto &s:law) if(s.first>p.from+BRepExt::tolerance::length && s.first<p.from+p.length-BRepExt::tolerance::length)
                points.emplace_back((s.first-p.from)/p.length,s.second);
            points.emplace_back(1,radiusAt(p.from+p.length));
            NCollection_Array1<gp_Pnt2d> values(1,int(points.size()));
            for(size_t j=0;j<points.size();++j) values(int(j)+1)=points[j];
            maker.SetRadius(values,contour,int(i)+1);
        }
    };
    FaceHistory history(op,hashes,bound);
    if(limits.empty()) {
        for(int c=1;c<=all.NbContours();++c) setLaw(all,c);
        all.Build();
        if(!all.IsDone() || !BRepCheck_Analyzer(all.Shape()).IsValid()) fail(op,no_solution,"station fillet failed");
        history.collect(all,shape);
        return history.finish(store(kernel,all.Shape()));
    }
    BRep_Builder builder;
    TopoDS_Compound removed,added;
    builder.MakeCompound(removed); builder.MakeCompound(added);
    for(int c=1;c<=all.NbContours();++c) {
        BRepFilletAPI_MakeFillet maker(shape);configure(maker);
        maker.Add(radius,all.Edge(c,1)); setLaw(maker,1); maker.Build();
        if(!maker.IsDone() || !BRepCheck_Analyzer(maker.Shape()).IsValid()) fail(op,no_solution,"limited station fillet failed");
        // A limit interval must occupy one connected run of the chain. Halfspaces on a
        // winding spline can otherwise round several disjoint runs while returning a valid
        // body. Check the contour in its directed order, including the closed-chain seam.
        struct Stop { gp_Pnt point; gp_Vec direction; };
        std::vector<Stop> stops;
        for(size_t i=0;i<stopEdges.size();++i) {
            if(all.Contour(TopoDS::Edge(stopEdges[i]))!=c) continue;
            const auto edge=TopoDS::Edge(stopEdges[i]);
            const double at=projected(edge,gp_Pnt(limits[6*i],limits[6*i+1],limits[6*i+2]));
            gp_Vec direction(limits[6*i+3],limits[6*i+4],limits[6*i+5]);
            if(direction.Magnitude()<1e-12) fail(op,invalid_input,"limit direction has no length");
            stops.push_back({BRepAdaptor_Curve(edge).Value(at),direction.Normalized()});
        }
        if(stops.size()>2) fail(op,invalid_input,"at most two limits per chain are supported");
        std::vector<bool> kept;
        for(int j=1;j<=maker.NbEdges(1);++j) {
            const auto edge=maker.Edge(1,j);BRepAdaptor_Curve curve(edge);
            for(int sample=0;sample<257;++sample) {
                const double f=(sample+.5)/257;
                const double t=edge.Orientation()==TopAbs_REVERSED ?
                    curve.LastParameter()-(curve.LastParameter()-curve.FirstParameter())*f :
                    curve.FirstParameter()+(curve.LastParameter()-curve.FirstParameter())*f;
                const auto p=curve.Value(t);bool inside=true;
                for(const auto& stop:stops) inside &= gp_Vec(stop.point,p).Dot(stop.direction)>-BRepExt::tolerance::length;
                kept.push_back(inside);
            }
        }
        int transitions=0;
        for(size_t i=1;i<kept.size();++i) transitions += kept[i]!=kept[i-1];
        if(!kept.empty() && maker.Closed(1) && kept.front()!=kept.back()) ++transitions;
        if(transitions>2) fail(op,unsupported_geometry,"limited spline chain re-enters its stopping planes");
        if(!stops.empty() && std::none_of(kept.begin(),kept.end(),[](bool value){return value;}))
            fail(op,invalid_input,"limits leave no interval on the selected chain");
        auto clip=[&](TopoDS_Shape tool) {
            for(size_t i=0;i<stopEdges.size();++i) {
                if(all.Contour(TopoDS::Edge(stopEdges[i]))!=c) continue;
                const auto edge=TopoDS::Edge(stopEdges[i]);
                const double at=projected(edge,gp_Pnt(limits[6*i],limits[6*i+1],limits[6*i+2]));
                const gp_Pnt p=BRepAdaptor_Curve(edge).Value(at);
                const gp_Vec direction(limits[6*i+3],limits[6*i+4],limits[6*i+5]);
                if(direction.Magnitude()<1e-12) fail(op,invalid_input,"limit direction has no length");
                const auto plane=BRepBuilderAPI_MakeFace(gp_Pln(p,gp_Dir(direction))).Face();
                const auto half=BRepPrimAPI_MakeHalfSpace(plane,p.Translated(direction.Normalized())).Solid();
                Untouched<BRepAlgoAPI_Common> common(tool,half);
                if(!common.IsDone()) fail(op,no_solution,"limit intersection failed");
                tool=common.Shape();
            }
            return tool;
        };
        // Closely coincident spline supports can make an ordinary difference lose most
        // of the blend volume without reporting an error. Use the shared verified partition.
        const auto take=booleanCells({op,Select::Difference,{shape},{maker.Shape()}});
        const auto give=booleanCells({op,Select::Difference,{maker.Shape()},{shape}});
        builder.Add(removed,clip(take.shape)); builder.Add(added,clip(give.shape));
    }
    Untouched<BRepAlgoAPI_Cut> cut(shape,removed);
    if(!cut.IsDone()) fail(op,no_solution,"limited blend subtraction failed");
    Untouched<BRepAlgoAPI_Fuse> fuse(cut.Shape(),added);
    if(!fuse.IsDone() || !BRepCheck_Analyzer(fuse.Shape()).IsValid()) fail(op,no_solution,"limited blend is invalid");
    history.collectComposed(cut,fuse,shape);
    return history.finish(store(kernel,fuse.Shape()));
}
} // namespace BRepExt
