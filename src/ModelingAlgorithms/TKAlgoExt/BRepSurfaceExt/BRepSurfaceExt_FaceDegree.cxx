// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepToolsExt_Operations.hxx>
#include <BRepToolsExt_History.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Geom2d_Curve.hxx>
#include <TopoDS_Edge.hxx>
#include <NCollection_Array2.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <TopExp.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopoDS.hxx>
#include <cmath>

namespace loft {
EvolutionData raiseFaceDegrees(OcctKernel& k,uint32_t id,const std::vector<uint32_t>& selected,
                               int du,int dv,double tolerance,const std::vector<int>& hashes,int bound) {
    const char* op="raiseFaceDegrees";
    if(du<0 || dv<0 || du>24 || dv>24 || !std::isfinite(tolerance) || tolerance<=0)
        fail(op,invalid_input,"degree increments must be 0..24 and tolerance positive");
    // A second copy lets the operation work correctly even when called directly from C++.
    const auto& input=get(k,id);BRepBuilderAPI_Copy copy(input,true,false);
    auto output=copy.Shape();
    NCollection_IndexedMap<TopoDS_Shape,TopTools_ShapeMapHasher> members;
    TopExp::MapShapes(input,TopAbs_FACE,members);
    FaceHistory history(op,hashes,bound);
    for(auto ref:selected) {
        const auto& original=get(k,ref);
        if(original.ShapeType()!=TopAbs_FACE || !members.Contains(original))
            fail(op,invalid_input,"select faces of the body to raise");
        auto face=TopoDS::Face(copy.ModifiedShape(original));
        TopLoc_Location location;auto support=BRep_Tool::Surface(face,location);
        auto basis=support;
        while(auto trim=Handle(Geom_RectangularTrimmedSurface)::DownCast(basis)) basis=trim->BasisSurface();
        auto spline=Handle(Geom_BSplineSurface)::DownCast(basis);
        if(spline.IsNull()) {
            if(BRepAdaptor_Surface(face).GetType()!=GeomAbs_Plane)
                fail(op,unsupported_geometry,"exact face elevation supports planes and B-spline surfaces");
            double a,b,c,d;BRepTools::UVBounds(face,a,b,c,d);
            NCollection_Array2<gp_Pnt> poles(1,2,1,2);
            poles(1,1)=support->Value(a,c);poles(2,1)=support->Value(b,c);
            poles(1,2)=support->Value(a,d);poles(2,2)=support->Value(b,d);
            NCollection_Array1<double> uk(1,2),vk(1,2);uk(1)=a;uk(2)=b;vk(1)=c;vk(2)=d;
            NCollection_Array1<int> mult(1,2);mult.Init(2);
            spline=new Geom_BSplineSurface(poles,uk,vk,mult,mult,1,1);
        } else spline=Handle(Geom_BSplineSurface)::DownCast(spline->Copy());
        if(spline->UDegree()+du>25 || spline->VDegree()+dv>25)
            fail(op,unsupported_geometry,"degree exceeds OpenCascade's maximum of 25");
        struct Curve { TopoDS_Edge edge; Handle(Geom2d_Curve) curve,second; double first,last; };
        std::vector<Curve> curves;
        auto forward=TopoDS::Face(face.Oriented(TopAbs_FORWARD));
        for(TopExp_Explorer it(forward,TopAbs_EDGE);it.More();it.Next()) {
            auto edge=TopoDS::Edge(it.Current().Oriented(TopAbs_FORWARD));double first,last;
            auto curve=BRep_Tool::CurveOnSurface(edge,forward,first,last);
            if(curve.IsNull()) fail(op,no_solution,"face edge has no parameter curve");
            Handle(Geom2d_Curve) second;
            if(BRep_Tool::IsClosed(edge,forward)) {
                double a,b;second=BRep_Tool::CurveOnSurface(TopoDS::Edge(edge.Reversed()),forward,a,b);
                if(second.IsNull()) fail(op,no_solution,"seam has no second parameter curve");
            }
            curves.push_back({edge,curve,second,first,last});
        }
        spline->IncreaseDegree(spline->UDegree()+du,spline->VDegree()+dv);
        // Exact elevation preserves the parameterization, wires, p-curves, seams and holes.
        BRep_Builder builder;builder.UpdateFace(face,spline,location,BRep_Tool::Tolerance(face));
        builder.NaturalRestriction(face,false);
        for(const auto& curve:curves) {
            // Reattach p-curves to the elevated support; preserve both sides of a seam.
            if(curve.second.IsNull()) builder.UpdateEdge(curve.edge,curve.curve,face,BRep_Tool::Tolerance(curve.edge));
            else builder.UpdateEdge(curve.edge,curve.curve,curve.second,face,BRep_Tool::Tolerance(curve.edge));
            builder.Range(curve.edge,face,curve.first,curve.last);
        }
    }
    if(selected.empty()) fail(op,invalid_input,"select at least one face");
    if(!BRepCheck_Analyzer(output).IsValid()) fail(op,no_solution,"degree elevation produced invalid topology");
    for(const auto& face:members) if(history.tracks(face)) {
        NCollection_List<TopoDS_Shape> successor;successor.Append(copy.ModifiedShape(face));
        history.modified(face,successor);
    }
    return history.finish(store(k,output));
}
}
