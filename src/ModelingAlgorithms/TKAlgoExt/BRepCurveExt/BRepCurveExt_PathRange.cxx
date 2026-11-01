// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Tool.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Wire.hxx>
#include <cmath>
namespace loft {
EvolutionData trimPath(OcctKernel &kernel,uint32_t id,double start,double end) {
    const char *op="trimPath";
    const auto shape=get(kernel,id);
    if(shape.ShapeType()!=TopAbs_WIRE || !std::isfinite(start) || !std::isfinite(end) || start<0 || end>1 || start>=end)
        fail(op,invalid_input,"wire and increasing distances within [0,1] required");
    double total=0;
    for(BRepTools_WireExplorer it(TopoDS::Wire(shape));it.More();it.Next()) {
        BRepAdaptor_Curve c(it.Current()); total+=GCPnts_AbscissaPoint::Length(c);
    }
    if(total<=loft::tolerance::exactLength) fail(op,invalid_input,"path has no length");
    BRepBuilderAPI_MakeWire wire;
    double distance=0;
    for(BRepTools_WireExplorer it(TopoDS::Wire(shape));it.More();it.Next()) {
        const auto edge=it.Current(); BRepAdaptor_Curve c(edge);
        const double len=GCPnts_AbscissaPoint::Length(c);
        const double a=std::max(0.0,start*total-distance),b=std::min(len,end*total-distance);
        distance+=len;
        if(b-a<=loft::tolerance::exactLength) continue;
        if(a<=loft::tolerance::exactLength && len-b<=loft::tolerance::exactLength) {wire.Add(edge);continue;}
        const bool reverse=edge.Orientation()==TopAbs_REVERSED;
        const double begin=reverse ? c.LastParameter():c.FirstParameter();
        auto parameter=[&](double length) {
            if(length<=loft::tolerance::exactLength) return begin;
            if(len-length<=loft::tolerance::exactLength) return reverse ? c.FirstParameter():c.LastParameter();
            GCPnts_AbscissaPoint solve(1e-9,c,reverse ? -length:length,begin);
            if(!solve.IsDone()) fail(op,no_solution,"arc length inversion failed");
            return solve.Parameter();
        };
        double u=parameter(a),v=parameter(b),first,last;
        auto geometry=BRep_Tool::Curve(edge,first,last);
        if(geometry.IsNull()) fail(op,unsupported_geometry,"path edge has no 3D curve");
        auto trimmed=BRepBuilderAPI_MakeEdge(geometry,std::min(u,v),std::max(u,v)).Edge();
        if(reverse) trimmed.Reverse();
        wire.Add(trimmed);
    }
    if(!wire.IsDone()) fail(op,no_solution,"trimmed path is disconnected");
    return stored(kernel, wire.Wire());
}
}
