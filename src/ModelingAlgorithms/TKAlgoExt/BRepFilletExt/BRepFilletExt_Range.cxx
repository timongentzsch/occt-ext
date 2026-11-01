// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_History.hxx>
#include <BRepFilletExt_Untouched.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeHalfSpace.hxx>
#include <GProp_GProps.hxx>
#include <TopoDS.hxx>
#include <cmath>
#include <NCollection_Array1.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Pln.hxx>

namespace BRepExt {
EvolutionData filletRange(OcctKernel &kernel, uint32_t shapeId, uint32_t edgeId, double radius, double start,
                          double end, bool propagate, const std::vector<int> &hashes, int bound, const std::vector<double> &law) {
    const char *op = "filletRange";
    const auto &shape = get(kernel, shapeId);
    const auto selected = subShapes(kernel, op, shape, {edgeId}, TopAbs_EDGE);
    const auto edge = TopoDS::Edge(selected[0]);
    if (!std::isfinite(radius) || radius <= 0 || !std::isfinite(start) || !std::isfinite(end) || start < 0 ||
        end > 1 || start >= end)
        fail(op, invalid_input, "positive radius and ordered range within [0,1] required");
    BRepAdaptor_Curve rangeCurve(edge);
    const double rangeLength=GCPnts_AbscissaPoint::Length(rangeCurve,
        rangeCurve.FirstParameter()+(rangeCurve.LastParameter()-rangeCurve.FirstParameter())*start,
        rangeCurve.FirstParameter()+(rangeCurve.LastParameter()-rangeCurve.FirstParameter())*end);
    if(rangeLength<=BRepExt::tolerance::length) fail(op,invalid_input,"empty range: limits are within kernel resolution");
    BRepFilletAPI_MakeFillet maker(shape);
    maker.Add(radius, edge);
    const int contour = maker.Contour(edge);
    if(!law.empty()) {
        if(law.size()<4 || law.size()%2 || law.front()!=0 || law[law.size()-2]!=1)
            fail(op,invalid_input,"radius points must run from 0 to 1");
        const int size=int(law.size()/2);
        int index=1;while(index<=maker.NbEdges(contour) && !maker.Edge(contour,index).IsSame(edge)) ++index;
        if(index>maker.NbEdges(contour)) fail(op,no_solution,"selected edge is absent from contour");
        const bool reverse=maker.Edge(contour,index).Orientation()!=edge.Orientation();
        NCollection_Array1<gp_Pnt2d> points(1,size);
        for(int i=0;i<size;++i) {
            if(!std::isfinite(law[2*i]) || !std::isfinite(law[2*i+1]) || law[2*i+1]<=0 ||
               (i && law[2*i]<=law[2*i-2])) fail(op,invalid_input,"ordered finite parameters and positive radii required");
            points(reverse ? size-i:i+1)=gp_Pnt2d(reverse ? 1-law[2*i]:law[2*i],law[2*i+1]);
        }
        maker.SetRadius(points,contour,index);
    }
    BRepAdaptor_Curve curve(edge);
    const bool limited = start > 0 || end < 1 || (!propagate && maker.NbEdges(contour) > 1);
    maker.Build();
    if (!maker.IsDone() || !BRepCheck_Analyzer(maker.Shape()).IsValid())
        fail(op, no_solution, "fillet construction failed");
    FaceHistory history(op, hashes, bound);
    if (!limited) {
        history.collect(maker, shape);
        return history.finish(store(kernel, maker.Shape()));
    }
    // Isolate both sides of the change. Concave blends add material; convex
    // blends remove it. Mixed contours can do both in one operation.
    Untouched<BRepAlgoAPI_Cut> removed(shape, maker.Shape()), added(maker.Shape(), shape);
    if (!removed.IsDone() || !added.IsDone())
        fail(op, no_solution, "could not isolate fillet volume");
    double a = curve.FirstParameter(), b = curve.LastParameter();
    if (edge.Orientation() == TopAbs_REVERSED) std::swap(a, b);
    gp_Pnt p0, p1;
    gp_Vec t0, t1;
    curve.D1(a + (b-a)*start, p0, t0);
    curve.D1(a + (b-a)*end, p1, t1);
    if (t0.Magnitude() < 1e-12 || t1.Magnitude() < 1e-12)
        fail(op, unsupported_geometry, "a limit has a singular tangent");
    if (b < a) { t0.Reverse(); t1.Reverse(); }
    const gp_Dir n0(t0), n1(t1);
    auto inside = [&](const gp_Pnt &p, double tol) {
        return gp_Vec(p0,p).Dot(gp_Vec(n0)) > tol && gp_Vec(p1,p).Dot(gp_Vec(n1)) < -tol;
    };
    // A pair of normal stopping planes cannot describe a winding interval.
    // Refuse those inputs instead of accidentally rounding a second arc.
    for (int i=0; i<=256; ++i) {
        const double f=double(i)/256;
        const auto p=curve.Value(a+(b-a)*f);
        if (f > start+1e-8 && f < end-1e-8) {
            if (!inside(p,-BRepExt::tolerance::length)) fail(op, unsupported_geometry, "range crosses its normal stopping planes");
        } else if ((f < start-1e-8 || f > end+1e-8) && inside(p,BRepExt::tolerance::length))
            fail(op, unsupported_geometry, "range returns inside its normal stopping planes");
    }
    for (int e=1; e<=maker.NbEdges(contour); ++e) {
        if (maker.Edge(contour,e).IsSame(edge)) continue;
        BRepAdaptor_Curve other(maker.Edge(contour,e));
        for (int i=1; i<64; ++i)
            if (inside(other.Value(other.FirstParameter()+(other.LastParameter()-other.FirstParameter())*i/64),BRepExt::tolerance::length))
                fail(op, unsupported_geometry, "tangent continuation re-enters the limited interval");
    }
    auto clip = [&](TopoDS_Shape tool) {
        for (int i=0; i<2; ++i) {
            const auto point=i ? p1:p0;
            const auto normal=i ? n1:n0;
            const auto plane=BRepBuilderAPI_MakeFace(gp_Pln(point,normal)).Face();
            const auto witness=point.Translated(gp_Vec(normal)*(i ? -1:1));
            const auto half=BRepPrimAPI_MakeHalfSpace(plane,witness).Solid();
            BRepAlgoAPI_Common common(tool,half);
            if (!common.IsDone()) fail(op,no_solution,"fillet limit intersection failed");
            tool=common.Shape();
        }
        return tool;
    };
    Untouched<BRepAlgoAPI_Cut> cut(shape,clip(removed.Shape()));
    if (!cut.IsDone()) fail(op,no_solution,"limited fillet subtraction failed");
    BRepAlgoAPI_Fuse fuse(cut.Shape(),clip(added.Shape()));
    if (!fuse.IsDone() || !BRepCheck_Analyzer(fuse.Shape()).IsValid())
        fail(op,no_solution,"limited fillet is invalid");
    history.collectComposed(cut, fuse, shape);
    return history.finish(store(kernel,fuse.Shape()));
}
} // namespace BRepExt
