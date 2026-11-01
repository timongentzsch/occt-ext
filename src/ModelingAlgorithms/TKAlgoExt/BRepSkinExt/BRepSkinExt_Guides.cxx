// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_CubicThrough.hxx>
#include <BRepSurfaceExt_Tools.hxx>
#include <BRepAdaptor_CompCurve.hxx>
#include <BOPAlgo_ArgumentAnalyzer.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepOffsetAPI_MakeFilling.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <GProp_GProps.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Wire.hxx>
#include <cmath>
#include <gp_Pln.hxx>
#include <gp_Ax3.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRep_Tool.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <GeomConvert.hxx>
#include <array>
#include <functional>
#include <tuple>

namespace BRepExt {
namespace {
/// An edge's curve as a spline from 0 to 1, run the other way where `reversed`.
Handle(Geom_BSplineCurve) unitCurve(const TopoDS_Edge &edge, bool reversed) {
    double a, b;
    const auto curve = BRep_Tool::Curve(edge, a, b);
    Handle(Geom_BSplineCurve) spline = GeomConvert::CurveToBSplineCurve(new Geom_TrimmedCurve(curve, a, b));
    if ((edge.Orientation() == TopAbs_REVERSED) != reversed)
        spline->Reverse();
    NCollection_Array1<double> knots(1, spline->NbKnots());
    spline->Knots(knots);
    BSplCLib::Reparametrize(0, 1, knots);
    spline->SetKnots(knots);
    return spline;
}

/// A guided skin between two profiles with G2 asked at an end: the sum of the surface straight
/// across between the guides, and the quintic in the guides' direction that has each profile, the
/// support's plane and bend along it, less the same between the guides' ends. Along a guide the
/// last two cancel, so each guide is kept; along a profile the first and last cancel. The two supports
/// give the plane and the bend, the guides how fast the surface leaves. Between flat supports
/// this is the guide swept straight across.
TopoDS_Face guidedSum(OcctKernel &kernel, const char *op, const std::vector<uint32_t> &edges, uint32_t firstFace,
                      uint32_t lastFace, int firstOrder, int lastOrder, double tolerance) {
    const auto edgeAt = [&](size_t i) { return TopoDS::Edge(get(kernel, edges[i])); };
    const auto startOf = [](const TopoDS_Edge &e) {
        double a, b;
        const auto c = BRep_Tool::Curve(e, a, b);
        return c->Value(e.Orientation() == TopAbs_REVERSED ? b : a);
    };
    const auto endOf = [](const TopoDS_Edge &e) {
        double a, b;
        const auto c = BRep_Tool::Curve(e, a, b);
        return c->Value(e.Orientation() == TopAbs_REVERSED ? a : b);
    };
    const auto first = edgeAt(0), last = edgeAt(1), g0 = edgeAt(2), g1 = edgeAt(3);
    // Each guide from the first profile to the last, and each profile from the first guide to the second.
    const auto distanceTo = [](const gp_Pnt &p, const TopoDS_Edge &e) {
        BRepExtrema_DistShapeShape near(BRepBuilderAPI_MakeVertex(p).Vertex(), e);
        return near.Value();
    };
    const bool g0Back = distanceTo(startOf(g0), first) > distanceTo(endOf(g0), first);
    const bool g1Back = distanceTo(startOf(g1), first) > distanceTo(endOf(g1), first);
    const auto guideA = unitCurve(g0, g0Back), guideB = unitCurve(g1, g1Back);
    const bool firstBack = startOf(first).Distance(guideA->StartPoint()) > endOf(first).Distance(guideA->StartPoint());
    const bool lastBack = startOf(last).Distance(guideA->EndPoint()) > endOf(last).Distance(guideA->EndPoint());
    const auto rimA = unitCurve(first, firstBack), rimB = unitCurve(last, lastBack);
    // A profile at a share of its edge's parameter, which is where a support's jet is read; the
    // spline of a conic runs at another pace.
    const auto pointOn = [](const TopoDS_Edge &edge, bool back) {
        double a, b;
        const Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, a, b);
        const bool flipped = (edge.Orientation() == TopAbs_REVERSED) != back;
        return [curve, a, b, flipped](double u) { return curve->Value(a + (b - a) * (flipped ? 1 - u : u)); };
    };
    const auto rimPointA = pointOn(first, firstBack), rimPointB = pointOn(last, lastBack);
    const auto jetsOf = [&](uint32_t edgeId, uint32_t faceId, int order, bool back) -> std::function<SurfaceJet(double)> {
        if (order < 1)
            return {};
        const auto boundary = surfacing::boundaryOf(kernel, op, edgeId, faceId);
        const bool flipped = (boundary.edge.Orientation() == TopAbs_REVERSED) != back;
        return [boundary, flipped, op](double u) { return boundary.jet(op, flipped ? 1 - u : u); };
    };
    // How the surface leaves a profile at u: the guides' way there, in the support's plane, with
    // the support's bend in that direction.
    const auto leaving = [&](const std::function<SurfaceJet(double)> &jetAt, int order, double v, double u) {
        gp_Pnt ignored;
        gp_Vec a1, a2, b1, b2;
        guideA->D2(v, ignored, a1, a2);
        guideB->D2(v, ignored, b1, b2);
        gp_Vec d1 = a1 * (1 - u) + b1 * u, d2 = a2 * (1 - u) + b2 * u;
        if (order >= 1) {
            const auto jet = jetAt(u);
            d1 -= jet.n * d1.Dot(jet.n);
            if (order >= 2)
                d2 += jet.n * (jet.curvature(d1, d1) - d2.Dot(jet.n));
        }
        return std::pair{d1, d2};
    };
    const auto jetsA = jetsOf(edges[0], firstFace, firstOrder, firstBack), jetsB = jetsOf(edges[1], lastFace, lastOrder, lastBack);
    const auto splineOf = [](const std::vector<gp_Pnt> &values, const std::vector<double> &places) {
        const GeomAlgoExt::CubicThrough law(places, GeomAlgoExt::CubicThrough::Ends::Follow);
        const auto poles = law.through(values);
        NCollection_Array1<gp_Pnt> array(1, law.poles);
        for (int k = 1; k <= law.poles; ++k)
            array(k) = poles[size_t(k - 1)];
        return Handle(Geom_BSplineCurve)(new Geom_BSplineCurve(array, law.knots, law.mults, 3));
    };
    // The fields along each profile as splines, and a rational profile the same way, for the sum
    // of poles, to the tolerance: over a span of 1 of the guides' parameter these are lengths.
    Handle(Geom_BSplineCurve) a1, a2, b1, b2, onA = rimA, onB = rimB;
    for (int n = 16;; n *= 2) {
        if (n > 4096)
            fail(op, no_solution, "the supports vary too much along the profiles to follow");
        std::vector<double> places;
        std::vector<gp_Pnt> va1, va2, vb1, vb2, pa, pb;
        for (int k = 0; k <= n; ++k) {
            const double u = double(k) / n;
            places.push_back(u);
            const auto [d1, d2] = leaving(jetsA, firstOrder, 0, u);
            const auto [e1, e2] = leaving(jetsB, lastOrder, 1, u);
            va1.emplace_back(d1.XYZ());
            va2.emplace_back(d2.XYZ());
            vb1.emplace_back(e1.XYZ());
            vb2.emplace_back(e2.XYZ());
            pa.push_back(rimPointA(u));
            pb.push_back(rimPointB(u));
        }
        a1 = splineOf(va1, places);
        a2 = splineOf(va2, places);
        b1 = splineOf(vb1, places);
        b2 = splineOf(vb2, places);
        if (rimA->IsRational())
            onA = splineOf(pa, places);
        if (rimB->IsRational())
            onB = splineOf(pb, places);
        double worst = 0;
        for (int k = 0; k < n; ++k) {
            const double u = (k + 0.5) / n;
            const auto [d1, d2] = leaving(jetsA, firstOrder, 0, u);
            const auto [e1, e2] = leaving(jetsB, lastOrder, 1, u);
            worst = std::max({worst, a1->Value(u).Distance(gp_Pnt(d1.XYZ())), a2->Value(u).Distance(gp_Pnt(d2.XYZ())),
                              b1->Value(u).Distance(gp_Pnt(e1.XYZ())), b2->Value(u).Distance(gp_Pnt(e2.XYZ())),
                              onA->Value(u).Distance(rimPointA(u)), onB->Value(u).Distance(rimPointB(u))});
        }
        if (worst <= tolerance / 4)
            break;
    }
    // The quintic: six rows of poles along the profiles.
    std::vector<Handle(Geom_BSplineCurve)> along{Handle(Geom_BSplineCurve)::DownCast(onA->Copy()), a1, a2, b1, b2,
                                                 Handle(Geom_BSplineCurve)::DownCast(onB->Copy())};
    GeomAlgoExt::unify(along);
    const int nu = along[0]->NbPoles();
    NCollection_Array2<gp_Pnt> quintic(1, nu, 1, 6);
    for (int a = 1; a <= nu; ++a) {
        const gp_XYZ p0 = along[0]->Pole(a).XYZ(), d1 = along[1]->Pole(a).XYZ(), d2 = along[2]->Pole(a).XYZ(),
                     e1 = along[3]->Pole(a).XYZ(), e2 = along[4]->Pole(a).XYZ(), p1 = along[5]->Pole(a).XYZ();
        quintic(a, 1) = gp_Pnt(p0);
        quintic(a, 2) = gp_Pnt(p0 + d1 / 5);
        quintic(a, 3) = gp_Pnt(p0 + d1 * 2 / 5 + d2 / 20);
        quintic(a, 4) = gp_Pnt(p1 - e1 * 2 / 5 + e2 / 20);
        quintic(a, 5) = gp_Pnt(p1 - e1 / 5);
        quintic(a, 6) = gp_Pnt(p1);
    }
    NCollection_Array1<double> bezierKnots(1, 2);
    NCollection_Array1<int> quinticMults(1, 2), lineMults(1, 2);
    bezierKnots(1) = 0;
    bezierKnots(2) = 1;
    quinticMults(1) = quinticMults(2) = 6;
    lineMults(1) = lineMults(2) = 2;
    Handle(Geom_BSplineSurface) profilesSkin = new Geom_BSplineSurface(
        quintic, along[0]->Knots(), bezierKnots, along[0]->Multiplicities(), quinticMults, along[0]->Degree(), 5);
    // The same quintic between the ends of the profiles, straight across.
    NCollection_Array2<gp_Pnt> corners(1, 2, 1, 6);
    for (int b = 1; b <= 6; ++b) {
        corners(1, b) = quintic(1, b);
        corners(2, b) = quintic(nu, b);
    }
    Handle(Geom_BSplineSurface) crossings = new Geom_BSplineSurface(corners, bezierKnots, bezierKnots, lineMults, quinticMults, 1, 5);
    // The surface straight across between the guides.
    std::vector<Handle(Geom_BSplineCurve)> sides{guideA, guideB};
    GeomAlgoExt::unify(sides);
    const int nv = sides[0]->NbPoles();
    NCollection_Array2<gp_Pnt> straight(1, 2, 1, nv);
    for (int b = 1; b <= nv; ++b) {
        straight(1, b) = sides[0]->Pole(b);
        straight(2, b) = sides[1]->Pole(b);
    }
    Handle(Geom_BSplineSurface) guidesSkin =
        new Geom_BSplineSurface(straight, bezierKnots, sides[0]->Knots(), lineMults, sides[0]->Multiplicities(), 1, sides[0]->Degree());
    // On one set of knots, and summed.
    const std::vector<Handle(Geom_BSplineSurface)> skins{profilesSkin, crossings, guidesSkin};
    int du = 1, dv = 1;
    for (const auto &skin : skins) {
        du = std::max(du, skin->UDegree());
        dv = std::max(dv, skin->VDegree());
    }
    for (const auto &skin : skins)
        skin->IncreaseDegree(du, dv);
    for (const auto &from : skins)
        for (const auto &to : skins) {
            for (int k = 2; k < from->NbUKnots(); ++k)
                to->InsertUKnot(from->UKnot(k), from->UMultiplicity(k), 1e-12, false);
            for (int k = 2; k < from->NbVKnots(); ++k)
                to->InsertVKnot(from->VKnot(k), from->VMultiplicity(k), 1e-12, false);
        }
    for (int a = 1; a <= profilesSkin->NbUPoles(); ++a)
        for (int b = 1; b <= profilesSkin->NbVPoles(); ++b)
            profilesSkin->SetPole(a, b, gp_Pnt(profilesSkin->Pole(a, b).XYZ() + guidesSkin->Pole(a, b).XYZ() - crossings->Pole(a, b).XYZ()));
    const TopoDS_Face face = surfacing::makeFace(op, profilesSkin, tolerance);
    // What comes back is held to the supports apart from how it was made.
    for (const auto &[edgeId, faceId, order, v] : std::array{std::tuple{edges[0], firstFace, firstOrder, 0.0}, std::tuple{edges[1], lastFace, lastOrder, 1.0}}) {
        if (order < 1)
            continue;
        const auto reference = surfacing::boundaryOf(kernel, op, edgeId, faceId);
        TopoDS_Edge own;
        double nearest = Precision::Infinite();
        for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next()) {
            BRepAdaptor_Curve curve(TopoDS::Edge(it.Current()));
            const double d = curve.Value((curve.FirstParameter() + curve.LastParameter()) / 2).Distance(reference.point(0.5));
            if (d < nearest) {
                nearest = d;
                own = TopoDS::Edge(it.Current());
            }
        }
        (void)v;
        surfacing::measure(op, face, {surfacing::constraint(op, face, own, reference, order)}, tolerance, false);
    }
    return face;
}
} // namespace

EvolutionData skinGuides(OcctKernel &kernel, const std::vector<uint32_t> &edges, uint32_t firstFace,
                         uint32_t lastFace, int firstOrder, int lastOrder, double tolerance) {
    const char *op = "skinGuides";
    if (edges.size() < 4 || firstOrder < 0 || firstOrder > 2 || lastOrder < 0 || lastOrder > 2 ||
        !std::isfinite(tolerance) || tolerance <= 0)
        fail(op, invalid_input, "at least two profiles, two guides and G0-G2 end constraints are required");
    BRepOffsetAPI_MakeFilling fill(3, 50, 5, false, tolerance / 10, tolerance / 10, 1e-5, 1e-4, 8, 128);
    BRepBuilderAPI_MakeWire boundary;
    // Connectivity is checked independently from filling, which otherwise accepts
    // disconnected constraints.
    for (auto id : edges) {
        if (get(kernel, id).ShapeType() != TopAbs_EDGE)
            fail(op, invalid_input, "each profile and guide must be a single edge");
    }
    const int last=int(edges.size())-3;
    for (int i : {0, last+1, last, last+2})
        boundary.Add(TopoDS::Edge(get(kernel, edges[i])));
    if (!boundary.IsDone() || !boundary.Wire().Closed())
        fail(op, invalid_input, "profiles and guides must form one closed four-sided boundary");
    if (edges.size() == 4 && (firstOrder == 2 || lastOrder == 2)) {
        if ((firstOrder > 0 && !firstFace) || (lastOrder > 0 && !lastFace))
            fail(op, invalid_input, "G1/G2 requires an adjacent support face");
        return stored(kernel, guidedSum(kernel, op, edges, firstFace, lastFace, firstOrder, lastOrder, tolerance));
    }
    for (int i = 0; i < int(edges.size()); ++i) {
        const auto edge = TopoDS::Edge(get(kernel, edges[i]));
        const int order = i == 0 ? firstOrder : i == last ? lastOrder : 0;
        const uint32_t support = i == 0 ? firstFace : i == last ? lastFace : 0;
        if (order) {
            if (!support || get(kernel, support).ShapeType() != TopAbs_FACE)
                fail(op, invalid_input, "G1/G2 requires an adjacent support face");
            subShapes(kernel, op, get(kernel, support), {edges[i]}, TopAbs_EDGE);
            // This OCCT build forwards GeomAbs_Shape numerically to the plate
            // constraint, which accepts integer orders 0, 1, 2 (not enum G2=3).
            fill.Add(edge, TopoDS::Face(get(kernel, support)), static_cast<GeomAbs_Shape>(order));
        } else
            fill.Add(edge, GeomAbs_C0, i==0 || i>=last);
    }
    fill.Build();
    if (!fill.IsDone() || fill.G0Error() > tolerance || !BRepCheck_Analyzer(fill.Shape()).IsValid())
        fail(op, no_solution, "guided skin did not meet its boundary tolerance");
    if ((firstOrder > 0 || lastOrder > 0) && fill.G1Error() > 1e-4)
        fail(op, no_solution, "guided skin did not meet G1");
    if ((firstOrder > 1 || lastOrder > 1) && fill.G2Error() > 1e-3)
        fail(op, no_solution, "guided skin did not meet G2: " + std::to_string(fill.G2Error()));
    return stored(kernel, fill.Shape());
}
} // namespace BRepExt

#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <NCollection_Array2.hxx>
namespace BRepExt {
EvolutionData skinPeriodic(OcctKernel &kernel,const std::vector<uint32_t> &ids,double tolerance) {
    const char *op="skinPeriodic";
    if(ids.size()<3 || !std::isfinite(tolerance) || tolerance<=0)
        fail(op,invalid_input,"at least three sections and positive tolerance required");
    constexpr int samples=65;
    NCollection_Array2<gp_Pnt> grid(1,int(ids.size()),1,samples);
    std::vector<TopoDS_Edge> edges;
    for(size_t i=0;i<ids.size();++i) {
        const auto &shape=get(kernel,ids[i]);
        TopExp_Explorer it(shape,TopAbs_EDGE);
        if(!it.More()) fail(op,invalid_input,"empty profile");
        auto edge=TopoDS::Edge(it.Current());it.Next();
        if(it.More() || edge.Closed()) fail(op,unsupported_geometry,"periodic skin requires single open-edge sections");
        BRepAdaptor_Curve curve(edge);
        double a=curve.FirstParameter(),b=curve.LastParameter();
        if(edge.Orientation()==TopAbs_REVERSED) std::swap(a,b);
        if(i>0) {
            const auto previousStart=grid(int(i),1),previousEnd=grid(int(i),samples);
            const double same=previousStart.SquareDistance(curve.Value(a))+previousEnd.SquareDistance(curve.Value(b));
            const double reverse=previousStart.SquareDistance(curve.Value(b))+previousEnd.SquareDistance(curve.Value(a));
            if(reverse<same) {std::swap(a,b);edge.Reverse();}
        }
        edges.push_back(edge);
        for(int j=0;j<samples;++j) grid(int(i)+1,j+1)=curve.Value(a+(b-a)*j/(samples-1));
    }
    GeomAPI_PointsToBSplineSurface fit;
    fit.Interpolate(grid,Approx_IsoParametric,true);
    if(!fit.IsDone()) fail(op,no_solution,"periodic interpolation failed");
    auto surface=fit.Surface();
    if(!surface->IsUPeriodic()) fail(op,no_solution,"interpolation did not close tangentially");
    double u0,u1,v0,v1;surface->Bounds(u0,u1,v0,v1);
    // Check between interpolation sites as well. This is a sampled tolerance
    // check, not a certified error bound for arbitrarily oscillatory profiles.
    for(size_t i=0;i<edges.size();++i) {
        BRepAdaptor_Curve curve(edges[i]);double a=curve.FirstParameter(),b=curve.LastParameter();
        if(edges[i].Orientation()==TopAbs_REVERSED) std::swap(a,b);
        for(int j=0;j<2*samples-1;++j) {
            const double f=double(j)/(2*samples-2);
            const auto expected=curve.Value(a+(b-a)*f);
            if(expected.Distance(surface->Value(u0+(u1-u0)*i/ids.size(),v0+(v1-v0)*f))>tolerance)
                fail(op,no_solution,"periodic skin misses a section between samples; rebuild or simplify the section");
        }
    }
    auto face=BRepBuilderAPI_MakeFace(surface,tolerance).Face();
    if(!BRepCheck_Analyzer(face).IsValid()) fail(op,no_solution,"invalid periodic surface");
    return stored(kernel, face);
}
}
