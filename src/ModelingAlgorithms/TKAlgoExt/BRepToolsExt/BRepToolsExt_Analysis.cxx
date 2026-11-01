// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_SurfaceJet.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRep_Tool.hxx>
#include <Geom2d_Curve.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <algorithm>
#include <cmath>

namespace BRepExt {
std::vector<double> minimumDistance(OcctKernel &kernel, uint32_t first, uint32_t second) {
    BRepExtrema_DistShapeShape distance(get(kernel, first), get(kernel, second));
    if (!distance.IsDone() || distance.NbSolution() < 1)
        fail("minimumDistance", no_solution, "no closest points found");
    const auto a = distance.PointOnShape1(1), b = distance.PointOnShape2(1);
    return {distance.Value(), a.X(), a.Y(), a.Z(), b.X(), b.Y(), b.Z()};
}
std::vector<double> faceContinuity(OcctKernel &kernel, uint32_t first, uint32_t second, uint32_t edgeId,
                                   int samples, uint32_t otherEdge) {
    const char *op = "faceContinuity";
    if (get(kernel, first).ShapeType() != TopAbs_FACE || get(kernel, second).ShapeType() != TopAbs_FACE ||
        get(kernel, edgeId).ShapeType() != TopAbs_EDGE || samples < 2 || samples > 4096)
        fail(op, invalid_input, "two faces, a shared edge and 2-4096 samples are required");
    const auto a = TopoDS::Face(get(kernel, first)), b = TopoDS::Face(get(kernel, second));
    const auto edge = TopoDS::Edge(get(kernel, edgeId));
    if (a.IsSame(b))
        fail(op, invalid_input, "two distinct faces are required");
    const auto edgeB = otherEdge ? TopoDS::Edge(get(kernel, otherEdge)) : edge;
    const std::pair<TopoDS_Face, TopoDS_Edge> sides[] = {{a, edge}, {b, edgeB}};
    for (const auto &[face, ownedEdge] : sides) {
        bool found = false;
        for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next())
            if (it.Current().IsSame(ownedEdge))
                found = true;
        if (!found)
            fail(op, invalid_input, "the edge must belong to both faces");
    }
    BRepAdaptor_Curve curve(edge), curveB(edgeB);
    const double firstA = curve.FirstParameter(), lastA = curve.LastParameter();
    const double firstB = curveB.FirstParameter(), lastB = curveB.LastParameter();
    // Quarter samples disambiguate closed edges, whose start and end coincide.
    double forward = 0, backward = 0;
    for (const double fraction : {0.0, 0.25, 0.5, 0.75, 1.0}) {
        const auto p = curve.Value(firstA + (lastA - firstA) * fraction);
        forward += p.Distance(curveB.Value(firstB + (lastB - firstB) * fraction));
        backward += p.Distance(curveB.Value(lastB - (lastB - firstB) * fraction));
    }
    const bool reverse = forward > backward;
    const double lengthA = otherEdge ? GCPnts_AbscissaPoint::Length(curve) : 0;
    const double lengthB = otherEdge ? GCPnts_AbscissaPoint::Length(curveB) : 0;
    auto atLength = [&](BRepAdaptor_Curve &c, double length, double fraction) {
        GCPnts_AbscissaPoint point(c, length * fraction, c.FirstParameter());
        if (!point.IsDone()) fail(op, no_solution, "could not place an edge continuity sample");
        return point.Parameter();
    };
    std::vector<double> out;
    for (int i = 0; i < samples; ++i) {
        // Interior samples avoid undefined pole normals and corner derivatives at
        // endpoints.
        const double fraction = (i + 0.5) / samples;
        const double u = otherEdge ? atLength(curve, lengthA, fraction) :
                                     firstA + (lastA - firstA) * fraction;
        const double v = otherEdge ? atLength(curveB, lengthB, reverse ? 1 - fraction : fraction) : u;
        auto ja = surfaceJet(op, edge, a, u), jb = surfaceJet(op, edgeB, b, v);
        const double cosine = std::clamp(ja.n.Dot(jb.n), -1.0, 1.0);
        gp_Pnt p;
        gp_Vec t;
        curve.D1(u, p, t);
        t.Normalize();
        const auto [deviation, relative] = bendDifference(ja, jb, t);
        out.insert(out.end(), {p.X(), p.Y(), p.Z(), ja.p.Distance(jb.p), std::acos(cosine) * 180 / M_PI,
                               deviation, relative});
    }
    return out;
}
} // namespace BRepExt
