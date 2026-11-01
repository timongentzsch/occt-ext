// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepToolsExt_History.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAlgoAPI_Splitter.hxx>
#include <BRepBndLib.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgo_NormalProjection.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepLib.hxx>
#include <BRepTools.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <Geom2dAdaptor_Curve.hxx>
#include <Geom2dInt_GInter.hxx>
#include <Geom2d_Line.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomConvert.hxx>
#include <GeomConvert_CompCurveToBSplineCurve.hxx>
#include <GeomInt_IntSS.hxx>
#include <GeomLProp_SLProps.hxx>
#include <Geom_Plane.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <IntTools_CommonPrt.hxx>
#include <IntTools_EdgeEdge.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <NCollection_List.hxx>
#include <Precision.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <gp_Pln.hxx>
#include <algorithm>
#include <cmath>
#include <vector>

// Isoparam (edges on a face of a solid or sheet along the U or V direction), and Imprint with its Normal method, with what is hidden left out,
// and with its edges carried on to the edges of the face (Completion: Edge).

namespace BRepExt {
namespace {
EvolutionData splitBy(const char *op, OcctKernel &kernel, const TopoDS_Shape &shape,
                      const NCollection_List<TopoDS_Shape> &tools, const std::vector<int> &hashes, int bound,
                      double fuzzy = 0) {
    NCollection_List<TopoDS_Shape> arguments;
    arguments.Append(shape);
    BRepAlgoAPI_Splitter split;
    split.SetArguments(arguments);
    split.SetTools(tools);
    split.SetNonDestructive(true);
    if (fuzzy > 0)
        split.SetFuzzyValue(fuzzy);
    split.Build();
    if (!split.IsDone() || !BRepCheck_Analyzer(split.Shape()).IsValid())
        fail(op, no_solution, "the faces could not be split along the curves");
    FaceHistory history(op, hashes, bound);
    history.collect(split, shape);
    return history.finish(store(kernel, split.Shape()));
}
} // namespace

// Clip to the selected face before splitting its owner. Offset curves outside that face
// must not imprint another face of the body, or become detached curve objects.
EvolutionData splitFaceEdges(OcctKernel &kernel, uint32_t shapeId, uint32_t faceId,
                             const std::vector<uint32_t> &curveIds, double tolerance,
                             const std::vector<int> &hashes, int bound) {
    const char *op = "splitFaceEdges";
    const auto &shape = get(kernel, shapeId);
    const auto face = subShapes(kernel, op, shape, {faceId}, TopAbs_FACE)[0];
    if (curveIds.empty() || !std::isfinite(tolerance) || tolerance <= 0)
        fail(op, invalid_input, "curves and a positive tolerance are required");
    NCollection_List<TopoDS_Shape> tools;
    for (const auto id : curveIds) {
        NCollection_List<TopoDS_Shape> arguments, cutters;
        arguments.Append(face);
        cutters.Append(get(kernel, id));
        BRepAlgoAPI_Common clip;
        clip.SetArguments(arguments);
        clip.SetTools(cutters);
        clip.SetNonDestructive(true);
        clip.SetFuzzyValue(tolerance);
        clip.Build();
        if (!clip.IsDone())
            fail(op, no_solution, "the offset could not be clipped to the face");
        for (TopExp_Explorer it(clip.Shape(), TopAbs_EDGE); it.More(); it.Next())
            tools.Append(it.Current());
    }
    if (tools.IsEmpty())
        fail(op, no_solution, "the offset lies outside the selected face; no edges can be inserted");
    auto result = splitBy(op, kernel, shape, tools, hashes, bound);
    auto count = [](const TopoDS_Shape &s) {
        int n = 0;
        for (TopExp_Explorer it(s, TopAbs_FACE); it.More(); it.Next()) ++n;
        return n;
    };
    if (count(get(kernel, result.resultId)) <= count(shape)) {
        release(kernel, result.resultId);
        fail(op, no_solution, "the offset does not divide the selected face");
    }
    return result;
}

// Edges on a face of a body along lines of constant u or v of its surface, at `parameters`
// (shares of the face's extent in that direction): the face is split along them and stays
// part of its body. With `subdivide` each part gets a surface of its own, with its own control
// points: that is done for a sheet of one whole rectangular spline face, and refused otherwise.
EvolutionData isoparamEdges(OcctKernel &kernel, uint32_t shapeId, uint32_t faceId, bool alongU,
                            const std::vector<double> &parameters, bool subdivide, double tolerance,
                            const std::vector<int> &hashes, int bound) {
    const char *op = "isoparamEdges";
    const auto &shape = get(kernel, shapeId);
    const auto faces = subShapes(kernel, op, shape, {faceId}, TopAbs_FACE);
    const TopoDS_Face face = TopoDS::Face(faces[0]);
    if (parameters.empty() || !std::isfinite(tolerance) || tolerance <= 0)
        fail(op, invalid_input, "one or more parameters and a positive tolerance are required");
    std::vector<double> sorted = parameters;
    std::sort(sorted.begin(), sorted.end());
    for (const double f : sorted)
        if (!std::isfinite(f) || f <= 0 || f >= 1)
            fail(op, invalid_input, "parameters must lie strictly between 0 and 1");
    const auto surface = BRep_Tool::Surface(face);
    double u0, u1, v0, v1;
    BRepTools::UVBounds(face, u0, u1, v0, v1);
    if (subdivide) {
        int count = 0;
        for (TopExp_Explorer it(shape, TopAbs_FACE); it.More(); it.Next())
            ++count;
        auto base = surface;
        while (auto trimmed = Handle(Geom_RectangularTrimmedSurface)::DownCast(base))
            base = trimmed->BasisSurface();
        const auto spline = Handle(Geom_BSplineSurface)::DownCast(base);
        const auto rectangle = BRepBuilderAPI_MakeFace(surface, u0, u1, v0, v1, tolerance).Face();
        // Stored only to be measured, and taken out of the store again.
        const uint32_t whole = store(kernel, rectangle);
        const double all = surfaceArea(kernel, whole);
        release(kernel, whole);
        if (count != 1 || spline.IsNull() || spline->IsUPeriodic() || spline->IsVPeriodic() ||
            std::abs(surfaceArea(kernel, faceId) - all) > tolerance * std::max(1.0, all))
            fail(op, unsupported_geometry,
                 "Subdivide is for a sheet of one whole rectangular spline face; other faces get edges only");
        std::vector<double> cuts = {0};
        cuts.insert(cuts.end(), sorted.begin(), sorted.end());
        cuts.push_back(1);
        BRepBuilderAPI_Sewing sewing(tolerance);
        for (size_t i = 0; i + 1 < cuts.size(); ++i) {
            auto part = Handle(Geom_BSplineSurface)::DownCast(spline->Copy());
            const double a = alongU ? u0 + cuts[i] * (u1 - u0) : u0, b = alongU ? u0 + cuts[i + 1] * (u1 - u0) : u1,
                         c = alongU ? v0 : v0 + cuts[i] * (v1 - v0), d = alongU ? v1 : v0 + cuts[i + 1] * (v1 - v0);
            part->Segment(a, b, c, d);
            TopoDS_Face piece = BRepBuilderAPI_MakeFace(part, tolerance).Face();
            piece.Orientation(face.Orientation());
            sewing.Add(piece);
        }
        sewing.Perform();
        if (sewing.SewedShape().IsNull() || !BRepCheck_Analyzer(sewing.SewedShape()).IsValid())
            fail(op, no_solution, "the subdivided faces could not be joined");
        FaceHistory history(op, hashes, bound);
        if (history.tracks(face)) {
            NCollection_List<TopoDS_Shape> parts;
            for (TopExp_Explorer it(sewing.SewedShape(), TopAbs_FACE); it.More(); it.Next())
                parts.Append(it.Current());
            history.modified(face, parts);
        }
        return history.finish(store(kernel, sewing.SewedShape()));
    }
    // Of each line only what lies inside the face is drawn: the line is cut where it crosses
    // the face's edges, and of its pieces those are kept whose middle is in the face. Drawn
    // across the face's whole extent it would also divide whatever else of the body lies on
    // the same surface there, such as a face inside a hole of this one.
    BRepTopAdaptor_FClass2d inside(face, Precision::PConfusion());
    const double low = alongU ? v0 : u0, high = alongU ? v1 : u1;
    // Two lines, or a line and the end of the face's extent, that are nowhere further apart
    // than the tolerance are one line: nothing would lie between them.
    auto apart = [&](double f, double g) {
        double most = 0;
        for (const double share : {0.0, 0.25, 0.5, 0.75, 1.0}) {
            const double other = low + (high - low) * share;
            const double a = alongU ? u0 + f * (u1 - u0) : v0 + f * (v1 - v0),
                         b = alongU ? u0 + g * (u1 - u0) : v0 + g * (v1 - v0);
            most = std::max(most, (alongU ? surface->Value(a, other) : surface->Value(other, a))
                                      .Distance(alongU ? surface->Value(b, other) : surface->Value(other, b)));
        }
        return most > tolerance;
    };
    for (size_t i = 0; i < sorted.size(); ++i) {
        if (!apart(0, sorted[i]) || !apart(sorted[i], 1))
            fail(op, invalid_input,
                 "a line at " + std::to_string(sorted[i]) + " lies within the tolerance of the end of the face");
        if (i > 0 && !apart(sorted[i - 1], sorted[i]))
            fail(op, invalid_input, "two of the lines lie within the tolerance of each other");
    }
    NCollection_List<TopoDS_Shape> tools;
    for (const double f : sorted) {
        const double at = alongU ? u0 + f * (u1 - u0) : v0 + f * (v1 - v0);
        const auto curve = alongU ? surface->UIso(at) : surface->VIso(at);
        const Handle(Geom2d_Line) line = new Geom2d_Line(alongU ? gp_Pnt2d(at, 0) : gp_Pnt2d(0, at),
                                                          alongU ? gp_Dir2d(0, 1) : gp_Dir2d(1, 0));
        const Geom2dAdaptor_Curve whole(line, low, high);
        std::vector<double> stops = {low, high};
        for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next()) {
            const TopoDS_Edge edge = TopoDS::Edge(it.Current());
            double a, b;
            const auto onFace = BRep_Tool::CurveOnSurface(edge, face, a, b);
            if (onFace.IsNull())
                continue;
            Geom2dInt_GInter crossings(whole, Geom2dAdaptor_Curve(onFace, a, b), Precision::PConfusion(),
                                       Precision::PConfusion());
            if (!crossings.IsDone())
                fail(op, no_solution, "a line of the surface could not be cut at the face's edges");
            for (int k = 1; k <= crossings.NbPoints(); ++k)
                stops.push_back(crossings.Point(k).ParamOnFirst());
            // An edge that the line runs along bounds nothing of it that is inside.
            for (int k = 1; k <= crossings.NbSegments(); ++k) {
                stops.push_back(crossings.Segment(k).FirstPoint().ParamOnFirst());
                stops.push_back(crossings.Segment(k).LastPoint().ParamOnFirst());
            }
        }
        std::sort(stops.begin(), stops.end());
        bool any = false;
        for (size_t i = 0; i + 1 < stops.size(); ++i) {
            const double a = stops[i], b = stops[i + 1], mid = (a + b) / 2;
            if (b - a <= Precision::PConfusion() ||
                inside.Perform(alongU ? gp_Pnt2d(at, mid) : gp_Pnt2d(mid, at)) != TopAbs_IN)
                continue;
            if (curve->Value(a).Distance(curve->Value(mid)) + curve->Value(mid).Distance(curve->Value(b)) <=
                tolerance)
                continue;
            BRepBuilderAPI_MakeEdge edge(curve, a, b);
            if (!edge.IsDone())
                fail(op, no_solution, "a line of the surface could not be drawn");
            tools.Append(edge.Edge());
            any = true;
        }
        if (!any)
            fail(op, invalid_input,
                 "a line at " + std::to_string(f) +
                     " does not cross the inside of the face: it lies on the face's edge or outside it");
    }
    return splitBy(op, kernel, shape, tools, hashes, bound);
}

// Imprint: `method` 0 projects each curve along `direction` (`both`: and against it), 1 along
// the normals of the faces (onto the nearest of the body). `hideOccluded` keeps of a
// projection along a direction only what has nothing of the body between it and the curve.
// `completion` 0 leaves the projected edges as they are; 1 carries each loose end on, straight
// in the parameters of the face it lies on, until it leaves that face (it stretches projected
// edges until they contact other edges). A third completion, to the boundary, is not built.
EvolutionData imprintCurves(OcctKernel &kernel, uint32_t shapeId, const std::vector<uint32_t> &curveIds,
                            int method, Vec3 direction, bool both, bool hideOccluded, int completion,
                            const std::vector<int> &hashes, int bound) {
    const char *op = "imprintCurves";
    const auto &shape = get(kernel, shapeId);
    if (method < 0 || method > 1 || completion < 0 || completion > 1 || curveIds.empty())
        fail(op, invalid_input, "curves, a method and a completion are required");
    std::vector<TopoDS_Edge> projected;
    for (const auto id : curveIds) {
        if (method == 1) {
            // Onto the nearest of the body: what lies further from the curve than the curve ever
            // is from the body (at 32 places to each edge, and a hundredth more) is left out, so
            // that a curve on a face is not also drawn on the faces behind it.
            double furthest = 0;
            DistanceTo toShape(shape);
            for (TopExp_Explorer e(get(kernel, id), TopAbs_EDGE); e.More(); e.Next()) {
                BRepAdaptor_Curve c(TopoDS::Edge(e.Current()));
                for (int k = 0; k <= 32; ++k)
                    if (const auto gap = toShape.from(c.Value(c.FirstParameter() + (c.LastParameter() - c.FirstParameter()) * k / 32)))
                        furthest = std::max(furthest, *gap);
            }
            BRepAlgo_NormalProjection normal(shape);
            normal.SetMaxDistance(furthest * 1.01 + 10 * Precision::Confusion());
            normal.Add(get(kernel, id));
            normal.Build();
            if (!normal.IsDone())
                fail(op, no_solution, "a curve could not be projected along the faces' normals");
            for (TopExp_Explorer it(normal.Projection(), TopAbs_EDGE); it.More(); it.Next())
                // A point's worth of projection, where the curve only touches a face's edge,
                // divides nothing.
                if (GCPnts_AbscissaPoint::Length(BRepAdaptor_Curve(TopoDS::Edge(it.Current()))) >
                    10 * std::max(BRep_Tool::Tolerance(TopoDS::Edge(it.Current())), Precision::Confusion()))
                    projected.push_back(TopoDS::Edge(it.Current()));
        } else {
            const auto made = hideOccluded ? projectCurveVisible(kernel, id, shapeId, direction, both)
                                           : projectCurve(kernel, id, shapeId, direction, both);
            for (TopExp_Explorer it(get(kernel, made.resultId), TopAbs_EDGE); it.More(); it.Next())
                projected.push_back(TopoDS::Edge(it.Current()));
            // The projection is only a step on the way: it is not left in the store.
            release(kernel, made.resultId);
        }
    }
    if (projected.empty())
        fail(op, no_solution, "nothing of the curves falls on the body");
    NCollection_List<TopoDS_Shape> tools;
    for (const auto &edge : projected)
        tools.Append(edge);
    if (completion == 1) {
        // The ends of the projected edges, and how many edges each is an end of. Two ends are
        // one where they are no further apart than the tolerances their edges were made with.
        struct End {
            gp_Pnt at;
            gp_Vec leaving; // the edge's direction there, out of the edge
            double tolerance;
            int uses;
        };
        std::vector<End> ends;
        for (const auto &edge : projected) {
            BRepAdaptor_Curve c(edge);
            const double a = c.FirstParameter(), b = c.LastParameter();
            const double tolerance = std::max(BRep_Tool::Tolerance(edge), Precision::Confusion());
            for (const bool first : {true, false}) {
                gp_Pnt at;
                gp_Vec leaving;
                c.D1(first ? a : b, at, leaving);
                if (first)
                    leaving.Reverse();
                bool known = false;
                for (auto &end : ends)
                    if (end.at.Distance(at) <= end.tolerance + tolerance) {
                        ++end.uses;
                        known = true;
                    }
                if (!known)
                    ends.push_back({at, leaving, tolerance, 1});
            }
        }
        for (const auto &end : ends) {
            if (end.uses != 1)
                continue;
            for (TopExp_Explorer it(shape, TopAbs_FACE); it.More(); it.Next()) {
                const TopoDS_Face face = TopoDS::Face(it.Current());
                const auto surface = BRep_Tool::Surface(face);
                ShapeAnalysis_Surface on(surface);
                gp_Pnt2d uv = on.ValueOfUV(end.at, Precision::Confusion());
                if (on.Gap() > end.tolerance)
                    continue;
                // On a surface that closes on itself the point is taken where the face has it.
                double u0, u1, v0, v1;
                BRepTools::UVBounds(face, u0, u1, v0, v1);
                gp_Vec2d turns(0, 0);
                if (surface->IsUPeriodic())
                    turns.SetX(surface->UPeriod() * std::ceil((u0 - uv.X()) / surface->UPeriod() - 1e-9));
                if (surface->IsVPeriodic())
                    turns.SetY(surface->VPeriod() * std::ceil((v0 - uv.Y()) / surface->VPeriod() - 1e-9));
                uv.Translate(turns);
                BRepTopAdaptor_FClass2d inside(face, Precision::Confusion());
                // An end on the face's own boundary is where it should be.
                if (inside.Perform(uv) != TopAbs_IN)
                    continue;
                // The way on in the face's parameters: the one whose direction on the surface
                // is nearest the edge's own as it leaves its end.
                gp_Pnt here;
                gp_Vec du, dv;
                surface->D1(uv.X(), uv.Y(), here, du, dv);
                const double e = du.Dot(du), f = du.Dot(dv), g = dv.Dot(dv), det = e * g - f * f;
                if (end.leaving.Magnitude() < 1e-14 || det <= 1e-14 * e * g)
                    fail(op, no_solution,
                         "a projected edge has no direction in its face's parameters to be carried on along");
                const double p = end.leaving.Dot(du), q = end.leaving.Dot(dv);
                gp_Vec2d way((g * p - f * q) / det, (e * q - f * p) / det);
                way.Normalize();
                // The first edge of the face that the straight way meets: where a line from
                // the end crosses the edges' own curves in the face's parameters, a seam among
                // them. (On a face that closes on itself the way never leaves the face.)
                const double reach = 2 * std::hypot(u1 - u0, v1 - v0);
                const Handle(Geom2d_Line) line = new Geom2d_Line(uv, gp_Dir2d(way));
                const Geom2dAdaptor_Curve ray(line, 0, reach);
                double hi = 0;
                for (TopExp_Explorer e(face, TopAbs_EDGE); e.More(); e.Next()) {
                    double a, b;
                    const auto onFace = BRep_Tool::CurveOnSurface(TopoDS::Edge(e.Current()), face, a, b);
                    if (onFace.IsNull())
                        continue;
                    Geom2dInt_GInter crossings(ray, Geom2dAdaptor_Curve(onFace, a, b), Precision::PConfusion(),
                                               Precision::PConfusion());
                    if (!crossings.IsDone())
                        continue;
                    for (int k = 1; k <= crossings.NbPoints(); ++k) {
                        const double t = crossings.Point(k).ParamOnFirst();
                        if (t > Precision::PConfusion() && (hi == 0 || t < hi))
                            hi = t;
                    }
                }
                if (hi == 0)
                    fail(op, no_solution, "a projected edge could not be carried on to an edge of its face");
                // Along a direction the edge carried on is the projection of the curve carried
                // on: it stays in the plane of the edge's direction and the direction of
                // projection. Straight in the parameters of a face is that on a plane, along
                // or round a cylinder, round a sphere; beside the pole of a sphere it is
                // not, and a wrong edge is not drawn. The direction of an approximated
                // projected edge is known to about a ten-thousandth, so that much of the
                // carried length is allowed.
                if (method == 0) {
                    const gp_Vec across = end.leaving.Crossed(gp_Vec(direction.x, direction.y, direction.z));
                    double off = 0, length = 0;
                    gp_Pnt before = end.at;
                    for (int k = 1; k <= 16; ++k) {
                        const gp_Pnt2d q = uv.Translated(way * (hi * k / 16));
                        const gp_Pnt p = surface->Value(q.X(), q.Y());
                        length += before.Distance(p);
                        before = p;
                        if (across.Magnitude() > 1e-14)
                            off = std::max(off, std::abs(gp_Vec(end.at, p).Dot(across.Normalized())));
                    }
                    if (off > end.tolerance + 1e-4 * length)
                        fail(op, unsupported_geometry,
                             "Completion Edge: carried on straight in the parameters of its face a projected "
                             "edge would leave the projection (the face's parameters bend there, as beside "
                             "the pole of a sphere)");
                }
                BRepBuilderAPI_MakeEdge carried(line, surface, 0, hi);
                if (!carried.IsDone())
                    fail(op, no_solution, "a projected edge could not be carried on to an edge of its face");
                TopoDS_Edge edge = carried.Edge();
                BRepLib::BuildCurves3d(edge);
                tools.Append(edge);
                break;
            }
        }
    }
    return splitBy(op, kernel, shape, tools, hashes, bound);
}

// The Knife: the faces of a body split along curves that lie on them, and only those faces.
// Every edge of the curves has to lie on the body to the tolerance (at 16 places to each edge
// and its ends); an edge that only touches a face's boundary at a point divides nothing of it.
EvolutionData splitOnFaces(OcctKernel &kernel, uint32_t shapeId, const std::vector<uint32_t> &curveIds,
                           double tolerance, const std::vector<int> &hashes, int bound) {
    const char *op = "splitOnFaces";
    const auto &shape = get(kernel, shapeId);
    if (curveIds.empty() || !std::isfinite(tolerance) || tolerance <= 0)
        fail(op, invalid_input, "curves and a positive tolerance are required");
    NCollection_List<TopoDS_Shape> tools;
    DistanceTo toShape(shape);
    for (const auto id : curveIds)
        for (TopExp_Explorer e(get(kernel, id), TopAbs_EDGE); e.More(); e.Next()) {
            const TopoDS_Edge edge = TopoDS::Edge(e.Current());
            BRepAdaptor_Curve c(edge);
            for (int k = 0; k <= 16; ++k) {
                const gp_Pnt p = c.Value(c.FirstParameter() + (c.LastParameter() - c.FirstParameter()) * k / 16);
                const auto gap = toShape.from(p);
                if (!gap || *gap > tolerance)
                    fail(op, invalid_input,
                         "a curve does not lie on the body: it is " + std::to_string(gap ? *gap : -1.0) +
                             " from it somewhere");
            }
            tools.Append(edge);
        }
    if (tools.IsEmpty())
        fail(op, invalid_input, "the curves have no edges");
    return splitBy(op, kernel, shape, tools, hashes, bound);
}

// Complete Edge: "Extend the selected Edge until it intersects with another Edge". Each loose
// end of the selected edges, an end that no other edge of the body meets, runs on across its
// face until it meets an edge of the face, and the body is split along the runs. A run lies
// where the face's surface meets the plane of the edge's direction at its end and the
// surface's normal there: on a plane the straight line the edge leaves along, round a cylinder
// an ellipse through its axis. Two runs that cross both stop where they cross.
EvolutionData completeEdges(OcctKernel &kernel, uint32_t shapeId, const std::vector<uint32_t> &edgeIds,
                            double tolerance, const std::vector<int> &hashes, int bound) {
    const char *op = "completeEdges";
    const auto &shape = get(kernel, shapeId);
    if (edgeIds.empty() || !std::isfinite(tolerance) || tolerance <= 0)
        fail(op, invalid_input, "edges and a positive tolerance are required");
    const auto selected = subShapes(kernel, op, shape, edgeIds, TopAbs_EDGE);
    NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher> edgesAt,
        facesOf;
    TopExp::MapShapesAndUniqueAncestors(shape, TopAbs_VERTEX, TopAbs_EDGE, edgesAt);
    TopExp::MapShapesAndUniqueAncestors(shape, TopAbs_EDGE, TopAbs_FACE, facesOf);
    struct Run {
        TopoDS_Edge edge;
        TopoDS_Face face;
        Handle(Geom_Curve) way;
        double from, to;
        double near; // how far along it a meeting is still its own start
        // How far along the run `t` is.
        double travelled(double t) const { return GCPnts_AbscissaPoint::Length(GeomAdaptor_Curve(way), from, t); }
    };
    std::vector<Run> runs;
    double fuzzy = tolerance;
    for (const auto &each : selected) {
        const TopoDS_Edge edge = TopoDS::Edge(each);
        TopoDS_Vertex first, last;
        TopExp::Vertices(edge, first, last);
        if (BRep_Tool::Degenerated(edge) || first.IsSame(last))
            continue;
        BRepAdaptor_Curve along(edge);
        for (const auto &vertex : {first, last}) {
            if (edgesAt.FindFromKey(vertex).Extent() != 1)
                continue;
            const auto &faces = facesOf.FindFromKey(edge);
            if (faces.Extent() != 1)
                fail(op, invalid_input, "an edge with a loose end does not lie inside one face", {edge});
            const TopoDS_Face face = TopoDS::Face(faces.First());
            // The edge's direction at this end, out of the edge.
            const double t = BRep_Tool::Parameter(vertex, edge);
            gp_Pnt on;
            gp_Vec leaving;
            along.D1(t, on, leaving);
            if (std::abs(t - along.FirstParameter()) < std::abs(t - along.LastParameter()))
                leaving.Reverse();
            const gp_Pnt start = BRep_Tool::Pnt(vertex);
            const auto surface = BRep_Tool::Surface(face);
            const gp_Pnt2d uv = ShapeAnalysis_Surface(surface).ValueOfUV(start, tolerance);
            GeomLProp_SLProps at(surface, uv.X(), uv.Y(), 1, Precision::Confusion());
            const gp_Vec across = at.IsNormalDefined() ? leaving.Crossed(gp_Vec(at.Normal())) : gp_Vec();
            if (across.Magnitude() <= Precision::Confusion() * leaving.Magnitude())
                fail(op, unsupported_geometry, "the face has no normal at a loose end of an edge to run on along",
                     {edge});
            // The way on: of where the surface meets the plane, the branch through the end. The
            // meeting comes in pieces, cut where the surface's parameters close: a piece of a
            // closed curve is taken whole, and pieces of a spline are joined end to end.
            GeomInt_IntSS meet(surface, new Geom_Plane(gp_Pln(start, gp_Dir(across))), tolerance);
            const double reached = std::max(meet.IsDone() ? meet.TolReached3d() : 0.0, BRep_Tool::Tolerance(vertex));
            std::vector<Handle(Geom_Curve)> pieces;
            for (int i = 1; meet.IsDone() && i <= meet.NbLines(); ++i) {
                Handle(Geom_Curve) piece = meet.Line(i);
                const auto trimmed = Handle(Geom_TrimmedCurve)::DownCast(piece);
                if (!trimmed.IsNull() && trimmed->BasisCurve()->IsPeriodic())
                    piece = trimmed->BasisCurve();
                pieces.push_back(piece);
            }
            Handle(Geom_Curve) way;
            double from = 0, gap = Precision::Infinite();
            for (const auto &piece : pieces) {
                GeomAPI_ProjectPointOnCurve foot(start, piece);
                if (foot.NbPoints() > 0 && foot.LowerDistance() < gap) {
                    gap = foot.LowerDistance();
                    way = piece;
                    from = foot.LowerDistanceParameter();
                }
            }
            if (way.IsNull() || gap > tolerance + reached)
                fail(op, no_solution, "the face's surface could not be followed on from a loose end of an edge", {edge});
            if (!way->IsPeriodic() && !Precision::IsInfinite(way->LastParameter())) {
                const gp_Pnt at = way->Value(from);
                GeomConvert_CompCurveToBSplineCurve joined(GeomConvert::CurveToBSplineCurve(way));
                for (bool more = true; more;) {
                    more = false;
                    for (auto &piece : pieces)
                        if (!piece.IsNull() && piece != way && !piece->IsPeriodic() &&
                            joined.Add(GeomConvert::CurveToBSplineCurve(piece), tolerance + reached)) {
                            piece.Nullify();
                            more = true;
                        }
                }
                // Joined round to where they began, the pieces are a closed curve, gone round
                // from any point of it.
                const Handle(Geom_BSplineCurve) whole = joined.BSplineCurve();
                if (whole->StartPoint().Distance(whole->EndPoint()) <= tolerance + reached)
                    whole->SetPeriodic();
                way = whole;
                from = GeomAPI_ProjectPointOnCurve(at, way).LowerDistanceParameter();
            }
            fuzzy = std::max(fuzzy, gap);
            gp_Pnt p;
            gp_Vec d;
            way->D1(from, p, d);
            if (d.Dot(leaving) < 0) {
                from = way->ReversedParameter(from);
                way = way->Reversed();
                way->D1(from, p, d);
            }
            // As far as the way goes: once round if it closes, else across the face and more.
            double to = way->LastParameter();
            if (way->IsPeriodic())
                to = from + way->Period();
            else if (Precision::IsInfinite(to)) {
                Bnd_Box box;
                BRepBndLib::Add(face, box);
                to = from + 2 * std::sqrt(box.SquareExtent()) / d.Magnitude();
            }
            if (to - from <= Precision::PConfusion())
                fail(op, no_solution, "a loose end of an edge lies on the end of its face's surface", {edge});
            runs.push_back({edge, face, way, from, to, 10 * (tolerance + reached)});
        }
    }
    if (runs.empty())
        fail(op, invalid_input,
             "none of the selected edges has a loose end: other edges meet each of them at both ends", selected);
    // Where two edges meet, as the Boolean operations find it.
    auto meetings = [](const TopoDS_Edge &a, const TopoDS_Edge &b, double fuzz, auto &&take) {
        IntTools_EdgeEdge meeting(a, b);
        meeting.SetFuzzyValue(fuzz);
        meeting.Perform();
        if (!meeting.IsDone())
            return;
        for (const auto &part : meeting.CommonParts())
            take(part);
    };
    auto runEdge = [](const Run &run, double to) { return BRepBuilderAPI_MakeEdge(run.way, run.from, to).Edge(); };
    // Each run stops at the first edge of its face it meets. Seams are where the face's
    // parameters close, not edges to stop at.
    for (auto &run : runs) {
        const TopoDS_Edge ray = runEdge(run, run.to);
        double hit = Precision::Infinite();
        bool own = false;
        for (TopExp_Explorer it(run.face, TopAbs_EDGE); it.More(); it.Next()) {
            const TopoDS_Edge other = TopoDS::Edge(it.Current());
            if (BRep_Tool::Degenerated(other) || BRep_Tool::IsClosed(other, run.face))
                continue;
            meetings(ray, other, tolerance, [&](const IntTools_CommonPrt &part) {
                const double t = part.Type() == TopAbs_VERTEX ? part.VertexParameter1() : part.Range1().First();
                if (t < hit && run.travelled(t) > run.near) {
                    hit = t;
                    own = other.IsSame(run.edge);
                }
            });
        }
        if (Precision::IsInfinite(hit))
            fail(op, no_solution, "a loose end of an edge runs on without meeting an edge of its face", {run.edge});
        // A way that closes on itself comes round to the end it started from, along the edge
        // there: it would touch its own edge rather than meet another.
        if (own && run.way->IsPeriodic())
            fail(op, unsupported_geometry,
                 "a loose end of an edge runs round its face and back into the edge itself without meeting "
                 "another edge",
                 {run.edge});
        run.to = hit;
    }
    // Runs across one face that cross stop where they cross, both of them.
    std::vector<double> ends;
    for (const auto &run : runs)
        ends.push_back(run.to);
    for (size_t i = 0; i < runs.size(); ++i)
        for (size_t j = i + 1; j < runs.size(); ++j) {
            if (!runs[i].face.IsSame(runs[j].face))
                continue;
            meetings(runEdge(runs[i], runs[i].to), runEdge(runs[j], runs[j].to), tolerance,
                     [&](const IntTools_CommonPrt &part) {
                         if (part.Type() != TopAbs_VERTEX)
                             return;
                         const double a = part.VertexParameter1(), b = part.VertexParameter2();
                         if (runs[i].travelled(a) > runs[i].near)
                             ends[i] = std::min(ends[i], a);
                         if (runs[j].travelled(b) > runs[j].near)
                             ends[j] = std::min(ends[j], b);
                     });
        }
    NCollection_List<TopoDS_Shape> tools;
    for (size_t i = 0; i < runs.size(); ++i)
        tools.Append(runEdge(runs[i], ends[i]));
    // The runs begin where the surface and the plane were found to meet, which is as near the
    // ends of the edges as that was found.
    return splitBy(op, kernel, shape, tools, hashes, bound, fuzzy);
}
} // namespace BRepExt
