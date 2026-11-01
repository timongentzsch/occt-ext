// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepToolsExt_Checks.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepAlgoExt_Prism.hxx>
#include <BRepAlgoExt_FaceBlend.hxx>
#include <BRepFilletExt_BlendChecks.hxx>
#include <optional>
#include <algorithm>
#include <BRepFilletExt_Untouched.hxx>
#include <BRepToolsExt_History.hxx>
#include <BRepAlgoExt_Contours.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <gp_Cylinder.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepFilletAPI_MakeChamfer.hxx>
#include <BRepPrimAPI_MakeHalfSpace.hxx>
#include <BRep_Tool.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Solid.hxx>
#include <cmath>
#include <limits>
#include <gp_Pln.hxx>

namespace BRepExt {
namespace {
std::pair<gp_Vec, gp_Vec> chamferDirections(const TopoDS_Face& a, const TopoDS_Face& b) {
    const gp_Vec na(orientedPlaneAxis(a)), nb(orientedPlaneAxis(b));
    auto da = na * na.Dot(nb) - nb, db = nb * nb.Dot(na) - na;
    if (da.Magnitude() < 1e-12 || db.Magnitude() < 1e-12)
        fail("chamfer", invalid_input, "coplanar faces have no bevel");
    return {da.Normalized(), db.Normalized()};
}

// OCCT's spline chamfer can be valid with a cap contact that is not the requested offset.
// Certify both contacts for a planar cap against vertical extrusion/cylinder walls.
bool chamferHasSetbacks(const TopoDS_Shape& shape, BRepFilletAPI_MakeChamfer& maker, double first, double second, const TopoDS_Shape& reference) {
    Adjacency beside;
    TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, beside);
    for (int contour = 1; contour <= maker.NbContours(); ++contour) {
        ShapeSet caps, walls, outline, newFaces;
        for (int e = 1; e <= maker.NbEdges(contour); ++e) {
            const auto edge = maker.Edge(contour, e);
            for (const auto& made : maker.Generated(edge))
                if (made.ShapeType() == TopAbs_FACE) newFaces.Add(made);
            for (TopExp_Explorer vertex(edge, TopAbs_VERTEX); vertex.More(); vertex.Next())
                for (const auto& made : maker.Generated(vertex.Current()))
                    if (made.ShapeType() == TopAbs_FACE) newFaces.Add(made);
            if (!beside.Contains(edge) || beside.FindFromKey(edge).Extent() != 2) continue;
            auto cap = TopoDS::Face(beside.FindFromKey(edge).First());
            auto wall = TopoDS::Face(beside.FindFromKey(edge).Last());
            if (BRepAdaptor_Surface(cap).GetType() != GeomAbs_Plane) std::swap(cap, wall);
            BRepAdaptor_Surface flat(cap), side(wall);
            if (flat.GetType() != GeomAbs_Plane) continue;
            const auto normal = flat.Plane().Axis().Direction();
            bool vertical = side.GetType() == GeomAbs_SurfaceOfExtrusion && normal.IsParallel(side.Direction(), BRepExt::tolerance::angle);
            if (side.GetType() == GeomAbs_Cylinder) vertical = normal.IsParallel(side.Cylinder().Axis().Direction(), BRepExt::tolerance::angle);
            if (!vertical) continue;
            caps.Add(cap); walls.Add(wall); outline.Add(edge);
        }
        for (int c = 1; c <= caps.Extent(); ++c) {
            const auto cap = TopoDS::Face(caps(c));
            const auto plane = BRepAdaptor_Surface(cap).Plane();
            const bool capFirst = reference.IsNull() || reference.IsSame(cap);
            const double capDistance = capFirst ? first : second, wallDistance = capFirst ? second : first;
            const auto check = [&](const TopoDS_Shape& support, bool onCap) {
                for (const auto& successor : maker.Modified(support)) {
                    for (TopExp_Explorer e(successor, TopAbs_EDGE); e.More(); e.Next()) {
                        bool contact = false;
                        for (int f = 1; f <= newFaces.Extent(); ++f)
                            for (TopExp_Explorer boundary(newFaces(f), TopAbs_EDGE); boundary.More(); boundary.Next())
                                contact |= e.Current().IsSame(boundary.Current());
                        if (!contact || BRep_Tool::Degenerated(TopoDS::Edge(e.Current()))) continue;
                        BRepAdaptor_Curve curve(TopoDS::Edge(e.Current()));
                        for (int sample = 1; sample < 16; ++sample) {
                            const auto p = curve.Value(curve.FirstParameter() + (curve.LastParameter() - curve.FirstParameter()) * sample / 16);
                            double setback = plane.Distance(p);
                            if (onCap) {
                                setback = std::numeric_limits<double>::infinity();
                                const auto point = BRepBuilderAPI_MakeVertex(p).Vertex();
                                for (int old = 1; old <= outline.Extent(); ++old) {
                                    BRepExtrema_DistShapeShape nearest(point, outline(old));
                                    if (!nearest.IsDone()) return false;
                                    setback = std::min(setback, nearest.Value());
                                }
                            }
                            const double distance = onCap ? capDistance : wallDistance;
                            if (std::abs(setback - distance) > 1e-6 * distance) return false;
                        }
                    }
                }
                return true;
            };
            if (!check(cap, true)) return false;
            for (int w = 1; w <= walls.Extent(); ++w)
                if (!check(walls(w), false)) return false;
        }
    }
    return true;
}
} // namespace

EvolutionData chamfer(OcctKernel &kernel, uint32_t shapeId, const std::vector<uint32_t> &edgeIds,
                      double distance1, double distance2, uint32_t referenceId,
                      const std::vector<int> &hashes, int bound, double angle, bool propagate) {
    const char *op = "chamfer";
    BooleanResult faced;
    const auto &shape = get(kernel, shapeId);
    if (!std::isfinite(distance1) || !std::isfinite(distance2) || distance1 <= 0 || distance2 <= 0 ||
        !std::isfinite(angle) || angle < 0 || angle >= M_PI / 2)
        fail(op, invalid_input, "chamfer distances must be positive");
    if (!referenceId && (distance1 != distance2 || angle > 0))
        fail(op, invalid_input, "asymmetric chamfer requires a reference face");
    auto edges = subShapes(kernel, op, shape, edgeIds, TopAbs_EDGE);
    edges.erase(std::remove_if(edges.begin(), edges.end(), [&](const auto& edge) {
        return tangentBlendEdge(shape, TopoDS::Edge(edge));
    }), edges.end());
    if (edges.empty())
        fail(op, unsupported_geometry, "the selected edges are tangent: there is nothing to bevel", {}, FailureKind::TangentEdge);
    if (!propagate && edges.size() != 1)
        fail(op, unsupported_geometry, "isolated chamfer requires one selected edge");
    if (propagate && distance1 == distance2 && angle == 0 && !referenceId)
        edges = blendContours<BRepFilletAPI_MakeChamfer>(shape, edges);
    Adjacency
        adjacency;
    TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, adjacency);
    // The bevel of a flat face's whole outline against walls square to it, also past the
    // tightest bend of the outline (blend_face.cpp): where OpenCascade's builder cannot.
    // Sets `vanished` where nothing of the face would be left.
    BRepFilletAPI_MakeChamfer maker(shape);
    bool vanished = false, wrongSetbacks = false, nativeFailed = false;
    std::optional<Failure> faceRefusal, nativeRefusal;
    const auto byFace = [&]() -> bool {
        if (!propagate || distance1 != distance2 || angle > 0 || referenceId) return false;
        try {
            // The whole runs of the chosen edges, as the builder has found them.
            ShapeSet selectedRuns;
            for (const auto& edge : edges) selectedRuns.Add(edge);
            for (int contour = 1; contour <= maker.NbContours(); ++contour)
                for (int i = 1; i <= maker.NbEdges(contour); ++i) selectedRuns.Add(maker.Edge(contour, i));
            std::vector<TopoDS_Shape> runs;
            for (int i = 1; i <= selectedRuns.Extent(); ++i) runs.push_back(selectedRuns(i));
            faced = faceBlend(op, shape, runs, distance1, true);
            return true;
        } catch (const Failure &refused) {
            faceRefusal = refused;
            vanished = refused.kind == FailureKind::CapConsumed;
            return false;
        }
    };
    std::vector<std::pair<TopoDS_Face, TopoDS_Face>> supports;
    for (const auto &selected : edges) {
        const auto &adjacent = adjacency.FindFromKey(selected);
        if (adjacent.Extent() != 2)
            fail(op, invalid_input, "each chamfer edge must have two adjacent faces");
        auto a = TopoDS::Face(adjacent.First()), b = TopoDS::Face(adjacent.Last());
        if (referenceId) {
            const auto &ref = get(kernel, referenceId);
            if (b.IsSame(ref))
                std::swap(a, b);
            else if (!a.IsSame(ref))
                fail(op, invalid_input, "reference face must border every selected edge");
        }
        supports.emplace_back(a, b);
        try {
            if (angle > 0) maker.AddDA(distance1, angle, TopoDS::Edge(selected), a);
            else maker.Add(distance1, distance2, TopoDS::Edge(selected), a);
        } catch (const Standard_Failure& error) {
            nativeFailed = true;
            nativeRefusal.emplace(op, no_solution, std::string("OCCT chamfer: ") + error.what(), std::vector<int>{}, "CONSTRUCTION_FAILED", FailureKind::NativeBuilder);
        }
    }
    if (!propagate) {
        // A triangle extruded only over the selected straight edge has normal end caps.
        // It never visits the tangent contour, even when that contour closes around a rim.
        const auto edge = TopoDS::Edge(edges.front());
        BRepAdaptor_Curve curve(edge);
        const auto [a, b] = supports.front();
        BRepAdaptor_Surface sa(a), sb(b);
        if (curve.GetType() != GeomAbs_Line || sa.GetType() != GeomAbs_Plane || sb.GetType() != GeomAbs_Plane)
            fail(op, unsupported_geometry, "isolated chamfer requires a straight edge between planar supports");
        if (angle > 0) fail(op, unsupported_geometry, "isolated distance-angle chamfer is not built");
        const auto [da, db] = chamferDirections(a, b);
        const auto first = curve.Value(curve.FirstParameter()), last = curve.Value(curve.LastParameter());
        const auto middle = curve.Value((curve.FirstParameter() + curve.LastParameter()) / 2);
        if (BRepClass3d_SolidClassifier(shape, middle.Translated((da + db) * (1e-3 * std::min(distance1, distance2))), BRepExt::tolerance::length).State() != TopAbs_IN)
            fail(op, unsupported_geometry, "isolated concave chamfer is not built");
        BRepBuilderAPI_MakePolygon triangle(first, first.Translated(da * distance1), first.Translated(db * distance2), true);
        const auto section = BRepBuilderAPI_MakeFace(triangle.Wire()).Face();
        const auto tool = BRepPrimAPI_MakePrism(section, gp_Vec(first, last)).Shape();
        Untouched<BRepAlgoAPI_Cut> cut(shape, tool);
        if (!cut.IsDone() || !BRepExt::invalidity(cut.Shape(), shape).empty())
            fail(op, no_solution, "isolated chamfer subtraction failed");
        auto solids = TopExp_Explorer(cut.Shape(), TopAbs_SOLID);
        if (!solids.More()) fail(op, no_solution, "isolated chamfer consumes the body");
        const auto result = solids.Current();
        solids.Next();
        if (solids.More()) fail(op, no_solution, "isolated chamfer splits the body");
        FaceHistory history(op, hashes, bound);
        history.collect(cut, shape);
        return history.finish(store(kernel, result));
    }
    // The chamfer walker mismeasures the setback when a planar cap contour
    // changes from an extrusion wall to a B-spline support (the surface made by
    // rounding two spline walls). BRepCheck still accepts that wrong geometry.
    // Inspect the complete propagated contour, not just the user's seed edge.
    for (int contour = 1; contour <= maker.NbContours(); ++contour) {
        bool extrusion = false, spline = false;
        for (int i = 1; i <= maker.NbEdges(contour); ++i) {
            const auto& adjacent = adjacency.FindFromKey(maker.Edge(contour, i));
            if (adjacent.Extent() != 2) continue;
            BRepAdaptor_Surface a(TopoDS::Face(adjacent.First())),
                                b(TopoDS::Face(adjacent.Last()));
            if (a.GetType() != GeomAbs_Plane && b.GetType() != GeomAbs_Plane) continue;
            const auto side = a.GetType() == GeomAbs_Plane ? b.GetType() : a.GetType();
            extrusion |= side == GeomAbs_SurfaceOfExtrusion;
            spline |= side == GeomAbs_BSplineSurface;
        }
        if (extrusion && spline && byFace()) {
            FaceHistory history(op, hashes, bound);
            history.collect(faced, shape);
            return history.finish(store(kernel, faced.shape));
        }
        if (extrusion && spline)
            fail(op, unsupported_geometry,
                 "chamfer across an extrusion/B-spline support transition is unsupported: "
                 "OCCT does not preserve the requested setback");
    }
    try {
        if (!nativeFailed) maker.Build();
        if (!nativeFailed && maker.IsDone() && BRepExt::invalidity(maker.Shape(), shape).empty()) {
            wrongSetbacks = angle == 0 && !chamferHasSetbacks(shape, maker, distance1, distance2,
                                                          referenceId ? get(kernel, referenceId) : TopoDS_Shape());
            if (!wrongSetbacks) {
                FaceHistory history(op, hashes, bound);
                history.collect(maker, shape);
                return history.finish(store(kernel, solidIfSingle(maker.Shape())));
            }
        }
    } catch (const Standard_Failure& error) {
        nativeRefusal.emplace(op, no_solution, std::string("OCCT chamfer: ") + error.what(), std::vector<int>{}, "CONSTRUCTION_FAILED", FailureKind::NativeBuilder);
    }
    if (byFace()) {
        FaceHistory history(op, hashes, bound);
        history.collect(faced, shape);
        return history.finish(store(kernel, faced.shape));
    }
    if (wrongSetbacks)
        fail(op, no_solution, "OCCT chamfer contacts do not preserve the requested setbacks; " +
             (faceRefusal ? faceRefusal->message : std::string("the face bevel does not support these distances or reference")));
    if (vanished) fail(op, no_solution, "requested distance consumes the planar cap: top face vanishes");
    if (angle == 0 && distance1 == distance2 &&
        blendConsumesPrismCap(shape, edges, distance1))
        fail(op, no_solution, "requested distance consumes the planar cap: top face vanishes");
    if (angle > 0) {
        if (nativeRefusal) fail(op, no_solution, nativeRefusal->message, edges, FailureKind::NativeBuilder);
        fail(op, no_solution, "distance-angle chamfer builder failed", edges, FailureKind::NativeBuilder);
    }
    // Face-consuming chamfers of convex polyhedra are intersections of half
    // spaces. OCCT's blend builder rejects a vanishing support face; its boolean
    // builder can remove it. Restrict the fallback to convex planar polyhedra: an
    // infinite bevel plane is incorrect on concave bodies or curved faces. Check
    // every face against every vertex.
    for (TopExp_Explorer f(shape, TopAbs_FACE); f.More(); f.Next()) {
        const auto face = TopoDS::Face(f.Current());
        BRepAdaptor_Surface s(face);
        if (s.GetType() != GeomAbs_Plane)
            fail(op, no_solution, "builder failed; plane-cut fallback requires a convex polyhedron");
        const auto plane = s.Plane();
        const gp_Vec n(orientedPlaneAxis(face));
        for (TopExp_Explorer v(shape, TopAbs_VERTEX); v.More(); v.Next())
            if (gp_Vec(plane.Location(), BRep_Tool::Pnt(TopoDS::Vertex(v.Current()))).Dot(n) > BRepExt::tolerance::length)
                fail(op, no_solution,
                     "builder failed; plane-cut fallback cannot change concave "
                     "topology");
    }
    NCollection_List<TopoDS_Shape> tools, args;
    args.Append(shape);
    for (size_t i = 0; i < edges.size(); ++i) {
        BRepAdaptor_Curve c(TopoDS::Edge(edges[i]));
        if (c.GetType() != GeomAbs_Line)
            fail(op, no_solution, "plane-cut fallback requires straight edges");
        const auto p = c.Value((c.FirstParameter() + c.LastParameter()) / 2);
        const auto [a, b] = supports[i];
        BRepAdaptor_Surface sa(a), sb(b);
        const auto [da, db] = chamferDirections(a, b);
        const auto p1 = p.Translated(da * distance1),
                   p2 = p.Translated(db * distance2);
        gp_Vec n = gp_Vec(p1, p2).Crossed(gp_Vec(c.Line().Direction()));
        tools.Append(
            BRepPrimAPI_MakeHalfSpace(BRepBuilderAPI_MakeFace(gp_Pln(p1, gp_Dir(n))).Face(), p).Solid());
    }
    BRepAlgoAPI_Cut cut;
    cut.SetArguments(args);
    cut.SetTools(tools);
    cut.SetNonDestructive(true);
    cut.Build();
    if (!cut.IsDone() || !BRepExt::invalidity(cut.Shape(), shape).empty() ||
        !TopExp_Explorer(cut.Shape(), TopAbs_SOLID).More())
        fail(op, no_solution, "support-plane chamfer did not leave a valid solid");
    FaceHistory history(op, hashes, bound);
    history.collect(cut, shape);
    return history.finish(store(kernel, solidIfSingle(cut.Shape())));
}
} // namespace BRepExt

namespace BRepExt {
EvolutionData clipPlane(OcctKernel &kernel, uint32_t shapeId, Vec3 origin, Vec3 normal, bool positive,
                        const std::vector<int> &hashes, int bound) {
    const auto &shape = get(kernel, shapeId);
    const gp_Pnt p(origin.x, origin.y, origin.z);
    gp_Vec n(normal.x, normal.y, normal.z);
    if (n.Magnitude() < 1e-15)
        fail("clipPlane", invalid_input, "plane normal must be nonzero");
    n.Normalize();
    const auto half = BRepPrimAPI_MakeHalfSpace(BRepBuilderAPI_MakeFace(gp_Pln(p, gp_Dir(n))).Face(),
                                                p.Translated(n * (positive ? 1 : -1)))
                          .Solid();
    Untouched<BRepAlgoAPI_Common> clip(shape, half);
    if (!clip.IsDone() || !BRepExt::invalidity(clip.Shape(), shape).empty())
        fail("clipPlane", no_solution, "plane clipping failed");
    FaceHistory history("clipPlane", hashes, bound);
    history.collect(clip, shape);
    return history.finish(store(kernel, clip.Shape()));
}
} // namespace BRepExt
