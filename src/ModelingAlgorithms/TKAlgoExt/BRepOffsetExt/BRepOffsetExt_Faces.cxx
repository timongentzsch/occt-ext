// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// offsetFaces: push-pull.
//
// Three ways, in this order:
//
// 1. BRepOffset_MakeOffset with a global offset of zero, a per-face offset (SetOffsetOnFace)
//    on the moved faces and GeomAbs_Intersection joins, which makes the unmoved neighbours
//    extend or trim to the moved faces. It is the engine of BRepOffsetAPI_MakeThickSolid, used
//    directly because the API class has no per-face offsets. It keeps the topology: a move
//    that removes a face or opens a hole makes it return no shape. And it returns
//    plausible-looking wrong solids in three situations, so each is checked:
//      - a moved face is tangent to a face that does not move with it (never tried),
//      - a concave face is offset past its centre of curvature and comes out inverted,
//      - a face travels past the far side of the body and the result intersects itself.
//    Measured on the same cases: a non-zero global offset with zero set on the other faces
//    hangs on a plain box, and GeomAbs_Arc joins fail with "command not done".
// 2. The exact reconstructions of offset_reconstruct.cpp for a few analytic cases.
// 3. The swept solid and a Boolean (push_face.cpp), which is what moves a face that is tangent
//    to a neighbour and what takes a move through any change of topology.
//
// Grow, for a face that is pulled out to a wall that faces the
// same way: None goes on past the wall by itself; Fixed stops where the face is flush with the
// wall; Moving goes on from there with the wall, the two as one face.
//
// A face of a sheet is moved by the first way alone.

#include <GeomAlgoExt_Cancel.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepToolsExt_History.hxx>
#include <BRepOffsetExt_Reconstruct.hxx>
#include <BRepOffsetExt_Solid.hxx>
#include <BRepToolsExt_Operations.hxx>

#include <BRepAdaptor_Surface.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <Bnd_Box.hxx>
#include <gp_Ax3.hxx>
#include <BRepAlgoAPI_Defeaturing.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepFilletExt_Untouched.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <BRepOffset_MakeOffset.hxx>
#include <BRep_Tool.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomLProp_SLProps.hxx>
#include <Precision.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>

#include <algorithm>
#include <cmath>
#include <optional>

namespace BRepExt {

namespace {

const char *const op = "offsetFaces";

/// A moved face and the neighbour that stays although it is tangent to it, if there is one.
std::optional<std::pair<TopoDS_Shape, TopoDS_Shape>> tangentNeighbour(const TopoDS_Shape &body,
                                                                       const FaceSet &moved) {
    Adjacency
        facesOfEdge;
    TopExp::MapShapesAndAncestors(body, TopAbs_EDGE, TopAbs_FACE, facesOfEdge);
    for (int i = 1; i <= facesOfEdge.Extent(); i++) {
        const TopoDS_Edge &edge = TopoDS::Edge(facesOfEdge.FindKey(i));
        const NCollection_List<TopoDS_Shape> &faces = facesOfEdge(i);
        if (faces.Size() != 2)
            continue;
        const TopoDS_Face &a = TopoDS::Face(faces.First());
        const TopoDS_Face &b = TopoDS::Face(faces.Last());
        if (moved.Contains(a) == moved.Contains(b) || !tangentAlong(edge, a, b))
            continue;
        return moved.Contains(a) ? std::pair<TopoDS_Shape, TopoDS_Shape>(a, b)
                                 : std::pair<TopoDS_Shape, TopoDS_Shape>(b, a);
    }
    return {};
}

/// Outward normal of `face` at (u, v), or false where the surface has none (an
/// apex, a pole).
bool normalOf(const TopoDS_Face &face, const Handle(Geom_Surface) & surface, double u, double v,
              gp_Dir &normal) {
    GeomLProp_SLProps props(surface, u, v, 1, Precision::Confusion());
    if (!props.IsNormalDefined())
        return false;
    normal = props.Normal();
    if (face.Orientation() == TopAbs_REVERSED)
        normal.Reverse();
    return true;
}

/// Whether `after` is `before` moved by `distance` along its outward normal,
/// facing the same way.
bool landed(const TopoDS_Face &before, const TopoDS_Face &after, double distance) {
    const BRepAdaptor_Surface bounds(before);
    const double u = 0.5 * (bounds.FirstUParameter() + bounds.LastUParameter());
    const double v = 0.5 * (bounds.FirstVParameter() + bounds.LastVParameter());
    const Handle(Geom_Surface) source = BRep_Tool::Surface(before);
    gp_Dir normal;
    if (!normalOf(before, source, u, v, normal))
        return true;
    const gp_Pnt expected = source->Value(u, v).Translated(gp_Vec(normal) * distance);

    const Handle(Geom_Surface) target = BRep_Tool::Surface(after);
    const GeomAPI_ProjectPointOnSurf projection(expected, target);
    // A gross check: the failure it guards against puts the surface a whole
    // radius away.
    if (projection.NbPoints() == 0 || projection.LowerDistance() > 0.01 * std::abs(distance)) {
        return false;
    }
    double targetU, targetV;
    projection.LowerDistanceParameters(targetU, targetV);
    gp_Dir targetNormal;
    if (!normalOf(after, target, targetU, targetV, targetNormal))
        return true;
    return normal.Dot(targetNormal) > 0.0;
}

/// Whether `result` is a body, and has more than `body` after a move outwards, less after one
/// inwards: what comes back inside out or as something else altogether has not.
bool gains(const TopoDS_Shape &body, const TopoDS_Shape &result, double distance) {
    const double after = estimatedVolume(result);
    return after > 0.0 && (after - estimatedVolume(body)) * distance > 0.0;
}

/// The first way. A null shape, and why, where the engine has no result that holds.
TopoDS_Shape byEngine(const TopoDS_Shape &body, const std::vector<TopoDS_Shape> &moved, double distance,
                      double tolerance, Lineage &faces, std::string &why) {
    TopoDS_Shape result;
    BRepOffset_MakeOffset maker;
    try {
        maker.Initialize(body, 0.0, tolerance, BRepOffset_Skin, false, false, GeomAbs_Intersection);
        for (const TopoDS_Shape &face : moved)
            maker.SetOffsetOnFace(TopoDS::Face(face), distance);
        maker.MakeOffsetShape(progress(true)->Start());
        checkCancelled("offsetFaces");
    } catch (const Standard_Failure &error) {
        checkCancelled("offsetFaces");
        why = "the neighbouring faces cannot be extended or trimmed to the moved faces";
        const char *detail = error.what();
        if (detail && *detail)
            why += std::string(": ") + detail;
        return {};
    }
    if (!maker.IsDone() || maker.Shape().IsNull()) {
        why = "the neighbouring faces cannot be extended or trimmed to the moved faces; the move "
              "probably removes a face or opens a hole";
        return {};
    }
    result = maker.Shape();
    if ((body.ShapeType() == TopAbs_SOLID && result.ShapeType() != TopAbs_SOLID) ||
        !BRepCheck_Analyzer(result).IsValid()) {
        why = "the result is not valid";
        return {};
    }
    // Every face is rebuilt, also the unmoved ones, so each input face has its successor as
    // an image (BRepOffset files it under Generated) or is gone.
    faces.images(maker);
    faces.within(result);
    for (const TopoDS_Shape &face : moved) {
        if (faces.of(face).empty()) {
            why = "a moved face vanished";
            return {};
        }
        for (const TopoDS_Shape &successor : faces.of(face))
            if (!landed(TopoDS::Face(face), TopoDS::Face(successor), distance)) {
                why = "a moved face came out inverted; the distance exceeds its radius of curvature";
                return {};
            }
    }
    if (intersectsItself(result) || (body.ShapeType() == TopAbs_SOLID && !gains(body, result, distance))) {
        why = "the result intersects itself; the faces were moved past another part of the body";
        return {};
    }
    return result;
}

struct Moved {
    TopoDS_Shape shape;
    Lineage faces;
};

/// `moved` of `body` along their normals by `distance`, by the first of the three ways that
/// has a result. Fails, naming the faces, where none has.
Moved moveFaces(const TopoDS_Shape &body, const std::vector<TopoDS_Shape> &moved, double distance,
                double tolerance) {
    FaceSet set;
    for (const TopoDS_Shape &face : moved)
        set.Add(face);
    const auto tangent = tangentNeighbour(body, set);
    const bool solid = body.ShapeType() == TopAbs_SOLID;
    if (tangent && !solid)
        fail(op, unsupported_geometry,
             "a face to move is tangent to a neighbour that stays; in a sheet they cannot be made "
             "to meet again. Move every face of the tangent chain",
             {tangent->first, tangent->second});

    if (!solid) {
        // A collapsed sheet can pass BRepCheck with zero area and no usable normal.
        TopoDS_Compound selected;
        BRep_Builder().MakeCompound(selected);
        for (const TopoDS_Shape &face : moved)
            BRep_Builder().Add(selected, face);
        const OffsetFate fate = classifyForOffset(selected, distance);
        if (!fate.vanish.empty() || !fate.partial.empty())
            fail(op, no_solution, "the distance reaches the radius of curvature of a moved sheet face; its offset collapses or folds",
                 fate.partial.empty() ? fate.vanish : fate.partial);
    }
    std::string why;
    if (!tangent) {
        Lineage faces(body);
        const TopoDS_Shape result = byEngine(body, moved, distance, tolerance, faces, why);
        if (!result.IsNull())
            return {result, faces};
        if (!solid)
            fail(op, no_solution, why, moved);
        Lineage rebuilt(body);
        const TopoDS_Shape exact = reconstructOffset(body, moved, distance, tolerance, rebuilt);
        if (!exact.IsNull() && gains(body, exact, distance))
            return {exact, rebuilt};
    }
    Lineage faces(body);
    std::string refused;
    std::vector<TopoDS_Shape> at;
    const TopoDS_Shape swept = sweepFaces(body, moved, distance, tolerance, faces, refused, at);
    if (!swept.IsNull())
        return {swept, faces};
    if (tangent)
        fail(op, refused.empty() ? no_solution : unsupported_geometry,
             "a face to move is tangent to a neighbour that stays (a fillet, typically), and " +
                 (refused.empty() ? std::string("the solid that it sweeps could not be joined to the body")
                                  : refused),
             at.empty() ? std::vector<TopoDS_Shape>{tangent->first, tangent->second} : at);
    fail(op, no_solution, refused.empty() ? why : why + "; and " + refused, at.empty() ? moved : at);
}

/// Whether two sets of faces have a face or an edge in common.
bool touching(const std::vector<TopoDS_Shape> &a, const std::vector<TopoDS_Shape> &b) {
    ShapeSet edges;
    for (const TopoDS_Shape &face : a)
        TopExp::MapShapes(face, TopAbs_EDGE, edges);
    for (const TopoDS_Shape &face : b)
        for (TopExp_Explorer of(face, TopAbs_EDGE); of.More(); of.Next())
            if (edges.Contains(of.Current()))
                return true;
    return false;
}

/// A flat face's outward plane, or false for another face.
bool outwardPlaneOf(const TopoDS_Shape &shape, gp_Pnt &point, gp_Dir &normal) {
    const TopoDS_Face face = TopoDS::Face(shape);
    const BRepAdaptor_Surface surface(face);
    return surface.GetType() == GeomAbs_Plane &&
           outwardNormal(face, surface.FirstUParameter(), surface.FirstVParameter(), point, normal);
}

} // namespace

EvolutionData offsetFaces(OcctKernel &kernel, uint32_t solidId, const std::vector<uint32_t> &faceIds,
                          double distance, double tolerance, const std::vector<int> &inputFaceHashes,
                          int hashUpperBound, int grow) {
    TopoDS_Shape body = get(kernel, solidId);
    if (body.ShapeType() == TopAbs_COMPOUND && !singleSolid(body).IsNull())
        body = singleSolid(body);
    if (body.ShapeType() != TopAbs_SOLID && body.ShapeType() != TopAbs_SHELL && body.ShapeType() != TopAbs_FACE)
        fail(op, invalid_input, "the shape is not a solid or a sheet");
    if (!std::isfinite(distance) || std::abs(distance) <= Precision::Confusion())
        fail(op, invalid_input, "the distance is zero");
    if (!std::isfinite(tolerance) || tolerance <= 0.0)
        fail(op, invalid_input, "the tolerance must be positive");
    if (grow < 0 || grow > 2)
        fail(op, invalid_input, "grow is 0 (none), 1 (fixed) or 2 (moving)");

    FaceHistory history(op, inputFaceHashes, hashUpperBound);
    std::vector<TopoDS_Shape> moved = subShapes(kernel, op, body, faceIds, TopAbs_FACE);

    const TopoDS_Shape original = body;
    std::optional<Lineage> prepared;
    if (distance > 0 && moved.size() == 1 && body.ShapeType() == TopAbs_SOLID) {
        gp_Pnt origin;
        gp_Dir normal;
        FaceSet selected;
        selected.Add(moved[0]);
        const auto tangent = tangentNeighbour(body, selected);
        if (tangent && outwardPlaneOf(moved[0], origin, normal))
            for (TopExp_Explorer all(body, TopAbs_FACE); all.More(); all.Next()) {
                gp_Pnt point;
                gp_Dir facing;
                if (all.Current().IsSame(moved[0]) || !outwardPlaneOf(all.Current(), point, facing) ||
                    facing.Dot(normal) < 1.0 - alongSlack)
                    continue;
                const double ahead = gp_Vec(origin, point).Dot(gp_Vec(normal));
                if (ahead <= tolerance) continue;
                const TopoDS_Face blend = TopoDS::Face(tangent->second);
                const BRepAdaptor_Surface support(blend);
                if (ahead > distance + tolerance) {
                    if (support.GetType() == GeomAbs_Cylinder && blend.Orientation() == TopAbs_REVERSED) {
                        // Measure the actual trimmed blend in the move's frame; its radius
                        // alone is not its reach for a corner of an arbitrary angle.
                        gp_Trsf frame;
                        frame.SetTransformation(gp_Ax3(origin, normal));
                        Bnd_Box extent;
                        BRepBndLib::AddOptimal(BRepBuilderAPI_Transform(blend, frame, true).Shape(), extent, false, false);
                        if (extent.IsVoid() || ahead < distance + extent.CornerMax().Z() - tolerance)
                            fail(op, unsupported_geometry,
                                 "the tangent blend would run past the end of the stationary face before the riser reaches the forward wall; partial trimming of that blend is not supported",
                                 {moved[0], blend, all.Current()});
                    }
                    continue;
                }
                // A concave cylindrical blend fills the sharp step's corner. Once the
                // riser reaches the forward wall that corner is wholly covered by the
                // sweep. Reconstruct the sharp step first, rather than translate a nub
                // past the wall. Convex blends cannot be discarded this way.
                if (support.GetType() != GeomAbs_Cylinder || blend.Orientation() != TopAbs_REVERSED)
                    fail(op, unsupported_geometry,
                         "the moved face reaches a forward wall with a tangent blend that cannot be trimmed flush by this reconstruction",
                         {moved[0], blend, all.Current()});
                BRepAlgoAPI_Defeaturing remove;
                remove.SetShape(body);
                remove.AddFaceToRemove(blend);
                remove.SetToFillHistory(true);
                remove.Build(progress(true)->Start());
                checkCancelled(op);
                const TopoDS_Shape sharp = remove.IsDone() && !remove.HasErrors() ? singleSolid(remove.Shape()) : TopoDS_Shape();
                if (sharp.IsNull() || !remove.IsDeleted(blend) || !BRepCheck_Analyzer(sharp).IsValid())
                    fail(op, no_solution, "the concave blend could not be removed to trim the step flush with the forward wall", {blend});
                prepared.emplace(body);
                prepared->through(remove);
                moved = prepared->of(moved[0]);
                if (moved.size() != 1)
                    fail(op, no_solution, "the riser did not remain one face after removing its concave blend");
                body = sharp;
                break;
            }
    }
    std::optional<Moved> done;
    gp_Pnt origin;
    gp_Dir normal;
    if (grow != 0 && distance > 0 && moved.size() == 1 && body.ShapeType() == TopAbs_SOLID &&
        outwardPlaneOf(moved[0], origin, normal)) {
        // The walls that face the same way, ahead of the face and within its reach, nearest first.
        std::vector<std::pair<double, TopoDS_Shape>> walls;
        for (TopExp_Explorer all(body, TopAbs_FACE); all.More(); all.Next()) {
            gp_Pnt point;
            gp_Dir facing;
            if (all.Current().IsSame(moved[0]) || !outwardPlaneOf(all.Current(), point, facing) ||
                facing.Dot(normal) < 1.0 - alongSlack)
                continue;
            const double ahead = gp_Vec(origin, point).Dot(gp_Vec(normal));
            if (ahead > tolerance && ahead <= distance + tolerance)
                walls.emplace_back(ahead, all.Current());
        }
        std::sort(walls.begin(), walls.end(), [](const auto &a, const auto &b) { return a.first < b.first; });
        for (const auto &[ahead, wall] : walls) {
            std::optional<Moved> flush;
            try {
                flush = moveFaces(body, moved, ahead, tolerance);
            } catch (const Failure &) {
                continue;
            }
            // Flush with it and beside it: this is the wall that the face has got to.
            const auto &face = flush->faces.of(moved[0]);
            const auto &beside = flush->faces.of(wall);
            const bool one = !face.empty() && !beside.empty() && face[0].IsSame(beside[0]);
            if (!one && !touching(face, beside))
                continue;
            if (grow == 2 && distance - ahead > tolerance) {
                std::vector<TopoDS_Shape> both = face;
                if (!one)
                    both.insert(both.end(), beside.begin(), beside.end());
                const Moved rest = moveFaces(flush->shape, both, distance - ahead, tolerance);
                flush->faces.then(rest.faces);
                flush->shape = rest.shape;
            }
            done = flush;
            break;
        }
    }
    if (!done)
        done = moveFaces(body, moved, distance, tolerance);
    if (prepared) {
        // Discarding the blend is justified only when the completed outward move covers
        // all of the original material, including the filled concave corner.
        Untouched<BRepAlgoAPI_Cut> missing(original, done->shape);
        if (!missing.IsDone() || missing.HasErrors() || estimatedVolume(missing.Shape()) > tolerance * tolerance * tolerance)
            fail(op, no_solution, "the move would lose material at the removed concave blend; it cannot be trimmed flush");
        prepared->then(done->faces);
        prepared->modified(history);
    } else {
        done->faces.modified(history);
    }
    return history.finish(store(kernel, done->shape));
}

} // namespace BRepExt
