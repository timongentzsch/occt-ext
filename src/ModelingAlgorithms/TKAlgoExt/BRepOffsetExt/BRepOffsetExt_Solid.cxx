// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// The offset of a whole body, with the faces that it cannot keep taken out first.
//
// OpenCascade removes faces that an offset shrinks to nothing only where every face is flat
// (BRepOffset_MakeOffset.cxx, `myIsPlanar`): a fillet whose radius the distance reaches comes
// back turned inside out, as an invalid solid or not at all. So the faces are classified first
// (offset_parts.h), those that vanish are removed from a working copy by extending their
// neighbours (BRepAlgoAPI_Defeaturing), and that copy is offset. What is bent too tightly in a
// part of a face only cannot be removed like that and is refused, naming the face.

#include <GeomAlgoExt_Cancel.hxx>
#include <BRepOffsetExt_Solid.hxx>
#include <BRepToolsExt_Operations.hxx>

#include <BRepAlgoAPI_Defeaturing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepOffsetAPI_MakeOffsetShape.hxx>
#include <BRepOffset_MakeOffset.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Surface.hxx>
#include <ShapeBuild_ReShape.hxx>
#include <ShapeUpgrade_ShapeDivideContinuity.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>

#include <cmath>

namespace BRepExt {

TopoDS_Shape dividedAtC0(const TopoDS_Shape &body, Lineage &faces) {
    bool any = false;
    for (TopExp_Explorer all(body, TopAbs_FACE); all.More() && !any; all.Next())
        any = BRep_Tool::Surface(TopoDS::Face(all.Current()))->Continuity() == GeomAbs_C0;
    if (!any)
        return body;
    ShapeUpgrade_ShapeDivideContinuity divide(body);
    divide.SetSurfaceCriterion(GeomAbs_C1);
    divide.SetBoundaryCriterion(GeomAbs_C0);
    divide.SetPCurveCriterion(GeomAbs_C0);
    if (!divide.Perform())
        return body;
    faces.throughHistory(divide.GetContext()->History());
    return divide.Result();
}

TopoDS_Shape withoutVanishing(const char *op, const TopoDS_Shape &solid, double distance, Lineage &faces) {
    const OffsetFate fate = classifyForOffset(solid, distance);
    if (!fate.partial.empty())
        fail(op, unsupported_geometry,
             "a face is bent more tightly than the distance in a part of it; its offset would cross "
             "itself there",
             fate.partial);
    if (fate.vanish.empty())
        return solid;
    BRepAlgoAPI_Defeaturing remove;
    remove.SetShape(solid);
    for (const TopoDS_Shape &face : fate.vanish)
        remove.AddFaceToRemove(face);
    remove.SetToFillHistory(true);
    remove.Build(progress(true)->Start());
    checkCancelled(op);
    std::vector<TopoDS_Shape> left;
    if (remove.IsDone() && !remove.HasErrors())
        for (const TopoDS_Shape &face : fate.vanish)
            if (!remove.IsDeleted(face))
                left.push_back(face);
    const TopoDS_Shape without = remove.IsDone() && !remove.HasErrors() ? singleSolid(remove.Shape()) : TopoDS_Shape();
    if (without.IsNull() || !left.empty() || !BRepCheck_Analyzer(without).IsValid())
        fail(op, no_solution,
             "faces are bent more tightly than the distance, and the faces around them cannot be "
             "extended to close the body without them",
             left.empty() ? fate.vanish : left);
    faces.through(remove);
    return without;
}

TopoDS_Shape offsetBody(const TopoDS_Shape &work, double distance, double tolerance, int join,
                        const std::vector<TopoDS_Shape> &pushed, double pushedBy, Lineage &faces) {
    TopoDS_Shape result;
    try {
        if (pushed.empty()) {
            BRepOffsetAPI_MakeOffsetShape maker;
            // A drag or Esc can stop it in there: some bodies keep it a long time.
            maker.PerformByJoin(work, distance, tolerance, BRepOffset_Skin, false, false,
                                static_cast<GeomAbs_JoinType>(join), false, progress(true)->Start());
            checkCancelled("offset");
            if (!maker.IsDone())
                return {};
            result = maker.Shape();
            faces.images(maker);
        } else {
            BRepOffset_MakeOffset maker;
            // Faces with an offset of their own need joins by intersection.
            maker.Initialize(work, distance, tolerance, BRepOffset_Skin, false, false, GeomAbs_Intersection);
            for (const TopoDS_Shape &face : pushed)
                maker.SetOffsetOnFace(TopoDS::Face(face), pushedBy);
            maker.MakeOffsetShape(progress(true)->Start());
            checkCancelled("offset");
            if (!maker.IsDone())
                return {};
            result = maker.Shape();
            faces.images(maker);
        }
    } catch (const Standard_Failure &error) {
        checkCancelled("offset");
        const char *detail = error.what();
        // Only the new bounded normal retry is a named refusal. Other native errors
        // retain the null result: Hollow uses it to try its existing alternative.
        if (detail && std::string(detail).find("offset normal retry") != std::string::npos)
            fail("offset", no_solution, std::string("offset construction refused: ") + detail);
        return {};
    }
    result = singleSolid(result);
    if (result.IsNull() || !BRepCheck_Analyzer(result).IsValid())
        return {};
    // An offset past the middle of the body comes back inside out, or crossing itself.
    const double before = estimatedVolume(work), after = estimatedVolume(result);
    if (!(after > 0.0) || (pushed.empty() && (after - before) * distance <= 0.0) || intersectsItself(result))
        return {};
    faces.within(result);
    return result;
}

EvolutionData offsetSolid(OcctKernel &kernel, uint32_t shapeId, double distance, double tolerance, int join,
                          const std::vector<int> &hashes, int bound) {
    const char *const op = "offsetSolid";
    const TopoDS_Shape solid = singleSolid(get(kernel, shapeId));
    if (solid.IsNull() || !std::isfinite(distance) || distance == 0.0 || !std::isfinite(tolerance) ||
        tolerance <= 0.0 || join < 0 || join > 2)
        fail(op, invalid_input, "one solid, a distance that is not zero, a positive tolerance and a supported join are required");
    FaceHistory history(op, hashes, bound);
    Lineage faces(solid);
    const TopoDS_Shape work = withoutVanishing(op, dividedAtC0(solid, faces), distance, faces);
    const TopoDS_Shape result = offsetBody(work, distance, tolerance, join, {}, 0.0, faces);
    if (result.IsNull())
        fail(op, no_solution,
             distance < 0 ? "the body has no inside left at this distance, or its offset crosses itself"
                          : "the offset of the body crosses itself or could not be built");
    // Apart from the kernel's own word: nothing of the result is nearer to the input than asked.
    std::vector<TopoDS_Shape> skin;
    for (TopExp_Explorer all(result, TopAbs_FACE); all.More(); all.Next())
        skin.push_back(all.Current()); // Arc joins need not descend from an input face.
    const double nearest = clearance(skin, solid, {}, 0.0, !work.IsSame(solid));
    if (!std::isfinite(nearest))
        fail(op, no_solution, "the minimum clearance of the reconstructed offset could not be measured");
    if (nearest < std::abs(distance) - clearanceSlack * tolerance)
        fail(op, no_solution, "the offset comes nearer to the body than the distance; its reconstructed faces or joins do not preserve clearance");
    // As upstream's offsetWithHistory has it: an offset face is generated by the face it is of.
    faces.generated(history);
    for (TopExp_Explorer all(solid, TopAbs_FACE); all.More(); all.Next())
        if (history.tracks(all.Current()) && faces.of(all.Current()).empty())
            history.deleted(all.Current());
    return history.finish(store(kernel, result));
}

} // namespace BRepExt
