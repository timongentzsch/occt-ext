// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepToolsExt_Checks.hxx>

#include <BOPAlgo_Builder.hxx>
#include <BOPAlgo_BuilderFace.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_FindPlane.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepLib.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>
#include <BRep_Builder.hxx>
#include <Geom_Plane.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <gp_Pln.hxx>

#include <cmath>

namespace BRepExt {

EvolutionData regions(OcctKernel& kernel, const std::vector<uint32_t>& curveIds,
                      Vec3 origin, Vec3 normal, double tolerance) {
    const char* const op = "regions";
    if (curveIds.empty() || !std::isfinite(tolerance) || tolerance <= 0.0 || !finite(origin) || !finite(normal) ||
        gp_Vec(normal.x, normal.y, normal.z).Magnitude() <= gp::Resolution()) {
        fail(op, invalid_input, "curves, a finite plane and a positive tolerance are required");
    }
    BRep_Builder builder;
    TopoDS_Compound input;
    builder.MakeCompound(input);
    for (uint32_t id : curveIds) {
        const auto& curve = get(kernel, id);
        if (curve.ShapeType() != TopAbs_WIRE && curve.ShapeType() != TopAbs_EDGE) {
            fail(op, invalid_input, "regions accepts only edges and wires");
        }
        builder.Add(input, curve);
    }
    // P-curves and intersection tolerances belong to this operation, never to the inputs.
    const TopoDS_Shape copy = BRepBuilderAPI_Copy(input, true, false).Shape();
    BRepBuilderAPI_FindPlane found(copy, tolerance);
    const gp_Pln plane(gp_Pnt(origin.x, origin.y, origin.z), gp_Dir(normal.x, normal.y, normal.z));
    if (!found.Found() || !found.Plane()->Pln().Axis().Direction().IsParallel(plane.Axis().Direction(), BRepExt::tolerance::angle) ||
        plane.Distance(found.Plane()->Location()) > tolerance) {
        fail(op, invalid_input, "the curves must lie in the specified plane");
    }
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> edges;
    TopExp::MapShapes(copy, TopAbs_EDGE, edges);
    TopoDS_Shape split = copy;
    if (edges.Extent() > 1) {
        BOPAlgo_Builder intersections;
        for (const auto& edge : edges) intersections.AddArgument(edge);
        intersections.SetNonDestructive(true);
        intersections.SetFuzzyValue(tolerance);
        intersections.Perform();
        if (intersections.HasErrors()) fail(op, no_solution, "curve intersections failed");
        split = intersections.Shape();
    }
    edges.Clear();
    TopExp::MapShapes(split, TopAbs_EDGE, edges);
    const TopoDS_Face support = BRepBuilderAPI_MakeFace(plane).Face();
    NCollection_List<TopoDS_Shape> directed;
    for (const auto& edge : edges) {
        BRepLib::BuildPCurveForEdgeOnPlane(TopoDS::Edge(edge), support);
        directed.Append(edge.Oriented(TopAbs_FORWARD));
        directed.Append(edge.Oriented(TopAbs_REVERSED));
    }
    BOPAlgo_BuilderFace areas;
    areas.SetFace(support);
    areas.SetShapes(directed);
    areas.Perform();
    if (areas.HasErrors()) fail(op, no_solution, "closed regions could not be assembled");
    TopoDS_Compound result;
    builder.MakeCompound(result);
    int count = 0;
    for (const auto& area : areas.Areas()) {
        const auto& face = TopoDS::Face(area);
        if (BRepTopAdaptor_FClass2d(face, tolerance).PerformInfinitePoint() != TopAbs_OUT) continue;
        if (!BRepCheck_Analyzer(face).IsValid()) fail(op, no_solution, "a region is not a valid face");
        builder.Add(result, face);
        count++;
    }
    if (count == 0) fail(op, no_solution, "the curves bound no closed region");
    return stored(kernel, result);
}

} // namespace BRepExt
