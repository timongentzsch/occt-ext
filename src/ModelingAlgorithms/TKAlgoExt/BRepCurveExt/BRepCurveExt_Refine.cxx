// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepToolsExt_Operations.hxx>
#include <BRepToolsExt_Checks.hxx>

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Tool.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <GeomConvert.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Wire.hxx>

#include <cmath>

namespace BRepExt {

EvolutionData refineCurve(OcctKernel& kernel, uint32_t wireId, int mode, int amount, double tolerance) {
    const char* const op = "refineCurve";
    TopoDS_Shape shape = get(kernel, wireId);
    if (shape.ShapeType() != TopAbs_WIRE || mode < 0 || mode > 2 || amount < 1 || amount > 24 ||
        !std::isfinite(tolerance) || tolerance <= 0)
        fail(op, invalid_input, "a wire, refinement mode, positive amount and tolerance are required");
    if (mode == 2) {
        ShapeUpgrade_UnifySameDomain unify(shape, true, false, false);
        unify.SetLinearTolerance(tolerance);
        unify.Build();
        shape = unify.Shape();
    }
    BRepBuilderAPI_MakeWire wire;
    for (BRepTools_WireExplorer it(TopoDS::Wire(shape)); it.More(); it.Next()) {
        const TopoDS_Edge edge = it.Current();
        double first, last;
        const auto curve = BRep_Tool::Curve(edge, first, last);
        if (curve.IsNull()) fail(op, unsupported_geometry, "an edge has no 3D curve");
        if (mode == 2 && Handle(Geom_BSplineCurve)::DownCast(curve).IsNull()) {
            wire.Add(edge);
            continue;
        }
        // Convert just the trimmed, located geometry. Upstream elevation exposes the whole
        // underlying spline, which lengthens split edges and loses their reversed orientation.
        auto spline = GeomConvert::CurveToBSplineCurve(new Geom_TrimmedCurve(curve, first, last));
        if (mode == 0) {
            if (spline->Degree() + amount > Geom_BSplineCurve::MaxDegree())
                fail(op, invalid_input, "the resulting degree exceeds OCCT's maximum (25)");
            spline->IncreaseDegree(spline->Degree() + amount);
        } else if (mode == 1) {
            const int target = 2 * spline->NbPoles();
            if (target > 16384) fail(op, invalid_input, "subdivision would exceed 16384 control points");
            while (spline->NbPoles() < target) {
                int widest = 1;
                for (int i = 2; i < spline->NbKnots(); ++i)
                    if (spline->Knot(i + 1) - spline->Knot(i) > spline->Knot(widest + 1) - spline->Knot(widest)) widest = i;
                const int before = spline->NbPoles();
                spline->InsertKnot((spline->Knot(widest) + spline->Knot(widest + 1)) / 2, 1, 0.0);
                if (spline->NbPoles() == before) fail(op, no_solution, "knot spacing is below numerical precision");
            }
        } else {
            const double budget = tolerance / spline->NbPoles();
            for (int i = spline->NbKnots() - 1; i > 1; --i) {
                for (int mult = spline->Multiplicity(i) - 1; mult >= 0; --mult)
                    if (!spline->RemoveKnot(i, mult, budget)) break;
            }
        }
        TopoDS_Edge result = BRepBuilderAPI_MakeEdge(spline);
        result.Orientation(edge.Orientation());
        wire.Add(result);
    }
    if (!wire.IsDone() || !BRepCheck_Analyzer(wire.Wire()).IsValid())
        fail(op, no_solution, "refinement did not produce a valid connected wire");
    return stored(kernel, wire.Wire());
}

} // namespace BRepExt
