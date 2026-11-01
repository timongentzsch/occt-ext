// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>

#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepOffsetAPI_MakeFilling.hxx>
#include <BRep_Tool.hxx>
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Geom_BSplineSurface.hxx>
#include <NCollection_Array2.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>

#include <cmath>

namespace loft {

EvolutionData surfaceThroughGrid(OcctKernel &kernel, const std::vector<double> &coordinates, int rows,
                                 int columns, double tolerance) {
    const char *const op = "surfaceThroughGrid";
    if (rows < 2 || columns < 2 || rows > 1024 || columns > 1024 ||
        coordinates.size() != static_cast<size_t>(rows) * columns * 3 || !std::isfinite(tolerance) ||
        tolerance <= 0)
        fail(op, invalid_input, "a rectangular grid and a positive finite tolerance are required");
    NCollection_Array2<gp_Pnt> points(1, rows, 1, columns);
    for (int r = 0; r < rows; ++r)
        for (int c = 0; c < columns; ++c) {
            const size_t i = (static_cast<size_t>(r) * columns + c) * 3;
            if (!std::isfinite(coordinates[i]) || !std::isfinite(coordinates[i + 1]) ||
                !std::isfinite(coordinates[i + 2]))
                fail(op, invalid_input, "grid coordinates must be finite");
            points.SetValue(r + 1, c + 1, gp_Pnt(coordinates[i], coordinates[i + 1], coordinates[i + 2]));
        }
    // Uniform interpolation exposes the exact parameter of each grid point.
    // Checking values is linear in the grid size and does not build a global
    // Extrema per point.
    GeomAPI_PointsToBSplineSurface fit;
    fit.Interpolate(points, Approx_IsoParametric);
    if (!fit.IsDone())
        fail(op, no_solution, "grid interpolation failed");
    double u0, u1, v0, v1;
    fit.Surface()->Bounds(u0, u1, v0, v1);
    for (int r = 0; r < rows; ++r)
        for (int c = 0; c < columns; ++c)
            if (fit.Surface()
                    ->Value(u0 + (u1 - u0) * r / (rows - 1), v0 + (v1 - v0) * c / (columns - 1))
                    .Distance(points(r + 1, c + 1)) > tolerance)
                fail(op, no_solution, "the surface does not meet the grid tolerance");
    BRepBuilderAPI_MakeFace face(fit.Surface(), tolerance);
    if (!face.IsDone() || !BRepCheck_Analyzer(face.Face()).IsValid())
        fail(op, no_solution, "surface fitting produced an invalid face");
    return stored(kernel, face.Face());
}

EvolutionData fillBoundary(OcctKernel &kernel, uint32_t wireId, double tolerance) {
    const char *const op = "fillBoundary";
    const auto wire = BRepBuilderAPI_Copy(get(kernel, wireId), true, false).Shape();
    if (wire.ShapeType() != TopAbs_WIRE || !BRep_Tool::IsClosed(wire) || !std::isfinite(tolerance) ||
        tolerance <= 0)
        fail(op, invalid_input, "a closed boundary wire and positive finite tolerance are required");
    BRepOffsetAPI_MakeFilling fill(3, 30, 3, false, tolerance, tolerance, 0.01, 0.1, 8, 32);
    for (TopExp_Explorer it(wire, TopAbs_EDGE); it.More(); it.Next())
        fill.Add(TopoDS::Edge(it.Current()), GeomAbs_C0);
    fill.Build();
    if (!fill.IsDone() || !BRepCheck_Analyzer(fill.Shape()).IsValid())
        fail(op, no_solution, "boundary filling produced an invalid face");
    if (!std::isfinite(fill.G0Error()) || fill.G0Error() > tolerance)
        fail(op, no_solution, "the patch does not meet the boundary tolerance");
    return stored(kernel, fill.Shape());
}

} // namespace loft
