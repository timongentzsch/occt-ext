// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <Approx_Curve3d.hxx>
#include <BRepAdaptor_CompCurve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <Geom_BSplineCurve.hxx>
#include <NCollection_Array1.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Wire.hxx>
#include <algorithm>
#include <cmath>

namespace loft {
namespace {
Handle(Geom_BSplineCurve) fixedFit(const Handle(Adaptor3d_Curve) & source, int count, int degree) {
    // Clamped uniform spline interpolating at Greville abscissae. Collocation is
    // a square system, solved for all three coordinates with pivoting; no
    // normal-equation squaring.
    const bool periodic =
        source->Value(source->FirstParameter()).Distance(source->Value(source->LastParameter())) <= loft::tolerance::exactLength;
    const double shift = degree % 2 == 0 ? 0.5 : 0.0;
    const int spans = periodic ? count : count - degree;
    std::vector<double> knots;
    for (int i = 0; i <= spans; ++i)
        for (int j = 0; j < (i == 0 || i == spans ? degree + 1 : 1); ++j)
            knots.push_back(double(i) / spans);
    std::vector<std::vector<double>> matrix(count, std::vector<double>(count + 3));
    Handle(Geom_BSplineCurve) periodicBasis;
    if (periodic) {
        NCollection_Array1<gp_Pnt> poles(1, count);
        NCollection_Array1<double> knots(1, count + 1);
        NCollection_Array1<int> mult(1, count + 1);
        for (int i = 1; i <= count; ++i)
            poles(i) = gp_Pnt(0, 0, 0);
        for (int i = 1; i <= count + 1; ++i) {
            knots(i) = double(i - 1) / count;
            mult(i) = 1;
        }
        periodicBasis = new Geom_BSplineCurve(poles, knots, mult, degree, true);
        for (int col = 0; col < count; ++col) {
            periodicBasis->SetPole(col + 1, gp_Pnt(1, 0, 0));
            for (int row = 0; row < count; ++row)
                matrix[row][col] = periodicBasis->Value((row + shift) / count).X();
            periodicBasis->SetPole(col + 1, gp_Pnt(0, 0, 0));
        }
    }
    for (int row = 0; row < count; ++row) {
        if (periodic) {
            const auto p =
                source->Value(source->FirstParameter() +
                              (row + shift) / count * (source->LastParameter() - source->FirstParameter()));
            matrix[row][count] = p.X();
            matrix[row][count + 1] = p.Y();
            matrix[row][count + 2] = p.Z();
            continue;
        }
        double u = 0;
        for (int j = 1; j <= degree; ++j)
            u += knots[row + j] / degree;
        std::vector<double> basis(knots.size() - 1);
        if (row == count - 1)
            basis[count - 1] = 1;
        else
            for (size_t i = 0; i + 1 < knots.size(); ++i)
                if (u >= knots[i] && u < knots[i + 1])
                    basis[i] = 1;
        for (int d = 1; d <= degree; ++d) {
            std::vector<double> next(basis.size());
            for (int i = 0; i < count; ++i) {
                double a = knots[i + d] - knots[i], b = knots[i + d + 1] - knots[i + 1];
                if (a > 0)
                    next[i] += (u - knots[i]) / a * basis[i];
                if (b > 0)
                    next[i] += (knots[i + d + 1] - u) / b * basis[i + 1];
            }
            basis = std::move(next);
        }
        if (row == count - 1) {
            std::fill(basis.begin(), basis.end(), 0);
            basis[count - 1] = 1;
        }
        for (int i = 0; i < count; ++i)
            matrix[row][i] = basis[i];
        const auto p = source->Value(source->FirstParameter() +
                                     u * (source->LastParameter() - source->FirstParameter()));
        matrix[row][count] = p.X();
        matrix[row][count + 1] = p.Y();
        matrix[row][count + 2] = p.Z();
    }
    if (!solveForThree(matrix, count))
        fail("rebuildCurve", no_solution, "spline collocation is singular");
    NCollection_Array1<gp_Pnt> poles(1, count);
    NCollection_Array1<double> distinct(1, spans + 1);
    NCollection_Array1<int> mult(1, spans + 1);
    for (int i = 0; i < count; ++i)
        poles(i + 1) = gp_Pnt(matrix[i][count], matrix[i][count + 1], matrix[i][count + 2]);
    for (int i = 0; i <= spans; ++i) {
        distinct(i + 1) = double(i) / spans;
        mult(i + 1) = !periodic && (i == 0 || i == spans) ? degree + 1 : 1;
    }
    return new Geom_BSplineCurve(poles, distinct, mult, degree, periodic);
}
} // namespace
EvolutionData rebuildCurve(OcctKernel &kernel, uint32_t wireId, int mode, double tolerance, bool keepCorners,
                           int degree, int count) {
    const char *op = "rebuildCurve";
    const auto &shape = get(kernel, wireId);
    if (shape.ShapeType() != TopAbs_WIRE || mode < 0 || mode > 1 || !std::isfinite(tolerance) ||
        tolerance <= 0 || degree < 1 || degree > 9 || count <= degree || count > 256)
        fail(op, invalid_input,
             "rebuild requires a wire, positive tolerance, degree 1-9 and degree+1 "
             "to 256 points");
    std::vector<TopoDS_Wire> pieces;
    if (mode == 0 && keepCorners)
        for (BRepTools_WireExplorer it(TopoDS::Wire(shape)); it.More(); it.Next())
            pieces.push_back(BRepBuilderAPI_MakeWire(it.Current()).Wire());
    else
        pieces.push_back(TopoDS::Wire(shape));
    BRepBuilderAPI_MakeWire output;
    for (const auto &wire : pieces) {
        Handle(Adaptor3d_Curve) curve = new BRepAdaptor_CompCurve(wire, true);
        Handle(Geom_BSplineCurve) spline;
        if (mode == 0) {
            Approx_Curve3d fit(curve, tolerance, GeomAbs_C1, 256, degree);
            if (!fit.IsDone() || fit.MaxError() > tolerance)
                fail(op, no_solution, "refit could not meet the requested tolerance");
            spline = fit.Curve();
        } else
            spline = fixedFit(curve, count, degree);
        output.Add(BRepBuilderAPI_MakeEdge(spline).Edge());
    }
    if (!output.IsDone() || !BRepCheck_Analyzer(output.Wire()).IsValid())
        fail(op, no_solution, "rebuild produced an invalid wire");
    return stored(kernel, output.Wire());
}
} // namespace loft
