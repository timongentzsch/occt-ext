// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepToolsExt_Operations.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepToolsExt_Checks.hxx>

#include <BRepAdaptor_Curve.hxx>
#include <TopoDS.hxx>

#include <cmath>

namespace loft {

std::vector<double> curvature(OcctKernel& kernel, uint32_t edgeId, const std::vector<double>& parameters) {
    const char* const op = "curvature";
    const TopoDS_Shape& shape = get(kernel, edgeId);
    if (shape.ShapeType() != TopAbs_EDGE) fail(op, invalid_input, "curvature requires an edge");
    BRepAdaptor_Curve curve(TopoDS::Edge(shape));
    const bool reversed = shape.Orientation() == TopAbs_REVERSED;
    const double first = curve.FirstParameter(), last = curve.LastParameter();
    std::vector<double> result;
    result.reserve(parameters.size() * 9);
    for (double fraction : parameters) {
        if (!std::isfinite(fraction) || fraction < 0.0 || fraction > 1.0) {
            fail(op, invalid_input, "sample parameters must be between zero and one");
        }
        const double u = first + (reversed ? 1.0 - fraction : fraction) * (last - first);
        gp_Pnt point;
        gp_Vec d1, d2;
        curve.D2(u, point, d1, d2);
        const double speed = d1.Magnitude();
        if (!std::isfinite(speed) || speed == 0.0) fail(op, no_solution, "curvature is undefined at a singular point");
        gp_Vec tangent = d1 / speed;
        const gp_Vec bending = ((d2 - tangent * d2.Dot(tangent)) / speed) / speed;
        if (!finite(bending)) {
            fail(op, no_solution, "curvature is not finite");
        }
        if (reversed) tangent.Reverse();
        result.insert(result.end(), {point.X(), point.Y(), point.Z(), tangent.X(), tangent.Y(), tangent.Z(),
                                     bending.X(), bending.Y(), bending.Z()});
    }
    return result;
}

} // namespace loft
