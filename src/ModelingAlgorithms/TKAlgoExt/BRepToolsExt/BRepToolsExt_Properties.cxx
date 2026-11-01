// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepToolsExt_Session.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepAlgoExt_Properties.hxx>
#include <cmath>

namespace BRepExt {
double volume(OcctKernel &kernel, uint32_t shapeId) { return measure(get(kernel,shapeId), false); }
std::vector<double> volumeResult(OcctKernel& kernel,uint32_t shapeId) {
    bool converged;const double value=measure(get(kernel,shapeId),false,&converged);
    return {value,converged?1.0:0.0};
}
double volumeEstimate(OcctKernel &kernel, uint32_t shapeId) {
    const double value=estimatedVolume(get(kernel, shapeId));
    if (!std::isfinite(value))
        fail("volumeEstimate", no_solution, "property integration failed");
    return value;
}
double surfaceAreaEstimate(OcctKernel &kernel, uint32_t shapeId) {
    const double value=estimatedArea(get(kernel, shapeId));
    if (!std::isfinite(value))
        fail("surfaceAreaEstimate", no_solution, "property integration failed");
    return value;
}
double surfaceArea(OcctKernel &kernel, uint32_t shapeId) { return measure(get(kernel,shapeId),true); }
std::vector<double> surfaceAreaResult(OcctKernel& kernel,uint32_t shapeId) {
    bool converged;const double value=measure(get(kernel,shapeId),true,&converged);
    return {value,converged?1.0:0.0};
}
std::vector<double> surfaceCenter(OcctKernel& kernel,uint32_t shapeId) { return surfaceCentre(get(kernel,shapeId)); }

} // namespace BRepExt
