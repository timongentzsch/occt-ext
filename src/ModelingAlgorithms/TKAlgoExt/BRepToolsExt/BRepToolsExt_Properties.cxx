// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepToolsExt_Session.hxx>
#include <BRepToolsExt_Checks.hxx>
// Compile the authoritative numerical core in its original translation unit.
#include <BRepAlgoExt_Properties.cxx>

namespace loft {
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
std::vector<double> surfaceCenter(OcctKernel& kernel,uint32_t shapeId) {
    const auto shape=get(kernel,shapeId);
    double area=0;gp_XYZ moment(0,0,0);
    for(TopExp_Explorer it(shape,TopAbs_FACE);it.More();it.Next()) {
        const auto face=TopoDS::Face(it.Current());TrimProperties measured;
        if(trimmedSupport(face)) measured=trimmedFace(face);
        else {GProp_GProps props;BRepGProp::SurfaceProperties(face,props);measured.area=props.Mass();measured.center=props.CentreOfMass();}
        area+=measured.area;moment+=measured.center.XYZ()*measured.area;
    }
    if(!(area>0) || !std::isfinite(area)) {
        // With no area to weigh, the facade's centre as before.
        GProp_GProps props;BRepGProp::SurfaceProperties(shape,props);const auto center=props.CentreOfMass();
        return {center.X(),center.Y(),center.Z()};
    }
    const auto center=moment/area;
    return {center.X(),center.Y(),center.Z()};
}

} // namespace loft
