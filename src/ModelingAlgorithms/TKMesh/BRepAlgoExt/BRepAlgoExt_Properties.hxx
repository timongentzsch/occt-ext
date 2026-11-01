// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <TopoDS_Face.hxx>

namespace loft {
// Adaptive measurement. Each face has its own integration budget; callers needing
// an accuracy certificate must inspect `converged` as well as the returned value.
double faceMeasure(const TopoDS_Face&, bool surface, const gp_Pnt& origin, bool* converged, const double* areaEstimate = nullptr);
// Closed forms for recognised analytic rectangles, otherwise fixed OCCT quadrature.
double faceArea(const TopoDS_Face&);
// Fixed boundary quadrature used only to bracket Boolean volume error.
double faceMeasureBracket(const TopoDS_Face&, const gp_Pnt& origin);
// Fixed OCCT quadrature, with no error/convergence certificate. These estimates
// retain the existing heuristic acceptance policy; use faceMeasure/volumeResult
// when a numerical accuracy contract is required. Keep this distinct from the
// adaptive public measurement rather than silently changing solver decisions.
inline double estimatedVolume(const TopoDS_Shape& shape) {
    GProp_GProps properties;
    BRepGProp::VolumeProperties(shape, properties);
    return properties.Mass();
}
inline double estimatedArea(const TopoDS_Shape& shape) {
    GProp_GProps properties;
    BRepGProp::SurfaceProperties(shape, properties);
    return properties.Mass();
}
} // namespace loft
