// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
#include <Standard_Macro.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <vector>

namespace BRepExt
{
// Adaptive measurement. Each face has its own integration budget; callers needing
// an accuracy certificate must inspect `converged` as well as the returned value.
Standard_EXPORT double faceMeasure(const TopoDS_Face&,
                                   bool          surface,
                                   const gp_Pnt& origin,
                                   bool*         converged,
                                   const double* areaEstimate = nullptr);
// Closed forms for recognised analytic rectangles, otherwise fixed OCCT quadrature.
Standard_EXPORT double faceArea(const TopoDS_Face&);
// Signed volume (surface false) or area of a whole shape by the adaptive face measurement above.
// With `converged`, reports whether every face met its budget.
Standard_EXPORT double measure(const TopoDS_Shape&, bool surface, bool* converged = nullptr);
// The area centroid, integrated over each face's trim where the support is not an elementary
// surface.
Standard_EXPORT std::vector<double> surfaceCentre(const TopoDS_Shape&);
// Fixed boundary quadrature used only to bracket Boolean volume error.
Standard_EXPORT double faceMeasureBracket(const TopoDS_Face&, const gp_Pnt& origin);

// Fixed OCCT quadrature, with no error/convergence certificate. These estimates
// retain the existing heuristic acceptance policy; use faceMeasure/volumeResult
// when a numerical accuracy contract is required. Keep this distinct from the
// adaptive public measurement rather than silently changing solver decisions.
inline double estimatedVolume(const TopoDS_Shape& shape)
{
  GProp_GProps properties;
  BRepGProp::VolumeProperties(shape, properties);
  return properties.Mass();
}

inline double estimatedArea(const TopoDS_Shape& shape)
{
  GProp_GProps properties;
  BRepGProp::SurfaceProperties(shape, properties);
  return properties.Mass();
}
} // namespace BRepExt
