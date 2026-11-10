// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#pragma once

#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <Geom2d_Curve.hxx>
#include <optional>

namespace BRepExt
{
//! An analytic image of an affine isoparametric trim. The spatial range is
//! independent of the pcurve parameter. Error bounds the two-sided set distance
//! to the lifted pcurve, including outward rounding of the image endpoints.
struct TrimImage
{
  Handle(Geom_Curve) curve;
  double             first, last, error;
};

//! Elementary carriers only. An absent image leaves the general trim builder
//! responsible; it does not change the trimmed domain or its clearance test.
Standard_EXPORT std::optional<TrimImage> isoTrimImage(const Handle(Geom_Surface)& surface,
                                                      const Handle(Geom2d_Curve)& pcurve,
                                                      double                      first,
                                                      double                      last,
                                                      double                      tolerance);
} // namespace BRepExt
