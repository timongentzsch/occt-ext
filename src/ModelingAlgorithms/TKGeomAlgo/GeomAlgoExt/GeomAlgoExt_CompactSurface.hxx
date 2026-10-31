// Copyright (c) 2026 Timon Gentzsch
//
// This file is part of Open CASCADE Technology software library.
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult LICENSE_LGPL_21.txt for the complete license.

#pragma once
#include <Geom_BSplineSurface.hxx>
#include <GeomAlgoExt_Execution.hxx>
#include <optional>

namespace GeomAlgoExt {
struct CompactPatchSurface {
    occ::handle<Geom_BSplineSurface> surface;
    double displacement;
    int removals;
};

/// Optional private reduction of a nonrational, nonperiodic clamped surface. Same-parameter
/// displacement is bounded over its entire rectangle by outward-rounded interval differences of coefficients
/// in a common Bezier basis (including both knot conversions). No curve, trim, normal or continuity acceptance is implied:
/// the caller must rebuild the original rims and repeat all original physical gates.
/// No certified removal returns nullopt. A work limit may return a certified partial reduction.
Standard_EXPORT std::optional<CompactPatchSurface> compactSurface(
    const Message_ProgressRange &operation, const occ::handle<Geom_BSplineSurface> &original,
    double maximumDisplacement);
}
