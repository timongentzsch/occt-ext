// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once

#include <TopoDS_Shape.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <GeomAlgoExt_Resolution.hxx>

#include <algorithm>
#include <limits>

namespace loft::tolerance {

// Lengths in model units. Preserve the existing thresholds; callers' document tolerances
// still govern requested fits. Do not substitute these for solver steps or relative errors.
inline constexpr double length = 1e-7;       // OCCT's nominal topological resolution
inline constexpr double sewing = 1e-6;      // default joining / point-coincidence allowance
inline constexpr double exactLength = 1e-9; // exact curve endpoint / analytic coincidence
// Angles in radians (or the existing normalized dot/cross comparison).
inline constexpr double coarseAngle = 1e-6; // blend tangency / sweep angular allowance
inline constexpr double angularResolution = 1e-7; // full-period / directional recognition
inline constexpr double angle = 1e-9;
// Parametric coordinates, scaled by a span where the caller already does so.
inline constexpr double parameter = 1e-12;
inline constexpr double uv = 1e-9;           // UV projection / face classification resolution

// How finely a length that is the difference of two points computed at model coordinates can be
// resolved: each coordinate is good to some units in the last place (ulps) of its magnitude, not
// of the length. A point on a surface takes a few tens of roundings (de Boor's recurrence to
// degree 9, a rational weight, an offset along the normal) and a residual is a difference of two
// such points; 64 ulps bound that with room. Below 7e6 from the origin it stays under `length`.
inline constexpr double roundingUlps = 64;
inline double rounding(double magnitude) { return roundingUlps * std::numeric_limits<double>::epsilon() * magnitude; }
// A tolerance on such a length that is `share` of the size of what is built (a blend's radius),
// but never finer than the rounding of points within `magnitude` of the origin: the same
// relative accuracy at any scale, and one that can be reached wherever the model sits.
inline double resolved(double share, double size, double magnitude) {
    return std::max(share * size, rounding(magnitude + size));
}

// The cell Boolean's fuzzy ceiling and the tolerance its swept crossings are found within.
// Scaling them by model size refused Booleans these constants make, on small models and
// at a 1e-6 relative crossing; they stay absolute
// until a corpus covers scales.
inline constexpr double booleanCeiling = 1e-5;
inline constexpr double booleanCrossing = 1e-4;
inline constexpr double booleanSphereRelative = 1e-4;
inline constexpr double booleanSurfaceParameter = 1e-10;
inline constexpr double booleanQuadrature = 5e-8;
inline constexpr double booleanFuzzyRelative = 1e-5;

} // namespace loft::tolerance
