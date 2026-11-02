// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
#include <Geom_BSplineSurface.hxx>
#include <optional>

namespace BRepExt
{
struct CompactPatchSurface
{
  occ::handle<Geom_BSplineSurface> surface;
  double                           displacement;
  int                              removals;
};

/// Optional private reduction of a nonrational, nonperiodic clamped surface. Same-parameter
/// displacement is bounded over its entire rectangle by outward-rounded interval differences of
/// coefficients in a common Bezier basis (including both knot conversions). No curve, trim, normal
/// or continuity acceptance is implied: the caller must rebuild the original rims and repeat all
/// original physical gates. No certified removal returns nullopt. A work limit may return a
/// certified partial reduction.
std::optional<CompactPatchSurface> compactPatchSurface(
  const char*                             operation,
  const occ::handle<Geom_BSplineSurface>& original,
  double                                  maximumDisplacement);
} // namespace BRepExt
