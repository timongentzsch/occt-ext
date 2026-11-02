// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
#include <cmath>

namespace BRepExt
{
/// Whether every coordinate of a point or vector is a number.
template <typename Coordinates>
bool finite(const Coordinates& p)
{
  return std::isfinite(p.X()) && std::isfinite(p.Y()) && std::isfinite(p.Z());
}

} // namespace BRepExt
