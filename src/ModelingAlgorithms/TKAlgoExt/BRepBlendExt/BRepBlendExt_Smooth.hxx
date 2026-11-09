// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#pragma once
#include <algorithm>
#include <cmath>

namespace BRepExt::blending
{
// Compact C2 upper rounding of max(a,b). Derived by matching |x| and its
// first two derivatives at x=+-1. Its excess is at most 3 width / 16.
inline double smoothMaximum(double a, double b, double width)
{
  if (!(width > 0) || std::abs(a - b) >= width)
    return std::max(a, b);
  const double x = (a - b) / width;
  const double rounded = width * (3.0 / 8 + x * x * (3.0 / 4 - x * x / 8));
  return (a + b + rounded) / 2;
}

inline double smoothMinimum(double a, double b, double width)
{
  return -smoothMaximum(-a, -b, width);
}

// Quintic Hermite activation (Perlin, GPU Gems, 2004, section 5.4).
inline double smoothActivation(double x)
{
  if (x <= 0)
    return 0;
  if (x >= 1)
    return 1;
  return x * x * x * (10 + x * (-15 + 6 * x));
}
} // namespace BRepExt::blending
