// Copyright (c) 2026 Timon Gentzsch
// This file is part of Open CASCADE Technology software library.
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
// distribution for complete text of the license and disclaimer of any warranty.
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
