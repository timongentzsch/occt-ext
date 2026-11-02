// Copyright (c) 2026 Timon Gentzsch
//
// This file is part of Open CASCADE Technology software library.
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
// distribution for complete text of the license and disclaimer of any warranty.
//
#pragma once
#include <vector>

namespace GeomAlgoExt
{
//! The clamped, flat knot vector (each knot repeated by its multiplicity) of `spans` equal spans
//! of degree `degree` over [from, to].
inline std::vector<double> uniformKnots(int degree, int spans, double from = 0, double to = 1)
{
  std::vector<double> knots(degree + 1, from);
  for (int i = 1; i < spans; ++i)
    knots.push_back(from + (to - from) * i / spans);
  knots.insert(knots.end(), degree + 1, to);
  return knots;
}
} // namespace GeomAlgoExt
