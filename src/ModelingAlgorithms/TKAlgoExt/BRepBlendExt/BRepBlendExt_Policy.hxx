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
// Alternatively, this file may be used under the terms of Open CASCADE
// commercial license or contractual agreement.

#pragma once
#include <array>
#include <utility>
#include <cmath>
#include <GeomAlgoExt_Failure.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>

namespace BRepExt::blending
{
struct NativePolicy
{
  inline static constexpr double lawShare = 0;

  static bool tangentFor(double angle, double size)
  {
    return angle <= tolerance::coarseAngle || angle * size <= tolerance::length;
  }

  static double conicWeight(double rho, double wedge)
  {
    return std::sin(wedge / 2) * rho / (1 - rho);
  }

  static std::array<std::pair<double, double>, 8> g2Poles(double)
  {
    fail("blend", unsupported_geometry, "supply a curvature-continuous section policy");
  }

  inline static constexpr std::array<double, 3> g2Knots          = {};
  inline static constexpr std::array<int, 3>    g2Multiplicities = {};

  template <class Station, class Point, class Vector>
  static bool curvedG2(double,
                       const Station&,
                       const Point&,
                       const Vector*,
                       const double*,
                       std::vector<Point>&,
                       int&)
  {
    return false;
  }
};
} // namespace BRepExt::blending
