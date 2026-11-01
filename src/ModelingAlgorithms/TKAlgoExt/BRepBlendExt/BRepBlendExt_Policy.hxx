// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
#include <array>
#include <utility>
#include <cmath>
#include <GeomAlgoExt_Failure.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>
namespace BRepExt::blend2 {
struct NativePolicy {
    inline static constexpr double lawShare = 0;
    static bool tangentFor(double angle, double size) { return angle <= tolerance::coarseAngle || angle * size <= tolerance::length; }
    static double conicWeight(double rho, double wedge) { return std::sin(wedge / 2) * rho / (1 - rho); }
    static std::array<std::pair<double,double>,8> g2Poles(double) { fail("blend", unsupported_geometry, "supply a curvature-continuous section policy"); }
    inline static constexpr std::array<double,3> g2Knots = {};
    inline static constexpr std::array<int,3> g2Multiplicities = {};
    template<class Station, class Point, class Vector>
    static bool curvedG2(double, const Station&, const Point&, const Vector*, const double*, std::vector<Point>&, int&) { return false; }
};
}
