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

#include <BRepToolsExt_Operations.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <Geom_BSplineCurve.hxx>
#include <TopoDS_Edge.hxx>
#include <NCollection_Array1.hxx>
#include <math_Gauss.hxx>
#include <math_Matrix.hxx>
#include <math_Vector.hxx>
#include <cmath>

namespace BRepExt
{
EvolutionData interpolateClosed(OcctKernel& kernel, const std::vector<double>& coordinates)
{
  const char* op = "interpolateClosed";
  if (coordinates.size() < 9 || coordinates.size() % 3 != 0)
    fail(op, invalid_input, "at least three points are required");
  for (double v : coordinates)
    if (!std::isfinite(v))
      fail(op, invalid_input, "points must be finite");
  const int           n = static_cast<int>(coordinates.size() / 3);
  std::vector<gp_Pnt> points;
  for (int i = 0; i < n; ++i)
    points.emplace_back(coordinates[3 * i], coordinates[3 * i + 1], coordinates[3 * i + 2]);
  NCollection_Array1<double> knots(1, n + 1);
  NCollection_Array1<int>    mults(1, n + 1);
  knots(1) = 0;
  mults.Init(1);
  for (int i = 0; i < n; ++i)
  {
    const double span = points[i].Distance(points[(i + 1) % n]);
    if (!std::isfinite(span) || span < 1e-12)
      fail(op, invalid_input, "successive points must be distinct");
    knots(i + 2) = knots(i + 1) + span;
    if (!std::isfinite(knots(i + 2)) || knots(i + 2) <= knots(i + 1))
      fail(op, invalid_input, "chord-length knots are not representable");
  }
  const double period = knots(n + 1);
  auto         knot   = [&](int i) {
    const int turn = static_cast<int>(std::floor(static_cast<double>(i) / n));
    return knots(i - turn * n + 1) + turn * period;
  };
  math_Matrix matrix(1, n, 1, n, 0.0);
  for (int i = 0; i < n; ++i)
  {
    const double x     = knot(i);
    auto         basis = [&](auto&& self, int start, int degree) -> double {
      if (!degree)
        return x >= knot(start) && x < knot(start + 1) ? 1.0 : 0.0;
      return (x - knot(start)) / (knot(start + degree) - knot(start))
               * self(self, start, degree - 1)
             + (knot(start + degree + 1) - x) / (knot(start + degree + 1) - knot(start + 1))
                 * self(self, start + 1, degree - 1);
    };
    for (int j = 0; j < 3; ++j)
      matrix(i + 1, (i - 3 + j + n) % n + 1) += basis(basis, i - 3 + j, 3);
  }
  math_Gauss solve(matrix, 1e-14);
  if (!solve.IsDone())
    fail(op, no_solution, "periodic interpolation system is singular");
  NCollection_Array1<gp_Pnt> poles(1, n);
  for (int axis = 1; axis <= 3; ++axis)
  {
    math_Vector right(1, n), solution(1, n);
    for (int i = 0; i < n; ++i)
      right(i + 1) = points[i].Coord(axis);
    solve.Solve(right, solution);
    // OCCT's first pole is the basis on the three knots preceding the origin.
    // The cyclic basis index is therefore n-3, independent of the input points.
    for (int i = 0; i < n; ++i)
    {
      const double value = solution((i + n - 3) % n + 1);
      if (!std::isfinite(value))
        fail(op, no_solution, "interpolation poles are not finite");
      poles(i + 1).SetCoord(axis, value);
    }
  }
  Handle(Geom_BSplineCurve) curve = new Geom_BSplineCurve(poles, knots, mults, 3, true);
  for (int i = 0; i < n; ++i)
    if (curve->Value(knots(i + 1)).Distance(points[i]) > 1e-10 * (1 + period))
      fail(op, no_solution, "periodic interpolation residual exceeds tolerance");
  return stored(kernel, BRepBuilderAPI_MakeEdge(curve).Edge());
}
} // namespace BRepExt
