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

#include <BRepAdaptor_Curve2d.hxx>
#include <BRepTools.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_BezierCurve.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <algorithm>
#include <cmath>
#include <limits>

namespace GeomAlgoExt
{
// A conservative shortcut contract shared by properties and bounds. One wire,
// covering all four sides of its UV box. Lines are exact; positive-weight spline
// pole hulls certify constant coordinates without sampling between knots.
// An unrecognised boundary takes the general trimmed-face path.
inline bool parameterRectangle(const TopoDS_Face& face, double u0, double u1, double v0, double v1)
{
  int wires = 0;
  for (TopExp_Explorer it(face, TopAbs_WIRE); it.More(); it.Next())
    ++wires;
  if (wires != 1 || !(u1 > u0) || !(v1 > v0))
    return false;
  const auto allowance = [](double a, double b) {
    return 64 * std::numeric_limits<double>::epsilon() * std::max(std::abs(a), std::abs(b))
           + 1e-12 * (b - a);
  };
  const double tu = allowance(u0, u1), tv = allowance(v0, v1);
  unsigned     sides = 0;
  for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next())
  {
    BRepAdaptor_Curve2d curve(TopoDS::Edge(it.Current()), face);
    double              a = std::numeric_limits<double>::infinity(), b = a, c = -a, d = -a;
    const auto          add = [&](const gp_Pnt2d& p) {
      a = std::min(a, p.X());
      b = std::min(b, p.Y());
      c = std::max(c, p.X());
      d = std::max(d, p.Y());
    };
    if (curve.GetType() == GeomAbs_Line)
    {
      add(curve.Value(curve.FirstParameter()));
      add(curve.Value(curve.LastParameter()));
    }
    else if (curve.GetType() == GeomAbs_BSplineCurve)
    {
      const auto spline = curve.BSpline();
      for (int i = 1; i <= spline->NbPoles(); ++i)
      {
        if (!(spline->Weight(i) > 0))
          return false;
        add(spline->Pole(i));
      }
    }
    else if (curve.GetType() == GeomAbs_BezierCurve)
    {
      const auto spline = curve.Bezier();
      for (int i = 1; i <= spline->NbPoles(); ++i)
      {
        if (!(spline->Weight(i) > 0))
          return false;
        add(spline->Pole(i));
      }
    }
    else
      return false;
    unsigned side = 0;
    if (std::abs(a - u0) <= tu && std::abs(c - u0) <= tu)
      side |= 1;
    if (std::abs(a - u1) <= tu && std::abs(c - u1) <= tu)
      side |= 2;
    if (std::abs(b - v0) <= tv && std::abs(d - v0) <= tv)
      side |= 4;
    if (std::abs(b - v1) <= tv && std::abs(d - v1) <= tv)
      side |= 8;
    if (!side)
      return false;
    sides |= side;
  }
  return sides == 15;
}

inline bool parameterRectangle(const TopoDS_Face& face)
{
  double u0, u1, v0, v1;
  BRepTools::UVBounds(face, u0, u1, v0, v1);
  return parameterRectangle(face, u0, u1, v0, v1);
}
} // namespace GeomAlgoExt
