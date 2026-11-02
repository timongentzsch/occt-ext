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

#include <BRepToolsExt_Checks.hxx>
#include <BRepSurfaceExt_Tools.hxx>
#include <algorithm>
#include <cmath>

namespace BRepExt
{
EvolutionData raiseSurfaceDegree(OcctKernel& k,
                                 uint32_t    id,
                                 int         deltaU,
                                 int         deltaV,
                                 double      tolerance)
{
  const char* op = "raiseSurfaceDegree";
  if (deltaU < 0 || deltaV < 0 || deltaU > 25 || deltaV > 25 || !std::isfinite(tolerance)
      || tolerance <= 0)
    fail(op, invalid_input, "degree increments must be 0..25 and tolerance positive and finite");
  const auto face    = surfacing::faceOf(k, op, id);
  auto       surface = surfacing::editable(op, face, false);
  const int  pu = surface->UDegree() + deltaU, pv = surface->VDegree() + deltaV;
  if (pu > Geom_BSplineSurface::MaxDegree() || pv > Geom_BSplineSurface::MaxDegree())
    fail(op, unsupported_geometry, "the requested degree exceeds OpenCascade's maximum of 25");
  auto before = occ::down_cast<Geom_BSplineSurface>(surface->Copy());
  surface->IncreaseDegree(pu, pv);
  double a, b, c, d;
  surface->Bounds(a, b, c, d);
  for (int i = 0; i <= 64; ++i)
    for (int j = 0; j <= 64; ++j)
    {
      double u = a + (b - a) * i / 64, v = c + (d - c) * j / 64;
      if (surface->Value(u, v).Distance(before->Value(u, v)) > tolerance)
        fail(op, no_solution, "degree elevation changed the surface beyond tolerance");
    }
  return stored(k, surfacing::makeFace(op, surface, tolerance, face.Orientation()));
}
} // namespace BRepExt
