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
#include <GeomAlgoExt_Failure.hxx>
#include <BndLib_Add3dCurve.hxx>
#include <Bnd_Box.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <Geom_Curve.hxx>
#include <algorithm>
#include <cmath>
#include <gp_Ax3.hxx>
#include <gp_Pln.hxx>
#include <gp_Trsf.hxx>

namespace BRepExt
{
// Bound the actual trimmed curve in plane coordinates. Fixed sample counts miss
// excursions of high-degree splines between samples; OCCT's optimal curve
// bounds include those extrema.
inline void requirePlanarCurve(const char*               op,
                               const Handle(Geom_Curve)& curve,
                               double                    first,
                               double                    last,
                               const gp_Pln&             plane,
                               double                    tolerance)
{
  if (curve.IsNull())
    fail(op, invalid_input, "edge has no 3D curve");
  gp_Trsf local;
  local.SetTransformation(plane.Position());
  const auto        transformed = Handle(Geom_Curve)::DownCast(curve->Transformed(local));
  GeomAdaptor_Curve adaptor(transformed, first, last);
  Bnd_Box           box;
  BndLib_Add3dCurve::AddOptimal(adaptor, first, last, 0.0, box);
  double xmin, ymin, zmin, xmax, ymax, zmax;
  box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
  if (std::max(std::abs(zmin), std::abs(zmax)) > tolerance)
    fail(op, invalid_input, "curve does not lie in the supplied plane");
}
} // namespace BRepExt
