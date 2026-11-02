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
#include <Bnd_Box.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_BezierSurface.hxx>
#include <GeomConvert_BSplineSurfaceToBezierSurface.hxx>
#include <vector>

namespace GeomAlgoExt
{
inline std::vector<Handle(Geom_BezierSurface)> surfacePatches(
  const Handle(Geom_BSplineSurface)& surface,
  double                             u0,
  double                             u1,
  double                             v0,
  double                             v1)
{
  auto opened = Handle(Geom_BSplineSurface)::DownCast(surface->Copy());
  if (opened->IsUPeriodic())
    opened->SetUNotPeriodic();
  if (opened->IsVPeriodic())
    opened->SetVNotPeriodic();
  GeomConvert_BSplineSurfaceToBezierSurface convert(opened, u0, u1, v0, v1, 1e-10);
  std::vector<Handle(Geom_BezierSurface)>   result;
  for (int i = 1; i <= convert.NbUPatches(); ++i)
    for (int j = 1; j <= convert.NbVPatches(); ++j)
      result.push_back(convert.Patch(i, j));
  return result;
}

// OCCT's rational surfaces have positive weights, so their Bezier pole hulls
// enclose every point. This is shared by tight bounds and collision broad phase.
inline Bnd_Box bezierBounds(const Handle(Geom_BezierSurface)& patch)
{
  Bnd_Box box;
  for (int i = 1; i <= patch->NbUPoles(); ++i)
    for (int j = 1; j <= patch->NbVPoles(); ++j)
      box.Add(patch->Pole(i, j));
  return box;
}
} // namespace GeomAlgoExt
