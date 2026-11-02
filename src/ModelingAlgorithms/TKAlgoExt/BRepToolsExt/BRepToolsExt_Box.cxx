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

// A box by its three sizes, refused below what the kernel can tell apart.
#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <Precision.hxx>
#include <cmath>

namespace BRepExt
{
EvolutionData makeBox(OcctKernel& k, double x, double y, double z)
{
  if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z) || x <= Precision::Confusion()
      || y <= Precision::Confusion() || z <= Precision::Confusion())
    fail("makeBox", invalid_input, "box dimensions must exceed kernel resolution (1e-7)");
  BRepPrimAPI_MakeBox maker(x, y, z);
  maker.Build();
  if (!maker.IsDone())
    fail("makeBox", no_solution, "box dimensions cannot be represented by the kernel");
  return stored(k, maker.Shape());
}
} // namespace BRepExt
