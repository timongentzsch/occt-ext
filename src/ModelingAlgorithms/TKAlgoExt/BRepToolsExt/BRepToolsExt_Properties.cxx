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

#include <BRepToolsExt_Session.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepAlgoExt_Properties.hxx>
#include <cmath>

namespace BRepExt
{
double volume(OcctKernel& kernel, uint32_t shapeId)
{
  return measure(get(kernel, shapeId), false);
}

std::vector<double> volumeResult(OcctKernel& kernel, uint32_t shapeId)
{
  bool         converged;
  const double value = measure(get(kernel, shapeId), false, &converged);
  return {value, converged ? 1.0 : 0.0};
}

double volumeEstimate(OcctKernel& kernel, uint32_t shapeId)
{
  const double value = estimatedVolume(get(kernel, shapeId));
  if (!std::isfinite(value))
    fail("volumeEstimate", no_solution, "property integration failed");
  return value;
}

double surfaceAreaEstimate(OcctKernel& kernel, uint32_t shapeId)
{
  const double value = estimatedArea(get(kernel, shapeId));
  if (!std::isfinite(value))
    fail("surfaceAreaEstimate", no_solution, "property integration failed");
  return value;
}

double surfaceArea(OcctKernel& kernel, uint32_t shapeId)
{
  return measure(get(kernel, shapeId), true);
}

std::vector<double> surfaceAreaResult(OcctKernel& kernel, uint32_t shapeId)
{
  bool         converged;
  const double value = measure(get(kernel, shapeId), true, &converged);
  return {value, converged ? 1.0 : 0.0};
}

std::vector<double> surfaceCenter(OcctKernel& kernel, uint32_t shapeId)
{
  return surfaceCentre(get(kernel, shapeId));
}

} // namespace BRepExt
