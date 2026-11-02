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
#include <Standard_Macro.hxx>
#include <BRepToolsExt_Session.hxx>
#include <vector>

namespace BRepExt
{
struct PolySplineResult
{
  uint32_t resultId = 0;
  int      patches = 0, extraordinary = 0, refinement = 0, freeEdges = 0;
  int      samples = 0, seams = 0;
  double   deviation = 0, gap = 0, normalAngle = 0, minJacobian = 0;
  bool     regular = false, solid = false;
  // Original cage face for each output face (not an identity claim across
  // regeneration).
  std::vector<int> provenance;
};

Standard_EXPORT PolySplineResult polySpline(OcctKernel&,
                                            const std::vector<double>&,
                                            const std::vector<int>&,
                                            double tolerance,
                                            int    refinement,
                                            bool   merge,
                                            bool   rounded,
                                            bool   boundaryExactly,
                                            bool   requireG1,
                                            double angularTolerance);
} // namespace BRepExt
