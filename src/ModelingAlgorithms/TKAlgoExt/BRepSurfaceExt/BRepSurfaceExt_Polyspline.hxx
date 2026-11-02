// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
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

PolySplineResult polySpline(OcctKernel&,
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
