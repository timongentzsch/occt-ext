// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

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

//! A spline surface over a quad cage of points `xyz` (x, y, z in turn) and `quads` (four point
//! indices each), refined `refinement` times (0 to 3). Fails on a cage too large or irregular.
Standard_EXPORT PolySplineResult polySpline(OcctKernel&                kernel,
                                            const std::vector<double>& xyz,
                                            const std::vector<int>&    quads,
                                            double                     tolerance,
                                            int                        refinement,
                                            bool                       merge,
                                            bool                       rounded,
                                            bool                       boundaryExactly,
                                            bool                       requireG1,
                                            double                     angularTolerance);
} // namespace BRepExt
