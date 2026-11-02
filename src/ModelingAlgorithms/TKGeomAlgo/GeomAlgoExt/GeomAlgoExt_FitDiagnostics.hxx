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
#include <Geom_Surface.hxx>
#include <gp_Pnt.hxx>
#include <gp_XYZ.hxx>
#include <vector>

namespace BRepExt
{
// Fixed parameters are the input's correspondence to the result, not nearest rendered triangles.
struct FitTarget
{
  gp_Pnt position;
  gp_XYZ normal;
  double u, v;
  int    index;
  bool   hasNormal;
};

struct FitMiss
{
  int    index;
  double residual, normalDegrees;
  bool   hasNormal;
};

struct FitDiagnostics
{
  std::vector<FitMiss> constraints;
  double               maxResidual = 0, maxNormalDegrees = 0, bending = 0;
  double               minOrientedJacobian = 1e100;
  int                  samples             = 1089;
};

// OCCT D2 evaluation, independent of the fitting basis/conditions. Fairness is the sampled
// integral of squared normalized-parametric second derivatives divided by modelScale^2.
// The reference area vector sets both the orientation and the Jacobian's dimensionless scale.
// Samples are an operational check, not a certified global error or intersection bound.
Standard_EXPORT FitDiagnostics measureSurfaceFit(const char*                      operation,
                                                 const occ::handle<Geom_Surface>& surface,
                                                 const std::vector<FitTarget>&    targets,
                                                 double                           modelScale,
                                                 gp_XYZ                           referenceArea);
} // namespace BRepExt
