// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
#include <Geom_Surface.hxx>
#include <gp_Pnt.hxx>
#include <gp_XYZ.hxx>
#include <vector>

namespace BRepExt {
// Fixed parameters are the input's correspondence to the result, not nearest rendered triangles.
struct FitTarget {
    gp_Pnt position;
    gp_XYZ normal;
    double u, v;
    int index;
    bool hasNormal;
};
struct FitMiss {
    int index;
    double residual, normalDegrees;
    bool hasNormal;
};
struct FitDiagnostics {
    std::vector<FitMiss> constraints;
    double maxResidual = 0, maxNormalDegrees = 0, bending = 0;
    double minOrientedJacobian = 1e100;
    int samples = 1089;
};
// OCCT D2 evaluation, independent of the fitting basis/conditions. Fairness is the sampled
// integral of squared normalized-parametric second derivatives divided by modelScale^2.
// The reference area vector sets both the orientation and the Jacobian's dimensionless scale.
// Samples are an operational check, not a certified global error or intersection bound.
FitDiagnostics measureSurfaceFit(const char *operation, const occ::handle<Geom_Surface> &surface,
                                 const std::vector<FitTarget> &targets, double modelScale,
                                 gp_XYZ referenceArea);
} // namespace BRepExt
