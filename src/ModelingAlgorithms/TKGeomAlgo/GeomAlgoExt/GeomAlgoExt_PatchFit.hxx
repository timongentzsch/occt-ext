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

// A tensor-product B-spline surface whose poles are found by least squares: conditions sampled
// along curves are met as closely as they are asked to be, and the freedom they leave goes to
// the fairest surface, the one that bends least and whose bend changes least over its parameter
// rectangle.

#include <Standard_Macro.hxx>
#include <GeomAlgoExt_Band.hxx>

namespace BRepExt
{
using GeomAlgoExt::Band;
}

#include <Geom_BSplineSurface.hxx>
#include <gp_XYZ.hxx>

#include <functional>
#include <optional>
#include <vector>

namespace BRepExt
{

//! Degrees up to 5; flat knots; the pole (i, j) at `poles[i * countV() + j]`.
struct SplineNet
{
  int                 degreeU = 0, degreeV = 0;
  std::vector<double> knotsU, knotsV;
  std::vector<gp_XYZ> poles;

  int countU() const { return static_cast<int>(knotsU.size()) - degreeU - 1; }

  int countV() const { return static_cast<int>(knotsV.size()) - degreeV - 1; }
};

//! The basis functions of one direction that are not zero at a parameter, from the one numbered
//! `first`, with their derivatives up to the third.
constexpr int splineOrder = 6, splineDerivatives = 3;

struct SplineBasis
{
  int    first;
  double value[splineDerivatives + 1][splineOrder];
};

Standard_EXPORT SplineBasis splineBasis(const char*                operation,
                                        const std::vector<double>& knots,
                                        int                        degree,
                                        double                     t);

//! One scalar condition at (u, v): the sum over `count` terms of
//! `direction . d^(du+dv) S / du^du dv^dv` is to be `value`. Derivatives up to the second.
struct FitCondition
{
  struct Term
  {
    int    du, dv;
    gp_XYZ direction;
  };

  double u, v;
  Term   terms[5];
  int    count;
  double value;
  //! Against the energies, which count as one per unit of squared second derivative.
  double weight;
  //! A wish and not a condition: weighed against the energies and left at that, where a
  //! condition is met as closely as the net can.
  bool soft = false;
  //! How far from its value the condition may stay and count as met.
  double within = 0;
  //! Read and not asked: the fit does nothing for it, and goes on until it is met with the
  //! others. For places between those that are asked, where the surface is on its own.
  bool watch = false;
};

//! A wish at one place of a line of the net across which the surface is only once
//! differentiable: that its second derivative across the line be the same on both sides there.
//! The line is `u = at` with `acrossU`, else `v = at`; `along` is the other parameter.
struct FitJump
{
  bool   acrossU;
  double at, along;
  //! Against the energies, as for a condition.
  double weight;
};

struct FitEnergy
{
  //! Of the integral of the squared second and of the squared third derivatives.
  double bending, variation;
};

//! The energies of a fit as assembled, before any condition: the next fit of the same net from
//! the same poles, with the same `fixed`, energy, density, `along` and jumps, starts from them.
//! Empty before the first.
struct FitEnergies
{
  std::optional<Band> matrix;
  std::vector<double> constant;
  bool                separate = false;
};

//! Move the poles of `net` that are not `fixed` (one flag per pole) from where they are given.
//! With `along`, a direction of unit length, they move along it only: the net is then a height
//! over the plane that its poles were given in, a third of the unknowns and a ninth of the
//! work. `density` is what the energies count for at a parameter pair, one everywhere when it
//! is empty. The energies are of the parameters: a surface that is steep or stretched over them
//! pays for that. The conditions that are not `soft` are met until each is `within` its own
//! distance of its value, or the net can follow them no better. Returns which of the
//! conditions that are only read (`watch`) are not met in the end, by their place in
//! `conditions`. `jumps` are weighed against the energies like the wishes among the conditions.
//! Fails with a std::runtime_error when the system has no solution. `leastSquares` keeps only the
//! unique initial weighted fair fit for a height graph; its caller must correct the original
//! boundary and certify all delivered geometry. The default retains the multiplier correction used
//! by supported surfaces. `again`, where given, keeps the energies as assembled for the next fit of
//! the same net from the same poles with other conditions (see FitEnergies).
//! Optional pole prior: add influence times the mean free-pole energy diagonal times
//! sum |P - prior|^2. The prior has the same basis and ordering as net. It is applied
//! after cached energies, so iterative fits may keep one original reference.
Standard_EXPORT std::vector<int> fitNet(const char*                                  operation,
                                        SplineNet&                                   net,
                                        const std::vector<char>&                     fixed,
                                        const std::vector<FitCondition>&             conditions,
                                        const FitEnergy&                             energy,
                                        const std::function<double(double, double)>& density,
                                        const gp_XYZ*               along        = nullptr,
                                        const std::vector<FitJump>& jumps        = {},
                                        bool                        leastSquares = false,
                                        FitEnergies*                again        = nullptr,
                                        const std::vector<gp_XYZ>* polePrior    = nullptr,
                                        double                     poleInfluence = 0);

Standard_EXPORT occ::handle<Geom_BSplineSurface> surfaceOf(const SplineNet& net);

} // namespace BRepExt
