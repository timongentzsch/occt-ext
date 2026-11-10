// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#ifndef BRepBlendExt_EndBoundary_HeaderFile
#define BRepBlendExt_EndBoundary_HeaderFile

#include <BRepBlendExt_Internal.hxx>
#include <gp_Pln.hxx>

namespace BRepExt::blending::detail
{
//! A point of a section/plane intersection, on the original support surfaces.
//! Coordinates are (u0,v0,u1,v1,w), where w is the section parameter.
struct EndBoundaryPoint
{
  std::array<double, 5> parameters{};
  gp_Pnt                point, centre, contact[2];
  double                ballResidual = 0, planeResidual = 0;
};

struct EndBoundaryTrace
{
  enum class End
  {
    Contact,
    SurfaceBoundary,
    Closed,
    Singular,
    Stalled,
    InvalidSeed
  };
  std::vector<EndBoundaryPoint> points;
  End                           end                = End::InvalidSeed;
  int                           boundaryCoordinate = -1;
  //! All simultaneous contact/carrier restrictions, using the same residual allowance.
  std::vector<int> boundaryCoordinates;
};

//! Experimental continuation of one regular section/plane component from a contact root.
//! An interior seed requires a direction, for continuing across a shared support restriction.
//! The three offset-equality equations and section/plane equation define a curve in five
//! parameters. No opposite-contact end assumption is made. A supplied seed selects the
//! component; this does not discover other components or classify the exposed swept region.
//! Constant-radius conic/G2 sections on regular surfaces only. Sample spacing is a numerical
//! policy, not a curve approximation certificate. Each returned point solves the equations
//! to the existing station tolerance. This is not yet used by solid attachment.
Standard_EXPORT EndBoundaryTrace tracePlaneEnd(const Rail&    rail,
                                               const gp_Pln&  plane,
                                               const Station& seed,
                                               double         sectionParameter,
                                               double         spacing,
                                               const gp_Vec&  direction = gp_Vec());
} // namespace BRepExt::blending::detail
#endif
