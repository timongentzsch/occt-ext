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

#include <Geom_BezierCurve.hxx>
#include <Geom_Curve.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <array>

//! Bridge curves of geometric continuity G0 to G3.
//!
//! A bridge leaves a place of a curve (or a point with a direction) and arrives at another. At an
//! end of geometric continuity G^k its derivatives with respect to the parameter t of the bridge
//! are fixed by the end's arc-length jet (point A, unit tangent T, curvature vector r2 = kappa N,
//! third derivative r3 = -kappa^2 T + kappa' N + kappa tau B) and by up to three free scalars
//! alpha > 0, beta, gamma (Hohmeyer and Barsky, "Rational continuity: parametric, geometric and
//! Frenet frame continuity of rational curves", ACM TOG 8(4), 1989, equations 105 to 109):
//!
//!   B'(0)   = alpha T
//!   B''(0)  = alpha^2 r2 + beta T
//!   B'''(0) = alpha^3 r3 + 3 alpha beta r2 + gamma T
//!
//! The components of B'' and B''' normal to T are therefore fixed by the curve being continued and
//! by alpha and beta; only the T-components of the first three control points past the end remain
//! free. They are the three signed distances x1, x2, x3 along T of those control points from the
//! end, and are what this package calls tensions. A Bezier curve of degree n = k_a + k_b + 1 gives
//! each end its k + 1 control points without sharing, so all six tensions are free.
namespace BRepExt
{

//! What an end of a bridge continues: its point, the unit tangent T pointing into the bridge, the
//! curvature vector r2 and the vector r3 of the third derivative with respect to arc length.
//! `order` is the highest continuity the end can give (0: a point, 1: with its tangent, 2: and its
//! curvature, 3: and its change).
struct BridgeEnd
{
  gp_Pnt point;
  gp_Vec tangent;
  gp_Vec curvature;
  gp_Vec third;
  int    order = 0;
};

//! The end of a bridge that continues `curve` at parameter `u`, leaving against the curve's
//! direction when `against`. The derivatives are those of the curve's own parameter turned into arc
//! length; a curve with fewer than three derivatives gives an end of order 2 or 1, as `order` says.
//! Fails when the curve is singular at `u`.
Standard_EXPORT BridgeEnd bridgeEndOf(const Handle(Geom_Curve)& curve, double u, bool against);

//! The end of a bridge from the first three derivatives of a curve at a point, in any regular
//! parametrisation and in the direction the bridge leaves in.
Standard_EXPORT BridgeEnd bridgeEndFromDerivatives(const gp_Pnt& point,
                                                   const gp_Vec& d1,
                                                   const gp_Vec& d2,
                                                   const gp_Vec& d3);

//! The Bezier curve of degree k_first + k_second + 1 that joins `first` to `second` with G^k_first
//! at the first and G^k_second at the second, where k is the lesser of the order asked and the
//! order the end has. `tensions` are the signed distances x1, x2, x3 of the control points from
//! each end along its tangent, the first end's three then the second's; NaN takes the default x_i =
//! i L / n, L the distance between the ends and n the degree (the uniform parametrisation of a
//! chord). Fails when an x1 is not positive, the ends coincide, or the hodograph vanishes (a cusp).
//! `used` receives the six distances the curve was made with.
Standard_EXPORT Handle(Geom_BezierCurve) bridgeCurve(const BridgeEnd&             first,
                                                     int                          firstOrder,
                                                     const BridgeEnd&             second,
                                                     int                          secondOrder,
                                                     const std::array<double, 6>& tensions,
                                                     std::array<double, 6>&       used);

} // namespace BRepExt
