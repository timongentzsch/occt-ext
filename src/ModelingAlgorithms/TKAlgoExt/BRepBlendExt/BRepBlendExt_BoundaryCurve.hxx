// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#ifndef BRepBlendExt_BoundaryCurve_HeaderFile
#define BRepBlendExt_BoundaryCurve_HeaderFile

#include <BRepBlendExt_Internal.hxx>

namespace BRepExt::blending::detail
{
//! Enclosure of |C(t) - S(P(t))| for polynomial cubic Bezier curves C and P,
//! t in [0,1]. Evaluates the stored spline basis with outward-rounded Bernstein
//! arithmetic, including every knot cell intersecting P's convex hull. Infinity
//! means that the rational denominator or the finite parameter domain was not
//! certified. Does not measure distance to an implicit rolling-ball surface.
Standard_EXPORT double boundaryCurveGap(const Handle(Geom_BSplineSurface)& surface,
                                        const std::array<gp_Pnt, 4>&       curve,
                                        const std::array<gp_Pnt2d, 4>&     pcurve);

struct BoundaryCurve
{
  TopoDS_Edge edge;
  double      sheetGap = 0, planeGap = 0;
  int         spans = 0;
};

//! Fit one regular, ordered sheet/plane component supplied by continuation.
//! The guide may turn in either sheet coordinate and have both ends on the same
//! contact. Cubic Hermite spans are refined until their whole-span discrepancies
//! from both surface representations are enclosed by the length tolerance.
//! Component discovery, exposed-region selection and support transitions belong
//! to the boundary construction, not to this approximation routine.
Standard_EXPORT std::optional<BoundaryCurve> planeBoundaryCurve(const Sheet&                 sheet,
                                                                const TopoDS_Face&           face,
                                                                const std::vector<gp_Pnt2d>& guide,
                                                                const TopoDS_Vertex&         first,
                                                                const TopoDS_Vertex&         last);
} // namespace BRepExt::blending::detail
#endif
