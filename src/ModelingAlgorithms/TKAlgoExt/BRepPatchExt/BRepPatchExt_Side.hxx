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

// The sides of a Patch loop, and the two routines of it for a boundary that does not close
// (chaining its edges, and closing a gap between runs), defined in BRepPatchExt_Chain.cxx.

#include <Standard_Macro.hxx>
#include <GeomAlgoExt_SurfaceJet.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <algorithm>
#include <vector>

namespace BRepExt::patch_side
{

//! An edge of the loop or a curve inside it, run from `from` to `to`.
struct Side
{
  int                     input;
  TopoDS_Edge             edge;
  TopoDS_Face             support;
  int                     order;
  occ::handle<Geom_Curve> curve;
  double                  from, to, length, slack;

  double parameter(double fraction) const { return from + (to - from) * fraction; }

  gp_Pnt at(double fraction) const { return curve->Value(parameter(fraction)); }

  gp_Vec tangent(double fraction) const
  {
    gp_Pnt point;
    gp_Vec tangent;
    curve->D1(parameter(fraction), point, tangent);
    return to < from ? -tangent : tangent;
  }

  //! What the side asks at a place: its face's jet, eased near a corner (see Ease). The ease
  //! falls as 1 - (10 t^3 - 15 t^4 + 6 t^5), t the distance from the corner over its length: the
  //! quintic that goes from one to none with its first two derivatives none at both ends, so
  //! that the bend asked stays twice differentiable along the side and agrees at the corner to
  //! second order with what the side on the other side asks there.
  SurfaceJet jet(double fraction) const
  {
    SurfaceJet jet  = faceJet(fraction);
    double     bend = 0;
    for (int end = 0; end < 2; ++end)
      if (eases[end].length > 0)
      {
        const double t = std::min(1.0, at(fraction).Distance(at(end)) / eases[end].length);
        bend += eases[end].bend * (1 - t * t * t * (10 - 15 * t + 6 * t * t));
      }
    if (bend == 0)
      return jet;
    // Bent across the side by that much, and in no direction along it.
    const gp_Vec across = jet.n.Crossed(tangent(fraction)).Normalized();
    const double a = jet.u.Dot(across), b = jet.v.Dot(across);
    jet.uu += jet.n * (bend * a * a);
    jet.uv += jet.n * (bend * a * b);
    jet.vv += jet.n * (bend * b * b);
    return jet;
  }

  //! The jet of the face itself, uneased.
  SurfaceJet faceJet(double fraction) const
  {
    return surfaceJet("patch", edge, support, parameter(fraction));
  }

  bool eased() const { return eases[0].length > 0 || eases[1].length > 0; }

  //! The largest ease that reaches a place, as a difference of curvature, or 0. Within it the
  //! bend asked runs from the corner's to the face's, and a miss is a share of no less than
  //! that: where both are near flat (straight edges whose faces bend across them), a share of
  //! the bend asked alone would count a rounding as a miss.
  double easing(double fraction) const
  {
    double out = 0;
    for (int end = 0; end < 2; ++end)
      if (eases[end].length > 0 && at(fraction).Distance(at(end)) < eases[end].length)
        out = std::max(out, std::abs(eases[end].bend));
    return out;
  }

  //! The face that the edge was given with, at any order, as its surface and the edge's
  //! curve in that surface's parameters.
  occ::handle<Geom_Surface> held;
  occ::handle<Geom2d_Curve> onHeld;
  bool                      originalReversed = false;

  //! Meet the face at its own edge parameter curve, exactly as faceContinuity reads it.
  //! The support's curve may differ from the spatial edge, especially after a Boolean;
  //! projecting the edge onto the support hides that parameter discrepancy from the report.
  gp_Pnt meetAt(double t) const
  {
    if (held.IsNull())
      return curve->Value(t);
    const gp_Pnt2d uv = onHeld->Value(t);
    return held->Value(uv.X(), uv.Y());
  }

  gp_Pnt meet(double fraction) const { return meetAt(parameter(fraction)); }

  //! A gap between two runs of an open boundary: nothing is asked along it, and its curve is
  //! only the chord that stands for it until the surface is known (see closeGaps).
  bool free = false;

  //! Near a corner where the faces on both sides bend differently, what this side asks across
  //! itself eases from the bend a surface can have at the corner to its face's own: by `bend`
  //! at the corner, as a difference of curvature, to none at `length` from it (see
  //! agreeAtCorners). At the start and at the end of the side; none where `length` is 0.
  struct Ease
  {
    double bend = 0, length = 0;
  };

  Ease eases[2];
};

struct Rim
{
  occ::handle<Geom2d_BSplineCurve> onSurface;
  double                           deviation = 0;
};

//! The sides of an open boundary in turn, with a gap wherever one run ends and the next begins.
//! All the edges are made one loop.
Standard_EXPORT void chain(std::vector<Side>& sides, double tolerance);

//! Each gap's curve in the surface's parameters; sides that are not gaps keep
//! their rims.
Standard_EXPORT void closeGaps(const std::vector<Side>& sides, std::vector<Rim>& rims);

} // namespace BRepExt::patch_side
