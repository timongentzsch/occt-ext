// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once

// The sides of a Patch loop, and the two routines of it that the host supplies (chaining the
// edges of an open boundary, and closing a gap between runs), as they are not part of the
// library. The host defines them in the namespace below.

#include <GeomAlgoExt_SurfaceJet.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <vector>

namespace BRepExt::patch_host
{

/// An edge of the loop or a curve inside it, run from `from` to `to`.
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

  SurfaceJet jet(double fraction) const
  {
    return surfaceJet("patch", edge, support, parameter(fraction));
  }

  /// The face that the edge was given with, at any order, as its surface and the edge's
  /// curve in that surface's parameters.
  occ::handle<Geom_Surface> held;
  occ::handle<Geom2d_Curve> onHeld;
  bool                      originalReversed = false;

  /// Meet the face at its own edge parameter curve, exactly as faceContinuity reads it.
  /// The support's curve may differ from the spatial edge, especially after a Boolean;
  /// projecting the edge onto the support hides that parameter discrepancy from the report.
  gp_Pnt meetAt(double t) const
  {
    if (held.IsNull())
      return curve->Value(t);
    const gp_Pnt2d uv = onHeld->Value(t);
    return held->Value(uv.X(), uv.Y());
  }

  gp_Pnt meet(double fraction) const { return meetAt(parameter(fraction)); }

  /// A gap between two runs of an open boundary: nothing is asked along it, and its curve is
  /// only the chord that stands for it until the surface is known (see closeGaps).
  bool free = false;
};

struct Rim
{
  occ::handle<Geom2d_BSplineCurve> onSurface;
  double                           deviation = 0;
};

/// The sides of an open boundary in turn, with a gap wherever one run ends and the next begins.
/// Host-defined: all the edges are made one loop.
void chain(std::vector<Side>& sides, double tolerance);

/// Each gap's curve in the surface's parameters. Host-defined; sides that are not gaps keep
/// their rims.
void closeGaps(const std::vector<Side>& sides, std::vector<Rim>& rims);

} // namespace BRepExt::patch_host
