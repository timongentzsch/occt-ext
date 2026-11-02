// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Geometry as a file gives it back: BRepTools and BinTools normalise every direction they read
// and rebuild every frame, so a shape reads back bit for bit only where its directions and frames
// are fixed points of that. Built geometry need not be; these settle it.

#pragma once

#include <Geom2d_Curve.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Dir2d.hxx>
#include <cmath>
#include <type_traits>

namespace BRepExt
{

/// A curve or surface without the trimming OpenCascade wraps it in.
Handle(Geom_Curve)   basisOf(Handle(Geom_Curve) curve);
Handle(Geom2d_Curve) basisOf(Handle(Geom2d_Curve) curve);
Handle(Geom_Surface) basisOf(Handle(Geom_Surface) surface);

/// A direction as a file gives it back: BRepTools and BinTools read its coordinates and normalise
/// them again (SetCoord); this is that normalising's fixed point.
inline bool same(const gp_Dir& a, const gp_Dir& b)
{
  return a.X() == b.X() && a.Y() == b.Y() && a.Z() == b.Z();
}

inline bool same(const gp_Dir2d& a, const gp_Dir2d& b)
{
  return a.X() == b.X() && a.Y() == b.Y();
}

template <class Dir>
Dir readBack(Dir d)
{
  const auto read = [](const Dir& from) {
    Dir again = from;
    if constexpr (std::is_same_v<Dir, gp_Dir>)
      again.SetCoord(from.X(), from.Y(), from.Z());
    else
      again.SetCoord(from.X(), from.Y());
    return again;
  };
  for (int i = 0; i < 8; ++i)
  {
    const Dir again = read(d);
    if (same(again, d))
      return d;
    d = again;
  }
  // Normalising can alternate between two neighbours: then the nearest coordinates, a few
  // units in the last place away, that it keeps.
  for (int ulps = 1; ulps <= 8; ++ulps)
    for (int c = 0; c < (std::is_same_v<Dir, gp_Dir> ? 3 : 2); ++c)
      for (const double towards : {HUGE_VAL, -HUGE_VAL})
      {
        double xyz[3] = {d.X(), d.Y(), 0};
        if constexpr (std::is_same_v<Dir, gp_Dir>)
          xyz[2] = d.Z();
        for (int i = 0; i < ulps; ++i)
          xyz[c] = std::nextafter(xyz[c], towards);
        Dir nudged = d;
        if constexpr (std::is_same_v<Dir, gp_Dir>)
          nudged.SetCoord(xyz[0], xyz[1], xyz[2]);
        else
          nudged.SetCoord(xyz[0], xyz[1]);
        if (same(read(nudged), nudged))
          return nudged;
      }
  return d;
}

/// A frame as a file gives it back: OpenCascade rebuilds a surface's frame from its axis and x
/// direction, read back, re-deriving y and x; this is a fixed point of that rebuilding, with x
/// turned about the axis by at most a few 1e-13 where the frame itself is none, so that new
/// surfaces read back bit for bit, as a file's own analytic surfaces do.
gp_Ax3 fixedFrame(gp_Ax3 frame);

/// The same for a new edge's curves, in space and on its faces: lines' directions and conics'
/// frames (a conic's y comes back re-derived from its x).
void settle(const Handle(Geom_Curve)& stored);

void settle(const Handle(Geom2d_Curve)& stored);

/// Settles the geometry `result` has that `input` does not share: elementary surfaces' frames and
/// the curves of its edges. Geometry the input holds is left as it is, so the input is untouched.
void settleNew(const TopoDS_Shape& input, const TopoDS_Shape& result);

} // namespace BRepExt
