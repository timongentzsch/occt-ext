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

// The boundary of a fill that does not close.
//
// Chaining: the edges are joined into one loop by greedy growth from both ends of the chain, each
// step taking, among the free edges and the four pairs of ends (either end of the chain, either end
// of the edge), the pair nearest to each other. Greedy nearest-end growth is the standard
// heuristic for assembling unordered curves; it is exact where the ends that belong together are
// nearer each other than to any other end, and it takes the nearest where they are not (a
// branching or a near tie is a matter for the caller's tolerance, not for this routine). A gap
// wherever the end of one run is farther from the start of the next than the tolerance, the loop's
// closing included, is a side of its own that nothing is asked of.
//
// Closing a gap: the gap runs straight between the surface parameters of the two ends it joins,
// the preimage of the runs' ends on either side. A segment in the parameters is the simplest curve
// that joins two points of a chart and keeps the trimmed region a polygon-like domain; its image on
// the surface is whatever the surface makes of it, which the fill then measures.

#include <BRepPatchExt_Side.hxx>
#include <GeomAlgoExt_Failure.hxx>

#include <Geom2d_BSplineCurve.hxx>
#include <Geom_Line.hxx>
#include <NCollection_Array1.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt2d.hxx>

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <vector>

namespace BRepExt::patch_side
{

void chain(std::vector<Side>& sides, double tolerance)
{
  if (sides.empty())
    return;
  const int         originals = static_cast<int>(sides.size());
  std::deque<Side>  run{sides.front()};
  std::vector<Side> rest(sides.begin() + 1, sides.end());
  while (!rest.empty())
  {
    size_t best    = 0;
    double nearest = std::numeric_limits<double>::infinity();
    bool   front = false, far = false;
    for (size_t i = 0; i < rest.size(); ++i)
      for (const bool atFront : {false, true})
        for (const bool farEnd : {false, true})
        {
          const gp_Pnt end      = atFront ? run.front().at(0) : run.back().at(1);
          const double distance = end.Distance(rest[i].at(farEnd ? 1 : 0));
          if (distance < nearest)
          {
            nearest = distance;
            best    = i;
            front   = atFront;
            far     = farEnd;
          }
        }
    Side next = rest[best];
    rest.erase(rest.begin() + static_cast<std::ptrdiff_t>(best));
    // Onto the back the edge starts at the near end; onto the front it ends there.
    if (far != front)
      std::swap(next.from, next.to);
    if (front)
      run.push_front(next);
    else
      run.push_back(next);
  }

  std::vector<Side> loop;
  int               gaps = 0;
  for (size_t k = 0; k < run.size(); ++k)
  {
    const Side& side = run[k];
    const Side& next = run[(k + 1) % run.size()];
    loop.push_back(side);
    const gp_Pnt from     = side.at(1);
    const gp_Pnt to       = next.at(0);
    const double distance = from.Distance(to);
    if (distance <= std::max({tolerance, side.slack, next.slack}))
      continue;
    Side gap{originals + gaps++,
             {},
             {},
             0,
             new Geom_Line(from, gp_Dir(gp_Vec(from, to))),
             0,
             distance,
             distance,
             0,
             {},
             {}};
    gap.free = true;
    loop.push_back(gap);
  }
  sides = std::move(loop);
}

void closeGaps(const std::vector<Side>& sides, std::vector<Rim>& rims)
{
  const size_t count = sides.size();
  for (size_t i = 0; i < count; ++i)
  {
    if (!sides[i].free)
      continue;
    const size_t before = (i + count - 1) % count, after = (i + 1) % count;
    if (sides[before].free || sides[after].free || rims[before].onSurface.IsNull()
        || rims[after].onSurface.IsNull())
      fail("patch", unsupported_geometry, "a gap of the boundary is not between two edges");
    // Where the run before ends and the run after begins, in the surface's parameters.
    const gp_Pnt2d               start = rims[before].onSurface->Value(sides[before].to);
    const gp_Pnt2d               end   = rims[after].onSurface->Value(sides[after].from);
    NCollection_Array1<gp_Pnt2d> poles(1, 2);
    poles(1) = start;
    poles(2) = end;
    NCollection_Array1<double> knots(1, 2);
    knots(1) = sides[i].from;
    knots(2) = sides[i].to;
    NCollection_Array1<int> multiplicities(1, 2);
    multiplicities.Init(2);
    rims[i].onSurface = new Geom2d_BSplineCurve(poles, knots, multiplicities, 1);
    rims[i].deviation = 0;
  }
}

} // namespace BRepExt::patch_side
