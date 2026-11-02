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

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>

#include <BRepAdaptor_Surface.hxx>
#include <BRepTools.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>
#include <TopoDS.hxx>

#include <cmath>

namespace BRepExt
{

// How fair a face is, from the face alone: its principal curvatures are taken at the middles of
// a grid of its parameter rectangle, where they are inside the face, and summed over its area.
std::vector<double> faceFairness(OcctKernel& kernel, uint32_t faceId, int samples)
{
  const char* const op = "faceFairness";
  if (get(kernel, faceId).ShapeType() != TopAbs_FACE || samples < 4 || samples > 512)
    fail(op, invalid_input, "a face and 4 to 512 samples each way are required");
  const TopoDS_Face face = TopoDS::Face(get(kernel, faceId));
  double            u0, u1, v0, v1;
  BRepTools::UVBounds(face, u0, u1, v0, v1);
  BRepAdaptor_Surface     surface(face);
  BRepTopAdaptor_FClass2d classifier(face, BRepExt::tolerance::uv);

  struct Place
  {
    bool   in = false;
    gp_Pnt point;
    double k1 = 0, k2 = 0, area = 0;
  };

  std::vector<Place> grid(static_cast<size_t>(samples) * samples);
  const double       du = (u1 - u0) / samples, dv = (v1 - v0) / samples;
  double             area = 0, bending = 0, variation = 0, sharpest = 0;
  for (int i = 0; i < samples; ++i)
    for (int j = 0; j < samples; ++j)
    {
      const double u = u0 + (i + 0.5) * du, v = v0 + (j + 0.5) * dv;
      if (classifier.Perform(gp_Pnt2d(u, v)) != TopAbs_IN)
        continue;
      Place& place = grid[static_cast<size_t>(i) * samples + j];
      gp_Vec su, sv, suu, svv, suv;
      surface.D2(u, v, place.point, su, sv, suu, svv, suv);
      const gp_Vec normal  = su.Crossed(sv);
      const double stretch = normal.Magnitude();
      if (stretch < 1e-15)
        continue;
      const double E = su.Dot(su), F = su.Dot(sv), G = sv.Dot(sv);
      const double L = suu.Dot(normal) / stretch, M = suv.Dot(normal) / stretch,
                   N      = svv.Dot(normal) / stretch;
      const double mean   = (E * N - 2 * F * M + G * L) / (2 * stretch * stretch);
      const double gauss  = (L * N - M * M) / (stretch * stretch);
      const double spread = std::sqrt(std::max(0.0, mean * mean - gauss));
      place.in            = true;
      place.k1            = mean + spread;
      place.k2            = mean - spread;
      place.area          = stretch * du * dv;
      area += place.area;
      bending += (place.k1 * place.k1 + place.k2 * place.k2) * place.area;
      sharpest = std::max({sharpest, std::abs(place.k1), std::abs(place.k2)});
    }
  // How fast the curvatures change from each place to the next along u and along v.
  for (int i = 0; i < samples; ++i)
    for (int j = 0; j < samples; ++j)
    {
      const Place& place = grid[static_cast<size_t>(i) * samples + j];
      for (const auto [di, dj] : {std::pair{1, 0}, std::pair{0, 1}})
      {
        if (!place.in || i + di >= samples || j + dj >= samples)
          continue;
        const Place& next = grid[static_cast<size_t>(i + di) * samples + j + dj];
        const double step = place.point.Distance(next.point);
        if (!next.in || step < 1e-15)
          continue;
        const double d1 = (next.k1 - place.k1) / step, d2 = (next.k2 - place.k2) / step;
        variation += (d1 * d1 + d2 * d2) * (place.area + next.area) / 2;
      }
    }
  return {area, bending, variation * area, sharpest};
}

} // namespace BRepExt
