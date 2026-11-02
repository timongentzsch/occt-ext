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
#include <BRepToolsExt_Operations.hxx>
#include <Geom_Surface.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Pln.hxx>

namespace BRepExt
{
// The normalized UVN contract, shared by carried curves and B-rep deformation.
// Swap, then mirror, then scale/offset; N uses oriented face normals.
class DeformUVN
{
  const char*               op;
  gp_Pln                    plane;
  occ::handle<Geom_Surface> surface;
  Vec3                      scale, offset;
  bool                      swap, mirror;
  double                    a, b, c, d, u0, u1, v0, v1, up, out;

public:
  DeformUVN(const char*, const TopoDS_Face&, const TopoDS_Face&, Vec3, Vec3, bool, bool, bool);
  gp_Pnt local(const gp_Pnt&) const;
  gp_Pnt point(const gp_Pnt&) const;
  gp_Pnt world(const gp_Pnt&) const;
  gp_XY  uv(const gp_Pnt&) const;

  double orientation() const
  {
    return out * scale.x * scale.y * scale.z * (swap ? -1 : 1) * (mirror ? -1 : 1);
  }

  void checkChart(double, double, double, double) const;
};

struct DeformResult
{
  EvolutionData evolution;
  double        faceResidual = 0, edgeResidual = 0, pcurveResidual = 0, tangentDegrees = 0,
         minJacobian = 0;
  int spans          = 0;
};

DeformResult deformBody(OcctKernel&,
                        uint32_t,
                        uint32_t,
                        uint32_t,
                        Vec3,
                        Vec3,
                        bool,
                        bool,
                        bool,
                        bool,
                        double,
                        const std::vector<int>&,
                        int);
} // namespace BRepExt
