// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
#include <Standard_Macro.hxx>
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
  bool                      swap = false, mirror = false;
  double                    a = 0, b = 0, c = 0, d = 0, u0 = 0, u1 = 0, v0 = 0, v1 = 0, up = 0,
                            out = 0;

public:
  //! The map from the planar `reference` face to `target`. Fails naming `operation` on a
  //! reference that is not planar or a zero scale in u or v.
  DeformUVN(const char*        operation,
            const TopoDS_Face& reference,
            const TopoDS_Face& target,
            Vec3               scale,
            Vec3               offset,
            bool               swapUV,
            bool               mirror,
            bool               flipNormal);
  //! `point` in the reference plane's frame: x, y and height along its oriented normal.
  gp_Pnt local(const gp_Pnt& point) const;
  //! `point` carried onto the target and its normal. Fails where the normal is singular.
  gp_Pnt point(const gp_Pnt& point) const;
  //! The inverse of `local`.
  gp_Pnt world(const gp_Pnt& point) const;
  //! The target parameters `point` maps to: swapped, mirrored, then scaled and offset.
  gp_XY uv(const gp_Pnt& point) const;

  double orientation() const
  {
    return out * scale.x * scale.y * scale.z * (swap ? -1 : 1) * (mirror ? -1 : 1);
  }

  //! Fails where the parameter range crosses the seam of a periodic target.
  void checkChart(double uMin, double uMax, double vMin, double vMax) const;
};

struct DeformResult
{
  EvolutionData evolution;
  double        faceResidual = 0, edgeResidual = 0, pcurveResidual = 0, tangentDegrees = 0,
         minJacobian = 0;
  int spans          = 0;
};

//! `body` carried from the planar `reference` face onto `target` by the UVN map (DeformUVN),
//! with its residuals. Fails with `reblend`, above 16 faces or for a zero scale along the normal.
Standard_EXPORT DeformResult deformBody(OcctKernel&             kernel,
                                        uint32_t                body,
                                        uint32_t                reference,
                                        uint32_t                target,
                                        Vec3                    scale,
                                        Vec3                    offset,
                                        bool                    swapUV,
                                        bool                    mirror,
                                        bool                    flipNormal,
                                        bool                    reblend,
                                        double                  tolerance,
                                        const std::vector<int>& hashes,
                                        int                     bound);
} // namespace BRepExt
