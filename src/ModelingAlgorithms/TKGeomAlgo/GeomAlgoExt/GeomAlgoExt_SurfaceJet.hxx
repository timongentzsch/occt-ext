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

// A surface at one point to second order: where it is, its tangent plane and how it bends there.
// Continuity across an edge is read by comparing the jets of the faces on both sides, so what
// measures a joint (faceContinuity) and what builds one (patch) share this reading.

#include <GeomAlgoExt_Failure.hxx>

#include <BRepAdaptor_Surface.hxx>
#include <BRep_Tool.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom_Surface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace BRepExt
{

struct SurfaceJet
{
  const char* operation;
  gp_Pnt      p;
  gp_Vec      u, v, uu, vv, uv, n;

  //! The second fundamental form of two tangent vectors given in space.
  double curvature(gp_Vec a, gp_Vec b) const
  {
    const double E = u.SquareMagnitude(), F = u.Dot(v), G = v.SquareMagnitude(),
                 det = E * G - F * F;
    // Independent column scaling removes parameter units from the rank test.
    // Cross-product dual coordinates avoid cancellation in E*G-F*F.
    const double lu = std::hypot(u.X(), u.Y(), u.Z()),
                 lv = std::hypot(v.X(), v.Y(), v.Z());
    if (!(lu > 0) || !(lv > 0) || !std::isfinite(lu) || !std::isfinite(lv))
      fail(operation, no_solution, "singular surface differential");
    const gp_Vec x = u / lu, y = v / lv, z = x.Crossed(y);
    const double sine2 = z.SquareMagnitude();
    constexpr double rankTolerance = 8 * std::numeric_limits<double>::epsilon();
    if (!(sine2 > rankTolerance * rankTolerance))
      fail(operation, no_solution, "singular surface differential");
    auto coordinates = [&](gp_Vec t) {
      // Preserve the arithmetic of the well-scaled path. For extreme UV units,
      // use the dual of the normalized tangent columns instead.
      if (det > 1e-24 && std::isfinite(det) && sine2 > 1e-8)
        return std::pair<double, double>((G * t.Dot(u) - F * t.Dot(v)) / det,
                                         (E * t.Dot(v) - F * t.Dot(u)) / det);
      return std::pair<double, double>(t.Dot(y.Crossed(z)) / sine2 / lu,
                                       t.Dot(z.Crossed(x)) / sine2 / lv);
    };
    auto [au, av] = coordinates(a);
    auto [bu, bv] = coordinates(b);
    return n.Dot(uu * (au * bu) + uv * (au * bv + av * bu) + vv * (av * bv));
  }

  void orient(bool reversed)
  {
    const double lu = std::hypot(u.X(), u.Y(), u.Z()),
                 lv = std::hypot(v.X(), v.Y(), v.Z());
    if (!(lu > 0) || !(lv > 0) || !std::isfinite(lu) || !std::isfinite(lv))
      fail(operation, no_solution, "surface normal is singular");
    const gp_Vec scaled = (u / lu).Crossed(v / lv);
    const double sine = scaled.Magnitude();
    if (!(sine > 8 * std::numeric_limits<double>::epsilon()))
      fail(operation, no_solution, "surface normal is singular");
    n = u.Crossed(v);
    const double length = n.Magnitude();
    if (length >= 1e-15 && std::isfinite(length))
      n.Normalize();
    else
      n = scaled / sine;
    if (reversed)
      n.Reverse();
  }
};

//! Read a one-sided U jet on adjacent distinct knot spans. Keep V in its actual
//! containing span: a broad LocalD2 V range can select an extrapolated span.
//! The caller checks regularity and orients the returned raw derivatives.
inline void surfaceD2OnUSpan(const Handle(Geom_BSplineSurface)& surface,
                             double                             u,
                             double                             v,
                             int                                uFrom,
                             int                                uTo,
                             SurfaceJet&                        jet)
{
  if (surface.IsNull() || !std::isfinite(u) || !std::isfinite(v) || uFrom < 1
      || uFrom >= surface->NbUKnots() || uTo != uFrom + 1 || uTo > surface->NbUKnots()
      || u < surface->UKnot(uFrom) || u > surface->UKnot(uTo) || v < surface->VKnot(1)
      || v > surface->VKnot(surface->NbVKnots()))
    fail(jet.operation, invalid_input, "invalid one-sided surface span");
  int vFrom = 1;
  while (vFrom < surface->NbVKnots() - 1 && surface->VKnot(vFrom + 1) <= v)
    ++vFrom;
  surface->LocalD2(u, v, uFrom, uTo, vFrom, vFrom + 1, jet.p, jet.u, jet.v, jet.uu, jet.vv, jet.uv);
}

//! Cache the adaptor and pcurve once per boundary, including the face's placement.
struct BoundaryJet
{
  const char*          operation;
  BRepAdaptor_Surface  surface;
  Handle(Geom2d_Curve) pc;
  bool                 reversed;

  BoundaryJet(const char* op, const TopoDS_Edge& edge, const TopoDS_Face& face)
      : operation(op),
        surface(face),
        reversed(face.Orientation() == TopAbs_REVERSED)
  {
    double a, b;
    pc = BRep_Tool::CurveOnSurface(edge, face, a, b);
    if (pc.IsNull())
      fail(op, invalid_input, "shared edge has no surface parameter curve");
  }

  SurfaceJet at(double parameter) const
  {
    const auto uv = pc->Value(parameter);
    SurfaceJet j{operation, {}, {}, {}, {}, {}, {}, {}};
    surface.D2(uv.X(), uv.Y(), j.p, j.u, j.v, j.uu, j.vv, j.uv);
    j.orient(reversed);
    return j;
  }
};

inline SurfaceJet surfaceJet(const char*        operation,
                             const TopoDS_Edge& edge,
                             const TopoDS_Face& face,
                             double             parameter)
{
  return BoundaryJet(operation, edge, face).at(parameter);
}

//! The jet of a surface at a parameter pair.
inline SurfaceJet surfaceJet(const char*                      operation,
                             const occ::handle<Geom_Surface>& surface,
                             double                           u,
                             double                           v)
{
  SurfaceJet j{operation, {}, {}, {}, {}, {}, {}, {}};
  surface->D2(u, v, j.p, j.u, j.v, j.uu, j.vv, j.uv);
  j.orient(false);
  return j;
}

//! How differently two jets with one tangent plane bend, as Measure Continuity reads it: the
//! largest difference of normal curvature over all directions, and that as a share of the
//! largest absolute principal curvature of either complete form. A `sign` of -1 reads `b` from
//! the other side, for a surface whose normal is the opposite of its neighbour's.
inline std::pair<double, double> bendDifference(const SurfaceJet& a,
                                                const SurfaceJet& b,
                                                gp_Vec            tangent,
                                                double            sign = 1)
{
  const auto   across    = a.n.Crossed(tangent).Normalized();
  const double ax = a.curvature(tangent, tangent), ay = a.curvature(across, across),
               az = a.curvature(tangent, across), bx = b.curvature(tangent, tangent),
               by = b.curvature(across, across), bz = b.curvature(tangent, across);
  // The spectral norm of a symmetric 2x2 form, including its mixed term.
  const auto norm = [](double x, double y, double z) {
    return std::abs(x / 2 + y / 2) + std::hypot(x / 2 - y / 2, z);
  };
  const double deviation = norm(ax - sign * bx, ay - sign * by, az - sign * bz);
  const double scale = std::max(norm(ax, ay, az), norm(bx, by, bz));
  return {deviation, scale > 0 ? std::min(1.0, deviation / scale) : 0};
}

} // namespace BRepExt
