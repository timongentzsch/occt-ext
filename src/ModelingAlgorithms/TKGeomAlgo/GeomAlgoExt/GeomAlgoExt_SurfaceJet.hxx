// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

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
#include <utility>

namespace BRepExt {

struct SurfaceJet {
    const char *operation;
    gp_Pnt p;
    gp_Vec u, v, uu, vv, uv, n;
    /// The second fundamental form of two tangent vectors given in space.
    double curvature(gp_Vec a, gp_Vec b) const {
        const double E = u.SquareMagnitude(), F = u.Dot(v), G = v.SquareMagnitude(), det = E * G - F * F;
        if (det <= 1e-24)
            fail(operation, no_solution, "singular surface differential");
        auto coordinates = [&](gp_Vec t) {
            return std::pair<double, double>((G * t.Dot(u) - F * t.Dot(v)) / det,
                                             (E * t.Dot(v) - F * t.Dot(u)) / det);
        };
        auto [au, av] = coordinates(a);
        auto [bu, bv] = coordinates(b);
        return n.Dot(uu * (au * bu) + uv * (au * bv + av * bu) + vv * (av * bv));
    }
    void orient(bool reversed) {
        n = u.Crossed(v);
        if (n.Magnitude() < 1e-15)
            fail(operation, no_solution, "surface normal is singular");
        n.Normalize();
        if (reversed)
            n.Reverse();
    }
};

/// Read a one-sided U jet on adjacent distinct knot spans. Keep V in its actual
/// containing span: a broad LocalD2 V range can select an extrapolated span.
/// The caller checks regularity and orients the returned raw derivatives.
inline void surfaceD2OnUSpan(const Handle(Geom_BSplineSurface)& surface, double u, double v,
                             int uFrom, int uTo, SurfaceJet& jet) {
    if (surface.IsNull() || !std::isfinite(u) || !std::isfinite(v) || uFrom < 1 ||
        uFrom >= surface->NbUKnots() || uTo != uFrom + 1 || uTo > surface->NbUKnots() ||
        u < surface->UKnot(uFrom) || u > surface->UKnot(uTo) ||
        v < surface->VKnot(1) || v > surface->VKnot(surface->NbVKnots()))
        fail(jet.operation, invalid_input, "invalid one-sided surface span");
    int vFrom = 1;
    while (vFrom < surface->NbVKnots() - 1 && surface->VKnot(vFrom + 1) <= v)
        ++vFrom;
    surface->LocalD2(u, v, uFrom, uTo, vFrom, vFrom + 1,
                     jet.p, jet.u, jet.v, jet.uu, jet.vv, jet.uv);
}

/// Cache the adaptor and pcurve once per boundary, including the face's placement.
struct BoundaryJet {
    const char* operation;
    BRepAdaptor_Surface surface;
    Handle(Geom2d_Curve) pc;
    bool reversed;
    BoundaryJet(const char* op, const TopoDS_Edge& edge, const TopoDS_Face& face)
        : operation(op), surface(face), reversed(face.Orientation() == TopAbs_REVERSED) {
        double a, b;
        pc = BRep_Tool::CurveOnSurface(edge, face, a, b);
        if (pc.IsNull()) fail(op, invalid_input, "shared edge has no surface parameter curve");
    }
    SurfaceJet at(double parameter) const {
        const auto uv = pc->Value(parameter);
        SurfaceJet j{operation, {}, {}, {}, {}, {}, {}, {}};
        surface.D2(uv.X(), uv.Y(), j.p, j.u, j.v, j.uu, j.vv, j.uv);
        j.orient(reversed);
        return j;
    }
};

inline SurfaceJet surfaceJet(const char *operation, const TopoDS_Edge &edge, const TopoDS_Face &face,
                             double parameter) {
    return BoundaryJet(operation, edge, face).at(parameter);
}

/// The jet of a surface at a parameter pair.
inline SurfaceJet surfaceJet(const char *operation, const occ::handle<Geom_Surface> &surface, double u,
                             double v) {
    SurfaceJet j{operation, {}, {}, {}, {}, {}, {}, {}};
    surface->D2(u, v, j.p, j.u, j.v, j.uu, j.vv, j.uv);
    j.orient(false);
    return j;
}

/// How differently two jets with one tangent plane bend, as Measure Continuity reads it: the
/// largest difference of normal curvature over all directions, and that as a share of the
/// largest curvature along and across `tangent` on either side. A `sign` of -1 reads `b` from
/// the other side, for a surface whose normal is the opposite of its neighbour's.
inline std::pair<double, double> bendDifference(const SurfaceJet &a, const SurfaceJet &b, gp_Vec tangent,
                                                double sign = 1) {
    const auto across = a.n.Crossed(tangent).Normalized();
    const double aa = a.curvature(tangent, tangent) - sign * b.curvature(tangent, tangent);
    const double bb = a.curvature(across, across) - sign * b.curvature(across, across);
    const double ab = a.curvature(tangent, across) - sign * b.curvature(tangent, across);
    const double deviation = std::abs((aa + bb) / 2) + std::hypot((aa - bb) / 2, ab);
    const double scale =
        std::max({std::abs(a.curvature(tangent, tangent)), std::abs(a.curvature(across, across)),
                  std::abs(b.curvature(tangent, tangent)), std::abs(b.curvature(across, across))});
    return {deviation, scale > 1e-15 ? std::min(1.0, deviation / scale) : 0};
}

} // namespace BRepExt
