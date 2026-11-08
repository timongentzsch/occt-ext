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

// Blend: how the faces beside an edge turn along all of it, with the sign proven on every part of
// the edge rather than read at sample places.
//
// Along the edge both normals are normal to its tangent d, so h = (n0 x n1) . d / |d| =
// +-sin(theta), theta the angle between the faces' normals (oriented so that h > 0 is convex). With
// the faces' unnormalised normals X = S_u x S_v and the tangent d = d/dt S(u(t), v(t)), each a
// polynomial in t on a part of the edge where both parameter curves and both surfaces are single
// polynomial pieces, H = (X0 x X1) . d is a polynomial with h = H / (|X0| |X1| |d|): its Bernstein
// coefficients bound it with the correlation between the faces kept (faces tangent along the edge
// give coefficients at rounding, however curved the faces), and |X|^2 is bounded below by its least
// coefficient. Rational pieces compose homogeneously (B_i(U / W) W^p = C(p, i) U^i (W - U)^(p-i));
// positive factors drop. A part is halved until it is proven of the edge's sign (h >= -tau convex,
// h <= tau concave, tau the sine of the angle below which the faces are tangent for the blend's
// size, tolerance::tangentFor), proven tangent all over, or its middle shows the other sign (the
// edge is then refused by name); the coefficients reach the values quadratically [Lane and
// Riesenfeld 1981]. Where a parameter curve meets several patches every pair is bounded (each holds
// where it is the face). A plane's normal is constant; a quadric's or torus's is the gradient of
// its implicit equation at the edge's point, a polynomial in the point's homogeneous coordinates
// (that point lies within the edge's tolerance of the surface, so the normal is off by about the
// surface's curvature times that tolerance). The common parameter is a spline face's parameter
// curve's (the edge's), else the edge's curve as a rational B-spline (for a conic another
// parameter, which every quantity then shares). Certified for spline and Bezier surfaces and
// extrusions of spline curves with spline, Bezier or line parameter curves, planes, cylinders,
// cones, spheres and tori; for other representations (revolutions, offsets) the turn is read at the
// places k/16 of the edge (Turn::unproven).

#include <BRepBlendExt_Internal.hxx>
#include <GeomConvert.hxx>
#include <Geom2dConvert.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_BezierCurve.hxx>
#include <Geom2d_Line.hxx>
#include <Geom_BezierSurface.hxx>
#include <Geom_Plane.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Geom_SurfaceOfLinearExtrusion.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_Line.hxx>

#include <deque>
#include <map>

namespace BRepExt::blending::detail
{

namespace
{

//! Polynomials in Bernstein form on [0, 1]: the coefficients, degree size - 1. Their convex hull
//! encloses the polynomial, and halving the interval brings the coefficients to the values
//! quadratically [Farouki and Rajan, On the numerical condition of polynomials in Bernstein form,
//! 1987; Lane and Riesenfeld, Bounds on a polynomial, 1981]. Products and sums keep the
//! correlation between the terms that an interval bound of each would lose.
using Poly = std::vector<double>;

std::vector<double> binomials(int n)
{
  std::vector<double> c(n + 1, 1.0);
  for (int k = 1; k <= n; ++k)
    c[k] = c[k - 1] * (n - k + 1) / k;
  return c;
}

Poly times(const Poly& a, const Poly& b)
{
  const int  n = int(a.size()) - 1, m = int(b.size()) - 1;
  const auto cn = binomials(n), cm = binomials(m), cnm = binomials(n + m);
  Poly       c(n + m + 1, 0.0);
  for (int i = 0; i <= n; ++i)
    if (a[i] != 0)
      for (int j = 0; j <= m; ++j)
        c[i + j] += a[i] * cn[i] * b[j] * cm[j];
  for (int k = 0; k <= n + m; ++k)
    c[k] /= cnm[k];
  return c;
}

//! a + s b, raised to the higher degree (a polynomial times 1 = sum of the basis).
Poly plus(const Poly& a, const Poly& b, double s = 1)
{
  const size_t n  = std::max(a.size(), b.size());
  const Poly   ea = a.size() < n ? times(a, Poly(n - a.size() + 1, 1.0)) : a;
  const Poly   eb = b.size() < n ? times(b, Poly(n - b.size() + 1, 1.0)) : b;
  Poly         c(n);
  for (size_t k = 0; k < n; ++k)
    c[k] = ea[k] + s * eb[k];
  return c;
}

//! The derivative on [0, 1]: n (a_k+1 - a_k).
Poly derivative(const Poly& a)
{
  const int n = int(a.size()) - 1;
  if (n == 0)
    return {0.0};
  Poly d(n);
  for (int k = 0; k < n; ++k)
    d[k] = n * (a[k + 1] - a[k]);
  return d;
}

using P3 = std::array<Poly, 3>;

P3 cross(const P3& a, const P3& b)
{
  return {plus(times(a[1], b[2]), times(a[2], b[1]), -1),
          plus(times(a[2], b[0]), times(a[0], b[2]), -1),
          plus(times(a[0], b[1]), times(a[1], b[0]), -1)};
}

Poly dot(const P3& a, const P3& b)
{
  return plus(plus(times(a[0], b[0]), times(a[1], b[1])), times(a[2], b[2]));
}

//! One polynomial patch of a (rational) B-spline surface in Bernstein form over its own knot span:
//! homogeneous poles (w x, w y, w z, w), row i along u.
struct Patch
{
  int                                p = 0, q = 0;
  std::vector<std::array<double, 4>> poles;
  double                             u0 = 0, u1 = 1, v0 = 0, v1 = 1;

  const std::array<double, 4>& at(int i, int j) const { return poles[size_t(i) * (q + 1) + j]; }
};

//! The derivative patch along u (`alongU`) or v, in the patch's own parameter on [0, 1].
Patch derivativeOf(const Patch& patch, bool alongU)
{
  Patch     out            = patch;
  const int degree         = alongU ? patch.p : patch.q;
  const int rows           = alongU ? std::max(patch.p, 1) : patch.p + 1;
  const int cols           = alongU ? patch.q + 1 : std::max(patch.q, 1);
  (alongU ? out.p : out.q) = std::max(degree - 1, 0);
  out.poles.assign(size_t(rows) * cols, {0, 0, 0, 0});
  if (degree > 0)
    for (int i = 0; i < rows; ++i)
      for (int j = 0; j < cols; ++j)
        for (int k = 0; k < 4; ++k)
          out.poles[size_t(i) * cols + j][k] =
            degree * ((alongU ? patch.at(i + 1, j) : patch.at(i, j + 1))[k] - patch.at(i, j)[k]);
  return out;
}

//! The patch at (u(t), v(t)) for a (rational) parameter curve u = U / W, v = V / W in the patch's
//! own parameters: sum P_ij B_i(u) B_j(v) W^(p+q), homogeneous (a positive factor W^(p+q) on all
//! four coordinates, so the point is unchanged). B_i(U / W) W^p = C(p, i) U^i (W - U)^(p-i).
std::array<Poly, 4> composed(const Patch& patch, const Poly& U, const Poly& V, const Poly& W)
{
  const auto basis = [&](int degree, const Poly& X) {
    std::vector<Poly> power(degree + 1, Poly{1.0}), rest(degree + 1, Poly{1.0});
    const Poly        other = plus(W, X, -1);
    for (int k = 1; k <= degree; ++k)
    {
      power[k] = times(power[k - 1], X);
      rest[k]  = times(rest[k - 1], other);
    }
    const auto        c = binomials(degree);
    std::vector<Poly> out(degree + 1);
    for (int k = 0; k <= degree; ++k)
    {
      out[k] = times(power[k], rest[degree - k]);
      for (auto& x : out[k])
        x *= c[k];
    }
    return out;
  };
  const auto          bu = basis(patch.p, U), bv = basis(patch.q, V);
  std::array<Poly, 4> out;
  for (int i = 0; i <= patch.p; ++i)
    for (int j = 0; j <= patch.q; ++j)
    {
      const Poly term = times(bu[i], bv[j]);
      for (int k = 0; k < 4; ++k)
        out[k] = out[k].empty() ? Poly(term.size(), 0.0) : out[k];
      for (int k = 0; k < 4; ++k)
        for (size_t m = 0; m < term.size(); ++m)
          out[k][m] += patch.at(i, j)[k] * term[m];
    }
  return out;
}

//! A face's unnormalised normal X = S_u x S_v and its edge's tangent d/dt S(u(t), v(t)) along a
//! part, as polynomials with positive factors dropped: with S = A / w, S_u = (w A_u - w_u A) / w^2,
//! and the tangent (w A' - w' A) / w^2.
struct FacePart
{
  P3                  normal, tangent;
  std::array<Poly, 4> point; // homogeneous
};

FacePart facePart(const Patch& patch, const Poly& U, const Poly& V, const Poly& W, double outward)
{
  const auto A = composed(patch, U, V, W), Au = composed(derivativeOf(patch, true), U, V, W),
             Av = composed(derivativeOf(patch, false), U, V, W);
  P3         xu, xv, tangent;
  const Poly dw = derivative(A[3]);
  for (int k = 0; k < 3; ++k)
  {
    xu[k]      = plus(times(A[3], Au[k]), times(Au[3], A[k]), -1);
    xv[k]      = plus(times(A[3], Av[k]), times(Av[3], A[k]), -1);
    tangent[k] = plus(times(A[3], derivative(A[k])), times(dw, A[k]), -1);
  }
  FacePart out{cross(xu, xv), tangent, A};
  for (auto& c : out.normal)
    for (auto& x : c)
      x *= outward;
  return out;
}

Poly scaled(Poly a, double s)
{
  for (auto& x : a)
    x *= s;
  return a;
}

//! The value at the middle (de Casteljau).
double valueAt(Poly a)
{
  for (size_t n = a.size(); n > 1; --n)
    for (size_t k = 0; k + 1 < n; ++k)
      a[k] = (a[k] + a[k + 1]) / 2;
  return a[0];
}

double least(const Poly& a)
{
  return *std::min_element(a.begin(), a.end());
}

double most(const Poly& a)
{
  return *std::max_element(a.begin(), a.end());
}

} // namespace

Turn turnAlong(const char*        op,
               const TopoDS_Edge& edge,
               const TopoDS_Face& first,
               const TopoDS_Face& second,
               double             size)
{
  double                   f, l;
  const Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, f, l);
  if (curve.IsNull())
    fail(op, invalid_input, "a blended edge has no curve", {edge});
  const TopoDS_Face    faces[2] = {first, second};
  Handle(Geom2d_Curve) pcurves[2];
  Handle(Geom_Surface) surfaces[2];
  Side                 sides[2];
  for (int i = 0; i < 2; ++i)
  {
    double a, b;
    pcurves[i] = BRep_Tool::CurveOnSurface(edge, faces[i], a, b);
    if (pcurves[i].IsNull())
      fail(op, invalid_input, "an edge has no parameter curve on a face beside it", {edge});
    surfaces[i] = BRep_Tool::Surface(faces[i]);
    sides[i] =
      Side{faces[i], surfaces[i], faces[i].Orientation() == TopAbs_REVERSED ? -1.0 : 1.0, 0};
  }
  TopAbs_Orientation inFirst = TopAbs_FORWARD;
  for (TopExp_Explorer e(first, TopAbs_EDGE); e.More(); e.Next())
    if (e.Current().IsSame(edge))
      inFirst = e.Current().Orientation();
  const double orientation = inFirst == TopAbs_REVERSED ? -1 : 1;
  // The turn at which the faces are tangent for a blend of this size (tolerance::tangentFor).
  const double tau = std::sin(std::max(tol::coarseAngle, tol::length / size));
  // The turn h at t (convex where (n0 x n1) . d > 0, convexity in BRepBlendExt_Solver.cxx);
  // nothing where a face is singular.
  const auto turn = [&](double t) -> std::optional<double> {
    gp_Pnt at;
    gp_Vec d1;
    curve->D1(t, at, d1);
    if (!(d1.Magnitude() > limit::degenerate))
      return {};
    gp_Vec n[2];
    for (int i = 0; i < 2; ++i)
    {
      const gp_Pnt2d uv = pcurves[i]->Value(t);
      Jet            j;
      if (!jetAt(sides[i], uv.X(), uv.Y(), j))
        return {};
      n[i] = j.n;
    }
    return n[0].Crossed(n[1]).Dot(d1.Normalized()) * orientation;
  };
  // Which representations compose exactly: (rational) spline and Bezier surfaces and planes
  // (as their bilinear splines, which keep their parameters), and parameter curves that are
  // splines, Bezier curves or lines (whose conversion keeps their parameter).
  const auto basis2d = [](Handle(Geom2d_Curve) c) {
    while (const auto t = Handle(Geom2d_TrimmedCurve)::DownCast(c))
      c = t->BasisCurve();
    return c;
  };
  const auto basis3d = [](Handle(Geom_Surface) s) {
    while (const auto t = Handle(Geom_RectangularTrimmedSurface)::DownCast(s))
      s = t->BasisSurface();
    return s;
  };
  // How each face's normal is had along a part, all in one parameter: a spline or Bezier surface
  // through its parameter curve (spline, Bezier or line: their conversion keeps the parameter); a
  // plane's is constant; a quadric's or torus's is the gradient of its implicit equation at the
  // edge's point, a polynomial in the point's homogeneous coordinates (a point within the edge's
  // tolerance of the surface, so off by about its curvature times that tolerance).
  enum class Form
  {
    Spline,
    Plane,
    Implicit,
    Other
  };
  Form form[2];
  bool composes[2];
  for (int i = 0; i < 2; ++i)
  {
    const auto c = basis2d(pcurves[i]);
    composes[i]  = c->IsKind(STANDARD_TYPE(Geom2d_BSplineCurve))
                   || c->IsKind(STANDARD_TYPE(Geom2d_BezierCurve))
                   || c->IsKind(STANDARD_TYPE(Geom2d_Line));
    const auto s = basis3d(surfaces[i]);
    // An extrusion of a spline, Bezier curve or line is a spline in its own parameters.
    Handle(Geom_Curve) profile;
    if (const auto e = Handle(Geom_SurfaceOfLinearExtrusion)::DownCast(s))
    {
      profile = e->BasisCurve();
      while (const auto t = Handle(Geom_TrimmedCurve)::DownCast(profile))
        profile = t->BasisCurve();
    }
    const bool extrusion = !profile.IsNull()
                           && (profile->IsKind(STANDARD_TYPE(Geom_BSplineCurve))
                               || profile->IsKind(STANDARD_TYPE(Geom_BezierCurve))
                               || profile->IsKind(STANDARD_TYPE(Geom_Line)));
    form[i]              = s->IsKind(STANDARD_TYPE(Geom_Plane)) ? Form::Plane
                           : (s->IsKind(STANDARD_TYPE(Geom_BSplineSurface))
                              || s->IsKind(STANDARD_TYPE(Geom_BezierSurface)) || extrusion)
                             ? (composes[i] ? Form::Spline : Form::Other)
                           : (s->IsKind(STANDARD_TYPE(Geom_CylindricalSurface))
                              || s->IsKind(STANDARD_TYPE(Geom_SphericalSurface))
                              || s->IsKind(STANDARD_TYPE(Geom_ConicalSurface))
                              || s->IsKind(STANDARD_TYPE(Geom_ToroidalSurface)))
                             ? Form::Implicit
                             : Form::Other;
  }
  // The carrier of the common parameter: a spline face (its parameter curve's, the edge's), else
  // the edge's curve as a (rational) B-spline (for a conic another parameter, which all the
  // quantities then share).
  const int  carrierFace = form[0] == Form::Spline ? 0 : form[1] == Form::Spline ? 1 : -1;
  const bool certified   = form[0] != Form::Other && form[1] != Form::Other;
  // The gradient of a quadric's or torus's implicit equation at a homogeneous point (A, w), with
  // positive factors dropped: cylinder (P - O) - ((P - O).D) D, sphere P - C, cone Q - (R0 + z
  // tan a) tan a D with Q the part of P - O across the axis and z along it, torus (|P - O|^2 +
  // R^2 - r^2)(P - O) - 2 R^2 Q (derived from |Q|^2 = R^2, |P - C|^2 = R^2, |Q| = R0 + z tan a,
  // (|P - O|^2 + R^2 - r^2)^2 = 4 R^2 |Q|^2).
  const auto gradient = [&](int i, const std::array<Poly, 4>& A) -> P3 {
    const auto s = basis3d(surfaces[i]);
    gp_Ax3     frame;
    if (const auto c = Handle(Geom_ElementarySurface)::DownCast(s))
      frame = c->Position();
    const gp_XYZ O = frame.Location().XYZ(), D = frame.Direction().XYZ();
    P3           rel; // A - O w
    for (int k = 0; k < 3; ++k)
    {
      rel[k] = A[k];
      for (size_t m = 0; m < rel[k].size(); ++m)
        rel[k][m] -= O.Coord(k + 1) * A[3][m];
    }
    if (s->IsKind(STANDARD_TYPE(Geom_SphericalSurface)))
      return rel;
    const Poly z = plus(plus(scaled(rel[0], D.X()), scaled(rel[1], D.Y())), scaled(rel[2], D.Z()));
    P3         across;
    for (int k = 0; k < 3; ++k)
      across[k] = plus(rel[k], scaled(z, D.Coord(k + 1)), -1);
    if (s->IsKind(STANDARD_TYPE(Geom_CylindricalSurface)))
      return across;
    if (const auto cone = Handle(Geom_ConicalSurface)::DownCast(s))
    {
      const double t      = std::tan(cone->SemiAngle());
      const Poly   radius = plus(scaled(A[3], cone->RefRadius()), scaled(z, t));
      P3           out;
      for (int k = 0; k < 3; ++k)
        out[k] = plus(across[k], scaled(radius, t * D.Coord(k + 1)), -1);
      return out;
    }
    const auto   torus = Handle(Geom_ToroidalSurface)::DownCast(s);
    const double R = torus->MajorRadius(), r = torus->MinorRadius();
    const Poly   w2   = times(A[3], A[3]);
    const Poly   lead = plus(dot(rel, rel), scaled(w2, R * R - r * r));
    P3           out;
    for (int k = 0; k < 3; ++k)
      out[k] = plus(times(lead, rel[k]), scaled(times(w2, across[k]), 2 * R * R), -1);
    return out;
  };
  // The sign that turns an implicit gradient into the face's outward normal: their dot product at
  // the edge's middle (they are parallel and the gradient does not vanish on the face; a cosine
  // below one half would mean the surface is not the one its type says, and is not certified).
  double gradientSign[2] = {1, 1};
  bool   converted       = certified;
  for (int i = 0; i < 2 && converted; ++i)
    if (form[i] == Form::Implicit)
    {
      const double   t  = (f + l) / 2;
      const gp_Pnt2d uv = pcurves[i]->Value(t);
      Jet            j;
      gp_Pnt         p;
      curve->D0(t, p);
      std::array<Poly, 4> point = {Poly{p.X()}, Poly{p.Y()}, Poly{p.Z()}, Poly{1.0}};
      const P3            g     = gradient(i, point);
      const double        along = jetAt(sides[i], uv.X(), uv.Y(), j)
                                    ? g[0][0] * j.n.X() + g[1][0] * j.n.Y() + g[2][0] * j.n.Z()
                                    : 0;
      const double length = std::sqrt(g[0][0] * g[0][0] + g[1][0] * g[1][0] + g[2][0] * g[2][0]);
      if (!(std::abs(along) > 0.5 * length))
        converted = false;
      gradientSign[i] = along > 0 ? 1 : -1;
    }
  // Spline faces as B-splines (Bezier converted, which keeps the parameters); planes' normals.
  Handle(Geom_BSplineSurface) splines[2];
  gp_Vec                      planeNormal[2];
  Handle(Geom_BSplineCurve)   edgeSpline;
  if (converted)
    try
    {
      for (int i = 0; i < 2; ++i)
        if (form[i] == Form::Plane)
        {
          gp_Pnt p;
          gp_Vec du, dv;
          surfaces[i]->D1(0, 0, p, du, dv);
          planeNormal[i] = du.Crossed(dv).Normalized() * sides[i].outward;
        }
        else if (form[i] == Form::Spline)
        {
          const auto s = basis3d(surfaces[i]);
          if (s->IsKind(STANDARD_TYPE(Geom_SurfaceOfLinearExtrusion)))
          {
            // Infinite along the extrusion: over the parameter curve's span of v and more.
            const auto pc =
              Geom2dConvert::CurveToBSplineCurve(new Geom2d_TrimmedCurve(pcurves[i], f, l));
            double v0 = HUGE_VAL, v1 = -HUGE_VAL;
            for (int k = 1; k <= pc->NbPoles(); ++k)
            {
              v0 = std::min(v0, pc->Pole(k).Y());
              v1 = std::max(v1, pc->Pole(k).Y());
            }
            splines[i] = GeomConvert::SurfaceToBSplineSurface(
              new Geom_RectangularTrimmedSurface(s, v0 - 1, v1 + 1, false));
          }
          else
            splines[i] = GeomConvert::SurfaceToBSplineSurface(s);
          converted &= !splines[i].IsNull();
        }
      if (carrierFace < 0)
      {
        edgeSpline = GeomConvert::CurveToBSplineCurve(new Geom_TrimmedCurve(curve, f, l));
        converted &= !edgeSpline.IsNull();
      }
    }
    catch (const Standard_Failure&)
    {
      converted = false;
    }
  // The spline's patch over the knot spans [ku, ku+1] x [kv, kv+1], cached.
  std::map<std::array<int, 3>, Patch> patches;
  const auto                          patchOf = [&](int i, int ku, int kv) -> const Patch* {
    const std::array<int, 3> key = {i, ku, kv};
    if (const auto found = patches.find(key); found != patches.end())
      return &found->second;
    const auto& s = splines[i];
    Patch       patch;
    patch.u0                          = s->UKnot(ku);
    patch.u1                          = s->UKnot(ku + 1);
    patch.v0                          = s->VKnot(kv);
    patch.v1                          = s->VKnot(kv + 1);
    Handle(Geom_BSplineSurface) piece = Handle(Geom_BSplineSurface)::DownCast(s->Copy());
    piece->Segment(patch.u0, patch.u1, patch.v0, patch.v1);
    patch.p = piece->UDegree();
    patch.q = piece->VDegree();
    if (piece->NbUPoles() != patch.p + 1 || piece->NbVPoles() != patch.q + 1)
      return nullptr;
    for (int a = 1; a <= piece->NbUPoles(); ++a)
      for (int b = 1; b <= piece->NbVPoles(); ++b)
      {
        const double w = piece->Weight(a, b);
        const gp_XYZ x = piece->Pole(a, b).XYZ() * w;
        patch.poles.push_back({x.X(), x.Y(), x.Z(), w});
      }
    return &(patches[key] = patch);
  };
  // The knot spans of a spline's parameter that meet [lo, hi]; past its ends the end span, whose
  // polynomial runs on.
  const auto spans = [](int count, const auto& knot, double lo, double hi) {
    std::vector<int> out;
    for (int k = 1; k < count; ++k)
      if ((knot(k + 1) >= lo || k + 1 == count) && (knot(k) <= hi || k == 1)
          && knot(k + 1) > knot(k))
        out.push_back(k);
    return out;
  };
  // A spline face's patches along a part [a, b] of its parameter curve.
  const auto faceParts = [&](int i, double a, double b) -> std::optional<std::vector<FacePart>> {
    const auto pc = Geom2dConvert::CurveToBSplineCurve(new Geom2d_TrimmedCurve(pcurves[i], a, b));
    if (pc->NbKnots() != 2)
      return {};
    Poly   U, V, W;
    double u0 = HUGE_VAL, u1 = -HUGE_VAL, v0 = HUGE_VAL, v1 = -HUGE_VAL;
    for (int k = 1; k <= pc->NbPoles(); ++k)
    {
      const double w = pc->Weight(k);
      U.push_back(w * pc->Pole(k).X());
      V.push_back(w * pc->Pole(k).Y());
      W.push_back(w);
      u0 = std::min(u0, pc->Pole(k).X());
      u1 = std::max(u1, pc->Pole(k).X());
      v0 = std::min(v0, pc->Pole(k).Y());
      v1 = std::max(v1, pc->Pole(k).Y());
    }
    std::vector<FacePart> out;
    const auto&           s = splines[i];
    for (const int ku : spans(s->NbUKnots(), [&](int k) { return s->UKnot(k); }, u0, u1))
      for (const int kv : spans(s->NbVKnots(), [&](int k) { return s->VKnot(k); }, v0, v1))
      {
        const Patch* patch = patchOf(i, ku, kv);
        if (!patch)
          return {};
        // The patch's own parameters: (u - u0) / (u1 - u0), homogeneous.
        const Poly Us = scaled(plus(U, W, -patch->u0), 1 / (patch->u1 - patch->u0));
        const Poly Vs = scaled(plus(V, W, -patch->v0), 1 / (patch->v1 - patch->v0));
        out.push_back(facePart(*patch, Us, Vs, W, sides[i].outward));
      }
    return out;
  };

  // The turn's bounds over [a, b] (of the carrier's parameter) from the Bernstein coefficients of
  // H = (X0 x X1) . d for every combination of patches there: h = H / (|X0| |X1| |d|), the norms
  // bounded below by the square roots of their squares' least coefficients, so that h >= min(0,
  // least H / norms) and h <= max(0, most H / norms). {lowest, highest, h at the middle}, or
  // nothing where a norm is not bounded away from zero (a part nearer a face's singular point).
  struct Bounds
  {
    double lowest = 0, highest = 0, middle = 0;
  };

  const auto certify = [&](double a, double b) -> std::optional<Bounds> {
    // The carrier's options: the point (homogeneous) and tangent, and its face's normal.
    std::vector<FacePart> carriers;
    if (carrierFace >= 0)
    {
      const auto parts = faceParts(carrierFace, a, b);
      if (!parts)
        return {};
      carriers = *parts;
    }
    else
    {
      Handle(Geom_BSplineCurve) piece = Handle(Geom_BSplineCurve)::DownCast(edgeSpline->Copy());
      piece->Segment(a, b);
      if (piece->NbKnots() != 2)
        return {};
      FacePart edgePart;
      for (int k = 1; k <= piece->NbPoles(); ++k)
      {
        const double w = piece->Weight(k);
        for (int c = 0; c < 3; ++c)
          edgePart.point[c].push_back(w * piece->Pole(k).Coord(c + 1));
        edgePart.point[3].push_back(w);
      }
      const Poly dw = derivative(edgePart.point[3]);
      for (int c = 0; c < 3; ++c)
        edgePart.tangent[c] = plus(times(edgePart.point[3], derivative(edgePart.point[c])),
                                   times(dw, edgePart.point[c]),
                                   -1);
      carriers.push_back(edgePart);
    }
    const auto normalsOf = [&](int i, const FacePart& carrier) -> std::optional<std::vector<P3>> {
      if (i == carrierFace)
        return std::vector<P3>{carrier.normal};
      if (form[i] == Form::Plane)
        return std::vector<P3>{
          P3{Poly{planeNormal[i].X()}, Poly{planeNormal[i].Y()}, Poly{planeNormal[i].Z()}}};
      if (form[i] == Form::Implicit)
      {
        P3 g = gradient(i, carrier.point);
        for (auto& c : g)
          c = scaled(c, gradientSign[i]);
        return std::vector<P3>{g};
      }
      const auto parts = faceParts(i, a, b);
      if (!parts)
        return {};
      std::vector<P3> out;
      for (const auto& p : *parts)
        out.push_back(p.normal);
      return out;
    };
    Bounds out;
    bool   firstPair = true;
    for (const FacePart& carrier : carriers)
    {
      const auto n0 = normalsOf(0, carrier), n1 = normalsOf(1, carrier);
      if (!n0 || !n1)
        return {};
      for (const P3& x0 : *n0)
        for (const P3& x1 : *n1)
        {
          const Poly s0 = dot(x0, x0), s1 = dot(x1, x1), sd = dot(carrier.tangent, carrier.tangent);
          const double norms = least(s0) * least(s1) * least(sd);
          if (!(norms > 0) || !std::isfinite(norms))
            return {};
          const Poly G = scaled(dot(cross(x0, x1), carrier.tangent), orientation);
          out.lowest   = std::min(out.lowest, least(G) / std::sqrt(norms));
          out.highest  = std::max(out.highest, most(G) / std::sqrt(norms));
          if (firstPair)
            out.middle = valueAt(G) / std::sqrt(valueAt(s0) * valueAt(s1) * valueAt(sd));
          firstPair = false;
        }
    }
    return out;
  };
  // The parts between the knots of the carrier (a certificate holds within one polynomial piece).
  double              from = f, to = l;
  std::vector<double> breaks;
  const auto          knotsOf = [&](const auto& c, double lo, double hi) {
    if (c.IsNull())
      return;
    for (int i = 1; i <= c->NbKnots(); ++i)
      if (c->Knot(i) > lo && c->Knot(i) < hi)
        breaks.push_back(c->Knot(i));
  };
  if (converted && carrierFace < 0)
  {
    from = edgeSpline->FirstParameter();
    to   = edgeSpline->LastParameter();
    knotsOf(edgeSpline, from, to);
  }
  else
    for (int i = 0; i < 2; ++i)
      if (form[i] == Form::Spline)
        knotsOf(Handle(Geom2d_BSplineCurve)::DownCast(basis2d(pcurves[i])), f, l);
  breaks.push_back(from);
  breaks.push_back(to);
  std::sort(breaks.begin(), breaks.end());
  std::deque<std::pair<double, double>> parts;
  for (size_t k = 0; k + 1 < breaks.size(); ++k)
    if (breaks[k + 1] - breaks[k] > tol::parameterMerge * (to - from))
      parts.emplace_back(breaks[k], breaks[k + 1]);
  Turn         out;
  bool         regular = false;
  const double finest  = tol::parameterMerge * (to - from);
  const auto   read    = [&](double where, double h) {
    regular = true;
    if (std::abs(h) <= tau)
    {
      out.touches = true;
      return;
    }
    const int here = h > 0 ? -1 : 1;
    if (out.kind != 0 && here != out.kind)
      fail(op,
           unsupported_geometry,
           "the faces beside the edge turn from convex to concave along it (at "
             + std::to_string(where) + " of its length)",
           {edge});
    out.kind = here;
  };
  // The middle first: a convexity read there before decides as it did.
  if (const auto h = turn((f + l) / 2))
    read(0.5, *h);
  if (!converted)
  {
    // Not certified for these representations: read at the places k / 16 of the edge.
    for (int k = 1; k < 16; ++k)
      if (const auto h = turn(f + (l - f) * k / 16))
        read(k / 16.0, *h);
    out.unproven = 1;
  }
  while (converted && !parts.empty())
  {
    checkCancelled(op);
    const auto [a, b] = parts.front();
    parts.pop_front();
    std::optional<Bounds> bounds;
    try
    {
      bounds = certify(a, b);
    }
    catch (const Standard_Failure&)
    {
      bounds.reset();
    }
    if (bounds)
    {
      read(((a + b) / 2 - from) / (to - from), bounds->middle);
      // Tangent all over the part, or of the edge's sign (h >= -tau convex, h <= tau concave).
      if (bounds->lowest >= -tau && bounds->highest <= tau)
      {
        out.touches = true;
        continue;
      }
      if ((out.kind < 0 && bounds->lowest >= -tau) || (out.kind > 0 && bounds->highest <= tau))
        continue;
    }
    if (b - a <= finest)
    {
      // Unresolved at the parameter resolution: a face's singular point, or a turn that grazes
      // tau there; the places beside it decide (Turn::unproven).
      ++out.unproven;
      continue;
    }
    parts.emplace_back(a, (a + b) / 2);
    parts.emplace_back((a + b) / 2, b);
  }
  if (!regular)
    fail(op, invalid_input, "a face beside the edge is singular there", {edge});
  return out;
}

} // namespace BRepExt::blending::detail
