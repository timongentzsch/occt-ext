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

// Blend: The ball resting on both supports at a place along the edge, and the section through
// its feet.

#include <BRepBlendExt_Smooth.hxx>
#include <BRepBlendExt_Internal.hxx>
#include <Geom_BezierSurface.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_SurfaceOfRevolution.hxx>
#include <gp_Lin.hxx>
#include <Geom2dAdaptor_Curve.hxx>

namespace BRepExt::blending::detail
{

namespace
{
//! A cliff's frame across its edge at `w`: the unit tangent, and two unit normals to it, with
//! their derivatives along w (the normals from a fixed direction, made normal to the tangent).
struct CliffFrame
{
  gp_Pnt p;
  gp_Vec d1, d2, t, n, b, dn, db;
};

bool cliffFrame(const Side& side, double w, CliffFrame& f)
{
  side.cliffCurve->D2(w, f.p, f.d1, f.d2);
  const double speed = f.d1.Magnitude();
  if (!(speed > limit::degenerate))
    return false;
  f.t                  = f.d1 / speed;
  const gp_Vec  dt     = (f.d2 - f.t * f.t.Dot(f.d2)) / speed;
  const gp_Vec& ref    = side.cliffReference;
  const gp_Vec  m      = ref - f.t * f.t.Dot(ref);
  const double  length = m.Magnitude();
  if (!(length > tol::angle))
    return false;
  f.n             = m / length;
  const gp_Vec dm = -(f.t * dt.Dot(ref) + dt * f.t.Dot(ref));
  f.dn            = (dm - f.n * f.n.Dot(dm)) / length;
  f.b             = f.t.Crossed(f.n);
  f.db            = dt.Crossed(f.n) + f.t.Crossed(f.dn);
  return true;
}

} // namespace

double cliffAngle(const Side& side, double w, const gp_Pnt& towards)
{
  CliffFrame f;
  if (!cliffFrame(side, w, f))
    return 0;
  const gp_Vec to(f.p, towards);
  return std::atan2(to.Dot(f.b), to.Dot(f.n));
}

std::array<std::optional<std::pair<double, double>>, 2> runsOn(const Handle(Geom_Surface)& surface)
{
  std::array<std::optional<std::pair<double, double>>, 2> out;
  const auto spline = [](const Handle(Geom_Curve)& curve) {
    const auto basis = basisOf(curve);
    const auto b     = Handle(Geom_BSplineCurve)::DownCast(basis);
    return (!b.IsNull() && !b->IsPeriodic()) || !Handle(Geom_BezierCurve)::DownCast(basis).IsNull();
  };
  const auto basis = basisOf(surface);
  double     u0, u1, v0, v1;
  basis->Bounds(u0, u1, v0, v1);
  bool inU = false, inV = false;
  if (const auto e = Handle(Geom_SurfaceOfLinearExtrusion)::DownCast(basis))
    inU = spline(e->BasisCurve());
  else if (const auto r = Handle(Geom_SurfaceOfRevolution)::DownCast(basis))
    inV = spline(r->BasisCurve());
  else if (const auto b = Handle(Geom_BSplineSurface)::DownCast(basis))
  {
    inU = !b->IsUPeriodic();
    inV = !b->IsVPeriodic();
  }
  else if (!Handle(Geom_BezierSurface)::DownCast(basis).IsNull())
    inU = inV = true;
  if (inU)
    out[0] = std::make_pair(u0, u1);
  if (inV)
    out[1] = std::make_pair(v0, v1);
  return out;
}

std::optional<double> pastPoint(const Handle(Geom_Surface)& surface, double v)
{
  double u0, u1, v0, v1;
  surface->Bounds(u0, u1, v0, v1);
  const bool cone = surface->DynamicType() == STANDARD_TYPE(Geom_ConicalSurface);
  if (v >= v0 && v <= v1 && !cone)
    return {};
  const auto basis = basisOf(surface);
  if (const auto sphere = Handle(Geom_SphericalSurface)::DownCast(basis))
  {
    if (std::cos(v) >= 0)
      return {};
    return std::copysign(M_PI / 2, std::sin(v));
  }
  if (const auto c = Handle(Geom_ConicalSurface)::DownCast(basis))
  {
    const double apex = -c->RefRadius() / std::sin(c->SemiAngle());
    // Past the apex, on the other side of it from the cone's reference circle.
    return apex * v > apex * apex ? std::optional<double>(apex) : std::nullopt;
  }
  if (const auto revolution = Handle(Geom_SurfaceOfRevolution)::DownCast(basis))
  {
    const auto   profile = revolution->BasisCurve();
    const gp_Lin axis(revolution->Location(), revolution->Direction());
    for (const double end : {profile->FirstParameter(), profile->LastParameter()})
      if ((v - end) * (end - (profile->FirstParameter() + profile->LastParameter()) / 2) > 0
          && axis.Distance(profile->Value(end)) <= tol::length)
        return end;
  }
  return {};
}

void supportD2(const Handle(Geom_Surface)& surface,
               double                      u,
               double                      v,
               gp_Pnt&                     p,
               gp_Vec&                     du,
               gp_Vec&                     dv,
               gp_Vec&                     duu,
               gp_Vec&                     dvv,
               gp_Vec&                     duv)
{
  if (const auto point = pastPoint(surface, v))
  {
    supportD2(surface, u + M_PI, 2 * *point - v, p, du, dv, duu, dvv, duv);
    dv.Reverse();
    duv.Reverse();
    return;
  }
  double u0, u1, v0, v1;
  surface->Bounds(u0, u1, v0, v1);
  const auto ends =
    u >= u0 && u <= u1 && v >= v0 && v <= v1 ? decltype(runsOn(surface)){} : runsOn(surface);
  // Within parameter roundoff, retain the one-sided spline jet without
  // moving the evaluation point: snapping its position creates a residual
  // floor for a ball whose section lies just beyond the endpoint.
  const auto jetBase = [](double parameter, const std::pair<double, double>& end) {
    const double boundary = std::clamp(parameter, end.first, end.second);
    const double rounding = tol::rounding(std::max(
      {std::abs(parameter), std::abs(end.first), std::abs(end.second), end.second - end.first}));
    return std::abs(parameter - boundary) <= rounding ? parameter : boundary;
  };
  const double cu = ends[0] ? jetBase(u, *ends[0]) : u;
  const double cv = ends[1] ? jetBase(v, *ends[1]) : v;
  const double eu = u - cu, ev = v - cv;
  if (eu == 0 && ev == 0)
  {
    surface->D2(u, v, p, du, dv, duu, dvv, duv);
    return;
  }
  // S(cu + eu, cv + ev) = S(cu, cv) + eu Su + ev Sv along whichever runs past its end.
  gp_Vec duuu, dvvv, duuv, duvv;
  surface->D3(cu, cv, p, du, dv, duu, dvv, duv, duuu, dvvv, duuv, duvv);
  p.Translate(du * eu + dv * ev);
  if (eu != 0 && ev != 0)
  {
    duu = dvv = duv = gp_Vec(0, 0, 0);
  }
  else if (eu != 0)
  {
    dv += duv * eu;
    dvv += duvv * eu;
    duu = gp_Vec(0, 0, 0);
  }
  else
  {
    du += duv * ev;
    duu += duuv * ev;
    dvv = gp_Vec(0, 0, 0);
  }
}

gp_Pnt supportPoint(const Handle(Geom_Surface)& surface, double u, double v)
{
  double u0, u1, v0, v1;
  surface->Bounds(u0, u1, v0, v1);
  if (u >= u0 && u <= u1 && v >= v0 && v <= v1 && !pastPoint(surface, v))
    return surface->Value(u, v);
  gp_Pnt p;
  gp_Vec du, dv, duu, dvv, duv;
  supportD2(surface, u, v, p, du, dv, duu, dvv, duv);
  return p;
}

bool jetAt(const Side& side, double u, double v, Jet& j)
{
  if (side.cliffCurve)
  {
    // On a cliff: the foot on the edge at u, the normal at angle v about it.
    CliffFrame f;
    if (!cliffFrame(side, u, f))
      return false;
    const double c = std::cos(v), s = std::sin(v);
    j.p   = f.p;
    j.du  = f.d1;
    j.dv  = gp_Vec(0, 0, 0);
    j.duu = f.d2;
    j.dvv = j.duv = gp_Vec(0, 0, 0);
    j.n           = f.n * c + f.b * s;
    j.ndu         = f.dn * c + f.db * s;
    j.ndv         = f.b * c - f.n * s;
    return true;
  }
  if (side.polynomialSpan[0] || side.polynomialSpan[1])
  {
    const auto spline = Handle(Geom_BSplineSurface)::DownCast(basisOf(side.surface));
    if (spline.IsNull())
      return false;
    const int u0 = side.polynomialSpan[0] ? side.polynomialSpan[0] : 1;
    const int u1 = side.polynomialSpan[0] ? u0 + 1 : spline->NbUKnots();
    const int v0 = side.polynomialSpan[1] ? side.polynomialSpan[1] : 1;
    const int v1 = side.polynomialSpan[1] ? v0 + 1 : spline->NbVKnots();
    if (u0 < 1 || u1 > spline->NbUKnots() || v0 < 1 || v1 > spline->NbVKnots())
      return false;
    spline->LocalD2(u, v, u0, u1, v0, v1, j.p, j.du, j.dv, j.duu, j.dvv, j.duv);
  }
  else
    supportD2(side.surface, u, v, j.p, j.du, j.dv, j.duu, j.dvv, j.duv);
  const gp_Vec cross  = j.du.Crossed(j.dv);
  const double length = cross.Magnitude();
  if (!(length > limit::degenerate))
    return false;
  const gp_Vec n  = cross / length;
  const gp_Vec cu = j.duu.Crossed(j.dv) + j.du.Crossed(j.duv),
               cv = j.duv.Crossed(j.dv) + j.du.Crossed(j.dvv);
  // Past the point a surface closes to its parameters run the other way round (pastPoint).
  const double outward = pastPoint(side.surface, v) ? -side.outward : side.outward;
  j.n                  = n * outward;
  j.ndu                = (cu - n * cu.Dot(n)) / length * outward;
  j.ndv                = (cv - n * cv.Dot(n)) / length * outward;
  return true;
}

//! Where the centre of a ball of signed radii r0 and r1 would be at the edge's point E at `t` if
//! its faces were their tangent planes there: off each by its radius, in the plane across the edge
//! (Rail::byCentre).
static std::optional<gp_Pnt> centreAcross(const Rail&   rail,
                                          double        t,
                                          const gp_Pnt& E,
                                          double        r0,
                                          double        r1)
{
  const auto side = [&](int i) -> const Side& {
    return rail.parameterSupports ? (*rail.parameterSupports)[i] : rail.s[i];
  };
  gp_XYZ n[2];
  for (int i = 0; i < 2; ++i)
  {
    double     f, l;
    const auto pcurve = BRep_Tool::CurveOnSurface(rail.edge, side(i).face, f, l);
    if (pcurve.IsNull())
      return {};
    // Evaluate the endpoint span as a cached polynomial when carrying the
    // station origin beyond the edge, rather than extrapolating the de Boor
    // recurrence repeatedly through large signed coefficients.
    const gp_Pnt2d uv = rail.centreInterpolation
      && (t < pcurve->FirstParameter() || t > pcurve->LastParameter())
      ? Geom2dAdaptor_Curve(pcurve).Value(t) : pcurve->Value(t);
    Jet            j;
    if (!jetAt(side(i), uv.X(), uv.Y(), j))
      return {};
    n[i] = j.n.XYZ();
  }
  const double c = n[0].Dot(n[1]), det = 1 - c * c;
  if (!(det > tol::angle))
    return {};
  return gp_Pnt(E.XYZ() + n[0] * ((r0 - c * r1) / det) + n[1] * ((r1 - c * r0) / det));
}

//! The ball of signed radii r0 and r1 off the supports whose section's plane passes through the
//! edge's point at `t` (or the centre across it there, Rail::byCentre).
static bool ballAt(const Rail& rail, double t, Station& s, double r0, double r1)
{
  gp_Pnt E;
  gp_Vec d;
  rail.stationD1(t, E, d);
  if (!(d.Magnitude() > limit::degenerate))
    return false;
  // At an end where the faces are tangent the ball rests on both at the edge's point: the
  // limit of the balls before it, its path there along the edge.
  if (rail.pinchedAt(t))
  {
    Jet j[2];
    for (int i = 0; i < 2; ++i)
    {
      double     f, l;
      const auto pcurve = BRep_Tool::CurveOnSurface(rail.edge, rail.s[i].face, f, l);
      if (pcurve.IsNull())
        return false;
      const gp_Pnt2d uv = pcurve->Value(t);
      if (!jetAt(rail.s[i], uv.X(), uv.Y(), j[i]))
        return false;
      s.uv[i][0]  = uv.X();
      s.uv[i][1]  = uv.Y();
      s.foot[i]   = j[i].p;
      s.normal[i] = j[i].n;
    }
    s.t      = t;
    s.centre = gp_Pnt((j[0].p.XYZ() + j[0].n.XYZ() * r0 + j[1].p.XYZ() + j[1].n.XYZ() * r1) / 2);
    s.along  = d.Normalized();
    return true;
  }
  if (rail.byCentre && rail.spine.IsNull())
  {
    const auto across = centreAcross(rail, t, E, r0, r1);
    if (!across)
      return false;
    E = *across;
  }
  std::array<double, 4> x = {s.uv[0][0], s.uv[0][1], s.uv[1][0], s.uv[1][1]};
  Jet                   a, b;
  gp_XYZ                along;
  // The floor Newton's residual can reach (newton), a forward error bound of its evaluation, unit
  // roundoff u, gamma_k = k u / (1 - k u) [Higham, Accuracy and Stability of Numerical Algorithms,
  // 2002, 3.1 and 3.5]. Each face's first derivatives come with relative error eta <= gamma_k,
  // k their evaluation's operation count (relativeError; derived for a derivative without
  // cancellation among its terms). Then, derived: the computed X = S_u x S_v is off by at most
  // (2 eta + gamma_2) |S_u||S_v|, and the unit normal by 2 (2 eta + gamma_2) s with s = |S_u||S_v|
  // / |S_u x S_v| (|x/|x| - y/|y|| <= 2|x - y|/|x|); n0 x n1 then by e = e0 + e1 + sqrt(6) gamma_2,
  // its direction by 2 e / |n0 x n1|, and arm . along by |arm| 2 e / |n0 x n1| + gamma_3 (|arm| +
  // |c0| + |E|), the last for rounding arm and the product. Near tangent faces the first term is
  // the floor; no constant is chosen. Zero where the error of a face's derivatives is not
  // bounded, and where the section's width times the plane's uncertainty exceeds the solve's
  // tolerance.
  double     rounding = 0;
  const auto eta0 = relativeError(rail.s[0]), eta1 = relativeError(rail.s[1]);
  // The ball whose own section passes through the edge's point: the section of a rolling ball
  // lies in the plane through its centre across the centres' path, whose direction is n0 x n1
  // (the path lies on both offsets). The ball's place along the edge depends then only on the
  // ball's own path, not on how the edge bends (the plane across the edge stalls where the centre
  // reaches the edge's centre of curvature).
  const auto residual =
    [&](const std::array<double, 4>& y, std::array<std::array<double, 5>, 4>& m, double& norm) {
      if (!jetAt(rail.s[0], y[0], y[1], a) || !jetAt(rail.s[1], y[2], y[3], b))
        return false;
      const gp_XYZ c0 = a.p.XYZ() + a.n.XYZ() * r0, c1 = b.p.XYZ() + b.n.XYZ() * r1;
      const gp_XYZ f      = c0 - c1;
      const gp_XYZ cross  = a.n.XYZ().Crossed(b.n.XYZ());
      const double length = cross.Modulus();
      if (!(length > tol::angle))
        return false;
      along            = cross / length;
      const gp_XYZ arm = c0 - E.XYZ();
      const double g   = arm.Dot(along);
      if (eta0 && eta1)
      {
        const double u     = std::numeric_limits<double>::epsilon() / 2;
        const auto   gamma = [u](double k) { return k * u / (1 - k * u); };
        const auto   unit  = [&](const Jet& j, double eta) {
          const double cross = j.du.Crossed(j.dv).Magnitude();
          return 2 * (2 * eta + gamma(2)) * j.du.Magnitude() * j.dv.Magnitude() / cross;
        };
        const double e = unit(a, *eta0) + unit(b, *eta1) + std::sqrt(6.0) * gamma(2);
        // The plane's direction is known to 2 e / |n0 x n1|. That uncertainty turns the section
        // about the centre, moving it by its width times as much: the floor stands only where
        // that stays within the solve's tolerance (a narrow section between all but tangent
        // faces), not where the plane decides a wide section (faces turned against each other,
        // the ball held between them).
        const double uncertain = 2 * e / length;
        rounding = a.p.Distance(b.p) * uncertain <= rail.resolved(limit::solved)
                     ? arm.Modulus() * uncertain
                         + gamma(3) * (arm.Modulus() + c0.Modulus() + E.XYZ().Modulus())
                     : 0;
      }
      // The derivative of the unit direction from that of the cross product.
      const auto turn = [&](const gp_XYZ& dcross) {
        return arm.Dot(dcross - along * along.Dot(dcross)) / length;
      };
      const gp_XYZ ju0 = (a.du + a.ndu * r0).XYZ(), jv0 = (a.dv + a.ndv * r0).XYZ();
      const gp_XYZ ju1 = (b.du + b.ndu * r1).Reversed().XYZ(),
                   jv1 = (b.dv + b.ndv * r1).Reversed().XYZ();
      for (int i = 0; i < 3; ++i)
        m[i] = {ju0.Coord(i + 1),
                jv0.Coord(i + 1),
                ju1.Coord(i + 1),
                jv1.Coord(i + 1),
                -f.Coord(i + 1)};
      m[3] = {ju0.Dot(along) + turn(a.ndu.XYZ().Crossed(b.n.XYZ())),
              jv0.Dot(along) + turn(a.ndv.XYZ().Crossed(b.n.XYZ())),
              turn(a.n.XYZ().Crossed(b.ndu.XYZ())),
              turn(a.n.XYZ().Crossed(b.ndv.XYZ())),
              -g};
      norm = f.Modulus() + std::abs(g);
      return std::isfinite(norm);
    };
  if (!newton(x, rail.resolved(limit::solved), residual, &rounding))
    return false;
  s.t         = t;
  s.uv[0][0]  = x[0];
  s.uv[0][1]  = x[1];
  s.uv[1][0]  = x[2];
  s.uv[1][1]  = x[3];
  s.foot[0]   = a.p;
  s.foot[1]   = b.p;
  s.normal[0] = a.n;
  s.normal[1] = b.n;
  s.centre    = gp_Pnt((a.p.XYZ() + a.n.XYZ() * r0 + b.p.XYZ() + b.n.XYZ() * r1) / 2);
  // The way the ball rolls, as the edge runs.
  s.along = gp_Vec(along.Dot(d.XYZ()) < 0 ? along.Reversed() : along);
  return true;
}

//! A full round's ball (Law::Full): resting on both supports and on the face between them, its
//! centre in the plane across its path through the edge's point at `t`, its radius what that takes.
//! The path lies on the three offsets, across all three normals, so its direction is that of the
//! first support's normal crossed with the middle face's (the supports may be parallel).
static bool fullAt(const Rail& rail, double t, Station& s)
{
  gp_Pnt E;
  gp_Vec d;
  rail.curve.D1(t, E, d);
  if (!(d.Magnitude() > limit::degenerate))
    return false;
  const double          last = s.centre.Distance(s.middle);
  std::array<double, 7> x    = {s.uv[0][0],
                                s.uv[0][1],
                                s.uv[1][0],
                                s.uv[1][1],
                                s.uvMiddle[0],
                                s.uvMiddle[1],
                             last > limit::degenerate * rail.scale ? last : rail.scale};
  Jet                   a, b, m;
  gp_XYZ                along;
  const double          sigma = rail.sigma;
  const auto            residual =
    [&](const std::array<double, 7>& y, std::array<std::array<double, 8>, 7>& J, double& norm) {
      if (!(y[6] > 0) || !jetAt(rail.s[0], y[0], y[1], a) || !jetAt(rail.s[1], y[2], y[3], b)
          || !jetAt(rail.middle, y[4], y[5], m))
        return false;
      const double r  = sigma * y[6];
      const gp_XYZ cM = m.p.XYZ() + m.n.XYZ() * r;
      const gp_XYZ f0 = a.p.XYZ() + a.n.XYZ() * r - cM, f1 = b.p.XYZ() + b.n.XYZ() * r - cM;
      const gp_XYZ cross  = a.n.XYZ().Crossed(m.n.XYZ());
      const double length = cross.Modulus();
      if (!(length > tol::angle))
        return false;
      along             = cross / length;
      const gp_XYZ arm  = a.p.XYZ() + a.n.XYZ() * r - E.XYZ();
      const double g    = arm.Dot(along);
      const auto   turn = [&](const gp_XYZ& dcross) {
        return arm.Dot(dcross - along * along.Dot(dcross)) / length;
      };
      const gp_XYZ ua = (a.du + a.ndu * r).XYZ(), va = (a.dv + a.ndv * r).XYZ();
      const gp_XYZ ub = (b.du + b.ndu * r).XYZ(), vb = (b.dv + b.ndv * r).XYZ();
      const gp_XYZ um = (m.du + m.ndu * r).XYZ(), vm = (m.dv + m.ndv * r).XYZ();
      const gp_XYZ ra = (a.n.XYZ() - m.n.XYZ()) * sigma, rb = (b.n.XYZ() - m.n.XYZ()) * sigma;
      for (int i = 0; i < 3; ++i)
      {
        J[i]     = {ua.Coord(i + 1),
                    va.Coord(i + 1),
                    0,
                    0,
                    -um.Coord(i + 1),
                    -vm.Coord(i + 1),
                    ra.Coord(i + 1),
                    -f0.Coord(i + 1)};
        J[i + 3] = {0,
                    0,
                    ub.Coord(i + 1),
                    vb.Coord(i + 1),
                    -um.Coord(i + 1),
                    -vm.Coord(i + 1),
                    rb.Coord(i + 1),
                    -f1.Coord(i + 1)};
      }
      J[6] = {ua.Dot(along) + turn(a.ndu.XYZ().Crossed(m.n.XYZ())),
              va.Dot(along) + turn(a.ndv.XYZ().Crossed(m.n.XYZ())),
              0,
              0,
              turn(a.n.XYZ().Crossed(m.ndu.XYZ())),
              turn(a.n.XYZ().Crossed(m.ndv.XYZ())),
              sigma * a.n.XYZ().Dot(along),
              -g};
      norm = f0.Modulus() + f1.Modulus() + std::abs(g);
      return std::isfinite(norm);
    };
  if (!newton(x, rail.resolved(limit::solved), residual))
    return false;
  // On the material's side of the face between at the edge's point: across a hole the plane
  // through the edge's point holds a ball on the far side too.
  double     f, l;
  const auto onMiddle = BRep_Tool::CurveOnSurface(rail.edge, rail.middle.face, f, l);
  Jet        at;
  if (onMiddle.IsNull() || !jetAt(rail.middle, onMiddle->Value(t).X(), onMiddle->Value(t).Y(), at)
      || !(gp_Vec(E, m.p.Translated(m.n * (sigma * x[6]))).Dot(at.n) * sigma > 0))
    return false;
  s.t           = t;
  s.uv[0][0]    = x[0];
  s.uv[0][1]    = x[1];
  s.uv[1][0]    = x[2];
  s.uv[1][1]    = x[3];
  s.uvMiddle[0] = x[4];
  s.uvMiddle[1] = x[5];
  s.foot[0]     = a.p;
  s.foot[1]     = b.p;
  s.middle      = m.p;
  s.normal[0]   = a.n;
  s.normal[1]   = b.n;
  s.centre      = m.p.Translated(m.n * (sigma * x[6]));
  s.along       = gp_Vec(along.Dot(d.XYZ()) < 0 ? along.Reversed() : along);
  return true;
}

//! A chordal blend's ball at `t`.
static bool chordAt(const Rail& rail, double t, Station& s, double scale)
{
  // A chordal blend's ball is the one whose feet are the chord apart: its radius from the last
  // ball's (or the chord's), scaled by the chord asked over the chord found until they agree.
  // The chord grows with the radius at nearly its own rate (2 r sin(a / 2), the wedge a changing
  // only as the feet move over curved supports), so each step takes most of the error.
  const double chord = rail.s[0].distance * scale;
  const double last  = s.foot[0].Distance(s.centre);
  double       r     = last > limit::degenerate * rail.scale ? last : chord;
  for (int iteration = 0; iteration < 50; ++iteration)
  {
    Station trial = s;
    if (!ballAt(rail, t, trial, rail.sigma * r, rail.sigma * r))
      return false;
    const double found = trial.foot[0].Distance(trial.foot[1]);
    if (!(found > limit::degenerate * rail.scale))
      return false;
    if (std::abs(found - chord) < rail.resolved(limit::solved))
    {
      s = trial;
      return true;
    }
    r *= chord / found;
  }
  return false;
}

//! distances: setbacks, each foot on its face in the plane across the edge at its point, its
//! distance from that point. Each face is solved on its own: two conditions, its two parameters.
//! Where a foot is still at the edge (the first station), it starts from a setback into its face:
//! into the face from the edge is along the other face's normal, across the edge, for a concave
//! blend, and against it for a convex one.
static bool setbackAt(const Rail& rail, double t, Station& s, double d0, double d1)
{
  gp_Pnt E;
  gp_Vec d;
  rail.curve.D1(t, E, d);
  if (!(d.Magnitude() > limit::degenerate))
    return false;
  const gp_XYZ T            = d.XYZ() / d.Magnitude();
  double       distances[2] = {d0, d1};
  if (rail.angle > 0)
  {
    // A bevel at angle a off the reference face across an opening theta between the faces' tangent
    // planes at the edge: in the triangle of the edge's point and the feet, the setbacks are as the
    // sines of the angles opposite, d' = d sin a / sin(theta + a).
    gp_Vec n[2];
    for (int i = 0; i < 2; ++i)
    {
      double     f, l;
      const auto pcurve = BRep_Tool::CurveOnSurface(rail.edge, rail.s[i].face, f, l);
      Jet        j;
      if (pcurve.IsNull() || !jetAt(rail.s[i], pcurve->Value(t).X(), pcurve->Value(t).Y(), j))
        return false;
      n[i] = j.n;
    }
    const double opening = M_PI - n[0].Angle(n[1]), across = std::sin(opening + rail.angle);
    if (!(across > tol::angle))
      return false;
    distances[1 - rail.reference] = distances[rail.reference] * std::sin(rail.angle) / across;
  }
  Jet at[2], edge[2];
  for (int i = 0; i < 2; ++i)
  {
    double     f, l;
    const auto pcurve = rail.s[i].cliffCurve
                          ? Handle(Geom2d_Curve)()
                          : BRep_Tool::CurveOnSurface(rail.edge, rail.s[i].face, f, l);
    if (pcurve.IsNull() || !jetAt(rail.s[i], s.uv[i][0], s.uv[i][1], at[i])
        || !jetAt(rail.s[i], pcurve->Value(t).X(), pcurve->Value(t).Y(), edge[i]))
      return false;
  }
  for (int i = 0; i < 2; ++i)
  {
    const double r = distances[i];
    if (!(r > 0))
      return false;
    Jet j = at[i];
    // The setback lies across the edge on the face's side: two roots solve the distance, and
    // the one into the face is the chamfer's. The side is read at the edge's point, from both
    // faces' normals there: a foot a setback round a curved face has its normal turned by the
    // setback's angle, which where the faces meet nearly flat or nearly closed turns the
    // direction read at the feet over, towards the root on the face's surface beyond the edge.
    // The start from the edge keeps the feet's direction, on that side.
    const gp_XYZ other = at[1 - i].n.XYZ();
    gp_XYZ       into  = (other - j.n.XYZ() * j.n.XYZ().Dot(other)) * rail.sigma;
    into -= T * T.Dot(into);
    const gp_XYZ otherAtEdge = edge[1 - i].n.XYZ(), own = edge[i].n.XYZ();
    gp_XYZ       side = (otherAtEdge - own * own.Dot(otherAtEdge)) * rail.sigma;
    side -= T * T.Dot(side);
    if (!(into.Modulus() > tol::angle) || !(side.Modulus() > tol::angle))
      return false;
    into.Normalize();
    side.Normalize();
    if (into.Dot(side) < 0)
      into.Reverse();
    if (j.p.Distance(E) < r / 4 || (j.p.XYZ() - E.XYZ()).Dot(side) <= 0)
    {
      GeomAPI_ProjectPointOnSurf onto(gp_Pnt(E.XYZ() + into * r), rail.s[i].surface);
      if (onto.NbPoints() == 0)
        return false;
      double u, v;
      onto.LowerDistanceParameters(u, v);
      // In the same period as the edge's own parameters on the face.
      const auto& surface = rail.s[i].surface;
      s.uv[i][0] =
        surface->IsUPeriodic() ? BRepExt::nearPeriod(u, s.uv[i][0], surface->UPeriod()) : u;
      s.uv[i][1] =
        surface->IsVPeriodic() ? BRepExt::nearPeriod(v, s.uv[i][1], surface->VPeriod()) : v;
      if (!jetAt(rail.s[i], s.uv[i][0], s.uv[i][1], j))
        return false;
    }
    std::array<double, 2> x = {s.uv[i][0], s.uv[i][1]};
    const auto            residual =
      [&](const std::array<double, 2>& y, std::array<std::array<double, 3>, 2>& m, double& norm) {
        if (!jetAt(rail.s[i], y[0], y[1], j))
          return false;
        const gp_XYZ arm    = j.p.XYZ() - E.XYZ();
        const double across = arm.Dot(T), reach = (arm.SquareModulus() - r * r) / (2 * r);
        m[0] = {j.du.XYZ().Dot(T), j.dv.XYZ().Dot(T), -across};
        m[1] = {j.du.XYZ().Dot(arm) / r, j.dv.XYZ().Dot(arm) / r, -reach};
        norm = std::abs(across) + std::abs(reach);
        return std::isfinite(norm);
      };
    if (!newton(x, rail.resolved(limit::solved), residual) || (j.p.XYZ() - E.XYZ()).Dot(side) <= 0)
      return false;
    s.uv[i][0]  = x[0];
    s.uv[i][1]  = x[1];
    s.foot[i]   = j.p;
    s.normal[i] = j.n;
  }
  s.t      = t;
  s.centre = gp_Pnt((s.foot[0].XYZ() + s.foot[1].XYZ()) / 2);
  s.along  = gp_Vec(T);
  return true;
}

//! The angle between the faces' normals at the edge's point at `t`, or nothing where a face has no
//! normal there.
static std::optional<double> turnAt(const Rail& rail, double t)
{
  gp_Vec n[2];
  for (int i = 0; i < 2; ++i)
  {
    double     f, l;
    const auto pcurve = BRep_Tool::CurveOnSurface(rail.edge, rail.s[i].face, f, l);
    Jet        j;
    if (pcurve.IsNull() || !jetAt(rail.s[i], pcurve->Value(t).X(), pcurve->Value(t).Y(), j))
      return {};
    n[i] = j.n;
  }
  // From the sine and cosine both: an angle from its cosine alone is not resolved below
  // sqrt(2 epsilon), some 2e-8 (gp_Vec::Angle).
  return std::atan2(n[0].Crossed(n[1]).Magnitude(), n[0].Dot(n[1]));
}

std::optional<double> relativeError(const Side& side)
{
  if (side.cliffCurve)
    return {};
  const double u     = std::numeric_limits<double>::epsilon() / 2;
  const auto   gamma = [u](double k) { return k * u / (1 - k * u); };
  // Operation counts of the first derivatives: the basis functions' Cox-de Boor recursion over the
  // degrees in u and v (each level a difference, a division, a product and a sum: four), and the
  // derivative's difference and scale (two); elementary surfaces as products of at most three
  // factors with a faithfully rounded sine or cosine.
  const auto surface = basisOf(side.surface);
  if (const auto b = Handle(Geom_BSplineSurface)::DownCast(surface))
    return gamma(4 * (b->UDegree() + b->VDegree()) + 2);
  if (const auto b = Handle(Geom_BezierSurface)::DownCast(surface))
    return gamma(4 * (b->UDegree() + b->VDegree()) + 2);
  if (surface->IsKind(STANDARD_TYPE(Geom_ElementarySurface)))
    return gamma(4);
  return {};
}

bool Rail::tangentAt(double t, bool forSize) const
{
  const auto turn = turnAt(*this, t);
  return turn && (forSize ? SectionPolicy::tangentFor(*turn, scale) : *turn <= tol::angle);
}

double seedPlace(const Rail& rail)
{
  const double middle = (rail.first + rail.last) / 2;
  if (!rail.touching || !rail.tangentAt(middle, true))
    return middle;
  // The places k / n of the rail, as the edge's convexity is decided (BRepBlendExt_Solver.cxx).
  constexpr int n      = 16;
  double        widest = middle, most = 0;
  for (int k = 1; k < n; ++k)
  {
    const double t    = rail.first + (rail.last - rail.first) * k / n;
    const auto   turn = turnAt(rail, t);
    if (turn && *turn > most)
    {
      most   = *turn;
      widest = t;
    }
  }
  return widest;
}

bool contactAt(const Rail& rail, double t, Station& s, double scale)
{
  const Station seed = s;
  // A size that changes along the edge (Rail::sizes) scales both distances.
  const double f = rail.factor(t);
  if (!(f > 0))
    return false; // no ball where the size has fallen to nothing (past an end, Rail::factor)
  if (!(rail.law == Law::Full      ? fullAt(rail, t, s)
        : rail.law == Law::Chordal ? chordAt(rail, t, s, scale)
        : rail.law == Law::Chamfer
          ? setbackAt(rail, t, s, rail.s[0].distance * scale * f, rail.s[1].distance * scale * f)
          : ballAt(rail,
                   t,
                   s,
                   rail.sigma * rail.s[0].distance * scale * f,
                   rail.sigma * rail.s[1].distance * scale * f)))
    return false;
  // Near the point a surface closes to Newton's steps in u are long, over the point and on by whole
  // turns; the feet stay parameters a fit can follow: this side of it, in the turn they were sought
  // in.
  for (int i = 0; i < 2; ++i)
    if (!rail.s[i].cliffCurve)
    {
      const auto& surface = rail.s[i].surface;
      if (const auto point = pastPoint(surface, s.uv[i][1]))
      {
        s.uv[i][0] += M_PI;
        s.uv[i][1] = 2 * *point - s.uv[i][1];
      }
      if (surface->IsUPeriodic())
        s.uv[i][0] = BRepExt::nearPeriod(s.uv[i][0], seed.uv[i][0], surface->UPeriod());
    }
  return true;
}

//! The centre of the circle the section is an arc of: the ball's own, but for a full round, whose
//! radius changes along its path, where the ball touches the blend in a small circle across the
//! path (a canal surface's characteristic circle): the circle through its three feet. (For a
//! chordal round the radius barely changes, and the circle is within 5e-5 of the great circle.)
static std::optional<gp_Pnt> sectionCentre(const Rail& rail, const Station& s)
{
  if (rail.law != Law::Full)
    return s.centre;
  const gp_Vec a(s.middle, s.foot[0]), b(s.middle, s.foot[1]);
  const gp_Vec n  = a.Crossed(b);
  const double n2 = n.SquareMagnitude();
  if (!(n2 > limit::degenerate * std::pow(rail.scale, 4)))
    return {};
  return s.middle.Translated(
    (b.Crossed(n) * a.SquareMagnitude() + n.Crossed(a) * b.SquareMagnitude()) / (2 * n2));
}

// Rotate only within each seam's compact geometric neighborhood. Disjoint windows
// preserve the natural section between two crossings and at untouched chain ends.
static bool frameTangents(const Rail&                      rail,
                          const Station&                   s,
                          const std::vector<SectionFrame>& frames,
                          gp_Vec                           tangent[2],
                          bool                             framed[2])
{
  for (int side = 0; side < 2; ++side)
  {
    double total = 0;
    gp_Vec target;
    for (const auto& frame : frames)
      if (frame.side == side)
      {
        const double distance =
          railLength(rail, s.t, frame.t);
        if (!(frame.reach > limit::degenerate * rail.scale))
          return false;
        const double x      = std::clamp(distance / frame.reach, 0.0, 1.0);
        const double weight = 1 - x * x * x * (10 + x * (-15 + 6 * x));
        if (!(weight > 0))
          continue;
        gp_Vec into = frame.tangent - s.normal[side] * s.normal[side].Dot(frame.tangent);
        if (!std::isfinite(into.SquareMagnitude()) || !(into.Magnitude() > tol::angle))
          return false;
        into.Normalize();
        target += into * weight;
        total += weight;
      }
    if (total > 0)
    {
      framed[side]  = true;
      tangent[side] = tangent[side] * std::max(0.0, 1 - total) + target;
      if (!std::isfinite(tangent[side].SquareMagnitude())
          || !(tangent[side].Magnitude() > tol::angle))
        return false;
      tangent[side].Normalize();
    }
  }
  return true;
}

//! A rational quadratic split at its middle into two (Rail::halves): de Casteljau on the
//! homogeneous poles, so the curve is the same.
static void halve(Section& out)
{
  const double w  = out.weights[1];
  const gp_XYZ A0 = out.poles[0].XYZ(), A1 = out.poles[1].XYZ() * w, A2 = out.poles[2].XYZ();
  const double b = (1 + w) / 2, m = (1 + w) / 2;
  const gp_XYZ B1 = (A0 + A1) / 2, M = (A0 + A1 * 2 + A2) / 4, C1 = (A1 + A2) / 2;
  out.poles          = {out.poles[0], gp_Pnt(B1 / b), gp_Pnt(M / m), gp_Pnt(C1 / b), out.poles[2]};
  out.weights        = {1, b, m, b, 1};
  out.knots          = {0, 0.5, 1};
  out.multiplicities = {3, 2, 3};
}

bool sectionAt(const Rail& rail, const Station& s, Section& out)
{
  return sectionAt(rail, s, out, {});
}

bool sectionsBesideSupportKnot(const Rail&                      rail,
                               const Station&                   seed,
                               int                              side,
                               int                              parameter,
                               int                              knot,
                               Section                          out[2],
                               const std::vector<SectionFrame>& frames)
{
  if (side < 0 || side > 1 || parameter < 0 || parameter > 1 || rail.s[side].cliffCurve)
    return false;
  const auto spline = Handle(Geom_BSplineSurface)::DownCast(basisOf(rail.s[side].surface));
  if (spline.IsNull())
    return false;
  const int count = parameter ? spline->NbVKnots() : spline->NbUKnots();
  if (knot < 2 || knot >= count)
    return false;
  const int continuity = parameter ? spline->VDegree() - spline->VMultiplicity(knot)
                                   : spline->UDegree() - spline->UMultiplicity(knot);
  if (continuity < 2)
    return false;
  for (int end = 0; end < 2; ++end)
  {
    Rail branch = rail;
    branch.knotWindows.reset();
    // Continue all uses of this carrier, including the source-normal station
    // plane. Keeping that plane on the original branch leaves an independent
    // frame break even if the support-curvature scalar has been smoothed.
    const auto continueSide = [&](Side& support) {
      if (basisOf(support.surface) == spline)
        support.polynomialSpan[parameter] = knot - 1 + end;
    };
    for (auto& support : branch.s)
      continueSide(support);
    if (branch.parameterSupports)
      for (auto& support : *branch.parameterSupports)
        continueSide(support);
    Station station = seed;
    if (!contactAt(branch, seed.t, station) || !sectionAt(branch, station, out[end], frames))
      return false;
  }
  const auto& a = out[0];
  const auto& b = out[1];
  if (a.degree != b.degree || a.knots != b.knots || a.multiplicities != b.multiplicities
      || a.poles.size() != b.poles.size() || a.weights.size() != a.poles.size()
      || b.weights.size() != b.poles.size())
    return false;
  for (size_t j = 0; j < a.poles.size(); ++j)
    if (!(a.weights[j] > 0 && b.weights[j] > 0))
      return false;
  return true;
}

bool sectionAcrossSupportKnot(const Rail&                      rail,
                              const Station&                   seed,
                              int                              side,
                              int                              parameter,
                              int                              knot,
                              double                           position,
                              Section&                         out,
                              const std::vector<SectionFrame>& frames)
{
  Section sections[2];
  if (!std::isfinite(position) || position < 0 || position > 1
      || !sectionsBesideSupportKnot(rail, seed, side, parameter, knot, sections, frames))
    return false;
  const auto&  a     = sections[0];
  const auto&  b     = sections[1];
  const double blend = smoothActivation(position);
  Section      result = a;
  for (size_t j = 0; j < a.poles.size(); ++j)
  {
    const double wa = (1 - blend) * a.weights[j], wb = blend * b.weights[j];
    const double weight = wa + wb;
    // Anchored homogeneous interpolation retains small differences at a
    // translated model origin, and preserves positive section weights.
    result.poles[j] = a.poles[j].Translated(gp_Vec(a.poles[j], b.poles[j]) * (wb / weight));
    result.weights[j] = weight;
  }
  out = std::move(result);
  return true;
}

bool sectionAt(const Rail&                      rail,
               const Station&                   s,
               Section&                         out,
               const std::vector<SectionFrame>& frames)
{
  if (rail.law == Law::G2 && rail.knotWindows)
    for (const auto& window : *rail.knotWindows)
    {
      const double x = s.uv[window.side][window.parameter] - window.value;
      if (std::abs(x) < window.half)
        return sectionAcrossSupportKnot(rail, s, window.side, window.parameter, window.knot,
                                        (x + window.half) / (2 * window.half), out, frames);
    }
  if (rail.law == Law::Chamfer)
  {
    // The straight line between the feet.
    if (!(s.foot[0].Distance(s.foot[1]) > limit::degenerate * rail.scale))
      return false;
    out                = Section{};
    out.poles          = {s.foot[0], s.foot[1]};
    out.weights        = {1, 1};
    out.knots          = {0, 1};
    out.multiplicities = {2, 2};
    out.degree         = 1;
    return true;
  }
  // Where the ball rests on both faces at one place (Rail::pinch) the section is that place, its
  // weights those the sections before it tend to (the wedge between the faces a half turn).
  if (rail.pinchedAt(s.t) && (rail.law == Law::Conic || rail.law == Law::Chordal))
  {
    out                = Section{};
    out.poles          = {s.foot[0], s.foot[0], s.foot[0]};
    out.weights        = {1, rho::conicWeight(rail.rho, M_PI), 1};
    out.knots          = {0, 1};
    out.multiplicities = {3, 3};
    out.degree         = 2;
    if (rail.halves)
      halve(out);
    return true;
  }
  const auto found = sectionCentre(rail, s);
  if (!found)
    return false;
  const gp_Pnt centre = *found;
  if (rail.law == Law::Full)
  {
    // The ball's circle in two arcs, from the first foot to where it touches the face between
    // and on to the second foot (together up to a whole half turn and more).
    out                  = Section{};
    const gp_Pnt feet[3] = {s.foot[0], s.middle, s.foot[1]};
    out.poles            = {feet[0]};
    out.weights          = {1};
    for (int i = 0; i < 2; ++i)
    {
      const gp_Vec e0(centre, feet[i]), e1(centre, feet[i + 1]);
      const double r = e0.Magnitude(), angle = e0.Angle(e1);
      if (!(r > limit::degenerate * rail.scale) || !(angle > tol::angle)
          || !(angle < M_PI - tol::angle))
        return false;
      const gp_Vec bisector = (e0 / r + e1.Normalized()).Normalized();
      out.poles.push_back(centre.Translated(bisector * (r / std::cos(angle / 2))));
      out.poles.push_back(feet[i + 1]);
      out.weights.push_back(std::cos(angle / 2));
      out.weights.push_back(1);
    }
    out.knots          = {0, 0.5, 1};
    out.multiplicities = {3, 2, 3};
    out.degree         = 2;
    return true;
  }
  const gp_Vec e0(centre, s.foot[0]), e1(centre, s.foot[1]);
  if (!(e0.Magnitude() > limit::degenerate * rail.scale)
      || !(e1.Magnitude() > limit::degenerate * rail.scale))
    return false;
  const gp_Vec u0 = e0.Normalized(), u1 = e1.Normalized();
  gp_Vec       t0 = u1 - u0 * u1.Dot(u0), t1 = u0 - u1 * u0.Dot(u1);
  if (!(t0.Magnitude() > tol::angle) || !(t1.Magnitude() > tol::angle))
    return false;
  t0.Normalize();
  t1.Normalize();
  // The apex: where the tangents at the two feet meet, in the plane of the centre and feet.
  const gp_Vec between(s.foot[0], s.foot[1]);
  const double c = t0.Dot(t1), det = 1 - c * c;
  // Where the faces are all but tangent (a place Rail::touching tells of, or one between the
  // places convexity is decided at) the tangents are opposite to within sqrt(fineAngle) and the
  // lines through the feet nearly one: any point between the feet is within |F1 - F0|
  // sqrt(fineAngle) of both lines (1e-6 of the section's width), and the apex's limit as the
  // wedge opens to a half turn is their middle.
  const bool flat =
    !(det > tol::fineAngle) && c < 0 && (rail.law == Law::Conic || rail.law == Law::Chordal);
  if (!(det > tol::fineAngle) && !flat)
    return false;
  const double a = flat ? between.Magnitude() / 2 : (between.Dot(t0) - c * between.Dot(t1)) / det;
  const double b = flat ? a : (c * between.Dot(t0) - between.Dot(t1)) / det;
  if (!(a > 0 && b > 0))
    return false;
  const gp_Pnt apex = s.foot[0].Translated(t0 * a);
  out               = Section{};
  switch (rail.law)
  {
    case Law::Chamfer:
    case Law::Conic:
    case Law::Chordal:
      out.poles = {s.foot[0], apex, s.foot[1]};
      // The supports' wedge at the apex is between the tangents back to the feet.
      out.weights        = {1, rho::conicWeight(rail.rho, std::acos(std::clamp(c, -1.0, 1.0))), 1};
      out.knots          = {0, 1};
      out.multiplicities = {3, 3};
      out.degree         = 2;
      if (rail.conicWindow)
      {
        const auto& window = *rail.conicWindow;
        const double x = conicPhase(rail, s) / window.halfWidth;
        if (std::abs(x) < 1)
        {
          // Both shoulders lie on both contact tangent planes. Their affine
          // combination keeps the two contact positions and tangent planes exact.
          const gp_Pnt q = s.foot[window.cylinder];
          const double d = gp_Vec(q, window.planeOrigin).Dot(window.planeNormal)
                           / window.axis.Dot(window.planeNormal);
          const gp_Pnt flat = q.Translated(window.axis * d);
          const double w = (1 - x * x) * (1 - x * x);
          out.poles[1] = apex.Translated(gp_Vec(apex, flat) * w);
        }
      }
      if (rail.halves)
        halve(out);
      return true;
    case Law::G2: {

      // curved supports its section as measured (curvedG2).
      gp_Vec tangent[2] = {t0, t1};
      bool   framed[2]  = {false, false};
      if (!frameTangents(rail, s, frames, tangent, framed))
        return false;
      double bend[2]       = {0, 0};
      double activation[2] = {0, 0};
      for (int i = 0; i < 2; ++i)
      {
        // Past a surface's end (runsOn) the ruled run-on is outside the faces: the section keeps
        // the curvature at the end, as its zero there would make the section jump. A choice held

        double foot[2] = {s.uv[i][0], s.uv[i][1]};
        if (!rail.s[i].cliffCurve)
        {
          const auto ends = runsOn(rail.s[i].surface);
          for (int k = 0; k < 2; ++k)
            if (ends[k])
              foot[k] = std::clamp(foot[k], ends[k]->first, ends[k]->second);
        }
        Jet j;
        if (!jetAt(rail.s[i], foot[0], foot[1], j))
          return false;
        bend[i] = SurfaceJet{rail.op, j.p, j.du, j.dv, j.duu, j.dvv, j.duv, s.normal[i]}.curvature(
          tangent[i],
          tangent[i]);
        // A support is flat where its curvature moves the section by no more than the fit
        // tolerance (it shifts the third pole by (5/4) kappa h^2, h under the ball's radius):
        // a plane carried as a spline, whose curvature is rounding, gives the planar section.
        const double relativeBend =
          std::abs(bend[i]) * s.centre.SquareDistance(s.foot[i]) / fitTolerance(rail);
        activation[i] = smoothActivation(relativeBend - 1);
        // Each side must ease independently even while the other side is curved.
        bend[i] *= activation[i];
      }
      for (const auto& [along, across] : rho::g2Poles(rail.rho))
        out.poles.push_back(
          s.foot[0].Translated(gp_Vec(s.foot[0], apex) * along + gp_Vec(apex, s.foot[1]) * across));
      const bool activeFrame = framed[0] || framed[1];
      if (activeFrame)
      {
        // A frame may rotate a calibrated profile only where its curved
        // construction has the same flat limit. Certify the geometry instead
        // of silently switching between incompatible high-rho profiles.
        auto         flatLimit   = out.poles;
        const gp_Vec natural[2]  = {t0, t1};
        const double flatBend[2] = {0, 0};
        int          flatRegime  = 0;
        if (!SectionPolicy::curvedG2(rail.rho, s, apex, natural, flatBend, flatLimit, flatRegime))
          return false;
        for (size_t j = 0; j < flatLimit.size(); ++j)
          if (!(flatLimit[j].Distance(out.poles[j]) <= fitTolerance(rail)))
            fail(rail.op,
                 no_solution,
                 "the G2 section has incompatible flat-support limits",
                 {rail.edge},
                 FailureKind::Curvature);
      }
      const double activity =
        activeFrame ? 1 : activation[0] + activation[1] - activation[0] * activation[1];
      if (activity > 0)
      {
        const auto planar = out.poles;
        if (!SectionPolicy::curvedG2(rail.rho, s, apex, tangent, bend, out.poles, out.regime))
          return false;
        // Blend the entire profile, not just its curvature: host policies may have
        // a different flat limit. Zero endpoint jets of activity make this C2 too.
        if (activity < 1)
          for (size_t j = 0; j < out.poles.size(); ++j)
            out.poles[j] = planar[j].Translated(gp_Vec(planar[j], out.poles[j]) * activity);
      }
      out.weights.assign(8, 1);
      out.knots.assign(rho::g2Knots.begin(), rho::g2Knots.end());
      out.multiplicities.assign(rho::g2Multiplicities.begin(), rho::g2Multiplicities.end());
      out.degree = 5;
      return true;
    }
    case Law::Full: // above
      break;
  }
  return false;
}

std::optional<Ball> restingBall(const std::vector<Side>& faces,
                                const gp_Pnt&            seed,
                                double                   distance,
                                double                   sigma,
                                double                   scale)
{
  const int  count = int(faces.size());
  Ball       ball{seed, std::vector<gp_Pnt>(count), std::vector<gp_Vec>(count)};
  const auto feet = [&] {
    for (int i = 0; i < count; ++i)
    {
      GeomAPI_ProjectPointOnSurf onto(ball.centre, faces[i].surface);
      if (onto.NbPoints() == 0)
        return false;
      double u, v;
      onto.LowerDistanceParameters(u, v);
      Jet j;
      if (!jetAt(faces[i], u, v, j))
        return false;
      ball.feet[i]    = j.p;
      ball.normals[i] = j.n;
    }
    return true;
  };
  if (count < 3 || !feet())
    return {};
  for (int iteration = 0; iteration < 50; ++iteration)
  {
    math_Matrix m(1, count, 1, 3);
    math_Vector rhs(1, count);
    double      worst = 0;
    for (int i = 0; i < count; ++i)
    {
      for (int j = 0; j < 3; ++j)
        m(i + 1, j + 1) = ball.normals[i].Coord(j + 1);
      const double off = gp_Vec(ball.feet[i], ball.centre).Dot(ball.normals[i]) - sigma * distance;
      rhs(i + 1)       = -off;
      worst            = std::max(worst, std::abs(off));
    }
    if (worst < tol::resolved(limit::solved, scale, seed.XYZ().Modulus()))
      return ball;
    // Three faces: the step that zeroes all three; more: the one nearest doing so.
    math_Matrix square = count == 3 ? m : m.Transposed() * m;
    math_Vector right  = count == 3 ? rhs : m.Transposed() * rhs, step(1, 3);
    math_Gauss  solve(square);
    if (!solve.IsDone())
      return {};
    solve.Solve(right, step);
    // Least squares that no longer move: the ball rests on all the faces where it meets each
    // within the length resolution, as a vertex lies on a face (a ball off one by more is not
    // that face's: the corner parts, Blend::vertexBlend).
    if (count > 3 && step.Norm() < tol::resolved(limit::solved, scale, seed.XYZ().Modulus()))
      return worst <= limit::one(scale, seed.XYZ().Modulus()) ? std::optional(ball) : std::nullopt;
    ball.centre.Translate(gp_Vec(step(1), step(2), step(3)));
    if (!feet())
      return {};
  }
  return {};
}

} // namespace BRepExt::blending::detail
