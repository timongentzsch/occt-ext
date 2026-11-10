// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepOffsetExt_Trim.hxx>
#include <GeomAlgoExt_Interval.hxx>
#include <Geom2dAdaptor_Curve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <array>
#include <gp_Pln.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Cone.hxx>
#include <gp_Sphere.hxx>
#include <gp_Torus.hxx>
#include <gp_Circ.hxx>
#include <gp_Lin.hxx>
#include <gp_Lin2d.hxx>

namespace BRepExt
{
namespace
{
namespace I  = GeomAlgoExt::interval;
using Vector = std::array<I::Range, 3>;

Vector vector(const gp_XYZ& xyz)
{
  return {I::point(xyz.X()), I::point(xyz.Y()), I::point(xyz.Z())};
}

Vector add(Vector a, Vector b)
{
  for (int i = 0; i < 3; ++i)
    a[i] = I::add(a[i], b[i]);
  return a;
}

Vector scale(Vector a, I::Range b)
{
  for (auto& x : a)
    x = I::multiply(x, b);
  return a;
}

double magnitude(I::Range x)
{
  return I::valid(x) ? std::max(std::abs(x.low), std::abs(x.high))
                     : std::numeric_limits<double>::infinity();
}

//! Taylor enclosures, with an interval pi for range reduction. For |x| <= 4,
//! the omitted terms after degree 49/48 are bounded by the next absolute term
//! (the alternating tails decrease). 4^50/50! < epsilon^2 for binary64.
//! No accuracy assumption about libm sin/cos.
std::pair<I::Range, I::Range> sincos(double t)
{
  if (t == 0)
    return {I::point(0), I::point(1)};
  const double k = std::nearbyint(t / (2 * M_PI));
  const auto   x =
    I::subtract(I::point(t), I::multiply(I::point(2 * k), {I::roundDown(M_PI), I::roundUp(M_PI)}));
  if (magnitude(x) > 4)
    return {I::invalid(), I::invalid()};
  const auto square = I::multiply(x, x);
  auto       st = x, ct = I::point(1), s = st, c = ct;
  for (int j = 1; j <= 24; ++j)
  {
    st = I::dividePositive(I::multiply(st, I::multiply(square, I::point(-1))),
                           I::point((2.0 * j) * (2.0 * j + 1)));
    ct = I::dividePositive(I::multiply(ct, I::multiply(square, I::point(-1))),
                           I::point((2.0 * j - 1) * (2.0 * j)));
    s  = I::add(s, st);
    c  = I::add(c, ct);
  }
  const double sr = magnitude(I::dividePositive(I::multiply(st, square), I::point(50.0 * 51)));
  const double cr = magnitude(I::dividePositive(I::multiply(ct, square), I::point(49.0 * 50)));
  return {I::add(s, {-sr, sr}), I::add(c, {-cr, cr})};
}

//! S along an elementary iso is A+B*t or A+B*cos(t)+C*sin(t).
struct Coefficients
{
  Vector a, b, c;
};

std::optional<Coefficients> coefficients(const GeomAdaptor_Surface& surface,
                                         bool                       uIso,
                                         double                     fixed)
{
  gp_Ax3     axes;
  double     radius = 0, minor = 0, angle = 0;
  const auto type = surface.GetType();
  switch (type)
  {
    case GeomAbs_Plane:
      axes = surface.Plane().Position();
      break;
    case GeomAbs_Cylinder:
      axes   = surface.Cylinder().Position();
      radius = surface.Cylinder().Radius();
      break;
    case GeomAbs_Cone:
      axes   = surface.Cone().Position();
      radius = surface.Cone().RefRadius();
      angle  = surface.Cone().SemiAngle();
      break;
    case GeomAbs_Sphere:
      axes   = surface.Sphere().Position();
      radius = surface.Sphere().Radius();
      break;
    case GeomAbs_Torus:
      axes   = surface.Torus().Position();
      radius = surface.Torus().MajorRadius();
      minor  = surface.Torus().MinorRadius();
      break;
    default:
      return std::nullopt;
  }
  const auto   x = vector(axes.XDirection().XYZ()), y = vector(axes.YDirection().XYZ());
  const auto   z      = vector(axes.Direction().XYZ());
  auto         origin = vector(axes.Location().XYZ());
  const Vector zero   = vector(gp_XYZ(0, 0, 0));
  if (type == GeomAbs_Plane)
    return Coefficients{add(origin, scale(uIso ? x : y, I::point(fixed))), uIso ? y : x, zero};
  const auto [s, c] = (!uIso && (type == GeomAbs_Cylinder || type == GeomAbs_Cone))
                        ? std::pair<I::Range, I::Range>{I::point(0), I::point(1)}
                        : sincos(fixed);
  const auto radial = add(scale(x, c), scale(y, s));
  const auto r      = I::point(radius);
  if (uIso)
  {
    if (type == GeomAbs_Cylinder)
      return Coefficients{add(origin, scale(radial, r)), z, zero};
    if (type == GeomAbs_Cone)
    {
      const auto [sa, ca] = sincos(angle);
      return Coefficients{add(origin, scale(radial, r)),
                          add(scale(radial, sa), scale(z, ca)),
                          zero};
    }
    if (type == GeomAbs_Sphere)
      return Coefficients{origin, scale(radial, r), scale(z, r)};
    return Coefficients{add(origin, scale(radial, r)),
                        scale(radial, I::point(minor)),
                        scale(z, I::point(minor))};
  }
  auto at = r;
  if (type == GeomAbs_Cylinder)
    origin = add(origin, scale(z, I::point(fixed)));
  if (type == GeomAbs_Cone)
  {
    const auto [sa, ca] = sincos(angle);
    origin              = add(origin, scale(z, I::multiply(I::point(fixed), ca)));
    at                  = I::add(r, I::multiply(I::point(fixed), sa));
  }
  if (type == GeomAbs_Sphere)
  {
    origin = add(origin, scale(z, I::multiply(r, s)));
    at     = I::multiply(r, c);
  }
  if (type == GeomAbs_Torus)
  {
    origin = add(origin, scale(z, I::multiply(I::point(minor), s)));
    at     = I::add(r, I::multiply(I::point(minor), c));
  }
  return Coefficients{origin, scale(x, at), scale(y, at)};
}
} // namespace

std::optional<TrimImage> isoTrimImage(const Handle(Geom_Surface)& surface,
                                      const Handle(Geom2d_Curve)& pcurve,
                                      double                      first,
                                      double                      last,
                                      double                      tolerance)
{
  const Geom2dAdaptor_Curve on(pcurve, first, last);
  if (on.GetType() != GeomAbs_Line)
    return std::nullopt;
  const auto line      = on.Line();
  const auto direction = line.Direction();
  if (direction.X() != 0 && direction.Y() != 0)
    return std::nullopt;
  const bool   uIso  = direction.X() == 0;
  const double fixed = uIso ? line.Location().X() : line.Location().Y();
  const auto   exact = coefficients(GeomAdaptor_Surface(surface), uIso, fixed);
  if (!exact)
    return std::nullopt;
  const double offset   = uIso ? line.Location().Y() : line.Location().X();
  const double slope    = uIso ? direction.Y() : direction.X();
  const auto   endpoint = [&](double t) {
    return I::add(I::point(offset), I::multiply(I::point(slope), I::point(t)));
  };
  const auto a = endpoint(first), b = endpoint(last);
  if (!I::valid(a) || !I::valid(b))
    return std::nullopt;
  const I::Range range{std::min(a.low, b.low), std::max(a.high, b.high)};
  auto           curve = uIso ? surface->UIso(fixed) : surface->VIso(fixed);
  while (const auto trimmed = Handle(Geom_TrimmedCurve)::DownCast(curve))
    curve = trimmed->BasisCurve();
  const GeomAdaptor_Curve spatial(curve);
  Coefficients            actual;
  double                  extent = 1, speed = 0;
  if (spatial.GetType() == GeomAbs_Circle)
  {
    const auto circle = spatial.Circle();
    actual = {vector(circle.Location().XYZ()),
              scale(vector(circle.XAxis().Direction().XYZ()), I::point(circle.Radius())),
              scale(vector(circle.YAxis().Direction().XYZ()), I::point(circle.Radius()))};
  }
  else if (spatial.GetType() == GeomAbs_Line)
  {
    const auto axis = spatial.Line();
    actual          = {vector(axis.Location().XYZ()),
                       vector(axis.Direction().XYZ()),
                       vector(gp_XYZ(0, 0, 0))};
    extent          = magnitude(range);
  }
  else
    return std::nullopt;
  for (int i = 0; i < 3; ++i)
    speed = I::add(I::point(speed),
                   I::add(I::point(magnitude(actual.b[i])), I::point(magnitude(actual.c[i]))))
              .high;
  auto error = I::point(0);
  for (int i = 0; i < 3; ++i)
    error =
      I::add(error,
             I::add(I::point(magnitude(I::subtract(exact->a[i], actual.a[i]))),
                    I::add(I::multiply(I::point(extent),
                                       I::point(magnitude(I::subtract(exact->b[i], actual.b[i])))),
                           I::point(magnitude(I::subtract(exact->c[i], actual.c[i]))))));
  // The expanded image interval can add at most an endpoint enclosure's width.
  // Arc length bounds this two-sided set error, for opening exclusions as well
  // as lower distance bounds. All carriers still use their original UV trims.
  const auto width = [](I::Range x) { return I::subtract(I::point(x.high), I::point(x.low)).high; };
  error = I::add(error, I::multiply(I::point(speed), I::point(std::max(width(a), width(b)))));
  if (!I::valid(error) || error.high > tolerance || !(range.low < range.high))
    return std::nullopt;
  return TrimImage{curve, range.low, range.high, error.high};
}
} // namespace BRepExt
