// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepBlendExt_BoundaryCurve.hxx>
#include <GeomAlgoExt_Interval.hxx>

#include <functional>
#include <limits>
#include <map>

namespace BRepExt::blending::detail
{
namespace
{
namespace I               = GeomAlgoExt::interval;
using Polynomial          = std::vector<I::Range>;
constexpr double infinity = std::numeric_limits<double>::infinity();

// All coefficients, including the binomial ratios, enclose exact arithmetic on
// the stored doubles. No rounded conversion to Bezier surface patches is used.
Polynomial binomials(int n)
{
  Polynomial row(1, I::point(1));
  for (int k = 1; k <= n; ++k)
  {
    Polynomial next(k + 1, I::point(1));
    for (int j = 1; j < k; ++j)
      next[j] = I::add(row[j - 1], row[j]);
    row = std::move(next);
  }
  return row;
}

Polynomial product(const Polynomial& a, const Polynomial& b)
{
  const int  m = int(a.size()) - 1, n = int(b.size()) - 1;
  // These outward-rounded constants depend only on the two polynomial degrees.
  // Reusing them preserves the arithmetic enclosure and accumulation order.
  static thread_local std::map<std::pair<int, int>, Polynomial> factors;
  auto [entry, inserted] = factors.try_emplace({m, n});
  if (inserted)
  {
    const auto ca = binomials(m), cb = binomials(n), cc = binomials(m + n);
    for (int i = 0; i <= m; ++i)
      for (int j = 0; j <= n; ++j)
        entry->second.push_back(I::dividePositive(I::multiply(ca[i], cb[j]), cc[i + j]));
  }
  const auto& coefficient = entry->second;
  Polynomial out(m + n + 1, I::point(0));
  for (int i = 0; i <= m; ++i)
    for (int j = 0; j <= n; ++j)
      out[i + j] = I::add(out[i + j],
                          I::multiply(I::multiply(a[i], b[j]),
                                      coefficient[i * (n + 1) + j]));
  return out;
}

Polynomial scaled(Polynomial a, I::Range s)
{
  for (auto& x : a)
    x = I::multiply(x, s);
  return a;
}

void accumulate(Polynomial& a, const Polynomial& b)
{
  for (size_t i = 0; i < a.size(); ++i)
    a[i] = I::add(a[i], b[i]);
}

//! Cox-de Boor restricted to one knot cell, with its parameter a Bernstein
//! polynomial. The polynomial extension outside this cell is harmless: its
//! enclosure is used only for those t where the pcurve lies inside the cell.
std::vector<Polynomial> basis(const Polynomial&          u,
                              const std::vector<double>& knots,
                              int                        span,
                              int                        degree)
{
  std::vector<Polynomial> previous(1, Polynomial{I::point(1)});
  for (int d = 1; d <= degree; ++d)
  {
    std::vector<Polynomial> next(d + 1, Polynomial(3 * d + 1, I::point(0)));
    for (int j = 0; j <= d; ++j)
    {
      const int i = span - d + j;
      if (j > 0 && knots[i + d] > knots[i])
      {
        Polynomial factor = u;
        for (auto& x : factor)
          x = I::dividePositive(I::subtract(x, I::point(knots[i])),
                                I::subtract(I::point(knots[i + d]), I::point(knots[i])));
        accumulate(next[j], product(factor, previous[j - 1]));
      }
      if (j < d && knots[i + d + 1] > knots[i + 1])
      {
        Polynomial factor = u;
        for (auto& x : factor)
          x = I::dividePositive(I::subtract(I::point(knots[i + d + 1]), x),
                                I::subtract(I::point(knots[i + d + 1]), I::point(knots[i + 1])));
        accumulate(next[j], product(factor, previous[j]));
      }
    }
    previous = std::move(next);
  }
  return previous;
}

double normBound(const std::array<Polynomial, 3>& numerator, double denominator)
{
  if (!(denominator > 0))
    return infinity;
  I::Range sum = I::point(0);
  for (const auto& coordinate : numerator)
  {
    double high = 0;
    for (const auto x : coordinate)
    {
      if (!I::valid(x))
        return infinity;
      high = std::max({high, std::abs(x.low), std::abs(x.high)});
    }
    const auto quotient = I::dividePositive(I::point(high), I::point(denominator));
    sum                 = I::add(sum, I::multiply(quotient, quotient));
  }
  return I::valid(sum) ? I::roundUp(std::sqrt(sum.high)) : infinity;
}

//! The affine plane representation is bounded with the same stored-point
//! arithmetic as the spline representation, rather than sampled in space.
double planeGap(const gp_Ax3&                  frame,
                const std::array<gp_Pnt, 4>&   curve,
                const std::array<gp_Pnt2d, 4>& pc)
{
  std::array<Polynomial, 3> residual;
  for (int j = 0; j < 3; ++j)
    for (int i = 0; i < 4; ++i)
    {
      auto value =
        I::add(I::point(frame.Location().Coord(j + 1)),
               I::add(I::multiply(I::point(pc[i].X()), I::point(frame.XDirection().Coord(j + 1))),
                      I::multiply(I::point(pc[i].Y()), I::point(frame.YDirection().Coord(j + 1)))));
      residual[j].push_back(I::subtract(value, I::point(curve[i].Coord(j + 1))));
    }
  return normBound(residual, 1);
}
} // namespace

double boundaryCurveGap(const Handle(Geom_BSplineSurface)& surface,
                        const std::array<gp_Pnt, 4>&       curve,
                        const std::array<gp_Pnt2d, 4>&     pcurve)
{
  if (surface.IsNull() || surface->IsUPeriodic() || surface->IsVPeriodic())
    return infinity;
  std::vector<double> knots[2];
  std::vector<int>    cells[2];
  Polynomial          parameters[2];
  const int           degrees[2] = {surface->UDegree(), surface->VDegree()};
  double              u0, u1, v0, v1;
  surface->Bounds(u0, u1, v0, v1);
  const double bounds[2][2] = {{u0, u1}, {v0, v1}};
  for (int d = 0; d < 2; ++d)
  {
    const auto& values = d == 0 ? surface->UKnotSequence() : surface->VKnotSequence();
    for (int i = values.Lower(); i <= values.Upper(); ++i)
      knots[d].push_back(values(i));
    double lo = infinity, hi = -infinity;
    for (const auto& p : pcurve)
    {
      const double value = p.Coord(d + 1);
      if (!std::isfinite(value))
        return infinity;
      lo = std::min(lo, value);
      hi = std::max(hi, value);
      parameters[d].push_back(I::point(value));
    }
    if (lo < bounds[d][0] || hi > bounds[d][1])
      return infinity;
    for (int i = degrees[d]; i + degrees[d] + 1 < int(knots[d].size()); ++i)
      if (knots[d][i] < knots[d][i + 1] && knots[d][i] <= hi && knots[d][i + 1] >= lo)
        cells[d].push_back(i);
    if (cells[d].empty())
      return infinity;
  }
  const gp_Pnt origin = curve[0];
  double       answer = 0;
  for (const int u : cells[0])
  {
    const auto bu = basis(parameters[0], knots[0], u, degrees[0]);
    for (const int v : cells[1])
    {
      const auto                bv    = basis(parameters[1], knots[1], v, degrees[1]);
      const int                 count = 3 * (degrees[0] + degrees[1]) + 1;
      Polynomial                weight(count, I::point(0));
      std::array<Polynomial, 3> numerator;
      for (auto& coordinate : numerator)
        coordinate.assign(count, I::point(0));
      for (int i = 0; i <= degrees[0]; ++i)
        for (int j = 0; j <= degrees[1]; ++j)
        {
          const int  pi = u - degrees[0] + i + 1, pj = v - degrees[1] + j + 1;
          const auto w = scaled(product(bu[i], bv[j]), I::point(surface->Weight(pi, pj)));
          accumulate(weight, w);
          for (int k = 0; k < 3; ++k)
            accumulate(numerator[k],
                       scaled(w,
                              I::subtract(I::point(surface->Pole(pi, pj).Coord(k + 1)),
                                          I::point(origin.Coord(k + 1)))));
        }
      double low = infinity;
      for (const auto w : weight)
      {
        if (!I::valid(w))
          return infinity;
        low = std::min(low, w.low);
      }
      for (int k = 0; k < 3; ++k)
      {
        Polynomial coordinate;
        for (const auto& p : curve)
          coordinate.push_back(
            I::subtract(I::point(p.Coord(k + 1)), I::point(origin.Coord(k + 1))));
        // Multiplication by the constant one in degree three elevates exactly.
        numerator[k]           = product(numerator[k], Polynomial(4, I::point(1)));
        const auto represented = product(weight, coordinate);
        for (size_t i = 0; i < represented.size(); ++i)
          numerator[k][i] = I::subtract(numerator[k][i], represented[i]);
      }
      answer = std::max(answer, normBound(numerator, low));
    }
  }
  return answer;
}

std::optional<BoundaryCurve> planeBoundaryCurve(const Sheet&                 sheet,
                                                const TopoDS_Face&           face,
                                                const std::vector<gp_Pnt2d>& guide,
                                                const TopoDS_Vertex&         first,
                                                const TopoDS_Vertex&         last)
{
  const auto   spline = Handle(Geom_BSplineSurface)::DownCast(sheet.surface);
  const Placed one = placedSheet(sheet), two = placedFace(face);
  const auto   plane = Handle(Geom_Plane)::DownCast(two.surface);
  if (spline.IsNull() || plane.IsNull() || guide.size() < 2)
    return {};
  const gp_Pln pln = plane->Pln();
  const gp_Vec normal(pln.Axis().Direction());

  struct Point
  {
    gp_Pnt2d uv;
    gp_Pnt   p;
    gp_Vec2d rate;
    gp_Vec   tangent;
  };

  const auto point = [&](gp_Pnt2d uv, const gp_Vec2d& along) -> std::optional<Point> {
    const gp_Pnt2d seed       = uv;
    gp_Vec2d       transverse = along;
    double         u0, u1, v0, v1;
    spline->Bounds(u0, u1, v0, v1);
    if (uv.X() == u0 || uv.X() == u1)
      transverse = gp_Vec2d(1, 0);
    else if (uv.Y() == v0 || uv.Y() == v1)
      transverse = gp_Vec2d(0, 1);
    for (;;)
    {
      gp_Pnt p;
      gp_Vec du, dv;
      spline->D1(uv.X(), uv.Y(), p, du, dv);
      const double f = gp_Vec(pln.Location(), p).Dot(normal), a = du.Dot(normal),
                   b = dv.Dot(normal), g = gp_Vec2d(seed, uv).Dot(transverse),
                   determinant = a * transverse.Y() - b * transverse.X();
      if (!std::isfinite(determinant) || determinant == 0)
        return {};
      gp_Vec2d rate(b, -a);
      gp_Vec   tangent = du * rate.X() + dv * rate.Y();
      if (!(tangent.Magnitude() > 0))
        return {};
      if (rate.Dot(along) < 0)
      {
        rate.Reverse();
        tangent.Reverse();
      }
      const double speed = tangent.Magnitude();
      if (std::abs(f) <= std::max(limit::marched, tol::rounding(p.XYZ().Modulus()))
          && std::abs(g) <= tol::rounding(seed.XY().Modulus() * transverse.Magnitude()))
        return Point{uv, p, rate / speed, tangent / speed};
      const gp_Vec2d change((f * transverse.Y() - b * g) / determinant,
                            (a * g - f * transverse.X()) / determinant);
      const gp_Pnt2d next = uv.Translated(-change);
      if (next.IsEqual(uv, 0) || !std::isfinite(next.X()) || !std::isfinite(next.Y()))
        return {};
      // Newton must reduce the two defining equations. A failed correction is
      // unresolved; it cannot be turned into a geometric obstruction.
      const gp_Pnt q = spline->Value(next.X(), next.Y());
      if (std::abs(gp_Vec(pln.Location(), q).Dot(normal)) >= std::abs(f)
          && std::abs(f) > std::max(limit::marched, tol::rounding(p.XYZ().Modulus())))
        return {};
      uv = next;
    }
  };

  struct Span
  {
    std::array<gp_Pnt, 4>   p;
    std::array<gp_Pnt2d, 4> uv, plane;
    double                  length;
  };

  std::vector<Span>                                     spans;
  BoundaryCurve                                         result;
  const std::function<bool(const Point&, const Point&)> refine = [&](const Point& a,
                                                                     const Point& b) {
    const double length = a.p.Distance(b.p);
    if (!(length > 0))
      return false;
    Span span{
      {a.p,
       a.p.Translated(a.tangent * (length / 3)),
       b.p.Translated(b.tangent * (-length / 3)),
       b.p},
      {a.uv, a.uv.Translated(a.rate * (length / 3)), b.uv.Translated(b.rate * (-length / 3)), b.uv},
      {},
      length};
    for (int i = 0; i < 4; ++i)
    {
      double u, v;
      ElSLib::Parameters(pln, span.p[i], u, v);
      span.plane[i] = gp_Pnt2d(u, v);
    }
    const double gap  = boundaryCurveGap(spline, span.p, span.uv),
                 flat = planeGap(pln.Position(), span.p, span.plane);
    if (gap <= tol::length && flat <= tol::length)
    {
      spans.push_back(span);
      result.sheetGap = std::max(result.sheetGap, gap);
      result.planeGap = std::max(result.planeGap, flat);
      return true;
    }
    const gp_Pnt2d middle((a.uv.XY() + b.uv.XY()) / 2);
    if (middle.IsEqual(a.uv, 0) || middle.IsEqual(b.uv, 0))
      return false;
    const auto m = point(middle, gp_Vec2d(a.uv, b.uv));
    return m && refine(a, *m) && refine(*m, b);
  };
  std::vector<Point> points;
  for (size_t i = 0; i < guide.size(); ++i)
  {
    const gp_Vec2d along(guide[i == 0 ? i : i - 1], guide[i + 1 == guide.size() ? i : i + 1]);
    const auto     p = point(guide[i], along);
    if (!p)
      return {};
    points.push_back(*p);
  }
  if (points.front().p.Distance(BRep_Tool::Pnt(first)) > tol::length
      || points.back().p.Distance(BRep_Tool::Pnt(last)) > tol::length)
    return {};
  for (size_t i = 1; i < points.size(); ++i)
    if (!refine(points[i - 1], points[i]))
      return {};
  const int                    n = int(spans.size());
  NCollection_Array1<gp_Pnt>   poles(1, 3 * n + 1);
  NCollection_Array1<gp_Pnt2d> uv(1, 3 * n + 1), planar(1, 3 * n + 1);
  NCollection_Array1<double>   knots(1, n + 1);
  NCollection_Array1<int>      multiplicities(1, n + 1);
  double                       t = 0;
  for (int i = 0; i < n; ++i)
  {
    knots(i + 1)          = t;
    multiplicities(i + 1) = i == 0 ? 4 : 3;
    for (int j = 0; j < 4; ++j)
    {
      poles(3 * i + j + 1)  = spans[i].p[j];
      uv(3 * i + j + 1)     = spans[i].uv[j];
      planar(3 * i + j + 1) = spans[i].plane[j];
    }
    const double next = t + spans[i].length;
    if (!(next > t))
      return {};
    t = next;
  }
  knots(n + 1)                            = t;
  multiplicities(n + 1)                   = 4;
  const Handle(Geom_BSplineCurve)   curve = new Geom_BSplineCurve(poles, knots, multiplicities, 3);
  const Handle(Geom2d_BSplineCurve) pc    = new Geom2d_BSplineCurve(uv, knots, multiplicities, 3),
                                    pp = new Geom2d_BSplineCurve(planar, knots, multiplicities, 3);
  result.edge  = edgeBetween(curve,
                             first,
                             0,
                             last,
                             t,
                             {{Handle(Geom2d_Curve)(pc), one.basis, one.location},
                              {Handle(Geom2d_Curve)(pp), two.basis, two.location}},
                             tol::length);
  result.spans = n;
  return result;
}
} // namespace BRepExt::blending::detail
