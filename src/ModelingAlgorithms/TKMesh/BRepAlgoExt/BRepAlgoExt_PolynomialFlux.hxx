// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#ifndef _BRepAlgoExt_PolynomialFlux_HeaderFile
#define _BRepAlgoExt_PolynomialFlux_HeaderFile

#include <BRepAdaptor_Surface.hxx>
#include <BRepAdaptor_Curve2d.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_BezierCurve.hxx>
#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <vector>
#include <type_traits>

namespace BRepExt::polynomialFlux
{
// Keep the stored binary64 geometry exact on input. Extra precision is needed
// through evaluation and quadrature, not just in the final contour sum.
// Unevaluated sum of two binary64 numbers. Error-free addition (TwoSum)
// and multiplication (Dekker splitting) preserve the low part without FMA.
// The geometry and final result remain binary64 on every platform.
struct Real
{
  double high = 0, low = 0;
  Real() = default;

  Real(double x)
      : high(x)
  {
  }

  Real(double x, double y)
      : high(x),
        low(y)
  {
  }

  explicit operator double() const { return high + low; }

  static Real sum(double a, double b)
  {
    const double s = a + b, v = s - a;
    return {s, (a - (s - v)) + (b - v)};
  }

  // Renormalize a canonical expansion after its low-order corrections.
  static Real quickSum(double a, double b)
  {
    const double s = a + b;
    return {s, b - (s - a)};
  }

  static Real product(double a, double b)
  {
    constexpr double split    = 134217729.; // 2^27+1 splits 53-bit significands.
    constexpr double largest  = std::numeric_limits<double>::max() / split;
    constexpr double smallest = std::numeric_limits<double>::min() * split;
    const double     aa = std::abs(a), ab = std::abs(b);
    if (aa <= largest && ab <= largest && (aa >= smallest || a == 0) && (ab >= smallest || b == 0))
    {
      const double ca = split * a, cb = split * b;
      const double ah = ca - (ca - a), al = a - ah, bh = cb - (cb - b), bl = b - bh;
      const double p = a * b;
      return {p, ((ah * bh - p) + ah * bl + al * bh) + al * bl};
    }
    // Scale only when splitting could overflow or underflow.
    int          ea, eb;
    const double ma = std::frexp(a, &ea), mb = std::frexp(b, &eb);
    const double ca = 134217729. * ma, cb = 134217729. * mb;
    const double ah = ca - (ca - ma), al = ma - ah, bh = cb - (cb - mb), bl = mb - bh;
    const double p = ma * mb;
    return {a * b, std::ldexp(((ah * bh - p) + ah * bl + al * bh) + al * bl, ea + eb)};
  }

  friend Real operator+(Real a, Real b)
  {
    const auto s = sum(a.high, b.high), t = sum(a.low, b.low);
    const auto u = quickSum(s.high, s.low + t.high);
    return quickSum(u.high, u.low + t.low);
  }

  friend Real operator-(Real a) { return {-a.high, -a.low}; }

  friend Real operator-(Real a, Real b) { return a + -b; }

  friend Real operator*(Real a, Real b)
  {
    const auto p = product(a.high, b.high);
    return quickSum(p.high, p.low + a.high * b.low + a.low * b.high + a.low * b.low);
  }

  friend Real operator/(Real a, Real b)
  {
    const Real q(a.high / b.high);
    const Real r = a - b * q;
    return quickSum(q.high, (r.high + r.low) / b.high);
  }

  friend bool operator==(Real a, Real b) { return a.high == b.high && a.low == b.low; }

  friend bool operator<(Real a, Real b)
  {
    return a.high < b.high || (a.high == b.high && a.low < b.low);
  }

  friend bool operator>(Real a, Real b) { return b < a; }

  friend bool operator<=(Real a, Real b) { return !(b < a); }

  friend bool operator>=(Real a, Real b) { return !(a < b); }

  Real& operator+=(Real x) { return *this = *this + x; }

  Real& operator-=(Real x) { return *this = *this - x; }

  Real& operator*=(Real x) { return *this = *this * x; }

  Real& operator/=(Real x) { return *this = *this / x; }
};

inline Real absolute(Real x)
{
  return x < Real(0) ? -x : x;
}

using Vector = std::array<Real, 3>;

inline Vector difference(const Vector& a, const Vector& b)
{
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

using Rule = std::vector<std::pair<Real, Real>>;

inline const Rule& rule(int order)
{
  static thread_local std::map<int, Rule> rules;
  auto                                    it = rules.find(order);
  if (it == rules.end())
  {
    Rule rule;
    for (int i = 0; i < (order + 1) / 2; ++i)
    {
      Real x = std::cos(std::acos(-1.) * (i + .75) / (order + .5)), dp = 0;
      for (;;)
      {
        Real p = 1, previous = 0;
        for (int j = 1; j <= order; ++j)
        {
          const Real q = p;
          p            = ((2 * j - 1) * x * p - (j - 1) * previous) / j;
          previous     = q;
        }
        dp              = order * (x * p - previous) / (x * x - 1);
        const Real step = p / dp, next = x - step;
        if (next == x || absolute(step) <= Real(0x1p-104) * absolute(x))
          break;
        x = next;
      }
      const Real w = 2 / ((1 - x * x) * dp * dp);
      rule.emplace_back(x, w);
      if (2 * i + 1 < order)
        rule.emplace_back(-x, w);
    }
    it = rules.emplace(order, std::move(rule)).first;
  }
  return it->second;
}

// Legendre projection factors depend only on the quadrature order.
inline const std::vector<std::vector<Real>>& projection(int order)
{
  static thread_local std::map<int, std::vector<std::vector<Real>>> cache;
  auto                                                              it = cache.find(order);
  if (it == cache.end())
  {
    // Convert the finite Legendre projection to Chebyshev coefficients.
    // x*T_0=T_1; x*T_j=(T_{j-1}+T_{j+1})/2 for j>0.
    std::vector<std::vector<Real>> basis(order, std::vector<Real>(order));
    basis[0][0] = 1;
    if (order > 1)
      basis[1][1] = 1;
    for (int k = 1; k + 1 < order; ++k)
    {
      const Real alpha = Real(2 * k + 1) / (k + 1), beta = Real(k) / (k + 1);
      for (int j = 0; j <= k; ++j)
      {
        if (j == 0)
          basis[k + 1][1] += alpha * basis[k][0];
        else
        {
          basis[k + 1][j - 1] += alpha * basis[k][j] / 2;
          basis[k + 1][j + 1] += alpha * basis[k][j] / 2;
        }
        basis[k + 1][j] -= beta * basis[k - 1][j];
      }
    }
    std::vector<std::vector<Real>> factors;
    for (const auto& [x, w] : rule(order))
    {
      std::vector<Real> row;
      Real              p = 1, previous = 0;
      for (int k = 0; k < order; ++k)
      {
        row.push_back(w * p * (2 * k + 1) / 2);
        const Real next = ((2 * k + 1) * x * p - k * previous) / (k + 1);
        previous        = p;
        p               = next;
      }
      std::vector<Real> chebyshev(order);
      for (int k = 0; k < order; ++k)
        for (int j = 0; j <= k; ++j)
          chebyshev[j] += row[k] * basis[k][j];
      factors.push_back(std::move(chebyshev));
    }
    it = cache.emplace(order, std::move(factors)).first;
  }
  return it->second;
}

struct Curve
{
  int                       degree    = 1;
  int                       dimension = 3;
  bool                      periodic  = false;
  std::vector<Vector>       poles, derivatives;
  std::vector<Real>         knots;
  mutable std::vector<Real> inverseDifferences;

  template <class C>
  void read(const C& c)
  {
    degree = c->Degree();
    for (int i = 1; i <= c->NbPoles(); ++i)
    {
      const auto p = c->Pole(i);
      Vector     v{p.X(), p.Y(), 0};
      if constexpr (std::is_same_v<decltype(p), const gp_Pnt>)
        v[2] = p.Z();
      poles.push_back(v);
    }
  }

  template <class C>
  void spline(const C& c)
  {
    read(c);
    periodic = c->IsPeriodic();
    for (double k : c->KnotSequence())
      knots.push_back(k);
    const size_t n = poles.size(), needed = knots.size() - degree - 1;
    for (size_t i = n; i < needed; ++i)
      poles.push_back(poles[(i - n) % n]);
    inverseDifferences.resize(knots.size() * degree);
  }

  template <class C>
  void bezier(const C& c)
  {
    read(c);
    knots.assign(degree + 1, 0);
    knots.insert(knots.end(), degree + 1, 1);
    inverseDifferences.resize(knots.size() * degree);
  }

  void differentiate()
  {
    for (size_t i = 0; i + 1 < poles.size(); ++i)
    {
      auto       d   = difference(poles[i + 1], poles[i]);
      const Real den = knots[i + degree + 1] - knots[i + 1];
      for (auto& x : d)
        x = den == 0 ? 0 : degree * x / den;
      derivatives.push_back(d);
    }
  }

  explicit Curve(const Handle(Adaptor3d_Curve)& c)
  {
    if (c->GetType() == GeomAbs_BSplineCurve)
      spline(c->BSpline());
    else if (c->GetType() == GeomAbs_BezierCurve)
      bezier(c->Bezier());
    else
    {
      const auto l = c->Line();
      poles        = {{l.Location().X(), l.Location().Y(), l.Location().Z()},
                      {l.Direction().X(), l.Direction().Y(), l.Direction().Z()}};
    }
  }

  explicit Curve(const BRepAdaptor_Curve2d& c)
      : dimension(2)
  {
    if (c.GetType() == GeomAbs_BSplineCurve)
      spline(c.BSpline());
    else if (c.GetType() == GeomAbs_BezierCurve)
      bezier(c.Bezier());
    else
    {
      const auto l = c.Line();
      poles = {{l.Location().X(), l.Location().Y(), 0}, {l.Direction().X(), l.Direction().Y(), 0}};
    }
  }

  Vector direct(Real t, bool derivative = false) const
  {
    if (knots.empty())
    {
      if (derivative)
        return poles[1];
      return {poles[0][0] + t * poles[1][0],
              poles[0][1] + t * poles[1][1],
              poles[0][2] + t * poles[1][2]};
    }
    if (periodic)
    {
      const Real first = knots[degree], period = knots[poles.size()] - first;
      t -= Real(std::floor(double((t - first) / period))) * period;
    }
    const int   offset = derivative ? 1 : 0, p = degree - offset;
    const auto& points = derivative ? derivatives : poles;
    const auto  begin = knots.begin() + offset, end = knots.end() - offset;
    const int   span =
      std::max(p,
               std::min(int(points.size()) - 1, int(std::upper_bound(begin, end, t) - begin) - 1));
    std::vector<Vector> d(points.begin() + span - p, points.begin() + span + 1);
    for (int r = 1; r <= p; ++r)
      for (int j = p; j >= r; --j)
      {
        const int i       = span - p + j + offset;
        auto&     inverse = inverseDifferences[i * degree + p - r];
        if (inverse == 0)
          inverse = 1 / (knots[i + p - r + 1] - knots[i]);
        const Real a = (t - knots[i]) * inverse;
        for (int k = 0; k < dimension; ++k)
          d[j][k] = (1 - a) * d[j - 1][k] + a * d[j][k];
      }
    return d.back();
  }

  struct Span
  {
    Real                first, last;
    std::vector<Vector> coefficients, derivativeCoefficients;
  };

  std::vector<Span> spans;

  void prepareTrim(Real first, Real last)
  {
    if (knots.empty())
      return;
    std::vector<Real> cuts{first, last};
    const Real        period = knots[poles.size()] - knots[degree];
    for (Real k : knots)
    {
      if (periodic)
      {
        const int begin = int(std::floor(double((first - k) / period))) + 1;
        const int end   = int(std::ceil(double((last - k) / period)));
        for (int winding = begin; winding < end; ++winding)
          cuts.push_back(k + winding * period);
      }
      else if (k > first && k < last)
        cuts.push_back(k);
    }
    std::sort(cuts.begin(), cuts.end());
    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
    for (size_t i = 1; i < cuts.size(); ++i)
    {
      const int   n = degree + 1;
      Span        span{cuts[i - 1], cuts[i], std::vector<Vector>(n), std::vector<Vector>(degree)};
      const Real  mid = (span.first + span.last) / 2, half = (span.last - span.first) / 2;
      const auto& factors = projection(n);
      size_t      node    = 0;
      for (const auto& [x, w] : rule(n))
      {
        const auto v = direct(mid + half * x);
        for (int k = 0; k < n; ++k)
          for (int j = 0; j < 2; ++j)
            span.coefficients[k][j] += factors[node][k] * v[j];
        ++node;
      }
      // T'_j=j*U_{j-1}; the backward Chebyshev derivative recurrence.
      for (int k = degree - 1; k >= 0; --k)
        for (int c = 0; c < 2; ++c)
          span.derivativeCoefficients[k][c] =
            2 * (k + 1) * span.coefficients[k + 1][c]
            + (k + 2 < degree ? span.derivativeCoefficients[k + 2][c] : Real(0));
      for (int c = 0; c < 2; ++c)
        span.derivativeCoefficients[0][c] /= 2;
      for (auto& coefficient : span.derivativeCoefficients)
        for (int c = 0; c < 2; ++c)
          coefficient[c] /= half;
      spans.push_back(std::move(span));
    }
  }

  Vector evaluate(Real t, bool derivative = false) const
  {
    if (spans.empty())
      return direct(t, derivative);
    const auto  it   = std::upper_bound(spans.begin(), spans.end(), t, [](Real x, const Span& s) {
      return x < s.last;
    });
    const auto& span = it == spans.end() ? spans.back() : *it;
    const Real  x    = (2 * t - span.first - span.last) / (span.last - span.first);
    const auto& coefficients = derivative ? span.derivativeCoefficients : span.coefficients;
    Vector      value{}, next{};
    for (int k = int(coefficients.size()) - 1; k >= 0; --k)
    {
      const Real factor = (k == 0 ? x : 2 * x);
      for (int j = 0; j < 2; ++j)
      {
        const Real previous = value[j];
        value[j]            = coefficients[k][j] + factor * value[j] - next[j];
        next[j]             = previous;
      }
    }
    return value;
  }

  // The contour needs (u,v,u'), never v'. Evaluate that jet together.
  Vector trimJet(Real t) const
  {
    if (spans.empty())
    {
      const auto uv = direct(t);
      return {uv[0], uv[1], direct(t, true)[0]};
    }
    const auto  it   = std::upper_bound(spans.begin(), spans.end(), t, [](Real x, const Span& s) {
      return x < s.last;
    });
    const auto& span = it == spans.end() ? spans.back() : *it;
    const Real  x    = (2 * t - span.first - span.last) / (span.last - span.first);
    Vector      value{}, next{};
    for (int k = int(span.coefficients.size()) - 1; k >= 0; --k)
    {
      const Real factor = (k == 0 ? x : 2 * x);
      for (int j = 0; j < 3; ++j)
      {
        const Real coefficient = j < 2 ? span.coefficients[k][j]
                                 : k < int(span.derivativeCoefficients.size())
                                   ? span.derivativeCoefficients[k][0]
                                   : Real(0);
        const Real previous    = value[j];
        value[j]               = coefficient + factor * value[j] - next[j];
        next[j]                = previous;
      }
    }
    return value;
  }
};

struct Flux
{
  Curve  carrier;
  Vector direction;

  Flux(const BRepAdaptor_Surface& s, const gp_Pnt& origin)
      : carrier(s.AdaptorSurfaceOriginal().BasisCurve())
  {
    const auto d = s.AdaptorSurfaceOriginal().Direction();
    for (int k = 0; k < 3; ++k)
    {
      direction[k] = 0;
      for (int j = 0; j < 3; ++j)
        direction[k] += Real(s.Trsf().Value(k + 1, j + 1)) * d.Coord(j + 1);
    }
    for (auto& p : carrier.poles)
    {
      const auto old = p;
      for (int k = 0; k < 3; ++k)
      {
        p[k] = Real(s.Trsf().Value(k + 1, 4)) - origin.Coord(k + 1);
        for (int j = 0; j < 3; ++j)
          p[k] += Real(s.Trsf().Value(k + 1, j + 1)) * old[j];
      }
    }
    // Differentiate after placement and origin subtraction, retaining their precision.
    if (!carrier.knots.empty())
    {
      carrier.derivatives.clear();
      carrier.differentiate();
    }
    else
      for (int k = 0; k < 3; ++k)
      {
        carrier.poles[1][k] = 0;
        const auto line     = s.AdaptorSurfaceOriginal().BasisCurve()->Line();
        for (int j = 0; j < 3; ++j)
          carrier.poles[1][k] += Real(s.Trsf().Value(k + 1, j + 1)) * line.Direction().Coord(j + 1);
      }
  }

  Real direct(Real u) const
  {
    const auto p = carrier.evaluate(u), d = carrier.evaluate(u, true);
    Real       value = 0;
    for (int k = 0; k < 3; ++k)
      value +=
        p[k] * (d[(k + 1) % 3] * direction[(k + 2) % 3] - d[(k + 2) % 3] * direction[(k + 1) % 3]);
    return value / 3;
  }

  struct Span
  {
    Real              first, last;
    std::vector<Real> coefficients;
  };

  std::vector<Span> spans;

  // A degree-d flux has a finite Legendre expansion. Its coefficient
  // integrands have degree <=2d, so d+1 Gauss nodes compute them exactly
  // in real arithmetic; change basis to Chebyshev for Clenshaw evaluation.
  // Cache the scalar, rather than reevaluating both
  // three-dimensional carrier jets at every contour node.
  void prepare(const std::vector<double>& breaks)
  {
    const int n = 2 * carrier.degree - 1;
    for (size_t i = 1; i < breaks.size(); ++i)
    {
      Span        span{breaks[i - 1], breaks[i], std::vector<Real>(n)};
      const Real  mid = (span.first + span.last) / 2, half = (span.last - span.first) / 2;
      const auto& factors = projection(n);
      size_t      node    = 0;
      for (const auto& [x, w] : rule(n))
      {
        const Real value = direct(mid + half * x);
        for (int k = 0; k < n; ++k)
          span.coefficients[k] += factors[node][k] * value;
        ++node;
      }
      spans.push_back(std::move(span));
    }
  }

  Real operator()(Real u) const
  {
    const auto  it    = std::upper_bound(spans.begin(), spans.end(), u, [](Real x, const Span& s) {
      return x < s.last;
    });
    const auto& span  = it == spans.end() ? spans.back() : *it;
    const Real  x     = (2 * u - span.first - span.last) / (span.last - span.first);
    Real        value = 0, next = 0;
    for (int k = int(span.coefficients.size()) - 1; k >= 0; --k)
    {
      const Real factor   = (k == 0 ? x : 2 * x);
      const Real previous = value;
      value               = span.coefficients[k] + factor * value - next;
      next                = previous;
    }
    return value;
  }
};

template <class F>
Real integral(F&& f, Real a, Real b, int order)
{
  const Real mid = (a + b) / 2, half = (b - a) / 2;
  Real       total = 0;
  for (const auto& [x, w] : rule(order))
    total += w * f(mid + half * x);
  return half * total;
}
} // namespace BRepExt::polynomialFlux
#endif
