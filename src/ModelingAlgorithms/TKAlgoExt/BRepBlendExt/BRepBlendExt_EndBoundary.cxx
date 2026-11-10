// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepBlendExt_EndBoundary.hxx>

namespace BRepExt::blending::detail
{
namespace
{
using Vector = std::array<double, 5>;
using Matrix = std::array<Vector, 5>;

bool solve(Matrix a, Vector b, Vector& x)
{
  for (int k = 0; k < 5; ++k)
  {
    int p = k;
    for (int i = k + 1; i < 5; ++i)
      if (std::abs(a[i][k]) > std::abs(a[p][k]))
        p = i;
    if (!(std::abs(a[p][k]) > 0))
      return false;
    std::swap(a[p], a[k]);
    std::swap(b[p], b[k]);
    for (int i = k + 1; i < 5; ++i)
    {
      const double q = a[i][k] / a[k][k];
      for (int j = k + 1; j < 5; ++j)
        a[i][j] -= q * a[k][j];
      b[i] -= q * b[k];
    }
  }
  for (int i = 4; i >= 0; --i)
  {
    x[i] = b[i];
    for (int j = i + 1; j < 5; ++j)
      x[i] -= a[i][j] * x[j];
    x[i] /= a[i][i];
    if (!std::isfinite(x[i]))
      return false;
  }
  return true;
}

//! Raw differential geometry: no tangent extension beyond a finite carrier.
bool rawJet(const Side& side, double u, double v, Jet& j)
{
  side.surface->D2(u, v, j.p, j.du, j.dv, j.duu, j.dvv, j.duv);
  const gp_Vec cross  = j.du.Crossed(j.dv);
  const double length = cross.Magnitude();
  if (!(length > 0))
    return false;
  const gp_Vec n  = cross / length;
  const gp_Vec cu = j.duu.Crossed(j.dv) + j.du.Crossed(j.duv);
  const gp_Vec cv = j.duv.Crossed(j.dv) + j.du.Crossed(j.dvv);
  j.n             = n * side.outward;
  j.ndu           = (cu - n * cu.Dot(n)) * (side.outward / length);
  j.ndv           = (cv - n * cv.Dot(n)) * (side.outward / length);
  return true;
}

class Equations
{
public:
  const Rail&   rail;
  const gp_Pln& plane;
  double        station;
  double        tolerance;

  struct Value
  {
    EndBoundaryPoint p;
    Vector           f{};
    gp_Vec           offset[4];
    Vector           metric{};
  };

  bool value(const Vector& x, Value& out) const
  {
    Station s;
    s.t = station;
    Jet    j[2];
    gp_Pnt c[2];
    for (int i = 0; i < 2; ++i)
    {
      if (!rawJet(rail.s[i], x[2 * i], x[2 * i + 1], j[i]))
        return false;
      s.uv[i][0]            = x[2 * i];
      s.uv[i][1]            = x[2 * i + 1];
      s.foot[i]             = j[i].p;
      s.normal[i]           = j[i].n;
      const double r        = rail.sigma * rail.s[i].distance;
      c[i]                  = j[i].p.Translated(j[i].n * r);
      out.offset[2 * i]     = j[i].du + j[i].ndu * r;
      out.offset[2 * i + 1] = j[i].dv + j[i].ndv * r;
      out.metric[2 * i]     = j[i].du.Magnitude();
      out.metric[2 * i + 1] = j[i].dv.Magnitude();
      out.p.contact[i]      = j[i].p;
    }
    s.centre = gp_Pnt((c[0].XYZ() + c[1].XYZ()) / 2);
    s.along  = s.normal[0].Crossed(s.normal[1]);
    if (!(s.along.Magnitude() > 0))
      return false;
    s.along.Normalize();
    Section section;
    if (!sectionAt(rail, s, section))
      return false;
    const auto curve = sectionCurve(section);
    out.p.point      = curve->Value(x[4]);
    out.p.centre     = s.centre;
    out.p.parameters = x;
    const gp_Vec gap(c[1], c[0]);
    for (int i = 0; i < 3; ++i)
      out.f[i] = gap.Coord(i + 1);
    out.f[3]      = gp_Vec(plane.Location(), out.p.point).Dot(gp_Vec(plane.Axis().Direction()));
    out.metric[4] = rail.scale;
    out.p.ballResidual  = gap.Magnitude();
    out.p.planeResidual = std::abs(out.f[3]);
    return std::isfinite(out.p.ballResidual + out.p.planeResidual);
  }

  bool differential(const Vector& x, Value& v, Matrix& a, gp_Vec derivatives[5]) const
  {
    if (!value(x, v))
      return false;
    const gp_Vec normal(plane.Axis().Direction());
    for (int k = 0; k < 5; ++k)
    {
      if (!(v.metric[k] > 0))
        return false;
      for (int i = 0; i < 3; ++i)
        a[i][k] = k == 4 ? 0 : v.offset[k].Coord(i + 1) * (k < 2 ? 1 : -1);
      // Centred difference balances O(h^2) truncation with O(epsilon/h) rounding.
      // It guides continuation only; acceptance uses freshly evaluated equations.
      const double h = std::cbrt(std::numeric_limits<double>::epsilon()) * rail.scale / v.metric[k];
      Vector       left = x, right = x;
      left[k] -= h;
      right[k] += h;
      Value l, r;
      if (!value(left, l) || !value(right, r))
        return false;
      derivatives[k] = gp_Vec(l.p.point, r.p.point) / (2 * h);
      a[3][k]        = derivatives[k].Dot(normal);
    }
    return true;
  }

  bool correct(Vector& x, const Vector& anchor, const Vector& row) const
  {
    for (;;)
    {
      Value  v;
      Matrix a{};
      gp_Vec derivatives[5];
      if (!differential(x, v, a, derivatives))
        return false;
      double extra = 0;
      for (int k = 0; k < 5; ++k)
      {
        extra += row[k] * (x[k] - anchor[k]);
        a[4][k] = row[k];
      }
      const double error = v.p.ballResidual + v.p.planeResidual + std::abs(extra);
      if (error <= tolerance)
        return true;
      Vector b{}, delta{};
      for (int i = 0; i < 4; ++i)
        b[i] = -v.f[i];
      b[4] = -extra;
      for (int k = 0; k < 5; ++k)
        for (int i = 0; i < 5; ++i)
          a[i][k] /= v.metric[k];
      if (!solve(a, b, delta))
        return false;
      double fraction = 1;
      bool   improved = false;
      while (fraction > 0)
      {
        Vector trial = x;
        double move = 0, constraint = 0;
        for (int k = 0; k < 5; ++k)
        {
          trial[k] += fraction * delta[k] / v.metric[k];
          move += std::abs(trial[k] - x[k]);
          constraint += row[k] * (trial[k] - anchor[k]);
        }
        if (!(move > 0))
          break;
        Value next;
        if (value(trial, next)
            && next.p.ballResidual + next.p.planeResidual + std::abs(constraint) < error)
        {
          x        = trial;
          improved = true;
          break;
        }
        fraction /= 2;
      }
      if (!improved)
        return false;
    }
  }

  bool tangent(const Vector& x, Vector& rate, gp_Vec& direction, Value& v) const
  {
    Matrix a{};
    gp_Vec d[5];
    if (!differential(x, v, a, d))
      return false;
    // Try each coordinate chart. Maximising the free component of a scaled null
    // vector selects a nonsingular chart without privileging either contact.
    double best = -1;
    Vector chosen{};
    for (int chart = 0; chart < 5; ++chart)
    {
      Matrix b = a;
      Vector rhs{}, null{};
      rhs[4] = 1;
      for (int k = 0; k < 5; ++k)
      {
        for (int i = 0; i < 4; ++i)
          b[i][k] /= v.metric[k];
        b[4][k] = k == chart ? 1 : 0;
      }
      if (!solve(b, rhs, null))
        continue;
      double norm = 0;
      for (double z : null)
        norm += z * z;
      const double quality = 1 / std::sqrt(norm);
      if (quality > best)
      {
        best   = quality;
        chosen = null;
      }
    }
    if (!(best > 0))
      return false;
    direction = gp_Vec();
    for (int k = 0; k < 5; ++k)
    {
      rate[k] = chosen[k] / v.metric[k];
      direction += d[k] * rate[k];
    }
    const double speed = direction.Magnitude();
    if (!(speed > 0))
      return false;
    direction /= speed;
    for (double& z : rate)
      z /= speed;
    return true;
  }
};
} // namespace

EndBoundaryTrace tracePlaneEnd(const Rail&    rail,
                               const gp_Pln&  plane,
                               const Station& seed,
                               double         sectionParameter,
                               double         spacing,
                               const gp_Vec&  initialDirection)
{
  EndBoundaryTrace result;
  if (sectionParameter < 0 || sectionParameter > 1 || !(spacing > 0)
      || (sectionParameter > 0 && sectionParameter < 1 && !(initialDirection.Magnitude() > 0))
      || !rail.sizes.places.empty() || (rail.law != Law::Conic && rail.law != Law::G2)
      || rail.s[0].cliffCurve || rail.s[1].cliffCurve || rail.knotWindows)
    return result;
  Equations equations{rail, plane, seed.t, rail.resolved(limit::solved)};
  Vector    x{seed.uv[0][0], seed.uv[0][1], seed.uv[1][0], seed.uv[1][1], sectionParameter};
  Vector    row{};
  row[4]                  = rail.scale;
  const Vector seedAnchor = x;
  if (!equations.correct(x, seedAnchor, row))
    return result;
  Equations::Value initialValue;
  if (!equations.value(x, initialValue))
    return result;
  // A root on a polynomial extension is not a seed on this finite support.
  // Permit only the positional rounding already allowed by the event equations.
  for (int side = 0; side < 2; ++side)
  {
    const auto finite = runsOn(rail.s[side].surface);
    for (int k = 0; k < 2; ++k)
      if (finite[k])
      {
        const int    coordinate = 2 * side + k;
        const double outside =
          std::max(finite[k]->first - x[coordinate], x[coordinate] - finite[k]->second);
        if (outside * initialValue.metric[coordinate] > equations.tolerance)
          return result;
      }
  }
  const Vector initial = x;
  gp_Vec       previous;
  double       length = 0;
  for (;;)
  {
    Vector           rate{};
    gp_Vec           direction;
    Equations::Value v;
    if (!equations.tangent(x, rate, direction, v))
    {
      result.end = EndBoundaryTrace::End::Singular;
      return result;
    }
    const bool reverse =
      result.points.empty()
        ? (initialDirection.Magnitude() > 0 ? direction.Dot(initialDirection) < 0
                                            : (sectionParameter == 0 ? rate[4] < 0 : rate[4] > 0))
        : direction.Dot(previous) < 0;
    if (reverse)
    {
      for (double& z : rate)
        z = -z;
      direction.Reverse();
    }
    result.points.push_back(v.p);
    Vector next{};
    double step = spacing;
    for (;;)
    {
      Vector predicted = x;
      row              = {};
      for (int k = 0; k < 5; ++k)
      {
        predicted[k] += step * rate[k];
        row[k] = rate[k] * v.metric[k] * v.metric[k];
      }
      next                       = predicted;
      const bool       corrected = equations.correct(next, predicted, row);
      Equations::Value end;
      if (corrected && equations.value(next, end)
          && end.p.point.Distance(v.p.point.Translated(direction * step)) <= step / 8)
        break;
      step /= 2;
      if (!(step > equations.tolerance))
      {
        result.end = EndBoundaryTrace::End::Stalled;
        return result;
      }
    }
    // The first crossed boundary is an event equation, solved jointly with the four
    // geometric equations. Both contacts may label the same component's endpoints.
    int        boundary      = -1;
    double     boundaryValue = 0, fraction = 2;
    const auto crossing = [&](int k, double bound) {
      if ((x[k] - bound) * (next[k] - bound) >= 0 || x[k] == bound)
        return;
      const double f = (bound - x[k]) / (next[k] - x[k]);
      if (f < fraction)
      {
        fraction      = f;
        boundary      = k;
        boundaryValue = bound;
      }
    };
    crossing(4, 0);
    crossing(4, 1);
    for (int side = 0; side < 2; ++side)
    {
      const auto finite = runsOn(rail.s[side].surface);
      for (int k = 0; k < 2; ++k)
        if (finite[k])
        {
          crossing(2 * side + k, finite[k]->first);
          crossing(2 * side + k, finite[k]->second);
        }
    }
    if (boundary >= 0)
    {
      Vector anchor = next;
      for (int k = 0; k < 5; ++k)
        next[k] = x[k] + fraction * (next[k] - x[k]);
      anchor           = next;
      anchor[boundary] = boundaryValue;
      row              = {};
      row[boundary]    = v.metric[boundary];
      if (!equations.correct(next, anchor, row))
      {
        result.end = EndBoundaryTrace::End::Stalled;
        return result;
      }
      Equations::Value end;
      if (!equations.value(next, end))
      {
        result.end = EndBoundaryTrace::End::Stalled;
        return result;
      }
      result.points.push_back(end.p);
      result.boundaryCoordinate = boundary;
      const auto active         = [&](int k, double bound) {
        if (std::abs(next[k] - bound) * end.metric[k] <= equations.tolerance)
        {
          if (k == 4)
            result.contactSide = int(bound);
          if (std::find(result.boundaryCoordinates.begin(), result.boundaryCoordinates.end(), k)
              == result.boundaryCoordinates.end())
            result.boundaryCoordinates.push_back(k);
        }
      };
      active(4, 0);
      active(4, 1);
      for (int side = 0; side < 2; ++side)
      {
        const auto finite = runsOn(rail.s[side].surface);
        for (int k = 0; k < 2; ++k)
          if (finite[k])
          {
            active(2 * side + k, finite[k]->first);
            active(2 * side + k, finite[k]->second);
          }
      }
      const bool onContact =
        std::find(result.boundaryCoordinates.begin(), result.boundaryCoordinates.end(), 4)
        != result.boundaryCoordinates.end();
      // The corrected seed may already satisfy a finite-carrier restriction. Its
      // first crossing is then the same event, not a curve returning to a contact.
      // Compare parameter states in the same positional metric used by correction;
      // distinct charts at a spatial self-intersection are not merged.
      bool sameSeed = result.points.size() == 2;
      for (int k = 0; k < 5 && sameSeed; ++k)
        sameSeed = std::abs(next[k] - initial[k]) * end.metric[k] <= equations.tolerance;
      result.end = onContact && !sameSeed ? EndBoundaryTrace::End::Contact
                                          : EndBoundaryTrace::End::SurfaceBoundary;
      return result;
    }
    length += step;
    if (length > 2 * spacing)
    {
      double separation = 0;
      for (int k = 0; k < 5; ++k)
        separation += std::pow((next[k] - initial[k]) * v.metric[k], 2);
      if (std::sqrt(separation) <= equations.tolerance)
      {
        result.end = EndBoundaryTrace::End::Closed;
        return result;
      }
    }
    x        = next;
    previous = direction;
  }
}
} // namespace BRepExt::blending::detail
