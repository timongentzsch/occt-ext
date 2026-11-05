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
#include <BRepOffsetExt_Focal.hxx>
#include <GeomAlgoExt_Cancel.hxx>
#include <GeomAlgoExt_Interval.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_BezierSurface.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Precision.hxx>

#include <array>
#include <deque>
#include <vector>

namespace BRepExt
{
namespace
{
using namespace GeomAlgoExt::interval;
using Vec = std::array<Range, 3>;

Range scaled(Range a, double b)
{
  return multiply(a, point(b));
}

Range square(Range a)
{
  if (!valid(a))
    return invalid();
  const Range r = multiply(a, a);
  return {a.low <= 0 && a.high >= 0 ? 0.0 : std::max(0.0, r.low), r.high};
}

Range root(Range a)
{
  if (!valid(a) || a.high < 0)
    return invalid();
  return {a.low <= 0 ? 0.0 : roundDown(std::sqrt(a.low)), roundUp(std::sqrt(a.high))};
}

Range dot(const Vec& a, const Vec& b)
{
  Range r = point(0);
  for (int k = 0; k < 3; ++k)
    r = add(r, multiply(a[k], b[k]));
  return r;
}

Range norm2(const Vec& a)
{
  Range r = point(0);
  for (Range x : a)
    r = add(r, square(x));
  return r;
}

Vec cross(const Vec& a, const Vec& b)
{
  Vec r;
  for (int k = 0; k < 3; ++k)
    r[k] =
      subtract(multiply(a[(k + 1) % 3], b[(k + 2) % 3]), multiply(a[(k + 2) % 3], b[(k + 1) % 3]));
  return r;
}

// Cox-de Boor on one open knot span. All arithmetic includes rounding,
// including knot differences. No rounded knot insertion precedes the enclosure.
std::vector<Range> basis(const std::vector<double>& knots, int degree, int span, Range t)
{
  std::vector<Range> previous(knots.size() - 1, point(0));
  previous[span] = point(1);
  for (int p = 1; p <= degree; ++p)
  {
    std::vector<Range> next(knots.size() - p - 1, point(0));
    for (int i = std::max(0, span - p); i <= span && i < int(next.size()); ++i)
    {
      if (knots[i + p] > knots[i])
        next[i] = add(next[i],
                      multiply(previous[i],
                               dividePositive(subtract(t, point(knots[i])),
                                              subtract(point(knots[i + p]), point(knots[i])))));
      if (knots[i + p + 1] > knots[i + 1])
        next[i] =
          add(next[i],
              multiply(previous[i + 1],
                       dividePositive(subtract(point(knots[i + p + 1]), t),
                                      subtract(point(knots[i + p + 1]), point(knots[i + 1])))));
      next[i] = intersect(next[i], {0, 1});
    }
    previous = std::move(next);
  }
  return previous;
}

struct Net
{
  int                               p, q, nu, nv;
  std::vector<double>               u, v;
  std::vector<std::array<Range, 4>> poles;
};

bool readNet(Handle(Geom_Surface) surface, const gp_Trsf& location, Net& net)
{
  while (const auto trim = Handle(Geom_RectangularTrimmedSurface)::DownCast(surface))
    surface = trim->BasisSurface();
  const auto spline = Handle(Geom_BSplineSurface)::DownCast(surface);
  const auto bezier = Handle(Geom_BezierSurface)::DownCast(surface);
  if (spline.IsNull() && bezier.IsNull())
    return false;
  if (!spline.IsNull() && (spline->IsUPeriodic() || spline->IsVPeriodic()))
    return false;
  net.p  = spline.IsNull() ? bezier->UDegree() : spline->UDegree();
  net.q  = spline.IsNull() ? bezier->VDegree() : spline->VDegree();
  net.nu = spline.IsNull() ? bezier->NbUPoles() : spline->NbUPoles();
  net.nv = spline.IsNull() ? bezier->NbVPoles() : spline->NbVPoles();
  if (!spline.IsNull())
  {
    for (double k : spline->UKnotSequence())
      net.u.push_back(k);
    for (double k : spline->VKnotSequence())
      net.v.push_back(k);
  }
  else
  {
    net.u.insert(net.u.end(), net.p + 1, 0);
    net.u.insert(net.u.end(), net.p + 1, 1);
    net.v.insert(net.v.end(), net.q + 1, 0);
    net.v.insert(net.v.end(), net.q + 1, 1);
  }
  const gp_Pnt origin = spline.IsNull() ? bezier->Pole(1, 1) : spline->Pole(1, 1);
  for (int i = 1; i <= net.nu; ++i)
    for (int j = 1; j <= net.nv; ++j)
    {
      const gp_Pnt pole = spline.IsNull() ? bezier->Pole(i, j) : spline->Pole(i, j);
      const double w    = spline.IsNull() ? bezier->Weight(i, j) : spline->Weight(i, j);
      if (!(w > 0) || !std::isfinite(w))
        return false;
      std::array<Range, 4> h;
      for (int k = 0; k < 3; ++k)
      {
        Range coordinate = point(0);
        for (int l = 0; l < 3; ++l)
          coordinate =
            add(coordinate,
                multiply(point(location.HVectorialPart().Value(k + 1, l + 1)),
                         subtract(point(pole.Coord(l + 1)), point(origin.Coord(l + 1)))));
        h[k] = multiply(scaled(coordinate, location.ScaleFactor()), point(w));
      }
      h[3] = point(w);
      net.poles.push_back(h);
    }
  return true;
}

struct Cell
{
  double u0, u1, v0, v1;
  int    uspan, vspan, depth = 0;
};

Net derivative(const Net& source, bool alongU)
{
  Net         out    = source;
  const auto& knots  = alongU ? source.u : source.v;
  const int   degree = alongU ? source.p : source.q;
  if (degree == 0)
  {
    for (auto& pole : out.poles)
      pole.fill(point(0));
    return out;
  }
  if (alongU)
  {
    --out.p;
    --out.nu;
    out.u.erase(out.u.begin());
    out.u.pop_back();
  }
  else
  {
    --out.q;
    --out.nv;
    out.v.erase(out.v.begin());
    out.v.pop_back();
  }
  out.poles.clear();
  for (int i = 0; i < out.nu; ++i)
    for (int j = 0; j < out.nv; ++j)
    {
      const int   index = alongU ? i : j;
      const Range width = subtract(point(knots[index + degree + 1]), point(knots[index + 1]));
      std::array<Range, 4> pole;
      for (int k = 0; k < 4; ++k)
        pole[k] =
          knots[index + degree + 1] > knots[index + 1]
            ? scaled(dividePositive(
                       subtract(
                         source.poles[(i + (alongU ? 1 : 0)) * source.nv + j + (alongU ? 0 : 1)][k],
                         source.poles[i * source.nv + j][k]),
                       width),
                     degree)
            : point(
                0); // Empty-support basis at a repeated knot contributes nothing on an open span.
      out.poles.push_back(pole);
    }
  return out;
}

FocalSweep enclose(const std::array<Net, 6>& jets, const Cell& cell, double distance)
{
  std::array<std::array<Range, 4>, 6> h;
  for (int d = 0; d < 6; ++d)
  {
    const Net& net   = jets[d];
    const int  uspan = cell.uspan - (jets[0].p - net.p);
    const int  vspan = cell.vspan - (jets[0].q - net.q);
    const auto u     = basis(net.u, net.p, uspan, {cell.u0, cell.u1});
    const auto v     = basis(net.v, net.q, vspan, {cell.v0, cell.v1});
    h[d].fill(point(0));
    std::array<Range, 4> hull;
    for (auto& r : hull)
      r = {1e300, -1e300};
    for (int i = std::max(0, uspan - net.p); i <= uspan && i < net.nu; ++i)
      for (int j = std::max(0, vspan - net.q); j <= vspan && j < net.nv; ++j)
      {
        const Range b = multiply(u[i], v[j]);
        for (int k = 0; k < 4; ++k)
        {
          const Range pole = net.poles[i * net.nv + j][k];
          h[d][k]          = add(h[d][k], multiply(b, pole));
          hull[k].low      = std::min(hull[k].low, pole.low);
          hull[k].high     = std::max(hull[k].high, pole.high);
        }
      }
    // On this span the basis is nonnegative and sums to one. The hull and the
    // interval evaluation are independent enclosures; their intersection is tighter.
    for (int k = 0; k < 4; ++k)
      h[d][k] = intersect(h[d][k], hull[k]);
  }
  if (!valid(h[0][3]) || h[0][3].low <= 0)
    return FocalSweep::Unknown;
  std::array<Vec, 6> s;
  for (int k = 0; k < 3; ++k)
  {
    s[0][k] = dividePositive(h[0][k], h[0][3]);
    for (int d = 1; d <= 2; ++d)
      s[d][k] = dividePositive(subtract(h[d][k], multiply(s[0][k], h[d][3])), h[0][3]);
    for (int d = 3; d < 6; ++d)
    {
      const int a = d == 5 ? 2 : 1, b = d == 3 ? 1 : 2;
      s[d][k] = dividePositive(subtract(subtract(subtract(h[d][k], multiply(s[0][k], h[d][3])),
                                                 multiply(s[a][k], h[b][3])),
                                        multiply(s[b][k], h[a][3])),
                               h[0][3]);
    }
  }
  const Vec   normal = cross(s[1], s[2]);
  const Range length = root(norm2(normal));
  if (!valid(length) || length.low <= 0)
    return FocalSweep::Unknown;
  // I - d II is positive definite iff both 1 - d*k_i are positive.
  // Sylvester's criterion avoids subtracting nearly equal principal-curvature roots.
  const Range e =
    subtract(norm2(s[1]), scaled(dividePositive(dot(s[3], normal), length), distance));
  const Range f =
    subtract(dot(s[1], s[2]), scaled(dividePositive(dot(s[4], normal), length), distance));
  const Range g =
    subtract(norm2(s[2]), scaled(dividePositive(dot(s[5], normal), length), distance));
  const Range determinant = subtract(multiply(e, g), square(f));
  if (valid(e) && valid(g) && valid(determinant))
  {
    if (e.low > 0 && determinant.low > 0)
      return FocalSweep::Regular;
    if (e.high < 0 || g.high < 0 || determinant.high < 0)
      return FocalSweep::Crosses;
  }
  return FocalSweep::Unknown;
}
} // namespace

FocalSweep focalSweep(const TopoDS_Face& face, double distance)
{
  if (!std::isfinite(distance))
    return FocalSweep::Unknown;
  Net                         net;
  TopLoc_Location             location;
  const Handle(Geom_Surface)& surface = BRep_Tool::Surface(face, location);
  // Transform relative poles with interval arithmetic: Surface(face) would first
  // round the transformed poles in an ordinary floating-point copy.
  if (!readNet(surface, location.Transformation(), net))
    return FocalSweep::Unknown;
  std::array<Net, 6> jets;
  jets[0] = net;
  jets[1] = derivative(net, true);
  jets[2] = derivative(net, false);
  jets[3] = derivative(jets[1], true);
  jets[4] = derivative(jets[1], false);
  jets[5] = derivative(jets[2], false);
  if (face.Orientation() == TopAbs_REVERSED)
    distance = -distance;
  double u0, u1, v0, v1;
  BRepTools::UVBounds(face, u0, u1, v0, v1);
  if (!(u0 >= net.u[net.p] && u1 <= net.u[net.nu] && v0 >= net.v[net.q] && v1 <= net.v[net.nv]))
    return FocalSweep::Unknown;
  BRepTopAdaptor_FClass2d inside(face, Precision::Confusion());
  std::deque<Cell>        pending;
  for (int i = net.p; i < net.nu; ++i)
    for (int j = net.q; j < net.nv; ++j)
    {
      Cell c{std::max(u0, net.u[i]),
             std::min(u1, net.u[i + 1]),
             std::max(v0, net.v[j]),
             std::min(v1, net.v[j + 1]),
             i,
             j};
      if (c.u0 < c.u1 && c.v0 < c.v1)
        pending.push_back(c);
    }
  if (pending.empty())
    return FocalSweep::Unknown;
  // Work limits only turn a proof into Unknown. They never justify a refusal.
  bool unknown = false;
  int  visited = 0;
  while (!pending.empty() && visited++ < 4096)
  {
    checkCancelled("thicken");
    const Cell c = pending.front();
    pending.pop_front();
    const auto result = enclose(jets, c, distance);
    if (result == FocalSweep::Regular)
      continue;
    const double um = c.u0 + (c.u1 - c.u0) / 2, vm = c.v0 + (c.v1 - c.v0) / 2;
    if (result == FocalSweep::Crosses && inside.Perform(gp_Pnt2d(um, vm)) == TopAbs_IN)
      return FocalSweep::Crosses;
    // Also examine a point interval: a strict witness need not wait until the whole cell crosses.
    Cell pointCell = c;
    pointCell.u0 = pointCell.u1 = um;
    pointCell.v0 = pointCell.v1 = vm;
    if (inside.Perform(gp_Pnt2d(um, vm)) == TopAbs_IN
        && enclose(jets, pointCell, distance) == FocalSweep::Crosses)
      return FocalSweep::Crosses;
    if (c.depth >= 24 || !(um > c.u0 && um < c.u1 && vm > c.v0 && vm < c.v1))
    {
      unknown = true;
      continue;
    }
    Cell lo = c, hi = c;
    lo.depth++;
    hi.depth++;
    if ((c.u1 - c.u0) / (net.u[c.uspan + 1] - net.u[c.uspan])
        >= (c.v1 - c.v0) / (net.v[c.vspan + 1] - net.v[c.vspan]))
    {
      lo.u1 = um;
      hi.u0 = um;
    }
    else
    {
      lo.v1 = vm;
      hi.v0 = vm;
    }
    pending.push_back(hi);
    pending.push_back(lo);
  }
  return unknown || !pending.empty() ? FocalSweep::Unknown : FocalSweep::Regular;
}
} // namespace BRepExt
