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

#include <BRepSurfaceExt_FaceSeparation.hxx>
#include <BRepSurfaceExt_FaceSeparationAnalytic.hxx>

#include <BRepToolsExt_Session.hxx>
#include <GeomAlgoExt_Bernstein.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepOffsetExt_Parts.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BOPAlgo_ArgumentAnalyzer.hxx>
#include <BOPAlgo_CheckResult.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BndLib.hxx>
#include <BndLib_AddSurface.hxx>
#include <Bnd_Box.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <GeomConvert.hxx>
#include <GeomConvert_BSplineSurfaceToBezierSurface.hxx>
#include <Geom2d_Line.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <Geom_ElementarySurface.hxx>
#include <Geom_Plane.hxx>
#include <Geom_SurfaceOfLinearExtrusion.hxx>
#include <Geom_SurfaceOfRevolution.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <gp_Circ.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_BezierSurface.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Standard_Failure.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <NCollection_List.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <set>
#include <vector>

namespace BRepExt
{
namespace
{

// Work allowed to one pre-pass, in box tests and Bernstein coefficients: past it, what is not
// yet proven is left to the check.
constexpr long WORK = 2'000'000;
// How small halving makes the pieces of faces it tries to prove apart, or clear of a plane,
// before it gives up: a share of the smaller face (what a face is cut into Bezier patches does
// not matter, nor the scale it is drawn at), and never below a share of the gap the proof is to
// tell (pieces that close and that small are within two gaps of each other).
constexpr double FINEST_OF_FACE = 1.0 / 1024, FINEST_OF_GAP = 0.25;
// Refinements of a chart.
constexpr int CHART_ROUNDS = 40, CHART_CELLS = 8192;
// Bezier patches a rational face needs for a proof of it clear of itself to be worth trying.
constexpr size_t SELF_PATCHES = 16;
// How far from their shared edges two faces meeting there may come within their tolerances
// of each other, in their tolerances, still counted as meeting along them. Two planes at an
// angle t come within a gap of each other up to gap / sin t from their line, which the bound
// over-estimates about twofold at a right angle: eight admits creases of about 15 degrees and
// more. That is how close they come, not where they cross: their points in common are bounded
// apart from that, by how far the two sides' curves along the edge part (the bound's jump, the
// distance of the plane's edge poles), and that must be within the shared edge's own tolerance,
// where they are on it. Faces whose tolerances are loose (imported ones) cross there unseen no
// further from the edge than their own geometry is off it.
constexpr double NEAR_EDGE = 8;

// How far from a shared edge two faces may cross and still be on it: within its tolerance.
double edgeTube(const TopoDS_Shape& edge)
{
  return std::max(Precision::Confusion(), BRep_Tool::Tolerance(TopoDS::Edge(edge)));
}

// ---- Rational Bezier patches ----

struct Box
{
  gp_XYZ lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};

  void add(const gp_XYZ& p)
  {
    lo.SetCoord(std::min(lo.X(), p.X()), std::min(lo.Y(), p.Y()), std::min(lo.Z(), p.Z()));
    hi.SetCoord(std::max(hi.X(), p.X()), std::max(hi.Y(), p.Y()), std::max(hi.Z(), p.Z()));
  }

  void add(const Box& b)
  {
    add(b.lo);
    add(b.hi);
  }

  // Whether the boxes are further apart than `gap` along some axis.
  bool apart(const Box& b, double gap) const
  {
    for (int k = 1; k <= 3; ++k)
      if (b.lo.Coord(k) - hi.Coord(k) > gap || lo.Coord(k) - b.hi.Coord(k) > gap)
        return true;
    return false;
  }

  double size() const { return (hi - lo).Modulus(); }
};

// Poles (i, j) at i * (m + 1) + j, with positive weights. Positive weights keep the patch in the
// hull of its poles, and its corner poles are on it.
struct Patch
{
  int                 n = 0, m = 0;
  std::vector<gp_XYZ> poles;
  std::vector<double> weights;

  const gp_XYZ& pole(int i, int j) const { return poles[i * (m + 1) + j]; }

  double weight(int i, int j) const { return weights[i * (m + 1) + j]; }

  Box box() const
  {
    Box b;
    for (const auto& p : poles)
      b.add(p);
    return b;
  }
};

Patch patchOf(const Handle(Geom_BezierSurface)& bezier)
{
  Patch p;
  p.n = bezier->NbUPoles() - 1;
  p.m = bezier->NbVPoles() - 1;
  for (int i = 1; i <= p.n + 1; ++i)
    for (int j = 1; j <= p.m + 1; ++j)
    {
      p.poles.push_back(bezier->Pole(i, j).XYZ());
      p.weights.push_back(bezier->Weight(i, j));
    }
  return p;
}

// Re-indexed: u and v swapped, or one of them run backwards.
Patch transposed(const Patch& a)
{
  Patch p{a.m, a.n, {}, {}};
  for (int i = 0; i <= p.n; ++i)
    for (int j = 0; j <= p.m; ++j)
    {
      p.poles.push_back(a.pole(j, i));
      p.weights.push_back(a.weight(j, i));
    }
  return p;
}

Patch flipped(const Patch& a, bool alongU)
{
  Patch p{a.n, a.m, {}, {}};
  for (int i = 0; i <= p.n; ++i)
    for (int j = 0; j <= p.m; ++j)
    {
      const int ii = alongU ? a.n - i : i, jj = alongU ? j : a.m - j;
      p.poles.push_back(a.pole(ii, jj));
      p.weights.push_back(a.weight(ii, jj));
    }
  return p;
}

// De Casteljau on homogeneous poles, at `t` along u or v: the two parts, exactly sharing the cut.
std::pair<Patch, Patch> split(const Patch& a, bool alongU, double t)
{
  const int           count = alongU ? a.n : a.m, other = alongU ? a.m : a.n;
  Patch               lo = a, hi = a;
  std::vector<gp_XYZ> h(count + 1);
  std::vector<double> w(count + 1);
  for (int k = 0; k <= other; ++k)
  {
    const auto at = [&](int i) { return alongU ? i * (a.m + 1) + k : k * (a.m + 1) + i; };
    for (int i = 0; i <= count; ++i)
    {
      w[i] = a.weights[at(i)];
      h[i] = a.poles[at(i)] * w[i];
    }
    for (int level = 0; level <= count; ++level)
    {
      lo.weights[at(level)]         = w[0];
      lo.poles[at(level)]           = h[0] / w[0];
      hi.weights[at(count - level)] = w[count - level];
      hi.poles[at(count - level)]   = h[count - level] / w[count - level];
      for (int i = 0; i < count - level; ++i)
      {
        h[i] = h[i] * (1 - t) + h[i + 1] * t;
        w[i] = w[i] * (1 - t) + w[i + 1] * t;
      }
    }
  }
  return {lo, hi};
}

// The part of `a` over [from, to] of its parameter along u or v.
Patch piece(const Patch& a, bool alongU, double from, double to)
{
  if (from <= 0 && to >= 1)
    return a;
  auto head = to >= 1 ? a : split(a, alongU, to).first;
  return from <= 0 ? head : split(head, alongU, from / to).second;
}

// Whether the patch's poles spread more along u than along v.
bool spreadsAlongU(const Patch& a)
{
  double u = 0, v = 0;
  for (int i = 0; i <= a.n; ++i)
    for (int j = 0; j <= a.m; ++j)
    {
      if (i > 0)
        u = std::max(u, (a.pole(i, j) - a.pole(i - 1, j)).SquareModulus());
      if (j > 0)
        v = std::max(v, (a.pole(i, j) - a.pole(i, j - 1)).SquareModulus());
    }
  return u >= v;
}

// Halved across the way its poles spread most.
std::pair<Patch, Patch> halve(const Patch& a)
{
  return split(a, spreadsAlongU(a), .5);
}

// Whether the patches' pole hulls are further apart than `gap` along the line between their
// pole centres: what separates two thin strips at a crease, which axis boxes need not.
bool apartAlongCentres(const Patch& a, const Patch& b, double gap)
{
  gp_XYZ ca, cb;
  for (const auto& p : a.poles)
    ca += p;
  for (const auto& p : b.poles)
    cb += p;
  gp_XYZ       d      = cb / double(b.poles.size()) - ca / double(a.poles.size());
  const double length = d.Modulus();
  if (!(length > 0))
    return false;
  d /= length;
  double most = -1e300, least = 1e300;
  for (const auto& p : a.poles)
    most = std::max(most, d.Dot(p));
  for (const auto& p : b.poles)
    least = std::min(least, d.Dot(p));
  return least - most > gap;
}

// Whether two patches are further apart than `gap`, by halving the larger until their pole
// boxes are; false where the larger is no more than `finest` across first, or the work runs out.
bool apart(const Patch& a, const Patch& b, double gap, double finest, long& work)
{
  struct Pair
  {
    Patch a, b;
  };

  std::vector<Pair> stack{{a, b}};
  while (!stack.empty())
  {
    auto pair = std::move(stack.back());
    stack.pop_back();
    if (--work < 0)
      return false;
    const auto boxA = pair.a.box(), boxB = pair.b.box();
    if (boxA.apart(boxB, gap) || apartAlongCentres(pair.a, pair.b, gap))
      continue;
    const bool first = boxA.size() >= boxB.size();
    if ((first ? boxA : boxB).size() <= finest)
      return false;
    auto [one, two] = halve(first ? pair.a : pair.b);
    if (first)
    {
      stack.push_back({std::move(one), pair.b});
      stack.push_back({std::move(two), std::move(pair.b)});
    }
    else
    {
      stack.push_back({pair.a, std::move(one)});
      stack.push_back({std::move(pair.a), std::move(two)});
    }
  }
  return true;
}

bool rational(const Patch& p)
{
  const auto [lo, hi] = std::minmax_element(p.weights.begin(), p.weights.end());
  return *hi - *lo > 1e-14 * *hi;
}

// The numerators of the patch's derivatives along u and v, both of degrees (N, M): the
// derivatives are these over w^2, with w between `wMin` and `wMax`.
struct Derivatives
{
  int                 N = 0, M = 0;
  std::vector<gp_XYZ> u, v;
  double              wMin = 1, wMax = 1;
};

Derivatives derivatives(const Patch& p, long& work)
{
  Derivatives         d;
  const int           n = p.n, m = p.m;
  const bool          weighted = rational(p);
  std::vector<gp_XYZ> h(p.poles.size());
  std::vector<double> w = weighted ? p.weights : std::vector<double>(p.poles.size(), 1.0);
  for (size_t k = 0; k < h.size(); ++k)
    h[k] = p.poles[k] * w[k];
  std::vector<gp_XYZ> hu(n * (m + 1)), hv((n + 1) * m);
  std::vector<double> wu(n * (m + 1)), wv((n + 1) * m);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j <= m; ++j)
    {
      hu[i * (m + 1) + j] = (h[(i + 1) * (m + 1) + j] - h[i * (m + 1) + j]) * n;
      wu[i * (m + 1) + j] = (w[(i + 1) * (m + 1) + j] - w[i * (m + 1) + j]) * n;
    }
  for (int i = 0; i <= n; ++i)
    for (int j = 0; j < m; ++j)
    {
      hv[i * m + j] = (h[i * (m + 1) + j + 1] - h[i * (m + 1) + j]) * m;
      wv[i * m + j] = (w[i * (m + 1) + j + 1] - w[i * (m + 1) + j]) * m;
    }
  if (!weighted)
  {
    d.N = n;
    d.M = m;
    d.u = GeomAlgoExt::elevate(hu, n - 1, m, n, m);
    d.v = GeomAlgoExt::elevate(hv, n, m - 1, n, m);
    work -= long(d.u.size()) * 2;
    return d;
  }
  // (h / w)' = (h' w - h w') / w^2.
  auto u  = GeomAlgoExt::product(hu, n - 1, m, w, n, m),
       u2 = GeomAlgoExt::product(h, n, m, wu, n - 1, m);
  auto v  = GeomAlgoExt::product(hv, n, m - 1, w, n, m),
       v2 = GeomAlgoExt::product(h, n, m, wv, n, m - 1);
  for (size_t k = 0; k < u.size(); ++k)
    u[k] -= u2[k];
  for (size_t k = 0; k < v.size(); ++k)
    v[k] -= v2[k];
  d.N                 = 2 * n;
  d.M                 = 2 * m;
  d.u                 = GeomAlgoExt::elevate(u, 2 * n - 1, 2 * m, d.N, d.M);
  d.v                 = GeomAlgoExt::elevate(v, 2 * n, 2 * m - 1, d.N, d.M);
  const auto [lo, hi] = std::minmax_element(p.weights.begin(), p.weights.end());
  d.wMin              = *lo;
  d.wMax              = *hi;
  work -= long(n * m) * long(n * m) * 8 + long(d.u.size()) * 2;
  return d;
}

// A tensor Bernstein polynomial of degrees (n, m) cut at t along u or v by de Casteljau.
template <class T>
std::pair<std::vector<T>, std::vector<T>> splitCoefficients(const std::vector<T>& a,
                                                            int                   n,
                                                            int                   m,
                                                            bool                  alongU,
                                                            double                t)
{
  const int      count = alongU ? n : m, other = alongU ? m : n;
  std::vector<T> lo(a.size()), hi(a.size()), h(count + 1);
  for (int k = 0; k <= other; ++k)
  {
    const auto at = [&](int i) { return alongU ? i * (m + 1) + k : k * (m + 1) + i; };
    for (int i = 0; i <= count; ++i)
      h[i] = a[at(i)];
    for (int level = 0; level <= count; ++level)
    {
      lo[at(level)]         = h[0];
      hi[at(count - level)] = h[count - level];
      for (int i = 0; i < count - level; ++i)
        h[i] = h[i] * (1 - t) + h[i + 1] * t;
    }
  }
  return {lo, hi};
}

// The derivatives of the two halves of a patch at the middle, from its own: the numerators
// restricted, those across the cut at half the speed, at a fraction of the cost of new ones.
std::pair<Derivatives, Derivatives> halves(const Derivatives& d,
                                           const Patch&       lo,
                                           const Patch&       hi,
                                           bool               alongU)
{
  auto [ulo, uhi] = splitCoefficients(d.u, d.N, d.M, alongU, .5);
  auto [vlo, vhi] = splitCoefficients(d.v, d.N, d.M, alongU, .5);
  for (auto* part : {&ulo, &uhi, &vlo, &vhi})
    if ((part == &ulo || part == &uhi) == alongU)
      for (auto& x : *part)
        x *= .5;
  const auto weights = [&](const Patch& p, Derivatives out) {
    if (d.wMin != 1 || d.wMax != 1 || rational(p))
    {
      const auto [a, b] = std::minmax_element(p.weights.begin(), p.weights.end());
      out.wMin          = *a;
      out.wMax          = *b;
    }
    return out;
  };
  return {weights(lo, {d.N, d.M, std::move(ulo), std::move(vlo), 1, 1}),
          weights(hi, {d.N, d.M, std::move(uhi), std::move(vhi), 1, 1})};
}

// How far apart the curves along two patches' sides can be: `a` along its side `sideA`, `b`
// along `sideB` (0 u = 0, 1 u = 1, 2 v = 0, 3 v = 1), both run the same way. The difference
// of two rational curves is (ha wb - hb wa) / (wa wb), bounded coefficient by coefficient.
double mismatch(const Patch& a, int sideA, const Patch& b, int sideB)
{
  const auto curve = [](const Patch& p, int side, std::vector<gp_XYZ>& h, std::vector<double>& w) {
    const bool alongV = side < 2;
    const int  count  = alongV ? p.m : p.n;
    for (int k = 0; k <= count; ++k)
    {
      const int i = alongV ? (side == 0 ? 0 : p.n) : k, j = alongV ? k : (side == 2 ? 0 : p.m);
      w.push_back(p.weight(i, j));
      h.push_back(p.pole(i, j) * p.weight(i, j));
    }
    return count;
  };
  std::vector<gp_XYZ> ha, hb;
  std::vector<double> wa, wb;
  const int           da = curve(a, sideA, ha, wa), db = curve(b, sideB, hb, wb);
  auto                one = GeomAlgoExt::product(ha, 0, da, wb, 0, db),
       two                = GeomAlgoExt::product(hb, 0, db, wa, 0, da);
  double largest          = 0;
  for (size_t k = 0; k < one.size(); ++k)
    largest = std::max(largest, (one[k] - two[k]).Modulus());
  return largest
         / (*std::min_element(wa.begin(), wa.end()) * *std::min_element(wb.begin(), wb.end()));
}

// ---- Faces ----

// An edge occurrence on a side of its face's parameter rectangle (0 u0, 1 u1, 2 v0, 3 v1):
// the coordinate along the side is scale * t + offset at the edge's parameter t.
struct Side
{
  TopoDS_Edge edge;
  int         side  = -1;
  double      scale = 0, offset = 0, first = 0, last = 0;
};

struct FaceData
{
  TopoDS_Face face;

  enum Kind
  {
    Spline,
    Plane,
    Hull,
    Unknown
  } kind = Unknown;

  double tolerance = 0;
  Box    box;
  gp_Pln plane;
  double u0 = 0, u1 = 0, v0 = 0, v1 = 0;
  // Bezier patches by u span then v span, and the spans' breaks in u and v.
  int                 nu = 0, nv = 0;
  std::vector<Patch>  patches;
  std::vector<double> uBreaks, vBreaks;
  // Spline faces whose every edge runs along a whole side of the parameter rectangle.
  bool              rectangle = false;
  std::vector<Side> sides;
  ShapeSet          edges, vertices;
};

std::optional<Side> sideOf(const TopoDS_Edge& edge, const FaceData& f)
{
  double first, last;
  auto   curve = BRep_Tool::CurveOnSurface(edge, f.face, first, last);
  while (!curve.IsNull() && curve->IsKind(STANDARD_TYPE(Geom2d_TrimmedCurve)))
    curve = Handle(Geom2d_TrimmedCurve)::DownCast(curve)->BasisCurve();
  const auto line = Handle(Geom2d_Line)::DownCast(curve);
  if (line.IsNull())
    return std::nullopt;
  const gp_Pnt2d at    = line->Location();
  const gp_Dir2d along = line->Direction();
  const double   du = 1e-9 * (1 + std::abs(f.u1 - f.u0)), dv = 1e-9 * (1 + std::abs(f.v1 - f.v0));
  Side           s{edge, -1, 0, 0, first, last};
  if (std::abs(along.X()) < 1e-12)
  {
    if (std::abs(at.X() - f.u0) <= du)
      s.side = 0;
    else if (std::abs(at.X() - f.u1) <= du)
      s.side = 1;
    s.scale  = along.Y();
    s.offset = at.Y();
  }
  else if (std::abs(along.Y()) < 1e-12)
  {
    if (std::abs(at.Y() - f.v0) <= dv)
      s.side = 2;
    else if (std::abs(at.Y() - f.v1) <= dv)
      s.side = 3;
    s.scale  = along.X();
    s.offset = at.X();
  }
  if (s.side < 0)
    return std::nullopt;
  // The whole side: the rectangle's faces are what charts put side by side.
  const double a = s.scale * first + s.offset, b = s.scale * last + s.offset;
  const double lo = s.side < 2 ? f.v0 : f.u0, hi = s.side < 2 ? f.v1 : f.u1,
               d = s.side < 2 ? dv : du;
  if (std::abs(std::min(a, b) - lo) > d || std::abs(std::max(a, b) - hi) > d)
    return std::nullopt;
  return s;
}

void bezierGrid(FaceData& f, const Handle(Geom_BSplineSurface)& surface)
{
  auto opened = Handle(Geom_BSplineSurface)::DownCast(surface->Copy());
  if (opened->IsUPeriodic())
    opened->SetUNotPeriodic();
  if (opened->IsVPeriodic())
    opened->SetVNotPeriodic();
  GeomConvert_BSplineSurfaceToBezierSurface convert(opened, f.u0, f.u1, f.v0, f.v1, 1e-10);
  f.nu = convert.NbUPatches();
  f.nv = convert.NbVPatches();
  NCollection_Array1<double> uk(1, f.nu + 1), vk(1, f.nv + 1);
  convert.UKnots(uk);
  convert.VKnots(vk);
  for (int i = 1; i <= f.nu + 1; ++i)
    f.uBreaks.push_back(uk(i));
  for (int j = 1; j <= f.nv + 1; ++j)
    f.vBreaks.push_back(vk(j));
  for (int i = 1; i <= f.nu; ++i)
    for (int j = 1; j <= f.nv; ++j)
      f.patches.push_back(patchOf(convert.Patch(i, j)));
}

// A box that encloses `curve` over [t0, t1]: in closed form for a line, a circle or an ellipse,
// the hull of its poles for a Bezier or B-spline curve; void for any other curve.
Bnd_Box enclosing(Handle(Geom_Curve) curve, double t0, double t1)
{
  Bnd_Box box;
  while (curve->IsKind(STANDARD_TYPE(Geom_TrimmedCurve)))
    curve = Handle(Geom_TrimmedCurve)::DownCast(curve)->BasisCurve();
  const GeomAdaptor_Curve adaptor(curve);
  switch (adaptor.GetType())
  {
    case GeomAbs_Line:
      box.Add(curve->Value(t0));
      box.Add(curve->Value(t1));
      break;
    case GeomAbs_Circle:
      BndLib::Add(adaptor.Circle(), t0, t1, 0, box);
      break;
    case GeomAbs_Ellipse:
      BndLib::Add(adaptor.Ellipse(), t0, t1, 0, box);
      break;
    case GeomAbs_BezierCurve:
    case GeomAbs_BSplineCurve: {
      // All its poles where it is periodic: [t0, t1] may run past its first period.
      auto spline = GeomConvert::CurveToBSplineCurve(curve);
      if (!spline->IsPeriodic() && (t0 > spline->FirstParameter() || t1 < spline->LastParameter()))
        spline->Segment(std::max(t0, spline->FirstParameter()),
                        std::min(t1, spline->LastParameter()));
      for (int i = 1; i <= spline->NbPoles(); ++i)
        box.Add(spline->Pole(i));
      break;
    }
    default:
      break;
  }
  return box;
}

// A box that encloses the face's surface over its parameter bounds, where one is known: OCCT's
// for an elementary surface, which it works out in closed form; for a surface of extrusion or of
// revolution, its basis curve's (enclosing) swept along the direction or round the axis, which
// encloses it as the sweep is linear in the curve's points. Void for other surfaces, whose OCCT
// boxes are samples of them, or of their basis curves, enlarged by an estimate
// (GeomBndLib_OtherSurface, GeomBndLib_BSplineCurve): no bound.
Bnd_Box enclosing(const TopoDS_Face& face, double u0, double u1, double v0, double v1)
{
  TopLoc_Location location;
  auto            surface = BRep_Tool::Surface(face, location);
  while (surface->IsKind(STANDARD_TYPE(Geom_RectangularTrimmedSurface)))
    surface = Handle(Geom_RectangularTrimmedSurface)::DownCast(surface)->BasisSurface();
  Bnd_Box local, box;
  if (surface->IsKind(STANDARD_TYPE(Geom_ElementarySurface))
      && !surface->IsKind(STANDARD_TYPE(Geom_Plane)))
    BndLib_AddSurface::Add(GeomAdaptor_Surface(surface), u0, u1, v0, v1, 0, local);
  else if (const auto extrusion = Handle(Geom_SurfaceOfLinearExtrusion)::DownCast(surface);
           !extrusion.IsNull())
  {
    const Bnd_Box curve = enclosing(extrusion->BasisCurve(), u0, u1);
    if (curve.IsVoid())
      return box;
    const gp_XYZ along                  = extrusion->Direction().XYZ();
    const auto [x0, x1, y0, y1, z0, z1] = curve.Get();
    for (const double x : {x0, x1})
      for (const double y : {y0, y1})
        for (const double z : {z0, z1})
          for (const double v : {v0, v1})
            local.Add(gp_Pnt(gp_XYZ(x, y, z) + along * v));
  }
  else if (const auto revolution = Handle(Geom_SurfaceOfRevolution)::DownCast(surface);
           !revolution.IsNull())
  {
    const Bnd_Box curve = enclosing(revolution->BasisCurve(), v0, v1);
    if (curve.IsVoid())
      return box;
    const gp_Ax1 axis                   = revolution->Axis();
    const auto [x0, x1, y0, y1, z0, z1] = curve.Get();
    for (const double x : {x0, x1})
      for (const double y : {y0, y1})
        for (const double z : {z0, z1})
        {
          // The corner's arc about the axis over [u0, u1], a point where it is on it.
          const gp_Pnt corner(x, y, z);
          const gp_XYZ from = axis.Location().XYZ(), up = axis.Direction().XYZ();
          const gp_Pnt centre(from + up * up.Dot(corner.XYZ() - from));
          const double radius = centre.Distance(corner);
          if (radius <= Precision::Confusion())
            local.Add(corner);
          else
            BndLib::Add(gp_Circ(gp_Ax2(centre, axis.Direction(), gp_Vec(centre, corner)), radius),
                        u0,
                        u1,
                        0,
                        local);
        }
  }
  if (local.IsVoid() || local.IsOpen())
    return box;
  const auto [x0, x1, y0, y1, z0, z1] = local.Get();
  for (const double x : {x0, x1})
    for (const double y : {y0, y1})
      for (const double z : {z0, z1})
        box.Add(gp_Pnt(x, y, z).Transformed(location.Transformation()));
  return box;
}

FaceData faceData(const TopoDS_Face& face)
{
  FaceData f;
  f.face      = face;
  f.tolerance = tolerance::largest(face);
  TopExp::MapShapes(face, TopAbs_EDGE, f.edges);
  TopExp::MapShapes(face, TopAbs_VERTEX, f.vertices);
  if (!std::isfinite(f.tolerance))
    return f;
  try
  {
    BRepTools::UVBounds(face, f.u0, f.u1, f.v0, f.v1);
    BRepAdaptor_Surface surface(face);
    const auto          type = surface.GetType();
    // The face's own surface, its location applied.
    auto own = BRep_Tool::Surface(face);
    while (own->IsKind(STANDARD_TYPE(Geom_RectangularTrimmedSurface)))
      own = Handle(Geom_RectangularTrimmedSurface)::DownCast(own)->BasisSurface();
    const auto spline = Handle(Geom_BSplineSurface)::DownCast(own);
    const auto bezier = Handle(Geom_BezierSurface)::DownCast(own);
    if (type == GeomAbs_Plane)
    {
      // The face lies in its plane's rectangle of its parameter bounds.
      f.kind  = FaceData::Plane;
      f.plane = surface.Plane();
      for (double u : {f.u0, f.u1})
        for (double v : {f.v0, f.v1})
          f.box.add(surface.Value(u, v).XYZ());
    }
    else if (!spline.IsNull() || !bezier.IsNull())
    {
      f.kind = FaceData::Spline;
      bezierGrid(f, !spline.IsNull() ? spline : GeomConvert::SurfaceToBSplineSurface(bezier));
      unsigned covered = 0;
      bool     lines   = true;
      for (TopExp_Explorer it(face, TopAbs_EDGE); it.More() && lines; it.Next())
      {
        const auto side = sideOf(TopoDS::Edge(it.Current()), f);
        if (!side)
          lines = false;
        else
        {
          covered |= 1u << side->side;
          f.sides.push_back(*side);
        }
      }
      f.rectangle = lines && covered == 15;
    }
    else
    {
      // Any other surface only by a box that encloses it: the check intersects analytic
      // surfaces in closed form, where proofs gain it little. None for surfaces that have
      // no such box, which stay Unknown, for the check.
      const Bnd_Box box = enclosing(face, f.u0, f.u1, f.v0, f.v1);
      if (!box.IsVoid())
      {
        f.kind = FaceData::Hull;
        double x0, y0, z0, x1, y1, z1;
        box.Get(x0, y0, z0, x1, y1, z1);
        f.box.add(gp_XYZ(x0, y0, z0));
        f.box.add(gp_XYZ(x1, y1, z1));
      }
    }
  }
  catch (const Standard_Failure&)
  {
    f.kind = FaceData::Unknown;
    return f;
  }
  for (const auto& p : f.patches)
    f.box.add(p.box());
  return f;
}

// Two faces' tolerances, and an allowance for the rounding of their poles' arithmetic.
double gapOf(const FaceData& a, const FaceData& b)
{
  return a.tolerance + b.tolerance + 1e-12 * (a.box.size() + b.box.size());
}

// How small pieces of the two faces get before a proof gives up (FINEST_OF_FACE).
double finestOf(const FaceData& a, const FaceData& b, double gap)
{
  return std::max(FINEST_OF_FACE * std::min(a.box.size(), b.box.size()), FINEST_OF_GAP * gap);
}

// ---- Charts: faces put side by side in one rectangle of cells ----

struct Cell
{
  Patch                      patch;
  int                        face = 0;
  std::optional<Derivatives> derivatives;
  // Unique to the cell as it stands: a halved cell's halves are new ones.
  long id = next();

  static long next()
  {
    static long count = 0;
    return ++count;
  }
};

struct Chart
{
  int  columns = 0, rows = 0;
  bool wrapS = false, wrapT = false;
  // The columns at which the second face starts: glue lines with no derivative continuity.
  int               split = -1;
  std::vector<Cell> cells; // column by column

  Cell& at(int c, int r) { return cells[c * rows + r]; }

  const Cell& at(int c, int r) const { return cells[c * rows + r]; }

  int column(int c) const { return wrapS ? (c % columns + columns) % columns : c; }

  int row(int r) const { return wrapT ? (r % rows + rows) % rows : r; }
};

// A cell halved at its middle along s or t, its derivatives carried over where it has them.
std::pair<Cell, Cell> halvesOf(const Cell& cell, bool alongS)
{
  auto [a, b] = split(cell.patch, alongS, .5);
  Cell one{std::move(a), cell.face, {}, Cell::next()},
    two{std::move(b), cell.face, {}, Cell::next()};
  if (cell.derivatives)
  {
    auto [da, db]   = halves(*cell.derivatives, one.patch, two.patch, alongS);
    one.derivatives = std::move(da);
    two.derivatives = std::move(db);
  }
  return {std::move(one), std::move(two)};
}

// How wide a column is: the longest chord across one of its cells.
double width(const Chart& chart, int c)
{
  double out = 0;
  for (int r = 0; r < chart.rows; ++r)
  {
    const auto& p = chart.at(c, r).patch;
    for (int j = 0; j <= p.m; ++j)
      out = std::max(out, (p.pole(p.n, j) - p.pole(0, j)).Modulus());
  }
  return out;
}

// How tall a row is beside column c: the longest chord along one of its cells there.
double height(const Chart& chart, int c, int r)
{
  const auto& p   = chart.at(c, r).patch;
  double      out = 0;
  for (int i = 0; i <= p.n; ++i)
    out = std::max(out, (p.pole(i, p.m) - p.pole(i, 0)).Modulus());
  return out;
}

// Halves the rows, or the columns, that `marked` picks: cells stay a chart however wide each
// row or column is, as its own parameter runs across it.
void refine(Chart& chart, bool alongS, const std::vector<char>& marked)
{
  Chart out = chart;
  out.cells.clear();
  if (alongS)
  {
    out.columns = 0;
    for (int c = 0; c < chart.columns; ++c)
    {
      if (c == chart.split)
        out.split = out.columns;
      if (!marked[c])
      {
        for (int r = 0; r < chart.rows; ++r)
          out.cells.push_back(chart.at(c, r));
        ++out.columns;
        continue;
      }
      std::vector<Cell> left, right;
      for (int r = 0; r < chart.rows; ++r)
      {
        auto [a, b] = halvesOf(chart.at(c, r), true);
        left.push_back(std::move(a));
        right.push_back(std::move(b));
      }
      out.cells.insert(out.cells.end(), left.begin(), left.end());
      out.cells.insert(out.cells.end(), right.begin(), right.end());
      out.columns += 2;
    }
  }
  else
  {
    out.rows = 0;
    for (int r = 0; r < chart.rows; ++r)
      out.rows += marked[r] ? 2 : 1;
    for (int c = 0; c < chart.columns; ++c)
      for (int r = 0; r < chart.rows; ++r)
      {
        if (!marked[r])
        {
          out.cells.push_back(chart.at(c, r));
          continue;
        }
        auto [a, b] = halvesOf(chart.at(c, r), false);
        out.cells.push_back(std::move(a));
        out.cells.push_back(std::move(b));
      }
  }
  chart = std::move(out);
}

void refine(Chart& chart, bool alongS)
{
  refine(chart, alongS, std::vector<char>(alongS ? chart.columns : chart.rows, 1));
}

// Whether the cells of the block at column c, row r (two by two, or fewer at an open edge)
// map one to one, with the bound that proves it. A segment between two of its points has
// G(b) - G(a) = integral of M J (b - a), plus at most a jump at the one vertical and the one
// horizontal cell boundary it may cross; with the symmetric part of M J at least lambda over
// every cell, (b - a) . (G(b) - G(a)) >= lambda |b - a|^2 - |M| jump |b - a|.
struct BlockBound
{
  bool   proven  = false;
  bool   refineS = true;
  double lambda = 0, norm = 0, jump = 0, speedS = 0;
};

BlockBound block(Chart& chart, int c, int r, long& work)
{
  BlockBound                      out;
  std::vector<std::array<int, 3>> cells; // chart column, row, and the offset along s
  const int                       cs = chart.wrapS || c + 1 < chart.columns ? 2 : 1,
            rs                       = chart.wrapT || r + 1 < chart.rows ? 2 : 1;
  for (int dc = 0; dc < cs; ++dc)
    for (int dr = 0; dr < rs; ++dr)
      cells.push_back({chart.column(c + dc), chart.row(r + dr), dc});
  gp_XYZ meanS, meanT;
  double terms = 0;
  for (const auto& [cc, rr, dc] : cells)
  {
    auto& cell = chart.at(cc, rr);
    if (!cell.derivatives)
      cell.derivatives = derivatives(cell.patch, work);
    const double w2 = cell.derivatives->wMax * cell.derivatives->wMax;
    for (const auto& d : cell.derivatives->u)
      meanS += d / w2;
    for (const auto& d : cell.derivatives->v)
      meanT += d / w2;
    terms += double(cell.derivatives->u.size());
  }
  meanS /= terms;
  meanT /= terms;

  // M = (T^t T)^-1 T^t: the mean derivatives to the unit axes, whatever the crease between.
  const double a = meanS.Dot(meanS), b = meanS.Dot(meanT), d = meanT.Dot(meanT),
               det = a * d - b * b;
  if (!(det > 1e-24 * a * d))
    return out;
  const gp_XYZ m0 = (meanS * d - meanT * b) / det, m1 = (meanT * a - meanS * b) / det;
  // |M|, its largest singular value.
  {
    const double p = m0.SquareModulus(), q = m0.Dot(m1), r = m1.SquareModulus();
    out.norm = std::sqrt((p + r) / 2 + std::sqrt((p - r) * (p - r) / 4 + q * q));
  }
  // What spoils M J being the identity: the s derivatives' share along t, which halving along s
  // shrinks, and the t derivatives' along s, their turning along t, which halving along t does.
  double leakS = 0, leakT = 0, lambda = 1e300;
  for (const auto& [cc, rr, dc] : cells)
  {
    const auto& cell  = chart.at(cc, rr);
    const auto& dv    = *cell.derivatives;
    double      least = 1e300, scale = 0;
    for (size_t k = 0; k < dv.u.size(); ++k)
    {
      const double c00 = m0.Dot(dv.u[k]), c10 = m1.Dot(dv.u[k]), c01 = m0.Dot(dv.v[k]),
                   c11 = m1.Dot(dv.v[k]);
      const double h01 = (c01 + c10) / 2, mid = (c00 + c11) / 2, half = (c00 - c11) / 2;
      leakS       = std::max(leakS, std::abs(c10));
      leakT       = std::max(leakT, std::abs(c01));
      out.refineS = leakS >= leakT;
      least       = std::min(least, mid - std::sqrt(half * half + h01 * h01));
      scale       = std::max(scale, std::abs(c00) + std::abs(c11) + std::abs(h01));
      out.speedS  = std::max(out.speedS, dv.u[k].Modulus() / (dv.wMin * dv.wMin));
    }
    work -= long(dv.u.size());
    if (!(least > 1e-10 * scale))
      return out;
    lambda = std::min(lambda, least / (dv.wMax * dv.wMax));
  }
  // The jumps across the block's inner boundaries: zero in exact arithmetic within a face.
  double jumpS = 0, jumpT = 0;
  if (cs == 2)
    for (int dr = 0; dr < rs; ++dr)
      jumpS = std::max(jumpS,
                       mismatch(chart.at(chart.column(c), chart.row(r + dr)).patch,
                                1,
                                chart.at(chart.column(c + 1), chart.row(r + dr)).patch,
                                0));
  if (rs == 2)
    for (int dc = 0; dc < cs; ++dc)
      jumpT = std::max(jumpT,
                       mismatch(chart.at(chart.column(c + dc), chart.row(r)).patch,
                                3,
                                chart.at(chart.column(c + dc), chart.row(r + 1)).patch,
                                2));
  out.jump   = jumpS + jumpT;
  out.lambda = lambda;
  out.proven = true;
  return out;
}

// Proves the chart's blocks (those `wanted` picks) one to one, refining it as they ask, and
// then its cells that are no block's neighbours apart. `nearGlue` bounds, for blocks across
// the glue line, how far from it the two sides may come within `gap` of each other, `onGlue`
// how far from it they may cross.
template <class Wanted, class Far>
bool proveChart(Chart& chart,
                Wanted wanted,
                Far    far,
                double gap,
                double finest,
                double nearGlue,
                double onGlue,
                long&  work)
{
  // A chart of two faces asks only for blocks across its glue line: it halves columns there.
  const bool glueOnly = chart.split >= 0;
  // Derivatives once, at the coarsest: halving carries them down at a fraction of the cost.
  for (auto& cell : chart.cells)
    if (!cell.derivatives)
      cell.derivatives = derivatives(cell.patch, work);
  // A wrapped direction needs three cells for no block to hold one twice.
  while (chart.wrapS && chart.columns < 3)
    refine(chart, true);
  while (chart.wrapT && chart.rows < 3)
    refine(chart, false);
  // Blocks proven, by their cells: a block whose cells no round halved stays proven.
  std::set<std::array<long, 4>> held;
  const auto                    key = [&](int c, int r) {
    std::array<long, 4> out{};
    const int           cs = chart.wrapS || c + 1 < chart.columns ? 2 : 1,
              rs           = chart.wrapT || r + 1 < chart.rows ? 2 : 1;
    for (int dc = 0; dc < cs; ++dc)
      for (int dr = 0; dr < rs; ++dr)
        out[dc * 2 + dr] = chart.at(chart.column(c + dc), chart.row(r + dr)).id;
    return out;
  };
  // How far the worst glue block's reach is over what is allowed, by round (infinite while a
  // glue block has no bound): where halving no longer brings it down, the faces meet too flat
  // for the bound, and the check decides.
  std::vector<double> over;
  for (int round = 0; round <= CHART_ROUNDS; ++round)
  {
    // The rows and columns of every block that does not hold, each to halve the way it asks.
    std::vector<char> rows(chart.rows, 0), columns(chart.columns, 0);
    bool              again = false;
    double            worst = 0;
    for (int c = 0; c < chart.columns; ++c)
      for (int r = 0; r < chart.rows; ++r)
      {
        if (!wanted(chart, c, r))
          continue;
        const auto id = key(c, r);
        if (held.count(id))
          continue;
        const auto bound = block(chart, c, r, work);
        if (work < 0)
          return false;
        bool       holds = bound.proven;
        const bool glue =
          chart.split >= 0
          && (chart.column(c + 1) == chart.split || (chart.wrapS && chart.column(c + 1) == 0));
        if (holds && glue)
        {
          // Points of the two sides within `gap` of each other are within r of each
          // other in the chart, so within r times the speed along s of the glue line.
          // Points in common, within the share of that the jump alone makes.
          const double reach = bound.norm * (gap + bound.jump) / bound.lambda * bound.speedS;
          const double cross = bound.norm * bound.jump / bound.lambda * bound.speedS;
          holds              = reach <= nearGlue && cross <= onGlue;
          worst              = std::max({worst, reach / nearGlue, cross / onGlue});
        }
        else if (glue)
        {
          // No bound at all yet: nothing to tell a stall by.
          worst = std::numeric_limits<double>::infinity();
        }
        if (holds)
        {
          held.insert(id);
          continue;
        }
        again = true;
        // A row four times as tall as the one beside it or more (a sliver, where the two
        // faces' breaks along their edge are a hair apart) keeps the derivatives of the
        // block jumping between them however both are halved: the taller is halved
        // alone, as the wider column beside a glue line is.
        if (chart.wrapT || r + 1 < chart.rows)
        {
          const int next = chart.row(r + 1), cs = chart.wrapS || c + 1 < chart.columns ? 2 : 1;
          double    here = 0, there = 0;
          for (int dc = 0; dc < cs; ++dc)
          {
            here  = std::max(here, height(chart, chart.column(c + dc), r));
            there = std::max(there, height(chart, chart.column(c + dc), next));
          }
          if (std::min(here, there) < .25 * std::max(here, there))
          {
            rows[here > there ? r : next] = 1;
            continue;
          }
        }
        if (bound.refineS)
        {
          columns[c] = 1;
          if (chart.wrapS || c + 1 < chart.columns)
            columns[chart.column(c + 1)] = 1;
        }
        else
        {
          rows[r] = 1;
          if (chart.wrapT || r + 1 < chart.rows)
            rows[chart.row(r + 1)] = 1;
        }
      }
    if (!again)
    {
      // Pairs whose boxes come within the gap, swept along x.
      const int        count = int(chart.cells.size());
      std::vector<Box> boxes;
      std::vector<int> order(count);
      for (int p = 0; p < count; ++p)
      {
        boxes.push_back(chart.cells[p].patch.box());
        order[p] = p;
      }
      std::sort(order.begin(), order.end(), [&](int p, int q) {
        return boxes[p].lo.X() < boxes[q].lo.X();
      });
      for (int a = 0; a < count; ++a)
        for (int b = a + 1; b < count && boxes[order[b]].lo.X() <= boxes[order[a]].hi.X() + gap;
             ++b)
        {
          const int p = std::min(order[a], order[b]), q = std::max(order[a], order[b]);
          if (--work < 0)
            return false;
          if (!far(chart, p, q) || boxes[p].apart(boxes[q], gap))
            continue;
          if (!apart(chart.cells[p].patch, chart.cells[q].patch, gap, finest, work))
            return false;
        }
      return true;
    }
    if (round == CHART_ROUNDS || chart.cells.size() * 2 > CHART_CELLS)
      return false;
    over.push_back(worst);
    if (over.size() > 2 && worst > 1 && std::isfinite(worst) && std::isfinite(over[over.size() - 3])
        && worst > .8 * over[over.size() - 3])
      return false;
    if (glueOnly && std::count(columns.begin(), columns.end(), 1))
    {
      // Of the two columns beside the glue line, the wider where one is much the wider:
      // the two sides bound it best when their columns are about as wide.
      const double left = width(chart, chart.split - 1), right = width(chart, chart.split);
      columns.assign(chart.columns, 0);
      columns[chart.split - 1] = left >= .75 * std::max(left, right);
      columns[chart.split]     = right >= .75 * std::max(left, right);
    }
    if (std::count(columns.begin(), columns.end(), 1))
      refine(chart, true, columns);
    if (std::count(rows.begin(), rows.end(), 1))
      refine(chart, false, rows);
  }
  return false;
}

bool neighbours(const Chart& chart, int p, int q)
{
  const int  cp = p / chart.rows, rp = p % chart.rows, cq = q / chart.rows, rq = q % chart.rows;
  const auto near = [](int a, int b, int count, bool wrap) {
    const int d = std::abs(a - b);
    return d <= 1 || (wrap && count - d <= 1);
  };
  return near(cp, cq, chart.columns, chart.wrapS) && near(rp, rq, chart.rows, chart.wrapT);
}

// The face's patches as cells, transposed and flipped so that what is asked for is its s and
// t, t increasing with `t` of the native coordinate along t.
Chart cellsOf(const FaceData& f, int index, bool transpose, bool flipS, bool flipT)
{
  Chart chart;
  chart.columns = transpose ? f.nv : f.nu;
  chart.rows    = transpose ? f.nu : f.nv;
  for (int c = 0; c < chart.columns; ++c)
    for (int r = 0; r < chart.rows; ++r)
    {
      const int cc = flipS ? chart.columns - 1 - c : c, rr = flipT ? chart.rows - 1 - r : r;
      const int i = transpose ? rr : cc, j = transpose ? cc : rr;
      auto      patch = f.patches[i * f.nv + j];
      if (transpose)
        patch = transposed(patch);
      if (flipS)
        patch = flipped(patch, true);
      if (flipT)
        patch = flipped(patch, false);
      chart.cells.push_back({std::move(patch), index, {}});
    }
  return chart;
}

// Whether the face's same edge runs along both of two opposite sides with the same
// coordinate along them: a seam the chart wraps across.
bool seam(const FaceData& f, int low, int high)
{
  for (const auto& a : f.sides)
    for (const auto& b : f.sides)
      if (a.side == low && b.side == high && a.edge.IsSame(b.edge)
          && std::abs(a.scale - b.scale) <= 1e-12 * std::abs(a.scale)
          && std::abs(a.offset - b.offset) <= 1e-9 * (1 + std::abs(a.offset)))
        return true;
  return false;
}

bool selfClear(const FaceData& f, long& work)
{
  if (f.kind != FaceData::Spline || !f.rectangle)
    return false;
  auto chart       = cellsOf(f, 0, false, false, false);
  chart.wrapS      = seam(f, 0, 1);
  chart.wrapT      = seam(f, 2, 3);
  const double gap = gapOf(f, f);
  return proveChart(
    chart,
    [](const Chart&, int, int) { return true; },
    [](const Chart& chart, int p, int q) { return p != q && !neighbours(chart, p, q); },
    gap,
    finestOf(f, f, gap),
    0,
    0,
    work);
}

// Two spline faces put side by side across their one shared edge: the first face's side
// along it becomes its last column's s = 1, the second's its first column's s = 0, both run
// with the edge's parameter, rows cut at the breaks of both.
bool besideClear(const FaceData& a, const FaceData& b, const TopoDS_Edge& edge, long& work)
{
  if (a.kind != FaceData::Spline || b.kind != FaceData::Spline || !a.rectangle || !b.rectangle)
    return false;
  const auto sideOn = [&](const FaceData& f) -> std::optional<Side> {
    std::optional<Side> out;
    for (const auto& s : f.sides)
      if (s.edge.IsSame(edge))
      {
        if (out)
          return std::nullopt; // a seam of its own: not a side to put a face beside
        out = s;
      }
    return out;
  };
  const auto sa = sideOn(a), sb = sideOn(b);
  if (!sa || !sb)
    return false;
  // Wrapping along t: both faces closed along the edge, with seams across it.
  const bool wrapT = seam(a, sa->side < 2 ? 2 : 0, sa->side < 2 ? 3 : 1);
  if (wrapT != seam(b, sb->side < 2 ? 2 : 0, sb->side < 2 ? 3 : 1))
    return false;
  // Native coordinates along each side, as the edge's parameter.
  std::vector<double> breaks;
  const auto          place = [&](const FaceData& f, const Side& s) {
    for (double x : s.side < 2 ? f.vBreaks : f.uBreaks)
      breaks.push_back((x - s.offset) / s.scale);
  };
  place(a, *sa);
  place(b, *sb);
  std::sort(breaks.begin(), breaks.end());
  std::vector<double> rows;
  for (double x : breaks)
    if (rows.empty() || x - rows.back() > 1e-12 * (1 + std::abs(x)))
      rows.push_back(x);
    else
      rows.back() = x;
  const auto chartOf = [&](const FaceData& f, const Side& s, int index, bool last) {
    // s along the native direction across the side, t along it.
    const bool transpose = s.side >= 2;
    const bool flipS     = last ? (s.side == 0 || s.side == 2) : (s.side == 1 || s.side == 3);
    const bool flipT     = s.scale < 0;
    auto       chart     = cellsOf(f, index, transpose, flipS, flipT);
    // Cut its rows at every break, in the edge's parameter.
    const auto&         native = s.side < 2 ? f.vBreaks : f.uBreaks;
    std::vector<double> edgeBreaks;
    for (double x : native)
      edgeBreaks.push_back((x - s.offset) / s.scale);
    if (flipT)
      std::reverse(edgeBreaks.begin(), edgeBreaks.end());
    Chart out;
    out.columns = chart.columns;
    out.rows    = int(rows.size()) - 1;
    for (int c = 0; c < chart.columns; ++c)
      for (int r = 0; r + 1 < int(rows.size()); ++r)
      {
        const double from = rows[r], to = rows[r + 1];
        int          own = 0;
        while (own + 2 < int(edgeBreaks.size())
               && edgeBreaks[own + 1] <= from + 1e-12 * (1 + std::abs(from)))
          ++own;
        const double lo = edgeBreaks[own], hi = edgeBreaks[own + 1];
        out.cells.push_back(
          {piece(chart.at(c, own).patch, false, (from - lo) / (hi - lo), (to - lo) / (hi - lo)),
           index,
           {}});
      }
    return out;
  };
  auto  left = chartOf(a, *sa, 0, true), right = chartOf(b, *sb, 1, false);
  Chart chart;
  chart.columns = left.columns + right.columns;
  chart.rows    = left.rows;
  chart.wrapT   = wrapT;
  chart.split   = left.columns;
  chart.cells   = std::move(left.cells);
  chart.cells.insert(chart.cells.end(), right.cells.begin(), right.cells.end());
  const double gap = gapOf(a, b);
  // Only blocks across the glue line and pairs of the two faces: each face's own are its own.
  return proveChart(
    chart,
    [](const Chart& chart, int c, int) { return chart.column(c + 1) == chart.split; },
    [](const Chart& chart, int p, int q) {
      const bool first = p / chart.rows < chart.split, second = q / chart.rows < chart.split;
      if (first == second)
        return false;
      const int  cp = p / chart.rows, cq = q / chart.rows, rp = p % chart.rows, rq = q % chart.rows;
      const int  d = std::abs(rp - rq);
      const bool beside =
        std::abs(cp - cq) == 1 && (d <= 1 || (chart.wrapT && chart.rows - d <= 1));
      return !beside;
    },
    gap,
    finestOf(a, b, gap),
    NEAR_EDGE * gap,
    edgeTube(edge),
    work);
}

// A face beside a planar face it shares edges with: the face leaves the plane monotonically
// from each shared edge and is clear of the plane, or of the planar face's box, elsewhere.
bool planeClear(const FaceData& plane, const FaceData& f, const ShapeSet& shared, long& work)
{
  if (plane.kind != FaceData::Plane || f.kind != FaceData::Spline || !f.rectangle)
    return false;
  // As for its own check: a rational face of few patches costs the check less than a proof.
  if (f.patches.size() < SELF_PATCHES && std::any_of(f.patches.begin(), f.patches.end(), rational))
    return false;
  unsigned sides = 0;
  for (int k = 1; k <= shared.Extent(); ++k)
  {
    if (shared(k).ShapeType() != TopAbs_EDGE)
      continue;
    int found = 0;
    for (const auto& s : f.sides)
      if (s.edge.IsSame(shared(k)))
      {
        sides |= 1u << s.side;
        ++found;
      }
    if (found != 1)
      return false;
  }
  // Every shared vertex on a shared edge.
  for (int k = 1; k <= shared.Extent(); ++k)
  {
    if (shared(k).ShapeType() != TopAbs_VERTEX)
      continue;
    bool onEdge = false;
    for (int e = 1; e <= shared.Extent() && !onEdge; ++e)
      if (shared(e).ShapeType() == TopAbs_EDGE)
        for (TopExp_Explorer it(shared(e), TopAbs_VERTEX); it.More(); it.Next())
          onEdge |= it.Current().IsSame(shared(k));
    if (!onEdge)
      return false;
  }
  const gp_XYZ origin = plane.plane.Location().XYZ(), normal = plane.plane.Axis().Direction().XYZ();
  double       crossing = 1e300;
  for (int k = 1; k <= shared.Extent(); ++k)
    if (shared(k).ShapeType() == TopAbs_EDGE)
      crossing = std::min(crossing, edgeTube(shared(k)));
  const double gap = gapOf(plane, f), finest = finestOf(plane, f, gap);

  struct Job
  {
    Patch    patch;
    unsigned touching;
  };

  std::vector<Job> jobs;
  for (int i = 0; i < f.nu; ++i)
    for (int j = 0; j < f.nv; ++j)
    {
      unsigned touching = 0;
      if (i == 0)
        touching |= 1;
      if (i == f.nu - 1)
        touching |= 2;
      if (j == 0)
        touching |= 4;
      if (j == f.nv - 1)
        touching |= 8;
      jobs.push_back({f.patches[i * f.nv + j], touching & sides});
    }
  while (!jobs.empty())
  {
    auto job = std::move(jobs.back());
    jobs.pop_back();
    if (--work < 0)
      return false;
    auto&               p = job.patch;
    std::vector<double> distance;
    for (const auto& pole : p.poles)
      distance.push_back(normal.Dot(pole - origin));
    const auto [lo, hi] = std::minmax_element(distance.begin(), distance.end());
    if (!job.touching)
    {
      const Box box = p.box();
      if (*lo > gap || *hi < -gap || box.apart(plane.box, gap))
        continue;
      if (box.size() <= finest)
        return false;
      auto [one, two] = halve(p);
      jobs.push_back({std::move(one), 0});
      jobs.push_back({std::move(two), 0});
      continue;
    }
    const bool single = (job.touching & (job.touching - 1)) == 0;
    if (single)
    {
      // Turned so that the shared side is t = 0.
      const int side = job.touching == 1 ? 0 : job.touching == 2 ? 1 : job.touching == 4 ? 2 : 3;
      Patch     q    = side < 2 ? transposed(p) : p;
      if (side == 1 || side == 3)
        q = flipped(q, false);
      const auto d = derivatives(q, work);
      // Distance d = H / w with H = sum w_ij d_ij; its t derivative's numerator is
      // H_t w - H w_t, here from the same homogeneous algebra.
      std::vector<double> h, w;
      const bool          weighted = rational(q);
      for (int i = 0; i <= q.n; ++i)
        for (int j = 0; j <= q.m; ++j)
        {
          w.push_back(weighted ? q.weight(i, j) : 1.0);
          h.push_back(w.back() * normal.Dot(q.pole(i, j) - origin));
        }
      std::vector<double> ht((q.n + 1) * q.m), wt((q.n + 1) * q.m);
      for (int i = 0; i <= q.n; ++i)
        for (int j = 0; j < q.m; ++j)
        {
          ht[i * q.m + j] = (h[i * (q.m + 1) + j + 1] - h[i * (q.m + 1) + j]) * q.m;
          wt[i * q.m + j] = (w[i * (q.m + 1) + j + 1] - w[i * (q.m + 1) + j]) * q.m;
        }
      auto rise    = GeomAlgoExt::product(ht, q.n, q.m - 1, w, q.n, q.m),
           fall    = GeomAlgoExt::product(h, q.n, q.m, wt, q.n, q.m - 1);
      double least = 1e300, most = -1e300;
      for (size_t k = 0; k < rise.size(); ++k)
      {
        least = std::min(least, rise[k] - fall[k]);
        most  = std::max(most, rise[k] - fall[k]);
      }
      const double slope = least > 0 ? least : most < 0 ? -most : 0;
      double       edge  = 0;
      for (int i = 0; i <= q.n; ++i)
        edge = std::max(edge, std::abs(normal.Dot(q.pole(i, 0) - origin)));
      if (slope > 0)
      {
        // Within gap of the plane only for t up to (gap + edge) wMax^2 / slope, so within
        // that times the speed along t of the shared edge; on it only for t up to
        // edge wMax^2 / slope.
        double speed = 0;
        for (const auto& v : d.v)
          speed = std::max(speed, v.Modulus());
        speed /= d.wMin * d.wMin;
        const double reach = (gap + edge) * d.wMax * d.wMax / slope * speed;
        const double cross = edge * d.wMax * d.wMax / slope * speed;
        if (reach <= NEAR_EDGE * gap && cross <= crossing)
          continue;
      }
    }
    // Halved: along the shared side keeps it, across it leaves one part off it.
    if (p.box().size() <= finest)
      return false;
    const bool alongU       = spreadsAlongU(p);
    auto [one, two]         = split(p, alongU, .5);
    const unsigned keepLow  = alongU ? (job.touching & ~2u) : (job.touching & ~8u);
    const unsigned keepHigh = alongU ? (job.touching & ~1u) : (job.touching & ~4u);
    jobs.push_back({std::move(one), keepLow});
    jobs.push_back({std::move(two), keepHigh});
  }
  return true;
}

bool apartFaces(const FaceData& a, const FaceData& b, long& work)
{
  // An Unknown face's box bounds nothing (it may be empty, or what was read before a failure).
  if (a.kind == FaceData::Unknown || b.kind == FaceData::Unknown)
    return false;
  const double gap = gapOf(a, b);
  if (a.box.apart(b.box, gap))
    return true;
  if (a.kind == FaceData::Hull || b.kind == FaceData::Hull)
    return false;
  if (a.kind == FaceData::Plane && b.kind == FaceData::Plane)
    return false;
  if (a.kind == FaceData::Plane || b.kind == FaceData::Plane)
  {
    // A planar face is only its box here; the other face is halved until clear of it.
    const auto&  flat   = a.kind == FaceData::Plane ? a : b;
    const auto&  curved = a.kind == FaceData::Plane ? b : a;
    const double finest = finestOf(a, b, gap);
    for (const auto& p : curved.patches)
    {
      std::vector<Patch> stack{p};
      while (!stack.empty())
      {
        auto c = std::move(stack.back());
        stack.pop_back();
        if (--work < 0)
          return false;
        const Box box = c.box();
        if (box.apart(flat.box, gap))
          continue;
        if (box.size() <= finest)
          return false;
        auto [one, two] = halve(c);
        stack.push_back(std::move(one));
        stack.push_back(std::move(two));
      }
    }
    return true;
  }
  const double finest = finestOf(a, b, gap);
  for (const auto& p : a.patches)
  {
    if (p.box().apart(b.box, gap))
      continue;
    for (const auto& q : b.patches)
      if (!apart(p, q, gap, finest, work))
        return false;
  }
  return true;
}

} // namespace

std::vector<TopoDS_Shape> unprovenParts(const TopoDS_Shape& shape, SeparationCount* count)
{
  ShapeSet faces, edges, faceEdges;
  TopExp::MapShapes(shape, TopAbs_FACE, faces);
  TopExp::MapShapes(shape, TopAbs_EDGE, edges);
  for (int i = 1; i <= faces.Extent(); ++i)
    TopExp::MapShapes(faces(i), TopAbs_EDGE, faceEdges);
  // Edges or vertices of no face: what only faces are proven for leaves them to the check.
  ShapeSet vertices, edgeVertices;
  TopExp::MapShapes(shape, TopAbs_VERTEX, vertices);
  for (int i = 1; i <= edges.Extent(); ++i)
    TopExp::MapShapes(edges(i), TopAbs_VERTEX, edgeVertices);
  if (edges.Extent() != faceEdges.Extent() || vertices.Extent() != edgeVertices.Extent())
  {
    if (count)
      *count = {};
    return {shape};
  }
  const int n = faces.Extent();
  // Worth trying only where the check costs much: spline faces of many spans, or a spline face
  // beside another face, whose intersection the check traces along their whole edge.
  {
    const auto spans = [](const TopoDS_Face& face) {
      auto surface = BRep_Tool::Surface(face);
      while (!surface.IsNull() && surface->IsKind(STANDARD_TYPE(Geom_RectangularTrimmedSurface)))
        surface = Handle(Geom_RectangularTrimmedSurface)::DownCast(surface)->BasisSurface();
      const auto spline = Handle(Geom_BSplineSurface)::DownCast(surface);
      if (!spline.IsNull())
        return (spline->NbUKnots() - 1) * (spline->NbVKnots() - 1);
      return Handle(Geom_BezierSurface)::DownCast(surface).IsNull() ? 0 : 1;
    };
    bool worth = false;
    for (int i = 1; i <= n && !worth; ++i)
      worth = spans(TopoDS::Face(faces(i))) >= int(SELF_PATCHES);
    if (!worth)
    {
      NCollection_IndexedDataMap<TopoDS_Shape,
                                 NCollection_List<TopoDS_Shape>,
                                 TopTools_ShapeMapHasher>
        owners;
      TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, owners);
      for (int e = 1; e <= owners.Extent() && !worth; ++e)
      {
        ShapeSet beside;
        bool     spline = false;
        for (const auto& face : owners(e))
        {
          beside.Add(face);
          spline |= spans(TopoDS::Face(face)) > 0;
        }
        worth = spline && beside.Extent() >= 2;
      }
    }
    if (!worth)
      worth = analyticWorth(shape);
    if (!worth)
    {
      if (count)
        *count = {};
      return {shape};
    }
  }
  std::vector<FaceData> data;
  for (int i = 1; i <= n; ++i)
    data.push_back(faceData(TopoDS::Face(faces(i))));
  long work = WORK;
  // Faces a part of the check still needs, joined where a pair between them is not proven.
  std::vector<char> needed(n, 0);
  std::vector<int>  part(n);
  for (int i = 0; i < n; ++i)
    part[i] = i;
  const auto root = [&](int i) {
    while (part[i] != i)
      i = part[i] = part[part[i]];
    return i;
  };
  SeparationCount counted;
  for (int i = 0; i < n; ++i)
  {
    ++counted.pairs;
    // A rational face of few patches costs the check little to meet with itself, and a proof
    // much, its derivatives being of twice its degrees: the check does that.
    const bool weighted = std::any_of(data[i].patches.begin(), data[i].patches.end(), rational);
    if (work > 0 && (!weighted || data[i].patches.size() >= SELF_PATCHES)
        && selfClear(data[i], work))
    {
      ++counted.proven;
      counted.provenPairs.push_back({i, i});
    }
    else
      needed[i] = 1;
    for (int j = i + 1; j < n; ++j)
    {
      ++counted.pairs;
      ShapeSet shared;
      for (int k = 1; k <= data[i].edges.Extent(); ++k)
        if (data[j].edges.Contains(data[i].edges(k)))
          shared.Add(data[i].edges(k));
      for (int k = 1; k <= data[i].vertices.Extent(); ++k)
        if (data[j].vertices.Contains(data[i].vertices(k)))
          shared.Add(data[i].vertices(k));
      int         sharedEdges = 0, ends = 0;
      TopoDS_Edge edge;
      for (int k = 1; k <= shared.Extent(); ++k)
        if (shared(k).ShapeType() == TopAbs_EDGE)
        {
          ++sharedEdges;
          edge = TopoDS::Edge(shared(k));
        }
      for (int k = 1; k <= shared.Extent() && sharedEdges == 1; ++k)
        if (shared(k).ShapeType() == TopAbs_VERTEX)
          for (TopExp_Explorer it(edge, TopAbs_VERTEX); it.More(); it.Next())
            if (it.Current().IsSame(shared(k)))
            {
              ++ends;
              break;
            }
      bool proven = false;
      if (work > 0)
      {
        if (shared.IsEmpty())
          proven = apartFaces(data[i], data[j], work);
        else if (data[i].kind == FaceData::Plane)
          proven = planeClear(data[i], data[j], shared, work);
        else if (data[j].kind == FaceData::Plane)
          proven = planeClear(data[j], data[i], shared, work);
        else if (sharedEdges == 1 && ends == shared.Extent() - 1)
          proven = besideClear(data[i], data[j], edge, work);
        if (!proven && !shared.IsEmpty())
          proven = analyticBeside(data[i].face, data[j].face, shared);
      }
      if (proven)
      {
        ++counted.proven;
        counted.provenPairs.push_back({i, j});
      }
      else
      {
        needed[i] = needed[j] = 1;
        part[root(i)]         = root(j);
      }
    }
  }
  if (count)
    *count = counted;
  std::vector<std::vector<int>> parts(n);
  for (int i = 0; i < n; ++i)
    if (needed[i])
      parts[root(i)].push_back(i);
  std::vector<TopoDS_Shape> out;
  for (const auto& faces_ : parts)
  {
    if (faces_.empty())
      continue;
    if (int(faces_.size()) == n)
      return {shape};
    if (faces_.size() == 1)
    {
      out.push_back(faces(faces_[0] + 1));
      continue;
    }
    TopoDS_Compound compound;
    BRep_Builder    builder;
    builder.MakeCompound(compound);
    for (int i : faces_)
      builder.Add(compound, faces(i + 1));
    out.push_back(compound);
  }
  return out;
}

namespace
{
// The self-interference check of a shape, all of it at once.
BOPAlgo_ArgumentAnalyzer analyzed(const TopoDS_Shape& shape)
{
  BOPAlgo_ArgumentAnalyzer analyzer;
  analyzer.SetShape1(shape);
  analyzer.ArgumentTypeMode() = false;
  analyzer.SmallEdgeMode()    = false;
  analyzer.SelfInterMode()    = true;
  analyzer.Perform();
  return analyzer;
}
} // namespace

std::vector<double> faceSeparationAudit(OcctKernel& kernel, uint32_t id)
{
  const auto&     shape = get(kernel, id);
  SeparationCount count;
  unprovenParts(shape, &count);
  ShapeSet faces;
  TopExp::MapShapes(shape, TopAbs_FACE, faces);
  int crossing = 0, aborted = 0;
  for (const auto& [i, j] : count.provenPairs)
  {
    const auto&  a      = faces(i + 1);
    const auto&  b      = faces(j + 1);
    TopoDS_Shape tested = a;
    if (i != j)
    {
      TopoDS_Compound compound;
      BRep_Builder    builder;
      builder.MakeCompound(compound);
      builder.Add(compound, a);
      builder.Add(compound, b);
      tested = compound;
    }
    const auto analyzer = analyzed(tested);
    ShapeSet   ofA, ofB;
    TopExp::MapShapes(a, ofA);
    TopExp::MapShapes(b, ofB);
    // Of the pair, what interferes across it: a sub-shape of one face only with one of the
    // other only; of a face by itself, anything.
    const auto only = [](const ShapeSet& in, const ShapeSet& out, const TopoDS_Shape& s) {
      return in.Contains(s) && !out.Contains(s);
    };
    bool across = false, failed = false;
    for (const auto& result : analyzer.GetCheckResult())
    {
      if (result.GetCheckStatus() == BOPAlgo_OperationAborted)
        failed = true;
      if (result.GetCheckStatus() != BOPAlgo_SelfIntersect)
        continue;
      const auto& faulty = result.GetFaultyShapes1();
      if (i == j)
        across = true;
      else if (faulty.Size() == 2)
      {
        const auto& s1 = faulty.First();
        const auto& s2 = faulty.Last();
        across |=
          (only(ofA, ofB, s1) && only(ofB, ofA, s2)) || (only(ofB, ofA, s1) && only(ofA, ofB, s2));
      }
    }
    crossing += across;
    aborted += failed;
  }
  // The verdict with the pre-pass (what the operations ask) and without it.
  return {double(count.provenPairs.size()),
          double(crossing),
          double(aborted),
          double(intersectsItself(shape)),
          double(analyzed(shape).HasFaulty())};
}

std::vector<double> faceSeparation(OcctKernel& kernel, uint32_t id)
{
  const auto&     shape = get(kernel, id);
  SeparationCount count;
  ShapeSet        faces, left;
  TopExp::MapShapes(shape, TopAbs_FACE, faces);
  const auto parts = unprovenParts(shape, &count);
  for (const auto& part : parts)
    TopExp::MapShapes(part, TopAbs_FACE, left);
  return {double(count.pairs),
          double(count.proven),
          double(faces.Extent()),
          double(left.Extent()),
          double(parts.size())};
}

} // namespace BRepExt
