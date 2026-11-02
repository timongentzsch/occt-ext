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

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepTools.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Line.hxx>
#include <Geom_Surface.hxx>
#include <Geom_SurfaceOfRevolution.hxx>
#include <Geom_BSplineSurface.hxx>
#include <IntCurvesFace_Intersector.hxx>
#include <NCollection_HArray1.hxx>
#include <NCollection_IndexedMap.hxx>
#include <Precision.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <gp_Ax2.hxx>
#include <gp_Lin.hxx>
#include <gp_Sphere.hxx>
#include <gp_Torus.hxx>
#include <gp_Trsf.hxx>
#include <map>
#include <memory>
#include <set>
#include <vector>

// The outline of a body seen along a direction: the boundary of its shadow. It is worked out
// without OpenCascade's hidden-line removal, which loses or doubles curves here (a box seen
// along a side, a torus seen side-on or steeply). Every edge of the body, and on every curved
// face the curve where the surface turns away from the eye (n . v = 0), is a candidate.
//
// Whether a candidate is outline can change only at certain places of it, the events: where its
// image crosses or touches another candidate's or its own, where another candidate ends on it or
// stops having the same image, where its own image turns back (a fold or a cusp), and where a
// face beside an edge turns edge-on. The events are found exactly (Newton on the two curves in
// the image), every candidate through an event's point is cut there, and between two events a
// candidate is asked once: it is kept where it has the body on one side and nothing on the other
// (a ray to each side against the exact faces) and its faces lie to one side of it in the image.
// Of kept stretches with the same image the one nearest the eye stays. The pieces therefore end
// on each other exactly, and nothing depends on a number of samples.
//
// Lengths are measured against the size of the body, and the work is done with the body moved
// to the origin, so that neither its size nor its place changes the answer.

namespace BRepExt
{
namespace
{

// A point of the image: what is seen of a point in space.
struct Image
{
  double x, y;
};

struct View
{
  gp_Pnt origin;
  gp_Vec x, y, n; // n towards the eye

  Image of(const gp_Pnt& p) const
  {
    const gp_Vec d(origin, p);
    return {d.Dot(x), d.Dot(y)};
  }

  Image of(const gp_Vec& d) const { return {d.Dot(x), d.Dot(y)}; }
};

double gap(Image a, Image b)
{
  return std::hypot(a.x - b.x, a.y - b.y);
}

// The lengths things are compared to, all of them shares of the body's size (the diagonal of
// its box).
struct Lengths
{
  double span;
  // Two candidates have the same image where they are this near in it: the accuracy the
  // silhouettes are traced to, times ten.
  double sameImage;
  // A silhouette that is not an exact curve lies this near the true one.
  double traced;
  // Two points are at the same depth to this.
  double sameDepth;
  // How far to each side of a candidate the rays are cast at most: far enough to clear the
  // tolerance the rays hit faces with, near enough to stay beside the candidate.
  double aside;
  // The tolerance the rays hit faces with.
  double rays;
  // A polyline of a candidate's image leaves it by no more than this (and by no more than a
  // twentieth of each of its pieces).
  double flat;
  // Below this two points of the image are one: the rounding of the arithmetic.
  double point;
  // Two events this near in the image are one event: twenty times the rounding, a hundredth
  // of what makes two images the same.
  double one;

  explicit Lengths(double size)
      : span(size),
        sameImage(2e-7 * size),
        traced(2e-8 * size),
        sameDepth(1e-8 * size),
        aside(1e-6 * size),
        rays(1e-8 * size),
        flat(2e-3 * size),
        point(1e-10 * size),
        one(2e-9 * size)
  {
  }
};

// Two candidates run the same way where the sine of the angle between their images is under this.
constexpr double sameWay = 1e-4;
// A face beside an edge is seen edge-on where the sine of its angle to the view is under this
// (the faces' normals are exact to about 1e-12; this is far above that and far below any angle
// a face is modelled at).
constexpr double edgeOn = BRepExt::tolerance::angularResolution;

// Looking along a straight prism, its two caps have the same shadow and its walls
// are edge-on. Recognise that exactly for line and full-circle boundaries; no rays
// or arrangement are needed. Other bodies take the general event path below.
TopoDS_Shape axialPrism(const TopoDS_Shape& shape, const View& view, const Lengths& len)
{
  std::vector<TopoDS_Face> caps;
  for (TopExp_Explorer it(shape, TopAbs_FACE); it.More(); it.Next())
  {
    const auto                face = TopoDS::Face(it.Current());
    const BRepAdaptor_Surface surface(face);
    if (surface.GetType() == GeomAbs_Plane)
    {
      const double facing = std::abs(gp_Vec(surface.Plane().Axis().Direction()).Dot(view.n));
      if (facing >= 1 - 1e-14)
        caps.push_back(face);
      else if (facing > BRepExt::tolerance::fineAngle)
        return {};
    }
    else if (surface.GetType() != GeomAbs_Cylinder
             || std::abs(gp_Vec(surface.Cylinder().Axis().Direction()).Dot(view.n)) < 1 - 1e-14)
    {
      return {};
    }
  }
  if (caps.size() != 2)
    return {};
  auto edgesOf = [](const TopoDS_Face& face) {
    std::vector<TopoDS_Edge> edges;
    for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next())
      edges.push_back(TopoDS::Edge(it.Current()));
    return edges;
  };
  const auto a = edgesOf(caps[0]), b = edgesOf(caps[1]);
  if (a.empty() || a.size() != b.size())
    return {};
  auto same = [&](const TopoDS_Edge& one, const TopoDS_Edge& two) {
    const BRepAdaptor_Curve first(one), second(two);
    if (first.GetType() != second.GetType())
      return false;
    if (first.GetType() == GeomAbs_Circle)
    {
      if (std::abs(first.LastParameter() - first.FirstParameter() - 2 * M_PI)
            > BRepExt::tolerance::fineAngle
          || std::abs(second.LastParameter() - second.FirstParameter() - 2 * M_PI)
               > BRepExt::tolerance::fineAngle)
        return false;
      return std::abs(first.Circle().Radius() - second.Circle().Radius()) <= len.point
             && gap(view.of(first.Circle().Location()), view.of(second.Circle().Location()))
                  <= len.point;
    }
    if (first.GetType() != GeomAbs_Line)
      return false;
    const Image a0 = view.of(first.Value(first.FirstParameter())),
                a1 = view.of(first.Value(first.LastParameter()));
    const Image b0 = view.of(second.Value(second.FirstParameter())),
                b1 = view.of(second.Value(second.LastParameter()));
    return (gap(a0, b0) <= len.point && gap(a1, b1) <= len.point)
           || (gap(a0, b1) <= len.point && gap(a1, b0) <= len.point);
  };
  std::vector<bool> used(b.size(), false);
  for (const auto& edge : a)
  {
    bool found = false;
    for (size_t j = 0; j < b.size() && !found; ++j)
      if (!used[j] && same(edge, b[j]))
      {
        used[j] = true;
        found   = true;
      }
    if (!found)
      return {};
  }
  const auto depth = [&](const TopoDS_Face& face) {
    return gp_Vec(view.origin, BRepAdaptor_Surface(face).Plane().Location()).Dot(view.n);
  };
  const auto&     front = depth(caps[0]) >= depth(caps[1]) ? a : b;
  BRep_Builder    builder;
  TopoDS_Compound result;
  builder.MakeCompound(result);
  for (const auto& edge : front)
    builder.Add(result, edge);
  return result;
}

struct OnFace
{
  Handle(ShapeAnalysis_Surface)            surface;
  std::shared_ptr<BRepTopAdaptor_FClass2d> inside;
};

// A face beside an edge: where the edge runs in its parameters, and whether the face lies to
// the left of the edge there (+1) or to its right (-1), seen against the surface's normal.
struct Beside
{
  Handle(Geom_Surface) surface;
  Handle(Geom2d_Curve) trace;
  double               left;
  bool                 flat;
};

struct Candidate
{
  Handle(Geom_Curve)  curve;
  double              a, b;
  int                 face;   // the face a silhouette has to lie on, or -1 for an edge
  std::vector<Beside> beside; // the faces of an edge
  bool                closed; // it ends where it starts: its ends are no events
  std::vector<double> at;     // the parameters of a polyline of its image
  std::vector<Image>  seen;   // its image there
  std::vector<double> slack;  // how far the image leaves each piece of the polyline, at most
  double              low[2], high[2]; // the box of the polyline, widened by the slack
  std::vector<double> events;          // where it may start or stop being outline
  std::vector<char>   state;           // between each two events: 1 kept, 0 not
};

Image imageOf(const Candidate& c, const View& view, double t)
{
  return view.of(c.curve->Value(t));
}

// The distance of q from the piece between a and b, and where along it (0 to 1) the nearest
// point lies.
double fromPiece(Image q, Image a, Image b, double& share)
{
  const double dx = b.x - a.x, dy = b.y - a.y, all = dx * dx + dy * dy;
  share = all > 0 ? std::clamp(((q.x - a.x) * dx + (q.y - a.y) * dy) / all, 0.0, 1.0) : 0.0;
  return std::hypot(q.x - a.x - share * dx, q.y - a.y - share * dy);
}

// The share of the image's speed in the curve's: 0 where the curve is seen end-on.
double still(const Candidate& c, const View& view, double t)
{
  gp_Pnt p;
  gp_Vec d;
  c.curve->D1(t, p, d);
  const Image  way = view.of(d);
  const double all = d.SquareMagnitude();
  return all < 1e-300 ? 1.0 : (way.x * way.x + way.y * way.y) / all;
}

template <typename F>
void roots(F f, double a, double fa, double b, double fb, int depth, std::vector<double>& out);

// The polyline of a candidate's image, as fine as the image is curved, and the places where
// the image turns back on itself. False for a candidate seen end-on: a point bounds nothing.
bool prepare(Candidate& c, const View& view, const Lengths& len, std::vector<double>& folds)
{
  const bool straight = c.curve->IsKind(STANDARD_TYPE(Geom_Line));
  int        pieces   = straight ? 1 : 16;
  if (const auto spline = Handle(Geom_BSplineCurve)::DownCast(c.curve); !spline.IsNull())
    pieces = std::max(pieces, 2 * spline->NbKnots());
  c.at.assign(1, c.a);
  c.seen.assign(1, imageOf(c, view, c.a));
  std::function<void(double, Image, double, Image, int)> add =
    [&](double t0, Image q0, double t1, Image q1, int depth) {
      double away = 0;
      Image  middle{};
      if (!straight)
        for (int k = 1; k <= 3; ++k)
        {
          const Image q = imageOf(c, view, t0 + (t1 - t0) * k / 4);
          double      share;
          away = std::max(away, fromPiece(q, q0, q1, share));
          if (k == 2)
            middle = q;
        }
      if (depth < 16 && away > std::min(len.flat, 0.05 * gap(q0, q1)) && away > len.point)
      {
        const double tm = (t0 + t1) / 2;
        add(t0, q0, tm, middle, depth + 1);
        add(tm, middle, t1, q1, depth + 1);
        return;
      }
      c.at.push_back(t1);
      c.seen.push_back(q1);
      c.slack.push_back(2 * away + len.point);
    };
  for (int i = 0; i < pieces; ++i)
  {
    const double t1 = i + 1 == pieces ? c.b : c.a + (c.b - c.a) * (i + 1) / pieces;
    const double t0 = c.at.back();
    const Image  q0 = c.seen.back();
    add(t0, q0, t1, imageOf(c, view, t1), 0);
  }
  c.low[0] = c.low[1] = std::numeric_limits<double>::infinity();
  c.high[0] = c.high[1] = -std::numeric_limits<double>::infinity();
  double length = 0, widest = 0;
  for (size_t i = 0; i < c.seen.size(); ++i)
  {
    if (i)
    {
      length += gap(c.seen[i], c.seen[i - 1]);
      widest = std::max(widest, c.slack[i - 1]);
    }
    c.low[0]  = std::min(c.low[0], c.seen[i].x);
    c.low[1]  = std::min(c.low[1], c.seen[i].y);
    c.high[0] = std::max(c.high[0], c.seen[i].x);
    c.high[1] = std::max(c.high[1], c.seen[i].y);
  }
  for (int k = 0; k < 2; ++k)
  {
    c.low[k] -= widest;
    c.high[k] += widest;
  }
  if (length <= 1e-7 * len.span)
    return false;
  // A circle projects to an ellipse unless its plane is exactly edge-on. An
  // ellipse never stops or turns back, so it needs no search for turning points.
  const auto circle  = Handle(Geom_Circle)::DownCast(c.curve);
  const bool canFold = !straight
                       && (circle.IsNull()
                           || std::abs(gp_Vec(circle->Axis().Direction()).Dot(view.n))
                                < BRepExt::tolerance::fineAngle);
  // Folds: a curve seen edge-on in part (a circle from its side) runs out and back along the
  // same image, and turns where its image stands still; a silhouette's image has cusps.
  if (canFold)
  {
    std::vector<double> value(c.at.size());
    for (size_t i = 0; i < c.at.size(); ++i)
      value[i] = still(c, view, c.at[i]);
    for (size_t i = 1; i + 1 < c.at.size(); ++i)
    {
      if (value[i] > value[i - 1] || value[i] > value[i + 1]
          || value[i] >= std::max(value[i - 1], value[i + 1]) - 1e-12)
        continue;
      double lo = c.at[i - 1], hi = c.at[i + 1];
      for (int k = 0; k < 100; ++k)
      {
        const double m1 = lo + (hi - lo) / 3, m2 = hi - (hi - lo) / 3;
        (still(c, view, m1) < still(c, view, m2) ? hi : lo) =
          still(c, view, m1) < still(c, view, m2) ? m2 : m1;
      }
      if (still(c, view, (lo + hi) / 2) >= 1e-12)
        continue;
      // The image runs one way before the fold and the other way after it: the fold is
      // where its direction changes sign, which is found more sharply than a least value.
      auto way = [&](double t) {
        gp_Pnt p;
        gp_Vec d;
        c.curve->D1(t, p, d);
        return view.of(d);
      };
      const Image before = way(c.at[i - 1]), after = way(c.at[i + 1]);
      double      fold = (lo + hi) / 2;
      if (before.x * after.x + before.y * after.y < 0)
      {
        double from = c.at[i - 1], to = c.at[i + 1];
        for (int k = 0; k < 80; ++k)
        {
          const double mid                                        = (from + to) / 2;
          const Image  here                                       = way(mid);
          (here.x * before.x + here.y * before.y > 0 ? from : to) = mid;
        }
        fold = (from + to) / 2;
      }
      folds.push_back(fold);
    }
  }
  if (canFold)
  {
    for (int component = 0; component < 2; ++component)
    {
      auto derivative = [&](double t) {
        gp_Pnt p;
        gp_Vec d;
        c.curve->D1(t, p, d);
        const Image q = view.of(d);
        return component == 0 ? q.x : q.y;
      };
      std::vector<double> found;
      for (size_t i = 0; i + 1 < c.at.size(); ++i)
        roots(derivative,
              c.at[i],
              derivative(c.at[i]),
              c.at[i + 1],
              derivative(c.at[i + 1]),
              0,
              found);
      for (const double t : found)
        if (still(c, view, t) < 1e-12)
          folds.push_back(t);
    }
  }
  return true;
}

// The parameter near `start` at which the image of `c` is nearest to `q`.
double foot(const Candidate& c, const View& view, Image q, double start, double low, double high)
{
  double     s     = std::clamp(start, low, high);
  const bool whole = low <= c.a && high >= c.b;
  for (int k = 0; k < 30; ++k)
  {
    gp_Pnt p;
    gp_Vec d1;
    c.curve->D1(s, p, d1);
    const Image  at = view.of(p), t = view.of(d1);
    const double tt = t.x * t.x + t.y * t.y;
    if (tt < 1e-300)
      break;
    const double move = ((q.x - at.x) * t.x + (q.y - at.y) * t.y) / tt;
    double       next = s + move;
    if (whole && c.closed && (next < c.a || next > c.b))
      next = c.a + std::fmod(std::fmod(next - c.a, c.b - c.a) + (c.b - c.a), c.b - c.a);
    next               = std::clamp(next, low, high);
    const bool settled = std::abs(next - s) <= 1e-15 * std::max(1.0, std::abs(s));
    s                  = next;
    if (settled)
      break;
  }
  return s;
}

double foot(const Candidate& c, const View& view, Image q, double start)
{
  return foot(c, view, q, start, c.a, c.b);
}

// The parameters at which the image of `c` is nearest to `q`, locally, and within `reach`.
std::vector<double> nearest(const Candidate& c, const View& view, Image q, double reach)
{
  std::vector<double> found;
  if (q.x < c.low[0] - reach || q.x > c.high[0] + reach || q.y < c.low[1] - reach
      || q.y > c.high[1] + reach)
    return found;
  for (size_t k = 0; k + 1 < c.at.size(); ++k)
  {
    double share;
    if (fromPiece(q, c.seen[k], c.seen[k + 1], share) > reach + c.slack[k])
      continue;
    const double s = foot(c, view, q, c.at[k] + (c.at[k + 1] - c.at[k]) * share);
    if (gap(imageOf(c, view, s), q) > reach)
      continue;
    bool known = false;
    for (const double other : found)
      known |= std::abs(other - s) <= BRepExt::tolerance::parameterMerge * (c.b - c.a);
    if (!known)
      found.push_back(s);
  }
  return found;
}

// The roots of f between a and b: its changes of sign, looked for more finely where it bends.
template <typename F>
void roots(F f, double a, double fa, double b, double fb, int depth, std::vector<double>& out)
{
  const double m = (a + b) / 2, fm = f(m);
  if (std::isnan(fa) || std::isnan(fb) || std::isnan(fm))
    return;
  const bool left = (fa >= 0) != (fm >= 0), right = (fm >= 0) != (fb >= 0);
  if (!left && !right)
  {
    const double least = std::min({std::abs(fa), std::abs(fm), std::abs(fb)});
    if (depth < 3 || (depth < 10 && std::abs(fm - (fa + fb) / 2) > 0.25 * least))
    {
      roots(f, a, fa, m, fm, depth + 1, out);
      roots(f, m, fm, b, fb, depth + 1, out);
    }
    return;
  }
  if (left && right && depth < 10)
  {
    roots(f, a, fa, m, fm, depth + 1, out);
    roots(f, m, fm, b, fb, depth + 1, out);
    return;
  }
  double lo = left ? a : m, hi = left ? m : b, flo = left ? fa : fm;
  for (int i = 0; i < 80 && hi - lo > 1e-16 * std::max(1.0, std::abs(lo)); ++i)
  {
    const double mid = (lo + hi) / 2, fmid = f(mid);
    if ((fmid >= 0) == (flo >= 0))
    {
      lo  = mid;
      flo = fmid;
    }
    else
    {
      hi = mid;
    }
  }
  out.push_back((lo + hi) / 2);
}

// A root of f between a and b, where it changes sign.
template <typename F>
double root(F f, double a, double b)
{
  double fa = f(a);
  for (int i = 0; i < 80; ++i)
  {
    const double m = (a + b) / 2, fm = f(m);
    if ((fm >= 0) == (fa >= 0))
    {
      a  = m;
      fa = fm;
    }
    else
    {
      b = m;
    }
  }
  return (a + b) / 2;
}

// How far the image of a at t is from the image of b, taken near the parameter `beside` of b
// (which is moved to the nearest one). Past an end of b, where b's image does not go on, it is
// no distance at all: the two do not have the same image there.
double apart(const Candidate& a,
             double           t,
             const Candidate& b,
             double&          beside,
             const View&      view,
             double           point,
             double           low,
             double           high)
{
  const Image q = imageOf(a, view, t);
  beside        = foot(b, view, q, beside, low, high);
  gp_Pnt p;
  gp_Vec d;
  b.curve->D1(beside, p, d);
  const Image on = view.of(p), way = view.of(d);
  if (!b.closed && (beside <= b.a || beside >= b.b))
  {
    const double speed = std::hypot(way.x, way.y);
    const double along = speed > 0 ? ((q.x - on.x) * way.x + (q.y - on.y) * way.y) / speed : 0;
    if ((beside <= b.a ? -along : along) > point)
      return std::numeric_limits<double>::infinity();
  }
  return gap(q, on);
}

struct Scene
{
  View                   view;
  Lengths                len{1};
  std::vector<Candidate> candidates;
  std::vector<OnFace>    faces;

  struct RayFace
  {
    Handle(IntCurvesFace_Intersector) ray;
    double                            low[2], high[2];
  };

  std::vector<RayFace> rays;
  size_t               lastHit = 0;
  std::vector<Image>   points; // the events' places in the image
  std::map<std::pair<long long, long long>, std::vector<size_t>>
    known; // ... by cell, to take each once

  void load(const TopoDS_Shape& shape)
  {
    for (TopExp_Explorer it(shape, TopAbs_FACE); it.More(); it.Next())
    {
      const auto face = TopoDS::Face(it.Current());
      RayFace    item;
      // A modelling tolerance may be much wider than the shadow feature. Rays
      // classify against the trims themselves, not a tolerance band round them.
      item.ray = new IntCurvesFace_Intersector(face, len.rays, true, false);
      Bnd_Box box;
      BRepBndLib::Add(face, box, false);
      double x0, y0, z0, x1, y1, z1;
      box.Get(x0, y0, z0, x1, y1, z1);
      item.low[0] = item.low[1] = 1e300;
      item.high[0] = item.high[1] = -1e300;
      for (const double x : {x0, x1})
        for (const double y : {y0, y1})
          for (const double z : {z0, z1})
          {
            const Image q = view.of(gp_Pnt(x, y, z));
            item.low[0]   = std::min(item.low[0], q.x);
            item.high[0]  = std::max(item.high[0], q.x);
            item.low[1]   = std::min(item.low[1], q.y);
            item.high[1]  = std::max(item.high[1], q.y);
          }
      rays.push_back(item);
    }
  }

  bool occupied(const gp_Pnt& p)
  {
    const Image  q = view.of(p);
    const gp_Lin line(p, gp_Dir(view.n));
    for (size_t step = 0; step < rays.size(); ++step)
    {
      const size_t   index = (lastHit + step) % rays.size();
      const RayFace& face  = rays[index];
      if (q.x < face.low[0] - len.rays || q.x > face.high[0] + len.rays
          || q.y < face.low[1] - len.rays || q.y > face.high[1] + len.rays)
        continue;
      face.ray->Perform(line, -2 * len.span, 2 * len.span);
      if (face.ray->NbPnt() > 0)
      {
        lastHit = index;
        return true;
      }
    }
    return false;
  }

  bool onFace(const Candidate& c, const gp_Pnt& p)
  {
    if (c.face < 0)
      return true;
    const auto&    face = faces[c.face];
    const gp_Pnt2d uv   = face.surface->ValueOfUV(p, len.rays);
    return face.inside->Perform(uv) != TopAbs_OUT;
  }

  // An event at a place of the image. The same event is found more than once, by different
  // arithmetic (two curves that touch where a face turns edge-on): places within `len.one`
  // of an earlier one are that one. Every candidate is then cut at the same places.
  void event(Image q)
  {
    const long long x = std::llround(q.x / len.one), y = std::llround(q.y / len.one);
    for (long long i = x - 1; i <= x + 1; ++i)
      for (long long j = y - 1; j <= y + 1; ++j)
      {
        const auto cell = known.find({i, j});
        if (cell == known.end())
          continue;
        for (const size_t k : cell->second)
          if (gap(points[k], q) <= len.one)
            return;
      }
    known[{x, y}].push_back(points.size());
    points.push_back(q);
  }

  void event(size_t index, double t)
  {
    Candidate&  c = candidates[index];
    const Image q = imageOf(c, view, t);
    for (const double known : c.events)
      if (gap(imageOf(c, view, known), q) <= len.one
          && gap(imageOf(c, view, (known + t) / 2), q) <= len.one)
        return;
    c.events.push_back(t);
    event(q);
  }

  // Whether candidate `index` bounds the shadow at its parameter t, seen by itself: `aside`
  // is how far to each side the rays are cast.
  bool bounds(size_t index, double t, double aside)
  {
    const Candidate& c = candidates[index];
    gp_Pnt           p;
    gp_Vec           d1;
    c.curve->D1(t, p, d1);
    const Image  way   = view.of(d1);
    const double speed = std::hypot(way.x, way.y);
    if (speed < 1e-9 * d1.Magnitude() || speed < 1e-300)
      return false;
    if (!onFace(c, p))
      return false;
    // An edge is outline only where its faces lie to one side of it in the image (or are
    // seen edge-on): where one lies to each side, the body goes on across it.
    bool toLeft = false, toRight = false;
    for (const Beside& face : c.beside)
    {
      const gp_Pnt2d uv = face.trace->Value(t);
      gp_Pnt         on;
      gp_Vec         du, dv;
      face.surface->D1(uv.X(), uv.Y(), on, du, dv);
      const gp_Vec inward = du.Crossed(dv).Crossed(d1) * face.left;
      const double size   = inward.Magnitude();
      if (size < 1e-300)
        continue;
      const Image  seenInward = view.of(inward);
      const double across     = (seenInward.x * -way.y + seenInward.y * way.x) / (speed * size);
      (across > edgeOn ? toLeft : toRight) |= std::abs(across) > edgeOn;
    }
    if (toLeft && toRight)
      return false;
    const gp_Vec side = (view.x * (-way.y / speed) + view.y * (way.x / speed)) * aside;
    return occupied(p.Translated(side)) != occupied(p.Translated(-side));
  }

  // The state of candidate `index` at its parameter t.
  bool keptAt(size_t index, double t) const
  {
    const Candidate& c    = candidates[index];
    const auto       next = std::upper_bound(c.events.begin(), c.events.end(), t);
    const long k = std::clamp(long(next - c.events.begin()) - 1, 0L, long(c.state.size()) - 1);
    return !c.state.empty() && c.state[size_t(k)] != 0;
  }

  // Whether another kept stretch with the same image, nearer the eye (or, at the same
  // depth, earlier in the list), takes the place of candidate `index` at t.
  bool covered(size_t index, double t) const
  {
    const Candidate& c = candidates[index];
    gp_Pnt           p;
    gp_Vec           d1;
    c.curve->D1(t, p, d1);
    const Image  q = view.of(p), way = view.of(d1);
    const double speed = std::hypot(way.x, way.y);
    const double depth = gp_Vec(view.origin, p).Dot(view.n);
    for (size_t other = 0; other < candidates.size(); ++other)
    {
      const Candidate& b = candidates[other];
      for (const double s : nearest(b, view, q, len.sameImage))
      {
        gp_Pnt bp;
        gp_Vec bd;
        b.curve->D1(s, bp, bd);
        if (other == index && bp.Distance(p) <= 10 * len.sameImage)
          continue;
        const Image  bway   = view.of(bd);
        const double bspeed = std::hypot(bway.x, bway.y);
        if (bspeed < 1e-300 || std::abs(way.x * bway.y - way.y * bway.x) > sameWay * speed * bspeed)
          continue;
        // A tangency is one point, not an overlap. Near it two distinct curves
        // can be inside sameImage; only curves that continue together can cover.
        bool together = true;
        for (const double side : {-1.0, 1.0})
        {
          const double test    = std::clamp(t + side * 0.02 * (c.b - c.a), c.a, c.b);
          const auto   matched = nearest(b, view, imageOf(c, view, test), len.sameImage);
          if (matched.empty())
            together = false;
        }
        if (!together || !keptAt(other, s))
          continue;
        const double bdepth = gp_Vec(view.origin, bp).Dot(view.n);
        if (bdepth > depth + len.sameDepth)
          return true;
        if (bdepth >= depth - len.sameDepth && (other < index || (other == index && s < t)))
          return true;
      }
    }
    return false;
  }

  // Where the images of a and b cross, by Newton from (u, s), within the given stretches.
  bool crossing(const Candidate& a,
                const Candidate& b,
                double&          u,
                double&          s,
                double           u0,
                double           u1,
                double           s0,
                double           s1) const
  {
    for (int k = 0; k < 40; ++k)
    {
      gp_Pnt p, bp;
      gp_Vec d, bd;
      a.curve->D1(u, p, d);
      b.curve->D1(s, bp, bd);
      const Image  fa = view.of(p), fb = view.of(bp), ta = view.of(d), tb = view.of(bd);
      const double fx = fa.x - fb.x, fy = fa.y - fb.y;
      const double sa = std::hypot(ta.x, ta.y), sb = std::hypot(tb.x, tb.y);
      const double det = -ta.x * tb.y + tb.x * ta.y;
      if (std::abs(det) < 1e-6 * sa * sb)
        return false; // they run the same way: no crossing to find
      const double du = (tb.x * fy - fx * tb.y) / det, ds = (ta.x * fy - ta.y * fx) / det;
      const double un = std::clamp(u - du, u0, u1), sn = std::clamp(s - ds, s0, s1);
      const double moved = std::abs(un - u) * sa + std::abs(sn - s) * sb;
      u                  = un;
      s                  = sn;
      if (moved <= 1e-14 * len.span)
        return gap(imageOf(a, view, u), imageOf(b, view, s)) <= len.point;
    }
    return false;
  }

  // What the pieces ka of candidate i and kb of candidate j do to each other.
  void meet(size_t i, size_t ka, size_t j, size_t kb)
  {
    const Candidate &a = candidates[i], &b = candidates[j];
    // The nearest points of the two pieces of the polylines.
    const Image a0 = a.seen[ka], a1 = a.seen[ka + 1], b0 = b.seen[kb], b1 = b.seen[kb + 1];
    double      la = 0, lb = 0, least;
    {
      const double rx = a1.x - a0.x, ry = a1.y - a0.y, sx = b1.x - b0.x, sy = b1.y - b0.y;
      const double det   = rx * sy - ry * sx;
      bool         cross = false;
      if (std::abs(det) > 1e-300)
      {
        la    = ((b0.x - a0.x) * sy - (b0.y - a0.y) * sx) / det;
        lb    = ((b0.x - a0.x) * ry - (b0.y - a0.y) * rx) / det;
        cross = la >= 0 && la <= 1 && lb >= 0 && lb <= 1;
      }
      if (cross)
      {
        least = 0;
      }
      else
      {
        double share;
        least = fromPiece(a0, b0, b1, share);
        la    = 0;
        lb    = share;
        if (const double d = fromPiece(a1, b0, b1, share); d < least)
        {
          least = d;
          la    = 1;
          lb    = share;
        }
        if (const double d = fromPiece(b0, a0, a1, share); d < least)
        {
          least = d;
          la    = share;
          lb    = 0;
        }
        if (const double d = fromPiece(b1, a0, a1, share); d < least)
        {
          least = d;
          la    = share;
          lb    = 1;
        }
      }
    }
    if (least > a.slack[ka] + b.slack[kb] + len.sameImage)
      return;
    const double t0 = a.at[ka], t1 = a.at[ka + 1], s0 = b.at[kb], s1 = b.at[kb + 1];
    const double ustart = t0 + (t1 - t0) * la, sstart = s0 + (s1 - s0) * lb;
    const double half = len.sameImage / 2;
    // How far the image of a at t is from the image of b; `beside` follows the nearest
    // parameter of b.
    // (Of a candidate against itself: of the other piece and what follows it, not of the
    // place itself.)
    const double low  = i == j ? (t1 + s0) / 2 : b.a,
                 high = i == j ? std::min(b.b, 2 * s1 - s0) : b.b;
    double beside     = sstart;
    auto   off = [&](double t) { return apart(a, t, b, beside, view, len.point, low, high); };
    // The same image all along this piece: nothing happens here.
    const bool in0 = off(t0) <= half, in1 = off(t1) <= half;
    if (in0 && in1)
      return;
    // The crossing may lie a little past the ends of the pieces of the polylines (but a
    // candidate does not cross itself where it merely runs on).
    const double more = 0.5;
    double       u = in0 ? t0 : in1 ? t1 : ustart, s = in0 || in1 ? beside : sstart;
    if (crossing(a,
                 b,
                 u,
                 s,
                 std::max(a.a, t0 - more * (t1 - t0)),
                 std::min(a.b, t1 + more * (t1 - t0)),
                 std::max(b.a, s0 - more * (s1 - s0)),
                 std::min(b.b, s1 + more * (s1 - s0)))
        && (i != j || s - u > (s0 - t1) / 4))
    {
      event(i, u);
      event(j, s);
      return;
    }
    if (in0 || in1)
    {
      contact(i, j, in0 ? t0 : t1, sstart, t1 - t0);
      return;
    }
    // They may touch between the ends of the piece.
    u = ustart;
    s = sstart;
    for (int k = 0; k < 60; ++k)
    {
      const double un = std::clamp(foot(a, view, imageOf(b, view, s), u), t0, t1);
      const double sn = std::clamp(foot(b, view, imageOf(a, view, un), s), s0, s1);
      const bool   settled =
        std::abs(un - u) + std::abs(sn - s) <= 1e-15 * (std::abs(u) + std::abs(s) + 1);
      u = un;
      s = sn;
      if (settled)
        break;
    }
    if (gap(imageOf(a, view, u), imageOf(b, view, s)) <= half)
      contact(i, j, u, s, t1 - t0);
  }

  // The images of candidates i and j are the same at the parameter `in` of i (near `near` of
  // j). Either they touch there, and the event is the place where they touch; or they run
  // together for a stretch, and the events are where they part. `step` is how far along i to
  // look for the ends of the stretch at a time.
  void contact(size_t i, size_t j, double in, double near, double step)
  {
    const Candidate &a = candidates[i], &b = candidates[j];
    const double     half = len.sameImage / 2;
    // At a trimmed branch's end, the endpoint itself is the event. A tolerance
    // band around its tangency would stop the neighbouring curve too early.
    if (!b.closed)
    {
      for (const double end : {b.a, b.b})
      {
        const Image  q = imageOf(b, view, end);
        const double t = foot(a, view, q, in);
        if (std::abs(t - in) <= 2 * step && gap(imageOf(a, view, t), q) <= len.sameImage)
        {
          event(i, t);
          event(j, end);
          return;
        }
      }
    }
    double beside = near;
    auto   off    = [&](double t) { return apart(a, t, b, beside, view, len.point, b.a, b.b); };
    // Where, to one side of `in`, the images part: found within two steps, or not at all
    // (the candidate ends first, or they run on together).
    auto parting = [&](double way, double& at) {
      double inside = in;
      beside        = near;
      for (int k = 0; k < 2; ++k)
      {
        const double out = std::clamp(inside + way * step, a.a, a.b);
        if (out == inside)
          return false;
        const double keep = beside;
        if (off(out) <= half)
        {
          inside = out;
          continue;
        }
        beside    = keep;
        double lo = inside, hi = out;
        for (int n = 0; n < 60; ++n)
        {
          const double mid    = (lo + hi) / 2;
          const double before = beside;
          if (off(mid) <= half)
          {
            lo = mid;
          }
          else
          {
            hi     = mid;
            beside = before;
          }
        }
        at = (lo + hi) / 2;
        return true;
      }
      return false;
    };
    double     low = in, high = in;
    const bool below = parting(-1, low), above = parting(1, high);
    if (below && above)
    {
      // A short stretch. Where two curves touch, they are apart by a quarter of what
      // they are at its ends when half way to them; where they are one curve for the
      // stretch, by nothing.
      beside              = near;
      const double middle = (low + high) / 2;
      const bool   touching =
        off((low + middle) / 2) > 0.02 * half || off((middle + high) / 2) > 0.02 * half;
      if (touching)
      {
        // The place where they run the same way.
        auto turned = [&](double t) {
          gp_Pnt p, bp;
          gp_Vec d, bd;
          a.curve->D1(t, p, d);
          beside = foot(b, view, view.of(p), beside);
          b.curve->D1(beside, bp, bd);
          const Image ta = view.of(d), tb = view.of(bd);
          return ta.x * tb.y - ta.y * tb.x;
        };
        beside             = near;
        const double first = turned(low), last = turned(high);
        double       at = middle;
        if ((first >= 0) != (last >= 0))
        {
          double lo = low, hi = high;
          for (int n = 0; n < 60; ++n)
          {
            const double mid                               = (lo + hi) / 2;
            ((turned(mid) >= 0) == (first >= 0) ? lo : hi) = mid;
          }
          at = (lo + hi) / 2;
        }
        else
        {
          double lo = low, hi = high;
          for (int n = 0; n < 80; ++n)
          {
            const double m1 = lo + (hi - lo) / 3, m2 = hi - (hi - lo) / 3;
            (off(m1) < off(m2) ? hi : lo) = off(m1) < off(m2) ? m2 : m1;
          }
          at = (lo + hi) / 2;
        }
        event(i, at);
        return;
      }
    }
    if (below)
      event(i, low);
    if (above)
      event(i, high);
  }
};

// The sine of the angle between a surface's normal and the plane square to the view, over a
// rectangle of its parameters.
struct Facing
{
  Handle(Geom_Surface) surface;
  gp_Vec               view;
  double               u1, u2, v1, v2;

  double raw(double u, double v) const
  {
    gp_Pnt p;
    gp_Vec du, dv;
    surface->D1(u, v, p, du, dv);
    const gp_Vec n    = du.Crossed(dv);
    const double size = n.Magnitude();
    return size <= 1e-10 * (du.SquareMagnitude() + dv.SquareMagnitude()) ? std::nan("")
                                                                         : n.Dot(view) / size;
  }

  // At a pole of the surface the normal is taken a little way in: NaN only where that has
  // none either.
  double operator()(double u, double v) const
  {
    const double value = raw(u, v);
    if (!std::isnan(value))
      return value;
    return raw(u + ((u1 + u2) / 2 - u) * 1e-6, v + ((v1 + v2) / 2 - v) * 1e-6);
  }
};

struct OnSurface
{
  gp_Pnt p;
  double u, v;
};

// A curve through points of a silhouette, made finer until, between each two of them, it is
// within `len.traced` of the surface and of the place where the surface turns away from the eye.
// `between` gives a point of the silhouette between two given ones.
Handle(Geom_Curve) fitted(
  const char*                                                                               op,
  std::vector<OnSurface>                                                                    points,
  bool                                                                                      closed,
  const Facing&                                                                             facing,
  const Handle(ShapeAnalysis_Surface)&                                                      project,
  const Lengths&                                                                            len,
  const std::function<bool(const OnSurface&, const OnSurface&, const gp_Pnt&, OnSurface&)>& between)
{
  const auto&  surface = facing.surface;
  const double hu = 1e-6 * (facing.u2 - facing.u1), hv = 1e-6 * (facing.v2 - facing.v1);
  // How far a point of the surface is from the silhouette, along the surface.
  auto toSilhouette = [&](double u, double v) {
    const double g  = facing(u, v);
    const double gu = (facing(u + hu, v) - facing(u - hu, v)) / (2 * hu);
    const double gv = (facing(u, v + hv) - facing(u, v - hv)) / (2 * hv);
    gp_Pnt       p;
    gp_Vec       du, dv;
    surface->D1(u, v, p, du, dv);
    const double e = du.Dot(du), f = du.Dot(dv), k = dv.Dot(dv), det = e * k - f * f;
    if (std::isnan(g) || std::isnan(gu) || std::isnan(gv) || det <= 0)
      return std::numeric_limits<double>::infinity();
    const double steep = (k * gu * gu - 2 * f * gu * gv + e * gv * gv) / det;
    return steep > 0 ? std::abs(g) / std::sqrt(steep) : std::numeric_limits<double>::infinity();
  };
  for (int round = 0; round < 32; ++round)
  {
    std::vector<OnSurface> distinct;
    for (const auto& p : points)
      if (distinct.empty()
          || distinct.back().p.Distance(p.p) > BRepExt::tolerance::exactLength * len.span)
        distinct.push_back(p);
    if (closed && distinct.size() > 1
        && distinct.front().p.Distance(distinct.back().p)
             <= BRepExt::tolerance::exactLength * len.span)
      distinct.pop_back();
    if (distinct.size() < (closed ? 3u : 2u))
      return nullptr;
    points                                    = distinct;
    const int                           count = int(points.size());
    Handle(NCollection_HArray1<gp_Pnt>) array = new NCollection_HArray1<gp_Pnt>(1, count);
    Handle(NCollection_HArray1<double>) knots =
      new NCollection_HArray1<double>(1, count + (closed ? 1 : 0));
    double along = 0;
    for (int i = 0; i < count; ++i)
    {
      array->SetValue(i + 1, points[size_t(i)].p);
      knots->SetValue(i + 1, along);
      along += points[size_t(i)].p.Distance(points[size_t((i + 1) % count)].p);
    }
    if (closed)
      knots->SetValue(count + 1, along);
    GeomAPI_Interpolate interpolate(array, knots, closed, 1e-12 * len.span);
    if (!closed)
    {
      auto tangent = [&](const OnSurface& at, const gp_Vec& towards) {
        const double gu = (facing(at.u + hu, at.v) - facing(at.u - hu, at.v)) / (2 * hu);
        const double gv = (facing(at.u, at.v + hv) - facing(at.u, at.v - hv)) / (2 * hv);
        gp_Pnt       p;
        gp_Vec       du, dv;
        surface->D1(at.u, at.v, p, du, dv);
        gp_Vec way = du * gv - dv * gu;
        if (way.Magnitude() <= 1e-300 || !std::isfinite(way.SquareMagnitude()))
          return towards;
        if (way.Dot(towards) < 0)
          way.Reverse();
        return way.Normalized();
      };
      interpolate.Load(
        tangent(points.front(), gp_Vec(points[0].p, points[1].p)),
        tangent(points.back(), gp_Vec(points[size_t(count - 2)].p, points.back().p)));
    }
    interpolate.Perform();
    if (!interpolate.IsDone())
      fail(op, no_solution, "the outline could not be drawn through its points on a face");
    Handle(Geom_Curve)     curve(interpolate.Curve());
    std::vector<OnSurface> finer;
    bool                   added = false;
    const int              spans = closed ? count : count - 1;
    for (int i = 0; i < spans; ++i)
    {
      const OnSurface &from = points[size_t(i)], &to = points[size_t((i + 1) % count)];
      finer.push_back(from);
      const double t    = (knots->Value(i + 1) + (i + 1 < count ? knots->Value(i + 2) : along)) / 2;
      const gp_Pnt p    = curve->Value(t);
      double       away = 0;
      const double first = knots->Value(i + 1), last = i + 1 < count ? knots->Value(i + 2) : along;
      for (const double share : {0.25, 0.5, 0.75})
      {
        const gp_Pnt2d uv =
          project->ValueOfUV(curve->Value(first + (last - first) * share), len.traced / 10);
        away = std::max(away, std::hypot(project->Gap(), toSilhouette(uv.X(), uv.Y())));
      }
      if (away <= len.traced)
        continue;
      OnSurface middle;
      if (!between(from, to, p, middle)
          || middle.p.Distance(from.p) <= BRepExt::tolerance::exactLength * len.span
          || middle.p.Distance(to.p) <= BRepExt::tolerance::exactLength * len.span)
      {
        fail(op, no_solution, "the adaptive surface silhouette could not be resolved to tolerance");
      }
      finer.push_back(middle);
      added = true;
    }
    if (!closed)
      finer.push_back(points.back());
    if (!added)
      return curve;
    if (finer.size() > 100000)
      break;
    points = finer;
  }
  fail(op, no_solution, "the adaptive surface silhouette could not be resolved to tolerance");
}

// The curves of a general surface where it turns away from the view. A grid over its
// parameters, as fine as its pieces and doubled while a cell hides a change of sign, gives how
// the curves run: each crossing of a grid line is found exactly and the crossings of a cell are
// joined. The curves drawn through them are then made finer until they are on the silhouette
// to tolerance.
void trace(const char*                          op,
           const Facing&                        facing,
           const Handle(ShapeAnalysis_Surface)& project,
           const Lengths&                       len,
           std::vector<Handle(Geom_Curve)>&     out,
           const std::vector<double>&           isoU = {},
           const std::vector<double>&           isoV = {})
{
  const auto&  surface = facing.surface;
  const double u1 = facing.u1, u2 = facing.u2, v1 = facing.v1, v2 = facing.v2;
  // Exact isoparametric branches are emitted separately. Divide them out of the
  // seed field so a crossing of two such branches is not rounded into a false turn.
  auto level = [&](double u, double v) {
    double factor = 1;
    for (const double root : isoU)
    {
      if (std::abs(u - root) < 1e-7 * (u2 - u1))
        u = root + 1e-7 * (u2 - u1);
      factor *= (u - root) / (u2 - u1);
    }
    for (const double root : isoV)
    {
      if (std::abs(v - root) < 1e-7 * (v2 - v1))
        v = root + 1e-7 * (v2 - v1);
      factor *= (v - root) / (v2 - v1);
    }
    return facing(u, v) / factor;
  };
  const GeomAdaptor_Surface pieces(surface, u1, u2, v1, v2);
  int                       nu = std::clamp(6 * pieces.NbUIntervals(GeomAbs_C2), 24, 96);
  int                       nv = std::clamp(6 * pieces.NbVIntervals(GeomAbs_C2), 24, 96);
  std::vector<double>       g;
  for (int attempt = 0;; ++attempt)
  {
    g.assign(size_t(nu + 1) * size_t(nv + 1), 0.0);
    auto   at  = [&](int i, int j) -> double& { return g[size_t(j) * size_t(nu + 1) + size_t(i)]; };
    auto   uOf = [&](int i) { return u1 + (u2 - u1) * i / nu; };
    auto   vOf = [&](int j) { return v1 + (v2 - v1) * j / nv; };
    double largest = 0;
    bool   plus = false, minus = false;
    for (int j = 0; j <= nv; ++j)
      for (int i = 0; i <= nu; ++i)
      {
        const double value = level(uOf(i), vOf(j));
        at(i, j)           = value;
        if (std::isnan(value))
          continue;
        largest                     = std::max(largest, std::abs(value));
        (value >= 0 ? plus : minus) = true;
      }
    // Seen edge-on all over (a wall along the view), the face has no curve of its own: its
    // edges are the outline.
    if (largest < 1e-9)
      return;
    // A cell whose corners agree, but where the surface is nearer to turning away than it
    // changes across the cell, is looked at in its middle and at the middles of its sides.
    bool hidden = false;
    for (int j = 0; j < nv && !hidden; ++j)
      for (int i = 0; i < nu && !hidden; ++i)
      {
        const double c[4] = {at(i, j), at(i + 1, j), at(i, j + 1), at(i + 1, j + 1)};
        if (std::isnan(c[0]) || std::isnan(c[1]) || std::isnan(c[2]) || std::isnan(c[3]))
          continue;
        const bool up = c[0] >= 0;
        if ((c[1] >= 0) != up || (c[2] >= 0) != up || (c[3] >= 0) != up)
          continue;
        const double low  = std::min({c[0], c[1], c[2], c[3]}),
                     high = std::max({c[0], c[1], c[2], c[3]});
        if (std::min(std::abs(low), std::abs(high)) > high - low)
          continue;
        const double um = (uOf(i) + uOf(i + 1)) / 2, vm = (vOf(j) + vOf(j + 1)) / 2;
        for (const auto& [u, v] : {std::pair{um, vm},
                                   std::pair{um, vOf(j)},
                                   std::pair{um, vOf(j + 1)},
                                   std::pair{uOf(i), vm},
                                   std::pair{uOf(i + 1), vm}})
        {
          const double value = level(u, v);
          hidden |= !std::isnan(value) && (value >= 0) != up;
        }
      }
    if (hidden && attempt < 3)
    {
      nu *= 2;
      nv *= 2;
      continue;
    }
    if (hidden)
      fail(op,
           no_solution,
           "the surface silhouette has unresolved branches after adaptive subdivision");
    if (!plus || !minus)
      return; // one sign all over: it does not turn away anywhere
    // Crossings on grid lines: 2 * (n + 1)^2 slots, even for those along u, odd along v.
    std::map<int, gp_Pnt2d> crossing;
    auto                    slot = [&](int i, int j, bool alongV) {
      return 2 * (j * (nu + 1) + i) + (alongV ? 1 : 0);
    };
    auto crosses = [&](double a, double b) {
      return !std::isnan(a) && !std::isnan(b) && (a >= 0) != (b >= 0);
    };
    for (int j = 0; j <= nv; ++j)
      for (int i = 0; i <= nu; ++i)
      {
        if (i < nu && crosses(at(i, j), at(i + 1, j)))
        {
          const double v = vOf(j);
          crossing[slot(i, j, false)] =
            gp_Pnt2d(root([&](double u) { return level(u, v); }, uOf(i), uOf(i + 1)), v);
        }
        if (j < nv && crosses(at(i, j), at(i, j + 1)))
        {
          const double u = uOf(i);
          crossing[slot(i, j, true)] =
            gp_Pnt2d(u, root([&](double v) { return level(u, v); }, vOf(j), vOf(j + 1)));
        }
      }
    std::map<int, std::vector<int>> next;
    auto                            join = [&](int a, int b) {
      next[a].push_back(b);
      next[b].push_back(a);
    };
    for (int j = 0; j < nv; ++j)
      for (int i = 0; i < nu; ++i)
      {
        const int bottom = slot(i, j, false), top = slot(i, j + 1, false), left = slot(i, j, true),
                  right = slot(i + 1, j, true);
        std::vector<int> around;
        for (const int s : {bottom, right, top, left})
          if (crossing.count(s))
            around.push_back(s);
        if (around.size() == 2)
        {
          join(around[0], around[1]);
        }
        else if (around.size() == 4)
        {
          // A saddle: the middle of the cell says which corners belong together.
          const double middle = level((uOf(i) + uOf(i + 1)) / 2, (vOf(j) + vOf(j + 1)) / 2);
          if ((middle >= 0) == (at(i, j) >= 0))
          {
            join(bottom, right);
            join(top, left);
          }
          else
          {
            join(bottom, left);
            join(top, right);
          }
        }
      }
    if (surface->IsUPeriodic())
      for (int j = 0; j < nv; ++j)
      {
        const int a = slot(0, j, true), b = slot(nu, j, true);
        if (crossing.count(a) && crossing.count(b) && next[a].size() == 1 && next[b].size() == 1)
          join(a, b);
      }
    if (surface->IsVPeriodic())
      for (int i = 0; i < nu; ++i)
      {
        const int a = slot(i, 0, false), b = slot(i, nv, false);
        if (crossing.count(a) && crossing.count(b) && next[a].size() == 1 && next[b].size() == 1)
          join(a, b);
      }
    // A point of the silhouette between two others: from the middle of their parameters,
    // down the slope of the facing to where it is nought.
    const double hu = 1e-6 * (u2 - u1), hv = 1e-6 * (v2 - v1);
    auto between = [&](const OnSurface& a, const OnSurface& b, const gp_Pnt&, OnSurface& found) {
      double bu = b.u, bv = b.v;
      if (surface->IsUPeriodic())
        bu = a.u + std::remainder(b.u - a.u, surface->UPeriod());
      if (surface->IsVPeriodic())
        bv = a.v + std::remainder(b.v - a.v, surface->VPeriod());
      double u = (a.u + bu) / 2, v = (a.v + bv) / 2;
      auto   wrap = [](double t, double low, double high) {
        return low + std::fmod(std::fmod(t - low, high - low) + high - low, high - low);
      };
      if (surface->IsUPeriodic())
        u = wrap(u, u1, u2);
      if (surface->IsVPeriodic())
        v = wrap(v, v1, v2);
      for (int k = 0; k < 20; ++k)
      {
        const double value = level(u, v);
        const double gu    = (level(u + hu, v) - level(u - hu, v)) / (2 * hu);
        const double gv    = (level(u, v + hv) - level(u, v - hv)) / (2 * hv);
        gp_Pnt       p;
        gp_Vec       du, dv;
        surface->D1(u, v, p, du, dv);
        const double e = du.Dot(du), f = du.Dot(dv), m = dv.Dot(dv), det = e * m - f * f;
        if (std::isnan(value) || std::isnan(gu) || std::isnan(gv) || det <= 0)
          return false;
        // The shortest step on the surface that takes the facing to nought.
        const double su = (m * gu - f * gv) / det, sv = (e * gv - f * gu) / det;
        const double steep = gu * su + gv * sv;
        if (!(steep > 0))
          return false;
        const double stepU = -value * su / steep, stepV = -value * sv / steep;
        u = std::clamp(u + stepU, u1, u2);
        v = std::clamp(v + stepV, v1, v2);
        if (std::abs(stepU) <= 1e-14 * (u2 - u1) && std::abs(stepV) <= 1e-14 * (v2 - v1))
          break;
      }
      found = {surface->Value(u, v), u, v};
      return true;
    };
    std::map<int, bool> used;
    auto                walk = [&](int start, bool closed) {
      std::vector<OnSurface> points;
      int                    current = start, previous = -1;
      while (true)
      {
        used[current]     = true;
        const gp_Pnt2d uv = crossing[current];
        points.push_back({surface->Value(uv.X(), uv.Y()), uv.X(), uv.Y()});
        int forward = -1;
        for (const int candidate : next[current])
          if (candidate != previous && !used[candidate])
            forward = candidate;
        if (forward < 0)
          break;
        previous = current;
        current  = forward;
      }
      const auto curve = fitted(op, points, closed, facing, project, len, between);
      if (!curve.IsNull())
        out.push_back(curve);
    };
    for (const auto& [s, list] : next)
      if (list.size() == 1 && !used[s])
        walk(s, false);
    for (const auto& [s, list] : next)
      if (!used[s])
        walk(s, true);
    return;
  }
}

// The silhouettes of one face: the curves on its surface where n . v = 0.
void silhouettes(const char*                         op,
                 const TopoDS_Face&                  face,
                 const OnFace&                       onFace,
                 const gp_Vec&                       view,
                 const Lengths&                      len,
                 std::vector<Handle(Geom_Curve)>&    out,
                 std::vector<std::array<double, 2>>& ranges)
{
  const BRepAdaptor_Surface adaptor(face);
  if (adaptor.GetType() == GeomAbs_Plane)
    return;
  const Handle(Geom_Surface) surface = BRep_Tool::Surface(face);
  double                     u1, u2, v1, v2;
  BRepTools::UVBounds(face, u1, u2, v1, v2);
  const Facing facing{surface, view, u1, u2, v1, v2};
  auto         add = [&](const Handle(Geom_Curve)& curve, double a, double b) {
    out.push_back(curve);
    ranges.push_back({a, b});
  };
  // The roots of the facing along u at one v, as lines of constant u.
  auto isoU = [&](double v) {
    constexpr int       steps = 16;
    std::vector<double> found;
    auto                f = [&](double u) { return facing.raw(u, v); };
    for (int i = 0; i < steps; ++i)
    {
      const double a = u1 + (u2 - u1) * i / steps, b = u1 + (u2 - u1) * (i + 1) / steps;
      roots(f, a, f(a), b, f(b), 0, found);
    }
    for (const double u : found)
      add(surface->UIso(u), v1, v2);
  };
  switch (adaptor.GetType())
  {
    case GeomAbs_Cylinder:
    case GeomAbs_Cone: {
      // The normal depends on u alone: the silhouettes are rulings. Along the axis of a
      // cylinder there are none.
      const double v       = (v1 + v2) / 2;
      double       largest = 0;
      for (int i = 0; i <= 16; ++i)
      {
        const double value = facing.raw(u1 + (u2 - u1) * i / 16, v);
        if (!std::isnan(value))
          largest = std::max(largest, std::abs(value));
      }
      if (largest > 1e-9)
        isoU(v);
      return;
    }
    case GeomAbs_SurfaceOfRevolution: {
      const auto revolution = Handle(Geom_SurfaceOfRevolution)::DownCast(surface);
      if (!revolution.IsNull()
          && view.Crossed(gp_Vec(revolution->Axis().Direction())).Magnitude()
               < BRepExt::tolerance::angle)
      {
        auto                f = [&](double v) { return facing.raw((u1 + u2) / 2, v); };
        std::vector<double> found;
        const int           count = std::max(16, 4 * adaptor.NbVIntervals(GeomAbs_C2));
        for (int i = 0; i < count; ++i)
        {
          const double a = v1 + (v2 - v1) * i / count, b = v1 + (v2 - v1) * (i + 1) / count;
          roots(f, a, f(a), b, f(b), 0, found);
        }
        for (const double v : found)
          add(surface->VIso(v), u1, u2);
        return;
      }
      std::vector<Handle(Geom_Curve)> traced;
      trace(op, facing, onFace.surface, len, traced);
      for (const auto& curve : traced)
        add(curve, curve->FirstParameter(), curve->LastParameter());
      return;
    }
    case GeomAbs_Sphere: {
      const gp_Sphere sphere = adaptor.Sphere();
      add(new Geom_Circle(gp_Ax2(sphere.Location(), gp_Dir(view)), sphere.Radius()), 0, 2 * M_PI);
      return;
    }
    case GeomAbs_Torus: {
      const gp_Torus torus = adaptor.Torus();
      const gp_Vec   axis(torus.Axis().Direction());
      const double   along = view.Dot(axis), across = (view - axis * along).Magnitude();
      if (across < BRepExt::tolerance::fineAngle)
      { // from above: the outer and the inner circle
        add(surface->VIso(0), u1, u2);
        add(surface->VIso(M_PI), u1, u2);
        return;
      }
      if (std::abs(along) < BRepExt::tolerance::fineAngle)
      { // side-on: the top and bottom circles and two of the tube
        add(surface->VIso(M_PI / 2), u1, u2);
        add(surface->VIso(3 * M_PI / 2), u1, u2);
        isoU(0);
        return;
      }
      if (std::abs(along) < 2e-4)
        fail(op,
             unsupported_geometry,
             "the torus is too nearly side-on to resolve its silhouette branches to tolerance");
      // Any other view. Along a silhouette the surface's normal runs once round the circle
      // of directions square to the view; the tube has each such normal n twice, on its
      // outer half (where the way out from the axis is that of n) and on its inner half
      // (where it is the opposite). So each of the two curves is given by the angle of its
      // normal, without a place where it could not be followed.
      const gp_Vec p1 = (axis - view * along).Normalized(), p2 = view.Crossed(p1);
      const gp_Pnt centre = torus.Location();
      const double major = torus.MajorRadius(), minor = torus.MinorRadius();
      for (const double sheet : {1.0, -1.0})
      {
        auto pointAt = [&](double angle) {
          const gp_Vec n     = p1 * std::cos(angle) + p2 * std::sin(angle);
          const gp_Vec level = n - axis * n.Dot(axis);
          const double size  = level.Magnitude();
          const gp_Vec out   = level * (sheet / size);
          const gp_Pnt p     = centre.Translated(out * major + n * minor);
          return OnSurface{p, angle, 0};
        };
        // Seen nearly side-on the curve swings half round the torus within a short
        // stretch of the angle: points are added until none is far from the next.
        std::vector<OnSurface> points;
        constexpr int          count = 64;
        for (int i = 0; i < count; ++i)
          points.push_back(pointAt(2 * M_PI * i / count));
        const double step = 2 * M_PI * (major + minor) / count;
        for (size_t i = 0; i < points.size() && points.size() < 100000;)
        {
          const double  from = points[i].u, to = i + 1 < points.size() ? points[i + 1].u : 2 * M_PI;
          const gp_Pnt& next = points[(i + 1) % points.size()].p;
          if (points[i].p.Distance(next) > step && to - from > 1e-13)
            points.insert(points.begin() + long(i) + 1, pointAt((from + to) / 2));
          else
            ++i;
        }
        auto between =
          [&](const OnSurface& a, const OnSurface& b, const gp_Pnt&, OnSurface& found) {
            const double to = b.u > a.u ? b.u : b.u + 2 * M_PI;
            found           = pointAt((a.u + to) / 2);
            found.u         = std::fmod(found.u, 2 * M_PI);
            return to - a.u > 1e-13;
          };
        const auto curve = fitted(op, points, true, facing, onFace.surface, len, between);
        if (!curve.IsNull())
          add(curve, curve->FirstParameter(), curve->LastParameter());
      }
      return;
    }
    default: {
      std::vector<double> uBranches, vBranches;
      const int           nu = std::max(24, 6 * adaptor.NbUIntervals(GeomAbs_C2));
      const int           nv = std::max(24, 6 * adaptor.NbVIntervals(GeomAbs_C2));
      auto branches          = [&](bool alongU, int count, int across, std::vector<double>& out) {
        const double low = alongU ? u1 : v1, high = alongU ? u2 : v2;
        auto         value = [&](double t, double share) {
          return alongU ? facing(t, v1 + (v2 - v1) * share) : facing(u1 + (u2 - u1) * share, t);
        };
        auto                f = [&](double t) { return value(t, 0.371); };
        std::vector<double> found;
        for (int i = 0; i < count; ++i)
        {
          const double a = low + (high - low) * i / count, b = low + (high - low) * (i + 1) / count;
          roots(f, a, f(a), b, f(b), 0, found);
        }
        for (const double t : found)
        {
          bool whole = true;
          for (int j = 0; j <= across; ++j)
            whole &= std::abs(value(t, double(j) / across)) < 1e-10;
          if (whole && (out.empty() || std::abs(t - out.back()) > 1e-8 * (high - low)))
            out.push_back(t);
        }
      };
      branches(true, nu, nv, uBranches);
      branches(false, nv, nu, vBranches);
      for (const double u : uBranches)
        add(surface->UIso(u), v1, v2);
      for (const double v : vBranches)
        add(surface->VIso(v), u1, u2);
      std::vector<Handle(Geom_Curve)> traced;
      trace(op, facing, onFace.surface, len, traced, uBranches, vBranches);
      for (const auto& curve : traced)
        add(curve, curve->FirstParameter(), curve->LastParameter());
    }
  }
}

} // namespace

EvolutionData outline(OcctKernel& kernel,
                      uint32_t    shapeId,
                      Vec3        origin,
                      Vec3        normal,
                      Vec3        x,
                      bool        in3d)
{
  const char* op    = "outline";
  const auto& given = get(kernel, shapeId);
  if (!TopExp_Explorer(given, TopAbs_FACE).More())
    fail(op, invalid_input, "outline requires a solid or sheet");
  const gp_Vec n(normal.x, normal.y, normal.z), towards(x.x, x.y, x.z);
  if (!(n.Magnitude() > 1e-15) || !(towards.Magnitude() > 1e-15)
      || !(n.Crossed(towards).Magnitude() > 1e-15))
    fail(op, invalid_input, "direction must be finite and nonzero");
  const gp_Ax2 axes(gp_Pnt(origin.x, origin.y, origin.z), gp_Dir(n), gp_Dir(towards));
  // The body is looked at where it would stand with the middle of its box at the origin.
  Bnd_Box bounds;
  // (From the faces themselves: a mesh the body may carry would give another box.)
  BRepBndLib::Add(given, bounds, false);
  double x0, y0, z0, x1, y1, z1;
  bounds.Get(x0, y0, z0, x1, y1, z1);
  const double size = std::sqrt(bounds.SquareExtent());
  if (!(size > 0))
    fail(op, invalid_input, "the body has no size");
  // OCCT's intersection algorithms also use Confusion internally. Work at one size
  // as well as one place so those internal lengths have the same relative meaning.
  constexpr double workingSpan = 10000;
  gp_Trsf          home;
  home.SetScale(gp_Pnt(0, 0, 0), size / workingSpan);
  home.SetTranslationPart(gp_Vec((x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2));
  const TopoDS_Shape shape = BRepBuilderAPI_Transform(given, home.Inverted(), true).Shape();
  Scene              scene;
  scene.view         = {gp_Pnt(0, 0, 0),
                        gp_Vec(axes.XDirection()),
                        gp_Vec(axes.YDirection()),
                        gp_Vec(axes.Direction())};
  scene.len          = Lengths(workingSpan);
  const Lengths& len = scene.len;
  if (!(len.span > 0))
    fail(op, invalid_input, "the body has no size");
  auto finish = [&](const TopoDS_Shape& result) {
    const TopoDS_Shape back = BRepBuilderAPI_Transform(result, home, true).Shape();
    if (!BRepCheck_Analyzer(back).IsValid())
      fail(op, no_solution, "result curves are invalid");
    EvolutionData out;
    out.resultId = store(kernel, back);
    if (in3d)
      return out;
    // The curves in space were stored only to be laid flat.
    try
    {
      const EvolutionData flat = projectPlane(kernel, out.resultId, origin, normal);
      release(kernel, out.resultId);
      return flat;
    }
    catch (...)
    {
      release(kernel, out.resultId);
      throw;
    }
  };
  const TopoDS_Shape axial = axialPrism(shape, scene.view, len);
  if (!axial.IsNull())
    return finish(axial);
  scene.load(shape);

  std::vector<std::vector<double>> folds;
  auto                             consider = [&](const Handle(Geom_Curve)&  curve,
                      double                     a,
                      double                     b,
                      int                        face,
                      const std::vector<Beside>& beside) {
    Candidate candidate{curve, a, b, face, beside, false, {}, {}, {}, {}, {}, {}, {}};
    candidate.closed = curve->Value(a).Distance(curve->Value(b)) <= len.point;
    std::vector<double> turns;
    if (prepare(candidate, scene.view, len, turns))
    {
      scene.candidates.push_back(candidate);
      folds.push_back(turns);
    }
  };
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> edges, faces;
  TopExp::MapShapes(shape, TopAbs_EDGE, edges);
  TopExp::MapShapes(shape, TopAbs_FACE, faces);
  std::vector<std::vector<Beside>> beside(size_t(edges.Extent()) + 1);
  for (int i = 1; i <= faces.Extent(); ++i)
  {
    // In a face taken forwards, it lies to the left of each of its edges as the edge is
    // oriented there, seen against the surface's normal.
    const auto face    = TopoDS::Face(faces(i).Oriented(TopAbs_FORWARD));
    const auto surface = BRep_Tool::Surface(face);
    const bool flat    = BRepAdaptor_Surface(face).GetType() == GeomAbs_Plane;
    for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next())
    {
      const auto edge = TopoDS::Edge(it.Current());
      double     a, b;
      const auto trace = BRep_Tool::CurveOnSurface(edge, face, a, b);
      if (!trace.IsNull())
        beside[size_t(edges.FindIndex(edge))].push_back(
          {surface, trace, edge.Orientation() == TopAbs_REVERSED ? -1.0 : 1.0, flat});
    }
  }
  for (int i = 1; i <= edges.Extent(); ++i)
  {
    const auto edge = TopoDS::Edge(edges(i));
    bool       seam = false;
    for (int j = 1; j <= faces.Extent() && !seam; ++j)
      seam = BRep_Tool::IsClosed(edge, TopoDS::Face(faces(j)));
    if (seam || BRep_Tool::Degenerated(edge))
      continue;
    double     a, b;
    const auto curve = BRep_Tool::Curve(edge, a, b);
    if (!curve.IsNull())
      consider(curve, a, b, -1, beside[size_t(i)]);
  }
  for (int i = 1; i <= faces.Extent(); ++i)
  {
    const auto                         face = TopoDS::Face(faces(i));
    const OnFace                       onFace{new ShapeAnalysis_Surface(BRep_Tool::Surface(face)),
                        std::make_shared<BRepTopAdaptor_FClass2d>(face, Precision::Confusion())};
    std::vector<Handle(Geom_Curve)>    found;
    std::vector<std::array<double, 2>> ranges;
    silhouettes(op, face, onFace, scene.view.n, len, found, ranges);
    if (found.empty())
      continue;
    scene.faces.push_back(onFace);
    for (size_t k = 0; k < found.size(); ++k)
      consider(found[k], ranges[k][0], ranges[k][1], int(scene.faces.size()) - 1, {});
  }

  // The events. First those a candidate has by itself: its ends, its folds, and where a
  // curved face beside an edge turns edge-on.
  auto& all = scene.candidates;
  for (size_t index = 0; index < all.size(); ++index)
  {
    const Candidate& c = all[index];
    if (!c.closed)
    {
      scene.event(imageOf(c, scene.view, c.a));
      scene.event(imageOf(c, scene.view, c.b));
    }
    for (const double t : folds[index])
      scene.event(index, t);
    for (const Beside& face : c.beside)
    {
      if (face.flat)
        continue;
      const GeomAdaptor_Surface kind(face.surface);
      if (kind.GetType() == GeomAbs_Cylinder
          && std::abs(gp_Vec(kind.Cylinder().Axis().Direction()).Dot(scene.view.n)) >= 1 - 1e-14)
        continue; // along its axis the cylinder is edge-on everywhere
      auto turned = [&](double t) {
        const gp_Pnt2d uv = face.trace->Value(t);
        gp_Pnt         p;
        gp_Vec         du, dv;
        face.surface->D1(uv.X(), uv.Y(), p, du, dv);
        const gp_Vec normal = du.Crossed(dv);
        const double size   = normal.Magnitude();
        return size < 1e-300 ? std::nan("") : normal.Dot(scene.view.n) / size;
      };
      std::vector<double>       found;
      const std::vector<double> at = c.at;
      for (size_t k = 0; k + 1 < at.size(); ++k)
        roots(turned, at[k], turned(at[k]), at[k + 1], turned(at[k + 1]), 0, found);
      for (const double t : found)
      {
        scene.event(index, t);
      }
    }
  }
  // Then where two images meet.
  for (size_t i = 0; i < all.size(); ++i)
    for (size_t j = i; j < all.size(); ++j)
    {
      const double pad = len.sameImage;
      if (all[i].low[0] > all[j].high[0] + pad || all[j].low[0] > all[i].high[0] + pad
          || all[i].low[1] > all[j].high[1] + pad || all[j].low[1] > all[i].high[1] + pad)
        continue;
      const size_t na = all[i].at.size() - 1, nb = all[j].at.size() - 1;
      for (size_t ka = 0; ka < na; ++ka)
      {
        const double reach = all[i].slack[ka] + pad;
        const double lowX  = std::min(all[i].seen[ka].x, all[i].seen[ka + 1].x) - reach,
                     highX = std::max(all[i].seen[ka].x, all[i].seen[ka + 1].x) + reach,
                     lowY  = std::min(all[i].seen[ka].y, all[i].seen[ka + 1].y) - reach,
                     highY = std::max(all[i].seen[ka].y, all[i].seen[ka + 1].y) + reach;
        if (lowX > all[j].high[0] || highX < all[j].low[0] || lowY > all[j].high[1]
            || highY < all[j].low[1])
          continue;
        for (size_t kb = i == j ? ka + 2 : 0; kb < nb; ++kb)
        {
          if (i == j && all[i].closed && ka == 0 && kb + 1 == nb)
            continue;
          const double more = all[j].slack[kb];
          const Image  b0 = all[j].seen[kb], b1 = all[j].seen[kb + 1];
          if (std::min(b0.x, b1.x) - more > highX || std::max(b0.x, b1.x) + more < lowX
              || std::min(b0.y, b1.y) - more > highY || std::max(b0.y, b1.y) + more < lowY)
            continue;
          scene.meet(i, ka, j, kb);
        }
      }
    }
  // Every candidate through the place of an event is cut there.
  for (size_t index = 0; index < all.size(); ++index)
  {
    Candidate&          c      = all[index];
    std::vector<double> direct = c.events;
    direct.push_back(c.a);
    direct.push_back(c.b);
    const auto isDirect = [&](double t) {
      return std::find(direct.begin(), direct.end(), t) != direct.end();
    };
    for (const Image q : scene.points)
      for (const double t : nearest(c, scene.view, q, len.sameImage))
        c.events.push_back(t);
    c.events.push_back(c.a);
    c.events.push_back(c.b);
    std::sort(c.events.begin(), c.events.end());
    // Two events are one where the candidate does not leave the place between them (a
    // candidate that crosses itself is at one place twice, with a way between).
    auto same = [&](double from, double to) {
      if (to - from <= 1e-13 * (c.b - c.a))
        return true;
      const Image q = imageOf(c, scene.view, from);
      return gap(imageOf(c, scene.view, to), q) <= len.one
             && gap(imageOf(c, scene.view, (from + to) / 2), q) <= len.one;
    };
    std::vector<double> apart;
    for (const double t : c.events)
    {
      if (!apart.empty() && same(apart.back(), t))
      {
        // The end of the candidate itself is the one that stays.
        if ((t >= c.b && apart.size() > 1) || (isDirect(t) && !isDirect(apart.back())))
          apart.back() = t;
        continue;
      }
      apart.push_back(t);
    }
    if (apart.size() < 2)
      apart = {c.a, c.b};
    c.events = apart;
  }
  // Between two events a candidate is outline or not.
  for (size_t index = 0; index < all.size(); ++index)
  {
    Candidate& c = all[index];
    c.state.assign(c.events.size() - 1, 0);
    for (size_t k = 0; k + 1 < c.events.size(); ++k)
    {
      const double from = c.events[k], to = c.events[k + 1], middle = (from + to) / 2;
      const Image  q = imageOf(c, scene.view, middle);
      const double length =
        gap(imageOf(c, scene.view, from), q) + gap(q, imageOf(c, scene.view, to));
      double aside = std::min(len.aside, std::max(length / 4, 3 * len.rays));
      // Near a tangent joint the distance between projected branches is quadratic
      // in their length along the curve. A ray must stay inside that interval too.
      for (size_t other = 0; other < all.size(); ++other)
      {
        if (other == index)
          continue;
        for (const double t : nearest(all[other], scene.view, q, 2 * aside))
        {
          const double distance = gap(imageOf(all[other], scene.view, t), q);
          if (distance > len.point)
            aside = std::min(aside, distance / 4);
        }
      }
      c.state[k] = scene.bounds(index, middle, aside) ? 1 : 0;
    }
  }
  std::vector<std::vector<char>> stays(all.size());
  for (size_t index = 0; index < all.size(); ++index)
  {
    const Candidate& c = all[index];
    stays[index].assign(c.state.size(), 0);
    for (size_t k = 0; k < c.state.size(); ++k)
      stays[index][k] =
        c.state[k] && !scene.covered(index, (c.events[k] + c.events[k + 1]) / 2) ? 1 : 0;
  }

  BRep_Builder    builder;
  TopoDS_Compound result;
  builder.MakeCompound(result);
  std::vector<Image> ends;
  for (size_t index = 0; index < all.size(); ++index)
  {
    const Candidate& c    = all[index];
    auto             emit = [&](double a, double b) {
      const gp_Pnt first = c.curve->Value(a);
      if (first.Distance(c.curve->Value(a + (b - a) / 4)) <= len.point
          && first.Distance(c.curve->Value((a + b) / 2)) <= len.point
          && first.Distance(c.curve->Value(a + 3 * (b - a) / 4)) <= len.point)
        return;
      BRepBuilderAPI_MakeEdge edge(c.curve, a, b);
      if (!edge.IsDone())
        fail(op, no_solution, "a silhouette edge could not be built: no outline is returned");
      builder.Add(result, edge.Edge());
      const Image from = imageOf(c, scene.view, a), to = imageOf(c, scene.view, b);
      if (gap(from, to) > len.point || !c.closed || a > c.a || b < c.b)
      {
        ends.push_back(from);
        ends.push_back(to);
      }
    };
    for (size_t k = 0; k < stays[index].size();)
    {
      if (!stays[index][k])
      {
        ++k;
        continue;
      }
      size_t last = k;
      while (last + 1 < stays[index].size() && stays[index][last + 1])
        ++last;
      emit(c.events[k], c.events[last + 1]);
      k = last + 1;
    }
  }
  if (!TopExp_Explorer(result, TopAbs_EDGE).More())
    fail(op, no_solution, "no curves were produced");
  // The boundary of a shadow is closed: every end of a piece has to be the end of another
  // piece (the ends are kept in pairs, piece by piece).
  for (size_t i = 0; i < ends.size(); ++i)
  {
    bool met = false;
    for (size_t j = 0; j < ends.size() && !met; ++j)
      met = j != i && j != (i ^ 1) && gap(ends[i], ends[j]) <= 5 * len.sameImage;
    if (!met)
    {
      fail(op,
           no_solution,
           "the projected silhouette pieces do not meet to tolerance: no outline is returned");
    }
  }
  return finish(result);
}
} // namespace BRepExt
