// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
// The blend of the whole outline of a flat face against walls that stand square on it. See
// BRepFilletExt_BlendFace.cxx for the surface that is built.
//
// How:
//
//   1. The body is cut off (a convex outline) or filled up (a concave one, a pocket's floor)
//      at the blend's distance from the face, by the Boolean on cells. What that leaves in
//      place of the face, "the foot", has the outline where the blend leaves the walls.
//   2. In the foot's plane, for each place of the outline: how far along its inward normal it
//      is still the nearest place of the outline ("its reach"), and which other place stops it.
//      A reach of the whole distance ends on the face's own plane; a shorter one ends in a
//      crease, over the outline's medial axis. Where what stops a place changes, the crease
//      forks or meets the plane: those places are found by bisection.
//   3. Over each edge of the outline one face of the blend: a cylinder or a plane over a
//      straight edge, a spline surface over any other; bounded by the edge, the sections at
//      its two ends, and the creases and the pieces of the face's plane along its top. A
//      crease is one edge between the two faces that meet in it (the same face twice where an
//      edge's own bend closes it).
//   4. What is left of the face: flat faces inside the pieces of the outline moved in.
//   5. The foot is taken out of the cut body and these faces are put in its place. No Boolean
//      has part in that: every edge is shared as it was built.

#include <BRepAlgoExt_FaceBlend.hxx>
#include <BRepFilletExt_BlendChecks.hxx>

#include <BRepAlgoExt_Prism.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepOffsetExt_Solid.hxx>
#include <GeomAlgoExt_SurfaceJet.hxx>

#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepLib.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepOffsetAPI_MakeOffset.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <Geom2dAPI_Interpolate.hxx>
#include <Geom2dAdaptor_Curve.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_Line.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Circle.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Line.hxx>
#include <Geom_Plane.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <NCollection_HArray1.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Pln.hxx>
#include <gp_Lin.hxx>
#include <math_Gauss.hxx>
#include <math_Matrix.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <map>
#include <set>
#include <tuple>

#include <BRepAlgoExt_Properties.hxx>

namespace BRepExt
{

namespace
{

//! The normal of a plane as its parametrisation has it: d/du x d/dv.
gp_Vec planeNormal(const gp_Pln& plane)
{
  return gp_Vec(plane.XAxis().Direction()).Crossed(gp_Vec(plane.YAxis().Direction()));
}

//! The volume of a body as the kernel measures it: curved faces integrated to 1e-8.
double held(const TopoDS_Shape& shape, bool bracket = false)
{
  Bnd_Box box;
  BRepBndLib::Add(shape, box, false);
  if (box.IsVoid() || box.IsOpen())
    return 0;
  const gp_Pnt origin((box.CornerMin().XYZ() + box.CornerMax().XYZ()) / 2);
  double       total = 0;
  for (TopExp_Explorer it(shape, TopAbs_FACE); it.More(); it.Next())
    total += bracket ? faceMeasureBracket(TopoDS::Face(it.Current()), origin)
                     : faceMeasure(TopoDS::Face(it.Current()), false, origin, nullptr);
  // A mirrored body's faces are integrated with the mirror's sign.
  return std::abs(total);
}

// ---- Tolerances and sizes, by name ---------------------------------------------------------

//! The degree of the spline surfaces along the outline; each span is fitted by itself.
constexpr int kDegree = 9;
//! A span of a spline surface is halved until it is this near the outline and the outline moved in.
constexpr double kFit = 2e-10;
//! Places along each edge at which its reach is first looked at.
constexpr int kScan = 96;
//! A crease is refined until it is this near the true crease between its points.
constexpr double kCrease = 5e-8;
//! Two places of the outline nearer than this share of the spacing of its polyline are not
//! asked which stops which; and nothing is looked for within five times that of an end of an
//! edge, where the edge that goes on is as near.
constexpr double kNear = 1e-2;

// Walls stand square on the face within tolerance::angle (radians): the blend is tangent to
// them within it.

//! What stops a place of the outline: nothing before the blend's whole distance (`by` -1), or
//! the place `t` of the edge `by` (the edge itself at its own place: its own bend).
struct Hit
{
  double d  = 0;
  int    by = -1;
  double t  = 0;
};

struct Frame
{
  gp_Pnt2d p;
  gp_Vec2d n, dn;    // the inward normal, and how it turns
  double   bend = 0; // curvature towards the inside
  gp_Vec2d along;    // unit tangent
  double   curveSpeed = 0;
};

struct Element
{
  TopoDS_Edge           edge;  // as its curve runs
  double                s = 1; // +1: the face is on the left of the curve as its parameter grows
  Geom2dAdaptor_Curve   curve;
  double                a = 0, b = 0;
  GeomAbs_CurveType     kind = GeomAbs_OtherCurve;
  std::vector<double>   ts;
  std::vector<gp_Pnt2d> points;
  double                spacing = 0;
  // Junctions with the edges before (at a) and after (at b) in parameter order.
  int    before = -1, after = -1;        // the edge that shares the vertex at a / at b
  double beforeAt = 0, afterAt = 0;      // its parameter there
  bool   sharpA = false, sharpB = false; // a convex corner there: the reach is nothing
  // The blend's face over it.
  Handle(Geom_Surface) surface;
  bool straight = false;      // a cylinder or a plane, with (u, v) from the edge's parameter
  bool flipped  = false;      // which way the cylinder's angle runs
  Handle(Geom2d_Curve) moved; // the outline moved in by the distance, in the plane
  Handle(Geom_Curve)   movedInSpace;

  Frame frame(double t) const
  {
    Frame    f;
    gp_Vec2d d1, d2;
    curve.D2(t, f.p, d1, d2);
    const double length = d1.Magnitude();
    f.curveSpeed        = length;
    f.along             = d1 / length;
    f.n                 = gp_Vec2d(-f.along.Y() * s, f.along.X() * s);
    const gp_Vec2d turn = (d2 - f.along * f.along.Dot(d2)) / length;
    f.dn                = gp_Vec2d(-turn.Y() * s, turn.X() * s);
    f.bend              = d2.Dot(f.n) / (length * length);
    return f;
  }
};

//! The matrix that takes values at the fitting nodes of [0, 1] to Bezier poles of kDegree.
struct Fit
{
  std::array<double, kDegree + 1>                          nodes;
  std::array<std::array<double, kDegree + 1>, kDegree + 1> inverse;

  Fit()
  {
    math_Matrix m(0, kDegree, 0, kDegree), inverted(0, kDegree, 0, kDegree);
    for (int k = 0; k <= kDegree; ++k)
    {
      nodes[k]        = (1 - std::cos(M_PI * k / kDegree)) / 2;
      double binomial = 1;
      for (int j = 0; j <= kDegree; ++j)
      {
        m(k, j)  = binomial * std::pow(nodes[k], j) * std::pow(1 - nodes[k], kDegree - j);
        binomial = binomial * (kDegree - j) / (j + 1);
      }
    }
    math_Gauss(m).Invert(inverted);
    for (int k = 0; k <= kDegree; ++k)
      for (int j = 0; j <= kDegree; ++j)
        inverse[k][j] = inverted(k, j);
  }
};

double bezierAt(const std::array<double, kDegree + 1>& poles, double s)
{
  std::array<double, kDegree + 1> w = poles;
  for (int level = kDegree; level > 0; --level)
    for (int i = 0; i < level; ++i)
      w[i] = w[i] * (1 - s) + w[i + 1] * s;
  return w[0];
}

struct Blend
{
  const char* op;
  double      r;
  bool        bevel;
  gp_Pln      foot; // the plane of the foot; the face's own plane is `r` from it against `w`
  gp_Vec      w;    // from the face's plane along the walls
  gp_Vec      X, Y;
  std::vector<Element> elements;
  double               scale = 1;

  // ---- The section ----
  //! How far from the face's plane the blend is, over a place `d` from the outline.
  double height(double d) const
  {
    if (bevel)
      return r - d;
    const double q = r - d;
    return r - std::sqrt(std::max(0.0, r * r - q * q));
  }

  gp_Pnt lift(const gp_Pnt2d& p, double h) const
  {
    return foot.Location().Translated(X * p.X() + Y * p.Y() + w * (h - r));
  }

  gp_Pnt at(const Element& e, double t, double d) const
  {
    const Frame f = e.frame(t);
    return lift(f.p.Translated(f.n * d), height(d));
  }

  //! The parameter across a spline face at the distance d from its edge.
  double across(double d) const
  {
    if (bevel)
      return d / r;
    if (d <= 0)
      return 0;
    const double h = height(d);
    return 1 / (1 + std::sqrt(std::max(0.0, h) / d));
  }

  //! The angle of the section's arc from the wall, at the distance d.
  double angle(double d) const { return std::acos(std::clamp(1 - d / r, -1.0, 1.0)); }

  gp_Pnt2d uv(const Element& e, double t, double d) const
  {
    if (!e.straight)
      return gp_Pnt2d(t, across(d));
    if (bevel)
      return gp_Pnt2d(t, d * M_SQRT2);
    return gp_Pnt2d(e.flipped ? M_PI / 2 - angle(d) : angle(d), t);
  }

  // ---- Reach ----
  //! How far along its normal the place `f` gets before the place `x` of `o` is as near.
  double stopped(const Frame& f, const Element& o, double x) const
  {
    const gp_Vec2d to(f.p, o.curve.Value(x));
    const double   acrossIt = to.Dot(f.n);
    // Two places too near each other say nothing that rounding has not eaten: what stops
    // a place from that near is a bend, which is asked by itself.
    if (to.SquareMagnitude() < kNear * kNear * o.spacing * o.spacing)
      return HUGE_VAL;
    return acrossIt > 1e-14 * scale ? to.SquareMagnitude() / (2 * acrossIt) : HUGE_VAL;
  }

  //! The place of `other` that stops the place `f` (of the edge `own` at `t`) first, near a
  //! guess, and how far out. It is where `stopped` is least: looked for as that, which holds
  //! where the two places run together and the equations of the meeting normals do not, and
  //! then made exact by those equations where they hold. A place does not stop itself: where
  //! the least is at the place itself, it is the edge's own bend that stops it, and this
  //! says no.
  bool partner(const Frame& f, int own, double t, int other, double& tp, double& d) const
  {
    const Element& o     = elements[other];
    const double   span  = o.b - o.a;
    double         least = o.a, most = o.b;
    if (own == other)
    {
      const double margin = kNear * span / double(std::max<size_t>(o.points.size(), 2) - 1);
      if (tp > t)
        least = t + margin;
      else
        most = t - margin;
      if (!(least < most))
        return false;
    }
    double step = 2 * span / double(std::max<size_t>(o.points.size(), 2) - 1);
    double lo = std::clamp(tp - step, least, most), hi = std::clamp(tp + step, least, most);
    double x = tp, value = HUGE_VAL;
    for (int widen = 0; widen < 40; ++widen)
    {
      double a = lo, b = hi;
      double x1 = b - 0.6180339887498949 * (b - a), x2 = a + 0.6180339887498949 * (b - a);
      double g1 = stopped(f, o, x1), g2 = stopped(f, o, x2);
      for (int k = 0; k < 70 && b - a > 1e-15 * span; ++k)
      {
        if (g1 <= g2)
        {
          b  = x2;
          x2 = x1;
          g2 = g1;
          x1 = b - 0.6180339887498949 * (b - a);
          g1 = stopped(f, o, x1);
        }
        else
        {
          a  = x1;
          x1 = x2;
          g1 = g2;
          x2 = a + 0.6180339887498949 * (b - a);
          g2 = stopped(f, o, x2);
        }
      }
      x     = (a + b) / 2;
      value = stopped(f, o, x);
      // At an end of the window that is not an end of what may be asked: look further.
      const double edge = 1e-6 * (hi - lo);
      if (x - lo < edge && lo > least)
      {
        hi = lo + (hi - lo) / 2;
        lo = std::max(least, lo - step);
      }
      else if (hi - x < edge && hi < most)
      {
        lo = hi - (hi - lo) / 2;
        hi = std::min(most, hi + step);
      }
      else
        break;
      step *= 2;
    }
    if (!std::isfinite(value) || value == HUGE_VAL)
      return false;
    if (own == other
        && (x - least < BRepExt::tolerance::parameterMerge * span
            || most - x < BRepExt::tolerance::parameterMerge * span))
      return false;
    if (f.p.Distance(o.curve.Value(x)) < 2 * kNear * o.spacing)
      return false;
    tp = x;
    d  = value;
    // Exactly: p + d n = q + d m, for the place q of the other edge with its normal m.
    double xt = tp, xd = d;
    for (int iteration = 0; iteration < 12; ++iteration)
    {
      const Frame    g   = o.frame(xt);
      const gp_Vec2d now = gp_Vec2d(g.p, f.p) + (f.n - g.n) * xd;
      const gp_Vec2d dx = g.along * (-g.curveSpeed) - g.dn * xd, dd = f.n - g.n;
      const double   det = dx.X() * dd.Y() - dx.Y() * dd.X();
      if (std::abs(det) < 1e-12 * std::max(1.0, g.curveSpeed))
        break;
      const double stepX = (-now.X() * dd.Y() + now.Y() * dd.X()) / det;
      const double stepD = (-dx.X() * now.Y() + dx.Y() * now.X()) / det;
      xt += stepX;
      xd += stepD;
      if (!(std::abs(xt - tp) < 1e-4 * span) || !(std::abs(xd - d) < 1e-6 * scale) || xt < least
          || xt > most)
        break;
      const Frame solved = o.frame(xt);
      if (gp_Vec2d(solved.p, f.p).Added((f.n - solved.n) * xd).Magnitude() < 1e-12 * scale)
      {
        tp = xt;
        d  = xd;
        break;
      }
    }
    return true;
  }

  //! The reach of the place `t` of the edge `own`, and what stops it.
  Hit reach(int own, double t) const
  {
    const Element& e = elements[own];
    const Frame    f = e.frame(t);
    Hit            best{r, -1, 0};
    if (f.bend > 1e-14 && 1 / f.bend < best.d)
      best = {1 / f.bend, own, t};
    // At an end where the outline goes on smoothly, the edge that goes on may bend more
    // tightly there than this one: its bend stops the place as well.
    for (const bool last : {false, true})
    {
      const int next = last ? e.after : e.before;
      if (next < 0 || next == own || (last ? e.sharpB : e.sharpA)
          || std::abs(t - (last ? e.b : e.a)) > BRepExt::tolerance::parameter * (e.b - e.a))
        continue;
      const double at   = last ? e.afterAt : e.beforeAt;
      const double bend = elements[next].frame(at).bend;
      if (bend > 1e-14 && 1 / bend < best.d)
        best = {1 / bend, next, at};
    }
    // The places of its own edge next to it, on either side: where the edge bends ever more
    // tightly one way, a place a little further that way stops this one before its own
    // bend does. The polyline is too coarse to say so.
    if (e.kind != GeomAbs_Line && e.kind != GeomAbs_Circle)
      for (const double way : {1.0, -1.0})
      {
        double tp = t + way * 3 * kNear * (e.b - e.a) / double(e.points.size() - 1), d = best.d;
        if (tp <= e.a || tp >= e.b)
          continue;
        if (partner(f, own, t, own, tp, d) && d < best.d - 1e-12 * scale)
          best = {std::max(d, 0.0), own, tp};
      }
    for (int other = 0; other < int(elements.size()); ++other)
    {
      const Element& o = elements[other];
      // A straight edge has no place that stops another of its own, and all the places
      // of an arc stop each other at its centre, which its bend has said already.
      if (other == own && (o.kind == GeomAbs_Line || o.kind == GeomAbs_Circle))
        continue;
      const int  count  = int(o.points.size());
      double     before = HUGE_VAL, here = HUGE_VAL;
      const auto g = [&](int k) {
        const gp_Vec2d to(f.p, o.points[k]);
        const double   across = to.Dot(f.n);
        if (across <= 1e-12 * scale)
          return HUGE_VAL;
        return to.SquareMagnitude() / (2 * across);
      };
      here = g(0);
      for (int k = 0; k < count; ++k)
      {
        const double next = k + 1 < count ? g(k + 1) : HUGE_VAL;
        const bool   low  = here <= before && here <= next && here < 1.05 * best.d + 2 * o.spacing;
        if (low && !(other == own && std::abs(o.ts[k] - t) < 1.5 * (o.b - o.a) / (count - 1)))
        {
          double tp = o.ts[k], d = std::min(here, r);
          // An arc's centre is where its own places stop each other; a place of
          // another edge that is as far (the end of a line going on from the arc)
          // has to be nearer by more than rounding to stand instead.
          const bool apex = best.by == own && best.t == t && e.kind == GeomAbs_Circle;
          if (partner(f, own, t, other, tp, d) && d < best.d - (apex ? 1e-7 : 1e-12) * scale)
          {
            best = {std::max(d, 0.0), other, std::clamp(tp, o.a, o.b)};
          }
        }
        before = here;
        here   = next;
      }
    }
    return best;
  }
};

//! A corner of the blend's faces.
struct Corners
{
  struct Corner
  {
    gp_Pnt        p;
    TopoDS_Vertex v;
  };

  std::vector<Corner> all;
  double              near = BRepExt::tolerance::sewing;

  int add(const gp_Pnt& p)
  {
    for (size_t i = 0; i < all.size(); ++i)
      if (all[i].p.Distance(p) < near)
        return int(i);
    all.push_back({p, BRepBuilderAPI_MakeVertex(p).Vertex()});
    return int(all.size()) - 1;
  }

  int add(const TopoDS_Vertex& v)
  {
    const gp_Pnt p = BRep_Tool::Pnt(v);
    for (size_t i = 0; i < all.size(); ++i)
      if (all[i].v.IsSame(v))
        return int(i);
    all.push_back({p, v});
    return int(all.size()) - 1;
  }
};

//! A piece of the top of an edge's face, between two places of the edge.
struct Piece
{
  enum Kind
  {
    Cap,
    Crease,
    Apex
  } kind = Cap;

  double from = 0, to = 0; // of the edge's parameter, from < to
  double dFrom = 0, dTo = 0;
  int    with     = -1; // a crease: the edge whose face it is shared with
  double withFrom = 0, withTo = 0;
  int    cornerFrom = -1, cornerTo = -1;
};

struct Maker : Blend
{
  Corners                         corners;
  BRep_Builder                    builder;
  std::vector<std::vector<Piece>> chains;

  struct Shared
  {
    TopoDS_Edge edge;
    int         builtBy = -1;
    bool        opposed = true; // the other face's edge runs against the one it was built from
  };

  std::map<std::tuple<int, int, int, int>, Shared> creases;

  struct Rise
  {
    TopoDS_Edge edge;
    int         top = -1;
    double      d   = 0;
  };

  NCollection_IndexedDataMap<TopoDS_Shape, int, TopTools_ShapeMapHasher> riseOf;
  std::vector<Rise>                                                      rises;

  struct CapEdge
  {
    TopoDS_Edge           edge;
    int                   from, to; // corners along the way the face lies to the left
    Handle(Geom2d_Curve)  inPlane;
    std::vector<gp_Pnt2d> polygon;
  };

  std::vector<CapEdge> capEdges;

  [[noreturn]] void refuse(const std::string& why) const { fail(op, no_solution, why); }

  [[noreturn]] void refuse(const std::string& why, const TopoDS_Shape& at) const
  {
    fail(op, no_solution, why, {at});
  }

  // ---- Step 2: the outline ----
  void addElements(const TopoDS_Face& face)
  {
    const Handle(Geom_Surface)       surface = BRep_Tool::Surface(face);
    std::vector<std::pair<int, int>> loops; // first element and count of each wire
    for (TopExp_Explorer wires(face, TopAbs_WIRE); wires.More(); wires.Next())
    {
      const int first = int(elements.size());
      for (BRepTools_WireExplorer walk(TopoDS::Wire(wires.Current()), face); walk.More();
           walk.Next())
      {
        Element           e;
        const TopoDS_Edge edge = walk.Current();
        if (BRep_Tool::Degenerated(edge))
          continue;
        e.s    = edge.Orientation() == TopAbs_REVERSED ? -1 : 1;
        e.edge = TopoDS::Edge(edge.Oriented(TopAbs_FORWARD));
        double     a, b;
        const auto curve = BRep_Tool::CurveOnSurface(e.edge, face, a, b);
        if (curve.IsNull())
          refuse("an edge of the outline has no curve in its plane", edge);
        e.curve.Load(curve, a, b);
        e.a    = a;
        e.b    = b;
        e.kind = e.curve.GetType();
        elements.push_back(e);
      }
      loops.emplace_back(first, int(elements.size()) - first);
    }
    // The wire explorer walks the way the face lies to the left; an edge that is reversed
    // in it is walked from b to a.
    for (const auto& [first, count] : loops)
    {
      for (int k = 0; k < count; ++k)
      {
        const int      i = first + k, j = first + (k + 1) % count;
        Element&       e     = elements[i];
        Element&       n     = elements[j];
        const double   leave = e.s > 0 ? e.b : e.a, arrive = n.s > 0 ? n.a : n.b;
        const gp_Vec2d out = e.frame(leave).along * e.s, in = n.frame(arrive).along * n.s;
        const double   turn = out.X() * in.Y() - out.Y() * in.X();
        if (turn < -BRepExt::tolerance::coarseAngle)
          refuse("the outline has a corner that turns away from the face: a blend around such a "
                 "corner is not built",
                 e.edge);
        const bool sharp                   = turn > BRepExt::tolerance::coarseAngle;
        (e.s > 0 ? e.after : e.before)     = j;
        (e.s > 0 ? e.afterAt : e.beforeAt) = arrive;
        (e.s > 0 ? e.sharpB : e.sharpA)    = sharp;
        (n.s > 0 ? n.before : n.after)     = i;
        (n.s > 0 ? n.beforeAt : n.afterAt) = leave;
        (n.s > 0 ? n.sharpA : n.sharpB)    = sharp;
      }
    }
    for (auto& e : elements)
    {
      double   length = 0;
      gp_Pnt2d last   = e.curve.Value(e.a);
      for (int k = 1; k <= 64; ++k)
      {
        const gp_Pnt2d p = e.curve.Value(e.a + (e.b - e.a) * k / 64);
        length += p.Distance(last);
        last = p;
      }
      const int count = std::clamp(int(length / (r / 24)) + 1, 16, 4096);
      e.spacing       = length / count;
      for (int k = 0; k <= count; ++k)
      {
        e.ts.push_back(e.a + (e.b - e.a) * k / count);
        e.points.push_back(e.curve.Value(e.ts.back()));
      }
    }
  }

  // ---- Step 2: the reach along every edge ----
  Hit atEnd(int own, bool last) const
  {
    const Element& e = elements[own];
    if (last ? e.sharpB : e.sharpA)
      return {0, last ? e.after : e.before, last ? e.afterAt : e.beforeAt};
    return reach(own, last ? e.b : e.a);
  }

  gp_Pnt2d placeOf(int own, double t, const Hit& hit) const
  {
    if (hit.by < 0)
      return gp_Pnt2d(1e300, 1e300);
    if (hit.by == own && hit.t == t)
      return elements[own].frame(t).p;
    return elements[hit.by].curve.Value(hit.t);
  }

  bool alike(int own, double ta, const Hit& a, double tb, const Hit& b) const
  {
    if (a.by != b.by)
      return false;
    if (a.by < 0)
      return true;
    const bool selfA = a.by == own && a.t == ta, selfB = b.by == own && b.t == tb;
    // Next to where an edge's bend closes its face, its own bend and the place beside it
    // that stops it are one stop: which is found there is rounding.
    if (selfA != selfB)
      return a.by == own && std::abs(a.d - b.d) < 1e-4 * scale;
    if (selfA)
      return true;
    const Element& e     = elements[own];
    const double   moved = e.curve.Value(ta).Distance(e.curve.Value(tb)) + std::abs(a.d - b.d);
    return placeOf(own, ta, a).Distance(placeOf(own, tb, b)) <= 64 * moved + 1e-6 * scale;
  }

  struct Event
  {
    double t;
    Hit    left, right;
    bool   leaf = false;
  };

  void resolve(int own, double ta, Hit a, double tb, Hit b, std::vector<Event>& events, int depth)
    const
  {
    const Element& e = elements[own];
    if (alike(own, ta, a, tb, b))
    {
      // Its own bend closes the edge's face where what stops it passes it by.
      if (a.by == own && b.by == own && a.t != ta && b.t != tb && a.t > ta && b.t < tb)
      {
        double lo = ta, hi = tb;
        for (int k = 0; k < 80; ++k)
        {
          const double x = lo + (hi - lo) * 0.381966, y = lo + (hi - lo) * 0.618034;
          if (e.frame(x).bend < e.frame(y).bend)
            lo = x;
          else
            hi = y;
        }
        const double t = (lo + hi) / 2, bend = e.frame(t).bend;
        if (bend > 1e-14)
          events.push_back({t, {1 / bend, own, t}, {1 / bend, own, t}, true});
      }
      return;
    }
    if (depth > 12)
      refuse("the outline is too intricate for the blend: what stops its places changes too often",
             e.edge);
    while (tb - ta > 1e-13 * (e.b - e.a))
    {
      const double tm   = (ta + tb) / 2;
      const Hit    m    = reach(own, tm);
      const bool   left = alike(own, ta, a, tm, m), right = alike(own, tm, m, tb, b);
      if (left && right)
        return; // it only moved fast
      if (left)
      {
        ta = tm;
        a  = m;
      }
      else if (right)
      {
        tb = tm;
        b  = m;
      }
      else
      {
        resolve(own, ta, a, tm, m, events, depth + 1);
        resolve(own, tm, m, tb, b, events, depth + 1);
        return;
      }
    }
    if (a.by == b.by && placeOf(own, ta, a).Distance(placeOf(own, tb, b)) < 1e-6 * scale)
      return;
    events.push_back({(ta + tb) / 2, a, b, false});
  }

  void chainOf(int own)
  {
    const Element&     e = elements[own];
    std::vector<Event> events;
    // Places looked at: evenly along the edge, and every place where the edge bends more
    // tightly than the blend, which a crease must reach however short it is.
    std::vector<double> ts;
    const double        inside = 5 * kNear * (e.b - e.a) / double(e.points.size() - 1);
    for (int k = 0; k <= kScan; ++k)
      ts.push_back(e.a + (e.b - e.a) * k / kScan);
    if (e.kind != GeomAbs_Line && e.kind != GeomAbs_Circle)
    {
      const int           fine = 16 * kScan;
      std::vector<double> bends(fine + 1);
      for (int k = 0; k <= fine; ++k)
        bends[k] = e.frame(e.a + (e.b - e.a) * k / fine).bend;
      for (int k = 1; k < fine; ++k)
      {
        if (bends[k] < bends[k - 1] || bends[k] < bends[k + 1])
          continue;
        double lo = e.a + (e.b - e.a) * (k - 1) / fine, hi = e.a + (e.b - e.a) * (k + 1) / fine;
        for (int step = 0; step < 60; ++step)
        {
          const double x = lo + (hi - lo) * 0.381966, y = lo + (hi - lo) * 0.618034;
          if (e.frame(x).bend < e.frame(y).bend)
            lo = x;
          else
            hi = y;
        }
        const double t = (lo + hi) / 2;
        if (e.frame(t).bend * r > 1 && t > e.a + inside && t < e.b - inside)
          ts.push_back(t);
      }
      std::sort(ts.begin(), ts.end());
    }
    std::vector<Hit> hits(ts.size());
    for (size_t k = 0; k < ts.size(); ++k)
      hits[k] = k == 0               ? atEnd(own, false)
                : k + 1 == ts.size() ? atEnd(own, true)
                                     : reach(own, ts[k]);
    for (size_t k = 0; k + 1 < ts.size(); ++k)
    {
      // The ends are known by rule; what is between is asked a little inside them.
      double ta = ts[k], tb = ts[k + 1];
      Hit    a = hits[k], b = hits[k + 1];
      if (k == 0)
      {
        ta = e.a + inside;
        a  = reach(own, ta);
      }
      if (k + 2 == ts.size())
      {
        tb = e.b - inside;
        b  = reach(own, tb);
      }
      if (tb <= ta)
        continue;
      resolve(own, ta, a, tb, b, events, 0);
    }
    std::sort(events.begin(), events.end(), [](const Event& x, const Event& y) {
      return x.t < y.t;
    });
    // Events found at one place (where rounding makes what stops it flicker) are one: what
    // stops the places before it and what stops those after.
    {
      std::vector<Event> merged;
      for (const auto& event : events)
      {
        if (!merged.empty()
            && event.t - merged.back().t < BRepExt::tolerance::parameterMerge * (e.b - e.a)
            && !event.leaf && !merged.back().leaf)
        {
          merged.back().right = event.right;
          continue;
        }
        merged.push_back(event);
      }
      events.clear();
      for (const auto& event : merged)
        if (event.leaf || !alike(own, event.t, event.left, event.t, event.right))
          events.push_back(event);
    }
    // The pieces between the events.
    // Where a crease comes up to the face's plane, the outline moved in crosses itself:
    // that place exactly, by Newton's method from where bisection left it.
    for (auto& event : events)
    {
      if (event.leaf || (event.left.by < 0) == (event.right.by < 0))
        continue;
      Hit&           stop = event.left.by < 0 ? event.right : event.left;
      const Element& o    = elements[stop.by];
      double         t = event.t, tp = stop.t;
      for (int iteration = 0; iteration < 12; ++iteration)
      {
        const Frame    f = e.frame(t), g = o.frame(tp);
        const gp_Vec2d now = gp_Vec2d(g.p, f.p) + (f.n - g.n) * r;
        const gp_Vec2d dt  = f.along * f.curveSpeed + f.dn * r,
                       dp  = (g.along * g.curveSpeed + g.dn * r) * -1.0;
        const double det   = dt.X() * dp.Y() - dt.Y() * dp.X();
        if (std::abs(det) < 1e-300)
          break;
        const double stepT = (-now.X() * dp.Y() + now.Y() * dp.X()) / det,
                     stepP = (-dt.X() * now.Y() + dt.Y() * now.X()) / det;
        t += stepT;
        tp += stepP;
        if (std::abs(t - event.t) > 1e-6 * (e.b - e.a)
            || std::abs(tp - stop.t) > 1e-6 * (o.b - o.a))
          break;
        const Frame solved = e.frame(t), other = o.frame(tp);
        if (gp_Vec2d(other.p, solved.p).Added((solved.n - other.n) * r).Magnitude() < 1e-12 * scale)
        {
          event.t = t;
          stop.t  = tp;
          stop.d  = r;
          break;
        }
      }
    }
    // A stretch whose places stop each other from its two ends has, between them, the
    // place where the edge's bend closes its face, also where it is too short for two of
    // the places looked at to have fallen in it.
    {
      std::vector<double> at{e.a};
      for (const auto& event : events)
        at.push_back(event.t);
      at.push_back(e.b);
      std::vector<Event> leaves;
      for (size_t k = 0; k + 1 < at.size(); ++k)
      {
        if ((k > 0 && events[k - 1].leaf) || (k < events.size() && events[k].leaf))
          continue;
        const double middle = (at[k] + at[k + 1]) / 2;
        const Hit    hit    = reach(own, middle);
        if (hit.by != own || hit.t == middle || hit.t <= at[k] || hit.t >= at[k + 1])
          continue;
        double lo = std::min(middle, hit.t), hi = std::max(middle, hit.t);
        for (int step = 0; step < 80; ++step)
        {
          const double x = lo + (hi - lo) * 0.381966, y = lo + (hi - lo) * 0.618034;
          if (e.frame(x).bend < e.frame(y).bend)
            lo = x;
          else
            hi = y;
        }
        const double t = (lo + hi) / 2, bend = e.frame(t).bend;
        if (bend > 1e-14)
          leaves.push_back({t, {1 / bend, own, t}, {1 / bend, own, t}, true});
      }
      for (const auto& leaf : leaves)
        events.push_back(leaf);
      std::sort(events.begin(), events.end(), [](const Event& x, const Event& y) {
        return x.t < y.t;
      });
    }
    std::vector<double> breaks{e.a};
    std::vector<double> reachAt{atEnd(own, false).d};
    for (const auto& event : events)
    {
      breaks.push_back(event.t);
      reachAt.push_back(std::min(event.left.d, event.right.d));
    }
    breaks.push_back(e.b);
    reachAt.push_back(atEnd(own, true).d);
    std::vector<Piece>& chain = chains[own];
    for (size_t k = 0; k + 1 < breaks.size(); ++k)
    {
      Piece piece;
      piece.from = breaks[k];
      piece.to   = breaks[k + 1];
      if (piece.to - piece.from < BRepExt::tolerance::parameter * (e.b - e.a))
        continue;
      double middle = (piece.from + piece.to) / 2;
      Hit    hit    = reach(own, middle);
      // A curve that is not an arc is stopped by its own bend at single places only; the
      // middle of a piece that is one is asked again a little to either side.
      for (const double share : {0.3, 0.7, 0.1, 0.9})
      {
        if (e.kind == GeomAbs_Circle || hit.by != own || hit.t != middle)
          break;
        middle = piece.from + (piece.to - piece.from) * share;
        hit    = reach(own, middle);
      }
      piece.dFrom = reachAt[k];
      piece.dTo   = reachAt[k + 1];
      if (hit.by < 0)
      {
        piece.kind  = Piece::Cap;
        piece.dFrom = piece.dTo = r;
      }
      else if (hit.by == own && hit.t == middle)
      {
        if (e.kind != GeomAbs_Circle)
          refuse("an edge of the outline is stopped by its own bend along a stretch, as only an "
                 "arc can be",
                 e.edge);
        piece.kind  = Piece::Apex;
        piece.dFrom = piece.dTo = hit.d;
      }
      else
      {
        piece.kind = Piece::Crease;
        piece.with = hit.by;
        // The partner's places at both ends, by following it out from the middle.
        const auto follow =
          [&](double to, bool leafEnd, bool cornerEnd, double cornerAt, double& where, double& d) {
            if (leafEnd)
            {
              where = to;
              return;
            }
            if (cornerEnd)
            {
              where = cornerAt;
              d     = 0;
              return;
            }
            double    tp = hit.t, dist = hit.d;
            const int steps = 24;
            for (int step = 1; step < steps; ++step)
            {
              const double t = middle + (to - middle) * step / steps;
              // Where it cannot be asked (two places that run together), the place
              // before stands; the end is come to below.
              double nextT = tp, nextD = dist;
              if (partner(e.frame(t), own, t, hit.by, nextT, nextD))
              {
                tp   = nextT;
                dist = nextD;
              }
            }
            // The end itself may be a place where nothing can be asked (three places at
            // one distance): it is come to by ever smaller steps.
            double left = (to - middle) / steps;
            for (int step = 0; step < 40; ++step)
            {
              left /= 2;
              double nextT = tp, nextD = dist;
              if (!partner(e.frame(to - left), own, to - left, hit.by, nextT, nextD))
                break;
              tp   = nextT;
              dist = nextD;
            }
            where = std::clamp(tp, elements[hit.by].a, elements[hit.by].b);
            // Where the outline goes on smoothly into the edge on the other side of the
            // crease, the crease ends over that corner, at both edges' ends.
            const bool last = to == e.b;
            if ((to == e.a || last) && (last ? e.after : e.before) == hit.by
                && !(last ? e.sharpB : e.sharpA))
            {
              const double corner = last ? e.afterAt : e.beforeAt;
              if (std::abs(where - corner) < 1e-3 * (elements[hit.by].b - elements[hit.by].a))
                where = corner;
            }
          };
        const bool leafFrom   = k > 0 && events[k - 1].leaf,
                   leafTo     = k < events.size() && events[k].leaf;
        const bool cornerFrom = k == 0 && e.sharpA, cornerTo = k + 2 == breaks.size() && e.sharpB;
        follow(piece.from, leafFrom, cornerFrom, e.beforeAt, piece.withFrom, piece.dFrom);
        follow(piece.to, leafTo, cornerTo, e.afterAt, piece.withTo, piece.dTo);
        if (leafFrom)
          piece.withFrom = piece.from;
        if (leafTo)
          piece.withTo = piece.to;
      }
      piece.cornerFrom = corners.add(at(e, piece.from, piece.dFrom));
      piece.cornerTo   = corners.add(at(e, piece.to, piece.dTo));
      // A piece that begins and ends in one corner is nothing, unless it goes round.
      if (piece.kind != Piece::Apex && piece.cornerFrom == piece.cornerTo
          && at(e, middle, hit.d).Distance(corners.all[piece.cornerFrom].p) < corners.near)
        continue;
      chain.push_back(piece);
    }
    if (chain.empty())
      refuse("an edge of the outline has no face in the blend", e.edge);
  }

  // ---- Step 3: the surfaces ----
  void surfaceOf(Element& e)
  {
    const gp_Vec2d moveIn = e.frame((e.a + e.b) / 2).n;
    if (e.kind == GeomAbs_Line)
    {
      const gp_Lin2d line   = e.curve.Line();
      const gp_Vec   along  = X * line.Direction().X() + Y * line.Direction().Y();
      const gp_Vec   inward = X * moveIn.X() + Y * moveIn.Y();
      e.straight            = true;
      const gp_Pnt2d origin = line.Location();
      if (bevel)
      {
        const gp_Vec rise = (inward - w).Normalized();
        e.surface =
          new Geom_Plane(gp_Ax3(lift(origin, r), gp_Dir(along.Crossed(rise)), gp_Dir(along)));
      }
      else
      {
        e.flipped = w.Dot(along.Crossed(inward)) < 0;
        e.surface = new Geom_CylindricalSurface(gp_Ax3(lift(origin.Translated(moveIn * r), r),
                                                       gp_Dir(along),
                                                       gp_Dir(e.flipped ? -w : -inward)),
                                                r);
      }
      e.moved        = new Geom2d_Line(origin.Translated(moveIn * r), line.Direction());
      e.movedInSpace = new Geom_Line(lift(origin.Translated(moveIn * r), 0), gp_Dir(along));
      return;
    }
    // Spans: where the curve is smooth, halved until a polynomial of kDegree fits both the
    // outline and the outline moved in.
    static const Fit    fit;
    std::vector<double> breaks;
    {
      const int                  count = e.curve.NbIntervals(GeomAbs_CN);
      NCollection_Array1<double> intervals(1, count + 1);
      e.curve.Intervals(intervals, GeomAbs_CN);
      for (int k = 1; k <= count + 1; ++k)
        breaks.push_back(std::clamp(intervals(k), e.a, e.b));
      if (e.kind == GeomAbs_Circle || e.kind == GeomAbs_Ellipse)
      {
        const int pieces = int(std::ceil((e.b - e.a) / (M_PI / 4)));
        breaks.clear();
        for (int k = 0; k <= pieces; ++k)
          breaks.push_back(e.a + (e.b - e.a) * k / pieces);
      }
    }
    using Poles      = std::array<std::array<double, kDegree + 1>, 4>;
    const auto exact = [&](double t, double out[4]) {
      const Frame f = e.frame(t);
      out[0]        = f.p.X();
      out[1]        = f.p.Y();
      out[2]        = f.p.X() + r * f.n.X();
      out[3]        = f.p.Y() + r * f.n.Y();
    };
    std::vector<double>                            knots;
    std::vector<Poles>                             spans;
    const std::function<void(double, double, int)> span = [&](double from, double to, int depth) {
      Poles  poles{};
      double values[kDegree + 1][4];
      for (int k = 0; k <= kDegree; ++k)
        exact(from + (to - from) * fit.nodes[k], values[k]);
      for (int c = 0; c < 4; ++c)
        for (int j = 0; j <= kDegree; ++j)
        {
          double sum = 0;
          for (int k = 0; k <= kDegree; ++k)
            sum += fit.inverse[j][k] * values[k][c];
          poles[c][j] = sum;
        }
      double worst = 0;
      for (int k = 0; k < kDegree; ++k)
      {
        const double s = (fit.nodes[k] + fit.nodes[k + 1]) / 2;
        double       value[4];
        exact(from + (to - from) * s, value);
        for (int c = 0; c < 4; ++c)
          worst = std::max(worst, std::abs(bezierAt(poles[c], s) - value[c]));
      }
      if (worst > kFit * scale && depth < 12)
      {
        span(from, (from + to) / 2, depth + 1);
        span((from + to) / 2, to, depth + 1);
        return;
      }
      knots.push_back(from);
      spans.push_back(poles);
    };
    for (size_t k = 0; k + 1 < breaks.size(); ++k)
      if (breaks[k + 1] - breaks[k] > BRepExt::tolerance::parameter * (e.b - e.a))
        span(breaks[k], breaks[k + 1], 0);
    knots.push_back(e.b);
    const int count = int(spans.size()), poles = count * kDegree + 1, rows = bevel ? 2 : 3;
    NCollection_Array1<double> uKnots(1, count + 1), vKnots(1, 2);
    NCollection_Array1<int>    uMults(1, count + 1), vMults(1, 2);
    for (int k = 0; k <= count; ++k)
    {
      uKnots(k + 1) = knots[k];
      uMults(k + 1) = k == 0 || k == count ? kDegree + 1 : kDegree;
    }
    vKnots(1) = 0;
    vKnots(2) = 1;
    vMults(1) = vMults(2) = rows;
    NCollection_Array2<gp_Pnt>   grid(1, poles, 1, rows);
    NCollection_Array2<double>   weights(1, poles, 1, rows);
    NCollection_Array1<gp_Pnt2d> movedPoles(1, poles);
    for (int k = 0; k < count; ++k)
      for (int j = 0; j <= kDegree; ++j)
      {
        const int      index = k * kDegree + j + 1;
        const gp_Pnt2d on(spans[k][0][j], spans[k][1][j]), in(spans[k][2][j], spans[k][3][j]);
        movedPoles(index) = in;
        grid(index, 1)    = lift(on, r);
        weights(index, 1) = 1;
        if (!bevel)
        {
          grid(index, 2)    = lift(on, 0);
          weights(index, 2) = M_SQRT1_2;
        }
        grid(index, rows)    = lift(in, 0);
        weights(index, rows) = 1;
      }
    e.surface =
      new Geom_BSplineSurface(grid, weights, uKnots, vKnots, uMults, vMults, kDegree, rows - 1);
    e.moved        = new Geom2d_BSplineCurve(movedPoles, uKnots, uMults, kDegree);
    e.movedInSpace = Handle(Geom_BSplineSurface)::DownCast(e.surface)->VIso(1);
  }

  // ---- Step 3: edges ----
  //! An edge on `curve` between two corners, whose tolerances are raised to cover how far
  //! the curve's ends are from them: a corner is found from two sides, a little apart.
  TopoDS_Edge joined(const Handle(Geom_Curve)& curve,
                     const TopoDS_Vertex&      from,
                     const TopoDS_Vertex&      to,
                     double                    a,
                     double                    b,
                     const char*               what)
  {
    const double gapA = curve->Value(a).Distance(BRep_Tool::Pnt(from)),
                 gapB = curve->Value(b).Distance(BRep_Tool::Pnt(to));
    if (std::max(gapA, gapB) > 1e-4 * scale)
      refuse(what);
    TopoDS_Edge edge;
    builder.MakeEdge(edge, curve, BRepExt::tolerance::length);
    builder.Add(edge, from.Oriented(TopAbs_FORWARD));
    builder.Add(edge, to.Oriented(TopAbs_REVERSED));
    builder.Range(edge, a, b);
    builder.UpdateVertex(from, std::max(BRep_Tool::Tolerance(from), 1.5 * gapA));
    builder.UpdateVertex(to, std::max(BRep_Tool::Tolerance(to), 1.5 * gapB));
    return edge;
  }

  //! Give an edge the largest tolerance that its curve in a face asks for.
  void fitTolerance(const TopoDS_Edge&          edge,
                    const Handle(Geom2d_Curve)& onFace,
                    const Handle(Geom_Surface)& surface)
  {
    double     a, b;
    const auto curve = BRep_Tool::Curve(edge, a, b);
    if (curve.IsNull())
      return;
    double worst = 0;
    // At a multiple of the places BRepCheck_Edge looks at.
    for (int k = 0; k <= 176; ++k)
    {
      const double   t  = a + (b - a) * k / 176;
      const gp_Pnt2d uv = onFace->Value(t);
      worst             = std::max(worst, curve->Value(t).Distance(surface->Value(uv.X(), uv.Y())));
    }
    builder.UpdateEdge(edge, std::max(BRepExt::tolerance::length, 2 * worst));
  }

  void fitVertices(const TopoDS_Edge& edge)
  {
    const double tolerance = BRep_Tool::Tolerance(edge);
    double       a, b;
    const auto   curve = BRep_Tool::Curve(edge, a, b);
    for (TopExp_Explorer it(edge, TopAbs_VERTEX); it.More(); it.Next())
    {
      const auto vertex = TopoDS::Vertex(it.Current());
      double     off    = 0;
      if (!curve.IsNull())
      {
        const gp_Pnt p = BRep_Tool::Pnt(vertex);
        off            = std::min(p.Distance(curve->Value(a)), p.Distance(curve->Value(b)));
      }
      builder.UpdateVertex(vertex, std::max(tolerance, 1.5 * off));
    }
  }

  Handle(Geom2d_Curve) through(const std::vector<gp_Pnt2d>& points,
                               const std::vector<double>&   params,
                               bool                         cuspFirst = false,
                               bool                         cuspLast  = false) const
  {
    Handle(NCollection_HArray1<gp_Pnt2d>) array =
      new NCollection_HArray1<gp_Pnt2d>(1, int(points.size()));
    Handle(NCollection_HArray1<double>) at = new NCollection_HArray1<double>(1, int(points.size()));
    for (size_t k = 0; k < points.size(); ++k)
    {
      array->SetValue(int(k) + 1, points[k]);
      at->SetValue(int(k) + 1, params[k]);
    }
    Geom2dAPI_Interpolate interpolate(array, at, false, 1e-12);
    interpolate.Perform();
    if (!interpolate.IsDone())
      refuse("a curve of the blend could not be put through its points");
    const auto curve = interpolate.Curve();
    // At a self-crease's cusp the two outline parameters merge. The boundary must
    // stay on its own side of the split through that parameter. Interpolation can
    // overshoot that coordinate and cross the split before the vertex.
    // The convex hull of the poles keeps that coordinate on the prescribed side.
    for (int k = 1; k <= curve->NbPoles(); ++k)
    {
      gp_Pnt2d pole = curve->Pole(k);
      if (cuspFirst)
        pole.SetX(points[1].X() > points[0].X() ? std::max(pole.X(), points[0].X())
                                                : std::min(pole.X(), points[0].X()));
      if (cuspLast)
        pole.SetX(points.back().X() > points[points.size() - 2].X()
                    ? std::min(pole.X(), points.back().X())
                    : std::max(pole.X(), points.back().X()));
      if (cuspFirst)
        pole.SetY(std::max(pole.Y(), points.front().Y()));
      if (cuspLast)
        pole.SetY(std::max(pole.Y(), points.back().Y()));
      curve->SetPole(k, pole);
    }
    return Handle(Geom2d_Curve)(curve);
  }

  //! The section at an end of an edge, from the foot up to where the edge's face ends
  //! there; shared with the edge that goes on from that end.
  int riseAt(int own, bool last)
  {
    const Element& e = elements[own];
    const double   t = last ? e.b : e.a;
    TopoDS_Vertex  first, second;
    TopExp::Vertices(e.edge, first, second);
    const TopoDS_Vertex base = last ? second : first;
    if (riseOf.Contains(base))
      return riseOf.FindFromKey(base);
    Rise      rise;
    const Hit hit = atEnd(own, last);
    rise.d        = hit.d;
    if (hit.d < 1e-9 * scale)
    {
      rise.top = corners.add(base);
    }
    else
    {
      rise.top                  = corners.add(at(e, t, hit.d));
      const Frame        f      = e.frame(t);
      const gp_Vec       inward = X * f.n.X() + Y * f.n.Y();
      Handle(Geom_Curve) curve;
      double             end;
      if (bevel)
      {
        curve = new Geom_Line(lift(f.p, r), gp_Dir(inward - w));
        end   = hit.d * M_SQRT2;
      }
      else
      {
        curve = new Geom_Circle(
          gp_Ax2(lift(f.p.Translated(f.n * r), r), gp_Dir(inward.Crossed(w)), gp_Dir(-inward)),
          r);
        end = angle(hit.d);
      }
      rise.edge = joined(curve,
                         base,
                         corners.all[rise.top].v,
                         0,
                         end,
                         "a section of the blend could not be made");
    }
    rises.push_back(rise);
    riseOf.Add(base, int(rises.size()) - 1);
    return int(rises.size()) - 1;
  }

  //! The curve of a section in an edge's face.
  Handle(Geom2d_Curve) riseIn(const Element& e, double t, const Rise& rise) const
  {
    double a, b;
    BRep_Tool::Range(rise.edge, a, b);
    if (e.straight)
    {
      if (bevel)
        return new Geom2d_Line(gp_Pnt2d(t, 0), gp_Dir2d(0, 1));
      return new Geom2d_Line(gp_Pnt2d(e.flipped ? M_PI / 2 : 0, t),
                             gp_Dir2d(e.flipped ? -1 : 1, 0));
    }
    std::vector<gp_Pnt2d> points;
    std::vector<double>   params;
    const int             count = 64;
    for (int k = 0; k <= count; ++k)
    {
      const double x = a + (b - a) * k / count;
      const double d = bevel ? x / M_SQRT2 : r * (1 - std::cos(x));
      points.emplace_back(t, across(d));
      params.push_back(x);
    }
    return through(points, params);
  }

  //! The crease that a piece of `own` shares with another edge's face: made once.
  Shared& creaseOf(int own, const Piece& piece)
  {
    const auto key   = std::make_tuple(std::min(own, piece.with),
                                     std::max(own, piece.with),
                                     std::min(piece.cornerFrom, piece.cornerTo),
                                     std::max(piece.cornerFrom, piece.cornerTo));
    const auto found = creases.find(key);
    if (found != creases.end())
      return found->second;
    const Element& e = elements[own];
    const Element& o = elements[piece.with];

    // A crease leaves a wall along the section there, and ends where an edge's own bend
    // closes its face in a point where the face has no width: its curve in a face has no
    // speed at the one, the crease itself none at the other, by the edge's parameter. So it
    // is run by a length that counts both: how far it goes in space and how far its two
    // places on the outline go.
    struct Place
    {
      double t, with, d;
    };

    const Hit middle = reach(own, (piece.from + piece.to) / 2);
    if (middle.by != piece.with)
      refuse("a crease of the blend changes what it lies between", e.edge);
    //! The places at these parameters, from the middle outwards, each from the one before.
    const auto placesAt = [&](const std::vector<double>& ts) {
      const int          count = int(ts.size()) - 1;
      std::vector<Place> places(count + 1);
      places[0]     = {piece.from, piece.withFrom, piece.dFrom};
      places[count] = {piece.to, piece.withTo, piece.dTo};
      int centre    = 1;
      for (int k = 1; k < count; ++k)
        if (std::abs(ts[k] - (piece.from + piece.to) / 2)
            < std::abs(ts[centre] - (piece.from + piece.to) / 2))
          centre = k;
      for (int direction : {1, -1})
      {
        double tp = middle.t, d = middle.d;
        for (int k = direction > 0 ? centre : centre - 1; k > 0 && k < count; k += direction)
        {
          double nextT = tp, nextD = d;
          if (partner(e.frame(ts[k]), own, ts[k], piece.with, nextT, nextD))
          {
            tp = nextT;
            d  = nextD;
          }
          else
          {
            // Next to where its own bend closes the face the two places run together.
            const Place& end   = direction > 0 ? places[count] : places[0];
            const double share = (ts[k] - ts[k - direction]) / (end.t - ts[k - direction]);
            tp += (end.with - tp) * share;
            const gp_Vec2d to(e.frame(ts[k]).p, o.frame(tp).p);
            const double   acrossIt = to.Dot(e.frame(ts[k]).n);
            d = acrossIt > 1e-12 * scale ? to.SquareMagnitude() / (2 * acrossIt)
                                         : d + (end.d - d) * share;
          }
          places[k] = {ts[k], tp, d};
        }
      }
      return places;
    };
    const auto lengths = [&](const std::vector<Place>& places, std::vector<gp_Pnt>& points) {
      std::vector<double> run(places.size(), 0);
      points.clear();
      for (const auto& place : places)
        points.push_back(at(e, place.t, place.d));
      points.front() = corners.all[piece.cornerFrom].p;
      points.back()  = corners.all[piece.cornerTo].p;
      for (size_t k = 1; k < places.size(); ++k)
        run[k] = run[k - 1] + points[k].Distance(points[k - 1])
                 + e.curve.Value(places[k].t).Distance(e.curve.Value(places[k - 1].t))
                 + o.curve.Value(places[k].with).Distance(o.curve.Value(places[k - 1].with));
      return run;
    };
    std::vector<Place>   places;
    std::vector<gp_Pnt>  points;
    std::vector<double>  run;
    Handle(Geom_Curve)   curve;
    Handle(Geom2d_Curve) mine, theirs;
    for (int count = 16; count <= 512; count *= 2)
    {
      std::vector<double> ts(count + 1);
      for (int k = 0; k <= count; ++k)
        ts[k] = piece.from + (piece.to - piece.from) * k / count;
      // Evenly by that length: twice over, from places even by the parameter.
      for (int pass = 0; pass < 3; ++pass)
      {
        places = placesAt(ts);
        run    = lengths(places, points);
        if (pass == 2)
          break;
        std::vector<double> even(count + 1);
        even[0]     = piece.from;
        even[count] = piece.to;
        size_t at   = 1;
        for (int k = 1; k < count; ++k)
        {
          const double want = run.back() * k / count;
          while (at + 1 < run.size() && run[at] < want)
            ++at;
          const double share =
            run[at] > run[at - 1] ? (want - run[at - 1]) / (run[at] - run[at - 1]) : 0;
          even[k] = ts[at - 1] + (ts[at] - ts[at - 1]) * share;
        }
        ts = even;
      }
      Handle(NCollection_HArray1<gp_Pnt>) array = new NCollection_HArray1<gp_Pnt>(1, count + 1);
      Handle(NCollection_HArray1<double>) at3   = new NCollection_HArray1<double>(1, count + 1);
      std::vector<gp_Pnt2d>               inMine, inTheirs;
      for (int k = 0; k <= count; ++k)
      {
        if (k && !(run[k] > run[k - 1]))
          refuse("a crease of the blend has no length between two of its points", e.edge);
        array->SetValue(k + 1, points[k]);
        at3->SetValue(k + 1, run[k]);
        inMine.push_back(uv(e, places[k].t, places[k].d));
        inTheirs.push_back(uv(o, places[k].with, places[k].d));
      }
      GeomAPI_Interpolate interpolate(array, at3, false, 1e-12);
      interpolate.Perform();
      if (!interpolate.IsDone())
        refuse("a crease of the blend could not be put through its points", e.edge);
      curve  = interpolate.Curve();
      mine   = through(inMine,
                     run,
                     piece.with == own && piece.withFrom == piece.from,
                     piece.with == own && piece.withTo == piece.to);
      theirs = through(inTheirs,
                       run,
                       piece.with == own && piece.withFrom == piece.from,
                       piece.with == own && piece.withTo == piece.to);
      // Near enough to both faces between its points?
      double worst = 0;
      for (int k = 0; k < count; ++k)
      {
        const double   x  = (run[k] + run[k + 1]) / 2;
        const gp_Pnt   on = curve->Value(x);
        const gp_Pnt2d a = mine->Value(x), b = theirs->Value(x);
        worst = std::max({worst,
                          on.Distance(e.surface->Value(a.X(), a.Y())),
                          on.Distance(o.surface->Value(b.X(), b.Y()))});
      }
      if (worst < kCrease * scale)
        break;
    }
    Shared shared;
    shared.edge      = joined(curve,
                         corners.all[piece.cornerFrom].v,
                         corners.all[piece.cornerTo].v,
                         0,
                         run.back(),
                         "a crease of the blend could not be made");
    shared.builtBy   = own;
    shared.opposed   = piece.withTo < piece.withFrom;
    creaseIn[0][key] = mine;
    creaseIn[1][key] = theirs;
    return creases.emplace(key, shared).first->second;
  }

  std::map<std::tuple<int, int, int, int>, Handle(Geom2d_Curve)> creaseIn[2];

  //! The faces of the blend over one edge of the outline: one, or where the edge's own bend
  //! closes its face in a crease, one to each side of every such crease. A face cannot have
  //! one edge twice running there and back from a corner that nothing else comes to, so the
  //! face is divided from the end of the crease down to the corner at the edge's start; the
  //! pieces lie in one surface and are smooth across the cut.
  std::vector<TopoDS_Face> facesOf(int own)
  {
    Element&            e     = elements[own];
    const auto&         chain = chains[own];
    std::vector<size_t> cuts; // the first piece after each such crease's end
    for (size_t k = 0; k + 1 < chain.size(); ++k)
      if (chain[k].kind == Piece::Crease && chain[k].with == own
          && chain[k + 1].kind == Piece::Crease && chain[k + 1].with == own
          && chain[k].cornerTo == chain[k + 1].cornerFrom && chain[k].withTo == chain[k].to)
        cuts.push_back(k + 1);
    const int    riseA = riseAt(own, false), riseB = riseAt(own, true);
    const size_t count = cuts.size() + 1;
    const bool   seam  = count == 1 && riseA == riseB && !rises[riseA].edge.IsNull();

    struct Cut
    {
      TopoDS_Edge          slant, rise;
      Handle(Geom2d_Curve) slantIn, riseIn;
    };

    std::vector<Cut> made;
    if (count > 1)
    {
      const auto spline = Handle(Geom_BSplineSurface)::DownCast(e.surface);
      if (spline.IsNull())
        refuse("a straight edge of the outline is closed by its own bend", e.edge);
      // The cuts run from the start of the edge low along it, each lower than the one
      // before, and then straight up to their crease: under everything the face has.
      double least = HUGE_VAL;
      for (int k = 1; k < kScan; ++k)
      {
        const double t = e.a + (e.b - e.a) * k / kScan;
        least          = std::min(least, across(reach(own, t).d) / (t - e.a));
      }
      for (const size_t cut : cuts)
        least = std::min(least, across(chain[cut].dFrom) / (chain[cut].from - e.a));
      TopoDS_Vertex first, second;
      TopExp::Vertices(e.edge, first, second);
      for (size_t k = 0; k < cuts.size(); ++k)
      {
        const Piece& after = chain[cuts[k]];
        const double t = after.from, low = least / (2.0 * double(k + 1)) * (t - e.a),
                     high = across(after.dFrom);
        const gp_Pnt2d      from(e.a, 0), turn(t, low);
        const TopoDS_Vertex corner = BRepBuilderAPI_MakeVertex(spline->Value(t, low)).Vertex();
        Cut                 cut;
        const gp_Vec2d      along(from, turn);
        cut.slantIn = new Geom2d_Line(from, gp_Dir2d(along));
        // Its curve in space from its curve in the face, to a tolerance.
        BRepBuilderAPI_MakeEdge slant(cut.slantIn, spline, first, corner, 0, along.Magnitude());
        if (!slant.IsDone()
            || !BRepLib::BuildCurve3d(slant.Edge(),
                                      BRepExt::tolerance::exactLength,
                                      GeomAbs_C1,
                                      14,
                                      256))
          refuse("a face of the blend could not be divided", e.edge);
        cut.slant  = slant.Edge();
        cut.rise   = joined(spline->UIso(t),
                          corner,
                          corners.all[after.cornerFrom].v,
                          low,
                          high,
                          "a face of the blend could not be divided");
        cut.riseIn = new Geom2d_Line(gp_Pnt2d(t, 0), gp_Dir2d(0, 1));
        made.push_back(cut);
      }
    }
    std::vector<TopoDS_Face> faces;
    for (size_t part = 0; part < count; ++part)
    {
      // Each piece has the surface for itself, so that an edge between two of them has
      // its curve in each.
      const Handle(Geom_Surface) surface =
        count == 1 ? e.surface : Handle(Geom_Surface)::DownCast(e.surface->Copy());
      TopoDS_Face face;
      builder.MakeFace(face, surface, BRepExt::tolerance::length);
      TopoDS_Wire wire;
      builder.MakeWire(wire);
      const auto put =
        [&](const TopoDS_Edge& edge, const Handle(Geom2d_Curve)& onFace, bool forward) {
          builder.UpdateEdge(edge, onFace, face, BRepExt::tolerance::length);
          fitTolerance(edge, onFace, surface);
          builder.Add(wire, edge.Oriented(forward ? TopAbs_FORWARD : TopAbs_REVERSED));
        };
      const size_t from = part == 0 ? 0 : cuts[part - 1],
                   to   = part + 1 == count ? chain.size() : cuts[part];
      if (part + 1 == count)
      {
        // The foot edge, as its curve runs, and up at its end.
        Handle(Geom2d_Curve) base;
        if (!e.straight || bevel)
          base = new Geom2d_Line(gp_Pnt2d(0, 0), gp_Dir2d(1, 0));
        else
          base = new Geom2d_Line(gp_Pnt2d(e.flipped ? M_PI / 2 : 0, 0), gp_Dir2d(0, 1));
        put(e.edge, base, true);
        if (seam)
        {
          // One edge closed on itself: the section at its one corner is the face's seam.
          const auto atB = riseIn(e, e.b, rises[riseA]), atA = riseIn(e, e.a, rises[riseA]);
          builder.UpdateEdge(rises[riseA].edge, atB, atA, face, BRepExt::tolerance::length);
          fitTolerance(rises[riseA].edge, atB, surface);
          fitTolerance(rises[riseA].edge, atA, surface);
          builder.Add(wire, rises[riseA].edge.Oriented(TopAbs_FORWARD));
        }
        else if (!rises[riseB].edge.IsNull())
          put(rises[riseB].edge, riseIn(e, e.b, rises[riseB]), true);
      }
      else
      {
        put(made[part].slant, made[part].slantIn, true);
        put(made[part].rise, made[part].riseIn, true);
      }
      for (size_t index = to; index-- > from;)
      {
        const Piece* piece = &chain[index];
        if (piece->kind == Piece::Cap)
        {
          const TopoDS_Edge capEdge =
            joined(e.movedInSpace,
                   corners.all[piece->cornerFrom].v,
                   corners.all[piece->cornerTo].v,
                   piece->from,
                   piece->to,
                   "an edge of what is left of the face could not be made");
          Handle(Geom2d_Curve) onFace;
          if (!e.straight)
            onFace = new Geom2d_Line(gp_Pnt2d(0, 1), gp_Dir2d(1, 0));
          else if (bevel)
            onFace = new Geom2d_Line(gp_Pnt2d(0, r * M_SQRT2), gp_Dir2d(1, 0));
          else
            onFace = new Geom2d_Line(gp_Pnt2d(e.flipped ? 0 : M_PI / 2, 0), gp_Dir2d(0, 1));
          put(capEdge, onFace, false);
          CapEdge cap;
          cap.edge    = capEdge;
          cap.from    = e.s > 0 ? piece->cornerFrom : piece->cornerTo;
          cap.to      = e.s > 0 ? piece->cornerTo : piece->cornerFrom;
          cap.inPlane = e.moved;
          for (int k = 0; k <= 16; ++k)
          {
            const double t = e.s > 0 ? piece->from + (piece->to - piece->from) * k / 16
                                     : piece->to + (piece->from - piece->to) * k / 16;
            cap.polygon.push_back(e.moved->Value(t));
          }
          if (e.s < 0)
            cap.edge.Reverse();
          capEdges.push_back(cap);
        }
        else if (piece->kind == Piece::Apex)
        {
          // All the sections of an arc end in one point.
          const TopoDS_Edge point =
            degenerateEdge(new Geom2d_Line(gp_Pnt2d(0, across(piece->dFrom)), gp_Dir2d(1, 0)),
                           surface,
                           TopLoc_Location(),
                           corners.all[piece->cornerFrom].v,
                           piece->from,
                           piece->to,
                           BRepExt::tolerance::length);
          builder.Add(wire, point.Oriented(TopAbs_REVERSED));
        }
        else
        {
          Shared&    shared = creaseOf(own, *piece);
          const auto key    = std::make_tuple(std::min(own, piece->with),
                                           std::max(own, piece->with),
                                           std::min(piece->cornerFrom, piece->cornerTo),
                                           std::max(piece->cornerFrom, piece->cornerTo));
          // Each face goes back along its top against its own edge. Of a crease
          // between two faces of one edge, the first to come here made it.
          const bool mine = piece->with == own ? closing.insert(key).second : shared.builtBy == own;
          put(shared.edge, creaseIn[mine ? 0 : 1][key], !mine && shared.opposed);
        }
      }
      if (part == 0)
      {
        if (seam)
          builder.Add(wire, rises[riseA].edge.Oriented(TopAbs_REVERSED));
        else if (!rises[riseA].edge.IsNull())
          put(rises[riseA].edge, riseIn(e, e.a, rises[riseA]), false);
      }
      else
      {
        put(made[part - 1].rise, made[part - 1].riseIn, false);
        put(made[part - 1].slant, made[part - 1].slantIn, false);
      }
      wire.Closed(true);
      builder.Add(face, wire);
      // The wire goes along the edge, up, back and down: the face is to its left where the
      // surface's first parameter is along the edge and its second across. A cylinder has
      // them the other way about, and is to the right unless its angle runs backwards.
      if (e.straight && !bevel && !e.flipped)
      {
        TopoDS_Face other;
        builder.MakeFace(other, surface, BRepExt::tolerance::length);
        builder.Add(other, wire.Reversed());
        face = other;
      }
      faces.push_back(face);
    }
    return faces;
  }

  std::set<std::tuple<int, int, int, int>> closing;

  //! What is left of the face: flat faces inside the pieces of the outline moved in.
  std::vector<TopoDS_Face> caps(const gp_Vec& outward)
  {
    struct Loop
    {
      TopoDS_Wire           wire;
      std::vector<gp_Pnt2d> polygon;
      double                area = 0;
      std::vector<int>      edges;
    };

    std::vector<Loop> loops;
    std::vector<bool> used(capEdges.size(), false);
    for (size_t start = 0; start < capEdges.size(); ++start)
    {
      if (used[start])
        continue;
      Loop loop;
      builder.MakeWire(loop.wire);
      size_t now = start;
      for (size_t guard = 0; guard <= capEdges.size(); ++guard)
      {
        used[now] = true;
        loop.edges.push_back(int(now));
        builder.Add(loop.wire, capEdges[now].edge);
        for (const auto& p : capEdges[now].polygon)
          loop.polygon.push_back(p);
        if (capEdges[now].to == capEdges[start].from)
          break;
        size_t next = capEdges.size();
        for (size_t k = 0; k < capEdges.size(); ++k)
          if (!used[k] && capEdges[k].from == capEdges[now].to)
            next = k;
        if (next == capEdges.size())
          refuse("what is left of the face does not close");
        now = next;
      }
      loop.wire.Closed(true);
      for (size_t k = 0; k < loop.polygon.size(); ++k)
      {
        const auto& p = loop.polygon[k];
        const auto& q = loop.polygon[(k + 1) % loop.polygon.size()];
        loop.area += (p.X() * q.Y() - q.X() * p.Y()) / 2;
      }
      loops.push_back(loop);
    }
    const auto inside = [](const std::vector<gp_Pnt2d>& polygon, const gp_Pnt2d& p) {
      bool in = false;
      for (size_t k = 0, j = polygon.size() - 1; k < polygon.size(); j = k++)
      {
        const auto& a = polygon[k];
        const auto& b = polygon[j];
        if ((a.Y() > p.Y()) != (b.Y() > p.Y())
            && p.X() < (b.X() - a.X()) * (p.Y() - a.Y()) / (b.Y() - a.Y()) + a.X())
          in = !in;
      }
      return in;
    };
    gp_Ax3 position = foot.Position();
    position.SetLocation(foot.Location().Translated(w * -r));
    const Handle(Geom_Plane) plane = new Geom_Plane(position);
    gp_Pnt                   origin;
    gp_Vec                   du, dv;
    plane->D1(0, 0, origin, du, dv);
    const bool               reversed = du.Crossed(dv).Dot(outward) < 0;
    std::vector<TopoDS_Face> faces;
    for (const auto& outer : loops)
    {
      if (outer.area <= 0)
        continue;
      TopoDS_Face face;
      builder.MakeFace(face, plane, BRepExt::tolerance::length);
      const auto take = [&](const Loop& loop) {
        for (int k : loop.edges)
        {
          builder.UpdateEdge(capEdges[k].edge,
                             capEdges[k].inPlane,
                             face,
                             BRepExt::tolerance::length);
          fitTolerance(capEdges[k].edge, capEdges[k].inPlane, plane);
        }
        builder.Add(face, loop.wire);
      };
      take(outer);
      for (const auto& hole : loops)
        if (hole.area < 0 && inside(outer.polygon, hole.polygon.front()))
          take(hole);
      if (reversed)
        face.Reverse();
      faces.push_back(face);
    }
    return faces;
  }
};

//! A closing concave round is a machining envelope, not the distance profile. The retained
//! cap pieces are the trimmed inward offset: their projection bounds the valid ball centres.
//! For a convex polygon of centres the envelope consists exactly of cylinders on its sides
//! and spheres at its vertices. Above the centres it continues as the vertical tool. Without
//! that continuation the unreachable corner would end in a sharp horizontal shelf at r.
BooleanResult closingBall(Maker&              maker,
                          const TopoDS_Shape& shape,
                          const TopoDS_Face&  floor,
                          const gp_Vec&       up,
                          const Adjacency&    adjacent)
{
  const char* op      = maker.op;
  const auto  tangent = [&](const TopoDS_Edge& edge, const NCollection_List<TopoDS_Shape>& faces) {
    if (faces.Extent() != 2)
      return false;
    BRepAdaptor_Curve curve(edge);
    const double      cosine = std::cos(0.5 * M_PI / 180);
    const BoundaryJet aAt(op, edge, TopoDS::Face(faces.First())),
      bAt(op, edge, TopoDS::Face(faces.Last()));
    for (int i = 0; i < 257; ++i)
    {
      const double t =
        curve.FirstParameter() + (curve.LastParameter() - curve.FirstParameter()) * (i + 0.5) / 257;
      const auto a = aAt.at(t);
      const auto b = bAt.at(t);
      if (a.n.Dot(b.n) < cosine)
        return false;
    }
    return true;
  };
  if (maker.bevel)
    fail(op, unsupported_geometry, "a closing concave bevel has no rolling-ball envelope: refused");

  // This builder takes a single convex polygon of centres. Curved retained offsets, holes
  // and disconnected regions need other envelopes; do not substitute a crease for them.
  for (int i = 0; i < int(maker.elements.size()); ++i)
  {
    auto& e = maker.elements[i];
    for (const auto& piece : maker.chains[i])
    {
      if (piece.kind != Piece::Cap)
        continue;
      if (e.kind != GeomAbs_Line)
        fail(op,
             unsupported_geometry,
             "closing concave round requires straight trimmed-offset centre paths: refused");
      maker.surfaceOf(e);
      Maker::CapEdge cap;
      cap.from    = e.s > 0 ? piece.cornerFrom : piece.cornerTo;
      cap.to      = e.s > 0 ? piece.cornerTo : piece.cornerFrom;
      cap.inPlane = e.moved;
      cap.edge    = maker.joined(e.movedInSpace,
                              maker.corners.all[piece.cornerFrom].v,
                              maker.corners.all[piece.cornerTo].v,
                              piece.from,
                              piece.to,
                              "a trimmed-offset centre path could not be built");
      if (e.s < 0)
        cap.edge.Reverse();
      cap.polygon = {e.moved->Value(e.s > 0 ? piece.from : piece.to),
                     e.moved->Value(e.s > 0 ? piece.to : piece.from)};
      maker.capEdges.push_back(cap);
    }
  }
  const auto caps = maker.caps(up);
  if (caps.size() != 1)
    fail(op,
         unsupported_geometry,
         "closing concave round requires one convex trimmed-offset region: refused");
  int wires = 0;
  for (TopExp_Explorer it(caps.front(), TopAbs_WIRE); it.More(); it.Next())
    ++wires;
  if (wires != 1)
    fail(op,
         unsupported_geometry,
         "closing concave round with centre-path holes is not built: refused");
  TopoDS_Wire         centreWire = BRepTools::OuterWire(caps.front());
  std::vector<gp_Pnt> polygon;
  double              perimeter = 0;
  for (BRepTools_WireExplorer it(centreWire, caps.front()); it.More(); it.Next())
  {
    polygon.push_back(BRep_Tool::Pnt(it.CurrentVertex()));
    BRepAdaptor_Curve curve(it.Current());
    perimeter += curve.Value(curve.FirstParameter()).Distance(curve.Value(curve.LastParameter()));
  }
  if (polygon.size() < 3)
    fail(op, no_solution, "the trimmed-offset polygon has fewer than three sides");
  const double tol = 1e-7 * maker.scale;
  for (size_t i = 0; i < polygon.size(); ++i)
  {
    const gp_Vec edge(polygon[i], polygon[(i + 1) % polygon.size()]);
    for (const auto& p : polygon)
      if (edge.Crossed(gp_Vec(polygon[i], p)).Dot(up) < -tol * edge.Magnitude())
        fail(op,
             unsupported_geometry,
             "closing concave round with a nonconvex centre region is not built: refused");
  }

  // The vertical continuation must stop at a common planar mouth. Certify the complete
  // walls, not just their first r: a tapered wall or an overhang must not be filled by this
  // prismatic tool. Every boundary is a level curve or a straight vertical generator.
  const gp_Pnt origin = BRepAdaptor_Surface(floor).Plane().Location();
  const auto   height = [&](const gp_Pnt& p) { return gp_Vec(origin, p).Dot(up); };
  double       depth  = 0;
  ShapeSet     walls;
  for (TopExp_Explorer it(floor, TopAbs_EDGE); it.More(); it.Next())
    for (const auto& face : adjacent.FindFromKey(it.Current()))
      if (!face.IsSame(floor))
        walls.Add(face);
  for (const auto& wall : walls)
  {
    BRepAdaptor_Surface surface(TopoDS::Face(wall));
    const bool          vertical =
      surface.GetType() == GeomAbs_Plane
                 ? std::abs(planeNormal(surface.Plane()).Dot(up)) < BRepExt::tolerance::angle
                 : surface.GetType() == GeomAbs_Cylinder
            && surface.Cylinder().Axis().Direction().IsParallel(gp_Dir(up),
                                                                BRepExt::tolerance::angle);
    if (!vertical)
      fail(op,
           unsupported_geometry,
           "closing concave round requires planar or cylindrical vertical prism walls: refused");
    double high = 0;
    for (TopExp_Explorer v(wall, TopAbs_VERTEX); v.More(); v.Next())
      high = std::max(high, height(BRep_Tool::Pnt(TopoDS::Vertex(v.Current()))));
    if (!(high > maker.r + tol) || (depth > 0 && std::abs(depth - high) > tol))
      fail(op,
           unsupported_geometry,
           "closing concave round requires walls ending at one planar mouth above the ball "
           "centres: refused");
    depth = high;
    for (TopExp_Explorer edge(wall, TopAbs_EDGE); edge.More(); edge.Next())
    {
      BRepAdaptor_Curve curve(TopoDS::Edge(edge.Current()));
      const double      a = curve.FirstParameter(), b = curve.LastParameter();
      const auto        p = curve.Value(a), q = curve.Value(b);
      const double      ha = height(p), hb = height(q);
      const bool        level = (std::abs(ha) < tol && std::abs(hb) < tol)
                         || (std::abs(ha - depth) < tol && std::abs(hb - depth) < tol);
      if (level)
      {
        if (curve.GetType() != GeomAbs_Line
            && !(curve.GetType() == GeomAbs_Circle
                 && curve.Circle().Axis().Direction().IsParallel(gp_Dir(up),
                                                                 BRepExt::tolerance::angle)))
          fail(op,
               unsupported_geometry,
               "closing concave round requires straight or circular prism-level wall boundaries: "
               "refused");
        if (std::abs(ha - depth) < tol)
        {
          bool mouth = false;
          for (const auto& neighbour : adjacent.FindFromKey(edge.Current()))
          {
            if (neighbour.IsSame(wall))
              continue;
            BRepAdaptor_Surface cap(TopoDS::Face(neighbour));
            if (cap.GetType() != GeomAbs_Plane)
              continue;
            gp_Vec normal = planeNormal(cap.Plane());
            if (neighbour.Orientation() == TopAbs_REVERSED)
              normal.Reverse();
            mouth |=
              normal.Dot(up) > 0 && normal.Crossed(up).Magnitude() < BRepExt::tolerance::angle;
          }
          if (!mouth)
            fail(op,
                 unsupported_geometry,
                 "closing concave round requires an open planar mouth; a ceiling or overhang is "
                 "not built: refused");
        }
      }
      else if (curve.GetType() != GeomAbs_Line || gp_Vec(p, q).Crossed(up).Magnitude() > tol
               || std::abs(std::min(ha, hb)) > tol || std::abs(std::max(ha, hb) - depth) > tol)
        fail(op,
             unsupported_geometry,
             "closing concave round has a wall that is not a full vertical prism: refused");
    }
  }

  // The whole centre region, at the ball's height, must be at least r from retained
  // material. This also catches a remote overhang above an otherwise open prism mouth.
  gp_Trsf toCentres;
  toCentres.SetTranslation(up * maker.r);
  const TopoDS_Shape centres = BRepBuilderAPI_Transform(caps.front(), toCentres, true).Shape();
  BRepExtrema_DistShapeShape clearance(shape, centres);
  if (!clearance.IsDone())
    fail(op,
         no_solution,
         "closing concave round could not certify the ball-centre region's clearance");
  if (clearance.Value() < maker.r - tol)
    fail(op,
         unsupported_geometry,
         "closing concave round's resting ball intersects retained material: refused");

  // Dilate the trimmed offset, then round the bottom of that tool. At each polygon vertex
  // the dilating arc has radius r; rounding it by r makes a sphere, not a spindle torus.
  BRepOffsetAPI_MakeOffset offset(caps.front(), GeomAbs_Arc, false);
  offset.Perform(maker.r);
  if (!offset.IsDone())
    fail(op, no_solution, "the trimmed-offset tool could not be dilated");
  TopoDS_Wire outline;
  for (TopExp_Explorer it(offset.Shape(), TopAbs_WIRE); it.More(); it.Next())
  {
    if (!outline.IsNull())
      fail(op, no_solution, "the dilated centre polygon split into several outlines");
    outline = TopoDS::Wire(it.Current());
  }
  if (outline.IsNull())
    fail(op, no_solution, "the dilated centre polygon has no outline");
  BRepBuilderAPI_MakeFace footprint(BRepAdaptor_Surface(floor).Plane(), outline, true);
  if (!footprint.IsDone())
    fail(op, no_solution, "the ball tool's footprint could not be closed");
  const double             toolHeight = depth + maker.r;
  const TopoDS_Shape       prism = BRepPrimAPI_MakePrism(footprint.Face(), up * toolHeight).Shape();
  BRepFilletAPI_MakeFillet round(prism);
  for (TopExp_Explorer it(prism, TopAbs_EDGE); it.More(); it.Next())
  {
    BRepAdaptor_Curve curve(TopoDS::Edge(it.Current()));
    if (std::abs(height(curve.Value(curve.FirstParameter()))) < tol
        && std::abs(height(curve.Value(curve.LastParameter()))) < tol)
      round.Add(maker.r, TopoDS::Edge(it.Current()));
  }
  round.Build();
  if (!round.IsDone() || !BRepCheck_Analyzer(round.Shape()).IsValid())
    fail(op, no_solution, "the trimmed-offset ball envelope could not be built as a valid tool");
  const TopoDS_Shape tool = round.Shape();
  const double       r = maker.r, area = faceMeasure(caps.front(), true, gp_Pnt(), nullptr);
  const double expected = area * toolHeight + perimeter * (M_PI * r * r / 4 + (toolHeight - r) * r)
                          + M_PI * (2 * r * r * r / 3 + (toolHeight - r) * r * r);
  if (std::abs(held(tool) - expected) > 1e-7 * expected)
    fail(op,
         no_solution,
         "the ball envelope does not hold the convex centre polygon's swept volume");
  // No approximate corner patch is admitted here. All curved surfaces must be the exact
  // spheres/cylinders of radius r, and their shared edges below the mouth must be tangent.
  for (TopExp_Explorer it(tool, TopAbs_FACE); it.More(); it.Next())
  {
    BRepAdaptor_Surface surface(TopoDS::Face(it.Current()));
    if (surface.GetType() == GeomAbs_Plane)
      continue;
    const double radius = surface.GetType() == GeomAbs_Sphere     ? surface.Sphere().Radius()
                          : surface.GetType() == GeomAbs_Cylinder ? surface.Cylinder().Radius()
                                                                  : -1;
    if (std::abs(radius - r) > tol)
      fail(op,
           no_solution,
           "the ball envelope contains a surface other than a sphere or cylinder of the requested "
           "radius");
    bool onPath = false;
    if (surface.GetType() == GeomAbs_Sphere)
    {
      const auto centre = surface.Sphere().Location();
      for (const auto& vertex : polygon)
        onPath |= centre.Distance(vertex.Translated(up * r)) < tol;
    }
    else
    {
      const gp_Lin axis(surface.Cylinder().Axis());
      if (axis.Direction().IsParallel(gp_Dir(up), BRepExt::tolerance::angle))
      {
        for (const auto& vertex : polygon)
          onPath |= axis.Distance(vertex) < tol;
      }
      else if (std::abs(height(axis.Location()) - r) < tol)
      {
        for (size_t i = 0; i < polygon.size(); ++i)
          onPath |= axis.Distance(polygon[i].Translated(up * r)) < tol
                    && axis.Distance(polygon[(i + 1) % polygon.size()].Translated(up * r)) < tol;
      }
    }
    if (!onPath)
      fail(op, no_solution, "the ball envelope's centres are not on the trimmed-offset path");
  }
  Adjacency toolFaces;
  TopExp::MapShapesAndAncestors(tool, TopAbs_EDGE, TopAbs_FACE, toolFaces);
  for (int i = 1; i <= toolFaces.Extent(); ++i)
  {
    const auto edge = TopoDS::Edge(toolFaces.FindKey(i));
    if (BRep_Tool::Degenerated(edge))
      continue;
    BRepAdaptor_Curve curve(edge);
    if (height(curve.Value((curve.FirstParameter() + curve.LastParameter()) / 2)) >= depth - tol)
      continue;
    const auto& faces = toolFaces.FindFromIndex(i);
    if (!tangent(edge, faces))
      fail(op, no_solution, "the trimmed-offset ball envelope has a nontangent shared edge");
  }
  BooleanRequest request;
  request.operation = op;
  request.select    = Select::Difference;
  request.sources   = {BRepPrimAPI_MakePrism(floor, up * depth).Shape()};
  request.tools     = {tool};
  const auto fill   = booleanCells(request);
  request.select    = Select::Union;
  request.sources   = {shape};
  request.tools     = {fill.shape};
  auto result       = booleanCells(request);
  if (singleSolid(result.shape).IsNull() || !BRepCheck_Analyzer(result.shape).IsValid()
      || held(result.shape) <= held(shape))
    fail(op, no_solution, "the closing concave round did not give one valid filled body");
  clearance.LoadS1(result.shape);
  clearance.Perform();
  if (!clearance.IsDone() || clearance.Value() < r - tol)
    fail(op, no_solution, "the closing concave round's result does not clear the resting ball");
  std::vector<TopoDS_Shape> selected;
  for (TopExp_Explorer e(floor, TopAbs_EDGE); e.More(); e.Next())
    selected.push_back(e.Current());
  auto      roles = blendBoundaryRoles(shape, selected, result.shape, result);
  Adjacency boundaries;
  TopExp::MapShapesAndUniqueAncestors(result.shape, TopAbs_EDGE, TopAbs_FACE, boundaries);
  for (int i = 1; i <= boundaries.Extent(); ++i)
  {
    const auto edge = TopoDS::Edge(boundaries.FindKey(i));
    if (BRep_Tool::Degenerated(edge))
      continue;
    BRepAdaptor_Curve curve(edge);
    if (height(curve.Value((curve.FirstParameter() + curve.LastParameter()) / 2)) >= depth - tol)
      roles.sharp.Add(edge); // the construction's mouth is a permitted sharp stop
  }
  certifyBlendBoundaries(op, result.shape, roles, maker.r, .5, 257);
  return result;
}

} // namespace

BooleanResult faceBlend(const char*                      op,
                        const TopoDS_Shape&              shape,
                        const std::vector<TopoDS_Shape>& edges,
                        double                           distance,
                        bool                             bevel)
{
  if (!(distance > 0) || !std::isfinite(distance) || edges.empty())
    fail(op, invalid_input, "edges and a positive distance are required");
  ShapeSet asked;
  for (const auto& edge : edges)
    asked.Add(edge);
  // The flat face whose outline these edges are, all of them.
  TopoDS_Face top;
  for (TopExp_Explorer it(shape, TopAbs_FACE); it.More() && top.IsNull(); it.Next())
  {
    const auto face = TopoDS::Face(it.Current());
    if (BRepAdaptor_Surface(face).GetType() != GeomAbs_Plane)
      continue;
    ShapeSet own;
    TopExp::MapShapes(face, TopAbs_EDGE, own);
    bool all = own.Extent() == asked.Extent();
    for (int k = 1; k <= own.Extent() && all; ++k)
      all = asked.Contains(own(k));
    if (all)
      top = face;
  }
  if (top.IsNull())
    fail(op, unsupported_geometry, "the edges are not the whole outline of one flat face");
  const gp_Pln plane = BRepAdaptor_Surface(top).Plane();
  // The surface's own normal (a plane's axis is the other way in a mirrored frame).
  gp_Vec outward = planeNormal(plane);
  if (top.Orientation() == TopAbs_REVERSED)
    outward.Reverse();
  // Which way the walls go from the face: into the body at a convex outline.
  Adjacency facesOfEdge;
  TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, facesOfEdge);
  double side = 0;
  for (int k = 1; k <= asked.Extent(); ++k)
  {
    const auto edge = TopoDS::Edge(asked(k));
    if (BRep_Tool::Degenerated(edge))
      continue;
    TopoDS_Face wall;
    for (const auto& face : facesOfEdge.FindFromKey(edge))
      if (!face.IsSame(top))
        wall = TopoDS::Face(face);
    if (wall.IsNull())
      fail(op, unsupported_geometry, "an edge of the outline has no wall", {edge});
    BRepAdaptor_Curve curve(edge);
    for (double at : {0.2, 0.5, 0.8})
    {
      const double t =
        curve.FirstParameter() + (curve.LastParameter() - curve.FirstParameter()) * at;
      gp_Pnt p;
      gp_Vec tangent, into;
      curve.D1(t, p, tangent);
      if (!supportInward(wall, edge, t, tangent, into)
          || into.Crossed(outward).Magnitude() > BRepExt::tolerance::angle * into.Magnitude())
        fail(op, unsupported_geometry, "a wall does not stand square on the face", {edge});
      const double way = into.Dot(outward) > 0 ? 1 : -1;
      if (side != 0 && way != side)
        fail(op,
             unsupported_geometry,
             "walls go both ways from the face: its outline is convex in part and concave in part",
             {edge});
      side = way;
    }
  }
  Maker maker;
  maker.op    = op;
  maker.r     = distance;
  maker.bevel = bevel;
  maker.w     = outward * side;
  Bnd_Box box;
  BRepBndLib::Add(shape, box, false);
  const double diagonal = std::sqrt(box.SquareExtent());
  maker.scale           = std::max(1.0, distance);
  // Two corners of the blend nearer than this (a share of the blend's size) are one.
  maker.corners.near = BRepExt::tolerance::sewing * std::max(1.0, distance);

  // Step 1: the body cut off or filled up at the distance.
  const double   area = faceMeasure(top, true, gp_Pnt(), nullptr), before = held(shape);
  const gp_Pnt   footAt = plane.Location().Translated(maker.w * distance);
  BooleanRequest request;
  request.operation = op;
  request.sources   = {shape};
  if (side < 0)
  {
    const gp_Vec x(plane.XAxis().Direction()), y = outward.Crossed(x);
    const double reach =
      2 * diagonal
      + gp_Vec(gp_Pnt((box.CornerMin().XYZ() + box.CornerMax().XYZ()) / 2), footAt).Magnitude();
    request.select = Select::Difference;
    request.tools  = {
      BRepPrimAPI_MakeBox(
        gp_Ax2(footAt.Translated(x * -reach + y * -reach), gp_Dir(outward), gp_Dir(x)),
        2 * reach,
        2 * reach,
        reach)
        .Shape()};
  }
  else
  {
    request.select = Select::Union;
    request.tools  = {BRepPrimAPI_MakePrism(top, outward * distance).Shape()};
  }
  auto               cutOff = booleanCells(request);
  const TopoDS_Shape tub    = cutOff.shape;
  const double       inTub  = held(tub);
  // Between the foot and the face the body has to be the prism over the face and no more.
  if (tub.IsNull() || !BRepCheck_Analyzer(tub).IsValid()
      || std::abs(std::abs(before - inTub) - area * distance)
           > 1e-6 * std::max(1.0, before) + 1e-4 * area * distance)
    fail(op,
         unsupported_geometry,
         "the body is not a prism over the face for as far as the blend goes down its walls");
  TopoDS_Face foot;
  for (TopExp_Explorer it(tub, TopAbs_FACE); it.More(); it.Next())
  {
    const auto          face = TopoDS::Face(it.Current());
    BRepAdaptor_Surface surface(face);
    if (surface.GetType() != GeomAbs_Plane)
      continue;
    const gp_Pln own = surface.Plane();
    if (own.Distance(footAt) > BRepExt::tolerance::length
        || !own.Axis().Direction().IsParallel(gp_Dir(outward), BRepExt::tolerance::angle))
      continue;
    gp_Vec normal = planeNormal(own);
    if (face.Orientation() == TopAbs_REVERSED)
      normal.Reverse();
    if (normal.Dot(outward) < 0)
      continue;
    if (!foot.IsNull())
      fail(op, unsupported_geometry, "the face is in several pieces at the foot of the blend");
    foot = face;
  }
  if (foot.IsNull())
    fail(op, no_solution, "the foot of the blend was not found");
  if (std::abs(faceMeasure(foot, true, gp_Pnt(), nullptr) - area) > 1e-6 * std::max(1.0, area))
    fail(op,
         unsupported_geometry,
         "the walls do not stand square on the face for as far as the blend goes down them");
  {
    const Handle(Geom_Plane) own = Handle(Geom_Plane)::DownCast(BRep_Tool::Surface(foot));
    if (own.IsNull())
      fail(op, no_solution, "the foot of the blend is not a plane");
    maker.foot = own->Pln();
    gp_Pnt origin;
    own->D1(0, 0, origin, maker.X, maker.Y);
  }

  // Steps 2 to 4.
  maker.addElements(foot);
  maker.chains.resize(maker.elements.size());
  for (auto& e : maker.elements)
  {
    TopoDS_Vertex first, second;
    TopExp::Vertices(e.edge, first, second);
    maker.corners.add(first);
    maker.corners.add(second);
  }
  for (int k = 0; k < int(maker.elements.size()); ++k)
    maker.chainOf(k);
  bool left = false;
  for (const auto& chain : maker.chains)
    for (const auto& piece : chain)
      left |= piece.kind == Piece::Cap;
  if (!left)
    fail(op,
         no_solution,
         "requested distance consumes the planar cap: top face vanishes",
         {},
         FailureKind::CapConsumed);
  // A concave closing outline uses the envelope of the ball on the trimmed offset. Convex
  // closing outlines keep the distance profile and its creases.
  if (side > 0)
    for (const auto& chain : maker.chains)
      for (const auto& piece : chain)
        if (piece.kind != Piece::Cap)
          try
          {
            return closingBall(maker, shape, top, outward, facesOfEdge);
          }
          catch (Failure& failure)
          {
            if (failure.kind == FailureKind::Construction)
              failure.kind = FailureKind::ClosingConcave;
            throw;
          }
  for (auto& e : maker.elements)
    maker.surfaceOf(e);
  std::vector<TopoDS_Face> lid;
  for (int k = 0; k < int(maker.elements.size()); ++k)
  {
    const Element& e = maker.elements[k];
    // Outwards is towards the corner that the blend takes off, or away from the one it
    // fills: the same for every face over one edge.
    const double   t   = (e.a + e.b) / 2;
    const Hit      hit = maker.reach(k, t);
    const gp_Pnt2d uv  = maker.uv(e, t, hit.d / 2);
    gp_Pnt         p;
    gp_Vec         du, dv;
    e.surface->D1(uv.X(), uv.Y(), p, du, dv);
    const gp_Vec corner(p, maker.lift(e.frame(t).p, 0));
    const bool   reversed = du.Crossed(dv).Dot(corner) * (side < 0 ? 1 : -1) < 0;
    for (auto face : maker.facesOf(k))
    {
      if (reversed)
        face.Reverse();
      lid.push_back(face);
    }
  }
  const auto caps = maker.caps(outward);
  if (caps.empty())
    fail(op,
         no_solution,
         "requested distance consumes the planar cap: top face vanishes",
         {},
         FailureKind::CapConsumed);
  for (const auto& corner : maker.corners.all)
    maker.builder.UpdateVertex(corner.v, BRepExt::tolerance::length);
  for (const auto& face : lid)
    for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next())
      maker.fitVertices(TopoDS::Edge(it.Current()));

  // Step 5: the foot out, the blend in.
  BRep_Builder builder;
  TopoDS_Shell shell;
  builder.MakeShell(shell);
  for (TopExp_Explorer it(tub, TopAbs_FACE); it.More(); it.Next())
    if (!it.Current().IsSame(foot))
      builder.Add(shell, it.Current());
  for (const auto& face : lid)
    builder.Add(shell, face);
  for (const auto& face : caps)
    builder.Add(shell, face);
  shell.Closed(true);
  BRepBuilderAPI_MakeSolid solid(shell);
  if (!solid.IsDone())
    fail(op, no_solution, "the blend did not close the body");
  const TopoDS_Shape made = solid.Solid();
  if (const auto what = invalidity(made); !what.empty())
    fail(op,
         no_solution,
         "the blend of the face's outline did not give a valid body (shape kind and status:" + what
           + ")");
  const double after = held(made, true);
  const double low = std::min(before, inTub), high = std::max(before, inTub);
  if (!(after > low) || !(after < high))
    fail(op, no_solution, "the blend of the face's outline does not hold what it should");
  BooleanResult result;
  result.shape = made;
  result.faces = cutOff.faces;
  NCollection_List<TopoDS_Shape> remains;
  for (const auto& face : caps)
    remains.Append(face);
  if (result.faces.IsBound(top))
    result.faces.UnBind(top);
  result.faces.Bind(top, remains);
  result.rungs = cutOff.rungs;
  if (!bevel)
  {
    std::vector<TopoDS_Shape> generated(lid.begin(), lid.end());
    auto                      roles = blendBoundaryRoles(shape, edges, made, result, generated);
    for (const auto& crease : maker.creases)
      roles.sharp.Add(crease.second.edge);
    certifyBlendBoundaries(op, made, roles, distance);
  }
  return result;
}

} // namespace BRepExt
