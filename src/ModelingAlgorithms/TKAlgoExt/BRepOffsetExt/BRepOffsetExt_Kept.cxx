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

// The offset that keeps the topology: each face on its own offset surface, each edge and vertex
// where those surfaces meet. The reference modeller's face offset keeps it (the same number of
// faces before and after in 536 of its 537 results), which is the intersection join.
//
// Faces. BRepOffset::Surface gives the offset of a surface in the surface's own parameters: a plane
// is translated, a cylinder, sphere, torus or cone keeps its axes and (u, v), anything else is a
// Geom_OffsetSurface of the same (u, v). So a point (u, v) of a face goes to (u, v) of its offset.
//
// Edges. Where two faces are tangent along an edge and move by the same distance, they meet along
// the offset of the edge: at a point c of it both have the normal n, so c + d n is on both offset
// surfaces, with the normal n on both. The edge's curves on the two faces stay as they are; this
// is measured, not assumed (`kinds`): the two faces' points of each parameter of the edge must
// agree within the tolerance after the move. Everywhere else (a sharp edge, or a face moved beside
// one that stays) the new edge is where the two surfaces cross. It is found section by section: in
// the plane through the point c(t) of the old edge normal to it, the point on both surfaces, by
// Newton's method on the two surfaces' parameters from the old edge's curves on them. The new
// edge keeps the old edge's parameter t, its curves are interpolated through the sections, and
// the sections are doubled until the curves between two of them are within a tenth of the
// tolerance of the section there.
//
// Vertices. A vertex is on every face around it. Where one of its edges keeps its curves on its
// faces (a tangent edge, a seam, an edge between two faces that stay), the vertex is where that
// curve, on the moved faces, meets another face around the vertex, the one it crosses most
// steeply; where the edge has no such face, the vertex goes with the edge. Else it is where three
// of its faces cross, the three whose normals are furthest from a common plane. Either way it must
// then lie on all its faces, or the move would change the topology there (four planes of a
// pyramid's apex, a face that shrinks to nothing).
//
// Certificate, apart from how the result was built: a valid solid; points inside each face at the
// distance from its source surface, on the side moved to; points along each edge at their
// distances from both faces' sources; the faces clear of each other; the volume grown by a move
// outwards and shrunk by one inwards; and where every face moves, nothing nearer the body than the
// distance. Each distance is to the unbounded surface (Newton on the foot point), as a face that
// was extended to meet its neighbours lies beyond its old bounds.

#include <BRepOffsetExt_Kept.hxx>
#include <BRepAlgoExt_Geometry.hxx>
#include <BRepAlgoExt_Properties.hxx>
#include <GeomAlgoExt_Cancel.hxx>

#include <BRepAdaptor_Surface.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <Geom_OffsetSurface.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopoDS_Compound.hxx>
#include <BRepSurfaceExt_FaceSeparation.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BOPAlgo_ArgumentAnalyzer.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>
#include <BRepTools.hxx>
#include <GeomConvert_ApproxSurface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <NCollection_Map.hxx>
#include <BRepOffset.hxx>
#include <BRepTools_Modification.hxx>
#include <BRepTools_Modifier.hxx>
#include <BRep_Builder.hxx>
#include <BRep_CurveRepresentation.hxx>
#include <BRep_PointRepresentation.hxx>
#include <BRep_TEdge.hxx>
#include <BRep_TVertex.hxx>
#include <BRep_Tool.hxx>
#include <Geom2dAPI_Interpolate.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <GeomAPI_IntSS.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomProjLib.hxx>
#include <Geom_Conic.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Line.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Geom2dAdaptor_Curve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_Surface.hxx>
#include <NCollection_DataMap.hxx>
#include <Precision.hxx>
#include <Standard_ConstructionError.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <math_Gauss.hxx>
#include <math_Matrix.hxx>
#include <math_Vector.hxx>

#include <algorithm>
#include <cmath>
#include <functional>

namespace BRepExt
{

namespace
{

const char* const op = "offsetFaces";

//! Newton's method on `z` for eval(z, f, J) = 0, whose rows are all lengths: true once |f| is
//! within `accuracy`, or where it stalls within `stalled`. A step that does not reduce
//! |f| is halved until it does (at most 30 times), so that a start further out than Newton's own
//! basin still comes in.
bool newton(math_Vector&                                                       z,
            const std::function<void(const math_Vector&, math_Vector&, math_Matrix&)>& eval,
            double                                                             accuracy,
            double                                                             stalled)
{
  const int   n = z.Length();
  math_Vector f(1, n), step(1, n), trial(1, n), g(1, n);
  math_Matrix jacobian(1, n, 1, n), unused(1, n, 1, n);
  eval(z, f, jacobian);
  for (int i = 0; i < 50; ++i)
  {
    const double size = f.Norm();
    if (size <= accuracy)
      return true;
    math_Gauss lu(jacobian);
    if (!lu.IsDone())
      return false;
    lu.Solve(f, step);
    bool lower = false;
    for (int half = 0; half <= 30 && !lower; ++half)
    {
      trial = z - step;
      eval(trial, g, unused);
      lower = g.Norm() < size;
      step *= 0.5;
    }
    // Stalled: the surfaces come no nearer than this, as two that only touch there. Within
    // `stalled` that is still a point of both.
    if (!lower)
      return size <= stalled;
    z = trial;
    eval(z, f, jacobian);
  }
  return f.Norm() <= stalled;
}

//! Where Newton stalls on two surfaces that only touch: within half the tolerance (the callers ask
//! a thousandth of it as `accuracy`) a point is still on both, leaving the other half to the
//! interpolation (a tenth) and the input's own edges.
double stalledWithin(double accuracy)
{
  return 500.0 * accuracy;
}

void put(math_Vector& f, int row, const gp_Vec& v)
{
  f(row)     = v.X();
  f(row + 1) = v.Y();
  f(row + 2) = v.Z();
}

void put(math_Matrix& j, int row, int column, const gp_Vec& v)
{
  j(row, column)     = v.X();
  j(row + 1, column) = v.Y();
  j(row + 2, column) = v.Z();
}

//! A step of an iteration toward a point near `p` small enough to stop at: far below any
//! tolerance, and above the rounding of an offset surface's derivatives (its second ones come
//! from the third of its basis), which can keep a step at 1e-11 of the size.
double settled(const gp_Pnt& p)
{
  return 1e-10 * (1.0 + p.XYZ().Modulus());
}

//! The point of the surface nearest to `p` from (u, v) on, the surface taken unbounded.
bool foot(const Handle(Geom_Surface)& surface, const gp_Pnt& p, double& u, double& v, gp_Pnt& at)
{
  for (int i = 0; i < 100; ++i)
  {
    gp_Vec su, sv, suu, svv, suv;
    surface->D2(u, v, at, su, sv, suu, svv, suv);
    const gp_Vec r(p, at);
    // Newton's step where the Hessian of the squared distance is positive definite, else the
    // Gauss-Newton step (its first-order part), which always descends.
    double a = su.Dot(su) + r.Dot(suu), b = su.Dot(sv) + r.Dot(suv), c = sv.Dot(sv) + r.Dot(svv);
    if (a <= 0.0 || a * c - b * b <= 0.0)
    {
      a = su.Dot(su);
      b = su.Dot(sv);
      c = sv.Dot(sv);
    }
    const double g = r.Dot(su), h = r.Dot(sv), det = a * c - b * b;
    if (det <= 0.0)
      return false;
    const double du = (c * g - b * h) / det, dv = (a * h - b * g) / det;
    u -= du;
    v -= dv;
    if ((su * du + sv * dv).Magnitude() <= settled(p))
    {
      at = surface->Value(u, v);
      return true;
    }
  }
  return false;
}

//! The point of the surface nearest to `p`: the foot point from (u, v) on, or where the
//! iteration has no footing (a pole, where the surface has no normal), the nearest point of the
//! surface's own domain.
bool nearest(const Handle(Geom_Surface)& surface, const gp_Pnt& p, double& u, double& v, gp_Pnt& at)
{
  const double u0 = u, v0 = v;
  if (foot(surface, p, u, v, at))
    return true;
  GeomAPI_ProjectPointOnSurf projection(p, surface);
  if (projection.NbPoints() == 0)
  {
    u = u0;
    v = v0;
    return false;
  }
  projection.LowerDistanceParameters(u, v);
  at = projection.NearestPoint();
  return true;
}

//! A face's surface after the move, and how far it moves (zero: it stays).
struct Moved
{
  Handle(Geom_Surface) surface;
  double               by          = 0.0;
  bool                 approximate = false;
};

//! An edge of the result: its curve over [first, last] of the old edge's parameter, and its
//! curves on the faces where they are new (none where it keeps them).
struct NewEdge
{
  Handle(Geom_Curve)                                                                 curve;
  double                                                                             first = 0.0;
  double                                                                             last  = 0.0;
  double                                                                             tolerance = 0.0;
  bool                                                                               crossing = false;
  NCollection_DataMap<TopoDS_Shape, Handle(Geom2d_Curve), TopTools_ShapeMapHasher> onFace;
};

struct NewVertex
{
  gp_Pnt point;
  double tolerance = 0.0;
};

//! The modification that BRepTools_Modifier applies: everything is worked out beforehand.
class Plan : public BRepTools_Modification
{
public:
  NCollection_DataMap<TopoDS_Shape, Moved, TopTools_ShapeMapHasher>     faces;
  NCollection_DataMap<TopoDS_Shape, NewEdge, TopTools_ShapeMapHasher>   edges;
  NCollection_DataMap<TopoDS_Shape, NewVertex, TopTools_ShapeMapHasher> vertices;

  bool NewSurface(const TopoDS_Face&    face,
                  Handle(Geom_Surface)& surface,
                  TopLoc_Location&      location,
                  double&               tolerance,
                  bool&                 reverseWires,
                  bool&                 reverseFace) override
  {
    const Moved* moved = faces.Seek(face);
    if (!moved)
      return false;
    surface      = moved->surface;
    location     = TopLoc_Location();
    tolerance    = BRep_Tool::Tolerance(face);
    reverseWires = false;
    reverseFace  = false;
    return true;
  }

  bool NewCurve(const TopoDS_Edge&  edge,
                Handle(Geom_Curve)& curve,
                TopLoc_Location&    location,
                double&             tolerance) override
  {
    const NewEdge* made = edges.Seek(edge);
    if (!made || made->curve.IsNull())
      return false;
    curve     = made->curve;
    location  = TopLoc_Location();
    tolerance = made->tolerance;
    return true;
  }

  bool NewPoint(const TopoDS_Vertex& vertex, gp_Pnt& point, double& tolerance) override
  {
    const NewVertex* made = vertices.Seek(vertex);
    if (!made)
      return false;
    point     = made->point;
    tolerance = made->tolerance;
    return true;
  }

  bool NewCurve2d(const TopoDS_Edge& edge,
                  const TopoDS_Face& face,
                  const TopoDS_Edge&,
                  const TopoDS_Face&,
                  Handle(Geom2d_Curve)& curve,
                  double&               tolerance) override
  {
    const NewEdge* made = edges.Seek(edge);
    if (!made)
      return false;
    const Handle(Geom2d_Curve)* onFace = made->onFace.Seek(face);
    double                      first, last;
    curve     = onFace ? *onFace : BRep_Tool::CurveOnSurface(edge, face, first, last);
    // Not shared with the input's edge: the result's geometry is its own.
    if (!curve.IsNull() && !onFace)
      curve = Handle(Geom2d_Curve)::DownCast(curve->Copy());
    tolerance = made->tolerance;
    return !curve.IsNull();
  }

  bool NewParameter(const TopoDS_Vertex& vertex,
                    const TopoDS_Edge&   edge,
                    double&              parameter,
                    double&              tolerance) override
  {
    const NewEdge* made = edges.Seek(edge);
    if (!made)
      return false;
    double first, last;
    BRep_Tool::Range(edge, first, last);
    const double old = BRep_Tool::Parameter(vertex, edge);
    parameter        = std::abs(old - first) <= std::abs(old - last) ? made->first : made->last;
    const NewVertex* moved = vertices.Seek(vertex);
    tolerance              = moved ? moved->tolerance : BRep_Tool::Tolerance(vertex);
    return true;
  }

  GeomAbs_Shape Continuity(const TopoDS_Edge& edge,
                           const TopoDS_Face& one,
                           const TopoDS_Face& other,
                           const TopoDS_Edge&,
                           const TopoDS_Face&,
                           const TopoDS_Face&) override
  {
    const NewEdge* made = edges.Seek(edge);
    return made && made->crossing ? GeomAbs_C0 : BRep_Tool::Continuity(edge, one, other);
  }
};

//! The two faces of an edge (the same face twice for a seam), or false for a free edge or one of
//! more than two faces.
bool facesOf(const NCollection_List<TopoDS_Shape>& list, TopoDS_Face& one, TopoDS_Face& other)
{
  if (list.Size() == 2)
  {
    one   = TopoDS::Face(list.First());
    other = TopoDS::Face(list.Last());
    return true;
  }
  if (list.Size() == 1)
  {
    one = other = TopoDS::Face(list.First());
    return true;
  }
  return false;
}

//! The point of `face` after the move at the parameter t of its curve of `edge`, with its
//! derivative in t.
gp_Pnt along(const Plan& plan, const TopoDS_Edge& edge, const TopoDS_Face& face, double t, gp_Vec* d)
{
  double                     first, last;
  const Handle(Geom2d_Curve) curve = BRep_Tool::CurveOnSurface(edge, face, first, last);
  gp_Pnt2d                   uv;
  gp_Vec2d                   duv;
  curve->D1(t, uv, duv);
  gp_Pnt p;
  gp_Vec su, sv;
  plan.faces.Find(face).surface->D1(uv.X(), uv.Y(), p, su, sv);
  if (d)
    *d = su * duv.X() + sv * duv.Y();
  return p;
}

//! The section of the crossing of the surfaces of `one` and `other` in the plane through c(t)
//! normal to `curve`, from the parameters `z` (u, v on each) on.
bool section(const Handle(Geom_Curve)&   curve,
             double                      t,
             const Handle(Geom_Surface)& one,
             const Handle(Geom_Surface)& other,
             math_Vector&                z,
             double                      accuracy)
{
  gp_Pnt c;
  gp_Vec tangent;
  curve->D1(t, c, tangent);
  if (tangent.Magnitude() <= gp::Resolution())
    return false;
  tangent.Normalize();
  return newton(
    z,
    [&](const math_Vector& at, math_Vector& f, math_Matrix& j) {
      gp_Pnt p, q;
      gp_Vec pu, pv, qu, qv;
      one->D1(at(1), at(2), p, pu, pv);
      other->D1(at(3), at(4), q, qu, qv);
      put(f, 1, gp_Vec(q, p));
      f(4) = gp_Vec(c, p).Dot(tangent);
      put(j, 1, 1, pu);
      put(j, 1, 2, pv);
      put(j, 1, 3, -qu);
      put(j, 1, 4, -qv);
      j(4, 1) = pu.Dot(tangent);
      j(4, 2) = pv.Dot(tangent);
      j(4, 3) = j(4, 4) = 0.0;
    },
    accuracy,
    stalledWithin(accuracy));
}

//! The parameter near `t` of the plane through c(t) normal to `curve` that holds `p`.
bool stationOf(const Handle(Geom_Curve)& curve, const gp_Pnt& p, double& t)
{
  for (int i = 0; i < 50; ++i)
  {
    gp_Pnt c;
    gp_Vec d1, d2;
    curve->D2(t, c, d1, d2);
    const double g = gp_Vec(c, p).Dot(d1), dg = gp_Vec(c, p).Dot(d2) - d1.Dot(d1);
    if (dg == 0.0)
      return false;
    t -= g / dg;
    if (std::abs(g / dg) * d1.Magnitude() <= settled(p))
      return true;
  }
  return false;
}

//! The parameter near `t` of the point of a curve q nearest to `p`; `q` gives a point and its
//! derivative.
bool nearestOn(const std::function<gp_Pnt(double, gp_Vec&)>& q, const gp_Pnt& p, double& t)
{
  for (int i = 0; i < 50; ++i)
  {
    gp_Vec       d;
    const gp_Pnt at   = q(t, d);
    const double size = d.Dot(d);
    if (size <= gp::Resolution())
      return false;
    const double step = gp_Vec(p, at).Dot(d) / size;
    t -= step;
    if (std::abs(step) * std::sqrt(size) <= settled(p))
      return true;
  }
  return false;
}

//! Cubic B-splines that follow on, each starting where the one before ends, as one: C0 where
//! they meet. `Point` is gp_Pnt or gp_Pnt2d, `Curve` the matching B-spline.
template <class Curve, class Point>
Handle(Curve) joined(const std::vector<Handle(Curve)>& pieces)
{
  if (pieces.size() == 1)
    return pieces.front();
  std::vector<Point>  poles;
  std::vector<double> knots;
  std::vector<int>    mults;
  for (size_t k = 0; k < pieces.size(); ++k)
  {
    const Handle(Curve)& piece = pieces[k];
    for (int i = k == 0 ? 1 : 2; i <= piece->NbPoles(); ++i)
      poles.push_back(piece->Pole(i));
    for (int i = k == 0 ? 1 : 2; i <= piece->NbKnots(); ++i)
    {
      const bool end = i == piece->NbKnots() && k + 1 < pieces.size();
      knots.push_back(piece->Knot(i));
      mults.push_back(end ? 3 : piece->Multiplicity(i));
    }
  }
  NCollection_Array1<Point>  p(1, static_cast<int>(poles.size()));
  NCollection_Array1<double> kn(1, static_cast<int>(knots.size()));
  NCollection_Array1<int>    m(1, static_cast<int>(mults.size()));
  for (size_t i = 0; i < poles.size(); ++i)
    p(static_cast<int>(i) + 1) = poles[i];
  for (size_t i = 0; i < knots.size(); ++i)
  {
    kn(static_cast<int>(i) + 1) = knots[i];
    m(static_cast<int>(i) + 1)  = mults[i];
  }
  return new Curve(p, kn, m, 3);
}

//! Samples of [first, last] as interpolated, a midpoint added between two samples wherever the
//! curves are further than `bound` from the exact point there: `at(t, point, uvs)` is the exact
//! point and its parameters on each face, false where there is none; `lift(f, uv)` the point of
//! face f at uv. Interpolated apart between `breaks`, where the curve is less smooth than a cubic
//! spline, and joined there. The curve and, with `faceCurves`, the curves on the faces.
bool interpolate(double                                                        first,
                 double                                                        last,
                 const std::vector<double>&                                    breaks,
                 const std::function<bool(double, gp_Pnt&, std::vector<gp_Pnt2d>&)>& at,
                 const std::function<gp_Pnt(size_t, const gp_Pnt2d&)>&         lift,
                 size_t                                                        faces,
                 bool                                                          faceCurves,
                 double                                                        bound,
                 Handle(Geom_Curve)&                                           curve,
                 std::vector<Handle(Geom2d_Curve)>&                            onFaces,
                 double&                                                       error)
{
  struct Sample
  {
    double                t;
    gp_Pnt                point;
    std::vector<gp_Pnt2d> uvs;
  };
  // Breaks nearer each other or an end than a millionth of the range are one.
  const double        apart = 1e-6 * (last - first);
  std::vector<double> ends{first};
  for (const double b : breaks)
    if (b > ends.back() + apart && b < last - apart)
      ends.push_back(b);
  ends.push_back(last);
  std::vector<Handle(Geom_BSplineCurve)>                pieces;
  std::vector<std::vector<Handle(Geom2d_BSplineCurve)>> pieces2(faceCurves ? faces : 0);
  error = 0.0;
  for (size_t piece = 0; piece + 1 < ends.size(); ++piece)
  {
    const double        a = ends[piece], b = ends[piece + 1];
    std::vector<Sample> samples;
    for (int i = 0; i <= 8; ++i)
    {
      Sample s{a + (b - a) * i / 8, {}, {}};
      if (!at(s.t, s.point, s.uvs))
      {
        return false;
      }
      samples.push_back(s);
    }
    Handle(Geom_BSplineCurve)              fitted;
    std::vector<Handle(Geom2d_BSplineCurve)> fitted2;
    double                                 worst = 0.0;
    // Judgement: 16 halvings of the first spacing or 4096 samples is where the refinement gives
    // up: a curve that a cubic still cannot follow there has a kink between the breaks.
    for (int round = 0;; ++round)
    {
      checkCancelled(op);
      if (round == 16 || samples.size() >= 4096)
      {
        return false;
      }
      const int n      = static_cast<int>(samples.size());
      auto      values = new NCollection_HArray1<double>(1, n);
      auto      p3     = new NCollection_HArray1<gp_Pnt>(1, n);
      for (int i = 0; i < n; ++i)
      {
        values->SetValue(i + 1, samples[i].t);
        p3->SetValue(i + 1, samples[i].point);
      }
      GeomAPI_Interpolate fit(p3, values, false, Precision::Confusion() * 1e-3);
      fit.Perform();
      if (!fit.IsDone())
        return false;
      fitted = fit.Curve();
      fitted2.clear();
      for (size_t f = 0; f < pieces2.size(); ++f)
      {
        auto p2 = new NCollection_HArray1<gp_Pnt2d>(1, n);
        for (int i = 0; i < n; ++i)
          p2->SetValue(i + 1, samples[i].uvs[f]);
        Geom2dAPI_Interpolate fit2(p2, values, false, 0.0);
        fit2.Perform();
        if (!fit2.IsDone())
          return false;
        fitted2.push_back(fit2.Curve());
      }
      // The midpoints: exact, and as interpolated.
      worst = 0.0;
      std::vector<Sample> refined;
      for (int i = 0; i < n; ++i)
      {
        refined.push_back(samples[i]);
        if (i + 1 == n)
          break;
        Sample mid{0.5 * (samples[i].t + samples[i + 1].t), {}, {}};
        if (!at(mid.t, mid.point, mid.uvs))
        {
          return false;
        }
        double off = mid.point.Distance(fitted->Value(mid.t));
        // Each face's curve against that face's own point: how far the faces are apart is the
        // edge's, not the interpolation's.
        for (size_t f = 0; f < fitted2.size(); ++f)
          off = std::max(off, lift(f, mid.uvs[f]).Distance(lift(f, fitted2[f]->Value(mid.t))));
        worst = std::max(worst, off);
        // A spline is global: halving one span moves its neighbours' error a little, so a span
        // just within the bound is halved too (a quarter of it; halving cuts a cubic's error 16
        // times), or the excess walks along the curve a span a round.
        if (off > 0.25 * bound)
          refined.push_back(mid);
      }
      if (worst <= bound)
        break;
      samples.swap(refined);
    }
    error = std::max(error, worst);
    pieces.push_back(fitted);
    for (size_t f = 0; f < fitted2.size(); ++f)
      pieces2[f].push_back(fitted2[f]);
  }
  curve = joined<Geom_BSplineCurve, gp_Pnt>(pieces);
  onFaces.clear();
  for (const auto& on : pieces2)
    onFaces.push_back(joined<Geom2d_BSplineCurve, gp_Pnt2d>(on));
  return true;
}

//! The parameters in [first, last] where `uv(t)`, parameters on `surface`, cross a line of the
//! surface along which it is less than C2 (less than `continuity`): a sign change of u - k or
//! v - k between `n` samples, then bisection.
void knotCrossings(const Handle(Geom_Surface)&                   surface,
                   const std::function<bool(double, gp_Pnt2d&)>& uv,
                   double                                        first,
                   double                                        last,
                   std::vector<double>&                          out,
                   GeomAbs_Shape                                 continuity = GeomAbs_C2)
{
  const GeomAdaptor_Surface adaptor(surface);
  std::vector<std::pair<int, double>> lines;
  if (adaptor.GetType() == GeomAbs_BSplineSurface || adaptor.GetType() == GeomAbs_OffsetSurface
      || adaptor.GetType() == GeomAbs_SurfaceOfExtrusion
      || adaptor.GetType() == GeomAbs_SurfaceOfRevolution)
  {
    const int nu = adaptor.NbUIntervals(continuity), nv = adaptor.NbVIntervals(continuity);
    NCollection_Array1<double> us(1, nu + 1), vs(1, nv + 1);
    adaptor.UIntervals(us, continuity);
    adaptor.VIntervals(vs, continuity);
    for (int i = 2; i <= nu; ++i)
      lines.emplace_back(0, us(i));
    for (int i = 2; i <= nv; ++i)
      lines.emplace_back(1, vs(i));
  }
  if (lines.empty())
    return;
  constexpr int         n = 32;
  std::vector<gp_Pnt2d> at(n + 1);
  for (int i = 0; i <= n; ++i)
    if (!uv(first + (last - first) * i / n, at[i]))
      return;
  for (const auto& [which, k] : lines)
    for (int i = 0; i < n; ++i)
    {
      const double s0 = at[i].Coord(which + 1) - k, s1 = at[i + 1].Coord(which + 1) - k;
      if (s0 * s1 >= 0.0)
        continue;
      double a = first + (last - first) * i / n, b = first + (last - first) * (i + 1) / n, sa = s0;
      for (int it = 0; it < 60 && b - a > 1e-15 * (std::abs(a) + std::abs(b) + 1.0); ++it)
      {
        const double m = 0.5 * (a + b);
        gp_Pnt2d     w;
        if (!uv(m, w))
          break;
        const double sm = w.Coord(which + 1) - k;
        if (sm * sa <= 0.0)
          b = m;
        else
        {
          a  = m;
          sa = sm;
        }
      }
      out.push_back(0.5 * (a + b));
    }
}

//! The parameters in (first, last) where a curve is less than C2 (`continuity` C2), or its
//! derivative is (C3).
template <class Adaptor>
void curveBreaks(const Adaptor& adaptor, GeomAbs_Shape continuity, std::vector<double>& out)
{
  const int                  n = adaptor.NbIntervals(continuity);
  NCollection_Array1<double> t(1, n + 1);
  adaptor.Intervals(t, continuity);
  for (int i = 2; i <= n; ++i)
    out.push_back(t(i));
}

//! The largest angle between the outward normals of `one` and `other` along `edge`: sampled 8
//! times in each piece of the edge where both curves on the faces are polynomial and cross no knot
//! line of either surface, and each local maximum of the samples refined by golden-section
//! search between its two neighbours. A judgement, not a bound: a bound needs the normals' rate of
//! turn along the edge, that is the faces' curvatures, bounded over it, and the hull bound on
//! curvature was three orders above the sampled curvature on blend sheets. Its measured margin is
//! in operations.md.
double largestAngle(const TopoDS_Edge& edge, const TopoDS_Face& one, const TopoDS_Face& other)
{
  double                     first, last, c0, c1;
  const Handle(Geom2d_Curve) on1 = BRep_Tool::CurveOnSurface(edge, one, c0, c1);
  const Handle(Geom2d_Curve) on2 = BRep_Tool::CurveOnSurface(edge, other, c0, c1);
  BRep_Tool::Range(edge, first, last);
  const auto angle = [&](double t) {
    gp_Vec n[2];
    int    k = 0;
    for (const auto& [face, on] : {std::pair(one, on1), std::pair(other, on2)})
    {
      const gp_Pnt2d uv = on->Value(t);
      gp_Pnt         p;
      gp_Vec         su, sv;
      BRep_Tool::Surface(face)->D1(uv.X(), uv.Y(), p, su, sv);
      n[k] = su.Crossed(sv);
      if (face.Orientation() == TopAbs_REVERSED)
        n[k].Reverse();
      ++k;
    }
    return n[0].Magnitude() > gp::Resolution() && n[1].Magnitude() > gp::Resolution()
             ? n[0].Angle(n[1])
             : 0.0;
  };
  std::vector<double> pieces{first, last};
  for (const auto& [face, on] : {std::pair(one, on1), std::pair(other, on2)})
  {
    curveBreaks(Geom2dAdaptor_Curve(on, first, last), GeomAbs_CN, pieces);
    knotCrossings(
      BRep_Tool::Surface(face),
      [&](double t, gp_Pnt2d& uv) {
        uv = on->Value(t);
        return true;
      },
      first,
      last,
      pieces,
      GeomAbs_CN);
  }
  std::sort(pieces.begin(), pieces.end());
  std::vector<double> ts;
  for (size_t i = 0; i + 1 < pieces.size(); ++i)
    for (int k = 0; k < 8; ++k)
      if (pieces[i + 1] > pieces[i])
        ts.push_back(pieces[i] + (pieces[i + 1] - pieces[i]) * k / 8);
  ts.push_back(last);
  std::vector<double> sampled(ts.size());
  for (size_t k = 0; k < ts.size(); ++k)
    sampled[k] = angle(ts[k]);
  double       largest = *std::max_element(sampled.begin(), sampled.end());
  const double golden  = 0.5 * (std::sqrt(5.0) - 1.0);
  for (size_t k = 0; k < ts.size(); ++k)
  {
    if ((k > 0 && sampled[k - 1] > sampled[k]) || (k + 1 < ts.size() && sampled[k + 1] > sampled[k]))
      continue;
    double a = ts[k > 0 ? k - 1 : k], b = ts[k + 1 < ts.size() ? k + 1 : k];
    double x = b - golden * (b - a), y = a + golden * (b - a);
    double fx = angle(x), fy = angle(y);
    for (int i = 0; i < 30; ++i)
    {
      if (fx > fy)
      {
        b  = y;
        y  = x;
        fy = fx;
        x  = b - golden * (b - a);
        fx = angle(x);
      }
      else
      {
        a  = x;
        x  = y;
        fx = fy;
        y  = a + golden * (b - a);
        fy = angle(y);
      }
    }
    largest = std::max({largest, fx, fy});
  }
  return largest;
}

//! Where two elementary surfaces (plane, cylinder, cone, sphere) cross, the crossing is a line or a
//! conic in closed form (IntAna, through GeomAPI_IntSS): the branch through `through`, turned to
//! run the way of `along`. `first` and `last` are the parameters of `from` and `to` on it, round
//! `through` (once round for a `closed` edge); `on1`, `on2` its curves on the surfaces, shifted by
//! whole periods to `uv1`, `uv2` at `through`. False where there is no such branch.
bool exactCrossing(const Handle(Geom_Surface)& s1,
                   const Handle(Geom_Surface)& s2,
                   const gp_Pnt&               through,
                   const gp_Vec&               along,
                   const gp_Pnt2d&             uv1,
                   const gp_Pnt2d&             uv2,
                   const gp_Pnt&               from,
                   const gp_Pnt&               to,
                   bool                        closed,
                   double                      near,
                   Handle(Geom_Curve)&         curve,
                   double&                     first,
                   double&                     last,
                   Handle(Geom2d_Curve)&       on1,
                   Handle(Geom2d_Curve)&       on2)
{
  const auto elementary = [](const Handle(Geom_Surface)& s) {
    const GeomAbs_SurfaceType type = GeomAdaptor_Surface(s).GetType();
    return type == GeomAbs_Plane || type == GeomAbs_Cylinder || type == GeomAbs_Cone
           || type == GeomAbs_Sphere;
  };
  if (!elementary(s1) || !elementary(s2))
    return false;
  GeomAPI_IntSS crossing(s1, s2, Precision::Confusion());
  if (!crossing.IsDone())
    return false;
  for (int i = 1; i <= crossing.NbLines(); ++i)
  {
    Handle(Geom_Curve) c = crossing.Line(i);
    if (const Handle(Geom_TrimmedCurve) trimmed = Handle(Geom_TrimmedCurve)::DownCast(c))
      c = trimmed->BasisCurve();
    if (!c->IsKind(STANDARD_TYPE(Geom_Line)) && !c->IsKind(STANDARD_TYPE(Geom_Conic)))
      continue;
    const auto parameter = [&](const gp_Pnt& q, double& u) {
      GeomAPI_ProjectPointOnCurve on(q, c);
      if (on.NbPoints() == 0 || on.LowerDistance() > near)
        return false;
      u = on.LowerDistanceParameter();
      return true;
    };
    double mid;
    if (!parameter(through, mid))
      continue;
    gp_Pnt p;
    gp_Vec d;
    c->D1(mid, p, d);
    if (d.Dot(along) < 0.0)
    {
      c = c->Reversed();
      if (!parameter(through, mid))
        return false;
    }
    if (!parameter(from, first) || !parameter(to, last))
      return false;
    if (c->IsPeriodic())
    {
      const double period = c->Period();
      first = mid - std::fmod(std::fmod(mid - first, period) + period, period);
      last  = closed ? first + period : mid + std::fmod(std::fmod(last - mid, period) + period, period);
    }
    if (!(first < mid && mid < last))
      return false;
    Handle(Geom2d_Curve) curves[2];
    const Handle(Geom_Surface) surfaces[2] = {s1, s2};
    const gp_Pnt2d             uvs[2]      = {uv1, uv2};
    for (int k = 0; k < 2; ++k)
    {
      curves[k] = GeomProjLib::Curve2d(c, first, last, surfaces[k]);
      if (curves[k].IsNull())
        return false;
      const gp_Pnt2d w = curves[k]->Value(mid);
      gp_Vec2d       shift(0.0, 0.0);
      if (surfaces[k]->IsUPeriodic())
        shift.SetX(surfaces[k]->UPeriod()
                   * std::round((uvs[k].X() - w.X()) / surfaces[k]->UPeriod()));
      if (surfaces[k]->IsVPeriodic())
        shift.SetY(surfaces[k]->VPeriod()
                   * std::round((uvs[k].Y() - w.Y()) / surfaces[k]->VPeriod()));
      curves[k]->Translate(shift);
    }
    curve = c;
    on1   = curves[0];
    on2   = curves[1];
    return true;
  }
  return false;
}

//! The curve of points q(t) in closed form where the old edge's curve c is a line or a circle and
//! q is c moved with its parameter kept: a line translated, or a circle about the same axis with
//! another radius and height (an edge along a plane, cylinder, cone, sphere or torus round its
//! axis moves so). Null unless q is that curve within `accuracy` at 17 parameters of
//! [first, last]. Display meshes and Booleans then see the same line or circle the input had, not
//! a spline through it.
Handle(Geom_Curve) movedLineOrCircle(const TopoDS_Edge&                     edge,
                                     const std::function<gp_Pnt(double)>& q,
                                     double                               first,
                                     double                               last,
                                     double                               accuracy)
{
  double             f, l;
  Handle(Geom_Curve) old = BRep_Tool::Curve(edge, f, l);
  while (!old.IsNull() && old->IsKind(STANDARD_TYPE(Geom_TrimmedCurve)))
    old = Handle(Geom_TrimmedCurve)::DownCast(old)->BasisCurve();
  Handle(Geom_Curve) moved;
  if (const Handle(Geom_Line) line = Handle(Geom_Line)::DownCast(old))
    moved = Handle(Geom_Curve)::DownCast(
      line->Translated(gp_Vec(line->Value(first), q(first))));
  else if (const Handle(Geom_Circle) circle = Handle(Geom_Circle)::DownCast(old))
  {
    const gp_Ax2 axes   = circle->Position();
    const gp_Vec offset(axes.Location(), q(first));
    const double height = offset.Dot(gp_Vec(axes.Direction()));
    const double radius = (offset - gp_Vec(axes.Direction()) * height).Magnitude();
    if (radius <= gp::Resolution())
      return {};
    gp_Ax2 moved2 = axes;
    moved2.SetLocation(axes.Location().Translated(gp_Vec(axes.Direction()) * height));
    moved = new Geom_Circle(moved2, radius);
  }
  if (moved.IsNull())
    return {};
  for (int k = 0; k <= 16; ++k)
  {
    const double t = first + (last - first) * k / 16;
    if (moved->Value(t).Distance(q(t)) > accuracy)
      return {};
  }
  return moved;
}

//! `made`, a face of `shape` on `exact`, put on a spline over its box in its parameters within
//! `bound` of `exact` (GeomConvert_ApproxSurface), so that its edges' curves on it stay: they, the
//! vertices' points on them and the edges' records of how smoothly two faces meet across them are
//! moved onto the spline, and the edges' and vertices' tolerances grow by the error. False, the face
//! untouched, where no spline is within `bound`.
bool asSpline(const TopoDS_Shape&         shape,
              const TopoDS_Face&          made,
              const Handle(Geom_Surface)& exact,
              double                      bound)
{
  BRep_Builder builder;
  double       u0, u1, v0, v1;
  BRepTools::UVBounds(made, u0, u1, v0, v1);
  Handle(Geom_BSplineSurface) spline;
  double                      error = HUGE_VAL;
  try
  {
    // Judgement: degree 9 and 1000 segments each way at most (AdvApp2Var's own scale).
    GeomConvert_ApproxSurface approx(new Geom_RectangularTrimmedSurface(exact, u0, u1, v0, v1),
                                     bound,
                                     GeomAbs_C0,
                                     GeomAbs_C0,
                                     9,
                                     9,
                                     1000,
                                     0);
    if (approx.IsDone() && approx.HasResult())
    {
      spline = approx.Surface();
      error  = approx.MaxError();
    }
  }
  catch (const Standard_Failure&)
  {
  }
  if (spline.IsNull() || !(error <= bound))
    return false;
  // The edges' curves on the face are kept by its surface: read them off before the surface
  // changes, and put them back on the new one (a seam's two).
  struct OnFace
  {
    TopoDS_Edge          edge;
    Handle(Geom2d_Curve) one, two;
    double               first, last;
  };
  std::vector<OnFace> curves;
  for (TopExp_Explorer e(made, TopAbs_EDGE); e.More(); e.Next())
  {
    const TopoDS_Edge edge = TopoDS::Edge(e.Current());
    OnFace            on{edge, {}, {}, 0.0, 0.0};
    on.one = BRep_Tool::CurveOnSurface(TopoDS::Edge(edge.Oriented(TopAbs_FORWARD)),
                                       made,
                                       on.first,
                                       on.last);
    if (BRep_Tool::IsClosed(edge, made))
      on.two = BRep_Tool::CurveOnSurface(TopoDS::Edge(edge.Oriented(TopAbs_REVERSED)),
                                         made,
                                         on.first,
                                         on.last);
    if (std::none_of(curves.begin(), curves.end(), [&](const OnFace& o) {
          return o.edge.IsSame(edge);
        }))
      curves.push_back(on);
  }
  // The old curves go with the old surface, which a saved file would otherwise still carry.
  for (const OnFace& on : curves)
    if (on.two.IsNull())
      builder.UpdateEdge(on.edge, Handle(Geom2d_Curve)(), made, BRep_Tool::Tolerance(on.edge));
    else
      builder.UpdateEdge(on.edge,
                         Handle(Geom2d_Curve)(),
                         Handle(Geom2d_Curve)(),
                         made,
                         BRep_Tool::Tolerance(on.edge));
  // So do the vertices' points on its curves (where an edge has no curve of its own, a pole),
  // and the edges' records of how smoothly two faces meet across them.
  const Handle(Geom_Surface) old = BRep_Tool::Surface(made);
  struct Smooth
  {
    TopoDS_Edge   edge;
    GeomAbs_Shape continuity;
  };
  std::vector<Smooth> smooth;
  for (const OnFace& on : curves)
  {
    auto& reps = Handle(BRep_TEdge)::DownCast(on.edge.TShape())->ChangeCurves();
    for (NCollection_List<Handle(BRep_CurveRepresentation)>::Iterator it(reps); it.More();)
      if (it.Value()->IsRegularity()
          && (it.Value()->Surface() == old || it.Value()->Surface2() == old))
      {
        smooth.push_back({on.edge, it.Value()->Continuity()});
        reps.Remove(it);
      }
      else
        it.Next();
  }
  for (TopExp_Explorer v(made, TopAbs_VERTEX); v.More(); v.Next())
  {
    auto& points = Handle(BRep_TVertex)::DownCast(v.Current().TShape())->ChangePoints();
    for (NCollection_List<Handle(BRep_PointRepresentation)>::Iterator it(points); it.More();)
      if ((it.Value()->IsPointOnCurveOnSurface() || it.Value()->IsPointOnSurface())
          && it.Value()->Surface() == old)
        points.Remove(it);
      else
        it.Next();
  }
  builder.UpdateFace(made, spline, TopLoc_Location(), BRep_Tool::Tolerance(made));
  for (const OnFace& on : curves)
  {
    const double edgeTolerance = BRep_Tool::Tolerance(on.edge) + error;
    if (on.two.IsNull())
      builder.UpdateEdge(on.edge, on.one, made, edgeTolerance);
    else
      builder.UpdateEdge(on.edge, on.one, on.two, made, edgeTolerance);
    builder.Range(on.edge, made, on.first, on.last);
    if (BRep_Tool::Degenerated(on.edge))
      for (TopExp_Explorer v(on.edge, TopAbs_VERTEX); v.More(); v.Next())
        builder.UpdateVertex(TopoDS::Vertex(v.Current()),
                             v.Current().Orientation() == TopAbs_FORWARD ? on.first : on.last,
                             on.edge,
                             made,
                             BRep_Tool::Tolerance(TopoDS::Vertex(v.Current())));
  }
  for (TopExp_Explorer v(made, TopAbs_VERTEX); v.More(); v.Next())
    builder.UpdateVertex(TopoDS::Vertex(v.Current()),
                         BRep_Tool::Tolerance(TopoDS::Vertex(v.Current())) + error);
  for (const Smooth& at : smooth)
  {
    // The face across the edge, as the result now has it.
    TopoDS_Face across = made;
    for (TopExp_Explorer f(shape, TopAbs_FACE); f.More(); f.Next())
      if (!f.Current().IsSame(made))
        for (TopExp_Explorer e(f.Current(), TopAbs_EDGE); e.More(); e.Next())
          if (e.Current().IsSame(at.edge))
            across = TopoDS::Face(f.Current());
    builder.Continuity(at.edge, made, across, at.continuity);
  }
  return true;
}

//! Whether the faces of `result` cross each other: BOPAlgo's self-interference check on what the
//! separation pre-pass (BRepSurfaceExt_FaceSeparation.hxx) cannot prove clear, the pre-pass run on
//! proxies. An offset surface, which the pre-pass does not read, has a spline proxy within `eps` of
//! it (asSpline) on a copy, its face's tolerance grown by `eps`: two proxies proven further apart
//! than their tolerances have their exact faces apart too. That a proxy is one to one says nothing
//! of its exact face, so an offset face goes to the check by itself as well where the pre-pass
//! has nothing else of it to check. The result keeps its exact faces.
bool crossesItself(const TopoDS_Shape& result, double eps)
{
  const TopoDS_Shape proxy = BRepBuilderAPI_Copy(result, true, false).Shape();
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> exact, proxies;
  TopExp::MapShapes(result, TopAbs_FACE, exact);
  TopExp::MapShapes(proxy, TopAbs_FACE, proxies);
  if (exact.Extent() != proxies.Extent())
    return intersectsItself(result);
  BRep_Builder      builder;
  std::vector<bool> offset(exact.Extent() + 1, false);
  for (int i = 1; i <= proxies.Extent(); ++i)
  {
    const TopoDS_Face    face    = TopoDS::Face(proxies(i));
    Handle(Geom_Surface) surface = BRep_Tool::Surface(face);
    while (surface->IsKind(STANDARD_TYPE(Geom_RectangularTrimmedSurface)))
      surface = Handle(Geom_RectangularTrimmedSurface)::DownCast(surface)->BasisSurface();
    if (!surface->IsKind(STANDARD_TYPE(Geom_OffsetSurface)))
      continue;
    if (!asSpline(proxy, face, BRep_Tool::Surface(face), eps))
      return intersectsItself(result);
    builder.UpdateFace(face, BRep_Tool::Tolerance(face) + eps);
    offset[i] = true;
  }
  const auto check = [](const TopoDS_Shape& part) {
    BOPAlgo_ArgumentAnalyzer analyzer;
    analyzer.SetShape1(part);
    analyzer.ArgumentTypeMode() = false;
    analyzer.SmallEdgeMode()    = false;
    analyzer.SelfInterMode()    = true;
    analyzer.Perform();
    return analyzer.HasFaulty();
  };
  std::vector<bool> checked(exact.Extent() + 1, false);
  for (const TopoDS_Shape& part : unprovenParts(proxy))
  {
    TopoDS_Compound counterpart;
    builder.MakeCompound(counterpart);
    for (TopExp_Explorer f(part, TopAbs_FACE); f.More(); f.Next())
    {
      const int i = proxies.FindIndex(f.Current());
      if (i == 0)
        return intersectsItself(result);
      builder.Add(counterpart, exact(i));
      checked[i] = true;
    }
    if (check(part.ShapeType() == TopAbs_FACE ? exact(proxies.FindIndex(part)) : counterpart))
      return true;
  }
  for (int i = 1; i <= exact.Extent(); ++i)
    if (offset[i] && !checked[i] && check(exact(i)))
      return true;
  return false;
}

} // namespace

TopoDS_Shape keptOffset(const TopoDS_Shape& body,
                        const FaceSet&      moved,
                        double              distance,
                        double              tolerance,
                        Lineage&            lineage,
                        std::string&        why)
{
  Handle(Plan) plan = new Plan;
  const double accuracy = 1e-3 * tolerance;

  // Faces.
  for (TopExp_Explorer all(body, TopAbs_FACE); all.More(); all.Next())
  {
    const TopoDS_Face& face = TopoDS::Face(all.Current());
    if (plan->faces.IsBound(face))
      continue;
    Moved m;
    m.surface = BRep_Tool::Surface(face);
    if (moved.Contains(face))
    {
      m.by = distance;
      BRepOffset_Status status;
      // Checked for continuity (allowC0 false) as the binary BREP reader builds an offset surface
      // (BinTools_SurfaceSet), so that what is saved loads again. A basis that the check finds
      // not G1 at a knot (a blend's sheet joined parametrically C0 where its parts meet) has its
      // offset built unchecked and, once the face's extent is known, replaced by a spline within
      // a tenth of the tolerance (below), which loads.
      const Handle(Geom_Surface) basis =
        BRepOffset::CollapseSingularities(m.surface, face, Precision::Confusion());
      const double along = face.Orientation() == TopAbs_REVERSED ? -distance : distance;
      try
      {
        m.surface = BRepOffset::Surface(basis, along, status, false);
      }
      catch (const Standard_ConstructionError&)
      {
        m.surface     = BRepOffset::Surface(basis, along, status, true);
        m.approximate = true;
      }
      if (status != BRepOffset_Good || m.surface.IsNull())
      {
        why = "a face is bent more tightly than the distance; its offset surface turns inside out";
        return {};
      }
    }
    plan->faces.Bind(face, m);
  }

  Adjacency facesOfEdge, edgesOfVertex, facesOfVertex;
  TopExp::MapShapesAndUniqueAncestors(body, TopAbs_EDGE, TopAbs_FACE, facesOfEdge);
  TopExp::MapShapesAndUniqueAncestors(body, TopAbs_VERTEX, TopAbs_EDGE, edgesOfVertex);
  TopExp::MapShapesAndUniqueAncestors(body, TopAbs_VERTEX, TopAbs_FACE, facesOfVertex);

  // Edges that keep their curves on their faces: the faces' points of each parameter agree after
  // the move, within the tolerance or the edge's own, whichever is the larger: two faces tangent
  // along an edge only to an angle theta (an input that is G1 to a few 1e-5 rad, its edges' tolerances
  // that wide) come apart by |d| theta more, and the result is held to no more than the input
  // held. A degenerate edge (a pole) must stay a point. `spread` is |d| theta_max (largestAngle).
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>         keeps;
  NCollection_DataMap<TopoDS_Shape, double, TopTools_ShapeMapHasher> spread;
  for (int i = 1; i <= facesOfEdge.Extent(); ++i)
  {
    const TopoDS_Edge& edge = TopoDS::Edge(facesOfEdge.FindKey(i));
    TopoDS_Face        one, other;
    if (!facesOf(facesOfEdge(i), one, other))
    {
      why = "an edge is not between two faces";
      return {};
    }
    double first, last;
    BRep_Tool::Range(edge, first, last);
    // How far apart the two faces' points of the edge are before the move, and how much further
    // after it.
    double gap = 0.0;
    for (int k = 0; k <= 8; ++k)
    {
      const double t = first + (last - first) * k / 8;
      double       c0, c1;
      const gp_Pnt2d w1 = BRep_Tool::CurveOnSurface(edge, one, c0, c1)->Value(t);
      const gp_Pnt2d w2 = BRep_Tool::CurveOnSurface(edge, other, c0, c1)->Value(t);
      const double   before =
        BRep_Tool::Surface(one)->Value(w1.X(), w1.Y()).Distance(BRep_Tool::Surface(other)->Value(w2.X(), w2.Y()));
      gap = std::max(gap,
                     along(*plan, edge, one, t, nullptr).Distance(along(*plan, edge, other, t, nullptr))
                       - before);
    }
    if (BRep_Tool::Degenerated(edge))
    {
      for (int k = 1; k <= 8; ++k)
        gap = std::max(gap, along(*plan, edge, one, first, nullptr)
                              .Distance(along(*plan, edge, one, first + (last - first) * k / 8, nullptr)));
      if (gap > tolerance)
      {
        why = "a pole of a face would open into a curve";
        return {};
      }
    }
    if (gap <= std::max(tolerance, BRep_Tool::Tolerance(edge)) || one.IsSame(other))
    {
      keeps.Add(edge);
      spread.Bind(edge, one.IsSame(other) ? 0.0 : std::abs(distance) * largestAngle(edge, one, other));
    }
  }

  // Vertices.
  for (int i = 1; i <= facesOfVertex.Extent(); ++i)
  {
    const TopoDS_Vertex&                  vertex = TopoDS::Vertex(facesOfVertex.FindKey(i));
    const NCollection_List<TopoDS_Shape>& around = facesOfVertex(i);
    bool                                  any    = false;
    for (const TopoDS_Shape& face : around)
      any = any || plan->faces.Find(face).by != 0.0;
    if (!any)
      continue;
    gp_Pnt       point;
    bool         found = false;
    // On an edge that keeps its curves, crossed with the face it crosses most steeply of those that
    // meet another at a sharp edge here: a face tangent to the others at the vertex has their
    // normal there and adds nothing. Where all are tangent, the vertex moves along that normal.
    NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> sharp;
    for (const TopoDS_Shape& e : edgesOfVertex.FindFromKey(vertex))
      if (!keeps.Contains(e))
        for (const TopoDS_Shape& f : facesOfEdge.FindFromKey(e))
          sharp.Add(f);
    for (const TopoDS_Shape& e : edgesOfVertex.FindFromKey(vertex))
    {
      const TopoDS_Edge& edge = TopoDS::Edge(e);
      if (!keeps.Contains(edge) || BRep_Tool::Degenerated(edge))
        continue;
      TopoDS_Face one, other;
      facesOf(facesOfEdge.FindFromKey(edge), one, other);
      const double t0 = BRep_Tool::Parameter(vertex, edge);
      gp_Vec       d;
      along(*plan, edge, one, t0, &d);
      TopoDS_Face steepest;
      double      steep = 0.0;
      for (const TopoDS_Shape& f : around)
      {
        if (f.IsSame(one) || f.IsSame(other) || !sharp.Contains(f))
          continue;
        const TopoDS_Face face = TopoDS::Face(f);
        const gp_Pnt2d    uv   = BRep_Tool::Parameters(vertex, face);
        gp_Pnt            p;
        gp_Vec            su, sv;
        plan->faces.Find(face).surface->D1(uv.X(), uv.Y(), p, su, sv);
        const gp_Vec normal = su.Crossed(sv);
        if (normal.Magnitude() <= gp::Resolution() || d.Magnitude() <= gp::Resolution())
          continue;
        const double s = std::abs(normal.Normalized().Dot(d.Normalized()));
        if (s > steep)
        {
          steep    = s;
          steepest = face;
        }
      }
      if (steepest.IsNull())
      {
        point = along(*plan, edge, one, t0, nullptr);
        found = true;
        break;
      }
      const gp_Pnt2d       uv      = BRep_Tool::Parameters(vertex, steepest);
      const Handle(Geom_Surface)& crossed = plan->faces.Find(steepest).surface;
      math_Vector          z(1, 3);
      z(1) = t0;
      z(2) = uv.X();
      z(3) = uv.Y();
      found = newton(
        z,
        [&](const math_Vector& at, math_Vector& f, math_Matrix& j) {
          gp_Vec       dq;
          const gp_Pnt q = along(*plan, edge, one, at(1), &dq);
          gp_Pnt       s;
          gp_Vec       su, sv;
          crossed->D1(at(2), at(3), s, su, sv);
          put(f, 1, gp_Vec(s, q));
          put(j, 1, 1, dq);
          put(j, 1, 2, -su);
          put(j, 1, 3, -sv);
        },
        accuracy,
        stalledWithin(accuracy));
      if (found)
        point = along(*plan, edge, one, z(1), nullptr);
      break;
    }
    // Between two edges where the same two faces cross: in the section of one of them there.
    if (!found && around.Size() == 2)
      for (const TopoDS_Shape& e : edgesOfVertex.FindFromKey(vertex))
      {
        const TopoDS_Edge&       edge = TopoDS::Edge(e);
        double                   f0, l0;
        const Handle(Geom_Curve) old = BRep_Tool::Curve(edge, f0, l0);
        if (old.IsNull())
          continue;
        const TopoDS_Face one   = TopoDS::Face(around.First());
        const TopoDS_Face other = TopoDS::Face(around.Last());
        const gp_Pnt2d    w1 = BRep_Tool::Parameters(vertex, one), w2 = BRep_Tool::Parameters(vertex, other);
        math_Vector       z(1, 4);
        z(1)  = w1.X();
        z(2)  = w1.Y();
        z(3)  = w2.X();
        z(4)  = w2.Y();
        found = section(old,
                        BRep_Tool::Parameter(vertex, edge),
                        plan->faces.Find(one).surface,
                        plan->faces.Find(other).surface,
                        z,
                        accuracy);
        if (found)
          point = plan->faces.Find(one).surface->Value(z(1), z(2));
        break;
      }
    // Else where three faces cross.
    if (!found && around.Size() >= 3)
    {
      std::vector<TopoDS_Face> list;
      std::vector<gp_Dir>      normals;
      for (const TopoDS_Shape& f : around)
      {
        const TopoDS_Face face = TopoDS::Face(f);
        const gp_Pnt2d    uv   = BRep_Tool::Parameters(vertex, face);
        gp_Pnt            p;
        gp_Vec            su, sv;
        plan->faces.Find(face).surface->D1(uv.X(), uv.Y(), p, su, sv);
        if (su.Crossed(sv).Magnitude() <= gp::Resolution())
          continue;
        list.push_back(face);
        normals.emplace_back(su.Crossed(sv));
      }
      double best = 0.0;
      int    a = -1, b = -1, c = -1;
      for (size_t x = 0; x < list.size(); ++x)
        for (size_t y = x + 1; y < list.size(); ++y)
          for (size_t w = y + 1; w < list.size(); ++w)
          {
            const double det =
              std::abs(gp_Vec(normals[x]).Dot(gp_Vec(normals[y]).Crossed(gp_Vec(normals[w]))));
            if (det > best)
            {
              best = det;
              a    = static_cast<int>(x);
              b    = static_cast<int>(y);
              c    = static_cast<int>(w);
            }
          }
      if (a >= 0)
      {
        const TopoDS_Face           three[3] = {list[a], list[b], list[c]};
        const Handle(Geom_Surface)* s[3];
        math_Vector                 z(1, 6);
        for (int k = 0; k < 3; ++k)
        {
          s[k]              = &plan->faces.Find(three[k]).surface;
          const gp_Pnt2d uv = BRep_Tool::Parameters(vertex, three[k]);
          z(2 * k + 1)      = uv.X();
          z(2 * k + 2)      = uv.Y();
        }
        found = newton(
          z,
          [&](const math_Vector& at, math_Vector& f, math_Matrix& j) {
            gp_Pnt p[3];
            gp_Vec su[3], sv[3];
            for (int k = 0; k < 3; ++k)
              (*s[k])->D1(at(2 * k + 1), at(2 * k + 2), p[k], su[k], sv[k]);
            j.Init(0.0);
            put(f, 1, gp_Vec(p[1], p[0]));
            put(f, 4, gp_Vec(p[2], p[0]));
            for (int row : {1, 4})
            {
              put(j, row, 1, su[0]);
              put(j, row, 2, sv[0]);
            }
            put(j, 1, 3, -su[1]);
            put(j, 1, 4, -sv[1]);
            put(j, 4, 5, -su[2]);
            put(j, 4, 6, -sv[2]);
          },
          accuracy,
          stalledWithin(accuracy));
        if (found)
          point = (*s[0])->Value(z(1), z(2));
      }
    }
    if (!found)
    {
      why = "a vertex has no place where its faces meet after the move";
      return {};
    }
    // On every face around it, within the tolerance or what the input's vertex and the spread of
    // its tangent edges allow.
    double carried = 0.0;
    for (const TopoDS_Shape& e : edgesOfVertex.FindFromKey(vertex))
      if (const double* d = spread.Seek(e))
        carried = std::max(carried, BRep_Tool::Tolerance(vertex) + *d);
    double off = 0.0;
    for (const TopoDS_Shape& f : around)
    {
      const TopoDS_Face face = TopoDS::Face(f);
      gp_Pnt2d          uv   = BRep_Tool::Parameters(vertex, face);
      double            u = uv.X(), v = uv.Y();
      gp_Pnt            at;
      if (!nearest(plan->faces.Find(face).surface, point, u, v, at)
          || at.Distance(point) > std::max(tolerance, carried))
      {
        why = "the faces round a vertex do not meet in one point after the move; the topology "
              "would change there";
        return {};
      }
      off = std::max(off, at.Distance(point));
    }
    // No finer than the input's: its edges meet within that.
    plan->vertices.Bind(vertex, {point, std::max(BRep_Tool::Tolerance(vertex), 1.001 * off)});
  }

  // Edges.
  for (int i = 1; i <= facesOfEdge.Extent(); ++i)
  {
    checkCancelled(op);
    const TopoDS_Edge& edge = TopoDS::Edge(facesOfEdge.FindKey(i));
    TopoDS_Face        one, other;
    facesOf(facesOfEdge(i), one, other);
    const Moved& a = plan->faces.Find(one);
    const Moved& b = plan->faces.Find(other);
    TopoDS_Vertex start, end;
    TopExp::Vertices(edge, start, end);
    const bool startMoves = !start.IsNull() && plan->vertices.IsBound(start);
    const bool endMoves   = !end.IsNull() && plan->vertices.IsBound(end);
    if (a.by == 0.0 && b.by == 0.0 && !startMoves && !endMoves)
      continue;
    double first, last;
    BRep_Tool::Range(edge, first, last);
    NewEdge made;
    made.first = first;
    made.last  = last;
    if (BRep_Tool::Degenerated(edge))
    {
      made.tolerance = BRep_Tool::Tolerance(edge);
      plan->edges.Bind(edge, made);
      continue;
    }
    const bool closed = !start.IsNull() && start.IsSame(end);
    if (keeps.Contains(edge))
    {
      // The curve of the first face, moved; the vertices' parameters on it.
      const auto q = [&](double t, gp_Vec& d) { return along(*plan, edge, one, t, &d); };
      if (startMoves && !nearestOn(q, plan->vertices.Find(start).point, made.first))
      {
        why = "a vertex is not on the curve of its edge after the move";
        return {};
      }
      if (endMoves && !nearestOn(q, plan->vertices.Find(end).point, made.last))
      {
        why = "a vertex is not on the curve of its edge after the move";
        return {};
      }
      if (closed)
        made.last = made.first + (last - first);
      // The curves on the faces stay unless the edge now runs past their ends; then they are
      // interpolated too, through their own continuation.
      double                     c0, c1;
      const Handle(Geom2d_Curve) onOne   = BRep_Tool::CurveOnSurface(edge, one, c0, c1);
      const Handle(Geom2d_Curve) onOther = BRep_Tool::CurveOnSurface(edge, other, c0, c1);
      const auto                 within  = [&](const Handle(Geom2d_Curve)& c) {
        return c->IsPeriodic()
               || (c->FirstParameter() <= made.first + Precision::PConfusion()
                   && made.last <= c->LastParameter() + Precision::PConfusion());
      };
      const bool refit = !within(onOne) || !within(onOther);
      if (refit && one.IsSame(other))
      {
        why = "a seam runs past the ends of its curves on its face after the move";
        return {};
      }
      // Breaks where a curve on a face, or the moved surface along it, is less than C2.
      std::vector<double> breaks;
      for (const auto& [face, on] : {std::pair(one, onOne), std::pair(other, onOther)})
      {
        curveBreaks(Geom2dAdaptor_Curve(on), GeomAbs_C2, breaks);
        knotCrossings(
          plan->faces.Find(face).surface,
          [&](double t, gp_Pnt2d& uv) {
            uv = on->Value(t);
            return true;
          },
          made.first,
          made.last,
          breaks);
      }
      std::sort(breaks.begin(), breaks.end());
      std::vector<Handle(Geom2d_Curve)> curves;
      double                            error = 0.0, gap = 0.0;
      if (!refit)
        made.curve = movedLineOrCircle(
          edge,
          [&](double t) { return along(*plan, edge, one, t, nullptr); },
          made.first,
          made.last,
          accuracy);
      if (!made.curve.IsNull())
        for (int k = 0; k <= 16; ++k)
        {
          const double t = made.first + (made.last - made.first) * k / 16;
          gap = std::max(gap, made.curve->Value(t).Distance(along(*plan, edge, other, t, nullptr)));
        }
      const bool fitted = !made.curve.IsNull() || interpolate(
        made.first,
        made.last,
        breaks,
        [&](double t, gp_Pnt& p, std::vector<gp_Pnt2d>& uvs) {
          p   = along(*plan, edge, one, t, nullptr);
          gap = std::max(gap, p.Distance(along(*plan, edge, other, t, nullptr)));
          uvs = {onOne->Value(t), onOther->Value(t)};
          return true;
        },
        [&](size_t f, const gp_Pnt2d& uv) {
          return plan->faces.Find(f == 0 ? one : other).surface->Value(uv.X(), uv.Y());
        },
        2,
        refit,
        0.1 * tolerance,
        made.curve,
        curves,
        error);
      if (!fitted)
      {
        why = "the offset of a tangent edge could not be fitted within the tolerance";
        return {};
      }
      if (refit)
      {
        made.onFace.Bind(one, curves[0]);
        made.onFace.Bind(other, curves[1]);
      }
      made.tolerance =
        std::max({BRep_Tool::Tolerance(edge) + spread.Find(edge), error, gap});
    }
    else
    {
      // Where the two surfaces cross, section by section.
      double                   f0, l0;
      const Handle(Geom_Curve) old = BRep_Tool::Curve(edge, f0, l0);
      if (old.IsNull())
      {
        why = "an edge has no curve";
        return {};
      }
      const auto stationAt = [&](const TopoDS_Vertex& vertex, double& t) {
        return stationOf(old, plan->vertices.IsBound(vertex) ? plan->vertices.Find(vertex).point
                                                             : BRep_Tool::Pnt(vertex),
                         t);
      };
      if (!start.IsNull() && !stationAt(start, made.first))
      {
        why = "a vertex is off the section planes of its edge";
        return {};
      }
      if (!end.IsNull())
      {
        if (closed)
          made.last = made.first + (last - first);
        else if (!stationAt(end, made.last))
        {
          why = "a vertex is off the section planes of its edge";
          return {};
        }
      }
      if (!(made.last > made.first))
      {
        why = "an edge shrinks to nothing";
        return {};
      }
      double                     c0, c1;
      const Handle(Geom2d_Curve) onOne   = BRep_Tool::CurveOnSurface(edge, one, c0, c1);
      const Handle(Geom2d_Curve) onOther = BRep_Tool::CurveOnSurface(edge, other, c0, c1);
      const auto station = [&](double t, math_Vector& z) {
        const gp_Pnt2d w1 = onOne->Value(t), w2 = onOther->Value(t);
        z(1)              = w1.X();
        z(2)              = w1.Y();
        z(3)              = w2.X();
        z(4)              = w2.Y();
        return section(old, t, a.surface, b.surface, z, accuracy);
      };
      // Two elementary surfaces cross along a line or a conic: that, exactly.
      std::vector<Handle(Geom2d_Curve)> curves;
      double                            error = 0.0;
      {
        const double tm = 0.5 * (made.first + made.last);
        math_Vector  z(1, 4);
        gp_Pnt       c;
        gp_Vec       along;
        old->D1(tm, c, along);
        // Into copies: a crossing that fails part way leaves the edge's parameters as they were.
        Handle(Geom2d_Curve) on1, on2;
        Handle(Geom_Curve)   exact;
        double               e0 = made.first, e1 = made.last;
        if (station(tm, z)
            && exactCrossing(a.surface,
                             b.surface,
                             a.surface->Value(z(1), z(2)),
                             along,
                             gp_Pnt2d(z(1), z(2)),
                             gp_Pnt2d(z(3), z(4)),
                             plan->vertices.IsBound(start) ? plan->vertices.Find(start).point
                                                           : BRep_Tool::Pnt(start),
                             plan->vertices.IsBound(end) ? plan->vertices.Find(end).point
                                                         : BRep_Tool::Pnt(end),
                             closed,
                             tolerance,
                             exact,
                             e0,
                             e1,
                             on1,
                             on2))
        {
          made.curve = exact;
          made.first = e0;
          made.last  = e1;
          curves     = {on1, on2};
        }
      }
      // Else by sections. Breaks where the old edge is less than C3 (its tangent, which sets the
      // section planes, less than C2), and where the crossing runs over a line of either surface
      // where it is less than C2.
      std::vector<double> breaks;
      if (curves.empty())
        curveBreaks(GeomAdaptor_Curve(old), GeomAbs_C3, breaks);
      for (int f = 0; f < 2 && curves.empty(); ++f)
        knotCrossings(
          (f == 0 ? a : b).surface,
          [&](double t, gp_Pnt2d& uv) {
            math_Vector z(1, 4);
            if (!station(t, z))
              return false;
            uv = gp_Pnt2d(z(2 * f + 1), z(2 * f + 2));
            return true;
          },
          made.first,
          made.last,
          breaks);
      std::sort(breaks.begin(), breaks.end());
      const bool                        fitted = !curves.empty() || interpolate(
        made.first,
        made.last,
        breaks,
        [&](double t, gp_Pnt& p, std::vector<gp_Pnt2d>& uvs) {
          math_Vector z(1, 4);
          if (!station(t, z))
            return false;
          p   = a.surface->Value(z(1), z(2));
          uvs = {gp_Pnt2d(z(1), z(2)), gp_Pnt2d(z(3), z(4))};
          return true;
        },
        [&](size_t f, const gp_Pnt2d& uv) {
          return (f == 0 ? a : b).surface->Value(uv.X(), uv.Y());
        },
        2,
        true,
        0.1 * tolerance,
        made.curve,
        curves,
        error);
      if (!fitted)
      {
        why = "the crossing of two faces after the move could not be followed along their edge";
        return {};
      }
      made.crossing = true;
      made.onFace.Bind(one, curves[0]);
      made.onFace.Bind(other, curves[1]);
      // How far the curves on the faces are from the curve, measured.
      double off = error;
      for (int k = 0; k <= 64; ++k)
      {
        const double t = made.first + (made.last - made.first) * k / 64;
        const gp_Pnt p = made.curve->Value(t);
        const gp_Pnt2d w1 = curves[0]->Value(t), w2 = curves[1]->Value(t);
        off = std::max({off,
                        p.Distance(a.surface->Value(w1.X(), w1.Y())),
                        p.Distance(b.surface->Value(w2.X(), w2.Y()))});
      }
      made.tolerance = std::max(BRep_Tool::Tolerance(edge), off);
    }
    plan->edges.Bind(edge, made);
  }
  // A vertex's tolerance covers the ends of its edges.
  for (NCollection_DataMap<TopoDS_Shape, NewEdge, TopTools_ShapeMapHasher>::Iterator it(plan->edges);
       it.More();
       it.Next())
  {
    const NewEdge& made = it.Value();
    if (made.curve.IsNull())
      continue;
    TopoDS_Vertex start, end;
    TopExp::Vertices(TopoDS::Edge(it.Key()), start, end);
    for (const auto& [vertex, t] : {std::pair(start, made.first), std::pair(end, made.last)})
      if (NewVertex* v = plan->vertices.ChangeSeek(vertex))
        v->tolerance = std::max({v->tolerance,
                                 made.tolerance,
                                 1.001 * v->point.Distance(made.curve->Value(t))});
  }

  BRepTools_Modifier modifier(body, plan);
  if (!modifier.IsDone())
  {
    why = "the faces could not be rebuilt on their offset surfaces";
    return {};
  }
  TopoDS_Shape result = modifier.ModifiedShape(body);
  BRep_Builder builder;
  for (NCollection_DataMap<TopoDS_Shape, NewEdge, TopTools_ShapeMapHasher>::Iterator it(plan->edges);
       it.More();
       it.Next())
  {
    const TopoDS_Edge edge = TopoDS::Edge(modifier.ModifiedShape(it.Key()));
    builder.Range(edge, it.Value().first, it.Value().last);
    builder.SameRange(edge, true);
    builder.SameParameter(edge, true);
  }
  // Faces whose offset surface would not load again: a spline over the face's own box.
  for (TopExp_Explorer all(body, TopAbs_FACE); all.More(); all.Next())
    if (plan->faces.Find(all.Current()).approximate
        && !asSpline(result,
                     TopoDS::Face(modifier.ModifiedShape(all.Current())),
                     plan->faces.Find(all.Current()).surface,
                     0.1 * tolerance))
    {
      why = "the offset of a face joined at a crease could not be made a spline within the "
            "tolerance";
      return {};
    }
  for (TopExp_Explorer all(body, TopAbs_FACE); all.More(); all.Next())
    lineage.set(all.Current(), {modifier.ModifiedShape(all.Current())});
  lineage.within(result);

  // The certificate.
  if (!BRepCheck_Analyzer(result).IsValid())
  {
    why = "the faces moved with their topology kept do not make a valid shape";
    return {};
  }
  // Within the tolerance, or a result edge's own where it carries an input's (above).
  const auto atDistance = [&](const TopoDS_Face& source,
                              const gp_Pnt&      p,
                              double             u,
                              double             v,
                              double             carried = 0.0) {
    const double within = std::max(tolerance, carried);
    const Moved& m = plan->faces.Find(source);
    const Handle(Geom_Surface) surface = BRep_Tool::Surface(source);
    gp_Pnt       at;
    if (!nearest(surface, p, u, v, at))
    {
      return false;
    }
    gp_Vec su, sv;
    gp_Pnt s;
    surface->D1(u, v, s, su, sv);
    gp_Vec normal = su.Crossed(sv);
    if (source.Orientation() == TopAbs_REVERSED)
      normal.Reverse();
    // At a pole S_u x S_v vanishes against the derivatives' own size and its direction is
    // rounding (GeomLProp's limit there can point either way): the distance alone is checked.
    // Judgement: a millionth of |S_u|^2 + |S_v|^2 counts as vanished.
    const bool   pole    = normal.Magnitude() <= 1e-6 * (su.SquareMagnitude() + sv.SquareMagnitude());
    const double signed_ = pole ? std::copysign(at.Distance(p), m.by)
                                : gp_Vec(at, p).Dot(normal.Normalized());
    return std::abs(signed_ - m.by) <= within && std::abs(at.Distance(p) - std::abs(m.by)) <= within;
  };
  for (TopExp_Explorer all(body, TopAbs_FACE); all.More(); all.Next())
  {
    const TopoDS_Face source = TopoDS::Face(all.Current());
    const TopoDS_Face made   = TopoDS::Face(modifier.ModifiedShape(source));
    // Points inside the face on a 5 x 5 grid of its parameter box: a point (u, v) of the moved
    // surface, from whose (u, v) the foot on the source is found.
    const BRepAdaptor_Surface  box(made, false);
    BRepTopAdaptor_FClass2d    inside(made, Precision::PConfusion());
    const Handle(Geom_Surface) movedSurface = BRep_Tool::Surface(made);
    const double u0 = box.FirstUParameter(), u1 = box.LastUParameter(),
                 v0 = box.FirstVParameter(), v1 = box.LastVParameter();
    for (int i = 0; i < 5; ++i)
      for (int j = 0; j < 5; ++j)
      {
        const double u = u0 + (u1 - u0) * (i + 0.5) / 5, v = v0 + (v1 - v0) * (j + 0.5) / 5;
        if (inside.Perform(gp_Pnt2d(u, v)) != TopAbs_IN)
          continue;
        if (!atDistance(source, movedSurface->Value(u, v), u, v))
        {
          why = "a face of the result is not at the distance from its source";
          return {};
        }
      }
  }
  for (int i = 1; i <= facesOfEdge.Extent(); ++i)
  {
    const TopoDS_Edge& edge = TopoDS::Edge(facesOfEdge.FindKey(i));
    if (BRep_Tool::Degenerated(edge))
      continue;
    const TopoDS_Edge made = TopoDS::Edge(modifier.ModifiedShape(edge));
    double            first, last;
    const Handle(Geom_Curve) curve = BRep_Tool::Curve(made, first, last);
    for (const TopoDS_Shape& f : facesOfEdge(i))
    {
      const TopoDS_Face source = TopoDS::Face(f);
      double            c0, c1;
      const Handle(Geom2d_Curve) onFace =
        BRep_Tool::CurveOnSurface(made, TopoDS::Face(modifier.ModifiedShape(source)), c0, c1);
      for (int k = 0; k <= 16 && !onFace.IsNull(); ++k)
      {
        const double   t  = first + (last - first) * k / 16;
        const gp_Pnt2d uv = onFace->Value(t);
        if (!atDistance(source,
                        curve->Value(t),
                        uv.X(),
                        uv.Y(),
                        keeps.Contains(edge) ? BRep_Tool::Tolerance(made) : 0.0))
        {
          why = "an edge of the result is not at the distance from both its faces' sources";
          return {};
        }
      }
    }
  }
  if (crossesItself(result, 0.01 * tolerance))
  {
    why = "the faces moved with their topology kept cross each other";
    return {};
  }
  if (body.ShapeType() == TopAbs_SOLID)
  {
    const double before = estimatedVolume(body), after = estimatedVolume(result);
    if (!(after > 0.0) || (after - before) * distance <= 0.0)
    {
      why = "the faces moved with their topology kept enclose the wrong side";
      return {};
    }
  }
  bool every = true;
  for (TopExp_Explorer all(body, TopAbs_FACE); all.More() && every; all.Next())
    every = moved.Contains(all.Current());
  if (every)
  {
    std::vector<TopoDS_Shape> skin;
    for (TopExp_Explorer all(result, TopAbs_FACE); all.More(); all.Next())
      skin.push_back(all.Current());
    if (clearance(skin, body) < std::abs(distance) - clearanceSlack * tolerance)
    {
      why = "the offset comes nearer to the body than the distance";
      return {};
    }
  }
  return result;
}

} // namespace BRepExt
