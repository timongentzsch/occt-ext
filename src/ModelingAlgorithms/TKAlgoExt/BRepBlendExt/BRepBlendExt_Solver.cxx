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

// Blend (BRepBlendExt_Internal.hxx): a request's pieces, their joints and splits, and the blends
// attached.

#include <GeomAlgoExt_Cancel.hxx>
#include <BRepBlendExt_Internal.hxx>
#include <BRepAlgoExt_Cells.hxx>
#include <BRepAlgoExt_FaceBlend.hxx>
#include <BRepAlgoExt_Prism.hxx>
#include <BRepAlgoExt_FillHoles.hxx>
#include <BRepAlgoExt_RoundFace.hxx>
#include <BRepAlgoAPI_Defeaturing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <BRepAlgoAPI_Splitter.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <ElCLib.hxx>
#include <Extrema_ExtPC.hxx>

namespace BRepExt::blend2
{

using namespace detail;

static Blended blendStage(const Request& request, bool continuation);

static Blended blendCopy(const Request&                   request,
                         BRepBuilderAPI_Copy&             copy,
                         const std::vector<TopoDS_Shape>& byCentre = {});

//! How the faces beside an edge turn at its middle: -1 where the edge is convex, 1 where it is
//! concave, 0 where they are tangent for a blend of `size` (tolerance::tangentFor). A face lies to
//! the left of its boundary seen along its outward normal, so with d the edge's direction in the
//! first face the edge is convex where (n0 x n1) . d > 0.
static int convexity(const char*        op,
                     const TopoDS_Edge& edge,
                     const TopoDS_Face& first,
                     const TopoDS_Face& second,
                     double             size)
{
  const BRepAdaptor_Curve curve(edge);
  const double            middle = (curve.FirstParameter() + curve.LastParameter()) / 2;
  gp_Pnt                  at;
  gp_Vec                  d;
  curve.D1(middle, at, d);
  gp_Vec            normals[2];
  const TopoDS_Face faces[2] = {first, second};
  for (int side = 0; side < 2; ++side)
  {
    double         f, l;
    const gp_Pnt2d uv = BRep_Tool::CurveOnSurface(edge, faces[side], f, l)->Value(middle);
    const Side     s{faces[side],
                 BRep_Tool::Surface(faces[side]),
                 faces[side].Orientation() == TopAbs_REVERSED ? -1.0 : 1.0,
                 0};
    Jet            j;
    if (!jetAt(s, uv.X(), uv.Y(), j))
      fail(op, invalid_input, "a face beside the edge is singular there", {edge});
    normals[side] = j.n;
  }
  if (SectionPolicy::tangentFor(normals[0].Angle(normals[1]), size))
    return 0;
  TopAbs_Orientation inFirst = TopAbs_FORWARD;
  for (TopExp_Explorer e(first, TopAbs_EDGE); e.More(); e.Next())
    if (e.Current().IsSame(edge))
      inFirst = e.Current().Orientation();
  return normals[0].Crossed(normals[1]).Dot(d) * (inFirst == TopAbs_REVERSED ? -1 : 1) > 0 ? -1 : 1;
}

//! How far an edge reaches from the origin: the largest coordinate of its box (Rail::magnitude).
static double magnitudeOf(const TopoDS_Edge& edge)
{
  Bnd_Box box;
  BRepBndLib::Add(edge, box);
  if (box.IsVoid())
    return 0;
  double x0, y0, z0, x1, y1, z1;
  box.Get(x0, y0, z0, x1, y1, z1);
  return std::max(
    {std::abs(x0), std::abs(y0), std::abs(z0), std::abs(x1), std::abs(y1), std::abs(z1)});
}

//! The side as a cliff on `edge`: the ball rests on the edge, its angle about it measured from
//! the axis most across the edge.
static bool restOn(Side& side, const TopoDS_Edge& edge)
{
  double     f, l;
  const auto curve = BRep_Tool::Curve(edge, f, l);
  if (curve.IsNull())
    return false;
  gp_Pnt p;
  gp_Vec d;
  curve->D1((f + l) / 2, p, d);
  const double x = std::abs(d.X()), y = std::abs(d.Y()), h = std::abs(d.Z());
  side.cliff          = edge;
  side.cliffCurve     = curve;
  side.cliffReference = x <= y && x <= h ? gp_Vec(1, 0, 0)
                        : y <= h         ? gp_Vec(0, 1, 0)
                                         : gp_Vec(0, 0, 1);
  return true;
}

//! The edge of `face` first met going across it from `edge` at its parameter `t`, along the face's
//! parameters at right angles to the edge's curve on it, into the face: empty where none is met.
static TopoDS_Edge firstAcross(const TopoDS_Face& face,
                               const TopoDS_Edge& edge,
                               double             t,
                               Inside&            inside)
{
  double     f, l;
  const auto onFace = BRep_Tool::CurveOnSurface(edge, face, f, l);
  if (onFace.IsNull())
    return {};
  gp_Pnt2d at;
  gp_Vec2d d;
  onFace->D1(t, at, d);
  if (!(d.Magnitude() > gp::Resolution()))
    return {};
  gp_Dir2d normal(-d.Y(), d.X());
  double   u0, u1, v0, v1;
  BRepTools::UVBounds(face, u0, u1, v0, v1);
  const double step = 1e-6 * std::max(u1 - u0, v1 - v0);
  if (inside.out(at.Translated(gp_Vec2d(normal) * step)))
    normal.Reverse();
  const Handle(Geom2d_Line) ray = new Geom2d_Line(at, normal);
  TopoDS_Edge               first;
  double                    nearest = Precision::Infinite();
  for (TopExp_Explorer e(face, TopAbs_EDGE); e.More(); e.Next())
  {
    const TopoDS_Edge other = TopoDS::Edge(e.Current());
    if (other.IsSame(edge))
      continue;
    double     a, b;
    const auto curve = BRep_Tool::CurveOnSurface(other, face, a, b);
    if (curve.IsNull())
      continue;
    Geom2dAPI_InterCurveCurve meet(
      new Geom2d_TrimmedCurve(ray, 0, 4 * std::hypot(u1 - u0, v1 - v0)),
      new Geom2d_TrimmedCurve(curve, a, b),
      tol::parameter);
    for (int i = 1; i <= meet.NbPoints(); ++i)
    {
      const double away = at.Distance(meet.Point(i));
      if (away > step && away < nearest)
      {
        nearest = away;
        first   = other;
      }
    }
  }
  return first;
}

//! The far edge of a face the blend's contact runs off along the whole edge: the edge first met
//! across it from the blended one, the same along it; whether the face is that strip alone (four
//! edges, the far one joined to each end of the blended one by one), and whether the far edge ends
//! with the strip at one end of the blend at least (the face's edge from the blended edge's vertex
//! there reaches it). `inside` classifies points of the face (Blend::inside keeps one per face).
struct Far
{
  TopoDS_Edge edge;
  bool        whole = false, ends = false;
};

static std::optional<Far> farEdge(const Rail& rail, int side, Inside& inside)
{
  const TopoDS_Face& face = rail.s[side].face;
  Far                out;
  double             f, l;
  BRep_Tool::Range(rail.edge, f, l);
  for (const double share : {0.1, 0.5, 0.9})
  {
    const TopoDS_Edge across = firstAcross(face, rail.edge, f + (l - f) * share, inside);
    if (across.IsNull() || (!out.edge.IsNull() && !out.edge.IsSame(across)))
      return {};
    out.edge = across;
  }
  if (BRep_Tool::Degenerated(out.edge))
    return {};
  std::vector<TopoDS_Edge> edges;
  int                      wires = 0;
  for (TopoDS_Iterator w(face); w.More(); w.Next(), ++wires)
    for (TopoDS_Iterator e(w.Value()); e.More(); e.Next())
      edges.push_back(TopoDS::Edge(e.Value()));
  TopoDS_Vertex a, z, f0, f1;
  TopExp::Vertices(rail.edge, a, z);
  TopExp::Vertices(out.edge, f0, f1);
  out.whole = wires == 1 && edges.size() == 4;
  for (const auto& e : edges)
  {
    if (e.IsSame(rail.edge) || e.IsSame(out.edge))
      continue;
    TopoDS_Vertex p, q;
    TopExp::Vertices(e, p, q);
    const bool atRail = p.IsSame(a) || p.IsSame(z) || q.IsSame(a) || q.IsSame(z);
    const bool atFar  = p.IsSame(f0) || p.IsSame(f1) || q.IsSame(f0) || q.IsSame(f1);
    out.whole         = out.whole && atRail;
    out.ends          = out.ends || (atRail && atFar);
  }
  return out;
}

//! Where a blend's contact runs off its face `side` along the whole edge, the face is narrower
//! than the blend: a strip between the edge and its far edge. The blend replaces the strip; where
//! the strip goes on past both ends, it keeps the round and trims it by the face beyond (a notch).
//! The Round only, which is the ball's own section: no law is tangent to an edge.
static bool cliffOf(Rail& rail, int side, const Far& far)
{
  if (rail.law != Law::Conic || std::abs(rail.rho - 0.5) > limit::equal
      || std::abs(rail.s[0].distance - rail.s[1].distance) > limit::equal * rail.scale || !far.ends)
    return false;
  rail.s[side].cliffWhole = far.whole;
  return restOn(rail.s[side], far.edge);
}

//! The edges asked, with those tangent to them added: a blend goes on to the next edge where it
//! meets the previous one smoothly (tangent for a blend of `size`, tolerance::tangentFor). An edge
//! that goes on straight between
//! other faces is a sharp continuation, not added. Where `mayAdd` is not empty, only edges in it
//! are added.
static std::vector<TopoDS_Shape> withTangent(const TopoDS_Shape&              body,
                                             std::vector<TopoDS_Shape>        edges,
                                             double                           size,
                                             const std::vector<TopoDS_Shape>& mayAdd = {})
{
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    edgesOfVertex, facesOfEdge;
  TopExp::MapShapesAndUniqueAncestors(body, TopAbs_VERTEX, TopAbs_EDGE, edgesOfVertex);
  TopExp::MapShapesAndUniqueAncestors(body, TopAbs_EDGE, TopAbs_FACE, facesOfEdge);
  const auto has = [](const std::vector<TopoDS_Shape>& list, const TopoDS_Shape& shape) {
    return std::any_of(list.begin(), list.end(), [&](const TopoDS_Shape& e) {
      return e.IsSame(shape);
    });
  };
  // The normals of the faces beside an edge at a vertex of it.
  const auto normals =
    [&](const TopoDS_Shape& edge, const TopoDS_Vertex& vertex, std::vector<gp_Vec>& out) {
      if (!facesOfEdge.Contains(edge))
        return false;
      for (const auto& face : facesOfEdge.FindFromKey(edge))
      {
        gp_Vec n;
        if (!normalAt(TopoDS::Face(face), TopoDS::Edge(edge), vertex, n))
          return false;
        out.push_back(n);
      }
      return true;
    };
  const auto smooth =
    [&](const TopoDS_Shape& edge, const TopoDS_Shape& other, const TopoDS_Vertex& vertex) {
      std::vector<gp_Vec> mine, theirs;
      if (!normals(edge, vertex, mine) || !normals(other, vertex, theirs))
        return false;
      return std::all_of(theirs.begin(), theirs.end(), [&](const gp_Vec& n) {
        return std::any_of(mine.begin(), mine.end(), [&](const gp_Vec& m) {
          return m.Dot(n) > 0 && SectionPolicy::tangentFor(m.Angle(n), size);
        });
      });
    };
  for (size_t i = 0; i < edges.size(); ++i)
    for (TopExp_Explorer v(edges[i], TopAbs_VERTEX); v.More(); v.Next())
    {
      const TopoDS_Vertex vertex = TopoDS::Vertex(v.Current());
      gp_Vec              mine, theirs;
      if (!edgesOfVertex.Contains(vertex) || !leaving(TopoDS::Edge(edges[i]), vertex, mine))
        continue;
      for (const auto& other : edgesOfVertex.FindFromKey(vertex))
        if (!has(edges, other) && (mayAdd.empty() || has(mayAdd, other))
            && leaving(TopoDS::Edge(other), vertex, theirs) && tangent(mine, theirs)
            && smooth(edges[i], other, vertex))
          edges.push_back(other);
    }
  return edges;
}

//! Convex and concave edges asked together: those of the minority convexity at the vertices of
//! mixed convexity go first (an edge is of the majority where most of the edges at the vertex,
//! blended or not, are of its kind), so that the later balls roll on over the earlier blends: a
//! convex ball over a concave blend (a torus about its axis), a boss's concave ball round the
//! convex round. The
//! concave edges go first where the vertices tie or disagree. Empty where the edges asked are all
//! of one kind. `opposite` tries the other order after the preferred one refuses; each stage
//! keeps its ordinary geometry and validity checks.
static std::optional<Blended> minorityFirst(const Request& request,
                                            bool           opposite,
                                            bool           continuation)
{
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    facesOfEdge, edgesOfVertex;
  TopExp::MapShapesAndUniqueAncestors(request.body, TopAbs_EDGE, TopAbs_FACE, facesOfEdge);
  TopExp::MapShapesAndUniqueAncestors(request.body, TopAbs_VERTEX, TopAbs_EDGE, edgesOfVertex);
  const auto kindOf = [&](const TopoDS_Shape& shape) {
    if (!facesOfEdge.Contains(shape) || facesOfEdge.FindFromKey(shape).Extent() != 2)
      return 0;
    const auto& faces = facesOfEdge.FindFromKey(shape);
    return convexity(request.operation,
                     TopoDS::Edge(shape),
                     TopoDS::Face(faces.First()),
                     TopoDS::Face(faces.Last()),
                     std::max(request.first, request.second));
  };
  std::vector<TopoDS_Shape> convex, concave;
  for (const auto& shape : request.edges)
  {
    const int kind = kindOf(shape);
    if (kind < 0)
      convex.push_back(shape);
    if (kind > 0)
      concave.push_back(shape);
  }
  if (convex.empty() || concave.empty())
    return {};
  // The minority at each vertex where an asked convex edge meets an asked concave one.
  int convexFirst = 0, concaveFirstVotes = 0; // votes of the vertices
  for (const auto& a : convex)
    for (TopExp_Explorer v(a, TopAbs_VERTEX); v.More(); v.Next())
    {
      const auto& around = edgesOfVertex.FindFromKey(v.Current());
      const bool  mixed  = std::any_of(concave.begin(), concave.end(), [&](const TopoDS_Shape& b) {
        return std::any_of(around.begin(), around.end(), [&](const TopoDS_Shape& e) {
          return e.IsSame(b);
        });
      });
      if (!mixed)
        continue;
      int convexAt = 0, concaveAt = 0;
      for (const auto& e : around)
      {
        const int kind = kindOf(e);
        convexAt += kind < 0;
        concaveAt += kind > 0;
      }
      if (convexAt < concaveAt)
        ++convexFirst;
      else
        ++concaveFirstVotes; // the concave minority, or a tie
    }
  const bool convexGoFirst         = (convexFirst > 0 && concaveFirstVotes == 0) != opposite;
  std::vector<TopoDS_Shape>& early = convexGoFirst ? convex : concave;
  std::vector<TopoDS_Shape>& late  = convexGoFirst ? concave : convex;
  Request                    first = request;
  first.edges                      = early;
  const Blended made               = blendStage(first, continuation);
  const auto    image              = [&](const TopoDS_Shape& shape) {
    std::vector<TopoDS_Shape> out;
    for (const auto& next : made.history->Modified(shape))
      out.push_back(next);
    if (out.empty() && !made.history->IsRemoved(shape))
      out.push_back(shape);
    return out;
  };
  Request second = request;
  second.body    = made.shape;
  second.edges.clear();
  for (const auto& shape : late)
  {
    // An asked edge the first blends took whole (running on along it where it goes on straight
    // from theirs) would be dropped unasked: that order does not make what was asked.
    const auto next = image(shape);
    if (next.empty())
      fail(request.operation,
           unsupported_geometry,
           "the blends of one kind use up an asked edge of the other kind",
           {shape});
    second.edges.insert(second.edges.end(), next.begin(), next.end());
  }
  // The later blends run on over what the earlier ones made, and along edges tangent to them in
  // the body as it was; not along edges that the earlier blends left tangent to them.
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> old;
  for (const auto& shape : withTangent(request.body, late, std::max(request.first, request.second)))
    for (const auto& next : image(shape))
      old.Add(next);
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> images;
  for (TopExp_Explorer e(request.body, TopAbs_EDGE); e.More(); e.Next())
    for (const auto& next : image(e.Current()))
      images.Add(next);
  for (TopExp_Explorer e(made.shape, TopAbs_EDGE); e.More(); e.Next())
    if (old.Contains(e.Current()) || !images.Contains(e.Current()))
      second.mayAdd.push_back(e.Current());
  if (!request.reference.IsNull())
  {
    const auto faces = image(request.reference);
    second.reference = faces.empty() ? TopoDS_Shape() : faces.front();
  }
  Blended then = blendStage(second, continuation);
  Blended out;
  out.shape   = then.shape;
  out.history = made.history;
  out.history->Merge(then.history);
  for (const auto& [before, middle] : made.unchanged)
    for (const auto& [again, after] : then.unchanged)
      if (again.IsSame(middle))
        out.unchanged.emplace_back(before, after);
  out.rails      = made.rails + then.rails;
  out.propagated = made.propagated + then.propagated;
  out.analytic   = made.analytic + then.analytic;
  for (const auto& [before, after] : out.unchanged)
    (before.ShapeType() == TopAbs_FACE ? out.keptFaces : out.keptEdges)++;
  return out;
}

//! Whether a ball of radius `r` beside `edge` folds on `cap`: somewhere the edge bends towards
//! the cap's inside, read with the cap's outward normal along the edge as the cap's loop runs it,
//! by 1/r or more (its geodesic curvature), so that the contact there, its offset by r, folds back.
static bool foldsOn(const TopoDS_Edge& edge, const TopoDS_Face& cap, double r)
{
  for (TopExp_Explorer e(cap, TopAbs_EDGE); e.More(); e.Next())
  {
    if (!e.Current().IsSame(edge))
      continue;
    const TopoDS_Edge       inCap = TopoDS::Edge(e.Current());
    const BRepAdaptor_Curve curve(inCap);
    const BoundaryJet       jet("filletV2", inCap, cap);
    const double            sign    = inCap.Orientation() == TopAbs_REVERSED ? -1 : 1;
    constexpr int           samples = 16;
    for (int i = 0; i <= samples; ++i)
    {
      const double t =
        curve.FirstParameter() + (curve.LastParameter() - curve.FirstParameter()) * i / samples;
      gp_Pnt at;
      gp_Vec d1, d2;
      curve.D2(t, at, d1, d2);
      if (!(d1.Magnitude() > limit::degenerate))
        continue;
      const gp_Vec along  = d1.Normalized() * sign;
      const gp_Vec inside = jet.at(t).n.Crossed(along);
      const gp_Vec bend   = (d2 - d1.Normalized() * d2.Dot(d1.Normalized())) / d1.SquareMagnitude();
      if (bend.Dot(inside) * r >= 1 - limit::agree)
        return true;
    }
  }
  return false;
}

//! An existing round that the blends asked cannot roll over, because its radius is the smaller
//! (the contact on the face beside it folds, or no ball rests on it): the round is taken off, the
//! edges asked are blended on the sharp body, and the round is put back at its own radius along
//! them (an edge rounded 1 under a top rounded 3 ends in a rolling ball of 1 between the two 3s).
//! `folded`
//! are the edges that could not be blended; empty where no round beside them is smaller.
//! What one stage made of these shapes: each modified, kept or removed.
static std::vector<TopoDS_Shape> onward(const Handle(BRepTools_History)& stage,
                                        const std::vector<TopoDS_Shape>& from)
{
  std::vector<TopoDS_Shape> to;
  for (const auto& shape : from)
  {
    if (stage->IsRemoved(shape))
      continue;
    if (stage->Modified(shape).IsEmpty())
      to.push_back(shape);
    for (const auto& next : stage->Modified(shape))
      to.push_back(next);
  }
  return to;
}

//! What became of each shape of `body` through the histories one after the other: its images
//! carried through every stage (modified, kept or removed), and what each stage generated from
//! them carried through the stages after.
static Handle(BRepTools_History) composed(const TopoDS_Shape&                           body,
                                          const std::vector<Handle(BRepTools_History)>& stages)
{
  Handle(BRepTools_History)                                     out = new BRepTools_History;
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> parts;
  TopExp::MapShapes(body, parts);
  for (int i = 1; i <= parts.Extent(); ++i)
  {
    const TopoDS_Shape& part = parts(i);
    if (!BRepTools_History::IsSupportedType(part))
      continue;
    std::vector<TopoDS_Shape> images = {part}, generated;
    for (const auto& stage : stages)
    {
      std::vector<TopoDS_Shape> made;
      for (const auto& image : images)
        for (const auto& g : stage->Generated(image))
          made.push_back(g);
      generated = onward(stage, generated);
      generated.insert(generated.end(), made.begin(), made.end());
      images = onward(stage, images);
    }
    if (images.empty())
      out->Remove(part);
    for (const auto& image : images)
      if (!image.IsSame(part))
        out->AddModified(part, image);
    for (const auto& g : generated)
      out->AddGenerated(part, g);
  }
  return out;
}

static std::optional<Blended> replayOverrun(const Request&                   request,
                                            const std::vector<TopoDS_Shape>& folded)
{
  const char*  op    = request.operation;
  const double reach = std::max(request.first, request.second);
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    facesOfEdge;
  TopExp::MapShapesAndUniqueAncestors(request.body, TopAbs_EDGE, TopAbs_FACE, facesOfEdge);
  std::vector<TopoDS_Face> rounds;
  std::vector<RoundFace>   recognised;
  for (int i = 1; i <= facesOfEdge.Extent(); ++i)
  {
    if (std::none_of(folded.begin(), folded.end(), [&](const TopoDS_Shape& e) {
          return e.IsSame(facesOfEdge.FindKey(i));
        }))
      continue;
    const auto& faces = facesOfEdge.FindFromIndex(i);
    if (faces.Extent() != 2)
      continue;
    for (const auto& face : faces)
    {
      const auto  round  = recognizeRound(request.body, TopoDS::Face(face), tol::length);
      const auto& other  = faces.First().IsSame(face) ? faces.Last() : faces.First();
      const bool  across = std::none_of(round.supports.begin(),
                                       round.supports.end(),
                                       [&](const TopoDS_Face& s) { return s.IsSame(other); });
      if (!round.valid || !(round.radius() < reach) || !across
          || !foldsOn(TopoDS::Edge(facesOfEdge.FindKey(i)), TopoDS::Face(other), reach))
        continue;
      if (std::none_of(rounds.begin(), rounds.end(), [&](const TopoDS_Face& f) {
            return f.IsSame(face);
          }))
      {
        rounds.push_back(TopoDS::Face(face));
        recognised.push_back(round);
      }
    }
  }
  if (rounds.empty())
    return {};
  // Defeaturing adds curves to the edges it is given: it is given a copy, so that the input is
  // left as it is; the copy's history goes first.
  BRepBuilderAPI_Copy       copy(request.body, true, false);
  Handle(BRepTools_History) lifted = new BRepTools_History;
  NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>
    original; // copy -> input
  {
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> parts;
    TopExp::MapShapes(request.body, parts);
    for (int i = 1; i <= parts.Extent(); ++i)
      if (BRepTools_History::IsSupportedType(parts(i)))
      {
        lifted->AddModified(parts(i), copy.ModifiedShape(parts(i)));
        original.Bind(copy.ModifiedShape(parts(i)), parts(i));
      }
  }
  BRepAlgoAPI_Defeaturing lift;
  lift.SetShape(copy.Shape());
  for (const auto& face : rounds)
    lift.AddFaceToRemove(copy.ModifiedShape(face));
  lift.SetRunParallel(false);
  lift.Build(progress(true)->Start());
  checkCancelled(op);
  if (!lift.IsDone() || lift.Shape().IsNull())
    return {};
  const Handle(BRepTools_History) defeatured = lift.History();
  lifted->Merge(defeatured);
  const auto through = [](const Handle(BRepTools_History)& history, const TopoDS_Shape& shape) {
    std::vector<TopoDS_Shape> out;
    for (const auto& next : history->Modified(shape))
      out.push_back(next);
    if (out.empty() && !history->IsRemoved(shape))
      out.push_back(shape);
    return out;
  };
  Request first = request;
  first.body    = lift.Shape();
  first.replays = request.replays + 1;
  // The chain as asked runs on through the rounds' arcs; the rounds taken off leave sharp
  // corners that tangency no longer crosses, so the sharp body is given the whole chain.
  // The chain's edges along a round taken off go with it.
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> alongRounds;
  for (const auto& face : rounds)
    TopExp::MapShapes(face, TopAbs_EDGE, alongRounds);
  first.edges.clear();
  for (const auto& edge : withTangent(request.body, request.edges, reach, request.mayAdd))
    if (!alongRounds.Contains(edge))
      for (const auto& next : through(lifted, edge))
        first.edges.push_back(next);
  first.mayAdd.clear();
  for (const auto& edge : request.mayAdd)
    for (const auto& next : through(lifted, edge))
      first.mayAdd.push_back(next);
  if (!request.reference.IsNull())
  {
    const auto faces = through(lifted, request.reference);
    first.reference  = faces.empty() ? TopoDS_Shape() : faces.front();
  }
  if (first.edges.empty())
    return {};
  // Where the replay cannot be made, the fold is what is reported, on the input's own shapes.
  try
  {
    const Blended made = blend(first);
    // Each round's edge on the sharp body: the one its two supports' images now share.
    NCollection_IndexedDataMap<TopoDS_Shape,
                               NCollection_List<TopoDS_Shape>,
                               TopTools_ShapeMapHasher>
      sharpFaces;
    TopExp::MapShapesAndUniqueAncestors(first.body, TopAbs_EDGE, TopAbs_FACE, sharpFaces);
    Blended out;
    // Each stage's history, the input's shapes carried through them all at the end
    // (composed here: BRepTools_History::Merge loses a face that one stage modifies and
    // another generates from).
    std::vector<Handle(BRepTools_History)> stages = {lifted, made.history};
    // From the sharp body on, through each stage: where its edges are now.
    Handle(BRepTools_History) fromSharp = new BRepTools_History;
    fromSharp->Merge(made.history);
    TopoDS_Shape shape = made.shape;
    int          rails = made.rails, propagated = made.propagated, analytic = made.analytic;

    // A round taken off, the faces put back, and the stages after them.
    struct Back
    {
      TopoDS_Face               round;
      std::vector<TopoDS_Shape> faces;
      size_t                    after;
    };

    std::vector<Back>                                  putBack;
    std::vector<std::pair<TopoDS_Shape, TopoDS_Shape>> kept;
    for (const auto& [before, after] : made.unchanged)
      if (const auto* input = original.Seek(before);
          input && defeatured->Modified(before).IsEmpty() && !defeatured->IsRemoved(before))
        kept.emplace_back(*input, after);
    for (size_t r = 0; r < rounds.size(); ++r)
    {
      std::vector<TopoDS_Shape> a, b, sharp;
      for (const auto& s : through(lifted, recognised[r].supports[0]))
        a.push_back(s);
      for (const auto& s : through(lifted, recognised[r].supports[1]))
        b.push_back(s);
      const auto holds = [](const std::vector<TopoDS_Shape>& in, const TopoDS_Shape& s) {
        return std::any_of(in.begin(), in.end(), [&](const TopoDS_Shape& x) {
          return x.IsSame(s);
        });
      };
      for (int i = 1; i <= sharpFaces.Extent(); ++i)
      {
        const auto& faces = sharpFaces.FindFromIndex(i);
        if (faces.Extent() == 2
            && ((holds(a, faces.First()) && holds(b, faces.Last()))
                || (holds(b, faces.First()) && holds(a, faces.Last()))))
          sharp.push_back(sharpFaces.FindKey(i));
      }
      if (sharp.empty())
        fail(op,
             unsupported_geometry,
             "a round the blend runs over left no edge between its supports to put it back on",
             {rounds[r]});
      Request again;
      again.operation = op;
      again.replays   = request.replays + 1;
      again.body      = shape;
      again.first     = recognised[r].radius();
      again.law       = Law::Conic;
      again.rho       = 0.5;
      for (const auto& edge : sharp)
        for (const auto& next : through(fromSharp, edge))
          again.edges.push_back(next);
      const Blended then = blend(again);
      // The round taken off goes on as the faces put back along its edge.
      std::vector<TopoDS_Shape> back;
      for (const auto& edge : again.edges)
        for (const auto& g : then.history->Generated(edge))
          if (g.ShapeType() == TopAbs_FACE)
            back.push_back(g);
      stages.push_back(then.history);
      fromSharp->Merge(then.history);
      putBack.push_back({rounds[r], back, stages.size()});
      std::vector<std::pair<TopoDS_Shape, TopoDS_Shape>> still;
      for (const auto& [before, middle] : kept)
        for (const auto& [put, after] : then.unchanged)
          if (put.IsSame(middle))
            still.emplace_back(before, after);
      kept  = std::move(still);
      shape = then.shape;
      rails += then.rails;
      propagated += then.propagated;
      analytic += then.analytic;
    }
    out.history = composed(request.body, stages);
    // Taken off, then put back: the round goes on as the faces put back along its edge.
    for (const auto& back : putBack)
    {
      auto faces = back.faces;
      for (size_t k = back.after; k < stages.size(); ++k)
        faces = onward(stages[k], faces);
      if (faces.empty())
        continue;
      out.history->ReplaceModified(back.round, faces.front());
      for (size_t i = 1; i < faces.size(); ++i)
        out.history->AddModified(back.round, faces[i]);
    }
    out.shape      = shape;
    out.rails      = rails;
    out.propagated = propagated;
    out.analytic   = analytic;
    for (const auto& [before, after] : kept)
    {
      out.unchanged.emplace_back(before, after);
      (before.ShapeType() == TopAbs_FACE ? out.keptFaces : out.keptEdges)++;
    }
    return out;
  }
  catch (const Failure& failure)
  {
    if (failure.reason == std::string("cancelled"))
      throw;
    return {};
  }
}

//! A feature taken off a copy of the body for the blends: the body without it, the copy's shapes
//! carried there (the feature's faces kept, as they are put back), the feature's solids to put
//! back, its faces in the copy, and (where the solids were cut from the bodies) how the copy's
//! feature faces went into them.
struct Lifted
{
  TopoDS_Shape              shape;
  Handle(BRepTools_History) history = new BRepTools_History;
  std::vector<TopoDS_Shape> tools, faces;
  Handle(BRepTools_History) intoTools;
  bool                      continuedBoss = false;
};

//! The copy's feature faces whose features meet the rest only in hole loops of their neighbours:
//! those loops filled (fillHoles), and each feature's solid its faces and the loops capped on
//! their neighbours' surfaces (holeCap), used against them. Its faces face into it for a hole or a
//! pocket, out of it for a boss; `kind` -1 asks for the first (a convex blend), 1 the second.
//! A prismatic boss may also trim a convex blend with its walls continued inside the body.
static std::optional<Lifted> liftLoops(const TopoDS_Shape&              body,
                                       const std::vector<TopoDS_Shape>& lift,
                                       int                              kind,
                                       const char*                      op)
{
  const Filled filled = fillHoles(body, lift);
  if (!filled.rest.empty())
    return {};
  Lifted out;
  out.shape = filled.shape;
  out.faces = lift;
  for (TopExp_Explorer f(body, TopAbs_FACE); f.More(); f.Next())
    if (std::none_of(lift.begin(), lift.end(), [&](const TopoDS_Shape& l) {
          return l.IsSame(f.Current());
        }))
      for (const auto& now : filled.history->Modified(f.Current()))
        out.history->AddModified(f.Current(), now);
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    facesOfEdge;
  TopExp::MapShapesAndUniqueAncestors(body, TopAbs_EDGE, TopAbs_FACE, facesOfEdge);
  BRep_Builder      builder;
  std::vector<bool> done(lift.size(), false);
  int               polarity = 0;
  for (size_t start = 0; start < lift.size(); ++start)
  {
    if (done[start])
      continue;
    TopoDS_Shell shell;
    builder.MakeShell(shell);
    NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> capped;
    std::vector<TopoDS_Shape>                              component;
    std::vector<TopoDS_Face>                               caps;
    std::vector<TopoDS_Face>                               walls;
    std::deque<size_t>                                     next = {start};
    done[start]                                                 = true;
    while (!next.empty())
    {
      const TopoDS_Shape face = lift[next.front()];
      next.pop_front();
      builder.Add(shell, face);
      component.push_back(face);
      for (TopExp_Explorer e(face, TopAbs_EDGE); e.More(); e.Next())
      {
        if (BRep_Tool::Degenerated(TopoDS::Edge(e.Current())))
          continue;
        for (const auto& other : facesOfEdge.FindFromKey(e.Current()))
        {
          bool member = false;
          for (size_t j = 0; j < lift.size(); ++j)
            if (lift[j].IsSame(other))
            {
              member = true;
              if (!done[j])
              {
                done[j] = true;
                next.push_back(j);
              }
            }
          if (member)
            continue;
          const TopoDS_Face host = TopoDS::Face(other.Oriented(TopAbs_FORWARD));
          for (TopoDS_Iterator w(host); w.More(); w.Next())
          {
            bool holds = false;
            for (TopoDS_Iterator x(w.Value()); x.More() && !holds; x.Next())
              holds = x.Value().IsSame(e.Current());
            if (!holds || !capped.Add(w.Value()))
              continue;
            const auto cap =
              TopoDS::Face(holeCap(host, w.Value())
                             .Oriented(other.Orientation() == TopAbs_REVERSED ? TopAbs_FORWARD
                                                                              : TopAbs_REVERSED));
            builder.Add(shell, cap);
            caps.push_back(cap);
            for (TopExp_Explorer x(cap, TopAbs_EDGE); x.More(); x.Next())
              for (const auto& wall : facesOfEdge.FindFromKey(x.Current()))
                if (std::any_of(lift.begin(), lift.end(), [&](const TopoDS_Shape& l) {
                      return l.IsSame(wall);
                    }))
                  walls.push_back(TopoDS::Face(wall));
          }
        }
      }
    }
    shell.Closed(true);
    TopoDS_Solid solid;
    builder.MakeSolid(solid);
    builder.Add(solid, shell);
    GProp_GProps mass;
    BRepGProp::VolumeProperties(solid, mass);
    const int sign = mass.Mass() < 0 ? -1 : 1;
    if (polarity != 0 && polarity != sign)
      return {};
    polarity = sign;
    if (sign != kind)
    {
      // A boss meeting a convex round trims it with the continuation of
      // its walls below the host, rather than its original flat base.
      // Certify a prismatic attachment before extending its footprint.
      if (kind >= 0 || mass.Mass() <= 0 || caps.size() != 1 || walls.empty())
        return {};
      const BRepAdaptor_Surface base(caps.front());
      if (base.GetType() != GeomAbs_Plane)
        return {};
      gp_Vec inward(base.Plane().Axis().Direction());
      if (caps.front().Orientation() == TopAbs_REVERSED)
        inward.Reverse();
      for (const auto& wall : walls)
      {
        const BRepAdaptor_Surface s(wall);
        if (s.GetType() == GeomAbs_Plane)
        {
          if (std::abs(inward.Dot(gp_Vec(s.Plane().Axis().Direction()))) > tol::angle)
            return {};
        }
        else if (s.GetType() == GeomAbs_Cylinder)
        {
          if (inward.Crossed(gp_Vec(s.Cylinder().Axis().Direction())).Magnitude() > tol::angle)
            return {};
        }
        else if (s.GetType() == GeomAbs_SurfaceOfExtrusion)
        {
          if (inward.Crossed(gp_Vec(s.Direction())).Magnitude() > tol::angle)
            return {};
        }
        else
          return {};
      }
      Bnd_Box bounds;
      BRepBndLib::Add(body, bounds);
      BRepPrimAPI_MakePrism prism(caps.front(), inward * std::sqrt(bounds.SquareExtent()), true);
      if (!prism.IsDone())
        return {};
      BooleanRequest join;
      join.operation      = op;
      join.select         = Select::Union;
      join.sources        = {solid};
      join.tools          = {prism.Shape()};
      join.journal        = true;
      join.cancellable    = true;
      const auto extended = booleanCells(join);
      // Extension may restore material only inside the original solid.
      join.select  = Select::Intersection;
      join.sources = {extended.shape};
      join.tools   = {body};
      auto clipped = booleanCells(join);
      simplifyBooleanResult(clipped, join, tol::angle);
      if (singleSolid(clipped.shape).IsNull())
        return {};
      solid                = TopoDS::Solid(singleSolid(clipped.shape));
      out.continuedBoss    = true;
      const auto continued = composed(body, {extended.history, clipped.history});
      if (out.intoTools.IsNull())
        out.intoTools = new BRepTools_History;
      NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> parts;
      for (const auto& face : component)
        TopExp::MapShapes(face, parts);
      for (int i = 1; i <= parts.Extent(); ++i)
      {
        const auto& part = parts(i);
        if (!BRepTools_History::IsSupportedType(part))
          continue;
        if (continued->IsRemoved(part))
          out.intoTools->Remove(part);
        for (const auto& now : continued->Modified(part))
          out.intoTools->AddModified(part, now);
      }
    }
    // A hole's faces turned to face out of it: each face, as classifiers do not read a solid
    // used reversed.
    if (mass.Mass() < 0)
    {
      TopoDS_Shell turned;
      builder.MakeShell(turned);
      for (TopoDS_Iterator f(shell); f.More(); f.Next())
        builder.Add(turned, f.Value().Reversed());
      turned.Closed(true);
      builder.MakeSolid(solid);
      builder.Add(solid, turned);
    }
    out.tools.push_back(solid);
  }
  return out;
}

//! The copy's faces `lift` removed as Delete Face removes them where no hole loop bounds them: the
//! faces beside them extended over the gap (BRepAlgoAPI_Defeaturing, patched to extend each only
//! across its edges with them). What that filled, the body without them less the body, is the
//! hole put back from a convex blend (`kind` -1); what it took, the body less the body without
//! them, the boss put back to a concave one. Empty where the faces cannot be removed or the
//! feature is of the other kind.
static std::optional<Lifted> liftFaces(const char*                      op,
                                       const TopoDS_Shape&              body,
                                       const std::vector<TopoDS_Shape>& lift,
                                       int                              kind)
{
  BRepAlgoAPI_Defeaturing remove;
  remove.SetShape(body);
  for (const auto& face : lift)
    remove.AddFaceToRemove(face);
  remove.SetRunParallel(false);
  remove.SetToFillHistory(true);
  remove.Build(progress(true)->Start());
  checkCancelled(op);
  if (!remove.IsDone() || remove.HasErrors() || remove.Shape().IsNull())
    return {};
  for (const auto& face : lift)
    if (!remove.IsDeleted(face))
      return {};
  Lifted out;
  out.shape = remove.Shape();
  out.faces = lift;
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> parts;
  TopExp::MapShapes(body, parts);
  for (int i = 1; i <= parts.Extent(); ++i)
  {
    if (!BRepTools_History::IsSupportedType(parts(i))
        || std::any_of(lift.begin(), lift.end(), [&](const TopoDS_Shape& l) {
             return l.IsSame(parts(i));
           }))
      continue;
    if (remove.IsDeleted(parts(i)))
      out.history->Remove(parts(i));
    for (const auto& now : remove.Modified(parts(i)))
      out.history->AddModified(parts(i), now);
  }
  BooleanRequest cut;
  cut.operation               = op;
  cut.select                  = Select::Difference;
  cut.sources                 = {kind < 0 ? out.shape : body};
  cut.tools                   = {kind < 0 ? body : out.shape};
  cut.journal                 = true;
  cut.cancellable             = true;
  const BooleanResult feature = booleanCells(cut);
  // The other difference must hold nothing: a feature of one kind only.
  cut.sources.swap(cut.tools);
  const BooleanResult other = booleanCells(cut);
  if (feature.shape.IsNull() || !TopExp_Explorer(feature.shape, TopAbs_SOLID).More()
      || (!other.shape.IsNull() && TopExp_Explorer(other.shape, TopAbs_SOLID).More()))
    return {};
  out.tools     = {feature.shape};
  out.intoTools = feature.history;
  return out;
}

//! A blend that runs into faces of the body beside none of its edges (`reached`, the input's), or
//! rolls on past a face its edge ends on to the next piece of the same edge (a hole drilled
//! through it, narrower than the blend): the feature is lifted off, the edges asked blended on the
//! body without it, and the feature put back with the kernel's Boolean, its faces trimming the
//! blend: a hole or a pocket cut from a convex blend, a boss joined to a concave one. A feature
//! that meets the rest only in hole loops of its neighbours is lifted by filling the loops

//! This makes the rounds of a top edge running into a hole and a pocket, and a round over a hole
//! drilled through its edge, to 4e-15 mm. Empty otherwise, and for a
//! feature of the other kind unless its walls can be continued: a prismatic boss meeting a
//! convex blend is extended inside the original body. Other extensions remain unsupported.
static std::optional<Blended> liftFeatures(const Request&                   request,
                                           const std::vector<TopoDS_Shape>& reached)
{
  const char*  op   = request.operation;
  const double size = std::max(request.first, request.second);
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    facesOfEdge;
  TopExp::MapShapesAndUniqueAncestors(request.body, TopAbs_EDGE, TopAbs_FACE, facesOfEdge);
  // Every edge asked convex, or every one concave.
  int kind = 0;
  for (const auto& edge : request.edges)
  {
    if (!facesOfEdge.Contains(edge) || facesOfEdge.FindFromKey(edge).Extent() != 2)
      return {};
    const auto& faces = facesOfEdge.FindFromKey(edge);
    const int   k     = convexity(op,
                            TopoDS::Edge(edge),
                            TopoDS::Face(faces.First()),
                            TopoDS::Face(faces.Last()),
                            size);
    if (k == 0 || (kind != 0 && k != kind))
      return {};
    kind = k;
  }
  // No face beside a blended edge may be lifted: the blend rests on it.
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> supports;
  for (const auto& edge : withTangent(request.body, request.edges, size, request.mayAdd))
    for (const auto& face : facesOfEdge.FindFromKey(edge))
      supports.Add(face);
  // Nor any face of another shell: a cavity the blend reaches would be opened by it.
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> beside;
  for (TopExp_Explorer shell(request.body, TopAbs_SHELL); shell.More(); shell.Next())
  {
    bool holds = false;
    for (TopExp_Explorer f(shell.Current(), TopAbs_FACE); f.More() && !holds; f.Next())
      holds = supports.Contains(f.Current());
    if (holds)
      for (TopExp_Explorer f(shell.Current(), TopAbs_FACE); f.More(); f.Next())
        beside.Add(f.Current());
  }
  const auto clear = [&](const std::vector<TopoDS_Shape>& faces) {
    return std::none_of(faces.begin(), faces.end(), [&](const TopoDS_Shape& f) {
      return supports.Contains(f) || !beside.Contains(f);
    });
  };
  // On a copy, as the blends are made: what is put back must not be the input's.
  BRepBuilderAPI_Copy copy(request.body, true, false);
  // The faces in the copy, as its shells use them.
  const auto inCopy = [&](const std::vector<TopoDS_Shape>& faces) {
    std::vector<TopoDS_Shape> out;
    for (TopExp_Explorer f(copy.Shape(), TopAbs_FACE); f.More(); f.Next())
      if (std::any_of(faces.begin(), faces.end(), [&](const TopoDS_Shape& face) {
            return copy.ModifiedShape(face).IsSame(f.Current());
          }))
        out.push_back(f.Current());
    return out;
  };
  std::optional<Lifted>     lifted;
  std::vector<TopoDS_Shape> features = featuresOf(request.body, reached);
  if (clear(features))
    lifted = liftLoops(copy.Shape(), inCopy(features), kind, op);
  if (!lifted && clear(reached))
  {
    features = reached;
    lifted   = liftFaces(op, copy.Shape(), inCopy(features), kind);
  }
  // With the faces joined to those by edges, short of the supports: a pocket open across the edge
  // whose back wall the blend does not reach, without which its other walls cannot be removed.
  if (!lifted && clear(reached))
  {
    features = reached;
    for (size_t i = 0; i < features.size(); ++i)
      for (TopExp_Explorer e(features[i], TopAbs_EDGE); e.More(); e.Next())
        for (const auto& face : facesOfEdge.FindFromKey(e.Current()))
          if (!supports.Contains(face)
              && std::none_of(features.begin(), features.end(), [&](const TopoDS_Shape& f) {
                   return f.IsSame(face);
                 }))
            features.push_back(face);
    if (features.size() > reached.size() && clear(features))
      lifted = liftFaces(op, copy.Shape(), inCopy(features), kind);
  }
  if (!lifted)
    return {};
  // The input's shapes into the copy, then the feature lifted off.
  Handle(BRepTools_History) lifting = new BRepTools_History;
  NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>
    original; // copy -> input
  {
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> parts;
    TopExp::MapShapes(request.body, parts);
    for (int i = 1; i <= parts.Extent(); ++i)
      if (BRepTools_History::IsSupportedType(parts(i)))
      {
        lifting->AddModified(parts(i), copy.ModifiedShape(parts(i)));
        original.Bind(copy.ModifiedShape(parts(i)), parts(i));
      }
  }
  lifting->Merge(lifted->history);
  const auto through = [&](const TopoDS_Shape& shape) { return onward(lifting, {shape}); };
  Request    first   = request;
  first.body         = lifted->shape;
  first.replays      = request.replays + 1;
  first.edges.clear();
  first.mayAdd.clear();
  // An edge the lift joined (the pieces of an edge a hole cut) is asked once.
  const auto add = [](std::vector<TopoDS_Shape>& to, const TopoDS_Shape& shape) {
    if (std::none_of(to.begin(), to.end(), [&](const TopoDS_Shape& s) { return s.IsSame(shape); }))
      to.push_back(shape);
  };
  for (const auto& edge : request.edges)
    for (const auto& next : through(edge))
      add(first.edges, next);
  for (const auto& edge : request.mayAdd)
    for (const auto& next : through(edge))
      add(first.mayAdd, next);
  if (!request.reference.IsNull())
  {
    const auto faces = through(request.reference);
    first.reference  = faces.empty() ? TopoDS_Shape() : faces.front();
  }
  try
  {
    const Blended  made = blend(first);
    BooleanRequest back;
    back.operation   = op;
    back.select      = kind < 0 && !lifted->continuedBoss ? Select::Difference : Select::Union;
    back.sources     = {made.shape};
    back.tools       = lifted->tools;
    back.journal     = true; // its history, for the faces' history through it
    back.cancellable = true;
    const BooleanResult put = booleanCells(back);
    if (put.shape.IsNull() || singleSolid(put.shape).IsNull())
      return {};
    Blended out;
    out.shape   = singleSolid(put.shape);
    out.history = composed(request.body, {lifting, made.history, put.history});
    // Where the feature's solid was cut from the bodies, its faces go on through that cut.
    if (!lifted->intoTools.IsNull())
      for (const auto& face : features)
      {
        const auto images =
          onward(put.history, onward(lifted->intoTools, {copy.ModifiedShape(face)}));
        if (images.empty())
        {
          if (!out.history->IsRemoved(face))
            out.history->Remove(face);
          continue;
        }
        out.history->ReplaceModified(face, images.front());
        for (size_t i = 1; i < images.size(); ++i)
          out.history->AddModified(face, images[i]);
      }
    const auto kept = [&](const TopoDS_Shape& shape) {
      return !put.history->IsRemoved(shape) && put.history->Modified(shape).IsEmpty();
    };
    for (const auto& [before, after] : made.unchanged)
      if (const auto* input = original.Seek(before); input && kept(after))
        out.unchanged.emplace_back(*input, after);
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> result;
    TopExp::MapShapes(out.shape, TopAbs_FACE, result);
    for (const auto& face : features)
      if (result.Contains(copy.ModifiedShape(face)) && kept(copy.ModifiedShape(face)))
        out.unchanged.emplace_back(face, copy.ModifiedShape(face));
    for (const auto& [before, after] : out.unchanged)
      (before.ShapeType() == TopAbs_FACE ? out.keptFaces : out.keptEdges)++;
    out.rails      = made.rails;
    out.propagated = made.propagated;
    out.analytic   = made.analytic;
    return out;
  }
  catch (const Failure& failure)
  {
    if (failure.reason == std::string("cancelled"))
      throw;
    return {};
  }
}

// A coplanar support seam is a topological partition, not a blend stop.
// Retry on a private copy with only such seams removed; preserve every asked
// edge and carry the result and display history back through the normalization.
static std::optional<Blended> joinedPlanarSupports(const Request& request)
{
  BRepBuilderAPI_Copy          copy(request.body, true, false);
  ShapeUpgrade_UnifySameDomain unify(copy.Shape(), false, true, false);
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    faces;
  TopExp::MapShapesAndUniqueAncestors(copy.Shape(), TopAbs_EDGE, TopAbs_FACE, faces);
  bool candidate = false;
  for (int i = 1; i <= faces.Extent(); ++i)
  {
    const auto& pair   = faces.FindFromIndex(i);
    bool        planar = pair.Extent() == 2;
    for (const auto& face : pair)
      planar = planar && BRepAdaptor_Surface(TopoDS::Face(face)).GetType() == GeomAbs_Plane;
    if (planar)
      candidate = true;
    else
      unify.KeepShape(faces.FindKey(i));
  }
  if (!candidate)
    return {};
  for (const auto* selected : {&request.edges, &request.mayAdd})
    for (const auto& edge : *selected)
      unify.KeepShape(copy.ModifiedShape(edge));
  unify.SetLinearTolerance(tol::length);
  unify.SetAngularTolerance(tol::angle);
  unify.Build();
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> before, after;
  TopExp::MapShapes(copy.Shape(), TopAbs_FACE, before);
  TopExp::MapShapes(unify.Shape(), TopAbs_FACE, after);
  if (after.Extent() >= before.Extent() || !BRepCheck_Analyzer(unify.Shape()).IsValid())
    return {};
  Handle(BRepTools_History)                                     copied = new BRepTools_History;
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> parts;
  TopExp::MapShapes(request.body, parts);
  for (int i = 1; i <= parts.Extent(); ++i)
    if (BRepTools_History::IsSupportedType(parts(i)))
      copied->AddModified(parts(i), copy.ModifiedShape(parts(i)));
  const auto normalization = composed(request.body, {copied, unify.History()});
  Request    next          = request;
  next.body                = unify.Shape();
  next.replays++;
  next.edges  = onward(normalization, request.edges);
  next.mayAdd = onward(normalization, request.mayAdd);
  if (!request.reference.IsNull())
  {
    const auto mapped = onward(normalization, {request.reference});
    if (mapped.size() != 1)
      return {};
    next.reference = mapped.front();
  }
  try
  {
    auto made    = blend(next);
    made.history = composed(request.body, {normalization, made.history});
    std::vector<std::pair<TopoDS_Shape, TopoDS_Shape>> kept;
    for (int i = 1; i <= parts.Extent(); ++i)
    {
      const auto& part = parts(i);
      if (!BRepTools_History::IsSupportedType(part))
        continue;
      const auto copiedPart = copy.ModifiedShape(part);
      if (!unify.History()->Modified(copiedPart).IsEmpty()
          || unify.History()->IsRemoved(copiedPart))
        continue;
      for (const auto& middle : onward(normalization, {part}))
        for (const auto& [old, now] : made.unchanged)
          if (old.IsSame(middle))
            kept.emplace_back(part, now);
    }
    made.unchanged = std::move(kept);
    made.keptFaces = made.keptEdges = 0;
    for (const auto& [old, now] : made.unchanged)
      (old.ShapeType() == TopAbs_FACE ? made.keptFaces : made.keptEdges)++;
    return made;
  }
  catch (const Failure& failure)
  {
    if (failure.reason == std::string("cancelled"))
      throw;
    return {};
  }
}

//! Rounds of the two edges of a face that is only the strip between them, of half its width: their
//! contacts meet on the strip, and the full round is made as two faces that meet where they touch
//! it. Empty where the request is not that.
static std::optional<Blended> meetingRounds(const Request& request)
{
  if (request.law != Law::Conic || std::abs(request.rho - 0.5) > limit::equal
      || request.edges.size() != 2
      || (request.second > 0
          && std::abs(request.second - request.first) > limit::equal * request.first))
    return {};
  Request full = request;
  full.law     = Law::Full;
  full.meet    = request.first;
  try
  {
    return blend(full);
  }
  catch (const Failure& failure)
  {
    if (failure.reason == std::string("cancelled"))
      throw;
    return {};
  }
}

// A disconnected planar cap has medial-axis creases between its islands. Use
// the existing exact distance envelope instead of smoothing those folds below
// that envelope. Connected caps retain their ordinary rolling-ball construction.
static std::optional<Blended> islandEnvelope(const Request& request)
{
  if ((request.law != Law::Chamfer
       && (request.law != Law::Conic || std::abs(request.rho - .5) > limit::equal))
      || (request.second > 0
          && std::abs(request.second - request.first) > limit::equal * request.first))
    return {};
  TopoDS_Shape cap;
  for (TopExp_Explorer f(request.body, TopAbs_FACE); f.More(); f.Next())
  {
    if (BRepAdaptor_Surface(TopoDS::Face(f.Current())).GetType() != GeomAbs_Plane)
      continue;
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> rim;
    TopExp::MapShapes(f.Current(), TopAbs_EDGE, rim);
    if (rim.Extent() == int(request.edges.size())
        && std::all_of(request.edges.begin(), request.edges.end(), [&](const auto& edge) {
             return rim.Contains(edge);
           }))
    {
      cap = f.Current();
      break;
    }
  }
  if (cap.IsNull() || !blendExceedsConvexRadius(request.body, request.edges, request.first))
    return {};
  BRepBuilderAPI_Copy       copy(request.body, true, false);
  std::vector<TopoDS_Shape> edges;
  for (const auto& edge : request.edges)
    edges.push_back(copy.ModifiedShape(edge));
  try
  {
    const auto made =
      faceBlend(request.operation, copy.Shape(), edges, request.first, request.law == Law::Chamfer);
    if (made.Modified(copy.ModifiedShape(cap)).Extent() < 2)
      return {};
    Blended out;
    out.shape = made.shape;
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> inputs, final;
    TopExp::MapShapes(request.body, inputs);
    TopExp::MapShapes(out.shape, final);
    for (int i = 1; i <= inputs.Extent(); ++i)
    {
      const auto& before = inputs(i);
      if (!BRepTools_History::IsSupportedType(before))
        continue;
      const auto copied   = copy.ModifiedShape(before);
      bool       retained = false;
      if (before.ShapeType() == TopAbs_FACE)
      {
        if (const auto* successors = made.faces.Seek(copied))
        {
          for (const auto& after : *successors)
            if (final.Contains(after))
            {
              out.history->AddModified(before, after);
              retained = true;
            }
        }
        else if (final.Contains(copied))
        {
          out.history->AddModified(before, copied);
          retained = true;
        }
      }
      else if (final.Contains(copied))
      {
        out.history->AddModified(before, copied);
        retained = true;
      }
      if (!retained)
        out.history->Remove(before);
    }
    // This constructor changes shared boundary vertices while closing its
    // sheets. Keep display reuse conservative until those copies are proved.
    out.rails = int(request.edges.size());
    return out;
  }
  catch (const Failure& failure)
  {
    if (failure.reason == std::string("cancelled"))
      throw;
    return {};
  }
  catch (const Standard_Failure&)
  {
    return {};
  }
}

Blended blend(const Request& request)
{
  return blendStage(request, true);
}

static Blended blendStage(const Request& request, bool continuation)
{
  const char* op = request.operation;
  if ((request.law != Law::Full && !(request.first > 0)) || request.second < 0)
    fail(op, invalid_input, "the distances must be positive");
  if (request.law != Law::Chamfer && !(request.rho > 0 && request.rho < 1))
    fail(op, invalid_input, "rho is between 0 and 1");
  if (request.angle != 0
      && (request.law != Law::Chamfer || request.reference.IsNull()
          || !(request.angle > 0 && request.angle < M_PI / 2)))
    fail(op, invalid_input, "a chamfer's angle lies off a reference face, below a right angle");
  if (request.law == Law::Chordal && request.second > 0
      && std::abs(request.second - request.first) > limit::equal * request.first)
    fail(op, invalid_input, "a chordal blend has one chord");
  // A size that changes along the edges (Request::radii) is made by the blend of its own edges
  // alone: the other routes are for one size.
  const bool varies = std::any_of(request.radii.begin(),
                                  request.radii.end(),
                                  [](const auto& points) { return !points.empty(); });
  if (varies)
  {
    if (request.radii.size() > request.edges.size())
      fail(op, invalid_input, "a size law for an edge not asked");
    if (request.law == Law::Full || request.law == Law::Chordal)
      fail(op, invalid_input, "a full or chordal round takes its size from the faces");
    if (request.law != Law::Chamfer && request.second > 0
        && std::abs(request.second - request.first) > limit::equal * request.first)
      fail(op,
           unsupported_geometry,
           "a blend whose size changes along the edge has one distance off both faces");
  }
  Request own = request;
  if (varies)
    own.mayAdd = request.edges;
  if (!varies)
    if (auto islands = islandEnvelope(request))
      return *islands;
  // The blends are made on a copy without its meshes, so that the result shares nothing with the
  // input (whose meshes and edges' samples belong to its own neighbours); what they leave alone
  // is the copy's.
  if (request.law != Law::Full && !varies)
  {
    std::exception_ptr original;
    // Preserve the established two construction orders before permitting
    // newly required support continuations. Every attempt keeps its gates.
    const std::array<std::pair<bool, bool>, 4> attempts = {
      {{false, false}, {true, false}, {false, true}, {true, true}}};
    for (const auto& [opposite, complete] : attempts)
    {
      if (complete && !continuation)
        break;
      try
      {
        if (auto both = minorityFirst(request, opposite, complete))
          return *both;
        break; // This is not a mixed-convexity request.
      }
      catch (const Failure& failure)
      {
        if (failure.reason != no_solution && failure.reason != unsupported_geometry)
          throw;
        if (!original)
          original = std::current_exception();
      }
      catch (const Standard_Failure&)
      {
        if (!original)
          throw;
      }
    }
    if (original)
      std::rethrow_exception(original);
  }
  BRepBuilderAPI_Copy copy(request.body, true, false);
  try
  {
    return blendCopy(own, copy);
  }
  catch (Failure& failure)
  {
    // The faults are the copy's sub-shapes: name the input's, by shape.
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> parts;
    TopExp::MapShapes(request.body, parts);
    NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher> input; // copy -> input
    for (int i = 1; i <= parts.Extent(); ++i)
      input.Bind(copy.ModifiedShape(parts(i)), parts(i));
    for (size_t i = 0; i < failure.shapes.size(); ++i)
      if (const TopoDS_Shape* found = input.Seek(failure.shapes[i]))
      {
        failure.shapes[i] = *found;
        if (i < failure.at.size())
          failure.at[i] = failureHash(*found);
      }
    if (varies)
      throw;
    if (continuation && dynamic_cast<const ContinuationRequired*>(&failure)
        && failure.shapes.size() == 2 && !request.mayAdd.empty()
        && request.replays < limit::replays)
    {
      const auto& next = failure.shapes[1];
      if (std::none_of(request.mayAdd.begin(), request.mayAdd.end(), [&](const TopoDS_Shape& edge) {
            return edge.IsSame(next);
          }))
      {
        Request continued = request;
        continued.mayAdd.push_back(next);
        ++continued.replays;
        return blendStage(continued, true);
      }
    }
    // A blend that runs into a hole beyond its faces: the hole lifted off and put back.
    if (const auto* collision = dynamic_cast<const Collision*>(&failure);
        collision && request.replays < limit::replays)
    {
      std::vector<TopoDS_Shape> reached;
      for (TopExp_Explorer f(request.body, TopAbs_FACE); f.More(); f.Next())
        if (std::any_of(collision->faces.begin(),
                        collision->faces.end(),
                        [&](const TopoDS_Shape& face) {
                          return face.IsSame(copy.ModifiedShape(f.Current()));
                        }))
          reached.push_back(f.Current());
      if (auto lifted = liftFeatures(request, reached))
        return *lifted;
    }
    // A rail that cannot be made beside a smaller round (its contact folds, or no ball rests
    // on the round): the round is put back after the blend instead. Overflows and material
    // bounds are what they say whatever the faces are.
    if (failure.kind != FailureKind::Overflow && failure.kind != FailureKind::MaterialBound
        && failure.reason != std::string("cancelled") && request.law != Law::Full)
    {
      if (request.replays < limit::replays)
        if (auto joined = joinedPlanarSupports(request))
          return *joined;
      if (auto met = meetingRounds(request))
        return *met;
      if (request.replays < limit::replays)
        if (auto replayed = replayOverrun(request, failure.shapes))
          return *replayed;
    }
    throw;
  }
}

//! A request's blends, stage by stage: its pieces, how their ends meet, their sheets and splits,
//! and the blends attached.
struct Blend
{
  Blend(const Request&            request_,
        BRepBuilderAPI_Copy&      copy_,
        std::vector<TopoDS_Shape> byCentre_ = {})
      : request(request_),
        copy(copy_),
        op(request_.operation),
        body(copy_.Shape()),
        byCentre(std::move(byCentre_))
  {
    TopExp::MapShapesAndUniqueAncestors(body, TopAbs_EDGE, TopAbs_FACE, facesOfEdge);
    TopExp::MapShapesAndUniqueAncestors(body, TopAbs_VERTEX, TopAbs_EDGE, edgesOfVertex);
  }

  //! The pieces: one per edge asked and per tangent edge added, with its supports, distances and
  //! kind.
  void rails();
  //! A full round's one piece.
  void fullRail(const std::vector<TopoDS_Shape>& asked);
  //! How each piece's ends meet the others': a joint, a mitre or a corner.
  void joints();
  //! The corner where more than three blended edges meet.
  void vertexBlend(const TopoDS_Vertex& vertex, const std::vector<std::pair<int, int>>& meeting);
  //! Each rail carried on past an end that a face stops.
  void carries();
  //! Each piece's stations and sheet, split where a contact crosses a smooth edge onto the next
  //! face.
  void sheets();
  //! The loops where a contact on a cap runs back across a joint into the next piece.
  void crossFolds(const std::vector<std::pair<int, std::vector<Fold>>>& folds);
  //! One kind of sheet along each chain: exact where every piece of it can be, else fitted.
  void chains();
  //! The parts of one split edge on one sheet.
  void joinSplits();
  //! Which sides pair at each joint.
  void pairJoints();
  //! The edges of each chain whose blends end in different balls at a joint, where the balls
  //! found through the edge's points fold back (Rail::byCentre); empty where every joint agrees.
  std::vector<TopoDS_Shape> foldedJoints() const;
  //! The blends attached, with the history from the input.
  Blended result();
  void    splitMeeting(Blended& out);
  //! Which points of a face's surface lie off it: made once for each face, as a face of many
  //! holes takes long to make it for.
  Inside& inside(const TopoDS_Face& face);
  //! farEdge of a rail's side, worked out once for each edge and face: the sheets are made again
  //! for each attempt, and the far edge depends on the edge and the face alone.
  const std::optional<Far>& farOf(const Rail& rail, int side);
  //! Whether piece k's rail is one whose balls may be found through their centres instead
  //! (Rail::byCentre): a round's or a G2 blend's, off faces (not a cliff), not yet so.
  bool mayCentre(int k) const;
  //! The balls of piece k's chain (its tangent edges and split parts, which share sections at
  //! their joints) found through their centres from now on: sheets() is made again.
  void centreChain(int k);

  const Request&       request;
  BRepBuilderAPI_Copy& copy;
  const char*          op;
  const TopoDS_Shape   body;
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
                            facesOfEdge, edgesOfVertex;
  std::vector<TopoDS_Shape> edges; // asked and added
  TopoDS_Shape              reference;
  std::vector<Piece>        pieces;
  //! Where a contact on a flat cap runs back but crosses itself nowhere within its piece: the
  //! places, for the loop across a joint into the next piece (crossFolds).
  std::vector<std::pair<int, std::vector<Fold>>> open;

  //! The refusal of a fold the blend finds no loop for, on piece k: its contact's fold, as
  //! checkFold says.
  [[noreturn]] void refuseFold(int k, const char* why)
  {
    Dips dips{pieces[k].rail, pieces[k].fitted};
    checkFold(dips);
    fail(op, unsupported_geometry, why, {pieces[k].rail.edge}, FailureKind::Curvature);
  }

  std::vector<Corner>                 corners;
  std::vector<std::optional<Station>> seeds; // [piece]: where a split part's stations start
  std::vector<TopoDS_Shape>
    byCentre; // edges whose balls are found through their centres (foldedJoints, centreChain)
  std::vector<TopoDS_Edge>
       ending; // a chain's last edges, wholly under the blend of the edge before (absorbEnds)
  void absorbEnds();
  NCollection_DataMap<TopoDS_Shape, int, TopTools_ShapeMapHasher> insideOf; // face -> insides
  std::deque<Inside>                                              insides;
  std::deque<std::tuple<TopoDS_Shape, TopoDS_Shape, std::optional<Far>>> fars; // edge, face
  bool again = false; // sheets() to be made again (centreChain)

  bool centred(const TopoDS_Shape& edge) const
  {
    return std::any_of(byCentre.begin(), byCentre.end(), [&](const TopoDS_Shape& e) {
      return e.IsSame(edge);
    });
  }
};

bool Blend::mayCentre(int k) const
{
  const Rail& rail = pieces[k].rail;
  return !rail.byCentre && (rail.law == Law::Conic || rail.law == Law::G2) && !rail.s[0].cliffCurve
         && !rail.s[1].cliffCurve;
}

//! A chain's last edge that a contact is off the whole way, past the edge of its face where the
//! face that stops the chain begins (a short edge, a sphere's seam or a patch's edge just before
//! the stop): no blend of its own, the one before runs on over it to that stop. It is made again
//! so, the ball being the same whichever edge it is found along.
void Blend::absorbEnds()
{
  for (const auto& edge : ending)
  {
    const auto found = std::find_if(pieces.begin(), pieces.end(), [&](const Piece& p) {
      return p.rail.edge.IsSame(edge);
    });
    if (found == pieces.end())
      continue;
    const int     b       = int(found - pieces.begin());
    const int     stopEnd = found->joint[0].meets == Meets::Stop ? 0 : 1;
    const auto    joint   = found->joint[1 - stopEnd];
    TopoDS_Vertex ends[2];
    TopExp::Vertices(TopoDS::Edge(edge.Oriented(TopAbs_FORWARD)), ends[0], ends[1]);
    Piece& a                = pieces[joint.piece];
    a.joint[joint.end]      = {};
    a.stopVertex[joint.end] = ends[stopEnd];
    a.swallowed.push_back(TopoDS::Edge(edge));
    a.swallowed.insert(a.swallowed.end(), found->swallowed.begin(), found->swallowed.end());
    // On past its vertex by the length of the edge it runs on over, and on from there as to a stop.
    gp_Pnt p;
    gp_Vec d;
    Rail&  rail = a.rail;
    rail.curve.D1(joint.end == 0 ? rail.first : rail.last, p, d);
    const double over =
      GCPnts_AbscissaPoint::Length(BRepAdaptor_Curve(TopoDS::Edge(edge))) / d.Magnitude();
    const double more = 4 * rail.scale / d.Magnitude();
    if (joint.end == 0)
    {
      rail.first -= over;
      rail.from = rail.first - more;
    }
    else
    {
      rail.last += over;
      rail.to = rail.last + more;
    }
    pieces.erase(found);
    for (auto& piece : pieces)
      for (auto& j : piece.joint)
        if (j.meets != Meets::Stop && j.piece > b)
          --j.piece;
    for (auto& corner : corners)
      for (auto& [p2, e2] : corner.ends)
        if (p2 > b)
          --p2;
  }
}

void Blend::centreChain(int k)
{
  std::vector<int>  todo = {k};
  std::vector<bool> seen(pieces.size(), false);
  seen[k] = true;
  while (!todo.empty())
  {
    const int p = todo.back();
    todo.pop_back();
    if (!centred(pieces[p].rail.edge))
      byCentre.push_back(pieces[p].rail.edge);
    for (const auto& joint : pieces[p].joint)
      if ((joint.meets == Meets::Joint || joint.meets == Meets::Split) && !seen[joint.piece])
      {
        seen[joint.piece] = true;
        todo.push_back(joint.piece);
      }
  }
  again = true;
}

Inside& Blend::inside(const TopoDS_Face& face)
{
  if (const int* at = insideOf.Seek(face))
    return insides[*at];
  insideOf.Bind(face, int(insides.size()));
  return insides.emplace_back(face);
}

const std::optional<Far>& Blend::farOf(const Rail& rail, int side)
{
  const TopoDS_Face& face = rail.s[side].face;
  for (const auto& [edge, of, far] : fars)
    if (edge.IsSame(rail.edge) && of.IsSame(face))
      return far;
  return std::get<2>(fars.emplace_back(rail.edge, face, farEdge(rail, side, inside(face))));
}

//! Refused where the fitted blend crosses itself between its contacts: the envelope of the balls
//! along a path that bends more tightly than the ball (a pipe surface c + r (cos a N + sin a B)
//! moves along the path at |c'| (1 - kappa r cos a), which changes sign there), so its sections
//! run back against the way the balls roll. Read off the sheet: S_u against the station's rolling
//! direction, at 65 x 7 places; where it turns against it, the surface has folded. The reference modeller trims
//! the blend there, which v2 does not yet.
static void checkCrossing(const Rail& rail, Stations& stations, const Fitted& fitted)
{
  if (fitted.blend.IsNull())
    return;
  // Where a contact runs back, its loop is patched or refused as a fold of the contact
  // (capFolds, checkFold): this is the crossing between contacts that both roll on.
  Dips dips{rail, fitted};
  if (runsBack(dips))
    return;
  constexpr int along = 64, across = 8;
  for (int i = 0; i <= along; ++i)
  {
    const double u = double(i) / along;
    Station      s;
    if (!stations.at(rail.t(u), s))
      continue;
    for (int j = 1; j < across; ++j)
    {
      gp_Pnt p;
      gp_Vec du, dv;
      fitted.blend->D1(u, double(j) / across, p, du, dv);
      if (du.Dot(s.along) < -limit::agree * du.Magnitude() - 1e-3 * du.Magnitude())
        fail(rail.op,
             unsupported_geometry,
             "the blend crosses itself between its contacts, where the balls' path bends more "
             "tightly than the ball, which v2 does not trim yet",
             {rail.edge},
             FailureKind::Curvature);
    }
  }
}

//! Whether a conic section of the rail turns by more than 120 degrees somewhere along its edge
//! (Rail::halves): the angle between the faces' normals at the edge's points, which the section
//! turns by. Sampled at 17 places.
static bool turnsWide(const Rail& rail)
{
  if (rail.law != Law::Conic && rail.law != Law::Chordal)
    return false;
  constexpr int samples = 16;
  for (int k = 0; k <= samples; ++k)
  {
    const double t = rail.first + (rail.last - rail.first) * k / samples;
    gp_Vec       n[2];
    bool         found = true;
    for (int i = 0; i < 2 && found; ++i)
    {
      double     f, l;
      const auto pcurve = BRep_Tool::CurveOnSurface(rail.edge, rail.s[i].face, f, l);
      Jet        j;
      found = !pcurve.IsNull() && !rail.s[i].cliffCurve
              && jetAt(rail.s[i], pcurve->Value(t).X(), pcurve->Value(t).Y(), j);
      if (found)
        n[i] = j.n;
    }
    if (found && n[0].Angle(n[1]) > 2 * M_PI / 3)
      return true;
  }
  return false;
}

//! The rail's size law (Rail::sizes) from the (place, size) pairs asked for its edge: the factor on
//! its distances by the edge's parameter, positive all along; round an edge that closes, the size
//! must end as it begins. Its scale is the largest.
static void sizesOf(Rail& rail, const std::vector<std::pair<double, double>>& points, double first)
{
  const char* op = rail.op;
  const int   n  = int(points.size());
  if (n < 2 || points.front().first != 0 || points.back().first != 1)
    fail(op,
         invalid_input,
         "a size law runs from the edge's start (0) to its end (1)",
         {rail.edge});
  for (int i = 0; i < n; ++i)
  {
    const auto [place, size] = points[i];
    if (!std::isfinite(place) || !std::isfinite(size) || !(size > 0)
        || (i > 0 && !(place > points[i - 1].first)))
      fail(op, invalid_input, "a size law takes rising places and positive sizes", {rail.edge});
    rail.sizes.places.push_back(rail.first + (rail.last - rail.first) * place);
    rail.sizes.values.push_back(size / first);
  }
  TopoDS_Vertex a, b;
  TopExp::Vertices(rail.edge, a, b);
  rail.sizes.periodic = !a.IsNull() && a.IsSame(b);
  if (rail.sizes.periodic
      && std::abs(points.front().second - points.back().second) > limit::equal * first)
    fail(op,
         invalid_input,
         "the size round a closed edge ends where it begins: "
           + std::to_string(points.front().second) + " and " + std::to_string(points.back().second)
           + " are at one point",
         {rail.edge},
         FailureKind::RadiusEndpoints);
  rail.sizes.solve();
  constexpr int samples = 64;
  double        most    = 0;
  for (int k = 0; k <= samples; ++k)
  {
    const double f = rail.factor(rail.first + (rail.last - rail.first) * k / samples);
    if (!(f > 0))
      fail(op,
           unsupported_geometry,
           "the size law falls to zero between its places",
           {rail.edge},
           FailureKind::RadiusLaw);
    most = std::max(most, f);
  }
  rail.scale *= most;
}

void Blend::rails()
{

  std::vector<TopoDS_Shape> asked, mayAdd;
  for (const auto& edge : request.edges)
    asked.push_back(copy.ModifiedShape(edge));
  for (const auto& edge : request.mayAdd)
    mayAdd.push_back(copy.ModifiedShape(edge));
  if (request.law == Law::Full)
    return fullRail(asked);
  edges     = withTangent(body, asked, std::max(request.first, request.second), mayAdd);
  reference = request.reference.IsNull() ? TopoDS_Shape() : copy.ModifiedShape(request.reference);
  pieces.reserve(edges.size() * 8);
  for (const auto& shape : edges)
  {
    const TopoDS_Edge edge = TopoDS::Edge(shape);
    if (std::any_of(pieces.begin(), pieces.end(), [&](const Piece& p) {
          return p.rail.edge.IsSame(edge);
        }))
      continue;
    // A seam has its face on both sides and an edge between tangent faces has nothing to

    if (facesOfEdge.Contains(edge) && facesOfEdge.FindFromKey(edge).Extent() == 1
        && BRep_Tool::IsClosed(edge, TopoDS::Face(facesOfEdge.FindFromKey(edge).First())))
      continue;
    if (!facesOfEdge.Contains(edge) || facesOfEdge.FindFromKey(edge).Extent() != 2)
      fail(op, invalid_input, "a blended edge must lie between two faces", {edge});
    pieces.emplace_back();
    Rail& rail = pieces.back().rail;
    rail.op    = op;
    rail.edge  = edge;
    rail.curve.Initialize(edge);
    rail.law = request.law;
    rail.rho = request.law == Law::Chordal ? 0.5 : request.rho;
    int i    = 0;
    for (const auto& face : facesOfEdge.FindFromKey(edge))
    {
      rail.s[i].face    = TopoDS::Face(face);
      rail.s[i].surface = BRep_Tool::Surface(rail.s[i].face);
      rail.s[i].outward = rail.s[i].face.Orientation() == TopAbs_REVERSED ? -1 : 1;
      ++i;
    }
    const double second        = request.second > 0 ? request.second : request.first;
    const bool   firstOnSecond = !reference.IsNull() && rail.s[1].face.IsSame(reference);
    rail.s[0].distance         = firstOnSecond ? second : request.first;
    rail.s[1].distance         = firstOnSecond ? request.first : second;
    rail.angle                 = request.angle;
    rail.reference             = firstOnSecond ? 1 : 0;
    rail.scale                 = std::max(rail.s[0].distance, rail.s[1].distance);
    rail.magnitude             = magnitudeOf(edge);
    rail.first                 = rail.curve.FirstParameter();
    rail.last                  = rail.curve.LastParameter();
    for (size_t k = 0; k < asked.size() && k < request.radii.size(); ++k)
      if (asked[k].IsSame(edge) && !request.radii[k].empty())
        sizesOf(rail, request.radii[k], request.first);
    const int kind = convexity(op, edge, rail.s[0].face, rail.s[1].face, rail.scale);
    if (kind == 0)
    {
      pieces.pop_back();
      continue;
    }
    rail.sigma  = kind;
    rail.halves = turnsWide(rail);
    TopoDS_Vertex start, end;
    TopExp::Vertices(edge, start, end);
    rail.closed   = !start.IsNull() && start.IsSame(end);
    rail.byCentre = std::any_of(byCentre.begin(), byCentre.end(), [&](const TopoDS_Shape& e) {
      return e.IsSame(edge);
    });
  }
}

void Blend::fullRail(const std::vector<TopoDS_Shape>& asked)
{
  if (asked.size() != 2)
    fail(op, invalid_input, "a full round takes the two edges of one face");
  TopoDS_Edge rails[2];
  TopoDS_Face beside[2], between;
  for (int i = 0; i < 2; ++i)
  {
    rails[i] = TopoDS::Edge(asked[i]);
    if (!facesOfEdge.Contains(rails[i]) || facesOfEdge.FindFromKey(rails[i]).Extent() != 2)
      fail(op, invalid_input, "a blended edge must lie between two faces", {rails[i]});
  }
  for (const auto& face : facesOfEdge.FindFromKey(rails[0]))
    for (const auto& other : facesOfEdge.FindFromKey(rails[1]))
      if (face.IsSame(other))
        between = TopoDS::Face(face);
  if (between.IsNull())
    fail(op, invalid_input, "a full round's two edges must bound one face", {rails[0], rails[1]});
  for (int i = 0; i < 2; ++i)
    for (const auto& face : facesOfEdge.FindFromKey(rails[i]))
      if (!face.IsSame(between))
        beside[i] = TopoDS::Face(face);
  const auto sideOf = [](const TopoDS_Face& face) {
    Side side;
    side.face    = face;
    side.surface = BRep_Tool::Surface(face);
    side.outward = face.Orientation() == TopAbs_REVERSED ? -1 : 1;
    return side;
  };
  edges = {rails[0]};
  pieces.reserve(8);
  pieces.emplace_back();
  Rail& rail = pieces.back().rail;
  rail.op    = op;
  rail.edge  = rails[0];
  rail.curve.Initialize(rails[0]);
  rail.law       = Law::Full;
  rail.magnitude = magnitudeOf(rails[0]);
  rail.s[0]      = sideOf(beside[0]);
  rail.s[1] = rail.middle = sideOf(between);
  rail.first              = rail.curve.FirstParameter();
  rail.last               = rail.curve.LastParameter();
  TopoDS_Vertex start, end;
  TopExp::Vertices(rails[0], start, end);
  rail.closed    = !start.IsNull() && start.IsSame(end);
  Inside     onStrip(rail.s[1].face);
  const auto far = farEdge(rail, 1, onStrip);
  if (!far || !far->edge.IsSame(rails[1]) || !far->whole)
    fail(op,
         unsupported_geometry,
         "a full round's face must be only the strip between its two edges",
         {rails[0], rails[1]});
  // Half the width across, which the ball's radius is near: for tolerances and the carry.
  BRepExtrema_DistShapeShape width(rails[0], rails[1]);
  if (!width.IsDone() || !(width.Value() > tol::length))
    fail(op, no_solution, "a full round's two edges meet", {rails[0], rails[1]});
  rail.scale     = width.Value() / 2;
  const int kind = convexity(op, rails[0], beside[0], between, rail.scale);
  if (kind == 0 || convexity(op, rails[1], between, beside[1], rail.scale) != kind)
    fail(op,
         unsupported_geometry,
         "a full round's two edges must both be convex or both concave",
         {rails[0], rails[1]});
  rail.sigma = kind;
  // The notch of the face between, whose face beyond the ball rests on (Rail::middle).
  rail.s[1]          = sideOf(beside[1]);
  rail.s[1].face     = between;
  rail.s[1].trim     = beside[1];
  rail.s[1].trimEdge = rails[1];
  rail.s[0].distance = rail.s[1].distance = rail.scale;
}

void Blend::joints()
{
  // Where two blended edges meet tangent to each other, one blend runs on from one into the
  // other; where more meet, or at an angle, the corner is blended (milestone 3).
  NCollection_IndexedDataMap<TopoDS_Shape,
                             std::vector<std::pair<int, int>>,
                             TopTools_ShapeMapHasher>
                                            endsAt;
  std::vector<std::array<TopoDS_Vertex, 2>> endOf(pieces.size());
  for (int k = 0; k < int(pieces.size()); ++k)
  {
    auto& ends = endOf[k];
    TopExp::Vertices(TopoDS::Edge(pieces[k].rail.edge.Oriented(TopAbs_FORWARD)), ends[0], ends[1]);
    for (int end = 0; end < 2; ++end)
    {
      if (!endsAt.Contains(ends[end]))
        endsAt.Add(ends[end], {});
      endsAt.ChangeFromKey(ends[end]).push_back({k, end});
    }
  }
  for (int v = 1; v <= endsAt.Extent(); ++v)
  {
    const auto& meeting = endsAt(v);
    if (meeting.size() == 1)
      continue;
    const auto [a, ea]          = meeting[0];
    const TopoDS_Vertex& vertex = endOf[a][ea];
    if (meeting.size() == 3)
    {
      Corner corner{vertex, meeting};
      for (const auto& [k, end] : meeting)
        pieces[k].joint[end] = {Meets::Corner, int(corners.size()), 0, false};
      corners.push_back(corner);
      continue;
    }
    if (meeting.size() > 3)
    {
      vertexBlend(vertex, meeting);
      continue;
    }
    const auto [b, eb] = meeting.back();
    gp_Vec da, db;
    if (!leaving(pieces[a].rail.edge, vertex, da, ea)
        || !leaving(pieces[b].rail.edge, vertex, db, eb))
      fail(op,
           unsupported_geometry,
           "a blended edge has no direction where it meets another",
           {pieces[a].rail.edge});
    const Meets meets = tangent(da, db) ? Meets::Joint : Meets::Miter;
    if (meets == Meets::Miter && a == b)
      fail(op,
           unsupported_geometry,
           "an edge that closes on itself at an angle is blended with its corner, which v2 does "
           "not do yet",
           {pieces[a].rail.edge});
    pieces[a].joint[ea] = {meets, b, eb, false};
    pieces[b].joint[eb] = {meets, a, ea, false};
  }
}

// Where more than three blended edges meet: the rolling ball's corner, where the body is what the
// balls inside it fill (its opening by the ball). Where one ball rests on all the faces round the
// vertex, the corner is that ball's cap about its apex. Else, between planes, the corner parts as
// the region of the balls' centres does (the faces' planes offset by the radius): each corner of
// that region is a ball resting on three faces or more, faces next to each other round the vertex,
// and two such balls that rest on the same two faces not next to each other are joined by the round
// between those faces, which meet only at the vertex: along the line where their planes meet, from
// one ball to the other. Four faces part so into two balls and a round between the faces across
// from each other, as

void Blend::vertexBlend(const TopoDS_Vertex&                    vertex,
                        const std::vector<std::pair<int, int>>& meeting)
{
  const TopoDS_Edge lead = pieces[meeting[0].first].rail.edge;
  const double      r    = pieces[meeting[0].first].rail.s[0].distance,
               sigma     = pieces[meeting[0].first].rail.sigma;
  // Balls whose centres are that near are one (limit::one); a face is met by a ball where its
  // offset is that near the centre.
  const double scale = pieces[meeting[0].first].rail.scale,
               near  = limit::one(scale, BRep_Tool::Pnt(vertex).XYZ().Modulus());
  const auto refuse  = [&](const char* why) { fail(op, unsupported_geometry, why, {lead}); };
  for (const auto& [k, end] : meeting)
  {
    const Rail& rail  = pieces[k].rail;
    const bool  round = rail.law == Law::Conic && std::abs(rail.rho - 0.5) < limit::equal
                       && std::abs(rail.s[0].distance - r) < limit::equal * scale
                       && std::abs(rail.s[1].distance - r) < limit::equal * scale;
    if (!round || rail.sigma != sigma)
      refuse("where more than three blended edges meet, v2 blends the corner of rounds of one "
             "radius, all convex or all concave");
  }
  // The faces round the vertex in turn: after[i] is the end of the blend between ring[i] and ring[i
  // + 1].
  const int                        count = int(meeting.size());
  std::vector<TopoDS_Face>         ring  = {pieces[meeting[0].first].rail.s[0].face};
  std::vector<std::pair<int, int>> after;
  std::vector<bool>                used(count, false);
  for (bool found = true; found && int(after.size()) < count;)
  {
    found = false;
    for (int m = 0; m < count && !found; ++m)
      for (int side = 0; side < 2 && !found && !used[m]; ++side)
        if (pieces[meeting[m].first].rail.s[side].face.IsSame(ring.back()))
        {
          used[m] = found = true;
          after.push_back(meeting[m]);
          ring.push_back(pieces[meeting[m].first].rail.s[1 - side].face);
        }
  }
  const char* ringless = "the faces at a vertex of more than three blended edges do not go round "
                         "it in turn, each between two of the blends";
  if (int(after.size()) != count || !ring.back().IsSame(ring.front()))
    refuse(ringless);
  ring.pop_back();
  for (int i = 0; i < count; ++i)
    for (int j = i + 1; j < count; ++j)
      if (ring[i].IsSame(ring[j]))
        refuse(ringless);
  std::vector<Side> sides(count);
  for (int i = 0; i < count; ++i)
    sides[i] = Side{ring[i],
                    BRep_Tool::Surface(ring[i]),
                    ring[i].Orientation() == TopAbs_REVERSED ? -1.0 : 1.0,
                    r};
  const auto corner = [&](const std::vector<std::pair<int, int>>& ends) {
    for (const auto& [k, end] : ends)
      pieces[k].joint[end] = {Meets::Corner, int(corners.size()), 0, false};
    corners.push_back({vertex, ends});
  };
  const gp_Pnt        V = BRep_Tool::Pnt(vertex);
  std::vector<gp_Vec> n(count);
  for (int i = 0; i < count; ++i)
    if (!normalAt(ring[i], pieces[after[i].first].rail.edge, vertex, n[i]))
      refuse("a face at a vertex of more than three blended edges has no normal there");
  // Faces tangent to each other at the vertex: the ball rests on them all there, and each blend
  // closes to the vertex (Rail::pinch).
  if (std::all_of(n.begin(), n.end(), [&](const gp_Vec& v) {
        return v.Angle(n[0]) <= tol::coarseAngle;
      }))
  {
    for (const auto& [k, end] : after)
      pieces[k].rail.pinch[end] = true;
    corner(after);
    corners.back().pinch = true;
    return;
  }
  if (restingBall(sides, V, r, sigma, scale))
    return corner(after);
  // The corners of the centres' region: where three of the offset planes meet, clear of the
  // others (on the material's side of each for a convex corner, the air's for a concave one).
  for (int i = 0; i < count; ++i)
    if (GeomAdaptor_Surface(sides[i].surface).GetType() != GeomAbs_Plane)
      refuse("where more than three blended edges meet between faces that are not all planes and "
             "no one ball rests on them all, v2 does not blend the corner yet");

  struct Cell
  {
    gp_Pnt           centre;
    std::vector<int> faces; // ascending: in turn round the vertex
  };

  std::vector<Cell> cells;
  for (int a = 0; a < count; ++a)
    for (int b = a + 1; b < count; ++b)
      for (int c = b + 1; c < count; ++c)
      {
        std::array<std::array<double, 4>, 3> m{};
        const int                            rows[3] = {a, b, c};
        for (int i = 0; i < 3; ++i)
        {
          for (int j = 0; j < 3; ++j)
            m[i][j] = n[rows[i]].Coord(j + 1);
          m[i][3] = sigma * r;
        }
        std::array<double, 3> y{};
        if (!solve<3>(m, y))
          continue;
        const gp_Vec     offset(y[0], y[1], y[2]);
        std::vector<int> tight;
        bool             clear = true;
        for (int i = 0; i < count && clear; ++i)
        {
          const double off = sigma * n[i].Dot(offset) - r;
          if (off < -near)
            clear = false;
          else if (off <= near)
            tight.push_back(i);
        }
        if (!clear)
          continue;
        const gp_Pnt centre = V.Translated(offset);
        const auto   same   = std::find_if(cells.begin(), cells.end(), [&](const Cell& cell) {
          return cell.centre.Distance(centre) <= near;
        });
        if (same == cells.end())
          cells.push_back({centre, tight});
        else
        {
          same->faces.insert(same->faces.end(), tight.begin(), tight.end());
          std::sort(same->faces.begin(), same->faces.end());
          same->faces.erase(std::unique(same->faces.begin(), same->faces.end()), same->faces.end());
        }
      }
  // Each blend's end in one corner, each pair of faces not next to each other in two or none,
  // and one corner more than such pairs: the corners part the faces' ring.
  const char* parted = "the corner where more than three blended edges meet does not part into "
                       "balls on faces next to each other";
  const auto  pairOf = [](int a, int b) { return std::make_pair(std::min(a, b), std::max(a, b)); };
  const auto  next   = [&](const std::pair<int, int>& pair) {
    return pair.second == pair.first + 1 || (pair.first == 0 && pair.second == count - 1);
  };
  std::map<std::pair<int, int>, std::vector<int>> sharedBy;
  for (int c = 0; c < int(cells.size()); ++c)
  {
    const auto& f = cells[c].faces;
    if (f.size() < 3)
      refuse(parted);
    for (size_t m = 0; m < f.size(); ++m)
      sharedBy[pairOf(f[m], f[(m + 1) % f.size()])].push_back(c);
  }
  int diagonals = 0;
  for (const auto& [pair, users] : sharedBy)
  {
    if (users.size() != (next(pair) ? 1u : 2u))
      refuse(parted);
    diagonals += !next(pair);
  }
  for (int i = 0; i < count; ++i)
    if (!sharedBy.count(pairOf(i, (i + 1) % count)))
      refuse(parted);
  if (int(cells.size()) != diagonals + 1)
    refuse(parted);
  // The round between each such pair: its rail the line where their planes meet, from the
  // place across from one ball's centre to the other's.
  std::map<std::pair<int, int>, int> roundOf;
  for (auto& [pair, users] : sharedBy)
  {
    if (next(pair))
      continue;
    const auto [i, j] = pair;
    gp_Vec d          = n[i].Crossed(n[j]);
    if (!(d.Magnitude() > tol::angle))
      refuse(parted);
    d.Normalize();
    double ta = gp_Vec(V, cells[users[0]].centre).Dot(d),
           tb = gp_Vec(V, cells[users[1]].centre).Dot(d);
    if (ta > tb)
    {
      std::swap(users[0], users[1]);
      std::swap(ta, tb);
    }
    const Handle(Geom_Line) line = new Geom_Line(V, gp_Dir(d));
    std::vector<OnSurface>  on;
    for (const int f : {i, j})
    {
      const auto [surface, location] = surfaceOf(ring[f]);
      const gp_Pln  plane            = GeomAdaptor_Surface(surface).Plane();
      const gp_Trsf back             = location.Inverted().Transformation();
      double        u0, v0, u1, v1;
      ElSLib::Parameters(plane, line->Value(ta).Transformed(back), u0, v0);
      ElSLib::Parameters(plane, line->Value(tb).Transformed(back), u1, v1);
      on.push_back({segment2d({u0, v0}, {u1, v1}, ta, tb), surface, location});
    }
    TopoDS_Vertex va, vb;
    BRep_Builder  builder;
    builder.MakeVertex(va, line->Value(ta), tol::length);
    builder.MakeVertex(vb, line->Value(tb), tol::length);
    if (pieces.size() >= pieces.capacity())
      refuse(parted);
    Piece piece;
    Rail& rail = piece.rail;
    rail.op    = op;
    rail.edge  = edgeBetween(line, va, ta, vb, tb, on, tol::length);
    rail.curve.Initialize(rail.edge);
    rail.s[0]     = sides[i];
    rail.s[1]     = sides[j];
    rail.sigma    = sigma;
    rail.scale    = scale;
    rail.first    = ta;
    rail.last     = tb;
    piece.apex    = vertex;
    roundOf[pair] = int(pieces.size());
    pieces.push_back(piece);
  }
  for (int c = 0; c < int(cells.size()); ++c)
  {
    const auto&                      f = cells[c].faces;
    std::vector<std::pair<int, int>> ends;
    for (size_t m = 0; m < f.size(); ++m)
    {
      const auto pair = pairOf(f[m], f[(m + 1) % f.size()]);
      if (next(pair))
        ends.push_back(after[pair.second == pair.first + 1 ? pair.first : count - 1]);
      else
        ends.push_back({roundOf[pair], sharedBy[pair][0] == c ? 0 : 1});
    }
    corner(ends);
  }
}

void Blend::carries()
{
  // Each piece's rail runs on past an end that a face stops, so that the blend reaches a face
  // that leans away from it; at a joint it ends with the edge.
  for (auto& piece : pieces)
  {
    Rail&      rail = piece.rail;
    const auto run  = [&](int end) {
      if (piece.joint[end].meets == Meets::Joint || piece.joint[end].meets == Meets::Corner)
        return 0.0;
      gp_Pnt p;
      gp_Vec d;
      rail.curve.D1(end == 0 ? rail.first : rail.last, p, d);
      // The blend's size at this end, where it changes along the edge (Rail::sizes).
      const double size = rail.sizes.places.empty()
                             ? rail.scale
                             : std::max(rail.s[0].distance, rail.s[1].distance)
                                * rail.factor(end == 0 ? rail.first : rail.last);
      // Twice the blend's size on, so that it reaches a face that leans away; a mitre's
      // blends may meet a leaning face further beyond the vertex: twice as far there.
      double length = (piece.joint[end].meets == Meets::Miter ? 4 : 2) * size;
      // At a notch (Attach::notch), where an edge of a support goes on straight from the
      // blended one, the blend may run on along it to the face at its far end: that far on.
      // Where it leaves at an angle to the rail's run on, the contact on that support, its
      // offset off the edge, crosses it that offset over the angle's tangent past the
      // vertex: that far on, within the edge. So too where a face stops the blend
      // (Attach::stops), whose contacts end where they cross the supports' edges there.
      // A chamfer's offset is its setback; a ball's is its radius over the half angle's
      // tangent between the faces.
      if (piece.joint[end].meets == Meets::Miter || piece.joint[end].meets == Meets::Stop)
      {
        TopoDS_Vertex ends[2];
        TopExp::Vertices(TopoDS::Edge(rail.edge.Oriented(TopAbs_FORWARD)), ends[0], ends[1]);
        gp_Vec mine, theirs, n0, n1;
        if (leaving(rail.edge, ends[end], mine, end)
            && normalAt(rail.s[0].face, rail.edge, ends[end], n0)
            && normalAt(rail.s[1].face, rail.edge, ends[end], n1))
        {
          const double ballOffset = size * std::tan(n0.Angle(n1) / 2);
          for (const auto& other : edgesOfVertex.FindFromKey(ends[end]))
          {
            if (other.IsSame(rail.edge) || !leaving(TopoDS::Edge(other), ends[end], theirs))
              continue;
            int on = -1;
            for (const auto& face : facesOfEdge.FindFromKey(other))
              for (int side = 0; side < 2; ++side)
                if (face.IsSame(rail.s[side].face))
                  on = side;
            if (on < 0)
              continue;
            const double offset =
              rail.law == Law::Chamfer ? rail.s[on].distance * size / rail.scale : ballOffset;
            const double reach =
              GCPnts_AbscissaPoint::Length(BRepAdaptor_Curve(TopoDS::Edge(other))) + 2 * size;
            const double angle = theirs.Angle(-mine);
            if (tangent(mine, theirs))
              length = std::max(length, reach);
            else if (angle < M_PI / 2)
              length = std::max(length, std::min(reach, offset / std::tan(angle) + 2 * size));
            // A mitre's partner on the same face: the contacts, offsets d and d' off the
            // two edges at an angle theta between them, cross (d' + d cos theta) / sin theta
            // past the vertex along this one, far beyond it at a narrow re-entrant notch.
            const auto& joint = piece.joint[end];
            if (joint.meets != Meets::Miter || !other.IsSame(pieces[joint.piece].rail.edge))
              continue;
            const Rail& partner = pieces[joint.piece].rail;
            gp_Vec      m0, m1;
            if (rail.law == Law::Chamfer
                || !normalAt(partner.s[0].face, partner.edge, ends[end], m0)
                || !normalAt(partner.s[1].face, partner.edge, ends[end], m1))
              continue;
            const double theta = mine.Angle(theirs);
            if (!(std::sin(theta) > tol::angle))
              continue;
            const double theirOffset = partner.scale * std::tan(m0.Angle(m1) / 2);
            length =
              std::max(length,
                       (theirOffset + offset * std::cos(theta)) / std::sin(theta) + 2 * size);
          }
        }
      }
      // Along a curved edge no further than it turns by a quarter of a turn: the stations'
      // planes across the rail turn with it, and past the end they need not meet the ball's
      // path (a split piece between two planes on an arc rolls straight on). The march ends
      // the carry where the ball strays from the rail (Stations::reached); an L corner's
      // ball rolling on round a concave blend goes a sixth of a turn to the face beyond.
      gp_Vec d2;
      rail.curve.D2(end == 0 ? rail.first : rail.last, p, d, d2);
      const double curvature = d.Crossed(d2).Magnitude() / std::pow(d.Magnitude(), 3);
      if (curvature > 0)
        length = std::min(length, M_PI / 2 / curvature);
      return length / d.Magnitude();
    };
    rail.from = rail.first - run(0);
    rail.to   = rail.last + run(1);
  }
}

//! A chordal blend whose ball keeps one size along the edge (between planes, or round a rim
//! whose wedge does not change) is that round: made as one, exact where a round is.
static void roundOfChord(Rail& rail, Stations& stations)
{
  constexpr int samples = 9;
  double        least = Precision::Infinite(), most = 0;
  for (int k = 0; k < samples; ++k)
  {
    Station s;
    if (!stations.at(rail.t(double(k) / (samples - 1)), s))
      return;
    const double r = s.centre.Distance(s.foot[0]);
    least          = std::min(least, r);
    most           = std::max(most, r);
  }
  if (most - least > limit::agree * most)
    return;
  rail.law           = Law::Conic;
  rail.s[0].distance = rail.s[1].distance = rail.scale = (least + most) / 2;
}

//! Whether a full round has the radius `r` all along (Request::meet).
static void meetsAt(const Rail& rail, Stations& stations, double r)
{
  constexpr int samples = 9;
  for (int k = 0; k < samples; ++k)
  {
    Station s;
    if (!stations.at(rail.t(double(k) / (samples - 1)), s)
        || std::abs(s.centre.Distance(s.foot[0]) - r)
             > tol::resolved(limit::together, r, rail.magnitude))
      fail(rail.op,
           no_solution,
           "the rounds' contacts on the face between do not meet along it",
           {rail.edge});
  }
}

//! Whether a contact may cross a seam of its face between its ends: a seam along a line of the
//! face's parameters, at one of them, is crossed only where the contact's poles reach past its
//! place in some period inside their range (a contact from seam to seam, round a hole, does not
//! cross it between); any other seam may be.
static bool acrossSeam(const Handle(Geom2d_BSplineCurve)& contact,
                       const TopoDS_Edge&                 seam,
                       const TopoDS_Face&                 face,
                       const Handle(Geom_Surface)&        surface)
{
  double     f, l;
  const auto pcurve = BRep_Tool::CurveOnSurface(seam, face, f, l);
  const auto line   = Handle(Geom2d_Line)::DownCast(pcurve);
  if (line.IsNull())
    return true;
  const gp_Dir2d d      = line->Direction();
  const int      across = std::abs(d.X()) < tol::angle ? 0 : std::abs(d.Y()) < tol::angle ? 1 : -1;
  if (across < 0)
    return true;
  const double period = across == 0 ? (surface->IsUPeriodic() ? surface->UPeriod() : 0)
                                    : (surface->IsVPeriodic() ? surface->VPeriod() : 0);
  if (!(period > 0))
    return true;
  double lo = HUGE_VAL, hi = -HUGE_VAL;
  for (int i = 1; i <= contact->NbPoles(); ++i)
  {
    const double x = across == 0 ? contact->Pole(i).X() : contact->Pole(i).Y();
    lo             = std::min(lo, x);
    hi             = std::max(hi, x);
  }
  const double at    = across == 0 ? line->Location().X() : line->Location().Y(),
               slack = 1e-9 * period;
  const double first = at + std::ceil((lo + slack - at) / period) * period;
  return first < hi - slack;
}

void Blend::sheets()
{
  // Sections in two spans along a whole chain, or along none of it: its sections share their
  // form at its joints.
  for (bool changed = true; changed;)
  {
    changed = false;
    for (auto& piece : pieces)
      for (const auto& joint : piece.joint)
        if ((joint.meets == Meets::Joint || joint.meets == Meets::Split)
            && pieces[joint.piece].rail.halves != piece.rail.halves)
        {
          pieces[joint.piece].rail.halves = piece.rail.halves = true;
          changed                                             = true;
        }
  }
  // Each piece's blend, exact where every piece of its chain can be; a chain's sections then
  // share one parameter at each joint.
  // Where a contact leaves its face across a smooth edge onto the next face, the piece is split
  // there and its second part rolls on that face (the support sequence).
  seeds.assign(pieces.size(), std::nullopt);
  std::vector<std::pair<int, std::vector<Fold>>> folds;
  // Whether a contact that runs back is left to the caller: taken by neither a fold nor a loop
  // across a joint.
  const auto capFoldsOf = [&](int k, Dips& dips) {
    std::vector<Fold> across;
    auto              found = capFolds(dips, &across);
    const bool        taken = !found.empty() || !across.empty();
    if (!found.empty())
      folds.emplace_back(k, std::move(found));
    if (!across.empty())
      open.emplace_back(k, std::move(across));
    return !taken;
  };
  // Where a contact on face `s2` leaves it over `edge` at `at`: whether the edge is convex there
  // (towards the blended edge, in the face, runs under the face beyond).
  const auto convexOver = [&](const Rail&        rail,
                              int                s2,
                              const TopoDS_Edge& edge,
                              const Meeting&     at,
                              const Handle(Geom2d_Curve)&,
                              double) {
    for (const auto& face : facesOfEdge.FindFromKey(edge))
    {
      if (face.IsSame(rail.s[s2].face))
        continue;
      double     f, l;
      const auto pc = BRep_Tool::CurveOnSurface(edge, TopoDS::Face(face), f, l);
      if (pc.IsNull())
        return false;
      Side beyond    = rail.s[s2];
      beyond.face    = TopoDS::Face(face);
      beyond.surface = BRep_Tool::Surface(beyond.face);
      beyond.outward = face.Orientation() == TopAbs_REVERSED ? -1 : 1;
      beyond.cliffCurve.Nullify();
      const gp_Pnt2d there = pc->Value(at.second);
      Jet            b;
      if (!jetAt(beyond, there.X(), there.Y(), b))
        return false;
      gp_Pnt inside;
      rail.curve.D0(rail.t(at.first), inside);
      return gp_Vec(b.p, inside).Dot(b.n) < 0;
    }
    return false;
  };
  for (size_t k = 0; k < pieces.size(); ++k)
  {
    if (pieces.size() >= pieces.capacity())
      fail(op,
           unsupported_geometry,
           "the blend crosses too many faces along one edge",
           {pieces[k].rail.edge});
    Piece& piece = pieces[k];
    if (centred(piece.rail.edge))
      piece.rail.byCentre = true;
    // Where the planes through the edge's points cross before they reach it, the balls found
    // through those points do not go on along it, and through their centres they do (the
    // same balls, placed by another parameter: Rail::byCentre); the piece's chain is made so
    // (centreChain). Through the edge's points the march stalls there, or near there the
    // balls' place along the edge turns back (a contact runs back) or changes too sharply
    // to fit. Only where through the centres the blend is fitted and no contact runs back:
    // else the fold or the failure is the blend's own, and stands.
    const auto centreMends = [&] {
      if (!mayCentre(int(k)))
        return false;
      Rail rail     = piece.rail;
      rail.byCentre = true;
      Stations stations(rail);
      try
      {
        stations.march(seeds[k]);
        const Fitted fitted = fit(rail, stations);
        Dips         centred{rail, fitted};
        return !runsBack(centred);
      }
      catch (const Failure& failure)
      {
        if (failure.reason == std::string("cancelled"))
          throw;
        return false;
      }
    };
    for (int attempt = 0; attempt < 2; ++attempt)
    {
      Stations stations(piece.rail);
      // No ball's plane passes through the edge's points somewhere along it (past a bend
      // tighter than the ball, where the contact on a cap folds): the balls are found through
      // their centres instead (Rail::byCentre), the same balls.
      const auto byCentre = [&] {
        const bool ball = piece.rail.law == Law::Conic || piece.rail.law == Law::G2;
        if (piece.rail.byCentre || !ball || piece.rail.s[0].cliffCurve
            || piece.rail.s[1].cliffCurve)
          return false;
        piece.rail.byCentre = true;
        --attempt;
        return true;
      };
      try
      {
        stations.march(seeds[k]);
      }
      catch (const Failure&)
      {
        if (!stations.stalled() || !mayCentre(int(k)))
          throw;
        centreChain(int(k));
        return;
      }
      if (piece.rail.law == Law::Chordal)
        roundOfChord(piece.rail, stations);
      if (request.meet > 0)
        meetsAt(piece.rail, stations, request.meet);
      // A carry that stopped short ends where the ball last went on with the rail.
      if (const auto reached = stations.reached(0))
        piece.rail.from = std::max(piece.rail.from, *reached);
      if (const auto reached = stations.reached(1))
        piece.rail.to = std::min(piece.rail.to, *reached);
      // An exact sheet needs no fit: its contacts are its own.
      // A conic's corner patch takes its weights from fitted sheets: none extruded there (a
      // chamfer's corner takes none).
      const bool corner =
        piece.rail.law != Law::Chamfer
        && (piece.joint[0].meets == Meets::Corner || piece.joint[1].meets == Meets::Corner);
      auto exact = exactSheet(piece.rail, stations);
      if (!exact && !corner)
        exact = extrudedSheet(piece.rail, stations);
      const auto contacts = exact ? exactContacts(piece.rail, stations, *exact) : std::nullopt;
      // Locate support crossings from the ball's contacts before fitting the G2 section.
      // A finite spline support's tangent continuation is only G1: its curvature jumps,
      // so fitting a curvature-dependent section across it before splitting cannot converge.
      // The circular pilot uses the same balls and contact parameters. It is always replaced
      // by the requested G2 section after crossing discovery, and never reaches attachment.
      const bool contactPilot = !contacts && attempt == 0 && piece.rail.law == Law::G2
                                && !piece.rail.s[0].cliffCurve && !piece.rail.s[1].cliffCurve
                                && piece.rail.s[0].trim.IsNull() && piece.rail.s[1].trim.IsNull();
      if (contacts)
      {
        piece.fitted = *contacts;
        piece.sheet  = *exact;
      }
      else
      {
        try
        {
          piece.fitted = fit(piece.rail, stations, contactPilot);
        }
        catch (const Failure& failure)
        {
          if (failure.reason == std::string("cancelled") || !centreMends())
            throw;
          centreChain(int(k));
          return;
        }
        piece.sheet = exact ? *exact : Sheet{piece.fitted.blend};
        if (!exact)
          checkCrossing(piece.rail, stations, piece.fitted);
      }
      const bool onCliff = piece.rail.s[0].cliffCurve || piece.rail.s[1].cliffCurve
                           || !piece.rail.s[0].trim.IsNull() || !piece.rail.s[1].trim.IsNull();
      if (attempt == 1 || onCliff)
      {
        if (!contacts)
        {
          Dips dips{piece.rail, piece.fitted};
          if (onCliff || capFoldsOf(int(k), dips))
          {
            if (runsBack(dips) && centreMends())
            {
              centreChain(int(k));
              return;
            }
            checkFold(dips);
          }
        }
        break;
      }
      // The first place where a contact changes between inside its face and outside, and
      // the edge of that face it crosses there: leaving it, or entering it from the face it
      // started on.
      const double  u0 = piece.uAt(0), u1 = piece.uAt(1);
      constexpr int samples = 256;
      double        at = 0, overAt = 0;
      int           side     = -1;
      bool          entering = false, notched = false;
      TopoDS_Edge   over;
      // A contact wholly off its face: on both, the blend rests on neither and does not fit;
      // on one, the ball may roll over onto the next face, which the blend does not do yet.
      int off = 0, offSide = -1;
      for (int s2 = 0; s2 < 2; ++s2)
      {
        Inside& onFace = inside(piece.rail.s[s2].face);
        // On the face's far edge all along counts as off it: the ball rests on the face and that
        // edge at once (a side exactly as high as the ball is wide), and the blend uses the face
        // up as where the side is narrower.
        const auto&       far = farOf(piece.rail, s2);
        BRepAdaptor_Curve farCurve;
        if (far)
          farCurve.Initialize(far->edge);
        const auto onFar = [&](double u) {
          if (!far)
            return false;
          const gp_Pnt2d      uv = piece.fitted.contact[s2]->Value(u);
          const gp_Pnt        p  = footAt(piece.rail.s[s2], uv.X(), uv.Y());
          const Extrema_ExtPC to(p, farCurve);
          if (!to.IsDone() || to.NbExt() == 0)
            return false;
          double nearest = HUGE_VAL;
          for (int i = 1; i <= to.NbExt(); ++i)
            nearest = std::min(nearest, to.SquareDistance(i));
          return std::sqrt(nearest) <= piece.rail.resolved(limit::together);
        };
        const auto out = [&](double u) { return onFace.out(piece.fitted.contact[s2]->Value(u)); };
        bool       wholly = true, along = true;
        for (int i = 0; i <= samples && (wholly || along); i += samples / 16)
        {
          wholly = wholly && out(u0 + (u1 - u0) * i / samples);
          along  = along && onFar(u0 + (u1 - u0) * i / samples);
        }
        if (!wholly && !along)
          continue;
        ++off;
        offSide = s2;
        // A flat face all of whose edges are blended: the blends use the cap up.
        if (GeomAdaptor_Surface(piece.rail.s[s2].surface).GetType() != GeomAbs_Plane)
          continue;
        bool all = true;
        for (TopExp_Explorer e(piece.rail.s[s2].face, TopAbs_EDGE); e.More() && all; e.Next())
          all = std::any_of(pieces.begin(), pieces.end(), [&](const Piece& other) {
            return other.rail.edge.IsSame(e.Current());
          });
        if (all)
          fail(op,
               no_solution,
               "the blend consumes the cap: the face beside its edge vanishes",
               {piece.rail.edge},
               FailureKind::CapConsumed);
      }
      if (off == 2)
        fail(op,
             no_solution,
             "the blend does not fit: it is wider than both faces beside its edge",
             {piece.rail.edge},
             FailureKind::Overflow);
      if (off == 1)
      {
        // A face narrower than the blend: the ball rests on its far edge (a cliff).
        // Between stops at both ends, or round a closed edge (a band narrower than the blend).
        const bool stopped =
          piece.joint[0].meets == Meets::Stop && piece.joint[1].meets == Meets::Stop;
        // Not on an edge blended too: it is not there for the ball to rest on, nor to trim by.
        const auto blended = [&](const TopoDS_Edge& edge) {
          return std::any_of(pieces.begin(), pieces.end(), [&](const Piece& other) {
            return other.rail.edge.IsSame(edge);
          });
        };
        const auto& far = farOf(piece.rail, offSide);
        if (far && !blended(far->edge) && (stopped || piece.rail.closed)
            && cliffOf(piece.rail, offSide, *far))
        {
          --attempt;
          continue;
        }
        // Else the face beyond the far edge trims the blend (a notch), between stops or round a
        // closed edge.
        if (far && !blended(far->edge) && (stopped || piece.rail.closed))
        {
          TopoDS_Face beyond;
          for (const auto& face : facesOfEdge.FindFromKey(far->edge))
            if (!face.IsSame(piece.rail.s[offSide].face))
              beyond = TopoDS::Face(face);
          if (!beyond.IsNull() && !beyond.IsSame(piece.rail.s[1 - offSide].face))
          {
            piece.rail.s[offSide].trim       = beyond;
            piece.rail.s[offSide].trimEdge   = far->edge;
            piece.rail.s[offSide].cliffWhole = far->whole;
            --attempt;
            continue;
          }
        }
        // A part of a split, wholly past the place its contact leaves the face: that is seen
        // along the whole edge, which came first.
        if (piece.rail.parameterSupports)
          off = 0;
        else if ((piece.joint[0].meets == Meets::Stop) != (piece.joint[1].meets == Meets::Stop)
                 && piece.joint[piece.joint[0].meets == Meets::Stop ? 1 : 0].meets == Meets::Joint)
        {
          ending.push_back(TopoDS::Edge(piece.rail.edge));
          again = true;
          return;
        }
        else
          fail(op,
               unsupported_geometry,
               "the blend's contact runs off the face beside the edge, onto faces v2 does not roll "
               "it over yet",
               {piece.rail.edge});
      }
      for (int s2 = 0; s2 < 2 && side < 0; ++s2)
      {
        Inside&    onFace = inside(piece.rail.s[s2].face);
        const auto out = [&](double u) { return onFace.out(piece.fitted.contact[s2]->Value(u)); };
        bool       startsOut = out(u0);
        double     previous  = u0;
        bool       trimmed   = false;
        for (int i = 1; i <= samples && side < 0; ++i)
        {
          trimmed        = false;
          const double u = u0 + (u1 - u0) * i / samples;
          // ON is a tolerance band, not a side of the boundary. Retain
          // the last definite side so a delayed exit cannot discard
          // the geometric crossing while traversing that band.
          const auto state = onFace.classify(piece.fitted.contact[s2]->Value(u));
          if (state == TopAbs_ON)
            continue;
          if ((state == TopAbs_OUT) == startsOut)
          {
            previous = u;
            continue;
          }
          // Keep the scan's forward lookahead, with its lower bound
          // still on the definite side before the boundary band.
          const double past = std::min(u1, u + (u1 - u0) / samples);
          for (TopExp_Explorer e(piece.rail.s[s2].face, TopAbs_EDGE); e.More(); e.Next())
          {
            // A seam has the face on both sides: crossing it, the contact stays on it.
            if (e.Current().IsSame(piece.rail.edge)
                || BRep_Tool::IsClosed(TopoDS::Edge(e.Current()), piece.rail.s[s2].face))
              continue;
            // Over another blended edge the contact runs past the corner or mitre where
            // the other blend's contact trims it: not a face it rolls on.
            if (std::any_of(pieces.begin(), pieces.end(), [&](const Piece& other) {
                  return other.rail.edge.IsSame(e.Current());
                }))
            {
              if (crossing(piece.fitted.contact[s2],
                           previous,
                           past,
                           previous,
                           TopoDS::Edge(e.Current()),
                           piece.rail.s[s2].face))
              {
                trimmed = true;
                break;
              }
              continue;
            }
            const auto found = crossing(piece.fitted.contact[s2],
                                        previous,
                                        past,
                                        previous,
                                        TopoDS::Edge(e.Current()),
                                        piece.rail.s[s2].face);
            if (!found)
              continue;
            // At an end of the piece: the place a split already made it end.
            const auto placeAt = [&](double v) {
              const gp_Pnt2d uv = piece.fitted.contact[s2]->Value(v);
              return supportPoint(piece.rail.s[s2].surface, uv.X(), uv.Y());
            };
            if (placeAt(found->first).Distance(placeAt(u0)) <= piece.rail.resolved(limit::together)
                || placeAt(found->first).Distance(placeAt(u1))
                     <= piece.rail.resolved(limit::together))
            {
              trimmed = true;
              break;
            }
            // Over an edge where the faces are not tangent the contact meets what ends the
            // blend: near an end of the piece, that end's stopping face trims it there.
            bool smooth = false;
            for (const auto& next : facesOfEdge.FindFromKey(e.Current()))
            {
              if (next.IsSame(piece.rail.s[s2].face))
                continue;
              double     f, l;
              const auto pc =
                BRep_Tool::CurveOnSurface(TopoDS::Edge(e.Current()), TopoDS::Face(next), f, l);
              if (pc.IsNull())
                continue;
              Side beyond          = piece.rail.s[s2];
              beyond.face          = TopoDS::Face(next);
              beyond.surface       = BRep_Tool::Surface(beyond.face);
              beyond.outward       = next.Orientation() == TopAbs_REVERSED ? -1 : 1;
              const gp_Pnt2d here  = pc->Value(found->second),
                             there = piece.fitted.contact[s2]->Value(found->first);
              Jet a, b;
              smooth = jetAt(beyond, here.X(), here.Y(), a)
                       && jetAt(piece.rail.s[s2], there.X(), there.Y(), b)
                       && SectionPolicy::tangentFor(a.n.Angle(b.n), piece.rail.scale);
            }
            if (!smooth)
            {
              // Into a hole of its face, a feature that can be lifted off and put back
              // (liftFeatures): the blend runs into its walls.
              for (TopoDS_Iterator w(TopoDS::Face(piece.rail.s[s2].face.Oriented(TopAbs_FORWARD)));
                   w.More();
                   w.Next())
              {
                bool holds = false;
                for (TopoDS_Iterator x(w.Value()); x.More() && !holds; x.Next())
                  holds = x.Value().IsSame(e.Current());
                if (!holds
                    || !holeOf(w.Value(),
                               TopoDS::Face(piece.rail.s[s2].face.Oriented(TopAbs_FORWARD))))
                  continue;
                std::vector<TopoDS_Shape> walls;
                for (const auto& next : facesOfEdge.FindFromKey(e.Current()))
                  if (!next.IsSame(piece.rail.s[s2].face))
                    walls.push_back(next);
                // Only a feature that is put back as it is: a hole (its rim convex)
                // beside a convex blend, a boss (its foot concave) beside a concave one.
                if (walls.size() != 1
                    || convexity(op,
                                 TopoDS::Edge(e.Current()),
                                 piece.rail.s[s2].face,
                                 TopoDS::Face(walls.front()),
                                 piece.rail.scale)
                         != (piece.rail.sigma < 0 ? -1 : 1))
                  continue;
                const auto feature = featuresOf(body, walls);
                if (std::none_of(feature.begin(), feature.end(), [&](const TopoDS_Shape& f) {
                      return std::any_of(pieces.begin(), pieces.end(), [&](const Piece& p) {
                        return f.IsSame(p.rail.s[0].face) || f.IsSame(p.rail.s[1].face);
                      });
                    }))
                  throw Collision(op,
                                  "the blend's contact runs into a hole or a boss in its face",
                                  piece.rail.edge,
                                  walls);
              }
              // An edge of the face at an end vertex of the blended edge is what stops the
              // blend there, however far from the end the contact crosses it (at a wedge of
              // angle w, its offset over tan w); else near an end, within three radii.
              const double reach = 3 * piece.rail.scale;
              gp_Pnt       p0, p1, pc;
              piece.rail.curve.D0(piece.rail.first, p0);
              piece.rail.curve.D0(piece.rail.last, p1);
              piece.rail.curve.D0(piece.rail.t(found->first), pc);
              // Or over the edge onto a face that stops the blend at an end of its edge
              // (a top narrowing to an oblique end): that stop trims it there,
              // wherever along the edge the contact reaches it.
              bool          ontoStop = false;
              TopoDS_Vertex ends[2];
              TopExp::Vertices(TopoDS::Edge(piece.rail.edge.Oriented(TopAbs_FORWARD)),
                               ends[0],
                               ends[1]);
              const bool atEnd =
                std::any_of(std::begin(ends), std::end(ends), [&](const TopoDS_Vertex& v) {
                  const auto& at = edgesOfVertex.FindFromKey(v);
                  return std::any_of(at.begin(), at.end(), [&](const TopoDS_Shape& edge) {
                    return edge.IsSame(e.Current());
                  });
                });
              for (const auto& next : facesOfEdge.FindFromKey(e.Current()))
                for (int end = 0; end < 2 && !next.IsSame(piece.rail.s[0].face)
                                  && !next.IsSame(piece.rail.s[1].face);
                     ++end)
                {
                  const TopoDS_Vertex& v = ends[end];
                  if (piece.joint[end].meets != Meets::Stop || v.IsNull()
                      || !edgesOfVertex.Contains(v))
                    continue;
                  for (const auto& at : edgesOfVertex.FindFromKey(v))
                    for (const auto& face : facesOfEdge.FindFromKey(at))
                      ontoStop = ontoStop || face.IsSame(next);
                }
              if (atEnd || ontoStop || pc.Distance(p0) <= reach || pc.Distance(p1) <= reach)
              {
                trimmed = true;
                break;
              }
              // Partway over a convex edge the blend goes on as it is on the face carried

              if (!convexOver(piece.rail,
                              s2,
                              TopoDS::Edge(e.Current()),
                              *found,
                              piece.fitted.contact[s2],
                              previous))
                fail(op,
                     unsupported_geometry,
                     "a blend's contact runs off its face partway along the edge, which v2 does "
                     "not do yet",
                     {piece.rail.edge});
              notched = true;
            }
            side     = s2;
            at       = found->first;
            overAt   = found->second;
            over     = TopoDS::Edge(e.Current());
            entering = startsOut;
            break;
          }
          if (side < 0 && !trimmed)
            fail(op,
                 unsupported_geometry,
                 "a blend's contact leaves its face where v2 does not find the edge it crosses",
                 {piece.rail.edge});
          if (trimmed)
          {
            // A stop at one end may be followed by a support crossing at the other.
            // Continue from this side of the ignored transition, rather than dropping
            // the rest of the contact after its first entry into the face.
            startsOut = out(u);
            previous  = u;
          }
        }
      }
      // Over a seam of its face the contact stays on the face, a period on, which the test
      // above does not see. Where it crosses one between the piece's ends (a rim whose
      // dihedral turns along it, so that its sections lean off the seam), the face's new
      // boundary has its vertex there on the seam: the piece is split there as where it
      // crosses onto the next face, the part after rolling on the same face a period over.
      bool     overSeam = false;
      gp_Vec2d seamShift; // the period the part off the face's own parameters is carried by
      if (side < 0 && !contacts)
        for (int s2 = 0; s2 < 2 && side < 0; ++s2)
        {
          const TopoDS_Face&          face    = piece.rail.s[s2].face;
          const Handle(Geom_Surface)& surface = piece.rail.s[s2].surface;
          // Only a face on a closed surface has a seam.
          if (!surface->IsUClosed() && !surface->IsVClosed())
            continue;
          const auto pointAt = [&](double u) {
            const gp_Pnt2d uv = piece.fitted.contact[s2]->Value(u);
            return surface->Value(uv.X(), uv.Y());
          };
          for (TopExp_Explorer e(face, TopAbs_EDGE); e.More() && side < 0; e.Next())
          {
            const TopoDS_Edge seam = TopoDS::Edge(e.Current());
            if (!BRep_Tool::IsClosed(seam, face)
                || !acrossSeam(piece.fitted.contact[s2], seam, face, surface)
                || poleInStrip(piece, s2, seam))
              continue;
            const auto found =
              crossing(piece.fitted.contact[s2], u0, u1, (u0 + u1) / 2, seam, face);
            if (!found)
              continue;
            const gp_Pnt there = pointAt(found->first);
            const double near  = piece.rail.resolved(limit::together);
            if (there.Distance(pointAt(u0)) <= near || there.Distance(pointAt(u1)) <= near)
              continue;
            // The part that lies off the face's own parameters, before the crossing or
            // after it, is the one carried a period over onto them.
            double ua, ub, va, vb;
            BRepTools::UVBounds(face, ua, ub, va, vb);
            // Each parameter's slack a share of its own range: an angle's and a length's
            // ranges differ by the face's size.
            const double su = 1e-9 * (ub - ua), sv = 1e-9 * (vb - va);
            const auto   within = [&](const gp_Pnt2d& p) {
              return p.X() >= ua - su && p.X() <= ub + su && p.Y() >= va - sv && p.Y() <= vb + sv;
            };
            const gp_Pnt2d ahead =
              piece.fitted.contact[s2]->Value(found->first + (u1 - found->first) * 1e-3);
            const gp_Pnt2d behind =
              piece.fitted.contact[s2]->Value(found->first - (found->first - u0) * 1e-3);
            if (within(ahead) == within(behind))
              continue;
            const gp_Pnt2d off      = within(ahead) ? behind : ahead;
            const double   tu       = surface->IsUPeriodic() ? surface->UPeriod() : 0,
                         tv         = surface->IsVPeriodic() ? surface->VPeriod() : 0;
            const gp_Vec2d shifts[] = {{tu, 0}, {-tu, 0}, {0, tv}, {0, -tv}};
            for (const auto& shift : shifts)
              if (shift.SquareMagnitude() > 0 && within(off.Translated(shift)))
              {
                seamShift = shift;
                side      = s2;
                at        = found->first;
                overAt    = found->second;
                over      = seam;
                overSeam  = true;
                entering  = within(ahead);
                break;
              }
            if (side < 0)
              fail(op,
                   unsupported_geometry,
                   "a blend's contact crosses its face's seam where v2 does not find the face's "
                   "period beyond",
                   {piece.rail.edge});
          }
        }
      if (side < 0)
      {
        if (contactPilot)
        {
          piece.fitted = fit(piece.rail, stations);
          piece.sheet  = exact ? *exact : Sheet{piece.fitted.blend};
        }
        // An exact sheet's contacts are lines or circles at constant speed: they cannot fold.
        if (!contacts)
        {
          Dips dips{piece.rail, piece.fitted};
          if (capFoldsOf(int(k), dips))
          {
            if (runsBack(dips) && centreMends())
            {
              centreChain(int(k));
              return;
            }
            checkFold(dips);
          }
        }
        break;
      }
      // The next face across that edge, smooth with this one there.
      TopoDS_Face next = overSeam ? piece.rail.s[side].face : TopoDS_Face();
      for (const auto& face : facesOfEdge.FindFromKey(over))
        if (!face.IsSame(piece.rail.s[side].face))
          next = TopoDS::Face(face);
      const double t = piece.rail.t(at);
      Station      station;
      if (next.IsNull() || !stations.at(t, station))
        fail(op, unsupported_geometry, "a blend's contact runs off its face", {piece.rail.edge});
      station.t = t;
      // The edge's part before the crossing and its part after; the one off the face the
      // edge lies on rolls on the next face, from the station there moved onto it.
      Piece first = piece, second = piece;
      if (!first.rail.parameterSupports)
        first.rail.parameterSupports = std::array<Side, 2>{piece.rail.s[0], piece.rail.s[1]};
      second.rail.parameterSupports = first.rail.parameterSupports;
      first.rail.last = first.rail.to = t;
      second.rail.first = second.rail.from = t;
      Piece&  moved                        = entering ? first : second;
      Station onNext                       = station;
      if (notched)
      {
        moved.rail.s[side].trim       = next;
        moved.rail.s[side].trimEdge   = over;
        moved.rail.s[side].cliffWhole = false;
      }
      else
      {
        moved.rail.s[side].face    = next;
        moved.rail.s[side].surface = BRep_Tool::Surface(next);
        moved.rail.s[side].outward = next.Orientation() == TopAbs_REVERSED ? -1 : 1;
        GeomAPI_ProjectPointOnSurf onto(station.foot[side], moved.rail.s[side].surface);
        if (onto.NbPoints() == 0 || onto.LowerDistance() > piece.rail.resolved(limit::together))
          fail(op,
               unsupported_geometry,
               "a blend's contact runs off its face onto one it does not meet",
               {piece.rail.edge});
        onto.LowerDistanceParameters(onNext.uv[side][0], onNext.uv[side][1]);
        // Across its own seam, a period over (seamShift).
        if (overSeam)
        {
          onNext.uv[side][0] = station.uv[side][0] + seamShift.X();
          onNext.uv[side][1] = station.uv[side][1] + seamShift.Y();
        }
        else
        {
          // In the period of the next face's own parameters where the contact crosses onto
          // it: two faces on one periodic surface meet at the seam's place a period apart.
          double     f, l;
          const auto onOver = BRep_Tool::CurveOnSurface(over, next, f, l);
          if (!onOver.IsNull())
          {
            const gp_Pnt2d near = nearPeriod(gp_Pnt2d(onNext.uv[side][0], onNext.uv[side][1]),
                                             onOver->Value(overAt),
                                             moved.rail.s[side].surface);
            onNext.uv[side][0]  = near.X();
            onNext.uv[side][1]  = near.Y();
          }
        }
      }
      Jet before, after;
      if (!jetAt(piece.rail.s[side], station.uv[side][0], station.uv[side][1], before)
          || !jetAt(moved.rail.s[side], onNext.uv[side][0], onNext.uv[side][1], after)
          || !SectionPolicy::tangentFor(before.n.Angle(after.n), piece.rail.scale))
        fail(op,
             unsupported_geometry,
             "a blend's contact runs over an edge where its faces are not tangent, which v2 does "
             "not do yet",
             {piece.rail.edge});
      // `piece` refers into `pieces`, which must not grow past what it holds.
      if (pieces.size() >= pieces.capacity())
        fail(op,
             unsupported_geometry,
             "the blend crosses too many faces along one edge",
             {piece.rail.edge});
      const int  index = int(pieces.size());
      const auto outer = piece.joint[1];
      first.joint[1]   = {Meets::Split, index, 0, false, side, over};
      second.joint[0]  = {Meets::Split, int(k), 1, false, side, over};
      seeds[k]         = entering ? onNext : station;
      seeds.push_back(entering ? station : onNext);
      piece = first;
      // What met the edge's far end meets the second part's: on a closed edge, the first
      // part's own start.
      if (outer.meets == Meets::Joint || outer.meets == Meets::Miter || outer.meets == Meets::Split)
      {
        pieces[outer.piece].joint[outer.end].piece = index;
        pieces[outer.piece].joint[outer.end].end   = 1;
      }
      if (outer.meets == Meets::Corner)
        for (auto& [p2, e2] : corners[outer.piece].ends)
          if (p2 == int(k) && e2 == 1)
            p2 = index;
      pieces.push_back(second);
      // The first part is looked at again, as the second will be: its other contact may cross
      // an edge too (the other face's seam just before a joint at a cone's, an edge of the
      // other support just after a joint of crossed extrusions).
      attempt = -1;
    }
  }
  crossFolds(folds);
  // Where a contact on a flat cap folds over itself, the edge's blend is cut at the two ends of
  // each loop, where the balls are one, and the loop is patched (Attach::fold): a piece on the
  // one sheet between each loop and the next.
  for (const auto& [edge, found] : folds)
  {
    int        k = edge;
    const bool closed =
      pieces[edge].joint[1].meets == Meets::Joint && pieces[edge].joint[1].piece == edge;
    for (const Fold& fold : found)
    {
      if (pieces.size() >= pieces.capacity())
        fail(op,
             unsupported_geometry,
             "the blend folds too often along one edge",
             {pieces[k].rail.edge});
      Piece first = pieces[k], second = pieces[k];
      first.rail.last   = first.rail.t(fold.from);
      second.rail.first = second.rail.t(fold.to);
      // Its own copy of the sheet, so that an edge both parts meet on (a closed edge's seam) has
      // a curve on each face's surface, not two on one.
      second.fitted.blend  = Handle(Geom_BSplineSurface)::DownCast(first.fitted.blend->Copy());
      second.sheet.surface = second.fitted.blend;
      const int  index     = int(pieces.size());
      const auto outer     = pieces[k].joint[1];
      first.joint[1]       = {Meets::Fold, index, 0, false, fold.side, {}};
      second.joint[0]      = {Meets::Fold, k, 1, false, fold.side, {}};
      pieces[k]            = first;
      pieces.push_back(second);
      seeds.push_back(std::nullopt);
      // What met the edge's far end meets the second part's (the first's start, on a closed edge).
      if (outer.meets == Meets::Joint || outer.meets == Meets::Miter || outer.meets == Meets::Split
          || outer.meets == Meets::Fold)
      {
        pieces[outer.piece].joint[outer.end].piece = index;
        pieces[outer.piece].joint[outer.end].end   = 1;
      }
      if (outer.meets == Meets::Corner)
        for (auto& [p2, e2] : corners[outer.piece].ends)
          if (p2 == k && e2 == 1)
            p2 = index;
      k = index;
    }
    // On a closed edge the last part runs on into the first across the edge's vertex, where the
    // contact on the face with a seam there crosses it: one blend face (joinSplits), as the
    // ball rolls on through the vertex.
    if (!closed)
      continue;
    TopoDS_Vertex vertex =
      TopExp::FirstVertex(TopoDS::Edge(pieces[edge].rail.edge.Oriented(TopAbs_FORWARD)));
    for (int side = 0; side < 2; ++side)
    {
      const TopoDS_Face& face = pieces[edge].rail.s[side].face;
      for (const auto& other : edgesOfVertex.FindFromKey(vertex))
      {
        if (other.IsSame(pieces[edge].rail.edge) || !BRep_Tool::IsClosed(TopoDS::Edge(other), face))
          continue;
        pieces[k].joint[1]    = {Meets::Split, edge, 0, false, side, TopoDS::Edge(other)};
        pieces[edge].joint[0] = {Meets::Split, k, 1, false, side, TopoDS::Edge(other)};
      }
    }
  }
}

void Blend::crossFolds(const std::vector<std::pair<int, std::vector<Fold>>>& folds)
{
  // A loop that runs on across a joint of a chain (a closed outline split into edges, its
  // tightest bend at their vertex): the contact of the piece where it runs back crosses the next
  // piece's contact on the same cap. Both pieces are cut where they cross, the narrowest such
  // crossing round each place, and the loop between is patched as one within a piece is.
  struct Cross
  {
    int    a, ea, b, eb, side, sideB;
    double ua, ub;
  };

  std::vector<Cross> crosses;
  const auto         length = [](const Piece& p, double u0, double u1) {
    const double t0 = p.rail.t(u0), t1 = p.rail.t(u1);
    return GCPnts_AbscissaPoint::Length(p.rail.curve, std::min(t0, t1), std::max(t0, t1));
  };
  // Whether a place of piece k lies in a crossing's loop.
  const auto inLoop = [](const Cross& c, int k, double at) {
    return (c.a == k && (c.ea == 1 ? at > c.ua : at < c.ua))
           || (c.b == k && (c.eb == 1 ? at > c.ub : at < c.ub));
  };
  for (const auto& [k, places] : open)
    for (const Fold& place : places)
    {
      if (std::any_of(crosses.begin(), crosses.end(), [&](const Cross& c) {
            return inLoop(c, k, place.from);
          }))
        continue;
      const Piece&         one = pieces[k];
      std::optional<Cross> best;
      double               shortest = HUGE_VAL;
      for (int e = 0; e < 2; ++e)
      {
        const auto& joint = one.joint[e];
        if (joint.meets != Meets::Joint || joint.piece == k)
          continue;
        const Piece& two = pieces[joint.piece];
        // The other piece's contact on the same cap: the one that ends where this one does.
        const gp_Pnt2d here  = one.fitted.contact[place.side]->Value(one.uAt(e));
        int            sideB = -1;
        for (int s2 = 0; s2 < 2; ++s2)
          if (!two.fitted.contact[s2].IsNull()
              && two.rail.s[s2].face.IsSame(one.rail.s[place.side].face)
              && (sideB < 0
                  || two.fitted.contact[s2]->Value(two.uAt(joint.end)).Distance(here)
                       < two.fitted.contact[sideB]->Value(two.uAt(joint.end)).Distance(here)))
            sideB = s2;
        if (sideB < 0)
          continue;
        const Geom2dAPI_InterCurveCurve crossing(
          new Geom2d_TrimmedCurve(one.fitted.contact[place.side], one.uAt(0), one.uAt(1)),
          new Geom2d_TrimmedCurve(two.fitted.contact[sideB], two.uAt(0), two.uAt(1)),
          tol::parameter);
        for (int i = 1; i <= crossing.NbPoints(); ++i)
        {
          const auto&  point = crossing.Intersector().Point(i);
          const double ua = point.ParamOnFirst(), ub = point.ParamOnSecond();
          // Not where the two contacts meet at the joint itself.
          const auto& cap = one.rail.s[place.side].surface;
          if (cap->Value(point.Value().X(), point.Value().Y())
                  .Distance(cap->Value(here.X(), here.Y()))
                <= one.rail.resolved(limit::together)
              || !(e == 1 ? ua < place.from : ua > place.from))
            continue;
          const double loop = length(one, ua, one.uAt(e)) + length(two, two.uAt(joint.end), ub);
          if (loop < shortest)
          {
            shortest = loop;
            best     = Cross{int(k), e, joint.piece, joint.end, place.side, sideB, ua, ub};
          }
        }
      }
      if (!best)
        refuseFold(k, "a contact on the cap folds where v2 finds no loop");
      crosses.push_back(*best);
    }
  // Each place where a contact runs back lies in one loop, and no loop overlaps one within a piece.
  for (const auto& [k, places] : open)
    for (const Fold& place : places)
      if (std::none_of(crosses.begin(), crosses.end(), [&](const Cross& c) {
            return inLoop(c, k, place.from);
          }))
        refuseFold(k, "a contact on the cap folds where v2 finds no loop");
  for (const auto& c : crosses)
    for (const auto& [k, found] : folds)
      for (const Fold& fold : found)
        if (inLoop(c, k, fold.from) || inLoop(c, k, fold.to))
          refuseFold(k, "loops where a contact folds overlap");
  for (const auto& c : crosses)
  {
    if (pieces[c.a].joint[c.ea].meets != Meets::Joint
        || pieces[c.b].joint[c.eb].meets != Meets::Joint)
      refuseFold(c.a, "loops where a contact folds overlap");
    for (const auto& [k, e, at] : {std::tuple{c.a, c.ea, c.ua}, std::tuple{c.b, c.eb, c.ub}})
    {
      Rail& rail                        = pieces[k].rail;
      (e == 1 ? rail.last : rail.first) = rail.t(at);
      if (!(rail.last > rail.first))
        refuseFold(k, "the loop where a contact folds takes a whole edge");
    }
    pieces[c.a].joint[c.ea] = {Meets::Fold, c.b, c.eb, c.side != c.sideB, c.side, {}};
    pieces[c.b].joint[c.eb] = {Meets::Fold, c.a, c.ea, c.side != c.sideB, c.sideB, {}};
  }
}

void Blend::chains()
{
  std::vector<int> chain(pieces.size(), -1);
  for (int k = 0; k < int(pieces.size()); ++k)
  {
    if (chain[k] >= 0)
      continue;
    std::vector<int> members, todo = {k};
    chain[k] = k;
    while (!todo.empty())
    {
      const int p = todo.back();
      todo.pop_back();
      members.push_back(p);
      for (const auto& joint : pieces[p].joint)
        if ((joint.meets == Meets::Joint || joint.meets == Meets::Split) && chain[joint.piece] < 0)
        {
          chain[joint.piece] = k;
          todo.push_back(joint.piece);
        }
    }
    const bool fitted =
      std::any_of(members.begin(), members.end(), [&](int p) { return !pieces[p].exact(); });
    if (fitted)
      for (int p : members)
      {
        // A piece whose sheet was exact is fitted now, to share its chain's parameter.
        if (pieces[p].fitted.blend.IsNull())
        {
          Stations stations(pieces[p].rail);
          stations.march(seeds[p]);
          pieces[p].fitted = fit(pieces[p].rail, stations);
        }
        pieces[p].sheet = Sheet{pieces[p].fitted.blend};
      }
  }
}

void Blend::joinSplits()
{
  joinSplitPieces(op, pieces, seeds);
}

void Blend::pairJoints()
{
  // Which sides pair at each joint: the contacts that end at the same place.
  for (auto& piece : pieces)
    for (int end = 0; end < 2; ++end)
    {
      auto& joint = piece.joint[end];
      if (joint.meets != Meets::Joint)
        continue;
      const Piece& other = pieces[joint.piece];
      const gp_Pnt here  = piece.sheet.contact(0)->Value(piece.sheet.rail(piece.uAt(end)));
      const auto   there = [&](int side) {
        return other.sheet.contact(side)->Value(other.sheet.rail(other.uAt(joint.end)));
      };
      joint.swapped = here.Distance(there(1)) < here.Distance(there(0));
    }
}

std::vector<TopoDS_Shape> Blend::foldedJoints() const
{
  std::vector<int> bad;
  for (int k = 0; k < int(pieces.size()); ++k)
    for (int end = 0; end < 2; ++end)
    {
      const auto& joint = pieces[k].joint[end];
      if (joint.meets != Meets::Joint)
        continue;
      const Piece& piece = pieces[k];
      const Piece& other = pieces[joint.piece];
      const auto   at    = [](const Piece& p, int side, int e) {
        return p.sheet.contact(side)->Value(p.sheet.rail(p.uAt(e)));
      };
      for (int side = 0; side < 2; ++side)
        if (at(piece, side, end).Distance(at(other, side ^ int(joint.swapped), joint.end))
            > piece.rail.resolved(limit::together))
          bad.push_back(k);
    }
  // Each such chain whole, through its joints and splits; only balls of a law the centres give.
  std::vector<TopoDS_Shape> edges;
  std::vector<bool>         seen(pieces.size(), false);
  while (!bad.empty())
  {
    const int k = bad.back();
    bad.pop_back();
    if (seen[k])
      continue;
    seen[k]          = true;
    const Rail& rail = pieces[k].rail;
    const bool  ball = rail.law == Law::Conic || rail.law == Law::G2;
    if (rail.byCentre || !ball || rail.s[0].cliffCurve || rail.s[1].cliffCurve)
      return {};
    if (std::none_of(edges.begin(), edges.end(), [&](const TopoDS_Shape& e) {
          return e.IsSame(rail.edge);
        }))
      edges.push_back(rail.edge);
    for (const auto& joint : pieces[k].joint)
      if (joint.meets == Meets::Joint || joint.meets == Meets::Split)
        bad.push_back(joint.piece);
  }
  return edges;
}

//! The full round of meeting rounds (Request::meet) split where it touches the face between: the
//! line of its sheet through the ball's foot there, between the face's boundary.
void Blend::splitMeeting(Blended& out)
{
  const Piece& piece = pieces.front();
  Stations     stations(piece.rail);
  stations.march();
  Station s;
  if (!stations.at((piece.rail.first + piece.rail.last) / 2, s))
    fail(op, no_solution, "the full round has no ball in its middle", {piece.rail.edge});
  const Sheet&               sheet = piece.sheet;
  GeomAPI_ProjectPointOnSurf onSheet(s.middle, sheet.surface);
  if (onSheet.NbPoints() == 0)
    fail(op, no_solution, "the full round does not touch the face between", {piece.rail.edge});
  double u, v;
  onSheet.LowerDistanceParameters(u, v);
  const double across = sheet.swapped ? u : v;
  TopoDS_Face  face;
  for (TopExp_Explorer f(out.shape, TopAbs_FACE); f.More() && face.IsNull(); f.Next())
    if (BRep_Tool::Surface(TopoDS::Face(f.Current())) == sheet.surface)
      face = TopoDS::Face(f.Current());
  if (face.IsNull())
    fail(op, no_solution, "the full round's face is missing", {piece.rail.edge});
  // Where the line meets the face's boundary, along the rail.
  const Handle(Geom2d_Line) line =
    new Geom2d_Line(sheet.uv(0, across), sheet.swapped ? gp_Dir2d(0, 1) : gp_Dir2d(1, 0));
  std::vector<double> ends;
  for (TopExp_Explorer e(face, TopAbs_EDGE); e.More(); e.Next())
  {
    double     f, l;
    const auto pcurve = BRep_Tool::CurveOnSurface(TopoDS::Edge(e.Current()), face, f, l);
    if (pcurve.IsNull())
      continue;
    Geom2dAPI_InterCurveCurve meet(new Geom2d_TrimmedCurve(line, -1e6, 1e6),
                                   new Geom2d_TrimmedCurve(pcurve, f, l),
                                   tol::parameter);
    for (int i = 1; i <= meet.NbPoints(); ++i)
      ends.push_back(ElCLib::Parameter(line->Lin2d(), meet.Point(i)));
  }
  if (ends.size() < 2)
    fail(op,
         no_solution,
         "the line where the full round touches the face between does not cross it",
         {piece.rail.edge});
  const auto [lo, hi] = std::minmax_element(ends.begin(), ends.end());
  const auto  curve   = sheet.swapped ? sheet.surface->UIso(across) : sheet.surface->VIso(across);
  TopoDS_Edge edge    = BRepBuilderAPI_MakeEdge(curve, *lo, *hi);
  BRep_Builder().UpdateEdge(edge, line, face, tol::length);
  BRepAlgoAPI_Splitter           split;
  NCollection_List<TopoDS_Shape> arguments, tools;
  arguments.Append(out.shape);
  tools.Append(edge);
  split.SetArguments(arguments);
  split.SetTools(tools);
  split.Build();
  if (!split.IsDone())
    fail(op,
         no_solution,
         "the full round cannot be split where it touches the face between",
         {piece.rail.edge});
  TopoDS_Shape solid;
  for (TopExp_Explorer x(split.Shape(), TopAbs_SOLID); x.More(); x.Next())
    solid = x.Current();
  if (solid.IsNull())
    fail(op, no_solution, "splitting the full round gave no solid", {piece.rail.edge});
  out.history->Merge(split.History());
  out.shape = solid;
}

Blended Blend::result()
{
  Blended                         out;
  const Handle(BRepTools_History) made = new BRepTools_History;
  out.shape = pieces.empty() ? body : attach(body, pieces, corners, made);
  // The input's history: into the copy, then through the blends.
  NCollection_List<TopoDS_Shape> inputs;
  inputs.Append(request.body);
  out.history = new BRepTools_History(inputs, copy);
  out.history->Merge(made);
  if (request.meet > 0)
    splitMeeting(out);
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> result;
  TopExp::MapShapes(out.shape, TopAbs_FACE, result);
  TopExp::MapShapes(out.shape, TopAbs_EDGE, result);
  for (const auto type : {TopAbs_FACE, TopAbs_EDGE})
  {
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> parts;
    TopExp::MapShapes(request.body, type, parts);
    for (int i = 1; i <= parts.Extent(); ++i)
    {
      const TopoDS_Shape inCopy = copy.ModifiedShape(parts(i));
      if (!result.Contains(inCopy))
        continue;
      out.unchanged.emplace_back(parts(i), inCopy);
      (type == TopAbs_FACE ? out.keptFaces : out.keptEdges)++;
    }
  }
  out.rails      = int(pieces.size());
  out.propagated = int(edges.size() - request.edges.size());
  for (const auto& piece : pieces)
    out.analytic += piece.exact();
  return out;
}

static Blended blendCopy(const Request&                   request,
                         BRepBuilderAPI_Copy&             copy,
                         const std::vector<TopoDS_Shape>& byCentre)
{
  Blend b(request, copy, byCentre);
  // Each stage a phase of the kernel's profile, where one is asked for (diagnostics).
  const auto stage = [&](const char* name, auto&& run) {
    checkCancelled(request.operation);
    BooleanPhase phase(name);
    run();
  };
  stage("fillet2.rails", [&] { b.rails(); });
  stage("fillet2.joints", [&] { b.joints(); });
  stage("fillet2.carries", [&] { b.carries(); });
  stage("fillet2.sheets", [&] {
    // Made again from the pieces as they were, with each chain centreChain marked, until
    // none is marked anew.
    const std::vector<Piece>  pieces  = b.pieces;
    const std::vector<Corner> corners = b.corners;
    // A fold set aside for a loop across a joint came first along the chain: where a later part
    // fails, the fold is what is refused.
    try
    {
      for (b.sheets(); b.again; b.sheets())
      {
        b.again = false;
        b.open.clear();
        b.pieces.assign(pieces.begin(), pieces.end());
        b.corners = corners;
        b.absorbEnds();
      }
    }
    catch (const Failure& failure)
    {
      if (failure.reason == std::string("cancelled") || b.open.empty())
        throw;
      b.refuseFold(b.open.front().first, failure.message.c_str());
    }
  });
  stage("fillet2.chains", [&] { b.chains(); });
  stage("fillet2.joinSplits", [&] { b.joinSplits(); });
  stage("fillet2.pairJoints", [&] { b.pairJoints(); });
  // Two blends of a chain that end in different balls at their joint: the balls whose planes
  // pass through the edge's points fold back along it (an acute groove on a tight bend), and the
  // chain's balls are found through their centres instead, the same balls (Rail::byCentre).
  if (byCentre.empty())
    if (auto folded = b.foldedJoints(); !folded.empty())
    {
      for (const auto& edge : b.byCentre)
        if (std::none_of(folded.begin(), folded.end(), [&](const TopoDS_Shape& e) {
              return e.IsSame(edge);
            }))
          folded.push_back(edge);
      return blendCopy(request, copy, folded);
    }
  BooleanPhase phase("fillet2.attach");
  return b.result();
}

} // namespace BRepExt::blend2
