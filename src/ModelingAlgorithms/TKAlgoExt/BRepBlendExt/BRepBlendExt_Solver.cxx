// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Blend (BRepBlendExt_Internal.hxx): a request's blends, made by the stages on a copy of the body
// (BRepBlendExt_Stages.cxx), and the routes around them: features lifted off and put back,
// overruns replayed, one kind of edge before the other.

#include <GeomAlgoExt_Cancel.hxx>
#include <BRepBlendExt_Stages.hxx>
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

namespace BRepExt::blending
{

using namespace detail;

static Blended blendStage(const Request& request, bool continuation);

//! Carry edge intent only along modifications of that edge, never along generated
//! faces or unrelated generated edges. Retention is resolved on the final topology.
std::vector<Join> carriedJoins(const std::vector<Join>& joins,
                               const Handle(BRepTools_History)& history,
                               const TopoDS_Shape& result)
{
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> retained;
  TopExp::MapShapes(result, TopAbs_EDGE, retained);
  std::vector<Join> out;
  for (const auto& join : joins)
  {
    const auto& modified = history->Modified(join.edge);
    if (modified.IsEmpty())
    {
      if (retained.Contains(join.edge))
        out.push_back(join);
    }
    else
      for (const auto& edge : modified)
        if (edge.ShapeType() == TopAbs_EDGE && retained.Contains(edge))
          out.push_back({edge, join.role, join.order});
  }
  return out;
}

//! How the faces beside an edge turn along it (turnAlong): -1 where the edge is convex, 1 where it
//! is concave, 0 where they are tangent for a blend of `size` all along; `touches` tells whether
//! they are tangent somewhere (Rail::touching).
int convexity(const char*        op,
              const TopoDS_Edge& edge,
              const TopoDS_Face& first,
              const TopoDS_Face& second,
              double             size,
              bool*              touches)
{
  const Turn turn = turnAlong(op, edge, first, second, size);
  if (touches)
    *touches = turn.touches;
  return turn.kind;
}


//! The edges asked, with those tangent to them added: a blend goes on to the next edge where it
//! meets the previous one smoothly (tangent for a blend of `size`, tolerance::tangentFor). An edge
//! that goes on straight between
//! other faces is a sharp continuation, not added. Where `mayAdd` is not empty, only edges in it
//! are added.
std::vector<TopoDS_Shape> withTangent(const TopoDS_Shape&              body,
                                      std::vector<TopoDS_Shape>        edges,
                                      double                           size,
                                      const std::vector<TopoDS_Shape>& mayAdd)
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
  out.joins = carriedJoins(made.joins, then.history, then.shape);
  out.joins.insert(out.joins.end(), then.joins.begin(), then.joins.end());
  out.history = made.history;
  out.history->Merge(then.history);
  for (const auto& [before, middle] : made.unchanged)
    for (const auto& [again, after] : then.unchanged)
      if (again.IsSame(middle))
        out.unchanged.emplace_back(before, after);
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
    const BoundaryJet       jet("blendEdges", inCap, cap);
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
    out.joins = made.joins;

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
      again.fitWork = request.fitWork;
      again.reads = request.reads;
      again.refineContacts = request.refineContacts;
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
      out.joins = carriedJoins(out.joins, then.history, then.shape);
      out.joins.insert(out.joins.end(), then.joins.begin(), then.joins.end());
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
    out.shape = shape;
    for (const auto& [before, after] : kept)
      out.unchanged.emplace_back(before, after);
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
      join.operation       = op;
      join.select          = Select::Union;
      join.sources         = {solid};
      join.tools           = {prism.Shape()};
      join.subShapeHistory = true;
      join.cancellable     = true;
      const auto extended  = booleanCells(join);
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
  cut.subShapeHistory         = true;
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
    back.operation       = op;
    back.select          = kind < 0 && !lifted->continuedBoss ? Select::Difference : Select::Union;
    back.sources         = {made.shape};
    back.tools           = lifted->tools;
    back.subShapeHistory = true; // its history, for the faces' history through it
    back.cancellable     = true;
    const BooleanResult put = booleanCells(back);
    if (put.shape.IsNull() || singleSolid(put.shape).IsNull())
      return {};
    Blended out;
    out.shape   = singleSolid(put.shape);
    out.joins = carriedJoins(made.joins, put.history, out.shape);
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
  if (request.fitWork != FitWork::Deferred)
    return blendStage(request, true);
  // Each stage after the first changes one input of the request. Where the stage before it read
  // none of that input (StageReads), it is that stage again, step by step, and fails as it did.
  StageReads reads;
  Request    bounded = request;
  bounded.fitWork    = FitWork::Bounded;
  bounded.reads      = &reads;
  // The first stage's refusal is the one reported, unless it names a construction that failed and
  // a later stage measured the obstruction on its sheet (FailureKind::Curvature: the sheet crosses
  // itself, a contact folds): then the last such measurement is reported.
  std::exception_ptr original, measured;
  bool               obstruction = false;
  try
  {
    return blendStage(bounded, true);
  }
  catch (const Failure& failure)
  {
    if (failure.reason != no_solution && failure.reason != unsupported_geometry)
      throw;
    original    = std::current_exception();
    obstruction = failure.kind == FailureKind::Curvature;
  }
  catch (const Standard_Failure&)
  {
    original = std::current_exception();
  }
  Request extended = request;
  extended.fitWork = FitWork::Extended;
  extended.reads   = &reads;
  const auto stage = [&](bool run) -> std::optional<Blended> {
    if (!run)
      return std::nullopt;
    reads = {};
    try
    {
      return blendStage(extended, true);
    }
    catch (const Failure& failure)
    {
      if (failure.reason == std::string("cancelled"))
        throw;
      if (failure.kind == FailureKind::Curvature)
        measured = std::current_exception();
    }
    catch (const Standard_Failure&)
    {
    }
    return std::nullopt;
  };
  if (auto made = stage(reads.fitBound))
    return *made;
  extended.refineContacts = true;
  if (auto made = stage(reads.refines))
    return *made;
  // A place whose spine cannot be traced (it throws) is one the stage may fail at differently.
  const auto changes = [](const std::function<bool()>& at) {
    try
    {
      return at();
    }
    catch (const Failure& failure)
    {
      if (failure.reason == std::string("cancelled"))
        throw;
    }
    catch (const Standard_Failure&)
    {
    }
    return true;
  };
  extended.intersectSpines = true;
  if (auto made = stage(std::any_of(reads.spines.begin(), reads.spines.end(), changes)))
    return *made;
  std::rethrow_exception(obstruction || !measured ? original : measured);
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

} // namespace BRepExt::blending
