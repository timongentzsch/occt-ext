// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepBlendExt_Attach.hxx>
#include <BRepBlendExt_EndBoundary.hxx>
#include <BRepBlendExt_BoundaryCurve.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_TVertex.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>

namespace BRepExt::blending::detail
{
namespace
{
Station boundarySeed(const Piece& piece, double u)
{
  Station seed;
  for (int side = 0; side < 2; ++side)
  {
    const auto uv = piece.fitted.contact[side]->Value(u);
    seed.uv[side][0] = uv.X();
    seed.uv[side][1] = uv.Y();
  }
  return seed;
}
}

void Attach::stopPlaces(int k)
{
  const Piece& piece = pieces[k];
  const Rail& rail = piece.rail;
  // Continuation currently describes regular constant conic and G2 sections.
  if ((rail.law != Law::Conic && rail.law != Law::G2) || !rail.sizes.places.empty())
    return;
  for (const auto& side : rail.s)
    if (side.cliffCurve || !side.trim.IsNull())
      return;
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> visited;
  for (int end = 0; end < 2; ++end)
  {
    if (piece.joint[end].meets != Meets::Stop)
      continue;
    const auto edge = TopoDS::Edge(rail.edge.Oriented(TopAbs_FORWARD));
    const auto vertex = end ? TopExp::LastVertex(edge) : TopExp::FirstVertex(edge);
    if (!edgesOfVertex.Contains(vertex))
      continue;
    for (const auto& incident : edgesOfVertex.FindFromKey(vertex))
      for (const auto& shape : facesOfEdge.FindFromKey(incident))
      {
        if (shape.IsSame(rail.s[0].face) || shape.IsSame(rail.s[1].face) || !visited.Add(shape))
          continue;
        const auto face = TopoDS::Face(shape);
        const auto plane = Handle(Geom_Plane)::DownCast(placedFace(face).surface);
        if (plane.IsNull())
          continue;
        for (int side = 0; side < 2; ++side)
        {
          GeomAPI_IntCS intersection(piece.sheet.contact(side), plane);
          for (int i = 1; intersection.IsDone() && i <= intersection.NbPoints(); ++i)
          {
            double a, b, t;
            intersection.Parameters(i, a, b, t);
            const double u = (t - piece.sheet.alpha) / piece.sheet.beta;
            if (u < 0 || u > 1)
              continue;
            const auto trace = tracePlaneEnd(rail, plane->Pln(), boundarySeed(piece, u),
                                              side, rail.scale / 20);
            // Incidence is determined by the exact section equations, not by a missing
            // opposite corner. An ordinary opposite-contact component keeps its curves.
            if (trace.end == EndBoundaryTrace::End::Contact && trace.points.size() > 1
                && trace.points.back().parameters[4] == side)
            {
              tracedBoundary(k);
              return;
            }
          }
        }
      }
  }
}

void Attach::tracedBoundary(int k)
{
  const auto completed = completeCarrier(pieces[k]);
  if (!completed)
    fail(op, no_solution, "the retained end component has no complete regular carrier", {pieces[k].rail.edge});
  pieces[k] = *completed;
  const Piece& piece = pieces[k];
  const Rail& rail = piece.rail;
  const Sheet& sheet = piece.sheet;
  struct Port
  {
    gp_Pnt2d uv;
    TopoDS_Vertex vertex;
  };
  struct Curve
  {
    TopoDS_Edge edge;
    size_t from, to;
    TopoDS_Face owner;
  };
  std::vector<Port> ports;
  std::vector<Curve> curves;
  for (TopExp_Explorer f(body, TopAbs_FACE); f.More(); f.Next())
  {
    const auto face = TopoDS::Face(f.Current());
    if (face.IsSame(rail.s[0].face) || face.IsSame(rail.s[1].face))
      continue;
    const auto plane = Handle(Geom_Plane)::DownCast(placedFace(face).surface);
    if (plane.IsNull())
      continue;
    // Roots are labelled by their contact and their place on its retained carrier.
    struct Root { double u; int side; bool used = false; };
    std::vector<Root> roots;
    for (int side = 0; side < 2; ++side)
    {
      GeomAPI_IntCS meet(sheet.contact(side), plane);
      if (!meet.IsDone())
        fail(op, no_solution, "a contact-plane intersection is unresolved", {rail.edge});
      for (int i = 1; i <= meet.NbPoints(); ++i)
      {
        double a, b, t;
        meet.Parameters(i, a, b, t);
        const double u = (t - sheet.alpha) / sheet.beta;
        if (u >= 0 && u <= 1)
          roots.push_back({u, side});
      }
    }
    for (size_t index = 0; index < roots.size(); ++index)
    {
      auto& root = roots[index];
      if (root.used)
        continue;
      const auto trace = tracePlaneEnd(rail, plane->Pln(), boundarySeed(piece, root.u),
                                        root.side, rail.scale / 20);
      if (trace.end != EndBoundaryTrace::End::Contact || trace.points.size() < 2)
        fail(op, no_solution, "a retained planar end has unresolved terminal incidence", {rail.edge});
      std::vector<gp_Pnt2d> guide;
      for (const auto& point : trace.points)
      {
        const double w = std::clamp(point.parameters[4], 0., 1.);
        GeomAPI_ProjectPointOnCurve projection(point.point, sheet.surface->VIso(w));
        if (!projection.NbPoints())
          fail(op, no_solution, "a section-plane component misses its fitted sheet", {rail.edge});
        guide.emplace_back(projection.LowerDistanceParameter(), w);
      }
      const int lastSide = int(trace.points.back().parameters[4]);
      size_t mate = roots.size();
      for (size_t j = 0; j < roots.size(); ++j)
        if (j != index && roots[j].side == lastSide
            && sheet.contact(lastSide)->Value(sheet.rail(roots[j].u)).Distance(trace.points.back().point)
                 <= tolerance)
        {
          if (mate != roots.size())
            fail(op, no_solution, "a contact event has ambiguous root incidence", {rail.edge});
          mate = j;
        }
      if (mate == roots.size() || roots[mate].used)
        fail(op, no_solution, "a planar end has no matching terminal contact root", {rail.edge});
      root.used = roots[mate].used = true;
      const auto first = BRepBuilderAPI_MakeVertex(trace.points.front().point).Vertex();
      const auto last = BRepBuilderAPI_MakeVertex(trace.points.back().point).Vertex();
      const auto fitted = planeBoundaryCurve(sheet, face, guide, first, last);
      if (!fitted)
        fail(op, no_solution, "a retained planar boundary could not be enclosed at length tolerance", {rail.edge});
      const size_t at = ports.size();
      ports.push_back({guide.front(), first});
      ports.push_back({guide.back(), last});
      curves.push_back({fitted->edge, at, at + 1, face});
    }
  }
  for (int side = 0; side < 2; ++side)
  {
    contactIntervals[k][side].clear();
    std::vector<size_t> events;
    for (size_t i = 0; i < ports.size(); ++i)
      if (ports[i].uv.Y() == side)
        events.push_back(i);
    std::sort(events.begin(), events.end(), [&](size_t a, size_t b) { return ports[a].uv.X() < ports[b].uv.X(); });
    for (size_t i = 1; i < events.size(); ++i)
    {
      const size_t a = events[i - 1], b = events[i];
      const double middle = (ports[a].uv.X() + ports[b].uv.X()) / 2;
      const auto point = sheet.surface->Value(middle, side);
      GeomAPI_ProjectPointOnSurf projection(point, rail.s[side].surface);
      double u, v;
      projection.LowerDistanceParameters(u, v);
      BRepClass_FaceClassifier state(rail.s[side].face, gp_Pnt2d(u, v), tolerance);
      if (state.State() != TopAbs_IN && state.State() != TopAbs_ON)
        continue;
      ContactInterval interval;
      interval.ends = {{{ports[a].vertex, ports[a].uv.X()}, {ports[b].vertex, ports[b].uv.X()}}};
      const auto support = placedFace(rail.s[side].face);
      interval.edge = edgeBetween(sheet.contact(side), ports[a].vertex, ports[a].uv.X(),
                                   ports[b].vertex, ports[b].uv.X(),
                                   {{sheet.contactUV(side), sheet.surface, {}},
                                    {Handle(Geom2d_Curve)(piece.fitted.contact[side]), support.basis, support.location}}, tolerance);
      contactIntervals[k][side].push_back(interval);
      curves.push_back({interval.edge, a, b, rail.s[side].face});
    }
  }
  if (curves.empty())
    fail(op, no_solution, "the retained boundary has no arcs", {rail.edge});
  BoundaryGraph graph;
  std::vector<bool> used(curves.size(), false);
  for (size_t first = 0; first < curves.size(); ++first)
  {
    if (used[first])
      continue;
    size_t current = first, start = curves[first].from, at = start;
    // The halfspace kept by the original face orients the first arc of each loop.
    const auto& arc = curves[first];
    const auto plane = Handle(Geom_Plane)::DownCast(placedFace(arc.owner).surface);
    if (!plane.IsNull())
    {
      double f, l;
      const auto pc = BRep_Tool::CurveOnSurface(arc.edge, sheet.surface, {}, f, l);
      const auto uv = pc->Value((f + l) / 2);
      gp_Pnt p; gp_Vec du, dv, tangent, pu, pv;
      sheet.surface->D1(uv.X(), uv.Y(), p, du, dv);
      BRepAdaptor_Curve(arc.edge).D1((f + l) / 2, p, tangent);
      plane->D1(0, 0, p, pu, pv);
      auto normal = pu.Crossed(pv);
      if (arc.owner.Orientation() == TopAbs_REVERSED)
        normal.Reverse();
      if (normal.Dot(du.Crossed(dv).Crossed(tangent)) > 0)
        start = at = arc.to;
    }
    do
    {
      const auto& arc = curves[current];
      if (used[current])
        fail(op, no_solution, "the end boundary repeats before closing", {rail.edge});
      used[current] = true;
      const bool forward = arc.from == at;
      const size_t to = forward ? arc.to : arc.from;
      graph.arcs.push_back({forward ? arc.edge : TopoDS::Edge(arc.edge.Reversed()), at, to});
      at = to;
      if (at == start)
        break;
      size_t next = curves.size();
      for (size_t i = 0; i < curves.size(); ++i)
        if (!used[i] && (curves[i].from == at || curves[i].to == at))
        {
          if (next != curves.size())
            fail(op, no_solution, "the end boundary has branching incidence", {rail.edge});
          next = i;
        }
      if (next == curves.size())
        fail(op, no_solution, "the end boundary has an open event", {rail.edge});
      current = next;
    } while (true);
  }
  retainedBoundaries[k] = std::move(graph);
}

namespace
{
TopoDS_Edge boundaryPart(const TopoDS_Edge& original, const TopoDS_Vertex& a, double first,
                         const TopoDS_Vertex& b, double last)
{
  BRep_Builder builder;
  auto edge = TopoDS::Edge(original.Oriented(TopAbs_FORWARD).EmptyCopied());
  builder.Add(edge, a.Oriented(TopAbs_FORWARD));
  builder.Add(edge, b.Oriented(TopAbs_REVERSED));
  builder.Range(edge, first, last);
  for (const auto& event : {std::pair{a, first}, std::pair{b, last}})
    if (Handle(BRep_TVertex)::DownCast(event.first.TShape())->Points().IsEmpty())
      builder.UpdateVertex(event.first, event.second, edge, tol::length);
  return edge;
}
}

void Attach::installBoundaries()
{
  if (std::none_of(retainedBoundaries.begin(), retainedBoundaries.end(), [](const auto& b) { return bool(b); }))
    return;
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> oldEdges, oldFaces;
  TopExp::MapShapes(body, TopAbs_EDGE, oldEdges);
  TopExp::MapShapes(body, TopAbs_FACE, oldFaces);
  struct Cut { double t; TopoDS_Vertex vertex; };
  std::map<int, std::vector<Cut>> eventsOnEdge;
  const auto addCut = [&](int index, double t, TopoDS_Vertex vertex) {
    auto& list = eventsOnEdge[index];
    for (const auto& event : list)
      if (event.vertex.IsSame(vertex))
        return event.vertex;
    for (const auto& event : list)
      if (BRep_Tool::Pnt(event.vertex).Distance(BRep_Tool::Pnt(vertex)) <= tolerance)
        return event.vertex;
    list.push_back({t, vertex});
    return vertex;
  };
  struct Use { TopoDS_Edge edge; int owner, piece; };
  std::vector<Use> boundary;
  std::vector<bool> reversed(n, false);
  for (int k = 0; k < n; ++k)
  {
    if (!retainedBoundaries[k])
      continue;
    const auto& piece = pieces[k];
    gp_Pnt p; gp_Vec du, dv, outward;
    piece.sheet.surface->D1(.5, .5, p, du, dv);
    for (int side = 0; side < 2; ++side)
    {
      const auto uv = piece.fitted.contact[side]->Value(.5);
      Jet jet;
      if (!jetAt(piece.rail.s[side], uv.X(), uv.Y(), jet))
        fail(op, no_solution, "a retained contact has no regular support normal", {piece.rail.edge});
      outward += jet.n;
      contactIntervals[k][side].clear();
    }
    reversed[k] = du.Crossed(dv).Dot(outward) < 0;
    BoundaryGraph next;
    for (const auto& arc : retainedBoundaries[k]->arcs)
    {
      const auto edge = arc.edge;
      int owner = 0;
      Handle(Geom2d_Curve) pc;
      double first, last;
      for (int fi = 1; fi <= oldFaces.Extent(); ++fi)
      {
        bool stored = false;
        const auto curve = BRep_Tool::CurveOnSurface(edge, TopoDS::Face(oldFaces(fi)), first, last, &stored);
        if (!stored)
          continue;
        if (owner)
          fail(op, no_solution, "a boundary arc has multiple original face owners", {edge});
        owner = fi;
        pc = curve;
      }
      if (!owner)
        fail(op, no_solution, "a boundary arc has no original face owner", {edge});
      const auto face = TopoDS::Face(oldFaces(owner));
      const auto surface = placedFace(face).surface;
      double ef, el;
      const auto curve = BRep_Tool::Curve(edge, ef, el);
      const auto forwardEdge = TopoDS::Edge(edge.Oriented(TopAbs_FORWARD));
      std::vector<Cut> events = {{ef, TopExp::FirstVertex(forwardEdge)}, {el, TopExp::LastVertex(forwardEdge)}};
      NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> seams;
      for (TopExp_Explorer x(face, TopAbs_EDGE); x.More(); x.Next())
      {
        const auto seam = TopoDS::Edge(x.Current());
        if (!BRep_Tool::IsClosed(seam, face) || !seams.Add(seam))
          continue;
        double sf, sl;
        const auto sc = BRep_Tool::CurveOnSurface(seam, face, sf, sl);
        if (sc.IsNull())
          continue;
        for (const auto& meeting : meetings(pc, sc, surface))
        {
          if (!(meeting.first > ef + tol::parameter && meeting.first < el - tol::parameter
                && meeting.second >= sf && meeting.second <= sl))
            continue;
          double a, b;
          const auto original = BRep_Tool::Curve(seam, a, b);
          const auto point = curve->Value(meeting.first);
          if (point.Distance(original->Value(meeting.second)) > tolerance)
            fail(op, no_solution, "a periodic event misses its original seam", {edge, seam});
          auto vertex = BRepBuilderAPI_MakeVertex(point).Vertex();
          vertex = addCut(oldEdges.FindIndex(seam), meeting.second, vertex);
          if (std::none_of(events.begin(), events.end(), [&](const Cut& e) { return e.vertex.IsSame(vertex); }))
            events.push_back({meeting.first, vertex});
        }
      }
      std::sort(events.begin(), events.end(), [](const Cut& a, const Cut& b) { return a.t < b.t; });
      std::vector<TopoDS_Edge> parts;
      for (size_t i = 1; i < events.size(); ++i)
      {
        auto part = boundaryPart(edge, events[i - 1].vertex, events[i - 1].t, events[i].vertex, events[i].t);
        double a, b;
        const auto on = BRep_Tool::CurveOnSurface(part, face, a, b);
        const double middle = (a + b) / 2;
        GeomAPI_ProjectPointOnSurf projection(curve->Value(middle), surface);
        double u, v;
        projection.LowerDistanceParameters(u, v);
        builder.UpdateEdge(part, aligned(on, middle, gp_Pnt2d(u, v), surface), face, tolerance);
        parts.push_back(part);
        int side = -1;
        for (int s = 0; s < 2; ++s)
          if (face.IsSame(piece.rail.s[s].face))
            side = s;
        if (side >= 0)
        {
          ContactInterval interval;
          interval.ends = {{{events[i - 1].vertex, events[i - 1].t}, {events[i].vertex, events[i].t}}};
          interval.edge = part;
          contactIntervals[k][side].push_back(interval);
        }
        else
          notches.push_back({k, part, face});
      }
      if (edge.Orientation() == TopAbs_REVERSED)
      {
        std::reverse(parts.begin(), parts.end());
        for (auto& part : parts)
          part.Reverse();
      }
      size_t from = arc.from;
      for (size_t i = 0; i < parts.size(); ++i)
      {
        // New seam ports are distinct from every pre-existing event port.
        const size_t to = i + 1 == parts.size() ? arc.to
            : retainedBoundaries[k]->arcs.size() * 2 + next.arcs.size();
        next.arcs.push_back({parts[i], from, to});
        boundary.push_back({parts[i], owner, k});
        from = to;
      }
    }
    for (int side = 0; side < 2; ++side)
      std::sort(contactIntervals[k][side].begin(), contactIntervals[k][side].end(),
                [](const auto& a, const auto& b) { return a.ends[0].parameter < b.ends[0].parameter; });
    retainedBoundaries[k] = std::move(next);
  }
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> vertices;
  for (const auto& use : boundary)
    TopExp::MapShapes(use.edge, TopAbs_VERTEX, vertices);
  for (const auto& shape : vertices)
  {
    const auto vertex = TopoDS::Vertex(shape);
    bool found = false;
    for (int ei = 1; ei <= oldEdges.Extent(); ++ei)
    {
      const auto edge = TopoDS::Edge(oldEdges(ei));
      if (BRep_Tool::Degenerated(edge))
        continue;
      double f, l;
      const auto curve = BRep_Tool::Curve(edge, f, l);
      if (curve.IsNull())
        continue;
      const auto at = onCurve(BRep_Tool::Pnt(vertex), new Geom_TrimmedCurve(curve, f, l));
      if (!at || at->distance > tolerance)
        continue;
      if (!addCut(ei, at->parameter, vertex).IsSame(vertex))
        fail(op, no_solution, "coincident boundary events do not share their vertex", {edge});
      found = true;
    }
    if (!found)
      fail(op, no_solution, "a boundary event misses the original face wires", {vertex});
  }
  std::map<int, std::vector<TopoDS_Edge>> split;
  for (int ei = 1; ei <= oldEdges.Extent(); ++ei)
  {
    const auto edge = TopoDS::Edge(oldEdges(ei).Oriented(TopAbs_FORWARD));
    double f, l;
    BRep_Tool::Range(edge, f, l);
    auto events = eventsOnEdge[ei];
    if (events.empty())
    {
      split[ei] = {edge};
      continue;
    }
    events.push_back({f, TopExp::FirstVertex(edge)});
    events.push_back({l, TopExp::LastVertex(edge)});
    std::sort(events.begin(), events.end(), [](const Cut& a, const Cut& b) { return a.t < b.t; });
    for (size_t i = 1; i < events.size(); ++i)
    {
      if (events[i].t == events[i - 1].t)
        fail(op, no_solution, "an original edge has duplicate boundary events", {edge});
      const auto part = boundaryPart(edge, events[i - 1].vertex, events[i - 1].t, events[i].vertex, events[i].t);
      split[ei].push_back(part);
      boundaryHistory.push_back({edge, part});
    }
  }
  for (int fi = 1; fi <= oldFaces.Extent(); ++fi)
  {
    const auto face = TopoDS::Face(oldFaces(fi));
    bool changed = std::any_of(boundary.begin(), boundary.end(), [&](const Use& use) { return use.owner == fi; });
    for (TopExp_Explorer e(face, TopAbs_EDGE); !changed && e.More(); e.Next())
      changed = !eventsOnEdge[oldEdges.FindIndex(e.Current())].empty();
    if (!changed)
      continue;
    const auto forwardFace = TopoDS::Face(face.Oriented(TopAbs_FORWARD));
    const auto surface = placedFace(face).surface;
    struct Node { TopoDS_Vertex vertex; gp_Pnt2d uv; };
    struct Arc { TopoDS_Edge edge; size_t from, to; bool old, kept; };
    std::vector<Node> nodes;
    std::vector<Arc> arcs;
    for (TopExp_Explorer wi(forwardFace, TopAbs_WIRE); wi.More(); wi.Next())
    {
      std::vector<TopoDS_Edge> uses;
      for (BRepTools_WireExplorer walk(TopoDS::Wire(wi.Current()), forwardFace); walk.More(); walk.Next())
        uses.push_back(walk.Current());
      const size_t base = nodes.size();
      for (const auto& use : uses)
      {
        double f, l;
        const auto pc = BRep_Tool::CurveOnSurface(use, forwardFace, f, l);
        nodes.push_back({TopExp::FirstVertex(use, true), pc->Value(use.Orientation() == TopAbs_FORWARD ? f : l)});
      }
      for (size_t i = 0; i < uses.size(); ++i)
      {
        const auto use = uses[i];
        double f, l;
        const auto pc = BRep_Tool::CurveOnSurface(use, forwardFace, f, l);
        auto parts = split.at(oldEdges.FindIndex(use));
        if (use.Orientation() == TopAbs_REVERSED)
        {
          std::reverse(parts.begin(), parts.end());
          for (auto& part : parts)
            part.Reverse();
        }
        size_t at = base + i;
        bool removed = false;
        for (int k = 0; k < n; ++k)
          removed = removed || (retainedBoundaries[k] && ownEdge(pieces[k], use));
        for (size_t j = 0; j < parts.size(); ++j)
        {
          size_t to = base + (i + 1) % uses.size();
          if (j + 1 < parts.size())
          {
            double a, b;
            BRep_Tool::Range(parts[j], a, b);
            to = nodes.size();
            nodes.push_back({TopExp::LastVertex(parts[j], true),
                             pc->Value(parts[j].Orientation() == TopAbs_FORWARD ? b : a)});
          }
          arcs.push_back({parts[j], at, to, true, !removed});
          at = to;
        }
      }
    }
    const auto port = [&](const TopoDS_Vertex& vertex, const gp_Pnt2d& uv) {
      size_t found = nodes.size();
      double distance = Precision::Infinite();
      for (size_t i = 0; i < nodes.size(); ++i)
        if (nodes[i].vertex.IsSame(vertex) && uv.Distance(nodes[i].uv) < distance)
        {
          distance = uv.Distance(nodes[i].uv);
          found = i;
        }
      if (found == nodes.size())
        fail(op, no_solution, "a shared boundary vertex has no original parameter port", {vertex});
      gp_Pnt p; gp_Vec du, dv;
      surface->D1(uv.X(), uv.Y(), p, du, dv);
      const auto delta = uv.XY() - nodes[found].uv.XY();
      if ((du * delta.X() + dv * delta.Y()).Magnitude() > tolerance)
        fail(op, no_solution, "a boundary event misses its original parameter port", {vertex});
      return found;
    };
    for (const auto& use : boundary)
      if (use.owner == fi)
      {
        auto edge = use.edge;
        if (!reversed[use.piece])
          edge.Reverse();
        if (face.Orientation() == TopAbs_REVERSED)
          edge.Reverse();
        double f, l;
        const auto pc = BRep_Tool::CurveOnSurface(edge, forwardFace, f, l);
        const bool forward = edge.Orientation() == TopAbs_FORWARD;
        const size_t a = port(TopExp::FirstVertex(edge, true), pc->Value(forward ? f : l));
        const size_t b = port(TopExp::LastVertex(edge, true), pc->Value(forward ? l : f));
        arcs.push_back({edge, a, b, false, true});
      }
    for (auto& arc : arcs)
      if (arc.old)
        for (const auto& bridge : arcs)
          if (!bridge.old && (arc.from == bridge.from || arc.to == bridge.to))
            arc.kept = false;
    bool progress = true;
    while (progress)
    {
      progress = false;
      std::vector<int> incoming(nodes.size()), outgoing(nodes.size());
      for (const auto& arc : arcs)
        if (arc.kept)
        {
          ++outgoing[arc.from];
          ++incoming[arc.to];
        }
      for (auto& arc : arcs)
        if (arc.kept && arc.old && (!incoming[arc.from] || !outgoing[arc.to]))
        {
          arc.kept = false;
          progress = true;
        }
    }
    BoundaryGraph graph;
    for (const auto& arc : arcs)
      if (arc.kept)
        graph.arcs.push_back({arc.edge, arc.from, arc.to});
    const auto loops = graph.loops();
    if (!loops)
      fail(op, no_solution, "the retained original face has unresolved boundary incidence", {face});
    auto out = TopoDS::Face(face.EmptyCopied());
    for (const auto& loop : *loops)
    {
      TopoDS_Wire wire;
      builder.MakeWire(wire);
      for (const auto& edge : loop)
        builder.Add(wire, edge);
      wire.Closed(true);
      builder.Add(out, wire);
    }
    editOf(edits, face);
    retainedFaces.push_back({face, {out}});
  }
}
} // namespace BRepExt::blending::detail
