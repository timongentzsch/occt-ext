// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepBlendExt_Attach.hxx>
#include <BRepBlendExt_EndBoundary.hxx>
#include <BRepBlendExt_BoundaryCurve.hxx>
#include <GeomAlgoExt_Timing.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_TVertex.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>

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
  BooleanPhase phase("fillet2.endIncidence");
  if (retainedBoundaries[k])
    return;
  const Piece& piece = pieces[k];
  const Rail& rail = piece.rail;
  // Continuation currently describes regular constant conic and G2 sections.
  if ((rail.law != Law::Conic && rail.law != Law::G2) || !rail.sizes.places.empty())
    return;
  for (const auto& side : rail.s)
    if (side.cliffCurve || !side.trim.IsNull())
      return;
  for (int end = 0; end < 2; ++end)
  {
    if (piece.joint[end].meets != Meets::Stop)
      continue;
    const auto edge = TopoDS::Edge(rail.edge.Oriented(TopAbs_FORWARD));
    auto vertex = end ? TopExp::LastVertex(edge) : TopExp::FirstVertex(edge);
    if (!piece.stopVertex[end].IsNull())
      vertex = piece.stopVertex[end];
    if (!edgesOfVertex.Contains(vertex))
      continue;
    OneStop source;
    source.from[0] = source.from[1] = vertex;
    source.middle = piece.uAt(0) / 2 + piece.uAt(1) / 2;
    for (const auto& incident : edgesOfVertex.FindFromKey(vertex))
      if (!ownEdge(piece, incident))
        for (int side = 0; side < 2; ++side)
          if (onFace(incident, rail.s[side].face))
            source.beside[side] = TopoDS::Edge(incident);
    if (source.beside[0].IsNull() || source.beside[1].IsNull())
      continue;
    const auto a = faceAcross(source.beside[0], rail.s[0].face);
    const auto b = faceAcross(source.beside[1], rail.s[1].face);
    if (!a.IsNull() && a.IsSame(b))
    {
      source.stop = a;
      // These are the same source events subsequently used to construct the
      // ordinary end's vertices, curves and ordered graph occurrences.
      if (simpleStop(k, end, source))
        continue;
    }
    StopEnd plan{k, end, vertex};
    plan.middle = source.middle;
    plan.carried[0] = piece.sheet.rail(end == 0 ? 0 : plan.middle);
    plan.carried[1] = piece.sheet.rail(end == 0 ? plan.middle : 1);
    for (int side = 0; side < 2; ++side)
      walkSupport(plan, side);
    stopWalks.insert_or_assign(std::pair{k, end}, plan);
    // The source walk owns the finite-face component. An open contact can be
    // paired with the other contact, or the section can return to that same
    // contact. Its exact trace distinguishes those incidences before attachment.
    const int open = plan.walk[0].open.IsNull() ? 1 : 0, side = 1 - open;
    if (plan.walk[open].open.IsNull() || !plan.walk[side].open.IsNull()
        || plan.walk[side].crossed.IsNull() || plan.walk[0].continued
        || plan.walk[1].continued || !plan.walk[side].level.IsNull())
      continue;
    const auto face = plan.walk[side].faces.back();
    const auto plane = Handle(Geom_Plane)::DownCast(placedFace(face).surface);
    if (plane.IsNull())
      continue;
    const auto trace = tracePlaneEnd(rail, plane->Pln(),
                                      boundarySeed(piece, plan.walk[side].at.first),
                                      side, rail.scale / 20);
    if (trace.end == EndBoundaryTrace::End::Contact && trace.points.size() > 1
        && trace.contactSide == side)
    {
      tracedBoundary(k);
      return;
    }
  }
}

void Attach::tracedBoundary(int k)
{
  BooleanPhase phase("fillet2.endBoundary");
  const auto   component = [&] {
    BooleanPhase carrier("fillet2.endCarrier");
    return carrierComponent(pieces[k], body);
  }();
  if (!component)
    fail(op, no_solution, "the retained end component has no complete regular carrier", {pieces[k].rail.edge});
  // The completed component replaces the connected source chain, including
  // chart splits already found by the section stage. Its branches can revisit a
  // source face, so chart ownership cannot be recovered from a face label alone.
  std::vector<int> indices = {k};
  for (size_t i = 0; i < indices.size(); ++i)
    for (const auto& joint : pieces[indices[i]].joint)
      if ((joint.meets == Meets::Joint || joint.meets == Meets::Split)
          && std::find(indices.begin(), indices.end(), joint.piece) == indices.end())
        indices.push_back(joint.piece);
  std::vector<TopoDS_Edge> owned;
  for (int index : indices)
  {
    const auto retain = [&](const TopoDS_Edge& edge) {
      if (std::none_of(owned.begin(), owned.end(), [&](const auto& e) { return e.IsSame(edge); }))
        owned.push_back(edge);
    };
    retain(pieces[index].rail.edge);
    for (const auto& edge : pieces[index].swallowed)
      retain(edge);
    stopWalks.erase({index, 0});
    stopWalks.erase({index, 1});
  }
  contactCrossings.erase(
    std::remove_if(contactCrossings.begin(),
                   contactCrossings.end(),
                   [&](const auto& event) {
                     return std::find(indices.begin(), indices.end(), event.piece) != indices.end();
                   }),
    contactCrossings.end());
  while (indices.size() > component->charts.size())
  {
    const int unused = indices.back();
    indices.pop_back();
    retainedBoundaries[unused] = BoundaryGraph{};
    for (auto& side : contactIntervals[unused])
      side.clear();
    pieces[unused].joint[0] = {};
    pieces[unused].joint[1] = {};
    pieces[unused].swallowed.clear();
  }
  while (indices.size() < component->charts.size())
  {
    indices.push_back(int(pieces.size()));
    pieces.emplace_back();
  }
  for (size_t i = 0; i < component->charts.size(); ++i)
  {
    auto& piece    = pieces[indices[i]];
    piece          = component->charts[i].piece;
    piece.joint[0] = {};
    piece.joint[1] = {};
    piece.swallowed.clear();
    if (i == 0)
      for (const auto& edge : owned)
        if (!edge.IsSame(piece.rail.edge))
          piece.swallowed.push_back(edge);
  }
  n = int(pieces.size());
  retainedBoundaries.resize(n);
  contactIntervals.resize(n);
  section.resize(n);
  stopAt.resize(n);
  blends.resize(n);
  chain.resize(n);

  struct Port
  {
    gp_Pnt2d      uv;
    TopoDS_Vertex vertex;
  };

  struct Curve
  {
    TopoDS_Edge edge;
    size_t      from, to;
    TopoDS_Face owner;
    bool        forward;
  };

  struct Chart
  {
    std::vector<Port>  ports;
    std::vector<Curve> curves;
  };

  std::vector<Chart> charts(indices.size());
  const auto         port = [&](size_t c, const gp_Pnt2d& uv, const TopoDS_Vertex& vertex) {
    auto& ports = charts[c].ports;
    for (size_t i = 0; i < ports.size(); ++i)
      if (ports[i].vertex.IsSame(vertex))
        return i;
    ports.push_back({uv, vertex});
    return ports.size() - 1;
  };

  struct JoinEvent
  {
    double                w;
    TopoDS_Vertex         vertex;
    std::array<size_t, 2> ports;
  };

  std::vector<std::vector<JoinEvent>> joins(component->joins.size());
  const auto joinEvent = [&](size_t j, double w, const gp_Pnt& point) -> JoinEvent {
    auto& events = joins[j];
    for (const auto& event : events)
      if (BRep_Tool::Pnt(event.vertex).Distance(point) <= tolerance)
        return event;
    JoinEvent event{w, BRepBuilderAPI_MakeVertex(point).Vertex(), {}};
    for (int s = 0; s < 2; ++s)
    {
      const auto [c, end] = component->joins[j].ends[s];
      const auto& sheet   = pieces[indices[c]].sheet;
      event.ports[s]      = port(c, gp_Pnt2d(sheet.rail(end), w), event.vertex);
    }
    events.push_back(event);
    return event;
  };
  for (size_t j = 0; j < joins.size(); ++j)
  {
    const auto [c, end] = component->joins[j].ends[0];
    for (int side = 0; side < 2; ++side)
      joinEvent(j, side, component->charts[c].ends[end].foot[side]);
  }
  const auto arcForward =
    [&](const Piece& piece, const TopoDS_Edge& edge, const TopoDS_Face& face) {
      const auto plane = Handle(Geom_Plane)::DownCast(placedFace(face).surface);
      double     f, l;
      const auto pc = BRep_Tool::CurveOnSurface(edge, piece.sheet.surface, {}, f, l);
      const auto uv = pc->Value((f + l) / 2);
      gp_Pnt     p;
      gp_Vec     du, dv, tangent, pu, pv;
      piece.sheet.surface->D1(uv.X(), uv.Y(), p, du, dv);
      BRepAdaptor_Curve(edge).D1((f + l) / 2, p, tangent);
      plane->D1(0, 0, p, pu, pv);
      auto normal = pu.Crossed(pv);
      if (face.Orientation() == TopAbs_REVERSED)
        normal.Reverse();
      return normal.Dot(du.Crossed(dv).Crossed(tangent)) <= 0;
    };
  for (TopExp_Explorer f(body, TopAbs_FACE); f.More(); f.Next())
  {
    const auto face = TopoDS::Face(f.Current());
    const auto plane = Handle(Geom_Plane)::DownCast(placedFace(face).surface);
    if (plane.IsNull())
      continue;

    struct Root
    {
      double u;
      int    side;
      bool   used = false;
    };

    std::vector<std::vector<Root>> roots(charts.size());
    for (size_t c = 0; c < charts.size(); ++c)
    {
      const auto& piece = pieces[indices[c]];
      if (face.IsSame(piece.rail.s[0].face) || face.IsSame(piece.rail.s[1].face))
        continue;
      for (int side = 0; side < 2; ++side)
      {
        GeomAPI_IntCS meet(piece.sheet.contact(side), plane);
        if (!meet.IsDone())
          fail(op, no_solution, "a contact-plane intersection is unresolved", {piece.rail.edge});
        for (int i = 1; i <= meet.NbPoints(); ++i)
        {
          double a, b, t;
          meet.Parameters(i, a, b, t);
          const double u = (t - piece.sheet.alpha) / piece.sheet.beta;
          if (u >= 0 && u <= 1)
            roots[c].push_back({u, side});
        }
      }
    }
    for (size_t startChart = 0; startChart < charts.size(); ++startChart)
      for (size_t ri = 0; ri < roots[startChart].size(); ++ri)
      {
        auto& root = roots[startChart][ri];
        if (root.used)
          continue;
        size_t                              c       = startChart;
        const auto&                         initial = pieces[indices[c]];
        Station                             seed    = boundarySeed(initial, root.u);
        double                              w       = root.side;
        gp_Vec                              direction;
        TopoDS_Vertex                       incoming;
        int                                 incomingEnd = -1;
        std::vector<std::pair<int, double>> visited;
        for (;;)
        {
          const auto& piece = pieces[indices[c]];
          const auto& rail  = piece.rail;
          const auto& sheet = piece.sheet;
          const auto trace = tracePlaneEnd(rail, plane->Pln(), seed, w, rail.scale / 20, direction);
          if (trace.points.size() < 2
              || (trace.end != EndBoundaryTrace::End::Contact
                  && trace.end != EndBoundaryTrace::End::SurfaceBoundary))
            fail(op,
                 no_solution,
                 "a retained planar end has unresolved terminal incidence",
                 {rail.edge});
          std::vector<gp_Pnt2d> guide;
          for (const auto& point : trace.points)
          {
            const double                across = std::clamp(point.parameters[4], 0., 1.);
            GeomAPI_ProjectPointOnCurve projection(point.point, sheet.surface->VIso(across));
            if (!projection.NbPoints())
              fail(op,
                   no_solution,
                   "a section-plane component misses its fitted sheet",
                   {rail.edge});
            guide.emplace_back(projection.LowerDistanceParameter(), across);
          }
          if (incomingEnd >= 0)
            guide.front().SetX(sheet.rail(incomingEnd));
          else
            guide.front().SetY(root.side);
          const auto    first = incoming.IsNull()
                                  ? BRepBuilderAPI_MakeVertex(trace.points.front().point).Vertex()
                                  : incoming;
          size_t        mate  = 0;
          int           outgoingEnd = -1, join = -1;
          TopoDS_Vertex last;
          if (trace.end == EndBoundaryTrace::End::Contact)
          {
            const int lastSide = trace.contactSide;
            guide.back().SetY(lastSide);
            mate = roots[c].size();
            for (size_t j = 0; j < roots[c].size(); ++j)
              if (!(c == startChart && j == ri) && roots[c][j].side == lastSide
                  && sheet.contact(lastSide)
                         ->Value(sheet.rail(roots[c][j].u))
                         .Distance(trace.points.back().point)
                       <= tolerance)
              {
                if (mate != roots[c].size())
                  fail(op,
                       no_solution,
                       "a contact event has ambiguous root incidence",
                       {rail.edge});
                mate = j;
              }
            if (mate == roots[c].size() || roots[c][mate].used)
              fail(op,
                   no_solution,
                   "a planar end has no matching terminal contact root",
                   {rail.edge});
            last = BRepBuilderAPI_MakeVertex(trace.points.back().point).Vertex();
          }
          else
          {
            for (int end = 0; end < 2; ++end)
              if (trace.points.back().centre.Distance(component->charts[c].ends[end].centre)
                  <= rail.resolved(limit::solved))
                outgoingEnd = end;
            if (outgoingEnd < 0 || component->charts[c].joins[outgoingEnd] < 0)
              fail(op,
                   no_solution,
                   "a planar boundary reaches an unjoined carrier restriction",
                   {rail.edge});
            join = component->charts[c].joins[outgoingEnd];
            guide.back().SetX(sheet.rail(outgoingEnd));
            last =
              joinEvent(join, trace.points.back().parameters[4], trace.points.back().point).vertex;
          }
          const auto fitted = [&] {
            BooleanPhase curve("fillet2.endCurve");
            return planeBoundaryCurve(sheet, face, guide, first, last);
          }();
          if (!fitted)
            fail(op,
                 no_solution,
                 "a retained planar boundary could not be enclosed at length tolerance",
                 {rail.edge});
          const size_t a = port(c, guide.front(), first), b = port(c, guide.back(), last);
          charts[c].curves.push_back(
            {fitted->edge, a, b, face, arcForward(piece, fitted->edge, face)});
          if (join < 0)
          {
            root.used = roots[c][mate].used = true;
            break;
          }
          const double nextW = trace.points.back().parameters[4];
          for (const auto& previous : visited)
            if (previous.first == join && previous.second == nextW)
              fail(op,
                   no_solution,
                   "a planar boundary repeats a carrier event before reaching a contact",
                   {rail.edge});
          visited.push_back({join, nextW});
          const auto& ends = component->joins[join].ends;
          const auto  next =
            ends[0].first == int(c) && ends[0].second == outgoingEnd ? ends[1] : ends[0];
          c           = next.first;
          incomingEnd = next.second;
          incoming    = last;
          seed        = component->charts[c].ends[incomingEnd];
          w           = nextW;
          direction =
            gp_Vec(trace.points[trace.points.size() - 2].point, trace.points.back().point);
        }
      }
  }
  for (size_t c = 0; c < charts.size(); ++c)
  {
    const auto& piece = pieces[indices[c]];
    const auto& rail  = piece.rail;
    const auto& sheet = piece.sheet;
    auto&       ports = charts[c].ports;
    for (int side = 0; side < 2; ++side)
    {
      contactIntervals[indices[c]][side].clear();
      std::vector<size_t> events;
      for (size_t i = 0; i < ports.size(); ++i)
        if (ports[i].uv.Y() == side)
          events.push_back(i);
      std::sort(events.begin(), events.end(), [&](size_t a, size_t b) {
        return ports[a].uv.X() < ports[b].uv.X();
      });
      for (size_t i = 1; i < events.size(); ++i)
      {
        const size_t               a = events[i - 1], b = events[i];
        const double               middle = (ports[a].uv.X() + ports[b].uv.X()) / 2;
        const auto                 point  = sheet.surface->Value(middle, side);
        GeomAPI_ProjectPointOnSurf projection(point, rail.s[side].surface);
        double                     u, v;
        projection.LowerDistanceParameters(u, v);
        BRepClass_FaceClassifier state(rail.s[side].face, gp_Pnt2d(u, v), tolerance);
        if (state.State() != TopAbs_IN && state.State() != TopAbs_ON)
          continue;
        ContactInterval interval;
        interval.ends = {{{ports[a].vertex, ports[a].uv.X()}, {ports[b].vertex, ports[b].uv.X()}}};
        const auto support = placedFace(rail.s[side].face);
        interval.edge      = edgeBetween(
          sheet.contact(side),
          ports[a].vertex,
          ports[a].uv.X(),
          ports[b].vertex,
          ports[b].uv.X(),
          {{sheet.contactUV(side), sheet.surface, {}},
           {Handle(Geom2d_Curve)(piece.fitted.contact[side]), support.basis, support.location}},
          tolerance);
        contactIntervals[indices[c]][side].push_back(interval);
        charts[c].curves.push_back({interval.edge, a, b, rail.s[side].face, side == 0});
      }
    }
  }
  const auto demand = [&](size_t c, size_t p) {
    int in = 0, out = 0;
    for (const auto& arc : charts[c].curves)
    {
      in += (arc.forward ? arc.to : arc.from) == p;
      out += (arc.forward ? arc.from : arc.to) == p;
    }
    return in == 1 && out == 0 ? 1 : in == 0 && out == 1 ? -1 : 0;
  };
  for (size_t j = 0; j < joins.size(); ++j)
  {
    auto events = joins[j];
    std::sort(events.begin(), events.end(), [](const auto& a, const auto& b) { return a.w < b.w; });
    const auto& ends = component->joins[j].ends;
    for (size_t i = 1; i < events.size(); ++i)
    {
      const auto& a = events[i - 1];
      const auto& b = events[i];
      bool        kept[2];
      for (int s = 0; s < 2; ++s)
      {
        const int sense = ends[s].second ? 1 : -1;
        kept[s] =
          demand(ends[s].first, a.ports[s]) == sense && demand(ends[s].first, b.ports[s]) == -sense;
      }
      if (kept[0] != kept[1])
        fail(op,
             no_solution,
             "adjacent carrier charts disagree on retained section incidence",
             {pieces[k].rail.edge});
      if (!kept[0])
        continue;
      std::array<Sheet, 2>  sheets;
      std::array<double, 2> parameters;
      for (int s = 0; s < 2; ++s)
      {
        sheets[s]     = pieces[indices[ends[s].first]].sheet;
        parameters[s] = sheets[s].rail(ends[s].second);
      }
      const auto shared = sharedSectionCurve(sheets, parameters, a.w, b.w, a.vertex, b.vertex);
      if (!shared)
        fail(op,
             no_solution,
             "a shared carrier section could not be enclosed on both sheets",
             {pieces[k].rail.edge});
      for (int s = 0; s < 2; ++s)
        charts[ends[s].first].curves.push_back(
          {shared->edge, a.ports[s], b.ports[s], {}, ends[s].second == 1});
    }
  }
  for (size_t c = 0; c < charts.size(); ++c)
  {
    const auto& curves = charts[c].curves;
    const auto& rail   = pieces[indices[c]].rail;
    if (curves.empty())
      fail(op, no_solution, "the retained boundary has no arcs", {rail.edge});
    BoundaryGraph     graph;
    std::vector<bool> used(curves.size(), false);
    for (size_t first = 0; first < curves.size(); ++first)
    {
      if (used[first])
        continue;
      size_t current = first, start = curves[first].forward ? curves[first].from : curves[first].to,
             at = start;
      do
      {
        const auto& arc = curves[current];
        if (used[current])
          fail(op, no_solution, "the end boundary repeats before closing", {rail.edge});
        used[current]        = true;
        const bool   forward = arc.from == at;
        const size_t to      = forward ? arc.to : arc.from;
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
    retainedBoundaries[indices[c]] = std::move(graph);
  }
  // A chart can meet the retained region in several connected components.
  // Nesting in its parameter domain assigns holes to their enclosing outer
  // cycle; disconnected outer cycles each own a face of the same sheet.
  for (int index : indices)
  {
    const auto graph = *retainedBoundaries[index];
    const auto loops = graph.loops();
    if (!loops || loops->size() <= 1)
      continue;
    const auto& sheet = pieces[index].sheet;
    const bool  ccw =
      sheet.swapped ? sheet.across[1] < sheet.across[0] : sheet.across[1] > sheet.across[0];
    std::vector<TopoDS_Face> domains;
    std::vector<gp_Pnt2d>    samples;
    std::vector<bool>        outer;
    for (const auto& loop : *loops)
    {
      TopoDS_Wire wire;
      builder.MakeWire(wire);
      for (const auto& edge : loop)
        builder.Add(wire, ccw ? edge : TopoDS::Edge(edge.Reversed()));
      wire.Closed(true);
      TopoDS_Face face;
      builder.MakeFace(face, sheet.surface, tolerance);
      builder.Add(face, wire);
      const auto infinite = BRepTopAdaptor_FClass2d(face, tol::parameter).PerformInfinitePoint();
      outer.push_back(infinite == TopAbs_OUT);
      if (!outer.back())
      {
        builder.MakeFace(face, sheet.surface, tolerance);
        builder.Add(face, wire.Reversed());
      }
      domains.push_back(face);
      double     first, last;
      const auto pc = BRep_Tool::CurveOnSurface(loop.front(), sheet.surface, {}, first, last);
      samples.push_back(pc->Value((first + last) / 2));
    }
    std::vector<std::vector<bool>> inside(loops->size(), std::vector<bool>(loops->size(), false));
    std::vector<int>               depth(loops->size(), 0);
    for (size_t i = 0; i < loops->size(); ++i)
      for (size_t j = 0; j < loops->size(); ++j)
        if (i != j)
        {
          inside[i][j] =
            BRepTopAdaptor_FClass2d(domains[j], tol::parameter).Perform(samples[i], false)
            == TopAbs_IN;
          depth[i] += inside[i][j];
        }
    std::vector<BoundaryGraph> components;
    for (size_t i = 0; i < loops->size(); ++i)
    {
      if (outer[i] != (depth[i] % 2 == 0))
        fail(op,
             no_solution,
             "retained cycles disagree with parameter-domain nesting",
             {pieces[index].rail.edge});
      if (!outer[i])
        continue;
      BoundaryGraph part;
      for (size_t j = 0; j < loops->size(); ++j)
        if (i == j || (!outer[j] && inside[j][i] && depth[j] == depth[i] + 1))
          for (const auto& edge : (*loops)[j])
            for (const auto& arc : graph.arcs)
              if (arc.edge.IsEqual(edge))
                part.arcs.push_back(arc);
      components.push_back(std::move(part));
    }
    for (size_t i = 0; i < components.size(); ++i)
    {
      int target = index;
      if (i)
      {
        target      = int(pieces.size());
        Piece piece = pieces[index];
        piece.swallowed.clear();
        pieces.push_back(std::move(piece));
        retainedBoundaries.emplace_back();
      }
      retainedBoundaries[target] = std::move(components[i]);
    }
  }
  n = int(pieces.size());
  contactIntervals.resize(n);
  section.resize(n);
  stopAt.resize(n);
  blends.resize(n);
  chain.resize(n);
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
  BooleanPhase phase("fillet2.endInstall");
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
      {
        int uses = 0;
        for (int other = 0; other < n; ++other)
          if (retainedBoundaries[other])
            for (const auto& shared : retainedBoundaries[other]->arcs)
              if (shared.edge.IsSame(edge))
                ++uses;
        if (uses != 2)
          fail(op,
               no_solution,
               "an internal boundary arc does not join two carrier charts",
               {edge});
        next.arcs.push_back(arc);
        continue;
      }
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
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> interiorVertices;
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
    {
      int owner = 0, incidence = 0;
      for (const auto& use : boundary)
        if (TopExp::FirstVertex(use.edge).IsSame(vertex)
            || TopExp::LastVertex(use.edge).IsSame(vertex))
        {
          if (owner && owner != use.owner)
            fail(op, no_solution, "an interior boundary event changes original face", {vertex});
          owner = use.owner;
          ++incidence;
        }
      if (!owner || incidence != 2)
        fail(op, no_solution, "an interior boundary event has unresolved incidence", {vertex});
      GeomAPI_ProjectPointOnSurf projection(BRep_Tool::Pnt(vertex),
                                            placedFace(TopoDS::Face(oldFaces(owner))).surface);
      if (!projection.NbPoints() || projection.LowerDistance() > tolerance)
        fail(op, no_solution, "an interior boundary event misses its original surface", {vertex});
      interiorVertices.Add(vertex);
    }
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
            const auto vertex = TopExp::LastVertex(parts[j], true);
            const auto inherited = pc->Value(parts[j].Orientation() == TopAbs_FORWARD ? b : a);
            // The new event is an exact surface point. The inherited edge pcurve
            // supplies its chart period, not its geometric position: that pcurve
            // may have the larger approximation tolerance of the original edge.
            GeomAPI_ProjectPointOnSurf projection(BRep_Tool::Pnt(vertex), surface);
            if (!projection.NbPoints() || projection.LowerDistance() > tolerance)
              fail(op, no_solution, "a boundary event misses its original surface", {vertex});
            double u, v;
            projection.LowerDistanceParameters(u, v);
            nodes.push_back({vertex, nearPeriod(gp_Pnt2d(u, v), inherited, surface)});
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
      {
        if (!interiorVertices.Contains(vertex))
          fail(op,
               no_solution,
               "a shared boundary vertex has no original parameter port",
               {vertex});
        nodes.push_back({vertex, uv});
      }
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
    if (loops->empty())
    {
      editOf(edits, face);
      retainedFaces.push_back({face, {}});
      continue;
    }
    // Arcs above are occurrences in the forward face. Add them in that frame:
    // TopoDS_Builder::Add otherwise reverses each wire for a reversed parent.
    auto out = TopoDS::Face(forwardFace.EmptyCopied());
    for (const auto& loop : *loops)
    {
      TopoDS_Wire wire;
      builder.MakeWire(wire);
      for (const auto& edge : loop)
        builder.Add(wire, edge);
      wire.Closed(true);
      builder.Add(out, wire);
    }
    out.Orientation(face.Orientation());
    editOf(edits, face);
    retainedFaces.push_back({face, {out}});
  }
}
} // namespace BRepExt::blending::detail
