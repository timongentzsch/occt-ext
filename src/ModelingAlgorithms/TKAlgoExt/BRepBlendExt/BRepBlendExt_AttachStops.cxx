// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Blend: The faces that stop a blend at an end, and the notches between two blends.

#include <BRepBlendExt_Attach.hxx>

namespace BRepExt::blending::detail
{

const Attach::TrimCurve& Attach::trimOf(int k, int side)
{
  const auto key = std::make_pair(k, side);
  if (const auto it = trimCurves.find(key); it != trimCurves.end())
    return it->second;
  const Piece& piece = pieces[k];
  const Sheet& sheet = piece.sheet;
  // A full round's sections end on the face beyond, touching it: the contact is the edge there.
  if (piece.rail.law == Law::Full)
    return trimCurves.emplace(key, TrimCurve{sheet.contact(side), sheet.contactUV(side), {}})
      .first->second;
  const auto surface = BRep_Tool::Surface(piece.rail.s[side].trim);
  // Along the sheet's rail as its parameter rises, which interpolation needs.
  const bool   part           = refines(piece.uAt(0) != 0 || piece.uAt(1) != 1);
  const double s0             = std::min(sheet.rail(part ? piece.uAt(0) : 0),
                                                    sheet.rail(part ? piece.uAt(1) : 1)),
               s1             = std::max(sheet.rail(part ? piece.uAt(0) : 0),
                                         sheet.rail(part ? piece.uAt(1) : 1));
  const double lo             = std::min(sheet.across[0], sheet.across[1]),
               hi             = std::max(sheet.across[0], sheet.across[1]);
  constexpr int         count = 257;
  std::vector<gp_Pnt>   points;
  std::vector<gp_Pnt2d> uvs;
  std::vector<double>   at;
  TrimCurve             out;
  int                   missed = 0;
  for (int i = 0; i < count; ++i)
  {
    const double             rail    = s0 + (s1 - s0) * i / (count - 1);
    const Handle(Geom_Curve) section = new Geom_TrimmedCurve(sheet.section(rail), lo, hi);
    // Coming from the other contact, where the section first meets the face's surface.
    GeomAPI_IntCS         meet(section, surface);
    std::optional<double> v;
    // The first from the other contact, then the one nearest the last (a face leaning over the
    // rail meets a section twice past the place it is tangent to the sheet).
    const double from = out.across.empty() ? sheet.across[1 - side] : out.across.back();
    if (meet.IsDone())
      for (int j = 1; j <= meet.NbPoints(); ++j)
      {
        double pu, pv, w;
        meet.Parameters(j, pu, pv, w);
        if (!v || std::abs(w - from) < std::abs(*v - from))
          v = w;
      }
    // Where the face reaches the sections only along one run of the rail (a face leaning over
    // the rail takes the sections away past the blend's end), the edge is that run's.
    if (!v)
    {
      if (!points.empty())
        missed = i;
      continue;
    }
    if (missed)
      fail(op,
           no_solution,
           "the face beyond a narrow face does not cross the blend's section",
           {piece.rail.edge});
    points.push_back(section->Value(*v));
    uvs.push_back(sheet.uv(rail, *v));
    at.push_back(rail);
    out.across.push_back(*v);
  }
  if (points.size() < 8)
    fail(op,
         no_solution,
         "the face beyond a narrow face does not cross the blend's section",
         {piece.rail.edge});
  const int                    made = int(points.size());
  NCollection_Array1<gp_Pnt>   pointArray(1, made);
  NCollection_Array1<gp_Pnt2d> uvArray(1, made);
  NCollection_Array1<double>   atArray(1, made);
  for (int i = 0; i < made; ++i)
  {
    pointArray(i + 1) = points[i];
    uvArray(i + 1)    = uvs[i];
    atArray(i + 1)    = at[i];
  }
  GeomAPI_Interpolate   curve(new NCollection_HArray1<gp_Pnt>(pointArray),
                            new NCollection_HArray1<double>(atArray),
                            false,
                            tol::length);
  Geom2dAPI_Interpolate onSheet(new NCollection_HArray1<gp_Pnt2d>(uvArray),
                                new NCollection_HArray1<double>(atArray),
                                false,
                                tol::parameter);
  curve.Perform();
  onSheet.Perform();
  if (!curve.IsDone() || !onSheet.IsDone())
    fail(op,
         no_solution,
         "the blend's edge on the face beyond a narrow face cannot be made",
         {piece.rail.edge});
  out.curve   = curve.Curve();
  out.onSheet = onSheet.Curve();
  return trimCurves.emplace(key, out).first->second;
}

//! Where a notch's narrow face goes on past the face that stops the blend: the surface of that face
//! is taken to where the blend's edge on the face beyond meets it, up the face beyond to the far
//! edge, and across the strip to the vertex.
void Attach::capNotch(int                  k,
                      int                  side,
                      int                  end,
                      const TopoDS_Vertex& vertex,
                      const TopoDS_Face&   stop,
                      gp_Pnt2d&            cornerOnStop)
{
  const Piece& piece   = pieces[k];
  const Rail&  rail    = piece.rail;
  const Side&  narrow  = rail.s[side];
  const auto   placed  = BRep_Tool::Surface(stop);
  const gp_Pnt V       = BRep_Tool::Pnt(vertex);
  const auto   nearest = [&](const Handle(Geom_Curve)& curve, const char* why) {
    GeomAPI_IntCS                            meet(curve, placed);
    std::optional<std::pair<gp_Pnt, double>> best;
    if (meet.IsDone())
      for (int i = 1; i <= meet.NbPoints(); ++i)
      {
        double pu, pv, w;
        meet.Parameters(i, pu, pv, w);
        if (!best || meet.Point(i).Distance(V) < best->first.Distance(V))
          best = {meet.Point(i), w};
      }
    if (!best)
      fail(op, unsupported_geometry, why, {rail.edge});
    return *best;
  };
  const auto [Q, sQ] = nearest(
    trimOf(k, side).curve,
    "the blend's edge on the face beyond a narrow face does not reach the face that stops it");
  double     ef, el;
  const auto far =
    nearest(BRep_Tool::Curve(narrow.trimEdge, ef, el),
            "the face that stops a blend at a notch does not cross the narrow face's far edge");
  builder.MakeVertex(corner[k][side][end], Q, tolerance);
  u[k][side][end] = (sQ - piece.sheet.alpha) / piece.sheet.beta;
  TopoDS_Vertex F;
  builder.MakeVertex(F, far.first, tolerance);
  trimFar[{k, side, end}] = F;
  // The cap's two edges, each the segment on both faces' surfaces.
  const auto segment =
    [&](const TopoDS_Vertex& from, const TopoDS_Vertex& to, const TopoDS_Face& face) {
      const gp_Pnt p0 = BRep_Tool::Pnt(from), p1 = BRep_Tool::Pnt(to);
      const double length = p0.Distance(p1);
      if (!(length > limit::degenerate * rail.scale))
        fail(op, no_solution, "the cap of a blend at a notch collapses", {rail.edge});
      const Handle(Geom_Curve) line(GC_MakeSegment(p0, p1).Value());
      const auto               other = BRep_Tool::Surface(face);
      for (const auto& surface : {placed, other})
      {
        GeomAPI_ProjectPointOnSurf mid(line->Value(length / 2), surface);
        if (mid.NbPoints() == 0 || mid.LowerDistance() > rail.resolved(limit::together))
          fail(op,
               unsupported_geometry,
               "the face that stops a blend at a notch does not cross the faces beside it in lines",
               {rail.edge});
      }
      const auto onStop  = GeomProjLib::Curve2d(line, 0, length, placed),
                 onOther = GeomProjLib::Curve2d(line, 0, length, other);
      if (onStop.IsNull() || onOther.IsNull())
        fail(op,
             no_solution,
             "the cap of a blend at a notch cannot be laid on its faces",
             {rail.edge});
      const auto [stopBasis, stopPlace]   = surfaceOf(stop);
      const auto [otherBasis, otherPlace] = surfaceOf(face);
      const TopoDS_Edge edge =
        edgeBetween(line,
                    from,
                    0,
                    to,
                    length,
                    {{onStop, stopBasis, stopPlace}, {onOther, otherBasis, otherPlace}},
                    tolerance);
      editOf(edits, stop).bridges.push_back(edge);
      editOf(edits, face).bridges.push_back(edge);
      continued.push_back({edge, {stop, face}});
      return onStop;
    };
  cornerOnStop = segment(corner[k][side][end], F, narrow.trim)->Value(0);
  segment(F, vertex, narrow.face);
}

//! A stop face a notch's blend uses up: a triangle of the vertex, the notch's narrow face's far
//! vertex and a third, whose edge from the vertex on the other support and edge from the far
//! vertex on the face beyond the narrow face end at that third vertex, where the support and the
//! face beyond meet in an edge. The other support's contact crosses that edge, and there the blend
//! closes: the face beyond, which a notch's blend is cut by, takes the blend's end whole (the
//! round meets it along the notch's edge, which runs from the blend's start to that crossing,
//! turning back along the rail where the face leans over it, with no section of the blend at that
//! end), and the stop face is gone with the edges it had (a back face leaning over a front edge).
bool Attach::consumed(int                  k,
                      int                  end,
                      const TopoDS_Vertex& vertex,
                      const TopoDS_Face&   stop,
                      const TopoDS_Edge (&beside)[2])
{
  const Piece& piece = pieces[k];
  const Rail&  rail  = piece.rail;
  const Sheet& sheet = piece.sheet;
  const int    notch = rail.s[0].trim.IsNull() ? 1 : 0, plain = 1 - notch;
  const Side&  narrow = rail.s[notch];
  if (narrow.trim.IsNull() || rail.law == Law::Full || !rail.s[plain].trim.IsNull()
      || rail.s[plain].cliffCurve || narrow.cliffCurve)
    return false;
  const auto otherEnd = [&](const TopoDS_Edge& edge, const TopoDS_Vertex& from) {
    TopoDS_Vertex a, z;
    TopExp::Vertices(edge, a, z);
    return a.IsSame(from) ? z : a;
  };
  const auto edgeAt = [&](const TopoDS_Vertex& at,
                          const TopoDS_Face&   one,
                          const TopoDS_Face&   two,
                          const TopoDS_Edge&   not1,
                          const TopoDS_Edge&   not2) {
    TopoDS_Edge found;
    for (const auto& shape : edgesOfVertex.FindFromKey(at))
      if (!shape.IsSame(not1) && !shape.IsSame(not2) && onFace(shape, one) && onFace(shape, two))
        found = TopoDS::Edge(shape);
    return found;
  };
  const TopoDS_Vertex far = otherEnd(beside[notch], vertex);
  TopoDS_Vertex       t0, t1;
  TopExp::Vertices(narrow.trimEdge, t0, t1);
  if (!far.IsSame(t0) && !far.IsSame(t1))
    return false;
  const TopoDS_Edge down = edgeAt(far, narrow.trim, stop, narrow.trimEdge, beside[notch]);
  if (down.IsNull())
    return false;
  const TopoDS_Vertex third = otherEnd(down, far);
  if (!third.IsSame(otherEnd(beside[plain], vertex)))
    return false;
  const TopoDS_Edge meet = edgeAt(third, rail.s[plain].face, narrow.trim, down, beside[plain]);
  if (meet.IsNull())
    return false;
  const double middle  = piece.uAt(0) / 2 + piece.uAt(1) / 2;
  const auto   onPlain = crossing(piece.fitted.contact[plain],
                                end == 0 ? 0 : middle,
                                end == 0 ? middle : 1,
                                piece.uAt(end),
                                meet,
                                rail.s[plain].face);
  if (!onPlain)
    return false;
  // The blend closes to that point on the edge between the faces: its corners are one.
  u[k][0][end] = u[k][1][end] = onPlain->first;
  builder.MakeVertex(corner[k][0][end],
                     sheet.contact(plain)->Value(sheet.rail(onPlain->first)),
                     tolerance);
  corner[k][1][end] = corner[k][0][end];
  cut(meet, third, corner[k][0][end], onPlain->second);
  for (const TopoDS_Edge& gone : {beside[0], beside[1], down})
  {
    for (const auto& face : facesOfEdge.FindFromKey(gone))
      editOf(edits, face).replace.Bind(gone, TopoDS_Edge());
    swept.push_back(gone);
  }
  swept.push_back(far);
  swept.push_back(third);
  removed.push_back(stop);
  stopAt[k][end]   = stop;
  rolled[{k, end}] = stop;
  return true;
}

//! A notch's edge on the face beyond, for a blend that closes on that face: the sheet's meeting
//! with the face from the notch's first corner to the second, by the arc that lies across the
//! sheet between its contacts. The face must hold all of it on the stop face's inner side.
TopoDS_Edge Attach::rolledEdge(int k, int side, const TopoDS_Face& stop)
{
  const Piece&       piece       = pieces[k];
  const Rail&        rail        = piece.rail;
  const Sheet&       sheet       = piece.sheet;
  const TopoDS_Face& beyond      = rail.s[side].trim;
  const auto [surface, location] = surfaceOf(beyond);
  const SurfaceMeeting meeting(sheet.surface, BRep_Tool::Surface(beyond), tolerance);
  const gp_Pnt a = BRep_Tool::Pnt(corner[k][side][0]), b = BRep_Tool::Pnt(corner[k][side][1]);
  const double lo = std::min(sheet.across[0], sheet.across[1]),
               hi = std::max(sheet.across[0], sheet.across[1]);
  for (int i = 1; i <= meeting.count(); ++i)
  {
    const Handle(Geom_Curve) curve = meeting.curve(i);
    const auto               atA = onCurve(a, curve), atB = onCurve(b, curve);
    if (!atA || !atB || atA->distance > rail.resolved(limit::together)
        || atB->distance > rail.resolved(limit::together))
      continue;
    // The ways from one corner to the other, the other way round a closed curve too.
    const double                           period = curve->IsPeriodic() ? curve->Period() : 0;
    std::vector<std::pair<double, double>> ways   = {{atA->parameter, atB->parameter}};
    if (period > 0)
      ways = {{atA->parameter, atB->parameter + (atB->parameter > atA->parameter ? 0 : period)},
              {atB->parameter, atA->parameter + (atA->parameter > atB->parameter ? 0 : period)}};
    for (const auto& [t0, t1] : ways)
    {
      // Across the sheet between its contacts all the way, and on the stop face's inner side.
      bool within = true;
      for (int j = 1; j < 16 && within; ++j)
      {
        const gp_Pnt               p = curve->Value(t0 + (t1 - t0) * j / 16);
        GeomAPI_ProjectPointOnSurf on(p, sheet.surface);
        if (on.NbPoints() == 0)
        {
          within = false;
          break;
        }
        double pu, pv;
        on.LowerDistanceParameters(pu, pv);
        const double across = sheet.swapped ? pu : pv;
        within              = across > lo + tol::parameter && across < hi - tol::parameter;
      }
      if (!within)
        continue;
      const auto stopSurface = BRep_Tool::Surface(stop);
      const auto plane       = Handle(Geom_Plane)::DownCast(stopSurface);
      if (plane.IsNull())
        fail(op,
             unsupported_geometry,
             "a blend that closes on the face beyond its narrow face ends on a stop face that is "
             "not flat",
             {rail.edge});
      gp_Vec outward(plane->Axis().Direction());
      if (stop.Orientation() == TopAbs_REVERSED)
        outward.Reverse();
      for (int j = 0; j <= 16; ++j)
        if (gp_Vec(plane->Location(), curve->Value(t0 + (t1 - t0) * j / 16)).Dot(outward)
            > rail.resolved(limit::together))
          fail(op,
               unsupported_geometry,
               "a blend that closes on the face beyond its narrow face runs past the face that "
               "stops it, which the blend does not do yet",
               {rail.edge});
      const auto onSheet  = meeting.on(i, 1, std::min(t0, t1), std::max(t0, t1));
      const auto onBeyond = meeting.on(i, 2, std::min(t0, t1), std::max(t0, t1));
      if (onSheet.IsNull() || onBeyond.IsNull())
        fail(op,
             no_solution,
             "the blend's edge on the face beyond a narrow face cannot be laid on it",
             {rail.edge});
      const gp_Pnt2d cornerOnSheet = sheet.uv(sheet.rail(u[k][side][0]), sheet.across[side]);
      return edgeBetween(curve,
                         corner[k][side][0],
                         t0,
                         corner[k][side][1],
                         t1,
                         {{aligned(onSheet, t0, cornerOnSheet, sheet.surface), sheet.surface, {}},
                          {onBeyond, surface, location}},
                         tolerance);
    }
  }
  fail(op,
       no_solution,
       "the blend's edge on the face beyond a narrow face does not run between its corners",
       {rail.edge});
}

//! A face across both sides stops a piece at one end (its vertex there): the edges it shares with
//! the supports are cut back to where the contacts cross them, and the blend ends in its section
//! on that face, or its meeting with it. At a vertex of more edges, or where a contact sweeps over
//! short edges of its support, several faces stop it (stops).
void Attach::stopFace(int k, int end, const TopoDS_Vertex& vertex, OneStop& one)
{
  const Piece&         piece  = pieces[k];
  const Rail&          rail   = piece.rail;
  const auto&          edges  = edgesOfVertex.FindFromKey(vertex);
  TopoDS_Face&         stop   = one.stop;
  TopoDS_Edge(&beside)[2]     = one.beside;
  TopoDS_Vertex(&from)[2]     = one.from;
  TopoDS_Edge sliverEdge[2];
  TopoDS_Face sliverOf[2];
  bool continuedSupport[2] = {false, false};
  for (const auto& shape : edges)
  {
    if (ownEdge(piece, shape))
      continue;
    for (const auto& face : facesOfEdge.FindFromKey(shape))
    {
      int side = -1;
      for (int s = 0; s < 2; ++s)
        if (face.IsSame(rail.s[s].face))
          side = s;
      if (side >= 0)
      {
        beside[side] = TopoDS::Edge(shape);
        continue;
      }
      // A sliver of a support under the blend, beside an edge it runs on over (Blend::absorbEnds):
      // it goes, and the edge from its far corner on the support's own face stands in for its edge.
      bool sliver = false;
      for (const auto& over : piece.swallowed)
        if (onFace(over, face))
          for (int s = 0; s < 2 && !sliver; ++s)
            if (onFace(over, rail.s[1 - s].face) && !onFace(over, rail.s[s].face))
            {
              sliverOf[s]   = TopoDS::Face(face);
              sliverEdge[s] = TopoDS::Edge(shape);
              sliver        = true;
            }
      if (!sliver && rail.parameterSupports)
        for (int side = 0; side < 2 && !sliver; ++side)
          if (piece.terminalContact[end][side]
              && face.IsSame((*rail.parameterSupports)[side].face)
              && !face.IsSame(rail.s[side].face))
          {
            sliverOf[side] = TopoDS::Face(face);
            sliverEdge[side] = TopoDS::Edge(shape);
            continuedSupport[side] = true;
            sliver = true;
          }
      if (sliver)
        continue;
      if (!stop.IsNull() && !stop.IsSame(face))
        fail(op,
             unsupported_geometry,
             "an edge that ends among several faces is blended with the corner, which the blend "
             "does not do yet",
             {rail.edge});
      stop = TopoDS::Face(face);
    }
  }
  from[0] = from[1] = vertex;
  for (int s = 0; s < 2; ++s)
    if (!sliverOf[s].IsNull() && !stop.IsNull())
    {
      TopoDS_Vertex a, z;
      TopExp::Vertices(sliverEdge[s], a, z);
      const TopoDS_Vertex far = a.IsSame(vertex) ? z : a;
      for (const auto& shape : edgesOfVertex.FindFromKey(far))
        if (onFace(shape, rail.s[s].face) && onFace(shape, stop))
        {
          beside[s] = TopoDS::Edge(shape);
          from[s]   = far;
        }
      if (from[s].IsSame(vertex))
        fail(op,
             unsupported_geometry,
             "an edge that ends among several faces is blended with the corner, which the blend "
             "does not do yet",
             {rail.edge});
      if (continuedSupport[s])
      {
        const auto& terminal = *piece.terminalContact[end][s];
        if (!terminal.support.IsSame(rail.s[s].face) || !terminal.arc.IsSame(beside[s]))
          fail(op, unsupported_geometry, "terminal contact lost its restriction arc", {rail.edge});
        for (const auto& face : facesOfEdge.FindFromKey(sliverEdge[s]))
          editOf(edits, face).replace.Bind(sliverEdge[s], TopoDS_Edge());
        swept.push_back(sliverEdge[s]);
      }
      else
      {
        for (TopExp_Explorer e(sliverOf[s], TopAbs_EDGE); e.More(); e.Next())
        {
          if (facesOfEdge.Contains(e.Current()))
            for (const auto& face : facesOfEdge.FindFromKey(e.Current()))
              editOf(edits, face).replace.Bind(e.Current(), TopoDS_Edge());
          swept.push_back(e.Current());
        }
        removed.push_back(sliverOf[s]);
      }
      swept.push_back(far);
    }
  if (stop.IsNull() || beside[0].IsNull() || beside[1].IsNull())
    fail(op,
         unsupported_geometry,
         "the faces at an end of the edge are not one face across both sides",
         {rail.edge});
}

// A contact that leaves its support for the stop face before the end's half of the piece (a
// support that narrows towards the stop face, as a triangular top to its oblique end) crosses
// the edge between them in the other half.
std::optional<Meeting> Attach::earlyCrossing(int k, int end, int side, const OneStop& s) const
{
  const Piece&         piece  = pieces[k];
  const Rail&          rail   = piece.rail;
  const Sheet&         sheet  = piece.sheet;
  const double         middle = s.middle;
  const TopoDS_Edge(&beside)[2] = s.beside;
  auto found = crossing(piece.fitted.contact[side],
                        end == 0 ? middle : 0,
                        end == 0 ? 1 : middle,
                        piece.uAt(end),
                        beside[side],
                        rail.s[side].face);
  // Or further on along the parts of a split, one sheet: a seam crossed by the other contact
  // before this one reaches the stop face.
  for (int m = k; !found && piece.joint[1 - end].meets == Meets::Split
                  && pieces[m].joint[1 - end].meets == Meets::Split;)
  {
    m = pieces[m].joint[1 - end].piece;
    if (m == k || !pieces[m].rail.s[side].face.IsSame(rail.s[side].face))
      break;
    found = crossing(pieces[m].fitted.contact[side],
                     0,
                     1,
                     pieces[m].uAt(end),
                     beside[side],
                     rail.s[side].face);
    if (found)
      found->first = (pieces[m].sheet.rail(found->first) - sheet.alpha) / sheet.beta;
  }
  return found;
}

bool Attach::simpleStop(int k, int end, const OneStop& s) const
{
  const Piece&         piece  = pieces[k];
  const Rail&          rail   = piece.rail;
  const double         middle = s.middle;
  const int            cliff  = s.cliff;
  const TopoDS_Face&   stop   = s.stop;
  const TopoDS_Edge(&beside)[2] = s.beside;
  const TopoDS_Vertex(&from)[2] = s.from;
  const auto early = [&](int side) { return earlyCrossing(k, end, side, s); };
  bool simple = true;
  for (int side = 0; side < 2 && simple; ++side)
    simple = side == cliff || !rail.s[side].trim.IsNull()
             || crossing(piece.fitted.contact[side],
                         end == 0 ? 0 : middle,
                         end == 0 ? middle : 1,
                         piece.uAt(end),
                         beside[side],
                         rail.s[side].face)
             || beyond(piece.fitted.contact[side],
                       end == 0 ? 0 : middle,
                       end == 0 ? middle : 1,
                       piece.uAt(end),
                       beside[side],
                       from[side],
                       rail.s[side].face,
                       stop,
                       tolerance)
             || continuation(piece.fitted.contact[side],
                             end == 0 ? 0 : middle,
                             end == 0 ? middle : 1,
                             piece.uAt(end),
                             beside[side],
                             from[side],
                             rail.s[side].face,
                             stop,
                             tolerance)
             || early(side);
  return simple;
}

void Attach::notchStopCorner(int                  k,
                             int                  end,
                             int                  side,
                             const TopoDS_Vertex& vertex,
                             const OneStop&       s,
                             gp_Pnt2d&            cornerOnStop)
{
  const Piece&         piece  = pieces[k];
  const Rail&          rail   = piece.rail;
  const Sheet&         sheet  = piece.sheet;
  const TopoDS_Face&   stop   = s.stop;
  const TopoDS_Edge(&beside)[2] = s.beside;
  // A notch: the narrow face's edge on the stop face goes with the face, and the face
  // beyond's edge on the stop face from the far edge's vertex is cut back to where the
  // blend's edge on the face beyond meets the stop face: the corner.
  const TopoDS_Face& beyond = rail.s[side].trim;
  TopoDS_Vertex      a, z, f0, f1;
  TopExp::Vertices(beside[side], a, z);
  TopExp::Vertices(rail.s[side].trimEdge, f0, f1);
  const TopoDS_Vertex across = a.IsSame(vertex) ? z : a;
  const auto&         trim   = trimOf(k, side);
  if (!across.IsSame(f0) && !across.IsSame(f1))
  {
    capNotch(k, side, end, vertex, stop, cornerOnStop);
    return;
  }
  TopoDS_Edge down;
  for (const auto& shape : edgesOfVertex.FindFromKey(across))
    if (onFace(shape, beyond) && onFace(shape, stop))
      down = TopoDS::Edge(shape);
  if (down.IsNull())
    fail(op,
         unsupported_geometry,
         "the face beyond a narrow face does not meet the face that stops the blend",
         {rail.edge});
  GeomAPI_IntCS                            meet(trim.curve, BRep_Tool::Surface(stop));
  const gp_Pnt                             V = BRep_Tool::Pnt(vertex);
  std::optional<std::pair<gp_Pnt, double>> best;
  if (meet.IsDone())
    for (int i = 1; i <= meet.NbPoints(); ++i)
    {
      double pu, pv, w;
      meet.Parameters(i, pu, pv, w);
      if (!best || meet.Point(i).Distance(V) < best->first.Distance(V))
        best = {meet.Point(i), w};
    }
  if (!best)
    fail(
      op,
      no_solution,
      "the blend's edge on the face beyond a narrow face does not reach the face that stops it",
      {rail.edge});
  double                      f, l;
  const auto                  curve = BRep_Tool::Curve(down, f, l);
  GeomAPI_ProjectPointOnCurve on(best->first, curve, f, l);
  if (on.NbPoints() == 0 || on.LowerDistance() > rail.resolved(limit::together))
    fail(op,
         no_solution,
         "the blend's edge on the face beyond a narrow face misses that face's edge on the "
         "stop face",
         {rail.edge});
  builder.MakeVertex(corner[k][side][end], best->first, tolerance);
  u[k][side][end] = (best->second - sheet.alpha) / sheet.beta;
  cut(down, across, corner[k][side][end], on.LowerDistanceParameter());
  editOf(edits, stop).replace.Bind(beside[side], TopoDS_Edge());
  if (!pieces[k].rail.s[side].cliffWhole)
    editOf(edits, rail.s[side].face).replace.Bind(beside[side], TopoDS_Edge());
  swept.push_back(beside[side]);
  double g0, l0;
  cornerOnStop =
    BRep_Tool::CurveOnSurface(down, stop, g0, l0)->Value(on.LowerDistanceParameter());
}

void Attach::cliffStopCorner(int                  k,
                             int                  end,
                             int                  side,
                             const TopoDS_Vertex& vertex,
                             const OneStop&       s,
                             gp_Pnt2d&            cornerOnStop)
{
  const Piece&         piece  = pieces[k];
  const Rail&          rail   = piece.rail;
  const Sheet&         sheet  = piece.sheet;
  const TopoDS_Face&   stop   = s.stop;
  const TopoDS_Edge(&beside)[2] = s.beside;
  const Side&   narrow = rail.s[side];
  TopoDS_Vertex a, z, f0, f1;
  TopExp::Vertices(beside[side], a, z);
  TopExp::Vertices(narrow.cliff, f0, f1);
  const TopoDS_Vertex across = a.IsSame(vertex) ? z : a;
  if (across.IsSame(f0) || across.IsSame(f1))
  {
    // On a cliff the corner is the far edge's own vertex: the narrow face's edge on the
    // stop face goes, and the blend's end runs to that vertex.
    corner[k][side][end] = across;
    GeomAPI_ProjectPointOnCurve on(BRep_Tool::Pnt(corner[k][side][end]), sheet.contact(side));
    if (on.NbPoints() == 0 || on.LowerDistance() > rail.resolved(limit::together))
      fail(op,
           no_solution,
           "the blend's contact misses the end of the face's far edge it rests on",
           {rail.edge});
    u[k][side][end] = (on.LowerDistanceParameter() - sheet.alpha) / sheet.beta;
    editOf(edits, stop).replace.Bind(beside[side], TopoDS_Edge());
    if (!narrow.cliffWhole)
      editOf(edits, narrow.face).replace.Bind(beside[side], TopoDS_Edge());
    swept.push_back(beside[side]);
    double f0, l0;
    cornerOnStop = BRep_Tool::CurveOnSurface(beside[side], stop, f0, l0)
                             ->Value(BRep_Tool::Parameter(corner[k][side][end], beside[side]));
    return;
  }
  // The narrow face goes on past the stop: the stop face's surface crosses the strip from
  // the vertex to the far edge, and is carried over it to cap the blend's end there (the

  if (narrow.cliffWhole)
    fail(op, no_solution, "the blend's end on a cliff misses the far edge's end", {rail.edge});
  const auto            placed = BRep_Tool::Surface(stop);
  GeomAPI_IntCS         meet(narrow.cliffCurve, placed);
  const gp_Pnt          V = BRep_Tool::Pnt(vertex);
  std::optional<gp_Pnt> F;
  if (meet.IsDone())
    for (int i = 1; i <= meet.NbPoints(); ++i)
      if (!F || meet.Point(i).Distance(V) < F->Distance(V))
        F = meet.Point(i);
  if (!F || !(F->Distance(V) > limit::degenerate * rail.scale))
    fail(op,
         unsupported_geometry,
         "the face that stops a blend on a cliff does not cross the cliff's edge",
         {rail.edge});
  builder.MakeVertex(corner[k][side][end], *F, tolerance);
  GeomAPI_ProjectPointOnCurve on(*F, sheet.contact(side));
  if (on.NbPoints() == 0 || on.LowerDistance() > rail.resolved(limit::together))
    fail(op,
         no_solution,
         "the blend's contact misses where the face that stops it crosses the cliff's edge",
         {rail.edge});
  u[k][side][end] = (on.LowerDistanceParameter() - sheet.alpha) / sheet.beta;
  // The cap's edge: the segment from there to the vertex, on both faces' surfaces.
  const Handle(Geom_Curve) line(GC_MakeSegment(*F, V).Value());
  const double             length       = F->Distance(V);
  const auto               narrowPlaced = BRep_Tool::Surface(narrow.face);
  for (const auto& surface : {placed, narrowPlaced})
  {
    GeomAPI_ProjectPointOnSurf mid(line->Value(length / 2), surface);
    if (mid.NbPoints() == 0 || mid.LowerDistance() > rail.resolved(limit::together))
      fail(op,
           unsupported_geometry,
           "the face that stops a blend on a cliff does not cross the narrow face in a line",
           {rail.edge});
  }
  const auto onStop   = GeomProjLib::Curve2d(line, 0, length, placed),
             onNarrow = GeomProjLib::Curve2d(line, 0, length, narrowPlaced);
  if (onStop.IsNull() || onNarrow.IsNull())
    fail(op,
         no_solution,
         "the cap of a blend on a cliff cannot be laid on its faces",
         {rail.edge});
  const auto [stopBasis, stopPlace]     = surfaceOf(stop);
  const auto [narrowBasis, narrowPlace] = surfaceOf(narrow.face);
  const TopoDS_Edge cap =
    edgeBetween(line,
                corner[k][side][end],
                0,
                vertex,
                length,
                {{onStop, stopBasis, stopPlace}, {onNarrow, narrowBasis, narrowPlace}},
                tolerance);
  editOf(edits, stop).bridges.push_back(cap);
  editOf(edits, narrow.face).bridges.push_back(cap);
  continued.push_back({cap, {stop, narrow.face}});
  cornerOnStop = onStop->Value(0);
}

void Attach::continuedStopCorner(int                  k,
                                 int                  end,
                                 int                  side,
                                 const TopoDS_Vertex& vertex,
                                 const OneStop&       s,
                                 Continuation&        more,
                                 gp_Pnt2d&            cornerOnStop)
{
  const Rail&          rail   = pieces[k].rail;
  const Sheet&         sheet  = pieces[k].sheet;
  const TopoDS_Face&   stop   = s.stop;
  // On the stop face's surface itself where the meeting it was found on, true to the
  // intersector's tolerance, leaves it further off than a corner may be.
  double     place         = more.u;
  const auto onStopSurface = [&](double at) {
    GeomAPI_ProjectPointOnSurf project(sheet.contact(side)->Value(sheet.rail(at)),
                                       BRep_Tool::Surface(stop));
    return project.NbPoints() ? project.LowerDistance() : HUGE_VAL;
  };
  if (onStopSurface(place) > tolerance)
  {
    const double at    = sheet.rail(more.u),
                 reach = 1e-3 * std::abs(sheet.rail(1) - sheet.rail(0));
    GeomAPI_IntCS onStop(new Geom_TrimmedCurve(sheet.contact(side), at - reach, at + reach),
                         BRep_Tool::Surface(stop));
    for (int i = 1; onStop.IsDone() && i <= onStop.NbPoints(); ++i)
    {
      double su, sv, w;
      onStop.Parameters(i, su, sv, w);
      const double q = (w - sheet.alpha) / sheet.beta;
      if (i == 1 || std::abs(q - more.u) < std::abs(place - more.u))
        place = q;
    }
  }
  u[k][side][end]                   = place;
  const gp_Pnt                there = sheet.contact(side)->Value(sheet.rail(place));
  GeomAPI_ProjectPointOnCurve onMeeting(there, more.curve);
  if (place != more.u && onMeeting.NbPoints() > 0)
    more.from = onMeeting.LowerDistanceParameter();
  // A hair over the gap: the check reads the vertex against its tolerance to rounding.
  builder.MakeVertex(
    corner[k][side][end],
    there,
    std::max(tolerance, 1.01 * there.Distance(more.curve->Value(more.from))));
  // The new piece of the edge between the support and the stop face, from the
  // corner back to the vertex, in both faces' wires.
  const auto [supportSurface, supportPlace] = surfaceOf(rail.s[side].face);
  const auto [stopSurface, stopPlace]       = surfaceOf(stop);
  const TopoDS_Edge piece                   = edgeBetween(
    more.curve,
    corner[k][side][end],
    more.from,
    vertex,
    more.to,
    {{more.onSupport, supportSurface, supportPlace}, {more.onStop, stopSurface, stopPlace}},
    tolerance);
  editOf(edits, rail.s[side].face).bridges.push_back(piece);
  editOf(edits, stop).bridges.push_back(piece);
  continued.push_back({piece, {rail.s[side].face, stop}});
  cornerOnStop = more.onStop->Value(more.from);
}

void Attach::stopCorner(int                  k,
                        int                  end,
                        int                  side,
                        const TopoDS_Vertex& vertex,
                        const OneStop&       s,
                        gp_Pnt2d&            cornerOnStop)
{
  const Piece&         piece  = pieces[k];
  const Rail&          rail   = piece.rail;
  const Sheet&         sheet  = piece.sheet;
  const double         middle = s.middle;
  const TopoDS_Face&   stop   = s.stop;
  const TopoDS_Edge(&beside)[2] = s.beside;
  const TopoDS_Vertex(&from)[2] = s.from;
  if (!rail.s[side].trim.IsNull())
  {
    // A notch: the narrow face's edge on the stop face goes with the face, and the face
    // beyond's edge on the stop face from the far edge's vertex is cut back to where the
    notchStopCorner(k, end, side, vertex, s, cornerOnStop);
    return;
  }
  if (side == s.cliff)
  {
    cliffStopCorner(k, end, side, vertex, s, cornerOnStop);
    return;
  }
  auto found = crossing(piece.fitted.contact[side],
                        end == 0 ? 0 : middle,
                        end == 0 ? middle : 1,
                        piece.uAt(end),
                        beside[side],
                        rail.s[side].face);
  if (!found)
    found = beyond(piece.fitted.contact[side],
                   end == 0 ? 0 : middle,
                   end == 0 ? middle : 1,
                   piece.uAt(end),
                   beside[side],
                   from[side],
                   rail.s[side].face,
                   stop,
                   tolerance);
  std::optional<Continuation> more;
  if (!found)
    more = continuation(piece.fitted.contact[side],
                        end == 0 ? 0 : middle,
                        end == 0 ? middle : 1,
                        piece.uAt(end),
                        beside[side],
                        from[side],
                        rail.s[side].face,
                        stop,
                        tolerance);
  if (!found && !more)
    found = earlyCrossing(k, end, side, s);
  if (!found)
  {
    if (!more)
      fail(op,
           unsupported_geometry,
           "the blend runs on past the face that stops it, where its supports do not meet the "
           "stop face's surface",
           {rail.edge});
    continuedStopCorner(k, end, side, vertex, s, *more, cornerOnStop);
    return;
  }
  u[k][side][end] = found->first;
  builder.MakeVertex(corner[k][side][end],
                     sheet.contact(side)->Value(sheet.rail(found->first)),
                     tolerance);
  cut(beside[side], from[side], corner[k][side][end], found->second);
  {
    double f0, l0;
    cornerOnStop =
      BRep_Tool::CurveOnSurface(beside[side], stop, f0, l0)->Value(found->second);
  }
}

void Attach::stopSection(int k, int end, const TopoDS_Face& stop, const gp_Pnt2d& cornerOnStop)
{
  const Piece&         piece  = pieces[k];
  const Rail&          rail   = piece.rail;
  const Sheet&         sheet  = piece.sheet;
  const auto [surface, location] = surfaceOf(stop);
  const auto     placed          = BRep_Tool::Surface(stop);
  const gp_Pnt2d cornerOnSheet   = sheet.uv(sheet.rail(u[k][0][end]), sheet.across[0]);
  // The sheet's own section where it lies in the face, else the two surfaces' meeting.
  const double along = (u[k][0][end] + u[k][1][end]) / 2;
  const auto   curve = sheet.section(sheet.rail(along));
  bool         iso   = true;
  for (int i = 0; i <= 8 && iso; ++i)
  {
    GeomAPI_ProjectPointOnSurf on(
      curve->Value(sheet.across[0] + (sheet.across[1] - sheet.across[0]) * i / 8),
      placed);
    iso = on.NbPoints() > 0 && on.LowerDistance() < tolerance;
  }
  TopoDS_Edge& edge = section[k][end];
  // A notch's side ends where the face beyond meets the section, short of the contact.
  double ends[2] = {sheet.across[0], sheet.across[1]};
  for (int side = 0; side < 2; ++side)
    if (!rail.s[side].trim.IsNull())
      ends[side] = trimAcross(k, side, sheet.rail(along));
  if (iso)
  {
    edge = edgeBetween(sheet.section(sheet.rail(along)),
                       corner[k][0][end],
                       ends[0],
                       corner[k][1][end],
                       ends[1],
                       {},
                       tol::length);
    builder.UpdateEdge(edge,
                       sheet.sectionUV(sheet.rail(along)),
                       sheet.surface,
                       TopLoc_Location(),
                       tolerance);
    double f, l;
    BRep_Tool::Range(edge, f, l);
    const auto plane  = Handle(Geom_Plane)::DownCast(placed);
    const auto onStop = plane.IsNull() ? GeomProjLib::Curve2d(curve, f, l, placed)
                                       : GeomAPI::To2d(curve, plane->Pln());
    if (onStop.IsNull())
      fail(op,
           no_solution,
           "the blend's end cannot be laid on the face that stops it",
           {rail.edge});
    builder.UpdateEdge(edge,
                       aligned(onStop, ends[0], cornerOnStop, placed),
                       surface,
                       location,
                       tolerance);
    builder.Range(edge, f, l);
  }
  else
  {
    if ((!rail.s[0].trim.IsNull() || !rail.s[1].trim.IsNull()) && rail.law != Law::Full)
      fail(op,
           unsupported_geometry,
           "a notch's blend ends on a face it does not meet in a section, which the blend does not "
           "make yet",
           {rail.edge});
    auto met = meetingEdge(placedSheet(sheet),
                           cornerOnSheet,
                           placedFace(stop),
                           cornerOnStop,
                           corner[k][0][end],
                           corner[k][1][end],
                           rail.scale);
    if (refines(!met))
      met = planeMeeting(sheet,
                         stop,
                         cornerOnSheet,
                         sheet.uv(sheet.rail(u[k][1][end]), sheet.across[1]),
                         corner[k][0][end],
                         corner[k][1][end]);
    if (!met)
      fail(op,
           no_solution,
           "the blend's end on the face that stops it does not reach both contacts",
           {rail.edge});
    edge = *met;
  }
  editOf(edits, stop).bridges.push_back(edge);
}

void Attach::stop(int k, int end, const TopoDS_Vertex& vertex)
{
  const Piece& piece = pieces[k];
  const Rail&  rail  = piece.rail;
  const auto&  edges = edgesOfVertex.FindFromKey(vertex);
  if (edges.Extent() > 3)
  {
    stops(k, end, vertex);
    return;
  }
  if (edges.Extent() != 3)
    fail(op,
         unsupported_geometry,
         "where blended edges meet other edges, the corner is blended too, which the blend does "
         "not do yet",
         {rail.edge});
  OneStop s;
  stopFace(k, end, vertex, s);
  s.middle = piece.uAt(0) / 2 + piece.uAt(1) / 2;
  s.cliff  = rail.s[0].cliffCurve ? 0 : rail.s[1].cliffCurve ? 1 : -1;
  if (!simpleStop(k, end, s))
  {
    if (!consumed(k, end, vertex, s.stop, s.beside))
      stops(k, end, vertex);
    return;
  }
  stopAt[k][end] = s.stop;
  gp_Pnt2d cornerOnStop2d[2];
  for (int side = 0; side < 2; ++side)
    stopCorner(k, end, side, vertex, s, cornerOnStop2d[side]);
  stopSection(k, end, s.stop, cornerOnStop2d[0]);
}

//! Two blends of one kind whose edges meet at a vertex of four faces with none in common (a notch).
//! The one of the edges at the vertex that lies between a face of each blend is in question; on
//! it, the blend whose contact reaches
//! the edge before it reaches the other blend's sheet wins: the edge is cut back to where that
//! contact crosses it, the winner's end runs from there along its sheet on the loser's face (which
//! grows to it where needed) to where the loser's contact meets the winner's sheet, and the
//! loser's contact ends there. Between those two places both ends run where the sheets meet.
//! A contact that never crosses the edge, which goes on straight from its own blended edge on
//! its face, wins it whole: the edge goes, and its blend runs on along it, trimmed by the loser's
//! notch (whose wall stands parallel to the other blended edge).
//! Everything is found before anything is changed: false where the blends do not end so.
bool Attach::notch(int ka, int ea, int kb, int eb)
{
  const int     k[2] = {ka, kb}, end[2] = {ea, eb};
  TopoDS_Vertex vertex, other;
  TopExp::Vertices(TopoDS::Edge(pieces[ka].rail.edge.Oriented(TopAbs_FORWARD)), vertex, other);
  if (ea == 1)
    vertex = other;
  const auto& atVertex = edgesOfVertex.FindFromKey(vertex);
  if (atVertex.Extent() != 4)
    return false;

  // On each other edge: the faces of both blends it lies between, and who wins there.
  struct Notched
  {
    TopoDS_Edge   edge;
    int           side[2] = {-1, -1}; // each blend's face on it
    int           winner  = -1;
    Meeting       cross; // the winner's contact across the edge (or, running on, across onWinner)
    double        q = 0; // the loser's contact where it meets the winner's sheet
    TopoDS_Vertex X, Q;
    std::vector<TopoDS_Edge> along; // the winner's end from its corner to Q
    // Where the winner runs on along the edge: the edge's far vertex, the face across there,
    // the edges of the winner's and the loser's faces to it, and the sheet's meeting with the
    // latter (its parameter there), past which it grows; the winner's corner Y is on the first.
    bool          runs = false;
    TopoDS_Vertex far, Y;
    TopoDS_Face   stop;
    TopoDS_Edge   onWinner, onLoser;
    double        onLoserAt = 0;
  } at[2];

  int count = 0;
  for (const auto& shape : atVertex)
  {
    if (railEdge(shape))
      continue;
    if (count == 2)
      return false;
    Notched& e = at[count++];
    e.edge     = TopoDS::Edge(shape);
    for (int p = 0; p < 2; ++p)
      for (int s = 0; s < 2; ++s)
        if (onFace(shape, pieces[k[p]].rail.s[s].face))
          e.side[p] = s;
    if (e.side[0] < 0 || e.side[1] < 0)
      return false;
    // Each blend's contact on its face there crosses the edge; the edge is cut back to the
    // crossing farther from the vertex, whose blend wins there, and the other contact ends
    // where it meets the winner's sheet (the nearer crossing lies in what the blends replace).
    std::optional<Meeting> crosses[2];
    double                 reach[2] = {-1, -1};
    const gp_Pnt           V        = BRep_Tool::Pnt(vertex);
    for (int p = 0; p < 2; ++p)
    {
      const Piece& piece  = pieces[k[p]];
      const auto [u0, u1] = halfOf(piece, end[p]);
      crosses[p]          = crossing(piece.fitted.contact[e.side[p]],
                            u0,
                            u1,
                            piece.uAt(end[p]),
                            e.edge,
                            piece.rail.s[e.side[p]].face);
      if (crosses[p])
        reach[p] =
          piece.sheet.contact(e.side[p])->Value(piece.sheet.rail(crosses[p]->first)).Distance(V);
    }
    // A contact that runs along the edge from its own blended edge does not cross it.
    for (int p = 0; p < 2; ++p)
    {
      gp_Vec mine, theirs;
      if (!crosses[p] && leaving(pieces[k[p]].rail.edge, vertex, mine, end[p])
          && leaving(e.edge, vertex, theirs) && tangent(mine, theirs))
      {
        if (e.runs)
          return false;
        e.runs   = true;
        e.winner = p;
      }
    }
    if (!e.runs)
    {
      if (reach[0] < 0 && reach[1] < 0)
        return false;
      e.winner = reach[0] > reach[1] ? 0 : 1;
    }
    const auto meets =
      meetsSheet(k[1 - e.winner], end[1 - e.winner], e.side[1 - e.winner], k[e.winner]);
    if (!meets)
      return false;
    e.q = *meets;
    if (!e.runs)
    {
      e.cross = *crosses[e.winner];
      continue;
    }
    // The edge's far vertex, with three edges: it, and one of each face to the face across.
    const int          w      = e.winner;
    const Piece&       winner = pieces[k[w]];
    const TopoDS_Face& mine   = winner.rail.s[e.side[w]].face;
    const TopoDS_Face& theirs = pieces[k[1 - w]].rail.s[e.side[1 - w]].face;
    TopoDS_Vertex      a, z;
    TopExp::Vertices(e.edge, a, z);
    e.far             = a.IsSame(vertex) ? z : a;
    const auto& atFar = edgesOfVertex.FindFromKey(e.far);
    if (atFar.Extent() != 3)
      return false;
    for (const auto& shape : atFar)
    {
      if (shape.IsSame(e.edge) || railEdge(shape))
        continue;
      if (onFace(shape, mine))
        e.onWinner = TopoDS::Edge(shape);
      else if (onFace(shape, theirs))
        e.onLoser = TopoDS::Edge(shape);
    }
    if (e.onWinner.IsNull() || e.onLoser.IsNull())
      return false;
    for (const auto& face : facesOfEdge.FindFromKey(e.onWinner))
      if (!face.IsSame(mine))
        e.stop = TopoDS::Face(face);
    if (e.stop.IsNull() || !onFace(e.onLoser, e.stop))
      return false;
    // The winner's contact across its face's edge there, and its sheet across the loser's,
    // that edge carried on past the vertex where needed.
    const auto [u0, u1] = halfOf(winner, end[w]);
    const auto across =
      crossing(winner.fitted.contact[e.side[w]], u0, u1, winner.uAt(end[w]), e.onWinner, mine);
    if (!across)
      return false;
    e.cross = *across;
    double        f, l;
    GeomAPI_IntCS meet(basisOf(BRep_Tool::Curve(e.onLoser, f, l)), winner.sheet.surface);
    if (!meet.IsDone() || meet.NbPoints() == 0)
      return false;
    const gp_Pnt F    = BRep_Tool::Pnt(e.far);
    int          best = 1;
    for (int i = 2; i <= meet.NbPoints(); ++i)
      if (meet.Point(i).Distance(F) < meet.Point(best).Distance(F))
        best = i;
    double su, sv;
    meet.Parameters(best, su, sv, e.onLoserAt);
    if (meet.Point(best).Distance(F) > 4 * winner.rail.scale)
      return false;
  }
  if (count != 2 || at[0].side[0] == at[1].side[0] || at[0].side[1] == at[1].side[1])
    return false;
  // The edges of the ends, found before anything changes.
  for (auto& e : at)
  {
    const int    w = e.winner, l = 1 - w;
    const Piece& winner = pieces[k[w]];
    const Piece& loser  = pieces[k[l]];
    builder.MakeVertex(e.Q,
                       loser.sheet.contact(e.side[l])->Value(loser.sheet.rail(e.q)),
                       tolerance);
    const TopoDS_Face& face = loser.rail.s[e.side[l]].face;
    double             f, l2;
    const gp_Pnt2d     onContact =
      winner.sheet.uv(winner.sheet.rail(e.cross.first), winner.sheet.across[e.side[w]]);
    if (!e.runs)
    {
      builder.MakeVertex(e.X,
                         winner.sheet.contact(e.side[w])->Value(winner.sheet.rail(e.cross.first)),
                         tolerance);
      const gp_Pnt2d xOnFace =
        BRep_Tool::CurveOnSurface(e.edge, face, f, l2)->Value(e.cross.second);
      const auto met = meetingEdge(placedSheet(winner.sheet),
                                   onContact,
                                   placedFace(face),
                                   xOnFace,
                                   e.X,
                                   e.Q,
                                   winner.rail.scale);
      if (!met)
        return false;
      e.along = {*met};
      continue;
    }
    // Running on: across the face at the far end from Y to X, then down the loser's face to Q.
    builder.MakeVertex(e.Y,
                       winner.sheet.contact(e.side[w])->Value(winner.sheet.rail(e.cross.first)),
                       tolerance);
    builder.MakeVertex(e.X, BRepAdaptor_Curve(e.onLoser).Value(e.onLoserAt), tolerance);
    const gp_Pnt2d yOnStop =
      BRep_Tool::CurveOnSurface(e.onWinner, e.stop, f, l2)->Value(e.cross.second);
    const gp_Pnt2d xOnFace = BRep_Tool::CurveOnSurface(e.onLoser, face, f, l2)->Value(e.onLoserAt);
    const auto     top     = meetingEdge(placedSheet(winner.sheet),
                                 onContact,
                                 placedFace(e.stop),
                                 yOnStop,
                                 e.Y,
                                 e.X,
                                 winner.rail.scale);
    if (!top)
      return false;
    const auto down = meetingEdge(placedSheet(winner.sheet),
                                  onSheet(k[w], end[w], BRep_Tool::Pnt(e.X)),
                                  placedFace(face),
                                  xOnFace,
                                  e.X,
                                  e.Q,
                                  winner.rail.scale);
    if (!down)
      return false;
    e.along = {*top, *down};
  }
  // Where the sheets meet, between the two places the losers' contacts end: across b's sheet
  // where b lost on both edges, else across a's.
  const auto onLoser = [&](const Notched& e, int p) {
    const Piece& piece = pieces[k[p]];
    return p == e.winner ? onSheet(k[p], end[p], BRep_Tool::Pnt(e.Q))
                         : piece.sheet.uv(piece.sheet.rail(e.q), piece.sheet.across[e.side[p]]);
  };
  const int    across       = at[0].winner == 0 && at[1].winner == 0 ? 1 : 0;
  const Piece& one          = pieces[k[across]];
  const Piece& two          = pieces[k[1 - across]];
  const bool   fittedSheets = !conical(one.sheet) || !conical(two.sheet);
  const auto   middle       = fittedSheets ? marchedMeeting(one.sheet,
                                                    onLoser(at[0], across),
                                                    onLoser(at[1], across),
                                                    two.sheet,
                                                    onLoser(at[0], 1 - across),
                                                    onLoser(at[1], 1 - across),
                                                    at[0].Q,
                                                    at[1].Q)
                                           : meetingEdge(placedSheet(one.sheet),
                                                 onLoser(at[0], across),
                                                 placedSheet(two.sheet),
                                                 onLoser(at[0], 1 - across),
                                                 at[0].Q,
                                                 at[1].Q,
                                                 one.rail.scale);
  if (!middle)
    return false;
  // Now the corners, the cut edges, the faces' new edges and each blend's end.
  for (auto& e : at)
  {
    const int          w = e.winner, l = 1 - w;
    const TopoDS_Face& face         = pieces[k[l]].rail.s[e.side[l]].face;
    corner[k[w]][e.side[w]][end[w]] = e.runs ? e.Y : e.X;
    u[k[w]][e.side[w]][end[w]]      = e.cross.first;
    corner[k[l]][e.side[l]][end[l]] = e.Q;
    u[k[l]][e.side[l]][end[l]]      = e.q;
    editOf(edits, face).bridges.push_back(e.along.back());
    notches.push_back({k[w], e.along.back(), face});
    if (!e.runs)
    {
      cut(e.edge, vertex, e.X, e.cross.second);
      continue;
    }
    // The edge and its far vertex go; the face there is cut back to Y and grows to X.
    for (const auto& f : facesOfEdge.FindFromKey(e.edge))
      editOf(edits, f).replace.Bind(e.edge, TopoDS_Edge());
    swept.push_back(e.edge);
    swept.push_back(e.far);
    cut(e.onWinner, e.far, e.Y, e.cross.second);
    cut(e.onLoser, e.far, e.X, e.onLoserAt);
    editOf(edits, e.stop).bridges.push_back(e.along.front());
    notches.push_back({k[w], e.along.front(), e.stop});
  }
  for (int p = 0; p < 2; ++p)
  {
    // From side 0's corner to side 1's: side 0's own edge where the blend won there, the
    // sheets' meeting, side 1's own edge where it won there.
    const Notched&           on0 = at[0].side[p] == 0 ? at[0] : at[1];
    const Notched&           on1 = at[0].side[p] == 0 ? at[1] : at[0];
    std::vector<TopoDS_Edge> ends;
    if (on0.winner == p)
      ends.insert(ends.end(), on0.along.begin(), on0.along.end());
    ends.push_back(*middle);
    if (on1.winner == p)
      ends.insert(ends.end(), on1.along.rbegin(), on1.along.rend());
    section[k[p]][end[p]] = *middle;
    if (ends.size() > 1)
      chain[k[p]][end[p]] = ends;
  }
  return true;
}

//! A stop where the blend's end is not one face's section: at a vertex of more edges, or where a
//! contact sweeps over short edges of its support before it crosses one. Out from the vertex
//! along each support's boundary, the edges the contact sweeps over go, and the faces across
//! them stop the blend, in turn with the faces around the vertex between the two supports' first
//! ones (and those around each swept vertex). The edges between consecutive stop faces are cut
//! where the sheet crosses them, the crossed edges at the contacts, and the blend's end runs along
//! the sheet on each stop face in turn.
TopoDS_Face Attach::faceAcross(const TopoDS_Shape& edge, const TopoDS_Shape& face) const
{
  for (const auto& other : facesOfEdge.FindFromKey(edge))
    if (!other.IsSame(face))
      return TopoDS::Face(other);
  return TopoDS_Face();
}

void Attach::facesAround(const Rail&               rail,
                         const TopoDS_Vertex&      at,
                         const TopoDS_Face&        first,
                         const TopoDS_Edge&        firstEdge,
                         const TopoDS_Edge&        lastEdge,
                         std::vector<TopoDS_Face>& faces,
                         std::vector<TopoDS_Edge>& edges)
{
  const auto support = [&](const TopoDS_Shape& face) {
    return face.IsSame(rail.s[0].face) || face.IsSame(rail.s[1].face);
  };
  const auto refuse = [&](const char* why) { fail(op, unsupported_geometry, why, {rail.edge}); };
  // Where an edge, as `face` uses it (a seam's two uses run on its two curves), reaches `at` in
  // the face's parameters.
  const auto reach = [&](const TopoDS_Edge&   use,
                         const TopoDS_Face&   face,
                         const TopoDS_Vertex& at) -> std::optional<gp_Pnt2d> {
    TopoDS_Vertex a, z;
    TopExp::Vertices(use, a, z);
    double     f, l;
    const auto pcurve = BRep_Tool::CurveOnSurface(use, face, f, l);
    if (BRep_Tool::Degenerated(use) || a.IsSame(z) || pcurve.IsNull())
      return {};
    return pcurve->Value(BRep_Tool::Parameter(at, use));
  };
  const auto wrong = [&] {
    refuse("the faces around an end of the edge do not run from one support to the other");
  };
  TopoDS_Face current  = first;
  TopoDS_Edge previous = firstEdge;
  // Each edge at the vertex is passed at most once on each of its sides.
  const int uses = 2 * edgesOfVertex.FindFromKey(at).Extent();
  for (int step = 0;; ++step)
  {
    const auto from = reach(previous, current, at);
    if (!from || step > uses)
      wrong();
    std::optional<std::pair<double, TopoDS_Edge>> next;
    for (const auto& shape : edgesOfVertex.FindFromKey(at))
    {
      if (!onFace(shape, current))
        continue;
      const TopoDS_Edge edge = TopoDS::Edge(shape);
      const bool        seam = BRep_Tool::IsClosed(edge, current);
      for (const auto orientation : {TopAbs_FORWARD, TopAbs_REVERSED})
      {
        const TopoDS_Edge use = TopoDS::Edge(edge.Oriented(orientation));
        if (use.IsSame(previous) && (!seam || use.Orientation() == previous.Orientation()))
          continue;
        if (const auto there = reach(use, current, at);
            there && (!next || there->Distance(*from) < next->first))
          next = {there->Distance(*from), use};
        if (!seam)
          break;
      }
    }
    if (!next)
      wrong();
    const TopoDS_Edge& edge = next->second;
    // Both uses end within their edge's tolerance of their curves and the vertex's of
    // them, so they meet in the face's parameters within what twice those span there.
    const double apart = 2
                         * (BRep_Tool::Tolerance(at)
                            + std::max(BRep_Tool::Tolerance(previous), BRep_Tool::Tolerance(edge)));
    const BRepAdaptor_Surface surface(current, false);
    if (!(next->first <= std::hypot(surface.UResolution(apart), surface.VResolution(apart))))
      wrong();
    if (edge.IsSame(lastEdge))
      return;
    if (railEdge(edge))
      wrong();
    const bool seam = BRep_Tool::IsClosed(edge, current);
    current         = seam ? current : faceAcross(edge, current);
    if (current.IsNull() || support(current))
      wrong();
    previous = seam ? TopoDS::Edge(edge.Reversed()) : edge;
    edges.push_back(previous);
    faces.push_back(current);
  }
}

bool Attach::endsOnFace(const Piece&       piece,
                        int                side,
                        int                end,
                        double             middle,
                        const StopWalk&    w,
                        const TopoDS_Face& own) const
{
  const Rail& rail = piece.rail;
  const double entry  = w.passed.empty() ? middle : w.passed.back().at.first;
  bool         leaves = false, rollsOn = false;
  for (TopExp_Explorer x(own, TopAbs_EDGE); x.More() && !leaves && !rollsOn; x.Next())
  {
    const TopoDS_Edge& e = TopoDS::Edge(x.Current());
    if (railEdge(e))
    {
      rollsOn = !ownEdge(piece, e) && onFace(e, rail.s[1 - side].face);
      continue;
    }
    if (!w.passed.empty() && e.IsSame(w.passed.back().edge))
      continue;
    leaves = crossing(piece.fitted.contact[side],
                      end == 0 ? 0 : entry,
                      end == 0 ? entry : 1,
                      piece.uAt(end),
                      e,
                      own)
               .has_value();
  }
  return !leaves && !rollsOn;
}

void Attach::walkSupport(StopEnd& st, int side)
{
  const int            end    = st.end;
  const TopoDS_Vertex& vertex = st.vertex;
  const Piece&         piece  = pieces[st.k];
  const Rail&          rail   = piece.rail;
  const Sheet&         sheet  = piece.sheet;
  const double         middle = st.middle;
  const double (&carried)[2]  = st.carried;
  const auto support          = [&](const TopoDS_Shape& face) {
    return face.IsSame(rail.s[0].face) || face.IsSame(rail.s[1].face);
  };
  const auto across = [&](const TopoDS_Shape& edge, const TopoDS_Shape& face) {
    return faceAcross(edge, face);
  };
  const auto refuse = [&](const char* why) { fail(op, unsupported_geometry, why, {rail.edge}); };
  // Whether a side's contact, where it is carried at this end, meets a face's surface: a face
  // ends the blend only where its surface cuts the whole section, from one contact to the other.
  const auto meets = [&](int side, const TopoDS_Face& face) {
    GeomAPI_IntCS meet(new Geom_TrimmedCurve(sheet.contact(side),
                                             std::min(carried[0], carried[1]),
                                             std::max(carried[0], carried[1])),
                       BRep_Tool::Surface(face));
    return meet.IsDone() && meet.NbPoints() > 0;
  };
  StopWalk&   w   = st.walk[side];
  TopoDS_Face own = rail.s[side].face;
  TopoDS_Edge edge;
  for (const auto& shape : edgesOfVertex.FindFromKey(vertex))
    if (!ownEdge(piece, shape) && onFace(shape, own))
      edge = TopoDS::Edge(shape);
  if (edge.IsNull() || railEdge(edge))
    refuse(
      "where blended edges meet other edges, the corner is blended too, which the blend does not "
      "do yet");
  TopoDS_Vertex from        = vertex;
  size_t        facesBefore = 0, betweenBefore = 0;
  // Each edge of the support is swept over at most once.
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> ofOwn;
  TopExp::MapShapes(own, TopAbs_EDGE, ofOwn);
  for (int step = 0, onOwn = 0;; ++step, ++onOwn)
  {
    if (onOwn > ofOwn.Extent())
      refuse("the blend's contact does not reach the faces that stop it");
    const TopoDS_Face beyond = across(edge, own);
    // A split may have reached the next edge between the same supports.
    // There is no stopping face there; continue only that required rail.
    if (!beyond.IsNull() && beyond.IsSame(rail.s[1 - side].face) && step == 0
        && piece.joint[1 - end].meets == Meets::Split && !railEdge(edge))
    {
      gp_Vec mine, next;
      if (leaving(rail.edge, vertex, mine) && leaving(edge, vertex, next) && tangent(mine, next))
        throw ContinuationRequired(op, rail.edge, edge);
    }
    if (beyond.IsNull() || support(beyond))
      refuse("the faces at an end of the edge do not stop it");
    if (w.faces.empty())
      w.faces.push_back(beyond);
    const auto found = crossing(piece.fitted.contact[side],
                                end == 0 ? 0 : middle,
                                end == 0 ? middle : 1,
                                piece.uAt(end),
                                edge,
                                own);
    // A sharp stop surface can meet the contact past this edge's vertex. Searching
    // the whole carried contact is necessary: a stop may consume over half its span.
    if (!found && step == 0 && meets(1 - side, beyond))
    {
      if (auto more = continuation(piece.fitted.contact[side],
                                   0,
                                   1,
                                   piece.uAt(end),
                                   edge,
                                   vertex,
                                   own,
                                   beyond,
                                   tolerance))
      {
        w.faces     = {beyond};
        w.up        = edge;
        w.at        = {more->u, more->from, {}};
        w.continued = std::move(more);
        break;
      }
    }
    // The support goes on smoothly past the vertex into the face across this edge (two
    // faces of one wall), and the contact runs on beside it: the face across that face's
    // other edge at the vertex stops the blend, its surface carried to where the contact
    // meets it, as a single stop face's is (continuation).
    if (!found && step == 0 && smoothAt(own, beyond, edge, vertex))
    {
      TopoDS_Edge up;
      for (const auto& shape : edgesOfVertex.FindFromKey(vertex))
        if (!shape.IsSame(edge) && onFace(shape, beyond))
          up = TopoDS::Edge(shape);
      const TopoDS_Face stopFace = up.IsNull() ? TopoDS_Face() : across(up, beyond);
      if (!stopFace.IsNull() && !support(stopFace) && meets(1 - side, stopFace))
        if (auto more = continuation(piece.fitted.contact[side],
                                     end == 0 ? 0 : middle,
                                     end == 0 ? middle : 1,
                                     piece.uAt(end),
                                     up,
                                     vertex,
                                     own,
                                     stopFace,
                                     tolerance))
        {
          w.faces     = {stopFace};
          w.up        = up;
          w.at        = {more->u, more->from, {}};
          w.continued = std::move(more);
          break;
        }
    }
    if (found)
    {
      // The face across the crossed edge goes on from the support (a niche's floor level
      // with it): the blend does not stop on it, but on the face before, whose surface
      // the crossed edge lies in; from the vertex to the crossing the edge is then that
      // face's with the face across.
      if (step > 0 && smoothAt(own, beyond, edge, from))
      {
        // A face of the support's own surface, where the contact crosses onto it short of
        // the face before (the crossed edge does not lie in that face's surface): the blend's
        // sheet lies on it as on the support, and the contact runs on over it to the faces
        // that stop it, as over the support (the faces around the vertex up to it are those).
        // The two faces' surfaces are one where they give the same points at the contact's
        // parameters from the crossing to the end of the carried contact (copies of one
        // surface, as a Boolean leaves the parts of a face).
        const auto&  contact    = piece.fitted.contact[side];
        const auto   ownSurface = BRep_Tool::Surface(own), nextSurface = BRep_Tool::Surface(beyond);
        const double away = end == 0 ? 0 : 1;
        bool         one  = true;
        for (int i = 0; i <= 16 && one; ++i)
        {
          const gp_Pnt2d uv = contact->Value(found->first + (away - found->first) * i / 16);
          one = ownSurface->Value(uv.X(), uv.Y()).Distance(nextSurface->Value(uv.X(), uv.Y()))
                <= tol::length;
        }
        const gp_Pnt2d             onOwnUV = contact->Value(found->first);
        const gp_Pnt               at      = ownSurface->Value(onOwnUV.X(), onOwnUV.Y());
        GeomAPI_ProjectPointOnSurf onBefore(at, BRep_Tool::Surface(w.faces[facesBefore - 1]));
        if (one && w.faces.size() > facesBefore && w.faces.back().IsSame(beyond)
            && (onBefore.NbPoints() == 0
                || onBefore.LowerDistance() > rail.resolved(limit::together)))
        {
          w.passed.push_back({edge, from, *found, beyond});
          own  = beyond;
          edge = w.between.back();
          w.faces.pop_back();
          w.between.pop_back();
          w.from.pop_back();
          ofOwn.Clear();
          TopExp::MapShapes(own, TopAbs_EDGE, ofOwn);
          onOwn = -1;
          continue;
        }
        w.faces.resize(facesBefore);
        w.between.resize(betweenBefore);
        w.from.resize(betweenBefore);
        w.level = beyond;
      }
      w.crossed     = edge;
      w.crossedFrom = from;
      w.at          = *found;
      break;
    }
    // Swept over: on to the next edge of the support at its other end.
    TopoDS_Vertex a, z;
    TopExp::Vertices(edge, a, z);
    const TopoDS_Vertex at = a.IsSame(from) ? z : a;
    TopoDS_Edge         next;
    for (const auto& shape : edgesOfVertex.FindFromKey(at))
      if (!shape.IsSame(edge) && onFace(shape, own))
        next = TopoDS::Edge(shape);
    // Past the stop faces to another edge asked between the same two faces (the next piece
    // of an edge a hole is drilled through): the stop faces do not stop the ball, which
    // rolls on over them to that piece (liftFeatures).
    if (!next.IsNull() && railEdge(next) && onFace(next, rail.s[1 - side].face))
      throw Collision(op,
                      "the blend's contact does not reach the faces that stop it",
                      rail.edge,
                      std::vector<TopoDS_Shape>(w.faces.begin(), w.faces.end()));
    // The contact ends on the face: it crosses none of the face's edges (but the rail's and
    // the one it came onto the face over) where it is carried, and no other piece of the
    // blended edges lies on the face beside the other support (over which the ball rolls on,
    // above). The face that stops the other contact then stops it too, its surface carried
    // to it.
    if (endsOnFace(piece, side, end, middle, w, own))
    {
      w.open     = own;
      w.openEdge = edge;
      w.openFrom = from;
      break;
    }
    w.swept.push_back(edge);
    w.sweptAt.push_back(at);
    if (next.IsNull() || railEdge(next))
      refuse("the blend's contact does not reach the faces that stop it");
    const TopoDS_Face nextBeyond = across(next, own);
    if (nextBeyond.IsNull() || support(nextBeyond))
      refuse("the faces at an end of the edge do not stop it");
    const size_t before = w.between.size();
    facesBefore         = w.faces.size();
    betweenBefore       = before;
    facesAround(rail, at, w.faces.back(), edge, next, w.faces, w.between);
    w.from.insert(w.from.end(), w.between.size() - before, at);
    edge = next;
    from = at;
  }
}

void Attach::stopFaces(StopEnd& st)
{
  const TopoDS_Vertex& vertex = st.vertex;
  const Piece&         piece  = pieces[st.k];
  const Rail&          rail   = piece.rail;
  const auto&          walk   = st.walk;
  const auto refuse = [&](const char* why) { fail(op, unsupported_geometry, why, {rail.edge}); };
  // All the stop faces in order from side 0's crossing to side 1's.
  std::vector<TopoDS_Face>& faces = st.faces;
  faces.assign(walk[0].faces.rbegin(), walk[0].faces.rend());
  std::vector<TopoDS_Edge>& between = st.between;
  between.assign(walk[0].between.rbegin(), walk[0].between.rend());
  // Side 0's walk ran the other way: each seam among its edges is used from its other side.
  for (size_t i = 0; i < between.size(); ++i)
    if (BRep_Tool::IsClosed(between[i], faces[i]))
      between[i].Reverse();
  std::vector<TopoDS_Vertex>& from = st.from;
  from.assign(walk[0].from.rbegin(), walk[0].from.rend());
  {
    std::vector<TopoDS_Face> ring;
    std::vector<TopoDS_Edge> edges;
    TopoDS_Edge              beside[2];
    for (const auto& shape : edgesOfVertex.FindFromKey(vertex))
      for (int side = 0; side < 2; ++side)
        if (!ownEdge(piece, shape) && onFace(shape, rail.s[side].face))
          beside[side] = TopoDS::Edge(shape);
    if (walk[0].continued)
      beside[0] = walk[0].up;
    facesAround(rail, vertex, walk[0].faces.front(), beside[0], beside[1], ring, edges);
    faces.insert(faces.end(), ring.begin(), ring.end());
    between.insert(between.end(), edges.begin(), edges.end());
    from.insert(from.end(), edges.size(), vertex);
  }
  faces.insert(faces.end(), walk[1].faces.begin() + 1, walk[1].faces.end());
  between.insert(between.end(), walk[1].between.begin(), walk[1].between.end());
  from.insert(from.end(), walk[1].from.begin(), walk[1].from.end());
  if (faces.size() != between.size() + 1)
    refuse("the faces that stop the blend do not run from one support to the other");
}

void Attach::openCorner(StopEnd& st, int o, const Handle(Geom_Surface)& stopSurface)
{
  const int    end                    = st.end;
  const Piece& piece                  = pieces[st.k];
  const Rail&  rail                   = piece.rail;
  const Sheet& sheet                  = piece.sheet;
  const double (&carried)[2]          = st.carried;
  auto&        walk                   = st.walk;
  const auto refuse = [&](const char* why) { fail(op, unsupported_geometry, why, {rail.edge}); };
  // Where the open contact meets the stop face's surface: first past where it came onto its
  // face, on that face.
  const double entry = walk[o].passed.empty() ? piece.uAt(end) : walk[o].passed.back().at.first;
  std::optional<double> uOpen;
  {
    GeomAPI_IntCS meet(new Geom_TrimmedCurve(sheet.contact(o),
                                             std::min(carried[0], carried[1]),
                                             std::max(carried[0], carried[1])),
                       stopSurface);
    for (int p = 1; meet.IsDone() && p <= meet.NbPoints(); ++p)
    {
      double su, sv, w;
      meet.Parameters(p, su, sv, w);
      const double at = (w - sheet.alpha) / sheet.beta;
      if ((end == 0 ? at < entry : at > entry)
          && (!uOpen || std::abs(at - entry) < std::abs(*uOpen - entry)))
        uOpen = at;
    }
  }
  if (!uOpen || Inside(walk[o].open).out(piece.fitted.contact[o]->Value(*uOpen)))
    refuse("the faces at an end of the edge do not stop it");
  // Its contact in that face's parameters: by whole periods near the face's edge before it.
  {
    double     f, l;
    const auto pcurve = BRep_Tool::CurveOnSurface(walk[o].openEdge, walk[o].open, f, l);
    if (pcurve.IsNull())
      refuse("the faces at an end of the edge do not stop it");
    const gp_Pnt2d uv = piece.fitted.contact[o]->Value(*uOpen);
    const gp_Pnt2d near =
      nearPeriod(uv,
                 pcurve->Value(BRep_Tool::Parameter(walk[o].openFrom, walk[o].openEdge)),
                 BRep_Tool::Surface(walk[o].open));
    walk[o].at = {*uOpen, 0, gp_Vec2d(near, uv)};
  }
}

void Attach::openStop(StopEnd& st)
{
  const int                   k       = st.k;
  const int                   end     = st.end;
  const Piece&                piece   = pieces[k];
  const Rail&                 rail    = piece.rail;
  const Sheet&                sheet   = piece.sheet;
  auto&                       walk    = st.walk;
  std::vector<TopoDS_Face>&   faces   = st.faces;
  std::vector<TopoDS_Edge>&   between = st.between;
  std::vector<TopoDS_Vertex>& from    = st.from;
  const auto refuse = [&](const char* why) { fail(op, unsupported_geometry, why, {rail.edge}); };
  // One contact ends on its face: the face the other crosses onto stops the blend, if its
  // surface meets this contact too, and grows to the blend. The faces from it round to the
  // open contact's face are cut by its surface: its new edges with them run from its own
  // vertex on the first to the open contact's corner, and what lies under the blend goes.
  const int o = walk[0].open.IsNull() ? 1 : 0, s = 1 - o;
  if (!walk[s].open.IsNull() || walk[s].crossed.IsNull() || walk[0].continued || walk[1].continued
      || !walk[s].level.IsNull())
    refuse("the faces at an end of the edge do not stop it");
  if (s == 1)
  {
    std::reverse(faces.begin(), faces.end());
    std::reverse(between.begin(), between.end());
    std::reverse(from.begin(), from.end());
  }
  faces.push_back(walk[o].open);
  between.push_back(walk[o].openEdge);
  from.push_back(walk[o].openFrom);
  for (size_t i = 0; i < between.size(); ++i)
    if (BRep_Tool::IsClosed(between[i], faces[i]) || BRep_Tool::IsClosed(between[i], faces[i + 1]))
      refuse("the faces that stop the blend do not run from one support to the other");
  const TopoDS_Face& stopFace    = faces.front();
  const auto         stopSurface = placedFace(stopFace).surface;
  openCorner(st, o, stopSurface);
  for (int side = 0; side < 2; ++side)
  {
    u[k][side][end] = walk[side].at.first;
    builder.MakeVertex(corner[k][side][end],
                       sheet.contact(side)->Value(sheet.rail(walk[side].at.first)),
                       tolerance);
  }
  const auto onSurface = [](const gp_Pnt& p, const Placed& placed, gp_Pnt2d& uv) {
    GeomAPI_ProjectPointOnSurf project(p, placed.surface);
    if (project.NbPoints() == 0 || project.LowerDistance() > tol::length)
      return false;
    double pu, pv;
    project.LowerDistanceParameters(pu, pv);
    uv.SetCoord(pu, pv);
    return true;
  };
  // The section on the stop face, from side 0's corner to side 1's.
  gp_Pnt2d onStop;
  if (!onSurface(BRep_Tool::Pnt(corner[k][0][end]), placedFace(stopFace), onStop))
    refuse("the faces at an end of the edge do not stop it");
  auto ended = meetingEdge(placedSheet(sheet),
                           sheet.uv(sheet.rail(u[k][0][end]), sheet.across[0]),
                           placedFace(stopFace),
                           onStop,
                           corner[k][0][end],
                           corner[k][1][end],
                           rail.scale);
  if (refines(!ended))
    ended = planeMeeting(sheet,
                         stopFace,
                         sheet.uv(sheet.rail(u[k][0][end]), sheet.across[0]),
                         sheet.uv(sheet.rail(u[k][1][end]), sheet.across[1]),
                         corner[k][0][end],
                         corner[k][1][end]);
  if (!ended)
    fail(op,
         no_solution,
         "the blend's end on the face that stops it does not reach both contacts",
         {rail.edge});
  // The stop face's new edges: from its vertex on the first edge between, through where its
  // surface crosses each further edge between (each cut back to it), to the open corner.
  std::vector<TopoDS_Vertex> through = {from[0]};
  for (size_t i = 1; i < between.size(); ++i)
  {
    double     f, l;
    const auto curve = BRep_Tool::Curve(between[i], f, l);
    if (curve.IsNull())
      refuse("the faces at an end of the edge do not stop it");
    GeomAPI_IntCS                            meet(new Geom_TrimmedCurve(curve, f, l), stopSurface);
    const gp_Pnt                             near = BRep_Tool::Pnt(from[i]);
    std::optional<std::pair<gp_Pnt, double>> crossed;
    for (int p = 1; meet.IsDone() && p <= meet.NbPoints(); ++p)
    {
      double su, sv, w;
      meet.Parameters(p, su, sv, w);
      if (!crossed || meet.Point(p).Distance(near) < crossed->first.Distance(near))
        crossed = {meet.Point(p), w};
    }
    if (!crossed)
      refuse("the faces at an end of the edge do not stop it");
    TopoDS_Vertex v;
    builder.MakeVertex(v, crossed->first, tolerance);
    cut(between[i], from[i], v, crossed->second);
    through.push_back(v);
  }
  through.push_back(corner[k][o][end]);
  for (size_t i = 0; i + 1 < through.size(); ++i)
  {
    gp_Pnt2d     a, b;
    const gp_Pnt p = BRep_Tool::Pnt(through[i]);
    if (!onSurface(p, placedFace(stopFace), a) || !onSurface(p, placedFace(faces[i + 1]), b))
      refuse("the faces at an end of the edge do not stop it");
    const auto grown = meetingEdge(placedFace(stopFace),
                                   a,
                                   placedFace(faces[i + 1]),
                                   b,
                                   through[i],
                                   through[i + 1],
                                   rail.scale);
    if (!grown)
      fail(op,
           no_solution,
           "the face that stops the blend does not grow to the contact that ends before it",
           {rail.edge});
    editOf(edits, stopFace).bridges.push_back(*grown);
    editOf(edits, faces[i + 1]).bridges.push_back(*grown);
    continued.push_back({*grown, {stopFace, faces[i + 1]}});
  }
  // The crossed edge, the edges each contact crossed onto faces of its surface, and what the
  // contacts swept over; the stop face's vertex on the first edge between stays.
  cut(walk[s].crossed, walk[s].crossedFrom, corner[k][s][end], walk[s].at.second);
  for (int side = 0; side < 2; ++side)
  {
    const auto& passed = walk[side].passed;
    for (size_t j = 0; j < passed.size(); ++j)
    {
      TopoDS_Vertex at;
      builder.MakeVertex(at, sheet.contact(side)->Value(sheet.rail(passed[j].at.first)), tolerance);
      cut(passed[j].edge, passed[j].from, at, passed[j].at.second);
      passages[{k, side, end}].push_back(
        {at,
         passed[j].at.first,
         passed[j].face,
         j + 1 < passed.size() ? passed[j + 1].at.shift : walk[side].at.shift});
    }
    for (const auto& edge : walk[side].swept)
    {
      for (const auto& face : facesOfEdge.FindFromKey(edge))
        editOf(edits, face).replace.Bind(edge, TopoDS_Edge());
      swept.push_back(edge);
    }
    for (const auto& v : walk[side].sweptAt)
      if (!v.IsSame(from[0]))
        swept.push_back(v);
  }
  editOf(edits, stopFace).bridges.push_back(*ended);
  stopAt[k][end]  = stopFace;
  section[k][end] = *ended;
}

void Attach::stopPlaces(StopEnd& st)
{
  const int                         k       = st.k;
  const int                         end     = st.end;
  const Piece&                      piece   = pieces[k];
  const Rail&                       rail    = piece.rail;
  const Sheet&                      sheet   = piece.sheet;
  const auto&                       walk    = st.walk;
  const std::vector<TopoDS_Face>&   faces   = st.faces;
  const std::vector<TopoDS_Edge>&   between = st.between;
  const std::vector<TopoDS_Vertex>& from    = st.from;
  const auto refuse = [&](const char* why) { fail(op, unsupported_geometry, why, {rail.edge}); };
  // The places: the two corners, and where the sheet crosses each edge between stop faces.
  std::vector<TopoDS_Vertex>& at    = st.at;
  std::vector<double>&        param = st.param;
  at.resize(faces.size() + 1);
  param.resize(faces.size() + 1);
  for (int side = 0; side < 2; ++side)
  {
    u[k][side][end] = walk[side].at.first;
    builder.MakeVertex(corner[k][side][end],
                       sheet.contact(side)->Value(sheet.rail(walk[side].at.first)),
                       tolerance);
    at[side == 0 ? 0 : at.size() - 1]    = corner[k][side][end];
    param[side == 0 ? 0 : at.size() - 1] = walk[side].at.second;
  }
  for (size_t i = 0; i < between.size(); ++i)
  {
    double     f, l;
    const auto curve = BRep_Tool::Curve(between[i], f, l);
    if (curve.IsNull())
      refuse("the blend does not cross an edge between the faces that stop it");
    // Where the sheet crosses the edge, nearest the end the blend replaces; or just past that
    // end, within the blend's reach, where the faces on either side grow to it.
    const gp_Pnt                             end0   = BRep_Tool::Pnt(from[i]);
    const double                             atFrom = BRep_Tool::Parameter(from[i], between[i]);
    std::optional<std::pair<gp_Pnt, double>> crossed;
    // As far past the end as the blend is carried, and twice its size more.
    gp_Pnt first, last, carriedFirst, carriedLast;
    rail.curve.D0(rail.first, first);
    rail.curve.D0(rail.last, last);
    rail.curve.D0(rail.from, carriedFirst);
    rail.curve.D0(rail.to, carriedLast);
    const double reach =
      std::max(first.Distance(carriedFirst), last.Distance(carriedLast)) + 2 * rail.scale;
    for (const bool past : {false, true})
    {
      GeomAPI_IntCS meet(past ? basisOf(curve)
                              : Handle(Geom_Curve)(new Geom_TrimmedCurve(curve, f, l)),
                         sheet.surface);
      for (int p = 1; meet.IsDone() && p <= meet.NbPoints(); ++p)
      {
        double su, sv, w;
        meet.Parameters(p, su, sv, w);
        if (past
            && !((w - atFrom) * (atFrom - (f + l) / 2) > 0
                 && meet.Point(p).Distance(end0) <= reach))
          continue;
        if (!crossed || meet.Point(p).Distance(end0) < crossed->first.Distance(end0))
          crossed = {meet.Point(p), w};
      }
      if (crossed)
        break;
    }
    if (!crossed)
      refuse("the blend does not cross an edge between the faces that stop it");
    builder.MakeVertex(at[i + 1], crossed->first, tolerance);
    param[i + 1] = crossed->second;
  }
  // A stop face that the end only touches, where a corner falls on an edge between stop faces
  // (a support's edge lands on a stop face's edge): no piece of the end runs on it, and the
  // crossing there is the corner.
  std::vector<bool>& touched = st.touched;
  touched.assign(faces.size(), false);
  for (size_t i = 0; i < faces.size(); ++i)
  {
    touched[i] =
      BRep_Tool::Pnt(at[i]).Distance(BRep_Tool::Pnt(at[i + 1])) <= rail.resolved(limit::together);
    if (!touched[i])
      continue;
    // One vertex for the place: the corner where it is one.
    if (i + 1 == faces.size())
      at[i] = at[i + 1];
    else
      at[i + 1] = at[i];
  }
  if (std::all_of(touched.begin(), touched.end(), [](bool t) { return t; }))
    refuse("the blend's end does not run on the faces that stop it");
}

void Attach::stopEdges(StopEnd& st)
{
  const int                         k       = st.k;
  const int                         end     = st.end;
  const Piece&                      piece   = pieces[k];
  const Rail&                       rail    = piece.rail;
  const Sheet&                      sheet   = piece.sheet;
  const auto&                       walk    = st.walk;
  const std::vector<TopoDS_Face>&   faces   = st.faces;
  const std::vector<TopoDS_Edge>&   between = st.between;
  const std::vector<TopoDS_Vertex>& at      = st.at;
  const std::vector<double>&        param   = st.param;
  const std::vector<bool>&          touched = st.touched;
  const auto refuse = [&](const char* why) { fail(op, unsupported_geometry, why, {rail.edge}); };
  // Along the sheet on each stop face, from one place to the next.
  std::vector<TopoDS_Edge>& made = st.made;
  made.resize(faces.size());
  for (size_t i = 0; i < faces.size(); ++i)
  {
    if (touched[i])
      continue;
    // Where it starts on the face: on the edge the contact crossed, or where it meets the
    // stop face's surface carried past the vertex.
    gp_Pnt2d start;
    if (i == 0 && walk[0].continued)
      start = walk[0].continued->onStop->Value(walk[0].continued->from);
    else
    {
      const TopoDS_Edge& on = i == 0 ? walk[0].crossed : between[i - 1];
      double             f, l;
      const auto         pcurve = BRep_Tool::CurveOnSurface(on, faces[i], f, l);
      if (pcurve.IsNull())
        refuse("the faces that stop the blend do not run from one support to the other");
      start = pcurve->Value(param[i]);
    }
    gp_Pnt2d onSheet;
    if (i == 0)
      onSheet = sheet.uv(sheet.rail(u[k][0][end]), sheet.across[0]);
    else
    {
      GeomAPI_ProjectPointOnSurf project(BRep_Tool::Pnt(at[i]), sheet.surface);
      double                     pu, pv;
      project.LowerDistanceParameters(pu, pv);
      onSheet =
        nearPeriod(gp_Pnt2d(pu, pv),
                   sheet.uv(sheet.rail(piece.uAt(end)), (sheet.across[0] + sheet.across[1]) / 2),
                   sheet.surface);
    }
    const auto met = meetingEdge(placedSheet(sheet),
                                 onSheet,
                                 placedFace(faces[i]),
                                 start,
                                 at[i],
                                 at[i + 1],
                                 rail.scale);
    if (!met)
      fail(op,
           no_solution,
           "the blend's end on the faces that stop it does not reach both contacts",
           {rail.edge});
    made[i] = *met;
  }
}

void Attach::stopCuts(StopEnd& st)
{
  const int                         k       = st.k;
  const int                         end     = st.end;
  const Piece&                      piece   = pieces[k];
  const Rail&                       rail    = piece.rail;
  const Sheet&                      sheet   = piece.sheet;
  const TopoDS_Vertex&              vertex  = st.vertex;
  const auto&                       walk    = st.walk;
  const std::vector<TopoDS_Face>&   faces   = st.faces;
  const std::vector<TopoDS_Edge>&   between = st.between;
  const std::vector<TopoDS_Vertex>& from    = st.from;
  const std::vector<TopoDS_Vertex>& at      = st.at;
  const std::vector<double>&        param   = st.param;
  const std::vector<bool>&          touched = st.touched;
  const std::vector<TopoDS_Edge>&   made    = st.made;
  const auto refuse = [&](const char* why) { fail(op, unsupported_geometry, why, {rail.edge}); };
  // Now the cuts, what the contacts swept over, and the faces' new edges.
  for (int side = 0; side < 2; ++side)
  {
    if (const auto& more = walk[side].continued)
    {
      // The new piece of the edge between the support and the stop face, from the corner
      // back to the vertex, in both faces' wires.
      const TopoDS_Face& own                    = rail.s[side].face;
      const TopoDS_Face& stopFace               = walk[side].faces.front();
      const auto [supportSurface, supportPlace] = surfaceOf(own);
      const auto [stopSurface, stopPlace]       = surfaceOf(stopFace);
      const TopoDS_Edge run                     = edgeBetween(
        more->curve,
        corner[k][side][end],
        more->from,
        vertex,
        more->to,
        {{more->onSupport, supportSurface, supportPlace}, {more->onStop, stopSurface, stopPlace}},
        tolerance);
      editOf(edits, own).bridges.push_back(run);
      editOf(edits, stopFace).bridges.push_back(run);
      continued.push_back({run, {own, stopFace}});
    }
    else
      cut(walk[side].crossed, walk[side].crossedFrom, corner[k][side][end], walk[side].at.second);
    // Each edge the contact crosses onto a face of its support's surface is cut back to the
    // crossing; the contact's part on that face lies in its parameters as the next crossing's.
    const auto& passed = walk[side].passed;
    for (size_t j = 0; j < passed.size(); ++j)
    {
      TopoDS_Vertex at;
      builder.MakeVertex(at, sheet.contact(side)->Value(sheet.rail(passed[j].at.first)), tolerance);
      cut(passed[j].edge, passed[j].from, at, passed[j].at.second);
      passages[{k, side, end}].push_back(
        {at,
         passed[j].at.first,
         passed[j].face,
         j + 1 < passed.size() ? passed[j + 1].at.shift : walk[side].at.shift});
    }
    if (const TopoDS_Face& level = walk[side].level; !level.IsNull())
    {
      // The crossed edge from the vertex the blend replaces to the corner, now between the
      // face level with the support and the stop face.
      const TopoDS_Face& stopFace    = side == 0 ? faces.front() : faces.back();
      const TopoDS_Edge& crossedEdge = walk[side].crossed;
      double             f, l, f1, l1, f2, l2;
      const auto         curve   = BRep_Tool::Curve(crossedEdge, f, l);
      const auto         onLevel = BRep_Tool::CurveOnSurface(crossedEdge, level, f1, l1);
      const auto         onStop  = BRep_Tool::CurveOnSurface(crossedEdge, stopFace, f2, l2);
      if (curve.IsNull() || onLevel.IsNull() || onStop.IsNull())
        refuse("the edge the contact crosses onto a face level with its support does not lie on "
               "the face that stops the blend");
      const auto [levelSurface, levelPlace] = surfaceOf(level);
      const auto [stopSurface, stopPlace]   = surfaceOf(stopFace);
      const TopoDS_Edge part =
        edgeBetween(curve,
                    walk[side].crossedFrom,
                    BRep_Tool::Parameter(walk[side].crossedFrom, crossedEdge),
                    corner[k][side][end],
                    walk[side].at.second,
                    {{onLevel, levelSurface, levelPlace}, {onStop, stopSurface, stopPlace}},
                    tolerance);
      editOf(edits, level).bridges.push_back(part);
      editOf(edits, stopFace).bridges.push_back(part);
      continued.push_back({part, {level, stopFace}});
    }
    for (const auto& edge : walk[side].swept)
    {
      for (const auto& face : facesOfEdge.FindFromKey(edge))
        editOf(edits, face).replace.Bind(edge, TopoDS_Edge());
      swept.push_back(edge);
    }
    swept.insert(swept.end(), walk[side].sweptAt.begin(), walk[side].sweptAt.end());
  }
  for (size_t i = 0; i < between.size(); ++i)
    cut(between[i], from[i], at[i + 1], param[i + 1]);
  std::vector<TopoDS_Edge> ends;
  for (size_t i = 0; i < faces.size(); ++i)
  {
    if (touched[i])
      continue;
    editOf(edits, faces[i]).bridges.push_back(made[i]);
    if (!ends.empty())
      notches.push_back({k, made[i], faces[i]});
    else
      stopAt[k][end] = faces[i];
    ends.push_back(made[i]);
  }
  section[k][end] = ends.front();
  if (ends.size() > 1)
    chain[k][end] = ends;
}

void Attach::stops(int k, int end, const TopoDS_Vertex& vertex)
{
  const Piece& piece = pieces[k];
  StopEnd      st{k, end, vertex};
  st.middle     = piece.uAt(0) / 2 + piece.uAt(1) / 2;
  st.carried[0] = piece.sheet.rail(end == 0 ? 0 : st.middle);
  st.carried[1] = piece.sheet.rail(end == 0 ? st.middle : 1);
  for (int side = 0; side < 2; ++side)
    walkSupport(st, side);
  // All the stop faces in order from side 0's crossing to side 1's.
  stopFaces(st);
  if (!st.walk[0].open.IsNull() || !st.walk[1].open.IsNull())
  {
    openStop(st);
    return;
  }
  // The places: the two corners, and where the sheet crosses each edge between stop faces.
  stopPlaces(st);
  // Along the sheet on each stop face, from one place to the next.
  stopEdges(st);
  // Now the cuts, what the contacts swept over, and the faces' new edges.
  stopCuts(st);
}

} // namespace BRepExt::blending::detail
