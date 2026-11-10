// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Blend: The blends attached to the body, stage by stage; the chains' ends and the assembly.

#include <BRepBlendExt_Attach.hxx>

namespace BRepExt::blending::detail
{

Attach::Attach(TopoDS_Shape                     body_,
               const std::vector<Piece>&        pieces_,
               const std::vector<Corner>&       cornersAt_,
               const Handle(BRepTools_History)& history_,
               bool refineContacts_,
               StageReads* reads_)
    : refineContacts(refineContacts_),
      reads(reads_),
      retainedBoundaries(pieces_.size()),
      body(std::move(body_)),
      pieces(pieces_),
      cornersAt(cornersAt_),
      history(history_),
      op(pieces_.front().rail.op),
      n(int(pieces_.size())),
      contactIntervals(n),
      section(n),
      stopAt(n),
      patches(cornersAt_.size()),
      parts(cornersAt_.size()),
      blends(n),
      chain(n)
{
  for (auto& sides : contactIntervals)
    for (auto& intervals : sides)
      intervals.emplace_back();
  TopExp::MapShapesAndUniqueAncestors(body, TopAbs_VERTEX, TopAbs_EDGE, edgesOfVertex);
  TopExp::MapShapesAndUniqueAncestors(body, TopAbs_EDGE, TopAbs_FACE, facesOfEdge);
  for (const auto& piece : pieces)
    for (const auto& joint : piece.joint)
      if (joint.meets == Meets::Split)
        crossed.Add(joint.crossed);
}

void Attach::splitNotchCut(int k, int end)
{
  const Piece& piece = pieces[k];
  const Rail&  rail  = piece.rail;
  const auto&  joint = piece.joint[end];
  const Piece& other = pieces[joint.piece];
  // A far edge the contact runs off over partway: its part along the face that stays
  // keeps the face's end there; the rest is inside the blend.
  const TopoDS_Edge& edge  = joint.crossed;
  const Piece&       plain = splitNotch(k, joint.side) ? other : piece;
  const int          away  = splitNotch(k, joint.side) ? 1 - joint.end : 1 - end;
  gp_Pnt             start;
  plain.rail.curve.D0(away == 0 ? plain.rail.first : plain.rail.last, start);
  TopoDS_Vertex a, z;
  TopExp::Vertices(edge, a, z);
  const TopoDS_Vertex keep =
    BRep_Tool::Pnt(a).Distance(start) < BRep_Tool::Pnt(z).Distance(start) ? a : z;
  double                      f, l;
  const auto                  curve = BRep_Tool::Curve(edge, f, l);
  GeomAPI_ProjectPointOnCurve on(BRep_Tool::Pnt(contactEnd(k, joint.side, end).vertex), curve, f, l);
  if (on.NbPoints() == 0 || on.LowerDistance() > rail.resolved(limit::together))
    fail(op, no_solution, "a contact misses the far edge it runs off over", {rail.edge});
  // A bounded notch can remove an interior interval, leaving two edge components.
  // Sequential cuts from the same original endpoint cannot represent that set.
  const bool  hereNotched = splitNotch(k, joint.side);
  const int   notchIndex  = hereNotched ? k : joint.piece;
  const int   notchSide   = hereNotched ? joint.side : joint.side ^ int(joint.swapped);
  const auto& notchPiece  = pieces[notchIndex];
  const auto& beginJoint  = notchPiece.joint[0];
  const auto& finish      = notchPiece.joint[1];
  int         crossings   = 0;
  for (const auto& p : pieces)
    crossings += p.joint[1].meets == Meets::Split && p.joint[1].crossed.IsSame(edge);
  const bool paired =
    refines(crossings == 2 && beginJoint.meets == Meets::Split && finish.meets == Meets::Split
            && beginJoint.side == notchSide && finish.side == notchSide
            && beginJoint.crossed.IsSame(edge) && finish.crossed.IsSame(edge)
            && pieces[beginJoint.piece].rail.s[notchSide ^ int(beginJoint.swapped)].trim.IsNull()
            && pieces[finish.piece].rail.s[notchSide ^ int(finish.swapped)].trim.IsNull()
            && !cuts.Contains(edge));
  if (!paired)
    cut(edge, keep.IsSame(a) ? z : a, contactEnd(k, joint.side, end).vertex, on.LowerDistanceParameter());
  else
  {
    const double crossing = on.LowerDistanceParameter();
    if (const auto found = seamCrossed.find(edge.TShape().get()); found == seamCrossed.end())
      seamCrossed[edge.TShape().get()] = {contactEnd(k, joint.side, end).vertex, crossing};
    else
    {
      const auto& [other, there] = found->second;
      const double low = std::min(crossing, there), high = std::max(crossing, there);
      if (!(low > f + tol::parameter && high < l - tol::parameter && high - low > tol::parameter))
        fail(op,
             unsupported_geometry,
             "a contact interval does not lie inside its crossed edge",
             {edge});
      const bool    lower = crossing < there;
      TopoDS_Vertex e0, e1;
      TopExp::Vertices(TopoDS::Edge(edge.Oriented(TopAbs_FORWARD)), e0, e1);
      TopoDS_Compound parts;
      builder.MakeCompound(parts);
      builder.Add(parts, trimmed(edge, e1, lower ? contactEnd(k, joint.side, end).vertex : other, low));
      builder.Add(parts, trimmed(edge, e0, lower ? other : contactEnd(k, joint.side, end).vertex, high));
      for (const auto& face : facesOfEdge.FindFromKey(edge))
        editOf(edits, face).replace.Bind(edge, parts);
    }
  }
}

void Attach::splitCrossingCut(int k, int end)
{
  const Piece& piece = pieces[k];
  const Rail&  rail  = piece.rail;
  const auto&  joint = piece.joint[end];
  // The edge the contact crosses onto the next face is cut back to it, from its
  // end at the blended edges.
  const TopoDS_Edge& edge = joint.crossed;
  TopoDS_Vertex      a, z;
  TopExp::Vertices(edge, a, z);
  const auto atRail = [&](const TopoDS_Vertex& v) {
    return std::any_of(pieces.begin(), pieces.end(), [&](const Piece& p) {
      for (TopExp_Explorer x(p.rail.edge, TopAbs_VERTEX); x.More(); x.Next())
        if (x.Current().IsSame(v))
          return true;
      return false;
    });
  };
  // Where both its ends are at blended edges (a hole's seam whose two rims are both
  // rounded), from the one at this edge.
  const auto atOwn = [&](const TopoDS_Vertex& v) {
    for (TopExp_Explorer x(rail.edge, TopAbs_VERTEX); x.More(); x.Next())
      if (x.Current().IsSame(v))
        return true;
    return false;
  };
  double                      f, l;
  const auto                  curve = BRep_Tool::Curve(edge, f, l);
  GeomAPI_ProjectPointOnCurve on(BRep_Tool::Pnt(contactEnd(k, joint.side, end).vertex), curve, f, l);
  if (on.NbPoints() == 0 || on.LowerDistance() > rail.resolved(limit::together))
    fail(op, no_solution, "a contact misses the edge it crosses onto the next face", {rail.edge});
  const double crossing = on.LowerDistanceParameter();
  // A sphere's seam the contact crosses twice (across the cap, from rim to rim) keeps
  // its part between: each crossing cuts it from the end it is nearer.
  int crossings = 0;
  for (const auto& p : pieces)
    crossings += p.joint[1].meets == Meets::Split && p.joint[1].crossed.IsSame(edge);
  // Two crossings retain the interval between them when the edge reaches
  // the blended chain. Its two ends need not both belong to that chain.
  double otherCrossing = crossing;
  if (refines(crossings == 2 && (atRail(a) || atRail(z))))
    for (int m = 0; m < n; ++m)
    {
      const auto& joint = pieces[m].joint[1];
      if (joint.meets != Meets::Split || !joint.crossed.IsSame(edge) || (m == k && end == 1)
          || (joint.piece == k && joint.end == end))
        continue;
      const auto&  q     = pieces[m];
      const gp_Pnt point = q.sheet.contact(joint.side)->Value(q.sheet.rail(q.uAt(1)));
      GeomAPI_ProjectPointOnCurve other(point, curve, f, l);
      if (other.NbPoints() && other.LowerDistance() <= rail.resolved(limit::together)
          && std::abs(other.LowerDistanceParameter() - crossing) > tol::parameter)
        otherCrossing = other.LowerDistanceParameter();
    }
  const TopoDS_Vertex from = otherCrossing != crossing ? (crossing < otherCrossing ? a : z)
                             : atRail(a) && atRail(z)
                               ? (crossings > 1           ? (crossing - f < l - crossing ? a : z)
                                  : atOwn(z) && !atOwn(a) ? z
                                                          : a)
                             : atRail(a) ? a
                             : atRail(z) ? z
                                         : TopoDS_Vertex();
  if (from.IsNull() && crossings == 2)
  {
    // An edge (or a seam) the contact crosses and comes back over, in the strip the
    // blend replaces between the crossings and on the faces past them: it is left in two.
    if (const auto found = seamCrossed.find(edge.TShape().get()); found == seamCrossed.end())
      seamCrossed[edge.TShape().get()] = {contactEnd(k, joint.side, end).vertex, crossing};
    else
    {
      const auto& [other, there] = found->second;
      const bool      lower      = crossing < there;
      TopoDS_Compound parts;
      builder.MakeCompound(parts);
      builder.Add(
        parts,
        trimmed(edge, z, lower ? contactEnd(k, joint.side, end).vertex : other, std::min(crossing, there)));
      builder.Add(
        parts,
        trimmed(edge, a, lower ? other : contactEnd(k, joint.side, end).vertex, std::max(crossing, there)));
      for (const auto& face : facesOfEdge.FindFromKey(edge))
        editOf(edits, face).replace.Bind(edge, parts);
    }
    return;
  }
  if (from.IsNull())
    fail(op,
         unsupported_geometry,
         "a contact crosses onto the next face over an edge that does not reach the "
         "blended edges, which the blend does not do yet",
         {rail.edge});
  cut(edge, from, contactEnd(k, joint.side, end).vertex, crossing);
}

void Attach::jointEdgeCuts(int k, int end, const TopoDS_Vertex& vertex)
{
  const Piece& piece = pieces[k];
  const Rail&  rail  = piece.rail;
  const auto&  joint = piece.joint[end];
  const Piece& other = pieces[joint.piece];
  // The other edges at the vertex lie between faces of the supports and end there:
  // each is cut back to the corner on its side.
  for (const auto& shape : edgesOfVertex.FindFromKey(vertex))
  {
    if (railEdge(shape) || crossed.Contains(shape))
      continue;
    const auto consumed = [&](const Rail& r) {
      return std::any_of(r.consumedSupports.begin(),
                         r.consumedSupports.end(),
                         [&](const TopoDS_Face& face) { return onFace(shape, face); });
    };
    if (consumed(rail) || consumed(other.rail))
    {
      swept.push_back(shape);
      continue;
    }
    // A narrow face a cliff or a notch replaces whole takes its edges with it.
    if (std::any_of(std::begin(rail.s), std::end(rail.s), [&](const Side& s) {
          return (s.cliffCurve || !s.trim.IsNull()) && s.cliffWhole && onFace(shape, s.face);
        }))
    {
      swept.push_back(shape);
      continue;
    }
    int side = -1;
    for (int s = 0; s < 2; ++s)
      if (onFace(shape, rail.s[s].face) || onFace(shape, other.rail.s[s ^ int(joint.swapped)].face))
        side = s;
    if (side < 0)
      fail(op,
           unsupported_geometry,
           "where a chain of blended edges meets another face, the corner is blended too, "
           "which the blend does not do yet",
           {rail.edge});
    if (poleInStrip(piece, side, TopoDS::Edge(shape)) && sweepToPole(TopoDS::Edge(shape), vertex))
      continue;
    double                      f, l;
    const auto                  curve = BRep_Tool::Curve(TopoDS::Edge(shape), f, l);
    GeomAPI_ProjectPointOnCurve on(BRep_Tool::Pnt(contactEnd(k, side, end).vertex), curve, f, l);
    if (on.NbPoints() == 0 || on.LowerDistance() > rail.resolved(limit::together))
    {
      if (sweepToPole(TopoDS::Edge(shape), vertex))
        continue;
      // A body pinched at the joint (a hole tangent to a support there, which meets it in
      // a line of four faces): no blend passes it.
      for (const auto& near : edgesOfVertex.FindFromKey(vertex))
        if (facesOfEdge.Contains(near) && facesOfEdge.FindFromKey(near).Extent() > 2)
          fail(op,
               invalid_input,
               "the body is pinched where the edge's blend would pass: more than two faces "
               "meet along an edge at its vertex (a hole tangent to a face beside it)",
               {rail.edge, near});
      fail(op,
           unsupported_geometry,
           "the contact of a chain's blend leaves the edge between its faces there, which "
           "the blend does not do yet",
           {rail.edge});
    }
    cut(shape, vertex, contactEnd(k, side, end).vertex, on.LowerDistanceParameter());
  }
}

void Attach::jointSection(int k, int end)
{
  const Piece& piece = pieces[k];
  const Rail&  rail  = piece.rail;
  const auto&  joint = piece.joint[end];
  const Piece& other = pieces[joint.piece];
  const Sheet& sheet = piece.sheet;
  const double here = piece.uAt(end), there = other.uAt(joint.end);
  // Where a closed blend's side is a notch, its section ends on the face beyond.
  double reach[2] = {sheet.across[0], sheet.across[1]};
  if (joint.piece == k)
    for (int side = 0; side < 2; ++side)
      if (!rail.s[side].trim.IsNull())
      {
        const gp_Pnt2d at = trimOf(k, side).onSheet->Value(sheet.rail(here));
        reach[side]       = sheet.swapped ? at.X() : at.Y();
      }
  const TopoDS_Edge x = edgeBetween(sheet.section(sheet.rail(here)),
                                    contactEnd(k, 0, end).vertex,
                                    reach[0],
                                    contactEnd(k, 1, end).vertex,
                                    reach[1],
                                    {},
                                    tol::length);
  double            lo, hi;
  BRep_Tool::Range(x, lo, hi);
  if (joint.piece == k || other.sheet.surface == sheet.surface)
  {
    // The seam, of a closed edge's sheet or of the one sheet a closed chain of split
    // parts is joined into: the forward use is where the face's boundary runs from
    // the first contact to the second, at the end of the rail when counter-clockwise.
    const bool ccw =
      sheet.swapped ? sheet.across[1] < sheet.across[0] : sheet.across[1] > sheet.across[0];
    const bool   forwardAtEnd = (sheet.across[0] < sheet.across[1]) == ccw;
    const double atHere = sheet.rail(here), atThere = other.sheet.rail(there);
    builder.UpdateEdge(x,
                       sheet.sectionUV(forwardAtEnd ? atThere : atHere),
                       sheet.sectionUV(forwardAtEnd ? atHere : atThere),
                       sheet.surface,
                       TopLoc_Location(),
                       tolerance);
  }
  else
  {
    builder.UpdateEdge(x,
                       sheet.sectionUV(sheet.rail(here)),
                       sheet.surface,
                       TopLoc_Location(),
                       tolerance);
    // On the other sheet the section runs between the same corners at the same
    // rate: the sections of a chain share their law and their parameter.
    const int atLo = sheet.across[0] < sheet.across[1] ? 0 : 1;
    builder.UpdateEdge(
      x,
      segment2d(
        other.sheet.uv(other.sheet.rail(there), other.sheet.across[atLo ^ int(joint.swapped)]),
        other.sheet.uv(other.sheet.rail(there),
                       other.sheet.across[(1 - atLo) ^ int(joint.swapped)]),
        lo,
        hi),
      other.sheet.surface,
      TopLoc_Location(),
      tolerance);
  }
  builder.Range(x, lo, hi);
  section[k][end] = section[joint.piece][joint.end] = x;
}

void Attach::jointEnd(int k, int end, const TopoDS_Vertex& endVertex)
{
  const Piece& piece = pieces[k];
  const Rail&  rail  = piece.rail;
  const auto&  joint = piece.joint[end];
  const Piece& other = pieces[joint.piece];
  const Sheet& sheet = piece.sheet;
  const double here = piece.uAt(end), there = other.uAt(joint.end);
  // One corner on each side, where both pieces' contacts end.
  for (int side = 0; side < 2; ++side)
  {
    const int to = side ^ int(joint.swapped);
    // Round a closed edge whose narrow face the blend uses up (a notch), the corner is on
    // the blend's edge on the face beyond.
    const bool   notched = joint.piece == k && !rail.s[side].trim.IsNull();
    const gp_Pnt point   = notched ? trimOf(k, side).curve->Value(sheet.rail(here))
                                   : sheet.contact(side)->Value(sheet.rail(here));
    const double gap =
      point.Distance(notched ? trimOf(k, side).curve->Value(sheet.rail(there))
                             : other.sheet.contact(to)->Value(other.sheet.rail(there)));
    if (gap > rail.resolved(limit::together))
      fail(op,
           unsupported_geometry,
           "the blends of tangent edges do not meet in one section there, which the blend "
           "does not do yet",
           {rail.edge});
    TopoDS_Vertex vertex;
    builder.MakeVertex(vertex, point, std::max(tolerance, gap));
    contactEnd(k, side, end).vertex = contactEnd(joint.piece, to, joint.end).vertex = vertex;
    contactEnd(k, side, end).parameter                                           = here;
    contactEnd(joint.piece, to, joint.end).parameter                             = there;
  }
  if (joint.meets == Meets::Split
      && (splitNotch(k, joint.side) || splitNotch(joint.piece, joint.side ^ int(joint.swapped))))
    splitNotchCut(k, end);
  else if (joint.meets == Meets::Split)
    splitCrossingCut(k, end);
  else
    jointEdgeCuts(k, end, endVertex);
  // The section between them: the sheet's own there, laid on the other sheet too,
  // or on its own sheet's other end where the piece closes on itself.
  // A split's parts are one face: no section between them.
  if (joint.meets == Meets::Split)
    return;
  jointSection(k, end);
}

void Attach::chainEnd(int k, int end)
{
  if (retainedBoundaries[k])
    return;
  const Piece&  piece = pieces[k];
  const Rail&   rail  = piece.rail;
  TopoDS_Vertex ends[2];
  TopExp::Vertices(TopoDS::Edge(rail.edge.Oriented(TopAbs_FORWARD)), ends[0], ends[1]);
  if (!piece.stopVertex[end].IsNull())
    ends[end] = piece.stopVertex[end];
  const auto& joint = piece.joint[end];
  if (joint.meets == Meets::Miter || joint.meets == Meets::Corner)
    return; // below
  if (joint.meets == Meets::Fold)
  {
    if (joint.piece > k || (joint.piece == k && joint.end > end))
      fold(k, end, joint.piece, joint.end);
    return;
  }
  if (joint.meets == Meets::Joint || joint.meets == Meets::Split)
  {
    if (joint.piece < k || (joint.piece == k && joint.end < end))
      return; // made from its other side
    jointEnd(k, end, ends[end]);
    return;
  }
  // A face stops the edge here.
  stop(k, end, ends[end]);
}

void Attach::chainEnds()
{
  for (int k = 0; k < n; ++k)
    stopPlaces(k);
  installBoundaries();
  for (int k = 0; k < n; ++k)
    for (int end = 0; end < 2; ++end)
      chainEnd(k, end);
}

TopoDS_Shape Attach::assemble()
{
  // The changed faces, and the solids with them and the blend faces; other shells as they were.
  std::vector<std::pair<TopoDS_Face, std::vector<TopoDS_Face>>> faces = retainedFaces;
  for (const auto& [oldFace, replacements] : retainedFaces)
    for (const auto& face : replacements)
      history->AddModified(oldFace, face);
  for (const auto& [original, replacement] : boundaryHistory)
    history->AddModified(original, replacement);
  for (const auto& edit : edits)
  {
    if (std::any_of(removed.begin(), removed.end(), [&](const TopoDS_Face& gone) {
          return gone.IsSame(edit.face);
        }))
      continue;
    EdgeMap replace = edit.replace;
    for (TopExp_Explorer e(edit.face, TopAbs_EDGE); e.More(); e.Next())
      if (!replace.IsBound(e.Current()) && cuts.Contains(e.Current()))
        replace.Bind(e.Current(), cuts.FindFromKey(e.Current()));
    const auto retained = std::find_if(faces.begin(), faces.end(), [&](const auto& pair) {
      return pair.first.IsSame(edit.face);
    });
    if (retained != faces.end())
    {
      if (!replace.IsEmpty() || !edit.bridges.empty())
      {
        std::vector<TopoDS_Face> updated;
        for (const auto& face : retained->second)
        {
          const auto parts = rebuiltFaces(face, replace, edit.bridges);
          updated.insert(updated.end(), parts.begin(), parts.end());
        }
        retained->second = std::move(updated);
      }
      continue;
    }
    faces.push_back({edit.face, rebuiltFaces(edit.face, replace, edit.bridges)});
    for (const auto& now : faces.back().second)
      history->AddModified(edit.face, now);
  }
  // What became of the rest: the cut edges are modified; each blended edge is gone, with the
  // vertices it ended at, and gives its blend, contacts and sections; a corner's vertex gives
  // its patch.
  for (int i = 1; i <= cuts.Extent(); ++i)
    history->AddModified(cuts.FindKey(i), cuts.FindFromIndex(i));
  for (const auto& shape : swept)
    if (!history->IsRemoved(shape))
      history->Remove(shape);
  for (const auto& face : removed)
    if (!history->IsRemoved(face))
      history->Remove(face);
  for (int k = 0; k < n; ++k)
  {
    // A round between faces that meet only at a vertex comes of that vertex.
    if (!pieces[k].apex.IsNull())
    {
      history->AddGenerated(pieces[k].apex, blends[k]);
      for (int side = 0; side < 2; ++side)
        history->AddGenerated(pieces[k].apex, contactEdge(k, side));
      continue;
    }
    const TopoDS_Shape& edge = pieces[k].rail.edge;
    if (!blends[k].IsNull())
      history->AddGenerated(edge, blends[k]);
    for (int side = 0; side < 2; ++side)
      for (const auto& made : sideEdges(k, side))
        if (!made.IsNull())
          history->AddGenerated(edge, made);
    for (int end = 0; end < 2; ++end)
      if (!section[k][end].IsNull() && pieces[k].joint[end].meets != Meets::Corner)
        history->AddGenerated(edge, section[k][end]);
    for (int end = 0; end < 2; ++end)
      if (pieces[k].joint[end].meets != Meets::Corner || parts[pieces[k].joint[end].piece].empty())
        for (const auto& part : chain[k][end])
          history->AddGenerated(edge, part);
    if (!history->IsRemoved(edge))
      history->Remove(edge);
    for (TopExp_Explorer v(edge, TopAbs_VERTEX); v.More(); v.Next())
      if (!history->IsRemoved(v.Current()))
        history->Remove(v.Current());
    for (const auto& over : pieces[k].swallowed)
      if (!blends[k].IsNull())
        history->AddGenerated(over, blends[k]);
  }
  for (const auto& f : folded)
  {
    if (!f.patch.IsNull())
      history->AddGenerated(pieces[f.a].rail.edge, f.patch);
    for (const auto& run : f.run)
      history->AddGenerated(pieces[f.a].rail.edge, run.edge);
  }
  for (const auto& roll : rolls)
  {
    history->AddGenerated(roll.vertex, roll.patch);
    history->AddGenerated(roll.vertex, roll.arc);
  }
  for (size_t c = 0; c < cornersAt.size(); ++c)
  {
    if (!patches[c].IsNull())
      history->AddGenerated(cornersAt[c].vertex, patches[c]);
    for (const auto& part : parts[c])
      history->AddGenerated(cornersAt[c].vertex, part);
    for (const auto& [k, end] : cornersAt[c].ends)
      if (parts[c].empty())
        history->AddGenerated(cornersAt[c].vertex, section[k][end]);
      else
        for (const auto& half : chain[k][end])
          history->AddGenerated(cornersAt[c].vertex, half);
  }
  std::vector<TopoDS_Shape> solids;
  for (TopExp_Explorer s(body, TopAbs_SOLID); s.More(); s.Next())
    solids.push_back(s.Current());
  for (const auto& solid : solids)
  {
    const auto holds = [&](const TopoDS_Shape& shape, const TopoDS_Shape& face) {
      for (TopExp_Explorer f(shape, TopAbs_FACE); f.More(); f.Next())
        if (f.Current().IsSame(face))
          return true;
      return false;
    };
    std::vector<int> own;
    for (int k = 0; k < n; ++k)
      if (holds(solid, pieces[k].rail.s[0].face))
        own.push_back(k);
    if (own.empty())
      continue;
    TopoDS_Solid out;
    builder.MakeSolid(out);
    for (TopoDS_Iterator shell(solid); shell.More(); shell.Next())
    {
      std::vector<int> here;
      for (int k : own)
        if (holds(shell.Value(), pieces[k].rail.s[0].face))
          here.push_back(k);
      if (here.empty())
      {
        builder.Add(out, shell.Value());
        continue;
      }
      TopoDS_Shell next;
      builder.MakeShell(next);
      for (TopExp_Explorer f(shell.Value(), TopAbs_FACE); f.More(); f.Next())
      {
        const TopoDS_Shape& face = f.Current();
        if (std::any_of(removed.begin(), removed.end(), [&](const TopoDS_Face& gone) {
              return gone.IsSame(face);
            }))
          continue;
        const auto changed = std::find_if(faces.begin(), faces.end(), [&](const auto& edit) {
          return edit.first.IsSame(face);
        });
        if (changed == faces.end())
          builder.Add(next, face);
        else
          for (const auto& now : changed->second)
            builder.Add(next, now);
      }
      NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> added;
      for (int k : here)
        if (!blends[k].IsNull() && added.Add(blends[k]))
          builder.Add(next, blends[k]);
      for (size_t c = 0; c < cornersAt.size(); ++c)
        if (std::find(here.begin(), here.end(), cornersAt[c].ends[0].first) != here.end())
        {
          if (!patches[c].IsNull())
            builder.Add(next, patches[c]);
          for (const auto& part : parts[c])
            builder.Add(next, part);
        }
      for (const auto& f : folded)
        if (!f.patch.IsNull() && std::find(here.begin(), here.end(), f.a) != here.end())
          builder.Add(next, f.patch);
      for (const auto& roll : rolls)
        if (std::find(here.begin(), here.end(), roll.a) != here.end())
          builder.Add(next, roll.patch);
      next.Closed(true);
      builder.Add(out, next);
    }
    body = withSolid(body, solid, out);
  }
  return body;
}

//! Common refinement of the final sheet, face ownership and contact-map domains.
//! Finite isolated crossings on the witnessed far edge bound uniformly owned fragments.
//! The surface is unchanged; Segment restricts copies of its fitted pcurves exactly.
//! Outer tails remain available to the stopping-face continuation. Overlap segments and
//! coincident crossings on distinct edges need a more general arrangement and are refused.
static void refineOwnership(std::vector<Piece>& pieces, std::vector<Corner>& corners)
{
  const int count = int(pieces.size());
  for (int k = 0; k < count; ++k)
  {
    checkCancelled(pieces[k].rail.op);
    const Piece original = pieces[k];
    const double low = original.uAt(0), high = original.uAt(1);
    if (!(high > low + tol::parameter))
      continue;
    struct Event { double u; int side; TopoDS_Edge edge; };
    std::vector<Event> events;
    for (int side = 0; side < 2; ++side)
    {
      const auto& owner = original.rail.s[side];
      if (owner.trim.IsNull() || owner.trimEdge.IsNull())
        continue;
      GeomAPI_IntCS crossing(original.sheet.contact(side), BRep_Tool::Surface(owner.trim));
      if (!crossing.IsDone() || crossing.NbSegments() != 0)
        fail(original.rail.op, unsupported_geometry, "contact ownership needs isolated boundary events", {original.rail.edge});
      double first, last;
      const auto curve = BRep_Tool::Curve(owner.trimEdge, first, last);
      for (int j = 1; j <= crossing.NbPoints(); ++j)
      {
        double pu, pv, w;
        crossing.Parameters(j, pu, pv, w);
        const double u = (w - original.sheet.alpha) / original.sheet.beta;
        if (!(u > low + tol::parameter && u < high - tol::parameter))
          continue;
        GeomAPI_ProjectPointOnCurve on(crossing.Point(j), curve, first, last);
        if (on.NbPoints() != 0 && on.LowerDistance() <= original.rail.resolved(limit::together))
          events.push_back({u, side, owner.trimEdge});
      }
    }
    if (events.empty())
      continue;
    std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) { return a.u < b.u; });
    for (size_t j = 1; j < events.size(); ++j)
      if (events[j].u - events[j-1].u <= tol::parameter
          && (events[j].side != events[j-1].side || !events[j].edge.IsSame(events[j-1].edge)))
        fail(original.rail.op, unsupported_geometry, "contact ownership has simultaneous boundary events", {original.rail.edge});
    events.erase(std::unique(events.begin(), events.end(), [](const Event& a, const Event& b) {
      return std::abs(a.u - b.u) <= tol::parameter;
    }), events.end());
    std::vector<double> bounds = {low};
    for (const auto& event : events)
      bounds.push_back(event.u);
    bounds.push_back(high);
    std::vector<int> indices = {k};
    for (size_t j = 0; j < events.size(); ++j)
    {
      indices.push_back(int(pieces.size()));
      pieces.push_back(original);
    }
    const auto outer = original.joint[1];
    for (size_t j = 0; j < indices.size(); ++j)
    {
      auto& part = pieces[indices[j]];
      part.rail.first = part.rail.t(bounds[j]);
      part.rail.last = part.rail.t(bounds[j+1]);
      for (int side = 0; side < 2; ++side)
      {
        if (!original.rail.s[side].trim.IsNull())
        {
          const auto point = original.sheet.contact(side)->Value(original.sheet.rail((bounds[j]+bounds[j+1])/2));
          BRepClass_FaceClassifier on(original.rail.s[side].face, point, original.rail.resolved(limit::together));
          if (on.State() == TopAbs_IN || on.State() == TopAbs_ON)
          {
            part.rail.s[side].trim.Nullify();
            part.rail.s[side].trimEdge.Nullify();
            part.rail.s[side].cliffWhole = false;
          }
        }
        const auto map = Handle(Geom2d_BSplineCurve)::DownCast(original.fitted.contact[side]->Copy());
        map->Segment(j == 0 ? map->FirstParameter() : bounds[j],
                     j+1 == indices.size() ? map->LastParameter() : bounds[j+1]);
        part.fitted.contact[side] = map;
      }
      if (j > 0)
        part.joint[0] = {Meets::Split, indices[j-1], 1, false, events[j-1].side, events[j-1].edge};
      if (j+1 < indices.size())
        part.joint[1] = {Meets::Split, indices[j+1], 0, false, events[j].side, events[j].edge};
    }
    if (outer.meets == Meets::Joint || outer.meets == Meets::Miter || outer.meets == Meets::Split)
    {
      pieces[outer.piece].joint[outer.end].piece = indices.back();
      pieces[outer.piece].joint[outer.end].end = 1;
    }
    if (outer.meets == Meets::Corner)
      for (auto& entry : corners[outer.piece].ends)
        if (entry.first == k && entry.second == 1)
          entry.first = indices.back();
  }
}

TopoDS_Shape attach(TopoDS_Shape                     body,
                    const std::vector<Piece>&        pieces,
                    const std::vector<Corner>&       cornersAt,
                    const Handle(BRepTools_History)& history,
                    std::vector<Join>* joins,
                    bool refineContacts,
                    StageReads* reads)
{
  std::vector<Piece> fragments;
  std::vector<Corner> fragmentCorners;
  // Only a contact on a face its support trims has an owner to change (refineOwnership).
  const bool owned = std::any_of(pieces.begin(), pieces.end(), [](const Piece& piece) {
    return std::any_of(std::begin(piece.rail.s), std::end(piece.rail.s), [](const Side& side) {
      return !side.trim.IsNull() && !side.trimEdge.IsNull();
    });
  });
  if (owned && !refineContacts && reads)
    reads->refines = true;
  if (refineContacts)
  {
    fragments = pieces;
    fragmentCorners = cornersAt;
    refineOwnership(fragments, fragmentCorners);
  }
  Attach a(std::move(body), refineContacts ? fragments : pieces,
           refineContacts ? fragmentCorners : cornersAt, history, refineContacts, reads);
  // Each stage may be cancelled before it starts (Patch also checks within its own).
  for (const auto stage : {&Attach::chainEnds,
                           &Attach::miters,
                           &Attach::cornerPatches,
                           &Attach::ridges,
                           &Attach::contacts,
                           &Attach::blendFaces,
                           &Attach::foldPatches,
                           &Attach::rollPatches,
                           &Attach::continuity,
                           &Attach::tolerances,
                           &Attach::collisions})
  {
    checkCancelled(a.op);
    (a.*stage)();
  }
  TopoDS_Shape result = a.assemble();
  if (joins)
  {
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> retained;
    TopExp::MapShapes(result, TopAbs_EDGE, retained);
    for (const auto& join : a.joins)
      if (retained.Contains(join.edge))
        joins->push_back(join);
  }
  return result;
}

} // namespace BRepExt::blending::detail
