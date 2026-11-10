// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Blend: The contacts, the blend faces, their continuity, tolerances and collisions.

#include <BRepBlendExt_Attach.hxx>
#include <BRepBlendExt_Separation.hxx>

namespace BRepExt::blending::detail
{

//! `edge` oriented to run as `along` does, by their directions where it starts (round a closed
//! edge its start is a vertex of its own).
static TopoDS_Edge sameWay(const TopoDS_Edge& edge, const TopoDS_Edge& along)
{
  BRepAdaptor_Curve own(TopoDS::Edge(edge.Oriented(TopAbs_FORWARD)));
  gp_Pnt            start;
  gp_Vec            direction;
  own.D1(own.FirstParameter(), start, direction);
  double       f, l;
  const auto   curve = BRep_Tool::Curve(along, f, l);
  const double w     = GeomAPI_ProjectPointOnCurve(start, curve, f, l).LowerDistanceParameter();
  gp_Pnt       on;
  gp_Vec       other;
  BRepAdaptor_Curve(TopoDS::Edge(along.Oriented(TopAbs_FORWARD))).D1(w, on, other);
  return direction.Dot(other) > 0 ? edge : TopoDS::Edge(edge.Reversed());
}

void Attach::ridges()
{
  for (int k = 0; k < n; ++k)
    for (int m = k + 1; m < n; ++m)
      for (int sk = 0; sk < 2; ++sk)
        for (int sm = 0; sm < 2; ++sm)
        {
          const Piece& one = pieces[k];
          const Piece& two = pieces[m];
          const Side&  a   = one.rail.s[sk];
          const Side&  b   = two.rail.s[sm];
          if (!a.face.IsSame(b.face) || a.cliffCurve || b.cliffCurve || !a.trim.IsNull()
              || !b.trim.IsNull())
            continue;
          if (!GeomLib_IsPlanarSurface(a.surface, tol::length).IsPlanar())
            continue;
          // The corners' parameters are solved to a tolerance: one may lie a little beyond the
          // fit's ends.
          const auto within = [](const Handle(Geom2d_BSplineCurve)& contact, double t) {
            return std::clamp(t, contact->FirstParameter(), contact->LastParameter());
          };
          const double ka = within(one.fitted.contact[sk], std::min(contactEnd(k, sk, 0).parameter, contactEnd(k, sk, 1).parameter)),
                       kb = within(one.fitted.contact[sk], std::max(contactEnd(k, sk, 0).parameter, contactEnd(k, sk, 1).parameter));
          const double ma = within(two.fitted.contact[sm], std::min(contactEnd(m, sm, 0).parameter, contactEnd(m, sm, 1).parameter)),
                       mb = within(two.fitted.contact[sm], std::max(contactEnd(m, sm, 0).parameter, contactEnd(m, sm, 1).parameter));
          // A blend with no length crosses nothing (its corners' balls are one).
          if (!(kb - ka > tol::parameterMerge) || !(mb - ma > tol::parameterMerge))
            continue;
          const Geom2dAPI_InterCurveCurve crossing(
            new Geom2d_TrimmedCurve(one.fitted.contact[sk], ka, kb),
            new Geom2d_TrimmedCurve(two.fitted.contact[sm], ma, mb),
            tol::parameter);
          // Contacts that run together are one (a face used up between the blends), not a neck.
          if (crossing.NbSegments() > 0)
            continue;
          // Not where the two contacts end together (a joint, a fold's crossing, a corner).
          const double near =
            std::max(one.rail.resolved(limit::together), two.rail.resolved(limit::together));
          std::vector<std::pair<double, double>> at;
          for (int i = 1; i <= crossing.NbPoints(); ++i)
          {
            const auto&  point   = crossing.Intersector().Point(i);
            const gp_Pnt there   = a.surface->Value(point.Value().X(), point.Value().Y());
            bool         corner_ = false;
            for (const auto& v :
                 {contactEnd(k, sk, 0).vertex, contactEnd(k, sk, 1).vertex, contactEnd(m, sm, 0).vertex, contactEnd(m, sm, 1).vertex})
              corner_ |= !v.IsNull() && BRep_Tool::Pnt(v).Distance(there) <= near;
            if (!corner_)
              at.emplace_back(point.ParamOnFirst(), point.ParamOnSecond());
          }
          if (at.empty())
            continue;
          // They part again past each crossing: the crossings pair up along either contact.
          // An odd count leaves a pair across another piece (the crease runs over a joint).
          if (at.size() % 2 != 0)
            fail(op,
                 unsupported_geometry,
                 "the blends of two parts of an outline cross over their cap across a joint of "
                 "their edges, which the blend does not do yet",
                 {one.rail.edge});
          std::sort(at.begin(), at.end());
          for (size_t i = 0; i < at.size(); i += 2)
          {
            Ridge r{k,
                    sk,
                    m,
                    sm,
                    {at[i].first, at[i + 1].first},
                    {at[i].second, at[i + 1].second},
                    {},
                    {}};
            for (int e = 0; e < 2; ++e)
            {
              const gp_Pnt2d uv = one.fitted.contact[sk]->Value(r.uk[e]);
              builder.MakeVertex(r.at[e], a.surface->Value(uv.X(), uv.Y()), tolerance);
            }
            const auto onOne = [&](double uk) {
              return one.sheet.uv(one.sheet.rail(uk), one.sheet.across[sk]);
            };
            const auto onTwo = [&](double um) {
              return two.sheet.uv(two.sheet.rail(um), two.sheet.across[sm]);
            };
            const auto edge = crossedMeeting(one.sheet,
                                             onOne(r.uk[0]),
                                             onOne(r.uk[1]),
                                             two.sheet,
                                             onTwo(r.um[0]),
                                             onTwo(r.um[1]),
                                             r.at[0],
                                             r.at[1]);
            if (!edge)
              fail(op,
                   no_solution,
                   "the crease where the blends of two parts of an outline cross cannot be made",
                   {one.rail.edge});
            r.edge = *edge;
            ridged.push_back(r);
          }
        }
}

void Attach::notchContact(int k, int side)
{
  const Piece& piece = pieces[k];
  const Sheet& sheet = piece.sheet;
  const Side&  notch = piece.rail.s[side];
  // A notch: the blend's edge on the face beyond, which takes it in the far edge's
  // place; the narrow face goes whole, with its far edge.
  const auto&  trim = trimOf(k, side);
  const double s0 = sheet.rail(contactEnd(k, side, 0).parameter), s1 = sheet.rail(contactEnd(k, side, 1).parameter);
  const auto   closed =
    std::find_if(rolled.begin(), rolled.end(), [&](const auto& r) { return r.first.first == k; });
  if (closed != rolled.end())
    contactEdge(k, side) = rolledEdge(k, side, closed->second);
  else
  {
    checkTrim(k, side, s0, s1);
    const auto [surface, location] = surfaceOf(notch.trim);
    const auto onBeyond            = GeomProjLib::Curve2d(trim.curve,
                                                          std::min(s0, s1),
                                                          std::max(s0, s1),
                                                          BRep_Tool::Surface(notch.trim));
    if (onBeyond.IsNull())
      fail(op,
           no_solution,
           "the blend's edge on the face beyond a narrow face cannot be laid on it",
           {piece.rail.edge});
    contactEdge(k, side) =
      edgeBetween(trim.curve,
                  contactEnd(k, side, 0).vertex,
                  s0,
                  contactEnd(k, side, 1).vertex,
                  s1,
                  {{trim.onSheet, sheet.surface, {}}, {onBeyond, surface, location}},
                  tolerance);
  }
  if (piece.rail.law == Law::Full || (piece.rail.closed && notch.cliffWhole))
  {
    // A full round's contact takes the far edge's place in the face beyond, running
    // as it ran (round a closed edge the contact is a loop of its own); so does a
    // closed blend's that uses up its narrow face (a chamfer rolling off a cylinder's
    // side onto its other cap).
    editOf(edits, notch.trim)
      .replace.Bind(notch.trimEdge, sameWay(contactEdge(k, side), notch.trimEdge));
    removed.push_back(notch.face);
    swept.push_back(notch.trimEdge);
    return;
  }
  Edit& edit = editOf(edits, notch.trim);
  edit.bridges.push_back(contactEdge(k, side));
  // From a split partway: the far edge is cut back to the crossing, the rest stays.
  if (splitNotch(k, side))
    return;
  edit.replace.Bind(notch.trimEdge, TopoDS_Edge());
  const auto F0 = trimFar.find({k, side, 0}), F1 = trimFar.find({k, side, 1});
  if (notch.cliffWhole)
  {
    removed.push_back(notch.face);
    swept.push_back(notch.trimEdge);
    return;
  }
  // The narrow face goes on past an end: it loses the blended edge and the strip, and
  // the far edge keeps its parts past the ends where the face goes on, in both faces.
  Edit& onNarrow = editOf(edits, notch.face);
  onNarrow.replace.Bind(notch.trimEdge, TopoDS_Edge());
  onNarrow.replace.Bind(piece.rail.edge, TopoDS_Edge());
  const TopoDS_Edge forward = TopoDS::Edge(notch.trimEdge.Oriented(TopAbs_FORWARD));
  TopoDS_Vertex     v0, v1;
  TopExp::Vertices(forward, v0, v1);
  double     ff, fl;
  const auto farCurve = BRep_Tool::Curve(forward, ff, fl);
  for (const auto& at : {F0, F1})
  {
    if (at == trimFar.end())
      continue;
    const double w = GeomAPI_ProjectPointOnCurve(BRep_Tool::Pnt(at->second), farCurve, ff, fl)
                       .LowerDistanceParameter();
    // The part on the side of the vertex the blend does not reach.
    const double w0 =
      GeomAPI_ProjectPointOnCurve(BRep_Tool::Pnt(contactEnd(k, side, 0).vertex), farCurve, ff, fl)
        .LowerDistanceParameter();
    const double w1 =
      GeomAPI_ProjectPointOnCurve(BRep_Tool::Pnt(contactEnd(k, side, 1).vertex), farCurve, ff, fl)
        .LowerDistanceParameter();
    const double      inside = (w0 + w1) / 2;
    const TopoDS_Edge part =
      inside > w ? trimmed(forward, v1, at->second, w) : trimmed(forward, v0, at->second, w);
    onNarrow.bridges.push_back(part);
    edit.bridges.push_back(part);
  }
  swept.push_back(notch.trimEdge);
}

void Attach::cliffContact(int k, int side)
{
  const Piece& piece = pieces[k];
  const Sheet& sheet = piece.sheet;
  const Side&  cliff = piece.rail.s[side];
  // A cliff: its far edge, now between the blend and the face beyond, which keeps
  // it; the narrow face goes whole.
  const TopoDS_Face beyond       = beyondCliff(cliff);
  const auto [surface, location] = surfaceOf(beyond);
  const double s0 = sheet.rail(contactEnd(k, side, 0).parameter), s1 = sheet.rail(contactEnd(k, side, 1).parameter);
  const auto   onBeyond = GeomProjLib::Curve2d(sheet.contact(side),
                                               std::min(s0, s1),
                                               std::max(s0, s1),
                                               BRep_Tool::Surface(beyond));
  if (beyond.IsNull() || onBeyond.IsNull())
    fail(op,
         no_solution,
         "the blend's contact on the face's far edge cannot be laid on the face beyond",
         {piece.rail.edge});
  contactEdge(k, side) =
    edgeBetween(sheet.contact(side),
                contactEnd(k, side, 0).vertex,
                s0,
                contactEnd(k, side, 1).vertex,
                s1,
                {{sheet.contactUV(side), sheet.surface, {}}, {onBeyond, surface, location}},
                tolerance);
  if (!cliff.cliffWhole)
  {
    // The narrow face goes on past the blend: the far edge keeps its parts outside
    // the blend's in both faces, the blended edge goes from the narrow one, and the
    // face beyond takes the contact between those parts.
    double f, l;
    BRep_Tool::Range(cliff.cliff, f, l);
    const auto at = [&](const TopoDS_Vertex& v) {
      return GeomAPI_ProjectPointOnCurve(BRep_Tool::Pnt(v), cliff.cliffCurve, f, l)
        .LowerDistanceParameter();
    };
    const double         w0 = at(contactEnd(k, side, 0).vertex), w1 = at(contactEnd(k, side, 1).vertex);
    const TopoDS_Vertex& low  = w0 < w1 ? contactEnd(k, side, 0).vertex : contactEnd(k, side, 1).vertex;
    const TopoDS_Vertex& high = w0 < w1 ? contactEnd(k, side, 1).vertex : contactEnd(k, side, 0).vertex;
    TopoDS_Vertex        v0, v1;
    TopExp::Vertices(TopoDS::Edge(cliff.cliff.Oriented(TopAbs_FORWARD)), v0, v1);
    Edit& onNarrow = editOf(edits, cliff.face);
    Edit& onBeyond = editOf(edits, beyond);
    onNarrow.replace.Bind(cliff.cliff, TopoDS_Edge());
    onBeyond.replace.Bind(cliff.cliff, TopoDS_Edge());
    onNarrow.replace.Bind(piece.rail.edge, TopoDS_Edge());
    const TopoDS_Edge forward = TopoDS::Edge(cliff.cliff.Oriented(TopAbs_FORWARD));
    if (!low.IsSame(v0))
    {
      const TopoDS_Edge part = trimmed(forward, v1, low, std::min(w0, w1));
      onNarrow.bridges.push_back(part);
      onBeyond.bridges.push_back(part);
    }
    if (!high.IsSame(v1))
    {
      const TopoDS_Edge part = trimmed(forward, v0, high, std::max(w0, w1));
      onNarrow.bridges.push_back(part);
      onBeyond.bridges.push_back(part);
    }
    onBeyond.bridges.push_back(contactEdge(k, side));
    return;
  }
  editOf(edits, beyond).replace.Bind(cliff.cliff, sameWay(contactEdge(k, side), cliff.cliff));
  removed.push_back(cliff.face);
}

void Attach::passedContact(int k, int side)
{
  const Piece& piece = pieces[k];
  const Sheet& sheet = piece.sheet;
  if (piece.joint[0].meets == Meets::Split || piece.joint[1].meets == Meets::Split
      || std::any_of(ridged.begin(), ridged.end(), [&](const Ridge& r) {
           return (r.k == k && r.sk == side) || (r.m == k && r.sm == side);
         }))
    fail(op,
         unsupported_geometry,
         "a contact crosses onto another face of its support beside a split or a ridge, "
         "which the blend does not do yet",
         {piece.rail.edge});
  const Handle(Geom2d_Curve) onSupport = alongRail(piece.fitted.contact[side], sheet);

  struct Place
  {
    double        u;
    TopoDS_Vertex at;
  };

  // The places in order from the start's corner to the end's, and the face after each.
  std::vector<Place>                            places = {{contactEnd(k, side, 0).parameter, contactEnd(k, side, 0).vertex}};
  std::vector<std::pair<TopoDS_Face, gp_Vec2d>> on;
  const auto&                                   start = passages[{k, side, 0}];
  for (auto p = start.rbegin(); p != start.rend(); ++p)
  {
    on.push_back({p->face, p->shift});
    places.push_back({p->u, p->at});
  }
  on.push_back({piece.rail.s[side].face, gp_Vec2d()});
  for (const auto& p : passages[{k, side, 1}])
  {
    places.push_back({p.u, p.at});
    on.push_back({p.face, p.shift});
  }
  places.push_back({contactEnd(k, side, 1).parameter, contactEnd(k, side, 1).vertex});
  std::vector<TopoDS_Edge>& along = runs[{k, side}];
  for (size_t i = 0; i < on.size(); ++i)
  {
    const double s0 = sheet.rail(places[i].u), s1 = sheet.rail(places[i + 1].u);
    const auto& [face, shift] = on[i];
    const Handle(Geom2d_Curve) laid =
      shift.SquareMagnitude() == 0 ? onSupport
                                   : Handle(Geom2d_Curve)::DownCast(onSupport->Translated(-shift));
    const bool own              = face.IsSame(piece.rail.s[side].face);
    const auto [surface, place] = surfaceOf(face);
    if (own && piece.rail.parameterSupports && !contactWithinFace(face, laid, s0, s1))
      fail(op,
           unsupported_geometry,
           "a retained blend contact lies outside its finite support",
           {piece.rail.edge});
    const TopoDS_Edge part =
      edgeBetween(sheet.contact(side),
                  places[i].at,
                  s0,
                  places[i + 1].at,
                  s1,
                  {{sheet.contactUV(side), sheet.surface, {}}, {laid, surface, place}},
                  tolerance);
    Edit& edit = editOf(edits, face);
    if (own && onFace(piece.rail.edge, edit.face) && !edit.replace.IsBound(piece.rail.edge))
      edit.replace.Bind(piece.rail.edge, part);
    else
      edit.bridges.push_back(part);
    along.push_back(part);
    if (own)
      contactEdge(k, side) = part;
    else
      partFace.Bind(part, face);
  }
}

bool Attach::ridgedContact(int                         k,
                           int                         side,
                           const Handle(Geom_Surface)& support,
                           const TopLoc_Location&      location)
{
  const Piece&                                  piece = pieces[k];
  const Sheet&                                  sheet = piece.sheet;
  std::vector<std::pair<double, TopoDS_Vertex>> stops;
  std::vector<std::pair<double, TopoDS_Edge>>   creases;
  for (const auto& r : ridged)
    for (const auto& [p, sp, us] : {std::tuple{r.k, r.sk, r.uk}, std::tuple{r.m, r.sm, r.um}})
      if (p == k && sp == side)
      {
        creases.emplace_back(std::min(us[0], us[1]), r.edge);
        stops.emplace_back(us[0], r.at[0]);
        stops.emplace_back(us[1], r.at[1]);
      }
  if (stops.empty())
    return false;
  if (piece.joint[0].meets == Meets::Split || piece.joint[1].meets == Meets::Split)
    fail(op,
         unsupported_geometry,
         "the blends of two parts of an outline cross over their cap beside a split, which "
         "the blend does not do yet",
         {piece.rail.edge});
  const bool rising = contactEnd(k, side, 0).parameter < contactEnd(k, side, 1).parameter;
  const auto order  = [&](const auto& x, const auto& y) {
    return rising ? x.first < y.first : x.first > y.first;
  };
  std::sort(stops.begin(), stops.end(), order);
  std::sort(creases.begin(), creases.end(), order);
  stops.insert(stops.begin(), {contactEnd(k, side, 0).parameter, contactEnd(k, side, 0).vertex});
  stops.push_back({contactEnd(k, side, 1).parameter, contactEnd(k, side, 1).vertex});
  // Every part of the edge's contact on the face is a bridge then: the parts of the
  // cap they bound are islands of it.
  Edit& edit = editOf(edits, piece.rail.s[side].face);
  if (edit.replace.IsBound(piece.rail.edge) && !edit.replace.Find(piece.rail.edge).IsNull())
    edit.bridges.push_back(TopoDS::Edge(edit.replace.Find(piece.rail.edge)));
  edit.replace.Bind(piece.rail.edge, TopoDS_Edge());
  std::vector<TopoDS_Edge>& along = runs[{k, side}];
  for (size_t i = 0; i + 1 < stops.size(); i += 2)
  {
    const double s0 = sheet.rail(stops[i].first), s1 = sheet.rail(stops[i + 1].first);
    if (piece.rail.parameterSupports
        && !contactWithinFace(piece.rail.s[side].face,
                              alongRail(piece.fitted.contact[side], sheet),
                              s0,
                              s1))
      fail(op,
           unsupported_geometry,
           "a retained blend contact lies outside its finite support",
           {piece.rail.edge});
    const TopoDS_Edge part =
      edgeBetween(sheet.contact(side),
                  stops[i].second,
                  s0,
                  stops[i + 1].second,
                  s1,
                  {{sheet.contactUV(side), sheet.surface, {}},
                   {alongRail(piece.fitted.contact[side], sheet), support, location}},
                  tolerance);
    edit.bridges.push_back(part);
    along.push_back(part);
    if (i / 2 < creases.size())
      along.push_back(creases[i / 2].second);
  }
  contactEdge(k, side) = along.front();
  return true;
}

void Attach::plainContact(int k, int side,
                          const Handle(Geom_Surface)& support,
                          const TopLoc_Location& location)
{
  for (auto& interval : contactIntervals[k][side])
    plainInterval(k, side, interval, support, location);
}

void Attach::plainInterval(int k, int side, ContactInterval& interval,
                           const Handle(Geom_Surface)& support,
                           const TopLoc_Location& location)
{
  const Piece& piece = pieces[k];
  const Sheet& sheet = piece.sheet;
  // The parts of a split edge that stay on this face share one contact, as their face
  // is one: from the first part's start to the last one on the face, its curve on the
  // face the parts' curves joined.
  // A seam of the face, which it crosses on itself, ends its edge there all the same.
  std::vector<int> run = {k};
  while (
    pieces[run.back()].joint[1].meets == Meets::Split
    && pieces[pieces[run.back()].joint[1].piece].rail.s[side].face.IsSame(piece.rail.s[side].face)
    && !BRep_Tool::IsClosed(pieces[run.back()].joint[1].crossed, piece.rail.s[side].face)
    && pieces[pieces[run.back()].joint[1].piece].rail.s[side].trim.IsNull())
    run.push_back(pieces[run.back()].joint[1].piece);
  const int    last = run.back();
  TopoDS_Edge& edge = interval.edge;
  // On an exact sheet the contact is a line or a circle: laid on the support exactly where
  // OpenCascade can, else the fit's feet.
  Handle(Geom2d_Curve) onSupport = alongRail(piece.fitted.contact[side], sheet);
  if (run.size() > 1)
  {
    std::vector<Handle(Geom2d_BSplineCurve)> parts;
    for (int m : run)
      parts.push_back(alongRail(pieces[m].fitted.contact[side], pieces[m].sheet));
    onSupport = endToEnd(parts, piece.rail.resolved(limit::together));
    if (onSupport.IsNull())
      fail(op,
           no_solution,
           "the contacts of a split blend's parts do not join on their face",
           {piece.rail.edge});
  }
  if (piece.exact())
  {
    const double s0 = sheet.rail(interval.ends[0].parameter), s1 = sheet.rail(interval.ends[1].parameter);
    const auto   exact =
      GeomProjLib::Curve2d(sheet.contact(side), s0, s1, piece.rail.s[side].surface);
    if (!exact.IsNull())
      onSupport = aligned(exact, s0, onSupport->Value(s0), piece.rail.s[side].surface);
  }
  const double s0 = sheet.rail(interval.ends[0].parameter), s1 = pieces[last].sheet.rail((last == k ? interval.ends[1] : contactEnd(last, side, 1)).parameter);
  if (piece.rail.parameterSupports
      && !contactWithinFace(piece.rail.s[side].face, onSupport, s0, s1))
    fail(op,
         unsupported_geometry,
         "a retained blend contact lies outside its finite support",
         {piece.rail.edge});
  // Where another blend's contact on this face is this one, the face between them is
  // used up (a box rounded to half its width): the two blends meet along it.
  if (const auto same = sharedContact(piece.rail.s[side].face,
                                      interval.ends[0].vertex,
                                      (last == k ? interval.ends[1] : contactEnd(last, side, 1)).vertex,
                                      sheet.contact(side)->Value((s0 + s1) / 2),
                                      piece.rail.scale);
      !same.IsNull())
  {
    double     f, l;
    const auto curve   = BRep_Tool::Curve(same, f, l);
    const auto onSheet = laidOn(curve, f, l, sheet.surface);
    if (onSheet.IsNull())
      fail(op,
           no_solution,
           "a contact two blends share cannot be laid on the second",
           {piece.rail.edge});
    builder.UpdateEdge(same, onSheet, sheet.surface, TopLoc_Location(), tolerance);
    builder.Range(same, f, l);
    interval.edge = same;
    for (int m : run)
      if (m != k)
        contactEdge(m, side) = same;
    usedUp.Add(piece.rail.s[side].face);
    return;
  }
  edge = edgeBetween(sheet.contact(side),
                     interval.ends[0].vertex,
                     s0,
                     (last == k ? interval.ends[1] : contactEnd(last, side, 1)).vertex,
                     s1,
                     {{sheet.contactUV(side), sheet.surface, {}}, {onSupport, support, location}},
                     tolerance);
  for (int m : run)
    if (m != k)
      contactEdge(m, side) = edge;
  // In a face beside the edge it takes the edge's place; in a face the contact crossed onto,
  // or beside a part of the edge already taken, it closes the gap it leaves.
  Edit& edit = editOf(edits, piece.rail.s[side].face);
  if (onFace(piece.rail.edge, edit.face) && !edit.replace.IsBound(piece.rail.edge))
    edit.replace.Bind(piece.rail.edge, edge);
  else
    edit.bridges.push_back(edge);
}

void Attach::sideContact(int k, int side)
{
  if (!contactEdge(k, side).IsNull())
    return; // one edge for parts of a split already made
  if (vanished(k))
    return;
  const Piece& piece = pieces[k];
  // A part that a split's earlier part runs on into on this face is made with that one's.
  if (const auto& before = piece.joint[0];
      before.meets == Meets::Split
      && pieces[before.piece].rail.s[side].face.IsSame(piece.rail.s[side].face)
      && !BRep_Tool::IsClosed(before.crossed, piece.rail.s[side].face)
      && piece.rail.s[side].trim.IsNull() && pieces[before.piece].rail.s[side].trim.IsNull())
    return;
  if (!piece.rail.s[side].trim.IsNull())
  {
    // A notch: the blend's edge on the face beyond, which takes it in the far edge's
    // place; the narrow face goes whole, with its far edge.
    notchContact(k, side);
    return;
  }
  if (piece.rail.s[side].cliffCurve)
  {
    // A cliff: its far edge, now between the blend and the face beyond, which keeps
    // it; the narrow face goes whole.
    cliffContact(k, side);
    return;
  }
  const auto [support, location] = surfaceOf(piece.rail.s[side].face);
  // A contact that crosses onto other faces of its support's surface at an end: one edge
  // on each face, from corner to corner, as one run along the side.
  if (passages.count({k, side, 0}) || passages.count({k, side, 1}))
  {
    passedContact(k, side);
    return;
  }
  // Where ridges take the contact's place: its parts between them, each an edge of the
  // support, which bridge the cap's gaps (two islands of it where a neck closes), and the
  // ridges, as one run along the side.
  if (ridgedContact(k, side, support, location))
    return;
  plainContact(k, side, support, location);
}

void Attach::contacts()
{
  for (const auto& piece : pieces)
    for (const auto& face : piece.rail.consumedSupports)
    {
      removed.push_back(face);
      for (TopExp_Explorer e(face, TopAbs_EDGE); e.More(); e.Next())
      {
        if (railEdge(e.Current()))
          continue;
        for (const auto& neighbor : facesOfEdge.FindFromKey(e.Current()))
          if (!neighbor.IsSame(face))
            editOf(edits, neighbor).replace.Bind(e.Current(), TopoDS_Edge());
        swept.push_back(e.Current());
      }
    }
  // Edges a blend runs on over (Blend::absorbEnds) go.
  for (const auto& piece : pieces)
    for (const auto& edge : piece.swallowed)
    {
      for (const auto& face : facesOfEdge.FindFromKey(edge))
        editOf(edits, face).replace.Bind(edge, TopoDS_Edge());
      swept.push_back(edge);
    }
  // The contacts: on each sheet, the lines at `across`; on the supports, the fit's feet. Each
  // takes the edge's place in its support.
  for (int k = 0; k < n; ++k)
    for (int side = 0; side < 2; ++side)
      sideContact(k, side);
  // A face used up between the blends goes: each of its edges is a blend's that has no length,
  // or whose contact on it another blend shares.
  for (int k = 0; k < n; ++k)
    if (vanished(k))
      for (const auto& side : pieces[k].rail.s)
        usedUp.Add(side.face);
  for (int i = 1; i <= usedUp.Extent(); ++i)
  {
    const TopoDS_Face& face = TopoDS::Face(usedUp.FindKey(i));
    for (TopExp_Explorer e(face, TopAbs_EDGE); e.More(); e.Next())
    {
      bool gone = false;
      for (int k = 0; k < n && !gone; ++k)
      {
        if (!pieces[k].rail.edge.IsSame(e.Current()))
          continue;
        const int side  = pieces[k].rail.s[0].face.IsSame(face) ? 0 : 1;
        int       users = 0;
        for (int m = 0; m < n; ++m)
          for (int t = 0; t < 2; ++t)
            users += !contactEdge(k, side).IsNull() && contactEdge(m, t).IsSame(contactEdge(k, side));
        gone = vanished(k) || users == 2;
      }
      if (!gone)
        fail(op,
             unsupported_geometry,
             "the blends use up part of a face between them, which the blend does not make yet",
             {pieces.front().rail.edge});
    }
    removed.push_back(face);
  }
}

BoundaryGraph Attach::boundaryGraph(const std::vector<int>& group) const
{
  const int  k = group.front(), last = group.back();
  if (retainedBoundaries[k])
    return *retainedBoundaries[k];
  const auto from = [](const TopoDS_Edge& edge, const TopoDS_Vertex& start) {
    return TopExp::FirstVertex(edge).IsSame(start) ? edge : TopoDS::Edge(edge.Reversed());
  };
  // Along the first contact, across the far end, back along the second, across the near end.
  std::vector<TopoDS_Edge> loop;
  const auto               once = [&](const TopoDS_Edge& edge) {
    if (loop.empty() || !loop.back().IsSame(edge))
      loop.push_back(edge);
  };
  // An end of several edges runs from side 0's corner to side 1's: forward at the far end,
  // back at the near one.
  const auto across = [&](int m, int end, const TopoDS_Vertex& start) {
    // A blend that closes on the face beyond its narrow face has no section there: its
    // contacts meet at a corner (Attach::consumed).
    if (section[m][end].IsNull() && chain[m][end].empty())
      return;
    // A degenerate end (Rail::pinch) runs across the sheet as the end does.
    if (BRep_Tool::Degenerated(section[m][end]))
    {
      const bool rising = pieces[m].sheet.across[0] < pieces[m].sheet.across[1];
      loop.push_back((end == 1) == rising ? section[m][end]
                                          : TopoDS::Edge(section[m][end].Reversed()));
      return;
    }
    if (chain[m][end].empty())
    {
      loop.push_back(from(section[m][end], start));
      return;
    }
    TopoDS_Vertex at    = start;
    const auto&   edges = chain[m][end];
    for (size_t i = 0; i < edges.size(); ++i)
    {
      const TopoDS_Edge& edge = edges[end == 1 ? i : edges.size() - 1 - i];
      loop.push_back(from(edge, at));
      at = TopExp::LastVertex(loop.back(), true);
    }
  };
  // A side's edges in turn, from the corner it starts at.
  const auto along = [&](int m, int side, const TopoDS_Vertex& start) {
    const auto edges = sideEdges(m, side);
    if (edges.size() == 1)
    {
      once(side == 0 ? from(edges.front(), start) : TopoDS::Edge(edges.front().Reversed()));
      return;
    }
    TopoDS_Vertex at = start;
    for (size_t i = 0; i < edges.size(); ++i)
    {
      once(from(edges[side == 0 ? i : edges.size() - 1 - i], at));
      at = TopExp::LastVertex(loop.back(), true);
    }
  };
  for (int m : group)
    along(m, 0, contactEnd(m, 0, 0).vertex);
  across(last, 1, contactEnd(last, 0, 1).vertex);
  for (auto m = group.rbegin(); m != group.rend(); ++m)
    along(*m, 1, contactEnd(*m, 1, 1).vertex);
  // A seam's first vertex is that of both its uses, so the closed edge's are set apart.
  const size_t far = std::find_if(loop.begin(),
                                  loop.end(),
                                  [&](const TopoDS_Edge& e) { return e.IsSame(section[last][1]); })
                     - loop.begin();
  if (section[k][0].IsSame(section[last][1]) && chain[k][0].empty())
    loop.push_back(TopoDS::Edge(loop[far].Reversed()));
  else
    across(k, 0, contactEnd(k, 1, 0).vertex);
  // The existing end constructors supply a single cycle in occurrence order.
  // Keep their ports separate even when a seam or pinch identifies vertices.
  BoundaryGraph graph;
  for (size_t i = 0; i < loop.size(); ++i)
    graph.arcs.push_back({loop[i], i, (i + 1) % loop.size()});
  return graph;
}

void Attach::blendFaces()
{
  // Faces bounded by the retained cycles, facing out as the supports do.
  for (int k = 0; k < n; ++k)
  {
    // One face for the parts of an edge that a split joins, made at the first of them; none
    // for a blend that has no length.
    if (pieces[k].joint[0].meets == Meets::Split || vanished(k))
      continue;
    std::vector<int> group = {k};
    while (pieces[group.back()].joint[1].meets == Meets::Split)
      group.push_back(pieces[group.back()].joint[1].piece);
    const Piece& piece = pieces[k];
    const Sheet& sheet = piece.sheet;
    builder.MakeFace(blends[k], sheet.surface, tolerance);
    const bool ccw =
      sheet.swapped ? sheet.across[1] < sheet.across[0] : sheet.across[1] > sheet.across[0];
    const auto loops = boundaryGraph(group).loops();
    if (!loops)
      fail(op,
           no_solution,
           "the retained blend boundary has open or ambiguous incidence",
           {piece.rail.edge});
    for (const auto& loop : *loops)
    {
      TopoDS_Wire wire;
      builder.MakeWire(wire);
      for (const auto& edge : loop)
        builder.Add(wire, ccw ? edge : TopoDS::Edge(edge.Reversed()));
      wire.Closed(true);
      builder.Add(blends[k], wire);
    }
    gp_Pnt         p;
    gp_Vec         du, dv, outward;
    const double   middle = piece.uAt(0) / 2 + piece.uAt(1) / 2;
    const gp_Pnt2d at     = sheet.uv(sheet.rail(middle), (sheet.across[0] + sheet.across[1]) / 2);
    sheet.surface->D1(at.X(), at.Y(), p, du, dv);
    for (int side = 0; side < 2; ++side)
    {
      const gp_Pnt2d foot = piece.fitted.contact[side]->Value(middle);
      Jet            j;
      if (jetAt(piece.rail.s[side], foot.X(), foot.Y(), j))
        outward += j.n;
    }
    // A full round's supports may face apart: the face between, which its middle touches.
    if (piece.rail.law == Law::Full)
    {
      GeomAPI_ProjectPointOnSurf on(p, piece.rail.middle.surface);
      double                     mu, mv;
      Jet                        j;
      if (on.NbPoints() == 0)
        fail(op,
             no_solution,
             "a full round does not touch the face between its edges",
             {piece.rail.edge});
      on.LowerDistanceParameters(mu, mv);
      if (jetAt(piece.rail.middle, mu, mv, j))
        outward = j.n;
    }
    if (du.Crossed(dv).Dot(outward) < 0)
      blends[k].Reverse();
    for (int m : group)
      blends[m] = blends[k];
  }
}

void Attach::continuity()
{
  // The continuity of each new edge between its faces, as the blends are made: what OpenCascade's
  // own fillets record, and what algorithms that take tangent neighbours read. Contacts are G2
  // for the G2 law and G1 for conics; joints, seams and corner sides likewise where the law
  // gives it. Mitres and ends on stopping faces are sharp and carry none.
  // Across a joint where a support's second form jumps, the ball's path and so the sheet's
  // curvature jump with it: that joint is only G1 (the limit is the one blending::accepted holds
  // the rest to).
  const auto supportBreaks = [&](int k, int end) {
    const Piece& a     = pieces[k];
    const auto&  joint = a.joint[end];
    const Piece& b     = pieces[joint.piece];
    for (int side = 0; side < 2; ++side)
    {
      const int other = joint.swapped ? 1 - side : side;
      if (a.rail.s[side].cliffCurve || b.rail.s[other].cliffCurve)
        continue;
      const gp_Pnt2d ua = a.fitted.contact[side]->Value(end),
                     ub = b.fitted.contact[other]->Value(joint.end);
      Jet ja, jb;
      if (!jetAt(a.rail.s[side], ua.X(), ua.Y(), ja) || !jetAt(b.rail.s[other], ub.X(), ub.Y(), jb))
        return true;
      const SurfaceJet here{op, ja.p, ja.du, ja.dv, ja.duu, ja.dvv, ja.duv, ja.n},
        there{op, jb.p, jb.du, jb.dv, jb.duu, jb.dvv, jb.duv, jb.n};
      if (bendDifference(here, there, ja.du.Normalized()).first * a.rail.scale > limit::seamBend)
        return true;
    }
    return false;
  };
  const auto record = [&](const TopoDS_Shape& edge, JoinRole role, int order) {
    if (!edge.IsNull())
      joins.push_back({edge, role, order});
  };
  const auto mark = [&](const TopoDS_Edge& edge, const TopoDS_Face& a,
                        const TopoDS_Face& b, JoinRole role, GeomAbs_Shape smooth) {
    record(edge, role, smooth == GeomAbs_G2 ? 2 : 1);
    builder.Continuity(edge, a, b, smooth);
  };
  for (int k = 0; k < n; ++k)
  {
    const Piece& piece = pieces[k];
    const bool bevel = piece.rail.law == Law::Chamfer;
    const GeomAbs_Shape smooth = piece.rail.law == Law::G2 ? GeomAbs_G2 : GeomAbs_G1;
    if (vanished(k))
    {
      if (!bevel)
        mark(section[k][0], patches[piece.joint[0].piece],
             patches[piece.joint[1].piece], JoinRole::Corner, smooth);
      else
        record(section[k][0], JoinRole::Corner, 0);
      continue;
    }
    for (int side = 0; side < 2; ++side)
    {
      if (usedUp.Contains(piece.rail.s[side].face))
      {
        for (int m = 0; m < n; ++m)
          for (int t = 0; t < 2; ++t)
            if (m != k && contactEdge(m, t).IsSame(contactEdge(k, side)))
            {
              if (!bevel)
                mark(contactEdge(k, side), blends[k], blends[m], JoinRole::UsedSupport, smooth);
              else
                record(contactEdge(k, side), JoinRole::UsedSupport, 0);
            }
        continue;
      }
      for (const auto& edge : sideEdges(k, side))
      {
        if (std::any_of(ridged.begin(), ridged.end(), [&](const Ridge& r) {
              return r.edge.IsSame(edge);
            }))
          continue;
        if (!bevel && !piece.rail.s[side].cliffCurve
            && (piece.rail.s[side].trim.IsNull() || piece.rail.law == Law::Full))
          mark(edge, blends[k], partFace.IsBound(edge) ? partFace.Find(edge)
                 : piece.rail.s[side].trim.IsNull()
                 ? piece.rail.s[side].face : piece.rail.s[side].trim,
               JoinRole::Contact, smooth);
        else
          record(edge, JoinRole::Trim, 0);
      }
    }
    for (int end = 0; end < 2; ++end)
    {
      const auto& joint = piece.joint[end];
      if (bevel)
      {
        record(section[k][end], joint.meets == Meets::Stop ? JoinRole::Stop : JoinRole::Miter, 0);
        continue;
      }
      if (joint.meets == Meets::Corner && !parts[joint.piece].empty())
        continue; // with the corner's parts, below
      if (joint.meets == Meets::Joint)
      {
        const bool breaks = smooth == GeomAbs_G2 && supportBreaks(k, end);
        mark(section[k][end], blends[k], blends[joint.piece],
             breaks ? JoinRole::SupportBreak : JoinRole::Joint, breaks ? GeomAbs_G1 : smooth);
      }
      else if (joint.meets == Meets::Fold && !foldAt(k, end).crease.empty())
      {
        // Once, at its first piece: the crease G0, the section beyond it within the one sheet.
        const Folded& f = foldAt(k, end);
        if (f.a == k && f.ea == end)
        {
          record(f.crease[0], JoinRole::Fold, 0);
          for (size_t i = 1; i < f.crease.size(); ++i)
            mark(f.crease[i], blends[f.a], blends[f.b], JoinRole::Joint, smooth);
        }
        continue;
      }
      else if (joint.meets == Meets::Fold)
        mark(section[k][end], blends[k], foldAt(k, end).patch, JoinRole::Fold, GeomAbs_G1);
      else if (const Roll* roll = rollAt(k, end))
        mark(section[k][end], blends[k], roll->patch, JoinRole::Roll, GeomAbs_G1);
      else if (joint.meets == Meets::Corner && !patches[joint.piece].IsNull())
        mark(section[k][end], blends[k], patches[joint.piece], JoinRole::Corner, smooth);
      else
        record(section[k][end], joint.meets == Meets::Stop ? JoinRole::Stop
                             : joint.meets == Meets::Split ? JoinRole::Split : JoinRole::Miter, 0);
      for (const auto& edge : chain[k][end])
        record(edge, JoinRole::Stop, 0);
    }
  }
  for (const auto& [edge, k, part] : partEnds)
    mark(edge, blends[k], part, JoinRole::Corner, GeomAbs_G1);
  for (const auto& [edge, a, b] : partInner)
    mark(edge, a, b, JoinRole::Corner, GeomAbs_G1);
  for (const auto& f : folded)
    for (const auto& run : f.run)
      record(run.edge, JoinRole::Fold, 1);
  for (const auto& roll : rolls)
    record(roll.arc, JoinRole::Roll, 1);
  for (const auto& ridge : ridged)
    record(ridge.edge, JoinRole::Ridge, 0);
  for (const auto& notch : notches)
    record(notch.edge, JoinRole::Stop, 0);
}

void Attach::tolerances()
{
  // Tolerances of what was made (Jackson): each new edge's own gap between its curves, and each
  // corner's distance to the ends of the edges at it.
  for (int i = 1; i <= cuts.Extent(); ++i)
    for (const auto& face : facesOfEdge.FindFromKey(cuts.FindKey(i)))
      editOf(edits, face);
  std::vector<std::pair<TopoDS_Edge, std::vector<TopoDS_Face>>> made = continued;
  for (const auto& notch : notches)
    made.push_back({notch.edge, {blends[notch.blend], notch.face}});
  for (int k = 0; k < n; ++k)
  {
    if (vanished(k))
    {
      made.push_back(
        {section[k][0], {patches[pieces[k].joint[0].piece], patches[pieces[k].joint[1].piece]}});
      continue;
    }
    for (int side = 0; side < 2; ++side)
    {
      const Side& s = pieces[k].rail.s[side];
      for (const auto& edge : sideEdges(k, side))
        if (std::none_of(ridged.begin(), ridged.end(), [&](const Ridge& r) {
              return r.edge.IsSame(edge);
            }))
          made.push_back({edge,
                          {blends[k],
                           partFace.IsBound(edge) ? partFace.Find(edge)
                           : !s.trim.IsNull() ? s.trim
                           : s.cliffCurve   ? beyondCliff(s)
                                            : s.face}});
    }
    for (int end = 0; end < 2; ++end)
    {
      const auto& joint = pieces[k].joint[end];
      if (joint.meets == Meets::Stop)
      {
        if (!section[k][end].IsNull())
          made.push_back({section[k][end], {blends[k], stopAt[k][end]}});
      }
      else if (joint.meets == Meets::Corner && !parts[joint.piece].empty())
        continue; // with the corner's parts, below
      else if (joint.meets == Meets::Corner)
        made.push_back({section[k][end],
                        patches[joint.piece].IsNull()
                          ? std::vector<TopoDS_Face>{blends[k]}
                          : std::vector<TopoDS_Face>{blends[k], patches[joint.piece]}});
      else if (joint.meets == Meets::Fold && !foldAt(k, end).crease.empty())
      {
        const Folded& f = foldAt(k, end);
        if (f.a == k && f.ea == end)
          for (const auto& edge : f.crease)
            made.push_back({edge, {blends[f.a], blends[f.b]}});
      }
      else if (joint.meets == Meets::Fold)
        made.push_back({section[k][end], {blends[k], foldAt(k, end).patch}});
      else if (const Roll* roll = rollAt(k, end))
        made.push_back({section[k][end], {blends[k], roll->patch}});
      else if (joint.meets != Meets::Split && (joint.piece > k || (joint.piece == k && end == 0)))
        made.push_back({section[k][end], {blends[k], blends[joint.piece]}});
    }
  }
  for (const auto& [edge, k, part] : partEnds)
    made.push_back({edge, {blends[k], part}});
  for (const auto& [edge, a, b] : partInner)
    made.push_back({edge, {a, b}});
  for (int i = 1; i <= cuts.Extent(); ++i)
  {
    std::vector<TopoDS_Face> faces;
    for (const auto& face : facesOfEdge.FindFromKey(cuts.FindKey(i)))
      faces.push_back(TopoDS::Face(face));
    made.push_back({TopoDS::Edge(cuts.FindFromIndex(i)), faces});
  }
  for (const auto& f : folded)
    for (const auto& run : f.run)
      made.push_back({run.edge, {run.side->face, f.patch}});
  for (const auto& r : ridged)
    made.push_back({r.edge, {blends[r.k], blends[r.m]}});
  for (const auto& roll : rolls)
    made.push_back({roll.arc, {roll.patch, pieces[roll.a].rail.s[roll.sa].face}});
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> corners;
  for (const auto& piece : contactIntervals)
    for (const auto& side : piece)
      for (const auto& interval : side)
        for (const auto& event : interval.ends)
          corners.Add(event.vertex);
  for (const auto& r : ridged)
    for (const auto& vertex : r.at)
      corners.Add(vertex);
  ShapeAnalysis_Edge analysis;
  for (const auto& [edge, faces] : made)
  {
    if (BRep_Tool::Degenerated(edge))
      continue; // one point: a blend's end where it closes (Rail::pinch)
    // Every curve the edge stores, whatever surface and place it is kept for.
    for (const auto& kept : Handle(BRep_TEdge)::DownCast(edge.TShape())->Curves())
    {
      if (kept->IsCurve3D())
        settle(kept->Curve3D());
      if (kept->IsCurveOnSurface())
        settle(kept->PCurve());
      if (kept->IsCurveOnClosedSurface())
        settle(kept->PCurve2());
    }
    double deviation = 0;
    analysis.CheckSameParameter(edge, deviation);
    if (deviation > BRep_Tool::Tolerance(edge))
      builder.UpdateEdge(edge, deviation);
    for (TopExp_Explorer v(edge, TopAbs_VERTEX); v.More(); v.Next())
    {
      const TopoDS_Vertex vertex = TopoDS::Vertex(v.Current());
      if (!corners.Contains(vertex))
        continue;
      const double t     = BRep_Tool::Parameter(vertex, edge);
      const gp_Pnt point = BRep_Tool::Pnt(vertex);
      double       gap =
        std::max(BRep_Tool::Tolerance(edge), BRepAdaptor_Curve(edge).Value(t).Distance(point));
      for (const auto& face : faces)
        gap = std::max(gap, BRepAdaptor_Curve(edge, face).Value(t).Distance(point));
      if (gap > BRep_Tool::Tolerance(vertex))
        builder.UpdateVertex(vertex, gap);
    }
  }
}

namespace
{
// The box of a face over its own parameter range. BRepBndLib boxes a face that is not planar by its
// whole surface, so a torus patch gets the box of the entire torus and collisions() would test its
// exact distance to every plane that the entire torus reaches. The boxes of quadrics and the torus
// over a parameter range bound the surface there (BndLib), which is all the pre-filter needs.
void faceBox(const TopoDS_Face& face, Bnd_Box& box)
{
  const BRepAdaptor_Surface surface(face, true);
  switch (surface.GetType())
  {
    case GeomAbs_Cylinder:
    case GeomAbs_Cone:
    case GeomAbs_Sphere:
    case GeomAbs_Torus:
      BndLib_AddSurface::Add(surface, 0., box);
      break;
    default:
      BRepBndLib::Add(face, box);
  }
}
//! A positive-weight tensor spline lies in the convex hull of its poles. Separating
//! that hull from the other face's entire supporting plane therefore separates the
//! faces. Unlike world-axis boxes this test is unchanged by a rigid motion.
bool separatedByPlane(const TopoDS_Face& splineFace,
                      const TopoDS_Face& planeFace,
                      double tolerance)
{
  const GeomAdaptor_Surface plane(BRep_Tool::Surface(planeFace));
  if (plane.GetType() != GeomAbs_Plane)
    return false;
  const auto spline = Handle(Geom_BSplineSurface)::DownCast(BRep_Tool::Surface(splineFace));
  if (spline.IsNull())
    return false;
  double u0, u1, v0, v1, su0, su1, sv0, sv1;
  BRepTools::UVBounds(splineFace, u0, u1, v0, v1);
  spline->Bounds(su0, su1, sv0, sv1);
  if (!(u0 >= su0 && u1 <= su1 && v0 >= sv0 && v1 <= sv1))
    return false; // Extrapolation beyond the spline domain need not obey its convex hull.
  // A face's stored edge/vertex geometry can deviate from its support by its
  // own tolerance. Include both boundary envelopes before skipping a shape query.
  tolerance += tol::largest(splineFace) + tol::largest(planeFace);
  const gp_Pln support = plane.Plane();
  const gp_Vec normal(support.Axis().Direction());
  double low = Precision::Infinite(), high = -Precision::Infinite();
  for (int i = 1; i <= spline->NbUPoles(); ++i)
    for (int j = 1; j <= spline->NbVPoles(); ++j)
    {
      if (!(spline->Weight(i, j) > 0))
        return false;
      const gp_Pnt pole = spline->Pole(i, j);
      const gp_Vec delta(support.Location(), pole);
      const double distance = delta.Dot(normal);
      // Resolve subtraction and dot-product roundoff at the coordinates used.
      const double error = tol::rounding(pole.XYZ().Modulus()
                                         + support.Location().XYZ().Modulus());
      if (!std::isfinite(distance) || !std::isfinite(error))
        return false;
      low = std::min(low, distance - error);
      high = std::max(high, distance + error);
      if (low <= tolerance && high >= -tolerance)
        return false;
    }
  return low > tolerance || high < -tolerance;
}
} // namespace

void Attach::collisions()
{
  // Nothing else of the body may meet a blend: a face beside none of them that one reaches
  // would be cut through (a cavity's wall, a far side of a thin part). Refused with every such
  // face, for blend() to lift off the features they are of where it can (liftFeatures).
  std::vector<TopoDS_Shape> reached;
  {
    std::vector<TopoDS_Face> made;
    std::copy_if(blends.begin(), blends.end(), std::back_inserter(made), [](const TopoDS_Face& f) {
      return !f.IsNull();
    });
    std::copy_if(patches.begin(),
                 patches.end(),
                 std::back_inserter(made),
                 [](const TopoDS_Face& f) { return !f.IsNull(); });
    for (const auto& faces : parts)
      made.insert(made.end(), faces.begin(), faces.end());
    for (const auto& f : folded)
      if (!f.patch.IsNull())
        made.push_back(f.patch);
    for (const auto& roll : rolls)
      made.push_back(roll.patch);
    NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> beside;
    for (const auto& edit : edits)
      beside.Add(edit.face);
    for (const auto& piece : pieces)
      for (const auto& side : piece.rail.s)
        beside.Add(side.face);
    // Two faces meet where they come within both their tolerances, as their edges and vertices
    // do: each face lies within its own of its geometry.
    double most = 0;
    for (TopExp_Explorer f(body, TopAbs_FACE); f.More(); f.Next())
      most = std::max(most, tol::own(f.Current()));
    for (const auto& face : made)
    {
      Bnd_Box box;
      faceBox(face, box);
      box.Enlarge(tol::own(face) + most);
      for (TopExp_Explorer f(body, TopAbs_FACE); f.More(); f.Next())
      {
        if (beside.Contains(f.Current())
            || std::any_of(reached.begin(), reached.end(), [&](const TopoDS_Shape& r) {
                 return r.IsSame(f.Current());
               }))
          continue;
        Bnd_Box other;
        faceBox(TopoDS::Face(f.Current()), other);
        if (box.IsOut(other))
          continue;
        const double tolerance = tol::own(face) + tol::own(f.Current());
        if (separatedByPlane(face, TopoDS::Face(f.Current()), tolerance)
            || separatedByPlane(TopoDS::Face(f.Current()), face, tolerance)
            || separated(face, TopoDS::Face(f.Current()), tolerance))
          continue;
        BRepExtrema_DistShapeShape distance(face, f.Current());
        if (distance.IsDone() && distance.Value() <= tol::own(face) + tol::own(f.Current()))
          reached.push_back(f.Current());
      }
    }
  }
  if (!reached.empty())
    throw Collision(op,
                    "the blend reaches further than the faces beside its edges and runs into a "
                    "nonadjacent face of the body",
                    pieces.front().rail.edge,
                    reached);
}

} // namespace BRepExt::blending::detail
