// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Blend: Mitres between two blends, and balls rolled round a re-entrant vertex.

#include <BRepBlendExt_Attach.hxx>

namespace BRepExt::blending::detail
{

void Attach::miters()
{
  // Two edges that meet at an angle, with a face in common and a third edge between their
  // other faces: each blend ends where it meets the other (the mitre), the contacts on the
  // common face end where they cross, and the third edge is cut back to where both other
  // contacts cross it, or, where they cross it at different places, to the farther one.
  for (int k = 0; k < n; ++k)
    for (int ea = 0; ea < 2; ++ea)
    {
      const auto& joint = pieces[k].joint[ea];
      if (joint.meets != Meets::Miter || joint.piece < k)
        continue;
      const int     j = joint.piece, eb = joint.end;
      const Piece&  a = pieces[k];
      const Piece&  b = pieces[j];
      TopoDS_Vertex vertex, other;
      TopExp::Vertices(TopoDS::Edge(a.rail.edge.Oriented(TopAbs_FORWARD)), vertex, other);
      if (ea == 1)
        vertex = other;
      int sa = -1, sb = -1;
      for (int s = 0; s < 2; ++s)
        for (int t = 0; t < 2; ++t)
          if (a.rail.s[s].face.IsSame(b.rail.s[t].face))
          {
            sa = s;
            sb = t;
          }
      if (sa < 0 && a.rail.sigma == b.rail.sigma)
      {
        if (!notch(k, ea, j, eb) && !notch(j, eb, k, ea))
          fail(op,
               unsupported_geometry,
               "two blended edges that meet at an angle without a face in common do not end on "
               "each other's faces, which the blend does not do yet",
               {a.rail.edge});
        continue;
      }
      if (sa < 0 || a.rail.sigma != b.rail.sigma)
        fail(op,
             unsupported_geometry,
             "one convex and one concave blended edge that meet at an angle are blended with their "
             "corner, which the blend does not do yet",
             {a.rail.edge});
      const int   oa = 1 - sa, ob = 1 - sb;
      TopoDS_Edge third;
      for (const auto& shape : edgesOfVertex.FindFromKey(vertex))
      {
        if (shape.IsSame(a.rail.edge) || shape.IsSame(b.rail.edge))
          continue;
        if (!third.IsNull() || !onFace(shape, a.rail.s[oa].face)
            || !onFace(shape, b.rail.s[ob].face))
          fail(op,
               unsupported_geometry,
               "where two blended edges meet among more faces, the corner is blended too, which "
               "the blend does not do yet",
               {a.rail.edge});
        third = TopoDS::Edge(shape);
      }
      if (third.IsNull())
        fail(op,
             unsupported_geometry,
             "two blended edges that meet at an angle have no edge between their other faces",
             {a.rail.edge});
      const auto [a0, a1] = halfOf(a, ea);
      const auto [b0, b1] = halfOf(b, eb);
      const auto q        = commonCrossing(k, ea, sa);
      const auto pa = crossing(a.fitted.contact[oa], a0, a1, a.uAt(ea), third, a.rail.s[oa].face);
      const auto pb = crossing(b.fitted.contact[ob], b0, b1, b.uAt(eb), third, b.rail.s[ob].face);
      if (!q || !pa || !pb)
        fail(op,
             unsupported_geometry,
             "the blends of two edges that meet at an angle do not reach each other within their "
             "faces",
             {a.rail.edge});
      if (runsAcross(k, ea, false))
      {
        rollRound(k, ea, j, eb, sa, sb, vertex, third, *pa, *pb);
        continue;
      }
      const gp_Pnt Q    = a.sheet.contact(sa)->Value(a.sheet.rail(q->first)),
                   P    = a.sheet.contact(oa)->Value(a.sheet.rail(pa->first));
      const gp_Pnt Pb   = b.sheet.contact(ob)->Value(b.sheet.rail(pb->first));
      const double gapQ = Q.Distance(b.sheet.contact(sb)->Value(b.sheet.rail(q->second)));
      const double gapP = P.Distance(Pb);
      if (gapQ > a.rail.resolved(limit::together))
        fail(op,
             no_solution,
             "the contacts of two blends that meet at an angle do not cross at one place on the "
             "face they share",
             {a.rail.edge});
      TopoDS_Vertex vq;
      builder.MakeVertex(vq, Q, std::max(tolerance, gapQ));
      contactEnd(k, sa, ea).vertex = contactEnd(j, sb, eb).vertex = vq;
      contactEnd(k, sa, ea).parameter                          = q->first;
      contactEnd(j, sb, eb).parameter                          = q->second;
      // Between cylinders and tori the meeting is a conic OpenCascade finds exactly; between
      // fitted sheets it is marched, as they touch at the inner corner. Chamfers do not touch
      // there: they cross, and OpenCascade intersects them.
      const gp_Pnt2d aQ = a.sheet.uv(a.sheet.rail(q->first), a.sheet.across[sa]),
                     bQ = b.sheet.uv(b.sheet.rail(q->second), b.sheet.across[sb]);
      const bool fittedSheets =
        (!conical(a.sheet) || !conical(b.sheet)) && a.rail.law != Law::Chamfer;
      const auto meeting = [&](const Sheet&         one,
                               const gp_Pnt2d&      oneFirst,
                               const gp_Pnt2d&      oneSecond,
                               const Sheet&         two,
                               const gp_Pnt2d&      twoFirst,
                               const gp_Pnt2d&      twoSecond,
                               const TopoDS_Vertex& first) {
        auto met = fittedSheets
                     ? marchedMeeting(one, oneFirst, oneSecond, two, twoFirst, twoSecond, first, vq)
                     : meetingEdge(placedSheet(one),
                                   oneFirst,
                                   placedSheet(two),
                                   twoFirst,
                                   first,
                                   vq,
                                   a.rail.scale);
        // A looping contact can reach this same 3D corner at several parameter pairs.
        // The nearest pair need not be on the incident intersection branch. Only a pair
        // on this same corner whose whole connecting edge passes the marcher's spatial
        // and parameter-curve certificates can replace it; the contact trims follow it.
        if (!met && fittedSheets)
        {
          for (const auto& alternative : meetings(between(a.fitted.contact[sa], a0, a1),
                                                  between(b.fitted.contact[sb], b0, b1),
                                                  a.rail.s[sa].surface))
          {
            if (std::abs(alternative.first - q->first) <= tol::parameter
                && std::abs(alternative.second - q->second) <= tol::parameter)
              continue;
            const gp_Pnt2d aEnd = a.sheet.uv(a.sheet.rail(alternative.first), a.sheet.across[sa]),
                           bEnd = b.sheet.uv(b.sheet.rail(alternative.second), b.sheet.across[sb]);
            if (!(a.sheet.surface->Value(aEnd.X(), aEnd.Y()).Distance(Q) <= tolerance
                  && b.sheet.surface->Value(bEnd.X(), bEnd.Y()).Distance(Q) <= tolerance))
              continue;
            met = marchedMeeting(one,
                                 oneFirst,
                                 &one == &a.sheet ? aEnd : bEnd,
                                 two,
                                 twoFirst,
                                 &two == &a.sheet ? aEnd : bEnd,
                                 first,
                                 vq);
            if (!met)
              continue;
            contactEnd(k, sa, ea).parameter = alternative.first;
            contactEnd(j, sb, eb).parameter = alternative.second;
            break;
          }
        }
        if (!met)
          fail(op,
               no_solution,
               "the blends of two edges that meet at an angle do not meet each other between their "
               "contacts",
               {a.rail.edge});
        return *met;
      };
      if (gapP <= a.rail.resolved(limit::together))
      {
        TopoDS_Vertex vp;
        builder.MakeVertex(vp, P, std::max(tolerance, gapP));
        contactEnd(k, oa, ea).vertex = contactEnd(j, ob, eb).vertex = vp;
        contactEnd(k, oa, ea).parameter                          = pa->first;
        contactEnd(j, ob, eb).parameter                          = pb->first;
        cut(third, vertex, vp, pa->second);
        const gp_Pnt2d aP = a.sheet.uv(a.sheet.rail(pa->first), a.sheet.across[oa]),
                       bP = b.sheet.uv(b.sheet.rail(pb->first), b.sheet.across[ob]);
        section[k][ea] = section[j][eb] = meeting(a.sheet, aP, aQ, b.sheet, bP, bQ, vp);
        continue;
      }
      // The other contacts cross the third edge at different places (supports that lean
      // differently at the vertex): as at a notch (Attach::notch), the blend whose contact
      // crosses it farther from the vertex wins it. The edge is cut back to that crossing,
      // the winner's end runs from there along its sheet on the loser's face to where the
      // loser's contact meets the winner's sheet, and the loser's contact ends there; from
      // there to the face they share the sheets meet.
      const int     p[2] = {k, j}, at[2] = {ea, eb}, own[2] = {sa, sb}, far[2] = {oa, ob};
      const Meeting crosses[2] = {*pa, *pb};
      const gp_Pnt  V          = BRep_Tool::Pnt(vertex);
      const int     w = P.Distance(V) > Pb.Distance(V) ? 0 : 1, l = 1 - w;
      const Piece&  winner = pieces[p[w]];
      const Piece&  loser  = pieces[p[l]];
      const auto    r      = meetsSheet(p[l], at[l], far[l], p[w]);
      if (!r)
        fail(op,
             unsupported_geometry,
             "the blends of two edges that meet at an angle do not reach each other within their "
             "faces",
             {a.rail.edge});
      const TopoDS_Face& face = loser.rail.s[far[l]].face;
      TopoDS_Vertex      vx, vr;
      builder.MakeVertex(vx,
                         winner.sheet.contact(far[w])->Value(winner.sheet.rail(crosses[w].first)),
                         tolerance);
      builder.MakeVertex(vr, loser.sheet.contact(far[l])->Value(loser.sheet.rail(*r)), tolerance);
      double         f, l2;
      const gp_Pnt2d xOnFace =
        BRep_Tool::CurveOnSurface(third, face, f, l2)->Value(crosses[w].second);
      const auto along = meetingEdge(
        placedSheet(winner.sheet),
        winner.sheet.uv(winner.sheet.rail(crosses[w].first), winner.sheet.across[far[w]]),
        placedFace(face),
        xOnFace,
        vx,
        vr,
        winner.rail.scale);
      if (!along)
        fail(
          op,
          no_solution,
          "the blend that wins the edge between two mitred blends does not meet the other's face",
          {a.rail.edge});
      const gp_Pnt2d rOnLoser = loser.sheet.uv(loser.sheet.rail(*r), loser.sheet.across[far[l]]);
      const auto     middle   = meeting(loser.sheet,
                                  rOnLoser,
                                  l == 0 ? aQ : bQ,
                                  winner.sheet,
                                  onSheet(p[w], at[w], BRep_Tool::Pnt(vr)),
                                  w == 0 ? aQ : bQ,
                                  vr);
      contactEnd(p[w], far[w], at[w]).vertex = vx;
      contactEnd(p[w], far[w], at[w]).parameter      = crosses[w].first;
      contactEnd(p[l], far[l], at[l]).vertex = vr;
      contactEnd(p[l], far[l], at[l]).parameter      = *r;
      cut(third, vertex, vx, crosses[w].second);
      editOf(edits, face).bridges.push_back(*along);
      notches.push_back({p[w], *along, face});
      section[k][ea] = section[j][eb] = middle;
      chain[p[w]][at[w]]              = own[w] == 0 ? std::vector<TopoDS_Edge>{middle, *along}
                                                    : std::vector<TopoDS_Edge>{*along, middle};
    }
}

//! Two mitred blends whose contacts, carried on past a re-entrant vertex, would run across each
//! other before they cross (Attach::miters): the ball rolls round the vertex from one blend's end
//! to the other's, as it does round any obstacle it cannot pass. The other face reaches the third
//! edge at X, where both must reach it; every ball between touches
//! the face they share and passes through X: their centres lie on a circle about the shared
//! face's normal through X, and they sweep a torus that closes to a point at X. The third edge is
//! cut back to X, and the shared face takes the arc of the balls' contacts between the blends.
//! For the round on a flat shared face, where the balls' envelope is that torus; elsewhere the
//! blend declines.
void Attach::rollRound(int                  k,
                       int                  ea,
                       int                  j,
                       int                  eb,
                       int                  sa,
                       int                  sb,
                       const TopoDS_Vertex& vertex,
                       const TopoDS_Edge&   third,
                       const Meeting&       pa,
                       const Meeting&       pb)
{
  const Piece&       a  = pieces[k];
  const Piece&       b  = pieces[j];
  const int          oa = 1 - sa, ob = 1 - sb;
  const TopoDS_Face& common = a.rail.s[sa].face;
  if (!round(a.rail) || !round(b.rail)
      || std::abs(a.rail.s[0].distance - b.rail.s[0].distance) > limit::equal * a.rail.scale)
    fail(op,
         unsupported_geometry,
         "two blends that cannot be mitred at a re-entrant vertex roll round it only as rounds of "
         "one radius",
         {a.rail.edge});
  if (GeomAdaptor_Surface(BRep_Tool::Surface(common)).GetType() != GeomAbs_Plane)
    fail(op,
         unsupported_geometry,
         "two blends that cannot be mitred at a re-entrant vertex roll round it only on a flat "
         "face, which the blend does not do yet",
         {a.rail.edge});
  const gp_Pnt P = a.sheet.contact(oa)->Value(a.sheet.rail(pa.first));
  if (P.Distance(b.sheet.contact(ob)->Value(b.sheet.rail(pb.first)))
      > a.rail.resolved(limit::together))
    fail(op,
         unsupported_geometry,
         "two blends that cannot be mitred at a re-entrant vertex reach the edge between them at "
         "different places, which the blend does not roll round",
         {a.rail.edge});
  // Each end ball: its foot on the shared face, and its centre.
  const double r = a.rail.s[0].distance;
  const auto   ball =
    [&](const Piece& piece, int side, double at, gp_Pnt& foot, gp_Pnt& centre, gp_Vec& up) {
      const gp_Pnt2d uv = piece.fitted.contact[side]->Value(at);
      Jet            jet;
      if (!jetAt(piece.rail.s[side], uv.X(), uv.Y(), jet))
        fail(op,
             no_solution,
             "the ball has no normal on the face where it rolls round a vertex",
             {piece.rail.edge});
      foot   = piece.sheet.contact(side)->Value(piece.sheet.rail(at));
      centre = jet.p.Translated(jet.n * (piece.rail.sigma * r));
      up     = jet.n * -piece.rail.sigma; // from the centres to the face
    };
  gp_Pnt Ta, Tb, Ca, Cb;
  gp_Vec up, upB;
  ball(a, sa, pa.first, Ta, Ca, up);
  ball(b, sb, pb.first, Tb, Cb, upB);
  // The circle of centres: about the line through X along the face's normal, in the plane of the
  // centres.
  const gp_Pnt X0 = P.Translated(up * gp_Vec(P, Ca).Dot(up));
  const gp_Vec toA(X0, Ca), toB(X0, Cb);
  const double rho = toA.Magnitude();
  if (!(rho > limit::degenerate * a.rail.scale)
      || std::abs(toB.Magnitude() - rho) > a.rail.resolved(limit::together))
    fail(op,
         no_solution,
         "the balls at two blends' ends do not roll round the vertex between them",
         {a.rail.edge});
  const gp_Ax3 frame(X0, gp_Dir(up), gp_Dir(toA));
  const double ub =
    std::atan2(toB.Dot(gp_Vec(frame.YDirection())), toB.Dot(gp_Vec(frame.XDirection())));
  // X on each meridian: at the centres' distance r from it, h above their plane.
  const double h  = gp_Vec(X0, P).Dot(up);
  double       vX = std::atan2(h, -rho);
  if (vX < M_PI / 2)
    vX += 2 * M_PI;
  Roll roll{k,
            ea,
            j,
            eb,
            sa,
            sb,
            vertex,
            {},
            {},
            new Geom_ToroidalSurface(frame, rho, r),
            0,
            ub,
            vX,
            {}};
  builder.MakeVertex(roll.X, P, tolerance);
  TopoDS_Vertex ta, tb;
  builder.MakeVertex(ta, Ta, tolerance);
  builder.MakeVertex(tb, Tb, tolerance);
  contactEnd(k, oa, ea).vertex = contactEnd(j, ob, eb).vertex = roll.X;
  contactEnd(k, sa, ea).vertex                     = ta;
  contactEnd(j, sb, eb).vertex                     = tb;
  contactEnd(k, oa, ea).parameter = contactEnd(k, sa, ea).parameter = pa.first;
  contactEnd(j, ob, eb).parameter = contactEnd(j, sb, eb).parameter = pb.first;
  cut(third, vertex, roll.X, pa.second);
  for (const auto& [m, end, at] : {std::tuple{k, ea, pa.first}, std::tuple{j, eb, pb.first}})
  {
    const Sheet& sheet = pieces[m].sheet;
    section[m][end]    = isoSection(sheet, sheet.rail(at), contactEnd(m, 0, end).vertex, contactEnd(m, 1, end).vertex);
    builder.UpdateEdge(section[m][end],
                       sheet.sectionUV(sheet.rail(at)),
                       sheet.surface,
                       TopLoc_Location(),
                       tolerance);
  }
  // The balls' contacts on the shared face, between the blends' ends: on the torus where v is a
  // quarter turn, its parameter the angle about the axis.
  const Handle(Geom_Circle) arc =
    new Geom_Circle(gp_Ax2(X0.Translated(up * r), gp_Dir(up), gp_Dir(toA)), rho);
  const auto [surface, location] = surfaceOf(common);
  const auto onCommon =
    GeomProjLib::Curve2d(arc, std::min(0.0, ub), std::max(0.0, ub), BRep_Tool::Surface(common));
  if (onCommon.IsNull())
    fail(op,
         no_solution,
         "the balls rolled round a vertex do not lie on the face they share",
         {a.rail.edge});
  roll.arc = edgeBetween(arc,
                         ta,
                         0,
                         tb,
                         ub,
                         {{segment2d(gp_Pnt2d(std::min(0.0, ub), M_PI / 2),
                                     gp_Pnt2d(std::max(0.0, ub), M_PI / 2),
                                     std::min(0.0, ub),
                                     std::max(0.0, ub)),
                           roll.torus,
                           {}},
                          {onCommon, surface, location}},
                         tolerance);
  editOf(edits, common).bridges.push_back(roll.arc);
  rolls.push_back(roll);
}

void Attach::rollPatches()
{
  for (auto& roll : rolls)
  {
    const auto&  torus = roll.torus;
    const gp_Ax3 frame = torus->Position();
    const double rho   = torus->MajorRadius();
    // A section, the meridian at u: its v at each place from the angle about the ball's centre.
    const auto meridian = [&](const TopoDS_Edge& edge, double u) -> Handle(Geom2d_Curve) {
      double       lo, hi;
      const auto   curve = BRep_Tool::Curve(edge, lo, hi);
      const gp_Vec out =
        gp_Vec(frame.XDirection()) * std::cos(u) + gp_Vec(frame.YDirection()) * std::sin(u);
      const gp_Pnt                 C     = frame.Location().Translated(out * rho);
      constexpr int                count = 33;
      NCollection_Array1<gp_Pnt2d> values(1, count);
      NCollection_Array1<double>   at(1, count);
      for (int i = 1; i <= count; ++i)
      {
        const double t = lo + (hi - lo) * (i - 1) / (count - 1);
        const gp_Vec to(C, curve->Value(t));
        double       v = std::atan2(to.Dot(gp_Vec(frame.Direction())), to.Dot(out));
        if (v < M_PI / 2 - tol::angle)
          v += 2 * M_PI;
        at(i)     = t;
        values(i) = gp_Pnt2d(u, v);
      }
      Geom2dAPI_Interpolate fit(new NCollection_HArray1<gp_Pnt2d>(values),
                                new NCollection_HArray1<double>(at),
                                false,
                                tol::parameter);
      fit.Perform();
      return fit.IsDone() ? Handle(Geom2d_Curve)(fit.Curve()) : Handle(Geom2d_Curve)();
    };
    const TopoDS_Edge& sa = section[roll.a][roll.ea];
    const TopoDS_Edge& sb = section[roll.b][roll.eb];
    const auto         pa = meridian(sa, roll.ua), pb = meridian(sb, roll.ub);
    if (pa.IsNull() || pb.IsNull())
      fail(op,
           no_solution,
           "a blend's end does not lie on the ball rolled round its vertex",
           {pieces[roll.a].rail.edge});
    for (const auto& [edge, pc] : {std::pair{sa, pa}, std::pair{sb, pb}})
    {
      double f, l;
      BRep_Tool::Range(edge, f, l);
      builder.UpdateEdge(edge, pc, torus, TopLoc_Location(), tolerance);
      builder.Range(edge, f, l);
    }
    // Round the loop: X along the first section to the shared face, the arc, back down the
    // second, and the point.
    const TopoDS_Edge point =
      degenerateEdge(segment2d(gp_Pnt2d(roll.ub, roll.vX), gp_Pnt2d(roll.ua, roll.vX), 0, 1),
                     torus,
                     TopLoc_Location(),
                     roll.X,
                     0,
                     1,
                     tolerance);
    TopoDS_Wire wire;
    builder.MakeWire(wire);
    const auto from = [](const TopoDS_Edge& edge, const TopoDS_Vertex& start) {
      return TopExp::FirstVertex(edge).IsSame(start) ? edge : TopoDS::Edge(edge.Reversed());
    };
    builder.Add(wire, from(sa, roll.X));
    builder.Add(wire, from(roll.arc, contactEnd(roll.a, roll.sa, roll.ea).vertex));
    builder.Add(wire, from(sb, contactEnd(roll.b, roll.sb, roll.eb).vertex));
    builder.Add(wire, point);
    wire.Closed(true);
    TopoDS_Face face;
    builder.MakeFace(face, torus, tolerance);
    builder.Add(face, wire);
    // Each edge's real gap to the torus is its tolerance (Jackson), as the fit leaves it.
    ShapeAnalysis_Edge analysis;
    for (const TopoDS_Edge& edge : {sa, roll.arc, sb})
    {
      double deviation = 0;
      analysis.CheckSameParameter(edge, face, deviation);
      if (deviation > pieces[roll.a].rail.resolved(limit::together))
        fail(op,
             no_solution,
             "a blend's end does not lie on the ball rolled round its vertex",
             {pieces[roll.a].rail.edge});
      if (deviation > BRep_Tool::Tolerance(edge))
        builder.UpdateEdge(edge, deviation);
      builder.SameParameter(edge, true);
    }
    if (!BRepCheck_Analyzer(face).IsValid())
    {
      // The other way round the loop.
      TopoDS_Wire back = TopoDS::Wire(wire.Reversed());
      builder.MakeFace(face, torus, tolerance);
      builder.Add(face, back);
      if (!BRepCheck_Analyzer(face).IsValid())
        fail(op,
             no_solution,
             "the ball rolled round a re-entrant vertex makes no valid face",
             {pieces[roll.a].rail.edge});
    }
    // Facing out as the blends do: across a shared edge the two faces run it opposite ways.
    const auto use = [](const TopoDS_Face& in, const TopoDS_Edge& edge) {
      for (TopExp_Explorer e(in, TopAbs_EDGE); e.More(); e.Next())
        if (e.Current().IsSame(edge))
          return e.Current().Orientation();
      return TopAbs_EXTERNAL;
    };
    if (use(face, sa) == use(blends[roll.a], sa))
      face.Reverse();
    roll.patch = face;
  }
}

} // namespace BRepExt::blending::detail
