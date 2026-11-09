// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Blend: The corners where three or more blends meet.

#include <BRepBlendExt_Attach.hxx>

namespace BRepExt::blending::detail
{

void Attach::pinchCorner(const Corner& at)
{
  const Piece& first = pieces[at.ends[0].first];
  // Every blend closes to the vertex: its end there is that point, a degenerate edge
  // across its sheet, and the corner has no patch.
  const TopoDS_Vertex point = footAt(BRep_Tool::Pnt(at.vertex), first.rail.scale);
  for (const auto& [k, end] : at.ends)
  {
    const Sheet& sheet = pieces[k].sheet;
    const double along = pieces[k].uAt(end);
    if (sheet.contact(0)->Value(sheet.rail(along)).Distance(BRep_Tool::Pnt(point))
          > pieces[k].rail.resolved(limit::together)
        || sheet.contact(1)->Value(sheet.rail(along)).Distance(BRep_Tool::Pnt(point))
             > pieces[k].rail.resolved(limit::together))
      fail(op,
           no_solution,
           "a blend between faces tangent at its end does not close to the vertex",
           {pieces[k].rail.edge});
    corner[k][0][end] = corner[k][1][end] = point;
    u[k][0][end] = u[k][1][end] = along;
    section[k][end]             = degenerateEdge(sheet.sectionUV(sheet.rail(along)),
                                                 sheet.surface,
                                                 TopLoc_Location(),
                                                 point,
                                                 std::min(sheet.across[0], sheet.across[1]),
                                                 std::max(sheet.across[0], sheet.across[1]),
                                                 tolerance);
  }
}

void Attach::ballSections(const Corner&                     at,
                          const std::vector<TopoDS_Face>&   faces,
                          const std::vector<TopoDS_Vertex>& foot)
{
  // Each blend ends at the section through the ball's centre, from foot to foot.
  for (const auto& [k, end] : at.ends)
  {
    const Piece& piece = pieces[k];
    double       ends[2];
    for (int side = 0; side < 2; ++side)
    {
      int i = 0;
      while (!faces[i].IsSame(piece.rail.s[side].face))
        ++i;
      corner[k][side][end] = foot[i];
      GeomAPI_ProjectPointOnCurve along(BRep_Tool::Pnt(foot[i]), piece.sheet.contact(side));
      if (along.NbPoints() == 0 || along.LowerDistance() > piece.rail.resolved(limit::together))
        fail(op, no_solution, "a blend's contact misses the corner's ball", {piece.rail.edge});
      ends[side] = (along.LowerDistanceParameter() - piece.sheet.alpha) / piece.sheet.beta;
    }
    if (railStationSeparation(piece.rail, ends[0], ends[1]) > piece.rail.resolved(limit::together))
      fail(op,
           no_solution,
           "a blend's contacts reach the corner's ball at different places",
           {piece.rail.edge});
    u[k][0][end] = ends[0];
    u[k][1][end] = ends[1];
    // Where the corner at the blend's other end has the same ball, the blend has no length:
    // the corners share its section.
    if (const int other = 1 - end; pieces[k].joint[other].meets == Meets::Corner
                                   && !section[k][other].IsNull()
                                   && corner[k][0][other].IsSame(corner[k][0][end])
                                   && corner[k][1][other].IsSame(corner[k][1][end]))
    {
      section[k][end] = section[k][other];
      continue;
    }
    const double along = (ends[0] + ends[1]) / 2;
    section[k][end] =
      isoSection(piece.sheet, piece.sheet.rail(along), corner[k][0][end], corner[k][1][end]);
    builder.UpdateEdge(section[k][end],
                       piece.sheet.sectionUV(piece.sheet.rail(along)),
                       piece.sheet.surface,
                       TopLoc_Location(),
                       tolerance);
  }
}

gp_Ax3 Attach::ballFrame(const Corner&                   at,
                         const std::vector<TopoDS_Face>& faces,
                         const gp_Pnt&                   ball,
                         const std::vector<gp_Vec>&      normal,
                         double                          sigma,
                         bool&                           square,
                         int&                            poleFoot)
{
  const size_t count = faces.size();
  // The sphere's frame. Where the edges are square to each other, as OpenCascade's own corners
  // are made, its pole is the ball's foot on the face across from the first edge and its x
  // axis an end of that edge's section: that section is its equator and the other two are
  // meridians, every side an isoline. Elsewhere the poles and seam are kept away from the patch.
  gp_Vec e[3];
  square = count == 3;
  for (int i = 0; i < 3 && square; ++i)
    square &= leaving(pieces[at.ends[i].first].rail.edge, at.vertex, e[i], at.ends[i].second);
  square = square && std::abs(e[0].Dot(e[1])) < tol::fineAngle
           && std::abs(e[0].Dot(e[2])) < tol::fineAngle
           && std::abs(e[1].Dot(e[2])) < tol::fineAngle;
  gp_Ax3 frame;
  poleFoot = -1;
  if (square)
  {
    const Piece& p0         = pieces[at.ends[0].first];
    int          equator[2] = {-1, -1};
    for (int i = 0; i < 3; ++i)
    {
      if (faces[i].IsSame(p0.rail.s[0].face))
        equator[0] = i;
      else if (faces[i].IsSame(p0.rail.s[1].face))
        equator[1] = i;
      else
        poleFoot = i;
    }
    // From the centre to each foot is the face's normal, exactly: against it at a convex
    // corner, along it at a concave one.
    const gp_Dir up(normal[poleFoot] * -sigma);
    // x towards the equator foot from which the other is a quarter turn ahead.
    gp_Dir x(normal[equator[0]] * -sigma);
    if (gp_Vec(up).Crossed(gp_Vec(x)).Dot(normal[equator[1]] * -sigma) < 0)
      x = gp_Dir(normal[equator[1]] * -sigma);
    frame = fixedFrame(gp_Ax3(ball, up, x));
  }
  else
  {
    gp_Vec sum;
    for (const auto& n : normal)
      sum += n;
    const gp_Vec middle = sum.Normalized() * -sigma;
    gp_Vec       pole   = middle.Crossed(normal[0]);
    if (pole.Magnitude() < tol::coarseAngle)
      pole = middle.Crossed(normal[1]);
    frame = fixedFrame(gp_Ax3(ball, gp_Dir(pole), gp_Dir(middle.Reversed())));
  }
  return frame;
}

void Attach::ballPatchFace(const Corner&                     at,
                           size_t                            c,
                           const std::vector<TopoDS_Vertex>& foot,
                           const gp_Ax3&                     frame,
                           bool                              square,
                           int                               poleFoot)
{
  const Piece&                        first = pieces[at.ends[0].first];
  const double                        d = first.rail.s[0].distance, sigma = first.rail.sigma;
  const size_t                        count  = foot.size();
  const Handle(Geom_SphericalSurface) sphere = new Geom_SphericalSurface(frame, d);
  TopoDS_Face&                        patch  = patches[c];
  builder.MakeFace(patch, sphere, tolerance);
  // Its boundary: the blends' sections, from foot to foot, each laid on the sphere and moved
  // by whole turns to where its middle is.
  const auto uvOf = [&](const gp_Pnt& point) {
    double su, sv;
    ElSLib::Parameters(gp_Sphere(frame, d), point, su, sv);
    if (square && su > 1.5 * M_PI)
      su -= 2 * M_PI;
    return gp_Pnt2d(su, sv);
  };
  std::vector<TopoDS_Edge> loop;
  TopoDS_Vertex            at0 = foot[0];
  std::vector<bool>        used(count, false);
  for (size_t step = 0; step < count; ++step)
    for (size_t i = 0; i < count; ++i)
    {
      if (used[i])
        continue;
      const TopoDS_Edge& edge = section[at.ends[i].first][at.ends[i].second];
      TopoDS_Vertex      a, z;
      TopExp::Vertices(edge, a, z);
      if (!a.IsSame(at0) && !z.IsSame(at0))
        continue;
      double f, l;
      BRep_Tool::Range(edge, f, l);
      const auto curve = BRep_Tool::Curve(edge, f, l);
      const auto uv    = laidOn(curve, f, l, sphere);
      if (uv.IsNull())
        fail(op,
             no_solution,
             "a blend's end cannot be laid on the corner's ball",
             {first.rail.edge});
      builder.UpdateEdge(edge,
                         aligned(uv, (f + l) / 2, uvOf(curve->Value((f + l) / 2)), sphere),
                         sphere,
                         TopLoc_Location(),
                         tolerance);
      builder.Range(edge, f, l);
      loop.push_back(a.IsSame(at0) ? edge : TopoDS::Edge(edge.Reversed()));
      at0     = a.IsSame(at0) ? z : a;
      used[i] = true;
      break;
    }
  if (loop.size() != count)
    fail(op,
         no_solution,
         "the ends of the blends do not close around their corner",
         {first.rail.edge});
  // Where two meridians meet at the pole, the side between them in the parameters is the
  // pole itself: a degenerate edge. Elsewhere two sections end at one foot, their curves on
  // the sphere within that place's tolerance of each other (a fitted sheet's section lies
  // on the ball to its fit).
  const auto endOf = [&](const TopoDS_Edge& edge, bool last) {
    double     f, l;
    const auto uv       = BRep_Tool::CurveOnSurface(edge, patch, f, l);
    const bool reversed = edge.Orientation() == TopAbs_REVERSED;
    return uv->Value(last != reversed ? l : f);
  };
  std::vector<TopoDS_Edge> closed;
  for (size_t i = 0; i < loop.size(); ++i)
  {
    closed.push_back(loop[i]);
    const gp_Pnt2d from = endOf(loop[i], true), to = endOf(loop[(i + 1) % loop.size()], false);
    if (from.Distance(to) < tol::uv)
      continue;
    const bool atPole =
      square && poleFoot >= 0 && TopExp::LastVertex(loop[i], true).IsSame(foot[poleFoot]);
    if (!atPole)
    {
      const double gap =
        d * std::hypot((to.X() - from.X()) * std::cos(from.Y()), to.Y() - from.Y());
      if (gap > first.rail.resolved(limit::together))
        fail(op,
             no_solution,
             "the ends of the blends do not close around their corner's ball",
             {first.rail.edge});
      continue;
    }
    closed.push_back(degenerateEdge(new Geom2d_Line(from, gp_Dir2d(gp_Vec2d(from, to))),
                                    sphere,
                                    TopLoc_Location(),
                                    foot[poleFoot],
                                    0,
                                    from.Distance(to),
                                    tolerance));
  }
  // Counter-clockwise in the sphere's parameters: the area the loop's parameter curves enclose.
  double area = 0;
  for (const auto& edge : closed)
  {
    double     f, l;
    const auto uv       = BRep_Tool::CurveOnSurface(edge, patch, f, l);
    gp_Pnt2d   previous = uv->Value(edge.Orientation() == TopAbs_REVERSED ? l : f);
    for (int i = 1; i <= 16; ++i)
    {
      const double t =
        edge.Orientation() == TopAbs_REVERSED ? l - (l - f) * i / 16 : f + (l - f) * i / 16;
      const gp_Pnt2d next = uv->Value(t);
      area += previous.X() * next.Y() - next.X() * previous.Y();
      previous = next;
    }
  }
  TopoDS_Wire wire;
  builder.MakeWire(wire);
  for (const auto& edge : closed)
    builder.Add(wire, area > 0 ? edge : TopoDS::Edge(edge.Reversed()));
  wire.Closed(true);
  builder.Add(patch, wire);
  // Facing out of the body: away from the centre at a convex corner, towards it at a concave.
  if (sigma > 0)
    patch.Reverse();
}

void Attach::cornerPatch(size_t c)
{
  const auto&  at    = cornersAt[c];
  const Piece& first = pieces[at.ends[0].first];
  if (at.pinch)
  {
    pinchCorner(at);
    return;
  }
  std::vector<TopoDS_Face> faces;
  for (const auto& [k, end] : at.ends)
    for (const auto& side : pieces[k].rail.s)
      if (std::none_of(faces.begin(), faces.end(), [&](const TopoDS_Face& f) {
            return f.IsSame(side.face);
          }))
        faces.push_back(side.face);
  // One law and one side throughout, each face between two of the blends.
  bool alike = faces.size() == at.ends.size(), even = true;
  for (const auto& [k, end] : at.ends)
  {
    const Rail& r = pieces[k].rail;
    alike &= r.law == first.rail.law && std::abs(r.rho - first.rail.rho) < limit::equal
             && r.sigma == first.rail.sigma;
    even &= std::abs(r.s[0].distance - first.rail.s[0].distance) < limit::equal * r.scale
            && std::abs(r.s[1].distance - first.rail.s[0].distance) < limit::equal * r.scale;
  }
  if (!alike)
    fail(op,
         unsupported_geometry,
         "where three blended edges meet, the blend makes the corner of one law so far",
         {first.rail.edge});
  const bool round =
    even && first.rail.law == Law::Conic && std::abs(first.rail.rho - 0.5) < limit::equal;
  if (!round)
  {
    if (at.ends.size() > 3 && first.rail.law == Law::Conic)
    {
      splitCorner(at, c);
      return;
    }
    if (at.ends.size() != 3)
      fail(op,
           unsupported_geometry,
           "where more than three G2, conic or chamfer blends meet, the corner is blended with "
           "them, which the blend does not do yet",
           {first.rail.edge});
    patches[c] = lawCorner(at, c);
    return;
  }
  // The ball resting on all the faces, from the corner's vertex.
  const double      d = first.rail.s[0].distance, sigma = first.rail.sigma;
  const size_t      count = faces.size();
  std::vector<Side> sides(count);
  for (size_t i = 0; i < count; ++i)
  {
    sides[i].face    = faces[i];
    sides[i].surface = BRep_Tool::Surface(faces[i]);
    sides[i].outward = faces[i].Orientation() == TopAbs_REVERSED ? -1 : 1;
  }
  const auto rests = restingBall(sides, BRep_Tool::Pnt(at.vertex), d, sigma, first.rail.scale);
  if (!rests)
    fail(op,
         no_solution,
         "no ball of this size rests on the faces at the corner",
         {first.rail.edge});
  const gp_Pnt               ball   = rests->centre;
  const std::vector<gp_Vec>& normal = rests->normals;
  std::vector<TopoDS_Vertex> foot(count);
  for (size_t i = 0; i < count; ++i)
    foot[i] = footAt(rests->feet[i], first.rail.scale);
  // Each blend ends at the section through the ball's centre, from foot to foot.
  ballSections(at, faces, foot);
  // The sphere's frame. Where the edges are square to each other, as OpenCascade's own corners
  // are made, its pole is the ball's foot on the face across from the first edge and its x
  // axis an end of that edge's section: that section is its equator and the other two are
  // meridians, every side an isoline. Elsewhere the poles and seam are kept away from the patch.
  bool         square   = false;
  int          poleFoot = -1;
  const gp_Ax3 frame    = ballFrame(at, faces, ball, normal, sigma, square, poleFoot);
  ballPatchFace(at, c, foot, frame, square, poleFoot);
}

void Attach::cornerPatches()
{
  // Three or more edges that meet: the rolling ball's corner. At equal distances the round's
  // corner is the ball itself where one rests on all the faces (three always have one; more need
  // not).
  for (size_t c = 0; c < cornersAt.size(); ++c)
    cornerPatch(c);
}

// The corner of three blends of one law between three planes at any angles: the vertex
// blend for box-type corners (Várady & Rockwood; Hartmann), a tensor product of the section's
// own law. In the basis of the edges' directions e away from the corner V, with Q_ij the
// corner on the face of edges i and j where their contacts cross, Q_ij = V + p_ij e_i + p_ji
// e_j, X the section's coordinate along one side (its poles x), Z(u) = X(1 - u) and H a step
// from 0 to 1 with no slope or bend at its ends:
//   S(u, v) = V + X(u) (p01 + (p02 - p01) H(v)) e0
//               + (Z(u) p12 + (1 - Z(u)) X(v) p10) e1 + (Z(u) p21 + (1 - Z(u)) X(1 - v) p20) e2.
// Its sides lie on the blends: at u = 0, (s1, s2) is blend 0's section and s0 runs along its
// rail; at v = 0 and v = 1 the same for blends 1 and 2; u = 1 collapses to Q12. Across each
// side its first and second derivatives run along that blend's rail, so it is G2 to the G2
// blends (X ends in three zero poles) and G1 to conics. A conic's weight depends on the angle
// of its faces: its corner is made where the three blends share one weight and the step H
// is not needed. A chamfer's section is the line between its feet, X = (1, 0): the corner is

std::pair<TopoDS_Vertex, double> Attach::contactsCross(const Corner& at, int i, int j)
{
  // Where the contacts of the blends at ends i and j cross on the face they share: the corner of
  // both, with the gap between their sheets' points there.
  const auto& [ki, ei] = at.ends[i];
  const auto& [kj, ej] = at.ends[j];
  const Piece& a       = pieces[ki];
  const Piece& b       = pieces[kj];
  int          si = -1, sj = -1;
  for (int s = 0; s < 2; ++s)
    for (int t = 0; t < 2; ++t)
      if (a.rail.s[s].face.IsSame(b.rail.s[t].face))
      {
        si = s;
        sj = t;
      }
  if (si < 0)
    fail(op,
         unsupported_geometry,
         "two blends at a corner have no face in common",
         {a.rail.edge});
  const auto [a0, a1] = halfOf(a, ei);
  const auto [b0, b1] = halfOf(b, ej);
  const auto at2      = nearest(meetings(between(a.fitted.contact[si], a0, a1),
                                    between(b.fitted.contact[sj], b0, b1),
                                    a.rail.s[si].surface),
                           a.uAt(ei));
  if (!at2)
    fail(op,
         no_solution,
         "the contacts of two blends at a corner do not cross on their face",
         {a.rail.edge});
  const gp_Pnt here = a.sheet.contact(si)->Value(a.sheet.rail(at2->first));
  const double gap  = here.Distance(b.sheet.contact(sj)->Value(b.sheet.rail(at2->second)));
  if (gap > a.rail.resolved(limit::together))
    fail(op,
         no_solution,
         "the contacts of two blends at a corner cross at different places",
         {a.rail.edge});
  TopoDS_Vertex Q;
  builder.MakeVertex(Q, here, std::max(tolerance, gap));
  corner[ki][si][ei] = corner[kj][sj][ej] = Q;
  u[ki][si][ei]                           = at2->first;
  u[kj][sj][ej]                           = at2->second;
  return {Q, gap};
}

TopoDS_Face Attach::lawCorner(const Corner& at, size_t c)
{
  const Piece& first = pieces[at.ends[0].first];
  const gp_Pnt V     = BRep_Tool::Pnt(at.vertex);
  gp_Vec       e[3];
  for (int i = 0; i < 3; ++i)
  {
    const auto& [k, end] = at.ends[i];
    if (!leaving(pieces[k].rail.edge, at.vertex, e[i], end))
      fail(op, no_solution, "a blended edge has no direction at its corner", {pieces[k].rail.edge});
  }
  // The patch's corners on the faces, where each pair's contacts cross.
  double        p[3][3] = {};
  TopoDS_Vertex Q[3][3];
  double crossingGap = 0;
  bool planarContacts = true;
  for (int i = 0; i < 3; ++i)
    for (int j = i + 1; j < 3; ++j)
    {
      const Piece& a          = pieces[at.ends[i].first];
      const auto [vertex, gap] = contactsCross(at, i, j);
      crossingGap             = std::max(crossingGap, gap);
      Q[i][j] = Q[j][i] = vertex;
      const gp_Pnt here = BRep_Tool::Pnt(vertex);
      // Its coordinates along e_i and e_j.
      const gp_Vec q(V, here);
      const double ii = e[i].Dot(e[i]), ij = e[i].Dot(e[j]), jj = e[j].Dot(e[j]),
                   det = ii * jj - ij * ij;
      p[i][j]          = (q.Dot(e[i]) * jj - q.Dot(e[j]) * ij) / det;
      p[j][i]          = (q.Dot(e[j]) * ii - q.Dot(e[i]) * ij) / det;
      if ((q - e[i] * p[i][j] - e[j] * p[j][i]).Magnitude() > a.rail.resolved(limit::together))
      {
        planarContacts = false;
        if (first.rail.law == Law::Chamfer)
          fail(op, no_solution, "a corner of two blends' contacts is off their face", {a.rail.edge});
      }
    }
  // The affine tensor requires every contact crossing in its incident tangent
  // plane. Retain that exact construction whenever its geometric premise holds;
  // otherwise use the carried sheets and the already solved crossings. The
  // allowance is the same positional bound used by the tensor construction.
  if (!planarContacts)
    return tracedCorner(at, c, Q, crossingGap);
  // The law: its coordinate along a side and its weights, and the step.
  std::vector<double> x, w, h, knots;
  std::vector<int>    mults;
  int                 degree;
  if (first.rail.law == Law::G2)
  {
    for (const auto& [along, across] : rho::g2Poles(first.rail.rho))
      x.push_back(1 - along);
    w.assign(x.size(), 1);
    h = {0, 0, 0, 0, 1, 1, 1, 1};
    knots.assign(rho::g2Knots.begin(), rho::g2Knots.end());
    mults.assign(rho::g2Multiplicities.begin(), rho::g2Multiplicities.end());
    degree = 5;
  }
  else if (first.rail.law == Law::Chamfer)
  {
    x      = {1, 0};
    w      = {1, 1};
    h      = {0, 1};
    knots  = {0, 1};
    mults  = {2, 2};
    degree = 1;
  }
  else
  {
    // Each blend's own middle weight, from its sheet: one for all three, and no step.
    double weight = 0;
    for (const auto& [k, end] : at.ends)
    {
      const auto& sheet = pieces[k].fitted.blend;
      if (sheet.IsNull())
        fail(op,
             unsupported_geometry,
             "a conic's corner on an exact sheet, which the blend does not make yet",
             {first.rail.edge});
      const double mine =
        sheet->Weight(sheet->NbUPoles() / 2 + 1, 2) / sheet->Weight(sheet->NbUPoles() / 2 + 1, 1);
      if (weight > 0 && std::abs(mine - weight) > limit::agree)
        fail(op,
             unsupported_geometry,
             "a conic's corner where its faces meet at different angles, which the blend does not "
             "make yet",
             {first.rail.edge});
      weight = mine;
    }
    if (std::abs(p[0][1] - p[0][2]) > limit::agree * first.rail.scale)
      fail(op,
           unsupported_geometry,
           "a conic's corner whose blends reach it unevenly, which the blend does not make yet",
           {first.rail.edge});
    x      = {1, 0, 0};
    w      = {1, weight, 1};
    h      = {0, 0, 1};
    knots  = {0, 1};
    mults  = {3, 3};
    degree = 2;
  }
  const int                  count = int(x.size());
  NCollection_Array2<gp_Pnt> poles(1, count, 1, count);
  NCollection_Array2<double> weights(1, count, 1, count);
  for (int i = 0; i < count; ++i)
    for (int j = 0; j < count; ++j)
    {
      const double z = x[count - 1 - i], beta = p[0][1] + (p[0][2] - p[0][1]) * h[j];
      poles(i + 1, j + 1) =
        V.Translated(e[0] * (x[i] * beta) + e[1] * (z * p[1][2] + (1 - z) * x[j] * p[1][0])
                     + e[2] * (z * p[2][1] + (1 - z) * x[count - 1 - j] * p[2][0]));
      weights(i + 1, j + 1) = w[i] * w[j];
    }
  NCollection_Array1<double> uk(1, int(knots.size()));
  NCollection_Array1<int>    um(1, int(knots.size()));
  for (size_t i = 0; i < knots.size(); ++i)
  {
    uk(int(i) + 1) = knots[i];
    um(int(i) + 1) = mults[i];
  }
  const Handle(Geom_BSplineSurface) patchSurface =
    new Geom_BSplineSurface(poles, weights, uk, uk, um, um, degree, degree);
  // Its corners are the contacts' corners.
  const std::array<std::array<double, 2>, 3> at01 = {
    {{0, 0}, {0, 1}, {1, 0}}}; // Q01, Q02, Q12 in (u, v)
  const TopoDS_Vertex* qs[3] = {&Q[0][1], &Q[0][2], &Q[1][2]};
  for (int i = 0; i < 3; ++i)
    if (patchSurface->Value(at01[i][0], at01[i][1]).Distance(BRep_Tool::Pnt(*qs[i]))
        > first.rail.resolved(limit::together))
      fail(op,
           no_solution,
           "the corner's patch misses a corner of its blends' contacts",
           {first.rail.edge});
  TopoDS_Face face;
  builder.MakeFace(face, patchSurface, tolerance);

  // Its sides: blend 0's end on u = 0 (Q01 to Q02), blend 1's on v = 0 (Q01 to Q12), blend 2's on
  // v = 1 (Q02 to Q12), each laid on its blend's sheet too; and the collapsed u = 1 at Q12.
  const struct
  {
    Handle(Geom_Curve)   curve;
    gp_Pnt2d             from, to;
    const TopoDS_Vertex *start, *stop;
  } sides[3] = {{patchSurface->UIso(0), {0, 0}, {0, 1}, &Q[0][1], &Q[0][2]},
                {patchSurface->VIso(0), {0, 0}, {1, 0}, &Q[0][1], &Q[1][2]},
                {patchSurface->VIso(1), {0, 1}, {1, 1}, &Q[0][2], &Q[1][2]}};

  TopoDS_Edge edges[3];
  for (int i = 0; i < 3; ++i)
  {
    const auto& [k, end] = at.ends[i];
    const Piece&   piece = pieces[k];
    const int      side0 = corner[k][0][end].IsSame(*sides[i].start) ? 0 : 1;
    const gp_Pnt2d startOnSheet =
      piece.sheet.uv(piece.sheet.rail(u[k][side0][end]), piece.sheet.across[side0]);
    // Where the side is the blend's own section at its end (both corners at one place along
    // it, the patch's law there the sheet's), it is the line of the sheet's parameters between
    // them, at the side's own rate; else it is projected onto the sheet.
    Handle(Geom2d_Curve) onSheet = segment2d(
      startOnSheet,
      piece.sheet.uv(piece.sheet.rail(u[k][1 - side0][end]), piece.sheet.across[1 - side0]));
    for (int j = 1; j < 8; ++j)
    {
      const double   t  = j / 8.0;
      const gp_Pnt2d uv = onSheet->Value(t);
      if (piece.sheet.surface->Value(uv.X(), uv.Y()).Distance(sides[i].curve->Value(t)) > tolerance)
      {
        onSheet = GeomProjLib::Curve2d(sides[i].curve, 0, 1, piece.sheet.surface);
        break;
      }
    }
    if (onSheet.IsNull())
      fail(op, no_solution, "a corner's side cannot be laid on its blend", {piece.rail.edge});
    // A chamfer's end runs straight across its sheet only where the sheet is flat there.
    if (piece.rail.law == Law::Chamfer)
      for (const double t : {0.25, 0.5, 0.75})
        if (piece.sheet.surface->Value(onSheet->Value(t).X(), onSheet->Value(t).Y())
              .Distance(sides[i].curve->Value(t))
            > piece.rail.resolved(limit::together))
          fail(op,
               unsupported_geometry,
               "a corner of three chamfers whose ends do not lie flat, which the blend does not "
               "make yet",
               {piece.rail.edge});
    edges[i] = edgeBetween(
      sides[i].curve,
      *sides[i].start,
      0,
      *sides[i].stop,
      1,
      {OnSurface{Handle(Geom2d_Curve)(segment2d(sides[i].from, sides[i].to)),
                 Handle(Geom_Surface)(patchSurface),
                 {}},
       OnSurface{Handle(Geom2d_Curve)(aligned(onSheet, 0, startOnSheet, piece.sheet.surface)),
                 Handle(Geom_Surface)(piece.sheet.surface),
                 {}}},
      tolerance);
    section[k][end] = edges[i];
  }
  const TopoDS_Edge collapsed = degenerateEdge(new Geom2d_Line(gp_Pnt2d(1, 0), gp_Dir2d(0, 1)),
                                               patchSurface,
                                               TopLoc_Location(),
                                               Q[1][2],
                                               0,
                                               1,
                                               tolerance);
  // Counter-clockwise: along v = 0, up the collapsed side, back along v = 1, down u = 0.
  TopoDS_Wire wire;
  builder.MakeWire(wire);
  builder.Add(wire, edges[1]);
  builder.Add(wire, collapsed);
  builder.Add(wire, TopoDS::Edge(edges[2].Reversed()));
  builder.Add(wire, TopoDS::Edge(edges[0].Reversed()));
  wire.Closed(true);
  builder.Add(face, wire);
  // Facing out of the body, as the blends' supports do.
  gp_Vec outward;
  for (const auto& [k, end] : at.ends)
    for (int side = 0; side < 2; ++side)
    {
      const gp_Pnt2d foot = pieces[k].fitted.contact[side]->Value(pieces[k].uAt(end));
      Jet            j;
      if (jetAt(pieces[k].rail.s[side], foot.X(), foot.Y(), j))
        outward += j.n;
    }
  gp_Pnt middle;
  gp_Vec du, dv;
  patchSurface->D1(0.5, 0.5, middle, du, dv);
  if (du.Crossed(dv).Dot(outward) < 0)
    face.Reverse();
  (void)c;
  return face;
}


//! A corner bounded by the carried sheets' actual end traces. The supports are
//! retained as the first- and second-order boundary data of the constrained fill.
TopoDS_Face Attach::tracedCorner(const Corner& at, size_t c,
                               const TopoDS_Vertex (&Q)[3][3], double crossingGap)
{
  const Piece& first = pieces[at.ends[0].first];
  // Each end trace connects the two independently solved contacts in the carried
  // sheet's parameter domain. It is not projected from a proposed corner surface.
  const TopoDS_Vertex starts[3] = {Q[0][1], Q[1][2], Q[0][2]};
  const TopoDS_Vertex stops[3]  = {Q[0][2], Q[0][1], Q[1][2]};
  std::vector<PatchSide> boundary;
  std::vector<TopoDS_Edge> edges;
  const int order = first.rail.law == Law::G2 ? 2 : 1;
  // A spatial trace and a fitted corner are two approximations. Their errors
  // add, so share the existing allowance after reserving the endpoint mismatch.
  const double approximation =
    (first.rail.resolved(limit::together) - crossingGap) / 2;
  for (int i = 0; i < 3; ++i)
  {
    const auto& [k, end] = at.ends[i];
    const Piece& piece = pieces[k];
    const int    side  = corner[k][0][end].IsSame(starts[i]) ? 0 : 1;
    const gp_Pnt2d from =
      piece.sheet.uv(piece.sheet.rail(u[k][side][end]), piece.sheet.across[side]);
    const gp_Pnt2d to =
      piece.sheet.uv(piece.sheet.rail(u[k][1 - side][end]), piece.sheet.across[1 - side]);
    const Handle(Geom2d_Curve) trace = segment2d(from, to);
    BRepBuilderAPI_MakeEdge make(trace, piece.sheet.surface, starts[i], stops[i], 0, 1);
    if (!make.IsDone())
      fail(op, no_solution, "a corner's end trace cannot be bounded on its sheet", {piece.rail.edge});
    TopoDS_Edge edge = make.Edge();
    if (!BRepLib::BuildCurve3d(edge, approximation, GeomAbs_C2, 9))
      fail(op, no_solution, "a corner's carried end trace cannot be approximated", {piece.rail.edge});
    TopoDS_Face support;
    builder.MakeFace(support, piece.sheet.surface, tolerance);
    boundary.push_back({edge, support, order});
    edges.push_back(edge);
    section[k][end] = edge;
  }
  TopoDS_Face made = patchLoop(boundary, approximation, limit::patchDegrees);
  const Handle(Geom_Surface) surface = BRep_Tool::Surface(made);
  // Keep the shared edges and the exact pcurves on the carried sheets. Only their
  // parameter curves on the new corner are obtained from the constrained fit.
  for (const auto& edge : edges)
  {
    double lo, hi;
    const auto curve = BRep_Tool::Curve(edge, lo, hi);
    const auto onPatch = GeomProjLib::Curve2d(curve, lo, hi, surface);
    if (onPatch.IsNull())
      fail(op, no_solution, "a carried end trace cannot be laid on its corner", {first.rail.edge});
    builder.UpdateEdge(edge, onPatch, surface, TopLoc_Location(), tolerance);
    builder.Range(edge, lo, hi);
  }
  BRepBuilderAPI_MakeWire wire;
  wire.Add(edges[0]);
  wire.Add(edges[2]);
  wire.Add(edges[1]);
  if (!wire.IsDone())
    fail(op, no_solution, "the carried end traces do not close", {first.rail.edge});
  BRepBuilderAPI_MakeFace face(surface, wire.Wire(), true);
  if (!face.IsDone())
    fail(op, no_solution, "the corner cannot be bounded by its carried traces", {first.rail.edge});
  TopoDS_Face result = face.Face();
  gp_Vec outward;
  for (const auto& [k, end] : at.ends)
    for (int side = 0; side < 2; ++side)
    {
      const gp_Pnt2d foot = pieces[k].fitted.contact[side]->Value(u[k][side][end]);
      Jet j;
      if (jetAt(pieces[k].rail.s[side], foot.X(), foot.Y(), j))
        outward += j.n;
    }
  double uf, ul, vf, vl;
  surface->Bounds(uf, ul, vf, vl);
  gp_Pnt point;
  gp_Vec du, dv;
  surface->D1((uf + ul) / 2, (vf + vl) / 2, point, du, dv);
  if (du.Crossed(dv).Dot(outward) < 0)
    result.Reverse();
  (void)c;
  return result;
}

// The corner of more than three conic blends between planes, split about a centre C into one
// four-sided face at each place Q_i where the contacts of blends i and i + 1 cross, bounded by a
// half of either blend's end and by two inner curves from the ends' splits M_i to C: the central
// split of an n-sided region (Charrot and Gregory, A pentagonal surface patch for computer aided
// geometric design, 1984; Chiyokura and Kimura, Design of solids with free-form surfaces, 1983).
// Every face is a rational Bezier patch of degree 6, made directly (no fit), and every seam is G1
// by construction, to rounding:
// - The ends. Between planes a sheet is its section's extrusion along its edge, S(s, r) =
//   C(s) + r d. Blend i's end runs from Q_{i-1} to Q_i at r(s) = R(s) / W(s), W the section's
//   denominator and R the cubic with the ends' places and slopes such that the end leaves Q_{i-1}
//   along d_{i-1} and arrives at Q_i along d_{i+1}: a rational cubic on the sheet, exactly. At Q_i
//   both those directions and both sheets' tangent planes are the face's (derived).
// - Across an end. Where a patch's second row of homogeneous poles is its first's moved along d
//   at the same weights, H_j1 = H_j0 + w_j0 l_j (d, 0), its derivative across that side is
//   6 sum(B_j w_j0 l_j) / W d: along the sheet's extrusion, so G1 to the sheet whatever the
//   l_j (derived). At Q_i the two rules give one pole, as either end's tangent there is the
//   other side's direction. (G2 to both sheets is not possible at Q_i: they bend differently
//   along the face there.)
// - Each half end is put in the form whose two first and two last weights are 1 (a Moebius change
//   of parameter, then a cubic factor of numerator and denominator), and split where the halves'
//   speeds agree: the weights along the inner curves are then 1 and the halves join C1 at M_i.
// - Across an inner curve, of tangent T, between the patch a it ends at u = 1 and the patch b it
//   ends at v = 1: D_a + D_b = 2 cos(2 pi / n) t T, the G1 condition of a split with n patches
//   about one point (Loop, Smooth spline surfaces over irregular meshes, 1994; Peters, Patching
//   Catmull-Clark meshes, 2000). Linear in the poles, it holds at the first two coefficients
//   by the rules above (with the offsets next to M_i set so) and is solved for at the others
//   together with the poles that remain, kept fair: the least thin-plate energy of the nets. Its
//   last coefficient makes the curves' tangents at C an affine image of the regular n-gon's.
//! The central split's construction (Attach::splitCorner) in the corner's frame at V: each blend's
//! end as a homogeneous cubic on its sheet, split at M_i into two halves of degree 6; the patches'
//! poles, fixed or solved for under the G1 conditions across the inner curves; the centre where
//! the corner bends least.
struct Attach::CentralSplit
{
  using Homogeneous = std::array<double, 6>; // w (x - V), w s, w r, w: on the corner and a sheet
  static constexpr int degree = 6;

  //! A pole: one of the unknowns, else fixed at `at`.
  struct Ref
  {
    int    var = -1;
    gp_XYZ at;
  };

  //! A term of a fairness stencil or of a condition across an inner curve.
  struct Term
  {
    double coefficient;
    Ref    ref;
  };

  const char*                            op     = nullptr;
  const char*                            notYet = nullptr;
  int                                    count  = 0;
  double                                 scale  = 0;
  gp_XYZ                                 V;
  std::vector<TopoDS_Edge>               edges; // the blends', in turn round it
  std::vector<TopoDS_Vertex>             X;     // X[i]: where blends i and i + 1 cross
  std::vector<gp_Vec>                    d;
  std::vector<Handle(Geom_BSplineCurve)> profile;
  std::vector<std::vector<Homogeneous>>  cubics;
  std::vector<gp_Vec> outwardAt; // out of the face blends i and i + 1 share, at X[i]
  //! Each blend's end: half[i][0] from X[i - 1] to M_i, half[i][1] from X[i] to M_i.
  std::vector<std::array<std::vector<Homogeneous>, 2>> half;
  double                                               beta = 0;
  std::vector<double>                                  legRow, legCol, ell;
  int                                                  curveBase = 0, centre = 0, unknowns = 0;
  std::vector<std::vector<Term>>                       stencils, rows;
  std::vector<math_Vector>                             solution;

  static gp_XYZ point(const Homogeneous& h) { return gp_XYZ(h[0], h[1], h[2]) / h[5]; }

  //! de Casteljau's split of a homogeneous Bezier curve at t.
  static std::pair<std::vector<Homogeneous>, std::vector<Homogeneous>> split(
    std::vector<Homogeneous> b,
    double                   t)
  {
    std::vector<Homogeneous> left, right(b.size());
    for (size_t level = 0; level < b.size(); ++level)
    {
      left.push_back(b[0]);
      right[b.size() - 1 - level] = b[b.size() - 1 - level];
      for (size_t m = 0; m + 1 < b.size() - level; ++m)
        for (int x = 0; x < 6; ++x)
          b[m][x] = (1 - t) * b[m][x] + t * b[m + 1][x];
    }
    return std::make_pair(left, right);
  }

  //! A cubic half in the form of weights 1, 1, ..., 1, 1: its weights' Moebius standard form, then
  //! the cubic factor (1, 2 - w1, 2 - w2, 1) of numerator and denominator. Empty where a weight is
  //! not positive.
  static std::optional<std::vector<Homogeneous>> sextic(std::vector<Homogeneous> b)
  {
    if (!(b[0][5] > 0 && b[3][5] > 0))
      return {};
    const double ratio = std::cbrt(b[0][5] / b[3][5]);
    double       f     = 1 / b[0][5];
    for (auto& pole : b)
    {
      for (auto& x : pole)
        x *= f;
      f *= ratio;
    }
    const double             g[4]  = {1, 2 - b[1][5], 2 - b[2][5], 1};
    constexpr double         c3[4] = {1, 3, 3, 1};
    constexpr double         c6[7] = {1, 6, 15, 20, 15, 6, 1};
    std::vector<Homogeneous> product(degree + 1, Homogeneous{});
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        for (int x = 0; x < 6; ++x)
          product[i + j][x] += c3[i] * c3[j] / c6[i + j] * g[j] * b[i][x];
    for (const auto& pole : product)
      if (!(pole[5] > 0))
        return {};
    return product;
  }

  //! Blend i's end split where its halves' speeds at M_i agree: the first's grows from nothing as
  //! the split moves along, the second's shrinks to nothing.
  void splitEnd(int i)
  {
    const auto halves = [&](double t) {
      const auto [left, right] = split(cubics[i], t);
      return std::make_pair(sextic(left), sextic(right));
    };
    const auto excess = [&](double t) -> std::optional<double> {
      const auto [left, right] = halves(t);
      if (!left || !right)
        return {};
      return (point((*left)[degree]) - point((*left)[degree - 1])).Modulus()
             - (point((*right)[1]) - point((*right)[0])).Modulus();
    };
    const TopoDS_Edge& edge = edges[i];
    const char*        why  = "the end of a blend at a corner cannot be split evenly";
    double             lo = 0.02, hi = 0.98;
    const auto         low = excess(lo), high = excess(hi);
    if (!low || !high || !(*low < 0) || !(*high > 0))
      fail(op, no_solution, why, {edge});
    for (int step = 0; step < 60; ++step)
    {
      const double mid = (lo + hi) / 2;
      const auto   e   = excess(mid);
      if (!e)
        fail(op, no_solution, why, {edge});
      (*e < 0 ? lo : hi) = mid;
    }
    const auto [left, right] = halves((lo + hi) / 2);
    half[i][0]               = *left;
    half[i][1]               = std::vector<Homogeneous>(right->rbegin(), right->rend());
  }

  gp_XYZ splitPoint(int i) const { return point(half[i][0][degree]); }

  //! The offsets across each half end, along its blend's edge: l0 at the crossing (the other
  //! end's first leg), l_i at M_i, and l_i (1 - cos(2 pi / n) / 6) next to M_i, where the
  //! condition across the inner curve then holds at its second coefficient; between, from one to
  //! the other in even steps. l_i is a sixth of the chord from M_i to the splits' mean: the speed
  //! of the inner curve's first leg that of a straight curve to the centre's likely place.
  void legs()
  {
    beta = 2 * std::cos(2 * M_PI / count);
    legRow.resize(count);
    legCol.resize(count);
    ell.resize(count);
    for (int i = 0; i < count; ++i)
    {
      const int    next = (i + 1) % count;
      const gp_Vec alongRow(point(half[i][1][1]) - point(half[i][1][0]));
      const gp_Vec alongCol(point(half[next][0][1]) - point(half[next][0][0]));
      legRow[i] = alongRow.Dot(d[next]);
      legCol[i] = alongCol.Dot(d[i]);
      if (alongRow.Crossed(d[next]).Magnitude() > tol::fineAngle * alongRow.Magnitude()
          || alongCol.Crossed(d[i]).Magnitude() > tol::fineAngle * alongCol.Magnitude())
        fail(op, unsupported_geometry, notYet, {edges[0]});
    }
    gp_XYZ centred;
    for (int i = 0; i < count; ++i)
      centred += splitPoint(i) / count;
    for (int i = 0; i < count; ++i)
    {
      // Both patches beside blend i's end lie on one side of it: towards the corner.
      if (!(legRow[(i + count - 1) % count] * legCol[i] > 0))
        fail(op, no_solution, "the ends of the blends at a corner turn away from it", {edges[i]});
      ell[i] = (legCol[i] > 0 ? 1 : -1) * (centred - splitPoint(i)).Modulus() / degree;
    }
  }

  double offset(double l0, double l, int j) const
  {
    const double inner = l * (1 - beta / (2 * degree));
    return j <= 1            ? l0
           : j == degree     ? l
           : j == degree - 1 ? inner
                             : l0 + (inner - l0) * (j - 1) / 4;
  }

  //! Patch q at X[q]: u along blend q + 1's end to M_{q+1} (its side v = 0), v along blend q's to
  //! M_q (u = 0); the inner curve q + 1 is its side u = 1, curve q its side v = 1. A pole is fixed
  //! or one of the unknowns: per patch its 16 inner ones, per curve 4, and the centre.
  Ref onCurve(int i, int t) const
  {
    if (t <= 1)
      return {-1, splitPoint(i) + d[i].XYZ() * (t * ell[i])};
    return {t == degree ? centre : curveBase + 4 * i + t - 2, {}};
  }

  Ref pole(int q, int j, int k) const
  {
    const int next = (q + 1) % count;
    if (k == degree)
      return onCurve(q, j);
    if (j == degree)
      return onCurve(next, k);
    const auto& A = half[next][0];
    const auto& B = half[q][1];
    if (k == 0)
      return {-1, point(A[j])};
    if (j == 0)
      return {-1, point(B[k])};
    if (k == 1)
      return {-1, point(A[j]) + d[next].XYZ() * offset(legRow[q], ell[next], j)};
    if (j == 1)
      return {-1, point(B[k]) + d[q].XYZ() * offset(legCol[q], ell[q], k)};
    return {16 * q + 4 * (j - 2) + k - 2, {}};
  }

  double weight(int q, int j, int k) const
  {
    if (k <= 1 && j < degree)
      return half[(q + 1) % count][0][j][5];
    if (j <= 1 && k < degree)
      return half[q][1][k][5];
    return 1.0;
  }

  //! The fairness, the thin-plate energy of the nets: their second differences along u and v and
  //! their twists. The conditions across the inner curves, coefficient by coefficient.
  void equations()
  {
    curveBase = 16 * count;
    centre    = 20 * count;
    unknowns  = 20 * count + 1;
    for (int q = 0; q < count; ++q)
      for (int j = 0; j <= degree; ++j)
        for (int k = 0; k <= degree; ++k)
        {
          if (j > 0 && j < degree)
            stencils.push_back(
              {{1, pole(q, j - 1, k)}, {-2, pole(q, j, k)}, {1, pole(q, j + 1, k)}});
          if (k > 0 && k < degree)
            stencils.push_back(
              {{1, pole(q, j, k - 1)}, {-2, pole(q, j, k)}, {1, pole(q, j, k + 1)}});
          if (j < degree && k < degree)
            stencils.push_back({{M_SQRT2, pole(q, j + 1, k + 1)},
                                {-M_SQRT2, pole(q, j + 1, k)},
                                {-M_SQRT2, pole(q, j, k + 1)},
                                {M_SQRT2, pole(q, j, k)}});
        }
    for (int i = 0; i < count; ++i)
    {
      const int a = (i + count - 1) % count, b = i;
      for (int t = 0; t <= degree; ++t)
        rows.push_back({{2 * degree - beta * t, onCurve(i, t)},
                        {beta * t, t > 0 ? onCurve(i, t - 1) : Ref{}},
                        {-double(degree), pole(a, degree - 1, t)},
                        {-double(degree), pole(b, t, degree - 1)}});
    }
  }

  //! The unknowns as an affine function of the centre c, and the centre where the nets' energy is
  //! least.
  gp_XYZ solve()
  {
    // One system for the three coordinates, [H A'; A 0] [x; y] = [-g; b], with one more row that
    // holds the centre at a point c; its constraints are dependent (the centre's take two of them),
    // so it is solved by singular values. Its solution is x0 + g c, g the same for every
    // coordinate.
    std::vector<std::vector<Term>> conditions;
    std::copy_if(rows.begin(), rows.end(), std::back_inserter(conditions), [](const auto& row) {
      return std::any_of(row.begin(), row.end(), [](const Term& term) {
        return term.ref.var >= 0;
      });
    });
    const int                size = unknowns + int(conditions.size()) + 1;
    math_Matrix              K(1, size, 1, size, 0.0);
    std::vector<math_Vector> rhs(4, math_Vector(1, size, 0.0));
    const auto               fixedPart = [&](int row, double c, const gp_XYZ& at) {
      for (int x = 0; x < 3; ++x)
        rhs[x](row) -= c * at.Coord(x + 1);
    };
    for (const auto& stencil : stencils)
      for (const auto& p : stencil)
        if (p.ref.var >= 0)
          for (const auto& q : stencil)
            if (q.ref.var >= 0)
              K(p.ref.var + 1, q.ref.var + 1) += p.coefficient * q.coefficient;
            else
              fixedPart(p.ref.var + 1, p.coefficient * q.coefficient, q.ref.at);
    for (size_t r = 0; r < conditions.size(); ++r)
    {
      const int row = unknowns + int(r) + 1;
      for (const auto& term : conditions[r])
        if (term.ref.var >= 0)
        {
          K(row, term.ref.var + 1) += term.coefficient;
          K(term.ref.var + 1, row) += term.coefficient;
        }
        else
          fixedPart(row, term.coefficient, term.ref.at);
    }
    K(size, centre + 1) = K(centre + 1, size) = 1;
    rhs[3](size)                              = 1;
    math_SVD svd(K);
    solution.assign(4, math_Vector(1, size));
    for (int x = 0; x < 4; ++x)
      svd.Solve(rhs[x], solution[x], 1e-12);
    const math_Vector& g = solution[3];
    // Where the thin-plate energy of the nets is least: x' H x - 2 x' r is least along x0 + g c.
    gp_XYZ fair;
    {
      double gHg = 0;
      gp_XYZ gR;
      for (int a = 1; a <= unknowns; ++a)
      {
        double Hg = 0;
        gp_XYZ Hx;
        for (int b = 1; b <= unknowns; ++b)
        {
          Hg += K(a, b) * g(b);
          Hx += gp_XYZ(solution[0](b), solution[1](b), solution[2](b)) * K(a, b);
        }
        gHg += g(a) * Hg;
        gR += (gp_XYZ(rhs[0](a), rhs[1](a), rhs[2](a)) - Hx) * g(a);
      }
      fair = gR / gHg;
    }
    return fair;
  }

  gp_XYZ valueFor(const Ref& ref, const gp_XYZ& c) const
  {
    const math_Vector& g = solution[3];
    if (ref.var < 0)
      return ref.at;
    const int x = ref.var + 1;
    return gp_XYZ(solution[0](x), solution[1](x), solution[2](x)) + c * g(x);
  }

  // Per point the homogeneous point and its derivatives (u, v, uu, uv, vv), split into the part
  // that does not move with the centre and the factor that does, and the weight's; summed over
  // u first, then v.
  struct Node
  {
    gp_XYZ h[6];
    double gamma[6], w[6], area;
  };

  struct Homogeneous5
  {
    gp_XYZ h;
    double gamma = 0, w = 0;
  };

  //! The surface's points and derivatives at the Gauss points of every patch, as affine functions
  //! of the centre.
  std::vector<Node> gaussNodes() const
  {
    const math_Vector& g = solution[3];
    // A composite Gauss rule, 4 x 4 squares of 8 x 8 points per patch. The energy's integrand
    // varies sharply where a patch is narrow (one Gauss rule of 10 x 10 points misses it by a
    // quarter); on the uneven family the centre this rule finds lies within 0.7 % of the corner's
    // size of the one 8 x 8 squares of 6 x 6 points find (measured).
    constexpr int                           cells = 4, nodes = 8, points = cells * nodes;
    const math_ComputeGaussPointsAndWeights gauss(nodes);
    const math_Vector                       gx = gauss.Points(), gw = gauss.Weights();
    // Bernstein values and derivatives at the points on [0, 1], and their weights.
    std::vector<std::array<std::array<double, degree + 1>, 3>> basis(points);
    std::vector<double>                                        share(points);
    for (int a = 0; a < points; ++a)
    {
      const double t = (a / nodes + (gx(gx.Lower() + a % nodes) + 1) / 2) / cells;
      share[a]       = gw(gw.Lower() + a % nodes) / (2 * cells);
      const auto b   = [&](int m, int i) {
        if (i < 0 || i > m)
          return 0.0;
        double c = 1;
        for (int k = 1; k <= i; ++k)
          c = c * (m - i + k) / k;
        return c * std::pow(t, i) * std::pow(1 - t, m - i);
      };
      for (int j = 0; j <= degree; ++j)
      {
        basis[a][0][j] = b(degree, j);
        basis[a][1][j] = degree * (b(degree - 1, j - 1) - b(degree - 1, j));
        basis[a][2][j] = degree * (degree - 1)
                         * (b(degree - 2, j - 2) - 2 * b(degree - 2, j - 1) + b(degree - 2, j));
      }
    }

    std::vector<Node> at(size_t(count) * points * points);
    constexpr int     du[6] = {0, 1, 0, 2, 1, 0}, dv[6] = {0, 0, 1, 0, 1, 2};
    for (int q = 0; q < count; ++q)
    {
      Homogeneous5 net[degree + 1][degree + 1];
      for (int j = 0; j <= degree; ++j)
        for (int k = 0; k <= degree; ++k)
        {
          const Ref    ref = pole(q, j, k);
          const double w   = weight(q, j, k);
          net[j][k] = {valueFor(ref, gp_XYZ()) * w, ref.var < 0 ? 0.0 : g(ref.var + 1) * w, w};
        }
      std::vector<std::array<std::array<Homogeneous5, degree + 1>, 3>> row(points);
      for (int a = 0; a < points; ++a)
        for (int order = 0; order < 3; ++order)
          for (int k = 0; k <= degree; ++k)
          {
            Homogeneous5& sum = row[a][order][k];
            for (int j = 0; j <= degree; ++j)
            {
              const double c = basis[a][order][j];
              sum.h += net[j][k].h * c;
              sum.gamma += net[j][k].gamma * c;
              sum.w += net[j][k].w * c;
            }
          }
      for (int a = 0; a < points; ++a)
        for (int b = 0; b < points; ++b)
        {
          Node& node = at[(size_t(q) * points + a) * points + b];
          for (int m = 0; m < 6; ++m)
            for (int k = 0; k <= degree; ++k)
            {
              const double        c   = basis[b][dv[m]][k];
              const Homogeneous5& sum = row[a][du[m]][k];
              node.h[m] += sum.h * c;
              node.gamma[m] += sum.gamma * c;
              node.w[m] += sum.w * c;
            }
          node.area = share[a] * share[b];
        }
    }
    return at;
  }

  // The centre itself where the corner bends least: the least integral of k1^2 + k2^2 over the
  // patches, the bending energy of the surface (not of its nets). The nets' energy is a stand-in
  // for it only where the patches are parametrised evenly; uneven corners have patches of very
  // different shapes, and it pulls the centre towards the long ones. Every other pole stays the
  // nets' fairest for that centre. The energy is a sum of squares of the shape operator's entries
  // in an orthonormal frame at Gauss points, minimised from the nets' fairest centre by
  // Levenberg-Marquardt steps until no step lowers it.
  gp_XYZ centreOf(const gp_XYZ& fair) const
  {
    const std::vector<Node> at        = gaussNodes();
    const auto              residuals = [&](const gp_XYZ& c, std::vector<double>& r) {
      r.resize(3 * at.size());
      for (size_t i = 0; i < at.size(); ++i)
      {
        const Node& node = at[i];
        gp_XYZ      H[6];
        for (int m = 0; m < 6; ++m)
          H[m] = node.h[m] + c * node.gamma[m];
        const double* W   = node.w;
        const gp_XYZ  S   = H[0] / W[0];
        const gp_XYZ  Su  = (H[1] - S * W[1]) / W[0];
        const gp_XYZ  Sv  = (H[2] - S * W[2]) / W[0];
        const gp_XYZ  Suu = (H[3] - Su * (2 * W[1]) - S * W[3]) / W[0];
        const gp_XYZ  Suv = (H[4] - Su * W[2] - Sv * W[1] - S * W[4]) / W[0];
        const gp_XYZ  Svv = (H[5] - Sv * (2 * W[2]) - S * W[5]) / W[0];
        const gp_XYZ  N   = Su.Crossed(Sv);
        const double  J   = N.Modulus();
        const gp_XYZ  n   = N / J;
        const double  L = Suu.Dot(n), M = Suv.Dot(n), Nn = Svv.Dot(n);
        // e1 = alpha Su, e2 = beta Su + gamma Sv, orthonormal.
        const double alpha = 1 / Su.Modulus();
        const double along = Sv.Dot(Su) * alpha;
        const double gamma = 1 / std::sqrt(std::max(Sv.SquareModulus() - along * along, 0.0));
        const double beta  = -gamma * along * alpha;
        const double s     = std::sqrt(node.area * J);
        r[3 * i]           = s * alpha * alpha * L;
        r[3 * i + 1]       = s * M_SQRT2 * alpha * (beta * L + gamma * M);
        r[3 * i + 2]       = s * (beta * beta * L + 2 * beta * gamma * M + gamma * gamma * Nn);
      }
    };
    const auto energy = [](const std::vector<double>& r) {
      double e = 0;
      for (const double x : r)
        e += x * x;
      return e;
    };
    // The gradient 2 J' r, J the residuals' Jacobian by differences.
    const double h        = 1e-7 * scale;
    const auto   gradient = [&](const gp_XYZ& c, const std::vector<double>& r, math_Matrix* JtJ) {
      std::vector<double>                plus;
      math_Vector                        grad(1, 3, 0.0);
      std::vector<std::array<double, 3>> jac(r.size());
      for (int x = 0; x < 3; ++x)
      {
        gp_XYZ step;
        step.SetCoord(x + 1, h);
        residuals(c + step, plus);
        for (size_t i = 0; i < r.size(); ++i)
        {
          jac[i][x] = (plus[i] - r[i]) / h;
          grad(x + 1) += 2 * jac[i][x] * r[i];
        }
      }
      if (JtJ)
        for (size_t i = 0; i < r.size(); ++i)
          for (int x = 0; x < 3; ++x)
            for (int y = 0; y < 3; ++y)
              (*JtJ)(x + 1, y + 1) += 2 * jac[i][x] * jac[i][y];
      return grad;
    };
    // Quasi-Newton steps (BFGS) from the Gauss-Newton matrix 2 J'J: the residuals are not small,
    // so that matrix alone misses the energy's curvature and its steps converge slowly; the
    // secant updates supply what it misses [Dennis and Schnabel 1983, 9.2]. A step is halved until
    // it lowers the energy; none does by more than its rounding where the energy is least.
    gp_XYZ              c = fair;
    std::vector<double> r, trial;
    residuals(c, r);
    double      e = energy(r);
    math_Matrix B(1, 3, 1, 3, 0.0);
    math_Vector grad = gradient(c, r, &B);
    for (;;)
    {
      math_Vector step(1, 3);
      math_SVD(B).Solve(-grad, step, 1e-14);
      gp_XYZ move(step(1), step(2), step(3));
      double f;
      for (;;)
      {
        if (!(move.Modulus() >= 1e-12 * scale))
          return c;
        residuals(c + move, trial);
        f = energy(trial);
        if (f < e)
          break;
        move /= 2;
      }
      const bool settled = e - f <= 1e-12 * e;
      c += move;
      r.swap(trial);
      e = f;
      if (settled)
        return c;
      const math_Vector next = gradient(c, r, nullptr);
      const math_Vector y    = next - grad;
      math_Vector       sv(1, 3);
      sv(1)           = move.X();
      sv(2)           = move.Y();
      sv(3)           = move.Z();
      const double sy = sv * y;
      if (sy > 0)
      {
        const math_Vector Bs  = B * sv;
        const double      sBs = sv * Bs;
        for (int x = 1; x <= 3; ++x)
          for (int z = 1; z <= 3; ++z)
            B(x, z) += y(x) * y(z) / sy - Bs(x) * Bs(z) / sBs;
      }
      grad = next;
    }
  }

  //! The conditions hold: the first two by construction, the others as solved.
  void checkJoins(const gp_XYZ& centrePoint) const
  {
    const auto   value = [&](const Ref& ref) { return valueFor(ref, centrePoint); };
    const double exact = 1e-9 * scale;
    for (const auto& row : rows)
    {
      gp_XYZ sum;
      for (const auto& term : row)
        sum += value(term.ref) * term.coefficient;
      if (sum.Modulus() > exact)
        fail(op, no_solution, "the corner's split does not join smoothly", {edges[0]});
    }
  }
};

void Attach::splitCornerEnds(const Corner& at, CentralSplit& cs)
{
  using Homogeneous       = CentralSplit::Homogeneous;
  const int     count     = cs.count;
  const gp_XYZ& V         = cs.V;
  const char*   notYet    = cs.notYet;
  const auto&   X         = cs.X;
  const auto&   d         = cs.d;
  const auto&   profile   = cs.profile;
  auto&         cubics    = cs.cubics;
  auto&         outwardAt = cs.outwardAt;
  for (int i = 0; i < count; ++i)
  {
    const auto& [k, end] = at.ends[i];
    const Piece& piece   = pieces[k];
    const int    prev = (i + count - 1) % count, next = (i + 1) % count;
    const int    a = corner[k][0][end].IsSame(X[prev]) ? 0 : 1;
    if (!corner[k][a][end].IsSame(X[prev]) || !corner[k][1 - a][end].IsSame(X[i]))
      fail(op,
           no_solution,
           "the contacts of the blends at a corner do not cross in turn round it",
           {piece.rail.edge});
    {
      const gp_Pnt2d foot = piece.fitted.contact[1 - a]->Value(u[k][1 - a][end]);
      Jet            j;
      if (!jetAt(piece.rail.s[1 - a], foot.X(), foot.Y(), j))
        fail(op, unsupported_geometry, notYet, {piece.rail.edge});
      outwardAt[i] = j.n;
    }
    const gp_Pnt2d from = piece.sheet.uv(piece.sheet.rail(u[k][a][end]), piece.sheet.across[a]);
    const gp_Pnt2d to =
      piece.sheet.uv(piece.sheet.rail(u[k][1 - a][end]), piece.sheet.across[1 - a]);
    const auto& conic = profile[i];
    const bool  reversed =
      std::abs(from.X() - conic->LastParameter()) < std::abs(from.X() - conic->FirstParameter());
    gp_XYZ p[3];
    double w[3];
    for (int m = 0; m < 3; ++m)
    {
      p[m] = conic->Pole(reversed ? 3 - m : m + 1).XYZ() - V;
      w[m] = conic->Weight(reversed ? 3 - m : m + 1);
    }
    // The section's slopes at its ends; the rail's slopes there that turn them along the
    // neighbours' edges.
    const gp_Vec start(gp_Vec(p[1] - p[0]) * (2 * w[1] / w[0]));
    const gp_Vec stop(gp_Vec(p[2] - p[1]) * (2 * w[1] / w[2]));
    const auto   slope = [&](const gp_Vec& c, const gp_Vec& along) {
      const gp_Vec normal = d[i].Crossed(along);
      const double rho    = -c.Crossed(along).Dot(normal) / normal.SquareMagnitude();
      const gp_Vec turned = c + d[i] * rho;
      if (turned.Crossed(along).Magnitude() > tol::fineAngle * turned.Magnitude())
        fail(op, unsupported_geometry, notYet, {piece.rail.edge});
      return rho;
    };
    const double             rhoA = slope(start, d[prev]), rhoB = slope(stop, d[next]);
    const double             R0 = from.Y() * w[0], R3 = to.Y() * w[2];
    const double             R[4] = {R0,
                                     R0 + (rhoA * w[0] + from.Y() * 2 * (w[1] - w[0])) / 3,
                                     R3 - (rhoB * w[2] + to.Y() * 2 * (w[2] - w[1])) / 3,
                                     R3};
    const double             s0 = from.X(), s1 = to.X();
    const gp_XYZ             hp[4] = {p[0] * w[0],
                                      (p[0] * w[0] + p[1] * (2 * w[1])) / 3,
                                      (p[1] * (2 * w[1]) + p[2] * w[2]) / 3,
                                      p[2] * w[2]};
    const double             hw[4] = {w[0], (w[0] + 2 * w[1]) / 3, (2 * w[1] + w[2]) / 3, w[2]};
    const double             hs[4] = {s0 * w[0],
                                      (2 * s0 * w[1] + s1 * w[0]) / 3,
                                      (s0 * w[2] + 2 * s1 * w[1]) / 3,
                                      s1 * w[2]};
    std::vector<Homogeneous> cubic(4);
    for (int m = 0; m < 4; ++m)
    {
      const gp_XYZ x = hp[m] + d[i].XYZ() * R[m];
      cubic[m]       = {x.X(), x.Y(), x.Z(), hs[m], R[m], hw[m]};
    }
    // Its ends at the crossings' vertices, which both sheets meet to within their gap.
    for (const auto& [m, vertex] : {std::make_pair(0, X[prev]), std::make_pair(3, X[i])})
    {
      const gp_XYZ x = (BRep_Tool::Pnt(vertex).XYZ() - V) * cubic[m][5];
      cubic[m][0]    = x.X();
      cubic[m][1]    = x.Y();
      cubic[m][2]    = x.Z();
    }
    cubics[i] = cubic;
  }
}

void Attach::splitCornerFaces(const Corner&       at,
                              size_t              c,
                              const CentralSplit& cs,
                              const gp_XYZ&       centrePoint)
{
  constexpr int degree    = CentralSplit::degree;
  const Piece&  first     = pieces[at.ends[0].first];
  const int     count     = cs.count;
  const gp_XYZ& V         = cs.V;
  const auto&   X         = cs.X;
  const auto&   half      = cs.half;
  const auto&   outwardAt = cs.outwardAt;
  const auto    point     = CentralSplit::point;
  const auto    pole      = [&](int q, int j, int k) { return cs.pole(q, j, k); };
  const auto    weight    = [&](int q, int j, int k) { return cs.weight(q, j, k); };
  const auto    onCurve   = [&](int i, int t) { return cs.onCurve(i, t); };
  const auto    value = [&](const CentralSplit::Ref& ref) { return cs.valueFor(ref, centrePoint); };
  // The patches, the ends' halves and the inner curves, all on shared edges.
  NCollection_Array1<double> knots(1, 2);
  NCollection_Array1<int>    mults(1, 2);
  knots(1) = 0;
  knots(2) = 1;
  mults(1) = mults(2) = degree + 1;
  std::vector<Handle(Geom_BSplineSurface)> surface(count);
  for (int q = 0; q < count; ++q)
  {
    NCollection_Array2<gp_Pnt> poles(1, degree + 1, 1, degree + 1);
    NCollection_Array2<double> weights(1, degree + 1, 1, degree + 1);
    for (int j = 0; j <= degree; ++j)
      for (int k = 0; k <= degree; ++k)
      {
        poles(j + 1, k + 1)   = gp_Pnt(value(pole(q, j, k)) + V);
        weights(j + 1, k + 1) = weight(q, j, k);
      }
    surface[q] =
      new Geom_BSplineSurface(poles, weights, knots, knots, mults, mults, degree, degree);
  }
  const auto curve3d = [&](const std::vector<gp_XYZ>& at, const std::vector<double>& w) {
    NCollection_Array1<gp_Pnt> poles(1, degree + 1);
    NCollection_Array1<double> weights(1, degree + 1);
    for (int m = 0; m <= degree; ++m)
    {
      poles(m + 1)   = gp_Pnt(at[m] + V);
      weights(m + 1) = w[m];
    }
    return Handle(Geom_Curve)(new Geom_BSplineCurve(poles, weights, knots, mults, degree));
  };
  std::vector<TopoDS_Vertex> M(count);
  for (int i = 0; i < count; ++i)
    builder.MakeVertex(M[i], gp_Pnt(point(half[i][0][degree]) + V), tolerance);
  TopoDS_Vertex C;
  builder.MakeVertex(C, gp_Pnt(value(onCurve(0, degree)) + V), tolerance);
  std::vector<std::array<TopoDS_Edge, 2>> ends(count);  // [blend][half]
  std::vector<TopoDS_Edge>                inner(count); // inner curve i, from M_i to C
  for (int i = 0; i < count; ++i)
  {
    const int    prev  = (i + count - 1) % count;
    const Sheet& sheet = pieces[at.ends[i].first].sheet;
    for (int h = 0; h < 2; ++h)
    {
      std::vector<gp_XYZ>          at3(degree + 1);
      std::vector<double>          w(degree + 1);
      NCollection_Array1<gp_Pnt2d> uv(1, degree + 1);
      NCollection_Array1<double>   uvWeights(1, degree + 1);
      for (int m = 0; m <= degree; ++m)
      {
        const auto& pole = half[i][h][m];
        at3[m]           = point(pole);
        w[m]             = pole[5];
        uv(m + 1)        = gp_Pnt2d(pole[3] / pole[5], pole[4] / pole[5]);
        uvWeights(m + 1) = pole[5];
      }
      // On its patch: patch prev's side v = 0 for the first half, patch i's side u = 0 for the
      // second.
      const Handle(Geom2d_Curve) onPatch =
        segment2d(gp_Pnt2d(0, 0), h == 0 ? gp_Pnt2d(1, 0) : gp_Pnt2d(0, 1));
      ends[i][h] =
        edgeBetween(curve3d(at3, w),
                    h == 0 ? X[prev] : X[i],
                    0,
                    M[i],
                    1,
                    {OnSurface{onPatch, Handle(Geom_Surface)(surface[h == 0 ? prev : i]), {}},
                     OnSurface{new Geom2d_BSplineCurve(uv, uvWeights, knots, mults, degree),
                               sheet.surface,
                               {}}},
                    tolerance);
    }
    std::vector<gp_XYZ> at3(degree + 1);
    for (int t = 0; t <= degree; ++t)
      at3[t] = value(onCurve(i, t));
    inner[i] = edgeBetween(
      curve3d(at3, std::vector<double>(degree + 1, 1.0)),
      M[i],
      0,
      C,
      1,
      {OnSurface{segment2d(gp_Pnt2d(1, 0), gp_Pnt2d(1, 1)),
                 Handle(Geom_Surface)(surface[prev]),
                 {}},
       OnSurface{segment2d(gp_Pnt2d(0, 1), gp_Pnt2d(1, 1)), Handle(Geom_Surface)(surface[i]), {}}},
      tolerance);
  }
  std::vector<TopoDS_Face> faces(count);
  for (int q = 0; q < count; ++q)
  {
    const int next = (q + 1) % count;
    builder.MakeFace(faces[q], surface[q], tolerance);
    // Counter-clockwise: along v = 0, up u = 1, back along v = 1, down u = 0.
    TopoDS_Wire wire;
    builder.MakeWire(wire);
    builder.Add(wire, ends[next][0]);
    builder.Add(wire, inner[next]);
    builder.Add(wire, TopoDS::Edge(inner[q].Reversed()));
    builder.Add(wire, TopoDS::Edge(ends[q][1].Reversed()));
    wire.Closed(true);
    builder.Add(faces[q], wire);
    // Facing out of the body, as the face it is tangent to at its crossing does: there its sides
    // leave along both blends' edges, in that face.
    gp_Pnt at;
    gp_Vec du, dv;
    surface[q]->D1(0, 0, at, du, dv);
    if (du.Crossed(dv).Dot(outwardAt[q]) < 0)
      faces[q].Reverse();
    if (faces[q].Orientation() != faces[0].Orientation())
      fail(op, no_solution, "the faces of a corner's split do not face one way", {first.rail.edge});
  }
  // Each blend's end: its two halves from side 0's corner to side 1's.
  for (int i = 0; i < count; ++i)
  {
    const auto& [k, end] = at.ends[i];
    const int  prev      = (i + count - 1) % count;
    const bool fromPrev  = corner[k][0][end].IsSame(X[prev]);
    chain[k][end]        = fromPrev ? std::vector<TopoDS_Edge>{ends[i][0], ends[i][1]}
                                    : std::vector<TopoDS_Edge>{ends[i][1], ends[i][0]};
    section[k][end]      = chain[k][end].front();
    partEnds.push_back({ends[i][0], k, faces[prev]});
    partEnds.push_back({ends[i][1], k, faces[i]});
    partInner.push_back({inner[i], faces[prev], faces[i]});
  }
  parts[c] = faces;
}

void Attach::splitCorner(const Corner& at, size_t c)
{
  const Piece& first = pieces[at.ends[0].first];
  CentralSplit cs;
  cs.op           = op;
  cs.count        = int(at.ends.size());
  cs.scale        = first.rail.scale;
  cs.V            = BRep_Tool::Pnt(at.vertex).XYZ();
  cs.notYet       = "where more than three conic blends meet, the blend makes the corner "
                    "only between planes so far";
  const int count = cs.count;
  cs.X.resize(count);
  for (int i = 0; i < count; ++i)
    cs.X[i] = contactsCross(at, i, (i + 1) % count).first;
  cs.d.resize(count);
  cs.profile.resize(count);
  for (int i = 0; i < count; ++i)
  {
    const Piece& piece     = pieces[at.ends[i].first];
    const auto   extrusion = Handle(Geom_SurfaceOfLinearExtrusion)::DownCast(piece.sheet.surface);
    if (!extrusion.IsNull())
      cs.profile[i] = Handle(Geom_BSplineCurve)::DownCast(extrusion->BasisCurve());
    if (!piece.sheet.swapped || cs.profile[i].IsNull() || cs.profile[i]->Degree() != 2
        || cs.profile[i]->NbPoles() != 3)
      fail(op, unsupported_geometry, cs.notYet, {piece.rail.edge});
    cs.d[i] = gp_Vec(extrusion->Direction());
  }
  for (int i = 0; i < count; ++i)
    cs.edges.push_back(pieces[at.ends[i].first].rail.edge);
  cs.half.resize(count);
  cs.cubics.resize(count);
  cs.outwardAt.resize(count);
  splitCornerEnds(at, cs);
  for (int i = 0; i < count; ++i)
    cs.splitEnd(i);
  cs.legs();
  cs.equations();
  const gp_XYZ fair = cs.solve();
  // The centre itself where the corner bends least (centreOf); every other pole stays the nets'
  // fairest for that centre.
  const gp_XYZ centrePoint = cs.centreOf(fair);
  cs.checkJoins(centrePoint);
  splitCornerFaces(at, c, cs, centrePoint);
}

} // namespace BRepExt::blending::detail
