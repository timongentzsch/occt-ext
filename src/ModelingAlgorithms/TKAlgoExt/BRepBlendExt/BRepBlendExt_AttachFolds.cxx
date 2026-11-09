// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Blend: The loops where a contact on a cap folds over itself, and their patches.

#include <BRepBlendExt_Attach.hxx>

namespace BRepExt::blending::detail
{

//! The loop where the contact on the cap folds over itself, between piece `a`'s end `ea` and
//! piece `b`'s end `eb`: both end at the place where the contact crosses itself (the one ball
//! there), each in its own section of that ball down to the other contact, and the other contact's
//! run between them is an edge of its face, which the patch over the loop meets. Within one edge
//! the run is the fit's, on either part's sheet alike; across a joint of a chain it runs on each
//! piece's sheet to the joint's place, where the edges at the joint's vertex between the faces it
//! runs on end, cut back to it as at any joint.
void Attach::fold(int a, int ea, int b, int eb)
{
  if (pieces[a].joint[ea].crease)
    return crossingFold(a, ea, b, eb);
  if (pieces[a].rail.law == Law::Chamfer)
    return creaseFold(a, ea, b, eb);
  const Piece& one = pieces[a];
  const Piece& two = pieces[b];
  const int    cap = one.joint[ea].side, capB = two.joint[eb].side;
  const int    other = 1 - cap, otherB = 1 - capB;
  const double ua = one.uAt(ea), ub = two.uAt(eb);
  const gp_Pnt X   = one.sheet.contact(cap)->Value(one.sheet.rail(ua));
  const double gap = X.Distance(two.sheet.contact(capB)->Value(two.sheet.rail(ub)));
  if (gap > one.rail.resolved(limit::together))
    fail(op,
         no_solution,
         "the contact on the cap does not cross itself where the blend folds",
         {one.rail.edge});
  builder.MakeVertex(corner[a][cap][ea], X, std::max(tolerance, gap));
  corner[b][capB][eb] = corner[a][cap][ea];
  builder.MakeVertex(corner[a][other][ea],
                     one.sheet.contact(other)->Value(one.sheet.rail(ua)),
                     tolerance);
  builder.MakeVertex(corner[b][otherB][eb],
                     two.sheet.contact(otherB)->Value(two.sheet.rail(ub)),
                     tolerance);
  for (int side = 0; side < 2; ++side)
  {
    u[a][side][ea] = ua;
    u[b][side][eb] = ub;
  }
  for (const auto& [k, end, at] : {std::tuple{a, ea, ua}, std::tuple{b, eb, ub}})
  {
    const Sheet& sheet = pieces[k].sheet;
    TopoDS_Edge  x     = isoSection(sheet, sheet.rail(at), corner[k][0][end], corner[k][1][end]);
    double       lo, hi;
    BRep_Tool::Range(x, lo, hi);
    builder.UpdateEdge(x,
                       sheet.sectionUV(sheet.rail(at)),
                       sheet.surface,
                       TopLoc_Location(),
                       tolerance);
    builder.Range(x, lo, hi);
    section[k][end] = x;
  }
  Folded     f{a, ea, b, eb, {}, {}};
  const auto runOn = [&](const Piece&         p,
                         int                  side,
                         const TopoDS_Vertex& from,
                         double               uFrom,
                         const TopoDS_Vertex& to,
                         double               uTo) {
    const Sheet&       sheet       = p.sheet;
    const TopoDS_Face& face        = p.rail.s[side].face;
    const auto [support, location] = surfaceOf(face);
    const TopoDS_Edge edge =
      edgeBetween(sheet.contact(side),
                  from,
                  sheet.rail(uFrom),
                  to,
                  sheet.rail(uTo),
                  {{alongRail(p.fitted.contact[side], sheet), support, location}},
                  tolerance);
    editOf(edits, face).bridges.push_back(edge);
    f.run.push_back({edge, &p.rail.s[side], TopExp::FirstVertex(edge).IsSame(from)});
  };
  if (one.rail.edge.IsSame(two.rail.edge))
    runOn(one, other, corner[a][other][ea], ua, corner[b][otherB][eb], ub);
  else
  {
    TopoDS_Vertex ends[2], endsB[2];
    TopExp::Vertices(TopoDS::Edge(one.rail.edge.Oriented(TopAbs_FORWARD)), ends[0], ends[1]);
    TopExp::Vertices(TopoDS::Edge(two.rail.edge.Oriented(TopAbs_FORWARD)), endsB[0], endsB[1]);
    const TopoDS_Vertex& V = ends[ea];
    if (!V.IsSame(endsB[eb]))
      fail(op,
           no_solution,
           "the loop where the blend folds crosses no joint of its edges",
           {one.rail.edge});
    const double ja    = one.rail.u(BRep_Tool::Parameter(V, one.rail.edge)),
                 jb    = two.rail.u(BRep_Tool::Parameter(V, two.rail.edge));
    const gp_Pnt P     = one.sheet.contact(other)->Value(one.sheet.rail(ja));
    const double apart = P.Distance(two.sheet.contact(otherB)->Value(two.sheet.rail(jb)));
    if (apart > one.rail.resolved(limit::together))
      fail(op,
           unsupported_geometry,
           "the blends of tangent edges do not meet in one section there, which the blend does not "
           "do yet",
           {one.rail.edge});
    TopoDS_Vertex J;
    builder.MakeVertex(J, P, std::max(tolerance, apart));
    runOn(one, other, corner[a][other][ea], ua, J, ja);
    runOn(two, otherB, J, jb, corner[b][otherB][eb], ub);
    for (const auto& shape : edgesOfVertex.FindFromKey(V))
    {
      if (railEdge(shape))
        continue;
      if (!onFace(shape, one.rail.s[other].face) && !onFace(shape, two.rail.s[otherB].face))
        fail(op,
             unsupported_geometry,
             "where the loop of a fold crosses a joint another edge meets the cap there, which the "
             "blend does not do yet",
             {one.rail.edge});
      double                      lo, hi;
      const auto                  curve = BRep_Tool::Curve(TopoDS::Edge(shape), lo, hi);
      GeomAPI_ProjectPointOnCurve on(P, curve, lo, hi);
      if (on.NbPoints() == 0 || on.LowerDistance() > one.rail.resolved(limit::together))
        fail(
          op,
          unsupported_geometry,
          "the contact of a chain's blend leaves the edge between its faces there, which the blend "
          "does not do yet",
          {one.rail.edge});
      cut(shape, V, J, on.LowerDistanceParameter());
    }
  }
  folded.push_back(f);
}

//! Where a chamfer's contact on the cap crosses itself at X, its two sheets (one surface, the
//! sections before the loop and after it) cross each other along a crease from X, and every
//! section of the loop between lies inside the material the others cut away, above that crease:
//! the chamfer is the outer envelope of its sections, which on the cap is the contact trimmed at
//! X (the points the setback away from the edge and no nearer, the distance set). The crease is
//! the curve S(ra, va) = S(rb, vb) of the sheet S(rail, across), traced from X (ra and rb the
//! sections of piece a's end and piece b's start, va = vb on the cap) as the gap rb - ra between
//! the two sections closes, to where they are one section (the swallowtail point V of the
//! sheet). Beyond V the sheet folds no more: the two faces meet along that section, from V to the
//! other contact, and the loop's sections are taken by both faces up to it. Each face is bounded
//! there by the crease and that section, an end of two edges (chain).
void Attach::creaseFold(int a, int ea, int b, int eb)
{
  const Piece& one = pieces[a];
  const Piece& two = pieces[b];
  const int    cap = one.joint[ea].side, other = 1 - cap;
  if (two.joint[eb].side != cap)
    fail(op, no_solution, "the two sides of a chamfer's fold are not on one cap", {one.rail.edge});
  const Sheet& sa = one.sheet;
  const Sheet& sb = two.sheet;
  const double c = sa.across[cap], o = sa.across[other];
  if (std::abs(sb.across[cap] - c) > tol::parameter
      || std::abs(sb.across[other] - o) > tol::parameter)
    fail(op,
         no_solution,
         "the two sides of a chamfer's fold are not on one sheet",
         {one.rail.edge});
  const auto at = [](const Sheet& sheet, double r, double v, gp_Pnt& p, gp_Vec& pr, gp_Vec& pv) {
    const gp_Pnt2d q = sheet.uv(r, v);
    gp_Vec         du, dv;
    sheet.surface->D1(q.X(), q.Y(), p, du, dv);
    pr = sheet.swapped ? dv : du;
    pv = sheet.swapped ? du : dv;
  };
  // Both branches are traced on piece a's surface: the pieces' sheets are one surface, each
  // with its own parameters for the fit's u (Sheet::rail), which carry b's back to its own.
  double       ra = sa.rail(one.uAt(ea)), rb = sa.rail(two.uAt(eb)), va = c, vb = c;
  const gp_Pnt X   = sa.contact(cap)->Value(ra);
  const double gap = X.Distance(sb.contact(cap)->Value(sb.rail(two.uAt(eb))));
  const auto   onB = [&](double r) { return sb.rail((r - sa.alpha) / sa.beta); };
  if (gap > one.rail.resolved(limit::together))
    fail(op,
         no_solution,
         "the contact on the cap does not cross itself where the blend folds",
         {one.rail.edge});
  // The gap closes from its value at X to nothing: each step a share of what is left, more finely
  // as the sections close (the crease turns into the sheet's cuspidal edge there), until the two
  // sections lie within the length tolerance of each other.
  const double first = rb - ra;

  struct Place
  {
    double ra, va, rb, vb;
    gp_Pnt p;
  };

  std::vector<Place>    path = {{ra, va, rb, vb, X}};
  std::array<double, 4> seed{}; // the first step: ra, va, vb and the gap
  const double          length = tolerance;
  double                delta = first, share = 1.0 / 16;
  // Two branches of S(ra, va) = S(rb, vb) leave X: the crease into the sheet, along the line the
  // two sheets' tangent planes there share, towards the other contact, and the loop's tail beyond
  // the cap. The first step goes along that line into the sheet, a sixteenth of the section's
  // length, each foot moved in its own sheet's tangent plane (least squares on its two tangents).
  {
    gp_Pnt pa, pb;
    gp_Vec ar, av, br, bv;
    at(sa, ra, c, pa, ar, av);
    at(sa, rb, c, pb, br, bv);
    gp_Vec d = ar.Crossed(av).Crossed(br.Crossed(bv));
    if (!(d.Magnitude() > tol::angle * ar.Crossed(av).Magnitude() * br.Crossed(bv).Magnitude()))
      fail(op,
           no_solution,
           "the chamfer's sheets touch where its contact crosses itself",
           {one.rail.edge});
    d.Normalize();
    if (d.Dot(av) * (o - c) < 0)
      d.Reverse();
    gp_Pnt far;
    gp_Vec fr, fv;
    at(sa, ra, o, far, fr, fv);
    const gp_Vec h  = d * (X.Distance(far) / 16);
    const auto   in = [&](const gp_Vec& tr, const gp_Vec& tv, double& dr, double& dv) {
      const double a11 = tr.Dot(tr), a12 = tr.Dot(tv), a22 = tv.Dot(tv),
                   det = a11 * a22 - a12 * a12;
      if (!(std::abs(det) > 0))
        return false;
      dr = (h.Dot(tr) * a22 - h.Dot(tv) * a12) / det;
      dv = (h.Dot(tv) * a11 - h.Dot(tr) * a12) / det;
      return true;
    };
    double dra, dva, drb, dvb;
    if (!in(ar, av, dra, dva) || !in(br, bv, drb, dvb))
      fail(op,
           no_solution,
           "the chamfer's sheets are singular where its contact crosses itself",
           {one.rail.edge});
    seed = {ra + dra, va + dva, vb + dvb, (rb + drb) - (ra + dra)};
  }
  for (int guard = 0;; ++guard)
  {
    checkCancelled(op);
    gp_Pnt pa, pb;
    gp_Vec ar, av, br, bv;
    at(sa, path.back().ra, path.back().va, pa, ar, av);
    if (std::abs(delta) * std::max(ar.Magnitude(), length) <= length || guard > 4096)
      break;
    const double target = path.size() == 1 ? seed[3] : delta * (1 - share);
    // From the place before, along the step the last two places took (the first from X along
    // the crease's tangent).
    double x[3] = {path.back().ra, path.back().va, path.back().vb};
    if (path.size() == 1)
    {
      x[0] = seed[0];
      x[1] = seed[1];
      x[2] = seed[2];
    }
    else
    {
      const Place& p0 = path[path.size() - 2];
      const Place& p1 = path.back();
      const double f  = (target - delta) / ((p1.rb - p1.ra) - (p0.rb - p0.ra));
      x[0] += (p1.ra - p0.ra) * f;
      x[1] += (p1.va - p0.va) * f;
      x[2] += (p1.vb - p0.vb) * f;
    }
    bool solved = false;
    for (int it = 0; it < 40; ++it)
    {
      at(sa, x[0], x[1], pa, ar, av);
      at(sa, x[0] + target, x[2], pb, br, bv);
      const gp_Vec F(pa, pb);
      if (F.Magnitude() <= 1e-3 * length)
      {
        solved = true;
        break;
      }
      // d(pb - pa) / d(ra, va, vb), the gap held.
      const gp_Vec cols[3] = {br - ar, -av, bv};
      math_Matrix  J(1, 3, 1, 3);
      math_Vector  rhs(1, 3);
      for (int r = 0; r < 3; ++r)
      {
        for (int k = 0; k < 3; ++k)
          J(r + 1, k + 1) = cols[k].Coord(r + 1);
        rhs(r + 1) = -F.Coord(r + 1);
      }
      math_SVD svd(J);
      if (!svd.IsDone())
        break;
      math_Vector dx(1, 3);
      svd.Solve(rhs, dx);
      for (int k = 0; k < 3; ++k)
        x[k] += dx(k + 1);
      if (!std::isfinite(x[0] + x[1] + x[2]))
        break;
    }
    // On the sheet between the cap and the other contact, and both sections between X's.
    const bool inside = solved && (x[1] - c) * (o - c) >= -tol::parameter
                        && (x[2] - c) * (o - c) >= -tol::parameter && (x[1] - o) * (o - c) < 0
                        && (x[2] - o) * (o - c) < 0;
    if (!inside)
    {
      if (path.size() == 1)
      {
        seed[0] = ra + (seed[0] - ra) / 2;
        seed[1] = va + (seed[1] - va) / 2;
        seed[2] = vb + (seed[2] - vb) / 2;
        seed[3] = first + (seed[3] - first) / 2;
      }
      share /= 2;
      if (share < 1e-6)
        fail(op,
             no_solution,
             "the crease where the chamfer's sheets cross cannot be traced to where they part",
             {one.rail.edge});
      continue;
    }
    path.push_back({x[0], x[1], x[0] + target, x[2], pa.Translated(gp_Vec(pa, pb) / 2)});
    delta = target;
    share = std::min(0.25, share * 1.5);
  }
  if (path.size() < 4)
    fail(op,
         no_solution,
         "the crease where the chamfer's sheets cross is too short to make",
         {one.rail.edge});
  // The crease's end: where the sections are one (within the tolerance), its place on both.
  const Place&                 end = path.back();
  const double                 rm = (end.ra + end.rb) / 2, vm = (end.va + end.vb) / 2;
  const int                    count = int(path.size());
  NCollection_Array1<gp_Pnt>   points(1, count);
  NCollection_Array1<gp_Pnt2d> onA(1, count), onBs(1, count);
  NCollection_Array1<double>   params(1, count);
  double                       along = 0;
  for (int i = 0; i < count; ++i)
  {
    if (i > 0)
      along += path[i].p.Distance(path[i - 1].p);
    params(i + 1) = along;
    points(i + 1) = path[i].p;
    onA(i + 1)    = sa.uv(path[i].ra, path[i].va);
    onBs(i + 1)   = sb.uv(onB(path[i].rb), path[i].vb);
  }
  // Repeated places at the end (steps under the resolution) are merged: interpolation needs
  // rising parameters.
  std::vector<int> keep;
  for (int i = 1; i <= count; ++i)
    if (keep.empty() || params(i) - params(keep.back()) > 1e-3 * length || i == count)
    {
      if (!keep.empty() && i == count && params(i) - params(keep.back()) <= 1e-3 * length)
        keep.back() = i;
      else
        keep.push_back(i);
    }
  const int m = int(keep.size());
  if (m < 4)
    fail(op,
         no_solution,
         "the crease where the chamfer's sheets cross is too short to make",
         {one.rail.edge});
  const Handle(NCollection_HArray1<gp_Pnt>)   p3  = new NCollection_HArray1<gp_Pnt>(1, m);
  const Handle(NCollection_HArray1<gp_Pnt2d>) pa2 = new NCollection_HArray1<gp_Pnt2d>(1, m),
                                              pb2 = new NCollection_HArray1<gp_Pnt2d>(1, m);
  const Handle(NCollection_HArray1<double>) t     = new NCollection_HArray1<double>(1, m);
  for (int i = 1; i <= m; ++i)
  {
    p3->SetValue(i, points(keep[i - 1]));
    pa2->SetValue(i, onA(keep[i - 1]));
    pb2->SetValue(i, onBs(keep[i - 1]));
    t->SetValue(i, params(keep[i - 1]));
  }
  Handle(Geom_Curve)   curve;
  Handle(Geom2d_Curve) curveA, curveB;
  try
  {
    GeomAPI_Interpolate   i3(p3, t, false, tol::parameter);
    Geom2dAPI_Interpolate ia(pa2, t, false, tol::parameter), ib(pb2, t, false, tol::parameter);
    i3.Perform();
    ia.Perform();
    ib.Perform();
    if (!i3.IsDone() || !ia.IsDone() || !ib.IsDone())
      fail(op,
           no_solution,
           "the crease where the chamfer's sheets cross cannot be interpolated",
           {one.rail.edge});
    curve  = i3.Curve();
    curveA = ia.Curve();
    curveB = ib.Curve();
  }
  catch (const Standard_Failure&)
  {
    fail(op,
         no_solution,
         "the crease where the chamfer's sheets cross cannot be interpolated",
         {one.rail.edge});
  }
  const double tEnd = t->Value(m);
  // X, V and W (where the section through V meets the other contact).
  builder.MakeVertex(corner[a][cap][ea], X, std::max(tolerance, gap));
  corner[b][cap][eb] = corner[a][cap][ea];
  gp_Pnt V, W;
  gp_Vec d1, d2;
  at(sa, rm, vm, V, d1, d2);
  at(sa, rm, o, W, d1, d2);
  TopoDS_Vertex vertexV;
  builder.MakeVertex(vertexV, V, tolerance);
  builder.MakeVertex(corner[a][other][ea], W, tolerance);
  corner[b][other][eb] = corner[a][other][ea];
  u[a][cap][ea]        = one.uAt(ea);
  u[b][cap][eb]        = two.uAt(eb);
  u[a][other][ea]      = (rm - sa.alpha) / sa.beta;
  u[b][other][eb]      = (rm - sa.alpha) / sa.beta;
  const TopoDS_Edge creaseEdge =
    edgeBetween(curve,
                corner[a][cap][ea],
                0,
                vertexV,
                tEnd,
                {{curveA, sa.surface, TopLoc_Location()}, {curveB, sb.surface, TopLoc_Location()}},
                tolerance);
  const TopoDS_Edge beyond = edgeBetween(sa.section(rm),
                                         vertexV,
                                         vm,
                                         corner[a][other][ea],
                                         o,
                                         {{sa.sectionUV(rm), sa.surface, TopLoc_Location()},
                                          {sb.sectionUV(onB(rm)), sb.surface, TopLoc_Location()}},
                                         tolerance);
  // An end of two edges, from side 0's corner to side 1's.
  const std::vector<TopoDS_Edge> ends = cap == 0 ? std::vector<TopoDS_Edge>{creaseEdge, beyond}
                                                 : std::vector<TopoDS_Edge>{beyond, creaseEdge};
  chain[a][ea] = chain[b][eb] = ends;
  section[a][ea] = section[b][eb] = ends.front();
  Folded f{a, ea, b, eb, {}, {}, {creaseEdge, beyond}};
  folded.push_back(f);
}

//! Where a round crosses itself between its contacts (Fold::crease): piece `a`'s end `ea` and
//! piece `b`'s end `eb` are the sheet before the fold and after it, each on its own copy of the one
//! surface. The blend is the boundary of the balls' union, so the two meet along the crease where
//! the sheet crosses itself, between its swallowtail points P0 and P1, which lies inside the
//! sheet; from each of those to the contact beside it they meet along the section through it, where
//! the sheet is one. Each face's end is so
//! three edges from side 0's contact to side 1's: the section through P0, the crease, the section
//! through P1 (chain).
void Attach::crossingFold(int a, int ea, int b, int eb)
{
  const Piece&  one    = pieces[a];
  const Piece&  two    = pieces[b];
  const Crease& crease = *one.joint[ea].crease;
  const Sheet&  sa     = one.sheet;
  const Sheet&  sb     = two.sheet;
  if (sa.swapped || sb.swapped || sa.alpha != 0 || sa.beta != 1 || sb.alpha != 0 || sb.beta != 1
      || sa.across[0] != 0 || sa.across[1] != 1 || sb.across[0] != 0 || sb.across[1] != 1)
    fail(op,
         no_solution,
         "the two sides of a blend's crease are not on its fitted sheet",
         {one.rail.edge});
  // The sections through the crease's ends, and the contacts' corners there.
  TopoDS_Vertex ends[2];
  TopoDS_Edge   sections[2];
  for (int side = 0; side < 2; ++side)
  {
    const gp_Pnt2d& end = crease.ends[side];
    const double    r   = end.X();
    builder.MakeVertex(ends[side], crease.curve->Value(side == 0 ? crease.first : crease.last),
                       tolerance);
    builder.MakeVertex(corner[a][side][ea], sa.contact(side)->Value(r), tolerance);
    corner[b][side][eb] = corner[a][side][ea];
    u[a][side][ea]      = r;
    u[b][side][eb]      = r;
    const TopoDS_Vertex& from = side == 0 ? corner[a][0][ea] : ends[1];
    const TopoDS_Vertex& to   = side == 0 ? ends[0] : corner[a][1][ea];
    sections[side]            = edgeBetween(sa.section(r),
                                 from,
                                 side == 0 ? sa.across[0] : end.Y(),
                                 to,
                                 side == 0 ? end.Y() : sa.across[1],
                                 {{sa.sectionUV(r), sa.surface, TopLoc_Location()},
                                  {sb.sectionUV(r), sb.surface, TopLoc_Location()}},
                                 tolerance);
  }
  const TopoDS_Edge creaseEdge =
    edgeBetween(crease.curve,
                ends[0],
                crease.first,
                ends[1],
                crease.last,
                {{Handle(Geom2d_Curve)(crease.onA), sa.surface, TopLoc_Location()},
                 {Handle(Geom2d_Curve)(crease.onB), sb.surface, TopLoc_Location()}},
                tolerance);
  chain[a][ea] = chain[b][eb] = {sections[0], creaseEdge, sections[1]};
  section[a][ea] = section[b][eb] = sections[0];
  folded.push_back(Folded{a, ea, b, eb, {}, {}, {creaseEdge, sections[0], sections[1]}});
}

//! The bicubic through a grid of points at uniform parameters in [0, 1] both ways, its knots
//! at the inner parameters (de Boor's not-a-knot ends); empty where it is singular.
//! (GeomAPI_PointsToBSplineSurface's interpolation sets end tangents that ripple back from the
//! ends across the first spans by up to the spacing itself.)
static Handle(Geom_BSplineSurface) interpolated(const NCollection_Array2<gp_Pnt>& grid)
{
  constexpr int degree = 3;
  const auto    basis  = [](int                         count,
                        NCollection_Array1<double>& knots,
                        NCollection_Array1<int>&    mults,
                        NCollection_Array1<double>& flat,
                        NCollection_Array1<double>& sites) {
    for (int i = 1; i <= count; ++i)
      sites(i) = double(i - 1) / (count - 1);
    knots(1) = 0;
    mults(1) = degree + 1;
    for (int i = 2; i < count - 2; ++i)
    {
      knots(i) = sites(i + 1);
      mults(i) = 1;
    }
    knots(count - 2) = 1;
    mults(count - 2) = degree + 1;
    BSplCLib::KnotSequence(knots, mults, flat);
  };
  const int nu = grid.ColLength(), nv = grid.RowLength();
  if (nu < degree + 1 || nv < degree + 1)
    return {};
  NCollection_Array1<double> uKnots(1, nu - 2), vKnots(1, nv - 2), uFlat(1, nu + degree + 1),
    vFlat(1, nv + degree + 1), uSites(1, nu), vSites(1, nv);
  NCollection_Array1<int> uMults(1, nu - 2), vMults(1, nv - 2), uOrders(1, nu), vOrders(1, nv);
  basis(nu, uKnots, uMults, uFlat, uSites);
  basis(nv, vKnots, vMults, vFlat, vSites);
  uOrders.Init(0);
  vOrders.Init(0);
  NCollection_Array2<gp_Pnt> poles(1, nu, 1, nv);
  // Along v for each row, then along u for each column of those poles.
  for (int i = 1; i <= nu; ++i)
  {
    NCollection_Array1<gp_Pnt> row(1, nv);
    for (int j = 1; j <= nv; ++j)
      row(j) = grid(grid.LowerRow() + i - 1, grid.LowerCol() + j - 1);
    int singular = 0;
    BSplCLib::Interpolate(degree, vFlat, vSites, vOrders, row, singular);
    if (singular != 0)
      return {};
    for (int j = 1; j <= nv; ++j)
      poles(i, j) = row(j);
  }
  for (int j = 1; j <= nv; ++j)
  {
    NCollection_Array1<gp_Pnt> column(1, nu);
    for (int i = 1; i <= nu; ++i)
      column(i) = poles(i, j);
    int singular = 0;
    BSplCLib::Interpolate(degree, uFlat, uSites, uOrders, column, singular);
    if (singular != 0)
      return {};
    for (int i = 1; i <= nu; ++i)
      poles(i, j) = column(i);
  }
  return new Geom_BSplineSurface(poles, uKnots, vKnots, uMults, vMults, degree, degree);
}

//! A round's patch over the loop where its contact on the cap folded: the ball where the contact
//! crosses itself, whose arcs both sections are and which touches both blends along them (so it is
//! G1 to them as it is), carried out radially to meet the other contact, which runs outside it,
//! with that face's tangent plane. Over u from the place on the cap (0, where the surface closes
//! to a point) to the other contact (1), and v from the first section (0) to the second (1): with
//! d(u, v) the direction from the centre turning from the cap's place to the contact at v, and
//! g(v) how far the contact lies beyond the ball, the point is C + (r + h1(u) g(v) + h2(u) k(v)) d,
//! h1 = 3u^2 - 2u^3 and h2 = u^2 (u - 1), k chosen so the radial slope at u = 1 lies in that face.
//! The contact is where the ball touches that face at both sections, so g and its slope vanish
//! there and the sections stay the ball's. Interpolated at a grid; empty where it does not hold.
Handle(Geom_BSplineSurface) Attach::ballPatch(const Folded& f, int nv)
{
  const Piece&   one  = pieces[f.a];
  const Rail&    rail = one.rail;
  const int      cap  = one.joint[f.ea].side;
  const double   r    = rail.s[cap].distance;
  Jet            foot;
  const gp_Pnt2d onCap = one.fitted.contact[cap]->Value(one.uAt(f.ea));
  if (!jetAt(rail.s[cap], onCap.X(), onCap.Y(), foot))
    return {};
  const gp_Pnt  C  = foot.p.Translated(foot.n * (rail.sigma * r));
  const gp_Vec  up = gp_Vec(C, BRep_Tool::Pnt(corner[f.a][cap][f.ea])) / r;
  constexpr int nu = 33;
  std::vector<std::vector<std::pair<double, double>>> onEdge;
  const auto                                          places = runPlaces(f, nv, onEdge);
  NCollection_Array2<gp_Pnt>                          grid(1, nu, 1, nv);
  for (int j = 1; j <= nv; ++j)
  {
    const Run& run = f.run[places[j - 1].edge];
    double     t0, t1;
    const auto curve     = BRep_Tool::Curve(run.edge, t0, t1);
    const auto onSupport = BRep_Tool::CurveOnSurface(run.edge, run.side->face, t0, t1);
    if (curve.IsNull() || onSupport.IsNull())
      return {};
    const double   t  = places[j - 1].t;
    const gp_Pnt   Q  = curve->Value(t);
    const gp_Pnt2d uv = onSupport->Value(t);
    Jet            at;
    if (!jetAt(*run.side, uv.X(), uv.Y(), at))
      return {};
    const double far   = C.Distance(Q);
    const gp_Vec to    = gp_Vec(C, Q) / far;
    const double angle = up.Angle(to);
    if (!(angle > tol::angle) || !(angle < M_PI - tol::angle))
      return {};
    // The direction's turn at u = 1, and the radial slope there that lies in the face.
    const gp_Vec turn = (up * -angle + to * (angle * std::cos(angle))) / std::sin(angle);
    if (!(std::abs(at.n.Dot(to)) > tol::angle))
      return {};
    const double k = -far * at.n.Dot(turn) / at.n.Dot(to);
    for (int i = 1; i <= nu; ++i)
    {
      const double u = double(i - 1) / (nu - 1);
      const gp_Vec d =
        (up * std::sin((1 - u) * angle) + to * std::sin(u * angle)) / std::sin(angle);
      const double h1 = u * u * (3 - 2 * u), h2 = u * u * (u - 1);
      grid(i, j) = C.Translated(d * (r + h1 * (far - r) + h2 * k));
    }
  }
  return interpolated(grid);
}

//! The face laid over a fold cut out of the rolling-ball blend:
//! for a round the ball there (ballPatch), else Patch's face, G2 (for a G2 blend, where Patch can
//! make it) or G1 to both blends along their sections and to the face the other contact runs on;
//! laid on the edges already made so that it shares them.
void Attach::foldPatches()
{
  for (auto& f : folded)
  {
    const TopoDS_Edge&       sa   = section[f.a][f.ea];
    const TopoDS_Edge&       sb   = section[f.b][f.eb];
    const Rail&              rail = pieces[f.a].rail;
    std::vector<TopoDS_Edge> loop = {sa};
    for (const auto& run : f.run)
      loop.push_back(run.edge);
    loop.push_back(sb);
    // Sampled more densely until the loop's edges lie on it within the patch tolerance.
    if (!f.crease.empty())
      continue; // a chamfer's sheets meet in their crease: no patch
    if (round(rail))
      for (int nv = 65; nv <= 257 && f.patch.IsNull(); nv = 2 * nv - 1)
        if (const auto ball = ballPatch(f, nv))
          if (auto face = ballFace(f, ball, nv))
            f.patch = *face;
    if (!f.patch.IsNull())
      continue;
    TopoDS_Face made;
    // As smooth as the blends where Patch can make it: a G2 blend's patch G2 to them and to
    // the face, else G1 (the joints record G1 either way).
    for (int order = rail.law == Law::G2 ? 2 : 1; order >= 1 && made.IsNull(); --order)
    {
      try
      {
        std::vector<PatchSide> boundary = {{sa, blends[f.a], order}};
        for (const auto& run : f.run)
          boundary.push_back({run.edge, run.side->face, order});
        boundary.push_back({sb, blends[f.b], order});
        made = patchLoop(boundary, rail.resolved(limit::patch), limit::patchDegrees);
      }
      catch (const Failure& failure)
      {
        if (failure.reason == std::string("cancelled"))
          throw;
        if (order == 1)
          fail(op,
               no_solution,
               std::string("the loop where the blend folds cannot be patched: ") + failure.what(),
               {pieces[f.a].rail.edge});
      }
    }
    const Handle(Geom_Surface) surface = BRep_Tool::Surface(made);
    for (const TopoDS_Edge& edge : loop)
    {
      double     lo, hi;
      const auto curve   = BRep_Tool::Curve(edge, lo, hi);
      const auto onPatch = GeomProjLib::Curve2d(curve, lo, hi, surface);
      if (onPatch.IsNull())
        fail(op,
             no_solution,
             "the edges of the loop where the blend folds cannot be laid on its patch",
             {pieces[f.a].rail.edge});
      builder.UpdateEdge(edge, onPatch, surface, TopLoc_Location(), tolerance);
      builder.Range(edge, lo, hi);
    }
    BRepBuilderAPI_MakeWire wire;
    for (const TopoDS_Edge& edge : loop)
      wire.Add(edge);
    if (!wire.IsDone())
      fail(op,
           no_solution,
           "the edges of the loop where the blend folds do not close",
           {pieces[f.a].rail.edge});
    BRepBuilderAPI_MakeFace face(surface, wire.Wire(), true);
    if (!face.IsDone())
      fail(op,
           no_solution,
           "the patch where the blend folds cannot be bounded by its edges",
           {pieces[f.a].rail.edge});
    f.patch = face.Face();
    // Facing out as the blends do: across a shared edge the two faces run it opposite ways.
    const auto use = [](const TopoDS_Face& in, const TopoDS_Edge& edge) {
      for (TopExp_Explorer e(in, TopAbs_EDGE); e.More(); e.Next())
        if (e.Current().IsSame(edge))
          return e.Current().Orientation();
      return TopAbs_EXTERNAL;
    };
    if (use(f.patch, sa) == use(blends[f.a], sa))
      f.patch.Reverse();
  }
}

//! The ball patch's face on the loop's own edges: the sections along v = 0 and v = 1 (their curves
//! on it from the angle at the centre, as their points lie on the ball), the other contact along
//! u = 1 at its own parameter, and the place on the cap a degenerate edge along u = 0. Empty
//! where the edges do not lie on it (the tolerances stage records how near they do).
std::optional<TopoDS_Face> Attach::ballFace(Folded&                            f,
                                            const Handle(Geom_BSplineSurface)& surface,
                                            int                                nv)
{
  const Piece& one = pieces[f.a];
  const int    cap = one.joint[f.ea].side;
  double       u0, u1, v0, v1;
  surface->Bounds(u0, u1, v0, v1);
  const gp_Pnt         X         = BRep_Tool::Pnt(corner[f.a][cap][f.ea]);
  const TopoDS_Vertex& vx        = corner[f.a][cap][f.ea];
  const auto           onSection = [&](const TopoDS_Edge& edge, double v) -> Handle(Geom2d_Curve) {
    double       lo, hi;
    const auto   curve = BRep_Tool::Curve(edge, lo, hi);
    const gp_Pnt far   = BRep_Tool::Pnt(
      TopExp::FirstVertex(edge).IsSame(vx) ? TopExp::LastVertex(edge) : TopExp::FirstVertex(edge));
    Jet            foot;
    const gp_Pnt2d onCap = one.fitted.contact[cap]->Value(one.uAt(f.ea));
    if (!jetAt(one.rail.s[cap], onCap.X(), onCap.Y(), foot))
      return {};
    const gp_Pnt  C = foot.p.Translated(foot.n * (one.rail.sigma * one.rail.s[cap].distance));
    const gp_Vec  up(C, X), down(C, far);
    const double  whole = up.Angle(down);
    constexpr int count = 33;
    NCollection_Array1<gp_Pnt2d> values(1, count);
    NCollection_Array1<double>   at(1, count);
    for (int i = 1; i <= count; ++i)
    {
      const double t = lo + (hi - lo) * (i - 1) / (count - 1);
      at(i)          = t;
      values(i) = gp_Pnt2d(u0 + (u1 - u0) * up.Angle(gp_Vec(C, curve->Value(t))) / whole, v);
    }
    const Handle(NCollection_HArray1<gp_Pnt2d>) points = new NCollection_HArray1<gp_Pnt2d>(values);
    const Handle(NCollection_HArray1<double>)   parameters = new NCollection_HArray1<double>(at);
    Geom2dAPI_Interpolate interpolate(points, parameters, false, tol::parameter);
    interpolate.Perform();
    return interpolate.IsDone() ? Handle(Geom2d_Curve)(interpolate.Curve())
                                          : Handle(Geom2d_Curve)();
  };
  const TopoDS_Edge& sa = section[f.a][f.ea];
  const TopoDS_Edge& sb = section[f.b][f.eb];
  const auto         pa = onSection(sa, v0), pb = onSection(sb, v1);
  if (pa.IsNull() || pb.IsNull())
    return {};
  // The run along u = 1, each edge over its share of v by length (as ballPatch sampled it).
  std::vector<std::vector<std::pair<double, double>>> onEdge;
  runPlaces(f, nv, onEdge);
  std::vector<std::pair<TopoDS_Edge, Handle(Geom2d_Curve)>> laid = {{sa, pa}};
  for (size_t i = 0; i < f.run.size(); ++i)
  {
    const auto& at = onEdge[i];
    const auto  v  = [&](size_t k) { return v0 + (v1 - v0) * at[k].second; };
    if (at.size() == 2)
    {
      laid.push_back({f.run[i].edge,
                      segment2d(gp_Pnt2d(u1, v(0)), gp_Pnt2d(u1, v(1)), at[0].first, at[1].first)});
      continue;
    }
    NCollection_Array1<gp_Pnt2d> values(1, int(at.size()));
    NCollection_Array1<double>   parameters(1, int(at.size()));
    for (size_t k = 0; k < at.size(); ++k)
    {
      values(int(k) + 1)     = gp_Pnt2d(u1, v(k));
      parameters(int(k) + 1) = at[k].first;
    }
    try
    {
      Geom2dAPI_Interpolate interpolate(new NCollection_HArray1<gp_Pnt2d>(values),
                                        new NCollection_HArray1<double>(parameters),
                                        false,
                                        tol::parameter);
      interpolate.Perform();
      if (!interpolate.IsDone())
        return {};
      laid.emplace_back(f.run[i].edge, interpolate.Curve());
    }
    catch (const Standard_Failure&)
    {
      return {};
    }
  }
  laid.push_back({sb, pb});
  for (const auto& [edge, pc] : laid)
  {
    double a, b;
    BRep_Tool::Range(edge, a, b);
    builder.UpdateEdge(edge, pc, surface, TopLoc_Location(), tolerance);
    builder.Range(edge, a, b);
  }
  // The place on the cap, where the surface closes to a point.
  const TopoDS_Edge degenerate = degenerateEdge(segment2d(gp_Pnt2d(u0, v1), gp_Pnt2d(u0, v0), 0, 1),
                                                surface,
                                                TopLoc_Location(),
                                                vx,
                                                0,
                                                1,
                                                tolerance);
  // Round the loop: X along the first section to A, the contact to B, back up the second, and
  // the point.
  TopoDS_Wire wire;
  builder.MakeWire(wire);
  builder.Add(wire, TopExp::FirstVertex(sa).IsSame(vx) ? sa : TopoDS::Edge(sa.Reversed()));
  for (const auto& run : f.run)
    builder.Add(wire, run.forward ? run.edge : TopoDS::Edge(run.edge.Reversed()));
  builder.Add(wire, TopExp::LastVertex(sb).IsSame(vx) ? sb : TopoDS::Edge(sb.Reversed()));
  builder.Add(wire, degenerate);
  wire.Closed(true);
  TopoDS_Face face;
  builder.MakeFace(face, surface, tolerance);
  builder.Add(face, wire);
  // Each edge's real gap to the surface is its tolerance (Jackson), as the fit leaves it.
  ShapeAnalysis_Edge analysis;
  for (const auto& [edge, pc] : laid)
  {
    double deviation = 0;
    analysis.CheckSameParameter(edge, face, deviation);
    if (deviation > one.rail.resolved(limit::patch))
      return {};
    if (deviation > BRep_Tool::Tolerance(edge))
      builder.UpdateEdge(edge, deviation);
    builder.SameParameter(edge, true);
  }
  if (!BRepCheck_Analyzer(face).IsValid())
    return {};
  // Facing out as the blends do: across a shared edge the two faces run it opposite ways.
  const auto use = [](const TopoDS_Face& in, const TopoDS_Edge& edge) {
    for (TopExp_Explorer e(in, TopAbs_EDGE); e.More(); e.Next())
      if (e.Current().IsSame(edge))
        return e.Current().Orientation();
    return TopAbs_EXTERNAL;
  };
  if (use(face, sa) == use(blends[f.a], sa))
    face.Reverse();
  return face;
}

} // namespace BRepExt::blending::detail
