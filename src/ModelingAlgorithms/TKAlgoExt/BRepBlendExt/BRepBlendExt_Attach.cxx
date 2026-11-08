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

// Blend: The blends attached to the body, stage by stage.

#include <BRepBlendExt_Internal.hxx>
#include <BRepAlgoExt_PatchLoop.hxx>
#include <GC_MakeSegment.hxx>
#include <math_SVD.hxx>
#include <math_ComputeGaussPointsAndWeights.hxx>
#include <math_Vector.hxx>
#include <math_Matrix.hxx>
#include <Geom2dAPI_Interpolate.hxx>
#include <GeomAPI_Interpolate.hxx>

#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepLib.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BndLib_AddSurface.hxx>
#include <GeomLib_IsPlanarSurface.hxx>
#include <ShapeAnalysis_Edge.hxx>
#include <numeric>

namespace BRepExt::blending::detail
{

//! The blends of chains of edges attached to the body: each support cut back to its
//! contacts, each stopping face to the blend's end section, the other edges at a chain's inner
//! vertices to its corners, and the blend faces put in between, all on shared edges, so that
//! nothing but the faces beside the edges is touched. The faces it changed are listed with what
//! each became.
struct Attach
{
  Attach(TopoDS_Shape                     body,
         const std::vector<Piece>&        pieces,
         const std::vector<Corner>&       cornersAt,
         const Handle(BRepTools_History)& history,
         bool refineContacts);

  const bool refineContacts;

  // The stages, in order.
  //! Each end of each piece: a stop face cut back to the blend's end section, a joint or a
  //! split shared with the next piece.
  void         chainEnds();
  void         miters();
  void         cornerPatches();
  void         ridges();
  void         contacts();
  void         blendFaces();
  void         continuity();
  std::vector<Join> joins;
  void         tolerances();
  void         collisions();
  TopoDS_Shape assemble();

  void fold(int a, int ea, int b, int eb);
  void foldPatches();
  struct Folded;
  Handle(Geom_BSplineSurface) ballPatch(const Folded& f, int nv);
  void                             creaseFold(int a, int ea, int b, int eb);
  std::optional<TopoDS_Face>  ballFace(Folded&                            f,
                                       const Handle(Geom_BSplineSurface)& surface,
                                       int                                nv);
  TopoDS_Face                 lawCorner(const Corner& at, size_t c);
  std::pair<TopoDS_Vertex, double> contactsCross(const Corner& at, int i, int j);
  void                        splitCorner(const Corner& at, size_t c);
  TopoDS_Face                 tracedCorner(const Corner& at, size_t c,
                                         const TopoDS_Vertex (&Q)[3][3], double crossingGap);
  bool                        notch(int a, int ea, int b, int eb);
  void                        rollRound(int                  k,
                                        int                  ea,
                                        int                  j,
                                        int                  eb,
                                        int                  sa,
                                        int                  sb,
                                        const TopoDS_Vertex& vertex,
                                        const TopoDS_Edge&   third,
                                        const Meeting&       pa,
                                        const Meeting&       pb);
  void                        rollPatches();

  //! Whether a rail's section is the ball's own circle: the round (a conic of rho 1/2 at equal
  //! distances).
  static bool round(const Rail& rail)
  {
    return rail.law == Law::Conic && std::abs(rail.rho - 0.5) < limit::equal
           && std::abs(rail.s[0].distance - rail.s[1].distance) < limit::equal * rail.scale;
  }

  void        stop(int k, int end, const TopoDS_Vertex& vertex);
  bool        consumed(int                  k,
                       int                  end,
                       const TopoDS_Vertex& vertex,
                       const TopoDS_Face&   stop,
                       const TopoDS_Edge (&beside)[2]);
  TopoDS_Edge rolledEdge(int k, int side, const TopoDS_Face& stop);
  std::map<std::pair<int, int>, TopoDS_Face>
       rolled; // [piece][end] -> the stop face a notch's blend uses up there (consumed)
  void stops(int k, int end, const TopoDS_Vertex& vertex);

  bool railEdge(const TopoDS_Shape& edge) const
  {
    return std::any_of(pieces.begin(), pieces.end(), [&](const Piece& p) {
      return ownEdge(p, edge);
    });
  }

  static bool ownEdge(const Piece& piece, const TopoDS_Shape& edge)
  {
    return piece.rail.edge.IsSame(edge)
           || std::any_of(piece.swallowed.begin(),
                          piece.swallowed.end(),
                          [&](const TopoDS_Edge& e) { return e.IsSame(edge); });
  }

  bool onFace(const TopoDS_Shape& edge, const TopoDS_Shape& face) const
  {
    if (!facesOfEdge.Contains(edge))
      return false; // a round between faces that meet only at a vertex (Piece::apex)
    for (const auto& f : facesOfEdge.FindFromKey(edge))
      if (f.IsSame(face))
        return true;
    return false;
  }

  //! A seam from the vertex to a pole in the strip the blend replaces (poleInStrip) goes, with the
  //! pole's degenerate edge: the contact closes by itself.
  bool sweepToPole(const TopoDS_Edge& seam, const TopoDS_Vertex& vertex)
  {
    const TopoDS_Face face = TopoDS::Face(facesOfEdge.FindFromKey(seam).First());
    TopoDS_Vertex     a, z;
    TopExp::Vertices(seam, a, z);
    if (!BRep_Tool::IsClosed(seam, face))
      return false;
    const TopoDS_Vertex& pole = a.IsSame(vertex) ? z : a;
    for (const auto& other : edgesOfVertex.FindFromKey(pole))
    {
      if (!BRep_Tool::Degenerated(TopoDS::Edge(other)))
        continue;
      auto& edit = editOf(edits, face);
      edit.replace.Bind(seam, TopoDS_Shape());
      edit.replace.Bind(other, TopoDS_Shape());
      swept.push_back(seam);
      swept.push_back(other);
      swept.push_back(pole);
      return true;
    }
    return false;
  }

  void cut(const TopoDS_Shape& edge, const TopoDS_Vertex& from, const TopoDS_Vertex& to, double at)
  {
    const TopoDS_Edge current = TopoDS::Edge(cuts.Contains(edge) ? cuts.FindFromKey(edge) : edge);
    const TopoDS_Edge out     = trimmed(current, from, to, at);
    // Cut back from both ends to each other or past, the edge is gone: the blends at its ends
    // meet on it (the full round, which the blend does not make) or overlap.
    double f, l;
    BRep_Tool::Range(out, f, l);
    if (!(l - f > tol::parameter * std::max(1.0, std::abs(f) + std::abs(l))))
      fail(op,
           unsupported_geometry,
           "the blends at the two ends of an edge between them meet on it or reach past each "
           "other, which the blend does not make",
           {edge});
    if (cuts.Contains(edge))
      cuts.ChangeFromKey(edge) = out;
    else
      cuts.Add(edge, out);
  }

  //! The half of a piece's fit at one end.
  static std::pair<double, double> halfOf(const Piece& p, int end)
  {
    const double middle = p.uAt(0) / 2 + p.uAt(1) / 2;
    return end == 0 ? std::make_pair(0.0, middle) : std::make_pair(middle, 1.0);
  }

  //! Where piece k's contact on its face `side`, in the half at `end`, meets piece `other`'s
  //! sheet: the fit's u, nearest the end.
  std::optional<double> meetsSheet(int k, int end, int side, int other) const
  {
    const Piece& mine   = pieces[k];
    const auto [u0, u1] = halfOf(mine, end);
    GeomAPI_IntCS meet(new Geom_TrimmedCurve(mine.sheet.contact(side),
                                             mine.sheet.rail(std::min(u0, u1)),
                                             mine.sheet.rail(std::max(u0, u1))),
                       pieces[other].sheet.surface);
    if (!meet.IsDone())
      return {};
    std::optional<double> best;
    for (int i = 1; i <= meet.NbPoints(); ++i)
    {
      double su, sv, w;
      meet.Parameters(i, su, sv, w);
      const double q = (w - mine.sheet.alpha) / mine.sheet.beta;
      if (!best || std::abs(q - mine.uAt(end)) < std::abs(*best - mine.uAt(end)))
        best = q;
    }
    return best;
  }

  //! Where piece k's contact on its face `side` crosses, in the half at `end`, the contact of
  //! the piece it is mitred against there on that face (the mitre's inner end), nearest the end;
  //! none where the end is no mitre on that face.
  std::optional<Meeting> commonCrossing(int k, int end, int side) const
  {
    const auto& joint = pieces[k].joint[end];
    if (joint.meets != Meets::Miter)
      return {};
    const Piece& a = pieces[k];
    const Piece& b = pieces[joint.piece];
    for (int t = 0; t < 2; ++t)
    {
      if (!a.rail.s[side].face.IsSame(b.rail.s[t].face))
        continue;
      const auto [a0, a1] = halfOf(a, end);
      const auto [b0, b1] = halfOf(b, joint.end);
      return nearest(meetings(between(a.fitted.contact[side], a0, a1),
                              between(b.fitted.contact[t], b0, b1),
                              a.rail.s[side].surface),
                     a.uAt(end));
    }
    return {};
  }

  //! Whether the mitre at piece k's end cannot be made: its contacts on the face the two blends
  //! share, carried on past the vertex to where they cross (a re-entrant corner), run over a
  //! contact as it is kept before they get there, the other's or their own. The ball rolls round
  //! the vertex instead (rollRound). A contact is kept
  //! between its ends, or out to where a mitre at an end crosses it, where that mitre is made
  //! (judged with the contacts kept to their ends: `alone`).
  bool runsAcross(int k, int ea, bool alone) const
  {
    const auto& joint = pieces[k].joint[ea];
    if (joint.meets != Meets::Miter)
      return false;
    const int    j = joint.piece, eb = joint.end;
    const Piece& a  = pieces[k];
    const Piece& b  = pieces[j];
    int          sa = -1, sb = -1;
    for (int s = 0; s < 2; ++s)
      for (int t = 0; t < 2; ++t)
        if (a.rail.s[s].face.IsSame(b.rail.s[t].face))
        {
          sa = s;
          sb = t;
        }
    if (sa < 0)
      return false;
    const auto q = commonCrossing(k, ea, sa);
    if (!q)
      return false;
    const auto past = [](double x, double end, double to, int at) {
      return (at == 1 ? to > end : to < end) && x > std::min(end, to) && x < std::max(end, to);
    };
    const auto along = [&](double x, int p, int side) {
      double ends[2];
      for (int end = 0; end < 2; ++end)
      {
        const auto there = alone ? std::nullopt : commonCrossing(p, end, side);
        ends[end]        = there && !runsAcross(p, end, true) ? there->first : pieces[p].uAt(end);
      }
      return x >= std::min(ends[0], ends[1]) && x <= std::max(ends[0], ends[1]);
    };
    for (const auto& m : meetings(a.fitted.contact[sa], b.fitted.contact[sb], a.rail.s[sa].surface))
      if (std::abs(m.first - q->first) > tol::parameterMerge
          && ((past(m.first, a.uAt(ea), q->first, ea) && along(m.second, j, sb))
              || (past(m.second, b.uAt(eb), q->second, eb) && along(m.first, k, sa))))
        return true;
    // Nor over its own contact: a run-on that comes back across the wall it leaves.
    for (const auto& [p, side, end, to] :
         {std::tuple{k, sa, ea, q->first}, std::tuple{j, sb, eb, q->second}})
    {
      const Piece&              piece = pieces[p];
      Geom2dAPI_InterCurveCurve self(piece.fitted.contact[side], tol::parameter);
      for (int i = 1; i <= self.NbPoints(); ++i)
      {
        const auto&  point = self.Intersector().Point(i);
        const double u0 = point.ParamOnFirst(), u1 = point.ParamOnSecond();
        if ((past(u0, piece.uAt(end), to, end) && along(u1, p, side))
            || (past(u1, piece.uAt(end), to, end) && along(u0, p, side)))
          return true;
      }
    }
    return false;
  }

  //! A point of piece k's sheet in its parameters, the period taken nearest its end `end`.
  gp_Pnt2d onSheet(int k, int end, const gp_Pnt& point) const
  {
    const Piece&               piece = pieces[k];
    GeomAPI_ProjectPointOnSurf on(point, piece.sheet.surface);
    double                     pu, pv;
    on.LowerDistanceParameters(pu, pv);
    // The projection is good to OpenCascade's confusion, which a small blend's share is not:
    // a point on the sheet is found to rounding by Gauss-Newton from there, while it gains.
    for (int step = 0; step < 4; ++step)
    {
      gp_Pnt at;
      gp_Vec du, dv;
      piece.sheet.surface->D1(pu, pv, at, du, dv);
      const gp_Vec off(at, point);
      const double a = du.Dot(du), b = du.Dot(dv), c = dv.Dot(dv), det = a * c - b * b;
      if (!(det > 0))
        break;
      const double su = (c * off.Dot(du) - b * off.Dot(dv)) / det,
                   sv = (a * off.Dot(dv) - b * off.Dot(du)) / det;
      if (!(piece.sheet.surface->Value(pu + su, pv + sv).Distance(point) < at.Distance(point)))
        break;
      pu += su;
      pv += sv;
    }
    return nearPeriod(gp_Pnt2d(pu, pv),
                      piece.sheet.uv(piece.sheet.rail(piece.uAt(end)),
                                     (piece.sheet.across[0] + piece.sheet.across[1]) / 2),
                      piece.sheet.surface);
  }

  TopoDS_Shape               body;
  const std::vector<Piece>&  pieces;
  const std::vector<Corner>& cornersAt;
  Handle(BRepTools_History)  history;
  const char*                op;
  const double               tolerance = tol::length;
  BRep_Builder               builder;
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
                                                           edgesOfVertex, facesOfEdge;
  const int                                                n;
  std::vector<std::array<std::array<TopoDS_Vertex, 2>, 2>> corner; // [piece][side][end]
  std::vector<std::array<std::array<double, 2>, 2>>        u;
  std::vector<std::array<TopoDS_Edge, 2>>                  section; // [piece][end]
  std::vector<std::array<TopoDS_Edge, 2>>                  contact; // [piece][side]
  std::vector<std::pair<TopoDS_Edge, std::vector<TopoDS_Face>>>
                                          continued; // edges carried past a stop's vertex
  std::vector<std::array<TopoDS_Face, 2>> stopAt;    // [piece][end]
  NCollection_IndexedDataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>
    cuts; // the body's edges cut back at one end or both, the same in every face; in the order they
          // were cut
  Edits edits;
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>
                           crossed; // edges contacts cross onto the next face
  std::vector<TopoDS_Face> patches; // [corner]
  //! A corner of several faces (splitCorner): its faces, each blend's end at it with the face it
  //! bounds there, and the edges between its faces.
  std::vector<std::vector<TopoDS_Face>>                          parts; // [corner]
  std::vector<std::tuple<TopoDS_Edge, int, TopoDS_Face>>         partEnds;  // edge, piece, part
  std::vector<std::tuple<TopoDS_Edge, TopoDS_Face, TopoDS_Face>> partInner; // edge, parts
  //! The corners' feet, one vertex for each place (within the length resolution): where the
  //! balls of two corners are one (a face used up between them, as a cube's round to half its
  //! side is a ball) they share them.
  std::vector<TopoDS_Vertex> feet;

  TopoDS_Vertex footAt(const gp_Pnt& p, double scale)
  {
    for (const auto& v : feet)
      if (BRep_Tool::Pnt(v).Distance(p) <= limit::one(scale, p.XYZ().Modulus()))
        return v;
    TopoDS_Vertex v;
    builder.MakeVertex(v, p, tolerance);
    return feet.emplace_back(v);
  }

  //! A blend between two corners whose balls are one: it has no length and no face, and the
  //! corners share its section.
  bool vanished(int k) const
  {
    return pieces[k].joint[0].meets == Meets::Corner && pieces[k].joint[1].meets == Meets::Corner
           && !section[k][0].IsNull() && section[k][0].IsSame(section[k][1]);
  }

  //! Faces used up between blends whose contacts on them are one, or that have no length.
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> usedUp;

  //! Another blend's contact on `face` from `from` to `to` through `middle`, else null.
  TopoDS_Edge sharedContact(const TopoDS_Face&   face,
                            const TopoDS_Vertex& from,
                            const TopoDS_Vertex& to,
                            const gp_Pnt&        middle,
                            double               scale) const
  {
    for (int m = 0; m < n; ++m)
      for (int t = 0; t < 2; ++t)
      {
        const TopoDS_Edge& e = contact[m][t];
        if (e.IsNull() || !pieces[m].rail.s[t].face.IsSame(face))
          continue;
        TopoDS_Vertex a, z;
        TopExp::Vertices(e, a, z);
        if (!((a.IsSame(from) && z.IsSame(to)) || (a.IsSame(to) && z.IsSame(from))))
          continue;
        double f, l;
        // A contact that has no 3D curve (a degenerate edge at a pole or a seam) cannot hold the
        // point.
        const Handle(Geom_Curve) curve = BRep_Tool::Curve(e, f, l);
        if (curve.IsNull())
          continue;
        GeomAPI_ProjectPointOnCurve on(middle, curve, f, l);
        if (on.NbPoints() > 0 && on.LowerDistance() <= limit::one(scale, middle.XYZ().Modulus()))
          return e;
      }
    return {};
  }

  //! Where a contact on a flat cap folded over itself: the piece ending at the loop on one side
  //! (a at its end ea) and the one on the other (b at eb), the other contact's run round the loop
  //! from a's section to b's (one edge, or one on each piece's face where the loop runs across a
  //! joint), and the patch over the loop.
  struct Run
  {
    TopoDS_Edge edge;
    const Side* side;    // the support it lies on
    bool        forward; // from a's section towards b's
  };

  struct Folded
  {
    int              a, ea, b, eb;
    std::vector<Run> run;
    TopoDS_Face      patch;
    //! A chamfer's: the crease where its two sheets cross, from X, and the section beyond it to
    //! the other contact, where the two pieces' faces meet on the one sheet (creaseFold).
    std::vector<TopoDS_Edge> crease;
  };

  std::vector<Folded> folded;

  //! A ball rolled round a re-entrant vertex between two blends that cannot be mitred there
  //! (rollRound): the pieces and their ends, their sides on the face they share, the vertex,
  //! the place X on the third edge where both other contacts cross it, the contact on the shared
  //! face between the blends' ends, and the surface the ball sweeps (a torus about the shared
  //! face's normal through X) over u0..u1 about its axis, closing to X at vX.
  struct Roll
  {
    int                          a, ea, b, eb, sa, sb;
    TopoDS_Vertex                vertex, X;
    TopoDS_Edge                  arc;
    Handle(Geom_ToroidalSurface) torus;
    double                       ua, ub, vX;
    TopoDS_Face                  patch;
  };

  std::vector<Roll> rolls;

  //! The roll a piece's end ends in, if any.
  const Roll* rollAt(int k, int end) const
  {
    for (const auto& r : rolls)
      if ((r.a == k && r.ea == end) || (r.b == k && r.eb == end))
        return &r;
    return nullptr;
  }

  //! The loop a piece's end meets.
  const Folded& foldAt(int k, int end) const
  {
    for (const auto& f : folded)
      if ((f.a == k && f.ea == end) || (f.b == k && f.eb == end))
        return f;
    fail(op, no_solution, "a blend's end at a fold has no loop", {pieces[k].rail.edge});
  }

  //! A loop's run at `count` equal steps of its length from a's section towards b's (the patch's
  //! v, from 0 to 1): each place's edge and parameter there, with each edge's own parameter and v
  //! at its ends and at the places on it, rising along the edge (its curve on the patch).
  struct RunPlace
  {
    size_t edge;
    double t;
  };

  std::vector<RunPlace> runPlaces(const Folded&                                        f,
                                  int                                                  count,
                                  std::vector<std::vector<std::pair<double, double>>>& onEdge) const
  {
    std::vector<double> length, start;
    double              total = 0;
    for (const auto& run : f.run)
    {
      start.push_back(total);
      length.push_back(GCPnts_AbscissaPoint::Length(BRepAdaptor_Curve(run.edge)));
      total += length.back();
    }
    onEdge.assign(f.run.size(), {});
    for (size_t i = 0; i < f.run.size(); ++i)
    {
      double t0, t1;
      BRep_Tool::Range(f.run[i].edge, t0, t1);
      const double v0 = start[i] / total, v1 = (start[i] + length[i]) / total;
      onEdge[i].push_back({t0, f.run[i].forward ? v0 : v1});
      onEdge[i].push_back({t1, f.run[i].forward ? v1 : v0});
    }
    std::vector<RunPlace> out;
    for (int j = 0; j < count; ++j)
    {
      const double at = total * j / (count - 1);
      size_t       i  = 0;
      while (i + 1 < f.run.size() && at > start[i] + length[i])
        ++i;
      double t0, t1;
      BRep_Tool::Range(f.run[i].edge, t0, t1);
      const double            along = std::clamp(at - start[i], 0.0, length[i]);
      const BRepAdaptor_Curve curve(f.run[i].edge);
      double                  t = f.run[i].forward ? t0 : t1;
      if (along > 0 && along < length[i])
      {
        GCPnts_AbscissaPoint place(curve, f.run[i].forward ? along : -along, t);
        if (place.IsDone())
          t = place.Parameter();
      }
      else if (along >= length[i])
        t = f.run[i].forward ? t1 : t0;
      out.push_back({i, t});
      // On the edge's curve on the patch, the places apart from its ends (a place merged
      // with an end is that end).
      const double apart = tol::parameterMerge * (t1 - t0);
      if (t > t0 + apart && t < t1 - apart)
        onEdge[i].push_back({t, double(j) / (count - 1)});
    }
    for (auto& places : onEdge)
      std::sort(places.begin(), places.end());
    return out;
  }

  //! Where the blends of two parts of an outline cross each other over a flat cap (a neck narrower
  //! than both: the cap's offset is islands): between each two places where their contacts on it
  //! cross, the two sheets meet in a crease, which takes both contacts' place there.
  struct Ridge
  {
    int           k, sk, m, sm; // the pieces and their sides on the cap
    double        uk[2], um[2];
    TopoDS_Vertex at[2];
    TopoDS_Edge   edge;
  };

  std::vector<Ridge> ridged;
  //! A side whose contact ridges cross: its contact's parts and the ridges, from its corner at the
  //! start to the one at the end.
  std::map<std::pair<int, int>, std::vector<TopoDS_Edge>> runs;

  //! The edges along a side, from its start's corner to its end's.
  std::vector<TopoDS_Edge> sideEdges(int k, int side) const
  {
    const auto it = runs.find({k, side});
    return it != runs.end() ? it->second : std::vector<TopoDS_Edge>{contact[k][side]};
  }

  std::vector<TopoDS_Face> blends; // [piece]
  //! An end made of several edges, from side 0's corner to side 1's (a notch's), else empty.
  std::vector<std::array<std::vector<TopoDS_Edge>, 2>> chain; // [piece][end]

  //! Each notch's edge where a blend's end runs on the other blend's face.
  struct Notch
  {
    int         blend;
    TopoDS_Edge edge;
    TopoDS_Face face;
  };

  std::vector<Notch>        notches;
  std::vector<TopoDS_Shape> swept; // edges and vertices of the body the blends swept over
  std::map<const TopoDS_TShape*, std::pair<TopoDS_Vertex, double>>
                           seamCrossed; // a seam's first crossing, until its second
  std::vector<TopoDS_Face> removed;     // faces a blend replaces whole (a cliff's narrow face)

  //! A notch's edge on the face that trims the blend: where each of the sheet's sections meets
  //! that face's surface, along the sheet's whole rail (its parameter that of the sheet's rail).
  struct TrimCurve
  {
    Handle(Geom_Curve)   curve;
    Handle(Geom2d_Curve) onSheet;
    std::vector<double>  across; // the section's own parameter there, at the samples
  };

  std::map<std::pair<int, int>, TrimCurve> trimCurves;
  //! Where a notch's narrow face goes on past an end: the far edge's vertex made where the stop
  //! face's surface crosses it [piece][side][end] (null where the strip ends there too).
  std::map<std::tuple<int, int, int>, TopoDS_Vertex> trimFar;
  const TrimCurve&                                   trimOf(int k, int side);

  //! Refused where the notch's edge does not run over the rail parameters s0 to s1: the face
  //! beyond does not reach the blend's sections there.
  void checkTrim(int k, int side, double s0, double s1)
  {
    const auto&  trim  = trimOf(k, side);
    const double slack = tol::parameter * (1 + std::abs(s0) + std::abs(s1));
    if (std::min(s0, s1) < trim.curve->FirstParameter() - slack
        || std::max(s0, s1) > trim.curve->LastParameter() + slack)
      fail(op,
           no_solution,
           "the face beyond a narrow face does not cross the blend's section",
           {pieces[k].rail.edge});
  }

  void capNotch(int                  k,
                int                  side,
                int                  end,
                const TopoDS_Vertex& vertex,
                const TopoDS_Face&   stop,
                gp_Pnt2d&            cornerOnStop);

  //! The section's own parameter where the notch's edge meets it at the sheet's rail parameter s.
  double trimAcross(int k, int side, double s)
  {
    const auto& trim = trimOf(k, side);
    checkTrim(k, side, s, s);
    const gp_Pnt2d uv = trim.onSheet->Value(s);
    return pieces[k].sheet.swapped ? uv.X() : uv.Y();
  }

  //! Whether piece k's side is notched from where a split found its contact running off the face
  //! over the far edge partway: the narrow face stays, the edge cut back to the crossing.
  bool splitNotch(int k, int side) const
  {
    const Side& s = pieces[k].rail.s[side];
    if (s.trim.IsNull())
      return false;
    return std::any_of(std::begin(pieces[k].joint),
                       std::end(pieces[k].joint),
                       [&](const Piece::Joint& j) {
                         return j.meets == Meets::Split && j.side == side
                                && j.crossed.IsSame(s.trimEdge);
                       });
  }

  //! A cliff's face beyond its edge: the one beside it that the blend does not replace.
  TopoDS_Face beyondCliff(const Side& side) const
  {
    for (const auto& face : facesOfEdge.FindFromKey(side.cliff))
      if (!face.IsSame(side.face))
        return TopoDS::Face(face);
    return TopoDS_Face();
  }
};

Attach::Attach(TopoDS_Shape                     body_,
               const std::vector<Piece>&        pieces_,
               const std::vector<Corner>&       cornersAt_,
               const Handle(BRepTools_History)& history_,
               bool refineContacts_)
    : refineContacts(refineContacts_),
      body(std::move(body_)),
      pieces(pieces_),
      cornersAt(cornersAt_),
      history(history_),
      op(pieces_.front().rail.op),
      n(int(pieces_.size())),
      corner(n),
      u(n),
      section(n),
      contact(n),
      stopAt(n),
      patches(cornersAt_.size()),
      parts(cornersAt_.size()),
      blends(n),
      chain(n)
{
  TopExp::MapShapesAndUniqueAncestors(body, TopAbs_VERTEX, TopAbs_EDGE, edgesOfVertex);
  TopExp::MapShapesAndUniqueAncestors(body, TopAbs_EDGE, TopAbs_FACE, facesOfEdge);
  for (const auto& piece : pieces)
    for (const auto& joint : piece.joint)
      if (joint.meets == Meets::Split)
        crossed.Add(joint.crossed);
}

void Attach::chainEnds()
{
  for (int k = 0; k < n; ++k)
    for (int end = 0; end < 2; ++end)
    {
      const Piece&  piece = pieces[k];
      const Rail&   rail  = piece.rail;
      const Sheet&  sheet = piece.sheet;
      TopoDS_Vertex ends[2];
      TopExp::Vertices(TopoDS::Edge(rail.edge.Oriented(TopAbs_FORWARD)), ends[0], ends[1]);
      if (!piece.stopVertex[end].IsNull())
        ends[end] = piece.stopVertex[end];
      const auto& joint = piece.joint[end];
      if (joint.meets == Meets::Miter || joint.meets == Meets::Corner)
        continue; // below
      if (joint.meets == Meets::Fold)
      {
        if (joint.piece > k || (joint.piece == k && joint.end > end))
          fold(k, end, joint.piece, joint.end);
        continue;
      }
      if (joint.meets == Meets::Joint || joint.meets == Meets::Split)
      {
        if (joint.piece < k || (joint.piece == k && joint.end < end))
          continue; // made from its other side
        const Piece& other = pieces[joint.piece];
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
          corner[k][side][end] = corner[joint.piece][to][joint.end] = vertex;
          u[k][side][end]                                           = here;
          u[joint.piece][to][joint.end]                             = there;
        }
        if (joint.meets == Meets::Split
            && (splitNotch(k, joint.side)
                || splitNotch(joint.piece, joint.side ^ int(joint.swapped))))
        {
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
          GeomAPI_ProjectPointOnCurve on(BRep_Tool::Pnt(corner[k][joint.side][end]), curve, f, l);
          if (on.NbPoints() == 0 || on.LowerDistance() > rail.resolved(limit::together))
            fail(op, no_solution, "a contact misses the far edge it runs off over", {rail.edge});
          // A bounded notch can remove an interior interval, leaving two edge components.
          // Sequential cuts from the same original endpoint cannot represent that set.
          const bool hereNotched = splitNotch(k, joint.side);
          const int notchIndex = hereNotched ? k : joint.piece;
          const int notchSide = hereNotched ? joint.side : joint.side ^ int(joint.swapped);
          const auto& notchPiece = pieces[notchIndex];
          const auto& beginJoint = notchPiece.joint[0];
          const auto& finish = notchPiece.joint[1];
          int crossings = 0;
          for (const auto& p : pieces)
            crossings += p.joint[1].meets == Meets::Split && p.joint[1].crossed.IsSame(edge);
          const bool paired = refineContacts && crossings == 2
            && beginJoint.meets == Meets::Split && finish.meets == Meets::Split
            && beginJoint.side == notchSide && finish.side == notchSide
            && beginJoint.crossed.IsSame(edge) && finish.crossed.IsSame(edge)
            && pieces[beginJoint.piece].rail.s[notchSide ^ int(beginJoint.swapped)].trim.IsNull()
            && pieces[finish.piece].rail.s[notchSide ^ int(finish.swapped)].trim.IsNull()
            && !cuts.Contains(edge);
          if (!paired)
            cut(edge, keep.IsSame(a) ? z : a, corner[k][joint.side][end], on.LowerDistanceParameter());
          else
          {
            const double crossing = on.LowerDistanceParameter();
            if (const auto found = seamCrossed.find(edge.TShape().get()); found == seamCrossed.end())
              seamCrossed[edge.TShape().get()] = {corner[k][joint.side][end], crossing};
            else
            {
              const auto& [other, there] = found->second;
              const double low = std::min(crossing, there), high = std::max(crossing, there);
              if (!(low > f + tol::parameter && high < l - tol::parameter
                    && high - low > tol::parameter))
                fail(op, unsupported_geometry, "a contact interval does not lie inside its crossed edge", {edge});
              const bool lower = crossing < there;
              TopoDS_Vertex e0, e1;
              TopExp::Vertices(TopoDS::Edge(edge.Oriented(TopAbs_FORWARD)), e0, e1);
              TopoDS_Compound parts;
              builder.MakeCompound(parts);
              builder.Add(parts, trimmed(edge, e1, lower ? corner[k][joint.side][end] : other, low));
              builder.Add(parts, trimmed(edge, e0, lower ? other : corner[k][joint.side][end], high));
              for (const auto& face : facesOfEdge.FindFromKey(edge))
                editOf(edits, face).replace.Bind(edge, parts);
            }
          }
        }
        else if (joint.meets == Meets::Split)
        {
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
          GeomAPI_ProjectPointOnCurve on(BRep_Tool::Pnt(corner[k][joint.side][end]), curve, f, l);
          if (on.NbPoints() == 0 || on.LowerDistance() > rail.resolved(limit::together))
            fail(op,
                 no_solution,
                 "a contact misses the edge it crosses onto the next face",
                 {rail.edge});
          const double crossing = on.LowerDistanceParameter();
          // A sphere's seam the contact crosses twice (across the cap, from rim to rim) keeps
          // its part between: each crossing cuts it from the end it is nearer.
          int crossings = 0;
          for (const auto& p : pieces)
            crossings += p.joint[1].meets == Meets::Split && p.joint[1].crossed.IsSame(edge);
          // Two crossings retain the interval between them when the edge reaches
          // the blended chain. Its two ends need not both belong to that chain.
          double otherCrossing = crossing;
          if (refineContacts && crossings == 2 && (atRail(a) || atRail(z)))
            for (int m = 0; m < n; ++m)
            {
              const auto& joint = pieces[m].joint[1];
              if (joint.meets != Meets::Split || !joint.crossed.IsSame(edge)
                  || (m == k && end == 1) || (joint.piece == k && joint.end == end))
                continue;
              const auto& q = pieces[m];
              const gp_Pnt point = q.sheet.contact(joint.side)->Value(q.sheet.rail(q.uAt(1)));
              GeomAPI_ProjectPointOnCurve other(point, curve, f, l);
              if (other.NbPoints() && other.LowerDistance() <= rail.resolved(limit::together)
                  && std::abs(other.LowerDistanceParameter() - crossing) > tol::parameter)
                otherCrossing = other.LowerDistanceParameter();
            }
          const TopoDS_Vertex from = refineContacts && otherCrossing != crossing
                                       ? (crossing < otherCrossing ? a : z)
                                     : atRail(a) && atRail(z)
                                       ? (crossings > 1 ? (crossing - f < l - crossing ? a : z)
                                          : atOwn(z) && !atOwn(a) ? z : a)
                                     : atRail(a) ? a
                                     : atRail(z) ? z
                                                 : TopoDS_Vertex();
          if (from.IsNull() && crossings == 2)
          {
            // An edge (or a seam) the contact crosses and comes back over, in the strip the
            // blend replaces between the crossings and on the faces past them: it is left in two.
            if (const auto found = seamCrossed.find(edge.TShape().get());
                found == seamCrossed.end())
              seamCrossed[edge.TShape().get()] = {corner[k][joint.side][end], crossing};
            else
            {
              const auto& [other, there] = found->second;
              const bool      lower      = crossing < there;
              TopoDS_Compound parts;
              builder.MakeCompound(parts);
              builder.Add(parts,
                          trimmed(edge,
                                  z,
                                  lower ? corner[k][joint.side][end] : other,
                                  std::min(crossing, there)));
              builder.Add(parts,
                          trimmed(edge,
                                  a,
                                  lower ? other : corner[k][joint.side][end],
                                  std::max(crossing, there)));
              for (const auto& face : facesOfEdge.FindFromKey(edge))
                editOf(edits, face).replace.Bind(edge, parts);
            }
            continue;
          }
          if (from.IsNull())
            fail(op,
                 unsupported_geometry,
                 "a contact crosses onto the next face over an edge that does not reach the "
                 "blended edges, which the blend does not do yet",
                 {rail.edge});
          cut(edge, from, corner[k][joint.side][end], crossing);
        }
        else
        {
          // The other edges at the vertex lie between faces of the supports and end there:
          // each is cut back to the corner on its side.
          for (const auto& shape : edgesOfVertex.FindFromKey(ends[end]))
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
                  return (s.cliffCurve || !s.trim.IsNull()) && s.cliffWhole
                         && onFace(shape, s.face);
                }))
            {
              swept.push_back(shape);
              continue;
            }
            int side = -1;
            for (int s = 0; s < 2; ++s)
              if (onFace(shape, rail.s[s].face)
                  || onFace(shape, other.rail.s[s ^ int(joint.swapped)].face))
                side = s;
            if (side < 0)
              fail(op,
                   unsupported_geometry,
                   "where a chain of blended edges meets another face, the corner is blended too, "
                   "which the blend does not do yet",
                   {rail.edge});
            if (poleInStrip(piece, side, TopoDS::Edge(shape))
                && sweepToPole(TopoDS::Edge(shape), ends[end]))
              continue;
            double                      f, l;
            const auto                  curve = BRep_Tool::Curve(TopoDS::Edge(shape), f, l);
            GeomAPI_ProjectPointOnCurve on(BRep_Tool::Pnt(corner[k][side][end]), curve, f, l);
            if (on.NbPoints() == 0 || on.LowerDistance() > rail.resolved(limit::together))
            {
              if (sweepToPole(TopoDS::Edge(shape), ends[end]))
                continue;
              // A body pinched at the joint (a hole tangent to a support there, which meets it in
              // a line of four faces): no blend passes it (the reference modeller refuses it too).
              for (const auto& near : edgesOfVertex.FindFromKey(ends[end]))
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
            cut(shape, ends[end], corner[k][side][end], on.LowerDistanceParameter());
          }
        }
        // The section between them: the sheet's own there, laid on the other sheet too,
        // or on its own sheet's other end where the piece closes on itself.
        // A split's parts are one face: no section between them.
        if (joint.meets == Meets::Split)
          continue;
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
                                          corner[k][0][end],
                                          reach[0],
                                          corner[k][1][end],
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
            segment2d(other.sheet.uv(other.sheet.rail(there),
                                     other.sheet.across[atLo ^ int(joint.swapped)]),
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
        continue;
      }
      // A face stops the edge here.
      stop(k, end, ends[end]);
    }
}

//! The loop where the contact on the cap folds over itself, between piece `a`'s end `ea` and
//! piece `b`'s end `eb`: both end at the place where the contact crosses itself (the one ball
//! there), each in its own section of that ball down to the other contact, and the other contact's
//! run between them is an edge of its face, which the patch over the loop meets. Within one edge
//! the run is the fit's, on either part's sheet alike; across a joint of a chain it runs on each
//! piece's sheet to the joint's place, where the edges at the joint's vertex between the faces it
//! runs on end, cut back to it as at any joint.
void Attach::fold(int a, int ea, int b, int eb)
{
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
  const double s0             = std::min(sheet.rail(refineContacts ? piece.uAt(0) : 0),
                                                    sheet.rail(refineContacts ? piece.uAt(1) : 1)),
               s1             = std::max(sheet.rail(refineContacts ? piece.uAt(0) : 0),
                                         sheet.rail(refineContacts ? piece.uAt(1) : 1));
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
      laid.push_back({f.run[i].edge, interpolate.Curve()});
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
void Attach::stop(int k, int end, const TopoDS_Vertex& vertex)
{
  const Piece& piece = pieces[k];
  const Rail&  rail  = piece.rail;
  const Sheet& sheet = piece.sheet;
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
  TopoDS_Face stop;
  TopoDS_Edge beside[2], sliverEdge[2];
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
  TopoDS_Vertex from[2] = {vertex, vertex};
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
  const double middle = piece.uAt(0) / 2 + piece.uAt(1) / 2;
  const int    cliff  = rail.s[0].cliffCurve ? 0 : rail.s[1].cliffCurve ? 1 : -1;
  // A contact that leaves its support for the stop face before the end's half of the piece (a
  // support that narrows towards the stop face, as a triangular top to its oblique end) crosses
  // the edge between them in the other half.
  const auto early = [&](int side) {
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
  };
  {
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
    if (!simple)
    {
      if (!consumed(k, end, vertex, stop, beside))
        stops(k, end, vertex);
      return;
    }
  }
  stopAt[k][end] = stop;
  gp_Pnt2d cornerOnStop2d[2];
  for (int side = 0; side < 2; ++side)
  {
    if (!rail.s[side].trim.IsNull())
    {
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
        capNotch(k, side, end, vertex, stop, cornerOnStop2d[side]);
        continue;
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
      cornerOnStop2d[side] =
        BRep_Tool::CurveOnSurface(down, stop, g0, l0)->Value(on.LowerDistanceParameter());
      continue;
    }
    if (side == cliff)
    {
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
        cornerOnStop2d[side] = BRep_Tool::CurveOnSurface(beside[side], stop, f0, l0)
                                 ->Value(BRep_Tool::Parameter(corner[k][side][end], beside[side]));
        continue;
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
      cornerOnStop2d[side] = onStop->Value(0);
      continue;
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
      found = early(side);
    if (!found)
    {
      if (!more)
        fail(op,
             unsupported_geometry,
             "the blend runs on past the face that stops it, where its supports do not meet the "
             "stop face's surface",
             {rail.edge});
      // On the stop face's surface itself where the meeting it was found on, true to the
      // intersector's tolerance, leaves it further off than a corner may be.
      double     place         = more->u;
      const auto onStopSurface = [&](double at) {
        GeomAPI_ProjectPointOnSurf project(sheet.contact(side)->Value(sheet.rail(at)),
                                           BRep_Tool::Surface(stop));
        return project.NbPoints() ? project.LowerDistance() : HUGE_VAL;
      };
      if (onStopSurface(place) > tolerance)
      {
        const double at    = sheet.rail(more->u),
                     reach = 1e-3 * std::abs(sheet.rail(1) - sheet.rail(0));
        GeomAPI_IntCS onStop(new Geom_TrimmedCurve(sheet.contact(side), at - reach, at + reach),
                             BRep_Tool::Surface(stop));
        for (int i = 1; onStop.IsDone() && i <= onStop.NbPoints(); ++i)
        {
          double su, sv, w;
          onStop.Parameters(i, su, sv, w);
          const double q = (w - sheet.alpha) / sheet.beta;
          if (i == 1 || std::abs(q - more->u) < std::abs(place - more->u))
            place = q;
        }
      }
      u[k][side][end]                   = place;
      const gp_Pnt                there = sheet.contact(side)->Value(sheet.rail(place));
      GeomAPI_ProjectPointOnCurve onMeeting(there, more->curve);
      if (place != more->u && onMeeting.NbPoints() > 0)
        more->from = onMeeting.LowerDistanceParameter();
      // A hair over the gap: the check reads the vertex against its tolerance to rounding.
      builder.MakeVertex(
        corner[k][side][end],
        there,
        std::max(tolerance, 1.01 * there.Distance(more->curve->Value(more->from))));
      // The new piece of the edge between the support and the stop face, from the
      // corner back to the vertex, in both faces' wires.
      const auto [supportSurface, supportPlace] = surfaceOf(rail.s[side].face);
      const auto [stopSurface, stopPlace]       = surfaceOf(stop);
      const TopoDS_Edge piece                   = edgeBetween(
        more->curve,
        corner[k][side][end],
        more->from,
        vertex,
        more->to,
        {{more->onSupport, supportSurface, supportPlace}, {more->onStop, stopSurface, stopPlace}},
        tolerance);
      editOf(edits, rail.s[side].face).bridges.push_back(piece);
      editOf(edits, stop).bridges.push_back(piece);
      continued.push_back({piece, {rail.s[side].face, stop}});
      cornerOnStop2d[side] = more->onStop->Value(more->from);
      continue;
    }
    u[k][side][end] = found->first;
    builder.MakeVertex(corner[k][side][end],
                       sheet.contact(side)->Value(sheet.rail(found->first)),
                       tolerance);
    cut(beside[side], from[side], corner[k][side][end], found->second);
    {
      double f0, l0;
      cornerOnStop2d[side] =
        BRep_Tool::CurveOnSurface(beside[side], stop, f0, l0)->Value(found->second);
    }
  }
  const auto [surface, location] = surfaceOf(stop);
  const auto     placed          = BRep_Tool::Surface(stop);
  const gp_Pnt2d cornerOnStop    = cornerOnStop2d[0];
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
    if (!met && refineContacts)
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
      corner[k][sa][ea] = corner[j][sb][eb] = vq;
      u[k][sa][ea]                          = q->first;
      u[j][sb][eb]                          = q->second;
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
            u[k][sa][ea] = alternative.first;
            u[j][sb][eb] = alternative.second;
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
        corner[k][oa][ea] = corner[j][ob][eb] = vp;
        u[k][oa][ea]                          = pa->first;
        u[j][ob][eb]                          = pb->first;
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
      corner[p[w]][far[w]][at[w]] = vx;
      u[p[w]][far[w]][at[w]]      = crosses[w].first;
      corner[p[l]][far[l]][at[l]] = vr;
      u[p[l]][far[l]][at[l]]      = *r;
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
  corner[k][oa][ea] = corner[j][ob][eb] = roll.X;
  corner[k][sa][ea]                     = ta;
  corner[j][sb][eb]                     = tb;
  u[k][oa][ea] = u[k][sa][ea] = pa.first;
  u[j][ob][eb] = u[j][sb][eb] = pb.first;
  cut(third, vertex, roll.X, pa.second);
  for (const auto& [m, end, at] : {std::tuple{k, ea, pa.first}, std::tuple{j, eb, pb.first}})
  {
    const Sheet& sheet = pieces[m].sheet;
    section[m][end]    = isoSection(sheet, sheet.rail(at), corner[m][0][end], corner[m][1][end]);
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
    builder.Add(wire, from(roll.arc, corner[roll.a][roll.sa][roll.ea]));
    builder.Add(wire, from(sb, corner[roll.b][roll.sb][roll.eb]));
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
void Attach::stops(int k, int end, const TopoDS_Vertex& vertex)
{
  const Piece& piece   = pieces[k];
  const Rail&  rail    = piece.rail;
  const Sheet& sheet   = piece.sheet;
  const auto   support = [&](const TopoDS_Shape& face) {
    return face.IsSame(rail.s[0].face) || face.IsSame(rail.s[1].face);
  };
  const auto across = [&](const TopoDS_Shape& edge, const TopoDS_Shape& face) {
    for (const auto& other : facesOfEdge.FindFromKey(edge))
      if (!other.IsSame(face))
        return TopoDS::Face(other);
    return TopoDS_Face();
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
  // The faces around `at` from `first`, the face beyond `firstEdge`, to the face beyond `lastEdge`,
  // and the edges between them. The next edge is the one the face's boundary passes `at` with
  // after the previous one: where it reaches `at` in the face's parameters. A periodic face's
  // seam there leads into the face again, from its other side (a hole whose seam runs from the
  // vertex). Each edge between is kept as the face after it uses it.
  const auto around = [&](const TopoDS_Vertex&      at,
                          const TopoDS_Face&        first,
                          const TopoDS_Edge&        firstEdge,
                          const TopoDS_Edge&        lastEdge,
                          std::vector<TopoDS_Face>& faces,
                          std::vector<TopoDS_Edge>& edges) {
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
      const double apart =
        2
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
      current         = seam ? current : across(edge, current);
      if (current.IsNull() || support(current))
        wrong();
      previous = seam ? TopoDS::Edge(edge.Reversed()) : edge;
      edges.push_back(previous);
      faces.push_back(current);
    }
  };
  const double middle = piece.uAt(0) / 2 + piece.uAt(1) / 2;

  struct Walk
  {
    std::vector<TopoDS_Face>   faces;   // out from the vertex
    std::vector<TopoDS_Edge>   between; // between consecutive faces
    std::vector<TopoDS_Vertex> from;    // each one's end that the blend replaces
    std::vector<TopoDS_Edge>   swept;
    std::vector<TopoDS_Vertex> sweptAt;
    TopoDS_Edge                crossed;
    TopoDS_Vertex              crossedFrom;
    Meeting                    at;
    // Where the support goes on past the vertex: the edge at the vertex between the face it
    // goes on into and the stop face, and where the contact meets the stop face's surface.
    TopoDS_Edge                 up;
    std::optional<Continuation> continued;
    // Where the crossed edge's face across goes on from the support: that face.
    TopoDS_Face level;
  } walk[2];

  for (int side = 0; side < 2; ++side)
  {
    Walk&              w   = walk[side];
    const TopoDS_Face& own = rail.s[side].face;
    TopoDS_Edge        edge;
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
    for (int step = 0;; ++step)
    {
      if (step > ofOwn.Extent())
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
      if (!found && step == 0)
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
        if (!stopFace.IsNull() && !support(stopFace))
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
      w.swept.push_back(edge);
      TopoDS_Vertex a, z;
      TopExp::Vertices(edge, a, z);
      const TopoDS_Vertex at = a.IsSame(from) ? z : a;
      w.sweptAt.push_back(at);
      TopoDS_Edge next;
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
      if (next.IsNull() || railEdge(next))
        refuse("the blend's contact does not reach the faces that stop it");
      const TopoDS_Face nextBeyond = across(next, own);
      if (nextBeyond.IsNull() || support(nextBeyond))
        refuse("the faces at an end of the edge do not stop it");
      const size_t before = w.between.size();
      facesBefore         = w.faces.size();
      betweenBefore       = before;
      around(at, w.faces.back(), edge, next, w.faces, w.between);
      w.from.insert(w.from.end(), w.between.size() - before, at);
      edge = next;
      from = at;
    }
  }
  // All the stop faces in order from side 0's crossing to side 1's.
  std::vector<TopoDS_Face> faces(walk[0].faces.rbegin(), walk[0].faces.rend());
  std::vector<TopoDS_Edge> between(walk[0].between.rbegin(), walk[0].between.rend());
  // Side 0's walk ran the other way: each seam among its edges is used from its other side.
  for (size_t i = 0; i < between.size(); ++i)
    if (BRep_Tool::IsClosed(between[i], faces[i]))
      between[i].Reverse();
  std::vector<TopoDS_Vertex> from(walk[0].from.rbegin(), walk[0].from.rend());
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
    around(vertex, walk[0].faces.front(), beside[0], beside[1], ring, edges);
    faces.insert(faces.end(), ring.begin(), ring.end());
    between.insert(between.end(), edges.begin(), edges.end());
    from.insert(from.end(), edges.size(), vertex);
  }
  faces.insert(faces.end(), walk[1].faces.begin() + 1, walk[1].faces.end());
  between.insert(between.end(), walk[1].between.begin(), walk[1].between.end());
  from.insert(from.end(), walk[1].from.begin(), walk[1].from.end());
  if (faces.size() != between.size() + 1)
    refuse("the faces that stop the blend do not run from one support to the other");
  // The places: the two corners, and where the sheet crosses each edge between stop faces.
  std::vector<TopoDS_Vertex> at(faces.size() + 1);
  std::vector<double>        param(faces.size() + 1);
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
  std::vector<bool> touched(faces.size(), false);
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
  // Along the sheet on each stop face, from one place to the next.
  std::vector<TopoDS_Edge> made(faces.size());
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

void Attach::cornerPatches()
{
  // Three or more edges that meet: the rolling ball's corner. At equal distances the round's
  // corner is the ball itself where one rests on all the faces (three always have one; more need
  // not).
  for (size_t c = 0; c < cornersAt.size(); ++c)
  {
    const auto&  at    = cornersAt[c];
    const Piece& first = pieces[at.ends[0].first];
    if (at.pinch)
    {
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
      continue;
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
        continue;
      }
      if (at.ends.size() != 3)
        fail(op,
             unsupported_geometry,
             "where more than three G2, conic or chamfer blends meet, the corner is blended with "
             "them, which the blend does not do yet",
             {first.rail.edge});
      patches[c] = lawCorner(at, c);
      continue;
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
      if (railStationSeparation(piece.rail, ends[0], ends[1])
          > piece.rail.resolved(limit::together))
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
    // The sphere's frame. Where the edges are square to each other, as OpenCascade's own corners
    // are made, its pole is the ball's foot on the face across from the first edge and its x
    // axis an end of that edge's section: that section is its equator and the other two are
    // meridians, every side an isoline. Elsewhere the poles and seam are kept away from the patch.
    gp_Vec e[3];
    bool   square = count == 3;
    for (int i = 0; i < 3 && square; ++i)
      square &= leaving(pieces[at.ends[i].first].rail.edge, at.vertex, e[i], at.ends[i].second);
    square = square && std::abs(e[0].Dot(e[1])) < tol::fineAngle
             && std::abs(e[0].Dot(e[2])) < tol::fineAngle
             && std::abs(e[1].Dot(e[2])) < tol::fineAngle;
    gp_Ax3 frame;
    int    poleFoot = -1;
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
void Attach::splitCorner(const Corner& at, size_t c)
{
  using Homogeneous = std::array<double, 6>; // w (x - V), w s, w r, w: on the corner and a sheet
  constexpr int degree = 6;
  const Piece&  first  = pieces[at.ends[0].first];
  const int     count  = int(at.ends.size());
  const double  scale  = first.rail.scale;
  const gp_XYZ  V      = BRep_Tool::Pnt(at.vertex).XYZ();
  const char*   notYet = "where more than three conic blends meet, the blend makes the corner "
                         "only between planes so far";
  std::vector<TopoDS_Vertex> X(count); // X[i]: where the contacts of blends i and i + 1 cross
  for (int i = 0; i < count; ++i)
    X[i] = contactsCross(at, i, (i + 1) % count).first;
  std::vector<gp_Vec>                    d(count);
  std::vector<Handle(Geom_BSplineCurve)> profile(count);
  for (int i = 0; i < count; ++i)
  {
    const Piece& piece = pieces[at.ends[i].first];
    const auto   extrusion =
      Handle(Geom_SurfaceOfLinearExtrusion)::DownCast(piece.sheet.surface);
    if (!extrusion.IsNull())
      profile[i] = Handle(Geom_BSplineCurve)::DownCast(extrusion->BasisCurve());
    if (!piece.sheet.swapped || profile[i].IsNull() || profile[i]->Degree() != 2
        || profile[i]->NbPoles() != 3)
      fail(op, unsupported_geometry, notYet, {piece.rail.edge});
    d[i] = gp_Vec(extrusion->Direction());
  }
  const auto point = [&](const Homogeneous& h) {
    return gp_XYZ(h[0], h[1], h[2]) / h[5];
  };
  // de Casteljau's split of a homogeneous Bezier curve at t.
  const auto split = [](std::vector<Homogeneous> b, double t) {
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
  };
  // A cubic half in the form of weights 1, 1, ..., 1, 1: its weights' Moebius standard form, then
  // the cubic factor (1, 2 - w1, 2 - w2, 1) of numerator and denominator. Empty where a weight is
  // not positive.
  const auto sextic = [&](std::vector<Homogeneous> b) -> std::optional<std::vector<Homogeneous>> {
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
    const double             g[4]    = {1, 2 - b[1][5], 2 - b[2][5], 1};
    constexpr double         c3[4]   = {1, 3, 3, 1};
    constexpr double         c6[7]   = {1, 6, 15, 20, 15, 6, 1};
    std::vector<Homogeneous> product(degree + 1, Homogeneous{});
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        for (int x = 0; x < 6; ++x)
          product[i + j][x] += c3[i] * c3[j] / c6[i + j] * g[j] * b[i][x];
    for (const auto& pole : product)
      if (!(pole[5] > 0))
        return {};
    return product;
  };
  // Each blend's end: half[i][0] from X[i - 1] to M_i, half[i][1] from X[i] to M_i.
  std::vector<std::array<std::vector<Homogeneous>, 2>> half(count);
  std::vector<std::vector<Homogeneous>>                 cubics(count);
  std::vector<gp_Vec> outwardAt(count); // out of the face blends i and i + 1 share, at X[i]
  for (int i = 0; i < count; ++i)
  {
    const auto& [k, end] = at.ends[i];
    const Piece& piece   = pieces[k];
    const int    prev = (i + count - 1) % count, next = (i + 1) % count;
    const int    a    = corner[k][0][end].IsSame(X[prev]) ? 0 : 1;
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
    const auto&  conic    = profile[i];
    const bool   reversed = std::abs(from.X() - conic->LastParameter())
                          < std::abs(from.X() - conic->FirstParameter());
    gp_XYZ       p[3];
    double       w[3];
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
    const double rhoA = slope(start, d[prev]), rhoB = slope(stop, d[next]);
    const double R0 = from.Y() * w[0], R3 = to.Y() * w[2];
    const double R[4] = {R0,
                         R0 + (rhoA * w[0] + from.Y() * 2 * (w[1] - w[0])) / 3,
                         R3 - (rhoB * w[2] + to.Y() * 2 * (w[2] - w[1])) / 3,
                         R3};
    const double s0 = from.X(), s1 = to.X();
    const gp_XYZ hp[4] = {p[0] * w[0],
                          (p[0] * w[0] + p[1] * (2 * w[1])) / 3,
                          (p[1] * (2 * w[1]) + p[2] * w[2]) / 3,
                          p[2] * w[2]};
    const double hw[4] = {w[0], (w[0] + 2 * w[1]) / 3, (2 * w[1] + w[2]) / 3, w[2]};
    const double hs[4] = {s0 * w[0],
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
  // Blend i's end split where its halves' speeds at M_i agree: the first's grows from nothing as
  // the split moves along, the second's shrinks to nothing.
  const auto splitEnd = [&](int i) {
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
    const TopoDS_Edge& edge = pieces[at.ends[i].first].rail.edge;
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
  };
  const auto splitPoint = [&](int i) { return point(half[i][0][degree]); };
  for (int i = 0; i < count; ++i)
    splitEnd(i);
  // The offsets across each half end, along its blend's edge: l0 at the crossing (the other
  // end's first leg), l_i at M_i, and l_i (1 - cos(2 pi / n) / 6) next to M_i, where the
  // condition across the inner curve then holds at its second coefficient; between, from one to
  // the other in even steps. l_i is a sixth of the chord from M_i to the splits' mean: the speed
  // of the inner curve's first leg that of a straight curve to the centre's likely place.
  const double        beta = 2 * std::cos(2 * M_PI / count);
  std::vector<double> legRow(count), legCol(count), ell(count);
  for (int i = 0; i < count; ++i)
  {
    const int    next = (i + 1) % count;
    const gp_Vec alongRow(point(half[i][1][1]) - point(half[i][1][0]));
    const gp_Vec alongCol(point(half[next][0][1]) - point(half[next][0][0]));
    legRow[i] = alongRow.Dot(d[next]);
    legCol[i] = alongCol.Dot(d[i]);
    if (alongRow.Crossed(d[next]).Magnitude() > tol::fineAngle * alongRow.Magnitude()
        || alongCol.Crossed(d[i]).Magnitude() > tol::fineAngle * alongCol.Magnitude())
      fail(op, unsupported_geometry, notYet, {first.rail.edge});
  }
  gp_XYZ centred;
  for (int i = 0; i < count; ++i)
    centred += splitPoint(i) / count;
  for (int i = 0; i < count; ++i)
  {
    // Both patches beside blend i's end lie on one side of it: towards the corner.
    if (!(legRow[(i + count - 1) % count] * legCol[i] > 0))
      fail(op,
           no_solution,
           "the ends of the blends at a corner turn away from it",
           {pieces[at.ends[i].first].rail.edge});
    ell[i] = (legCol[i] > 0 ? 1 : -1) * (centred - splitPoint(i)).Modulus() / degree;
  }
  const auto offset = [&](double l0, double l, int j) {
    const double inner = l * (1 - beta / (2 * degree));
    return j <= 1             ? l0
           : j == degree      ? l
           : j == degree - 1  ? inner
                              : l0 + (inner - l0) * (j - 1) / 4;
  };
  // Patch q at X[q]: u along blend q + 1's end to M_{q+1} (its side v = 0), v along blend q's to
  // M_q (u = 0); the inner curve q + 1 is its side u = 1, curve q its side v = 1. A pole is fixed
  // or one of the unknowns: per patch its 16 inner ones, per curve 4, and the centre.
  struct Ref
  {
    int    var = -1;
    gp_XYZ at;
  };
  const int  curveBase = 16 * count, centre = 20 * count, unknowns = 20 * count + 1;
  const auto onCurve   = [&](int i, int t) -> Ref {
    if (t <= 1)
      return {-1, splitPoint(i) + d[i].XYZ() * (t * ell[i])};
    return {t == degree ? centre : curveBase + 4 * i + t - 2, {}};
  };
  const auto pole = [&](int q, int j, int k) -> Ref {
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
  };
  const auto weight = [&](int q, int j, int k) {
    if (k <= 1 && j < degree)
      return half[(q + 1) % count][0][j][5];
    if (j <= 1 && k < degree)
      return half[q][1][k][5];
    return 1.0;
  };
  // The fairness, the thin-plate energy of the nets: their second differences along u and v and
  // their twists. The conditions across the inner curves, coefficient by coefficient.
  struct Term
  {
    double coefficient;
    Ref    ref;
  };
  std::vector<std::vector<Term>> stencils;
  for (int q = 0; q < count; ++q)
    for (int j = 0; j <= degree; ++j)
      for (int k = 0; k <= degree; ++k)
      {
        if (j > 0 && j < degree)
          stencils.push_back({{1, pole(q, j - 1, k)}, {-2, pole(q, j, k)}, {1, pole(q, j + 1, k)}});
        if (k > 0 && k < degree)
          stencils.push_back({{1, pole(q, j, k - 1)}, {-2, pole(q, j, k)}, {1, pole(q, j, k + 1)}});
        if (j < degree && k < degree)
          stencils.push_back({{M_SQRT2, pole(q, j + 1, k + 1)},
                              {-M_SQRT2, pole(q, j + 1, k)},
                              {-M_SQRT2, pole(q, j, k + 1)},
                              {M_SQRT2, pole(q, j, k)}});
      }
  std::vector<std::vector<Term>> rows;
  for (int i = 0; i < count; ++i)
  {
    const int a = (i + count - 1) % count, b = i;
    for (int t = 0; t <= degree; ++t)
      rows.push_back({{2 * degree - beta * t, onCurve(i, t)},
                      {beta * t, t > 0 ? onCurve(i, t - 1) : Ref{}},
                      {-double(degree), pole(a, degree - 1, t)},
                      {-double(degree), pole(b, t, degree - 1)}});
  }
  // One system for the three coordinates, [H A'; A 0] [x; y] = [-g; b], with one more row that
  // holds the centre at a point c; its constraints are dependent (the centre's take two of them),
  // so it is solved by singular values. Its solution is x0 + g c, g the same for every coordinate.
  std::vector<std::vector<Term>> conditions;
  std::copy_if(rows.begin(), rows.end(), std::back_inserter(conditions), [](const auto& row) {
    return std::any_of(row.begin(), row.end(), [](const Term& term) { return term.ref.var >= 0; });
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
  rhs[3](size)                         = 1;
  math_SVD                 svd(K);
  std::vector<math_Vector> solution(4, math_Vector(1, size));
  for (int x = 0; x < 4; ++x)
    svd.Solve(rhs[x], solution[x], 1e-12);
  const math_Vector& g        = solution[3];
  const auto         valueFor = [&](const Ref& ref, const gp_XYZ& c) {
    if (ref.var < 0)
      return ref.at;
    const int x = ref.var + 1;
    return gp_XYZ(solution[0](x), solution[1](x), solution[2](x)) + c * g(x);
  };
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
  // The centre itself where the corner bends least: the least integral of k1^2 + k2^2 over the
  // patches, the bending energy of the surface (not of its nets). The nets' energy is a stand-in
  // for it only where the patches are parametrised evenly; uneven corners have patches of very
  // different shapes, and it pulls the centre towards the long ones. Every other pole stays the
  // nets' fairest for that centre. The energy is a sum of squares of the shape operator's entries
  // in an orthonormal frame at Gauss points, minimised from the nets' fairest centre by
  // Levenberg-Marquardt steps until no step lowers it.
  const gp_XYZ centrePoint = [&] {
    // A composite Gauss rule, 4 x 4 squares of 8 x 8 points per patch. The energy's integrand
    // varies sharply where a patch is narrow (one Gauss rule of 10 x 10 points misses it by a
    // quarter); on the uneven family the centre this rule finds lies within 0.7 % of the corner's
    // size of the one 8 x 8 squares of 6 x 6 points find (measured).
    constexpr int cells = 4, nodes = 8, points = cells * nodes;
    const math_ComputeGaussPointsAndWeights gauss(nodes);
    const math_Vector gx = gauss.Points(), gw = gauss.Weights();
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
          net[j][k]        = {valueFor(ref, gp_XYZ()) * w, ref.var < 0 ? 0.0 : g(ref.var + 1) * w, w};
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
    const auto residuals = [&](const gp_XYZ& c, std::vector<double>& r) {
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
      std::vector<double> plus;
      math_Vector         grad(1, 3, 0.0);
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
      sv(1) = move.X();
      sv(2) = move.Y();
      sv(3) = move.Z();
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
  }();
  const auto value = [&](const Ref& ref) { return valueFor(ref, centrePoint); };
  // The conditions hold: the first two by construction, the others as solved.
  const double exact = 1e-9 * scale;
  for (const auto& row : rows)
  {
    gp_XYZ sum;
    for (const auto& term : row)
      sum += value(term.ref) * term.coefficient;
    if (sum.Modulus() > exact)
      fail(op, no_solution, "the corner's split does not join smoothly", {first.rail.edge});
  }
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
  std::vector<std::array<TopoDS_Edge, 2>> ends(count); // [blend][half]
  std::vector<TopoDS_Edge>                inner(count);  // inner curve i, from M_i to C
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
      ends[i][h] = edgeBetween(
        curve3d(at3, w),
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
    const int prev       = (i + count - 1) % count;
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
          const double ka = within(one.fitted.contact[sk], std::min(u[k][sk][0], u[k][sk][1])),
                       kb = within(one.fitted.contact[sk], std::max(u[k][sk][0], u[k][sk][1]));
          const double ma = within(two.fitted.contact[sm], std::min(u[m][sm][0], u[m][sm][1])),
                       mb = within(two.fitted.contact[sm], std::max(u[m][sm][0], u[m][sm][1]));
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
                 {corner[k][sk][0], corner[k][sk][1], corner[m][sm][0], corner[m][sm][1]})
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
    {
      if (!contact[k][side].IsNull())
        continue; // one edge for parts of a split already made
      if (vanished(k))
        continue;
      const Piece& piece = pieces[k];
      // A part that a split's earlier part runs on into on this face is made with that one's.
      if (const auto& before = piece.joint[0];
          before.meets == Meets::Split
          && pieces[before.piece].rail.s[side].face.IsSame(piece.rail.s[side].face)
          && !BRep_Tool::IsClosed(before.crossed, piece.rail.s[side].face)
          && piece.rail.s[side].trim.IsNull() && pieces[before.piece].rail.s[side].trim.IsNull())
        continue;
      const Sheet& sheet = piece.sheet;
      if (const Side& notch = piece.rail.s[side]; !notch.trim.IsNull())
      {
        // A notch: the blend's edge on the face beyond, which takes it in the far edge's
        // place; the narrow face goes whole, with its far edge.
        const auto&  trim = trimOf(k, side);
        const double s0 = sheet.rail(u[k][side][0]), s1 = sheet.rail(u[k][side][1]);
        const auto   closed = std::find_if(rolled.begin(), rolled.end(), [&](const auto& r) {
          return r.first.first == k;
        });
        if (closed != rolled.end())
          contact[k][side] = rolledEdge(k, side, closed->second);
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
          contact[k][side] =
            edgeBetween(trim.curve,
                        corner[k][side][0],
                        s0,
                        corner[k][side][1],
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
            .replace.Bind(notch.trimEdge, sameWay(contact[k][side], notch.trimEdge));
          removed.push_back(notch.face);
          swept.push_back(notch.trimEdge);
          continue;
        }
        Edit& edit = editOf(edits, notch.trim);
        edit.bridges.push_back(contact[k][side]);
        // From a split partway: the far edge is cut back to the crossing, the rest stays.
        if (splitNotch(k, side))
          continue;
        edit.replace.Bind(notch.trimEdge, TopoDS_Edge());
        const auto F0 = trimFar.find({k, side, 0}), F1 = trimFar.find({k, side, 1});
        if (notch.cliffWhole)
        {
          removed.push_back(notch.face);
          swept.push_back(notch.trimEdge);
          continue;
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
            GeomAPI_ProjectPointOnCurve(BRep_Tool::Pnt(corner[k][side][0]), farCurve, ff, fl)
              .LowerDistanceParameter();
          const double w1 =
            GeomAPI_ProjectPointOnCurve(BRep_Tool::Pnt(corner[k][side][1]), farCurve, ff, fl)
              .LowerDistanceParameter();
          const double      inside = (w0 + w1) / 2;
          const TopoDS_Edge part =
            inside > w ? trimmed(forward, v1, at->second, w) : trimmed(forward, v0, at->second, w);
          onNarrow.bridges.push_back(part);
          edit.bridges.push_back(part);
        }
        swept.push_back(notch.trimEdge);
        continue;
      }
      if (const Side& cliff = piece.rail.s[side]; cliff.cliffCurve)
      {
        // A cliff: its far edge, now between the blend and the face beyond, which keeps
        // it; the narrow face goes whole.
        const TopoDS_Face beyond       = beyondCliff(cliff);
        const auto [surface, location] = surfaceOf(beyond);
        const double s0 = sheet.rail(u[k][side][0]), s1 = sheet.rail(u[k][side][1]);
        const auto   onBeyond = GeomProjLib::Curve2d(sheet.contact(side),
                                                   std::min(s0, s1),
                                                   std::max(s0, s1),
                                                   BRep_Tool::Surface(beyond));
        if (beyond.IsNull() || onBeyond.IsNull())
          fail(op,
               no_solution,
               "the blend's contact on the face's far edge cannot be laid on the face beyond",
               {piece.rail.edge});
        contact[k][side] =
          edgeBetween(sheet.contact(side),
                      corner[k][side][0],
                      s0,
                      corner[k][side][1],
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
          const double         w0 = at(corner[k][side][0]), w1 = at(corner[k][side][1]);
          const TopoDS_Vertex& low  = w0 < w1 ? corner[k][side][0] : corner[k][side][1];
          const TopoDS_Vertex& high = w0 < w1 ? corner[k][side][1] : corner[k][side][0];
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
          onBeyond.bridges.push_back(contact[k][side]);
          continue;
        }
        editOf(edits, beyond).replace.Bind(cliff.cliff, sameWay(contact[k][side], cliff.cliff));
        removed.push_back(cliff.face);
        continue;
      }
      const auto [support, location] = surfaceOf(piece.rail.s[side].face);
      // Where ridges take the contact's place: its parts between them, each an edge of the
      // support, which bridge the cap's gaps (two islands of it where a neck closes), and the
      // ridges, as one run along the side.
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
      if (!stops.empty())
      {
        if (piece.joint[0].meets == Meets::Split || piece.joint[1].meets == Meets::Split)
          fail(op,
               unsupported_geometry,
               "the blends of two parts of an outline cross over their cap beside a split, which "
               "the blend does not do yet",
               {piece.rail.edge});
        const bool rising = u[k][side][0] < u[k][side][1];
        const auto order  = [&](const auto& x, const auto& y) {
          return rising ? x.first < y.first : x.first > y.first;
        };
        std::sort(stops.begin(), stops.end(), order);
        std::sort(creases.begin(), creases.end(), order);
        stops.insert(stops.begin(), {u[k][side][0], corner[k][side][0]});
        stops.push_back({u[k][side][1], corner[k][side][1]});
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
        contact[k][side] = along.front();
        continue;
      }
      // The parts of a split edge that stay on this face share one contact, as their face
      // is one: from the first part's start to the last one on the face, its curve on the
      // face the parts' curves joined.
      // A seam of the face, which it crosses on itself, ends its edge there all the same.
      std::vector<int> run = {k};
      while (pieces[run.back()].joint[1].meets == Meets::Split
             && pieces[pieces[run.back()].joint[1].piece].rail.s[side].face.IsSame(
               piece.rail.s[side].face)
             && !BRep_Tool::IsClosed(pieces[run.back()].joint[1].crossed, piece.rail.s[side].face)
             && pieces[pieces[run.back()].joint[1].piece].rail.s[side].trim.IsNull())
        run.push_back(pieces[run.back()].joint[1].piece);
      const int    last = run.back();
      TopoDS_Edge& edge = contact[k][side];
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
        const double s0 = sheet.rail(u[k][side][0]), s1 = sheet.rail(u[k][side][1]);
        const auto   exact =
          GeomProjLib::Curve2d(sheet.contact(side), s0, s1, piece.rail.s[side].surface);
        if (!exact.IsNull())
          onSupport = aligned(exact, s0, onSupport->Value(s0), piece.rail.s[side].surface);
      }
      const double s0 = sheet.rail(u[k][side][0]), s1 = pieces[last].sheet.rail(u[last][side][1]);
      if (piece.rail.parameterSupports
          && !contactWithinFace(piece.rail.s[side].face, onSupport, s0, s1))
        fail(op,
             unsupported_geometry,
             "a retained blend contact lies outside its finite support",
             {piece.rail.edge});
      // Where another blend's contact on this face is this one, the face between them is
      // used up (a box rounded to half its width): the two blends meet along it.
      if (const auto same = sharedContact(piece.rail.s[side].face,
                                          corner[k][side][0],
                                          corner[last][side][1],
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
        for (int m : run)
          contact[m][side] = same;
        usedUp.Add(piece.rail.s[side].face);
        continue;
      }
      edge =
        edgeBetween(sheet.contact(side),
                    corner[k][side][0],
                    s0,
                    corner[last][side][1],
                    s1,
                    {{sheet.contactUV(side), sheet.surface, {}}, {onSupport, support, location}},
                    tolerance);
      for (int m : run)
        contact[m][side] = edge;
      // In a face beside the edge it takes the edge's place; in a face the contact crossed onto,
      // or beside a part of the edge already taken, it closes the gap it leaves.
      Edit& edit = editOf(edits, piece.rail.s[side].face);
      if (onFace(piece.rail.edge, edit.face) && !edit.replace.IsBound(piece.rail.edge))
        edit.replace.Bind(piece.rail.edge, edge);
      else
        edit.bridges.push_back(edge);
    }
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
            users += !contact[k][side].IsNull() && contact[m][t].IsSame(contact[k][side]);
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

void Attach::blendFaces()
{
  // The blend faces: around the first contact, the far end, the second contact and the near
  // end, counter-clockwise in their parameters and facing out as the supports do.
  for (int k = 0; k < n; ++k)
  {
    // One face for the parts of an edge that a split joins, made at the first of them; none
    // for a blend that has no length.
    if (pieces[k].joint[0].meets == Meets::Split || vanished(k))
      continue;
    std::vector<int> group = {k};
    while (pieces[group.back()].joint[1].meets == Meets::Split)
      group.push_back(pieces[group.back()].joint[1].piece);
    const int    last  = group.back();
    const Piece& piece = pieces[k];
    const Sheet& sheet = piece.sheet;
    builder.MakeFace(blends[k], sheet.surface, tolerance);
    const bool ccw =
      sheet.swapped ? sheet.across[1] < sheet.across[0] : sheet.across[1] > sheet.across[0];
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
      along(m, 0, corner[m][0][0]);
    across(last, 1, corner[last][0][1]);
    for (auto m = group.rbegin(); m != group.rend(); ++m)
      along(*m, 1, corner[*m][1][1]);
    // A seam's first vertex is that of both its uses, so the closed edge's are set apart.
    const size_t far =
      std::find_if(loop.begin(),
                   loop.end(),
                   [&](const TopoDS_Edge& e) { return e.IsSame(section[last][1]); })
      - loop.begin();
    if (section[k][0].IsSame(section[last][1]) && chain[k][0].empty())
      loop.push_back(TopoDS::Edge(loop[far].Reversed()));
    else
      across(k, 0, corner[k][1][0]);
    TopoDS_Wire wire;
    builder.MakeWire(wire);
    for (const auto& edge : loop)
      builder.Add(wire, ccw ? edge : TopoDS::Edge(edge.Reversed()));
    wire.Closed(true);
    builder.Add(blends[k], wire);
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
            if (m != k && contact[m][t].IsSame(contact[k][side]))
            {
              if (!bevel)
                mark(contact[k][side], blends[k], blends[m], JoinRole::UsedSupport, smooth);
              else
                record(contact[k][side], JoinRole::UsedSupport, 0);
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
          mark(edge, blends[k], piece.rail.s[side].trim.IsNull()
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
          mark(f.crease[1], blends[f.a], blends[f.b], JoinRole::Joint, smooth);
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
                           !s.trim.IsNull() ? s.trim
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
  for (const auto& piece : corner)
    for (const auto& side : piece)
      for (const auto& vertex : side)
        corners.Add(vertex);
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
            || separatedByPlane(TopoDS::Face(f.Current()), face, tolerance))
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

TopoDS_Shape Attach::assemble()
{
  // The changed faces, and the solids with them and the blend faces; other shells as they were.
  std::vector<std::pair<TopoDS_Face, std::vector<TopoDS_Face>>> faces;
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
        history->AddGenerated(pieces[k].apex, contact[k][side]);
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
                    bool refineContacts)
{
  std::vector<Piece> fragments;
  std::vector<Corner> fragmentCorners;
  if (refineContacts)
  {
    fragments = pieces;
    fragmentCorners = cornersAt;
    refineOwnership(fragments, fragmentCorners);
  }
  Attach a(std::move(body), refineContacts ? fragments : pieces,
           refineContacts ? fragmentCorners : cornersAt, history, refineContacts);
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
