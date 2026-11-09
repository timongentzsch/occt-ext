// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once

// Blend: the attachment of the blends to the body (Attach), shared by its stage files.

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
         bool refineContacts,
         StageReads* reads);

  const bool refineContacts;
  StageReads* const reads;
  //! Whether the contacts are refined at a place where refining them `changes` what is attached;
  //! where it does, a stage that does not refine them notes it (StageReads::refines).
  bool refines(bool changes) const
  {
    if (changes && !refineContacts && reads)
      reads->refines = true;
    return changes && refineContacts;
  }

  // The stages, in order.
  //! Each end of each piece: a stop face cut back to the blend's end section, a joint or a
  //! split shared with the next piece.
  void         chainEnds();
  //! The steps of chainEnds: one end of one piece, by how it meets.
  void chainEnd(int k, int end);
  //! A joint or a split shared with the next piece: the corners, the cuts and the section.
  void jointEnd(int k, int end, const TopoDS_Vertex& endVertex);
  //! A far edge a split's contact runs off over partway: cut back to the crossing.
  void splitNotchCut(int k, int end);
  //! The edge a split's contact crosses onto the next face: cut back to it.
  void splitCrossingCut(int k, int end);
  //! The other edges at a joint's vertex: each cut back to the corner on its side.
  void jointEdgeCuts(int k, int end, const TopoDS_Vertex& vertex);
  //! The section between the two pieces at a joint.
  void              jointSection(int k, int end);
  void              miters();
  void              cornerPatches();
  //! The steps of cornerPatches: one corner, by kind.
  void cornerPatch(size_t c);
  //! Every blend closes to the vertex: a degenerate end each, and no patch.
  void pinchCorner(const Corner& at);
  //! Each blend's end at the section through the corner ball's centre, from foot to foot.
  void ballSections(const Corner&                     at,
                    const std::vector<TopoDS_Face>&   faces,
                    const std::vector<TopoDS_Vertex>& foot);
  //! The corner ball's frame: its pole at a foot where the edges are square, else off the patch.
  gp_Ax3 ballFrame(const Corner&                   at,
                   const std::vector<TopoDS_Face>& faces,
                   const gp_Pnt&                   ball,
                   const std::vector<gp_Vec>&      normal,
                   double                          sigma,
                   bool&                           square,
                   int&                            poleFoot);
  //! The ball's face, bounded by the blends' sections laid on it.
  void              ballPatchFace(const Corner&                     at,
                                  size_t                            c,
                                  const std::vector<TopoDS_Vertex>& foot,
                                  const gp_Ax3&                     frame,
                                  bool                              square,
                                  int                               poleFoot);
  void              ridges();
  void              contacts();
  void              blendFaces();
  void              continuity();
  std::vector<Join> joins;
  void         tolerances();
  void         collisions();
  TopoDS_Shape assemble();

  //! The steps of contacts: one side's contact edge, by kind.
  void sideContact(int k, int side);
  //! A notch's: the blend's edge on the face beyond, in the far edge's place.
  void notchContact(int k, int side);
  //! A cliff's: on the far edge, now between the blend and the face beyond.
  void cliffContact(int k, int side);
  //! One crossing onto other faces of its support's surface: one edge on each face.
  void passedContact(int k, int side);
  //! Where ridges take the contact's place: its parts between them; false where there are none.
  bool ridgedContact(int                         k,
                     int                         side,
                     const Handle(Geom_Surface)& support,
                     const TopLoc_Location&      location);
  //! The contact on its support, shared by a split's parts on the face, or with another blend.
  void plainContact(int                         k,
                    int                         side,
                    const Handle(Geom_Surface)& support,
                    const TopLoc_Location&      location);

  void fold(int a, int ea, int b, int eb);
  void foldPatches();
  struct Folded;
  Handle(Geom_BSplineSurface) ballPatch(const Folded& f, int nv);
  void                             creaseFold(int a, int ea, int b, int eb);
  void                             crossingFold(int a, int ea, int b, int eb);
  std::optional<TopoDS_Face>  ballFace(Folded&                            f,
                                       const Handle(Geom_BSplineSurface)& surface,
                                       int                                nv);
  TopoDS_Face                 lawCorner(const Corner& at, size_t c);
  std::pair<TopoDS_Vertex, double> contactsCross(const Corner& at, int i, int j);
  void                             splitCorner(const Corner& at, size_t c);
  //! The central split's construction (defined beside splitCorner).
  struct CentralSplit;
  //! Each blend's end at the corner as a homogeneous cubic on its sheet, from crossing to crossing.
  void splitCornerEnds(const Corner& at, CentralSplit& cs);
  //! The split's patches, the ends' halves and the inner curves, all on shared edges.
  void                        splitCornerFaces(const Corner&       at,
                                               size_t              c,
                                               const CentralSplit& cs,
                                               const gp_XYZ&       centrePoint);
  TopoDS_Face                 tracedCorner(const Corner& at,
                                           size_t        c,
                                           const TopoDS_Vertex (&Q)[3][3],
                                           double crossingGap);
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

  void stop(int k, int end, const TopoDS_Vertex& vertex);

  //! A stop at one face (stop): that face, the edges beside the end on each support and the ends
  //! of those the blend replaces (a sliver's far corner), the piece's middle and its cliff side.
  struct OneStop
  {
    TopoDS_Face   stop;
    TopoDS_Edge   beside[2];
    TopoDS_Vertex from[2];
    double        middle = 0;
    int           cliff  = -1;
  };

  //! The steps of stop: the face across both sides, with the slivers it swallows.
  void stopFace(int k, int end, const TopoDS_Vertex& vertex, OneStop& one);
  //! Where a contact leaves its support for the stop face before the end's half of the piece.
  std::optional<Meeting> earlyCrossing(int k, int end, int side, const OneStop& s) const;
  //! Whether each contact meets the stop face simply: else several faces stop the blend.
  bool simpleStop(int k, int end, const OneStop& s) const;
  //! The corner of one side on the stop face, and its place in the stop face's parameters.
  void stopCorner(int                  k,
                  int                  end,
                  int                  side,
                  const TopoDS_Vertex& vertex,
                  const OneStop&       s,
                  gp_Pnt2d&            cornerOnStop);
  //! A notch's corner: where the blend's edge on the face beyond meets the stop face.
  void notchStopCorner(int                  k,
                       int                  end,
                       int                  side,
                       const TopoDS_Vertex& vertex,
                       const OneStop&       s,
                       gp_Pnt2d&            cornerOnStop);
  //! A cliff's corner: the far edge's own vertex, or where the stop face crosses the cliff.
  void cliffStopCorner(int                  k,
                       int                  end,
                       int                  side,
                       const TopoDS_Vertex& vertex,
                       const OneStop&       s,
                       gp_Pnt2d&            cornerOnStop);
  //! A corner where the contact meets the stop face's surface carried past the vertex.
  void continuedStopCorner(int                  k,
                           int                  end,
                           int                  side,
                           const TopoDS_Vertex& vertex,
                           const OneStop&       s,
                           Continuation&        more,
                           gp_Pnt2d&            cornerOnStop);
  //! The blend's end on the stop face: the sheet's section where it lies in it, else the meeting.
  void stopSection(int k, int end, const TopoDS_Face& stop, const gp_Pnt2d& cornerOnStop);

  bool        consumed(int                  k,
                       int                  end,
                       const TopoDS_Vertex& vertex,
                       const TopoDS_Face&   stop,
                       const TopoDS_Edge (&beside)[2]);
  TopoDS_Edge rolledEdge(int k, int side, const TopoDS_Face& stop);
  std::map<std::pair<int, int>, TopoDS_Face>
       rolled; // [piece][end] -> the stop face a notch's blend uses up there (consumed)
  void stops(int k, int end, const TopoDS_Vertex& vertex);

  //! A face of a support's own surface a contact crosses onto before the faces that stop it:
  //! the edge it crosses, that edge's end the blend replaces, where it crosses, the face.
  struct StopPassage
  {
    TopoDS_Edge   edge;
    TopoDS_Vertex from;
    Meeting       at;
    TopoDS_Face   face;
  };

  //! The walk out from a stop's vertex along one support to the faces that stop its contact.
  struct StopWalk
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
    // Faces of the support's own surface the contact crosses onto before it reaches the faces
    // that stop it, in the order it crosses them, each with the edge it crosses onto it, that
    // edge's end the blend replaces and where it crosses.
    std::vector<StopPassage> passed;
    // Where the contact ends on its support (or a face of its surface) without crossing an edge
    // of it: that face, and the edge across which the face before it lies, from its end the
    // blend replaces.
    TopoDS_Face   open;
    TopoDS_Edge   openEdge;
    TopoDS_Vertex openFrom;
  };

  //! A stop at several faces (stops): the piece and its end, the walks along both supports, the
  //! stop faces in order from side 0's crossing to side 1's with the edges between them and the
  //! ends of those the blend replaces, the places on them and the edges the end runs along.
  struct StopEnd
  {
    int                        k, end;
    TopoDS_Vertex              vertex;
    double                     middle     = 0;
    double                     carried[2] = {0, 0};
    StopWalk                   walk[2];
    std::vector<TopoDS_Face>   faces;
    std::vector<TopoDS_Edge>   between;
    std::vector<TopoDS_Vertex> from;
    std::vector<TopoDS_Vertex> at;
    std::vector<double>        param;
    std::vector<bool>          touched;
    std::vector<TopoDS_Edge>   made;
  };

  //! The face across `edge` from `face`, else null.
  TopoDS_Face faceAcross(const TopoDS_Shape& edge, const TopoDS_Shape& face) const;
  //! The faces around `at` from `first`, the face beyond `firstEdge`, to the face beyond
  //! `lastEdge`, and the edges between them.
  void facesAround(const Rail&               rail,
                   const TopoDS_Vertex&      at,
                   const TopoDS_Face&        first,
                   const TopoDS_Edge&        firstEdge,
                   const TopoDS_Edge&        lastEdge,
                   std::vector<TopoDS_Face>& faces,
                   std::vector<TopoDS_Edge>& edges);
  //! The steps of stops: out along one support from the vertex to the faces that stop its contact.
  void walkSupport(StopEnd& st, int side);
  //! Whether a contact ends on `own` where it is carried: it crosses none of its edges (but the
  //! rail's and the one it came onto the face over) and the ball rolls on over no other piece.
  bool endsOnFace(const Piece&       piece,
                  int                side,
                  int                end,
                  double             middle,
                  const StopWalk&    w,
                  const TopoDS_Face& own) const;
  //! All the stop faces in order from side 0's crossing to side 1's.
  void stopFaces(StopEnd& st);
  //! One contact ends on its face: the face the other crosses onto stops the blend and grows.
  //! Where the open contact o meets the stop face's surface, in its own face's parameters.
  void openCorner(StopEnd& st, int o, const Handle(Geom_Surface)& stopSurface);
  void openStop(StopEnd& st);
  //! The two corners, and where the sheet crosses each edge between stop faces.
  void stopPlaces(StopEnd& st);
  //! The end along the sheet on each stop face, from one place to the next.
  void stopEdges(StopEnd& st);
  //! The cuts, what the contacts swept over, and the faces' new edges.
  void stopCuts(StopEnd& st);

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
  //! Where a side's contact, short of the faces that stop its end, crosses onto other faces of
  //! its support's own surface, in the order it crosses them: the place, the fit's u there, the
  //! face it crosses onto and the whole periods its contact is moved by to lie in that face's
  //! parameters. [piece, side, end]
  struct Passed
  {
    TopoDS_Vertex at;
    double        u;
    TopoDS_Face   face;
    gp_Vec2d      shift;
  };
  std::map<std::array<int, 3>, std::vector<Passed>> passages;
  //! The face each part of such a contact lies on, where it is not the support's own.
  NCollection_DataMap<TopoDS_Shape, TopoDS_Face, TopTools_ShapeMapHasher> partFace;

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

} // namespace BRepExt::blending::detail
