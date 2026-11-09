// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Blend (BRepBlendExt_Internal.hxx): a request's blends stage by stage (struct Blend): its
// pieces, how their ends meet, their sheets and splits, and the blends attached.

#include <GeomAlgoExt_Cancel.hxx>
#include <BRepBlendExt_Stages.hxx>
#include <BRepBlendExt_Internal.hxx>
#include <BRepAlgoExt_Cells.hxx>
#include <BRepAlgoExt_FaceBlend.hxx>
#include <BRepAlgoExt_Prism.hxx>
#include <BRepAlgoExt_FillHoles.hxx>
#include <BRepAlgoExt_RoundFace.hxx>
#include <BRepAlgoAPI_Defeaturing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <BRepAlgoAPI_Splitter.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <ElCLib.hxx>
#include <Extrema_ExtPC.hxx>

namespace BRepExt::blending
{

using namespace detail;

//! How far an edge reaches from the origin: the largest coordinate of its box (Rail::magnitude).
static double magnitudeOf(const TopoDS_Edge& edge)
{
  Bnd_Box box;
  BRepBndLib::Add(edge, box);
  if (box.IsVoid())
    return 0;
  double x0, y0, z0, x1, y1, z1;
  box.Get(x0, y0, z0, x1, y1, z1);
  return std::max(
    {std::abs(x0), std::abs(y0), std::abs(z0), std::abs(x1), std::abs(y1), std::abs(z1)});
}

//! The side as a cliff on `edge`: the ball rests on the edge, its angle about it measured from
//! the axis most across the edge.
static bool restOn(Side& side, const TopoDS_Edge& edge)
{
  double     f, l;
  const auto curve = BRep_Tool::Curve(edge, f, l);
  if (curve.IsNull())
    return false;
  gp_Pnt p;
  gp_Vec d;
  curve->D1((f + l) / 2, p, d);
  const double x = std::abs(d.X()), y = std::abs(d.Y()), h = std::abs(d.Z());
  side.cliff          = edge;
  side.cliffCurve     = curve;
  side.cliffReference = x <= y && x <= h ? gp_Vec(1, 0, 0)
                        : y <= h         ? gp_Vec(0, 1, 0)
                                         : gp_Vec(0, 0, 1);
  return true;
}

//! The edge of `face` first met going across it from `edge` at its parameter `t`, along the face's
//! parameters at right angles to the edge's curve on it, into the face: empty where none is met.
static TopoDS_Edge firstAcross(const TopoDS_Face& face,
                               const TopoDS_Edge& edge,
                               double             t,
                               Inside&            inside)
{
  double     f, l;
  const auto onFace = BRep_Tool::CurveOnSurface(edge, face, f, l);
  if (onFace.IsNull())
    return {};
  gp_Pnt2d at;
  gp_Vec2d d;
  onFace->D1(t, at, d);
  if (!(d.Magnitude() > gp::Resolution()))
    return {};
  gp_Dir2d normal(-d.Y(), d.X());
  double   u0, u1, v0, v1;
  BRepTools::UVBounds(face, u0, u1, v0, v1);
  const double step = 1e-6 * std::max(u1 - u0, v1 - v0);
  if (inside.out(at.Translated(gp_Vec2d(normal) * step)))
    normal.Reverse();
  const Handle(Geom2d_Line) ray = new Geom2d_Line(at, normal);
  TopoDS_Edge               first;
  double                    nearest = Precision::Infinite();
  for (TopExp_Explorer e(face, TopAbs_EDGE); e.More(); e.Next())
  {
    const TopoDS_Edge other = TopoDS::Edge(e.Current());
    if (other.IsSame(edge))
      continue;
    double     a, b;
    const auto curve = BRep_Tool::CurveOnSurface(other, face, a, b);
    if (curve.IsNull())
      continue;
    Geom2dAPI_InterCurveCurve meet(
      new Geom2d_TrimmedCurve(ray, 0, 4 * std::hypot(u1 - u0, v1 - v0)),
      new Geom2d_TrimmedCurve(curve, a, b),
      tol::parameter);
    for (int i = 1; i <= meet.NbPoints(); ++i)
    {
      const double away = at.Distance(meet.Point(i));
      if (away > step && away < nearest)
      {
        nearest = away;
        first   = other;
      }
    }
  }
  return first;
}

//! The far edge of a face the blend's contact runs off along the whole edge: the edge first met
//! across it from the blended one, the same along it; whether the face is that strip alone (four
//! edges, the far one joined to each end of the blended one by one), and whether the far edge ends
//! with the strip at one end of the blend at least (the face's edge from the blended edge's vertex
//! there reaches it). `inside` classifies points of the face (Blend::inside keeps one per face).
struct Far
{
  TopoDS_Edge edge;
  bool        whole = false, ends = false;
};

static std::optional<Far> farEdge(const Rail& rail, int side, Inside& inside)
{
  const TopoDS_Face& face = rail.s[side].face;
  Far                out;
  double             f, l;
  BRep_Tool::Range(rail.edge, f, l);
  for (const double share : {0.1, 0.5, 0.9})
  {
    const TopoDS_Edge across = firstAcross(face, rail.edge, f + (l - f) * share, inside);
    if (across.IsNull() || (!out.edge.IsNull() && !out.edge.IsSame(across)))
      return {};
    out.edge = across;
  }
  if (BRep_Tool::Degenerated(out.edge))
    return {};
  std::vector<TopoDS_Edge> edges;
  int                      wires = 0;
  for (TopoDS_Iterator w(face); w.More(); w.Next(), ++wires)
    for (TopoDS_Iterator e(w.Value()); e.More(); e.Next())
      edges.push_back(TopoDS::Edge(e.Value()));
  TopoDS_Vertex a, z, f0, f1;
  TopExp::Vertices(rail.edge, a, z);
  TopExp::Vertices(out.edge, f0, f1);
  out.whole = wires == 1 && edges.size() == 4;
  for (const auto& e : edges)
  {
    if (e.IsSame(rail.edge) || e.IsSame(out.edge))
      continue;
    TopoDS_Vertex p, q;
    TopExp::Vertices(e, p, q);
    const bool atRail = p.IsSame(a) || p.IsSame(z) || q.IsSame(a) || q.IsSame(z);
    const bool atFar  = p.IsSame(f0) || p.IsSame(f1) || q.IsSame(f0) || q.IsSame(f1);
    out.whole         = out.whole && atRail;
    out.ends          = out.ends || (atRail && atFar);
  }
  return out;
}

//! Where a blend's contact runs off its face `side` along the whole edge, the face is narrower
//! than the blend: a strip between the edge and its far edge. The blend replaces the strip; where
//! the strip goes on past both ends, it keeps the round and trims it by the face beyond (a notch).
//! The Round only, which is the ball's own section: no law is tangent to an edge.
static bool cliffOf(Rail& rail, int side, const Far& far)
{
  if (rail.law != Law::Conic || std::abs(rail.rho - 0.5) > limit::equal
      || std::abs(rail.s[0].distance - rail.s[1].distance) > limit::equal * rail.scale || !far.ends)
    return false;
  rail.s[side].cliffWhole = far.whole;
  return restOn(rail.s[side], far.edge);
}

//! A request's blends, stage by stage: its pieces, how their ends meet, their sheets and splits,
//! and the blends attached.
struct Blend
{
  Blend(const Request&            request_,
        BRepBuilderAPI_Copy&      copy_,
        std::vector<TopoDS_Shape> byCentre_ = {})
      : request(request_),
        copy(copy_),
        op(request_.operation),
        body(copy_.Shape()),
        byCentre(std::move(byCentre_))
  {
    TopExp::MapShapesAndUniqueAncestors(body, TopAbs_EDGE, TopAbs_FACE, facesOfEdge);
    TopExp::MapShapesAndUniqueAncestors(body, TopAbs_VERTEX, TopAbs_EDGE, edgesOfVertex);
  }

  //! The pieces: one per edge asked and per tangent edge added, with its supports, distances and
  //! kind.
  void rails();
  //! A full round's one piece.
  void fullRail(const std::vector<TopoDS_Shape>& asked);
  //! How each piece's ends meet the others': a joint, a mitre or a corner.
  void joints();
  //! The corner where more than three blended edges meet.
  void vertexBlend(const TopoDS_Vertex& vertex, const std::vector<std::pair<int, int>>& meeting);
  //! Each rail carried on past an end that a face stops.
  void carries();
  //! Each piece's stations and sheet, split where a contact crosses a smooth edge onto the next
  //! face.
  void sheets();
  bool continueWholeSupport(int k, int side, Stations& stations);

  //! The loops found on each piece, and each piece's creases where its round crosses itself.
  using Folds   = std::vector<std::pair<int, std::vector<Fold>>>;
  using Creases = std::vector<std::vector<std::shared_ptr<const Crease>>>;
  //! What a step of a piece's sheet leads to: on, the same attempt again, or sheets() ends.
  enum class Next
  {
    Go,
    Retry,
    Stop
  };

  //! Where a contact first changes between inside its face and outside: its side, the fit's u
  //! there, the edge it crosses (and that edge's parameter), whether it enters the face, whether
  //! the face beyond trims it (a notch), and over a seam the period it is carried by.
  struct Crossing
  {
    int         side = -1;
    double      at = 0, overAt = 0;
    bool        entering = false, notched = false;
    TopoDS_Edge over;
    bool        overSeam = false;
    gp_Vec2d    seamShift;
  };

  //! A piece's sheet as fitted: exact where it can be, the exact contacts, and whether the fit is
  //! the circular pilot that finds the support crossings of a G2 section first.
  struct SheetFit
  {
    std::optional<Sheet>  exact;
    std::optional<Fitted> contacts;
    bool                  contactPilot = false;
  };

  //! The steps of sheets(): the sections' spans shared along each chain.
  void shareHalves();
  //! Piece k's sheet, split where a contact crosses an edge; false where sheets() ends.
  bool sheetOf(int k, Folds& folds, Creases& creased);
  //! The rail's windows, chord, meeting radius and reach made ready for the fit.
  void prepareRail(Piece& piece, Stations& stations);
  //! The sheet fitted (or exact); false where the chain is found through its centres instead.
  bool fitSheet(int                                         k,
                int                                         attempt,
                Stations&                                   stations,
                SheetFit&                                   out,
                std::vector<std::shared_ptr<const Crease>>& creases);
  //! A contact wholly off its face: a cliff or a notch tried again, or a chain's end absorbed.
  Next offFace(int k, Stations& stations, double u0, double u1);
  //! The first edge of its face a contact crosses between the piece's ends.
  void crossingOverEdge(int k, double u0, double u1, Crossing& c);
  //! Whether the faces either side of `over` are tangent where the contact crosses it.
  bool smoothOver(const Piece& piece, int s2, const TopoDS_Edge& over, const Meeting& found) const;
  //! Refused where the contact runs over `over` into a hole or a boss that stays.
  void refuseFeature(const Piece& piece, int s2, const TopoDS_Edge& over) const;
  //! Whether what stops the blend at an end trims the contact where it crosses `over`.
  bool stopsNearEnd(const Piece& piece, const TopoDS_Edge& over, const Meeting& found) const;
  //! The first seam of its face a contact crosses, where it crosses no edge.
  void crossingOverSeam(int k, double u0, double u1, Crossing& c);
  //! The piece split at the crossing, the part off the face rolling on the next face.
  void splitAt(int k, Stations& stations, const Crossing& c);
  //! The folds of piece k's contacts on a cap; whether a contact that runs back is left over.
  bool capFoldsOf(int k, Dips& dips, Folds& folds);
  //! The folds settled; false where the chain is found through its centres instead.
  bool settleFolds(int k, Folds& folds, bool onCliff);
  //! Whether through the centres piece k is fitted and no contact runs back.
  bool centreMends(int k);
  //! Whether the edge a contact on face s2 leaves over at `at` is convex there.
  bool convexOver(const Rail& rail, int s2, const TopoDS_Edge& edge, const Meeting& at) const;
  //! Each round's creases added to its folds, in order along it.
  void creaseFolds(Folds& folds, Creases& creased);
  //! Each edge cut at its folds' ends, a piece on the one sheet between each loop and the next.
  void cutAtFolds(const Folds& folds);
  //! The loops where a contact on a cap runs back across a joint into the next piece.
  void crossFolds(const std::vector<std::pair<int, std::vector<Fold>>>& folds);
  //! One kind of sheet along each chain: exact where every piece of it can be, else fitted.
  void chains();
  //! The parts of one split edge on one sheet.
  void joinSplits();
  //! Which sides pair at each joint.
  void pairJoints();
  //! The edges of each chain whose blends end in different balls at a joint, where the balls
  //! found through the edge's points fold back (Rail::byCentre); empty where every joint agrees.
  std::vector<TopoDS_Shape> foldedJoints() const;
  //! The blends attached, with the history from the input.
  Blended result();
  void    splitMeeting(Blended& out);
  //! Which points of a face's surface lie off it: made once for each face, as a face of many
  //! holes takes long to make it for.
  Inside& inside(const TopoDS_Face& face);
  //! farEdge of a rail's side, worked out once for each edge and face: the sheets are made again
  //! for each attempt, and the far edge depends on the edge and the face alone.
  const std::optional<Far>& farOf(const Rail& rail, int side);
  //! Whether piece k's rail is one whose balls may be found through their centres instead
  //! (Rail::byCentre): a round's or a G2 blend's, off faces (not a cliff), not yet so.
  bool mayCentre(int k) const;
  //! The balls of piece k's chain (its tangent edges and split parts, which share sections at
  //! their joints) found through their centres from now on: sheets() is made again.
  void centreChain(int k);

  const Request&       request;
  BRepBuilderAPI_Copy& copy;
  const char*          op;
  const TopoDS_Shape   body;
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
                            facesOfEdge, edgesOfVertex;
  std::vector<TopoDS_Shape> edges; // asked and added
  TopoDS_Shape              reference;
  std::vector<Piece>        pieces;
  //! Where a contact on a flat cap runs back but crosses itself nowhere within its piece: the
  //! places, for the loop across a joint into the next piece (crossFolds).
  std::vector<std::pair<int, std::vector<Fold>>> open;

  //! The refusal of a fold the blend finds no loop for, on piece k: its contact's fold, as
  //! checkFold says.
  [[noreturn]] void refuseFold(int k, const char* why)
  {
    Dips dips{pieces[k].rail, pieces[k].fitted};
    checkFold(dips);
    fail(op, unsupported_geometry, why, {pieces[k].rail.edge}, FailureKind::Curvature);
  }

  std::vector<Corner>                 corners;
  std::vector<std::optional<Station>> seeds; // [piece]: where a split part's stations start
  std::vector<TopoDS_Shape>
    byCentre; // edges whose balls are found through their centres (foldedJoints, centreChain)
  std::vector<TopoDS_Edge>
       ending; // a chain's last edges, wholly under the blend of the edge before (absorbEnds)
  void absorbEnds();
  NCollection_DataMap<TopoDS_Shape, int, TopTools_ShapeMapHasher> insideOf; // face -> insides
  std::deque<Inside>                                              insides;
  std::deque<std::tuple<TopoDS_Shape, TopoDS_Shape, std::optional<Far>>> fars; // edge, face
  bool again = false; // sheets() to be made again (centreChain)

  bool centred(const TopoDS_Shape& edge) const
  {
    return std::any_of(byCentre.begin(), byCentre.end(), [&](const TopoDS_Shape& e) {
      return e.IsSame(edge);
    });
  }
};

bool Blend::mayCentre(int k) const
{
  const Rail& rail = pieces[k].rail;
  return !rail.byCentre && (rail.law == Law::Conic || rail.law == Law::G2) && !rail.s[0].cliffCurve
         && !rail.s[1].cliffCurve;
}

//! A chain's last edge that a contact is off the whole way, past the edge of its face where the
//! face that stops the chain begins (a short edge, a sphere's seam or a patch's edge just before
//! the stop): no blend of its own, the one before runs on over it to that stop. It is made again
//! so, the ball being the same whichever edge it is found along.
void Blend::absorbEnds()
{
  for (const auto& edge : ending)
  {
    const auto found = std::find_if(pieces.begin(), pieces.end(), [&](const Piece& p) {
      return p.rail.edge.IsSame(edge);
    });
    if (found == pieces.end())
      continue;
    const int     b       = int(found - pieces.begin());
    const int     stopEnd = found->joint[0].meets == Meets::Stop ? 0 : 1;
    const auto    joint   = found->joint[1 - stopEnd];
    TopoDS_Vertex ends[2];
    TopExp::Vertices(TopoDS::Edge(edge.Oriented(TopAbs_FORWARD)), ends[0], ends[1]);
    Piece& a                = pieces[joint.piece];
    a.joint[joint.end]      = {};
    a.stopVertex[joint.end] = ends[stopEnd];
    a.swallowed.push_back(TopoDS::Edge(edge));
    a.swallowed.insert(a.swallowed.end(), found->swallowed.begin(), found->swallowed.end());
    // On past its vertex by the length of the edge it runs on over, and on from there as to a stop.
    gp_Pnt p;
    gp_Vec d;
    Rail&  rail = a.rail;
    rail.curve.D1(joint.end == 0 ? rail.first : rail.last, p, d);
    const double over =
      GCPnts_AbscissaPoint::Length(BRepAdaptor_Curve(TopoDS::Edge(edge))) / d.Magnitude();
    const double more = 4 * rail.scale / d.Magnitude();
    if (joint.end == 0)
    {
      rail.first -= over;
      rail.from = rail.first - more;
    }
    else
    {
      rail.last += over;
      rail.to = rail.last + more;
    }
    pieces.erase(found);
    for (auto& piece : pieces)
      for (auto& j : piece.joint)
        if (j.meets != Meets::Stop && j.piece > b)
          --j.piece;
    for (auto& corner : corners)
      for (auto& [p2, e2] : corner.ends)
        if (p2 > b)
          --p2;
  }
}

void Blend::centreChain(int k)
{
  std::vector<int>  todo = {k};
  std::vector<bool> seen(pieces.size(), false);
  seen[k] = true;
  while (!todo.empty())
  {
    const int p = todo.back();
    todo.pop_back();
    if (!centred(pieces[p].rail.edge))
      byCentre.push_back(pieces[p].rail.edge);
    for (const auto& joint : pieces[p].joint)
      if ((joint.meets == Meets::Joint || joint.meets == Meets::Split) && !seen[joint.piece])
      {
        seen[joint.piece] = true;
        todo.push_back(joint.piece);
      }
  }
  again = true;
}

Inside& Blend::inside(const TopoDS_Face& face)
{
  if (const int* at = insideOf.Seek(face))
    return insides[*at];
  insideOf.Bind(face, int(insides.size()));
  return insides.emplace_back(face);
}

const std::optional<Far>& Blend::farOf(const Rail& rail, int side)
{
  const TopoDS_Face& face = rail.s[side].face;
  for (const auto& [edge, of, far] : fars)
    if (edge.IsSame(rail.edge) && of.IsSame(face))
      return far;
  return std::get<2>(fars.emplace_back(rail.edge, face, farEdge(rail, side, inside(face))));
}

//! Refused where the fitted blend crosses itself between its contacts: the envelope of the balls
//! along a path that bends more tightly than the ball (a pipe surface c + r (cos a N + sin a B)
//! moves along the path at |c'| (1 - kappa r cos a), which changes sign there), so its sections
//! run back against the way the balls roll. Read off the sheet: S_u against the station's rolling
//! direction, at 65 x 7 places; where it turns against it, the surface has folded. A round's sheet
//! is the balls' envelope, and the blend the boundary of the balls' union, which is the sheet
//! trimmed where it crosses itself (Farouki and Neff, *Analytic properties of plane offset
//! curves*, CAGD 7 (1990), on the trimmed offset): its creases are returned (crossings). Other
//! laws' sections are not the balls', and are refused there.
static std::vector<std::shared_ptr<const Crease>> checkCrossing(const Rail& rail,
                                                                Stations&   stations,
                                                                const Fitted& fitted)
{
  if (fitted.blend.IsNull())
    return {};
  // Where a contact runs back, its loop is patched or refused as a fold of the contact
  // (capFolds, checkFold): this is the crossing between contacts that both roll on.
  Dips dips{rail, fitted};
  if (runsBack(dips))
    return {};
  const bool round = rail.law == Law::Conic && std::abs(rail.rho - 0.5) < limit::equal
                     && std::abs(rail.s[0].distance - rail.s[1].distance)
                          < limit::equal * rail.scale;
  constexpr int along = 64, across = 8;
  for (int i = 0; i <= along; ++i)
  {
    const double u = double(i) / along;
    Station      s;
    if (!stations.at(rail.t(u), s))
      continue;
    for (int j = 1; j < across; ++j)
    {
      gp_Pnt p;
      gp_Vec du, dv;
      fitted.blend->D1(u, double(j) / across, p, du, dv);
      if (du.Dot(s.along) < -limit::agree * du.Magnitude() - 1e-3 * du.Magnitude())
      {
        // A crease that cannot be traced within the piece is refused as before.
        if (round)
          try
          {
            return crossings(rail, stations, fitted);
          }
          catch (const Standard_Failure&)
          {
          }
        fail(rail.op,
             unsupported_geometry,
             "the blend crosses itself between its contacts, where the balls' path bends more "
             "tightly than the ball, which the blend does not trim yet",
             {rail.edge},
             FailureKind::Curvature);
      }
    }
  }
  return {};
}

//! Whether a conic section of the rail turns by more than 120 degrees somewhere along its edge
//! (Rail::halves): the angle between the faces' normals at the edge's points, which the section
//! turns by. Sampled at 17 places.
static bool turnsWide(const Rail& rail)
{
  if (rail.law != Law::Conic && rail.law != Law::Chordal)
    return false;
  constexpr int samples = 16;
  for (int k = 0; k <= samples; ++k)
  {
    const double t = rail.first + (rail.last - rail.first) * k / samples;
    gp_Vec       n[2];
    bool         found = true;
    for (int i = 0; i < 2 && found; ++i)
    {
      double     f, l;
      const auto pcurve = BRep_Tool::CurveOnSurface(rail.edge, rail.s[i].face, f, l);
      Jet        j;
      found = !pcurve.IsNull() && !rail.s[i].cliffCurve
              && jetAt(rail.s[i], pcurve->Value(t).X(), pcurve->Value(t).Y(), j);
      if (found)
        n[i] = j.n;
    }
    if (found && n[0].Angle(n[1]) > 2 * M_PI / 3)
      return true;
  }
  return false;
}

//! The rail's size law (Rail::sizes) from the (place, size) pairs asked for its edge: the factor on
//! its distances by the edge's parameter, positive all along; round an edge that closes, the size
//! must end as it begins. Its scale is the largest.
static void sizesOf(Rail& rail, const std::vector<std::pair<double, double>>& points, double first)
{
  const char* op = rail.op;
  const int   n  = int(points.size());
  if (n < 2 || points.front().first != 0 || points.back().first != 1)
    fail(op,
         invalid_input,
         "a size law runs from the edge's start (0) to its end (1)",
         {rail.edge});
  for (int i = 0; i < n; ++i)
  {
    const auto [place, size] = points[i];
    if (!std::isfinite(place) || !std::isfinite(size) || !(size > 0)
        || (i > 0 && !(place > points[i - 1].first)))
      fail(op, invalid_input, "a size law takes rising places and positive sizes", {rail.edge});
    rail.sizes.places.push_back(rail.first + (rail.last - rail.first) * place);
    rail.sizes.values.push_back(size / first);
  }
  TopoDS_Vertex a, b;
  TopExp::Vertices(rail.edge, a, b);
  rail.sizes.periodic = !a.IsNull() && a.IsSame(b);
  if (rail.sizes.periodic
      && std::abs(points.front().second - points.back().second) > limit::equal * first)
    fail(op,
         invalid_input,
         "the size round a closed edge ends where it begins: "
           + std::to_string(points.front().second) + " and " + std::to_string(points.back().second)
           + " are at one point",
         {rail.edge},
         FailureKind::RadiusEndpoints);
  rail.sizes.solve();
  const auto bounds = rail.sizes.bounds();
  if (!GeomAlgoExt::interval::valid(bounds) || !(bounds.low > 0))
    fail(op,
         unsupported_geometry,
         "the size law falls to zero between its places or its sign is unresolved",
         {rail.edge},
         FailureKind::RadiusLaw);
  rail.scale *= bounds.high;
}

void Blend::rails()
{

  std::vector<TopoDS_Shape> asked, mayAdd;
  for (const auto& edge : request.edges)
    asked.push_back(copy.ModifiedShape(edge));
  for (const auto& edge : request.mayAdd)
    mayAdd.push_back(copy.ModifiedShape(edge));
  if (request.law == Law::Full)
    return fullRail(asked);
  edges     = withTangent(body, asked, std::max(request.first, request.second), mayAdd);
  reference = request.reference.IsNull() ? TopoDS_Shape() : copy.ModifiedShape(request.reference);
  pieces.reserve(edges.size() * 8);
  for (const auto& shape : edges)
  {
    const TopoDS_Edge edge = TopoDS::Edge(shape);
    if (std::any_of(pieces.begin(), pieces.end(), [&](const Piece& p) {
          return p.rail.edge.IsSame(edge);
        }))
      continue;
    // A seam has its face on both sides and an edge between tangent faces has nothing to

    if (facesOfEdge.Contains(edge) && facesOfEdge.FindFromKey(edge).Extent() == 1
        && BRep_Tool::IsClosed(edge, TopoDS::Face(facesOfEdge.FindFromKey(edge).First())))
      continue;
    if (!facesOfEdge.Contains(edge) || facesOfEdge.FindFromKey(edge).Extent() != 2)
      fail(op, invalid_input, "a blended edge must lie between two faces", {edge});
    pieces.emplace_back();
    Rail& rail = pieces.back().rail;
    rail.op    = op;
    rail.boundedFit = request.fitWork == FitWork::Bounded;
    rail.reads      = request.reads;
    rail.edge  = edge;
    rail.curve.Initialize(edge);
    rail.law = request.law;
    rail.rho = request.law == Law::Chordal ? 0.5 : request.rho;
    int i    = 0;
    for (const auto& face : facesOfEdge.FindFromKey(edge))
    {
      rail.s[i].face    = TopoDS::Face(face);
      rail.s[i].surface = BRep_Tool::Surface(rail.s[i].face);
      rail.s[i].outward = rail.s[i].face.Orientation() == TopAbs_REVERSED ? -1 : 1;
      ++i;
    }
    const double second        = request.second > 0 ? request.second : request.first;
    const bool   firstOnSecond = !reference.IsNull() && rail.s[1].face.IsSame(reference);
    rail.s[0].distance         = firstOnSecond ? second : request.first;
    rail.s[1].distance         = firstOnSecond ? request.first : second;
    rail.angle                 = request.angle;
    rail.reference             = firstOnSecond ? 1 : 0;
    rail.scale                 = std::max(rail.s[0].distance, rail.s[1].distance);
    rail.magnitude             = magnitudeOf(edge);
    rail.first                 = rail.curve.FirstParameter();
    rail.last                  = rail.curve.LastParameter();
    for (size_t k = 0; k < asked.size() && k < request.radii.size(); ++k)
      if (asked[k].IsSame(edge) && !request.radii[k].empty())
        sizesOf(rail, request.radii[k], request.first);
    const int kind =
      convexity(op, edge, rail.s[0].face, rail.s[1].face, rail.scale, &rail.touching);
    if (kind == 0)
    {
      pieces.pop_back();
      continue;
    }
    rail.sigma  = kind;
    rail.halves = turnsWide(rail);
    TopoDS_Vertex start, end;
    TopExp::Vertices(edge, start, end);
    rail.closed   = !start.IsNull() && start.IsSame(end);
    rail.byCentre = std::any_of(byCentre.begin(), byCentre.end(), [&](const TopoDS_Shape& e) {
      return e.IsSame(edge);
    });
  }
}

void Blend::fullRail(const std::vector<TopoDS_Shape>& asked)
{
  if (asked.size() != 2)
    fail(op, invalid_input, "a full round takes the two edges of one face");
  TopoDS_Edge rails[2];
  TopoDS_Face beside[2], between;
  for (int i = 0; i < 2; ++i)
  {
    rails[i] = TopoDS::Edge(asked[i]);
    if (!facesOfEdge.Contains(rails[i]) || facesOfEdge.FindFromKey(rails[i]).Extent() != 2)
      fail(op, invalid_input, "a blended edge must lie between two faces", {rails[i]});
  }
  for (const auto& face : facesOfEdge.FindFromKey(rails[0]))
    for (const auto& other : facesOfEdge.FindFromKey(rails[1]))
      if (face.IsSame(other))
        between = TopoDS::Face(face);
  if (between.IsNull())
    fail(op, invalid_input, "a full round's two edges must bound one face", {rails[0], rails[1]});
  for (int i = 0; i < 2; ++i)
    for (const auto& face : facesOfEdge.FindFromKey(rails[i]))
      if (!face.IsSame(between))
        beside[i] = TopoDS::Face(face);
  const auto sideOf = [](const TopoDS_Face& face) {
    Side side;
    side.face    = face;
    side.surface = BRep_Tool::Surface(face);
    side.outward = face.Orientation() == TopAbs_REVERSED ? -1 : 1;
    return side;
  };
  edges = {rails[0]};
  pieces.reserve(8);
  pieces.emplace_back();
  Rail& rail = pieces.back().rail;
  rail.op    = op;
  rail.boundedFit = request.fitWork == FitWork::Bounded;
  rail.reads      = request.reads;
  rail.edge  = rails[0];
  rail.curve.Initialize(rails[0]);
  rail.law       = Law::Full;
  rail.magnitude = magnitudeOf(rails[0]);
  rail.s[0]      = sideOf(beside[0]);
  rail.s[1] = rail.middle = sideOf(between);
  rail.first              = rail.curve.FirstParameter();
  rail.last               = rail.curve.LastParameter();
  TopoDS_Vertex start, end;
  TopExp::Vertices(rails[0], start, end);
  rail.closed    = !start.IsNull() && start.IsSame(end);
  Inside     onStrip(rail.s[1].face);
  const auto far = farEdge(rail, 1, onStrip);
  if (!far || !far->edge.IsSame(rails[1]) || !far->whole)
    fail(op,
         unsupported_geometry,
         "a full round's face must be only the strip between its two edges",
         {rails[0], rails[1]});
  // Half the width across, which the ball's radius is near: for tolerances and the carry.
  BRepExtrema_DistShapeShape width(rails[0], rails[1]);
  if (!width.IsDone() || !(width.Value() > tol::length))
    fail(op, no_solution, "a full round's two edges meet", {rails[0], rails[1]});
  rail.scale     = width.Value() / 2;
  const int kind = convexity(op, rails[0], beside[0], between, rail.scale);
  if (kind == 0 || convexity(op, rails[1], between, beside[1], rail.scale) != kind)
    fail(op,
         unsupported_geometry,
         "a full round's two edges must both be convex or both concave",
         {rails[0], rails[1]});
  rail.sigma = kind;
  // The notch of the face between, whose face beyond the ball rests on (Rail::middle).
  rail.s[1]          = sideOf(beside[1]);
  rail.s[1].face     = between;
  rail.s[1].trim     = beside[1];
  rail.s[1].trimEdge = rails[1];
  rail.s[0].distance = rail.s[1].distance = rail.scale;
}

void Blend::joints()
{
  // Where two blended edges meet tangent to each other, one blend runs on from one into the
  // other; where more meet, or at an angle, the corner is blended (milestone 3).
  NCollection_IndexedDataMap<TopoDS_Shape,
                             std::vector<std::pair<int, int>>,
                             TopTools_ShapeMapHasher>
                                            endsAt;
  std::vector<std::array<TopoDS_Vertex, 2>> endOf(pieces.size());
  for (int k = 0; k < int(pieces.size()); ++k)
  {
    auto& ends = endOf[k];
    TopExp::Vertices(TopoDS::Edge(pieces[k].rail.edge.Oriented(TopAbs_FORWARD)), ends[0], ends[1]);
    for (int end = 0; end < 2; ++end)
    {
      if (!endsAt.Contains(ends[end]))
        endsAt.Add(ends[end], {});
      endsAt.ChangeFromKey(ends[end]).push_back({k, end});
    }
  }
  for (int v = 1; v <= endsAt.Extent(); ++v)
  {
    const auto& meeting = endsAt(v);
    if (meeting.size() == 1)
      continue;
    const auto [a, ea]          = meeting[0];
    const TopoDS_Vertex& vertex = endOf[a][ea];
    if (meeting.size() == 3)
    {
      Corner corner{vertex, meeting};
      for (const auto& [k, end] : meeting)
        pieces[k].joint[end] = {Meets::Corner, int(corners.size()), 0, false};
      corners.push_back(corner);
      continue;
    }
    if (meeting.size() > 3)
    {
      vertexBlend(vertex, meeting);
      continue;
    }
    const auto [b, eb] = meeting.back();
    gp_Vec da, db;
    if (!leaving(pieces[a].rail.edge, vertex, da, ea)
        || !leaving(pieces[b].rail.edge, vertex, db, eb))
      fail(op,
           unsupported_geometry,
           "a blended edge has no direction where it meets another",
           {pieces[a].rail.edge});
    const Meets meets = tangent(da, db) ? Meets::Joint : Meets::Miter;
    if (meets == Meets::Miter && a == b)
      fail(op,
           unsupported_geometry,
           "an edge that closes on itself at an angle is blended with its corner, which the blend "
           "does not do yet",
           {pieces[a].rail.edge});
    pieces[a].joint[ea] = {meets, b, eb, false};
    pieces[b].joint[eb] = {meets, a, ea, false};
  }
}

// Where more than three blended edges meet: the rolling ball's corner, where the body is what the
// balls inside it fill (its opening by the ball). Where one ball rests on all the faces round the
// vertex, the corner is that ball's cap about its apex. Else, between planes, the corner parts as
// the region of the balls' centres does (the faces' planes offset by the radius): each corner of
// that region is a ball resting on three faces or more, faces next to each other round the vertex,
// and two such balls that rest on the same two faces not next to each other are joined by the round
// between those faces, which meet only at the vertex: along the line where their planes meet, from
// one ball to the other. Four faces part so into two balls and a round between the faces across
// from each other, as

void Blend::vertexBlend(const TopoDS_Vertex&                    vertex,
                        const std::vector<std::pair<int, int>>& meeting)
{
  const TopoDS_Edge lead = pieces[meeting[0].first].rail.edge;
  const double      r    = pieces[meeting[0].first].rail.s[0].distance,
               sigma     = pieces[meeting[0].first].rail.sigma;
  // Balls whose centres are that near are one (limit::one); a face is met by a ball where its
  // offset is that near the centre.
  const double scale = pieces[meeting[0].first].rail.scale,
               near  = limit::one(scale, BRep_Tool::Pnt(vertex).XYZ().Modulus());
  const auto refuse  = [&](const char* why) { fail(op, unsupported_geometry, why, {lead}); };
  // Rounds of one radius make the rolling ball's corner (below); conics of one shape, all convex
  // or all concave, one corner bounded by their ends in turn (Attach::splitCorner).
  const Rail& leadRail = pieces[meeting[0].first].rail;
  bool        rounds = true, conics = true;
  for (const auto& [k, end] : meeting)
  {
    const Rail& rail = pieces[k].rail;
    rounds &= rail.law == Law::Conic && std::abs(rail.rho - 0.5) < limit::equal
              && std::abs(rail.s[0].distance - r) < limit::equal * scale
              && std::abs(rail.s[1].distance - r) < limit::equal * scale && rail.sigma == sigma;
    conics &= rail.law == Law::Conic && std::abs(rail.rho - leadRail.rho) < limit::equal
              && rail.sigma == sigma;
  }
  if (!rounds && !conics)
    refuse("where more than three blended edges meet, only the corner of rounds of one "
           "radius, all convex or all concave, is blended");
  // The faces round the vertex in turn: after[i] is the end of the blend between ring[i] and ring[i
  // + 1].
  const int                        count = int(meeting.size());
  std::vector<TopoDS_Face>         ring  = {pieces[meeting[0].first].rail.s[0].face};
  std::vector<std::pair<int, int>> after;
  std::vector<bool>                used(count, false);
  for (bool found = true; found && int(after.size()) < count;)
  {
    found = false;
    for (int m = 0; m < count && !found; ++m)
      for (int side = 0; side < 2 && !found && !used[m]; ++side)
        if (pieces[meeting[m].first].rail.s[side].face.IsSame(ring.back()))
        {
          used[m] = found = true;
          after.push_back(meeting[m]);
          ring.push_back(pieces[meeting[m].first].rail.s[1 - side].face);
        }
  }
  const char* ringless = "the faces at a vertex of more than three blended edges do not go round "
                         "it in turn, each between two of the blends";
  if (int(after.size()) != count || !ring.back().IsSame(ring.front()))
    refuse(ringless);
  ring.pop_back();
  for (int i = 0; i < count; ++i)
    for (int j = i + 1; j < count; ++j)
      if (ring[i].IsSame(ring[j]))
        refuse(ringless);
  if (!rounds)
  {
    for (const auto& [k, end] : after)
      pieces[k].joint[end] = {Meets::Corner, int(corners.size()), 0, false};
    corners.push_back({vertex, after});
    return;
  }
  std::vector<Side> sides(count);
  for (int i = 0; i < count; ++i)
    sides[i] = Side{ring[i],
                    BRep_Tool::Surface(ring[i]),
                    ring[i].Orientation() == TopAbs_REVERSED ? -1.0 : 1.0,
                    r};
  const auto corner = [&](const std::vector<std::pair<int, int>>& ends) {
    for (const auto& [k, end] : ends)
      pieces[k].joint[end] = {Meets::Corner, int(corners.size()), 0, false};
    corners.push_back({vertex, ends});
  };
  const gp_Pnt        V = BRep_Tool::Pnt(vertex);
  std::vector<gp_Vec> n(count);
  for (int i = 0; i < count; ++i)
    if (!normalAt(ring[i], pieces[after[i].first].rail.edge, vertex, n[i]))
      refuse("a face at a vertex of more than three blended edges has no normal there");
  // Faces tangent to each other at the vertex: the ball rests on them all there, and each blend
  // closes to the vertex (Rail::pinch).
  if (std::all_of(n.begin(), n.end(), [&](const gp_Vec& v) {
        return v.Angle(n[0]) <= tol::coarseAngle;
      }))
  {
    for (const auto& [k, end] : after)
      pieces[k].rail.pinch[end] = true;
    corner(after);
    corners.back().pinch = true;
    return;
  }
  if (restingBall(sides, V, r, sigma, scale))
    return corner(after);
  // The corners of the centres' region: where three of the offset planes meet, clear of the
  // others (on the material's side of each for a convex corner, the air's for a concave one).
  for (int i = 0; i < count; ++i)
    if (GeomAdaptor_Surface(sides[i].surface).GetType() != GeomAbs_Plane)
      refuse("where more than three blended edges meet between faces that are not all planes and "
             "no one ball rests on them all, the blend does not blend the corner yet");

  struct Cell
  {
    gp_Pnt           centre;
    std::vector<int> faces; // ascending: in turn round the vertex
  };

  std::vector<Cell> cells;
  for (int a = 0; a < count; ++a)
    for (int b = a + 1; b < count; ++b)
      for (int c = b + 1; c < count; ++c)
      {
        std::array<std::array<double, 4>, 3> m{};
        const int                            rows[3] = {a, b, c};
        for (int i = 0; i < 3; ++i)
        {
          for (int j = 0; j < 3; ++j)
            m[i][j] = n[rows[i]].Coord(j + 1);
          m[i][3] = sigma * r;
        }
        std::array<double, 3> y{};
        if (!solve<3>(m, y))
          continue;
        const gp_Vec     offset(y[0], y[1], y[2]);
        std::vector<int> tight;
        bool             clear = true;
        for (int i = 0; i < count && clear; ++i)
        {
          const double off = sigma * n[i].Dot(offset) - r;
          if (off < -near)
            clear = false;
          else if (off <= near)
            tight.push_back(i);
        }
        if (!clear)
          continue;
        const gp_Pnt centre = V.Translated(offset);
        const auto   same   = std::find_if(cells.begin(), cells.end(), [&](const Cell& cell) {
          return cell.centre.Distance(centre) <= near;
        });
        if (same == cells.end())
          cells.push_back({centre, tight});
        else
        {
          same->faces.insert(same->faces.end(), tight.begin(), tight.end());
          std::sort(same->faces.begin(), same->faces.end());
          same->faces.erase(std::unique(same->faces.begin(), same->faces.end()), same->faces.end());
        }
      }
  // Each blend's end in one corner, each pair of faces not next to each other in two or none,
  // and one corner more than such pairs: the corners part the faces' ring.
  const char* parted = "the corner where more than three blended edges meet does not part into "
                       "balls on faces next to each other";
  const auto  pairOf = [](int a, int b) { return std::make_pair(std::min(a, b), std::max(a, b)); };
  const auto  next   = [&](const std::pair<int, int>& pair) {
    return pair.second == pair.first + 1 || (pair.first == 0 && pair.second == count - 1);
  };
  std::map<std::pair<int, int>, std::vector<int>> sharedBy;
  for (int c = 0; c < int(cells.size()); ++c)
  {
    const auto& f = cells[c].faces;
    if (f.size() < 3)
      refuse(parted);
    for (size_t m = 0; m < f.size(); ++m)
      sharedBy[pairOf(f[m], f[(m + 1) % f.size()])].push_back(c);
  }
  int diagonals = 0;
  for (const auto& [pair, users] : sharedBy)
  {
    if (users.size() != (next(pair) ? 1u : 2u))
      refuse(parted);
    diagonals += !next(pair);
  }
  for (int i = 0; i < count; ++i)
    if (!sharedBy.count(pairOf(i, (i + 1) % count)))
      refuse(parted);
  if (int(cells.size()) != diagonals + 1)
    refuse(parted);
  // The round between each such pair: its rail the line where their planes meet, from the
  // place across from one ball's centre to the other's.
  std::map<std::pair<int, int>, int> roundOf;
  for (auto& [pair, users] : sharedBy)
  {
    if (next(pair))
      continue;
    const auto [i, j] = pair;
    gp_Vec d          = n[i].Crossed(n[j]);
    if (!(d.Magnitude() > tol::angle))
      refuse(parted);
    d.Normalize();
    double ta = gp_Vec(V, cells[users[0]].centre).Dot(d),
           tb = gp_Vec(V, cells[users[1]].centre).Dot(d);
    if (ta > tb)
    {
      std::swap(users[0], users[1]);
      std::swap(ta, tb);
    }
    const Handle(Geom_Line) line = new Geom_Line(V, gp_Dir(d));
    std::vector<OnSurface>  on;
    for (const int f : {i, j})
    {
      const auto [surface, location] = surfaceOf(ring[f]);
      const gp_Pln  plane            = GeomAdaptor_Surface(surface).Plane();
      const gp_Trsf back             = location.Inverted().Transformation();
      double        u0, v0, u1, v1;
      ElSLib::Parameters(plane, line->Value(ta).Transformed(back), u0, v0);
      ElSLib::Parameters(plane, line->Value(tb).Transformed(back), u1, v1);
      on.push_back({segment2d({u0, v0}, {u1, v1}, ta, tb), surface, location});
    }
    TopoDS_Vertex va, vb;
    BRep_Builder  builder;
    builder.MakeVertex(va, line->Value(ta), tol::length);
    builder.MakeVertex(vb, line->Value(tb), tol::length);
    if (pieces.size() >= pieces.capacity())
      refuse(parted);
    Piece piece;
    Rail& rail = piece.rail;
    rail.op    = op;
    rail.boundedFit = request.fitWork == FitWork::Bounded;
    rail.reads      = request.reads;
    rail.edge  = edgeBetween(line, va, ta, vb, tb, on, tol::length);
    rail.curve.Initialize(rail.edge);
    rail.s[0]     = sides[i];
    rail.s[1]     = sides[j];
    rail.sigma    = sigma;
    rail.scale    = scale;
    rail.first    = ta;
    rail.last     = tb;
    piece.apex    = vertex;
    roundOf[pair] = int(pieces.size());
    pieces.push_back(piece);
  }
  for (int c = 0; c < int(cells.size()); ++c)
  {
    const auto&                      f = cells[c].faces;
    std::vector<std::pair<int, int>> ends;
    for (size_t m = 0; m < f.size(); ++m)
    {
      const auto pair = pairOf(f[m], f[(m + 1) % f.size()]);
      if (next(pair))
        ends.push_back(after[pair.second == pair.first + 1 ? pair.first : count - 1]);
      else
        ends.push_back({roundOf[pair], sharedBy[pair][0] == c ? 0 : 1});
    }
    corner(ends);
  }
}

void Blend::carries()
{
  // Each piece's rail runs on past an end that a face stops, so that the blend reaches a face
  // that leans away from it; at a joint it ends with the edge.
  for (auto& piece : pieces)
  {
    Rail&      rail = piece.rail;
    const auto run  = [&](int end) {
      if (piece.joint[end].meets == Meets::Joint || piece.joint[end].meets == Meets::Corner)
        return 0.0;
      gp_Pnt p;
      gp_Vec d;
      rail.curve.D1(end == 0 ? rail.first : rail.last, p, d);
      // The blend's size at this end, where it changes along the edge (Rail::sizes).
      const double size = rail.sizes.places.empty()
                             ? rail.scale
                             : std::max(rail.s[0].distance, rail.s[1].distance)
                                * rail.factor(end == 0 ? rail.first : rail.last);
      // Twice the blend's size on, so that it reaches a face that leans away; a mitre's
      // blends may meet a leaning face further beyond the vertex: twice as far there.
      double length = (piece.joint[end].meets == Meets::Miter ? 4 : 2) * size;
      // At a notch (Attach::notch), where an edge of a support goes on straight from the
      // blended one, the blend may run on along it to the face at its far end: that far on.
      // Where it leaves at an angle to the rail's run on, the contact on that support, its
      // offset off the edge, crosses it that offset over the angle's tangent past the
      // vertex: that far on, within the edge. So too where a face stops the blend
      // (Attach::stops), whose contacts end where they cross the supports' edges there.
      // A chamfer's offset is its setback; a ball's is its radius over the half angle's
      // tangent between the faces.
      if (piece.joint[end].meets == Meets::Miter || piece.joint[end].meets == Meets::Stop)
      {
        TopoDS_Vertex ends[2];
        TopExp::Vertices(TopoDS::Edge(rail.edge.Oriented(TopAbs_FORWARD)), ends[0], ends[1]);
        gp_Vec mine, theirs, n0, n1;
        if (leaving(rail.edge, ends[end], mine, end)
            && normalAt(rail.s[0].face, rail.edge, ends[end], n0)
            && normalAt(rail.s[1].face, rail.edge, ends[end], n1))
        {
          const double ballOffset = size * std::tan(n0.Angle(n1) / 2);
          for (const auto& other : edgesOfVertex.FindFromKey(ends[end]))
          {
            if (other.IsSame(rail.edge) || !leaving(TopoDS::Edge(other), ends[end], theirs))
              continue;
            int on = -1;
            for (const auto& face : facesOfEdge.FindFromKey(other))
              for (int side = 0; side < 2; ++side)
                if (face.IsSame(rail.s[side].face))
                  on = side;
            if (on < 0)
              continue;
            const double offset =
              rail.law == Law::Chamfer ? rail.s[on].distance * size / rail.scale : ballOffset;
            const double reach =
              GCPnts_AbscissaPoint::Length(BRepAdaptor_Curve(TopoDS::Edge(other))) + 2 * size;
            const double angle = theirs.Angle(-mine);
            if (tangent(mine, theirs))
              length = std::max(length, reach);
            else if (angle < M_PI / 2)
              length = std::max(length, std::min(reach, offset / std::tan(angle) + 2 * size));
            // A mitre's partner on the same face: the contacts, offsets d and d' off the
            // two edges at an angle theta between them, cross (d' + d cos theta) / sin theta
            // past the vertex along this one, far beyond it at a narrow re-entrant notch.
            const auto& joint = piece.joint[end];
            if (joint.meets != Meets::Miter || !other.IsSame(pieces[joint.piece].rail.edge))
              continue;
            const Rail& partner = pieces[joint.piece].rail;
            gp_Vec      m0, m1;
            if (rail.law == Law::Chamfer
                || !normalAt(partner.s[0].face, partner.edge, ends[end], m0)
                || !normalAt(partner.s[1].face, partner.edge, ends[end], m1))
              continue;
            const double theta = mine.Angle(theirs);
            if (!(std::sin(theta) > tol::angle))
              continue;
            const double theirOffset = partner.scale * std::tan(m0.Angle(m1) / 2);
            length =
              std::max(length,
                       (theirOffset + offset * std::cos(theta)) / std::sin(theta) + 2 * size);
          }
        }
      }
      // Along a curved edge no further than it turns by a quarter of a turn: the stations'
      // planes across the rail turn with it, and past the end they need not meet the ball's
      // path (a split piece between two planes on an arc rolls straight on). The march ends
      // the carry where the ball strays from the rail (Stations::reached); an L corner's
      // ball rolling on round a concave blend goes a sixth of a turn to the face beyond.
      gp_Vec d2;
      rail.curve.D2(end == 0 ? rail.first : rail.last, p, d, d2);
      const double curvature = d.Crossed(d2).Magnitude() / std::pow(d.Magnitude(), 3);
      if (curvature > 0)
        length = std::min(length, M_PI / 2 / curvature);
      return length / d.Magnitude();
    };
    rail.from = rail.first - run(0);
    rail.to   = rail.last + run(1);
  }
}

//! Continue an already identified terminal branch to its restriction arcs. The temporary
//! tangent carrier is corrected on both exact offsets at every step; it is not an envelope.
//! Whether continueTerminal may carry a rail's ends on: a conic's of one size, split onto
//! other faces, off faces.
static bool continuesTerminal(const Rail& rail)
{
  return rail.parameterSupports && !rail.closed && rail.sizes.places.empty()
         && rail.law == Law::Conic && !rail.s[0].cliffCurve && !rail.s[1].cliffCurve;
}

static void continueTerminal(Piece& piece, Stations& stations,
                             const NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>& facesOfEdge,
                             const NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>& edgesOfVertex)
{
  Rail& rail = piece.rail;
  if (!continuesTerminal(rail))
    return;
  TopoDS_Vertex vertices[2];
  TopExp::Vertices(TopoDS::Edge(rail.edge.Oriented(TopAbs_FORWARD)), vertices[0], vertices[1]);
  for (int end = 0; end < 2; ++end)
  {
    if (piece.joint[end].meets != Meets::Stop || vertices[end].IsNull())
      continue;
    TopoDS_Edge arcs[2];
    TopoDS_Face stop;
    for (const auto& edge : edgesOfVertex.FindFromKey(vertices[end]))
      for (const auto& face : facesOfEdge.FindFromKey(edge))
      {
        if (face.IsSame((*rail.parameterSupports)[0].face)
            || face.IsSame((*rail.parameterSupports)[1].face)
            || face.IsSame(rail.s[0].face) || face.IsSame(rail.s[1].face))
          continue;
        TopoDS_Edge found[2];
        for (int side = 0; side < 2; ++side)
          for (TopExp_Explorer e(rail.s[side].face, TopAbs_EDGE); e.More(); e.Next())
            for (const auto& adjacent : facesOfEdge.FindFromKey(e.Current()))
              if (adjacent.IsSame(face))
                found[side] = TopoDS::Edge(e.Current());
        if (!found[0].IsNull() && !found[1].IsNull())
        {
          stop = TopoDS::Face(face);
          arcs[0] = found[0]; arcs[1] = found[1];
        }
      }
    if (stop.IsNull())
      continue;
    const auto surface = BRep_Tool::Surface(stop);
    const auto distance = [&](const Station& s, int side) {
      GeomAPI_ProjectPointOnSurf onto(s.foot[side], surface);
      if (onto.NbPoints() == 0) return std::numeric_limits<double>::quiet_NaN();
      double u,v; onto.LowerDistanceParameters(u,v);
      gp_Pnt p; gp_Vec du,dv; surface->D1(u,v,p,du,dv);
      const gp_Vec n = du.Crossed(dv);
      return n.Magnitude() > 0 ? gp_Vec(p,s.foot[side]).Dot(n.Normalized())
                              : std::numeric_limits<double>::quiet_NaN();
    };
    const double start = end == 0 ? rail.from : rail.to;
    Station current;
    if (!stations.at(start,current)) continue;
    double before[2] = {distance(current,0),distance(current,1)};
    // March in length units; preserve orientation of the original parameter.
    const double sign = end == 0 ? -1 : 1;
    const double step = rail.scale / 16;
    const auto original = rail.terminalCarrier[end];
    std::vector<Station> trace{current};
    bool crossed[2] = {rail.s[0].face.IsSame((*rail.parameterSupports)[0].face),
                       rail.s[1].face.IsSame((*rail.parameterSupports)[1].face)};
    std::optional<Piece::TerminalContact> witnesses[2];
    // A finite exploration policy; reaching it is not a geometric nonexistence claim.
    for (int count = 0; count < 256 && !(crossed[0] && crossed[1]); ++count)
    {
      gp_Vec tangent = current.along;
      gp_Pnt p; gp_Vec oldDirection;
      rail.stationD1(current.t,p,oldDirection);
      if (tangent.Dot(oldDirection) < 0) tangent.Reverse();
      tangent.Normalize();
      bool solved = false;
      Station next;
      double h = step;
      for (int cut = 0; cut < 24; ++cut, h /= 2)
      {
        const double t = current.t + sign*h;
        const gp_Pnt origin = current.centre.Translated(tangent * (-current.t));
        const Handle(Geom_Curve) line = new Geom_Line(origin,gp_Dir(tangent));
        // The finite restriction keeps the carrier's active side explicit.
        rail.terminalCarrier[end] = new Geom_TrimmedCurve(line,
                                     std::min(start,t),std::max(start,t));
        next = current;
        if (contactAt(rail,t,next) && next.along.Dot(current.along) > 0
            && next.centre.Distance(current.centre) <= 2*h)
        { solved = true; break; }
      }
      if (!solved) break;
      double after[2] = {distance(next,0),distance(next,1)};
      for (int side = 0; side < 2; ++side)
        if (!crossed[side] && before[side]*after[side] <= 0)
        {
          Station a = current, b = next, root;
          double da = before[side];
          for (int iteration = 0; iteration < 64; ++iteration)
          {
            const double t = (a.t+b.t)/2;
            root = a;
            if (!contactAt(rail,t,root)) break;
            const double d = distance(root,side);
            if (std::abs(d) <= rail.resolved(limit::solved)) break;
            if (da*d <= 0) b=root;
            else { a=root; da=d; }
          }
          double f,l; const auto arc = BRep_Tool::Curve(arcs[side],f,l);
          GeomAPI_ProjectPointOnCurve on(root.foot[side],arc,f,l);
          crossed[side] = std::abs(distance(root,side)) <= rail.resolved(limit::solved)
                          && on.NbPoints() && on.LowerDistance() <= rail.resolved(limit::together);
          if (crossed[side])
            witnesses[side] = Piece::TerminalContact{
              rail.s[side].face, arcs[side], on.LowerDistanceParameter(), root.t};
        }
      trace.push_back(next); current = next;
      before[0]=after[0];before[1]=after[1];
    }
    rail.terminalCarrier[end] = original;
    if (!(crossed[0] && crossed[1]) || trace.size() < 3) continue;
    if (end == 0) std::reverse(trace.begin(),trace.end());
    Handle(NCollection_HArray1<gp_Pnt>) points = new NCollection_HArray1<gp_Pnt>(1,int(trace.size()));
    Handle(NCollection_HArray1<double>) parameters = new NCollection_HArray1<double>(1,int(trace.size()));
    for (int i=0;i<int(trace.size());++i)
    { points->SetValue(i+1,trace[i].centre); parameters->SetValue(i+1,trace[i].t); }
    GeomAPI_Interpolate fit(points,parameters,false,rail.resolved(limit::solved));
    fit.Perform();
    if (!fit.IsDone()) continue;
    rail.terminalCarrier[end] = fit.Curve();
    for (int side = 0; side < 2; ++side)
      piece.terminalContact[end][side] = witnesses[side];
    if (end == 0) rail.from = trace.front().t; else rail.to = trace.back().t;
  }
}

//! A chordal blend whose ball keeps one size along the edge (between planes, or round a rim
//! whose wedge does not change) is that round: made as one, exact where a round is.
static void roundOfChord(Rail& rail, Stations& stations)
{
  constexpr int samples = 9;
  double        least = Precision::Infinite(), most = 0;
  for (int k = 0; k < samples; ++k)
  {
    Station s;
    if (!stations.at(rail.t(double(k) / (samples - 1)), s))
      return;
    const double r = s.centre.Distance(s.foot[0]);
    least          = std::min(least, r);
    most           = std::max(most, r);
  }
  if (most - least > limit::agree * most)
    return;
  rail.law           = Law::Conic;
  rail.s[0].distance = rail.s[1].distance = rail.scale = (least + most) / 2;
}

//! Whether a full round has the radius `r` all along (Request::meet).
static void meetsAt(const Rail& rail, Stations& stations, double r)
{
  constexpr int samples = 9;
  for (int k = 0; k < samples; ++k)
  {
    Station s;
    if (!stations.at(rail.t(double(k) / (samples - 1)), s)
        || std::abs(s.centre.Distance(s.foot[0]) - r)
             > tol::resolved(limit::together, r, rail.magnitude))
      fail(rail.op,
           no_solution,
           "the rounds' contacts on the face between do not meet along it",
           {rail.edge});
  }
}

//! Whether a contact may cross a seam of its face between its ends: a seam along a line of the
//! face's parameters, at one of them, is crossed only where the contact's poles reach past its
//! place in some period inside their range (a contact from seam to seam, round a hole, does not
//! cross it between); any other seam may be.
static bool acrossSeam(const Handle(Geom2d_BSplineCurve)& contact,
                       const TopoDS_Edge&                 seam,
                       const TopoDS_Face&                 face,
                       const Handle(Geom_Surface)&        surface)
{
  double     f, l;
  const auto pcurve = BRep_Tool::CurveOnSurface(seam, face, f, l);
  const auto line   = Handle(Geom2d_Line)::DownCast(pcurve);
  if (line.IsNull())
    return true;
  const gp_Dir2d d      = line->Direction();
  const int      across = std::abs(d.X()) < tol::angle ? 0 : std::abs(d.Y()) < tol::angle ? 1 : -1;
  if (across < 0)
    return true;
  const double period = across == 0 ? (surface->IsUPeriodic() ? surface->UPeriod() : 0)
                                    : (surface->IsVPeriodic() ? surface->VPeriod() : 0);
  if (!(period > 0))
    return true;
  double lo = HUGE_VAL, hi = -HUGE_VAL;
  for (int i = 1; i <= contact->NbPoles(); ++i)
  {
    const double x = across == 0 ? contact->Pole(i).X() : contact->Pole(i).Y();
    lo             = std::min(lo, x);
    hi             = std::max(hi, x);
  }
  const double at    = across == 0 ? line->Location().X() : line->Location().Y(),
               slack = 1e-9 * period;
  const double first = at + std::ceil((lo + slack - at) / period) * period;
  return first < hi - slack;
}

bool Blend::continueWholeSupport(int k, int side, Stations& stations)
{
  Piece& piece = pieces[k];
  if (piece.rail.parameterSupports)
    return false;
  const double u0 = piece.uAt(0), u1 = piece.uAt(1);
  // A contact may have crossed a smooth support boundary over the whole rail.
  // There is then no longitudinal crossing at which to split the piece. Keep
  // the original station frame, transfer its middle ball to the adjacent
  // support and march the same offset equations on that support instead.
  bool        transferred = false;
  const auto& beyondEdge  = farOf(piece.rail, side);
  if (beyondEdge)
    for (const auto& adjacent : facesOfEdge.FindFromKey(beyondEdge->edge))
    {
      if (adjacent.IsSame(piece.rail.s[side].face) || adjacent.IsSame(piece.rail.s[1 - side].face))
        continue;
      Side next    = piece.rail.s[side];
      next.face    = TopoDS::Face(adjacent);
      next.surface = BRep_Tool::Surface(next.face);
      next.outward = adjacent.Orientation() == TopAbs_REVERSED ? -1 : 1;
      double     first, last, a, b;
      const auto here =
        BRep_Tool::CurveOnSurface(beyondEdge->edge, piece.rail.s[side].face, first, last);
      const auto there = BRep_Tool::CurveOnSurface(beyondEdge->edge, next.face, a, b);
      if (here.IsNull() || there.IsNull())
        continue;
      bool smooth = true;
      for (double share : {0.1, 0.5, 0.9})
      {
        const double   t  = first + (last - first) * share;
        const gp_Pnt2d uv = here->Value(t), xy = there->Value(t);
        Jet            before, after;
        smooth = smooth && jetAt(piece.rail.s[side], uv.X(), uv.Y(), before)
                 && jetAt(next, xy.X(), xy.Y(), after)
                 && SectionPolicy::tangentFor(before.n.Angle(after.n), piece.rail.scale);
      }
      if (!smooth)
        continue;
      Station seed;
      if (!stations.at(piece.rail.t((u0 + u1) / 2), seed))
        continue;
      GeomAPI_ProjectPointOnSurf onto(seed.foot[side], next.surface);
      if (onto.NbPoints() == 0)
        continue;
      onto.LowerDistanceParameters(seed.uv[side][0], seed.uv[side][1]);
      // Projection is only a Newton seed. The offset equations must put the
      // contact on the finite destination, rather than its carrier extension.
      Rail transferredRail = piece.rail;
      transferredRail.parameterSupports = std::array<Side, 2>{piece.rail.s[0], piece.rail.s[1]};
      transferredRail.s[side] = next;
      if (!contactAt(transferredRail, seed.t, seed)
          || inside(next.face).out(gp_Pnt2d(seed.uv[side][0], seed.uv[side][1])))
        continue;
      bool consumed = true;
      for (TopExp_Explorer edge(piece.rail.s[side].face, TopAbs_EDGE); edge.More(); edge.Next())
        consumed =
          consumed
          && (edge.Current().IsSame(piece.rail.edge) || edge.Current().IsSame(beyondEdge->edge));
      if (consumed)
        piece.rail.consumedSupports.push_back(piece.rail.s[side].face);
      piece.rail.parameterSupports = std::array<Side, 2>{piece.rail.s[0], piece.rail.s[1]};
      piece.rail.s[side]           = next;
      piece.rail.wholeSupportContinuation = true;
      seeds[k]                     = seed;
      transferred                  = true;
      break;
    }
  return transferred;
}

namespace
{
//! Where a piece's contacts are tested against their faces: at this many equal steps of its fit
//! (chosen: fine enough to find each edge a contact crosses; every 16th for a whole contact).
constexpr int contactSamples = 256;
} // namespace

void Blend::shareHalves()
{
  // Sections in two spans along a whole chain, or along none of it: its sections share their
  // form at its joints.
  for (bool changed = true; changed;)
  {
    changed = false;
    for (auto& piece : pieces)
      for (const auto& joint : piece.joint)
        if ((joint.meets == Meets::Joint || joint.meets == Meets::Split)
            && pieces[joint.piece].rail.halves != piece.rail.halves)
        {
          pieces[joint.piece].rail.halves = piece.rail.halves = true;
          changed                                             = true;
        }
  }
}

bool Blend::capFoldsOf(int k, Dips& dips, Folds& folds)
{
  std::vector<Fold> across;
  auto              found = capFolds(dips, &across);
  const bool        taken = !found.empty() || !across.empty();
  if (!found.empty())
    folds.emplace_back(k, std::move(found));
  if (!across.empty())
    open.emplace_back(k, std::move(across));
  return !taken;
}

bool Blend::convexOver(const Rail& rail, int s2, const TopoDS_Edge& edge, const Meeting& at) const
{
  for (const auto& face : facesOfEdge.FindFromKey(edge))
  {
    if (face.IsSame(rail.s[s2].face))
      continue;
    double     f, l;
    const auto pc = BRep_Tool::CurveOnSurface(edge, TopoDS::Face(face), f, l);
    if (pc.IsNull())
      return false;
    Side beyond    = rail.s[s2];
    beyond.face    = TopoDS::Face(face);
    beyond.surface = BRep_Tool::Surface(beyond.face);
    beyond.outward = face.Orientation() == TopAbs_REVERSED ? -1 : 1;
    beyond.cliffCurve.Nullify();
    const gp_Pnt2d there = pc->Value(at.second);
    Jet            b;
    if (!jetAt(beyond, there.X(), there.Y(), b))
      return false;
    gp_Pnt inside;
    rail.curve.D0(rail.t(at.first), inside);
    return gp_Vec(b.p, inside).Dot(b.n) < 0;
  }
  return false;
}

bool Blend::centreMends(int k)
{
  const Piece& piece = pieces[k];
  if (!mayCentre(k))
    return false;
  Rail rail     = piece.rail;
  rail.byCentre = true;
  Stations stations(rail);
  try
  {
    stations.march(seeds[k]);
    const Fitted fitted = fit(rail, stations);
    Dips         centred{rail, fitted};
    return !runsBack(centred);
  }
  catch (const Failure& failure)
  {
    if (failure.reason == std::string("cancelled"))
      throw;
    return false;
  }
}

bool Blend::settleFolds(int k, Folds& folds, bool onCliff)
{
  Dips dips{pieces[k].rail, pieces[k].fitted};
  if (onCliff || capFoldsOf(k, dips, folds))
  {
    if (runsBack(dips) && centreMends(k))
    {
      centreChain(k);
      return false;
    }
    checkFold(dips);
  }
  return true;
}

void Blend::prepareRail(Piece& piece, Stations& stations)
{
  if (request.intersectSpines)
    prepareConicWindow(piece.rail);
  if (piece.rail.law == Law::Chordal)
    roundOfChord(piece.rail, stations);
  if (request.meet > 0)
    meetsAt(piece.rail, stations, request.meet);
  // A carry that stopped short ends where the ball last went on with the rail.
  if (const auto reached = stations.reached(0))
    piece.rail.from = std::max(piece.rail.from, *reached);
  if (const auto reached = stations.reached(1))
    piece.rail.to = std::min(piece.rail.to, *reached);
  if (request.intersectSpines)
    continueTerminal(piece, stations, facesOfEdge, edgesOfVertex);
  else if (request.reads && continuesTerminal(piece.rail))
    request.reads->spines.push_back([] { return true; });
}

bool Blend::fitSheet(int                                         k,
                     int                                         attempt,
                     Stations&                                   stations,
                     SheetFit&                                   out,
                     std::vector<std::shared_ptr<const Crease>>& creases)
{
  Piece& piece = pieces[k];
  // An exact sheet needs no fit: its contacts are its own.
  // A conic's corner patch takes its weights from fitted sheets: none extruded there (a
  // chamfer's corner takes none). The corner of more than three conics other than rounds is
  // split from their extruded sheets (Attach::splitCorner).
  const bool round   = piece.rail.law == Law::Conic && std::abs(piece.rail.rho - 0.5) < limit::equal
                       && std::abs(piece.rail.s[0].distance - piece.rail.s[1].distance)
                            < limit::equal * piece.rail.scale;
  const auto patched = [&](int end) {
    return piece.joint[end].meets == Meets::Corner
           && (round || corners[piece.joint[end].piece].ends.size() <= 3);
  };
  const bool corner = piece.rail.law != Law::Chamfer && (patched(0) || patched(1));
  auto&      exact  = out.exact;
  exact             = exactSheet(piece.rail, stations);
  if (!exact && !corner)
    exact = extrudedSheet(piece.rail, stations);
  out.contacts         = exact ? exactContacts(piece.rail, stations, *exact) : std::nullopt;
  const auto& contacts = out.contacts;
  // Locate support crossings from the ball's contacts before fitting the G2 section.
  // A finite spline support's tangent continuation is only G1: its curvature jumps,
  // so fitting a curvature-dependent section across it before splitting cannot converge.
  // The circular pilot uses the same balls and contact parameters. It is always replaced
  // by the requested G2 section after crossing discovery, and never reaches attachment.
  out.contactPilot = !contacts && attempt == 0 && piece.rail.law == Law::G2
                     && !piece.rail.s[0].cliffCurve && !piece.rail.s[1].cliffCurve
                     && piece.rail.s[0].trim.IsNull() && piece.rail.s[1].trim.IsNull();
  if (contacts)
  {
    piece.fitted = *contacts;
    piece.sheet  = *exact;
  }
  else
  {
    try
    {
      piece.fitted = fit(piece.rail, stations, out.contactPilot);
    }
    catch (const Failure& failure)
    {
      if (failure.reason == std::string("cancelled") || !centreMends(k))
        throw;
      centreChain(k);
      return false;
    }
    piece.sheet = exact ? *exact : Sheet{piece.fitted.blend};
    if (!exact)
      creases = checkCrossing(piece.rail, stations, piece.fitted);
  }
  return true;
}

Blend::Next Blend::offFace(int k, Stations& stations, double u0, double u1)
{
  Piece& piece = pieces[k];
  // A contact wholly off its face: on both, the blend rests on neither and does not fit;
  // on one, the ball may roll over onto the next face, which the blend does not do yet.
  int off = 0, offSide = -1;
  for (int s2 = 0; s2 < 2; ++s2)
  {
    Inside& onFace = inside(piece.rail.s[s2].face);
    // On the face's far edge all along counts as off it: the ball rests on the face and that
    // edge at once (a side exactly as high as the ball is wide), and the blend uses the face
    // up as where the side is narrower.
    const auto&       far = farOf(piece.rail, s2);
    BRepAdaptor_Curve farCurve;
    if (far)
      farCurve.Initialize(far->edge);
    const auto onFar = [&](double u) {
      if (!far)
        return false;
      const gp_Pnt2d      uv = piece.fitted.contact[s2]->Value(u);
      const gp_Pnt        p  = footAt(piece.rail.s[s2], uv.X(), uv.Y());
      const Extrema_ExtPC to(p, farCurve);
      if (!to.IsDone() || to.NbExt() == 0)
        return false;
      double nearest = HUGE_VAL;
      for (int i = 1; i <= to.NbExt(); ++i)
        nearest = std::min(nearest, to.SquareDistance(i));
      return std::sqrt(nearest) <= piece.rail.resolved(limit::together);
    };
    const auto out    = [&](double u) { return onFace.out(piece.fitted.contact[s2]->Value(u)); };
    bool       wholly = true, along = true;
    for (int i = 0; i <= contactSamples && (wholly || along); i += contactSamples / 16)
    {
      wholly = wholly && out(u0 + (u1 - u0) * i / contactSamples);
      along  = along && onFar(u0 + (u1 - u0) * i / contactSamples);
    }
    if (!wholly && !along)
      continue;
    ++off;
    offSide = s2;
    // A flat face all of whose edges are blended: the blends use the cap up.
    if (GeomAdaptor_Surface(piece.rail.s[s2].surface).GetType() != GeomAbs_Plane)
      continue;
    bool all = true;
    for (TopExp_Explorer e(piece.rail.s[s2].face, TopAbs_EDGE); e.More() && all; e.Next())
      all = std::any_of(pieces.begin(), pieces.end(), [&](const Piece& other) {
        return other.rail.edge.IsSame(e.Current());
      });
    if (all)
      fail(op,
           no_solution,
           "the blend consumes the cap: the face beside its edge vanishes",
           {piece.rail.edge},
           FailureKind::CapConsumed);
  }
  if (off == 2)
    fail(op,
         no_solution,
         "the blend does not fit: it is wider than both faces beside its edge",
         {piece.rail.edge},
         FailureKind::Overflow);
  if (off == 1)
  {
    // A face narrower than the blend: the ball rests on its far edge (a cliff).
    // Between stops at both ends, or round a closed edge (a band narrower than the blend).
    const bool stopped = piece.joint[0].meets == Meets::Stop && piece.joint[1].meets == Meets::Stop;
    // Not on an edge blended too: it is not there for the ball to rest on, nor to trim by.
    const auto blended = [&](const TopoDS_Edge& edge) {
      return std::any_of(pieces.begin(), pieces.end(), [&](const Piece& other) {
        return other.rail.edge.IsSame(edge);
      });
    };
    const auto& far = farOf(piece.rail, offSide);
    if (far && !blended(far->edge) && (stopped || piece.rail.closed)
        && cliffOf(piece.rail, offSide, *far))
    {
      return Next::Retry;
    }
    // Else the face beyond the far edge trims the blend (a notch), between stops or round a
    // closed edge.
    if (far && !blended(far->edge) && (stopped || piece.rail.closed))
    {
      TopoDS_Face beyond;
      for (const auto& face : facesOfEdge.FindFromKey(far->edge))
        if (!face.IsSame(piece.rail.s[offSide].face))
          beyond = TopoDS::Face(face);
      if (!beyond.IsNull() && !beyond.IsSame(piece.rail.s[1 - offSide].face))
      {
        piece.rail.s[offSide].trim       = beyond;
        piece.rail.s[offSide].trimEdge   = far->edge;
        piece.rail.s[offSide].cliffWhole = far->whole;
        return Next::Retry;
      }
    }
    // A part of a split, wholly past the place its contact leaves the face: that is seen
    // along the whole edge, which came first.
    if (piece.rail.parameterSupports)
      off = 0;
    else if ((piece.joint[0].meets == Meets::Stop) != (piece.joint[1].meets == Meets::Stop)
             && piece.joint[piece.joint[0].meets == Meets::Stop ? 1 : 0].meets == Meets::Joint)
    {
      ending.push_back(TopoDS::Edge(piece.rail.edge));
      again = true;
      return Next::Stop;
    }
    else
    {
      if (continueWholeSupport(k, offSide, stations))
      {
        return Next::Retry;
      }
      fail(op,
           unsupported_geometry,
           "the blend's contact runs off the face beside the edge, onto faces the blend does "
           "not roll it over yet",
           {piece.rail.edge});
    }
  }
  return Next::Go;
}

bool Blend::smoothOver(const Piece&       piece,
                       int                s2,
                       const TopoDS_Edge& over,
                       const Meeting&     found) const
{
  // Over an edge where the faces are not tangent the contact meets what ends the
  // blend: near an end of the piece, that end's stopping face trims it there.
  bool smooth = false;
  for (const auto& next : facesOfEdge.FindFromKey(over))
  {
    if (next.IsSame(piece.rail.s[s2].face))
      continue;
    double     f, l;
    const auto pc = BRep_Tool::CurveOnSurface(over, TopoDS::Face(next), f, l);
    if (pc.IsNull())
      continue;
    Side beyond          = piece.rail.s[s2];
    beyond.face          = TopoDS::Face(next);
    beyond.surface       = BRep_Tool::Surface(beyond.face);
    beyond.outward       = next.Orientation() == TopAbs_REVERSED ? -1 : 1;
    const gp_Pnt2d here  = pc->Value(found.second),
                   there = piece.fitted.contact[s2]->Value(found.first);
    Jet            a, b;
    smooth = jetAt(beyond, here.X(), here.Y(), a)
             && jetAt(piece.rail.s[s2], there.X(), there.Y(), b)
             && SectionPolicy::tangentFor(a.n.Angle(b.n), piece.rail.scale);
  }
  return smooth;
}

void Blend::refuseFeature(const Piece& piece, int s2, const TopoDS_Edge& over) const
{
  // Into a hole of its face, a feature that can be lifted off and put back
  // (liftFeatures): the blend runs into its walls.
  for (TopoDS_Iterator w(TopoDS::Face(piece.rail.s[s2].face.Oriented(TopAbs_FORWARD))); w.More();
       w.Next())
  {
    bool holds = false;
    for (TopoDS_Iterator x(w.Value()); x.More() && !holds; x.Next())
      holds = x.Value().IsSame(over);
    if (!holds || !holeOf(w.Value(), TopoDS::Face(piece.rail.s[s2].face.Oriented(TopAbs_FORWARD))))
      continue;
    std::vector<TopoDS_Shape> walls;
    for (const auto& next : facesOfEdge.FindFromKey(over))
      if (!next.IsSame(piece.rail.s[s2].face))
        walls.push_back(next);
    // Only a feature that is put back as it is: a hole (its rim convex)
    // beside a convex blend, a boss (its foot concave) beside a concave one.
    if (walls.size() != 1
        || convexity(op, over, piece.rail.s[s2].face, TopoDS::Face(walls.front()), piece.rail.scale)
             != (piece.rail.sigma < 0 ? -1 : 1))
      continue;
    const auto feature = featuresOf(body, walls);
    if (std::none_of(feature.begin(), feature.end(), [&](const TopoDS_Shape& f) {
          return std::any_of(pieces.begin(), pieces.end(), [&](const Piece& p) {
            return f.IsSame(p.rail.s[0].face) || f.IsSame(p.rail.s[1].face);
          });
        }))
      throw Collision(op,
                      "the blend's contact runs into a hole or a boss in its face",
                      piece.rail.edge,
                      walls);
  }
}

bool Blend::stopsNearEnd(const Piece& piece, const TopoDS_Edge& over, const Meeting& found) const
{
  // An edge of the face at an end vertex of the blended edge is what stops the
  // blend there, however far from the end the contact crosses it (at a wedge of
  // angle w, its offset over tan w); else near an end, within three radii.
  const double reach = 3 * piece.rail.scale;
  gp_Pnt       p0, p1, pc;
  piece.rail.curve.D0(piece.rail.first, p0);
  piece.rail.curve.D0(piece.rail.last, p1);
  piece.rail.curve.D0(piece.rail.t(found.first), pc);
  // Or over the edge onto a face that stops the blend at an end of its edge
  // (a top narrowing to an oblique end): that stop trims it there,
  // wherever along the edge the contact reaches it.
  bool          ontoStop = false;
  TopoDS_Vertex ends[2];
  TopExp::Vertices(TopoDS::Edge(piece.rail.edge.Oriented(TopAbs_FORWARD)), ends[0], ends[1]);
  const bool atEnd = std::any_of(std::begin(ends), std::end(ends), [&](const TopoDS_Vertex& v) {
    const auto& at = edgesOfVertex.FindFromKey(v);
    return std::any_of(at.begin(), at.end(), [&](const TopoDS_Shape& edge) {
      return edge.IsSame(over);
    });
  });
  for (const auto& next : facesOfEdge.FindFromKey(over))
    for (int end = 0;
         end < 2 && !next.IsSame(piece.rail.s[0].face) && !next.IsSame(piece.rail.s[1].face);
         ++end)
    {
      const TopoDS_Vertex& v = ends[end];
      if (piece.joint[end].meets != Meets::Stop || v.IsNull() || !edgesOfVertex.Contains(v))
        continue;
      for (const auto& at : edgesOfVertex.FindFromKey(v))
        for (const auto& face : facesOfEdge.FindFromKey(at))
          ontoStop = ontoStop || face.IsSame(next);
    }
  return atEnd || ontoStop || pc.Distance(p0) <= reach || pc.Distance(p1) <= reach;
}

void Blend::crossingOverEdge(int k, double u0, double u1, Crossing& c)
{
  Piece&                           piece    = pieces[k];
  int&                             side     = c.side;
  double&                          at       = c.at;
  double&                          overAt   = c.overAt;
  bool&                            entering = c.entering;
  bool&                            notched  = c.notched;
  TopoDS_Edge&                     over     = c.over;
  std::optional<std::vector<Fold>> loops; // where a contact crosses itself, found once
  for (int s2 = 0; s2 < 2 && side < 0; ++s2)
  {
    Inside&    onFace    = inside(piece.rail.s[s2].face);
    const auto out       = [&](double u) { return onFace.out(piece.fitted.contact[s2]->Value(u)); };
    bool       startsOut = out(u0);
    double     previous  = u0;
    bool       trimmed   = false;
    for (int i = 1; i <= contactSamples && side < 0; ++i)
    {
      trimmed        = false;
      const double u = u0 + (u1 - u0) * i / contactSamples;
      // ON is a tolerance band, not a side of the boundary. Retain
      // the last definite side so a delayed exit cannot discard
      // the geometric crossing while traversing that band.
      const auto state = onFace.classify(piece.fitted.contact[s2]->Value(u));
      if (state == TopAbs_ON)
        continue;
      if ((state == TopAbs_OUT) == startsOut)
      {
        previous = u;
        continue;
      }
      // Keep the scan's forward lookahead, with its lower bound
      // still on the definite side before the boundary band.
      const double past = std::min(u1, u + (u1 - u0) / contactSamples);
      for (TopExp_Explorer e(piece.rail.s[s2].face, TopAbs_EDGE); e.More(); e.Next())
      {
        // A seam has the face on both sides: crossing it, the contact stays on it.
        if (e.Current().IsSame(piece.rail.edge)
            || BRep_Tool::IsClosed(TopoDS::Edge(e.Current()), piece.rail.s[s2].face))
          continue;
        // Over another blended edge the contact runs past the corner or mitre where
        // the other blend's contact trims it: not a face it rolls on.
        if (std::any_of(pieces.begin(), pieces.end(), [&](const Piece& other) {
              return other.rail.edge.IsSame(e.Current());
            }))
        {
          if (crossing(piece.fitted.contact[s2],
                       previous,
                       past,
                       previous,
                       TopoDS::Edge(e.Current()),
                       piece.rail.s[s2].face))
          {
            trimmed = true;
            break;
          }
          continue;
        }
        const auto found = crossing(piece.fitted.contact[s2],
                                    previous,
                                    past,
                                    previous,
                                    TopoDS::Edge(e.Current()),
                                    piece.rail.s[s2].face);
        if (!found)
          continue;
        // At an end of the piece: the place a split already made it end.
        const auto placeAt = [&](double v) {
          const gp_Pnt2d uv = piece.fitted.contact[s2]->Value(v);
          return supportPoint(piece.rail.s[s2].surface, uv.X(), uv.Y());
        };
        if (placeAt(found->first).Distance(placeAt(u0)) <= piece.rail.resolved(limit::together)
            || placeAt(found->first).Distance(placeAt(u1)) <= piece.rail.resolved(limit::together))
        {
          trimmed = true;
          break;
        }
        if (!smoothOver(piece, s2, TopoDS::Edge(e.Current()), *found))
        {
          refuseFeature(piece, s2, TopoDS::Edge(e.Current()));
          if (stopsNearEnd(piece, TopoDS::Edge(e.Current()), *found))
          {
            trimmed = true;
            break;
          }
          // Partway over a convex edge the blend goes on as it is on the face carried

          if (!convexOver(piece.rail, s2, TopoDS::Edge(e.Current()), *found))
            fail(op,
                 unsupported_geometry,
                 "a blend's contact runs off its face partway along the edge, which the blend "
                 "does not do yet",
                 {piece.rail.edge});
          notched = true;
        }
        side     = s2;
        at       = found->first;
        overAt   = found->second;
        over     = TopoDS::Edge(e.Current());
        entering = startsOut;
        break;
      }
      if (side < 0 && !trimmed)
      {
        // Inside the loop where the contact crosses itself it may leave its face over the
        // blended edge itself (the loop's tail past the edge): that part of the contact is
        // cut away with the loop (capFolds), and the scan goes on past the loop's end.
        if (!loops)
        {
          Dips dips{piece.rail, piece.fitted};
          loops = capFolds(dips);
        }
        const auto loop = std::find_if(loops->begin(), loops->end(), [&](const Fold& f) {
          return f.side == s2 && f.from < u && u < f.to;
        });
        if (loop == loops->end())
          fail(op,
               unsupported_geometry,
               "a blend's contact leaves its face where the blend does not find the edge it "
               "crosses",
               {piece.rail.edge});
        while (i < contactSamples && u0 + (u1 - u0) * i / contactSamples <= loop->to)
          ++i;
        previous  = u0 + (u1 - u0) * i / contactSamples;
        startsOut = out(previous);
        continue;
      }
      if (trimmed)
      {
        // A stop at one end may be followed by a support crossing at the other.
        // Continue from this side of the ignored transition, rather than dropping
        // the rest of the contact after its first entry into the face.
        startsOut = out(u);
        previous  = u;
      }
    }
  }
}

void Blend::crossingOverSeam(int k, double u0, double u1, Crossing& c)
{
  const Piece& piece     = pieces[k];
  int&         side      = c.side;
  double&      at        = c.at;
  double&      overAt    = c.overAt;
  bool&        entering  = c.entering;
  TopoDS_Edge& over      = c.over;
  bool&        overSeam  = c.overSeam;
  gp_Vec2d&    seamShift = c.seamShift;
  for (int s2 = 0; s2 < 2 && side < 0; ++s2)
  {
    const TopoDS_Face&          face    = piece.rail.s[s2].face;
    const Handle(Geom_Surface)& surface = piece.rail.s[s2].surface;
    // Only a face on a closed surface has a seam.
    if (!surface->IsUClosed() && !surface->IsVClosed())
      continue;
    const auto pointAt = [&](double u) {
      const gp_Pnt2d uv = piece.fitted.contact[s2]->Value(u);
      return surface->Value(uv.X(), uv.Y());
    };
    for (TopExp_Explorer e(face, TopAbs_EDGE); e.More() && side < 0; e.Next())
    {
      const TopoDS_Edge seam = TopoDS::Edge(e.Current());
      if (!BRep_Tool::IsClosed(seam, face)
          || !acrossSeam(piece.fitted.contact[s2], seam, face, surface)
          || poleInStrip(piece, s2, seam))
        continue;
      const auto found = crossing(piece.fitted.contact[s2], u0, u1, (u0 + u1) / 2, seam, face);
      if (!found)
        continue;
      const gp_Pnt there = pointAt(found->first);
      const double near  = piece.rail.resolved(limit::together);
      if (there.Distance(pointAt(u0)) <= near || there.Distance(pointAt(u1)) <= near)
        continue;
      // The part that lies off the face's own parameters, before the crossing or
      // after it, is the one carried a period over onto them.
      double ua, ub, va, vb;
      BRepTools::UVBounds(face, ua, ub, va, vb);
      // Each parameter's slack a share of its own range: an angle's and a length's
      // ranges differ by the face's size.
      const double su = 1e-9 * (ub - ua), sv = 1e-9 * (vb - va);
      const auto   within = [&](const gp_Pnt2d& p) {
        return p.X() >= ua - su && p.X() <= ub + su && p.Y() >= va - sv && p.Y() <= vb + sv;
      };
      const gp_Pnt2d ahead =
        piece.fitted.contact[s2]->Value(found->first + (u1 - found->first) * 1e-3);
      const gp_Pnt2d behind =
        piece.fitted.contact[s2]->Value(found->first - (found->first - u0) * 1e-3);
      if (within(ahead) == within(behind))
        continue;
      const gp_Pnt2d off      = within(ahead) ? behind : ahead;
      const double   tu       = surface->IsUPeriodic() ? surface->UPeriod() : 0,
                     tv       = surface->IsVPeriodic() ? surface->VPeriod() : 0;
      const gp_Vec2d shifts[] = {{tu, 0}, {-tu, 0}, {0, tv}, {0, -tv}};
      for (const auto& shift : shifts)
        if (shift.SquareMagnitude() > 0 && within(off.Translated(shift)))
        {
          seamShift = shift;
          side      = s2;
          at        = found->first;
          overAt    = found->second;
          over      = seam;
          overSeam  = true;
          entering  = within(ahead);
          break;
        }
      if (side < 0)
        fail(op,
             unsupported_geometry,
             "a blend's contact crosses its face's seam where the blend does not find the "
             "face's period beyond",
             {piece.rail.edge});
    }
  }
}

void Blend::splitAt(int k, Stations& stations, const Crossing& c)
{
  Piece&             piece     = pieces[k];
  const int          side      = c.side;
  const double       at        = c.at;
  const double       overAt    = c.overAt;
  const bool         entering  = c.entering;
  const bool         notched   = c.notched;
  const TopoDS_Edge& over      = c.over;
  const bool         overSeam  = c.overSeam;
  const gp_Vec2d&    seamShift = c.seamShift;
  // The next face across that edge, smooth with this one there.
  TopoDS_Face next = overSeam ? piece.rail.s[side].face : TopoDS_Face();
  for (const auto& face : facesOfEdge.FindFromKey(over))
    if (!face.IsSame(piece.rail.s[side].face))
      next = TopoDS::Face(face);
  const double t = piece.rail.t(at);
  Station      station;
  if (next.IsNull() || !stations.at(t, station))
    fail(op, unsupported_geometry, "a blend's contact runs off its face", {piece.rail.edge});
  station.t = t;
  // The edge's part before the crossing and its part after; the one off the face the
  // edge lies on rolls on the next face, from the station there moved onto it.
  Piece first = piece, second = piece;
  if (!first.rail.parameterSupports)
    first.rail.parameterSupports = std::array<Side, 2>{piece.rail.s[0], piece.rail.s[1]};
  second.rail.parameterSupports = first.rail.parameterSupports;
  first.rail.last = first.rail.to = t;
  second.rail.first = second.rail.from = t;
  Piece&  moved                        = entering ? first : second;
  Station onNext                       = station;
  if (notched)
  {
    moved.rail.s[side].trim       = next;
    moved.rail.s[side].trimEdge   = over;
    moved.rail.s[side].cliffWhole = false;
  }
  else
  {
    moved.rail.s[side].face    = next;
    moved.rail.s[side].surface = BRep_Tool::Surface(next);
    moved.rail.s[side].outward = next.Orientation() == TopAbs_REVERSED ? -1 : 1;
    GeomAPI_ProjectPointOnSurf onto(station.foot[side], moved.rail.s[side].surface);
    if (onto.NbPoints() == 0 || onto.LowerDistance() > piece.rail.resolved(limit::together))
      fail(op,
           unsupported_geometry,
           "a blend's contact runs off its face onto one it does not meet",
           {piece.rail.edge});
    onto.LowerDistanceParameters(onNext.uv[side][0], onNext.uv[side][1]);
    // Across its own seam, a period over (seamShift).
    if (overSeam)
    {
      onNext.uv[side][0] = station.uv[side][0] + seamShift.X();
      onNext.uv[side][1] = station.uv[side][1] + seamShift.Y();
    }
    else
    {
      // In the period of the next face's own parameters where the contact crosses onto
      // it: two faces on one periodic surface meet at the seam's place a period apart.
      double     f, l;
      const auto onOver = BRep_Tool::CurveOnSurface(over, next, f, l);
      if (!onOver.IsNull())
      {
        const gp_Pnt2d near = nearPeriod(gp_Pnt2d(onNext.uv[side][0], onNext.uv[side][1]),
                                         onOver->Value(overAt),
                                         moved.rail.s[side].surface);
        onNext.uv[side][0]  = near.X();
        onNext.uv[side][1]  = near.Y();
      }
    }
  }
  Jet before, after;
  if (!jetAt(piece.rail.s[side], station.uv[side][0], station.uv[side][1], before)
      || !jetAt(moved.rail.s[side], onNext.uv[side][0], onNext.uv[side][1], after)
      || !SectionPolicy::tangentFor(before.n.Angle(after.n), piece.rail.scale))
    fail(op,
         unsupported_geometry,
         "a blend's contact runs over an edge where its faces are not tangent, which the blend "
         "does not do yet",
         {piece.rail.edge});
  // `piece` refers into `pieces`, which must not grow past what it holds.
  if (pieces.size() >= pieces.capacity())
    fail(op,
         unsupported_geometry,
         "the blend crosses too many faces along one edge",
         {piece.rail.edge});
  const int  index = int(pieces.size());
  const auto outer = piece.joint[1];
  first.joint[1]   = {Meets::Split, index, 0, false, side, over};
  second.joint[0]  = {Meets::Split, k, 1, false, side, over};
  seeds[k]         = entering ? onNext : station;
  seeds.push_back(entering ? station : onNext);
  piece = first;
  // What met the edge's far end meets the second part's: on a closed edge, the first
  // part's own start.
  if (outer.meets == Meets::Joint || outer.meets == Meets::Miter || outer.meets == Meets::Split)
  {
    pieces[outer.piece].joint[outer.end].piece = index;
    pieces[outer.piece].joint[outer.end].end   = 1;
  }
  if (outer.meets == Meets::Corner)
    for (auto& [p2, e2] : corners[outer.piece].ends)
      if (p2 == k && e2 == 1)
        p2 = index;
  pieces.push_back(second);
}

bool Blend::sheetOf(int k, Folds& folds, Creases& creased)
{
  if (pieces.size() >= pieces.capacity())
    fail(op,
         unsupported_geometry,
         "the blend crosses too many faces along one edge",
         {pieces[k].rail.edge});
  Piece& piece = pieces[k];
  if (centred(piece.rail.edge))
    piece.rail.byCentre = true;
  // An open rail is walked along its spine only where its own stations fold there
  // (stationsFold): elsewhere they are the same balls, by the edge's parameter.
  if (request.intersectSpines && piece.rail.spine.IsNull() && intersectSpine(piece.rail)
      && piece.rail.spineArc && stationsFold(piece.rail) != true)
  {
    piece.rail.spine.Nullify();
    piece.rail.spineArc.reset();
  }
  else if (!request.intersectSpines && request.reads && piece.rail.spine.IsNull())
    request.reads->spines.push_back([rail = piece.rail]() mutable {
      return intersectSpine(rail) && !(rail.spineArc && stationsFold(rail) != true);
    });
  for (int attempt = 0; attempt < 2; ++attempt)
  {
    creased.resize(pieces.size());
    creased[k].clear();
    Stations stations(piece.rail);
    // No ball's plane passes through the edge's points somewhere along it (past a bend
    // tighter than the ball, where the contact on a cap folds): the balls are found through
    // their centres instead (Rail::byCentre), the same balls.
    try
    {
      stations.march(seeds[k]);
    }
    catch (const Failure&)
    {
      if (!stations.stalled() || !mayCentre(k))
        throw;
      centreChain(k);
      return false;
    }
    prepareRail(piece, stations);
    SheetFit fitting;
    if (!fitSheet(k, attempt, stations, fitting, creased[k]))
      return false;
    const bool onCliff = piece.rail.s[0].cliffCurve || piece.rail.s[1].cliffCurve
                         || !piece.rail.s[0].trim.IsNull() || !piece.rail.s[1].trim.IsNull();
    if (attempt == 1 || onCliff)
    {
      if (!fitting.contacts && !settleFolds(k, folds, onCliff))
        return false;
      break;
    }
    // The first place where a contact changes between inside its face and outside, and the edge
    // of that face it crosses there: leaving it, or entering it from the face it started on.
    const double u0 = piece.uAt(0), u1 = piece.uAt(1);
    const Next   off = offFace(k, stations, u0, u1);
    if (off == Next::Retry)
    {
      --attempt;
      continue;
    }
    if (off == Next::Stop)
      return false;
    Crossing crossed;
    crossingOverEdge(k, u0, u1, crossed);
    // Over a seam of its face the contact stays on the face, a period on, which the test
    // above does not see. Where it crosses one between the piece's ends (a rim whose
    // dihedral turns along it, so that its sections lean off the seam), the face's new
    // boundary has its vertex there on the seam: the piece is split there as where it
    // crosses onto the next face, the part after rolling on the same face a period over.
    if (crossed.side < 0 && !fitting.contacts)
      crossingOverSeam(k, u0, u1, crossed);
    if (crossed.side < 0)
    {
      if (fitting.contactPilot)
      {
        piece.fitted = fit(piece.rail, stations);
        piece.sheet  = fitting.exact ? *fitting.exact : Sheet{piece.fitted.blend};
      }
      // An exact sheet's contacts are lines or circles at constant speed: they cannot fold.
      if (!fitting.contacts && !settleFolds(k, folds, false))
        return false;
      break;
    }
    splitAt(k, stations, crossed);
    // The first part is looked at again, as the second will be: its other contact may cross
    // an edge too (the other face's seam just before a joint at a cone's, an edge of the
    // other support just after a joint of crossed extrusions).
    attempt = -1;
  }
  return true;
}

void Blend::creaseFolds(Folds& folds, Creases& creased)
{
  // Where a round crosses itself between its contacts, the edge's blend is cut where the sheets
  // before and after the fold meet (Attach::crossingFold): its parts end at the sections through
  // the crease's ends, the crease between them.
  creased.resize(pieces.size());
  for (int k = 0; k < int(pieces.size()); ++k)
  {
    if (creased[k].empty())
      continue;
    auto found =
      std::find_if(folds.begin(), folds.end(), [&](const auto& f) { return f.first == k; });
    if (found == folds.end())
      found = folds.insert(folds.end(), {k, {}});
    for (const auto& crease : creased[k])
    {
      const double a = crease->ends[0].X(), b = crease->ends[1].X();
      found->second.push_back(Fold{-1, std::min(a, b), std::max(a, b), crease});
    }
    std::sort(found->second.begin(), found->second.end(), [](const Fold& x, const Fold& y) {
      return x.from < y.from;
    });
    for (size_t i = 1; i < found->second.size(); ++i)
      if (found->second[i].from <= found->second[i - 1].to)
        refuseFold(k, "loops where a contact folds overlap");
  }
}

void Blend::cutAtFolds(const Folds& folds)
{
  // Where a contact on a flat cap folds over itself, the edge's blend is cut at the two ends of
  // each loop, where the balls are one, and the loop is patched (Attach::fold): a piece on the
  // one sheet between each loop and the next.
  for (const auto& [edge, found] : folds)
  {
    int        k = edge;
    const bool closed =
      pieces[edge].joint[1].meets == Meets::Joint && pieces[edge].joint[1].piece == edge;
    for (const Fold& fold : found)
    {
      if (pieces.size() >= pieces.capacity())
        fail(op,
             unsupported_geometry,
             "the blend folds too often along one edge",
             {pieces[k].rail.edge});
      Piece first = pieces[k], second = pieces[k];
      first.rail.last   = first.rail.t(fold.from);
      second.rail.first = second.rail.t(fold.to);
      // Its own copy of the sheet, so that an edge both parts meet on (a closed edge's seam) has
      // a curve on each face's surface, not two on one.
      second.fitted.blend  = Handle(Geom_BSplineSurface)::DownCast(first.fitted.blend->Copy());
      second.sheet.surface = second.fitted.blend;
      const int  index     = int(pieces.size());
      const auto outer     = pieces[k].joint[1];
      first.joint[1]       = {Meets::Fold, index, 0, false, fold.side, {}, fold.crease};
      second.joint[0]      = {Meets::Fold, k, 1, false, fold.side, {}, fold.crease};
      pieces[k]            = first;
      pieces.push_back(second);
      seeds.push_back(std::nullopt);
      // What met the edge's far end meets the second part's (the first's start, on a closed
      // edge).
      if (outer.meets == Meets::Joint || outer.meets == Meets::Miter || outer.meets == Meets::Split
          || outer.meets == Meets::Fold)
      {
        pieces[outer.piece].joint[outer.end].piece = index;
        pieces[outer.piece].joint[outer.end].end   = 1;
      }
      if (outer.meets == Meets::Corner)
        for (auto& [p2, e2] : corners[outer.piece].ends)
          if (p2 == k && e2 == 1)
            p2 = index;
      k = index;
    }
    // On a closed edge the last part runs on into the first across the edge's vertex, where the
    // contact on the face with a seam there crosses it: one blend face (joinSplits), as the
    // ball rolls on through the vertex.
    if (!closed)
      continue;
    TopoDS_Vertex vertex =
      TopExp::FirstVertex(TopoDS::Edge(pieces[edge].rail.edge.Oriented(TopAbs_FORWARD)));
    for (int side = 0; side < 2; ++side)
    {
      const TopoDS_Face& face = pieces[edge].rail.s[side].face;
      for (const auto& other : edgesOfVertex.FindFromKey(vertex))
      {
        if (other.IsSame(pieces[edge].rail.edge) || !BRep_Tool::IsClosed(TopoDS::Edge(other), face))
          continue;
        pieces[k].joint[1]    = {Meets::Split, edge, 0, false, side, TopoDS::Edge(other)};
        pieces[edge].joint[0] = {Meets::Split, k, 1, false, side, TopoDS::Edge(other)};
      }
    }
  }
}

void Blend::sheets()
{
  // Sections in two spans along a whole chain, or along none of it: its sections share their
  // form at its joints.
  shareHalves();
  // Each piece's blend, exact where every piece of its chain can be; a chain's sections then
  // share one parameter at each joint.
  // Where a contact leaves its face across a smooth edge onto the next face, the piece is split
  // there and its second part rolls on that face (the support sequence).
  seeds.assign(pieces.size(), std::nullopt);
  Folds folds;
  // Each piece's creases where its round crosses itself between its contacts, from its last fit.
  Creases creased;
  for (size_t k = 0; k < pieces.size(); ++k)
    if (!sheetOf(int(k), folds, creased))
      return;
  creaseFolds(folds, creased);
  crossFolds(folds);
  cutAtFolds(folds);
}

void Blend::crossFolds(const std::vector<std::pair<int, std::vector<Fold>>>& folds)
{
  // A loop that runs on across a joint of a chain (a closed outline split into edges, its
  // tightest bend at their vertex): the contact of the piece where it runs back crosses the next
  // piece's contact on the same cap. Both pieces are cut where they cross, the narrowest such
  // crossing round each place, and the loop between is patched as one within a piece is.
  struct Cross
  {
    int    a, ea, b, eb, side, sideB;
    double ua, ub;
  };

  std::vector<Cross> crosses;
  const auto         length = [](const Piece& p, double u0, double u1) {
    const double t0 = p.rail.t(u0), t1 = p.rail.t(u1);
    return railLength(p.rail, t0, t1);
  };
  // Whether a place of piece k lies in a crossing's loop.
  const auto inLoop = [](const Cross& c, int k, double at) {
    return (c.a == k && (c.ea == 1 ? at > c.ua : at < c.ua))
           || (c.b == k && (c.eb == 1 ? at > c.ub : at < c.ub));
  };
  for (const auto& [k, places] : open)
    for (const Fold& place : places)
    {
      if (std::any_of(crosses.begin(), crosses.end(), [&](const Cross& c) {
            return inLoop(c, k, place.from);
          }))
        continue;
      const Piece&         one = pieces[k];
      std::optional<Cross> best;
      double               shortest = HUGE_VAL;
      for (int e = 0; e < 2; ++e)
      {
        const auto& joint = one.joint[e];
        if (joint.meets != Meets::Joint || joint.piece == k)
          continue;
        const Piece& two = pieces[joint.piece];
        // The other piece's contact on the same cap: the one that ends where this one does.
        const gp_Pnt2d here  = one.fitted.contact[place.side]->Value(one.uAt(e));
        int            sideB = -1;
        for (int s2 = 0; s2 < 2; ++s2)
          if (!two.fitted.contact[s2].IsNull()
              && two.rail.s[s2].face.IsSame(one.rail.s[place.side].face)
              && (sideB < 0
                  || two.fitted.contact[s2]->Value(two.uAt(joint.end)).Distance(here)
                       < two.fitted.contact[sideB]->Value(two.uAt(joint.end)).Distance(here)))
            sideB = s2;
        if (sideB < 0)
          continue;
        const Geom2dAPI_InterCurveCurve crossing(
          new Geom2d_TrimmedCurve(one.fitted.contact[place.side], one.uAt(0), one.uAt(1)),
          new Geom2d_TrimmedCurve(two.fitted.contact[sideB], two.uAt(0), two.uAt(1)),
          tol::parameter);
        for (int i = 1; i <= crossing.NbPoints(); ++i)
        {
          const auto&  point = crossing.Intersector().Point(i);
          const double ua = point.ParamOnFirst(), ub = point.ParamOnSecond();
          // Not where the two contacts meet at the joint itself.
          const auto& cap = one.rail.s[place.side].surface;
          if (cap->Value(point.Value().X(), point.Value().Y())
                  .Distance(cap->Value(here.X(), here.Y()))
                <= one.rail.resolved(limit::together)
              || !(e == 1 ? ua < place.from : ua > place.from))
            continue;
          const double loop = length(one, ua, one.uAt(e)) + length(two, two.uAt(joint.end), ub);
          if (loop < shortest)
          {
            shortest = loop;
            best     = Cross{int(k), e, joint.piece, joint.end, place.side, sideB, ua, ub};
          }
        }
      }
      if (!best)
        refuseFold(k, "a contact on the cap folds where the blend finds no loop");
      crosses.push_back(*best);
    }
  // Each place where a contact runs back lies in one loop, and no loop overlaps one within a piece.
  for (const auto& [k, places] : open)
    for (const Fold& place : places)
      if (std::none_of(crosses.begin(), crosses.end(), [&](const Cross& c) {
            return inLoop(c, k, place.from);
          }))
        refuseFold(k, "a contact on the cap folds where the blend finds no loop");
  for (const auto& c : crosses)
    for (const auto& [k, found] : folds)
      for (const Fold& fold : found)
        if (inLoop(c, k, fold.from) || inLoop(c, k, fold.to))
          refuseFold(k, "loops where a contact folds overlap");
  for (const auto& c : crosses)
  {
    if (pieces[c.a].joint[c.ea].meets != Meets::Joint
        || pieces[c.b].joint[c.eb].meets != Meets::Joint)
      refuseFold(c.a, "loops where a contact folds overlap");
    for (const auto& [k, e, at] : {std::tuple{c.a, c.ea, c.ua}, std::tuple{c.b, c.eb, c.ub}})
    {
      Rail& rail                        = pieces[k].rail;
      (e == 1 ? rail.last : rail.first) = rail.t(at);
      if (!(rail.last > rail.first))
        refuseFold(k, "the loop where a contact folds takes a whole edge");
    }
    pieces[c.a].joint[c.ea] = {Meets::Fold, c.b, c.eb, c.side != c.sideB, c.side, {}};
    pieces[c.b].joint[c.eb] = {Meets::Fold, c.a, c.ea, c.side != c.sideB, c.sideB, {}};
  }
}

void Blend::chains()
{
  std::vector<int> chain(pieces.size(), -1);
  for (int k = 0; k < int(pieces.size()); ++k)
  {
    if (chain[k] >= 0)
      continue;
    std::vector<int> members, todo = {k};
    chain[k] = k;
    while (!todo.empty())
    {
      const int p = todo.back();
      todo.pop_back();
      members.push_back(p);
      for (const auto& joint : pieces[p].joint)
        if ((joint.meets == Meets::Joint || joint.meets == Meets::Split) && chain[joint.piece] < 0)
        {
          chain[joint.piece] = k;
          todo.push_back(joint.piece);
        }
    }
    const bool fitted =
      std::any_of(members.begin(), members.end(), [&](int p) { return !pieces[p].exact(); });
    if (fitted)
      for (int p : members)
      {
        // A piece whose sheet was exact is fitted now, to share its chain's parameter.
        if (pieces[p].fitted.blend.IsNull())
        {
          Stations stations(pieces[p].rail);
          stations.march(seeds[p]);
          pieces[p].fitted = fit(pieces[p].rail, stations);
        }
        pieces[p].sheet = Sheet{pieces[p].fitted.blend};
      }
  }
}

void Blend::joinSplits()
{
  joinSplitPieces(op, pieces, seeds);
}

void Blend::pairJoints()
{
  // Which sides pair at each joint: the contacts that end at the same place.
  for (auto& piece : pieces)
    for (int end = 0; end < 2; ++end)
    {
      auto& joint = piece.joint[end];
      if (joint.meets != Meets::Joint)
        continue;
      const Piece& other = pieces[joint.piece];
      const gp_Pnt here  = piece.sheet.contact(0)->Value(piece.sheet.rail(piece.uAt(end)));
      const auto   there = [&](int side) {
        return other.sheet.contact(side)->Value(other.sheet.rail(other.uAt(joint.end)));
      };
      joint.swapped = here.Distance(there(1)) < here.Distance(there(0));
    }
}

std::vector<TopoDS_Shape> Blend::foldedJoints() const
{
  std::vector<int> bad;
  for (int k = 0; k < int(pieces.size()); ++k)
    for (int end = 0; end < 2; ++end)
    {
      const auto& joint = pieces[k].joint[end];
      if (joint.meets != Meets::Joint)
        continue;
      const Piece& piece = pieces[k];
      const Piece& other = pieces[joint.piece];
      const auto   at    = [](const Piece& p, int side, int e) {
        return p.sheet.contact(side)->Value(p.sheet.rail(p.uAt(e)));
      };
      for (int side = 0; side < 2; ++side)
        if (at(piece, side, end).Distance(at(other, side ^ int(joint.swapped), joint.end))
            > piece.rail.resolved(limit::together))
          bad.push_back(k);
    }
  // Each such chain whole, through its joints and splits; only balls of a law the centres give.
  std::vector<TopoDS_Shape> edges;
  std::vector<bool>         seen(pieces.size(), false);
  while (!bad.empty())
  {
    const int k = bad.back();
    bad.pop_back();
    if (seen[k])
      continue;
    seen[k]          = true;
    const Rail& rail = pieces[k].rail;
    const bool  ball = rail.law == Law::Conic || rail.law == Law::G2;
    if (rail.byCentre || !ball || rail.s[0].cliffCurve || rail.s[1].cliffCurve)
      return {};
    if (std::none_of(edges.begin(), edges.end(), [&](const TopoDS_Shape& e) {
          return e.IsSame(rail.edge);
        }))
      edges.push_back(rail.edge);
    for (const auto& joint : pieces[k].joint)
      if (joint.meets == Meets::Joint || joint.meets == Meets::Split)
        bad.push_back(joint.piece);
  }
  return edges;
}

//! The full round of meeting rounds (Request::meet) split where it touches the face between: the
//! line of its sheet through the ball's foot there, between the face's boundary.
void Blend::splitMeeting(Blended& out)
{
  const Piece& piece = pieces.front();
  Stations     stations(piece.rail);
  stations.march();
  Station s;
  if (!stations.at((piece.rail.first + piece.rail.last) / 2, s))
    fail(op, no_solution, "the full round has no ball in its middle", {piece.rail.edge});
  const Sheet&               sheet = piece.sheet;
  GeomAPI_ProjectPointOnSurf onSheet(s.middle, sheet.surface);
  if (onSheet.NbPoints() == 0)
    fail(op, no_solution, "the full round does not touch the face between", {piece.rail.edge});
  double u, v;
  onSheet.LowerDistanceParameters(u, v);
  const double across = sheet.swapped ? u : v;
  TopoDS_Face  face;
  for (TopExp_Explorer f(out.shape, TopAbs_FACE); f.More() && face.IsNull(); f.Next())
    if (BRep_Tool::Surface(TopoDS::Face(f.Current())) == sheet.surface)
      face = TopoDS::Face(f.Current());
  if (face.IsNull())
    fail(op, no_solution, "the full round's face is missing", {piece.rail.edge});
  // Where the line meets the face's boundary, along the rail.
  const Handle(Geom2d_Line) line =
    new Geom2d_Line(sheet.uv(0, across), sheet.swapped ? gp_Dir2d(0, 1) : gp_Dir2d(1, 0));
  std::vector<double> ends;
  for (TopExp_Explorer e(face, TopAbs_EDGE); e.More(); e.Next())
  {
    double     f, l;
    const auto pcurve = BRep_Tool::CurveOnSurface(TopoDS::Edge(e.Current()), face, f, l);
    if (pcurve.IsNull())
      continue;
    Geom2dAPI_InterCurveCurve meet(new Geom2d_TrimmedCurve(line, -1e6, 1e6),
                                   new Geom2d_TrimmedCurve(pcurve, f, l),
                                   tol::parameter);
    for (int i = 1; i <= meet.NbPoints(); ++i)
      ends.push_back(ElCLib::Parameter(line->Lin2d(), meet.Point(i)));
  }
  if (ends.size() < 2)
    fail(op,
         no_solution,
         "the line where the full round touches the face between does not cross it",
         {piece.rail.edge});
  const auto [lo, hi] = std::minmax_element(ends.begin(), ends.end());
  const auto  curve   = sheet.swapped ? sheet.surface->UIso(across) : sheet.surface->VIso(across);
  TopoDS_Edge edge    = BRepBuilderAPI_MakeEdge(curve, *lo, *hi);
  BRep_Builder().UpdateEdge(edge, line, face, tol::length);
  BRepAlgoAPI_Splitter           split;
  NCollection_List<TopoDS_Shape> arguments, tools;
  arguments.Append(out.shape);
  tools.Append(edge);
  split.SetArguments(arguments);
  split.SetTools(tools);
  split.Build();
  if (!split.IsDone())
    fail(op,
         no_solution,
         "the full round cannot be split where it touches the face between",
         {piece.rail.edge});
  TopoDS_Shape solid;
  for (TopExp_Explorer x(split.Shape(), TopAbs_SOLID); x.More(); x.Next())
    solid = x.Current();
  if (solid.IsNull())
    fail(op, no_solution, "splitting the full round gave no solid", {piece.rail.edge});
  out.joins = carriedJoins(out.joins, split.History(), solid);
  for (const auto& made : split.Modified(edge))
    if (made.ShapeType() == TopAbs_EDGE)
      out.joins.push_back({made, JoinRole::UsedSupport, 1});
  out.history->Merge(split.History());
  out.shape = solid;
}

Blended Blend::result()
{
  Blended                         out;
  const Handle(BRepTools_History) made = new BRepTools_History;
  out.shape = pieces.empty() ? body
                             : attach(body, pieces, corners, made, &out.joins,
                                      request.refineContacts, request.reads);
  // The input's history: into the copy, then through the blends.
  NCollection_List<TopoDS_Shape> inputs;
  inputs.Append(request.body);
  out.history = new BRepTools_History(inputs, copy);
  out.history->Merge(made);
  if (request.meet > 0)
    splitMeeting(out);
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> result;
  TopExp::MapShapes(out.shape, TopAbs_FACE, result);
  TopExp::MapShapes(out.shape, TopAbs_EDGE, result);
  for (const auto type : {TopAbs_FACE, TopAbs_EDGE})
  {
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> parts;
    TopExp::MapShapes(request.body, type, parts);
    for (int i = 1; i <= parts.Extent(); ++i)
    {
      const TopoDS_Shape inCopy = copy.ModifiedShape(parts(i));
      if (!result.Contains(inCopy))
        continue;
      out.unchanged.emplace_back(parts(i), inCopy);
    }
  }
  return out;
}

Blended blendCopy(const Request&                   request,
                  BRepBuilderAPI_Copy&             copy,
                  const std::vector<TopoDS_Shape>& byCentre)
{
  Blend b(request, copy, byCentre);
  // Each stage a phase of the kernel's profile, where one is asked for (diagnostics).
  const auto stage = [&](const char* name, auto&& run) {
    checkCancelled(request.operation);
    BooleanPhase phase(name);
    run();
  };
  stage("fillet2.rails", [&] { b.rails(); });
  stage("fillet2.joints", [&] { b.joints(); });
  stage("fillet2.carries", [&] { b.carries(); });
  stage("fillet2.sheets", [&] {
    // Made again from the pieces as they were, with each chain centreChain marked, until
    // none is marked anew.
    const std::vector<Piece>  pieces  = b.pieces;
    const std::vector<Corner> corners = b.corners;
    // A fold set aside for a loop across a joint came first along the chain: where a later part
    // fails, the fold is what is refused.
    try
    {
      for (b.sheets(); b.again; b.sheets())
      {
        b.again = false;
        b.open.clear();
        b.pieces.assign(pieces.begin(), pieces.end());
        b.corners = corners;
        b.absorbEnds();
      }
    }
    catch (const Failure& failure)
    {
      if (failure.reason == std::string("cancelled") || b.open.empty())
        throw;
      b.refuseFold(b.open.front().first, failure.message.c_str());
    }
  });
  stage("fillet2.chains", [&] { b.chains(); });
  stage("fillet2.joinSplits", [&] { b.joinSplits(); });
  stage("fillet2.pairJoints", [&] { b.pairJoints(); });
  // Two blends of a chain that end in different balls at their joint: the balls whose planes
  // pass through the edge's points fold back along it (an acute groove on a tight bend), and the
  // chain's balls are found through their centres instead, the same balls (Rail::byCentre).
  if (byCentre.empty())
    if (auto folded = b.foldedJoints(); !folded.empty())
    {
      for (const auto& edge : b.byCentre)
        if (std::none_of(folded.begin(), folded.end(), [&](const TopoDS_Shape& e) {
              return e.IsSame(edge);
            }))
          folded.push_back(edge);
      return blendCopy(request, copy, folded);
    }
  BooleanPhase phase("fillet2.attach");
  return b.result();
}

} // namespace BRepExt::blending
