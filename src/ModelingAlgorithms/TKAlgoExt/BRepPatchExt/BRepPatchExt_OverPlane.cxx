// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Any loop: a surface over the loop's mean plane, trimmed by the loop. Its net is finer where the
// rim turns fast, and may change its bend across lines through the corners where the faces along
// the rim do.

#include <BRepPatchExt_Fill.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>

#include <GeomAlgoExt_Cancel.hxx>
#include <GeomAlgoExt_PatchRim.hxx>

#include <BRepAdaptor_Surface.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BSplCLib.hxx>
#include <GeomConvert.hxx>
#include <GeomConvert_CompCurveToBSplineCurve.hxx>
#include <Geom_BezierCurve.hxx>
#include <ShapeAnalysis_Surface.hxx>

#include <limits>
#include <tuple>

namespace BRepExt::patch_fill
{
namespace
{

//! How much less the surface moved onto its curves must bend than the one fitted at the
//! tolerance to be kept instead.
constexpr double clearlyFairer = 0.75;

//! How far the rim may turn within one span of the coarsest net, in radians, before the net
//! is made finer there: its own direction, or the tangent plane of the face it runs on.
constexpr double turnPerSpan = 0.2;

//! The places along a side at which the turning of its face is read.
constexpr int turningSteps = 48;

//! The nets that are tried with the neighbours' bend wished for, before the surface is made
//! tangent and no more.
constexpr int wishNets = 3;

//! The steps along a curve at which it is looked for the lines of the net that it crosses.
constexpr int crossingSteps = 256;

//! Two lines of the net across which it may change its bend count as well apart when this many
//! of them, as far apart, would cross the loop: about the spans of the third net.
constexpr double apartSpans = 16;

//! A miss of no more than this many times its limit is nearly met.
constexpr double nearly = 2;

//! What a step of the second derivative across a line of lower continuity counts for against
//! the bend, per length of the line and with the loop's size as the unit; the places to a span
//! along the line at which it is read; and how far from its corner, as a share of the loop's
//! size, a step counts in full.
constexpr double jumpWeight = 30, jumpReach = 0.3;

constexpr int jumpPlaces = 3;

//! A planar section guide carried along a rim, with support-derived Hermite extensions.
//! This is the translation member of the curve-network Boolean-sum construction.
//! It is only a candidate: all original rims, jets, guides and interior checks still apply.
std::optional<Built> sectionReference(const Asked& asked, const Frame& frame)
{
  if (asked.guides.size() != 1 || !asked.points.empty() || asked.sides.size() != 4)
    return {};
  std::vector<const Side*> rims;
  for (const Side& side : asked.sides)
    if (!side.free)
      rims.push_back(&side);
  if (rims.size() != 2)
    return {};
  if (rims[1]->order > rims[0]->order)
    std::swap(rims[0], rims[1]);
  const Side &rim = *rims[0], &foot = *rims[1], &guide = asked.guides.front();
  auto        profile = plain(guide, guide.length, false, asked.limits.g0);
  auto        rail    = plain(rim, rim.length, false, asked.limits.g0);
  if (profile.IsNull() || rail.IsNull())
    return {};
  std::vector<gp_Pnt> samples;
  for (int i = 1; i <= profile->NbPoles(); ++i)
    samples.push_back(profile->Pole(i));
  const Frame section = meanPlane(samples, guide.length);
  // The section must be planar; reject rather than moving any supplied guide.
  for (int i = 1; i <= profile->NbPoles(); ++i)
    if (std::abs(gp_Vec(section.origin, profile->Pole(i)).Dot(section.z)) > asked.limits.g0)
      return {};
  const auto intersect = [&](const Side& side) -> std::optional<double> {
    const auto height = [&](double f) { return gp_Vec(section.origin, side.at(f)).Dot(section.z); };
    double     a = 0, b = 1, ha = height(a), hb = height(b);
    if (ha * hb >= 0)
      return {};
    for (int i = 0; i < std::numeric_limits<double>::digits; ++i)
    {
      const double m = (a + b) / 2, h = height(m);
      if ((h >= 0) == (ha >= 0))
      {
        a  = m;
        ha = h;
      }
      else
        b = m;
    }
    return (a + b) / 2;
  };
  const auto ra = intersect(rim), rb = intersect(foot);
  if (!ra || !rb)
    return {};
  const gp_Pnt a = rim.meet(*ra), b = foot.meet(*rb);
  if (a.Distance(profile->Value(profile->FirstParameter()))
        + b.Distance(profile->Value(profile->LastParameter()))
      > a.Distance(profile->Value(profile->LastParameter()))
          + b.Distance(profile->Value(profile->FirstParameter())))
    profile->Reverse();
  gp_Pnt g0, g1;
  gp_Vec d0, d1, e0, e1;
  profile->D2(profile->FirstParameter(), g0, d0, e0);
  profile->D2(profile->LastParameter(), g1, d1, e1);
  const double h0 = a.Distance(g0), h1 = b.Distance(g1);
  if (h0 <= asked.limits.g0 || h1 <= asked.limits.g0)
    return {};
  const auto jet = [&](const Side&   side,
                       double        f,
                       const gp_Pnt& toward,
                       const gp_Vec& fallback) -> std::pair<gp_Vec, gp_Vec> {
    if (side.support.IsNull() || side.order == 0)
      return {fallback, gp_Vec(0, 0, 0)};
    const SurfaceJet j = side.jet(f);
    gp_Vec           t = j.n.Crossed(section.z);
    t.Normalize();
    if (t.Dot(gp_Vec(side.meet(f), toward)) < 0)
      t.Reverse();
    const double E = j.u.SquareMagnitude(), F = j.u.Dot(j.v), G = j.v.SquareMagnitude();
    const double u = (G * t.Dot(j.u) - F * t.Dot(j.v)) / (E * G - F * F),
                 v = (E * t.Dot(j.v) - F * t.Dot(j.u)) / (E * G - F * F);
    gp_Vec e       = j.uu * (u * u) + j.uv * (2 * u * v) + j.vv * (v * v);
    // Differentiate the plane-intersection constraint, preserving the normal acceleration.
    const gp_Vec w = section.z - j.n * section.z.Dot(j.n);
    e -= w * (e.Dot(section.z) / w.SquareMagnitude());
    return {t, e};
  };
  const auto [da, ea] = jet(rim, *ra, g0, d0);
  auto [db, eb]       = jet(foot, *rb, g1, -d1);
  if (foot.order > 0)
    db.Reverse(); // The section arrives at the second rim.
  else
    db = d1;
  const auto hermite = [](const gp_Pnt& p,
                          const gp_Vec& d,
                          const gp_Vec& e,
                          const gp_Pnt& q,
                          const gp_Vec& f,
                          const gp_Vec& g,
                          double        h) {
    NCollection_Array1<gp_Pnt> poles(1, 6);
    poles(1)   = p;
    poles(2)   = p.Translated(d * (h / 5));
    poles(3)   = p.Translated(d * (2 * h / 5) + e * (h * h / 20));
    poles(4)   = q.Translated(f * (-2 * h / 5) + g * (h * h / 20));
    poles(5)   = q.Translated(f * (-h / 5));
    poles(6)   = q;
    auto curve = GeomConvert::CurveToBSplineCurve(new Geom_BezierCurve(poles));
    NCollection_Array1<double> knots = curve->Knots();
    BSplCLib::Reparametrize(0, h, knots);
    curve->SetKnots(knots);
    return curve;
  };
  profile->IncreaseDegree(5);
  GeomConvert_CompCurveToBSplineCurve joined(hermite(a, da, ea, g0, d0, e0, h0));
  if (!joined.Add(profile, asked.limits.g0 * 1e-6, true, false, 3)
      || !joined.Add(hermite(g1, d1, e1, b, db, eb, h1), asked.limits.g0 * 1e-6, true, false, 3))
    return {};
  const auto                 cross = joined.BSplineCurve();
  NCollection_Array2<gp_Pnt> poles(1, rail->NbPoles(), 1, cross->NbPoles());
  for (int i = 1; i <= rail->NbPoles(); ++i)
    for (int j = 1; j <= cross->NbPoles(); ++j)
      poles(i, j) = gp_Pnt(rail->Pole(i).XYZ() + cross->Pole(j).XYZ() - a.XYZ());
  const occ::handle<Geom_BSplineSurface> surface = new Geom_BSplineSurface(poles,
                                                                           rail->Knots(),
                                                                           cross->Knots(),
                                                                           rail->Multiplicities(),
                                                                           cross->Multiplicities(),
                                                                           rail->Degree(),
                                                                           cross->Degree());
  ShapeAnalysis_Surface                  project(surface);
  const auto                             layout = [&](const gp_Pnt& p) {
    return project.ValueOfUV(p, asked.limits.g0 * 1e-3).XY();
  };
  std::vector<Rim> traces(asked.sides.size());
  for (size_t i = 0; i < asked.sides.size(); ++i)
    if (!asked.sides[i].free)
      traces[i] = preimage(surface, asked.sides[i], frame, asked.limits.g0, asked.size, layout);
  patch_side::closeGaps(asked.sides, traces);
  std::vector<gp_XY> loop;
  for (size_t i = 0; i < asked.sides.size(); ++i)
    loop.push_back(traces[i].onSurface->Value(asked.sides[i].from).XY());
  Built built{
    3,
    {},
    measureAll(surface,
               asked,
               [&](const Side& side, double f, const gp_Pnt&) -> std::optional<Place> {
                 const auto i = asked.boundaryIndex(side);
                 return i ? std::optional<Place>(
                              Place{traces[*i].onSurface->Value(side.parameter(f)).XY(), true})
                          : std::nullopt;
               })};
  if (!built.measured.within)
    return {};
  built.face     = trimmed(surface, asked, traces, enclosed(loop) < 0);
  built.crossing = Crossing{true, rim.order, foot.order};
  for (size_t i = 0; i < asked.sides.size(); ++i)
  {
    const Side& side = asked.sides[i];
    if (!side.free)
      built.spatialBoundary.push_back(
        {side.input, traces[i].onSurface, side.from, side.to, side.originalReversed});
    built.gaps.push_back(side.free ? traces[i].onSurface : occ::handle<Geom2d_Curve>());
  }
  measureSpatialBoundary(built, asked, false);
  if (!built.measured.within
      || !regularMargin(surface, [](double, double, double, double) { return true; })
      || !BRepCheck_Analyzer(built.face).IsValid()
      || insideFault(built.face, asked.points.empty() ? rimBend(asked) : 0, built.crossing))
    return {};
  return built;
}

} // namespace

Built overPlane(const Asked& asked, const Frame& flat, int direction)
{
  // Corners where the faces on both sides are tangent and bend differently (a rim that runs
  // from a plane onto a round): there the surface has to change its bend, and a net that
  // is twice differentiable everywhere can only smear that over a span.
  std::vector<gp_Pnt> turns;
  std::vector<gp_Vec> along;
  for (size_t i = 0; i < asked.sides.size(); ++i)
  {
    const Side &a = asked.sides[i], &b = asked.sides[(i + 1) % asked.sides.size()];
    if (a.order < 1 || b.order < 1 || a.support.IsSame(b.support)
        || cornerConflict(a, b) <= asked.limits.g2)
      continue;
    turns.push_back(a.at(1));
    along.push_back(a.tangent(1));
    along.push_back(b.tangent(0));
  }
  // The net's lines through such a corner carry the change best when the rim crosses both
  // of them well: the net is turned in its plane to where the rim that is nearest to one of
  // its directions is furthest from it.
  // And two such lines that run close beside each other leave the surface next to nothing
  // to hold it between them: the corners are to lie well apart both ways.
  Frame frame = flat;
  if (!turns.empty())
  {
    const std::vector<gp_Pnt> outline = polygon(asked.sides);
    double                    best = -1, angle = 0;
    for (int step = 0; step < 90; ++step)
    {
      const double turn  = step * M_PI / 180;
      Frame        tried = flat;
      tried.x            = flat.x * std::cos(turn) + flat.y * std::sin(turn);
      tried.y            = flat.y * std::cos(turn) - flat.x * std::sin(turn);
      double worst       = 1;
      for (const gp_Vec& way : along)
        worst =
          std::min(worst, std::abs(std::sin(2 * std::atan2(way.Dot(tried.y), way.Dot(tried.x)))));
      gp_XY low = tried.project(outline.front()), high = low;
      for (const gp_Pnt& point : outline)
      {
        const gp_XY at = tried.project(point);
        low.SetCoord(std::min(low.X(), at.X()), std::min(low.Y(), at.Y()));
        high.SetCoord(std::max(high.X(), at.X()), std::max(high.Y(), at.Y()));
      }
      for (size_t i = 0; i < turns.size(); ++i)
        for (size_t k = i + 1; k < turns.size(); ++k)
        {
          const gp_XY apart = tried.project(turns[i]) - tried.project(turns[k]);
          worst             = std::min({worst,
                                        std::abs(apart.X()) / (high.X() - low.X()) * apartSpans,
                                        std::abs(apart.Y()) / (high.Y() - low.Y()) * apartSpans});
        }
      if (worst > best + 1e-9)
      {
        best  = worst;
        angle = turn;
      }
    }
    angle += direction * M_PI / 9;
    frame.x = flat.x * std::cos(angle) + flat.y * std::sin(angle);
    frame.y = flat.y * std::cos(angle) - flat.x * std::sin(angle);
  }
  std::vector<double> breaksU, breaksV;
  for (const gp_Pnt& turn : turns)
  {
    breaksU.push_back(frame.project(turn).X());
    breaksV.push_back(frame.project(turn).Y());
  }
  std::vector<gp_XY> loop;
  for (const gp_Pnt& point : polygon(asked.sides))
    loop.push_back(frame.project(point));
  // A loop that is refused over the mean plane, and whose rim's faces have normals that no plane
  // has all within `steepest` of its own, is a height over no plane at all: it is filled over a
  // chart of its own instead (see overChart). Every other loop is filled as it was, and one that
  // the chart does not fill either is refused as it was.
  const auto refused = [&](const std::string& refusal) -> Built {
    // The plane's check reads the faces between the ends; whether any plane would do is read
    // at the ends as well, where a rim's normals close their sweep (the nose of a wing).
    std::vector<gp_Vec> normals;
    for (const Side& side : asked.sides)
      for (int i = 0; direction == 0 && side.order > 0 && i <= gaussPlaces; ++i)
        normals.push_back(
          side.jet(std::clamp(static_cast<double>(i) / gaussPlaces, cornerInset, 1 - cornerInset))
            .n);
    if (direction == 0 && graphMargin(normals, {frame.x, frame.y, frame.z}) < steepest)
    {
      std::optional<Built> charted;
      try
      {
        charted = overChart(asked, flat);
      }
      catch (const Failure& failure)
      {
        if (failure.reason == "cancelled")
          throw;
      }
      catch (const Standard_Failure&)
      {
      }
      if (charted)
        return *charted;
    }
    fail(op, unsupported_geometry, refusal);
  };
  double                     area = 0;
  std::optional<std::string> crossed;
  try
  {
    area = enclosed(loop);
  }
  catch (const Failure& failure)
  {
    crossed = failure.message;
  }
  if (crossed)
    return refused(*crossed);
  const Inside within(loop, 0);
  // A fixed graph avoids spending fairness on tangential parameter drift. Keep the
  // spatial fit for steep supports, where its movable parameters approximate shape
  // energy better. The graph also makes the support's G2 conditions exactly linear.
  constexpr double graphNormal   = 0.8;
  bool             graphSupports = true;
  for (const Side& side : asked.sides)
    for (int i = 0; side.order > 0 && i < 16; ++i)
    {
      const double normal = std::abs(side.jet((i + 0.5) / 16).n.Dot(frame.z));
      graphSupports       = graphSupports && normal >= graphNormal;
      if (normal < steepest)
        return refused(
          "a face to be smooth with stands too steeply against the boundary's mean plane");
    }
  for (const Side& guide : asked.guides)
    for (int i = 0; i <= 16; ++i)
      if (!inside(loop, frame.project(guide.at(i / 16.0)), outsideSlack * asked.size))
        fail(op, invalid_input, "a guide leaves the boundary");
  for (const gp_Pnt& point : asked.points)
    if (!inside(loop, frame.project(point), outsideSlack * asked.size))
      fail(op, invalid_input, "a point lies outside the boundary");

  if (!asked.guides.empty())
  {
    Asked unguided = asked;
    unguided.guides.clear();
    const bool           gapped    = hasGaps(asked);
    std::optional<Built> candidate = gapped ? std::nullopt : spherical(unguided, flat);
    if (!candidate && !gapped)
      candidate = ruledRound(unguided, flat);
    if (!candidate && !gapped && unguided.points.empty()
        && squarish(unguided.sides, asked.limits.g1))
      candidate = foursided(unguided, flat);
    Built base = candidate && candidate->measured.within ? *candidate : overPlane(unguided, flat);
    if (base.measured.within)
    {
      const auto          surface = BRep_Tool::Surface(base.face);
      std::vector<double> misses;
      for (const Side& guide : asked.guides)
        misses.push_back(measure(surface, guide, [&](const Side&, double, const gp_Pnt& point) {
                           return std::optional<Place>(Place{frame.project(point), false});
                         }).g0);
      if (std::all_of(misses.begin(), misses.end(), [&](double miss) {
            return miss <= asked.limits.g0;
          }))
      {
        base.measured.values.insert(base.measured.values.begin() + asked.sides.size() * 4,
                                    misses.begin(),
                                    misses.end());
        return base;
      }
    }
  }

  gp_XY low = loop.front(), high = loop.front();
  for (const gp_XY& point : loop)
  {
    low.SetCoord(std::min(low.X(), point.X()), std::min(low.Y(), point.Y()));
    high.SetCoord(std::max(high.X(), point.X()), std::max(high.Y(), point.Y()));
  }
  const double width = high.X() - low.X(), depth = high.Y() - low.Y(),
               larger = std::max(width, depth);
  // Each way by its own extent: a tenth of the longer one on every side would leave a thin
  // strip a small part of its own surface.
  const gp_XY margin(marginShare * width, marginShare * depth);
  int         order = 0;
  for (const Side& side : asked.sides)
    order = std::max(order, side.order);
  // Degree 5 for tangency and guides; cubic for free boundaries without guides.
  // A net of degree 3 follows a rim's tangent planes as a cubic follows a curve, and the
  // corner of three rounds, a ball's eighth, took it 25 x 25 poles and came out rippled (a
  // largest bend of 0.56 for 0.25); of degree 5 it takes 13 x 13 and bends 0.31 at most.
  int        degree = order > 0 || !asked.guides.empty() ? 5 : 3;
  const auto typeOf = [](const Side& side) { return BRepAdaptor_Surface(side.support).GetType(); };
  const bool supported = std::all_of(asked.sides.begin(), asked.sides.end(), [](const Side& side) {
    return !side.support.IsNull();
  });
  const bool exchange =
    order > 0 && asked.guides.empty() && supported
    && ((asked.sides.size() == 1 && typeOf(asked.sides.front()) == GeomAbs_Cone)
        || (std::any_of(asked.sides.begin(),
                        asked.sides.end(),
                        [&](const Side& side) { return typeOf(side) == GeomAbs_Plane; })
            && std::any_of(asked.sides.begin(), asked.sides.end(), [&](const Side& side) {
                 return typeOf(side) != GeomAbs_Plane;
               })));

  // The fairest surface tangent to a ball is not the ball: of all that leave the rim in the
  // neighbours' planes, the one that bends least is flatter. Where the neighbours agree
  // about the bend all round the rim, their bend is wished for as well, and what comes of
  // that is taken where it is within the limits and one to keep; where they do not agree,
  // or it is not, the surface is tangent and no more.
  const double bend    = rimBend(asked);
  const bool   fairRim = order == 0 && asked.guides.empty() && asked.points.empty();
  // Where curves are only met, the fit is asked to come as near them as a fair surface over
  // its net does, a share of the loop's size and nothing of the tolerance, and the surface is
  // moved onto them afterwards: brought nearer by the fit, a coarse net ripples (z = sin 2a
  // on a radius of 6 bends 0.44 for 0.056 at 2e-5 and 0.06 at 5e-5). Where a face is gone
  // over into, the fit is asked for the tolerance itself: its net of degree 5 follows a rim
  // on a face to a millionth with few poles, and moved there afterwards took three times as
  // many.
  bool nearlyQuadratic = false;
  if (order == 0 && (!asked.guides.empty() || !asked.points.empty()))
  {
    // Recognition only: a quadratic needs neither the general guide fit's degree five
    // nor its extra poles. The real net is seeded below in the recognised degree.
    SplineNet probe{2,
                    2,
                    clampedKnots(low.X() - margin.X(), high.X() + margin.X(), 1, 2),
                    clampedKnots(low.Y() - margin.Y(), high.Y() + margin.Y(), 1, 2),
                    std::vector<gp_XYZ>(9)};
    if (polynomialHeight(asked, frame, probe))
      degree = 2;
    else if (!asked.guides.empty())
      nearlyQuadratic = polynomialHeight(asked, frame, probe, true);
  }
  Limits reach            = asked.limits;
  double inputUncertainty = 0;
  for (const Side& side : asked.sides)
    for (int k = 0; k < places; ++k)
    {
      const double f   = (k + .5) / places;
      inputUncertainty = std::max(inputUncertainty, side.at(f).Distance(side.meet(f)));
    }
  // Leave extra angle headroom when sewing must reparameterize a coarse input edge.
  if (inputUncertainty > 10 * asked.limits.g0)
    reach.g1 *= .5;
  GraphFit   graphMode = GraphFit::transverse;
  const auto build     = [&](bool wish, int nets) {
    std::optional<Built>      best;
    int                       tried = 0, unkept = 0;
    double                    placeOnly     = HUGE_VAL;
    int                       previousSpans = 0;
    occ::handle<Geom_Surface> previous;
    bool                      previousWithin = false, previousSpatial = false, settledOn = false;
    std::vector<Built>        meeting;
    // The corners of two sides that are only to be met at which a surface of the ladder left the
    // sides nearer the corner than the even places of its correction reach (see below).
    std::vector<char> forced(asked.sides.size(), 0);
    for (int rung : ladder)
    {
      const int spans = std::max(rung,
                                 asked.controls.quality == PatchQuality::Max    ? 16
                                     : asked.controls.quality == PatchQuality::High ? 10
                                                                                    : 0);
      if (spans == previousSpans)
        continue;
      previousSpans = spans;
      // The last is for nets that only meet curves: a height of degree 3, which that many
      // spans cost no more than the net before it costs one that is free in space.
      if (nets-- <= 0 || (spans == std::end(ladder)[-1] && order > 0))
        break;
      auto along = [&](double extent) {
        // The complete-G2 retry refines both directions, including a narrow
        // direction whose aspect-ratio share never exceeds the coarse minimum.
        // Preserve established G0/G1 and transverse constructions.
        const int minimum =
          graphMode == GraphFit::complete
                ? static_cast<int>(std::lround(static_cast<double>(fewestSpans) * spans / ladder[0]))
                : fewestSpans;
        return std::max(minimum, static_cast<int>(std::lround(spans * extent / larger)));
      };
      const int spansU = along(width), spansV = along(depth);
      // Finer where the rim turns faster than the even net follows: by `turnPerSpan` of
      // turning to a span on the coarsest net, and in step with the net after that.
      std::vector<Finer> finerU, finerV;
      for (const Side& side : asked.sides)
        for (int i = 0; side.order > 0 && i < turningSteps; ++i)
        {
          const double f = (i + 0.5) / turningSteps, step = 0.25 / turningSteps;
          const double turning =
            side.jet(f - step).n.Angle(side.jet(f + step).n) / (2 * step * side.length);
          const double span = turnPerSpan * ladder[0] / spans / std::max(turning, 1e-12);
          if (span >= larger / spans)
            continue;
          const gp_XY at = frame.project(side.at(f));
          finerU.push_back({at.X(), span});
          finerV.push_back({at.Y(), span});
        }
      // A narrow opening needs longitudinal resolution near its ends too:
      // the G2 boundary layer can be as short as the opening is wide, even
      // when the boundary tangent hardly turns. Refine only the failed-G2
      // retry, and only where its even spans exceed that physical scale.
      if (graphMode == GraphFit::complete)
      {
        const double span = std::min(width, depth) / fewestSpans;
        const auto   ends = [&](double from, double to, int count, std::vector<Finer>& finer) {
          if ((to - from) / count <= span)
            return;
          for (double end : {from, to})
            for (int i = -fewestSpans; i <= fewestSpans; ++i)
              finer.push_back({end + i * span, span});
        };
        ends(low.X(), high.X(), spansU, finerU);
        ends(low.Y(), high.Y(), spansV, finerV);
      }
      SplineNet net{
        degree,
        degree,
        gradedKnots(low.X() - margin.X(), high.X() + margin.X(), spansU, degree, finerU, breaksU),
        gradedKnots(low.Y() - margin.Y(), high.Y() + margin.Y(), spansV, degree, finerV, breaksV),
            {}};
      // The net lies in the mean plane, each pole where the plane's own parameters put it,
      // and is free to rise from it: a height over the plane.
      const int  nu = net.countU(), nv = net.countV();
      const auto mean = [](const std::vector<double>& knots, int degree, int i) {
        double sum = 0;
        for (int k = 1; k <= degree; ++k)
          sum += knots[i + k];
        return sum / degree;
      };
      net.poles.resize(static_cast<size_t>(nu) * nv);
      for (int i = 0; i < nu; ++i)
        for (int j = 0; j < nv; ++j)
          net.poles[static_cast<size_t>(i) * nv + j] =
            (frame.origin.XYZ() + frame.x.XYZ() * mean(net.knotsU, degree, i)
             + frame.y.XYZ() * mean(net.knotsV, degree, j));
      const gp_XYZ up = frame.z.XYZ();
      const bool   seeded =
        order == 0 && (reach.g0 > asked.limits.g0 || !asked.guides.empty() || !asked.points.empty())
        && polynomialHeight(asked, frame, net, reach.g0 > asked.limits.g0);
      std::vector<Target> targets;
      for (const auto* curves : {&asked.sides, &asked.guides})
        for (const Side& side : *curves)
        {
          if (side.free)
            continue;
          // The same number of places in every span of the net that the curve passes
          // through, however briefly: a span it only clips has poles of its own, and
          // without a place there nothing but the energy would say where they go.
          std::vector<double> cuts{0, 1};
          gp_XY               before = frame.project(side.at(0));
          for (int i = 1; i <= crossingSteps; ++i)
          {
            const gp_XY here = frame.project(side.at(static_cast<double>(i) / crossingSteps));
            for (int way = 1; way <= 2; ++way)
            {
              const std::vector<double>& knots = way == 1 ? net.knotsU : net.knotsV;
              const double               a = before.Coord(way), b = here.Coord(way);
              for (auto knot = std::upper_bound(knots.begin(), knots.end(), std::min(a, b));
                   knot != knots.end() && *knot < std::max(a, b);
                   ++knot)
                cuts.push_back((i - 1 + (*knot - a) / (b - a)) / crossingSteps);
            }
            before = here;
          }
          std::sort(cuts.begin(), cuts.end());
          // Next to a gap, or where the next edge lies on the same face, nothing else is
          // asked of that end: the face has its say there as along the rest.
          const auto index = asked.boundaryIndex(side);
          const auto quiet = [&](size_t other) {
            const Side& near = asked.sides[other % asked.sides.size()];
            return near.free || near.support.IsSame(side.support);
          };
          const bool lone       = index && hasGaps(asked);
          const bool quietStart = lone && quiet(*index + asked.sides.size() - 1),
                     quietEnd   = lone && quiet(*index + 1);
          auto place            = [&](double f, double share, bool end, bool watch = false) {
            if (end && ((f == 0 && quietStart) || (f == 1 && quietEnd)))
            {
              f   = std::clamp(f, cornerInset, 1 - cornerInset);
              end = false;
            }
            const gp_Pnt point = side.meet(f);
            // Tangency and bend are not asked for at the ends, where two sides meet.
            const bool smooth = side.order > 0 && !end;
            Target     target{point,
                          frame.project(point),
                          smooth ? side.order : 0,
                                             {},
                                             {},
                          share * side.length / asked.size,
                          smooth && wish && side.order == 1};
            if (smooth)
            {
              target.jet          = side.jet(f);
              target.easing       = side.easing(f);
              const gp_Vec along3 = side.tangent(f);
              const gp_XY  along2(along3.Dot(frame.x), along3.Dot(frame.y));
              target.across = gp_XY(-along2.Y(), along2.X()) / along2.Modulus();
            }
            target.watch = watch;
            targets.push_back(target);
          };
          const double least = 1.0 / (crossingSteps * 16);
          place(0, least, true);
          place(1, least, true);
          for (size_t i = 0; i + 1 < cuts.size(); ++i)
          {
            const double from = cuts[i], to = cuts[i + 1];
            if (to - from < BRepExt::tolerance::parameterMerge)
              continue;
            for (int k = 0; k < placesPerSpan; ++k)
            {
              place(from + (to - from) * (k + 0.5) / placesPerSpan,
                    std::max(least, (to - from) / placesPerSpan),
                    false);
              // And between the places, where the surface is only read.
              place(from + (to - from) * (k + 1.0) / placesPerSpan,
                    std::max(least, (to - from) / placesPerSpan),
                    false,
                    true);
            }
          }
        }
      for (const gp_Pnt& point : asked.points)
        targets.push_back({point, frame.project(point), 0, {}, {}, 1.0 / spans});
      // Where two sides that are only to be met meet at a corner, a surface that holds both is
      // tangent there to the plane of their two tangents, whatever it is elsewhere. Within the
      // tolerance it need not be, and that plane is asked at the corner only once a surface of the
      // ladder has left the sides near the corner (see below). Not where they run on, where the two
      // give no plane.
      for (size_t i = 0; i < asked.sides.size(); ++i)
      {
        const Side &a = asked.sides[i], &b = asked.sides[(i + 1) % asked.sides.size()];
        if (!forced[i] || a.free || b.free || a.order > 0 || b.order > 0)
          continue;
        const gp_Vec ta = a.tangent(1).Normalized(), tb = b.tangent(0).Normalized();
        const gp_Vec normal = ta.Crossed(tb);
        if (normal.Magnitude() < std::sin(asked.limits.g1 * M_PI / 180))
          continue;
        const gp_Pnt corner = b.meet(0);
        Target       target{corner, frame.project(corner), 1, {}, {}, 1.0 / spans};
        SurfaceJet   jet{op, corner, ta, tb, {}, {}, {}, normal.Normalized()};
        target.jet = jet;
        targets.push_back(target);
      }
      const std::vector<char> fixed(net.poles.size(), 0);
      // Across a line through a corner where the rim's faces change their bend, the surface
      // may change its own: at the corner, where the rim needs it, and less and less the
      // farther along the line from it, so that the step stays at the rim and the surface
      // inside is of one bend across the line.
      std::vector<FitJump> jumps;
      for (const gp_Pnt& turn : turns)
      {
        const gp_XY corner = frame.project(turn);
        for (const bool acrossU : {true, false})
        {
          const std::vector<double>& knots = acrossU ? net.knotsV : net.knotsU;
          const double               here  = acrossU ? corner.Y() : corner.X();
          for (size_t k = degree; k + degree + 1 < knots.size(); ++k)
            for (int part = 0; knots[k + 1] > knots[k] && part < jumpPlaces; ++part)
            {
              const double along = knots[k] + (knots[k + 1] - knots[k]) * (part + 0.5) / jumpPlaces;
              const double away  = std::abs(along - here) / (jumpReach * asked.size);
              jumps.push_back({acrossU,
                               acrossU ? corner.X() : corner.Y(),
                               along,
                               jumpWeight * bendingWeight * asked.size * std::min(1.0, away * away)
                                 * (knots[k + 1] - knots[k]) / jumpPlaces});
            }
        }
      }
      // Fair inside the loop. Outside it the surface is only kept from running wild: counted
      // in full there, it would be held to the loop from both sides, like a plate clamped
      // along a line, and kink across it.
      const bool heightGraph =
        order > 0 && graphSupports
        && std::all_of(targets.begin(), targets.end(), [&](const Target& target) {
             return target.order == 0 || std::abs(target.jet->n.Dot(frame.z)) >= graphNormal;
           });
      const bool spatial =
        !heightGraph && order > 0 && asked.guides.empty() && nu * nv <= freePoles;
      // A net free in space and a height over the plane are fair by different measures and
      // come nearer different surfaces: the one is not a finer net of the other, and where
      // nets free in space have met the limits, a height is not tried after them.
      if (!previous.IsNull() && spatial != previousSpatial && best && best->measured.within)
        break;
      occ::handle<Geom_BSplineSurface> surface =
        seeded
              ? surfaceOf(net)
              : fit(
              net,
              fixed,
              targets,
              asked.size,
              reach,
              {},
              [&](double u, double v) { return within(gp_XY(u, v)) ? 1 : outsideWeight; },
              spatial ? nullptr : &up,
              jumps,
              order > 0 || asked.controls.open || !asked.guides.empty() || !asked.points.empty()
                    ? &frame
                    : nullptr,
              asked.flatness,
              fairRim,
              heightGraph ? graphMode : GraphFit::none);
      // The fit leaves of a curve what its net cannot follow fairly; the surface is then
      // moved onto the curves by the least movement of its poles, with knots added where
      // it takes them.
      // Not where the fit has left the tangent planes or the bends outside their limits: a
      // finer net is then tried in any case.
      // A boundary-only height graph is fitted by its weighted fair objective and corrected onto
      // the original rim: multiplier corrections amplify nearly dependent conditions.
      // Supported surfaces keep their original tangent/bend fit and correction gate.
      const bool moved =
        (fairRim || order > 0 || reach.g0 > asked.limits.g0)
        && measureAll(surface,
                          {asked.sides,
                       asked.guides,
                       asked.points,
                       {HUGE_VAL, asked.limits.g1, asked.limits.g2},
                       asked.size,
                       asked.flatness,
                       asked.controls},
                      [&](const Side&, double, const gp_Pnt& point) {
                        return std::optional<Place>(Place{frame.project(point), false});
                      })
             .within;
      // What the correction leaves at its own places where it had no room to make the net
      // finer, the net being as fine as it may make it already.
      double left = 0;
      if (moved)
      {
        const size_t knots = net.knotsU.size() + net.knotsV.size();
        // Per side, the places near its corners that the correction holds as well: those the
        // measure reads there (see measuredFractions) where the surface it leaves misses them.
        std::vector<std::vector<double>> nearCorners(asked.sides.size());
        const auto                       onto = [&](double span) {
          std::vector<RimPlace> onto;
          for (const auto* curves : {&asked.sides, &asked.guides})
            for (const Side& side : *curves)
            {
              if (side.free)
                continue;
              const int count = std::clamp(static_cast<int>(rimPlacesPerSpan * side.length / span),
                                           fewestRimPlaces,
                                           mostRimPlaces);
              std::vector<double> fractions;
              for (int i = 0; i <= count; ++i)
                fractions.push_back(static_cast<double>(i) / count);
              if (curves == &asked.sides)
                fractions.insert(fractions.end(),
                                 nearCorners[&side - asked.sides.data()].begin(),
                                 nearCorners[&side - asked.sides.data()].end());
              for (const double f : fractions)
              {
                const gp_Pnt point = side.meet(f);
                gp_XY        across(0, 0);
                if (side.order > 0 && f > 0 && f < 1)
                {
                  const gp_Vec along3 = side.tangent(f);
                  const gp_XY  along2(along3.Dot(frame.x), along3.Dot(frame.y));
                  across = gp_XY(-along2.Y(), along2.X()) / along2.Modulus();
                }
                onto.push_back({point, frame.project(point), across, side.order > 1});
              }
            }
          for (const gp_Pnt& point : asked.points)
            onto.push_back({point, frame.project(point), gp_XY(0, 0), false});
          return onto;
        };
        const auto correct = [&] {
          return correctRim(
            op,
            net,
            onto,
            up,
            std::min(rimShare * asked.limits.g0, exchange ? exchangeWithin : HUGE_VAL),
            std::min(enoughShare * asked.limits.g0, exchange ? 2 * exchangeWithin : HUGE_VAL),
            std::min(mostRimSpans,
                     (order > 0 ? rimFiner : rimFinerMet) * std::max(spansU, spansV)));
        };
        double worst = correct();
        // The correction's even places reach a corner of two sides that are only to be met no
        // nearer than its step: where the surface leaves a side nearer the corner than that by
        // more than the tolerance, those places are held as well and the correction made again,
        // and finer nets are asked the corner's tangent plane (see above).
        {
          const occ::handle<Geom_BSplineSurface> corrected = surfaceOf(net);
          bool                                   more      = false;
          const size_t                           count     = asked.sides.size();
          for (size_t s = 0; s < count; ++s)
          {
            const Side& side = asked.sides[s];
            if (side.free || side.order > 0)
              continue;
            for (const double f : measuredFractions())
            {
              const bool start = f < 0.5 / places, end = f > 1 - 0.5 / places;
              // The corner at that end, numbered by the side before it.
              const size_t corner = start ? (s + count - 1) % count : s;
              if ((!start && !end) || asked.sides[start ? corner : (s + 1) % count].order > 0)
                continue;
              const gp_Pnt point = side.meet(f);
              const gp_XY  uv    = nearest(corrected, frame.project(point), point);
              if (corrected->Value(uv.X(), uv.Y()).Distance(point) > asked.limits.g0)
              {
                nearCorners[s].push_back(f);
                forced[corner] = 1;
                more           = true;
              }
            }
          }
          if (more)
            worst = correct();
        }
        if (net.knotsU.size() + net.knotsV.size() == knots)
          left = worst;
        surface = surfaceOf(net);
      }
      std::vector<Rim> rims;
      for (const Side& side : asked.sides)
        rims.push_back(side.free
                         ? Rim{}
                         : preimage(surface,
                                    side,
                                    frame,
                                    std::min(asked.limits.g0, exchange ? exchangeWithin : HUGE_VAL),
                                    asked.size));
      closeGaps(asked.sides, rims);
      // A side is read where the face holds it; a guide and a point at the surface's nearest.
      Built built{3,
                      {},
                  measureAll(surface, asked, [&](const Side& side, double f, const gp_Pnt& point) {
                    const auto s = asked.boundaryIndex(side);
                    return std::optional<Place>(
                      s ? Place{rims[*s].onSurface->Value(side.parameter(f)).XY(), true}
                        : Place{frame.project(point), false});
                  })};
      for (size_t s = 0; s < rims.size(); ++s)
        rims[s].deviation =
          std::max(rims[s].deviation, built.measured.values[asked.sides[s].input * 4]);
      // What is not within the limits is kept for what it says of the miss.
      built.face = trimmed(surface, asked, rims, area < 0);
      if (fairRim)
        for (size_t i = 0; i < rims.size(); ++i)
          if (!asked.sides[i].free)
            built.spatialBoundary.push_back({asked.sides[i].input,
                                             rims[i].onSurface,
                                             asked.sides[i].from,
                                             asked.sides[i].to,
                                             asked.sides[i].originalReversed});
      for (size_t i = 0; hasGaps(asked) && i < rims.size(); ++i)
        built.gaps.push_back(asked.sides[i].free ? rims[i].onSurface : occ::handle<Geom2d_Curve>());
      // The actual trim must meet the original spatial edge too. Let the existing finite
      // ladder react to a missed trace instead of accepting a coarse preimage and refusing
      // only after placement. The final world-coordinate measurement remains authoritative.
      measureSpatialBoundary(built, asked, false);
      // Within every limit along its edges is not yet a patch: see BRepPatchExt_Inside.hxx. A net
      // that meets the edges with a surface that is not one to keep is too coarse for what
      // is asked along them, and what it cannot follow shows inside the loop; a finer one
      // may do. Points inside may ask for any bend, so with those only a fold is a fault.
      if (built.measured.within)
        if (const auto fault = BRepCheck_Analyzer(built.face).IsValid()
                                     ? insideFault(built.face, asked.points.empty() ? bend : 0)
                                     : std::optional<std::string>("it is not a valid face"))
        {
          built.measured.within = false;
          built.measured.unkept = true;
          built.measured.worst  = *fault;
          built.measured.excess = unkeptExcess;
        }
      const double change    = previous.IsNull() || spatial != previousSpatial
                                     ? HUGE_VAL
                                     : apart(previous,
                                      surface,
                                      // Twice to a span of the finer net each way: the
                                      // difference of the two is a spline of its spans.
                                      2 * spansU,
                                      2 * spansV,
                                      low,
                                      high,
                                      within);
      const bool   wasWithin = previousWithin;
      previous               = surface;
      previousWithin         = built.measured.within;
      previousSpatial        = spatial;
      // A net that meets every limit is not yet the surface asked for: that is the fairest of
      // all that meet them, which the nets of the ladder come nearer to as they become finer,
      // and a coarse net can meet the edges and still stand far from it inside (a seat of
      // eight sides: 9.2 off on six spans, 0.17 on every finer net). Two nets in a row that
      // lie within `settled` of each other have come as near it as finer ones would take
      // them: then the coarser of them that meets the limits is kept. Until then the finest
      // that meets them is, should the ladder or the nets allowed end first.
      const double before = best ? best->measured.excess : HUGE_VAL;
      if (built.measured.within)
        meeting.push_back(built);
      if (built.measured.within || (best && best->measured.within))
      {
        if (change <= settled * asked.size)
        {
          if (built.measured.within && !wasWithin)
            best = built;
          settledOn = true;
          break;
        }
        if (built.measured.within)
          best = built;
        continue;
      }
      if (built.measured.excess < before)
        best = built;
      // A few finer nets for a surface that is not one to keep, and no more: after those
      // it is the surface asked for that is not one to keep.
      if (built.measured.unkept)
      {
        if (++unkept >= mostUnkept)
          break;
        continue;
      }
      // A net that the fit left within its tangent planes and bends, and that was moved
      // onto its curves, is read at the places it was moved onto, and between them. Where
      // the fitted net was already as fine as the correction may make it, and the
      // correction's solve still leaves its own places further off than the tolerance, what
      // is left is what the fit left and the solve could not take, a share of the size that
      // finer nets of the ladder, fitted finer than their correction may go as well, do not
      // bring below it: the corner of rounds 2 / 3 / 5 at a hundred times its size (nets of
      // 25 x 20 and 31 x 26 spans for corrections bounded at 12 and 20) is left 3.7e-6 off by
      // the first and 6.0e-6 by the second (the limit 1e-6). Where the correction could make
      // the net finer, a finer net may meet them (a guide across a hole in a plate is met on
      // the third, at 32 x 40 spans).
      if (!fairRim && moved && left > asked.limits.g0)
      {
        best->measured.worst +=
          ". The bounded rim correction leaves the curves further off than the tolerance";
        break;
      }
      if (moved)
      {
        if (!fairRim && built.measured.excess > stillFalling * placeOnly)
        {
          best->measured.worst += ". A finer net came no nearer";
          break;
        }
        placeOnly = std::min(placeOnly, built.measured.excess);
      }
      // The coarsest nets miss by what they cannot follow at all, and by no rule, unless
      // it is by the very same.
      // But not, where curves are only met and the nets are cheap, while the best of them is
      // within twice its limits: with each net the fit stops where its own places are met,
      // and the misses of a ladder are not a falling series (one closed curve 0.5 sin 3a at a
      // tolerance of 1e-5 misses by 1.35, 4.5, 1.72 and 3.64 times that on the first four
      // nets and is met on the fifth).
      // Full Hessian errors need not fall monotonically as the narrow direction
      // gains spans. Exhaust its existing finite ladder rather than stop on one rise.
      if (!fairRim && graphMode != GraphFit::complete
          && ((++tried > unruly && built.measured.excess > stillFalling * before
               && (order > 0 || best->measured.excess > nearly))
              || std::abs(built.measured.excess - before) < same * before))
      {
        best->measured.worst += ". A finer net came no nearer";
        break;
      }
    }
    // Where no two nets in a row have come to agree, the ladder has not shown which of the
    // surfaces it made is near the fairest that meets the limits; of those that meet them, the
    // one that bends least is nearest it by the energy it is the least of.
    if (!settledOn && meeting.size() > 1)
    {
      double least = HUGE_VAL;
      for (const Built& candidate : meeting)
        if (const double bending = bendingInside(candidate.face); bending < least)
        {
          least = bending;
          best  = candidate;
        }
    }
    return *best;
  };
  if (nearlyQuadratic)
  {
    // Interpolating parabolic points does not make an exact parabola. Use the nearby
    // quadratic as a start on a cubic net, then actually meet the interpolated curves
    // by the same measured rim correction. It is not degree-two recognition.
    degree             = 3;
    reach              = Limits{fitReach * asked.size, asked.limits.g1, asked.limits.g2};
    const Built guided = build(false, static_cast<int>(std::size(ladder)));
    if (guided.measured.within)
      return guided;
    degree = 5;
    reach  = asked.limits;
  }
  // G1 on independently designed spline sheets leaves transverse curvature free.
  // Borrowing each sheet's G2 as a soft target can manufacture alternating bends in
  // an otherwise convex hole. Minimise the requested energy without those targets.
  // Analytic round supports retain their shape wishes (for example spherical caps).
  const bool splineRim = std::all_of(asked.sides.begin(), asked.sides.end(), [](const Side& side) {
    return !side.support.IsNull()
           && BRepAdaptor_Surface(side.support).GetType() == GeomAbs_BSplineSurface;
  });
  bool       wishable =
    !splineRim && std::any_of(asked.sides.begin(), asked.sides.end(), [](const Side& side) {
      return side.order == 1;
    });
  for (size_t i = 0; wishable && i < asked.sides.size(); ++i)
  {
    const Side &a = asked.sides[i], &b = asked.sides[(i + 1) % asked.sides.size()];
    if (a.order < 1 || b.order < 1 || a.support.IsSame(b.support))
      continue;
    const SurfaceJet before = a.jet(1 - cornerInset), after = b.jet(cornerInset);
    Miss             miss;
    std::tie(miss.g2, miss.g2Relative) = bendDifference(before,
                                                        after,
                                                        a.tangent(1 - cornerInset).Normalized(),
                                                        before.n.Dot(after.n) < 0 ? -1 : 1);
    wishable                           = miss.within(asked.limits, 2, asked.size);
  }
  if (wishable)
  {
    const Built wished = build(true, wishNets);
    if (wished.measured.within)
      return wished;
  }
  const Built strict = build(false, direction == 0 ? static_cast<int>(std::size(ladder)) : 1);
  if (!strict.measured.within && order == 2 && graphSupports)
  {
    // Approximate position/normal constraints do not guarantee the two Hessian
    // components omitted by the transverse row. Retry the same accuracy contract
    // with all three components, preserving already delivered constructions.
    graphMode = GraphFit::complete;
    try
    {
      const Built complete = build(false, direction == 0 ? static_cast<int>(std::size(ladder)) : 1);
      if (complete.measured.within)
        return complete;
    }
    catch (const Failure& failure)
    {
      if (failure.reason != no_solution)
        throw;
    }
    graphMode = GraphFit::transverse;
  }
  if (!strict.measured.within && !asked.guides.empty())
  {
    try
    {
      if (const auto reference = sectionReference(asked, frame))
        return *reference;
    }
    catch (const Failure& failure)
    {
      if (failure.reason == "cancelled")
        throw;
    }
    catch (const Standard_Failure&)
    {
    }
  }
  if (order > 0)
  {
    // A near miss can occur between the fit's own places even when every asked and
    // watched condition was met. Retry once with extra angular headroom; the request
    // and every delivered contact/interior/validity check keep their original limits.
    if (order == 1 && direction == 0 && !strict.measured.within && !strict.measured.unkept
        && strict.measured.excess < 2 && asked.guides.empty() && asked.points.empty())
    {
      reach.g1 *= .5;
      try
      {
        const Built tighter = build(false, static_cast<int>(std::size(ladder)));
        if (tighter.measured.within)
          return tighter;
      }
      catch (const Failure& failure)
      {
        if (failure.reason != no_solution)
          throw;
      }
    }
    if (direction == 0 && strict.measured.within && !turns.empty() && asked.guides.empty()
        && asked.points.empty())
    {
      // Repeated-radius rounds retain their clearance frame. A turned retry on
      // those supports added CPU time without a candidate worth keeping.
      std::vector<TopoDS_Face> rounds;
      std::vector<double>      radii;
      bool                     repeated = false;
      for (const Side& side : asked.sides)
      {
        if (side.support.IsNull() || BRepAdaptor_Surface(side.support).GetType() != GeomAbs_Cylinder
            || std::any_of(rounds.begin(), rounds.end(), [&](const TopoDS_Face& face) {
                 return face.IsSame(side.support);
               }))
          continue;
        const double radius = BRepAdaptor_Surface(side.support).Cylinder().Radius();
        repeated = repeated || std::any_of(radii.begin(), radii.end(), [&](double other) {
                     return std::abs(other - radius) < 1e-9 * std::max(other, radius);
                   });
        rounds.push_back(side.support);
        radii.push_back(radius);
      }
      const InsideCurvature before = insideCurvature(strict.face);
      if (!repeated && before.wrongSign() > .34 * bend)
      {
        std::optional<Built> turned;
        try
        {
          turned = overPlane(asked, flat, 1);
        }
        catch (const Failure&)
        {
          return strict;
        }
        if (turned->measured.within)
        {
          const InsideCurvature after = insideCurvature(turned->face);
          if (after.wrongSign() < before.wrongSign() && after.tightest() <= before.tightest()
              && insideStep(turned->face) <= insideStep(strict.face))
            return *turned;
        }
      }
    }
    return strict;
  }
  // Where curves are only met there is a second surface: fitted only as near as a fair one
  // over its net comes, a share of the loop's size, and moved onto the curves afterwards.
  // Brought nearer by the fit, a coarse net ripples to meet a curve (z = sin 2a on a radius
  // of 6 bends 0.25 at 1e-6 for the 0.056 of the fair surface); moved, it does not (0.060).
  // But for other loops the fit at the tolerance is the fairer (a C of eight sides: 0.06,
  // moved 0.45). Of the two that meet the limits the fairer is kept, the one fitted at the
  // tolerance unless the other bends clearly less. Not for a loop with a corner where its
  // sides meet at an angle: no surface follows that without a crease, and moving one onto
  // it only makes its net finer there (a polygon of 16 sides in space: 99 x 99 poles).
  for (size_t i = 0; asked.guides.empty() && i < asked.sides.size(); ++i)
  {
    const Side &a = asked.sides[i], &b = asked.sides[(i + 1) % asked.sides.size()];
    if (a.tangent(1).Angle(b.tangent(0)) * 180 / M_PI > asked.limits.g1)
      return strict;
  }
  reach       = Limits{fitReach * asked.size, asked.limits.g1, asked.limits.g2};
  Built loose = build(false, static_cast<int>(std::size(ladder)));
  if (!loose.measured.within && !strict.measured.within && degree == 3 && asked.guides.empty()
      && asked.points.empty())
  {
    degree = 5;
    loose  = build(false, static_cast<int>(std::size(ladder)));
  }
  if (!loose.measured.within)
    return strict;
  if (!strict.measured.within)
    return loose;
  return tightestInside(loose.face) < clearlyFairer * tightestInside(strict.face) ? loose : strict;
}

} // namespace BRepExt::patch_fill
