// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Over a chart of the loop's own: a fill that is a height over no plane, fitted over the loop
// spanned by a disc and laid flat (see BRepPatchExt_Chart.hxx).

#include <BRepPatchExt_Fill.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>

#include <BRepAlgoExt_PatchFlow.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <GeomAlgoExt_Cancel.hxx>
#include <GeomAlgoExt_Phase.hxx>

#include <BRepAdaptor_Surface.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepTools.hxx>
#include <Geom2dAPI_Interpolate.hxx>
#include <Geom2d_BSplineCurve.hxx>

#include <array>

namespace BRepExt::patch_fill
{
namespace
{

//! How far beside an edge, as a share of the diagonal of its face's parameters, the face is
//! looked for on either side of it. Chosen.
constexpr double sideProbe = 1e-4;

//! The places along a guide that the loop's chart is drawn through (see proxyChart).
constexpr int guidePins = 32;

//! The places the loop is spanned by round it for its chart (see proxyChart): the disc inside
//! has about a tenth of the square of that many nodes.
constexpr int chartPlaces = 160;

//! Which way along N x T, N the face's normal and T the side's tangent, a fill leaves a side it
//! is to be tangent to: away from the face, which lies to one side of the edge in the face's own
//! parameters. Zero where nothing is asked or that side is not found.
double leaving(const Side& side)
{
  if (side.free || side.order < 1 || side.support.IsNull())
    return 0;
  double     first, last;
  const auto onFace = BRep_Tool::CurveOnSurface(side.edge, side.support, first, last);
  if (onFace.IsNull())
    return 0;
  gp_Pnt2d uv;
  gp_Vec2d along;
  onFace->D1(side.parameter(0.5), uv, along);
  if (!(along.Magnitude() > 0))
    return 0;
  double u0, u1, v0, v1;
  BRepTools::UVBounds(side.support, u0, u1, v0, v1);
  const gp_Vec2d left =
    gp_Vec2d(-along.Y(), along.X()).Normalized() * (sideProbe * std::hypot(u1 - u0, v1 - v0));
  gp_Pnt point;
  gp_Vec su, sv;
  BRepAdaptor_Surface(side.support).D1(uv.X(), uv.Y(), point, su, sv);
  const gp_Vec across = side.faceJet(0.5).n.Crossed(side.tangent(0.5));
  for (const double sign : {1.0, -1.0})
  {
    const gp_Vec2d way = left * sign;
    if (BRepClass_FaceClassifier(side.support, uv.Translated(way), Precision::PConfusion()).State()
        == TopAbs_IN)
      return (su * way.X() + sv * way.Y()).Dot(across) > 0 ? -1 : 1;
  }
  return 0;
}

//! The fractions of a side at `pieces` even steps of its length, from its start to its end.
std::vector<double> evenSteps(const Side& side, int pieces)
{
  const int           fine = 16 * pieces;
  std::vector<double> length(fine + 1, 0);
  gp_Pnt              before = side.at(0);
  for (int i = 1; i <= fine; ++i)
  {
    const gp_Pnt here = side.at(static_cast<double>(i) / fine);
    length[i]         = length[i - 1] + before.Distance(here);
    before            = here;
  }
  std::vector<double> out{0};
  for (int k = 1; k < pieces; ++k)
  {
    const double want = length[fine] * k / pieces;
    const size_t i    = std::upper_bound(length.begin(), length.end(), want) - length.begin() - 1;
    const double part =
      length[i + 1] > length[i] ? (want - length[i]) / (length[i + 1] - length[i]) : 0;
    out.push_back((i + part) / fine);
  }
  out.push_back(1);
  return out;
}

//! A curve in the chart through `places` at `fractions`, for a curve parameter from `from` to
//! `to`.
std::function<gp_XY(double)> laidThrough(const std::vector<gp_XY>&  places,
                                         const std::vector<double>& fractions,
                                         double                     from,
                                         double                     to,
                                         double                     size)
{
  const int                                        count = static_cast<int>(places.size());
  const occ::handle<NCollection_HArray1<gp_Pnt2d>> points =
    new NCollection_HArray1<gp_Pnt2d>(1, count);
  const occ::handle<NCollection_HArray1<double>> at = new NCollection_HArray1<double>(1, count);
  for (int i = 0; i < count; ++i)
  {
    points->SetValue(i + 1, gp_Pnt2d(places[i]));
    at->SetValue(i + 1, fractions[i]);
  }
  Geom2dAPI_Interpolate through(points, at, false, BRepExt::tolerance::exactLength * size);
  through.Perform();
  if (!through.IsDone())
    fail(op, no_solution, "a side of the loop has no curve in its chart");
  const occ::handle<Geom2d_BSplineCurve> curve = through.Curve();
  return [curve, from, to](double t) { return curve->Value((t - from) / (to - from)).XY(); };
}

//! The Favard constant that bounds the best approximation of degree p splines on even knots of span
//! h: no more than K (h / pi)^(p + 1 - k) |f^(p + 1)| for the k-th derivative, K between 1 and
//! pi / 2 (Favard 1937; de Boor, A Practical Guide to Splines, 2001, XII). The larger is taken.
constexpr double favard = M_PI / 2;

//! A cost forecast using the higher derivatives of the osculating circle as a surrogate.
//! Those derivatives are not bounded by curvature for a general surface or chart. This is
//! a computation policy, not an approximation theorem or a necessary pole-count bound.
double netSpacing(double bend, int order, int degree, const Limits& limits)
{
  if (!(bend > 0))
    return HUGE_VAL;
  double span =
    M_PI * std::pow(fitShare * limits.g0 / (favard * std::pow(bend, degree)), 1.0 / (degree + 1));
  if (order > 0)
    span = std::min(
      span,
      M_PI * std::pow(fitShare * std::sin(limits.g1 * M_PI / 180) / favard, 1.0 / degree) / bend);
  if (order > 1)
    span =
      std::min(span, M_PI * std::pow(fitShare * limits.g2 / favard, 1.0 / (degree - 1)) / bend);
  return span;
}

//! The span each node of `chart` asks of a net (netSpacing): inside, by how the disc bends
//! there; along the rim, by how the side's curve and its face bend, to the order the side asks.
//! `steps` are the fractions of each side at the rim's nodes, from `firstNode` of the side on.
std::vector<double> nodeSpacing(const Asked&                            asked,
                                const ConformalChart&                   chart,
                                const std::vector<std::vector<double>>& steps,
                                const std::vector<size_t>&              firstNode,
                                int                                     degree)
{
  const size_t              count = asked.sides.size();
  const std::vector<double> bends = nodeBends(chart);
  std::vector<double>       spacing(chart.point.size());
  for (size_t i = 0; i < spacing.size(); ++i)
    spacing[i] = std::min(asked.size, netSpacing(bends[i], 0, degree, asked.limits));
  for (size_t s = 0; s < count; ++s)
  {
    const Side& side = asked.sides[s];
    for (size_t j = 0; j < steps[s].size(); ++j)
    {
      const size_t node = j + 1 < steps[s].size() ? firstNode[s] + j : firstNode[(s + 1) % count];
      const double f    = steps[s][j];
      // The side's bend from the circle through its places beside this one, at even steps of its
      // length: its parameter may slow to nothing (a trailing edge) where its derivatives say
      // nothing of the bend, and a corner, at a side's end, is no bend.
      double curve = 0;
      if (steps[s].size() >= 3)
      {
        const size_t middle = std::clamp<size_t>(j, 1, steps[s].size() - 2);
        const gp_Pnt before = side.at(steps[s][middle - 1]), here = side.at(steps[s][middle]),
                     after = side.at(steps[s][middle + 1]);
        const double product =
          before.Distance(here) * here.Distance(after) * before.Distance(after);
        if (product > 0)
          curve = 2 * gp_Vec(here, before).Crossed(gp_Vec(here, after)).Magnitude() / product;
      }
      const double face =
        side.order > 0 ? sharpest(side.faceJet(std::clamp(f, cornerInset, 1 - cornerInset))) : 0;
      spacing[node] =
        std::min({spacing[node],
                  netSpacing(std::max({curve, face, bends[node]}), 0, degree, asked.limits),
                  netSpacing(std::max(face, bends[node]), side.order, degree, asked.limits)});
    }
  }
  return spacing;
}

//! Whether the spans a net over the chart needs at least make more poles than a net free in
//! space is let have: each way, over the chart's extent, the span each node asks as the chart
//! stretches it there, the finest across the other way.
bool tooManyPoles(const ConformalChart& chart, const std::vector<double>& spacing, int degree)
{
  const std::vector<double> stretch = nodeStretch(chart);
  constexpr int             bins    = 256;
  const auto                needed  = [&](int way, double from, double to) {
    std::vector<double> span(bins, HUGE_VAL);
    for (size_t i = 0; i < chart.flat.size(); ++i)
    {
      const int bin =
        std::clamp(static_cast<int>((chart.flat[i].Coord(way) - from) / (to - from) * bins),
                   0,
                   bins - 1);
      span[bin] = std::min(span[bin], spacing[i] * stretch[i]);
    }
    double spans = 0;
    for (const double one : span)
      if (std::isfinite(one))
        spans += (to - from) / bins / one;
    return std::max(fewestSpans, static_cast<int>(std::ceil(spans)));
  };
  gp_XY low = chart.flat.front(), high = low;
  for (const gp_XY& at : chart.flat)
  {
    low.SetCoord(std::min(low.X(), at.X()), std::min(low.Y(), at.Y()));
    high.SetCoord(std::max(high.X(), at.X()), std::max(high.Y(), at.Y()));
  }
  const gp_XY margin = (high - low) * marginShare;
  const int   spansU = needed(1, low.X() - margin.X(), high.X() + margin.X()),
            spansV   = needed(2, low.Y() - margin.Y(), high.Y() + margin.Y());
  if ((spansU + degree) * (spansV + degree) > freePoles)
    return true;
  return false;
}

//! The fill over a chart of the loop's own once the chart is made: the sides, guides and points
//! laid into it, and a net that starts as the coarsest of the ladder over the chart.
class DiscFill
{
public:
  DiscFill(const Asked&                            asked,
           const Frame&                            frame,
           const ConformalChart&                   chart,
           const std::vector<std::vector<double>>& steps,
           const std::vector<size_t>&              firstNode,
           const std::vector<double>&              sign,
           size_t                                  rimNodes,
           int                                     degree)
      : asked(asked),
        frame(frame),
        chart(chart),
        count(asked.sides.size()),
        sign(sign),
        degree(degree),
        layout(count),
        loop(chart.flat.begin(), chart.flat.begin() + rimNodes),
        within(loop, 0)
  {
    for (size_t s = 0; s < count; ++s)
    {
      std::vector<gp_XY> places;
      for (size_t k = 0; k + 1 < steps[s].size(); ++k)
        places.push_back(chart.flat[firstNode[s] + k]);
      places.push_back(chart.flat[firstNode[(s + 1) % count]]);
      layout[s] = laidThrough(places, steps[s], asked.sides[s].from, asked.sides[s].to, asked.size);
    }
    for (const Side& guide : asked.guides)
    {
      constexpr int       guidePlaces = 64;
      std::vector<gp_XY>  places;
      std::vector<double> fractions;
      for (int i = 0; i <= guidePlaces; ++i)
      {
        fractions.push_back(static_cast<double>(i) / guidePlaces);
        places.push_back(onDisc(guide.at(fractions.back())));
      }
      guideLayout.push_back(laidThrough(places, fractions, 0, 1, asked.size));
    }
    for (const gp_Pnt& point : asked.points)
      pointPlaces.push_back(onDisc(point));

    low  = loop.front();
    high = low;
    for (const gp_XY& at : chart.flat)
    {
      low.SetCoord(std::min(low.X(), at.X()), std::min(low.Y(), at.Y()));
      high.SetCoord(std::max(high.X(), at.X()), std::max(high.Y(), at.Y()));
    }
    const double width = high.X() - low.X(), depth = high.Y() - low.Y(),
                 larger = std::max(width, depth);
    const gp_XY margin(marginShare * width, marginShare * depth);

    // The net starts as the coarsest of the ladder over the chart and is made finer span by span
    // where a place it is held to is missed, as the Adjacent fit is (see adjacentOver): the nose
    // of a wing asks for spans far shorter than the rest of its chart does.
    const auto even = [&](std::vector<double>& breaks, double from, double to, double extent) {
      const int spans =
        std::max(fewestSpans, static_cast<int>(std::lround(ladder[0] * extent / larger)));
      for (int i = 0; i <= spans; ++i)
        breaks.push_back(from + (to - from) * i / spans);
    };
    even(breaksU, low.X() - margin.X(), high.X() + margin.X(), width);
    even(breaksV, low.Y() - margin.Y(), high.Y() + margin.Y(), depth);
  }

  DiscFill(const DiscFill&)            = delete;
  DiscFill& operator=(const DiscFill&) = delete;

  //! The fill, made finer span by span where a place it is held to is missed, as the Adjacent
  //! fit is (see adjacentOver): the nose of a wing asks for spans far shorter than the rest of
  //! its chart does. Nothing where no net meets every limit.
  std::optional<Built> run();

private:
  //! Where a guide or a point lies on the disc in space; fails where that is outside the loop.
  gp_XY onDisc(const gp_Pnt& point) const;

  //! The knots of a net of `degree` with these breaks, clamped at both ends.
  std::vector<double> knots(const std::vector<double>& breaks) const;

  //! `surface` trimmed by the loop, measured along its edges and checked to be regular, valid
  //! and one to keep where the face uses it.
  Built certify(const occ::handle<Geom_BSplineSurface>& surface) const;

  //! The net of the current breaks, flat in the frame's plane.
  SplineNet flatNet() const;

  //! The shortest span of the current breaks.
  double smallestSpan() const;

  //! What the fit is held to along the sides and guides and at the points, about
  //! `placesPerSpan` places to the shortest span; per target, in `along`, the way along its side
  //! in space where it asks a tangent plane.
  void placeTargets(std::vector<Target>& targets,
                    std::vector<gp_Vec>& along,
                    double               smallest) const;

  //! The breaks made finer in every span where the fit misses a place it was held to, or a side
  //! is missed between them: halved, up to the finest net of the ladder each way. Whether to go
  //! on: not where no span was, nor where a finer net would have more poles than allowed.
  bool refine(const occ::handle<Geom_BSplineSurface>& surface,
              const std::vector<Target>&              targets,
              const std::vector<gp_Vec>&              along,
              const Built&                            built,
              bool                                    split,
              double                                  smallest,
              int&                                    unkept);

  const Asked&                              asked;
  const Frame&                              frame;
  const ConformalChart&                     chart;
  const size_t                              count;
  const std::vector<double>&                sign;
  const int                                 degree;
  std::vector<std::function<gp_XY(double)>> layout;
  const std::vector<gp_XY>                  loop;
  const Inside                              within;
  std::vector<std::function<gp_XY(double)>> guideLayout;
  std::vector<gp_XY>                        pointPlaces;
  gp_XY                                     low, high;
  std::vector<double>                       breaksU, breaksV;
};

std::optional<Built> DiscFill::run()
{
  std::optional<Built>      best;
  occ::handle<Geom_Surface> previous;
  int                       idle = 0, unkept = 0;
  for (double before = HUGE_VAL;;)
  {
    SplineNet           net      = flatNet();
    const double        smallest = smallestSpan();
    std::vector<Target> targets;
    std::vector<gp_Vec> along;
    placeTargets(targets, along, smallest);
    const occ::handle<Geom_BSplineSurface> surface = fit(
      net,
      std::vector<char>(net.poles.size(), 0),
      targets,
      asked.size,
      asked.limits,
      {},
      [&](double u, double v) { return within(gp_XY(u, v)) ? 1 : outsideWeight; },
      nullptr,
      {},
      &frame);
    Built        built  = certify(surface);
    const double change = previous.IsNull() ? HUGE_VAL
                                            : apart(previous,
                                                    surface,
                                                    2 * static_cast<int>(breaksU.size() - 1),
                                                    2 * static_cast<int>(breaksV.size() - 1),
                                                    low,
                                                    high,
                                                    within);
    previous            = surface;
    // As over the mean plane: a net that meets every limit is kept once the next finer one lies
    // within `settled` of it, and until then the finest that meets them.
    const bool split = built.measured.within || (best && best->measured.within);
    if (split)
    {
      if (change <= settled * asked.size)
      {
        if (!best->measured.within)
          best = std::move(built);
        break;
      }
      if (built.measured.within)
        best = std::move(built);
    }
    else if (!best || built.measured.excess < best->measured.excess)
      best = built;
    // Nets that come no nearer, more than the coarsest few may: what is left is not the net's.
    if (split)
    {
    }
    else if (built.measured.excess > stillFalling * before)
    {
      if (++idle > unruly)
      {
        best->measured.worst += ". A finer net came no nearer";
        break;
      }
    }
    else
    {
      idle   = 0;
      before = built.measured.excess;
    }
    if (!refine(surface, targets, along, built, split, smallest, unkept))
      break;
  }
  if (!best || !best->measured.within)
    return std::nullopt;
  return best;
}

gp_XY DiscFill::onDisc(const gp_Pnt& point) const
{
  double     distance = 0;
  const auto at       = onChart(chart, point, distance);
  if (!inside(loop, at, outsideSlack * asked.size))
    fail(op, invalid_input, "a guide or a point lies outside the boundary");
  return at;
}

std::vector<double> DiscFill::knots(const std::vector<double>& breaks) const
{
  std::vector<double> out(degree, breaks.front());
  out.insert(out.end(), breaks.begin(), breaks.end());
  out.insert(out.end(), degree, breaks.back());
  return out;
}

Built DiscFill::certify(const occ::handle<Geom_BSplineSurface>& surface) const
{
  std::vector<Rim> rims;
  for (size_t s = 0; s < count; ++s)
    rims.push_back(
      asked.sides[s].free
        ? Rim{}
        : preimage(surface, asked.sides[s], frame, asked.limits.g0, asked.size, {}, layout[s]));
  closeGaps(asked.sides, rims);
  std::vector<gp_XY> around;
  for (size_t s = 0; s < count; ++s)
    for (int k = 0; k < 64; ++k)
      around.push_back(rims[s].onSurface->Value(asked.sides[s].parameter(k / 64.0)).XY());
  double turning = 0;
  for (size_t i = 0; i < around.size(); ++i)
    turning += around[i].Crossed(around[(i + 1) % around.size()]);
  Built built{3,
              {},
              measureAll(surface,
                         asked,
                         [&](const Side& side, double f, const gp_Pnt&) -> std::optional<Place> {
                           if (const auto s = asked.boundaryIndex(side))
                             return Place{rims[*s].onSurface->Value(side.parameter(f)).XY(), true};
                           // A guide is read at the surface's nearest from its chart place.
                           for (size_t g = 0; g < asked.guides.size(); ++g)
                             if (&side == &asked.guides[g])
                               return Place{guideLayout[g](f), false};
                           return std::nullopt;
                         })};
  for (size_t s = 0; s < count; ++s)
    rims[s].deviation =
      std::max(rims[s].deviation, built.measured.values[asked.sides[s].input * 4]);
  built.face = trimmed(surface, asked, rims, turning < 0);
  for (size_t s = 0; s < count; ++s)
  {
    const Side& side = asked.sides[s];
    if (!side.free)
      built.spatialBoundary.push_back(
        {side.input, rims[s].onSurface, side.from, side.to, side.originalReversed});
    if (hasGaps(asked))
      built.gaps.push_back(side.free ? rims[s].onSurface : occ::handle<Geom2d_Curve>());
  }
  measureSpatialBoundary(built, asked, false);
  if (built.measured.within)
  {
    // The cells the face uses: those that hold a place inside the loop or of its rim.
    const auto used = [&](double u0, double u1, double v0, double v1) {
      for (const gp_XY& at : around)
        if (u0 <= at.X() && at.X() <= u1 && v0 <= at.Y() && at.Y() <= v1)
          return true;
      for (int i = 0; i <= 4; ++i)
        for (int j = 0; j <= 4; ++j)
          if (within(gp_XY(u0 + (u1 - u0) * i / 4, v0 + (v1 - v0) * j / 4)))
            return true;
      return false;
    };
    std::optional<std::string> fault;
    if (!regularMargin(surface, used))
      fault = "it is not shown regular over the loop's chart";
    else if (!BRepCheck_Analyzer(built.face).IsValid())
      fault = "it is not a valid face";
    else
      fault = insideFault(built.face, asked.points.empty() ? rimBend(asked) : 0);
    if (fault)
    {
      built.measured.within = false;
      built.measured.unkept = true;
      built.measured.worst  = *fault;
      built.measured.excess = unkeptExcess;
    }
  }
  return built;
}

SplineNet DiscFill::flatNet() const
{
  SplineNet  net{degree, degree, knots(breaksU), knots(breaksV), {}};
  const int  nu = net.countU(), nv = net.countV();
  const auto mean = [&](const std::vector<double>& knots, int i) {
    double sum = 0;
    for (int k = 1; k <= degree; ++k)
      sum += knots[i + k];
    return sum / degree;
  };
  for (int i = 0; i < nu; ++i)
    for (int j = 0; j < nv; ++j)
      net.poles.push_back(frame.origin.XYZ() + frame.x.XYZ() * mean(net.knotsU, i)
                          + frame.y.XYZ() * mean(net.knotsV, j));
  return net;
}

double DiscFill::smallestSpan() const
{
  double smallest = HUGE_VAL;
  for (const auto* breaks : {&breaksU, &breaksV})
    for (size_t k = 1; k < breaks->size(); ++k)
      smallest = std::min(smallest, (*breaks)[k] - (*breaks)[k - 1]);
  return smallest;
}

void DiscFill::placeTargets(std::vector<Target>& targets,
                            std::vector<gp_Vec>& along,
                            double               smallest) const
{
  for (size_t s = 0; s < count; ++s)
  {
    const Side& side = asked.sides[s];
    if (side.free)
      continue;
    const int places =
      std::clamp(static_cast<int>(std::ceil(placesPerSpan * side.length / smallest)),
                 fewestRimPlaces,
                 mostRimPlaces);
    const double step = 1.0 / places;
    // Next to a gap the face has its say at the end as along the rest (see overPlane).
    const auto quiet = [&](size_t other) {
      const Side& near = asked.sides[other % count];
      return near.free || near.support.IsSame(side.support);
    };
    const bool lone       = hasGaps(asked);
    const bool quietStart = lone && quiet(s + count - 1), quietEnd = lone && quiet(s + 1);
    const auto add = [&](double f, bool end, bool watch) {
      if (end && ((f == 0 && quietStart) || (f == 1 && quietEnd)))
      {
        f   = std::clamp(f, cornerInset, 1 - cornerInset);
        end = false;
      }
      const bool smooth = side.order > 0 && !end;
      Target     target{side.meet(f),
                    layout[s](side.parameter(f)),
                    smooth ? side.order : 0,
                        {},
                        {},
                    step * side.length / asked.size};
      target.watch = watch;
      if (smooth)
      {
        target.jet           = side.jet(f);
        target.easing        = side.easing(f);
        const double h       = std::min({step / 4, f, 1 - f});
        const gp_XY  tangent = layout[s](side.parameter(f + h)) - layout[s](side.parameter(f - h));
        target.across        = gp_XY(-tangent.Y(), tangent.X()) / tangent.Modulus();
        if (sign[s] != 0)
          target.leave = side.faceJet(f).n.Crossed(side.tangent(f)).Normalized() * sign[s];
      }
      targets.push_back(target);
      along.push_back(smooth ? side.tangent(f).Normalized() : gp_Vec(0, 0, 0));
    };
    add(0, true, false);
    add(1, true, false);
    for (int k = 0; k < places; ++k)
    {
      add((k + 0.5) * step, false, false);
      if (k + 1 < places)
        add((k + 1) * step, false, true);
    }
  }
  for (size_t g = 0; g < asked.guides.size(); ++g)
  {
    const Side& guide = asked.guides[g];
    const int   places =
      std::clamp(static_cast<int>(std::ceil(placesPerSpan * guide.length / smallest)),
                 fewestRimPlaces,
                 mostRimPlaces);
    for (int i = 0; i <= 2 * places; ++i)
    {
      const double f = static_cast<double>(i) / (2 * places);
      targets.push_back(
        {guide.meet(f), guideLayout[g](f), 0, {}, {}, guide.length / places / asked.size});
      targets.back().watch = i % 2 != 0;
      along.push_back(gp_Vec(0, 0, 0));
    }
  }
  for (size_t i = 0; i < asked.points.size(); ++i)
  {
    targets.push_back({asked.points[i], pointPlaces[i], 0, {}, {}, 1});
    along.push_back(gp_Vec(0, 0, 0));
  }
}

bool DiscFill::refine(const occ::handle<Geom_BSplineSurface>& surface,
                      const std::vector<Target>&              targets,
                      const std::vector<gp_Vec>&              along,
                      const Built&                            built,
                      bool                                    split,
                      double                                  smallest,
                      int&                                    unkept)
{
  // Finer in every span where the fit misses a place it was held to, or a side is missed
  // between them: halved, up to the finest net of the ladder each way.
  std::vector<char> splitU(breaksU.size() - 1, split), splitV(breaksV.size() - 1, split);
  const auto        mark = [&](const gp_XY& uv) {
    const auto span = [](const std::vector<double>& breaks, double at) {
      return std::clamp(
        static_cast<int>(std::upper_bound(breaks.begin(), breaks.end(), at) - breaks.begin()) - 1,
        0,
        static_cast<int>(breaks.size()) - 2);
    };
    splitU[span(breaksU, uv.X())] = splitV[span(breaksV, uv.Y())] = 1;
  };
  for (size_t t = 0; !split && t < targets.size(); ++t)
  {
    const Target& target = targets[t];
    const double  share  = target.watch ? watchShare : fitShare;
    bool          missed =
      surface->Value(target.uv.X(), target.uv.Y()).Distance(target.point) > share * asked.limits.g0;
    if (!missed && target.order > 0)
    {
      const SurfaceJet mine   = surfaceJet(op, surface, target.uv.X(), target.uv.Y());
      const double     cosine = mine.n.Dot(target.jet->n);
      missed = std::acos(std::min(1.0, std::abs(cosine))) * 180 / M_PI > share * asked.limits.g1;
      if (!missed && target.order > 1)
      {
        const auto [bend, relative] =
          bendDifference(*target.jet, mine, along[t], cosine < 0 ? -1 : 1);
        missed = relative > share * asked.limits.g2
                 && bend * asked.size > share * asked.limits.g2 * flatShare;
      }
    }
    if (missed)
      mark(target.uv);
  }
  for (size_t s = 0; s < count && !split; ++s)
  {
    const Side& side = asked.sides[s];
    if (side.free)
      continue;
    Miss miss;
    miss.g0         = built.measured.values[side.input * 4];
    miss.g1         = built.measured.values[side.input * 4 + 1];
    miss.g2Relative = built.measured.values[side.input * 4 + 2];
    miss.g2         = built.measured.values[side.input * 4 + 3];
    if (miss.within(asked.limits, side.order, asked.size))
      continue;
    const int places =
      std::clamp(static_cast<int>(std::ceil(4 * placesPerSpan * side.length / smallest)),
                 fewestRimPlaces,
                 4 * mostRimPlaces);
    for (int i = 0; i <= places; ++i)
      mark(layout[s](side.parameter(static_cast<double>(i) / places)));
  }
  // A net that meets the edges with a surface that is not one to keep is too coarse for what
  // is asked along them, wherever that shows: a few finer nets, and no more.
  if (!split && built.measured.unkept && ++unkept <= mostUnkept)
  {
    std::fill(splitU.begin(), splitU.end(), 1);
    std::fill(splitV.begin(), splitV.end(), 1);
  }
  bool                               finer = false;
  std::array<std::vector<double>, 2> next;
  for (int way = 0; way < 2; ++way)
  {
    const auto& breaks = way ? breaksV : breaksU;
    const auto& marked = way ? splitV : splitU;
    next[way]          = {breaks.front()};
    int spans          = static_cast<int>(breaks.size()) - 1;
    for (size_t i = 0; i + 1 < breaks.size(); ++i)
    {
      if (marked[i] && spans < std::end(ladder)[-1])
      {
        next[way].push_back((breaks[i] + breaks[i + 1]) / 2);
        ++spans;
        finer = true;
      }
      next[way].push_back(breaks[i + 1]);
    }
  }
  // A net free in space is let have no more poles than over the mean plane (see freePoles).
  if ((next[0].size() + degree - 1) * (next[1].size() + degree - 1) > freePoles)
    return false;
  breaksU = next[0];
  breaksV = next[1];
  return finer;
}

} // namespace

std::optional<Built> overChart(const Asked& asked, const Frame& frame)
{
  PatchPhase   phase("chart");
  const size_t count     = asked.sides.size();
  double       perimeter = 0;
  for (const Side& side : asked.sides)
    perimeter += side.length;
  std::vector<ChartRim>            rim;
  std::vector<std::vector<double>> steps(count);
  std::vector<size_t>              firstNode(count);
  std::vector<double>              sign(count);
  for (size_t s = 0; s < count; ++s)
  {
    const Side& side = asked.sides[s];
    const int   pieces =
      std::max(1, static_cast<int>(std::lround(chartPlaces * side.length / perimeter)));
    steps[s]     = evenSteps(side, pieces);
    firstNode[s] = rim.size();
    sign[s]      = leaving(side);
    for (int k = 0; k < pieces; ++k)
    {
      const double f = (steps[s][k] + steps[s][k + 1]) / 2;
      rim.push_back({side.at(steps[s][k]),
                     sign[s] == 0
                       ? gp_Vec(0, 0, 0)
                       : side.faceJet(f).n.Crossed(side.tangent(f)).Normalized() * sign[s]});
    }
  }
  std::vector<gp_Pnt> through = asked.points;
  for (const Side& guide : asked.guides)
    for (int i = 0; i <= guidePins; ++i)
      through.push_back(guide.at(static_cast<double>(i) / guidePins));
  int order = 0;
  for (const Side& side : asked.sides)
    order = std::max(order, side.order);
  const int                 degree  = order > 0 || !asked.guides.empty() ? 5 : 3;
  const ConformalChart      chart   = proxyChart(rim, through);
  const std::vector<double> spacing = nodeSpacing(asked, chart, steps, firstNode, degree);
  // The spans a net over the chart needs at least: each way, over the chart's extent, the span
  // each node asks as the chart stretches it there, the finest across the other way. Where those
  // make more poles than a net free in space is let have, no net over the chart is tried.
  if (tooManyPoles(chart, spacing, degree))
    return std::nullopt;
  DiscFill fill(asked, frame, chart, steps, firstNode, sign, rim.size(), degree);
  return fill.run();
}

} // namespace BRepExt::patch_fill
