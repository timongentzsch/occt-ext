// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Adjacent over a conformal chart: the faces round the loop flattened together with it, and the
// fill fitted over that chart, finer span by span where what it is held to is not met.

#include <BRepPatchExt_Fill.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>

#include <BRepAlgoExt_PatchFlow.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <GeomAlgoExt_Cancel.hxx>
#include <GeomAlgoExt_PatchCompact.hxx>
#include <GeomAlgoExt_PatchRim.hxx>
#include <GeomAlgoExt_Phase.hxx>

#include <BRepCheck_Analyzer.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>

namespace BRepExt::patch_fill
{
namespace
{

//! The triangulation that is flattened: at least this many triangles across the loop's
//! narrowest part (its area over its length round).
constexpr double trianglesAcross = 4;

//! How far the triangulation may stand off the face, as a share of the tolerance: where the
//! face bends its triangles are as small as that makes them, and the chart is as near the
//! face's own as its triangles are to the face.
constexpr double chartDeflection = 0.1;

//! The loop's places in the chart, in turn.
std::vector<gp_XY> boundaryLoop(const ConformalChart& chart)
{
  std::vector<gp_XY> loop;
  for (int b : chart.boundary)
    loop.push_back(chart.flat[b]);
  return loop;
}

// Where the result is read inside the loop: every node and the middle of every triangle,
// each with the Natural fill's normal there, which a fill that does not fold keeps to.
struct Check
{
  gp_XY  flat;
  gp_Vec normal;
};

std::vector<Check> interiorChecks(const ConformalChart&            chart,
                                  const occ::handle<Geom_Surface>& source)
{
  std::vector<Check> checks;
  const auto         check = [&](const gp_XY& flat, const gp_XY& uv) {
    gp_Pnt point;
    gp_Vec su, sv;
    source->D1(uv.X(), uv.Y(), point, su, sv);
    checks.push_back({flat, su.Crossed(sv)});
  };
  for (size_t i = 0; i < chart.flat.size(); ++i)
    check(chart.flat[i], chart.source[i]);
  for (const auto& t : chart.triangles)
    check((chart.flat[t[0]] + chart.flat[t[1]] + chart.flat[t[2]]) / 3,
          (chart.source[t[0]] + chart.source[t[1]] + chart.source[t[2]]) / 3);
  return checks;
}

//! The Adjacent fit over the chart: a net made finer span by span where what it is held to is
//! not met, from the coarsest of the ladder.
class ChartFill
{
public:
  ChartFill(const Asked&                                  asked,
            const Frame&                                  frame,
            const Built&                                  natural,
            const ConformalChart&                         chart,
            const std::vector<occ::handle<Geom2d_Curve>>& held,
            double                                        area)
      : asked(asked),
        frame(frame),
        natural(natural),
        chart(chart),
        held(held),
        area(area),
        loop(boundaryLoop(chart)),
        onBoundary(chart.flat.size(), 0),
        within(loop, 0),
        source(BRep_Tool::Surface(natural.face)),
        sourceProject(source),
        checks(interiorChecks(chart, source))
  {
    for (int b : chart.boundary)
      onBoundary[b] = 1;
    gp_XY low = chart.flat.front(), high = low;
    for (const gp_XY& at : chart.flat)
    {
      low.SetCoord(std::min(low.X(), at.X()), std::min(low.Y(), at.Y()));
      high.SetCoord(std::max(high.X(), at.X()), std::max(high.Y(), at.Y()));
    }
    const double width = high.X() - low.X(), depth = high.Y() - low.Y(),
                 larger = std::max(width, depth);
    const gp_XY margin(marginShare * width, marginShare * depth);
    for (const Side& side : asked.sides)
      order = std::max(order, side.order);
    // As over the mean plane: degree 5 for tangency, cubic where edges are only met.
    degree = order > 0 ? 5 : 3;
    // The net starts as the coarsest of the ladder and is made finer span by span where what
    // it is held to is not met: the knots follow the loop's geometry and the tolerance.
    const auto even = [&](std::vector<double>& breaks, double from, double to, double extent) {
      const int spans =
        std::max(fewestSpans, static_cast<int>(std::lround(ladder[0] * extent / larger)));
      for (int i = 0; i <= spans; ++i)
        breaks.push_back(from + (to - from) * i / spans);
    };
    even(breaksU, low.X() - margin.X(), high.X() + margin.X(), width);
    even(breaksV, low.Y() - margin.Y(), high.Y() + margin.Y(), depth);
  }

  ChartFill(const ChartFill&)            = delete;
  ChartFill& operator=(const ChartFill&) = delete;

  //! The fill, the coarsest net first.
  Built run();

private:
  //! A fill read inside the loop: the result, and where and whether it folds over.
  struct Checked
  {
    Built             built;
    std::vector<char> folded;
    bool              folds;
  };

  //! The places of side `s` in the chart, by its curve's parameter.
  std::function<gp_XY(double)> layout(size_t s) const;

  //! Where a point inside the loop lies in the chart, through the Natural fill's parameters.
  //! Interior constraints belong to the Natural fill's chart as well. Keep that
  //! correspondence fixed while fitting; reprojecting onto each candidate would
  //! allow a constraint to slide away from the chart that defines the flow.
  gp_XY interiorPlace(const gp_Pnt& point) const;

  //! The knots of a net of `degree` with these breaks, clamped at both ends.
  std::vector<double> knots(const std::vector<double>& breaks) const;

  //! `surface` trimmed by the loop, measured along its edges and read inside it for folds and
  //! for how square its U and V are; `spatial` keeps the edges' curves in space for measuring.
  Checked certify(const occ::handle<Geom_BSplineSurface>& surface, bool spatial = false) const;

  //! The net of the current breaks, flat in the frame's plane.
  SplineNet flatNet() const;

  //! The shortest span of the current breaks.
  double smallestSpan() const;

  //! The Natural fill inside, wished for at the node nearest the middle of every part of
  //! every span, `placesPerSpan` parts each way: as many wishes as the net can follow.
  void wishInterior(std::vector<Target>& targets) const;

  //! The places along the sides that the fit is held to, about `placesPerSpan` to the shortest
  //! span.
  void placeSides(std::vector<Target>& targets, double smallest) const;

  //! The places along the guides, and the points.
  void placeGuides(std::vector<Target>& targets, double smallest) const;

  //! Where curves and points are only met: the fit moved onto them.
  void correctPositions(SplineNet& net, occ::handle<Geom_BSplineSurface>& surface) const;

  //! The breaks made finer in every span where a place along an edge is missed, the surface
  //! folds or the shape wished for is not followed; whether any span was.
  bool refine(const occ::handle<Geom_BSplineSurface>& surface,
              const std::vector<Target>&              targets,
              const Built&                            built,
              const std::vector<char>&                folded,
              bool                                    folds,
              double                                  wishedWithin,
              double                                  smallest,
              int&                                    unkept);

  //! Where curves are only met: the best fill with fewer poles, where that keeps every check.
  void compactBest(std::optional<Built>& best, const std::vector<Target>& bestTargets) const;

  const Asked&                                  asked;
  const Frame&                                  frame;
  const Built&                                  natural;
  const ConformalChart&                         chart;
  const std::vector<occ::handle<Geom2d_Curve>>& held;
  const double                                  area;
  const std::vector<gp_XY>                      loop;
  std::vector<char>                             onBoundary;
  const Inside                                  within;
  int                                           order = 0, degree = 0;
  const occ::handle<Geom_Surface>               source;
  mutable ShapeAnalysis_Surface                 sourceProject;
  const std::vector<Check>                      checks;
  std::vector<double>                           breaksU, breaksV;
};

Built ChartFill::run()
{
  std::optional<Built> best;
  std::vector<Target>  bestTargets;
  int                  idle = 0, unkept = 0, iterations = 0;
  for (double before = HUGE_VAL;;)
  {
    ++iterations;
    SplineNet    net      = flatNet();
    const double smallest = smallestSpan();
    // The Natural fill inside, wished for at the node nearest the middle of every part of
    // every span, `placesPerSpan` parts each way: as many wishes as the net can follow.
    std::vector<Target> targets;
    wishInterior(targets);
    placeSides(targets, smallest);
    placeGuides(targets, smallest);
    occ::handle<Geom_BSplineSurface> surface = fit(
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
    // As in the Natural fit, the fair objective can leave a small positional
    // residual. Use the shared geometric correction rather than relaxing the
    // delivery tolerance. Supported derivatives remain the responsibility of
    // the constrained fit; this correction applies only to positional inputs.
    correctPositions(net, surface);
    auto        checked = certify(surface);
    Built&      built   = checked.built;
    const auto& folded  = checked.folded;
    const bool  folds   = checked.folds;
    // The fill's shape is wished for, not asked: a net that meets everything asked is kept
    // and is made finer only where the shape is followed by less than the tolerance, so
    // that the chart is followed and U and V come out square. That is never a refusal.
    // It is asked no more than a fair surface over its net comes to a curve (see
    // fitReach): a finer net beyond that costs more than the flow gains.
    const double wishedWithin =
      asked.controls.open ? std::max(asked.limits.g0, fitReach * asked.size) : asked.limits.g0;
    bool followed = true;
    for (const Target& target : targets)
      followed =
        followed
        && (!target.soft
            || surface->Value(target.uv.X(), target.uv.Y()).Distance(target.point) <= wishedWithin);
    built.shapeTarget = {followed ? 1.0 : 0.0, static_cast<double>(iterations)};
    if (built.measured.within)
    {
      best        = built;
      bestTargets = targets;
      if (followed)
        break;
    }
    else if (!best || (!best->measured.within && built.measured.excess < best->measured.excess))
    {
      best = built;
    }
    // Nets that come no nearer, more than the coarsest few may: what is left is not the net's.
    if (built.measured.within)
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
    if (!refine(surface, targets, built, folded, folds, wishedWithin, smallest, unkept))
      break;
  }
  compactBest(best, bestTargets);
  best->shapeTarget[1] = static_cast<double>(iterations);
  return *best;
}

std::function<gp_XY(double)> ChartFill::layout(size_t s) const
{
  return std::function<gp_XY(double)>(
    [&, s](double t) { return chart.map(held[s]->Value(t).XY()); });
}

gp_XY ChartFill::interiorPlace(const gp_Pnt& point) const
{
  return chart.map(sourceProject.ValueOfUV(point, BRepExt::tolerance::exactLength).XY());
}

std::vector<double> ChartFill::knots(const std::vector<double>& breaks) const
{
  std::vector<double> out(degree, breaks.front());
  out.insert(out.end(), breaks.begin(), breaks.end());
  out.insert(out.end(), degree, breaks.back());
  return out;
}

ChartFill::Checked ChartFill::certify(const occ::handle<Geom_BSplineSurface>& surface,
                                      bool                                    spatial) const
{
  std::vector<Rim>   rims;
  std::vector<gp_XY> around;
  for (size_t s = 0; s < asked.sides.size(); ++s)
    rims.push_back(
      asked.sides[s].free
        ? Rim{}
        : preimage(surface, asked.sides[s], frame, asked.limits.g0, asked.size, {}, layout(s)));
  closeGaps(asked.sides, rims);
  for (size_t s = 0; s < asked.sides.size(); ++s)
    for (int k = 0; k < 16; ++k)
      around.push_back(rims[s].onSurface->Value(asked.sides[s].parameter(k / 16.0)).XY());
  Built built{3,
              {},
              measureAll(surface,
                         asked,
                         [&](const Side& side, double f, const gp_Pnt&) -> std::optional<Place> {
                           if (const auto s = asked.boundaryIndex(side))
                             return Place{rims[*s].onSurface->Value(side.parameter(f)).XY(), true};
                           // Guides are measured independently by spatial projection.
                           return std::nullopt;
                         })};
  for (size_t s = 0; s < rims.size(); ++s)
    rims[s].deviation =
      std::max(rims[s].deviation, built.measured.values[asked.sides[s].input * 4]);
  built.face = trimmed(surface, asked, rims, enclosed(around) < 0);
  if (spatial || !natural.spatialBoundary.empty())
    for (size_t i = 0; i < rims.size(); ++i)
      if (!asked.sides[i].free)
        built.spatialBoundary.push_back({asked.sides[i].input,
                                         rims[i].onSurface,
                                         asked.sides[i].from,
                                         asked.sides[i].to,
                                         asked.sides[i].originalReversed});
  // Inside the loop: no fold, and how square U and V are.
  std::vector<double> angles;
  double              least = HUGE_VAL, most = 0;
  std::vector<char>   folded(checks.size(), 0);
  bool                folds = false;
  for (size_t i = 0; i < checks.size(); ++i)
  {
    gp_Pnt p;
    gp_Vec su, sv;
    surface->D1(checks[i].flat.X(), checks[i].flat.Y(), p, su, sv);
    const double scale = su.Magnitude() * sv.Magnitude();
    folded[i]          = !(scale > 0) || su.Crossed(sv).Dot(checks[i].normal) <= 0;
    folds              = folds || folded[i];
    if (folded[i])
      continue;
    angles.push_back(std::asin(std::min(1.0, std::abs(su.Dot(sv)) / scale)) * 180 / M_PI);
    least = std::min(least, su.Magnitude() / sv.Magnitude());
    most  = std::max(most, su.Magnitude() / sv.Magnitude());
  }
  if (folds)
  {
    built.measured.within = false;
    built.measured.unkept = true;
    built.measured.worst  = "the conformal surface folds over inside the loop";
    built.measured.excess = std::max(built.measured.excess, unkeptExcess);
  }
  else if (!angles.empty())
  {
    double squares = 0;
    for (double angle : angles)
      squares += angle * angle;
    std::sort(angles.begin(), angles.end());
    built.squareness = {std::sqrt(squares / angles.size()),
                        angles[static_cast<size_t>(0.95 * (angles.size() - 1))],
                        angles.back(),
                        least,
                        most};
  }
  if (built.measured.within)
    if (const auto fault = BRepCheck_Analyzer(built.face).IsValid()
                             ? insideFault(built.face, asked.points.empty() ? rimBend(asked) : 0)
                             : std::optional<std::string>("it is not a valid face"))
    {
      built.measured.within = false;
      built.measured.unkept = true;
      built.measured.worst  = *fault;
      built.measured.excess = unkeptExcess;
    }
  return Checked{std::move(built), std::move(folded), folds};
}

SplineNet ChartFill::flatNet() const
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

double ChartFill::smallestSpan() const
{
  double smallest = HUGE_VAL;
  for (const auto* breaks : {&breaksU, &breaksV})
    for (size_t i = 1; i < breaks->size(); ++i)
      smallest = std::min(smallest, (*breaks)[i] - (*breaks)[i - 1]);
  return smallest;
}

void ChartFill::wishInterior(std::vector<Target>& targets) const
{
  const auto part = [&](const std::vector<double>& breaks, double at) {
    const int span = std::clamp(
      static_cast<int>(std::upper_bound(breaks.begin(), breaks.end(), at) - breaks.begin()) - 1,
      0,
      static_cast<int>(breaks.size()) - 2);
    const double share = (at - breaks[span]) / (breaks[span + 1] - breaks[span]);
    const int    sub   = std::clamp(static_cast<int>(share * placesPerSpan), 0, placesPerSpan - 1);
    return std::pair{span * placesPerSpan + sub,
                     breaks[span]
                       + (breaks[span + 1] - breaks[span]) * (sub + 0.5) / placesPerSpan};
  };
  std::unordered_map<long long, std::pair<int, double>> nearest;
  const long long columns = static_cast<long long>(breaksV.size()) * placesPerSpan;
  for (size_t i = 0; i < chart.flat.size(); ++i)
  {
    if (onBoundary[i])
      continue;
    const auto [a, middleU] = part(breaksU, chart.flat[i].X());
    const auto [b, middleV] = part(breaksV, chart.flat[i].Y());
    const double distance   = (chart.flat[i] - gp_XY(middleU, middleV)).SquareModulus();
    auto& kept = nearest.try_emplace(a * columns + b, static_cast<int>(i), HUGE_VAL).first->second;
    if (distance < kept.second)
      kept = {static_cast<int>(i), distance};
  }
  std::vector<int> chosen;
  for (const auto& [cell, kept] : nearest)
    chosen.push_back(kept.first);
  std::sort(chosen.begin(), chosen.end());
  // Each stands for its share of the loop's area, as a length.
  const double share = std::sqrt(area / std::max<size_t>(1, chosen.size())) / asked.size;
  for (int i : chosen)
  {
    targets.push_back({chart.point[i], chart.flat[i], 0, {}, {}, share});
    targets.back().soft = true;
  }
}

void ChartFill::placeSides(std::vector<Target>& targets, double smallest) const
{
  for (size_t s = 0; s < asked.sides.size(); ++s)
  {
    const Side& side = asked.sides[s];
    if (side.free)
      continue;
    const auto place = layout(s);
    const int  count =
      std::clamp(static_cast<int>(std::ceil(placesPerSpan * side.length / smallest)),
                 fewestRimPlaces,
                 mostRimPlaces);
    const double step = 1.0 / count;
    // Next to a gap the face has its say at the end as along the rest (see overPlane).
    const auto quiet = [&](size_t other) {
      const Side& near = asked.sides[other % asked.sides.size()];
      return near.free || near.support.IsSame(side.support);
    };
    const bool lone       = hasGaps(asked);
    const bool quietStart = lone && quiet(s + asked.sides.size() - 1),
               quietEnd   = lone && quiet(s + 1);
    auto add              = [&](double f, bool end, bool watch) {
      if (end && ((f == 0 && quietStart) || (f == 1 && quietEnd)))
      {
        f   = std::clamp(f, cornerInset, 1 - cornerInset);
        end = false;
      }
      const bool smooth = side.order > 0 && !end;
      // Follow the Natural fill at its own edge correspondence. Any
      // positional correction below targets the original spatial curve.
      const gp_XY uv = held[s]->Value(side.parameter(f)).XY();
      Target      target{source->Value(uv.X(), uv.Y()),
                    place(side.parameter(f)),
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
        const gp_XY  tangent = place(side.parameter(f + h)) - place(side.parameter(f - h));
        target.across        = gp_XY(-tangent.Y(), tangent.X()) / tangent.Modulus();
      }
      targets.push_back(target);
    };
    add(0, true, false);
    add(1, true, false);
    for (int k = 0; k < count; ++k)
    {
      add((k + 0.5) * step, false, false);
      if (k + 1 < count)
        add((k + 1) * step, false, true);
    }
  }
}

void ChartFill::placeGuides(std::vector<Target>& targets, double smallest) const
{
  for (const Side& guide : asked.guides)
  {
    const int count =
      std::clamp(static_cast<int>(std::ceil(placesPerSpan * guide.length / smallest)),
                 fewestRimPlaces,
                 mostRimPlaces);
    for (int i = 0; i <= 2 * count; ++i)
    {
      const gp_Pnt point = guide.meet(static_cast<double>(i) / (2 * count));
      targets.push_back(
        {point, interiorPlace(point), 0, {}, {}, guide.length / count / asked.size});
      targets.back().watch = i % 2 != 0;
    }
  }
  for (const gp_Pnt& point : asked.points)
    targets.push_back({point, interiorPlace(point), 0, {}, {}, 1});
}

void ChartFill::correctPositions(SplineNet& net, occ::handle<Geom_BSplineSurface>& surface) const
{
  if (order == 0 && (!asked.guides.empty() || !asked.points.empty()))
  {
    const auto onto = [&](double span) {
      std::vector<RimPlace> places;
      for (size_t s = 0; s < asked.sides.size(); ++s)
      {
        const Side& side = asked.sides[s];
        if (side.free)
          continue;
        const auto place = layout(s);
        const int  count = std::clamp(static_cast<int>(rimPlacesPerSpan * side.length / span),
                                     fewestRimPlaces,
                                     mostRimPlaces);
        for (int i = 0; i <= count; ++i)
        {
          const double f = static_cast<double>(i) / count;
          places.push_back({side.meet(f), place(side.parameter(f)), {}, false});
        }
      }
      for (const Side& guide : asked.guides)
      {
        const int count = std::clamp(static_cast<int>(rimPlacesPerSpan * guide.length / span),
                                     fewestRimPlaces,
                                     mostRimPlaces);
        for (int i = 0; i <= count; ++i)
        {
          const gp_Pnt point = guide.meet(static_cast<double>(i) / count);
          places.push_back({point, interiorPlace(point), {}, false});
        }
      }
      for (const gp_Pnt& point : asked.points)
        places.push_back({point, interiorPlace(point), {}, false});
      return places;
    };
    net = netOf(surface);
    // The correction is held to the net's own scale, as over the mean plane: a coarse
    // net that is far off is made finer by the fit, and not by thousands of knots here.
    correctRim(
      op,
      net,
      onto,
      frame.z.XYZ(),
      rimShare * asked.limits.g0,
      enoughShare * asked.limits.g0,
      asked.controls.open
        ? std::min(mostRimSpans,
                   rimFinerMet * static_cast<int>(std::max(breaksU.size(), breaksV.size()) - 1))
        : mostRimSpans);
    surface = surfaceOf(net);
  }
}

bool ChartFill::refine(const occ::handle<Geom_BSplineSurface>& surface,
                       const std::vector<Target>&              targets,
                       const Built&                            built,
                       const std::vector<char>&                folded,
                       bool                                    folds,
                       double                                  wishedWithin,
                       double                                  smallest,
                       int&                                    unkept)
{
  std::vector<char> splitU(breaksU.size() - 1, 0), splitV(breaksV.size() - 1, 0);
  const auto        mark = [&](const gp_XY& uv) {
    const auto span = [](const std::vector<double>& breaks, double at) {
      return std::clamp(
        static_cast<int>(std::upper_bound(breaks.begin(), breaks.end(), at) - breaks.begin()) - 1,
        0,
        static_cast<int>(breaks.size()) - 2);
    };
    splitU[span(breaksU, uv.X())] = splitV[span(breaksV, uv.Y())] = 1;
  };
  for (const Target& target : targets)
  {
    bool missed = surface->Value(target.uv.X(), target.uv.Y()).Distance(target.point)
                  > (target.soft ? wishedWithin : fitShare * asked.limits.g0);
    if (!built.measured.within && !missed && target.order > 0)
    {
      const auto   mine   = surfaceJet(op, surface, target.uv.X(), target.uv.Y());
      const double cosine = mine.n.Dot(target.jet->n);
      missed = std::acos(std::min(1., std::abs(cosine))) * 180 / M_PI > fitShare * asked.limits.g1;
      if (!missed && target.order > 1)
      {
        const auto [bend, relative] = bendDifference(*target.jet,
                                                     mine,
                                                     target.jet->u / target.jet->u.Magnitude(),
                                                     cosine < 0 ? -1 : 1);
        missed                      = relative > fitShare * asked.limits.g2
                 && bend * asked.size > fitShare * asked.limits.g2 * flatShare;
      }
    }
    if (missed)
      mark(target.uv);
  }
  // The complete certification reads between fit rows as well, with the same
  // easing and flat-curvature alternatives as final acceptance.
  for (size_t s = 0; s < asked.sides.size(); ++s)
  {
    const Side& side = asked.sides[s];
    if (side.free || built.measured.within)
      continue;
    const auto place = layout(s);
    const Miss miss  = measure(surface, side, [&](const Side&, double f, const gp_Pnt&) {
      return std::optional<Place>(Place{place(side.parameter(f)), true});
    });
    if (miss.within(asked.limits, side.order, asked.size))
      continue;
    const int count =
      std::clamp(static_cast<int>(std::ceil(4 * placesPerSpan * side.length / smallest)),
                 fewestRimPlaces,
                 4 * mostRimPlaces);
    for (int i = 0; i <= count; ++i)
      mark(place(side.parameter(static_cast<double>(i) / count)));
  }
  for (size_t i = 0; i < checks.size(); ++i)
    if (folded[i])
      mark(checks[i].flat);
  // A net that meets the edges with a surface that is not one to keep is too coarse for
  // what is asked along them, wherever that shows: a few finer nets, and no more.
  if (built.measured.unkept && !folds && ++unkept <= mostUnkept)
  {
    std::fill(splitU.begin(), splitU.end(), 1);
    std::fill(splitV.begin(), splitV.end(), 1);
  }
  bool finer = false;
  for (auto [breaks, split] : {std::pair{&breaksU, &splitU}, std::pair{&breaksV, &splitV}})
  {
    std::vector<double> next{breaks->front()};
    int                 spans = static_cast<int>(breaks->size()) - 1;
    for (size_t i = 0; i + 1 < breaks->size(); ++i)
    {
      if ((*split)[i] && spans < std::end(ladder)[-1])
      {
        next.push_back(((*breaks)[i] + (*breaks)[i + 1]) / 2);
        ++spans;
        finer = true;
      }
      next.push_back((*breaks)[i + 1]);
    }
    *breaks = next;
  }
  return finer;
}

void ChartFill::compactBest(std::optional<Built>&      best,
                            const std::vector<Target>& bestTargets) const
{
  if (order == 0)
    measureSpatialBoundary(*best, asked, false);
  if (order == 0 && best->measured.within && !best->measured.unkept)
  {
    double gap = 0;
    for (const Side& side : asked.sides)
      gap = std::max(gap, best->measured.values[side.input * 4]);
    const double allowance = std::nextafter(asked.limits.g0 - gap, 0.0);
    const auto   original  = occ::down_cast<Geom_BSplineSurface>(BRep_Tool::Surface(best->face));
    if (std::isfinite(gap) && gap >= 0 && !original.IsNull()
        && std::all_of(best->measured.values.begin(),
                       best->measured.values.end(),
                       [](double value) { return std::isfinite(value); }))
      if (const auto compact = compactPatchSurface(op, original, allowance))
      {
        try
        {
          auto reduced = certify(compact->surface, true);
          measureSpatialBoundary(reduced.built, asked, false);
          bool followed = true;
          for (const Target& target : bestTargets)
          {
            if (!target.soft)
              continue;
            const double before =
              original->Value(target.uv.X(), target.uv.Y()).Distance(target.point);
            const double after =
              compact->surface->Value(target.uv.X(), target.uv.Y()).Distance(target.point);
            followed = followed && std::isfinite(before) && std::isfinite(after)
                       && after <= std::max(asked.limits.g0, before);
          }
          const bool finite = std::all_of(reduced.built.squareness.begin(),
                                          reduced.built.squareness.end(),
                                          [](double value) { return std::isfinite(value); });
          if (reduced.built.measured.within && !reduced.built.measured.unkept && followed && finite)
          {
            reduced.built.shapeTarget = best->shapeTarget;
            best                      = std::move(reduced.built);
          }
        }
        catch (const Failure& failure)
        {
          if (failure.reason != no_solution)
            throw;
        }
        catch (const Standard_Failure&)
        {
          // The certified original remains available if its reduced trim cannot be rebuilt.
        }
      }
  }
}

//! The Adjacent fit over the chart.
Built adjacentOver(const Asked&                                  asked,
                   const Frame&                                  frame,
                   const Built&                                  natural,
                   const ConformalChart&                         chart,
                   const std::vector<occ::handle<Geom2d_Curve>>& held,
                   double                                        area)
{
  ChartFill fill(asked, frame, natural, chart, held, area);
  return fill.run();
}

} // namespace

Built adjacent(const Asked& asked, const Frame& frame, const Built& natural)
{
  PatchPhase phase("adjacent");
  double     perimeter = 0;
  for (const Side& side : asked.sides)
    perimeter += side.length;
  std::vector<gp_XY> outline;
  for (const gp_Pnt& point : polygon(asked.sides))
    outline.push_back(frame.project(point));
  const double area = std::abs(enclosed(outline));
  // Each edge as the Natural face holds it, in that face's parameters, and so in the chart.
  std::vector<occ::handle<Geom2d_Curve>> held;
  for (const Side& side : asked.sides)
  {
    if (side.free)
    {
      held.push_back(natural.gaps[&side - asked.sides.data()]);
      continue;
    }
    occ::handle<Geom2d_Curve> found;
    for (TopExp_Explorer it(natural.face, TopAbs_EDGE); it.More() && found.IsNull(); it.Next())
    {
      double          first, last;
      TopLoc_Location location;
      const auto      edge = TopoDS::Edge(it.Current());
      if (BRep_Tool::Curve(edge, location, first, last) == side.curve && location.IsIdentity()
          && first == std::min(side.from, side.to) && last == std::max(side.from, side.to))
        found = BRep_Tool::CurveOnSurface(edge, natural.face, first, last);
    }
    if (found.IsNull())
      fail(op,
           no_solution,
           "the Natural fill does not hold boundary edge " + std::to_string(side.input + 1));
    held.push_back(found);
  }
  const ConformalChart chart = conformalChart(natural.face,
                                              area / perimeter / trianglesAcross,
                                              chartDeflection * asked.limits.g0);
  return adjacentOver(asked, frame, natural, chart, held, area);
}

} // namespace BRepExt::patch_fill
