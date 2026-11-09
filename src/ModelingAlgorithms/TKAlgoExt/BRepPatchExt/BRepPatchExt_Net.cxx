// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// The nets a fill is fitted on: the conditions a fit is held to, the fit itself, and the knots of
// its net.

#include <BRepPatchExt_Fill.hxx>

#include <math_Gauss.hxx>
#include <math_Matrix.hxx>
#include <math_Vector.hxx>

#include <array>

namespace BRepExt::patch_fill
{
namespace
{

//! How near none the derivative across a side along the axis of the flow asked must come for
//! the condition to count as met. Chosen.
constexpr double flowWithin = 1e-7;

//! Where curves are only met approximately, a quadratic is near enough to start from where it
//! misses them by less than this share of the loop's size. Chosen.
constexpr double approximateShare = 1e-3;

//! Two breaks of a net closer than this share of an even span count as one. Chosen.
constexpr double sameBreak = 1e-6;

//! What the change of the bend counts for against the bend, with the loop's size as the unit
//! of length: in full where a bend is asked or wished for along the rim, and by these shares
//! where the rim is only tangent to its faces, and where curves are only met. With more of
//! it the surface is pulled flat between its edges; with none, a point inside the loop stands
//! in a spike. Measured for the tangent rim on the corner of three rounds of 4, whose own
//! surface is a ball, as the distance from that ball at tolerances of 1e-3 and 1e-6: with
//! 0.03, 0.066 and 0.065; with 0.003, 0.024 and 0.068; with none, 0.010 and 0.069. For curves
//! that are only met, 0.03 is the first author's, from a round hole in a sphere; with 0.003
//! an L of six sides on a curved height is missed by 2.6e-6 on every net where 0.03 meets a
//! millionth on the fourth.
constexpr double variationTangent = 0.003, variationMet = 0.03;

//! The fits of a net whose bend is asked along an edge, the first of them without the bend:
//! what is asked is taken about the surface before, and over a steep rim (a ball cut at a
//! quarter of its radius) the third still left the bend 6 to 8 % off where the fifth meets it.
constexpr int askedPasses = 5;

//! The same where the bend is only wished for: nearer than the third is not needed, and on a
//! ball's cap the third is nearer the ball than the fifth.
constexpr int wishedPasses = 3;

//! How often the places that a fit misses between those it was asked are added to them.
constexpr int placeRounds = 3;

//! The speed a surface over a chart is to leave its rim at, across it, and the band about it
//! that is met: over a chart that keeps the rim's lengths (proxyChart), a surface that the chart
//! is conformal to leaves the rim as fast as it runs along it, at one; the band keeps it from
//! running off to nothing or far off, which meet every tangent plane and are no parametrisation.
//! The band is a judgement, not derived: the chart is the proxy's, not the fill's.
constexpr double leaveSpeed = 1, leaveBand = 0.5;

//! The conditions that say so. For a height graph in the normalized mean-plane frame,
//! tangency determines its spatial derivatives exactly, making G2 linear as well. Otherwise
//! G2 is taken about `current`, and left out while there is none.
std::vector<FitCondition> conditions(const std::vector<Target>&              targets,
                                     const occ::handle<Geom_BSplineSurface>& current,
                                     double                                  size,
                                     const Limits&                           limits,
                                     std::vector<int>&                       owner,
                                     GraphFit                                graph)
{
  std::vector<FitCondition> out;
  owner.clear();
  for (const Target& target : targets)
  {
    // Whose every condition is: all that this place adds below.
    struct Own
    {
      std::vector<int>&                owner;
      const std::vector<FitCondition>& out;
      int                              target;

      ~Own() { owner.resize(out.size(), target); }
    } own{owner, out, static_cast<int>(&target - targets.data())};

    const double share = target.watch ? watchShare : fitShare;
    const double near = share * limits.g0, flat = share * std::sin(limits.g1 * M_PI / 180);
    for (int c = 1; c <= 3; ++c)
    {
      gp_XYZ axis(0, 0, 0);
      axis.SetCoord(c, 1);
      out.push_back({target.uv.X(),
                     target.uv.Y(),
                     {{0, 0, axis}},
                     1,
                     axis.Dot(target.point.XYZ()),
                     target.share / (size * size),
                     target.soft,
                     near,
                     target.watch});
    }
    if (target.perpendicular)
    {
      const gp_XYZ direction = target.perpendicular->Normalized().XYZ();
      out.push_back({target.uv.X(),
                     target.uv.Y(),
                     {{1, 0, direction * target.across.X()}, {0, 1, direction * target.across.Y()}},
                     2,
                     0,
                     target.share * tangentWeight,
                     false,
                     flowWithin});
    }
    if (target.leave)
    {
      const gp_XYZ direction = target.leave->XYZ();
      out.push_back({target.uv.X(),
                     target.uv.Y(),
                     {{1, 0, direction * target.across.X()}, {0, 1, direction * target.across.Y()}},
                     2,
                     leaveSpeed,
                     target.share * tangentWeight,
                     false,
                     leaveBand,
                     !target.leaveAsked});
    }
    if (target.order < 1)
      continue;
    const SurfaceJet& jet = *target.jet;
    const gp_XYZ      n   = jet.n.XYZ();
    const double      du = target.across.X(), dv = target.across.Y();
    // One tangent plane with the neighbour: both directions of the surface in it. Along
    // the edge that follows from holding the edge, and is asked all the same: measured on
    // a hole across a round, the direction across alone left the fit short of a millionth
    // on every net, and both met it on the second.
    out.push_back({target.uv.X(),
                   target.uv.Y(),
                   {{1, 0, n}},
                   1,
                   0,
                   target.share * tangentWeight,
                   false,
                   flat,
                   target.watch});
    out.push_back({target.uv.X(),
                   target.uv.Y(),
                   {{0, 1, n}},
                   1,
                   0,
                   target.share * tangentWeight,
                   false,
                   flat,
                   target.watch});
    const bool wished = target.order < 2 && target.wish;
    if ((target.order < 2 && !wished)
        || (graph == GraphFit::none && (current.IsNull() || target.watch)))
      continue;
    if (graph != GraphFit::none)
    {
      // The support normal fixes the graph's first derivatives, so evaluating
      // its second fundamental form is exact and linear before fitting. Start
      // with transverse curvature; if measured G2 still fails, the complete
      // form also constrains tangent/twist errors from approximate rim fits.
      const double               diagonal = std::sqrt(0.5);
      const std::array<gp_XY, 3> directions =
        graph == GraphFit::complete ? std::array<gp_XY, 3>{{{1, 0}, {0, 1}, {diagonal, diagonal}}}
                                    : std::array<gp_XY, 3>{{{du, dv}, {}, {}}};
      for (int i = 0; i < (graph == GraphFit::complete ? 3 : 1); ++i)
      {
        const double u = directions[i].X(), v = directions[i].Y();
        const gp_Vec tangent(u, v, -(n.X() * u + n.Y() * v) / n.Z());
        const double bend = jet.curvature(tangent, tangent);
        out.push_back({target.uv.X(),
                       target.uv.Y(),
                       {{2, 0, n * (u * u)}, {1, 1, n * (2 * u * v)}, {0, 2, n * (v * v)}},
                       3,
                       bend,
                       target.share * size * size * bendWeight,
                       wished,
                       share * limits.g2
                         * std::max({std::abs(bend),
                                     target.easing * tangent.SquareMagnitude(),
                                     tangent.SquareMagnitude() * flatShare / size}),
                       target.watch});
      }
      continue;
    }
    // With T the way across in space, the bend across is n . S_dd, to equal II(T, T) of the
    // neighbour; about the T of `current` that is 2 II(T0, T) - II(T0, T0).
    gp_Pnt point;
    gp_Vec su, sv;
    current->D1(target.uv.X(), target.uv.Y(), point, su, sv);
    const gp_Vec t0 = su * du + sv * dv;
    const gp_Vec e1 = jet.u.Normalized(), e2 = jet.n.Crossed(e1);
    const gp_XYZ w = (e1 * jet.curvature(t0, e1) + e2 * jet.curvature(t0, e2)).XYZ();
    out.push_back({target.uv.X(),
                   target.uv.Y(),
                   {{2, 0, n * (du * du)},
                    {1, 1, n * (2 * du * dv)},
                    {0, 2, n * (dv * dv)},
                    {1, 0, w * (-2 * du)},
                    {0, 1, w * (-2 * dv)}},
                   5,
                   -jet.curvature(t0, t0),
                   target.share * size * size * bendWeight,
                   wished,
                   fitShare * limits.g2 * t0.SquareMagnitude()
                     * std::max({std::abs(jet.curvature(t0, t0)) / t0.SquareMagnitude(),
                                 target.easing,
                                 flatShare / size})});
  }
  return out;
}

//! The most spans a net is given across one direction, however much is asked to be finer.
constexpr int mostSpans = 72;

} // namespace

occ::handle<Geom_BSplineSurface> fit(SplineNet&                                   net,
                                     const std::vector<char>&                     fixed,
                                     std::vector<Target>                          targets,
                                     double                                       size,
                                     const Limits&                                limits,
                                     occ::handle<Geom_BSplineSurface>             surface,
                                     const std::function<double(double, double)>& density,
                                     const gp_XYZ*                                along,
                                     const std::vector<FitJump>&                  jumps,
                                     const Frame*                                 frame,
                                     double                                       flatness,
                                     bool                                         leastSquares,
                                     GraphFit                                     graph)
{
  if (frame)
  {
    // Solve in the loop's orthonormal frame, including the unknowns and the energy's
    // right hand side. Stopping on world-coordinate components makes a rotated request
    // stop on another pass; subtracting large world positions also loses the small miss.
    const auto vector = [&](const gp_Vec& v) {
      return gp_Vec(v.Dot(frame->x), v.Dot(frame->y), v.Dot(frame->z));
    };
    // Coordinates in units of the size are rounded to a grid: what rotating them into the frame
    // rounds differently (a few ulps) then rounds to the same point, except within that much of
    // the grid's halfway points. The grid is of the size, never of the model's units, and finer
    // than a tenth of the relative tolerance.
    const double quantum    = std::min(1 / frameSteps, 0.1 * limits.g0 / size);
    const auto   coordinate = [&](double v) { return std::round(v / (size * quantum)) * quantum; };
    const auto   scaled     = [&](const gp_Vec& v) {
      const gp_Vec q = vector(v);
      return gp_Vec(coordinate(q.X()), coordinate(q.Y()), coordinate(q.Z()));
    };
    for (gp_XYZ& pole : net.poles)
      pole = scaled(gp_Vec(frame->origin, gp_Pnt(pole))).XYZ();
    for (auto* knots : {&net.knotsU, &net.knotsV})
      for (double& knot : *knots)
        knot = coordinate(knot);
    for (Target& target : targets)
    {
      if (target.perpendicular)
        target.perpendicular = vector(*target.perpendicular);
      if (target.leave)
        target.leave = vector(*target.leave);
      target.point = gp_Pnt(scaled(gp_Vec(frame->origin, target.point)).XYZ());
      target.uv.SetCoord(coordinate(target.uv.X()), coordinate(target.uv.Y()));
      // A bend, as the jets' are, in the units of the size.
      target.easing *= size;
      if (target.jet)
      {
        auto& j = *target.jet;
        j.p     = gp_Pnt(scaled(gp_Vec(frame->origin, j.p)).XYZ());
        for (gp_Vec* v : {&j.u, &j.v, &j.uu, &j.vv, &j.uv})
          *v = scaled(*v);
        j.n = vector(j.n);
        for (int c = 1; c <= 3; ++c)
          j.n.SetCoord(c, std::round(j.n.Coord(c) * frameSteps) / frameSteps);
        j.n.Normalize();
      }
    }
    std::vector<FitJump> localJumps = jumps;
    for (auto& jump : localJumps)
    {
      jump.at    = coordinate(jump.at);
      jump.along = coordinate(jump.along);
      jump.weight /= size * size;
    }
    const gp_XYZ localAlong = along ? vector(gp_Vec(*along)).XYZ() : gp_XYZ(0, 0, 0);
    std::function<double(double, double)> localDensity;
    if (density)
      localDensity = [&](double u, double v) { return density(u * size, v * size); };
    fit(net,
        fixed,
        std::move(targets),
        1,
        {limits.g0 / size, limits.g1, limits.g2},
        surface.IsNull() ? occ::handle<Geom_BSplineSurface>{} : surfaceOf(net),
        localDensity,
        along ? &localAlong : nullptr,
        localJumps,
        nullptr,
        flatness,
        leastSquares,
        graph);
    for (auto* knots : {&net.knotsU, &net.knotsV})
      for (double& knot : *knots)
        knot *= size;
    for (gp_XYZ& pole : net.poles)
      pole = frame->origin.XYZ()
             + (frame->x.XYZ() * pole.X() + frame->y.XYZ() * pole.Y() + frame->z.XYZ() * pole.Z())
                 * size;
    return surfaceOf(net);
  }
  const bool      bends   = std::any_of(targets.begin(), targets.end(), [](const Target& target) {
    return target.order > 1 || target.wish;
  });
  const bool      tangent = std::any_of(targets.begin(), targets.end(), [](const Target& target) {
    return target.order > 0;
  });
  const FitEnergy energy{(flatness < 0 ? 1 : 1 - flatness) * bendingWeight,
                         (flatness < 0 ? (bends     ? 1
                                          : tangent ? variationTangent
                                                    : variationMet)
                                       : flatness)
                           * bendingWeight * size * size};
  const bool      asked = std::any_of(targets.begin(), targets.end(), [](const Target& target) {
    return target.order > 1;
  });
  // Only free spatial fits need to iterate their nonlinear bend conditions. A graph's
  // tangents are fixed by the support normal, so its first solve already asks for G2.
  for (int pass = 0; pass < (graph != GraphFit::none ? 1
                             : asked                 ? askedPasses
                             : bends                 ? wishedPasses
                                                     : 1);
       ++pass)
  {
    // A place between those asked that the surface misses on its own is asked as well,
    // and the net fitted again from where it stood: places are added where the energy
    // alone does not keep the surface on a curve, and nowhere else, since every place
    // asked is one more thing that the surface has to follow whatever it costs it inside.
    const std::vector<gp_XYZ> start = net.poles;
    FitEnergies               energies;
    for (int round = 0;; ++round)
    {
      std::vector<int>       owner;
      const std::vector<int> missed =
        fitNet(op,
               net,
               fixed,
               conditions(targets, surface, size, limits, owner, graph),
               energy,
               density,
               along,
               jumps,
               leastSquares,
               &energies);
      if (missed.empty() || round == placeRounds)
        break;
      for (int condition : missed)
      {
        targets[owner[condition]].watch      = false;
        targets[owner[condition]].leaveAsked = true;
      }
      net.poles = start;
    }
    surface = surfaceOf(net);
  }
  return surface;
}

std::vector<double> clampedKnots(double from, double to, int spans, int degree)
{
  std::vector<double> knots(degree + 1, from);
  for (int i = 1; i < spans; ++i)
    knots.push_back(from + (to - from) * i / spans);
  knots.insert(knots.end(), degree + 1, to);
  return knots;
}

bool polynomialHeight(const Asked& asked, const Frame& frame, SplineNet& net, bool approximate)
{
  if (asked.flatness >= 0)
    return false;
  const bool interior = !asked.guides.empty() || !asked.points.empty();
  const auto terms    = [&](const gp_Pnt& point, int degree) {
    const gp_XY         uv = frame.project(point) / asked.size;
    std::vector<double> out{1};
    gp_XY               power(1, 0);
    for (int k = 1; k <= degree; ++k)
    {
      power.SetCoord(power.X() * uv.X() - power.Y() * uv.Y(),
                     power.X() * uv.Y() + power.Y() * uv.X());
      out.push_back(power.X());
      out.push_back(power.Y());
    }
    if (interior)
      out.push_back(uv.SquareModulus());
    return out;
  };
  const auto choose = [](int n, int k) {
    double value = 1;
    for (int j = 1; j <= k; ++j)
      value *= static_cast<double>(n - j + 1) / j;
    return value;
  };
  const auto powers = [&](const std::vector<double>& knots, int degree, int i) {
    std::array<double, 6> out{1, 0, 0, 0, 0, 0};
    for (int j = 1; j <= degree; ++j)
      for (int k = j; k > 0; --k)
        out[k] += out[k - 1] * knots[i + j] / asked.size;
    for (int k = 1; k <= degree; ++k)
      out[k] /= choose(degree, k);
    return out;
  };
  for (int degree = 2; degree <= (interior ? 2 : std::min(net.degreeU, net.degreeV)); ++degree)
  {
    const int   count = 1 + 2 * degree + static_cast<int>(interior);
    math_Matrix matrix(1, count, 1, count, 0.0);
    math_Vector value(1, count, 0.0), coefficient(1, count);
    const auto  add = [&](const gp_Pnt& at, double weight) {
      const auto term = terms(at, degree);
      for (int a = 0; a < count; ++a)
      {
        value(a + 1) += weight * term[a] * frame.height(at);
        for (int b = 0; b < count; ++b)
          matrix(a + 1, b + 1) += weight * term[a] * term[b];
      }
    };
    for (const auto* curves : {&asked.sides, &asked.guides})
      for (const Side& side : *curves)
        for (int i = 0; i < places; ++i)
          add(side.meet((i + .5) / places), side.length / places);
    for (const gp_Pnt& at : asked.points)
      add(at, asked.size);
    math_Gauss solve(matrix);
    if (!solve.IsDone())
      continue;
    solve.Solve(value, coefficient);
    double     miss    = 0;
    const auto measure = [&](const gp_Pnt& at) {
      const auto term   = terms(at, degree);
      double     height = 0;
      for (int k = 0; k < count; ++k)
        height += coefficient(k + 1) * term[k];
      miss = std::max(miss, std::abs(height - frame.height(at)));
    };
    for (const auto* curves : {&asked.sides, &asked.guides})
      for (const Side& side : *curves)
        for (int i = 0; i <= places; ++i)
          measure(side.meet(static_cast<double>(i) / places));
    for (const gp_Pnt& at : asked.points)
      measure(at);
    // An interior constraint must identify the quadratic within the actual tolerance.
    // Approximate harmonic free rims are corrected afterwards, as before.
    if (miss >= (interior && !approximate ? .1 * asked.limits.g0 : approximateShare * asked.size))
      continue;
    for (int i = 0; i < net.countU(); ++i)
      for (int j = 0; j < net.countV(); ++j)
      {
        const auto u = powers(net.knotsU, net.degreeU, i), v = powers(net.knotsV, net.degreeV, j);
        double     height = coefficient(1);
        for (int m = 1; m <= degree; ++m)
          for (int k = 0; k <= m; ++k)
          {
            const double term = choose(m, k) * u[m - k] * v[k], sign = k % 4 < 2 ? 1 : -1;
            height += coefficient(2 * m + k % 2) * sign * term;
          }
        if (interior)
          height += coefficient(count) * (u[2] + v[2]);
        net.poles[static_cast<size_t>(i) * net.countV() + j] += frame.z.XYZ() * height;
      }
    return true;
  }
  return false;
}

std::vector<double> gradedKnots(double                     from,
                                double                     to,
                                int                        spans,
                                int                        degree,
                                const std::vector<Finer>&  finer,
                                const std::vector<double>& breaks)
{
  // What is wanted is summed over this many cells, which is as finely as a knot is placed.
  constexpr int       cells = 4096;
  constexpr double    snap  = 0.35;
  std::vector<double> wanted(cells, spans / (to - from));
  for (const Finer& place : finer)
  {
    const int first =
                std::clamp(static_cast<int>((place.at - place.span - from) / (to - from) * cells),
                           0,
                           cells - 1),
              last =
                std::clamp(static_cast<int>((place.at + place.span - from) / (to - from) * cells),
                           0,
                           cells - 1);
    for (int i = first; i <= last; ++i)
      wanted[i] = std::max(wanted[i], 1 / place.span);
  }
  double total = 0;
  for (double density : wanted)
    total += density * (to - from) / cells;
  const int count =
    std::clamp(static_cast<int>(std::lround(total)), spans, std::max(spans, mostSpans));
  std::vector<double> even{from};
  double              sum = 0;
  for (int i = 0; i < cells && static_cast<int>(even.size()) < count; ++i)
  {
    const double step = wanted[i] * (to - from) / cells * count / total;
    while (static_cast<int>(even.size()) < count && sum + step >= static_cast<double>(even.size()))
      even.push_back(from
                     + (to - from) * (i + (static_cast<double>(even.size()) - sum) / step) / cells);
    sum += step;
  }
  even.push_back(to);
  std::vector<double> at = breaks;
  std::sort(at.begin(), at.end());
  std::vector<std::pair<double, int>> inner;
  for (double value : at)
    if (value > from && value < to
        && (inner.empty() || value - inner.back().first > sameBreak * (to - from) / count))
      inner.emplace_back(value, degree - 1);
  const size_t broken = inner.size();
  for (size_t i = 1; i + 1 < even.size(); ++i)
  {
    const double span  = std::min(even[i] - even[i - 1], even[i + 1] - even[i]);
    bool         taken = false;
    for (size_t k = 0; k < broken; ++k)
      taken = taken || std::abs(inner[k].first - even[i]) < snap * span;
    if (!taken)
      inner.emplace_back(even[i], 1);
  }
  std::sort(inner.begin(), inner.end());
  std::vector<double> knots(degree + 1, from);
  for (const auto& [knot, times] : inner)
    knots.insert(knots.end(), times, knot);
  knots.insert(knots.end(), degree + 1, to);
  return knots;
}

} // namespace BRepExt::patch_fill
