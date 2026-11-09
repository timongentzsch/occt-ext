// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// How a fill is measured: its misses along each edge and at the corners, the way Measure Continuity
// reads them, inside its loop, and once its face is sewn.

#include <BRepPatchExt_Fill.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>

#include <BRepToolsExt_GeomUtil.hxx>
#include <GeomAlgoExt_Cancel.hxx>
#include <GeomAlgoExt_PatchRim.hxx>
#include <GeomAlgoExt_Phase.hxx>

#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomConvert.hxx>
#include <GeomConvert_ApproxCurve.hxx>
#include <GeomConvert_BSplineCurveToBezierCurve.hxx>
#include <GeomConvert_BSplineSurfaceToBezierSurface.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_BezierSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <math_Gauss.hxx>
#include <math_Matrix.hxx>
#include <math_Vector.hxx>

#include <array>
#include <cstdio>
#include <map>
#include <memory>
#include <tuple>

namespace BRepExt::patch_fill
{
namespace
{

//! The parameters of the surface's point nearest `point`, looked for from `uv` within the
//! surface's parameter bounds: steps of Gauss and Newton, and along the border once a step
//! would leave it, until the step falls to the rounding of the parameters. A coordinate held at
//! the border is let go again when the free step leads back inside. Nothing where the steps do
//! not settle; the global search is then the one to ask.
std::optional<gp_XY> nearestWithin(const occ::handle<Geom_Surface>& surface,
                                   gp_XY                            uv,
                                   const gp_Pnt&                    point)
{
  double low[2], high[2];
  surface->Bounds(low[0], high[0], low[1], high[1]);
  const bool periodic[2] = {surface->IsUPeriodic(), surface->IsVPeriodic()};
  for (int step = 0; step < 32; ++step)
  {
    gp_Pnt at;
    gp_Vec su, sv;
    surface->D1(uv.X(), uv.Y(), at, su, sv);
    const gp_Vec to(at, point);
    const double E = su.Dot(su), F = su.Dot(sv), G = sv.Dot(sv), det = E * G - F * F;
    if (!(E > 0) || !(G > 0) || !(det > 1e-12 * E * G))
      return std::nullopt;
    double move[2] = {(G * to.Dot(su) - F * to.Dot(sv)) / det,
                      (E * to.Dot(sv) - F * to.Dot(su)) / det};
    // A coordinate on its border whose free step leaves it is held there, and the other
    // takes the step along that border alone.
    bool held[2] = {false, false};
    for (int c = 0; c < 2; ++c)
    {
      const double here = c ? uv.Y() : uv.X();
      held[c]           = !periodic[c]
                && ((here <= low[c] && move[c] < 0) || (here >= high[c] && move[c] > 0)
                    || here + move[c] < low[c] || here + move[c] > high[c]);
    }
    if (held[0] && !held[1])
      move[1] = to.Dot(sv) / G;
    if (held[1] && !held[0])
      move[0] = to.Dot(su) / E;
    double next[2] = {uv.X() + (held[0] ? 0 : move[0]), uv.Y() + (held[1] ? 0 : move[1])};
    for (int c = 0; c < 2; ++c)
    {
      if (held[c])
        next[c] = std::clamp(c ? uv.Y() + move[1] : uv.X() + move[0], low[c], high[c]);
      else if (!periodic[c])
        next[c] = std::clamp(next[c], low[c], high[c]);
    }
    const double du = next[0] - uv.X(), dv = next[1] - uv.Y();
    uv.SetCoord(next[0], next[1]);
    if (du * du * E + dv * dv * G <= 1e-28 * (E + G) * (1 + uv.SquareModulus()))
      return uv;
  }
  return std::nullopt;
}

//! What rounding makes of a curvature read from a surface's jet, as a share of its largest
//! principal curvature: a curvature that is nothing comes out as up to about this much of the ones
//! beside it. A judgement, not derived: a few thousand units in the last place of a double (2^-52)
//! for the sums of the jet's derivatives, far below any curvature asked and far above rounding.
constexpr double jetRounding = 0x1p-40;

//! What the faces on the two sides of a corner say of the one surface that is to leave both
//! tangentially. Along its edge a face gives the surface its bend along the edge and its twist
//! (the turn of the tangent plane along the edge); at a corner that is four numbers for the
//! three of one bend, and where the four do not fit, the surface has to change its bend across a
//! line through the corner. `conflict` is by what share of their curvature they do not: zero
//! where they fit, as the two rounds of a ball's corner do, and a plane and a round tangent to it
//! do not. Where everything the four numbers ask is no more than rounding makes of the faces'
//! curvatures (`jetRounding` of the larger of the two faces' own principal curvatures, uneased),
//! it is none: what is asked is rounding, and a share of rounding over rounding says nothing; it
//! was read against an absolute 1e-15, above which rounding came at one size and not another. With
//! `bends`, each face's bend across its edge is asked as well, counting for as much less as a bend
//! counts for against a tangent plane (see tangentWeight), so that it decides the one bend only
//! where the four numbers leave it open (an edge that runs straight on); `across` is then by how
//! much each face's bend across its edge differs from that one bend, a's and b's, each in its own
//! face's sense.
struct CornerBend
{
  double conflict  = 0;
  double across[2] = {0, 0};
};

CornerBend cornerBend(const Side& a, const Side& b, bool bends)
{
  const SurfaceJet ja = a.jet(1 - cornerInset), jb = b.jet(cornerInset);
  const double     sign = ja.n.Dot(jb.n) < 0 ? -1 : 1;
  const gp_Vec     ta = a.tangent(1 - cornerInset).Normalized(), na = ja.n.Crossed(ta);
  gp_Vec           tb = b.tangent(cornerInset);
  double           c = tb.Dot(ta), s = tb.Dot(na);
  const double     length = std::hypot(c, s);
  if (length < 1e-12)
    return {};
  c /= length;
  s /= length;
  tb              = ta * c + na * s;
  const gp_Vec nb = na * c - ta * s;
  // The bend as L, M, N over (ta, na).
  const double rows[6][3] = {{1, 0, 0},
                             {0, 1, 0},
                             {c * c, 2 * c * s, s * s},
                             {-c * s, c * c - s * s, c * s},
                             {0, 0, 1},
                             {s * s, -2 * c * s, c * c}};
  const double asked[6]   = {ja.curvature(ta, ta),
                             ja.curvature(ta, na),
                             sign * jb.curvature(tb, tb),
                             sign * jb.curvature(tb, nb),
                             ja.curvature(na, na),
                             sign * jb.curvature(nb, nb)};
  const double weight     = bendWeight / tangentWeight;
  math_Matrix  normal(1, 3, 1, 3, 0.0);
  math_Vector  right(1, 3, 0.0), bend(1, 3, 0.0);
  for (int r = 0; r < (bends ? 6 : 4); ++r)
  {
    const double w = r < 4 ? 1 : weight * weight;
    for (int i = 0; i < 3; ++i)
    {
      right(i + 1) += w * rows[r][i] * asked[r];
      for (int k = 0; k < 3; ++k)
        normal(i + 1, k + 1) += w * rows[r][i] * rows[r][k];
    }
  }
  // Where the edge runs straight on, nothing says how the surface bends across it.
  for (int i = 1; i <= 3; ++i)
    normal(i, i) += 1e-9;
  math_Gauss solve(normal);
  if (!solve.IsDone())
    return {};
  solve.Solve(right, bend);
  CornerBend out;
  double     scale = 0;
  for (int r = 0; r < 4; ++r)
  {
    out.conflict = std::max(
      out.conflict,
      std::abs(rows[r][0] * bend(1) + rows[r][1] * bend(2) + rows[r][2] * bend(3) - asked[r]));
    scale = std::max(scale, std::abs(asked[r]));
  }
  // Rounding comes from the faces' own derivatives, whatever an ease has taken off them.
  const double rounding =
    jetRounding * std::max(sharpest(a.faceJet(1 - cornerInset)), sharpest(b.faceJet(cornerInset)));
  out.conflict = scale > rounding ? out.conflict / scale : 0;
  for (int side = 0; bends && side < 2; ++side)
  {
    const double* row = rows[4 + side];
    out.across[side]  = (row[0] * bend(1) + row[1] * bend(2) + row[2] * bend(3) - asked[4 + side])
                       * (side ? sign : 1);
  }
  return out;
}

} // namespace

std::string number(double value)
{
  char text[32];
  std::snprintf(text, sizeof text, "%.3g", value);
  return text;
}

const std::vector<double>& measuredFractions()
{
  static const std::vector<double> out = [] {
    std::vector<double> toEnd, all;
    for (double t = 0.25 / places; t >= cornerInset; t /= 2)
      toEnd.push_back(t);
    all.insert(all.end(), toEnd.rbegin(), toEnd.rend());
    for (int i = 0; i < places; ++i)
      all.push_back((i + 0.5) / places);
    for (const double t : toEnd)
      all.push_back(1 - t);
    return all;
  }();
  return out;
}

Miss measure(const occ::handle<Geom_Surface>& surface, const Side& side, const Known& known)
{
  ShapeAnalysis_Surface      project(surface);
  Miss                       miss;
  gp_Pnt2d                   uv;
  const std::vector<double>& reads = measuredFractions();
  const int                  count = static_cast<int>(reads.size());
  for (int i = 0; i < count + 2; ++i)
  {
    const bool   end      = i >= count;
    const double fraction = end ? i - count : reads[i];
    const gp_Pnt point    = side.meet(fraction);
    const auto   start    = known ? known(side, fraction, point) : std::nullopt;
    if (start && start->own)
    {
      uv = gp_Pnt2d(start->uv);
    }
    else
    {
      // From the place before, or the one known, within the surface's border: OpenCascade's
      // next step gives up on the border, where a side that is the surface's own edge has
      // all its places, and searches the whole surface again for each.
      const auto local =
        start || i > 0 ? nearestWithin(surface, start ? start->uv : uv.XY(), point) : std::nullopt;
      uv = local ? gp_Pnt2d(*local)
           : start
             ? project.NextValueOfUV(gp_Pnt2d(start->uv), point, BRepExt::tolerance::exactLength)
           : i == 0 ? project.ValueOfUV(point, BRepExt::tolerance::exactLength)
                    : project.NextValueOfUV(uv, point, BRepExt::tolerance::exactLength);
      if (start
          && surface->Value(start->uv.X(), start->uv.Y()).Distance(point)
               < surface->Value(uv.X(), uv.Y()).Distance(point))
        uv = gp_Pnt2d(start->uv);
    }
    miss.g0 = std::max(miss.g0, surface->Value(uv.X(), uv.Y()).Distance(point));
    // Not at the ends: there the neighbours on both sides have their say.
    if (side.order < 1 || end)
      continue;
    const SurfaceJet mine = surfaceJet(op, surface, uv.X(), uv.Y()), theirs = side.jet(fraction);
    const double     cosine = theirs.n.Dot(mine.n);
    miss.g1 = std::max(miss.g1, std::acos(std::min(1.0, std::abs(cosine))) * 180 / M_PI);
    // The bend is read up to a corner too, but not within an ease that the corner's faces are
    // shown to ask for (see agreeAtCorners): there it is read from the first place on, as the
    // ease's own allowance was set.
    const bool nearEnd = fraction < 0.5 / places || fraction > 1 - 0.5 / places;
    if (side.order < 2 || (nearEnd && side.easing(fraction) > 0))
      continue;
    auto [bend, relative] =
      bendDifference(theirs, mine, side.tangent(fraction).Normalized(), cosine < 0 ? -1 : 1);
    if (const double easing = side.easing(fraction); easing > 0)
      relative = std::min(relative, bend / easing);
    miss.g2                             = std::max(miss.g2, bend);
    miss.g2Relative                     = std::max(miss.g2Relative, relative);
    const auto [faceBend, faceRelative] = side.eased()
                                            ? bendDifference(side.faceJet(fraction),
                                                             mine,
                                                             side.tangent(fraction).Normalized(),
                                                             cosine < 0 ? -1 : 1)
                                            : std::pair<double, double>(bend, relative);
    miss.faceG2                         = std::max(miss.faceG2, faceBend);
    miss.faceG2Relative                 = std::max(miss.faceG2Relative, faceRelative);
  }
  return miss;
}

double sharpest(const SurfaceJet& jet)
{
  const gp_Vec a = jet.u.Normalized(), b = jet.n.Crossed(a);
  const double aa = jet.curvature(a, a), bb = jet.curvature(b, b), ab = jet.curvature(a, b);
  return std::abs((aa + bb) / 2) + std::hypot((aa - bb) / 2, ab);
}

std::string agreeAtCorners(std::vector<Side>& sides,
                           const Limits&      limits,
                           double             size,
                           double             shortest)
{
  std::string eased;
  for (size_t i = 0; i < sides.size(); ++i)
  {
    Side &a = sides[i], &b = sides[(i + 1) % sides.size()];
    if (a.free || b.free)
      continue;
    // A side that is only to be met leaves the corner in a direction of its own, and the
    // surface with it: where that is not in the plane of the face next to it, the surface
    // cannot be tangent to the face there.
    if ((a.order > 0) != (b.order > 0))
    {
      const bool       first = a.order > 0;
      const SurfaceJet face  = first ? a.jet(1 - cornerInset) : b.jet(cornerInset);
      const gp_Vec     leaves =
        (first ? b.tangent(cornerInset) : a.tangent(1 - cornerInset)).Normalized();
      const double angle = std::asin(std::min(1.0, std::abs(face.n.Dot(leaves)))) * 180 / M_PI;
      if (angle > limits.g1)
        fail(op,
             unsupported_geometry,
             "boundary edges " + std::to_string(a.input + 1) + " and " + std::to_string(b.input + 1)
               + " meet at a corner where the one that is only to be met leaves the face of the "
                 "other by "
               + number(angle) + " degrees: no surface holds it and is tangent to that face there");
      continue;
    }
    if (a.order < 1 || b.order < 1 || a.support.IsSame(b.support))
      continue;
    const SurfaceJet  before = a.jet(1 - cornerInset), after = b.jet(cornerInset);
    const std::string where = "boundary edges " + std::to_string(a.input + 1) + " and "
                              + std::to_string(b.input + 1)
                              + " meet at a corner where their faces ";
    const double cosine = before.n.Dot(after.n);
    const double angle  = std::acos(std::min(1.0, std::abs(cosine))) * 180 / M_PI;
    if (angle > limits.g1)
      fail(op,
           unsupported_geometry,
           where + "are " + number(angle) + " degrees apart: no surface is tangent to both there");
    if (a.order < 2 || b.order < 2)
      continue;
    Miss miss;
    std::tie(miss.g2, miss.g2Relative) =
      bendDifference(before, after, a.tangent(1 - cornerInset).Normalized(), cosine < 0 ? -1 : 1);
    if (miss.within(limits, 2, size))
      continue;
    if (eased.empty())
      eased = where + "bend differently, by a share of " + number(miss.g2Relative)
              + " of the curvature: no surface has the curvature of both there, and none was "
                "found that eases to the one it can have near the corner and meets both faces "
                "beyond (G1 can be had)";
    const CornerBend corner = cornerBend(a, b, true);
    for (int k = 0; k < 2; ++k)
    {
      Side&        side = k ? b : a;
      const double bend = corner.across[k];
      if (bend != 0)
        side.eases[1 - k] = {bend, std::min(shortest, side.length / 2)};
    }
  }
  return eased;
}

bool hasGaps(const Asked& asked)
{
  return std::any_of(asked.sides.begin(), asked.sides.end(), [](const Side& side) {
    return side.free;
  });
}

std::optional<std::string> guideEndHint(const Asked& asked, double size)
{
  for (const Side& guide : asked.guides)
    for (double fraction : {0.0, 1.0})
    {
      const gp_Pnt point = guide.at(fraction);
      for (const Side& side : asked.sides)
      {
        if (side.order < 1)
          continue;
        GeomAPI_ProjectPointOnCurve project(point,
                                            side.curve,
                                            std::min(side.from, side.to),
                                            std::max(side.from, side.to));
        if (!project.NbPoints() || project.LowerDistance() > asked.limits.g0)
          continue;
        const double t = project.LowerDistanceParameter();
        const double f =
          std::clamp((t - side.from) / (side.to - side.from), cornerInset, 1 - cornerInset);
        const SurfaceJet jet     = side.jet(f);
        const gp_Vec     tangent = guide.tangent(fraction).Normalized();
        const double angle = std::asin(std::min(1.0, std::abs(jet.n.Dot(tangent)))) * 180 / M_PI;
        const std::string where = "guide " + std::to_string(guide.input + 1)
                                  + " ends on boundary edge " + std::to_string(side.input + 1);
        if (angle > asked.limits.g1)
          return where + " with tangent plane angle " + number(angle) + " degrees";
        if (side.order < 2)
          continue;
        gp_Pnt at;
        gp_Vec first, second;
        guide.curve->D2(guide.parameter(fraction), at, first, second);
        const double mine       = second.Dot(jet.n) / first.SquareMagnitude(),
                     theirs     = jet.curvature(tangent, tangent);
        const double difference = std::abs(mine - theirs),
                     scale      = std::max({std::abs(mine), std::abs(theirs), flatShare / size});
        if (difference > asked.limits.g2 * scale)
          return where + " with normal-curvature difference share " + number(difference / scale);
      }
    }
  return std::nullopt;
}

double rimBend(const Asked& asked)
{
  double bend  = 2 / asked.size;
  auto   along = [&](const Side& side, bool support) {
    // Between the ends, never at them: at a corner that is a pole of its face the face
    // has no differential to read.
    for (int i = 0; i < 32; ++i)
    {
      const double fraction = (i + 0.5) / 32;
      gp_Pnt       at;
      gp_Vec       first, second;
      side.curve->D2(side.parameter(fraction), at, first, second);
      const double speed = first.Magnitude();
      if (speed > 1e-12)
        bend = std::max(bend, first.Crossed(second).Magnitude() / (speed * speed * speed));
      if (!support || side.support.IsNull() || side.order < 1)
        continue;
      const auto   jet = side.jet(fraction);
      gp_Vec       a = jet.u.Normalized(), b = jet.n.Crossed(a);
      const double aa = jet.curvature(a, a), bb = jet.curvature(b, b), ab = jet.curvature(a, b);
      bend = std::max(bend, std::abs((aa + bb) / 2) + std::hypot((aa - bb) / 2, ab));
    }
  };
  for (const auto& side : asked.sides)
    along(side, true);
  for (const auto& guide : asked.guides)
    along(guide, false);
  return bend;
}

Measured measureAll(const occ::handle<Geom_Surface>& surface,
                    const Asked&                     asked,
                    const Known&                     known)
{
  PatchPhase phase("measurement");
  Measured   out;
  out.values.assign(asked.sides.size() * 4 + asked.guides.size() + asked.points.size(), 0);
  auto note = [&](bool ok, double excess, const std::string& what) {
    if (!ok && (out.within || excess > out.excess))
    {
      out.within = false;
      out.worst  = what;
    }
    out.excess = std::max(out.excess, excess);
  };
  const Limits& limits = asked.limits;
  for (const Side& side : asked.sides)
  {
    if (side.free)
      continue;
    const Miss   miss      = measure(surface, side, known);
    const double numbers[] = {miss.g0, miss.g1, miss.faceG2Relative, miss.faceG2};
    std::copy_n(numbers, 4, out.values.begin() + side.input * 4);
    std::string said = "boundary edge " + std::to_string(side.input + 1) + " is missed by "
                       + number(miss.g0) + " (the limit is " + number(limits.g0) + ")";
    if (side.order > 0)
      said += ", its tangent plane by " + number(miss.g1) + " degrees (" + number(limits.g1) + ")";
    if (side.order > 1)
      said +=
        ", its curvature by a share of " + number(miss.g2Relative) + " (" + number(limits.g2) + ")";
    // Only a miss in position is one that a larger tolerance could meet.
    if (miss.g0 > limits.g0)
      said += ". A larger tolerance may be met";
    note(miss.within(limits, side.order, asked.size),
         miss.excess(limits, side.order, asked.size),
         said);
  }
  size_t at = asked.sides.size() * 4;
  for (const Side& guide : asked.guides)
  {
    const double miss            = measure(surface, guide, known).g0;
    out.values[at + guide.input] = miss;
    note(miss <= limits.g0,
         miss / limits.g0,
         "guide " + std::to_string(guide.input + 1) + " is missed by " + number(miss)
           + " (the limit is " + number(limits.g0) + ")");
  }
  at += asked.guides.size();
  ShapeAnalysis_Surface project(surface);
  for (size_t i = 0; i < asked.points.size(); ++i)
  {
    const gp_Pnt2d uv   = project.ValueOfUV(asked.points[i], BRepExt::tolerance::exactLength);
    const double   miss = surface->Value(uv.X(), uv.Y()).Distance(asked.points[i]);
    out.values[at + i]  = miss;
    note(miss <= limits.g0,
         miss / limits.g0,
         "point " + std::to_string(i + 1) + " is missed by " + number(miss) + " (the limit is "
           + number(limits.g0) + ")");
  }
  return out;
}

void measureSpatialBoundary(Built& built, const Asked& input, bool originalDirection)
{
  if (built.measured.within && !built.measured.unkept && !built.spatialBoundary.empty())
  {
    const auto surface = BRep_Tool::Surface(built.face);
    for (const auto& rim : built.spatialBoundary)
    {
      checkCancelled(op);
      if (rim.input < 0 || static_cast<size_t>(rim.input) >= input.sides.size()
          || static_cast<size_t>(rim.input) >= built.measured.values.size() / 4)
        fail(op, no_solution, "the delivered spatial boundary has an invalid original edge");
      const auto side = std::find_if(input.sides.begin(), input.sides.end(), [&](const Side& item) {
        return item.input == rim.input;
      });
      if (side == input.sides.end() || rim.curve.IsNull())
        fail(op, no_solution, "the delivered spatial boundary has no original edge");
      double miss = 0;
      // Dense delivered-trace measurement, alongside the existing held-PC gate;
      // this is not an analytic bound between the measured parameters.
      for (int i = 0; i <= 512; ++i)
      {
        if ((i & 31) == 0)
          checkCancelled(op);
        const double   fraction = i / 512.0;
        const gp_Pnt2d uv       = rim.curve->Value(rim.from + fraction * (rim.to - rim.from));
        if (!std::isfinite(uv.X()) || !std::isfinite(uv.Y()))
          fail(op, no_solution, "the delivered spatial boundary parameter is not finite");
        const double gap =
          surface->Value(uv.X(), uv.Y())
            .Distance(
              side->at(originalDirection && rim.originalReversed ? 1 - fraction : fraction));
        if (!std::isfinite(gap))
          fail(op, no_solution, "the delivered spatial boundary is not finite");
        miss = std::max(miss, gap);
      }
      built.measured.values[rim.input * 4] = std::max(built.measured.values[rim.input * 4], miss);
      if (miss > input.limits.g0)
      {
        built.measured.within = false;
        const double excess   = miss / input.limits.g0;
        if (excess > built.measured.excess)
        {
          built.measured.worst = "original spatial boundary edge " + std::to_string(rim.input + 1)
                                 + " is missed by " + number(miss) + " (the limit is "
                                 + number(input.limits.g0) + ")";
          built.measured.excess = excess;
        }
      }
    }
  }
}

void measureTrimmedInterior(Built& built, const Asked& asked)
{
  if (asked.points.empty() && asked.guides.empty())
    return;
  const auto                 interior    = trimmedInteriorCells(built.face);
  const auto                 faceSurface = BRep_Tool::Surface(built.face);
  ShapeAnalysis_Surface      projection(faceSurface);
  std::optional<gp_Pnt2d>    previousUV;
  BRepExtrema_DistShapeShape query;
  query.LoadS2(built.face);

  struct Witness
  {
    double distance;
    gp_XY  uv;
  };

  std::map<std::array<double, 3>, Witness> witnesses;
  const auto distance = [&](const gp_Pnt& point, gp_XY* uvOut = nullptr) {
    checkCancelled(op);
    const std::array<double, 3> key{point.X(), point.Y(), point.Z()};
    if (!std::all_of(key.begin(), key.end(), [](double x) { return std::isfinite(x); }))
      fail(op, no_solution, "an interior constraint has a nonfinite point");
    if (const auto found = witnesses.find(key); found != witnesses.end())
    {
      if (uvOut)
        *uvOut = found->second.uv;
      return found->second.distance;
    }
    // A local projection that supplies an in-trim face witness is sufficient
    // for an upper distance bound. Fall back to the global trimmed query whenever
    // it cannot meet the request; it never supplies a negative feasibility verdict.
    try
    {
      const auto local =
        previousUV ? nearestWithin(faceSurface, previousUV->XY(), point) : std::nullopt;
      const auto   uv = local ? gp_Pnt2d(*local) : projection.ValueOfUV(point, asked.limits.g0 / 8);
      const double d  = point.Distance(faceSurface->Value(uv.X(), uv.Y()));
      if (std::isfinite(d) && d <= asked.limits.g0 && interior(uv.X(), uv.X(), uv.Y(), uv.Y()))
      {
        previousUV = uv;
        witnesses.emplace(key, Witness{d, uv.XY()});
        if (uvOut)
          *uvOut = uv.XY();
        return d;
      }
    }
    catch (const Standard_Failure&)
    {
    }
    query.LoadS1(BRepBuilderAPI_MakeVertex(point).Vertex());
    query.Perform();
    if (!query.IsDone() || query.NbSolution() == 0 || !std::isfinite(query.Value()))
      fail(op, no_solution, "an interior constraint has no resolved trimmed-face distance");
    gp_XY uv(HUGE_VAL, HUGE_VAL);
    if (query.SupportTypeShape2(1) == BRepExtrema_IsInFace)
    {
      double u, v;
      query.ParOnFaceS2(1, u, v);
      uv = gp_XY(u, v);
    }
    if (uvOut)
      *uvOut = uv;
    // Distance to the returned face point is an upper witness, not merely the
    // optimizer's minimum-distance estimate. Reuse only exactly equal coordinates.
    const double d = point.Distance(query.PointOnShape2(1));
    witnesses.emplace(key, Witness{d, uv});
    return d;
  };
  bool convexPlane = BRepAdaptor_Surface(built.face).GetType() == GeomAbs_Plane;
  int  wires       = 0;
  for (TopExp_Explorer it(built.face, TopAbs_WIRE); it.More(); it.Next())
    ++wires;
  convexPlane = convexPlane && wires == 1;
  std::vector<gp_XY> polygon;
  for (BRepTools_WireExplorer it(BRepTools::OuterWire(built.face), built.face);
       convexPlane && it.More();
       it.Next())
  {
    BRepAdaptor_Curve edge(it.Current());
    convexPlane = edge.GetType() == GeomAbs_Line;
    double     a, b;
    const auto pc = BRep_Tool::CurveOnSurface(it.Current(), built.face, a, b);
    if (pc.IsNull())
    {
      convexPlane = false;
      break;
    }
    polygon.push_back(pc->Value(it.Current().Orientation() == TopAbs_REVERSED ? b : a).XY());
  }
  double orientation = 0;
  for (size_t i = 0; convexPlane && i < polygon.size(); ++i)
  {
    const double turn =
      (polygon[(i + 1) % polygon.size()] - polygon[i])
        .Crossed(polygon[(i + 2) % polygon.size()] - polygon[(i + 1) % polygon.size()]);
    if (orientation * turn < 0)
      convexPlane = false;
    if (turn != 0)
      orientation = turn;
  }
  convexPlane        = convexPlane && orientation != 0;
  const auto carrier = occ::down_cast<Geom_BSplineSurface>(BRep_Tool::Surface(built.face));
  std::unique_ptr<GeomConvert_BSplineSurfaceToBezierSurface> patches;
  std::vector<double>                                        us, vs;
  if (!carrier.IsNull() && !carrier->IsURational() && !carrier->IsVRational())
  {
    patches = std::make_unique<GeomConvert_BSplineSurfaceToBezierSurface>(carrier);
    NCollection_Array1<double> u(1, patches->NbUPatches() + 1), v(1, patches->NbVPatches() + 1);
    patches->UKnots(u);
    patches->VKnots(v);
    for (double x : u)
      us.push_back(x);
    for (double x : v)
      vs.push_back(x);
  }
  const auto choose = [](int n, int k) {
    double c = 1;
    for (int i = 1; i <= k; ++i)
      c *= double(n - i + 1) / i;
    return c;
  };
  const auto pairedCover = [&](const Handle(Geom_BSplineCurve)& part, double allowance) {
    if (carrier.IsNull() || carrier->IsURational() || carrier->IsVRational() || part->IsRational())
      return false;
    gp_XY a, b;
    if (distance(part->StartPoint(), &a) > allowance || distance(part->EndPoint(), &b) > allowance
        || !std::isfinite(a.X() + a.Y() + b.X() + b.Y()))
      return false;
    const double u0 = std::min(a.X(), b.X()), u1 = std::max(a.X(), b.X()),
                 v0 = std::min(a.Y(), b.Y()), v1 = std::max(a.Y(), b.Y());
    if (!interior(u0, u1, v0, v1))
      return false;
    try
    {
      // A chord can cross several polynomial cells. Partition its certificate
      // at those crossings and the guide's own knots, instead of forcing the
      // geometric subdivision to converge on each non-dyadic knot parameter.
      const gp_XY         origin = a, direction = b - a;
      const double        first = part->FirstParameter(), last = part->LastParameter();
      std::vector<double> cuts{0., 1.};
      for (int axis = 1; axis <= 2; ++axis)
        if (direction.Coord(axis) != 0)
          for (double knot : axis == 1 ? us : vs)
          {
            const double t = (knot - origin.Coord(axis)) / direction.Coord(axis);
            if (t > 0 && t < 1)
              cuts.push_back(t);
          }
      for (double knot : part->Knots())
      {
        const double t = (knot - first) / (last - first);
        if (t > 0 && t < 1)
          cuts.push_back(t);
      }
      std::sort(cuts.begin(), cuts.end());
      cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
      const auto polynomialCover = [&](const gp_XY& a, const gp_XY& b, double first, double last) {
        const double u0 = std::min(a.X(), b.X()), u1 = std::max(a.X(), b.X()),
                     v0 = std::min(a.Y(), b.Y()), v1 = std::max(a.Y(), b.Y());
        Handle(Geom_BezierCurve) trace;
        // Restrict to one polynomial cell before taking its diagonal. Degenerate
        // UV rectangles are exact isoparametric curves instead.
        if (u1 == u0 || v1 == v0)
        {
          auto         iso   = u1 == u0 ? carrier->UIso(u0) : carrier->VIso(v0);
          const double first = u1 == u0 ? v0 : u0, last = u1 == u0 ? v1 : u1;
          if (first == last)
            return false;
          GeomConvert_BSplineCurveToBezierCurve spans(
            GeomConvert::CurveToBSplineCurve(new Geom_TrimmedCurve(iso, first, last)));
          if (spans.NbArcs() != 1)
            return false;
          trace = spans.Arc(1);
          if ((u1 == u0 ? a.Y() > b.Y() : a.X() > b.X()))
            trace->Reverse();
        }
        else
        {
          const auto span = [](const std::vector<double>& k, double a, double b) {
            const int i = std::clamp(int(std::upper_bound(k.begin(), k.end(), a) - k.begin()) - 1,
                                     0,
                                     int(k.size()) - 2);
            return a >= k[i] && b <= k[i + 1] ? i : -1;
          };
          const int iu = span(us, u0, u1), iv = span(vs, v0, v1);
          if (iu < 0 || iv < 0)
            return false;
          auto patch = occ::down_cast<Geom_BezierSurface>(patches->Patch(iu + 1, iv + 1)->Copy());
          patch->Segment((u0 - us[iu]) / (us[iu + 1] - us[iu]),
                         (u1 - us[iu]) / (us[iu + 1] - us[iu]),
                         (v0 - vs[iv]) / (vs[iv + 1] - vs[iv]),
                         (v1 - vs[iv]) / (vs[iv + 1] - vs[iv]));
          const int                  p = patch->UDegree(), q = patch->VDegree(), degree = p + q;
          NCollection_Array1<gp_Pnt> poles(1, degree + 1);
          for (int k = 0; k <= degree; ++k)
          {
            gp_XYZ value(0, 0, 0);
            for (int i = std::max(0, k - q); i <= std::min(p, k); ++i)
            {
              const int j = k - i;
              value +=
                patch->Pole(a.X() > b.X() ? p - i + 1 : i + 1, a.Y() > b.Y() ? q - j + 1 : j + 1)
                  .XYZ()
                * (choose(p, i) * choose(q, j) / choose(degree, k));
            }
            poles(k + 1) = gp_Pnt(value);
          }
          trace = new Geom_BezierCurve(poles);
        }
        GeomConvert_BSplineCurveToBezierCurve guide(part, first, last, 0.);
        if (guide.NbArcs() != 1)
          return false;
        auto      c      = guide.Arc(1);
        const int degree = std::max(c->Degree(), trace->Degree());
        c->Increase(degree);
        trace->Increase(degree);
        double bound = 0;
        for (int i = 1; i <= degree + 1; ++i)
          bound = std::max(bound, c->Pole(i).Distance(trace->Pole(i)));
        return std::nextafter(bound, HUGE_VAL) <= allowance;
      };
      for (size_t i = 1; i < cuts.size(); ++i)
        if (!polynomialCover(origin + direction * cuts[i - 1],
                             origin + direction * cuts[i],
                             first + (last - first) * cuts[i - 1],
                             first + (last - first) * cuts[i]))
          return false;
      return true;
    }
    catch (const Standard_Failure&)
    {
      return false;
    }
  };
  size_t at = asked.sides.size() * 4;
  for (const Side& guide : asked.guides)
  {
    const double first = std::min(guide.from, guide.to), last = std::max(guide.from, guide.to);
    const Handle(Geom_Curve)  cut = new Geom_TrimmedCurve(guide.curve, first, last);
    Handle(Geom_BSplineCurve) curve;
    double                    approximation = 0;
    if (const auto type = GeomAdaptor_Curve(cut).GetType();
        type == GeomAbs_OtherCurve || type == GeomAbs_OffsetCurve)
    {
      // Keep the approximation's reported error in the distance budget; unlike
      // CurveToBSplineCurve's fixed offset allowance, this follows the request.
      GeomConvert_ApproxCurve fit(cut, asked.limits.g0 / 8, GeomAbs_C1, 1000, 14);
      if (!fit.IsDone() || !std::isfinite(fit.MaxError()))
        fail(op, no_solution, "a guide has no resolved spline approximation for trim measurement");
      curve         = fit.Curve();
      approximation = fit.MaxError();
    }
    else
      curve = GeomConvert::CurveToBSplineCurve(cut);
    double     observed = 0;
    const auto covered  = [&](auto&& self, const Handle(Geom_BSplineCurve)& part) -> bool {
      checkCancelled(op);
      const double a = part->FirstParameter(), b = part->LastParameter(), m = a + (b - a) / 2;
      const gp_Pnt point = part->Value(m);
      const double miss  = distance(point) + approximation;
      observed           = std::max(observed, miss);
      if (miss > asked.limits.g0)
        return false;
      double radius      = 0;
      bool   convexCover = convexPlane;
      for (int i = 1; i <= part->NbPoles(); ++i)
      {
        if (!(part->Weight(i) > 0))
          return false;
        radius = std::max(radius, point.Distance(part->Pole(i)));
        if (convexCover && distance(part->Pole(i)) + approximation > asked.limits.g0)
          convexCover = false;
      }
      if (convexCover || miss + std::nextafter(radius, HUGE_VAL) <= asked.limits.g0
          || pairedCover(part, asked.limits.g0 - approximation))
        return true;
      if (m == a || m == b)
        return false;
      auto left = occ::down_cast<Geom_BSplineCurve>(part->Copy());
      left->Segment(a, m);
      auto right = occ::down_cast<Geom_BSplineCurve>(part->Copy());
      right->Segment(m, b);
      return self(self, left) && self(self, right);
    };
    if (curve.IsNull() || !covered(covered, curve))
      fail(op,
           no_solution,
           "guide " + std::to_string(guide.input + 1)
             + " is not within tolerance of the delivered trimmed face");
    built.measured.values[at + guide.input] =
      std::max(built.measured.values[at + guide.input], observed);
  }
  at += asked.guides.size();
  for (size_t i = 0; i < asked.points.size(); ++i)
  {
    const double miss             = distance(asked.points[i]);
    built.measured.values[at + i] = std::max(built.measured.values[at + i], miss);
    if (miss > asked.limits.g0)
      fail(op,
           no_solution,
           "point " + std::to_string(i + 1)
             + " is outside the tolerance of the delivered trimmed face");
  }
}

double cornerConflict(const Side& a, const Side& b)
{
  return cornerBend(a, b, false).conflict;
}

void verifySewn(Built& built, const Asked& asked)
{
  PatchPhase          phase("sewing");
  std::vector<double> delivered(built.measured.values.size() / 4, -1);
  // Sewing may replace a p-curve even while the surface and its spatial edge stay put.
  // The preimage's residual alone cannot predict that delivered representation.
  BRepBuilderAPI_Sewing    sewing(asked.limits.g0);
  std::vector<TopoDS_Face> supports;
  for (const Side& side : asked.sides)
    if (!side.support.IsNull()
        && std::none_of(supports.begin(), supports.end(), [&](const TopoDS_Face& face) {
             return face.IsSame(side.support);
           }))
    {
      supports.push_back(side.support);
      sewing.Add(side.support);
    }
  if (supports.empty())
    return;
  sewing.Add(built.face);
  sewing.Perform();
  TopoDS_Face sewnPatch;
  for (TopExp_Explorer it(sewing.SewedShape(), TopAbs_FACE); it.More(); it.Next())
  {
    const auto face = TopoDS::Face(it.Current());
    if (BRep_Tool::Surface(face) == BRep_Tool::Surface(built.face))
      sewnPatch = face;
  }
  if (!sewnPatch.IsNull())
  {
    for (const TopoDS_Face& original : supports)
    {
      TopoDS_Face     support;
      TopLoc_Location location;
      const auto      held = BRep_Tool::Surface(original, location);
      for (TopExp_Explorer it(sewing.SewedShape(), TopAbs_FACE); it.More(); it.Next())
      {
        const auto      face = TopoDS::Face(it.Current());
        TopLoc_Location loc;
        if (BRep_Tool::Surface(face, loc) == held && loc == location)
          support = face;
      }
      if (support.IsNull())
        continue;
      for (TopExp_Explorer it(sewnPatch, TopAbs_EDGE); it.More(); it.Next())
      {
        const auto edge   = TopoDS::Edge(it.Current());
        bool       shared = false;
        for (TopExp_Explorer other(support, TopAbs_EDGE); other.More(); other.Next())
          shared = shared || edge.IsSame(other.Current());
        if (!shared || BRep_Tool::Degenerated(edge))
          continue;
        const BRepAdaptor_Curve curve(edge);
        const double            from = curve.FirstParameter(), to = curve.LastParameter();
        double                  gap = 0;
        for (int k = 0; k < places; ++k)
        {
          const double t = from + (to - from) * (k + .5) / places;
          const auto   a = surfaceJet(op, edge, sewnPatch, t);
          const auto   b = surfaceJet(op, edge, support, t);
          gap            = std::max(gap, a.p.Distance(b.p));
        }
        double allowance = asked.limits.g0;
        for (const Side& side : asked.sides)
        {
          if (!side.support.IsSame(original)
              || side.at(.5).Distance(curve.Value((from + to) / 2)) > 1e-5 * asked.size)
            continue;
          double geometric = 0, parameter = 0;
          if (gap > allowance)
          {
            for (int k = 0; k < places; ++k)
            {
              const double f     = (k + .5) / places;
              const gp_Pnt point = side.at(f);
              const gp_XY  uv =
                nearest(side.held, side.onHeld->Value(side.parameter(f)).XY(), point);
              geometric = std::max(geometric, point.Distance(side.held->Value(uv.X(), uv.Y())));
              parameter = std::max(parameter, point.Distance(side.meet(f)));
            }
            // An inaccurate spatial edge is not an allowance for a bad joint. Only
            // existing tangential reparameterisation of an edge on its support is.
            if (geometric <= asked.limits.g0)
              allowance = std::max(allowance, parameter);
          }
          // Report the delivered representation's gap, even when the input edge
          // already has a tangential parameter discrepancy.
          delivered[side.input] = std::max(delivered[side.input], gap);
        }
        if (gap > allowance)
          fail(op,
               no_solution,
               "the joint after sewing is missed by " + number(gap)
                 + " (the limit including the input edge's tangential parameter uncertainty is "
                 + number(allowance) + "): its edge parameter curves do not agree");
      }
    }
  }
  // Construction has already passed its spatial-boundary checks. The public G0 report
  // describes the delivered seam, not the maximum of two different measurements. Keep
  // the construction reading only where no sewn seam could be identified.
  for (size_t i = 0; i < delivered.size(); ++i)
    if (delivered[i] >= 0)
      built.measured.values[i * 4] = delivered[i];
}

} // namespace BRepExt::patch_fill
