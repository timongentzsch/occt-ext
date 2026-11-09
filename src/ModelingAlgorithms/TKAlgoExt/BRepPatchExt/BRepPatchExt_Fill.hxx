// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once

// The internals of Patch shared between its files: what a loop asks, what a fill built is and
// how it measured, and the stages one file hands to another. Not installed; see
// BRepPatchExt_Patch.cxx for the order in which the stages are tried.

#include <BRepPatchExt_Chart.hxx>
#include <BRepPatchExt_Inside.hxx>
#include <BRepPatchExt_Side.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <GeomAlgoExt_PatchFit.hxx>
#include <GeomAlgoExt_SurfaceJet.hxx>

#include <Geom_BSplineCurve.hxx>
#include <gp_Trsf.hxx>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace BRepExt::patch_fill
{

const char* const op = "patch";

//! How many places along an edge the result is measured at: those of Measure Continuity.
constexpr int places = 129;

//! Below this share of a G2 limit, as a difference of curvature times the loop's size, two
//! bends count as the same though one of them is flat and every difference is all of it.
constexpr double flatShare = 1e-3;

//! What the fairness of a trimmed surface counts for outside its loop.
constexpr double outsideWeight = 1e-3;

//! The steps to the unit of the grid that a fit in a loop's frame rounds its coordinates, in
//! units of the loop's size, and its normals to: about the square root of the rounding of a
//! double (2.2e-16), far above what a rotation rounds differently and far below any tolerance.
constexpr double frameSteps = 1e8;

//! Where a surface goes over into planes and other faces at once (or into a cone alone), it is
//! moved onto its curves, and its edges' curves in its parameters are made, to a tenth of
//! OpenCascade's least distance, Precision::Confusion(), when the tolerance asks less: its
//! edges then keep that least distance as their tolerance (see `trimmed`), as the faces they
//! are sewn to do, and STEP gives the joint back as it was (at half of 1e-6 the cap of a cone
//! read 2.0e-7 before STEP and 2.8e-3 after, the corner 2 / 3 / 5 5.4e-7 and 1.7e-6). It is
//! OpenCascade's precision, which is absolute, and not a share of the model's size. A surface
//! within twice that is left as it is.
constexpr double exchangeWithin = Precision::Confusion() / 10;

//! The steepest a trimmed surface may stand against the loop's mean plane, as the cosine: 80
//! degrees. A surface over a plane is a height above it, and the slope of that height is the
//! tangent of this angle; at 80 degrees it is 5.7 and the net's spans, even across the plane,
//! are six times as long on the surface.
constexpr double steepest = 0.17;

//! How far outside the loop's polygon a guide or a point may lie and still count as inside, as
//! a share of the loop's size: the polygon has 32 corners per side and cuts off that much of a
//! curved side (a full circle of 32 chords is 0.5 % of its radius inside the circle).
constexpr double outsideSlack = 0.01;

struct Limits
{
  double g0, g1, g2;
};

//! The largest misses along an edge: a length, degrees, a difference of curvature and that as
//! a share (see bendDifference).
struct Miss
{
  double g0 = 0, g1 = 0, g2 = 0, g2Relative = 0;
  //! The same two against the face itself, where the side eases near a corner (see
  //! Side::Ease): what is reported, while the side is held to what it asks.
  double faceG2 = 0, faceG2Relative = 0;

  bool within(const Limits& limits, int order, double size) const
  {
    return g0 <= limits.g0 && (order < 1 || g1 <= limits.g1)
           && (order < 2 || g2Relative <= limits.g2 || g2 * size <= limits.g2 * flatShare);
  }

  //! The largest of the misses as a multiple of its limit.
  double excess(const Limits& limits, int order, double size) const
  {
    return std::max(
      {g0 / limits.g0,
       order < 1 ? 0 : g1 / limits.g1,
       order < 2 ? 0 : std::min(g2Relative / limits.g2, g2 * size / (limits.g2 * flatShare))});
  }
};

using patch_side::Side;

using patch_side::chain;

//! The plane a loop lies nearest to, by the area its polygon encloses seen from every side.
struct Frame
{
  gp_Pnt origin;
  gp_Vec x, y, z;

  gp_XY project(const gp_Pnt& point) const
  {
    const gp_Vec to(origin, point);
    return {to.Dot(x), to.Dot(y)};
  }

  double height(const gp_Pnt& point) const { return gp_Vec(origin, point).Dot(z); }
};

std::vector<gp_Pnt> polygon(const std::vector<Side>& sides);

Frame meanPlane(const std::vector<gp_Pnt>& points, double perimeter);

//! The area a polygon encloses, negative when it runs clockwise. Fails where it crosses itself.
double enclosed(const std::vector<gp_XY>& loop);

//! Whether `point` is no further than `slack` from the side of a polygon from `from` to `to`.
bool near(const gp_XY& from, const gp_XY& to, const gp_XY& point, double slack);

//! Where the side of a polygon from `from` to `to` crosses the line at height `y`, if it does.
std::optional<double> crossing(const gp_XY& from, const gp_XY& to, double y);

//! Whether a point is inside a polygon, or no further outside it than `slack`: the polygon
//! stands for curves, which a point on them may be a little outside of.
bool inside(const std::vector<gp_XY>& loop, const gp_XY& point, double slack);

//! `inside` for many points of one polygon, as a fit reads its density over a grid: per
//! height, the polygon's crossings are found once and sorted, and only sides that come near a
//! point are measured. The same answers as `inside`.
class Inside
{
public:
  Inside(const std::vector<gp_XY>& loop, double slack)
      : loop(loop),
        slack(slack)
  {
    for (const gp_XY& point : loop)
      reach = std::max({reach, std::abs(point.X()), std::abs(point.Y())});
    // Beyond this from its box no side's measured distance is within `slack`.
    reach = slack + 1e-9 * (1 + reach);
  }

  bool operator()(const gp_XY& point) const
  {
    auto [at, made] = rows.try_emplace(point.Y());
    Row& row        = at->second;
    if (made)
      for (size_t i = 0, j = loop.size() - 1; i < loop.size(); j = i++)
      {
        if (const auto x = crossing(loop[j], loop[i], point.Y()))
          row.crossings.push_back(*x);
        if (std::min(loop[i].Y(), loop[j].Y()) - reach <= point.Y()
            && point.Y() <= std::max(loop[i].Y(), loop[j].Y()) + reach)
          row.sides.push_back(j);
      }
    if (made)
      std::sort(row.crossings.begin(), row.crossings.end());
    for (size_t j : row.sides)
    {
      const gp_XY &from = loop[j], &to = loop[(j + 1) % loop.size()];
      if (std::min(from.X(), to.X()) - reach <= point.X()
          && point.X() <= std::max(from.X(), to.X()) + reach && near(from, to, point, slack))
        return true;
    }
    const auto beyond = std::upper_bound(row.crossings.begin(), row.crossings.end(), point.X());
    return (row.crossings.end() - beyond) % 2 == 1;
  }

private:
  struct Row
  {
    std::vector<double> crossings;
    std::vector<size_t> sides;
  };

  const std::vector<gp_XY>&               loop;
  double                                  slack, reach = 0;
  mutable std::unordered_map<double, Row> rows;
};

//! What the surface is to meet at one place: a point, and with an order the tangent plane and
//! the bend of `jet` there. `across` is a direction of the surface's parameters that leaves
//! the curve the place is on.
struct Target
{
  gp_Pnt                    point;
  gp_XY                     uv;
  int                       order;
  std::optional<SurfaceJet> jet;
  gp_XY                     across;
  //! The share of the loop's size that the place stands for.
  double share;
  //! With order 1: the neighbour's bend is wished for as well, and not asked.
  bool wish = false;
  //! Read and not asked: see FitCondition.
  bool                  watch         = false;
  std::optional<gp_Vec> perpendicular = {};
  //! A wish and not a condition: its point is weighed against the energies and not met.
  bool soft = false;
  //! Near a corner, the least bend a miss of the bend asked is a share of (see Side::easing).
  double easing = 0;
  //! Over a chart whose lengths along the rim are the rim's own: the way the surface leaves the
  //! curve in space, a unit vector, along which its derivative `across` is to be within a band
  //! about one (see overChart). Read and not asked until it is missed: a band and no wish.
  std::optional<gp_Vec> leave      = {};
  bool                  leaveAsked = false;
};

//! What tangency and bend count for against position, in the units of the loop's size. Where
//! not all can be had, at a corner that two neighbours bend differently into, the lower order
//! is kept and the higher one gives.
constexpr double tangentWeight = 1e-3, bendWeight = 1e-6;

//! What the energy of bending counts for against a condition of weight one. The conditions are
//! met by multipliers, whatever it is; it decides how near the first solution is to them, and
//! with that how many passes the multipliers need, and must stay well above the rounding of
//! the sums it is added to (1e-16 of the conditions' weights).
constexpr double bendingWeight = 1e-11;

//! How near its value a condition of the fit is brought, as a share of the limit that the
//! finished surface is held to. The fit reads at its own places and the surface is measured
//! between them, where it is a little further off: a quarter leaves room for that.
constexpr double fitShare = 0.25;

//! The same for the places between, which the fit only reads: half, since the surface is
//! measured at yet other places, and a miss that is at its largest between two places of the
//! fit is found there.
constexpr double watchShare = 0.5;

enum class GraphFit
{
  none,
  transverse,
  complete
};

//! Fit the free poles: once, or where a bend is asked or wished for several times over, the
//! bend each time about the surface before. With `along`, the poles are free along it only.
occ::handle<Geom_BSplineSurface> fit(SplineNet&                                   net,
                                     const std::vector<char>&                     fixed,
                                     std::vector<Target>                          targets,
                                     double                                       size,
                                     const Limits&                                limits,
                                     occ::handle<Geom_BSplineSurface>             surface = {},
                                     const std::function<double(double, double)>& density = {},
                                     const gp_XYZ*                                along   = nullptr,
                                     const std::vector<FitJump>&                  jumps   = {},
                                     const Frame*                                 frame   = nullptr,
                                     double                                       flatness = -1,
                                     bool     leastSquares                                 = false,
                                     GraphFit graph = GraphFit::none);

std::string number(double value);

//! Where on the surface a place of a side is, for the constructions that know. `own` says that
//! this is where the face itself holds the place, by the curve of its edge in the surface's
//! parameters: the place is then read there and nowhere else, as Measure Continuity reads it.
//! Otherwise it is where to start looking for the surface's nearest point, which has then not
//! to be looked for over a surface that may fold outside the loop.
struct Place
{
  gp_XY uv;
  bool  own;
};

using Known = std::function<std::optional<Place>(const Side&, double fraction, const gp_Pnt&)>;

//! How far from a corner, as a share of a side, the face along the side is read for what it
//! asks at the corner: off the corner itself, where the faces of both sides meet and a face
//! may have no one direction.
constexpr double cornerInset = 1e-6;

//! The fractions of a side that the measure reads, in order: the places of Measure Continuity,
//! and between the first and the last of them and the side's ends places that halve their way to
//! each end down to `cornerInset`, so that what a surface does up to a corner is read as well. How
//! near is at the corner is a judgement, a fraction of the side and so the same at every size: read
//! down to 1e-5 or 1e-7 instead, the Patch tests came out the same (see
//! dox/extensions/operations/BRepPatchExt.md).
const std::vector<double>& measuredFractions();

//! How far a surface is from what a side asks, read from the surface alone.
Miss measure(const occ::handle<Geom_Surface>& surface, const Side& side, const Known& known);

//! The largest principal curvature of a surface at a place, from its jet.
double sharpest(const SurfaceJet& jet);

//! One surface has one tangent plane and one bend at a corner. Where the faces on both sides
//! of a corner ask for different tangent planes, no patch can give both: say so, and where.
//! Where they ask for different bends, a surface that is twice differentiable at the corner has
//! the one bend there that the edges and the tangent planes along them leave it (cornerBend),
//! and no other: the bend each side asks across itself eases to that one near the corner (see
//! Side::Ease), and what the face itself is missed by there is reported, not refused. Over how
//! long: with the faces' bends met along both sides beyond the ease, the surface's bend has to
//! turn from the one side's to the other's across the corner at every distance d from it, by as
//! much as they differ, so that it changes by that over a length of about d: as fast as the net
//! has spans for only where d is no shorter than a span. An ease is therefore no shorter than
//! `shortest`, a span of the coarsest net (on the corners of three arms of rounds, eases of 0.3,
//! 3.9 and 7.8 against that span's 24 left every net short of the bends asked or bulging inside,
//! 16 and 24 were met), and no longer than half its side. Returns what is to be said where no
//! fill eases so and meets the faces elsewhere (three rounds round a corner at a millionth,
//! unequal rounds and their planes round a spherical cut): the first such corner; empty where
//! there is none.
std::string agreeAtCorners(std::vector<Side>& sides,
                           const Limits&      limits,
                           double             size,
                           double             shortest);

struct Asked;
bool hasGaps(const Asked& asked);

struct Asked
{
  std::vector<Side>   sides, guides;
  std::vector<gp_Pnt> points;
  Limits              limits;
  double              size;
  double              flatness = -1;
  PatchControls       controls;

  // Guides and boundaries occupy separate arrays. Pointer subtraction between
  // them is undefined, even when its result is only used to test membership.
  std::optional<size_t> boundaryIndex(const Side& side) const
  {
    for (size_t i = 0; i < sides.size(); ++i)
      if (&side == &sides[i])
        return i;
    return std::nullopt;
  }
};

bool hasGaps(const Asked& asked);

//! Work on a private, consistently directed spline rim in its intrinsic frame. Equivalent
//! pole reversal and rigid placement must be removed before sampling and net decisions,
//! rather than rounding only the eventual linear system. Original curves are measured later.
std::optional<gp_Trsf> localSplineRim(Asked& asked, bool* copyBudgetMiss = nullptr);

//! A failed world-space copy may instead be represented relative to an exact anchor.
//! TwoDiff certifies each subtraction; no source-pole displacement is spent here. Keep
//! the anchor as a face placement, so it never gets rounded back into the fitted poles.
std::optional<gp_Trsf> rebasedSplineRim(Asked& asked, gp_Trsf& anchor);

//! Endpoint jet information for an already rejected fit; never a feasibility veto.
std::optional<std::string> guideEndHint(const Asked& asked, double size);

//! The tightest bend on the rim of what is asked: of its curves and guides, of the faces
//! that it is to run on from, and of a ball as large as the loop.
double rimBend(const Asked& asked);

//! The misses in the order of the input: per boundary edge four numbers, then per guide and
//! per point the distance. `worst` says what the first miss beyond the limits is.
struct Measured
{
  std::vector<double> values;
  bool                within = true;
  std::string         worst;
  //! The largest miss as a multiple of its limit.
  double excess = 0;
  //! Within the limits along the edges, and not a surface to keep: `worst` says why.
  bool unkept = false;
};

Measured measureAll(const occ::handle<Geom_Surface>& surface,
                    const Asked&                     asked,
                    const Known&                     known = {});

//! A surface with the face made of it and how it measures.
struct Built
{
  int         construction;
  TopoDS_Face face;
  Measured    measured;
  //! For a surface over a conformal chart, how square its U and V come out inside the loop,
  //! a goal and not a limit: the root mean square, the 95th percentile and the largest angle
  //! from square in degrees, then the least and largest ratio |S_u| / |S_v|.
  std::vector<double> squareness = {};
  //! Whether the soft interior-shape target was met, and how many nets were tried.
  std::vector<double> shapeTarget = {};

  // The actual delivered boundary traces of the new free height-graph construction.
  // UV is unchanged by restoring its world placement; certify original spatial curves
  // separately from a support's potentially different owned parameter curve.
  struct SpatialRim
  {
    int                       input;
    occ::handle<Geom2d_Curve> curve;
    double                    from, to;
    bool                      originalReversed;
  };

  std::vector<SpatialRim> spatialBoundary = {};
  //! Per side, the curve of a gap in the surface's parameters; null for every other side.
  std::vector<occ::handle<Geom2d_Curve>> gaps = {};
  //! For a quad of two rims across two gaps, how it crosses from rim to rim (see insideFault).
  std::optional<Crossing> crossing = {};
};

// The local construction keeps fraction direction; canonical restoration may reverse it.
// Measure the actual delivered rim, not an unconstrained surface projection.
void measureSpatialBoundary(Built& built, const Asked& input, bool originalDirection);

//! Distance to the delivered face, including its trim. A guide's positive-weight
//! spline hull is subdivided until its whole image has a tolerance-ball cover.
void measureTrimmedInterior(Built& built, const Asked& asked);

std::vector<double> clampedKnots(double from, double to, int spans, int degree);

//! A harmonic height polynomial for a free rim, or a quadratic fixed by interior positions.
//! Seed the exact polynomial rather than driving tiny interpolation errors into the fit.
//! The later correction and measurement hold the same limits as for a general fitted surface.
bool polynomialHeight(const Asked& asked,
                      const Frame& frame,
                      SplineNet&   net,
                      bool         approximate = false);

//! A place where the net is to be finer than it is elsewhere: spans no longer than `span`
//! from `at - span` to `at + span`.
struct Finer
{
  double at, span;
};

//! The knots of `spans` even spans, made finer where `finer` says and no coarser anywhere:
//! each knot where the number of spans wanted up to there is whole. Then `breaks` among them
//! at one less than the degree times: the surface is there only once differentiable across,
//! and free to change its bend. A knot that is nearer to a break than `snap` of a span gives
//! way to it, so that no span is left a sliver.
std::vector<double> gradedKnots(double                     from,
                                double                     to,
                                int                        spans,
                                int                        degree,
                                const std::vector<Finer>&  finer,
                                const std::vector<double>& breaks);

//! By what share of their curvature the faces on the two sides of a corner disagree about the
//! one surface that is to leave both tangentially (see cornerBend).
double cornerConflict(const Side& a, const Side& b);

std::optional<TopoDS_Face> planar(const Asked& asked, const Frame& frame);

//! A side as a polynomial B-spline over [0, length] with the side's own run of parameters, or
//! nothing where no such curve comes within the tolerance.
occ::handle<Geom_BSplineCurve> plain(const Side& side,
                                     double      length,
                                     bool        turned,
                                     double      tolerance);

//! How near each other, as a share of the loop's size, the surfaces of two nets in a row of the
//! ladder are to count as having come as near the fairest surface
//! as finer nets would take them. A share of the size, so that the stop is the same at every size
//! and pose; the length tolerance, which is absolute, would not be, and nothing inside the loop is
//! asked to any tolerance. The value is a judgement: between the nets of the tests, these distances
//! spread from 1e-5 to 1 of the size with no gap between settled and unsettled pairs; an
//! eight-sided seat surface, two domed tops and a hull bottom that a coarse net filled 0.6 to 9 off
//! their masters stop on the same nets for any value from 3.1e-4 to 2.4e-3 (see
//! dox/extensions/operations/BRepPatchExt.md).
constexpr double settled = 1e-3;

//! How far apart the surfaces of two nets are: from places of the coarser, `countU` x `countV`
//! of them evenly over the parameters from `low` to `high` that are `inside`, to the nearest
//! place of the finer.
double apart(const occ::handle<Geom_Surface>&         coarser,
             const occ::handle<Geom_Surface>&         finer,
             int                                      countU,
             int                                      countV,
             const gp_XY&                             low,
             const gp_XY&                             high,
             const std::function<bool(const gp_XY&)>& inside);

//! A finer net is not tried when the one before it brought the worst miss down by less than
//! to this share of what it was: the miss has then stopped falling, and what is left of the
//! ladder would be spent on a net that comes no nearer.
constexpr double stillFalling = 0.75;

//! How far one side turns into the next at their corner, in degrees: none where it goes straight
//! on, 180 where it turns back.
double turning(const Side& from, const Side& to);

//! Whether four sides meet at corners, so that they can be the edges of one surface.
bool squarish(const std::vector<Side>& sides, double g1Degrees);

//! One side of a quad: an edge, or edges in loop order that meet without a corner.
struct Run
{
  std::vector<const Side*> parts;
  //! Where each part starts, as a share of the run's length.
  std::vector<double> starts;
  double              length = 0;

  //! The part at a share of the run, and how far along that part.
  std::pair<const Side*, double> at(double fraction) const
  {
    size_t k = 0;
    while (k + 1 < parts.size() && fraction > starts[k + 1])
      ++k;
    const double end = k + 1 < parts.size() ? starts[k + 1] : 1;
    return {parts[k], std::clamp((fraction - starts[k]) / (end - starts[k]), 0.0, 1.0)};
  }

  gp_Pnt point(double fraction) const
  {
    const auto [part, f] = at(fraction);
    return part->at(f);
  }

  gp_Vec tangent(double fraction) const
  {
    const auto [part, f] = at(fraction);
    return part->tangent(f);
  }
};

//! The four sides of a quad, the first holding the first edge. Four edges are four sides.
//! Of more edges, those that meet without a corner (within the angular tolerance) make one
//! side: there must be four corners then.
std::optional<std::vector<Run>> quadRuns(const std::vector<Side>& sides, double g1Degrees);

//! Whether two edges are certainly farther apart than `gap`, from the pole boxes of their
//! Bezier pieces: positive weights keep each piece inside its poles. Pieces whose boxes come
//! nearer than `gap` are halved, the longer first, a few times over; false says only that it
//! could not be shown so, and the edges are then measured.
bool apart(const TopoDS_Edge& a, const TopoDS_Edge& b, double gap);

std::optional<Built> foursided(const Asked& asked, const Frame& frame);

//! An edge of a trimmed face: the curve of a side in the surface's parameters, run by the
//! side's own parameter, and the furthest the surface is there from the side at that parameter.
using patch_side::Rim;

//! Where a side lies on a surface: for every parameter of the side's curve, the parameters of
//! the surface's point nearest to the curve's. A face whose edge has this as its curve in the
//! surface's parameters is as near its edge at every parameter as its surface is to the curve
//! at all, which is what continuity across the edge is read by.
//!
//! It starts as the side laid into the mean plane, which is where the fit asked the surface to
//! be, and is then moved onto the true places pole by pole (each pole by what its own place on
//! the curve is off: progressive iteration, which converges for a B-spline's basis), with knots
//! added while the curve between those places is off by more than a share of the tolerance.
//! `parameterLayout` gives the place for a curve parameter in a chart that is not linear in
//! space: the curve then starts from those places, not from its poles laid into the chart.
Rim preimage(const occ::handle<Geom_BSplineSurface>&    surface,
             const Side&                                side,
             const Frame&                               frame,
             double                                     tolerance,
             double                                     size,
             const std::function<gp_XY(const gp_Pnt&)>& layout          = {},
             const std::function<gp_XY(double)>&        parameterLayout = {});

using patch_side::closeGaps;

//! The face of `surface` inside the loop. Each edge is the side's own curve, with `rims` as
//! its curve in the surface's parameters.
TopoDS_Face trimmed(const occ::handle<Geom_BSplineSurface>& surface,
                    const Asked&                            asked,
                    const std::vector<Rim>&                 rims,
                    bool                                    clockwise);

//! Check the p-curves that sewing actually delivers, including its same-parameter repair.
void verifySewn(Built& built, const Asked& asked);

//! Trim an exact continuation and measure its delivered parameter curves, rather than
//! merely measuring the untrimmed surface at nearest points.
std::optional<Built> continued(const occ::handle<Geom_BSplineSurface>&    surface,
                               const Asked&                               asked,
                               const Frame&                               frame,
                               const std::function<gp_XY(const gp_Pnt&)>& layout);

//! A G0 loop may have four logical sides split into tangent edge chains. Seed a
//! Coons surface from those chains, including degrees above the fitter's limit.
//! `plain` bounds rational approximation and Coons averages mismatched corner ends;
//! every original edge still supplies its own rim and measured constraint. This path
//! never calls fitNet and only follows a regular fit that has missed its limits.
std::optional<Built> chainedQuad(const Asked& asked, const Frame& frame);

//! If every edge supplies the very same spline surface, its continuation is a candidate
//! before a variational fill. This uses the supplied geometry, never a guessed interior
//! reconstructed from G0 curves. Different supports and explicit shape controls still fit.
std::optional<Built> splineContinuation(const Asked& asked, const Frame&);

//! A toroidal neighbour has an exact rational continuation too. Trim the analytic torus
//! in its angular parameters before conversion, so the rational rectangle has no seam
//! inside the requested window. Its B-spline parameters differ from the angles.
std::optional<Built> toroidal(const Asked& asked, const Frame& frame);

//! Continue a spherical neighbour exactly, using stereographic parameters. Unlike longitude
//! and latitude this rational biquadratic rectangle has no pole or seam inside the cap.
//! Only a rim wholly on the same ball qualifies; the delivered joint is still measured.
std::optional<Built> spherical(const Asked& asked, const Frame& flat);

//! Continue two tangent planes and the cylinder between them as one ruled NURBS.
//! Its only curvature breaks are the two existing plane/round joins, along the round's axis.
std::optional<Built> ruledRound(const Asked& asked, const Frame& frame);

//! How far the surface reaches past the loop, as a share of the loop's extent that way: what is
//! asked along the loop is then asked inside the surface and not on its own edge.
constexpr double marginShare = 0.1;

//! The spans across the longer extent of the loop, net by net, each about 1.6 times the one
//! before: a miss that falls with the fourth power of the span falls sixfold from one to the
//! next.
constexpr int ladder[] = {6, 10, 16, 24, 36, 54};

//! The fewest spans a direction is given, however short the loop is that way.
constexpr int fewestSpans = 3;

//! How near the curves a fit is asked to come where they are only met, as a share of the
//! loop's size: about as near as a fair surface over the coarsest net comes (z = sin 2a on a
//! radius of 6: 1.7e-5 over six spans, which is 1.4e-6 of the size 12).
constexpr double fitReach = 1e-4;

//! How near the curves the surface is moved after the fit, as a share of the tolerance: half,
//! which leaves the other half to what the curve of an edge in the surface's parameters adds
//! and to the places between those that are looked at.
constexpr double rimShare = 0.5;

//! A surface that is within this share of the tolerance at every place is left as it is: it
//! is within the limit wherever it is read, and moved, its net would be made finer for nothing
//! (a window in the inside of a torus went from 10 x 11 to 17 x 17 poles).
constexpr double enoughShare = 0.9;

//! The places along a curve at which the surface is moved onto it: this many to a span of the
//! net as fitted, and no fewer than this many to a curve (64 short sides with 96 each took
//! minutes).
constexpr int rimPlacesPerSpan = 6, fewestRimPlaces = 16, mostRimPlaces = 4096;

//! The most spans either way of a net whose surface is moved onto its curves, and how many
//! times finer than it was fitted a net may become that goes over into a face: there the
//! fit has come near already, and what is left is what its passes did not reach; and where
//! curves are only met, where the fit was not asked for the tolerance (a hexagon round a
//! cube's corner asked for 1e-9 took 26 s on the way to the most).
constexpr int mostRimSpans = 192, rimFiner = 2, rimFinerMet = 4;

//! Places along a side per span of the net that it crosses. Three met everything with the
//! fewest poles and rippled; two is what was kept.
constexpr int placesPerSpan = 2;

//! The most poles of a net that is free in space. A net that is to go over smoothly into a
//! face does better so than as a height over the mean plane: its parameters can then run with
//! the surface where it is steep, and the energy, which is of the parameters, is nearer to
//! one of the shape (a ball's cap cut at 6 of 10 comes out 0.13 off the ball, as a height 1.3).
//! It has three times the unknowns and takes 27 times the work: 25 x 25 poles free are about
//! what a height of 57 x 57 takes.
constexpr int freePoles = 640;

//! What a surface that meets its edges and is not one to keep counts as missing by, in units
//! of the limits: no more than any that misses an edge, so that where no net gives a patch it
//! is said that the edges can be met and what is wrong with the surface that meets them.
constexpr double unkeptExcess = 1;

//! The nets that may meet the edges with a surface that is not one to keep before it is given
//! up: the two coarsest and one more.
constexpr int mostUnkept = 3;

//! Within this share of the miss before, a miss counts as the same: what is in the way is then
//! not the net.
constexpr double same = 0.02;

//! The coarsest nets, whose misses say nothing yet of how the finer ones will do.
constexpr int unruly = 2;

//! The places along a side at which its face's normal is read for whether any plane is one the
//! fill could be a height over (see graphMargin), its ends included.
constexpr int gaussPlaces = 64;

//! A fill over a chart of the loop's own, for a loop whose rim's faces stand so that it is a
//! height over no plane (see overPlane). The loop is spanned by a disc of triangles that leaves
//! every side as its face asks and is laid flat conformally (proxyChart); a net free in space is
//! fitted over that chart, held to the sides at their places in it as over the mean plane: with
//! its tangent plane the face's and its bend across the face's where they are asked, which are
//! conditions on the derivatives in space and none on how the chart lies. It is trimmed by the
//! sides' curves in its parameters and kept only where S_u x S_v is shown to keep away from zero
//! inside the loop (regularMargin). Nothing where no net over the chart meets what is asked.
std::optional<Built> overChart(const Asked& asked, const Frame& frame);

Built overPlane(const Asked& asked, const Frame& flat, int direction = 0);

//! Next/Adjacent with Quad sided off. Over a conformal chart of a surface, the
//! surface's U and V meet at right angles; so the Natural fill is laid flat conformally, and
//! a new surface is fitted over that chart: one that meets the same edges as asked and keeps
//! near the Natural fill inside, as a wish and not a limit. How square U and V
//! come out is read at every node and in every triangle and reported (see Built::squareness).
Built adjacent(const Asked& asked, const Frame& frame, const Built& natural);

} // namespace BRepExt::patch_fill
