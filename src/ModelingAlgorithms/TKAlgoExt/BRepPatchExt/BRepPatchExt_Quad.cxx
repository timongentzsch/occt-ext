// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Four sides: an untrimmed surface with the sides as its own edges, the Coons patch of the four
// curves, its inner poles then moved to meet G1 and G2; and a quad chained from runs of sides.

#include <BRepPatchExt_Fill.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>

#include <GeomAlgoExt_Cancel.hxx>
#include <GeomAlgoExt_PatchRim.hxx>
#include <GeomAlgoExt_Phase.hxx>

#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BSplCLib.hxx>
#include <Bnd_Box.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomConvert.hxx>
#include <GeomConvert_ApproxCurve.hxx>
#include <GeomConvert_BSplineCurveToBezierCurve.hxx>
#include <GeomConvert_CompCurveToBSplineCurve.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_Line.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <ShapeAnalysis_Surface.hxx>

#include <array>

namespace BRepExt::patch_fill
{
namespace
{

//! Where two guides come nearest is narrowed down until the step of their parameters, which
//! run from 0 to 1, is below this: some thousands of roundings of a double. Chosen.
constexpr double finestStep = 1e-12;

//! What share of the tolerance a side that has no polynomial form (an arc, a rational spline)
//! may be off the spline that stands for it in a four-sided surface, and the most pieces of
//! degree 5 that spline is given: a half circle is within 1e-8 of its radius with 16.
constexpr double plainShare = 0.1;

constexpr int plainSegments = 64;

std::vector<double> flatKnots(const occ::handle<Geom_BSplineCurve>& curve)
{
  const NCollection_Array1<double>& sequence = curve->KnotSequence();
  return {sequence.begin(), sequence.end()};
}

//! The Coons patch of four curves that agree in degree and knots pairwise: `bottom` and `top`
//! along u, `left` and `right` along v. A linear blend is exact in a B-spline basis when taken
//! at the Greville abscissae.
SplineNet coons(const occ::handle<Geom_BSplineCurve>& bottom,
                const occ::handle<Geom_BSplineCurve>& right,
                const occ::handle<Geom_BSplineCurve>& top,
                const occ::handle<Geom_BSplineCurve>& left)
{
  SplineNet net{bottom->Degree(), left->Degree(), flatKnots(bottom), flatKnots(left), {}};
  const int nu = net.countU(), nv = net.countV();
  auto      greville = [](const std::vector<double>& knots, int degree, int i) {
    double sum = 0;
    for (int k = 1; k <= degree; ++k)
      sum += knots[i + k];
    return sum / degree / knots.back();
  };
  auto         mid = [](const gp_Pnt& a, const gp_Pnt& b) { return (a.XYZ() + b.XYZ()) / 2; };
  const gp_XYZ p00 = mid(bottom->Pole(1), left->Pole(1)),
               p10 = mid(bottom->Pole(nu), right->Pole(1)), p01 = mid(top->Pole(1), left->Pole(nv)),
               p11 = mid(top->Pole(nu), right->Pole(nv));
  net.poles.resize(static_cast<size_t>(nu) * nv);
  for (int i = 0; i < nu; ++i)
    for (int j = 0; j < nv; ++j)
    {
      const double a = greville(net.knotsU, net.degreeU, i),
                   b = greville(net.knotsV, net.degreeV, j);
      net.poles[static_cast<size_t>(i) * nv + j] =
        bottom->Pole(i + 1).XYZ() * (1 - b) + top->Pole(i + 1).XYZ() * b
        + left->Pole(j + 1).XYZ() * (1 - a) + right->Pole(j + 1).XYZ() * a
        - (p00 * ((1 - a) * (1 - b)) + p10 * (a * (1 - b)) + p01 * ((1 - a) * b) + p11 * (a * b));
    }
  return net;
}

//! The sharpest and the widest corner of a four-sided surface, as the cosine of what the
//! corner lacks of a straight line or of a turn back: 15 and 165 degrees. A surface's two
//! directions stand at the corner's angle to each other there, and at 15 degrees its area is
//! still a quarter of what a right angle gives it (the sine); nearer to a line than that, the
//! sides are better taken as one loop and trimmed.
constexpr double sharpestCorner = 0.966;

//! Whether two sides meet at a corner of a quad: one that turns by more than the angular
//! tolerance, and turns back by less: four curves make a quad at turns from under 4 to over 176
//! degrees.
bool cornered(const Side& from, const Side& to, double g1Degrees)
{
  const double angle = turning(from, to);
  return angle > g1Degrees && angle < 180 - g1Degrees;
}

//! Places along a side of a four-sided net per pole of it that way: four, so that between
//! two knots of the border, which the surface follows exactly, the tangent plane is asked
//! more often than the net has freedom to wiggle between.
constexpr int placesPerPole = 4;

//! The most poles each way that the Coons patch of four sides is built with: the sides' own,
//! which a curve with some hundred poles of its own may ask for.
constexpr int mostCoonsPoles = 1024;

//! The most poles each way of a four-sided net whose inner poles are fitted: about what the
//! finest net of a trimmed surface has.
constexpr int mostFittedPoles = 64;

} // namespace

occ::handle<Geom_BSplineCurve> plain(const Side& side, double length, bool turned, double tolerance)
{
  const double            first = std::min(side.from, side.to), last = std::max(side.from, side.to);
  occ::handle<Geom_Curve> trimmed = new Geom_TrimmedCurve(side.curve, first, last);
  occ::handle<Geom_BSplineCurve> curve;
  const auto                     spline = occ::down_cast<Geom_BSplineCurve>(side.curve);
  const auto                     bezier = occ::down_cast<Geom_BezierCurve>(side.curve);
  if (side.curve->IsKind(STANDARD_TYPE(Geom_Line)) || (!spline.IsNull() && !spline->IsRational())
      || (!bezier.IsNull() && !bezier->IsRational()))
  {
    curve = GeomConvert::CurveToBSplineCurve(trimmed);
  }
  else
  {
    // Arcs and rational splines: their exact B-spline has weights, or other parameters.
    GeomConvert_ApproxCurve approximation(trimmed,
                                          tolerance * plainShare,
                                          GeomAbs_C2,
                                          plainSegments,
                                          5);
    if (!approximation.HasResult() || approximation.MaxError() > tolerance * plainShare)
      return {};
    curve = approximation.Curve();
  }
  if (curve->IsPeriodic())
    curve->SetNotPeriodic();
  if ((side.to < side.from) != turned)
    curve->Reverse();
  NCollection_Array1<double> knots = curve->Knots();
  BSplCLib::Reparametrize(0, length, knots);
  curve->SetKnots(knots);
  return curve;
}

double apart(const occ::handle<Geom_Surface>&         coarser,
             const occ::handle<Geom_Surface>&         finer,
             int                                      countU,
             int                                      countV,
             const gp_XY&                             low,
             const gp_XY&                             high,
             const std::function<bool(const gp_XY&)>& inside)
{
  double farthest = 0;
  for (int i = 0; i < countU; ++i)
    for (int j = 0; j < countV; ++j)
    {
      const gp_XY at(low.X() + (high.X() - low.X()) * (i + 0.5) / countU,
                     low.Y() + (high.Y() - low.Y()) * (j + 0.5) / countV);
      if (!inside(at))
        continue;
      const gp_Pnt point = coarser->Value(at.X(), at.Y());
      const gp_XY  uv    = nearest(finer, at, point);
      farthest           = std::max(farthest, point.Distance(finer->Value(uv.X(), uv.Y())));
    }
  return farthest;
}

double turning(const Side& from, const Side& to)
{
  return std::acos(
           std::clamp(from.tangent(1).Normalized().Dot(to.tangent(0).Normalized()), -1.0, 1.0))
         * 180 / M_PI;
}

bool squarish(const std::vector<Side>& sides, double g1Degrees)
{
  if (sides.size() != 4)
    return false;
  for (size_t i = 0; i < 4; ++i)
    if (!cornered(sides[i], sides[(i + 1) % 4], g1Degrees))
      return false;
  return true;
}

std::optional<std::vector<Run>> quadRuns(const std::vector<Side>& sides, double g1Degrees)
{
  const size_t n = sides.size();
  if (n < 4)
    return std::nullopt;
  std::vector<size_t> corners;
  for (size_t i = 0; i < n; ++i)
  {
    const Side &from = sides[i], &to = sides[(i + 1) % n];
    if (n == 4 ? !cornered(from, to, g1Degrees) : turning(from, to) >= 180 - g1Degrees)
      return std::nullopt;
    if (n == 4 || turning(from, to) > g1Degrees)
      corners.push_back((i + 1) % n);
  }
  if (corners.size() != 4)
    return std::nullopt;
  std::sort(corners.begin(), corners.end());
  // The side that holds the first edge is the first: the last when that edge goes on from it.
  std::rotate(corners.begin(), corners.begin() + (corners.front() == 0 ? 0 : 3), corners.end());
  std::vector<Run> runs(4);
  for (size_t r = 0; r < 4; ++r)
  {
    Run& run = runs[r];
    for (size_t i = corners[r]; run.parts.empty() || i != corners[(r + 1) % 4]; i = (i + 1) % n)
    {
      run.starts.push_back(run.length);
      run.parts.push_back(&sides[i]);
      run.length += sides[i].length;
    }
    for (double& start : run.starts)
      start /= run.length;
  }
  return runs;
}

bool apart(const TopoDS_Edge& a, const TopoDS_Edge& b, double gap)
{
  const auto pieces = [](const TopoDS_Edge& edge) {
    std::vector<occ::handle<Geom_BezierCurve>> out;
    double                                     first, last;
    TopLoc_Location                            location;
    const auto curve = BRep_Tool::Curve(edge, location, first, last);
    if (curve.IsNull())
      return out;
    auto spline = GeomConvert::CurveToBSplineCurve(new Geom_TrimmedCurve(curve, first, last));
    if (!location.IsIdentity())
      spline->Transform(location.Transformation());
    GeomConvert_BSplineCurveToBezierCurve convert(spline);
    for (int i = 1; i <= convert.NbArcs(); ++i)
      out.push_back(convert.Arc(i));
    return out;
  };
  const auto box = [](const occ::handle<Geom_BezierCurve>& piece) {
    Bnd_Box out;
    for (int i = 1; i <= piece->NbPoles(); ++i)
      out.Add(piece->Pole(i));
    return out;
  };
  const auto halves = [](const occ::handle<Geom_BezierCurve>& piece) {
    std::array<occ::handle<Geom_BezierCurve>, 2> out;
    for (int half = 0; half < 2; ++half)
    {
      out[half] = occ::down_cast<Geom_BezierCurve>(piece->Copy());
      out[half]->Segment(half * 0.5, (half + 1) * 0.5);
    }
    return out;
  };

  struct Pair
  {
    occ::handle<Geom_BezierCurve> a, b;
    int                           depth;
  };

  std::vector<Pair> pairs;
  const auto        ofA = pieces(a), ofB = pieces(b);
  if (ofA.empty() || ofB.empty())
    return false;
  for (const auto& x : ofA)
    for (const auto& y : ofB)
      pairs.push_back({x, y, 0});
  for (int visits = 0; !pairs.empty(); ++visits)
  {
    const Pair pair = pairs.back();
    pairs.pop_back();
    const Bnd_Box boxA = box(pair.a), boxB = box(pair.b);
    if (boxA.Distance(boxB) > gap)
      continue;
    if (pair.depth >= 12 || visits >= 4096)
      return false;
    const bool first = boxA.SquareExtent() >= boxB.SquareExtent();
    for (const auto& half : halves(first ? pair.a : pair.b))
      pairs.push_back(first ? Pair{half, pair.b, pair.depth + 1}
                            : Pair{pair.a, half, pair.depth + 1});
  }
  return true;
}

namespace
{

// A guide from side to side runs straight across the rectangle, through the place of
// each guide it crosses: the nearest point of a surface that does not yet rise to it
// would drag the parameters along with the rise. Others go by that nearest point.
struct Across
{
  bool                                  laid = false;
  std::vector<std::pair<double, gp_XY>> stops;
};

//! The place of a laid guide at `f`, along the straight pieces between its `stops`.
gp_XY laidAt(std::vector<std::pair<double, gp_XY>> stops, double f)
{
  std::sort(stops.begin(), stops.end(), [](const auto& a, const auto& b) {
    return a.first < b.first;
  });
  size_t i = 0;
  while (i + 2 < stops.size() && f > stops[i + 1].first)
    ++i;
  const double t = (f - stops[i].first) / std::max(stops[i + 1].first - stops[i].first, vanishing);
  return gp_XY(stops[i].second + (stops[i + 1].second - stops[i].second) * t);
}

//! The even places along two guides at which where they come nearest is first looked for,
//! before it is narrowed down. Chosen: about as many as an edge is measured at (`places`).
constexpr int crossingPlaces = 129;

//! Four sides: the runs of the loop as the four sides of an untrimmed surface, the Coons patch
//! of their curves, its inner poles then moved to meet G1 and G2, net by finer net.
class QuadFill
{
public:
  QuadFill(const Asked& asked, const Frame& frame, std::vector<Run> found)
      : asked(asked),
        frame(frame),
        runs(std::move(found))
  {
    // A side made only of gaps of an open boundary: nothing is asked along it.
    for (size_t s = 0; s < 4; ++s)
      gap[s] = std::all_of(runs[s].parts.begin(), runs[s].parts.end(), [](const Side* part) {
        return part->free;
      });
    for (int pair = 0; pair < 2; ++pair)
      for (const Side* a : runs[pair].parts)
        for (const Side* b : runs[pair + 2].parts)
          if (!a->free && !b->free && !apart(a->edge, b->edge, asked.limits.g0)
              && BRepExtrema_DistShapeShape(a->edge, b->edge).Value() <= asked.limits.g0)
            fail(op, invalid_input, "the boundary crosses itself");
    for (const Side& side : asked.sides)
      order = std::max(order, side.order);
    fitting = order > 0 || asked.flatness >= 0 || asked.controls.flow > PatchFlow::Natural
              || !asked.guides.empty() || !asked.points.empty();
    lu = (runs[0].length + runs[2].length) / 2;
    lv = (runs[1].length + runs[3].length) / 2;
  }

  //! The surface, or nothing where the sides cannot be its edges.
  std::optional<Built> run();

private:
  //! Each run as one curve, round the loop; false where an edge cannot be made a plain curve.
  bool joinSides();

  //! Opposite curves brought to one degree and one set of knots; false where they cannot be.
  bool matchSides();

  //! A place of a side is on the surface's own edge, as far along it.
  gp_XY along(size_t s, double f) const;

  //! Where an edge is in its side, for a place of it to be found on the surface.
  gp_XY place(const Side& side, double f) const;

  //! The way the isocurves across a side are to leave it, as the one direction in the
  //! tangent plane that the derivative across is to have nothing of: square to the side for
  //! Normal. For Adjacent they turn from the side before to the side after, which are those
  //! isocurves at the corners; for Next they head for the same share of the opposite side,
  //! turned at the corners into those same sides. Square to the side on a rectangle, both.
  gp_Vec axis(size_t s, double f) const;

  //! The nets, from the Coons patch of the curves, each with twice the spans of the one before.
  std::optional<Built> rounds();

  //! The poles a fit keeps: those of the sides that are not gaps, and the corners.
  std::vector<char> fixedPoles(const SplineNet& net, int nu, int nv) const;

  //! What a fit is held to along the sides, `placesPerPole` places to a pole.
  std::vector<Target> sideTargets(int nu, int nv) const;

  //! What a fit is held to along the guides and at the points inside.
  void interiorTargets(std::vector<Target>&                    targets,
                       const SplineNet&                        net,
                       const std::vector<char>&                fixed,
                       const occ::handle<Geom_BSplineSurface>& surface,
                       int                                     nu,
                       int                                     nv) const;

  //! Whether a place in the parameters is on the rectangle's edge, within the tolerance.
  bool onSide(const gp_XY& uv) const;

  //! The guides laid straight across the rectangle where they run from side to side.
  std::vector<Across> layAcross(ShapeAnalysis_Surface& projection) const;

  //! Two laid guides that cross pass through one place of the rectangle.
  void crossLaid(std::vector<Across>& across) const;

  //! A guide that ends on one so laid, at a side or at another such guide, is laid to there.
  void layEnds(std::vector<Across>& across, ShapeAnalysis_Surface& projection) const;

  //! The face of `surface`, measured along the sides and for the flow across them.
  Built measured(const occ::handle<Geom_BSplineSurface>& surface) const;

  //! Twice the spans for the next round.
  void finer();

  const Asked&                   asked;
  const Frame&                   frame;
  std::vector<Run>               runs;
  bool                           gap[4];
  int                            order   = 0;
  bool                           fitting = false;
  double                         lu = 0, lv = 0;
  occ::handle<Geom_BSplineCurve> curves[4];
};

std::optional<Built> QuadFill::run()
{
  if (!joinSides() || !matchSides())
    return {};
  // Quality budgets are minimum span counts.
  const int minimum = asked.controls.quality == PatchQuality::Max    ? 16
                      : asked.controls.quality == PatchQuality::High ? 10
                                                                     : 1;
  for (int i = 0; i < 2; ++i)
    for (int k = 1; k < minimum; ++k)
      for (int c : {i, i + 2})
        curves[c]->InsertKnot((i ? lv : lu) * k / minimum);
  return rounds();
}

bool QuadFill::joinSides()
{
  for (int i = 0; i < 4; ++i)
  {
    // Round the loop the top and the left run backwards.
    if (runs[i].parts.size() == 1)
    {
      curves[i] = plain(*runs[i].parts[0], i % 2 ? lv : lu, i >= 2, asked.limits.g0);
      if (curves[i].IsNull())
        return false;
      continue;
    }
    // The edges joined end to end, each scaled to go on at the speed the one before ends
    // with: where they meet without a corner, the curve is then once differentiable there
    // and so is the net across, which a fit would otherwise be free to crease. Each edge
    // keeps its own parameter within its stretch, which `starts` then gives.
    std::optional<GeomConvert_CompCurveToBSplineCurve> joined;
    std::vector<double>                                junctions;
    for (const Side* part : runs[i].parts)
    {
      const auto piece = plain(*part, part->length, false, asked.limits.g0);
      if (piece.IsNull())
        return false;
      if (!joined)
      {
        joined.emplace(piece);
        continue;
      }
      junctions.push_back(joined->BSplineCurve()->LastParameter());
      if (!joined->Add(piece, asked.limits.g0, true, true))
        return false;
    }
    curves[i]          = joined->BSplineCurve();
    const double first = curves[i]->FirstParameter(), span = curves[i]->LastParameter() - first;
    for (size_t k = 0; k < junctions.size(); ++k)
    {
      runs[i].starts[k + 1] = (junctions[k] - first) / span;
      for (int knot = 2; knot < curves[i]->NbKnots(); ++knot)
        if (std::abs(curves[i]->Knot(knot) - junctions[k])
            <= BRepExt::tolerance::parameterMerge * span)
        {
          // Kept as it is where the edges do not go on from each other closely enough.
          curves[i]->RemoveKnot(knot, curves[i]->Degree() - 1, asked.limits.g0 * plainShare);
          break;
        }
    }
    if (i >= 2)
      curves[i]->Reverse();
    NCollection_Array1<double> knots = curves[i]->Knots();
    BSplCLib::Reparametrize(0, i % 2 ? lv : lu, knots);
    curves[i]->SetKnots(knots);
  }
  return true;
}

bool QuadFill::matchSides()
{
  for (int i = 0; i < 2; ++i)
  {
    auto &a = curves[i], &b = curves[i + 2];
    // Asked for a quad, a side above the net's highest degree, or one with more poles at the
    // net's degree than a fitted net may have, is brought down to it as an arc is, within the
    // same share of the tolerance and over the same parameters.
    const int least = !fitting ? 1 : order < 2 ? 3 : 5;
    for (auto* curve : {&a, &b})
      if (asked.controls.quad == PatchQuad::On
          && ((*curve)->Degree() > 5
              || (fitting
                  && (*curve)->NbPoles()
                         + (std::max(least, (*curve)->Degree()) - (*curve)->Degree())
                             * ((*curve)->NbKnots() - 1)
                       > mostFittedPoles)))
      {
        GeomConvert_ApproxCurve lower(*curve,
                                      asked.limits.g0 * plainShare,
                                      GeomAbs_C2,
                                      plainSegments,
                                      5);
        if (!lower.HasResult() || lower.MaxError() > asked.limits.g0 * plainShare)
          return false;
        *curve = lower.Curve();
      }
    const int degree = std::max({a->Degree(), b->Degree(), least});
    if (degree > 5)
      return false;
    a->IncreaseDegree(degree);
    b->IncreaseDegree(degree);
    const double same = BRepExt::tolerance::parameterMerge * (i ? lv : lu);
    a->InsertKnots(b->Knots(), b->Multiplicities(), same, false);
    b->InsertKnots(a->Knots(), a->Multiplicities(), same, false);
    if (a->NbPoles() != b->NbPoles())
      return false;
  }
  return true;
}

gp_XY QuadFill::along(size_t s, double f) const
{
  return s == 0   ? gp_XY(f * lu, 0)
         : s == 1 ? gp_XY(lu, f * lv)
         : s == 2 ? gp_XY((1 - f) * lu, lv)
                  : gp_XY(0, (1 - f) * lv);
}

gp_XY QuadFill::place(const Side& side, double f) const
{
  for (size_t s = 0; s < 4; ++s)
    for (size_t k = 0; k < runs[s].parts.size(); ++k)
      if (runs[s].parts[k] == &side)
      {
        const double end = k + 1 < runs[s].parts.size() ? runs[s].starts[k + 1] : 1;
        return along(s, runs[s].starts[k] + f * (end - runs[s].starts[k]));
      }
  return gp_XY();
}

gp_Vec QuadFill::axis(size_t s, double f) const
{
  const auto [part, g] = runs[s].at(f);
  if (asked.controls.flow == PatchFlow::Normal)
    return part->tangent(g);
  const gp_Vec tangent = part->tangent(g).Normalized();
  gp_Vec       before  = -runs[(s + 3) % 4].tangent(1).Normalized(),
         after         = runs[(s + 1) % 4].tangent(0).Normalized();
  const auto chord     = [&](double t) {
    return gp_Vec(runs[s].point(t), runs[(s + 2) % 4].point(1 - t)).Normalized();
  };
  // A gap is no side to turn to: next to one the isocurves leave as for Normal, square to
  // the side, for Adjacent, and along the chord for Next.
  const bool fromGap = gap[(s + 3) % 4], toGap = gap[(s + 1) % 4];
  if (asked.controls.flow == PatchFlow::Adjacent && fromGap && toGap)
    return part->tangent(g);
  const auto square = [&](const gp_Vec& way, double t) {
    const gp_Vec along = runs[s].tangent(t).Normalized(), out = way - along * way.Dot(along);
    return out.Magnitude() > vanishing ? out.Normalized() : out;
  };
  if (fromGap)
    before = asked.controls.flow == PatchFlow::Next ? chord(0) : square(after, 0);
  if (toGap)
    after = asked.controls.flow == PatchFlow::Next ? chord(1) : square(before, 1);
  gp_Vec way = before * (1 - f) + after * f;
  if (asked.controls.flow == PatchFlow::Next)
    way = chord(f) + (before - chord(0)) * (1 - f) + (after - chord(1)) * f;
  if (part->order > 0)
  {
    const gp_Vec n = part->jet(g).n;
    way -= n * way.Dot(n);
  }
  if (way.Magnitude() <= BRepExt::tolerance::exactLength * asked.size
      || std::abs(way.Normalized().Dot(tangent)) > sharpestCorner)
    fail(op,
         unsupported_geometry,
         "the UV flow asked for runs along boundary edge " + std::to_string(part->input + 1));
  way.Normalize();
  return tangent - way * tangent.Dot(way);
}

std::optional<Built> QuadFill::rounds()
{
  std::optional<Built> best;
  for (int round = 0; round < 5; ++round)
  {
    const int most = fitting ? mostFittedPoles : mostCoonsPoles;
    if (curves[0]->NbPoles() > most || curves[1]->NbPoles() > most)
      break;
    SplineNet                        net = coons(curves[0], curves[1], curves[2], curves[3]);
    const int                        nu = net.countU(), nv = net.countV();
    occ::handle<Geom_BSplineSurface> surface = surfaceOf(net);
    if (fitting)
    {
      // The poles of a side that is a gap are free like the inner ones, but for the corners.
      const std::vector<char> fixed   = fixedPoles(net, nu, nv);
      std::vector<Target>     targets = sideTargets(nu, nv);
      interiorTargets(targets, net, fixed, surface, nu, nv);
      // As for supported trimmed fits, solve the coordinates near the loop's origin.
      // Subtracting the energy of large world-coordinate poles loses the small interior
      // correction, even when the fixed boundary still passes its continuity checks.
      surface = fit(net,
                    fixed,
                    targets,
                    asked.size,
                    asked.limits,
                    surface,
                    {},
                    nullptr,
                    {},
                    &frame,
                    asked.flatness);
    }
    Built        built  = measured(surface);
    const double before = best ? best->measured.excess : HUGE_VAL;
    if (built.measured.excess < before)
      best = built;
    // A finer net that comes no nearer is the last, as for a trimmed surface.
    if (built.measured.within || !fitting || built.measured.excess > stillFalling * before)
      break;
    // Twice the spans for the next round.
    finer();
  }
  return best;
}

std::vector<char> QuadFill::fixedPoles(const SplineNet& net, int nu, int nv) const
{
  std::vector<char> fixed(net.poles.size(), 0);
  for (int i = 0; i < nu; ++i)
    for (int j = 0; j < nv; ++j)
      fixed[static_cast<size_t>(i) * nv + j] = ((i == 0 || i == nu - 1) && (j == 0 || j == nv - 1))
                                               || (j == 0 && !gap[0]) || (i == nu - 1 && !gap[1])
                                               || (j == nv - 1 && !gap[2]) || (i == 0 && !gap[3]);
  return fixed;
}

std::vector<Target> QuadFill::sideTargets(int nu, int nv) const
{
  std::vector<Target> targets;
  for (size_t s = 0; s < 4; ++s)
  {
    const int count = placesPerPole * (s % 2 ? nv : nu);
    for (int k = 0; k < count; ++k)
    {
      const double f       = (k + 0.5) / count;
      const auto [part, g] = runs[s].at(f);
      if (part->free || (part->order < 1 && asked.controls.flow == PatchFlow::Natural))
        continue;
      Target target{part->at(g),
                    along(s, f),
                    part->order,
                    part->order ? std::optional<SurfaceJet>(part->jet(g)) : std::nullopt,
                    s % 2 ? gp_XY(1, 0) : gp_XY(0, 1),
                    runs[s].length / count / asked.size};
      if (asked.controls.flow != PatchFlow::Natural)
        target.perpendicular = axis(s, f);
      if (part->order)
        target.easing = part->easing(g);
      targets.push_back(target);
    }
  }
  return targets;
}

void QuadFill::interiorTargets(std::vector<Target>&                    targets,
                               const SplineNet&                        net,
                               const std::vector<char>&                fixed,
                               const occ::handle<Geom_BSplineSurface>& surface,
                               int                                     nu,
                               int                                     nv) const
{
  ShapeAnalysis_Surface projection(surface);
  const auto            guideTarget = [&](const gp_Pnt& point, double share) {
    gp_XY uv = projection.ValueOfUV(point, BRepExt::tolerance::exactLength).XY();
    // A guide that ends on a side comes out on it to within the tolerance, either way:
    // the parameters here are lengths.
    const double slack = asked.limits.g0;
    if (uv.X() < -slack || uv.X() > lu + slack || uv.Y() < -slack || uv.Y() > lv + slack)
      fail(op, invalid_input, "an interior constraint lies outside the quad");
    uv.SetCoord(std::clamp(uv.X(), 0.0, lu), std::clamp(uv.Y(), 0.0, lv));
    targets.push_back({point, uv, 0, {}, {}, share});
  };
  const std::vector<Across> across = layAcross(projection);
  for (size_t g = 0; g < asked.guides.size(); ++g)
  {
    const Side& guide = asked.guides[g];
    const int   count = std::max(129, placesPerPole * std::max(nu, nv));
    for (int k = 0; k <= count; ++k)
    {
      const double f = static_cast<double>(k) / count, share = guide.length / count / asked.size;
      if (across[g].laid)
        targets.push_back({guide.at(f), laidAt(across[g].stops, f), 0, {}, {}, share});
      else
        guideTarget(guide.at(f), share);
    }
  }
  // A point goes where the surface that holds the laid guides passes nearest it, for the
  // same reason: the rectangle's own surface does not rise to it.
  const bool anyLaid =
    std::any_of(across.begin(), across.end(), [](const Across& a) { return a.laid; });
  if (anyLaid && !asked.points.empty())
  {
    SplineNet             first = net;
    const auto            held  = fit(first,
                          fixed,
                          targets,
                          asked.size,
                          asked.limits,
                          surface,
                                      {},
                          nullptr,
                                      {},
                          &frame,
                          asked.flatness);
    ShapeAnalysis_Surface nearestOn(held);
    for (const gp_Pnt& point : asked.points)
    {
      const gp_XY uv = nearestOn.ValueOfUV(point, BRepExt::tolerance::exactLength).XY();
      targets.push_back(
        {point, gp_XY(std::clamp(uv.X(), 0.0, lu), std::clamp(uv.Y(), 0.0, lv)), 0, {}, {}, 1});
    }
  }
  else
    for (const gp_Pnt& point : asked.points)
      guideTarget(point, 1);
}

bool QuadFill::onSide(const gp_XY& uv) const
{
  const double slack = asked.limits.g0;
  return std::abs(uv.X()) <= slack || std::abs(uv.X() - lu) <= slack || std::abs(uv.Y()) <= slack
         || std::abs(uv.Y() - lv) <= slack;
}

std::vector<Across> QuadFill::layAcross(ShapeAnalysis_Surface& projection) const
{
  std::vector<Across> across(asked.guides.size());
  for (size_t g = 0; g < asked.guides.size(); ++g)
  {
    const gp_XY a =
                  projection.ValueOfUV(asked.guides[g].at(0), BRepExt::tolerance::exactLength).XY(),
                b =
                  projection.ValueOfUV(asked.guides[g].at(1), BRepExt::tolerance::exactLength).XY();
    if (onSide(a) && onSide(b))
      across[g] = {true, {{0.0, a}, {1.0, b}}};
  }
  crossLaid(across);
  layEnds(across, projection);
  return across;
}

void QuadFill::crossLaid(std::vector<Across>& across) const
{
  for (size_t g = 0; g < asked.guides.size(); ++g)
    for (size_t h = g + 1; h < asked.guides.size(); ++h)
    {
      if (!across[g].laid || !across[h].laid)
        continue;
      double least = HUGE_VAL, fg = 0, fh = 0;
      for (int i = 0; i <= crossingPlaces; ++i)
        for (int j = 0; j <= crossingPlaces; ++j)
        {
          const double x = static_cast<double>(i) / crossingPlaces,
                       y = static_cast<double>(j) / crossingPlaces;
          if (const double d = asked.guides[g].at(x).SquareDistance(asked.guides[h].at(y));
              d < least)
          {
            least = d;
            fg    = x;
            fh    = y;
          }
        }
      // Closer, each way round the nearest places, to well below the tolerance.
      for (double step = 1.0 / crossingPlaces; step > finestStep; step /= 8)
      {
        const double x0 = fg, y0 = fh;
        for (int i = -8; i <= 8; ++i)
          for (int j = -8; j <= 8; ++j)
          {
            const double x = std::clamp(x0 + i * step / 8, 0.0, 1.0),
                         y = std::clamp(y0 + j * step / 8, 0.0, 1.0);
            if (const double d = asked.guides[g].at(x).SquareDistance(asked.guides[h].at(y));
                d < least)
            {
              least = d;
              fg    = x;
              fh    = y;
            }
          }
      }
      // Within the tolerance they cross there, and a conflict at that place is theirs.
      if (std::sqrt(least) > asked.limits.g0)
        continue;
      const gp_XY a = across[g].stops[0].second, b = across[g].stops.back().second,
                  c = across[h].stops[0].second, d = across[h].stops.back().second;
      const double det = (b - a).Crossed(d - c);
      if (std::abs(det)
          <= BRepExt::tolerance::parameterMerge * (b - a).Modulus() * (d - c).Modulus())
        continue;
      const gp_XY cross = a + (b - a) * ((c - a).Crossed(d - c) / det);
      across[g].stops.push_back({fg, cross});
      across[h].stops.push_back({fh, cross});
    }
}

void QuadFill::layEnds(std::vector<Across>& across, ShapeAnalysis_Surface& projection) const
{
  const size_t laidFirst = asked.guides.size();
  for (size_t g = 0; g < laidFirst; ++g)
  {
    if (across[g].laid)
      continue;
    const auto end = [&](double f) -> std::optional<gp_XY> {
      const gp_Pnt at = asked.guides[g].at(f);
      const gp_XY  uv = projection.ValueOfUV(at, BRepExt::tolerance::exactLength).XY();
      if (onSide(uv))
        return uv;
      for (size_t h = 0; h < laidFirst; ++h)
      {
        if (h == g || !across[h].laid || across[h].stops.size() < 2)
          continue;
        double least = HUGE_VAL, foot = 0;
        for (int i = 0; i <= crossingPlaces; ++i)
        {
          const double x = static_cast<double>(i) / crossingPlaces;
          if (const double d = at.SquareDistance(asked.guides[h].at(x)); d < least)
          {
            least = d;
            foot  = x;
          }
        }
        for (double step = 1.0 / crossingPlaces; step > finestStep; step /= 8)
        {
          const double x0 = foot;
          for (int i = -8; i <= 8; ++i)
          {
            const double x = std::clamp(x0 + i * step / 8, 0.0, 1.0);
            if (const double d = at.SquareDistance(asked.guides[h].at(x)); d < least)
            {
              least = d;
              foot  = x;
            }
          }
        }
        if (std::sqrt(least) <= asked.limits.g0)
          return laidAt(across[h].stops, foot);
      }
      return std::nullopt;
    };
    const auto a = end(0), b = end(1);
    if (a && b)
      across[g] = {true, {{0.0, *a}, {1.0, *b}}};
  }
}

Built QuadFill::measured(const occ::handle<Geom_BSplineSurface>& surface) const
{
  Built built{fitting ? 2 : 1,
              {},
              measureAll(surface,
                         asked,
                         [&](const Side& side, double f, const gp_Pnt&) -> std::optional<Place> {
                           return std::any_of(
                                    asked.sides.begin(),
                                    asked.sides.end(),
                                    [&](const Side& boundary) { return &side == &boundary; })
                                    ? std::optional<Place>(Place{place(side, f), false})
                                    : std::nullopt;
                         })};
  for (size_t side = 0; asked.controls.flow != PatchFlow::Natural && side < 4; ++side)
    for (int k = 0; k < places && !gap[side]; ++k)
    {
      const double     f      = (k + .5) / places;
      const gp_XY      uv     = along(side, f);
      const SurfaceJet jet    = surfaceJet(op, surface, uv.X(), uv.Y());
      const gp_Vec     across = side % 2 ? jet.u : jet.v;
      const double     angle =
        std::asin(std::min(1.0, std::abs(across.Normalized().Dot(axis(side, f).Normalized()))))
        * 180 / M_PI;
      if (angle > asked.limits.g1)
      {
        built.measured.within = false;
        built.measured.excess = std::max(built.measured.excess, angle / asked.limits.g1);
        built.measured.worst  = "boundary UV flow is missed by " + number(angle) + " degrees";
      }
    }
  built.face = BRepBuilderAPI_MakeFace(surface, Precision::Confusion()).Face();
  // The orders a rim asks along all of it.
  const auto least = [&](size_t s) {
    int lowest = 2;
    for (const Side* part : runs[s].parts)
      lowest = std::min(lowest, part->order);
    return lowest;
  };
  if (gap[1] && gap[3] && !gap[0] && !gap[2])
    built.crossing = Crossing{true, least(0), least(2)};
  else if (gap[0] && gap[2] && !gap[1] && !gap[3])
    built.crossing = Crossing{false, least(3), least(1)};
  return built;
}

void QuadFill::finer()
{
  for (int i = 0; i < 2; ++i)
  {
    const std::vector<double> knots = flatKnots(curves[i]);
    for (size_t k = 0; k + 1 < knots.size(); ++k)
      if (knots[k + 1] > knots[k])
        for (int c : {i, i + 2})
          curves[c]->InsertKnot((knots[k] + knots[k + 1]) / 2);
  }
}

} // namespace

std::optional<Built> foursided(const Asked& asked, const Frame& frame)
{
  PatchPhase phase("quad");
  const auto found = quadRuns(asked.sides, asked.limits.g1);
  if (!found)
    return {};
  std::vector<Run> runs = *found;
  QuadFill         quad(asked, frame, std::move(runs));
  return quad.run();
}

std::optional<Built> chainedQuad(const Asked& asked, const Frame& frame)
{
  if (asked.sides.size() <= 4 || !asked.guides.empty() || !asked.points.empty()
      || asked.flatness >= 0 || asked.controls.flow != PatchFlow::Natural
      || asked.controls.quad >= PatchQuad::Off || asked.controls.quality != PatchQuality::Auto
      || std::any_of(asked.sides.begin(), asked.sides.end(), [](const Side& side) {
           return side.order != 0;
         }))
    return std::nullopt;
  std::vector<size_t> starts;
  for (size_t i = 0; i < asked.sides.size(); ++i)
  {
    const double turn = asked.sides[i].tangent(1).Normalized().Dot(
      asked.sides[(i + 1) % asked.sides.size()].tangent(0).Normalized());
    if (std::acos(std::clamp(turn, -1.0, 1.0)) * 180 / M_PI <= asked.limits.g1)
      continue;
    if (std::abs(turn) > sharpestCorner)
      return std::nullopt;
    starts.push_back((i + 1) % asked.sides.size());
  }
  if (starts.size() != 4)
    return std::nullopt;
  std::sort(starts.begin(), starts.end());
  occ::handle<Geom_BSplineCurve> curves[4];
  double                         lengths[4] = {};
  for (size_t group = 0; group < 4; ++group)
  {
    const size_t first = starts[group], stop = starts[(group + 1) % 4];
    std::optional<GeomConvert_CompCurveToBSplineCurve> joined;
    for (size_t i = first;; i = (i + 1) % asked.sides.size())
    {
      if (i == stop && i != first)
        break;
      const Side& side  = asked.sides[i];
      auto        piece = plain(side, side.length, false, asked.limits.g0);
      if (piece.IsNull())
        return std::nullopt;
      if (!joined)
        joined.emplace(piece);
      else if (!joined->Add(piece,
                            asked.limits.g0,
                            true,
                            true,
                            std::max(piece->Degree(), joined->BSplineCurve()->Degree())))
        return std::nullopt;
      lengths[group] += side.length;
    }
    curves[group] = joined->BSplineCurve();
    if (group >= 2)
      curves[group]->Reverse();
  }
  const double u = (lengths[0] + lengths[2]) / 2, v = (lengths[1] + lengths[3]) / 2;
  for (int i = 0; i < 4; ++i)
  {
    NCollection_Array1<double> knots = curves[i]->Knots();
    BSplCLib::Reparametrize(0, i % 2 ? v : u, knots);
    curves[i]->SetKnots(knots);
  }
  for (int i = 0; i < 2; ++i)
  {
    auto &    a = curves[i], &b = curves[i + 2];
    const int degree = std::max(a->Degree(), b->Degree());
    a->IncreaseDegree(degree);
    b->IncreaseDegree(degree);
    const double same = BRepExt::tolerance::parameterMerge * (i ? v : u);
    a->InsertKnots(b->Knots(), b->Multiplicities(), same, false);
    b->InsertKnots(a->Knots(), a->Multiplicities(), same, false);
    if (a->NbPoles() > mostCoonsPoles || b->NbPoles() > mostCoonsPoles)
      return std::nullopt;
  }
  const auto surface = surfaceOf(coons(curves[0], curves[1], curves[2], curves[3]));
  const auto layout  = [&](const gp_Pnt& point) {
    double distance = HUGE_VAL;
    gp_XY  uv;
    for (int i = 0; i < 4; ++i)
    {
      GeomAPI_ProjectPointOnCurve projected(point, curves[i]);
      if (projected.NbPoints() == 0 || projected.LowerDistance() >= distance)
        continue;
      distance       = projected.LowerDistance();
      const double t = projected.LowerDistanceParameter();
      uv = i == 0 ? gp_XY(t, 0) : i == 1 ? gp_XY(u, t) : i == 2 ? gp_XY(t, v) : gp_XY(0, t);
    }
    return uv;
  };
  return continued(surface, asked, frame, layout);
}

} // namespace BRepExt::patch_fill
