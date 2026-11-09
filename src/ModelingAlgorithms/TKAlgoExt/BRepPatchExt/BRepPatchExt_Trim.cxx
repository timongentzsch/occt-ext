// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// The rim on the surface: the loop's sides laid into the surface's parameters, the trimmed face
// they bound, and a rim of one spline worked on in a frame of its own.

#include <BRepPatchExt_Fill.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>

#include <BRepToolsExt_GeomUtil.hxx>
#include <GeomAlgoExt_Cancel.hxx>
#include <GeomAlgoExt_PatchRim.hxx>

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepLib.hxx>
#include <BRep_Builder.hxx>
#include <BSplCLib.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomConvert.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_Line.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <TopoDS_Wire.hxx>

#include <limits>

namespace BRepExt::patch_fill
{
namespace
{

// Shared pole access keeps private Bezier and nonperiodic B-spline copies on the same
// frame/budget path. Positive unchanged weights retain the convex-hull displacement bound.
struct RimPoles
{
  occ::handle<Geom_BezierCurve>  bezier;
  occ::handle<Geom_BSplineCurve> spline;

  explicit RimPoles(const occ::handle<Geom_Curve>& curve)
      : bezier(occ::down_cast<Geom_BezierCurve>(curve)),
        spline(occ::down_cast<Geom_BSplineCurve>(curve))
  {
  }

  bool supported() const
  {
    if (bezier.IsNull() && (spline.IsNull() || spline->IsPeriodic()))
      return false;
    for (int i = 1; i <= count(); ++i)
    {
      if ((i & 63) == 0)
        checkCancelled(op);
      const double weight = !bezier.IsNull() ? bezier->Weight(i) : spline->Weight(i);
      if (!(weight > 0) || !std::isfinite(weight))
        return false;
    }
    return true;
  }

  int count() const { return !bezier.IsNull() ? bezier->NbPoles() : spline->NbPoles(); }

  gp_Pnt pole(int i) const { return !bezier.IsNull() ? bezier->Pole(i) : spline->Pole(i); }

  void set(int i, const gp_Pnt& point) const
  {
    if (!bezier.IsNull())
      bezier->SetPole(i, point);
    else
      spline->SetPole(i, point);
  }
};

//! What share of the tolerance the curve of an edge in the surface's parameters may add to
//! the distance of the surface from the edge, which no such curve can do anything about.
constexpr double preimageShare = 0.05;

//! The most poles such a curve is given before it is taken as it is.
constexpr int preimagePoles = 1000;

//! The most passes in which its poles are moved onto their places: each takes a third or more
//! off what is left, so that these leave a millionth of it.
constexpr int preimagePasses = 40;

//! How much more than the deviation that was measured an edge is given as its tolerance: the
//! deviation is read at places, and between them it may be a little more.
constexpr double edgeHeadroom = 1.1;

} // namespace

std::optional<gp_Trsf> localSplineRim(Asked& asked, bool* copyBudgetMiss)
{
  if (copyBudgetMiss)
    *copyBudgetMiss = false;
  if (asked.sides.empty() || !asked.guides.empty() || !asked.points.empty()
      || std::any_of(asked.sides.begin(), asked.sides.end(), [](const Side& side) {
           const RimPoles poles(side.curve);
           if (!side.support.IsNull() || !poles.supported())
             return true;
           // Positive spline bases bound pole displacement on their natural domain,
           // not on an extrapolated edge interval.
           return !poles.spline.IsNull()
                  && (!std::isfinite(side.from) || !std::isfinite(side.to)
                      || !std::isfinite(side.curve->FirstParameter())
                      || !std::isfinite(side.curve->LastParameter())
                      || std::min(side.from, side.to) < side.curve->FirstParameter()
                      || std::max(side.from, side.to) > side.curve->LastParameter());
         }))
    return std::nullopt;
  auto sides = asked.sides;
  for (Side& side : sides)
  {
    auto curve = occ::down_cast<Geom_Curve>(side.curve->Copy());
    if (const auto bezier = occ::down_cast<Geom_BezierCurve>(curve); !bezier.IsNull())
    {
      const double from = std::min(side.from, side.to), to = std::max(side.from, side.to);
      if (from != bezier->FirstParameter() || to != bezier->LastParameter())
        bezier->Segment(from, to);
      if (side.to < side.from)
        bezier->Reverse();
      side.from = bezier->FirstParameter();
      side.to   = bezier->LastParameter();
    }
    // A B-spline keeps its complete basis, weights and active interval. Direction
    // lives in from/to; do not reflect knots or subtract parameter origins.
    side.curve  = curve;
    side.length = GCPnts_AbscissaPoint::Length(GeomAdaptor_Curve(curve),
                                               std::min(side.from, side.to),
                                               std::max(side.from, side.to));
  }
  gp_XYZ centre(0, 0, 0);
  double perimeter = 0;
  size_t first     = 0;
  for (size_t i = 0; i < sides.size(); ++i)
  {
    const Side& side = sides[i];
    if (!(side.length > 0) || !std::isfinite(side.length))
      return std::nullopt;
    centre += side.at(.5).XYZ() * side.length;
    perimeter += side.length;
    if (side.length > sides[first].length * (1 + 1e-12))
      first = i;
  }
  if (!(perimeter > 0) || !std::isfinite(perimeter))
    return std::nullopt;
  // A tied anchor has no unique intrinsic ordering; retain the original path.
  for (size_t i = 0; i < sides.size(); ++i)
    if (i != first && std::abs(sides[i].length - sides[first].length) <= 1e-12 * perimeter)
      return std::nullopt;
  centre /= perimeter;
  if (!finite(centre))
    return std::nullopt;
  const gp_Pnt middle(centre);
  const int    selected = sides[first].input;
  bool         reverse = false, directed = false;
  for (double f : {0.0, .25})
  {
    const double a = sides[first].at(f).SquareDistance(middle);
    const double b = sides[first].at(1 - f).SquareDistance(middle);
    if (std::abs(a - b) > 1e-12 * perimeter * perimeter)
    {
      reverse  = b > a;
      directed = true;
      break;
    }
  }
  if (!directed)
    return std::nullopt;
  if (reverse)
  {
    std::reverse(sides.begin(), sides.end());
    for (Side& side : sides)
    {
      if (!occ::down_cast<Geom_BezierCurve>(side.curve).IsNull())
      {
        const double from = side.curve->ReversedParameter(side.to);
        const double to   = side.curve->ReversedParameter(side.from);
        side.curve->Reverse();
        side.from = from;
        side.to   = to;
      }
      else
      {
        std::swap(side.from, side.to);
      }
      side.originalReversed = !side.originalReversed;
    }
  }
  const auto start = std::find_if(sides.begin(), sides.end(), [&](const Side& side) {
    return side.input == selected;
  });
  std::rotate(sides.begin(), start, sides.end());
  const Frame frame = meanPlane(polygon(sides), perimeter);
  gp_Trsf     world;
  world.SetValues(frame.x.X(),
                  frame.y.X(),
                  frame.z.X(),
                  frame.origin.X(),
                  frame.x.Y(),
                  frame.y.Y(),
                  frame.z.Y(),
                  frame.origin.Y(),
                  frame.x.Z(),
                  frame.y.Z(),
                  frame.z.Z(),
                  frame.origin.Z());
  const double allowance = asked.limits.g0 * 1e-6;
  if (!(allowance > 0) || !std::isfinite(allowance))
    return std::nullopt;
  const double step = std::min(std::ldexp(1.0, std::ilogb(perimeter) - 40), allowance);
  if (!(step > 0) || !std::isfinite(step))
    return std::nullopt;
  double error = 0;
  double magnitude =
    std::max({std::abs(frame.origin.X()), std::abs(frame.origin.Y()), std::abs(frame.origin.Z())});
  for (Side& side : sides)
  {
    const RimPoles curve(side.curve);
    for (int i = 1; i <= curve.count(); ++i)
    {
      if ((i & 63) == 0)
        checkCancelled(op);
      const gp_Vec at(frame.origin, curve.pole(i));
      gp_Pnt       point(at.Dot(frame.x), at.Dot(frame.y), at.Dot(frame.z));
      for (int axis = 1; axis <= 3; ++axis)
        magnitude =
          std::max({magnitude, std::abs(curve.pole(i).Coord(axis)), std::abs(point.Coord(axis))});
      for (int axis = 1; axis <= 3; ++axis)
        point.SetCoord(axis, std::round(point.Coord(axis) / step) * step);
      const double displacement = curve.pole(i).Distance(point.Transformed(world));
      if (!std::isfinite(displacement))
        return std::nullopt;
      error = std::max(error, displacement);
      curve.set(i, point);
    }
    side.edge = BRepBuilderAPI_MakeEdge(side.curve,
                                        std::min(side.from, side.to),
                                        std::max(side.from, side.to))
                  .Edge();
    if (side.to < side.from)
      side.edge.Reverse();
    side.length = GCPnts_AbscissaPoint::Length(GeomAdaptor_Curve(side.curve),
                                               std::min(side.from, side.to),
                                               std::max(side.from, side.to));
  }
  // Positive unchanged spline weights make the maximum pole displacement a bound on the
  // entire curve. If world-coordinate resolution cannot hold this private-copy
  // budget, leave the original construction alone rather than perturbing its rim.
  // Round-tripping at large coordinates can hide displacement below one world ulp.
  // Reserve the dot products/transform rounding before accepting the private copy.
  const double roundoff = 32 * std::numeric_limits<double>::epsilon() * magnitude;
  if (error + roundoff > allowance)
  {
    if (copyBudgetMiss)
      *copyBudgetMiss = true;
    return std::nullopt;
  }
  asked.sides = std::move(sides);
  return world;
}

std::optional<gp_Trsf> rebasedSplineRim(Asked& asked, gp_Trsf& anchor)
{
  Asked          local = asked;
  const RimPoles first(local.sides.front().curve);
  if (!first.supported())
    return std::nullopt;
  const gp_Pnt origin = first.pole(1);
  for (Side& side : local.sides)
  {
    if (!RimPoles(side.curve).supported())
      return std::nullopt;
    const auto     copy = occ::down_cast<Geom_Curve>(side.curve->Copy());
    const RimPoles curve(copy);
    for (int i = 1; i <= curve.count(); ++i)
    {
      if ((i & 63) == 0)
        checkCancelled(op);
      gp_Pnt point = curve.pole(i);
      for (int axis = 1; axis <= 3; ++axis)
      {
        const double a = point.Coord(axis), b = origin.Coord(axis), difference = a - b;
        const double virtualB = a - difference, virtualA = difference + virtualB;
        const double residual = (a - virtualA) + (virtualB - b);
        if (!std::isfinite(difference) || !std::isfinite(residual) || residual != 0)
          return std::nullopt;
        point.SetCoord(axis, difference);
      }
      curve.set(i, point);
    }
    side.curve = copy;
  }
  const auto placement = localSplineRim(local);
  if (!placement)
    return std::nullopt;
  anchor.SetTranslation(gp_Vec(origin.XYZ()));
  asked = std::move(local);
  return placement;
}

Rim preimage(const occ::handle<Geom_BSplineSurface>&    surface,
             const Side&                                side,
             const Frame&                               frame,
             double                                     tolerance,
             double                                     size,
             const std::function<gp_XY(const gp_Pnt&)>& layout,
             const std::function<gp_XY(double)>&        parameterLayout)
{
  const auto   project = [&](const gp_Pnt& p) { return layout ? layout(p) : frame.project(p); };
  const double first = std::min(side.from, side.to), last = std::max(side.from, side.to);
  occ::handle<Geom2d_BSplineCurve> flat;
  // A line, a Bezier curve and a B-spline keep their parameters as B-splines, and laid into
  // a plane they are their own poles laid into it.
  occ::handle<Geom_Curve> basis = side.curve;
  while (const auto cut = occ::down_cast<Geom_TrimmedCurve>(basis))
    basis = cut->BasisCurve();
  if (basis->IsKind(STANDARD_TYPE(Geom_Line)) || basis->IsKind(STANDARD_TYPE(Geom_BSplineCurve))
      || basis->IsKind(STANDARD_TYPE(Geom_BezierCurve)))
  {
    const auto curve =
      GeomConvert::CurveToBSplineCurve(new Geom_TrimmedCurve(side.curve, first, last));
    if (curve->IsPeriodic())
      curve->SetNotPeriodic();
    bool same = std::abs(curve->FirstParameter() - first) + std::abs(curve->LastParameter() - last)
                < BRepExt::tolerance::parameter * std::max(1.0, std::abs(first) + std::abs(last));
    for (int k = 0; same && k <= 4; ++k)
    {
      const double t = first + (last - first) * k / 4;
      same =
        curve->Value(t).Distance(side.curve->Value(t)) < BRepExt::tolerance::exactLength * size;
    }
    // One with weights goes the way of an arc: its poles are not found as a plain
    // curve's are.
    if (same && !curve->IsRational() && !parameterLayout)
    {
      NCollection_Array1<gp_Pnt2d> poles(1, curve->NbPoles());
      for (int k = 1; k <= curve->NbPoles(); ++k)
        poles.SetValue(k, gp_Pnt2d(project(curve->Pole(k))));
      flat =
        new Geom2d_BSplineCurve(poles, curve->Knots(), curve->Multiplicities(), curve->Degree());
    }
  }
  if (flat.IsNull())
  {
    // An arc and the like run by an angle, which no B-spline does: a curve of degree 5
    // near it, to begin with, its poles at the places of the curve that they stand for.
    // Of degree 5, since twice the knots then bring it 64 times nearer: a cubic of the
    // thousand poles allowed is 2e-9 of an arc's size off, which is more than a
    // millionth on a model a thousand across.
    const int                    count = 16, degree = 5;
    const std::vector<double>    knots = clampedKnots(first, last, count, degree);
    NCollection_Array1<gp_Pnt2d> poles(1, count + degree);
    for (int k = 0; k < count + degree; ++k)
    {
      double sum = 0;
      for (int i = 1; i <= degree; ++i)
        sum += knots[k + i];
      poles.SetValue(k + 1,
                     gp_Pnt2d(parameterLayout ? parameterLayout(sum / degree)
                                              : project(side.curve->Value(sum / degree))));
    }
    NCollection_Array1<double> distinct(1, count + 1);
    NCollection_Array1<int>    times(1, count + 1);
    for (int k = 0; k <= count; ++k)
    {
      distinct.SetValue(k + 1, first + (last - first) * k / count);
      times.SetValue(k + 1, k == 0 || k == count ? degree + 1 : 1);
    }
    flat = new Geom2d_BSplineCurve(poles, distinct, times, degree);
  }
  // A same-parameter support curve can be only C1 at its own knots. An everywhere-C4
  // preimage cannot follow those small changes of parameter at an absolute tolerance.
  occ::handle<Geom2d_Curve> onHeld = side.onHeld;
  while (const auto cut = occ::down_cast<Geom2d_TrimmedCurve>(onHeld))
    onHeld = cut->BasisCurve();
  if (const auto spline = occ::down_cast<Geom2d_BSplineCurve>(onHeld); !spline.IsNull())
  {
    flat->IncreaseDegree(std::max(flat->Degree(), spline->Degree()));
    for (int k = 1; k <= spline->NbKnots(); ++k)
    {
      const double t = spline->Knot(k);
      if (t <= first || t >= last)
        continue;
      NCollection_Array1<double> knots(1, 1);
      knots.SetValue(1, t);
      NCollection_Array1<int> times(1, 1);
      times.SetValue(1, flat->Degree() - spline->Degree() + spline->Multiplicity(k));
      flat->InsertKnots(knots, times, 1e-10 * std::max(1.0, last - first), false);
    }
  }
  const auto off = [&](double t, const gp_Pnt2d& uv) {
    return surface->Value(uv.X(), uv.Y()).Distance(side.curve->Value(t));
  };
  Rim    rim{flat, 0};
  double before = HUGE_VAL;
  for (;;)
  {
    const int                        count = flat->NbPoles(), degree = flat->Degree();
    const NCollection_Array1<double> sequence = flat->KnotSequence();
    // Each pole's own place on the curve: the mean of its knots.
    std::vector<double> at(count);
    std::vector<gp_XY>  wanted(count), miss(count);
    for (int i = 0; i < count; ++i)
    {
      double sum = 0;
      for (int k = 1; k <= degree; ++k)
        sum += sequence(sequence.Lower() + i + k);
      at[i]     = std::clamp(sum / degree, first, last);
      wanted[i] = nearest(surface,
                          parameterLayout ? parameterLayout(at[i]) : flat->Value(at[i]).XY(),
                          side.meetAt(at[i]));
    }
    // The curve through those places, which are as many as it has poles and lie one under
    // each: its poles solve a banded system. Where that has no solution (knots as often
    // as the order, which no curve here has) each pole is moved by what its own place is
    // off instead, over and over.
    NCollection_Array1<gp_Pnt2d> through(1, count);
    NCollection_Array1<double>   parameters(1, count);
    NCollection_Array1<int>      contacts(1, count);
    for (int i = 0; i < count; ++i)
    {
      through.SetValue(i + 1, gp_Pnt2d(wanted[i]));
      parameters.SetValue(i + 1, at[i]);
      contacts.SetValue(i + 1, 0);
    }
    // Not where the curve is at its places as it stands (over a net that is a height
    // above the mean plane the side laid into that plane is there already), and not kept
    // where the curve it gives is further from them than the one before.
    const auto furthest = [&] {
      double worst = 0;
      for (int i = 0; i < count; ++i)
        worst =
          std::max(worst,
                   surface->Value(wanted[i].X(), wanted[i].Y())
                     .Distance(surface->Value(flat->Value(at[i]).X(), flat->Value(at[i]).Y())));
      return worst;
    };
    const double asGiven  = furthest();
    int          unsolved = asGiven < preimageShare * tolerance / 10;
    if (!unsolved)
    {
      const NCollection_Array1<gp_Pnt2d> given = flat->Poles();
      BSplCLib::Interpolate(degree, sequence, parameters, contacts, through, unsolved);
      for (int i = 0; !unsolved && i < count; ++i)
        flat->SetPole(i + 1, through(i + 1));
      if (!unsolved && furthest() > asGiven)
      {
        unsolved = 1;
        for (int i = 0; i < count; ++i)
          flat->SetPole(i + 1, given(i + 1));
      }
    }
    for (int pass = 0; unsolved && pass < preimagePasses; ++pass)
    {
      double worst = 0;
      for (int i = 0; i < count; ++i)
      {
        miss[i] = wanted[i] - flat->Value(at[i]).XY();
        worst =
          std::max(worst,
                   surface->Value(wanted[i].X(), wanted[i].Y())
                     .Distance(surface->Value(flat->Value(at[i]).X(), flat->Value(at[i]).Y())));
      }
      if (worst < preimageShare * tolerance / 10)
        break;
      for (int i = 0; i < count; ++i)
        flat->SetPole(i + 1, gp_Pnt2d(flat->Pole(i + 1).XY() + miss[i]));
    }
    // Between those places: how far the curve is from where it should be, in space.
    double beside = 0;
    rim.deviation = std::max(off(first, flat->Value(first)), off(last, flat->Value(last)));
    std::vector<double> middles;
    for (int k = sequence.Lower() + degree; k < sequence.Upper() - degree; ++k)
    {
      const double from = sequence(k), to = sequence(k + 1);
      if (to <= from)
        continue;
      middles.push_back((from + to) / 2);
      for (int part = 1; part <= 31; ++part)
      {
        const double   t  = from + (to - from) * part / 32;
        const gp_Pnt2d uv = flat->Value(t);
        const gp_XY    best =
          nearest(surface, parameterLayout ? parameterLayout(t) : uv.XY(), side.meetAt(t));
        beside =
          std::max(beside,
                   surface->Value(uv.X(), uv.Y()).Distance(surface->Value(best.X(), best.Y())));
        rim.deviation = std::max(rim.deviation, off(t, uv));
      }
    }
    // Twice the knots bring a curve that is being followed at all sixteen times nearer;
    // one that comes no nearer by half is not following (the surface has no one nearest
    // point there), and more knots would not make it.

    if (beside <= preimageShare * tolerance || beside > before / 2
        || count + static_cast<int>(middles.size()) > preimagePoles)
      break;
    before = beside;
    for (double middle : middles)
      flat->InsertKnot(middle);
  }
  return rim;
}

TopoDS_Face trimmed(const occ::handle<Geom_BSplineSurface>& surface,
                    const Asked&                            asked,
                    const std::vector<Rim>&                 rims,
                    bool                                    clockwise)
{
  BRep_Builder             builder;
  const size_t             count = asked.sides.size();
  std::vector<TopoDS_Edge> edges(count);
  std::vector<double>      gaps(count);
  for (size_t i = 0; i < count; ++i)
  {
    const Side& side = asked.sides[i];
    // What the edge's curve is off the surface where the face holds it, with room for the
    // places between those it was read at: an edge that is off its own face by more than
    // the tolerance is that far off the patch too, which meets the face. OpenCascade reads
    // no edge more tightly than its own least distance.
    gaps[i] = std::max(Precision::Confusion(), rims[i].deviation * edgeHeadroom);
    if (side.free)
    {
      // A gap is what the surface makes of its line in the parameters.
      builder.MakeEdge(edges[i]);
      builder.UpdateEdge(edges[i], rims[i].onSurface, surface, TopLoc_Location(), gaps[i]);
      builder.Range(edges[i], side.from, side.to);
      if (!BRepLib::BuildCurve3d(edges[i], asked.limits.g0))
        fail(op, no_solution, "a gap of the boundary has no curve on the surface");
      continue;
    }
    builder.MakeEdge(edges[i], side.curve, gaps[i]);
    builder.Range(edges[i], std::min(side.from, side.to), std::max(side.from, side.to));
    builder.UpdateEdge(edges[i], rims[i].onSurface, surface, TopLoc_Location(), gaps[i]);
  }
  for (size_t i = 0; i < count; ++i)
  {
    const size_t  before = (i + count - 1) % count;
    const Side &  side = asked.sides[i], &last = asked.sides[before];
    TopoDS_Vertex vertex;
    builder.MakeVertex(vertex,
                       side.at(0),
                       side.at(0).Distance(last.at(1))
                         + std::max({gaps[i],
                                     gaps[before],
                                     rims[i].deviation * edgeHeadroom,
                                     rims[before].deviation * edgeHeadroom}));
    // A side that runs against its curve starts at the curve's end.
    builder.Add(edges[i], vertex.Oriented(side.from < side.to ? TopAbs_FORWARD : TopAbs_REVERSED));
    builder.Add(edges[before],
                vertex.Oriented(last.from < last.to ? TopAbs_REVERSED : TopAbs_FORWARD));
  }
  TopoDS_Wire wire;
  builder.MakeWire(wire);
  for (size_t i = 0; i < count; ++i)
    builder.Add(wire,
                edges[i].Oriented((asked.sides[i].from < asked.sides[i].to) != clockwise
                                    ? TopAbs_FORWARD
                                    : TopAbs_REVERSED));
  wire.Closed(true);
  TopoDS_Face face;
  builder.MakeFace(face, surface, Precision::Confusion());
  builder.Add(face, wire);
  return face;
}

} // namespace BRepExt::patch_fill
