// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepPatchExt_Fill.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>

#include <BRepAlgoExt_PatchLoop.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <GeomAlgoExt_Cancel.hxx>
#include <GeomAlgoExt_PatchRim.hxx>
#include <GeomAlgoExt_Phase.hxx>

#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepFill_Filling.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>

// Patch: one face inside a closed loop of edges. Each edge may ask the face to go over smoothly
// into a face it belongs to (G1: one tangent plane; G2: one bend as well), and curves and points
// inside the loop may be given for the face to pass through.
//
// Constructions, the first that applies and is measured to meet what was asked:
//   - a plane, where the loop and everything asked of it are flat;
//   - an exact rational continuation of a torus, a ball or a plane/round/plane strip;
//   - four sides: an untrimmed surface with the sides as its own edges, the Coons patch of the
//     four curves, its inner poles then moved to meet G1 and G2;
//   - any loop: a surface over the loop's mean plane, trimmed by the loop. Its net is finer
//     where the rim turns fast, and may change its bend across lines through the corners where
//     the faces along the rim do; each edge of the face lies on the surface where the surface is
//     nearest to it, parameter by parameter.
// The last two are fitted by patch_fit. `method` Plate is OpenCascade's plate instead, kept to
// measure against.
//
// Whatever was built is measured again from outside, at the places Measure Continuity reads and
// for a trimmed face in the way it reads them, and looked at inside its loop (patch_inside);
// no surface that misses its constraints or is not one to keep is returned. Sewn readings
// include the input edges' existing tangential parameter uncertainty, and are reported as
// delivered.

namespace BRepExt
{
namespace
{

using namespace patch_fill;

Side sideOf(const TopoDS_Shape& shape, const TopoDS_Shape& face, int order, int input)
{
  if (shape.ShapeType() != TopAbs_EDGE || BRep_Tool::Degenerated(TopoDS::Edge(shape)))
    fail(op, invalid_input, "the boundary and the guides must be edges");
  if (order < 0 || order > 2)
    fail(op, invalid_input, "continuity is 0 (G0), 1 (G1) or 2 (G2)");
  Side side{input, TopoDS::Edge(shape), {}, !face.IsNull() ? order : 0, {}, 0, 0, 0, 0, {}, {}};
  if (order > 0 && face.IsNull())
    fail(op, invalid_input, "G1 and G2 need the face that the edge is to be smooth with");
  if (!face.IsNull())
  {
    if (face.ShapeType() != TopAbs_FACE)
      fail(op, invalid_input, "a support must be a face");
    // The edge as its face holds it: a seam has two parameter curves, told apart by that.
    bool found = false;
    for (TopExp_Explorer it(face, TopAbs_EDGE); it.More() && !found; it.Next())
      if (it.Current().IsSame(shape))
      {
        side.edge = TopoDS::Edge(it.Current());
        found     = true;
      }
    // An edge laid on the face's surface is held to it too (a patch for another operation,
    // whose new edges are not yet the face's).
    if (!found && side.order > 0)
    {
      double first, last;
      found = !BRep_Tool::CurveOnSurface(side.edge, TopoDS::Face(face), first, last).IsNull();
    }
    if (!found && side.order > 0)
      fail(op, invalid_input, "an edge must belong to the face it is to be smooth with");
    if (found)
    {
      side.support = TopoDS::Face(face);
      double first, last;
      side.onHeld = BRep_Tool::CurveOnSurface(side.edge, TopoDS::Face(face), first, last);
      if (!side.onHeld.IsNull())
        side.held = BRep_Tool::Surface(TopoDS::Face(face));
    }
  }
  side.curve = BRep_Tool::Curve(side.edge, side.from, side.to);
  if (side.curve.IsNull())
    fail(op, invalid_input, "an edge has no curve in space");
  side.length = GCPnts_AbscissaPoint::Length(GeomAdaptor_Curve(side.curve), side.from, side.to);
  for (TopExp_Explorer it(side.edge, TopAbs_VERTEX); it.More(); it.Next())
    side.slack = std::max(side.slack, BRep_Tool::Tolerance(TopoDS::Vertex(it.Current())));
  return side;
}

Side sideOf(const OcctKernel& kernel, uint32_t id, uint32_t supportId, int order, int input)
{
  return sideOf(get(kernel, id), supportId ? get(kernel, supportId) : TopoDS_Shape(), order, input);
}

//! Put the sides in the order and the directions in which they run round the loop. With `open`
//! a boundary that does not close says so and is not refused.
bool arrange(std::vector<Side>& sides, double tolerance, bool open = false)
{
  auto refuse = [](double gap, const gp_Pnt& at) {
    fail(op,
         invalid_input,
         "the boundary does not close: its edges are " + std::to_string(gap) + " apart at ("
           + std::to_string(at.X()) + ", " + std::to_string(at.Y()) + ", " + std::to_string(at.Z())
           + ")");
  };
  for (size_t done = 1; done < sides.size(); ++done)
  {
    const Side&  last    = sides[done - 1];
    const gp_Pnt end     = last.at(1);
    size_t       next    = done;
    double       nearest = HUGE_VAL;
    bool         turned  = false;
    for (size_t i = done; i < sides.size(); ++i)
      for (bool far : {false, true})
      {
        const double distance = end.Distance(sides[i].at(far ? 1 : 0));
        if (distance < nearest)
        {
          nearest = distance;
          next    = i;
          turned  = far;
        }
      }
    if (nearest > std::max({tolerance, last.slack, sides[next].slack}))
    {
      if (open)
        return false;
      refuse(nearest, end);
    }
    std::swap(sides[done], sides[next]);
    if (turned)
      std::swap(sides[done].from, sides[done].to);
  }
  const double gap = sides.back().at(1).Distance(sides.front().at(0));
  if (gap > std::max({tolerance, sides.back().slack, sides.front().slack}))
  {
    if (open)
      return false;
    refuse(gap, sides.back().at(1));
  }
  return true;
}

Built plate(const Asked& asked)
{
  int order = 0;
  for (const Side& side : asked.sides)
    order = std::max(order, side.order);
  // As tight as it is asked anywhere in this kernel, whatever the limits here: it is measured
  // like the others afterwards.
  BRepFill_Filling fill(order + 2 > 3 ? order + 2 : 3,
                        30,
                        3,
                        false,
                        1e-7 * asked.size,
                        1e-7 * asked.size,
                        1e-5,
                        1e-4,
                        8,
                        64);
  for (const Side& side : asked.sides)
    // The plate takes the order as a number, 0 to 2, not as the enumeration's G2.
    if (side.order > 0)
      fill.Add(side.edge, side.support, static_cast<GeomAbs_Shape>(side.order));
    else
      fill.Add(side.edge, GeomAbs_C0);
  for (const Side& guide : asked.guides)
    fill.Add(guide.edge, GeomAbs_C0, false);
  for (const gp_Pnt& point : asked.points)
    fill.Add(point);
  fill.Build();
  if (!fill.IsDone())
    fail(op, no_solution, "the plate found no surface");
  return {4, fill.Face(), measureAll(BRep_Tool::Surface(fill.Face()), asked)};
}

//! The face of the loop `asked` holds once its sides are in turn, held to every check Patch
//! makes; fails with what it misses otherwise.
Built fillFace(Asked asked, PatchMethod method, double flatness, const PatchControls& controls)
{
  const bool gapped    = hasGaps(asked);
  asked.controls       = controls;
  Asked input          = asked;
  bool  copyBudgetMiss = false;
  auto  placement =
    method != PatchMethod::Plate ? localSplineRim(asked, &copyBudgetMiss) : std::nullopt;
  gp_Trsf anchor;
  bool    anchored = false;
  if (!placement && copyBudgetMiss)
  {
    placement = rebasedSplineRim(asked, anchor);
    anchored  = placement.has_value();
  }
  bool canonical = placement.has_value();
  // Free spline interpolation can differ by a few ulps across machines. Quantise a
  // private copy before any net/refinement decisions, not just the linear solve. The
  // original curves are still measured below; the lattice is at most a millionth
  // of the requested positional tolerance and uses a binary scale independent of ulps.
  for (auto* sides : {&asked.sides, &asked.guides})
    for (Side& side : *sides)
    {
      if (!side.held.IsNull())
        continue;
      auto curve = occ::down_cast<Geom_BSplineCurve>(side.curve);
      if (curve.IsNull() || !curve->IsPeriodic()
          || side.at(0).Distance(side.at(1)) > asked.limits.g0)
        continue;
      const double parameterStep =
        std::ldexp(1.0, std::ilogb(std::max(1.0, std::abs(curve->LastParameter()))) - 38);
      bool separated = true;
      for (int i = 2; i <= curve->NbKnots(); ++i)
        separated = separated && curve->Knot(i) - curve->Knot(i - 1) > 16 * parameterStep;
      if (!separated)
        continue; // Never merge a narrow knot interval by rounding.
      curve        = occ::down_cast<Geom_BSplineCurve>(curve->Copy());
      double scale = 1;
      for (int i = 1; i <= curve->NbPoles(); ++i)
        for (int axis = 1; axis <= 3; ++axis)
          scale = std::max(scale, std::abs(curve->Pole(i).Coord(axis)));
      const double step = std::min(std::ldexp(1.0, std::ilogb(scale) - 38), asked.limits.g0 * 1e-6);
      for (int i = 1; i <= curve->NbPoles(); ++i)
      {
        gp_Pnt point = curve->Pole(i);
        for (int axis = 1; axis <= 3; ++axis)
          point.SetCoord(axis, std::round(point.Coord(axis) / step) * step);
        curve->SetPole(i, point);
      }
      for (int i = 1; i <= curve->NbKnots(); ++i)
        curve->SetKnot(i, std::round(curve->Knot(i) / parameterStep) * parameterStep);
      side.from   = std::round(side.from / parameterStep) * parameterStep;
      side.to     = std::round(side.to / parameterStep) * parameterStep;
      side.curve  = curve;
      side.length = GCPnts_AbscissaPoint::Length(GeomAdaptor_Curve(curve),
                                                 std::min(side.from, side.to),
                                                 std::max(side.from, side.to));
      canonical   = true;
    }
  double perimeter = 0;
  for (const Side& side : asked.sides)
    perimeter += side.length;
  const std::vector<gp_Pnt> outline = polygon(asked.sides);
  const Frame               frame   = meanPlane(outline, perimeter);
  for (const gp_Pnt& point : outline)
    asked.size = std::max(asked.size, frame.origin.Distance(point) * 2);
  const std::string eased =
    agreeAtCorners(asked.sides, asked.limits, asked.size, asked.size / ladder[0]);
  // The delivered face is held to what the sides ask near their corners as well.
  for (size_t i = 0; i < input.sides.size(); ++i)
    std::copy_n(asked.sides[i].eases, 2, input.sides[i].eases);

  // A plane in the mean-plane coordinates has orthogonal UV everywhere, so all
  // nonquad flow choices described as Normal can be satisfied exactly there.
  const auto plane = method == PatchMethod::Auto && controls.quad != PatchQuad::On
                       ? planar(asked, frame)
                       : std::nullopt;
  // With Quad sided off, Next and Adjacent both ask for Normal flow: the Natural
  // fill, then the same over its conformal chart (see adjacent). With it on, they are
  // directions across the quad's sides (see foursided).
  const bool conformal = controls.flow >= (gapped ? PatchFlow::Normal : PatchFlow::Next)
                         && controls.quad == PatchQuad::Off && !plane;
  if (controls.flow >= PatchFlow::Next && controls.quad != PatchQuad::On && !plane && !conformal)
    fail(op,
         unsupported_geometry,
         "Next and Adjacent UV flow are not implemented for this boundary; use Natural flow");
  std::optional<Built> built;
  if (plane)
  {
    built = Built{0, *plane, measureAll(BRep_Tool::Surface(*plane), asked)};
  }
  else if (method == PatchMethod::Plate)
  {
    built = plate(asked);
  }
  else if (controls.quad == PatchQuad::On)
  {
    if (!quadRuns(asked.sides, asked.limits.g1))
      fail(op,
           unsupported_geometry,
           "Quad sided requires four sides with regular corners: four edges, or edges that meet "
           "without a corner");
    built = foursided(asked, frame);
    if (!built || !built->measured.within)
      fail(op,
           no_solution,
           built ? built->measured.worst : "the four sides cannot form an untrimmed quad");
  }
  else if (controls.flow > PatchFlow::Natural && !conformal)
  {
    fail(op,
         unsupported_geometry,
         "UV flow controls currently require four regular sides and Quad sided");
  }
  else if (flatness >= 0 || controls.quad == PatchQuad::Off
           || controls.quality > PatchQuality::Auto)
  {
    built = overPlane(asked, frame);
  }
  else if (const auto face = method == PatchMethod::Auto ? planar(asked, frame) : std::nullopt)
  {
    built = Built{0, *face, measureAll(BRep_Tool::Surface(*face), asked)};
  }
  else
  {
    built = splineContinuation(asked, frame);
    if (!built)
      built = toroidal(asked, frame);
    if (!built)
      built = spherical(asked, frame);
    if (!built)
      built = ruledRound(asked, frame);
    std::optional<Built> rotatingQuad;
    if (!built && asked.guides.empty() && asked.points.empty()
        && squarish(asked.sides, asked.limits.g1))
    {
      built = foursided(asked, frame);
      // Retain the established graph preference where a quad turns rapidly.
      // This is candidate ordering, not a fold proof. A regular rotating quad
      // remains available when the preferred trimmed construction cannot deliver.
      if (built && built->measured.within && insideFault(built->face, 0, {}, false))
      {
        if (!insideFault(built->face, 0))
          rotatingQuad = built;
        built.reset();
      }
    }
    if (!built || !built->measured.within)
    {
      try
      {
        built = overPlane(asked, frame);
      }
      catch (const Failure& failure)
      {
        if (failure.reason == "cancelled" || !rotatingQuad)
          throw;
        built = rotatingQuad;
      }
      if (!built->measured.within && rotatingQuad)
        built = rotatingQuad;
    }
  }
  if (!built->measured.within && method != PatchMethod::Plate && !gapped)
    if (const auto quad = chainedQuad(asked, frame))
      built = quad;
  if (conformal && built->measured.within && !built->measured.unkept)
    built = adjacent(asked, frame, *built);
  if (placement && built->measured.within && !built->measured.unkept)
  {
    built->face = TopoDS::Face(BRepBuilderAPI_Transform(built->face, *placement, true).Shape());
    if (anchored)
      built->face.Move(TopLoc_Location(anchor));
    const double size = asked.size;
    asked             = input;
    asked.size        = size;
  }
  if (canonical && built->measured.within && !built->measured.unkept)
  {
    Asked delivered    = input;
    delivered.size     = asked.size;
    const auto surface = BRep_Tool::Surface(built->face);
    built->measured    = measureAll(surface, delivered);
    if (!built->measured.within && !built->measured.unkept)
    {
      // The generic projector's absolute resolution can exceed the requested
      // physical tolerance. Refine its UV seed without moving the already-fitted
      // surface. Successful original measurements and all construction stay intact.
      try
      {
        ShapeAnalysis_Surface project(surface);
        bool                  finite = true;
        const Measured        refined =
          measureAll(surface,
                     delivered,
                     [&](const Side&, double, const gp_Pnt& point) -> std::optional<Place> {
                       const gp_XY seed =
                         project.ValueOfUV(point, BRepExt::tolerance::exactLength).XY();
                       const gp_XY uv = nearest(surface, seed, point);
                       if (!std::isfinite(uv.X()) || !std::isfinite(uv.Y())
                           || !std::isfinite(surface->Value(uv.X(), uv.Y()).Distance(point)))
                       {
                         finite = false;
                         return std::nullopt;
                       }
                       return Place{uv, false};
                     });
        if (finite && refined.within && !refined.unkept && refined.excess < built->measured.excess
            && std::all_of(refined.values.begin(), refined.values.end(), [](double value) {
                 return std::isfinite(value);
               }))
          built->measured = refined;
      }
      catch (const Standard_Failure&)
      {
        // Retain the original typed refusal when projection is not regular.
      }
    }
  }
  measureSpatialBoundary(*built, input, true);
  // Where a corner was eased and nothing met, the corner is the reason.
  if (!eased.empty() && (built->measured.unkept || !built->measured.within))
    fail(op, unsupported_geometry, eased);
  if (built->measured.unkept)
    fail(op,
         no_solution,
         "the surface that meets the edges is not one to keep: " + built->measured.worst);
  if (!built->measured.within)
  {
    std::string reason = "no Patch candidate meets what was asked: " + built->measured.worst;
    if (const auto hint = guideEndHint(input, asked.size))
      reason += "; " + *hint;
    fail(op, no_solution, reason);
  }
  if (!BRepCheck_Analyzer(built->face).IsValid())
    fail(op, no_solution, "the patch is not a valid face");
  // Within every limit along its edges is not yet a patch: see BRepPatchExt_Inside.hxx. Points
  // inside may ask for any bend, so with those only a fold is a fault. A plane, or the one surface
  // that four curves have between them with nothing fitted, bends as its edges make it: there only
  // a fold is a fault.
  const bool fitted = built->construction >= 2;
  if (const auto fault = insideFault(built->face,
                                     fitted && asked.points.empty() ? rimBend(asked) : 0,
                                     built->crossing))
  {
    if (!eased.empty())
      fail(op, unsupported_geometry, eased);
    fail(op, no_solution, "the surface that meets the edges is not one to keep: " + *fault);
  }

  if (method != PatchMethod::Plate)
    verifySewn(*built, asked);
  if (!controls.trim)
  {
    const auto surface = BRep_Tool::Surface(built->face);
    if (occ::down_cast<Geom_BSplineSurface>(surface).IsNull())
      fail(op, unsupported_geometry, "Untrimmed output requires a bounded spline surface");
    built->face = BRepBuilderAPI_MakeFace(surface, Precision::Confusion()).Face();
    if (const auto fault = insideFault(built->face,
                                       fitted && asked.points.empty() ? rimBend(asked) : 0,
                                       built->crossing))
      fail(op, no_solution, "the untrimmed rectangle is not one to keep: " + *fault);
  }
  if (!BRepCheck_Analyzer(built->face).IsValid())
    fail(op, no_solution, "the untrimmed result is not a valid face");
  measureTrimmedInterior(*built, input);
  return *built;
}

//! The face of the loop `asked` holds, held to every check Patch makes; fails with what it
//! misses otherwise.
Built patchFace(Asked asked, PatchMethod method, double flatness, PatchControls controls)
{
  if (!arrange(asked.sides, asked.limits.g0, controls.open && method != PatchMethod::Plate))
    chain(asked.sides, asked.limits.g0);
  const bool gapped = hasGaps(asked);
  // What `open` changes is for a boundary that does not close. One that does is the fill it
  // always was, whoever asks for it.
  controls.open = gapped;
  // An open boundary asked to be a quad, whose runs and gaps are four sides with corners, is the
  // untrimmed quad with the gaps' poles free (see foursided). Where that is not a face to keep,
  // it is filled as any open boundary is: every flow but Natural is then the refit over the
  // conformal chart.
  if (gapped && controls.quad == PatchQuad::On && quadRuns(asked.sides, asked.limits.g1))
    try
    {
      return fillFace(asked, method, flatness, controls);
    }
    catch (const Failure& failure)
    {
      if (failure.reason == "cancelled")
        throw;
    }
  if (gapped)
    controls.quad = PatchQuad::Off;
  return fillFace(asked, method, flatness, controls);
}

PatchResult patchAsAsked(OcctKernel&                  kernel,
                         const std::vector<uint32_t>& edges,
                         const std::vector<uint32_t>& supports,
                         const std::vector<int>&      orders,
                         const std::vector<uint32_t>& guides,
                         const std::vector<double>&   points,
                         double                       tolerance,
                         double                       g1Degrees,
                         double                       g2Share,
                         PatchMethod                  method,
                         double                       flatness,
                         const PatchControls&         controls)
{
  PatchPhase phase("total");
  if (edges.empty() || edges.size() > 64 || supports.size() != edges.size()
      || orders.size() != edges.size() || points.size() % 3 != 0 || method < PatchMethod::Auto
      || method > PatchMethod::Plate)
    fail(op,
         invalid_input,
         "one to 64 boundary edges, each with a support (or 0) and an order, are required");
  if (!(tolerance > 0) || !(g1Degrees > 0) || !(g2Share > 0)
      || !std::isfinite(tolerance + g1Degrees + g2Share))
    fail(op, invalid_input, "the tolerances must be positive and finite");
  if (!std::isfinite(flatness) || (flatness != -1 && (flatness < 0 || flatness > 1))
      || (flatness >= 0 && method == PatchMethod::Plate))
    fail(op, invalid_input, "Flatness must be in 0..1 and uses the variational solver");
  if (controls.quad < PatchQuad::Auto || controls.quad > PatchQuad::On
      || controls.flow < PatchFlow::Natural || controls.flow > PatchFlow::Adjacent
      || controls.quality < PatchQuality::Auto || controls.quality > PatchQuality::Max)
    fail(op, invalid_input, "invalid parameterization, flow or quality control");
  if (method == PatchMethod::Plate
      && (controls.quad != PatchQuad::Auto || controls.flow != PatchFlow::Natural
          || controls.quality != PatchQuality::Auto || !controls.trim))
    fail(op, invalid_input, "parameterization, flow, quality and untrim controls use this solver");
  Asked asked{{}, {}, {}, {tolerance, g1Degrees, g2Share}, 0, flatness, controls};
  for (size_t i = 0; i < edges.size(); ++i)
    asked.sides.push_back(sideOf(kernel, edges[i], supports[i], orders[i], static_cast<int>(i)));
  for (size_t i = 0; i < guides.size(); ++i)
    asked.guides.push_back(sideOf(kernel, guides[i], 0, 0, static_cast<int>(i)));
  for (size_t i = 0; i < points.size(); i += 3)
  {
    if (!std::isfinite(points[i] + points[i + 1] + points[i + 2]))
      fail(op, invalid_input, "points must be finite");
    asked.points.emplace_back(points[i], points[i + 1], points[i + 2]);
  }
  Built       built = patchFace(asked, method, flatness, controls);
  PatchResult result;
  result.evolution.resultId = store(kernel, built.face);
  int degreeU = 0, degreeV = 0, polesU = 0, polesV = 0;
  if (const auto spline = occ::down_cast<Geom_BSplineSurface>(BRep_Tool::Surface(built.face));
      !spline.IsNull())
  {
    degreeU = spline->UDegree();
    degreeV = spline->VDegree();
    polesU  = spline->NbUPoles();
    polesV  = spline->NbVPoles();
  }
  result.values = {static_cast<double>(built.construction),
                   static_cast<double>(degreeU),
                   static_cast<double>(degreeV),
                   static_cast<double>(polesU),
                   static_cast<double>(polesV)};
  result.values.insert(result.values.end(),
                       built.measured.values.begin(),
                       built.measured.values.end());
  // The caller's numbers are its own edges': gaps stand between those and the guides.
  const size_t gaps =
    built.measured.values.size() - edges.size() * 4 - guides.size() - points.size() / 3;
  result.values.erase(result.values.begin() + 5 + edges.size() * 4,
                      result.values.begin() + 5 + edges.size() * 4 + gaps);
  result.values.insert(result.values.end(), built.squareness.begin(), built.squareness.end());
  result.values.insert(result.values.end(), built.shapeTarget.begin(), built.shapeTarget.end());
  return result;
}
} // namespace

PatchResult patch(OcctKernel&                  kernel,
                  const std::vector<uint32_t>& edges,
                  const std::vector<uint32_t>& supports,
                  const std::vector<int>&      orders,
                  const std::vector<uint32_t>& guides,
                  const std::vector<double>&   points,
                  double                       tolerance,
                  double                       g1Degrees,
                  double                       g2Share,
                  PatchMethod                  method,
                  double                       flatness,
                  const PatchControls&         controls)
{
  try
  {
    return patchAsAsked(kernel,
                        edges,
                        supports,
                        orders,
                        guides,
                        points,
                        tolerance,
                        g1Degrees,
                        g2Share,
                        method,
                        flatness,
                        controls);
  }
  catch (const Failure& failure)
  {
    // A surface that goes over smoothly into a face along an edge meets that edge: where
    // an edge that was only to be met came with its face and no surface was found, the
    // one that is tangent to that face is looked for, and is what was asked if it is
    // found. What is said of the edge is then its tangent plane's miss as well.
    std::vector<int> smoother = orders;
    bool             any      = false;
    for (size_t i = 0; i < smoother.size() && i < supports.size(); ++i)
      if (smoother[i] == 0 && supports[i])
      {
        smoother[i] = 1;
        any         = true;
      }
    if (!any || failure.reason == invalid_input)
      throw;
    try
    {
      return patchAsAsked(kernel,
                          edges,
                          supports,
                          smoother,
                          guides,
                          points,
                          tolerance,
                          g1Degrees,
                          g2Share,
                          method,
                          flatness,
                          controls);
    }
    catch (const Failure&)
    {
    }
    throw failure;
  }
}

TopoDS_Face patchLoop(const std::vector<PatchSide>& sides, double tolerance, double g1Degrees)
{
  Asked asked{{}, {}, {}, {tolerance, g1Degrees, 0.05}, 0, -1, PatchControls{}};
  for (size_t i = 0; i < sides.size(); ++i)
    asked.sides.push_back(
      sideOf(sides[i].edge, sides[i].support, sides[i].order, static_cast<int>(i)));
  return patchFace(asked, PatchMethod::Auto, -1, PatchControls{}).face;
}

} // namespace BRepExt
