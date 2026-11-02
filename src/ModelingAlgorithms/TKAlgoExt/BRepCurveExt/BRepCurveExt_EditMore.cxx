// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepSurfaceExt_DeformUvn.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepAdaptor_CompCurve.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <BSplCLib.hxx>
#include <GeomAPI_ExtremaCurveCurve.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomConvert_CompCurveToBSplineCurve.hxx>
#include <GeomProjLib.hxx>
#include <Geom_Plane.hxx>
#include <GeomConvert.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_Surface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <IntCurvesFace_ShapeIntersector.hxx>
#include <NCollection_Array1.hxx>
#include <Precision.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <algorithm>
#include <cmath>
#include <gp_Lin.hxx>
#include <gp_Pln.hxx>
#include <vector>

// More of the curve tools: Deform with heights (Scale N, Offset N, Mirror, Normal), Rebuild
// with a weight, splines cut to the part an edge uses, and a projection with what is hidden
// left out.

namespace BRepExt
{

// Deform Curve in UVN: a curve that lies on or above a flat reference face is
// carried to a target face. A point at (x, y) of the reference and h above it goes to the
// place of the target at the same shares (u, v) of its parameters, and from there
// h * scale N + offset N along the target's normal. Scale and offset in U and V act on the
// shares. `mirror` turns u into 1 - u (the manual: "Mirror the deformed curves in Z axis";
// that this is the mirror meant is this kernel's reading), `flipNormal` turns the heights to
// the other side of the target. Each edge becomes a spline through mapped points, refined
// until the sampled parameter-correspondence residual is within `tolerance`.
EvolutionData deformCurveUVN(OcctKernel& kernel,
                             uint32_t    wireId,
                             uint32_t    referenceId,
                             uint32_t    targetId,
                             Vec3        scale,
                             Vec3        offset,
                             bool        swapUV,
                             bool        mirror,
                             bool        flipNormal,
                             double      tolerance)
{
  const char* op = "deformCurveUVN";
  if (!std::isfinite(tolerance) || tolerance <= 0)
    fail(op, invalid_input, "tolerance must be positive");
  if (get(kernel, wireId).ShapeType() != TopAbs_WIRE
      || get(kernel, referenceId).ShapeType() != TopAbs_FACE
      || get(kernel, targetId).ShapeType() != TopAbs_FACE)
    fail(op, invalid_input, "a curve, a reference face and a target face are required");
  const auto reference = TopoDS::Face(get(kernel, referenceId)),
             target    = TopoDS::Face(get(kernel, targetId));
  const DeformUVN map(op, reference, target, scale, offset, swapUV, mirror, flipNormal);
  auto            mapped = [&](const gp_Pnt& point) { return map.point(point); };
  // Each edge becomes one spline. Where the edge's own curve is smooth (between the places
  // where it is less than twice differentiable) the carried curve is a spline of degree 5
  // with the edge's own parameters. Its knots are the edge's knots, three times each, so
  // that it may be as smooth there as the edge and no smoother; it is put through carried
  // points, as many as it has control points, at the mean places of its knots (the Greville
  // abscissae). Then the carried curve is held against it at seven places between every two
  // knots, at equal parameters; this is a sampled residual, not a global error bound;
  // where that is over the tolerance a knot is added in the middle and the spline drawn
  // again. It is given up where knots would come nearer than a 1024th of the shortest span.
  constexpr int           degree = 5;
  BRepBuilderAPI_MakeWire result;
  for (BRepTools_WireExplorer it(TopoDS::Wire(get(kernel, wireId))); it.More(); it.Next())
  {
    BRepAdaptor_Curve          curve(it.Current());
    const bool                 back    = it.Current().Orientation() == TopAbs_REVERSED;
    auto                       carried = [&](double t) { return mapped(curve.Value(t)); };
    NCollection_Array1<double> smooth(1, curve.NbIntervals(GeomAbs_C2) + 1),
      spans(1, curve.NbIntervals(GeomAbs_CN) + 1);
    curve.Intervals(smooth, GeomAbs_C2);
    curve.Intervals(spans, GeomAbs_CN);
    GeomConvert_CompCurveToBSplineCurve joined;
    for (int piece = smooth.Lower(); piece < smooth.Upper(); ++piece)
    {
      const double from = smooth(piece), to = smooth(piece + 1);
      // The knots between the two ends, and how often each is taken.
      std::vector<std::pair<double, int>> inner;
      double                              finest = to - from, last = from;
      for (int k = spans.Lower(); k <= spans.Upper(); ++k)
        if (spans(k) > from + BRepExt::tolerance::parameter * (to - from)
            && spans(k) < to - BRepExt::tolerance::parameter * (to - from))
        {
          inner.push_back({spans(k), degree - 2});
          finest = std::min(finest, spans(k) - last);
          last   = spans(k);
        }
      finest = std::min(finest, to - last) / 1024;
      Handle(Geom_BSplineCurve) made;
      while (made.IsNull())
      {
        std::vector<double> flat(degree + 1, from);
        for (const auto& [knot, times] : inner)
          flat.insert(flat.end(), times, knot);
        flat.insert(flat.end(), degree + 1, to);
        const int                  count = int(flat.size()) - degree - 1;
        NCollection_Array1<double> flatKnots(1, int(flat.size())), sites(1, count);
        NCollection_Array1<int>    contact(1, count);
        NCollection_Array1<gp_Pnt> poles(1, count);
        for (size_t i = 0; i < flat.size(); ++i)
          flatKnots(int(i) + 1) = flat[i];
        for (int i = 0; i < count; ++i)
        {
          double sum = 0;
          for (int k = 1; k <= degree; ++k)
            sum += flat[i + k];
          sites(i + 1)   = std::min(to, std::max(from, sum / degree));
          contact(i + 1) = 0;
          poles(i + 1)   = carried(sites(i + 1));
        }
        if (poles(1).Distance(poles(count)) <= Precision::Confusion()
            && poles(1).Distance(poles((count + 1) / 2)) <= Precision::Confusion())
          fail(op, no_solution, "a stretch of the curve is carried to a single point");
        int singular = 0;
        BSplCLib::Interpolate(degree, flatKnots, sites, contact, poles, singular);
        if (singular != 0)
          fail(op, no_solution, "the deformed curve could not be drawn through its points");
        NCollection_Array1<double> knots(1, int(inner.size()) + 2);
        NCollection_Array1<int>    mults(1, int(inner.size()) + 2);
        knots(1) = from;
        mults(1) = degree + 1;
        for (size_t i = 0; i < inner.size(); ++i)
        {
          knots(int(i) + 2) = inner[i].first;
          mults(int(i) + 2) = inner[i].second;
        }
        knots(knots.Upper())                = to;
        mults(mults.Upper())                = degree + 1;
        Handle(Geom_BSplineCurve) candidate = new Geom_BSplineCurve(poles, knots, mults, degree);
        std::vector<std::pair<double, int>> finer;
        bool                                near = true;
        for (int i = knots.Lower(); i < knots.Upper(); ++i)
        {
          bool here = true;
          for (int k = 0; k < 7 && here; ++k)
          {
            const double t = knots(i) + (knots(i + 1) - knots(i)) * (k + 0.5) / 7;
            here           = carried(t).Distance(candidate->Value(t)) <= tolerance;
          }
          if (!here)
          {
            if (knots(i + 1) - knots(i) < finest)
              fail(op,
                   no_solution,
                   "the deformed curve could not be drawn to tolerance: the target bends too "
                   "sharply under it");
            finer.push_back({(knots(i) + knots(i + 1)) / 2, 1});
            near = false;
          }
          if (i + 1 < knots.Upper())
            finer.push_back(inner[i - knots.Lower()]);
        }
        if (near)
          made = candidate;
        else
          inner = std::move(finer);
      }
      // The first piece is taken as it is; a later one joins on where the last ended.
      if (!joined.Add(made, tolerance, true, false, 0))
        fail(op, no_solution, "the deformed pieces of an edge do not join up");
    }
    Handle(Geom_BSplineCurve) whole = joined.BSplineCurve();
    if (whole.IsNull())
      fail(op, no_solution, "the curve is carried to a single point");
    if (back)
      whole->Reverse();
    result.Add(BRepBuilderAPI_MakeEdge(whole).Edge());
    if (!result.IsDone())
      fail(op, no_solution, "the deformed pieces do not join up");
  }
  if (!result.IsDone() || !BRepCheck_Analyzer(result.Wire()).IsValid())
    fail(op, no_solution, "deformed curve is invalid");
  return stored(kernel, result.Wire());
}

// Rebuild with a weight: a spline of the degree and number of control points asked, with
// evenly spaced knots, that is at once near the curve and smooth. It minimises
//   weight * (mean squared distance at 8 samples per control point)
//   + (1 - weight) * (mean squared second difference of the control points)
// with both ends held. At weight 1 it is the nearest such spline in the least-squares sense;
// lower weights let it leave the curve to become smoother. The weight is a convention of this
// code, not derived.
EvolutionData rebuildCurveWeighted(OcctKernel& kernel,
                                   uint32_t    wireId,
                                   int         degree,
                                   int         count,
                                   double      weight)
{
  const char* op    = "rebuildCurveWeighted";
  const auto& shape = get(kernel, wireId);
  if (shape.ShapeType() != TopAbs_WIRE || degree < 1 || degree > 9 || count <= degree || count > 256
      || !std::isfinite(weight) || weight <= 0 || weight > 1)
    fail(op,
         invalid_input,
         "rebuild requires a wire, degree 1-9, degree+1 to 256 points and a weight above 0, at "
         "most 1");
  BRepAdaptor_CompCurve source(TopoDS::Wire(shape), true);
  const double          first = source.FirstParameter(), last = source.LastParameter();
  if (source.Value(first).Distance(source.Value(last)) <= BRepExt::tolerance::exactLength)
    fail(op,
         unsupported_geometry,
         "weighted Explicit Control supports open curves only; use Refit or Points for a closed "
         "curve");
  const int                  spans = count - degree;
  NCollection_Array1<gp_Pnt> poles(1, count);
  NCollection_Array1<double> knots(1, spans + 1);
  NCollection_Array1<int>    mults(1, spans + 1);
  for (int i = 1; i <= count; ++i)
    poles(i) = gp_Pnt(0, 0, 0);
  for (int i = 0; i <= spans; ++i)
  {
    knots(i + 1) = double(i) / spans;
    mults(i + 1) = i == 0 || i == spans ? degree + 1 : 1;
  }
  Handle(Geom_BSplineCurve) basis   = new Geom_BSplineCurve(poles, knots, mults, degree);
  const int                 samples = 8 * count;
  // B(row, col): the share of control point col in the curve at sample row.
  std::vector<std::vector<double>> B(samples + 1, std::vector<double>(count, 0.0));
  for (int col = 0; col < count; ++col)
  {
    basis->SetPole(col + 1, gp_Pnt(1, 0, 0));
    for (int row = 0; row <= samples; ++row)
      B[row][col] = basis->Value(double(row) / samples).X();
    basis->SetPole(col + 1, gp_Pnt(0, 0, 0));
  }
  std::vector<gp_XYZ> at(samples + 1);
  for (int row = 0; row <= samples; ++row)
    at[row] = source.Value(first + (last - first) * row / samples).XYZ();
  // Normal equations over the inner control points; the two ends are known.
  const int                        n = count - 2;
  std::vector<std::vector<double>> m(std::max(n, 0), std::vector<double>(n + 3, 0.0));
  const gp_XYZ                     start = at.front(), end = at.back();
  const double near = weight / (samples + 1), smooth = (1 - weight) / std::max(1, count - 2);
  auto         known = [&](int col) { return col == 0 ? start : end; };
  for (int row = 0; row <= samples; ++row)
  {
    gp_XYZ rest = at[row] - start * B[row][0] - end * B[row][count - 1];
    for (int i = 0; i < n; ++i)
    {
      const double bi = B[row][i + 1];
      if (bi == 0)
        continue;
      for (int j = 0; j < n; ++j)
        m[i][j] += near * bi * B[row][j + 1];
      for (int k = 0; k < 3; ++k)
        m[i][n + k] += near * bi * rest.Coord(k + 1);
    }
  }
  // Second differences P[k-1] - 2 P[k] + P[k+1] for k = 1 .. count - 2.
  for (int k = 1; k + 1 < count; ++k)
  {
    const int    cols[3]   = {k - 1, k, k + 1};
    const double factor[3] = {1, -2, 1};
    gp_XYZ       rest(0, 0, 0);
    for (int e = 0; e < 3; ++e)
      if (cols[e] == 0 || cols[e] == count - 1)
        rest -= known(cols[e]) * factor[e];
    for (int e = 0; e < 3; ++e)
    {
      if (cols[e] == 0 || cols[e] == count - 1)
        continue;
      for (int f = 0; f < 3; ++f)
        if (cols[f] != 0 && cols[f] != count - 1)
          m[cols[e] - 1][cols[f] - 1] += smooth * factor[e] * factor[f];
      for (int c = 0; c < 3; ++c)
        m[cols[e] - 1][n + c] += smooth * factor[e] * rest.Coord(c + 1);
    }
  }
  if (!solveForThree(m, n))
    fail(op, no_solution, "the rebuild's equations are singular");
  poles(1)     = gp_Pnt(start);
  poles(count) = gp_Pnt(end);
  for (int i = 0; i < n; ++i)
    poles(i + 2) = gp_Pnt(m[i][n], m[i][n + 1], m[i][n + 2]);
  Handle(Geom_BSplineCurve) spline = new Geom_BSplineCurve(poles, knots, mults, degree);
  BRepBuilderAPI_MakeWire   output(BRepBuilderAPI_MakeEdge(spline).Edge());
  if (!output.IsDone() || !BRepCheck_Analyzer(output.Wire()).IsValid())
    fail(op, no_solution, "rebuild produced an invalid wire");
  return stored(kernel, output.Wire());
}

// Each spline edge that uses a part of its curve only gets a spline of exactly that part: the
// same shape, with control points of its own.
EvolutionData untrimCurve(OcctKernel& kernel, uint32_t wireId)
{
  const char* op    = "untrimCurve";
  const auto& shape = get(kernel, wireId);
  if (shape.ShapeType() != TopAbs_WIRE)
    fail(op, invalid_input, "a curve is required");
  BRepBuilderAPI_MakeWire output;
  for (BRepTools_WireExplorer it(TopoDS::Wire(shape)); it.More(); it.Next())
  {
    TopoDS_Edge edge = it.Current();
    double      a, b;
    auto        curve = BRep_Tool::Curve(edge, a, b);
    while (!curve.IsNull() && curve->IsKind(STANDARD_TYPE(Geom_TrimmedCurve)))
      curve = Handle(Geom_TrimmedCurve)::DownCast(curve)->BasisCurve();
    const auto spline = Handle(Geom_BSplineCurve)::DownCast(curve);
    if (!spline.IsNull() && !spline->IsPeriodic()
        && (a > spline->FirstParameter() + BRepExt::tolerance::parameter
            || b < spline->LastParameter() - BRepExt::tolerance::parameter))
    {
      auto part = Handle(Geom_BSplineCurve)::DownCast(spline->Copy());
      part->Segment(a, b);
      TopoDS_Edge cut = BRepBuilderAPI_MakeEdge(part).Edge();
      cut.Orientation(edge.Orientation());
      edge = cut;
    }
    output.Add(edge);
    if (!output.IsDone())
      fail(op, no_solution, "the pieces of the curve do not join up");
  }
  return stored(kernel, output.Wire());
}

// The projection of a curve onto a body with what is hidden left out: of each projected piece
// those stretches are kept that have nothing of the body between them and the curve they
// come from.
//
// Everything of the body that a ray from the source can meet lies on the projected curves
// themselves (they are the cut of the body with the sheet of all those rays). So along a
// projected curve what is seen can change only where its ray passes an end of a projected
// curve, a place where one of them turns back along the direction of projection, or a
// crossing of two of them. Each curve is divided at exactly those places, and each division
// asked once, at its middle, by a ray back to the source.
EvolutionData projectCurveVisible(OcctKernel& kernel,
                                  uint32_t    wire,
                                  uint32_t    body,
                                  Vec3        dir,
                                  bool        both)
{
  const char*        op        = "projectCurveVisible";
  const auto         all       = projectCurve(kernel, wire, body, dir, both);
  const TopoDS_Shape projected = get(kernel, all.resultId), source = get(kernel, wire),
                     target = get(kernel, body);
  // The whole projection is only a step on the way: it is not left in the store.
  release(kernel, all.resultId);
  const gp_Dir travel(dir.x, dir.y, dir.z);
  Bnd_Box      bounds;
  BRepBndLib::Add(source, bounds);
  BRepBndLib::Add(target, bounds);
  const double extent = std::sqrt(bounds.SquareExtent()) * 2 + 1;
  // Rays meet faces to OpenCascade's own measure of one point, and stop ten times that short
  // of the projected point, so as not to hit the face it is on; at the curve's end they go on
  // ten times that past it, so that a curve lying on a face hides what is behind that face.
  // What is hidden by less than that much of the body (beside the horizon of a round face)
  // is kept.
  const double                   same = Precision::Confusion(), clear = 10 * same;
  IntCurvesFace_ShapeIntersector rays;
  rays.Load(target, same);
  auto visible = [&](const gp_Pnt& q) {
    const TopoDS_Edge line = BRepBuilderAPI_MakeEdge(gp_Lin(q, travel), -extent, extent).Edge();
    BRepExtrema_DistShapeShape to(line, source);
    if (!to.IsDone() || to.NbSolution() == 0)
      return false;
    // How far back along the direction the source is.
    const double back = gp_Vec(q, to.PointOnShape1(1)).Dot(gp_Vec(travel));
    if (std::abs(back) <= 2 * clear)
      return true;
    const double lo = back < 0 ? back - clear : clear, hi = back < 0 ? -clear : back + clear;
    rays.Perform(gp_Lin(q, travel), lo, hi);
    return rays.NbPnt() == 0;
  };

  struct Piece
  {
    Handle(Geom_Curve)  curve, flat; // as it is, and seen along the direction
    double              a, b, tolerance;
    Bnd_Box             box;
    std::vector<double> cuts;
  };

  struct Corner
  {
    gp_Pnt at;
    double tolerance;
  };

  std::vector<Piece>       pieces;
  std::vector<Corner>      corners;
  const Handle(Geom_Plane) screen    = new Geom_Plane(gp_Pnt(0, 0, 0), travel);
  auto                     flattened = [&](const gp_Pnt& p) {
    const gp_Vec r(gp_Pnt(0, 0, 0), p);
    return gp_Pnt((r - gp_Vec(travel) * r.Dot(gp_Vec(travel))).XYZ());
  };
  for (TopExp_Explorer it(projected, TopAbs_EDGE); it.More(); it.Next())
  {
    const TopoDS_Edge edge = TopoDS::Edge(it.Current());
    Piece             piece;
    piece.curve = BRep_Tool::Curve(edge, piece.a, piece.b);
    if (piece.curve.IsNull())
      continue;
    piece.tolerance = BRep_Tool::Tolerance(edge);
    BRepBndLib::Add(edge, piece.box);
    piece.box.Enlarge(piece.tolerance);
    for (TopExp_Explorer v(edge, TopAbs_VERTEX); v.More(); v.Next())
      corners.push_back({BRep_Tool::Pnt(TopoDS::Vertex(v.Current())),
                         BRep_Tool::Tolerance(TopoDS::Vertex(v.Current()))});
    // Where the curve turns back along the direction its tangent is along it: the part of
    // the tangent across the direction passes through nothing and comes out reversed. It
    // is looked for between places a few to each span of the curve, more for a higher
    // degree, and closed in on by halving.
    GeomAdaptor_Curve shape(piece.curve, piece.a, piece.b);
    if (shape.GetType() != GeomAbs_Line)
    {
      auto across = [&](double t) {
        gp_Pnt p;
        gp_Vec d;
        piece.curve->D1(t, p, d);
        return d.Crossed(gp_Vec(travel));
      };
      const bool spline =
        shape.GetType() == GeomAbs_BSplineCurve || shape.GetType() == GeomAbs_BezierCurve;
      const int                  steps = spline ? 4 * (shape.Degree() + 1) : 32;
      NCollection_Array1<double> spans(1, shape.NbIntervals(GeomAbs_CN) + 1);
      shape.Intervals(spans, GeomAbs_CN);
      for (int k = spans.Lower(); k < spans.Upper(); ++k)
      {
        gp_Vec before = across(spans(k));
        for (int i = 1; i <= steps; ++i)
        {
          double lo        = spans(k) + (spans(k + 1) - spans(k)) * (i - 1) / steps,
                 hi        = spans(k) + (spans(k + 1) - spans(k)) * i / steps;
          const gp_Vec low = before, after = across(hi);
          before = after;
          if (low.Dot(after) >= 0)
            continue;
          for (int n = 0; n < 60; ++n)
          {
            const double mid                     = (lo + hi) / 2;
            (across(mid).Dot(low) > 0 ? lo : hi) = mid;
          }
          piece.cuts.push_back((lo + hi) / 2);
          corners.push_back({piece.curve->Value((lo + hi) / 2), same});
        }
      }
    }
    // A curve that runs along the direction is one point seen along it: it is seen or
    // hidden as a whole.
    const gp_Pnt first   = flattened(piece.curve->Value(piece.a));
    bool         upright = true;
    for (int i = 1; i <= 4 && upright; ++i)
      upright = flattened(piece.curve->Value(piece.a + (piece.b - piece.a) * i / 4)).Distance(first)
                <= piece.tolerance + same;
    if (!upright)
      piece.flat = GeomProjLib::ProjectOnPlane(new Geom_TrimmedCurve(piece.curve, piece.a, piece.b),
                                               screen,
                                               travel,
                                               true);
    pieces.push_back(piece);
  }
  for (size_t i = 0; i < pieces.size(); ++i)
  {
    Piece& piece = pieces[i];
    if (piece.flat.IsNull())
      continue;
    // Where the ray through a corner passes this curve.
    for (const auto& corner : corners)
    {
      GeomAPI_ProjectPointOnCurve feet(flattened(corner.at), piece.flat, piece.a, piece.b);
      for (int k = 1; k <= feet.NbPoints(); ++k)
        if (feet.Distance(k) <= corner.tolerance + piece.tolerance)
          piece.cuts.push_back(feet.Parameter(k));
    }
    // Where another curve crosses it (curves of two bodies that reach into each other).
    for (size_t j = i + 1; j < pieces.size(); ++j)
    {
      Piece& other = pieces[j];
      if (piece.box.IsOut(other.box))
        continue;
      GeomAPI_ExtremaCurveCurve nearest(piece.curve,
                                        other.curve,
                                        piece.a,
                                        piece.b,
                                        other.a,
                                        other.b);
      for (int k = 1; k <= nearest.NbExtrema(); ++k)
      {
        if (nearest.Distance(k) > piece.tolerance + other.tolerance)
          continue;
        double u, v;
        nearest.Parameters(k, u, v);
        piece.cuts.push_back(u);
        other.cuts.push_back(v);
      }
    }
  }
  BRep_Builder    builder;
  TopoDS_Compound result;
  builder.MakeCompound(result);
  for (auto& piece : pieces)
  {
    const auto&  curve = piece.curve;
    const double a = piece.a, b = piece.b;
    auto         emit = [&](double from, double to) {
      if (to - from <= BRepExt::tolerance::parameterMerge * (b - a)
          || curve->Value(from).Distance(curve->Value((from + to) / 2)) <= Precision::Confusion())
        return;
      BRepBuilderAPI_MakeEdge edge(curve, from, to);
      if (edge.IsDone())
        builder.Add(result, edge.Edge());
    };
    std::vector<double> stops = {a};
    std::sort(piece.cuts.begin(), piece.cuts.end());
    for (const double cut : piece.cuts)
      if (cut - stops.back() > BRepExt::tolerance::parameterMerge * (b - a)
          && b - cut > BRepExt::tolerance::parameterMerge * (b - a))
        stops.push_back(cut);
    stops.push_back(b);
    double from = a;
    bool   open = false;
    for (size_t i = 0; i + 1 < stops.size(); ++i)
    {
      const bool seen = visible(curve->Value((stops[i] + stops[i + 1]) / 2));
      if (seen && !open)
        from = stops[i];
      if (!seen && open)
        emit(from, stops[i]);
      open = seen;
    }
    if (open)
      emit(from, b);
  }
  if (!TopExp_Explorer(result, TopAbs_EDGE).More())
    fail(op, no_solution, "nothing of the projection is in sight of the curve");
  return stored(kernel, result);
}
} // namespace BRepExt
