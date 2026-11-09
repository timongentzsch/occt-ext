// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepSurfaceExt_DeformUvn.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_History.hxx>
#include <GeomAlgoExt_UniformKnots.hxx>
#include <GeomAlgoExt_PatchFit.hxx>
#include <BRepSurfaceExt_Tools.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Check.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_GTransform.hxx>
#include <gp_GTrsf.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepLib.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BSplCLib.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <GProp_GProps.hxx>
#include <NCollection_Array1.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

namespace BRepExt
{
namespace
{
//! How many poles a target may have along each direction: what a cubic needs for its bends at the
//! tolerance, or the 21 any target may have. A C2 cubic through points of spacing h is off by
//! 5 h^4 max|f(4)| / 384 inside; ends that follow the cubic through four points have a slope off by
//! h^3 max|f(4)| / 4, which moves the first span by 4/27 of h times that. Together about
//! h^4 max|f(4)| / 20. On an arc of radius R and a turn of d apart that is R d^4 / 20, so a turn T
//! over a length L needs T (L / (20 T tolerance))^(1/4) spans. A grid doubled until it fits may
//! have twice as many, with three more poles.
std::pair<int, int> targetBudget(const occ::handle<Geom_BSplineSurface>& s, double tolerance)
{
  const auto along = [&](bool u) {
    const int lines = u ? s->NbVPoles() : s->NbUPoles(), count = u ? s->NbUPoles() : s->NbVPoles();
    double    most = 0;
    for (int line = 1; line <= lines; ++line)
    {
      const auto pole = [&](int k) { return u ? s->Pole(k, line) : s->Pole(line, k); };
      double     turn = 0, length = 0;
      for (int k = 1; k < count; ++k)
      {
        const gp_Vec step(pole(k), pole(k + 1));
        length += step.Magnitude();
        if (k > 1 && step.Magnitude() > gp::Resolution()
            && pole(k - 1).Distance(pole(k)) > gp::Resolution())
          turn += gp_Vec(pole(k - 1), pole(k)).Angle(step);
      }
      if (turn > 0 && length > 0)
        most = std::max(most, turn * std::pow(length / (20 * turn * tolerance), 0.25));
    }
    return std::max(21, int(std::ceil(2 * most)) + 3);
  };
  return {along(true), along(false)};
}

constexpr auto op     = "deformBody";
constexpr int  degree = 5, maxSpans = 16;

struct Side
{
  int   edge, side;
  gp_XY first, last;
};

struct Face
{
  TopoDS_Face               shape;
  occ::handle<Geom_Surface> source;
  double                    a, b, c, d;
  std::vector<Side>         sides;

  gp_Pnt point(double u, double v) const { return source->Value(a + u * (b - a), c + v * (d - c)); }
};

struct Edge
{
  TopoDS_Edge             shape;
  occ::handle<Geom_Curve> source;
  double                  a, b;
  int                     first, last;
};

std::vector<double> sites(const std::vector<double>& knots)
{
  std::vector<double> out(knots.size() - degree - 1);
  for (size_t i = 0; i < out.size(); ++i)
  {
    for (int k = 1; k <= degree; ++k)
      out[i] += knots[i + k] / degree;
  }
  out.front() = 0;
  out.back()  = 1;
  return out;
}

occ::handle<Geom_BSplineCurve> edgeFit(const Edge& edge, const DeformUVN& map, int spans)
{
  auto                       flat = GeomAlgoExt::uniformKnots(degree, spans), at = sites(flat);
  int                        n = int(at.size());
  NCollection_Array1<double> knots(1, spans + 1), sequence(1, int(flat.size())), parameters(1, n);
  NCollection_Array1<int>    mults(1, spans + 1), contacts(1, n);
  NCollection_Array1<gp_Pnt> poles(1, n);
  for (size_t i = 0; i < flat.size(); ++i)
    sequence(int(i) + 1) = flat[i];
  for (int i = 0; i < n; ++i)
  {
    parameters(i + 1) = at[i];
    contacts(i + 1)   = 0;
    poles(i + 1)      = map.point(edge.source->Value(edge.a + at[i] * (edge.b - edge.a)));
  }
  int singular = 0;
  BSplCLib::Interpolate(degree, sequence, parameters, contacts, poles, singular);
  if (singular)
    fail(op, no_solution, "shared edge fit is singular");
  for (int i = 0; i <= spans; ++i)
  {
    knots(i + 1) = double(i) / spans;
    mults(i + 1) = i == 0 || i == spans ? degree + 1 : 1;
  }
  return new Geom_BSplineCurve(poles, knots, mults, degree);
}

int poleIndex(int side, int i, int n)
{
  return side == 0 ? i : side == 1 ? (n - 1) * n + i : side == 2 ? i * n : i * n + n - 1;
}

occ::handle<Geom_BSplineSurface> faceFit(const Face&                                        face,
                                         const DeformUVN&                                   map,
                                         const std::vector<occ::handle<Geom_BSplineCurve>>& curves,
                                         int                                                spans,
                                         double tolerance)
{
  SplineNet net;
  net.degreeU = net.degreeV = degree;
  net.knotsU = net.knotsV = GeomAlgoExt::uniformKnots(degree, spans);
  auto      at            = sites(net.knotsU);
  const int n             = int(at.size());
  net.poles.resize(n * n);
  std::vector<char> fixed(n * n, false);
  // Seed interior poles in world coordinates. They are then fitted to mapped evaluations,
  // never accepted merely because the original control points were warped.
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      net.poles[i * n + j] = map.point(face.point(at[i], at[j])).XYZ();
  for (const auto& side : face.sides)
  {
    auto       curve = curves[side.edge - 1];
    const bool back =
      side.side < 2 ? side.last.Y() < side.first.Y() : side.last.X() < side.first.X();
    for (int i = 0; i < n; ++i)
    {
      const int p  = poleIndex(side.side, i, n);
      net.poles[p] = curve->Pole(back ? n - i : i + 1).XYZ();
      fixed[p]     = true;
    }
  }
  std::vector<FitCondition> conditions;
  for (int i = 1; i < n - 1; ++i)
    for (int j = 1; j < n - 1; ++j)
    {
      auto p = map.point(face.point(at[i], at[j]));
      for (int axis = 1; axis <= 3; ++axis)
      {
        gp_XYZ direction(0, 0, 0);
        direction.SetCoord(axis, 1);
        conditions.push_back(
          {at[i], at[j], {{0, 0, direction}}, 1, p.Coord(axis), 1, false, tolerance * 0.001});
      }
    }
  fitNet(op, net, fixed, conditions, {0, 0}, {});
  return surfaceOf(net);
}

void dimensions(const DeformUVN& map, const std::vector<Face>& faces, gp_Pnt& low, gp_Pnt& high)
{
  low.SetCoord(1e100, 1e100, 1e100);
  high.SetCoord(-1e100, -1e100, -1e100);
  for (const auto& face : faces)
    for (int i = 0; i <= 16; ++i)
      for (int j = 0; j <= 16; ++j)
      {
        auto p = map.local(face.point(double(i) / 16, double(j) / 16));
        for (int axis = 1; axis <= 3; ++axis)
        {
          low.SetCoord(axis, std::min(low.Coord(axis), p.Coord(axis)));
          high.SetCoord(axis, std::max(high.Coord(axis), p.Coord(axis)));
        }
      }
}

double checkMap(const DeformUVN& map, const gp_Pnt& low, const gp_Pnt& high)
{
  // Numerical map Jacobian on an enclosing reference-frame grid. This is a bounded
  // operational check, not a certified global injectivity/error proof.
  const double step    = std::max(1e-6, low.Distance(high) * 1e-5);
  double       minimum = std::numeric_limits<double>::max();
  const double sign    = map.orientation() > 0 ? 1 : -1;
  double       umin = 1e100, umax = -1e100, vmin = 1e100, vmax = -1e100;
  for (int i = 0; i <= 12; ++i)
    for (int j = 0; j <= 12; ++j)
      for (int k = 0; k <= 4; ++k)
      {
        gp_Pnt p(low.X() + (high.X() - low.X()) * i / 12,
                 low.Y() + (high.Y() - low.Y()) * j / 12,
                 low.Z() + (high.Z() - low.Z()) * k / 4);
        auto   uv = map.uv(map.world(p));
        umin      = std::min(umin, uv.X());
        umax      = std::max(umax, uv.X());
        vmin      = std::min(vmin, uv.Y());
        vmax      = std::max(vmax, uv.Y());
        gp_Vec derivatives[3];
        for (int axis = 1; axis <= 3; ++axis)
        {
          auto before = p, after = p;
          before.SetCoord(axis, p.Coord(axis) - step);
          after.SetCoord(axis, p.Coord(axis) + step);
          derivatives[axis - 1] =
            gp_Vec(map.point(map.world(before)), map.point(map.world(after))) / (2 * step);
        }
        double       determinant = derivatives[0].Dot(derivatives[1].Crossed(derivatives[2]));
        const double product =
          derivatives[0].Magnitude() * derivatives[1].Magnitude() * derivatives[2].Magnitude();
        if (!std::isfinite(determinant) || product < 1e-14
            || std::abs(determinant) < product * 1e-5)
          fail(op, no_solution, "UVN map collapses or reaches a normal-offset focal distance");
        if (sign * determinant <= 0)
          fail(op, no_solution, "UVN map folds or reverses beyond a normal-offset focal distance");
        minimum = std::min(minimum, std::abs(determinant) / product);
      }
  map.checkChart(umin, umax, vmin, vmax);
  return minimum;
}

occ::handle<Geom2d_BSplineCurve> pcurve(const Side& side)
{
  NCollection_Array1<gp_Pnt2d> poles(1, 2);
  poles(1) = gp_Pnt2d(side.first);
  poles(2) = gp_Pnt2d(side.last);
  NCollection_Array1<double> knots(1, 2);
  knots(1) = 0;
  knots(2) = 1;
  NCollection_Array1<int> mults(1, 2);
  mults(1) = mults(2) = 2;
  return new Geom2d_BSplineCurve(poles, knots, mults, 1);
}
} // namespace

DeformResult deformBody(OcctKernel&             kernel,
                        uint32_t                bodyId,
                        uint32_t                referenceId,
                        uint32_t                targetId,
                        Vec3                    scale,
                        Vec3                    offset,
                        bool                    swapUV,
                        bool                    mirror,
                        bool                    flipNormal,
                        bool                    reblend,
                        double                  tolerance,
                        const std::vector<int>& hashes,
                        int                     bound)
{
  if (reblend)
    fail(op,
         unsupported_geometry,
         "Reblend unavailable: V1 cannot regenerate fillets or their radii");
  if (!std::isfinite(tolerance) || tolerance < Precision::Confusion())
    fail(op,
         invalid_input,
         "fit tolerance must be finite and at least the kernel precision (1e-7)");
  auto input = get(kernel, bodyId);
  if (input.ShapeType() != TopAbs_FACE && input.ShapeType() != TopAbs_SHELL
      && input.ShapeType() != TopAbs_SOLID)
    fail(op,
         unsupported_geometry,
         "select a whole rectangular sheet, sheet shell, or cuboid solid");
  if (!BRepCheck_Analyzer(input).IsValid())
    fail(op, invalid_input, "moving B-rep is invalid");
  auto reference        = surfacing::faceOf(kernel, op, referenceId),
       target           = surfacing::faceOf(kernel, op, targetId);
  const auto targetType = BRepAdaptor_Surface(target).GetType();
  if (targetType != GeomAbs_Plane && targetType != GeomAbs_Cylinder
      && targetType != GeomAbs_BSplineSurface)
    fail(op,
         unsupported_geometry,
         "V1 targets are planes, open cylinder charts or regular B-splines");
  if (targetType == GeomAbs_BSplineSurface)
  {
    auto spline                   = occ::down_cast<Geom_BSplineSurface>(BRep_Tool::Surface(target));
    const auto [uBudget, vBudget] = targetBudget(spline, tolerance);
    if (spline->NbUPoles() > uBudget || spline->NbVPoles() > vBudget)
      fail(op, unsupported_geometry, "target has more poles than its bends need at this tolerance");
  }
  DeformUVN map(op, reference, target, scale, offset, swapUV, mirror, flipNormal);
  if (scale.z == 0)
    fail(op, invalid_input, "Scale N must be nonzero; collapse is unsupported");
  ShapeSet edgeMap, vertexMap, faceMap;
  TopExp::MapShapes(input, TopAbs_FACE, faceMap);
  TopExp::MapShapes(input, TopAbs_EDGE, edgeMap);
  TopExp::MapShapes(input, TopAbs_VERTEX, vertexMap);
  if (faceMap.Extent() > 16)
    fail(op, unsupported_geometry, "V1 is limited to 16 rectangular sheet faces");
  std::vector<Edge> edges;
  for (const auto& shape : edgeMap)
  {
    auto edge = TopoDS::Edge(shape);
    if (BRep_Tool::Degenerated(edge))
      fail(op, unsupported_geometry, "degenerate edges are unsupported");
    double a, b;
    auto   curve = BRep_Tool::Curve(edge, a, b);
    if (curve.IsNull() || b <= a)
      fail(op, unsupported_geometry, "regular edges with 3D curves required");
    TopoDS_Vertex first, last;
    TopExp::Vertices(edge, first, last);
    if (first.IsNull() || last.IsNull() || first.IsSame(last))
      fail(op, unsupported_geometry, "closed edges require an explicit cut");
    edges.push_back({edge, curve, a, b, vertexMap.FindIndex(first), vertexMap.FindIndex(last)});
  }
  std::vector<Face> faces;
  for (const auto& shape : faceMap)
  {
    auto                face   = TopoDS::Face(shape);
    auto                source = BRep_Tool::Surface(face);
    BRepAdaptor_Surface adapt(face);
    if (adapt.GetType() != GeomAbs_Plane && adapt.GetType() != GeomAbs_BSplineSurface)
      fail(op,
           unsupported_geometry,
           "moving faces must be planar or rectangular B-splines",
           {face});
    if (adapt.GetType() == GeomAbs_BSplineSurface)
    {
      auto spline = surfacing::editable(op, face);
      if (spline->NbUPoles() > 21 || spline->NbVPoles() > 21)
        fail(op,
             unsupported_geometry,
             "moving surface exceeds the V1 21 by 21 pole budget",
             {face});
    }
    if (source->IsUPeriodic() || source->IsVPeriodic() || source->Continuity() < GeomAbs_C2)
      fail(op,
           unsupported_geometry,
           "moving periodic faces or continuity breaks require splitting",
           {face});
    Face f{face, source, 0, 0, 0, 0, {}};
    BRepTools::UVBounds(face, f.a, f.b, f.c, f.d);
    if (f.b <= f.a || f.d <= f.c)
      fail(op, unsupported_geometry, "degenerate moving face domain", {face});
    int wires = 0, mask = 0;
    for (TopExp_Explorer it(face, TopAbs_WIRE); it.More(); it.Next())
      ++wires;
    for (TopExp_Explorer it(face, TopAbs_EDGE); it.More(); it.Next())
    {
      auto edge = TopoDS::Edge(it.Current());
      int  side = surfacing::boundarySide(op, face, edge);
      if (mask & (1 << side))
        fail(op, unsupported_geometry, "split rectangular edges require separate charts", {face});
      mask |= 1 << side;
      double a, b;
      auto   pc         = BRep_Tool::CurveOnSurface(edge, face, a, b);
      auto   normalized = [&](double t) {
        auto p = pc->Value(t);
        return gp_XY((p.X() - f.a) / (f.b - f.a), (p.Y() - f.c) / (f.d - f.c));
      };
      auto        start = normalized(a), end = normalized(b);
      const auto& e = edges[edgeMap.FindIndex(edge) - 1];
      for (int i = 0; i <= 32; ++i)
      {
        const double t  = double(i) / 32;
        auto         uv = start + (end - start) * t;
        if ((normalized(a + t * (b - a)) - uv).Modulus() > 1e-9
            || f.point(uv.X(), uv.Y()).Distance(e.source->Value(e.a + t * (e.b - e.a)))
                 > tolerance * 0.1)
          fail(op,
               unsupported_geometry,
               "boundary pcurves must have affine parameter correspondence",
               {face, edge});
      }
      f.sides.push_back({edgeMap.FindIndex(edge), side, start, end});
    }
    if (wires != 1 || f.sides.size() != 4 || mask != 15)
      fail(op,
           unsupported_geometry,
           "trimmed faces and holes are unsupported; use whole rectangles",
           {face});
    faces.push_back(f);
  }
  gp_Pnt low, high;
  dimensions(map, faces, low, high);
  const bool solid = input.ShapeType() == TopAbs_SOLID;
  if (solid)
  {
    if (faces.size() != 6 || edges.size() != 12 || vertexMap.Extent() != 8)
      fail(op, unsupported_geometry, "whole solids are restricted to six-face planar cuboids");
    for (const auto& f : faces)
      if (BRepAdaptor_Surface(f.shape).GetType() != GeomAbs_Plane)
        fail(op,
             unsupported_geometry,
             "whole solids are restricted to planar cuboids; fillets are unsupported");
    for (const auto& vertex : vertexMap)
    {
      auto p = map.local(BRep_Tool::Pnt(TopoDS::Vertex(vertex)));
      for (int axis = 1; axis <= 3; ++axis)
        if (std::min(std::abs(p.Coord(axis) - low.Coord(axis)),
                     std::abs(p.Coord(axis) - high.Coord(axis)))
            > tolerance)
          fail(op, unsupported_geometry, "cuboid axes must align with the reference plane frame");
    }
  }
  DeformResult result;
  result.minJacobian = checkMap(map, low, high);
  if (BRepAdaptor_Surface(target).GetType() == GeomAbs_Plane)
  {
    // Exact affine case. The standard modifier preserves the shared topology graph;
    // no sampled approximation of surface control points is necessary here.
    gp_GTrsf   transform;
    const auto origin = map.point(gp_Pnt(0, 0, 0));
    for (int column = 1; column <= 3; ++column)
    {
      gp_Pnt axis(0, 0, 0);
      axis.SetCoord(column, 1);
      auto delta = gp_Vec(origin, map.point(axis));
      for (int row = 1; row <= 3; ++row)
        transform.SetValue(row, column, delta.Coord(row));
    }
    for (int row = 1; row <= 3; ++row)
      transform.SetValue(row, 4, origin.Coord(row));
    BRepBuilderAPI_GTransform maker(input, transform, true);
    if (!maker.IsDone() || !BRepCheck_Analyzer(maker.Shape(), true).IsValid())
      fail(op, no_solution, "affine deformation produced an invalid B-rep");
    auto output = maker.Shape();
    if (solid)
    {
      auto body = TopoDS::Solid(output);
      if (!BRepLib::OrientClosedSolid(body))
        fail(op, no_solution, "affine solid orientation failed");
      output = body;
    }
    BRepAlgoAPI_Check interference(output, false, true);
    if (interference.HasErrors() || !interference.IsValid())
      fail(op, no_solution, "affine deformation has self-interference");
    for (const auto& face : faces)
    {
      const auto changed = TopoDS::Face(maker.ModifiedShape(face.shape));
      auto       surface = BRep_Tool::Surface(changed);
      double     a, b, c, d;
      BRepTools::UVBounds(changed, a, b, c, d);
      for (int i = 0; i <= 32; ++i)
        for (int j = 0; j <= 32; ++j)
        {
          const double u = double(i) / 32, v = double(j) / 32;
          result.faceResidual = std::max(
            result.faceResidual,
            map.point(face.point(u, v)).Distance(surface->Value(a + u * (b - a), c + v * (d - c))));
        }
      for (TopExp_Explorer it(changed, TopAbs_EDGE); it.More(); it.Next())
      {
        const auto edge = TopoDS::Edge(it.Current());
        double     first, last, lo, hi;
        auto       curve = BRep_Tool::Curve(edge, first, last);
        auto       pc    = BRep_Tool::CurveOnSurface(edge, changed, lo, hi);
        if (pc.IsNull() || curve.IsNull())
          fail(op, no_solution, "affine result lost a boundary curve");
        for (int i = 0; i <= 64; ++i)
        {
          const double t        = double(i) / 64;
          auto         uv       = pc->Value(lo + t * (hi - lo));
          result.pcurveResidual = std::max(
            result.pcurveResidual,
            curve->Value(first + t * (last - first)).Distance(surface->Value(uv.X(), uv.Y())));
        }
      }
    }
    for (const auto& edge : edges)
    {
      auto   changed = TopoDS::Edge(maker.ModifiedShape(edge.shape));
      double a, b;
      auto   curve = BRep_Tool::Curve(changed, a, b);
      for (int i = 0; i <= 64; ++i)
      {
        const double t      = double(i) / 64;
        result.edgeResidual = std::max(result.edgeResidual,
                                       map.point(edge.source->Value(edge.a + t * (edge.b - edge.a)))
                                         .Distance(curve->Value(a + t * (b - a))));
      }
    }
    if (result.faceResidual > tolerance || result.edgeResidual > tolerance
        || result.pcurveResidual > tolerance)
      fail(op, no_solution, "affine numeric residual exceeds the document tolerance");
    FaceHistory history(op, hashes, bound);
    history.collect(maker, input);
    result.evolution = history.finish(store(kernel, output));
    // Zero spans identifies the exact affine construction in saved diagnostics.
    return result;
  }
  std::vector<occ::handle<Geom_BSplineCurve>>   curves;
  std::vector<occ::handle<Geom_BSplineSurface>> surfaces;
  bool                                          accepted = false;
  for (int spans = 1; spans <= maxSpans && !accepted; spans *= 2)
  {
    curves.clear();
    surfaces.clear();
    result.faceResidual = result.edgeResidual = result.pcurveResidual = result.tangentDegrees = 0;
    result.spans = spans;
    for (const auto& edge : edges)
      curves.push_back(edgeFit(edge, map, spans));
    for (const auto& face : faces)
      surfaces.push_back(faceFit(face, map, curves, spans, tolerance));
    // Offset sample lattice, different from fit sites. Bound resources and fail promptly.
    const int samples = std::max(32, spans * 8);
    for (size_t e = 0; e < edges.size(); ++e)
      for (int i = 0; i <= samples; ++i)
      {
        const double t = double(i) / samples;
        result.edgeResidual =
          std::max(result.edgeResidual,
                   map.point(edges[e].source->Value(edges[e].a + t * (edges[e].b - edges[e].a)))
                     .Distance(curves[e]->Value(t)));
      }
    for (size_t f = 0; f < faces.size(); ++f)
    {
      const auto& face    = faces[f];
      const auto& surface = surfaces[f];
      for (int i = 0; i <= samples; ++i)
        for (int j = 0; j <= samples; ++j)
        {
          const double u = double(i) / samples, v = double(j) / samples;
          result.faceResidual =
            std::max(result.faceResidual,
                     map.point(face.point(u, v)).Distance(surface->Value(u, v)));
          if (i == 0 || j == 0 || i == samples || j == samples)
            continue;
          constexpr double step = 1e-5;
          gp_Pnt           p;
          gp_Vec           du, dv;
          surface->D1(u, v, p, du, dv);
          auto       n  = du.Crossed(dv);
          const auto mu = gp_Vec(map.point(face.point(u - step, v)),
                                 map.point(face.point(u + step, v))),
                     mv = gp_Vec(map.point(face.point(u, v - step)),
                                 map.point(face.point(u, v + step)));
          auto expected = mu.Crossed(mv);
          if (n.Magnitude() < 1e-14 || expected.Magnitude() < 1e-14 || n.Dot(expected) <= 0)
            fail(op, no_solution, "fitted face folds or has singular derivatives");
          result.tangentDegrees = std::max(result.tangentDegrees, n.Angle(expected) * 180 / M_PI);
        }
      for (const auto& side : face.sides)
        for (int i = 0; i <= samples; ++i)
        {
          const double t  = double(i) / samples;
          auto         uv = side.first + (side.last - side.first) * t;
          result.pcurveResidual =
            std::max(result.pcurveResidual,
                     surface->Value(uv.X(), uv.Y()).Distance(curves[side.edge - 1]->Value(t)));
        }
    }
    accepted = result.faceResidual <= tolerance * 0.5 && result.edgeResidual <= tolerance * 0.5
               && result.pcurveResidual <= tolerance * 0.1 && result.tangentDegrees <= 0.1;
  }
  if (!accepted)
  {
    char detail[256];
    std::snprintf(detail,
                  sizeof detail,
                  "fit budget exhausted (16 spans): face %.9g, edge %.9g, pcurve %.9g, tangent "
                  "%.9g degrees; tolerance %.9g, angle limit 0.1",
                  result.faceResidual,
                  result.edgeResidual,
                  result.pcurveResidual,
                  result.tangentDegrees,
                  tolerance);
    fail(op, no_solution, detail);
  }
  // The graph is reconstructed directly: each input edge and vertex has exactly one
  // successor. No sewing, projection onto ambiguous branches, or tolerance inflation.
  BRep_Builder               builder;
  std::vector<TopoDS_Vertex> vertices;
  for (const auto& vertex : vertexMap)
    vertices.push_back(
      BRepBuilderAPI_MakeVertex(map.point(BRep_Tool::Pnt(TopoDS::Vertex(vertex)))).Vertex());
  std::vector<TopoDS_Edge> madeEdges;
  for (size_t e = 0; e < edges.size(); ++e)
  {
    const auto& edge = edges[e];
    auto        made =
      BRepBuilderAPI_MakeEdge(curves[e], vertices[edge.first - 1], vertices[edge.last - 1], 0, 1);
    if (!made.IsDone())
      fail(op, no_solution, "shared edge construction failed");
    madeEdges.push_back(made.Edge());
  }
  TopoDS_Shell shell;
  builder.MakeShell(shell);
  std::vector<TopoDS_Face> madeFaces;
  FaceHistory              history(op, hashes, bound);
  for (size_t f = 0; f < faces.size(); ++f)
  {
    const auto& face = faces[f];
    TopoDS_Face made;
    builder.MakeFace(made, surfaces[f], Precision::Confusion());
    TopoDS_Wire wire;
    builder.MakeWire(wire);
    for (const auto& side : face.sides)
    {
      auto edge = madeEdges[side.edge - 1];
      builder.UpdateEdge(edge, pcurve(side), made, Precision::Confusion());
      builder.Range(edge, made, 0, 1);
    }
    auto forward = face.shape;
    forward.Orientation(TopAbs_FORWARD);
    for (TopExp_Explorer wi(forward, TopAbs_WIRE); wi.More(); wi.Next())
      for (BRepTools_WireExplorer it(TopoDS::Wire(wi.Current()), forward); it.More(); it.Next())
      {
        auto edge = madeEdges[edgeMap.FindIndex(it.Current()) - 1];
        edge.Orientation(it.Current().Orientation());
        builder.Add(wire, edge);
      }
    wire.Closed(true);
    builder.Add(made, wire);
    made.Orientation(face.shape.Orientation());
    madeFaces.push_back(made);
    builder.Add(shell, made);
    NCollection_List<TopoDS_Shape> successor;
    successor.Append(made);
    history.modified(face.shape, successor);
  }
  for (auto& edge : madeEdges)
  {
    builder.SameRange(edge, true);
    builder.SameParameter(edge, true);
  }
  TopoDS_Shape output;
  if (solid)
  {
    shell.Closed(true);
    TopoDS_Solid body;
    builder.MakeSolid(body);
    builder.Add(body, shell);
    if (!BRepLib::OrientClosedSolid(body))
      fail(op, no_solution, "deformed shell cannot bound a solid");
    GProp_GProps props;
    BRepGProp::VolumeProperties(body, props);
    if (props.Mass() <= tolerance * tolerance * tolerance)
      fail(op, no_solution, "deformed solid has no positive volume");
    output = body;
  }
  else
    output =
      input.ShapeType() == TopAbs_FACE ? TopoDS_Shape(madeFaces.front()) : TopoDS_Shape(shell);
  if (!BRepCheck_Analyzer(output, true).IsValid())
    fail(op, no_solution, "deformed B-rep is invalid");
  BRepAlgoAPI_Check interference(output, false, true);
  if (interference.HasErrors() || !interference.IsValid())
    fail(op, no_solution, "deformed B-rep self-intersects or its interference check failed");
  ShapeSet finalEdges, finalVertices;
  TopExp::MapShapes(output, TopAbs_EDGE, finalEdges);
  TopExp::MapShapes(output, TopAbs_VERTEX, finalVertices);
  if (finalEdges.Extent() != edgeMap.Extent() || finalVertices.Extent() != vertexMap.Extent())
    fail(op, no_solution, "shared topology was not preserved");
  result.evolution = history.finish(store(kernel, output));
  return result;
}
} // namespace BRepExt
