// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepSurfaceExt_Tools.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <Geom_BSplineCurve.hxx>
#include <NCollection_HArray1.hxx>

#include <NCollection_Array2.hxx>
#include <algorithm>
#include <cmath>
#include <cstdio>

namespace BRepExt
{
namespace
{
std::vector<gp_Pnt> section(const char*                op,
                            const surfacing::Boundary& a,
                            const surfacing::Boundary& b,
                            double                     f,
                            int                        orderA,
                            int                        orderB,
                            double                     tensionA,
                            double                     tensionB,
                            int                        degree)
{
  const auto   p = a.point(f), q = b.point(f);
  const gp_Vec chord(p, q);
  if (chord.Magnitude() < BRepExt::tolerance::length)
    fail(op,
         unsupported_geometry,
         "the boundaries touch or intersect; bridge requires a separated strip");
  const auto end = [&](const surfacing::Boundary& side, int order, double tension, bool last) {
    auto   j = side.jet(op, f);
    gp_Pnt point;
    gp_Vec tangent;
    side.curve->D1(side.first + f * (side.last - side.first), point, tangent);
    gp_Vec d = chord, dd(0, 0, 0);
    if (order > 0)
    {
      d = j.n.Crossed(tangent).Normalized();
      if (d.Dot(chord) < 0)
        d.Reverse();
      d *= chord.Magnitude() * tension;
    }
    if (order > 1)
      dd = j.n * j.curvature(d, d);
    std::vector<gp_Pnt> poles;
    auto                base = last ? q : p;
    const double        sign = last ? -1 : 1;
    poles.push_back(base);
    poles.emplace_back(base.XYZ() + d.XYZ() * (sign / degree));
    if (degree >= 5)
      poles.emplace_back(base.XYZ() + d.XYZ() * (sign * 2 / degree)
                         + dd.XYZ() / (degree * (degree - 1)));
    if (degree >= 7)
      poles.emplace_back(base.XYZ() + d.XYZ() * (sign * 3 / degree)
                         + dd.XYZ() * (3.0 / (degree * (degree - 1))));
    return poles;
  };
  auto left = end(a, orderA, tensionA, false), right = end(b, orderB, tensionB, true);
  std::vector<gp_Pnt> out = left;
  for (auto it = right.rbegin(); it != right.rend(); ++it)
    out.push_back(*it);
  return out;
}
} // namespace

EvolutionData bridgeSurface(OcctKernel& k,
                            uint32_t    firstEdge,
                            uint32_t    firstFace,
                            uint32_t    secondEdge,
                            uint32_t    secondFace,
                            int         orderA,
                            int         orderB,
                            double      tensionA,
                            double      tensionB,
                            double      tolerance,
                            BridgeShape shape)
{
  const char* op = "bridgeSurface";
  if (orderA < 0 || orderA > 2 || orderB < 0 || orderB > 2 || !std::isfinite(tensionA)
      || !std::isfinite(tensionB) || tensionA <= 0 || tensionB <= 0 || tensionA > 10
      || tensionB > 10 || !std::isfinite(tolerance) || tolerance <= 0)
    fail(op,
         invalid_input,
         "G0/G1/G2, positive tensions up to ten and a positive finite tolerance are required");
  if (shape < BridgeShape::Smooth || shape > BridgeShape::Chamfer)
    fail(op, invalid_input, "shape must be smooth or chamfer");
  if (shape == BridgeShape::Chamfer)
  {
    orderA = 0;
    orderB = 0;
  }
  auto a = surfacing::boundaryOf(k, op, firstEdge, firstFace),
       b = surfacing::boundaryOf(k, op, secondEdge, secondFace);
  if (a.edge.IsSame(b.edge) || a.point(0).Distance(a.point(1)) < BRepExt::tolerance::length
      || b.point(0).Distance(b.point(1)) < BRepExt::tolerance::length)
    fail(op, unsupported_geometry, "two distinct open boundary edges are required");
  if (a.point(0).Distance(b.point(1)) + a.point(1).Distance(b.point(0))
      < a.point(0).Distance(b.point(0)) + a.point(1).Distance(b.point(1)))
    std::swap(b.first, b.last);
  const int degree = 2 * std::max(orderA, orderB) + 3;
  for (int samples = 17; samples <= 257; samples = 2 * samples - 1)
  {
    std::vector<occ::handle<Geom_BSplineCurve>> curves;
    for (int row = 0; row <= degree; ++row)
    {
      occ::handle<NCollection_HArray1<gp_Pnt>> points = new NCollection_HArray1<gp_Pnt>(1, samples);
      occ::handle<NCollection_HArray1<double>> parameters =
        new NCollection_HArray1<double>(1, samples);
      for (int i = 0; i < samples; ++i)
      {
        double f = double(i) / (samples - 1);
        points->SetValue(i + 1,
                         section(op, a, b, f, orderA, orderB, tensionA, tensionB, degree)[row]);
        parameters->SetValue(i + 1, f);
      }
      GeomAPI_Interpolate interpolate(points, parameters, false, 1e-10);
      interpolate.Perform();
      if (!interpolate.IsDone())
        fail(op, no_solution, "longitudinal interpolation of the Hermite blend failed");
      curves.push_back(interpolate.Curve());
    }
    // Interpolation with shared explicit parameters gives every row the same knots.
    const auto                 u = curves.front();
    NCollection_Array2<gp_Pnt> poles(1, u->NbPoles(), 1, degree + 1);
    for (int j = 0; j <= degree; ++j)
    {
      if (curves[j]->NbPoles() != u->NbPoles() || curves[j]->Degree() != u->Degree())
        fail(op, no_solution, "Hermite rows have incompatible knots");
      for (int i = 1; i <= u->NbPoles(); ++i)
        poles.SetValue(i, j + 1, curves[j]->Pole(i));
    }
    NCollection_Array1<double> vk(1, 2);
    vk(1) = 0;
    vk(2) = 1;
    NCollection_Array1<int> vm(1, 2);
    vm.Init(degree + 1);
    occ::handle<Geom_BSplineSurface> surface =
      new Geom_BSplineSurface(poles, u->Knots(), vk, u->Multiplicities(), vm, u->Degree(), degree);
    // Between all interpolation sites, test the exact Hermite control curves. This checks
    // the full surface through the convex hull property, not only its boundary positions.
    double gap = 0;
    for (int i = 0; i < 4 * (samples - 1); ++i)
    {
      double f     = (i + .5) / (4 * (samples - 1));
      auto   exact = section(op, a, b, f, orderA, orderB, tensionA, tensionB, degree);
      for (int j = 0; j <= degree; ++j)
        gap = std::max(gap, curves[j]->Value(f).Distance(exact[j]));
    }
    if (gap > tolerance * .1 && samples < 257)
      continue;
    if (gap > tolerance * .1)
    {
      char text[180];
      std::snprintf(text,
                    sizeof text,
                    "Hermite control curves miss by %.9g (limit %.9g) after 257 longitudinal sites",
                    gap,
                    tolerance * .1);
      fail(op, no_solution, text);
    }
    const std::vector<surfacing::Constraint> targets{{a, {0, 0}, {1, 0}, orderA},
                                                     {b, {0, 1}, {1, 1}, orderB}};
    const auto normal = surfaceJet(op, surface, .5, 0).n.Dot(a.jet(op, .5).n);
    auto       made =
      surfacing::makeFace(op, surface, tolerance, normal < 0 ? TopAbs_REVERSED : TopAbs_FORWARD);
    surfacing::measure(op, made, targets, tolerance);
    return stored(k, made);
  }
  fail(op, no_solution, "Hermite refinement exhausted");
}
} // namespace BRepExt
