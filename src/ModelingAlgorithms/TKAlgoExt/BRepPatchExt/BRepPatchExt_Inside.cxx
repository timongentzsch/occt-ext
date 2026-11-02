// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_Phase.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepPatchExt_Inside.hxx>

#include <BRepAdaptor_Surface.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <TopoDS.hxx>
#include <TopExp_Explorer.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_Line.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <Geom2d_BezierCurve.hxx>
#include <Geom2dConvert.hxx>
#include <Geom2dConvert_BSplineCurveToBezierCurve.hxx>
#include <optional>
#include <BRepTopAdaptor_FClass2d.hxx>
#include <GeomLProp_SLProps.hxx>
#include <Geom_BSplineSurface.hxx>
#include <BRep_Tool.hxx>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace BRepExt
{
namespace
{
/// The places a face is looked at, each way across its parameters.
constexpr int grid = 40;
/// How many times more tightly than its rim a patch may bend inside. A patch that fills a
/// hole fairly bends about as its rim does (a cap of a ball exactly so; the corner between
/// three rounds was measured at up to two and a half times, a hole whose rim is a little off
/// its face at three and a half); the spikes and bulges that this is to catch bend seven to
/// a hundred times more. Between the two, four is a judgement.
constexpr double tightest = 4;

// An inexpensive trim for choosing a second fitted net. The final acceptance still uses
// the actual trim in insideFault. Patch has one wire, but parity also handles inner wires.
class SurveyTrim
{
  std::vector<std::vector<gp_XY>> wires;

public:
  explicit SurveyTrim(const TopoDS_Face& face)
  {
    constexpr int samples = 64;
    for (TopExp_Explorer it(face, TopAbs_WIRE); it.More(); it.Next())
    {
      std::vector<gp_XY> wire;
      for (BRepTools_WireExplorer edge(TopoDS::Wire(it.Current()), face); edge.More(); edge.Next())
      {
        double     from, to;
        const auto curve = BRep_Tool::CurveOnSurface(edge.Current(), face, from, to);
        if (edge.Current().Orientation() == TopAbs_REVERSED)
          std::swap(from, to);
        for (int k = 0; k < samples; ++k)
          wire.push_back(curve->Value(from + (to - from) * k / samples).XY());
      }
      wires.push_back(std::move(wire));
    }
  }

  bool contains(double u, double v) const
  {
    bool inside = false;
    for (const auto& wire : wires)
      for (size_t i = 0; i < wire.size(); ++i)
      {
        const gp_XY &a = wire[i], &b = wire[(i + 1) % wire.size()];
        if ((a.Y() > v) != (b.Y() > v)
            && u < (b.X() - a.X()) * (v - a.Y()) / (b.Y() - a.Y()) + a.X())
          inside = !inside;
      }
    return inside;
  }
};

// Exact straight convex UV trims use half-planes. Curved trims use positive-weight
// Bezier control-hull tubes to certify polygon parity away from the boundary.
// Uncertain points and unsupported curves use OCCT, on the same acceptance grid.
class FinalTrim
{
  BRepTopAdaptor_FClass2d exact;
  std::vector<gp_XY>      polygon;

  struct Segment
  {
    gp_XY  a, b;
    double bound;
  };

  std::vector<Segment> segments;
  double               guard = 0;

  void approximate(const TopoDS_Face& face)
  {
    double u0, u1, v0, v1;
    BRepTools::UVBounds(face, u0, u1, v0, v1);
    const double goal = 1e-4 * std::max(u1 - u0, v1 - v0);
    const auto   append =
      [&](auto&& self, const occ::handle<Geom2d_BezierCurve>& curve, int depth) -> bool {
      const gp_XY  a = curve->StartPoint().XY(), b = curve->EndPoint().XY(), ab = b - a;
      const double squared = ab.SquareModulus();
      if (squared <= 1e-24)
        return false;
      double bound = 0;
      for (int k = 1; k <= curve->NbPoles(); ++k)
      {
        if (curve->Weight(k) <= 0)
          return false;
        const gp_XY to = curve->Pole(k).XY() - a;
        bound = std::max(bound, (to - ab * std::clamp(to.Dot(ab) / squared, 0.0, 1.0)).Modulus());
      }
      if (bound <= goal || depth == 12)
      {
        // Positive rational Bezier weights keep the entire boundary within
        // this control-hull tube, including excursions past either endpoint.
        segments.push_back(
          {a, b, bound + BRepExt::tolerance::uv + 1e-10 * std::max(u1 - u0, v1 - v0)});
        return true;
      }
      auto left  = occ::down_cast<Geom2d_BezierCurve>(curve->Copy());
      auto right = occ::down_cast<Geom2d_BezierCurve>(curve->Copy());
      left->Segment(0, .5);
      right->Segment(.5, 1);
      return self(self, left, depth + 1) && self(self, right, depth + 1);
    };
    try
    {
      for (TopExp_Explorer it(face, TopAbs_WIRE); it.More(); it.Next())
        for (BRepTools_WireExplorer edge(TopoDS::Wire(it.Current()), face); edge.More();
             edge.Next())
        {
          double     from, to;
          const auto curve = BRep_Tool::CurveOnSurface(edge.Current(), face, from, to);
          if (curve.IsNull())
          {
            segments.clear();
            return;
          }
          auto spline =
            Geom2dConvert::CurveToBSplineCurve(new Geom2d_TrimmedCurve(curve, from, to));
          if (edge.Current().Orientation() == TopAbs_REVERSED)
            spline->Reverse();
          Geom2dConvert_BSplineCurveToBezierCurve spans(spline);
          for (int k = 1; k <= spans.NbArcs(); ++k)
            if (!append(append, spans.Arc(k), 0))
            {
              segments.clear();
              return;
            }
        }
    }
    catch (const Standard_Failure&)
    {
      segments.clear();
    }
  }

public:
  explicit FinalTrim(const TopoDS_Face& face)
      : exact(face, BRepExt::tolerance::uv)
  {
    approximate(face);
    int wires = 0;
    for (TopExp_Explorer it(face, TopAbs_WIRE); it.More(); it.Next())
    {
      if (++wires > 1)
      {
        polygon.clear();
        return;
      }
      for (BRepTools_WireExplorer edge(TopoDS::Wire(it.Current()), face); edge.More(); edge.Next())
      {
        double from, to;
        auto   curve = BRep_Tool::CurveOnSurface(edge.Current(), face, from, to);
        if (curve.IsNull())
        {
          polygon.clear();
          return;
        }
        if (edge.Current().Orientation() == TopAbs_REVERSED)
          std::swap(from, to);
        const gp_XY  a = curve->Value(from).XY(), b = curve->Value(to).XY(), ab = b - a;
        const double length = ab.Modulus();
        if (length <= 1e-12)
        {
          polygon.clear();
          return;
        }
        guard = std::max(guard, length * 1e-10);
        while (const auto cut = occ::down_cast<Geom2d_TrimmedCurve>(curve))
          curve = cut->BasisCurve();
        if (const auto spline = occ::down_cast<Geom2d_BSplineCurve>(curve); !spline.IsNull())
        {
          for (int k = 1; k <= spline->NbPoles(); ++k)
          {
            const gp_XY  offset = spline->Pole(k).XY() - a;
            const double along  = offset.Dot(ab) / (length * length);
            if (std::abs(offset.Crossed(ab)) / length > length * 1e-12 || along < -1e-12
                || along > 1 + 1e-12 || spline->Weight(k) <= 0)
            {
              polygon.clear();
              return;
            }
          }
        }
        else if (occ::down_cast<Geom2d_Line>(curve).IsNull())
        {
          polygon.clear();
          return;
        }
        polygon.push_back(a);
      }
    }
    if (polygon.size() < 3)
    {
      polygon.clear();
      return;
    }
    double area = 0;
    for (size_t i = 0; i < polygon.size(); ++i)
      area += polygon[i].Crossed(polygon[(i + 1) % polygon.size()]);
    if (area < 0)
      std::reverse(polygon.begin(), polygon.end());
    for (size_t i = 0; i < polygon.size(); ++i)
    {
      const gp_XY a = polygon[i], ab = polygon[(i + 1) % polygon.size()] - a;
      for (const gp_XY& p : polygon)
        if (ab.Crossed(p - a) < -guard * ab.Modulus())
        {
          polygon.clear();
          return;
        }
    }
  }

  bool contains(double u, double v)
  {
    if (!polygon.empty())
    {
      double nearest = HUGE_VAL;
      for (size_t i = 0; i < polygon.size(); ++i)
      {
        const gp_XY a = polygon[i], ab = polygon[(i + 1) % polygon.size()] - a;
        nearest = std::min(nearest, ab.Crossed(gp_XY(u, v) - a) / ab.Modulus());
      }
      if (nearest > guard)
        return true;
      if (nearest < -guard)
        return false;
    }
    if (!segments.empty())
    {
      bool inside = false, uncertain = false;
      for (const auto& segment : segments)
      {
        const gp_XY ab = segment.b - segment.a, to = gp_XY(u, v) - segment.a;
        if ((to - ab * std::clamp(to.Dot(ab) / ab.SquareModulus(), 0.0, 1.0)).Modulus()
            <= segment.bound)
          uncertain = true;
        if ((segment.a.Y() > v) != (segment.b.Y() > v)
            && u < segment.a.X() + ab.X() * (v - segment.a.Y()) / ab.Y())
          inside = !inside;
      }
      // Outside every tube the polygon and the true trim have the same winding.
      // Near a boundary keep OCCT's exact classifier and its original tolerance.
      if (!uncertain)
        return inside;
    }
    return exact.Perform(gp_Pnt2d(u, v)) == TopAbs_IN;
  }
};
} // namespace

double InsideCurvature::wrongSign() const
{
  return std::min(std::max(0.0, high), std::max(0.0, -low));
}

double InsideCurvature::tightest() const
{
  return std::max(std::abs(low), std::abs(high));
}

InsideCurvature insideCurvature(const TopoDS_Face& face)
{
  const auto surface = BRep_Tool::Surface(face);
  double     u0, u1, v0, v1;
  BRepTools::UVBounds(face, u0, u1, v0, v1);
  SurveyTrim      inside(face);
  InsideCurvature range;
  constexpr int   samples = 100;
  for (int i = 0; i < samples; ++i)
    for (int j = 0; j < samples; ++j)
    {
      const double u = u0 + (u1 - u0) * (i + .5) / samples, v = v0 + (v1 - v0) * (j + .5) / samples;
      if (!inside.contains(u, v))
        continue;
      GeomLProp_SLProps at(surface, u, v, 2, 1e-12);
      if (!at.IsCurvatureDefined())
        continue;
      range.low  = std::min(range.low, at.MinCurvature());
      range.high = std::max(range.high, at.MaxCurvature());
    }
  return range;
}

double insideStep(const TopoDS_Face& face)
{
  const auto surface = occ::down_cast<Geom_BSplineSurface>(BRep_Tool::Surface(face));
  if (surface.IsNull())
    return 0;
  double u0, u1, v0, v1;
  BRepTools::UVBounds(face, u0, u1, v0, v1);
  SurveyTrim    inside(face);
  const auto    in      = [&](double u, double v) { return inside.contains(u, v); };
  double        worst   = 0;
  constexpr int samples = 300;
  for (const bool acrossU : {true, false})
  {
    const int count  = acrossU ? surface->NbUKnots() : surface->NbVKnots(),
              degree = acrossU ? surface->UDegree() : surface->VDegree();
    for (int k = 2; k < count; ++k)
    {
      if ((acrossU ? surface->UMultiplicity(k) : surface->VMultiplicity(k)) < degree - 1)
        continue;
      const double knot = acrossU ? surface->UKnot(k) : surface->VKnot(k);
      for (int part = 0; part < samples; ++part)
      {
        const double along =
          (acrossU ? v0 : u0) + ((acrossU ? v1 - v0 : u1 - u0) * (part + .5) / samples);
        const double u = acrossU ? knot : along, v = acrossU ? along : knot, hu = .03 * (u1 - u0),
                     hv = .03 * (v1 - v0);
        if (!in(u, v) || !in(u - hu, v) || !in(u + hu, v) || !in(u, v - hv) || !in(u, v + hv))
          continue;
        const auto bend = [&](double offset) {
          gp_Pnt p;
          gp_Vec su, sv, suu, svv, suv;
          surface
            ->D2(u + (acrossU ? offset : 0), v + (acrossU ? 0 : offset), p, su, sv, suu, svv, suv);
          const gp_Vec n = su.Crossed(sv).Normalized();
          return acrossU ? n.Dot(suu) / su.SquareMagnitude() : n.Dot(svv) / sv.SquareMagnitude();
        };
        worst = std::max(worst, std::abs(bend(-1e-9 * (acrossU ? u1 - u0 : v1 - v0)) - bend(0)));
      }
    }
  }
  return worst;
}

double tightestInside(const TopoDS_Face& face)
{
  const auto surface = BRep_Tool::Surface(face);
  double     u0, u1, v0, v1;
  BRepTools::UVBounds(face, u0, u1, v0, v1);
  FinalTrim inside(face);
  double    bend = 0;
  for (int i = 0; i < grid; ++i)
    for (int j = 0; j < grid; ++j)
    {
      const double u = u0 + (u1 - u0) * (i + .5) / grid, v = v0 + (v1 - v0) * (j + .5) / grid;
      if (!inside.contains(u, v))
        continue;
      GeomLProp_SLProps at(surface, u, v, 2, 1e-12);
      if (at.IsCurvatureDefined())
        bend = std::max({bend, std::abs(at.MaxCurvature()), std::abs(at.MinCurvature())});
    }
  return bend;
}

std::optional<std::string> insideFault(const TopoDS_Face& face, double rimBend)
{
  PatchPhase phase("interior-check");
  const auto surface = BRep_Tool::Surface(face);
  double     u0, u1, v0, v1;
  BRepTools::UVBounds(face, u0, u1, v0, v1);
  FinalTrim           inside(face);
  std::vector<gp_Vec> normals(grid * grid);
  std::vector<char>   seen(grid * grid, 0);
  double              bend = 0;
  for (int i = 0; i < grid; ++i)
    for (int j = 0; j < grid; ++j)
    {
      const double u = u0 + (u1 - u0) * (i + .5) / grid, v = v0 + (v1 - v0) * (j + .5) / grid;
      if (!inside.contains(u, v))
        continue;
      GeomLProp_SLProps at(surface, u, v, 2, 1e-12);
      if (!at.IsNormalDefined())
        return "its surface has no direction at a place inside the loop";
      normals[i * grid + j] = gp_Vec(at.Normal());
      seen[i * grid + j]    = 1;
      if (at.IsCurvatureDefined())
        bend = std::max({bend, std::abs(at.MaxCurvature()), std::abs(at.MinCurvature())});
    }
  // Folded over: two places side by side that face opposite ways.
  for (int i = 0; i < grid; ++i)
    for (int j = 0; j < grid; ++j)
    {
      if (!seen[i * grid + j])
        continue;
      for (const auto [di, dj] : {std::pair{1, 0}, std::pair{0, 1}})
      {
        if (i + di >= grid || j + dj >= grid || !seen[(i + di) * grid + j + dj])
          continue;
        if (normals[i * grid + j].Dot(normals[(i + di) * grid + j + dj]) < 0)
          return "it folds over inside the loop";
      }
    }
  if (rimBend > 0 && bend > tightest * rimBend)
  {
    char text[160];
    std::snprintf(text,
                  sizeof text,
                  "it bends %.0f times more tightly inside the loop than anything on its rim does",
                  bend / rimBend);
    return text;
  }
  return std::nullopt;
}
} // namespace BRepExt
