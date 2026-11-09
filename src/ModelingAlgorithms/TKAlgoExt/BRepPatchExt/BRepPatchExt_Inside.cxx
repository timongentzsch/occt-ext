// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_Phase.hxx>
#include <GeomAlgoExt_SurfaceJet.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepPatchExt_Inside.hxx>

#include <BRepAdaptor_Surface.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <TopoDS.hxx>
#include <TopExp_Explorer.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom2dAdaptor_Curve.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_Line.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <Geom2d_BezierCurve.hxx>
#include <Geom2dConvert.hxx>
#include <Geom2dConvert_BSplineCurveToBezierCurve.hxx>
#include <optional>
#include <memory>
#include <BRepTopAdaptor_FClass2d.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <GeomLProp_SLProps.hxx>
#include <Geom_BSplineSurface.hxx>
#include <NCollection_Array2.hxx>
#include <PLib.hxx>
#include <BRep_Tool.hxx>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace BRepExt
{
namespace
{
//! The places a face is looked at, each way across its parameters.
constexpr int grid = 40;
//! How many times more tightly than its rim a patch may bend inside. A patch that fills a
//! hole fairly bends about as its rim does (a cap of a ball exactly so; the corner between
//! three rounds was measured at up to two and a half times, a hole whose rim is a little off
//! its face at three and a half); the spikes and bulges that this is to catch bend seven to
//! a hundred times more. Between the two, four is a judgement.
constexpr double tightest = 4;

// A face that crosses from one rim to another over the gaps of an open boundary must turn from
// the tangent plane of the one to that of the other, and may have to bend far more tightly than
// either rim to do so: an elbow from the end of a tube onto the end of a plate does. What the
// rims ask of it there is read from the bridge between them: along each crossing of the face, the
// Hermite curve from the one rim to the other that leaves and meets each square to it in its
// tangent plane, at the speed of the chord between them, and with each face's bend that way
// where the face is to be met to G2 (a curve of degree first + last + 1, a straight line between
// two G0 rims). It holds nothing but what the rims ask (their places, tangent planes and bends,
// which the face has been measured to meet within the limits) and is the same curve in any pose
// and at any scale. It knows nothing of the face's inside, so a bulge, a ripple or a spike that
// leaves the rims as they are leaves it as it is. That the face may bend four times as tightly as
// the bridge does near the same share of the crossing is the judgement above, not a theorem:
// over 76 tube-to-plate and wall-to-plate elbows of three heights, tilts and widths, at rest and
// turned, that met their rims, the face bent at most 2.1 times as tightly as the bridge or its
// rims (against the rims alone, up to 18 times).

//! How much earlier or later along a crossing, as a share of its length, the face may turn than
//! its bridge does. Measured over the same elbows: a face at most 5.7 times as tight as the bridge
//! at the same share, 3.1 times within 0.05 and 2.1 times within this. A judgement.
constexpr double bridgeSlack = 0.1;

//! The places a crossing and its bridge are read at.
constexpr int crossingPlaces = 400;

//! The bridge along the crossing of `surface` at `along`, from the rim at `from` to that at `to`
//! (across parameters): its curvature at shares of its length, in turn.
std::vector<std::pair<double, double>> bridgeBends(const occ::handle<Geom_Surface>& surface,
                                                   const Crossing&                  crossing,
                                                   double                           along,
                                                   double                           from,
                                                   double                           to)
{
  // At each end: the place, the way square to the rim in the tangent plane that the face
  // crosses by, the normal and the face's normal curvature that way.
  struct End
  {
    gp_Pnt place;
    gp_Vec way, normal;
    double bend;
  };

  const auto end = [&](double across) -> std::optional<End> {
    const double u = crossing.acrossV ? along : across, v = crossing.acrossV ? across : along;
    gp_Pnt       place;
    gp_Vec       su, sv, suu, svv, suv;
    surface->D2(u, v, place, su, sv, suu, svv, suv);
    // No bridge from a rim where the face has no tangent plane: the rims alone judge there.
    if (su.Crossed(sv).Magnitude() <= 1e-12 * su.Magnitude() * sv.Magnitude())
      return std::nullopt;
    const gp_Vec normal = su.Crossed(sv).Normalized();
    gp_Vec       way    = normal.Crossed(crossing.acrossV ? su : sv).Normalized();
    if (way.Dot(crossing.acrossV ? sv : su) < 0)
      way.Reverse();
    // The way in the parameters, x su + y sv, for the second fundamental form.
    const double e = su.Dot(su), f = su.Dot(sv), g = sv.Dot(sv), det = e * g - f * f,
                 a = way.Dot(su), b = way.Dot(sv), x = (g * a - f * b) / det,
                 y = (e * b - f * a) / det;
    return End{place,
               way,
               normal,
               x * x * suu.Dot(normal) + 2 * x * y * suv.Dot(normal) + y * y * svv.Dot(normal)};
  };
  const auto start = end(from), finish = end(to);
  if (!start || !finish)
    return {};
  const End    ends[2]   = {*start, *finish};
  const int    orders[2] = {std::clamp(crossing.first, 0, 2), std::clamp(crossing.last, 0, 2)};
  const double chord     = ends[0].place.Distance(ends[1].place);
  NCollection_Array2<double> first(1, 3, 0, orders[0]), last(1, 3, 0, orders[1]);
  for (int side = 0; side < 2; ++side)
  {
    NCollection_Array2<double>& asked  = side ? last : first;
    const End&                  at     = ends[side];
    const gp_XYZ                jet[3] = {at.place.XYZ(),
                                          at.way.XYZ() * chord,
                                          at.normal.XYZ() * (at.bend * chord * chord)};
    for (int k = 0; k <= orders[side]; ++k)
      for (int axis = 1; axis <= 3; ++axis)
        asked(axis, k) = jet[k].Coord(axis);
  }
  const int                              degree = orders[0] + orders[1] + 1;
  NCollection_Array1<double>             coefficients(0, 3 * (degree + 1) - 1);
  std::vector<std::pair<double, double>> out;
  if (!PLib::HermiteInterpolate(3, 0, 1, orders[0], orders[1], first, last, coefficients))
    return out;
  double length = 0, before = 0;
  for (int k = 0; k <= crossingPlaces; ++k)
  {
    double values[9];
    PLib::EvalPolynomial(static_cast<double>(k) / crossingPlaces,
                         2,
                         degree,
                         3,
                         coefficients(0),
                         values[0]);
    const gp_Vec tangent(values[3], values[4], values[5]), second(values[6], values[7], values[8]);
    const double speed = tangent.Magnitude();
    if (k)
      length += (speed + before) / 2 / crossingPlaces;
    before = speed;
    out.push_back(
      {length, speed > 0 ? tangent.Crossed(second).Magnitude() / (speed * speed * speed) : 0});
  }
  if (!(length > 0))
    return {};
  for (auto& place : out)
    place.first /= length;
  return out;
}

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
  TopoDS_Face             ownedFace;
  BRepTopAdaptor_FClass2d exact;
  std::vector<gp_XY>      polygon;

  struct Segment
  {
    gp_XY                           a, b;
    double                          bound;
    occ::handle<Geom2d_BezierCurve> curve;
  };

  std::vector<Segment>                 segments;
  std::vector<std::pair<gp_XY, gp_XY>> joins;
  bool                                 closed         = true;
  double                               guard          = 0;
  double                               enclosureGuard = 0;

  void approximate(const TopoDS_Face& face)
  {
    double u0, u1, v0, v1;
    BRepTools::UVBounds(face, u0, u1, v0, v1);
    const double goal = 1e-4 * std::max(u1 - u0, v1 - v0);
    enclosureGuard    = BRepExt::tolerance::uv + 1e-10 * std::max(u1 - u0, v1 - v0);
    const auto append =
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
        segments.push_back({a, b, bound + enclosureGuard, curve});
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
      {
        const size_t firstSegment = segments.size();
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
          const auto type = Geom2dAdaptor_Curve(curve, from, to).GetType();
          // Offset conversion can approximate at a fixed allowance. It cannot
          // supply a proven enclosure of the actual owned trim.
          if (type == GeomAbs_OffsetCurve || type == GeomAbs_OtherCurve)
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
        const size_t lastSegment = segments.size();
        for (size_t k = firstSegment; k < lastSegment; ++k)
        {
          const gp_XY a = segments[k].b;
          const gp_XY b = segments[k + 1 == lastSegment ? firstSegment : k + 1].a;
          if ((b - a).Modulus() > enclosureGuard)
            closed = false;
          joins.emplace_back(a, b);
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
      : ownedFace(face),
        exact(face, BRepExt::tolerance::uv)
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

  bool boundaryMayIntersect(double u0, double u1, double v0, double v1) const
  {
    if (segments.empty())
      return true;
    for (const auto& segment : segments)
    {
      double from = 0, to = 1;
      for (int axis = 1; axis <= 2 && from <= to; ++axis)
      {
        const double low  = (axis == 1 ? u0 : v0) - segment.bound,
                     high = (axis == 1 ? u1 : v1) + segment.bound, origin = segment.a.Coord(axis),
                     direction = segment.b.Coord(axis) - origin;
        if (direction == 0)
        {
          if (origin < low || origin > high)
            from = 2;
        }
        else
        {
          const double a = (low - origin) / direction, b = (high - origin) / direction;
          from = std::max(from, std::nextafter(std::min(a, b), -HUGE_VAL));
          to   = std::min(to, std::nextafter(std::max(a, b), HUGE_VAL));
        }
      }
      if (from <= to)
        return true;
    }
    // A gap is local uncertainty, even though the stored pcurves do not close.
    // Exclude its entire coordinate box, rather than inventing a closing curve.
    for (const auto& join : joins)
      if (std::max(join.first.X(), join.second.X()) + enclosureGuard >= u0
          && std::min(join.first.X(), join.second.X()) - enclosureGuard <= u1
          && std::max(join.first.Y(), join.second.Y()) + enclosureGuard >= v0
          && std::min(join.first.Y(), join.second.Y()) - enclosureGuard <= v1)
        return true;
    return false;
  }

  std::optional<bool> pointMembership(double u, double v) const
  {
    if (segments.empty())
      return std::nullopt;
    if (!closed)
    {
      // A point needs only the delivered face's own classifier. A nonzero
      // cell separately excludes all boundary and gap enclosures before using
      // this point as its seed. No artificial closing curve is introduced.
      const BRepClass_FaceClassifier classifier(ownedFace, gp_Pnt2d(u, v), BRepExt::tolerance::uv);
      if (classifier.State() == TopAbs_IN)
        return true;
      if (classifier.State() == TopAbs_OUT)
        return false;
      return std::nullopt;
    }
    const gp_XY point(u, v);
    const auto  parity = [&](auto&&                                 self,
                            const occ::handle<Geom2d_BezierCurve>& curve) -> std::optional<bool> {
      const gp_XY  a = curve->StartPoint().XY(), b = curve->EndPoint().XY(), ab = b - a;
      const double squared = ab.SquareModulus();
      const auto   gap     = [&](const gp_XY& p) {
        return (p - a - ab * (squared > 0 ? std::clamp((p - a).Dot(ab) / squared, 0., 1.) : 0.))
          .Modulus();
      };
      double bound = 0;
      for (int k = 1; k <= curve->NbPoles(); ++k)
        bound = std::max(bound, gap(curve->Pole(k).XY()));
      if (gap(point) > std::nextafter(bound + enclosureGuard, HUGE_VAL))
        return (a.Y() > v) != (b.Y() > v) && u < a.X() + ab.X() * (v - a.Y()) / ab.Y();
      if (bound <= enclosureGuard)
        return std::nullopt;
      auto left  = occ::down_cast<Geom2d_BezierCurve>(curve->Copy());
      auto right = occ::down_cast<Geom2d_BezierCurve>(curve->Copy());
      left->Segment(0, .5);
      right->Segment(.5, 1);
      const auto l = self(self, left), r = self(self, right);
      if (!l || !r)
        return std::nullopt;
      return *l != *r;
    };
    bool inside = false;
    for (const auto& segment : segments)
    {
      const gp_XY  ab = segment.b - segment.a, to = point - segment.a;
      const double squared = ab.SquareModulus();
      const double distance =
        (to - ab * (squared > 0 ? std::clamp(to.Dot(ab) / squared, 0., 1.) : 0.)).Modulus();
      const auto crossing =
        distance > std::nextafter(segment.bound, HUGE_VAL)
          ? std::optional<bool>((segment.a.Y() > v) != (segment.b.Y() > v)
                                && u < segment.a.X() + ab.X() * (v - segment.a.Y()) / ab.Y())
          : parity(parity, segment.curve);
      if (!crossing)
        return std::nullopt;
      inside = inside != *crossing;
    }
    for (const auto& join : joins)
    {
      const gp_XY  a = join.first, ab = join.second - a, to = point - a;
      const double squared = ab.SquareModulus();
      if ((to - ab * (squared > 0 ? std::clamp(to.Dot(ab) / squared, 0., 1.) : 0.)).Modulus()
          <= enclosureGuard)
        return std::nullopt;
      if ((a.Y() > v) != (join.second.Y() > v) && u < a.X() + ab.X() * (v - a.Y()) / ab.Y())
        inside = !inside;
    }
    return inside;
  }

  bool containsCell(double u0, double u1, double v0, double v1)
  {
    if ((u0 != u1 || v0 != v1) && boundaryMayIntersect(u0, u1, v0, v1))
      return false;
    return pointMembership(u0 + (u1 - u0) / 2, v0 + (v1 - v0) / 2).value_or(false);
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

std::function<bool(double, double, double, double)> trimmedInteriorCells(const TopoDS_Face& face)
{
  auto trim = std::make_shared<FinalTrim>(face);
  return [trim](double u0, double u1, double v0, double v1) {
    return trim->containsCell(u0, u1, v0, v1);
  };
}

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

double bendingInside(const TopoDS_Face& face)
{
  const auto surface = BRep_Tool::Surface(face);
  double     u0, u1, v0, v1;
  BRepTools::UVBounds(face, u0, u1, v0, v1);
  FinalTrim    inside(face);
  const double cell    = (u1 - u0) * (v1 - v0) / (grid * grid);
  double       bending = 0;
  for (int i = 0; i < grid; ++i)
    for (int j = 0; j < grid; ++j)
    {
      const double u = u0 + (u1 - u0) * (i + .5) / grid, v = v0 + (v1 - v0) * (j + .5) / grid;
      if (!inside.contains(u, v))
        continue;
      GeomLProp_SLProps at(surface, u, v, 2, 1e-12);
      if (!at.IsCurvatureDefined())
        continue;
      const double k1 = at.MaxCurvature(), k2 = at.MinCurvature();
      bending += (k1 * k1 + k2 * k2) * at.D1U().Crossed(at.D1V()).Magnitude() * cell;
    }
  return bending;
}

std::optional<std::string> insideFault(const TopoDS_Face&             face,
                                       double                         rimBend,
                                       const std::optional<Crossing>& crossing,
                                       bool                           resolveNormalRotation)
{
  PatchPhase phase("interior-check");
  const auto surface = BRep_Tool::Surface(face);
  double     u0, u1, v0, v1;
  BRepTools::UVBounds(face, u0, u1, v0, v1);
  FinalTrim           inside(face);
  std::vector<gp_Vec> normals(grid * grid);
  std::vector<char>   seen(grid * grid, 0);
  // Across a crossing, the bend the bridge has near the same share of the crossing at each place
  // looked at, crossing by crossing (see bridgeBends).
  std::vector<double> bridged(grid * grid, 0);
  for (int c = 0; crossing && rimBend > 0 && c < grid; ++c)
  {
    const double along = crossing->acrossV ? u0 + (u1 - u0) * (c + .5) / grid
                                           : v0 + (v1 - v0) * (c + .5) / grid,
                 from = crossing->acrossV ? v0 : u0, to = crossing->acrossV ? v1 : u1;
    const auto bridge = bridgeBends(surface, *crossing, along, from, to);
    if (bridge.empty())
      continue;
    // The face's own crossing, by its length.
    std::vector<double> share(crossingPlaces + 1, 0);
    gp_Pnt              before;
    for (int k = 0; k <= crossingPlaces; ++k)
    {
      const double across = from + (to - from) * k / crossingPlaces;
      const gp_Pnt place =
        surface->Value(crossing->acrossV ? along : across, crossing->acrossV ? across : along);
      share[k] = k ? share[k - 1] + place.Distance(before) : 0;
      before   = place;
    }
    for (int k = 0; k < grid; ++k)
    {
      const double at   = share[(2 * k + 1) * crossingPlaces / (2 * grid)] / share.back();
      double       most = 0;
      for (const auto& [where, bend] : bridge)
        if (std::abs(where - at) <= bridgeSlack)
          most = std::max(most, bend);
      bridged[crossing->acrossV ? c * grid + k : k * grid + c] = most;
    }
  }
  double tighter = 0;
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
      {
        const double here = std::max(std::abs(at.MaxCurvature()), std::abs(at.MinCurvature()));
        tighter           = std::max(tighter, here / std::max(rimBend, bridged[i * grid + j]));
      }
    }
  // Opposite normals at separated places do not establish a fold. Subdivide
  // their parameter segment until the rotation is resolved or rank is lost.
  // This remains a sampled regularity check, not a global injectivity proof.
  const auto resolve =
    [&](auto&& self, double ua, double va, const gp_Vec& na, double ub, double vb, const gp_Vec& nb)
    -> bool {
    if (na.Dot(nb) >= 0)
      return true;
    if (!resolveNormalRotation)
      return false;
    const double u = ua + (ub - ua) / 2, v = va + (vb - va) / 2;
    if ((u == ua && v == va) || (u == ub && v == vb))
      return false;
    if (!inside.contains(u, v))
      return true;
    gp_Vec normal;
    try
    {
      normal = surfaceJet("patch", surface, u, v).n;
    }
    catch (const Failure&)
    {
      return false;
    }
    catch (const Standard_Failure&)
    {
      return false;
    }
    return self(self, ua, va, na, u, v, normal) && self(self, u, v, normal, ub, vb, nb);
  };
  for (int i = 0; i < grid; ++i)
    for (int j = 0; j < grid; ++j)
    {
      if (!seen[i * grid + j])
        continue;
      for (const auto [di, dj] : {std::pair{1, 0}, std::pair{0, 1}})
      {
        if (i + di >= grid || j + dj >= grid || !seen[(i + di) * grid + j + dj])
          continue;
        const double u = u0 + (u1 - u0) * (i + .5) / grid, v = v0 + (v1 - v0) * (j + .5) / grid;
        if (!resolve(resolve,
                     u,
                     v,
                     normals[i * grid + j],
                     u0 + (u1 - u0) * (i + di + .5) / grid,
                     v0 + (v1 - v0) * (j + dj + .5) / grid,
                     normals[(i + di) * grid + j + dj]))
          return "its differential regularity is unresolved inside the loop";
      }
    }
  if (rimBend > 0 && tighter > tightest)
  {
    char text[160];
    std::snprintf(text,
                  sizeof text,
                  crossing
                    ? "it bends %.0f times more tightly inside the loop than its rims and the "
                      "bridge between them do"
                    : "it bends %.0f times more tightly inside the loop than anything on its "
                      "rim does",
                  tighter);
    return text;
  }
  return std::nullopt;
}
} // namespace BRepExt
