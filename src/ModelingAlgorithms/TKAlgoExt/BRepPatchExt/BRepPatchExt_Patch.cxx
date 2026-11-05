// Copyright (c) 2026 Timon Gentzsch
//
// This file is part of Open CASCADE Technology software library.
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
// distribution for complete text of the license and disclaimer of any warranty.
//
// Alternatively, this file may be used under the terms of Open CASCADE
// commercial license or contractual agreement.

#include <GeomAlgoExt_Phase.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <GeomAlgoExt_Cancel.hxx>
#include <BRepPatchExt_Side.hxx>
#include <BRepAlgoExt_PatchFlow.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>
#include <BRepAlgoExt_Tolerance.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepAlgoExt_PatchLoop.hxx>
#include <GeomAlgoExt_PatchFit.hxx>
#include <GeomAlgoExt_PatchCompact.hxx>
#include <BRepPatchExt_Inside.hxx>
#include <BRepPatchExt_Chart.hxx>
#include <GeomAlgoExt_PatchRim.hxx>
#include <GeomAlgoExt_SurfaceJet.hxx>

#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepLib.hxx>
#include <BRepTools.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <Bnd_Box.hxx>
#include <BRepFill_Filling.hxx>
#include <BRep_Builder.hxx>
#include <BSplCLib.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <Geom2dConvert.hxx>
#include <Geom2dAPI_Interpolate.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomConvert.hxx>
#include <GeomConvert_BSplineCurveToBezierCurve.hxx>
#include <GeomConvert_ApproxCurve.hxx>
#include <GeomConvert_CompCurveToBSplineCurve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_Line.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Plane.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Pln.hxx>
#include <gp_Trsf.hxx>
#include <math_Gauss.hxx>
#include <math_Matrix.hxx>
#include <math_Vector.hxx>

#include <algorithm>
#include <deque>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <limits>
#include <optional>
#include <tuple>
#include <unordered_map>

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
// The last two are fitted by patch_fit. `method` 2 is OpenCascade's plate instead, kept to
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

//! The corners that a side is given in the polygon that stands for the loop.
constexpr int polygonCorners = 32;

std::vector<gp_Pnt> polygon(const std::vector<Side>& sides)
{
  std::vector<gp_Pnt> points;
  for (const Side& side : sides)
    for (int i = 0; i < polygonCorners; ++i)
      points.push_back(side.at(static_cast<double>(i) / polygonCorners));
  return points;
}

//! The least area a loop encloses seen along its mean direction, as a share of the square of
//! its length (a circle has 0.08): below it the loop is a line there and back.
constexpr double leastArea = 1e-6;

//! Below this share of its size, the difference between the longer and the shorter axis of
//! the area a loop encloses says nothing of a direction: the 32 chords that stand for a side
//! are that far from its curve.
constexpr double roundLoop = 1e-3;

//! Two corners are as far from the middle when they differ by less than this share: the
//! rounding of the sums they are measured from, with room.
constexpr double sameCorner = BRepExt::tolerance::exactLength;

Frame meanPlane(const std::vector<gp_Pnt>& points, double perimeter)
{
  gp_XYZ centre(0, 0, 0), area(0, 0, 0);
  for (const gp_Pnt& point : points)
    centre += point.XYZ() / static_cast<double>(points.size());
  for (size_t i = 0; i < points.size(); ++i)
    area += (points[i].XYZ() - centre).Crossed(points[(i + 1) % points.size()].XYZ() - centre) / 2;
  if (area.Modulus() < leastArea * perimeter * perimeter)
    fail(op, unsupported_geometry, "the boundary encloses no area seen from any side");
  const gp_Dir z(area);
  // The two directions in the plane are the loop's own, so that a model turned in space gets
  // the net it had, turned with it: the longer axis of the area the loop encloses; where
  // that has none (a square, a loop that is the same after a third of a turn) the way its
  // longest side runs, which lays the net along the sides of a square; and for one closed
  // curve the way to its one corner. Any two directions do to measure in.
  const gp_Vec e1(std::abs(z.X()) < 0.9 ? z.Crossed(gp_Dir(1, 0, 0)) : z.Crossed(gp_Dir(0, 1, 0)));
  const gp_Vec e2 = gp_Vec(z).Crossed(e1);
  std::vector<gp_XY> flat;
  for (const gp_Pnt& point : points)
    flat.emplace_back(gp_Vec(gp_Pnt(centre), point).Dot(e1), gp_Vec(gp_Pnt(centre), point).Dot(e2));
  double a = 0, cx = 0, cy = 0, xx = 0, yy = 0, xy = 0;
  for (size_t i = 0; i < flat.size(); ++i)
  {
    const gp_XY &p = flat[i], &q = flat[(i + 1) % flat.size()];
    const double cross = p.Crossed(q);
    a += cross / 2;
    cx += cross * (p.X() + q.X()) / 6;
    cy += cross * (p.Y() + q.Y()) / 6;
    xx += cross * (p.X() * p.X() + p.X() * q.X() + q.X() * q.X()) / 12;
    yy += cross * (p.Y() * p.Y() + p.Y() * q.Y() + q.Y() * q.Y()) / 12;
    xy += cross * (p.X() * q.Y() + 2 * p.X() * p.Y() + 2 * q.X() * q.Y() + q.X() * p.Y()) / 24;
  }
  cx /= a;
  cy /= a;
  xx -= a * cx * cx;
  yy -= a * cy * cy;
  xy -= a * cx * cy;
  // The farthest corner from the middle of the area, and of corners as far (as on a circle)
  // the one between the longest sides: nothing of the order the edges were given in, unless
  // the sides are as long too (a loop the same after a third of a turn), where the first of
  // them round the loop is kept and not the one their lengths' rounding favours at this size.
  const size_t sides = points.size() / polygonCorners;
  const auto   run   = [&](size_t side) {
    double length = 0;
    for (int i = 0; i < polygonCorners; ++i)
      length += points[side * polygonCorners + i].Distance(
        points[(side * polygonCorners + i + 1) % points.size()]);
    return length;
  };
  gp_XY  corner(1, 0);
  double far = -1, beside = 0;
  for (size_t s = 0; s < sides; ++s)
  {
    const gp_XY  at       = flat[s * polygonCorners] - gp_XY(cx, cy);
    const double distance = at.Modulus(), lengths = run(s) + run((s + sides - 1) % sides);
    if (distance > far * (1 + sameCorner)
        || (distance > far * (1 - sameCorner) && lengths > beside * (1 + sameCorner)))
    {
      far    = distance;
      beside = lengths;
      corner = at;
    }
  }
  // Across the way to the corner, for one closed curve: no better than along it, and the
  // way such a curve was laid out before the directions were the loop's own, when it began
  // on the world's x axis.
  gp_XY  axis(-corner.Y(), corner.X());
  double longest = 0;
  for (size_t s = 0; s < sides; ++s)
  {
    const gp_XY chord = flat[(s + 1) % sides * polygonCorners] - flat[s * polygonCorners];
    if (chord.Modulus() > longest * (1 + sameCorner) && chord.Modulus() > sameCorner * far)
    {
      longest = chord.Modulus();
      axis    = chord.Dot(corner) < 0 ? chord.Reversed() : chord;
    }
  }
  axis.Normalize();
  if (std::hypot(xx - yy, 2 * xy) > roundLoop * (xx + yy))
  {
    const double turn = std::atan2(2 * xy, xx - yy) / 2;
    axis.SetCoord(std::cos(turn), std::sin(turn));
    if (axis.Dot(corner) < 0)
      axis.Reverse();
  }
  const gp_Vec x = (e1 * axis.X() + e2 * axis.Y()).Normalized();
  return {gp_Pnt(centre), x, gp_Vec(z).Crossed(x), gp_Vec(z)};
}

//! Whether the segments a-b and c-d cross. A point counts as on one side of a segment only where
//! the sign of the determinant that says so is beyond what rounding its two products and their
//! difference can change (Shewchuk, "Adaptive precision floating-point arithmetic and fast robust
//! geometric predicates", 1997: the bound ccwerrboundA = (3 + 16 e) e, e = 2^-53, of orient2d);
//! where it is not, the point is on the segment's line as far as the numbers tell, and the
//! segments are not taken to cross there. Collinear pieces of one straight side, whose
//! determinants are rounding alone, otherwise crossed or not by the size of the model.
bool crosses(const gp_XY& a, const gp_XY& b, const gp_XY& c, const gp_XY& d)
{
  auto side = [](const gp_XY& p, const gp_XY& q, const gp_XY& r) {
    constexpr double e       = std::numeric_limits<double>::epsilon() / 2;
    const double     left    = (q.X() - p.X()) * (r.Y() - p.Y()),
                 right       = (q.Y() - p.Y()) * (r.X() - p.X());
    const double determinant = left - right,
                 bound       = (3 + 16 * e) * e * (std::abs(left) + std::abs(right));
    return determinant > bound ? 1 : determinant < -bound ? -1 : 0;
  };
  return side(a, b, c) * side(a, b, d) < 0 && side(c, d, a) * side(c, d, b) < 0;
}

//! The area a polygon encloses, negative when it runs clockwise. Fails where it crosses itself.
double enclosed(const std::vector<gp_XY>& loop)
{
  double area = 0;
  for (size_t i = 0; i < loop.size(); ++i)
  {
    const gp_XY &a = loop[i], &b = loop[(i + 1) % loop.size()];
    area += a.Crossed(b) / 2;
    for (size_t j = i + 2; j < loop.size() && (i > 0 || j + 1 < loop.size()); ++j)
      if (crosses(a, b, loop[j], loop[(j + 1) % loop.size()]))
        fail(op,
             unsupported_geometry,
             "the boundary crosses itself, or folds over itself seen along its mean direction");
  }
  return area;
}

//! Whether `point` is no further than `slack` from the side of a polygon from `from` to `to`.
bool near(const gp_XY& from, const gp_XY& to, const gp_XY& point, double slack)
{
  const gp_XY edge = to - from, away = point - from;
  return (away - edge * std::clamp(away.Dot(edge) / edge.SquareModulus(), 0.0, 1.0)).Modulus()
         <= slack;
}

//! Where the side of a polygon from `from` to `to` crosses the line at height `y`, if it does.
std::optional<double> crossing(const gp_XY& from, const gp_XY& to, double y)
{
  if ((to.Y() > y) == (from.Y() > y))
    return std::nullopt;
  return from.X() + (to.X() - from.X()) * (y - from.Y()) / (to.Y() - from.Y());
}

//! Whether a point is inside a polygon, or no further outside it than `slack`: the polygon
//! stands for curves, which a point on them may be a little outside of.
bool inside(const std::vector<gp_XY>& loop, const gp_XY& point, double slack)
{
  bool in = false;
  for (size_t i = 0, j = loop.size() - 1; i < loop.size(); j = i++)
  {
    if (near(loop[j], loop[i], point, slack))
      return true;
    if (const auto x = crossing(loop[j], loop[i], point.Y()); x && point.X() < *x)
      in = !in;
  }
  return in;
}

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

//! What the change of the bend counts for against the bend, with the loop's size as the unit
//! of length: in full where a bend is asked or wished for along the rim, and by these shares
//! where the rim is only tangent to its faces, and where curves are only met. With more of
//! it the surface is pulled flat between its edges; with none, a point inside the loop stands
//! in a spike. Measured for the tangent rim on the corner of three rounds of 4, whose own
//! surface is a ball, as the distance from that ball at tolerances of 1e-3 and 1e-6: with
//! 0.03, 0.066 and 0.065; with 0.003, 0.024 and 0.068; with none, 0.010 and 0.069. For curves
//! that are only met, 0.03 is the first author's, from a round hole in a sphere; with 0.003
//! an L of six sides on a curved height is missed by 2.6e-6 on every net where 0.03 meets a
//! millionth on the fourth.
constexpr double variationTangent = 0.003, variationMet = 0.03;

//! The fits of a net whose bend is asked along an edge, the first of them without the bend:
//! what is asked is taken about the surface before, and over a steep rim (a ball cut at a
//! quarter of its radius) the third still left the bend 6 to 8 % off where the fifth meets it.
constexpr int askedPasses = 5;

//! The same where the bend is only wished for: nearer than the third is not needed, and on a
//! ball's cap the third is nearer the ball than the fifth.
constexpr int wishedPasses = 3;

//! How often the places that a fit misses between those it was asked are added to them.
constexpr int placeRounds = 3;

//! How near its value a condition of the fit is brought, as a share of the limit that the
//! finished surface is held to. The fit reads at its own places and the surface is measured
//! between them, where it is a little further off: a quarter leaves room for that.
constexpr double fitShare = 0.25;

//! The same for the places between, which the fit only reads: half, since the surface is
//! measured at yet other places, and a miss that is at its largest between two places of the
//! fit is found there.
constexpr double watchShare = 0.5;

//! The speed a surface over a chart is to leave its rim at, across it, and the band about it
//! that is met: over a chart that keeps the rim's lengths (proxyChart), a surface that the chart
//! is conformal to leaves the rim as fast as it runs along it, at one; the band keeps it from
//! running off to nothing or far off, which meet every tangent plane and are no parametrisation.
//! The band is a judgement, not derived: the chart is the proxy's, not the fill's.
constexpr double leaveSpeed = 1, leaveBand = 0.5;

enum class GraphFit
{
  none,
  transverse,
  complete
};

//! The conditions that say so. For a height graph in the normalized mean-plane frame,
//! tangency determines its spatial derivatives exactly, making G2 linear as well. Otherwise
//! G2 is taken about `current`, and left out while there is none.
std::vector<FitCondition> conditions(const std::vector<Target>&              targets,
                                     const occ::handle<Geom_BSplineSurface>& current,
                                     double                                  size,
                                     const Limits&                           limits,
                                     std::vector<int>&                       owner,
                                     GraphFit                                graph)
{
  std::vector<FitCondition> out;
  owner.clear();
  for (const Target& target : targets)
  {
    // Whose every condition is: all that this place adds below.
    struct Own
    {
      std::vector<int>&                owner;
      const std::vector<FitCondition>& out;
      int                              target;

      ~Own() { owner.resize(out.size(), target); }
    } own{owner, out, static_cast<int>(&target - targets.data())};

    const double share = target.watch ? watchShare : fitShare;
    const double near = share * limits.g0, flat = share * std::sin(limits.g1 * M_PI / 180);
    for (int c = 1; c <= 3; ++c)
    {
      gp_XYZ axis(0, 0, 0);
      axis.SetCoord(c, 1);
      out.push_back({target.uv.X(),
                     target.uv.Y(),
                     {{0, 0, axis}},
                     1,
                     axis.Dot(target.point.XYZ()),
                     target.share / (size * size),
                     target.soft,
                     near,
                     target.watch});
    }
    if (target.perpendicular)
    {
      const gp_XYZ direction = target.perpendicular->Normalized().XYZ();
      out.push_back({target.uv.X(),
                     target.uv.Y(),
                     {{1, 0, direction * target.across.X()}, {0, 1, direction * target.across.Y()}},
                     2,
                     0,
                     target.share * tangentWeight,
                     false,
                     1e-7});
    }
    if (target.leave)
    {
      const gp_XYZ direction = target.leave->XYZ();
      out.push_back({target.uv.X(),
                     target.uv.Y(),
                     {{1, 0, direction * target.across.X()}, {0, 1, direction * target.across.Y()}},
                     2,
                     leaveSpeed,
                     target.share * tangentWeight,
                     false,
                     leaveBand,
                     !target.leaveAsked});
    }
    if (target.order < 1)
      continue;
    const SurfaceJet& jet = *target.jet;
    const gp_XYZ      n   = jet.n.XYZ();
    const double      du = target.across.X(), dv = target.across.Y();
    // One tangent plane with the neighbour: both directions of the surface in it. Along
    // the edge that follows from holding the edge, and is asked all the same: measured on
    // a hole across a round, the direction across alone left the fit short of a millionth
    // on every net, and both met it on the second.
    out.push_back({target.uv.X(),
                   target.uv.Y(),
                   {{1, 0, n}},
                   1,
                   0,
                   target.share * tangentWeight,
                   false,
                   flat,
                   target.watch});
    out.push_back({target.uv.X(),
                   target.uv.Y(),
                   {{0, 1, n}},
                   1,
                   0,
                   target.share * tangentWeight,
                   false,
                   flat,
                   target.watch});
    const bool wished = target.order < 2 && target.wish;
    if ((target.order < 2 && !wished)
        || (graph == GraphFit::none && (current.IsNull() || target.watch)))
      continue;
    if (graph != GraphFit::none)
    {
      // The support normal fixes the graph's first derivatives, so evaluating
      // its second fundamental form is exact and linear before fitting. Start
      // with transverse curvature; if measured G2 still fails, the complete
      // form also constrains tangent/twist errors from approximate rim fits.
      const double               diagonal = std::sqrt(0.5);
      const std::array<gp_XY, 3> directions =
        graph == GraphFit::complete ? std::array<gp_XY, 3>{{{1, 0}, {0, 1}, {diagonal, diagonal}}}
                                    : std::array<gp_XY, 3>{{{du, dv}, {}, {}}};
      for (int i = 0; i < (graph == GraphFit::complete ? 3 : 1); ++i)
      {
        const double u = directions[i].X(), v = directions[i].Y();
        const gp_Vec tangent(u, v, -(n.X() * u + n.Y() * v) / n.Z());
        const double bend = jet.curvature(tangent, tangent);
        out.push_back({target.uv.X(),
                       target.uv.Y(),
                       {{2, 0, n * (u * u)}, {1, 1, n * (2 * u * v)}, {0, 2, n * (v * v)}},
                       3,
                       bend,
                       target.share * size * size * bendWeight,
                       wished,
                       share * limits.g2
                         * std::max({std::abs(bend),
                                     target.easing * tangent.SquareMagnitude(),
                                     tangent.SquareMagnitude() * flatShare / size}),
                       target.watch});
      }
      continue;
    }
    // With T the way across in space, the bend across is n . S_dd, to equal II(T, T) of the
    // neighbour; about the T of `current` that is 2 II(T0, T) - II(T0, T0).
    gp_Pnt point;
    gp_Vec su, sv;
    current->D1(target.uv.X(), target.uv.Y(), point, su, sv);
    const gp_Vec t0 = su * du + sv * dv;
    const gp_Vec e1 = jet.u.Normalized(), e2 = jet.n.Crossed(e1);
    const gp_XYZ w = (e1 * jet.curvature(t0, e1) + e2 * jet.curvature(t0, e2)).XYZ();
    out.push_back({target.uv.X(),
                   target.uv.Y(),
                   {{2, 0, n * (du * du)},
                    {1, 1, n * (2 * du * dv)},
                    {0, 2, n * (dv * dv)},
                    {1, 0, w * (-2 * du)},
                    {0, 1, w * (-2 * dv)}},
                   5,
                   -jet.curvature(t0, t0),
                   target.share * size * size * bendWeight,
                   wished,
                   fitShare * limits.g2 * t0.SquareMagnitude()
                     * std::max({std::abs(jet.curvature(t0, t0)) / t0.SquareMagnitude(),
                                 target.easing,
                                 flatShare / size})});
  }
  return out;
}

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
                                     GraphFit graph = GraphFit::none)
{
  if (frame)
  {
    // Solve in the loop's orthonormal frame, including the unknowns and the energy's
    // right hand side. Stopping on world-coordinate components makes a rotated request
    // stop on another pass; subtracting large world positions also loses the small miss.
    const auto vector = [&](const gp_Vec& v) {
      return gp_Vec(v.Dot(frame->x), v.Dot(frame->y), v.Dot(frame->z));
    };
    // Coordinates in units of the size are rounded to a grid: what rotating them into the frame
    // rounds differently (a few ulps) then rounds to the same point, except within that much of
    // the grid's halfway points. The grid is of the size, never of the model's units, and finer
    // than a tenth of the relative tolerance.
    const double quantum    = std::min(1 / frameSteps, 0.1 * limits.g0 / size);
    const auto   coordinate = [&](double v) { return std::round(v / (size * quantum)) * quantum; };
    const auto   scaled     = [&](const gp_Vec& v) {
      const gp_Vec q = vector(v);
      return gp_Vec(coordinate(q.X()), coordinate(q.Y()), coordinate(q.Z()));
    };
    for (gp_XYZ& pole : net.poles)
      pole = scaled(gp_Vec(frame->origin, gp_Pnt(pole))).XYZ();
    for (auto* knots : {&net.knotsU, &net.knotsV})
      for (double& knot : *knots)
        knot = coordinate(knot);
    for (Target& target : targets)
    {
      if (target.perpendicular)
        target.perpendicular = vector(*target.perpendicular);
      if (target.leave)
        target.leave = vector(*target.leave);
      target.point = gp_Pnt(scaled(gp_Vec(frame->origin, target.point)).XYZ());
      target.uv.SetCoord(coordinate(target.uv.X()), coordinate(target.uv.Y()));
      // A bend, as the jets' are, in the units of the size.
      target.easing *= size;
      if (target.jet)
      {
        auto& j = *target.jet;
        j.p     = gp_Pnt(scaled(gp_Vec(frame->origin, j.p)).XYZ());
        for (gp_Vec* v : {&j.u, &j.v, &j.uu, &j.vv, &j.uv})
          *v = scaled(*v);
        j.n = vector(j.n);
        for (int c = 1; c <= 3; ++c)
          j.n.SetCoord(c, std::round(j.n.Coord(c) * frameSteps) / frameSteps);
        j.n.Normalize();
      }
    }
    std::vector<FitJump> localJumps = jumps;
    for (auto& jump : localJumps)
    {
      jump.at    = coordinate(jump.at);
      jump.along = coordinate(jump.along);
      jump.weight /= size * size;
    }
    const gp_XYZ localAlong = along ? vector(gp_Vec(*along)).XYZ() : gp_XYZ(0, 0, 0);
    std::function<double(double, double)> localDensity;
    if (density)
      localDensity = [&](double u, double v) { return density(u * size, v * size); };
    fit(net,
        fixed,
        std::move(targets),
        1,
        {limits.g0 / size, limits.g1, limits.g2},
        surface.IsNull() ? occ::handle<Geom_BSplineSurface>{} : surfaceOf(net),
        localDensity,
        along ? &localAlong : nullptr,
        localJumps,
        nullptr,
        flatness,
        leastSquares,
        graph);
    for (auto* knots : {&net.knotsU, &net.knotsV})
      for (double& knot : *knots)
        knot *= size;
    for (gp_XYZ& pole : net.poles)
      pole = frame->origin.XYZ()
             + (frame->x.XYZ() * pole.X() + frame->y.XYZ() * pole.Y() + frame->z.XYZ() * pole.Z())
                 * size;
    return surfaceOf(net);
  }
  const bool      bends   = std::any_of(targets.begin(), targets.end(), [](const Target& target) {
    return target.order > 1 || target.wish;
  });
  const bool      tangent = std::any_of(targets.begin(), targets.end(), [](const Target& target) {
    return target.order > 0;
  });
  const FitEnergy energy{(flatness < 0 ? 1 : 1 - flatness) * bendingWeight,
                         (flatness < 0 ? (bends     ? 1
                                          : tangent ? variationTangent
                                                    : variationMet)
                                       : flatness)
                           * bendingWeight * size * size};
  const bool      asked = std::any_of(targets.begin(), targets.end(), [](const Target& target) {
    return target.order > 1;
  });
  // Only free spatial fits need to iterate their nonlinear bend conditions. A graph's
  // tangents are fixed by the support normal, so its first solve already asks for G2.
  for (int pass = 0; pass < (graph != GraphFit::none ? 1
                             : asked                 ? askedPasses
                             : bends                 ? wishedPasses
                                                     : 1);
       ++pass)
  {
    // A place between those asked that the surface misses on its own is asked as well,
    // and the net fitted again from where it stood: places are added where the energy
    // alone does not keep the surface on a curve, and nowhere else, since every place
    // asked is one more thing that the surface has to follow whatever it costs it inside.
    const std::vector<gp_XYZ> start = net.poles;
    FitEnergies               energies;
    for (int round = 0;; ++round)
    {
      std::vector<int>       owner;
      const std::vector<int> missed =
        fitNet(op,
               net,
               fixed,
               conditions(targets, surface, size, limits, owner, graph),
               energy,
               density,
               along,
               jumps,
               leastSquares,
               &energies);
      if (missed.empty() || round == placeRounds)
        break;
      for (int condition : missed)
      {
        targets[owner[condition]].watch      = false;
        targets[owner[condition]].leaveAsked = true;
      }
      net.poles = start;
    }
    surface = surfaceOf(net);
  }
  return surface;
}

std::string number(double value)
{
  char text[32];
  std::snprintf(text, sizeof text, "%.3g", value);
  return text;
}

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

//! How far from a corner, as a share of a side, the face along the side is read for what it
//! asks at the corner: off the corner itself, where the faces of both sides meet and a face
//! may have no one direction.
constexpr double cornerInset = 1e-6;

//! The fractions of a side that the measure reads, in order: the places of Measure Continuity,
//! and between the first and the last of them and the side's ends places that halve their way to
//! each end down to `cornerInset`, so that what a surface does up to a corner is read as well. How
//! near is at the corner is a judgement, a fraction of the side and so the same at every size: read
//! down to 1e-5 or 1e-7 instead, the Patch tests came out the same (see
//! dox/extensions/operations.md).
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

//! How far a surface is from what a side asks, read from the surface alone.
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

//! The largest principal curvature of a surface at a place, from its jet.
double sharpest(const SurfaceJet& jet)
{
  const gp_Vec a = jet.u.Normalized(), b = jet.n.Crossed(a);
  const double aa = jet.curvature(a, a), bb = jet.curvature(b, b), ab = jet.curvature(a, b);
  return std::abs((aa + bb) / 2) + std::hypot((aa - bb) / 2, ab);
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

bool hasGaps(const Asked& asked)
{
  return std::any_of(asked.sides.begin(), asked.sides.end(), [](const Side& side) {
    return side.free;
  });
}

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

//! Work on a private, consistently directed spline rim in its intrinsic frame. Equivalent
//! pole reversal and rigid placement must be removed before sampling and net decisions,
//! rather than rounding only the eventual linear system. Original curves are measured later.
std::optional<gp_Trsf> localSplineRim(Asked& asked, bool* copyBudgetMiss = nullptr)
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

//! A failed world-space copy may instead be represented relative to an exact anchor.
//! TwoDiff certifies each subtraction; no source-pole displacement is spent here. Keep
//! the anchor as a face placement, so it never gets rounded back into the fitted poles.
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

//! A guide ending on a supported edge must leave in that face's tangent plane, and at G2
//! have its normal curvature. Test this before any fit, using the same jet as the joint.
void agreeAtGuideEnds(const Asked& asked)
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
          fail(op,
               unsupported_geometry,
               where + " but leaves its tangent plane by " + number(angle)
                 + " degrees: no surface holds the guide and is tangent to that face there");
        if (side.order < 2)
          continue;
        gp_Pnt at;
        gp_Vec first, second;
        guide.curve->D2(guide.parameter(fraction), at, first, second);
        const double mine       = second.Dot(jet.n) / first.SquareMagnitude(),
                     theirs     = jet.curvature(tangent, tangent);
        const double difference = std::abs(mine - theirs),
                     scale = std::max({std::abs(mine), std::abs(theirs), flatShare / asked.size});
        if (difference > asked.limits.g2 * scale)
          fail(op,
               unsupported_geometry,
               where + " but has a different normal curvature, by a share of "
                 + number(difference / scale) + ": G2 cannot hold there (G1 may be met)");
      }
    }
}

//! The tightest bend on the rim of what is asked: of its curves and guides, of the faces
//! that it is to run on from, and of a ball as large as the loop.
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
                    const Known&                     known = {})
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

std::vector<double> clampedKnots(double from, double to, int spans, int degree)
{
  std::vector<double> knots(degree + 1, from);
  for (int i = 1; i < spans; ++i)
    knots.push_back(from + (to - from) * i / spans);
  knots.insert(knots.end(), degree + 1, to);
  return knots;
}

//! A harmonic height polynomial for a free rim, or a quadratic fixed by interior positions.
//! Seed the exact polynomial rather than driving tiny interpolation errors into the fit.
//! The later correction and measurement hold the same limits as for a general fitted surface.
bool polynomialHeight(const Asked& asked,
                      const Frame& frame,
                      SplineNet&   net,
                      bool         approximate = false)
{
  if (asked.flatness >= 0)
    return false;
  const bool interior = !asked.guides.empty() || !asked.points.empty();
  const auto terms    = [&](const gp_Pnt& point, int degree) {
    const gp_XY         uv = frame.project(point) / asked.size;
    std::vector<double> out{1};
    gp_XY               power(1, 0);
    for (int k = 1; k <= degree; ++k)
    {
      power.SetCoord(power.X() * uv.X() - power.Y() * uv.Y(),
                     power.X() * uv.Y() + power.Y() * uv.X());
      out.push_back(power.X());
      out.push_back(power.Y());
    }
    if (interior)
      out.push_back(uv.SquareModulus());
    return out;
  };
  const auto choose = [](int n, int k) {
    double value = 1;
    for (int j = 1; j <= k; ++j)
      value *= static_cast<double>(n - j + 1) / j;
    return value;
  };
  const auto powers = [&](const std::vector<double>& knots, int degree, int i) {
    std::array<double, 6> out{1, 0, 0, 0, 0, 0};
    for (int j = 1; j <= degree; ++j)
      for (int k = j; k > 0; --k)
        out[k] += out[k - 1] * knots[i + j] / asked.size;
    for (int k = 1; k <= degree; ++k)
      out[k] /= choose(degree, k);
    return out;
  };
  for (int degree = 2; degree <= (interior ? 2 : std::min(net.degreeU, net.degreeV)); ++degree)
  {
    const int   count = 1 + 2 * degree + static_cast<int>(interior);
    math_Matrix matrix(1, count, 1, count, 0.0);
    math_Vector value(1, count, 0.0), coefficient(1, count);
    const auto  add = [&](const gp_Pnt& at, double weight) {
      const auto term = terms(at, degree);
      for (int a = 0; a < count; ++a)
      {
        value(a + 1) += weight * term[a] * frame.height(at);
        for (int b = 0; b < count; ++b)
          matrix(a + 1, b + 1) += weight * term[a] * term[b];
      }
    };
    for (const auto* curves : {&asked.sides, &asked.guides})
      for (const Side& side : *curves)
        for (int i = 0; i < places; ++i)
          add(side.meet((i + .5) / places), side.length / places);
    for (const gp_Pnt& at : asked.points)
      add(at, asked.size);
    math_Gauss solve(matrix);
    if (!solve.IsDone())
      continue;
    solve.Solve(value, coefficient);
    double     miss    = 0;
    const auto measure = [&](const gp_Pnt& at) {
      const auto term   = terms(at, degree);
      double     height = 0;
      for (int k = 0; k < count; ++k)
        height += coefficient(k + 1) * term[k];
      miss = std::max(miss, std::abs(height - frame.height(at)));
    };
    for (const auto* curves : {&asked.sides, &asked.guides})
      for (const Side& side : *curves)
        for (int i = 0; i <= places; ++i)
          measure(side.meet(static_cast<double>(i) / places));
    for (const gp_Pnt& at : asked.points)
      measure(at);
    // An interior constraint must identify the quadratic within the actual tolerance.
    // Approximate harmonic free rims are corrected afterwards, as before.
    if (miss >= (interior && !approximate ? .1 * asked.limits.g0 : 1e-3 * asked.size))
      continue;
    for (int i = 0; i < net.countU(); ++i)
      for (int j = 0; j < net.countV(); ++j)
      {
        const auto u = powers(net.knotsU, net.degreeU, i), v = powers(net.knotsV, net.degreeV, j);
        double     height = coefficient(1);
        for (int m = 1; m <= degree; ++m)
          for (int k = 0; k <= m; ++k)
          {
            const double term = choose(m, k) * u[m - k] * v[k], sign = k % 4 < 2 ? 1 : -1;
            height += coefficient(2 * m + k % 2) * sign * term;
          }
        if (interior)
          height += coefficient(count) * (u[2] + v[2]);
        net.poles[static_cast<size_t>(i) * net.countV() + j] += frame.z.XYZ() * height;
      }
    return true;
  }
  return false;
}

//! The most spans a net is given across one direction, however much is asked to be finer.
constexpr int mostSpans = 72;

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
                                const std::vector<double>& breaks)
{
  // What is wanted is summed over this many cells, which is as finely as a knot is placed.
  constexpr int       cells = 4096;
  constexpr double    snap  = 0.35;
  std::vector<double> wanted(cells, spans / (to - from));
  for (const Finer& place : finer)
  {
    const int first =
                std::clamp(static_cast<int>((place.at - place.span - from) / (to - from) * cells),
                           0,
                           cells - 1),
              last =
                std::clamp(static_cast<int>((place.at + place.span - from) / (to - from) * cells),
                           0,
                           cells - 1);
    for (int i = first; i <= last; ++i)
      wanted[i] = std::max(wanted[i], 1 / place.span);
  }
  double total = 0;
  for (double density : wanted)
    total += density * (to - from) / cells;
  const int count =
    std::clamp(static_cast<int>(std::lround(total)), spans, std::max(spans, mostSpans));
  std::vector<double> even{from};
  double              sum = 0;
  for (int i = 0; i < cells && static_cast<int>(even.size()) < count; ++i)
  {
    const double step = wanted[i] * (to - from) / cells * count / total;
    while (static_cast<int>(even.size()) < count && sum + step >= static_cast<double>(even.size()))
      even.push_back(from
                     + (to - from) * (i + (static_cast<double>(even.size()) - sum) / step) / cells);
    sum += step;
  }
  even.push_back(to);
  std::vector<double> at = breaks;
  std::sort(at.begin(), at.end());
  std::vector<std::pair<double, int>> inner;
  for (double value : at)
    if (value > from && value < to
        && (inner.empty() || value - inner.back().first > 1e-6 * (to - from) / count))
      inner.emplace_back(value, degree - 1);
  const size_t broken = inner.size();
  for (size_t i = 1; i + 1 < even.size(); ++i)
  {
    const double span  = std::min(even[i] - even[i - 1], even[i + 1] - even[i]);
    bool         taken = false;
    for (size_t k = 0; k < broken; ++k)
      taken = taken || std::abs(inner[k].first - even[i]) < snap * span;
    if (!taken)
      inner.emplace_back(even[i], 1);
  }
  std::sort(inner.begin(), inner.end());
  std::vector<double> knots(degree + 1, from);
  for (const auto& [knot, times] : inner)
    knots.insert(knots.end(), times, knot);
  knots.insert(knots.end(), degree + 1, to);
  return knots;
}

//! By what share of their curvature the faces on the two sides of a corner disagree about the
//! one surface that is to leave both tangentially (see cornerBend).
double cornerConflict(const Side& a, const Side& b)
{
  return cornerBend(a, b, false).conflict;
}

// --- A plane ---------------------------------------------------------------------------------

std::optional<TopoDS_Face> planar(const Asked& asked, const Frame& frame)
{
  const double flat = asked.limits.g0, turn = std::cos(asked.limits.g1 * M_PI / 180 / 10);
  for (const auto* curves : {&asked.sides, &asked.guides})
    for (const Side& side : *curves)
      for (int i = 0; i <= 32; ++i)
      {
        if (std::abs(frame.height(side.at(i / 32.0))) > flat)
          return {};
        if (side.order < 1)
          continue;
        const SurfaceJet jet   = side.jet(i / 32.0);
        const gp_Vec     along = side.tangent(i / 32.0).Normalized(), across = jet.n.Crossed(along);
        if (std::abs(jet.n.Dot(frame.z)) < turn
            || (side.order > 1
                && std::abs(jet.curvature(across, across)) * asked.size
                     > asked.limits.g2 * flatShare))
          return {};
      }
  for (const gp_Pnt& point : asked.points)
    if (std::abs(frame.height(point)) > flat)
      return {};
  std::vector<gp_XY> loop;
  for (const gp_Pnt& point : polygon(asked.sides))
    loop.push_back(frame.project(point));
  enclosed(loop);
  // In the plane, but beside the loop, is not on the patch: as for the other constructions.
  for (const Side& guide : asked.guides)
    for (int i = 0; i <= 16; ++i)
      if (!inside(loop, frame.project(guide.at(i / 16.0)), outsideSlack * asked.size))
        fail(op, invalid_input, "a guide leaves the boundary");
  for (const gp_Pnt& point : asked.points)
    if (!inside(loop, frame.project(point), outsideSlack * asked.size))
      fail(op, invalid_input, "a point lies outside the boundary");
  BRepBuilderAPI_MakeWire wire;
  for (const Side& side : asked.sides)
    wire.Add(BRepBuilderAPI_MakeEdge(side.curve,
                                     std::min(side.from, side.to),
                                     std::max(side.from, side.to)));
  if (!wire.IsDone())
    return {};
  BRepBuilderAPI_MakeFace face(gp_Pln(frame.origin, gp_Dir(frame.z)), wire.Wire(), true);
  if (!face.IsDone())
    return {};
  return face.Face();
}

// --- Four sides ------------------------------------------------------------------------------

//! What share of the tolerance a side that has no polynomial form (an arc, a rational spline)
//! may be off the spline that stands for it in a four-sided surface, and the most pieces of
//! degree 5 that spline is given: a half circle is within 1e-8 of its radius with 16.
constexpr double plainShare    = 0.1;
constexpr int    plainSegments = 64;

//! A side as a polynomial B-spline over [0, length] with the side's own run of parameters, or
//! nothing where no such curve comes within the tolerance.
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

std::vector<double> flatKnots(const occ::handle<Geom_BSplineCurve>& curve)
{
  const NCollection_Array1<double>& sequence = curve->KnotSequence();
  return {sequence.begin(), sequence.end()};
}

//! How near each other, as a share of the loop's size, the surfaces of two nets in a row of the
//! ladder are to count as having come as near the fairest surface
//! as finer nets would take them. A share of the size, so that the stop is the same at every size
//! and pose; the length tolerance, which is absolute, would not be, and nothing inside the loop is
//! asked to any tolerance. The value is a judgement: between the nets of the tests, these distances
//! spread from 1e-5 to 1 of the size with no gap between settled and unsettled pairs; a seat, two
//! mouse tops and a hull bottom that a coarse net filled 0.6 to 9 off their masters stop on
//! the same nets for any value from 3.1e-4 to 2.4e-3 (see dox/extensions/operations.md).
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

//! A finer net is not tried when the one before it brought the worst miss down by less than
//! to this share of what it was: the miss has then stopped falling, and what is left of the
//! ladder would be spent on a net that comes no nearer.
constexpr double stillFalling = 0.75;

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

//! How far one side turns into the next at their corner, in degrees: none where it goes straight
//! on, 180 where it turns back.
double turning(const Side& from, const Side& to)
{
  return std::acos(
           std::clamp(from.tangent(1).Normalized().Dot(to.tangent(0).Normalized()), -1.0, 1.0))
         * 180 / M_PI;
}

//! Whether two sides meet at a corner of a quad: one that turns by more than the angular
//! tolerance, and turns back by less: four curves make a quad at turns from under 4 to over 176
//! degrees.
bool cornered(const Side& from, const Side& to, double g1Degrees)
{
  const double angle = turning(from, to);
  return angle > g1Degrees && angle < 180 - g1Degrees;
}

//! Whether four sides meet at corners, so that they can be the edges of one surface.
bool squarish(const std::vector<Side>& sides, double g1Degrees)
{
  if (sides.size() != 4)
    return false;
  for (size_t i = 0; i < 4; ++i)
    if (!cornered(sides[i], sides[(i + 1) % 4], g1Degrees))
      return false;
  return true;
}

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

//! Whether two edges are certainly farther apart than `gap`, from the pole boxes of their
//! Bezier pieces: positive weights keep each piece inside its poles. Pieces whose boxes come
//! nearer than `gap` are halved, the longer first, a few times over; false says only that it
//! could not be shown so, and the edges are then measured.
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

std::optional<Built> foursided(const Asked& asked, const Frame& frame)
{
  PatchPhase phase("quad");
  const auto found = quadRuns(asked.sides, asked.limits.g1);
  if (!found)
    return {};
  std::vector<Run> runs = *found;
  // A side made only of gaps of an open boundary: nothing is asked along it.
  bool gap[4];
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
  int order = 0;
  for (const Side& side : asked.sides)
    order = std::max(order, side.order);
  const bool fitting = order > 0 || asked.flatness >= 0 || asked.controls.flow > 0
                       || !asked.guides.empty() || !asked.points.empty();
  const double lu = (runs[0].length + runs[2].length) / 2,
               lv = (runs[1].length + runs[3].length) / 2;
  occ::handle<Geom_BSplineCurve> curves[4];
  for (int i = 0; i < 4; ++i)
  {
    // Round the loop the top and the left run backwards.
    if (runs[i].parts.size() == 1)
    {
      curves[i] = plain(*runs[i].parts[0], i % 2 ? lv : lu, i >= 2, asked.limits.g0);
      if (curves[i].IsNull())
        return {};
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
        return {};
      if (!joined)
      {
        joined.emplace(piece);
        continue;
      }
      junctions.push_back(joined->BSplineCurve()->LastParameter());
      if (!joined->Add(piece, asked.limits.g0, true, true))
        return {};
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
  for (int i = 0; i < 2; ++i)
  {
    auto &a = curves[i], &b = curves[i + 2];
    // Asked for a quad, a side above the net's highest degree, or one with more poles at the
    // net's degree than a fitted net may have, is brought down to it as an arc is, within the
    // same share of the tolerance and over the same parameters.
    const int least = !fitting ? 1 : order < 2 ? 3 : 5;
    for (auto* curve : {&a, &b})
      if (asked.controls.quad == 1
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
          return {};
        *curve = lower.Curve();
      }
    const int degree = std::max({a->Degree(), b->Degree(), least});
    if (degree > 5)
      return {};
    a->IncreaseDegree(degree);
    b->IncreaseDegree(degree);
    const double same = BRepExt::tolerance::parameterMerge * (i ? lv : lu);
    a->InsertKnots(b->Knots(), b->Multiplicities(), same, false);
    b->InsertKnots(a->Knots(), a->Multiplicities(), same, false);
    if (a->NbPoles() != b->NbPoles())
      return {};
  }
  // Quality budgets are minimum span counts.
  const int minimum = asked.controls.quality == 2 ? 16 : asked.controls.quality == 1 ? 10 : 1;
  for (int i = 0; i < 2; ++i)
    for (int k = 1; k < minimum; ++k)
      for (int c : {i, i + 2})
        curves[c]->InsertKnot((i ? lv : lu) * k / minimum);
  // A place of a side is on the surface's own edge, as far along it.
  const auto along = [&](size_t s, double f) {
    return s == 0   ? gp_XY(f * lu, 0)
           : s == 1 ? gp_XY(lu, f * lv)
           : s == 2 ? gp_XY((1 - f) * lu, lv)
                    : gp_XY(0, (1 - f) * lv);
  };
  // Where an edge is in its side, for a place of it to be found on the surface.
  const auto place = [&](const Side& side, double f) {
    for (size_t s = 0; s < 4; ++s)
      for (size_t k = 0; k < runs[s].parts.size(); ++k)
        if (runs[s].parts[k] == &side)
        {
          const double end = k + 1 < runs[s].parts.size() ? runs[s].starts[k + 1] : 1;
          return along(s, runs[s].starts[k] + f * (end - runs[s].starts[k]));
        }
    return gp_XY();
  };
  // The way the isocurves across a side are to leave it, as the one direction in the
  // tangent plane that the derivative across is to have nothing of: square to the side for
  // Normal. For Adjacent they turn from the side before to the side after, which are those
  // isocurves at the corners; for Next they head for the same share of the opposite side,
  // turned at the corners into those same sides. Square to the side on a rectangle, both.
  const auto axis = [&](size_t s, double f) -> gp_Vec {
    const auto [part, g] = runs[s].at(f);
    if (asked.controls.flow == 1)
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
    if (asked.controls.flow == 3 && fromGap && toGap)
      return part->tangent(g);
    const auto square = [&](const gp_Vec& way, double t) {
      const gp_Vec along = runs[s].tangent(t).Normalized(), out = way - along * way.Dot(along);
      return out.Magnitude() > 1e-12 ? out.Normalized() : out;
    };
    if (fromGap)
      before = asked.controls.flow == 2 ? chord(0) : square(after, 0);
    if (toGap)
      after = asked.controls.flow == 2 ? chord(1) : square(before, 1);
    gp_Vec way = before * (1 - f) + after * f;
    if (asked.controls.flow == 2)
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
  };
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
      std::vector<char> fixed(net.poles.size(), 0);
      for (int i = 0; i < nu; ++i)
        for (int j = 0; j < nv; ++j)
          fixed[static_cast<size_t>(i) * nv + j] =
            ((i == 0 || i == nu - 1) && (j == 0 || j == nv - 1)) || (j == 0 && !gap[0])
            || (i == nu - 1 && !gap[1]) || (j == nv - 1 && !gap[2]) || (i == 0 && !gap[3]);
      std::vector<Target> targets;
      for (size_t s = 0; s < 4; ++s)
      {
        const int count = placesPerPole * (s % 2 ? nv : nu);
        for (int k = 0; k < count; ++k)
        {
          const double f       = (k + 0.5) / count;
          const auto [part, g] = runs[s].at(f);
          if (part->free || (part->order < 1 && !asked.controls.flow))
            continue;
          Target target{part->at(g),
                        along(s, f),
                        part->order,
                        part->order ? std::optional<SurfaceJet>(part->jet(g)) : std::nullopt,
                        s % 2 ? gp_XY(1, 0) : gp_XY(0, 1),
                        runs[s].length / count / asked.size};
          if (asked.controls.flow)
            target.perpendicular = axis(s, f);
          if (part->order)
            target.easing = part->easing(g);
          targets.push_back(target);
        }
      }
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

      // A guide from side to side runs straight across the rectangle, through the place of
      // each guide it crosses: the nearest point of a surface that does not yet rise to it
      // would drag the parameters along with the rise. Others go by that nearest point.
      struct Across
      {
        bool                                  laid = false;
        std::vector<std::pair<double, gp_XY>> stops;
      };

      std::vector<Across> across(asked.guides.size());
      const auto          onSide = [&](const gp_XY& uv) {
        const double slack = asked.limits.g0;
        return std::abs(uv.X()) <= slack || std::abs(uv.X() - lu) <= slack
               || std::abs(uv.Y()) <= slack || std::abs(uv.Y() - lv) <= slack;
      };
      for (size_t g = 0; g < asked.guides.size(); ++g)
      {
        const gp_XY
          a = projection.ValueOfUV(asked.guides[g].at(0), BRepExt::tolerance::exactLength).XY(),
          b = projection.ValueOfUV(asked.guides[g].at(1), BRepExt::tolerance::exactLength).XY();
        if (onSide(a) && onSide(b))
          across[g] = {true, {{0.0, a}, {1.0, b}}};
      }
      constexpr int crossingPlaces = 129;
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
          for (double step = 1.0 / crossingPlaces; step > 1e-12; step /= 8)
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
      const auto laidAt = [](std::vector<std::pair<double, gp_XY>> stops, double f) {
        std::sort(stops.begin(), stops.end(), [](const auto& a, const auto& b) {
          return a.first < b.first;
        });
        size_t i = 0;
        while (i + 2 < stops.size() && f > stops[i + 1].first)
          ++i;
        const double t =
          (f - stops[i].first) / std::max(stops[i + 1].first - stops[i].first, 1e-12);
        return gp_XY(stops[i].second + (stops[i + 1].second - stops[i].second) * t);
      };
      // A guide that ends on one so laid, at a side or at another such guide, is laid to there.
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
            for (double step = 1.0 / crossingPlaces; step > 1e-12; step /= 8)
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
      for (size_t g = 0; g < asked.guides.size(); ++g)
      {
        const Side& guide = asked.guides[g];
        const int   count = std::max(129, placesPerPole * std::max(nu, nv));
        for (int k = 0; k <= count; ++k)
        {
          const double f     = static_cast<double>(k) / count,
                       share = guide.length / count / asked.size;
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
    for (size_t side = 0; asked.controls.flow && side < 4; ++side)
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
      int order = 2;
      for (const Side* part : runs[s].parts)
        order = std::min(order, part->order);
      return order;
    };
    if (gap[1] && gap[3] && !gap[0] && !gap[2])
      built.crossing = Crossing{true, least(0), least(2)};
    else if (gap[0] && gap[2] && !gap[1] && !gap[3])
      built.crossing = Crossing{false, least(3), least(1)};
    const double before = best ? best->measured.excess : HUGE_VAL;
    if (built.measured.excess < before)
      best = built;
    // A finer net that comes no nearer is the last, as for a trimmed surface.
    if (built.measured.within || !fitting || built.measured.excess > stillFalling * before)
      break;
    // Twice the spans for the next round.
    for (int i = 0; i < 2; ++i)
    {
      const std::vector<double> knots = flatKnots(curves[i]);
      for (size_t k = 0; k + 1 < knots.size(); ++k)
        if (knots[k + 1] > knots[k])
          for (int c : {i, i + 2})
            curves[c]->InsertKnot((knots[k] + knots[k + 1]) / 2);
    }
  }
  return best;
}

// --- Any loop --------------------------------------------------------------------------------

//! An edge of a trimmed face: the curve of a side in the surface's parameters, run by the
//! side's own parameter, and the furthest the surface is there from the side at that parameter.
using patch_side::Rim;

//! What share of the tolerance the curve of an edge in the surface's parameters may add to
//! the distance of the surface from the edge, which no such curve can do anything about.
constexpr double preimageShare = 0.05;

//! The most poles such a curve is given before it is taken as it is.
constexpr int preimagePoles = 1000;

//! The most passes in which its poles are moved onto their places: each takes a third or more
//! off what is left, so that these leave a millionth of it.
constexpr int preimagePasses = 40;

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
             const std::function<gp_XY(double)>&        parameterLayout = {})
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

using patch_side::closeGaps;

//! How much more than the deviation that was measured an edge is given as its tolerance: the
//! deviation is read at places, and between them it may be a little more.
constexpr double edgeHeadroom = 1.1;

//! The face of `surface` inside the loop. Each edge is the side's own curve, with `rims` as
//! its curve in the surface's parameters.
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

//! Check the p-curves that sewing actually delivers, including its same-parameter repair.
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

//! Trim an exact continuation and measure its delivered parameter curves, rather than
//! merely measuring the untrimmed surface at nearest points.
std::optional<Built> continued(const occ::handle<Geom_BSplineSurface>&    surface,
                               const Asked&                               asked,
                               const Frame&                               frame,
                               const std::function<gp_XY(const gp_Pnt&)>& layout)
{
  std::vector<gp_XY> loop;
  for (const gp_Pnt& point : polygon(asked.sides))
    loop.push_back(layout(point));
  // The surface is exact: its edges' curves are made as near as for an exchange.
  std::vector<Rim> rims;
  for (const Side& side : asked.sides)
    rims.push_back(preimage(surface,
                            side,
                            frame,
                            std::min(asked.limits.g0, exchangeWithin),
                            asked.size,
                            layout));
  Built built{
    3,
    {},
    measureAll(surface,
               asked,
               [&](const Side& side, double f, const gp_Pnt&) -> std::optional<Place> {
                 const auto i = asked.boundaryIndex(side);
                 return i ? std::optional<Place>(
                              Place{rims[*i].onSurface->Value(side.parameter(f)).XY(), true})
                          : std::nullopt;
               })};
  if (!built.measured.within)
    return std::nullopt;
  built.face = trimmed(surface, asked, rims, enclosed(loop) < 0);
  if (!BRepCheck_Analyzer(built.face).IsValid())
    return std::nullopt;
  return built;
}

//! A G0 loop may have four logical sides split into tangent edge chains. Seed a
//! Coons surface from those chains, including degrees above the fitter's limit.
//! `plain` bounds rational approximation and Coons averages mismatched corner ends;
//! every original edge still supplies its own rim and measured constraint. This path
//! never calls fitNet and only follows a regular fit that has missed its limits.
std::optional<Built> chainedQuad(const Asked& asked, const Frame& frame)
{
  if (asked.sides.size() <= 4 || !asked.guides.empty() || !asked.points.empty()
      || asked.flatness >= 0 || asked.controls.flow || asked.controls.quad >= 0
      || asked.controls.quality
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

//! If every edge supplies the very same spline surface, its continuation is a candidate
//! before a variational fill. This uses the supplied geometry, never a guessed interior
//! reconstructed from G0 curves. Different supports and explicit shape controls still fit.
std::optional<Built> splineContinuation(const Asked& asked, const Frame&)
{
  if (!asked.guides.empty() || !asked.points.empty())
    return std::nullopt;
  const auto surface = occ::down_cast<Geom_BSplineSurface>(asked.sides.front().held);
  if (surface.IsNull())
    return std::nullopt;
  for (const Side& side : asked.sides)
  {
    const auto other = occ::down_cast<Geom_BSplineSurface>(side.held);
    if (other.IsNull() || other->UDegree() != surface->UDegree()
        || other->VDegree() != surface->VDegree() || other->NbUPoles() != surface->NbUPoles()
        || other->NbVPoles() != surface->NbVPoles() || other->NbUKnots() != surface->NbUKnots()
        || other->NbVKnots() != surface->NbVKnots()
        || other->IsUPeriodic() != surface->IsUPeriodic()
        || other->IsVPeriodic() != surface->IsVPeriodic())
      return std::nullopt;
    if (other == surface)
      continue;
    for (int i = 1; i <= surface->NbUKnots(); ++i)
      if (surface->UKnot(i) != other->UKnot(i)
          || surface->UMultiplicity(i) != other->UMultiplicity(i))
        return std::nullopt;
    for (int j = 1; j <= surface->NbVKnots(); ++j)
      if (surface->VKnot(j) != other->VKnot(j)
          || surface->VMultiplicity(j) != other->VMultiplicity(j))
        return std::nullopt;
    for (int i = 1; i <= surface->NbUPoles(); ++i)
      for (int j = 1; j <= surface->NbVPoles(); ++j)
        if (surface->Pole(i, j).Distance(other->Pole(i, j)) > 1e-12 * asked.size
            || surface->Weight(i, j) != other->Weight(i, j))
          return std::nullopt;
  }
  // The supplied p-curves already use these parameters. Re-projecting and fitting them
  // adds error and work; private exact conversions retain the original edge parameter.
  std::vector<Rim>   rims;
  std::vector<gp_XY> loop;
  for (const Side& side : asked.sides)
  {
    const double first = std::min(side.from, side.to), last = std::max(side.from, side.to);
    const auto   cut = new Geom2d_TrimmedCurve(side.onHeld, first, last);
    rims.push_back({Geom2dConvert::CurveToBSplineCurve(cut), 0});
    for (int k = 0; k < 32; ++k)
      loop.push_back(side.onHeld->Value(side.parameter(k / 32.0)).XY());
  }
  Built built{3, {}, measureAll(surface, asked, [&](const Side& side, double f, const gp_Pnt&) {
                return std::optional<Place>(
                  Place{side.onHeld->Value(side.parameter(f)).XY(), true});
              })};
  if (!built.measured.within)
    return std::nullopt;
  built.face = trimmed(surface, asked, rims, enclosed(loop) < 0);
  if (!BRepCheck_Analyzer(built.face).IsValid())
    return std::nullopt;
  return built;
}

//! A toroidal neighbour has an exact rational continuation too. Trim the analytic torus
//! in its angular parameters before conversion, so the rational rectangle has no seam
//! inside the requested window. Its B-spline parameters differ from the angles.
std::optional<Built> toroidal(const Asked& asked, const Frame& frame)
{
  if (!asked.guides.empty() || !asked.points.empty())
    return std::nullopt;
  const Side& first = asked.sides.front();
  if (first.support.IsNull() || BRepAdaptor_Surface(first.support).GetType() != GeomAbs_Torus)
    return std::nullopt;
  const auto torus = BRepAdaptor_Surface(first.support).Torus();
  for (const Side& side : asked.sides)
  {
    if (side.support.IsNull() || BRepAdaptor_Surface(side.support).GetType() != GeomAbs_Torus)
      return std::nullopt;
    const auto here = BRepAdaptor_Surface(side.support).Torus();
    if (here.Location().Distance(torus.Location()) > 1e-12 * asked.size
        || std::abs(here.MajorRadius() - torus.MajorRadius()) > 1e-12 * asked.size
        || std::abs(here.MinorRadius() - torus.MinorRadius()) > 1e-12 * asked.size
        || !here.Axis().IsParallel(torus.Axis(), BRepExt::tolerance::fineAngle))
      return std::nullopt;
  }
  ShapeAnalysis_Surface angular(first.held);
  const gp_XY           anchor = angular.ValueOfUV(first.meet(.5), 1e-10).XY();
  gp_XY                 low = anchor, high = anchor;
  for (const Side& side : asked.sides)
    for (int k = 0; k <= places; ++k)
    {
      const gp_XY uv = angular.ValueOfUV(side.meet(static_cast<double>(k) / places), 1e-10).XY();
      const gp_XY at(anchor.X() + std::remainder(uv.X() - anchor.X(), 2 * M_PI),
                     anchor.Y() + std::remainder(uv.Y() - anchor.Y(), 2 * M_PI));
      low.SetCoord(std::min(low.X(), at.X()), std::min(low.Y(), at.Y()));
      high.SetCoord(std::max(high.X(), at.X()), std::max(high.Y(), at.Y()));
    }
  const gp_XY margin = (high - low) * .1;
  low -= margin;
  high += margin;
  if (high.X() - low.X() >= 2 * M_PI || high.Y() - low.Y() >= 2 * M_PI)
    return std::nullopt;
  const auto surface = GeomConvert::SurfaceToBSplineSurface(
    new Geom_RectangularTrimmedSurface(first.held, low.X(), high.X(), low.Y(), high.Y()));
  ShapeAnalysis_Surface projection(surface);
  const auto layout = [&](const gp_Pnt& point) { return projection.ValueOfUV(point, 1e-10).XY(); };
  return continued(surface, asked, frame, layout);
}

//! Continue a spherical neighbour exactly, using stereographic parameters. Unlike longitude
//! and latitude this rational biquadratic rectangle has no pole or seam inside the cap.
//! Only a rim wholly on the same ball qualifies; the delivered joint is still measured.
std::optional<Built> spherical(const Asked& asked, const Frame& flat)
{
  if (!asked.guides.empty() || !asked.points.empty())
    return std::nullopt;
  std::optional<gp_Sphere> ball;
  std::vector<gp_Cylinder> cylinders;
  for (const Side& side : asked.sides)
  {
    if (side.support.IsNull())
      return std::nullopt;
    BRepAdaptor_Surface support(side.support);
    if (support.GetType() == GeomAbs_Cylinder)
    {
      cylinders.push_back(support.Cylinder());
      continue;
    }
    if (support.GetType() != GeomAbs_Sphere)
      return std::nullopt;
    const gp_Sphere here = support.Sphere();
    if (ball
        && (ball->Location().Distance(here.Location()) > 1e-12 * asked.size
            || std::abs(ball->Radius() - here.Radius()) > 1e-12 * asked.size))
      return std::nullopt;
    ball = here;
  }
  if (!ball)
  {
    // Equal rounds about intersecting axes admit the same ball between them. Recover
    // its centre from the axes, then verify the complete rim against that ball.
    if (cylinders.size() < 2)
      return std::nullopt;
    math_Matrix  matrix(1, 3, 1, 3, 0.0);
    math_Vector  value(1, 3, 0.0), centre(1, 3);
    const gp_Pnt origin = cylinders.front().Location();
    for (const auto& cylinder : cylinders)
    {
      const gp_XYZ d = cylinder.Axis().Direction().XYZ(),
                   q = cylinder.Location().XYZ() - origin.XYZ();
      for (int i = 1; i <= 3; ++i)
        for (int j = 1; j <= 3; ++j)
        {
          const double m = (i == j ? 1 : 0) - d.Coord(i) * d.Coord(j);
          matrix(i, j) += m;
          value(i) += m * q.Coord(j);
        }
    }
    math_Gauss solve(matrix);
    if (!solve.IsDone())
      return std::nullopt;
    solve.Solve(value, centre);
    ball = gp_Sphere(
      gp_Ax3(gp_Pnt(origin.XYZ() + gp_XYZ(centre(1), centre(2), centre(3))), gp_Dir(0, 0, 1)),
      cylinders.front().Radius());
  }
  for (const auto& cylinder : cylinders)
    if (std::abs(cylinder.Radius() - ball->Radius()) > 1e-12 * asked.size
        || gp_Vec(cylinder.Location(), ball->Location())
               .Crossed(gp_Vec(cylinder.Axis().Direction()))
               .Magnitude()
             > 1e-12 * asked.size)
      return std::nullopt;
  const double radius = ball->Radius();
  gp_Vec       z(ball->Location(), flat.origin);
  if (z.Magnitude() < 1e-9 * radius)
    return std::nullopt;
  z.Normalize();
  const gp_Vec x = (flat.x - z * flat.x.Dot(z)).Normalized(), y = z.Crossed(x);
  const Frame  frame{ball->Location(), x, y, z};
  const auto   layout = [&](const gp_Pnt& point) {
    const gp_Vec to(frame.origin, point);
    return gp_XY(to.Dot(x), to.Dot(y)) * (2 * radius / (radius + to.Dot(z)));
  };
  std::vector<gp_XY> loop;
  for (const gp_Pnt& point : polygon(asked.sides))
  {
    if (radius + frame.height(point) < 0.1 * radius)
      return std::nullopt;
    loop.push_back(layout(point));
  }
  gp_XY low = loop.front(), high = low;
  for (const gp_XY& at : loop)
  {
    low.SetCoord(std::min(low.X(), at.X()), std::min(low.Y(), at.Y()));
    high.SetCoord(std::max(high.X(), at.X()), std::max(high.Y(), at.Y()));
  }
  const gp_XY margin = (high - low) * 0.1;
  low -= margin;
  high += margin;
  // Two quadratic spans each way keep every homogeneous weight positive even on a
  // deep cap. The knot at zero does not interrupt the geometry of the exact sphere.
  const double u[]  = {low.X(), low.X() / 2, 0, high.X() / 2, high.X()},
               v[]  = {low.Y(), low.Y() / 2, 0, high.Y() / 2, high.Y()},
               uu[] = {low.X() * low.X(), 0, 0, 0, high.X() * high.X()},
               vv[] = {low.Y() * low.Y(), 0, 0, 0, high.Y() * high.Y()};
  if (low.X() >= 0 || high.X() <= 0 || low.Y() >= 0 || high.Y() <= 0)
    return std::nullopt;
  NCollection_Array2<gp_Pnt> poles(1, 5, 1, 5);
  NCollection_Array2<double> weights(1, 5, 1, 5);
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j)
    {
      const double q = (uu[i] + vv[j]) / (4 * radius * radius), w = 1 + q;
      weights.SetValue(i + 1, j + 1, w);
      poles.SetValue(
        i + 1,
        j + 1,
        gp_Pnt(frame.origin.XYZ()
               + (x.XYZ() * u[i] + y.XYZ() * v[j] + z.XYZ() * (radius * (1 - q))) / w));
    }
  NCollection_Array1<double> us(1, 3), vs(1, 3);
  us.SetValue(1, low.X());
  us.SetValue(2, 0);
  us.SetValue(3, high.X());
  vs.SetValue(1, low.Y());
  vs.SetValue(2, 0);
  vs.SetValue(3, high.Y());
  NCollection_Array1<int> mult(1, 3);
  mult.SetValue(1, 3);
  mult.SetValue(2, 2);
  mult.SetValue(3, 3);
  const occ::handle<Geom_BSplineSurface> surface =
    new Geom_BSplineSurface(poles, weights, us, vs, mult, mult, 2, 2);
  return continued(surface, asked, frame, layout);
}

//! Continue two tangent planes and the cylinder between them as one ruled NURBS.
//! Its only curvature breaks are the two existing plane/round joins, along the round's axis.
std::optional<Built> ruledRound(const Asked& asked, const Frame& frame)
{
  if (!asked.guides.empty() || !asked.points.empty())
    return std::nullopt;
  std::optional<gp_Cylinder> round;
  std::vector<gp_Pln>        planes;
  for (const Side& side : asked.sides)
  {
    if (side.held.IsNull())
      return std::nullopt;
    // G0 still comes with the same geometric support, held separately from tangency.
    const auto cylinder = occ::down_cast<Geom_CylindricalSurface>(side.held);
    const auto plane    = occ::down_cast<Geom_Plane>(side.held);
    if (!cylinder.IsNull())
    {
      const gp_Cylinder here = cylinder->Cylinder();
      if (round
          && (std::abs(here.Radius() - round->Radius()) > 1e-12 * asked.size
              || !here.Axis().IsCoaxial(round->Axis(),
                                        BRepExt::tolerance::fineAngle,
                                        1e-12 * asked.size)))
        return std::nullopt;
      round = here;
    }
    else if (!plane.IsNull())
    {
      const gp_Pln here = plane->Pln();
      if (std::none_of(planes.begin(), planes.end(), [&](const gp_Pln& p) {
            return p.Axis().IsParallel(here.Axis(), BRepExt::tolerance::fineAngle)
                   && p.Distance(here.Location()) < 1e-12 * asked.size;
          }))
        planes.push_back(here);
    }
    else
      return std::nullopt;
  }
  if (!round || planes.size() != 2)
    return std::nullopt;
  const gp_Pnt centre = round->Location();
  const gp_Vec axis(round->Axis().Direction());
  const double radius = round->Radius();
  gp_Vec       radial[2];
  for (int i = 0; i < 2; ++i)
  {
    const gp_Vec n(planes[i].Axis().Direction());
    const double distance = gp_Vec(centre, planes[i].Location()).Dot(n);
    if (std::abs(n.Dot(axis)) > BRepExt::tolerance::fineAngle
        || std::abs(std::abs(distance) - radius) > 1e-12 * asked.size)
      return std::nullopt;
    radial[i] = n * distance;
  }
  const double cosine = radial[0].Dot(radial[1]) / (radius * radius);
  if (cosine < -0.9 || cosine > 0.9)
    return std::nullopt;
  const gp_Vec t0                   = (radial[1] - radial[0] * cosine).Normalized(),
               t1                   = (radial[1] * cosine - radial[0]).Normalized();
  const gp_XYZ               c      = centre.XYZ();
  const double               length = asked.size, weight = std::sqrt((1 + cosine) / 2);
  const gp_XYZ               profile[] = {c + radial[0].XYZ() - t0.XYZ() * length,
                                          c + radial[0].XYZ() - t0.XYZ() * (length / 2),
                                          c + radial[0].XYZ(),
                                          c + (radial[0].XYZ() + radial[1].XYZ()) / (1 + cosine),
                                          c + radial[1].XYZ(),
                                          c + radial[1].XYZ() + t1.XYZ() * (length / 2),
                                          c + radial[1].XYZ() + t1.XYZ() * length};
  NCollection_Array1<gp_Pnt> cp(1, 7);
  NCollection_Array1<double> cw(1, 7);
  for (int i = 0; i < 7; ++i)
  {
    cp.SetValue(i + 1, gp_Pnt(profile[i]));
    cw.SetValue(i + 1, i == 3 ? weight : 1);
  }
  const double               span = length / (2 * radius * std::sqrt((1 - cosine) / 2));
  NCollection_Array1<double> us(1, 4);
  us.SetValue(1, 0);
  us.SetValue(2, span);
  us.SetValue(3, span + 1);
  us.SetValue(4, 2 * span + 1);
  NCollection_Array1<int> um(1, 4);
  um.SetValue(1, 3);
  um.SetValue(2, 2);
  um.SetValue(3, 2);
  um.SetValue(4, 3);
  const occ::handle<Geom_BSplineCurve> curve = new Geom_BSplineCurve(cp, cw, us, um, 2);
  double                               low = HUGE_VAL, high = -HUGE_VAL;
  for (const auto& point : polygon(asked.sides))
  {
    const double v = gp_Vec(centre, point).Dot(axis);
    low            = std::min(low, v);
    high           = std::max(high, v);
  }
  low -= .1 * asked.size;
  high += .1 * asked.size;
  NCollection_Array2<gp_Pnt> poles(1, 7, 1, 2);
  NCollection_Array2<double> weights(1, 7, 1, 2);
  for (int i = 0; i < 7; ++i)
    for (int j = 0; j < 2; ++j)
    {
      poles.SetValue(i + 1, j + 1, gp_Pnt(profile[i] + axis.XYZ() * (j ? high : low)));
      weights.SetValue(i + 1, j + 1, cw(i + 1));
    }
  NCollection_Array1<double> vs(1, 2);
  vs.SetValue(1, low);
  vs.SetValue(2, high);
  NCollection_Array1<int> vm(1, 2);
  vm.Init(2);
  const occ::handle<Geom_BSplineSurface> surface =
    new Geom_BSplineSurface(poles, weights, us, vs, um, vm, 2, 1);
  const auto layout = [&](const gp_Pnt& point) {
    const double                v = gp_Vec(centre, point).Dot(axis);
    GeomAPI_ProjectPointOnCurve project(point.Translated(-axis * v), curve);
    return gp_XY(project.LowerDistanceParameter(), v);
  };
  return continued(surface, asked, frame, layout);
}

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

//! How much less the surface moved onto its curves must bend than the one fitted at the
//! tolerance to be kept instead.
constexpr double clearlyFairer = 0.75;

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

//! How far the rim may turn within one span of the coarsest net, in radians, before the net
//! is made finer there: its own direction, or the tangent plane of the face it runs on.
constexpr double turnPerSpan = 0.2;

//! The most poles of a net that is free in space. A net that is to go over smoothly into a
//! face does better so than as a height over the mean plane: its parameters can then run with
//! the surface where it is steep, and the energy, which is of the parameters, is nearer to
//! one of the shape (a ball's cap cut at 6 of 10 comes out 0.13 off the ball, as a height 1.3).
//! It has three times the unknowns and takes 27 times the work: 25 x 25 poles free are about
//! what a height of 57 x 57 takes.
constexpr int freePoles = 640;

//! The places along a side at which the turning of its face is read.
constexpr int turningSteps = 48;

//! Whether the neighbours' bend is wished for also where two of them disagree about it at a
//! corner of the rim.
constexpr bool wishAcrossCorners = false;

//! The nets that are tried with the neighbours' bend wished for, before the surface is made
//! tangent and no more.
constexpr int wishNets = 3;

//! The steps along a curve at which it is looked for the lines of the net that it crosses.
constexpr int crossingSteps = 256;

//! Two lines of the net across which it may change its bend count as well apart when this many
//! of them, as far apart, would cross the loop: about the spans of the third net.
constexpr double apartSpans = 16;

//! What a surface that meets its edges and is not one to keep counts as missing by, in units
//! of the limits: no more than any that misses an edge, so that where no net gives a patch it
//! is said that the edges can be met and what is wrong with the surface that meets them.
constexpr double unkeptExcess = 1;

//! The nets that may meet the edges with a surface that is not one to keep before it is given
//! up: the two coarsest and one more.
constexpr int mostUnkept = 3;

//! A miss of no more than this many times its limit is nearly met.
constexpr double nearly = 2;

//! What a step of the second derivative across a line of lower continuity counts for against
//! the bend, per length of the line and with the loop's size as the unit; the places to a span
//! along the line at which it is read; and how far from its corner, as a share of the loop's
//! size, a step counts in full.
constexpr double jumpWeight = 30, jumpReach = 0.3;
constexpr int    jumpPlaces = 3;

//! Within this share of the miss before, a miss counts as the same: what is in the way is then
//! not the net.
constexpr double same = 0.02;

//! The coarsest nets, whose misses say nothing yet of how the finer ones will do.
constexpr int unruly = 2;

// --- Over a chart of the loop's own ----------------------------------------------------------

//! The places along a side at which its face's normal is read for whether any plane is one the
//! fill could be a height over (see graphMargin), its ends included.
constexpr int gaussPlaces = 64;

//! The places along a guide that the loop's chart is drawn through (see proxyChart).
constexpr int guidePins = 32;

//! The places the loop is spanned by round it for its chart (see proxyChart): the disc inside
//! has about a tenth of the square of that many nodes.
constexpr int chartPlaces = 160;

//! Which way along N x T, N the face's normal and T the side's tangent, a fill leaves a side it
//! is to be tangent to: away from the face, which lies to one side of the edge in the face's own
//! parameters. Zero where nothing is asked or that side is not found.
double leaving(const Side& side)
{
  if (side.free || side.order < 1 || side.support.IsNull())
    return 0;
  double     first, last;
  const auto onFace = BRep_Tool::CurveOnSurface(side.edge, side.support, first, last);
  if (onFace.IsNull())
    return 0;
  gp_Pnt2d uv;
  gp_Vec2d along;
  onFace->D1(side.parameter(0.5), uv, along);
  if (!(along.Magnitude() > 0))
    return 0;
  double u0, u1, v0, v1;
  BRepTools::UVBounds(side.support, u0, u1, v0, v1);
  const gp_Vec2d left =
    gp_Vec2d(-along.Y(), along.X()).Normalized() * (1e-4 * std::hypot(u1 - u0, v1 - v0));
  gp_Pnt point;
  gp_Vec su, sv;
  BRepAdaptor_Surface(side.support).D1(uv.X(), uv.Y(), point, su, sv);
  const gp_Vec across = side.faceJet(0.5).n.Crossed(side.tangent(0.5));
  for (const double sign : {1.0, -1.0})
  {
    const gp_Vec2d way = left * sign;
    if (BRepClass_FaceClassifier(side.support, uv.Translated(way), Precision::PConfusion()).State()
        == TopAbs_IN)
      return (su * way.X() + sv * way.Y()).Dot(across) > 0 ? -1 : 1;
  }
  return 0;
}

//! The fractions of a side at `pieces` even steps of its length, from its start to its end.
std::vector<double> evenSteps(const Side& side, int pieces)
{
  const int           fine = 16 * pieces;
  std::vector<double> length(fine + 1, 0);
  gp_Pnt              before = side.at(0);
  for (int i = 1; i <= fine; ++i)
  {
    const gp_Pnt here = side.at(static_cast<double>(i) / fine);
    length[i]         = length[i - 1] + before.Distance(here);
    before            = here;
  }
  std::vector<double> out{0};
  for (int k = 1; k < pieces; ++k)
  {
    const double want = length[fine] * k / pieces;
    const size_t i    = std::upper_bound(length.begin(), length.end(), want) - length.begin() - 1;
    const double part =
      length[i + 1] > length[i] ? (want - length[i]) / (length[i + 1] - length[i]) : 0;
    out.push_back((i + part) / fine);
  }
  out.push_back(1);
  return out;
}

//! A curve in the chart through `places` at `fractions`, for a curve parameter from `from` to
//! `to`.
std::function<gp_XY(double)> laidThrough(const std::vector<gp_XY>&  places,
                                         const std::vector<double>& fractions,
                                         double                     from,
                                         double                     to,
                                         double                     size)
{
  const int                                        count = static_cast<int>(places.size());
  const occ::handle<NCollection_HArray1<gp_Pnt2d>> points =
    new NCollection_HArray1<gp_Pnt2d>(1, count);
  const occ::handle<NCollection_HArray1<double>> at = new NCollection_HArray1<double>(1, count);
  for (int i = 0; i < count; ++i)
  {
    points->SetValue(i + 1, gp_Pnt2d(places[i]));
    at->SetValue(i + 1, fractions[i]);
  }
  Geom2dAPI_Interpolate through(points, at, false, BRepExt::tolerance::exactLength * size);
  through.Perform();
  if (!through.IsDone())
    fail(op, no_solution, "a side of the loop has no curve in its chart");
  const occ::handle<Geom2d_BSplineCurve> curve = through.Curve();
  return [curve, from, to](double t) { return curve->Value((t - from) / (to - from)).XY(); };
}

//! The Favard constant that bounds the best approximation of degree p splines on even knots of span
//! h: no more than K (h / pi)^(p + 1 - k) |f^(p + 1)| for the k-th derivative, K between 1 and
//! pi / 2 (Favard 1937; de Boor, A Practical Guide to Splines, 2001, XII). The larger is taken.
constexpr double favard = M_PI / 2;

//! A cost forecast using the higher derivatives of the osculating circle as a surrogate.
//! Those derivatives are not bounded by curvature for a general surface or chart. This is
//! a computation policy, not an approximation theorem or a necessary pole-count bound.
double netSpacing(double bend, int order, int degree, const Limits& limits)
{
  if (!(bend > 0))
    return HUGE_VAL;
  double span =
    M_PI * std::pow(fitShare * limits.g0 / (favard * std::pow(bend, degree)), 1.0 / (degree + 1));
  if (order > 0)
    span = std::min(
      span,
      M_PI * std::pow(fitShare * std::sin(limits.g1 * M_PI / 180) / favard, 1.0 / degree) / bend);
  if (order > 1)
    span =
      std::min(span, M_PI * std::pow(fitShare * limits.g2 / favard, 1.0 / (degree - 1)) / bend);
  return span;
}

//! A fill over a chart of the loop's own, for a loop whose rim's faces stand so that it is a
//! height over no plane (see overPlane). The loop is spanned by a disc of triangles that leaves
//! every side as its face asks and is laid flat conformally (proxyChart); a net free in space is
//! fitted over that chart, held to the sides at their places in it as over the mean plane: with
//! its tangent plane the face's and its bend across the face's where they are asked, which are
//! conditions on the derivatives in space and none on how the chart lies. It is trimmed by the
//! sides' curves in its parameters and kept only where S_u x S_v is shown to keep away from zero
//! inside the loop (regularMargin). Nothing where no net over the chart meets what is asked.
std::optional<Built> overChart(const Asked& asked, const Frame& frame)
{
  PatchPhase   phase("chart");
  const size_t count     = asked.sides.size();
  double       perimeter = 0;
  for (const Side& side : asked.sides)
    perimeter += side.length;
  std::vector<ChartRim>            rim;
  std::vector<std::vector<double>> steps(count);
  std::vector<size_t>              firstNode(count);
  std::vector<double>              sign(count);
  for (size_t s = 0; s < count; ++s)
  {
    const Side& side = asked.sides[s];
    const int   pieces =
      std::max(1, static_cast<int>(std::lround(chartPlaces * side.length / perimeter)));
    steps[s]     = evenSteps(side, pieces);
    firstNode[s] = rim.size();
    sign[s]      = leaving(side);
    for (int k = 0; k < pieces; ++k)
    {
      const double f = (steps[s][k] + steps[s][k + 1]) / 2;
      rim.push_back({side.at(steps[s][k]),
                     sign[s] == 0
                       ? gp_Vec(0, 0, 0)
                       : side.faceJet(f).n.Crossed(side.tangent(f)).Normalized() * sign[s]});
    }
  }
  std::vector<gp_Pnt> through = asked.points;
  for (const Side& guide : asked.guides)
    for (int i = 0; i <= guidePins; ++i)
      through.push_back(guide.at(static_cast<double>(i) / guidePins));
  int order = 0;
  for (const Side& side : asked.sides)
    order = std::max(order, side.order);
  const int            degree = order > 0 || !asked.guides.empty() ? 5 : 3;
  const ConformalChart chart  = proxyChart(rim, through);
  // The span each node asks of a net (netSpacing): inside, by how the disc bends there; along the
  // rim, by how the side's curve and its face bend, to the order the side asks.
  const std::vector<double> bends = nodeBends(chart);
  std::vector<double>       spacing(chart.point.size());
  for (size_t i = 0; i < spacing.size(); ++i)
    spacing[i] = std::min(asked.size, netSpacing(bends[i], 0, degree, asked.limits));
  for (size_t s = 0; s < count; ++s)
  {
    const Side& side = asked.sides[s];
    for (size_t j = 0; j < steps[s].size(); ++j)
    {
      const size_t node = j + 1 < steps[s].size() ? firstNode[s] + j : firstNode[(s + 1) % count];
      const double f    = steps[s][j];
      // The side's bend from the circle through its places beside this one, at even steps of its
      // length: its parameter may slow to nothing (a trailing edge) where its derivatives say
      // nothing of the bend, and a corner, at a side's end, is no bend.
      double curve = 0;
      if (steps[s].size() >= 3)
      {
        const size_t middle = std::clamp<size_t>(j, 1, steps[s].size() - 2);
        const gp_Pnt before = side.at(steps[s][middle - 1]), here = side.at(steps[s][middle]),
                     after = side.at(steps[s][middle + 1]);
        const double product =
          before.Distance(here) * here.Distance(after) * before.Distance(after);
        if (product > 0)
          curve = 2 * gp_Vec(here, before).Crossed(gp_Vec(here, after)).Magnitude() / product;
      }
      const double face =
        side.order > 0 ? sharpest(side.faceJet(std::clamp(f, cornerInset, 1 - cornerInset))) : 0;
      spacing[node] =
        std::min({spacing[node],
                  netSpacing(std::max({curve, face, bends[node]}), 0, degree, asked.limits),
                  netSpacing(std::max(face, bends[node]), side.order, degree, asked.limits)});
    }
  }
  // The spans a net over the chart needs at least: each way, over the chart's extent, the span
  // each node asks as the chart stretches it there, the finest across the other way. Where those
  // make more poles than a net free in space is let have, no net over the chart is tried.
  {
    const std::vector<double> stretch = nodeStretch(chart);
    constexpr int             bins    = 256;
    const auto                needed  = [&](int way, double from, double to) {
      std::vector<double> span(bins, HUGE_VAL);
      for (size_t i = 0; i < chart.flat.size(); ++i)
      {
        const int bin =
          std::clamp(static_cast<int>((chart.flat[i].Coord(way) - from) / (to - from) * bins),
                     0,
                     bins - 1);
        span[bin] = std::min(span[bin], spacing[i] * stretch[i]);
      }
      double spans = 0;
      for (const double one : span)
        if (std::isfinite(one))
          spans += (to - from) / bins / one;
      return std::max(fewestSpans, static_cast<int>(std::ceil(spans)));
    };
    gp_XY low = chart.flat.front(), high = low;
    for (const gp_XY& at : chart.flat)
    {
      low.SetCoord(std::min(low.X(), at.X()), std::min(low.Y(), at.Y()));
      high.SetCoord(std::max(high.X(), at.X()), std::max(high.Y(), at.Y()));
    }
    const gp_XY margin = (high - low) * marginShare;
    const int   spansU = needed(1, low.X() - margin.X(), high.X() + margin.X()),
              spansV   = needed(2, low.Y() - margin.Y(), high.Y() + margin.Y());
    if ((spansU + degree) * (spansV + degree) > freePoles)
      return std::nullopt;
  }
  std::vector<std::function<gp_XY(double)>> layout(count);
  for (size_t s = 0; s < count; ++s)
  {
    std::vector<gp_XY> places;
    for (size_t k = 0; k + 1 < steps[s].size(); ++k)
      places.push_back(chart.flat[firstNode[s] + k]);
    places.push_back(chart.flat[firstNode[(s + 1) % count]]);
    layout[s] = laidThrough(places, steps[s], asked.sides[s].from, asked.sides[s].to, asked.size);
  }
  std::vector<gp_XY> loop(chart.flat.begin(), chart.flat.begin() + rim.size());
  const Inside       within(loop, 0);
  // Guides and points where they lie on the disc in space.
  const auto onDisc = [&](const gp_Pnt& point) {
    double     distance = 0;
    const auto at       = onChart(chart, point, distance);
    if (!inside(loop, at, outsideSlack * asked.size))
      fail(op, invalid_input, "a guide or a point lies outside the boundary");
    return at;
  };
  std::vector<std::function<gp_XY(double)>> guideLayout;
  for (const Side& guide : asked.guides)
  {
    constexpr int       guidePlaces = 64;
    std::vector<gp_XY>  places;
    std::vector<double> fractions;
    for (int i = 0; i <= guidePlaces; ++i)
    {
      fractions.push_back(static_cast<double>(i) / guidePlaces);
      places.push_back(onDisc(guide.at(fractions.back())));
    }
    guideLayout.push_back(laidThrough(places, fractions, 0, 1, asked.size));
  }
  std::vector<gp_XY> pointPlaces;
  for (const gp_Pnt& point : asked.points)
    pointPlaces.push_back(onDisc(point));

  gp_XY low = loop.front(), high = low;
  for (const gp_XY& at : chart.flat)
  {
    low.SetCoord(std::min(low.X(), at.X()), std::min(low.Y(), at.Y()));
    high.SetCoord(std::max(high.X(), at.X()), std::max(high.Y(), at.Y()));
  }
  const double width = high.X() - low.X(), depth = high.Y() - low.Y(),
               larger = std::max(width, depth);
  const gp_XY margin(marginShare * width, marginShare * depth);

  const auto certify = [&](const occ::handle<Geom_BSplineSurface>& surface) {
    std::vector<Rim> rims;
    for (size_t s = 0; s < count; ++s)
      rims.push_back(
        asked.sides[s].free
          ? Rim{}
          : preimage(surface, asked.sides[s], frame, asked.limits.g0, asked.size, {}, layout[s]));
    closeGaps(asked.sides, rims);
    std::vector<gp_XY> around;
    for (size_t s = 0; s < count; ++s)
      for (int k = 0; k < 64; ++k)
        around.push_back(rims[s].onSurface->Value(asked.sides[s].parameter(k / 64.0)).XY());
    double turning = 0;
    for (size_t i = 0; i < around.size(); ++i)
      turning += around[i].Crossed(around[(i + 1) % around.size()]);
    Built built{3,
                {},
                measureAll(surface,
                           asked,
                           [&](const Side& side, double f, const gp_Pnt&) -> std::optional<Place> {
                             if (const auto s = asked.boundaryIndex(side))
                               return Place{rims[*s].onSurface->Value(side.parameter(f)).XY(),
                                            true};
                             // A guide is read at the surface's nearest from its chart place.
                             for (size_t g = 0; g < asked.guides.size(); ++g)
                               if (&side == &asked.guides[g])
                                 return Place{guideLayout[g](f), false};
                             return std::nullopt;
                           })};
    for (size_t s = 0; s < count; ++s)
      rims[s].deviation =
        std::max(rims[s].deviation, built.measured.values[asked.sides[s].input * 4]);
    built.face = trimmed(surface, asked, rims, turning < 0);
    for (size_t s = 0; s < count; ++s)
    {
      const Side& side = asked.sides[s];
      if (!side.free)
        built.spatialBoundary.push_back(
          {side.input, rims[s].onSurface, side.from, side.to, side.originalReversed});
      if (hasGaps(asked))
        built.gaps.push_back(side.free ? rims[s].onSurface : occ::handle<Geom2d_Curve>());
    }
    measureSpatialBoundary(built, asked, false);
    if (built.measured.within)
    {
      // The cells the face uses: those that hold a place inside the loop or of its rim.
      const auto used = [&](double u0, double u1, double v0, double v1) {
        for (const gp_XY& at : around)
          if (u0 <= at.X() && at.X() <= u1 && v0 <= at.Y() && at.Y() <= v1)
            return true;
        for (int i = 0; i <= 4; ++i)
          for (int j = 0; j <= 4; ++j)
            if (within(gp_XY(u0 + (u1 - u0) * i / 4, v0 + (v1 - v0) * j / 4)))
              return true;
        return false;
      };
      std::optional<std::string> fault;
      if (!regularMargin(surface, used))
        fault = "it is not shown regular over the loop's chart";
      else if (!BRepCheck_Analyzer(built.face).IsValid())
        fault = "it is not a valid face";
      else
        fault = insideFault(built.face, asked.points.empty() ? rimBend(asked) : 0);
      if (fault)
      {
        built.measured.within = false;
        built.measured.unkept = true;
        built.measured.worst  = *fault;
        built.measured.excess = unkeptExcess;
      }
    }
    return built;
  };

  // The net starts as the coarsest of the ladder over the chart and is made finer span by span
  // where a place it is held to is missed, as the Adjacent fit is (see adjacentOver): the nose
  // of a wing asks for spans far shorter than the rest of its chart does.
  std::vector<double> breaksU, breaksV;
  const auto even = [&](std::vector<double>& breaks, double from, double to, double extent) {
    const int spans =
      std::max(fewestSpans, static_cast<int>(std::lround(ladder[0] * extent / larger)));
    for (int i = 0; i <= spans; ++i)
      breaks.push_back(from + (to - from) * i / spans);
  };
  even(breaksU, low.X() - margin.X(), high.X() + margin.X(), width);
  even(breaksV, low.Y() - margin.Y(), high.Y() + margin.Y(), depth);
  const auto knots = [&](const std::vector<double>& breaks) {
    std::vector<double> out(degree, breaks.front());
    out.insert(out.end(), breaks.begin(), breaks.end());
    out.insert(out.end(), degree, breaks.back());
    return out;
  };
  std::optional<Built>      best;
  occ::handle<Geom_Surface> previous;
  int                       idle = 0, unkept = 0;
  for (double before = HUGE_VAL;;)
  {
    SplineNet  net{degree, degree, knots(breaksU), knots(breaksV), {}};
    const int  nu = net.countU(), nv = net.countV();
    const auto mean = [&](const std::vector<double>& knots, int i) {
      double sum = 0;
      for (int k = 1; k <= degree; ++k)
        sum += knots[i + k];
      return sum / degree;
    };
    for (int i = 0; i < nu; ++i)
      for (int j = 0; j < nv; ++j)
        net.poles.push_back(frame.origin.XYZ() + frame.x.XYZ() * mean(net.knotsU, i)
                            + frame.y.XYZ() * mean(net.knotsV, j));
    double smallest = HUGE_VAL;
    for (const auto* breaks : {&breaksU, &breaksV})
      for (size_t k = 1; k < breaks->size(); ++k)
        smallest = std::min(smallest, (*breaks)[k] - (*breaks)[k - 1]);
    std::vector<Target> targets;
    // Per target, the way along its side in space where it asks a tangent plane.
    std::vector<gp_Vec> along;
    for (size_t s = 0; s < count; ++s)
    {
      const Side& side = asked.sides[s];
      if (side.free)
        continue;
      const int places =
        std::clamp(static_cast<int>(std::ceil(placesPerSpan * side.length / smallest)),
                   fewestRimPlaces,
                   mostRimPlaces);
      const double step = 1.0 / places;
      // Next to a gap the face has its say at the end as along the rest (see overPlane).
      const auto quiet = [&](size_t other) {
        const Side& near = asked.sides[other % count];
        return near.free || near.support.IsSame(side.support);
      };
      const bool lone       = hasGaps(asked);
      const bool quietStart = lone && quiet(s + count - 1), quietEnd = lone && quiet(s + 1);
      const auto add = [&](double f, bool end, bool watch) {
        if (end && ((f == 0 && quietStart) || (f == 1 && quietEnd)))
        {
          f   = std::clamp(f, cornerInset, 1 - cornerInset);
          end = false;
        }
        const bool smooth = side.order > 0 && !end;
        Target     target{side.meet(f),
                      layout[s](side.parameter(f)),
                      smooth ? side.order : 0,
                          {},
                          {},
                      step * side.length / asked.size};
        target.watch = watch;
        if (smooth)
        {
          target.jet          = side.jet(f);
          target.easing       = side.easing(f);
          const double h      = std::min({step / 4, f, 1 - f});
          const gp_XY tangent = layout[s](side.parameter(f + h)) - layout[s](side.parameter(f - h));
          target.across       = gp_XY(-tangent.Y(), tangent.X()) / tangent.Modulus();
          if (sign[s] != 0)
            target.leave = side.faceJet(f).n.Crossed(side.tangent(f)).Normalized() * sign[s];
        }
        targets.push_back(target);
        along.push_back(smooth ? side.tangent(f).Normalized() : gp_Vec(0, 0, 0));
      };
      add(0, true, false);
      add(1, true, false);
      for (int k = 0; k < places; ++k)
      {
        add((k + 0.5) * step, false, false);
        if (k + 1 < places)
          add((k + 1) * step, false, true);
      }
    }
    for (size_t g = 0; g < asked.guides.size(); ++g)
    {
      const Side& guide = asked.guides[g];
      const int   places =
        std::clamp(static_cast<int>(std::ceil(placesPerSpan * guide.length / smallest)),
                   fewestRimPlaces,
                   mostRimPlaces);
      for (int i = 0; i <= 2 * places; ++i)
      {
        const double f = static_cast<double>(i) / (2 * places);
        targets.push_back(
          {guide.meet(f), guideLayout[g](f), 0, {}, {}, guide.length / places / asked.size});
        targets.back().watch = i % 2 != 0;
        along.push_back(gp_Vec(0, 0, 0));
      }
    }
    for (size_t i = 0; i < asked.points.size(); ++i)
    {
      targets.push_back({asked.points[i], pointPlaces[i], 0, {}, {}, 1});
      along.push_back(gp_Vec(0, 0, 0));
    }
    const occ::handle<Geom_BSplineSurface> surface = fit(
      net,
      std::vector<char>(net.poles.size(), 0),
      targets,
      asked.size,
      asked.limits,
      {},
      [&](double u, double v) { return within(gp_XY(u, v)) ? 1 : outsideWeight; },
      nullptr,
      {},
      &frame);
    Built        built  = certify(surface);
    const double change = previous.IsNull() ? HUGE_VAL
                                            : apart(previous,
                                                    surface,
                                                    2 * static_cast<int>(breaksU.size() - 1),
                                                    2 * static_cast<int>(breaksV.size() - 1),
                                                    low,
                                                    high,
                                                    within);
    previous            = surface;
    // As over the mean plane: a net that meets every limit is kept once the next finer one lies
    // within `settled` of it, and until then the finest that meets them.
    const bool split = built.measured.within || (best && best->measured.within);
    if (split)
    {
      if (change <= settled * asked.size)
      {
        if (!best->measured.within)
          best = std::move(built);
        break;
      }
      if (built.measured.within)
        best = std::move(built);
    }
    else if (!best || built.measured.excess < best->measured.excess)
      best = built;
    // Nets that come no nearer, more than the coarsest few may: what is left is not the net's.
    if (split)
    {
    }
    else if (built.measured.excess > stillFalling * before)
    {
      if (++idle > unruly)
      {
        best->measured.worst += ". A finer net came no nearer";
        break;
      }
    }
    else
    {
      idle   = 0;
      before = built.measured.excess;
    }
    // Finer in every span where the fit misses a place it was held to, or a side is missed
    // between them: halved, up to the finest net of the ladder each way.
    std::vector<char> splitU(breaksU.size() - 1, split), splitV(breaksV.size() - 1, split);
    const auto        mark = [&](const gp_XY& uv) {
      const auto span = [](const std::vector<double>& breaks, double at) {
        return std::clamp(
          static_cast<int>(std::upper_bound(breaks.begin(), breaks.end(), at) - breaks.begin()) - 1,
          0,
          static_cast<int>(breaks.size()) - 2);
      };
      splitU[span(breaksU, uv.X())] = splitV[span(breaksV, uv.Y())] = 1;
    };
    for (size_t t = 0; !split && t < targets.size(); ++t)
    {
      const Target& target = targets[t];
      const double  share  = target.watch ? watchShare : fitShare;
      bool          missed = surface->Value(target.uv.X(), target.uv.Y()).Distance(target.point)
                    > share * asked.limits.g0;
      if (!missed && target.order > 0)
      {
        const SurfaceJet mine   = surfaceJet(op, surface, target.uv.X(), target.uv.Y());
        const double     cosine = mine.n.Dot(target.jet->n);
        missed = std::acos(std::min(1.0, std::abs(cosine))) * 180 / M_PI > share * asked.limits.g1;
        if (!missed && target.order > 1)
        {
          const auto [bend, relative] =
            bendDifference(*target.jet, mine, along[t], cosine < 0 ? -1 : 1);
          missed = relative > share * asked.limits.g2
                   && bend * asked.size > share * asked.limits.g2 * flatShare;
        }
      }
      if (missed)
        mark(target.uv);
    }
    for (size_t s = 0; s < count && !split; ++s)
    {
      const Side& side = asked.sides[s];
      if (side.free)
        continue;
      Miss miss;
      miss.g0         = built.measured.values[side.input * 4];
      miss.g1         = built.measured.values[side.input * 4 + 1];
      miss.g2Relative = built.measured.values[side.input * 4 + 2];
      miss.g2         = built.measured.values[side.input * 4 + 3];
      if (miss.within(asked.limits, side.order, asked.size))
        continue;
      const int places =
        std::clamp(static_cast<int>(std::ceil(4 * placesPerSpan * side.length / smallest)),
                   fewestRimPlaces,
                   4 * mostRimPlaces);
      for (int i = 0; i <= places; ++i)
        mark(layout[s](side.parameter(static_cast<double>(i) / places)));
    }
    // A net that meets the edges with a surface that is not one to keep is too coarse for what
    // is asked along them, wherever that shows: a few finer nets, and no more.
    if (!split && built.measured.unkept && ++unkept <= mostUnkept)
    {
      std::fill(splitU.begin(), splitU.end(), 1);
      std::fill(splitV.begin(), splitV.end(), 1);
    }
    bool                               finer = false;
    std::array<std::vector<double>, 2> next;
    for (int way = 0; way < 2; ++way)
    {
      const auto& breaks = way ? breaksV : breaksU;
      const auto& marked = way ? splitV : splitU;
      next[way]          = {breaks.front()};
      int spans          = static_cast<int>(breaks.size()) - 1;
      for (size_t i = 0; i + 1 < breaks.size(); ++i)
      {
        if (marked[i] && spans < std::end(ladder)[-1])
        {
          next[way].push_back((breaks[i] + breaks[i + 1]) / 2);
          ++spans;
          finer = true;
        }
        next[way].push_back(breaks[i + 1]);
      }
    }
    // A net free in space is let have no more poles than over the mean plane (see freePoles).
    if ((next[0].size() + degree - 1) * (next[1].size() + degree - 1) > freePoles)
      break;
    breaksU = next[0];
    breaksV = next[1];
    if (!finer)
      break;
  }
  if (!best || !best->measured.within)
    return std::nullopt;
  return best;
}

//! A planar section guide carried along a rim, with support-derived Hermite extensions.
//! This is the translation member of the curve-network Boolean-sum construction.
//! It is only a candidate: all original rims, jets, guides and interior checks still apply.
std::optional<Built> sectionReference(const Asked& asked, const Frame& frame)
{
  if (asked.guides.size() != 1 || !asked.points.empty() || asked.sides.size() != 4)
    return {};
  std::vector<const Side*> rims;
  for (const Side& side : asked.sides)
    if (!side.free)
      rims.push_back(&side);
  if (rims.size() != 2)
    return {};
  if (rims[1]->order > rims[0]->order)
    std::swap(rims[0], rims[1]);
  const Side &rim = *rims[0], &foot = *rims[1], &guide = asked.guides.front();
  auto        profile = plain(guide, guide.length, false, asked.limits.g0);
  auto        rail    = plain(rim, rim.length, false, asked.limits.g0);
  if (profile.IsNull() || rail.IsNull())
    return {};
  std::vector<gp_Pnt> samples;
  for (int i = 1; i <= profile->NbPoles(); ++i)
    samples.push_back(profile->Pole(i));
  const Frame section = meanPlane(samples, guide.length);
  // The section must be planar; reject rather than moving any supplied guide.
  for (int i = 1; i <= profile->NbPoles(); ++i)
    if (std::abs(gp_Vec(section.origin, profile->Pole(i)).Dot(section.z)) > asked.limits.g0)
      return {};
  const auto intersect = [&](const Side& side) -> std::optional<double> {
    const auto height = [&](double f) { return gp_Vec(section.origin, side.at(f)).Dot(section.z); };
    double     a = 0, b = 1, ha = height(a), hb = height(b);
    if (ha * hb >= 0)
      return {};
    for (int i = 0; i < std::numeric_limits<double>::digits; ++i)
    {
      const double m = (a + b) / 2, h = height(m);
      if ((h >= 0) == (ha >= 0))
      {
        a  = m;
        ha = h;
      }
      else
        b = m;
    }
    return (a + b) / 2;
  };
  const auto ra = intersect(rim), rb = intersect(foot);
  if (!ra || !rb)
    return {};
  const gp_Pnt a = rim.meet(*ra), b = foot.meet(*rb);
  if (a.Distance(profile->Value(profile->FirstParameter()))
        + b.Distance(profile->Value(profile->LastParameter()))
      > a.Distance(profile->Value(profile->LastParameter()))
          + b.Distance(profile->Value(profile->FirstParameter())))
    profile->Reverse();
  gp_Pnt g0, g1;
  gp_Vec d0, d1, e0, e1;
  profile->D2(profile->FirstParameter(), g0, d0, e0);
  profile->D2(profile->LastParameter(), g1, d1, e1);
  const double h0 = a.Distance(g0), h1 = b.Distance(g1);
  if (h0 <= asked.limits.g0 || h1 <= asked.limits.g0)
    return {};
  const auto jet = [&](const Side&   side,
                       double        f,
                       const gp_Pnt& toward,
                       const gp_Vec& fallback) -> std::pair<gp_Vec, gp_Vec> {
    if (side.support.IsNull() || side.order == 0)
      return {fallback, gp_Vec(0, 0, 0)};
    const SurfaceJet j = side.jet(f);
    gp_Vec           t = j.n.Crossed(section.z);
    t.Normalize();
    if (t.Dot(gp_Vec(side.meet(f), toward)) < 0)
      t.Reverse();
    const double E = j.u.SquareMagnitude(), F = j.u.Dot(j.v), G = j.v.SquareMagnitude();
    const double u = (G * t.Dot(j.u) - F * t.Dot(j.v)) / (E * G - F * F),
                 v = (E * t.Dot(j.v) - F * t.Dot(j.u)) / (E * G - F * F);
    gp_Vec       e = j.uu * (u * u) + j.uv * (2 * u * v) + j.vv * (v * v);
    // Differentiate the plane-intersection constraint, preserving the normal acceleration.
    const gp_Vec w = section.z - j.n * section.z.Dot(j.n);
    e -= w * (e.Dot(section.z) / w.SquareMagnitude());
    return {t, e};
  };
  const auto [da, ea] = jet(rim, *ra, g0, d0);
  auto [db, eb]       = jet(foot, *rb, g1, -d1);
  if (foot.order > 0)
    db.Reverse(); // The section arrives at the second rim.
  else
    db = d1;
  const auto hermite = [](const gp_Pnt& p,
                          const gp_Vec& d,
                          const gp_Vec& e,
                          const gp_Pnt& q,
                          const gp_Vec& f,
                          const gp_Vec& g,
                          double        h) {
    NCollection_Array1<gp_Pnt> poles(1, 6);
    poles(1)   = p;
    poles(2)   = p.Translated(d * (h / 5));
    poles(3)   = p.Translated(d * (2 * h / 5) + e * (h * h / 20));
    poles(4)   = q.Translated(f * (-2 * h / 5) + g * (h * h / 20));
    poles(5)   = q.Translated(f * (-h / 5));
    poles(6)   = q;
    auto curve = GeomConvert::CurveToBSplineCurve(new Geom_BezierCurve(poles));
    NCollection_Array1<double> knots = curve->Knots();
    BSplCLib::Reparametrize(0, h, knots);
    curve->SetKnots(knots);
    return curve;
  };
  profile->IncreaseDegree(5);
  GeomConvert_CompCurveToBSplineCurve joined(hermite(a, da, ea, g0, d0, e0, h0));
  if (!joined.Add(profile, asked.limits.g0 * 1e-6, true, false, 3)
      || !joined.Add(hermite(g1, d1, e1, b, db, eb, h1), asked.limits.g0 * 1e-6, true, false, 3))
    return {};
  const auto                 cross = joined.BSplineCurve();
  NCollection_Array2<gp_Pnt> poles(1, rail->NbPoles(), 1, cross->NbPoles());
  for (int i = 1; i <= rail->NbPoles(); ++i)
    for (int j = 1; j <= cross->NbPoles(); ++j)
      poles(i, j) = gp_Pnt(rail->Pole(i).XYZ() + cross->Pole(j).XYZ() - a.XYZ());
  const occ::handle<Geom_BSplineSurface> surface = new Geom_BSplineSurface(poles,
                                                                           rail->Knots(),
                                                                           cross->Knots(),
                                                                           rail->Multiplicities(),
                                                                           cross->Multiplicities(),
                                                                           rail->Degree(),
                                                                           cross->Degree());
  ShapeAnalysis_Surface                  project(surface);
  const auto                             layout = [&](const gp_Pnt& p) {
    return project.ValueOfUV(p, asked.limits.g0 * 1e-3).XY();
  };
  std::vector<Rim> traces(asked.sides.size());
  for (size_t i = 0; i < asked.sides.size(); ++i)
    if (!asked.sides[i].free)
      traces[i] = preimage(surface, asked.sides[i], frame, asked.limits.g0, asked.size, layout);
  patch_side::closeGaps(asked.sides, traces);
  std::vector<gp_XY> loop;
  for (size_t i = 0; i < asked.sides.size(); ++i)
    loop.push_back(traces[i].onSurface->Value(asked.sides[i].from).XY());
  Built built{
    3,
    {},
    measureAll(surface,
               asked,
               [&](const Side& side, double f, const gp_Pnt&) -> std::optional<Place> {
                 const auto i = asked.boundaryIndex(side);
                 return i ? std::optional<Place>(
                              Place{traces[*i].onSurface->Value(side.parameter(f)).XY(), true})
                          : std::nullopt;
               })};
  if (!built.measured.within)
    return {};
  built.face     = trimmed(surface, asked, traces, enclosed(loop) < 0);
  built.crossing = Crossing{true, rim.order, foot.order};
  for (size_t i = 0; i < asked.sides.size(); ++i)
  {
    const Side& side = asked.sides[i];
    if (!side.free)
      built.spatialBoundary.push_back(
        {side.input, traces[i].onSurface, side.from, side.to, side.originalReversed});
    built.gaps.push_back(side.free ? traces[i].onSurface : occ::handle<Geom2d_Curve>());
  }
  measureSpatialBoundary(built, asked, false);
  if (!built.measured.within
      || !regularMargin(surface, [](double, double, double, double) { return true; })
      || !BRepCheck_Analyzer(built.face).IsValid()
      || insideFault(built.face, asked.points.empty() ? rimBend(asked) : 0, built.crossing))
    return {};
  return built;
}

Built overPlane(const Asked& asked, const Frame& flat, int direction = 0)
{
  // Corners where the faces on both sides are tangent and bend differently (a rim that runs
  // from a plane onto a round): there the surface has to change its bend, and a net that
  // is twice differentiable everywhere can only smear that over a span.
  std::vector<gp_Pnt> turns;
  std::vector<gp_Vec> along;
  for (size_t i = 0; i < asked.sides.size(); ++i)
  {
    const Side &a = asked.sides[i], &b = asked.sides[(i + 1) % asked.sides.size()];
    if (a.order < 1 || b.order < 1 || a.support.IsSame(b.support)
        || cornerConflict(a, b) <= asked.limits.g2)
      continue;
    turns.push_back(a.at(1));
    along.push_back(a.tangent(1));
    along.push_back(b.tangent(0));
  }
  // The net's lines through such a corner carry the change best when the rim crosses both
  // of them well: the net is turned in its plane to where the rim that is nearest to one of
  // its directions is furthest from it.
  // And two such lines that run close beside each other leave the surface next to nothing
  // to hold it between them: the corners are to lie well apart both ways.
  Frame frame = flat;
  if (!turns.empty())
  {
    const std::vector<gp_Pnt> outline = polygon(asked.sides);
    double                    best = -1, angle = 0;
    for (int step = 0; step < 90; ++step)
    {
      const double turn  = step * M_PI / 180;
      Frame        tried = flat;
      tried.x            = flat.x * std::cos(turn) + flat.y * std::sin(turn);
      tried.y            = flat.y * std::cos(turn) - flat.x * std::sin(turn);
      double worst       = 1;
      for (const gp_Vec& way : along)
        worst =
          std::min(worst, std::abs(std::sin(2 * std::atan2(way.Dot(tried.y), way.Dot(tried.x)))));
      gp_XY low = tried.project(outline.front()), high = low;
      for (const gp_Pnt& point : outline)
      {
        const gp_XY at = tried.project(point);
        low.SetCoord(std::min(low.X(), at.X()), std::min(low.Y(), at.Y()));
        high.SetCoord(std::max(high.X(), at.X()), std::max(high.Y(), at.Y()));
      }
      for (size_t i = 0; i < turns.size(); ++i)
        for (size_t k = i + 1; k < turns.size(); ++k)
        {
          const gp_XY apart = tried.project(turns[i]) - tried.project(turns[k]);
          worst             = std::min({worst,
                                        std::abs(apart.X()) / (high.X() - low.X()) * apartSpans,
                                        std::abs(apart.Y()) / (high.Y() - low.Y()) * apartSpans});
        }
      if (worst > best + 1e-9)
      {
        best  = worst;
        angle = turn;
      }
    }
    angle += direction * M_PI / 9;
    frame.x = flat.x * std::cos(angle) + flat.y * std::sin(angle);
    frame.y = flat.y * std::cos(angle) - flat.x * std::sin(angle);
  }
  std::vector<double> breaksU, breaksV;
  for (const gp_Pnt& turn : turns)
  {
    breaksU.push_back(frame.project(turn).X());
    breaksV.push_back(frame.project(turn).Y());
  }
  std::vector<gp_XY> loop;
  for (const gp_Pnt& point : polygon(asked.sides))
    loop.push_back(frame.project(point));
  // A loop that is refused over the mean plane, and whose rim's faces have normals that no plane
  // has all within `steepest` of its own, is a height over no plane at all: it is filled over a
  // chart of its own instead (see overChart). Every other loop is filled as it was, and one that
  // the chart does not fill either is refused as it was.
  const auto refused = [&](const std::string& refusal) -> Built {
    // The plane's check reads the faces between the ends; whether any plane would do is read
    // at the ends as well, where a rim's normals close their sweep (the nose of a wing).
    std::vector<gp_Vec> normals;
    for (const Side& side : asked.sides)
      for (int i = 0; direction == 0 && side.order > 0 && i <= gaussPlaces; ++i)
        normals.push_back(
          side.jet(std::clamp(static_cast<double>(i) / gaussPlaces, cornerInset, 1 - cornerInset))
            .n);
    if (direction == 0 && graphMargin(normals, {frame.x, frame.y, frame.z}) < steepest)
    {
      std::optional<Built> charted;
      try
      {
        charted = overChart(asked, flat);
      }
      catch (const Failure& failure)
      {
        if (failure.reason == "cancelled")
          throw;
      }
      catch (const Standard_Failure&)
      {
      }
      if (charted)
        return *charted;
    }
    fail(op, unsupported_geometry, refusal);
  };
  double                     area = 0;
  std::optional<std::string> crossed;
  try
  {
    area = enclosed(loop);
  }
  catch (const Failure& failure)
  {
    crossed = failure.message;
  }
  if (crossed)
    return refused(*crossed);
  const Inside within(loop, 0);
  // A fixed graph avoids spending fairness on tangential parameter drift. Keep the
  // spatial fit for steep supports, where its movable parameters approximate shape
  // energy better. The graph also makes the support's G2 conditions exactly linear.
  constexpr double graphNormal   = 0.8;
  bool             graphSupports = true;
  for (const Side& side : asked.sides)
    for (int i = 0; side.order > 0 && i < 16; ++i)
    {
      const double normal = std::abs(side.jet((i + 0.5) / 16).n.Dot(frame.z));
      graphSupports       = graphSupports && normal >= graphNormal;
      if (normal < steepest)
        return refused(
          "a face to be smooth with stands too steeply against the boundary's mean plane");
    }
  for (const Side& guide : asked.guides)
    for (int i = 0; i <= 16; ++i)
      if (!inside(loop, frame.project(guide.at(i / 16.0)), outsideSlack * asked.size))
        fail(op, invalid_input, "a guide leaves the boundary");
  for (const gp_Pnt& point : asked.points)
    if (!inside(loop, frame.project(point), outsideSlack * asked.size))
      fail(op, invalid_input, "a point lies outside the boundary");

  if (!asked.guides.empty())
  {
    Asked unguided = asked;
    unguided.guides.clear();
    const bool           gapped    = hasGaps(asked);
    std::optional<Built> candidate = gapped ? std::nullopt : spherical(unguided, flat);
    if (!candidate && !gapped)
      candidate = ruledRound(unguided, flat);
    if (!candidate && !gapped && unguided.points.empty()
        && squarish(unguided.sides, asked.limits.g1))
      candidate = foursided(unguided, flat);
    Built base = candidate && candidate->measured.within ? *candidate : overPlane(unguided, flat);
    if (base.measured.within)
    {
      const auto          surface = BRep_Tool::Surface(base.face);
      std::vector<double> misses;
      for (const Side& guide : asked.guides)
        misses.push_back(measure(surface, guide, [&](const Side&, double, const gp_Pnt& point) {
                           return std::optional<Place>(Place{frame.project(point), false});
                         }).g0);
      if (std::all_of(misses.begin(), misses.end(), [&](double miss) {
            return miss <= asked.limits.g0;
          }))
      {
        base.measured.values.insert(base.measured.values.begin() + asked.sides.size() * 4,
                                    misses.begin(),
                                    misses.end());
        return base;
      }
    }
  }

  gp_XY low = loop.front(), high = loop.front();
  for (const gp_XY& point : loop)
  {
    low.SetCoord(std::min(low.X(), point.X()), std::min(low.Y(), point.Y()));
    high.SetCoord(std::max(high.X(), point.X()), std::max(high.Y(), point.Y()));
  }
  const double width = high.X() - low.X(), depth = high.Y() - low.Y(),
               larger = std::max(width, depth);
  // Each way by its own extent: a tenth of the longer one on every side would leave a thin
  // strip a small part of its own surface.
  const gp_XY margin(marginShare * width, marginShare * depth);
  int         order = 0;
  for (const Side& side : asked.sides)
    order = std::max(order, side.order);
  // Degree 5 for tangency and guides; cubic for free boundaries without guides.
  // A net of degree 3 follows a rim's tangent planes as a cubic follows a curve, and the
  // corner of three rounds, a ball's eighth, took it 25 x 25 poles and came out rippled (a
  // largest bend of 0.56 for 0.25); of degree 5 it takes 13 x 13 and bends 0.31 at most.
  int        degree = order > 0 || !asked.guides.empty() ? 5 : 3;
  const auto typeOf = [](const Side& side) { return BRepAdaptor_Surface(side.support).GetType(); };
  const bool supported = std::all_of(asked.sides.begin(), asked.sides.end(), [](const Side& side) {
    return !side.support.IsNull();
  });
  const bool exchange =
    order > 0 && asked.guides.empty() && supported
    && ((asked.sides.size() == 1 && typeOf(asked.sides.front()) == GeomAbs_Cone)
        || (std::any_of(asked.sides.begin(),
                        asked.sides.end(),
                        [&](const Side& side) { return typeOf(side) == GeomAbs_Plane; })
            && std::any_of(asked.sides.begin(), asked.sides.end(), [&](const Side& side) {
                 return typeOf(side) != GeomAbs_Plane;
               })));

  // The fairest surface tangent to a ball is not the ball: of all that leave the rim in the
  // neighbours' planes, the one that bends least is flatter. Where the neighbours agree
  // about the bend all round the rim, their bend is wished for as well, and what comes of
  // that is taken where it is within the limits and one to keep; where they do not agree,
  // or it is not, the surface is tangent and no more.
  const double bend    = rimBend(asked);
  const bool   fairRim = order == 0 && asked.guides.empty() && asked.points.empty();
  // Where curves are only met, the fit is asked to come as near them as a fair surface over
  // its net does, a share of the loop's size and nothing of the tolerance, and the surface is
  // moved onto them afterwards: brought nearer by the fit, a coarse net ripples (z = sin 2a
  // on a radius of 6 bends 0.44 for 0.056 at 2e-5 and 0.06 at 5e-5). Where a face is gone
  // over into, the fit is asked for the tolerance itself: its net of degree 5 follows a rim
  // on a face to a millionth with few poles, and moved there afterwards took three times as
  // many.
  bool nearlyQuadratic = false;
  if (order == 0 && (!asked.guides.empty() || !asked.points.empty()))
  {
    // Recognition only: a quadratic needs neither the general guide fit's degree five
    // nor its extra poles. The real net is seeded below in the recognised degree.
    SplineNet probe{2,
                    2,
                    clampedKnots(low.X() - margin.X(), high.X() + margin.X(), 1, 2),
                    clampedKnots(low.Y() - margin.Y(), high.Y() + margin.Y(), 1, 2),
                    std::vector<gp_XYZ>(9)};
    if (polynomialHeight(asked, frame, probe))
      degree = 2;
    else if (!asked.guides.empty())
      nearlyQuadratic = polynomialHeight(asked, frame, probe, true);
  }
  Limits reach            = asked.limits;
  double inputUncertainty = 0;
  for (const Side& side : asked.sides)
    for (int k = 0; k < places; ++k)
    {
      const double f   = (k + .5) / places;
      inputUncertainty = std::max(inputUncertainty, side.at(f).Distance(side.meet(f)));
    }
  // Leave extra angle headroom when sewing must reparameterize a coarse input edge.
  if (inputUncertainty > 10 * asked.limits.g0)
    reach.g1 *= .5;
  GraphFit   graphMode = GraphFit::transverse;
  const auto build     = [&](bool wish, int nets) {
    std::optional<Built>      best;
    int                       tried = 0, unkept = 0;
    double                    placeOnly     = HUGE_VAL;
    int                       previousSpans = 0;
    occ::handle<Geom_Surface> previous;
    bool                      previousWithin = false, previousSpatial = false, settledOn = false;
    std::vector<Built>        meeting;
    // The corners of two sides that are only to be met at which a surface of the ladder left the
    // sides nearer the corner than the even places of its correction reach (see below).
    std::vector<char> forced(asked.sides.size(), 0);
    for (int rung : ladder)
    {
      const int spans = std::max(rung,
                                 asked.controls.quality == 2   ? 16
                                     : asked.controls.quality == 1 ? 10
                                                                   : 0);
      if (spans == previousSpans)
        continue;
      previousSpans = spans;
      // The last is for nets that only meet curves: a height of degree 3, which that many
      // spans cost no more than the net before it costs one that is free in space.
      if (nets-- <= 0 || (spans == std::end(ladder)[-1] && order > 0))
        break;
      auto along = [&](double extent) {
        // The complete-G2 retry refines both directions, including a narrow
        // direction whose aspect-ratio share never exceeds the coarse minimum.
        // Preserve established G0/G1 and transverse constructions.
        const int minimum =
          graphMode == GraphFit::complete
                ? static_cast<int>(std::lround(static_cast<double>(fewestSpans) * spans / ladder[0]))
                : fewestSpans;
        return std::max(minimum, static_cast<int>(std::lround(spans * extent / larger)));
      };
      const int spansU = along(width), spansV = along(depth);
      // Finer where the rim turns faster than the even net follows: by `turnPerSpan` of
      // turning to a span on the coarsest net, and in step with the net after that.
      std::vector<Finer> finerU, finerV;
      for (const Side& side : asked.sides)
        for (int i = 0; side.order > 0 && i < turningSteps; ++i)
        {
          const double f = (i + 0.5) / turningSteps, step = 0.25 / turningSteps;
          const double turning =
            side.jet(f - step).n.Angle(side.jet(f + step).n) / (2 * step * side.length);
          const double span = turnPerSpan * ladder[0] / spans / std::max(turning, 1e-12);
          if (span >= larger / spans)
            continue;
          const gp_XY at = frame.project(side.at(f));
          finerU.push_back({at.X(), span});
          finerV.push_back({at.Y(), span});
        }
      // A narrow opening needs longitudinal resolution near its ends too:
      // the G2 boundary layer can be as short as the opening is wide, even
      // when the boundary tangent hardly turns. Refine only the failed-G2
      // retry, and only where its even spans exceed that physical scale.
      if (graphMode == GraphFit::complete)
      {
        const double span = std::min(width, depth) / fewestSpans;
        const auto   ends = [&](double from, double to, int count, std::vector<Finer>& finer) {
          if ((to - from) / count <= span)
            return;
          for (double end : {from, to})
            for (int i = -fewestSpans; i <= fewestSpans; ++i)
              finer.push_back({end + i * span, span});
        };
        ends(low.X(), high.X(), spansU, finerU);
        ends(low.Y(), high.Y(), spansV, finerV);
      }
      SplineNet net{
        degree,
        degree,
        gradedKnots(low.X() - margin.X(), high.X() + margin.X(), spansU, degree, finerU, breaksU),
        gradedKnots(low.Y() - margin.Y(), high.Y() + margin.Y(), spansV, degree, finerV, breaksV),
            {}};
      // The net lies in the mean plane, each pole where the plane's own parameters put it,
      // and is free to rise from it: a height over the plane.
      const int  nu = net.countU(), nv = net.countV();
      const auto mean = [](const std::vector<double>& knots, int degree, int i) {
        double sum = 0;
        for (int k = 1; k <= degree; ++k)
          sum += knots[i + k];
        return sum / degree;
      };
      net.poles.resize(static_cast<size_t>(nu) * nv);
      for (int i = 0; i < nu; ++i)
        for (int j = 0; j < nv; ++j)
          net.poles[static_cast<size_t>(i) * nv + j] =
            (frame.origin.XYZ() + frame.x.XYZ() * mean(net.knotsU, degree, i)
             + frame.y.XYZ() * mean(net.knotsV, degree, j));
      const gp_XYZ up = frame.z.XYZ();
      const bool   seeded =
        order == 0 && (reach.g0 > asked.limits.g0 || !asked.guides.empty() || !asked.points.empty())
        && polynomialHeight(asked, frame, net, reach.g0 > asked.limits.g0);
      std::vector<Target> targets;
      for (const auto* curves : {&asked.sides, &asked.guides})
        for (const Side& side : *curves)
        {
          if (side.free)
            continue;
          // The same number of places in every span of the net that the curve passes
          // through, however briefly: a span it only clips has poles of its own, and
          // without a place there nothing but the energy would say where they go.
          std::vector<double> cuts{0, 1};
          gp_XY               before = frame.project(side.at(0));
          for (int i = 1; i <= crossingSteps; ++i)
          {
            const gp_XY here = frame.project(side.at(static_cast<double>(i) / crossingSteps));
            for (int way = 1; way <= 2; ++way)
            {
              const std::vector<double>& knots = way == 1 ? net.knotsU : net.knotsV;
              const double               a = before.Coord(way), b = here.Coord(way);
              for (auto knot = std::upper_bound(knots.begin(), knots.end(), std::min(a, b));
                   knot != knots.end() && *knot < std::max(a, b);
                   ++knot)
                cuts.push_back((i - 1 + (*knot - a) / (b - a)) / crossingSteps);
            }
            before = here;
          }
          std::sort(cuts.begin(), cuts.end());
          // Next to a gap, or where the next edge lies on the same face, nothing else is
          // asked of that end: the face has its say there as along the rest.
          const auto index = asked.boundaryIndex(side);
          const auto quiet = [&](size_t other) {
            const Side& near = asked.sides[other % asked.sides.size()];
            return near.free || near.support.IsSame(side.support);
          };
          const bool lone       = index && hasGaps(asked);
          const bool quietStart = lone && quiet(*index + asked.sides.size() - 1),
                     quietEnd   = lone && quiet(*index + 1);
          auto place            = [&](double f, double share, bool end, bool watch = false) {
            if (end && ((f == 0 && quietStart) || (f == 1 && quietEnd)))
            {
              f   = std::clamp(f, cornerInset, 1 - cornerInset);
              end = false;
            }
            const gp_Pnt point = side.meet(f);
            // Tangency and bend are not asked for at the ends, where two sides meet.
            const bool smooth = side.order > 0 && !end;
            Target     target{point,
                          frame.project(point),
                          smooth ? side.order : 0,
                                             {},
                                             {},
                          share * side.length / asked.size,
                          smooth && wish && side.order == 1};
            if (smooth)
            {
              target.jet          = side.jet(f);
              target.easing       = side.easing(f);
              const gp_Vec along3 = side.tangent(f);
              const gp_XY  along2(along3.Dot(frame.x), along3.Dot(frame.y));
              target.across = gp_XY(-along2.Y(), along2.X()) / along2.Modulus();
            }
            target.watch = watch;
            targets.push_back(target);
          };
          const double least = 1.0 / (crossingSteps * 16);
          place(0, least, true);
          place(1, least, true);
          for (size_t i = 0; i + 1 < cuts.size(); ++i)
          {
            const double from = cuts[i], to = cuts[i + 1];
            if (to - from < BRepExt::tolerance::parameterMerge)
              continue;
            for (int k = 0; k < placesPerSpan; ++k)
            {
              place(from + (to - from) * (k + 0.5) / placesPerSpan,
                    std::max(least, (to - from) / placesPerSpan),
                    false);
              // And between the places, where the surface is only read.
              place(from + (to - from) * (k + 1.0) / placesPerSpan,
                    std::max(least, (to - from) / placesPerSpan),
                    false,
                    true);
            }
          }
        }
      for (const gp_Pnt& point : asked.points)
        targets.push_back({point, frame.project(point), 0, {}, {}, 1.0 / spans});
      // Where two sides that are only to be met meet at a corner, a surface that holds both is
      // tangent there to the plane of their two tangents, whatever it is elsewhere. Within the
      // tolerance it need not be, and that plane is asked at the corner only once a surface of the
      // ladder has left the sides near the corner (see below). Not where they run on, where the two
      // give no plane.
      for (size_t i = 0; i < asked.sides.size(); ++i)
      {
        const Side &a = asked.sides[i], &b = asked.sides[(i + 1) % asked.sides.size()];
        if (!forced[i] || a.free || b.free || a.order > 0 || b.order > 0)
          continue;
        const gp_Vec ta = a.tangent(1).Normalized(), tb = b.tangent(0).Normalized();
        const gp_Vec normal = ta.Crossed(tb);
        if (normal.Magnitude() < std::sin(asked.limits.g1 * M_PI / 180))
          continue;
        const gp_Pnt corner = b.meet(0);
        Target       target{corner, frame.project(corner), 1, {}, {}, 1.0 / spans};
        SurfaceJet   jet{op, corner, ta, tb, {}, {}, {}, normal.Normalized()};
        target.jet = jet;
        targets.push_back(target);
      }
      const std::vector<char> fixed(net.poles.size(), 0);
      // Across a line through a corner where the rim's faces change their bend, the surface
      // may change its own: at the corner, where the rim needs it, and less and less the
      // farther along the line from it, so that the step stays at the rim and the surface
      // inside is of one bend across the line.
      std::vector<FitJump> jumps;
      for (const gp_Pnt& turn : turns)
      {
        const gp_XY corner = frame.project(turn);
        for (const bool acrossU : {true, false})
        {
          const std::vector<double>& knots = acrossU ? net.knotsV : net.knotsU;
          const double               here  = acrossU ? corner.Y() : corner.X();
          for (size_t k = degree; k + degree + 1 < knots.size(); ++k)
            for (int part = 0; knots[k + 1] > knots[k] && part < jumpPlaces; ++part)
            {
              const double along = knots[k] + (knots[k + 1] - knots[k]) * (part + 0.5) / jumpPlaces;
              const double away  = std::abs(along - here) / (jumpReach * asked.size);
              jumps.push_back({acrossU,
                               acrossU ? corner.X() : corner.Y(),
                               along,
                               jumpWeight * bendingWeight * asked.size * std::min(1.0, away * away)
                                 * (knots[k + 1] - knots[k]) / jumpPlaces});
            }
        }
      }
      // Fair inside the loop. Outside it the surface is only kept from running wild: counted
      // in full there, it would be held to the loop from both sides, like a plate clamped
      // along a line, and kink across it.
      const bool heightGraph =
        order > 0 && graphSupports
        && std::all_of(targets.begin(), targets.end(), [&](const Target& target) {
             return target.order == 0 || std::abs(target.jet->n.Dot(frame.z)) >= graphNormal;
           });
      const bool spatial =
        !heightGraph && order > 0 && asked.guides.empty() && nu * nv <= freePoles;
      // A net free in space and a height over the plane are fair by different measures and
      // come nearer different surfaces: the one is not a finer net of the other, and where
      // nets free in space have met the limits, a height is not tried after them.
      if (!previous.IsNull() && spatial != previousSpatial && best && best->measured.within)
        break;
      occ::handle<Geom_BSplineSurface> surface =
        seeded
              ? surfaceOf(net)
              : fit(
              net,
              fixed,
              targets,
              asked.size,
              reach,
              {},
              [&](double u, double v) { return within(gp_XY(u, v)) ? 1 : outsideWeight; },
              spatial ? nullptr : &up,
              jumps,
              order > 0 || asked.controls.open || !asked.guides.empty() || !asked.points.empty()
                    ? &frame
                    : nullptr,
              asked.flatness,
              fairRim,
              heightGraph ? graphMode : GraphFit::none);
      // The fit leaves of a curve what its net cannot follow fairly; the surface is then
      // moved onto the curves by the least movement of its poles, with knots added where
      // it takes them.
      // Not where the fit has left the tangent planes or the bends outside their limits: a
      // finer net is then tried in any case.
      // A boundary-only height graph is fitted by its weighted fair objective and corrected onto
      // the original rim: multiplier corrections amplify nearly dependent conditions.
      // Supported surfaces keep their original tangent/bend fit and correction gate.
      const bool moved =
        (fairRim || order > 0 || reach.g0 > asked.limits.g0)
        && measureAll(surface,
                          {asked.sides,
                       asked.guides,
                       asked.points,
                       {HUGE_VAL, asked.limits.g1, asked.limits.g2},
                       asked.size,
                       asked.flatness,
                       asked.controls},
                      [&](const Side&, double, const gp_Pnt& point) {
                        return std::optional<Place>(Place{frame.project(point), false});
                      })
             .within;
      // What the correction leaves at its own places where it had no room to make the net
      // finer, the net being as fine as it may make it already.
      double left = 0;
      if (moved)
      {
        const size_t knots = net.knotsU.size() + net.knotsV.size();
        // Per side, the places near its corners that the correction holds as well: those the
        // measure reads there (see measuredFractions) where the surface it leaves misses them.
        std::vector<std::vector<double>> nearCorners(asked.sides.size());
        const auto                       onto = [&](double span) {
          std::vector<RimPlace> onto;
          for (const auto* curves : {&asked.sides, &asked.guides})
            for (const Side& side : *curves)
            {
              if (side.free)
                continue;
              const int count = std::clamp(static_cast<int>(rimPlacesPerSpan * side.length / span),
                                           fewestRimPlaces,
                                           mostRimPlaces);
              std::vector<double> fractions;
              for (int i = 0; i <= count; ++i)
                fractions.push_back(static_cast<double>(i) / count);
              if (curves == &asked.sides)
                fractions.insert(fractions.end(),
                                 nearCorners[&side - asked.sides.data()].begin(),
                                 nearCorners[&side - asked.sides.data()].end());
              for (const double f : fractions)
              {
                const gp_Pnt point = side.meet(f);
                gp_XY        across(0, 0);
                if (side.order > 0 && f > 0 && f < 1)
                {
                  const gp_Vec along3 = side.tangent(f);
                  const gp_XY  along2(along3.Dot(frame.x), along3.Dot(frame.y));
                  across = gp_XY(-along2.Y(), along2.X()) / along2.Modulus();
                }
                onto.push_back({point, frame.project(point), across, side.order > 1});
              }
            }
          for (const gp_Pnt& point : asked.points)
            onto.push_back({point, frame.project(point), gp_XY(0, 0), false});
          return onto;
        };
        const auto correct = [&] {
          return correctRim(
            op,
            net,
            onto,
            up,
            std::min(rimShare * asked.limits.g0, exchange ? exchangeWithin : HUGE_VAL),
            std::min(enoughShare * asked.limits.g0, exchange ? 2 * exchangeWithin : HUGE_VAL),
            std::min(mostRimSpans,
                     (order > 0 ? rimFiner : rimFinerMet) * std::max(spansU, spansV)));
        };
        double worst = correct();
        // The correction's even places reach a corner of two sides that are only to be met no
        // nearer than its step: where the surface leaves a side nearer the corner than that by
        // more than the tolerance, those places are held as well and the correction made again,
        // and finer nets are asked the corner's tangent plane (see above).
        {
          const occ::handle<Geom_BSplineSurface> corrected = surfaceOf(net);
          bool                                   more      = false;
          const size_t                           count     = asked.sides.size();
          for (size_t s = 0; s < count; ++s)
          {
            const Side& side = asked.sides[s];
            if (side.free || side.order > 0)
              continue;
            for (const double f : measuredFractions())
            {
              const bool start = f < 0.5 / places, end = f > 1 - 0.5 / places;
              // The corner at that end, numbered by the side before it.
              const size_t corner = start ? (s + count - 1) % count : s;
              if ((!start && !end) || asked.sides[start ? corner : (s + 1) % count].order > 0)
                continue;
              const gp_Pnt point = side.meet(f);
              const gp_XY  uv    = nearest(corrected, frame.project(point), point);
              if (corrected->Value(uv.X(), uv.Y()).Distance(point) > asked.limits.g0)
              {
                nearCorners[s].push_back(f);
                forced[corner] = 1;
                more           = true;
              }
            }
          }
          if (more)
            worst = correct();
        }
        if (net.knotsU.size() + net.knotsV.size() == knots)
          left = worst;
        surface = surfaceOf(net);
      }
      std::vector<Rim> rims;
      for (const Side& side : asked.sides)
        rims.push_back(side.free
                         ? Rim{}
                         : preimage(surface,
                                    side,
                                    frame,
                                    std::min(asked.limits.g0, exchange ? exchangeWithin : HUGE_VAL),
                                    asked.size));
      closeGaps(asked.sides, rims);
      // A side is read where the face holds it; a guide and a point at the surface's nearest.
      Built built{3,
                      {},
                  measureAll(surface, asked, [&](const Side& side, double f, const gp_Pnt& point) {
                    const auto s = asked.boundaryIndex(side);
                    return std::optional<Place>(
                      s ? Place{rims[*s].onSurface->Value(side.parameter(f)).XY(), true}
                        : Place{frame.project(point), false});
                  })};
      for (size_t s = 0; s < rims.size(); ++s)
        rims[s].deviation =
          std::max(rims[s].deviation, built.measured.values[asked.sides[s].input * 4]);
      // What is not within the limits is kept for what it says of the miss.
      built.face = trimmed(surface, asked, rims, area < 0);
      if (fairRim)
        for (size_t i = 0; i < rims.size(); ++i)
          if (!asked.sides[i].free)
            built.spatialBoundary.push_back({asked.sides[i].input,
                                             rims[i].onSurface,
                                             asked.sides[i].from,
                                             asked.sides[i].to,
                                             asked.sides[i].originalReversed});
      for (size_t i = 0; hasGaps(asked) && i < rims.size(); ++i)
        built.gaps.push_back(asked.sides[i].free ? rims[i].onSurface : occ::handle<Geom2d_Curve>());
      // The actual trim must meet the original spatial edge too. Let the existing finite
      // ladder react to a missed trace instead of accepting a coarse preimage and refusing
      // only after placement. The final world-coordinate measurement remains authoritative.
      measureSpatialBoundary(built, asked, false);
      // Within every limit along its edges is not yet a patch: see BRepPatchExt_Inside.hxx. A net
      // that meets the edges with a surface that is not one to keep is too coarse for what
      // is asked along them, and what it cannot follow shows inside the loop; a finer one
      // may do. Points inside may ask for any bend, so with those only a fold is a fault.
      if (built.measured.within)
        if (const auto fault = BRepCheck_Analyzer(built.face).IsValid()
                                     ? insideFault(built.face, asked.points.empty() ? bend : 0)
                                     : std::optional<std::string>("it is not a valid face"))
        {
          built.measured.within = false;
          built.measured.unkept = true;
          built.measured.worst  = *fault;
          built.measured.excess = unkeptExcess;
        }
      const double change    = previous.IsNull() || spatial != previousSpatial
                                     ? HUGE_VAL
                                     : apart(previous,
                                      surface,
                                      // Twice to a span of the finer net each way: the
                                      // difference of the two is a spline of its spans.
                                      2 * spansU,
                                      2 * spansV,
                                      low,
                                      high,
                                      within);
      const bool   wasWithin = previousWithin;
      previous               = surface;
      previousWithin         = built.measured.within;
      previousSpatial        = spatial;
      // A net that meets every limit is not yet the surface asked for: that is the fairest of
      // all that meet them, which the nets of the ladder come nearer to as they become finer,
      // and a coarse net can meet the edges and still stand far from it inside (a seat of
      // eight sides: 9.2 off on six spans, 0.17 on every finer net). Two nets in a row that
      // lie within `settled` of each other have come as near it as finer ones would take
      // them: then the coarser of them that meets the limits is kept. Until then the finest
      // that meets them is, should the ladder or the nets allowed end first.
      const double before = best ? best->measured.excess : HUGE_VAL;
      if (built.measured.within)
        meeting.push_back(built);
      if (built.measured.within || (best && best->measured.within))
      {
        if (change <= settled * asked.size)
        {
          if (built.measured.within && !wasWithin)
            best = built;
          settledOn = true;
          break;
        }
        if (built.measured.within)
          best = built;
        continue;
      }
      if (built.measured.excess < before)
        best = built;
      // A few finer nets for a surface that is not one to keep, and no more: after those
      // it is the surface asked for that is not one to keep.
      if (built.measured.unkept)
      {
        if (++unkept >= mostUnkept)
          break;
        continue;
      }
      // A net that the fit left within its tangent planes and bends, and that was moved
      // onto its curves, is read at the places it was moved onto, and between them. Where
      // the fitted net was already as fine as the correction may make it, and the
      // correction's solve still leaves its own places further off than the tolerance, what
      // is left is what the fit left and the solve could not take, a share of the size that
      // finer nets of the ladder, fitted finer than their correction may go as well, do not
      // bring below it: the corner of rounds 2 / 3 / 5 at a hundred times its size (nets of
      // 25 x 20 and 31 x 26 spans for corrections bounded at 12 and 20) is left 3.7e-6 off by
      // the first and 6.0e-6 by the second (the limit 1e-6). Where the correction could make
      // the net finer, a finer net may meet them (a guide across a hole in a plate is met on
      // the third, at 32 x 40 spans).
      if (!fairRim && moved && left > asked.limits.g0)
      {
        best->measured.worst +=
          ". The bounded rim correction leaves the curves further off than the tolerance";
        break;
      }
      if (moved)
      {
        if (!fairRim && built.measured.excess > stillFalling * placeOnly)
        {
          best->measured.worst += ". A finer net came no nearer";
          break;
        }
        placeOnly = std::min(placeOnly, built.measured.excess);
      }
      // The coarsest nets miss by what they cannot follow at all, and by no rule, unless
      // it is by the very same.
      // But not, where curves are only met and the nets are cheap, while the best of them is
      // within twice its limits: with each net the fit stops where its own places are met,
      // and the misses of a ladder are not a falling series (one closed curve 0.5 sin 3a at a
      // tolerance of 1e-5 misses by 1.35, 4.5, 1.72 and 3.64 times that on the first four
      // nets and is met on the fifth).
      // Full Hessian errors need not fall monotonically as the narrow direction
      // gains spans. Exhaust its existing finite ladder rather than stop on one rise.
      if (!fairRim && graphMode != GraphFit::complete
          && ((++tried > unruly && built.measured.excess > stillFalling * before
               && (order > 0 || best->measured.excess > nearly))
              || std::abs(built.measured.excess - before) < same * before))
      {
        best->measured.worst += ". A finer net came no nearer";
        break;
      }
    }
    // Where no two nets in a row have come to agree, the ladder has not shown which of the
    // surfaces it made is near the fairest that meets the limits; of those that meet them, the
    // one that bends least is nearest it by the energy it is the least of.
    if (!settledOn && meeting.size() > 1)
    {
      double least = HUGE_VAL;
      for (const Built& candidate : meeting)
        if (const double bending = bendingInside(candidate.face); bending < least)
        {
          least = bending;
          best  = candidate;
        }
    }
    return *best;
  };
  if (nearlyQuadratic)
  {
    // Interpolating parabolic points does not make an exact parabola. Use the nearby
    // quadratic as a start on a cubic net, then actually meet the interpolated curves
    // by the same measured rim correction. It is not degree-two recognition.
    degree             = 3;
    reach              = Limits{fitReach * asked.size, asked.limits.g1, asked.limits.g2};
    const Built guided = build(false, static_cast<int>(std::size(ladder)));
    if (guided.measured.within)
      return guided;
    degree = 5;
    reach  = asked.limits;
  }
  // G1 on independently designed spline sheets leaves transverse curvature free.
  // Borrowing each sheet's G2 as a soft target can manufacture alternating bends in
  // an otherwise convex hole. Minimise the requested energy without those targets.
  // Analytic round supports retain their shape wishes (for example spherical caps).
  const bool splineRim = std::all_of(asked.sides.begin(), asked.sides.end(), [](const Side& side) {
    return !side.support.IsNull()
           && BRepAdaptor_Surface(side.support).GetType() == GeomAbs_BSplineSurface;
  });
  bool       wishable =
    !splineRim && std::any_of(asked.sides.begin(), asked.sides.end(), [](const Side& side) {
      return side.order == 1;
    });
  for (size_t i = 0; wishable && i < asked.sides.size(); ++i)
  {
    const Side &a = asked.sides[i], &b = asked.sides[(i + 1) % asked.sides.size()];
    if (a.order < 1 || b.order < 1 || a.support.IsSame(b.support))
      continue;
    const SurfaceJet before = a.jet(1 - cornerInset), after = b.jet(cornerInset);
    Miss             miss;
    std::tie(miss.g2, miss.g2Relative) = bendDifference(before,
                                                        after,
                                                        a.tangent(1 - cornerInset).Normalized(),
                                                        before.n.Dot(after.n) < 0 ? -1 : 1);
    wishable = miss.within(asked.limits, 2, asked.size) || wishAcrossCorners;
  }
  if (wishable)
  {
    const Built wished = build(true, wishNets);
    if (wished.measured.within)
      return wished;
  }
  const Built strict = build(false, direction == 0 ? static_cast<int>(std::size(ladder)) : 1);
  if (!strict.measured.within && order == 2 && graphSupports)
  {
    // Approximate position/normal constraints do not guarantee the two Hessian
    // components omitted by the transverse row. Retry the same accuracy contract
    // with all three components, preserving already delivered constructions.
    graphMode = GraphFit::complete;
    try
    {
      const Built complete = build(false, direction == 0 ? static_cast<int>(std::size(ladder)) : 1);
      if (complete.measured.within)
        return complete;
    }
    catch (const Failure& failure)
    {
      if (failure.reason != no_solution)
        throw;
    }
    graphMode = GraphFit::transverse;
  }
  if (!strict.measured.within && !asked.guides.empty())
  {
    try
    {
      if (const auto reference = sectionReference(asked, frame))
        return *reference;
    }
    catch (const Failure& failure)
    {
      if (failure.reason == "cancelled")
        throw;
    }
    catch (const Standard_Failure&)
    {
    }
  }
  if (order > 0)
  {
    // A near miss can occur between the fit's own places even when every asked and
    // watched condition was met. Retry once with extra angular headroom; the request
    // and every delivered contact/interior/validity check keep their original limits.
    if (order == 1 && direction == 0 && !strict.measured.within && !strict.measured.unkept
        && strict.measured.excess < 2 && asked.guides.empty() && asked.points.empty())
    {
      reach.g1 *= .5;
      try
      {
        const Built tighter = build(false, static_cast<int>(std::size(ladder)));
        if (tighter.measured.within)
          return tighter;
      }
      catch (const Failure& failure)
      {
        if (failure.reason != no_solution)
          throw;
      }
    }
    if (direction == 0 && strict.measured.within && !turns.empty() && asked.guides.empty()
        && asked.points.empty())
    {
      // Repeated-radius rounds retain their clearance frame. A turned retry on
      // those supports added CPU time without a candidate worth keeping.
      std::vector<TopoDS_Face> rounds;
      std::vector<double>      radii;
      bool                     repeated = false;
      for (const Side& side : asked.sides)
      {
        if (side.support.IsNull() || BRepAdaptor_Surface(side.support).GetType() != GeomAbs_Cylinder
            || std::any_of(rounds.begin(), rounds.end(), [&](const TopoDS_Face& face) {
                 return face.IsSame(side.support);
               }))
          continue;
        const double radius = BRepAdaptor_Surface(side.support).Cylinder().Radius();
        repeated = repeated || std::any_of(radii.begin(), radii.end(), [&](double other) {
                     return std::abs(other - radius) < 1e-9 * std::max(other, radius);
                   });
        rounds.push_back(side.support);
        radii.push_back(radius);
      }
      const InsideCurvature before = insideCurvature(strict.face);
      if (!repeated && before.wrongSign() > .34 * bend)
      {
        std::optional<Built> turned;
        try
        {
          turned = overPlane(asked, flat, 1);
        }
        catch (const Failure&)
        {
          return strict;
        }
        if (turned->measured.within)
        {
          const InsideCurvature after = insideCurvature(turned->face);
          if (after.wrongSign() < before.wrongSign() && after.tightest() <= before.tightest()
              && insideStep(turned->face) <= insideStep(strict.face))
            return *turned;
        }
      }
    }
    return strict;
  }
  // Where curves are only met there is a second surface: fitted only as near as a fair one
  // over its net comes, a share of the loop's size, and moved onto the curves afterwards.
  // Brought nearer by the fit, a coarse net ripples to meet a curve (z = sin 2a on a radius
  // of 6 bends 0.25 at 1e-6 for the 0.056 of the fair surface); moved, it does not (0.060).
  // But for other loops the fit at the tolerance is the fairer (a C of eight sides: 0.06,
  // moved 0.45). Of the two that meet the limits the fairer is kept, the one fitted at the
  // tolerance unless the other bends clearly less. Not for a loop with a corner where its
  // sides meet at an angle: no surface follows that without a crease, and moving one onto
  // it only makes its net finer there (a polygon of 16 sides in space: 99 x 99 poles).
  for (size_t i = 0; asked.guides.empty() && i < asked.sides.size(); ++i)
  {
    const Side &a = asked.sides[i], &b = asked.sides[(i + 1) % asked.sides.size()];
    if (a.tangent(1).Angle(b.tangent(0)) * 180 / M_PI > asked.limits.g1)
      return strict;
  }
  reach       = Limits{fitReach * asked.size, asked.limits.g1, asked.limits.g2};
  Built loose = build(false, static_cast<int>(std::size(ladder)));
  if (!loose.measured.within && !strict.measured.within && degree == 3 && asked.guides.empty()
      && asked.points.empty())
  {
    degree = 5;
    loose  = build(false, static_cast<int>(std::size(ladder)));
  }
  if (!loose.measured.within)
    return strict;
  if (!strict.measured.within)
    return loose;
  return tightestInside(loose.face) < clearlyFairer * tightestInside(strict.face) ? loose : strict;
}

// --- Adjacent over a conformal chart ---------------------------------------------------------

//! The triangulation that is flattened: at least this many triangles across the loop's
//! narrowest part (its area over its length round).
constexpr double trianglesAcross = 4;

//! How far the triangulation may stand off the face, as a share of the tolerance: where the
//! face bends its triangles are as small as that makes them, and the chart is as near the
//! face's own as its triangles are to the face.
constexpr double chartDeflection = 0.1;

//! The Adjacent fit over the chart: a net made finer span by span where what it is held to is
//! not met, from the coarsest of the ladder.
Built adjacentOver(const Asked&                                  asked,
                   const Frame&                                  frame,
                   const Built&                                  natural,
                   const ConformalChart&                         chart,
                   const std::vector<occ::handle<Geom2d_Curve>>& held,
                   double                                        area)
{
  const auto layout = [&](size_t s) {
    return std::function<gp_XY(double)>(
      [&, s](double t) { return chart.map(held[s]->Value(t).XY()); });
  };

  gp_XY low = chart.flat.front(), high = low;
  for (const gp_XY& at : chart.flat)
  {
    low.SetCoord(std::min(low.X(), at.X()), std::min(low.Y(), at.Y()));
    high.SetCoord(std::max(high.X(), at.X()), std::max(high.Y(), at.Y()));
  }
  std::vector<gp_XY> loop;
  std::vector<char>  onBoundary(chart.flat.size(), 0);
  for (int b : chart.boundary)
  {
    loop.push_back(chart.flat[b]);
    onBoundary[b] = 1;
  }
  const Inside within(loop, 0);
  const double width = high.X() - low.X(), depth = high.Y() - low.Y(),
               larger = std::max(width, depth);
  const gp_XY margin(marginShare * width, marginShare * depth);
  int         order = 0;
  for (const Side& side : asked.sides)
    order = std::max(order, side.order);
  // As over the mean plane: degree 5 for tangency, cubic where edges are only met.
  const int  degree = order > 0 ? 5 : 3;
  const auto source = BRep_Tool::Surface(natural.face);
  // Interior constraints belong to the Natural fill's chart as well. Keep that
  // correspondence fixed while fitting; reprojecting onto each candidate would
  // allow a constraint to slide away from the chart that defines the flow.
  ShapeAnalysis_Surface sourceProject(source);
  const auto            interiorPlace = [&](const gp_Pnt& point) {
    return chart.map(sourceProject.ValueOfUV(point, BRepExt::tolerance::exactLength).XY());
  };

  // Where the result is read inside the loop: every node and the middle of every triangle,
  // each with the Natural fill's normal there, which a fill that does not fold keeps to.
  struct Check
  {
    gp_XY  flat;
    gp_Vec normal;
  };

  std::vector<Check> checks;
  const auto         check = [&](const gp_XY& flat, const gp_XY& uv) {
    gp_Pnt point;
    gp_Vec su, sv;
    source->D1(uv.X(), uv.Y(), point, su, sv);
    checks.push_back({flat, su.Crossed(sv)});
  };
  for (size_t i = 0; i < chart.flat.size(); ++i)
    check(chart.flat[i], chart.source[i]);
  for (const auto& t : chart.triangles)
    check((chart.flat[t[0]] + chart.flat[t[1]] + chart.flat[t[2]]) / 3,
          (chart.source[t[0]] + chart.source[t[1]] + chart.source[t[2]]) / 3);

  // The net starts as the coarsest of the ladder and is made finer span by span where what
  // it is held to is not met: the knots follow the loop's geometry and the tolerance.
  std::vector<double> breaksU, breaksV;
  const auto even = [&](std::vector<double>& breaks, double from, double to, double extent) {
    const int spans =
      std::max(fewestSpans, static_cast<int>(std::lround(ladder[0] * extent / larger)));
    for (int i = 0; i <= spans; ++i)
      breaks.push_back(from + (to - from) * i / spans);
  };
  even(breaksU, low.X() - margin.X(), high.X() + margin.X(), width);
  even(breaksV, low.Y() - margin.Y(), high.Y() + margin.Y(), depth);
  const auto knots = [&](const std::vector<double>& breaks) {
    std::vector<double> out(degree, breaks.front());
    out.insert(out.end(), breaks.begin(), breaks.end());
    out.insert(out.end(), degree, breaks.back());
    return out;
  };

  struct Checked
  {
    Built             built;
    std::vector<char> folded;
    bool              folds;
  };

  const auto certify = [&](const occ::handle<Geom_BSplineSurface>& surface, bool spatial = false) {
    std::vector<Rim>   rims;
    std::vector<gp_XY> around;
    for (size_t s = 0; s < asked.sides.size(); ++s)
      rims.push_back(
        asked.sides[s].free
          ? Rim{}
          : preimage(surface, asked.sides[s], frame, asked.limits.g0, asked.size, {}, layout(s)));
    closeGaps(asked.sides, rims);
    for (size_t s = 0; s < asked.sides.size(); ++s)
      for (int k = 0; k < 16; ++k)
        around.push_back(rims[s].onSurface->Value(asked.sides[s].parameter(k / 16.0)).XY());
    Built built{3,
                {},
                measureAll(surface,
                           asked,
                           [&](const Side& side, double f, const gp_Pnt&) -> std::optional<Place> {
                             if (const auto s = asked.boundaryIndex(side))
                               return Place{rims[*s].onSurface->Value(side.parameter(f)).XY(),
                                            true};
                             // Guides are measured independently by spatial projection.
                             return std::nullopt;
                           })};
    for (size_t s = 0; s < rims.size(); ++s)
      rims[s].deviation =
        std::max(rims[s].deviation, built.measured.values[asked.sides[s].input * 4]);
    built.face = trimmed(surface, asked, rims, enclosed(around) < 0);
    if (spatial || !natural.spatialBoundary.empty())
      for (size_t i = 0; i < rims.size(); ++i)
        if (!asked.sides[i].free)
          built.spatialBoundary.push_back({asked.sides[i].input,
                                           rims[i].onSurface,
                                           asked.sides[i].from,
                                           asked.sides[i].to,
                                           asked.sides[i].originalReversed});
    // Inside the loop: no fold, and how square U and V are.
    std::vector<double> angles;
    double              least = HUGE_VAL, most = 0;
    std::vector<char>   folded(checks.size(), 0);
    bool                folds = false;
    for (size_t i = 0; i < checks.size(); ++i)
    {
      gp_Pnt p;
      gp_Vec su, sv;
      surface->D1(checks[i].flat.X(), checks[i].flat.Y(), p, su, sv);
      const double scale = su.Magnitude() * sv.Magnitude();
      folded[i]          = !(scale > 0) || su.Crossed(sv).Dot(checks[i].normal) <= 0;
      folds              = folds || folded[i];
      if (folded[i])
        continue;
      angles.push_back(std::asin(std::min(1.0, std::abs(su.Dot(sv)) / scale)) * 180 / M_PI);
      least = std::min(least, su.Magnitude() / sv.Magnitude());
      most  = std::max(most, su.Magnitude() / sv.Magnitude());
    }
    if (folds)
    {
      built.measured.within = false;
      built.measured.unkept = true;
      built.measured.worst  = "the conformal surface folds over inside the loop";
      built.measured.excess = std::max(built.measured.excess, unkeptExcess);
    }
    else if (!angles.empty())
    {
      double squares = 0;
      for (double angle : angles)
        squares += angle * angle;
      std::sort(angles.begin(), angles.end());
      built.squareness = {std::sqrt(squares / angles.size()),
                          angles[static_cast<size_t>(0.95 * (angles.size() - 1))],
                          angles.back(),
                          least,
                          most};
    }
    if (built.measured.within)
      if (const auto fault = BRepCheck_Analyzer(built.face).IsValid()
                               ? insideFault(built.face, asked.points.empty() ? rimBend(asked) : 0)
                               : std::optional<std::string>("it is not a valid face"))
      {
        built.measured.within = false;
        built.measured.unkept = true;
        built.measured.worst  = *fault;
        built.measured.excess = unkeptExcess;
      }
    return Checked{std::move(built), std::move(folded), folds};
  };
  std::optional<Built> best;
  std::vector<Target>  bestTargets;
  int                  idle = 0, unkept = 0, iterations = 0;
  for (double before = HUGE_VAL;;)
  {
    ++iterations;
    SplineNet net{degree, degree, knots(breaksU), knots(breaksV), {}};
    const int nu = net.countU(), nv = net.countV();
    double    smallest = HUGE_VAL;
    for (const auto* breaks : {&breaksU, &breaksV})
      for (size_t i = 1; i < breaks->size(); ++i)
        smallest = std::min(smallest, (*breaks)[i] - (*breaks)[i - 1]);
    const auto mean = [&](const std::vector<double>& knots, int i) {
      double sum = 0;
      for (int k = 1; k <= degree; ++k)
        sum += knots[i + k];
      return sum / degree;
    };
    for (int i = 0; i < nu; ++i)
      for (int j = 0; j < nv; ++j)
        net.poles.push_back(frame.origin.XYZ() + frame.x.XYZ() * mean(net.knotsU, i)
                            + frame.y.XYZ() * mean(net.knotsV, j));
    // The Natural fill inside, wished for at the node nearest the middle of every part of
    // every span, `placesPerSpan` parts each way: as many wishes as the net can follow.
    std::vector<Target> targets;
    {
      const auto part = [&](const std::vector<double>& breaks, double at) {
        const int span = std::clamp(
          static_cast<int>(std::upper_bound(breaks.begin(), breaks.end(), at) - breaks.begin()) - 1,
          0,
          static_cast<int>(breaks.size()) - 2);
        const double share = (at - breaks[span]) / (breaks[span + 1] - breaks[span]);
        const int sub = std::clamp(static_cast<int>(share * placesPerSpan), 0, placesPerSpan - 1);
        return std::pair{span * placesPerSpan + sub,
                         breaks[span]
                           + (breaks[span + 1] - breaks[span]) * (sub + 0.5) / placesPerSpan};
      };
      std::unordered_map<long long, std::pair<int, double>> nearest;
      const long long columns = static_cast<long long>(breaksV.size()) * placesPerSpan;
      for (size_t i = 0; i < chart.flat.size(); ++i)
      {
        if (onBoundary[i])
          continue;
        const auto [a, middleU] = part(breaksU, chart.flat[i].X());
        const auto [b, middleV] = part(breaksV, chart.flat[i].Y());
        const double distance   = (chart.flat[i] - gp_XY(middleU, middleV)).SquareModulus();
        auto&        kept =
          nearest.try_emplace(a * columns + b, static_cast<int>(i), HUGE_VAL).first->second;
        if (distance < kept.second)
          kept = {static_cast<int>(i), distance};
      }
      std::vector<int> chosen;
      for (const auto& [cell, kept] : nearest)
        chosen.push_back(kept.first);
      std::sort(chosen.begin(), chosen.end());
      // Each stands for its share of the loop's area, as a length.
      const double share = std::sqrt(area / std::max<size_t>(1, chosen.size())) / asked.size;
      for (int i : chosen)
      {
        targets.push_back({chart.point[i], chart.flat[i], 0, {}, {}, share});
        targets.back().soft = true;
      }
    }
    for (size_t s = 0; s < asked.sides.size(); ++s)
    {
      const Side& side = asked.sides[s];
      if (side.free)
        continue;
      const auto place = layout(s);
      const int  count =
        std::clamp(static_cast<int>(std::ceil(placesPerSpan * side.length / smallest)),
                   fewestRimPlaces,
                   mostRimPlaces);
      const double step = 1.0 / count;
      // Next to a gap the face has its say at the end as along the rest (see overPlane).
      const auto quiet = [&](size_t other) {
        const Side& near = asked.sides[other % asked.sides.size()];
        return near.free || near.support.IsSame(side.support);
      };
      const bool lone       = hasGaps(asked);
      const bool quietStart = lone && quiet(s + asked.sides.size() - 1),
                 quietEnd   = lone && quiet(s + 1);
      auto add              = [&](double f, bool end, bool watch) {
        if (end && ((f == 0 && quietStart) || (f == 1 && quietEnd)))
        {
          f   = std::clamp(f, cornerInset, 1 - cornerInset);
          end = false;
        }
        const bool smooth = side.order > 0 && !end;
        // Follow the Natural fill at its own edge correspondence. Any
        // positional correction below targets the original spatial curve.
        const gp_XY uv = held[s]->Value(side.parameter(f)).XY();
        Target      target{source->Value(uv.X(), uv.Y()),
                      place(side.parameter(f)),
                      smooth ? side.order : 0,
                                        {},
                                        {},
                      step * side.length / asked.size};
        target.watch = watch;
        if (smooth)
        {
          target.jet           = side.jet(f);
          target.easing        = side.easing(f);
          const double h       = std::min({step / 4, f, 1 - f});
          const gp_XY  tangent = place(side.parameter(f + h)) - place(side.parameter(f - h));
          target.across        = gp_XY(-tangent.Y(), tangent.X()) / tangent.Modulus();
        }
        targets.push_back(target);
      };
      add(0, true, false);
      add(1, true, false);
      for (int k = 0; k < count; ++k)
      {
        add((k + 0.5) * step, false, false);
        if (k + 1 < count)
          add((k + 1) * step, false, true);
      }
    }
    for (const Side& guide : asked.guides)
    {
      const int count =
        std::clamp(static_cast<int>(std::ceil(placesPerSpan * guide.length / smallest)),
                   fewestRimPlaces,
                   mostRimPlaces);
      for (int i = 0; i <= 2 * count; ++i)
      {
        const gp_Pnt point = guide.meet(static_cast<double>(i) / (2 * count));
        targets.push_back(
          {point, interiorPlace(point), 0, {}, {}, guide.length / count / asked.size});
        targets.back().watch = i % 2 != 0;
      }
    }
    for (const gp_Pnt& point : asked.points)
      targets.push_back({point, interiorPlace(point), 0, {}, {}, 1});
    occ::handle<Geom_BSplineSurface> surface = fit(
      net,
      std::vector<char>(net.poles.size(), 0),
      targets,
      asked.size,
      asked.limits,
      {},
      [&](double u, double v) { return within(gp_XY(u, v)) ? 1 : outsideWeight; },
      nullptr,
      {},
      &frame);
    // As in the Natural fit, the fair objective can leave a small positional
    // residual. Use the shared geometric correction rather than relaxing the
    // delivery tolerance. Supported derivatives remain the responsibility of
    // the constrained fit; this correction applies only to positional inputs.
    if (order == 0 && (!asked.guides.empty() || !asked.points.empty()))
    {
      const auto onto = [&](double span) {
        std::vector<RimPlace> places;
        for (size_t s = 0; s < asked.sides.size(); ++s)
        {
          const Side& side = asked.sides[s];
          if (side.free)
            continue;
          const auto place = layout(s);
          const int  count = std::clamp(static_cast<int>(rimPlacesPerSpan * side.length / span),
                                       fewestRimPlaces,
                                       mostRimPlaces);
          for (int i = 0; i <= count; ++i)
          {
            const double f = static_cast<double>(i) / count;
            places.push_back({side.meet(f), place(side.parameter(f)), {}, false});
          }
        }
        for (const Side& guide : asked.guides)
        {
          const int count = std::clamp(static_cast<int>(rimPlacesPerSpan * guide.length / span),
                                       fewestRimPlaces,
                                       mostRimPlaces);
          for (int i = 0; i <= count; ++i)
          {
            const gp_Pnt point = guide.meet(static_cast<double>(i) / count);
            places.push_back({point, interiorPlace(point), {}, false});
          }
        }
        for (const gp_Pnt& point : asked.points)
          places.push_back({point, interiorPlace(point), {}, false});
        return places;
      };
      net = netOf(surface);
      // The correction is held to the net's own scale, as over the mean plane: a coarse
      // net that is far off is made finer by the fit, and not by thousands of knots here.
      correctRim(
        op,
        net,
        onto,
        frame.z.XYZ(),
        rimShare * asked.limits.g0,
        enoughShare * asked.limits.g0,
        asked.controls.open
          ? std::min(mostRimSpans,
                     rimFinerMet * static_cast<int>(std::max(breaksU.size(), breaksV.size()) - 1))
          : mostRimSpans);
      surface = surfaceOf(net);
    }
    auto        checked = certify(surface);
    Built&      built   = checked.built;
    const auto& folded  = checked.folded;
    const bool  folds   = checked.folds;
    // The fill's shape is wished for, not asked: a net that meets everything asked is kept
    // and is made finer only where the shape is followed by less than the tolerance, so
    // that the chart is followed and U and V come out square. That is never a refusal.
    // It is asked no more than a fair surface over its net comes to a curve (see
    // fitReach): a finer net beyond that costs more than the flow gains.
    const double wishedWithin =
      asked.controls.open ? std::max(asked.limits.g0, fitReach * asked.size) : asked.limits.g0;
    bool followed = true;
    for (const Target& target : targets)
      followed =
        followed
        && (!target.soft
            || surface->Value(target.uv.X(), target.uv.Y()).Distance(target.point) <= wishedWithin);
    built.shapeTarget = {followed ? 1.0 : 0.0, static_cast<double>(iterations)};
    if (built.measured.within)
    {
      best        = built;
      bestTargets = targets;
      if (followed)
        break;
    }
    else if (!best || (!best->measured.within && built.measured.excess < best->measured.excess))
    {
      best = built;
    }
    // Nets that come no nearer, more than the coarsest few may: what is left is not the net's.
    if (built.measured.within)
    {
    }
    else if (built.measured.excess > stillFalling * before)
    {
      if (++idle > unruly)
      {
        best->measured.worst += ". A finer net came no nearer";
        break;
      }
    }
    else
    {
      idle   = 0;
      before = built.measured.excess;
    }
    // Finer in every span where a place along an edge is missed, the surface folds or the
    // shape wished for is not followed: halved, up to the finest net of the ladder that the
    // Natural fill is given.
    std::vector<char> splitU(breaksU.size() - 1, 0), splitV(breaksV.size() - 1, 0);
    const auto        mark = [&](const gp_XY& uv) {
      const auto span = [](const std::vector<double>& breaks, double at) {
        return std::clamp(
          static_cast<int>(std::upper_bound(breaks.begin(), breaks.end(), at) - breaks.begin()) - 1,
          0,
          static_cast<int>(breaks.size()) - 2);
      };
      splitU[span(breaksU, uv.X())] = splitV[span(breaksV, uv.Y())] = 1;
    };
    for (const Target& target : targets)
      if (surface->Value(target.uv.X(), target.uv.Y()).Distance(target.point)
          > (target.soft ? wishedWithin : fitShare * asked.limits.g0))
        mark(target.uv);
    // Certification reads between fit rows. A miss there must refine the
    // crossed spans too, even when every row used by the solver was met.
    for (size_t s = 0; s < asked.sides.size(); ++s)
    {
      const Side& side = asked.sides[s];
      if (side.free || built.measured.values[side.input * 4] <= asked.limits.g0)
        continue;
      const auto place = layout(s);
      const int  count =
        std::clamp(static_cast<int>(std::ceil(4 * placesPerSpan * side.length / smallest)),
                   fewestRimPlaces,
                   4 * mostRimPlaces);
      for (int i = 0; i <= count; ++i)
        mark(place(side.parameter(static_cast<double>(i) / count)));
    }
    for (size_t i = 0; i < checks.size(); ++i)
      if (folded[i])
        mark(checks[i].flat);
    // A net that meets the edges with a surface that is not one to keep is too coarse for
    // what is asked along them, wherever that shows: a few finer nets, and no more.
    if (built.measured.unkept && !folds && ++unkept <= mostUnkept)
    {
      std::fill(splitU.begin(), splitU.end(), 1);
      std::fill(splitV.begin(), splitV.end(), 1);
    }
    bool finer = false;
    for (auto [breaks, split] : {std::pair{&breaksU, &splitU}, std::pair{&breaksV, &splitV}})
    {
      std::vector<double> next{breaks->front()};
      int                 spans = static_cast<int>(breaks->size()) - 1;
      for (size_t i = 0; i + 1 < breaks->size(); ++i)
      {
        if ((*split)[i] && spans < std::end(ladder)[-1])
        {
          next.push_back(((*breaks)[i] + (*breaks)[i + 1]) / 2);
          ++spans;
          finer = true;
        }
        next.push_back((*breaks)[i + 1]);
      }
      *breaks = next;
    }
    if (!finer)
      break;
  }
  if (order == 0)
    measureSpatialBoundary(*best, asked, false);
  if (order == 0 && best->measured.within && !best->measured.unkept)
  {
    double gap = 0;
    for (const Side& side : asked.sides)
      gap = std::max(gap, best->measured.values[side.input * 4]);
    const double allowance = std::nextafter(asked.limits.g0 - gap, 0.0);
    const auto   original  = occ::down_cast<Geom_BSplineSurface>(BRep_Tool::Surface(best->face));
    if (std::isfinite(gap) && gap >= 0 && !original.IsNull()
        && std::all_of(best->measured.values.begin(),
                       best->measured.values.end(),
                       [](double value) { return std::isfinite(value); }))
      if (const auto compact = compactPatchSurface(op, original, allowance))
      {
        try
        {
          auto reduced = certify(compact->surface, true);
          measureSpatialBoundary(reduced.built, asked, false);
          bool followed = true;
          for (const Target& target : bestTargets)
          {
            if (!target.soft)
              continue;
            const double before =
              original->Value(target.uv.X(), target.uv.Y()).Distance(target.point);
            const double after =
              compact->surface->Value(target.uv.X(), target.uv.Y()).Distance(target.point);
            followed = followed && std::isfinite(before) && std::isfinite(after)
                       && after <= std::max(asked.limits.g0, before);
          }
          const bool finite = std::all_of(reduced.built.squareness.begin(),
                                          reduced.built.squareness.end(),
                                          [](double value) { return std::isfinite(value); });
          if (reduced.built.measured.within && !reduced.built.measured.unkept && followed && finite)
          {
            reduced.built.shapeTarget = best->shapeTarget;
            best                      = std::move(reduced.built);
          }
        }
        catch (const Failure& failure)
        {
          if (failure.reason != no_solution)
            throw;
        }
        catch (const Standard_Failure&)
        {
          // The certified original remains available if its reduced trim cannot be rebuilt.
        }
      }
  }
  best->shapeTarget[1] = static_cast<double>(iterations);
  return *best;
}

//! Next/Adjacent with Quad sided off. Over a conformal chart of a surface, the
//! surface's U and V meet at right angles; so the Natural fill is laid flat conformally, and
//! a new surface is fitted over that chart: one that meets the same edges as asked and keeps
//! near the Natural fill inside, as a wish and not a limit. How square U and V
//! come out is read at every node and in every triangle and reported (see Built::squareness).
Built adjacent(const Asked& asked, const Frame& frame, const Built& natural)
{
  PatchPhase phase("adjacent");
  double     perimeter = 0;
  for (const Side& side : asked.sides)
    perimeter += side.length;
  std::vector<gp_XY> outline;
  for (const gp_Pnt& point : polygon(asked.sides))
    outline.push_back(frame.project(point));
  const double area = std::abs(enclosed(outline));
  // Each edge as the Natural face holds it, in that face's parameters, and so in the chart.
  std::vector<occ::handle<Geom2d_Curve>> held;
  for (const Side& side : asked.sides)
  {
    if (side.free)
    {
      held.push_back(natural.gaps[&side - asked.sides.data()]);
      continue;
    }
    occ::handle<Geom2d_Curve> found;
    for (TopExp_Explorer it(natural.face, TopAbs_EDGE); it.More() && found.IsNull(); it.Next())
    {
      double          first, last;
      TopLoc_Location location;
      const auto      edge = TopoDS::Edge(it.Current());
      if (BRep_Tool::Curve(edge, location, first, last) == side.curve && location.IsIdentity()
          && first == std::min(side.from, side.to) && last == std::max(side.from, side.to))
        found = BRep_Tool::CurveOnSurface(edge, natural.face, first, last);
    }
    if (found.IsNull())
      fail(op,
           no_solution,
           "the Natural fill does not hold boundary edge " + std::to_string(side.input + 1));
    held.push_back(found);
  }
  const ConformalChart chart = conformalChart(natural.face,
                                              area / perimeter / trianglesAcross,
                                              chartDeflection * asked.limits.g0);
  return adjacentOver(asked, frame, natural, chart, held, area);
}

// --- OpenCascade's plate ---------------------------------------------------------------------

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
Built fillFace(Asked asked, int method, double flatness, const PatchControls& controls)
{
  const bool gapped      = hasGaps(asked);
  asked.controls         = controls;
  Asked   input          = asked;
  bool    copyBudgetMiss = false;
  auto    placement      = method != 2 ? localSplineRim(asked, &copyBudgetMiss) : std::nullopt;
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
  agreeAtGuideEnds(asked);

  // A plane in the mean-plane coordinates has orthogonal UV everywhere, so all
  // nonquad flow choices described as Normal can be satisfied exactly there.
  const auto plane = method == 0 && controls.quad != 1 ? planar(asked, frame) : std::nullopt;
  // With Quad sided off, Next and Adjacent both ask for Normal flow: the Natural
  // fill, then the same over its conformal chart (see adjacent). With it on, they are
  // directions across the quad's sides (see foursided).
  const bool conformal = controls.flow >= (gapped ? 1 : 2) && controls.quad == 0 && !plane;
  if (controls.flow >= 2 && controls.quad != 1 && !plane && !conformal)
    fail(op,
         unsupported_geometry,
         "Next and Adjacent UV flow are not implemented for this boundary; select Natural");
  std::optional<Built> built;
  if (plane)
  {
    built = Built{0, *plane, measureAll(BRep_Tool::Surface(*plane), asked)};
  }
  else if (method == 2)
  {
    built = plate(asked);
  }
  else if (controls.quad == 1)
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
  else if (controls.flow > 0 && !conformal)
  {
    fail(op,
         unsupported_geometry,
         "UV flow controls currently require four regular sides and Quad sided");
  }
  else if (flatness >= 0 || controls.quad == 0 || controls.quality > 0)
  {
    built = overPlane(asked, frame);
  }
  else if (const auto face = method == 0 ? planar(asked, frame) : std::nullopt)
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
  if (!built->measured.within && method != 2 && !gapped)
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
    fail(op, no_solution, "no Patch candidate meets what was asked: " + built->measured.worst);
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

  if (method != 2)
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
  return *built;
}

//! The face of the loop `asked` holds, held to every check Patch makes; fails with what it
//! misses otherwise.
Built patchFace(Asked asked, int method, double flatness, PatchControls controls)
{
  if (!arrange(asked.sides, asked.limits.g0, controls.open && method != 2))
    chain(asked.sides, asked.limits.g0);
  const bool gapped = hasGaps(asked);
  // What `open` changes is for a boundary that does not close. One that does is the fill it
  // always was, whoever asks for it.
  controls.open = gapped;
  // An open boundary asked to be a quad, whose runs and gaps are four sides with corners, is the
  // untrimmed quad with the gaps' poles free (see foursided). Where that is not a face to keep,
  // it is filled as any open boundary is: every flow but Natural is then the refit over the
  // conformal chart.
  if (gapped && controls.quad == 1 && quadRuns(asked.sides, asked.limits.g1))
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
    controls.quad = 0;
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
                         int                          method,
                         double                       flatness,
                         const PatchControls&         controls)
{
  PatchPhase phase("total");
  if (edges.empty() || edges.size() > 64 || supports.size() != edges.size()
      || orders.size() != edges.size() || points.size() % 3 != 0 || method < 0 || method > 2)
    fail(op,
         invalid_input,
         "one to 64 boundary edges, each with a support (or 0) and an order, are required");
  if (!(tolerance > 0) || !(g1Degrees > 0) || !(g2Share > 0)
      || !std::isfinite(tolerance + g1Degrees + g2Share))
    fail(op, invalid_input, "the tolerances must be positive and finite");
  if (!std::isfinite(flatness) || (flatness != -1 && (flatness < 0 || flatness > 1))
      || (flatness >= 0 && method == 2))
    fail(op, invalid_input, "Flatness must be in 0..1 and uses the variational solver");
  if (controls.quad < -1 || controls.quad > 1 || controls.flow < 0 || controls.flow > 3
      || controls.quality < 0 || controls.quality > 2)
    fail(op, invalid_input, "invalid parameterization, flow or quality control");
  if (method == 2 && (controls.quad != -1 || controls.flow || controls.quality || !controls.trim))
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
                  int                          method,
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
  return patchFace(asked, 0, -1, PatchControls{}).face;
}

} // namespace BRepExt
