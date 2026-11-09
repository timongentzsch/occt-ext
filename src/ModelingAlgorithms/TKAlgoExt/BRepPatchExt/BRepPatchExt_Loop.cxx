// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// The loop seen flat: the plane it lies nearest to, the polygon that stands for it, and whether
// points and segments lie inside that polygon.

#include <BRepPatchExt_Fill.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>

#include <GeomAlgoExt_Cancel.hxx>

#include <limits>

namespace BRepExt::patch_fill
{
namespace
{

//! The corners that a side is given in the polygon that stands for the loop.
constexpr int polygonCorners = 32;

//! Below this share of its size, the difference between the longer and the shorter axis of
//! the area a loop encloses says nothing of a direction: the 32 chords that stand for a side
//! are that far from its curve.
constexpr double roundLoop = 1e-3;

//! Two corners are as far from the middle when they differ by less than this share: the
//! rounding of the sums they are measured from, with room.
constexpr double sameCorner = BRepExt::tolerance::exactLength;

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

} // namespace

std::vector<gp_Pnt> polygon(const std::vector<Side>& sides)
{
  std::vector<gp_Pnt> points;
  for (const Side& side : sides)
    for (int i = 0; i < polygonCorners; ++i)
      points.push_back(side.at(static_cast<double>(i) / polygonCorners));
  return points;
}

Frame meanPlane(const std::vector<gp_Pnt>& points, double perimeter)
{
  gp_XYZ centre(0, 0, 0), area(0, 0, 0);
  for (const gp_Pnt& point : points)
    centre += point.XYZ() / static_cast<double>(points.size());
  for (size_t i = 0; i < points.size(); ++i)
    area += (points[i].XYZ() - centre).Crossed(points[(i + 1) % points.size()].XYZ() - centre) / 2;
  // A small aspect ratio is not degeneracy. Only reject an area unresolved
  // at the arithmetic scale of the centered polygon's accumulated cross products.
  const double roundoff =
    8 * points.size() * std::numeric_limits<double>::epsilon() * perimeter * perimeter;
  if (!(area.Modulus() > roundoff))
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

bool near(const gp_XY& from, const gp_XY& to, const gp_XY& point, double slack)
{
  const gp_XY edge = to - from, away = point - from;
  return (away - edge * std::clamp(away.Dot(edge) / edge.SquareModulus(), 0.0, 1.0)).Modulus()
         <= slack;
}

std::optional<double> crossing(const gp_XY& from, const gp_XY& to, double y)
{
  if ((to.Y() > y) == (from.Y() > y))
    return std::nullopt;
  return from.X() + (to.X() - from.X()) * (y - from.Y()) / (to.Y() - from.Y());
}

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

} // namespace BRepExt::patch_fill
