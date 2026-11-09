// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Blend: a certificate that two faces are farther apart than a tolerance.
//
// A face lies on its surface inside the parameter box of its edges' curves on it, and its edges
// and vertices lie within its largest tolerance of that patch. The patch is converted exactly into
// a rational B-spline with positive weights, which lies in the convex hull of its poles (Piegl and
// Tiller, The NURBS Book, 2nd ed. 1997, P4.4 and its tensor-product form). Two hulls separated
// along any direction by more than the tolerance and both growths separate the faces. Pairs of
// patch pieces that are not are halved, the larger first (Johnson and Cohen, "A framework for
// efficient minimum distance computations", ICRA 1998: the bounding-volume hierarchy on control
// hulls).

#include <BRepBlendExt_Separation.hxx>
#include <BRepBlendExt_Internal.hxx>

#include <BRep_Tool.hxx>
#include <BndLib_Add2dCurve.hxx>
#include <Bnd_Box2d.hxx>
#include <Extrema_GenLocateExtSS.hxx>
#include <Geom2dAdaptor_Curve.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <GeomConvert.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_BezierSurface.hxx>
#include <Geom_Circle.hxx>
#include <Geom_ElementarySurface.hxx>
#include <Geom_Ellipse.hxx>
#include <Geom_Line.hxx>
#include <Geom_OffsetSurface.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_SurfaceOfLinearExtrusion.hxx>
#include <Geom_SurfaceOfRevolution.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Standard_Failure.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>

#include <array>
#include <optional>
#include <vector>

namespace BRepExt::blending::detail
{
namespace
{
//! A face's patch as an exact rational B-spline, and how far the face may lie off it.
struct Cover
{
  Handle(Geom_BSplineSurface) surface;
  double                      grow = 0;
};

//! The parameter box of a face: the union of the boxes of its edges' curves on it. Only curves
//! whose box is bounded by their poles or exactly (lines, circles, B-splines, Beziers) are taken.
bool parameterBox(const TopoDS_Face& face, double& u0, double& u1, double& v0, double& v1)
{
  Bnd_Box2d box;
  for (TopExp_Explorer e(face, TopAbs_EDGE); e.More(); e.Next())
  {
    double     first, last;
    const auto curve = BRep_Tool::CurveOnSurface(TopoDS::Edge(e.Current()), face, first, last);
    if (curve.IsNull())
      return false;
    switch (Geom2dAdaptor_Curve(curve).GetType())
    {
      case GeomAbs_Line:
      case GeomAbs_Circle:
      case GeomAbs_BezierCurve:
      case GeomAbs_BSplineCurve:
        BndLib_Add2dCurve::Add(curve, first, last, 0., box);
        break;
      default:
        return false;
    }
  }
  if (box.IsVoid())
    return false;
  box.Get(u0, v0, u1, v1);
  return std::isfinite(u0) && std::isfinite(u1) && std::isfinite(v0) && std::isfinite(v1)
         && u0 < u1 && v0 < v1;
}

//! Whether GeomConvert turns a curve into a B-spline exactly (its conic and spline branches).
bool convertsExactly(Handle(Geom_Curve) curve)
{
  while (const auto trimmed = Handle(Geom_TrimmedCurve)::DownCast(curve))
    curve = trimmed->BasisCurve();
  return curve->IsKind(STANDARD_TYPE(Geom_Line)) || curve->IsKind(STANDARD_TYPE(Geom_Circle))
         || curve->IsKind(STANDARD_TYPE(Geom_Ellipse))
         || curve->IsKind(STANDARD_TYPE(Geom_BezierCurve))
         || curve->IsKind(STANDARD_TYPE(Geom_BSplineCurve));
}

//! The face's patch converted exactly, or nothing where no exact conversion is known. An offset
//! surface's points lie within its offset of its basis, so the basis is converted and the offset
//! added to the growth.
std::optional<Cover> coverOf(const TopoDS_Face& face)
{
  Cover        cover;
  double       u0, u1, v0, v1;
  Handle(Geom_Surface) surface = BRep_Tool::Surface(face);
  if (surface.IsNull() || !parameterBox(face, u0, u1, v0, v1))
    return std::nullopt;
  cover.grow = tol::largest(face);
  for (;;)
  {
    if (const auto trimmed = Handle(Geom_RectangularTrimmedSurface)::DownCast(surface))
      surface = trimmed->BasisSurface();
    else if (const auto offset = Handle(Geom_OffsetSurface)::DownCast(surface))
    {
      cover.grow += std::abs(offset->Offset());
      surface = offset->BasisSurface();
    }
    else
      break;
  }
  if (const auto r = Handle(Geom_SurfaceOfRevolution)::DownCast(surface))
  {
    if (!convertsExactly(r->BasisCurve()))
      return std::nullopt;
  }
  else if (const auto e = Handle(Geom_SurfaceOfLinearExtrusion)::DownCast(surface))
  {
    if (!convertsExactly(e->BasisCurve()))
      return std::nullopt;
  }
  else if (!surface->IsKind(STANDARD_TYPE(Geom_ElementarySurface))
           && !surface->IsKind(STANDARD_TYPE(Geom_BSplineSurface))
           && !surface->IsKind(STANDARD_TYPE(Geom_BezierSurface)))
    return std::nullopt;
  // A periodic direction is covered by one period; a bounded one only inside its bounds, since an
  // extrapolated patch need not lie in any hull of the surface's poles.
  double su0, su1, sv0, sv1;
  surface->Bounds(su0, su1, sv0, sv1);
  if (surface->IsUPeriodic())
    u1 = std::min(u1, u0 + surface->UPeriod());
  else if (u0 < su0 || u1 > su1)
    return std::nullopt;
  if (surface->IsVPeriodic())
    v1 = std::min(v1, v0 + surface->VPeriod());
  else if (v0 < sv0 || v1 > sv1)
    return std::nullopt;
  try
  {
    cover.surface = GeomConvert::SurfaceToBSplineSurface(
      new Geom_RectangularTrimmedSurface(surface, u0, u1, v0, v1));
  }
  catch (const Standard_Failure&)
  {
    return std::nullopt;
  }
  if (cover.surface.IsNull())
    return std::nullopt;
  if (const auto* weights = cover.surface->Weights())
    for (int i = weights->LowerRow(); i <= weights->UpperRow(); ++i)
      for (int j = weights->LowerCol(); j <= weights->UpperCol(); ++j)
        if (!((*weights)(i, j) > 0))
          return std::nullopt;
  return cover;
}

//! A rectangle of a cover's parameters, the poles of the cover's segment over it and their box.
struct Cell
{
  double                       u0, u1, v0, v1;
  Handle(Geom_BSplineSurface)  piece;
  std::array<double, 6>        box; // x, y, z low, then x, y, z high
  double                       diagonal, magnitude;
  gp_Pnt                       middle;
  gp_Vec                       normal;
};

std::optional<Cell> cellOf(const Handle(Geom_BSplineSurface)& surface,
                           double                             u0,
                           double                             u1,
                           double                             v0,
                           double                             v1)
{
  Cell cell{u0, u1, v0, v1, Handle(Geom_BSplineSurface)::DownCast(surface->Copy()), {}, 0, 0, {}, {}};
  double su0, su1, sv0, sv1;
  surface->Bounds(su0, su1, sv0, sv1);
  try
  {
    // Always from the whole cover, so that the poles carry the rounding of one segmentation.
    if (u0 != su0 || u1 != su1 || v0 != sv0 || v1 != sv1)
      cell.piece->Segment(u0, u1, v0, v1);
  }
  catch (const Standard_Failure&)
  {
    return std::nullopt;
  }
  const auto& poles = cell.piece->Poles();
  cell.box          = {poles(poles.LowerRow(), poles.LowerCol()).X(),
                       poles(poles.LowerRow(), poles.LowerCol()).Y(),
                       poles(poles.LowerRow(), poles.LowerCol()).Z(),
                       0,
                       0,
                       0};
  cell.box[3] = cell.box[0], cell.box[4] = cell.box[1], cell.box[5] = cell.box[2];
  for (int i = poles.LowerRow(); i <= poles.UpperRow(); ++i)
    for (int j = poles.LowerCol(); j <= poles.UpperCol(); ++j)
      for (int k = 0; k < 3; ++k)
      {
        const double c = poles(i, j).Coord(k + 1);
        cell.box[k]    = std::min(cell.box[k], c);
        cell.box[k + 3] = std::max(cell.box[k + 3], c);
      }
  for (int k = 0; k < 3; ++k)
  {
    cell.diagonal += (cell.box[k + 3] - cell.box[k]) * (cell.box[k + 3] - cell.box[k]);
    cell.magnitude = std::max({cell.magnitude, std::abs(cell.box[k]), std::abs(cell.box[k + 3])});
  }
  cell.diagonal = std::sqrt(cell.diagonal);
  gp_Vec du, dv;
  surface->D1((u0 + u1) / 2, (v0 + v1) / 2, cell.middle, du, dv);
  cell.normal = du.Crossed(dv);
  if (!std::isfinite(cell.diagonal) || !std::isfinite(cell.magnitude))
    return std::nullopt;
  return cell;
}

//! The interval of the cell's poles along a direction.
std::pair<double, double> along(const Cell& cell, const gp_XYZ& direction)
{
  const auto& poles = cell.piece->Poles();
  double      low = Precision::Infinite(), high = -Precision::Infinite();
  for (int i = poles.LowerRow(); i <= poles.UpperRow(); ++i)
    for (int j = poles.LowerCol(); j <= poles.UpperCol(); ++j)
    {
      const double d = poles(i, j).XYZ().Dot(direction);
      low            = std::min(low, d);
      high           = std::max(high, d);
    }
  return {low, high};
}

//! A lower bound of the distance between two cells' hulls: their boxes' distance, and the gaps
//! of their poles along the line between their middles and along each one's normal there.
double hullGap(const Cell& x, const Cell& y)
{
  double square = 0;
  for (int k = 0; k < 3; ++k)
  {
    const double gap = std::max(x.box[k] - y.box[k + 3], y.box[k] - x.box[k + 3]);
    if (gap > 0)
      square += gap * gap;
  }
  double best = std::sqrt(square);
  for (gp_Vec axis : {gp_Vec(x.middle, y.middle), x.normal, y.normal})
  {
    const double length = axis.Magnitude();
    if (!(length > 0) || !std::isfinite(length))
      continue;
    const gp_XYZ direction = axis.XYZ() / length;
    const auto [xl, xh]    = along(x, direction);
    const auto [yl, yh]    = along(y, direction);
    best                   = std::max(best, std::max(yl - xh, xl - yh));
  }
  return best;
}

//! The halves of a cell, cut across its longer control polygon at the middle parameter.
std::optional<std::pair<Cell, Cell>> halves(const Handle(Geom_BSplineSurface)& surface,
                                            const Cell&                        cell)
{
  const auto& poles = cell.piece->Poles();
  const int   row = (poles.LowerRow() + poles.UpperRow()) / 2,
            col   = (poles.LowerCol() + poles.UpperCol()) / 2;
  double lu = 0, lv = 0;
  for (int i = poles.LowerRow(); i < poles.UpperRow(); ++i)
    lu += poles(i, col).Distance(poles(i + 1, col));
  for (int j = poles.LowerCol(); j < poles.UpperCol(); ++j)
    lv += poles(row, j).Distance(poles(row, j + 1));
  std::optional<Cell> first, second;
  if (lu >= lv)
  {
    const double m = (cell.u0 + cell.u1) / 2;
    first          = cellOf(surface, cell.u0, m, cell.v0, cell.v1);
    second         = cellOf(surface, m, cell.u1, cell.v0, cell.v1);
  }
  else
  {
    const double m = (cell.v0 + cell.v1) / 2;
    first          = cellOf(surface, cell.u0, cell.u1, cell.v0, m);
    second         = cellOf(surface, cell.u0, cell.u1, m, cell.v1);
  }
  if (!first || !second)
    return std::nullopt;
  return std::make_pair(std::move(*first), std::move(*second));
}
} // namespace

bool separated(const TopoDS_Face& a, const TopoDS_Face& b, double tolerance)
{
  const auto ca = coverOf(a), cb = coverOf(b);
  if (!ca || !cb)
    return false;
  const auto whole = [](const Handle(Geom_BSplineSurface)& s) {
    double u0, u1, v0, v1;
    s->Bounds(u0, u1, v0, v1);
    return cellOf(s, u0, u1, v0, v1);
  };
  auto x = whole(ca->surface), y = whole(cb->surface);
  if (!x || !y)
    return false;
  const GeomAdaptor_Surface sa(ca->surface), sb(cb->surface);
  std::vector<std::pair<Cell, Cell>> pending;
  pending.emplace_back(std::move(*x), std::move(*y));
  while (!pending.empty())
  {
    auto [p, q] = std::move(pending.back());
    pending.pop_back();
    // Each pole's projection is rounded within its coordinates' rounding.
    const double gap =
      hullGap(p, q) - ca->grow - cb->grow - tol::rounding(p.magnitude + q.magnitude);
    if (gap > tolerance)
      continue;
    // Two cover points within the tolerance: no hull of these covers can separate them.
    if (p.middle.Distance(q.middle) <= tolerance)
      return false;
    if (!(gap > 0))
    {
      Extrema_GenLocateExtSS near(sa,
                                  sb,
                                  (p.u0 + p.u1) / 2,
                                  (p.v0 + p.v1) / 2,
                                  (q.u0 + q.u1) / 2,
                                  (q.v0 + q.v1) / 2,
                                  tol::uv,
                                  tol::uv);
      if (near.IsDone() && near.SquareDistance() <= tolerance * tolerance)
        return false;
    }
    // Cells within the tolerance across: their points come within three tolerances, which no
    // finer cells can resolve.
    if (std::max(p.diagonal, q.diagonal) <= tolerance)
      return false;
    const bool split = p.diagonal >= q.diagonal;
    auto       parts = halves(split ? ca->surface : cb->surface, split ? p : q);
    if (!parts)
      return false;
    if (split)
    {
      pending.emplace_back(std::move(parts->first), q);
      pending.emplace_back(std::move(parts->second), std::move(q));
    }
    else
    {
      pending.emplace_back(p, std::move(parts->first));
      pending.emplace_back(std::move(p), std::move(parts->second));
    }
  }
  return true;
}
} // namespace BRepExt::blending::detail
