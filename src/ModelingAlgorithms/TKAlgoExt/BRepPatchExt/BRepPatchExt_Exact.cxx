// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Exact fills: a plane, and the continuations of a torus, a ball, a plane/round/plane strip and
// a spline face the loop already lies on.

#include <BRepPatchExt_Fill.hxx>
#include <GeomAlgoExt_ModelTolerance.hxx>

#include <GeomAlgoExt_Cancel.hxx>

#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <Geom2dConvert.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomConvert.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Plane.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <gp_Pln.hxx>
#include <math_Gauss.hxx>
#include <math_Matrix.hxx>
#include <math_Vector.hxx>

namespace BRepExt::patch_fill
{

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

} // namespace BRepExt::patch_fill
