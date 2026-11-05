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

// Split blend pieces: inherited support breaks, compatible section frames and joined sheets.
#include <BRepBlendExt_Internal.hxx>
#include <GeomAlgoExt_Cancel.hxx>
#include <GeomAlgoExt_Hermite.hxx>
#include <set>

namespace BRepExt::blending::detail
{
// Express the sections in a common basis before comparing or joining their poles.
static bool commonSections(const std::vector<Handle(Geom_BSplineSurface)>& surfaces)
{
  if (surfaces.empty())
    return false;
  std::map<double, int> knots;
  const auto&           first = surfaces.front();
  for (const auto& surface : surfaces)
  {
    if (surface->VDegree() != first->VDegree() || surface->IsVPeriodic()
        || std::abs(surface->VKnot(1)) > tol::parameter
        || std::abs(surface->VKnot(surface->NbVKnots()) - 1) > tol::parameter)
      return false;
    for (int i = 1; i <= surface->NbVKnots(); ++i)
      knots[surface->VKnot(i)] = std::max(knots[surface->VKnot(i)], surface->VMultiplicity(i));
  }
  for (const auto& surface : surfaces)
    for (const auto& [knot, multiplicity] : knots)
      surface->InsertVKnot(knot, multiplicity, tol::parameter, false);
  return true;
}

// A compact C2 displacement with unit longitudinal derivative at the join.
// Its physical reach is independent of the tessellation/fitting knot spacing.
struct CouplingWindow
{
  double              reach = 0, width = 0;
  std::vector<double> derivative, elevation;
};

static CouplingWindow couplingWindow(const Handle(Geom_BSplineSurface)& surface,
                                     const Rail&                        rail,
                                     bool                               last,
                                     double                             requested = 0,
                                     bool                               lift      = false)
{
  CouplingWindow result;
  const double   length = railLength(rail, rail.from, rail.to);
  const double   radius = std::max(std::abs(rail.s[0].distance), std::abs(rail.s[1].distance));
  if (!std::isfinite(radius) || !std::isfinite(length) || !std::isfinite(requested))
    return result;
  const double reach = std::min(requested > 0 ? requested : radius, length / 2);
  if (!(reach > limit::degenerate * rail.scale) || surface->UDegree() < 5)
    return result;
  const auto at = railUAtDistance(rail, last ? rail.to : rail.from, last ? -reach : reach);
  if (!at)
    return result;
  const double boundary = *at, width = last ? 1 - boundary : boundary;
  if (!(boundary > tol::parameter && boundary < 1 - tol::parameter) || !std::isfinite(width))
    return result;
  surface->InsertUKnot(boundary, surface->UDegree() - 2, tol::parameter, false);
  std::vector<double> flat;
  for (int i = 1; i <= surface->NbUKnots(); ++i)
    flat.insert(flat.end(), surface->UMultiplicity(i), surface->UKnot(i));
  const auto sites = GeomAlgoExt::hermiteSites(surface->UDegree(), flat, 2, 2);
  if (sites.empty())
    return result;
  std::vector<gp_Vec> derivative(sites.size()), elevation(sites.size());
  for (size_t i = 0; i < sites.size(); ++i)
  {
    const auto& site = sites[i];
    double      h = 0, e = 0;
    if (site.order == 0)
    {
      const double x = last ? (site.parameter - boundary) / width : 1 - site.parameter / width;
      if (x > 0 && x < 1)
      {
        h = (last ? 1 : -1) * width * x * x * x * (-4 + x * (7 - 3 * x));
        e = x * x * x * (10 + x * (-15 + 6 * x));
      }
      else if (x >= 1)
        e = 1;
    }
    else if (site.order == 1 && std::abs(site.parameter - (last ? 1 : 0)) <= tol::parameter)
      h = 1;
    derivative[i] = gp_Vec(h, 0, 0);
    elevation[i]  = gp_Vec(e, 0, 0);
  }
  if (!GeomAlgoExt::interpolateHermite(surface->UDegree(), flat, sites, derivative)
      || (lift && !GeomAlgoExt::interpolateHermite(surface->UDegree(), flat, sites, elevation)))
    return result;
  result.reach = reach;
  result.width = width;
  for (const auto& v : derivative)
    result.derivative.push_back(v.X());
  if (lift)
    for (const auto& v : elevation)
      result.elevation.push_back(v.X());
  return result;
}

// A regularized section law must retain the sheet orientation and immersion.
// Its acceptance budget permits at most one ball-curvature unit of added shape
// operator (Frobenius norm, hence also an operator-norm bound). This bounds the
// correction, not the total support curvature or equivalence to a native law.
static bool regularizedSheet(const char*                        op,
                             const Handle(Geom_BSplineSurface)& moved,
                             const Handle(Geom_BSplineSurface)& original,
                             const Rail&                        rail,
                             double                             lawGap)
{
  const double radius = std::max(std::abs(rail.s[0].distance), std::abs(rail.s[1].distance));
  if (!std::isfinite(radius) || !(radius > limit::degenerate * rail.scale))
    return false;
  const gp_Vec axes[3] = {gp_Vec(1, 0, 0), gp_Vec(0, 1, 0), gp_Vec(0, 0, 1)};
  for (int i = 2; i <= moved->NbUKnots(); ++i)
  {
    checkCancelled(op);
    for (int j = 2; j <= moved->NbVKnots(); ++j)
      for (double qu : {0.211324865405187, 0.788675134594813})
        for (double qv : {0.211324865405187, 0.788675134594813})
        {
          const double u = moved->UKnot(i - 1) + (moved->UKnot(i) - moved->UKnot(i - 1)) * qu,
                       v = moved->VKnot(j - 1) + (moved->VKnot(j) - moved->VKnot(j - 1)) * qv;
          SurfaceJet a{op}, b{op};
          moved->D2(u, v, a.p, a.u, a.v, a.uu, a.vv, a.uv);
          original->D2(u, v, b.p, b.u, b.v, b.uu, b.vv, b.uv);
          a.n                 = a.u.Crossed(a.v);
          b.n                 = b.u.Crossed(b.v);
          const double metric = a.u.Magnitude() * a.v.Magnitude();
          if (!std::isfinite(metric) || !std::isfinite(a.n.SquareMagnitude())
              || !std::isfinite(b.n.SquareMagnitude()) || !(metric > 0)
              || !(a.n.Magnitude() > tol::coarseAngle * metric) || !(a.n.Dot(b.n) > 0))
            return false;
          a.n.Normalize();
          b.n.Normalize();
          // The sheet stays on its section law (the fit left it within fitTolerance of it). Where
          // a support's curvature jumps the law is two, `lawGap` apart at the join, and holding
          // G1 there moves the sheet: by no more than they differ (measured at most 0.88 of it

          if (!(std::abs(gp_Vec(b.p, a.p).Dot(b.n)) <= limit::section * radius + lawGap))
            return false;
          double change = 0;
          for (int x = 0; x < 3; ++x)
            for (int y = x; y < 3; ++y)
            {
              const double d = a.curvature(axes[x], axes[y]) - b.curvature(axes[x], axes[y]);
              change += (x == y ? 1 : 2) * d * d;
            }
          if (!(change * radius * radius <= 1))
            return false;
        }
  }
  return true;
}

// Couple the transverse derivatives along the shared section. The inherited
// support seam may carry a curvature jump, but the blend must have one tangent
// plane all the way across it. The original shared-row mismatch must fit the
// geometric error budget; the subsequent correction fixes endpoint jets so G2
// contact to each support and neighboring chain joints are preserved.
static bool coupleSplitNormals(const char*                  op,
                               Handle(Geom_BSplineSurface)& a,
                               Handle(Geom_BSplineSurface)& b,
                               double                       tolerance,
                               const Rail&                  railA,
                               const Rail&                  railB,
                               bool                         localLaw,
                               const gp_Pnt2d&              supportUV,
                               double                       lawGap,
                               double                       certified = tol::angularResolution)
{
  if (a->IsURational() || a->IsVRational() || b->IsURational() || b->IsVRational()
      || !commonSections({a, b}))
  {
    return false;
  }
  // Match the shared section before fitting derivatives: replacing it afterward
  // would invalidate the tangent-plane solve when the U span is short.
  for (int j = 1; j <= a->NbVPoles(); ++j)
  {
    const double gap = a->Pole(a->NbUPoles(), j).Distance(b->Pole(1, j));
    if (!(gap <= tolerance))
    {
      return false;
    }
    const gp_Pnt shared((a->Pole(a->NbUPoles(), j).XYZ() + b->Pole(1, j).XYZ()) / 2);
    a->SetPole(a->NbUPoles(), j, shared);
    b->SetPole(1, j, shared);
  }
  bool compatible = true;
  for (int span = 2; span <= a->NbVKnots(); ++span)
    for (int sample = 1; sample <= 16; ++sample)
    {
      const double v = a->VKnot(span - 1) + (a->VKnot(span) - a->VKnot(span - 1)) * sample / 17;
      gp_Pnt       pa, pb;
      gp_Vec       au, av, bu, bv;
      a->D1(1, v, pa, au, av);
      b->D1(0, v, pb, bu, bv);
      compatible = compatible && au.Crossed(av).Angle(bu.Crossed(bv)) <= tol::angularResolution;
    }
  if (compatible)
    return true;
  const auto          originalA    = Handle(Geom_BSplineSurface)::DownCast(a->Copy());
  const auto          originalB    = Handle(Geom_BSplineSurface)::DownCast(b->Copy());
  constexpr int       initialSpans = 16, maxRefinements = 1024, maxIterations = 20;
  std::vector<double> refinements;
  if (!localLaw)
    for (int i = 1; i < initialSpans; ++i)
      refinements.push_back(double(i) / initialSpans);
  for (int iteration = 0; iteration < maxIterations; ++iteration)
  {
    checkCancelled(op);
    a = Handle(Geom_BSplineSurface)::DownCast(originalA->Copy());
    b = Handle(Geom_BSplineSurface)::DownCast(originalB->Copy());
    for (double v : refinements)
    {
      a->InsertVKnot(v, 1, tol::parameter, false);
      b->InsertVKnot(v, 1, tol::parameter, false);
    }
    if (localLaw)
    {
      // Keep exact contact jets through order two. Only their adjacent
      // basis-support intervals shrink when the angular certificate asks
      // for more resolution; the interior derivative field stays exact.
      const double boundary = std::ldexp(.25, -iteration);
      a->InsertVKnot(boundary, a->VDegree() - 2, tol::parameter, false);
      a->InsertVKnot(1 - boundary, a->VDegree() - 2, tol::parameter, false);
      b->InsertVKnot(boundary, b->VDegree() - 2, tol::parameter, false);
      b->InsertVKnot(1 - boundary, b->VDegree() - 2, tol::parameter, false);
    }
    const int           count = a->NbVPoles(), degree = a->VDegree();
    std::vector<double> flat;
    for (int i = 1; i <= a->NbVKnots(); ++i)
      for (int j = 0; j < a->VMultiplicity(i); ++j)
        flat.push_back(a->VKnot(i));
    const auto sites = GeomAlgoExt::hermiteSites(degree, flat, 2, 2);
    if (sites.empty())
      return false;
    std::vector<gp_Vec> correctionA(count), correctionB(count);
    double              inwardSign = 1;
    if (localLaw)
    {
      gp_Pnt p;
      gp_Vec du, dv;
      originalA->D1(1, 0, p, du, dv);
      Jet support;
      if (!std::isfinite(supportUV.X()) || !std::isfinite(supportUV.Y())
          || !jetAt(railA.s[0], supportUV.X(), supportUV.Y(), support)
          || !(p.Distance(support.p) <= fitTolerance(railA)))
      {
        return false;
      }
      const double orientation = du.Crossed(dv).Dot(support.n);
      if (!std::isfinite(orientation) || orientation == 0)
        return false;
      inwardSign = orientation > 0 ? -1 : 1;
    }
    double rateA = 1, rateB = 1;
    if (localLaw)
    {
      gp_Pnt ratePoint;
      gp_Vec railSpeedA, railSpeedB;
      railA.curve.D1(railA.to, ratePoint, railSpeedA);
      railB.curve.D1(railB.from, ratePoint, railSpeedB);
      rateA = railSpeedA.Magnitude() * std::abs(railA.to - railA.from);
      rateB = railSpeedB.Magnitude() * std::abs(railB.to - railB.from);
      if (!(std::isfinite(rateA) && std::isfinite(rateB) && rateA > 0 && rateB > 0))
        return false;
    }
    if (localLaw)
    {
      // Arithmetic mean in physical rail arclength. These are exact
      // spline coefficient operations, including tangential components.
      // Away from preserved endpoint basis-support intervals, both
      // corrected derivatives are proportional. Their tangent planes then
      // agree for any common displacement of the section. The angular
      // certificate below also covers those contact boundary intervals.
      const double factorA = a->UDegree() / (a->UKnot(a->NbUKnots()) - a->UKnot(a->NbUKnots() - 1)),
                   factorB = b->UDegree() / (b->UKnot(2) - b->UKnot(1));
      for (int j = 3; j < count - 3; ++j)
      {
        const gp_Vec da =
          gp_Vec(a->Pole(a->NbUPoles() - 1, j + 1), a->Pole(a->NbUPoles(), j + 1)) * factorA;
        const gp_Vec db = gp_Vec(b->Pole(1, j + 1), b->Pole(2, j + 1)) * factorB;
        correctionA[j]  = (db * (rateA / rateB) - da) * .5;
        correctionB[j]  = correctionA[j] * (-rateB / rateA);
      }
    }
    else
    {
      for (int i = 3; i < count - 3; ++i)
      {
        const double v = sites[i].parameter;
        gp_Pnt       pa, pb;
        gp_Vec       au, av, bu, bv;
        originalA->D1(1, v, pa, au, av);
        originalB->D1(0, v, pb, bu, bv);
        gp_Vec n = au.Crossed(av).Normalized() + bu.Crossed(bv).Normalized();
        if (!(n.Magnitude() > tol::angle))
          return false;
        n.Normalize();
        correctionA[i] = n * (-n.Dot(au));
        correctionB[i] = n * (-n.Dot(bu));
      }
      if (!GeomAlgoExt::interpolateHermite(degree, flat, sites, correctionA)
          || !GeomAlgoExt::interpolateHermite(degree, flat, sites, correctionB))
        return false;
    }
    NCollection_Array1<double> fieldKnots(1, int(flat.size()));
    for (size_t i = 0; i < flat.size(); ++i)
      fieldKnots(int(i) + 1) = flat[i];
    const auto fieldAt = [&](const std::vector<gp_Vec>& field, double v, int order) {
      math_Matrix basis(1, order + 1, 1, degree + 1);
      int         first = 0;
      BSplCLib::EvalBsplineBasis(order, degree + 1, fieldKnots, v, first, basis);
      gp_Vec value;
      for (int k = 1; k <= degree + 1; ++k)
      {
        const int j = first + k - 2;
        if (j >= 0 && j < count)
          value += field[j] * basis(order + 1, k);
      }
      return value;
    };
    double              reach = 0, required = 0;
    std::vector<gp_Vec> release(count);
    if (localLaw)
    {
      double slope = 0;
      for (int span = 2; span <= a->NbVKnots(); ++span)
        for (int sample = 0; sample <= 16; ++sample)
        {
          const double v = a->VKnot(span - 1) + (a->VKnot(span) - a->VKnot(span - 1)) * sample / 16;
          gp_Pnt       p;
          gp_Vec       au, av, bu, bv;
          originalA->D1(1, v, p, au, av);
          originalB->D1(0, v, p, bu, bv);
          const double sa =
            std::abs(au.Crossed(av).Normalized().Dot(fieldAt(correctionA, v, 0))) / rateA;
          const double sb =
            std::abs(bu.Crossed(bv).Normalized().Dot(fieldAt(correctionB, v, 0))) / rateB;
          if (!std::isfinite(sa) || !std::isfinite(sb))
            return false;
          slope = std::max({slope, sa, sb});
        }
      const double radius = std::max({std::abs(railA.s[0].distance),
                                      std::abs(railA.s[1].distance),
                                      std::abs(railB.s[0].distance),
                                      std::abs(railB.s[1].distance)});
      // |q''| <= 3.940234, |e''| <= 5.773503, |q|/e <= .4.
      // Only normal derivative deficit enters this longitudinal geometric
      // curvature initializer. Tangential corrections remain in the actual
      // surface and in its displacement hull and full world certificate,
      // which includes mixed/transverse curvature and the changed metric.
      const double least = (3.9402339529696953 + .4 * 5.773502691896258) * slope * radius;
      // That is the least reach; the curvature it adds falls as the reach grows, so the
      // widest the sections can bear: the displacement, at most (16/81 + .4) slope reach by the
      // bounds below, within what two constructions share, and no wider than the ball.
      // (A deficit of 1e-7 would otherwise leave a window of 1e-7 mm, a span nothing can resolve.)
      reach =
        std::max(least,
                 slope > 0 ? std::min(radius, limit::together * radius / ((16.0 / 81 + .4) * slope))
                           : radius);
      // At the least reach the sections cannot hold within what two constructions share: the
      // law's own kink asks for that move, by the same bounds.
      if (reach == least)
        required = (16.0 / 81 + .4) * slope * reach;
      if (!(std::isfinite(reach) && reach > 0))
        return false;
    }
    const auto windowA = couplingWindow(a, railA, true, reach, localLaw);
    const auto windowB = couplingWindow(b, railB, false, reach, localLaw);
    if (windowA.derivative.empty() || windowB.derivative.empty())
    {
      return false;
    }
    if (localLaw)
    {
      const double widthA = windowA.width, widthB = windowB.width;
      for (int j = 3; j < count - 3; ++j)
      {
        const double v = sites[j].parameter;
        gp_Pnt       p;
        gp_Vec       da, dv, db, vb;
        originalA->D1(1, v, p, da, dv);
        originalB->D1(0, v, p, db, vb);
        const gp_Vec n = (da.Crossed(dv).Normalized() + db.Crossed(vb).Normalized()).Normalized(),
                     m = n * inwardSign;
        release[j]     = m
                     * (.4
                        * std::max({0.0,
                                    widthA * fieldAt(correctionA, v, 0).Dot(m),
                                    -widthB * fieldAt(correctionB, v, 0).Dot(m)}));
      }
      if (!GeomAlgoExt::interpolateHermite(degree, flat, sites, release))
        return false;
      // Convex-hull bound for the full parameter-matched displacement:
      // max|q| = 16/81 and 0 <= e <= 1. This is not a normal-distance
      // or material-containment certificate.
      for (int j = 0; j < count; ++j)
      {
        const double da =
          (16.0 / 81) * windowA.width * correctionA[j].Magnitude() + release[j].Magnitude();
        const double db =
          (16.0 / 81) * windowB.width * correctionB[j].Magnitude() + release[j].Magnitude();
        if (!std::isfinite(da) || !std::isfinite(db))
          return false;
      }
    }

    for (int j = 1; j <= count; ++j)
    {
      const gp_Vec da = correctionA[j - 1], db = correctionB[j - 1];
      for (int i = 1; i <= a->NbUPoles(); ++i)
        a->SetPole(i,
                   j,
                   a->Pole(i, j).Translated(
                     da * windowA.derivative[i - 1]
                     + (localLaw ? release[j - 1] * windowA.elevation[i - 1] : gp_Vec())));
      for (int i = 1; i <= b->NbUPoles(); ++i)
        b->SetPole(i,
                   j,
                   b->Pole(i, j).Translated(
                     db * windowB.derivative[i - 1]
                     + (localLaw ? release[j - 1] * windowB.elevation[i - 1] : gp_Vec())));
    }

    std::vector<double> split;
    for (int span = 2; span <= a->NbVKnots(); ++span)
    {
      bool bad = false;
      for (int i = 1; i < 17; ++i)
      {
        double v = a->VKnot(span - 1) + (a->VKnot(span) - a->VKnot(span - 1)) * double(i) / 17;
        gp_Pnt pa, pb;
        gp_Vec au, av, bu, bv;
        a->D1(1, v, pa, au, av);
        b->D1(0, v, pb, bu, bv);
        const double angle = au.Crossed(av).Angle(bu.Crossed(bv));

        bad = bad || !(angle <= certified);
      }
      if (bad)
        split.push_back((a->VKnot(span - 1) + a->VKnot(span)) / 2);
    }
    if (split.empty())
      return regularizedSheet(op, a, originalA, railA, lawGap + required)
             && regularizedSheet(op, b, originalB, railB, lawGap + required);
    if (localLaw)
      continue;
    refinements.insert(refinements.end(), split.begin(), split.end());
    if (refinements.size() > maxRefinements)
      return false;
  }
  return false;
}

// One-sided derivatives of exact sections distinguish inherited jet incompatibility
// from derivative error in the fitted sheets. Physical arclength steps and their
// halving keep the trigger independent of the rail's parameterization.
static bool sectionJetBreak(const Rail& a,
                            const Rail& b,
                            Stations&   atA,
                            Stations&   atB,
                            double      radius)
{
  const double length = std::min({railLength(a, a.first, a.last),
                                  railLength(b, b.first, b.last),
                                  radius});
  if (!(length > limit::degenerate * std::max(a.scale, b.scale)))
    return false;
  const auto curveAt = [&](const Rail& rail,
                           Stations&   stations,
                           double      endpoint,
                           double      distance) -> Handle(Geom_BSplineCurve) {
    double t = endpoint;
    if (distance != 0)
    {
      const auto along = railParameterAt(rail, endpoint, distance);
      if (!along)
        return {};
      t = *along;
    }
    Station station;
    Section section;
    if (!stations.at(t, station) || !sectionAt(rail, station, section))
      return {};
    return sectionCurve(section);
  };
  const auto endA = curveAt(a, atA, a.last, 0), endB = curveAt(b, atB, b.first, 0);
  if (endA.IsNull() || endB.IsNull())
    return false;
  double previous = -1;
  int    stable   = 0;
  for (int refinement = 0; refinement < 6; ++refinement)
  {
    const double h  = length * .005 / std::pow(2.0, refinement);
    const auto   a1 = curveAt(a, atA, a.last, -h), a2 = curveAt(a, atA, a.last, -2 * h),
               b1 = curveAt(b, atB, b.first, h), b2 = curveAt(b, atB, b.first, 2 * h);
    if (a1.IsNull() || a2.IsNull() || b1.IsNull() || b2.IsNull())
      return false;
    double worst = 0;
    for (int i = 1; i < 16; ++i)
    {
      const double v = double(i) / 16;
      gp_Pnt       pa, pb;
      gp_Vec       av, bv;
      endA->D1(v, pa, av);
      endB->D1(v, pb, bv);
      const gp_Vec au = (gp_Vec(a1->Value(v), pa) * 4 - gp_Vec(a2->Value(v), pa)) / (2 * h);
      const gp_Vec bu = (gp_Vec(pb, b1->Value(v)) * 4 - gp_Vec(pb, b2->Value(v))) / (2 * h);
      const gp_Vec na = au.Crossed(av), nb = bu.Crossed(bv);
      if (!std::isfinite(na.SquareMagnitude()) || !std::isfinite(nb.SquareMagnitude()))
        fail(a.op, no_solution, "an inherited section has a non-finite jet");
      if (na.Magnitude() <= limit::degenerate * a.scale
          || nb.Magnitude() <= limit::degenerate * b.scale)
        return false;
      const double angle = na.Angle(nb);
      if (!std::isfinite(angle))
        fail(a.op, no_solution, "an inherited section has a non-finite normal");
      worst = std::max(worst, angle);
    }
    stable =
      previous >= 0 && std::abs(worst - previous) <= std::max(tol::angularResolution, .05 * worst)
        ? stable + 1
        : 0;
    if (stable >= 2)
      return std::min(worst, previous) > 2 * tol::angularResolution;
    previous = worst;
  }
  return false;
}

struct SplitRequirements
{
  std::map<int, std::vector<SectionFrame>> frames;
  std::set<size_t>                         framedJoins;
  std::set<size_t>                         localJoins;
};

static SplitRequirements splitFrames(const char*                                op,
                                     const std::vector<Piece>&                  pieces,
                                     const std::vector<std::optional<Station>>& seeds,
                                     const std::vector<int>&                    group)
{
  SplitRequirements requirements;
  auto&             frames = requirements.frames;
  for (size_t m = 1; m < group.size(); ++m)
  {
    const auto& a = pieces[group[m - 1]];
    const auto& b = pieces[group[m]];
    if (a.rail.law != Law::G2 || b.rail.law != Law::G2)
      continue;
    const auto& sa = a.fitted.blend;
    const auto& sb = b.fitted.blend;
    if (sa->NbVPoles() != sb->NbVPoles())
      continue;
    const double tol = std::min(fitTolerance(a.rail), fitTolerance(b.rail));
    if (sa->Pole(sa->NbUPoles(), 1).Distance(sb->Pole(1, 1)) > tol
        || sa->Pole(sa->NbUPoles(), sa->NbVPoles()).Distance(sb->Pole(1, sb->NbVPoles())) > tol)
      continue;
    const auto& joint = a.joint[1];
    const int   side  = joint.side;
    double      f, l;
    const auto  curve = BRep_Tool::Curve(joint.crossed, f, l);
    if (curve.IsNull())
      continue;
    const gp_Pnt                foot = sa->Pole(sa->NbUPoles(), side ? sa->NbVPoles() : 1);
    GeomAPI_ProjectPointOnCurve on(foot, curve, f, l);
    if (on.NbPoints() == 0)
      continue;
    gp_Pnt point;
    gp_Vec tangent;
    curve->D1(on.LowerDistanceParameter(), point, tangent);
    Stations stationsA(a.rail), stationsB(b.rail);
    stationsA.march(seeds[group[m - 1]]);
    stationsB.march(seeds[group[m]]);
    Station stationA, stationB;
    if (!stationsA.at(a.rail.last, stationA) || !stationsB.at(b.rail.first, stationB))
      continue;
    Jet jetA, jetB;
    if (!jetAt(a.rail.s[side], stationA.uv[side][0], stationA.uv[side][1], jetA)
        || !jetAt(b.rail.s[side], stationB.uv[side][0], stationB.uv[side][1], jetB))
      continue;
    const gp_Vec normalA = stationA.normal[side], normalB = stationB.normal[side];
    if (normalA.Angle(normalB) > tol::angularResolution)
      continue;
    gp_Vec transverse = normalA.Crossed(tangent);
    if (transverse.Magnitude() <= tol::angle)
      continue;
    transverse.Normalize();
    const SurfaceJet
      supportA{a.rail.op, jetA.p, jetA.du, jetA.dv, jetA.duu, jetA.dvv, jetA.duv, normalA};
    const SurfaceJet
      supportB{b.rail.op, jetB.p, jetB.du, jetB.dv, jetB.duu, jetB.dvv, jetB.duv, normalB};
    const gp_Vec along = transverse.Crossed(normalA).Normalized();
    const double aa    = supportA.curvature(along, along) - supportB.curvature(along, along),
                 bb    = supportA.curvature(transverse, transverse)
                      - supportB.curvature(transverse, transverse),
                 ab = supportA.curvature(along, transverse) - supportB.curvature(along, transverse);
    const double jump   = std::abs((aa + bb) / 2) + std::hypot((aa - bb) / 2, ab);
    const double radius = std::max(stationA.centre.Distance(stationA.foot[side]),
                                   stationB.centre.Distance(stationB.foot[side]));
    if (!std::isfinite(jump) || !std::isfinite(radius))
      fail(op, no_solution, "a support crossing has a non-finite curvature jet");
    const bool continuousCurvature = jump * radius * radius <= tol;
    if (continuousCurvature && !sectionJetBreak(a.rail, b.rail, stationsA, stationsB, radius))
      continue;
    requirements.framedJoins.insert(m);
    // A third-order support break needs longitudinal normal coupling, but
    // does not require replacing the calibrated section family by a frame.
    if (continuousCurvature)
    {
      requirements.localJoins.insert(m);
      continue;
    }
    // Compact geometric windows cannot overlap with a frame at the other
    // end. Ball radius and inherited support-curvature radius set the geometric reach.
    const double lengthA = railLength(a.rail, a.rail.first, a.rail.last);
    const double lengthB = railLength(b.rail, b.rail.first, b.rail.last);
    gp_Pnt       contact;
    gp_Vec       alongContact, acrossContact;
    sa->D1(1, side ? 1 : 0, contact, alongContact, acrossContact);
    if (alongContact.Magnitude() <= limit::degenerate * a.rail.scale
        || alongContact.Normalized().Crossed(tangent.Normalized()).Magnitude() <= tol::coarseAngle)
      fail(op, no_solution, "a support crossing is tangent to the blend contact");
    // Lift the seam line to one oriented tangent at the joint, then transport
    // that orientation continuously through both windows.
    if (tangent.Dot(side ? -acrossContact : acrossContact) < 0)
      tangent.Reverse();
    // A continuous-curvature crossing has no finite II-jump radius. Limit
    // both windows by the shorter adjacent support span, so a remote flat
    // extension cannot change this local construction.
    const double reach =
      continuousCurvature ? std::min(lengthA, lengthB) : std::max(radius, 1 / jump);
    frames[group[m - 1]].push_back({a.rail.last, std::min(reach, lengthA), side, tangent});
    frames[group[m]].push_back({b.rail.first, std::min(reach, lengthB), side, tangent});
  }
  for (auto& [index, at] : frames)
    for (auto& frame : at)
    {
      if (std::count_if(at.begin(),
                        at.end(),
                        [&](const SectionFrame& other) { return other.side == frame.side; })
          > 1)
      {
        const auto& rail = pieces[index].rail;
        frame.reach      = std::min(frame.reach,
                               railLength(rail, rail.first, rail.last) / 2);
      }
    }
  return requirements;
}

//! A G2 sheet's knot where the law itself kinks (a support's C2 knot, whose curvature's slope
//! jumps; the section law changing form) is coupled as a split's join is: the sheet is cut there,
//! its two parts' derivatives across brought together in a local window, and it is joined again. A
//! kink within the G1 limit is left as it is, and so is one coupling cannot take: sheetJoins
//! refuses it.
static bool coupleInnerKinks(const char* op, Piece& piece)
{
  bool coupled = false;
  if (piece.rail.law != Law::G2 || piece.fitted.blend.IsNull() || piece.exact())
    return coupled;
  for (int i = 2; i < piece.fitted.blend->NbUKnots(); ++i)
  {
    const auto& surface = piece.fitted.blend;
    if (surface->IsURational() || surface->IsVRational()
        || surface->UMultiplicity(i) < surface->UDegree() - 1)
      continue;
    const double w = surface->UKnot(i);
    if (w <= piece.uAt(0) || w >= piece.uAt(1))
      continue;
    double tilt = 0;
    for (int j = 2; j <= surface->NbVKnots(); ++j)
      for (const double q : {0.25, 0.5, 0.75})
      {
        const double v = surface->VKnot(j - 1) + (surface->VKnot(j) - surface->VKnot(j - 1)) * q;
        SurfaceJet   a{op}, b{op};
        surfaceD2OnUSpan(surface, w, v, i - 1, i, a);
        surfaceD2OnUSpan(surface, w, v, i, i + 1, b);
        a.orient(false);
        b.orient(false);
        tilt = std::max(tilt, a.n.Angle(b.n));
      }
    if (!(tilt * 180 / M_PI > limit::seamDegrees))
      continue;
    // Its two parts, each over 0 to 1, on their own parts of the rail.
    const auto part = [&](double lo, double hi) {
      const auto cut = Handle(Geom_BSplineSurface)::DownCast(surface->Copy());
      cut->Segment(lo, hi, surface->VKnot(1), surface->VKnot(surface->NbVKnots()));
      NCollection_Array1<double> knots(1, cut->NbUKnots());
      for (int k = 1; k <= cut->NbUKnots(); ++k)
        knots(k) = (cut->UKnot(k) - lo) / (hi - lo);
      cut->SetUKnots(knots);
      return cut;
    };
    auto         a = part(0, w), b = part(w, 1);
    Rail         railA = piece.rail, railB = piece.rail;
    const double t = piece.rail.t(w);
    railA.last = railA.to = t;
    railB.first = railB.from = t;
    if (!coupleSplitNormals(op,
                            a,
                            b,
                            fitTolerance(piece.rail),
                            railA,
                            railB,
                            true,
                            piece.fitted.contact[0]->Value(w),
                            0,
                            limit::seamDegrees * M_PI / 180 / 2)
        || a->NbVPoles() != b->NbVPoles())
      continue;
    // Joined again, one row at the join.
    const int                  nu = a->NbUPoles() + b->NbUPoles() - 1, nv = a->NbVPoles();
    NCollection_Array2<gp_Pnt> poles(1, nu, 1, nv);
    for (int v = 1; v <= nv; ++v)
    {
      for (int u = 1; u <= a->NbUPoles(); ++u)
        poles(u, v) = a->Pole(u, v);
      for (int u = 2; u <= b->NbUPoles(); ++u)
        poles(a->NbUPoles() + u - 1, v) = b->Pole(u, v);
    }
    const int                  ka = a->NbUKnots(), kb = b->NbUKnots(), degree = a->UDegree();
    NCollection_Array1<double> uKnots(1, ka + kb - 1), vKnots(1, a->NbVKnots());
    NCollection_Array1<int>    uMults(1, ka + kb - 1), vMults(1, a->NbVKnots());
    for (int k = 1; k <= ka; ++k)
    {
      uKnots(k) = a->UKnot(k) * w;
      uMults(k) = k == ka ? degree : a->UMultiplicity(k);
    }
    for (int k = 2; k <= kb; ++k)
    {
      uKnots(ka + k - 1) = w + b->UKnot(k) * (1 - w);
      uMults(ka + k - 1) = b->UMultiplicity(k);
    }
    for (int k = 1; k <= a->NbVKnots(); ++k)
    {
      vKnots(k) = a->VKnot(k);
      vMults(k) = a->VMultiplicity(k);
    }
    piece.fitted.blend =
      new Geom_BSplineSurface(poles, uKnots, vKnots, uMults, vMults, degree, a->VDegree());
    piece.sheet = Sheet{piece.fitted.blend};
    coupled     = true;
  }
  return coupled;
}

//! How much the G2 law's curvature along the rail jumps at a place of a sheet: the exact sections'
//! one-sided second differences there. The law reads a support's curvature, so a support that is
//! C^k gives a sheet of C^(k-2) at best (a cubic's knot, C2, gives a jump here), which no fit can
//! take away. Physical steps, as sectionJetBreak's.
static double lawBreak(const Rail& a,
                       Stations&   atA,
                       double      ta,
                       const Rail& b,
                       Stations&   atB,
                       double      tb,
                       double      radius)
{
  const auto curveAt = [&](const Rail& rail,
                           Stations&   stations,
                           double      from,
                           double      distance) -> Handle(Geom_BSplineCurve) {
    double t = from;
    if (distance != 0)
    {
      const auto along = railParameterAt(rail, from, distance);
      if (!along)
        return {};
      t = *along;
    }
    Station station;
    Section section;
    if (!stations.at(t, station) || !sectionAt(rail, station, section))
      return {};
    return sectionCurve(section);
  };
  const double h  = radius * .005;
  const auto   a0 = curveAt(a, atA, ta, 0), a1 = curveAt(a, atA, ta, -h),
             a2 = curveAt(a, atA, ta, -2 * h);
  const auto b0 = curveAt(b, atB, tb, 0), b1 = curveAt(b, atB, tb, h),
             b2 = curveAt(b, atB, tb, 2 * h);
  if (a0.IsNull() || a1.IsNull() || a2.IsNull() || b0.IsNull() || b1.IsNull() || b2.IsNull())
    return 0;
  double bent = 0;
  for (int i = 1; i < 16; ++i)
  {
    const double v = double(i) / 16;
    gp_Pnt       pa, pb;
    gp_Vec       av, bv;
    a0->D1(v, pa, av);
    b0->D1(v, pb, bv);
    const gp_Pnt qa = a1->Value(v), qb = b1->Value(v), ra = a2->Value(v), rb = b2->Value(v);
    const gp_Vec au(qa, pa), bu(pb, qb);
    const gp_Vec auu = gp_Vec(qa, pa) - gp_Vec(ra, qa), buu = gp_Vec(qb, rb) - gp_Vec(pb, qb);
    const gp_Vec na = au.Crossed(av), nb = bu.Crossed(bv);
    if (!(na.Magnitude() > limit::degenerate) || !(nb.Magnitude() > limit::degenerate))
      continue;
    const gp_Vec n = na.Normalized();
    bent =
      std::max(bent,
               std::abs(n.Dot(auu) / au.SquareMagnitude() - n.Dot(buu) / bu.SquareMagnitude()));
  }
  return bent;
}

//! Inside a G2 blend's face, each join of its sheet where its U derivatives may break (a knot below
//! C2: where the parts of a split were joined, or a support's knot): held to what its edges are
//! held to, past what the law itself breaks there (lawBreak). Over the G2 matrix a join's break is
//! 0.1 to 6 times the law's, median 1, and the next is over 400 times: limit::lawShare, measured,
//! not derived.
static void sheetJoins(const char*                                op,
                       const std::vector<Piece>&                  pieces,
                       const std::vector<std::optional<Station>>& seeds)
{
  NCollection_Map<Standard_Address> seen;
  for (int k = 0; k < int(pieces.size()); ++k)
  {
    if (pieces[k].rail.law != Law::G2)
      continue;
    const auto surface = Handle(Geom_BSplineSurface)::DownCast(pieces[k].sheet.surface);
    if (surface.IsNull() || pieces[k].sheet.swapped || !seen.Add(surface.get()))
      continue;
    std::vector<int> group;
    for (int m = 0; m < int(pieces.size()); ++m)
      if (pieces[m].sheet.surface.get() == surface.get())
        group.push_back(m);
    const auto holding = [&](double w, bool ahead) -> std::optional<std::pair<int, double>> {
      for (int m : group)
      {
        const double at = (w - pieces[m].sheet.alpha) / pieces[m].sheet.beta, u0 = pieces[m].uAt(0),
                     u1 = pieces[m].uAt(1);
        if (ahead ? at >= u0 && at < u1 : at > u0 && at <= u1)
          return std::pair{m, at};
      }
      return {};
    };
    for (int i = 2; i < surface->NbUKnots(); ++i)
    {
      if (surface->UMultiplicity(i) < surface->UDegree() - 1)
        continue;
      const double w      = surface->UKnot(i);
      const auto   before = holding(w, false), after = holding(w, true);
      if (!before || !after)
        continue;
      double tilt = 0, bent = 0;
      for (int j = 2; j <= surface->NbVKnots(); ++j)
        for (const double q : {0.25, 0.5, 0.75})
        {
          const double v = surface->VKnot(j - 1) + (surface->VKnot(j) - surface->VKnot(j - 1)) * q;
          SurfaceJet   a{op}, b{op};
          surfaceD2OnUSpan(surface, w, v, i - 1, i, a);
          surfaceD2OnUSpan(surface, w, v, i, i + 1, b);
          a.orient(false);
          b.orient(false);
          tilt = std::max(tilt, a.n.Angle(b.n));
          bent = std::max(bent, bendDifference(a, b, a.v.Normalized()).first);
        }
      const Piece& p     = pieces[before->first];
      const Piece& q     = pieces[after->first];
      const double scale = p.rail.scale, degrees = 180 / M_PI;
      if (tilt * degrees <= limit::seamDegrees && bent * scale <= limit::seamBend)
        continue;
      // The crease a support itself has there, read on either side of the join, is the only
      // excuse for one here.
      double       crease = 0;
      const double aside =
        limit::knotMerge * std::min(w - surface->UKnot(i - 1), surface->UKnot(i + 1) - w);
      for (int side = 0; side < 2; ++side)
      {
        if (p.rail.s[side].cliffCurve || q.rail.s[side].cliffCurve)
          continue;
        const gp_Pnt2d uvp = p.fitted.contact[side]->Value(before->second - aside / p.sheet.beta),
                       uvq = q.fitted.contact[side]->Value(after->second + aside / q.sheet.beta);
        Jet jp, jq;
        if (jetAt(p.rail.s[side], uvp.X(), uvp.Y(), jp)
            && jetAt(q.rail.s[side], uvq.X(), uvq.Y(), jq))
          crease = std::max(crease, jp.n.Angle(jq.n));
      }
      if (tilt * degrees > limit::seamDegrees && tilt > crease)
        fail(op,
             no_solution,
             "a G2 blend is not tangent continuous inside its face",
             {p.rail.edge},
             FailureKind::Curvature);
      if (bent * scale <= limit::seamBend || !seeds[before->first] || !seeds[after->first])
        continue;
      Stations atP(p.rail), atQ(q.rail);
      atP.march(seeds[before->first]);
      atQ.march(seeds[after->first]);
      const double lawBent = lawBreak(p.rail,
                                      atP,
                                      p.rail.t(before->second),
                                      q.rail,
                                      atQ,
                                      q.rail.t(after->second),
                                      scale);
      gp_Pnt       place;
      gp_Vec       alongU, acrossU;
      surface->D1(w,
                  (surface->VKnot(1) + surface->VKnot(surface->NbVKnots())) / 2,
                  place,
                  alongU,
                  acrossU);
      const double span =
        alongU.Magnitude() * std::min(w - surface->UKnot(i - 1), surface->UKnot(i + 1) - w);
      if (bent * scale > limit::seamBend + limit::lawShare * lawBent * scale
          && !unresolved(bent, span, scale))
        fail(op,
             no_solution,
             "a G2 blend is not curvature continuous inside its face",
             {p.rail.edge},
             FailureKind::Curvature);
    }
  }
}

void joinSplitPieces(const char*                                op,
                     std::vector<Piece>&                        pieces,
                     const std::vector<std::optional<Station>>& seeds)
{
  // Keep the single-face representation of one edge split by support crossings.
  // Framed inherited curvature breaks retain G2 contacts and a common join normal.
  // Other groups retain their established section law and join construction.

  for (int k = 0; k < int(pieces.size()); ++k)
  {
    checkCancelled(op);
    if (pieces[k].joint[0].meets == Meets::Split || pieces[k].joint[1].meets != Meets::Split)
      continue;
    std::vector<int> group = {k};
    while (pieces[group.back()].joint[1].meets == Meets::Split)
      group.push_back(pieces[group.back()].joint[1].piece);
    const auto  requirements = splitFrames(op, pieces, seeds, group);
    const auto& frames       = requirements.frames;
    for (const auto& [index, at] : frames)
    {
      Stations stations(pieces[index].rail);
      stations.march(seeds[index]);
      const Fitted terminalReference = pieces[index].fitted;
      pieces[index].fitted = fit(pieces[index].rail, stations, false, at, &terminalReference);
      pieces[index].sheet  = Sheet{pieces[index].fitted.blend};
      Dips dips{pieces[index].rail, pieces[index].fitted};
      checkFold(dips);
    }
    const bool framedGroup = !requirements.framedJoins.empty();
    if (framedGroup)
      for (size_t m = 1; m < group.size(); ++m)
      {
        if (requirements.framedJoins.find(m) == requirements.framedJoins.end())
          continue;
        const bool localLaw = requirements.localJoins.count(m) > 0;
        gp_Pnt2d   supportUV;
        if (localLaw)
        {
          const auto& contact = pieces[group[m - 1]].fitted.contact[0];
          if (contact.IsNull())
            fail(op, no_solution, "a split contact has no support parameter curve");
          supportUV = contact->Value(1);
        }
        // How far apart the two sides' own laws put the section at the join.
        double lawGap = 0;
        if (!localLaw)
        {
          const Rail& railA = pieces[group[m - 1]].rail;
          const Rail& railB = pieces[group[m]].rail;
          Stations    atA(railA), atB(railB);
          atA.march(seeds[group[m - 1]]);
          atB.march(seeds[group[m]]);
          Station x, y;
          Section xs, ys;
          if (atA.at(railA.last, x) && atB.at(railB.first, y) && sectionAt(railA, x, xs)
              && sectionAt(railB, y, ys) && xs.poles.size() == ys.poles.size())
            for (size_t i = 0; i < xs.poles.size(); ++i)
              lawGap = std::max(lawGap, xs.poles[i].Distance(ys.poles[i]));
        }
        const bool glued = coupleSplitNormals(
          op,
          pieces[group[m - 1]].fitted.blend,
          pieces[group[m]].fitted.blend,
          std::min(fitTolerance(pieces[group[m - 1]].rail), fitTolerance(pieces[group[m]].rail)),
          pieces[group[m - 1]].rail,
          pieces[group[m]].rail,
          localLaw,
          supportUV,
          lawGap);
        if (!glued)
          fail(op, no_solution, "the coupled jets at a split could not be fitted");
      }
    bool coupled = false;
    for (int index : group)
      coupled = coupleInnerKinks(op, pieces[index]) || coupled;
    if (framedGroup || coupled)
    {
      std::vector<Handle(Geom_BSplineSurface)> sections;
      for (int index : group)
        sections.push_back(pieces[index].fitted.blend);
      if (!commonSections(sections))
        fail(op, no_solution, "the section bases of split blends do not agree");
    }
    const auto&                                         head   = pieces[k].fitted.blend;
    const int                                           degree = head->UDegree();
    std::vector<double>                                 knots;
    std::vector<int>                                    mults;
    std::vector<std::vector<std::pair<gp_Pnt, double>>> rows;
    for (size_t m = 0; m < group.size(); ++m)
    {
      const auto& part         = pieces[group[m]].fitted.blend;
      const auto  incompatible = [&] {
        fail(op,
             no_solution,
             "the parts of a blend split onto the next face do not join",
              {pieces[k].rail.edge});
      };
      // Sharing one row is valid only for clamped, equally parameterised sections.
      // Never replace a section by another merely because its pole count agrees.
      if (part->IsUPeriodic() || part->IsVPeriodic() || part->UDegree() != degree
          || part->VDegree() != head->VDegree() || part->NbVPoles() != head->NbVPoles()
          || part->NbVKnots() != head->NbVKnots() || std::abs(part->UKnot(1)) > tol::parameter
          || std::abs(part->UKnot(part->NbUKnots()) - 1) > tol::parameter
          || part->UMultiplicity(1) != degree + 1
          || part->UMultiplicity(part->NbUKnots()) != degree + 1)
        incompatible();
      for (int j = 1; j <= head->NbVKnots(); ++j)
        if (part->VMultiplicity(j) != head->VMultiplicity(j)
            || std::abs(part->VKnot(j) - head->VKnot(j)) > tol::parameter)
          incompatible();
      // Rational sections may differ by one common homogeneous scale. Apply it
      // to the whole part, preserving its geometry and the shared row's weights.
      const double weightScale = m ? rows.back()[0].second / part->Weight(1, 1) : 1;
      if (m)
      {
        const double tolerance =
          std::min(fitTolerance(pieces[group[m - 1]].rail), fitTolerance(pieces[group[m]].rail));
        double poleGap = 0, weightChange = 0;
        for (int j = 1; j <= part->NbVPoles(); ++j)
        {
          const auto&  previous = rows.back()[j - 1];
          const double weight   = part->Weight(1, j) * weightScale;
          if (!(previous.second > 0 && weight > 0))
            incompatible();
          poleGap = std::max(poleGap, previous.first.Distance(part->Pole(1, j)));
          weightChange =
            std::max(weightChange,
                     std::abs(previous.second - weight) / std::min(previous.second, weight));
        }
        // Replacing the first row affects only the first U span. The
        // original surface lies in its active pole hull. Pole gap
        // plus relative weight change * hull diameter bounds
        // the geometric change. A dimensionless weight threshold alone
        // rejects harmless rounding in a small section.
        Bnd_Box hull;
        for (int i = 1; i <= degree + 1; ++i)
          for (int j = 1; j <= part->NbVPoles(); ++j)
            hull.Add(part->Pole(i, j));
        hull.SetGap(0);
        if (poleGap + weightChange * std::sqrt(hull.SquareExtent()) > tolerance)
          incompatible();
      }
      for (int i = 1; i <= part->NbUKnots(); ++i)
      {
        if (m > 0 && i == 1)
        {
          mults.back() = degree;
          continue;
        }
        knots.push_back(part->UKnot(i) + double(m));
        mults.push_back(part->UMultiplicity(i));
      }
      for (int i = 1; i <= part->NbUPoles(); ++i)
      {
        std::vector<std::pair<gp_Pnt, double>> row;
        for (int j = 1; j <= part->NbVPoles(); ++j)
          row.push_back({part->Pole(i, j), part->Weight(i, j) * weightScale});
        if (m > 0 && i == 1)
        {
          continue;
        }
        rows.push_back(row);
      }
    }
    NCollection_Array2<gp_Pnt> poles(1, int(rows.size()), 1, head->NbVPoles());
    NCollection_Array2<double> weights(1, int(rows.size()), 1, head->NbVPoles());
    for (int i = 0; i < int(rows.size()); ++i)
      for (int j = 0; j < head->NbVPoles(); ++j)
      {
        poles(i + 1, j + 1)   = rows[i][j].first;
        weights(i + 1, j + 1) = rows[i][j].second;
      }
    NCollection_Array1<double> uk(1, int(knots.size())), vk(1, head->NbVKnots());
    NCollection_Array1<int>    um(1, int(knots.size())), vm(1, head->NbVKnots());
    for (int i = 0; i < int(knots.size()); ++i)
    {
      uk(i + 1) = knots[i];
      um(i + 1) = mults[i];
    }
    for (int i = 1; i <= head->NbVKnots(); ++i)
    {
      vk(i) = head->VKnot(i);
      vm(i) = head->VMultiplicity(i);
    }
    Handle(Geom_BSplineSurface) joined =
      new Geom_BSplineSurface(poles, weights, uk, vk, um, vm, degree, head->VDegree());
    std::vector<double> starts, widths;
    for (size_t m = 0; m < group.size(); ++m)
    {
      starts.push_back(double(m));
      widths.push_back(1);
    }
    // A join is stored C0 in the parameter however smooth it is, and OpenCascade's offset refuses a
    // surface so stored. With spans in proportion to each part's share of the rail's parameter, on
    // which the sections go on from each other, its derivatives across are the same where the
    // law's are: stored C1 where the poles then carry on to within a hundredth of the fit's
    // tolerance and it is tangent to begin with. One unit a part on the average, as they were.
    double along = 0;
    for (int index : group)
      along += pieces[index].rail.to - pieces[index].rail.from;
    std::vector<double> scaledStarts, scaledWidths;
    for (int index : group)
    {
      scaledStarts.push_back(scaledStarts.empty() ? 0 : scaledStarts.back() + scaledWidths.back());
      scaledWidths.push_back((pieces[index].rail.to - pieces[index].rail.from)
                             * double(group.size()) / along);
    }
    NCollection_Array1<double> scaledKnots(uk);
    for (int i = 1; i <= uk.Upper(); ++i)
    {
      const int m    = std::min(int(uk(i)), int(group.size()) - 1);
      scaledKnots(i) = scaledStarts[m] + (uk(i) - m) * scaledWidths[m];
    }
    const Handle(Geom_BSplineSurface) smooth =
      new Geom_BSplineSurface(poles, weights, scaledKnots, vk, um, vm, degree, head->VDegree());
    double tolerance = HUGE_VAL;
    for (int index : group)
      tolerance = std::min(tolerance, fitTolerance(pieces[index].rail));
    bool stored = false;
    for (int i = 2; i < smooth->NbUKnots(); ++i)
    {
      if (smooth->UMultiplicity(i) < degree)
        continue;
      double tilt = 0;
      for (int j = 2; j <= smooth->NbVKnots(); ++j)
        for (const double q : {0.25, 0.5, 0.75})
        {
          const double v = smooth->VKnot(j - 1) + (smooth->VKnot(j) - smooth->VKnot(j - 1)) * q;
          SurfaceJet   a{op}, b{op};
          surfaceD2OnUSpan(smooth, smooth->UKnot(i), v, i - 1, i, a);
          surfaceD2OnUSpan(smooth, smooth->UKnot(i), v, i, i + 1, b);
          a.orient(false);
          b.orient(false);
          tilt = std::max(tilt, a.n.Angle(b.n));
        }
      if (tilt * 180 / M_PI <= limit::seamDegrees
          && smooth->RemoveUKnot(i, degree - 1, tolerance / 100))
        stored = true;
    }
    if (stored)
    {
      joined = smooth;
      starts = scaledStarts;
      widths = scaledWidths;
    }
    if (!requirements.framedJoins.empty())
      for (int ui = 2; ui <= joined->NbUKnots(); ++ui)
      {
        checkCancelled(op);
        for (int vi = 2; vi <= joined->NbVKnots(); ++vi)
          for (double fu : {0.211324865405187, 0.788675134594813})
            for (double fv : {0.211324865405187, 0.788675134594813})
            {
              const double u =
                joined->UKnot(ui - 1) + (joined->UKnot(ui) - joined->UKnot(ui - 1)) * fu;
              const double v =
                joined->VKnot(vi - 1) + (joined->VKnot(vi) - joined->VKnot(vi - 1)) * fv;
              gp_Pnt point;
              gp_Vec du, dv;
              joined->D1(u, v, point, du, dv);
              const double metric = du.Magnitude() * dv.Magnitude();
              const double normal = du.Crossed(dv).Magnitude();
              if (!std::isfinite(metric) || !std::isfinite(normal) || !(metric > 0)
                  || !(normal > tol::coarseAngle * metric))
                fail(op, no_solution, "a split blend surface is singular");
            }
      }
    // Every framed join must retain the constructed common tangent plane.
    // Unframed groups keep the established section law and its proof contract.
    if (framedGroup)
      for (size_t m = 1; m < group.size(); ++m)
        if (requirements.framedJoins.find(m) != requirements.framedJoins.end())
          for (int span = 2; span <= joined->NbVKnots(); ++span)
            for (int sample = 1; sample <= 16; ++sample)
            {
              const double v = joined->VKnot(span - 1)
                               + (joined->VKnot(span) - joined->VKnot(span - 1)) * sample / 17;
              const double u = starts[m];
              SurfaceJet   left{op}, right{op};
              int          split = 2;
              while (split < joined->NbUKnots() && joined->UKnot(split) < u)
                ++split;
              surfaceD2OnUSpan(joined, u, v, split - 1, split, left);
              surfaceD2OnUSpan(joined, u, v, split, split + 1, right);
              left.orient(false);
              right.orient(false);
              if (!(left.n.Angle(right.n) <= tol::angularResolution))
                fail(op, no_solution, "the normal of split blends does not agree");
            }
    for (size_t m = 0; m < group.size(); ++m)
      pieces[group[m]].sheet = Sheet{joined, false, starts[m], widths[m], {0, 1}};
  }
  for (auto& piece : pieces)
    if (piece.joint[0].meets != Meets::Split && piece.joint[1].meets != Meets::Split)
      coupleInnerKinks(op, piece);
  sheetJoins(op, pieces, seeds);
}

} // namespace BRepExt::blending::detail
