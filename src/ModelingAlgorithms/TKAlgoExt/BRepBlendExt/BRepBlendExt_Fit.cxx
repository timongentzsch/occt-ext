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

// Blend: The blend approximated from its exact sections: knots at the
// data's breaks, Greville interpolation refined to the tolerance, and the fold
// check.

#include <BRepBlendExt_Internal.hxx>
#include <GeomAlgoExt_Hermite.hxx>

#include <Geom2dAdaptor_Curve.hxx>
#include <GeomLib_IsPlanarSurface.hxx>

namespace BRepExt::blend2::detail
{

std::vector<double> flatten(const std::map<double, int>& knots)
{
  std::vector<double> flat;
  for (const auto& [u, m] : knots)
    flat.insert(flat.end(), m, u);
  return flat;
}

static void addBreak(std::map<double, int>& breaks, double u, int m)
{
  if (!(u > tol::parameterMerge && u < 1 - tol::parameterMerge) || m <= 1)
    return;
  // A rail knot and a support crossing can describe the same event with
  // roundoff-different parameters. Keep its strongest continuity break.
  const double tolerance = tol::parameter;
  const auto   near      = breaks.lower_bound(u - tolerance);
  auto&        slot =
    near != breaks.end() && std::abs(near->first - u) <= tolerance ? near->second : breaks[u];
  slot = std::max(slot, std::min(m, kDegree));
}

std::map<double, int> breaksOf(const Rail& rail, Stations& stations, bool sought)
{
  const int  order  = rail.law == Law::G2 ? 2 : 1;
  const auto needed = [&](int degree, int multiplicity) {
    return kDegree - std::max(0, degree - multiplicity - order);
  };
  std::map<double, int> breaks;
  // The edge's knots: each station is the ball whose section's plane passes through the edge's
  // point there (ballAt), so the stations are as smooth as the edge's points, not its tangent:
  // a C1 edge (OpenCascade's approximated meetings, degree 8 with knots of multiplicity 7) makes
  // the blend C1 there, not C0.
  if (rail.curve.GetType() == GeomAbs_BSplineCurve)
  {
    const auto curve = rail.curve.BSpline();
    for (int k = 2; k < curve->NbKnots(); ++k)
      addBreak(breaks,
               rail.u(curve->Knot(k)),
               kDegree - std::max(0, curve->Degree() - curve->Multiplicity(k)));
  }

  // The supports' knot lines, each as a parameter (0 for u, 1 for v), a value and what it needs.
  struct Line
  {
    int    parameter;
    double value;
    int    multiplicity;
  };

  for (int side = 0; side < 2; ++side)
  {
    std::vector<Line> lines;
    if (rail.s[side].cliffCurve)
      continue; // its feet run along an edge, not across a face's knots
    const Handle(Geom_Surface) surface    = basisOf(rail.s[side].surface);
    const auto                 curveLines = [&](const Handle(Geom_Curve)& curve, int parameter) {
      if (const auto spline = Handle(Geom_BSplineCurve)::DownCast(basisOf(curve)))
        for (int k = 2; k < spline->NbKnots(); ++k)
          lines.push_back(
            {parameter, spline->Knot(k), needed(spline->Degree(), spline->Multiplicity(k))});
    };
    if (const auto spline = Handle(Geom_BSplineSurface)::DownCast(surface))
    {
      for (int k = 2; k < spline->NbUKnots(); ++k)
        lines.push_back({0, spline->UKnot(k), needed(spline->UDegree(), spline->UMultiplicity(k))});
      for (int k = 2; k < spline->NbVKnots(); ++k)
        lines.push_back({1, spline->VKnot(k), needed(spline->VDegree(), spline->VMultiplicity(k))});
    }
    else if (const auto extrusion = Handle(Geom_SurfaceOfLinearExtrusion)::DownCast(surface))
      curveLines(extrusion->BasisCurve(), 0);
    else if (const auto revolution = Handle(Geom_SurfaceOfRevolution)::DownCast(surface))
      curveLines(revolution->BasisCurve(), 1);
    // Where the surface ends and runs on along its tangent (runsOn), it is only G1: a break
    // where a ball's foot passes it, as marched (a scan asks for stations of its own, so a fit
    // that needs no such break asks for none).
    const auto ends = runsOn(rail.s[side].surface);
    const auto past = [&](int parameter, double value, int beyond) {
      for (const auto& [t, station] : stations.known())
        if ((station.uv[side][parameter] - value) * beyond > 0)
          return true;
      return false;
    };
    for (int parameter = 0; parameter < 2; ++parameter)
      if (ends[parameter])
      {
        if (sought || past(parameter, ends[parameter]->first, -1))
          lines.push_back({parameter, ends[parameter]->first, needed(1, 0)});
        if (sought || past(parameter, ends[parameter]->second, 1))
          lines.push_back({parameter, ends[parameter]->second, needed(1, 0)});
      }
    lines.erase(
      std::remove_if(lines.begin(), lines.end(), [](const Line& l) { return l.multiplicity <= 1; }),
      lines.end());
    if (lines.empty())
      continue;
    // Scan the feet along the rail and bisect each crossing to the parameter's resolution.
    constexpr int scan   = 256;
    const auto    footAt = [&](double u, Station& s) {
      if (!stations.at(rail.t(u), s))
        fail(rail.op,
             no_solution,
             "the ball cannot roll the length of the edge at this size",
                {rail.edge},
             FailureKind::SupportContact);
      return s.uv[side];
    };
    Station previous, current;
    footAt(0, previous);
    for (int k = 1; k <= scan; ++k)
    {
      footAt(double(k) / scan, current);
      for (const auto& line : lines)
      {
        const bool before = previous.uv[side][line.parameter] > line.value;
        if (before == (current.uv[side][line.parameter] > line.value))
          continue;
        double low = double(k - 1) / scan, high = double(k) / scan;
        for (int i = 0; i < 60 && high - low > limit::degenerate; ++i)
        {
          const double middle = (low + high) / 2;
          Station      at;
          footAt(middle, at);
          (at.uv[side][line.parameter] > line.value) == before ? low = middle : high = middle;
        }
        addBreak(breaks, (low + high) / 2, line.multiplicity);
      }
      previous = current;
    }
  }
  return breaks;
}

double fitTolerance(const Rail& rail)
{
  // The length tolerance, and for a blend smaller than a millimetre that share of its size, so
  // that a small blend is as true to its law as a large one; never finer than the rounding of
  // its points where it sits (tolerance::resolved).
  return std::max(tol::length * std::min(1.0, rail.scale) * (rail.law == Law::G2 ? 0.1 : 1),
                  tol::rounding(rail.magnitude + rail.scale));
}

Fitted fit(const Rail& rail, Stations& stations, bool contactPilot)
{
  return fit(rail, stations, contactPilot, {});
}

Fitted fit(const Rail&                      rail,
           Stations&                        stations,
           bool                             contactPilot,
           const std::vector<SectionFrame>& frames,
           const Fitted*                    terminalReference)
{
  const double          tolerance = fitTolerance(rail);
  std::map<double, int> knots     = {{0, kDegree + 1}, {1, kDegree + 1}};
  for (int k = 1; k < 8; ++k)
    knots[k / 8.0] = 1;
  auto naturalBreaks = breaksOf(rail, stations);
  // A compact frame is C2 at its window boundary. Preserve that actual law
  // break instead of asking a globally C4 fit to approximate it by refinement.
  for (const auto& frame : frames)
  {
    const auto u =
      railUAtDistance(rail, frame.t, frame.t == rail.last ? -frame.reach : frame.reach);
    if (!u)
      fail(rail.op, no_solution, "a section frame has no finite window boundary", {rail.edge});
    addBreak(naturalBreaks, *u, kDegree - 2);
  }
  // Only optional uniform knots yield to a nearby natural break. Removing a
  // previously inserted natural break could erase its higher multiplicity.
  for (auto it = std::next(knots.begin()); it != knots.end() && it->first < 1;)
  {
    const auto near = naturalBreaks.lower_bound(it->first - limit::knotMerge);
    if (near != naturalBreaks.end() && std::abs(near->first - it->first) <= limit::knotMerge)
      it = knots.erase(it);
    else
      ++it;
  }
  for (const auto& [u, m] : naturalBreaks)
    knots[u] = std::max(knots[u], m);
  Rail sectionRail = rail;
  if (contactPilot)
  {
    sectionRail.law = Law::Conic;
    sectionRail.rho = 0.5;
  }
  const auto exact = [&](double u, Station& s, Section& section) {
    if (!stations.at(rail.t(u), s))
      fail(rail.op,
           no_solution,
           "the ball cannot rest on both faces beside the edge there",
           {rail.edge},
           FailureKind::SupportContact);
    if (!sectionAt(sectionRail, s, section, frames))
      fail(rail.op,
           no_solution,
           "the section collapses between the faces beside the edge",
           {rail.edge});
  };
  Section reference;
  {
    Station s;
    exact(0.5, s, reference);
  }
  const int S = int(reference.poles.size());
  // Per site: the section's homogeneous poles (x, y, z, w) and the feet's
  // parameters.
  const int  width = 4 * S + 4;
  const auto row   = [&](double u) {
    Station s;
    Section section;
    exact(u, s, section);
    std::vector<double> r(width);
    for (int j = 0; j < S; ++j)
    {
      const double w = section.weights[j];
      r[4 * j]       = section.poles[j].X() * w;
      r[4 * j + 1]   = section.poles[j].Y() * w;
      r[4 * j + 2]   = section.poles[j].Z() * w;
      r[4 * j + 3]   = w;
    }
    for (int i = 0; i < 2; ++i)
    {
      r[4 * S + 2 * i]     = s.uv[i][0];
      r[4 * S + 2 * i + 1] = s.uv[i][1];
    }
    return r;
  };
  // Compact frame windows leave their exterior and their boundary unchanged
  // through order two. Preserve that existing terminal fit contract when an
  // interior re-fit changes the global interpolation basis. Test the actual
  // carried fitting endpoints, which need not equal rail.first/last.
  const auto unchangedTerminal = [&](double t) {
    if (!terminalReference || frames.empty())
      return false;
    for (const auto& frame : frames)
    {
      const double distance =
        GCPnts_AbscissaPoint::Length(rail.curve, std::min(t, frame.t), std::max(t, frame.t));
      if (!std::isfinite(distance) || !std::isfinite(frame.reach)
          || !(frame.reach > limit::degenerate * rail.scale))
        fail(rail.op, no_solution, "a section frame has no finite terminal distance", {rail.edge});
      if (!(distance >= frame.reach))
        return false;
    }
    return true;
  };
  const bool keepFirst = unchangedTerminal(rail.from), keepLast = unchangedTerminal(rail.to);
  std::vector<double> referenceFlat;
  if (keepFirst || keepLast)
  {
    const auto& base = terminalReference->blend;
    if (base->NbVPoles() != S || base->VDegree() != reference.degree
        || base->NbVKnots() != int(reference.knots.size()))
      fail(rail.op, no_solution, "a terminal jet has a different section basis", {rail.edge});
    for (int i = 1; i <= base->NbVKnots(); ++i)
      if (base->VKnot(i) != reference.knots[i - 1]
          || base->VMultiplicity(i) != reference.multiplicities[i - 1])
        fail(rail.op, no_solution, "a terminal jet has a different section basis", {rail.edge});
    for (int i = 1; i <= base->NbUKnots(); ++i)
      referenceFlat.insert(referenceFlat.end(), base->UMultiplicity(i), base->UKnot(i));
  }
  const auto terminalRow = [&](double u, int order) {
    const auto&                base = terminalReference->blend;
    NCollection_Array1<double> flat(1, int(referenceFlat.size()));
    for (size_t i = 0; i < referenceFlat.size(); ++i)
      flat(int(i) + 1) = referenceFlat[i];
    math_Matrix values(1, order + 1, 1, base->UDegree() + 1);
    int         first = 0;
    BSplCLib::EvalBsplineBasis(order, base->UDegree() + 1, flat, u, first, values);
    std::vector<double> row(width, 0);
    for (int j = 1; j <= S; ++j)
      for (int k = 1; k <= base->UDegree() + 1; ++k)
      {
        const int i = first + k - 1;
        if (i < 1 || i > base->NbUPoles())
          continue;
        const double h    = values(order + 1, k) * base->Weight(i, j);
        const auto&  pole = base->Pole(i, j);
        row[4 * (j - 1)] += h * pole.X();
        row[4 * (j - 1) + 1] += h * pole.Y();
        row[4 * (j - 1) + 2] += h * pole.Z();
        row[4 * (j - 1) + 3] += h;
      }
    for (int side = 0; side < 2; ++side)
    {
      gp_Pnt2d p;
      gp_Vec2d d, dd;
      terminalReference->contact[side]->D2(u, p, d, dd);
      row[4 * S + 2 * side]     = order == 0 ? p.X() : order == 1 ? d.X() : dd.X();
      row[4 * S + 2 * side + 1] = order == 0 ? p.Y() : order == 1 ? d.Y() : dd.Y();
    }
    if (!std::all_of(row.begin(), row.end(), [](double x) { return std::isfinite(x); }))
      fail(rail.op, no_solution, "a terminal fit jet is non-finite", {rail.edge});
    return row;
  };
  // Breaks found while refining, on the knot they fall on where they fall on one.
  const auto raise = [&](const std::map<double, int>& breaks) {
    bool more = false;
    for (const auto& [u, m] : breaks)
    {
      const auto near = knots.lower_bound(u - limit::knotMerge);
      int&       slot =
        near != knots.end() && near->first <= u + limit::knotMerge ? near->second : knots[u];
      if (slot < m)
      {
        slot = m;
        more = true;
      }
    }
    return more;
  };
  std::vector<double> worst;
  int                 bendRounds = 0;
  for (int refinement = 0; refinement < 16; ++refinement)
  {
    checkCancelled(rail.op);
    const auto flat  = flatten(knots);
    const int  count = int(flat.size()) - kDegree - 1;
    if (count > 4000)
      break;
    NCollection_Array1<double> flatKnots(1, int(flat.size()));
    for (size_t k = 0; k < flat.size(); ++k)
      flatKnots(int(k) + 1) = flat[k];
    std::vector<double> sites(count);
    for (int i = 0; i < count; ++i)
    {
      double g = 0;
      for (int j = 1; j <= kDegree; ++j)
        g += flat[i + j];
      sites[i] = g / kDegree;
    }
    std::vector<GeomAlgoExt::HermiteSite> jetSites;
    if (keepFirst || keepLast)
    {
      jetSites = GeomAlgoExt::hermiteSites(kDegree, flat, keepFirst ? 2 : 0, keepLast ? 2 : 0);
      if (jetSites.empty())
        fail(rail.op, no_solution, "invalid terminal Hermite sites", {rail.edge});
    }
    // Interpolation at the sites, all rows at once: the collocation matrix is
    // banded and, at Greville value sites (or terminal Hermite jets),
    // OpenCascade's banded elimination without pivoting solves it in linear
    // time.
    const auto basisAt = [&](double u, int& firstIndex, math_Matrix& values) {
      return BSplCLib::EvalBsplineBasis(0, kDegree + 1, flatKnots, u, firstIndex, values);
    };
    NCollection_Array1<double> parameters(1, count);
    NCollection_Array1<int>    orders(1, count);
    std::vector<double>        poles(size_t(count) * width);
    for (int i = 0; i < count; ++i)
    {
      const GeomAlgoExt::HermiteSite site =
        jetSites.empty() ? GeomAlgoExt::HermiteSite{sites[i], 0} : jetSites[i];
      parameters(i + 1)   = site.parameter;
      orders(i + 1)       = site.order;
      const bool terminal = (keepFirst && site.parameter == 0) || (keepLast && site.parameter == 1);
      const auto r = terminal ? terminalRow(site.parameter, site.order) : row(site.parameter);
      std::copy(r.begin(), r.end(), poles.begin() + size_t(i) * width);
    }
    int singular = 0;
    BSplCLib::Interpolate(kDegree, flatKnots, parameters, orders, width, poles[0], singular);
    if (singular != 0)
      fail(rail.op, no_solution, "the blend's interpolation is singular", {rail.edge});
    std::vector<std::vector<double>> coefficients(width, std::vector<double>(count));
    for (int c = 0; c < width; ++c)
      for (int i = 0; i < count; ++i)
        coefficients[c][i] = poles[size_t(i) * width + c];
    const auto evaluate = [&](double u) {
      math_Matrix values(1, 1, 1, kDegree + 1);
      int         firstIndex = 0;
      basisAt(u, firstIndex, values);
      std::vector<double> r(width, 0.0);
      for (int c = 0; c < width; ++c)
        for (int k = 1; k <= kDegree + 1; ++k)
        {
          const int i = firstIndex + k - 2;
          if (i >= 0 && i < count)
            r[c] += values(1, k) * coefficients[c][i];
        }
      return r;
    };
    // A G2 blend's normal curvature along each contact, against its support's
    // along the feet: the fit is G2 to its supports only so far as its contact
    // rows are, which the positions alone do not bound (beside a C1 knot where
    // an edge's parameter kinks, they bend apart by 3e-4 / r within the length
    // tolerance).
    const auto alongBend = [&](double u, int side) {
      math_Matrix values(1, 3, 1, kDegree + 1);
      int         firstIndex = 0;
      BSplCLib::EvalBsplineBasis(2, kDegree + 1, flatKnots, u, firstIndex, values);
      const auto row = [&](int c, int d) {
        double sum = 0;
        for (int k = 1; k <= kDegree + 1; ++k)
        {
          const int i = firstIndex + k - 2;
          if (i >= 0 && i < count)
            sum += values(d + 1, k) * coefficients[c][i];
        }
        return sum;
      };
      const int j = side == 0 ? 0 : S - 1;
      gp_XYZ    X[3];
      double    w[3];
      for (int d = 0; d < 3; ++d)
      {
        X[d] = gp_XYZ(row(4 * j, d), row(4 * j + 1, d), row(4 * j + 2, d));
        w[d] = row(4 * j + 3, d);
      }
      const gp_XYZ P = X[0] / w[0], P1 = (X[1] - P * w[1]) / w[0],
                   P2 = (X[2] - P1 * (2 * w[1]) - P * w[2]) / w[0];
      const double U = row(4 * S + 2 * side, 0), V = row(4 * S + 2 * side + 1, 0);
      const double U1 = row(4 * S + 2 * side, 1), V1 = row(4 * S + 2 * side + 1, 1);
      Jet          jet;
      if (!jetAt(rail.s[side], U, V, jet))
        return 0.0;
      const gp_XYZ n = jet.n.XYZ(), c1 = (jet.du * U1 + jet.dv * V1).XYZ();
      const gp_XYZ c2 = (jet.duu * (U1 * U1) + jet.duv * (2 * U1 * V1) + jet.dvv * (V1 * V1)).XYZ();
      if (!(P1.SquareModulus() > 0 && c1.SquareModulus() > 0))
        return 0.0;
      return std::abs(n.Dot(P2) / P1.SquareModulus() - n.Dot(c2) / c1.SquareModulus()) * rail.scale;
    };
    // Never in the last refinements, which the positions keep.
    const bool g2 = rail.law == Law::G2 && bendRounds < limit::bendRefinements && refinement < 8;
    // The error between the sites: section poles (homogeneous, a bound on the
    // section for equal weights), and the feet mapped back onto their supports.
    std::vector<double> refine;
    bool                bent       = false;
    double              worstError = 0;
    for (int i = 0; i + 1 < count; ++i)
      for (double f : {0.25, 0.5, 0.75})
      {
        const double u   = sites[i] + (sites[i + 1] - sites[i]) * f;
        const auto   got = evaluate(u), want = row(u);
        double       error = 0;
        for (int j = 0; j < S; ++j)
        {
          const double wg = got[4 * j + 3], ww = want[4 * j + 3];
          if (!(wg > 0))
          {
            error = HUGE_VAL;
            break;
          }
          const gp_Pnt pg(got[4 * j] / wg, got[4 * j + 1] / wg, got[4 * j + 2] / wg);
          const gp_Pnt pw(want[4 * j] / ww, want[4 * j + 1] / ww, want[4 * j + 2] / ww);
          error = std::max(error, pg.Distance(pw) + std::abs(wg - ww) * rail.scale);
        }
        for (int side = 0; side < 2; ++side)
        {
          const gp_Pnt foot =
            footAt(rail.s[side], got[4 * S + 2 * side], got[4 * S + 2 * side + 1]);
          const gp_Pnt footWant =
            footAt(rail.s[side], want[4 * S + 2 * side], want[4 * S + 2 * side + 1]);
          error = std::max(error, foot.Distance(footWant));
        }
        worstError = std::max(worstError, error);
        bool split = error > tolerance;
        if (!split && g2)
          for (int side = 0; side < 2 && !split; ++side)
            if (!rail.s[side].cliffCurve && alongBend(u, side) > limit::bend)
              split = bent = true;
        if (split)
        {
          const auto right = knots.upper_bound(u);
          refine.push_back((std::prev(right)->first + right->first) / 2);
        }
      }
    // Refined for the bend alone, a bounded number of times: the positions
    // decide the fit.
    if (bent)
      ++bendRounds;
    if (refine.empty())
    {
      // The surface and the curves, from the coefficients.
      std::vector<double> distinct;
      std::vector<int>    multiplicity;
      for (const auto& [u, m] : knots)
      {
        distinct.push_back(u);
        multiplicity.push_back(m);
      }
      NCollection_Array1<double> uKnots(1, int(distinct.size()));
      NCollection_Array1<int>    uMults(1, int(distinct.size()));
      for (size_t k = 0; k < distinct.size(); ++k)
      {
        uKnots(int(k) + 1) = distinct[k];
        uMults(int(k) + 1) = multiplicity[k];
      }
      NCollection_Array1<double> vKnots(1, int(reference.knots.size()));
      NCollection_Array1<int>    vMults(1, int(reference.knots.size()));
      for (size_t k = 0; k < reference.knots.size(); ++k)
      {
        vKnots(int(k) + 1) = reference.knots[k];
        vMults(int(k) + 1) = reference.multiplicities[k];
      }
      NCollection_Array2<gp_Pnt> poles(1, count, 1, S);
      NCollection_Array2<double> weights(1, count, 1, S);
      bool                       rational = false;
      for (int i = 0; i < count; ++i)
        for (int j = 0; j < S; ++j)
        {
          const double w = coefficients[4 * j + 3][i];
          if (!(w > 0))
            fail(rail.op, no_solution, "the blend's weights do not stay positive", {rail.edge});
          rational |= std::abs(w - 1) > limit::degenerate;
          poles(i + 1, j + 1)   = gp_Pnt(coefficients[4 * j][i] / w,
                                       coefficients[4 * j + 1][i] / w,
                                       coefficients[4 * j + 2][i] / w);
          weights(i + 1, j + 1) = w;
        }
      Fitted out;
      out.blend = rational ? new Geom_BSplineSurface(poles,
                                                     weights,
                                                     uKnots,
                                                     vKnots,
                                                     uMults,
                                                     vMults,
                                                     kDegree,
                                                     reference.degree)
                           : new Geom_BSplineSurface(poles,
                                                     uKnots,
                                                     vKnots,
                                                     uMults,
                                                     vMults,
                                                     kDegree,
                                                     reference.degree);
      for (int side = 0; side < 2; ++side)
      {
        NCollection_Array1<gp_Pnt2d> contact(1, count);
        for (int i = 0; i < count; ++i)
          contact(i + 1) =
            gp_Pnt2d(coefficients[4 * S + 2 * side][i], coefficients[4 * S + 2 * side + 1][i]);
        out.contact[side] = new Geom2d_BSplineCurve(contact, uKnots, uMults, kDegree);
      }
      return out;
    }
    // Breaks the fit only halves its error towards (de Boor) are looked for once a smooth fit
    // should have converged: a foot's short excursion past a support's end between the march's
    // stations, and where a G2 section's measured law changes form, a handle or the depth point's
    // derivative at its limit (curvedG2): C0 in the poles there.
    if (refinement == 2)
    {
      if (raise(breaksOf(rail, stations, true)))
        continue;
    }
    if (refinement == 3 && !contactPilot && rail.law == Law::G2)
    {
      const auto regime = [&](double u) {
        Station s;
        Section section;
        exact(u, s, section);
        return section.regime;
      };
      constexpr int         scan   = 128;
      int                   before = regime(0);
      std::map<double, int> changes;
      for (int k = 1; k <= scan; ++k)
      {
        const int now = regime(double(k) / scan);
        if (now != before)
        {
          double low = double(k - 1) / scan, high = double(k) / scan;
          for (int i = 0; i < 60 && high - low > limit::degenerate; ++i)
          {
            const double middle                     = (low + high) / 2;
            (regime(middle) == before ? low : high) = middle;
          }
          addBreak(changes, (low + high) / 2, kDegree);
        }
        before = now;
      }
      if (raise(changes))
        continue;
    }
    // Halving the spans takes a smooth blend's error down some 64 times each
    // (degree 5). An error still a thousand tolerances that does not halve in
    // three refinements is a jump in the data (a contact that folds); near the
    // tolerance a slow tail goes on.
    worst.push_back(worstError);
    if (worst.size() > 3 && worstError > 1e3 * tolerance
        && worstError > worst[worst.size() - 4] / 2)
      break;
    for (double u : refine)
      knots.emplace(u, 1);
  }
  fail(rail.op, no_solution, "the blend surface cannot meet the length tolerance", {rail.edge});
}

//! A contact's speed the way the ball rolls at the fit's u: its foot's velocity along the direction
//! across both supports' normals at the contacts, taken the way the edge runs there (as a station's
//! is, ballAt). Not the edge's direction itself: a contact far from the edge (in an acute wedge)
//! moves across it, down a leaning wall, without running back.
//! Through adaptors, which keep the span last evaluated: a thousand and more evaluations each.
struct Speed
{
  const Rail&         rail;
  int                 side;
  Geom2dAdaptor_Curve contact, other;
  GeomAdaptor_Surface surface;

  Speed(const Rail& r, const Fitted& fitted, int s)
      : rail(r),
        side(s),
        contact(fitted.contact[s]),
        other(fitted.contact[1 - s]),
        surface(r.s[s].surface)
  {
  }

  double operator()(double u) const
  {
    gp_Pnt2d uv;
    gp_Vec2d duv;
    contact.D1(u, uv, duv);
    gp_Pnt foot, at;
    gp_Vec su, sv, d;
    surface.D1(uv.X(), uv.Y(), foot, su, sv);
    rail.curve.D1(rail.t(u), at, d);
    // Where the supports are tangent at the feet (a crease's end) the ball rolls along the edge.
    Jet            j[2];
    const gp_Pnt2d there = other.Value(u);
    gp_Vec         along = d;
    if (jetAt(rail.s[side], uv.X(), uv.Y(), j[side])
        && jetAt(rail.s[1 - side], there.X(), there.Y(), j[1 - side])
        && j[0].n.Crossed(j[1].n).Magnitude() > tol::coarseAngle)
      along = j[0].n.Crossed(j[1].n) * (j[0].n.Crossed(j[1].n).Dot(d) < 0 ? -1.0 : 1.0);
    if (!(along.Magnitude() > limit::degenerate))
      return 0.0;
    return (su * duv.X() + sv * duv.Y()).Dot(along.Normalized());
  }
};

//! The places the speed is read at between u0 and u1: within each knot span of the fitted contact
//! (where its own refinement put them, as finely as the data asks: a fold however short on however
//! long an edge), 2 (kDegree - 1) places, two more than a span's speed (of degree kDegree - 1 on a
//! flat cap) can have turns.
static std::vector<double> sites(const Fitted& fitted, int side, double u0, double u1)
{
  const auto&         contact = fitted.contact[side];
  std::vector<double> knots   = {u0};
  for (int i = 1; i <= contact->NbKnots(); ++i)
    if (contact->Knot(i) > u0 && contact->Knot(i) < u1)
      knots.push_back(contact->Knot(i));
  knots.push_back(u1);
  constexpr int       perSpan = 2 * (kDegree - 1);
  std::vector<double> out;
  for (size_t i = 0; i + 1 < knots.size(); ++i)
    for (int j = 0; j < perSpan; ++j)
      out.push_back(knots[i] + (knots[i + 1] - knots[i]) * j / perSpan);
  out.push_back(u1);
  return out;
}

//! Where a contact's speed along the edge dips, and how low: each strict least of the samples,
//! refined by golden-section search to the true least, and each end where it is running back
//! (checkFold).
static std::vector<std::pair<double, double>> dips(const Rail& rail, const Fitted& fitted, int side)
{
  const double        u0 = rail.u(rail.first), u1 = rail.u(rail.last);
  const Speed         speed(rail, fitted, side);
  const auto          at      = sites(fitted, side, u0, u1);
  const int           samples = int(at.size()) - 1;
  std::vector<double> v(samples + 1);
  for (int i = 0; i <= samples; ++i)
    v[i] = speed(at[i]);
  std::vector<std::pair<double, double>> out;
  if (v.front() <= 0 && v[1] >= v[0])
    out.emplace_back(v.front(), u0);
  if (v.back() <= 0 && v[samples - 1] >= v[samples])
    out.emplace_back(v.back(), u1);
  for (int i = 1; i < samples; ++i)
  {
    // A strict least of the samples; a steady speed has none.
    if (!(v[i] <= v[i - 1] && v[i] <= v[i + 1] && (v[i] < v[i - 1] || v[i] < v[i + 1])))
      continue;
    double                    lo = at[i - 1], hi = at[i + 1];
    std::pair<double, double> least  = {v[i], at[i]};
    constexpr double          golden = 0.6180339887498949;
    double a = hi - golden * (hi - lo), b = lo + golden * (hi - lo), fa = speed(a), fb = speed(b);
    for (int k = 0; k < 60 && hi - lo > limit::degenerate; ++k)
    {
      if (fa < fb)
      {
        hi = b;
        b  = a;
        fb = fa;
        a  = hi - golden * (hi - lo);
        fa = speed(a);
      }
      else
      {
        lo = a;
        a  = b;
        fa = fb;
        b  = lo + golden * (hi - lo);
        fb = speed(b);
      }
      if (fa < least.first)
        least = {fa, a};
      if (fb < least.first)
        least = {fb, b};
    }
    out.push_back(least);
  }
  return out;
}

const std::vector<std::pair<double, double>>& Dips::of(int side)
{
  if (!found[side])
    found[side] = dips(rail, fitted, side);
  return *found[side];
}

//! Whether a contact runs back anywhere along the edge.
static bool runsBack(Dips& dips, int side)
{
  const auto& found = dips.of(side);
  return std::any_of(found.begin(), found.end(), [](const auto& d) { return d.first <= 0; });
}

bool runsBack(Dips& dips)
{
  for (int side = 0; side < 2; ++side)
    if (!dips.rail.s[side].cliffCurve && dips.rail.s[side].trim.IsNull() && runsBack(dips, side))
      return true;
  return false;
}

void checkFold(Dips& dips)
{
  const Rail& rail = dips.rail;
  for (int side = 0; side < 2; ++side)
  {
    if (rail.s[side].cliffCurve || !rail.s[side].trim.IsNull())
      continue; // a cliff's contact is its edge; a notch's is trimmed away
    if (!runsBack(dips, side))
      continue;
    // Where a contact folds, the blend's own swallowtail is to be trimmed, which the blend does
    // not do yet; past the faces' bound nothing can be made.
    const bool flat = GeomAdaptor_Surface(rail.s[side].surface).GetType() == GeomAbs_Plane;
    fail(rail.op,
         unsupported_geometry,
         flat ? "a contact on the planar cap beside the edge folds at this radius, which v2 does "
                "not trim yet"
              : "a contact folds where the radius is beyond a face's convex radius, which v2 does "
                "not trim yet",
         {rail.edge},
         FailureKind::Curvature);
  }
}

std::vector<Fold> capFolds(Dips& dips, std::vector<Fold>* open)
{
  const Rail&   rail   = dips.rail;
  const Fitted& fitted = dips.fitted;
  const double  u0 = rail.u(rail.first), u1 = rail.u(rail.last);
  for (int side = 0; side < 2; ++side)
  {
    // A flat cap, whatever its surface's kind (a plane scaled unevenly is a spline): only a plane
    // is patched (a cliff is not patched either).
    if (rail.s[side].cliffCurve
        || !GeomLib_IsPlanarSurface(rail.s[side].surface, tol::length).IsPlanar())
      continue;
    if (!runsBack(dips, side))
      continue;
    // The other contact rolls on: only the cap's folds.
    if (runsBack(dips, 1 - side))
      return {};
    const Speed speed(rail, fitted, side);
    const auto  spacing = sites(fitted, side, u0, u1);
    // Each place the contact runs back lies in a loop where it crosses itself, the narrowest
    // round it; the loops apart. The run back lies between two cusps, where the speed turns;
    // the loop's ends are where the contact before the first crosses the contact after the
    // second, looked for in reaches about the cusps that double from the run back's own
    // length: found from the run back's own geometry, whatever the edge's length.
    const auto cusp = [&](double at, int way) -> std::optional<double> {
      // The sites past the place in turn, to the first where the speed is positive again;
      // the turn between, by bisection to the parameter's resolution.
      std::vector<double> ahead;
      for (double u : spacing)
        if ((u - at) * way > 0)
          ahead.push_back(u);
      if (way < 0)
        std::reverse(ahead.begin(), ahead.end());
      double inside = at;
      for (double u : ahead)
      {
        if (speed(u) > 0)
        {
          double lo = inside, hi = u;
          for (int k = 0; k < 60 && std::abs(hi - lo) > tol::parameter; ++k)
            (speed((lo + hi) / 2) > 0 ? hi : lo) = (lo + hi) / 2;
          return (lo + hi) / 2;
        }
        inside = u;
      }
      return std::nullopt;
    };
    std::vector<Fold> out;
    for (const auto& [least, at] : dips.of(side))
    {
      if (least > 0)
        continue;
      std::optional<Fold> best;
      const auto          first = cusp(at, -1), last = cusp(at, 1);
      if (first && last)
        for (double reach = std::max(*last - *first, tol::parameter); !best; reach *= 2)
        {
          const double a0 = std::max(u0, *first - reach), b1 = std::min(u1, *last + reach);
          const Geom2dAPI_InterCurveCurve crossing(
            new Geom2d_TrimmedCurve(fitted.contact[side], a0, *first),
            new Geom2d_TrimmedCurve(fitted.contact[side], *last, b1),
            tol::parameter);
          for (int i = 1; i <= crossing.NbPoints(); ++i)
          {
            const auto&  point = crossing.Intersector().Point(i);
            const double a = point.ParamOnFirst(), b = point.ParamOnSecond();
            if (!best || b - a < best->to - best->from)
              best = Fold{side, a, b};
          }
          if (a0 <= u0 && b1 >= u1)
            break;
        }
      if (!best)
      {
        // No crossing within the piece: the loop may run on across an end into the next
        // piece of the chain, which the caller asks.
        if (!open)
          return {};
        open->push_back(Fold{side, at, at});
        continue;
      }
      if (std::none_of(out.begin(), out.end(), [&](const Fold& f) {
            return f.from == best->from && f.to == best->to;
          }))
        out.push_back(*best);
    }
    std::sort(out.begin(), out.end(), [](const Fold& x, const Fold& y) { return x.from < y.from; });
    for (size_t i = 1; i < out.size(); ++i)
      if (out[i].from <= out[i - 1].to)
        return {};
    return out;
  }
  return {};
}

} // namespace BRepExt::blend2::detail
