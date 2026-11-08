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
#include <GeomConvert_CompCurveToBSplineCurve.hxx>

namespace BRepExt::blending::detail
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

std::map<double, int> breaksOf(const Rail& rail, Stations& stations, bool sought, bool centre)
{
  const int  order  = !centre && rail.law == Law::G2 ? 2 : 1;
  const auto needed = [&](int degree, int multiplicity) {
    return kDegree - std::max(0, degree - multiplicity - order);
  };
  std::map<double, int> breaks;
  // The edge's knots: each station is the ball whose section's plane passes through the edge's
  // point there (ballAt), so the stations are as smooth as the edge's points, not its tangent:
  // a C1 edge (OpenCascade's approximated meetings, degree 8 with knots of multiplicity 7) makes
  // the blend C1 there, not C0.
  if (!rail.centreParameter && rail.curve.GetType() == GeomAbs_BSplineCurve)
  {
    const auto curve = rail.curve.BSpline();
    for (int k = 2; k < curve->NbKnots(); ++k)
      addBreak(breaks,
               rail.u(curve->Knot(k)),
               kDegree - std::max(0, curve->Degree() - curve->Multiplicity(k)));
  }

  // The size law's places (Rail::sizes): a cubic spline is twice differentiable there and the
  // stations follow its value, so the blend is too, as at an edge's knots.
  for (size_t k = 1; k + 1 < rail.sizes.places.size(); ++k)
    addBreak(breaks, rail.u(rail.sizes.places[k]), kDegree - 2);

  // The supports' knot lines, each as a parameter (0 for u, 1 for v), a value and what it needs.
  struct Line
  {
    int    parameter;
    double value;
    int    multiplicity;
  };

  // A centre-based station plane also follows normals evaluated on the
  // source edge. Its data can cross a support boundary before either contact.
  const int paths = centre && rail.byCentre && rail.spine.IsNull() ? 2 : 1;
  for (int path = 0; path < paths; ++path)
  for (int side = 0; side < 2; ++side)
  {
    const Side& support = path && rail.parameterSupports ? (*rail.parameterSupports)[side]
                                                        : rail.s[side];
    std::vector<Line> lines;
    if (rail.s[side].cliffCurve)
      continue; // its feet run along an edge, not across a face's knots
    const Handle(Geom_Surface) surface    = basisOf(support.surface);
    const auto                 curveLines = [&](const Handle(Geom_Curve)& curve, int parameter) {
      if (const auto spline = Handle(Geom_BSplineCurve)::DownCast(basisOf(curve)))
        for (int k = 2; k < spline->NbKnots(); ++k)
          lines.push_back(
            {parameter, spline->Knot(k), needed(spline->Degree(), spline->Multiplicity(k))});
    };
    if (const auto spline = Handle(Geom_BSplineSurface)::DownCast(surface))
    {
      // Within a knot window the section is C-infinity in the support's coordinate and the
      // transition's ends are C2 (supportKnotWindows); the feet stay C1 at the knot.
      const auto knotLine = [&](int parameter, int k, int degree, int multiplicity) {
        const double value = parameter ? spline->VKnot(k) : spline->UKnot(k);
        if (!path && rail.knotWindows)
          for (const auto& window : *rail.knotWindows)
            if (window.side == side && window.parameter == parameter && window.knot == k)
            {
              for (double x : {-window.half, 0.0, window.half})
                lines.push_back({parameter, value + x, kDegree - 2});
              return;
            }
        lines.push_back({parameter, value, needed(degree, multiplicity)});
      };
      for (int k = 2; k < spline->NbUKnots(); ++k)
        knotLine(0, k, spline->UDegree(), spline->UMultiplicity(k));
      for (int k = 2; k < spline->NbVKnots(); ++k)
        knotLine(1, k, spline->VDegree(), spline->VMultiplicity(k));
    }
    else if (const auto extrusion = Handle(Geom_SurfaceOfLinearExtrusion)::DownCast(surface))
      curveLines(extrusion->BasisCurve(), 0);
    else if (const auto revolution = Handle(Geom_SurfaceOfRevolution)::DownCast(surface))
      curveLines(revolution->BasisCurve(), 1);
    // Where the surface ends and runs on along its tangent (runsOn), it is only G1: a break
    // where a ball's foot passes it, as marched (a scan asks for stations of its own, so a fit
    // that needs no such break asks for none).
    const auto ends = runsOn(support.surface);
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
    double first, last;
    const auto pcurve = path ? BRep_Tool::CurveOnSurface(rail.edge, support.face, first, last)
                             : Handle(Geom2d_Curve)();
    if (path && pcurve.IsNull())
      fail(rail.op, invalid_input, "the station plane has no support parameter curve", {rail.edge});
    const auto    footAt = [&](double u, Station& s) {
      if (path)
      {
        const auto uv = pcurve->Value(rail.t(u));
        s.uv[side][0] = uv.X();
        s.uv[side][1] = uv.Y();
      }
      else if (!stations.at(rail.t(u), s))
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

std::shared_ptr<const std::vector<Rail::KnotWindow>> supportKnotWindows(const Rail& rail,
                                                                        Stations&   stations)
{
  if (rail.law != Law::G2)
    return {};
  // The section is the continued spans' P_a + S(z) (P_b - P_a), S the quintic activation and
  // z = (x + h) / 2h in the support's coordinate x from the knot. On either half it departs
  // from the exact law's own branch by S(z) D(x), D = P_b - P_a, which vanishes at the knot,
  // where the spans agree through order two (derived bounds, the constants maximised over z):
  // - with |D(x)| <= |x| m (the poles, and the section's curvature at the foot, which is the
  //   contact's across the blend), the departure is at most c m h, c = max r S((1 - r) / 2),
  //   the maximiser solving 4 - 15 r + 20 r^3 - 9 r^5 = 0;
  // - the feet, on spans agreeing through order two, part as |D(x)| <= |x|^3 m; the normal
  //   curvature along the contact then departs by at most k m h / v^2, v the feet's speed in x
  //   and k = max |(S(z) (2z - 1)^3)''| / 4 over z in [0, 1/2].
  // The poles bound the section (homogeneous, as the fit measures them). The window is the
  // widest that keeps the section within its limit, less the fit's tolerance, and the contact
  // within the fit's own limit on its curvature (limit::bend); its departures across and along
  // add.
  constexpr double spread       = 0.070552878102991504; // c
  constexpr double bending      = 0.32818934415354306;  // k
  const double     sectionShare = limit::section * rail.scale - fitTolerance(rail);
  const double     bendShare    = limit::bend / rail.scale;
  if (!(sectionShare > 0))
    return {};
  auto       windows = std::make_shared<std::vector<Rail::KnotWindow>>();
  const auto station = [&](double u, Station& s) {
    if (!stations.at(rail.t(u), s))
      fail(rail.op, no_solution, "the ball cannot roll the length of the edge at this size",
           {rail.edge}, FailureKind::SupportContact);
  };
  for (int side = 0; side < 2; ++side)
  {
    if (rail.s[side].cliffCurve)
      continue;
    const auto spline = Handle(Geom_BSplineSurface)::DownCast(basisOf(rail.s[side].surface));
    if (spline.IsNull())
      continue;
    for (int parameter = 0; parameter < 2; ++parameter)
    {
      const int count = parameter ? spline->NbVKnots() : spline->NbUKnots();
      const int degree = parameter ? spline->VDegree() : spline->UDegree();
      const auto knot = [&](int k) { return parameter ? spline->VKnot(k) : spline->UKnot(k); };
      for (int k = 2; k < count; ++k)
      {
        if (degree - (parameter ? spline->VMultiplicity(k) : spline->UMultiplicity(k)) != 2)
          continue;
        const double value = knot(k);
        // Both continued spans must stay within their neighbours, and neighbouring windows apart.
        const double widest = std::min(value - knot(k - 1), knot(k + 1) - value) / 2;
        const auto   coordinate = [&](double u) {
          Station s;
          station(u, s);
          return s.uv[side][parameter] - value;
        };
        constexpr int scan   = 256;
        double        before = coordinate(0);
        for (int n = 1; n <= scan; ++n)
        {
          const double after = coordinate(double(n) / scan);
          if ((before > 0) == (after > 0))
          {
            before = after;
            continue;
          }
          double low = double(n - 1) / scan, high = double(n) / scan;
          for (int i = 0; i < 60 && high - low > limit::degenerate; ++i)
          {
            const double middle = (low + high) / 2;
            ((coordinate(middle) > 0) == (before > 0) ? low : high) = middle;
          }
          before = after;
          const double event = (low + high) / 2;
          // The rate of the coordinate along the fit's u, to place the probes.
          const double step = std::min({1e-6, event, 1 - event});
          if (!(step > 0))
            continue;
          const double rate = (coordinate(event + step) - coordinate(event - step)) / (2 * step);
          if (!(std::abs(rate) > limit::degenerate))
            continue;
          // The window's ratio of departure to its limit, from both spans at x = -h and h.
          const auto ratio = [&](double h) {
            double  worst = 0, across = 0, along = 0;
            gp_Pnt  feet[2];
            double  places[2];
            for (int end = 0; end < 2; ++end)
            {
              const double u = event + (end ? h : -h) / rate;
              if (!(u > 0 && u < 1))
                return HUGE_VAL;
              Station s;
              station(u, s);
              const double x = s.uv[side][parameter] - value;
              Section      pair[2];
              if (!(std::abs(x) > 0)
                  || !sectionsBesideSupportKnot(rail, s, side, parameter, k, pair))
                return HUGE_VAL;
              double poles = 0;
              for (size_t j = 0; j < pair[0].poles.size(); ++j)
                poles = std::max(poles,
                                 pair[0].poles[j].Distance(pair[1].poles[j])
                                   + std::abs(pair[0].weights[j] - pair[1].weights[j])
                                       * rail.scale);
              double bend[2];
              for (int branch = 0; branch < 2; ++branch)
              {
                const auto curve = sectionCurve(pair[branch]);
                gp_Pnt     p;
                gp_Vec     d1, d2;
                curve->D2(side ? curve->LastParameter() : curve->FirstParameter(), p, d1, d2);
                bend[branch] = d1.Crossed(d2).Magnitude() / std::pow(d1.Magnitude(), 3);
              }
              const size_t foot = side ? pair[0].poles.size() - 1 : 0;
              const double apart =
                std::abs(gp_Vec(pair[0].poles[foot], pair[1].poles[foot]).Dot(s.normal[side]));
              feet[end]   = s.foot[side];
              places[end] = x;
              along       = std::max(along, apart / std::pow(std::abs(x), 3));
              const double scaled = spread * h / std::abs(x);
              worst  = std::max(worst, scaled * poles / sectionShare);
              across = std::max(across, scaled * std::abs(bend[0] - bend[1]));
            }
            const double speed = feet[0].Distance(feet[1]) / std::abs(places[1] - places[0]);
            if (!(speed > 0))
              return HUGE_VAL;
            return std::max(worst, (across + bending * along * h / (speed * speed)) / bendShare);
          };
          // At the knot the window's middle departs by |D(0)| / 2 however narrow it is. The
          // spans agree there unless the foot is on the support's extension past its edge,
          // which no polynomial span continues: then the exact law keeps its break.
          {
            Station s;
            station(event, s);
            Section pair[2];
            if (!sectionsBesideSupportKnot(rail, s, side, parameter, k, pair))
              continue;
            double apart = 0;
            for (size_t j = 0; j < pair[0].poles.size(); ++j)
              apart = std::max(apart, pair[0].poles[j].Distance(pair[1].poles[j]));
            if (!(apart / 2 <= sectionShare))
              continue;
          }
          double h = widest;
          if (!(ratio(h) <= 1))
          {
            // The departure grows with h: bisect, in proportion, down to where it meets the limit.
            double low = widest * 1e-9, high = widest;
            for (int i = 0; i < 48 && high > low * (1 + 1e-3); ++i)
            {
              const double middle = std::sqrt(low * high);
              (ratio(middle) <= 1 ? low : high) = middle;
            }
            h = low;
          }
          windows->push_back({side, parameter, k, value, h});
        }
      }
    }
  }
  if (windows->empty())
    return {};
  return windows;
}

double fitTolerance(const Rail& rail)
{
  // The length tolerance, and for a blend smaller than a millimetre that share of its size, so
  // that a small blend is as true to its law as a large one; never finer than the rounding of
  // its points where it sits (tolerance::resolved).
  return std::max(tol::length * std::min(1.0, rail.scale) * (rail.law == Law::G2 ? 0.1 : 1),
                  tol::rounding(rail.magnitude + rail.scale));
}

Fitted fit(const Rail&                      rail,
           Stations&                        stations,
           bool                             contactPilot,
           const std::vector<SectionFrame>& frames,
           const Fitted*                    terminalReference,
           bool                             reparameterizedTerminal)
{
  const double          tolerance = fitTolerance(rail);
  std::map<double, int> knots     = {{0, kDegree + 1}, {1, kDegree + 1}};
  for (int k = 1; k < 8; ++k)
    knots[k / 8.0] = 1;
  auto naturalBreaks = breaksOf(rail, stations);
  if (rail.conicWindow)
    for (double target : {-rail.conicWindow->halfWidth, rail.conicWindow->halfWidth})
    {
      const auto value = [&](double u) {
        Station station;
        if (!stations.at(rail.t(u), station))
          fail(rail.op, no_solution, "a conic window has no contact station", {rail.edge});
        return conicPhase(rail, station) - target;
      };
      double before = value(0);
      for (int k = 1; k <= 256; ++k)
      {
        const double u = k / 256.0, after = value(u);
        if (before * after < 0 && std::abs(after - before) < M_PI)
        {
          double lo = (k - 1) / 256.0, hi = u;
          for (int j = 0; j < 48; ++j)
          {
            const double mid = (lo + hi) / 2;
            if (value(mid) * before > 0)
              lo = mid;
            else
              hi = mid;
          }
          addBreak(naturalBreaks, (lo + hi) / 2, kDegree - 1);
        }
        before = after;
      }
    }

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
    if (terminalReference && reparameterizedTerminal)
    {
      const auto& basis = terminalReference->blend;
      auto curve = sectionCurve(section);
      curve->IncreaseDegree(basis->VDegree());
      curve->InsertKnots(basis->VKnots(), basis->VMultiplicities(), tol::parameter, false);
      section.degree = curve->Degree();
      section.poles.clear(); section.weights.clear();
      section.knots.clear(); section.multiplicities.clear();
      for (int i = 1; i <= curve->NbPoles(); ++i)
      {
        section.poles.push_back(curve->Pole(i));
        section.weights.push_back(curve->Weight(i));
      }
      for (int i = 1; i <= curve->NbKnots(); ++i)
      {
        section.knots.push_back(curve->Knot(i));
        section.multiplicities.push_back(curve->Multiplicity(i));
      }
    }
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
    if (terminalReference && reparameterizedTerminal)
      return true;
    if (!terminalReference || frames.empty())
      return false;
    for (const auto& frame : frames)
    {
      const double distance =
        railLength(rail, t, frame.t);
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
  const auto referenceRow = [&](double u, int order) {
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
  // Chain rule for the old normalized edge coordinate h as a function
  // of normalized centre length u: (F o h)'' = F'' h'^2 + F' h''.
  const auto terminalRow = [&](double u, int order) {
    auto result = referenceRow(u, order);
    if (!reparameterizedTerminal || !order)
      return result;
    gp_Pnt p;
    gp_Vec d, dd;
    const auto& curve = rail.centreParameter->curve;
    const double t = u == 0 ? rail.from : rail.to;
    int span = 1;
    while (span + 1 < curve->NbKnots()
           && (u == 0 ? curve->Knot(span + 1) <= t : curve->Knot(span + 1) < t))
      ++span;
    curve->LocalD2(t, span, span + 1, p, d, dd);
    const double speed = d.Magnitude();
    const double length = rail.centreParameter->length(rail.to)
                          - rail.centreParameter->length(rail.from);
    const double rate = length / (speed * (rail.to - rail.from));
    const double acceleration = -length * length * d.Dot(dd)
      / (std::pow(speed, 4) * (rail.to - rail.from));
    const auto first = order == 2 ? referenceRow(u, 1) : std::vector<double>();
    for (size_t i = 0; i < result.size(); ++i)
      result[i] = order == 1 ? result[i] * rate
                           : result[i] * rate * rate + first[i] * acceleration;
    return result;
  };
  std::vector<double> worst;
  int                 bendRounds = 0;
  // Local refinement can keep converging in a short interval long after sixteen rounds.
  // The pole budget bounds the work; require a new representable knot at every ordinary
  // refinement instead of imposing a second, unrelated iteration budget.
  for (int refinement = 0; !rail.boundedFit || refinement < 16; ++refinement)
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
      if (raiseFitBreaks(knots, breaksOf(rail, stations, true)))
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
      if (raiseFitBreaks(knots, changes))
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
    const auto previousKnots = knots.size();
    for (double u : refine)
      knots.emplace(u, 1);
    if (!rail.boundedFit && knots.size() == previousKnots)
    {
      // Break discovery (2), regime discovery (3), and bend testing (<8) can still
      // change the next fit without an inserted knot. These are finite events.
      // Round 8 has no bend check; after it no scheduled transition remains.
      if (refinement < 8)
        continue;
      break;
    }
  }
  fail(rail.op, no_solution, "the blend surface cannot meet the length tolerance", {rail.edge});
}

//! A contact's speed the way the ball rolls at the fit's u: its foot's velocity along the direction
//! across both supports' normals at the contacts, taken the way the edge runs there (as a station's
//! is, ballAt); a chamfer's, along the edge (its sections lie across the edge). Not the edge's
//! direction itself for a ball: a contact far from the edge (in an acute wedge)
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
    rail.stationD1(rail.t(u), at, d);
    // Where the supports are tangent at the feet (a crease's end) the ball rolls along the edge.
    // A chamfer's section lies in the plane across the edge at its point (setbackAt): the way it
    // goes on is the edge's own, whatever the faces' normals at its feet.
    Jet            j[2];
    const gp_Pnt2d there = other.Value(u);
    gp_Vec         along = d;
    if (rail.law != Law::Chamfer && jetAt(rail.s[side], uv.X(), uv.Y(), j[side])
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
    // A run back is the contact moving against the way the balls roll between two places where
    // it stops: only where that stretch of the contact has a length (more than the places two
    // constructions must share) is it a fold of the contact. Where the fit's own parameter stalls
    // (a stationary point of the sheet's parameterisation, at which both contacts stop together
    // and go on the same way), the speed's least is a rounding of nothing, on no length.
    if (least.first <= 0)
    {
      const auto zero = [&](int way) {
        double inside = least.second, outside = inside;
        for (int k = i + way; k >= 0 && k <= samples; k += way)
          if (speed(at[k]) > 0)
          {
            outside = at[k];
            break;
          }
        if (outside == inside)
          return inside;
        for (int k = 0; k < 60 && std::abs(outside - inside) > limit::degenerate; ++k)
          (speed((inside + outside) / 2) > 0 ? outside : inside) = (inside + outside) / 2;
        return inside;
      };
      const double from = zero(-1), to = zero(1);
      double       length = 0;
      gp_Pnt       previous;
      for (int k = 0; k <= 16; ++k)
      {
        const gp_Pnt2d uv = fitted.contact[side]->Value(from + (to - from) * k / 16);
        const gp_Pnt   p  = rail.s[side].surface->Value(uv.X(), uv.Y());
        if (k > 0)
          length += p.Distance(previous);
        previous = p;
      }
      if (length <= rail.resolved(limit::together))
        continue;
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
    fail(
      rail.op,
      unsupported_geometry,
      flat
        ? "a contact on the planar cap beside the edge folds at this radius, which the blend does "
          "not trim yet"
        : "a contact folds where the radius is beyond a face's convex radius, which the blend does "
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
    // The contact crosses itself in the face's parameters, whatever its surface (a cliff's
    // contact is its edge, and is not patched).
    if (rail.s[side].cliffCurve)
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

} // namespace BRepExt::blending::detail

namespace BRepExt::blending::detail
{

double CentreParameter::length(double t) const
{
  const auto   it = std::upper_bound(knots.begin(), knots.end(), t);
  const size_t i =
    it == knots.begin() ? 0 : std::min(size_t(it - knots.begin() - 1), knots.size() - 2);
  if (t == knots[i])
    return lengths[i];
  if (t == knots[i + 1])
    return lengths[i + 1];
  GeomAdaptor_Curve adaptor(curve, std::min(knots[i], t), std::max(knots[i + 1], t));
  const double      delta =
    GCPnts_AbscissaPoint::Length(adaptor, std::min(knots[i], t), std::max(knots[i], t), tolerance);
  return lengths[i] + (t < knots[i] ? -delta : delta);
}

double CentreParameter::parameter(double distance) const
{
  const auto   it = std::upper_bound(lengths.begin(), lengths.end(), distance);
  const size_t i =
    it == lengths.begin() ? 0 : std::min(size_t(it - lengths.begin() - 1), lengths.size() - 2);
  if (distance == lengths[i])
    return knots[i];
  if (distance == lengths[i + 1])
    return knots[i + 1];
  GeomAdaptor_Curve    adaptor(curve,
                               distance < lengths.front() ? curve->FirstParameter() : knots[i],
                               distance > lengths.back() ? curve->LastParameter() : knots[i + 1]);
  GCPnts_AbscissaPoint along(tolerance, adaptor, distance - lengths[i], knots[i]);
  if (!along.IsDone())
    fail("blendEdges", no_solution, "the centre-path length cannot be inverted");
  return along.Parameter();
}

//! Fit only the station centres, retaining the original edge-coordinate breaks.
//! Arc length is subsequently integrated span by span; no polygonal length is used.
std::shared_ptr<CentreParameter> centreParameter(
  const Rail& rail, const std::function<gp_Pnt(double)>& pointAt,
  const std::vector<double>& joins, const std::map<double, int>& events)
{
  if (!rail.spine.IsNull() && !rail.spineArc)
  {
    const auto source = Handle(Geom_BSplineCurve)::DownCast(rail.spine);
    if (source.IsNull())
      fail(rail.op, no_solution, "the station carrier has no spline representation");
    const double rate = (rail.spineLast - rail.spineFirst)
                        / (rail.spineParameterLast - rail.spineParameterFirst);
    const double begin = source->FirstParameter(), period = source->LastParameter()-begin;
    const auto mapped = [&](double t) {
      return rail.spineFirst + (t-rail.spineParameterFirst)*rate;
    };
    std::vector<double> cuts{rail.from, rail.to};
    const double a = mapped(rail.from), b = mapped(rail.to);
    for (int turn = int(std::ceil((std::min(a,b)-begin)/period));
         turn <= int(std::floor((std::max(a,b)-begin)/period)); ++turn)
    {
      const double t = rail.spineParameterFirst+(begin+turn*period-rail.spineFirst)/rate;
      if (t > rail.from && t < rail.to) cuts.push_back(t);
    }
    std::sort(cuts.begin(), cuts.end());
    GeomConvert_CompCurveToBSplineCurve joined;
    for (size_t i=1; i<cuts.size(); ++i)
    {
      const double lo=cuts[i-1], hi=cuts[i];
      const double turn=std::floor((mapped((lo+hi)/2)-begin)/period);
      const double p=mapped(lo)-turn*period, q=mapped(hi)-turn*period;
      auto part=Handle(Geom_BSplineCurve)::DownCast(source->Copy());
      part->Segment(std::min(p,q),std::max(p,q));
      if (rate<0) part->Reverse();
      NCollection_Array1<double> knots(part->Knots());
      BSplCLib::Reparametrize(lo,hi,knots);
      part->SetKnots(knots);
      if (!joined.Add(part,rail.resolved(limit::solved),true,false,part->Degree()))
        fail(rail.op,no_solution,"the station carrier circuit does not join");
    }
    auto out=std::make_shared<CentreParameter>();
    out->curve=joined.BSplineCurve();
    out->tolerance=rail.resolved(limit::solved);
    GeomAdaptor_Curve adaptor(out->curve);
    double length=0;
    for (int k=1; k<=out->curve->NbKnots(); ++k)
    {
      if (k>1) length+=GCPnts_AbscissaPoint::Length(adaptor,out->curve->Knot(k-1),
                                                  out->curve->Knot(k),out->tolerance);
      out->knots.push_back(out->curve->Knot(k));
      out->lengths.push_back(length);
    }
    return out;
  }
  std::map<double, int> knots{{0, kDegree + 1}, {1, kDegree + 1}};
  for (int k = 1; k < 8; ++k)
    knots[k / 8.0] = 1;
  if (rail.curve.GetType() == GeomAbs_BSplineCurve)
  {
    const auto curve = rail.curve.BSpline();
    for (int k = 2; k < curve->NbKnots(); ++k)
      addBreak(knots,
               rail.u(curve->Knot(k)),
               kDegree - std::max(0, curve->Degree() - curve->Multiplicity(k)));
  }
  for (const auto& [t, multiplicity] : events)
    addBreak(knots, rail.u(t), multiplicity);
  // A support change shares the centre, but need not share its tangent.
  // Retain that C0 event instead of approximating it by ever shorter smooth spans.
  for (double t : joins)
    addBreak(knots, rail.u(t), kDegree);
  const double tolerance = rail.resolved(limit::solved);
  const auto exact = [&](double u) { return pointAt(rail.t(u)); };
  const gp_Pnt origin = exact(0);
  for (;;)
  {
    const auto flat = flatten(knots);
    const int  n    = int(flat.size()) - kDegree - 1;
    NCollection_Array1<double> full(1, int(flat.size())), parameters(1, n);
    NCollection_Array1<int>    orders(1, n);
    std::vector<double>        poles(3 * n), sites(n);
    for (size_t i = 0; i < flat.size(); ++i)
      full(int(i) + 1) = flat[i];
    for (int i = 0; i < n; ++i)
    {
      double u = 0;
      for (int j = 1; j <= kDegree; ++j)
        u += flat[i + j] / kDegree;
      sites[i]          = u;
      parameters(i + 1) = u;
      orders(i + 1)     = 0;
      const auto p      = exact(u).XYZ() - origin.XYZ();
      for (int c = 0; c < 3; ++c)
        poles[3 * i + c] = p.Coord(c + 1);
    }
    int singular = 0;
    BSplCLib::Interpolate(kDegree, full, parameters, orders, 3, poles[0], singular);
    if (singular)
      fail(rail.op, no_solution, "the centre-path interpolation is singular", {rail.edge});
    NCollection_Array1<gp_Pnt> points(1, n);
    for (int i = 0; i < n; ++i)
      points(i + 1) =
        gp_Pnt(origin.XYZ() + gp_XYZ(poles[3 * i], poles[3 * i + 1], poles[3 * i + 2]));
    NCollection_Array1<double> distinct(1, int(knots.size()));
    NCollection_Array1<int>    mult(1, int(knots.size()));
    int                        j = 1;
    for (const auto& [u, m] : knots)
    {
      distinct(j) = rail.t(u);
      mult(j++)   = m;
    }
    Handle(Geom_BSplineCurve) curve = new Geom_BSplineCurve(points, distinct, mult, kDegree);
    std::vector<double>       refine;
    for (int i = 0; i + 1 < n; ++i)
      for (double f : {.25, .5, .75})
      {
        const double u = sites[i] + (sites[i + 1] - sites[i]) * f;
        if (curve->Value(rail.t(u)).Distance(exact(u)) > tolerance)
        {
          const auto next = knots.upper_bound(u);
          refine.push_back((std::prev(next)->first + next->first) / 2);
        }
      }
    if (refine.empty())
    {
      auto out       = std::make_shared<CentreParameter>();
      out->curve     = curve;
      out->tolerance = tolerance;
      GeomAdaptor_Curve adaptor(curve);
      double            length = 0;
      for (int k = 1; k <= curve->NbKnots(); ++k)
      {
        if (k > 1)
          length +=
            GCPnts_AbscissaPoint::Length(adaptor, curve->Knot(k - 1), curve->Knot(k), tolerance);
        out->knots.push_back(curve->Knot(k));
        out->lengths.push_back(length);
      }
      return out;
    }
    bool changed = false;
    for (double u : refine)
      changed |= knots.emplace(u, 1).second;
    if (!changed)
      fail(rail.op, no_solution, "the centre-path fit cannot refine", {rail.edge});
  }
}

bool defectiveRailKnot(const Rail& rail, const Fitted& fitted)
{
  if (rail.curve.GetType() != GeomAbs_BSplineCurve)
    return false;
  const auto curve = rail.curve.BSpline();
  const auto surface = fitted.blend;
  for (int k = 2; k < surface->NbUKnots(); ++k)
  {
    if (surface->UMultiplicity(k) < surface->UDegree() - 1)
      continue;
    const double u         = surface->UKnot(k);
    bool         inherited = false;
    for (int j = 2; j < curve->NbKnots(); ++j)
      inherited |= std::abs(u - rail.u(curve->Knot(j))) <= tol::parameter;
    if (!inherited)
      continue;
    for (int j = 1; j < surface->NbVKnots(); ++j)
      for (int q = 0; q <= 18; ++q)
      {
        const double v = surface->VKnot(j) + (surface->VKnot(j + 1) - surface->VKnot(j)) * q / 18.;
        SurfaceJet   a{rail.op}, b{rail.op};
        surfaceD2OnUSpan(surface, u, v, k - 1, k, a);
        surfaceD2OnUSpan(surface, u, v, k, k + 1, b);
        a.orient(false);
        b.orient(false);
        if (bendDifference(a, b, a.v.Normalized()).first * rail.scale > limit::seamBend)
          return true;
      }
  }
  return false;
}

Fitted fit(Rail& rail, Stations& stations, bool contactPilot)
{
  rail.knotWindows = contactPilot ? nullptr : supportKnotWindows(rail, stations);
  Fitted original = fit(rail, stations, contactPilot, {});
  // Stations walked by the spine's arc length (Rail::spineArc) are already by the centres' path
  // length, which is what a centre parameter gives.
  if (contactPilot || rail.law != Law::G2 || rail.centreParameter || rail.parameterSupports
      || rail.spineArc || !rail.sizes.places.empty()
      || rail.s[0].cliffCurve || rail.s[1].cliffCurve
      || (rail.spine.IsNull() && !defectiveRailKnot(rail, original)))
    return original;
  try
  {
    Rail centreRail = rail;
    centreRail.centreInterpolation = true;
    // Preserve delivering stations unless centre fitting extends the source
    // parameter curve: only that domain requires stable endpoint evaluation.
    bool extrapolated = false;
    if (rail.byCentre && rail.spine.IsNull())
      for (int side = 0; side < 2; ++side)
      {
        double first, last;
        const auto pcurve = BRep_Tool::CurveOnSurface(rail.edge, rail.s[side].face, first, last);
        if (!pcurve.IsNull())
          extrapolated |= rail.from < pcurve->FirstParameter() || rail.to > pcurve->LastParameter();
      }
    std::unique_ptr<Stations> stableStations;
    Stations* centreStations = &stations;
    if (extrapolated)
    {
      stableStations = std::make_unique<Stations>(centreRail);
      Station anchor;
      if (!stations.at((rail.from + rail.to) / 2, anchor))
        fail(rail.op, no_solution, "the delivering centre path has no interior station", {rail.edge});
      stableStations->march(anchor);
      centreStations = stableStations.get();
    }
    rail.centreParameter = centreParameter(rail, [&](double t) {
      Station s;
      if (!centreStations->at(t, s))
        fail(rail.op, no_solution, "the centre path cannot be fitted", {rail.edge});
      return s.centre;
    });
    Fitted candidate = fit(rail, stations, false, {});
    // A change of parameter can require more interpolation coefficients even
    // when it removes an artificial curvature break. Positional and contact
    // checks in fit remain the acceptance conditions for this construction.
    return candidate;
  }
  catch (const Failure& failure)
  {
    rail.centreParameter.reset();
    if (failure.reason == std::string("cancelled"))
      throw;
    return original;
  }
}
} // namespace BRepExt::blending::detail
