// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Fillet v2: Exact sheets (cylinders and tori where the blend is one).

#include <BRepBlendExt_Internal.hxx>

namespace loft::blend2::detail {

Handle(Geom_BSplineCurve) sectionCurve(const Section& section) {
    const int count = int(section.poles.size());
    NCollection_Array1<gp_Pnt> poles(1, count);
    NCollection_Array1<double> weights(1, count);
    for (int i = 0; i < count; ++i) {
        poles(i + 1) = section.poles[i];
        weights(i + 1) = section.weights[i];
    }
    NCollection_Array1<double> knots(1, int(section.knots.size()));
    NCollection_Array1<int> multiplicities(1, int(section.knots.size()));
    for (size_t i = 0; i < section.knots.size(); ++i) {
        knots(int(i) + 1) = section.knots[i];
        multiplicities(int(i) + 1) = section.multiplicities[i];
    }
    return new Geom_BSplineCurve(poles, weights, knots, multiplicities, section.degree);
}

std::optional<Sheet> exactSheet(const Rail& rail, Stations& stations) {
    // A full round of one size along a straight edge is a cylinder too, its arc through the foot
    // on the face between.
    const bool full = rail.law == Law::Full;
    // A chamfer round a circle between faces about its axis is a cone: every section the same
    // line in the plane through the axis.
    const bool cone = rail.law == Law::Chamfer && rail.curve.GetType() == GeomAbs_Circle;
    if (!full && !cone && (rail.law != Law::Conic || std::abs(rail.rho - 0.5) > limit::equal ||
                  std::abs(rail.s[0].distance - rail.s[1].distance) > limit::equal * rail.scale))
        return {};
    const double tolerance = rail.resolved(limit::exactSheet);
    constexpr int samples = 9;
    std::vector<Station> at(samples);
    for (int k = 0; k < samples; ++k)
        if (!stations.at(rail.t(double(k) / (samples - 1)), at[k])) return {};
    const auto& first = at.front();
    const auto& last = at.back();
    const double r = full ? first.centre.Distance(first.foot[0]) : rail.s[0].distance;
    const auto frame = [&](const gp_Pnt& centre, const gp_Dir& axis, const gp_Pnt& towards) {
        const gp_Vec x = gp_Vec(centre, towards) - gp_Vec(axis) * gp_Vec(centre, towards).Dot(gp_Vec(axis));
        return fixedFrame(gp_Ax3(centre, axis, gp_Dir(x)));
    };
    // Every foot on its contact line at its place along the rail.
    const auto exact = [&](const Sheet& sheet) {
        for (int k = 0; k < samples; ++k)
            for (int side = 0; side < 2; ++side) {
                const gp_Pnt2d p = sheet.uv(sheet.rail(double(k) / (samples - 1)), sheet.across[side]);
                if (sheet.surface->Value(p.X(), p.Y()).Distance(at[k].foot[side]) > tolerance) return false;
            }
        return true;
    };
    if (cone) {
        // About the edge's axis, turning the way the rail runs; its reference circle through the
        // foot nearer the axis's start, so that the line climbs to the other one.
        const gp_Circ edge = rail.curve.Circle();
        gp_Dir axis = edge.Axis().Direction();
        const gp_Pnt O = edge.Location();
        if (gp_Vec(O, first.foot[0]).Crossed(gp_Vec(O, at[1].foot[0])).Dot(gp_Vec(axis)) < 0) axis.Reverse();
        double z[2], rho[2];
        for (int side = 0; side < 2; ++side) {
            const gp_Vec v(O, first.foot[side]);
            z[side] = v.Dot(gp_Vec(axis));
            rho[side] = (v - gp_Vec(axis) * z[side]).Magnitude();
        }
        const int low = z[0] < z[1] ? 0 : 1, high = 1 - low;
        const double rise = z[high] - z[low], spread = rho[high] - rho[low], length = std::hypot(rise, spread);
        // A flat ring or a cylinder is not a cone: fitted.
        if (!(rise > tol::angle * length) || !(std::abs(spread) > tol::angle * length)) return {};
        gp_Ax3 ax = frame(O.Translated(gp_Vec(axis) * z[low]), axis, first.foot[0]);
        double turn = 2 * M_PI, start = 0;
        if (!rail.closed) {
            const gp_Vec w(ax.Location(), last.foot[0]);
            turn = std::atan2(w.Dot(gp_Vec(ax.YDirection())), w.Dot(gp_Vec(ax.XDirection())));
            if (turn < 0) turn += 2 * M_PI;
            // The seam on the far side from the arc.
            ax = fixedFrame(ax.Rotated(gp_Ax1(ax.Location(), axis), turn / 2 + M_PI));
            start = M_PI - turn / 2;
        }
        Sheet sheet{new Geom_ConicalSurface(ax, std::asin(spread / length), rho[low]), false, start, turn, {}};
        sheet.across[low] = 0;
        sheet.across[high] = length;
        if (exact(sheet)) return sheet;
        return {};
    }
    if (!rail.closed && first.centre.Distance(last.centre) > tolerance) {
        // A cylinder about the line of the centres, its seam on the far side from the blend: its
        // x axis away from the middle of the feet, so that their angles lie about a half turn,
        // in a positive range, with no seam near the blend or what meets it.
        // Its axis the edge's own line where it is one (exact, as OpenCascade's are), else the
        // line of the centres.
        gp_Dir axis(gp_Vec(first.centre, last.centre));
        if (rail.curve.GetType() == GeomAbs_Line) {
            const gp_Dir along = rail.curve.Line().Direction();
            axis = along.Dot(axis) > 0 ? along : along.Reversed();
        }
        const gp_Vec middle = full ? gp_Vec(first.centre, first.middle)
                                   : gp_Vec(first.centre, first.foot[0]).Normalized() + gp_Vec(first.centre, first.foot[1]).Normalized();
        const gp_Ax3 ax = frame(first.centre, axis, first.centre.Translated(middle.Reversed()));
        const auto angleOf = [&](const gp_Pnt& foot) {
            const gp_Vec v(ax.Location(), foot);
            double a = std::atan2(v.Dot(gp_Vec(ax.YDirection())), v.Dot(gp_Vec(ax.XDirection())));
            return a < 0 ? a + 2 * M_PI : a;
        };
        Sheet sheet{new Geom_CylindricalSurface(ax, r), true, 0, first.centre.Distance(last.centre),
                    {angleOf(first.foot[0]), angleOf(first.foot[1])}};
        if (exact(sheet)) return sheet;
    }
    // A torus about the circle of the centres, u from the first centre and v about the tube.
    gce_MakeCirc circle(first.centre, at[samples / 2].centre, at[samples - 2].centre);
    if (!circle.IsDone()) return {};
    const gp_Circ spine = circle.Value();
    // Its axis the edge's own circle's where it is one, exactly.
    gp_Dir axis = spine.Axis().Direction();
    if (rail.curve.GetType() == GeomAbs_Circle) {
        const gp_Dir own = rail.curve.Circle().Axis().Direction();
        if (std::abs(std::abs(own.Dot(axis)) - 1) < tol::angle) axis = own.Dot(axis) > 0 ? own : own.Reversed();
    }
    // Turning the way the rail runs.
    if (gp_Vec(spine.Location(), first.centre).Crossed(gp_Vec(spine.Location(), at[1].centre)).Dot(gp_Vec(axis)) < 0) axis.Reverse();
    gp_Ax3 ax = frame(spine.Location(), axis, first.centre);
    double turn = 2 * M_PI, start = 0;
    if (!rail.closed) {
        const gp_Vec w(spine.Location(), last.centre);
        turn = std::atan2(w.Dot(gp_Vec(ax.YDirection())), w.Dot(gp_Vec(ax.XDirection())));
        if (turn < 0) turn += 2 * M_PI;
        // The seam on the far side from the arc: the x axis opposite its middle.
        ax = fixedFrame(ax.Rotated(gp_Ax1(spine.Location(), axis), turn / 2 + M_PI));
        start = M_PI - turn / 2;
    }
    // The tube's angle about the first centre, from the direction out of the spine's centre.
    const gp_Vec out = gp_Vec(spine.Location(), first.centre).Normalized();
    const auto tube = [&](const gp_Pnt& foot) {
        const gp_Vec v(first.centre, foot);
        return std::atan2(v.Dot(gp_Vec(axis)), v.Dot(out));
    };
    double v0 = tube(first.foot[0]), v1 = tube(first.foot[1]);
    if (full) {
        // The way round the tube that passes the foot on the face between.
        const auto ahead = [&](double v) { return std::fmod(std::fmod(v - v0, 2 * M_PI) + 2 * M_PI, 2 * M_PI); };
        v1 = ahead(tube(first.middle)) < ahead(v1) ? v0 + ahead(v1) : v0 + ahead(v1) - 2 * M_PI;
    } else if (std::abs(v1 - v0) > M_PI) (v1 < v0 ? v1 : v0) += 2 * M_PI;
    // Its range in the tube's angle from zero up, as OpenCascade keeps it.
    if (std::min(v0, v1) < 0) {
        v0 += 2 * M_PI;
        v1 += 2 * M_PI;
    }
    Sheet sheet{new Geom_ToroidalSurface(ax, spine.Radius(), r), false, start, turn, {v0, v1}};
    if (exact(sheet)) return sheet;
    return {};
}

Handle(Geom2d_BSplineCurve) alongRail(const Handle(Geom2d_BSplineCurve)& curve, const Sheet& sheet) {
    const auto copy = Handle(Geom2d_BSplineCurve)::DownCast(curve->Copy());
    NCollection_Array1<double> knots(1, copy->NbKnots());
    for (int k = 1; k <= copy->NbKnots(); ++k) knots(k) = sheet.rail(copy->Knot(k));
    copy->SetKnots(knots);
    return copy;
}

template <class Curve, class Point>
static Handle(Curve) joined(std::vector<Handle(Curve)> parts, double tolerance) {
    int degree = 0;
    for (const auto& part : parts) degree = std::max(degree, part->Degree());
    std::vector<Point> poles;
    std::vector<double> weights, knots;
    std::vector<int> mults;
    for (size_t m = 0; m < parts.size(); ++m) {
        auto part = Handle(Curve)::DownCast(parts[m]->Copy());
        part->IncreaseDegree(degree);
        // Weights scaled so that the shared pole keeps one weight: the curve is the same.
        const double scale = m > 0 ? weights.back() / part->Weight(1) : 1;
        if (m > 0 && (poles.back().Distance(part->Pole(1)) > tolerance || std::abs(knots.back() - part->Knot(1)) > tol::parameter))
            return {};
        for (int i = m > 0 ? 2 : 1; i <= part->NbPoles(); ++i) {
            poles.push_back(part->Pole(i));
            weights.push_back(part->Weight(i) * scale);
        }
        for (int i = m > 0 ? 2 : 1; i <= part->NbKnots(); ++i) {
            knots.push_back(part->Knot(i));
            mults.push_back(part->Multiplicity(i));
        }
        if (m > 0) mults[mults.size() - part->NbKnots()] = degree;
    }
    NCollection_Array1<Point> p(1, int(poles.size()));
    NCollection_Array1<double> w(1, int(weights.size())), k(1, int(knots.size()));
    NCollection_Array1<int> u(1, int(mults.size()));
    for (size_t i = 0; i < poles.size(); ++i) {
        p(int(i) + 1) = poles[i];
        w(int(i) + 1) = weights[i];
    }
    for (size_t i = 0; i < knots.size(); ++i) {
        k(int(i) + 1) = knots[i];
        u(int(i) + 1) = mults[i];
    }
    return new Curve(p, w, k, u, degree);
}

Handle(Geom2d_BSplineCurve) endToEnd(std::vector<Handle(Geom2d_BSplineCurve)> parts, double tolerance) {
    return joined<Geom2d_BSplineCurve, gp_Pnt2d>(std::move(parts), tolerance);
}

Handle(Geom_BSplineCurve) endToEnd(std::vector<Handle(Geom_BSplineCurve)> parts, double tolerance) {
    return joined<Geom_BSplineCurve, gp_Pnt>(std::move(parts), tolerance);
}

std::optional<Sheet> extrudedSheet(const Rail& rail, Stations& stations) {
    if (rail.curve.GetType() != GeomAbs_Line) return {};
    for (const auto& side : rail.s)
        if (side.cliffCurve || GeomAdaptor_Surface(side.surface).GetType() != GeomAbs_Plane) return {};
    if (rail.law == Law::Full && GeomAdaptor_Surface(rail.middle.surface).GetType() != GeomAbs_Plane) return {};
    const double middle = (rail.first + rail.last) / 2;
    Station s;
    Section section;
    if (!stations.at(middle, s) || !sectionAt(rail, s, section)) return {};
    gp_Pnt p;
    gp_Vec d;
    rail.curve.D1(middle, p, d);
    const double speed = d.Magnitude();
    if (!(speed > limit::degenerate)) return {};
    const double tolerance = rail.resolved(limit::exactSheet);
    // The feet on their contacts at nine places along the rail.
    const auto exact = [&](const Sheet& sheet) {
        for (int k = 0; k <= 8; ++k) {
            Station at;
            if (!stations.at(rail.t(k / 8.0), at)) return false;
            for (int side = 0; side < 2; ++side) {
                const gp_Pnt2d uv = sheet.uv(sheet.rail(k / 8.0), sheet.across[side]);
                if (sheet.surface->Value(uv.X(), uv.Y()).Distance(at.foot[side]) > tolerance) return false;
            }
        }
        return true;
    };
    // A chamfer's section is the line between the feet: the blend is the plane through it and
    // the edge's direction, x along the rail from the middle section and y across from foot to foot.
    if (rail.law == Law::Chamfer) {
        const gp_Vec across(s.foot[0], s.foot[1]);
        if (!(across.Magnitude() > limit::degenerate)) return {};
        const gp_Dir x = readBack(gp_Dir(d));
        const gp_Vec y = across - gp_Vec(x) * across.Dot(gp_Vec(x));
        if (!(y.Magnitude() > limit::degenerate)) return {};
        const gp_Ax3 frame(s.foot[0], gp_Dir(gp_Vec(x).Crossed(y)), x);
        Sheet sheet{new Geom_Plane(frame), false, (rail.from - middle) * speed, (rail.to - rail.from) * speed, {0, y.Magnitude()}};
        if (exact(sheet)) return sheet;
        return {};
    }
    const auto curve = sectionCurve(section);
    // Along the rail from the middle section, so that the surface's v is the rail's length there.
    Sheet sheet{new Geom_SurfaceOfLinearExtrusion(curve, readBack(gp_Dir(d))), true, (rail.from - middle) * speed, (rail.to - rail.from) * speed,
                {curve->FirstParameter(), curve->LastParameter()}};
    if (exact(sheet)) return sheet;
    return {};
}

/// A curve of a support's parameters at the fit's rate u, where its own parameter is the sheet's
/// rail parameter alpha + beta u. A circle's or a conic's rational spline does not run at its
/// angle's rate, so the curve is interpolated at the rate instead: quintic, its spans halved until
/// it is within `tolerance` of the curve between its sites (empty past 1024 spans).
static Handle(Geom2d_BSplineCurve) atRate(const Handle(Geom2d_Curve)& curve, double alpha, double beta, double tolerance) {
    return interpolated([&](double u) { return curve->Value(alpha + beta * u); }, 0, 1, tolerance);
}

Handle(Geom2d_BSplineCurve) interpolated(const std::function<gp_Pnt2d(double)>& at, double first, double last, double tolerance) {
    constexpr int degree = 5;
    for (int spans = 8; spans <= 1024; spans *= 2) {
        const int count = spans + degree;
        NCollection_Array1<double> knots(1, spans + 1);
        NCollection_Array1<int> multiplicities(1, spans + 1);
        for (int i = 0; i <= spans; ++i) {
            knots(i + 1) = first + (last - first) * (double(i) / spans);
            multiplicities(i + 1) = i == 0 || i == spans ? degree + 1 : 1;
        }
        NCollection_Array1<double> flat(1, count + degree + 1);
        BSplCLib::KnotSequence(knots, multiplicities, flat);
        NCollection_Array1<double> sites(1, count);
        BSplCLib::BuildSchoenbergPoints(degree, flat, sites);
        NCollection_Array1<int> orders(1, count);
        orders.Init(0);
        NCollection_Array1<gp_Pnt2d> poles(1, count);
        for (int i = 1; i <= count; ++i) poles(i) = at(sites(i));
        int singular = 0;
        BSplCLib::Interpolate(degree, flat, sites, orders, poles, singular);
        if (singular != 0) return {};
        const Handle(Geom2d_BSplineCurve) spline = new Geom2d_BSplineCurve(poles, knots, multiplicities, degree);
        double worst = 0;
        for (int i = 0; i < spans; ++i)
            for (const double f : {0.25, 0.5, 0.75}) {
                const double u = first + (last - first) * ((i + f) / spans);
                worst = std::max(worst, spline->Value(u).Distance(at(u)));
            }
        if (worst <= tolerance) return spline;
    }
    return {};
}

std::optional<Fitted> exactContacts(const Rail& rail, Stations& stations, const Sheet& sheet) {
    Fitted out;
    const double s0 = sheet.rail(0), s1 = sheet.rail(1);
    for (int side = 0; side < 2; ++side) {
        // On a cliff the contact is in the edge's parameter and the angle about it (Side): through
        // the stations' own, which run along the edge at a constant angle on an exact sheet.
        if (rail.s[side].cliffCurve) {
            constexpr int count = 9;
            NCollection_Array1<gp_Pnt2d> values(1, count);
            NCollection_Array1<double> at(1, count);
            const auto& edge = rail.s[side].cliffCurve;
            for (int i = 0; i < count; ++i) {
                Station s;
                if (!stations.at(rail.t(i / double(count - 1)), s)) return {};
                double w = s.uv[side][0], angle = s.uv[side][1];
                // Round a closed edge, on from the last place in its own period; the angle likewise.
                if (i > 0) {
                    const gp_Pnt2d last = values(i);
                    if (edge->IsPeriodic()) w = loft::nearPeriod(w, last.X(), edge->Period());
                    angle = loft::nearPeriod(angle, last.Y(), 2 * M_PI);
                }
                values(i + 1) = gp_Pnt2d(w, angle);
                at(i + 1) = i / double(count - 1);
            }
            Geom2dAPI_Interpolate interpolate(new NCollection_HArray1<gp_Pnt2d>(values), new NCollection_HArray1<double>(at), false, tol::parameter);
            interpolate.Perform();
            if (!interpolate.IsDone()) return {};
            const Handle(Geom2d_BSplineCurve) spline = interpolate.Curve();
            for (int i = 0; i + 1 < count; ++i) {
                const double u = (i + 0.5) / (count - 1);
                Station s;
                Jet j;
                const gp_Pnt2d laid = spline->Value(u);
                if (!stations.at(rail.t(u), s) || !jetAt(rail.s[side], laid.X(), laid.Y(), j) || j.p.Distance(s.foot[side]) > rail.resolved(limit::exactSheet))
                    return {};
            }
            out.contact[side] = spline;
            continue;
        }
        const auto onSupport = GeomProjLib::Curve2d(sheet.contact(side), std::min(s0, s1), std::max(s0, s1), rail.s[side].surface);
        if (onSupport.IsNull()) return {};
        if (sheet.beta < 0) return {};
        // From the sheet's rail parameter to the fit's u: a line's spline runs at its own rate, so
        // its knots are rescaled; a circle's (or a conic's) does not, so it is laid at the rate.
        Handle(Geom2d_BSplineCurve) spline;
        if (basisOf(onSupport)->IsKind(STANDARD_TYPE(Geom2d_Line))) {
            spline = Geom2dConvert::CurveToBSplineCurve(new Geom2d_TrimmedCurve(onSupport, std::min(s0, s1), std::max(s0, s1)));
            if (spline.IsNull()) return {};
            NCollection_Array1<double> knots(1, spline->NbKnots());
            for (int i = 1; i <= spline->NbKnots(); ++i) knots(i) = (spline->Knot(i) - sheet.alpha) / sheet.beta;
            spline->SetKnots(knots);
        } else {
            spline = atRate(onSupport, sheet.alpha, sheet.beta, rail.resolved(limit::exactSheet) / 4);
            if (spline.IsNull()) return {};
        }
        // The ball's feet are on it, in the support's own period.
        for (int i = 0; i <= 8; ++i) {
            const double u = i / 8.0;
            Station s;
            if (!stations.at(rail.t(u), s)) return {};
            const gp_Pnt2d foot(s.uv[side][0], s.uv[side][1]);
            gp_Pnt2d laid = spline->Value(u);
            if (i == 0) {
                const gp_Vec2d shift(laid, nearPeriod(laid, foot, rail.s[side].surface));
                if (shift.SquareMagnitude() > 0) spline->Translate(shift);
                laid = spline->Value(u);
            }
            if (rail.s[side].surface->Value(laid.X(), laid.Y()).Distance(s.foot[side]) > rail.resolved(limit::exactSheet)) return {};
        }
        out.contact[side] = spline;
    }
    // A contact that runs back against the edge folds the blend over itself (the ball larger than
    // the face is round there), which the fitted contacts' checks find and treat: an exact one is
    // taken only where both run on with the edge.
    for (int side = 0; side < 2; ++side) {
        if (rail.s[side].cliffCurve || !rail.s[side].trim.IsNull()) continue;
        for (const double u : {0.25, 0.5, 0.75}) {
            gp_Pnt2d uv;
            gp_Vec2d duv;
            out.contact[side]->D1(u, uv, duv);
            gp_Pnt foot, at;
            gp_Vec su, sv, along;
            rail.s[side].surface->D1(uv.X(), uv.Y(), foot, su, sv);
            rail.curve.D1(rail.t(u), at, along);
            const gp_Vec run = su * duv.X() + sv * duv.Y();
            if (!(along.Magnitude() > limit::degenerate) || !(run.Dot(along.Normalized()) > 1e-3 * run.Magnitude())) return {};
        }
    }
    return out;
}

} // namespace loft::blend2::detail
