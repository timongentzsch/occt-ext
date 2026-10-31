// Copyright (c) 2026 Timon Gentzsch
//
// This file is part of Open CASCADE Technology software library.
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult LICENSE_LGPL_21.txt for the complete license.


#include <GeomAlgoExt_RimFit.hxx>
#include <GeomAlgoExt_Band.hxx>
#include <GeomAlgoExt_Execution.hxx>

#include <NCollection_Array1.hxx>

#include <algorithm>
#include <cmath>
#include <functional>
#include <cstdio>

namespace GeomAlgoExt {
namespace {

/// A span shorter than this share of the mean one is beside a knot of lower continuity, and
/// says nothing of how finely the net follows a curve.
constexpr double sliver = 0.2;

/// The solves on one net, the places looked for again on the surface between them: the first
/// takes all that a linear step can, the next what moving the places left.
constexpr int solves = 3;

/// What a pole's own movement counts for against a place it moves the surface at, whose
/// basis functions sum to one: a pole with less than this share in every place along the curve
/// is moved little, and where the curve is still missed the span there is halved instead.
/// Without it a pole in a corner of a span that the curve only clips was moved far to bring
/// the curve a little nearer, and z = sin 2a on a radius of 6 rippled to a bend of 0.43.
constexpr double damping = 0.03;

/// A solve that leaves more than this share of the miss before it has done what the net can.
constexpr double stillFalling = 0.5;

/// Where the surface stands more steeply than this against the direction its poles move
/// in, as the cosine, a place is left as it is: 80 degrees, as for the net itself.
constexpr double steepest = 0.17;

/// One thing asked of the movement, over the poles that it touches.
struct Row {
    std::vector<int> index;
    std::vector<double> coefficient;
    double value;
};

int spans(const std::vector<double> &knots) {
    int count = 0;
    for (size_t i = 0; i + 1 < knots.size(); ++i)
        count += knots[i + 1] > knots[i];
    return count;
}

/// The knots on both sides of a parameter.
std::pair<double, double> spanAt(const std::vector<double> &knots, int degree, double t) {
    const auto last = knots.end() - degree - 1;
    auto above = std::upper_bound(knots.begin() + degree, last, t);
    if (above == last && *(above - 1) >= *above)
        --above;
    return {*(above - 1), *above};
}

} // namespace

gp_XY nearest(const occ::handle<Geom_Surface> &surface, gp_XY uv, const gp_Pnt &point) {
    double u0, u1, v0, v1;
    surface->Bounds(u0, u1, v0, v1);
    for (int step = 0; step < 8; ++step) {
        gp_Pnt at;
        gp_Vec su, sv;
        surface->D1(uv.X(), uv.Y(), at, su, sv);
        const gp_Vec to(at, point);
        const double E = su.Dot(su), F = su.Dot(sv), G = sv.Dot(sv), det = E * G - F * F;
        if (det <= 1e-12 * E * G)
            break;
        const gp_XY move((G * to.Dot(su) - F * to.Dot(sv)) / det, (E * to.Dot(sv) - F * to.Dot(su)) / det);
        uv.SetCoord(surface->IsUPeriodic() ? uv.X() + move.X() : std::clamp(uv.X() + move.X(), u0, u1),
                    surface->IsVPeriodic() ? uv.Y() + move.Y() : std::clamp(uv.Y() + move.Y(), v0, v1));
        if (move.X() * move.X() * E + move.Y() * move.Y() * G < 1e-28 * (E + G))
            break;
    }
    return uv;
}

SplineNet netOf(const occ::handle<Geom_BSplineSurface> &surface) {
    SplineNet net{surface->UDegree(), surface->VDegree(), {}, {}, {}};
    const NCollection_Array1<double> &us = surface->UKnotSequence(), &vs = surface->VKnotSequence();
    net.knotsU.assign(us.begin(), us.end());
    net.knotsV.assign(vs.begin(), vs.end());
    const int nu = surface->NbUPoles(), nv = surface->NbVPoles();
    net.poles.resize(static_cast<size_t>(nu) * nv);
    for (int i = 0; i < nu; ++i)
        for (int j = 0; j < nv; ++j)
            net.poles[static_cast<size_t>(i) * nv + j] = surface->Pole(i + 1, j + 1).XYZ();
    return net;
}

double refineRim(const char *op, SplineNet &net, const std::function<std::vector<RimPlace>(double)> &placesFor,
                  const gp_XYZ &up, double within, double enough, int mostSpans, const Message_ProgressRange& range, StageObserver observer) {
    FitPhase phase(observer, "refinement");
    occ::handle<Geom_BSplineSurface> surface = surfaceOf(net);
    // The shortest span of the net either way, slivers beside a knot of lower continuity apart.
    const auto shortest = [&] {
        double least = HUGE_VAL;
        for (const auto *knots : {&net.knotsU, &net.knotsV}) {
            const double mean = (knots->back() - knots->front()) / spans(*knots);
            for (size_t i = 0; i + 1 < knots->size(); ++i)
                if ((*knots)[i + 1] - (*knots)[i] > sliver * mean)
                    least = std::min(least, (*knots)[i + 1] - (*knots)[i]);
        }
        return least;
    };
    std::vector<RimPlace> places = placesFor(shortest());
    // How far the surface is to rise along `up` at each place, for its nearest point to come
    // to the place: the miss is along the surface's normal there.
    std::vector<double> rise, miss;
    const auto look = [&] {
        double worst = 0;
        rise.resize(places.size());
        miss.resize(places.size());
        for (size_t k = 0; k < places.size(); ++k) {
            RimPlace &place = places[k];
            place.uv = nearest(surface, place.uv, place.point);
            gp_Pnt at;
            gp_Vec su, sv;
            surface->D1(place.uv.X(), place.uv.Y(), at, su, sv);
            const gp_Vec to(at, place.point);
            miss[k] = to.Magnitude();
            gp_Vec normal = su.Crossed(sv);
            const double length = normal.Magnitude(), lean = length > 0 ? normal.XYZ().Dot(up) / length : 0;
            rise[k] = std::abs(lean) > steepest ? to.Dot(normal) / length / lean : 0;
            worst = std::max(worst, miss[k]);
        }
        return worst;
    };
    double worst = look();
    if (worst <= enough)
        return worst;
    while (worst > within) {
        checkCancelled(range);
        const int nv = net.countV(), pu = net.degreeU, pv = net.degreeV;
        const int size = net.countU() * nv;
        for (int solve = 0; solve < solves && worst > within; ++solve) {
            // Per place how far it is to rise, and where the tangent plane is to stay, that
            // the rise does not change across the curve.
            std::vector<Row> rows;
            rows.reserve(places.size() * 2);
            for (size_t k = 0; k < places.size(); ++k) {
                const RimPlace &place = places[k];
                const SplineBasis bu = splineBasisAt(op, net.knotsU, pu, place.uv.X()),
                                  bv = splineBasisAt(op, net.knotsV, pv, place.uv.Y());
                const bool held = place.across.SquareModulus() > 0;
                const auto [u0, u1] = spanAt(net.knotsU, pu, place.uv.X());
                const auto [v0, v1] = spanAt(net.knotsV, pv, place.uv.Y());
                const double du = place.across.X() * (u1 - u0), dv = place.across.Y() * (v1 - v0);
                const bool bent = held && place.bend;
                Row at{{}, {}, rise[k]}, tilt{{}, {}, 0}, bend{{}, {}, 0};
                for (int a = 0; a <= pu; ++a)
                    for (int b = 0; b <= pv; ++b) {
                        const int pole = (bu.first + a) * nv + bv.first + b;
                        at.index.push_back(pole);
                        at.coefficient.push_back(bu.value[0][a] * bv.value[0][b]);
                        if (held) {
                            tilt.index.push_back(pole);
                            tilt.coefficient.push_back(du * bu.value[1][a] * bv.value[0][b] +
                                                       dv * bu.value[0][a] * bv.value[1][b]);
                        }
                        if (bent) {
                            bend.index.push_back(pole);
                            bend.coefficient.push_back(du * du * bu.value[2][a] * bv.value[0][b] +
                                                       2 * du * dv * bu.value[1][a] * bv.value[1][b] +
                                                       dv * dv * bu.value[0][a] * bv.value[2][b]);
                        }
                    }
                rows.push_back(std::move(at));
                if (held)
                    rows.push_back(std::move(tilt));
                if (bent)
                    rows.push_back(std::move(bend));
            }
            // The least movement of the poles that meets the rows as nearly as the net can, by
            // least squares. Damped: each pole's own movement counts as a row of its own, so
            // that a pole that hardly reaches the curve stays where it is, and the span is made
            // finer instead. Solved exactly, by the normal equations over the poles the rows
            // touch, a band along the curve in the order of the net: the other poles stay.
            std::vector<int> unknown(size, -1);
            for (const Row &row : rows)
                for (int pole : row.index)
                    unknown[pole] = 0;
            int count = 0;
            for (int &index : unknown)
                if (index == 0) index = count++;
            int half = 0;
            for (const Row &row : rows) {
                const auto [low, high] = std::minmax_element(row.index.begin(), row.index.end());
                half = std::max(half, unknown[*high] - unknown[*low]);
            }
            Band normal(count, half, true);
            std::vector<double> right(count, 0);
            for (int i = 0; i < count; ++i)
                normal.add(i, i, damping * damping);
            for (size_t r = 0; r < rows.size(); ++r) {
                if ((r & 1023) == 0) checkCancelled(range);
                const Row &row = rows[r];
                for (size_t a = 0; a < row.index.size(); ++a) {
                    right[unknown[row.index[a]]] += row.coefficient[a] * row.value;
                    for (size_t b = 0; b <= a; ++b)
                        normal.add(unknown[row.index[a]], unknown[row.index[b]], row.coefficient[a] * row.coefficient[b]);
                }
            }
            checkCancelled(range);
            if (!normal.factor())
                fail(op, no_solution, "the movement onto the curves cannot be solved: a place is not a number");
            normal.solve(right);
            std::vector<double> x(size, 0);
            for (int i = 0; i < size; ++i)
                if (unknown[i] >= 0) x[i] = right[unknown[i]];
            for (size_t p = 0; p < net.poles.size(); ++p)
                net.poles[p] += up * x[p];
            surface = surfaceOf(net);
            const double before = worst;
            worst = look();
            if (worst > stillFalling * before)
                break;
        }
        if (worst <= within)
            break;
        // What is left is what this net cannot follow: halve the spans it is left in.
        std::vector<double> addU, addV;
        const bool roomU = spans(net.knotsU) < mostSpans, roomV = spans(net.knotsV) < mostSpans;
        for (size_t k = 0; k < places.size(); ++k) {
            if (miss[k] <= within)
                continue;
            const auto [u0, u1] = spanAt(net.knotsU, pu, places[k].uv.X());
            const auto [v0, v1] = spanAt(net.knotsV, pv, places[k].uv.Y());
            if (roomU)
                addU.push_back((u0 + u1) / 2);
            if (roomV)
                addV.push_back((v0 + v1) / 2);
        }
        if (addU.empty() && addV.empty())
            break;
        for (auto *add : {&addU, &addV}) {
            std::sort(add->begin(), add->end());
            add->erase(std::unique(add->begin(), add->end()), add->end());
        }
        for (double knot : addU)
            surface->InsertUKnot(knot, 1, 0, true);
        for (double knot : addV)
            surface->InsertVKnot(knot, 1, 0, true);
        net = netOf(surface);
        // As many places to a span as before: fewer than the net has poles along the curve
        // would be met each, with the surface anywhere between them.
        places = placesFor(shortest());
        worst = look();
    }
    return worst;
}

} // namespace GeomAlgoExt
