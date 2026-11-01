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
#pragma once

#include <Geom_BSplineCurve.hxx>
#include <NCollection_Array1.hxx>
#include <gp_Pnt.hxx>
#include <gp_XYZ.hxx>

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

// What the skins of skins are made of.

namespace GeomAlgoExt {

/// The cubic spline through values at rising places, C2 with a knot at each place: the unique such
/// spline once its two end conditions are fixed. `Free` ends have no second derivative (the natural
/// spline, which among all C2 interpolants minimises the integral of the squared second
/// derivative). `Follow` ends take the slope of the cubic through the first and last four values,
/// which makes it good to the fourth power of the spacing for a smooth function. `Closed` makes the
/// first and last value one, and the slope and second derivative meet there too. What it gives for
/// each set of values is its poles.
struct CubicThrough {
    enum class Ends { Free, Follow, Closed };
    int poles;
    NCollection_Array1<double> knots;
    NCollection_Array1<int> mults;
    std::vector<double> h;
    Ends ends;
    explicit CubicThrough(const std::vector<double> &places, Ends endCondition = Ends::Free)
        : poles(int(places.size()) + 2), knots(1, int(places.size())), mults(1, int(places.size())), ends(endCondition) {
        const int n = int(places.size());
        for (int i = 0; i < n; ++i) {
            knots(i + 1) = places[size_t(i)];
            mults(i + 1) = i == 0 || i == n - 1 ? 4 : 1;
            if (i + 1 < n)
                h.push_back(places[size_t(i + 1)] - places[size_t(i)]);
        }
    }
    std::vector<gp_Pnt> through(const std::vector<gp_Pnt> &values) const {
        const size_t n = values.size();
        const auto slope = [&](size_t k) { return (values[k + 1].XYZ() - values[k].XYZ()) / h[k]; };
        // The slopes at the places, by the tridiagonal system of a C2 spline:
        // h[k] m[k-1] + 2 (h[k-1] + h[k]) m[k] + h[k-1] m[k+1] = 3 (h[k] d[k-1] + h[k-1] d[k]).
        std::vector<gp_XYZ> m(n);
        if (ends == Ends::Closed) {
            // The same round the loop, one unknown fewer: Thomas's algorithm with the corner terms
            // taken out by the Sherman-Morrison formula.
            const size_t c = n - 1;
            const auto before = [&](size_t k) { return h[(k + c - 1) % c]; };
            std::vector<double> diagonal(c), upper(c), lower(c);
            std::vector<gp_XYZ> right(c);
            for (size_t k = 0; k < c; ++k) {
                lower[k] = h[k];
                diagonal[k] = 2 * (before(k) + h[k]);
                upper[k] = before(k);
                right[k] = (slope((k + c - 1) % c) * h[k] + slope(k) * before(k)) * 3;
            }
            if (c == 1) {
                m[0] = m[1] = right[0] / (2 * diagonal[0]);
                m[0] = m[1] = gp_XYZ(0, 0, 0);
            } else {
                const double gamma = -diagonal[0];
                std::vector<double> d = diagonal, z(c, 0.0), f(c);
                d[0] -= gamma;
                d[c - 1] -= lower[0] * upper[c - 1] / gamma;
                z[0] = gamma;
                z[c - 1] = upper[c - 1];
                std::vector<gp_XYZ> y = right;
                f[0] = upper[0] / d[0];
                y[0] = y[0] / d[0];
                z[0] /= d[0];
                for (size_t k = 1; k < c; ++k) {
                    const double pivot = d[k] - lower[k] * f[k - 1];
                    f[k] = upper[k] / pivot;
                    y[k] = (y[k] - y[k - 1] * lower[k]) / pivot;
                    z[k] = (z[k] - z[k - 1] * lower[k]) / pivot;
                }
                for (size_t k = c - 1; k-- > 0;) {
                    y[k] -= y[k + 1] * f[k];
                    z[k] -= z[k + 1] * f[k];
                }
                const double share = 1 / (1 + z[0] + z[c - 1] * lower[0] / gamma);
                const gp_XYZ fix = (y[0] + y[c - 1] * (lower[0] / gamma)) * share;
                for (size_t k = 0; k < c; ++k)
                    m[k] = y[k] - fix * z[k];
                m[c] = m[0];
            }
        } else {
            const bool follow = ends == Ends::Follow && n >= 4;
            // The slope at a place of the cubic through four values, in Newton's divided differences.
            const auto endSlope = [&](size_t a, size_t b, size_t c, size_t d) {
                const auto place = [&](size_t i) { return knots(int(i) + 1); };
                const gp_XYZ ab = (values[b].XYZ() - values[a].XYZ()) / (place(b) - place(a));
                const gp_XYZ bc = (values[c].XYZ() - values[b].XYZ()) / (place(c) - place(b));
                const gp_XYZ cd = (values[d].XYZ() - values[c].XYZ()) / (place(d) - place(c));
                const gp_XYZ abc = (bc - ab) / (place(c) - place(a));
                const gp_XYZ bcd = (cd - bc) / (place(d) - place(b));
                const gp_XYZ abcd = (bcd - abc) / (place(d) - place(a));
                return ab + abc * (place(a) - place(b)) + abcd * ((place(a) - place(b)) * (place(a) - place(c)));
            };
            std::vector<double> lower(n), main(n), upper(n);
            std::vector<gp_XYZ> right(n);
            main[0] = 1;
            upper[0] = follow ? 0 : 0.5;
            right[0] = follow ? endSlope(0, 1, 2, 3) : slope(0) * 1.5;
            for (size_t k = 1; k + 1 < n; ++k) {
                lower[k] = h[k];
                main[k] = 2 * (h[k - 1] + h[k]);
                upper[k] = h[k - 1];
                right[k] = (slope(k - 1) * h[k] + slope(k) * h[k - 1]) * 3;
            }
            main[n - 1] = 1;
            lower[n - 1] = follow ? 0 : 0.5;
            right[n - 1] = follow ? endSlope(n - 1, n - 2, n - 3, n - 4) : slope(n - 2) * 1.5;
            for (size_t k = 1; k < n; ++k) {
                const double f = lower[k] / main[k - 1];
                main[k] -= f * upper[k - 1];
                right[k] -= right[k - 1] * f;
            }
            m[n - 1] = right[n - 1] / main[n - 1];
            for (size_t k = n - 1; k-- > 0;)
                m[k] = (right[k] - m[k + 1] * upper[k]) / main[k];
        }
        // Its poles, each from the inner Bezier points of a span and the pole before.
        std::vector<gp_Pnt> out(static_cast<size_t>(poles));
        out[0] = values[0];
        out[1] = gp_Pnt(values[0].XYZ() + m[0] * (h[0] / 3));
        for (size_t k = 0; k + 1 < n; ++k) {
            const gp_XYZ inner = values[k + 1].XYZ() - m[k + 1] * (h[k] / 3);
            const double previous = k > 0 ? h[k - 1] : 0, next = k + 2 < n ? h[k + 1] : 0;
            out[k + 2] = gp_Pnt((inner * (previous + h[k] + next) - out[k + 1].XYZ() * next) / (previous + h[k]));
        }
        out[n + 1] = values[n - 1];
        return out;
    }
};

/// Curves of one degree with the same knots, by raising degrees and inserting knots, which moves
/// none of them.
inline void unify(std::vector<Handle(Geom_BSplineCurve)> &curves) {
    int degree = 1;
    for (const auto &c : curves)
        degree = std::max(degree, c->Degree());
    std::vector<std::pair<double, int>> knots;
    for (const auto &c : curves) {
        c->IncreaseDegree(degree);
        for (int k = 1; k <= c->NbKnots(); ++k) {
            auto same = std::find_if(knots.begin(), knots.end(),
                                     [&](const auto &e) { return std::abs(e.first - c->Knot(k)) <= 1e-12; });
            if (same == knots.end())
                knots.emplace_back(c->Knot(k), c->Multiplicity(k));
            else
                same->second = std::max(same->second, c->Multiplicity(k));
        }
    }
    for (const auto &c : curves)
        for (const auto &[knot, multiplicity] : knots)
            c->InsertKnot(knot, std::min(multiplicity, degree), 1e-12, false);
}

} // namespace GeomAlgoExt
