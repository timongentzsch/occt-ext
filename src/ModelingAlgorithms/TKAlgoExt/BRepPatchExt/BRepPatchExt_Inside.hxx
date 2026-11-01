// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
// What a patch is like inside its loop. The misses along its edges say nothing of that: a
// surface can meet every edge and still stand up as a spike between them, or fold over.
#include <TopoDS_Face.hxx>
#include <optional>
#include <string>

namespace loft {
/// What is wrong with the inside of a face, if anything is: that it folds over, or that it
/// bends many times more tightly than anything on its rim does (`rimBend`, a curvature; zero
/// or less leaves that unjudged).
std::optional<std::string> insideFault(const TopoDS_Face &face, double rimBend);

struct InsideCurvature {
    double low = 0, high = 0;
    double wrongSign() const;
    double tightest() const;
};
InsideCurvature insideCurvature(const TopoDS_Face &face);
/// Largest normal-curvature step across a knot below C2, away from the rim.
double insideStep(const TopoDS_Face &face);

/// The largest principal curvature of a face at the places it is looked at inside its loop.
double tightestInside(const TopoDS_Face &face);
} // namespace loft
