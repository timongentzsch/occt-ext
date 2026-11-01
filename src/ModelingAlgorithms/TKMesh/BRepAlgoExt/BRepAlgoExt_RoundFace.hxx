// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once

// A face that is a constant round: a cylinder or torus tangent to two supports along its edges,
// as Measure's fillet dimension, Refillet and fillet v2's re-blending recognise it.

#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Torus.hxx>

#include <vector>

namespace BRepExt {

struct RoundFace {
    bool valid = false;
    gp_Cylinder cylinder;
    gp_Torus torus;
    bool toroidal = false;
    std::vector<TopoDS_Face> supports;
    /// The round's radius: the cylinder's, or the torus's minor one.
    double radius() const { return toroidal ? torus.MinorRadius() : cylinder.Radius(); }
};

/// `face` of `shape` as a round: valid where it is a cylinder (not closed round) or a torus that
/// is tangent, to `tolerance`, along every edge it shares with exactly two neighbours, its supports.
RoundFace recognizeRound(const TopoDS_Shape& shape, const TopoDS_Face& face, double tolerance);

} // namespace BRepExt
