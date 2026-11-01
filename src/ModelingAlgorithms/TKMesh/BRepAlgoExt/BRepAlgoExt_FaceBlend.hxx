// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once

// The round or bevel of the whole outline of a flat face against the walls that stand square
// on it, past the tightest bend of that outline (blend_face.cpp).

#include <BRepAlgoExt_Cells.hxx>

#include <vector>

namespace loft {

/// `shape` with every edge of one flat face blended by `distance`: a round (the ball's arc) or,
/// with `bevel`, a chamfer of that setback on both sides. `edges` must be all the edges of one
/// flat face, whose walls stand square on it for at least `distance`.
///
/// The surface is the one each place of the outline makes by itself, cut where two of them
/// meet: over a point at the distance d from the outline it lies sqrt(r^2 - (r - d)^2) - r
/// (a bevel: d - r) from the face's plane, for d < r. Where the outline bends more tightly
/// than the blend those pieces meet in creases over the outline's medial axis.
/// A closing concave round instead uses the ball envelope on the trimmed inward offset.
/// Currently built for one convex polygon of straight centre paths and planar/cylindrical
/// vertical prism walls ending at one planar mouth; other closing concave geometry refuses.
///
/// Fails with "the planar cap: top face vanishes" in its text where nothing of the face is
/// left, and otherwise with what it cannot build.
BooleanResult faceBlend(const char* operation, const TopoDS_Shape& shape, const std::vector<TopoDS_Shape>& edges,
                        double distance, bool bevel);

} // namespace loft
