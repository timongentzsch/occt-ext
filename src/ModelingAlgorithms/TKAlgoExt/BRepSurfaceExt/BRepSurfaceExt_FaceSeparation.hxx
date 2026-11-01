// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once

#include <TopoDS_Shape.hxx>
#include <array>
#include <vector>

namespace BRepExt {

// What a self-interference check of a shape (BOPAlgo_ArgumentAnalyzer's SelfInterMode, OCCT's
// CheckerSI) can leave out: a conservative pre-pass over the shape's faces, which proves pairs
// of faces clear of each other and faces clear of themselves, and leaves everything else to the
// check. The proofs bound the exact surfaces through the Bezier poles of their supports (pole
// hulls, and Bernstein coefficients of the derivatives), never samples of them:
// - apart: two faces' pole hulls, halved as needed, further apart than their tolerances;
// - clear of itself: a spline face's patches that are not neighbours apart, and every 2 x 2
//   block of them one to one, by a bound on the Jacobian over all coefficients of its
//   derivatives (across its seams too, for a closed face);
// - meeting only along a shared edge: two spline faces put side by side in one parameter
//   rectangle across the edge, the blocks across it one to one by the same bound, the rest of
//   their patches apart, and the two within their tolerances of each other only near the edge
//   (NEAR_EDGE times them: a crease of about 15 degrees or more); a spline face beside a planar
//   face it shares edges with leaves the plane monotonically from them and is apart from the
//   planar face elsewhere.
// Planar faces are not proven clear of themselves, nor anything about a face the pre-pass does
// not read as a rectangle of Bezier patches with straight parameter-line edges. What is not
// proven within its budget (near-touching, tangent or shallow creases, singular or trimmed
// faces, free edges) goes to the check, which then sees it as it saw it before.
struct SeparationCount {
    int pairs = 0, proven = 0;
    // The pairs proven, as indices of the shape's faces in TopExp::MapShapes order (i, i for a
    // face clear of itself).
    std::vector<std::array<int, 2>> provenPairs;
};
// What the check still has to see: a part for each set of faces joined by pairs not proven,
// a face alone where only it is not proven clear of itself; `shape` itself where that is all
// of its faces; nothing where all is proven. A face of few patches is left to the check to
// meet with itself, which costs it little.
std::vector<TopoDS_Shape> unprovenParts(const TopoDS_Shape& shape, SeparationCount* count = nullptr);

} // namespace BRepExt
