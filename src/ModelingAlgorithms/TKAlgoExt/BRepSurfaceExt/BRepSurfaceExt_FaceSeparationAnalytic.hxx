// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once

// The analytic branch of the self-interference pre-pass (face_separation.h): faces of planes
// beside faces of surfaces of revolution about an axis normal to them, proven to meet only along
// the edges they share (face_separation_analytic.cpp).

#include <NCollection_IndexedMap.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>

namespace loft {

/// Whether two faces that have `shared` in common (edges and vertices) meet only there and
/// come within their tolerances of each other only beside it. False where this cannot tell.
bool analyticBeside(const TopoDS_Face& a, const TopoDS_Face& b,
                    const NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>& shared);

/// Whether `shape` has two faces beside each other that analyticBeside is for.
bool analyticWorth(const TopoDS_Shape& shape);

} // namespace loft
