// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
// A chart for UV flow that keeps right angles: a face laid flat conformally. Over such a chart
// a surface's U and V meet at right angles wherever the surface follows the face.
#include <TopoDS_Face.hxx>
#include <gp_Pnt.hxx>
#include <gp_XY.hxx>

#include <array>
#include <vector>

namespace BRepExt {

/// A face's triangulation laid flat by Boundary First Flattening with the boundary's lengths
/// kept (Sawhney and Crane 2017, scale zero along the boundary): angles are kept everywhere,
/// lengths along the rim, and inside the area changes only as the face's curvature makes it.
/// The chart is turned to the principal axes of its own area, centred on it.
struct ConformalChart {
    /// Per node: the face's parameters, the place in space and the place in the chart.
    std::vector<gp_XY> source, flat;
    std::vector<gp_Pnt> point;
    std::vector<std::array<int, 3>> triangles;
    /// The boundary nodes in turn, anticlockwise in the chart.
    std::vector<int> boundary;

    /// Where face parameters lie in the chart: by the triangle that holds them, or linearly
    /// on from the nearest one for parameters just outside (a curve between two boundary
    /// nodes bulges out of their chord).
    gp_XY map(const gp_XY &uv) const;

    // Triangles by cell of an even grid over the face parameters, for `map`.
    gp_XY low, cell;
    int columns = 0, rows = 0;
    std::vector<std::vector<int>> cells;
};

/// Triangulate a face with edges no longer than `spacing` and within `deflection` of the face
/// in space, and lay it flat. Fails with no_solution where the face is not one disc or the
/// flattening folds a triangle over.
ConformalChart conformalChart(const TopoDS_Face &face, double spacing, double deflection);

} // namespace BRepExt
