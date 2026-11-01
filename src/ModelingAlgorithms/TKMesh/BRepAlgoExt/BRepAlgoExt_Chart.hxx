// Copyright (c) 2026 Timon Gentzsch
//
// This file is part of Open CASCADE Technology software library.
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult LICENSE_LGPL_21.txt for the complete license.

#pragma once
#include <Standard.hxx>
#include <GeomAlgoExt_Execution.hxx>
// A piecewise-linear UV-flow chart approximating a conformal flattening.
// Harmonic boundary interpolation does not guarantee exact conformality.
#include <TopoDS_Face.hxx>
#include <gp_Pnt.hxx>
#include <gp_XY.hxx>

#include <array>
#include <vector>

namespace BRepAlgoExt {

/// A Boundary First Flattening chart (Sawhney and Crane 2017), with zero boundary
/// log scale, a least-change boundary closure and harmonic coordinate extension.
/// Closure can alter boundary lengths; exact discrete conformality is not promised.
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
    Standard_EXPORT gp_XY map(const gp_XY &uv) const;

    // Triangles by cell of an even grid over the face parameters, for `map`.
    gp_XY low, cell;
    int columns = 0, rows = 0;
    std::vector<std::vector<int>> cells;
};

/// Build a sampled triangulation/chart driven by spacing and deflection targets.
/// Finite curvature sampling and the node cap do not certify global 3D deflection.
/// Fails with no_solution on an unsupported disc or an inverted output triangle.
Standard_EXPORT ConformalChart makeConformalChart(const TopoDS_Face &face, double spacing, double deflection, double parameterResolution, GeomAlgoExt::StageObserver observer);

} // namespace BRepAlgoExt
