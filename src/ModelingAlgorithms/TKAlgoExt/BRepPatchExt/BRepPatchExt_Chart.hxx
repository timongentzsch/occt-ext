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
// Alternatively, this file may be used under the terms of Open CASCADE
// commercial license or contractual agreement.

#pragma once

// A chart of a loop of its own, for a fill that is a height over no plane: the loop spanned by a
// disc of triangles that leaves it as the faces along it ask, laid flat conformally. And what a
// fill over such a chart is held to that a height is not: that no plane is one it could be a
// height over, which is when it is built so, and that it is regular where it is used.

#include <Standard_Macro.hxx>
#include <BRepAlgoExt_PatchFlow.hxx>

#include <Geom_BSplineSurface.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <functional>
#include <optional>
#include <vector>

namespace BRepExt
{

//! A place on a loop and the way from it to the next. `leave` is the way the fill is to leave
//! that piece of the loop, a unit vector in space, where a face it is to be tangent to says so;
//! zero where nothing says.
struct ChartRim
{
  gp_Pnt point;
  gp_Vec leave;
};

//! The loop `rim`, closed, spanned by a disc of triangles and laid flat (see layFlat). The disc
//! is an even lattice on the unit circle, the loop on its rim at even steps of its length, and
//! its places in space are the fairest that leave every piece of the loop along its `leave`
//! (see BRepPatchExt_Chart.cxx). The chart's first nodes are the loop's places, in turn; its
//! `source` is the chart it was last laid out in. It is drawn to pass through the points
//! `through` as well. Fails with no_solution where the flattening folds a triangle over.
Standard_EXPORT ConformalChart proxyChart(const std::vector<ChartRim>& rim,
                                          const std::vector<gp_Pnt>&   through);

//! Where in the chart the point of the chart's triangles nearest to `point` lies, and how far
//! that is from it in space.
Standard_EXPORT gp_XY onChart(const ConformalChart& chart, const gp_Pnt& point, double& distance);

//! How near some plane comes to having every one of `normals` within its own normal's cone: the
//! largest, over the directions d, of the least |n . d|. A surface with these normals along its
//! rim is a height over the plane across d with no slope above tan(acos(that)) there, and over
//! none with less. `axes` are three orthonormal directions the search is laid out in, so that
//! the answer turns with them.
Standard_EXPORT double graphMargin(const std::vector<gp_Vec>& normals, const gp_Vec (&axes)[3]);

//! Whether S_u x S_v keeps away from zero on every cell of `surface`'s knot spans that `used`
//! names by its parameter box (u0, u1, v0, v1): per cell the polynomial S_u x S_v in Bernstein
//! form, whose coefficients all lie on one side of a plane through the origin where it does,
//! split in four where they do not, a few times. Returns the least cosine, over the cells, of
//! the coefficients against the cell's mean direction; nothing where a cell is not shown regular.
Standard_EXPORT std::optional<double> regularMargin(
  const occ::handle<Geom_BSplineSurface>&                    surface,
  const std::function<bool(double, double, double, double)>& used);

} // namespace BRepExt
