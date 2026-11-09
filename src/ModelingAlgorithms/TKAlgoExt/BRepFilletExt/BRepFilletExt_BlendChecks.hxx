// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once

#include <Standard_Macro.hxx>
#include <BRepAlgoExt_Cells.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>

namespace BRepExt
{

// Roles come from construction/history, never from the angle of the delivered seam.
struct BlendBoundaryRoles
{
  ShapeSet                                      supports, stops, blends, sharp;
  double                                        constructionTolerance = 0;
  std::vector<std::pair<TopoDS_Face, ShapeSet>> contacts;
};

//! Whether faces `a` and `b` lie on one surface: the same support, or equal planes or cylinders.
Standard_EXPORT bool sameBlendSurface(const TopoDS_Face& a, const TopoDS_Face& b);
//! Whether the faces either side of `edge` meet tangentially (a seam or degenerate edge does).
//! Fails where the edge does not lie between two faces.
Standard_EXPORT bool tangentBlendEdge(const TopoDS_Shape& shape, const TopoDS_Edge& edge);
//! The role of each boundary of the blend faces of `result`, read from how they were made.
Standard_EXPORT BlendBoundaryRoles
                     blendBoundaryRoles(const TopoDS_Shape&              before,
                                        const std::vector<TopoDS_Shape>& selected,
                                        const TopoDS_Shape&              result,
                                        const BooleanResult&             history,
                                        const std::vector<TopoDS_Shape>& generated = {});
//! Fails naming `op` and the edge where a boundary of `result` leaves its curve or, at any of
//! `samples` places, its faces' normals part by more than `angleDegrees` where its role is G1.
Standard_EXPORT void certifyBlendBoundaries(const char*               op,
                                            const TopoDS_Shape&       result,
                                            const BlendBoundaryRoles& roles,
                                            double                    scale,
                                            double                    angleDegrees = .5,
                                            int                       samples      = 65);
} // namespace BRepExt
