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

Standard_EXPORT bool sameBlendSurface(const TopoDS_Face& a, const TopoDS_Face& b);
Standard_EXPORT bool tangentBlendEdge(const TopoDS_Shape& shape, const TopoDS_Edge& edge);
Standard_EXPORT BlendBoundaryRoles
                     blendBoundaryRoles(const TopoDS_Shape&              before,
                                        const std::vector<TopoDS_Shape>& selected,
                                        const TopoDS_Shape&              result,
                                        const BooleanResult&             history,
                                        const std::vector<TopoDS_Shape>& generated = {});
Standard_EXPORT void certifyBlendBoundaries(const char*               op,
                                            const TopoDS_Shape&       result,
                                            const BlendBoundaryRoles& roles,
                                            double                    scale,
                                            double                    angleDegrees = .5,
                                            int                       samples      = 65);
} // namespace BRepExt
