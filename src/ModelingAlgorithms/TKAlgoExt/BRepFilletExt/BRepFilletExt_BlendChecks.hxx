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

#include <BRepAlgoExt_Cells.hxx>
#include <BRepToolsExt_GeomUtil.hxx>
#include <BRepToolsExt_Checks.hxx>
#include <NCollection_List.hxx>
#include <functional>
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

bool               sameBlendSurface(const TopoDS_Face& a, const TopoDS_Face& b);
bool               tangentBlendEdge(const TopoDS_Shape& shape, const TopoDS_Edge& edge);
BlendBoundaryRoles blendBoundaryRoles(const TopoDS_Shape&              before,
                                      const std::vector<TopoDS_Shape>& selected,
                                      const TopoDS_Shape&              result,
                                      const BooleanResult&             history,
                                      const std::vector<TopoDS_Shape>& generated = {});
void               markMiterBoundaries(const TopoDS_Shape&              before,
                                       const std::vector<TopoDS_Shape>& selected,
                                       const TopoDS_Shape&              result,
                                       BlendBoundaryRoles&              roles);
void               certifyBlendBoundaries(const char*               op,
                                          const TopoDS_Shape&       result,
                                          const BlendBoundaryRoles& roles,
                                          double                    scale,
                                          double                    angleDegrees = .5,
                                          int                       samples      = 65);
//! The faces a face of `before` became in `after`, as the algorithm that made it reports them.
using BlendModified = std::function<const NCollection_List<TopoDS_Shape>&(const TopoDS_Shape&)>;
//! Fails if a face the blend made runs into a face of the body that is not beside it.
void checkBlendClear(const TopoDS_Shape&              before,
                     const TopoDS_Shape&              after,
                     const std::vector<TopoDS_Shape>& selected,
                     const BlendModified&             modified);
} // namespace BRepExt
