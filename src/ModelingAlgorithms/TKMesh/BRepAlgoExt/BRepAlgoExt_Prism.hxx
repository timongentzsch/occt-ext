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
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Vec.hxx>
#include <vector>

namespace BRepExt
{
// A sufficient, conservative cap-collapse test for a right prism. Only used
// after native construction failed; it cannot change a successful native blend.
bool blendConsumesPrismCap(const TopoDS_Shape&, const std::vector<TopoDS_Shape>&, double);
bool blendExceedsConvexRadius(const TopoDS_Shape&,
                              const std::vector<TopoDS_Shape>&,
                              double,
                              const std::vector<TopoDS_Shape>* supports = nullptr);
// The way into `face` from the place `t` of its edge, across the edge and along the surface.
bool supportInward(const TopoDS_Face& face,
                   const TopoDS_Edge& edge,
                   double             t,
                   const gp_Vec&      tangent,
                   gp_Vec&            inward);
} // namespace BRepExt
