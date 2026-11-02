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
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>

#include <vector>

namespace BRepExt
{

//! One edge of a loop that Patch is to fill for another operation: the face it is to be smooth
//! with to `order` (0, 1 or 2 for G0 to G2; a null face for G0).
struct PatchSide
{
  TopoDS_Edge edge;
  TopoDS_Face support;
  int         order = 0;
};

//! The face Patch makes inside the closed loop of `sides`, within `tolerance` and `g1Degrees`
//! along them (Patch's own checks, its own choice of construction). Fails as Patch does where
//! it makes none to keep.
Standard_EXPORT TopoDS_Face patchLoop(const std::vector<PatchSide>& sides,
                                      double                        tolerance,
                                      double                        g1Degrees);

} // namespace BRepExt
