// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

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
