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
#include <BRepTools_History.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopoDS_Shape.hxx>
#include <vector>
#include <BRepBuilderAPI_Sewing.hxx>
#include <functional>
#include <initializer_list>

namespace BRepExt
{

//! What `sewing` made of `added`, a shape given to it or a sub-shape of one: itself where kept.
TopoDS_Shape sewn(const BRepBuilderAPI_Sewing& sewing, const TopoDS_Shape& added);

//! Records in `history`, as modifications, what `sewing` made of each entity of `types` in
//! `sources` that `result` delivers. `added` names the shape given to the sewing for a source
//! entity (the entity itself by default). Returns the entities with no product in `result`.
std::vector<TopoDS_Shape> sewingHistory(
  const BRepBuilderAPI_Sewing&                            sewing,
  const std::vector<TopoDS_Shape>&                        sources,
  const TopoDS_Shape&                                     result,
  std::initializer_list<TopAbs_ShapeEnum>                 types,
  BRepTools_History&                                      history,
  const std::function<TopoDS_Shape(const TopoDS_Shape&)>& added = nullptr);

//! The solid bounded by the shell `sewing` made, oriented outward, where that shell is closed; a
//! null shape where the sewing left a sheet. Fails as `reason` where the closed shell cannot be
//! oriented or bounds no volume above `minimumVolume`.
TopoDS_Shape closedSolid(const char*                  op,
                         const BRepBuilderAPI_Sewing& sewing,
                         const char*                  reason,
                         double                       minimumVolume);
} // namespace BRepExt
