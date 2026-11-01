// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
#include <BRepTools_History.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopoDS_Shape.hxx>
#include <vector>
#include <BRepBuilderAPI_Sewing.hxx>
#include <functional>
#include <initializer_list>

namespace BRepExt {

/// What `sewing` made of `added`, a shape given to it or a sub-shape of one: itself where kept.
TopoDS_Shape sewn(const BRepBuilderAPI_Sewing& sewing, const TopoDS_Shape& added);

/// Records in `history`, as modifications, what `sewing` made of each entity of `types` in
/// `sources` that `result` delivers. `added` names the shape given to the sewing for a source
/// entity (the entity itself by default). Returns the entities with no product in `result`.
std::vector<TopoDS_Shape> sewingHistory(
    const BRepBuilderAPI_Sewing& sewing, const std::vector<TopoDS_Shape>& sources,
    const TopoDS_Shape& result, std::initializer_list<TopAbs_ShapeEnum> types,
    BRepTools_History& history,
    const std::function<TopoDS_Shape(const TopoDS_Shape&)>& added = nullptr);

/// The solid bounded by the shell `sewing` made, oriented outward, where that shell is closed; a
/// null shape where the sewing left a sheet. Fails as `reason` where the closed shell cannot be
/// oriented or bounds no volume above `minimumVolume`.
TopoDS_Shape closedSolid(const char* op, const BRepBuilderAPI_Sewing& sewing, const char* reason,
                         double minimumVolume);
}
