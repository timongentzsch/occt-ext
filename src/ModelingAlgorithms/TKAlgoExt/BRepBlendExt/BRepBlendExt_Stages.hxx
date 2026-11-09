// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once

// Blend: what the solver's entry (BRepBlendExt_Solver.cxx) and its stages
// (BRepBlendExt_Stages.cxx) share.

#include <BRepBlendExt_Internal.hxx>

namespace BRepExt::blending
{

//! The joins carried through a stage's history: along modifications of their edge only.
std::vector<Join> carriedJoins(const std::vector<Join>&         joins,
                               const Handle(BRepTools_History)& history,
                               const TopoDS_Shape&              result);

//! How the faces beside an edge turn along it: -1 convex, 1 concave, 0 tangent for `size`.
int convexity(const char*        op,
              const TopoDS_Edge& edge,
              const TopoDS_Face& first,
              const TopoDS_Face& second,
              double             size,
              bool*              touches = nullptr);

//! The edges asked, with those tangent to them added (only those in `mayAdd` where it is given).
std::vector<TopoDS_Shape> withTangent(const TopoDS_Shape&              body,
                                      std::vector<TopoDS_Shape>        edges,
                                      double                           size,
                                      const std::vector<TopoDS_Shape>& mayAdd = {});

//! The blends of a request on a copy of its body, stage by stage (struct Blend); the edges whose
//! balls are found through their centres in `byCentre`.
Blended blendCopy(const Request&                   request,
                  BRepBuilderAPI_Copy&             copy,
                  const std::vector<TopoDS_Shape>& byCentre = {});

} // namespace BRepExt::blending
