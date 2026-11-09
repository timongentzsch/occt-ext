// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once

// Features that meet the rest of a body only in hole loops of their neighbours (a hole, a pocket,
// a boss), taken off by dropping those loops: Delete Face's fast path, and how
// the blend lifts off a feature its blend runs into (BRepBlendExt_Solver.cxx).

#include <Standard_Macro.hxx>
#include <BRepTools_History.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>

#include <vector>

namespace BRepExt
{

//! Whether `wire` of `owner` (forward) is a hole of it: it winds against the outer bound in the
//! face's parameters, as ShapeAnalysis::OuterWire tells them apart.
Standard_EXPORT bool holeOf(const TopoDS_Shape& wire, const TopoDS_Face& owner);

//! The face that fills a hole loop `wire` of `owner` (forward): on the owner's surface, the loop
//! its outer bound (the owner's wire reversed).
Standard_EXPORT TopoDS_Face holeCap(const TopoDS_Face& owner, const TopoDS_Shape& wire);

//! What fillHoles made of a body: the body without the features it could take off, its
//! history, and the requested faces it left.
struct Filled
{
  TopoDS_Shape              shape;
  Handle(BRepTools_History) history;
  std::vector<TopoDS_Shape> rest;
};

//! Removes those features among `faces` that border the rest of `shape` only along whole hole
//! loops of their neighbours: each such neighbour loses those loops and nothing else changes.
//!
//! A feature here is what one removal has to take together: the requested faces joined by
//! their edges and by the loops they share. It qualifies when every edge of it is shared with
//! another of its faces, is a seam or degenerate, or lies on a loop of a face outside it that
//! is a hole of that face (it winds against the outer bound in the face's parameters, as
//! ShapeAnalysis::OuterWire tells them apart) and borders only faces of the feature. Then its
//! faces and those loops make a closed surface, and taking them off leaves the body closed.
//! What the filled loops then span (holeCap) must be clear of the rest of the body; where one
//! meets another face or another filled loop within their tolerances, the feature is left to
//! BRepAlgoAPI_Defeaturing, which intersects them.
Standard_EXPORT Filled fillHoles(const TopoDS_Shape& shape, const std::vector<TopoDS_Shape>& faces);

//! The features `faces` are of: they and the faces met across their edges, again from those,
//! except across an edge on a hole loop of the face beyond, which holds the feature. A face of
//! no such feature takes in what it is joined to, up to the whole body.
Standard_EXPORT std::vector<TopoDS_Shape> featuresOf(const TopoDS_Shape&              shape,
                                                     const std::vector<TopoDS_Shape>& faces);

} // namespace BRepExt
