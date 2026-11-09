// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
// What a patch is like inside its loop. The misses along its edges say nothing of that: a
// surface can meet every edge and still stand up as a spike between them, or fold over.
#include <Standard_Macro.hxx>
#include <TopoDS_Face.hxx>
#include <optional>
#include <functional>
#include <string>

namespace BRepExt
{
//! A four-sided face of which two opposite sides are rims and the other two gaps of an open
//! boundary: it crosses from the rim at its lowest V (U, if not `acrossV`) to the rim at its
//! highest, asked to meet their faces to the orders `first` and `last` (0 to 2).
struct Crossing
{
  bool acrossV = false;
  int  first = 0, last = 0;
};

//! What is wrong with the inside of a face, if anything is: that it folds over, or that it
//! bends many times more tightly than anything on its rim does (`rimBend`, a curvature; zero
//! or less leaves that unjudged). Across a `crossing` it may bend as tightly as the bridge
//! between its rims does as well (see BRepPatchExt_Inside.cxx). The unresolved coarse
//! normal-turn mode is for candidate ordering only; final acceptance resolves it.
Standard_EXPORT std::optional<std::string> insideFault(
  const TopoDS_Face&             face,
  double                         rimBend,
  const std::optional<Crossing>& crossing = std::nullopt,
  bool                           resolveNormalRotation = true);

//! A cell wholly inside the trim; unresolved boundary overlap returns false.
Standard_EXPORT std::function<bool(double, double, double, double)> trimmedInteriorCells(
  const TopoDS_Face& face);

struct InsideCurvature
{
  double low = 0, high = 0;
  double wrongSign() const;
  double tightest() const;
};

Standard_EXPORT InsideCurvature insideCurvature(const TopoDS_Face& face);
//! Largest normal-curvature step across a knot below C2, away from the rim.
Standard_EXPORT double insideStep(const TopoDS_Face& face);

//! The largest principal curvature of a face at the places it is looked at inside its loop.
Standard_EXPORT double tightestInside(const TopoDS_Face& face);

//! The bending energy of a face inside its loop, the integral of k1^2 + k2^2 over its area, read
//! at the same places.
Standard_EXPORT double bendingInside(const TopoDS_Face& face);
} // namespace BRepExt
