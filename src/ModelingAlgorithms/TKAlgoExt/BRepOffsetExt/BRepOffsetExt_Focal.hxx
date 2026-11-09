// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
#pragma once
#include <TopoDS_Face.hxx>

namespace BRepExt
{
//! A local normal sweep proof, separate from global self-intersection.
//! Unknown includes unsupported supports, singular parameters and work-budget exhaustion.
//! A loose enclosure is never a crossing witness.
enum class FocalSweep
{
  Regular,
  Crosses,
  Unknown
};
//! Whether moving `face` along its normal by `distance` folds it, proven by interval bounds.
Standard_EXPORT FocalSweep focalSweep(const TopoDS_Face& face, double distance);
} // namespace BRepExt
