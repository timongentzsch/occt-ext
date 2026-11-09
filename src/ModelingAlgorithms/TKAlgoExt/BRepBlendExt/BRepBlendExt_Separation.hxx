// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Blend: a certificate that two faces are farther apart than a tolerance.

#pragma once
#include <Standard_Macro.hxx>
#include <TopoDS_Face.hxx>

namespace BRepExt::blending::detail
{
//! Whether every point of `a`, its edges and its vertices is proven farther than `tolerance`
//! from every point of `b`, its edges and its vertices (BRepBlendExt_Separation.cxx). False
//! means not proven: the faces may still be apart, and an exact distance has to decide.
Standard_EXPORT bool separated(const TopoDS_Face& a, const TopoDS_Face& b, double tolerance);
} // namespace BRepExt::blending::detail
