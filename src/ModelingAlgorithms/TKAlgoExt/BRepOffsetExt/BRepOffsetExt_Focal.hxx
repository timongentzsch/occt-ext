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
Standard_EXPORT FocalSweep focalSweep(const TopoDS_Face& face, double distance);
} // namespace BRepExt
