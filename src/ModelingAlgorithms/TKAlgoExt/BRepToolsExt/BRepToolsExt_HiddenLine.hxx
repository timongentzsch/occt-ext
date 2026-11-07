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
#include <BRepToolsExt_Operations.hxx>

namespace BRepExt
{
// A projected curve retains its source body and line category. Coordinates are in the view's
// plane; only their approximation to a polyline depends on the drawing tolerance.
struct DrawingCurve
{
  int                 body = 0, kind = 0;
  std::vector<double> points;
};

Standard_EXPORT std::vector<DrawingCurve> hiddenLine(OcctKernel&,
                                                     const std::vector<uint32_t>&,
                                                     Vec3                       normal,
                                                     Vec3                       x,
                                                     const std::vector<double>& options);
} // namespace BRepExt
