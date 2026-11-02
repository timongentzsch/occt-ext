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

// The offset of a whole body as shells, hollows and solid offsets share it
// (BRepOffsetExt_Solid.cxx).

#include <Standard_Macro.hxx>
#include <BRepBooleanExt_Accept.hxx>
#include <BRepOffsetExt_Parts.hxx>

namespace BRepExt
{

//! `solid` without the faces that an offset by `distance` shrinks to nothing, their neighbours
//! extended over them; `solid` itself when there are none. `faces` is carried through. Fails,
//! naming the faces, where one is too tight in part only or cannot be removed.
Standard_EXPORT TopoDS_Shape withoutVanishing(const char*         op,
                                              const TopoDS_Shape& solid,
                                              double              distance,
                                              Lineage&            faces);

//! `work` offset by `distance` (negative: inwards), with `pushed` moved by `pushedBy` instead.
//! `faces` becomes the images in the result. A null shape where OpenCascade's offset fails or
//! what it returns is not one valid solid that keeps clear of itself.
Standard_EXPORT TopoDS_Shape offsetBody(const TopoDS_Shape&              work,
                                        double                           distance,
                                        double                           tolerance,
                                        int                              join,
                                        const std::vector<TopoDS_Shape>& pushed,
                                        double                           pushedBy,
                                        Lineage&                         faces);

} // namespace BRepExt
