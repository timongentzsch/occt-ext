// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once

// The offset of a whole body as shells, hollows and solid offsets share it
// (BRepOffsetExt_Solid.cxx).

#include <Standard_Macro.hxx>
#include <BRepBooleanExt_Accept.hxx>
#include <BRepOffsetExt_Parts.hxx>
#include <BRepToolsExt_Operations.hxx>

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
                                        OffsetJoin                       join,
                                        const std::vector<TopoDS_Shape>& pushed,
                                        double                           pushedBy,
                                        Lineage&                         faces);

} // namespace BRepExt
