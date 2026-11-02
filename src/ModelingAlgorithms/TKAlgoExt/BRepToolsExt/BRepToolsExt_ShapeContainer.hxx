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
#include <BRep_Builder.hxx>
#include <TopoDS_Shape.hxx>

namespace BRepExt
{
//! Gathering a caller-owned shape in a temporary container must not change its
//! administrative Free flag, including when OCCT rejects the insertion.
inline void addPreservingFree(const BRep_Builder& builder,
                              TopoDS_Shape&       container,
                              TopoDS_Shape        child)
{
  struct Restore
  {
    TopoDS_Shape shape;
    bool         free;

    ~Restore() { shape.Free(free); }
  } restore{child, child.Free()};

  builder.Add(container, child);
}
} // namespace BRepExt
