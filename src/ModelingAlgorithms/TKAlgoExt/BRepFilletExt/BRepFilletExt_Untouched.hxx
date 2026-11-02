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
// A Boolean of two shapes that leaves both as they were. OpenCascade otherwise adds p-curves
// to its arguments and may raise their tolerances, also when it fails; and the shapes of a
// document are shared by its earlier states, by previews and by commands that are refused.
#include <NCollection_List.hxx>
#include <TopoDS_Shape.hxx>

namespace BRepExt
{
template <class Operation>
struct Untouched : Operation
{
  Untouched(const TopoDS_Shape& object, const TopoDS_Shape& tool)
  {
    NCollection_List<TopoDS_Shape> objects, tools;
    objects.Append(object);
    tools.Append(tool);
    this->SetArguments(objects);
    this->SetTools(tools);
    this->SetNonDestructive(true);
    this->Build();
  }
};
} // namespace BRepExt
