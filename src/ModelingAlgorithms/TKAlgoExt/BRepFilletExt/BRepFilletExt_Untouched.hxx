// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
// A Boolean of two shapes that leaves both as they were. OpenCascade otherwise adds p-curves
// to its arguments and may raise their tolerances, also when it fails; and the shapes of a
// document are shared by its earlier states, by previews and by commands that are refused.
#include <NCollection_List.hxx>
#include <TopoDS_Shape.hxx>

namespace loft {
template <class Operation> struct Untouched : Operation {
    Untouched(const TopoDS_Shape& object, const TopoDS_Shape& tool) {
        NCollection_List<TopoDS_Shape> objects, tools;
        objects.Append(object);
        tools.Append(tool);
        this->SetArguments(objects);
        this->SetTools(tools);
        this->SetNonDestructive(true);
        this->Build();
    }
};
} // namespace loft
