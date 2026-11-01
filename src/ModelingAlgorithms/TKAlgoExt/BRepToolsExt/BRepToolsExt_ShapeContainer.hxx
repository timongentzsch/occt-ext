// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
#include <BRep_Builder.hxx>
#include <TopoDS_Shape.hxx>

namespace BRepExt {
/// Gathering a caller-owned shape in a temporary container must not change its
/// administrative Free flag, including when OCCT rejects the insertion.
inline void addPreservingFree(const BRep_Builder& builder, TopoDS_Shape& container,
                              TopoDS_Shape child) {
    struct Restore {
        TopoDS_Shape shape;
        bool free;
        ~Restore() { shape.Free(free); }
    } restore{child, child.Free()};
    builder.Add(container, child);
}
} // namespace BRepExt
