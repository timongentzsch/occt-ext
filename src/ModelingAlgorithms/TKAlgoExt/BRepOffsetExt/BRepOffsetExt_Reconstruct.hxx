// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
#include <Standard_Macro.hxx>
#include <BRepOffsetExt_Parts.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <NCollection_Map.hxx>
#include <string>

namespace BRepExt
{
using FaceSet = NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>;

// Bounded topology-changing fallback; no shape is returned unless its support
// geometry satisfies a recognized exact construction. `faces` is set to what each face of
// the solid has become.
Standard_EXPORT TopoDS_Shape reconstructOffset(const TopoDS_Shape&              solid,
                                               const std::vector<TopoDS_Shape>& moved,
                                               double                           distance,
                                               double                           tolerance,
                                               Lineage&                         faces);

//! Push-pull through any change of topology, by a solid of its own for what the moved faces
//! sweep and a Boolean with the body (BRepOffsetExt_Faces.cxx). A null shape, and in `why` the
//! reason, where the swept solid would not have the body's own faces for its sides.
Standard_EXPORT TopoDS_Shape sweepFaces(const TopoDS_Shape&              solid,
                                        const std::vector<TopoDS_Shape>& moved,
                                        double                           distance,
                                        double                           tolerance,
                                        Lineage&                         faces,
                                        std::string&                     why,
                                        std::vector<TopoDS_Shape>&       at);
} // namespace BRepExt
