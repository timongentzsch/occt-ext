// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once

#include <Standard_Macro.hxx>
#include <BRepOffsetExt_Reconstruct.hxx>

#include <string>

namespace BRepExt
{

//! `body` with the faces of `moved` moved by `distance` along their outward normals and its
//! topology kept where possible: each face on its own offset surface, the edges and vertices
//! where those surfaces meet. After that construction refuses a separated tangent junction,
//! a single-shell solid may split that junction and close it with a ruled contact strip.
//! `faces` includes strips descended from each incident face. A null shape, and in `why` the
//! reason, where neither construction passes the geometric and topology checks.
Standard_EXPORT TopoDS_Shape keptOffset(const TopoDS_Shape& body,
                                        const FaceSet&      moved,
                                        double              distance,
                                        double              tolerance,
                                        Lineage&            faces,
                                        std::string&        why);

} // namespace BRepExt
