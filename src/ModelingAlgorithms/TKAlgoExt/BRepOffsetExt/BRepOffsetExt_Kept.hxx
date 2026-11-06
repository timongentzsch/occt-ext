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
