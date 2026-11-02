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

#include <BRepOffsetExt_Parts.hxx>

#include <gp_Pln.hxx>

namespace BRepExt
{

//! `faces` of `solid`, walls along `pull`, each replaced by the ruled surface that makes
//! `angle` with `pull` and passes through the wall where the neutral plane crosses it
//! (BRepSkinExt_DraftIsocline.cxx). `lineage` is carried through. Fails, naming the face, where a
//! face is no such wall, where the taper folds it, and where the draft would leave a step.
TopoDS_Shape draftWalls(const TopoDS_Shape&              solid,
                        const std::vector<TopoDS_Shape>& faces,
                        double                           angle,
                        const gp_Dir&                    pull,
                        const gp_Pln&                    neutral,
                        double                           tolerance,
                        Lineage&                         lineage);

} // namespace BRepExt
