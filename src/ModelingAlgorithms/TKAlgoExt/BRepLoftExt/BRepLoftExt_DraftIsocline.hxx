// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once

#include <BRepOffsetExt_Parts.hxx>

#include <gp_Pln.hxx>

namespace loft {

/// `faces` of `solid`, walls along `pull`, each replaced by the ruled surface that makes
/// `angle` with `pull` and passes through the wall where the neutral plane crosses it
/// (draft_isocline.cpp). `lineage` is carried through. Fails, naming the face, where a face is
/// no such wall, where the taper folds it, and where the draft would leave a step.
TopoDS_Shape draftWalls(const TopoDS_Shape &solid, const std::vector<TopoDS_Shape> &faces, double angle,
                        const gp_Dir &pull, const gp_Pln &neutral, double tolerance, Lineage &lineage);

} // namespace loft
