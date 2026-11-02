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

// The round or bevel of the whole outline of a flat face against the walls that stand square
// on it, past the tightest bend of that outline (BRepFilletExt_BlendFace.cxx).

#include <BRepAlgoExt_Cells.hxx>

#include <vector>

namespace BRepExt
{

//! `shape` with every edge of one flat face blended by `distance`: a round (the ball's arc) or,
//! with `bevel`, a chamfer of that setback on both sides. `edges` must be all the edges of one
//! flat face, whose walls stand square on it for at least `distance`.
//!
//! The surface is the one each place of the outline makes by itself, cut where two of them
//! meet: over a point at the distance d from the outline it lies sqrt(r^2 - (r - d)^2) - r
//! (a bevel: d - r) from the face's plane, for d < r. Where the outline bends more tightly
//! than the blend those pieces meet in creases over the outline's medial axis.
//! A closing concave round instead uses the ball envelope on the trimmed inward offset.
//! Currently built for one convex polygon of straight centre paths and planar/cylindrical
//! vertical prism walls ending at one planar mouth; other closing concave geometry refuses.
//!
//! Fails with "the planar cap: top face vanishes" in its text where nothing of the face is
//! left, and otherwise with what it cannot build.
BooleanResult faceBlend(const char*                      operation,
                        const TopoDS_Shape&              shape,
                        const std::vector<TopoDS_Shape>& edges,
                        double                           distance,
                        bool                             bevel);

} // namespace BRepExt
