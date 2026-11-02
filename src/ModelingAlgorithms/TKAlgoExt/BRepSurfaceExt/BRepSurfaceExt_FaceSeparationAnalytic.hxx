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

// The analytic branch of the self-interference pre-pass (BRepSurfaceExt_FaceSeparation.hxx): faces
// of planes beside faces of surfaces of revolution about an axis normal to them, proven to meet
// only along the edges they share (BRepSurfaceExt_FaceSeparationAnalytic.cxx).

#include <NCollection_IndexedMap.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>

namespace BRepExt
{

//! Whether two faces that have `shared` in common (edges and vertices) meet only there and
//! come within their tolerances of each other only beside it. False where this cannot tell.
bool analyticBeside(const TopoDS_Face&                                                   a,
                    const TopoDS_Face&                                                   b,
                    const NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>& shared);

//! Whether `shape` has two faces beside each other that analyticBeside is for.
bool analyticWorth(const TopoDS_Shape& shape);

} // namespace BRepExt
