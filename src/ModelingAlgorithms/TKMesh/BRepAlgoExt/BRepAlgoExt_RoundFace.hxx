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

// A face that is a constant round: a cylinder or torus tangent to two supports along its edges,
// as Measure's fillet dimension, Refillet and the blend's re-blending recognise it.

#include <Standard_Macro.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Torus.hxx>

#include <vector>

namespace BRepExt
{

struct RoundFace
{
  bool                     valid = false;
  gp_Cylinder              cylinder;
  gp_Torus                 torus;
  bool                     toroidal = false;
  std::vector<TopoDS_Face> supports;

  //! The round's radius: the cylinder's, or the torus's minor one.
  double radius() const { return toroidal ? torus.MinorRadius() : cylinder.Radius(); }
};

//! `face` of `shape` as a round: valid where it is a cylinder (not closed round) or a torus that
//! is tangent, to `tolerance`, along every edge it shares with exactly two neighbours, its
//! supports.
Standard_EXPORT RoundFace recognizeRound(const TopoDS_Shape& shape,
                                         const TopoDS_Face&  face,
                                         double              tolerance);

} // namespace BRepExt
