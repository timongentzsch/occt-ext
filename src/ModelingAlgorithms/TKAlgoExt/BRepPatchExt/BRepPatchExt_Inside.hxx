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
// What a patch is like inside its loop. The misses along its edges say nothing of that: a
// surface can meet every edge and still stand up as a spike between them, or fold over.
#include <Standard_Macro.hxx>
#include <TopoDS_Face.hxx>
#include <optional>
#include <string>

namespace BRepExt
{
//! What is wrong with the inside of a face, if anything is: that it folds over, or that it
//! bends many times more tightly than anything on its rim does (`rimBend`, a curvature; zero
//! or less leaves that unjudged).
Standard_EXPORT std::optional<std::string> insideFault(const TopoDS_Face& face, double rimBend);

struct InsideCurvature
{
  double low = 0, high = 0;
  double wrongSign() const;
  double tightest() const;
};

Standard_EXPORT InsideCurvature insideCurvature(const TopoDS_Face& face);
//! Largest normal-curvature step across a knot below C2, away from the rim.
Standard_EXPORT double insideStep(const TopoDS_Face& face);

//! The largest principal curvature of a face at the places it is looked at inside its loop.
Standard_EXPORT double tightestInside(const TopoDS_Face& face);
} // namespace BRepExt
