// Copyright (c) 2026 Timon Gentzsch
//
// This file is part of Open CASCADE Technology software library.
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult LICENSE_LGPL_21.txt for the complete license.

#pragma once
#include <Standard.hxx>
#include <GeomAlgoExt_Execution.hxx>
// What a patch is like inside its loop. The misses along its edges say nothing of that: a
// surface can meet every edge and still stand up as a spike between them, or fold over.
#include <TopoDS_Face.hxx>
#include <optional>
#include <string>

namespace BRepAlgoExt {
/// What is wrong with the inside of a face, if anything is: that it folds over, or that it
/// bends many times more tightly than anything on its rim does (`rimBend`, a curvature; zero
/// or less leaves that unjudged).
Standard_EXPORT std::optional<std::string> interiorFault(const TopoDS_Face &face, double rimBend, double uvTolerance, double maximumBendRatio, GeomAlgoExt::StageObserver observer);

struct InteriorCurvature {
    double low = 0, high = 0;
    Standard_EXPORT double wrongSign() const;
    Standard_EXPORT double tightest() const;
};
Standard_EXPORT InteriorCurvature surveyCurvature(const TopoDS_Face &face);
/// Largest normal-curvature step across a knot below C2, away from the rim.
Standard_EXPORT double surveyStep(const TopoDS_Face &face);

/// The largest principal curvature of a face at the places it is looked at inside its loop.
Standard_EXPORT double maximumInteriorCurvature(const TopoDS_Face &face, double uvTolerance);
} // namespace BRepAlgoExt
