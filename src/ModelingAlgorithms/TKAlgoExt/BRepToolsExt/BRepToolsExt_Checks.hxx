// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once

// Typed failures and input checks shared by every operation.

#include <Standard_Macro.hxx>
#include <BRepToolsExt_Session.hxx>
#include <GeomAlgoExt_Failure.hxx>

#include <TopAbs_ShapeEnum.hxx>
#include <TopoDS_Shape.hxx>

#include <stdexcept>
#include <string>
#include <vector>

namespace BRepExt
{

//! What OpenCascade's checker finds wrong with `shape`, by kind of sub-shape and its own status,
//! each also as found within each shape that holds it ("4:27 5in4:28"); empty when valid.
Standard_EXPORT std::string invalidity(const TopoDS_Shape& shape);

//! The same for the result of an operation on `earlier`, looking only at what the operation made:
//! the faces that `earlier` does not hold, then each solid and shell as a whole. Faces that the
//! result shares with `earlier` keep the verdict they had.
Standard_EXPORT std::string invalidity(const TopoDS_Shape& shape, const TopoDS_Shape& earlier);

//! The sub-shapes of `owner` behind `ids`, as `owner` holds them (its orientation, not the
//! handle's). Fails with invalid_input when an id is not a `type` sub-shape of `owner`.
Standard_EXPORT std::vector<TopoDS_Shape> subShapes(const OcctKernel&            kernel,
                                                    const char*                  operation,
                                                    const TopoDS_Shape&          owner,
                                                    const std::vector<uint32_t>& ids,
                                                    TopAbs_ShapeEnum             type);

} // namespace BRepExt
