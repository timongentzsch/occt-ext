// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
#include <TopoDS_Shape.hxx>
namespace BRepExt::tolerance {
// Stored values, without a global application floor.
double own(const TopoDS_Shape& shape);
double largest(const TopoDS_Shape& shape);
}
