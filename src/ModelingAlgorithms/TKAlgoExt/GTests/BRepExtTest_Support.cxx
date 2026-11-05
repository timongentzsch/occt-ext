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

// The host seam of BRepToolsExt_Session.hxx for the tests: shapes live in a vector and a handle
// is the index plus one.

#include "BRepExtTest_Support.hxx"

#include <BRepPatchExt_Side.hxx>
#include <GeomAlgoExt_Failure.hxx>

#include <ctime>

namespace
{
std::vector<TopoDS_Shape> shapes;
alignas(8) char token[8];
} // namespace

OcctKernel& BRepExtTest::session()
{
  return *reinterpret_cast<OcctKernel*>(token);
}

size_t BRepExtTest::storedShapeCount() { return shapes.size(); }

namespace BRepExt
{
uint32_t store(OcctKernel&, const TopoDS_Shape& shape)
{
  shapes.push_back(shape);
  return static_cast<uint32_t>(shapes.size());
}

const TopoDS_Shape& get(const OcctKernel&, uint32_t id)
{
  if (id == 0 || id > shapes.size())
    fail("shape", invalid_input, "Invalid shape reference");
  return shapes[id - 1];
}

void release(OcctKernel&, uint32_t) {}

double booleanNow()
{
  return 1000. * std::clock() / CLOCKS_PER_SEC;
}

bool booleanTimingEnabled()
{
  return false;
}

void booleanPhase(const char*, double) {}

} // namespace BRepExt
