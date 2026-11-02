// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

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

namespace patch_host
{
void chain(std::vector<Side>&, double)
{
  fail("patch", unsupported_geometry, "open boundaries need a host");
}

void closeGaps(const std::vector<Side>&, std::vector<Rim>&) {}
} // namespace patch_host
} // namespace BRepExt
