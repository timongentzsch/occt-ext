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

#include <GeomAlgoExt_Failure.hxx>
#include <TopTools_ShapeMapHasher.hxx>

namespace BRepExt
{
const char* failureKindName(FailureKind kind)
{
  switch (kind)
  {
    case FailureKind::Construction:
      return "construction";
    case FailureKind::NativeBuilder:
      return "native_builder";
    case FailureKind::TangentEdge:
      return "tangent_edge";
    case FailureKind::UnsupportedRail:
      return "unsupported_rail";
    case FailureKind::InvalidTopology:
      return "invalid_topology";
    case FailureKind::CapConsumed:
      return "cap_consumed";
    case FailureKind::ClosingConcave:
      return "closing_concave";
    case FailureKind::Overflow:
      return "overflow";
    case FailureKind::SupportContact:
      return "support_contact";
    case FailureKind::Tangency:
      return "tangency";
    case FailureKind::Curvature:
      return "curvature";
    case FailureKind::RadiusEndpoints:
      return "radius_endpoints";
    case FailureKind::RadiusLaw:
      return "radius_law";
    case FailureKind::MaterialBound:
      return "material_bound";
  }
  return "construction";
}

void fail(const char* operation, const char* reason, const std::string& message)
{
  fail(operation, reason, message, {});
}

int failureHash(const TopoDS_Shape& shape)
{
  return static_cast<int>(TopTools_ShapeMapHasher{}(shape) % static_cast<size_t>(failureHashBound));
}

void fail(const char*                      operation,
          const char*                      reason,
          const std::string&               message,
          const std::vector<TopoDS_Shape>& at,
          const char*                      code,
          FailureKind                      kind)
{
  std::vector<int> hashes;
  for (const TopoDS_Shape& shape : at)
    hashes.push_back(failureHash(shape));
  Failure failure(operation, reason, message, std::move(hashes), code, kind);
  failure.shapes = at;
  throw failure;
}

} // namespace BRepExt
