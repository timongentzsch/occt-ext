// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

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
