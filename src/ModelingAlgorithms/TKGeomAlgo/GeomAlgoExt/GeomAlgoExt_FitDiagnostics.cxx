// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_FitDiagnostics.hxx>
#include <GeomAlgoExt_SurfaceJet.hxx>
#include <algorithm>
#include <cmath>

namespace BRepExt
{
FitDiagnostics measureSurfaceFit(const char*                      op,
                                 const occ::handle<Geom_Surface>& surface,
                                 const std::vector<FitTarget>&    targets,
                                 double                           scale,
                                 gp_XYZ                           referenceArea)
{
  if (!std::isfinite(scale) || scale <= 0 || referenceArea.SquareModulus() <= 1e-24)
    fail(op, invalid_input, "fit diagnostics require a finite scale and regular reference area");
  FitDiagnostics out;
  for (const auto& target : targets)
  {
    const auto   jet = surfaceJet(op, surface, target.u, target.v);
    const double gap = jet.p.Distance(target.position);
    const double angle =
      target.hasNormal
        ? std::acos(std::clamp(jet.n.XYZ().Dot(target.normal), -1.0, 1.0)) * 180 / M_PI
        : 0;
    if (!std::isfinite(gap) || !std::isfinite(angle))
      fail(op, no_solution, "nonfinite input residual or normal angle");
    out.constraints.push_back({target.index, gap, angle, target.hasNormal});
    out.maxResidual      = std::max(out.maxResidual, gap);
    out.maxNormalDegrees = std::max(out.maxNormalDegrees, angle);
  }
  double u0, u1, v0, v1;
  surface->Bounds(u0, u1, v0, v1);
  const double du = u1 - u0, dv = v1 - v0;
  for (int i = 0; i <= 32; ++i)
    for (int j = 0; j <= 32; ++j)
    {
      const auto   jet = surfaceJet(op, surface, u0 + du * i / 32, v0 + dv * j / 32);
      const double jac =
        jet.u.Crossed(jet.v).XYZ().Dot(referenceArea) * du * dv / referenceArea.SquareModulus();
      if (!std::isfinite(jac))
        fail(op, no_solution, "nonfinite surface differential");
      out.minOrientedJacobian = std::min(out.minOrientedJacobian, jac);
      const double energy     = (jet.uu.SquareMagnitude() * std::pow(du, 4)
                             + 2 * jet.uv.SquareMagnitude() * du * du * dv * dv
                             + jet.vv.SquareMagnitude() * std::pow(dv, 4))
                            / (scale * scale);
      if (!std::isfinite(energy))
        fail(op, no_solution, "nonfinite surface derivatives");
      const double weight = (i == 0 || i == 32 ? 0.5 : 1) * (j == 0 || j == 32 ? 0.5 : 1);
      out.bending += weight * energy / (32 * 32);
    }
  return out;
}
} // namespace BRepExt
