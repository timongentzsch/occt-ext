// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
#include <Standard_Macro.hxx>
#include <BRepToolsExt_Session.hxx>
#include <BRepBlendExt_Request.hxx>

namespace BRepExt::blending
{
//! The blend of a stored body's edges (Request). Where `counts` is not empty, edges[i] takes
//! counts[i] (place, size) pairs of `sizes`, flat, in turn (Request::radii; none: one size);
//! `angle` a chamfer's (Request::angle); `reference` 0 for none.
Standard_EXPORT Blended blendStored(const OcctKernel&            kernel,
                                    const char*                  operation,
                                    uint32_t                     shape,
                                    const std::vector<uint32_t>& edges,
                                    double                       first,
                                    double                       second,
                                    Law                          law,
                                    double                       rho,
                                    uint32_t                     reference,
                                    const std::vector<int>&      counts      = {},
                                    const std::vector<double>&   sizes       = {},
                                    double                       angle       = 0,
                                    bool                         deferredFit = true);
} // namespace BRepExt::blending
