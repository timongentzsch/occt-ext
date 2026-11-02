// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
#include <BRepToolsExt_Session.hxx>
#include <BRepBlendExt_Request.hxx>

namespace BRepExt::blend2
{
Blended blendStored(const OcctKernel&            kernel,
                    const char*                  operation,
                    uint32_t                     shape,
                    const std::vector<uint32_t>& edges,
                    double                       first,
                    double                       second,
                    int                          law,
                    double                       rho,
                    uint32_t                     reference);
}
