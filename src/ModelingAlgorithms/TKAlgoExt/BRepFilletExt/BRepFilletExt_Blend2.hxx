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
