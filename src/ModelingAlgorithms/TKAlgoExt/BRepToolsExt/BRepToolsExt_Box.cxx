// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// A box by its three sizes, refused below what the kernel can tell apart.
#include <BRepToolsExt_Checks.hxx>
#include <BRepToolsExt_Operations.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <Precision.hxx>
#include <cmath>

namespace loft {
EvolutionData makeBox(OcctKernel& k,double x,double y,double z) {
    if(!std::isfinite(x)||!std::isfinite(y)||!std::isfinite(z)||
       x<=Precision::Confusion()||y<=Precision::Confusion()||z<=Precision::Confusion())
        fail("makeBox",invalid_input,"box dimensions must exceed kernel resolution (1e-7)");
    BRepPrimAPI_MakeBox maker(x,y,z);maker.Build();
    if(!maker.IsDone()) fail("makeBox",no_solution,"box dimensions cannot be represented by the kernel");
    return stored(k, maker.Shape());
}
} // namespace loft
