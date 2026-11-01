// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
#include <GeomAlgoExt_Failure.hxx>
namespace loft {
extern "C" int sessionCancelled();
inline void checkCancelled(const char* operation, bool enabled = true) {
    if (enabled && sessionCancelled()) fail(operation, "cancelled", "Operation cancelled");
}
}
