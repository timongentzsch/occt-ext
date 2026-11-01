// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once


// Optional phase observer supplied by the application; native default is a no-op.
extern "C" void loft_patch_phase(const char *name, int begin);

namespace loft {
struct PatchPhase {
    const char *name;
    explicit PatchPhase(const char *name) : name(name) { loft_patch_phase(name, 1); }
    void finish() { if (name) { loft_patch_phase(name, 0); name = nullptr; } }
    ~PatchPhase() { finish(); }
    PatchPhase(const PatchPhase &) = delete;
    PatchPhase &operator=(const PatchPhase &) = delete;
};
} // namespace loft
