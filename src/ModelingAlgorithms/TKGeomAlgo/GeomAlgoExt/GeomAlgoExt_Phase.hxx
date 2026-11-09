// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once

// Optional phase observer supplied by the embedding program; native default is a no-op.
extern "C" void BRepExt_PatchPhase(const char* name, int begin);

namespace BRepExt
{
struct PatchPhase
{
  const char* name;

  explicit PatchPhase(const char* name)
      : name(name)
  {
    BRepExt_PatchPhase(name, 1);
  }

  void finish()
  {
    if (name)
    {
      BRepExt_PatchPhase(name, 0);
      name = nullptr;
    }
  }

  ~PatchPhase() { finish(); }

  PatchPhase(const PatchPhase&)            = delete;
  PatchPhase& operator=(const PatchPhase&) = delete;
};
} // namespace BRepExt
