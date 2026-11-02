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
