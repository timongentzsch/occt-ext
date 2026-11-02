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

// Opt-in host clock and phase sink. The host supplies executing-thread CPU milliseconds;
// without the hook the existing rung timers keep their elapsed-time clock. No logging by default.
namespace BRepExt
{
double booleanNow();
bool   booleanTimingEnabled();
void   booleanPhase(const char* name, double milliseconds);

struct BooleanPhase
{
  const char* name;
  bool        enabled;
  double      start;

  explicit BooleanPhase(const char* phase)
      : name(phase),
        enabled(booleanTimingEnabled()),
        start(enabled ? booleanNow() : 0)
  {
  }

  ~BooleanPhase()
  {
    if (enabled)
      booleanPhase(name, booleanNow() - start);
  }
};
} // namespace BRepExt
