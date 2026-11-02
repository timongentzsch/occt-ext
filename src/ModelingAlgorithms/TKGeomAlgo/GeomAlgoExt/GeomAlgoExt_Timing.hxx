// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

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
