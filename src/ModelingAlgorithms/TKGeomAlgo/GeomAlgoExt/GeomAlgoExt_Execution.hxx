// Copyright (c) 2026 Timon Gentzsch
//
// This file is part of Open CASCADE Technology software library.
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult LICENSE_LGPL_21.txt for the complete license.

#pragma once
#include <Standard_Failure.hxx>
#include <Message_ProgressRange.hxx>
#include <string>

//! Stable diagnostic payload; application front ends may translate the reason.
class GeomAlgoExt_Failure : public Standard_Failure {
public:
  GeomAlgoExt_Failure(const char* theReason, const char* theMessage)
    : Standard_Failure(theMessage), myReason(theReason) {}
  const char* Reason() const { return myReason.c_str(); }
private:
  std::string myReason;
};
namespace GeomAlgoExt {
inline constexpr const char* no_solution = "no_solution";
[[noreturn]] inline void fail(const char*, const char* reason, const char* message) {
  throw GeomAlgoExt_Failure(reason, message);
}
//! Optional generic stage observer; no clock or host runtime is required.
using StageObserver = void (*)(const char*, bool);
class FitPhase {
  const char* myName;
  StageObserver myObserver;
public:
  FitPhase(StageObserver observer, const char* name) : myName(name), myObserver(observer) {
    if (myObserver) myObserver(myName, true);
  }
  void finish() { if (myName) { if (myObserver) myObserver(myName, false); myName = nullptr; } }
  ~FitPhase() { finish(); }
  FitPhase(const FitPhase&) = delete;
  FitPhase& operator=(const FitPhase&) = delete;
};
inline void checkCancelled(const Message_ProgressRange& range) {
  if (range.UserBreak()) throw GeomAlgoExt_Failure("cancelled", "Operation cancelled");
}
}
