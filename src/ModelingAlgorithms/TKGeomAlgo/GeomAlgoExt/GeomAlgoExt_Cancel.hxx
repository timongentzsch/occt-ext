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
#include <GeomAlgoExt_Failure.hxx>
#include <Message_ProgressIndicator.hxx>
#include <Message_ProgressScope.hxx>

namespace BRepExt
{
extern "C" int sessionCancelled();
//! Host hook receiving the name of each progress scope OCCT opens (diagnostics only).
extern "C" void BRepExt_ReportProgress(const char* name);

inline void checkCancelled(const char* operation, bool enabled = true)
{
  if (enabled && sessionCancelled())
    fail(operation, "cancelled", "Operation cancelled");
}

//! A progress indicator that stops an OCCT algorithm when the host's cancel flag is raised.
class CancelProgress : public Message_ProgressIndicator
{
protected:
  bool UserBreak() override { return sessionCancelled(); }

  void Show(const Message_ProgressScope& scope, const bool) override
  {
    BRepExt_ReportProgress(scope.Name());
  }
};

inline Handle(Message_ProgressIndicator) progress(bool enabled)
{
  if (!enabled)
    return {};
  return new CancelProgress;
}
} // namespace BRepExt
