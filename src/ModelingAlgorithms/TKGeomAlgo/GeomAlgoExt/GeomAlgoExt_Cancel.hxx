// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
#include <GeomAlgoExt_Failure.hxx>
#include <Message_ProgressIndicator.hxx>
#include <Message_ProgressScope.hxx>
namespace loft {
extern "C" int sessionCancelled();
/// Host hook receiving the name of each progress scope OCCT opens (diagnostics only).
extern "C" void loft_report_progress(const char* name);
inline void checkCancelled(const char* operation, bool enabled = true) {
    if (enabled && sessionCancelled()) fail(operation, "cancelled", "Operation cancelled");
}

/// A progress indicator that stops an OCCT algorithm when the host's cancel flag is raised.
class CancelProgress : public Message_ProgressIndicator {
  protected:
    bool UserBreak() override { return sessionCancelled(); }
    void Show(const Message_ProgressScope& scope, const bool) override { loft_report_progress(scope.Name()); }
};
inline Handle(Message_ProgressIndicator) progress(bool enabled) {
    if (!enabled) return {};
    return new CancelProgress;
}
}
