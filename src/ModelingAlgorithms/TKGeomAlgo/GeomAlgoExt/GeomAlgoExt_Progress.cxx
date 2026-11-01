// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_Progress.hxx>
#include <GeomAlgoExt_Cancel.hxx>
#include <GeomAlgoExt_Phase.hxx>
namespace { thread_local Message_ProgressRange current; }
namespace BRepExt {
ProgressScope::ProgressScope(const Message_ProgressRange& range) : previous(current) { current = range; }
ProgressScope::~ProgressScope() { current = previous; }
extern "C" int sessionCancelled() { return current.UserBreak(); }
}
extern "C" void BRepExt_PatchPhase(const char*, int) {}
extern "C" void BRepExt_ReportProgress(const char*) {}
