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

#include <GeomAlgoExt_Progress.hxx>
#include <GeomAlgoExt_Cancel.hxx>
#include <GeomAlgoExt_Phase.hxx>

namespace
{
thread_local Message_ProgressRange current;
}

namespace BRepExt
{
ProgressScope::ProgressScope(const Message_ProgressRange& range)
    : previous(current)
{
  current = range;
}

ProgressScope::~ProgressScope()
{
  current = previous;
}

extern "C" int sessionCancelled()
{
  return current.UserBreak();
}
} // namespace BRepExt

extern "C" void BRepExt_PatchPhase(const char*, int) {}

extern "C" void BRepExt_ReportProgress(const char*) {}
