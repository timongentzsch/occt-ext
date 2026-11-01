// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
#include <Message_ProgressRange.hxx>
namespace loft {
// Scope a native entry call; existing algorithms poll the direct link-time hook.
class ProgressScope {
    Message_ProgressRange previous;
public:
    ProgressScope(const Message_ProgressRange& range);
    ~ProgressScope();
    ProgressScope(const ProgressScope&) = delete;
    ProgressScope& operator=(const ProgressScope&) = delete;
};
}
