#!/usr/bin/env bash
# Copyright (c) 2026 Timon Gentzsch
# SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
set -euo pipefail
build="$(cd "${1:?Usage: utilities.sh BUILD_DIRECTORY}" && pwd)"
here="$(cd "$(dirname "$0")" && pwd)"
mapfile -t archives < <(find "$build" -name 'libTK*.a' -print | sort)
"${CXX:-c++}" -std=c++17 -O3 -I"$build/include/opencascade" "$here/utilities.cxx" \
  -Wl,--start-group "${archives[@]}" -Wl,--end-group -pthread -ldl -o "$build/extensions-utilities"
"$build/extensions-utilities"
