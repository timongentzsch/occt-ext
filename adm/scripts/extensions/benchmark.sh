#!/usr/bin/env bash
# Copyright (c) 2026 Timon Gentzsch
# SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
# Build this revision with native CMake first, then pass its build directory.
set -euo pipefail
build="$(cd "${1:?Usage: benchmark.sh BUILD_DIRECTORY}" && pwd)"
here="$(cd "$(dirname "$0")" && pwd)"
mapfile -t archives < <(find "$build" -name 'libTK*.a' -print | sort)
[[ ${#archives[@]} -gt 0 ]] || { echo 'No static OCCT archives in build directory' >&2; exit 1; }
"${CXX:-c++}" -std=c++17 -O3 -I"$build/include/opencascade" "$here/benchmark.cxx" \
  -Wl,--start-group "${archives[@]}" -Wl,--end-group -pthread -ldl -o "$build/extensions-benchmark"
"$build/extensions-benchmark"
