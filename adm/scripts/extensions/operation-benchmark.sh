#!/usr/bin/env bash
# Copyright (c) 2026 Timon Gentzsch
# SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
# Whole-operation comparison of TKAlgoExt with stock OCCT (see operation-benchmark.cxx).
#   operation-benchmark.sh NATIVE_BUILD_DIRECTORY
# NATIVE_BUILD_DIRECTORY is a static Release build of this tree with the toolkits TKAlgoExt,
# TKGeomAlgo, TKMesh, TKFillet, TKOffset, TKHLR and TKBO built (their dependencies come along).
# Requires a C++17 compiler. Prints one row per operation, input and implementation.
set -euo pipefail
build="${1:?usage: operation-benchmark.sh NATIVE_BUILD_DIRECTORY}"
here="$(cd "$(dirname "$0")" && pwd)"
mapfile -t libs < <(find "$build" -name 'libTK*.a' -print | sort)
out="$(mktemp -d)/operation-benchmark"
c++ -std=c++17 -O3 -I"$build/include/opencascade" "$here/operation-benchmark.cxx" \
  -Wl,--start-group "${libs[@]}" -Wl,--end-group -pthread -ldl -o "$out"
"$out"
