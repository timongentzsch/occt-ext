#!/usr/bin/env bash
# Copyright (c) 2026 Timon Gentzsch
# SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
# Configure/build native Release static TKGeomAlgo in this revision and stock b8f597c6.
set -euo pipefail
candidate="$(cd "${1:?Usage: spline-fit-benchmark.sh FORK_BUILD STOCK_BASE_BUILD}" && pwd)"
stock="$(cd "${2:?Supply a separate stock b8f597c6 Release build}" && pwd)"
here="$(cd "$(dirname "$0")" && pwd)"
for variant in candidate stock; do
 build="${!variant}"
 flags=();[[ "$variant" != stock ]] || flags+=(-DOCCT_EXTENSION_STOCK_ONLY)
 mapfile -t archives < <(find "$build" -name 'libTK*.a' -print | sort)
 "${CXX:-c++}" -std=c++17 -O3 "${flags[@]}" -I"$build/include/opencascade" "$here/spline-fit-benchmark.cxx" \
  -Wl,--start-group "${archives[@]}" -Wl,--end-group -pthread -ldl -o "$build/spline-fit-benchmark-$variant"
 "$build/spline-fit-benchmark-$variant"
done
