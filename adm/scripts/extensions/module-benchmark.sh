#!/usr/bin/env bash
# Copyright (c) 2026 Timon Gentzsch
# SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
# Runs one module benchmark against this tree and against a stock OCCT build, to compare them.
#   module-benchmark.sh FAMILY FORK_BUILD STOCK_BASE_BUILD
# FAMILY is spline-fit, patch-chart or fit-aux (the file FAMILY-benchmark.cxx). Both builds are
# native static Release builds; the stock one is of the base commit b8f597c6.
set -euo pipefail

family="${1:?Usage: module-benchmark.sh FAMILY FORK_BUILD STOCK_BASE_BUILD}"
case "$family" in
  spline-fit | patch-chart | fit-aux) ;;
  *)
    echo 'Unknown benchmark family' >&2
    exit 1
    ;;
esac
candidate="$(cd "${2:?Supply native static Release fork build}" && pwd)"
stock="$(cd "${3:?Supply separate stock b8f597c6 Release build}" && pwd)"
here="$(cd "$(dirname "$0")" && pwd)"

for variant in candidate stock; do
  build="${!variant}"
  flags=()
  [[ "$variant" != stock ]] || flags+=(-DOCCT_EXTENSION_STOCK_ONLY)
  mapfile -t archives < <(find "$build" -name 'libTK*.a' -print | sort)
  "${CXX:-c++}" -std=c++17 -O3 "${flags[@]}" -I"$build/include/opencascade" \
    "$here/$family-benchmark.cxx" \
    -Wl,--start-group "${archives[@]}" -Wl,--end-group -pthread -ldl \
    -o "$build/$family-benchmark-$variant"
  "$build/$family-benchmark-$variant"
done
