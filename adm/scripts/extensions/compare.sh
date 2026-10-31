#!/usr/bin/env bash
# Copyright (c) 2026 Timon Gentzsch
# SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
set -euo pipefail
base="${1:?Usage: compare.sh BASE CANDIDATE OUTPUT_DIRECTORY}"
candidate="${2:?Missing candidate revision}"
mkdir -p "${3:?Missing output directory}"
out="$(cd "$3" && pwd)"
here="$(cd "$(dirname "$0")" && pwd)"
repo="$(git -C "$here" rev-parse --show-toplevel)"
for label in base candidate; do
  revision="$base"
  [[ "$label" != candidate ]] || revision="$candidate"
  git -C "$repo" worktree add --detach "$out/$label-source" "$revision"
  cmake -S "$out/$label-source" -B "$out/$label-build" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release -DBUILD_LIBRARY_TYPE=Static \
    -DBUILD_MODULE_Visualization=OFF -DBUILD_MODULE_DataExchange=OFF \
    -DBUILD_MODULE_ApplicationFramework=OFF -DBUILD_MODULE_Draw=OFF \
    -DUSE_TBB=OFF -DUSE_TK=OFF -DBUILD_USE_PCH=ON
  cmake --build "$out/$label-build" --parallel "${JOBS:-6}"
  "$here/benchmark.sh" "$out/$label-build" > "$out/$label-cpu.jsonl"
done
# Keep source/build trees for inspection. The caller can remove them with git worktree remove.
