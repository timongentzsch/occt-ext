#!/usr/bin/env python3
# Copyright (c) 2026 Timon Gentzsch
# SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0
"""Draw the output of operation-benchmark as a grouped bar chart (SVG, no dependencies).

    operation-benchmark.sh BUILD > results.txt
    operation-chart.py results.txt dox/extensions/operation-benchmark.svg

One group per input; two bars, stock and extension, of the median CPU time over ten runs.
"stock, then its result checked" adds the BRepCheck pass that the branch performs on every result.
A bar is hatched and labelled when the result is refused or invalid, since its time then
measures a call that did not do the work.
"""
import re
import sys
from xml.sax.saxutils import escape

rows = []
BARS = (("stock", "#6b7a8f"), ("stock+check", "#c9a227"), ("extension", "#2f8f5b"))
for line in open(sys.argv[1], encoding="utf-8"):
    parts = [p.strip() for p in line.split("|")]
    if len(parts) < 6 or parts[0] in ("operation", ""):
        continue
    op, inp, who, runs, cpu, quality = parts[0], parts[1], parts[2], parts[3], parts[4], "|".join(parts[5:])
    median = float(cpu.split("[")[0])
    ok = int(runs.split("/")[0])
    bad = ok == 0 or "INVALID" in quality or "REFUSED" in quality
    # A result that is valid but wrong by more than 1% of the input's scale is also marked.
    m = re.search(r"\|dV\|\s*([0-9.e+-]+)", quality)
    wrong = bool(m and float(m.group(1)) > 1.0)
    rows.append((f"{op}: {inp}", who, median, bad or wrong))

groups = []
for label, who, median, flagged in rows:
    if not groups or groups[-1][0] != label:
        groups.append((label, {}))
    groups[-1][1][who] = (median, flagged)

W, left, bar, gap, top = 980, 330, 16, 26, 66
H = top + len(groups) * (3 * bar + gap) + 40
peak = max(m for _, g in groups for m, _ in g.values())
scale = (W - left - 90) / peak
out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="Helvetica, Arial, sans-serif" font-size="12">',
       '<defs><pattern id="hatch" width="6" height="6" patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><rect width="6" height="6" fill="#d98a8a"/><line x1="0" y1="0" x2="0" y2="6" stroke="#8b1a1a" stroke-width="2"/></pattern></defs>',
       f'<rect width="{W}" height="{H}" fill="#ffffff"/>',
       f'<text x="{left}" y="22" font-size="14" font-weight="bold">CPU time per call, median of ten runs (ms, lower is better)</text>',
       f'<rect x="{left}" y="34" width="12" height="12" fill="#6b7a8f"/><text x="{left+16}" y="45">stock OCCT</text>',
       f'<rect x="{left+100}" y="34" width="12" height="12" fill="#c9a227"/><text x="{left+116}" y="45">stock, then its result checked</text>',
       f'<rect x="{left+310}" y="34" width="12" height="12" fill="#2f8f5b"/><text x="{left+326}" y="45">this branch</text>',
       f'<rect x="{left+410}" y="34" width="12" height="12" fill="url(#hatch)"/><text x="{left+426}" y="45">no valid result, or a wrong one</text>']
y = top
for label, g in groups:
    out.append(f'<text x="{left-10}" y="{y + bar + 4}" text-anchor="end">{escape(label)}</text>')
    for k, (who, color) in enumerate(BARS):
        if who not in g:
            continue
        median, flagged = g[who]
        width = max(median * scale, 1)
        fill = "url(#hatch)" if flagged else color
        yy = y + k * bar
        out.append(f'<rect x="{left}" y="{yy}" width="{width:.1f}" height="{bar-2}" fill="{fill}"/>')
        out.append(f'<text x="{left + width + 6:.1f}" y="{yy + bar - 5}">{median:.2f}</text>')
    y += 3 * bar + gap
out.append('</svg>')
open(sys.argv[2], "w", encoding="utf-8").write("\n".join(out) + "\n")
