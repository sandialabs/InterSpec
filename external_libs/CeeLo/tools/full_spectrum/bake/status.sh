#!/usr/bin/env bash
# CeeLo: a Monte Carlo photon-transport library for computing gamma-ray
# and X-ray detector efficiency - developed as part of InterSpec.
#
# Copyright 2026 National Technology & Engineering Solutions of Sandia, LLC
# (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
# Government retains certain rights in this software.
# For questions contact William Johnson via email at wcjohns@sandia.gov, or
# alternative emails of interspec@sandia.gov.
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

# Where a bake stands, from its files alone: safe to run at any time, from any
# shell, while a bake is running or after it ended.  all.sh writes this report
# to $WORK/bake_report.txt when it finishes.
#
#   tools/full_spectrum/bake/status.sh
#
# Per work area: steps done and the last one, gate verdict, and size; then the
# environment campaign's rounds, the published payloads, disk space, and the
# bake processes still running.
. "$(dirname "$0")/common.sh"

echo "Bake status $(date '+%Y-%m-%d %H:%M:%S'): WORK=$WORK RESULTS=$RESULTS SMOKE=$SMOKE PRECISION=$PRECISION"
echo "Source: $(git rev-parse --short HEAD 2>/dev/null || echo '?')$(git diff --quiet HEAD 2>/dev/null || echo ' (modified tree)')"
echo

if [ -f "$WORK/all/stages.txt" ]; then
    echo "all.sh stages:"
    sed 's/^/  /' "$WORK/all/stages.txt"
    echo
fi

"$PY" - "$WORK" <<'EOF'
import json, os, sys
from pathlib import Path

work = Path(sys.argv[1])

def gate_verdict(area: Path) -> str:
    report = area / "gates.json"
    if not report.exists():
        return "-"
    data = json.loads(report.read_text())
    passed = data.get("passed", data.get("pass"))
    if passed is None:
        return "?"
    failed = [c.get("name", "?") for c in data.get("checks", []) if c.get("passed") is False]
    return "PASS" if passed else f"FAIL ({len(failed) or len(data.get('failed', []))} checks)"

def size(area: Path) -> str:
    total = sum(f.stat().st_size for f in area.rglob("*") if f.is_file())
    return f"{total / 1e9:6.2f} GB"

print(f"{'work area':24s} {'steps':>5s}  {'last step done':38s} {'gates':16s} {'size':>9s}")
for area in sorted(p for p in work.iterdir() if (p / "done").is_dir()):
    done = sorted((area / "done").iterdir(), key=lambda f: f.stat().st_mtime)
    last = f"{done[-1].name} at {done[-1].read_text().strip()}" if done else "-"
    print(f"{area.name:24s} {len(done):5d}  {last:38s} {gate_verdict(area):16s} {size(area):>9s}")
    for out in (work / f"{area.name}.out",):
        if out.exists():
            lines = [l for l in out.read_text(errors="replace").splitlines()
                     if "_DONE" in l or "ERROR" in l or "FAILED" in l]
            if lines:
                print(f"{'':24s}        {lines[-1].strip()[:100]}")

for area in sorted(work.glob("env_*")):
    rounds = sorted(area.glob("round*/manifest.json"),
                    key=lambda p: int(p.parent.name[5:]))
    if not rounds:
        continue
    print(f"\n{area.name}: {len(rounds)} round(s) recorded")
    for manifest in rounds:
        state = json.loads(manifest.read_text())["campaign_state"]
        line = f"  {manifest.parent.name}: stopped={state['stopped']}"
        acceptance = manifest.parent / "acceptance.json"
        if acceptance.exists():
            selection = json.loads(acceptance.read_text()).get("selection", {})
            worst = [(name, s.get("max_yield_relative_error"), s.get("max_joint_tv"))
                     for name, s in selection.items() if isinstance(s, dict)]
            worst = [w for w in worst if w[1] is not None]
            if worst:
                # The campaign's targets: yield 0.5%, joint TV 1%.
                within = sum(1 for _, err, tv in worst
                             if err <= 0.005 and (tv is None or tv <= 0.01))
                name, err, tv = max(worst, key=lambda w: w[1])
                line += (f"; {within}/{len(worst)} components within 0.5% yield / 1% TV"
                         f"; worst yield error {100 * err:.2f}% ({name})")
        print(line)
EOF

echo
echo "Published payloads ($RESULTS):"
if [ -d "$RESULTS" ]; then
    for d in "$RESULTS"/*/; do
        [ -d "$d" ] || continue
        printf '  %-36s %s\n' "$(basename "$d")" "$(ls -lT "$d" | awk 'NR > 1 { print $6, $7, $8 }' | sort | tail -1)"
    done
fi

echo
echo "Disk:"
df -Ph "$WORK" | awk 'NR == 2 { print "  " $6 ": " $4 " free of " $2 }'
if find "$WORK" -name '._*' -print -quit 2>/dev/null | grep -q .; then
    echo "  WARNING: AppleDouble ._ files under $WORK"
fi

echo
echo "Bake processes running:"
ps -axo pid,ni,etime,command | grep -E 'bake/|ceelo_|env_kernel_bake|env_scatter_fit|full_spectrum/' \
    | grep -v -E 'grep|status\.sh' | cut -c1-150 | sed 's/^/  /' || true
