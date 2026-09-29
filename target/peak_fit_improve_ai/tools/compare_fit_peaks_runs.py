#!/usr/bin/env python3
"""Compare two fit_peaks_corpus_eval run directories.

Usage: compare_fit_peaks_runs.py RUN_A RUN_B [--peaks] [--top N]

Prints the aggregate deltas, the problems whose raw cost moved most (both directions), new and
fixed mechanical failures, and, with --peaks, the reference peaks whose verdict changed.
"""
import csv
import os
import sys


def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return 0.0


def main(argv):
    if len(argv) < 3:
        print(__doc__)
        return 2
    run_a, run_b = argv[1], argv[2]
    show_peaks = "--peaks" in argv
    top = 15
    if "--top" in argv:
        top = int(argv[argv.index("--top") + 1])

    a = {r["id"]: r for r in read_tsv(os.path.join(run_a, "per_problem.tsv"))}
    b = {r["id"]: r for r in read_tsv(os.path.join(run_b, "per_problem.tsv"))}
    common = sorted(set(a) & set(b))
    only_a = sorted(set(a) - set(b))
    only_b = sorted(set(b) - set(a))

    cols = ["raw_cost", "cost_missed", "cost_extra", "cost_share", "cost_family", "cost_extent",
            "cost_area", "cost_mean", "cost_failure", "missed_strong", "missed_moderate", "missed_weak",
            "extra_significant", "extra_weak", "ghost", "share_disagree", "pairs", "family_disagree",
            "rois_compared", "extent_gt1fwhm_sides", "mechanical_failure", "nondeterministic",
            "dev_check_failures", "wall_seconds"]

    print(f"# {run_a} -> {run_b}")
    print(f"problems: {len(common)} common, {len(only_a)} only in A, {len(only_b)} only in B\n")
    print("| term | A | B | delta |")
    print("|---|---:|---:|---:|")
    for c in cols:
        sa = sum(fnum(a[i].get(c)) for i in common)
        sb = sum(fnum(b[i].get(c)) for i in common)
        print(f"| {c} | {sa:.3f} | {sb:.3f} | {sb - sa:+.3f} |")

    deltas = [(fnum(b[i]["raw_cost"]) - fnum(a[i]["raw_cost"]), i) for i in common]
    deltas.sort()
    print(f"\n## Improved most (top {top})")
    for d, i in deltas[:top]:
        if d < 0:
            print(f"- {i}: {fnum(a[i]['raw_cost']):.3f} -> {fnum(b[i]['raw_cost']):.3f} ({d:+.3f})")
    print(f"\n## Regressed most (top {top})")
    for d, i in reversed(deltas[-top:]):
        if d > 0:
            print(f"- {i}: {fnum(a[i]['raw_cost']):.3f} -> {fnum(b[i]['raw_cost']):.3f} ({d:+.3f})")

    new_fail = [i for i in common if a[i]["mechanical_failure"] == "0" and b[i]["mechanical_failure"] == "1"]
    fixed_fail = [i for i in common if a[i]["mechanical_failure"] == "1" and b[i]["mechanical_failure"] == "0"]
    if new_fail:
        print("\n## New mechanical failures\n- " + "\n- ".join(f"{i}: {b[i]['error']}" for i in new_fail))
    if fixed_fail:
        print("\n## Fixed mechanical failures\n- " + "\n- ".join(fixed_fail))

    if show_peaks:
        def truth_rows(run):
            rows = {}
            for r in read_tsv(os.path.join(run, "per_peak.tsv")):
                if r["set"] == "truth":
                    rows[(r["id"], round(fnum(r["energy"]), 2))] = r
            return rows
        pa, pb = truth_rows(run_a), truth_rows(run_b)
        changed = []
        for key in sorted(set(pa) & set(pb)):
            va, vb = pa[key]["verdict"], pb[key]["verdict"]
            if va != vb:
                changed.append((key, va, vb, pa[key]["z_det"], pa[key]["source"]))
        print(f"\n## Reference peaks whose verdict changed ({len(changed)})")
        for (pid, e), va, vb, z, src in changed:
            print(f"- {pid} {e:.2f} keV {src} z={fnum(z):.1f}: {va} -> {vb}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
