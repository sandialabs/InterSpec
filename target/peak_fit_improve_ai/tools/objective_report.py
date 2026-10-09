#!/usr/bin/env python3
"""Tables from a peak_fit_objective_eval run (summary.tsv / paired.tsv); standard library only.

  objective_report.py RUN_DIR                      # per-problem table, every objective side by side
  objective_report.py RUN_DIR --by area,cont_level # aggregate (median over problems) by dimensions
  objective_report.py RUN_DIR --filter det=HPGe --filter suite=primary
  objective_report.py RUN_DIR --cols rel_bias_truth,pull_mean,pull_rms,cov1,fail_rate
  objective_report.py RUN_DIR --paired             # paired-difference (bad-fit tail) counts per objective
  objective_report.py RUN_A --compare RUN_B        # same objective, two runs (e.g. before/after a change)

Aggregation over problems uses the median (robust to the odd pathological configuration); counts
(n, n_ok) are summed and fail_rate is recomputed from them.  --zbins groups by the truth z
(S/sqrt(S+B)) class instead of / as well as by a dimension: use --by z_class.
"""
import argparse
import csv
import math
import os
import sys
from collections import OrderedDict, defaultdict

DEFAULT_COLS = ['rel_bias_truth', 'rel_bias_pseudo', 'pull_mean', 'pull_rms', 'cov1',
                'unc_ratio', 'tail_rate', 'fail_rate', 'fwhm_rel_bias', 'cpu_median']

SHORT = {
    'rel_bias_truth': 'bias', 'rel_bias_pseudo': 'biasP', 'rel_bias_truth_se': 'biasSE',
    'pull_mean': 'pullM', 'pull_rms': 'pullR', 'cov1': 'cov1', 'cov2': 'cov2', 'unc_ratio': 'uncR',
    'tail_rate': 'tail', 'fail_rate': 'fail', 'fwhm_rel_bias': 'fwhmB', 'fwhm_pull_mean': 'fwPullM',
    'fwhm_pull_rms': 'fwPullR', 'mean_bias_sigma': 'meanB', 'mean_pull_mean': 'mPullM',
    'mean_pull_rms': 'mPullR', 'cpu_median': 'cpuMed', 'cpu_p95': 'cpu95', 'cpu_mean': 'cpuAvg',
}


def fnum(s):
    try:
        v = float(s)
    except (TypeError, ValueError):
        return math.nan
    return v


def z_class(z):
    if math.isnan(z):
        return 'z?'
    for hi, name in ((2, 'z<2'), (5, 'z2-5'), (20, 'z5-20'), (100, 'z20-100')):
        if z < hi:
            return name
    return 'z>100'


def load_summary(run_dir):
    rows = []
    with open(os.path.join(run_dir, 'summary.tsv')) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            r['z_class'] = z_class(fnum(r.get('z')))
            rows.append(r)
    return rows


def median(vals):
    v = sorted(x for x in vals if not math.isnan(x))
    if not v:
        return math.nan
    n = len(v)
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def fmt(v, col):
    if isinstance(v, str):
        return v
    if math.isnan(v):
        return 'nan'
    if col.startswith('cpu'):
        return '%.2gms' % (1000.0 * v)
    if col in ('n', 'n_ok'):
        return '%d' % v
    if abs(v) >= 100:
        return '%.0f' % v
    return '%.3f' % v


def matches(r, filters):
    for k, v in filters:
        if r.get(k) != v:
            return False
    return True


def aggregate(rows, by, cols):
    """{group_key: {objective: {col: value}}}"""
    groups = OrderedDict()
    for r in rows:
        key = tuple(r.get(b, '') for b in by)
        groups.setdefault(key, defaultdict(list))[r['objective']].append(r)
    out = OrderedDict()
    for key, by_obj in groups.items():
        out[key] = {}
        for obj, rs in by_obj.items():
            vals = {}
            n = sum(fnum(r['n']) for r in rs)
            n_ok = sum(fnum(r['n_ok']) for r in rs)
            for c in cols:
                if c == 'fail_rate':
                    vals[c] = (1.0 - n_ok / n) if n > 0 else math.nan
                elif c in ('n', 'n_ok'):
                    vals[c] = n if c == 'n' else n_ok
                else:
                    vals[c] = median(fnum(r.get(c)) for r in rs)
            vals['_nprob'] = len(rs)
            out[key][obj] = vals
    return out


def print_table(agg, by, cols, objectives):
    header = list(by) + ['#'] + ['%s:%s' % (o, SHORT.get(c, c)) for c in cols for o in objectives]
    lines = [header]
    for key, by_obj in agg.items():
        nprob = max((v['_nprob'] for v in by_obj.values()), default=0)
        line = list(key) + [str(nprob)]
        for c in cols:
            for o in objectives:
                line.append(fmt(by_obj[o][c], c) if o in by_obj else '-')
        lines.append(line)
    widths = [max(len(l[i]) for l in lines) for i in range(len(header))]
    for l in lines:
        print('  '.join(s.rjust(w) for s, w in zip(l, widths)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run_dir')
    ap.add_argument('--by', default='', help='comma separated dimensions to aggregate by (default: problem)')
    ap.add_argument('--filter', action='append', default=[], help='dim=value (repeatable)')
    ap.add_argument('--cols', default=','.join(DEFAULT_COLS))
    ap.add_argument('--objectives', default='', help='comma separated; default all, chi2 first')
    ap.add_argument('--paired', action='store_true', help='summarize paired.tsv instead')
    ap.add_argument('--compare', default='', help='a second run dir; objectives get suffixes A/B')
    args = ap.parse_args()

    filters = [tuple(f.split('=', 1)) for f in args.filter]
    cols = [c for c in args.cols.split(',') if c]

    if args.paired:
        counts = defaultdict(lambda: [0, 0])
        with open(os.path.join(args.run_dir, 'paired.tsv')) as f:
            for r in csv.DictReader(f, delimiter='\t'):
                both_ok = r['chi2_status'] == 'ok' and r['obj_status'] == 'ok'
                counts[r['objective']][0 if both_ok else 1] += 1
        print('objective  |diff|>3sd  status-differs')
        for o, (a, b) in sorted(counts.items()):
            print('%-10s %9d %15d' % (o, a, b))
        return

    rows = [r for r in load_summary(args.run_dir) if matches(r, filters)]
    if args.compare:
        rows_b = [r for r in load_summary(args.compare) if matches(r, filters)]
        for r in rows:
            r['objective'] += '/A'
        for r in rows_b:
            r['objective'] += '/B'
        rows += rows_b

    if not rows:
        print('No rows match', file=sys.stderr)
        return

    objectives = [o for o in args.objectives.split(',') if o]
    if not objectives:
        seen = OrderedDict()
        for r in rows:
            seen[r['objective']] = 1
        objectives = sorted(seen, key=lambda o: (not o.startswith('chi2'), o))

    by = [b for b in args.by.split(',') if b] or ['problem', 'peak']
    agg = aggregate(rows, by, cols)
    print_table(agg, by, cols, objectives)


if __name__ == '__main__':
    main()
