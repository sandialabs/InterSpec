#!/usr/bin/env python3
"""search_cmp.py A B [--lost N] : compares two `fit_peaks_corpus_eval --search-only` results.

A and B are either run directories, or tags of tools/search_all.sh (every ${tag}_<set>_<dwell> pair
under $FPR_WORK/runs is compared).  Per set: truth photopeaks found by class (fg = the source+background
spectrum, rec = the background spectrum after background-peak recovery), unexplained search peaks (on
no truth photopeak, 511 keV or escape peak), and the totals; then the N truth lines B lost and gained
with the largest truth z.
"""
import csv, os, sys, glob, argparse

def load(run):
    truth, peaks = {}, {}
    with open(os.path.join(run, 'search_truth.tsv')) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            truth[(r['problem'], r['set'], round(float(r['energy']), 2))] = r
    with open(os.path.join(run, 'search_peaks.tsv')) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            peaks.setdefault(r['set'], []).append(r)
    return truth, peaks

def counts(truth, peaks):
    c = {}
    for (prob, s, e), r in truth.items():
        if s not in ('fg', 'bg_recovered'):
            continue
        key = ('fg' if s == 'fg' else 'rec') + '_' + r['class']
        n, f = c.get(key, (0, 0))
        c[key] = (n + 1, f + int(r['found']))
    for s in ('fg', 'recovered'):
        c[s + '_unexplained'] = (sum(1 for p in peaks.get(s, []) if p['verdict'] == 'unexplained'), 0)
        c[s + '_peaks'] = (sum(1 for p in peaks.get(s, []) if p['verdict'] != 'below_min_energy'), 0)
    return c

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('a'); ap.add_argument('b')
    ap.add_argument('--lost', type=int, default=15, help='lost/gained lines to list')
    ap.add_argument('--min-class', default='moderate', help='list lines of this class and above')
    args = ap.parse_args()
    work = os.environ.get('FPR_WORK', os.path.expanduser('~/fit_peaks_work'))
    if os.path.isdir(args.a) and os.path.isdir(args.b):
        pairs = [(os.path.basename(args.b.rstrip('/')), args.a, args.b)]
    else:
        pairs = []
        for run_b in sorted(glob.glob(os.path.join(work, 'runs', args.b + '_*'))):
            if not os.path.isdir(run_b):
                continue
            suffix = os.path.basename(run_b)[len(args.b) + 1:]
            run_a = os.path.join(work, 'runs', args.a + '_' + suffix)
            if os.path.isdir(run_a) and os.path.exists(os.path.join(run_b, 'search_truth.tsv')):
                pairs.append((suffix, run_a, run_b))
    classes = ['strong', 'moderate'] if args.min_class == 'moderate' else ['strong']
    cols = ['fg_strong', 'fg_moderate', 'fg_weak', 'fg_unexplained', 'rec_strong', 'rec_moderate', 'recovered_unexplained']
    print('%-14s' % 'set' + ''.join('%22s' % c for c in cols))
    totals = {}
    lost_all, gained_all = [], []
    for name, run_a, run_b in pairs:
        ta, pa = load(run_a); tb, pb = load(run_b)
        ca, cb = counts(ta, pa), counts(tb, pb)
        line = '%-14s' % name
        for c in cols:
            na, fa = ca.get(c, (0, 0)); nb, fb = cb.get(c, (0, 0))
            if c.endswith('unexplained'):
                line += '%22s' % ('%d -> %d' % (na, nb))
                totals[c] = (totals.get(c, (0, 0, 0, 0))[0] + na, totals.get(c, (0, 0, 0, 0))[1] + nb, 0, 0)
            else:
                line += '%22s' % ('%d -> %d /%d' % (fa, fb, nb))
                t = totals.get(c, (0, 0, 0, 0))
                totals[c] = (t[0] + fa, t[1] + fb, t[2] + nb, 0)
        print(line)
        for key, r in tb.items():
            if key[1] not in ('fg', 'bg_recovered') or r['class'] not in classes or key not in ta:
                continue
            fa, fb = int(ta[key]['found']), int(r['found'])
            row = (float(r['z']), name, key[0], key[1], float(r['energy']), r['class'],
                   ta[key]['found_marg_z'], ta[key]['found_det_z'], r['found_marg_z'], r['found_det_z'])
            if fa and not fb:
                lost_all.append(row)
            elif fb and not fa:
                gained_all.append(row)
    if len(pairs) > 1:
        line = '%-14s' % 'TOTAL'
        for c in cols:
            t = totals.get(c, (0, 0, 0, 0))
            line += '%22s' % (('%d -> %d' % (t[0], t[1])) if c.endswith('unexplained') else ('%d -> %d /%d' % (t[0], t[1], t[2])))
        print(line)
    for label, rows in (('LOST', lost_all), ('GAINED', gained_all)):
        print('\n%s (%d, %s and above):' % (label, len(rows), args.min_class))
        for row in sorted(rows, reverse=True)[:args.lost]:
            print('  z=%6.1f %-12s %-22s %-12s %9.2f keV %-8s  A: marg %s det %s   B: marg %s det %s' % row)

if __name__ == '__main__':
    main()
