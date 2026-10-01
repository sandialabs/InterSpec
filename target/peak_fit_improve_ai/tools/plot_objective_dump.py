#!/usr/bin/env python3
"""Plots a peak_fit_objective_eval `--dump` file: the channel data, its noise-free expectation, and
each objective's fitted model, with the per-channel residuals below.  Needs numpy + matplotlib.

  plot_objective_dump.py DUMP.tsv [--out PNG] [--log]
"""
import argparse
import csv
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dump')
    ap.add_argument('--out', default='')
    ap.add_argument('--log', action='store_true')
    args = ap.parse_args()

    with open(args.dump) as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    lo = np.array([float(r['energy_lower']) for r in rows])
    hi = np.array([float(r['energy_upper']) for r in rows])
    x = 0.5 * (lo + hi)
    data = np.array([float(r['data']) for r in rows])
    expected = np.array([float(r['expected']) for r in rows])
    models = [k for k in rows[0] if k.startswith('model_')]

    fig, (ax, axr) = plt.subplots(2, 1, figsize=(10, 7), sharex=True,
                                  gridspec_kw={'height_ratios': [3, 1]})
    ax.step(x, data, where='mid', color='0.3', lw=0.8, label='data')
    ax.plot(x, expected, color='0.6', ls=':', lw=1.2, label='noise-free expectation')
    for m in models:
        y = np.array([float(r[m]) for r in rows])
        ax.plot(x, y, lw=1.3, label=m[6:])
        axr.plot(x, (data - y) / np.sqrt(np.maximum(y, 1.0)), lw=1.0, label=m[6:])
    if args.log:
        ax.set_yscale('log')
    ax.set_ylabel('counts / channel')
    ax.legend(fontsize=8)
    ax.set_title(os.path.basename(args.dump), fontsize=9)
    axr.axhline(0, color='0.5', lw=0.6)
    axr.set_ylabel('(data-model)/sqrt(model)')
    axr.set_xlabel('energy (keV)')
    fig.tight_layout()
    out = args.out or os.path.splitext(args.dump)[0] + '.png'
    fig.savefig(out, dpi=110)
    print(out)


if __name__ == '__main__':
    main()
