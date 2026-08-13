#!/usr/bin/env python3
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

"""Compare 4M-event MC validation runs (with auto-enabled biasing) against
64M-event GEANT4 reference data for all configs with references."""
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
# CeeLo's own results, committed alongside the G4 references, are the default record.
# Pass `--mc-dir build/examples` to compare a fresh dev run instead.
DEFAULT_MC_DIR = os.path.join(HERE, '..', 'tests', 'data', 'ceelo_reference')
REF_DIR = os.path.join(HERE, '..', 'tests', 'data', 'geant4_reference')

# (config, energy_keV) -> reason; documented in DESIGN.md Known Limitations
SKIP = {
    (7, 100): 'Pb K X-ray-dominated total behind 2mm Pb (tau~6.3, exponentially '
              'sensitive); documented exclusion',
    # cfg 8 59/100/200 keV are NO LONGER skipped (June 26 2026): the prior "low-E
    # Marinelli deficit" was a STALE reference. Regenerating the cfg-8 G4 reference
    # at 32M on the current geometry (see the reference CSV header) brought MC into
    # <=0.7% total / <=1.1% FEP agreement for 100-2614 keV; 59 keV FEP carries a
    # genuine -2.1% near-peak-scatter residual (covered by the FEP tolerance below).
}

# config -> (reference csv, documented max |FEP| %, max |total| %)
CONFIGS = {
    1: ('nai_3x3_10cm_multi.csv', 1.0, 0.4),
    2: ('nai_3x3_al1mm_10cm_multi.csv', 1.0, 1.6),
    3: ('labr3_2x2_al05mm_5cm_multi.csv', 1.5, 1.1),
    5: ('czt_1x1x05cm_5cm_multi.csv', 7.0, 0.8),   # e- escape >1 MeV documented
    6: ('nai_3x3_offaxis45_15cm_multi.csv', 0.9, 0.6),
    7: ('nai_3x3_al1mm_pb2mm_15cm_multi.csv', 1.0, 1.1),
    11: ('nai_3x3_fe05cm_shield_10cm_multi.csv', 0.6, 1.2),
    12: ('nai_3x3_ss304box_cellulose_15cm_multi.csv', 2.1, 1.9),
    8: ('nai_3x3_al05mm_marinelli_water_multi.csv', 2.5, 1.0),  # 32M ref regen Jun26; FEP bound by 59 keV -2.1% near-peak-scatter residual
}


def read_multi(path):
    data = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('energy'):
                continue
            p = line.split(',')
            data[round(float(p[0]))] = (float(p[1]), float(p[2]),
                                        float(p[3]), float(p[4]))
    return data


# --- Cascade true-coincidence summing gate ---------------------------------------
# Surfaces the same engine-vs-GEANT4 cascade-summing observables as the C++ ctest
# gate (tests/test_cascade_summing.cpp) into this dashboard. The G4 reference (bands +
# per-decay areas) and the engine areas are two committed CSVs produced from ONE source
# region ("alcyl"): the reference is tests/data/geant4_reference/cascade_summing_multi.csv;
# the engine areas are our_cascade_multi.csv, written by `cascade_observables` (read from
# the same --mc-dir as the efficiency configs). Unlike the +/-%-tolerance efficiency
# configs, these are wide low-stats RATIO bands, so the gate is band-based.
CASCADE_REF = 'cascade_summing_multi.csv'
CASCADE_MC = 'our_cascade_multi.csv'


def read_cascade(path, has_bands):
    """Parse a cascade CSV keyed by (nuclide, estimator, name).

    Reference rows -> (g4, g4_sig, rlo, rhi); MC rows -> (area, area_unc)."""
    rows = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('nuclide'):
                continue
            p = line.split(',')
            key = (p[0], p[1], p[2])
            if has_bands:  # nuclide,estimator,name,lo,hi,peak_keV,g4,g4_sig,rlo,rhi,...
                rows[key] = (float(p[6]), float(p[7]), float(p[8]), float(p[9]))
            else:          # nuclide,estimator,name,area,area_unc
                rows[key] = (float(p[3]), float(p[4]))
    return rows


def compare_cascade(mc_dir):
    """Gate engine cascade-summing areas vs the baked G4 references. Returns #failures."""
    ref_path = os.path.join(REF_DIR, CASCADE_REF)
    mc_path = os.path.join(mc_dir, CASCADE_MC)
    if not os.path.exists(ref_path):
        return 0  # no cascade reference committed -> nothing to gate
    print('\n=== Cascade summing (engine vs GEANT4, per-decay areas) ===')
    if not os.path.exists(mc_path):
        print(f'  MISSING {mc_path} (run `cascade_observables --out {mc_path}` to gate)')
        return 0
    ref = read_cascade(ref_path, has_bands=True)
    mc = read_cascade(mc_path, has_bands=False)
    print(f"{'nuclide':>7} {'est':>4} {'name':>10} {'ratio':>7} {'z':>6} "
          f"{'band':>11}  flag")
    fails = 0
    for key in sorted(ref):
        if key not in mc:
            continue
        g4, g4_sig, rlo, rhi = ref[key]
        area, area_unc = mc[key]
        ratio = area / g4 if g4 > 0 else float('nan')
        sig = math.sqrt(area_unc**2 + g4_sig**2)
        z = (area - g4) / sig if sig > 0 else 0.0
        flag = ''
        if not (rlo < ratio < rhi):
            flag = ' FAIL'
            fails += 1
        nuc, est, name = key
        print(f'{nuc:>7} {est:>4} {name:>10} {ratio:>7.3f} {z:>6.1f} '
              f'{rlo:>4.2f}-{rhi:<4.2f}  {flag}')
    return fails


def main():
    mc_dir = DEFAULT_MC_DIR
    argv = sys.argv[1:]
    if '--mc-dir' in argv:
        mc_dir = argv[argv.index('--mc-dir') + 1]
    fails = 0
    for cfg, (ref_name, tol_fep, tol_tot) in sorted(CONFIGS.items()):
        mc_path = os.path.join(mc_dir, f'our_{cfg}_multi.csv')
        ref_path = os.path.join(REF_DIR, ref_name)
        if not os.path.exists(mc_path):
            print(f'config {cfg}: MISSING {mc_path}')
            fails += 1
            continue
        mc, ref = read_multi(mc_path), read_multi(ref_path)
        print(f'\n=== Config {cfg} (tolerances: FEP {tol_fep}%, total {tol_tot}%) ===')
        print(f"{'E':>6} {'FEP%':>7} {'z_fep':>6} {'TOT%':>7} {'z_tot':>6}  flag")
        for e in sorted(mc):
            if e not in ref:
                continue
            if (cfg, e) in SKIP:
                print(f'{e:>6}     ---    ---     ---    ---  SKIP ({SKIP[(cfg, e)]})')
                continue
            m, r = mc[e], ref[e]
            d_fep = 100.0 * (m[0] - r[0]) / r[0] if r[0] > 0 else float('nan')
            d_tot = 100.0 * (m[2] - r[2]) / r[2] if r[2] > 0 else float('nan')
            s_fep = math.sqrt(m[1]**2 + r[1]**2)
            s_tot = math.sqrt(m[3]**2 + r[3]**2)
            z_fep = abs(m[0] - r[0]) / s_fep if s_fep > 0 else 0.0
            z_tot = abs(m[2] - r[2]) / s_tot if s_tot > 0 else 0.0
            # tolerance + 3-sigma statistical allowance
            allow_fep = tol_fep + 300.0 * s_fep / r[0] if r[0] > 0 else 1e9
            allow_tot = tol_tot + 300.0 * s_tot / r[2] if r[2] > 0 else 1e9
            flag = ''
            if abs(d_fep) > allow_fep:
                flag += ' FEP_FAIL'
                fails += 1
            if abs(d_tot) > allow_tot:
                flag += ' TOT_FAIL'
                fails += 1
            print(f'{e:>6} {d_fep:>7.2f} {z_fep:>6.1f} {d_tot:>7.2f} {z_tot:>6.1f} {flag}')
    fails += compare_cascade(mc_dir)
    print(f'\n{"PASS" if fails == 0 else f"{fails} FAILURES"}')
    return 0 if fails == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
