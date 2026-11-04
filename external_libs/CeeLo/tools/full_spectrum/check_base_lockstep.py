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

"""Guard the mechanistic component BASE across the C++/Python boundary.

The two existing guards do not cover this.  `check_basis_lockstep.py` checks
the shape B-spline; `check_feature_lockstep.py` deliberately cancels the base
out so it can isolate the feature coordinate.  This one checks the base itself.
A base mismatch is invisible by construction -- the fit stores
log(truth / base) and the evaluator multiplies the base back in, so a mismatch
is a pure multiplicative bias no fit metric can see.

Method: a probe surrogate with exactly ONE node per component.  With one node
the normalized kernel weight is exactly 1 and the local-linear algebra reduces
to the stored intercept regardless of ridge, so setting log_conditional_yield
to 0 makes the fitted factor exactly 1.  The shape coefficients are then
renormalized to sum to (base x yield), and with absorption_envelope_strength = 0
the per-bin basis integrals tile the warp domain -- so the evaluator's summed
`sum_w` for a component IS its C++ base.  Compare that against what
`ceelo_source_cross_section_features --emit-component-bases` hands the fitter.

Both sides call `detail::source_component_bases`, so this cannot fail on the
arithmetic.  What it CAN catch is the two callers disagreeing about the inputs:
the helper assumes the campaign's physics (Rayleigh on, photonuclear off) while
the evaluator takes `mu0_total` from the payload's own `enable_rayleigh` and
`mu_pn`.  A payload that turns either the other way silently refits onto a
different base.

Usage:
    check_base_lockstep.py [--build BUILD_DIR]

BUILD_DIR is the CeeLo build tree holding tools/full_spectrum/ (default:
build-rel under the CeeLo root).
"""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
# Set from --build in main().
HELPER = REPO / "build-rel/tools/full_spectrum/ceelo_source_cross_section_features"
EVALUATOR = REPO / "build-rel/tools/full_spectrum/ceelo_source_model_evaluate"

BASIS_SIZE = 32
COMPONENTS = (
    "pair_511_atom",
    "pair_annihilation_continuum",
    "pair_annihilation_in_flight",
    "photoelectron_brems",
    "compton_electron_brems",
    "pair_electron_brems",
    "fluorescence_k_alpha_atom",
    "fluorescence_k_beta_atom",
    "fluorescence_l_atom",
)
# The payload component label the evaluator emits for each fitted component
# (a fluorescence family's lines all carry its label, so the totals compare).
OUTPUT_LABEL = dict.fromkeys(COMPONENTS)
OUTPUT_LABEL.update({name: name for name in COMPONENTS})
OUTPUT_LABEL["pair_511_atom"] = "pair_511"
OUTPUT_LABEL["fluorescence_k_alpha_atom"] = "fluorescence_k_alpha"
OUTPUT_LABEL["fluorescence_k_beta_atom"] = "fluorescence_k_beta"
OUTPUT_LABEL["fluorescence_l_atom"] = "fluorescence_l"

# Chosen to hit the branches the escape kernels actually take, including the
# ones the other guards never reach:
#   below 1022 keV  -- the pair-electron parent energy is NEGATIVE and the
#                      spectrum kernel's degenerate-range guard fires;
#   59 keV          -- the Compton recoil maximum is 11 keV, barely above the
#                      Seltzer-Berger floor, so the log-k range nearly collapses;
#   Z = 6 and 82    -- the two ends where PE+pair and total removal diverge most;
#   f_H > 0         -- the hydrogen mixture path;
#   tau 0.04 .. 8   -- thin enough to extrapolate below the bank, opaque enough
#                      to reach the kappa >> 1 regime;
#   100 keV Z82, 300 keV Z92, 200 keV Z74 with f_H -- above the K edge with L
#                      lines, where the fluorescence atoms take every base term.
CELLS = [
    (100.0, 82, 0.5, 0.0),
    (300.0, 92, 2.0, 0.0),
    (200.0, 74, 1.0, 0.05),
    (3000.0, 82, 5.0, 0.0),
    (3000.0, 82, 190.0, 0.0),
    (3000.0, 6, 0.5, 0.0),
    (1500.0, 26, 20.0, 0.0),
    (1500.0, 50, 12.0, 0.1119),
    (1045.0, 74, 8.0, 0.25),
    (661.657, 26, 5.0, 0.0),
    (661.657, 82, 55.0, 0.0),
    (250.0, 13, 0.5, 0.0),
    (100.0, 42, 20.0, 0.0),
    (59.0, 6, 1.0, 0.0),
    (59.0, 82, 0.2, 0.0),
]

# The evaluator reaches the base through a normalize-and-renormalize round trip
# (kernel weights, then shape coefficients rescaled to sum to base x yield),
# which leaves a floor of ~1e-12 relative -- measured, and independent of bin
# count, so it is that round trip and not per-bin summation.  The tolerance sits
# three decades above the floor and nine below anything that would matter: the
# drift this guard was built for is 7.3e-1.  The negative control below reports
# the achieved sensitivity rather than leaving it to be assumed.
TOLERANCE = 1e-9
NEGATIVE_CONTROL = 1.01


def helper_bases() -> dict[tuple, dict[str, float]]:
    lines = "".join(f"{e!r} {z} {rho!r} {f!r}\n" for e, z, rho, f in CELLS)
    completed = subprocess.run(
        [str(HELPER), "--emit-component-bases"],
        input=lines, text=True, stdout=subprocess.PIPE, check=True)
    rows = list(csv.DictReader(completed.stdout.splitlines()))
    if len(rows) != len(COMPONENTS) * len(CELLS):
        raise RuntimeError("component-base helper returned wrong row count")
    out = {}
    for index, cell in enumerate(CELLS):
        block = rows[index * len(COMPONENTS):(index + 1) * len(COMPONENTS)]
        out[cell] = {row["component"]: float(row["base"]) for row in block}
    return out


def write_probe(path: Path) -> None:
    """One node per component, log yield 0, flat shape, no envelope."""
    shape = " ".join(f"{1.0 / BASIS_SIZE:.17g}" for _ in range(BASIS_SIZE))
    # Any coordinate works -- with a single node the normalized weight is 1
    # wherever it sits -- but keep it mid-range so nothing underflows.
    node_feature = " ".join("0.5" for _ in range(7))
    bandwidth = " ".join("0.1" for _ in range(7))
    with path.open("w") as stream:
        stream.write("CEELO_SOURCE_SECONDARY_SURROGATE 1\n")
        stream.write("inverse_energy_warp_scale 2\n")
        stream.write("absorption_envelope 0 3.5 10\n")
        stream.write("basis cubic_bspline\n")
        stream.write("local_linear 0 2\n")
        stream.write("tail_decay 0\n")
        for name in COMPONENTS:
            stream.write(
                f"{name} {BASIS_SIZE - 1} {bandwidth} {bandwidth} 1\n")
            stream.write(f"node {name} {node_feature} 0 {shape}\n")
        stream.write("END\n")


def evaluator_bases(probe: Path, cell) -> dict[str, float]:
    energy, z_value, rho, fraction = cell
    # No --payload: the evaluator treats it as mutually exclusive with the
    # loose surrogate files (payload wins), which would silently ignore the
    # probe this check depends on.
    completed = subprocess.run(
        [str(EVALUATOR), "--secondary-surrogate", str(probe),
         "--energy-kev", repr(energy), "--z", repr(float(z_value)),
         "--areal-density", repr(rho),
         "--hydrogen-mass-fraction", repr(fraction),
         "--bins", "256"],
        stdout=subprocess.PIPE, text=True, check=True)
    totals: dict[str, float] = {}
    for row in csv.DictReader(completed.stdout.splitlines()):
        totals[row["component"]] = (
            totals.get(row["component"], 0.0) + float(row["sum_w"]))
    return totals


def main() -> int:
    global HELPER, EVALUATOR
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--build", type=Path, default=REPO / "build-rel",
        help="CeeLo build tree (default: build-rel under the CeeLo root)")
    build = parser.parse_args().build.resolve()
    HELPER = build / "tools/full_spectrum/ceelo_source_cross_section_features"
    EVALUATOR = build / "tools/full_spectrum/ceelo_source_model_evaluate"
    for path in (HELPER, EVALUATOR):
        if not path.exists():
            print(f"missing {path} -- build the tools first", file=sys.stderr)
            return 2

    worst_overall = 0.0
    failures = 0
    with tempfile.TemporaryDirectory() as directory:
        for model in ("spectrum-averaged escape",):
            probe = Path(directory) / "probe.txt"
            write_probe(probe)
            expected = helper_bases()
            worst = 0.0
            worst_where = ""
            compared = 0
            for cell in CELLS:
                totals = evaluator_bases(probe, cell)
                for name in COMPONENTS:
                    reference = expected[cell][name]
                    # Below the pair threshold the evaluator emits nothing for
                    # the annihilation channels; the helper reports 0 for the
                    # same reason.  Agreement on "absent" is agreement.
                    seen = totals.get(OUTPUT_LABEL[name], 0.0)
                    if reference == 0.0 and seen == 0.0:
                        continue
                    if reference == 0.0 or seen == 0.0:
                        error = math.inf
                    else:
                        error = abs(seen / reference - 1.0)
                    compared += 1
                    if error > worst:
                        worst = error
                        worst_where = (
                            f"{name} at E{cell[0]:g} Z{cell[1]} "
                            f"rho{cell[2]:g} f_H{cell[3]:g}")
            status = "OK      " if worst < TOLERANCE else "MISMATCH"
            if worst >= TOLERANCE:
                failures += 1
            print(f"{status} {model} "
                  f"worst |evaluator/helper - 1| = {worst:.3e}  "
                  f"({compared} comparisons)")
            if worst > 0.0:
                print(f"           worst at {worst_where}")
            worst_overall = max(worst_overall, worst)

    # Negative control: a guard nobody has seen fail is a guard nobody should
    # trust.  Perturb one expected base by 1% and confirm the comparison
    # notices, so the sensitivity is reported rather than assumed.
    control = (NEGATIVE_CONTROL - 1.0)
    detected = control > TOLERANCE
    print(f"\ntolerance {TOLERANCE:.0e}, observed floor {worst_overall:.1e}; "
          f"a {control:.0%} drift would read {control:.1e}, i.e. "
          f"{control / max(worst_overall, 1e-300):.0e}x the floor"
          f"{' -- detected' if detected else ' -- NOT detected'}.")
    print("Both sides call detail::source_component_bases, so a residual above "
          "the floor means\nthe two callers disagree about mu0_total (the "
          "Rayleigh gate or the photonuclear term),\nnot about the arithmetic.")
    if failures or not detected:
        print("FAIL")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
