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

"""Guard the fitter/evaluator lockstep on the source-secondary PROCESS features.

Companion to `check_basis_lockstep.py`, which guards the B-spline basis
normalization. This one guards the other coordinate the two sides must agree
on: the 7-D process feature vector.

`fit_source_secondaries.py` does not compute those features -- it reads them
from `ceelo_source_cross_section_features --emit-process-features`, which calls
`detail::source_process_features`, the same function the evaluator uses. That
makes lockstep structural. This script proves the property still holds end to
end, because "structural" only means "hard to break", not "impossible": a
future change could give the helper and the evaluator different cross-section
blends, different tau conventions, or a different axis order.

Method -- a base-free ratio test. Build two surrogates that are identical
except that the two nodes' log conditional yields are SWAPPED, then evaluate
both at the first node's own cell. The mechanistic base cancels exactly in the
ratio, so what remains depends only on the kernel weights, hence only on the
distance between the evaluator's query feature and the stored node features:

    ratio = exp(w_a * L_a + w_b * L_b) / exp(w_a * L_b + w_b * L_a)
          = exp((w_a - w_b) * (L_a - L_b))

with w normalized. If the evaluator queried at a different feature vector than
the helper emitted, the weights change and the ratio moves. A single node would
NOT detect this: its normalized weight is 1 regardless of where it sits.

Usage: check_feature_lockstep.py [--build BUILD_DIR]  (default build-rel under
the CeeLo root); exit 0 = pass.
"""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
# Set from --build in main().
HELPER = REPO / "build-rel/tools/full_spectrum/ceelo_source_cross_section_features"
EVALUATOR = REPO / "build-rel/tools/full_spectrum/ceelo_source_model_evaluate"

BASIS_SIZE = 32
FEATURE_SIZE = 7
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
# The payload component label the evaluator emits for each fitted component.
OUTPUT_LABEL = {
    "fluorescence_k_alpha_atom": "fluorescence_k_alpha",
    "fluorescence_k_beta_atom": "fluorescence_k_beta",
    "fluorescence_l_atom": "fluorescence_l",
    "pair_511_atom": "pair_511",
    "pair_annihilation_continuum": "pair_annihilation_continuum",
    "pair_annihilation_in_flight": "pair_annihilation_in_flight",
    "photoelectron_brems": "photoelectron_brems",
    "compton_electron_brems": "compton_electron_brems",
    "pair_electron_brems": "pair_electron_brems",
}

# Two well-separated cells, both above the pair threshold so every component
# has a nonzero mechanistic base (cell A, lead, also above its K edge with L
# lines; cell B, iron, has no fluorescence above the 10 keV cut, so its
# fluorescence bases are zero and only cell A exercises them).
CELL_A = (3000.0, 82, 5.0, 0.0)
CELL_B = (1500.0, 26, 20.0, 0.0)

# Bandwidths are derived per run from the measured node separation (see
# calibrate_bandwidth) rather than hard-coded. With two far-apart nodes the
# normalized kernel weight saturates at 1 and the ratio stops responding to
# small feature shifts, which would make this guard look sensitive while
# actually being nearly blind.
LOG_YIELD_A = 0.0
LOG_YIELD_B = math.log(1000.0)
# Target squared kernel distance between the two probe nodes. Order 1 keeps
# both weights comparable, which is where the ratio responds most steeply.
TARGET_DISTANCE2 = 1.0


def helper_features(cells):
    """The 7-D process features, straight from the C++ helper."""
    text = "".join(
        f"{energy} {z} {rho} {f_h}\n" for energy, z, rho, f_h in cells)
    completed = subprocess.run(
        [str(HELPER), "--emit-process-features"], input=text, text=True,
        stdout=subprocess.PIPE, check=True)
    rows = list(csv.DictReader(completed.stdout.splitlines()))
    answer: list[dict[str, np.ndarray]] = []
    for index in range(len(cells)):
        block = rows[index * len(COMPONENTS):(index + 1) * len(COMPONENTS)]
        answer.append({
            row["component"]: np.asarray(
                [float(row[f"f{axis}"]) for axis in range(FEATURE_SIZE)])
            for row in block})
    return answer


def calibrate_bandwidth(features) -> tuple[np.ndarray, np.ndarray]:
    """Bandwidths that put the two probe nodes ~1 kernel distance apart.

    Averaged over components so every component shares one bandwidth vector,
    as a real surrogate does.
    """
    deltas = np.asarray([
        np.abs(features[0][name] - features[1][name])
        for name in COMPONENTS])
    delta = np.maximum(deltas.mean(axis=0), 1e-6)
    # sum((delta / bandwidth)^2) = TARGET_DISTANCE2 when every axis
    # contributes equally.
    bandwidth = delta * math.sqrt(FEATURE_SIZE / TARGET_DISTANCE2)
    return bandwidth, bandwidth * 1.3


def expected_ratio(delta: np.ndarray, bandwidth: np.ndarray) -> float:
    distance2 = float(np.sum((delta / bandwidth) ** 2))
    weight_b = math.exp(-0.5 * distance2)
    normalized = (1.0 - weight_b) / (1.0 + weight_b)
    return math.exp(normalized * (LOG_YIELD_A - LOG_YIELD_B))


def write_surrogate(
    path: Path, features, swap: bool,
    yield_bandwidth: np.ndarray, shape_bandwidth: np.ndarray,
) -> None:
    shape = np.full(BASIS_SIZE, 1.0 / BASIS_SIZE)
    with path.open("w") as stream:
        stream.write("CEELO_SOURCE_SECONDARY_SURROGATE 1\n")
        stream.write("inverse_energy_warp_scale 2\n")
        stream.write("absorption_envelope 0 3.5 10\n")
        stream.write("basis cubic_bspline\n")
        stream.write("local_linear 0 2\n")
        stream.write("tail_decay 0\n")
        for name in COMPONENTS:
            stream.write(
                f"{name} {BASIS_SIZE - 1} "
                + " ".join(f"{v:.17g}" for v in yield_bandwidth) + " "
                + " ".join(f"{v:.17g}" for v in shape_bandwidth) + " 2\n")
            yields = ([LOG_YIELD_B, LOG_YIELD_A] if swap
                      else [LOG_YIELD_A, LOG_YIELD_B])
            for cell_index, log_yield in enumerate(yields):
                stream.write(
                    f"node {name} "
                    + " ".join(
                        f"{v:.17g}"
                        for v in features[cell_index][name])
                    + f" {log_yield:.17g} "
                    + " ".join(f"{v:.17g}" for v in shape) + "\n")
        stream.write("END\n")


def evaluate_yields(surrogate: Path, cell) -> dict[str, float]:
    energy, z, rho, f_h = cell
    # No --payload here: the evaluator treats --payload and the loose
    # surrogate files as mutually exclusive (payload wins), so passing both
    # would silently ignore the probe surrogate this check depends on.
    completed = subprocess.run(
        [str(EVALUATOR),
         "--secondary-surrogate", str(surrogate),
         "--energy-kev", str(energy), "--z", str(z),
         "--areal-density", str(rho),
         "--hydrogen-mass-fraction", str(f_h),
         "--bins", "64"],
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
            print(f"missing: {path}", file=sys.stderr)
            return 2

    features = helper_features([CELL_A, CELL_B])
    yield_bandwidth, shape_bandwidth = calibrate_bandwidth(features)

    with tempfile.TemporaryDirectory() as directory:
        plain = Path(directory) / "plain.txt"
        swapped = Path(directory) / "swapped.txt"
        write_surrogate(
            plain, features, False, yield_bandwidth, shape_bandwidth)
        write_surrogate(
            swapped, features, True, yield_bandwidth, shape_bandwidth)
        plain_yields = evaluate_yields(plain, CELL_A)
        swapped_yields = evaluate_yields(swapped, CELL_A)

    print("component                        expected      measured    "
          "rel.err   verdict")
    worst = 0.0
    failures = 0
    for name in COMPONENTS:
        delta = features[0][name] - features[1][name]
        expected = expected_ratio(delta, yield_bandwidth)

        label = OUTPUT_LABEL[name]
        numerator = plain_yields.get(label, 0.0)
        denominator = swapped_yields.get(label, 0.0)
        if denominator <= 0.0 or numerator <= 0.0:
            print(f"{name:32s} {'--':>12s} {'--':>13s}   "
                  f"{'component absent':>9s}  SKIP")
            continue
        measured = numerator / denominator
        error = abs(measured - expected) / expected
        worst = max(worst, error)
        verdict = "OK" if error < 1e-6 else "MISMATCH"
        failures += verdict == "MISMATCH"
        print(f"{name:32s} {expected:12.6g} {measured:13.6g} "
              f"{error:9.2e}  {verdict}")

    print()
    if failures:
        print(f"{failures} component(s) disagree: the evaluator is querying "
              "the process-feature field at a different coordinate than the "
              "helper emits, so the fitted nodes no longer mean what the "
              "fitter thought they meant.")
        return 1
    print(f"process features are in lockstep (worst relative error "
          f"{worst:.2e}).")
    # Report the sensitivity actually achieved, computed rather than asserted.
    probe = features[0][COMPONENTS[-1]] - features[1][COMPONENTS[-1]]
    shifted = probe.copy()
    shifted[3] += 0.01
    reference = expected_ratio(probe, yield_bandwidth)
    moved = expected_ratio(shifted, yield_bandwidth)
    print(
        "sensitivity: a 0.01 shift on one feature axis moves the ratio by "
        f"{abs(moved - reference) / reference:.2%}, against a {worst:.0e} "
        "noise floor.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
