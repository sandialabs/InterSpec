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

"""The detector bake's bank placement and arrival scoring helpers."""

import io
import sys
import unittest
from pathlib import Path

# The bake tools import their siblings (bank_io) as top-level modules.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import detector_bank_energies as placement  # noqa: E402
import score_arrival_bank as arrival  # noqa: E402

HPGE_COAX_LANDMARK_KEV = 55.000005689471955
LOWE = [35.0, 37.5, 47.5, 50.0, 57.5, 60.0, 65.0, 72.5]


class LandmarkScaling(unittest.TestCase):
    def test_the_reference_detector_keeps_its_energies_exactly(self):
        energies = LOWE + [41.0, 56.3, 120.0, 661.657]
        self.assertEqual(
            placement.landmark_scaled(energies, HPGE_COAX_LANDMARK_KEV),
            energies)

    def test_low_energies_follow_the_landmark(self):
        scaled = placement.landmark_scaled([35.0, 65.0, 99.0, 100.0, 340.0],
                                           63.44)
        self.assertAlmostEqual(scaled[0], 40.35)
        self.assertAlmostEqual(scaled[1], 74.95)
        self.assertAlmostEqual(scaled[2], 114.2)
        self.assertEqual(scaled[3:], [100.0, 340.0])
        for energy in scaled:
            self.assertAlmostEqual(energy * 20.0, round(energy * 20.0))


class HeldOutPlacement(unittest.TestCase):
    def test_energies_on_the_training_grid_move_off_it(self):
        training = [40.0, 59.54, 70.0, 80.0, 122.0, 1000.0]
        placed, moves = placement.held_out([41.0, 69.5, 120.0, 1000.0, 500.0],
                                           training)
        self.assertEqual(placed[0], 41.0)
        self.assertEqual(placed[4], 500.0)
        self.assertEqual([before for before, _ in moves], [69.5, 120.0, 1000.0])
        for energy in placed:
            for train in training:
                self.assertGreaterEqual(placement.separation(energy, train),
                                        placement.MIN_SEPARATION)
        # The smallest move that clears, upward first: 1000 -> 1020; with 1023
        # also on the grid it goes down instead, to 980.
        self.assertEqual(placed[3], 1020.0)
        placed, _ = placement.held_out([1000.0], training + [1023.0])
        self.assertEqual(placed, [980.0])

    def test_a_crowded_grid_is_an_error(self):
        training = [100.0 * 1.01 ** k for k in range(-30, 31)]
        with self.assertRaises(ValueError):
            placement.held_out([100.0], training)


def arrival_cell(energy: float, fep: float, zero: float, continuum: list,
                 digits: int) -> io.StringIO:
    """An arrival CSV of `events` = 1000 with a uniform continuum grid."""
    bins = len(continuum)
    header = ("cell_kind,energy_keV,Z,tau,areal_density_g_cm2,cos_theta,"
              "detector,component,bin_lo_keV,bin_hi_keV,sum_w,sum_w2,events,"
              "seed,analytic_direct,incidence,mu_lo,mu_hi,"
              "mean_projected_area_cm2,mean_projected_area_sigma_cm2\n")
    tail = ",1000,1,,iso,-1,1,80,0.1\n"
    rows = [header,
            f"detector_arrival,{energy},,,,,d,zero_deposit,0,0,{zero},{zero}{tail}",
            f"detector_arrival,{energy},,,,,d,full_energy,{energy},{energy},"
            f"{fep},{fep}{tail}"]
    for b, value in enumerate(continuum):
        lo = f"{energy * b / bins:.{digits}g}"
        hi = f"{energy * (b + 1) / bins:.{digits}g}"
        rows.append(f"detector_arrival,{energy},,,,,d,single_compton,{lo},{hi},"
                    f"{value},{value}{tail}")
    return io.StringIO("".join(rows))


class ArrivalScoring(unittest.TestCase):
    def test_bins_match_whatever_the_printed_precision(self):
        energy = 47.3
        continuum = [float(10 + (b % 7)) for b in range(64)]
        truth = arrival.read_cell(arrival_cell(energy, 300, 400, continuum, 12))
        model = arrival.read_cell(arrival_cell(energy, 300, 400, continuum, 17))
        bins = len(truth["bins"])
        self.assertEqual(bins, 64)
        for with_fep in (True, False):
            a = arrival.spectrum(truth, energy, bins, with_fep)
            b = arrival.spectrum(model, energy, bins, with_fep)
            self.assertEqual(arrival.total_variation(a, b), 0.0)
        self.assertEqual(arrival.valley(truth, energy, bins),
                         arrival.valley(model, energy, bins))

    def test_the_valley_window_sits_between_edge_and_line(self):
        lo, hi = arrival.valley_window(661.657)
        edge = 661.657 * (2 * 661.657 / 510.99895) / (1 + 2 * 661.657 / 510.99895)
        self.assertAlmostEqual(lo, edge + 0.25 * (661.657 - edge))
        self.assertAlmostEqual(hi, 661.657 * 0.98)

    def test_a_fixed_angle_cell_is_refused(self):
        text = arrival_cell(100.0, 1, 1, [1.0] * 8, 12).getvalue()
        text = text.replace("detector_arrival,", "detector,")
        with self.assertRaises(ValueError):
            arrival.read_cell(io.StringIO(text))


if __name__ == "__main__":
    unittest.main()
