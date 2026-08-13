/* CeeLo: a Monte Carlo photon-transport library for computing gamma-ray
 and X-ray detector efficiency - developed as part of InterSpec.

 Copyright 2026 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov.

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

/// @file cascade_correction.cpp
/// @brief Demonstrates true coincidence summing (TCS) cascade corrections
///        for Co-60 in a 3"×3" NaI detector at various source-face distances.
///
/// Co-60 emits two cascade gammas (1173.2 keV and 1332.5 keV) in prompt
/// coincidence.  At close geometry, both gammas may deposit energy in the
/// same detector resolving time, shifting or removing events from the
/// individual photopeaks.
///
/// Reference:
///   Debertin & Helmer, "Gamma and X-Ray Spectrometry with Semiconductor
///   Detectors" (1988), §6.3; PLAN.md §8.
///
/// Expected results (PLAN.md §8.7):
///   At ~1 cm:  large correction (~30–40%), C_out ≈ 0.60–0.70
///   At ~25 cm: small correction (< 2%),    C_out > 0.98

#include "efficiency/EfficiencyCalculator.h"
#include "materials/Material.h"

#include <iostream>
#include <iomanip>
#include <map>
#include <vector>

int main() {
    using namespace ceelo;

    Material nai = make_NaI();

    // Source distances to evaluate.
    const std::vector<double> distances = {1.0, 5.0, 10.0, 25.0, 50.0};
    const uint64_t N_per_energy = 50'000;

    std::cout << "\n  Co-60 cascade correction for 3\"×3\" NaI\n"
              << "  " << N_per_energy << " events per gamma per distance\n\n"
              << std::left
              << std::setw(12) << "Distance"
              << std::setw(16) << "ε_total(1173)"
              << std::setw(16) << "ε_total(1332)"
              << std::setw(16) << "C_out(1173)"
              << std::setw(16) << "C_out(1332)"
              << "\n"
              << std::string(76, '-') << "\n";

    for (double d : distances) {
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -d));

        auto res1173 = calc.compute(1173.2, N_per_energy, 0);
        auto res1332 = calc.compute(1332.5, N_per_energy, 0);

        std::map<double, EfficiencyResult> cache;
        cache[1173.2] = res1173;
        cache[1332.5] = res1332;

        // For 1173 keV primary: coincident with 1332.5 keV (f = 1.0)
        auto corr1173 = EfficiencyCalculator::cascade_correction(
            res1173.full_energy_peak_efficiency,
            res1173.total_efficiency,
            {{1332.5, 1.0}}, cache);

        // For 1332 keV primary: coincident with 1173.2 keV (f = 1.0)
        auto corr1332 = EfficiencyCalculator::cascade_correction(
            res1332.full_energy_peak_efficiency,
            res1332.total_efficiency,
            {{1173.2, 1.0}}, cache);

        std::cout << std::left
                  << std::setw(6) << d << " cm    "
                  << std::fixed << std::setprecision(4)
                  << std::setw(16) << res1173.total_efficiency
                  << std::setw(16) << res1332.total_efficiency
                  << std::setw(16) << corr1173.summing_out_factor
                  << std::setw(16) << corr1332.summing_out_factor
                  << "\n";
    }

    std::cout << "\nExpected C_out(1173) at 1 cm:  ~0.60–0.70  (literature: ~30–40% correction)\n"
              << "Expected C_out(1173) at 25 cm: ~0.98–1.00  (literature: <2% correction)\n";

    return 0;
}
