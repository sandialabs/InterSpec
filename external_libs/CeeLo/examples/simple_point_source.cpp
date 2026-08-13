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

/// @file simple_point_source.cpp
/// @brief Demonstrates computing ε_FEP and ε_total for a point source.
///
/// Simulates Cs-137 (661.7 keV) and the two Co-60 lines (1173.2 and 1332.5 keV)
/// for a 3"×3" NaI detector with an on-axis point source at z = -10 cm.
///
/// Printed output includes a comparison to published values from:
///   Knoll (2010), "Radiation Detection and Measurement" 4th ed., Fig. 10-10.

#include "efficiency/EfficiencyCalculator.h"
#include "materials/Material.h"

#include <iostream>
#include <iomanip>
#include <vector>

int main() {
    using namespace ceelo;

    // 3"×3" NaI(Tl): radius = 3.81 cm (1.5"), length = 7.62 cm (3")
    Material nai = make_NaI();

    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));  // 10 cm from front face

    // Simulate 100k events per energy using all available hardware threads.
    const uint64_t N = 100'000;

    struct TestCase {
        const char* isotope;
        double      energy_keV;
        double      ref_fep;   // literature ε_FEP (Knoll Fig. 10-10)
    };

    // Reference: Knoll (2010), 3"×3" NaI, on-axis, source-face distance = 10 cm.
    std::vector<TestCase> cases = {
        {"Cs-137",  661.7, 0.013},   // ≈ 1.2–1.8%
        {"Co-60",  1173.2, 0.007},   // ≈ 0.6–1.0%
        {"Co-60",  1332.5, 0.006},   // ≈ 0.5–0.8%
    };

    std::cout << "\n  3\"×3\" NaI, on-axis point source at z = -10 cm\n"
              << "  " << N << " events per energy\n\n"
              << std::left
              << std::setw(10) << "Isotope"
              << std::setw(12) << "E (keV)"
              << std::setw(16) << "ε_FEP"
              << std::setw(16) << "ε_total"
              << std::setw(14) << "Lit. ε_FEP"
              << "\n"
              << std::string(68, '-') << "\n";

    for (const auto& tc : cases) {
        auto res = calc.compute(tc.energy_keV, N, 0 /* auto threads */);

        std::cout << std::left
                  << std::setw(10) << tc.isotope
                  << std::setw(12) << std::fixed << std::setprecision(1) << tc.energy_keV
                  << std::setw(16) << std::scientific << std::setprecision(3)
                  << res.full_energy_peak_efficiency
                  << "± " << std::setw(12) << res.fep_uncertainty
                  << std::setw(16) << res.total_efficiency
                  << "~" << tc.ref_fep
                  << "\n";
    }

    // Export GDML and macro for GEANT4 validation (optional).
    // calc.export_geant4_gdml("nai_3x3_point.gdml");
    // calc.export_geant4_macro("nai_3x3_point_662.mac", 661.7, N);

    std::cout << "\nDone.\n";
    return 0;
}
