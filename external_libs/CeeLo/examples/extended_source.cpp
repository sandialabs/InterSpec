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

/// @file extended_source.cpp
/// @brief Demonstrates efficiency calculation for extended cylindrical and
///        rectangular sources.
///
/// Compares a point source against a thin coaxial cylindrical source
/// (converges to point source limit) and a larger cylindrical source.
/// Also demonstrates the rectangular source API.

#include "efficiency/EfficiencyCalculator.h"
#include "materials/Material.h"

#include <iostream>
#include <iomanip>

int main() {
    using namespace ceelo;

    Material nai = make_NaI();
    const uint64_t N = 50'000;

    std::cout << "\n  3\"×3\" NaI, various source configurations at 662 keV\n"
              << "  " << N << " events\n\n";

    // -----------------------------------------------------------------------
    // 1. Point source reference at z = -10 cm
    // -----------------------------------------------------------------------
    EfficiencyCalculator calc_pt;
    calc_pt.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_pt.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
    auto res_pt = calc_pt.compute(662.0, N, 0);

    // -----------------------------------------------------------------------
    // 2. Tiny cylindrical source (should converge to point)
    // -----------------------------------------------------------------------
    EfficiencyCalculator calc_tiny;
    calc_tiny.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_tiny.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -10.0), 0.001, 0.001);
    auto res_tiny = calc_tiny.compute(662.0, N, 0);

    // -----------------------------------------------------------------------
    // 3. Larger coaxial cylindrical source (R=2 cm, h=2 cm, centered at -10 cm)
    // -----------------------------------------------------------------------
    EfficiencyCalculator calc_cyl;
    calc_cyl.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_cyl.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -10.0), 2.0, 2.0);
    auto res_cyl = calc_cyl.compute(662.0, N, 0);

    // -----------------------------------------------------------------------
    // 4. Rectangular (box) source (4×4×4 cm, centered at -10 cm)
    // -----------------------------------------------------------------------
    EfficiencyCalculator calc_rect;
    calc_rect.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc_rect.set_rectangular_source(Eigen::Vector3d(0.0, 0.0, -10.0),
                                     Eigen::Vector3d(2.0, 2.0, 2.0));
    auto res_rect = calc_rect.compute(662.0, N, 0);

    auto print_row = [](const char* label, const EfficiencyResult& r) {
        std::cout << std::left << std::setw(32) << label
                  << std::fixed << std::setprecision(5)
                  << "ε_FEP=" << r.full_energy_peak_efficiency
                  << "  ε_total=" << r.total_efficiency
                  << "\n";
    };

    print_row("Point source (z=-10 cm)", res_pt);
    print_row("Tiny cyl source (R=0.001 cm)", res_tiny);
    print_row("Cyl source (R=2, h=2, z=-10)", res_cyl);
    print_row("Box source (2×2×2 cm, z=-10)", res_rect);

    // The tiny cylinder should agree with the point source within ~15%.
    double ratio = res_tiny.total_efficiency / res_pt.total_efficiency;
    std::cout << "\nTiny cylinder / point source ratio: " << ratio
              << " (expect 1.0 ± 0.15)\n";

    // Extended source has lower efficiency since some emission points are farther away.
    std::cout << "Extended cyl / point source ratio: "
              << res_cyl.total_efficiency / res_pt.total_efficiency
              << " (expect < 1.0 since volume spans z = [-12, -8] cm)\n";

    return 0;
}
