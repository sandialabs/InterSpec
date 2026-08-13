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

// Gamma-gamma angular-correlation (W(theta)) demonstration for compute_cascade.
//
// Builds Co-60 cascades from the enriched SandiaDecay data twice -- once with the
// angular correlation enabled (a2/a4 carried onto the coincidence links) and once
// forced isotropic -- and compares the resulting true-coincidence-summing factors
// on a close 3"x3" NaI geometry, where the effect is largest. Reports the 1173 /
// 1332 keV summing-out factors and the 2505 keV sum peak (FullRealization), with
// 1-sigma uncertainties and the z-score of the W-on vs W-off difference.
//
// Co-60 a2 ~ 0.10 > 0 aligns the two gammas, modestly changing the joint
// detection probability vs isotropic. This isolates the *feature's* effect; an
// absolute validation is against GEANT4 run with /process/had/deex/correlatedGamma
// true (harness flag --correlated-gamma).
//
// Usage: cascade_wtheta_demo [events] [distance_cm]

#include "efficiency/EfficiencyCalculator.h"
#include "materials/Material.h"
#include "cascade/SandiaDecayCascade.h"
#include "SandiaDecay.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace ceelo;
using namespace ceelo::cascade_adapter;

namespace {
const SandiaDecay::SandiaDecayDataBase& db() {
    static SandiaDecay::SandiaDecayDataBase database(SANDIA_DECAY_XML_PATH);
    return database;
}
}

int main(int argc, char** argv) {
    const uint64_t events = (argc > 1) ? std::strtoull(argv[1], nullptr, 10) : 4000000;
    const double d_cm = (argc > 2) ? std::atof(argv[2]) : 2.0;

    CascadeOptions opts_corr;   // angular_correlations = true (default)
    CascadeOptions opts_iso;
    opts_iso.angular_correlations = false;

    const auto casc_corr = build_cascades(db(), "Co60", opts_corr);
    const auto casc_iso  = build_cascades(db(), "Co60", opts_iso);

    Material nai = make_NaI();
    auto make_calc = [&]() {
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -d_cm));
        return calc;
    };

    std::printf("Co-60 W(theta) demo: 3\"x3\" NaI, point source %.1f cm on-axis, %llu events\n",
                d_cm, (unsigned long long)events);
    std::printf("%-10s %-10s  %12s %12s   %7s\n",
                "peak", "method", "factor(W)", "factor(iso)", "z");

    auto run_pair = [&](double peak, CascadeMethod method, const char* mname) {
        auto run_one = [&](const std::vector<DecayCascade>& c) {
            EfficiencyCalculator calc = make_calc();
            CascadeConfig cfg;
            cfg.cascades = c;
            cfg.peaks.push_back(PeakWindow{peak, 1.5});
            cfg.num_events = events;
            cfg.method = method;
            return calc.compute_cascade(cfg).peaks.at(0);
        };
        const auto w = run_one(casc_corr);
        const auto i = run_one(casc_iso);
        const double dz = w.summing_factor - i.summing_factor;
        const double sz = std::sqrt(w.summing_factor_unc * w.summing_factor_unc +
                                    i.summing_factor_unc * i.summing_factor_unc);
        const double z = (sz > 0) ? dz / sz : 0.0;
        std::printf("%-10.1f %-10s  %8.5f+-%.5f %8.5f+-%.5f  %6.2f\n",
                    peak, mname, w.summing_factor, w.summing_factor_unc,
                    i.summing_factor, i.summing_factor_unc, z);
    };

    run_pair(1173.2, CascadeMethod::Conditional, "Cond");
    run_pair(1332.5, CascadeMethod::Conditional, "Cond");
    // The 2505 keV sum peak only exists with summing-in (FullRealization).
    run_pair(2505.7, CascadeMethod::FullRealization, "FullReal");
    return 0;
}
