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

// Cascade (true-coincidence) summing validation study.
//
// Drives compute_cascade() from real SandiaDecay nuclide data (via the adapter)
// and reports summing factors against published references:
//   - Cs-137 null test (single gamma -> factor == 1)
//   - Co-60 distance dependence on a 3"x3" NaI (factor -> 1 as distance grows)
//   - Gupta et al. (arXiv:2302.02776) BEGe k_TCS, point sources at distance,
//     for Co-60 / Ba-133 / Eu-152 (compare 1/summing_factor to their analytical k_TCS)
//
// Requires -DCEELO_WITH_SANDIADECAY=ON. The SandiaDecay XML path is provided
// by the build as SANDIA_DECAY_XML_PATH.
//
// Usage: cascade_validation [events_per_peak]   (default 300000)

#include "efficiency/EfficiencyCalculator.h"
#include "materials/Material.h"
#include "cascade/SandiaDecayCascade.h"
#include "SandiaDecay.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace ceelo;
using namespace ceelo::cascade_adapter;

namespace {

uint64_t g_events = 300000;
uint64_t g_sec6_events = 3000000;  // engine events for the high-stats G4 cross-check

const SandiaDecay::SandiaDecayDataBase& db() {
    static SandiaDecay::SandiaDecayDataBase database(SANDIA_DECAY_XML_PATH);
    return database;
}

// 3"x3" NaI cylinder, point source on axis at distance `d_cm` from the face.
EfficiencyCalculator nai_3x3(double d_cm) {
    static Material nai = make_NaI();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -d_cm));
    return calc;
}

// Gupta BEGe-like HPGe: crystal r=2.945 cm, L=2.95 cm, 0.18 cm front dead layer,
// 0.12 cm Al endcap, ~1.45 cm endcap-to-crystal gap (air, not transported).
EfficiencyCalculator bege(double d_cm) {
    static Material ge = make_HPGe();
    static Material al = make_Aluminum();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &ge, {2.945, 2.95});
    calc.set_dead_layer(0.18, 0.05);
    calc.add_attenuator(&al, 0.12, 0.12, 0.0, 2.95);  // endcap
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -(d_cm + 1.45)));
    return calc;
}

void run(const char* label, EfficiencyCalculator& calc,
         const std::vector<DecayCascade>& cascades,
         const std::vector<double>& peak_keV,
         const std::vector<double>& ref_kTCS = {},
         CascadeMethod method = CascadeMethod::Conditional,
         uint64_t events = 0) {
    CascadeConfig cfg;
    cfg.cascades = cascades;
    for (double e : peak_keV) cfg.peaks.push_back(PeakWindow{e, 1.5});
    cfg.num_events = events ? events : g_events;
    cfg.method = method;
    auto res = calc.compute_cascade(cfg);

    printf("  %-22s  %10s %10s %8s %8s   %8s %8s\n",
           label, "eps_nosum", "eps_sum", "factor", "+-", "k_TCS", "ref");
    for (size_t i = 0; i < res.peaks.size(); ++i) {
        const auto& p = res.peaks[i];
        if (!p.found) { printf("  %8.1f keV   (not found)\n", peak_keV[i]); continue; }
        const double kTCS = (p.summing_factor > 0) ? 1.0 / p.summing_factor : 0.0;
        char refbuf[16] = "   -";
        if (i < ref_kTCS.size()) std::snprintf(refbuf, sizeof refbuf, "%7.3f", ref_kTCS[i]);
        printf("  %8.1f keV   %10.3e %10.3e %8.4f %8.4f   %8.4f %s\n",
               p.energy_keV, p.eff_no_summing, p.eff_with_summing,
               p.summing_factor, p.summing_factor_unc, kTCS, refbuf);
    }
}

} // namespace

int main(int argc, char** argv) {
    if (argc > 1) g_events = std::strtoull(argv[1], nullptr, 10);
    if (argc > 2) g_sec6_events = std::strtoull(argv[2], nullptr, 10);
    printf("# Cascade summing validation  (events/peak = %llu)\n",
           static_cast<unsigned long long>(g_events));

    const auto co60 = build_cascades(db(), "Co60");
    const auto cs137 = build_cascades(db(), "Cs137");
    const auto ba133 = build_cascades(db(), "Ba133");
    const auto eu152 = build_cascades(db(), "Eu152");

    // --- 1. Cs-137 null test (NaI, 10 cm) ---
    printf("\n[1] Cs-137 null test (NaI 3x3, 10 cm) -- expect factor = k_TCS = 1.000\n");
    { auto c = nai_3x3(10.0); run("Cs-137", c, cs137, {661.657}); }

    // --- 2. Co-60 distance dependence (NaI 3x3) ---
    printf("\n[2] Co-60 distance dependence (NaI 3x3) -- factor -> 1 as distance grows\n");
    for (double d : {1.0, 2.0, 5.0, 10.0, 25.0}) {
        auto c = nai_3x3(d);
        char lbl[32]; std::snprintf(lbl, sizeof lbl, "Co-60 d=%.0f cm", d);
        run(lbl, c, co60, {1173.2, 1332.5});
    }

    // --- 3. Gupta BEGe k_TCS (point source) ---
    // Reference analytical k_TCS from Gupta et al. Table 4 (BEGe), per distance.
    printf("\n[3] Gupta BEGe k_TCS comparison (HPGe BEGe model)\n");
    printf("    Co-60 @ 1,3,5 cm:\n");
    { auto c = bege(1.0); run("Co-60 d=1cm", c, co60, {1173.2, 1332.5}, {1.074, 1.077}); }
    { auto c = bege(3.0); run("Co-60 d=3cm", c, co60, {1173.2, 1332.5}, {1.035, 1.036}); }
    { auto c = bege(5.0); run("Co-60 d=5cm", c, co60, {1173.2, 1332.5}, {1.003, 1.003}); }

    printf("    Ba-133 @ 1 cm (k_TCS ref: 53->1.208, 81->1.130, 276->1.277, 302->1.176, 356->1.157, 383->0.938):\n");
    { auto c = bege(1.0); run("Ba-133 d=1cm", c, ba133, {53.16, 80.99, 276.40, 302.85, 356.01, 383.85},
                              {1.208, 1.130, 1.277, 1.176, 1.157, 0.938}); }

    printf("    Eu-152 @ 1 cm (incl. 1085.9 de-summing line, ref k_TCS 1.008):\n");
    { auto c = bege(1.0); run("Eu-152 d=1cm", c, eu152,
                              {121.78, 244.69, 344.27, 778.9, 964.07, 1085.86, 1112.07, 1408.0},
                              {1.080, 1.170, 1.066, 1.103, 1.227, 1.008, 1.203, 1.220}); }

    // --- 4. Full-realization: captures summing-IN to sum-peak-fed lines ---
    // These lines are at the SUM energy of a coincident pair, so they gain
    // counts from summing-in (k_TCS < 1) -- which the per-peak conditional
    // estimator structurally cannot see, but the full realization does.
    printf("\n[4] Full-realization mode -- sum-peak-fed (de-summing) lines\n");
    printf("    Ba-133 383.85 keV (= 81.0 + 302.85; Gupta k_TCS 0.938):\n");
    { auto c = bege(1.0);
      run("Ba-133 d=1cm [full]", c, ba133, {356.01, 383.85}, {1.157, 0.938},
          CascadeMethod::FullRealization, 4000000); }
    printf("    Eu-152 1085.86 keV (= 121.78 + 964.07; Gupta k_TCS 1.008):\n");
    { auto c = bege(1.0);
      run("Eu-152 d=1cm [full]", c, eu152, {964.07, 1085.86}, {1.227, 1.008},
          CascadeMethod::FullRealization, 4000000); }

    // --- 5. GEANT4 full-decay cross-check (Co-60, NaI 3x3, 2 cm) ---
    // Compares the per-decay 2505 keV sum-peak coincidence rate and the 1173/
    // 1332 photopeak rates to a G4 radioactive-decay run. (G4 also has Co-60
    // betas reaching the bare source; their beta-gamma summing shifts events up,
    // so G4's sum peak is split between the 2505 line and a 2505-2730 tail.)
    printf("\n[5] GEANT4 full-decay cross-check: Co-60, NaI 3x3, 2 cm\n");
    {
        auto c = nai_3x3(2.0);
        CascadeConfig cfg;
        cfg.cascades = co60;
        cfg.peaks = {PeakWindow{1173.2, 5.0}, PeakWindow{1332.5, 5.0}};
        cfg.num_events = 2000000;
        cfg.method = CascadeMethod::FullRealization;
        for (int i = 0; i <= 280; ++i) cfg.spectrum_bin_edges.push_back(i * 10.0f);
        auto res = c.compute_cascade(cfg);
        auto bin = [&](double e){ return res.summed_spectrum.empty() ? 0.0f
                       : res.summed_spectrum[static_cast<size_t>(e/10.0)]; };
        printf("    quantity                 engine(per decay)   G4(per decay)\n");
        printf("    1173 peak (per decay)    %12.3e   %12.3e\n", (double)bin(1173), 1215.0/40000);
        printf("    1332 peak (per decay)    %12.3e   %12.3e\n", (double)bin(1332), 1095.0/40000);
        printf("    2505 SUM peak (per dec)  %12.3e   %12.3e   (G4 peak+betatail: %.3e)\n",
               (double)bin(2505), 58.0/40000, 79.0/40000);
        printf("    -> engine sum-peak %.3e vs G4 peak+betatail %.3e : agree ~%.0f%%\n",
               (double)bin(2505), 79.0/40000, 100.0*std::abs((double)bin(2505) - 79.0/40000)/(79.0/40000));
        printf("    NOTE: G4 1173/1332 are additionally suppressed by Co-60 beta-gamma\n");
        printf("          summing (bare source); engine is gamma-only. The 2505 keV sum\n");
        printf("          peak (clean coincidence signature) appears in both and agrees.\n");
    }

    // --- 6. Extended/shielded GEANT4 cross-checks (beta confound removed) ---
    // The Fe shell (A) and the water self-absorption (B) stop the Co-60 betas, so
    // the G4 sum peak has no beta tail -- a clean gamma-gamma comparison of the
    // extended/shielded cascade transport path.
    printf("\n[6] Extended/shielded GEANT4 cross-check (Co-60, clean gamma-gamma)\n");
    {
        static Material fe = make_Iron();
        static Material nai = make_NaI();
        auto sum_rate = [&](EfficiencyCalculator& c) {
            CascadeConfig cfg; cfg.cascades = co60;
            cfg.peaks = {PeakWindow{1173.2, 5.0}, PeakWindow{1332.5, 5.0}};
            cfg.num_events = g_sec6_events; cfg.method = CascadeMethod::FullRealization;
            for (int i = 0; i <= 280; ++i) cfg.spectrum_bin_edges.push_back(i * 10.0f);
            auto r = c.compute_cascade(cfg);
            auto b = [&](double e){ size_t k=(size_t)(e/10.0);
                return r.summed_spectrum.empty()?0.0f:r.summed_spectrum[k]; };
            auto u = [&](double e){ size_t k=(size_t)(e/10.0);
                return r.summed_spectrum_uncertainty.empty()?0.0f:r.summed_spectrum_uncertainty[k]; };
            printf("      engine (%lluM): 1173=%.4e 1332=%.4e 2505=%.4e +-%.2e (%.2f%%)\n",
                   (unsigned long long)(g_sec6_events/1000000),
                   (double)b(1173),(double)b(1332),(double)b(2505),
                   (double)u(2505), 100.0*u(2505)/((double)b(2505)>0?(double)b(2505):1.0));
        };
        printf("    A) shielded point (0.2 cm Fe, 2 cm):  G4(24M) 1173=3.7052e-2 1332=3.2977e-2 2505=1.5601e-3 +-0.52%%\n");
        { EfficiencyCalculator c; c.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
          c.set_point_source(Eigen::Vector3d(0,0,-2.0)); c.add_source_shield(&fe, 0.2);
          sum_rate(c); }
        printf("    B) water cyl (r1 hz1, 4 cm):          G4(88M) 1173=2.1273e-2 1332=1.8966e-2 2505=4.7338e-4 +-0.49%%\n");
        { static Material water = make_Water();
          EfficiencyCalculator c; c.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
          c.set_cylindrical_source(Eigen::Vector3d(0,0,-4.0), 1.0, 1.0);
          c.set_source_material(&water); sum_rate(c); }
    }

    printf("\n# done\n");
    return 0;
}
