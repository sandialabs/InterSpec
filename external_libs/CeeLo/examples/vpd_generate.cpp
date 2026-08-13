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

// vpd_generate: per-detector virtual interaction depth delta(E) fit + export,
// validated across a detector zoo (NaI sizes, HPGe, CZT) for breadth/robustness.
//
// For each detector it:
//   1. fits delta(E) over an (energy x distance) MC grid and writes vpd_<name>.txt
//      (exp-of-log coefficients, consumable by InterSpec's p0/p1/p2 columns);
//   2. validates the delta(E)-corrected efficiency against FRESH direct MC at a
//      distance ladder x off-axis angles, reporting %error and z-scores, and the
//      delta=0 naive-front-face baseline for comparison;
//   3. quantifies model (a) [delta + on-axis solid angle] vs the model (b)
//      angular gap [eps(theta)/eps(0)] that a full angular response would close;
//   4. runs the Phase-4 fit-far->predict-close stress.
//
// Usage: vpd_generate [events_per_point]   (default 400000)
//
// All MC uncertainties + z-scores are reported (project convention). This is pure
// orchestration over EfficiencyCalculator::compute(); no transport code changes.

#include "efficiency/EfficiencyCalculator.h"
#include "io/SolidAngle.h"
#include "io/VirtualDepthFit.h"
#include "materials/Material.h"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace ceelo;

namespace {

constexpr double kPi = 3.14159265358979323846;
uint64_t g_events = 400000;

struct DetSpec {
    std::string name;
    Material (*make)();
    DetectorShape shape;
    std::vector<double> dims; // Cylinder {R, halflen}; Box {hx, hy, len}
    std::vector<double> energies_keV;
};

// FEP efficiency + 1-sigma at an off-axis point source (theta about +x, measured
// from the -z axis toward the face, matching the efficiency-grid study convention).
std::pair<double, double> fep_offaxis(EfficiencyCalculator& calc, double E,
                                      double d, double theta_deg, uint64_t nev) {
    const double th = theta_deg * kPi / 180.0;
    calc.set_point_source(Eigen::Vector3d(d * std::sin(th), 0.0, -d * std::cos(th)));
    const EfficiencyResult r = calc.compute(E, nev, 0);
    return {r.full_energy_peak_efficiency, r.fep_uncertainty};
}

} // namespace

int main(int argc, char** argv) {
    if (argc > 1) g_events = std::strtoull(argv[1], nullptr, 10);

    // Detector zoo. HPGe modeled as a solid Ge cylinder (no coaxial bore -- the
    // disk-face assumption needs a solid face); CZT as a thin box.
    const std::vector<double> nai_E = {122, 344, 662, 1173, 1332};
    const std::vector<double> ge_E = {60, 122, 344, 662, 1173, 1332};
    const std::vector<double> czt_E = {60, 122, 344, 662, 1000, 1332};
    std::vector<DetSpec> zoo = {
        {"NaI_2x2", make_NaI, DetectorShape::Cylinder, {2.54, 2.54}, nai_E},
        {"NaI_3x3", make_NaI, DetectorShape::Cylinder, {3.81, 3.81}, nai_E},
        {"NaI_4x4", make_NaI, DetectorShape::Cylinder, {5.08, 5.08}, nai_E},
        {"HPGe_small", make_HPGe, DetectorShape::Cylinder, {2.5, 2.0}, ge_E},
        {"HPGe_large", make_HPGe, DetectorShape::Cylinder, {3.5, 3.0}, ge_E},
        {"CZT_1x1x0.5", make_CZT, DetectorShape::Box, {0.5, 0.5, 0.5}, czt_E},
    };

    // Fit grid + independent validation grid.
    const std::vector<double> fit_ladder = {1, 2, 3, 5, 8, 12, 20, 30};
    const std::vector<double> val_dist = {1, 2, 5, 10, 25, 50};
    const std::vector<double> val_angles = {0, 15, 30, 45, 60};

    printf("# vpd_generate -- virtual depth delta(E) fit + validation\n");
    printf("# events/point = %llu\n", static_cast<unsigned long long>(g_events));

    for (const DetSpec& spec : zoo) {
        Material mat = spec.make();
        DetectorDescriptor desc =
            make_descriptor(spec.name, spec.shape, mat, spec.dims);

        printf("\n========================================================\n");
        printf("Detector: %s  (%s, R=%.3f cm, R%s, len=%.2f cm)\n",
               spec.name.c_str(), desc.material_name.c_str(),
               desc.crystal_radius_cm,
               spec.shape == DetectorShape::Box ? "_eff(box)" : "(cyl)",
               desc.crystal_length_cm);

        // ---- 1. Fit delta(E) and export ----
        VpdFitConfig cfg;
        cfg.energies_keV = spec.energies_keV;
        cfg.distance_ladder_cm = fit_ladder;
        cfg.events_per_point = g_events;
        cfg.poly_order = 3;
        VpdFit fit = fit_virtual_depth(desc, mat, spec.dims, cfg);

        const std::string fname = "vpd_" + spec.name + ".txt";
        save_vpd(fit, fname);
        printf("  [fit] wrote %s   rms_log_residual=%.4f\n", fname.c_str(),
               fit.rms_log_residual);
        printf("  %8s %8s %12s %10s %10s\n", "E[keV]", "delta", "intrinsic",
               "CV", "maxres%");
        for (const VpdPoint& p : fit.points)
            printf("  %8.0f %8.3f %12.4e %10.4f %10.2f\n", p.energy_keV, p.delta_cm,
                   p.intrinsic_mean, p.cv_residual, 100.0 * p.max_abs_res);

        // ---- 2. Validation vs fresh direct MC (on-axis), + naive delta=0 ----
        // Reference intrinsic for the naive baseline: anchor at d=10 cm, delta=0.
        EfficiencyCalculator calc;
        calc.set_detector(spec.shape, &mat, spec.dims);
        printf("\n  [validate on-axis] eps_pred(delta) vs fresh MC; "
               "naive=delta0 anchored @10cm\n");
        printf("  %8s %6s %12s %12s %9s %8s %7s %9s\n", "E[keV]", "d[cm]",
               "eps_pred", "eps_mc", "mc_unc", "err%", "z", "naive%");
        for (double E : spec.energies_keV) {
            const double delta_E = fit.delta_at(E);
            // Anchor naive model at 10 cm (delta=0): intrinsic0 = eps10 / f(10).
            auto [eps10, u10] = fep_offaxis(calc, E, 10.0, 0.0, g_events);
            (void)u10;
            const double intr0 = eps10 / disk_solid_angle_fraction(10.0, desc.crystal_radius_cm);
            // Recover intrinsic(E) from the fit's per-energy plateau (grid
            // energies == validation energies here).
            double intrinsic_E = 0.0;
            for (const VpdPoint& p : fit.points)
                if (std::abs(p.energy_keV - E) < 1e-6) intrinsic_E = p.intrinsic_mean;
            for (double d : val_dist) {
                auto [eps_mc, unc] = fep_offaxis(calc, E, d, 0.0, g_events);
                // Exported model: eps = intrinsic(E) * f_Omega(d + delta(E)).
                const double eps_pred = intrinsic_E *
                    disk_solid_angle_fraction(d + delta_E, desc.crystal_radius_cm);
                const double naive = intr0 *
                    disk_solid_angle_fraction(d, desc.crystal_radius_cm);
                const double err = (eps_mc > 0) ? 100.0 * (eps_pred - eps_mc) / eps_mc : 0.0;
                const double z = (unc > 0) ? (eps_pred - eps_mc) / unc : 0.0;
                const double nerr = (eps_mc > 0) ? 100.0 * (naive - eps_mc) / eps_mc : 0.0;
                printf("  %8.0f %6.0f %12.4e %12.4e %9.2e %8.2f %7.2f %9.2f\n",
                       E, d, eps_pred, eps_mc, unc, err, z, nerr);
            }
        }

        // ---- 3. Model comparison: off-axis angular gap (model b would close) ----
        printf("\n  [model gap] eps(theta)/eps(0) @10cm "
               "(model-a ignores this; ~ the model-b upside)\n");
        printf("  %8s", "E[keV]");
        for (double a : val_angles) printf("   %4.0fdeg", a);
        printf("\n");
        for (double E : spec.energies_keV) {
            auto [e0, u0] = fep_offaxis(calc, E, 10.0, 0.0, g_events);
            (void)u0;
            printf("  %8.0f", E);
            for (double a : val_angles) {
                auto [e, ue] = fep_offaxis(calc, E, 10.0, a, g_events);
                (void)ue;
                printf("   %6.3f", e0 > 0 ? e / e0 : 0.0);
            }
            printf("\n");
        }

        // ---- 4. Phase-4 stress: fit FAR, predict CLOSE ----
        printf("\n  [stress] fit FAR {8,12,20,30}, predict CLOSE {1,2,3} vs MC\n");
        printf("  %8s %6s %10s %7s\n", "E[keV]", "d[cm]", "err%(far)", "err0%");
        const std::vector<double> far_d = {8, 12, 20, 30};
        const std::vector<double> close_d = {1, 2, 3};
        for (double E : spec.energies_keV) {
            std::vector<double> fe, fu;
            for (double d : far_d) {
                auto [eps, unc] = fep_offaxis(calc, E, d, 0.0, g_events);
                fe.push_back(eps);
                fu.push_back(unc);
            }
            VpdPoint far_fit = fit_delta_one_energy(desc.crystal_radius_cm, far_d, fe, fu);
            for (double d : close_d) {
                auto [eps_mc, unc] = fep_offaxis(calc, E, d, 0.0, g_events);
                (void)unc;
                const double pred = far_fit.intrinsic_mean *
                    disk_solid_angle_fraction(d + far_fit.delta_cm, desc.crystal_radius_cm);
                // delta=0 (naive) far->close for contrast.
                const double intr0 = far_fit.intrinsic_mean; // far plateau already ~intrinsic
                const double naive = intr0 *
                    disk_solid_angle_fraction(d, desc.crystal_radius_cm);
                const double e1 = (eps_mc > 0) ? 100.0 * (pred - eps_mc) / eps_mc : 0.0;
                const double e0 = (eps_mc > 0) ? 100.0 * (naive - eps_mc) / eps_mc : 0.0;
                printf("  %8.0f %6.0f %10.2f %7.2f\n", E, d, e1, e0);
            }
        }
    }

    printf("\n# done\n");
    return 0;
}
