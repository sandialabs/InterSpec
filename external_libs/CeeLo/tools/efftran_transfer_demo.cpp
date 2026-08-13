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

// efftran_transfer_demo -- EFFTRAN-style efficiency transfer, validated against
// direct Monte Carlo (mirrors EFFTRAN's et.bat workflow).
//
// Runs CeeLo MC ONCE at a reference distance on-axis to build an anchor
// efficiency curve, then transfers it -- via the interaction-weighted effective
// solid angle ratio K(E,pos)/K(E,ref) (io/EfficiencyTransfer.h) -- to a ladder
// of distances and off-axis angles. At every target it also runs a DIRECT MC
// and prints the transfer-vs-MC comparison (eps, MC uncertainty, z-score,
// fractional error), quantifying the ~1% transfer envelope. It also emits a
// make_transfer_response DetectorResponse and cross-checks it against the
// standalone helper.
//
// Usage:
//   efftran_transfer_demo [preset] [precision] [seed]
//     preset     nai3x3 | czt        (default nai3x3)
//     precision  per-point MC FEP target   (default 0.01)
//     seed       base MC seed              (default 1)

#include "efficiency/EfficiencyCalculator.h"
#include "geometry/Geometry.h"
#include "io/DetectorResponse.h"
#include "io/EfficiencyTransfer.h"
#include "io/ResponseGenerator.h"
#include "io/ResponseKernel.h"
#include "materials/Material.h"

#include <Eigen/Core>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

using namespace ceelo;

namespace {

GeometryDescriptor preset_descriptor(const std::string& name,
                                     std::vector<double>& energies_keV) {
    GeometryDescriptor gd;
    if (name == "czt") {
        // Bare 1x1x0.5 cm CZT box.
        gd.shape = DetectorShape::Box;
        gd.dimensions_cm = {0.5, 0.5, 0.5};
        gd.materials = {MaterialSpec::from(make_CZT())};
        gd.crystal_material_index = 0;
        energies_keV = {59.5, 122.0, 356.0, 662.0};
    } else {
        // NaI 3"x3" + 0.5 mm Al can.
        gd.shape = DetectorShape::Cylinder;
        gd.dimensions_cm = {3.81, 7.62};
        gd.materials = {MaterialSpec::from(make_NaI()),
                        MaterialSpec::from(make_Aluminum())};
        gd.crystal_material_index = 0;
        LayerSpec can;
        can.material_index = 1;
        can.front_thickness_cm = 0.05;
        can.side_thickness_cm = 0.05;
        can.z_end_cm = 7.62;
        gd.layers.push_back(can);
        energies_keV = {59.5, 122.0, 662.0, 1332.0};
    }
    return gd;
}

EfficiencyResult run_mc(EfficiencyCalculator& calc, const Eigen::Vector3d& src,
                        double energy_keV, double precision, uint64_t seed) {
    calc.set_point_source(src);
    SimulationConfig cfg;
    cfg.energy_keV = energy_keV;
    cfg.termination.target_fep_rel_precision = precision;
    cfg.termination.target_total_rel_precision = precision;
    cfg.termination.max_events = 8000000;
    cfg.termination.min_events = 40000;
    cfg.seed = seed;
    return calc.compute(cfg);
}

}  // namespace

int main(int argc, char** argv) {
    const std::string preset = argc > 1 ? argv[1] : "nai3x3";
    const double precision = argc > 2 ? std::atof(argv[2]) : 0.01;
    const uint64_t seed = argc > 3 ? std::strtoull(argv[3], nullptr, 10) : 1;

    std::vector<double> energies;
    const GeometryDescriptor gd = preset_descriptor(preset, energies);
    const double a = gd.transverse_half_extent();
    const double d_ref = 10.0 * a;  // far-field on-axis reference

    // Configure the MC engine and build the ray-trace geometry from the SAME
    // descriptor (owned-material vectors must outlive their consumers).
    std::vector<std::unique_ptr<Material>> owned_calc, owned_geom;
    EfficiencyCalculator calc;
    ResponseGenerator::configure_calculator(calc, gd, owned_calc);
    const Geometry geom = gd.build_geometry(owned_geom);

    std::printf("# EFFTRAN transfer demo  preset=%s  a=%.3f cm  d_ref=%.2f cm "
                "(%.0f a)  precision=%.3f\n",
                preset.c_str(), a, d_ref, d_ref / a, precision);

    // ---- 1. MC anchor curve at the reference position ----------------------
    AnchorCurve fep_anchor, tot_anchor;
    std::printf("# Anchor (MC @ %.2f cm on-axis):\n", d_ref);
    std::printf("#   E_keV      eps_fep   frac_sig    eps_tot   frac_sig\n");
    for (size_t i = 0; i < energies.size(); ++i) {
        const EfficiencyResult r =
            run_mc(calc, {0.0, 0.0, -d_ref}, energies[i], precision, seed + i);
        fep_anchor.energies_keV.push_back(energies[i]);
        fep_anchor.eff.push_back(r.full_energy_peak_efficiency);
        fep_anchor.frac_sigma.push_back(r.fep_uncertainty /
                                        std::max(r.full_energy_peak_efficiency, 1e-30));
        tot_anchor.energies_keV.push_back(energies[i]);
        tot_anchor.eff.push_back(r.total_efficiency);
        tot_anchor.frac_sigma.push_back(r.total_uncertainty /
                                        std::max(r.total_efficiency, 1e-30));
        std::printf("#  %7.1f  %.4e  %.5f   %.4e  %.5f\n", energies[i],
                    r.full_energy_peak_efficiency, fep_anchor.frac_sigma.back(),
                    r.total_efficiency, tot_anchor.frac_sigma.back());
    }

    EfficiencyTransfer et_fep(geom, {0.0, 0.0, -d_ref}, fep_anchor,
                              MuChoice::Total);
    EfficiencyTransfer et_tot(geom, {0.0, 0.0, -d_ref}, tot_anchor,
                              MuChoice::NoRayleigh);

    // ---- 2. transfer vs direct MC over a distance x angle ladder -----------
    const std::vector<double> d_over_a{1.0, 1.5, 2.0, 3.0, 5.0, 20.0};
    const std::vector<double> cos_thetas{1.0, 0.5};
    std::printf("\n# quantity=FEP\n");
    std::printf("# %7s %8s %6s %8s  %11s %11s %10s %7s %8s\n", "E_keV", "d_cm",
                "d/a", "cos_th", "eps_xfer", "eps_MC", "mc_unc", "z", "err%");
    double worst_onaxis_far = 0.0, worst_offaxis = 0.0;
    uint64_t sd = seed + 1000;
    for (double dr : d_over_a) {
        const double d = dr * a;
        for (double ct : cos_thetas) {
            const Eigen::Vector3d tgt = source_position(d, ct);
            for (size_t i = 0; i < energies.size(); ++i) {
                const double E = energies[i];
                const double eps_x = et_fep.eps_at(E, tgt);
                const EfficiencyResult r = run_mc(calc, tgt, E, precision, sd++);
                const double mc = r.full_energy_peak_efficiency;
                const double unc = r.fep_uncertainty;
                const double z = (eps_x - mc) / std::max(unc, 1e-30);
                const double err = 100.0 * (eps_x - mc) / std::max(mc, 1e-30);
                std::printf("  %7.1f %8.2f %6.2f %8.2f  %.4e %.4e %.3e %7.2f %8.2f\n",
                            E, d, dr, ct, eps_x, mc, unc, z, err);
                if (ct > 0.99 && dr >= 2.0)
                    worst_onaxis_far = std::max(worst_onaxis_far, std::fabs(err));
                if (ct < 0.99)
                    worst_offaxis = std::max(worst_offaxis, std::fabs(err));
            }
        }
    }
    std::printf("# worst |err| on-axis (d>=2a): %.2f%%   worst |err| off-axis: "
                "%.2f%%\n", worst_onaxis_far, worst_offaxis);

    // ---- 3. cross-check make_transfer_response vs the standalone helper -----
    auto resp = make_transfer_response(gd, fep_anchor, {0.0, 0.0, -d_ref},
                                       &tot_anchor);
    const Eigen::Vector3d chk = source_position(3.0 * a, 1.0);
    std::printf("\n# make_transfer_response vs standalone (on-axis, d=3a):\n");
    std::printf("#   E_keV   resp.eps_fep   standalone   rel_diff%%\n");
    for (double E : energies) {
        const EffResult rr = resp->eps_fep(E, 0.0, 0.0, 3.0 * a);
        const double sa = et_fep.eps_at(E, chk);
        std::printf("#  %7.1f   %.4e   %.4e   %8.3f\n", E, rr.value, sa,
                    100.0 * (rr.value - sa) / std::max(sa, 1e-30));
    }

    return 0;
}
