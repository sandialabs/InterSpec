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

// GEANT4 validation harness for the spherical source geometries (solid ball,
// hollow shell, and shielded variants).  Three geometries, each at a NEAR and a
// FAR (~50 cm) distance:
//
//   G-A  bare self-attenuating Thorium sphere (Th-232 chain energies)
//   G-B  trace source in a soil sphere with two spherical shield layers
//   G-C  void-center spherical shell trace source + outer shield
//
// For each (geometry x distance) it:
//   * exports a vacuum-world GDML and a GPS volume-source macro per energy
//   * runs the MC `compute()` (precision-targeted) and prints FEP/total +/- 1sigma
//
// Then run GEANT4 with the generated artifacts and compare (z-scores), e.g.:
//   source <your-geant4-install>/bin/geant4.sh
//   ceelo_g4val sphere_GA_near.gdml sphere_GA_near_583keV.mac out.csv --histogram
// Both sides score the SAME full-energy window: this harness pins
// ceelo::kDefaultFepWindowKeV explicitly and ceelo_g4val takes it from the same
// header, so the two cannot drift.  (It used to claim 1.5 keV here while pinning
// nothing, which is how the DESIGN.md G-A/B/C numbers came to be quoted at a
// window the code no longer used.)
//
// Vacuum world (no air gap) matches the MC, which does not transport source->det air.

#include "efficiency/EfficiencyCalculator.h"
#include "materials/Material.h"
#include "physics/FepWindow.h"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

using namespace ceelo;

namespace {

// Representative Th-232 chain gamma energies (keV).
const std::vector<double> kEnergies = {238.6, 583.2, 911.2, 2614.5};

// Precision worth paying for: the comparison's sigma is dominated by the GEANT4
//  side, which at 16M isotropic events gives ~0.2% on FEP at NEAR and only ~1.2%
//  at FAR (Omega/4pi ~ 0.10 vs ~1.5e-3).  Matching 0.2% NEAR puts the combined
//  sigma at ~0.28%, so a 1% discrepancy is z ~ 3.5; pushing CeeLo to 0.1% would
//  cost 4x the runtime to move that to 0.22%.  FAR is a gross-error check, so a
//  looser target is the honest allocation.
SimulationConfig precision_config(double energy_keV, double target = 0.002) {
    SimulationConfig config;
    config.energy_keV = energy_keV;
    // FEP only.  It is the priority metric, and being the smaller efficiency it
    //  is also the binding one - `total` converges first and comes along for
    //  free, so asking for both would only cost wall time.
    config.termination.target_fep_rel_precision = target;
    config.termination.target_total_rel_precision = 0.0;
    config.termination.max_events = 200000000;
    config.termination.min_events = 20000;
    // A cored source gives up the electron-containment fast path (single_material
    //  is false), so every source electron takes the full Moliere walk and these
    //  runs are slower.  The old 60 s cap silently truncated below the target and
    //  the stop reason was never printed, so an under-converged number looked
    //  exactly like a converged one.  Both fixed: bigger budget, and the stop
    //  reason goes in the output.
    config.termination.max_wall_seconds = 900.0;
    config.num_threads = 0;
    config.batch_size = 20000;
    return config;
}

const char* stop_reason_name(StopReason r) {
    switch (r) {
    case StopReason::FepPrecision:   return "fep_prec";
    case StopReason::TotalPrecision: return "tot_prec";
    case StopReason::MaxEvents:      return "MAX_EVENTS";
    case StopReason::WallTime:       return "WALL_TIME";
    default:                         return "other";
    }
}

/// Same six columns ceelo_g4val writes, so the two sides can be compared by one
/// reader (see write_our_csv_batch in benchmark_mc_configs.cpp).
void write_multi_csv(const std::string& filename,
                     const std::vector<double>& energies,
                     const std::vector<EfficiencyResult>& results) {
    std::ofstream f(filename);
    f << "# CeeLo simulation results\n"
      << "energy_keV,fep_efficiency,fep_uncertainty,"
         "total_efficiency,total_uncertainty,num_events\n";
    for (size_t i = 0; i < results.size(); ++i) {
        const auto& r = results[i];
        f << std::fixed << std::setprecision(4) << energies[i] << ","
          << std::scientific << std::setprecision(6)
          << r.full_energy_peak_efficiency << "," << r.fep_uncertainty << ","
          << r.total_efficiency << "," << r.total_uncertainty << ","
          << r.num_events_simulated << "\n";
    }
}

// Build the detector + source for one geometry/distance. `tag` selects geometry.
void configure(EfficiencyCalculator& calc, const Material* nai, const Material* th,
               const Material* soil, const Material* al, const Material* fe,
               const std::string& geom, double center_z) {
    // Pin the window explicitly rather than inheriting the default: ceelo_g4val
    //  takes it from the same header, so pinning here is what guarantees the two
    //  sides score the same peak.
    calc.set_fep_window_keV(kDefaultFepWindowKeV);
    // Source-electron transport is OFF by default and every benchmark config in
    //  bench_configs.h turns it on explicitly - this harness never did.  Above
    //  the 1022 keV pair-production threshold the source/shield secondary
    //  channels (annihilation gammas, bremsstrahlung, electron escape) are worth
    //  several percent of TOTAL efficiency, so without this the sphere family
    //  reported a high-energy total deficit that looked like a physics defect.
    //  See DESIGN.md, config 11 at 3000 keV: -3.05% (z = 17) with it disabled.
    calc.enable_source_electron_transport(true);
    calc.set_detector(nai, CylinderDims{3.81, 7.62});  // 3"x3" NaI
    Eigen::Vector3d c(0.0, 0.0, center_z);
    if (geom == "GA") {
        // Bare self-attenuating thorium sphere, R=2 cm.
        calc.set_spherical_source(c, 2.0);
        calc.set_source_material(th);
    } else if (geom == "GB") {
        // Trace source in soil sphere R=3 cm + Al(0.2) + Fe(0.5) spherical shields.
        calc.set_spherical_source(c, 3.0);
        calc.set_source_material(soil);
        calc.add_source_shield(al, 0.2);
        calc.add_source_shield(fe, 0.5);
    } else if (geom == "GC") {
        // Void-center soil shell [2,3] cm + Fe(0.5) outer shield.
        calc.set_spherical_source(c, 3.0, Eigen::Matrix3d::Identity(), 2.0);
        calc.set_source_material(soil);
        calc.add_source_shield(fe, 0.5);
    } else if (geom == "GD") {
        // G-C's shell with the centre FILLED by an attenuating core.  Iron, not
        //  lead: photoelectric absorption in source material and shields emits no
        //  characteristic K X-ray (DESIGN.md, Known Limitations), so a high-Z core
        //  would bias TOTAL low and make it unjudgeable - the same channel that
        //  puts the G-A thorium sphere 22% low at 238.6 keV.  Iron's K X-rays are
        //  6.4/7.1 keV and reabsorb within microns.  It still shadows hard: tau
        //  over the 4 cm central chord runs 3.3 at 238.6 keV down to 1.2 at 2614,
        //  so this tests whether the core attenuation is RIGHT, not just present.
        calc.set_spherical_source(c, 3.0, Eigen::Matrix3d::Identity(), 2.0);
        calc.set_source_material(soil);
        calc.add_source_core(fe, 2.0);          // fills [0,2] exactly
        calc.add_source_shield(fe, 0.5);
    } else if (geom == "GE") {
        // Export control, numerator: a shell cored with its OWN material.  This
        //  is optically identical to a solid ball but emits only from [2,3], so
        //  eps(GE)/eps(GEsolid) is NOT 1 - and both codes must agree on it.
        //  A GEANT4 ratio of 1 would mean GPS confinement leaked into the core;
        //  a ratio matching a void-centre shell would mean the core is missing
        //  from the GDML.  See the G-E note in the validation write-up.
        calc.set_spherical_source(c, 3.0, Eigen::Matrix3d::Identity(), 2.0);
        calc.set_source_material(soil);
        calc.add_source_core(soil, 2.0);
    } else if (geom == "GEsolid") {
        // Export control, denominator.
        calc.set_spherical_source(c, 3.0);
        calc.set_source_material(soil);
    } else if (geom == "GEadd") {
        // Additivity control: four 0.5 cm cores describe the same scene as G-E's
        //  single 2 cm one, so this ratio IS 1 by construction, in both codes.
        //  It tests the nesting/placement of a multi-layer core stack.
        calc.set_spherical_source(c, 3.0, Eigen::Matrix3d::Identity(), 2.0);
        calc.set_source_material(soil);
        for (int i = 0; i < 4; ++i) calc.add_source_core(soil, 0.5);
    } else if (geom == "GF") {
        // Closed-cavity CYLINDER + iron core.  The only geometry carrying
        //  cyl_inner_half_length, and the one the export used to get wrong: the
        //  source solid was written as a through-bore tube regardless, so a
        //  nested stack came out a pipe.  Identity rotation - the export emits no
        //  <rotation>, so nothing here may be tilted.
        calc.set_cylindrical_source(c, 3.0, 3.0, Eigen::Matrix3d::Identity(),
                                    2.0, 2.0);
        calc.set_source_material(soil);
        calc.add_source_core(fe, 2.0, 2.0);
    } else {  // GG
        // Rectangular shell + per-axis iron core, the third export path.
        calc.set_rectangular_source(c, Eigen::Vector3d(3.0, 3.0, 3.0),
                                    Eigen::Matrix3d::Identity(),
                                    Eigen::Vector3d(2.0, 2.0, 2.0));
        calc.set_source_material(soil);
        calc.add_source_core(fe, 2.0, 2.0, 2.0);
    }
}

} // namespace

int main(int argc, char** argv) {
    const uint64_t g4_events = (argc > 1) ? std::stoull(argv[1]) : 16000000ull;
    // argv[2] overrides the NEAR precision target (FAR gets 2x it).
    const double near_target = (argc > 2) ? std::stod(argv[2]) : 0.002;
    // argv[3] is an optional comma-separated list of "<geom>_<dist>" tags to run
    //  (e.g. "GD_near,GF_near").  The full sweep at a tight FEP target is hours
    //  of CPU, and most of it is geometries a given change does not touch.
    std::vector<std::string> only;
    if (argc > 3) {
        std::string spec = argv[3], tok;
        for (char ch : spec + ",") {
            if (ch == ',') { if (!tok.empty()) only.push_back(tok); tok.clear(); }
            else tok += ch;
        }
    }

    static Material nai  = make_NaI();
    static Material th   = Material("Thorium", 11.72, {{90, 1.0}});  // Th metal
    static Material soil = make_Soil();
    static Material al   = make_Aluminum();
    static Material fe   = make_Iron();

    struct Case { std::string geom; std::string dist; double center_z; };
    const std::vector<Case> cases = {
        {"GA", "near", -3.0},  {"GA", "far", -50.0},
        {"GB", "near", -5.0},  {"GB", "far", -50.0},
        {"GC", "near", -5.0},  {"GC", "far", -50.0},
        // Cored geometries.  NEAR is the quantitative gate (Omega/4pi ~ 0.10, so
        //  16M isotropic G4 events give ~0.2% on FEP); FAR is a gross-error check
        //  only (~1.5e-3, so ~1.2%) - DESIGN.md already records that a
        //  quantitative 50 cm cross-check needs a dedicated high-stats run.
        {"GD", "near", -5.0},  {"GD", "far", -50.0},
        {"GF", "near", -5.0},  {"GF", "far", -50.0},
        {"GG", "near", -5.0},  {"GG", "far", -50.0},
        // Export controls: NEAR only - they test the GDML, and distance adds
        //  nothing to that.
        {"GE",      "near", -5.0},
        {"GEsolid", "near", -5.0},
        {"GEadd",   "near", -5.0},
    };

    printf("# Sphere G4 validation — MC efficiencies (3\"x3\" NaI, vacuum world)\n");
    printf("# G4 macro events per energy: %llu\n", (unsigned long long)g4_events);
    printf("# FEP window: %.3f keV half-width (both sides)\n", kDefaultFepWindowKeV);
    printf("# MC precision target: %.3f%% near, %.3f%% far\n",
           100.0 * near_target, 200.0 * near_target);
    printf("# %-8s %-4s %8s  %14s  %14s  %s\n",
           "geom", "dist", "E[keV]", "FEP+/-1sig", "total+/-1sig", "stop");

    for (const auto& cs : cases) {
        // One GDML per geometry/distance (geometry is energy-independent).
        std::string base = "sphere_" + cs.geom + "_" + cs.dist;
        if (!only.empty()
            && std::find(only.begin(), only.end(), cs.geom + "_" + cs.dist)
                   == only.end())
            continue;
        {
            EfficiencyCalculator calc;
            configure(calc, &nai, &th, &soil, &al, &fe, cs.geom, cs.center_z);
            calc.export_geant4_gdml(base + ".gdml", /*vacuum_world=*/true);
        }

        // One CSV per case, in the same six-column schema ceelo_g4val writes, so
        //  one reader parses both sides and z-scores are mechanical rather than
        //  hand-computed.  (The sphere family has never had this.)
        std::vector<EfficiencyResult> results;
        for (double E : kEnergies) {
            EfficiencyCalculator calc;
            configure(calc, &nai, &th, &soil, &al, &fe, cs.geom, cs.center_z);
            char mac[256];
            std::snprintf(mac, sizeof(mac), "%s_%.0fkeV.mac", base.c_str(), E);
            calc.export_geant4_macro(mac, E, g4_events);

            const double target = (cs.dist == "far") ? 2.0 * near_target
                                                     : near_target;
            auto r = calc.compute(precision_config(E, target));
            printf("  %-8s %-4s %8.1f  %.4e+/-%.1e  %.4e+/-%.1e  %s\n",
                   cs.geom.c_str(), cs.dist.c_str(), E,
                   r.full_energy_peak_efficiency, r.fep_uncertainty,
                   r.total_efficiency, r.total_uncertainty,
                   stop_reason_name(r.stop_reason));
            std::fflush(stdout);
            results.push_back(r);
        }
        write_multi_csv(base + "_multi.csv", kEnergies, results);
    }
    printf("# GDML + .mac + *_multi.csv written to the current directory.\n");
    return 0;
}
