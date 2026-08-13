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
// FEP_G4 = (counts in [E-1.5,E+1.5]) / N ; total_G4 = (sum all bins) / N.
//
// Vacuum world (no air gap) matches the MC, which does not transport source->det air.

#include "efficiency/EfficiencyCalculator.h"
#include "materials/Material.h"

#include <cstdio>
#include <string>
#include <vector>

using namespace ceelo;

namespace {

// Representative Th-232 chain gamma energies (keV).
const std::vector<double> kEnergies = {238.6, 583.2, 911.2, 2614.5};

SimulationConfig precision_config(double energy_keV, double target = 0.005) {
    SimulationConfig config;
    config.energy_keV = energy_keV;
    config.termination.target_fep_rel_precision = target;
    config.termination.target_total_rel_precision = target;
    config.termination.max_events = 200000000;
    config.termination.min_events = 20000;
    config.termination.max_wall_seconds = 60.0;
    config.num_threads = 0;
    config.batch_size = 20000;
    return config;
}

// Build the detector + source for one geometry/distance. `tag` selects geometry.
void configure(EfficiencyCalculator& calc, const Material* nai, const Material* th,
               const Material* soil, const Material* al, const Material* fe,
               const std::string& geom, double center_z) {
    calc.set_detector(DetectorShape::Cylinder, nai, {3.81, 7.62});  // 3"x3" NaI
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
    } else {  // GC
        // Void-center soil shell [2,3] cm + Fe(0.5) outer shield.
        calc.set_spherical_source(c, 3.0, Eigen::Matrix3d::Identity(), 2.0);
        calc.set_source_material(soil);
        calc.add_source_shield(fe, 0.5);
    }
}

} // namespace

int main(int argc, char** argv) {
    const uint64_t g4_events = (argc > 1) ? std::stoull(argv[1]) : 16000000ull;

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
    };

    printf("# Sphere G4 validation — MC efficiencies (3\"x3\" NaI, vacuum world)\n");
    printf("# G4 macro events per energy: %llu\n", (unsigned long long)g4_events);
    printf("# %-4s %-4s %8s  %14s  %14s\n",
           "geom", "dist", "E[keV]", "FEP+/-1sig", "total+/-1sig");

    for (const auto& cs : cases) {
        // One GDML per geometry/distance (geometry is energy-independent).
        std::string base = "sphere_" + cs.geom + "_" + cs.dist;
        {
            EfficiencyCalculator calc;
            configure(calc, &nai, &th, &soil, &al, &fe, cs.geom, cs.center_z);
            calc.export_geant4_gdml(base + ".gdml", /*vacuum_world=*/true);
        }
        for (double E : kEnergies) {
            EfficiencyCalculator calc;
            configure(calc, &nai, &th, &soil, &al, &fe, cs.geom, cs.center_z);
            char mac[256];
            std::snprintf(mac, sizeof(mac), "%s_%.0fkeV.mac", base.c_str(), E);
            calc.export_geant4_macro(mac, E, g4_events);

            auto r = calc.compute(precision_config(E));
            printf("  %-4s %-4s %8.1f  %.4e+/-%.1e  %.4e+/-%.1e\n",
                   cs.geom.c_str(), cs.dist.c_str(), E,
                   r.full_energy_peak_efficiency, r.fep_uncertainty,
                   r.total_efficiency, r.total_uncertainty);
            std::fflush(stdout);
        }
    }
    printf("# GDML + .mac files written to the current directory.\n");
    return 0;
}
