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

/// @file generate_g4_benchmarks.cpp
/// @brief Generate GDML + macro files for all §9.3 GEANT4 benchmark configurations,
///        and print our own MC efficiency results in a comparable CSV format.
///
/// After running this program, feed each config into the GEANT4 harness:
///   ./ceelo_g4val detector_<N>.gdml run_<N>_<E>keV.mac g4_<N>_<E>keV.csv
/// Then compare with:
///   python tools/geant4_validation/compare_results.py our_<N>_<E>keV.csv g4_<N>_<E>keV.csv

#include "efficiency/EfficiencyCalculator.h"
#include "materials/Material.h"
#include "export/Geant4Export.h"

#include <Eigen/Core>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

using namespace ceelo;

static const uint64_t N_MC      = 500'000;  // events per energy for our MC
static const uint64_t N_G4      = 500'000;  // events per energy written to macro
static const unsigned N_THREADS = 0;        // 0 = auto

/// Write a one-row CSV in compare_results.py format.
static void write_our_csv(const std::string& filename,
                          double energy_keV,
                          const EfficiencyResult& r)
{
    std::ofstream f(filename);
    f << "# CeeLo simulation results\n"
      << "energy_keV,fep_efficiency,fep_uncertainty,"
         "total_efficiency,total_uncertainty,num_events\n"
      << std::fixed << std::setprecision(4) << energy_keV << ","
      << std::scientific << std::setprecision(6)
      << r.full_energy_peak_efficiency << "," << r.fep_uncertainty << ","
      << r.total_efficiency << "," << r.total_uncertainty << ","
      << r.num_events_simulated << "\n";
}

/// Write a multi-row CSV for a batch of energies.
static void write_our_csv_batch(const std::string& filename,
                                const std::vector<double>& energies_keV,
                                const std::vector<EfficiencyResult>& results)
{
    std::ofstream f(filename);
    f << "# CeeLo simulation results\n"
      << "energy_keV,fep_efficiency,fep_uncertainty,"
         "total_efficiency,total_uncertainty,num_events\n";
    for (size_t i = 0; i < results.size(); ++i) {
        const auto& r = results[i];
        f << std::fixed << std::setprecision(4) << energies_keV[i] << ","
          << std::scientific << std::setprecision(6)
          << r.full_energy_peak_efficiency << "," << r.fep_uncertainty << ","
          << r.total_efficiency << "," << r.total_uncertainty << ","
          << r.num_events_simulated << "\n";
    }
}

// ---------------------------------------------------------------------------
// Energy grid: §9.3 energies (skip below threshold for each detector)
// ---------------------------------------------------------------------------
static const std::vector<double> ENERGIES_NAI = {
    30, 50, 80, 100, 200, 300, 500, 662, 800, 1000, 1173, 1332, 2000, 3000
};
static const std::vector<double> ENERGIES_HPGE = {
    40, 60, 80, 100, 200, 300, 500, 662, 1000, 1173, 1332, 2000, 3000, 5000, 7000, 10000
};
static const std::vector<double> ENERGIES_CZT = {
    30, 50, 80, 100, 200, 300, 500, 662, 800, 1000, 1500
};

int main()
{
    // -----------------------------------------------------------------------
    // Config 1: 3"×3" NaI cylinder, no shielding, on-axis point, 10 cm
    // -----------------------------------------------------------------------
    {
        const int cfg = 1;
        std::cout << "\n=== Config " << cfg
                  << ": 3\"x3\" NaI, no shielding, on-axis 10 cm ===\n";

        Material nai = make_NaI();
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));

        // Export GDML (same geometry for all energies in this config).
        calc.export_geant4_gdml("detector_1.gdml");

        auto results = calc.compute_batch(ENERGIES_NAI, N_MC, N_THREADS);
        write_our_csv_batch("our_1_multi.csv", ENERGIES_NAI, results);

        // Write one macro per energy for the G4 harness.
        for (double E : ENERGIES_NAI) {
            std::string mac = "run_1_" + std::to_string(static_cast<int>(E)) + "keV.mac";
            calc.export_geant4_macro(mac, E, N_G4);
        }

        std::cout << "  GDML:    detector_1.gdml\n"
                  << "  Our CSV: our_1_multi.csv\n"
                  << "  Macros:  run_1_<E>keV.mac\n";
        for (size_t i = 0; i < ENERGIES_NAI.size(); ++i) {
            std::cout << "    " << std::setw(7) << ENERGIES_NAI[i] << " keV:"
                      << "  ε_FEP=" << std::scientific << std::setprecision(3)
                      << results[i].full_energy_peak_efficiency
                      << "  ε_tot=" << results[i].total_efficiency << "\n";
        }
    }

    // -----------------------------------------------------------------------
    // Config 2: 3"×3" NaI, 1 mm Al endcap, on-axis point, 10 cm
    // -----------------------------------------------------------------------
    {
        const int cfg = 2;
        std::cout << "\n=== Config " << cfg
                  << ": 3\"x3\" NaI, 1mm Al endcap, on-axis 10 cm ===\n";

        Material nai = make_NaI();
        Material al  = make_Aluminum();

        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        // 1 mm Al front + side attenuator covering the detector face
        calc.add_attenuator(&al, 0.1, 0.1, 0.0, 7.62);
        calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));

        calc.export_geant4_gdml("detector_2.gdml");

        auto results = calc.compute_batch(ENERGIES_NAI, N_MC, N_THREADS);
        write_our_csv_batch("our_2_multi.csv", ENERGIES_NAI, results);

        for (double E : ENERGIES_NAI) {
            std::string mac = "run_2_" + std::to_string(static_cast<int>(E)) + "keV.mac";
            calc.export_geant4_macro(mac, E, N_G4);
        }

        std::cout << "  GDML:    detector_2.gdml\n"
                  << "  Our CSV: our_2_multi.csv\n";
    }

    // -----------------------------------------------------------------------
    // Config 3: 2"×2" LaBr3, 0.5 mm Al, on-axis point, 5 cm
    // -----------------------------------------------------------------------
    {
        const int cfg = 3;
        std::cout << "\n=== Config " << cfg
                  << ": 2\"x2\" LaBr3, 0.5mm Al, on-axis 5 cm ===\n";

        Material labr3 = make_LaBr3();
        Material al    = make_Aluminum();

        EfficiencyCalculator calc;
        // 2"×2": R=2.54 cm (1"), L=5.08 cm (2")
        calc.set_detector(DetectorShape::Cylinder, &labr3, {2.54, 5.08});
        calc.add_attenuator(&al, 0.05, 0.05, 0.0, 5.08);
        calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -5.0));

        calc.export_geant4_gdml("detector_3.gdml");

        auto results = calc.compute_batch(ENERGIES_NAI, N_MC, N_THREADS);
        write_our_csv_batch("our_3_multi.csv", ENERGIES_NAI, results);

        for (double E : ENERGIES_NAI) {
            std::string mac = "run_3_" + std::to_string(static_cast<int>(E)) + "keV.mac";
            calc.export_geant4_macro(mac, E, N_G4);
        }

        std::cout << "  GDML:    detector_3.gdml\n"
                  << "  Our CSV: our_3_multi.csv\n";
    }

    // -----------------------------------------------------------------------
    // Config 6: 3"×3" NaI, no shielding, off-axis point (45°), 15 cm
    //   Source at 45° from axis: pos = 15*sin(45°), 0, -15*cos(45°) = (10.6, 0, -10.6)
    // -----------------------------------------------------------------------
    {
        const int cfg = 6;
        std::cout << "\n=== Config " << cfg
                  << ": 3\"x3\" NaI, no shielding, off-axis 45deg, 15 cm ===\n";

        Material nai = make_NaI();

        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});

        const double d = 15.0;
        const double theta = 45.0 * M_PI / 180.0;
        calc.set_point_source(Eigen::Vector3d(d * std::sin(theta), 0.0, -d * std::cos(theta)));

        calc.export_geant4_gdml("detector_6.gdml");

        const std::vector<double> energies_6 = {100, 200, 300, 500, 662, 1000, 1173, 1332, 2000};
        auto results = calc.compute_batch(energies_6, N_MC, N_THREADS);
        write_our_csv_batch("our_6_multi.csv", energies_6, results);

        for (double E : energies_6) {
            std::string mac = "run_6_" + std::to_string(static_cast<int>(E)) + "keV.mac";
            calc.export_geant4_macro(mac, E, N_G4);
        }

        std::cout << "  GDML:    detector_6.gdml\n"
                  << "  Our CSV: our_6_multi.csv\n";
    }

    // -----------------------------------------------------------------------
    // Config 7: 3"×3" NaI, 1mm Al + 2mm Pb shielding, on-axis, 15 cm
    // -----------------------------------------------------------------------
    {
        const int cfg = 7;
        std::cout << "\n=== Config " << cfg
                  << ": 3\"x3\" NaI, 1mm Al + 2mm Pb, on-axis 15 cm ===\n";

        Material nai = make_NaI();
        Material al  = make_Aluminum();
        Material pb  = make_Lead();

        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc.add_attenuator(&al, 0.1, 0.1, 0.0, 7.62);   // innermost: 1mm Al
        calc.add_attenuator(&pb, 0.2, 0.2, 0.0, 7.62);   // outer: 2mm Pb
        calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -15.0));

        calc.export_geant4_gdml("detector_7.gdml");

        const std::vector<double> energies_7 = {100, 200, 300, 500, 662, 1000, 1173, 1332, 2000, 3000};
        auto results = calc.compute_batch(energies_7, N_MC, N_THREADS);
        write_our_csv_batch("our_7_multi.csv", energies_7, results);

        for (double E : energies_7) {
            std::string mac = "run_7_" + std::to_string(static_cast<int>(E)) + "keV.mac";
            calc.export_geant4_macro(mac, E, N_G4);
        }

        std::cout << "  GDML:    detector_7.gdml\n"
                  << "  Our CSV: our_7_multi.csv\n";
    }

    // -----------------------------------------------------------------------
    // Config 5: 1cm×1cm×5mm CZT, no shielding, on-axis, 5 cm
    //   Box detector: half_x=0.5, half_y=0.5, length=0.5 cm
    // -----------------------------------------------------------------------
    {
        const int cfg = 5;
        std::cout << "\n=== Config " << cfg
                  << ": 1x1x0.5cm CZT box, no shielding, on-axis 5 cm ===\n";

        Material czt = make_CZT();

        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Box, &czt, {0.5, 0.5, 0.5});
        calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -5.0));

        calc.export_geant4_gdml("detector_5.gdml");

        auto results = calc.compute_batch(ENERGIES_CZT, N_MC, N_THREADS);
        write_our_csv_batch("our_5_multi.csv", ENERGIES_CZT, results);

        for (double E : ENERGIES_CZT) {
            std::string mac = "run_5_" + std::to_string(static_cast<int>(E)) + "keV.mac";
            calc.export_geant4_macro(mac, E, N_G4);
        }

        std::cout << "  GDML:    detector_5.gdml\n"
                  << "  Our CSV: our_5_multi.csv\n";
    }

    std::cout << "\n=== All benchmark files generated. ===\n"
              << "Now run the GEANT4 harness for each config:\n"
              << "  ./ceelo_g4val detector_1.gdml run_1_662keV.mac g4_1_662keV.csv\n"
              << "Then compare:\n"
              << "  python ../tools/geant4_validation/compare_results.py our_1_multi.csv g4_1_multi.csv\n"
              << "\nNote: GEANT4 produces one CSV per run/macro; combine multiple energies manually\n"
              << "      or run the batch comparison script.\n\n";

    return 0;
}
