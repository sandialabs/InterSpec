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

/// @file generate_all_spectra.cpp
/// @brief Generate MC spectra and GDML+macro files for all 6 validation configs.
///
/// For each config, runs our MC code at representative energies and writes:
///   - mc_cfg{N}_{E}keV.csv          histogram (RealTime + Energy (keV),Counts)
///   - mc_cfg{N}_multi.csv           multi-energy efficiency table
///   - detector_{N}.gdml             GDML geometry for G4
///   - run_{N}_{E}keV.mac            G4 macro (isotropic, high stats)
///
/// Usage:
///   generate_all_spectra [num_events_per_energy]
///
/// Build:
///   cmake --build build -j8 --target generate_all_spectra

#include "efficiency/EfficiencyCalculator.h"
#include "materials/Material.h"
#include "transport/PhotonTransport.h"
#include "geometry/Geometry.h"
#include "geometry/SourceGeometry.h"
#include "export/Geant4Export.h"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <cstdint>
#include <cmath>
#include <random>
#include <chrono>
#include <thread>
#include <mutex>

using namespace ceelo;

// ---------------------------------------------------------------------------
// Config descriptor
// ---------------------------------------------------------------------------
struct ConfigDesc {
    int id;
    std::string label;
    std::vector<double> energies_keV;
};

// ---------------------------------------------------------------------------
// Histogram accumulator
// ---------------------------------------------------------------------------
struct HistogramResult {
    std::vector<uint64_t> counts;
    double bin_width_keV = 1.0;
    uint64_t n_fep = 0;
    uint64_t n_any = 0;
    uint64_t n_total = 0;
};

// ---------------------------------------------------------------------------
// Sample a random direction uniformly on the unit sphere (isotropic 4pi).
// ---------------------------------------------------------------------------
static Eigen::Vector3d sample_isotropic_direction(
    std::mt19937_64& rng, std::uniform_real_distribution<double>& uniform)
{
    double cos_theta = 2.0 * uniform(rng) - 1.0;
    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
    double phi = 2.0 * M_PI * uniform(rng);
    return Eigen::Vector3d(sin_theta * std::cos(phi),
                           sin_theta * std::sin(phi),
                           cos_theta);
}

// ---------------------------------------------------------------------------
// Run a histogram using transport_photon directly (true integer counts).
// When use_cone_bias=true, samples within a cone toward the detector.
// When use_cone_bias=false, samples isotropically into 4pi.
// Stores UNweighted histogram counts (one count per event with deposit).
// ---------------------------------------------------------------------------
static HistogramResult run_histogram(
    const EfficiencyCalculator& calc,
    double energy_keV,
    uint64_t num_events,
    bool use_cone_bias = true,
    double bin_width_keV = 1.0)
{
    const Geometry& geom = calc.geometry();
    TransportConfig config;
    config.enable_electron_csda = true;

    int n_bins = static_cast<int>(std::ceil((energy_keV + 10.0) / bin_width_keV));

    HistogramResult hist;
    hist.counts.resize(n_bins, 0);
    hist.bin_width_keV = bin_width_keV;
    hist.n_total = num_events;

    constexpr double kFepTol = 1.5; // keV

    unsigned n_threads = std::thread::hardware_concurrency();
    if (n_threads < 1) n_threads = 1;
    uint64_t per_thread = num_events / n_threads;

    // Get source position and compute cone parameters (used only if cone biasing)
    Eigen::Vector3d source_pos = calc.source_position();
    double cos_theta_max = -1.0; // full sphere default
    Eigen::Vector3d cone_axis(0, 0, 1);

    if (use_cone_bias) {
        double det_r = geom.outer_bounding_radius();
        auto [z_min, z_max] = geom.outer_z_extent();

        // Cone axis: source toward detector center
        double det_center_z = 0.5 * (z_min + z_max);
        Eigen::Vector3d det_center(0.0, 0.0, det_center_z);
        Eigen::Vector3d to_det = det_center - source_pos;
        double dist = to_det.norm();
        cone_axis = to_det / dist;

        // Cone half-angle: must subtend the ENTIRE detector bounding volume.
        // Use distance to nearest face (not center) for the angle calculation,
        // matching EfficiencyCalculator::compute_cone_half_angle().
        double sx = source_pos.x(), sy = source_pos.y(), sz = source_pos.z();
        double lateral_offset = std::sqrt(sx * sx + sy * sy);
        double apparent_r = det_r + lateral_offset;

        double dz = std::abs(sz - z_min); // default: source in front
        if (sz > z_max) dz = std::abs(sz - z_max); // source behind
        if (sz >= z_min && sz <= z_max) dz = 0.0; // source inside

        double half_angle = (dz < 1e-6) ? M_PI
            : std::min(std::atan2(apparent_r, dz) * 1.05, M_PI);
        cos_theta_max = std::cos(half_angle);
    }

    std::mutex mtx;
    std::vector<std::thread> threads;

    for (unsigned t = 0; t < n_threads; ++t) {
        uint64_t n = (t == n_threads - 1) ? (num_events - per_thread * (n_threads - 1)) : per_thread;
        threads.emplace_back([&, n, t]() {
            std::mt19937_64 rng(42 + t * 997);
            std::uniform_real_distribution<double> uniform(0.0, 1.0);

            std::vector<uint64_t> local_counts(n_bins, 0);
            uint64_t local_fep = 0, local_any = 0;

            // Build orthonormal basis for cone (only used if cone biasing)
            Eigen::Vector3d w = cone_axis;
            Eigen::Vector3d u, v;
            if (use_cone_bias) {
                Eigen::Vector3d ex(1, 0, 0), ey(0, 1, 0);
                if (std::abs(w.x()) < 0.9)
                    u = ex.cross(w).normalized();
                else
                    u = ey.cross(w).normalized();
                v = w.cross(u);
            }

            for (uint64_t i = 0; i < n; ++i) {
                Eigen::Vector3d dir;
                if (use_cone_bias) {
                    double cos_theta = 1.0 - uniform(rng) * (1.0 - cos_theta_max);
                    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
                    double phi = 2.0 * M_PI * uniform(rng);
                    dir = sin_theta * std::cos(phi) * u
                        + sin_theta * std::sin(phi) * v
                        + cos_theta * w;
                } else {
                    dir = sample_isotropic_direction(rng, uniform);
                }

                auto result = transport_photon(source_pos, dir, energy_keV,
                                               geom, config, rng);

                double edep = result.energy_deposited_scoring;
                if (edep > 0.0) {
                    local_any++;
                    if (std::abs(edep - energy_keV) < kFepTol) {
                        local_fep++;
                    }
                    int bin = static_cast<int>(edep / bin_width_keV);
                    if (bin >= 0 && bin < n_bins) {
                        local_counts[bin]++;
                    }
                }
            }

            std::lock_guard<std::mutex> lk(mtx);
            for (int i = 0; i < n_bins; ++i) {
                hist.counts[i] += local_counts[i];
            }
            hist.n_fep += local_fep;
            hist.n_any += local_any;
        });
    }

    for (auto& th : threads) th.join();
    return hist;
}

// ---------------------------------------------------------------------------
// Run a histogram for Marinelli beaker: samples positions inside the L-shaped
// volume, transports through source material + beaker wall, then detector.
// ---------------------------------------------------------------------------
static HistogramResult run_histogram_marinelli(
    EfficiencyCalculator& calc,
    double energy_keV,
    uint64_t num_events,
    double bin_width_keV = 1.0)
{
    // Use compute() which includes all Marinelli environmental fixes
    // (miss-bounce, re-entry bounce, escaped secondaries, PP in water,
    // full beaker wall transport).
    int n_bins = static_cast<int>(std::ceil((energy_keV + 10.0) / bin_width_keV));
    std::vector<float> bin_edges(n_bins + 1);
    for (int i = 0; i <= n_bins; ++i)
        bin_edges[i] = static_cast<float>(i * bin_width_keV);

    auto result = calc.compute(energy_keV, num_events, 0, bin_edges);

    HistogramResult hist;
    hist.bin_width_keV = bin_width_keV;
    uint64_t N = result.num_events_simulated;
    hist.n_total = N;

    // Convert pulse_height_distribution (weighted probabilities) back to
    // effective counts for the histogram writer.
    // Use num_events_simulated (not num_events) for consistency since
    // compute() may run more events than requested for precision batching.
    hist.counts.resize(n_bins, 0);
    for (int i = 0; i < n_bins; ++i) {
        if (i < static_cast<int>(result.pulse_height_distribution.size()))
            hist.counts[i] = static_cast<uint64_t>(
                result.pulse_height_distribution[i] * N + 0.5);
    }

    // Derive FEP/total counts from efficiencies
    hist.n_fep = static_cast<uint64_t>(
        result.full_energy_peak_efficiency * N + 0.5);
    hist.n_any = static_cast<uint64_t>(
        result.total_efficiency * N + 0.5);

    return hist;
}

// ---------------------------------------------------------------------------
// Write histogram CSV in the standard format
// ---------------------------------------------------------------------------
static void write_histogram_csv(const std::string& filename,
                                 const HistogramResult& hist,
                                 double energy_keV,
                                 const std::string& config_label,
                                 double omega_frac)
{
    // RealTime: N_events / (1E6 * omega_frac) for cone-biased
    constexpr double kActivity = 1.0e6;
    double real_time_s = static_cast<double>(hist.n_total) / (kActivity * omega_frac);

    double eps_fep = static_cast<double>(hist.n_fep) / hist.n_total * omega_frac;
    double eps_any = static_cast<double>(hist.n_any) / hist.n_total * omega_frac;
    double sig_fep = std::sqrt(eps_fep * (1.0 - eps_fep / omega_frac) / hist.n_total) * omega_frac;
    double sig_any = std::sqrt(eps_any * (1.0 - eps_any / omega_frac) / hist.n_total) * omega_frac;

    std::ofstream csv(filename);
    csv << "# MC energy deposit histogram\n"
        << "# " << config_label << "\n"
        << "# Primary energy: " << energy_keV << " keV\n"
        << "# Events: " << hist.n_total << "\n"
        << "# Bin width: " << hist.bin_width_keV << " keV\n"
        << "# Absolute FEP: " << eps_fep << " +/- " << sig_fep << "\n"
        << "# Absolute total: " << eps_any << " +/- " << sig_any << "\n"
        << "RealTime: " << std::fixed << std::setprecision(3) << real_time_s << " s\n"
        << "Energy (keV),Counts\n";

    int n_bins = static_cast<int>(hist.counts.size());
    for (int i = 0; i < n_bins; ++i) {
        if (hist.counts[i] > 0) {
            double lower_keV = i * hist.bin_width_keV;
            csv << std::fixed << std::setprecision(1) << lower_keV
                << "," << hist.counts[i] << "\n";
        }
    }
}

// ---------------------------------------------------------------------------
// Write multi-energy efficiency CSV
// ---------------------------------------------------------------------------
static void write_multi_csv(const std::string& filename,
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
// Run one config: generate GDML, macros, MC histograms, MC efficiencies
//
// use_cone_bias: true  = point source with no source shielding (cone sample)
//                false = extended source or source with shielding (isotropic 4pi)
// ---------------------------------------------------------------------------
static void run_config(EfficiencyCalculator& calc,
                        int cfg_id,
                        const std::string& label,
                        const std::vector<double>& energies_keV,
                        const std::vector<double>& histogram_energies_keV,
                        uint64_t n_mc_multi,
                        uint64_t n_mc_histo,
                        uint64_t n_g4,
                        bool use_cone_bias = true)
{
    std::string cfg_str = std::to_string(cfg_id);

    std::cout << "\n=== Config " << cfg_id << ": " << label << " ===" << std::endl;

    // 1) Export GDML (vacuum world so G4 doesn't include air attenuation —
    //    matches MC default of AirAttenuation::None)
    std::string gdml = "detector_" + cfg_str + ".gdml";
    calc.export_geant4_gdml(gdml, /*vacuum_world=*/true);
    std::cout << "  GDML: " << gdml << std::endl;

    // 2) Generate G4 macros (isotropic, high stats) for efficiency AND histogram energies
    {
        std::set<int> done_energies;
        auto gen_mac = [&](double E) {
            int e_int = static_cast<int>(E);
            if (done_energies.count(e_int)) return;
            done_energies.insert(e_int);
            std::string mac = "run_" + cfg_str + "_" + std::to_string(e_int) + "keV.mac";
            calc.export_geant4_macro(mac, E, n_g4);
        };
        for (double E : energies_keV) gen_mac(E);
        for (double E : histogram_energies_keV) gen_mac(E);
        std::cout << "  G4 macros: run_" << cfg_str << "_*keV.mac ("
                  << done_energies.size() << " energies, " << n_g4 << " events each)" << std::endl;
    }

    // 3) Run our MC for multi-energy efficiency curve
    std::cout << "  MC efficiency (" << n_mc_multi << " events/energy, "
              << energies_keV.size() << " energies)..." << std::flush;
    auto t_eff_0 = std::chrono::steady_clock::now();
    auto results = calc.compute_batch(energies_keV, n_mc_multi, 0 /*auto threads*/);
    auto t_eff_1 = std::chrono::steady_clock::now();
    double eff_elapsed = std::chrono::duration<double>(t_eff_1 - t_eff_0).count();
    std::string multi_csv = "mc_cfg" + cfg_str + "_multi.csv";
    write_multi_csv(multi_csv, energies_keV, results);
    std::cout << " done (" << std::fixed << std::setprecision(1)
              << eff_elapsed << " s) -> " << multi_csv << std::endl;

    // 4) Run our MC histograms at representative energies
    // Compute omega_frac for RealTime calculation
    double omega_frac = 1.0; // isotropic default
    if (use_cone_bias) {
        Eigen::Vector3d source_pos = calc.source_position();
        const Geometry& geom = calc.geometry();
        double det_r = geom.outer_bounding_radius();
        auto [z_min, z_max] = geom.outer_z_extent();

        // Use distance to nearest face (matching run_histogram cone and
        // EfficiencyCalculator::compute_cone_half_angle).
        double sx = source_pos.x(), sy = source_pos.y(), sz = source_pos.z();
        double lateral_offset = std::sqrt(sx * sx + sy * sy);
        double apparent_r = det_r + lateral_offset;

        double dz = std::abs(sz - z_min);
        if (sz > z_max) dz = std::abs(sz - z_max);
        if (sz >= z_min && sz <= z_max) dz = 0.0;

        double half_angle = (dz < 1e-6) ? M_PI
            : std::min(std::atan2(apparent_r, dz) * 1.05, M_PI);
        double cos_theta_max_local = std::cos(half_angle);
        omega_frac = (1.0 - cos_theta_max_local) / 2.0;
    }

    bool has_source_effects = calc.source_geometry().has_source_effects();
    bool use_compute_path = has_source_effects;  // Marinelli, extended sources, source shields
    std::string bias_label = use_compute_path ? "compute" : (use_cone_bias ? "cone" : "isotropic");
    for (double E : histogram_energies_keV) {
        int e_int = static_cast<int>(E);
        std::cout << "  Histogram " << e_int << " keV (" << n_mc_histo
                  << " events, " << bias_label << ")..." << std::flush;

        auto t0 = std::chrono::steady_clock::now();
        HistogramResult hist;
        if (use_compute_path) {
            // Use compute() for any config with source material/shielding.
            // This ensures source transport (attenuation, scattering in source
            // material and shields) is properly included in the histogram.
            hist = run_histogram_marinelli(calc, E, n_mc_histo);
        } else {
            hist = run_histogram(calc, E, n_mc_histo, use_cone_bias);
        }
        auto t1 = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(t1 - t0).count();

        std::string histo_file = "mc_cfg" + cfg_str + "_" + std::to_string(e_int) + "keV.csv";
        write_histogram_csv(histo_file, hist, E, label, omega_frac);

        double eps_fep = static_cast<double>(hist.n_fep) / hist.n_total * omega_frac;
        std::cout << " done (" << std::fixed << std::setprecision(1) << elapsed
                  << " s, FEP=" << std::scientific << std::setprecision(3) << eps_fep << ")" << std::endl;
    }
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char** argv) {
    uint64_t n_mc_histo = 4'000'000; // events per histogram energy
    uint64_t n_mc_multi = 100'000;   // events per energy for efficiency curve (cone-biased → plenty of stats)
    uint64_t n_g4 = 100'000'000;     // events per G4 run (high stats)

    if (argc >= 2) n_mc_histo = static_cast<uint64_t>(std::atoll(argv[1]));
    if (argc >= 3) n_g4 = static_cast<uint64_t>(std::atoll(argv[2]));
    if (argc >= 4) n_mc_multi = static_cast<uint64_t>(std::atoll(argv[3]));

    std::cout << "MC histogram events: " << n_mc_histo << "\n"
              << "MC efficiency events: " << n_mc_multi << " per energy\n"
              << "G4 macro events: " << n_g4 << " per energy\n" << std::flush;

    // Track which configs use cone biasing (for G4 run script generation).
    // Point sources with no source shielding use cone bias; extended sources
    // or sources with shielding use isotropic 4pi emission.
    std::set<int> cone_bias_configs;

    // Histogram energies: representative subset for spectral comparison
    const std::vector<double> HISTO_ENERGIES = {59, 122, 185, 375, 662, 1173, 2614};
    const std::vector<double> HISTO_CZT = {59, 122, 185, 375, 662};

    // Full energy grids for efficiency curves
    const std::vector<double> E_NAI = {30, 50, 80, 100, 200, 300, 500, 662, 800, 1000, 1173, 1332, 2000, 3000};
    const std::vector<double> E_CZT = {30, 50, 80, 100, 200, 300, 500, 662, 800, 1000, 1500};

    // -----------------------------------------------------------------------
    // Config 1: 3"x3" NaI, bare, on-axis 10 cm
    // -----------------------------------------------------------------------
    {
        Material nai = make_NaI();
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));

        cone_bias_configs.insert(1);
        run_config(calc, 1, "3\"x3\" NaI bare, on-axis 10cm",
                   E_NAI, HISTO_ENERGIES, n_mc_multi, n_mc_histo, n_g4, true);
    }

    // -----------------------------------------------------------------------
    // Config 2: 3"x3" NaI, 1mm Al, on-axis 10 cm
    // -----------------------------------------------------------------------
    {
        Material nai = make_NaI();
        Material al = make_Aluminum();
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc.add_attenuator(&al, 0.1, 0.1, 0.0, 7.62);
        calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));

        cone_bias_configs.insert(2);
        run_config(calc, 2, "3\"x3\" NaI + 1mm Al, on-axis 10cm",
                   E_NAI, HISTO_ENERGIES, n_mc_multi, n_mc_histo, n_g4, true);
    }

    // -----------------------------------------------------------------------
    // Config 3: 2"x2" LaBr3, 0.5mm Al, on-axis 5 cm
    // -----------------------------------------------------------------------
    {
        Material labr3 = make_LaBr3();
        Material al = make_Aluminum();
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &labr3, {2.54, 5.08});
        calc.add_attenuator(&al, 0.05, 0.05, 0.0, 5.08);
        calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -5.0));

        cone_bias_configs.insert(3);
        run_config(calc, 3, "2\"x2\" LaBr3 + 0.5mm Al, on-axis 5cm",
                   E_NAI, HISTO_ENERGIES, n_mc_multi, n_mc_histo, n_g4, true);
    }

    // -----------------------------------------------------------------------
    // Config 5: 1cm x 1cm x 0.5cm CZT, bare, on-axis 5 cm
    // -----------------------------------------------------------------------
    {
        Material czt = make_CZT();
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Box, &czt, {0.5, 0.5, 0.5});
        calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -5.0));

        cone_bias_configs.insert(5);
        run_config(calc, 5, "1x1x0.5cm CZT bare, on-axis 5cm",
                   E_CZT, HISTO_CZT, n_mc_multi, n_mc_histo, n_g4, true);
    }

    // -----------------------------------------------------------------------
    // Config 6: 3"x3" NaI, bare, off-axis 45deg, 15 cm
    // -----------------------------------------------------------------------
    {
        Material nai = make_NaI();
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        const double d = 15.0;
        const double theta = 45.0 * M_PI / 180.0;
        calc.set_point_source(Eigen::Vector3d(d * std::sin(theta), 0.0, -d * std::cos(theta)));

        const std::vector<double> E_6 = {100, 200, 300, 500, 662, 1000, 1173, 1332, 2000};
        const std::vector<double> H_6 = {122, 375, 662, 1173};
        cone_bias_configs.insert(6);
        run_config(calc, 6, "3\"x3\" NaI bare, off-axis 45deg 15cm",
                   E_6, H_6, n_mc_multi, n_mc_histo, n_g4, true);
    }

    // -----------------------------------------------------------------------
    // Config 7: 3"x3" NaI, 1mm Al + 2mm Pb, on-axis 15 cm
    // Note: Al and Pb are *detector attenuators* (not source shielding),
    // so cone biasing is still appropriate.
    // -----------------------------------------------------------------------
    {
        Material nai = make_NaI();
        Material al = make_Aluminum();
        Material pb = make_Lead();
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc.add_attenuator(&al, 0.1, 0.1, 0.0, 7.62);
        calc.add_attenuator(&pb, 0.2, 0.2, 0.0, 7.62);
        calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -15.0));

        const std::vector<double> E_7 = {100, 200, 300, 500, 662, 1000, 1173, 1332, 2000, 3000};
        const std::vector<double> H_7 = {185, 375, 662, 1173};
        cone_bias_configs.insert(7);
        run_config(calc, 7, "3\"x3\" NaI + 1mm Al + 2mm Pb, on-axis 15cm",
                   E_7, H_7, n_mc_multi, n_mc_histo, n_g4, true);
    }

    // -----------------------------------------------------------------------
    // Config 8: 3"x3" NaI + 0.5mm Al housing, Marinelli beaker (water, 2mm PE)
    // -----------------------------------------------------------------------
    {
        Material nai = make_NaI();
        Material al = make_Aluminum();
        Material water = make_Water();
        Material pe = make_Polyethylene();
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        // 0.5mm aluminum detector housing (industry standard for NaI)
        calc.add_attenuator(&al, 0.05, 0.05, 0.0, 7.62);
        // Marinelli: well_r=4.3cm, well_depth=6cm, outer_r=7.5cm,
        //            fill_height=4cm, endcap_to_beaker=0.5cm
        //            2mm polyethylene beaker wall, water sample
        calc.set_marinelli_beaker(
            /*well_inner_radius=*/4.3,
            /*well_depth=*/6.0,
            /*outer_radius=*/7.5,
            /*fill_height=*/4.0,
            /*endcap_to_beaker=*/0.5,
            &water, &pe, 0.2);
        calc.enable_source_electron_transport(true);

        const std::vector<double> E_8 = {59, 100, 200, 300, 500, 662, 1000, 1173, 1332, 2614};
        const std::vector<double> H_8 = {59, 200, 662, 1173, 2614};
        // No cone bias for Marinelli (source surrounds detector)
        run_config(calc, 8, "3\"x3\" NaI + 0.5mm Al, Marinelli beaker (water, 2mm PE)",
                   E_8, H_8, n_mc_multi, n_mc_histo, n_g4, false);
    }

    // Config 9: 3"x3" NaI + 0.5mm Al housing, Marinelli beaker (soil, 2mm Pyrex)
    // -----------------------------------------------------------------------
    {
        Material nai = make_NaI();
        Material al = make_Aluminum();
        Material soil = make_Soil();
        Material pyrex = make_Pyrex();
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc.add_attenuator(&al, 0.05, 0.05, 0.0, 7.62);
        calc.set_marinelli_beaker(
            /*well_inner_radius=*/4.3,
            /*well_depth=*/6.0,
            /*outer_radius=*/7.5,
            /*fill_height=*/4.0,
            /*endcap_to_beaker=*/0.5,
            &soil, &pyrex, 0.2);
        calc.enable_source_electron_transport(true);

        const std::vector<double> E_9 = {59, 100, 200, 300, 500, 662, 1000, 1173, 1332, 2614};
        const std::vector<double> H_9 = {59, 200, 662, 1173, 2614};
        // No cone bias for Marinelli (source surrounds detector)
        run_config(calc, 9, "3\"x3\" NaI + 0.5mm Al, Marinelli beaker (soil, 2mm Pyrex)",
                   E_9, H_9, n_mc_multi, n_mc_histo, n_g4, false);
    }

    // -----------------------------------------------------------------------
    // Config 10: In-situ HPGe with W collimator, soil cylinder source
    // Realistic coaxial HPGe: 60mm dia x 60mm length, bore hole, dead layer,
    // Al endcap, W collimator extending 3 inches (7.62 cm) past crystal front.
    // Source: soil cylinder (r=500cm, 30cm deep) at 100cm below detector.
    // -----------------------------------------------------------------------
    {
        Material hpge = make_HPGe();
        Material al = make_Aluminum();
        Material w = make_Tungsten();
        Material soil = make_Soil();
        EfficiencyCalculator calc;

        // HPGe crystal: 3cm radius x 6cm length (typical 50% relative efficiency)
        calc.set_detector(DetectorShape::Cylinder, &hpge, {3.0, 6.0});

        // Bore hole: 0.5cm radius, 4cm depth from back face
        calc.set_bore_hole(0.5, 4.0);

        // Dead layer: 0.7mm Li contact (front), 0.03mm outer contact (side)
        calc.set_dead_layer(0.07, 0.003);

        // Al endcap: 1mm front, 1mm side
        calc.add_attenuator(&al, 0.1, 0.1, 0.0, 6.0);

        // W collimator: 2cm side thickness, extends 3 inches (7.62cm) past crystal front.
        // z_start = -7.62 (extends forward toward source), z_end = 6.0 (crystal back).
        // Side-only tube (no front/back disk) via add_collimator.
        calc.add_collimator(&w, 2.0, -7.62, 6.0);

        // Soil source: cylinder r=500cm, 30cm deep, centered 115cm below detector
        // (detector front at z=0; soil surface at z=-100cm, soil extends to z=-130cm)
        // Center of source cylinder at z=-115cm, half_length = 15cm
        calc.set_cylindrical_source(
            Eigen::Vector3d(0.0, 0.0, -115.0),  // center
            500.0,     // radius (5 meters)
            15.0,      // half-length (30cm deep)
            Eigen::Matrix3d::Identity());
        calc.set_source_material(&soil);

        const std::vector<double> E_10 = {59, 100, 200, 300, 500, 662, 1000, 1173, 1332, 2614};
        const std::vector<double> H_10 = {59, 200, 662, 1173, 2614};
        // No cone bias: extended source with source material
        run_config(calc, 10, "HPGe + W collimator, in-situ soil (r=5m, d=1m)",
                   E_10, H_10, n_mc_multi, n_mc_histo, n_g4, false);
    }

    // -----------------------------------------------------------------------
    // Config 11: 3"x3" NaI, point source with 0.5cm Fe source shielding
    // -----------------------------------------------------------------------
    {
        Material nai = make_NaI();
        Material fe = make_Iron();
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
        calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
        calc.add_source_shield(&fe, 0.5);  // 0.5 cm iron source shielding

        const std::vector<double> E_11 = {100, 200, 300, 500, 662, 1000, 1173, 1332, 2000, 3000};
        const std::vector<double> H_11 = {200, 375, 662, 1173};
        // No cone bias: source shielding present in full mode
        run_config(calc, 11, "3\"x3\" NaI, point source + 0.5cm Fe shield",
                   E_11, H_11, n_mc_multi, n_mc_histo, n_g4, false);
    }

    // -----------------------------------------------------------------------
    // Config 12: 3"x3" NaI, 10x15x20cm steel box filled with cellulose source
    // Box is 2mm stainless steel, filled with cellulose (0.5 g/cm³).
    // Centered 15cm in front of detector.
    // -----------------------------------------------------------------------
    {
        Material nai = make_NaI();
        Material ss = make_StainlessSteel304();
        Material cellulose = make_Cellulose();
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});

        // Rectangular source: 10x15x20 cm box → half-dims 5x7.5x10 cm
        // Center the source 15cm + half-depth in front of detector.
        // Detector front at z=0, source center at z = -(15 + 10) = -25.
        calc.set_rectangular_source(
            Eigen::Vector3d(0.0, 0.0, -25.0),  // center
            Eigen::Vector3d(5.0, 7.5, 10.0),   // half-dims (10x15x20 cm)
            Eigen::Matrix3d::Identity());
        calc.set_source_material(&cellulose);

        // Steel box walls: 2mm stainless steel as source shielding
        calc.add_source_shield(&ss, 0.2);

        const std::vector<double> E_12 = {59, 100, 200, 300, 500, 662, 1000, 1173, 1332, 2614};
        const std::vector<double> H_12 = {59, 200, 662, 1173};
        // No cone bias: extended source with source material + shielding
        run_config(calc, 12, "3\"x3\" NaI, 10x15x20cm SS box + cellulose",
                   E_12, H_12, n_mc_multi, n_mc_histo, n_g4, false);
    }

    // -----------------------------------------------------------------------
    // Write G4 run script
    // -----------------------------------------------------------------------
    // Track which configs use cone biasing (point source, no source shielding)
    // vs isotropic 4pi (extended source or source with shielding).
    // Isotropic configs need higher G4 stats since most photons miss the detector.
    {
        std::ofstream sh("run_g4_all.sh");
        sh << "#!/bin/bash\n"
           << "# Run all GEANT4 simulations.\n"
           << "# Generated by generate_all_spectra\n"
           << "#\n"
           << "# Configs with point sources (no source shielding) use --cone-bias.\n"
           << "# Configs with extended sources or source shielding use isotropic 4pi\n"
           << "# with higher statistics.\n"
           << "#\n"
           << "# Usage: bash run_g4_all.sh\n"
           << "#\n"
           << "# Requires: ceelo_g4val built and on PATH, or point G4VAL at it:\n"
           << "#   G4VAL=/path/to/ceelo_g4val bash run_g4_all.sh\n"
           << "\n"
           << "set -e\n"
           << "\n"
           << "G4VAL=\"${G4VAL:-ceelo_g4val}\"\n"
           << "N_CONE_HISTO=4000000       # events for cone-biased histogram runs\n"
           << "N_CONE_EFF=2000000         # events for cone-biased efficiency runs\n"
           << "N_ISO_HISTO=100000000      # events for isotropic histogram runs\n"
           << "N_ISO_EFF=50000000         # events for isotropic efficiency runs\n"
           << "cd \"$(dirname \"$0\")\"\n"
           << "\n"
           << "echo \"=== GEANT4 Reference Runs ===\"\n"
           << "echo \"Started: $(date)\"\n"
           << "echo \"\"\n"
           << "\n"
           << "# Helper: create a temporary macro with specified event count\n"
           << "make_mac() {\n"
           << "    local src_mac=$1\n"
           << "    local n_events=$2\n"
           << "    local tmp_mac=\"tmp_${src_mac}\"\n"
           << "    if [ ! -f \"$src_mac\" ]; then\n"
           << "        echo \"  SKIP: $src_mac not found\"\n"
           << "        return 1\n"
           << "    fi\n"
           << "    sed \"s|/run/beamOn.*|/run/beamOn $n_events|\" \"$src_mac\" > \"$tmp_mac\"\n"
           << "    echo \"$tmp_mac\"\n"
           << "}\n\n";

        // Histogram runs descriptor
        struct G4Run {
            int cfg;
            std::vector<double> energies;
            bool cone_bias;
        };
        std::vector<G4Run> g4runs = {
            {1, {59, 122, 185, 375, 662, 1173, 2614},  true},
            {2, {59, 122, 185, 375, 662, 1173, 2614},  true},
            {3, {59, 122, 185, 375, 662, 1173, 2614},  true},
            {5, {59, 122, 185, 375, 662},               true},
            {6, {122, 375, 662, 1173},                  true},
            {7, {185, 375, 662, 1173},                  true},
            {8, {59, 200, 662, 1173, 2614},             false},
            {9, {59, 200, 662, 1173, 2614},             false},
            {10, {59, 200, 662, 1173, 2614},            false},
            {11, {200, 375, 662, 1173},                 false},
            {12, {59, 200, 662, 1173},                  false},
        };
        // Apply cone_bias settings from the configs we just ran
        // (future configs with source shielding or extended sources will set false)
        for (auto& run : g4runs) {
            run.cone_bias = cone_bias_configs.count(run.cfg) > 0;
        }

        sh << "echo \"=== Phase 1: Histogram spectrum runs ===\"\n\n";
        for (const auto& run : g4runs) {
            std::string cfg_str = std::to_string(run.cfg);
            std::string n_var = run.cone_bias ? "$N_CONE_HISTO" : "$N_ISO_HISTO";
            std::string bias_str = run.cone_bias ? " --cone-bias" : "";
            std::string bias_label = run.cone_bias ? "cone-biased" : "isotropic";

            sh << "echo \"--- Config " << run.cfg << " (" << bias_label << ") ---\"\n";
            for (double E : run.energies) {
                int e_int = static_cast<int>(E);
                std::string mac = "run_" + cfg_str + "_" + std::to_string(e_int) + "keV.mac";
                std::string gdml = "detector_" + cfg_str + ".gdml";
                std::string out = "g4_cfg" + cfg_str + "_" + std::to_string(e_int) + "keV.csv";

                sh << "MAC=$(make_mac " << mac << " " << n_var << ") && \\\n"
                   << "  echo \"  [$(date +%H:%M:%S)] Config " << run.cfg
                   << ", " << e_int << " keV, " << n_var << " events...\" && \\\n"
                   << "  $G4VAL " << gdml << " \"$MAC\" " << out
                   << " --histogram" << bias_str << " && \\\n"
                   << "  rm -f \"$MAC\" && \\\n"
                   << "  echo \"  Done: " << out << "\"\n";
            }
            sh << "\n";
        }

        sh << "echo \"Phase 1 complete: $(date)\"\n"
           << "echo \"\"\n\n";

        // Multi-energy efficiency runs
        struct MultiRun {
            int cfg;
            std::vector<double> energies;
            bool cone_bias;
        };
        std::vector<MultiRun> multi_runs = {
            {1, {30, 50, 80, 100, 200, 300, 500, 662, 800, 1000, 1173, 1332, 2000, 3000}, true},
            {2, {30, 50, 80, 100, 200, 300, 500, 662, 800, 1000, 1173, 1332, 2000, 3000}, true},
            {3, {30, 50, 80, 100, 200, 300, 500, 662, 800, 1000, 1173, 1332, 2000, 3000}, true},
            {5, {30, 50, 80, 100, 200, 300, 500, 662, 800, 1000, 1500}, true},
            {6, {100, 200, 300, 500, 662, 1000, 1173, 1332, 2000}, true},
            {7, {100, 200, 300, 500, 662, 1000, 1173, 1332, 2000, 3000}, true},
            {8, {59, 100, 200, 300, 500, 662, 1000, 1173, 1332, 2614}, false},
        };
        for (auto& mr : multi_runs) {
            mr.cone_bias = cone_bias_configs.count(mr.cfg) > 0;
        }

        sh << "echo \"=== Phase 2: Multi-energy efficiency runs ===\"\n\n";
        for (const auto& mr : multi_runs) {
            std::string cfg_str = std::to_string(mr.cfg);
            std::string n_var = mr.cone_bias ? "$N_CONE_EFF" : "$N_ISO_EFF";
            std::string bias_str = mr.cone_bias ? " --cone-bias" : "";
            std::string bias_label = mr.cone_bias ? "cone-biased" : "isotropic";

            sh << "echo \"--- Config " << mr.cfg << " (" << bias_label << ") ---\"\n";
            for (double E : mr.energies) {
                int e_int = static_cast<int>(E);
                std::string mac = "run_" + cfg_str + "_" + std::to_string(e_int) + "keV.mac";
                std::string gdml = "detector_" + cfg_str + ".gdml";
                std::string out = "g4_cfg" + cfg_str + "_" + std::to_string(e_int) + "keV.csv";

                sh << "MAC=$(make_mac " << mac << " " << n_var << ") && \\\n"
                   << "  $G4VAL " << gdml << " \"$MAC\" " << out
                   << bias_str << " && \\\n"
                   << "  rm -f \"$MAC\"\n";
            }

            // Combine into multi CSV
            std::string multi_out = "g4_cfg" + cfg_str + "_multi.csv";
            sh << "echo \"energy_keV,fep_efficiency,fep_uncertainty,"
               << "total_efficiency,total_uncertainty,num_events\" > " << multi_out << "\n";
            for (double E : mr.energies) {
                int e_int = static_cast<int>(E);
                std::string eff_csv = "g4_cfg" + cfg_str + "_" + std::to_string(e_int) + "keV.csv";
                sh << "tail -1 " << eff_csv << " >> " << multi_out << "\n";
            }
            sh << "echo \"  Combined: " << multi_out << "\"\n\n";
        }

        sh << "echo \"\"\n"
           << "echo \"=== All GEANT4 runs complete: $(date) ===\"\n";
        std::cout << "\n  G4 run script: run_g4_all.sh\n";
    }

    std::cout << "\n=== All MC spectra and G4 files generated ===\n";
    return 0;
}
