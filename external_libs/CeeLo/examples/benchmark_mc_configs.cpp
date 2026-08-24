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

#include "bench_configs.h"

#include "efficiency/EfficiencyCalculator.h"
#include "materials/Material.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace ceelo;

namespace {

constexpr uint64_t kDefaultEvents = 1'000'000;
constexpr unsigned kDefaultThreads = 0;
constexpr uint64_t kDefaultMacroEvents = 500'000;
constexpr uint64_t kPrecisionMaxEvents = 200'000'000;
constexpr double kPrecisionMaxSeconds = 600.0;

const char* stop_reason_name(StopReason r) {
    switch (r) {
        case StopReason::MaxEvents:      return "max_events";
        case StopReason::FepPrecision:   return "fep_precision";
        case StopReason::TotalPrecision: return "total_precision";
        case StopReason::WallTime:       return "wall_time";
    }
    return "unknown";
}

void print_usage(const char* argv0) {
    std::cout
        << "Usage: " << argv0 << " [options]\n\n"
        << "Options:\n"
        << "  --config <list>          Comma-separated config IDs to run\n"
        << "                           (default: 1,2,3,5,6,7; also supports 8,11,12)\n"
        << "  --events <N>             Events per energy for fixed-N runs (default: 1000000)\n"
        << "  --precision <p>          Run to relative precision p on FEP and total\n"
        << "                           (e.g. 0.005); overrides --events\n"
        << "  --threads <N>            Number of threads (0 = auto, default: 0)\n"
        << "  --energies <list>        Override energy grid (keV) for all selected configs\n"
        << "  --fom-csv <file>         Append per-run FOM rows to CSV file\n"
        << "  --label <name>           Mode label written to FOM CSV rows (default: baseline)\n"
        << "  --force-interaction      Enable forced first interaction in detector\n"
        << "  --fep-only               Run in FEP-only mode (no total efficiency)\n"
        << "  --mixture-alpha <a>      Mixture angular biasing alpha for source-effect\n"
        << "                           configs (0 = disabled)\n"
        << "  --two-stream <f>         Two-stream stratified estimator for source-effect\n"
        << "                           configs; f = direct-stream fraction in [0.05,0.95]\n"
        << "  --compton-bias <g>       Compton-angle mixture biasing gamma in [0,0.9]\n"
        << "                           (source-effect configs; combine with --two-stream)\n"
        << "  (source-escape model is COMPILE-TIME: -DCEELO_SOURCE_ESCAPE_MODEL=gate|tails|gs)\n"
        << "  --no-source-brems        Disable bremsstrahlung from electrons in source\n"
        << "  --no-source-electrons    Disable source-electron escape transport\n"
           "                           material/shielding (A/B quantification)\n"
        << "  --analog                 Force fully analog sampling (disable the\n"
        << "                           library's auto-enabled biasing)\n"
        << "  --histogram              Also write pulse-height spectra (1 keV bins)\n"
        << "                           to our_<cfg>_<E>keV_spec.csv (fixed-N mode)\n"
        << "  (default: library auto-enables biasing per geometry/energy)\n"
        << "  --skip-export            Skip GDML/macro export (recommended for profiling)\n"
        << "  --macro-events <N>       Events written into exported macros (default: 500000)\n"
        << "  --help                   Show this help\n";
}

bool parse_int_list(const std::string& text, std::vector<int>& out) {
    std::stringstream ss(text);
    std::string token;
    while (std::getline(ss, token, ',')) {
        if (token.empty()) {
            return false;
        }
        try {
            out.push_back(std::stoi(token));
        } catch (...) {
            return false;
        }
    }
    return !out.empty();
}

bool parse_double_list(const std::string& text, std::vector<double>& out) {
    std::stringstream ss(text);
    std::string token;
    while (std::getline(ss, token, ',')) {
        if (token.empty()) {
            return false;
        }
        try {
            out.push_back(std::stod(token));
        } catch (...) {
            return false;
        }
    }
    return !out.empty();
}

void write_our_csv_batch(const std::string& filename,
                         const std::vector<double>& energies_keV,
                         const std::vector<EfficiencyResult>& results) {
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

struct FomWriter {
    std::ofstream file;
    bool enabled = false;

    void open(const std::string& path, const std::string& /*label*/) {
        // Append so baseline and per-technique runs can share one file.
        bool write_header = true;
        {
            std::ifstream existing(path);
            write_header = !existing.good() || existing.peek() == EOF;
        }
        file.open(path, std::ios::app);
        enabled = file.is_open();
        if (enabled && write_header) {
            file << "config,energy_keV,mode,num_events,wall_s,"
                    "fep,fep_rel_unc,tot,tot_rel_unc,"
                    "fom_fep,fom_tot,stop_reason\n";
        }
    }

    void add_row(int cfg, double energy, const std::string& mode,
                 const EfficiencyResult& r) {
        if (!enabled) return;
        double fep_rel = (r.full_energy_peak_efficiency > 0.0)
            ? r.fep_uncertainty / r.full_energy_peak_efficiency : 0.0;
        double tot_rel = (r.total_efficiency > 0.0)
            ? r.total_uncertainty / r.total_efficiency : 0.0;
        double fom_fep = (fep_rel > 0.0 && r.wall_time_seconds > 0.0)
            ? 1.0 / (fep_rel * fep_rel * r.wall_time_seconds) : 0.0;
        double fom_tot = (tot_rel > 0.0 && r.wall_time_seconds > 0.0)
            ? 1.0 / (tot_rel * tot_rel * r.wall_time_seconds) : 0.0;
        file << cfg << ","
             << std::fixed << std::setprecision(1) << energy << ","
             << mode << ","
             << r.num_events_simulated << ","
             << std::fixed << std::setprecision(4) << r.wall_time_seconds << ","
             << std::scientific << std::setprecision(6)
             << r.full_energy_peak_efficiency << "," << fep_rel << ","
             << r.total_efficiency << "," << tot_rel << ","
             << fom_fep << "," << fom_tot << ","
             << stop_reason_name(r.stop_reason) << "\n";
        file.flush();
    }
};

bool run_config(int cfg,
                const std::vector<double>& energies,
                uint64_t n_mc,
                double precision,
                unsigned n_threads,
                bool skip_export,
                uint64_t macro_events,
                const std::string& label,
                const BiasingConfig* biasing,
                bool fep_only,
                bool histogram,
                bool no_source_brems,
                bool no_source_electrons,
                FomWriter& fom) {
    EfficiencyCalculator calc;
    bench::ConfigSetup setup;
    if (!bench::make_config(cfg, calc, setup)) {
        std::cerr << "Unsupported config: " << cfg << "\n";
        return false;
    }
    if (biasing) calc.set_biasing(*biasing);  // null = library auto-enable
    if (fep_only) calc.enable_fep_only_mode(true);
    if (no_source_brems) calc.set_source_brems(false);
    if (no_source_electrons) calc.enable_source_electron_transport(false);

    std::cout << "\n=== Config " << cfg << ": " << setup.description << " ===\n";

    // GDML/macro export only supported for the point-source detector-only
    // configs (source-effect configs are exported by generate_all_spectra).
    // cfg 20/21 (high-Z skin-escape validation, box + Pb/W wall) opt in here so
    // their GEANT4 inputs can be generated directly.
    // cfg 25/26 (sharp vs bulletized HPGe) are point-source detector-only
    // configs, so they export directly like configs 1-7.
    const bool can_export = (cfg <= 7) || cfg == 12 || cfg == 20 ||
                            cfg == 21 || cfg == 22 || cfg == 23 || cfg == 24 ||
                            cfg == 25 || cfg == 26;
    if (!skip_export && can_export) {
        // vacuum_world=true matches the validated generate_all_spectra
        // reference export (air in the source-detector gap would otherwise
        // perturb the comparison).
        calc.export_geant4_gdml("detector_" + std::to_string(cfg) + ".gdml",
                                /*vacuum_world=*/true);
    }

    std::vector<EfficiencyResult> results;
    if (precision > 0.0) {
        for (double e : energies) {
            SimulationConfig sim;
            sim.energy_keV = e;
            sim.num_threads = n_threads;
            sim.termination.target_fep_rel_precision = precision;
            sim.termination.target_total_rel_precision = precision;
            sim.termination.max_events = kPrecisionMaxEvents;
            sim.termination.max_wall_seconds = kPrecisionMaxSeconds;
            results.push_back(calc.compute(sim));
        }
        write_our_csv_batch("our_" + std::to_string(cfg) + "_multi.csv",
                            energies, results);
    } else if (histogram) {
        for (double e : energies) {
            std::vector<float> edges;
            for (float b = 0.0f; b <= static_cast<float>(e) + 10.0f; b += 1.0f) {
                edges.push_back(b);
            }
            auto r = calc.compute(e, n_mc, n_threads, edges);
            std::ostringstream name;
            name << "our_" << cfg << "_" << static_cast<int>(e)
                 << "keV_" << label << "_spec.csv";
            std::ofstream sf(name.str());
            sf << "# MC pulse-height spectrum (probability per 1 keV bin)\n"
               << "# Config " << cfg << ": " << setup.description << "\n"
               << "# Primary energy: " << e << " keV\n"
               << "# Events: " << r.num_events_simulated << "\n"
               << "Energy (keV),Probability\n";
            for (size_t i = 0; i < r.pulse_height_distribution.size(); ++i) {
                sf << edges[i] << "," << std::scientific << std::setprecision(6)
                   << r.pulse_height_distribution[i] << "\n";
            }
            results.push_back(std::move(r));
        }
        write_our_csv_batch("our_" + std::to_string(cfg) + "_multi.csv",
                            energies, results);
    } else {
        results = calc.compute_batch(energies, n_mc, n_threads);
        write_our_csv_batch("our_" + std::to_string(cfg) + "_multi.csv",
                            energies, results);
    }

    if (!skip_export && can_export) {
        for (double e : energies) {
            std::string mac = "run_" + std::to_string(cfg) + "_" +
                std::to_string(static_cast<int>(e)) + "keV.mac";
            calc.export_geant4_macro(mac, e, macro_events);
        }
    }

    for (size_t i = 0; i < energies.size(); ++i) {
        const auto& r = results[i];
        double fep_rel = (r.full_energy_peak_efficiency > 0.0)
            ? r.fep_uncertainty / r.full_energy_peak_efficiency : 0.0;
        std::cout << "  E=" << std::fixed << std::setprecision(0)
                  << std::setw(6) << energies[i]
                  << " keV  t=" << std::setprecision(3)
                  << r.wall_time_seconds << " s"
                  << "  N=" << std::setw(10) << r.num_events_simulated
                  << "  FEP=" << std::scientific << std::setprecision(6)
                  << r.full_energy_peak_efficiency
                  << " (" << std::fixed << std::setprecision(3)
                  << 100.0 * fep_rel << "%)"
                  << "  TOT=" << std::scientific << std::setprecision(6)
                  << r.total_efficiency
                  << "  [" << stop_reason_name(r.stop_reason) << "]\n";
        // Source-effect variance decomposition: share of efficiency and of
        // the estimator's second moment (~variance for small eps) carried by
        // the scattered-primary class.
        const auto& sd = r.src_diag;
        if (sd.n_u + sd.n_s > 0) {
            const double N = static_cast<double>(r.num_events_simulated);
            const double tot_u = sd.any_w_u / N, tot_s = sd.any_w_s / N;
            const double fep_u = sd.fep_w_u / N, fep_s = sd.fep_w_s / N;
            auto share = [](double a, double b) {
                return (a + b > 0.0) ? 100.0 * b / (a + b) : 0.0;
            };
            std::cout << "         src-classes: scattered "
                      << std::fixed << std::setprecision(1)
                      << share(static_cast<double>(sd.n_u),
                               static_cast<double>(sd.n_s))
                      << "% of events | TOT: eps_u=" << std::scientific
                      << std::setprecision(3) << tot_u
                      << " eps_s=" << tot_s
                      << " (s-share " << std::fixed << std::setprecision(1)
                      << share(tot_u, tot_s) << "% eff, "
                      << share(sd.any_w2_u, sd.any_w2_s) << "% w2)"
                      << " | FEP: eps_u=" << std::scientific
                      << std::setprecision(3) << fep_u
                      << " eps_s=" << fep_s
                      << " (s-share " << std::fixed << std::setprecision(1)
                      << share(fep_u, fep_s) << "% eff, "
                      << share(sd.fep_w2_u, sd.fep_w2_s) << "% w2)\n";
            // Absorbed-primary escaped-electron channel (dropped entirely
            // before June 2026; a subset of the s-class).
            if (sd.n_e_only > 0) {
                const double eps_e = sd.e_only_any_w / N;
                const double sig_e = std::sqrt(sd.e_only_any_w2) / N;
                const double tot = r.total_efficiency;
                std::cout << "         e-only channel: "
                          << sd.n_e_only << " events ("
                          << std::fixed << std::setprecision(2)
                          << static_cast<double>(sd.n_e_only_electrons) /
                                 sd.n_e_only
                          << " e/evt, mean escape KE "
                          << std::setprecision(0)
                          << sd.e_only_energy_keV / sd.n_e_only_electrons
                          << " keV) | eps_e=" << std::scientific
                          << std::setprecision(3) << eps_e
                          << " +- " << sig_e
                          << " (" << std::fixed << std::setprecision(3)
                          << (tot > 0.0 ? 100.0 * eps_e / tot : 0.0)
                          << "% of TOT)  fep_w=" << std::scientific
                          << std::setprecision(2) << sd.e_only_fep_w << "\n";
            }
        }
        if (r.pp_secondaries_processed > 0) {
            std::cout << "         source-PP gammas: "
                      << r.pp_secondaries_processed
                      << " (" << std::scientific << std::setprecision(3)
                      << static_cast<double>(r.pp_secondaries_processed)
                         / r.num_events_simulated
                      << "/event)  dep=" << r.pp_secondary_deposit_keV
                      << " keV  pp-only-any contrib="
                      << r.pp_only_any_weight / r.num_events_simulated
                      << "\n";
        }
        fom.add_row(cfg, energies[i], label, r);
    }
    return true;
}

} // namespace

int main(int argc, char** argv) {
    std::vector<int> configs = {1, 2, 3, 5, 6, 7};
    uint64_t n_mc = kDefaultEvents;
    double precision = 0.0;
    unsigned n_threads = kDefaultThreads;
    bool skip_export = false;
    uint64_t macro_events = kDefaultMacroEvents;
    std::vector<double> energy_override;
    std::string fom_csv;
    std::string label = "baseline";
    BiasingConfig biasing;
    bool biasing_explicit = false;
    bool fep_only = false;
    bool histogram = false;
    bool no_source_brems = false;
    bool no_source_electrons = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help") {
            print_usage(argv[0]);
            return 0;
        }
        if (arg == "--config" && i + 1 < argc) {
            std::vector<int> parsed;
            if (!parse_int_list(argv[++i], parsed)) {
                std::cerr << "Failed to parse --config list\n";
                return 1;
            }
            configs = parsed;
            continue;
        }
        if (arg == "--events" && i + 1 < argc) {
            n_mc = static_cast<uint64_t>(std::stoull(argv[++i]));
            continue;
        }
        if (arg == "--precision" && i + 1 < argc) {
            precision = std::stod(argv[++i]);
            if (precision <= 0.0 || precision >= 1.0) {
                std::cerr << "--precision must be in (0, 1)\n";
                return 1;
            }
            continue;
        }
        if (arg == "--threads" && i + 1 < argc) {
            n_threads = static_cast<unsigned>(std::stoul(argv[++i]));
            continue;
        }
        if (arg == "--energies" && i + 1 < argc) {
            if (!parse_double_list(argv[++i], energy_override)) {
                std::cerr << "Failed to parse --energies list\n";
                return 1;
            }
            continue;
        }
        if (arg == "--fom-csv" && i + 1 < argc) {
            fom_csv = argv[++i];
            continue;
        }
        if (arg == "--label" && i + 1 < argc) {
            label = argv[++i];
            continue;
        }
        if (arg == "--force-interaction") {
            biasing.force_detector_interaction = true;
            biasing_explicit = true;
            continue;
        }
        if (arg == "--analog") {
            biasing_explicit = true;
            continue;
        }
        if (arg == "--fep-only") {
            fep_only = true;
            continue;
        }
        if (arg == "--histogram") {
            histogram = true;
            continue;
        }
        if (arg == "--mixture-alpha" && i + 1 < argc) {
            biasing.mixture_cone_alpha = std::stod(argv[++i]);
            if (biasing.mixture_cone_alpha < 0.0 || biasing.mixture_cone_alpha > 1.0) {
                std::cerr << "--mixture-alpha must be in [0, 1]\n";
                return 1;
            }
            biasing_explicit = true;
            continue;
        }
        if (arg == "--two-stream" && i + 1 < argc) {
            biasing.two_stream = true;
            biasing.two_stream_direct_fraction = std::stod(argv[++i]);
            if (biasing.two_stream_direct_fraction < 0.05 ||
                biasing.two_stream_direct_fraction > 0.95) {
                std::cerr << "--two-stream fraction must be in [0.05, 0.95]\n";
                return 1;
            }
            biasing_explicit = true;
            continue;
        }
        if (arg == "--compton-bias" && i + 1 < argc) {
            biasing.compton_cone_fraction = std::stod(argv[++i]);
            if (biasing.compton_cone_fraction < 0.0 ||
                biasing.compton_cone_fraction > 0.9) {
                std::cerr << "--compton-bias must be in [0, 0.9]\n";
                return 1;
            }
            biasing_explicit = true;
            continue;
        }
        if (arg == "--no-source-electrons") {
            no_source_electrons = true;
            continue;
        }
        if (arg == "--no-source-brems") {
            no_source_brems = true;
            continue;
        }
        if (arg == "--skip-export") {
            skip_export = true;
            continue;
        }
        if (arg == "--macro-events" && i + 1 < argc) {
            macro_events = static_cast<uint64_t>(std::stoull(argv[++i]));
            continue;
        }
        std::cerr << "Unknown or incomplete option: " << arg << "\n";
        print_usage(argv[0]);
        return 1;
    }

    std::sort(configs.begin(), configs.end());
    configs.erase(std::unique(configs.begin(), configs.end()), configs.end());

    FomWriter fom;
    if (!fom_csv.empty()) {
        fom.open(fom_csv, label);
        if (!fom.enabled) {
            std::cerr << "Failed to open FOM CSV: " << fom_csv << "\n";
            return 1;
        }
    }

    std::cout << "CeeLo benchmark runner\n"
              << "  configs: ";
    for (size_t i = 0; i < configs.size(); ++i) {
        if (i) std::cout << ",";
        std::cout << configs[i];
    }
    std::cout << "\n  mode: ";
    if (precision > 0.0) {
        std::cout << "precision-targeted (" << precision * 100.0 << "% rel)";
    } else {
        std::cout << "fixed-N (" << n_mc << " events per energy)";
    }
    std::cout << "\n  label: " << label
              << "\n  threads: " << n_threads
              << "\n  skip export: " << (skip_export ? "true" : "false") << "\n";

    bool ok = true;
    for (int cfg : configs) {
        std::vector<double> energies =
            energy_override.empty() ? bench::default_energies_for_cfg(cfg)
                                    : energy_override;
        if (energies.empty()) {
            std::cerr << "No energies available for config " << cfg << "\n";
            ok = false;
            continue;
        }
        bool ran = run_config(cfg, energies, n_mc, precision, n_threads,
                              skip_export, macro_events, label,
                              biasing_explicit ? &biasing : nullptr,
                              fep_only, histogram, no_source_brems,
                              no_source_electrons, fom);
        ok = ok && ran;
    }

    std::cout << "\nDone.\n";
    return ok ? 0 : 1;
}
