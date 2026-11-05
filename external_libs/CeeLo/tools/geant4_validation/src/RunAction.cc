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

#include "RunAction.hh"

#include "G4Run.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <algorithm>

RunAction::RunAction(const std::string& output_csv,
                     const std::string& gdml_file,
                     const std::string& macro_file)
    : output_csv_(output_csv)
    , gdml_file_(gdml_file)
    , macro_file_(macro_file)
{}

void RunAction::BeginOfRunAction(const G4Run* run) {
    if (preserve_across_runs_) {
        // Accumulate events across sub-runs
        n_total_ += static_cast<uint64_t>(run->GetNumberOfEventToBeProcessed());
        return;
    }

    n_total_ = static_cast<uint64_t>(run->GetNumberOfEventToBeProcessed());
    w_fep_ = 0.0;
    w2_fep_ = 0.0;
    w_any_ = 0.0;
    w2_any_ = 0.0;
    sum_weights_ = 0.0;

    // Reset the sphere-flux sufficient statistics together.  In particular,
    // the scalar moments below are not reconstructed from cell uncertainties:
    // doing that would discard cross-cell covariance from one primary event.
    if (sphere_flux_enabled_) {
        std::fill(sphere_direct_.begin(), sphere_direct_.end(), 0.0);
        std::fill(sphere_scat_.begin(), sphere_scat_.end(), 0.0);
        std::fill(sphere_scat_ev_.begin(), sphere_scat_ev_.end(), 0.0);
        std::fill(sphere_scat_w2_.begin(), sphere_scat_w2_.end(), 0.0);
        std::fill(sphere_origin_.begin(), sphere_origin_.end(), 0.0);
        sphere_scat_total_n_ = 0;
        sphere_scat_total_s1_cm_ = 0.0L;
        sphere_scat_total_s2_cm2_ = 0.0L;
    }

    // Reset histogram data
    if (histogram_enabled_) {
        histogram_sumw_.assign(num_bins_, 0.0);
        histogram_sumw2_.assign(num_bins_, 0.0);
        if (component_labels_enabled_) {
            for (auto& values : component_histogram_sumw_)
                values.assign(num_bins_, 0.0);
            for (auto& values : component_histogram_sumw2_)
                values.assign(num_bins_, 0.0);
        }
        sum_electron_edep_MeV_ = 0.0;
        sum_gamma_edep_MeV_ = 0.0;
        sum_total_edep_MeV_ = 0.0;
        sum_deposit_weights_ = 0.0;
        n_events_with_deposit_ = 0;
    }

    // Reset fluence-at-detector scorer
    if (fluence_enabled_) {
        fluence_w_.assign(num_bins_, 0.0);
        fluence_w2_.assign(num_bins_, 0.0);
        unc_w_ = unc_w2_ = 0.0;
        scat_w_ = scat_w2_ = 0.0;
    }

    // Reset crystal entry diagnostics
    if (entry_diag_enabled_) {
        entry_energy_face_.assign(kEntryEnergyBins, 0);
        entry_energy_side_.assign(kEntryEnergyBins, 0);
        entry_z_.assign(kEntryZBins, 0);
        entry_face_count_ = 0;
        entry_side_count_ = 0;
        entry_back_count_ = 0;
        entry_face_energy_sum_ = 0.0;
        entry_side_energy_sum_ = 0.0;
        entry_back_energy_sum_ = 0.0;
        entry_total_ = 0;
    }

    // Reset external electron counters
    n_ext_electron_events_ = 0;
    n_ext_electron_with_deposit_ = 0;
    n_ext_electron_only_deposit_ = 0;
    n_src_u_ = 0; n_src_s_ = 0;
    src_any_w_u_ = 0.0; src_any_w2_u_ = 0.0;
    src_any_w_s_ = 0.0; src_any_w2_s_ = 0.0;
    src_fep_w_u_ = 0.0; src_fep_w2_u_ = 0.0;
    src_fep_w_s_ = 0.0; src_fep_w2_s_ = 0.0;

    // Reset diagnostic histograms
    if (diagnostics_enabled_) {
        diag_cs1_depth_.assign(kDiagDepthBins, 0);
        diag_cs1_angle_.assign(kDiagAngleBins, 0);
        diag_cs1_energy_.assign(kDiagEnergyBins, 0);
        diag_cs2_depth_.assign(kDiagDepthBins, 0);
        diag_cs2_angle_.assign(kDiagAngleBins, 0);
        diag_cs2_energy_.assign(kDiagEnergyBins, 0);
        diag_pe_depth_.assign(kDiagDepthBins, 0);
        diag_pe_energy_.assign(kDiagEnergyBins, 0);
        diag_pe_nscatter_.assign(kDiagScatterBins, 0);
        diag_scatter_mult_.assign(kDiagScatterBins, 0);
        diag_escape_energy_.assign(kDiagEnergyBins, 0);
        diag_fluor_escape_.assign(kDiagEnergyBins, 0);
        n_fluor_escaped_total_ = 0;
    }
}

void RunAction::RecordHistogram(double total_edep_MeV,
                                double electron_edep_MeV,
                                double gamma_edep_MeV,
                                double weight,
                                const std::vector<G4InteractionRecord>& interactions,
                                double primary_energy_MeV,
                                bool had_pair_production,
                                int pair_secondary_compton_count) {
    std::lock_guard<std::mutex> lk(mutex_);

    const int bin = std::min(
        num_bins_ - 1,
        std::max(0, static_cast<int>(total_edep_MeV / bin_width_MeV_)));
    histogram_sumw_[bin] += weight;
    histogram_sumw2_[bin] += weight * weight;

    if (component_labels_enabled_) {
        // zero, non-pair FEP, pair FEP, single escape, double escape,
        // annihilation-511, photoelectron escape, one Compton,
        // multiple Compton, pair+zero/one/multiple subsequent Comptons
        int component = 0;
        constexpr double line_tol_MeV = 0.0015;
        if (total_edep_MeV > 0.0) {
            int compton = 0;
            bool pair = had_pair_production;
            bool photoelectric = false;
            for (const auto& interaction : interactions) {
                if (interaction.type == G4InteractionRecord::kCompton)
                    ++compton;
                if (interaction.type
                    == G4InteractionRecord::kPhotoelectric)
                    photoelectric = true;
                if (interaction.type
                    == G4InteractionRecord::kPairProduction)
                    pair = true;
            }
            if (std::abs(total_edep_MeV - primary_energy_MeV)
                < line_tol_MeV) {
                component = pair ? 2 : 1;
            } else if (pair && primary_energy_MeV > 1.0219979
                       && std::abs(total_edep_MeV
                                   - (primary_energy_MeV - 0.51099895))
                              < line_tol_MeV) {
                component = 3;
            } else if (pair && primary_energy_MeV > 1.0219979
                       && std::abs(total_edep_MeV
                                   - (primary_energy_MeV - 1.0219979))
                              < line_tol_MeV) {
                component = 4;
            } else if (
                pair
                && std::abs(total_edep_MeV - 0.51099895)
                       < line_tol_MeV) {
                component = 5;
            } else {
                if (pair) {
                    component = pair_secondary_compton_count == 0
                        ? 9
                        : pair_secondary_compton_count == 1 ? 10 : 11;
                }
                else if (photoelectric && compton == 0)
                    component = 6;
                else
                    component = compton <= 1 ? 7 : 8;
            }
        }
        component_histogram_sumw_[component][bin] += weight;
        component_histogram_sumw2_[component][bin] += weight * weight;
    }

    if (total_edep_MeV > 0.0) {
        sum_total_edep_MeV_ += total_edep_MeV;
        sum_electron_edep_MeV_ += weight * electron_edep_MeV;
        sum_gamma_edep_MeV_ += weight * gamma_edep_MeV;
        sum_total_edep_MeV_ += (weight - 1.0) * total_edep_MeV;
        sum_deposit_weights_ += weight;
        n_events_with_deposit_++;
    }
}

void RunAction::RecordDiagnostics(
    const std::vector<G4InteractionRecord>& interactions,
    double edep_MeV,
    double primary_energy_MeV)
{
    std::lock_guard<std::mutex> lk(mutex_);

    int n_compton_in_crystal = 0;
    bool had_pe = false;

    for (const auto& rec : interactions) {
        // Depth: z position in mm -> cm. Detector front face at z=0.
        double depth_cm = rec.position.z() / 10.0;  // mm -> cm

        // Energy in keV
        double e_before_keV = rec.energy_before_MeV * 1000.0;
        double e_after_keV  = rec.energy_after_MeV * 1000.0;

        if (rec.type == G4InteractionRecord::kCompton) {
            if (rec.compton_index == 0) {
                // 1st Compton scatter
                int bd = static_cast<int>(depth_cm / kDiagDepthWidth);
                if (bd >= 0 && bd < kDiagDepthBins) diag_cs1_depth_[bd]++;

                double angle_deg = std::acos(std::clamp(rec.cos_theta, -1.0, 1.0))
                                   * 180.0 / M_PI;
                int ba = static_cast<int>(angle_deg / kDiagAngleWidth);
                if (ba >= 0 && ba < kDiagAngleBins) diag_cs1_angle_[ba]++;

                int be = static_cast<int>(e_after_keV / kDiagEnergyWidth);
                if (be >= 0 && be < kDiagEnergyBins) diag_cs1_energy_[be]++;

            } else if (rec.compton_index == 1) {
                // 2nd Compton scatter
                int bd = static_cast<int>(depth_cm / kDiagDepthWidth);
                if (bd >= 0 && bd < kDiagDepthBins) diag_cs2_depth_[bd]++;

                double angle_deg = std::acos(std::clamp(rec.cos_theta, -1.0, 1.0))
                                   * 180.0 / M_PI;
                int ba = static_cast<int>(angle_deg / kDiagAngleWidth);
                if (ba >= 0 && ba < kDiagAngleBins) diag_cs2_angle_[ba]++;

                int be = static_cast<int>(e_after_keV / kDiagEnergyWidth);
                if (be >= 0 && be < kDiagEnergyBins) diag_cs2_energy_[be]++;
            }
            n_compton_in_crystal++;

        } else if (rec.type == G4InteractionRecord::kPhotoelectric) {
            int bd = static_cast<int>(depth_cm / kDiagDepthWidth);
            if (bd >= 0 && bd < kDiagDepthBins) diag_pe_depth_[bd]++;

            int be = static_cast<int>(e_before_keV / kDiagEnergyWidth);
            if (be >= 0 && be < kDiagEnergyBins) diag_pe_energy_[be]++;

            int bs = std::min(rec.compton_index, kDiagScatterBins - 1);
            if (bs >= 0) diag_pe_nscatter_[bs]++;
            had_pe = true;
        }
    }

    // Scatter multiplicity: count Compton scatters for events with any deposit
    if (edep_MeV > 0.0) {
        int mult = std::min(n_compton_in_crystal, kDiagScatterBins - 1);
        diag_scatter_mult_[mult]++;
    }

    // Escape energy: primary gamma energy minus total deposit gives escape.
    // But more precisely, if the primary gamma escaped the crystal (wasn't
    // fully absorbed), its last energy is what escaped. For simplicity,
    // we look at whether there was a PE: if not, the last Compton scatter's
    // after-energy approximates the escape energy.
    if (edep_MeV > 0.0 && !had_pe && !interactions.empty()) {
        // Find the last Compton interaction
        for (int i = static_cast<int>(interactions.size()) - 1; i >= 0; --i) {
            if (interactions[i].type == G4InteractionRecord::kCompton) {
                double esc_keV = interactions[i].energy_after_MeV * 1000.0;
                int be = static_cast<int>(esc_keV / kDiagEnergyWidth);
                if (be >= 0 && be < kDiagEnergyBins) diag_escape_energy_[be]++;
                break;
            }
        }
    }
}

void RunAction::RecordFluorEscape(double escape_energy_keV) {
    std::lock_guard<std::mutex> lk(mutex_);
    int be = static_cast<int>(escape_energy_keV / kDiagEnergyWidth);
    if (be >= 0 && be < kDiagEnergyBins) diag_fluor_escape_[be]++;
    n_fluor_escaped_total_++;
}

void RunAction::EndOfRunAction(const G4Run* run) {
    // In precision-driven mode (PreserveAcrossRuns) the tallies accumulate over
    // many small beamOn batches, so the estimator denominator must be the running
    // total -- NOT this batch's event count, which would inflate the efficiency by
    // (total / batch) and report a zero uncertainty.
    uint64_t n_run = preserve_across_runs_
                         ? n_total_
                         : static_cast<uint64_t>(run->GetNumberOfEvent());
    if (n_run == 0) return;

    double N = static_cast<double>(n_run);

    // Efficiency = sum_of_weights / N
    // For unweighted (GPS) mode: weight=1 per event, so this = n_hits/N
    // For cone-bias mode: weight = cone_solid_angle_fraction per hit
    double eps_fep   = w_fep_ / N;
    double eps_total = w_any_ / N;

    // Uncertainty: sigma = sqrt( (sum_w2/N - (sum_w/N)^2) / N )
    // = sqrt( (sum_w2 - sum_w^2/N) / N^2 )
    double var_fep   = (w2_fep_ / N - eps_fep * eps_fep) / N;
    double var_total = (w2_any_ / N - eps_total * eps_total) / N;
    double sig_fep   = std::sqrt(std::max(0.0, var_fep));
    double sig_total = std::sqrt(std::max(0.0, var_total));

    // Read primary energy from the GPS macro.  Environment-export macros use
    // /gps/energy while some detector macros use the equivalent
    // /gps/ene/mono spelling.
    double energy_keV = 0.0;
    {
        std::ifstream mac(macro_file_);
        std::string line;
        while (std::getline(mac, line)) {
            std::istringstream iss(line);
            std::string cmd;
            double val = 0.0;
            std::string unit;
            if (!(iss >> cmd >> val >> unit)) continue;
            if (cmd != "/gps/ene/mono" && cmd != "/gps/energy") continue;
            if (unit == "MeV") energy_keV = val * 1000.0;
            else if (unit == "eV") energy_keV = val * 0.001;
            else energy_keV = val;
            break;
        }
    }

    std::ofstream csv(output_csv_);
    if (!csv) {
        std::cerr << "RunAction: cannot open output file: " << output_csv_ << "\n";
        return;
    }

    csv << "# ceelo_g4val reference data\n"
        << "# geometry: " << gdml_file_ << "\n"
        << "# macro: "    << macro_file_ << "\n"
        << "energy_keV,fep_efficiency,fep_uncertainty,"
           "total_efficiency,total_uncertainty,num_events\n"
        << energy_keV << ","
        << eps_fep    << "," << sig_fep   << ","
        << eps_total  << "," << sig_total << ","
        << n_run << "\n";

    std::cout << "\n=== ceelo_g4val results ===\n"
              << "  Energy:          " << energy_keV << " keV\n"
              << "  Events:          " << n_run << "\n"
              << "  FEP efficiency:  " << eps_fep   << " +/- " << sig_fep   << "\n"
              << "  Total efficiency:" << eps_total << " +/- " << sig_total << "\n"
              << "  Output:          " << output_csv_ << "\n";

    // u/s source-class decomposition (printed whenever any event was
    // classified s — i.e. there is material outside the crystal).
    if (n_src_s_ > 0) {
        auto eps = [N](double w) { return w / N; };
        auto sig = [N](double w2) { return std::sqrt(w2) / N; };
        std::cout << "\n=== Source-class decomposition (u/s) ===\n"
                  << "  events: u=" << n_src_u_ << "  s=" << n_src_s_ << "\n"
                  << std::scientific << std::setprecision(5)
                  << "  TOT: eps_u=" << eps(src_any_w_u_)
                  << " +/- " << sig(src_any_w2_u_)
                  << "  eps_s=" << eps(src_any_w_s_)
                  << " +/- " << sig(src_any_w2_s_) << "\n"
                  << "  FEP: eps_u=" << eps(src_fep_w_u_)
                  << " +/- " << sig(src_fep_w2_u_)
                  << "  eps_s=" << eps(src_fep_w_s_)
                  << " +/- " << sig(src_fep_w2_s_) << "\n";
    }

    // Print external electron entry statistics (always, if any were recorded)
    if (n_ext_electron_events_ > 0) {
        double frac_events = 100.0 * n_ext_electron_events_ / N;
        double frac_deposit = 100.0 * n_ext_electron_with_deposit_ / N;
        double frac_only = 100.0 * n_ext_electron_only_deposit_ / N;
        std::cout << "\n=== External Electron Entry ===\n"
                  << "  Events with e-/e+ entering crystal: "
                  << n_ext_electron_events_ << " (" << std::fixed
                  << std::setprecision(4) << frac_events << "% of total)\n"
                  << "  ... with energy deposit: "
                  << n_ext_electron_with_deposit_ << " (" << frac_deposit << "%)\n"
                  << "  ... deposit but NO gamma entered: "
                  << n_ext_electron_only_deposit_ << " (" << frac_only << "%)\n";
    }

    // Write histogram and summary if enabled
    if (histogram_enabled_) {
        WriteHistogramCSV(energy_keV, n_run);
        WriteSummaryCSV(energy_keV, n_run);
    }

    // Write diagnostic CSVs if enabled
    if (diagnostics_enabled_) {
        WriteDiagnosticCSVs(energy_keV, n_run);
    }

    // Write crystal entry diagnostic CSVs if enabled
    if (entry_diag_enabled_) {
        WriteEntryDiagCSVs(energy_keV, n_run);
    }

    // Write fluence-at-detector CSV if enabled
    if (fluence_enabled_) {
        WriteFluenceCSV(energy_keV, n_run);
    }
    if (sphere_flux_enabled_) {
        SetSphereEnergy(energy_keV);
        WriteSphereFluxCSV(n_run);
    }
}

void RunAction::WriteFluenceCSV(double energy_keV, uint64_t n_run) const {
    std::string flu_file = output_csv_;
    auto dot = flu_file.rfind('.');
    if (dot != std::string::npos)
        flu_file = flu_file.substr(0, dot) + "_fluence.csv";
    else
        flu_file += "_fluence.csv";

    std::ofstream csv(flu_file);
    if (!csv) {
        std::cerr << "RunAction: cannot open fluence file: " << flu_file << "\n";
        return;
    }

    const double N = (n_run > 0) ? static_cast<double>(n_run) : 1.0;
    // Per-source-gamma normalization: each event emits exactly one source gamma
    // (cone-biased events carry weight = cone solid-angle fraction). The unbiased
    // 4pi-per-source-gamma estimate of any tally is sum(w * I) / N.
    const double unc_per = unc_w_ / N;
    const double unc_per_sig = std::sqrt(unc_w2_) / N;
    const double scat_per = scat_w_ / N;
    const double scat_per_sig = std::sqrt(scat_w2_) / N;

    csv << "# Photon fluence incident on detector (air-scatter study)\n"
        << "# Primary energy (E0): " << energy_keV << " keV\n"
        << "# Events: " << n_run << "\n"
        << "# Cone omega_frac: " << cone_omega_frac_
        << " (0 = isotropic)\n"
        << "# Bin width: " << (bin_width_MeV_ * 1000.0) << " keV\n"
        << "# Normalization: fluence_per_src_gamma = weighted_counts / Events\n"
        << "# Uncollided (|E-E0|<1.5keV) per source gamma: "
        << std::scientific << std::setprecision(6) << unc_per
        << " +/- " << unc_per_sig << "\n"
        << "# Scattered (E<E0-1.5keV) per source gamma:  "
        << scat_per << " +/- " << scat_per_sig << "\n"
        << "# Scatter/uncollided ratio: "
        << (unc_per > 0 ? scat_per / unc_per : 0.0) << "\n"
        << "Energy (keV),fluence_per_src_gamma,uncertainty\n";

    csv << std::scientific << std::setprecision(6);
    for (int i = 0; i < num_bins_; ++i) {
        if (fluence_w_[i] > 0.0) {
            double lower_keV = i * bin_width_MeV_ * 1000.0;
            double per = fluence_w_[i] / N;
            double sig = std::sqrt(fluence_w2_[i]) / N;
            csv << std::fixed << std::setprecision(3) << lower_keV << ","
                << std::scientific << std::setprecision(6) << per << ","
                << sig << "\n";
        }
    }

    std::cout << "  Fluence:         " << flu_file << "\n"
              << "    uncollided/src=" << std::scientific << std::setprecision(4)
              << unc_per << " +/- " << unc_per_sig
              << "   scattered/src=" << scat_per << " +/- " << scat_per_sig << "\n";
}

void RunAction::WriteHistogramCSV(double energy_keV, uint64_t n_run) const {
    std::string histo_file = output_csv_;
    auto dot = histo_file.rfind('.');
    if (dot != std::string::npos)
        histo_file = histo_file.substr(0, dot) + "_histogram.csv";
    else
        histo_file += "_histogram.csv";

    std::ofstream csv(histo_file);
    if (!csv) {
        std::cerr << "RunAction: cannot open histogram file: " << histo_file << "\n";
        return;
    }

    // RealTime for equivalent 1E6 Bq source, 100% BR
    // Cone-biased: N_equiv_iso = n_run / omega_frac; Isotropic: N_equiv_iso = n_run
    constexpr double kActivity = 1.0e6; // Bq
    double n_equiv = (cone_omega_frac_ > 0.0)
        ? static_cast<double>(n_run) / cone_omega_frac_
        : static_cast<double>(n_run);
    double real_time_s = n_equiv / kActivity;

    csv << "# Energy deposit histogram\n"
        << "# Primary energy: " << energy_keV << " keV\n"
        << "# Events: " << n_run << "\n"
        << "# Bin width: " << (bin_width_MeV_ * 1000.0) << " keV\n"
        << "RealTime: " << std::fixed << std::setprecision(3) << real_time_s << " s\n"
        << "energy_lower_keV,component,sum_w,sum_w2\n";
    // Undo the RealTime line's fixed/3-decimal format: cone-weighted sums are
    //  small (sum_w2 ~ n w^2) and printed as 0.000 under it.
    csv << std::defaultfloat << std::setprecision(12);

    static constexpr const char* component_names[] = {
        "zero_deposit", "full_energy", "pair_full_energy",
        "single_escape", "double_escape", "annihilation_511",
        "photoelectric_escape", "single_compton", "multiple_compton",
        "pair_no_compton", "pair_single_compton",
        "pair_multiple_compton"};
    for (int i = 0; i < num_bins_; ++i) {
        if (histogram_sumw_[i] > 0.0) {
            double lower_keV = i * bin_width_MeV_ * 1000.0;
            csv << lower_keV << ",total," << histogram_sumw_[i]
                << ',' << histogram_sumw2_[i] << "\n";
            if (component_labels_enabled_) {
                for (int component = 0;
                     component < kNumSpectrumComponents; ++component) {
                    if (component_histogram_sumw_[component][i] > 0.0) {
                        csv << lower_keV << ','
                            << component_names[component] << ','
                            << component_histogram_sumw_[component][i] << ','
                            << component_histogram_sumw2_[component][i] << '\n';
                    }
                }
            }
        }
    }

    std::cout << "  Histogram:       " << histo_file << "\n";
}

void RunAction::WriteSummaryCSV(double energy_keV, uint64_t n_run) const {
    std::string summary_file = output_csv_;
    auto dot = summary_file.rfind('.');
    if (dot != std::string::npos)
        summary_file = summary_file.substr(0, dot) + "_summary.csv";
    else
        summary_file += "_summary.csv";

    std::ofstream csv(summary_file);
    if (!csv) {
        std::cerr << "RunAction: cannot open summary file: " << summary_file << "\n";
        return;
    }

    double w_total = (sum_weights_ > 0) ? sum_weights_ : static_cast<double>(n_run);
    double avg_total = (sum_deposit_weights_ > 0.0)
        ? (sum_total_edep_MeV_ * 1000.0 / sum_deposit_weights_) : 0.0;
    double avg_electron = (sum_deposit_weights_ > 0.0)
        ? (sum_electron_edep_MeV_ * 1000.0 / sum_deposit_weights_) : 0.0;
    double avg_gamma = (sum_deposit_weights_ > 0.0)
        ? (sum_gamma_edep_MeV_ * 1000.0 / sum_deposit_weights_) : 0.0;

    double frac_electron = (sum_total_edep_MeV_ > 0.0)
        ? (sum_electron_edep_MeV_ / sum_total_edep_MeV_) : 0.0;
    double frac_gamma = (sum_total_edep_MeV_ > 0.0)
        ? (sum_gamma_edep_MeV_ / sum_total_edep_MeV_) : 0.0;

    csv << "# Per-particle energy deposition summary\n"
        << "# Primary energy: " << energy_keV << " keV\n"
        << "# Events: " << n_run << "\n"
        << "# Events with deposit: " << n_events_with_deposit_ << "\n"
        << "particle_type,total_edep_keV,avg_edep_per_event_keV,fraction_of_total\n"
        << "electron," << (sum_electron_edep_MeV_ * 1000.0) << ","
        << avg_electron << "," << frac_electron << "\n"
        << "gamma," << (sum_gamma_edep_MeV_ * 1000.0) << ","
        << avg_gamma << "," << frac_gamma << "\n"
        << "total," << (sum_total_edep_MeV_ * 1000.0) << ","
        << avg_total << ",1.0\n";

    std::cout << "  Summary:         " << summary_file << "\n";
    std::cout << "  Edep breakdown:  "
              << "e-/e+ " << (frac_electron * 100.0) << "%, "
              << "gamma " << (frac_gamma * 100.0) << "%, "
              << "other " << ((1.0 - frac_electron - frac_gamma) * 100.0) << "%\n";
}

void RunAction::WriteDiagnosticCSVs(double energy_keV, uint64_t n_run) const {
    // Derive prefix from output filename
    std::string base = output_csv_;
    auto dot = base.rfind('.');
    if (dot != std::string::npos) base = base.substr(0, dot);

    std::string pfx = base + "_diag_";

    auto write_csv = [&](const std::string& suffix,
                          const std::string& x_label,
                          const std::vector<uint64_t>& counts,
                          double bin_min, double bin_width, int n_bins,
                          const std::string& description)
    {
        std::string filename = pfx + suffix + ".csv";
        std::ofstream csv(filename);
        if (!csv) {
            std::cerr << "RunAction: cannot open diagnostic file: " << filename << "\n";
            return;
        }
        csv << "# G4 diag, " << energy_keV << " keV, "
            << n_run << " events - " << description << "\n"
            << x_label << ",counts\n";
        for (int i = 0; i < n_bins; ++i) {
            if (counts[i] > 0) {
                double center = bin_min + (i + 0.5) * bin_width;
                csv << std::fixed << std::setprecision(2) << center
                    << "," << counts[i] << "\n";
            }
        }
        std::cout << "  Diagnostic: " << filename << "\n";
    };

    write_csv("cs1_depth",  "depth_cm", diag_cs1_depth_,
              0.0, kDiagDepthWidth, kDiagDepthBins, "1st Compton depth");
    write_csv("cs1_angle",  "angle_deg", diag_cs1_angle_,
              0.0, kDiagAngleWidth, kDiagAngleBins, "1st Compton angle");
    write_csv("cs1_energy", "energy_keV", diag_cs1_energy_,
              0.0, kDiagEnergyWidth, kDiagEnergyBins, "1st Compton energy");

    write_csv("cs2_depth",  "depth_cm", diag_cs2_depth_,
              0.0, kDiagDepthWidth, kDiagDepthBins, "2nd Compton depth");
    write_csv("cs2_angle",  "angle_deg", diag_cs2_angle_,
              0.0, kDiagAngleWidth, kDiagAngleBins, "2nd Compton angle");
    write_csv("cs2_energy", "energy_keV", diag_cs2_energy_,
              0.0, kDiagEnergyWidth, kDiagEnergyBins, "2nd Compton energy");

    write_csv("pe_depth",   "depth_cm", diag_pe_depth_,
              0.0, kDiagDepthWidth, kDiagDepthBins, "PE depth");
    write_csv("pe_energy",  "energy_keV", diag_pe_energy_,
              0.0, kDiagEnergyWidth, kDiagEnergyBins, "PE incoming energy");
    write_csv("pe_nscatter","n_scatters", diag_pe_nscatter_,
              0.0, 1.0, kDiagScatterBins, "# Compton before PE");

    write_csv("scatter_mult", "n_comptons", diag_scatter_mult_,
              0.0, 1.0, kDiagScatterBins, "Compton scatter multiplicity");
    write_csv("escape_energy","energy_keV", diag_escape_energy_,
              0.0, kDiagEnergyWidth, kDiagEnergyBins, "Escape photon energy");

    write_csv("fluor_escape","energy_keV", diag_fluor_escape_,
              0.0, kDiagEnergyWidth, kDiagEnergyBins, "Secondary gamma escape energy");
    std::cout << "  Total secondary gammas escaped: " << n_fluor_escaped_total_ << "\n";
}

void RunAction::RecordCrystalEntries(
    const std::vector<G4CrystalEntryRecord>& entries)
{
    std::lock_guard<std::mutex> lk(mutex_);

    for (const auto& rec : entries) {
        entry_total_++;

        int e_bin = static_cast<int>(rec.energy_keV / kEntryEnergyWidth);
        if (e_bin < 0) e_bin = 0;
        if (e_bin >= kEntryEnergyBins) e_bin = kEntryEnergyBins - 1;

        double z_cm = rec.position.z();  // already in cm from SteppingAction
        int z_bin = static_cast<int>((z_cm - kEntryZMin) / kEntryZBinWidth);

        switch (rec.surface) {
            case G4CrystalEntryRecord::kFace:
                entry_face_count_++;
                entry_face_energy_sum_ += rec.energy_keV;
                if (e_bin >= 0 && e_bin < kEntryEnergyBins)
                    entry_energy_face_[e_bin]++;
                break;
            case G4CrystalEntryRecord::kSide:
                entry_side_count_++;
                entry_side_energy_sum_ += rec.energy_keV;
                if (e_bin >= 0 && e_bin < kEntryEnergyBins)
                    entry_energy_side_[e_bin]++;
                break;
            case G4CrystalEntryRecord::kBack:
                entry_back_count_++;
                entry_back_energy_sum_ += rec.energy_keV;
                break;
        }

        if (z_bin >= 0 && z_bin < kEntryZBins)
            entry_z_[z_bin]++;
    }
}

void RunAction::WriteEntryDiagCSVs(double energy_keV, uint64_t n_run) const {
    std::string base = output_csv_;
    auto dot = base.rfind('.');
    if (dot != std::string::npos) base = base.substr(0, dot);

    int e_int = static_cast<int>(energy_keV);

    // Energy histogram
    {
        std::string fname = base + "_entry_energy.csv";
        std::ofstream f(fname);
        f << "# G4 crystal entry energy histogram\n"
          << "# Source energy: " << energy_keV << " keV\n"
          << "# Events: " << n_run << "\n"
          << "energy_keV,face_counts,side_counts\n";
        for (int i = 0; i < kEntryEnergyBins; ++i) {
            if (entry_energy_face_[i] > 0 || entry_energy_side_[i] > 0) {
                f << static_cast<int>(i * kEntryEnergyWidth) << ","
                  << entry_energy_face_[i] << "," << entry_energy_side_[i] << "\n";
            }
        }
        std::cout << "  Entry energy: " << fname << "\n";
    }

    // Z histogram
    {
        std::string fname = base + "_entry_z.csv";
        std::ofstream f(fname);
        f << "# G4 crystal entry z-position histogram\n"
          << "# Source energy: " << energy_keV << " keV\n"
          << "# Events: " << n_run << "\n"
          << "z_cm,counts\n";
        for (int i = 0; i < kEntryZBins; ++i) {
            if (entry_z_[i] > 0) {
                f << std::fixed << std::setprecision(3)
                  << (kEntryZMin + (i + 0.5) * kEntryZBinWidth) << ","
                  << entry_z_[i] << "\n";
            }
        }
        std::cout << "  Entry z-pos:  " << fname << "\n";
    }

    // Summary
    {
        std::string fname = base + "_entry_summary.csv";
        std::ofstream f(fname);
        f << "# G4 crystal entry summary\n"
          << "# Source energy: " << energy_keV << " keV\n"
          << "# Events: " << n_run << "\n"
          << "quantity,value\n"
          << "total_events," << n_run << "\n"
          << "total_entries," << entry_total_ << "\n"
          << "face_entries," << entry_face_count_ << "\n"
          << "side_entries," << entry_side_count_ << "\n"
          << "back_entries," << entry_back_count_ << "\n"
          << std::fixed << std::setprecision(2)
          << "face_mean_energy_keV,"
          << (entry_face_count_ > 0 ? entry_face_energy_sum_ / entry_face_count_ : 0.0) << "\n"
          << "side_mean_energy_keV,"
          << (entry_side_count_ > 0 ? entry_side_energy_sum_ / entry_side_count_ : 0.0) << "\n"
          << "back_mean_energy_keV,"
          << (entry_back_count_ > 0 ? entry_back_energy_sum_ / entry_back_count_ : 0.0) << "\n"
          << std::setprecision(4)
          << "face_fraction,"
          << (entry_total_ > 0 ? static_cast<double>(entry_face_count_) / entry_total_ : 0.0) << "\n"
          << "side_fraction,"
          << (entry_total_ > 0 ? static_cast<double>(entry_side_count_) / entry_total_ : 0.0) << "\n"
          << "back_fraction,"
          << (entry_total_ > 0 ? static_cast<double>(entry_back_count_) / entry_total_ : 0.0) << "\n";
        std::cout << "  Entry summary: " << fname << "\n";
    }

    // Print summary to console
    std::cout << "\n=== Crystal Entry Diagnostics ===\n"
              << "  Total entries: " << entry_total_ << "\n"
              << "  Face (z~0):  " << entry_face_count_
              << " (" << std::fixed << std::setprecision(1)
              << 100.0 * entry_face_count_ / std::max(entry_total_, uint64_t(1)) << "%)"
              << "  mean E=" << std::setprecision(1)
              << (entry_face_count_ > 0 ? entry_face_energy_sum_ / entry_face_count_ : 0.0) << " keV\n"
              << "  Side (r~R):  " << entry_side_count_
              << " (" << std::setprecision(1)
              << 100.0 * entry_side_count_ / std::max(entry_total_, uint64_t(1)) << "%)"
              << "  mean E=" << std::setprecision(1)
              << (entry_side_count_ > 0 ? entry_side_energy_sum_ / entry_side_count_ : 0.0) << " keV\n"
              << "  Back (z~L):  " << entry_back_count_
              << " (" << std::setprecision(1)
              << 100.0 * entry_back_count_ / std::max(entry_total_, uint64_t(1)) << "%)"
              << "  mean E=" << std::setprecision(1)
              << (entry_back_count_ > 0 ? entry_back_energy_sum_ / entry_back_count_ : 0.0) << " keV\n";
}

// --- Environment-scatter sphere flux scorer ---------------------------------

void RunAction::RecordSphereStep(double e_keV, double mu, double len_cm,
                                 bool scattered, int origin)
{
    if (!sphere_flux_enabled_ || len_cm <= 0.0) return;
    if (e_keV < kSphereEMin || e_keV > kSphereEMax) return;
    const double u = std::log(e_keV / kSphereEMin) /
                     std::log(kSphereEMax / kSphereEMin);
    int ie = static_cast<int>(u * kSphereNE);
    if (ie < 0) ie = 0;
    if (ie >= kSphereNE) ie = kSphereNE - 1;
    if (mu < -1.0) mu = -1.0;
    if (mu > 1.0) mu = 1.0;
    int im = static_cast<int>((mu + 1.0) * 0.5 * kSphereNMu);
    if (im < 0) im = 0;
    if (im >= kSphereNMu) im = kSphereNMu - 1;
    const int idx = ie * kSphereNMu + im;
    if (scattered) sphere_scat_ev_[idx] += len_cm;
    else           sphere_direct_[idx]  += len_cm;
    if (scattered && origin >= 0 && origin < kNumSphereOrigins)
        sphere_origin_[origin * kSphereNE + ie] += len_cm;
}

void RunAction::EndSphereEvent()
{
    if (!sphere_flux_enabled_) return;
    // Per-event summation before squaring: steps within one event are
    // correlated, so squaring individual steps would understate the variance.
    // Kahan summation limits roundoff while reducing the occupied (E,mu)
    // cells to the scalar event contribution.
    double event_total_cm = 0.0;
    double event_total_compensation_cm = 0.0;
    for (std::size_t i = 0; i < sphere_scat_ev_.size(); ++i) {
        const double v = sphere_scat_ev_[i];
        if (v == 0.0) continue;
        const double corrected = v - event_total_compensation_cm;
        const double updated = event_total_cm + corrected;
        event_total_compensation_cm = (updated - event_total_cm) - corrected;
        event_total_cm = updated;
        sphere_scat_[i]    += v;
        sphere_scat_w2_[i] += v * v;
        sphere_scat_ev_[i]  = 0.0;
    }
    ++sphere_scat_total_n_;
    const long double total = static_cast<long double>(event_total_cm);
    sphere_scat_total_s1_cm_ += total;
    sphere_scat_total_s2_cm2_ += total * total;
}

void RunAction::WriteSphereFluxCSV(uint64_t n_run) const
{
    if (!sphere_flux_enabled_ || n_run == 0) return;
    std::string path = output_csv_;
    const std::size_t dot = path.rfind('.');
    if (dot != std::string::npos) path = path.substr(0, dot);
    path += "_sphereflux.csv";

    std::ofstream os(path);
    if (!os) {
        G4cout << "ERROR: cannot write " << path << G4endl;
        return;
    }
    const double R = sphere_radius_cm_;
    const double volume = (4.0 / 3.0) * CLHEP::pi * R * R * R;

    const long double moment_n =
        static_cast<long double>(sphere_scat_total_n_);
    const long double volume_ld = static_cast<long double>(volume);
    const long double total_mean_flux =
        sphere_scat_total_n_ > 0
            ? sphere_scat_total_s1_cm_ / (moment_n * volume_ld)
            : 0.0L;
    long double centered_total = 0.0L;
    if (sphere_scat_total_n_ > 0) {
        centered_total = sphere_scat_total_s2_cm2_
            - sphere_scat_total_s1_cm_ * sphere_scat_total_s1_cm_ / moment_n;
        // Roundoff can produce a tiny negative value for a deterministic
        // tally.  Do not allow it to create a NaN in the explicit SE header.
        if (centered_total < 0.0L) centered_total = 0.0L;
    }
    const long double total_se_flux =
        sphere_scat_total_n_ > 1
            ? std::sqrt(centered_total /
                        (moment_n * (moment_n - 1.0L))) / volume_ld
            : 0.0L;

    os << std::scientific
       << std::setprecision(std::numeric_limits<long double>::max_digits10);
    os << "# GEANT4 environment-scatter sphere flux\n"
       << "# Track-length flux estimator phi = L / V inside a scoring sphere.\n"
       << "# ALL photon generations scored. Direct = primary before its first\n"
       << "# interaction; scattered = everything else. They are mutually\n"
       << "# exclusive outcomes of one primary and must not be summed.\n"
       << "# E0_keV: " << sphere_e0_keV_ << "\n"
       << "# sphere_radius_cm: " << R << "\n"
       << "# gdml_scoring_radius_cm: " << sphere_gdml_scoring_radius_cm_ << "\n"
       << "# sphere_volume_cm3: " << volume << "\n"
       << "# observable: scattered-track-length-flux-per-primary-v1\n"
       << "# physics_list: " << sphere_physics_list_ << "\n"
       << "# world_half_extent_m: " << sphere_world_half_extent_m_ << "\n"
       << "# primary_sampling: " << sphere_primary_sampling_ << "\n"
       << "# rng_engine: " << sphere_rng_engine_ << "\n"
       << "# rng_seed: " << sphere_rng_seed_ << "\n"
       << "# production_cut_mm: " << sphere_production_cut_mm_ << "\n"
       << "# gdml_sha256: " << sphere_gdml_sha256_ << "\n"
       << "# macro_sha256: " << sphere_macro_sha256_ << "\n"
       << "# geometry_tier: " << sphere_geometry_tier_ << "\n"
       << "# material: " << sphere_material_ << "\n"
       << "# distance_cm: " << sphere_distance_cm_ << "\n"
       << "# source_height_cm: " << sphere_source_height_cm_ << "\n"
       << "# detector_height_cm: " << sphere_detector_height_cm_ << "\n"
       << "# events: " << n_run << "\n"
       << "# scattered_total_moment_schema: per-event-track-length-v1\n"
       << "# scattered_total_N: " << sphere_scat_total_n_ << "\n"
       << "# scattered_total_S1_cm: " << sphere_scat_total_s1_cm_ << "\n"
       << "# scattered_total_S2_cm2: " << sphere_scat_total_s2_cm2_ << "\n"
       << "# scattered_total_mean_flux_cm-2: " << total_mean_flux << "\n"
       << "# scattered_total_standard_error_flux_cm-2: "
       << total_se_flux << "\n"
       << "# n_energy: " << kSphereNE << "\n"
       << "# energy_min_keV: " << kSphereEMin << "\n"
       << "# energy_max_keV: " << kSphereEMax << "\n"
       << "# n_mu: " << kSphereNMu << "\n"
       << "# mu_min: -1\n"
       << "# mu_max: 1\n"
       << "# mu_definition: cos(direction-of-travel,+z)\n"
       << "# origin classes below decompose the SCATTERED flux by how each\n"
       << "# photon was created: primary (and its whole scatter chain),\n"
       << "# bremsstrahlung, annihilation, other (fluorescence). One run, so\n"
       << "# no differencing of independent runs is involved.\n"
       << "energy_bin,mu_bin,e_low_keV,mu_low,phi_direct,phi_scattered,"
          "phi_scattered_sigma\n";

    const double n = static_cast<double>(n_run);
    for (int ie = 0; ie < kSphereNE; ++ie) {
        const double elo = kSphereEMin *
            std::pow(kSphereEMax / kSphereEMin, double(ie) / kSphereNE);
        for (int im = 0; im < kSphereNMu; ++im) {
            const int idx = ie * kSphereNMu + im;
            const double mlo = -1.0 + 2.0 * double(im) / kSphereNMu;
            const double phid = sphere_direct_[idx] / (volume * n);
            const double m    = sphere_scat_[idx] / n;
            const double var  = std::max(0.0, sphere_scat_w2_[idx] / n - m * m);
            os << ie << "," << im << "," << elo << "," << mlo << ","
               << phid << "," << m / volume << ","
               << std::sqrt(var / n) / volume << "\n";
        }
    }

    // Per-origin decomposition of the scattered flux, as a trailing block with
    // its own header so existing parsers (which read 7-field rows) ignore it.
    if (!sphere_origin_.empty()) {
        static const char* kNames[kNumSphereOrigins] =
            {"primary", "brems_air", "brems_soil", "annihil", "other"};
        os << "# ORIGIN DECOMPOSITION (scattered flux only)\n";
        os << "origin,energy_bin,e_low_keV,phi_scattered\n";
        double tot[kNumSphereOrigins] = {};
        for (int o = 0; o < kNumSphereOrigins; ++o) {
            for (int ie = 0; ie < kSphereNE; ++ie) {
                const double elo = kSphereEMin *
                    std::pow(kSphereEMax / kSphereEMin, double(ie) / kSphereNE);
                const double v = sphere_origin_[o * kSphereNE + ie] / (volume * n);
                tot[o] += v;
                os << kNames[o] << "," << ie << "," << elo << "," << v << "\n";
            }
        }
        double all = 0.0;
        for (int o = 0; o < kNumSphereOrigins; ++o) all += tot[o];
        os << "# origin totals:";
        for (int o = 0; o < kNumSphereOrigins; ++o)
            os << " " << kNames[o] << "=" << tot[o]
               << " (" << (all > 0.0 ? 100.0 * tot[o] / all : 0.0) << "%)";
        os << "\n";
        G4cout << "  scattered-flux origins:";
        for (int o = 0; o < kNumSphereOrigins; ++o)
            G4cout << " " << kNames[o] << " "
                   << (all > 0.0 ? 100.0 * tot[o] / all : 0.0) << "%";
        G4cout << G4endl;
    }
    G4cout << "Wrote sphere flux to " << path << G4endl;
}
