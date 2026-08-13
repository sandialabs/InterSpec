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

    // Reset histogram data
    if (histogram_enabled_) {
        histogram_.assign(num_bins_, 0);
        sum_electron_edep_MeV_ = 0.0;
        sum_gamma_edep_MeV_ = 0.0;
        sum_total_edep_MeV_ = 0.0;
        n_events_with_deposit_ = 0;
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
                                double weight) {
    std::lock_guard<std::mutex> lk(mutex_);

    // Fill histogram bin (raw counts -- unweighted for spectral shape)
    if (total_edep_MeV > 0.0) {
        int bin = static_cast<int>(total_edep_MeV / bin_width_MeV_);
        if (bin >= 0 && bin < num_bins_) {
            histogram_[bin]++;
        }

        sum_total_edep_MeV_ += total_edep_MeV;
        sum_electron_edep_MeV_ += electron_edep_MeV;
        sum_gamma_edep_MeV_ += gamma_edep_MeV;
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

    // Read primary energy from the GPS macro (first /gps/ene/mono line).
    double energy_keV = 0.0;
    {
        std::ifstream mac(macro_file_);
        std::string line;
        while (std::getline(mac, line)) {
            if (line.find("/gps/ene/mono") != std::string::npos) {
                std::istringstream iss(line);
                std::string cmd;
                double val;
                std::string unit;
                if (iss >> cmd >> val >> unit) {
                    energy_keV = (unit == "MeV") ? val * 1000.0 : val;
                }
                break;
            }
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
        << "Energy (keV),Counts\n";

    for (int i = 0; i < num_bins_; ++i) {
        if (histogram_[i] > 0) {
            double lower_keV = i * bin_width_MeV_ * 1000.0;
            csv << lower_keV << "," << histogram_[i] << "\n";
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
    double avg_total = (n_events_with_deposit_ > 0)
        ? (sum_total_edep_MeV_ * 1000.0 / n_events_with_deposit_) : 0.0;
    double avg_electron = (n_events_with_deposit_ > 0)
        ? (sum_electron_edep_MeV_ * 1000.0 / n_events_with_deposit_) : 0.0;
    double avg_gamma = (n_events_with_deposit_ > 0)
        ? (sum_gamma_edep_MeV_ * 1000.0 / n_events_with_deposit_) : 0.0;

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
