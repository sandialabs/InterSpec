#pragma once
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

#include "G4UserRunAction.hh"
#include "SteppingAction.hh"  // for G4InteractionRecord

#include <string>
#include <cstdint>
#include <cmath>
#include <mutex>
#include <vector>

/// Accumulates per-event tallies across all threads and writes results to a
/// CSV file at the end of the run.
///
/// Supports importance-weighted tallies for cone-bias mode:
///   efficiency = sum_of_weights_for_hits / num_events
///
/// When histogram mode is enabled, also accumulates:
///   - Fixed-bin energy histogram (1 keV bins, weighted)
///   - Per-particle-type energy breakdown (electron vs gamma)
///
/// When diagnostics mode is enabled, also accumulates:
///   - 11 per-interaction diagnostic histograms
class RunAction : public G4UserRunAction {
public:
    RunAction(const std::string& output_csv,
              const std::string& gdml_file,
              const std::string& macro_file);
    ~RunAction() override = default;

    void BeginOfRunAction(const G4Run* run) override;
    void EndOfRunAction(const G4Run* run) override;

    /// Called by EventAction (thread-safe).
    void RecordFEP(double /*edep_MeV*/, double weight) {
        std::lock_guard<std::mutex> lk(mutex_);
        w_fep_ += weight;
        w2_fep_ += weight * weight;
    }
    void RecordAny(double /*edep_MeV*/, double weight) {
        std::lock_guard<std::mutex> lk(mutex_);
        w_any_ += weight;
        w2_any_ += weight * weight;
    }

    /// Called every event to accumulate total weight (for normalization).
    void RecordWeight(double weight) {
        std::lock_guard<std::mutex> lk(mutex_);
        sum_weights_ += weight;
    }

    /// Called by EventAction when histogram mode is enabled.
    void RecordHistogram(double total_edep_MeV,
                         double electron_edep_MeV,
                         double gamma_edep_MeV,
                         double weight);

    /// Called by EventAction when diagnostics mode is enabled.
    void RecordDiagnostics(const std::vector<G4InteractionRecord>& interactions,
                           double edep_MeV,
                           double primary_energy_MeV);

    /// Called by EventAction: record a secondary gamma escaping the crystal.
    void RecordFluorEscape(double escape_energy_keV);

    /// Enable diagnostic histogram output.
    void EnableHistogram(bool on) { histogram_enabled_ = on; }

    /// Enable per-interaction diagnostic histograms.
    void EnableDiagnostics(bool on) { diagnostics_enabled_ = on; }

    /// Enable crystal entry diagnostics.
    void EnableEntryDiag(bool on) { entry_diag_enabled_ = on; }

    /// Called by EventAction: record crystal entry records for this event.
    void RecordCrystalEntries(const std::vector<G4CrystalEntryRecord>& entries);

    /// Called by EventAction: record external electron entry into crystal.
    /// @param had_deposit  Whether this event had any crystal energy deposit
    /// @param gamma_entered Whether a gamma also entered the crystal this event
    void RecordExternalElectron(bool had_deposit, bool gamma_entered) {
        std::lock_guard<std::mutex> lk(mutex_);
        n_ext_electron_events_++;
        if (had_deposit) n_ext_electron_with_deposit_++;
        if (had_deposit && !gamma_entered) n_ext_electron_only_deposit_++;
    }

    /// Called by EventAction every event: u/s source-class decomposition
    /// (s = primary interacted outside the crystal before first reaching
    /// it). Mirrors CeeLo's SourceEffectDiag for direct comparison.
    void RecordSrcClass(bool interacted, bool any, bool fep, double w) {
        std::lock_guard<std::mutex> lk(mutex_);
        if (interacted) {
            n_src_s_++;
            if (any) { src_any_w_s_ += w; src_any_w2_s_ += w * w; }
            if (fep) { src_fep_w_s_ += w; src_fep_w2_s_ += w * w; }
        } else {
            n_src_u_++;
            if (any) { src_any_w_u_ += w; src_any_w2_u_ += w * w; }
            if (fep) { src_fep_w_u_ += w; src_fep_w2_u_ += w * w; }
        }
    }

    /// Set cone solid angle fraction for RealTime header in histograms.
    /// If 0 (default), assumes isotropic (4pi) emission.
    void SetConeSolidAngleFraction(double f) { cone_omega_frac_ = f; }

    /// Set histogram bin width (default: 1 keV = 0.001 MeV).
    void SetBinWidth(double keV) {
        bin_width_MeV_ = keV * 0.001;
        num_bins_ = static_cast<int>(std::ceil(kMaxEnergyMeV / bin_width_MeV_));
    }

    /// Check if target FEP relative precision has been achieved.
    /// @param target  Target relative precision (e.g. 0.01 for 1%)
    /// @return true if precision achieved or events exceeded safety cap
    bool CheckFepPrecision(double target) const {
        std::lock_guard<std::mutex> lk(mutex_);
        if (n_total_ < 10000 || w_fep_ <= 0.0) return false;
        double N = static_cast<double>(n_total_);
        double eps = w_fep_ / N;
        double var = (w2_fep_ / N - eps * eps) / N;
        double sigma = std::sqrt(std::max(0.0, var));
        return (sigma / eps) <= target;
    }

    /// Get total events processed so far.
    uint64_t GetTotalEvents() const {
        std::lock_guard<std::mutex> lk(mutex_);
        return n_total_;
    }

    /// Accumulate tallies from a sub-run (for iterative beamOn).
    /// BeginOfRunAction resets tallies, so we need to preserve across calls.
    void PreserveAcrossRuns() { preserve_across_runs_ = true; }

private:
    std::string  output_csv_;
    std::string  gdml_file_;
    std::string  macro_file_;
    mutable std::mutex mutex_;
    uint64_t     n_total_ = 0;
    bool         preserve_across_runs_ = false;

    // Weighted tallies: efficiency = w_fep / N, uncertainty from w2
    double w_fep_ = 0.0;
    double w2_fep_ = 0.0;
    double w_any_ = 0.0;
    double w2_any_ = 0.0;
    double sum_weights_ = 0.0;

    // Histogram mode
    bool histogram_enabled_ = false;
    static constexpr double kMaxEnergyMeV = 3.0; // covers up to 3 MeV
    double bin_width_MeV_ = 0.001; // default 1 keV bins
    int num_bins_ = 3000;
    std::vector<uint64_t> histogram_;  // raw (unweighted) hit counts for spectral shape

    // Per-particle-type cumulative energy (MeV), weighted
    double sum_electron_edep_MeV_ = 0.0;
    double sum_gamma_edep_MeV_ = 0.0;
    double sum_total_edep_MeV_ = 0.0;
    uint64_t n_events_with_deposit_ = 0;

    double cone_omega_frac_ = 0.0; // 0 = isotropic

    void WriteHistogramCSV(double energy_keV, uint64_t n_run) const;
    void WriteSummaryCSV(double energy_keV, uint64_t n_run) const;

    // Diagnostics mode: 11 histograms
    bool diagnostics_enabled_ = false;

    // Histogram dimensions (must match MC side)
    static constexpr int kDiagDepthBins  = 80;
    static constexpr double kDiagDepthWidth = 0.1;  // cm
    static constexpr int kDiagAngleBins  = 180;
    static constexpr double kDiagAngleWidth = 1.0;  // deg
    static constexpr int kDiagEnergyBins = 3000;  // covers up to 3 MeV
    static constexpr double kDiagEnergyWidth = 1.0;   // keV
    static constexpr int kDiagScatterBins = 21;

    std::vector<uint64_t> diag_cs1_depth_;
    std::vector<uint64_t> diag_cs1_angle_;
    std::vector<uint64_t> diag_cs1_energy_;
    std::vector<uint64_t> diag_cs2_depth_;
    std::vector<uint64_t> diag_cs2_angle_;
    std::vector<uint64_t> diag_cs2_energy_;
    std::vector<uint64_t> diag_pe_depth_;
    std::vector<uint64_t> diag_pe_energy_;
    std::vector<uint64_t> diag_pe_nscatter_;
    std::vector<uint64_t> diag_scatter_mult_;
    std::vector<uint64_t> diag_escape_energy_;

    // Secondary photon (fluorescence) escape histogram
    std::vector<uint64_t> diag_fluor_escape_;
    uint64_t n_fluor_escaped_total_ = 0;

    void WriteDiagnosticCSVs(double energy_keV, uint64_t n_run) const;

    // Crystal entry diagnostics
    bool entry_diag_enabled_ = false;
    static constexpr int kEntryEnergyBins = 3000;  // 1 keV bins, covers up to 3 MeV
    static constexpr double kEntryEnergyWidth = 1.0; // keV
    static constexpr int kEntryZBins = 200;          // 0.05 cm bins, -1 to 9 cm
    static constexpr double kEntryZBinWidth = 0.05;  // cm
    static constexpr double kEntryZMin = -1.0;       // cm

    std::vector<uint64_t> entry_energy_face_;
    std::vector<uint64_t> entry_energy_side_;
    std::vector<uint64_t> entry_z_;

    uint64_t entry_face_count_ = 0;
    uint64_t entry_side_count_ = 0;
    uint64_t entry_back_count_ = 0;
    double entry_face_energy_sum_ = 0.0;
    double entry_side_energy_sum_ = 0.0;
    double entry_back_energy_sum_ = 0.0;
    uint64_t entry_total_ = 0;

    void WriteEntryDiagCSVs(double energy_keV, uint64_t n_run) const;

    // External electron entry counters
    uint64_t n_ext_electron_events_ = 0;        ///< Events with any e-/e+ boundary crossing into crystal
    uint64_t n_ext_electron_with_deposit_ = 0;   ///< Events with ext electron AND crystal deposit
    uint64_t n_ext_electron_only_deposit_ = 0;   ///< Events with ext electron AND deposit but NO gamma entered

    // u/s source-class decomposition (s = primary interacted outside the
    // crystal before first reaching it)
    uint64_t n_src_u_ = 0, n_src_s_ = 0;
    double src_any_w_u_ = 0.0, src_any_w2_u_ = 0.0;
    double src_any_w_s_ = 0.0, src_any_w2_s_ = 0.0;
    double src_fep_w_u_ = 0.0, src_fep_w2_u_ = 0.0;
    double src_fep_w_s_ = 0.0, src_fep_w2_s_ = 0.0;
};
