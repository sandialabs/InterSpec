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
#include <array>
#include <algorithm>
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
                         double weight,
                         const std::vector<G4InteractionRecord>& interactions,
                         double primary_energy_MeV,
                         bool had_pair_production,
                         int pair_secondary_compton_count);

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

    /// Enable the photon-fluence-at-detector scorer (air-scatter study).
    void EnableFluence(bool on) { fluence_enabled_ = on; }

    // --- Environment-scatter sphere flux scorer -------------------------------
    // Scores the photon TRACK LENGTH inside the scoring sphere, binned in
    // (E, mu), where mu is the cosine between the direction of travel and +z.
    // phi = L / V is the standard track-length flux estimator, and it is exactly
    // the quantity CeeLo's NextEventEstimator produces -- scoring the same
    // observable both sides is what makes the comparison apples-to-apples.
    // ALL generations are scored (not just the primary): environment scatter is
    // largely secondaries and later-generation descendants.
    void EnableSphereFlux(bool on, double radius_cm) {
        sphere_flux_enabled_ = on;
        sphere_radius_cm_ = radius_cm;
        if (on) {
            sphere_direct_.assign(kSphereNE * kSphereNMu, 0.0);
            sphere_scat_.assign(kSphereNE * kSphereNMu, 0.0);
            sphere_scat_ev_.assign(kSphereNE * kSphereNMu, 0.0);
            sphere_scat_w2_.assign(kSphereNE * kSphereNMu, 0.0);
            sphere_origin_.assign(kNumSphereOrigins * kSphereNE, 0.0);
            sphere_scat_total_n_ = 0;
            sphere_scat_total_s1_cm_ = 0.0L;
            sphere_scat_total_s2_cm2_ = 0.0L;
        }
    }
    bool SphereFluxEnabled() const { return sphere_flux_enabled_; }

    /// Record the run configuration determined by the process which actually
    /// ran the scorer, including the SHA-256 of its exact GDML and macro
    /// inputs.  These fields are emitted into the source-hashed CSV; the fit
    /// evidence generator does not accept manifest-only assertions.
    void SetSphereRunMetadata(const std::string& physics_list,
                              double world_half_extent_m,
                              const std::string& primary_sampling,
                              const std::string& rng_engine,
                              long rng_seed,
                              double production_cut_mm,
                              const std::string& gdml_sha256,
                              const std::string& macro_sha256,
                              const std::string& geometry_tier,
                              const std::string& material,
                              double distance_cm,
                              double source_height_cm,
                              double detector_height_cm,
                              double gdml_scoring_radius_cm) {
        sphere_physics_list_ = physics_list;
        sphere_world_half_extent_m_ = world_half_extent_m;
        sphere_primary_sampling_ = primary_sampling;
        sphere_rng_engine_ = rng_engine;
        sphere_rng_seed_ = rng_seed;
        sphere_production_cut_mm_ = production_cut_mm;
        sphere_gdml_sha256_ = gdml_sha256;
        sphere_macro_sha256_ = macro_sha256;
        sphere_geometry_tier_ = geometry_tier;
        sphere_material_ = material;
        sphere_distance_cm_ = distance_cm;
        sphere_source_height_cm_ = source_height_cm;
        sphere_detector_height_cm_ = detector_height_cm;
        sphere_gdml_scoring_radius_cm_ = gdml_scoring_radius_cm;
    }

    /// One step of `len_cm` inside the sphere at energy `e_keV`, travelling at
    /// cosine `mu` to +z. `scattered` is false only for the primary photon
    /// before its first interaction (matching CeeLo's direct/scattered split,
    /// which are mutually exclusive outcomes of one primary and must never be
    /// summed).
    /// `origin` classifies the photon by how it was CREATED, so one run
    /// decomposes the scattered flux instead of needing a second run with a
    /// process inactivated. Differencing two independent runs cannot resolve a
    /// few-percent channel: at 36M events each, the difference carries the sum
    /// of both variances and every individual bin is noise.
    ///
    /// A Compton- or Rayleigh-scattered photon keeps its track ID and creator in
    /// GEANT4, so each class automatically includes that photon's whole
    /// subsequent scatter chain -- which is exactly what CeeLo's channel A/B
    /// switches measure, making the two directly comparable.
    /// Brems is split by the material of the EMISSION vertex (air vs slab) --
    /// air-born electrons have metres-scale range, so where their photons are
    /// born decides whether creation-point emission is a valid approximation.
    enum SphereOrigin { kOriginPrimary = 0, kOriginBremsAir, kOriginBremsSoil,
                        kOriginAnnihil, kOriginOther, kNumSphereOrigins };
    void RecordSphereStep(double e_keV, double mu, double len_cm, bool scattered,
                          int origin = kOriginPrimary);
    const std::vector<double>& SphereOriginTotals() const { return sphere_origin_; }

    /// Fold the current event's scattered tally into the run sums (per-EVENT
    /// summation, because steps within one event are correlated -- the same
    /// reason CeeLo's FluenceTally buffers per history).
    void EndSphereEvent();

    void SetSphereEnergy(double e0_keV) { sphere_e0_keV_ = e0_keV; }

    /// Called by EventAction (fluence mode): record one gamma crossing into the
    /// detector with kinetic energy e_keV, event weight w, and whether it is
    /// uncollided (E ~= E0) vs air-scattered (E < E0).
    void RecordFluence(double e_keV, double w, bool uncollided) {
        std::lock_guard<std::mutex> lk(mutex_);
        int bin = static_cast<int>(e_keV / (bin_width_MeV_ * 1000.0));
        if (bin >= 0 && bin < num_bins_) {
            fluence_w_[bin]  += w;
            fluence_w2_[bin] += w * w;
        }
        if (uncollided) { unc_w_ += w; unc_w2_ += w * w; }
        else            { scat_w_ += w; scat_w2_ += w * w; }
    }

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
        num_bins_ = static_cast<int>(
            std::ceil(max_energy_MeV_ / bin_width_MeV_));
    }

    /// Histogram upper bound. FullSpecApprox uses 12 MeV for its requested
    /// domain and expands further when a higher-energy macro is supplied.
    void SetMaxEnergy(double MeV) {
        max_energy_MeV_ = std::max(12.0, MeV);
        num_bins_ = static_cast<int>(
            std::ceil(max_energy_MeV_ / bin_width_MeV_));
    }

    void EnableComponentLabels(bool on) { component_labels_enabled_ = on; }

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
    bool component_labels_enabled_ = false;
    double max_energy_MeV_ = 12.0;
    double bin_width_MeV_ = 0.001; // default 1 keV bins
    int num_bins_ = 12000;
    std::vector<double> histogram_sumw_;
    std::vector<double> histogram_sumw2_;
    static constexpr int kNumSpectrumComponents = 12;
    std::array<std::vector<double>, kNumSpectrumComponents>
        component_histogram_sumw_;
    std::array<std::vector<double>, kNumSpectrumComponents>
        component_histogram_sumw2_;

    // Per-particle-type cumulative energy (MeV), weighted
    double sum_electron_edep_MeV_ = 0.0;
    double sum_gamma_edep_MeV_ = 0.0;
    double sum_total_edep_MeV_ = 0.0;
    double sum_deposit_weights_ = 0.0;
    uint64_t n_events_with_deposit_ = 0;

    double cone_omega_frac_ = 0.0; // 0 = isotropic

    // Fluence-at-detector scorer (air-scatter study). Weighted photon fluence
    // crossing into the detector, binned in incident kinetic energy; per-source-
    // gamma normalization = bin / n_run. Uncollided (peak) vs scattered totals
    // are tracked separately. Uses the same bins as the pulse-height histogram.
    bool fluence_enabled_ = false;
    std::vector<double> fluence_w_;   // weighted counts per energy bin
    std::vector<double> fluence_w2_;  // weighted sum of squares per bin
    double unc_w_ = 0.0, unc_w2_ = 0.0;   // uncollided (E ~= E0) totals
    double scat_w_ = 0.0, scat_w2_ = 0.0; // air-scattered (E < E0) totals

    void WriteFluenceCSV(double energy_keV, uint64_t n_run) const;
    void WriteSphereFluxCSV(uint64_t n_run) const;

    // Environment-scatter sphere flux scorer.
    static constexpr int kSphereNE  = 48;    // log bins, 10 keV .. 3 MeV
    static constexpr int kSphereNMu = 32;    // uniform in mu over [-1, 1]
    static constexpr double kSphereEMin = 10.0;
    static constexpr double kSphereEMax = 3000.0;
    bool sphere_flux_enabled_ = false;
    double sphere_radius_cm_ = 0.0;
    double sphere_e0_keV_ = 0.0;
    std::vector<double> sphere_direct_;   // track length, direct (uncollided)
    std::vector<double> sphere_scat_;     // track length, scattered: run sum
    std::vector<double> sphere_scat_ev_;  // current event's scattered tally
    std::vector<double> sphere_scat_w2_;  // sum of per-event squares
    // Exact integrated-flux uncertainty needs the square of the event's WHOLE
    // scattered sphere track length, not a sum of independent-looking cell
    // variances.  These are raw per-event sufficient statistics in track-
    // length units; WriteSphereFluxCSV applies the deterministic sphere volume.
    uint64_t sphere_scat_total_n_ = 0;
    long double sphere_scat_total_s1_cm_ = 0.0L;
    long double sphere_scat_total_s2_cm2_ = 0.0L;
    std::string sphere_physics_list_;
    double sphere_world_half_extent_m_ = 0.0;
    std::string sphere_primary_sampling_;
    std::string sphere_rng_engine_;
    long sphere_rng_seed_ = 0;
    double sphere_production_cut_mm_ = 0.0;
    std::string sphere_gdml_sha256_;
    std::string sphere_macro_sha256_;
    std::string sphere_geometry_tier_;
    std::string sphere_material_;
    double sphere_distance_cm_ = 0.0;
    double sphere_source_height_cm_ = 0.0;
    double sphere_detector_height_cm_ = 0.0;
    double sphere_gdml_scoring_radius_cm_ = 0.0;
    /// Scattered track length per [origin][energy bin]; mu is summed out since
    /// this exists to size channels, not to resolve their angular shape.
    std::vector<double> sphere_origin_;

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
