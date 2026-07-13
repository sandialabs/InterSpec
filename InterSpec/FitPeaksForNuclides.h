/* InterSpec: an application to analyze spectral gamma radiation data.

 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
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

#ifndef FitPeaksForNuclides_h
#define FitPeaksForNuclides_h

#include <string>
#include <vector>
#include <memory>
#include <functional>

#include "InterSpec/PeakDef.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/PeakFitDetPrefs.h"
#include "InterSpec/DetectorPeakResponse.h"

namespace SpecUtils
{
  class Measurement;
}

namespace FitPeaksForNuclides
{

/** Enable/disable the verbose internal debug trace of `fit_peaks_for_nuclides` (rel-eff fit form/order,
 chi2/dof, gamma-cluster keep/drop decisions, etc.).  For development harnesses ONLY: it is a process-wide
 flag with no synchronization, so it must never be enabled during parallel GA optimization. */
void set_debug_printout( bool enable );

/** an updated implementation of `find_spectroscopic_extent(...)` - we will replace the old implementation after some more testing. */
std::pair<double,double> find_valid_energy_range( const std::shared_ptr<const SpecUtils::Measurement> &meas );


// Internal helpers exposed for unit tests and development harnesses; not part of the public API.
namespace detail
{
  /** A cheap, fit-free estimate of the local continuum: a straight line through averaged channel
   heights at the two edges of a region (the classic two-sideband estimator used for net-area
   determination).  Callers must pad the region beyond any expected peak (~1 FWHM past the
   outermost gamma) so the edge samples measure continuum rather than peak tail.
   */
  struct LocalContinuumEstimate
  {
    double coeffs[2] = { 0.0, 0.0 };   // linear continuum density, relative to reference_energy
    double reference_energy = 0.0;
    bool valid = false;

    // The sideband measurements the line was derived from (windows may have been relocated
    // outward past interfering unfit auto-search peaks; extents are the ones actually used).
    double lower_sideband_counts = 0.0;      // signal-subtracted counts, low-side window
    double upper_sideband_counts = 0.0;
    double lower_sideband_raw_counts = 0.0;  // raw counts (Poisson variance basis)
    double upper_sideband_raw_counts = 0.0;
    double lower_sideband_lo = 0.0, lower_sideband_hi = 0.0;  // low-side window extent, keV
    double upper_sideband_lo = 0.0, upper_sideband_hi = 0.0;  // high-side window extent, keV

    /** Integral of the estimated continuum density over [x0,x1], clamped to >= 0.
     Returns 0 when not valid. */
    double integral( const double x0, const double x1 ) const;

    /** z-score of (low-side density - high-side density) against the Poisson noise of the two
     sideband samples.  Positive when the continuum is higher below the region than above it -
     the signature of a step continuum.  Returns 0 when not valid. */
    double sideband_asymmetry_z() const;
  };//struct LocalContinuumEstimate

  /** Estimate the local continuum from sideband channel averages at `region_lower`/`region_upper`.
   `sideband_num_fwhm` sets how many FWHM of channels are averaged at each edge (minimum 2 channels).

   `predicted_signal`, when supplied, is the expected signal counts over an energy interval
   [x0,x1] (e.g., Gaussian tails of the cluster's own gammas); it is subtracted from each sideband
   sum so peak leakage does not bias the continuum estimate upward.

   `unfit_auto_peaks`, when supplied, veto sideband windows they overlap: the window slides one
   width further from the region (up to 3 tries) to find a clean sample.

   Result is not `valid` if the region is outside the spectrum, degenerate, so close to a spectrum
   edge that a sideband would extend past the first/last channel, or no uncontaminated sideband
   window exists on one side.
   */
  LocalContinuumEstimate estimate_local_continuum(
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const double region_lower,
    const double region_upper,
    const double fwhm,
    const double sideband_num_fwhm,
    const std::function<double(double,double)> &predicted_signal = std::function<double(double,double)>(),
    const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks
      = std::vector<std::shared_ptr<const PeakDef>>() );

  /** Result of the adaptive (data-driven) ROI extent determination. */
  struct AdaptiveExtentResult
  {
    double lower = 0.0, upper = 0.0;      // final ROI bounds
    double sideband_lower_kev = 0.0;      // accepted continuum sideband beyond the core, low side
    double sideband_upper_kev = 0.0;      // accepted continuum sideband beyond the core, high side
  };//struct AdaptiveExtentResult

  /** Determine a ROI's extent by data-driven sideband extension.

   A core region of `core_num_fwhm` x FWHM beyond the outermost expected gammas (plus a fixed skew
   allowance on the low side when `skew_type` is not NoSkew) is always included.  Each side is then
   extended in ~0.375-FWHM blocks while the newly added block stays statistically consistent with a
   linear continuum anchored just inside the already-accepted extent.  `extend_z` is the
   FAMILY-wise consistency z for a full side of extension: the per-block threshold is
   Bonferroni-split across the expected block count, so a genuinely flat continuum has the same
   chance of full extension regardless of how many blocks the cap allows (a fixed per-block z
   would false-stop ~28% of the time at z=2 with ~7 blocks/side).  The block z's denominator
   includes the Poisson noise of the block, the predicted tail leakage of the cluster's own
   gammas, and the estimation variance of the (extrapolated, leveraged) anchor line.  A cumulative
   drift guard catches slow curvature, and any unfit auto-search peak near the block vetoes
   further extension.  Extension stops at `max_num_fwhm` x FWHM beyond the outermost gammas.
   This replaces fixed +/- k x FWHM extents: ROIs shorten automatically next to Compton edges,
   backscatter humps and other peaks, and lengthen over clean flat continua - and the parameters
   are dimensionless, so they transfer across live-times and detector classes.

   @param gamma_energies   Expected gamma energies of the cluster (need not be sorted)
   @param gamma_amplitudes Expected counts of each gamma (parallel to gamma_energies); pass zeros
                           when amplitudes are unknown (tail-leakage prediction is then skipped)
   */
  AdaptiveExtentResult extend_roi_by_sidebands(
    const std::vector<double> &gamma_energies,
    const std::vector<double> &gamma_amplitudes,
    const double effective_fwhm,
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const std::function<double(double)> &fwhm_at_energy,
    const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks,
    const double core_num_fwhm,
    const double extend_z,
    const double max_num_fwhm,
    const PeakDef::SkewType skew_type,
    const double lowest_energy,
    const double highest_energy );

  /** Search for a "clean gap" between two peak groups: a contiguous window at least
   `clean_gap_num_fwhm` x FWHM wide, between the two anchor energies, where the predicted Gaussian
   tail contamination from BOTH groups is statistically negligible compared to the local continuum
   noise, tested at WINDOW level: S_pred over the window / sqrt(B_est over the window)
   < merge_tail_z.  (A former per-~0.25-FWHM-block form understated the window-level contamination
   by ~sqrt(block/window), biasing toward splitting.)  The local continuum estimate has both
   groups' predicted tails subtracted from its sideband samples, so strong close peaks no longer
   inflate B_est and spuriously pass the test.  If no clean window exists, the continuum between
   the peaks cannot be independently anchored and the ROIs should share one continuum (merge).
   This replaces an amplitude-relative tail-fraction merge test that ignored counting statistics:
   a 0.5% tail matters on a high-statistics spectrum and is invisible on a low-statistics one -
   the noise-relative form transfers across live-times.

   Returns true (and the least-contaminated window via the out-params) when a clean gap exists.
   When all amplitudes are zero/unknown the test degenerates to a pure gap-width check.
   */
  bool find_clean_gap_between(
    const std::vector<double> &left_energies,
    const std::vector<double> &left_amplitudes,
    const std::vector<double> &right_energies,
    const std::vector<double> &right_amplitudes,
    const double left_anchor,
    const double right_anchor,
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const std::function<double(double)> &fwhm_at_energy,
    const double merge_tail_z,
    const double clean_gap_num_fwhm,
    double *clean_win_lo,
    double *clean_win_hi );

  /** Choose Linear vs Quadratic continuum for a ROI by AICc over the ROI's continuum sidebands
   (the channels between the ROI bounds and the peak core [core_lo, core_hi], which are excluded).
   Fits both polynomial orders to the sideband channels by Poisson-weighted least squares and picks
   the penalized-chi2 winner; `aicc_penalty` is the kappa scale (2.0 = textbook AIC).  Returns
   Linear when there are too few sideband channels (< 8) to select on, or when curvature is not
   supported by the data.  Replaces a pure ROI-width-in-FWHM rule: whether the continuum actually
   curves is a property of the data, not of the window width.
   */
  PeakContinuum::OffsetType select_continuum_order_by_sidebands(
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const double roi_lower,
    const double roi_upper,
    const double core_lo,
    const double core_hi,
    const double aicc_penalty );
}//namespace detail


// Settings for the gamma clustering algorithm - different values may be used
// for the initial RelActManual stage vs subsequent RelActAuto refinement stages
struct GammaClusteringSettings
{
  double cluster_num_sigma;         // How many sigma to use for clustering gamma lines

  // Minimum Poisson detection significance z = S_est / sqrt(S_est + B_est) to keep a cluster,
  // where S_est is the expected peak counts and B_est the sideband-estimated continuum over the
  // cluster's CORE extent (outermost gammas +/- roi_core_num_fwhm x FWHM - the always-included
  // part of the ROI the fit will see), with the cluster's own predicted tails subtracted from the
  // sideband samples and the samples relocated away from interfering unfit auto-search peaks
  // (see detail::estimate_local_continuum).  Dimensionless, so - unlike the absolute-count
  // gates it replaces - the same value transfers across live-times and detector classes.
  // A fixed (non-configurable) minimum expected-count floor additionally protects the
  // Gaussian-statistics regime; see sm_keep_gate_min_est_counts in FitPeaksForNuclides.cpp.
  double keep_significance_z;

  // Adaptive ROI extent (see detail::extend_roi_by_sidebands): always-included core half-extent
  // beyond the outermost gamma, block-consistency z for data-driven sideband extension, and the
  // extension cap - all in FWHM/z units so they transfer across live-times and detector classes.
  double roi_core_num_fwhm;
  double roi_extend_z;
  double roi_max_num_fwhm;

  // Peak skew the eventual fit will use; extend_roi_by_sidebands adds a low-side core allowance
  // when not NoSkew.  Copied from PeakFitForNuclideConfig::skew_type (not GA-optimized).
  PeakDef::SkewType skew_type = PeakDef::SkewType::NoSkew;

  double max_fwhm_width;            // Maximum ROI width in FWHM before breaking up
  double min_fwhm_roi;              // Minimum ROI width in FWHM to keep

  // kappa for the per-ROI Linear-vs-Quadratic continuum AICc selection
  // (see detail::select_continuum_order_by_sidebands); replaces a width-in-FWHM threshold.
  double cont_order_aicc_penalty = 2.0;

  // Merge decision via the clean-gap test (see detail::find_clean_gap_between): overlapping
  // clusters stay separate only when a continuum-anchoring window exists between their dominant
  // gammas where the predicted tail contamination is < merge_tail_z x sqrt(local continuum).
  double merge_tail_z = 2.0;
  double merge_clean_gap_fwhm = 1.0;  // required clean-window width, in FWHM

  // Parameters for synthetic spectrum-based ROI breaking
  // Region around minimum/maximum to compute significance (in FWHM units)
  double break_check_fwhm_fraction = 0.5;

  // Threshold for considering a peak "significant" between breakpoints (sigma)
  // Must have a peak exceeding this between any two breakpoints (or ROI edge and breakpoint)
  double break_peak_significance_threshold = 2.0;

  // When multiple breakpoint candidates have similar significance, use depth_score as tiebreaker
  // This defines "similar" - candidates within this sigma of each other are considered tied
  double break_significance_tie_threshold = 0.5;

  // Step continuum decision thresholds
  // Minimum peak detection significance z = S_est / sqrt(S_est + B_est), with B_est from the
  // sideband continuum estimate, to consider a step continuum (a step only matters when the peak
  // towers over the continuum, which is inherently a significance statement - an absolute-count
  // gate would not transfer across live-times).
  double step_cont_min_peak_significance = 30.0;
  // Chi2 margin by which the step-continuum trial fit must beat the polynomial fit for the ROI to
  // get a step continuum (the trial pairs equal-parameter-count candidates - Linear vs FlatStep,
  // Quadratic vs LinearStep - so the AICc penalty terms cancel and the decision reduces to a
  // chi2 comparison with this tunable bias against the step).  Replaces the former left-vs-right
  // probe-window nsigma test, which self-vetoed on tight ROIs and read neighbor peaks as steps.
  double step_trial_chi2_margin = 4.0;
};


// Result from RelActAuto peak fitting
struct PeakFitResult
{
  RelActCalcAuto::RelActAutoSolution::Status status;
  std::string error_message;
  std::vector<std::string> warnings;  // warnings that don't prevent success

  // Peaks after combining overlapping peaks within ROIs.
  // Peaks that are close together (within 1.5 sigma) or where a smaller peak
  // is not statistically distinguishable from a larger peak's tail are merged
  // into a single peak with combined properties.
  // Note: even if fitting energy calibration was selected, these peaks are in the original spectrums energy cal.
  std::vector<PeakDef> fit_peaks;

  // Original uncombined peaks from the fit - preserves all individual peak information.
  // This is the raw output from RelActAuto before peak combination.
  // Note: even if fitting energy calibration was selected, these peaks are in the original spectrums energy cal.
  std::vector<PeakDef> uncombined_fit_peaks;

  // Peaks that are statistically observable in the spectrum.
  // Computed from fit_peaks by:
  // 1. Removing peaks with initial significance < threshold (using raw data area)
  // 2. Refitting each ROI with PeakFitLM::refitPeaksThatShareROI_LM (SmallRefinementOnly)
  // 3. Iteratively removing peaks with final significance < threshold and refitting
  // Note: even if fitting energy calibration was selected, these peaks are in the original spectrums energy cal.
  std::vector<PeakDef> observable_peaks;

  // Existing user peaks that should be removed from PeakModel when the result is accepted.
  // Behaviour depends on which option flags were set:
  //   - DoNotUseExistingRois: always empty (existing ROIs and their peaks are left untouched).
  //   - ExistingPeaksAsFreePeak: peaks that had a matching source (unconditionally replaced),
  //     plus bystander peaks whose ROI ended up in the solution.
  //   - Default (neither flag): user peaks whose source matches a fit source and whose energy
  //     falls within a fitted observable ROI (i.e. the fit replaced them).
  // These are the exact shared_ptr instances from the input user_peaks vector, so pointer
  // identity is preserved for PeakModel::removePeaks().
  std::vector<std::shared_ptr<const PeakDef>> original_peaks_to_remove;

  // The final RelActCalcAuto solution.  Note that `solution.m_foreground` and/or `solution.m_foreground`
  // may have a different energy calibration than input foreground/background, due to iteratively finding solution.
  // This means the various peak quantities of solution that are supposed to be in the original energy calibration,
  // would need translating (i.e., with `EnergyCal::translatePeaksForCalibrationChange(...)`) before using.
  RelActCalcAuto::RelActAutoSolution solution;  // only valid if status == Success
};//struct PeakFitResult


// Configuration parameters for peak fitting for nuclides
// These values will be optimized via genetic algorithm for different detector types
struct PeakFitForNuclideConfig
{
  /** Returns the default `PeakFitForNuclideConfig` for the detector type.
   */
  static const PeakFitForNuclideConfig &default_config( const PeakFitUtils::CoarseResolutionType det_type );


  // FWHM functional form
  DetectorPeakResponse::ResolutionFnctForm fwhm_functional_form = DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;

  // RelActManual parameters for initial relative efficiency estimation
  double rel_eff_manual_base_rel_eff_uncert = 0.0;
  double initial_nuc_match_cluster_num_sigma = 1.5;
  double manual_eff_cluster_num_sigma = 1.5;

  // The RelActManual equation form and order are selected per spectrum by an AICc ladder
  // (see manual_releff_aicc_penalty below), not configured here.

  // ROI clustering thresholds for manual RelEff stage.
  // Cluster keep-gate: minimum z = S_est/sqrt(S_est + B_est) (see GammaClusteringSettings::keep_significance_z).
  double manual_keep_significance_z = 2.0;
  double manual_rel_eff_sol_min_fwhm_roi = 1.0;
  double manual_rel_eff_sol_max_fwhm = 15.0;
  double manual_roi_core_num_fwhm = 1.5;  // Always-included ROI half-extent beyond outermost gamma

  // RelActAuto parameters
  RelActCalcAuto::FwhmForm fwhm_form = RelActCalcAuto::FwhmForm::Berstein_3;
  double rel_eff_auto_base_rel_eff_uncert = 0.1;  // BR uncertainty for RelActAuto

  // ROI clustering thresholds for auto RelEff refinement stage
  double auto_rel_eff_cluster_num_sigma = 2.0;  // Slightly wider clustering with better rel-eff
  double auto_keep_significance_z = 3.0;  // Cluster keep-gate z (see manual_keep_significance_z)
  double auto_roi_core_num_fwhm = 1.5;    // Always-included ROI half-extent beyond outermost gamma
  double auto_rel_eff_sol_max_fwhm = 12.0;  // Tighter constraint as solution improves
  double auto_rel_eff_sol_min_fwhm_roi = 1.25;

  // Adaptive ROI sideband extension (shared between manual and auto stages): block-consistency z
  // and extension cap in FWHM beyond the outermost gamma.  See detail::extend_roi_by_sidebands.
  double roi_extend_z = 2.0;
  double roi_max_num_fwhm = 4.0;

  // ROI merge/split decision (shared between stages): overlapping clusters stay separate only
  // when a clean continuum-anchoring gap exists between them.  See detail::find_clean_gap_between.
  double merge_tail_z = 2.0;
  double merge_clean_gap_fwhm = 1.0;

  // AICc complexity-penalty scales (kappa; 2.0 = textbook AIC) for the per-spectrum manual
  // rel-eff form/order ladder and the per-ROI continuum-order selection.  These two scalars
  // replace the former per-peak-count form/order genes and the width-threshold quadratic rule;
  // they also absorb the non-textbook scale of the judged/sideband chi2s they penalize.
  double manual_releff_aicc_penalty = 2.0;
  double cont_order_aicc_penalty = 2.0;


  // RelActAuto relative efficiency model parameters
  RelActCalc::RelEffEqnForm rel_eff_eqn_type = RelActCalc::RelEffEqnForm::LnXLnY;
  size_t rel_eff_eqn_order = 2;

  // Physical model shielding (empty vectors mean no shielding)
  std::vector<std::shared_ptr<RelActCalc::PhysicalModelShieldInput>> phys_model_self_atten;
  std::vector<std::shared_ptr<RelActCalc::PhysicalModelShieldInput>> phys_model_external_atten;

  // Desperation Physical Model shielding configuration
  // Atomic number for desperation external shielding (0 = don't use, 1-98 = valid elements)
  double desperation_phys_model_atomic_number = 26.0;  // Default: iron

  // Areal density starting value for desperation external shielding (in g/cm2)
  double desperation_phys_model_areal_density_g_per_cm2 = 1.0;  // Default: 1.0 g/cm2

  bool nucs_of_el_same_age = true;
  
  // Physical model options (only used when rel_eff_eqn_type == FramPhysicalModel)
  bool phys_model_use_hoerl = true;

  // Fields for RelActAuto options configuration - this is only used for optimiziation work - it is supersceeded by `FitSrcPeaksOptions::DoNotVaryEnergyCal
  bool fit_energy_cal = true;

  // Maximum acceptable chi2/dof for the initial RelActCalcManual rel-eff solution
  // (after retry, if any).  When the matched-peaks fit can't reach this quality the
  // source(s) are likely not actually present in the spectrum; the manual rel-eff
  // step is treated as failed and the caller falls back to `estimate_initial_rois_fallback`.
  // Set to a large value (e.g. 1e6) to effectively disable the cap.
  double initial_manual_rel_eff_max_chi2_dof = 25.0;

  // ROI significance threshold for iterative refinement, as an equivalent-z of a
  // likelihood-ratio test (Wilks): the chi2 improvement of adding the ROI's peaks over a
  // quadratic continuum-only null, referred to a chi2 distribution with dof = number of
  // peaks, converted to a normal quantile.  Replaces the former trio of thresholds
  // (min chi2-reduction, min single-peak significance, quad-continuum gate), which were
  // redundant parameterizations of this one test: a strong single peak gives a huge
  // delta-chi2, and the quadratic null already absorbs smooth continuum curvature.
  double roi_significance_z = 3.0;

  // Threshold for initial peak significance filter before refitting for observable_peaks.
  // Significance = S / sqrt(S + B) with S = peak_amplitude * 0.7607 (fraction of Gaussian
  // area within +/-1 FWHM) and B the FITTED continuum integral over the same +/-1 FWHM
  // window - using the fitted continuum rather than the gross data means a neighboring
  // peak's counts no longer dilute this peak's significance.
  double observable_peak_initial_significance_threshold = 2.25;

  // Threshold for final peak significance after refitting for observable_peaks.
  // Significance = peak_area / peak_area_uncertainty
  double observable_peak_final_significance_threshold = 2.0;

  // Step continuum decision parameters
  // Minimum peak detection significance z = S_est/sqrt(S_est + B_est) to consider step continuum
  double step_cont_min_peak_significance = 30.0;
  // Chi2 margin the step trial fit must beat the polynomial fit by (see
  // GammaClusteringSettings::step_trial_chi2_margin for the full description).
  double step_trial_chi2_margin = 4.0;

  // Peak skew type to apply during the RelActAuto fit - note this parameter should not be optimized, but rather
  //  something that might be over-rided according to the detector efficiency or user preferences.
  PeakDef::SkewType skew_type = PeakDef::SkewType::NoSkew;

  /** CSS color string for NORM background nuclide peaks (used when FitNormBkgrndPeaks is set).
   Non-empty default ensures Rel. Eff. chart data points render; override from the background
   ReferenceLineInfo color or ColorTheme at the call site. */
  std::string norm_css_color = "rgb(150,150,150)";


  // Get GammaClusteringSettings for manual RelEff stage
  GammaClusteringSettings get_manual_clustering_settings() const
  {
    GammaClusteringSettings settings;
    settings.cluster_num_sigma = manual_eff_cluster_num_sigma;
    settings.keep_significance_z = manual_keep_significance_z;
    settings.roi_core_num_fwhm = manual_roi_core_num_fwhm;
    settings.roi_extend_z = roi_extend_z;
    settings.roi_max_num_fwhm = roi_max_num_fwhm;
    settings.skew_type = skew_type;
    settings.max_fwhm_width = manual_rel_eff_sol_max_fwhm;
    settings.min_fwhm_roi = manual_rel_eff_sol_min_fwhm_roi;
    settings.cont_order_aicc_penalty = cont_order_aicc_penalty;
    settings.merge_tail_z = merge_tail_z;
    settings.merge_clean_gap_fwhm = merge_clean_gap_fwhm;
    settings.step_cont_min_peak_significance = step_cont_min_peak_significance;
    settings.step_trial_chi2_margin = step_trial_chi2_margin;
    return settings;
  }

  // Get GammaClusteringSettings for auto refinement stage
  GammaClusteringSettings get_auto_clustering_settings() const
  {
    GammaClusteringSettings settings;
    settings.cluster_num_sigma = auto_rel_eff_cluster_num_sigma;
    settings.keep_significance_z = auto_keep_significance_z;
    settings.roi_core_num_fwhm = auto_roi_core_num_fwhm;
    settings.roi_extend_z = roi_extend_z;
    settings.roi_max_num_fwhm = roi_max_num_fwhm;
    settings.skew_type = skew_type;
    settings.max_fwhm_width = auto_rel_eff_sol_max_fwhm;
    settings.min_fwhm_roi = auto_rel_eff_sol_min_fwhm_roi;
    settings.cont_order_aicc_penalty = cont_order_aicc_penalty;
    settings.merge_tail_z = merge_tail_z;
    settings.merge_clean_gap_fwhm = merge_clean_gap_fwhm;
    settings.step_cont_min_peak_significance = step_cont_min_peak_significance;
    settings.step_trial_chi2_margin = step_trial_chi2_margin;
    return settings;
  }

private:
  PeakFitForNuclideConfig(){};
};//struct PeakFitForNuclideConfig


/** Returns true if the source should be excluded from the peak_fit_improve
 GA's background-false-positive penalty.

 True for the canonical NORM nuclides (K40, Ra226, U235, U238, Th232), any
 nuclide whose decay-chain ancestors include U235, U238, or Th232, and a
 hand-curated extras list (initially U232, U233, F18) maintained alongside
 the implementation.

 False for elements (xrays) and reactions - they have no decay chain to
 test against.  Adjust the implementation if specific elements/reactions
 ever need exclusion.
 */
bool is_norm_like_for_ga( const RelActCalcAuto::SrcVariant &src );

/** Returns true if `energy_kev` is within `tolerance_kev` of any commonly-
 observed NORM gamma or NORM-element K-xray line.  Used by the GA's
 background-fit penalty to suppress false positives where a source's
 gamma happens to land on a real NORM peak in the background. */
bool is_near_strong_norm_gamma( double energy_kev, double tolerance_kev );


/** Options for fitting the peaks of nuclides.

 By default, when No FitSrcPeaksOptions are specified:
 - Existing ROIs containing only other-source peaks will not have peaks of the source(s) being fit added to them; if the photopeak of the source(s) being fit are inside an existing ROI, that photopeak will be ignored. If a photopeak is adjacent to an existing ROI, the existing ROI will not be modified, but the added ROI may slightly overlap (although it will be tried to make it not overlap - but if PeakFitForNuclideConfig::auto_rel_eff_sol_min_fwhm_roi dictates it must, then it will), but will not extend any closer than `0.5*PeakFitForNuclideConfig::auto_rel_eff_sol_min_fwhm_roi` to any peak mean in the existing ROI (e.g., new ROI wont cover the mean of any peak in the existing ROI).
 - If a peak/ROI for a source being fit (and only contains peaks for the source(s) being fit) is already present in the data, then these ROIs will be re-fit; their energy range and/or peak properties may become different.  If the fit to the source(s) determines that the ROI should not be present (i.e., it doesnt think the peak(s) are significant) then the original ROI will remain unaltered.
 - Existing ROIs containing both a source being fit, as well as other-source peaks (mixed ROIs) use the existing ROI bounds; other-source peaks are included as bystander floating peaks in the combined fit.  These bystander peaks will be included in the results (with the original bystander peak being in the collection of peaks that should be removed before adding the result peaks), and have the same sources associated with them.  It is possible the bystander peak will become insignificant in the fit, so it may just be in the collection of peaks to remove without a replacement, but the ROI will maintain the original bounds (even if the bystander peak disappears).
 - ROIs added for the sources being fit, will not overlap with each other - they will have at least one channel between them.
 - All peaks sharing a PeakContinuum with a removed peak are also removed and replaced together.
 
 When the DoNotUseExistingRois option is specified:
 - Any existing ROI will not be used, and any gammas from the current source(s) that fall within the ROI will not be considered, even if that ROI has a peak with a source being fitted (i.e., existing peaks/ROIs of the sources being fit will not be altered in any way).  Existing ROIs will not be combined with new ROIS. Like default, a new ROI may slightly overlap with an existing ROI, but not any closer than `0.5*PeakFitForNuclideConfig::auto_rel_eff_sol_min_fwhm_roi` to any peak mean in the existing ROI.
 
 When the ExistingPeaksAsFreePeak option is specified:
 - Potential peaks for the sources being fit, that are adjacent to, or within an existing ROI may be combined with the existing ROI (potentially, and likely altering its energy extent); the existing peaks will be treated as freely-floating peaks, and included in the results (and the set of peaks to delete).  If a new peak for a source being fit is not added to (or combined with) an existing ROI, the existing ROI will not be altered (either its energy range, or the peaks in it).
 
 DoNotUseExistingRois can not be used in combination with ExistingPeaksAsFreePeak.
 */
enum FitSrcPeaksOptions
{
  /** With this option, any existing ROI will not be used, and any gammas from the current source(s)
   that fall within the ROI will not be considered.  Without this option the peaks you have already
   fit, for the sources you are currently trying to fit peaks of, will be replaced.
   
   TODO: we should probably rename DoNotUseExistingRois to DoNotUseExistingRoisOfSourcesBeingFit
   */
  DoNotUseExistingRois = 0x01,

  /** Normally ROIs of source peak will try to be limited in energy range to mitigate effects of
   other nearby peaks (of sources you are not fitting); with this option, the nearby peaks may
   share the ROI will be left in as freely-floating peaks, and included in the results.
   */
  ExistingPeaksAsFreePeak = 0x02,
  
  /** The energy calibration stays fixed. */
  DoNotVaryEnergyCal = 0x04,
  
  /** The energy calibration is varied in the fit, but the ROI extent is not updated based on the fit calibration.

   Note: returned peaks are still in the original spectrum's energy calibration (they are built from
   the solve's spectrum-cal peak set, and the working foreground is never cal-advanced in this mode);
   the fitted calibration adjustment itself is simply discarded.
   */
  DoNotRefineEnergyCal = 0x08,
  
  /** Fit the NORM peaks, using a second relative efficiency curve. */
  FitNormBkgrndPeaks = 0x10,
  
  /** Fit the NORM peaks to assist in getting FWHM functional form and energy calibration right,
   but wont return in the solution peaks.
   */
  FitNormBkgrndPeaksDontUse = 0x20,
};

/** Function to fit all the observable peaks for one or more sources.

 This function performs the complete peak fitting workflow:
 1. Determines FWHM functional form from auto-search peaks or DRF
 2. Matches auto-search peaks to source nuclides
 3. Performs RelActManual fit to get initial relative efficiency estimate
 4. Clusters gamma lines into ROIs based on manual rel-eff
 5. Calls fit_peaks_for_nuclide_relactauto with the single configured RelEff curve type (config.rel_eff_eqn_type)
 6. Optionally retries with a Physical-Model "desperation" configuration if the initial fit is poor,
    keeping whichever result has the lower chi2/dof

 @param auto_search_peaks Initial peaks found by automatic search (user peaks merged with auto-detected)
 @param foreground Foreground spectrum to fit
 @param sources Vector of source nuclides to fit
 @param user_peaks The user's current peaks from PeakModel.  Used when DoNotUseExistingRois or
        ExistingPeaksAsFreePeak options are set to identify existing ROIs / add floating peaks.
        May be empty if neither option is set.
 @param background Background spectrum (can be nullptr)
 @param drf_input Detector response function (can be nullptr, will use generic if needed)
 @param options Options for how the fit should be done.
 @param config Configuration for peak fitting parameters
 @param peak_fit_prefs Peak fitting preferences (detector type, skew, FWHM method).
        The isHPGe flag is derived from prefs->m_det_type.  If nullptr, defaults are used.
 @return PeakFitResult with status, error message, fit peaks, and solution
 */

PeakFitResult fit_peaks_for_nuclides(
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::vector<RelActCalcAuto::NucInputInfo> &sources,
  const std::vector<std::shared_ptr<const PeakDef>> &user_peaks,
  std::shared_ptr<const SpecUtils::Measurement> background,
  const std::shared_ptr<const DetectorPeakResponse> &drf_input,
  const Wt::WFlags<FitSrcPeaksOptions> options,
  const PeakFitForNuclideConfig &config,
  const std::shared_ptr<const PeakFitDetPrefs> &peak_fit_prefs );

PeakFitResult fit_peaks_for_nuclides(
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::vector<RelActCalcAuto::SrcVariant> &sources,
  const std::vector<std::shared_ptr<const PeakDef>> &user_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &background,
  const std::shared_ptr<const DetectorPeakResponse> &drf_input,
  const Wt::WFlags<FitSrcPeaksOptions> options,
  const PeakFitForNuclideConfig &config,
  const std::shared_ptr<const PeakFitDetPrefs> &peak_fit_prefs );
  

// Helper function for estimating initial ROIs when no peaks are available
std::vector<RelActCalcAuto::RoiRange> estimate_initial_rois_without_peaks(
  const std::vector<RelActCalcAuto::NucInputInfo> &sources,
  const std::shared_ptr<const DetectorPeakResponse> &drf,
  const PeakFitUtils::CoarseResolutionType det_type,
  const DetectorPeakResponse::ResolutionFnctForm fwhmFnctnlForm,
  const std::vector<float> &fwhm_coefficients,
  const double lower_fwhm_energy,
  const double upper_fwhm_energy,
  const double min_valid_energy,
  const double max_valid_energy,
  const PeakFitForNuclideConfig &config );

/** Debug harness for fit_peaks_for_nuclides.
 Loads spectrum files, runs auto peak search, calls fit_peaks_for_nuclides,
 and prints diagnostics at each stage. Enable PERFORM_DEVELOPER_CHECKS for
 full internal trace output via the existing local_debug_printout flag.
 Called from main.cpp before server startup.
*/
int debug_fit_peaks_for_nuclides();

}//namespace FitPeaksForNuclides

#endif //FitPeaksForNuclides_h
