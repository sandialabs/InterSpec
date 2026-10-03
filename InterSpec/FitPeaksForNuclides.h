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

#include <set>
#include <string>
#include <vector>
#include <atomic>
#include <memory>
#include <cstdint>
#include <limits>
#include <optional>
#include <functional>

#include "SandiaDecay/SandiaDecay.h"

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

struct GammaClusteringSettings;  // defined below; the detail:: planner test seam takes it by reference

/** Outcome of one automatic boundary decision.  These values describe structural policy, not a
 fitted-peak quality classification. */
enum class AutomaticRoiDecision
{
  KeepSeparate,
  MergeInseparable,
  MergeInseparableWide,
  UnmodeledFeatureBlocked,
  ProtectedGeometry,
  R6LegacyBypass,
  /** A provisionally admitted source group in an over-wide overlap component was retained by the
   measured-data continuum/source/free-feature comparison. */
  SourceBridgeRetained,
  /** The local continuum-only model explained a provisional source group; it was rejected before
   atom admission, so it cannot bridge otherwise distinct ROIs. */
  SourceBridgeRejectedContinuum,
  /** A free data-peak explanation beat the requested-source-tied model; the provisional source
   group was rejected before atom admission and the found feature remains unmodeled evidence. */
  SourceBridgeRejectedFreeFeature,
  /** No core-safe partition existed, so the atom-safe layer retained the incumbent geometry
   (merged the pair or left it unchanged) rather than dropping a requested line. */
  InfeasiblePartition
};

/** Concise, reporter-ready evidence for an automatic ROI join/partition decision. */
struct AutomaticRoiDecisionDiagnostic
{
  AutomaticRoiDecision decision = AutomaticRoiDecision::KeepSeparate;
  std::string stage;
  std::string reason;
  double left_lower = 0.0;
  double left_upper = 0.0;
  double right_lower = 0.0;
  double right_upper = 0.0;
  double boundary_energy = 0.0;
  size_t boundary_channel = 0;
  size_t calibration_num_channels = 0;
  double separation_fwhm = 0.0;
  double observed_valley_counts = 0.0;
  double snip_valley_counts = 0.0;
  double modeled_tail_counts = 0.0;
  double modeled_tail_significance = 0.0;
  double unexplained_excess_significance = 0.0;
  double snip_mismatch_significance = 0.0;
  size_t left_sideband_channels = 0;
  size_t right_sideband_channels = 0;
  bool sidebands_adequate = false;
  bool unmodeled_core_blocked = false;
  bool used_global_continuum = false;
  double combined_width_fwhm = 0.0;
  double width_pressure = 0.0;
  double one_roi_aicc = 0.0;
  double two_roi_aicc = 0.0;
  /** Local H0/Hs/Hf evidence values for a provisional source group in an over-wide component.
   Unavailable hypotheses are NaN; these fields are meaningful only when
   `source_evidence_tested` is true. */
  bool source_evidence_tested = false;
  double source_null_aicc = 0.0;
  double source_tied_aicc = 0.0;
  double free_feature_aicc = 0.0;
  double source_likelihood_z = 0.0;
  double free_feature_energy = 0.0;
  /** Number of admitted atoms that the atom-safe partition assigned to a different child than
   their original group membership (spatial reassignment).  Purely informational. */
  size_t atoms_reassigned = 0;
  /** True when the atom-safe layer could not find a core-safe partition and retained the
   incumbent geometry (see AutomaticRoiDecision::InfeasiblePartition). */
  bool partition_infeasible = false;
};

const char *automatic_roi_decision_name( AutomaticRoiDecision decision );

/** Enable/disable the verbose internal debug trace of `fit_peaks_for_nuclides` (rel-eff fit form/order,
 chi2/dof, gamma-cluster keep/drop decisions, etc.).  For development harnesses ONLY: it is a process-wide
 flag with no synchronization, so it must never be enabled during parallel GA optimization. */
void set_debug_printout( bool enable );

/** Developer invariant checks (the FITPEAKS_DEV_CHECK macro in FitPeaksForNuclides.cpp) stay compiled
 in Release builds whenever PERFORM_DEVELOPER_CHECKS is on.  A failed check is always printed to stderr
 and appended to the fit result's `warnings` (prefixed "DevCheck:"); when this flag is set (development
 harnesses) it additionally throws std::logic_error, so the failure is scored as a mechanical failure
 instead of silently continuing.  Default false (GUI).  Process-wide and thread-safe. */
void set_dev_checks_throw( bool enable );
bool dev_checks_throw();

/** an updated implementation of `find_spectroscopic_extent(...)` - we will replace the old implementation after some more testing. */
std::pair<double,double> find_valid_energy_range( const std::shared_ptr<const SpecUtils::Measurement> &meas );


// Internal helpers exposed for unit tests and development harnesses; not part of the public API.
namespace detail
{
  // Forward declaration (full definition is GlobalContinuumEstimate, further below) so
  // find_clean_gap_between can accept an optional pointer to the shared global continuum.
  struct GlobalContinuumEstimate;

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
    /** (low-side density - high-side density) / low-side density: how big the step is relative to
     the continuum it sits on.  A step worth modelling is one a fit would visibly bend around; a
     0.02 step under a peak on a large foreign continuum is not, however significant it is.
     Returns 0 when not valid or when the low-side density is not positive. */
    double sideband_step_fraction() const;
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
      = std::vector<std::shared_ptr<const PeakDef>>() ,
    const double max_predicted_fraction = 0.0 );

  /** The chi2 `peaks` would leave against a quadratic continuum-only null over the channels starting
   at `energies` (lower edges, `nchannel + 1` of them) if they described `counts` exactly: their
   summed shape fit by a quadratic alone, weighted by max(counts, 1).  Measures how much a
   peaks-vs-quadratic test on these channels can see of the peaks - little in a narrow ROI, where the
   quadratic absorbs most of a peak.  Returns 0 on failure or with fewer than 4 channels.
   */
  double quadratic_null_power( const float * const energies,
                               const float * const counts,
                               const size_t nchannel,
                               const double ref_energy,
                               const std::vector<std::shared_ptr<const PeakDef>> &peaks );

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

  /** Search for a statistically unbridged boundary between two peak groups: a window at least
   `clean_gap_num_fwhm` x FWHM wide, between the two anchor energies, where the predicted Gaussian
   tail contamination from BOTH groups is statistically negligible compared to the local continuum
   noise, tested at WINDOW level: S_pred over the window / sqrt(B_est over the window)
   < merge_tail_z.  (A former per-~0.25-FWHM-block form understated the window-level contamination
   by ~sqrt(block/window), biasing toward splitting.)  The local continuum estimate has both
   groups' predicted tails subtracted from its sideband samples, so strong close peaks no longer
   inflate B_est and spuriously pass the test.  The eventual boundary decision also rejects an
   unexplained peak-like excess over the shared continuum.  Thus this is positive evidence for a
   lack of statistically significant peak content connecting the groups, not a requirement that
   the noisy raw spectrum contain a morphological local minimum.  If no such window exists, the
   continuum between the peaks cannot be independently anchored and the ROIs should share one
   continuum (merge).
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
    double *clean_win_hi,
    const GlobalContinuumEstimate *global_continuum = nullptr );

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
    const double aicc_penalty,
    const std::function<double(double,double)> &predicted_signal = nullptr,
    std::string *decision_note = nullptr );


  /** A single SNIP-based continuum estimate over the whole valid spectroscopic extent, shared by the
   clustering/gating decisions so they all reason about the SAME B(E) instead of each re-estimating a
   local two-sideband line (which is unreliable under broad low-resolution peaks on a structured
   Compton continuum).  Built once (see make_global_continuum) from the FWHM-window SNIP with
   per-detector-class parameters.  Every consumer MUST fall back to its prior local estimate when
   `valid()` is false, so an invalid provider reproduces the pre-R1-step2 behaviour exactly. */
  struct GlobalContinuumEstimate
  {
    std::shared_ptr<const SpecUtils::Measurement> snip;        // SNIP continuum (foreground binning)
    std::shared_ptr<const SpecUtils::Measurement> foreground;  // the data (for the variance bound)
    bool built = false;

    bool valid() const { return built && snip && foreground; }

    /** Integral of the SNIP continuum over [x0,x1] (counts), clamped >= 0; 0 if invalid or x1<=x0. */
    double integral( double x0, double x1 ) const;

    /** SNIP continuum density (counts/keV) at energy E; 0 if invalid. */
    double density_at( double E ) const;

    /** A conservative Poisson variance of the continuum over [x0,x1]: the DATA counts there (an upper
     bound, largest exactly where peaks make the SNIP least trustworthy).  0 if invalid. */
    double integral_variance( double x0, double x1 ) const;
  };

  /** Build a GlobalContinuumEstimate from `foreground` with the FWHM-window SNIP restricted to
   [restrict_lower_energy, restrict_upper_energy] (the valid spectroscopic extent).  SNIP parameters
   are selected by detector class: HPGe = 2.0xFWHM / order 2 / 3-ch presmooth / LLS on; else
   (NaI/LaBr/CZT) = 1.5xFWHM / order 2 / 7-ch presmooth / LLS off.  Returns an invalid estimate
   (`valid()==false`) on any failure, so callers transparently fall back to local estimation. */
  GlobalContinuumEstimate make_global_continuum(
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const std::function<double(double)> &fwhm_at_energy,
    PeakFitUtils::CoarseResolutionType det_type,
    double restrict_lower_energy,
    double restrict_upper_energy );

  /** Source-clean seed recovery is warranted only when at least two independent, significant
   source anchors would be lost by the current predicted-count keep gate. */

  /** Transactional source-clean acceptance: recover predicted anchors, preserve the incumbent's
   FWHM-distinct significant fitted source evidence, and improve the filtered data score.  A valid
   candidate may replace an incumbent for which no significant-ROI score was available. */

  /** Small-sample-corrected data-only information criterion used for common-channel challengers. */
  double data_only_aicc( double data_chi2, size_t num_data_rows,
                         size_t num_parameters, double parameter_penalty );

  /** The modeled content and proposed bounds on one side of an automatic ROI boundary. */
  struct AutomaticRoiGroup
  {
    double lower = 0.0;
    double upper = 0.0;
    std::vector<double> peak_energies;
    std::vector<double> peak_areas;
    size_t joined_groups = 1;
    bool protected_geometry = false;
  };

  struct AutomaticRoiPolicySettings
  {
    double merge_tail_z = 2.0;
    double merge_clean_gap_fwhm = 1.0;
    double continuum_aicc_penalty = 2.0;
    double peak_core_num_fwhm = 1.0;
    double max_width_fwhm = 12.0;
    // Optional admission rail for a proposed split.  Unlike force_partition_gap_fwhm, this
    // filters even AICc-preferred boundaries so dense resolved multiplets do not fan out merely
    // because two continua have more flexibility.
    double minimum_partition_gap_fwhm = 0.0;
    // Permit a tail-subtracted continuum-anchoring window to override the hard modeled-core
    // gap rail.  This is opt-in because core spacing alone cannot distinguish a visible valley
    // from a dense resolved multiplet, while the clean-window test can.
    bool allow_clean_gap_partition_override = false;
    // A final sparse-ROI challenger may instead use a residual valley: no statistically
    // significant unmodeled excess above the shared SNIP continuum and modeled Gaussian tails.
    // Zero disables it, preserving the stricter tail-clean window rule.
    double residual_valley_max_excess_z = 0.0;
    // Maximum children from one measured whole-component partition.  Two preserves the original
    // binary challenger; larger values are admitted only when clean-gap mode is explicitly on.
    size_t max_partition_children = 2;
    const GlobalContinuumEstimate *global_continuum = nullptr;
    double force_partition_gap_fwhm = 0.0;
    /** Permit an over-wide recovery component to continue past the ordinary overlapping-core
     short-circuit and seek a scored boundary.  False preserves the initial atom-safe policy. */
    bool allow_overwide_overlap_partition = false;
    std::string stage;
  };

  struct AutomaticRoiPolicyResult
  {
    AutomaticRoiDecision decision = AutomaticRoiDecision::KeepSeparate;
    double boundary_energy = 0.0;
    double exclusion_lower = 0.0;
    double exclusion_upper = 0.0;
    AutomaticRoiDecisionDiagnostic diagnostic;
  };

  enum class SourceClusterEvidenceDecision
  {
    RetainSource,
    RejectContinuumOnly,
    RejectFreeFeature,
    InsufficientEvidence
  };

  /** Result of the local, common-channel H0/Hs/Hf comparison used only for provisional source
   groups inside over-wide transitive overlap components.  H0 is continuum-only; Hs adds one
   locally scaled, fixed-ratio mixture of the requested source lines; Hf jointly adds every
   FWHM-distinct significant found peak outside all requested-line cores. */
  struct SourceClusterEvidenceResult
  {
    SourceClusterEvidenceDecision decision = SourceClusterEvidenceDecision::InsufficientEvidence;
    double null_aicc = std::numeric_limits<double>::quiet_NaN();
    double source_aicc = std::numeric_limits<double>::quiet_NaN();
    double free_feature_aicc = std::numeric_limits<double>::quiet_NaN();
    double source_likelihood_z = 0.0;
    double free_feature_energy = 0.0;
    std::string reason;
  };

  /** Transactionally classify one provisional source cluster on identical measured-data channels.
   No source names or shielding labels participate.  An unavailable/ill-conditioned comparison is
   conservative (`InsufficientEvidence`) and callers retain the provisional source group. */
  SourceClusterEvidenceResult evaluate_source_cluster_evidence(
    const std::vector<double> &source_energies,
    const std::vector<double> &source_areas,
    double lower_energy,
    double upper_energy,
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const std::function<double(double)> &fwhm_at_energy,
    const std::vector<std::shared_ptr<const PeakDef>> &found_peaks,
    double significance_z,
    double source_core_num_fwhm,
    double aicc_penalty );

  /** Decide whether adjacent automatic groups may share an ROI.  All statistical comparisons use
   the same foreground channels and the shared current-calibration SNIP estimate. */
  AutomaticRoiPolicyResult evaluate_automatic_roi_boundary(
    const AutomaticRoiGroup &left,
    const AutomaticRoiGroup &right,
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const GlobalContinuumEstimate *global_continuum,
    const std::function<double(double)> &fwhm_at_energy,
    const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks,
    const AutomaticRoiPolicySettings &settings );


  //=========================================================================================
  // Atom-safe automatic ROI partition layer.
  //
  // Every automatic (policy-mode) ROI split/combine operates on stable "atoms" - one per
  // admitted modeled gamma line (or line-like evidence) - carried WITH the ROI geometry rather
  // than reconstructed from a flat energy list by geometric containment.  The layer guarantees,
  // for each operation, that every admitted atom is represented exactly once afterward, no atom
  // is lost/duplicated/silently reassigned, each atom's core lies within its assigned ROI, the
  // resulting automatic ROIs are channel-disjoint, protected user/mixed geometry is untouched,
  // and unmodeled-exclusion regions are never split or merged through.  When no core-safe
  // partition exists it retains the incumbent geometry (merge or unchanged) instead of dropping
  // a side.  `use_automatic_roi_policy == false` (R6 legacy) paths never enter this layer.
  //=========================================================================================

  /** Kind of evidence an atom represents.  All kinds act as anchors for boundary decisions and
   are preserved exactly-once; the kind records provenance for diagnostics/tuning. */
  enum class RoiAtomKind
  {
    ModeledGamma,       // a requested/NORM/interferer source gamma line
    FoundPeakEvidence,  // a data-confirmed found+matched auto-search seed / user peak
    FloatingFeature     // an escape/511/floating-peak feature (no source)
  };

  /** The pipeline stage that first admitted an atom (diagnostics/tuning only). */
  enum class RoiAtomAdmission
  {
    InitialCluster, FallbackEstimate, NoPeakEstimate, FoundPeakSeed,
    UserPeak, RefinementCluster, R2Rescue, EscapeOr511
  };

  /** Stable identity + payload for one admitted modeled line.  IDs are unique per process
   (see next_roi_atom_id) and compared only within a single fit. */
  struct RoiAtom
  {
    uint64_t id = 0;
    double energy = 0.0;                 // keV
    double area = 0.0;                    // expected counts (0 => unknown)
    RoiAtomKind kind = RoiAtomKind::ModeledGamma;
    RelActCalcAuto::SrcVariant source{};  // monostate for evidence/floating atoms
    size_t rel_eff_curve_index = 0;
    RoiAtomAdmission admission = RoiAtomAdmission::InitialCluster;
  };

  /** Mint the next unique atom id (thread-safe; safe under GA parallelism). */
  uint64_t next_roi_atom_id();

  /** One materialized automatic component: channel-aligned bounds plus its exactly-once atoms. */
  struct AutomaticRoiComponent
  {
    double lower = 0.0;                   // == gamma_channel_lower(first_channel)
    double upper = 0.0;                   // == gamma_channel_upper(last_channel)
    size_t first_channel = 0;
    size_t last_channel = 0;
    std::vector<RoiAtom> atoms;           // sorted by energy; exactly-once ownership
    size_t joined_groups = 1;
    bool protected_geometry = false;
    PeakContinuum::OffsetType continuum_type = PeakContinuum::OffsetType::Linear;  // pass-through
    RelActCalcAuto::RoiRange::RangeLimitsType range_limits_type
        = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;                        // pass-through
  };

  /** Stage-independent geometric constraints governing materialization of a partition. */
  struct AutomaticRoiPartitionConstraints
  {
    double lowest_energy = 0.0;           // valid spectroscopic extent (widening clamp)
    double highest_energy = 0.0;
    double left_barrier = -std::numeric_limits<double>::infinity();  // may not widen below this
    double min_width_fwhm = 0.0;          // 0 => impose no minimum child width
    double peak_core_num_fwhm = 1.0;      // atom core half-width, in FWHM
  };

  enum class AutomaticRoiPartitionOutcome { KeptSeparate, Merged, Infeasible };

  struct AutomaticRoiPartitionResult
  {
    AutomaticRoiPartitionOutcome outcome = AutomaticRoiPartitionOutcome::Infeasible;
    std::vector<AutomaticRoiComponent> components;  // 2 (KeptSeparate) | 1 (Merged) | 0 (Infeasible)
    /** Atoms with no legal owner (protected-boundary straddle or spectrum edge only); each is
     accompanied by infeasible_reason and a diagnostic.  Empty in the normal case. */
    std::vector<RoiAtom> orphaned_atoms;
    std::string infeasible_reason;
    AutomaticRoiPolicyResult policy;      // underlying decision + diagnostic (unchanged oracle)
  };

  struct AutomaticRoiTransactionCheck
  {
    bool valid = false;
    std::string failure_reason;
  };

  /** Partition (or merge) one adjacent automatic pair into fully-materialized, channel-aligned
   components with spatially-assigned atoms, or report an explicit infeasible/merge fallback.
   Uses evaluate_automatic_roi_boundary as the decision oracle (unchanged) and owns all geometry
   materialization: core-safe channel boundary search, spatial atom assignment, min-width
   widening, exclusion-band carves, and protected-edge pinning.  Never drops an atom. */
  AutomaticRoiPartitionResult partition_automatic_roi_pair(
    const AutomaticRoiComponent &left,
    const AutomaticRoiComponent &right,
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const GlobalContinuumEstimate *global_continuum,
    const std::function<double(double)> &fwhm_at_energy,
    const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks,
    const AutomaticRoiPolicySettings &settings,
    const AutomaticRoiPartitionConstraints &constraints );

  struct AutomaticRoiReconcileResult
  {
    std::vector<AutomaticRoiComponent> components;   // channel-disjoint, sorted, atom-complete
    std::vector<RoiAtom> orphaned_atoms;             // aggregated; each carries a diagnostic
    bool valid = false;                              // whole-stage transaction validated
    std::string failure_reason;
  };

  /** Result of the measured-data whole-component partition search.  `components` is always a
   transactionally valid replacement when `valid`; `changed` says that the scored optimum has more
   than one child.  A declined or infeasible challenger returns the explicit merged incumbent. */
  struct AutomaticRoiComponentPartitionResult
  {
    std::vector<AutomaticRoiComponent> components;
    bool valid = false;
    bool changed = false;
    std::string failure_reason;
    AutomaticRoiDecisionDiagnostic diagnostic;
  };

  /** Jointly score every core-safe channel boundary producing two children from one over-wide
   transitive component.  Segment scores fit FWHM-distinct peaks plus a production continuum to
   the measured foreground; global AICc is evaluated on the identical union channels and includes
   the existing soft-width pressure.  The atom ledger is preserved exactly once or the incumbent
   merge is retained. */
  AutomaticRoiComponentPartitionResult partition_overwide_automatic_component(
    const std::vector<AutomaticRoiComponent> &component,
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const std::function<double(double)> &fwhm_at_energy,
    const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks,
    const AutomaticRoiPolicySettings &settings,
    const AutomaticRoiPartitionConstraints &constraints );

  /** The single policy-mode reconciliation driver: a left-fold over energy-sorted (possibly
   overlapping) components, folding each adjacent pair through partition_automatic_roi_pair and
   re-examining an enlarged component after a merge.  Validates the whole-stage transaction and
   works all-or-nothing on a copy, so a validation failure leaves `valid == false` and the caller
   retains its incumbent geometry. */
  AutomaticRoiReconcileResult reconcile_automatic_components(
    std::vector<AutomaticRoiComponent> components,
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const GlobalContinuumEstimate *global_continuum,
    const std::function<double(double)> &fwhm_at_energy,
    const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks,
    const AutomaticRoiPolicySettings &settings,
    const AutomaticRoiPartitionConstraints &constraints,
    std::vector<AutomaticRoiDecisionDiagnostic> *diagnostics );

  /** Verify a proposed replacement transaction preserves every invariant: atom-ID multiset
   (before == after together with reported orphans, each exactly once), sorted channel-disjoint
   components, atom energy + clamped-core containment, bit-identical protected bounds/metadata,
   and orphan reasons restricted to protected-straddle / spectrum-edge.  Cheap; always run. */
  AutomaticRoiTransactionCheck validate_automatic_roi_transaction(
    const std::vector<AutomaticRoiComponent> &before,
    const std::vector<AutomaticRoiComponent> &after,
    const std::vector<RoiAtom> &reported_orphans,
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const std::function<double(double)> &fwhm_at_energy,
    const double peak_core_num_fwhm );


  /** One requested source's in-range photon lines, pre-expanded by the caller so that
   find_strong_unmodeled_interferers() stays free of SandiaDecay/NuclideMixture dependencies and is
   unit-testable on synthetic input.  `energies` are already filtered to the valid energy range;
   `yields` are the parallel per-unit-activity intensities (used for single-line / doublet guards). */
  struct RequestedSourceGammas
  {
    RelActCalcAuto::SrcVariant source;
    std::vector<double> energies;
    std::vector<double> yields;
  };

  /** A detected strong line that interferes with a requested-source gamma but is not in the model. */
  struct InterfererCandidate
  {
    double energy = 0.0;                            // interfering line energy (keV)
    const SandiaDecay::Nuclide *nuclide = nullptr;  // co-fit nuclide; nullptr => add a floating peak
    double detection_z = 0.0;                       // area/uncert of the confirming auto-search peak
    bool from_background_search = false;            // false => foreground NORM-table path
  };

  /** One predicted line for the single-pass ROI planner test entry point. */
  struct PlannerLine
  {
    double energy = 0.0;           // keV
    double expected_counts = 0.0;  // predicted peak area, counts
  };

  /** One ROI produced by the planner test entry point. */
  struct PlannedRoiSummary
  {
    double lower = 0.0;
    double upper = 0.0;
    PeakContinuum::OffsetType continuum_type = PeakContinuum::OffsetType::Linear;
    std::vector<double> line_energies;   // the predicted lines assigned to this ROI
  };

  /** Physics envelope test for a claimed peak.

   Could `source_lines` (photon rates at unit activity, e.g. from a NuclideMixture) own a peak of
   `claimed_counts` at `claimed_energy`?  The source's strongest lines (by yield x `intrinsic_eff`)
   constrain its activity from the data: the data's upper limit at every line at least a tenth as
   strong as the claim (net counts within +-1 FWHM over the lower sideband, plus three sigma) bounds
   it from above, an observed peak (from `observed_peaks`, energy/area pairs) at one of the two
   principal lines bounds it from below, and the claimed peak's own size bounds it from above.  Iron and lead shields are scanned (lead's K-edge makes 90-120 keV lines attenuate
   more than 60-88 keV ones, iron has no edge in range); the claim passes if either material can
   explain it.  Principal lines are the two strongest of `source_gammas` (gamma rays only - x-rays
   are shared between nuclides, fluoresce from shields, and are attenuated hardest, so they never
   vouch for a thin shield); `source_lines` holds every photon (gammas and x-rays) for the upper
   bounds and the claim's own yield.  `claimed_fwhm` (0 = the model FWHM) widens the claim's own
   window for blended x-ray peaks.  Every bound depends on the unknown shielding,
   so a lead areal density x is scanned from 0 to `max_shield_g_cm2`: the claim is feasible at x when
   its own activity lower bound (`drf_slack` x claimed counts) does not exceed any line's upper bound,
   and the strong lines are mutually feasible there.  `worst_ratio` is the smallest
   claim-bound / tightest-limit ratio over the mutually-feasible x; above 1 no physical rel-eff curve
   lets the source own the peak (the Cs137 662 keV peak "matched" to Am241's 3.6e-6-branching 662.4
   keV line would need a 335 keV Am241 peak far above what the data allows, and the visible 59.5 keV
   peak rules out the thick shield that could hide it).  A claim within 1.5 sigma of the sum of two
   strong lines (a coincidence-sum peak) and a source whose strong lines cannot agree at any x are
   not judged (`judged` false, `worst_ratio` 0).  Lines within 1.5 sigma of the claimed energy count
   as the claimed peak itself. */
  struct SiblingAbsenceResult
  {
    bool on_source_line = false; // the claimed energy sits on at least one source line (own yield > 0)
    bool judged = false;
    double worst_ratio = 0.0;    // claim lower bound / tightest sibling limit at the best shielding
    double sibling_energy = 0.0; // the binding sibling line
    double required = 0.0;       // counts that sibling would need for the claim to hold
    double limit = 0.0;          // counts the data allows there
    double best_shield_g_cm2 = 0.0;
    float best_shield_z = 0.0f;   // shield material (atomic number) of the best scan point
    double own_yield = 0.0;       // diagnostics: the claim's own line yield ...
    double own_eff = 0.0;         // ... and the efficiency the check used there
    double sibling_eff = 0.0;     // efficiency at the binding sibling
  };

  /** Worst automated-search-peak displacement between two energy calibrations, in units of the
   FWHM: for each peak, where the channel holding it under `from_cal` would be labelled by
   `to_cal`.  Widths come from `fwhm_at` when supplied (the fitted resolution model), else the
   peak's own width - a noise spike's width is not a sensible yardstick for a calibration.  This is what bounds how far a fitted energy calibration may drift from the one the
   spectrum arrived with (see `sm_energy_cal_max_drift_fwhm`).  Peaks with no width are skipped;
   returns 0 when either calibration is unusable or no peak has a width.

   \param worst_energy If non-null, set to the energy of the peak that drifted the most.
   */
  double max_search_peak_drift_fwhm(
    const std::vector<std::shared_ptr<const PeakDef>> &peaks,
    const std::shared_ptr<const SpecUtils::EnergyCalibration> &from_cal,
    const std::shared_ptr<const SpecUtils::EnergyCalibration> &to_cal,
    const std::function<double(double)> &fwhm_at = nullptr,
    double *worst_energy = nullptr );

  /** Generic resolution curve (FWHM in keV) of a detector class, used as the shape prior of
   fit_fwhm_function_robust: NaI for Low/LowOrMedRes/Unknown, LaBr3 for LaBr/MedRes, CZT, HPGe. */
  double class_shape_fwhm( const PeakFitUtils::CoarseResolutionType det_type, const double energy );

  /** Fits the resolution function `form` to the automated-search peaks, anchored to a shape prior
   so the result is a physically shaped curve over the whole analysis range even from two or three
   clean peaks, and never a constant.

   The prior is `shape_prior` (a DRF's curve) when supplied and finite over the range, else the
   detector class's generic curve.  It is scaled by a significance-weighted median of the peaks'
   width ratios and kept inside the class width rails; peaks that disagree with the scaled prior
   (multiplets, backscatter and Compton-edge bumps, starved narrow peaks) do not vote; the
   sqrt-polynomial is then fit by correctly weighted least squares to the surviving peaks plus
   lightly weighted samples of the scaled prior outside their span, so the curve follows the data
   where there is data and the prior where there is none.  `lower_energy`/`upper_energy` are the
   range the curve is valid over: [max(analysis_floor_kev, spectrum min), spectrum max].
   `note` receives a one-line trace ("fwhm model: ...").  Throws when no Gaussian peak is usable.

   On a non-HPGe detector a peak reaching below the spectrum's first live channel does not vote,
   and when no remaining peak is significant enough to measure a width the prior's own scale is
   used; `scale_from_prior`, when given, is set to whether that happened.
   */
  void fit_fwhm_function_robust( const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
                                 const std::shared_ptr<const SpecUtils::Measurement> &foreground,
                                 const PeakFitUtils::CoarseResolutionType det_type,
                                 const double analysis_floor_kev,
                                 const std::function<double(double)> &shape_prior,
                                 const DetectorPeakResponse::ResolutionFnctForm form,
                                 std::vector<float> &coefficients,
                                 std::vector<float> &uncerts,
                                 double &lower_energy,
                                 double &upper_energy,
                                 std::string &note,
                                 bool *scale_from_prior = nullptr );

  SiblingAbsenceResult sibling_absence_check(
    const std::vector<SandiaDecay::EnergyRatePair> &source_lines,
    const std::vector<SandiaDecay::EnergyRatePair> &source_gammas,
    const double claimed_energy,
    const double claimed_counts,
    const double claimed_fwhm,
    const std::vector<std::pair<double,double>> &observed_peaks,
    const std::function<double(double)> &fwhm_at,
    const std::function<double(double)> &intrinsic_eff,
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const double lowest_energy,
    const double highest_energy,
    const double drf_slack,
    const double max_shield_g_cm2,
    const double max_eff_ratio,
    const bool robust_limits = false );

  /** Runs the single-pass ROI planner (GammaClusteringSettings::use_roi_plan) on a bare list of
   predicted lines - the seam unit tests use to exercise grouping, admission, sharing, extent and
   continuum decisions on synthetic spectra.  `fwhm_at` must be valid over [lowest, highest]. */
  std::vector<PlannedRoiSummary> plan_rois_for_lines(
    const std::vector<PlannerLine> &lines,
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const std::function<double(double)> &fwhm_at,
    const double lowest_energy,
    const double highest_energy,
    const GammaClusteringSettings &settings,
    const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks );

  /** R2 keep-gate classification seam for focused statistical tests.  Returns true only for a
   counts-floor-passing cluster in [0.7*keep_z, keep_z], i.e. one that the normal strict keep gate
   rejects but the bounded rescue pass may inspect. */

#if( PERFORM_DEVELOPER_CHECKS )
#endif

  /** Find strong foreground NORM lines NOT in the current model that sit within
   ~`sm_interferer_near_num_fwhm` FWHM of a requested-source gamma and are data-confirmed, so they
   can be considered for auto co-fitting (R6).

   A candidate is a strong-NORM-table line whose parent nuclide is not already modeled and is not
   itself a requested source, that is not explained by the source's own chain, and that is confirmed
   by a foreground auto-search peak within `sm_interferer_confirm_num_fwhm` FWHM at area/uncert >=
   `sm_interferer_min_detect_z`.  The currently active path emits attributable nuclide candidates.
   Ambient Cs137/Co60 scanning, a dedicated background search, and unattributable floating peaks
   remain deliberately disabled until the multi-source nuisance-model behavior is validated.

   If `warnings` is non-null, a human-readable note is appended for each interferer that was detected
   but deliberately NOT co-fit (e.g. an unresolvable single-line-source vs single-line-interferer
   doublet), so the caller can surface it.
   If `guard_energies` is non-null, it receives every data-confirmed interfering-line energy,
   including confirmed doublets that were deliberately not returned as candidates.  This lets the
   bounded rescue pass avoid those ranges without parsing warning text.

   `background`, `drf`, `peak_fit_prefs`, and `global_continuum` are reserved for the disabled
   background/residual-confirmation path and are not currently dereferenced. */
  std::vector<InterfererCandidate> find_strong_unmodeled_interferers(
    const std::vector<RequestedSourceGammas> &source_gammas,
    const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
    const std::function<double(double)> &fwhm_at_energy,
    const bool fit_norm_peaks,
    const double min_valid_energy,
    const double max_valid_energy,
    const std::shared_ptr<const SpecUtils::Measurement> &background,
    const std::shared_ptr<const DetectorPeakResponse> &drf,
    const std::shared_ptr<const PeakFitDetPrefs> &peak_fit_prefs,
    std::vector<std::string> *warnings = nullptr,
    const GlobalContinuumEstimate *global_continuum = nullptr,
    std::vector<double> *guard_energies = nullptr );


  /** Finds the `RelActCalcAuto::FloatingPeakResult` belonging to a bystander user-peak that was
   enrolled into the fit as a `FloatingPeak` at `enrolled_energy`, and marks it consumed.

   `RelActCalcAuto::FloatingPeakResult::energy` is a copy of the input `FloatingPeak::energy`, and
   bystanders are enrolled at exactly their peak mean, so this requires energy identity (to within
   `sm_enrolled_float_match_tol`) rather than a nearest-match over some window.  That matters
   because `m_floating_peaks` also holds floats this code injects for its own reasons - the 511 keV
   annihilation peak, escape peaks, and auto-detected interferers.  A nearest-match can bind an
   unmatched bystander to one of those (e.g. a source-less user peak at ~510.7 keV to the 511 keV
   float), after which the de-duplication pass erases the real 511 keV peak from the results.

   Each result is returned at most once - `consumed_results` is both read and updated - so two
   bystanders enrolled at the same energy bind to different results.

   Returns nullptr when the bystander's own float did not survive to the solve (e.g. it was dropped
   by `remove_floating_peaks_without_roi`); callers then retain the user's original peak. */
  const RelActCalcAuto::FloatingPeakResult *find_enrolled_float_result(
    const std::vector<RelActCalcAuto::FloatingPeakResult> &floating_results,
    std::set<const RelActCalcAuto::FloatingPeakResult *> &consumed_results,
    const double enrolled_energy );


  /** Builds the updated bystander peak from the fit's `FloatingPeakResult`, or indicates that the
   user's original peak should be retained instead.

   Shared by the `ExistingPeaksAsFreePeak` and default-mode reconciliation blocks of
   `fit_peaks_for_nuclide_relactauto` so their acceptance rules cannot drift apart.

   Returns `std::nullopt` when the solver did not actually determine this peak - a non-positive
   amplitude, or an amplitude whose uncertainty is missing, non-finite, or as large as the
   amplitude itself.  Retaining the original is the conservative choice there: swapping in an
   undetermined amplitude trades the user's measured peak for a meaningless one, and pairing the
   fit's new amplitude with the original's (small) uncertainty would report a confident peak that
   nothing ever measured.

   `mode_label` prefixes the `PERFORM_DEVELOPER_CHECKS` diagnostics, to identify the calling block. */
  std::optional<PeakDef> update_bystander_from_float_result(
    const PeakDef &orig_peak,
    const double orig_energy,
    const RelActCalcAuto::FloatingPeakResult &fpr,
    const std::shared_ptr<PeakContinuum> &roi_continuum,
    const char * const mode_label );
}//namespace detail


// Settings for the gamma clustering algorithm - different values may be used
// for the initial RelActManual stage vs subsequent RelActAuto refinement stages
struct GammaClusteringSettings
{
  // False only for an R6-enabled fit, whose source incumbent and nuisance transaction must retain
  // their complete legacy geometry behavior.  Not serialized or tuned.
  bool use_automatic_roi_policy = true;
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

  // Optional shared SNIP-based global continuum for gating B(E) estimates (R1 step 2).  Non-owning;
  // NULL => every consumer falls back to its local two-sideband estimate (pre-R1-step2 behaviour).
  // Lifetime is a synchronous stack frame in fit_peaks_for_nuclides.
  const detail::GlobalContinuumEstimate *global_continuum = nullptr;

  /** Whether `global_continuum` may supply the BACKGROUND the admission gate and the step
   statistics measure against, as opposed to the local sideband estimate.

   The planner needs a SNIP continuum for its curvature statistic whatever else is going on, but a
   SNIP estimate rides high inside a dense multiplet - it has no peak-free channels to dig down to -
   so using it as the gate's background there costs real lines (a Detective-X Pu239 fit lost seven
   of its 330-450 keV lines and another timed out).  Supplying the continuum for curvature while
   leaving the gate on its local estimate keeps both behaviours honest.
   */
  bool global_continuum_gates_admission = true;

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


  // Step continuum decision thresholds
  // Minimum peak detection significance z = S_est / sqrt(S_est + B_est), with B_est from the
  // sideband continuum estimate, to consider a step continuum (a step only matters when the peak
  // towers over the continuum, which is inherently a significance statement - an absolute-count
  // gate would not transfer across live-times).
  double step_cont_min_peak_significance = 40.0;
  // Chi2 margin by which the step-continuum trial fit must beat the polynomial fit for the ROI to
  // get a step continuum (the trial pairs equal-parameter-count candidates - Linear vs FlatStep,
  // Quadratic vs LinearStep - so the AICc penalty terms cancel and the decision reduces to a
  // chi2 comparison with this tunable bias against the step).  Replaces the former left-vs-right
  // probe-window nsigma test, which self-vetoed on tight ROIs and read neighbor peaks as steps.
  double step_trial_chi2_margin = 4.0;

  // Provenance stage stamped on atoms minted at the keep-gate (diagnostics only; not serialized or
  // tuned).  Refinement re-clustering sets this to RefinementCluster.
  detail::RoiAtomAdmission cluster_admission_stage = detail::RoiAtomAdmission::InitialCluster;

  // Single-pass ROI planner (see plan_rois_impl in FitPeaksForNuclides.cpp).  When enabled it
  // replaces the greedy cluster/merge/split cascade of cluster_gammas_to_rois with one decision
  // pass: dominant-anchored line groups, data confirmation by automated-search peaks, one
  // detectability gate, share/separate by separation bands with a clean-gap test in between,
  // sideband-driven extents clipped at neighbours, and per-ROI continuum order/step selection.
  bool use_roi_plan = true;
  // The reference fits share below ~3 FWHM separation and separate above ~4; in between the
  // clean-gap test decides (a hard cut at 3.4 scored worse on the reference corpus, 2026-09-05).
  double share_always_fwhm = 3.0;         // adjacent groups closer than this always share a ROI
  double separate_always_fwhm = 4.0;      // ... farther than this never share
  // When > 0, replaces the three rules above: adjacent groups separate when, cut at the valley of
  // their combined predicted signal, each side's predicted counts spilling across the cut are at
  // most this fraction of the other side's.  0 keeps the distance rules.
  double share_max_leak_fraction = 0.0;
  // With share_max_leak_fraction: separate only where the data come back down to the continuum - the
  // gross counts within a quarter FWHM of the cut may exceed the SNIP continuum there by at most this
  // fraction of it (or two sigma).  A cut on a peak's flank leaves the next ROI, which carries only
  // its own lines, to absorb that peak's tail: its continuum dives at the edge under an edge peak
  // (Tl201 135/167 keV, Br82 554/619 keV).  0 disables.
  double share_valley_max_excess_fraction = 0.0;
  // An ROI that channel alignment (with its one-channel gap to the previous ROI) leaves under three
  // channels takes the gap and any free room above it before it is dropped: dropping loses its lines
  // outright - Tl201's z=154 167 keV line on a 12.5 keV/channel NaI, squeezed to 164-189 keV.
  bool widen_roi_without_room = false;
  // Where two ROIs touch, the one-channel gap channel alignment leaves between them is the channel
  // holding their boundary, rather than the channel above the lower ROI's rounded-up edge.  The old
  // placement took one to two channels from the upper ROI alone, which on a 3 keV/channel NaI is up
  // to 0.8 FWHM of low-side sideband (a 43.5 keV Pu238 ROI started 0.3 FWHM below its line).
  bool roi_gap_at_boundary = false;
  // With share_max_leak_fraction, on a coarse binning: separate only where each side keeps at
  // least share_min_side_channels between its nearest visible line and the cut, and the component
  // being closed off spans at least share_min_roi_channels from its own lower boundary to the cut.
  // The leak cut sits about a FWHM from each line, which there is two or three channels: Pu238's
  // z=84 43.5 keV line, cut at the valleys 1.9 FWHM from its neighbours on a 2.9 keV/channel NaI,
  // got a three-channel ROI that fit worse than no peak.  The span test, not a larger per-side one,
  // decides which neighbour a cramped group joins; a per-side test alone glued Pu238's 15-132 keV
  // into one ROI through a weak escape group sitting against both of its cuts.  0 disables each.
  double share_min_side_channels = 0.0;
  double share_min_roi_channels = 0.0;
  // The extent of an ROI is built around the lines one could see (visible_line_min_z); with this set
  // every joined group's dominant line is part of that core, so a group admitted on the data's
  // evidence but predicted invisible cannot end up outside the ROI it was joined into (Pd103_Phantom's
  // z=10 39.8 keV line, joined to the Rh x-rays for room, got a 16-31 keV ROI).
  bool roi_core_covers_every_group = false;
  // With share_max_leak_fraction: separate only where the predicted signal has a valley - its value
  // at the cut at most this fraction of the smaller dominant peak's height.  0 disables.
  double share_max_valley_depth = 0.0;
  // Admit a group the search did not confirm only if, at one of its lines (with at least a tenth of the
  // group's largest predicted counts), a peak of the model width over a QUADRATIC continuum reaches this
  // significance in the data (detail::fixed_shape_peak_z).  A linear continuum under the smooth
  // low-energy scatter hump leaves room for Gaussians wherever lines are predicted there (Eu152's
  // 122 keV behind Pb, La140 at 128 keV); the quadratic follows the hump.  On 96 R500 spectra, z < 2
  // rejected half the ROIs whose line the truth has invisible, and 5 % of those it has visible (mostly
  // marginal lines).  0 disables.
  double admission_min_data_z = 0.0;

  /** With admission_min_data_z: the group's DOMINANT line (most predicted counts), when it can be
   judged, must itself show a peak at admission_min_data_z - it or a line within half a FWHM of it,
   which makes one peak with it; otherwise one minor line of the group passing the test admitted it.
   R500 W187_Sh's W K x-ray group was predicted at z 53 (the rel-eff misses the shield's x-ray
   attenuation) and showed z 0.0 at 72 keV, yet a minor line passed and the fit delivered a 44-85 keV
   ROI with a phantom on the scatter hump. */
  bool admission_dominant_line_required = false;
  // Admit a group the prediction gate rejected when one of its lines is a strong line of its source
  // (decay yield >= 5 % of the source's strongest in the analysis range) and the data show a peak of
  // the model width there at this significance (detail::fixed_shape_peak_z).  The prediction is what
  // fails for these: a rel-eff extrapolated to 1.8 MeV put Br76's clearly visible 1854 keV line at
  // z=1.1.  0 disables.
  double admission_data_evident_min_z = 0.0;
  // Predict each line's NaI/CsI iodine K x-ray escape peaks (PeakFitUtils::nai_iodine_escape_fractions)
  // as lines of their own source, so ROIs are planned around them and a search peak there counts as
  // accounted for; the solve models them with RelActCalcAuto::Options::iodine_escape_peaks.
  bool iodine_escape_peaks = false;
  double found_peak_match_num_fwhm = 0.5; // auto-search peak within this of a line confirms the group
  // Quadratic continuum: only for ROIs wide enough to show curvature AND carrying enough continuum
  // counts to measure it.  The reference fits are unambiguous on the second point - among their ROIs
  // at least 6 FWHM wide, the fraction with a quadratic/cubic continuum runs 1 % below 1000
  // continuum counts (hence the default), 15 % from 1000-5000, 42 % from 5000-20000 - and their
  // quadratic ROIs are also
  // the wide, crowded ones (median 8.4 FWHM and 4 peaks, against 6.3 FWHM and 1 peak for linear).
  // Past both gates the choice is made by fitting the whole ROI both ways (peak amplitudes free)
  // and comparing by AICc; <= 0 on either gate means "never quadratic".
  double quad_min_width_fwhm = 8.0;
  double quad_min_continuum_counts = 1000.0;
  // ... and a second, independent route to a quadratic: the SNIP continuum's own curvature across
  // the ROI.  A peak sitting near the crest of the broad backscatter/x-ray hump has a continuum
  // that is visibly convex under it, and no straight line fits that however the peaks are modelled.
  // Measured over thirds of the ROI: the middle third's continuum counts minus the mean of the two
  // outer thirds (zero for any straight line, whatever its slope), against the Poisson noise of the
  // continuum counts themselves; |z| past this is curvature.  <= 0 disables.
  double quad_min_curvature_z = 3.0;
  double found_peak_min_predicted_fraction = 0.25; // confirmation needs predicted counts >= this x found area
  // A step continuum is chosen on the SIZE of the step, not only its significance: the Compton step
  // under a peak is always physically present, so the question is whether it is big enough to bend
  // the fit.  The step fraction (low-side minus high-side continuum density, over the low-side
  // density) separates the reference fits better than the asymmetry z did - the z gate and the
  // fraction gate disagree with the reference on 96 and 86 of 1120 ROIs respectively (2026-09-07).
  double step_min_asym_z = 1.5;           // low-minus-high sideband z (0 = only require a downward step)
  // ... and the step must be measurably a step: its size against the Poisson noise of the two
  // flank samples (step_min_asym_z), not a fixed fraction of the continuum.  With the flanks
  // sampled correctly (see sm_step_core_num_fwhm) significance beats a fraction cut on the
  // reference set - 75 disagreements of 983 against 78 for the best fraction.
  double step_min_fraction = 0.0;
  bool step_use_chi2_trial = false;       // also require the LLS step-vs-polynomial chi2 trial to win
  double step_low_side_extra_fwhm = 0.5;  // extra low-side ROI room for step continua (the reference fits give it)
  // Physics envelope for found peaks (detail::sibling_absence_check): a found peak the source cannot
  // own without a stronger sibling line the data shows absent never confirms a group; the group under
  // it is rejected as swamped and the peak stays an obstacle.  <= 0 disables the test.
  double sibling_absence_max_ratio = 5.0;      // genuine contaminants score 10-1000; 2-4 is DRF-shape territory
  double sibling_absence_drf_slack = 0.4;      // generic-DRF shape / summing tolerance on the claim's activity bound
  double sibling_absence_shield_g_cm2 = 120.0; // heaviest lead areal density (g/cm2) the shield scan considers
  // A constraint line is only usable evidence when its detection efficiency is within this factor
  // of the claim's.  Past that the verdict is the efficiency SHAPE talking, not the physics: the
  // generic coaxial curve the check falls back on when no DRF is supplied has all but vanishing
  // efficiency at 20-35 keV, so on a planar or low-energy HPGe - where the K x-rays are among the
  // largest peaks in the spectrum - it concluded that an I123 27.4 keV x-ray peak of 93,000 counts
  // would require 1.5e8 counts in the 529 keV line, and rejected every K x-ray complex in the set.
  // <= 0 disables the guard (every line is usable).
  double sibling_absence_max_eff_ratio = 50.0;
  // How much room the data leaves a sibling line, and which sources may explain a peak, judged for
  // a scintillator: a sideband lying on the claimed peak or on another strong line of the source is
  // not continuum; the continuum under the sibling cannot exceed the data's own lowest density in its
  // window; and a peak is unexplained only if NO source with a line in the group can account for it.
  // Off keeps the original (HPGe-tuned) limits.
  bool sibling_absence_robust_limits = false;
  std::shared_ptr<const DetectorPeakResponse> sibling_check_drf;  // efficiency shape for the test; callers set it, null skips it
  // Data-detected admission ABOVE the spectroscopic extent (0 = off, the default; below the extent
  // it is always available, being the only route there): a group the prediction rejects but whose
  // predicted z is at least this is admitted when the spectrum itself shows the peak (net counts
  // over the core window at data_detect_min_data_z or better) and the physics envelope does not
  // rule it out, without the automated peak search having found it.  It was on while the rel-eff
  // curve could run away and the FWHM model could be 50 % wide; with the order cap and the robust
  // width model the ordinary admission path finds those lines, and switching it on now costs 20
  // raw on the corpus and quadruples the wall time (the ROIs it adds are mostly junk).  Kept as an
  // option for spectra with no usable prior information at all.
  double data_detect_min_predicted_z = 0.0;
  double data_detect_min_data_z = 3.0;      // net counts over the core window must reach this z (not the stage's keep z)

  /** Widest window, in FWHM, over which the SNIP estimate may supply the admission gate's
   background; wider than this the local sideband estimate is used instead.  Zero means no limit.

   SNIP finds the continuum by clipping features narrower than its window, so it needs peak-free
   channels to settle on.  Across a wide window of overlapping lines there are none and the estimate
   rides up on the peak shoulders; as the gate's background that inflates B, deflates z and rejects
   real lines - a Detective-X Pu239 fit lost seven of its crowded 330-450 keV lines exactly that
   way.  The local sideband estimate fails in the mirror way (it wants clean sidebands), so this is
   not "SNIP is worse", it is that the two err in opposite directions and the gate should use
   whichever one the window can support.
   */
  double snip_gate_max_window_fwhm = 0.0;

  /** Minimum significance a predicted line needs before it counts as VISIBLE, i.e. before it may
   set how far a ROI reaches and therefore whether two groups share one.  Zero keeps the old
   behaviour, where any line carrying a tenth of its group's brightest line was "visible".

   That fraction is relative to the group, not to the data: beside a strong line, a tenth of it can
   still be far under the continuum, and a chain of such lines walks a ROI across a valley the
   spectrum plainly shows.  On a uranium-ore NaI spectrum the planner put ONE region across
   927-2439 keV, joining the 1120, 1408, 1764 and 2204 keV peaks because unmeasurable lines between
   them left gaps of half a FWHM, where a hand fit uses five separate regions of 2-3 FWHM each.
   */
  double visible_line_min_z = 0.0;

  /** Refutation by the data: a line group predicted at `refute_min_predicted_z` or better whose
   own window shows less than `refute_max_data_fraction` of the predicted counts is not there,
   whatever the rel-eff says, and is rejected.

   The admission gate reads S/sqrt(S+B) from the PREDICTION, and the prediction comes from a
   rel-eff extrapolated from a handful of matched peaks, so it can be absurd: an Am241 fit
   predicted z=1492 at 59.5 keV, about two million counts, where the whole window holds three
   thousand.  A phantom that size anchors a 30 FWHM region, and because its predicted tail is
   subtracted from the sidebands it also zeroes the continuum estimate there.

   The comparison is against the window's GROSS counts, deliberately.  A peak cannot contain more
   counts than the window it sits in, whatever the continuum is doing, so this needs no continuum
   estimate and cannot be fooled by one.  Judging it against a local continuum instead cost 26 real
   strong peaks on one corpus: where the window sits on a steep Compton slope the sideband estimate
   exceeds the gross, the "net" goes hugely negative (-37726 counts under a z=63 peak) and a real
   line is refuted for a fault in the background model rather than any absence of signal.
   */
  double refute_min_predicted_z = 0.0;
  double refute_max_gross_multiple = 2.0;

  /** Fraction of a candidate sideband window's counts that MODELLED lines may explain before the
   window is rejected as measuring a peak rather than the continuum.  0 disables the test, which is
   the HPGe behaviour: there the peaks are narrow enough that a sideband placed against the region
   edge is already clean, and forcing the sample further out measured strictly worse (strong misses
   2 -> 4 and a mechanical failure on the Detective-X set).  On a scintillator the flank of a line
   two or three FWHM away is exactly what the sample lands on. */
  double sideband_max_predicted_fraction = 0.0;

  /** When two neighbouring components are closer than `roi_min_side_fwhm` allows, put a boundary
   between them instead of merging them into one region.  The boundary goes at the lowest point of
   the data between their outermost lines, so the two regions BUTT against each other rather than
   each keeping a clear sideband.  <= 0 keeps the merge.

   Hand fits do exactly this: over six NaI spectra the median gap between adjacent hand-drawn
   regions is 1.1 FWHM, two thirds are under 1.5, and several touch or overlap slightly - so most
   hand-drawn boundaries are ones a merge-on-small-gap rule can never reproduce.  The boundary is
   still kept `roi_touch_min_line_fwhm` away from any modelled line on either side, since a
   boundary sitting on a peak flank is what makes a continuum collapse.
   */
  double roi_touch_split_min_fwhm = 0.0;
  double roi_touch_min_line_fwhm = 0.6;
  // A ROI wider than this (the outermost lines plus a core on each side, in FWHM) is split at its
  // widest internal gap, repeatedly: one polynomial continuum does not hold across a 25 FWHM x-ray
  // region, and the reference fits break such regions into 7-9 FWHM pieces.  <= 0 disables.
  double max_shared_span_fwhm = 12.0;
  // Every ROI edge must sit at least this many FWHM from the outermost line inside it, so the fit
  // has continuum to anchor on.  Without it a split or a neighbour clip could leave a modelled line
  // a tenth of a FWHM from the edge, and the fitted continuum then collapsed toward zero while the
  // peak swallowed it (Ac225 81.5 keV, Pu239 658.9 keV, Br76 1224.5 keV: continuum ~0 at the edge
  // against 500, 7 and 2 counts of data).  A gap too small to give both neighbours this room is not
  // a place to split: the components stay in one ROI instead.
  double roi_min_side_fwhm = 1.0;
  // A continuum is linear (or a step) almost always: that is what the reference fits use throughout,
  // and a real continuum does slope.  A Constant continuum is used only where a slope cannot be
  // measured at all - BOTH sidebands agreeing within cont_constant_max_asym_z sigma AND fewer than
  // cont_constant_max_counts counts of continuum (about four counts per channel at the default,
  // ~18 % of this corpus's ROIs).  Measured on the 166-problem corpus: never 856.6 raw / 1379
  // reference peaks matched, under 100 counts 851.0 / 1380, under 300 counts 847.0 / 1383,
  // unconditional lower still but a third to a half of all ROIs - the gain is in weak peaks a free
  // slope would otherwise absorb, which is why it is capped at the low-statistics end.
  // Either value <= 0 disables it.
  double cont_constant_max_asym_z = 1.0;
  double cont_constant_max_counts = 100.0;
  // Where netting the sidebands of the lines' prediction left one of them empty, judge
  // cont_constant_max_asym_z from the RAW sideband counts: an over-predicted line nets both sides to
  // zero, which reads as "no slope" and forced a Constant continuum - R500 U233_Sh's 2614 keV ROI, whose
  // ~0 constant then sat under 8-10 counts/channel once the side cap trimmed it.  Elsewhere the netted
  // sidebands stand (raw ones keep the line's tails, and tilted sparse LaBr3 ROIs linear - Cu64
  // 1346 keV), and cont_constant_max_counts stays on the netted ROI continuum.
  bool cont_constant_from_raw_sidebands = false;
  // Obstacles (found peaks no source line explains) stay out of a ROI: its edge stops this many of
  // the obstacle's FWHM from the obstacle mean (3.5 sigma: the obstacle's Gaussian is 0.2 % of its
  // height there), or at the data valley when that would leave less than obstacle_min_side_fwhm of
  // sideband beyond the ROI's outermost line.  The extension test alone only refuses blocks
  // touching an obstacle's half-width, so a ROI could end on the flank of a neighbouring peak (a
  // Lu177m 233.9 keV ROI ended 0.3 FWHM below the Pb212 238.6 keV peak; the flank tilted its
  // continuum and the peak area came out 9 sigma high).  <= 0 disables.
  double obstacle_exclusion_fwhm = 1.5;
  double obstacle_min_side_fwhm = 1.0;
  // An obstacle within this many FWHM of one of a ROI's own lines (of at least a tenth of its dominant
  // line's predicted counts) is that line - displaced by the
  // energy calibration, or a peak the source's first estimate under-predicted - and does not clip the
  // ROI; 0 disables.  NGH Ca47_Sh's ROI was cut at 1303 keV "below the 1308.6 keV obstacle", its own
  // z=53 1297 keV line 0.16 FWHM away, and the line was lost; ~10 ROIs across the scintillator sets
  // were clipped at a search peak within 0.2 FWHM of their dominant line.
  double obstacle_own_line_fwhm = 0.0;
  // Below the spectroscopic extent (the peak search's turn-on estimate, which walks past the
  // low-energy fluorescence/backscatter band) the planner trusts no prediction: a line group whose
  // lines all lie below this energy is admitted only when the spectrum itself shows the peak
  // (data-detected admission) - the Ba133 K x-rays at 31/35 keV are real, Np237's hundreds of
  // L x-ray lines in the same band are not worth 470 s of solving.  Callers set it; 0 disables.
  // Such admission is limited to single-source groups of at most sub_extent_max_lines lines (K x-ray
  // clusters and lone gammas; the actinide L x-ray forests chain into 16 keV ROIs otherwise).
  double sub_extent_energy = 0.0;
  int sub_extent_max_lines = 6;
  double low_energy_abs_floor = 0.0;   // see PeakFitForNuclideConfig::low_energy_abs_floor
  bool low_energy_skip_threshold_ramp = false;   // see PeakFitForNuclideConfig
  // Keep an ROI's lower edge out of the detector turn-on: raise it to sub_extent_energy when its own
  // lines sit at least roi_min_side_fwhm above that.  A sideband reaching into the turn-on gave the
  // solve a rising edge no continuum follows, and it fit that edge with the TAILS of L x-rays below
  // every ROI at an extrapolated efficiency of e^16, zeroing every real activity (Tl201_Unsh with
  // its iodine-escape ROI at 17.6-55.7 keV).
  bool roi_floor_at_extent = false;
  // The live-time-scaled background for the data tests (see PeakFitForNuclideConfig::planner_net_data_tests);
  // callers set it, null leaves the tests on the gross spectrum.
  std::shared_ptr<const SpecUtils::Measurement> background;
  double background_scale = 0.0;
};


// Result from RelActAuto peak fitting
struct PeakFitResult
{
  RelActCalcAuto::RelActAutoSolution::Status status;
  std::string error_message;
  std::vector<std::string> warnings;  // warnings that don't prevent success

  // Structural decisions made while constructing automatic ROIs, in solve-calibration order.
  std::vector<AutomaticRoiDecisionDiagnostic> automatic_roi_diagnostics;

  // Human-readable decision trace of the single-pass ROI planner (one entry per planned ROI or
  // rejected line group), in the order the planner ran; empty unless use_roi_plan.
  std::vector<std::string> roi_plan_trace;

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

  /** The ROIs the planner produced, before the solve ran.  Scoring these alone answers "would a
   peak here have been fit at all", which is what most misses come down to, and costs milliseconds
   against ~10 s for the solve - so a planning change can be checked across every corpus in a
   minute instead of hours.  Populated whether or not `stop_after_plan` skips the solve. */
  std::vector<RelActCalcAuto::RoiRange> planned_rois;

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
  double manual_roi_core_num_fwhm = 2.5;  // Always-included ROI half-extent beyond outermost gamma

  // RelActAuto parameters
  /** The solve's resolution form: noise plus a power law with a varying log-log slope, which is
   positive and non-decreasing by construction (see RelActCalcAuto::FwhmForm::NoisePlusCurvedPower). */
  RelActCalcAuto::FwhmForm fwhm_form = RelActCalcAuto::FwhmForm::NoisePlusCurvedPower;

  /** How the solve treats the width function it is handed.

   The default lets the solve refine the resolution while it fits.  On a scintillator that freedom
   is not free: below ~300 keV the continuum families cannot follow the broad convex shoulder, and a
   peak that is allowed to grow WIDE is a cheap way to absorb it - fitted widths run ~30 % over the
   truth below 100 keV while being correct above 300 keV, and a peak fitted much too wide carries
   several times its true area.  `FixedToAllPeaksInSpectrum` instead holds the curve at the estimate
   made from the peaks in the spectrum (which, with a starting curve supplied, is the planner's own
   fit - see detail::fit_fwhm_function_robust), removing that degree of freedom entirely.
   */
  RelActCalcAuto::FwhmEstimationMethod fwhm_estimation_method
      = RelActCalcAuto::FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum;
  double rel_eff_auto_base_rel_eff_uncert = 0.1;  // BR uncertainty for RelActAuto
  // Smallest yield multiple the solve may give a peak range (see
  // RelActCalcAuto::Options::additional_br_min_yield_fraction); 0 lets it switch a line off.
  double rel_eff_auto_br_min_yield_fraction = 0.0;
  /** The solve's widths stay under this multiple of the planner's width model everywhere (see
   RelActCalcAuto::Options::fwhm_max_ratio_to_start; 0 = no limit).  Where a solve's low-energy
   widths ran to 1.6x the model or more - R500 25 of 104 solves, NGH 41 of 103 at 55-100 keV - the
   model matched the truth widths (median ratio 0.93-1.00) and the solve did not (1.55-1.67): x-ray
   blobs and scatter humps pull a free noise floor and flat slopes wide (R500 Yb169_Phantom ran
   ~17.5 keV at 10-60 keV against a true 8.7 keV at 60). */
  double rel_eff_auto_fwhm_max_ratio_to_model = 0.0;
  /** A first solve whose widths below 300 keV ran past this multiple of the planner's width model
   (itself no narrower than the solver's channel floor) is re-solved with the widths held within
   width_balloon_resolve_limit of the model (RelActCalcAuto::Options::fwhm_max_ratio_to_start), and
   the rest of the fit keeps that limit; 0 = never.  Limiting every solve instead
   (rel_eff_auto_fwhm_max_ratio_to_model) reshuffled the healthy ones through the retry and the
   refinement.  The worst balloons are a flat curve on a wide noise floor that models the whole scatter
   plateau with broad peaks over a sunken continuum: R500 Ta182_Sh solved ~55 keV at 55-264 keV (model
   10-17 keV), its 79-246 keV ROI three wide peaks on a step 600 counts under the data.  Measured at
   1.6 / 1.4 and left off (2026-09-26): 32 better / 28 worse of 127 changed spectra - the balloon also
   stands in for unresolved x-ray blends and the scatter hump, which narrowed widths then under-fit. */
  double width_balloon_resolve_ratio = 0.0;
  /** See width_balloon_resolve_ratio. */
  double width_balloon_resolve_limit = 1.4;
  /** The solve's narrowest allowed width in widths of the widest ROI channel (see
   RelActCalcAuto::Options::fwhm_channel_floor_factor; 1.25 is the solver's default). */
  double rel_eff_auto_fwhm_channel_floor_factor = 1.25;

  // ROI clustering thresholds for auto RelEff refinement stage
  double auto_rel_eff_cluster_num_sigma = 2.0;  // Slightly wider clustering with better rel-eff
  double auto_keep_significance_z = 3.0;  // Cluster keep-gate z (see manual_keep_significance_z)
  // The reference fits keep about 3 FWHM per side; a 1.5 FWHM core left the planner's ROIs
  // systematically 0.3-0.45 FWHM narrow on both sides.  1.5 -> 2.0 was worth ~60 raw on the corpus;
  // 2.0 / 2.5 / 3.0 are within its noise (1115 / 1119, and later 909.6 / 896.6 / 908.9), so 2.5 is
  // a choice inside a tie, not a measured optimum.
  double auto_roi_core_num_fwhm = 2.5;    // Always-included ROI half-extent beyond outermost gamma
  double auto_rel_eff_sol_max_fwhm = 12.0;  // Tighter constraint as solution improves
  double auto_rel_eff_sol_min_fwhm_roi = 1.25;
  // Permit the automatic late-geometry reconciler to score a core-safe partition of an
  // over-wide component.  Disabled by default while detector-specific evidence is gathered.
  bool auto_roi_partition_overwide = false;
  // Separately admit a partition proposed from a completed solver ROI.  This path has more
  // complete peak evidence than the early predicted-gamma reconciler, so it must be tuned and
  // validated independently rather than being an implicit consequence of the early policy.
  bool auto_roi_final_fitted_partition = false;
  // Minimum modeled-core separation required before a final automatic ROI can be partitioned.
  // Zero retains the established AICc-only behavior.
  double auto_roi_partition_min_gap_fwhm = 0.0;
  // When enabled, a global-SNIP-aware clean continuum window may admit a partition below the
  // core-gap rail.  Disabled by default so existing configurations retain their geometry.
  bool auto_roi_partition_allow_clean_gap_override = false;
  // Opt-in residual-SNIP valley admission for late sparse-ROI partitions. Zero disables it.
  // This is deliberately distinct from merge_tail_z: it never changes ordinary ROI merging.
  double auto_roi_partition_residual_valley_max_excess_z = 0.0;
  // Limit expensive post-solve local challengers per refinement.  The normal solver can revisit
  // geometry on later refinement passes; this prevents dense spectra from multiplying solves.
  size_t auto_roi_final_partition_max_proposals = 1;
  // Late partitions are intended for sparse, visibly separated fitted structures.  A dense
  // collapsed multiplet is both physically inappropriate to split and expensive to score.
  // Zero disables this atom-count admission rail.
  size_t auto_roi_final_partition_max_atoms = 12;
  // Optional lower width threshold for the late fitted-ROI challenger.  Zero retains the normal
  // auto_rel_eff_sol_max_fwhm admission threshold; a larger value targets only genuinely broad
  // fitted ROIs without spending the bounded proposal budget on ordinary borderline ranges.
  double auto_roi_final_partition_min_width_fwhm = 0.0;
  // Maximum children a clean-valley whole-component challenger may propose in one transaction.
  // Two retains the original binary behavior.
  size_t auto_roi_partition_max_children = 2;
  double auto_roi_partition_force_gap_fwhm = 0.0;

  // Adaptive ROI sideband extension (shared between manual and auto stages): block-consistency z
  // and extension cap in FWHM beyond the outermost gamma.  See detail::extend_roi_by_sidebands.
  double roi_extend_z = 2.0;
  double roi_max_num_fwhm = 3.5;   // extension cap per side; 3.5 scored best on the reference fits (2.5/3/3.5/4: 1328/1186/1183/1258)

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

  /** Leave out of the manual rel-eff stage the search peaks narrower than this fraction of the
   planner's width model at their energy; 0 disables.  A search peak that narrow is a noise spike or a
   weak line its narrow fit under-measures: of the scintillator search peaks under 0.5x the model, 303
   sit on no truth line and 47 do (those at z 3-5), against ~1500 real ones at 0.8-1.4x.  NGH Ca47_Sh's
   490 keV "peak" (FWHM 3.7 keV, 72 counts; the line holds 451) bent the rel-eff curve so the 1297 keV
   line (z 53) was predicted at 1/15 of the data, excluded as a contaminant, and its ROI clipped at
   its own search peak. */
  double manual_min_search_fwhm_ratio = 0.0;


  // RelActAuto relative efficiency model parameters
  // Absolute low-energy floor (keV) for planning source lines, below which no ROI is ever planned.
  // A DRF with a lower valid energy overrides it downward; 0 means "the built-in default" (20 keV
  // for HPGe, 25 otherwise).  The built-in defaults suit a coaxial detector, whose spectrum has
  // nothing usable that low - but a planar or low-energy HPGe puts real, strong lines at 10-18 keV
  // (a 50 % planar misses 99 of 100 strong truth peaks below 20 keV with the coaxial floor).  What
  // actually protects the fit is the data-alive walk in planning_low_energy_bound plus the
  // sub-extent admission rules, not this constant, so it can be set low for such detectors.
  double low_energy_abs_floor = 0.0;
  // ... and the data-alive bound is moved past the detector's threshold ramp (see detail
  // threshold_ramp_top): the few channels the discriminator cuts through climb from almost nothing to
  // the spectrum's level, an ROI or data-test window reaching into them fits the climb as a peak (a
  // CZT's 40-78 keV ROI, U233_Sh, and the planner's own data test read the ramp as a z=9 peak).
  bool low_energy_skip_threshold_ramp = false;
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

  /** Significance a requested-source peak must have in the incumbent before losing it makes the
   refinement pass unacceptable.  <= 0 uses the built-in max(8, 2*roi_significance_z).

   This guard is dangerous on a low-resolution spectrum, and the danger is self-reinforcing: where a
   wide ROI carries a broad continuum hump on a comb of Gaussians, every comb member has an enormous
   fitted significance, so the refinement that would dissolve the comb and split the ROI is seen to
   "remove strong peaks" and is rejected - the bad geometry manufactures the evidence that protects
   it.  Measured firing on ~25 % of problems across five NaI detectors.  Set very high to disable.
   */
  double refinement_anchor_guard_min_z = 0.0;

  /** Only requested-source peaks the automated search found in the data (a significant,
   photopeak-shaped search peak within half a FWHM) may anchor the refinement guard.  On a
   scintillator the solve's tied lines all carry large fitted significances, so without this a comb
   of Gaussians on the scatter hump vetoed the refinement that would dissolve it; with the guard off
   entirely, a solve that dropped U235's z=55 185.7 keV line was accepted instead. */
  bool refinement_anchor_guard_confirmed_only = false;

  /** With refinement_anchor_guard_confirmed_only: one anchor per search peak and source - the lines a
   source has under one found peak - held when the challenger keeps any of them significant, rather
   than one anchor per gamma line.  Blended members of a decay chain trade fitted significance among
   themselves whenever the ROIs change, so counted line by line a U233 spectrum had 143 "anchors",
   and the keep-80 % rule refused the refinement that split a 6.6 FWHM ROI into the two a hand fit
   uses.  An anchor is held while the source keeps half its significance under the peak (at most
   roi_significance_z is asked), so a marginal anchor drifting from z=3.6 to 3.4 is not a loss; a
   source that loses the peak outright still fails the guard. */
  bool refinement_anchor_guard_per_found_peak = false;

  /** When the first solve zeroes the activity of a requested nuclide, re-solve without the ROIs in
   the low-energy (x-ray / scatter-hump) region and, if that keeps every nuclide's activity, start
   the refinement from it.  See the note where it is applied in fit_peaks_for_nuclide_relactauto. */
  bool zero_activity_low_energy_retry = false;

  /** zero_activity_low_energy_retry judges a nuclide on the evidence the data show: its "principal
   anchor" - the peak the solve must still explain - is the match of its highest-yield line among all
   its significant matched search peaks, x-ray matches included, rather than among the clean matches
   only; and a lost principal does not make the nuclide "zeroed" while the solve fits an ROI of its
   own lines at least twice as significantly as the data show the principal.  R500 Pd103_Phantom: the
   search missed the 7000-count 22 keV Rh x-rays (on the turn-on), its weak 357 keV line (z 6.4) stood
   for Pd103, went unexplained while the solve fit the x-rays at z 25, and the retry dropped the x-ray
   ROI.  (Ranked by data z instead, R500 Th232_Sh's principal became its 182 keV backscatter hump,
   which no solve explains, and the re-solve that kept Th232 was refused.) */
  bool zero_activity_evidence_anchors = false;

  /** When zero_activity_low_energy_retry keeps its re-solve, the low-energy ROIs it dropped - for the
   solve's sake, the decay data's x-ray intensities or the scatter hump dragging the shared curve - are
   judged on their own data at the end (the final filter's local rescue, lines added strongest first)
   unless a final ROI covers them, so the strong peaks there are still delivered.  The retry lost SAM
   Xe133_Unsh's z=459 81 keV and z=403 31 keV lines (it zeroed the trace Xe133m, whose only anchor was a
   z=5 background line), NGH I123_Unsh's z=110 28 keV Te x-rays and SAM Au198_Unsh's z=44 70 keV Hg
   x-rays.  The lines are judged at the planner's widths (the retry's ran to many times the
   resolution: a 150 keV FWHM on SAM Au198_Unsh's 1-channel x-ray, its total 40x the data), each at
   z >= 5, and only those 1.5 FWHM clear of the detector turn-on (R500 Th232_Sh 19 keV and Lu177m_Sh
   27 keV were knee phantoms). */
  bool zero_activity_rescue_dropped_rois = false;

  /** Before the final filter deletes an ROI the solve left insignificant, refit that ROI's lines
   with free amplitudes over a quadratic continuum, and keep the ROI if they pass the same
   significance test.  The solve ties every line to one relative efficiency, so a region that
   curve mis-predicts can end with its visible peaks at a fraction of their counts (Br76 in a
   phantom: 1216 and 1854 keV deleted with the 559 keV region fitting worse than no peaks). */
  bool rescue_insignificant_rois = false;

  /** rescue_insignificant_rois adds the ROI's lines strongest first, each kept only if it improves
   the fit, instead of eliminating the weakest from all of them.  The local fit's chi2 holds the
   continuum at or above zero, so one line whose best amplitude needs a negative continuum (an edge
   line on a steep turn-on) spoils every set it is in, and the backward elimination then kept
   nothing: R500 I123_Phantom's z=130 159 keV line was dropped with its 98-186 keV ROI once a
   refinement change delivered that solve. */
  bool rescue_forward_selection = false;

  /** The local rescues (rescue_insignificant_rois, and final_filter_veto_min_lambda's confirmation) fit
   a LINEAR continuum under an ROI narrower than this many FWHM of its widest line, a quadratic
   otherwise (0 = always quadratic).  A quadratic absorbs 86 % of one peak's chi2 over 3 FWHM and 97 %
   over 2, so over a narrow ROI a real line cannot pass the local test (z 10-20 lines scored ROI z
   2.5-3: SAM At211_Unsh 570 keV, LaBr3 Ra226_Sh 1390, NGH Re188_Phantom 65 keV).  Applies to the
   confirmation only unless final_filter_rescue_linear. */
  double rescue_linear_below_num_fwhm = 0.0;

  /** Also use rescue_linear_below_num_fwhm in rescue_insignificant_rois (the final filter's rescue of
   ROIs the solve left insignificant), not only in final_filter_veto_min_lambda's confirmation of ROIs
   whose refit was significant.  Off: there a straight continuum under a curved slope or a Compton edge
   rescued phantoms (R500 Bi213_Phantom 50 keV on a smooth slope, SAM Ag110m_Unsh 1124 keV on the 1384
   keV Compton edge, R500 Np237_Sh a 180-239 keV ROI whose continuum dived to zero). */
  bool final_filter_rescue_linear = false;

  /** When the final filter's rescue of an ROI the solve left insignificant fails, try again with only
   the lines the automated search found a peak on (within sm_rescue_search_peak_num_fwhm of their
   FWHM), on the narrow-ROI linear continuum (rescue_linear_below_num_fwhm).  A quadratic absorbs
   ~97 % of one peak over 2 FWHM, so the final filter dropped NGH Pu238_Unsh's z=19 152.7 keV line
   (29 keV ROI, search peak at 154.9) once a solve sized it low; the search peak is the independent
   evidence a line on a smooth slope lacks - applied to every narrow ROI (final_filter_rescue_linear)
   the linear continuum rescued slope phantoms at the same local z as real lines.  Left off (2026-09-26):
   the search finds peaks on scatter humps and Compton shoulders too, and three of six slope phantoms
   came back (NGH At211_Sh 77 keV, Eu154_Unsh 90 keV, SAM Ag110m_Unsh 1118 keV) for one real line
   (R500 Lu176_Unsh 26 keV).  With the search peak also held to the resolution width
   (sm_rescue_search_peak_max_width_ratio), a full run changed only three spectra (+2 moderate truth
   lines): SAM As72_Unsh
   1213 keV better, R500 Lu176_Unsh's turn-on ROI mixed, NGH Eu154_Unsh's 90 keV (5x the escape line's
   area, on the scatter hump's rising edge) worse - still off. */
  bool final_filter_rescue_linear_at_search_peaks = false;

  /** Test each ROI's significance on the model the solve delivered for it.  The test (holding the
   peaks fixed, refitting the continuum, comparing to a continuum-only null) had two faults:
   - it took the peaks whose MEAN lies in the ROI, but the solve gives every ROI its own copy of each
     line whose tail reaches it, so a line near a neighbouring ROI was counted twice and a line
     centred just outside lost its tail; now the ROI's peaks are those on its own continuum.
   - a peak-CDF step continuum (FlatStepCDF / LinearStepCDF / BiLinearStepCDF) is built from the
     ROI's peaks, which fit_amp_and_offset only does for peaks it fits, so with every peak fixed the
     ROI was scored with no step; now the continuum is refit the way refit_roi_continuums made it.
   Either made adjacent or stepped ROIs read "worse than no peaks" when they fit well, which feeds
   the final filter, the zero-activity retry's health check and the refinement scores.  Off for
   HPGe until measured there. */
  bool roi_significance_delivered_model = false;

  /** An ROI the final filter keeps only on its strongest peak's significance although its peaks fit
   WORSE than the continuum-only null, where that test could see them - `lambda`, the chi2 their
   shape leaves against the null's quadratic, at least this (0 = off; needs
   roi_significance_delivered_model) - must pass the same test after the observable refit, or its
   refit peaks are dropped.  The strongest-peak path uses the solve's own amplitude, so it delivered
   ROIs the data do not support (NGH Sm153_Phantom 14-26 keV, worse by 76688).  Two measured traps:
   - a quadratic absorbs 86 % of one peak's chi2 in a 3 FWHM ROI and 97 % in a 2 FWHM one, so
     without the lambda guard "worse than no peaks" is no evidence (CZT Mo99_Sh's 1.5 FWHM 740 keV
     ROI); the ROIs the solve fit worse than no peaks held 160 strong truth lines against 41 extra
     peaks on the five scintillator sets;
   - judged on the SOLVE's peaks, those ROIs are mostly real lines the solve mis-sized or displaced,
     which the observable refit repairs: routing them to rescue_roi_locally lost R500 Yb169_Phantom's
     z=315 54.9 keV line and five other strong lines on the first 12 R500 exemplars.
   Only delivered peaks are judged, so the refinement's decisions are untouched. */
  double final_filter_veto_min_lambda = 0.0;

  /** Let the ROI refinement explore instead of stopping at its first challenger that scores worse:
   the loop keeps the best candidate so far, judges each challenger against it, plans each pass from
   the latest one whether or not it was the best, and delivers the best.  Stopping at the first worse
   challenger makes the result depend on luck: R500 I123_Phantom's first challenger scored 0.3 %
   worse, but planning from it reached a good fit of the z=130 159 keV line two passes later.  Off:
   with the default common-domain score the exploration only reshuffled R500 (+6 / -9 truth lines);
   it wants a score that compares any two candidates fairly (see refinement_whole_range_score). */
  bool refinement_keep_best = false;

  /** Judge refinement challengers over the whole analysis range, each ROI the final filter keeps
   modelled by its own peaks and every other channel by the shared SNIP continuum (see
   delivered_spectrum_chi2 in the source), rather than over the energy both candidates model; needs
   roi_significance_delivered_model.  It credits the ROIs a challenger adds (SAM Tl200_Sh's seven),
   but on NaI the SNIP continuum does not follow the broad scatter and backscatter structure, so an
   ROI putting peaks on it scores as progress.  With refinement_keep_best, over the five scintillator
   sets: +16 strong / +33 moderate truth lines and +17 extra peaks, but a visual review of the 102
   changed R500 fits called the result worse about twice as often as better (new ROIs on the Compton
   edge, scatter plateau, backscatter hump and turn-on; inverted step continua under x-ray clusters). */
  bool refinement_whole_range_score = false;

  /** With refinement_whole_range_score: the score's baseline is its own SNIP continuum with this
   half-window in FWHM (0 = the planner's 1.5 FWHM SNIP), and `refinement_score_snip_presmooth`
   channels of boxcar half-width.  The planner's SNIP clips NaI structure narrower than ~3 FWHM - the
   scatter plateau, the backscatter hump, the turn-on - so an ROI modelling it scored as progress;
   on R500 plot data a 0.75 FWHM window cut the baseline's chi2 per channel over those regions
   3-10x (La140_Sh plateau 294 -> 28, Np237_Unsh backscatter 197 -> 24).  A Compton edge is not
   helped (the order-2 clip rounds sharp corners), and a window much under 1 FWHM starts to keep
   part of a real narrow peak in the baseline. */
  double refinement_score_snip_window_fwhm = 0.0;
  int refinement_score_snip_presmooth = 3;

  /** With refinement_whole_range_score: an ROI's credit (the baseline's chi2 minus its model's over
   its channels) is capped at what its peaks gain over a quadratic continuum-only fit there.  Where
   a quadratic follows the data - broad structure, and the SNIP baseline's own low bias (it rides
   ~0.5-0.7 sigma per channel under plain NaI continuum) - an ROI then earns nothing for existing. */
  bool refinement_score_null_capped_credit = false;

  /** With refinement_score_null_capped_credit: the cap is what FREE amplitudes at the ROI's line
   positions gain over the quadratic continuum, not what the solve's own amplitudes do.  The solve's
   amplitudes can be far off a line the data show clearly (the ROIs it fits worse than no peaks),
   and capped by them R500 I123_Phantom's z=130 159 keV and Yb169_Phantom's z=116 177 keV lines
   earned no credit, so the refinement dropped them. */
  bool refinement_score_free_amplitude_cap = false;

  /** With refinement_whole_range_score: chi2 charged per free parameter of each delivered ROI (its
   continuum's parameters and the floating peaks in it; tied lines are not free), AIC-style (2),
   so an ROI around an invisible z 3-4 line, or a merge-versus-split near tie, does not win on noise
   (c83's review: empty ROIs at Cu67_Sh 390 and Bi213_Phantom 802 keV).  0 = none. */
  double refinement_score_complexity_charge = 0.0;

  /** With refinement_whole_range_score: an ROI counts as delivered only if it survives
   final_filter_veto_min_lambda's power veto, as the final filter judges it. */
  bool refinement_score_power_veto = false;

  /** Score refinement challengers over the energy both candidates model with each ROI's delivered
   model (see solution_chi2_over_segments in the source: the peaks on the ROI's own continuum, over
   exactly the shared channels) instead of the peaks centred in it over the continuum's whole range;
   needs roi_significance_delivered_model.  The corrected score still favours an ROI that fits
   exactly the shared range over a wider one holding a line at its edge (R500 Br76_Phantom's 559 keV)
   and cannot credit a line a challenger adds (Th228_Unsh's 2614 keV); its net effect over the
   scintillator corpora was noise (+5 of ~4000 truth lines; SAM -14). */
  bool refinement_delivered_segments = false;

  /** Each pass of the observable refit starts every peak from the solve's width, so the refit's
   +-15% width freedom stays relative to the solve's resolution curve.  Carried from pass to pass
   it compounded: a weak NaI line narrowed 95 -> 56 -> 42 -> 33 keV at 2 MeV, its ROI shrinking with
   it (the ROI bounds follow the peak width), ending as a spike a third of the detector resolution. */
  bool observable_width_from_solve = false;

  /** A flat-step continuum whose step the observable refit turned upside down - higher above the
   ROI's peaks than below them - is refit once on a straight continuum.  The step models a peak's own
   down-scattered photons, which only ever raise the continuum below the peak; the refit fits the step
   freely, and under a uranium x-ray cluster it ran the wrong way and sank ~200 counts below the
   live-time-scaled background (R500 U235_Unsh_0400 44-138 keV, 176 -> 378 counts across the ROI),
   with the broad peaks taking the difference. */
  bool observable_inverted_step_as_linear = false;

  /** When the observable refit finds several peaks of an ROI insignificant, drop only the least
   significant and refit the rest, rather than dropping them all.  Overlapping lines can each be
   insignificant while together they are overwhelming: five lines within 12 channels of a 12.5 keV
   per channel NaI (the 100 keV one holding 143000 counts +- 514000) were all dropped in one pass and
   the region lost every peak. */
  bool observable_backward_elimination = false;

  /** An observable refit that lost every peak of an ROI that arrived with a measured line "collapsed";
   with this, the ROI goes back to the peaks the SOLVE delivered, not to those the refit started this
   pass from.  After a first pass had refit an ROI's only line negative (and carried it, see
   observable_backward_elimination), the collapse kept that and delivered it at -1639 counts (NGH
   Th228_Unsh 776 keV; SAM Ra226_Sh 1572, Ra226_Unsh 1529 keV).  Dropping such an ROI instead lost NGH
   U233_Sh's real 730 keV line, which the final filter's confirmation rescues from the solve's peaks. */
  bool observable_collapse_restores_solve = false;

  /** Where an ROI holds a strong line of the background (its live-time-scaled area 5 sigma over the
   foreground there), do not judge the ROI on the observable refit: neither the final filter's
   confirmation (final_filter_veto_min_lambda) nor a valley split (observable_split_valley_fraction).
   The refit is of the gross spectrum, and cannot tell the background line's structure from the ROI's
   own peaks: LaBr3's intrinsic La138 1436 keV feature cost Cs134_Unsh its 1365 keV line (confirmation)
   and Sb124_Unsh / Tl200_Sh their 1325-1408 keV lines (splits).  Left off (2026-09-26): it restores those
   five LaBr3 lines, but NaI backgrounds hold strong lines in many ROIs (Pb x-rays, 186, 352, 662, 1460
   keV), and the solve's misfit there is often the very phantom the confirmation removes - NGH
   Pu239_Sh_100g's comb, Ca47_Sh's 493 keV and LaBr3 Np237_Sh's 86-131 keV phantoms came back.  A
   background-aware refit continuum is the real fix. */
  bool observable_skip_at_background_lines = false;

  /** Put the background's own lines into the observable refit, held fixed at their live-time-scaled
   areas, where a line shows in the foreground at this significance (scaled area over the square root
   of the foreground counts within +-1 FWHM of it); 0 disables.  The solve fits the background-subtracted
   spectrum, but the observable refit fits the gross foreground, so a background line inside an ROI that
   no model peak sits on was taken up by the refit's continuum and neighbouring peaks: LaBr3 Cs134_Unsh
   lost its 1365 keV line and Tl200_Sh its 1408 keV line beside the intrinsic La138 feature, and NGH
   Ca47_Sh its 489 keV line.  Only lines of the resolution's width (<= 1.4x the planner's model; the
   search also returns humps and blends) inside the ROI are held (a Gaussian falls short of the real tail
   of a line outside it), and a line within half a FWHM of a model peak is left out, so that peak keeps
   measuring the gross foreground as InterSpec's foreground peaks do (the Activity/Shielding fit subtracts
   the background's peaks itself).  The held lines are not delivered: the reported continuum carries
   them, as it carries neighbouring ROIs' tails.  Measured and left OFF (2026-09-26, c149-c150 review 1
   better / 19 worse in three batches): a polynomial continuum cannot carry a line's bump, so the drawn
   fit fell below the data at every held line (NGH's 662 keV Cs137, SAM/R500 K40 1460 keV) and the
   neighbouring peak overshot its apex (NGH I124 603 keV) - where the old refit's peak or continuum had
   followed the data.  Ready for use if the background's lines are ever delivered as peaks. */
  double observable_fixed_background_line_z = 0.0;

  /** Deliver the background's lines as peaks: those that show in the foreground at this significance
   (the same test as observable_fixed_background_line_z) and are of about the resolution's width join the
   model as ordinary peaks of the ROI that holds them - an ROI ending within a FWHM of one is taken past it,
   two ROIs either side of one are joined - so the observable refit fits them and they are drawn and
   returned, labelled by their parent (K40, Th232, Ra226, Cs137, ... - the nearest known background line)
   with user label "background" and norm_css_color.  A line on a model peak is left to that peak, which
   measures the gross foreground as before.  0 disables.  Holding the lines undelivered measured the
   neighbouring peaks better but left every background line undrawn, the fit far below the data there;
   ROI edges on the flanks of unfitted background lines (NGH's 662 keV Cs137, SAM's K40 1460 keV) were the
   most frequent severe continuum complaint of the 2026-09-26 reviews. */
  double observable_deliver_background_line_z = 0.0;

  /** With observable_deliver_background_line_z and no background spectrum, take as background lines the
   foreground's search peaks on the common field-background lines (K40 1461; the Th232 chain's 239, 583,
   911, 2615; the Ra226 chain's 352, 609, 1120, 1764 keV) that no model peak sits on. */
  bool observable_background_lines_without_background = false;

  /** Refit each ROI's observable peaks with the tails of neighbouring ROIs' peaks in the model (held
   fixed), then refit its continuum without them, as a hand fit's continuum carries a neighbour's
   tail.  Fit with only its own lines, an ROI that starts on a strong neighbour's flank inflated
   whatever line sat at its edge to absorb the tail, its continuum diving under that edge peak
   (Th232 62-108 keV, U235 135-225 keV, Ra223 308-489 keV on the R500). */
  bool observable_neighbour_tails = false;

  /** Keep an ROI's solved peaks when the observable refit sank its continuum at an edge below half the
   solve's continuum and below 60 % of the data there, where the solve's continuum was no more than
   20 % above the data.  A free linear continuum under several broad
   NaI peaks with no clean sideband is degenerate with them, and the refit traded the one for the
   other (Ra223 308-489 keV: continuum 5 counts under 220 of data, peaks inflated to fill it) while
   the solve, whose line amplitudes are tied together, had followed the data. */
  bool observable_keep_solve_on_collapse = false;

  /** Keep an ROI's solved peaks, unrefit, when it has fewer channels than its free continuum and
   peak parameters plus two.  With no degrees of freedom the refit returns every amplitude as
   insignificant: a 12.5 keV/channel NaI puts Xe133's 81 keV line (126000 counts, truth z=460) in a
   3-channel ROI, and the line was dropped. */
  bool observable_keep_solve_without_dof = false;

  /** In the observable refit, hold a Constant/Linear/Quadratic/Cubic continuum at the solve's values
   (drawn on the raw data with the solve's tied line amplitudes) and re-measure only the peaks.  A free
   continuum under a free-amplitude edge line is degenerate with it: the refit flattened the slope and
   grew the line to fill the gap (U235_Unsh_2000: a truth z=2.6 163 keV line at an ROI's low edge
   carried the Compton step under the 186 keV peak).  Step continua stay free (their step is built
   from the peaks). */
  bool observable_fixed_polynomial_continuum = false;

  /** In the observable refit, hold iodine escape companions (RelActCalcAuto::Options::iodine_escape_peaks)
   at the solve's amplitude, which is tied to their parent line's on the background-subtracted data.
   Freed on the raw spectrum, a companion at 70-90 keV took up the background's Pb K x-rays. */
  bool observable_hold_escape_peaks = false;

  /** Run the planner's data tests (detail::fixed_shape_peak_z: the data-peak admission gates, the
   refutation rescue) on the background-subtracted spectrum, and do not treat a search peak the
   background explains as an obstacle or as a peak that swamps a group.  The solve is fit to the net
   spectrum already; judged on the gross one, NGH's background Cs137 662 keV line zeroed I124's
   visible 722.8 keV line and swamped its 602.7 keV group. */
  bool planner_net_data_tests = false;

  /** When the observable refit's polynomial continuum sank more than 3 sigma of the data below the
   solve's at an ROI edge while a peak within a FWHM of that edge grew more than 2 sigma past its solved
   amplitude, re-measure the ROI's peaks by linear least squares against the solve's continuum.  The
   free continuum and the free edge amplitude are degenerate, and the refit traded one for the other
   (U235_Unsh_2000: a truth z=2.6 line at 163 keV carried the Compton step under 186 keV).  Holding
   every ROI's continuum (observable_fixed_polynomial_continuum) also exposed deficits the free
   continuum rightly followed. */
  bool observable_remeasure_on_edge_dive = false;

  /** In the observable refit, bound each peak's width around its own input width, rather than
   PeakFitLM's shared model, wherever the ROI's input widths differ by more than 20 %.  The shared
   model lets the highest-energy peak differ from the lowest by at most 20 %, so across a wide
   scintillator ROI, whose resolution can grow threefold from one end to the other, it forces equal
   widths (Lu177m_Unsh 85-452 keV: 15-42 keV FWHM refit as 21-25 keV at chi2/dof 600, and its 378
   and 418 keV lines dropped). */
  bool observable_independent_widths = false;

  /** In the observable refit, move an ROI edge only when the edge peak is gone (no peak within half a
   sigma of it), and judge the edges with the held escape companions as members.  Otherwise a refit
   that shifted the edge peak's mean by 0.1 keV re-drew the edge at 2.4 FWHM from it: Tl201_Unsh's
   x-ray ROI lost 19 keV of low sideband and its z=56 40 keV escape peak. */
  bool observable_edge_moves_on_removal_only = false;

  /** In the observable refit, no ROI side reaches more than this many FWHM beyond its outermost
   peak (0 = no limit): judged once, on the solve's peaks before the refit, in the planner's width
   model - which a solve whose widths ran away cannot inflate - with a peak combined from several
   lines reaching half its extra width further, to its outermost line.  The planner extends a side while the data stay consistent with a straight
   continuum, which over a sparse scintillator continuum is up to its cap: R500 sides ran a median
   1.7 FWHM and up to 4.3 (Ir192_Sh's 1061 keV ROI went 300 keV past the line) against 1.2 and 3.0
   for the same spectra fit by hand.  Capping the planner instead re-decided which refinement passes
   the solve accepted and lost more lines than it saved; trimming the delivered ROI does not. */
  double observable_max_side_fwhm = 0.0;

  /** In the observable refit, split an ROI between two consecutive peaks more than this many FWHM
   apart (0 = never), each part keeping its own continuum (and observable_max_side_fwhm).  A planner
   side extended into a weak neighbouring group, or a group of lines spread over its whole width,
   left two peaks 5 FWHM apart on one straight continuum (Br76_Phantom's 2358-2994 keV ROI held
   the 2391 and 2950 keV peaks); a hand fit splits NaI peaks more than ~3 FWHM apart. */
  double observable_split_gap_fwhm = 0.0;

  /** In the observable refit, also split an ROI at least observable_split_valley_min_roi_fwhm FWHM wide
   between two peaks at least 2 FWHM apart where the peaks' summed density falls below this fraction
   of the continuum (0 = never) - where the data come back down to the continuum between features.
   Heavily overlapping NaI lines leak into each other in the planner, so one ROI can hold several
   features with valleys between them and no gap wide enough for observable_split_gap_fwhm: R500
   Lu177m_Unsh 85-451 keV (14 FWHM, 5 groups; the data drop to the continuum at ~262 keV), Pu239_Unsh
   21-117 keV, Ag110m_Unsh 524-1014 keV.  Each part must keep a peak at observable_split_min_z, and the
   refit parts must fit their channels no worse than the whole ROI did, or it stays whole (SAM
   In111_Sh's 204-275 keV part put its step continuum far above the data). */
  double observable_split_valley_fraction = 0.0;
  double observable_split_valley_min_roi_fwhm = 6.0;
  /** observable_split_valley_fraction judges the valley by the DATA's excess over the refit continuum
   (where that excess is under 2 sigma), and between peaks only 1.2 FWHM apart, instead of by the
   model's peak density between peaks 2 FWHM apart: weak model lines fill the valleys of a chain of
   joined groups (R500 Lu177m_Unsh 85-451 keV: the data meet the continuum at ~262 and ~355 keV, but
   a z=4 272 keV line sits in the first).  Left off (2026-09-26): the refit continuum under such a joined
   ROI runs below the data throughout, so it still found no valley there, and it lost Pu239_Unsh's
   good 77 keV cut. */
  bool observable_split_valley_by_data = false;

  /** Model in the solve only the lines inside the energy span the ROIs cover (see
   RelActCalcAuto::Options::model_lines_outside_roi_span): a line below the lowest ROI rests on an
   extrapolated efficiency, and at the start of a solve that can predict orders of magnitude too much. */
  bool solve_lines_within_roi_span = false;

  /** Model the iodine K x-ray escape peaks of a NaI/CsI detector, 28.5 and 32.4 keV below each source
   line above the iodine K edge, at the fraction of the line's area the detector physics gives (9 %
   of a 60 keV line, 1.6 % at 122 keV); see RelActCalcAuto::Options::iodine_escape_peaks.  Only used
   when the detector type is PeakFitUtils::CoarseResolutionType::Low.  Unmodelled, a strong
   low-energy line's escape was either missed (Xe133's 52 keV, Cd109's 59 keV, Tl201's 41 keV) or
   absorbed by a neighbouring line (Am241's 31 keV escape pulled its weak 26/33 keV lines up). */
  bool iodine_escape_peaks = false;

  /** Model the pair-production single and double escape peaks (511 / 1022 keV below a strong line
   above 1.6 MeV) on detectors other than HPGe too, where HPGe always models them.  Only an escape
   the automated search found is modelled: a scintillator's broad escape peak overlaps neighbouring
   lines, and a free floating peak there with no data of its own would take their counts.
   Unmodelled, Tl208's 2103.5 keV single escape was missed at z=22-27 on SAM-Eagle Th228/U233 and
   z=9-15 on LaBr3, and Y88's 1325 keV single escape at z=10-21. */
  bool escape_peaks_non_hpge = false;

  // Threshold for initial peak significance filter before refitting for observable_peaks.
  // Significance = S / sqrt(S + B) with S = peak_amplitude * 0.7607 (fraction of Gaussian
  // area within +/-1 FWHM) and B the FITTED continuum integral over the same +/-1 FWHM
  // window - using the fitted continuum rather than the gross data means a neighboring
  // peak's counts no longer dilute this peak's significance.
  double observable_peak_initial_significance_threshold = 2.25;

  // Threshold for final peak significance of the fixed-shape (linear least-squares) observable fits
  // (measure_on_continuum, rescue_roi_locally); the observable refit uses observable_refit_min_z.
  // Significance = peak_area / peak_area_uncertainty
  double observable_peak_final_significance_threshold = 2.0;

  // Thresholds on the area / area-uncertainty of a peak fit with free shapes (automated-search peaks,
  // observable-refit peaks).  Until 2026-10 they were the same fields as the likelihood-ratio,
  // data-z and fixed-shape thresholds; they are kept apart to follow the peak fit's uncertainties.

  /** Minimum area/uncertainty of an automated-search peak used as evidence of a source line (manual
   rel-eff anchors, a refinement's found-peak confirmation, found-peak ROI re-seeding). */
  double found_peak_min_z = 3.0;

  /** A part of an observable ROI split at a valley must keep a peak of at least this area/uncertainty
   (see observable_split_valley_fraction). */
  double observable_split_min_z = 3.0;

  /** The observable refit drops peaks below this area/uncertainty. */
  double observable_refit_min_z = 2.0;

  /** Without a background spectrum (observable_background_lines_without_background), the
   area/uncertainty a foreground search peak on a common background line needs to be delivered (see
   observable_deliver_background_line_z). */
  double observable_background_line_fit_z = 5.0;

  // Single-pass ROI planner and its band quantities (see GammaClusteringSettings::use_roi_plan).
  bool use_roi_plan = true;   // HPGe default since 2026-09-06; the NaI defaults switch it off until tuned
  double share_always_fwhm = 3.0;
  double separate_always_fwhm = 4.0;
  double share_max_leak_fraction = 0.0;   // see GammaClusteringSettings::share_max_leak_fraction
  double share_valley_max_excess_fraction = 0.0;   // see GammaClusteringSettings
  bool widen_roi_without_room = false;              // see GammaClusteringSettings
  bool roi_gap_at_boundary = false;                 // see GammaClusteringSettings
  double share_min_side_channels = 0.0;             // see GammaClusteringSettings
  double share_min_roi_channels = 0.0;              // see GammaClusteringSettings
  bool roi_core_covers_every_group = false;         // see GammaClusteringSettings
  double share_max_valley_depth = 0.0;               // see GammaClusteringSettings
  double admission_min_data_z = 0.0;                 // see GammaClusteringSettings
  bool cont_constant_from_raw_sidebands = false;    // see GammaClusteringSettings
  bool admission_dominant_line_required = false;     // see GammaClusteringSettings
  double admission_data_evident_min_z = 0.0;         // see GammaClusteringSettings
  bool roi_floor_at_extent = false;                  // see GammaClusteringSettings
  double found_peak_match_num_fwhm = 0.5;
  double quad_min_width_fwhm = 8.0;
  double quad_min_continuum_counts = 1000.0;  // see GammaClusteringSettings::quad_min_continuum_counts
  double quad_min_curvature_z = 3.0;          // see GammaClusteringSettings::quad_min_curvature_z
  double found_peak_min_predicted_fraction = 0.25;
  double step_min_asym_z = 1.5;
  double step_min_fraction = 0.0;
  bool step_use_chi2_trial = false;
  double step_low_side_extra_fwhm = 0.5;
  double sibling_absence_max_ratio = 5.0;
  double sibling_absence_drf_slack = 0.4;
  double sibling_absence_shield_g_cm2 = 120.0;
  double sibling_absence_max_eff_ratio = 50.0; // see GammaClusteringSettings::sibling_absence_max_eff_ratio
  bool sibling_absence_robust_limits = false;   // see GammaClusteringSettings::sibling_absence_robust_limits
  // Tight ROIs seeded for matched automated-search peaks: only for peaks the fitted curve accounts
  // for (x-ray matches excepted - their yields are not trustworthy), sized by the width MODEL rather
  // than the search peak's own width.  Off keeps seeding every matched peak at its own width.
  bool seed_from_source_accounted_peaks = false;
  double data_detect_min_predicted_z = 0.0;
  double data_detect_min_data_z = 3.0;
  double snip_gate_max_window_fwhm = 0.0;  // see GammaClusteringSettings::snip_gate_max_window_fwhm
  double visible_line_min_z = 0.0;         // see GammaClusteringSettings::visible_line_min_z
  double refute_min_predicted_z = 0.0;     // see GammaClusteringSettings::refute_min_predicted_z
  double refute_max_gross_multiple = 2.0;
  double sideband_max_predicted_fraction = 0.0;  // see GammaClusteringSettings
  double roi_touch_split_min_fwhm = 0.0;   // see GammaClusteringSettings::roi_touch_split_min_fwhm
  double roi_touch_min_line_fwhm = 0.6;
  double max_shared_span_fwhm = 12.0;
  double roi_min_side_fwhm = 1.0;         // see GammaClusteringSettings::roi_min_side_fwhm
  double cont_constant_max_asym_z = 1.0;  // see GammaClusteringSettings::cont_constant_max_asym_z
  double cont_constant_max_counts = 100.0; // ... and the ROI must hold fewer counts than this (0 = never constant)
  double obstacle_exclusion_fwhm = 1.5;   // see GammaClusteringSettings::obstacle_exclusion_fwhm
  double obstacle_min_side_fwhm = 1.0;
  double obstacle_own_line_fwhm = 0.0;    // see GammaClusteringSettings::obstacle_own_line_fwhm
  int sub_extent_max_lines = 6;           // see GammaClusteringSettings::sub_extent_max_lines
  // Observable refit freedom: 0 = small (mean within 0.15 sigma, amplitude changes only slightly),
  // 1 = medium (mean within 0.5 sigma, amplitude free to move moderately).  The final per-peak
  // decision should rest on the data, not on the model amplitude the solve handed over: with the
  // small refinement, a weak line the rel-eff under-predicted (Eu152 564 keV at 17 counts against
  // ~70 in the data) can never reach the significance threshold.
  int observable_refit_level = 1;
  // The rel-eff curve's polynomial order is capped at (strong ROIs - rel_eff_order_strong_roi_margin),
  // where a strong ROI is one whose data show net counts at 5 sigma or better.  A curve of order n
  // carries n+1 shape terms, and against too few measured lines the solve trades the curve against
  // the activity: a two-line Pd103 fit with an order-2 curve landed on rel-eff ~ x^-7 with an
  // activity fifty thousand times too small.  Counting the activity as one of the parameters gives
  // margin 2, but the corpus prefers 1: the admission gate reads the curve, and a curve flattened
  // below what the physics supports predicts weak lines too low to plan (Gd153 151.7 keV, the Co60
  // 2505 keV sum peak).  <= 0 disables the cap.
  int rel_eff_order_strong_roi_margin = 1;
  // Also count, as evidence for the rel-eff order, the significant photopeak-shaped automated-search
  // peaks inside the ROIs (merged within 1.5 FWHM); the order cap then uses the larger count.  On a
  // scintillator one wide ROI can hold several well-measured lines, and counting only ROIs held a
  // spectrum with seven strong lines from 50 to 308 keV (Yb169 in a phantom) to a flat curve.
  bool rel_eff_order_count_found_peaks = false;
  // Use the rel-eff form/order the manual ladder chose (by AICc on the matched peaks) for the
  // RelActAuto solve instead of rel_eff_eqn_type/order, when at least this many sources share the
  // curve (0 disables).  A physical curve cannot bend by decades the way the polynomial did in the
  // four-source Trinitite fit (it zeroed Ba133 and the Eu152 122 keV peak) - but the physics
  // envelope now removes the contaminant matches that bent it, and the physical-model solve is far
  // slower and scored worse on the four-source Pu corpus problems, so this stays off.
  size_t auto_rel_eff_follow_manual_winner_min_sources = 0;   // off: the physical model made the four-source Pu fits worse and 10-30x slower (2026-09-06)
  // Energy-calibration freedom given to the RelActAuto solve (unless the caller passes
  // DoNotVaryEnergyCal): 0 = none, 1 = linear (offset + gain), 2 = non-linear (a free deviation
  // offset at every interior ROI).  Linear by default: per-ROI deviation offsets let each ROI slide
  // by more than a keV to absorb misassigned or NORM peaks (the four-source Trinitite fit fitted a
  // +-1.2 keV zig-zag), and the observable refit re-centres each peak anyway.
  int energy_cal_fit_type = 1;

  // Step continuum decision parameters
  // Minimum peak detection significance z = S_est/sqrt(S_est + B_est) to consider step continuum
  double step_cont_min_peak_significance = 40.0;
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
    settings.use_roi_plan = use_roi_plan;
    settings.share_always_fwhm = share_always_fwhm;
    settings.separate_always_fwhm = separate_always_fwhm;
    settings.share_max_leak_fraction = share_max_leak_fraction;
    settings.share_valley_max_excess_fraction = share_valley_max_excess_fraction;
    settings.widen_roi_without_room = widen_roi_without_room;
    settings.roi_gap_at_boundary = roi_gap_at_boundary;
    settings.share_min_side_channels = share_min_side_channels;
    settings.share_min_roi_channels = share_min_roi_channels;
    settings.roi_core_covers_every_group = roi_core_covers_every_group;
    settings.share_max_valley_depth = share_max_valley_depth;
    settings.admission_min_data_z = admission_min_data_z;
    settings.admission_dominant_line_required = admission_dominant_line_required;
    settings.cont_constant_from_raw_sidebands = cont_constant_from_raw_sidebands;
    settings.admission_data_evident_min_z = admission_data_evident_min_z;
    settings.iodine_escape_peaks = iodine_escape_peaks;
    settings.roi_floor_at_extent = roi_floor_at_extent;
    settings.found_peak_match_num_fwhm = found_peak_match_num_fwhm;
    settings.quad_min_width_fwhm = quad_min_width_fwhm;
    settings.quad_min_continuum_counts = quad_min_continuum_counts;
    settings.quad_min_curvature_z = quad_min_curvature_z;
    settings.found_peak_min_predicted_fraction = found_peak_min_predicted_fraction;
    settings.step_min_asym_z = step_min_asym_z;
    settings.step_min_fraction = step_min_fraction;
    settings.step_use_chi2_trial = step_use_chi2_trial;
    settings.step_low_side_extra_fwhm = step_low_side_extra_fwhm;
    settings.sibling_absence_max_ratio = sibling_absence_max_ratio;
    settings.sibling_absence_drf_slack = sibling_absence_drf_slack;
    settings.sibling_absence_shield_g_cm2 = sibling_absence_shield_g_cm2;
    settings.sibling_absence_max_eff_ratio = sibling_absence_max_eff_ratio;
    settings.sibling_absence_robust_limits = sibling_absence_robust_limits;
    settings.data_detect_min_predicted_z = data_detect_min_predicted_z;
    settings.data_detect_min_data_z = data_detect_min_data_z;
    settings.snip_gate_max_window_fwhm = snip_gate_max_window_fwhm;
    settings.visible_line_min_z = visible_line_min_z;
    settings.refute_min_predicted_z = refute_min_predicted_z;
    settings.refute_max_gross_multiple = refute_max_gross_multiple;
    settings.sideband_max_predicted_fraction = sideband_max_predicted_fraction;
    settings.roi_touch_split_min_fwhm = roi_touch_split_min_fwhm;
    settings.roi_touch_min_line_fwhm = roi_touch_min_line_fwhm;
    settings.max_shared_span_fwhm = max_shared_span_fwhm;
    settings.roi_min_side_fwhm = roi_min_side_fwhm;
    settings.cont_constant_max_asym_z = cont_constant_max_asym_z;
    settings.cont_constant_max_counts = cont_constant_max_counts;
    settings.obstacle_exclusion_fwhm = obstacle_exclusion_fwhm;
    settings.obstacle_min_side_fwhm = obstacle_min_side_fwhm;
    settings.obstacle_own_line_fwhm = obstacle_own_line_fwhm;
    settings.sub_extent_max_lines = sub_extent_max_lines;
    settings.low_energy_abs_floor = low_energy_abs_floor;
    settings.low_energy_skip_threshold_ramp = low_energy_skip_threshold_ramp;
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
    settings.use_roi_plan = use_roi_plan;
    settings.share_always_fwhm = share_always_fwhm;
    settings.separate_always_fwhm = separate_always_fwhm;
    settings.share_max_leak_fraction = share_max_leak_fraction;
    settings.share_valley_max_excess_fraction = share_valley_max_excess_fraction;
    settings.widen_roi_without_room = widen_roi_without_room;
    settings.roi_gap_at_boundary = roi_gap_at_boundary;
    settings.share_min_side_channels = share_min_side_channels;
    settings.share_min_roi_channels = share_min_roi_channels;
    settings.roi_core_covers_every_group = roi_core_covers_every_group;
    settings.share_max_valley_depth = share_max_valley_depth;
    settings.admission_min_data_z = admission_min_data_z;
    settings.admission_dominant_line_required = admission_dominant_line_required;
    settings.cont_constant_from_raw_sidebands = cont_constant_from_raw_sidebands;
    settings.admission_data_evident_min_z = admission_data_evident_min_z;
    settings.iodine_escape_peaks = iodine_escape_peaks;
    settings.roi_floor_at_extent = roi_floor_at_extent;
    settings.found_peak_match_num_fwhm = found_peak_match_num_fwhm;
    settings.quad_min_width_fwhm = quad_min_width_fwhm;
    settings.quad_min_continuum_counts = quad_min_continuum_counts;
    settings.quad_min_curvature_z = quad_min_curvature_z;
    settings.found_peak_min_predicted_fraction = found_peak_min_predicted_fraction;
    settings.step_min_asym_z = step_min_asym_z;
    settings.step_min_fraction = step_min_fraction;
    settings.step_use_chi2_trial = step_use_chi2_trial;
    settings.step_low_side_extra_fwhm = step_low_side_extra_fwhm;
    settings.sibling_absence_max_ratio = sibling_absence_max_ratio;
    settings.sibling_absence_drf_slack = sibling_absence_drf_slack;
    settings.sibling_absence_shield_g_cm2 = sibling_absence_shield_g_cm2;
    settings.sibling_absence_max_eff_ratio = sibling_absence_max_eff_ratio;
    settings.sibling_absence_robust_limits = sibling_absence_robust_limits;
    settings.data_detect_min_predicted_z = data_detect_min_predicted_z;
    settings.data_detect_min_data_z = data_detect_min_data_z;
    settings.snip_gate_max_window_fwhm = snip_gate_max_window_fwhm;
    settings.visible_line_min_z = visible_line_min_z;
    settings.refute_min_predicted_z = refute_min_predicted_z;
    settings.refute_max_gross_multiple = refute_max_gross_multiple;
    settings.sideband_max_predicted_fraction = sideband_max_predicted_fraction;
    settings.roi_touch_split_min_fwhm = roi_touch_split_min_fwhm;
    settings.roi_touch_min_line_fwhm = roi_touch_min_line_fwhm;
    settings.max_shared_span_fwhm = max_shared_span_fwhm;
    settings.roi_min_side_fwhm = roi_min_side_fwhm;
    settings.cont_constant_max_asym_z = cont_constant_max_asym_z;
    settings.cont_constant_max_counts = cont_constant_max_counts;
    settings.obstacle_exclusion_fwhm = obstacle_exclusion_fwhm;
    settings.obstacle_min_side_fwhm = obstacle_min_side_fwhm;
    settings.obstacle_own_line_fwhm = obstacle_own_line_fwhm;
    settings.sub_extent_max_lines = sub_extent_max_lines;
    settings.low_energy_abs_floor = low_energy_abs_floor;
    settings.low_energy_skip_threshold_ramp = low_energy_skip_threshold_ramp;
    return settings;
  }

  /** Return as soon as the ROIs are planned, skipping the solve entirely.  For evaluating planning
   changes across whole corpora quickly; the returned peak lists are empty. */
  bool stop_after_plan = false;

  /** Names of every scalar/enum/string field that can be read or set by name (development
   harnesses, parameter sweeps, config files).  The `phys_model_*` shielding vectors are excluded. */
  static std::vector<std::string> field_names();

  /** Sets one field from its string form: numbers for doubles/sizes, "true"/"false"/"1"/"0" for
   bools, and either the integer value or the symbolic name for enums.  Returns false (leaving the
   config unchanged) when the name is unknown or the value does not parse. */
  bool set_field( const std::string &name, const std::string &value );

  /** String form of one field (full precision for doubles; symbolic names for enums).
   Throws std::invalid_argument for an unknown name. */
  std::string get_field( const std::string &name ) const;

  /** Every field as a `name=value` entry, joined by `separator`. */
  std::string to_string( const std::string &separator = "\n" ) const;

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

  /** Disable automatic detection and nuisance co-fitting of strong unmodeled interferers.

   By default the fitter may transactionally add foreground-confirmed strong NORM nuisance
   nuclides near requested-source lines.  This option disables that entire automatic R6 path:
   no interferer discovery, warnings, nuisance curves, or nuisance ROIs are added.  Explicitly
   requested sources and the FitNormBkgrndPeaks options are unaffected.

   This is useful when a representative background spectrum is supplied, and for controlled
   tuning/evaluation runs that need to measure the requested-source fitter independently of the
   optional interferer model.
   */
  DisableAutoInterfererFit = 0x40,
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
 @param cancel_calc Optional cooperative cancel flag: when it becomes true, the current
        RelActAuto solve stops at its next check point and the fit returns (status UserCanceled
        or FailToSolveProblem).  Lets callers bound the wall time of a pathological fit.
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
  const std::shared_ptr<const PeakFitDetPrefs> &peak_fit_prefs,
  const std::shared_ptr<std::atomic_bool> cancel_calc = nullptr );

PeakFitResult fit_peaks_for_nuclides(
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::vector<RelActCalcAuto::SrcVariant> &sources,
  const std::vector<std::shared_ptr<const PeakDef>> &user_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &background,
  const std::shared_ptr<const DetectorPeakResponse> &drf_input,
  const Wt::WFlags<FitSrcPeaksOptions> options,
  const PeakFitForNuclideConfig &config,
  const std::shared_ptr<const PeakFitDetPrefs> &peak_fit_prefs,
  const std::shared_ptr<std::atomic_bool> cancel_calc = nullptr );
  

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
