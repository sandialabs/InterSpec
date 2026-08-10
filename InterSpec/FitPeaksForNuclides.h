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
#include <cstdint>
#include <limits>
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
    const double aicc_penalty );


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
  bool should_try_source_clean_recovery( size_t num_source_anchors,
                                         size_t num_preserved_anchors );

  /** Transactional source-clean acceptance: recover predicted anchors, preserve the incumbent's
   FWHM-distinct significant fitted source evidence, and improve the filtered data score.  A valid
   candidate may replace an incumbent for which no significant-ROI score was available. */
  bool should_accept_source_clean_challenger( bool solve_succeeded,
                                              size_t incumbent_preserved_anchors,
                                              size_t candidate_preserved_anchors,
                                              size_t incumbent_fitted_anchors,
                                              size_t candidate_fitted_anchors,
                                              double incumbent_score,
                                              double candidate_score );

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

  /** Exact-once snapshot assignment of an atom universe onto already channel-disjoint ROIs (the
   bridge at the two points where geometry lives in RelActCalcAuto::Options::rois, which cannot
   carry atoms).  Each atom goes to the single ROI containing its energy; ties resolve to the
   nearest-midpoint ROI then lowest index.  Atoms outside every ROI are returned in `unowned`
   (pre-existing floaters, not losses of this operation).  `rois` must be disjoint (dev-asserted).*/
  void assign_atoms_to_disjoint_rois(
    const std::vector<RoiAtom> &universe,
    const std::vector<RelActCalcAuto::RoiRange> &rois,
    std::vector<std::vector<RoiAtom>> &per_roi_atoms,
    std::vector<RoiAtom> &unowned_atoms );


  /** One accepted source-gamma group supplied to the R4 shadow boundary optimizer. */
  struct RoiBoundaryShadowGroup
  {
    double legacy_lower = 0.0;
    double legacy_upper = 0.0;
    std::vector<double> gamma_energies;
  };

  /** One interval in the shadow optimizer's proposed partition.  Shadow results never alter the
   production ROIs; they are diagnostics for paired evaluation and visual review. */
  struct RoiBoundaryShadowInterval
  {
    double lower = 0.0;
    double upper = 0.0;
    double legacy_lower = 0.0;
    double legacy_upper = 0.0;
    double width_fwhm = 0.0;
    size_t num_channels = 0;
    PeakContinuum::OffsetType continuum_type = PeakContinuum::OffsetType::Linear;
    double normalized_continuum_mismatch = 0.0;
    double interval_score = 0.0;
    PeakContinuum::OffsetType legacy_continuum_type = PeakContinuum::OffsetType::Linear;
    double legacy_normalized_continuum_mismatch = 0.0;
    double legacy_score = 0.0;
    size_t first_group = 0;
    size_t last_group = 0;
    std::vector<double> group_gamma_energies;
    std::vector<double> profile_energies;
    std::vector<double> profile_foreground;
    std::vector<double> profile_snip;
    std::vector<double> profile_continuum;
    std::vector<double> unmodeled_peak_energies;
    size_t unmodeled_peak_conflicts = 0;
    std::string reason;
  };

  struct RoiBoundaryShadowResult
  {
    bool valid = false;
    std::string stage;
    std::string fallback_reason;
    double legacy_total_score = 0.0;
    double proposed_total_score = 0.0;
    std::vector<RoiBoundaryShadowInterval> intervals;
  };

  /** Jointly optimize a non-overlapping ROI partition against the shared SNIP baseline.  This is
   shadow-only: callers may inspect the proposal, but production fitting retains its legacy ROIs. */
  RoiBoundaryShadowResult optimize_roi_boundaries_shadow(
    const std::vector<RoiBoundaryShadowGroup> &groups,
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const GlobalContinuumEstimate &global_continuum,
    const std::function<double(double)> &fwhm_at_energy,
    const std::vector<std::shared_ptr<const PeakDef>> &unfit_auto_peaks,
    double catastrophic_max_fwhm_width,
    double roi_core_num_fwhm = 1.5 );

  /** Drain shadow diagnostics produced on the calling thread since the previous drain. */
  std::vector<RoiBoundaryShadowResult> take_roi_boundary_shadow_diagnostics();

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

  /** R2 keep-gate classification seam for focused statistical tests.  Returns true only for a
   counts-floor-passing cluster in [0.7*keep_z, keep_z], i.e. one that the normal strict keep gate
   rejects but the bounded rescue pass may inspect. */
  bool is_marginal_keep_reject( double expected_counts, double significance, double keep_z );

#if( PERFORM_DEVELOPER_CHECKS )
  /** Developer-test controls for proving rescue causality and transactional exception fallback. */
  void set_bounded_rescue_enabled_for_test( bool enabled );
  void force_next_bounded_rescue_admission_failure_for_test();
  void force_next_bounded_rescue_evaluation_failure_for_test();
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

  // Provenance stage stamped on atoms minted at the keep-gate (diagnostics only; not serialized or
  // tuned).  Refinement re-clustering sets this to RefinementCluster.
  detail::RoiAtomAdmission cluster_admission_stage = detail::RoiAtomAdmission::InitialCluster;
};


// Result from RelActAuto peak fitting
struct PeakFitResult
{
  RelActCalcAuto::RelActAutoSolution::Status status;
  std::string error_message;
  std::vector<std::string> warnings;  // warnings that don't prevent success

  // Structural decisions made while constructing automatic ROIs, in solve-calibration order.
  std::vector<AutomaticRoiDecisionDiagnostic> automatic_roi_diagnostics;

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
