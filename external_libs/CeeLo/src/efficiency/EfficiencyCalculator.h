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

/// @file EfficiencyCalculator.h
/// @brief Main entry point for detector efficiency calculations.
///
/// Orchestrates the Monte Carlo simulation:
///   1. Emits photons from a source toward the detector
///   2. Traces each photon through attenuators to the detector
///   3. Transports the photon through the detector crystal
///   4. Tallies ε_FEP and ε_total
///   5. Computes statistical uncertainties
///
/// Supports point sources and extended sources (cylindrical, rectangular,
/// and Marinelli beaker).

#include "geometry/Geometry.h"
#include "geometry/SourceGeometry.h"
#include "transport/PhotonTransport.h"
#include "materials/Material.h"
#include "cascade/CascadeTypes.h"

#include <Eigen/Core>
#include <vector>
#include <cstdint>
#include <map>
#include <functional>
#include <chrono>
#include <utility>
#include <random>

namespace ceelo {

/// Reason why the simulation stopped.
enum class StopReason : uint8_t {
    MaxEvents,        ///< Reached max_events limit
    FepPrecision,     ///< Achieved target FEP relative precision
    TotalPrecision,   ///< Achieved target total efficiency relative precision
    WallTime,         ///< Exceeded max wall-clock time
    CpuTime           ///< Exceeded max summed per-thread CPU time
};

/// Snapshot passed to a ProgressCallback after a batch merge. The efficiency
/// fields mirror EfficiencyResult and use the same IS-variance estimator; all
/// four are zero until the first merged batch with nonzero summed weight.
struct Progress {
    uint64_t num_events = 0;                        ///< Total events simulated so far
    std::chrono::steady_clock::duration elapsed{};  ///< Wall-clock time since start
    double frac_complete = 0.0;      ///< Estimated fraction complete (0 to 1; exactly 1.0 on the final call)
    double fep_efficiency = 0.0;     ///< Current ε_FEP estimate
    double fep_uncertainty = 0.0;    ///< 1σ statistical uncertainty on ε_FEP
    double total_efficiency = 0.0;   ///< Current ε_total estimate
    double total_uncertainty = 0.0;  ///< 1σ statistical uncertainty on ε_total
};

/// Progress callback: invoked at most ~once per second during the run, plus
/// exactly once on completion with frac_complete = 1.0 and efficiency values
/// identical to the returned EfficiencyResult.
using ProgressCallback = std::function<void(const Progress&)>;

/// Termination criteria for a simulation.
/// Multiple criteria can be active simultaneously; whichever is met first stops the run.
struct TerminationCriteria {
    uint64_t max_events = 0;                 ///< 0 = no limit
    double target_fep_rel_precision = 0.0;   ///< e.g., 0.01 for 1%; 0 = disabled
    double target_total_rel_precision = 0.0; ///< e.g., 0.01 for 1%; 0 = disabled
    double max_wall_seconds = 0.0;           ///< 0 = no limit
    /// Cap on CPU time summed across worker threads (0 = no limit). Unlike
    /// max_wall_seconds this is invariant to machine load and to how many
    /// simulations run concurrently, so it is the preferred budget cap.
    double max_cpu_seconds = 0.0;
    uint64_t min_events = 10000;             ///< Minimum events before precision checks
};

/// Full simulation configuration
struct SimulationConfig {
    double energy_keV = 0.0;
    TerminationCriteria termination;
    unsigned num_threads = 0;              ///< Number of worker threads (0 = auto-detect from hardware)
    std::vector<float> energy_bin_edges;
    /// Progress callback, invoked at most once per second, plus exactly once
    /// on completion (payload identical to the returned EfficiencyResult).
    /// Called under an internal mutex — keep the callback lightweight.
    ProgressCallback progress_callback;
    /// Events per batch per thread. Each thread simulates this many events before merging
    /// results and checking termination. Smaller values give more responsive termination
    /// and progress updates; larger values reduce mutex contention.
    uint64_t batch_size = 10000;
    /// Master RNG seed. 0 (default) = seed from std::random_device (production
    /// behavior, non-deterministic). Nonzero = deterministic, reproducible run;
    /// combine with num_threads = 1 for a bit-exact result (multi-threaded runs
    /// stay reproducible in distribution but not bit-exact, since batches race).
    uint64_t seed = 0;
};

/// Diagnostics: decomposition of source-effect (full-mode) tallies by the
/// primary's source-transport outcome.
///   u = primary crossed the source geometry with ZERO interactions
///       (straight line, full energy, no secondaries) — the only class that
///       can meaningfully contribute FEP at E0;
///   s = at least one interaction (Compton/PE/PP/Rayleigh) — the
///       scattered/secondary class that dominates the wide-angle
///       scattered-in part of total efficiency.
/// Per-class event counts, summed event weights, and summed squared weights
/// for the FEP and any-deposit indicators. Dividing weight sums by
/// num_events_simulated gives the per-class efficiency terms. All zero
/// unless the run is full mode with source effects.
struct SourceEffectDiag {
    uint64_t n_u = 0;        ///< Events whose primary had zero source interactions
    uint64_t n_s = 0;        ///< Events whose primary interacted in source geometry
    double fep_w_u = 0.0, fep_w2_u = 0.0;
    double any_w_u = 0.0, any_w2_u = 0.0;
    double fep_w_s = 0.0, fep_w2_s = 0.0;
    double any_w_s = 0.0, any_w2_s = 0.0;

    /// Absorbed-primary escaped-electron channel: events whose primary was
    /// absorbed in the source geometry with no surviving secondary photons
    /// but with escaped recoil electrons. A subset of the s-class. (Before
    /// June 2026 these events were dropped by an early-out and contributed
    /// exactly zero, so this tally measures the fix's effect directly.)
    uint64_t n_e_only = 0;             ///< Events in the channel
    uint64_t n_e_only_electrons = 0;   ///< Escaped electrons on those events
    double e_only_energy_keV = 0.0;    ///< Summed electron KE at escape
    double e_only_any_w = 0.0, e_only_any_w2 = 0.0;  ///< Any-deposit tally
    /// FEP weight from this channel. A single escaped electron can never
    /// reach the FEP window (KE < E0), but rare multi-electron events can
    /// land in it as a continuum coincidence (measured ~1 event in 8k
    /// channel events for a 2 cm Pb shell at 662 keV). Expect ~0.
    double e_only_fep_w = 0.0;

    void merge(const SourceEffectDiag& o) {
        n_u += o.n_u;         n_s += o.n_s;
        fep_w_u += o.fep_w_u; fep_w2_u += o.fep_w2_u;
        any_w_u += o.any_w_u; any_w2_u += o.any_w2_u;
        fep_w_s += o.fep_w_s; fep_w2_s += o.fep_w2_s;
        any_w_s += o.any_w_s; any_w2_s += o.any_w2_s;
        n_e_only += o.n_e_only;
        n_e_only_electrons += o.n_e_only_electrons;
        e_only_energy_keV += o.e_only_energy_keV;
        e_only_any_w += o.e_only_any_w;
        e_only_any_w2 += o.e_only_any_w2;
        e_only_fep_w += o.e_only_fep_w;
    }
};

/// Result of a single-energy efficiency calculation.
struct EfficiencyResult {
    double full_energy_peak_efficiency;  ///< ε_FEP: P(deposit ALL energy)
    double fep_uncertainty;              ///< 1-sigma statistical uncertainty on ε_FEP

    double total_efficiency;             ///< ε_total: P(deposit ANY energy > 0)
    double total_uncertainty;            ///< 1-sigma uncertainty on ε_total

    uint64_t num_events_simulated;
    double wall_time_seconds;
    /// CPU time summed across worker threads (seconds). Load-invariant cost
    /// measure; ~= wall_time_seconds * effective thread count on an idle box.
    double cpu_time_seconds = 0.0;

    /// Fraction of primary photons that hit the max_compton_scatters cap and
    /// were force-absorbed (residual energy deposited locally). Diagnostic for
    /// the truncation bias: should be ~0 in the validated regime now that the
    /// default cap is 40. See A1 in the review report.
    double forced_absorption_fraction = 0.0;

    /// Why the simulation stopped.
    StopReason stop_reason = StopReason::MaxEvents;

    /// Optional: pulse-height distribution using caller-provided energy bins.
    std::vector<float> pulse_height_distribution;

    /// 1-sigma statistical uncertainty per bin, parallel to
    /// pulse_height_distribution: sigma_i = sqrt(Sum w_i^2) / sum_weights.
    /// Same length as pulse_height_distribution; empty if no bins requested.
    /// For analog (unweighted) sampling this recovers Poisson sqrt(count_i);
    /// for weighted sampling (cone, two-stream, Compton-mixture, forced
    /// interaction) per-bin sigma differs from sqrt(count_i) and cannot be
    /// reconstructed from the counts alone -- hence the companion Sum w^2.
    std::vector<float> pulse_height_uncertainty;

    /// Diagnostics: pair-production annihilation gammas generated in source
    /// material/shielding and processed against the detector, their summed
    /// (unweighted) scoring deposit, and the weighted count of events whose
    /// any-deposit status came ONLY from these secondaries (their summed
    /// event weight / num_events_simulated ≈ the channel's contribution to
    /// total efficiency). All zero unless source effects + PP occur.
    uint64_t pp_secondaries_processed = 0;
    double pp_secondary_deposit_keV = 0.0;
    double pp_only_any_weight = 0.0;

    /// Diagnostics: per-class (unscattered/scattered primary) decomposition
    /// of the source-effect tallies. See SourceEffectDiag.
    SourceEffectDiag src_diag;
};

/// Default energy-match tolerance (keV) for cascade-correction cache lookups.
/// Absorbs float round-trip / rounded-literal differences between the energies
/// InterSpec queries and the energies the MC ran at, while staying well below
/// typical gamma/x-ray line spacing. Callers may override per call.
inline constexpr double kCascadeEnergyMatchTolKeV = 0.05;

/// A gamma that is coincident with a primary gamma (provided by InterSpec).
struct CoincidentGamma {
    double energy_keV;             ///< energy of the coincident gamma
    double coincidence_fraction;   ///< P(this gamma emitted | primary gamma was emitted)
};

/// Result of cascade correction for a single primary gamma.
struct CascadeCorrectionResult {
    double summing_out_factor;  ///< C_out = Π_j [1 - f_j × ε_total(E_j)]
    double summing_in_term;     ///< C_in = Σ f_ab × ε_FEP(a) × ε_FEP(b) / ε_FEP(primary)
    double net_correction;      ///< C_net = C_out + C_in
    double corrected_fep;       ///< ε_FEP(primary) × C_net
    /// 1σ statistical uncertainty on corrected_fep, by first-order propagation
    /// of the cached per-energy FEP/total uncertainties and the primary FEP
    /// uncertainty. Assumes the contributing MC estimates are independent
    /// (separate runs); reusing one cache entry across terms makes this mildly
    /// conservative. Zero when no input uncertainties are supplied.
    double corrected_fep_uncertainty = 0.0;
    /// Energies requested (coincident gammas and summing-in pair members) that
    /// were NOT found in the efficiency cache within the match tolerance. A
    /// non-empty list signals a mis-keyed cache rather than "no summing" — the
    /// caller should treat it as an error, not silently trust the result.
    std::vector<double> unmatched_energies;
};

/// A pair of coincident gammas whose energies sum to (or near) a primary gamma
/// energy, contributing a summing-in effect to the primary photopeak.
struct SummingInPair {
    double energy_a_keV;    ///< Energy of first gamma in the pair (keV)
    double energy_b_keV;    ///< Energy of second gamma in the pair (keV)
    double joint_fraction;  ///< P(this pair emitted | primary gamma emitted)
};

// ===================== Gold-standard correlated-cascade MC =====================
//
// compute_cascade() evaluates true-coincidence summing for each requested peak
// by Monte-Carlo, with full transport physics (partial deposits, x-ray summing,
// 511 pairs, angular spread). For each peak it conditions on that peak's gamma
// being emitted, samples its coincident partners from the level-scheme / pairwise
// coincidence data, emits all of them from the source vertex, transports each,
// and compares the SUMMED crystal deposit (with summing) against the primary-
// alone deposit (no summing). The peak's own gamma is cone-biased toward the
// detector for efficiency; partners are emitted isotropically (unbiased). Point
// and extended sources are supported, with optional source material/shielding
// (each member is transported through the source geometry, including its 511/
// brems secondaries); extended or shielded sources emit isotropically (cone off).
//
// The with-summing score per history is X = I_sum + Σ_k w_k·P_k: the direct
// indicator I_sum (peak gamma + its sampled partners + coincident vacancy x-rays
// in the window) plus the analytic "sum-peak-fed" pair channels. Each channel k
// is a gamma pair (a,b) with E_a+E_b in the window where neither is the peak
// gamma; w_k is its per-primary-emission weight from the exact level-scheme joint
// probability minus the triple overlap the direct stream already samples:
//   w_k = [ Σ_B bw_B·P_B(t_a∧t_b)·pγ_a·pγ_b
//           − (B is the primary branch ? bw_A·P_A(t_a∧t_b∧t_p)·pγ_a·pγ_b·pγ_p : 0) ]
//         / (bw_A·pass_p·pγ_p),
// and P_k is the MC indicator that a (cone-biased, like the primary) + b + the
// pair's conditioned partners/vacancies land in the window. Cone-biasing photon a
// is unbiased because both pair photons are required to exceed 2·tol (see
// cascade_sum_pair_channels), so the window is unreachable when a misses.
//
// This path is fully opt-in and shares no code with simulate_thread(), so
// non-cascade results and performance are unaffected. FullRealization is the
// reference; see the CascadeMethod::Conditional note for the remaining
// approximations (partner independence, x-ray-fed / triple-fed sums not enumerated).

/// How compute_cascade() evaluates the summing.
enum class CascadeMethod {
    /// Per-peak conditional estimator: condition on the peak's gamma, sample its
    /// coincident partners, cone-bias the peak gamma toward the detector. Exact
    /// for the pairwise data and efficient (cone-biased, so best at far geometry).
    /// Partner emission probabilities come from the daughter level scheme
    /// (visit/reach DP, `cascade_level_pmate`) when available, so mutually-
    /// EXCLUSIVE transitions (alternative de-excitations of a shared level) are
    /// correctly not summed; only when no level scheme is present does it fall
    /// back to the pairwise/marginal `cascade_partner_prob`.
    ///
    /// Captures: (1) summing-OUT by the peak's coincident gammas and by the
    /// coincident daughter K/L vacancy x-rays (EC feed + internal conversion,
    /// via `cascade_level_vacancies` — the same vacancy model as FullRealization,
    /// so the default method is no longer γ-only); (2) summing-IN involving the
    /// peak's own partners; and (3) "sum-peak-fed" summing-IN — pairs (a,b) whose
    /// energies add into the window where NEITHER is the peak's gamma (e.g. Co-57
    /// 122.06+14.41 → 136.47), added analytically as weighted MC pair channels
    /// (`cascade_sum_pair_channels`) with the exact level-scheme joint probability.
    ///
    /// Remaining approximations (use FullRealization as the reference when these
    /// matter at close geometry): non-primary partners are sampled independently
    /// from exact conditional marginals rather than as one conditioned DAG path.
    /// In a 3+-gamma cascade that can co-emit two mutually exclusive downstream
    /// alternatives or lose sequential correlation; it is an approximation, not
    /// a physical decay realization. Sum-peak feeding
    /// enumerates gamma-gamma pairs only — x-ray-FED feeding (γ + coincident K
    /// x-rays adding into the window, e.g. Co-57 122 + two Fe Kα) and triple-fed
    /// sums are not enumerated (≈1–2% on Co-57 136 / Ba-133 at contact, ε²-small
    /// at far geometry).
    ///
    /// A window containing more than one positive-yield primary gamma cannot be
    /// conditioned on one member without biasing the line-area correction. Such
    /// ambiguous windows use an isotropic one-peak FullRealization sub-run with
    /// the same event and thread counts. This is unbiased, but can have much
    /// higher variance than cone-biased Conditional sampling at far geometry.
    Conditional,
    /// Full-realization estimator: per history, sample one complete decay
    /// realization (branch + correlated emitted set), emit all members
    /// isotropically, and tally every peak window plus the summed spectrum.
    /// For accepted level graphs, captures all represented channels (summing-
    /// out, summing-in including sum-peak-fed lines, and sum continuum). A
    /// rejected branch instead uses a marginal-preserving one-parent projection
    /// of its inconsistent pair links; affected peaks set
    /// `summing_model_complete=false` because that fallback cannot honor every
    /// tabulated pair simultaneously. Isotropic, so best at close geometry.
    FullRealization
};

/// Per-peak true-coincidence-summing result.
struct PeakCascadeResult {
    double energy_keV = 0.0;
    bool   found = false;               ///< a matching gamma member was located in `cascades`
    double eff_no_summing = 0.0;        ///< single-gamma FEP efficiency ε_FEP(E) (no coincidences)
    double eff_no_summing_unc = 0.0;
    double eff_with_summing = 0.0;      ///< FEP efficiency WITH coincidence summing
    double eff_with_summing_unc = 0.0;
    /// Summing factor = eff_with_summing / eff_no_summing. < 1 => net summing-out
    /// (peak loses counts); > 1 => net summing-in (peak gains). The literature
    /// TCS correction factor k_TCS is typically 1 / summing_factor.
    double summing_factor = 1.0;
    double summing_factor_unc = 0.0;
    /// False when this requested peak materially depends on a rejected branch.
    /// Conditional then omits invalid-branch sum-fed pairs; FullRealization uses
    /// a coherent one-parent Fréchet projection which preserves marginals but
    /// necessarily discards some inconsistent pair constraints. The numeric
    /// estimate remains available for compatibility in either case.
    bool summing_model_complete = true;
};

/// Snapshot passed to a CascadeProgressCallback during compute_cascade().
struct CascadeProgress {
    /// Histories completed so far. For Conditional this counts across peaks
    /// (total = num_events × number of found peaks); for FullRealization it
    /// is the overall history count (total = num_events).
    uint64_t num_events = 0;
    std::chrono::steady_clock::duration elapsed{};  ///< Wall-clock time since start
    double frac_complete = 0.0;  ///< Overall fraction complete (0 to 1; exactly 1.0 on the final call)
    /// Running per-peak snapshot, same order and struct as CascadeResult::peaks.
    /// Conditional: completed peaks carry their final values, the in-flight
    /// peak its running estimate, not-yet-started peaks zeros (found flag set);
    /// not-found peaks stay {energy, found = false}.
    std::vector<PeakCascadeResult> peaks;
};

/// Cascade progress callback: invoked at most ~once per second during the run,
/// plus exactly once on completion with frac_complete = 1.0 and peaks identical
/// to the returned CascadeResult::peaks. Called under an internal mutex — keep
/// the callback lightweight.
using CascadeProgressCallback = std::function<void(const CascadeProgress&)>;

/// Input to compute_cascade(): the nuclide's correlated emissions + the peaks
/// to report. Build `cascades` with the SandiaDecay adapter or by hand.
struct CascadeConfig {
    std::vector<DecayCascade> cascades;  ///< correlated decay branches
    std::vector<PeakWindow>   peaks;     ///< gamma peaks to report a corrected efficiency for
    uint64_t num_events = 1000000;       ///< MC histories (per peak for Conditional; total for FullRealization)
    unsigned num_threads = 0;            ///< 0 => hardware_concurrency()
    CascadeMethod method = CascadeMethod::Conditional;
    /// Lower-edge energy bins for the summed spectrum (FullRealization only).
    std::vector<float> spectrum_bin_edges;
    /// Progress callback, invoked at most once per second (and always on completion).
    CascadeProgressCallback progress_callback;
    /// Deposit internal-conversion (and coincident K-Auger) electron kinetic
    /// energy for converted transitions, in addition to the vacancy x-ray. Off
    /// by default (photon-only, bit-identical to the legacy path). When on, the
    /// electron is transported through the source geometry + air gap and
    /// deposited in the crystal only if it can physically reach it (contained /
    /// air-absorbed / off-solid-angle electrons contribute nothing), so distant
    /// external sources are unaffected. Matters for contact / well / 4π sources.
    bool enable_ic_electrons = false;
};

struct CascadeResult {
    std::vector<PeakCascadeResult> peaks;
    uint64_t num_events_per_peak = 0;
    /// Summed pulse-height spectrum (counts per parent decay per bin) and its
    /// 1σ. Populated only for CascadeMethod::FullRealization with spectrum bins.
    std::vector<float> summed_spectrum;
    std::vector<float> summed_spectrum_uncertainty;
    /// True when a populated spectrum used accepted level graphs for every
    /// positive-weight branch. False means at least one sampled branch used the
    /// marginal-preserving fallback forest. Remains true when no spectrum was
    /// requested or produced.
    bool summed_spectrum_model_complete = true;
};

/// How to handle air attenuation in the gap between source and detector.
enum class AirAttenuation : uint8_t {
    None,                ///< No air attenuation (default, current behavior)
    AnalyticNoScatter    ///< exp(-mu_no_rayleigh * air_distance) per event;
                         ///< excludes Rayleigh (isotropic cancellation argument)
};

/// Configuration for position-based importance sampling in extended sources.
///
/// Biases source position sampling toward regions that contribute the most
/// signal: near the detector axis (lateral) and/or near the detector-facing
/// surface (depth). Each event's weight is corrected by p_physical / p_biased
/// to preserve unbiased results.
///
/// Parameters are in detector-frame terms. Internally, the projection onto
/// the source local frame is handled automatically (accounting for rotation
/// and offset).
///
/// Only effective for Cylindrical and Rectangular source types.
struct PositionBiasConfig {
    double lambda_lateral_cm = 0.0;  ///< Lateral bias scale (cm); 0 = no lateral bias
    double lambda_depth_cm = 0.0;    ///< Depth bias scale (cm); 0 = no depth bias

    /// Solid-angle-matched lateral sampling for cylindrical sources.
    /// When > 0, samples r from p(r) ∝ r/(r²+h²) where h is this distance.
    /// With cone direction sampling, this produces approximately constant
    /// total weight for every hit event — desirable because it minimizes
    /// the variance from rare high-weight events. Overrides lambda_lateral.
    double solid_angle_h_cm = 0.0;

    /// Precomputed ln(1 + R²/h²) for solid-angle-matched sampling.
    /// Avoids recomputing log() per event. Set by compute_auto_bias_params().
    double log_solid_angle_base = 0.0;

    /// Full-surface detector sampling: sample a target point on the detector's
    /// outer surface (front disk + side cylinder + back disk) instead of cone
    /// sampling. Eliminates per-position cone computation and Rodrigues rotation.
    /// Every accepted event hits the detector; rejected events (cos_inc <= 0,
    /// i.e. surface faces away from source) are skipped cheaply.
    bool use_surface_sampling = false;
};

/// Event-biasing (variance reduction) configuration.
///
/// All techniques are mathematically unbiased: each one multiplies the
/// per-event score weight; the estimator denominator stays N (number of
/// events). Defaults preserve analog behavior.
struct BiasingConfig {
    /// Force the first interaction of the primary photon inside the detector
    /// geometry. Removes the pass-through branch; the history carries weight
    /// 1 - exp(-tau) <= 1. Biggest gains at high energy and for thin
    /// detectors where the analog interaction probability is small.
    /// Not applied to Marinelli sources (the un-collided pass-through
    /// branch feeds the re-entry physics there).
    bool force_detector_interaction = false;

    /// Mixture angular biasing for full-mode runs with source effects
    /// (source material / source shielding): emit with pdf
    ///   q(Omega) = alpha * (1/4pi) + (1-alpha) * (1/Omega_cone) inside the
    /// detector-subtending cone, weight w = (1/4pi)/q(Omega) <= 1/alpha.
    /// Unbiased even with source scatter because q > 0 over all 4pi.
    /// 0 = disabled (pure isotropic emission). Valid range (0, 1].
    /// Superseded by two_stream when that is enabled.
    double mixture_cone_alpha = 0.0;

    /// Two-stream stratified estimator for full-mode source-effect runs
    /// (non-Marinelli). Every tally partitions exactly by the primary's
    /// source-transport outcome:
    ///   u-term: primary crosses the source geometry with ZERO interactions
    ///           (sole meaningful FEP contributor + direct part of total);
    ///   s-term: primary interacted (the wide-angle scattered-in component
    ///           that dominates total efficiency under shielding).
    /// A deterministic round-robin fraction f of events runs the DIRECT
    /// stream: cone emission with weight (Omega/4pi) * T, where
    /// T = exp(-sum mu_tot*l) is the deterministic zero-interaction
    /// probability through the source geometry, then plain detector
    /// transport at the source energy. The remaining events run the
    /// SCATTER stream: isotropic emission + analog source transport,
    /// killed early if the primary did not interact (that term belongs to
    /// the direct stream), else the full legacy detector processing.
    /// Stream weights carry 1/f and 1/(1-f) so the standard single-pool
    /// estimator sum(w_i I_i)/N stays unbiased for every tally including
    /// the spectrum. Replaces mixture_cone_alpha (no 1/alpha weight
    /// outliers on the scattered-in component).
    bool two_stream = false;

    /// Fraction f of events assigned to the direct stream (quantized to
    /// 1/20). Clamped to [0.05, 0.95]: both streams must keep events since
    /// each owns one term of the estimator. Larger f favors FEP precision;
    /// smaller f favors total/spectrum precision.
    double two_stream_direct_fraction = 0.5;

    /// Compton-angle mixture biasing of the primary's first two Compton
    /// vertices in the source geometry (scatter stream / source-effect
    /// runs): with probability gamma the outgoing direction is drawn
    /// uniformly in the cone subtending the detector from the vertex,
    /// otherwise from the analog KN x S(x,Z) kernel; the event weight gains
    /// p_analog/q <= 1/(1-gamma) per biased vertex. Feeds the scattered-in
    /// component of total efficiency under thick shields. 0 = disabled.
    /// Valid range [0, 0.9]. See SourceGeometry::ComptonBiasConfig.
    double compton_cone_fraction = 0.0;
};

/// Source type
enum class SourceType : uint8_t {
    Point,
    Cylindrical,
    Rectangular,
    Marinelli,
    Spherical
};

/// Depth distribution for extended sources.
/// Uniform: activity uniformly distributed throughout the volume (default).
/// Exponential: activity follows A(d) = exp(-d/lambda) where d is depth from
///              the detector-facing surface and lambda is the relaxation length.
enum class DepthDistribution : uint8_t {
    Uniform,
    Exponential
};

/// Configuration for a cylindrical extended source.
///
/// The cylinder axis is along the source local z-axis.  The cylinder spans
/// [-half_length, +half_length] in z and has the given radius in the source
/// local frame.
///
/// The rotation matrix maps from the detector coordinate frame to the source
/// local frame:
///   local_vec = rotation × detector_vec
///   detector_vec = rotation.transpose() × local_vec
///
/// An identity rotation means the source axis is aligned with the detector axis
/// (coaxial configuration), which is the most common case.
struct CylindricalSourceConfig {
    Eigen::Vector3d center{0.0, 0.0, -10.0};  ///< Source center in detector frame (cm)
    double radius{1.0};                         ///< Source cylinder (outer) radius (cm)
    double inner_radius{0.0};                   ///< Inner bore radius (0 = solid); annular tube/pipe/ring when > 0
    double half_length{1.0};                    ///< Source half-length along its axis (cm)
    Eigen::Matrix3d rotation{Eigen::Matrix3d::Identity()};  ///< detector → source local
};

/// Configuration for a spherical extended source.
///
/// The active material occupies the shell [inner_radius, radius]; inner_radius = 0
/// is a solid ball, inner_radius > 0 a hollow spherical shell (non-attenuating
/// void center). The rotation is stored for API symmetry but is physically
/// irrelevant — a sphere is rotation-invariant.
struct SphericalSourceConfig {
    Eigen::Vector3d center{0.0, 0.0, -10.0};  ///< Source center in detector frame (cm)
    double inner_radius{0.0};                   ///< Inner void radius (0 = solid ball)
    double radius{1.0};                         ///< Outer radius (cm)
    Eigen::Matrix3d rotation{Eigen::Matrix3d::Identity()};  ///< detector → source local (unused)
};

/// Configuration for a rectangular (box) extended source.
///
/// The box spans [-hx, +hx] × [-hy, +hy] × [-hz, +hz] in the source local frame.
/// The rotation matrix maps from the detector coordinate frame to the source local frame.
///
/// A nonzero inner_half_dims makes the source a hollow box shell (crate,
/// container wall): the active material is the outer box minus the inner box
/// (both centered, same rotation); the inner box is an inactive,
/// non-attenuating void. All-zero = solid box.
struct RectangularSourceConfig {
    Eigen::Vector3d center{0.0, 0.0, -10.0};        ///< Source center in detector frame (cm)
    Eigen::Vector3d half_dims{1.0, 1.0, 1.0};        ///< Half-dimensions (hx, hy, hz) in cm
    Eigen::Vector3d inner_half_dims{0.0, 0.0, 0.0};  ///< Inner void half-dims (all-zero = solid)
    Eigen::Matrix3d rotation{Eigen::Matrix3d::Identity()};  ///< detector → source local
};

/// Configuration for a Marinelli beaker source.
///
/// A Marinelli beaker is a re-entrant sample container that surrounds the
/// detector, providing nearly 4π geometry coverage. The sample fills an
/// L-shaped volume: a disk in front of the detector plus an annular ring
/// around the detector sides.
///
/// All distances are measured relative to the **detector geometry front face**
/// (z_det_min = geometry_.outer_z_extent().first), which includes dead layers
/// and attenuator front thicknesses — NOT the crystal front face (z=0).
///
/// - endcap_to_beaker: air gap from detector geometry front to beaker well opening.
/// - well_depth: how far the well extends past the detector geometry front,
///               measured from z_det_min in the +z direction (over detector sides).
/// - fill_height: height of sample disk in front of the beaker well opening.
///
/// The detector must be fully configured (set_detector, set_dead_layer,
/// add_attenuator) BEFORE calling set_marinelli_beaker().
struct MarinelliBeakerConfig {
    double well_inner_radius = 0.0;   ///< Inner radius of re-entrant well (cm); must clear detector geometry
    double well_depth = 0.0;          ///< How far well extends past detector geometry front (cm)
    double outer_radius = 0.0;        ///< Outer radius of beaker body (cm)
    double fill_height = 0.0;         ///< Height of sample disk in front of beaker opening (cm)
    double endcap_to_beaker = 0.0;    ///< Air gap: detector geometry front to beaker well opening (cm)
    const Material* beaker_material = nullptr;  ///< Container wall material (e.g., polyethylene)
    double beaker_thickness = 0.0;    ///< Container wall thickness (cm)
};

/// Main calculator class.
class EfficiencyCalculator {
public:
    EfficiencyCalculator();

    // --- Geometry Configuration ---

    void set_detector(DetectorShape type, const Material* material,
                      const std::vector<double>& dimensions);
    void set_bore_hole(double bore_radius, double bore_depth);
    void set_dead_layer(double front_thickness, double side_thickness,
                        double back_thickness = 0.0);
    void add_attenuator(const Material* material,
                        double front_thickness, double side_thickness,
                        double z_start, double z_end);

    /// Add a collimator (tube-shaped, side walls only, no front/back disks).
    void add_collimator(const Material* material, double side_thickness,
                        double z_start, double z_end);

    // --- Source Configuration ---

    /// Configure a point source at the given position (detector frame, cm).
    void set_point_source(const Eigen::Vector3d& position);

    /// Configure a cylindrical extended source.
    ///
    /// **Orientation (maps to InterSpec's CylinderEndOn / CylinderSideOn).**
    /// The cylinder axis is the source local z-axis. The `rotation` (detector →
    /// source-local) sets the orientation relative to the detector axis (+z):
    ///   - **EndOn** (detector looks at a flat end): identity rotation — the
    ///     cylinder axis is parallel to the detector axis (the default, coaxial).
    ///   - **SideOn** (detector looks at the curved side): a 90° rotation taking
    ///     the cylinder's local z onto a detector transverse axis, e.g. about the
    ///     y-axis  R = [[0,0,1],[0,1,0],[-1,0,0]]  (local z ← detector x).
    ///
    /// @param center       Center of the source cylinder in the detector frame (cm).
    /// @param radius       Source cylinder (outer) radius (cm).
    /// @param half_length  Source half-length along the cylinder axis (cm).
    /// @param rotation     Rotation from detector frame to source local frame.
    ///                     Identity means source axis is coaxial with detector axis.
    /// @param inner_radius Inner bore radius (cm) for a hollow/annular cylinder
    ///                     (tube, pipe, ring); 0 = solid. The bore is an inactive,
    ///                     non-attenuating void.
    void set_cylindrical_source(const Eigen::Vector3d& center,
                                double radius, double half_length,
                                const Eigen::Matrix3d& rotation = Eigen::Matrix3d::Identity(),
                                double inner_radius = 0.0);

    /// Configure a spherical extended source (solid ball or hollow shell).
    ///
    /// @param center        Center of the source sphere in the detector frame (cm).
    /// @param radius        Outer radius (cm).
    /// @param rotation      Stored for API symmetry; physically irrelevant (a
    ///                      sphere is rotation-invariant).
    /// @param inner_radius  Inner void radius (cm) for a hollow spherical shell;
    ///                      0 = solid ball. The void center is non-attenuating.
    void set_spherical_source(const Eigen::Vector3d& center, double radius,
                              const Eigen::Matrix3d& rotation = Eigen::Matrix3d::Identity(),
                              double inner_radius = 0.0);

    /// Configure a rectangular (box) extended source.
    ///
    /// @param center    Center of the source box in the detector frame (cm).
    /// @param half_dims Half-dimensions (hx, hy, hz) in the source local frame (cm).
    /// @param rotation  Rotation from detector frame to source local frame.
    ///                  Identity means source box edges are aligned with detector axes.
    /// @param inner_half_dims  Inner void half-dimensions for a hollow box
    ///                  shell (crate, container wall); all-zero = solid box.
    ///                  The active material is the outer box minus the inner
    ///                  box (both centered, same rotation); the void is
    ///                  inactive and non-attenuating. Must satisfy
    ///                  0 <= inner < outer componentwise, or be all zero.
    void set_rectangular_source(const Eigen::Vector3d& center,
                                const Eigen::Vector3d& half_dims,
                                const Eigen::Matrix3d& rotation = Eigen::Matrix3d::Identity(),
                                const Eigen::Vector3d& inner_half_dims
                                    = Eigen::Vector3d::Zero());

    /// Configure a Marinelli beaker source geometry.
    ///
    /// All distances are relative to the detector geometry front face
    /// (z_det_min = geometry_.outer_z_extent().first), which includes
    /// dead layers and attenuator thicknesses — NOT z=0 (crystal front face).
    /// The detector must be fully configured before calling this method.
    ///
    /// @param well_inner_radius  Inner radius of re-entrant well (cm); must clear detector
    /// @param well_depth         How far well extends past detector geom front (cm)
    /// @param outer_radius       Outer radius of beaker body (cm)
    /// @param fill_height        Height of sample disk in front of beaker opening (cm)
    /// @param endcap_to_beaker   Air gap from detector geom front to beaker opening (cm)
    /// @param sample_material    Material filling the sample volume (e.g., water)
    /// @param beaker_material    Container wall material (e.g., polyethylene)
    /// @param beaker_thickness   Container wall thickness (cm)
    void set_marinelli_beaker(double well_inner_radius, double well_depth,
                              double outer_radius, double fill_height,
                              double endcap_to_beaker,
                              const Material* sample_material,
                              const Material* beaker_material, double beaker_thickness);

    // --- Source Material & Shielding ---

    /// Set the material that fills the source volume (for self-attenuation).
    /// Only meaningful for extended sources (cylindrical/rectangular/spherical).
    ///
    /// **Trace and self-attenuating sources (mapping to InterSpec).** CeeLo
    /// computes per-photon *efficiency*; the absolute activity scale is a
    /// downstream normalization (InterSpec multiplies efficiency by Bq×volume),
    /// so the only things that affect efficiency are the emission's *spatial
    /// distribution* and the *self-attenuating medium*:
    ///   - **Self-attenuating source** (the source IS the active matrix, e.g. U/Th
    ///     metal): set the extended geometry + `set_source_material(matrix)`. The
    ///     volume self-absorbs. InterSpec's per-element mass-fraction is part of
    ///     the activity normalization and does not change efficiency.
    ///   - **Trace source** (activity at a concentration in a host that may also be
    ///     the shield, e.g. soil/water contamination): set the extended geometry +
    ///     `set_source_material(host)`; choose the emission distribution —
    ///     **uniform** (InterSpec TotalActivity / ActivityPerCm3 / ActivityPerGram)
    ///     or **exponential depth** (InterSpec ExponentialDistribution, Bq/m²
    ///     surface contamination) via `set_exponential_depth_distribution(λ)`,
    ///     where λ maps to InterSpec's relaxation length. Concentration *units*
    ///     (per cm³/g/m²) only scale the absolute rate, not the efficiency.
    void set_source_material(const Material* mat);

    /// Add a shielding layer around the source (innermost first).
    /// For point sources: spherical shells (path = thickness).
    /// For extended sources: shells matching the source shape, with the same
    /// thickness in every dimension.
    void add_source_shield(const Material* mat, double thickness);

    /// Add a cylindrical-source shielding layer with independent radial and
    /// end-cap thicknesses (cm). The end thickness applies to both end caps.
    /// One (but not both) may be zero, e.g. (t_radial, 0) = side wall only.
    /// Only valid for cylindrical sources; call after set_cylindrical_source().
    void add_source_shield(const Material* mat, double t_radial, double t_end);

    /// Add a rectangular-source shielding layer with independent x/y/z
    /// thicknesses (cm), applied on both +/- faces of each axis. One or two
    /// may be zero (but not all three). Only valid for rectangular sources;
    /// call after set_rectangular_source().
    void add_source_shield(const Material* mat, double t_x, double t_y, double t_z);

    // --- Source Depth Distribution ---

    /// Set exponential depth distribution for extended sources.
    /// Activity is distributed as A(d) = exp(-d/relaxation_length) where d is
    /// depth from the surface nearest the detector (d=0 at the detector-facing
    /// surface, d increases away from detector).
    /// Only meaningful for Cylindrical and Rectangular source types.
    /// @param relaxation_length  Relaxation length lambda (cm). Must be > 0.
    void set_exponential_depth_distribution(double relaxation_length);

    /// Reset to uniform depth distribution (the default).
    void set_uniform_depth_distribution();

    // --- Transport Configuration ---

    void set_max_compton_scatters(int n) { transport_config_.max_compton_scatters = n; }
    void set_min_energy_keV(double e) { transport_config_.min_energy_keV = e; }
    void enable_rayleigh(bool on) { transport_config_.enable_rayleigh = on; }
    void enable_doppler_broadening(bool on) { transport_config_.enable_doppler_broadening = on; }
    void enable_electron_csda(bool on) { transport_config_.enable_electron_csda = on; }
    void enable_cone_sampling(bool on) { use_cone_sampling_ = on; }
    void enable_fep_only_mode(bool on) { fep_only_mode_ = on; }
    void set_disable_fep_early_kill(bool on) { transport_config_.disable_fep_early_kill = on; }

    /// Enable optional source electron transport. Compton recoil electrons
    /// from source material are tracked through source geometry to the
    /// detector via CSDA. Only effective in full mode (not FEP-only).
    /// Performance-gated: electrons are filtered by energy, geometric
    /// intersection with detector, and CSDA range budget. Default: disabled.
    void enable_source_electron_transport(bool on = true);
    void set_source_electron_threshold(double keV) { source_geometry_.set_source_electron_threshold(keV); }
    void set_source_electron_geom_check(bool on) { source_geometry_.set_source_electron_geom_check(on); }
    void set_disable_moliere(bool on) { transport_config_.disable_moliere = on; }
    void set_disable_brems(bool on) { transport_config_.disable_brems = on; }
    /// Disable bremsstrahlung from electrons stopping in source
    /// material/shielding (on by default; off switch for A/B quantification).
    void set_source_brems(bool on) { source_geometry_.set_source_brems(on); }
    void set_air_attenuation(AirAttenuation mode) { air_attenuation_ = mode; }

    // --- Event Biasing (variance reduction) ---

    /// Set event-biasing configuration explicitly. Overrides the
    /// auto-enabled defaults; pass a default-constructed BiasingConfig to
    /// force fully analog sampling (for A/B consistency tests).
    void set_biasing(const BiasingConfig& config) {
        biasing_ = config;
        biasing_manual_ = true;
    }
    const BiasingConfig& biasing() const { return biasing_; }

    /// The biasing configuration a compute() call would actually use:
    /// the manual config if set_biasing() was called, otherwise the
    /// auto-enable policy (measured on the validated benchmark configs and
    /// an out-of-sample stress set):
    /// - force_detector_interaction when the central-ray optical depth
    ///   through the detector is < 0.7 (high pass-through probability),
    ///   or < 4.0 when only total-efficiency precision is targeted.
    /// - two_stream for non-Marinelli full-mode source-effect runs whose
    ///   detector cone from the source center subtends omega/4pi < 0.15:
    ///   f = 0.5 when only FEP precision is targeted, else f = 0.25 with
    ///   compton_cone_fraction = 0.3 (total/spectrum requested).
    /// - mixture_cone_alpha is no longer auto-enabled (superseded by
    ///   two_stream); it remains available via set_biasing().
    BiasingConfig compute_effective_biasing(const SimulationConfig& config) const;

    // --- Position Biasing ---

    /// Enable position biasing with manually specified parameters.
    /// Overrides auto-computed parameters.
    void set_position_bias(const PositionBiasConfig& config);

    /// Enable position biasing with auto-computed parameters.
    /// Parameters are computed once at the start of each compute() call
    /// based on the source geometry, material, and photon energy.
    void enable_position_bias();

    /// Disable position biasing (return to uniform sampling).
    void disable_position_bias();

    // --- Simulation ---

    /// Run simulation with full configuration (flexible termination, progress callback).
    EfficiencyResult compute(const SimulationConfig& config);

    /// Run simulation for a single energy (convenience wrapper).
    EfficiencyResult compute(double energy_keV,
                             uint64_t num_events,
                             unsigned num_threads = 0,
                             const std::vector<float>& energy_bin_edges = {});

    /// Run simulation for multiple energies (batch mode).
    std::vector<EfficiencyResult> compute_batch(
                             const std::vector<double>& energies_keV,
                             uint64_t num_events_per_energy,
                             unsigned num_threads = 0,
                             const std::vector<float>& energy_bin_edges = {});

    // --- Gold-standard correlated-cascade MC ---

    /// Evaluate true-coincidence summing per peak via correlated-cascade MC.
    /// For each peak, returns the single-gamma FEP efficiency, the
    /// summing-corrected FEP efficiency, and their ratio (the summing factor),
    /// all with statistical uncertainties. Point and extended sources, with
    /// optional source material/shielding, are supported. Opt-in; does not
    /// affect compute()/simulate_thread().
    CascadeResult compute_cascade(const CascadeConfig& config) const;

    // --- Cascade Correction Helper ---

    /// Compute summing-out correction only.
    /// C_out = Π_j [1 - f_j × ε_total(E_j)]
    /// summing_in_term = 0, net_correction = summing_out_factor.
    /// @param primary_fep_uncertainty  1σ on primary_fep (propagated into
    ///                                 corrected_fep_uncertainty).
    /// @param match_tolerance_keV      Energies are matched to cache keys within
    ///                                 this absolute tolerance; any unmatched
    ///                                 energy is reported in the result.
    static CascadeCorrectionResult cascade_correction(
                             double primary_fep,
                             double primary_total,
                             const std::vector<CoincidentGamma>& coincident,
                             const std::map<double, EfficiencyResult>& efficiency_cache,
                             double primary_fep_uncertainty = 0.0,
                             double match_tolerance_keV = kCascadeEnergyMatchTolKeV);

    /// Compute summing-out AND summing-in corrections.
    /// Summing-in:  C_in = Σ_{pairs} f_ab × ε_FEP(a) × ε_FEP(b) / ε_FEP(primary)
    /// Net:         C_net = C_out + C_in
    /// Corrected:   ε_corrected = ε_FEP(primary) × C_net
    ///
    /// @param summing_in_pairs  Pairs of coincident gammas whose energies sum to
    ///                          the primary gamma energy.  Their FEP efficiencies
    ///                          are looked up in efficiency_cache.
    /// @param primary_fep_uncertainty  1σ on primary_fep (propagated into
    ///                                 corrected_fep_uncertainty).
    /// @param match_tolerance_keV      Cache-key match tolerance (keV); unmatched
    ///                                 energies are reported in the result.
    static CascadeCorrectionResult cascade_correction(
                             double primary_fep,
                             double primary_total,
                             const std::vector<CoincidentGamma>& coincident,
                             const std::map<double, EfficiencyResult>& efficiency_cache,
                             const std::vector<SummingInPair>& summing_in_pairs,
                             double primary_fep_uncertainty = 0.0,
                             double match_tolerance_keV = kCascadeEnergyMatchTolKeV);

    /// Conditional partner-emission probabilities P(member m gamma emitted |
    /// primary gamma emitted) derived from the daughter LEVEL SCHEME via the
    /// exact visit/reach dynamic program (same DP as `LevelWalk` in
    /// the cascade design study). Correctly gives p=0 for transitions
    /// that are mutually EXCLUSIVE with the primary (alternative de-excitations of
    /// a shared level), unlike the pairwise marginal fallback that spuriously sums
    /// them (Co-57 136/122, Ba-133, Eu-152...). Members with no matching
    /// level-scheme transition (x-rays / annihilation, emitted independently of
    /// the gamma path) keep their pairwise/marginal probability. Returns an EMPTY
    /// vector when `dc.level_scheme` is invalid or the primary is not a
    /// level-scheme gamma, signalling the caller to fall back to the pure pairwise
    /// construction. Public so it can be unit-tested directly (no MC / no state).
    static std::vector<double> cascade_level_pmate(const DecayCascade& dc,
                                                   std::size_t primary_index);

    /// One K/L-shell vacancy that fires coincident with a conditioned emission
    /// (the primary gamma, or a sum-feeding pair), with its per-primary-emission
    /// probability. Emitted as a single fluorescence line (or Auger) by the
    /// Conditional estimator, mirroring FullRealization's vacancy x-ray model.
    struct CascadeVacancyDraw {
        double prob = 0.0;     ///< P(this vacancy fires | the conditioned emission)
        bool is_L = false;     ///< false => K-shell vacancy; true => L-shell
        int l_subshell = -1;   ///< L subshell 0/1/2 (is_L); -1 = K or unresolved
        /// Transition energy (keV) of the internal conversion that made this
        /// vacancy, for depositing the IC electron KE (= trans_keV - binding).
        /// 0 for an electron-capture vacancy (no conversion electron).
        double trans_keV = 0.0;
        /// Alternatives sharing a non-negative group are selected with one
        /// categorical draw (none is the residual).  This enforces one EC shell
        /// per capture and one conversion shell per nuclear transition.
        int exclusive_group = -1;
        /// Gamma alternative of this conversion group, or -1 for EC/E0.  If it
        /// fired, no vacancy in this group may fire; otherwise `prob` is divided
        /// by the complementary gamma probability before the categorical draw.
        int gamma_member = -1;
        /// P(this vacancy | conditioned emission, selected weak outcome k).
        /// Empty for legacy cascades and EC vacancies (the latter are sampled
        /// directly from the selected WeakOutcome).  Keeping this conditional
        /// law prevents an EC-fed path from borrowing an IC transition that is
        /// reachable only from a competing beta+ outcome, and vice versa.
        std::vector<double> weak_outcome_prob;
        /// False for a shell-unresolved E0 conversion: transport the approximate
        /// electron but do not invent a K/L vacancy or fluorescence x-ray.
        bool produces_vacancy = true;
    };

    /// One sum-peak-feeding channel for a peak in the Conditional estimator: a
    /// pair of daughter gammas (a,b) whose energies add up into the peak window
    /// even though NEITHER is the peak's own gamma (e.g. Co-57 122.06+14.41 ->
    /// 136.47). `weight` is the analytic per-primary-emission expectation of the
    /// pair being emitted (from the exact level-scheme joint DP), with the triple
    /// overlap the primary conditional stream already samples subtracted out.
    struct CascadeSumPairChannel {
        double e_a_keV = 0.0, e_b_keV = 0.0;
        int member_a = -1, member_b = -1;  ///< member indices (W(theta) links)
        std::size_t branch = 0;            ///< index into the cascades vector
        double weight = 0.0;               ///< w_k per primary emission
        std::vector<double> partner_prob;  ///< P(member m | pair a&b), that branch
        /// Per-selected-weak-outcome version of partner_prob.  Rows correspond
        /// to WeakOutcomeLaw::outcomes; columns to DecayCascade::members.
        std::vector<std::vector<double>> weak_outcome_partner_prob;
        /// Posterior categorical weak outcome P(o | pair a&b emitted). Empty for
        /// legacy/host-created cascades without an explicit outcome law.
        std::vector<double> weak_outcome_posterior;
        std::vector<CascadeVacancyDraw> vacancies;  ///< vacancies joint with the pair
    };

    /// Vacancy x-ray draws that fire coincident with the primary gamma being
    /// emitted, from the daughter level scheme (EC feed + IC of the passed
    /// transitions), for the Conditional estimator's primary stream. Empty when
    /// `dc.level_scheme` is invalid or the primary is not a level-scheme gamma.
    /// Public for unit testing.
    static std::vector<CascadeVacancyDraw> cascade_level_vacancies(
        const DecayCascade& dc, std::size_t primary_index);

    /// Enumerate the sum-peak-feeding pair channels for `peak_keV` across all
    /// level-scheme-valid branches (the primary gamma is
    /// `cascades[primary_branch].members[primary_index]`). Each returned channel
    /// carries its analytic per-primary-emission weight, pair-conditioned partner
    /// probabilities, and coincident vacancies. Returns empty when the primary
    /// branch has no valid level scheme. Invalid-branch flat-member pairs are
    /// deliberately not synthesized from inconsistent pairwise metadata; the
    /// public result's `summing_model_complete` flag reports when that omission
    /// can affect a requested peak. Public for unit testing.
    static std::vector<CascadeSumPairChannel> cascade_sum_pair_channels(
        const std::vector<DecayCascade>& cascades, std::size_t primary_branch,
        std::size_t primary_index, double peak_keV, double tol_keV);

    // --- GEANT4 Export ---

    /// Export the detector geometry as a GDML file for GEANT4 validation.
    /// The active crystal volume is named "active_crystal".
    /// @param vacuum_world  If true, world volume uses vacuum instead of air.
    /// @see tools/geant4_validation/ for the matching GEANT4 harness.
    void export_geant4_gdml(const std::string& filename,
                             bool vacuum_world = false) const;

    /// Export a GEANT4 GPS macro file for a point-source simulation.
    /// Uses the currently configured source position (center for extended sources).
    void export_geant4_macro(const std::string& filename,
                             double energy_keV,
                             uint64_t num_events) const;

    // --- Access to geometry and source (for testing/tools) ---
    const Geometry& geometry() const { return geometry_; }
    /// Sample one unbiased source position (public for unit tests, e.g.
    /// hollow-shell rejection-sampling coverage). Consumes only `rng`.
    Eigen::Vector3d sample_source_position_for_test(std::mt19937_64& rng) const {
        return sample_source_position(rng);
    }
    Eigen::Vector3d source_position() const { return source_position_; }
    SourceType source_type() const { return source_type_; }
    const MarinelliBeakerConfig& marinelli_config() const { return marinelli_cfg_; }
    const SourceGeometry& source_geometry() const { return source_geometry_; }

    // Internal tally structure (public for SimulationState access)
    struct ThreadTally {
        uint64_t num_events = 0;
        uint64_t num_fep = 0;      ///< Full-energy peak detections
        uint64_t num_any = 0;      ///< Any-energy detections (ε_total)
        double sum_weights = 0.0;
        double sum_fep_weights = 0.0;
        double sum_any_weights = 0.0;
        double sum_fep_w_sq = 0.0;   ///< Sum of total_weight² for FEP events (IS variance)
        double sum_any_w_sq = 0.0;   ///< Sum of total_weight² for any-deposit events (IS variance)
        uint64_t num_forced_absorption = 0; ///< Primary photons force-absorbed at the
                                            ///< max_compton_scatters cap (diagnostic; see A1)
        std::vector<double> pulse_height_counts; ///< Weighted counts (Sum w) per bin
        std::vector<double> pulse_height_counts_sq; ///< Sum w^2 per bin (per-bin IS variance)

        // Temporary Marinelli re-entry diagnostics
        uint64_t dbg_initial_hit = 0;       ///< Initial crystal transports
        uint64_t dbg_miss_bounce_hit = 0;   ///< Miss-bounce recoveries that hit crystal
        uint64_t dbg_reentry_hit = 0;       ///< Post-crystal re-entries that hit crystal
        uint64_t dbg_secondary_hit = 0;     ///< Escaped secondaries that hit crystal
        uint64_t dbg_pp_secondary = 0;      ///< PP annihilation gammas processed
        double dbg_dep_initial = 0.0;       ///< Energy deposited from initial transport
        double dbg_dep_reentry = 0.0;       ///< Energy deposited from re-entries
        double dbg_dep_secondary = 0.0;     ///< Energy deposited from escaped secondaries
        double dbg_dep_pp = 0.0;            ///< Energy deposited from PP 511 keV gammas
        uint64_t dbg_pp_only_any = 0;       ///< Events where any-deposit came only from PP gammas
        double dbg_pp_only_any_w = 0.0;     ///< Summed event weight of those events

        // Re-entry pipeline funnel counters
        uint64_t dbg_n_escaped = 0;         ///< Photons that escape crystal
        uint64_t dbg_n_can_reenter = 0;     ///< Of escaped, compute_marinelli_reentry true
        uint64_t dbg_n_wall_pass = 0;       ///< Of reenterable, wall_transport passed
        uint64_t dbg_n_wall_scatter = 0;    ///< Of reenterable, wall_transport scattered
        uint64_t dbg_n_wall_absorb = 0;     ///< Of reenterable, wall_transport absorbed
        uint64_t dbg_n_water_survived = 0;  ///< Of wall-passed, source_photon survived
        uint64_t dbg_n_trace_hit = 0;       ///< Of water-survived, trace_ray hits crystal
        uint64_t dbg_n_trace_miss = 0;      ///< Of water-survived, trace_ray misses

        /// Source-effect variance decomposition by primary class (u/s).
        SourceEffectDiag dbg_src;
    };

private:
    Geometry geometry_;
    TransportConfig transport_config_;

    // Source configuration
    SourceType source_type_ = SourceType::Point;
    Eigen::Vector3d source_position_{0.0, 0.0, -10.0};
    CylindricalSourceConfig cyl_src_;
    RectangularSourceConfig rect_src_;
    SphericalSourceConfig sph_src_;
    MarinelliBeakerConfig marinelli_cfg_;

    // Options
    bool use_cone_sampling_ = true;
    bool fep_only_mode_ = false;
    bool source_electron_transport_ = false;
    AirAttenuation air_attenuation_ = AirAttenuation::None;

    // Event biasing (variance reduction)
    BiasingConfig biasing_;
    bool biasing_manual_ = false;  ///< true if user called set_biasing()

    /// Effective biasing for the current compute() call: the manual config
    /// if set_biasing() was called, otherwise auto-enabled per geometry,
    /// energy, and termination targets (see compute_effective_biasing()).
    BiasingConfig active_biasing_;

    /// Total optical depth through the detector geometry along the ray from
    /// the source center toward the detector front-face center, at the given
    /// energy. Returns +inf for Marinelli sources (no meaningful central ray).
    double central_ray_optical_depth(double energy_keV) const;

    // Depth distribution
    DepthDistribution depth_distribution_ = DepthDistribution::Uniform;
    double relaxation_length_ = 1.0;  ///< cm, only used when Exponential

    // Position biasing
    PositionBiasConfig position_bias_;
    bool position_bias_enabled_ = false;
    bool position_bias_manual_ = false;  ///< true if user called set_position_bias()

    // Source material and shielding
    SourceGeometry source_geometry_;

    ThreadTally simulate_thread(double energy_keV,
                                uint64_t num_events,
                                uint64_t seed,
                                const std::vector<float>& energy_bin_edges) const;

    // --- Correlated-cascade MC internals (compute_cascade) ---

    /// Accumulators for one peak's correlated-cascade MC. The per-history
    /// "with-summing" score X = I_sum + Σ_k w_k·P_k is a real number (the sum-fed
    /// pair channels carry analytic weights w_k), so it is accumulated in doubles.
    /// When there are no channels and no vacancies, X ∈ {0,1} and these reduce to
    /// the old integer n_sum / n_both counters exactly. The per-history primary
    /// cone-bias weight is constant, so it is applied once in compute_cascade().
    struct CascadePeakTally {
        uint64_t n = 0;        ///< histories
        uint64_t n_nosum = 0;  ///< primary-alone deposit in the peak window (I_nosum)
        double sum_x = 0.0;    ///< Σ X (with-summing score, incl. sum-fed channels)
        double sum_xx = 0.0;   ///< Σ X² (for the with-summing variance)
        double sum_xd = 0.0;   ///< Σ X·I_nosum (for the correlated ratio variance)
    };

    /// One thread of the per-peak correlated-cascade MC. Conditions on the
    /// primary gamma (cascades[primary_branch].members[primary_index]) being
    /// emitted, samples its coincident partners and vacancy x-rays, adds the
    /// analytic sum-peak-fed pair channels, transports all members, and tallies
    /// the windowed primary-alone and with-summing deposits. `progress_flush`,
    /// when set, is called periodically (time-throttled) with the thread-local
    /// tally so far; it never touches the RNG, so results are identical either way.
    CascadePeakTally cascade_peak_thread(
        const std::vector<DecayCascade>& cascades, std::size_t primary_branch,
        std::size_t primary_index, double peak_keV, double tol_keV,
        double cos_theta_max, bool coned,
        const std::vector<CascadeVacancyDraw>& prim_vacancies,
        const std::vector<CascadeSumPairChannel>& channels,
        uint64_t num_events, uint64_t seed, bool ic_electrons,
        const std::function<void(const CascadePeakTally&)>& progress_flush = {}) const;

    /// Convert an aggregated conditional-method tally into a PeakCascadeResult
    /// (binomial uncertainties + correlated-ratio summing-factor variance).
    /// Shared by the final result build and the progress snapshot.
    static PeakCascadeResult conditional_peak_result(const CascadePeakTally& agg,
                                                     double energy_keV, bool found,
                                                     double prim_w);

    /// Emit the fluorescence x-ray for one K/L vacancy from the vertex and return
    /// its scoring-volume deposit (0 for Auger / no data). Mirrors the K
    /// relaxation (incl. the Kalpha-induced L vacancy) and the L Coster-Kronig
    /// loop in cascade_full_thread -- KEEP IN SYNC with that reference path.
    double emit_vacancy_xray_deposit(int daughter_Z, bool is_L, int l_subshell,
                                     const Eigen::Vector3d& vertex,
                                     std::vector<PathSegment>& seg_buffer,
                                     std::uniform_real_distribution<double>& u,
                                     std::mt19937_64& rng,
                                     bool ic_electrons = false) const;

    /// Deposit an internal-conversion (or coincident K-Auger, binding_keV=0)
    /// electron of kinetic energy (trans_keV - binding_keV) born isotropically at
    /// the decay `vertex`, transported through the source geometry escape walk +
    /// the (charged-particle-degraded) air gap, then the detector CSDA. Returns
    /// the crystal deposit (keV, incl. bremsstrahlung). Contained-in-source,
    /// air-absorbed, or off-solid-angle electrons contribute 0 (so distant
    /// external sources are unaffected). RNG-NEUTRAL when `enabled` is false:
    /// returns 0 before any draw, so the legacy photon-only path is bit-identical.
    /// Shared by both cascade estimators; see the block comment at the call sites.
    double emit_ic_electron(double trans_keV, double binding_keV, bool enabled,
                            const Eigen::Vector3d& vertex,
                            std::vector<PathSegment>& seg_buffer,
                            std::uniform_real_distribution<double>& u,
                            std::mt19937_64& rng) const;

    /// Emit a cascade member from the source vertex and return its total
    /// scoring-volume deposit. When the source has material/shield effects the
    /// member is first transported through the source geometry (self-attenuation
    /// + Compton/Rayleigh scatter), and the survivor plus any source-generated
    /// secondary photons (511s from pair production, bremsstrahlung) are summed.
    /// Otherwise it goes straight to the detector. (Source-escaped recoil
    /// electrons are not yet included.)
    double transport_cascade_member(const Eigen::Vector3d& vertex,
                                    const Eigen::Vector3d& direction,
                                    double energy_keV,
                                    std::vector<PathSegment>& seg_buffer,
                                    std::mt19937_64& rng) const;

    /// Trace a photon from `start` (already past any source geometry) into the
    /// detector and return its scoring-volume deposit (0 if it misses).
    double transport_into_detector(const Eigen::Vector3d& start,
                                   const Eigen::Vector3d& direction,
                                   double energy_keV,
                                   std::vector<PathSegment>& seg_buffer,
                                   std::mt19937_64& rng) const;

    /// Per-peak target. `branch`/`member` retain the dominant primary used by the
    /// single-primary Conditional estimator; `primary_matches` contains every
    /// positive-absolute-yield gamma identity in the window for exact
    /// FullRealization normalization.
    struct CascadePeakTarget {
        struct PrimaryIdentity {
            std::size_t branch = 0;
            std::size_t member = 0;
        };
        std::size_t branch = 0;
        std::size_t member = 0;
        std::vector<PrimaryIdentity> primary_matches;
        double energy_keV = 0.0;
        double tol_keV = 1.5;
        bool found = false;
        bool summing_model_complete = true;
    };

    /// Accumulator for the full-realization estimator (one thread).
    struct CascadeFullTally {
        uint64_t n = 0;                       ///< histories
        std::vector<uint64_t> n_emit;         ///< per peak: peak gamma emitted
        std::vector<uint64_t> n_nosum;        ///< per peak: peak gamma alone in window
        std::vector<uint64_t> n_sum;          ///< per peak: summed deposit in window
        std::vector<double>   spectrum;       ///< summed-deposit histogram (counts)
        double branch_weight_sum = 0.0;       ///< Σ branch_weight over all branches (R)
    };

    /// One thread of the full-realization correlated-cascade MC.
    /// `progress_flush` as in cascade_peak_thread (RNG-neutral).
    CascadeFullTally cascade_full_thread(
        const std::vector<DecayCascade>& cascades,
        const std::vector<double>& branch_cdf, double branch_weight_sum,
        const std::vector<CascadePeakTarget>& targets,
        const std::vector<float>& spectrum_edges,
        uint64_t num_events, uint64_t seed, bool ic_electrons,
        const std::function<void(const CascadeFullTally&)>& progress_flush = {}) const;

    /// Convert one peak's aggregated full-realization counters into a
    /// PeakCascadeResult. Shared by the final result build and the progress
    /// snapshot.
    static PeakCascadeResult full_peak_result(uint64_t n_emit, uint64_t n_nosum,
                                              uint64_t n_sum,
                                              double energy_keV, bool found);

    /// Map each requested peak to all positive-yield gamma identities and the
    /// dominant identity used by single-primary Conditional sampling. A usable
    /// level scheme is ranked by its absolute gamma yield (branch_weight × exact
    /// weak-outcome reach × p_gamma); metadata-free cascades fall back to
    /// branch_weight × member intensity. Shared by both CascadeMethods.
    std::vector<CascadePeakTarget> cascade_locate_targets(
        const CascadeConfig& config) const;

    /// Sample a source position.
    /// For point sources, always returns source_position_.
    /// For extended sources, samples from the source volume according to
    /// the configured depth distribution (uniform or exponential).
    Eigen::Vector3d sample_source_position(std::mt19937_64& rng) const;

    /// Sample a depth from the truncated exponential distribution on [0, D].
    double sample_exponential_depth(double total_depth, double U) const;

    /// Sample a biased source position using importance sampling.
    /// Returns {position, weight} where weight = p_physical / p_biased.
    /// Weight is capped at 100.0 to prevent pathological outliers.
    std::pair<Eigen::Vector3d, double> sample_biased_position(
        std::mt19937_64& rng,
        const PositionBiasConfig& bias) const;

    /// Auto-compute position bias parameters from source geometry and energy.
    PositionBiasConfig compute_auto_bias_params(double energy_keV) const;

    // Cone sampling helpers
    double compute_raw_cone_half_angle(const Eigen::Vector3d& from_pos) const;
    double compute_cone_half_angle(const Eigen::Vector3d& from_pos) const;
    double compute_worst_case_cone() const;
    Eigen::Vector3d sample_cone_direction(const Eigen::Vector3d& from_pos,
                                          double cos_theta_max,
                                          std::mt19937_64& rng) const;

    /// Sample a direction by targeting a random point on the detector's outer
    /// surface. Returns {direction, weight} where weight = A * cos_inc / (4pi * d^2).
    /// Returns weight = 0 if the sampled surface faces away from the source.
    std::pair<Eigen::Vector3d, double> sample_surface_direction(
        const Eigen::Vector3d& src_pos, std::mt19937_64& rng) const;
};

} // namespace ceelo
