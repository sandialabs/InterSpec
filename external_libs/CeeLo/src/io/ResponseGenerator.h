#ifndef CEELO_IO_RESPONSE_GENERATOR_H
#define CEELO_IO_RESPONSE_GENERATOR_H
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

/// @file ResponseGenerator.h
/// @brief MC generation driver for the parameterized DetectorResponse.
///
/// Turns a GeometryDescriptor into a filled DetectorResponse by running
/// precision-targeted CeeLo MC at generator-chosen nodes (the July 2026
/// campaign recipe -- see DetectorResponse.h and, in the dev repo,
/// the detector-response parameterization campaign):
///
///   1. Energy backbone: dense on-axis far-field scan (log grid + K-edge
///      flanks) -> greedy node selection with segment endpoints forced ->
///      locally-smoothed node values (fit through denser data, not raw
///      noisy nodes).
///   2. Angular shape: cos-theta scans at a handful of energies (+ phi for
///      boxes) -> greedy shared cos-theta nodes -> the eta tensor is
///      backbone x lnE-interpolated shape (the S10-validated assembly).
///   3. Near field (general/contact profiles): direct (cos-theta, d) tensor
///      scan below ~4a -> tabulated ln N on that grid per scan energy,
///      PCHIP-interpolated -> measured 1% breakpoints.
///   4. eps_tot tier auto-choice: score kernel-exact and b(E) tiers against
///      the free eps_tot at every scanned node; ship the most compact tier
///      whose 95th-percentile error is <= 1%, else the eta_tot table.
///   5. Grounding (optional, separate call): Level-1 GLS fit of ln k on a
///      hat basis with per-point stat sigma + per-source certificate
///      covariance blocks (spec Eqs. 7a-7d).
///
/// All MC is deterministically seeded from options.base_seed; regeneration
/// with the same options is reproducible (in distribution across threads;
/// bit-exact with num_threads = 1).

#include "io/DetectorResponse.h"

#include <atomic>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace ceelo {

class Material;
class EfficiencyCalculator;

/// Per-MC-node cost/outcome record (generation accounting; one per run_node).
struct NodeStat {
    uint32_t stage = 0;        ///< scan family: 1 backbone, 2 angular, 3 near
    double energy_keV = 0.0;
    double cos_theta = 1.0;    ///< of the source position vs detector axis
    double d_cm = 0.0;         ///< source distance from the origin
    uint64_t events = 0;
    double cpu_s = 0.0;        ///< summed worker-thread CPU seconds
    double wall_s = 0.0;
    uint8_t stop = 0;          ///< numeric ceelo::StopReason value
    double fep_rel_prec = 0.0; ///< achieved sigma/eps (0 when eps == 0)
};

/// Aggregated per-generation MC cost statistics (see GenerationOptions::stats_out).
struct GenerationStats {
    std::vector<NodeStat> nodes;
    double cpu_s_by_stage[8] = {};  ///< indexed by NodeStat::stage (clamped)
    double total_cpu_s = 0.0;       ///< sum of node CPU seconds
    double total_wall_s = 0.0;      ///< sum of node wall seconds (serial nodes)
    uint64_t total_events = 0;

    void add(const NodeStat& ns) {
        nodes.push_back(ns);
        cpu_s_by_stage[ns.stage < 8 ? ns.stage : 7] += ns.cpu_s;
        total_cpu_s += ns.cpu_s;
        total_wall_s += ns.wall_s;
        total_events += ns.events;
    }
};

/// Per-node MC precision policy (the D0 speedup dial). Uniform reproduces the
/// historical single-target behaviour bit-for-bit; RelaxMild installs the
/// energy-graded map measured in the D0 policy memo (~2-2.5x cheaper high-E
/// nodes at a bounded p95 cost).
enum class PrecisionProfile : uint8_t {
    Uniform,    ///< node_fep_precision everywhere (default; unchanged)
    RelaxMild   ///< graded target precision by (stage, energy) -- see below
};

/// The `relax_mild` graded target-precision map (D0 policy memo): low-E rows
/// stay tight, stepping toward 3 MeV where the cap-limited nodes live. Applies
/// to STAGE 2 (angular) and STAGE 3 (near) rows; STAGE 1 backbone stays at
/// `base` (cheap + load-bearing). A pure function of (stage, E) -- deterministic.
std::function<double(uint32_t stage, double E_keV)>
relax_mild_precision_map(double base = 0.003);

/// Structured probe families for the closed-loop generator (D-b). Each value
/// doubles as the certificate Row::tag it produces; Random is the model-form
/// canary. See plan_structured_probes / structured_probe_bank.
enum class ProbeFamily : uint8_t {
    Random    = 0,  ///< quasi-random Halton probe (model-form canary)
    EnergyGap = 1,  ///< geom-mean E midpoint of adjacent backbone nodes, far, on-axis
    CtGap     = 2,  ///< cos-theta midpoint at a shape-E node, far
    ShapeEGap = 3,  ///< E between adjacent shape-E rows, mid-ct, far + d~0.6a
    NearCell  = 4   ///< near-field d-gap / ct-gap midpoint, near field
};

/// Bitmask over ProbeFamily (1 << family). kAllStructuredFamilies selects every
/// structured family (excludes Random); kAllProbeFamilies adds the canary.
using ProbeFamilyMask = uint32_t;
constexpr ProbeFamilyMask probe_family_bit(ProbeFamily f) {
    return ProbeFamilyMask{1} << static_cast<uint32_t>(f);
}
constexpr ProbeFamilyMask kAllStructuredFamilies =
    probe_family_bit(ProbeFamily::EnergyGap) | probe_family_bit(ProbeFamily::CtGap) |
    probe_family_bit(ProbeFamily::ShapeEGap) | probe_family_bit(ProbeFamily::NearCell);
constexpr ProbeFamilyMask kAllProbeFamilies =
    kAllStructuredFamilies | probe_family_bit(ProbeFamily::Random);

/// Closed-loop initial grid selector (D-b). Full = the D0 full grid + relax_mild
/// (safe default, ~2-2.5x, robust by construction). Coarse = the memo's naive
/// 5 shape-E x 6 ct start (fewer nodes -> higher initial p95) so the refinement
/// loop demonstrably lifts it back under the gate.
enum class InitialGrid : uint8_t { Full, Coarse };

/// Options controlling the MC spend and grid budgets.
struct GenerationOptions {
    ResponseProfile profile = ResponseProfile::General;

    /// Per-node MC FEP precision target. 0.01 "fast", 0.003 default,
    /// 0.001 "thorough" (a user-visible choice in InterSpec).
    double node_fep_precision = 0.003;

    /// Graded per-node precision. Uniform (default) => node_fep_precision
    /// everywhere, bit-identical to history. RelaxMild => generate() installs
    /// relax_mild_precision_map(node_fep_precision) into node_precision.
    PrecisionProfile precision_profile = PrecisionProfile::Uniform;

    /// Optional per-node target FEP precision as f(stage, E_keV). Null (default)
    /// uses the scalar node_fep_precision everywhere -- the exact historical
    /// path. Set explicitly to override, or leave null and pick a
    /// precision_profile. Consulted ONLY for generation nodes (stages 1-3);
    /// probe/certificate banks always run at their own uniform precision.
    std::function<double(uint32_t stage, double E_keV)> node_precision;
    uint64_t max_events_per_node = 20000000;
    /// Per-node CPU-time cap, summed across worker threads. A CPU (not wall)
    /// cap keeps the node budget invariant to machine load and to concurrent
    /// jobs; 80 CPU-s reproduces the old 8 wall-s cap at ~10 worker threads,
    /// so the same node population stays cap-limited. An internal wall
    /// backstop (10x this / thread count) still catches a wedged run.
    double max_cpu_seconds_per_node = 80.0;
    /// Minimum MC events per node before precision-based termination.
    uint64_t min_events_per_node = 20000;

    double e_min_keV = 35.0;
    double e_max_keV = 3000.0;

    int n_energy_nodes = 24;        ///< greedy node CAP (noise/tol-limited, see below)
    int n_cos_theta_nodes = 12;     ///< greedy node CAP (noise/tol-limited)
    int n_energy_scan = 40;         ///< dense backbone scan size
    int n_cos_theta_scan = 14;      ///< dense angular scan size
    int n_shape_energies = 9;       ///< energies carrying the angular/near scans
                                    ///< (log-spaced; denser here resolves the
                                    ///< angular shape's fast low-E variation)
    int n_phi_nodes = 4;            ///< boxes only (quadrant/octant)

    /// Node placement is noise-aware and tolerance-driven: the greedy adds a node
    /// only when its ln-eta interpolation error exceeds BOTH node_interp_tol and
    /// node_noise_k * the local MC sigma (resolving sub-noise structure just fits
    /// scatter). So n_energy_nodes / n_cos_theta_nodes above are CAPS - the target
    /// accuracy decides how many are actually used, and the caps can be generous.
    double node_interp_tol = 0.003; ///< target ln-eta interpolation accuracy
    double node_noise_k = 1.0;      ///< don't resolve below this * local MC sigma

    double far_distance_a = 10.0;   ///< far-field scan distance, units of a
    double cos_theta_min = 0.02;    ///< grazing cutoff of the theta grid

    /// EFFTRAN transfer mode: MC ONLY the on-axis energy backbone (stage 1) and
    /// let the geometry kernel K carry all angular/distance transfer -- the
    /// cheap "anchor once, transfer everywhere" path (~24-40 nodes vs hundreds).
    /// n_anchor_angles == 1 => angle-flat eta (skips the angular scan entirely);
    /// > 1 => a few forced cos-theta MC anchors build a coarse eta(E,theta) for
    /// off-axis accuracy. Either way the near-field MC tensor (stage 3) is
    /// skipped (far-field product) and a SigmaTransferModel is attached so
    /// off-axis/near queries report honest, inflated sigma.
    bool transfer_mode = false;
    int  n_anchor_angles = 1;

    uint64_t base_seed = 1;         ///< deterministic node-seed base (never 0)
    unsigned num_threads = 0;       ///< per-node MC threads (0 = auto)
    int kernel_n_rays = 2048;       ///< evaluation-time quadrature rays

    /// Half-width (keV) of the full-energy-peak window the MC scores FEP with.
    /// Applied to the EfficiencyCalculator (which keeps the transport early
    /// kill in step) and recorded in ResponseProvenance.  See
    /// physics/FepWindow.h.
    double fep_window_keV = kDefaultFepWindowKeV;

    // ---- closed-loop refinement (D-b) --------------------------------------
    /// Opt-in closed-loop generation: coarse/full start -> structured probe ->
    /// attribute -> incremental refine -> re-certify. FALSE (default) runs the
    /// fixed-grid generate() path VERBATIM (goldens bit-identical). When true,
    /// generate() dispatches to the loop; the response carries an accuracy
    /// certificate recording convergence (never a silent miss).
    bool closed_loop = false;
    double cert_tol = 0.012;        ///< p95 FEP-error target (D0 memo gate)
    int max_refine_iters = 3;       ///< closed-loop iteration cap
    double refine_cpu_budget_frac = 0.5;  ///< refine MC budget, x initial-grid CPU
    int n_cert_probes = 48;         ///< final-certificate random probes
    InitialGrid initial_grid = InitialGrid::Full;

    /// Optional cost accounting: when non-null, cleared at generation start
    /// and filled with one NodeStat per MC node (caller owns; must outlive
    /// the generate() call). Adds no measurable overhead when null.
    GenerationStats* stats_out = nullptr;

    std::string detector_name;

    /// Progress: fraction complete [0,1] + human-readable stage. Called from
    /// the generation thread between MC nodes (lightweight; no UI work).
    std::function<void(double, const std::string&)> progress;
    /// Cooperative cancellation, polled between MC nodes. A cancelled run
    /// throws GenerationCancelled.
    std::shared_ptr<std::atomic<bool>> cancel;
};

/// Thrown when GenerationOptions::cancel is set during a run.
struct GenerationCancelled : public std::exception {
    const char* what() const noexcept override { return "generation cancelled"; }
};

/// One MC-computed probe/validation point (also the golden-fixture row).
struct ProbePoint {
    double energy_keV = 0.0;
    double d_cm = 0.0;             ///< from the crystal-face origin
    double cos_theta = 1.0;
    double phi_deg = 0.0;
    double eps_fep = 0.0, fep_unc = 0.0;   ///< MC truth, absolute 1-sigma
    double eps_tot = 0.0, tot_unc = 0.0;
    uint64_t seed = 0;
    uint8_t tag = 0;   ///< ProbeFamily of the probe (0 = Random); set by the
                       ///< structured banks, 0 for the plain/adversarial banks
};

class ResponseGenerator {
public:
    /// Full characterization: minutes of MC (profile/precision dependent).
    /// Throws GenerationCancelled on cancel; std::runtime_error on bad input.
    static std::shared_ptr<DetectorResponse> generate(
        const GeometryDescriptor& descriptor, const GenerationOptions& options);

    /// Rough node count for the given options (UI time estimates: each node
    /// is one precision-targeted MC run, typically 0.1-8 s).
    static int estimated_node_count(const GeometryDescriptor& descriptor,
                                    const GenerationOptions& options);

    /// Level-1 grounding (spec Eqs. 7a-7d): fit ln k(E) on a hat basis to
    /// raw measured points via GLS with V = diag(stat^2) + per-source
    /// certificate blocks; fills response.grounding (incl. covariance and a
    /// copy of the points). `points[i].model_eff` may be 0, in which case it
    /// is computed here from the (ungrounded) response at the point's own
    /// geometry. `curve_derived` marks points sampled from a fitted legacy
    /// curve rather than raw peak fits (lower quality; flagged in the block).
    static void ground_to_points(DetectorResponse& response,
                                 std::vector<GroundingPoint> points,
                                 bool curve_derived);

    /// Fresh quasi-random MC validation bank (never used in any fit; Halton
    /// offset start_index picks a disjoint set -- never split by parity).
    static std::vector<ProbePoint> probe_bank(
        const GeometryDescriptor& descriptor, const GenerationOptions& options,
        int n_points, int start_index,
        double d_min_cm = 0.5, double d_max_cm = 100.0);

    /// Adversarial (worst-case-interpolation) MC validation bank: places points
    /// at the MIDPOINTS of the widest energy and cos-theta node gaps of an
    /// already-built response, where the stored PCHIP interpolation error peaks
    /// by construction. Unlike the quasi-random probe_bank (which can miss a bad
    /// corner between nodes), this deliberately stresses every interpolation gap,
    /// so it reports the true worst-case error a random bank under-estimates.
    /// Energy midpoints are geometric means within a K-edge segment (never across
    /// an edge); each (E-mid, ct-mid) pair is sampled at a near and a far
    /// distance. Widest gaps first, truncated to n_points.
    ///
    /// D-b: now a thin wrapper over structured_probe_bank with the EnergyGap,
    /// CtGap and NearCell families (the same interpolation-gap stresses it used
    /// to place by hand). Kept for signature compatibility.
    static std::vector<ProbePoint> adversarial_probe_bank(
        const GeometryDescriptor& descriptor, const GenerationOptions& options,
        const DetectorResponse& response, int n_points, int start_index,
        double d_min_cm = 0.5, double d_max_cm = 100.0);

    /// One planned structured probe (family coordinates + tag + per-probe MC
    /// precision), WITHOUT MC -- the pure, deterministic output of the family
    /// planner (unit-testable: counts per family, midpoints within gaps, K-edge
    /// not bridged). structured_probe_bank runs MC at these coordinates.
    struct StructuredProbe {
        double energy_keV = 0.0, d_cm = 0.0, cos_theta = 1.0, phi_deg = 0.0;
        ProbeFamily family = ProbeFamily::Random;
        double node_precision = 0.0;  ///< MC target for this probe (EnergyGap
                                      ///< runs at node precision; else 0.008)
    };

    /// Pure family planner (no MC): places `families` probe coordinates on the
    /// node grid of an already-built `response` -- EnergyGap/CtGap/ShapeEGap/
    /// NearCell at interpolation-gap midpoints, plus Random Halton canaries.
    /// Up to `n_per_family` per selected family; `start_index` offsets the
    /// deterministic seed/Halton sequence. Deterministic for fixed inputs.
    static std::vector<StructuredProbe> plan_structured_probes(
        const GeometryDescriptor& descriptor, const GenerationOptions& options,
        const DetectorResponse& response, ProbeFamilyMask families,
        int n_per_family, int start_index,
        double d_min_cm = 0.5, double d_max_cm = 100.0);

    /// Structured (worst-case-interpolation) MC validation bank -- EXTENDS
    /// adversarial_probe_bank into tagged families. Runs MC at the coordinates
    /// from plan_structured_probes (each at its planned precision) and fills
    /// ProbePoint::tag with the family. The closed-loop generator uses this both
    /// to probe (per iteration) and to certify (the full family).
    static std::vector<ProbePoint> structured_probe_bank(
        const GeometryDescriptor& descriptor, const GenerationOptions& options,
        const DetectorResponse& response, ProbeFamilyMask families,
        int n_per_family, int start_index,
        double d_min_cm = 0.5, double d_max_cm = 100.0);

    /// Attach an honest accuracy certificate to `response`: run a fresh probe
    /// bank (quasi-random, Halton offset `seed_offset` -- disjoint from the
    /// generation scans' 5000 and never parity-split) at a fixed uniform
    /// precision, score each probe's |model/mc - 1| via eps_fep_at/eps_total_at,
    /// and fill response.certificate (summary percentiles over converged probes
    /// + per-probe rows). Metadata about the response, not part of its content
    /// (excluded from content_hash). Deterministic for fixed options.
    /// `cert_families` (D-b) adds tagged structured-family probes to the random
    /// bank so the certificate's p95/max reflect the worst-case interpolation
    /// gaps, not just random coverage. 0 (default) = random-only (the D-a
    /// behaviour, unchanged).
    static void certify(DetectorResponse& response,
                        const GeometryDescriptor& descriptor,
                        const GenerationOptions& options,
                        int n_probes = 48, int seed_offset = 7000,
                        ProbeFamilyMask cert_families = 0);

    /// Configure `calc`'s DETECTOR side (crystal, bore, dead layer, attenuator
    /// layers, collimator) from a stored descriptor — the same mapping the
    /// generator itself uses per MC node. The instantiated materials are
    /// appended to `owned_materials`, which the caller must keep alive for the
    /// lifetime of `calc` (the calculator holds raw pointers). Source-side
    /// setup (set_*_source / add_source_shield / air) remains the caller's.
    /// For hosts (e.g. InterSpec) running bespoke scenes against a stored
    /// DetectorResponse::descriptor.
    static void configure_calculator(
        EfficiencyCalculator& calc, const GeometryDescriptor& descriptor,
        std::vector<std::unique_ptr<Material>>& owned_materials);
};

} // namespace ceelo

#endif // CEELO_IO_RESPONSE_GENERATOR_H
