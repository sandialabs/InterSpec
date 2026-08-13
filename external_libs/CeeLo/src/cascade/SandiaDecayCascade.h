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

// Optional SandiaDecay -> DecayCascade adapter.
//
// Builds the engine's data-agnostic DecayCascade list (CascadeTypes.h) from a
// SandiaDecay nuclide. Compiled only when -DCEELO_WITH_SANDIADECAY=ON so the
// core library never hard-depends on SandiaDecay. A host application that
// already links SandiaDecay (e.g. InterSpec) can use this directly; otherwise
// the host builds DecayCascade lists itself.
//
// Coincidence handling:
//  - gamma-gamma: copied verbatim from SandiaDecay's RadParticle::coincidences
//    (the only coincidences the data carries), with product indices remapped to
//    member indices.
//  - gamma-xray (transition-group model): each decay/IC x-ray of a transition is
//    linked as coincident with that transition's gammas, since atomic relaxation
//    accompanies the same decay event. This is approximate (near-exact for
//    electron-capture emitters, which dominate x-ray summing).
//  - beta+ : each positron is expanded into a 511 keV annihilation member
//    coincident with the branch gammas; the engine emits the back-to-back pair.

#include "cascade/CascadeTypes.h"
#include "cascade/BetaSpectrum.h"

#include <string>
#include <vector>

namespace SandiaDecay {
struct Nuclide;
class SandiaDecayDataBase;
}

namespace ceelo {
namespace cascade_adapter {

/// Options controlling translation of a nuclide's decay data into DecayCascades.
struct CascadeOptions {
    /// Source age (seconds). Controls decay-chain progeny in-growth/evolution.
    double age_seconds = 0.0;
    /// If true, add the parent in prompt/secular equilibrium with its short-lived
    /// daughters (so e.g. Cs-137's 661 keV gamma, emitted by Ba-137m, is present
    /// even at age 0). If false, start from a pure parent and let progeny grow in
    /// over age_seconds.
    bool prompt_equilibrium = true;
    /// Include decay/IC x-rays (transition-group coincidence model). Ignored for
    /// a branch when the vacancy-level model is active (see vacancy_xray_model).
    bool include_xrays = true;
    /// Use the vacancy-level K x-ray model when the enriched ICC/EC data is
    /// present: instead of flat per-line Xray members, the branch carries K-shell
    /// vacancy sources (electron capture + internal conversion) and the engine
    /// generates one K x-ray per fired vacancy (enforcing K-line exclusivity and
    /// the right per-decay multiplicity). Falls back to the transition-group
    /// model when the data is absent. Set false to force the legacy model.
    bool vacancy_xray_model = true;
    /// Expand positrons into a 511 keV annihilation member.
    bool include_annihilation = true;
    /// Carry gamma-gamma angular-correlation coefficients (a2/a4) from the
    /// enriched SandiaDecay data onto the cascade links, so the engine samples
    /// W(theta) instead of isotropic partner directions. When the data is
    /// un-enriched (a2=a4=0) this has no effect. Set false to force isotropic.
    bool angular_correlations = true;
    /// Drop members whose per-branch intensity is below this.
    double min_intensity = 0.0;
    /// Drop branches whose (relative) selection weight is below this.
    double min_branch_weight = 0.0;
};

/// Per-branch report on the level scheme's flux balance and on the decay-data
/// defects that had to be repaired to obtain it. Optional diagnostic output of
/// build_cascades(); nothing in the engine consumes it.
///
/// A decay enters its daughter's level scheme at most once, so the invariant is
///     total_feed  <=  1      and      total_feed == total_sink
/// where total_feed is the selected (possibly common-scaled) sum of levels with net out-flow
/// (the direct beta/EC/alpha
/// feeding) and total_sink sums those with net in-flow. The identity is exact:
/// summing (out - in) over the closed level graph telescopes to zero.
///
/// total_feed is BELOW 1 whenever some decays populate the daughter's ground (or
/// isomeric) state directly -- e.g. Th-234 -> Pa-234m feeds the 73.9 keV isomer
/// in ~72% of decays, so 0.286 is correct, not a defect.
///
/// `n_sinks` counts levels with net in-flow. It is NOT an error indicator on its
/// own: many valid schemes have more than one, including validated examples
/// (Am-241, Eu-152, and Co-57), because a cascade legitimately
/// terminates at a metastable level as well as the ground state. Read it
/// alongside raw_total_feed. Accepted values in (1, 1.05] are common-scaled to
/// total_feed == 1; larger values reject the scheme.
struct BranchFeedingAudit {
    std::string parent, child;
    double branch_weight = 0.0;
    double raw_total_feed = 0.0; ///< feed sum before the accepted <=5% repair
    double feeding_scale = 1.0;  ///< common transition/feed scale (<= 1)
    double total_feed = 0.0;   ///< selected sum after feeding_scale
    double total_sink = 0.0;   ///< selected sink sum after feeding_scale
    int n_sinks = 0;           ///< levels with net in-flow (multiple may be physical)
    int n_e0_repaired = 0;         ///< definitive E0 records repaired
    int n_intensity_capped = 0;    ///< transitions whose I_gamma*(1+alpha) exceeded 1
    int n_split_renormalized = 0;  ///< transitions whose outcome split summed past 1
    int n_graph_e0_repaired = 0;       ///< E0 repairs on matched DAG edges
    int n_graph_intensity_capped = 0;  ///< intensity caps on matched DAG edges
    int n_graph_split_renormalized = 0;///< split repairs on matched DAG edges
    int n_residual_e0_repaired = 0;       ///< E0 repairs on accepted residuals
    int n_residual_intensity_capped = 0;  ///< intensity caps on accepted residuals
    int n_residual_split_renormalized = 0;///< split repairs on accepted residuals
    int n_gammas_unmatched = 0;    ///< kept as members but absent from the level graph
    double unmatched_intensity = 0.0;  ///< intensity those unmatched gammas carry
    int n_residual_transitions = 0; ///< unmatched categorical transitions retained
    double residual_transition_occurrence_mass = 0.0; ///< sum(1-p_none); may exceed 1
    double residual_gamma_probability = 0.0; ///< sum absolute residual p_gamma
    int n_suppressed_duplicates = 0; ///< explicit enriched duplicate records omitted
    double matched_gamma_intensity = 0.0; ///< legacy raw kept intensity placed on the DAG
    double total_gamma_intensity = 0.0; ///< legacy raw kept matched + residual intensity
    double matched_transition_occurrence_mass = 0.0; ///< sum repaired occurrence weights with topology
    double total_transition_occurrence_mass = 0.0; ///< sum repaired weights, including E0; may exceed 1
    bool partial_exact_feed_compatible = false; ///< <99% occurrence graph passed feed replay
    bool partial_scheme_rejected = false; ///< partial graph failed replay/residual bound
    int n_invalid_topology_edges = 0; ///< non-descending candidate edges; any rejects
    bool scheme_valid = false;     ///< a usable level scheme was produced
    bool rejected = false;         ///< raw flux feed sum >1.05; distinct from partial rejection
};

/// Build the correlated DecayCascade list for a nuclide (and its in-grown
/// progeny at the requested age). When `audit` is non-null it receives one entry
/// per decay branch considered (see BranchFeedingAudit).
std::vector<DecayCascade> build_cascades(const SandiaDecay::Nuclide* parent,
                                         const CascadeOptions& opts = {});

/// Audit-producing overload. Kept separate so rebuilt source consumers retain
/// the historical two-argument call and link symbol. This change adds fields to
/// returned cascade structs and therefore does NOT preserve binary ABI for
/// objects compiled against older headers; such consumers must be rebuilt.
std::vector<DecayCascade> build_cascades(
    const SandiaDecay::Nuclide* parent,
    const CascadeOptions& opts,
    std::vector<BranchFeedingAudit>* audit);

/// Convenience: look the nuclide up by label ("Co60", "Cs-137", ...).
std::vector<DecayCascade> build_cascades(const SandiaDecay::SandiaDecayDataBase& db,
                                         const std::string& nuclide_label,
                                         const CascadeOptions& opts = {});

std::vector<DecayCascade> build_cascades(
    const SandiaDecay::SandiaDecayDataBase& db,
    const std::string& nuclide_label,
    const CascadeOptions& opts,
    std::vector<BranchFeedingAudit>* audit);

/// Build the direct radioactive source term for electron/brems source-escape
/// studies. Unlike build_cascades(), this intentionally uses only the named
/// parent's immediate transitions: no daughter in-growth, prompt equilibrium,
/// or coincidence realization is applied. Legacy beta-product vectors whose
/// conditional sum exceeds one are common-scaled per transition. If inconsistent
/// duplicate/split parent branches still make the aggregate beta+positron yield
/// exceed one by more than 1e-6, this fails closed with std::runtime_error rather
/// than silently scaling authoritative exact level-feed branches. Excess at or
/// below 1e-6 is treated as float roundoff and closed exactly to one. Also throws
/// std::runtime_error for a malformed host-constructed exact beta level-feed law.
RadioactiveEmissionSet build_radioactive_emissions(
    const SandiaDecay::Nuclide* parent);

/// Convenience lookup overload.
RadioactiveEmissionSet build_radioactive_emissions(
    const SandiaDecay::SandiaDecayDataBase& db,
    const std::string& nuclide_label);

} // namespace cascade_adapter
} // namespace ceelo
