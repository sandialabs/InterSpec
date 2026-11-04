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

#include "transport/PhotonTransport.h"
#include "transport/TransportUtils.h"
#include "transport/ComptonScatter.h"
#include "physics/ElectronCsda.h"
#include "geometry/Geometry.h"
#include "geometry/Cylinder.h"
#include "geometry/Box.h"
#include "geometry/RayTrace.h"
#include "materials/Material.h"
#include "cross_sections/CrossSectionData.h"

#include <cmath>
#include <cassert>
#include <algorithm>
#include <random>
#include <limits>

namespace ceelo {

namespace {
constexpr double kPi    = 3.14159265358979323846;
constexpr double kTwoPi = 2.0 * kPi;

// Polar angle (returns cosθ) of a pair-production lepton relative to the parent
// photon direction, sampled from the G4ModifiedTsai distribution: θ = u·mₑc²/E
// with u = −ln(ξ₁ξ₂) scaled by a1 (prob 0.25) or a2 = a1/3 (prob 0.75), capped
// at u = E·π/mₑc² so θ ≤ π.  E is the lepton TOTAL energy (keV).  The
// characteristic opening is ~mₑc²/E — small at MeV, but it spreads the e⁺/e⁻
// (and hence the positron annihilation site) off the photon axis.
double sample_pair_lepton_cos_theta(double E_total_keV, std::mt19937_64& rng) {
    constexpr double mc2 = 510.998950;
    constexpr double a1 = 1.6, a2 = a1 / 3.0, border = 0.25;
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    const double uMax = E_total_keV * kPi / mc2;
    double u;
    do {
        double sel = u01(rng);
        u = -std::log(u01(rng) * u01(rng));
        u *= (sel < border) ? a1 : a2;
    } while (u > uMax);
    return std::cos(u * mc2 / E_total_keV);
}

// Summarize one annihilation photon's sub-history for the opt-in pair
// diagnostics. Reads only values the sub-transport already produced, so it is
// pure bookkeeping: no random draws, no effect on any transport result.
void record_annihilation_photon(
    TransportResult::PairHistory::AnnihilationPhoton& record,
    double emitted_energy_keV,
    const TransportResult& sub)
{
    record.emitted_energy_keV = emitted_energy_keV;
    record.deposited_in_scoring_keV = sub.energy_deposited_scoring;
    record.escaped = sub.escaped;
    record.exit_energy_keV = sub.escaped ? sub.exit_energy_keV : 0.0;
    // "Scattered before escape" counts only Compton vertices inside the
    // scoring volume, matching the convention the detector-response campaign
    // uses for post-pair Compton counting.
    if (sub.escaped) {
        for (const auto& interaction : sub.interactions) {
            if (interaction.in_scoring
                && interaction.type == InteractionType::Compton) {
                record.scattered_before_escape = true;
                break;
            }
        }
    }
}
} // anonymous namespace


ElectronDepositResult walk_electron(const Eigen::Vector3d& position,
                                    const Eigen::Vector3d& direction,
                                    double KE_keV,
                                    const Geometry& geometry,
                                    const TransportConfig& config,
                                    std::mt19937_64& rng,
                                    bool is_positron)
{
    const auto& csda = ElectronCsda::instance();
    if (config.electron_transport == ElectronTransport::Detailed && !config.disable_moliere)
        return csda.deposited_in_scoring_detailed(
            position, direction, KE_keV, geometry, rng, config.disable_brems, is_positron,
            config.enable_impact_ionization ? ElectronCsda::impact_ionization_max_keV() : 0.0);
    return csda.deposited_in_scoring(position, direction, KE_keV, geometry, rng,
                                     config.disable_moliere, config.disable_brems,
                                     is_positron);
}


AnnihilationPair sample_annihilation_pair(bool in_flight, double residual_KE_keV,
                                          std::mt19937_64& rng)
{
    constexpr double kElectronMassKeV = 510.998950;
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    AnnihilationPair pair;
    if (in_flight) {
        const double Etot = 2.0 * kElectronMassKeV + residual_KE_keV;
        const double x    = 0.5 + (uniform(rng) - 0.5) * (residual_KE_keV / Etot);
        pair.energy_keV   = {x * Etot, Etot - x * Etot};
        pair.direction[0] = sample_isotropic_direction(rng);
        pair.direction[1] = sample_isotropic_direction(rng);
    } else {
        pair.energy_keV   = {kElectronMassKeV, kElectronMassKeV};
        pair.direction[0] = sample_isotropic_direction(rng);
        pair.direction[1] = -pair.direction[0];
    }
    return pair;
}


TransportResult transport_photon(
    Eigen::Vector3d position,
    Eigen::Vector3d direction,
    double energy_keV,
    const Geometry& geometry,
    const TransportConfig& config_in,
    std::mt19937_64& rng)
{
    // Forced first interaction applies to this history only: clear the flag
    // in the config used for the rest of transport so recursive calls
    // (fluorescence, brems, annihilation gammas) are never re-forced.
    const bool do_force = config_in.force_first_interaction;
    TransportConfig forced_cfg;
    if (do_force) {
        forced_cfg = config_in;
        forced_cfg.force_first_interaction = false;
    }
    const TransportConfig& config = do_force ? forced_cfg : config_in;

    TransportResult result{};
    result.energy_deposited_scoring = 0.0;
    result.energy_deposited_total = 0.0;
    result.any_interaction_in_scoring = false;
    result.num_compton_scatters = 0;
    result.forced_absorption = false;

    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    // FEP early termination: capture initial energy so we can detect when
    // enough energy has escaped scoring to make FEP impossible.
    const double initial_energy_keV = energy_keV;

    // Returns true if FEP is impossible: even if all remaining energy
    // (max_remaining_keV) deposits in scoring, the deficit exceeds tolerance.
    auto fep_killed = [&](double max_remaining_keV) -> bool {
        return config.fep_only_mode && !config.disable_fep_early_kill &&
               (initial_energy_keV - result.energy_deposited_scoring - max_remaining_keV)
                    > config.fep_window_keV;
    };

    // An electron-impact K x-ray among an electron's brems photons is
    // recorded like a fluorescence x-ray, so its escape is labelled.
    auto record_impact_xray = [&](const BremsPhoton& bp, const TransportResult& br) {
        if (bp.impact_Z == 0 || !config.record_interactions) return;
        FluorescenceRecord record{
            bp.impact_Z, true, true, bp.energy_keV,
            br.energy_deposited_scoring < bp.energy_keV * 0.99,
            br.energy_deposited_scoring};
        record.filling_subshell = bp.impact_filling_subshell;
        record.from_impact = true;
        result.fluorescence_records.push_back(record);
    };

    // Transport one relaxation x-ray (K or L line) isotropically from `from`
    // and fold its deposits and escapes into this history. Fluorescence may
    // cascade (e.g. Te Ka -> Cd K-shell PE -> Cd Ka) to max_fluorescence_depth.
    auto transport_relaxation_xray = [&](const Eigen::Vector3d& from,
                                         double x_keV) -> TransportResult {
        Eigen::Vector3d x_dir = sample_isotropic_direction(rng);
        TransportConfig fl_config = config;
        fl_config.fluorescence_depth = config.fluorescence_depth + 1;
        if (fl_config.fluorescence_depth >= config.max_fluorescence_depth)
            fl_config.enable_fluorescence = false;
        auto fl_res = transport_photon(from, x_dir, x_keV, geometry, fl_config, rng);
        result.energy_deposited_scoring += fl_res.energy_deposited_scoring;
        result.energy_deposited_total   += fl_res.energy_deposited_total;
        if (fl_res.any_interaction_in_scoring)
            result.any_interaction_in_scoring = true;
        if (fl_res.escaped && fl_res.exit_energy_keV > 1.0)
            result.escaped_secondaries.push_back({
                fl_res.exit_position, fl_res.exit_direction,
                fl_res.exit_energy_keV});
        result.escaped_secondaries.insert(
            result.escaped_secondaries.end(),
            fl_res.escaped_secondaries.begin(),
            fl_res.escaped_secondaries.end());
        return fl_res;
    };

    double energy_MeV = energy_keV * 1e-3;

    // Maximum number of transport steps to prevent infinite loops
    const int max_steps = 500;
    std::vector<PathSegment> segments;
    std::vector<PathSegment> next_segs;

    // ------------------------------------------------------------------
    // Forced first interaction (variance reduction).
    // Sample the first interaction point from the truncated exponential
    // over the total optical depth tau along the entry ray, and set
    // result.weight = 1 - exp(-tau) (the analog interaction probability).
    // The pass-through branch scores 0 in all tallies, so removing it and
    // carrying its probability as a weight is exactly unbiased.
    // ------------------------------------------------------------------
    bool forced_pending = false;
    if (do_force) {
        geometry.trace_ray(position, direction, segments);

        // Per-segment optical depth with the same Rayleigh/PP masking as
        // the main loop. In FEP-only mode, stop at the first non-scoring
        // material segment after a scoring one: the analog photon would be
        // killed at that boundary, so those segments cannot contribute.
        auto segment_tau = [&](const PathSegment& seg, double& mu, double& len) -> double {
            mu = 0.0;
            len = seg.t_end - std::max(seg.t_start, 0.0);
            if (!seg.material || len <= 0.0) return 0.0;
            MacroscopicXS xs = seg.material->macroscopic_xs(energy_MeV);
            if (!config.enable_rayleigh) xs.mu_rs = 0.0;
            if (!config.enable_pair_production || energy_keV < 1022.0) xs.mu_pp = 0.0;
            mu = xs.mu_total();
            return (mu > 0.0) ? mu * len : 0.0;
        };

        double tau_total = 0.0;
        bool seen_scoring = false;
        size_t num_ladder_segs = 0;
        for (const auto& seg : segments) {
            if (config.fep_only_mode && seen_scoring &&
                seg.material && !seg.is_scoring) {
                break;
            }
            if (seg.material && seg.is_scoring) seen_scoring = true;
            double mu, len;
            tau_total += segment_tau(seg, mu, len);
            ++num_ladder_segs;
        }

        if (tau_total > 1e-12) {
            const double w_fc = -std::expm1(-tau_total);  // 1 - exp(-tau)
            result.weight = w_fc;

            // Truncated-exponential optical depth in (0, tau_total)
            const double t_forced = -std::log1p(-uniform(rng) * w_fc);

            double tau_cum = 0.0;
            for (size_t i = 0; i < num_ladder_segs; ++i) {
                const auto& seg = segments[i];
                double mu, len;
                const double tau_seg = segment_tau(seg, mu, len);
                if (tau_seg <= 0.0) continue;
                const bool is_last = (tau_cum + tau_seg >= tau_total * (1.0 - 1e-12));
                if (t_forced <= tau_cum + tau_seg || is_last) {
                    double d = (t_forced - tau_cum) / mu;
                    // Keep strictly inside the segment for numerical safety
                    d = std::min(std::max(d, 0.0), len * (1.0 - 1e-12));
                    position += direction * (std::max(seg.t_start, 0.0) + d);
                    forced_pending = true;
                    break;
                }
                tau_cum += tau_seg;
            }
        }
        // If tau_total ~ 0 (no material along the ray), fall through to
        // normal transport: the photon escapes with weight 1.
    }

    for (int step = 0; step < max_steps; ++step) {
        // Minimum energy check
        if (energy_keV < config.min_energy_keV) {
            // Deposit remaining energy at current location
            // Check if we're in a scoring region
            geometry.trace_ray(position, direction, segments);
            if (!segments.empty() && segments[0].is_scoring && segments[0].t_start <= 1e-6) {
                result.energy_deposited_scoring += energy_keV;
                result.any_interaction_in_scoring = true;
            }
            result.energy_deposited_total += energy_keV;
            break;
        }

        // Trace ray from current position to find what material we're in
        // and how far to the next boundary
        geometry.trace_ray(position, direction, segments);

        if (segments.empty()) {
            // Photon has escaped the geometry
            result.escaped = true;
            result.exit_position = position;
            result.exit_direction = direction;
            result.exit_energy_keV = energy_keV;
            if (config.record_interactions) {
                result.escape_energy_keV = energy_keV;
            }
            break;
        }

        // Find the segment we're currently in (t_start closest to 0)
        const PathSegment* current_seg = nullptr;
        for (const auto& seg : segments) {
            if (seg.t_start < 1e-8 || (seg.t_end > 1e-8 && seg.t_start <= 1e-8)) {
                current_seg = &seg;
                break;
            }
        }

        if (!current_seg) {
            // We're not inside any material — we're in a gap (e.g., bore hole)
            // Find the next segment we'd enter
            for (const auto& seg : segments) {
                if (seg.t_start > 1e-8) {
                    // Advance to the start of the next segment
                    position += direction * (seg.t_start + 1e-9);
                    current_seg = &seg;
                    break;
                }
            }
            if (!current_seg) {
                // No more segments — photon escapes
                result.escaped = true;
                result.exit_position = position;
                result.exit_direction = direction;
                result.exit_energy_keV = energy_keV;
                break;
            }
            continue; // Re-trace from new position
        }

        if (!current_seg->material) {
            result.escaped = true;
            result.exit_position = position;
            result.exit_direction = direction;
            result.exit_energy_keV = energy_keV;
            break; // Vacuum — photon escapes
        }

        // Compute macroscopic cross-sections in the current material
        energy_MeV = energy_keV * 1e-3;
        MacroscopicXS xs = current_seg->material->macroscopic_xs(energy_MeV);

        if (!config.enable_rayleigh) xs.mu_rs = 0.0;
        if (!config.enable_pair_production || energy_keV < 1022.0) xs.mu_pp = 0.0;

        double mu_total = xs.mu_total();
        if (mu_total <= 0.0) {
            // No interactions possible — advance to boundary
            position += direction * (current_seg->t_end + 1e-9);
            continue;
        }

        // Sample interaction distance. A pending forced first interaction
        // occurs exactly at the current (pre-sampled) position.
        double s;
        if (forced_pending) {
            forced_pending = false;
            s = 0.0;
        } else {
            s = sample_interaction_distance(mu_total, rng);
        }

        // Distance to the end of the current segment
        double d_boundary = current_seg->t_end - std::max(current_seg->t_start, 0.0);
        // If we entered mid-segment, adjust
        if (current_seg->t_start < 0.0) {
            d_boundary = current_seg->t_end;
        }

        if (s >= d_boundary) {
            // Photon reaches the boundary without interacting.
            // In FEP-only mode: if leaving a scoring segment, kill the event.
            if (config.fep_only_mode && current_seg->is_scoring) {
                // Check if the next region is non-scoring or outside geometry
                Eigen::Vector3d next_pos = position + direction * (current_seg->t_end + 1e-9);
                geometry.trace_ray(next_pos, direction, next_segs);
                bool next_is_scoring = false;
                for (const auto& ns : next_segs) {
                    if (ns.t_start < 1e-8 || (ns.t_end > 1e-8 && ns.t_start <= 1e-8)) {
                        next_is_scoring = ns.is_scoring;
                        break;
                    }
                }
                if (!next_is_scoring) {
                    // Energy escaping scoring volume — kill for FEP-only
                    result.energy_deposited_scoring = 0.0;
                    result.any_interaction_in_scoring = false;
                    goto done;
                }
            }
            // Advance to just past the boundary
            position += direction * (current_seg->t_end + 1e-9);
            continue;
        }

        // Interaction occurs at distance s
        position += direction * s;

        // Select interaction type
        double xi = uniform(rng);
        InteractionType interaction = select_interaction(xs, xi);

        bool in_scoring = current_seg->is_scoring;

        switch (interaction) {

        // ----------------------------------------------------------------
        // Photoelectric absorption
        // ----------------------------------------------------------------
        case InteractionType::Photoelectric: {
            if (config.record_interactions) {
                result.interactions.push_back({
                    InteractionType::Photoelectric, position,
                    energy_keV, 0.0, 1.0, in_scoring,
                    result.num_compton_scatters});
            }
            double deposited = energy_keV; // default: deposit all energy
            double pe_kinetic_energy = energy_keV; // photoelectron KE for CSDA
            double auger_local = 0.0; // Auger energy deposited locally

            if (config.enable_fluorescence) {
                // Select which element was photoelectrically absorbed.
                // Use the material's per-element PE XS weights.
                int Z_pe = current_seg->material->select_element(energy_MeV, 0 /*PE*/, uniform(rng));
                const FluorescenceData* fl = CrossSectionData::instance().fluorescence(Z_pe);
                const LFluorescenceData* fl_L = CrossSectionData::instance().l_fluorescence(Z_pe);
                // Only elements with retained (>= 10 keV) L lines take the L
                // draws; for the rest an L vacancy stays local as before.
                const bool l_lines = has_l_lines(fl_L);

                // Photoabsorption outside the K shell: an L-subshell vacancy
                // (by the tabulated shares) radiates its retained L x-ray, the
                // photoelectron carries E - E_L. Returns true if the FEP early
                // kill fires.
                auto outer_shell_absorption = [&](double f_K) -> bool {
                    if (!l_lines) return false;
                    const int s_l = sample_l_photoelectric_subshell(*fl_L, energy_keV, f_K, rng);
                    if (s_l < 0 || energy_keV <= fl_L->sub[s_l].edge_keV) return false;
                    pe_kinetic_energy = energy_keV - fl_L->sub[s_l].edge_keV;
                    const LRelaxation l_relax = sample_l_relaxation(*fl_L, s_l, rng);
                    auger_local = l_relax.local_keV;
                    deposited = pe_kinetic_energy + auger_local;
                    if (l_relax.line_keV <= 0.0) return false;
                    transport_relaxation_xray(position, l_relax.line_keV);
                    return fep_killed(pe_kinetic_energy + auger_local);
                };

                if (fl && energy_keV > fl->k_edge_keV) {
                    // Probabilistic K-shell vs outer-shell PE selection.
                    // f_K = σ_K / σ_PE_total.
                    const auto& xs_data = CrossSectionData::instance();
                    double sigma_K = xs_data.sigma_K_photoelectric(Z_pe, energy_MeV);
                    double sigma_tot = xs_data.sigma_photoelectric(Z_pe, energy_MeV);
                    double f_K = (sigma_tot > 0.0) ? std::min(sigma_K / sigma_tot, 1.0) : 1.0;

                    if (uniform(rng) >= f_K) {
                        // Outer-shell PE: full energy deposited locally unless an
                        // L x-ray carries some of it away.
                        if (config.record_interactions) {
                            result.fluorescence_records.push_back({
                                Z_pe, false, false, 0.0, false, 0.0});
                        }
                        if (outer_shell_absorption(f_K)) {
                            energy_keV = 0.0;
                            goto done;
                        }
                        goto pe_deposit;
                    }

                    // K-shell absorption:
                    //   photoelectron KE = E_γ − E_K
                    pe_kinetic_energy = energy_keV - fl->k_edge_keV;
                    deposited = pe_kinetic_energy;

                    // K x-ray (Kα1, Kα2, Kβ, ... by intensity) or Auger electrons.
                    // The binding energy the x-ray does not carry goes to Auger
                    // electrons with sub-mm range: deposited locally. A K-alpha
                    // leaves an L vacancy, which may radiate a retained L line.
                    const KRelaxation relax = sample_k_relaxation(*fl, rng);
                    auger_local = relax.local_keV;
                    LRelaxation l_relax;
                    if (l_lines && relax.filling_subshell >= kEadlL1
                        && relax.filling_subshell <= kEadlL3) {
                        l_relax = sample_l_relaxation(
                            *fl_L, relax.filling_subshell - kEadlL1, rng);
                        auger_local -= l_relax.line_keV;
                    }
                    deposited += auger_local;
                    if (relax.line_keV > 0.0) {
                        // --- K-shell fluorescence X-ray ---
                        const double fl_energy = relax.line_keV;
                        auto fl_res = transport_relaxation_xray(position, fl_energy);

                        // Record fluorescence fate
                        if (config.record_interactions) {
                            bool escaped = (fl_res.energy_deposited_scoring <
                                            static_cast<double>(fl_energy) * 0.99);
                            FluorescenceRecord record{
                                Z_pe, true, true,
                                static_cast<double>(fl_energy),
                                escaped, fl_res.energy_deposited_scoring};
                            record.filling_subshell = relax.filling_subshell;
                            result.fluorescence_records.push_back(record);
                        }

                        // FEP early kill: fluorescence X-ray escaped scoring.
                        // Best case remaining = photoelectron KE + Auger (+ the
                        // L x-ray still to come), none deposited yet.
                        if (fep_killed(pe_kinetic_energy + auger_local + l_relax.line_keV)) {
                            energy_keV = 0.0;
                            goto done;
                        }
                        if (l_relax.line_keV > 0.0) {
                            transport_relaxation_xray(position, l_relax.line_keV);
                            if (fep_killed(pe_kinetic_energy + auger_local)) {
                                energy_keV = 0.0;
                                goto done;
                            }
                        }

                    } else {
                        // --- Auger electrons ---
                        // Binding energy released locally via Auger electrons
                        // (auger_local = E_K above): total deposited = E_γ.
                        if (config.record_interactions) {
                            result.fluorescence_records.push_back({
                                Z_pe, true, false, 0.0, false, 0.0});
                        }
                    }
                } else {
                    // Below the K edge or no K fluorescence data.
                    if (config.record_interactions && fl) {
                        result.fluorescence_records.push_back({
                            Z_pe, false, false, 0.0, false, 0.0});
                    }
                    if (outer_shell_absorption(0.0)) {
                        energy_keV = 0.0;
                        goto done;
                    }
                }
                // Outer-shell absorption in an element without retained L lines
                // (and M and beyond everywhere) leaves deposited = energy_keV and
                // pe_kinetic_energy = energy_keV.
            }

            pe_deposit:
            // CSDA electron tracking: photoelectron KE only (Auger deposited locally).
            if (config.enable_electron_csda && pe_kinetic_energy > 0.0) {
                Eigen::Vector3d pe_dir = sample_photoelectron_direction(
                    direction, pe_kinetic_energy, rng);
                auto e_result = walk_electron(position, pe_dir, pe_kinetic_energy,
                                              geometry, config, rng);
                result.energy_deposited_scoring += e_result.deposited_scoring_keV;
                if (e_result.deposited_scoring_keV > 0.0) result.any_interaction_in_scoring = true;

                // Track electron escape
                if (config.record_interactions) {
                    double e_escape = pe_kinetic_energy - e_result.deposited_scoring_keV;
                    if (e_escape > 0.01) result.electron_escape_keV += e_escape;
                }

                // Auger electrons deposit locally (sub-mm range)
                if (in_scoring && auger_local > 0.0) {
                    result.energy_deposited_scoring += auger_local;
                    result.any_interaction_in_scoring = true;
                }

                // Transport bremsstrahlung photons
                for (const auto& bp : e_result.brems_photons) {
                    auto br = transport_photon(bp.position, bp.direction, bp.energy_keV,
                                               geometry, config, rng);
                    record_impact_xray(bp, br);
                    result.energy_deposited_scoring += br.energy_deposited_scoring;
                    result.energy_deposited_total   += br.energy_deposited_total;
                    if (br.any_interaction_in_scoring) result.any_interaction_in_scoring = true;
                    if (br.escaped && br.exit_energy_keV > 1.0)
                        result.escaped_secondaries.push_back({
                            br.exit_position, br.exit_direction, br.exit_energy_keV});
                    result.escaped_secondaries.insert(
                        result.escaped_secondaries.end(),
                        br.escaped_secondaries.begin(), br.escaped_secondaries.end());
                }
            } else {
                if (in_scoring) {
                    result.energy_deposited_scoring += deposited;
                    result.any_interaction_in_scoring = true;
                }
            }
            // Total energy deposited in the system is conserved regardless of CSDA
            result.energy_deposited_total += deposited;
            energy_keV = 0.0;
            goto done; // photon absorbed
        }

        // ----------------------------------------------------------------
        // Compton scattering
        // ----------------------------------------------------------------
        case InteractionType::Compton: {
            // Count only in-crystal scatters: matches the header field doc
            // ("Number of Compton scatters in the crystal") and the
            // force-absorption cutoff below, which is gated on in_scoring.
            // Non-scoring scatters (attenuator/source-shield) must not push
            // the count toward the cap.
            if (in_scoring) result.num_compton_scatters++;

            // Check Compton scatter cutoff (only counts inside crystal)
            if (in_scoring && result.num_compton_scatters >= config.max_compton_scatters) {
                // Note: no interaction record for forced absorption (not a real interaction)
                // Force absorption: deposit ALL remaining energy
                result.energy_deposited_scoring += energy_keV;
                result.energy_deposited_total   += energy_keV;
                result.any_interaction_in_scoring = true;
                result.forced_absorption = true;
                energy_keV = 0.0;
                goto done;
            }

            // Select element for bound-electron Compton correction (S(x,Z))
            int Z_cs = current_seg->material->select_element(energy_MeV, 1 /*CS*/, uniform(rng));

            // Sample Compton scatter with bound-electron correction and
            // (optionally) Doppler broadening of the scattered energy.
            // Save incoming direction for recoil direction calculation
            Eigen::Vector3d dir_in = direction;
            double energy_in = energy_keV;
            auto scatter = sample_compton_scatter(energy_keV, direction, rng,
                                                  Z_cs, config.enable_doppler_broadening);

            if (config.record_interactions) {
                result.interactions.push_back({
                    InteractionType::Compton, position,
                    energy_in, scatter.scattered_energy_keV, scatter.cos_theta,
                    in_scoring, result.num_compton_scatters - 1});
            }

            // CSDA electron tracking: recoil electron direction from 4-momentum.
            // (Direction uses the free-electron momentum balance; the Doppler
            // momentum shift in direction is negligible for mm-range electrons.)
            if (config.enable_electron_csda && scatter.deposited_energy_keV > 0.0) {
                Eigen::Vector3d e_dir = compton_recoil_direction(
                    dir_in, scatter.new_direction, energy_in, scatter.scattered_energy_keV);
                auto e_result = walk_electron(position, e_dir, scatter.deposited_energy_keV,
                                              geometry, config, rng);
                result.energy_deposited_scoring += e_result.deposited_scoring_keV;
                if (in_scoring || e_result.deposited_scoring_keV > 0.0)
                    result.any_interaction_in_scoring = true;

                // Track electron escape
                if (config.record_interactions) {
                    double e_escape = scatter.deposited_energy_keV - e_result.deposited_scoring_keV;
                    if (e_escape > 0.01) result.electron_escape_keV += e_escape;
                }

                // Transport bremsstrahlung photons
                for (const auto& bp : e_result.brems_photons) {
                    auto br = transport_photon(bp.position, bp.direction, bp.energy_keV,
                                               geometry, config, rng);
                    record_impact_xray(bp, br);
                    result.energy_deposited_scoring += br.energy_deposited_scoring;
                    result.energy_deposited_total   += br.energy_deposited_total;
                    if (br.any_interaction_in_scoring) result.any_interaction_in_scoring = true;
                    if (br.escaped && br.exit_energy_keV > 1.0)
                        result.escaped_secondaries.push_back({
                            br.exit_position, br.exit_direction, br.exit_energy_keV});
                    result.escaped_secondaries.insert(
                        result.escaped_secondaries.end(),
                        br.escaped_secondaries.begin(), br.escaped_secondaries.end());
                }
            } else {
                if (in_scoring) {
                    result.energy_deposited_scoring += scatter.deposited_energy_keV;
                    result.any_interaction_in_scoring = true;
                }
            }
            result.energy_deposited_total += scatter.deposited_energy_keV;

            // Doppler: the struck subshell's binding energy U_i.  A K vacancy
            // relaxes like a photoelectric one when enable_compton_fluorescence
            // is set: its K x-ray is transported and only U_i - E_line stays
            // local.  Otherwise (and for outer shells, whose L/M lines are
            // deposited locally anyway) all of U_i is deposited on the spot.
            if (scatter.binding_deposit_keV > 0.0) {
                double local = scatter.binding_deposit_keV;
                const bool relax_vacancy =
                    config.enable_compton_fluorescence && config.enable_fluorescence;
                const FluorescenceData* fl = (relax_vacancy && scatter.shell_index == 0)
                    ? CrossSectionData::instance().fluorescence(Z_cs) : nullptr;
                const LFluorescenceData* fl_L = relax_vacancy
                    ? CrossSectionData::instance().l_fluorescence(Z_cs) : nullptr;
                const bool l_lines = has_l_lines(fl_L);
                if (fl && std::abs(scatter.binding_deposit_keV - fl->k_edge_keV)
                              < 0.05 * fl->k_edge_keV) {
                    const KRelaxation relax = sample_k_relaxation(*fl, rng);
                    if (relax.line_keV > 0.0
                        && relax.line_keV < scatter.binding_deposit_keV) {
                        local = scatter.binding_deposit_keV - relax.line_keV;
                        LRelaxation l_relax;
                        if (l_lines && relax.filling_subshell >= kEadlL1
                            && relax.filling_subshell <= kEadlL3) {
                            l_relax = sample_l_relaxation(
                                *fl_L, relax.filling_subshell - kEadlL1, rng);
                            local -= l_relax.line_keV;
                        }
                        auto fl_res = transport_relaxation_xray(position, relax.line_keV);
                        if (config.record_interactions) {
                            FluorescenceRecord record{
                                Z_cs, true, true, relax.line_keV,
                                fl_res.energy_deposited_scoring < relax.line_keV * 0.99,
                                fl_res.energy_deposited_scoring};
                            record.from_compton = true;
                            record.filling_subshell = relax.filling_subshell;
                            result.fluorescence_records.push_back(record);
                        }
                        if (l_relax.line_keV > 0.0)
                            transport_relaxation_xray(position, l_relax.line_keV);
                    }
                } else if (l_lines && scatter.shell_index >= 1 && scatter.shell_index <= 3) {
                    // L-shell Compton vacancy (xraylib shell order K, L1, L2, L3).
                    const int s_l = scatter.shell_index - 1;
                    const double edge = fl_L->sub[s_l].edge_keV;
                    if (std::abs(scatter.binding_deposit_keV - edge) < 0.05 * edge) {
                        const LRelaxation l_relax = sample_l_relaxation(*fl_L, s_l, rng);
                        if (l_relax.line_keV > 0.0
                            && l_relax.line_keV < scatter.binding_deposit_keV) {
                            local = scatter.binding_deposit_keV - l_relax.line_keV;
                            transport_relaxation_xray(position, l_relax.line_keV);
                        }
                    }
                }
                if (in_scoring) {
                    result.energy_deposited_scoring += local;
                    result.any_interaction_in_scoring = true;
                }
                result.energy_deposited_total += local;
            }

            // Update photon state
            energy_keV = scatter.scattered_energy_keV;
            direction  = scatter.new_direction;

            // FEP early kill: energy escaped via electron or brems escape.
            // Best case remaining = scattered photon energy.
            if (fep_killed(energy_keV)) {
                goto done;
            }
            break;
        }

        // ----------------------------------------------------------------
        // Rayleigh (coherent) scattering — form-factor sampling
        // ----------------------------------------------------------------
        case InteractionType::Rayleigh: {
            // Select element by Rayleigh XS contribution
            int Z_rs = current_seg->material->select_element(energy_MeV, 2 /*RS*/, uniform(rng));

            double cos_theta = sample_rayleigh_cos_theta(Z_rs, energy_keV, rng);

            if (config.record_interactions) {
                result.interactions.push_back({
                    InteractionType::Rayleigh, position,
                    energy_keV, energy_keV, cos_theta, in_scoring, -1});
            }

            double phi = kTwoPi * uniform(rng);
            direction = rotate_direction(direction, cos_theta, phi);
            // No energy deposition — Rayleigh is elastic
            break;
        }

        // ----------------------------------------------------------------
        // Pair production — track both 511 keV annihilation gammas
        // ----------------------------------------------------------------
        case InteractionType::PairProduction: {
            if (config.record_interactions) {
                result.interactions.push_back({
                    InteractionType::PairProduction, position,
                    energy_keV, 0.0, 1.0, in_scoring, -1});
            }
            constexpr double kPairThreshold = 1022.0; // keV

            // Kinetic energy of the e⁺/e⁻ pair: E − 1022 keV, split evenly.
            double deposited_local = energy_keV - kPairThreshold;
            if (deposited_local < 0.0) deposited_local = 0.0; // guard

            // CSDA tracking for the e⁺/e⁻ pair, with G4ModifiedTsai-sampled
            // lepton opening angles (sample_pair_lepton_cos_theta below).
            // Site of the back-to-back 511 keV annihilation pair.  Defaults to
            // the pair-production vertex; when the positron is CSDA-tracked it is
            // updated to the positron's rest point below (B1: the 511 source is
            // the positron endpoint, ~1.5 mm off the vertex at MeV in NaI, which
            // matters for single/double-escape-peak placement).
            Eigen::Vector3d annih_pos = position;
            // In-flight annihilation: when the positron annihilates before
            // stopping, the pair carries 2·mₑc² + residual KE (above 1022 keV,
            // off the clean 511/escape peaks) from the annihilation site.
            bool   pos_in_flight = false;
            double pos_residual_KE = 0.0;
            // A positron that walks out of the geometry annihilates outside it
            // (vacuum world): it emits no photon pair here.
            bool   pos_escaped = false;
            if (config.enable_electron_csda && deposited_local > 0.0) {
                double ke_each = 0.5 * deposited_local;
                // Pair opening angle (Koch-Motz / G4ModifiedTsai): each lepton is
                // emitted at ~mₑc²/E off the photon axis, on opposite azimuths
                // (transverse momentum of the absorbed photon is ~0).
                double E_lep   = ke_each + 510.998950;  // lepton total energy (keV)
                double cos_neg = sample_pair_lepton_cos_theta(E_lep, rng);
                double cos_pos = sample_pair_lepton_cos_theta(E_lep, rng);
                double phi     = kTwoPi * uniform(rng);
                Eigen::Vector3d dir_neg = rotate_direction(direction, cos_neg, phi);
                Eigen::Vector3d dir_pos = rotate_direction(direction, cos_pos, phi + kPi);
                auto e_neg = walk_electron(position, dir_neg, ke_each, geometry, config, rng);
                auto e_pos = walk_electron(position, dir_pos, ke_each, geometry, config, rng,
                                           /*is_positron=*/true);
                // The positron annihilates where it comes to rest (or in flight),
                // not at the creation vertex.  (The e⁻ does not annihilate.)
                annih_pos       = e_pos.stop_position;
                pos_in_flight   = e_pos.annihilated_in_flight;
                pos_residual_KE = e_pos.residual_KE_at_annih_keV;
                pos_escaped     = e_pos.escaped;
                double dep_pair = e_neg.deposited_scoring_keV + e_pos.deposited_scoring_keV;
                result.energy_deposited_scoring += dep_pair;
                if (dep_pair > 0.0) result.any_interaction_in_scoring = true;

                // Track electron escape.  In-flight residual KE left as photon
                // energy, not slowing-down deposit, so exclude it from "escape".
                if (config.record_interactions) {
                    double slowing = deposited_local - (pos_in_flight ? pos_residual_KE : 0.0);
                    double e_escape = slowing - dep_pair;
                    if (e_escape > 0.01) result.electron_escape_keV += e_escape;
                    if (!result.pair_history.occurred) {
                        auto& history = result.pair_history;
                        history.occurred = true;
                        history.in_scoring = in_scoring;
                        history.annihilated_in_flight = pos_in_flight;
                        history.positron_escaped = pos_escaped;
                        history.pair_kinetic_energy_keV = deposited_local;
                        history.lepton_deposit_in_scoring_keV = dep_pair;
                        double brems = 0.0;
                        for (const auto* e_res : {&e_neg, &e_pos})
                            for (const auto& bp : e_res->brems_photons)
                                brems += bp.energy_keV;
                        history.lepton_brems_energy_keV = brems;
                    }
                }

                // Transport bremsstrahlung photons from both e- and e+
                for (const auto* e_res : {&e_neg, &e_pos}) {
                    for (const auto& bp : e_res->brems_photons) {
                        auto br = transport_photon(bp.position, bp.direction, bp.energy_keV,
                                                   geometry, config, rng);
                        record_impact_xray(bp, br);
                        result.energy_deposited_scoring += br.energy_deposited_scoring;
                        result.energy_deposited_total   += br.energy_deposited_total;
                        if (br.any_interaction_in_scoring) result.any_interaction_in_scoring = true;
                        if (br.escaped && br.exit_energy_keV > 1.0)
                            result.escaped_secondaries.push_back({
                                br.exit_position, br.exit_direction, br.exit_energy_keV});
                        result.escaped_secondaries.insert(
                            result.escaped_secondaries.end(),
                            br.escaped_secondaries.begin(), br.escaped_secondaries.end());
                    }
                }
            } else {
                if (in_scoring && deposited_local > 0.0) {
                    result.energy_deposited_scoring += deposited_local;
                    result.any_interaction_in_scoring = true;
                }
            }
            // In-flight residual KE leaves as photon energy, not slowing-down
            // deposit, so don't also count it as locally deposited.
            result.energy_deposited_total += deposited_local
                                           - (pos_in_flight ? pos_residual_KE : 0.0);

            if (pos_escaped) {
                energy_keV = 0.0;
                goto done;
            }

            // Emit the two annihilation gammas (at rest or in flight).
            const AnnihilationPair pair =
                sample_annihilation_pair(pos_in_flight, pos_residual_KE, rng);
            const double ann_E1 = pair.energy_keV[0], ann_E2 = pair.energy_keV[1];
            const Eigen::Vector3d& ann_dir1 = pair.direction[0];
            const Eigen::Vector3d& ann_dir2 = pair.direction[1];

            // FEP early kill: if pair electron escape already blew the budget,
            // skip annihilation gamma transport entirely.
            if (fep_killed(ann_E1 + ann_E2)) {
                energy_keV = 0.0;
                goto done;
            }

            // Transport each annihilation gamma.
            // Keep the same config (fluorescence on, pair production off below threshold).
            // Each gamma is < 1022 keV at rest; an in-flight gamma can exceed it,
            // so pair production stays enabled by the recursive energy check.
            auto res1 = transport_photon(annih_pos, ann_dir1, ann_E1, geometry, config, rng);
            result.energy_deposited_scoring += res1.energy_deposited_scoring;
            result.energy_deposited_total   += res1.energy_deposited_total;
            if (res1.any_interaction_in_scoring)
                result.any_interaction_in_scoring = true;
            if (res1.escaped && res1.exit_energy_keV > 1.0)
                result.escaped_secondaries.push_back({
                    res1.exit_position, res1.exit_direction, res1.exit_energy_keV});
            result.escaped_secondaries.insert(
                result.escaped_secondaries.end(),
                res1.escaped_secondaries.begin(), res1.escaped_secondaries.end());
            if (config.record_interactions) {
                // Preserve the annihilation-photon history after the parent
                // pair-production record. Detector-response diagnostics use
                // this ordering to distinguish pair histories with zero, one,
                // or multiple subsequent Compton scatters. This is diagnostic
                // only and does not alter transport or the random stream.
                result.interactions.insert(
                    result.interactions.end(),
                    res1.interactions.begin(), res1.interactions.end());
                record_annihilation_photon(
                    result.pair_history.annihilation[0], ann_E1, res1);
            }

            // FEP early kill: first annihilation gamma escaped, check if second can save FEP.
            if (fep_killed(ann_E2)) {
                energy_keV = 0.0;
                goto done;
            }

            auto res2 = transport_photon(annih_pos, ann_dir2, ann_E2, geometry, config, rng);
            result.energy_deposited_scoring += res2.energy_deposited_scoring;
            result.energy_deposited_total   += res2.energy_deposited_total;
            if (res2.any_interaction_in_scoring)
                result.any_interaction_in_scoring = true;
            if (res2.escaped && res2.exit_energy_keV > 1.0)
                result.escaped_secondaries.push_back({
                    res2.exit_position, res2.exit_direction, res2.exit_energy_keV});
            result.escaped_secondaries.insert(
                result.escaped_secondaries.end(),
                res2.escaped_secondaries.begin(), res2.escaped_secondaries.end());
            if (config.record_interactions) {
                result.interactions.insert(
                    result.interactions.end(),
                    res2.interactions.begin(), res2.interactions.end());
                record_annihilation_photon(
                    result.pair_history.annihilation[1], ann_E2, res2);
            }

            energy_keV = 0.0;
            goto done;
        }

        } // switch

    } // for step

done:
    return result;
}

} // namespace ceelo
