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
} // anonymous namespace


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
    constexpr double kFepEscapeTol = 1.5; // keV — must match kFepTolerance in EfficiencyCalculator

    // Returns true if FEP is impossible: even if all remaining energy
    // (max_remaining_keV) deposits in scoring, the deficit exceeds tolerance.
    auto fep_killed = [&](double max_remaining_keV) -> bool {
        return config.fep_only_mode && !config.disable_fep_early_kill &&
               (initial_energy_keV - result.energy_deposited_scoring - max_remaining_keV) > kFepEscapeTol;
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

                if (fl && energy_keV > fl->k_edge_keV) {
                    // Probabilistic K-shell vs outer-shell PE selection.
                    // f_K = σ_K / σ_PE_total; outer shells deposit full energy locally.
                    const auto& xs_data = CrossSectionData::instance();
                    double sigma_K = xs_data.sigma_K_photoelectric(Z_pe, energy_MeV);
                    double sigma_tot = xs_data.sigma_photoelectric(Z_pe, energy_MeV);
                    double f_K = (sigma_tot > 0.0) ? std::min(sigma_K / sigma_tot, 1.0) : 1.0;

                    if (uniform(rng) >= f_K) {
                        // Outer-shell PE: full energy deposited locally (no K X-ray escape).
                        // deposited and pe_kinetic_energy remain = energy_keV (the defaults).
                        if (config.record_interactions) {
                            result.fluorescence_records.push_back({
                                Z_pe, false, false, 0.0, false, 0.0});
                        }
                        goto pe_deposit;
                    }

                    // K-shell absorption:
                    //   photoelectron KE = E_γ − E_K
                    pe_kinetic_energy = energy_keV - fl->k_edge_keV;
                    deposited = pe_kinetic_energy;

                    if (uniform(rng) < fl->fluorescence_yield) {
                        // --- K-shell fluorescence X-ray ---
                        // Select the transition line (Kα1, Kα2, Kβ, ...) by probability
                        double xi_line = uniform(rng), cum = 0.0;
                        float fl_energy = fl->line_energy_keV[0];
                        for (int li = 0; li < fl->num_lines; ++li) {
                            cum += fl->line_probability[li];
                            if (xi_line <= cum) {
                                fl_energy = fl->line_energy_keV[li];
                                break;
                            }
                        }

                        // L-shell Auger electrons from the K-vacancy fill:
                        // The transition (E_K → E_Kline) releases (E_K - E_Kline) as
                        // Auger electrons with sub-mm range → deposit locally.
                        auger_local = fl->k_edge_keV - static_cast<double>(fl_energy);
                        deposited += auger_local;

                        // Emit K-X-ray in an isotropic random direction
                        Eigen::Vector3d fl_dir = sample_isotropic_direction(rng);

                        // Allow secondary fluorescence cascade (e.g. Te Ka → Cd K-shell PE → Cd Ka)
                        // but cap recursion depth to prevent infinite loops.
                        TransportConfig fl_config = config;
                        fl_config.fluorescence_depth = config.fluorescence_depth + 1;
                        if (fl_config.fluorescence_depth >= config.max_fluorescence_depth)
                            fl_config.enable_fluorescence = false;

                        auto fl_res = transport_photon(
                            position, fl_dir, static_cast<double>(fl_energy),
                            geometry, fl_config, rng);

                        result.energy_deposited_scoring += fl_res.energy_deposited_scoring;
                        result.energy_deposited_total   += fl_res.energy_deposited_total;
                        if (fl_res.any_interaction_in_scoring)
                            result.any_interaction_in_scoring = true;

                        // Collect escaped fluorescence for environmental bounce
                        if (fl_res.escaped && fl_res.exit_energy_keV > 1.0)
                            result.escaped_secondaries.push_back({
                                fl_res.exit_position, fl_res.exit_direction,
                                fl_res.exit_energy_keV});
                        result.escaped_secondaries.insert(
                            result.escaped_secondaries.end(),
                            fl_res.escaped_secondaries.begin(),
                            fl_res.escaped_secondaries.end());

                        // Record fluorescence fate
                        if (config.record_interactions) {
                            bool escaped = (fl_res.energy_deposited_scoring <
                                            static_cast<double>(fl_energy) * 0.99);
                            result.fluorescence_records.push_back({
                                Z_pe, true, true,
                                static_cast<double>(fl_energy),
                                escaped, fl_res.energy_deposited_scoring});
                        }

                        // FEP early kill: fluorescence X-ray escaped scoring.
                        // Best case remaining = photoelectron KE + Auger (not yet deposited).
                        if (fep_killed(pe_kinetic_energy + auger_local)) {
                            energy_keV = 0.0;
                            goto done;
                        }

                    } else {
                        // --- Auger electrons ---
                        // Binding energy released locally via Auger electrons.
                        // Total deposited = (E_γ − E_K) + E_K = E_γ
                        auger_local = fl->k_edge_keV;
                        deposited += auger_local;

                        if (config.record_interactions) {
                            result.fluorescence_records.push_back({
                                Z_pe, true, false, 0.0, false, 0.0});
                        }
                    }
                } else {
                    // Below K-edge or no fluorescence data
                    if (config.record_interactions && fl) {
                        result.fluorescence_records.push_back({
                            Z_pe, false, false, 0.0, false, 0.0});
                    }
                }
                // If energy_keV <= k_edge_keV: only L-shell (or lower) absorption
                // possible — L X-rays have very short MFP, deposited locally.
                // deposited stays = energy_keV, pe_kinetic_energy = energy_keV (no change).
            }

            pe_deposit:
            // CSDA electron tracking: photoelectron KE only (Auger deposited locally).
            if (config.enable_electron_csda && pe_kinetic_energy > 0.0) {
                Eigen::Vector3d pe_dir = sample_photoelectron_direction(
                    direction, pe_kinetic_energy, rng);
                auto e_result = ElectronCsda::instance()
                    .deposited_in_scoring(position, pe_dir, pe_kinetic_energy, geometry, rng,
                                         config.disable_moliere, config.disable_brems);
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
                auto e_result = ElectronCsda::instance()
                    .deposited_in_scoring(position, e_dir,
                                         scatter.deposited_energy_keV, geometry, rng,
                                         config.disable_moliere, config.disable_brems);
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

            // Doppler: subshell binding energy is deposited locally (no
            // Compton-vacancy relaxation), keeping total deposit = E_in.
            if (scatter.binding_deposit_keV > 0.0) {
                if (in_scoring) {
                    result.energy_deposited_scoring += scatter.binding_deposit_keV;
                    result.any_interaction_in_scoring = true;
                }
                result.energy_deposited_total += scatter.binding_deposit_keV;
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
                auto e_neg = ElectronCsda::instance()
                    .deposited_in_scoring(position, dir_neg, ke_each, geometry, rng,
                                         config.disable_moliere, config.disable_brems);
                auto e_pos = ElectronCsda::instance()
                    .deposited_in_scoring(position, dir_pos, ke_each, geometry, rng,
                                         config.disable_moliere, config.disable_brems,
                                         /*is_positron=*/true);
                // The positron annihilates where it comes to rest (or in flight),
                // not at the creation vertex.  (The e⁻ does not annihilate.)
                annih_pos       = e_pos.stop_position;
                pos_in_flight   = e_pos.annihilated_in_flight;
                pos_residual_KE = e_pos.residual_KE_at_annih_keV;
                double dep_pair = e_neg.deposited_scoring_keV + e_pos.deposited_scoring_keV;
                result.energy_deposited_scoring += dep_pair;
                if (dep_pair > 0.0) result.any_interaction_in_scoring = true;

                // Track electron escape.  In-flight residual KE left as photon
                // energy, not slowing-down deposit, so exclude it from "escape".
                if (config.record_interactions) {
                    double slowing = deposited_local - (pos_in_flight ? pos_residual_KE : 0.0);
                    double e_escape = slowing - dep_pair;
                    if (e_escape > 0.01) result.electron_escape_keV += e_escape;
                }

                // Transport bremsstrahlung photons from both e- and e+
                for (const auto* e_res : {&e_neg, &e_pos}) {
                    for (const auto& bp : e_res->brems_photons) {
                        auto br = transport_photon(bp.position, bp.direction, bp.energy_keV,
                                                   geometry, config, rng);
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

            // Emit the two annihilation gammas.
            //  - At rest: two back-to-back 511 keV photons (isotropic axis).
            //  - In flight: total energy 2·mₑc² + residual KE, split by a sampled
            //    fraction and emitted in two independent isotropic directions
            //    (continuum, not a clean 511 pair).  Same approximation as the
            //    cascade β⁺ path; proper 2-body boost kinematics are not modeled.
            constexpr double kElectronMassKeV = 510.998950;
            double ann_E1, ann_E2;
            Eigen::Vector3d ann_dir1, ann_dir2;
            if (pos_in_flight) {
                double Etot = 2.0 * kElectronMassKeV + pos_residual_KE;
                double x    = 0.5 + (uniform(rng) - 0.5) * (pos_residual_KE / Etot);
                ann_E1   = x * Etot;
                ann_E2   = Etot - ann_E1;
                ann_dir1 = sample_isotropic_direction(rng);
                ann_dir2 = sample_isotropic_direction(rng);
            } else {
                ann_E1   = kElectronMassKeV;
                ann_E2   = kElectronMassKeV;
                ann_dir1 = sample_isotropic_direction(rng);
                ann_dir2 = -ann_dir1;
            }

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

            energy_keV = 0.0;
            goto done;
        }

        } // switch

    } // for step

done:
    return result;
}

} // namespace ceelo
