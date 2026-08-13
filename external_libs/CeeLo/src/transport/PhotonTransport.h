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

/// @file PhotonTransport.h
/// @brief Photon transport loop inside detector materials.
///
/// Tracks a photon through the detector crystal (and dead layer), handling:
///   - Interaction distance sampling
///   - Interaction type selection (PE, Compton, Rayleigh, PP)
///   - Compton scattering (Klein-Nishina sampling)
///   - Photoelectric absorption (full deposit for now; fluorescence in Phase 3)
///   - Configurable N_max Compton scatter cutoff
///   - Boundary crossing between scoring/non-scoring regions

#include <Eigen/Core>
#include <random>
#include <cstdint>
#include <vector>

namespace ceelo {

class Material;
class Geometry;
struct PathSegment;

/// Configuration for the transport simulation.
struct TransportConfig {
    int max_compton_scatters = 40;  ///< N_max: force absorption after this many in-crystal Compton scatters.
                                    ///< Set high enough that truncation effectively never fires in the
                                    ///< validated regime (≤3 MeV); the chain terminates naturally as
                                    ///< energy degrades and PE (~E^-3) rises. Force-absorption sets the
                                    ///< TransportResult::forced_absorption flag (diagnostic only).
    double min_energy_keV = 1.0;    ///< Stop tracking below this energy (deposit remaining locally)
    bool enable_rayleigh = true;    ///< Include Rayleigh scattering (form-factor sampling)
    bool enable_pair_production = true; ///< Include pair production with 511 keV tracking
    bool enable_fluorescence = true;    ///< Track K-shell fluorescence X-rays explicitly
    int max_fluorescence_depth = 3;     ///< Max recursive fluorescence cascade depth
    int fluorescence_depth = 0;         ///< Current fluorescence recursion depth (internal)

    /// Doppler-broaden the Compton scattered-photon energy using the bound
    /// electron's momentum (impulse approximation, PENELOPE-style analytic
    /// subshell profiles).  Broadens the Compton edge to match GEANT4's
    /// G4LowEPComptonModel.  The subshell binding energy U_i is deposited
    /// locally (no Compton-vacancy relaxation), so total energy is conserved
    /// exactly and FEP/total efficiency are unaffected.  No fep_only_mode
    /// gating: FEP-only runs full physics inside the crystal and local-U_i
    /// deposition keeps the FEP indicator intact.
    bool enable_doppler_broadening = true;

    /// Enable CSDA electron tracking at interaction sites.
    ///
    /// When true, photoelectrons, Compton recoil electrons, and pair e⁺/e⁻ are
    /// tracked through the geometry using the Continuous Slowing Down Approximation.
    /// This accounts for electrons escaping the scoring volume (e.g. near the crystal
    /// surface), which causes FEP efficiency to be overestimated at high energies.
    ///
    /// Improvement over "deposit locally":
    ///   200 keV: ~2%  |  662 keV: ~7–9%  |  1.5 MeV: ~12–15%  |  3 MeV: ~20%
    ///
    /// Residual GEANT4 discrepancy with CSDA (shell corrections, bremsstrahlung,
    /// multiple scattering not yet modelled):
    ///   < 200 keV: ~1–3%  |  662 keV: ~3–6%  |  > 1 MeV: ~5–10%
    ///
    /// TODO (future): full electron transport via Molière multiple scattering,
    ///   discrete bremsstrahlung photon generation, delta-ray production.
    bool enable_electron_csda = true;

    // --- Diagnostic flags for GEANT4 comparison ---

    /// Disable Molière multiple scattering in CSDA electron transport.
    /// All electrons use straight-line CSDA regardless of energy.
    bool disable_moliere = false;

    /// Disable bremsstrahlung emission in CSDA electron transport.
    /// Electrons lose energy only via collision stopping.
    bool disable_brems = false;

    /// Record per-interaction data for diagnostic histograms.
    /// Zero cost when false. Used by diag_interactions tool.
    bool record_interactions = false;

    /// FEP-only mode: kill event when any energy escapes the scoring volume.
    /// Used by EfficiencyCalculator to short-circuit non-FEP events.
    bool fep_only_mode = false;

    /// Disable FEP early-kill optimization (for benchmarking only).
    /// When true, FEP mode still kills on crystal exit but doesn't check
    /// energy deficit mid-transport.
    bool disable_fep_early_kill = false;

    /// Force the photon's first interaction to occur inside the detector
    /// geometry (variance reduction). The pass-through branch is removed and
    /// the history carries weight = 1 - exp(-tau) in TransportResult::weight,
    /// where tau is the total optical depth along the entry ray. The first
    /// interaction point is sampled from the truncated exponential across
    /// segments; everything downstream is analog. Secondaries are never
    /// re-forced (the flag is cleared before any recursive transport).
    /// In fep_only_mode the optical depth stops at the first non-scoring
    /// segment after scoring (the photon would be killed there anyway).
    bool force_first_interaction = false;
};

/// Interaction types
enum class InteractionType : uint8_t {
    Photoelectric,
    Compton,
    Rayleigh,
    PairProduction
};

/// Per-interaction diagnostic record (only populated when config.record_interactions == true).
struct InteractionRecord {
    InteractionType type;
    Eigen::Vector3d position;
    double energy_before_keV;   ///< Photon energy before this interaction
    double energy_after_keV;    ///< Photon energy after (0 for PE, same for Rayleigh)
    double cos_theta;           ///< Cosine of scattering angle (1.0 for PE)
    bool in_scoring;            ///< Was this interaction in the scoring volume?
    int compton_scatter_index;  ///< 0-based index of Compton scatters (or # prior Comptons for PE)
};

/// Record of a fluorescence event during PE absorption.
/// Only populated when config.record_interactions == true.
struct FluorescenceRecord {
    int element_Z = 0;              ///< Element where PE occurred
    bool was_k_shell = false;       ///< K-shell (true) or outer-shell (false) PE
    bool fluorescence_emitted = false; ///< Was a K-fluorescence X-ray emitted?
    double fluorescence_energy_keV = 0.0; ///< Energy of emitted X-ray (0 if Auger)
    bool fluorescence_escaped = false;    ///< Did the X-ray escape the scoring volume?
    double fluorescence_deposit_in_scoring = 0.0; ///< Energy deposited by X-ray in scoring (keV)
};

/// Result of transporting a single photon through the detector.
struct TransportResult {
    double energy_deposited_scoring;   ///< Total energy deposited in scoring (active) volume (keV)
    double energy_deposited_total;     ///< Total energy deposited everywhere (keV)
    bool any_interaction_in_scoring;   ///< Did the photon interact at all in scoring volume?
    int num_compton_scatters;          ///< Number of Compton scatters in the crystal
    bool forced_absorption;            ///< Was the photon force-absorbed at N_max?

    /// History weight multiplier accumulated by biasing inside transport
    /// (currently only forced first interaction; always in (0, 1]).
    /// The caller must multiply this into the event weight before tallying.
    double weight = 1.0;

    // Exit state (always populated when photon escapes the geometry)
    bool escaped = false;                       ///< Did the photon exit the detector geometry?
    Eigen::Vector3d exit_position{0,0,0};       ///< Position at escape (world coords)
    Eigen::Vector3d exit_direction{0,0,0};      ///< Direction at escape (normalized)
    double exit_energy_keV = 0.0;               ///< Photon energy at escape (keV)

    /// Per-interaction records (only when config.record_interactions == true).
    std::vector<InteractionRecord> interactions;
    /// Energy of photon at escape (0 if absorbed). Only set when record_interactions == true.
    double escape_energy_keV = 0.0;

    /// Total electron kinetic energy that escaped the scoring volume (keV).
    /// Only set when config.record_interactions == true and enable_electron_csda == true.
    double electron_escape_keV = 0.0;

    /// Fluorescence records (only when config.record_interactions == true).
    std::vector<FluorescenceRecord> fluorescence_records;

    /// Secondary photon that escaped the detector geometry.
    /// Used by the Marinelli environmental bounce loop to recover escaped
    /// fluorescence X-rays, bremsstrahlung, and 511 keV annihilation gammas.
    struct EscapedPhoton {
        Eigen::Vector3d position;
        Eigen::Vector3d direction;
        double energy_keV;
    };

    /// All secondary photons that escaped the crystal during this transport.
    /// Populated at every recursive transport_photon() call site (fluorescence,
    /// bremsstrahlung, pair annihilation). Zero-cost when empty.
    std::vector<EscapedPhoton> escaped_secondaries;
};

/// Transport a single photon through the detector geometry.
///
/// The photon enters the detector geometry at a known position and direction.
/// We track it through the path segments (which may include dead layer,
/// active crystal, and possibly bore hole gaps) until it:
///   - Is photoelectrically absorbed
///   - Escapes the geometry
///   - Is force-absorbed after N_max Compton scatters
///   - Drops below the minimum energy threshold
///
/// @param position      Current photon position (world coords)
/// @param direction     Current photon direction (normalized)
/// @param energy_keV    Photon energy in keV
/// @param geometry      The detector geometry
/// @param config        Transport configuration
/// @param rng           Random number generator
/// @return              Transport result with energy deposition info
TransportResult transport_photon(
    Eigen::Vector3d position,
    Eigen::Vector3d direction,
    double energy_keV,
    const Geometry& geometry,
    const TransportConfig& config,
    std::mt19937_64& rng);

} // namespace ceelo
