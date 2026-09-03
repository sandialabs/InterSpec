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

/// @file ElectronCsda.h
/// @brief CSDA (Continuous Slowing Down Approximation) electron transport.
///
/// Implements CSDA electron tracking for photoelectrons, Compton recoil electrons,
/// and pair-production e+/e- pairs.  The CSDA approximation models electron energy
/// loss as continuous (no straggling), propagating the electron in a straight line
/// until its range is exhausted.  This significantly improves FEP efficiency accuracy
/// compared to "deposit all KE locally at the interaction point":
///
///   Energy regime     | CSDA improvement | Residual GEANT4 error
///   < 200 keV        | < 2%             | ~1–3%  (shell corrections, binding)
///   200–1000 keV     | 7–15%            | ~4–8%  (Bethe-Bloch ~8% under at 1 MeV)
///   1–3 MeV          | 15–25%           | ~5–10% (Bethe-Bloch ~3% under at 10 MeV)
///
/// Molière multiple scattering (added for KE > 200 keV):
///   Uses the Highland condensed-history formula for per-step RMS projected angle:
///     θ₀ = (13.6 MeV / βcp) × √(x/X₀) × (1 + 0.038 ln(x/X₀))
///   20 sub-steps per electron track.  Corrects FEP overestimate at high energy
///   by allowing electrons to scatter out of the scoring volume.
///
/// References: EGS4 (Nelson 1985), Geant4 EM Physics (Ivanchenko 2010),
///             Highland (1975) Nucl.Instrum.Meth. 129 497.

#include "geometry/Geometry.h"
#include "materials/Material.h"

#include <Eigen/Core>
#include <random>
#include <array>
#include <vector>

namespace ceelo {

/// A bremsstrahlung photon emitted during electron CSDA transport.
struct BremsPhoton {
    Eigen::Vector3d position;
    Eigen::Vector3d direction;
    double energy_keV;
};

/// Result of CSDA electron tracking through detector geometry.
struct ElectronDepositResult {
    double deposited_scoring_keV = 0.0;
    std::vector<BremsPhoton> brems_photons;

    /// World-frame position where the electron's condensed-history walk
    /// terminated — the CSDA stop point if it came to rest, or the last
    /// boundary crossing if it left all geometry with residual energy.
    /// Defaults to the start position (used when the walk produces no path,
    /// e.g. an empty ray trace).  Callers that annihilate a positron emit
    /// the back-to-back 511 keV pair from here rather than the creation
    /// (pair-production) vertex.
    Eigen::Vector3d stop_position{0.0, 0.0, 0.0};

    /// Positron-only: set true when the positron annihilated IN FLIGHT during
    /// the walk (sampled per substep against the Heitler annihilation rate over
    /// the collision stopping power).  When true, `stop_position` is the
    /// annihilation site and `residual_KE_at_annih_keV` is the positron kinetic
    /// energy carried into the photon pair: the caller must emit two photons
    /// summing to 2·mₑc² + residual_KE_at_annih_keV from `stop_position` and
    /// must NOT emit the at-rest back-to-back 511 keV pair.  When false the
    /// positron stopped (annihilate at rest at `stop_position`).
    bool annihilated_in_flight = false;
    double residual_KE_at_annih_keV = 0.0;
};

class SourceGeometry;

/// Result of an electron condensed-history walk through source geometry
/// (source material + shielding shells).
struct ElectronSourceWalkResult {
    /// Brems photons emitted along the walk (positions inside the source
    /// geometry; caller transports them through the remaining material).
    std::vector<BremsPhoton> brems_photons;

    /// True if the electron exited the source geometry before stopping.
    bool escaped = false;
    Eigen::Vector3d exit_position{0.0, 0.0, 0.0};
    Eigen::Vector3d exit_direction{0.0, 0.0, 1.0};
    double exit_KE_keV = 0.0;
};

/// Minimum electron kinetic energy (keV) for Moliere condensed-history transport
/// with discrete bremsstrahlung emission. Below this threshold, straight-line
/// CSDA is used with no angular deflection and no bremsstrahlung.
///
/// This value must not be lowered below kSB_min_energy_keV without extending
/// the Seltzer-Berger brems tables to cover the new range.
static constexpr double kMoliereBremsThreshold_keV = 200.0;

/// Minimum electron energy (keV) covered by the Seltzer-Berger brems tables.
/// If kMoliereBremsThreshold_keV is lowered below this, the SB data in
/// element_data.cpp must be regenerated with an extended energy range.
static constexpr double kSB_min_energy_keV = 200.0;

static_assert(kMoliereBremsThreshold_keV >= kSB_min_energy_keV,
    "Seltzer-Berger brems tables only cover electron energies >= 200 keV. "
    "To lower the Moliere/brems threshold, regenerate element_data.cpp with "
    "SB data extending to lower energies.");

/// Singleton that provides pre-computed CSDA range tables for Z = 1..92.
///
/// Stopping power: NIST ESTAR collision stopping power read directly from the
/// generated fixtures (`estar_stopping_data.h`), with the positron obtained
/// from the electron via the analytic Bhabha/Møller F⁺/F⁻ ratio; ranges by
/// trapezoid integration of 1/S over a 200-point log grid.  (The original
/// Bethe-Bloch/ICRU-37 implementation this doc-block used to describe was
/// replaced in the stopping-power migration.)
///
/// Usage:
///   double R = ElectronCsda::instance().range_gcm2(53, 662.0); // I at 662 keV
class ElectronCsda {
public:
    static const ElectronCsda& instance();

    /// The source-escape multiple-scattering / gate model is selected at COMPILE
    /// TIME via the CEELO_ESCAPE_MODEL macro (CMake cache var
    /// CEELO_SOURCE_ESCAPE_MODEL = gate|tails|gs); see walk_in_source_geometry()
    /// and the kSourceEscapeModel constant in ElectronCsda.cpp. gate (0, default,
    /// validated) uses the empirical exit-state albedo gate; tails (1) and gs (2)
    /// replace it with first-principles MSC (straggling + Rutherford tail, differing
    /// in the Gaussian-core width: Highland-minus-hard vs soft-transport-moment).

    /// CSDA range for pure element Z at kinetic energy KE_keV (g/cm²).
    double range_gcm2(int Z, double KE_keV) const;

    /// CSDA range for a compound material using the Bragg range rule (g/cm²).
    /// R_compound = 1 / (Σᵢ wᵢ / Rᵢ(E))
    /// Exact when stopping powers have the same energy shape (good to ~5–10%).
    /// is_positron selects the positron (Berger-Seltzer F⁺) range table.
    double range_gcm2_material(const Material& mat, double KE_keV,
                               bool is_positron = false) const;

    /// Find the residual kinetic energy (keV) of an electron that started with
    /// KE_keV and traversed path_gcm2 (g/cm²) of material mat.
    /// Returns 0.0 if the electron stops before traversing path_gcm2.
    /// is_positron selects the positron (Berger-Seltzer F⁺) range table.
    double residual_energy_keV(const Material& mat,
                               double KE_keV,
                               double path_gcm2,
                               bool is_positron = false) const;

    /// Energy deposited (keV) in scoring regions as an electron tracks through
    /// the geometry starting at `position` in `direction` with kinetic energy KE_keV.
    ///
    /// For KE ≤ 200 keV: straight-line CSDA (range too short for significant
    /// angular deflection).  No bremsstrahlung emitted.
    ///
    /// For KE > 200 keV: Molière condensed-history multiple scattering walk
    /// (20 sub-steps, Highland formula for per-step RMS angle).  This corrects
    /// the overestimated FEP at high energy due to electrons deflecting out of
    /// the scoring volume.  Discrete bremsstrahlung photons are emitted when
    /// radiative energy loss exceeds 10 keV per step, and returned for
    /// subsequent photon transport by the caller.
    ///
    /// @param position   Electron start position (cm, world frame)
    /// @param direction  Initial electron direction (normalized)
    /// @param KE_keV     Electron kinetic energy (keV)
    /// @param geometry   Detector geometry
    /// @param rng        Random number generator (consumed by Molière scatter)
    /// @param is_positron  When true the walk uses the positron Berger-Seltzer
    ///                   collision stopping term and samples in-flight
    ///                   annihilation per substep (see ElectronDepositResult's
    ///                   annihilated_in_flight / residual_KE_at_annih_keV).
    ///                   Default false reproduces the electron path bit-for-bit.
    /// @return           ElectronDepositResult with scoring deposit and brems photons
    ElectronDepositResult deposited_in_scoring(
                                const Eigen::Vector3d& position,
                                const Eigen::Vector3d& direction,
                                double KE_keV,
                                const Geometry& geometry,
                                std::mt19937_64& rng,
                                bool disable_moliere = false,
                                bool disable_brems = false,
                                bool is_positron = false) const;

    /// Sample bremsstrahlung photons from an electron (or positron) of kinetic
    /// energy KE_keV slowing down entirely within an infinite medium of `mat`.
    ///
    /// Geometry-free companion to deposited_in_scoring(): same condensed-history
    /// substepping, Seltzer-Berger emission sampling, and Highland direction
    /// evolution, but every photon is reported at the creation position
    /// (valid when the electron range is much smaller than both the medium
    /// thickness and the distance to the detector, e.g. source shields).
    /// Returns no photons for KE ≤ 200 keV (below the SB table floor).
    std::vector<BremsPhoton> sample_brems_in_material(
                                const Material& mat,
                                const Eigen::Vector3d& position,
                                const Eigen::Vector3d& direction,
                                double KE_keV,
                                std::mt19937_64& rng) const;

    /// Condensed-history walk of an electron (or positron) through source
    /// geometry (source material + shields), emitting bremsstrahlung along
    /// the way and detecting escape through the outer boundary.
    ///
    /// Same physics as deposited_in_scoring() — adaptive substeps, CSDA
    /// collision loss, Seltzer-Berger brems, Highland multiple scattering —
    /// but walks SourceGeometry::trace_source_segments() instead of detector
    /// ray tracing, and scores nothing: the outputs are the brems photons and
    /// the escape state. The Molière direction diffusion is what suppresses
    /// straight-line escape overestimates in dense shields (detour factor).
    /// For KE ≤ 200 keV the walk degenerates to straight-line CSDA with no
    /// brems (range < 0.1 mm in typical shields).
    /// The skin-escape suppression is now applied REGIME-AWARE at the exit
    /// boundary (source_escape_survival_exit), one material-general gate for
    /// electrons AND positrons — there is no per-caller skip flag (the soft
    /// low-Z / range-limited escape it used to over-suppress is now kept
    /// automatically via the exit-energy window).
    ElectronSourceWalkResult walk_in_source_geometry(
                                const SourceGeometry& source,
                                const Material& birth_material,
                                const Eigen::Vector3d& position,
                                const Eigen::Vector3d& direction,
                                double KE_keV,
                                std::mt19937_64& rng) const;

    /// ICRU Report 49 mean excitation energy I(Z) in eV, Z = 1..92.
    static double mean_excitation_eV(int Z);

    /// NIST ESTAR electron collision stopping power in MeV cm²/g for Z=1..92
    /// at kinetic energy KE_keV. With is_positron=true, applies the Bhabha/Møller
    /// collision-term ratio to the ESTAR electron value.
    static double stopping_power_MeV_cm2_g(int Z, double A_g_mol, double KE_keV,
                                           bool is_positron = false);

    /// Standard atomic weight A (g/mol) for Z = 1..92.
    static double atomic_weight(int Z);

    /// Tsai radiation length (g/cm²) for a single element.
    /// Formula: X₀ = 716.4 × A / (Z × (Z+1) × ln(287 / √Z))
    /// Accuracy: ~1–2% vs PDG tabulated values.
    static double radiation_length_gcm2_element(int Z, double A_g_mol);

    /// Compound radiation length via Bragg additivity (g/cm²).
    /// 1/X₀_compound = Σᵢ wᵢ / X₀ᵢ
    static double radiation_length_gcm2(const Material& mat);

    /// Mass-fraction-weighted effective atomic number of a material.
    /// Used as the Z lever for the skin-escape detour/transmission model.
    static double effective_Z(const Material& mat);

    /// Practical (extrapolated) electron range in g/cm² — the projected
    /// penetration depth R_ex = R_csda × detour(KE, Z_eff), with the detour
    /// ratio < 1 capturing multiple-scattering path-lengthening (smaller at
    /// lower energy and higher Z). Used by walk_in_source_geometry()'s
    /// skin-escape transmission gate; the detour parameterization carries the
    /// energy/Z dependence that separates sub-MeV from MeV escape.
    double extrapolated_range_gcm2_material(const Material& mat,
                                            double KE_keV) const;

    /// Regime-aware skin-escape SURVIVAL probability ∈ [0,1], evaluated at the
    /// EXIT boundary from the particle's exit state (post-walk), the product of:
    ///   • a surface-emergence (1 − backscatter/albedo) factor keyed on the EXIT
    ///     energy and EXIT-material Z — the albedo applies only to skin-DIFFUSION
    ///     escape (still-energetic exit), and vanishes for range-LIMITED exit
    ///     (E_exit → 0) and for the forward-peaked MeV channel, via
    ///     skin_escape_window(); and
    ///   • a depth gate T_path(τ), τ = (g/cm² path to the outer boundary)/R_ex,
    ///     that collapses near τ≈1 (Tabata-style) for electrons born deeper in.
    /// Replaces the legacy birth-energy source_escape_transmission(); see
    /// walk_in_source_geometry(). One material-general gate for electrons AND
    /// positrons — no per-caller skip flag.
    static double source_escape_survival_exit(double exit_KE_keV, double Z_exit,
                                              double tau);

private:
    ElectronCsda();  // Builds range tables on first access.

    static constexpr int    kNGrid   = 200;       // Log-spaced grid points
    static constexpr double kEMin_keV = 1.0;      // Minimum electron KE
    static constexpr double kEMax_keV = 20000.0;  // Maximum electron KE (20 MeV)

    // range_table_[Z-1][i] = CSDA range (g/cm²) at energy_grid_[i], for Z = 1..92.
    // range_table_pos_ is the positron counterpart (Berger-Seltzer F⁺ stopping);
    // ~1–3% longer range at MeV energies.  Built once in the constructor.
    std::array<std::array<double, kNGrid>, 92> range_table_;
    std::array<std::array<double, kNGrid>, 92> range_table_pos_;
    std::array<double, kNGrid> energy_grid_keV_;   // Kinetic energy grid (keV)
    std::array<double, kNGrid> log_energy_grid_;   // ln(energy_grid_keV_)

    // Precomputed per-element Seltzer-Berger brems normalization integrals over
    // kappa = k/KE, evaluated on energy_grid_keV_ for Z = 1..92:
    //   sb_Jchi_[Z-1][i]   = ∫[kappa_min,1] chi(Z, E_i, kappa) dkappa
    //   sb_Jchiok_[Z-1][i] = ∫[kappa_min,1] chi(Z, E_i, kappa)/kappa dkappa
    // with kappa_min = kBremsThreshold_keV / E_i.  Because the compound chi is a
    // fixed radiation-weighted linear combination of per-element chi, the compound
    // integral is the same linear combination of these per-element integrals.  This
    // lets brems_pemit_corr() replace the per-substep 64-point integral (which
    // dominated runtime via sb_chi) with a small weighted sum of table lookups.
    std::array<std::array<double, kNGrid>, 92> sb_Jchi_;
    std::array<std::array<double, kNGrid>, 92> sb_Jchiok_;

    double interpolate_range(int Z, double KE_keV) const;

    // Linear interpolation of a per-element table over the uniform log-energy grid.
    double interpolate_log_grid(const std::array<double, kNGrid>& tbl,
                                double KE_keV) const;

    // Compound (Bragg-additive) CSDA range table for a material, evaluated on the
    // energy grid and cached per-thread.  Replaces the per-call harmonic sum over
    // the composition; combined with table inversion in residual_energy_keV() this
    // removes the former 50-iteration bisection (the dominant cost after the brems
    // cache).  Cached thread-locally (worker threads are spawned per compute(), so
    // the cache lifetime is naturally bounded) and guarded by a composition check.
    const std::array<double, kNGrid>& compound_range_table(const Material& mat,
                                                           bool is_positron = false) const;

    // Compound ESTAR/Tsai radiative ratio for a material, evaluated
    // on the energy grid and cached per-thread (same scheme as compound_range_table).
    // Replaces the per-substep element loop with a single interpolation.
    const std::array<double, kNGrid>& compound_radcorr_table(const Material& mat) const;

public:
    /// Bremsstrahlung emission-rate correction factor <k>_1/k / <k>_SB for a
    /// compound material at electron energy KE_keV, using the precomputed
    /// Seltzer-Berger integral tables.  Equivalent to the per-substep 64-point
    /// numerical integration it replaces (exact at grid energies, linearly
    /// interpolated between).  Returns 1.0 below the brems threshold.
    double brems_pemit_corr(const Material& mat, double KE_keV) const;
};

/// Sample the direction of a photoelectron emitted from K-shell absorption
/// using the Sauter relativistic dipole formula:
///   dP/d(cos θ) ∝ (1 − cos²θ) / (1 − β cos θ)⁴
///
/// Algorithm (Sauter 1931, see Geant4 G4SauterGavrilaAngularDistribution):
///  1. Sample cos θ from the proposal p(u) ∝ 1/(1 − β u)²  (analytically invertible).
///  2. Accept/reject with probability (1 − β²)(1 − u²) / (1 − β u)².
///  3. Sample φ uniformly in [0, 2π] and rotate around the photon direction axis.
///
/// @param photon_direction  Incoming photon direction (normalized)
/// @param KE_keV            Photoelectron kinetic energy (keV)
/// @param rng               Random number generator
/// @return                  Electron direction (normalized)
Eigen::Vector3d sample_photoelectron_direction(
    const Eigen::Vector3d& photon_direction,
    double KE_keV,
    std::mt19937_64& rng);

/// Compute the Compton recoil electron direction from 4-momentum conservation.
///
///   p_e = p_γ_in − p_γ_out
///   direction ∝ E_in × d_in − E_out × d_out
///
/// This is exact (no approximation) and avoids needing the azimuthal angle φ.
///
/// @param photon_dir_in   Incoming photon direction (normalized)
/// @param photon_dir_out  Outgoing photon direction (normalized)
/// @param energy_in_keV   Incoming photon energy (keV)
/// @param energy_out_keV  Outgoing (scattered) photon energy (keV)
/// @return                Recoil electron direction (normalized)
Eigen::Vector3d compton_recoil_direction(
    const Eigen::Vector3d& photon_dir_in,
    const Eigen::Vector3d& photon_dir_out,
    double energy_in_keV,
    double energy_out_keV);

} // namespace ceelo
