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

/// @file CrossSectionData.h
/// @brief Photon cross-section data access for all elements Z=1..92.
///
/// Cross-section data is compiled directly into the library as constant arrays.
/// Photon curves and angular factors come directly from EPICS2023 EPDL; K/L
/// relaxation comes directly from EPICS2023 EADL; bremsstrahlung comes from
/// public-domain NIST EPQ; and Compton profile support comes from xraylib 4.2.1.
/// Exact source locks and transformations are in tools/prepare_cross_sections/.
/// All data is immutable and thus thread-safe.
///
/// Storage format: log10(energy/MeV) vs log10(cross_section/barn), using
/// log-log linear interpolation.

#include <cstdint>
#include <cmath>
#include <cassert>
#include <random>

namespace ceelo {

struct PhotonProcessCurve;
struct PhotonRayleighValues;
struct PhotonAngularCurve;

/// Maximum atomic number we support data for.
static constexpr int kMaxZ = 92;

/// Seltzer-Berger bremsstrahlung spectral shape data.
/// chi(Z, T, kappa) = (beta^2 / Z^2) * kappa * dSigma/dkappa
/// Used for rejection sampling of brems photon energy.
extern const uint16_t kSB_n_kappa;       ///< Number of k/T fraction grid points (32)
extern const uint16_t kSB_n_energy;      ///< Number of electron energy grid points (27)
extern const float kSB_kappa[];          ///< k/T fraction values, ascending [kSB_n_kappa]
extern const float kSB_log_E_keV[];      ///< log10(electron KE / keV), ascending [kSB_n_energy]
extern const float kSB_chi_scale[kMaxZ]; ///< Per-element uint16 decode scales

/// Fluorescence line data for a single element.
/// Only K-shell fluorescence is tracked (L-shell and higher are deposited locally).
struct FluorescenceData {
    float k_edge_keV;           ///< K-shell binding energy (keV)
    float fluorescence_yield;   ///< omega_K: probability of X-ray vs. Auger

    // Radiative transition probabilities and energies (Kalpha1, Kalpha2, Kbeta)
    // Normalized so probabilities sum to 1 among the radiative transitions.
    uint8_t num_lines;
    const float* line_energy_keV;     ///< Transition energies (keV), num_lines entries
    const float* line_probability;    ///< Relative probabilities, num_lines entries
};

/// Radiative data for one L subshell (L1, L2 or L3).
struct LSubshellFluor {
    float fluorescence_yield;   ///< omega_Li: probability of an x-ray vs. Auger
    uint8_t num_lines;
    const float* line_energy_keV;  ///< L line energies (keV), num_lines entries
    const float* line_probability; ///< relative probabilities, num_lines entries
};

/// Per-subshell L-shell fluorescence data for a single element. The cascade
/// vacancy model resolves which L subshell a vacancy lands in (from the
/// per-transition IC subshell coefficients), then emits from that subshell's
/// line set. Only elements whose L lines clear the 10 keV x-ray cut have
/// populated subshells; lighter elements have all-zero (empty) subshells.
/// Generated separately from element_data.cpp (see relaxation_epics_data.cpp).
/// Relaxation data extends through Z=99 for radioactive-decay daughters even
/// though photon/electron transport tables stop at Z=92.
struct LFluorescenceData {
    float l3_edge_keV;          ///< L3-shell binding energy (keV)
    /// Coster-Kronig L-vacancy transfer yields {f12, f13, f23}: probability an
    /// L1 vacancy transfers (radiationlessly) to L2 / L3, and an L2 vacancy to
    /// L3, before radiating. Applied before line selection so emission shifts
    /// toward the higher subshells (notably L3 -> Lalpha) as in a full atomic
    /// deexcitation cascade. Values come directly from EPICS2023 EADL.
    float coster_kronig[3];
    LSubshellFluor sub[3];      ///< L1, L2, L3
};

/// Per-element support data that is not part of the compact photon tables.
/// Photon cross sections and angular factors are queried through CrossSectionData.
struct ElementData {
    uint8_t Z;

    /// Quantized Seltzer-Berger table. Decode coefficient i as
    /// sb_chi_quantized[i] * kSB_chi_scale[Z-1].
    const uint16_t* sb_chi_quantized;

    // Subshell data for Compton Doppler broadening (impulse approximation with
    // PENELOPE-style analytic one-parameter profiles).  Per occupied subshell:
    // occupancy (electrons), binding energy U_i (keV), and the per-electron
    // Compton profile at zero momentum J_i(0) (1/a.u., Biggs/EPDL via xraylib).
    // These fields are last so a stale generated file value-initializes them
    // (num_compton_shells = 0 disables Doppler for that element).
    uint8_t num_compton_shells;
    const float* shell_occupancy;     ///< electrons per subshell
    const float* shell_binding_keV;   ///< subshell binding energy U_i (keV)
    const float* shell_J0;            ///< per-electron J_i(p=0), 1/a.u.
};

/// Singleton-like accessor for all cross-section data.
/// All source data is immutable and compiled into the binary.
class CrossSectionData {
public:
    /// Get the singleton instance.
    static const CrossSectionData& instance();

    /// Get element data for atomic number Z (1-based).
    /// @pre Z >= 1 && Z <= kMaxZ
    const ElementData& element(int Z) const;

    /// Interpolate a single cross-section type at energy E (MeV) for element Z.
    /// Returns the cross-section in barns.
    /// Uses log-log linear interpolation.
    double sigma_photoelectric(int Z, double energy_MeV) const;
    double sigma_K_photoelectric(int Z, double energy_MeV) const;
    double sigma_compton(int Z, double energy_MeV) const;
    double sigma_rayleigh(int Z, double energy_MeV) const;
    double sigma_pair_production(int Z, double energy_MeV) const;

    /// Get all four partial cross-sections at once (more efficient — single binary search).
    /// Results in barns.
    struct PartialCrossSections {
        double sigma_pe;
        double sigma_cs;
        double sigma_rs;
        double sigma_pp;
        double sigma_total() const { return sigma_pe + sigma_cs + sigma_rs + sigma_pp; }
    };
    PartialCrossSections all_cross_sections(int Z, double energy_MeV) const;

    /// Evaluate the incoherent scattering function S(x, Z).
    /// x = momentum transfer parameter in inverse angstroms.
    /// x = (E/hc) * sin(theta/2) / (2*pi), but caller typically computes this.
    double scattering_function_S(int Z, double x) const;

    /// Sample a coherent (Rayleigh) momentum transfer x in [0, x_max] from the
    /// distribution p(x) ∝ F(x,Z)^2 (i.e. the form-factor part of the Rayleigh
    /// differential cross-section), by inverse transform of the precomputed
    /// cumulative ∫ F^2 d(x^2) table.  The caller converts x to cos(theta) and
    /// applies the (1+cos^2 theta)/2 polarization rejection.
    /// x is the momentum transfer in inverse angstroms; x_max = (E/hc) for the
    /// backscatter limit.
    double sample_rayleigh_x(int Z, double x_max, std::mt19937_64& rng) const;

    /// Get fluorescence data for element Z.
    /// Returns nullptr if the element has no K-shell fluorescence data.
    const FluorescenceData* fluorescence(int Z) const;

    /// Get effective L-shell fluorescence data for element Z (for the cascade
    /// vacancy model). Returns nullptr if no L fluorescence data is tabulated.
    const LFluorescenceData* l_fluorescence(int Z) const;

    /// Interpolate the Seltzer-Berger chi function for element Z at electron
    /// kinetic energy T_keV and reduced photon energy kappa = k/T.
    /// Uses bilinear interpolation in (log10(T_keV), kappa).
    double sb_chi(int Z, double T_keV, double kappa) const;

    /// Hoisted Seltzer-Berger chi evaluation for compound materials.
    /// The (T_keV, kappa) -> grid-bracket mapping uses the GLOBAL kSB_log_E_keV /
    /// kSB_kappa grids and is Z-INDEPENDENT, so a brems rejection loop can compute
    /// the energy bracket once per electron energy and the kappa bracket once per
    /// kappa, then loop only the per-element table fetch.  sb_chi_bracketed()
    /// reproduces sb_chi()'s exact kappa-first bilinear, so for the matching
    /// brackets the result is bit-for-bit identical to sb_chi(Z, T_keV, kappa).
    struct SBChiEBracket { int iE; double tE; };   ///< energy bracket (per emission)
    struct SBChiKBracket { int iK; double tK; };   ///< kappa bracket (per attempt)
    SBChiEBracket sb_chi_energy_bracket(double T_keV) const;
    SBChiKBracket sb_chi_kappa_bracket(double kappa) const;
    double sb_chi_bracketed(int Z, const SBChiEBracket& eb, const SBChiKBracket& kb) const;

    /// Atomic weight in g/mol for element Z.
    double atomic_weight(int Z) const;

private:
    CrossSectionData() = default;

    // Interpolation helpers: binary search + log-log linear interpolation.
    // log_x is log10 of the query point.
    static double interpolate_log_log(const PhotonProcessCurve& curve, double log_x);

    static double interpolate_rayleigh(int Z,
                                       const PhotonRayleighValues& values,
                                       double log_x);

    // Linear interpolation helper for scattering functions (not log-log in y).
    static double interpolate_log_lin(const PhotonAngularCurve& curve, double log_x);

};

} // namespace ceelo
