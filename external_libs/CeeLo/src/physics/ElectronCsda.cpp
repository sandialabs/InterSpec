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

#include "physics/ElectronCsda.h"
#include "physics/estar_stopping_data.h"
#include "cross_sections/CrossSectionData.h"
#include "geometry/Geometry.h"
#include "geometry/SourceGeometry.h"

#include <Eigen/Geometry>

#include <cmath>
#include <cassert>
#include <algorithm>
#include <stdexcept>
#include <random>
#include <cstdlib>

namespace ceelo {

namespace {

// ============================================================
// Physical constants
// ============================================================
constexpr double kMeC2_keV = 510.9989461;   // Electron rest energy (keV)
constexpr double kMeC2_MeV = 0.5109989461;  // Electron rest energy (MeV)
constexpr double kPi   = 3.14159265358979323846;
constexpr double kTwoPi = 2.0 * kPi;

// ============================================================
// Positron in-flight annihilation (Heitler).
//
// Total e⁺e⁻ annihilation cross-section PER ELECTRON (cm²) for a positron of
// kinetic energy T_keV — Heitler / G4eplusAnnihilation form.  (A parallel copy
// lives in EfficiencyCalculator.cpp for the cascade β⁺ source-side path; the two
// should be unified into one physics home.)  Used to sample whether a pair-
// production positron annihilates in flight before stopping in the crystal.
// ============================================================
double positron_annih_xsec_per_electron(double T_keV) {
    const double gamma = 1.0 + T_keV / kMeC2_keV;
    if (gamma <= 1.0 + 1.0e-6) return 0.0;
    const double g2 = gamma * gamma;
    const double s = std::sqrt(g2 - 1.0);
    constexpr double pi_re2 = 2.4946e-25;  // π·rₑ² (cm²)
    return (pi_re2 / (gamma + 1.0)) *
           (((g2 + 4.0 * gamma + 1.0) / (g2 - 1.0)) * std::log(gamma + s)
            - (gamma + 3.0) / s);
}

// Electron number density (1/cm³) of a material.
double electron_number_density(const Material& mat) {
    constexpr double NA = 6.02214076e23;
    const double rho = mat.density();
    double ne = 0.0;
    for (const auto& c : mat.composition())
        ne += rho * NA * c.mass_fraction * c.Z / ElectronCsda::atomic_weight(c.Z);
    return ne;
}

// ============================================================
// Hoisted compound Seltzer-Berger chi for brems rejection sampling.
//
// The rejection loop calls compound_chi(kappa) ~100x per emitted photon at a
// CONSTANT electron energy KE.  The SB grids (kSB_log_E_keV / kSB_kappa) are
// global, so the energy bracket and the per-element radiation weights
// w_rad = w_mass·Z(Z+1)/A are Z-/kappa-constant and precomputed once; each call
// then adds only the kappa bracket + per-element table fetch.  Bit-identical to
// the former per-call lambda: same w_rad values, same element accumulation order,
// same w_sum, and CrossSectionData::sb_chi_bracketed reproduces sb_chi's exact
// kappa-first bilinear.
// ============================================================
class CompoundBremsChi {
public:
    CompoundBremsChi(const CrossSectionData& xs, const Material& mat, double KE_keV)
        : xs_(xs)
        , eb_(xs.sb_chi_energy_bracket(KE_keV))
    {
        const auto& comps = mat.composition();
        ncomp_ = comps.size();
        assert(ncomp_ <= kMaxComp);
        for (size_t i = 0; i < ncomp_; ++i) {
            const int Z = comps[i].Z;
            z_[i]     = Z;
            w_rad_[i] = comps[i].mass_fraction * Z * (Z + 1.0)
                        / ElectronCsda::atomic_weight(Z);
            w_sum_   += w_rad_[i];
        }
    }

    double operator()(double kappa) const {
        const auto kb = xs_.sb_chi_kappa_bracket(kappa);
        double chi_val = 0.0;
        for (size_t i = 0; i < ncomp_; ++i)
            chi_val += w_rad_[i] * xs_.sb_chi_bracketed(z_[i], eb_, kb);
        return (w_sum_ > 0.0) ? chi_val / w_sum_ : 1.0;
    }

private:
    static constexpr size_t kMaxComp = 32;
    const CrossSectionData& xs_;
    CrossSectionData::SBChiEBracket eb_;
    std::array<int, kMaxComp>    z_{};
    std::array<double, kMaxComp> w_rad_{};
    size_t ncomp_ = 0;
    double w_sum_ = 0.0;
};

// ============================================================
// ICRU Report 49 (1993) Table 5.1: Mean excitation energies I(Z) in eV.
// Used in the Bethe-Bloch collision stopping power formula.
// Source: ICRU 49 / NIST PSTAR database.
// ============================================================
static const double kICRU49_I_eV[92] = {
    /* Z=1  H  */   19.2,
    /* Z=2  He */   41.8,
    /* Z=3  Li */   40.0,
    /* Z=4  Be */   63.7,
    /* Z=5  B  */   76.0,
    /* Z=6  C  */   81.0,
    /* Z=7  N  */   82.0,
    /* Z=8  O  */   95.0,
    /* Z=9  F  */  115.0,
    /* Z=10 Ne */  137.0,
    /* Z=11 Na */  149.0,
    /* Z=12 Mg */  156.0,
    /* Z=13 Al */  166.0,
    /* Z=14 Si */  173.0,
    /* Z=15 P  */  173.0,
    /* Z=16 S  */  180.0,
    /* Z=17 Cl */  174.0,
    /* Z=18 Ar */  188.0,
    /* Z=19 K  */  190.0,
    /* Z=20 Ca */  191.0,
    /* Z=21 Sc */  216.0,
    /* Z=22 Ti */  233.0,
    /* Z=23 V  */  245.0,
    /* Z=24 Cr */  257.0,
    /* Z=25 Mn */  272.0,
    /* Z=26 Fe */  286.0,
    /* Z=27 Co */  297.0,
    /* Z=28 Ni */  311.0,
    /* Z=29 Cu */  322.0,
    /* Z=30 Zn */  330.0,
    /* Z=31 Ga */  334.0,
    /* Z=32 Ge */  350.0,
    /* Z=33 As */  347.0,
    /* Z=34 Se */  348.0,
    /* Z=35 Br */  357.0,
    /* Z=36 Kr */  352.0,
    /* Z=37 Rb */  363.0,
    /* Z=38 Sr */  366.0,
    /* Z=39 Y  */  379.0,
    /* Z=40 Zr */  393.0,
    /* Z=41 Nb */  417.0,
    /* Z=42 Mo */  424.0,
    /* Z=43 Tc */  428.0,
    /* Z=44 Ru */  441.0,
    /* Z=45 Rh */  449.0,
    /* Z=46 Pd */  470.0,
    /* Z=47 Ag */  470.0,
    /* Z=48 Cd */  469.0,
    /* Z=49 In */  488.0,
    /* Z=50 Sn */  488.0,
    /* Z=51 Sb */  487.0,
    /* Z=52 Te */  485.0,
    /* Z=53 I  */  491.0,
    /* Z=54 Xe */  482.0,
    /* Z=55 Cs */  488.0,
    /* Z=56 Ba */  491.0,
    /* Z=57 La */  501.0,
    /* Z=58 Ce */  523.0,
    /* Z=59 Pr */  535.0,
    /* Z=60 Nd */  546.0,
    /* Z=61 Pm */  560.0,
    /* Z=62 Sm */  574.0,
    /* Z=63 Eu */  580.0,
    /* Z=64 Gd */  591.0,
    /* Z=65 Tb */  614.0,
    /* Z=66 Dy */  628.0,
    /* Z=67 Ho */  650.0,
    /* Z=68 Er */  658.0,
    /* Z=69 Tm */  674.0,
    /* Z=70 Yb */  684.0,
    /* Z=71 Lu */  694.0,
    /* Z=72 Hf */  705.0,
    /* Z=73 Ta */  718.0,
    /* Z=74 W  */  727.0,
    /* Z=75 Re */  736.0,
    /* Z=76 Os */  746.0,
    /* Z=77 Ir */  757.0,
    /* Z=78 Pt */  790.0,
    /* Z=79 Au */  790.0,
    /* Z=80 Hg */  800.0,
    /* Z=81 Tl */  810.0,
    /* Z=82 Pb */  823.0,
    /* Z=83 Bi */  823.0,
    /* Z=84 Po */  830.0,
    /* Z=85 At */  825.0,
    /* Z=86 Rn */  794.0,
    /* Z=87 Fr */  827.0,
    /* Z=88 Ra */  826.0,
    /* Z=89 Ac */  841.0,
    /* Z=90 Th */  847.0,
    /* Z=91 Pa */  878.0,
    /* Z=92 U  */  890.0,
};

// Standard atomic weights (g/mol) for Z = 1..92.
// Values from IUPAC 2021 Table of Atomic Weights (rounded).
static const double kAtomicWeight[92] = {
    /* Z=1  H  */   1.008,
    /* Z=2  He */   4.003,
    /* Z=3  Li */   6.941,
    /* Z=4  Be */   9.012,
    /* Z=5  B  */  10.811,
    /* Z=6  C  */  12.011,
    /* Z=7  N  */  14.007,
    /* Z=8  O  */  15.999,
    /* Z=9  F  */  18.998,
    /* Z=10 Ne */  20.180,
    /* Z=11 Na */  22.990,
    /* Z=12 Mg */  24.305,
    /* Z=13 Al */  26.982,
    /* Z=14 Si */  28.086,
    /* Z=15 P  */  30.974,
    /* Z=16 S  */  32.065,
    /* Z=17 Cl */  35.453,
    /* Z=18 Ar */  39.948,
    /* Z=19 K  */  39.098,
    /* Z=20 Ca */  40.078,
    /* Z=21 Sc */  44.956,
    /* Z=22 Ti */  47.867,
    /* Z=23 V  */  50.942,
    /* Z=24 Cr */  51.996,
    /* Z=25 Mn */  54.938,
    /* Z=26 Fe */  55.845,
    /* Z=27 Co */  58.933,
    /* Z=28 Ni */  58.693,
    /* Z=29 Cu */  63.546,
    /* Z=30 Zn */  65.380,
    /* Z=31 Ga */  69.723,
    /* Z=32 Ge */  72.630,
    /* Z=33 As */  74.922,
    /* Z=34 Se */  78.971,
    /* Z=35 Br */  79.904,
    /* Z=36 Kr */  83.798,
    /* Z=37 Rb */  85.468,
    /* Z=38 Sr */  87.620,
    /* Z=39 Y  */  88.906,
    /* Z=40 Zr */  91.224,
    /* Z=41 Nb */  92.906,
    /* Z=42 Mo */  95.950,
    /* Z=43 Tc */  98.000,
    /* Z=44 Ru */ 101.070,
    /* Z=45 Rh */ 102.906,
    /* Z=46 Pd */ 106.420,
    /* Z=47 Ag */ 107.868,
    /* Z=48 Cd */ 112.414,
    /* Z=49 In */ 114.818,
    /* Z=50 Sn */ 118.710,
    /* Z=51 Sb */ 121.760,
    /* Z=52 Te */ 127.600,
    /* Z=53 I  */ 126.904,
    /* Z=54 Xe */ 131.293,
    /* Z=55 Cs */ 132.905,
    /* Z=56 Ba */ 137.327,
    /* Z=57 La */ 138.905,
    /* Z=58 Ce */ 140.116,
    /* Z=59 Pr */ 140.908,
    /* Z=60 Nd */ 144.242,
    /* Z=61 Pm */ 145.000,
    /* Z=62 Sm */ 150.360,
    /* Z=63 Eu */ 151.964,
    /* Z=64 Gd */ 157.250,
    /* Z=65 Tb */ 158.925,
    /* Z=66 Dy */ 162.500,
    /* Z=67 Ho */ 164.930,
    /* Z=68 Er */ 167.259,
    /* Z=69 Tm */ 168.934,
    /* Z=70 Yb */ 173.045,
    /* Z=71 Lu */ 174.967,
    /* Z=72 Hf */ 178.490,
    /* Z=73 Ta */ 180.948,
    /* Z=74 W  */ 183.840,
    /* Z=75 Re */ 186.207,
    /* Z=76 Os */ 190.230,
    /* Z=77 Ir */ 192.217,
    /* Z=78 Pt */ 195.085,
    /* Z=79 Au */ 196.967,
    /* Z=80 Hg */ 200.592,
    /* Z=81 Tl */ 204.383,
    /* Z=82 Pb */ 207.200,
    /* Z=83 Bi */ 208.980,
    /* Z=84 Po */ 209.000,
    /* Z=85 At */ 210.000,
    /* Z=86 Rn */ 222.000,
    /* Z=87 Fr */ 223.000,
    /* Z=88 Ra */ 226.000,
    /* Z=89 Ac */ 227.000,
    /* Z=90 Th */ 232.038,
    /* Z=91 Pa */ 231.036,
    /* Z=92 U  */ 238.029,
};

/// Build a perpendicular basis for direction d (normalized).
/// Returns two vectors u, v such that (u, v, d) is orthonormal.
void build_perp_basis(const Eigen::Vector3d& d,
                      Eigen::Vector3d& u,
                      Eigen::Vector3d& v)
{
    // Pick the axis least aligned with d
    if (std::abs(d.x()) <= std::abs(d.y()) && std::abs(d.x()) <= std::abs(d.z())) {
        u = Eigen::Vector3d(0.0, -d.z(), d.y()).normalized();
    } else if (std::abs(d.y()) <= std::abs(d.z())) {
        u = Eigen::Vector3d(-d.z(), 0.0, d.x()).normalized();
    } else {
        u = Eigen::Vector3d(-d.y(), d.x(), 0.0).normalized();
    }
    v = d.cross(u);  // Already normalized since d and u are orthonormal
}

// Bohr energy-loss straggling constant: ξ_B = 4π r_e² (mₑc²)² N_A = 0.1569 MeV² cm²/g,
// giving the Gaussian energy-loss variance σ² = ξ_B × (Z/A) × Δx × (1 - β²/2)   [MeV²].
// NOTE: reserved — NOT currently used.  Energy loss in the condensed-history walk below is
// the deterministic CSDA value; no stochastic energy-loss straggling is sampled yet.
// TODO: try implementing per-step Gaussian energy-loss straggling with this constant
//   (one gauss() draw per step: ΔKE += N(0, σ); clamp ΔKE ≥ 0).  Expected to improve the
//   FEP-edge / escape-peak shape for thin detectors (CZT e⁻ escape > 1 MeV) and the
//   source-electron skin-escape energy spectrum — straggling lets near-threshold electrons
//   sometimes stop/escape against the deterministic prediction, which should reduce the
//   empirical skin-escape gate's calibration load (see source_escape_transmission / kAlbScale).
[[maybe_unused]] static constexpr double kBohrXi = 0.1569;  // MeV² cm²/g

/// Radiative stopping power correction factor for element Z at KE_keV.
/// Multiply KE_MeV/X0 by this to get S_rad consistent with NIST ESTAR.
static double radiative_correction_factor(int Z, double KE_keV)
{
    if (Z < 1 || Z > 92 || KE_keV <= 0.0) return 0.0;
    const double stopping = estar_data::radiative_stopping(Z, KE_keV);
    const double energy_MeV = KE_keV * 1.0e-3;
    const double A = ElectronCsda::atomic_weight(Z);
    const double X0 = ElectronCsda::radiation_length_gcm2_element(Z, A);
    return (energy_MeV > 0.0 && X0 > 0.0) ? stopping * X0 / energy_MeV : 0.0;
}

/// Compound radiative correction: mass-fraction-weighted average of per-element corrections.
static double radiative_correction_compound(const Material& mat, double KE_keV)
{
    // Weight by each element's contribution to S_rad (∝ w_i × Z_i(Z_i+1)/A_i)
    // For simplicity, use mass-fraction-weighted average of per-element X0-normalized
    // corrections (exact if all elements have similar S_rad shape).
    double sum_wf = 0.0;
    double sum_corr = 0.0;
    for (const auto& comp : mat.composition()) {
        double A  = ElectronCsda::atomic_weight(comp.Z);
        double X0_i = ElectronCsda::radiation_length_gcm2_element(comp.Z, A);
        double weight = comp.mass_fraction / X0_i;  // contribution to 1/X0
        double corr_i = radiative_correction_factor(comp.Z, KE_keV);
        sum_wf   += weight;
        sum_corr += weight * corr_i;
    }
    return (sum_wf > 0.0) ? sum_corr / sum_wf : 0.7;
}

} // anonymous namespace


// ============================================================
// ElectronCsda static accessors
// ============================================================

double ElectronCsda::mean_excitation_eV(int Z)
{
    if (Z < 1 || Z > 92) return 200.0; // fallback
    return kICRU49_I_eV[Z - 1];
}

double ElectronCsda::atomic_weight(int Z)
{
    if (Z < 1 || Z > 92) return static_cast<double>(Z) * 2.0; // rough fallback
    return kAtomicWeight[Z - 1];
}

/// Direct NIST ESTAR collision stopping power for electrons.
///
/// @param is_positron  When true, use the Bhabha F⁺(τ) collision term (positrons)
///   instead of the Møller F⁻(τ) (electrons). ESTAR tabulates electrons, so the
///   direct electron stopping power is multiplied by the analytic Bhabha/Møller
///   collision-term ratio for positron transport.
double ElectronCsda::stopping_power_MeV_cm2_g(int Z, double A_g_mol, double KE_keV,
                                              bool is_positron)
{
    if (KE_keV <= 0.0 || A_g_mol <= 0.0 || Z < 1 || Z > 92) return 1e30;

    // Direct all-element ESTAR collision values replace the former 12-element,
    // 15-energy correction surface. ESTAR's supported migration range begins at
    // 10 keV; below it, continue the first-interval log-log slope to keep the
    // existing 1 keV CSDA integration grid continuous.
    auto electron_estar = [Z](double energy_keV) {
        if (energy_keV >= 10.0)
            return estar_data::collision_stopping(Z, energy_keV);
        const double s10 = estar_data::collision_stopping(Z, 10.0);
        const double s12 = estar_data::collision_stopping(Z, 12.5);
        const double slope = std::log(s12 / s10) / std::log(12.5 / 10.0);
        return s10 * std::pow(energy_keV / 10.0, slope);
    };

    const double electron_stopping = electron_estar(KE_keV);
    if (!is_positron) return std::max(electron_stopping, 0.001);

    // ESTAR tabulates electrons. Preserve the prior positron treatment by
    // applying only the Bhabha/Møller collision-term ratio to the direct
    // electron value; material, shell, and density corrections cancel.
    const double tau = (KE_keV * 1.0e-3) / kMeC2_MeV;
    const double tau2 = tau * tau;
    const double beta2 = tau * (tau + 2.0) / ((tau + 1.0) * (tau + 1.0));
    const double I_ratio = (mean_excitation_eV(Z) * 1.0e-6) / kMeC2_MeV;
    const double arg = tau2 * (tau + 2.0) / (2.0 * I_ratio * I_ratio);
    if (beta2 <= 0.0 || arg <= 0.0) return electron_stopping;
    const double common = 0.5 * std::log(arg);
    const double f_minus = (1.0 - beta2)
        + (tau2 / 8.0 - (2.0 * tau + 1.0) * std::log(2.0))
          / ((tau + 1.0) * (tau + 1.0));
    const double t2 = tau + 2.0;
    const double f_plus = 2.0 * std::log(2.0)
        - (beta2 / 12.0) * (23.0 + 14.0 / t2 + 10.0 / (t2 * t2)
                            + 4.0 / (t2 * t2 * t2));
    const double denominator = common + f_minus;
    const double ratio = (denominator > 0.0) ? (common + f_plus) / denominator : 1.0;
    return std::max(electron_stopping * ratio, 0.001);
}


// ============================================================
// ElectronCsda construction & range table
// ============================================================

ElectronCsda::ElectronCsda()
{
    // Build log-spaced energy grid from kEMin_keV to kEMax_keV
    double log_min = std::log(kEMin_keV);
    double log_max = std::log(kEMax_keV);
    double log_step = (log_max - log_min) / (kNGrid - 1);

    for (int i = 0; i < kNGrid; ++i) {
        double lE = log_min + i * log_step;
        log_energy_grid_[i] = lE;
        energy_grid_keV_[i] = std::exp(lE);
    }

    // For each element Z = 1..92, build the CSDA range table via trapezoid integration
    // of 1/S(T) over the energy grid.
    //
    // R(T_i) = ∫₀^{T_i} dT' / S(T')
    //
    // Bootstrap: R(T_0) ≈ T_0 / S(T_0)  (assumes S constant from 0 to T_0).
    // This is conservative — the true range below T_0 is smaller since S increases
    // as T→0.  But T_0 = 1 keV electrons have R < 10⁻⁵ g/cm² in any material,
    // so the error is negligible.
    // Build the electron (is_pos=false) and positron (is_pos=true) range tables
    // with identical trapezoid integration; only the F±(τ) stopping term differs.
    for (int pass = 0; pass < 2; ++pass) {
        const bool is_pos = (pass == 1);
        auto& Rtab = is_pos ? range_table_pos_ : range_table_;
        for (int Zm1 = 0; Zm1 < 92; ++Zm1) {
            int Z = Zm1 + 1;
            double A = kAtomicWeight[Zm1];

            double S0 = stopping_power_MeV_cm2_g(Z, A, energy_grid_keV_[0], is_pos);
            // R(T_0) in g/cm² — trapezoid bootstrap from 0 to T_0
            Rtab[Zm1][0] = (energy_grid_keV_[0] * 1e-3) / S0;  // MeV / (MeV cm²/g)

            for (int i = 1; i < kNGrid; ++i) {
                double T0 = energy_grid_keV_[i - 1];
                double T1 = energy_grid_keV_[i];
                double dT = (T1 - T0) * 1e-3;  // keV → MeV
                double S1 = stopping_power_MeV_cm2_g(Z, A, T1, is_pos);
                double S_prev = stopping_power_MeV_cm2_g(Z, A, T0, is_pos);

                // Trapezoid rule: ∫ 1/S dT ≈ (ΔT/2)(1/S(T0) + 1/S(T1))
                double delta_R = 0.5 * dT * (1.0 / S_prev + 1.0 / S1);
                Rtab[Zm1][i] = Rtab[Zm1][i - 1] + delta_R;
            }
        }
    }

    // Build per-element Seltzer-Berger brems normalization integrals.
    // For each element Z and grid energy E, integrate chi(kappa) and chi(kappa)/kappa
    // over kappa = [kBremsThreshold/E, 1] using the SAME 64-point log-spaced
    // trapezoid scheme used previously at runtime, so the compound result matches
    // the old per-substep computation exactly at grid energies.
    constexpr double kBremsThreshold_keV = 10.0;  // matches deposited_in_scoring
    constexpr int    N_int = 64;
    const auto& xs_data = CrossSectionData::instance();
    for (int Zm1 = 0; Zm1 < 92; ++Zm1) {
        int Z = Zm1 + 1;
        for (int i = 0; i < kNGrid; ++i) {
            double E = energy_grid_keV_[i];
            if (E <= kBremsThreshold_keV) {
                sb_Jchi_[Zm1][i]   = 0.0;
                sb_Jchiok_[Zm1][i] = 0.0;
                continue;
            }
            double kappa_min = kBremsThreshold_keV / E;
            double log_kmin  = std::log(kappa_min);
            double dlog      = -log_kmin / N_int;   // positive step; log(1) = 0
            double prev_kap  = kappa_min;
            double prev_chi  = xs_data.sb_chi(Z, E, kappa_min);
            double I_chi = 0.0, I_chiok = 0.0;
            for (int j = 1; j <= N_int; ++j) {
                double kap_j = std::exp(log_kmin + j * dlog);
                double chi_j = xs_data.sb_chi(Z, E, kap_j);
                double dk    = kap_j - prev_kap;
                I_chi   += 0.5 * (prev_chi + chi_j) * dk;
                I_chiok += 0.5 * (prev_chi / prev_kap + chi_j / kap_j) * dk;
                prev_kap = kap_j;
                prev_chi = chi_j;
            }
            sb_Jchi_[Zm1][i]   = I_chi;
            sb_Jchiok_[Zm1][i] = I_chiok;
        }
    }
}

const ElectronCsda& ElectronCsda::instance()
{
    static ElectronCsda singleton;
    return singleton;
}


// ============================================================
// Range interpolation
// ============================================================

double ElectronCsda::interpolate_range(int Z, double KE_keV) const
{
    if (Z < 1 || Z > 92) return 0.0;
    int Zm1 = Z - 1;

    if (KE_keV <= kEMin_keV) return range_table_[Zm1][0];
    if (KE_keV >= kEMax_keV) return range_table_[Zm1][kNGrid - 1];

    // Binary search in log-energy space
    double lE = std::log(KE_keV);
    // Grid is uniformly spaced in log-E
    double log_min = log_energy_grid_[0];
    double log_step = log_energy_grid_[1] - log_energy_grid_[0];
    double frac = (lE - log_min) / log_step;
    int lo = static_cast<int>(frac);
    if (lo < 0) lo = 0;
    if (lo >= kNGrid - 1) lo = kNGrid - 2;
    int hi = lo + 1;

    double t = frac - lo;  // fractional position in [0,1]
    // Linear interpolation in log-range
    double r_lo = range_table_[Zm1][lo];
    double r_hi = range_table_[Zm1][hi];
    // Use linear interpolation in R (range is smooth and nearly linear over each step)
    return r_lo + t * (r_hi - r_lo);
}

double ElectronCsda::interpolate_log_grid(const std::array<double, kNGrid>& tbl,
                                          double KE_keV) const
{
    if (KE_keV <= kEMin_keV) return tbl[0];
    if (KE_keV >= kEMax_keV) return tbl[kNGrid - 1];

    double lE       = std::log(KE_keV);
    double log_min  = log_energy_grid_[0];
    double log_step = log_energy_grid_[1] - log_energy_grid_[0];
    double frac     = (lE - log_min) / log_step;
    int lo = static_cast<int>(frac);
    if (lo < 0) lo = 0;
    if (lo >= kNGrid - 1) lo = kNGrid - 2;
    double t = frac - lo;
    return tbl[lo] + t * (tbl[lo + 1] - tbl[lo]);
}

double ElectronCsda::brems_pemit_corr(const Material& mat, double KE_keV) const
{
    constexpr double k_min = 10.0;  // kBremsThreshold_keV (matches build above)
    if (KE_keV <= k_min) return 1.0;

    // Compound integrals via radiation-weighted sum of per-element tables.
    // w_rad = mass_fraction × Z(Z+1)/A  (same weighting as the old compound_chi).
    double I_chi = 0.0, I_chiok = 0.0, w_sum = 0.0;
    for (const auto& comp : mat.composition()) {
        double A_i = atomic_weight(comp.Z);
        if (A_i <= 0.0) continue;
        double w_rad = comp.mass_fraction * comp.Z * (comp.Z + 1.0) / A_i;
        I_chi   += w_rad * interpolate_log_grid(sb_Jchi_[comp.Z - 1], KE_keV);
        I_chiok += w_rad * interpolate_log_grid(sb_Jchiok_[comp.Z - 1], KE_keV);
        w_sum   += w_rad;
    }
    if (w_sum <= 0.0 || I_chiok <= 0.0) return 1.0;
    I_chi   /= w_sum;
    I_chiok /= w_sum;

    // <k>_1/k = (KE - k_min) / ln(KE/k_min);  <k>_SB = KE × I_chi / I_chi_over_k.
    double mean_k_1k = (KE_keV - k_min) / std::log(KE_keV / k_min);
    double mean_k_SB = KE_keV * I_chi / I_chiok;
    return (mean_k_SB > 0.0) ? (mean_k_1k / mean_k_SB) : 1.0;
}

double ElectronCsda::range_gcm2(int Z, double KE_keV) const
{
    return interpolate_range(Z, KE_keV);
}

const std::array<double, ElectronCsda::kNGrid>&
ElectronCsda::compound_range_table(const Material& mat, bool is_positron) const
{
    // Per-thread cache of compound (Bragg-additive) range tables.  The set of
    // materials in a run is tiny; worker threads are created per compute(), so a
    // thread_local cache has the right lifetime with no locking.  A small
    // composition signature guards against a Material address being reused for a
    // different material across runs.
    struct Entry {
        const Material* key;
        size_t n; int z0; double mf0; double density; bool positron;
        std::array<double, kNGrid> R;
    };
    static thread_local std::vector<Entry> cache = [] {
        std::vector<Entry> v; v.reserve(32); return v;   // reserve => stable refs
    }();

    const auto& comp = mat.composition();
    size_t n   = comp.size();
    int    z0  = n ? static_cast<int>(comp[0].Z) : 0;
    double mf0 = n ? comp[0].mass_fraction : 0.0;
    double den = mat.density();

    for (auto& e : cache) {
        if (e.key == &mat && e.n == n && e.z0 == z0 &&
            e.mf0 == mf0 && e.density == den && e.positron == is_positron) {
            return e.R;
        }
    }

    cache.emplace_back();
    Entry& e = cache.back();
    e.key = &mat; e.n = n; e.z0 = z0; e.mf0 = mf0; e.density = den;
    e.positron = is_positron;
    const auto& src_table = is_positron ? range_table_pos_ : range_table_;

    // Bragg additivity for ranges at each grid energy: 1/R_compound = Σᵢ wᵢ / Rᵢ.
    for (int i = 0; i < kNGrid; ++i) {
        double inv_range_sum = 0.0;
        for (const auto& c : comp) {
            double Ri = src_table[c.Z - 1][i];
            if (Ri > 1e-30) inv_range_sum += c.mass_fraction / Ri;
        }
        e.R[i] = (inv_range_sum > 0.0) ? 1.0 / inv_range_sum : 1e30;
    }
    return e.R;
}

double ElectronCsda::range_gcm2_material(const Material& mat, double KE_keV,
                                         bool is_positron) const
{
    if (KE_keV <= 0.0) return 0.0;
    // Bragg additivity for ranges: 1/R_compound = Σᵢ wᵢ / Rᵢ.
    // Evaluated from the per-thread compound range table (combine-on-grid then
    // interpolate); agrees with the former per-call harmonic sum to grid accuracy.
    return interpolate_log_grid(compound_range_table(mat, is_positron), KE_keV);
}

const std::array<double, ElectronCsda::kNGrid>&
ElectronCsda::compound_radcorr_table(const Material& mat) const
{
    // Per-thread cache of compound radiative-correction tables (deterministic in
    // (material, KE)); same lifetime/signature scheme as compound_range_table.
    struct Entry {
        const Material* key;
        size_t n; int z0; double mf0; double density;
        std::array<double, kNGrid> C;
    };
    static thread_local std::vector<Entry> cache = [] {
        std::vector<Entry> v; v.reserve(32); return v;
    }();

    const auto& comp = mat.composition();
    size_t n   = comp.size();
    int    z0  = n ? static_cast<int>(comp[0].Z) : 0;
    double mf0 = n ? comp[0].mass_fraction : 0.0;
    double den = mat.density();

    for (auto& e : cache) {
        if (e.key == &mat && e.n == n && e.z0 == z0 &&
            e.mf0 == mf0 && e.density == den) {
            return e.C;
        }
    }

    cache.emplace_back();
    Entry& e = cache.back();
    e.key = &mat; e.n = n; e.z0 = z0; e.mf0 = mf0; e.density = den;
    for (int i = 0; i < kNGrid; ++i) {
        e.C[i] = radiative_correction_compound(mat, energy_grid_keV_[i]);
    }
    return e.C;
}

double ElectronCsda::residual_energy_keV(const Material& mat,
                                         double KE_keV,
                                         double path_gcm2,
                                         bool is_positron) const
{
    if (path_gcm2 <= 0.0) return KE_keV;

    // Invert the compound range table directly instead of bisecting.  The table is
    // piecewise-linear in log-energy (same interpolation as interpolate_log_grid),
    // so its inverse is exact per segment: one binary search + one exp(), replacing
    // the former 50-iteration bisection that re-evaluated the compound range (and a
    // std::log per element) each step.
    const std::array<double, kNGrid>& Rtab = compound_range_table(mat, is_positron);
    double R0 = interpolate_log_grid(Rtab, KE_keV);
    double R_remaining = R0 - path_gcm2;
    if (R_remaining <= 0.0) return 0.0;
    if (R_remaining <= Rtab[0]) return kEMin_keV;

    // Largest grid index lo with Rtab[lo] <= R_remaining (Rtab increases with index).
    int lo = 0, hi = kNGrid - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) >> 1;
        if (Rtab[mid] <= R_remaining) lo = mid; else hi = mid;
    }
    double denom = Rtab[lo + 1] - Rtab[lo];
    double t = (denom > 0.0) ? (R_remaining - Rtab[lo]) / denom : 0.0;
    double log_min  = log_energy_grid_[0];
    double log_step = log_energy_grid_[1] - log_energy_grid_[0];
    double T = std::exp(log_min + (lo + t) * log_step);
    // Guard against tiny interpolation overshoot above the incident energy.
    return (T < KE_keV) ? T : KE_keV;
}


// ============================================================
// Radiation length (Tsai formula)
// ============================================================

double ElectronCsda::radiation_length_gcm2_element(int Z, double A_g_mol)
{
    // Tsai formula (PDG Eq. 34.26):
    //   X₀ = 716.4 × A / (Z × (Z+1) × ln(287 / √Z))  [g/cm²]
    if (Z < 1 || A_g_mol <= 0.0) return 1e30;
    double ln_term = std::log(287.0 / std::sqrt(static_cast<double>(Z)));
    if (ln_term <= 0.0) return 1e30;
    return 716.4 * A_g_mol
           / (static_cast<double>(Z) * static_cast<double>(Z + 1) * ln_term);
}

double ElectronCsda::radiation_length_gcm2(const Material& mat)
{
    // Bragg additivity: 1/X₀_compound = Σᵢ wᵢ / X₀ᵢ
    double inv_sum = 0.0;
    for (const auto& comp : mat.composition()) {
        double A  = atomic_weight(comp.Z);
        double X0 = radiation_length_gcm2_element(comp.Z, A);
        if (X0 > 0.0) {
            inv_sum += comp.mass_fraction / X0;
        }
    }
    return (inv_sum > 0.0) ? 1.0 / inv_sum : 1e30;
}


// ============================================================
// Source skin-escape transmission model (June 2026)
// ============================================================
//
// Calibration constants for the empirical electron skin-escape transmission
// gate in walk_in_source_geometry().  See that function for the physics
// rationale.  The detour parameterization (R_ex = R_csda × detour) and the
// transmission curve T(τ) are tuned jointly against the GEANT4 electron-entry
// rate at the two anchor points:
//   cfg 12 @662  (sub-MeV recoils, SS304 skin) — target channel ≈ 2.0e-5/event
//   cfg 11 @3000 (MeV recoils, Fe shell)        — preserve channel (×1.21 vs G4)
// The detour ratio < 1 is the projected penetration depth over the CSDA range;
// it shrinks at lower energy and higher Z (Tabata-Ito-Okabe extrapolated-range
// behaviour, NIM 103 85 1972), which is what makes the same geometric depth map
// to a larger τ at sub-MeV and so collapse T harder there.
namespace {
constexpr double kDetourEhalf_MeV = 0.35;  // energy where the MS deficit half-saturates
constexpr double kDetourQ         = 1.6;   // sharpness of the energy modulation
constexpr double kDetourCZ        = 0.0645;// deficit per unit Z_eff at low energy
constexpr double kDetourZpow      = 1.0;   // Z exponent of the deficit
constexpr double kDetourFloor     = 0.10;  // minimum detour ratio (high-Z/low-E clamp)
constexpr double kEscapeTau50     = 0.90;  // τ at 50% transmission (depth gate)
constexpr double kEscapeP         = 4.0;   // Fermi steepness of the τ collapse
// Surface-emergence (1 − albedo) factor A = 1 − albedo·w_skin(E_exit),
// albedo = min(kSkinAlbedoScale·max(0, η₀(Z) − kSkinEtaFloor), kSkinAlbedoCap),
// evaluated REGIME-AWARE at the EXIT energy + EXIT-material Z (see
// source_escape_survival_exit + skin_escape_window).  The albedo is the
// missing-backscatter correction to the Gaussian-core walk.  Two independent
// conditionings make one model material-general (no Fe-only constant, no
// positron-only flag):
//
//   (1) EXIT-ENERGY window w_skin(E_exit) — the REGIME axis.  The albedo is
//       physical only for skin-DIFFUSION escape (still-energetic exit); it is
//       absent for the forward-peaked MeV channel the Gaussian walk reproduces
//       (w_hi → 0 above ~1.5 MeV, preserving cfg 11's MeV electron entry) and is
//       removed for genuinely range-LIMITED exit (w_lo → 0 as E_exit → 0).
//
//   (2) EXIT-material Z floor — the BACKSCATTER-STRENGTH axis.  Empirically the
//       light-Al β⁺ escapers (Na-22) and the heavy-Fe recoil escapers (cfg 11/12)
//       exit in the SAME energy band (~200–450 keV, measured), so the energy
//       window alone cannot separate them — the true separator is Z.  The Gaussian
//       walk already captures the small/moderate-angle scattering that makes up the
//       light-Z backscatter; it misses only the LARGE-ANGLE tail, whose weight
//       grows above a light-Z baseline.  So the missed (correctable) fraction is
//       max(0, η₀(Z) − kSkinEtaFloor): ≈0 for Al/water (walk accurate → no
//       suppression → Na-22 P511 → 1.00), rising to the Fe-anchored value.
//
// kSkinEtaFloor sits just above η₀(Al)=0.153 (water/cellulose are well below);
// kSkinAlbedoScale is then re-anchored so the Fe deficit reproduces the previously
// validated Fe albedo (min(2.80·η₀(Fe),0.92)=0.783 ⇒ scale·(η₀(Fe)−floor)=0.783),
// keeping cfg 11/12 numerically unchanged.  The first-principles MSC tails (Task 3)
// ultimately remove the need for the empirical albedo entirely.
constexpr double kSkinAlbedoScale = 6.55;  // re-anchored to floored η₀ (Fe albedo preserved)
constexpr double kSkinAlbedoCap   = 0.92;  // albedo cap (A floor = 1 − kSkinAlbedoCap)
constexpr double kSkinEtaFloor    = 0.16;  // light-Z backscatter the walk already captures
constexpr double kSkinEoff_MeV    = 0.45;  // high-E turn-OFF center (preserve MeV channel)
constexpr double kSkinPoff        = 2.0;   // high-E turn-OFF sharpness
constexpr double kSkinEonset_MeV  = 0.080; // low-E turn-ON center (range-limited rolloff)
constexpr double kSkinPonset      = 2.0;   // low-E turn-ON sharpness

// Standard normal-incidence electron backscatter coefficient η₀(Z) (number
// fraction): the cubic Z-fit reproduced in electron-microscopy references
// (e.g. Goldstein et al., "Scanning Electron Microscopy and X-ray
// Microanalysis").  Monotone over Z = 1..92, ≈0 at low Z, saturating to ≈0.5
// at Z = 92 (η: C 0.06, Al 0.15, Fe 0.28, I 0.42, Pb 0.49).  Energy dependence
// is carried separately by ha(E); this is the low-energy plateau shape.
double backscatter_coeff_Z(double Z)
{
    double e = -0.0254 + 0.016 * Z - 1.86e-4 * Z * Z + 8.3e-7 * Z * Z * Z;
    if (e < 0.0)  e = 0.0;
    if (e > 0.55) e = 0.55;
    return e;
}

// Regime window in EXIT kinetic energy (keV): the fraction of the backscatter
// albedo that physically applies, given how the particle reaches the surface.
//   w_lo → 0 as E_exit → 0  : range-limited escape (stopping anyway; the
//                             Gaussian walk is RIGHT that it emerges, so no
//                             albedo suppression — fixes soft positrons / light
//                             low-Z escape that the old birth-E gate killed).
//   w_hi → 0 as E_exit → MeV: forward-peaked skin escape the Gaussian walk
//                             reproduces (preserve cfg 11's MeV channel).
// Product is a band-pass that isolates the sub-MeV skin-DIFFUSION escape the
// walk over-predicts (the only regime the albedo is physical for).
double skin_escape_window(double KE_keV)
{
    double E = KE_keV * 1e-3;  // MeV
    double w_lo = 1.0 / (1.0 + std::pow(kSkinEonset_MeV / std::max(E, 1e-6),
                                        kSkinPonset));
    double w_hi = 1.0 / (1.0 + std::pow(E / kSkinEoff_MeV, kSkinPoff));
    return w_lo * w_hi;
}

// ============================================================
// Source-escape model — COMPILE-TIME selected (CEELO_ESCAPE_MODEL, set by the
// CMake cache var CEELO_SOURCE_ESCAPE_MODEL = gate|tails|gs).
// ============================================================
//   0 GATE  : empirical exit-state (1-albedo) skin-escape gate [default, validated]
//   1 TAILS : straggling + Highland-B2 Gaussian (hard-variance subtracted) +
//             screened-Rutherford hard tail; no empirical gate
//   2 GS    : straggling + Goudsmit-Saunderson soft-moment Gaussian (theta0 from
//             the soft transport coefficient of the screened cross-section,
//             replacing the Highland formula) + screened-Rutherford hard tail;
//             no empirical gate
// TAILS and GS share the straggling + hard tail; they differ only in how the
// Gaussian core width is set (Highland-minus-hard vs soft-transport-moment).
#ifndef CEELO_ESCAPE_MODEL
#define CEELO_ESCAPE_MODEL 0
#endif
constexpr int kSourceEscapeModel = CEELO_ESCAPE_MODEL;
constexpr bool kEscapePrincipled = (kSourceEscapeModel != 0);  // TAILS or GS
constexpr bool kEscapeUseGS       = (kSourceEscapeModel == 2); // GS soft moment

// Bohr energy-loss straggling constant ξ_B = 4π r_e² (mₑc²)² N_A  [MeV² cm²/g].
constexpr double kBohrXi_MeV2cm2g = 0.1569;

// Large-angle cutoff for the explicit single-scatter tail: scatters beyond this
// polar angle are sampled discretely (the Gaussian Highland core handles the
// small-angle bulk; at this angle the Gaussian per-step contribution is small
// for the thin steps here, so the split avoids most double counting).
constexpr double kHardScatterMuCut = 0.940;  // cos(20°)

// Mean-square polar angle of ONE hard screened-Rutherford scatter (μ < μ_c),
// small-η limit: <θ²> = 2·ln(2/(1−μ_c)) / (1/(1−μ_c) − 1/2).  Used to subtract
// the explicit hard tail's variance from the Gaussian Highland core so that
// core + tail reproduces the true MS variance (no double counting — the tail
// must not pile moderate-angle scattering ON TOP of the full Gaussian, which
// over-diffuses and spuriously inflates box-geometry escape).
const double kHardThetaSq =
    2.0 * std::log(2.0 / (1.0 - kHardScatterMuCut))
    / (1.0 / (1.0 - kHardScatterMuCut) - 0.5);

// Moliere/Thomas-Fermi screening parameter η for screened-Rutherford elastic
// scattering: η = (ħc·Z^{1/3} / (2·pc·0.885·a0))².  The bracket constant is
// ħc/(2·0.885·a0) = 197327 keV·fm / (2·0.885·5.2918e4 fm) = 2.1067 keV.
double moliere_screening_eta(double Z, double pc_keV)
{
    if (pc_keV <= 0.0) return 1.0;
    double a = 2.1067 * std::cbrt(Z) / pc_keV;
    return a * a;
}

// Mean number of HARD (μ < μ_c) screened-Rutherford elastic scatters per g/cm²
// for a compound, at electron momentum pc (keV) and total energy E_tot (keV):
//   dσ/dΩ ∝ Z(Z+1) r_e² (mₑc²)² E_tot² / pc⁴ · 1/(1−μ+2η)²
//   σ_hard = π·Z(Z+1) r_e² (mₑc²/pc²)² E_tot² · [1/(1−μ_c+2η) − 1/(2+2η)]
// Σ_i (N_A w_i / A_i)·σ_hard,i.  (The 1/η forward divergence cancels because the
// tail is cut at μ_c, leaving a finite, η-insensitive hard-collision rate.)
double hard_elastic_mean_per_gcm2(const Material& mat, double pc_keV,
                                  double E_tot_keV, double mu_c)
{
    if (pc_keV <= 0.0) return 0.0;
    constexpr double kNA          = 6.02214076e23;     // 1/mol
    constexpr double kReSq_cm2    = 7.9408e-26;        // r_e² (cm²)
    constexpr double kMeC2_keV_   = 510.998950;
    const double kin = kMeC2_keV_ * kMeC2_keV_ * E_tot_keV * E_tot_keV
                     / (pc_keV * pc_keV * pc_keV * pc_keV);
    double sum = 0.0;
    for (const auto& c : mat.composition()) {
        double Z = static_cast<double>(c.Z);
        double A = ElectronCsda::atomic_weight(c.Z);
        if (A <= 0.0) continue;
        double eta   = moliere_screening_eta(Z, pc_keV);
        double tail  = 1.0 / (1.0 - mu_c + 2.0 * eta) - 1.0 / (2.0 + 2.0 * eta);
        if (tail < 0.0) tail = 0.0;
        double sigma = M_PI * Z * (Z + 1.0) * kReSq_cm2 * kin * tail;  // cm²
        sum += (kNA * c.mass_fraction / A) * sigma;                    // 1/(g/cm²)
    }
    return sum;
}

// SOFT transport first-moment coefficient G_1 (per g/cm²) for the screened-
// Rutherford cross-section restricted to μ > μ_c — the "soft" multiple scattering
// the Gaussian core represents (the hard tail μ < μ_c is sampled explicitly):
//   σ_1,soft = π Z(Z+1) r_e² (mₑc²/pc²)² E_tot² · I_soft(η,μ_c),
//   I_soft  = ∫_{μ_c}^1 (1−μ)/(1−μ+2η)² dμ
//           = ln((1−μ_c+2η)/(2η)) + 2η/(1−μ_c+2η) − 1.
// In GS mode the soft Gaussian-core projected variance is θ0² = G_1,soft·Δx —
// the Goudsmit-Saunderson first moment for the step, computed from the SAME
// cross-section as the hard tail (no Highland formula, no double counting).
double soft_transport_per_gcm2(const Material& mat, double pc_keV,
                               double E_tot_keV, double mu_c)
{
    if (pc_keV <= 0.0) return 0.0;
    constexpr double kNA        = 6.02214076e23;
    constexpr double kReSq_cm2  = 7.9408e-26;
    constexpr double kMeC2_keV_ = 510.998950;
    const double kin = kMeC2_keV_ * kMeC2_keV_ * E_tot_keV * E_tot_keV
                     / (pc_keV * pc_keV * pc_keV * pc_keV);
    double sum = 0.0;
    for (const auto& c : mat.composition()) {
        double Z = static_cast<double>(c.Z);
        double A = ElectronCsda::atomic_weight(c.Z);
        if (A <= 0.0) continue;
        double eta    = moliere_screening_eta(Z, pc_keV);
        double denom  = 1.0 - mu_c + 2.0 * eta;
        double I_soft = std::log(denom / (2.0 * eta)) + 2.0 * eta / denom - 1.0;
        if (I_soft < 0.0) I_soft = 0.0;
        double sigma1 = M_PI * Z * (Z + 1.0) * kReSq_cm2 * kin * I_soft;  // cm²
        sum += (kNA * c.mass_fraction / A) * sigma1;                      // 1/(g/cm²)
    }
    return sum;
}

// Sample μ = cosθ of one hard screened-Rutherford scatter restricted to μ < μ_c,
// from the unrestricted inverse-CDF μ = 1 − 2η ξ/(1+η−ξ) confined to ξ ∈ [ξ_c,1].
double sample_hard_mu(double eta, double mu_c, double xi01)
{
    double denom = 2.0 * eta + 1.0 - mu_c;
    double xi_c  = (denom > 0.0) ? (1.0 - mu_c) * (1.0 + eta) / denom : 0.0;
    double x     = xi_c + (1.0 - xi_c) * xi01;          // ∈ [ξ_c, 1]
    double mu    = 1.0 - 2.0 * eta * x / (1.0 + eta - x);
    if (mu < -1.0) mu = -1.0;
    if (mu >  1.0) mu =  1.0;
    return mu;
}
}  // namespace

double ElectronCsda::effective_Z(const Material& mat)
{
    double num = 0.0, den = 0.0;
    for (const auto& c : mat.composition()) {
        num += c.mass_fraction * static_cast<double>(c.Z);
        den += c.mass_fraction;
    }
    return (den > 0.0) ? num / den : 7.0;
}

double ElectronCsda::extrapolated_range_gcm2_material(const Material& mat,
                                                      double KE_keV) const
{
    double Rcsda = range_gcm2_material(mat, KE_keV);
    if (Rcsda <= 0.0) return 0.0;
    double Z     = effective_Z(mat);
    double E_MeV = KE_keV * 1e-3;
    // Energy modulation of the multiple-scattering detour deficit: ~1 at low
    // energy, → 0 at high energy (projected range approaches the CSDA range).
    double gE = 1.0 / (1.0 + std::pow(E_MeV / kDetourEhalf_MeV, kDetourQ));
    double deficit = kDetourCZ * std::pow(Z, kDetourZpow) * gE;
    double detour  = 1.0 - deficit;
    if (detour < kDetourFloor) detour = kDetourFloor;
    if (detour > 0.98)         detour = 0.98;
    return Rcsda * detour;
}

double ElectronCsda::source_escape_survival_exit(double exit_KE_keV,
                                                 double Z_exit,
                                                 double tau)
{
    // The empirical gate is a HEAVY-material over-prediction fix: it exists only
    // because the Gaussian-core walk misses the large-angle backscatter that turns
    // escapers back in dense, high-Z shields.  Both pieces below (the surface
    // albedo AND the depth gate) are that same fix.  For light materials the walk
    // is accurate, so NEITHER applies — fully analog escape.  Z is the EXIT
    // material's effective Z (the skin the particle reflects off).  Only the
    // large-angle backscatter the walk MISSES is correctable: the deficit of η₀(Z)
    // above the light-Z baseline the walk already captures.  η₀(Z) ≤ floor (Al,
    // water, cellulose) ⇒ survival 1 — this is what lets soft light-Z escape
    // (Na-22 Al positrons, water sources) match GEANT4 without the legacy
    // per-caller skip flag, and it removes the depth gate there too.
    double eta = backscatter_coeff_Z(Z_exit) - kSkinEtaFloor;
    if (eta <= 0.0) return 1.0;

    // Heavy material: depth gate T_path(τ) (Fermi/Tabata-style, 1 at τ=0, 0.5 at
    // τ=τ50, collapsing for τ ≫ τ50 — electrons born deep in the shield) times the
    // surface-emergence factor A = 1 − albedo.  The albedo is REGIME-AWARE via
    // w_skin(E_exit): zeroed for range-LIMITED exit (E_exit → 0, stopping anyway)
    // and for the forward-peaked MeV channel, isolating the sub-MeV over-escape.
    double T_path = (tau <= 0.0)
                      ? 1.0
                      : 1.0 / (1.0 + std::pow(tau / kEscapeTau50, kEscapeP));
    double albedo = std::min(kSkinAlbedoScale * eta, kSkinAlbedoCap)
                  * skin_escape_window(exit_KE_keV);
    double A      = 1.0 - albedo;
    if (A < 0.0) A = 0.0;

    return A * T_path;
}


// ============================================================
// CSDA energy deposition in scoring volume
// ============================================================

ElectronDepositResult ElectronCsda::deposited_in_scoring(
                                           const Eigen::Vector3d& position,
                                           const Eigen::Vector3d& direction,
                                           double KE_keV,
                                           const Geometry& geometry,
                                           std::mt19937_64& rng,
                                           bool disable_moliere,
                                           bool disable_brems,
                                           bool is_positron) const
{
    ElectronDepositResult result;
    result.stop_position = position;   // default: no walk → annihilate at vertex
    if (KE_keV < 0.01) return result;  // < 10 eV: negligible
    std::vector<PathSegment> segments;

    // ----------------------------------------------------------------
    // Low energy (KE ≤ 200 keV): straight-line CSDA.
    // At these energies the electron range is < 0.1 mm in typical
    // detector materials, so angular deflection is negligible.
    // No bremsstrahlung emitted (radiative losses negligible below 200 keV).
    // ----------------------------------------------------------------
    // Also use straight-line path for all energies if Molière is disabled.
    if (KE_keV <= kMoliereBremsThreshold_keV || disable_moliere) {
        geometry.trace_ray(position, direction, segments);
        if (segments.empty()) return result;

        double remaining_KE = KE_keV;

        for (const auto& seg : segments) {
            if (remaining_KE < 0.01) break;
            if (!seg.material) continue;
            if (seg.t_end < -1e-8) continue;

            double t_start_eff = std::max(0.0, seg.t_start);
            double seg_len     = seg.t_end - t_start_eff;
            if (seg_len <= 0.0) continue;

            double density   = seg.material->density();
            double path_gcm2 = density * seg_len;
            double range_g   = range_gcm2_material(*seg.material, remaining_KE, is_positron);

            if (path_gcm2 >= range_g) {
                if (seg.is_scoring) result.deposited_scoring_keV += remaining_KE;
                // Stops inside this segment: advance from the segment entry by the
                // residual CSDA range (g/cm² → cm via the segment density).
                result.stop_position =
                    position + direction * (t_start_eff + range_g / density);
                remaining_KE = 0.0;
                break;
            } else {
                double exit_KE  = residual_energy_keV(*seg.material, remaining_KE, path_gcm2, is_positron);
                double delta_KE = remaining_KE - exit_KE;
                if (delta_KE < 0.0) delta_KE = 0.0;
                if (seg.is_scoring) result.deposited_scoring_keV += delta_KE;
                // Crosses this segment; provisional stop point is its far face,
                // updated again if a later segment stops the electron.
                result.stop_position = position + direction * seg.t_end;
                remaining_KE = exit_KE;
            }
        }
        return result;
    }

    // ----------------------------------------------------------------
    // High energy (KE > 200 keV): Molière condensed-history walk.
    // Adaptive sub-steps (more at higher energy for finer brems sampling);
    // each step applies Highland multiple scattering and may emit a discrete
    // bremsstrahlung photon with NIST-corrected S_rad.  Energy loss per step is
    // the deterministic CSDA value — no stochastic energy-loss straggling is
    // applied yet (see the kBohrXi TODO above).
    // ----------------------------------------------------------------
    constexpr double m_e     = 511.0;  // electron rest mass (keV)
    constexpr double kBremsThreshold_keV = 10.0;  // minimum brems photon energy

    // Adaptive step count: 20 below 1 MeV, scaling up to ~60 at 10 MeV.
    // More steps at high energy gives finer brems sampling and avoids
    // overcorrection from single large brems photons.
    int n_steps = std::max(20, std::min(60, static_cast<int>(KE_keV / 50.0)));

    Eigen::Vector3d pos = position;
    Eigen::Vector3d dir = direction;
    double KE = KE_keV;
    int steps_remaining = n_steps;
    std::normal_distribution<double> gauss(0.0, 1.0);
    thread_local std::uniform_real_distribution<double> uniform(0.0, 1.0);

    while (steps_remaining > 0 && KE > 0.01) {
        geometry.trace_ray(pos, dir, segments);
        if (segments.empty()) break;

        // Find the segment we are currently in (t_start ≤ 0 < t_end)
        const PathSegment* cur = nullptr;
        for (const auto& seg : segments) {
            if (seg.t_start < 1e-8 || (seg.t_end > 1e-8 && seg.t_start <= 1e-8)) {
                cur = &seg;
                break;
            }
        }
        if (!cur) {
            // Not inside any material — advance to the next segment start
            for (const auto& seg : segments) {
                if (seg.t_start > 1e-8) {
                    pos += dir * (seg.t_start + 1e-9);
                    break;
                }
            }
            continue;
        }

        if (!cur->material) break;

        const Material& mat  = *cur->material;
        double density       = mat.density();
        double t_start_eff   = std::max(0.0, cur->t_start);
        double seg_len       = cur->t_end - t_start_eff;
        if (seg_len <= 0.0) { steps_remaining--; continue; }

        double seg_gcm2 = density * seg_len;
        double range_g  = range_gcm2_material(mat, KE, is_positron);

        // Step size: fraction of remaining range, or full segment if we reach boundary
        double planned_step_gcm2 = range_g / static_cast<double>(steps_remaining);
        double step_gcm2;
        bool reach_boundary = false;

        if (planned_step_gcm2 >= seg_gcm2) {
            step_gcm2 = seg_gcm2;
            reach_boundary = true;
        } else {
            step_gcm2 = planned_step_gcm2;
        }
        if (step_gcm2 <= 0.0) { steps_remaining--; continue; }

        // In-flight annihilation (positron only): probability of annihilating
        // over this step is n_e·σ_a(KE)·Δx, with Δx = step_gcm2/density.  Sampled
        // at the step-entry position/energy so the positron deposits nothing
        // further and the photon pair (handled by the caller) carries 2mₑc²+KE.
        // Covers the whole slowing-down because the walk runs to rest; positrons
        // born ≤200 keV (straight-line branch) skip this — their in-flight
        // fraction is tiny and the annihilation site is ~0.1 mm from the vertex.
        if (is_positron) {
            double dx_cm   = step_gcm2 / density;
            double P_annih = electron_number_density(mat)
                           * positron_annih_xsec_per_electron(KE) * dx_cm;
            if (P_annih > 0.0 && uniform(rng) < std::min(P_annih, 1.0)) {
                result.annihilated_in_flight   = true;
                result.residual_KE_at_annih_keV = KE;
                result.stop_position           = pos;  // step-entry site
                return result;
            }
        }

        // Advance position
        pos += dir * (step_gcm2 / density);

        // Push across boundary so next trace_ray finds the next volume
        if (reach_boundary) {
            pos += dir * 1e-6;
        }

        // Radiation length for this material (needed for brems and scattering)
        double X0 = radiation_length_gcm2(mat);

        // --- Bremsstrahlung emission (probabilistic, Seltzer-Berger) ---
        //
        // P_emit is derived from the 1/k envelope: N = (Δx/X0) × ln(KE/k_min) × rad_corr.
        // With Seltzer-Berger rejection sampling, the accepted photons have lower mean
        // energy <k>_SB < <k>_1/k.  To maintain the correct total radiative energy loss
        // (S_rad × Δx), we scale P_emit by <k>_1/k / <k>_SB.  This increases the
        // emission rate so that P_corrected × <k>_SB = P_old × <k>_1/k = S_rad × Δx.
        //
        double brems_keV = 0.0;
        if (X0 > 0.0 && !disable_brems && KE > kBremsThreshold_keV) {
            double k_min = kBremsThreshold_keV;
            double k_max = KE;
            double kappa_min = k_min / KE;
            const auto& xs_data = CrossSectionData::instance();

            // Helper: radiation-weight-averaged compound chi at given kappa
            // Hoisted compound SB chi (Z-independent energy bracket + per-element
            // radiation weights precomputed once; bit-identical to the old lambda).
            CompoundBremsChi compound_chi(xs_data, mat, KE);

            // P_emit correction <k>_1/k / <k>_SB, from precomputed per-element
            // Seltzer-Berger integral tables (replaces the former per-substep
            // 64-point numerical integration that dominated runtime via sb_chi).
            double P_emit_corr = brems_pemit_corr(mat, KE);

            double rad_corr = interpolate_log_grid(compound_radcorr_table(mat), KE);
            double P_emit = (step_gcm2 / X0)
                          * std::log(KE / kBremsThreshold_keV)
                          * rad_corr
                          * P_emit_corr;
            P_emit = std::min(P_emit, 1.0);

            if (uniform(rng) < P_emit) {
                // Rejection-sampling envelope ceiling: chi at the soft-photon limit
                // (kappa_min).  Computed only on the rare steps that actually emit a
                // photon, so it no longer contributes to the per-substep hot path.
                double chi_max = compound_chi(kappa_min);
                if (chi_max <= 0.0) chi_max = 1.0;

                // Rejection sampling: propose from 1/k, accept with chi(κ)/chi_max
                double E_brems = 0.0;
                bool accepted = false;
                for (int attempt = 0; attempt < 100; ++attempt) {
                    double xi = uniform(rng);
                    double k_proposed = k_min * std::pow(k_max / k_min, xi);
                    double kappa = k_proposed / KE;
                    double chi_val = compound_chi(kappa);

                    if (uniform(rng) * chi_max < chi_val) {
                        E_brems = std::min(k_proposed, KE);
                        accepted = true;
                        break;
                    }
                }
                if (accepted) {
                    result.brems_photons.push_back({pos, dir, E_brems});
                    brems_keV = E_brems;
                }
            }
        }

        // CSDA collision loss (computed from current KE, which includes collision
        // stopping only — brems is a separate radiative process).
        double KE_exit  = residual_energy_keV(mat, KE, step_gcm2, is_positron);
        double delta_KE = KE - KE_exit;
        if (delta_KE < 0.0) delta_KE = 0.0;

        // Collision energy is deposited locally (as heat).
        // Brems photon carries away energy that may escape the crystal.
        // To conserve energy: if brems was emitted, the collision deposit
        // must be reduced so that deposit + brems ≤ KE for this step.
        double deposit = delta_KE;
        if (brems_keV > 0.0) {
            // Total energy removed from electron this step = delta_KE + brems_keV.
            // Cap brems so it doesn't exceed remaining KE after collision:
            double max_brems = std::max(0.0, KE - delta_KE);
            if (brems_keV > max_brems) {
                // Reduce deposit to accommodate brems
                deposit = std::max(0.0, KE - brems_keV);
            }
        }
        if (cur->is_scoring) {
            result.deposited_scoring_keV += deposit;
        }
        KE = std::max(0.0, KE - deposit - brems_keV);
        if (KE <= 0.01) break;  // electron stopped

        // Highland-Molière scattering kick
        // θ₀ = (13.6 MeV / βcp) × √(x/X₀) × (1 + 0.038 ln(x/X₀))
        // βcp = √(KE × (KE + 2mₑ))   [all in keV]
        if (X0 > 0.0) {
            double x_X0  = step_gcm2 / X0;
            double bcp   = std::sqrt(KE * (KE + 2.0 * m_e));  // β×p×c in keV
            if (bcp > 1e-6 && x_X0 > 0.0) {
                double ln_x   = std::log(x_X0);
                double theta0 = (13600.0 / bcp) * std::sqrt(x_X0)
                                * (1.0 + 0.038 * ln_x);  // radians
                if (theta0 > 1e-10) {
                    // Sample two independent projected angles ~ N(0, θ₀)
                    double dx    = gauss(rng) * theta0;
                    double dy    = gauss(rng) * theta0;
                    double theta = std::sqrt(dx * dx + dy * dy);  // total deflection
                    if (theta > 1e-10) {
                        Eigen::Vector3d u, v;
                        build_perp_basis(dir, u, v);
                        Eigen::Vector3d dir_new =
                            std::cos(theta) * dir
                            + std::sin(theta) * (dx * u + dy * v) / theta;
                        dir = dir_new.normalized();
                    }
                }
            }
        }

        steps_remaining--;
    }

    // Terminal position of the condensed-history walk: the rest point if the
    // electron stopped (KE → 0), otherwise the last boundary it reached before
    // leaving all geometry.  Used as the positron annihilation site by callers.
    result.stop_position = pos;
    return result;
}


std::vector<BremsPhoton> ElectronCsda::sample_brems_in_material(
                                           const Material& mat,
                                           const Eigen::Vector3d& position,
                                           const Eigen::Vector3d& direction,
                                           double KE_keV,
                                           std::mt19937_64& rng) const
{
    std::vector<BremsPhoton> photons;
    if (KE_keV <= kMoliereBremsThreshold_keV) return photons;

    double X0 = radiation_length_gcm2(mat);
    if (X0 <= 0.0) return photons;

    constexpr double m_e = 511.0;  // electron rest mass (keV)
    constexpr double kBremsThreshold_keV = 10.0;  // minimum brems photon energy

    // Same condensed-history stepping as deposited_in_scoring(), minus the
    // geometry: the electron stops entirely inside `mat`, so each step is a
    // fixed fraction of the remaining CSDA range.
    int n_steps = std::max(20, std::min(60, static_cast<int>(KE_keV / 50.0)));

    Eigen::Vector3d dir = direction;
    double KE = KE_keV;
    int steps_remaining = n_steps;
    std::normal_distribution<double> gauss(0.0, 1.0);
    thread_local std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const auto& xs_data = CrossSectionData::instance();

    while (steps_remaining > 0 && KE > kBremsThreshold_keV) {
        double range_g   = range_gcm2_material(mat, KE);
        double step_gcm2 = range_g / static_cast<double>(steps_remaining);
        if (step_gcm2 <= 0.0) break;

        // --- Bremsstrahlung emission (identical to deposited_in_scoring) ---
        double kappa_min = kBremsThreshold_keV / KE;

        // Hoisted compound SB chi (Z-independent energy bracket + per-element
        // radiation weights precomputed once; bit-identical to the old lambda).
        CompoundBremsChi compound_chi(xs_data, mat, KE);

        double P_emit_corr = brems_pemit_corr(mat, KE);
        double rad_corr = interpolate_log_grid(compound_radcorr_table(mat), KE);
        double P_emit = (step_gcm2 / X0)
                      * std::log(KE / kBremsThreshold_keV)
                      * rad_corr
                      * P_emit_corr;
        P_emit = std::min(P_emit, 1.0);

        double brems_keV = 0.0;
        if (uniform(rng) < P_emit) {
            double chi_max = compound_chi(kappa_min);
            if (chi_max <= 0.0) chi_max = 1.0;

            for (int attempt = 0; attempt < 100; ++attempt) {
                double xi = uniform(rng);
                double k_proposed = kBremsThreshold_keV
                    * std::pow(KE / kBremsThreshold_keV, xi);
                double kappa = k_proposed / KE;
                double chi_val = compound_chi(kappa);

                if (uniform(rng) * chi_max < chi_val) {
                    double E_brems = std::min(k_proposed, KE);
                    // Creation-point approximation: report at the electron's
                    // birth position with the current (Highland-evolved)
                    // direction.
                    photons.push_back({position, dir, E_brems});
                    brems_keV = E_brems;
                    break;
                }
            }
        }

        // CSDA collision loss; conserve energy when a brems photon was emitted
        // (same bookkeeping as deposited_in_scoring).
        double KE_exit  = residual_energy_keV(mat, KE, step_gcm2);
        double delta_KE = KE - KE_exit;
        if (delta_KE < 0.0) delta_KE = 0.0;
        double deposit = delta_KE;
        if (brems_keV > 0.0 && brems_keV > KE - delta_KE) {
            deposit = std::max(0.0, KE - brems_keV);
        }
        KE = std::max(0.0, KE - deposit - brems_keV);
        if (KE <= 0.01) break;

        // Highland-Molière scattering kick on the electron direction so that
        // late-emitted photons get a realistic (diffused) direction.
        double x_X0 = step_gcm2 / X0;
        double bcp  = std::sqrt(KE * (KE + 2.0 * m_e));
        if (bcp > 1e-6 && x_X0 > 0.0) {
            double ln_x   = std::log(x_X0);
            double theta0 = (13600.0 / bcp) * std::sqrt(x_X0)
                            * (1.0 + 0.038 * ln_x);
            if (theta0 > 1e-10) {
                double dx    = gauss(rng) * theta0;
                double dy    = gauss(rng) * theta0;
                double theta = std::sqrt(dx * dx + dy * dy);
                if (theta > 1e-10) {
                    Eigen::Vector3d u, v;
                    build_perp_basis(dir, u, v);
                    Eigen::Vector3d dir_new =
                        std::cos(theta) * dir
                        + std::sin(theta) * (dx * u + dy * v) / theta;
                    dir = dir_new.normalized();
                }
            }
        }

        steps_remaining--;
    }

    return photons;
}


ElectronSourceWalkResult ElectronCsda::walk_in_source_geometry(
                                           const SourceGeometry& source,
                                           const Material& birth_material,
                                           const Eigen::Vector3d& position,
                                           const Eigen::Vector3d& direction,
                                           double KE_keV,
                                           std::mt19937_64& rng) const
{
    ElectronSourceWalkResult result;
    if (KE_keV < 0.01) return result;

    constexpr double m_e = 511.0;  // electron rest mass (keV)
    constexpr double kBremsThreshold_keV = 10.0;  // minimum brems photon energy

    // Source-escape model (COMPILE-TIME, kSourceEscapeModel): TAILS/GS add
    // straggling + Rutherford tails and bypass the empirical gate; GS additionally
    // sets the Gaussian core from the soft transport moment instead of Highland.
    constexpr bool principled = kEscapePrincipled;
    constexpr bool use_gs     = kEscapeUseGS;
    // Component A/B knobs (env, read once) — isolate straggling / Highland-B2 /
    // Rutherford tails for diagnosis and tuning of the TAILS/GS models.
    static const bool dbg_no_straggle = std::getenv("MCDET_NO_STRAGGLE") != nullptr;
    static const bool dbg_no_b2       = std::getenv("MCDET_NO_B2") != nullptr;
    static const bool dbg_no_tails    = std::getenv("MCDET_NO_TAILS") != nullptr;

    Eigen::Vector3d pos = position;
    Eigen::Vector3d dir = direction;
    double KE = KE_keV;

    // Reusable trace buffer (E2): trace_source_segments fills this in place
    // instead of heap-allocating a fresh vector per call. Shared across the
    // escape-tau trace, the low-energy straight-line trace, and the per-substep
    // Molière-walk traces below (each call clear()s it first; the escape-tau
    // result is fully consumed into tau before the walk loop begins, and the
    // low-energy branch returns before the loop, so there is no aliasing).
    std::vector<SourceGeometry::SourcePathSegment> seg_buf;

    // ----------------------------------------------------------------
    // Regime-aware skin-escape survival gate (analog absorption).
    //
    // The Gaussian-core Highland walk below carries no large-angle single-
    // scatter tails, no backscatter, and no energy-loss straggling, so it
    // overestimates escape of sub-MeV recoil electrons from the outer skin of a
    // dense shield ~5× vs GEANT4 (cfg 12 @662, June 2026 residual triage),
    // while the MeV-scale channel (cfg 11 @3000) is correct to ~1.2× and must
    // be preserved.  We restore that differential with a single analog survival
    // test, killing the escaper (it stops in the unscored source) with
    // probability 1 − T, where
    //   T = A(E_exit, Z_exit) · T_path(τ).
    //   • A(E_exit, Z_exit) is the surface-emergence (1 − backscatter) factor,
    //     evaluated at the EXIT state (not the birth energy).  Keying on exit
    //     energy is what makes the gate REGIME-AWARE: the albedo applies only to
    //     skin-DIFFUSION escape (still-energetic exit), and vanishes for range-
    //     LIMITED escape (E_exit→0, the particle is stopping anyway — no
    //     backscatter physics) and for the MeV channel.  This subsumes both the
    //     old Fe-only constant and the positron-only skip flag in one model.
    //   • τ = Σ_layer (g/cm² path to outer boundary)/R_ex (birth depth); T_path(τ)
    //     collapses near τ≈1 (Tabata-style) for electrons born deeper in (cfg 11's
    //     shell escapers, τ ≈ 0.5).  R_ex/R_csda shrinks at lower energy, so the
    //     same depth maps to a larger τ at sub-MeV.  Retained as-is.
    // Ordinary attenuation → unbiased and FEP-safe (no per-branch weights).
    // The test gates only the ESCAPE outcome, not the walk: it is applied at the
    // escape point below (using the exit state there) so that bremsstrahlung
    // emitted along the walk is kept even for escapers that are suppressed (the
    // walk is also run for brems when source-electron escape is disabled).
    // ----------------------------------------------------------------
    double tau = 0.0;
    {
        source.trace_source_segments(pos, dir, KE, seg_buf);
        for (const auto& seg : seg_buf) {
            if (!seg.material || seg.length <= 0.0) continue;
            double path_gcm2 = seg.length * seg.material->density();
            double Rex = extrapolated_range_gcm2_material(*seg.material, KE);
            if (Rex > 1e-12) tau += path_gcm2 / Rex;
        }
    }
    // Survival decision (true = the escape is physical / kept), evaluated from
    // the EXIT state (residual KE + exit-material Z) at the escape point.
    // Principled mode: the walk itself now produces backscatter (Rutherford
    // tails), so the empirical gate is bypassed — every walk-reported escape is
    // accepted.
    auto escape_survives = [&](double exit_KE, double Z_exit) -> bool {
        if (principled) return true;
        double s = source_escape_survival_exit(exit_KE, Z_exit, tau);
        if (s >= 1.0) return true;
        if (s <= 0.0) return false;
        thread_local std::uniform_real_distribution<double> u01(0.0, 1.0);
        return u01(rng) <= s;
    };

    // ----------------------------------------------------------------
    // Very low energy (KE ≤ 50 keV): straight-line CSDA, no brems.
    //
    // Unlike deposited_in_scoring (which switches to straight lines at
    // 200 keV — the path shape is irrelevant for deposits inside the
    // scoring crystal), the source-escape walk keeps Moliere direction
    // diffusion down to 50 keV: escape through the outer skin of a
    // dense shield is dominated by detour and diffusion, and the
    // straight-line model overestimated the sub-MeV escape channel
    // ~4.5x vs GEANT4 (config 12 @ 662 keV, June 2026 — it produced a
    // +1.8% total-efficiency excess). Below 50 keV escape no longer
    // counts (set_source_electron_threshold default) and brems is
    // negligible, so straight-line is safe there.
    // ----------------------------------------------------------------
    constexpr double kSourceWalkStraightLineFloor_keV = 50.0;
    if (KE <= kSourceWalkStraightLineFloor_keV) {
        source.trace_source_segments(pos, dir, KE, seg_buf);
        double dist = 0.0;
        const Material* exit_mat = &birth_material;  // exit-Z source for the gate
        for (const auto& seg : seg_buf) {
            if (!seg.material || seg.length <= 0.0) continue;
            double path_gcm2 = seg.length * seg.material->density();
            double range_g = range_gcm2_material(*seg.material, KE);
            if (path_gcm2 >= range_g) return result;  // stopped inside
            KE = residual_energy_keV(*seg.material, KE, path_gcm2);
            exit_mat = seg.material;
            dist += seg.length;
            if (KE < 1.0) return result;
        }
        // skin-escape gate (no brems below 50 keV); KE is the residual exit energy
        if (escape_survives(KE, effective_Z(*exit_mat))) {
            result.escaped = true;
            result.exit_position = pos + dir * (dist + 1e-7);
            result.exit_direction = dir;
            result.exit_KE_keV = KE;
        }
        return result;
    }

    // ----------------------------------------------------------------
    // High energy: Molière condensed-history walk through the source
    // segments, with Seltzer-Berger brems (same physics constants and
    // bookkeeping as deposited_in_scoring).
    // ----------------------------------------------------------------
    int n_steps = std::max(20, std::min(60, static_cast<int>(KE / 50.0)));
    int steps_remaining = n_steps;
    std::normal_distribution<double> gauss(0.0, 1.0);
    thread_local std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const auto& xs_data = CrossSectionData::instance();

    // Material the particle is currently traversing — the EXIT-Z source for the
    // skin-escape gate when the trace empties (the boundary it reflects off).
    const Material* last_material = &birth_material;

    // Phase D: cumulative path in radiation lengths, for the full-path Highland
    // log term (B2) — the per-step ln(x/X0) under-scatters the cumulative angle.
    double sum_x_X0 = 0.0;

    while (steps_remaining > 0 && KE > 0.01) {
        // Only the first (current) segment is needed here -- max_segments=1 stops
        // the trace after it, avoiding a full multi-shell re-trace every substep.
        source.trace_source_segments(pos, dir, KE, seg_buf, /*max_segments=*/1);
        if (seg_buf.empty() || !seg_buf[0].material || seg_buf[0].length <= 1e-10) {
            // Outside the source geometry — escaped (subject to the skin-escape
            // gate; brems already emitted along the walk is retained in result).
            // KE is the residual exit energy; gate on the exit-material Z.
            if (escape_survives(KE, effective_Z(*last_material))) {
                result.escaped = true;
                result.exit_position = pos;
                result.exit_direction = dir;
                result.exit_KE_keV = KE;
            }
            return result;
        }

        const Material& mat = *seg_buf[0].material;
        last_material = &mat;
        double density   = mat.density();
        double seg_gcm2  = seg_buf[0].length * density;
        double range_g   = range_gcm2_material(mat, KE);

        double planned_step_gcm2 = range_g / static_cast<double>(steps_remaining);
        double step_gcm2 = planned_step_gcm2;
        bool reach_boundary = false;
        if (planned_step_gcm2 >= seg_gcm2) {
            step_gcm2 = seg_gcm2;
            reach_boundary = true;
        }
        if (step_gcm2 <= 0.0) { steps_remaining--; continue; }

        pos += dir * (step_gcm2 / density);
        if (reach_boundary) pos += dir * 1e-6;

        double X0 = radiation_length_gcm2(mat);

        // --- Bremsstrahlung emission (identical to deposited_in_scoring) ---
        double brems_keV = 0.0;
        if (X0 > 0.0 && KE > kBremsThreshold_keV) {
            double kappa_min = kBremsThreshold_keV / KE;

            // Hoisted compound SB chi (Z-independent energy bracket + per-element
            // radiation weights precomputed once; bit-identical to the old lambda).
            CompoundBremsChi compound_chi(xs_data, mat, KE);

            double P_emit_corr = brems_pemit_corr(mat, KE);
            double rad_corr = interpolate_log_grid(compound_radcorr_table(mat), KE);
            double P_emit = (step_gcm2 / X0)
                          * std::log(KE / kBremsThreshold_keV)
                          * rad_corr
                          * P_emit_corr;
            P_emit = std::min(P_emit, 1.0);

            if (uniform(rng) < P_emit) {
                double chi_max = compound_chi(kappa_min);
                if (chi_max <= 0.0) chi_max = 1.0;

                for (int attempt = 0; attempt < 100; ++attempt) {
                    double xi = uniform(rng);
                    double k_proposed = kBremsThreshold_keV
                        * std::pow(KE / kBremsThreshold_keV, xi);
                    double kappa = k_proposed / KE;
                    double chi_val = compound_chi(kappa);

                    if (uniform(rng) * chi_max < chi_val) {
                        double E_brems = std::min(k_proposed, KE);
                        result.brems_photons.push_back({pos, dir, E_brems});
                        brems_keV = E_brems;
                        break;
                    }
                }
            }
        }

        // CSDA collision loss; conserve energy with the emitted brems photon.
        double KE_exit  = residual_energy_keV(mat, KE, step_gcm2);
        double delta_KE = KE - KE_exit;
        if (delta_KE < 0.0) delta_KE = 0.0;
        // Phase D1: per-step Gaussian energy-loss straggling (Bohr variance
        // σ² = ξ_B·(Z/A)·Δx·(1−β²/2)); lets near-threshold electrons stop/escape
        // against the deterministic CSDA value. Principled mode only.
        if (principled && !dbg_no_straggle && step_gcm2 > 1e-12) {
            double beta2  = KE * (KE + 2.0 * m_e) / ((KE + m_e) * (KE + m_e));
            double ZoverA = 0.0;
            for (const auto& c : mat.composition())
                ZoverA += c.mass_fraction * c.Z / atomic_weight(c.Z);
            double sigma2 = kBohrXi_MeV2cm2g * ZoverA * step_gcm2
                          * (1.0 - 0.5 * beta2);                 // MeV²
            if (sigma2 > 0.0) {
                delta_KE += gauss(rng) * std::sqrt(sigma2) * 1000.0;  // keV
                if (delta_KE < 0.0)  delta_KE = 0.0;
                if (delta_KE > KE)   delta_KE = KE;
            }
        }
        double deposit = delta_KE;
        if (brems_keV > 0.0 && brems_keV > KE - delta_KE) {
            deposit = std::max(0.0, KE - brems_keV);
        }
        KE = std::max(0.0, KE - deposit - brems_keV);
        if (KE <= 0.01) break;  // electron stopped inside

        // Highland-Molière scattering kick (Gaussian core) + explicit hard tail.
        if (X0 > 0.0) {
            double x_X0 = step_gcm2 / X0;
            sum_x_X0 += x_X0;
            double bcp  = std::sqrt(KE * (KE + 2.0 * m_e));

            // Phase D3 hard-scatter mean for this step (also feeds the Gaussian
            // variance subtraction below, so compute it first).
            double hard_mean = 0.0;
            if (principled && !dbg_no_tails) {
                hard_mean = hard_elastic_mean_per_gcm2(
                                mat, bcp, KE + m_e, kHardScatterMuCut) * step_gcm2;
                if (hard_mean > 50.0) hard_mean = 50.0;  // guard pathological
            }

            if (bcp > 1e-6 && x_X0 > 0.0) {
                double theta0;
                if (use_gs) {
                    // GS: soft-Gaussian projected variance θ0² = G_1,soft·Δx, the
                    // Goudsmit-Saunderson first moment of the screened cross-section
                    // restricted to μ > μ_c — computed from the SAME cross-section
                    // as the hard tail (no Highland, no double-count).
                    double g1soft = soft_transport_per_gcm2(
                                        mat, bcp, KE + m_e, kHardScatterMuCut);
                    theta0 = std::sqrt(std::max(0.0, g1soft * step_gcm2));
                } else {
                    // Highland; for TAILS, full-path log (B2) and subtract the hard
                    // tail's projected variance hard_mean·<θ²>/2 so core+tail match
                    // the Highland variance (no over-diffusion → no spurious box escape).
                    double ln_x = std::log((principled && !dbg_no_b2) ? sum_x_X0 : x_X0);
                    theta0 = (13600.0 / bcp) * std::sqrt(x_X0) * (1.0 + 0.038 * ln_x);
                    if (principled && !dbg_no_tails && hard_mean > 0.0 && theta0 > 0.0) {
                        double var_hard = 0.5 * hard_mean * kHardThetaSq;
                        double t0sq = theta0 * theta0 - var_hard;
                        theta0 = (t0sq > 0.0) ? std::sqrt(t0sq) : 0.0;
                    }
                }
                if (theta0 > 1e-10) {
                    double dx    = gauss(rng) * theta0;
                    double dy    = gauss(rng) * theta0;
                    double theta = std::sqrt(dx * dx + dy * dy);
                    if (theta > 1e-10) {
                        Eigen::Vector3d u, v;
                        build_perp_basis(dir, u, v);
                        Eigen::Vector3d dir_new =
                            std::cos(theta) * dir
                            + std::sin(theta) * (dx * u + dy * v) / theta;
                        dir = dir_new.normalized();
                    }
                }
            }

            // Phase D3: large-angle screened-Rutherford single-scatter tail — the
            // backscatter the Gaussian core misses (the physics the empirical
            // albedo gate stands in for). A Poisson number of hard (θ > cutoff)
            // elastic scatters this step, each from the screened Rutherford beyond
            // the cutoff; can fully reverse the direction. Principled only.
            if (hard_mean > 0.0) {
                std::poisson_distribution<int> pois(hard_mean);
                int n_hard = pois(rng);
                double eta = moliere_screening_eta(effective_Z(mat), bcp);
                for (int h = 0; h < n_hard; ++h) {
                    double mu = sample_hard_mu(eta, kHardScatterMuCut,
                                               uniform(rng));
                    double th = std::acos(std::max(-1.0, std::min(1.0, mu)));
                    double ph = 2.0 * M_PI * uniform(rng);
                    Eigen::Vector3d u, v;
                    build_perp_basis(dir, u, v);
                    dir = (std::cos(th) * dir
                           + std::sin(th) * (std::cos(ph) * u
                                             + std::sin(ph) * v)).normalized();
                }
            }
        }

        steps_remaining--;
    }

    return result;  // stopped inside the source geometry
}


// ============================================================
// Free functions
// ============================================================

Eigen::Vector3d sample_photoelectron_direction(
    const Eigen::Vector3d& photon_direction,
    double KE_keV,
    std::mt19937_64& rng)
{
    // Electron Lorentz factor and velocity
    double gamma_e = 1.0 + KE_keV / kMeC2_keV;
    double beta_e  = std::sqrt(1.0 - 1.0 / (gamma_e * gamma_e));

    // For low-energy electrons (β ≪ 1), direction is nearly isotropic (sin²θ pattern).
    // The Sauter algorithm still works but can be slow due to low acceptance.
    // Below 10 keV all electrons deposit locally (range < 1 μm), so direction
    // doesn't affect efficiency; sample isotropically for speed.
    thread_local std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double cos_theta;

    if (beta_e < 0.01) {
        // Non-relativistic: sin²θ distribution (Sauter in the β→0 limit)
        // Peak at θ = 90° (cos θ = 0). Sample by rejection.
        do {
            cos_theta = 2.0 * uniform(rng) - 1.0;
        } while (uniform(rng) > (1.0 - cos_theta * cos_theta));
    } else {
        // Relativistic Sauter distribution:
        //   f(u) ∝ (1 − u²) / (1 − β u)⁴
        //
        // Proposal: p(u) ∝ 1 / (1 − β u)²
        //   Inverse CDF: u = (ξ(1+β) − 1) / (1 + ξ β (1+β))
        // Acceptance probability: (1 − β²)(1 − u²) / (1 − β u)²  ∈ [0, 1]
        //   (Maximum = 1 at u = β, verified analytically)
        double beta = beta_e;
        double one_minus_beta2 = 1.0 - beta * beta;

        // Correct inverse CDF for p(u) ∝ 1/(1 − βu)²:
        //   CDF(u) = (1−β²)/2β × [1/(1−βu) − 1/(1+β)]
        //   Set CDF = ξ and solve: u = (β + 2ξ − 1) / (1 + β(2ξ − 1))
        // Verification: ξ=0 → u=−1, ξ=0.5 → u=β (median), ξ=1 → u=+1  ✓
        do {
            double xi = uniform(rng);
            double xi2 = 2.0 * xi - 1.0;  // ∈ [−1, +1]
            cos_theta = (beta + xi2) / (1.0 + beta * xi2);
            // Clamp to [-1, 1] for numerical safety
            cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
            // Acceptance probability: (1−β²)(1−u²)/(1−βu)² ∈ [0, 1]
            // Maximum = 1 at u=β  (verified analytically)
            double one_minus_bv = 1.0 - beta * cos_theta;
            double accept = one_minus_beta2 * (1.0 - cos_theta * cos_theta)
                            / (one_minus_bv * one_minus_bv);
            if (uniform(rng) <= accept) break;
        } while (true);
    }

    double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    double phi = kTwoPi * uniform(rng);

    // Build orthonormal basis around the photon direction
    Eigen::Vector3d u, v;
    build_perp_basis(photon_direction, u, v);

    return (sin_theta * std::cos(phi) * u
          + sin_theta * std::sin(phi) * v
          + cos_theta * photon_direction).normalized();
}


Eigen::Vector3d compton_recoil_direction(
    const Eigen::Vector3d& photon_dir_in,
    const Eigen::Vector3d& photon_dir_out,
    double energy_in_keV,
    double energy_out_keV)
{
    // From 4-momentum conservation: p_e = p_γ_in − p_γ_out
    // |p_γ| = E/c, so in units with c=1:  p_e ∝ E_in × d_in − E_out × d_out
    Eigen::Vector3d p_e = energy_in_keV * photon_dir_in
                        - energy_out_keV * photon_dir_out;

    double mag = p_e.norm();
    if (mag < 1e-30) {
        // Degenerate case (E_in ≈ E_out and d_in ≈ d_out): electron barely recoils.
        // Return the photon direction as a fallback.
        return photon_dir_in;
    }
    return p_e / mag;
}

} // namespace ceelo
