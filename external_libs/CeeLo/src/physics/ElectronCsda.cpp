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
#include "physics/BremsAngular.h"
#include "physics/estar_stopping_data.h"
#include "cross_sections/CrossSectionData.h"
#include "cross_sections/electron_elastic_data.h"
#include "cross_sections/electron_ionization_data.h"
#include "geometry/Geometry.h"
#include "geometry/SourceGeometry.h"

#include <Eigen/Geometry>

#include <cmath>
#include <cassert>
#include <algorithm>
#include <stdexcept>
#include <random>
#include <cstdlib>
#include <deque>
#include <string>
#include <limits>
#include <utility>

namespace ceelo {

static_assert(estar_data::kElementCount == static_cast<std::size_t>(kMaxElectronTableZ),
              "ESTAR tables must cover exactly the electron-table domain");

namespace {

// ============================================================
// Electron-side element domain
//
// The electron tables (ESTAR stopping, ICRU-49 I, Seltzer-Berger brems, and the
// range / brems-integral tables built from them) stop at uranium. Every element
// Z > kMaxElectronTableZ is treated as uranium throughout the electron physics:
// its tables AND the atomic weight, Z/A, Z(Z+1)/A weights, radiation length and
// effective Z used with them, so the per-material electron view stays
// self-consistent. electron_table_z() documents the size of the approximation.
// ============================================================

/// Map a photon-domain Z to the electron tables; throws for a Z no table covers
/// (this replaces the former silent I=200 eV / A=2Z / 0 / 1e30 fallbacks).
int checked_electron_z(int Z) {
    if (Z < 1 || Z > kMaxZ)
        throw std::out_of_range("ElectronCsda: Z=" + std::to_string(Z)
                                + " is outside 1.." + std::to_string(kMaxZ));
    return electron_table_z(Z);
}

/// Row index into a per-element electron table. Ze must already be mapped
/// through electron_table_z(); an unmapped Z fails here in Debug builds.
inline std::size_t electron_row(int Ze) {
    assert(Ze >= 1 && Ze <= kMaxElectronTableZ);
    return static_cast<std::size_t>(Ze - 1);
}

// ============================================================
// Physical constants
// ============================================================
constexpr double kMeC2_keV = 510.9989461;   // Electron rest energy (keV)
constexpr double kMeC2_MeV = 0.5109989461;  // Electron rest energy (MeV)
constexpr double kPi   = 3.14159265358979323846;
constexpr double kTwoPi = 2.0 * kPi;
constexpr double kFineStructure = 7.2973525693e-3;
constexpr double kBremsAvogadro = 6.02214076e23;
constexpr double kMillibarn_cm2 = 1.0e-27;

/// Exact capped Poisson draw optimized for the small optical depths used by
/// condensed charged-particle steps. libc++'s general Poisson sampler has
/// appreciable setup/rejection overhead when almost every draw is zero.
int sample_capped_poisson(
    double mean, int cap, std::mt19937_64& rng)
{
    if (!(mean > 0.0) || cap <= 0)
        return 0;
    if (mean >= 10.0) {
        std::poisson_distribution<int> poisson(mean);
        return std::min(poisson(rng), cap);
    }

    thread_local std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const double target = uniform(rng);
    double probability = std::exp(-mean);
    double cumulative = probability;
    int count = 0;
    while (target > cumulative && count < cap) {
        ++count;
        probability *= mean / static_cast<double>(count);
        cumulative += probability;
    }
    return count;
}

/// Positron/electron Seltzer-Berger differential-cross-section ratio used by
/// G4LivermoreBremsstrahlungModel.  It tends to one at high energy/soft photon
/// emission and suppresses hard positron bremsstrahlung.
double positron_brems_correction(int Z, double electron_KE_keV,
                                 double photon_energy_keV)
{
    if (!(electron_KE_keV > 0.0) || !(photon_energy_keV > 0.0)
        || photon_energy_keV >= electron_KE_keV)
        return 0.0;
    const auto inverse_beta = [](double kinetic_keV) {
        return (kinetic_keV + kMeC2_keV)
            / std::sqrt(kinetic_keV
                        * (kinetic_keV + 2.0 * kMeC2_keV));
    };
    const double exponent =
        kTwoPi * kFineStructure * static_cast<double>(Z)
        * (inverse_beta(electron_KE_keV)
           - inverse_beta(electron_KE_keV - photon_energy_keV));
    if (exponent < -12.0)
        return 0.0;
    return std::exp(exponent);
}

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
    for (const auto& c : mat.composition()) {
        const int Ze = electron_table_z(c.Z);
        ne += rho * NA * c.mass_fraction * Ze / ElectronCsda::atomic_weight(Ze);
    }
    return ne;
}

void build_perp_basis(const Eigen::Vector3d& d,
                      Eigen::Vector3d& u,
                      Eigen::Vector3d& v);

/// Sample the kinetic-energy fraction carried by a Møller/Bhabha delta
/// electron.  This is a unit-preserving port of
/// G4MollerBhabhaModel::SampleSecondaries.
double sample_delta_energy_keV(double primary_KE_keV,
                               double minimum_delta_keV,
                               bool primary_is_positron,
                               std::mt19937_64& rng)
{
    const double maximum_delta_keV = primary_is_positron
        ? primary_KE_keV : 0.5 * primary_KE_keV;
    if (!(minimum_delta_keV < maximum_delta_keV))
        return 0.0;

    const double xmin = minimum_delta_keV / primary_KE_keV;
    const double xmax = maximum_delta_keV / primary_KE_keV;
    const double gamma = 1.0 + primary_KE_keV / kMeC2_keV;
    const double gamma2 = gamma * gamma;
    const double beta2 = 1.0 - 1.0 / gamma2;
    thread_local std::uniform_real_distribution<double> uniform(0.0, 1.0);

    double x = 0.0;
    if (!primary_is_positron) {
        const double gg = (2.0 * gamma - 1.0) / gamma2;
        const double ymax = 1.0 - xmax;
        const double majorant =
            1.0 - gg * xmax
            + xmax * xmax
                * (1.0 - gg
                   + (1.0 - gg * ymax) / (ymax * ymax));
        for (int attempt = 0; attempt < 1000; ++attempt) {
            const double proposal = uniform(rng);
            x = xmin * xmax
                / (xmin * (1.0 - proposal)
                   + xmax * proposal);
            const double y = 1.0 - x;
            const double rejection =
                1.0 - gg * x
                + x * x
                    * (1.0 - gg + (1.0 - gg * y) / (y * y));
            if (majorant * uniform(rng) <= rejection)
                return x * primary_KE_keV;
        }
    } else {
        const double y = 1.0 / (1.0 + gamma);
        const double y2 = y * y;
        const double y12 = 1.0 - 2.0 * y;
        const double b1 = 2.0 - y2;
        const double b2 = y12 * (3.0 + y2);
        const double y122 = y12 * y12;
        const double b4 = y122 * y12;
        const double b3 = b4 + y122;
        const double xmax2 = xmax * xmax;
        const double majorant =
            1.0
            + (xmax2 * xmax2 * b4
               - xmin * xmin * xmin * b3
               + xmax2 * b2 - xmin * b1) * beta2;
        for (int attempt = 0; attempt < 1000; ++attempt) {
            const double proposal = uniform(rng);
            x = xmin * xmax
                / (xmin * (1.0 - proposal)
                   + xmax * proposal);
            const double x2 = x * x;
            const double rejection =
                1.0
                + (x2 * x2 * b4 - x * x2 * b3
                   + x2 * b2 - x * b1) * beta2;
            if (majorant * uniform(rng) <= rejection)
                return x * primary_KE_keV;
        }
    }
    return 0.0;
}

Eigen::Vector3d sample_delta_direction(
    const Eigen::Vector3d& primary_direction,
    double primary_KE_keV,
    double delta_KE_keV,
    std::mt19937_64& rng)
{
    const double delta_momentum =
        std::sqrt(delta_KE_keV
                  * (delta_KE_keV + 2.0 * kMeC2_keV));
    const double primary_momentum =
        std::sqrt(primary_KE_keV
                  * (primary_KE_keV + 2.0 * kMeC2_keV));
    double cosine =
        delta_KE_keV * (primary_KE_keV + 2.0 * kMeC2_keV)
        / (delta_momentum * primary_momentum);
    cosine = std::clamp(cosine, -1.0, 1.0);
    const double sine =
        std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
    thread_local std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const double phi = kTwoPi * uniform(rng);
    Eigen::Vector3d u, v;
    build_perp_basis(primary_direction, u, v);
    return (
        cosine * primary_direction
        + sine * (std::cos(phi) * u + std::sin(phi) * v)).normalized();
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
    CompoundBremsChi(const CrossSectionData& xs, const Material& mat,
                     double KE_keV, bool is_positron = false)
        : xs_(xs)
        , eb_(xs.sb_chi_energy_bracket(KE_keV))
        , KE_keV_(KE_keV)
        , is_positron_(is_positron)
    {
        const auto& comps = mat.composition();
        ncomp_ = comps.size();
        assert(ncomp_ <= kMaxComp);
        for (size_t i = 0; i < ncomp_; ++i) {
            const int Z = electron_table_z(comps[i].Z);
            z_[i]     = Z;
            // The tabulated chi carries the absolute per-atom normalization after
            // multiplication by Z²/beta².  Number density per gram supplies
            // mass_fraction/A; Z(Z+1)/A is a radiation-length approximation
            // and gives the wrong material trend, especially at low Z.
            w_rad_[i] = comps[i].mass_fraction * Z * Z
                        / ElectronCsda::atomic_weight(Z);
            w_sum_   += w_rad_[i];
        }
    }

    double operator()(double kappa) const {
        const auto kb = xs_.sb_chi_kappa_bracket(kappa);
        return evaluate(kappa, kb, /*apply_positron_correction=*/true);
    }

    double maximum(double kappa_min) const {
        // The positron Elwert factor is <= 1, so the electron chi envelope is
        // also a rigorous positron envelope. Avoid evaluating the expensive
        // correction at every grid node.
        double answer = evaluate(
            kappa_min, xs_.sb_chi_kappa_bracket(kappa_min),
            /*apply_positron_correction=*/false);
        for (int i = 0; i < kSB_n_kappa; ++i) {
            if (kSB_kappa[i] < kappa_min)
                continue;
            const CrossSectionData::SBChiKBracket kb =
                i == kSB_n_kappa - 1
                    ? CrossSectionData::SBChiKBracket{
                          kSB_n_kappa - 2, 1.0}
                    : CrossSectionData::SBChiKBracket{i, 0.0};
            answer = std::max(
                answer,
                evaluate(
                    kSB_kappa[i], kb,
                    /*apply_positron_correction=*/false));
        }
        // Protect against a shallow interpolation maximum between grid nodes.
        return 1.02 * answer;
    }

private:
    double evaluate(
        double kappa,
        const CrossSectionData::SBChiKBracket& kb,
        bool apply_positron_correction) const
    {
        double chi_val = 0.0;
        for (size_t i = 0; i < ncomp_; ++i) {
            double value = xs_.sb_chi_bracketed(z_[i], eb_, kb);
            if (is_positron_ && apply_positron_correction) {
                value *= positron_brems_correction(
                    z_[i], KE_keV_, kappa * KE_keV_);
            }
            chi_val += w_rad_[i] * value;
        }
        return (w_sum_ > 0.0) ? chi_val / w_sum_ : 1.0;
    }

    static constexpr size_t kMaxComp = 32;
    const CrossSectionData& xs_;
    CrossSectionData::SBChiEBracket eb_;
    std::array<int, kMaxComp>    z_{};
    std::array<double, kMaxComp> w_rad_{};
    size_t ncomp_ = 0;
    double w_sum_ = 0.0;
    double KE_keV_ = 0.0;
    bool is_positron_ = false;
};

// ============================================================
// ICRU Report 49 (1993) Table 5.1: Mean excitation energies I(Z) in eV.
// Used in the Bethe-Bloch collision stopping power formula.
// Source: ICRU 49 / NIST PSTAR database.
// ============================================================
static const double kICRU49_I_eV[kMaxElectronTableZ] = {
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

// Standard atomic weights (g/mol) for Z = 1..92, the electron-table domain.
// Values from IUPAC 2021 Table of Atomic Weights (rounded).
static const double kAtomicWeight[kMaxElectronTableZ] = {
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

/// Mean forward displacement divided by true condensed-history path length
/// for a Gaussian multiple-scattering step. Highland theta0 is the RMS of
/// either projected angle. Treating angular diffusion as Brownian through the
/// step gives <cos(theta(s))> = exp(-theta0² s/L), hence
///   <Delta z>/L = (1 - exp(-theta0²))/theta0².
///
/// Advancing by the full true path along the entry direction makes thin source
/// geometries leak electrons too early. The correction tends exactly to one
/// for a vanishing step and remains positive for large angular diffusion.
double highland_forward_displacement_fraction(double KE_keV,
                                              double step_gcm2,
                                              double X0_gcm2,
                                              double cumulative_x_X0)
{
    if (!(KE_keV > 0.0) || !(step_gcm2 > 0.0)
        || !(X0_gcm2 > 0.0))
        return 1.0;
    const double x_X0 = step_gcm2 / X0_gcm2;
    if (!(x_X0 > 0.0))
        return 1.0;
    // beta*p*c = T(T+2m)/(T+m), not p*c = sqrt(T(T+2m)).
    const double beta_pc =
        KE_keV * (KE_keV + 2.0 * kMeC2_keV)
        / (KE_keV + kMeC2_keV);
    if (!(beta_pc > 1e-12))
        return 1.0;
    // Highland's logarithmic correction describes the accumulated slab, not
    // each artificial condensed-history substep independently.  Reusing the
    // per-substep x/X0 underestimates both angular diffusion and detour, with
    // the largest bias for thin high-Z shields.
    const double log_path =
        std::max(cumulative_x_X0, x_X0);
    const double highland_log = 1.0 + 0.038 * std::log(log_path);
    if (!(highland_log > 0.0))
        return 1.0;
    const double theta0 =
        (13600.0 / beta_pc) * std::sqrt(x_X0) * highland_log;
    const double variance = theta0 * theta0;
    if (variance < 1e-8)
        return 1.0 - 0.5 * variance;
    return std::clamp(-std::expm1(-variance) / variance, 0.05, 1.0);
}

// sample_brems_direction_2bs lives in physics/BremsAngular.{h,cpp} (hoisted
// so the next-event estimator's angular EVALUATOR shares the same rejection
// function as this sampler; RNG stream unchanged). The unqualified calls
// below resolve to ceelo::sample_brems_direction_2bs from that header.

/// Apply momentum recoil after a brems photon is emitted.
void apply_brems_recoil(Eigen::Vector3d& electron_direction,
                        double electron_kinetic_energy_keV,
                        const Eigen::Vector3d& photon_direction,
                        double photon_energy_keV)
{
    const double momentum = std::sqrt(std::max(
        electron_kinetic_energy_keV
            * (electron_kinetic_energy_keV + 2.0 * kMeC2_keV),
        0.0));
    const Eigen::Vector3d outgoing_momentum =
        momentum * electron_direction
        - photon_energy_keV * photon_direction;
    if (outgoing_momentum.squaredNorm() > 1e-20)
        electron_direction = outgoing_momentum.normalized();
}

/// One Seltzer-Berger emission conditional on an emission occurring.
struct SampledBremsstrahlung {
    double energy_keV = 0.0;
    Eigen::Vector3d direction = Eigen::Vector3d::UnitZ();
};

SampledBremsstrahlung sample_single_bremsstrahlung(
    const Material& mat,
    Eigen::Vector3d& electron_direction,
    double available_keV,
    bool is_positron,
    std::mt19937_64& rng)
{
    constexpr double kPhotonCut_keV = 10.0;
    if (!(available_keV > kPhotonCut_keV))
        return {};

    thread_local std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const auto& xs_data = CrossSectionData::instance();
    const double kappa_min = kPhotonCut_keV / available_keV;
    CompoundBremsChi compound_chi(
        xs_data, mat, available_keV, is_positron);
    double chi_max = compound_chi.maximum(kappa_min);
    if (!(chi_max > 0.0))
        chi_max = 1.0;

    double photon_energy_keV = 0.0;
    for (int attempt = 0; attempt < 100; ++attempt) {
        const double proposed_keV =
            kPhotonCut_keV
            * std::pow(
                available_keV / kPhotonCut_keV, uniform(rng));
        const double kappa = proposed_keV / available_keV;
        if (uniform(rng) * chi_max < compound_chi(kappa)) {
            photon_energy_keV =
                std::min(proposed_keV, available_keV);
            break;
        }
    }
    if (!(photon_energy_keV > 0.0))
        return {};

    const Eigen::Vector3d photon_direction =
        sample_brems_direction_2bs(
            electron_direction, available_keV, photon_energy_keV,
            ElectronCsda::effective_Z(mat), rng);
    apply_brems_recoil(
        electron_direction, available_keV,
        photon_direction, photon_energy_keV);
    return {photon_energy_keV, photon_direction};
}

/// Sample every discrete bremsstrahlung emission in one condensed-history
/// step. A Bernoulli "one or zero" draw has expectation 1-exp(-lambda), which
/// undercounts the physical Poisson mean lambda by approximately lambda/2.
/// Sampling the multiplicity restores the exact tabulated macroscopic emission
/// rate while retaining CeeLo's inexpensive fixed charged-particle substeps.
template<typename PhotonSink>
double sample_step_bremsstrahlung(
    const ElectronCsda& csda,
    const Material& mat,
    Eigen::Vector3d& electron_direction,
    double electron_kinetic_energy_keV,
    double step_gcm2,
    bool is_positron,
    std::mt19937_64& rng,
    PhotonSink&& photon_sink)
{
    constexpr double kPhotonCut_keV = 10.0;
    if (!(electron_kinetic_energy_keV > kPhotonCut_keV)
        || !(step_gcm2 > 0.0))
        return 0.0;

    const double macro_xs = csda.brems_emission_rate_per_gcm2(
        mat, electron_kinetic_energy_keV, is_positron);
    const double mean_multiplicity = macro_xs * step_gcm2;
    if (!(mean_multiplicity > 0.0))
        return 0.0;

    // The ordinary source/detector substeps have lambda << 1.  This guard
    // prevents pathological work if the routine is ever called with a
    // grossly oversized step; energy conservation normally terminates first.
    const int multiplicity =
        sample_capped_poisson(mean_multiplicity, 64, rng);
    if (multiplicity == 0)
        return 0.0;

    double total_brems_keV = 0.0;

    for (int emission = 0; emission < multiplicity; ++emission) {
        const double available_keV =
            electron_kinetic_energy_keV - total_brems_keV;
        if (!(available_keV > kPhotonCut_keV))
            break;

        const SampledBremsstrahlung photon =
            sample_single_bremsstrahlung(
                mat, electron_direction, available_keV, is_positron, rng);
        if (!(photon.energy_keV > 0.0))
            continue;

        photon_sink(photon.direction, photon.energy_keV);
        total_brems_keV += photon.energy_keV;
    }

    return total_brems_keV;
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
    if (KE_keV <= 0.0) return 0.0;
    const int Ze = checked_electron_z(Z);
    const double stopping = estar_data::radiative_stopping(Ze, KE_keV);
    const double energy_MeV = KE_keV * 1.0e-3;
    const double A = ElectronCsda::atomic_weight(Ze);
    const double X0 = ElectronCsda::radiation_length_gcm2_element(Ze, A);
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
        const int Ze = electron_table_z(comp.Z);
        double A  = ElectronCsda::atomic_weight(Ze);
        double X0_i = ElectronCsda::radiation_length_gcm2_element(Ze, A);
        double weight = comp.mass_fraction / X0_i;  // contribution to 1/X0
        double corr_i = radiative_correction_factor(Ze, KE_keV);
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
    return kICRU49_I_eV[electron_row(checked_electron_z(Z))];
}

double ElectronCsda::atomic_weight(int Z)
{
    return kAtomicWeight[electron_row(checked_electron_z(Z))];
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
    if (KE_keV <= 0.0 || A_g_mol <= 0.0) return 1e30;
    Z = checked_electron_z(Z);  // Z 93-98 use uranium's ESTAR row and I

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

    // For each element Z = 1..92 (the electron-table domain), build the CSDA range
    // table via trapezoid integration
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
        for (int Zm1 = 0; Zm1 < kMaxElectronTableZ; ++Zm1) {
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

    // Build per-element electron Seltzer-Berger brems normalization integrals.
    // The more expensive positron Elwert-corrected counterparts are initialized
    // lazily, one element at a time, by ensure_positron_brems_integral().
    // For each element Z and grid energy E, integrate chi(kappa) and chi(kappa)/kappa
    // over kappa = [kBremsThreshold/E, 1] using the SAME 64-point log-spaced
    // trapezoid scheme used previously at runtime, so the compound result matches
    // the old per-substep computation exactly at grid energies.
    constexpr double kBremsThreshold_keV = 10.0;  // matches deposited_in_scoring
    constexpr int    N_int = 64;
    const auto& xs_data = CrossSectionData::instance();
    for (int Zm1 = 0; Zm1 < kMaxElectronTableZ; ++Zm1) {
        int Z = Zm1 + 1;
        for (int i = 0; i < kNGrid; ++i) {
            double E = energy_grid_keV_[i];
            if (E <= kBremsThreshold_keV) {
                sb_Jchi_[Zm1][i]   = 0.0;
                sb_Jchiok_[Zm1][i] = 0.0;
                sb_Jchiok_pos_[Zm1][i] = 0.0;
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
            sb_Jchiok_pos_[Zm1][i] = 0.0;
        }
    }
}

void ElectronCsda::ensure_positron_brems_integral(int Z) const
{
    const std::size_t row = electron_row(Z);  // Z is already an electron-table Z
    std::call_once(sb_Jchiok_pos_once_[row], [this, Z, row] {
        constexpr double kBremsThreshold_keV = 10.0;
        constexpr int N_int = 64;
        const auto& xs_data = CrossSectionData::instance();
        auto& table = sb_Jchiok_pos_[row];
        for (int i = 0; i < kNGrid; ++i) {
            const double E = energy_grid_keV_[i];
            if (E <= kBremsThreshold_keV) {
                table[i] = 0.0;
                continue;
            }
            const double kappa_min = kBremsThreshold_keV / E;
            const double log_kmin = std::log(kappa_min);
            const double dlog = -log_kmin / N_int;
            double previous_kappa = kappa_min;
            double previous_value =
                xs_data.sb_chi(Z, E, previous_kappa)
                * positron_brems_correction(
                    Z, E, previous_kappa * E);
            double integral = 0.0;
            for (int j = 1; j <= N_int; ++j) {
                const double kappa =
                    std::exp(log_kmin + j * dlog);
                const double value =
                    xs_data.sb_chi(Z, E, kappa)
                    * positron_brems_correction(Z, E, kappa * E);
                const double dk = kappa - previous_kappa;
                integral +=
                    0.5
                    * (previous_value / previous_kappa
                       + value / kappa)
                    * dk;
                previous_kappa = kappa;
                previous_value = value;
            }
            table[i] = integral;
        }
    });
}

const ElectronCsda& ElectronCsda::instance()
{
    static ElectronCsda singleton;
    return singleton;
}


// ============================================================
// Range interpolation
// ============================================================

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
        const int Ze = electron_table_z(comp.Z);
        double A_i = atomic_weight(Ze);
        if (A_i <= 0.0) continue;
        double w_rad = comp.mass_fraction * Ze * (Ze + 1.0) / A_i;
        I_chi   += w_rad * interpolate_log_grid(sb_Jchi_[electron_row(Ze)], KE_keV);
        I_chiok += w_rad * interpolate_log_grid(sb_Jchiok_[electron_row(Ze)], KE_keV);
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

double ElectronCsda::brems_emission_rate_per_gcm2(
    const Material& mat, double KE_keV, bool is_positron) const
{
    constexpr double kPhotonCut_keV = 10.0;
    if (!(KE_keV > kPhotonCut_keV))
        return 0.0;

    // G4LivermoreBremsstrahlungModel:
    //   sigma_atom = Z² millibarn / beta²
    //              * integral chi(T,kappa) d(log kappa).
    // Multiplication by N_A * mass_fraction/A converts the atomic cross
    // section to the requested macroscopic cross section per g/cm².
    const double total_energy = KE_keV + kMeC2_keV;
    const double inverse_beta_squared =
        total_energy * total_energy
        / (KE_keV * (KE_keV + 2.0 * kMeC2_keV));
    double rate = 0.0;
    for (const auto& comp : mat.composition()) {
        const int Z = electron_table_z(comp.Z);
        const double A = atomic_weight(Z);
        if (!(A > 0.0))
            continue;
        if (is_positron)
            ensure_positron_brems_integral(Z);
        const auto& integral_table =
            is_positron ? sb_Jchiok_pos_[electron_row(Z)]
                        : sb_Jchiok_[electron_row(Z)];
        const double integral =
            interpolate_log_grid(integral_table, KE_keV);
        rate += comp.mass_fraction * kBremsAvogadro / A
              * static_cast<double>(Z * Z) * kMillibarn_cm2
              * inverse_beta_squared * integral;
    }
    return std::max(0.0, rate);
}

double ElectronCsda::delta_electron_rate_per_gcm2(
    const Material& mat, double KE_keV, double minimum_delta_keV,
    bool is_positron) const
{
    const double maximum_delta_keV =
        is_positron ? KE_keV : 0.5 * KE_keV;
    if (!(minimum_delta_keV > 0.0)
        || !(minimum_delta_keV < maximum_delta_keV))
        return 0.0;

    const double xmin = minimum_delta_keV / KE_keV;
    const double xmax = maximum_delta_keV / KE_keV;
    const double tau = KE_keV / kMeC2_keV;
    const double gamma = tau + 1.0;
    const double gamma2 = gamma * gamma;
    const double beta2 = tau * (tau + 2.0) / gamma2;
    double reduced_cross_section = 0.0;

    if (!is_positron) {
        const double gg = (2.0 * gamma - 1.0) / gamma2;
        reduced_cross_section =
            ((xmax - xmin)
                 * (1.0 - gg + 1.0 / (xmin * xmax)
                    + 1.0
                        / ((1.0 - xmin) * (1.0 - xmax)))
             - gg
                 * std::log(
                     xmax * (1.0 - xmin)
                     / (xmin * (1.0 - xmax))))
            / beta2;
    } else {
        const double y = 1.0 / (1.0 + gamma);
        const double y2 = y * y;
        const double y12 = 1.0 - 2.0 * y;
        const double b1 = 2.0 - y2;
        const double b2 = y12 * (3.0 + y2);
        const double y122 = y12 * y12;
        const double b4 = y122 * y12;
        const double b3 = b4 + y122;
        reduced_cross_section =
            (xmax - xmin)
                * (1.0 / (beta2 * xmin * xmax) + b2
                   - 0.5 * b3 * (xmin + xmax)
                   + b4
                       * (xmin * xmin + xmin * xmax + xmax * xmax)
                       / 3.0)
            - b1 * std::log(xmax / xmin);
    }

    // G4's twopi_mc2_rcl2 / T, in cm² per target electron.
    constexpr double kClassicalElectronRadiusSq_cm2 =
        7.9407876e-26;
    const double cross_section_per_electron_cm2 =
        reduced_cross_section
        * kTwoPi * kMeC2_keV
        * kClassicalElectronRadiusSq_cm2 / KE_keV;
    double electrons_per_g = 0.0;
    for (const auto& comp : mat.composition()) {
        const int Ze = electron_table_z(comp.Z);
        electrons_per_g +=
            kBremsAvogadro * comp.mass_fraction * Ze / atomic_weight(Ze);
    }
    return std::max(
        0.0, cross_section_per_electron_cm2 * electrons_per_g);
}

double ElectronCsda::range_gcm2(
    int Z, double KE_keV, bool is_positron) const
{
    const std::size_t row = electron_row(checked_electron_z(Z));
    return interpolate_log_grid(
        is_positron ? range_table_pos_[row] : range_table_[row], KE_keV);
}

namespace {

/// Marker for a table that depends on nothing beyond the material itself.
struct NoExtra {
    friend bool operator==(NoExtra, NoExtra) { return true; }
};

/// Per-thread memo of a per-(material, extra) table over the energy grid.
///
/// Keyed on Material::cache_key() -- a per-OBJECT identity minted by every
/// constructor and assignment -- and NEVER on the Material address.  An address
/// aliases across a destroyed-and-rebuilt Material: same address, same density,
/// same first element gives a stale hit that silently serves a DIFFERENT
/// material's table, and a partial composition signature (count, first element,
/// density) narrows that window without closing it.  See Material::cache_key().
///
/// Route every new per-material grid table through here.  The keying rule then
/// lives in one place instead of being re-derived -- correctly or not -- at each
/// new cache, which is exactly how two of these caches were built wrong.
///
/// `extra` carries whatever else the table depends on (positron flag, delta-ray
/// cut, ...) and is compared with operator==.  `compute` fills a table that is
/// only inserted once it returns, so a throw leaves nothing cached (the
/// hand-rolled versions inserted first and filled after, which would have cached
/// a half-built table) and a table may still be built incrementally from its own
/// earlier elements.  `compute` must not recurse into the SAME memo -- it would
/// not find its own entry and would compute a duplicate.  Recursing into a
/// DIFFERENT one is fine and is what the depth tables do.
///
/// Each call site instantiates its own cache (the closure type is part of the
/// instantiation) -- one memo per table, which is the intent.
///
/// std::deque, not std::vector: these return a reference into the container and
/// callers hold it across further calls, so growth must not invalidate.  The old
/// hand-rolled caches papered over this with reserve(32) and a comment.
template <class Table, class Extra, class Compute>
const Table& material_grid_memo(const Material& mat, const Extra& extra,
                                Compute&& compute)
{
    struct Entry {
        uint64_t key;
        Extra extra;
        Table table;
    };
    // Worker threads are spawned per compute(), so the cache lifetime is
    // naturally bounded and no locking is needed.
    static thread_local std::deque<Entry> cache;

    const uint64_t key = mat.cache_key();
    for (auto& e : cache) {
        if (e.key == key && e.extra == extra) return e.table;
    }

    Table table{};
    compute(table);
    cache.push_back(Entry{key, extra, std::move(table)});
    return cache.back().table;
}

} // namespace

const std::array<double, ElectronCsda::kNGrid>&
ElectronCsda::compound_range_table(const Material& mat, bool is_positron) const
{
    return material_grid_memo<std::array<double, kNGrid>>(
        mat, is_positron,
        [&](std::array<double, kNGrid>& R) {
            const auto& comp = mat.composition();
            const auto& src_table = is_positron ? range_table_pos_ : range_table_;

            // Bragg additivity for ranges at each grid energy:
            //   1/R_compound = Σᵢ wᵢ / Rᵢ.
            for (int i = 0; i < kNGrid; ++i) {
                double inv_range_sum = 0.0;
                for (const auto& c : comp) {
                    double Ri = src_table[electron_row(electron_table_z(c.Z))][i];
                    if (Ri > 1e-30) inv_range_sum += c.mass_fraction / Ri;
                }
                R[i] = (inv_range_sum > 0.0) ? 1.0 / inv_range_sum : 1e30;
            }
        });
}

const std::array<double, ElectronCsda::kNGrid>&
ElectronCsda::compound_brems_depth_table(
    const Material& mat, bool is_positron) const
{
    return material_grid_memo<std::array<double, kNGrid>>(
        mat, is_positron,
        [&](std::array<double, kNGrid>& depth) {
            const auto& ranges = compound_range_table(mat, is_positron);
            double previous_rate =
                brems_emission_rate_per_gcm2(
                    mat, energy_grid_keV_[0], is_positron);
            depth[0] = previous_rate * ranges[0];
            for (int i = 1; i < kNGrid; ++i) {
                const double rate =
                    brems_emission_rate_per_gcm2(
                        mat, energy_grid_keV_[i], is_positron);
                const double delta_range =
                    std::max(0.0, ranges[i] - ranges[i - 1]);
                depth[i] =
                    depth[i - 1]
                    + 0.5 * (previous_rate + rate) * delta_range;
                previous_rate = rate;
            }
        });
}

const std::array<double, ElectronCsda::kNGrid>&
ElectronCsda::compound_delta_depth_table(
    const Material& mat,
    double minimum_delta_keV,
    bool is_positron) const
{
    // minimum_delta_keV is a genuine parameter of the table, not an identity
    // proxy, so it stays in the key alongside the material.
    return material_grid_memo<std::array<double, kNGrid>>(
        mat, std::make_pair(minimum_delta_keV, is_positron),
        [&](std::array<double, kNGrid>& depth) {
            const auto& ranges = compound_range_table(mat, is_positron);
            double previous_rate =
                delta_electron_rate_per_gcm2(
                    mat, energy_grid_keV_[0], minimum_delta_keV, is_positron);
            depth[0] = previous_rate * ranges[0];
            for (int i = 1; i < kNGrid; ++i) {
                const double rate =
                    delta_electron_rate_per_gcm2(
                        mat, energy_grid_keV_[i],
                        minimum_delta_keV, is_positron);
                const double delta_range =
                    std::max(0.0, ranges[i] - ranges[i - 1]);
                depth[i] =
                    depth[i - 1]
                    + 0.5 * (previous_rate + rate) * delta_range;
                previous_rate = rate;
            }
        });
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
    // Deterministic in (material, KE), so the material alone keys it.
    return material_grid_memo<std::array<double, kNGrid>>(
        mat, NoExtra{},
        [&](std::array<double, kNGrid>& C) {
            for (int i = 0; i < kNGrid; ++i) {
                C[i] = radiative_correction_compound(mat, energy_grid_keV_[i]);
            }
        });
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
    // Z is taken on the electron side (Z 93-98 -> U), matching the atomic weight
    // callers pass (ElectronCsda::atomic_weight is the uranium proxy there too).
    if (A_g_mol <= 0.0) return 1e30;
    Z = checked_electron_z(Z);
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
        const int Ze = electron_table_z(comp.Z);
        double A  = atomic_weight(Ze);
        double X0 = radiation_length_gcm2_element(Ze, A);
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
        const int Ze = electron_table_z(c.Z);
        double Z = static_cast<double>(Ze);
        double A = ElectronCsda::atomic_weight(Ze);
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
        const int Ze = electron_table_z(c.Z);
        double Z = static_cast<double>(Ze);
        double A = ElectronCsda::atomic_weight(Ze);
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
        num += c.mass_fraction * static_cast<double>(electron_table_z(c.Z));
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
        if (remaining_KE >= 0.01) {  // ran out of segments still moving
            result.escaped = true;
            result.exit_KE_keV = remaining_KE;
            result.exit_direction = direction;
        }
        return result;
    }

    // ----------------------------------------------------------------
    // High energy (KE > 200 keV): Molière condensed-history walk.
    // Adaptive sub-steps (more at higher energy for finer brems sampling);
    // each step applies Highland multiple scattering and may emit a discrete
    // bremsstrahlung photon at the absolute Seltzer-Berger rate.
    // Energy loss per step is
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
    double sum_x_X0 = 0.0;  // accumulated thickness for the Highland log term
    const Material* straggle_mat = nullptr;   // per-walk Z/A cache
    double straggle_ZoverA = 0.0;
    std::normal_distribution<double> gauss(0.0, 1.0);
    thread_local std::uniform_real_distribution<double> uniform(0.0, 1.0);

    int iter_guard = 20 * n_steps;  // hard cap: degenerate-geometry backstop
    bool left_geometry = false;
    while (steps_remaining > 0 && KE > 0.01 && iter_guard-- > 0) {
        geometry.trace_ray(pos, dir, segments);
        if (segments.empty()) { left_geometry = true; break; }

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
            bool advanced = false;
            for (const auto& seg : segments) {
                if (seg.t_start > 1e-8) {
                    pos += dir * (seg.t_start + 1e-9);
                    advanced = true;
                    break;
                }
            }
            if (!advanced) { left_geometry = true; break; }
            continue;
        }

        if (!cur->material) { left_geometry = true; break; }

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

        const Eigen::Vector3d step_start = pos;
        const Eigen::Vector3d step_direction = dir;
        const double step_cm = step_gcm2 / density;

        // The charged particle follows the entry direction during this
        // condensed-history step. Brems creation is sampled within this
        // segment below; it must not be placed systematically at the far end.
        pos += step_direction * step_cm;

        // Push across boundary so next trace_ray finds the next volume
        if (reach_boundary) {
            pos += dir * 1e-6;
        }

        // Radiation length for this material (needed for brems and scattering)
        double X0 = radiation_length_gcm2(mat);

        // --- Bremsstrahlung emission (absolute Seltzer-Berger rate) ---
        double brems_keV = 0.0;
        if (X0 > 0.0 && !disable_brems && KE > kBremsThreshold_keV) {
            // Emission point is sampled uniformly WITHIN the sub-step rather
            // than pinned at its far end, and the photon direction comes from
            // BremsAngular rather than being taken as the electron direction.
            // The zero-opening-angle alternative measures null against a
            // G4ModifiedTsai A/B (FEP -0.12% +- 0.13 @2614 keV, +0.08% +- 0.14
            // @3000 keV) -- the crystal is large compared to the ~m_e/E cone --
            // so this is a refinement, not a defect fix.
            brems_keV = sample_step_bremsstrahlung(
                *this, mat, dir, KE, step_gcm2, is_positron, rng,
                [&](const Eigen::Vector3d& photon_direction,
                    double photon_energy_keV) {
                    const Eigen::Vector3d emission_position =
                        step_start
                        + step_direction * (uniform(rng) * step_cm);
                    result.brems_photons.push_back(
                        {emission_position, photon_direction,
                         photon_energy_keV});
                });
        }

        // CSDA collision loss (computed from current KE, which includes collision
        // stopping only — brems is a separate radiative process).
        double KE_exit  = residual_energy_keV(mat, KE, step_gcm2, is_positron);
        double delta_KE = KE - KE_exit;
        if (delta_KE < 0.0) delta_KE = 0.0;

        // Per-step Gaussian energy-loss straggling (Bohr variance, same formula
        // as the principled source walk).  Contained electrons still deposit
        // their full KE — the walk continues until KE -> 0 (see the budget
        // guard below), so straggling moves only stop/escape endpoints, never
        // the contained-event total.  Without it the escape threshold is a
        // deterministic step function of birth depth, which under-predicts FEP
        // at >=2 MeV (measured +0.23% @2614, +0.50% @3000 on config 1).
        // The FINAL budgeted step is left deterministic: it is planned as the
        // full remaining range, so straggling it would leave residual KE ~half
        // the time and force costly walk extensions (measured +13% CPU @3 MeV)
        // for ~1/n_steps of the variance.  CEELO_NO_CRYSTAL_STRAGGLE=1
        // disables straggling entirely (A/B and debugging).
        static const bool kNoCrystalStraggle =
            std::getenv("CEELO_NO_CRYSTAL_STRAGGLE") != nullptr;
        static const bool kStraggleAllKE =
            std::getenv("CEELO_DIAG_STRAGGLE_ALL_KE") != nullptr;
        if (!kNoCrystalStraggle && steps_remaining > 1 && step_gcm2 > 1e-12 &&
            (KE > 500.0 || kStraggleAllKE)) {
            double beta2  = KE * (KE + 2.0 * m_e) / ((KE + m_e) * (KE + m_e));
            if (&mat != straggle_mat) {   // per-walk cache: material rarely changes
                straggle_ZoverA = 0.0;
                for (const auto& c : mat.composition()) {
                    const int Ze = electron_table_z(c.Z);
                    straggle_ZoverA += c.mass_fraction * Ze / atomic_weight(Ze);
                }
                straggle_mat = &mat;
            }
            double sigma2 = kBohrXi_MeV2cm2g * straggle_ZoverA * step_gcm2
                          * (1.0 - 0.5 * beta2);                 // MeV²
            if (sigma2 > 0.0) {
                delta_KE += gauss(rng) * std::sqrt(sigma2) * 1000.0;  // keV
                if (delta_KE < 0.0)  delta_KE = 0.0;
                if (delta_KE > KE)   delta_KE = KE;
            }
        }

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
        // The log term is evaluated on the ACCUMULATED path thickness, not the
        // per-sub-step thickness: the (1 + 0.038 ln x) correction is not
        // additive in variance, and applying it to x/n_steps under-scatters the
        // summed walk by ~0.038·ln(n_steps) (~13% in θ₀ for 20 steps).  Same
        // convention as the B2 source-geometry walk above.
        if (X0 > 0.0) {
            double x_X0  = step_gcm2 / X0;
            sum_x_X0    += x_X0;
            double bcp   = std::sqrt(KE * (KE + 2.0 * m_e));  // β×p×c in keV
            if (bcp > 1e-6 && x_X0 > 0.0) {
                double ln_x   = std::log(sum_x_X0);
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
        // A walk that still carries KE when the budget empties must keep
        // walking (the budget is a discretization aid, not a physics limit);
        // exiting here silently drops the residual KE.  With straggling the
        // final full-range step under-deposits ~half the time, so this guard
        // is what makes per-step straggling energy-conserving for contained
        // electrons.
        if (steps_remaining == 0 && KE > 0.01) steps_remaining = 1;
    }

    // Terminal position of the condensed-history walk: the rest point if the
    // electron stopped (KE → 0), otherwise the last boundary it reached before
    // leaving all geometry.  Used as the positron annihilation site by callers.
    result.stop_position = pos;
    if (left_geometry && KE > 0.01) {
        result.escaped = true;
        result.exit_KE_keV = KE;
        result.exit_direction = dir;
    }
    return result;
}


// ============================================================
// Detailed electron transport (spectrum fidelity)
// ============================================================
//
// A mixed (class II) elastic-scattering scheme in the PENELOPE mould, for every
// electron energy (no 200 keV split):
//   * Elastic scattering: per material, a Wentzel distribution
//       p(mu') ∝ 1/(A + mu')²,  mu' = (1 - cos theta)/2,
//     whose screening A and normalization reproduce the first two transport
//     cross sections (sigma_1, sigma_2) of the EPICS2023 EEDL elastic data
//     (electron_elastic_data.cpp) -- so the angular diffusion AND the
//     large-angle tail that backscatters electrons out of a surface are both
//     carried.  Collisions deflecting by more than mu'_c are sampled one by
//     one ("hard"), with mean free path lambda_h = max(lambda_W,
//     min(C1 lambda_1, C2 R)); the softer ones over each step become one
//     Gaussian deflection applied at a uniformly sampled point of the step
//     (the "random hinge"), with the soft first transport moment as variance.
//   * Energy loss (class II as well): collisions transferring more than
//     W_cc = 2 keV are sampled one by one from the Moller cross section, and
//     each delta ray is walked like any electron, so the ones born near a
//     surface can leave; the primary takes the binary-collision recoil. The
//     soft remainder (the ESTAR CSDA loss less the hard collisions' mean) is
//     drawn from PENELOPE's artificial distribution with the soft variance.
//     This straggling is what puts the early, energetic backscatter where
//     GEANT4 has it. Seltzer-Berger brems via sample_step_bremsstrahlung,
//     in-flight annihilation for positrons.
//   * Boundaries: a flight that reaches a region boundary stops there and a
//     new step is sampled beyond it (memoryless hard events make this exact).
//   * Fast path: an electron whose residual CSDA range is below
//     Geometry::safety_distance() cannot leave its region, so it deposits its
//     energy there at once (brems from the contained sampler).  Most
//     electrons are born deep inside a crystal, so they never walk.
// Positrons use the electron elastic moments (EEDL has none for positrons) and
// the Moller collisions (Bhabha differs only in the hard tail).
namespace {

/// <mu'> and <mu'^2> of the normalized Wentzel distribution with screening A.
struct WentzelMoments { double m1, m2; };
WentzelMoments wentzel_moments(double A)
{
    const double L = std::log1p(1.0 / A);
    const double K = A * (1.0 + A);
    return {K * (L + A / (1.0 + A) - 1.0),
            K * (1.0 - 2.0 * A * L + A / (1.0 + A))};
}

/// sigma_2 / sigma_1 of the Wentzel distribution: 3 as A -> 0, 1 as A -> inf.
double wentzel_ratio21(double A)
{
    const WentzelMoments m = wentzel_moments(A);
    return 3.0 * (1.0 - m.m2 / m.m1);
}

constexpr double kWentzelA_min = 1e-12;
constexpr double kWentzelA_max = 1e3;

/// The screening A whose Wentzel distribution has sigma_2/sigma_1 = R.
double fit_wentzel_screening(double R)
{
    double lo = std::log(kWentzelA_min), hi = std::log(kWentzelA_max);
    if (R >= wentzel_ratio21(kWentzelA_min)) return kWentzelA_min;
    if (R <= wentzel_ratio21(kWentzelA_max)) return kWentzelA_max;
    for (int it = 0; it < 80; ++it) {      // ratio decreases with A
        const double mid = 0.5 * (lo + hi);
        if (wentzel_ratio21(std::exp(mid)) > R) lo = mid; else hi = mid;
    }
    return std::exp(0.5 * (lo + hi));
}

/// Per-material elastic table on the EEDL energy grid: ln A and ln sigma_1
/// (per g/cm²), both linear in ln E between nodes.
struct ElasticTable {
    std::array<double, eedl_elastic::kEnergyCount> ln_A{};
    std::array<double, eedl_elastic::kEnergyCount> ln_G1{};
    double ln_E0 = 0.0, d_lnE = 1.0;
};

const ElasticTable& compound_elastic_table(const Material& mat)
{
    return material_grid_memo<ElasticTable>(mat, 0, [&](ElasticTable& t) {
        using namespace eedl_elastic;
        constexpr double kNA = 6.02214076e23;
        t.ln_E0 = std::log(kEnergy_keV[0]);
        t.d_lnE = std::log(kEnergy_keV[1]) - t.ln_E0;
        for (std::size_t k = 0; k < kEnergyCount; ++k) {
            // Compound transport moments are additive.
            double s1 = 0.0, s2 = 0.0;
            for (const auto& c : mat.composition()) {
                const int Ze = electron_table_z(c.Z);
                const double atoms_per_g =
                    kNA * c.mass_fraction / ElectronCsda::atomic_weight(Ze);
                const double sig1 = 1e-24 * kSigma1_barn[Ze - 1][k];
                s1 += atoms_per_g * sig1;
                s2 += atoms_per_g * sig1 * kRatio21[Ze - 1][k];
            }
            t.ln_G1[k] = std::log(std::max(s1, 1e-300));
            t.ln_A[k]  = std::log(fit_wentzel_screening(s1 > 0.0 ? s2 / s1 : 2.0));
        }
    });
}

/// Step parameters of the mixed scheme at kinetic energy KE.
struct ElasticStep {
    double lambda_h_g = 0.0;  ///< hard mean free path (g/cm²)
    double A = 1.0;           ///< Wentzel screening
    double P_hard = 1.0;      ///< fraction of collisions that are hard
    double G1_soft = 0.0;     ///< soft first transport coefficient (1/(g/cm²))
};

ElasticStep elastic_step(const ElasticTable& t, double KE_keV, double range_g,
                         double C1, double C2)
{
    const double x = (std::log(KE_keV) - t.ln_E0) / t.d_lnE;
    const int n = static_cast<int>(eedl_elastic::kEnergyCount);
    const int k = std::clamp(static_cast<int>(std::floor(x)), 0, n - 2);
    const double f = std::clamp(x - k, 0.0, 1.0);
    ElasticStep e;
    e.A = std::exp(t.ln_A[k] + f * (t.ln_A[k + 1] - t.ln_A[k]));
    const double G1 = std::exp(t.ln_G1[k] + f * (t.ln_G1[k + 1] - t.ln_G1[k]));
    const double A = e.A;
    const double sigma_W = G1 / (2.0 * wentzel_moments(A).m1);   // collisions per g/cm²
    const double lambda_W = 1.0 / sigma_W;
    e.lambda_h_g = std::max(lambda_W, std::min(C1 / G1, C2 * range_g));
    e.P_hard = std::min(1.0, lambda_W / e.lambda_h_g);
    const double mu_c = A * (1.0 - e.P_hard) / (A + e.P_hard);
    // First transport moment of the hard part: 2 sigma_W int_{mu_c}^1 mu' p(mu').
    const double G1_hard = 2.0 * sigma_W * A * (1.0 + A)
        * (std::log((1.0 + A) / (A + mu_c)) + A / (1.0 + A) - A / (A + mu_c));
    e.G1_soft = std::max(0.0, G1 - G1_hard);
    return e;
}

/// Turn `dir` by polar angle with cosine `cos_theta` at a uniform azimuth.
void deflect(Eigen::Vector3d& dir, double cos_theta, double phi)
{
    cos_theta = std::clamp(cos_theta, -1.0, 1.0);
    const double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    Eigen::Vector3d u, v;
    build_perp_basis(dir, u, v);
    dir = (cos_theta * dir + sin_theta * (std::cos(phi) * u + std::sin(phi) * v))
              .normalized();
}

double env_or(const char* name, double fallback)
{
    const char* text = std::getenv(name);
    return text ? std::atof(text) : fallback;
}

/// Moller collisions transferring W in [W_cc, T/2] (free-electron, per electron
/// of the medium): their count and summed transfer per g/cm², and a sampler.
struct MollerHard {
    double T = 0.0, W_cc = 0.0, k = 0.0;   // k = 2 pi r_e^2 m c^2 / beta^2 (keV cm²)
    double n_per_g = 0.0;                   // hard collisions per g/cm²
    double mean_W_per_g = 0.0;              // energy they carry per g/cm²
    double soft_var_per_g = 0.0;            // variance of the soft (W < W_cc) loss per g/cm²
};

MollerHard moller_hard(double T, double W_cc, double electrons_per_g)
{
    constexpr double m = 510.998950;
    constexpr double k0 = 2.549549e-22;   // 2 pi r_e^2 m c^2, keV cm²
    MollerHard h;
    h.T = T;
    h.W_cc = W_cc;
    const double beta2 = T * (T + 2.0 * m) / ((T + m) * (T + m));
    h.k = k0 / beta2;
    const double a = W_cc, b = 0.5 * T;
    h.soft_var_per_g = electrons_per_g * h.k * W_cc;
    if (b <= a) return h;
    // Integrals of the Moller DCS from a to b = T/2 (count, then W-weighted).
    const double Tm2 = (T + m) * (T + m);
    const double c4 = (2.0 * T + m) * m / Tm2;
    const double ln_ta = std::log((T - a) / a);   // ln((T-a)/a)
    const double ln_hb = std::log(b / (T - a));    // ln(b/(T-a)) <= 0
    const double n = 1.0 / a - 1.0 / (T - a) + (b - a) / Tm2 - c4 / T * ln_ta;
    const double w = (ln_ta + ln_hb) + (2.0 - T / (T - a) + ln_hb)
                   + 0.5 * (b * b - a * a) / Tm2 + c4 * ln_hb;
    h.n_per_g = electrons_per_g * h.k * n;
    h.mean_W_per_g = electrons_per_g * h.k * w;
    return h;
}

/// A non-negative soft energy loss with the given mean and variance
/// (PENELOPE's artificial distribution: truncated Gaussian, uniform, or a
/// zero/uniform mixture, as the spread grows against the mean).
double sample_soft_loss(double mean, double var, std::mt19937_64& rng)
{
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    if (!(var > 0.0) || !(mean > 0.0)) return std::max(0.0, mean);
    const double mean2 = mean * mean;
    if (mean2 > 9.0 * var) {
        std::normal_distribution<double> gauss(mean, std::sqrt(var));
        for (;;) {
            const double x = gauss(rng);
            if (x > 0.0 && x < 2.0 * mean) return x;
        }
    }
    if (mean2 > 3.0 * var) {
        const double half = std::sqrt(3.0 * var);
        return mean - half + 2.0 * half * uniform(rng);
    }
    // Zero with probability 1 - p, else uniform on (0, w0): mean p w0 / 2,
    // second moment p w0^2 / 3.
    const double w0 = 1.5 * (var + mean2) / mean;
    const double p = 2.0 * mean / w0;
    return uniform(rng) < p ? w0 * uniform(rng) : 0.0;
}

/// A transfer W from the restricted Moller distribution (1/W² proposal).
double sample_moller_W(const MollerHard& h, std::mt19937_64& rng)
{
    constexpr double m = 510.998950;
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const double T = h.T, a = h.W_cc, b = 0.5 * T;
    const double Tm2 = (T + m) * (T + m);
    const double c4 = (2.0 * T + m) * m / Tm2;
    for (;;) {
        const double W = 1.0 / (1.0 / a - uniform(rng) * (1.0 / a - 1.0 / b));
        const double g = 1.0 + W * W / ((T - W) * (T - W)) + W * W / Tm2 - c4 * W / (T - W);
        if (uniform(rng) * 2.25 <= g) return W;
    }
}

/// Per-material electron-impact K-shell ionization table: each element's
/// rate weight (atoms/g x 1e-24 cm²/barn x omega_K, so weight x sigma_K in
/// barn is K x-rays per g/cm²) and its CSDA K x-ray yield from each energy
/// of a log grid down to rest.
struct ImpactTable {
    static constexpr int kN = 96;
    struct Element {
        int Z = 0;
        double weight = 0.0;
        double threshold_keV = 0.0;  ///< EEDL's K-shell ionization threshold
        const FluorescenceData* fl = nullptr;
        std::array<double, kN> yield{};
    };
    std::vector<Element> elements;
    double ln_E0 = 0.0, d_lnE = 1.0;
};

const ImpactTable& compound_impact_table(const Material& mat, const ElectronCsda& csda)
{
    return material_grid_memo<ImpactTable>(mat, 0, [&](ImpactTable& t) {
        constexpr double kNA = 6.02214076e23;
        constexpr double kEmin_keV = 1.0, kEmax_keV = 20000.0;
        t.ln_E0 = std::log(kEmin_keV);
        t.d_lnE = std::log(kEmax_keV / kEmin_keV) / (ImpactTable::kN - 1);
        for (const auto& c : mat.composition()) {
            if (c.Z < 1 || c.Z > static_cast<int>(eedl_ionization::kElementCount)) continue;
            const FluorescenceData* fl = CrossSectionData::instance().fluorescence(c.Z);
            if (!fl || fl->num_lines == 0 || !(fl->fluorescence_yield > 0.0f)) continue;
            ImpactTable::Element e;
            e.Z = c.Z;
            e.threshold_keV = eedl_ionization::kKShell[c.Z - 1].energy_keV[0];
            e.fl = fl;
            e.weight = kNA * c.mass_fraction / ElectronCsda::atomic_weight(c.Z)
                     * 1e-24 * fl->fluorescence_yield;
            t.elements.push_back(e);
        }
        // Yield = integral of weight x sigma_K over the CSDA path, with
        // dR = R(E_b) - R(E_a) from the material's range table: midpoint rule on
        // eight linear sub-steps per grid interval, from rest upward.
        constexpr int kSub = 8;
        std::vector<double> acc(t.elements.size(), 0.0);
        double E_prev = 0.0, R_prev = 0.0;
        for (int k = 0; k < ImpactTable::kN; ++k) {
            const double E_k = std::exp(t.ln_E0 + k * t.d_lnE);
            double ea = E_prev, Ra = R_prev;
            for (int j = 1; j <= kSub; ++j) {
                const double eb = E_prev + (E_k - E_prev) * j / kSub;
                const double Rb = csda.range_gcm2_material(mat, eb);
                const double em = 0.5 * (ea + eb);
                for (std::size_t i = 0; i < t.elements.size(); ++i)
                    acc[i] += t.elements[i].weight
                            * eedl_ionization::k_shell_sigma_barn(t.elements[i].Z, em)
                            * std::max(0.0, Rb - Ra);
                ea = eb;
                Ra = Rb;
            }
            E_prev = E_k;
            R_prev = Ra;
            for (std::size_t i = 0; i < t.elements.size(); ++i)
                t.elements[i].yield[k] = acc[i];
        }
    });
}

/// An isotropic direction for an impact-ionization x-ray.
Eigen::Vector3d impact_xray_direction(std::mt19937_64& rng)
{
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const double c = 2.0 * uniform(rng) - 1.0;
    const double s = std::sqrt(std::max(0.0, 1.0 - c * c));
    const double phi = 6.283185307179586 * uniform(rng);
    return {s * std::cos(phi), s * std::sin(phi), c};
}

/// A yield column interpolated linearly in ln E (linearly to zero below 1 keV).
double impact_yield_at(const ImpactTable& t, const std::array<double, ImpactTable::kN>& y,
                       double E_keV)
{
    if (!(E_keV > 0.0)) return 0.0;
    const double x = (std::log(E_keV) - t.ln_E0) / t.d_lnE;
    if (x <= 0.0) return y[0] * E_keV / std::exp(t.ln_E0);
    if (x >= ImpactTable::kN - 1) return y[ImpactTable::kN - 1];
    const int i = static_cast<int>(x);
    const double f = x - i;
    return y[i] + f * (y[i + 1] - y[i]);
}

/// An element's slowing-down yield, zero at and below its K threshold: the log
/// grid would otherwise interpolate across the threshold and let an electron
/// just below it ionize the shell.
double element_impact_yield_at(const ImpactTable& t, const ImpactTable::Element& e,
                               double E_keV)
{
    return E_keV > e.threshold_keV ? impact_yield_at(t, e.yield, E_keV) : 0.0;
}

}  // namespace

double ElectronCsda::impact_xray_rate_per_gcm2(const Material& mat, double KE_keV) const
{
    const ImpactTable& t = compound_impact_table(mat, *this);
    double rate = 0.0;
    for (const auto& e : t.elements)
        rate += e.weight * eedl_ionization::k_shell_sigma_barn(e.Z, KE_keV);
    return rate;
}

double ElectronCsda::impact_xray_yield(const Material& mat, double KE_keV) const
{
    const ImpactTable& t = compound_impact_table(mat, *this);
    double yield = 0.0;
    for (const auto& e : t.elements) yield += element_impact_yield_at(t, e, KE_keV);
    return yield;
}

void ElectronCsda::sample_impact_xrays(const Material& mat, double KE_keV, double mean,
                                       bool slowing_down, std::mt19937_64& rng,
                                       std::vector<ImpactXray>& out) const
{
    if (!(mean > 0.0)) return;
    std::poisson_distribution<int> count(mean);
    int n = count(rng);
    if (n == 0) return;
    const ImpactTable& t = compound_impact_table(mat, *this);
    if (t.elements.empty()) return;
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    thread_local std::vector<double> share;
    share.assign(t.elements.size(), 0.0);
    double total = 0.0;
    for (std::size_t i = 0; i < t.elements.size(); ++i) {
        const auto& e = t.elements[i];
        share[i] = slowing_down
            ? element_impact_yield_at(t, e, KE_keV)
            : e.weight * eedl_ionization::k_shell_sigma_barn(e.Z, KE_keV);
        total += share[i];
    }
    if (!(total > 0.0)) return;
    for (; n > 0; --n) {
        double r = uniform(rng) * total;
        std::size_t i = 0;
        while (i + 1 < t.elements.size() && r >= share[i]) r -= share[i++];
        const FluorescenceData& fl = *t.elements[i].fl;
        const double xi = uniform(rng);
        double cum = 0.0;
        int chosen = 0;
        for (int li = 0; li < fl.num_lines; ++li) {
            cum += fl.line_probability[li];
            if (xi <= cum) { chosen = li; break; }
        }
        out.push_back({t.elements[i].Z, static_cast<double>(fl.line_energy_keV[chosen]),
                       fl.line_filling_subshell[chosen]});
    }
}

double ElectronCsda::impact_ionization_max_keV()
{
    static const double emax =
        env_or("CEELO_EWALK_EII_EMAX", std::numeric_limits<double>::infinity());
    return emax;
}

ElectronDepositResult ElectronCsda::deposited_in_scoring_detailed(
                                           const Eigen::Vector3d& position,
                                           const Eigen::Vector3d& direction,
                                           double KE_keV,
                                           const Geometry& geometry,
                                           std::mt19937_64& rng,
                                           bool disable_brems,
                                           bool is_positron,
                                           double impact_max_keV) const
{
    ElectronDepositResult result;
    result.stop_position = position;
    if (KE_keV < 0.01) return result;

    // Electron-impact K-shell ionization (impact_max_keV > 0): an ionization
    // tags part of the energy the stopping power already removes, so its
    // x-ray is emitted like a brems photon (BremsPhoton::impact_Z) without
    // changing the particle's energy; one born in a scoring region is paid for
    // out of the walk's scoring deposit at the end (`impact_debt`).
    const bool impact = impact_max_keV > 0.0;
    double impact_debt = 0.0;
    thread_local std::vector<ImpactXray> impact_buf;
    auto emit_impact_xrays = [&](const Eigen::Vector3d& from, const Eigen::Vector3d& along,
                                 double length_cm, bool scoring) {
        for (const ImpactXray& x : impact_buf) {
            std::uniform_real_distribution<double> u(0.0, 1.0);
            BremsPhoton photon{from + along * (u(rng) * length_cm),
                               impact_xray_direction(rng), x.line_keV};
            photon.impact_Z = static_cast<std::uint8_t>(x.Z);
            photon.impact_filling_subshell = x.filling_subshell;
            result.brems_photons.push_back(photon);
            if (scoring) impact_debt += x.line_keV;
        }
        impact_buf.clear();
    };

    // Step controls (PENELOPE's C1, C2) and the absorption energy; the
    // environment overrides exist for convergence A/B studies. Electron-beam
    // spectra against GEANT4 are unchanged between 0.01, 0.05 and 0.1 (bands
    // within ~2%); 0.2 starts to drift (full-energy +1-2%). 0.1 halves the
    // steps of 0.05.
    static const double kC1 = env_or("CEELO_EWALK_C1", 0.1);
    static const double kC2 = env_or("CEELO_EWALK_C2", 0.1);
    static const double kAbsorb_keV = env_or("CEELO_EWALK_EABS", 1.0);
    static const bool kNoFastPath = std::getenv("CEELO_EWALK_NO_FASTPATH") != nullptr;
    // Hard-collision threshold W_cc (keV); 0 turns the hard collisions off
    // (CSDA loss, Bohr straggling above 500 keV, no delta rays). Electron-beam
    // spectra against GEANT4 are the same at 1 and 2 keV; 2 is cheaper.
    static const double kMollerWcc = env_or("CEELO_EWALK_WCC", 2.0);
    constexpr double kTwoPi_ = 6.283185307179586;
    constexpr double m_e = 510.998950;

    thread_local std::vector<PathSegment> segments;
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::normal_distribution<double> gauss(0.0, 1.0);

    // Delta rays from hard collisions wait here; each is walked after the
    // particle that made it, starting in the region it was born in.
    struct Delta {
        Eigen::Vector3d pos, dir;
        double KE;
        const Material* mat;
        bool scoring;
    };
    thread_local std::vector<Delta> deltas;
    deltas.clear();

    Eigen::Vector3d pos = position;
    Eigen::Vector3d dir = direction;
    double KE = KE_keV;
    bool positron = is_positron;   // the primary only; delta rays are electrons
    bool primary = true;
    bool annihilated = false;

    // State of the current step: 0 = sample a new one, 1 = flying to the
    // hinge, 2 = flying from the hinge to the step end.
    int phase = 0;
    double flight_left_g = 0.0, after_hinge_g = 0.0, soft_mean = 0.0;
    bool step_hard = false;
    ElasticStep step;
    const Material* last_mat = nullptr;
    const ElasticTable* elastic = nullptr;
    double ZoverA = 0.0;
    bool left_geometry = false;

    // The region the particle is in, from the last ray trace. A trace is the
    // walk's dominant cost, so it is redone only when the region is unknown
    // (start, after a boundary or a gap) or the next flight might reach a
    // boundary -- a flight shorter than the safety distance cannot.
    bool need_region = true;
    const Material* region_mat = nullptr;
    bool region_scoring = false;
    double d_bound_cm = 0.0;
    bool d_bound_valid = false;
    auto locate = [&]() -> const PathSegment* {
        geometry.trace_ray(pos, dir, segments);
        for (const auto& seg : segments)
            if (seg.t_start < 1e-8 || (seg.t_end > 1e-8 && seg.t_start <= 1e-8))
                return &seg;
        return nullptr;
    };

    for (;;) {   // the primary, then its delta rays
    for (int guard = 0; guard < 200000; ++guard) {
        if (need_region) {
            const PathSegment* cur = locate();
            if (segments.empty()) { left_geometry = true; break; }
            if (!cur) {
                // A gap (bore, vacuum between layers): fly straight across it.
                bool advanced = false;
                for (const auto& seg : segments) {
                    if (seg.t_start > 1e-8) {
                        pos += dir * (seg.t_start + 1e-9);
                        advanced = true;
                        break;
                    }
                }
                if (!advanced) { left_geometry = true; break; }
                continue;
            }
            if (!cur->material) { left_geometry = true; break; }
            region_mat = cur->material;
            region_scoring = cur->is_scoring;
            d_bound_cm = cur->t_end - std::max(0.0, cur->t_start);
            d_bound_valid = true;
            need_region = false;
        }

        const Material& mat = *region_mat;
        const double rho = mat.density();
        if (KE <= kAbsorb_keV || !(rho > 0.0)) {
            if (region_scoring) result.deposited_scoring_keV += KE;
            KE = 0.0;
            break;
        }
        if (&mat != last_mat) {
            ZoverA = 0.0;
            for (const auto& c : mat.composition()) {
                const int Ze = electron_table_z(c.Z);
                ZoverA += c.mass_fraction * Ze / atomic_weight(Ze);
            }
            elastic = &compound_elastic_table(mat);
            last_mat = &mat;
        }
        const double range_g = range_gcm2_material(mat, KE, positron);
        double safety_cm = -1.0;   // computed on demand at this position

        if (phase == 0) {
            // Contained: cannot reach any boundary, so no walk is needed.
            if (!positron && !kNoFastPath
                && range_g / rho < (safety_cm = geometry.safety_distance(pos))) {
                double brems_E = 0.0;
                if (!disable_brems && KE > 10.0) {
                    for (const auto& bp : sample_brems_in_material(mat, pos, dir, KE, rng)) {
                        brems_E += bp.energy_keV;
                        result.brems_photons.push_back(bp);
                    }
                }
                if (region_scoring)
                    result.deposited_scoring_keV += std::max(0.0, KE - brems_E);
                if (impact) {
                    const double ke = std::min(KE, impact_max_keV);
                    impact_buf.clear();
                    sample_impact_xrays(mat, ke, impact_xray_yield(mat, ke), true, rng, impact_buf);
                    emit_impact_xrays(pos, dir, 0.0, region_scoring);
                }
                KE = 0.0;
                break;
            }
            step = elastic_step(*elastic, KE, range_g, kC1, kC2);
            const double s_hard = -std::log(1.0 - uniform(rng)) * step.lambda_h_g;
            const double s_cap  = kC2 * range_g;
            step_hard = s_hard < s_cap;
            const double s = step_hard ? s_hard : s_cap;
            const double tau = uniform(rng) * s;
            soft_mean = -std::expm1(-step.G1_soft * s);   // <1 - cos theta> of the soft part
            flight_left_g = tau;
            after_hinge_g = s - tau;
            phase = 1;
        }

        // Distance to the region boundary along this flight: from a trace,
        // unless the flight provably stays inside.
        if (!d_bound_valid) {
            if (safety_cm < 0.0) safety_cm = geometry.safety_distance(pos);
            if (flight_left_g / rho < safety_cm) {
                d_bound_cm = std::numeric_limits<double>::infinity();
            } else {
                const PathSegment* cur = locate();
                if (!cur || cur->material != region_mat || cur->is_scoring != region_scoring) {
                    need_region = true;   // sitting on a boundary: re-locate
                    continue;
                }
                d_bound_cm = cur->t_end - std::max(0.0, cur->t_start);
            }
        }
        d_bound_valid = false;   // one flight's worth

        // Fly along dir: to the next event, or to the boundary if nearer.
        const double d_bound_g = d_bound_cm * rho;
        const bool hit_boundary = d_bound_g < flight_left_g;
        const double path_g = hit_boundary ? d_bound_g : flight_left_g;
        const double path_cm = path_g / rho;
        const Eigen::Vector3d flight_start = pos;
        const Eigen::Vector3d flight_dir = dir;

        if (positron && path_cm > 0.0) {
            const double P_annih = electron_number_density(mat)
                                 * positron_annih_xsec_per_electron(KE) * path_cm;
            if (P_annih > 0.0 && uniform(rng) < std::min(P_annih, 1.0)) {
                result.annihilated_in_flight    = true;
                result.residual_KE_at_annih_keV = KE;
                result.stop_position            = pos;
                annihilated = true;
                KE = 0.0;
                break;
            }
        }

        // Bremsstrahlung over the flight (at most the residual range).
        const bool stops = path_g >= range_g;
        const double loss_path_g = stops ? range_g : path_g;
        double brems_keV = 0.0;
        if (!disable_brems && KE > 10.0 && loss_path_g > 0.0) {
            const double loss_path_cm = loss_path_g / rho;
            brems_keV = sample_step_bremsstrahlung(
                *this, mat, dir, KE, loss_path_g, positron, rng,
                [&](const Eigen::Vector3d& photon_direction, double photon_energy_keV) {
                    result.brems_photons.push_back(
                        {flight_start + flight_dir * (uniform(rng) * loss_path_cm),
                         photon_direction, photon_energy_keV});
                });
        }

        // Impact ionization over the flight: the CSDA yield of the rest of the
        // slowing-down if the particle stops in it, else the rate times the path.
        if (impact && loss_path_g > 0.0) {
            const double ke = std::min(KE, impact_max_keV);
            double mean = 0.0;
            if (stops) mean = impact_xray_yield(mat, ke);
            else if (KE <= impact_max_keV) mean = impact_xray_rate_per_gcm2(mat, KE) * loss_path_g;
            impact_buf.clear();
            sample_impact_xrays(mat, ke, mean, stops, rng, impact_buf);
            emit_impact_xrays(flight_start, flight_dir, loss_path_g / rho, region_scoring);
        }

        if (stops) {
            if (region_scoring)
                result.deposited_scoring_keV += std::max(0.0, KE - brems_keV);
            KE = 0.0;
            pos = flight_start + flight_dir * (range_g / rho);
            break;
        }

        // Energy loss. With hard collisions on (W_cc > 0), those transferring
        // more than W_cc are sampled one by one from the Moller cross section
        // and their delta rays walked like any electron (so the ones born
        // near a surface can leave); the soft remainder -- the CSDA loss less
        // the hard collisions' mean -- is drawn with the soft variance from
        // PENELOPE's artificial distribution. Otherwise: the CSDA mean,
        // straggled above 500 keV (Bohr; capped at the mean loss, where the
        // Gaussian form stops meaning anything).
        double dE = KE - residual_energy_keV(mat, KE, path_g, positron);
        double hard_keV = 0.0;
        if (kMollerWcc > 0.0 && dE > 0.0) {
            const MollerHard h = moller_hard(KE, kMollerWcc, 6.02214076e23 * ZoverA);
            dE = sample_soft_loss(dE - h.mean_W_per_g * path_g, h.soft_var_per_g * path_g, rng);
            if (h.n_per_g > 0.0) {
                std::poisson_distribution<int> count(h.n_per_g * path_g);
                for (int i = count(rng); i > 0; --i) {
                    const double W = sample_moller_W(h, rng);
                    if (dE + hard_keV + W > KE) break;
                    hard_keV += W;
                    const double phi = kTwoPi_ * uniform(rng);
                    Eigen::Vector3d delta_dir = flight_dir;
                    deflect(delta_dir,
                            std::sqrt(W * (KE + 2.0 * m_e) / (KE * (W + 2.0 * m_e))), phi);
                    deltas.push_back({flight_start + flight_dir * (uniform(rng) * path_cm),
                                      delta_dir, W, region_mat, region_scoring});
                    // The primary's recoil, applied at the flight end.
                    deflect(dir, std::sqrt((KE - W) * (KE + 2.0 * m_e)
                                           / (KE * (KE - W + 2.0 * m_e))),
                            phi + 0.5 * kTwoPi_);
                }
            }
        } else if (KE > 500.0 && dE > 0.0) {
            const double beta2 = KE * (KE + 2.0 * m_e) / ((KE + m_e) * (KE + m_e));
            const double sigma =
                1000.0 * std::sqrt(kBohrXi_MeV2cm2g * ZoverA * path_g * (1.0 - 0.5 * beta2));
            dE += std::min(sigma, dE) * gauss(rng);
        }
        dE = std::clamp(dE, 0.0, KE - hard_keV);
        double deposit = dE;
        if (brems_keV > KE - hard_keV - dE) deposit = std::max(0.0, KE - hard_keV - brems_keV);
        if (region_scoring) result.deposited_scoring_keV += deposit;
        KE = std::max(0.0, KE - hard_keV - deposit - brems_keV);
        pos = flight_start + flight_dir * path_cm;

        if (hit_boundary) {
            pos += flight_dir * 1e-9;   // into the next region
            phase = 0;
            need_region = true;
            continue;
        }
        if (phase == 1) {
            // The hinge: the step's soft deflection, then on to its end.
            if (soft_mean > 0.0) {
                const double theta = std::sqrt(-2.0 * soft_mean * std::log(1.0 - uniform(rng)));
                deflect(dir, std::cos(std::min(theta, 3.141592653589793)),
                        kTwoPi_ * uniform(rng));
            }
            flight_left_g = after_hinge_g;
            phase = 2;
            continue;
        }
        // Step end: a hard elastic collision, unless the step was capped.
        if (step_hard) {
            const double A = step.A;
            const double xi = (1.0 - step.P_hard) + step.P_hard * uniform(rng);
            const double mu = A * xi / (1.0 + A - xi);
            deflect(dir, 1.0 - 2.0 * mu, kTwoPi_ * uniform(rng));
        }
        phase = 0;
    }

    if (primary) {
        if (!annihilated) result.stop_position = pos;
        if (left_geometry && KE > 0.01) {
            result.escaped = true;
            result.exit_KE_keV = KE;
            result.exit_direction = dir;
        }
    } else if (left_geometry) {
        result.escaped_delta_keV += KE;
    }
    if (deltas.empty()) break;
    const Delta d = deltas.back();
    deltas.pop_back();
    pos = d.pos;
    dir = d.dir;
    KE = d.KE;
    positron = false;
    primary = false;
    phase = 0;
    left_geometry = false;
    need_region = false;
    region_mat = d.mat;
    region_scoring = d.scoring;
    d_bound_valid = false;
    }
    if (impact_debt > 0.0) {
        const double pay = std::min(impact_debt, result.deposited_scoring_keV);
        result.deposited_scoring_keV -= pay;
        result.unpaid_impact_keV = impact_debt - pay;
    }
    return result;
}


std::vector<BremsPhoton> ElectronCsda::sample_brems_in_material(
                                           const Material& mat,
                                           const Eigen::Vector3d& position,
                                           const Eigen::Vector3d& direction,
                                           double KE_keV,
                                           std::mt19937_64& rng,
                                           std::uint8_t charged_generation) const
{
    // Thin wrapper over the emission-recording core: same physics, same RNG
    // draw sequence (position accumulation adds no draws), so this is
    // bit-identical to a creation-point model -- every photon is
    // collapsed back to the creation position, which is what the source-
    // geometry callers assume (electron range << distance to detector).
    const std::vector<BremsEmission> emissions =
        sample_brems_emissions_in_material(
            mat, position, direction, KE_keV, rng,
            /*is_positron=*/false, charged_generation);
    std::vector<BremsPhoton> photons;
    photons.reserve(emissions.size());
    for (const BremsEmission& e : emissions)
        photons.push_back(
            {position, e.photon_dir, e.photon_keV, e.charged_generation});
    return photons;
}


std::vector<BremsEmission> ElectronCsda::sample_brems_emissions_in_material(
                                           const Material& mat,
                                           const Eigen::Vector3d& position,
                                           const Eigen::Vector3d& direction,
                                           double KE_keV,
                                           std::mt19937_64& rng,
                                           bool is_positron,
                                           std::uint8_t charged_generation) const
{
    std::vector<BremsEmission> photons;
    if (KE_keV <= kSourceBremsFloor_keV) return photons;

    double X0 = radiation_length_gcm2(mat);
    if (X0 <= 0.0) return photons;

    constexpr double m_e = 511.0;  // electron rest mass (keV)
    constexpr double kBremsThreshold_keV = 10.0;  // minimum brems photon energy

    // A charged particle known to stop inside this homogeneous material has no
    // boundary state to update.  Its brems and (for the primary generation)
    // delta-electron interactions are independent non-homogeneous Poisson
    // processes along CSDA depth.  Draw only their actual event locations
    // instead of taking 20--60 mostly empty condensed-history steps.
    //
    // B(T) = integral[0,T] sigma_brems(T') dR(T')
    // D(T) = integral[0,T] sigma_delta(T') dR(T')
    //
    // Exponential optical-flight clocks compete as T decreases. A brems
    // energy jump skips the corresponding interval of both optical-depth
    // curves, as it should: no charged-particle path was traversed there.
    constexpr double kDeltaTransportFloor_keV = 50.0;
    constexpr std::uint8_t kMaximumDeltaGeneration = 1;
    const auto& brems_depth =
        compound_brems_depth_table(mat, is_positron);
    const auto& range_table =
        compound_range_table(mat, is_positron);
    const bool transport_deltas =
        charged_generation < kMaximumDeltaGeneration;
    const auto* delta_depth = transport_deltas
        ? &compound_delta_depth_table(
              mat, kDeltaTransportFloor_keV,
              is_positron)
        : nullptr;

    auto energy_at_depth = [&](const auto& depth_table, double target) {
        if (target <= depth_table.front())
            return energy_grid_keV_.front();
        if (target >= depth_table.back())
            return energy_grid_keV_.back();
        const auto upper =
            std::upper_bound(
                depth_table.begin(), depth_table.end(), target);
        const int hi = static_cast<int>(upper - depth_table.begin());
        const int lo = std::max(0, hi - 1);
        const double span = depth_table[hi] - depth_table[lo];
        const double fraction =
            span > 0.0
                ? (target - depth_table[lo]) / span
                : 0.0;
        return std::exp(
            log_energy_grid_[lo]
            + fraction
                * (log_energy_grid_[hi] - log_energy_grid_[lo]));
    };

    Eigen::Vector3d dir = direction;
    Eigen::Vector3d pos = position;
    const double density = mat.density();
    double KE = KE_keV;
    std::exponential_distribution<double> optical_flight(1.0);
    std::normal_distribution<double> gauss(0.0, 1.0);
    double brems_remaining = optical_flight(rng);
    double delta_remaining =
        transport_deltas ? optical_flight(rng) : 0.0;

    // Actual multiplicities are small over the supported 12 MeV range. The
    // cap is only a guard against corrupt tables or a future out-of-range use.
    for (int event = 0;
         event < 256 && KE > kBremsThreshold_keV;
         ++event) {
        const double starting_brems_depth =
            interpolate_log_grid(brems_depth, KE);
        const double target_brems_depth =
            starting_brems_depth - brems_remaining;
        const double brems_event_KE =
            target_brems_depth > brems_depth.front()
                ? std::min(
                      KE,
                      energy_at_depth(
                          brems_depth, target_brems_depth))
                : -1.0;

        double starting_delta_depth = 0.0;
        double delta_event_KE = -1.0;
        if (transport_deltas) {
            starting_delta_depth =
                interpolate_log_grid(*delta_depth, KE);
            const double target_delta_depth =
                starting_delta_depth - delta_remaining;
            if (target_delta_depth > delta_depth->front()) {
                delta_event_KE =
                    std::min(
                        KE,
                        energy_at_depth(
                            *delta_depth, target_delta_depth));
            }
        }

        const bool brems_event = brems_event_KE >= delta_event_KE;
        const double interaction_KE =
            brems_event ? brems_event_KE : delta_event_KE;
        if (!(interaction_KE > kBremsThreshold_keV))
            break;

        const double path_gcm2 =
            std::max(
                0.0,
                interpolate_log_grid(range_table, KE)
                - interpolate_log_grid(range_table, interaction_KE));

        // Direction only affects emitted photons for a contained track. Apply
        // the same adaptive Highland diffusion as the old stepper, but only
        // over path intervals preceding an interaction.
        if (path_gcm2 > 0.0) {
            const double starting_range =
                interpolate_log_grid(range_table, KE);
            const int full_track_steps =
                std::max(
                    20,
                    std::min(60, static_cast<int>(KE / 50.0)));
            const int scatter_intervals =
                std::max(
                    1,
                    static_cast<int>(
                        std::ceil(
                            full_track_steps
                            * path_gcm2
                            / std::max(starting_range, 1e-30))));
            const double interval_gcm2 =
                path_gcm2 / scatter_intervals;
            double interval_KE = KE;
            for (int interval = 0;
                 interval < scatter_intervals
                     && interval_KE > interaction_KE;
                 ++interval) {
                // Advance along the CURRENT diffused direction, then apply
                // the interval's scattering kick (move-then-scatter).
                // Deterministic, so the RNG stream matches the creation-point
                // wrapper draw for draw.
                if (density > 0.0)
                    pos += dir * (interval_gcm2 / density);
                const double exit_KE =
                    interval + 1 == scatter_intervals
                        ? interaction_KE
                        : residual_energy_keV(
                              mat, interval_KE, interval_gcm2,
                              is_positron);
                const double bcp =
                    std::sqrt(
                        std::max(
                            exit_KE * (exit_KE + 2.0 * m_e),
                            0.0));
                const double x_X0 = interval_gcm2 / X0;
                if (bcp > 1e-6 && x_X0 > 0.0) {
                    const double theta0 =
                        (13600.0 / bcp) * std::sqrt(x_X0)
                        * (1.0 + 0.038 * std::log(x_X0));
                    if (theta0 > 1e-10) {
                        const double dx = gauss(rng) * theta0;
                        const double dy = gauss(rng) * theta0;
                        const double theta = std::hypot(dx, dy);
                        if (theta > 1e-10) {
                            Eigen::Vector3d u, v;
                            build_perp_basis(dir, u, v);
                            dir = (
                                std::cos(theta) * dir
                                + std::sin(theta)
                                    * (dx * u + dy * v) / theta)
                                      .normalized();
                        }
                    }
                }
                interval_KE = exit_KE;
            }
        }

        // Both clocks consume only the charged path down to this interaction.
        brems_remaining =
            std::max(
                0.0,
                brems_remaining
                - (starting_brems_depth
                   - interpolate_log_grid(
                       brems_depth, interaction_KE)));
        if (transport_deltas) {
            delta_remaining =
                std::max(
                    0.0,
                    delta_remaining
                    - (starting_delta_depth
                       - interpolate_log_grid(
                           *delta_depth, interaction_KE)));
        }
        KE = interaction_KE;

        if (brems_event) {
            // Capture the electron state BEFORE the emission:
            // sample_single_bremsstrahlung deflects `dir` via the recoil,
            // and the 2BS density the estimator evaluates is about the
            // pre-emission direction.
            const Eigen::Vector3d dir_at_emission = dir;
            const SampledBremsstrahlung photon =
                sample_single_bremsstrahlung(
                    mat, dir, KE, is_positron, rng);
            brems_remaining = optical_flight(rng);
            if (!(photon.energy_keV > 0.0))
                continue;
            photons.push_back(
                {pos, dir_at_emission, photon.direction, photon.energy_keV,
                 KE, charged_generation});
            KE = std::max(0.0, KE - photon.energy_keV);
            continue;
        }

        const double delta_KE_keV =
            sample_delta_energy_keV(
                KE, kDeltaTransportFloor_keV,
                is_positron, rng);
        delta_remaining = optical_flight(rng);
        if (!(delta_KE_keV > 0.0))
            continue;
        const Eigen::Vector3d delta_direction =
            sample_delta_direction(
                dir, KE, delta_KE_keV, rng);
        // Knock-on deltas are electrons regardless of the primary's charge;
        // they start from the primary's current along-walk position.
        auto delta_photons = sample_brems_emissions_in_material(
            mat, pos, delta_direction, delta_KE_keV,
            rng, /*is_positron=*/false, charged_generation + 1);
        for (auto& photon : delta_photons)
            photons.push_back(std::move(photon));
    }

    return photons;
}


ElectronSourceWalkResult ElectronCsda::walk_in_source_geometry(
                                           const SourceGeometry& source,
                                           const Material& birth_material,
                                           const Eigen::Vector3d& position,
                                           const Eigen::Vector3d& direction,
                                           double KE_keV,
                                           std::mt19937_64& rng,
                                           bool is_positron,
                                           std::uint8_t charged_generation) const
{
    ElectronSourceWalkResult result;
    result.stop_position = position;
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
            if (!seg.material || seg.length <= 0.0) continue;  // cavity adds no tau
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
    // Very low energy (KE ≤ 10 keV): straight-line CSDA, no brems.
    //
    // Unlike deposited_in_scoring (which switches to straight lines at
    // 200 keV — the path shape is irrelevant for deposits inside the
    // scoring crystal), the source-escape walk keeps Moliere direction
    // diffusion down to 50 keV: escape through the outer skin of a
    // dense shield is dominated by detour and diffusion. The former 50 keV
    // floor also suppressed observable brems from low-energy radioactive
    // conversion electrons, so the source walk now continues down to the
    // 10 keV Seltzer-Berger/tally floor.
    // ----------------------------------------------------------------
    if (KE <= kSourceBremsFloor_keV) {
        source.trace_source_segments(pos, dir, KE, seg_buf);
        double dist = 0.0;
        const Material* exit_mat = &birth_material;  // exit-Z source for the gate
        for (const auto& seg : seg_buf) {
            if (!seg.material || seg.length <= 0.0) {
                // An unfilled cavity: no energy loss, but it still has to be
                //  travelled, or exit_position lands short by the cavity chord.
                if (seg.length > 0.0) dist += seg.length;
                continue;
            }
            double path_gcm2 = seg.length * seg.material->density();
            double range_g =
                range_gcm2_material(*seg.material, KE, is_positron);
            if (path_gcm2 >= range_g) {
                result.traversed_gcm2 += range_g;
                result.stop_position =
                    pos + dir * (dist + range_g / seg.material->density());
                return result;  // stopped inside
            }
            result.traversed_gcm2 += path_gcm2;
            KE = residual_energy_keV(
                *seg.material, KE, path_gcm2, is_positron);
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
            result.stop_position = result.exit_position;
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

    // Electron-impact K x-rays (SourceGeometry::set_source_impact_ionization),
    // carried among the brems photons with BremsPhoton::impact_Z set.
    const bool impact = source.source_impact_ionization_enabled();
    const double impact_emax = impact_ionization_max_keV();
    auto push_impact_xray = [&](const Eigen::Vector3d& at, const ImpactXray& x) {
        BremsPhoton photon{at, impact_xray_direction(rng), x.line_keV, charged_generation};
        photon.impact_Z = static_cast<std::uint8_t>(x.Z);
        photon.impact_filling_subshell = x.filling_subshell;
        result.brems_photons.push_back(photon);
    };
    std::normal_distribution<double> gauss(0.0, 1.0);
    thread_local std::uniform_real_distribution<double> uniform(0.0, 1.0);

    // Material the particle is currently traversing — the EXIT-Z source for the
    // skin-escape gate when the trace empties (the boundary it reflects off).
    const Material* last_material = &birth_material;

    // Phase D: cumulative path in radiation lengths, for the full-path Highland
    // log term (B2) — the per-step ln(x/X0) under-scatters the cumulative angle.
    double sum_x_X0 = 0.0;

    // A cored source has at most one cavity, but a ray can cross it once per
    //  direction change; this only has to bound the loop.
    int void_crossings = 64;

    while (steps_remaining > 0 && KE > 0.01) {
        // Only the first (current) segment is needed here -- max_segments=1 stops
        // the trace after it, avoiding a full multi-shell re-trace every substep.
        source.trace_source_segments(pos, dir, KE, seg_buf, /*max_segments=*/1);

        // A null-material segment is an unfilled cavity INSIDE a cored source,
        //  not the world outside it.  Cross it - vacuum costs no energy and
        //  scatters nothing - and pick the walk up at the far wall.  Treating it
        //  as an exit declared electrons "escaped" at a point deep inside the
        //  source, with full residual KE, skipping the far shell wall and every
        //  shield beyond it.
        if (!seg_buf.empty() && !seg_buf[0].material && seg_buf[0].length > 1e-10) {
            pos += dir * (seg_buf[0].length + 1e-7);
            if (--void_crossings < 0) return result;   // cannot loop forever
            continue;
        }

        if (seg_buf.empty() || !seg_buf[0].material || seg_buf[0].length <= 1e-10) {
            // Outside the source geometry — escaped (subject to the skin-escape
            // gate; brems already emitted along the walk is retained in result).
            // KE is the residual exit energy; gate on the exit-material Z.
            if (escape_survives(KE, effective_Z(*last_material))) {
                result.escaped = true;
                result.exit_position = pos;
                result.exit_direction = dir;
                result.exit_KE_keV = KE;
                result.stop_position = pos;
            }
            return result;
        }

        const Material& mat = *seg_buf[0].material;
        last_material = &mat;
        double density   = mat.density();
        double seg_gcm2  = seg_buf[0].length * density;
        double range_g =
            range_gcm2_material(mat, KE, is_positron);

        double planned_step_gcm2 = range_g / static_cast<double>(steps_remaining);
        double step_gcm2 = planned_step_gcm2;
        bool reach_boundary = false;
        if (planned_step_gcm2 >= seg_gcm2) {
            step_gcm2 = seg_gcm2;
            reach_boundary = true;
        }
        if (step_gcm2 <= 0.0) { steps_remaining--; continue; }
        result.traversed_gcm2 += step_gcm2;

        if (is_positron) {
            const double path_cm = step_gcm2 / density;
            const double annihilation_probability =
                electron_number_density(mat)
                * positron_annih_xsec_per_electron(KE)
                * path_cm;
            if (annihilation_probability > 0.0
                && uniform(rng)
                    < std::min(annihilation_probability, 1.0)) {
                result.annihilated_in_flight = true;
                result.residual_KE_at_annih_keV = KE;
                result.stop_position = pos;
                return result;
            }
        }

        double X0 = radiation_length_gcm2(mat);
        const Eigen::Vector3d step_start = pos;
        const Eigen::Vector3d step_direction = dir;
        const double step_cm = step_gcm2 / density;
        const double cumulative_x_X0 =
            sum_x_X0 + ((X0 > 0.0) ? step_gcm2 / X0 : 0.0);
        const double displacement_fraction =
            highland_forward_displacement_fraction(
                KE, step_gcm2, X0, cumulative_x_X0);
        const double displacement_cm =
            step_cm * displacement_fraction;
        pos += step_direction * displacement_cm;
        // A condensed-history step is limited by true path. With scattering,
        // its forward displacement is shorter than that path, so consuming a
        // path numerically equal to the straight boundary distance does not
        // imply that the particle crossed the boundary. Re-trace from the
        // corrected endpoint instead.
        if (reach_boundary
            && displacement_cm >= step_cm * (1.0 - 1e-12)) {
            pos += step_direction * 1e-6;
        } else {
            reach_boundary = false;
        }

        // --- Bremsstrahlung emission (identical to deposited_in_scoring) ---
        double brems_keV = 0.0;
        if (X0 > 0.0 && KE > kBremsThreshold_keV) {
            brems_keV = sample_step_bremsstrahlung(
                *this, mat, dir, KE, step_gcm2, is_positron, rng,
                [&](const Eigen::Vector3d& photon_direction,
                    double photon_energy_keV) {
                        const Eigen::Vector3d emission_position =
                            step_start
                            + step_direction
                                * (uniform(rng) * displacement_cm);
                        result.brems_photons.push_back(
                            {emission_position, photon_direction,
                             photon_energy_keV, charged_generation});
                });
        }

        // Electron-impact K x-rays over the step (rate x path), born along it.
        if (impact && KE <= impact_emax) {
            std::vector<ImpactXray> xrays;
            sample_impact_xrays(mat, KE, impact_xray_rate_per_gcm2(mat, KE) * step_gcm2,
                                false, rng, xrays);
            for (const ImpactXray& x : xrays)
                push_impact_xray(step_start + step_direction * (uniform(rng) * displacement_cm), x);
        }

        // Explicit energetic Møller/Bhabha delta electrons.  The collision
        // stopping power already removes their energy from the requested
        // primary, but treating that loss as entirely local omits brems from
        // the secondary charged-particle track.  Sampling deltas above 50 keV
        // captures the HPGe-visible part without following the numerous
        // sub-threshold ionization electrons.  One secondary generation is
        // sufficient here; lower generations contribute below the present
        // statistical/HPGe resolution target.
        constexpr double kDeltaTransportFloor_keV = 50.0;
        constexpr std::uint8_t kMaximumDeltaGeneration = 1;
        const double maximum_delta_keV =
            is_positron ? KE : 0.5 * KE;
        if (charged_generation < kMaximumDeltaGeneration
            && maximum_delta_keV > kDeltaTransportFloor_keV) {
            const double delta_mean =
                delta_electron_rate_per_gcm2(
                    mat, KE, kDeltaTransportFloor_keV, is_positron)
                * step_gcm2;
            if (delta_mean > 0.0) {
                const int count =
                    sample_capped_poisson(delta_mean, 32, rng);
                for (int delta_index = 0;
                     delta_index < count; ++delta_index) {
                    const double delta_KE_keV =
                        sample_delta_energy_keV(
                            KE, kDeltaTransportFloor_keV,
                            is_positron, rng);
                    if (!(delta_KE_keV > 0.0))
                        continue;
                    const Eigen::Vector3d delta_position =
                        step_start
                        + step_direction
                            * (uniform(rng) * displacement_cm);
                    const Eigen::Vector3d delta_direction =
                        sample_delta_direction(
                            step_direction, KE, delta_KE_keV, rng);

                    // A knock-on whose full CSDA path is shorter than its
                    // clearance to the nearest material boundary cannot escape
                    // or enter another material. Avoid 20--60 geometry traces
                    // for that common case and use the same bulk SB sampler.
                    // Surface-reachable knock-ons retain the full geometry walk
                    // and its spatially resolved brems origins.
                    const double delta_range_cm =
                        mat.density() > 0.0
                            ? range_gcm2_material(
                                  mat, delta_KE_keV, /*is_positron=*/false)
                                / mat.density()
                            : 0.0;
                    const double material_clearance_cm =
                        source.min_distance_within_material(
                            delta_position, mat);
                    // Require another full range of clearance beyond strict
                    // charged-particle containment. Near-surface deltas cannot
                    // themselves escape, but resolving their brems creation
                    // positions still matters because those photons can.
                    constexpr double kBulkClearanceRanges = 2.0;
                    if (delta_range_cm > 0.0
                        && kBulkClearanceRanges * delta_range_cm
                            < material_clearance_cm) {
                        auto delta_photons = sample_brems_in_material(
                            mat, delta_position, delta_direction,
                            delta_KE_keV, rng, charged_generation + 1);
                        for (auto& photon : delta_photons)
                            result.brems_photons.push_back(std::move(photon));
                        if (impact) {
                            const double ke = std::min(delta_KE_keV, impact_emax);
                            std::vector<ImpactXray> xrays;
                            sample_impact_xrays(mat, ke, impact_xray_yield(mat, ke), true,
                                                rng, xrays);
                            for (const ImpactXray& x : xrays)
                                push_impact_xray(delta_position, x);
                        }
                    } else {
                        ElectronSourceWalkResult delta_walk =
                            walk_in_source_geometry(
                                source, mat, delta_position, delta_direction,
                                delta_KE_keV, rng, false,
                                charged_generation + 1);
                        for (auto& photon : delta_walk.brems_photons)
                            result.brems_photons.push_back(std::move(photon));
                    }
                }
            }
        }

        // CSDA collision loss; conserve energy with the emitted brems photon.
        double KE_exit =
            residual_energy_keV(mat, KE, step_gcm2, is_positron);
        double delta_KE = KE - KE_exit;
        if (delta_KE < 0.0) delta_KE = 0.0;
        // Phase D1: per-step Gaussian energy-loss straggling (Bohr variance
        // σ² = ξ_B·(Z/A)·Δx·(1−β²/2)); lets near-threshold electrons stop/escape
        // against the deterministic CSDA value. Principled mode only.
        if (principled && !dbg_no_straggle && step_gcm2 > 1e-12) {
            double beta2  = KE * (KE + 2.0 * m_e) / ((KE + m_e) * (KE + m_e));
            double ZoverA = 0.0;
            for (const auto& c : mat.composition()) {
                const int Ze = electron_table_z(c.Z);
                ZoverA += c.mass_fraction * Ze / atomic_weight(Ze);
            }
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
            double bcp =
                KE * (KE + 2.0 * m_e) / (KE + m_e);

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
                    double ln_x = std::log(
                        (principled && !dbg_no_b2)
                            ? sum_x_X0 : x_X0);
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

    result.stop_position = pos;
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
