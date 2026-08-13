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

/// @file BetaSpectrum.h
/// @brief Data-agnostic radioactive electron emissions and beta sampling.

#include <cstdint>
#include <random>
#include <string>
#include <vector>

namespace ceelo {

enum class BetaShape : std::uint8_t {
    Allowed,
    FirstForbidden,
    FirstUniqueForbidden,
    SecondForbidden,
    SecondUniqueForbidden,
    ThirdForbidden,
    ThirdUniqueForbidden,
    FourthForbidden
};

struct BetaBranch {
    double endpoint_keV = 0.0;
    double yield_per_decay = 0.0;
    int daughter_Z = 0;
    int daughter_A = 0;
    bool is_positron = false;
    BetaShape shape = BetaShape::Allowed;
};

enum class ConversionShell : std::uint8_t {
    K,
    L1,
    L2,
    L3,
    Outer
};

struct ConversionElectronLine {
    double energy_keV = 0.0;
    double yield_per_decay = 0.0;
    ConversionShell shell = ConversionShell::Outer;
};

struct RadioactivePhotonLine {
    double energy_keV = 0.0;
    double yield_per_decay = 0.0;
};

/// Direct emissions from one parent decay. This deliberately does not include
/// daughter in-growth or cascade coincidence state: it is the source term for
/// source-escape electron/brems validation.
struct RadioactiveEmissionSet {
    std::string nuclide;
    std::vector<BetaBranch> beta_branches;
    std::vector<ConversionElectronLine> conversion_electrons;
    std::vector<RadioactivePhotonLine> photons;

    double beta_yield() const;
    double conversion_electron_yield() const;
    double photon_yield() const;
};

/// Tabulated inverse-CDF sampler for an individual beta branch.
///
/// The shape is the relativistic allowed beta spectrum with finite-nuclear-size
/// Fermi factor. Unique-forbidden branches receive their leading-order shape
/// factor; non-unique forbidden branches fall back to allowed because their
/// nuclear matrix elements are not present in SandiaDecay.
class BetaSpectrumSampler {
public:
    explicit BetaSpectrumSampler(const BetaBranch& branch,
                                 std::size_t grid_points = 2049);

    double sample_keV(std::mt19937_64& rng) const;
    double normalized_density_per_keV(double kinetic_energy_keV) const;

    const BetaBranch& branch() const { return branch_; }
    bool valid() const { return normalization_ > 0.0; }

private:
    static double unnormalized_density(const BetaBranch& branch,
                                       double kinetic_energy_keV);

    BetaBranch branch_;
    std::vector<double> energy_grid_keV_;
    std::vector<double> cdf_;
    double normalization_ = 0.0;
};

} // namespace ceelo
