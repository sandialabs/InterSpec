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

#include "cascade/BetaSpectrum.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace ceelo {
namespace {

constexpr double kElectronMassKeV = 510.998950;
constexpr double kFineStructure = 7.2973525693e-3;
constexpr double kElectronComptonFm = 386.15926764;
constexpr double kPi = 3.14159265358979323846;

std::complex<double> log_gamma(std::complex<double> z)
{
    // Lanczos approximation of log(Gamma(z)).
    //
    // Reference: C. Lanczos, "A precision approximation of the gamma
    // function", J. SIAM Numer. Anal. Ser. B 1 (1964) 86-96.
    //
    // The g=7, n=9 coefficient set below is Godfrey's, published openly and
    // reproduced in Boost.Math (lanczos.hpp) and Press et al. among others; it
    // is a numerical constant of the published approximation, not code taken
    // from any implementation.  (Numerical Recipes, which also tabulates it, is
    // NOT redistributable -- nothing here derives from it.)
    static constexpr double coefficients[] = {
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    };
    if (z.real() < 0.5) {
        return std::log(kPi) - std::log(std::sin(kPi * z))
            - log_gamma(1.0 - z);
    }
    z -= 1.0;
    std::complex<double> x = coefficients[0];
    for (std::size_t i = 1; i < std::size(coefficients); ++i)
        x += coefficients[i] / (z + static_cast<double>(i));
    const std::complex<double> t = z + 7.5;
    return 0.5 * std::log(2.0 * kPi)
        + (z + 0.5) * std::log(t) - t + std::log(x);
}

double fermi_factor(const BetaBranch& branch, double total_energy_electron_mass,
                    double momentum_electron_mass)
{
    if (!(momentum_electron_mass > 0.0))
        return 0.0;
    const double az = kFineStructure
        * static_cast<double>(std::max(branch.daughter_Z, 0));
    if (!(az < 1.0))
        return 0.0;
    const double gamma = std::sqrt(1.0 - az * az);
    const double charge_sign = branch.is_positron ? -1.0 : 1.0;
    const double y = charge_sign * az * total_energy_electron_mass
        / momentum_electron_mass;
    const double radius =
        1.2 * std::cbrt(static_cast<double>(std::max(branch.daughter_A, 1)))
        / kElectronComptonFm;
    const std::complex<double> argument(gamma, y);
    double log_f = std::log(2.0 * (1.0 + gamma))
        + 2.0 * (gamma - 1.0)
            * std::log(std::max(2.0 * momentum_electron_mass * radius,
                                std::numeric_limits<double>::min()))
        + kPi * y
        + 2.0 * log_gamma(argument).real()
        - 2.0 * std::lgamma(2.0 * gamma + 1.0);
    log_f = std::clamp(log_f, -700.0, 700.0);
    return std::exp(log_f);
}

double forbidden_shape_factor(BetaShape shape, double p, double q)
{
    const double p2 = p * p;
    const double q2 = q * q;
    switch (shape) {
    case BetaShape::FirstUniqueForbidden:
        return p2 + q2;
    case BetaShape::SecondUniqueForbidden:
        return p2 * p2 + (10.0 / 3.0) * p2 * q2 + q2 * q2;
    case BetaShape::ThirdUniqueForbidden:
        return p2 * p2 * p2 + 7.0 * p2 * p2 * q2
            + 7.0 * p2 * q2 * q2 + q2 * q2 * q2;
    case BetaShape::Allowed:
    case BetaShape::FirstForbidden:
    case BetaShape::SecondForbidden:
    case BetaShape::ThirdForbidden:
    case BetaShape::FourthForbidden:
        return 1.0;
    }
    return 1.0;
}

double sum_beta_yield(const std::vector<BetaBranch>& branches)
{
    return std::accumulate(
        branches.begin(), branches.end(), 0.0,
        [](double sum, const BetaBranch& branch) {
            return sum + std::max(branch.yield_per_decay, 0.0);
        });
}

} // namespace

double RadioactiveEmissionSet::beta_yield() const
{
    return sum_beta_yield(beta_branches);
}

double RadioactiveEmissionSet::conversion_electron_yield() const
{
    return std::accumulate(
        conversion_electrons.begin(), conversion_electrons.end(), 0.0,
        [](double sum, const ConversionElectronLine& line) {
            return sum + std::max(line.yield_per_decay, 0.0);
        });
}

double RadioactiveEmissionSet::photon_yield() const
{
    return std::accumulate(
        photons.begin(), photons.end(), 0.0,
        [](double sum, const RadioactivePhotonLine& line) {
            return sum + std::max(line.yield_per_decay, 0.0);
        });
}

BetaSpectrumSampler::BetaSpectrumSampler(const BetaBranch& branch,
                                         std::size_t grid_points)
    : branch_(branch)
{
    if (!(branch.endpoint_keV > 0.0))
        throw std::invalid_argument("beta endpoint must be positive");
    grid_points = std::max<std::size_t>(grid_points, 65);
    energy_grid_keV_.resize(grid_points);
    cdf_.assign(grid_points, 0.0);

    double previous_density = 0.0;
    for (std::size_t i = 0; i < grid_points; ++i) {
        const double fraction =
            static_cast<double>(i) / static_cast<double>(grid_points - 1);
        const double energy = fraction * branch.endpoint_keV;
        energy_grid_keV_[i] = energy;
        const double density = unnormalized_density(branch, energy);
        if (i > 0) {
            const double width = energy - energy_grid_keV_[i - 1];
            cdf_[i] = cdf_[i - 1]
                + 0.5 * width * (previous_density + density);
        }
        previous_density = density;
    }
    normalization_ = cdf_.back();
    if (!(normalization_ > 0.0) || !std::isfinite(normalization_))
        throw std::runtime_error("beta spectrum normalization failed");
    for (double& value : cdf_)
        value /= normalization_;
    cdf_.back() = 1.0;
}

double BetaSpectrumSampler::unnormalized_density(
    const BetaBranch& branch, double kinetic_energy_keV)
{
    if (!(kinetic_energy_keV > 0.0)
        || !(kinetic_energy_keV < branch.endpoint_keV))
        return 0.0;
    const double W = 1.0 + kinetic_energy_keV / kElectronMassKeV;
    const double p = std::sqrt(std::max(W * W - 1.0, 0.0));
    const double q =
        (branch.endpoint_keV - kinetic_energy_keV) / kElectronMassKeV;
    const double fermi = fermi_factor(branch, W, p);
    const double shape = forbidden_shape_factor(branch.shape, p, q);
    const double density = p * W * q * q * fermi * shape;
    return std::isfinite(density) ? std::max(density, 0.0) : 0.0;
}

double BetaSpectrumSampler::normalized_density_per_keV(
    double kinetic_energy_keV) const
{
    return unnormalized_density(branch_, kinetic_energy_keV) / normalization_;
}

double BetaSpectrumSampler::sample_keV(std::mt19937_64& rng) const
{
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const double value = uniform(rng);
    const auto upper = std::lower_bound(cdf_.begin(), cdf_.end(), value);
    if (upper == cdf_.begin())
        return energy_grid_keV_.front();
    if (upper == cdf_.end())
        return energy_grid_keV_.back();
    const std::size_t hi =
        static_cast<std::size_t>(std::distance(cdf_.begin(), upper));
    const std::size_t lo = hi - 1;
    const double span = cdf_[hi] - cdf_[lo];
    const double fraction =
        span > 0.0 ? (value - cdf_[lo]) / span : 0.0;
    return energy_grid_keV_[lo]
        + fraction * (energy_grid_keV_[hi] - energy_grid_keV_[lo]);
}

} // namespace ceelo
