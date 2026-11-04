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

#include "physics/BremsAngular.h"

#include "transport/GaussLegendre.h"

#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>

namespace ceelo {

namespace {

constexpr double kMeC2_keV = 510.9989461;   // Electron rest energy (keV)
constexpr double kTwoPi = 2.0 * 3.14159265358979323846;

/// Everything the 2BS shape needs, computed once per (KE, k, Z).
struct TwoBSParams {
    bool valid = false;
    double ratio = 0.0;    // E_f / E_i (total energies)
    double ratio1 = 0.0;   // (1 + r)^2
    double ratio2 = 0.0;   // 1 + r^2
    double fz = 0.0;       // 8.116224e-5 * (Z (Z + 1))^(2/3)
    double y_max = 0.0;    // 2 beta (1 + beta) gamma^2
    double ceiling = 0.0;  // M = max(f(0), f(y_max))
};

/// The shared 2BS rejection function f(y). Sampler and evaluator both call
/// THIS -- any change here changes both in lockstep.
double rejection(const TwoBSParams& p, double y)
{
    const double one_plus_y_squared = (1.0 + y) * (1.0 + y);
    const double x = 4.0 * y * p.ratio / one_plus_y_squared;
    return 4.0 * x - p.ratio1
        - (p.ratio2 - x) * std::log(p.fz / one_plus_y_squared);
}

TwoBSParams make_params(double electron_kinetic_energy_keV,
                        double photon_energy_keV,
                        double atomic_number)
{
    TwoBSParams p;
    const double total_energy = electron_kinetic_energy_keV + kMeC2_keV;
    if (!(total_energy > kMeC2_keV)
        || !(photon_energy_keV > 0.0)
        || photon_energy_keV >= total_energy) {
        return p;
    }
    p.ratio = (total_energy - photon_energy_keV) / total_energy;
    p.ratio1 = (1.0 + p.ratio) * (1.0 + p.ratio);
    p.ratio2 = 1.0 + p.ratio * p.ratio;
    const double gamma = total_energy / kMeC2_keV;
    const double beta =
        std::sqrt(std::max((gamma - 1.0) * (gamma + 1.0), 0.0)) / gamma;
    const double z = std::clamp(atomic_number, 1.0, 92.0);
    p.fz = 0.00008116224 * std::pow(z * (z + 1.0), 2.0 / 3.0);
    p.y_max = 2.0 * beta * (1.0 + beta) * gamma * gamma;
    if (!(p.y_max > 0.0) || !(p.fz > 0.0))
        return p;
    p.ceiling = std::max(rejection(p, 0.0), rejection(p, p.y_max));
    if (!(p.ceiling > 0.0) || !std::isfinite(p.ceiling))
        return p;
    p.valid = true;
    return p;
}

/// min(max(f(y), 0), M): exactly what the accept step realizes.
double clamped_rejection(const TwoBSParams& p, double y)
{
    return std::min(std::max(rejection(p, y), 0.0), p.ceiling);
}

/// Build a perpendicular basis for direction d (normalized): (u, v, d)
/// orthonormal. Same construction as ElectronCsda.cpp's file-local helper.
void perp_basis(const Eigen::Vector3d& d, Eigen::Vector3d& u, Eigen::Vector3d& v)
{
    if (std::abs(d.x()) <= std::abs(d.y()) && std::abs(d.x()) <= std::abs(d.z())) {
        u = Eigen::Vector3d(0.0, -d.z(), d.y()).normalized();
    } else if (std::abs(d.y()) <= std::abs(d.z())) {
        u = Eigen::Vector3d(-d.z(), 0.0, d.x()).normalized();
    } else {
        u = Eigen::Vector3d(-d.y(), d.x(), 0.0).normalized();
    }
    v = d.cross(u);
}

} // namespace

double brems_2bs_ymax(double electron_kinetic_energy_keV)
{
    if (!(electron_kinetic_energy_keV > 0.0))
        return 0.0;
    const double gamma = 1.0 + electron_kinetic_energy_keV / kMeC2_keV;
    const double beta =
        std::sqrt(std::max((gamma - 1.0) * (gamma + 1.0), 0.0)) / gamma;
    return 2.0 * beta * (1.0 + beta) * gamma * gamma;
}

double brems_2bs_pdf_y_unnorm(double electron_kinetic_energy_keV,
                              double photon_energy_keV,
                              double atomic_number,
                              double y)
{
    const TwoBSParams p = make_params(
        electron_kinetic_energy_keV, photon_energy_keV, atomic_number);
    if (!p.valid || y < 0.0 || y > p.y_max)
        return 0.0;
    return clamped_rejection(p, y) / ((1.0 + y) * (1.0 + y));
}

double brems_2bs_norm(double electron_kinetic_energy_keV,
                      double photon_energy_keV,
                      double atomic_number)
{
    const TwoBSParams p = make_params(
        electron_kinetic_energy_keV, photon_energy_keV, atomic_number);
    if (!p.valid)
        return 0.0;
    // Exact cancellation: with y(q) = q y_max / (1 + y_max (1 - q)),
    // dy/dq = y_max (1 + y_max) / D^2 and (1 + y) = (1 + y_max)/D, so
    // [dy/dq] / (1+y)^2 = y_max / (1 + y_max) -- a constant. The q-space
    // integrand is just the clamped rejection function: bounded and smooth
    // except for the single kink where f crosses zero, which the composite
    // panels absorb.
    const double q_integral = gl::integrate32_composite(0.0, 1.0, 4,
        [&](double q) {
            const double y = q * p.y_max / (1.0 + p.y_max * (1.0 - q));
            return clamped_rejection(p, y);
        });
    return (p.y_max / (1.0 + p.y_max)) * q_integral;
}

double brems_2bs_pdf_sr(double electron_kinetic_energy_keV,
                        double photon_energy_keV,
                        double atomic_number,
                        double mu)
{
    const TwoBSParams p = make_params(
        electron_kinetic_energy_keV, photon_energy_keV, atomic_number);
    if (!p.valid || !(mu >= -1.0 && mu <= 1.0))
        return 0.0;
    const double y = p.y_max * (1.0 - mu) * 0.5;
    const double unnorm =
        clamped_rejection(p, y) / ((1.0 + y) * (1.0 + y));
    const double norm = brems_2bs_norm(
        electron_kinetic_energy_keV, photon_energy_keV, atomic_number);
    if (!(norm > 0.0))
        return 0.0;
    return unnorm * (p.y_max * 0.5) / norm / kTwoPi;
}

Eigen::Vector3d sample_brems_direction_2bs(
    const Eigen::Vector3d& electron_direction,
    double electron_kinetic_energy_keV,
    double photon_energy_keV,
    double atomic_number,
    std::mt19937_64& rng)
{
    const TwoBSParams p = make_params(
        electron_kinetic_energy_keV, photon_energy_keV, atomic_number);
    if (!p.valid)
        return electron_direction;

    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double y = 0.0;
    bool accepted = false;
    for (int attempt = 0; attempt < 10000; ++attempt) {
        const double q = uniform(rng);
        y = q * p.y_max / (1.0 + p.y_max * (1.0 - q));
        const double value = rejection(p, y);
        // Short-circuit preserved: no accept draw when f < 0 (keeps the RNG
        // stream identical to the pre-hoist implementation).
        if (value >= 0.0 && uniform(rng) * p.ceiling <= value) {
            accepted = true;
            break;
        }
    }
    if (!accepted)
        return electron_direction;

    const double cosine = std::clamp(1.0 - 2.0 * y / p.y_max, -1.0, 1.0);
    const double sine = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
    const double phi = kTwoPi * uniform(rng);
    Eigen::Vector3d u, v;
    perp_basis(electron_direction, u, v);
    return (
        cosine * electron_direction
        + sine * (std::cos(phi) * u + std::sin(phi) * v)
    ).normalized();
}

} // namespace ceelo
