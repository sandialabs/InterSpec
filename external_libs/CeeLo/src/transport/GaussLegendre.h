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

/// @file GaussLegendre.h
/// @brief Shared 32-point Gauss-Legendre nodes/weights and a composite
///        integrator over the substituted variable s = sin(theta/2).
///
/// The s-substitution (mu = 1 - 2 s^2, dmu = -4 s ds) is used by every angular
/// normalization in the transport layer. It makes the momentum transfer
/// x = (E/hc) * s LINEAR in the integration variable, removing the sqrt(1-mu)
/// cusp that S(x,Z) and F(x,Z) imprint on the mu-space integrand at mu -> 1
/// (which otherwise limits any quadrature to ~1e-4 at MeV energies).
///
/// Originally private to ComptonScatter.cpp; hoisted here when the Rayleigh
/// evaluator and the next-event estimator's cone integrals needed the same
/// machinery.

#include <cmath>

namespace ceelo {
namespace gl {

/// 32-point Gauss-Legendre on [-1, 1]: the positive half of the symmetric
/// node set (16 nodes, each used twice as +/-x).
constexpr int kGL32Half = 16;

constexpr double kGL32_x[kGL32Half] = {
    0.0483076656877383162, 0.1444719615827964934, 0.2392873622521370745,
    0.3318686022821276497, 0.4213512761306353454, 0.5068999089322293900,
    0.5877157572407623290, 0.6630442669302152010, 0.7321821187402896804,
    0.7944837959679424070, 0.8493676137325699701, 0.8963211557660521240,
    0.9349060759377396892, 0.9647622555875064308, 0.9856115115452683354,
    0.9972638618494815635};

constexpr double kGL32_w[kGL32Half] = {
    0.0965400885147278006, 0.0956387200792748594, 0.0938443990808045654,
    0.0911738786957638847, 0.0876520930044038111, 0.0833119242269467552,
    0.0781938957870703065, 0.0723457941088485062, 0.0658222227763618468,
    0.0586840934785355471, 0.0509980592623761762, 0.0428358980222266807,
    0.0342738629130214331, 0.0253920653092620595, 0.0162743947309056706,
    0.0070186100094700966};

/// 32-point Gauss-Legendre integral of `f` over [a, b].
template <typename F>
double integrate32(double a, double b, F&& f)
{
    const double c = 0.5 * (a + b), h = 0.5 * (b - a);
    double sum = 0.0;
    for (int i = 0; i < kGL32Half; ++i)
        sum += kGL32_w[i] * (f(c - h * kGL32_x[i]) + f(c + h * kGL32_x[i]));
    return sum * h;
}

/// Composite 32-point Gauss-Legendre over [a, b] in `panels` uniform panels.
template <typename F>
double integrate32_composite(double a, double b, int panels, F&& f)
{
    const double h = (b - a) / panels;
    double sum = 0.0;
    for (int k = 0; k < panels; ++k)
        sum += integrate32(a + k * h, a + (k + 1) * h, f);
    return sum;
}

/// Integrate an angular density g(mu) over mu in [mu_lo, mu_hi] using the
/// s = sin(theta/2) substitution, where s = sqrt((1-mu)/2):
///     integral g(mu) dmu = integral g(1 - 2 s^2) * 4 s ds.
/// `panels` composite panels are used across the s range.
template <typename G>
double integrate_mu_in_s(double mu_lo, double mu_hi, int panels, G&& g)
{
    // mu is decreasing in s, so the s range runs from s(mu_hi) to s(mu_lo).
    auto s_of_mu = [](double mu) {
        if (mu > 1.0) mu = 1.0;
        if (mu < -1.0) mu = -1.0;
        return std::sqrt(0.5 * (1.0 - mu));
    };
    const double s_a = s_of_mu(mu_hi);
    const double s_b = s_of_mu(mu_lo);
    if (s_b <= s_a) return 0.0;
    return integrate32_composite(s_a, s_b, panels, [&](double s) {
        return 4.0 * s * g(1.0 - 2.0 * s * s);
    });
}

} // namespace gl
} // namespace ceelo
