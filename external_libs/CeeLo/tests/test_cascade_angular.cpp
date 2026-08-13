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

// Unit tests for the gamma-gamma angular-correlation sampling helpers
// (src/cascade/AngularSampling.h): W(theta) cos-theta sampler + direction build.
//
// For the pdf f(c) = (1 + a2*P2(c) + a4*P4(c))/2 on c in [-1,1], Legendre
// orthogonality gives E[P2(c)] = a2/5 and E[P4(c)] = a4/9. We check the sample
// means against those, the isotropic limit, and the geometry of direction_at_angle.

#define BOOST_TEST_MODULE cascade_angular
#include <boost/test/included/unit_test.hpp>

#include "cascade/AngularSampling.h"

#include <cmath>
#include <random>

using namespace ceelo;

namespace {
std::mt19937_64 rng(12345);
std::uniform_real_distribution<double> U(0.0, 1.0);
}

BOOST_AUTO_TEST_CASE(w_sampler_moments_co60) {
    // Co-60 coefficients.
    const double a2 = 0.1020, a4 = 0.0091;
    const int N = 2'000'000;
    double sp2 = 0.0, sp4 = 0.0;
    for (int i = 0; i < N; ++i) {
        const double c = sample_cos_theta_W(a2, a4, U, rng);
        BOOST_REQUIRE(c >= -1.0 && c <= 1.0);
        sp2 += legendre_p2(c);
        sp4 += legendre_p4(c);
    }
    const double e2 = sp2 / N, e4 = sp4 / N;
    // Statistical error ~ 1/sqrt(N) ~ 7e-4; use a 4-sigma-ish window.
    BOOST_CHECK_CLOSE(e2, a2 / 5.0, /*pct*/ 6.0);   // a2/5 = 0.0204
    BOOST_CHECK_SMALL(e4 - a4 / 9.0, 3e-3);          // a4/9 = 0.00101
}

BOOST_AUTO_TEST_CASE(w_sampler_isotropic_limit) {
    const int N = 1'000'000;
    double sp2 = 0.0;
    for (int i = 0; i < N; ++i)
        sp2 += legendre_p2(sample_cos_theta_W(0.0, 0.0, U, rng));
    BOOST_CHECK_SMALL(sp2 / N, 3e-3);  // E[P2]=0 for uniform
}

BOOST_AUTO_TEST_CASE(w_sampler_strong_dipole) {
    // 0(E1)1(E1)0 -> a2 = 0.5: E[P2] should be 0.1.
    const int N = 1'000'000;
    double sp2 = 0.0;
    for (int i = 0; i < N; ++i)
        sp2 += legendre_p2(sample_cos_theta_W(0.5, 0.0, U, rng));
    BOOST_CHECK_CLOSE(sp2 / N, 0.5 / 5.0, 5.0);
}

BOOST_AUTO_TEST_CASE(direction_at_angle_geometry) {
    std::mt19937_64 r(999);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    for (int t = 0; t < 200; ++t) {
        // Random reference axis.
        const double rc = 2.0 * u(r) - 1.0, rs = std::sqrt(1.0 - rc * rc),
                     rp = kAngularTwoPi * u(r);
        Eigen::Vector3d ref(rs * std::cos(rp), rs * std::sin(rp), rc);
        const double cos_theta = 2.0 * u(r) - 1.0;
        const double phi = kAngularTwoPi * u(r);
        Eigen::Vector3d d = direction_at_angle(ref, cos_theta, phi);
        BOOST_CHECK_CLOSE(d.norm(), 1.0, 1e-6);
        BOOST_CHECK_SMALL(d.dot(ref) - cos_theta, 1e-9);  // angle to ref is theta
    }
}
