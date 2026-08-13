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

// Angular-correlation sampling helpers for the correlated-cascade MC.
//
// W(theta) = 1 + a2*P2(cos theta) + a4*P4(cos theta) is the gamma-gamma
// directional correlation; the coefficients a2/a4 are supplied per coincidence
// link (CascadeTypes.h, populated by the enriched SandiaDecay data). These
// helpers sample a partner direction relative to a reference direction. The pdf
// integrates to 1 over the sphere (P2,P4 average to zero), so sampling from it
// introduces no estimator weight.

#include <algorithm>
#include <cmath>
#include <random>

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace ceelo {

constexpr double kAngularTwoPi = 6.283185307179586;

inline double legendre_p2(double c) { return 0.5 * (3.0 * c * c - 1.0); }
inline double legendre_p4(double c) {
    const double c2 = c * c;
    return 0.125 * (35.0 * c2 * c2 - 30.0 * c2 + 3.0);
}

/// Sample cos(theta) in [-1,1] from f(c) = (1 + a2*P2(c) + a4*P4(c)) / 2 by
/// rejection against a constant envelope (max of the bounded quartic over a fine
/// scan, clamped to >= 0). Returns a uniform draw for negligible coefficients.
inline double sample_cos_theta_W(double a2, double a4,
                                 std::uniform_real_distribution<double>& u,
                                 std::mt19937_64& rng) {
    if (std::abs(a2) < 1e-12 && std::abs(a4) < 1e-12)
        return 2.0 * u(rng) - 1.0;
    auto W = [a2, a4](double c) {
        return std::max(0.0, 1.0 + a2 * legendre_p2(c) + a4 * legendre_p4(c));
    };
    double m = 0.0;
    const int kScan = 64;
    for (int i = 0; i <= kScan; ++i)
        m = std::max(m, W(-1.0 + 2.0 * i / kScan));
    m *= 1.0001;
    for (int tries = 0; tries < 1000; ++tries) {
        const double c = 2.0 * u(rng) - 1.0;
        if (u(rng) * m <= W(c))
            return c;
    }
    return 2.0 * u(rng) - 1.0;
}

/// Unit direction at polar angle theta (cos_theta given) about `ref`, uniform
/// azimuth phi.
inline Eigen::Vector3d direction_at_angle(const Eigen::Vector3d& ref,
                                          double cos_theta, double phi) {
    const double st = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    // Pick a seed axis not parallel to ref to build an orthonormal basis.
    Eigen::Vector3d a = (std::abs(ref.x()) <= std::abs(ref.y()) &&
                         std::abs(ref.x()) <= std::abs(ref.z()))
                            ? Eigen::Vector3d::UnitX()
                            : (std::abs(ref.y()) <= std::abs(ref.z())
                                   ? Eigen::Vector3d::UnitY()
                                   : Eigen::Vector3d::UnitZ());
    Eigen::Vector3d e1 = a - a.dot(ref) * ref;
    if (e1.norm() < 1e-12) {
        a = Eigen::Vector3d::UnitZ();
        e1 = a - a.dot(ref) * ref;
    }
    e1.normalize();
    const Eigen::Vector3d e2 = ref.cross(e1);
    return cos_theta * ref + st * (std::cos(phi) * e1 + std::sin(phi) * e2);
}

} // namespace ceelo
