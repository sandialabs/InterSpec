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

/// @file Pchip.cpp
/// @brief scipy-compatible PCHIP (see Pchip.h).

#include "io/Pchip.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>

namespace ceelo {

namespace {

inline double sign_of(double v) { return (v > 0.0) - (v < 0.0); }

// scipy PchipInterpolator._edge_case: shape-preserving three-point estimate.
double edge_case(double h0, double h1, double m0, double m1) {
    double d = ((2.0 * h0 + h1) * m0 - h0 * m1) / (h0 + h1);
    if (sign_of(d) != sign_of(m0)) {
        d = 0.0;
    } else if (sign_of(m0) != sign_of(m1) && std::abs(d) > 3.0 * std::abs(m0)) {
        d = 3.0 * m0;
    }
    return d;
}

} // namespace

Pchip::Pchip(std::vector<double> x, std::vector<double> y)
    : x_(std::move(x)), y_(std::move(y)) {
    const size_t n = x_.size();
    if (n < 2 || y_.size() != n)
        throw std::invalid_argument("Pchip: need >= 2 nodes, matching sizes");
    for (size_t i = 1; i < n; ++i)
        if (!(x_[i] > x_[i - 1]))
            throw std::invalid_argument("Pchip: x must be strictly ascending");

    d_.resize(n);
    if (n == 2) {
        const double m = (y_[1] - y_[0]) / (x_[1] - x_[0]);
        d_[0] = d_[1] = m;
        return;
    }

    std::vector<double> h(n - 1), m(n - 1);
    for (size_t i = 0; i + 1 < n; ++i) {
        h[i] = x_[i + 1] - x_[i];
        m[i] = (y_[i + 1] - y_[i]) / h[i];
    }
    // Interior: weighted harmonic mean where slopes agree in sign, else 0
    // (scipy PchipInterpolator._find_derivatives).
    for (size_t i = 1; i + 1 < n; ++i) {
        const double mk0 = m[i - 1], mk1 = m[i];
        if (sign_of(mk0) != sign_of(mk1) || mk0 == 0.0 || mk1 == 0.0) {
            d_[i] = 0.0;
        } else {
            const double w1 = 2.0 * h[i] + h[i - 1];
            const double w2 = h[i] + 2.0 * h[i - 1];
            d_[i] = (w1 + w2) / (w1 / mk0 + w2 / mk1);
        }
    }
    d_[0] = edge_case(h[0], h[1], m[0], m[1]);
    d_[n - 1] = edge_case(h[n - 2], h[n - 3], m[n - 2], m[n - 3]);
}

double Pchip::operator()(double xq) const {
    assert(valid());
    if (xq <= x_.front()) return y_.front();
    if (xq >= x_.back()) return y_.back();
    // Find interval i with x_[i] <= xq < x_[i+1].
    const size_t i = static_cast<size_t>(
        std::upper_bound(x_.begin(), x_.end(), xq) - x_.begin()) - 1;
    const double h = x_[i + 1] - x_[i];
    const double t = (xq - x_[i]) / h;
    const double t2 = t * t, t3 = t2 * t;
    const double h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    const double h10 = t3 - 2.0 * t2 + t;
    const double h01 = -2.0 * t3 + 3.0 * t2;
    const double h11 = t3 - t2;
    return h00 * y_[i] + h10 * h * d_[i] + h01 * y_[i + 1] + h11 * h * d_[i + 1];
}

std::vector<double> hat_basis(double x, const std::vector<double>& knots) {
    const size_t n = knots.size();
    std::vector<double> B(n, 0.0);
    if (n == 0) return B;
    if (n == 1) { B[0] = 1.0; return B; }
    if (x <= knots.front()) { B[0] = 1.0; return B; }
    if (x >= knots.back()) { B[n - 1] = 1.0; return B; }
    const size_t i = static_cast<size_t>(
        std::upper_bound(knots.begin(), knots.end(), x) - knots.begin()) - 1;
    const double t = (x - knots[i]) / (knots[i + 1] - knots[i]);
    B[i] = 1.0 - t;
    B[i + 1] = t;
    return B;
}

} // namespace ceelo
