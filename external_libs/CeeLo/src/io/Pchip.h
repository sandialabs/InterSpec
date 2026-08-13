#ifndef CEELO_IO_PCHIP_H
#define CEELO_IO_PCHIP_H
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

/// @file Pchip.h
/// @brief Monotone cubic Hermite interpolation (PCHIP), scipy-compatible.
///
/// Fritsch-Carlson monotone cubic interpolation with the derivative rule of
/// scipy.interpolate.PchipInterpolator (weighted harmonic mean at interior
/// points, shape-preserving one-sided formula at the ends). Chosen over
/// natural splines to forbid overshoot at the low-E efficiency turnover and
/// the pair-production onset, and over Akima for its monotonicity guarantee
/// (response_function_spec.md sec 3.1). Verified against scipy-generated
/// reference values in tests/test_detector_response_io.cpp.
///
/// ATTRIBUTION: the interpolant is Fritsch & Carlson, "Monotone piecewise
/// cubic interpolation", SIAM J. Numer. Anal. 17 (1980) 238-246. The specific
/// derivative rule (interior weighted harmonic mean; the three-point
/// shape-preserving end formula) follows scipy.interpolate.PchipInterpolator,
/// which is distributed under the 3-clause BSD licence -- reimplemented here in
/// C++ rather than translated line by line. See THIRD-PARTY-NOTICES.md.

#include <cstddef>
#include <vector>

namespace ceelo {

/// 1-D PCHIP through (x, y) nodes; x strictly ascending, size >= 2.
/// Evaluation outside [x.front(), x.back()] clamps to the end values
/// (log-log response curves must never be extrapolated -- spec sec 3.1).
class Pchip {
public:
    Pchip() = default;
    Pchip(std::vector<double> x, std::vector<double> y);

    double operator()(double x) const;

    bool valid() const { return x_.size() >= 2; }
    double x_min() const { return x_.empty() ? 0.0 : x_.front(); }
    double x_max() const { return x_.empty() ? 0.0 : x_.back(); }
    const std::vector<double>& x() const { return x_; }
    const std::vector<double>& y() const { return y_; }

private:
    std::vector<double> x_, y_, d_;  // nodes + endpoint-slopes
};

/// Piecewise-linear "hat" basis on a knot grid: B_j(x) is 1 at knot j, 0 at
/// its neighbors, clamped so the end basis functions extend flat beyond the
/// range. Linear in the coefficients -- exactly what the grounding-k(E)
/// covariance propagation Cov[ln k(E), ln k(E')] = B(E) C B(E')^T needs
/// (spec Eqs. 7b-8a; a PCHIP through the knots is NOT linear in the knot
/// values, so the uncertainty basis is the hats).
/// Returns the basis row vector B(x) of length knots.size().
std::vector<double> hat_basis(double x, const std::vector<double>& knots);

} // namespace ceelo

#endif // CEELO_IO_PCHIP_H
