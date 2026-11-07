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

/// @file FullSpectrumModel.h
/// @brief Jet-compatible source escape, detector response, and composition.
///
/// The implementation intentionally has no Ceres dependency. All operations
/// involving fitted values are templated on Scalar and use unqualified math
/// calls so ceres::Jet overloads participate through ADL.

#include "full_spectrum/SpectralMeasure.h"
#include "geometry/GeometryDescriptor.h"
#include "materials/Material.h"
#include "transport/ComptonScatter.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <mutex>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace ceelo::full_spectrum {

namespace detail {

template <typename Scalar>
inline Scalar scalar_exp(const Scalar& x) {
    using std::exp;
    return exp(x);
}

template <typename Scalar>
inline Scalar scalar_sqrt(const Scalar& x) {
    using std::sqrt;
    return sqrt(x);
}

template <typename Scalar>
inline Scalar scalar_log(const Scalar& x) {
    using std::log;
    return log(x);
}

/// Clamp a cosine into [-1, 1]. Callers routinely obtain cos(theta) from a
/// dot product, which can land a few ulps outside the domain; unclamped that
/// produced sqrt(negative)/acos(>1) and NaN'd an entire spectrum.
template <typename Scalar>
inline Scalar clamped_cosine(const Scalar& x) {
    return x < Scalar(-1) ? Scalar(-1) : (x > Scalar(1) ? Scalar(1) : x);
}

/// Branch-on-value minimum. The two arguments meet transversally wherever it
/// is used (chord selection), so the value is continuous and the derivative
/// jump is the genuine geometric one.
template <typename Scalar>
inline Scalar scalar_min(const Scalar& a, const Scalar& b) {
    return a < b ? a : b;
}

/// Gauss-Legendre nodes and weights on [-1, 1] by Newton iteration on the
/// Legendre polynomial. Deterministic, so the fitter and the evaluator share
/// one quadrature.
inline void gauss_legendre_nodes(
    std::size_t count, double* nodes, double* weights)
{
    constexpr double pi = 3.14159265358979323846;
    const std::size_t half = (count + 1) / 2;
    for (std::size_t i = 0; i < half; ++i) {
        double x = std::cos(pi * (static_cast<double>(i) + 0.75)
                            / (static_cast<double>(count) + 0.5));
        double derivative = 0.0;
        for (int iteration = 0; iteration < 100; ++iteration) {
            double p0 = 1.0, p1 = 0.0;
            for (std::size_t j = 0; j < count; ++j) {
                const double p2 = p1;
                p1 = p0;
                p0 = ((2.0 * static_cast<double>(j) + 1.0) * x * p1
                      - static_cast<double>(j) * p2)
                    / (static_cast<double>(j) + 1.0);
            }
            derivative = static_cast<double>(count) * (x * p0 - p1)
                / (x * x - 1.0);
            const double delta = p0 / derivative;
            x -= delta;
            if (std::abs(delta) < 1.0e-15)
                break;
        }
        nodes[i] = -x;
        nodes[count - 1 - i] = x;
        const double weight = 2.0 / ((1.0 - x * x)
                                     * derivative * derivative);
        weights[i] = weight;
        weights[count - 1 - i] = weight;
    }
}

/// Cached azimuth quadrature: the (node, weight) pairs already folded with
/// the phi = (pi/2) * node mapping and the cos(phi) projection weight, so a
/// caller only multiplies. Keyed by point count; built once per count.
struct AzimuthQuadrature {
    std::vector<double> cos_phi;
    std::vector<double> weight;
};

inline const AzimuthQuadrature& azimuth_quadrature_rule(std::size_t count)
{
    static std::map<std::size_t, AzimuthQuadrature> cache;
    static std::mutex guard;
    std::lock_guard<std::mutex> lock(guard);
    auto found = cache.find(count);
    if (found != cache.end())
        return found->second;
    constexpr double pi = 3.14159265358979323846;
    std::vector<double> nodes(count), weights(count);
    gauss_legendre_nodes(count, nodes.data(), weights.data());
    AzimuthQuadrature rule;
    rule.cos_phi.resize(count);
    rule.weight.resize(count);
    for (std::size_t index = 0; index < count; ++index) {
        const double phi = 0.5 * pi * nodes[index];
        rule.cos_phi[index] = std::cos(phi);
        rule.weight[index] =
            0.5 * pi * weights[index] * std::cos(phi) / 2.0;
    }
    return cache.emplace(count, std::move(rule)).first->second;
}

template <typename Scalar>
inline Scalar scalar_acos(const Scalar& x) {
    using std::acos;
    return acos(x);
}

template <typename Scalar>
inline Scalar scalar_erf(const Scalar& x) {
    using std::erf;
    return erf(x);
}

template <typename Scalar>
inline Scalar scalar_pow(const Scalar& x, double exponent) {
    using std::pow;
    return pow(x, exponent);
}

inline double binomial(int n, int k) {
    if (k < 0 || k > n)
        return 0.0;
    double result = 1.0;
    for (int i = 1; i <= k; ++i)
        result *= static_cast<double>(n - k + i) / static_cast<double>(i);
    return result;
}

/// Global Bernstein values are branch-free in the fitted coordinate, unlike
/// a knot-selected B-spline. This makes them suitable for a Jet-valued
/// incidence angle in the compact detector-rate surface.
template <typename Scalar, std::size_t N>
inline std::array<Scalar, N> bernstein_values(const Scalar& coordinate) {
    static_assert(N > 1, "Bernstein basis requires at least two terms");
    std::array<Scalar, N> answer{};
    for (std::size_t index = 0; index < N; ++index) {
        Scalar value = Scalar(binomial(
            static_cast<int>(N) - 1, static_cast<int>(index)));
        for (std::size_t power = 0; power < index; ++power)
            value *= coordinate;
        for (std::size_t power = index; power + 1 < N; ++power)
            value *= Scalar(1) - coordinate;
        answer[index] = value;
    }
    return answer;
}

/// Integral of the normalized Bernstein density on [lo, hi].
inline double bernstein_integral(int degree, int index, double lo, double hi) {
    if (!(hi > lo))
        return 0.0;
    lo = std::max(0.0, std::min(1.0, lo));
    hi = std::max(0.0, std::min(1.0, hi));
    if (!(hi > lo))
        return 0.0;

    // Four-point Gauss-Legendre integrates a degree-seven Bernstein basis
    // exactly, which is the highest Bernstein degree any payload uses.
    static constexpr double x[4] = {
        -0.8611363115940525752, -0.3399810435848562648,
         0.3399810435848562648,  0.8611363115940525752};
    static constexpr double w[4] = {
        0.3478548451374538574, 0.6521451548625461426,
        0.6521451548625461426, 0.3478548451374538574};
    const double mid = 0.5 * (lo + hi);
    const double half = 0.5 * (hi - lo);
    const double choose = binomial(degree, index);
    double answer = 0.0;
    for (int q = 0; q < 4; ++q) {
        const double t = mid + half * x[q];
        const double basis = (degree + 1.0) * choose
            * std::pow(t, index) * std::pow(1.0 - t, degree - index);
        answer += w[q] * basis;
    }
    return half * answer;
}

// Whether a detector atom sits at a constant offset from the incident energy
// (full energy, pair full energy, escape peaks, K escape), so that an incident
// energy BAND maps onto a band of the same width. The zero-deposit atom, the
// fixed 511 keV line and the absorbed K lines (a dead-layer K x-ray, 9.833 /
// 10.943 keV whatever the incident energy) do not.
inline bool detector_atom_tracks_incident(SpectralComponent component)
{
    switch (component) {
    case SpectralComponent::DetectorFullEnergy:
    case SpectralComponent::DetectorPairFullEnergy:
    case SpectralComponent::DetectorSingleEscape:
    case SpectralComponent::DetectorDoubleEscape:
    case SpectralComponent::DetectorKAlphaEscape:
    case SpectralComponent::DetectorKBetaEscape:
        return true;
    default:
        return false;
    }
}

inline double compton_scattered_energy_keV(
    double incident_energy_keV, double cosine)
{
    constexpr double electron_mass_keV = 510.99895;
    return incident_energy_keV
        / (1.0 + incident_energy_keV / electron_mass_keV
                     * (1.0 - cosine));
}

inline double inverse_energy_loss_warp(
    double energy_keV, double incident_energy_keV,
    double scale = 1.0)
{
    if (!(energy_keV > 0.0))
        return 1.0;
    constexpr double electron_mass_keV = 510.99895;
    const double loss = 0.5 * electron_mass_keV
        * (1.0 / energy_keV - 1.0 / incident_energy_keV);
    return std::max(0.0, std::min(1.0, loss / (scale + loss)));
}

template <typename Scalar, std::size_t N>
inline void append_bernstein_density(
    SpectralMeasure<Scalar>& result,
    const std::vector<double>& edges,
    double support_keV,
    const Scalar& total_yield,
    const std::array<double, N>& positive_coefficients,
    SpectralComponent component)
{
    if (!(support_keV > 0.0))
        return;
    double coefficient_sum = 0.0;
    for (double coefficient : positive_coefficients)
        coefficient_sum += coefficient;
    const int degree = static_cast<int>(N) - 1;

    for (std::size_t b = 0; b + 1 < edges.size(); ++b) {
        const double lo = std::max(0.0, edges[b] / support_keV);
        const double hi = std::min(1.0, edges[b + 1] / support_keV);
        if (!(hi > lo))
            continue;
        double fraction = 0.0;
        for (std::size_t k = 0; k < N; ++k) {
            fraction += positive_coefficients[k] / coefficient_sum
                * bernstein_integral(degree, static_cast<int>(k), lo, hi);
        }
        result.continuous.push_back(
            {edges[b], edges[b + 1], total_yield * fraction, component});
    }
}

/// Largest coefficient count in a single block of a block-partitioned detector
/// continuum, and the largest number of blocks any component uses.
constexpr std::size_t kDetectorMaxBlockSize = 96;
constexpr std::size_t kDetectorMaxBlocks = 4;
/// Largest Bernstein basis evaluated at run time (the detected-rate surface).
constexpr std::size_t kMaxRuntimeBernsteinSize = 32;

/// Bernstein basis of run-time size, evaluated by the branch-free degree
/// recurrence b_k^(n) = (1-x) b_k^(n-1) + x b_(k-1)^(n-1).
///
/// The detected-rate surface's angle coordinate depends on cos(theta), which a
/// caller may differentiate, so this must stay free of comparisons on the
/// coordinate; the recurrence is pure arithmetic and therefore ceres::Jet safe.
template <typename Scalar>
inline void bernstein_values_runtime(
    std::size_t size, const Scalar& coordinate, Scalar* answer)
{
    const Scalar complement = Scalar(1) - coordinate;
    answer[0] = Scalar(1);
    for (std::size_t degree = 1; degree < size; ++degree) {
        Scalar previous = Scalar(0);
        for (std::size_t index = 0; index < degree; ++index) {
            const Scalar current = answer[index];
            answer[index] = complement * current + coordinate * previous;
            previous = current;
        }
        answer[degree] = coordinate * previous;
    }
}

/// Append a positive basis expansion in the additive inverse-energy Compton
/// loss coordinate. Each coefficient is the integrated yield assigned to one
/// normalized Bernstein basis, so their sum is the component yield.
template <std::size_t N>
inline std::array<double, N> cubic_bspline_values(double coordinate)
{
    static_assert(N > 3, "cubic B-spline basis requires at least 4 terms");
    constexpr std::size_t degree = 3;
    constexpr std::size_t knot_count = N + degree + 1;
    std::array<double, knot_count> knots{};
    for (std::size_t index = 0; index < knot_count; ++index) {
        if (index <= degree)
            knots[index] = 0.0;
        else if (index >= N)
            knots[index] = 1.0;
        else
            knots[index] =
                static_cast<double>(index - degree)
                / static_cast<double>(N - degree);
    }
    std::array<double, N> answer{};
    if (coordinate >= 1.0) {
        answer[N - 1] = 1.0;
        return answer;
    }
    constexpr std::size_t work_size = N + degree;
    std::array<double, work_size> previous{};
    for (std::size_t index = 0; index + 1 < knot_count; ++index) {
        if (knots[index] <= coordinate
            && coordinate < knots[index + 1]) {
            previous[index] = 1.0;
        }
    }
    for (std::size_t order = 1; order <= degree; ++order) {
        std::array<double, work_size> current{};
        const std::size_t count = work_size - order;
        for (std::size_t index = 0; index < count; ++index) {
            const double left_denominator =
                knots[index + order] - knots[index];
            if (left_denominator > 0.0) {
                current[index] += previous[index]
                    * (coordinate - knots[index])
                    / left_denominator;
            }
            const double right_denominator =
                knots[index + order + 1] - knots[index + 1];
            if (right_denominator > 0.0) {
                current[index] += previous[index + 1]
                    * (knots[index + order + 1] - coordinate)
                    / right_denominator;
            }
        }
        previous = current;
    }
    std::copy_n(previous.begin(), N, answer.begin());
    return answer;
}

template <std::size_t N>
inline std::array<double, N> cubic_bspline_integrals(
    double coordinate_lo, double coordinate_hi)
{
    std::array<double, N> answer{};
    if (!(coordinate_hi > coordinate_lo))
        return answer;
    constexpr int degree = 3;
    constexpr int intervals = static_cast<int>(N) - degree;
    constexpr double inverse_sqrt_three = 0.57735026918962576451;
    double segment_lo = coordinate_lo;
    for (int boundary = 1; boundary <= intervals; ++boundary) {
        const double knot = static_cast<double>(boundary) / intervals;
        const double segment_hi = std::min(coordinate_hi, knot);
        if (segment_hi > segment_lo) {
            const double midpoint = 0.5 * (segment_lo + segment_hi);
            const double half_width = 0.5 * (segment_hi - segment_lo);
            const auto left = cubic_bspline_values<N>(
                midpoint - half_width * inverse_sqrt_three);
            const auto right = cubic_bspline_values<N>(
                midpoint + half_width * inverse_sqrt_three);
            for (std::size_t index = 0; index < N; ++index) {
                answer[index] += half_width
                    * (left[index] + right[index]);
            }
        }
        if (knot >= coordinate_hi)
            break;
        segment_lo = std::max(segment_lo, knot);
    }
    for (std::size_t index = 0; index < N; ++index) {
        const int knot_lo = static_cast<int>(index);
        const int knot_hi = static_cast<int>(index) + degree + 1;
        const auto knot_value = [&](int knot_index) {
            if (knot_index <= degree)
                return 0.0;
            if (knot_index >= static_cast<int>(N))
                return 1.0;
            return static_cast<double>(knot_index - degree)
                / intervals;
        };
        const double normalization =
            (knot_value(knot_hi) - knot_value(knot_lo))
            / (degree + 1);
        answer[index] /= normalization;
    }
    return answer;
}

template <typename Scalar, std::size_t N>
inline void append_warped_bernstein_density(
    SpectralMeasure<Scalar>& result,
    const std::vector<double>& edges,
    double incident_energy_keV,
    double inverse_energy_warp_scale,
    double absorption_envelope_strength,
    double absorption_envelope_power,
    double absorption_envelope_max_ratio,
    bool use_cubic_bspline,
    const Scalar& incident_photoelectric_tau,
    const std::array<Scalar, N>& coefficients,
    SpectralComponent component)
{
    const int degree = static_cast<int>(N) - 1;
    const auto basis_integrals =
        [&](double warp_lo, double warp_hi) {
            if (use_cubic_bspline)
                return cubic_bspline_integrals<N>(warp_lo, warp_hi);
            std::array<double, N> answer{};
            for (std::size_t k = 0; k < N; ++k) {
                answer[k] = bernstein_integral(
                    degree, static_cast<int>(k),
                    warp_lo, warp_hi);
            }
            return answer;
        };
    if (!(absorption_envelope_strength > 0.0)) {
        for (std::size_t b = 0; b + 1 < edges.size(); ++b) {
            const double energy_lo = std::max(0.0, edges[b]);
            const double energy_hi =
                std::min(edges[b + 1], incident_energy_keV);
            if (!(energy_hi > energy_lo))
                continue;
            const double warp_lo = inverse_energy_loss_warp(
                energy_hi, incident_energy_keV,
                inverse_energy_warp_scale);
            const double warp_hi = inverse_energy_loss_warp(
                energy_lo, incident_energy_keV,
                inverse_energy_warp_scale);
            if (!(warp_hi > warp_lo))
                continue;
            const auto fractions = basis_integrals(warp_lo, warp_hi);
            Scalar integral = Scalar(0);
            for (std::size_t k = 0; k < N; ++k) {
                integral += coefficients[k]
                    * fractions[k];
            }
            result.continuous.push_back(
                {edges[b], edges[b + 1], integral, component});
        }
        return;
    }
    std::array<Scalar, N> basis_normalization{};
    for (std::size_t b = 0; b + 1 < edges.size(); ++b) {
        const double energy_lo = std::max(0.0, edges[b]);
        const double energy_hi =
            std::min(edges[b + 1], incident_energy_keV);
        if (!(energy_hi > energy_lo))
            continue;
        const double warp_lo = inverse_energy_loss_warp(
            energy_hi, incident_energy_keV,
            inverse_energy_warp_scale);
        const double warp_hi = inverse_energy_loss_warp(
            energy_lo, incident_energy_keV,
            inverse_energy_warp_scale);
        const double center = 0.5 * (energy_lo + energy_hi);
        const double loss_factor = std::pow(std::min(
            absorption_envelope_max_ratio,
            incident_energy_keV / center),
            absorption_envelope_power) - 1.0;
        const Scalar envelope = scalar_exp(
            -absorption_envelope_strength
            * incident_photoelectric_tau * loss_factor);
        const auto fractions = basis_integrals(warp_lo, warp_hi);
        for (std::size_t k = 0; k < N; ++k) {
            basis_normalization[k] += envelope
                * fractions[k];
        }
    }
    for (std::size_t b = 0; b + 1 < edges.size(); ++b) {
        const double energy_lo = std::max(0.0, edges[b]);
        const double energy_hi =
            std::min(edges[b + 1], incident_energy_keV);
        if (!(energy_hi > energy_lo))
            continue;
        const double warp_lo = inverse_energy_loss_warp(
            energy_hi, incident_energy_keV,
            inverse_energy_warp_scale);
        const double warp_hi = inverse_energy_loss_warp(
            energy_lo, incident_energy_keV,
            inverse_energy_warp_scale);
        if (!(warp_hi > warp_lo))
            continue;
        const double center = 0.5 * (energy_lo + energy_hi);
        const double loss_factor = std::pow(std::min(
            absorption_envelope_max_ratio,
            incident_energy_keV / center),
            absorption_envelope_power) - 1.0;
        const Scalar envelope = scalar_exp(
            -absorption_envelope_strength
            * incident_photoelectric_tau * loss_factor);
        const auto fractions = basis_integrals(warp_lo, warp_hi);
        Scalar integral = Scalar(0);
        for (std::size_t k = 0; k < N; ++k) {
            integral += coefficients[k] * envelope
                * fractions[k]
                / (basis_normalization[k]
                   + Scalar(std::numeric_limits<double>::min()));
        }
        result.continuous.push_back(
            {edges[b], edges[b + 1], integral, component});
    }
}

/// Append a positive cubic B-spline expansion in logarithmic photon energy.
/// Each coefficient is the integrated yield assigned to one normalized basis.
/// This coordinate resolves the approximately 1/k bremsstrahlung spectrum
/// without storing or interpolating a histogram.
template <typename Scalar, std::size_t N>
inline void append_log_energy_bspline_density(
    SpectralMeasure<Scalar>& result,
    const std::vector<double>& edges,
    double minimum_energy_keV,
    double maximum_energy_keV,
    const std::array<Scalar, N>& coefficients,
    SpectralComponent component)
{
    if (!(maximum_energy_keV > minimum_energy_keV)
        || !(minimum_energy_keV > 0.0))
        return;
    const double inverse_log_span =
        1.0 / std::log(maximum_energy_keV / minimum_energy_keV);
    const auto coordinate = [&](double energy_keV) {
        if (energy_keV <= minimum_energy_keV)
            return 0.0;
        if (energy_keV >= maximum_energy_keV)
            return 1.0;
        return std::log(energy_keV / minimum_energy_keV)
            * inverse_log_span;
    };
    for (std::size_t b = 0; b + 1 < edges.size(); ++b) {
        const double energy_lo =
            std::max(minimum_energy_keV, edges[b]);
        const double energy_hi =
            std::min(maximum_energy_keV, edges[b + 1]);
        if (!(energy_hi > energy_lo))
            continue;
        const auto fractions = cubic_bspline_integrals<N>(
            coordinate(energy_lo), coordinate(energy_hi));
        Scalar integral = Scalar(0);
        for (std::size_t k = 0; k < N; ++k)
            integral += coefficients[k] * fractions[k];
        result.continuous.push_back(
            {edges[b], edges[b + 1], integral, component});
    }
}

/// Apply a precomputed smooth photon-absorption envelope to every logarithmic
/// B-spline basis and renormalize each basis over the supplied full-support
/// binning. This preserves the fitted component yield while moving the sharp,
/// material-dependent low-energy escape cutoff out of the learned shape.
template <typename Scalar, std::size_t N>
inline void append_attenuated_log_energy_bspline_density(
    SpectralMeasure<Scalar>& result,
    const std::vector<double>& edges,
    double minimum_energy_keV,
    double maximum_energy_keV,
    const std::vector<Scalar>& bin_escape_envelope,
    const std::array<Scalar, N>& coefficients,
    SpectralComponent component)
{
    if (bin_escape_envelope.size() + 1 != edges.size())
        throw std::invalid_argument(
            "absorption envelope must match output bins");
    if (!(maximum_energy_keV > minimum_energy_keV)
        || !(minimum_energy_keV > 0.0))
        return;
    const double inverse_log_span =
        1.0 / std::log(maximum_energy_keV / minimum_energy_keV);
    const auto coordinate = [&](double energy_keV) {
        if (energy_keV <= minimum_energy_keV)
            return 0.0;
        if (energy_keV >= maximum_energy_keV)
            return 1.0;
        return std::log(energy_keV / minimum_energy_keV)
            * inverse_log_span;
    };
    std::vector<std::array<double, N>> fractions(edges.size() - 1);
    std::array<Scalar, N> normalization{};
    for (std::size_t b = 0; b + 1 < edges.size(); ++b) {
        const double energy_lo =
            std::max(minimum_energy_keV, edges[b]);
        const double energy_hi =
            std::min(maximum_energy_keV, edges[b + 1]);
        if (!(energy_hi > energy_lo))
            continue;
        fractions[b] = cubic_bspline_integrals<N>(
            coordinate(energy_lo), coordinate(energy_hi));
        for (std::size_t k = 0; k < N; ++k) {
            normalization[k] +=
                bin_escape_envelope[b] * fractions[b][k];
        }
    }
    for (std::size_t b = 0; b + 1 < edges.size(); ++b) {
        const double energy_lo =
            std::max(minimum_energy_keV, edges[b]);
        const double energy_hi =
            std::min(maximum_energy_keV, edges[b + 1]);
        if (!(energy_hi > energy_lo))
            continue;
        Scalar integral = Scalar(0);
        for (std::size_t k = 0; k < N; ++k) {
            integral += coefficients[k] * bin_escape_envelope[b]
                * fractions[b][k]
                / (normalization[k]
                   + Scalar(std::numeric_limits<double>::min()));
        }
        result.continuous.push_back(
            {edges[b], edges[b + 1], integral, component});
    }
}

/// Append a positive cubic B-spline expansion on a finite linear-energy
/// interval. Used for the moving [511, 511+T] direct in-flight-annihilation
/// support, which is approximately uniform before photon attenuation.
template <typename Scalar, std::size_t N>
inline void append_linear_energy_bspline_density(
    SpectralMeasure<Scalar>& result,
    const std::vector<double>& edges,
    double minimum_energy_keV,
    double maximum_energy_keV,
    const std::array<Scalar, N>& coefficients,
    SpectralComponent component)
{
    if (!(maximum_energy_keV > minimum_energy_keV))
        return;
    const double inverse_span =
        1.0 / (maximum_energy_keV - minimum_energy_keV);
    for (std::size_t b = 0; b + 1 < edges.size(); ++b) {
        const double energy_lo =
            std::max(minimum_energy_keV, edges[b]);
        const double energy_hi =
            std::min(maximum_energy_keV, edges[b + 1]);
        if (!(energy_hi > energy_lo))
            continue;
        const double coordinate_lo =
            (energy_lo - minimum_energy_keV) * inverse_span;
        const double coordinate_hi =
            (energy_hi - minimum_energy_keV) * inverse_span;
        const auto fractions = cubic_bspline_integrals<N>(
            coordinate_lo, coordinate_hi);
        Scalar integral = Scalar(0);
        for (std::size_t k = 0; k < N; ++k)
            integral += coefficients[k] * fractions[k];
        result.continuous.push_back(
            {edges[b], edges[b + 1], integral, component});
    }
}

/// Runtime-size clamped cubic B-spline basis values.
///
/// Same knot layout and Cox-de Boor recursion as cubic_bspline_values<N>, with
/// the coefficient count supplied at run time so one deposited-energy interval
/// can be partitioned into independently sized blocks. Values are `double`
/// throughout: the coordinate depends only on incident energy and bin edges,
/// never on a fitted parameter, so no Jet ever reaches the knot search.
inline void cubic_bspline_values_runtime(
    std::size_t size, double coordinate, double* answer)
{
    constexpr std::size_t degree = 3;
    const std::size_t knot_count = size + degree + 1;
    std::array<double, kDetectorMaxBlockSize + degree + 1> knots{};
    for (std::size_t index = 0; index < knot_count; ++index) {
        if (index <= degree)
            knots[index] = 0.0;
        else if (index >= size)
            knots[index] = 1.0;
        else
            knots[index] = static_cast<double>(index - degree)
                / static_cast<double>(size - degree);
    }
    std::fill_n(answer, size, 0.0);
    if (coordinate >= 1.0) {
        answer[size - 1] = 1.0;
        return;
    }
    const std::size_t work_size = size + degree;
    std::array<double, kDetectorMaxBlockSize + degree> previous{};
    std::array<double, kDetectorMaxBlockSize + degree> current{};
    for (std::size_t index = 0; index + 1 < knot_count; ++index) {
        if (knots[index] <= coordinate && coordinate < knots[index + 1])
            previous[index] = 1.0;
    }
    for (std::size_t order = 1; order <= degree; ++order) {
        current.fill(0.0);
        const std::size_t count = work_size - order;
        for (std::size_t index = 0; index < count; ++index) {
            const double left = knots[index + order] - knots[index];
            if (left > 0.0) {
                current[index] += previous[index]
                    * (coordinate - knots[index]) / left;
            }
            const double right =
                knots[index + order + 1] - knots[index + 1];
            if (right > 0.0) {
                current[index] += previous[index + 1]
                    * (knots[index + order + 1] - coordinate) / right;
            }
        }
        previous = current;
    }
    std::copy_n(previous.begin(), size, answer);
}

/// Runtime-size version of cubic_bspline_integrals<N>.
inline void cubic_bspline_integrals_runtime(
    std::size_t size, double coordinate_lo, double coordinate_hi,
    double* answer)
{
    std::fill_n(answer, size, 0.0);
    if (!(coordinate_hi > coordinate_lo))
        return;
    constexpr int degree = 3;
    const int intervals = static_cast<int>(size) - degree;
    constexpr double inverse_sqrt_three = 0.57735026918962576451;
    std::array<double, kDetectorMaxBlockSize> left{};
    std::array<double, kDetectorMaxBlockSize> right{};
    double segment_lo = coordinate_lo;
    for (int boundary = 1; boundary <= intervals; ++boundary) {
        const double knot = static_cast<double>(boundary) / intervals;
        const double segment_hi = std::min(coordinate_hi, knot);
        if (segment_hi > segment_lo) {
            const double midpoint = 0.5 * (segment_lo + segment_hi);
            const double half_width = 0.5 * (segment_hi - segment_lo);
            cubic_bspline_values_runtime(
                size, midpoint - half_width * inverse_sqrt_three,
                left.data());
            cubic_bspline_values_runtime(
                size, midpoint + half_width * inverse_sqrt_three,
                right.data());
            for (std::size_t index = 0; index < size; ++index)
                answer[index] += half_width * (left[index] + right[index]);
        }
        if (knot >= coordinate_hi)
            break;
        segment_lo = std::max(segment_lo, knot);
    }
    for (std::size_t index = 0; index < size; ++index) {
        const auto knot_value = [&](int knot_index) {
            if (knot_index <= degree)
                return 0.0;
            if (knot_index >= static_cast<int>(size))
                return 1.0;
            return static_cast<double>(knot_index - degree) / intervals;
        };
        const double normalization =
            (knot_value(static_cast<int>(index) + degree + 1)
             - knot_value(static_cast<int>(index)))
            / (degree + 1);
        answer[index] /= normalization;
    }
}

/// Deposited-energy block partition for one detector continuum component.
///
/// Block boundaries are physical landmarks that move with incident energy; the
/// coefficient count per block is fixed per component, so coefficient k always
/// carries the same physical meaning and the incident-energy interpolation of
/// the coefficient vector stays well posed.
struct DetectorBlockLayout {
    std::size_t count = 0;
    std::array<double, kDetectorMaxBlocks + 1> bounds{};
    std::array<std::size_t, kDetectorMaxBlocks> sizes{};

    std::size_t total_size() const
    {
        std::size_t answer = 0;
        for (std::size_t index = 0; index < count; ++index)
            answer += sizes[index];
        return answer;
    }
};

/// Blocks for the ordinary Compton continua: [0, Compton edge] and
/// [Compton edge, Ein]. The free-electron edge is a near-step in the true
/// distribution; giving each side its own basis lets the expansion represent
/// that step instead of trading a below-edge deficit against an edge spike.
inline DetectorBlockLayout detector_compton_block_layout(
    double incident_energy_keV, std::size_t below_size,
    std::size_t above_size)
{
    constexpr double electron_mass_keV = 510.99895;
    const double edge = incident_energy_keV
        * (1.0 - 1.0
           / (1.0 + 2.0 * incident_energy_keV / electron_mass_keV));
    DetectorBlockLayout layout;
    layout.count = 2;
    layout.bounds[0] = 0.0;
    layout.bounds[1] = edge;
    layout.bounds[2] = incident_energy_keV;
    layout.sizes[0] = below_size;
    layout.sizes[1] = above_size;
    return layout;
}

/// Blocks for the pair-history continua, separated at the sorted landmarks
/// 0, 511, Ein-1022, Ein-511, Ein. Sorting keeps the partition well ordered
/// through the Ein = 1533 keV crossing of 511 and Ein-1022; the interval
/// between those two landmarks shrinks to zero width there and contributes
/// nothing, so the emitted spectrum stays continuous in incident energy.
inline DetectorBlockLayout detector_pair_block_layout(
    double incident_energy_keV,
    const std::array<std::size_t, 4>& sizes)
{
    constexpr double electron_mass_keV = 510.99895;
    std::array<double, 5> points{
        0.0,
        incident_energy_keV - 2.0 * electron_mass_keV,
        electron_mass_keV,
        incident_energy_keV - electron_mass_keV,
        incident_energy_keV};
    std::sort(points.begin(), points.end());
    DetectorBlockLayout layout;
    layout.count = 4;
    for (std::size_t index = 0; index < 5; ++index) {
        layout.bounds[index] =
            std::clamp(points[index], 0.0, incident_energy_keV);
    }
    for (std::size_t index = 0; index < 4; ++index)
        layout.sizes[index] = sizes[index];
    return layout;
}

/// Append a continuum expanded in independent per-block cubic B-splines.
///
/// Exactly one bin is emitted per requested edge interval that the component
/// covers; block contributions to the same bin are summed. Each block's basis
/// functions integrate to one over the block, so the sum of all coefficients
/// remains the component yield and the simplex normalization applied by the
/// caller is unchanged.
template <typename Scalar, std::size_t MaxN>
inline void append_detector_block_density(
    SpectralMeasure<Scalar>& result,
    const std::vector<double>& edges,
    const DetectorBlockLayout& layout,
    const std::array<Scalar, MaxN>& coefficients,
    SpectralComponent component)
{
    const double support = layout.bounds[layout.count];
    if (!(support > 0.0) || layout.count == 0)
        return;
    std::array<double, kDetectorMaxBlockSize> fractions{};
    for (std::size_t bin = 0; bin + 1 < edges.size(); ++bin) {
        const double energy_lo = std::max(0.0, edges[bin]);
        const double energy_hi = std::min(support, edges[bin + 1]);
        if (!(energy_hi > energy_lo))
            continue;
        Scalar integral = Scalar(0);
        std::size_t offset = 0;
        for (std::size_t block = 0; block < layout.count; ++block) {
            const std::size_t size = layout.sizes[block];
            const double low = layout.bounds[block];
            const double high = layout.bounds[block + 1];
            const double width = high - low;
            if (width > 1e-9 && offset + size <= MaxN) {
                const double overlap_lo = std::max(energy_lo, low);
                const double overlap_hi = std::min(energy_hi, high);
                if (overlap_hi > overlap_lo) {
                    cubic_bspline_integrals_runtime(
                        size, (overlap_lo - low) / width,
                        (overlap_hi - low) / width, fractions.data());
                    for (std::size_t index = 0; index < size; ++index) {
                        integral += coefficients[offset + index]
                            * fractions[index];
                    }
                }
            }
            offset += size;
        }
        result.continuous.push_back(
            {edges[bin], edges[bin + 1], integral, component});
    }
}

inline bool supported_detector_basis_size(std::size_t size)
{
    switch (size) {
    case 12:
    case 16:
    case 24:
    case 32:
    case 48:
    case 64:
    case 96:
        return true;
    default:
        return false;
    }
}

template <std::size_t N, typename Scalar, std::size_t MaxN>
inline std::array<Scalar, N> detector_coefficient_prefix(
    const std::array<Scalar, MaxN>& coefficients)
{
    static_assert(N <= MaxN);
    std::array<Scalar, N> answer{};
    std::copy_n(coefficients.begin(), N, answer.begin());
    return answer;
}

/// Dispatch a compact, per-component detector basis without heap allocation.
/// Atom-only components use scalar-yield interpolation and never enter here.
template <typename Scalar, std::size_t MaxN>
inline void append_detector_linear_density(
    SpectralMeasure<Scalar>& result,
    const std::vector<double>& edges,
    double minimum_energy_keV,
    double maximum_energy_keV,
    const std::array<Scalar, MaxN>& coefficients,
    std::size_t active_size,
    SpectralComponent component)
{
#define CEELO_APPEND_LINEAR_CASE(N)                                      \
    case N: {                                                           \
        const auto active = detector_coefficient_prefix<N>(coefficients); \
        append_linear_energy_bspline_density(                           \
            result, edges, minimum_energy_keV, maximum_energy_keV,      \
            active, component);                                        \
        return;                                                        \
    }
    switch (active_size) {
        CEELO_APPEND_LINEAR_CASE(12)
        CEELO_APPEND_LINEAR_CASE(16)
        CEELO_APPEND_LINEAR_CASE(24)
        CEELO_APPEND_LINEAR_CASE(32)
        CEELO_APPEND_LINEAR_CASE(48)
        CEELO_APPEND_LINEAR_CASE(64)
        CEELO_APPEND_LINEAR_CASE(96)
    default:
        throw std::invalid_argument("unsupported detector basis size");
    }
#undef CEELO_APPEND_LINEAR_CASE
}

} // namespace detail

struct PartialMassAttenuation {
    double pe = 0.0;
    double cs = 0.0;
    double rs = 0.0;
    double pp = 0.0;
};

/// Double-only elemental mass attenuation coefficients (cm^2/g).
PartialMassAttenuation element_mass_cross_sections(int Z, double energy_keV);

/// Double-only charged-particle properties used as smooth mechanistic features.
double element_csda_range_g_cm2(
    int Z, double kinetic_energy_keV, bool is_positron);
double element_inverse_radiation_length_cm2_g(int Z);

/// Fluorescence families of the source-escape model: K-alpha, K-beta, L.
constexpr int kFluorescenceFamilies = 3;
/// Fluorescence lines below this are neither modelled nor tallied (the photon
/// tables and the model start at 10 keV).
constexpr double kFluorescenceMinimumKeV = 10.0;

/// One fluorescence line of an element as the source-escape fluorescence
/// components use it: its family (0 K-alpha, 1 K-beta, 2 L) and its yield per
/// K vacancy (`per_k`) and per L1/L2/L3 vacancy (`per_l`). The yields follow
/// the relaxation the Monte Carlo samples: a K line by intensity times omega_K;
/// the L vacancy a K-alpha leaves; an L vacancy radiating with Coster-Kronig
/// transfer (lines at or above the 10 keV cut only).
struct FluorescenceLineYield {
    double energy_keV = 0.0;
    int family = 0;
    double per_k = 0.0;
    std::array<double, 3> per_l{};
};

/// An element's fluorescence lines (Z = 1..92), compressed: within a family,
/// lines under 0.5% of the family's yield fold into the nearest kept line and
/// lines within 0.05 keV merge. Cached; double-only.
struct ElementFluorescence {
    double k_edge_keV = 0.0;   ///< 0 without K fluorescence data
    double omega_K = 0.0;      ///< K fluorescence yield
    std::array<double, 3> l_edge_keV{};  ///< L1, L2, L3; 0 without L lines
    std::vector<FluorescenceLineYield> lines;
};
const ElementFluorescence& element_fluorescence(int Z);

/// {K, L1, L2, L3} shares of the element's photoelectric cross section at
/// `energy_keV` (the L shares only where they are tabulated, Z >= 69).
std::array<double, 4> element_photoelectric_vacancy_shares(int Z, double energy_keV);

/// The photoelectron's kinetic energy after photoabsorption at `energy_keV` in
/// pure element Z, as {energy, probability} for the three outcomes the source
/// walk samples: a K-shell absorption leaves E - E_K, a retained-L one
/// E - E_Ls (merged at the share-weighted mean binding), any other the whole
/// E. Probabilities sum to 1.
std::array<std::pair<double, double>, 3> element_photoelectron_energies(
    int Z, double energy_keV);

/// Electron-impact ionization K x-rays an electron of `kinetic_energy_keV`
/// makes slowing down to rest in pure element Z (EEDL K-shell cross section x
/// omega_K over the CSDA path; zero above Z = 92).
double element_impact_xray_yield(int Z, double kinetic_energy_keV);

constexpr std::size_t kMultipleScatterBasisSize = 32;
constexpr std::size_t kMultipleScatterYieldBasisSize = 6;

/// Widest scattered-RBF feature vector any source component uses. The "raw"
/// coordinate is 3-D (log incident energy, Z/92, log total optical depth); the
/// secondary surrogate's "process" coordinate is 7-D, appending the
/// parent-process optical depth, the descendant's own absorption depth, the
/// shield thickness in radiation lengths, and the hydrogen mass fraction.
/// Components carry their own `feature_dimension`.
constexpr std::size_t kSourceRbfMaxFeatures = 7;
constexpr std::size_t kDetectorResponseBasisSize = 96;
/// Serialized detector-response payload format version. Only this version is
/// read; it is not incremented before the first release.
constexpr std::uint32_t kDetectorPayloadVersion = 1;
constexpr std::size_t kDetectorRateEnergyBasisSize = 20;
constexpr std::size_t kDetectorRateAngleBasisSize = 12;
// Upper bounds for the stored detected-rate surface. The surface's actual
// dimensions travel in the payload; these only size the fixed arrays.
constexpr std::size_t kDetectorRateMaxEnergyBasisSize = 32;
constexpr std::size_t kDetectorRateMaxAngleBasisSize = 32;
constexpr std::size_t kDetectorRateSurfaceSize =
    kDetectorRateEnergyBasisSize * kDetectorRateAngleBasisSize;

struct MultipleScatterYieldNode {
    std::array<double, 2> feature{};
    std::array<double, kMultipleScatterYieldBasisSize> coefficients{};
};

struct MultipleScatterYieldModel {
    std::uint32_t degree = 5;
    std::array<double, 2> bandwidth{0.1, 0.1};
    double minimum_tau = 0.2;
    double maximum_tau = 12.0;
    std::vector<MultipleScatterYieldNode> nodes;
};

struct MultipleScatterRbfNode {
    /// Only the first `feature_dimension` entries of the owning component are
    /// populated or serialized; the remainder stay zero.
    std::array<double, kSourceRbfMaxFeatures> feature{};
    double log_conditional_yield = 0.0;
    std::array<double, kMultipleScatterBasisSize>
        shape_coefficients{};
};

struct MultipleScatterRbfComponent {
    using value_type = MultipleScatterRbfNode;
    std::array<double, kSourceRbfMaxFeatures>
        yield_bandwidth{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    std::array<double, kSourceRbfMaxFeatures>
        shape_bandwidth{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    /// Feature coordinate: 3 = the "raw" (log E, Z/92, log tau) coordinate
    /// used by the higher-Compton and photonuclear components; 7 = the
    /// secondary surrogate's process coordinate (see kSourceRbfMaxFeatures).
    std::uint32_t feature_dimension = 3;
    /// Nonzero replaces the normalized-Gaussian kernel MEAN (Nadaraya-Watson)
    /// with a ridge-stabilized local-LINEAR regression along
    /// `local_linear_feature_index`. A weighted average cannot follow a slope,
    /// so the kernel mean over-predicts the low side wherever a field varies
    /// steeply between training nodes. Zero selects the kernel mean.
    double local_linear_ridge = 0.0;
    std::uint32_t local_linear_feature_index = 2;
    /// Deepest optical depth this component's node cloud actually populates,
    /// derived from the nodes at load time (never serialized, so it cannot
    /// disagree with them). Zero disables the tail decay envelope.
    double tail_decay_tau = 0.0;
    /// Copied from SourceSecondarySurrogatePayload::tail_decay_kappa at load
    /// so evaluation depends only on the component. Zero disables the decay.
    double tail_decay_kappa = 0.0;
    std::vector<MultipleScatterRbfNode> nodes;

    /// Recompute `tail_decay_tau` from the node cloud. Feature index 2 stores
    /// log1p(tau)/log(25) in both the raw and process coordinates.
    void refresh_tail_decay_tau() {
        double maximum_feature = 0.0;
        for (const auto& node : nodes)
            maximum_feature = std::max(maximum_feature, node.feature[2]);
        tail_decay_tau =
            nodes.empty()
                ? 0.0
                : std::exp(maximum_feature * std::log(25.0)) - 1.0;
    }
};

/// One charged-particle payload component: a tensor-product clamped-uniform
/// cubic B-spline field for the log conditional yield plus a rank-R softmax
/// mixture of nonnegative archetype shapes whose log-odds are small tensor
/// fields. The Python fitter (tools/full_spectrum/
/// fit_charged_particle_source_structured.py) and this evaluator implement the
/// same definition; their parity gate is <= 1e-9 relative.
///
/// Records (axes = 3: log-KE warp, Z/92, warped range ratio):
///
///     <name>_domain 3 <lo0> <hi0> <lo1> <hi1> <lo2> <hi2>
///     <name>_yield_grid 3 <n0> <n1> <n2>
///     <name>_yield_coefficients <n0*n1*n2 values, last axis fastest>
///     <name>_shape_rank <R>
///     <name>_shape_basis <32 values>            (R times, each sums to 1)
///     <name>_weight_grid 3 <m0> <m1> <m2>       (only when R >= 2)
///     <name>_weight_log_odds <m0*m1*m2 values>  (R-1 times, only when R >= 2)
///
/// Constraints (the reader throws otherwise): 4 <= n_a, m_a <= 32; coefficient
/// counts match the grid products; 0 <= R <= 16; all values finite; archetypes
/// nonnegative.
///
/// Evaluation: per axis u = (q - lo) / (hi - lo), freeze-clamped to [0, 1], so
/// outside the trained domain every field is constant and the mechanistic base
/// carries the physical trend. The basis is the open-knot clamped-uniform cubic
/// B-spline of detail::cubic_bspline_values. conditional_yield =
/// exp(yield field). R == 0 is yield-only (uniform shape); R == 1 uses
/// archetype 0; R >= 2 mixes archetypes with softmax weights of the log-odds
/// fields (the last log-odds is fixed at 0). A component's coefficients are
/// base_yield * conditional_yield * shape.
struct ChargedParticleSplineField {
    std::vector<double> domain;               ///< (lo, hi) per axis
    std::vector<std::uint32_t> yield_grid;    ///< basis size per axis, 4..32
    std::vector<double> yield_coefficients;   ///< flattened, last axis fastest
    std::uint32_t shape_rank = 0;             ///< 0 = yield-only component
    std::vector<std::array<double, kMultipleScatterBasisSize>> shape_basis;
    std::vector<std::uint32_t> weight_grid;   ///< present when shape_rank >= 2
    std::vector<std::vector<double>> weight_log_odds;  ///< shape_rank - 1 fields
    bool enabled() const { return !yield_coefficients.empty(); }
};

constexpr std::size_t kChargedSplineMaxBasis = 32;

namespace detail {

/// Runtime-size Scalar-templated twin of cubic_bspline_values<N> with the
/// identical clamped-uniform open-knot convention. Jet-safe: the interval
/// membership branches compare the coordinate's value part only, and the
/// clamped cubic spline is C^2 in the coordinate, so one-sided local
/// polynomial derivatives agree at every knot; callers pre-freeze the
/// coordinate to [0, 1].
template <typename Scalar>
inline void charged_cubic_bspline_values_runtime(
    std::size_t size, const Scalar& coordinate, Scalar* answer)
{
    constexpr std::size_t degree = 3;
    const std::size_t knot_count = size + degree + 1;
    const auto knot = [&](std::size_t index) {
        if (index <= degree)
            return 0.0;
        if (index >= size)
            return 1.0;
        return static_cast<double>(index - degree)
            / static_cast<double>(size - degree);
    };
    for (std::size_t index = 0; index < size; ++index)
        answer[index] = Scalar(0);
    if (coordinate >= Scalar(1)) {
        answer[size - 1] = Scalar(1);
        return;
    }
    std::array<Scalar, kChargedSplineMaxBasis + degree> previous{};
    for (std::size_t index = 0; index + 1 < knot_count; ++index) {
        if (Scalar(knot(index)) <= coordinate
            && coordinate < Scalar(knot(index + 1)))
            previous[index] = Scalar(1);
    }
    for (std::size_t order = 1; order <= degree; ++order) {
        std::array<Scalar, kChargedSplineMaxBasis + degree> current{};
        const std::size_t count = size + degree - order;
        for (std::size_t index = 0; index < count; ++index) {
            const double left_denominator =
                knot(index + order) - knot(index);
            if (left_denominator > 0.0) {
                current[index] += previous[index]
                    * ((coordinate - Scalar(knot(index)))
                       / Scalar(left_denominator));
            }
            const double right_denominator =
                knot(index + order + 1) - knot(index + 1);
            if (right_denominator > 0.0) {
                current[index] += previous[index + 1]
                    * ((Scalar(knot(index + order + 1)) - coordinate)
                       / Scalar(right_denominator));
            }
        }
        previous = current;
    }
    for (std::size_t index = 0; index < size; ++index)
        answer[index] = previous[index];
}

} // namespace detail

struct DetectorResponseRbfNode {
    std::array<double, 3> feature{};
    double log_conditional_yield = 0.0;
    std::array<double, kDetectorResponseBasisSize> shape_coefficients{};
};

struct DetectorResponseRbfComponent {
    using value_type = DetectorResponseRbfNode;
    std::uint32_t basis_size = kDetectorResponseBasisSize;
    std::array<double, 3> yield_bandwidth{0.1, 0.1, 0.1};
    std::array<double, 3> shape_bandwidth{0.1, 0.1, 0.1};
    std::vector<DetectorResponseRbfNode> nodes;
    // Block partition of the deposited-energy support. When the first entry
    // is non-zero the coefficient vector is a concatenation of independent
    // per-block cubic B-spline expansions whose boundaries are the
    // component's moving physical landmarks, and the sizes must sum to
    // basis_size. Zero means no partition (the plain linear basis of the
    // photoelectric-escape continuum); the Compton and pair continua require
    // one.
    std::array<std::uint32_t, detail::kDetectorMaxBlocks> block_sizes{};

    bool uses_blocks() const { return block_sizes[0] != 0; }
};

/// One passive layer the rate-surface anchor attenuates through.
struct AnchorLayer {
    std::vector<MaterialComponent> composition;  ///< {Z, mass fraction}
    double density_g_cm3 = 0.0;
    double thickness_cm = 0.0;
};

/// Deterministic surface-entry anchor for the detected-rate surface. Below
/// ~150 keV the crystal is a surface absorber, so the detected fraction is the
/// passive transmission averaged over a parallel beam's ENTRY POINTS.
/// Averaging exp(-tau) rather than the path is what makes the yield rise
/// toward |cos| = 1 and produces the non-monotone angular shape that a
/// low-order Bernstein surface cannot interpolate from 10-keV-spaced nodes.
///
/// P0(E, cos theta) is normalized PER PHOTON CROSSING THE BEAM ENVELOPE -- the
/// cylinder the detector campaign aims its rays at
/// (Geometry::outer_bounding_radius / outer_z_extent) -- so it is on the same
/// footing as the truth. The entry surfaces are the face the beam meets
/// (front for cos theta >= 0, back otherwise), the curved side and, for back
/// incidence, the bore opening, whose photons reach active crystal with no
/// dead layer. Each face has its own passive stack, of compound materials,
/// as trace_ray() stacks the layers: a layer adds to the back only where it
/// reaches past the shell inside it.
///
/// The anchor carries NO fitted parameters. Its geometry fields are derived
/// from the payload's GeometryDescriptor (set_rate_surface_anchor_geometry)
/// and are never serialized on their own. It is added to the surface's log
/// odds, and the fitter subtracts the identical quantity from its fit target,
/// so the Bernstein surface only has to represent the smooth residual.
/// `enabled = false` reproduces the un-anchored model exactly.
struct DetectorRateSurfaceAnchor {
    bool enabled = false;
    /// The anchor is confined to the turn-on by w(E) = 1/(1+(E/blend)^power).
    double blend_keV = 150.0;
    double blend_power = 2.0;
    std::uint32_t azimuth_quadrature_points = 48;
    double radius_cm = 0.0;
    double length_cm = 0.0;
    double front_dead_cm = 0.0;
    double side_dead_cm = 0.0;
    double back_dead_cm = 0.0;
    /// Radius of the bore's opening in the back face; 0 when there is none.
    double bore_radius_cm = 0.0;
    std::vector<MaterialComponent> crystal_composition;
    double crystal_density_g_cm3 = 0.0;
    std::vector<AnchorLayer> front_layers;
    std::vector<AnchorLayer> side_layers;
    std::vector<AnchorLayer> back_layers;
    double envelope_radius_cm = 0.0;
    double envelope_length_cm = 0.0;
};

/// Fill `anchor`'s geometry fields from a cylindrical detector `geometry`
/// (blend and quadrature settings are left alone). Throws
/// std::invalid_argument for a box crystal or a descriptor that fails
/// GeometryDescriptor::problems().
void set_rate_surface_anchor_geometry(DetectorRateSurfaceAnchor& anchor,
                                      const GeometryDescriptor& geometry);

/// Front-entry optical depth at which the passive-layer landmark sits: the
/// front stack plus the front dead layer at normal incidence of the reference
/// hpge_coax detector at 55 keV, the landmark its feature map was selected
/// with.
inline constexpr double kPassiveLandmarkOpticalDepth = 0.770421;

/// The passive-layer landmark of a detector: the highest energy in
/// [10, 1000] keV at which the front entry's normal-incidence optical depth
/// is kPassiveLandmarkOpticalDepth (10 keV if it never gets that high).
/// Feeds the detector features' passive coordinate 1/(1+(E/landmark)^3).
double passive_landmark_keV(const DetectorRateSurfaceAnchor& anchor);

namespace detail {

inline double total_mass_attenuation(
    const std::vector<MaterialComponent>& composition, double energy_keV)
{
    double mu = 0.0;
    for (const MaterialComponent& component : composition) {
        const PartialMassAttenuation x = element_mass_cross_sections(
            static_cast<int>(component.Z), energy_keV);
        mu += component.mass_fraction * (x.pe + x.cs + x.rs + x.pp);
    }
    return mu;
}

inline double layer_stack_optical_depth(const std::vector<AnchorLayer>& layers,
                                        double energy_keV)
{
    double tau = 0.0;
    for (const AnchorLayer& layer : layers) {
        tau += total_mass_attenuation(layer.composition, energy_keV)
            * layer.density_g_cm3 * layer.thickness_cm;
    }
    return tau;
}

}  // namespace detail

struct DetectorRateSurface {
    bool enabled = false;
    double energy_warp_alpha = 3.0;
    // Angle coordinate reshaping. The polar angle theta/pi is passed through
    // x -> x*(x^2+eps^2)^((q-1)/2), renormalized to +/-1, on x = 2*theta/pi - 1.
    // q = 1 leaves the plain polar angle. Smaller q concentrates basis
    // resolution where the detected rate turns over; eps > 0 bounds the
    // derivative at cos(theta) = 0, which a plain |x|^q would make infinite.
    double angle_warp_exponent = 1.0;
    double angle_warp_epsilon = 0.03;
    std::uint32_t energy_size = kDetectorRateEnergyBasisSize;
    std::uint32_t angle_size = kDetectorRateAngleBasisSize;
    std::array<double,
               kDetectorRateMaxEnergyBasisSize
                   * kDetectorRateMaxAngleBasisSize>
        log_odds_coefficients{};
    DetectorRateSurfaceAnchor anchor;
};

/// Log odds of the detected-rate anchor, `w(E)*logit(P0(E, cos theta))`; see
/// DetectorRateSurfaceAnchor. Cross sections are looked up at the (double)
/// incident energy, so the autodiff type only ever touches the angular
/// factors.
template <typename Scalar>
Scalar rate_surface_anchor_log_odds(
    const DetectorRateSurfaceAnchor& anchor,
    double Ein_keV,
    const Scalar& cos_theta)
{
    if (!anchor.enabled)
        return Scalar(0);
    const double mu_crystal =
        detail::total_mass_attenuation(anchor.crystal_composition, Ein_keV)
        * anchor.crystal_density_g_cm3;
    const double tau_front =
        detail::layer_stack_optical_depth(anchor.front_layers, Ein_keV);
    const double tau_side =
        detail::layer_stack_optical_depth(anchor.side_layers, Ein_keV);
    const double tau_back =
        detail::layer_stack_optical_depth(anchor.back_layers, Ein_keV);
    // Clamp before squaring: a caller whose cosine came from a dot product
    // can hand us 1 + 1e-12, and 1 - cos^2 would then go negative and NaN
    // the entire measure.
    const Scalar shrunk = detail::clamped_cosine(cos_theta);
    // |cos| is smoothed just enough to keep the grazing-incidence derivative
    // finite; the face term it multiplies vanishes there anyway.
    constexpr double kGrazingFloor = 1.0e-6;
    const Scalar cosine = detail::scalar_sqrt(
        shrunk * shrunk + Scalar(kGrazingFloor * kGrazingFloor));
    const Scalar sine = detail::scalar_sqrt(
        Scalar(1) - shrunk * shrunk
        + Scalar(kGrazingFloor * kGrazingFloor));
    const double radius = anchor.radius_cm;
    const double length = anchor.length_cm;
    const Scalar chord = detail::scalar_min(
        Scalar(length) / cosine, Scalar(2.0 * radius) / sine);
    const Scalar absorption =
        Scalar(1) - detail::scalar_exp(Scalar(-mu_crystal) * chord);
    // The face the beam meets: front and back stacks differ, and only the
    // back carries the bore opening. The face weight is proportional to |cos|,
    // so this selection is continuous in value where the two meet.
    const bool back = cos_theta < Scalar(0);
    const double face_passive = back
        ? tau_back + mu_crystal * anchor.back_dead_cm
        : tau_front + mu_crystal * anchor.front_dead_cm;
    const double opening = back ? anchor.bore_radius_cm : 0.0;
    const Scalar face_probability =
        detail::scalar_exp(Scalar(-face_passive) / cosine) * absorption;
    // A photon entering the bore meets active crystal (the bore wall has no
    // dead layer) and crosses about a crystal radius of it.
    const Scalar bore_probability =
        detail::scalar_exp(Scalar(-tau_back) / cosine)
        * Scalar(1.0 - std::exp(-mu_crystal * radius));
    // Illuminated-azimuth average of exp(-tau/(sin*cos phi)) over
    // phi in [-pi/2, pi/2] with the cos(phi) projection weight.
    Scalar side_transmission = Scalar(0);
    const std::size_t points = anchor.azimuth_quadrature_points;
    // The rule depends only on the point count, so build it once: rebuilding
    // it per call would cost ~16 us of a ~17 us call and ~10% of a whole
    // detector evaluation.
    const auto& rule = detail::azimuth_quadrature_rule(points);
    const double side_passive = tau_side + mu_crystal * anchor.side_dead_cm;
    for (std::size_t index = 0; index < points; ++index) {
        const Scalar local = sine * Scalar(rule.cos_phi[index]);
        side_transmission += Scalar(rule.weight[index])
            * detail::scalar_exp(Scalar(-side_passive) / local);
    }
    constexpr double pi = 3.14159265358979323846;
    const Scalar face_area =
        Scalar(pi * (radius * radius - opening * opening)) * cosine;
    const Scalar bore_area = Scalar(pi * opening * opening) * cosine;
    const Scalar side_area = Scalar(2.0 * radius * length) * sine;
    const double envelope_radius = anchor.envelope_radius_cm;
    const Scalar envelope_area =
        Scalar(pi * envelope_radius * envelope_radius) * cosine
        + Scalar(2.0 * envelope_radius * anchor.envelope_length_cm) * sine;
    const Scalar probability =
        (face_area * face_probability + bore_area * bore_probability
         + side_area * side_transmission * absorption)
        / envelope_area;
    constexpr double kFloor = 1.0e-10;
    const Scalar clamped = probability < Scalar(kFloor)
        ? Scalar(kFloor)
        : (probability > Scalar(1.0 - 1e-9) ? Scalar(1.0 - 1e-9)
                                            : probability);
    const double blend = anchor.blend_keV > 0.0
        ? 1.0 / (1.0 + std::pow(Ein_keV / anchor.blend_keV,
                                anchor.blend_power))
        : 1.0;
    return Scalar(blend)
        * detail::scalar_log(clamped / (Scalar(1) - clamped));
}

struct MultipleScatterSurrogatePayload {
    std::uint32_t model_version = 1;
    double inverse_energy_warp_scale = 2.0;
    double absorption_envelope_strength = 0.0;
    double absorption_envelope_power = 3.5;
    double absorption_envelope_max_ratio = 10.0;
    bool use_cubic_bspline = false;
    bool merge_higher_rayleigh = false;
    MultipleScatterYieldModel higher_yield;
    MultipleScatterRbfComponent higher_compton;
    MultipleScatterRbfComponent higher_compton_rayleigh;
};

/// Positive process-resolved residuals for secondary photons born in source
/// shielding. Each component stores a smooth RBF field over
/// (log incident energy, effective Z, log total optical depth). The node
/// coefficients describe a normalized continuous shape; their common scale
/// is supplied by a mechanistic process-onset term and a fitted positive
/// conditional yield. The pair atom uses only the summed coefficient yield.
struct SourceSecondarySurrogatePayload {
    std::uint32_t model_version = 1;
    double inverse_energy_warp_scale = 2.0;
    double absorption_envelope_strength = 0.0;
    double absorption_envelope_power = 3.5;
    double absorption_envelope_max_ratio = 10.0;
    bool use_cubic_bspline = true;
    // The mechanistic bases fold in a closed-form escape probability for
    // components whose creation is otherwise depth-blind (see
    // source_component_bases): the steep escape falloff (355x-1528x between
    // tau = 0.2 and 8 for the bremsstrahlung and in-flight channels) lives in
    // the baseline, so the fitted field carries only a smooth residual.
    /// Mechanistic decay applied to the fitted conditional yield BEYOND the
    /// last optical depth the training bank populates:
    /// `yield *= exp(-kappa * (tau - tau_last_populated))`.
    ///
    /// Without it the kernel field extrapolates FLAT off the deepest node, so
    /// fully opaque shields kept spurious escaping yield (measured up to
    /// 7.4e-2 photons/primary at 250 keV through Pb at 50 g/cm^2, and +3625%
    /// at 59 keV through Al at 50 g/cm^2). The shape envelope
    /// (`absorption_envelope_*`) cannot do this: it is yield-preserving by
    /// construction and only redistributes the spectrum.
    ///
    /// Zero disables it. `tau_last_populated` is derived per component from
    /// its own node cloud
    /// at load time (see `MultipleScatterRbfComponent::tail_decay_tau`), so it
    /// costs no format field and can never disagree with the stored nodes.
    double tail_decay_kappa = 0.0;
    MultipleScatterRbfComponent pair_511_atom;
    MultipleScatterRbfComponent pair_annihilation_continuum;
    MultipleScatterRbfComponent pair_annihilation_in_flight;
    MultipleScatterRbfComponent photoelectron_brems;
    MultipleScatterRbfComponent compton_electron_brems;
    MultipleScatterRbfComponent pair_electron_brems;
    /// Fluorescence atoms (shield K-alpha, K-beta and L x-rays that escape
    /// unscattered): one fitted conditional yield per family, spread over the
    /// family's lines by their mechanistic bases.
    MultipleScatterRbfComponent fluorescence_k_alpha_atom;
    MultipleScatterRbfComponent fluorescence_k_beta_atom;
    MultipleScatterRbfComponent fluorescence_l_atom;

    /// Every component, in the order of detail::SourceComponentIndex, the
    /// payload file, the feature helper and the Python fitter's ALL_COMPONENTS.
    static constexpr std::size_t kComponentCount = 9;
    std::array<MultipleScatterRbfComponent*, kComponentCount> components() {
        return {&pair_511_atom, &pair_annihilation_continuum,
                &pair_annihilation_in_flight, &photoelectron_brems,
                &compton_electron_brems, &pair_electron_brems,
                &fluorescence_k_alpha_atom, &fluorescence_k_beta_atom,
                &fluorescence_l_atom};
    }
    std::array<const MultipleScatterRbfComponent*, kComponentCount> components() const {
        return {&pair_511_atom, &pair_annihilation_continuum,
                &pair_annihilation_in_flight, &photoelectron_brems,
                &compton_electron_brems, &pair_electron_brems,
                &fluorescence_k_alpha_atom, &fluorescence_k_beta_atom,
                &fluorescence_l_atom};
    }
};

struct PhotonuclearCrossSectionRow {
    double energy_keV = 0.0;
    std::array<double, 92> mass_attenuation_cm2_g{};
};

/// Compact Geant4-derived photonuclear feature table. Energy interpolation is
/// piecewise linear; effective Z uses the same smooth cardinal blend as the
/// electromagnetic cross sections. The emitted spectrum remains a functional
/// positive basis expansion, not a stored transport histogram.
struct PhotonuclearSourcePayload {
    std::uint32_t model_version = 1;
    /// Global yield scale on the mechanistic (1 - exp(-tau_pn)) baseline.
    /// Ignored when `secondary` carries nodes: the fitted conditional yield
    /// then absorbs the whole normalization.
    double escaped_yield_scale = 0.85;
    std::array<double, 8> shape{
        0.24894173205315184, 0.7489488311008247,
        0.0, 0.0, 0.0, 0.0, 0.002109436846023685, 0.0};
    std::vector<PhotonuclearCrossSectionRow> cross_sections;
    // ---- Fitted surrogate (inert while `secondary.nodes` is empty) ---------
    /// Fitted (E, Z, tau) conditional yield and warped-B-spline shape for the
    /// escaping photonuclear secondaries. A single fixed shape is ~17% median
    /// total-variation off across the Geant4 bank (p95 37%), driven mostly by
    /// Z and then by incident energy, which is what these nodes carry.
    MultipleScatterRbfComponent secondary;
    double inverse_energy_warp_scale = 1000.0;
    bool use_cubic_bspline = true;
    /// Ratio of the secondaries' mass attenuation to the primary's, used by
    /// the mechanistic escape factor below. The escaping photonuclear yield
    /// falls by three orders of magnitude between tau=0.2 and tau=8; folding
    /// that into the baseline leaves the fitted field varying by ~3x instead
    /// of ~1700x, which is what a kernel-mean interpolant can actually carry.
    double escape_kappa = 0.8;
};

namespace detail {

/// Mean over birth depth u in [0, 1] of the primary reaching u and the
/// secondary escaping the remaining path:
///     integral of exp(-tau u) exp(-kappa tau (1-u)) du
///   = exp(-kappa tau) (1 - exp(-(1-kappa) tau)) / ((1-kappa) tau).
/// Both limits (tau -> 0 and kappa -> 1) are 1 and exp(-tau) respectively;
/// the series guard keeps the value and its derivative finite at tau = 0,
/// which a bare source (zero areal density) hits exactly.
/// Mean over birth depth u in [0, 1] of the primary reaching u and the
/// secondary escaping the rest:
///     integral exp(-tau u) exp(-kappa tau (1-u)) du
///   = (exp(-tau) - exp(-kappa tau)) / (tau (kappa - 1)).
///
/// Written in THAT form rather than the algebraically equal
/// exp(-kappa tau) (1 - exp(-(1-kappa) tau)) / ((1-kappa) tau), because kappa
/// is genuinely much larger than 1 here: a bremsstrahlung child at the
/// eight-downscatter reference energy is far more absorbed than its primary,
/// so kappa reaches tens at low energy in high Z. The second form then
/// evaluates exp(+kappa tau), which overflows to inf and gives inf * 0 = NaN.
/// This one decays to exp(-tau) / (tau (kappa - 1)) instead. Same overflow
/// class as the frozen Chebyshev tail in evaluate_higher_yield.
///
/// One guard covers both removable singularities: kappa -> 1 and tau -> 0
/// both send a = tau (kappa - 1) to zero, where the value is
/// exp(-tau) (1 + a/2 + a^2/6).
template <typename Scalar>
inline Scalar secondary_escape_factor(const Scalar& tau, const Scalar& kappa)
{
    const Scalar a = tau * (kappa - Scalar(1));
    const Scalar decay = scalar_exp(-tau);
    if (a < Scalar(1e-6) && a > Scalar(-1e-6))
        return decay * (Scalar(1) + a / Scalar(2) + a * a / Scalar(6));
    return (decay - scalar_exp(-kappa * tau)) / a;
}

template <typename Scalar>
inline Scalar photonuclear_escape_factor(const Scalar& tau, double kappa)
{
    const double one_minus_kappa = 1.0 - kappa;
    const Scalar x = Scalar(one_minus_kappa) * tau;
    Scalar ratio;
    if (x < Scalar(1e-6) && x > Scalar(-1e-6)) {
        ratio = Scalar(1) - x / Scalar(2)
            + x * x / Scalar(6);
    } else {
        ratio = (Scalar(1) - scalar_exp(-x)) / x;
    }
    return scalar_exp(Scalar(-kappa) * tau) * ratio;
}

} // namespace detail

/// Cross sections for an equivalent element, evaluated with a differentiable
/// positive partition of unity over tabulated integer-Z elements.
class SmoothElementCrossSection {
public:
    explicit SmoothElementCrossSection(double blend_width_Z = 0.42)
        : blend_width_Z_(blend_width_Z) {}

    template <typename Scalar>
    std::array<Scalar, 4> evaluate(double energy_keV,
                                   const Scalar& effective_Z) const;

    /// Two-component mixture: hydrogen mass fraction `f_H` of exact Z=1 plus
    /// `1 - f_H` of the effective-Z blend. Exactly the pure blend at f_H = 0
    /// (the correction term is `f_H * (H - blend)`), smooth and positive in
    /// both fitted scalars. Hydrogen gets its own descriptor because no
    /// effective Z reproduces its electrons per gram.
    template <typename Scalar>
    std::array<Scalar, 4> evaluate(double energy_keV,
                                   const Scalar& effective_Z,
                                   const Scalar& hydrogen_mass_fraction) const;

    /// Smooth effective-Z blend of CeeLo's normalized
    /// Klein-Nishina × incoherent-scattering-function angular densities.
    template <typename Scalar>
    Scalar compton_angular_density(
        double energy_keV,
        const Scalar& effective_Z,
        double cosine) const;

    /// Mixture form: constituents are weighted by their share of Compton
    /// events (mass-fraction × incoherent cross section), which is the exact
    /// scattered-photon angular density of a homogeneous two-component
    /// material.
    template <typename Scalar>
    Scalar compton_angular_density(
        double energy_keV,
        const Scalar& effective_Z,
        const Scalar& hydrogen_mass_fraction,
        double cosine) const;

    /// Smooth effective-Z blend of the exact CeeLo shell-profile probability
    /// for a Doppler-broadened Compton photon to land in an energy bin at a
    /// fixed scattering angle.
    template <typename Scalar>
    Scalar compton_doppler_bin_probability(
        double energy_keV,
        const Scalar& effective_Z,
        double cosine,
        double lower_keV,
        double upper_keV) const;

    /// Mixture form, weighted like `compton_angular_density`.
    template <typename Scalar>
    Scalar compton_doppler_bin_probability(
        double energy_keV,
        const Scalar& effective_Z,
        const Scalar& hydrogen_mass_fraction,
        double cosine,
        double lower_keV,
        double upper_keV) const;

    template <typename Scalar>
    Scalar photonuclear_mass_attenuation(
        double energy_keV,
        const Scalar& effective_Z,
        const PhotonuclearSourcePayload& payload) const;

    /// Mixture form; the hydrogen column is the exact Z=1 table row.
    template <typename Scalar>
    Scalar photonuclear_mass_attenuation(
        double energy_keV,
        const Scalar& effective_Z,
        const Scalar& hydrogen_mass_fraction,
        const PhotonuclearSourcePayload& payload) const;

    double blend_width_Z() const { return blend_width_Z_; }

private:
    struct PhotonuclearInterpolation {
        std::size_t lower = 0;
        std::size_t upper = 0;
        double fraction = 0.0;
        bool below_table = false;
    };
    static PhotonuclearInterpolation photonuclear_interpolation(
        const PhotonuclearSourcePayload& payload, double energy_keV);

    /// Fraction of Compton events on the hydrogen component of the mixture.
    template <typename Scalar>
    Scalar hydrogen_compton_event_share(
        double energy_keV,
        const Scalar& effective_Z,
        const Scalar& hydrogen_mass_fraction) const;

    double blend_width_Z_;
};

struct SourceModelPayload {
    std::uint32_t model_version = 1;
    double min_energy_keV = 10.0;
    double max_energy_keV = 12000.0;
    double min_effective_Z = 1.0;
    double max_effective_Z = 92.0;
    double min_areal_density_g_cm2 = 0.0;
    double max_areal_density_g_cm2 = 60.0;
    double z_blend_width = 0.42;
    double single_escape_scale = 1.0;
    double double_escape_scale = 1.0;
    double multiple_escape_scale = 0.61;
    double compton_rayleigh_scale = 1.0;
    double double_rayleigh_scale = 1.0;
    double multiple_rayleigh_scale = 0.20;
    double elastic_escape_scale = 0.74;
    /// Survival of the Rayleigh-only ("elastic") line against Compton
    /// removal: exp(-elastic_compton_damping * tau_cs). A weighted log-space
    /// least squares over the 1794-cell training bank lands on 0.9994, i.e.
    /// the physically natural 1.0 (a Compton interaction removes the photon
    /// from the elastic line outright); with elastic_escape_scale = 0.978 that
    /// gives 2.7% median / +6.9% at tau > 6.
    double elastic_compton_damping = 1.0;
    double brems_yield_scale = 0.09;
    bool enable_rayleigh = true;
    bool enable_compton_doppler = true;
    bool enable_pair_annihilation = true;
    bool enable_source_electron_brems = false;
    /// Shield K and L x-rays (the fluorescence atoms; emitted only where the
    /// payload carries their fitted yields).
    bool enable_fluorescence = true;
    bool enable_photonuclear = false;
    std::array<double, 8> single_compton_shape{
        0.10, 0.20, 0.55, 1.00, 1.20, 0.90, 0.38, 0.08};
    std::array<double, 8> multiple_compton_shape{
        0.75, 1.15, 1.25, 1.05, 0.72, 0.43, 0.20, 0.06};
    std::array<double, 8> compton_rayleigh_shape{
        0.08, 0.24, 0.62, 1.05, 1.22, 0.88, 0.36, 0.07};
    std::array<double, 8> brems_shape{
        1.40, 1.10, 0.82, 0.58, 0.39, 0.24, 0.12, 0.04};
    MultipleScatterSurrogatePayload multiple_scatter_surrogate;
    SourceSecondarySurrogatePayload secondary_surrogate;
    PhotonuclearSourcePayload photonuclear;
};

class SourceGeomFluxModel {
public:
    explicit SourceGeomFluxModel(SourceModelPayload payload = {})
        : payload_(std::move(payload)),
          cross_sections_(payload_.z_blend_width) {}

    template <typename Scalar>
    SpectralMeasure<Scalar> evaluate(
        double E0_keV,
        const Scalar& areal_density_g_cm2,
        const Scalar& effective_Z,
        const std::vector<double>& output_bins_keV) const;

    /// Two-component shield: hydrogen mass fraction of exact Z=1 plus the
    /// effective-Z remainder (see SmoothElementCrossSection). The 4-argument
    /// form forwards f_H = 0.
    /// The fitted surrogate corrections do not yet carry an f_H coordinate;
    /// only the mechanistic terms respond to it (pilot-bank work quantifies
    /// the residual).
    template <typename Scalar>
    SpectralMeasure<Scalar> evaluate(
        double E0_keV,
        const Scalar& areal_density_g_cm2,
        const Scalar& effective_Z,
        const Scalar& hydrogen_mass_fraction,
        const std::vector<double>& output_bins_keV) const;

    const SourceModelPayload& payload() const { return payload_; }

    std::string serialize() const;
    static SourceGeomFluxModel deserialize(const std::string& text);

private:
    SourceModelPayload payload_;
    SmoothElementCrossSection cross_sections_;
};

/// Compact positive representation of photons escaping from a centered
/// monoenergetic electron or positron in an elemental-equivalent sphere.
///
/// Tensor B-spline fields (ChargedParticleSplineField) interpolate only smooth
/// corrections to mechanistic yield baselines. Continuum shapes are cubic
/// B-splines in log photon energy, which is well matched to bremsstrahlung's
/// approximately 1/k behavior. No source histogram is retained in the payload.
///
/// Serialized format (`CEELO_CHARGED_PARTICLE_SOURCE 1`), written with 17
/// significant digits:
///
///     CEELO_CHARGED_PARTICLE_SOURCE 1
///     validity <min_T> <max_T> <min_photon> <min_Z> <max_Z> <min_rho> <max_rho>
///     features <z_blend_width> <range_feature_scale> <photon_absorption_path_fraction>
///     <component> spline | <component> none      (one record per component)
///     ...                                        (spline fields, see below)
///     END
///
/// Components, in order: electron_brems, positron_brems, positron_containment,
/// positron_511_atom, positron_annihilation_continuum,
/// positron_annihilation_in_flight_direct,
/// positron_annihilation_in_flight_continuum. `none` marks a component with no
/// fitted field; it contributes nothing.
struct ChargedParticleSourcePayload {
    std::uint32_t model_version = 1;
    double min_kinetic_energy_keV = 10.0;
    double max_kinetic_energy_keV = 12000.0;
    double min_photon_energy_keV = 10.0;
    double min_effective_Z = 1.0;
    double max_effective_Z = 92.0;
    double min_areal_density_g_cm2 = 0.0;
    double max_areal_density_g_cm2 = 1000.0;
    double z_blend_width = 0.42;
    double range_feature_scale = 50.0;
    // Validation-selected effective path for the absorption shape envelope.
    // It is a transport-shape parameter, not a literal spherical chord.
    double photon_absorption_path_fraction = 3.0;
    // Fitted correction fields; a component is used iff enabled().
    ChargedParticleSplineField electron_brems_spline;
    ChargedParticleSplineField positron_brems_spline;
    ChargedParticleSplineField positron_containment_spline;
    ChargedParticleSplineField positron_511_atom_spline;
    ChargedParticleSplineField positron_annihilation_continuum_spline;
    ChargedParticleSplineField positron_annihilation_in_flight_direct_spline;
    ChargedParticleSplineField
        positron_annihilation_in_flight_continuum_spline;
    // Shield K-alpha, K-beta and L x-rays (yield-only atoms; the lines share a
    // family's yield by their mechanistic bases).
    ChargedParticleSplineField electron_fluorescence_k_alpha_spline;
    ChargedParticleSplineField electron_fluorescence_k_beta_spline;
    ChargedParticleSplineField electron_fluorescence_l_spline;
    ChargedParticleSplineField positron_fluorescence_k_alpha_spline;
    ChargedParticleSplineField positron_fluorescence_k_beta_spline;
    ChargedParticleSplineField positron_fluorescence_l_spline;
};

class ChargedParticleSourceFluxModel {
public:
    explicit ChargedParticleSourceFluxModel(
        ChargedParticleSourcePayload payload = {})
        : payload_(std::move(payload)),
          cross_sections_(payload_.z_blend_width) {}

    template <typename Scalar>
    SpectralMeasure<Scalar> evaluate(
        double kinetic_energy_keV,
        const Scalar& areal_density_g_cm2,
        const Scalar& effective_Z,
        bool is_positron,
        const std::vector<double>& output_bins_keV) const;

    /// Two-component shield with a hydrogen mass fraction: photon
    /// attenuation mixes exactly; CSDA ranges combine by Bragg additivity
    /// and radiation length by linear X0^-1 additivity. The 5-argument form
    /// forwards f_H = 0. Fitted correction fields carry no f_H coordinate
    /// yet.
    template <typename Scalar>
    SpectralMeasure<Scalar> evaluate(
        double kinetic_energy_keV,
        const Scalar& areal_density_g_cm2,
        const Scalar& effective_Z,
        const Scalar& hydrogen_mass_fraction,
        bool is_positron,
        const std::vector<double>& output_bins_keV) const;

    const ChargedParticleSourcePayload& payload() const {
        return payload_;
    }

    std::string serialize() const;
    static ChargedParticleSourceFluxModel deserialize(
        const std::string& text);

private:
    ChargedParticleSourcePayload payload_;
    SmoothElementCrossSection cross_sections_;
};

enum class DetectorFamily : std::uint8_t {
    NaI3x3,
    CoaxialHPGe
};

/// A point source's finite distance. The detector response is per photon on
/// the envelope of a parallel beam, so `predict`'s far-field acceptance (the
/// envelope's projected area over 4 pi d^2) overstates a point source at d:
/// its rays diverge, and photons entering near the rim leave through the side
/// before they interact (+5% at 1 m, +16% at 30 cm for hpge_coax at 662 keV).
/// A virtual detector disk corrects it: the response at distance d (on axis,
/// from the geometry's reference point -- crystal face or endcap front, as
/// every user distance in CeeLo; the depths are measured from it too) is the
/// far-field response times
///     f(d) = 2 d^2 / (s (s + x)),  x = d + delta,  s = sqrt(x^2 + R^2),
/// the solid angle of a disk of radius R at depth delta, normalized to the far
/// field. Both delta(E) and R(E) are measured per detector by Monte Carlo on a
/// ladder of distances (ceelo_detector_finite_distance): one pair for the atoms
/// that track the incident energy (full energy, escapes), one for everything
/// detected. R falls from near the crystal radius at 100 keV to ~0 above
/// 2 MeV, where the detector looks like a point at depth delta. Empty: none
/// measured, and `predict` is the far-field limit. The depths are measured on
/// the axis and hold to cos theta = 0.5 (hpge_coax, 662 keV, 30 cm - 3 m: within
/// 0.7%); below kFiniteDistanceMinCos `predict` keeps the far-field acceptance,
/// because from the side the source is nearer the crystal's surface than its
/// distance to the face centre says (on-axis depths read -16% side-on at 30 cm,
/// the far field -2%).
inline constexpr double kFiniteDistanceMinCos = 0.5;

struct FiniteDistanceDepth {
    std::vector<double> energies_keV;           ///< ascending
    std::vector<double> full_energy_cm;         ///< delta, tracking atoms
    std::vector<double> full_energy_radius_cm;  ///< R, tracking atoms
    std::vector<double> detected_cm;            ///< delta, detected total
    std::vector<double> detected_radius_cm;     ///< R, detected total
    bool enabled() const { return !energies_keV.empty(); }
};

namespace detail {

inline double finite_distance_factor(double distance_cm, double depth_cm,
                                     double radius_cm)
{
    const double x = distance_cm + depth_cm;
    const double s = std::sqrt(x * x + radius_cm * radius_cm);
    return 2.0 * distance_cm * distance_cm / (s * (s + x));
}

/// A table value (depth or radius) at `energy_keV`: linear in log energy
/// between the table's energies, constant beyond them.
inline double finite_distance_depth(const std::vector<double>& energies_keV,
                                    const std::vector<double>& depth_cm,
                                    double energy_keV)
{
    if (energy_keV <= energies_keV.front())
        return depth_cm.front();
    if (energy_keV >= energies_keV.back())
        return depth_cm.back();
    const auto it = std::upper_bound(energies_keV.begin(), energies_keV.end(),
                                     energy_keV);
    const std::size_t j = static_cast<std::size_t>(it - energies_keV.begin()) - 1;
    const double t = std::log(energy_keV / energies_keV[j])
        / std::log(energies_keV[j + 1] / energies_keV[j]);
    return (1.0 - t) * depth_cm[j] + t * depth_cm[j + 1];
}

} // namespace detail

struct DetectorModelPayload {
    std::uint32_t model_version = 1;
    double min_energy_keV = 10.0;
    double max_energy_keV = 12000.0;
    double min_cos_theta = -1.0;
    double max_cos_theta = 1.0;
    DetectorFamily family = DetectorFamily::NaI3x3;
    // The detector a fitted payload describes: a short label and its geometry
    // (the hand-set mechanistic pilots carry none). The rate-surface anchor,
    // the passive landmark and predict()'s acceptance derive from the
    // geometry; family, radius, length and density are set from it on load.
    std::string detector_name;
    std::optional<GeometryDescriptor> geometry;
    // Passive-layer landmark of the feature map's third coordinate,
    // 1/(1+(E/landmark)^3): passive_landmark_keV() of the geometry.
    double passive_landmark_keV = 55.0;
    FiniteDistanceDepth finite_distance;
    double radius_cm = 3.81;
    double length_cm = 7.62;
    double density_g_cm3 = 3.67;
    double back_incidence_scale = 0.08;
    double fep_base = 0.16;
    double fep_photoelectric_scale = 3.4;
    double fep_reinteraction_scale = 1.1;
    double single_escape_scale = 0.22;
    double double_escape_scale = 0.08;
    int resolution_half_window_channels = 160;
    std::array<double, 8> single_compton_shape{
        0.22, 0.60, 1.00, 1.22, 1.05, 0.63, 0.25, 0.05};
    std::array<double, 8> multiple_compton_shape{
        1.05, 1.24, 1.16, 0.91, 0.63, 0.38, 0.17, 0.04};
    // When populated, these positive process-resolved fields supersede the
    // original compact mechanistic pilot below. They interpolate functional
    // B-spline coefficients rather than stored response histograms.
    bool use_response_surrogate = false;
    DetectorResponseRbfComponent zero_deposit;
    DetectorResponseRbfComponent full_energy;
    // Hierarchical pair model. pair_fraction stores the odds of a detected
    // pulse having pair-production history. Pair outcomes are normalized only
    // within that parent population; full-energy atoms from pair and non-pair
    // histories remain separately fitted but coincide in observable energy.
    DetectorResponseRbfComponent pair_fraction;
    DetectorResponseRbfComponent pair_full_energy;
    DetectorResponseRbfComponent single_escape;
    DetectorResponseRbfComponent double_escape;
    DetectorResponseRbfComponent annihilation_511;
    // Crystal K-fluorescence escape lines: discrete atoms at
    // Ein - k_escape_alpha/beta_energy_keV, carved out of the smooth
    // photoelectric-escape continuum exactly like the pair escape peaks.
    // The energies are the probability-weighted K-line family means of the
    // crystal element (Ge: 9.833 / 10.943); the edge gates the atoms to
    // incident energies where K-shell absorption exists. Zeros disable the
    // atoms.
    double k_escape_alpha_energy_keV = 0.0;
    double k_escape_beta_energy_keV = 0.0;
    double k_escape_edge_keV = 0.0;
    // Upper end of each K component's fitted support, in the order
    // {alpha escape, beta escape, alpha absorbed, beta absorbed}. The RBF
    // yield fields are normalized Nadaraya-Watson averages, so above their
    // node cloud they level off at the nearest node's value instead of
    // decaying; without a gate the model emits a spurious satellite line
    // (~0.02-0.05% of the full-energy peak, and resolvable next to it) far
    // beyond the energies it was fitted at. The four supports genuinely
    // differ -- 320 / 175 / 280 / 150 keV for the fitted HPGe payload -- so one
    // shared limit would let the shorter-lived components extrapolate over
    // most of a decade. Zero disables the gate for a component.
    std::array<double, 4> k_escape_support_keV{};
    DetectorResponseRbfComponent k_alpha_escape;
    DetectorResponseRbfComponent k_beta_escape;
    // The complementary population: the primary photo-absorbs
    // in a dead layer and only its K x-ray reaches the active volume, where it
    // is absorbed, so the deposit is the LINE energy itself rather than
    // Ein - line. These two discrete lines carried 67% of the
    // photoelectric-escape continuum at 42.5 keV, which no nonnegative spline
    // over [0, Ein] can represent.
    DetectorResponseRbfComponent k_alpha_absorbed;
    DetectorResponseRbfComponent k_beta_absorbed;
    DetectorResponseRbfComponent detector_photoelectric_escape;
    DetectorResponseRbfComponent detector_single_compton;
    DetectorResponseRbfComponent detector_multiple_compton;
    DetectorResponseRbfComponent detector_pair_no_compton;
    DetectorResponseRbfComponent detector_pair_single_compton;
    DetectorResponseRbfComponent detector_pair_multiple_compton;
    // A small global functional surface replaces normalized RBF averaging for
    // the sharply varying total detection rate. Process-conditional rates and
    // all continuum shapes remain in the components above.
    DetectorRateSurface detected_rate_surface;
};

DetectorModelPayload nai3x3_payload();
DetectorModelPayload coaxial_hpge_payload();

/// Describe a fresh payload's detector: stores `geometry` and sets the crystal
/// radius, length and density, the rate-surface anchor's geometry and the
/// passive landmark from it -- what a fitted payload's reader does, except
/// that the reader keeps the stored landmark. Throws std::invalid_argument
/// for a box crystal or a descriptor that fails problems().
void set_detector_geometry(DetectorModelPayload& payload,
                           GeometryDescriptor geometry);

class DetectorPulseResponseModel {
public:
    explicit DetectorPulseResponseModel(
        DetectorModelPayload payload = nai3x3_payload())
        : payload_(std::move(payload)) {}

    template <typename Scalar>
    SpectralMeasure<Scalar> evaluate(
        double Ein_keV,
        const Scalar& cos_theta,
        const std::vector<double>& deposit_bins_keV) const;

    const DetectorModelPayload& payload() const { return payload_; }

    std::string serialize() const;
    static DetectorPulseResponseModel deserialize(const std::string& text);

private:
    PartialMassAttenuation mass_cross_sections(double energy_keV) const;

    DetectorModelPayload payload_;
};

struct EmissionLine {
    double energy_keV = 0.0;
    double photons_per_decay = 0.0;
};

/// Monoenergetic charged-particle source term. Continuous beta branches are
/// represented by a fixed positive quadrature of these lines; conversion
/// electrons remain discrete lines.
struct ChargedParticleLine {
    double kinetic_energy_keV = 0.0;
    double particles_per_decay = 0.0;
    bool is_positron = false;
};

template <typename Scalar>
struct ShieldParameters {
    Scalar areal_density_g_cm2 = Scalar(0);
    Scalar effective_Z = Scalar(13);
    /// Mass fraction of hydrogen, modeled as exact Z=1 alongside the
    /// effective-Z remainder. Hydrogen needs its own descriptor because no
    /// effective Z reproduces its electrons per gram (Compton strength) or
    /// its charged-particle stopping. Default 0 is the two-parameter shield.
    Scalar hydrogen_mass_fraction = Scalar(0);
};

template <typename Scalar>
struct ResolutionParameters {
    Scalar a0_keV2 = Scalar(0);
    Scalar a1_keV = Scalar(0);
    Scalar a2 = Scalar(0);
};

struct FullSpectrumPayload {
    std::uint32_t model_version = 1;
    std::string model_name = "FullSpecApprox-pilot";
    std::string validity = "10--12000 keV; elemental-equivalent spherical source";
    std::string training_hash = "untrained-mechanistic-pilot";
    std::string statistical_uncertainty = "not-estimated";
    std::vector<double> coefficient_standard_errors;
    std::string physics_toggles =
        "gamma-only;vacuum;no-fluorescence;no-beta-brems;no-pileup";
    SourceModelPayload source;
    DetectorModelPayload detector = nai3x3_payload();
    double distance_cm = 100.0;  ///< source to the detector geometry's reference point
    double incidence_cos_theta = 1.0;
    int source_quadrature_bins = 16;
};

/// Compose source escape and detector response and apply resolution last.
///
/// `activity` is the expected number of decays in the acquisition (Bq times
/// live time). The far-field geometric acceptance is derived from the detector
/// payload and distance. Air and room transport deliberately do not appear in
/// this v1 composition.
class FastSpectrumModel {
public:
    explicit FastSpectrumModel(FullSpectrumPayload payload = {})
        : payload_(std::move(payload)),
          source_(payload_.source),
          detector_(payload_.detector) {}

    FastSpectrumModel(
        FullSpectrumPayload payload,
        ChargedParticleSourceFluxModel charged_source)
        : payload_(std::move(payload)),
          source_(payload_.source),
          charged_source_(std::move(charged_source)),
          detector_(payload_.detector) {}

    template <typename Scalar>
    std::vector<Scalar> predict(
        const std::vector<EmissionLine>& lines,
        const Scalar& activity,
        const ShieldParameters<Scalar>& shield,
        const ResolutionParameters<Scalar>& resolution,
        const std::vector<double>& channels_keV) const;

    /// Compose nuclear photon lines, beta-spectrum quadrature nodes, and
    /// conversion-electron lines into one detector spectrum. All source terms
    /// are per parent decay and share the fitted activity/shield parameters.
    template <typename Scalar>
    std::vector<Scalar> predict(
        const std::vector<EmissionLine>& photon_lines,
        const std::vector<ChargedParticleLine>& charged_lines,
        const Scalar& activity,
        const ShieldParameters<Scalar>& shield,
        const ResolutionParameters<Scalar>& resolution,
        const std::vector<double>& channels_keV) const;

    const FullSpectrumPayload& payload() const { return payload_; }

    /// Compact, deterministic text payload. Raw campaign histograms are never
    /// serialized into this object.
    std::string serialize() const;
    static FastSpectrumModel deserialize(const std::string& text);

private:
    FullSpectrumPayload payload_;
    SourceGeomFluxModel source_;
    ChargedParticleSourceFluxModel charged_source_;
    DetectorPulseResponseModel detector_;
};

// ---- Template implementation ------------------------------------------------

template <typename Scalar>
std::array<Scalar, 4> SmoothElementCrossSection::evaluate(
    double energy_keV, const Scalar& effective_Z) const
{
    // Implemented as a positive soft cardinal interpolant. The narrow default
    // width reproduces integer-element cross sections while remaining smooth
    // throughout the fitted equivalent-Z range.
    std::array<Scalar, 4> numerator{
        Scalar(0), Scalar(0), Scalar(0), Scalar(0)};
    Scalar denominator = Scalar(0);
    const double inv_w2 = 1.0 / (blend_width_Z_ * blend_width_Z_);
    for (int Z = 1; Z <= 92; ++Z) {
        const Scalar dz = effective_Z - Scalar(Z);
        const Scalar weight = detail::scalar_exp(-dz * dz * inv_w2);
        const PartialMassAttenuation xs =
            element_mass_cross_sections(Z, energy_keV);
        numerator[0] += weight * xs.pe;
        numerator[1] += weight * xs.cs;
        numerator[2] += weight * xs.rs;
        numerator[3] += weight * xs.pp;
        denominator += weight;
    }
    const Scalar safe_denominator =
        denominator + Scalar(std::numeric_limits<double>::min());
    for (auto& value : numerator)
        value /= safe_denominator;
    return numerator;
}

template <typename Scalar>
Scalar SmoothElementCrossSection::compton_angular_density(
    double energy_keV,
    const Scalar& effective_Z,
    double cosine) const
{
    Scalar numerator = Scalar(0);
    Scalar denominator = Scalar(0);
    const double inv_w2 = 1.0 / (blend_width_Z_ * blend_width_Z_);
    for (int Z = 1; Z <= 92; ++Z) {
        const Scalar dz = effective_Z - Scalar(Z);
        const Scalar weight = detail::scalar_exp(-dz * dz * inv_w2);
        const double angular_density =
            ceelo::compton_angular_pdf_unnorm(energy_keV, Z, cosine)
            / ceelo::compton_angular_norm(energy_keV, Z);
        numerator += weight * angular_density;
        denominator += weight;
    }
    return numerator
        / (denominator + Scalar(std::numeric_limits<double>::min()));
}

template <typename Scalar>
Scalar SmoothElementCrossSection::compton_doppler_bin_probability(
    double energy_keV,
    const Scalar& effective_Z,
    double cosine,
    double lower_keV,
    double upper_keV) const
{
    Scalar numerator = Scalar(0);
    Scalar denominator = Scalar(0);
    const double inv_w2 = 1.0 / (blend_width_Z_ * blend_width_Z_);
    for (int Z = 1; Z <= 92; ++Z) {
        const Scalar dz = effective_Z - Scalar(Z);
        const Scalar weight = detail::scalar_exp(-dz * dz * inv_w2);
        numerator += weight * ceelo::doppler_energy_bin_probability(
            energy_keV, cosine, Z, lower_keV, upper_keV);
        denominator += weight;
    }
    return numerator
        / (denominator + Scalar(std::numeric_limits<double>::min()));
}

template <typename Scalar>
Scalar SmoothElementCrossSection::photonuclear_mass_attenuation(
    double energy_keV,
    const Scalar& effective_Z,
    const PhotonuclearSourcePayload& payload) const
{
    const PhotonuclearInterpolation interp =
        photonuclear_interpolation(payload, energy_keV);
    if (interp.below_table)
        return Scalar(0);
    Scalar numerator = Scalar(0);
    Scalar denominator = Scalar(0);
    const double inv_w2 = 1.0 / (blend_width_Z_ * blend_width_Z_);
    for (int Z = 1; Z <= 92; ++Z) {
        const Scalar dz = effective_Z - Scalar(Z);
        const Scalar weight = detail::scalar_exp(-dz * dz * inv_w2);
        const double lower_value =
            payload.cross_sections[interp.lower]
                .mass_attenuation_cm2_g[static_cast<std::size_t>(Z - 1)];
        const double upper_value =
            payload.cross_sections[interp.upper]
                .mass_attenuation_cm2_g[static_cast<std::size_t>(Z - 1)];
        const double value =
            lower_value
            + interp.fraction * (upper_value - lower_value);
        numerator += weight * value;
        denominator += weight;
    }
    return numerator
        / (denominator + Scalar(std::numeric_limits<double>::min()));
}

template <typename Scalar>
std::array<Scalar, 4> SmoothElementCrossSection::evaluate(
    double energy_keV,
    const Scalar& effective_Z,
    const Scalar& hydrogen_mass_fraction) const
{
    // Written as blend + f_H * (H - blend) so f_H = 0 reproduces the pure
    // effective-Z blend exactly.
    std::array<Scalar, 4> mixture = evaluate(energy_keV, effective_Z);
    const PartialMassAttenuation hydrogen =
        element_mass_cross_sections(1, energy_keV);
    const std::array<double, 4> hydrogen_values{
        hydrogen.pe, hydrogen.cs, hydrogen.rs, hydrogen.pp};
    for (std::size_t process = 0; process < mixture.size(); ++process) {
        mixture[process] += hydrogen_mass_fraction
            * (Scalar(hydrogen_values[process]) - mixture[process]);
    }
    return mixture;
}

template <typename Scalar>
Scalar SmoothElementCrossSection::hydrogen_compton_event_share(
    double energy_keV,
    const Scalar& effective_Z,
    const Scalar& hydrogen_mass_fraction) const
{
    const double hydrogen_cs =
        element_mass_cross_sections(1, energy_keV).cs;
    const Scalar blend_cs = evaluate(energy_keV, effective_Z)[1];
    const Scalar hydrogen_term =
        hydrogen_mass_fraction * Scalar(hydrogen_cs);
    return hydrogen_term
        / (hydrogen_term
           + (Scalar(1) - hydrogen_mass_fraction) * blend_cs
           + Scalar(std::numeric_limits<double>::min()));
}

template <typename Scalar>
Scalar SmoothElementCrossSection::compton_angular_density(
    double energy_keV,
    const Scalar& effective_Z,
    const Scalar& hydrogen_mass_fraction,
    double cosine) const
{
    const Scalar blend =
        compton_angular_density(energy_keV, effective_Z, cosine);
    const Scalar share = hydrogen_compton_event_share(
        energy_keV, effective_Z, hydrogen_mass_fraction);
    const double hydrogen_density =
        ceelo::compton_angular_pdf_unnorm(energy_keV, 1, cosine)
        / ceelo::compton_angular_norm(energy_keV, 1);
    return blend + share * (Scalar(hydrogen_density) - blend);
}

template <typename Scalar>
Scalar SmoothElementCrossSection::compton_doppler_bin_probability(
    double energy_keV,
    const Scalar& effective_Z,
    const Scalar& hydrogen_mass_fraction,
    double cosine,
    double lower_keV,
    double upper_keV) const
{
    const Scalar blend = compton_doppler_bin_probability(
        energy_keV, effective_Z, cosine, lower_keV, upper_keV);
    const Scalar share = hydrogen_compton_event_share(
        energy_keV, effective_Z, hydrogen_mass_fraction);
    const double hydrogen_probability =
        ceelo::doppler_energy_bin_probability(
            energy_keV, cosine, 1, lower_keV, upper_keV);
    return blend + share * (Scalar(hydrogen_probability) - blend);
}

template <typename Scalar>
Scalar SmoothElementCrossSection::photonuclear_mass_attenuation(
    double energy_keV,
    const Scalar& effective_Z,
    const Scalar& hydrogen_mass_fraction,
    const PhotonuclearSourcePayload& payload) const
{
    const Scalar blend =
        photonuclear_mass_attenuation(energy_keV, effective_Z, payload);
    const PhotonuclearInterpolation interp =
        photonuclear_interpolation(payload, energy_keV);
    if (interp.below_table)
        return blend;
    const double hydrogen_lower =
        payload.cross_sections[interp.lower].mass_attenuation_cm2_g[0];
    const double hydrogen_upper =
        payload.cross_sections[interp.upper].mass_attenuation_cm2_g[0];
    const double hydrogen_value =
        hydrogen_lower
        + interp.fraction * (hydrogen_upper - hydrogen_lower);
    return blend
        + hydrogen_mass_fraction * (Scalar(hydrogen_value) - blend);
}

namespace detail {

/// Eight-point Gauss-Legendre rule on [-1, 1], shared by every mechanistic
/// quadrature in this file, so the evaluator and the component-base helper use
/// the SAME literals -- `gauss_legendre_nodes()` regenerates the rule by
/// Newton iteration and would not reproduce the last ulp.
inline constexpr double kGaussLegendre8Nodes[8] = {
    -0.9602898564975362317, -0.7966664774136267396,
    -0.5255324099163289858, -0.1834346424956498049,
     0.1834346424956498049,  0.5255324099163289858,
     0.7966664774136267396,  0.9602898564975362317};
inline constexpr double kGaussLegendre8Weights[8] = {
    0.1012285362903762592, 0.2223810344533744705,
    0.3137066458778872873, 0.3626837833783619829,
    0.3626837833783619829, 0.3137066458778872873,
    0.2223810344533744705, 0.1012285362903762592};

/// Which mechanistic process creates the escaping secondary. Selects the
/// parent-attenuation and descendant-absorption axes of the 7-D "process"
/// feature vector.
enum class SourceSecondaryProcess {
    Pair511Atom,
    PairAnnihilationContinuum,
    PairAnnihilationInFlight,
    PhotoelectronBrems,
    ComptonElectronBrems,
    PairElectronBrems,
    FluorescenceKAlpha,
    FluorescenceKBeta,
    FluorescenceL,
};

/// Smooth, Scalar-valued inverse mass radiation length for the two-component
/// shield. Uses the same Gaussian effective-Z cardinal blend as the photon
/// cross sections, so it is differentiable in Z_eff, and the hydrogen mixture
/// is linear in X0^-1 (exact at f_H = 0).
template <typename Scalar>
Scalar smooth_inverse_radiation_length_cm2_g(
    const Scalar& effective_Z,
    const Scalar& hydrogen_mass_fraction,
    double z_blend_width)
{
    Scalar weight_sum = Scalar(0);
    Scalar inverse_sum = Scalar(0);
    const double inverse_blend_width2 =
        1.0 / (z_blend_width * z_blend_width);
    for (int Z = 1; Z <= 92; ++Z) {
        const Scalar delta_Z = effective_Z - Scalar(Z);
        const Scalar weight =
            scalar_exp(-delta_Z * delta_Z * inverse_blend_width2);
        weight_sum += weight;
        inverse_sum +=
            weight * element_inverse_radiation_length_cm2_g(Z);
    }
    Scalar inverse_x0 =
        inverse_sum
        / (weight_sum + Scalar(std::numeric_limits<double>::min()));
    inverse_x0 += hydrogen_mass_fraction
        * (Scalar(element_inverse_radiation_length_cm2_g(1)) - inverse_x0);
    return inverse_x0;
}

/// Seltzer-Berger production floor. Below it the emitted bremsstrahlung is
/// not tracked anywhere in this project, and in any shield thick enough to
/// matter a 10 keV photon escapes with probability ~0 regardless.
///
/// Not a tuning knob, and measured rather than asserted: moving it to 20 keV
/// scales the whole factor by 1.16 at tau = 0.2 and 1.18 at tau = 8, i.e. a
/// near-constant the fitted conditional yield absorbs exactly. Span and slope
/// are unchanged to three decimals.
constexpr double kBremsstrahlungFloorKeV = 10.0;

/// Escape probability of a broad-spectrum descendant, averaged over its own
/// emission spectrum rather than evaluated at one proxy energy:
///
///     escape(tau, [k_lo, k_hi])
///       = integral_{ln k_lo}^{ln k_hi} escape(tau, kappa(k)) d ln k
///         / (ln k_hi - ln k_lo)
///
/// Bremsstrahlung is emitted with a ~1/k number spectrum, which is FLAT in
/// ln k -- that is the whole reason for the variable. The emission weight is
/// then uniform and an unweighted Gauss-Legendre mean IS the quadrature:
/// GL weights sum to 2 on [-1, 1], so 0.5 * w_i is already the normalized
/// measure and the interval half-width cancels.
///
/// GL8 suffices: GL16 gives the same span and slope to three decimals over the
/// whole bank. The base is divided out at fit time and multiplied back in at
/// evaluation, so an identical rule on both sides makes any quadrature error a
/// smooth change of base rather than an accuracy loss.
///
/// Removal is PE + Compton + pair. A Compton scatter both lengthens the exit
/// chord and drops the photon toward the PE-dominated regime, so for a base
/// whose job is to carry the tau falloff it is removal -- which is already
/// what `pair_511_atom` does for its own child. Measured over 1122 bank cells,
/// the two changes are strictly complementary: median conditional-yield span
/// for photoelectron / compton-electron / pair-electron brems goes
/// 76.7 / 113.7 / 340.8 with a creation-only base, to 7.9 / 16.6 / 26.0 with a
/// single proxy-energy escape factor (PE+pair removal), to 5.8 / 5.5 / 29.9
/// spectrum-averaged with PE+pair removal -- and to 1.4 / 4.5 / 2.6 with both.
/// Neither change is worth much alone.
template <typename Scalar>
inline Scalar spectrum_averaged_escape_factor(
    const Scalar& tau_total,
    const Scalar& areal_density_g_cm2,
    const Scalar& effective_Z,
    const Scalar& hydrogen_mass_fraction,
    double child_min_keV,
    double child_max_keV,
    const SmoothElementCrossSection& cross_sections)
{
    // Both bounds are plain doubles -- functions of E0 alone -- so this branch
    // can never split a ceres::Jet derivative path. It fires for
    // pair_electron_brems at every E0 below 1022 keV, where the parent kinetic
    // energy is NEGATIVE, and for compton_electron_brems below ~54 keV, where
    // the recoil maximum falls under the floor. Returning 1 is exact: the
    // caller's turn-on is identically 0 for any parent under 200 keV, so the
    // base is 0 * 1 either way -- but 0 * NaN is NaN.
    if (!(child_max_keV > child_min_keV))
        return Scalar(1);
    const double log_low = std::log(child_min_keV);
    const double log_high = std::log(child_max_keV);
    const double log_middle = 0.5 * (log_high + log_low);
    const double log_half = 0.5 * (log_high - log_low);
    const Scalar safe_tau =
        tau_total + Scalar(std::numeric_limits<double>::min());
    Scalar average = Scalar(0);
    for (int node = 0; node < 8; ++node) {
        const double child_keV = std::exp(
            log_middle + log_half * kGaussLegendre8Nodes[node]);
        const std::array<Scalar, 4> mu_child = cross_sections.evaluate(
            child_keV, effective_Z, hydrogen_mass_fraction);
        const Scalar child_tau = areal_density_g_cm2
            * (mu_child[0] + mu_child[1] + mu_child[3]);
        average += Scalar(0.5 * kGaussLegendre8Weights[node])
            * secondary_escape_factor(tau_total, child_tau / safe_tau);
    }
    return average;
}

/// Index into `SourceComponentBases::base`, in the same order the surrogate
/// file, the feature helper and the Python fitter's ALL_COMPONENTS use.
enum SourceComponentIndex : std::size_t {
    kPair511Atom = 0,
    kPairAnnihilationContinuum,
    kPairAnnihilationInFlight,
    kPhotoelectronBrems,
    kComptonElectronBrems,
    kPairElectronBrems,
    kFluorescenceKAlpha,
    kFluorescenceKBeta,
    kFluorescenceL,
    kSourceComponentCount,
};
static_assert(kSourceComponentCount == SourceSecondarySurrogatePayload::kComponentCount,
              "SourceComponentIndex and SourceSecondarySurrogatePayload::components() disagree");

template <typename Scalar>
struct SourceComponentBases {
    /// Probability of a first pair conversion, shared by the three
    /// annihilation channels.
    Scalar pair_first_probability = Scalar(0);
    std::array<Scalar, kSourceComponentCount> base{};
    /// Per fluorescence family, the lines (energy, base) whose sum is that
    /// family's base: the evaluator spreads the family's fitted yield over
    /// them in proportion.
    std::array<std::vector<std::pair<double, Scalar>>, kFluorescenceFamilies>
        fluorescence_lines;
};

/// Vacancy sources the fluorescence bases include (a bit set; the helper's
/// --fluorescence-base-terms picks them for span studies): the primary's own
/// photoabsorption, a once-Compton-scattered primary photoabsorbed on its way
/// out, an absorbed annihilation photon, electron-impact ionization by the
/// primary's photo-, Compton and pair electrons, and those electrons'
/// bremsstrahlung photoabsorbed in the shield.
constexpr int kFluorescenceDirect = 1;
constexpr int kFluorescenceCompton = 2;
constexpr int kFluorescencePair = 4;
constexpr int kFluorescenceImpact = 8;
constexpr int kFluorescenceBrems = 16;
constexpr int kFluorescenceAllTerms = 31;

/// The effective-Z blend's elements that carry fluorescence, with their
/// normalized blend weights tapered smoothly (C1) to zero below 0.5% and
/// renormalized, so an integer effective Z emits exactly its own lines and a
/// fractional one its two neighbours'.
template <typename Scalar>
std::vector<std::pair<int, Scalar>> fluorescence_element_weights(
    const Scalar& effective_Z, double blend_width_Z)
{
    std::array<Scalar, 93> raw{};
    Scalar total = Scalar(0);
    const double inv_w2 = 1.0 / (blend_width_Z * blend_width_Z);
    for (int Z = 1; Z <= 92; ++Z) {
        const Scalar dz = effective_Z - Scalar(Z);
        raw[Z] = scalar_exp(-dz * dz * inv_w2);
        total += raw[Z];
    }
    std::vector<std::pair<int, Scalar>> out;
    Scalar kept = Scalar(0);
    for (int Z = 1; Z <= 92; ++Z) {
        const Scalar w = raw[Z] / (total + Scalar(std::numeric_limits<double>::min()));
        if (!(w > Scalar(0.005))) continue;
        Scalar taper = Scalar(1);
        if (w < Scalar(0.02)) {
            const Scalar t = (w - Scalar(0.005)) / Scalar(0.015);
            taper = t * t * (Scalar(3) - Scalar(2) * t);
        }
        out.emplace_back(Z, w * taper);
        kept += w * taper;
    }
    for (auto& element : out) element.second /= kept;
    return out;
}

/// Yield-weighted mean removal coefficient (PE + Compton + pair, cm^2/g) of a
/// fluorescence family's lines in the shield: the family's own absorption
/// depth per g/cm^2, for the process features.
template <typename Scalar>
Scalar fluorescence_family_removal(
    int family, const Scalar& effective_Z, const Scalar& hydrogen_mass_fraction,
    const SmoothElementCrossSection& cross_sections)
{
    Scalar numerator = Scalar(0), denominator = Scalar(0);
    for (const auto& [Z, weight] :
         fluorescence_element_weights(effective_Z, cross_sections.blend_width_Z())) {
        for (const auto& line : element_fluorescence(Z).lines) {
            if (line.family != family) continue;
            const double y = line.per_k + line.per_l[0] + line.per_l[1] + line.per_l[2];
            const auto mu = cross_sections.evaluate(
                line.energy_keV, effective_Z, hydrogen_mass_fraction);
            numerator += weight * Scalar(y) * (mu[0] + mu[1] + mu[3]);
            denominator += weight * Scalar(y);
        }
    }
    return denominator > Scalar(0) ? numerator / denominator : Scalar(0);
}

/// The six mechanistic component baselines.
///
/// This is the SINGLE implementation. `SourceGeomFluxModel::evaluate`
/// multiplies these back into the fitted conditional yield, and
/// `ceelo_source_cross_section_features --emit-component-bases` hands the same
/// numbers to the Python fitter, which divides them OUT at fit time. A second
/// copy of this arithmetic would be a purely multiplicative bias that no
/// metric in the fit report can see. Same reasoning as
/// `source_process_features`, one layer deeper.
///
/// `mu0_total` is supplied by the caller rather than recomputed, exactly as
/// `source_process_features` takes `tau_total`: it is what carries the
/// `enable_rayleigh` gate and the photonuclear term, neither of which is
/// recoverable from `mu0` alone.
template <typename Scalar>
SourceComponentBases<Scalar> source_component_bases(
    double E0_keV,
    const Scalar& areal_density_g_cm2,
    const Scalar& effective_Z,
    const Scalar& hydrogen_mass_fraction,
    const std::array<Scalar, 4>& mu0,
    const Scalar& mu0_total,
    const SmoothElementCrossSection& cross_sections,
    const SourceSecondarySurrogatePayload& options,
    int fluorescence_terms = kFluorescenceAllTerms)
{
    constexpr double electron_mass_keV = 510.99895;
    SourceComponentBases<Scalar> result;

    const auto mu511 = cross_sections.evaluate(
        electron_mass_keV, effective_Z, hydrogen_mass_fraction);
    const Scalar mu511_removal = mu511[0] + mu511[1] + mu511[3];
    const Scalar tau_total = areal_density_g_cm2 * mu0_total;

    // Escape probability of a secondary born at a depth drawn against the
    // primary's own attenuation. kappa = (child absorption depth) / (primary
    // total depth); both are already needed for the process features, so this
    // introduces no new physics input.
    const auto secondary_escape = [&](const Scalar& child_tau) {
        const Scalar kappa =
            child_tau
            / (tau_total + Scalar(std::numeric_limits<double>::min()));
        return secondary_escape_factor(tau_total, kappa);
    };

    // Same, for a descendant that is a broad spectrum rather than a line:
    // the closed form averaged over the descendant's own emission spectrum
    // in log k, with removal = PE + Compton + pair.
    const auto spectrum_escape = [&](double child_min_keV,
                                     double child_max_keV) {
        return spectrum_averaged_escape_factor(
            tau_total, areal_density_g_cm2, effective_Z,
            hydrogen_mass_fraction, child_min_keV, child_max_keV,
            cross_sections);
    };

    // Mechanistic first-pair baseline. The pair vertex lies at radius x
    // along the primary ray; each annihilation photon is marginally
    // isotropic, so average its no-energy-changing escape probability
    // over the remaining spherical chord. Rayleigh is treated as elastic
    // in this baseline. Its path-length correction, pair creation after
    // prior Compton/Rayleigh interactions, and other higher generations
    // are learned by the positive conditional residual.
    Scalar pair_atom_base = Scalar(0);
    for (int depth_index = 0; depth_index < 8; ++depth_index) {
        const double depth_fraction =
            0.5 * (kGaussLegendre8Nodes[depth_index] + 1.0);
        const Scalar pair_density =
            areal_density_g_cm2 * mu0[3]
            * scalar_exp(
                -areal_density_g_cm2 * mu0_total
                * depth_fraction);
        result.pair_first_probability +=
            Scalar(0.5 * kGaussLegendre8Weights[depth_index])
            * pair_density;
        Scalar child_escape_average = Scalar(0);
        for (int angle_index = 0; angle_index < 8; ++angle_index) {
            const double cosine = kGaussLegendre8Nodes[angle_index];
            const double transverse2 =
                depth_fraction * depth_fraction
                * (1.0 - cosine * cosine);
            const double exit_path_fraction =
                -depth_fraction * cosine
                + std::sqrt(std::max(0.0, 1.0 - transverse2));
            child_escape_average +=
                Scalar(0.5 * kGaussLegendre8Weights[angle_index])
                * scalar_exp(
                    -areal_density_g_cm2 * mu511_removal
                    * exit_path_fraction);
        }
        pair_atom_base +=
            Scalar(kGaussLegendre8Weights[depth_index])
            * pair_density * child_escape_average;
    }
    result.base[kPair511Atom] = pair_atom_base;

    // The at-rest annihilation continuum is a 511 keV photon that HAS
    // Compton-scattered, so its escape is the same closed form as
    // pair_511_atom's child.
    result.base[kPairAnnihilationContinuum] =
        Scalar(2)
        * result.pair_first_probability
        * (Scalar(1) - scalar_exp(
            -areal_density_g_cm2 * mu511[1]))
        * secondary_escape(
            areal_density_g_cm2 * mu511_removal);

    const double in_flight_turn_on =
        1.0 - std::exp(
            -std::max(
                0.0, E0_keV - 2.0 * electron_mass_keV)
            / 2000.0);
    // In-flight annihilation emits a forward-boosted PAIR spanning roughly
    // m_e to E0 + m_e, not a 511 keV line -- which is why the evaluator
    // appends this component's density at support E0 rather than at 511. The
    // upper bound is taken as E0: the extra m_e moves the slope by 0.003 and
    // would need its own justification. Span 44.9x creation-only, 8.9x at the
    // 511 line, 4.4x here.
    result.base[kPairAnnihilationInFlight] =
        Scalar(2.0 * in_flight_turn_on)
        * result.pair_first_probability
        * spectrum_escape(electron_mass_keV, E0_keV);

    const double backscatter_energy =
        E0_keV
        / (1.0 + 2.0 * E0_keV / electron_mass_keV);
    const auto electron_turn_on = [](double maximum_keV) {
        return 1.0 - std::exp(
            -std::max(0.0, maximum_keV - 200.0) / 650.0);
    };
    // Parent electron kinetic energies, one per creating interaction: the
    // photoelectron takes the photon less the binding of the shell it left
    // (per element of the blend, as the walk samples it), the Compton recoil
    // takes at most the Compton edge, the pair electron shares E0 - 2 m_e.
    // Each is the top of that channel's bremsstrahlung spectrum.
    const std::array<double, 3> parent_keV{{
        E0_keV,
        E0_keV - backscatter_energy,
        E0_keV - 2.0 * electron_mass_keV,
    }};
    std::array<Scalar, 3> brems_escape{
        {Scalar(1), Scalar(1), Scalar(1)}};
    for (std::size_t index = 1; index < 3; ++index) {
        // Skip the quadrature where the turn-on is identically zero: the
        // base is 0 either way, and this removes all 24 cross-section
        // evaluations below 200 keV.
        if (electron_turn_on(parent_keV[index]) > 0.0) {
            brems_escape[index] = spectrum_escape(
                kBremsstrahlungFloorKeV, parent_keV[index]);
        }
    }
    // The photoelectron's turn-on x escape, averaged over its energies: a
    // K-shell absorption in lead at 250 keV leaves 162 keV, under the walk's
    // 200 keV floor.
    Scalar photoelectron_turn_on_escape = Scalar(0);
    for (const auto& [Z, weight] :
         fluorescence_element_weights(effective_Z, cross_sections.blend_width_Z()))
        for (const auto& [energy_keV, probability] : element_photoelectron_energies(Z, E0_keV)) {
            const double turn_on = electron_turn_on(energy_keV);
            if (turn_on > 0.0 && probability > 0.0)
                photoelectron_turn_on_escape += weight * Scalar(probability * turn_on)
                    * spectrum_escape(kBremsstrahlungFloorKeV, energy_keV);
        }
    const Scalar tau_pe = areal_density_g_cm2 * mu0[0];
    const Scalar tau_cs = areal_density_g_cm2 * mu0[1];
    const Scalar tau_pp = areal_density_g_cm2 * mu0[3];
    result.base[kPhotoelectronBrems] =
        photoelectron_turn_on_escape
            * (Scalar(1) - scalar_exp(-tau_pe));
    result.base[kComptonElectronBrems] =
        electron_turn_on(parent_keV[1])
            * (Scalar(1) - scalar_exp(-tau_cs))
            * brems_escape[1];
    result.base[kPairElectronBrems] =
        electron_turn_on(parent_keV[2])
            * (Scalar(1) - scalar_exp(-tau_pp))
            * brems_escape[2];

    // Fluorescence atoms: shield K and L x-rays that escape unscattered. A
    // line's source per g/cm^2 at photon energy E is S(E) = sum over the
    // blend's elements of weight x mu_pe x (share_K x yield per K vacancy +
    // sum_s share_Ls x yield per L_s vacancy). Vacancy sources:
    //   D  the primary absorbed at depth x;
    //   C  the primary Compton-scattered once at x (angle c to the outward
    //      radius, which is the primary's direction) and absorbed on its way
    //      out -- the fraction S(E1)/removal(E1) of its removal -- with the
    //      x-ray leaving from the Compton vertex;
    //   P  an annihilation photon absorbed, the same way at 511 keV.
    // Each line escapes from depth x as the spherical-chord average with
    // removal = PE + Compton + pair at the line (a Rayleigh-only x-ray stays
    // in the atom). One photon per vacancy: weights w_i / 2.
    if (fluorescence_terms != 0) {
        const auto elements =
            fluorescence_element_weights(effective_Z, cross_sections.blend_width_Z());
        const Scalar non_h = Scalar(1) - hydrogen_mass_fraction;
        // Quadrature: an opaque line (Pb L x-rays at 50 g/cm^2: optical depth
        // ~6000) escapes only from a skin 1e-4 of the radius thick, so depth
        // uses Gauss-Legendre-8 segments refined geometrically toward the
        // surface, and the escape angle is split at its kink (c = 0, grazing).
        // Angle weights average over c in [-1, 1]; depth weights integrate
        // over [0, 1].
        constexpr std::array<double, 7> kSegment{
            0.0, 0.9, 0.99, 0.999, 0.9999, 0.99999, 1.0};
        constexpr int kDepth = 8 * 6, kAngle = 16;
        std::array<double, kDepth> depth{}, depth_weight{};
        for (int segment = 0; segment < 6; ++segment) {
            const double a = kSegment[segment], b = kSegment[segment + 1];
            for (int n = 0; n < 8; ++n) {
                depth[8 * segment + n] = a + 0.5 * (b - a) * (kGaussLegendre8Nodes[n] + 1.0);
                depth_weight[8 * segment + n] = 0.5 * (b - a) * kGaussLegendre8Weights[n];
            }
        }
        std::array<double, kAngle> cosine{}, angle_weight{};
        for (int n = 0; n < 8; ++n) {
            const double half = 0.5 * (kGaussLegendre8Nodes[n] + 1.0);
            cosine[n] = -half;
            cosine[8 + n] = half;
            angle_weight[n] = angle_weight[8 + n] = 0.25 * kGaussLegendre8Weights[n];
        }
        std::array<std::array<double, kAngle>, kDepth> chord{};
        std::array<Scalar, kDepth> primary{};
        for (int i = 0; i < kDepth; ++i) {
            for (int j = 0; j < kAngle; ++j) {
                const double c = cosine[j];
                chord[i][j] = -depth[i] * c
                    + std::sqrt(std::max(0.0, 1.0 - depth[i] * depth[i] * (1.0 - c * c)));
            }
            primary[i] = areal_density_g_cm2
                * scalar_exp(-areal_density_g_cm2 * mu0_total * depth[i]);
        }
        const auto escape_from = [&](int i, const Scalar& removal) {
            Scalar escape = Scalar(0);
            for (int j = 0; j < kAngle; ++j)
                escape += Scalar(angle_weight[j])
                    * scalar_exp(-areal_density_g_cm2 * removal * chord[i][j]);
            return escape;
        };
        // Once-scattered photons: energy, angular density (normalized over c)
        // and removal per angle node; the probability they are removed before
        // leaving, per depth and angle.
        const bool compton_term = (fluorescence_terms & kFluorescenceCompton) != 0;
        std::array<double, kAngle> e1{};
        std::array<Scalar, kAngle> angular{}, removal1{};
        std::array<std::array<Scalar, kAngle>, kDepth> removed1{};
        if (compton_term) {
            for (int j = 0; j < kAngle; ++j) {
                e1[j] = E0_keV / (1.0 + E0_keV / electron_mass_keV * (1.0 - cosine[j]));
                angular[j] = cross_sections.compton_angular_density(
                    E0_keV, effective_Z, hydrogen_mass_fraction, cosine[j]);
                const auto mu1 = cross_sections.evaluate(
                    e1[j], effective_Z, hydrogen_mass_fraction);
                removal1[j] = mu1[0] + mu1[1] + mu1[3];
                for (int i = 0; i < kDepth; ++i)
                    removed1[i][j] = Scalar(1)
                        - scalar_exp(-areal_density_g_cm2 * removal1[j] * chord[i][j]);
            }
        }
        // Annihilation photons: the probability one is removed before leaving.
        const bool pair_term = (fluorescence_terms & kFluorescencePair) != 0
            && E0_keV > 2.0 * electron_mass_keV;
        // Electrons of the primary's interactions: the photoelectron (~E0),
        // the mean Compton recoil and each pair lepton; their impact-ionization
        // x-rays leave from the interaction depth.
        const bool impact_term = (fluorescence_terms & kFluorescenceImpact) != 0;
        const bool brems_term = (fluorescence_terms & kFluorescenceBrems) != 0;
        double mean_recoil_keV = 0.0;
        if (impact_term || brems_term) {
            // Klein-Nishina x S of iron: the bound-electron correction moves the
            // mean recoil negligibly at the energies impact ionization matters.
            double norm = 0.0;
            for (int j = 0; j < kAngle; ++j) {
                const double c = cosine[j];
                const double density = ceelo::compton_angular_pdf_unnorm(E0_keV, 26, c);
                mean_recoil_keV += angle_weight[j] * density
                    * (E0_keV - E0_keV / (1.0 + E0_keV / electron_mass_keV * (1.0 - c)));
                norm += angle_weight[j] * density;
            }
            mean_recoil_keV = norm > 0.0 ? mean_recoil_keV / norm : 0.0;
        }
        const double pair_lepton_keV =
            std::max(0.0, 0.5 * (E0_keV - 2.0 * electron_mass_keV));
        std::array<Scalar, kDepth> removed511{};
        if (pair_term)
            for (int i = 0; i < kDepth; ++i)
                removed511[i] = Scalar(1) - escape_from(i, mu511_removal);
        for (const auto& [Z, weight] : elements) {
            const ElementFluorescence& element = element_fluorescence(Z);
            if (element.lines.empty()) continue;
            // Vacancy production per g/cm^2 of the shield, {K, L1, L2, L3},
            // by photoabsorption.
            const auto vacancies = [&](double energy_keV) {
                const double pe = element_mass_cross_sections(Z, energy_keV).pe;
                const auto share = element_photoelectric_vacancy_shares(Z, energy_keV);
                std::array<Scalar, 4> out{};
                for (int k = 0; k < 4; ++k) out[k] = non_h * weight * Scalar(pe * share[k]);
                return out;
            };
            auto v0 = vacancies(E0_keV);
            // The primary's own Compton collisions leave a vacancy in a shell in
            // proportion to its occupancy (Doppler sampling): K 2, L1 2, L2 2,
            // L3 4 of Z electrons, where the shell can be ionized.
            {
                const double cs = element_mass_cross_sections(Z, E0_keV).cs;
                const std::array<double, 4> occupancy{2.0, 2.0, 2.0, 4.0};
                for (int k = 0; k < 4; ++k)
                    if (Z > (k == 0 ? 2 : 4) && (k > 0 || E0_keV > element.k_edge_keV))
                        v0[k] += non_h * weight * Scalar(cs * occupancy[k] / Z);
            }
            std::array<std::array<Scalar, 4>, kAngle> v1{};
            if (compton_term)
                for (int j = 0; j < kAngle; ++j) v1[j] = vacancies(e1[j]);
            std::array<Scalar, 4> v511{};
            if (pair_term) v511 = vacancies(electron_mass_keV);
            // K vacancies per g/cm^2 from impact ionization (yields count K
            // x-rays, so divide omega_K back out).
            Scalar impact_k = Scalar(0);
            if (impact_term && element.omega_K > 0.0) {
                double y_pe = 0.0;
                for (const auto& [energy_keV, probability] :
                     element_photoelectron_energies(Z, E0_keV))
                    if (probability > 0.0)
                        y_pe += probability * element_impact_xray_yield(Z, energy_keV);
                const double y_cs = element_impact_xray_yield(Z, mean_recoil_keV);
                const double y_pp = 2.0 * element_impact_xray_yield(Z, pair_lepton_keV);
                impact_k = non_h * weight
                    * (mu0[0] * Scalar(y_pe) + mu0[1] * Scalar(y_cs) + mu0[3] * Scalar(y_pp))
                    / Scalar(element.omega_K);
            }
            // Bremsstrahlung of the primary's electrons absorbed in the shield:
            // a thick-target (Kramers) photon spectrum (R/X0)(1 - k/T)/k above
            // the element's K edge (or the 10 keV floor), with the brems bases'
            // turn-on (the Monte Carlo walks electrons above 200 keV only);
            // each photon absorbed before it leaves the depth it was born at
            // makes vacancies like the primary's photoabsorption.
            std::array<std::array<Scalar, 4>, kDepth> brems_vacancies{};
            if (brems_term) {
                const double inverse_x0 = element_inverse_radiation_length_cm2_g(Z);
                // {energy, parent process, probability}: the photoelectron at
                // each of its energies, the mean Compton recoil, a pair lepton.
                std::vector<std::tuple<double, int, double>> electrons;
                for (const auto& [energy_keV, probability] :
                     element_photoelectron_energies(Z, E0_keV))
                    if (probability > 0.0) electrons.emplace_back(energy_keV, 0, probability);
                electrons.emplace_back(mean_recoil_keV, 1, 1.0);
                electrons.emplace_back(pair_lepton_keV, 2, 1.0);
                for (const auto& [T, source_index, probability] : electrons) {
                    const double turn_on = 1.0 - std::exp(-std::max(0.0, T - 200.0) / 650.0);
                    const double k_lo = std::max(kFluorescenceMinimumKeV, element.k_edge_keV);
                    if (!(turn_on > 0.0) || !(T > k_lo)) continue;
                    const double photons_per_lnk = probability * turn_on
                        * element_csda_range_g_cm2(Z, T, false) * inverse_x0;
                    const Scalar parent = source_index == 0 ? mu0[0]
                        : (source_index == 1 ? mu0[1] : Scalar(2) * mu0[3]);
                    const double span = std::log(T / k_lo);
                    for (int q = 0; q < 8; ++q) {
                        const double k = k_lo * std::exp(0.5 * span * (kGaussLegendre8Nodes[q] + 1.0));
                        const double weight_q = 0.5 * span * kGaussLegendre8Weights[q]
                            * photons_per_lnk * (1.0 - k / T);
                        const auto mu_k = cross_sections.evaluate(k, effective_Z, hydrogen_mass_fraction);
                        const Scalar removal_k = mu_k[0] + mu_k[1] + mu_k[3];
                        const auto v_k = vacancies(k);
                        for (int i = 0; i < kDepth; ++i) {
                            const Scalar absorbed = (Scalar(1) - escape_from(i, removal_k))
                                / (removal_k + Scalar(std::numeric_limits<double>::min()));
                            for (int m = 0; m < 4; ++m)
                                brems_vacancies[i][m] += parent * Scalar(weight_q) * v_k[m] * absorbed;
                        }
                    }
                }
            }
            // Vacancy rates {K, L1, L2, L3} per g/cm^2 of primary path at each
            // depth, every source together; the lines then differ only by
            // their yields and their own escape.
            std::array<std::array<Scalar, 4>, kDepth> rate{};
            for (int i = 0; i < kDepth; ++i) {
                for (int m = 0; m < 4; ++m) {
                    Scalar r = Scalar(0);
                    if (fluorescence_terms & kFluorescenceDirect) r += v0[m];
                    if (compton_term) {
                        Scalar scattered = Scalar(0);
                        for (int j = 0; j < kAngle; ++j)
                            scattered += Scalar(2.0 * angle_weight[j]) * angular[j]
                                * removed1[i][j] * v1[j][m]
                                / (removal1[j] + Scalar(std::numeric_limits<double>::min()));
                        r += mu0[1] * scattered;
                    }
                    if (pair_term)
                        r += mu0[3] * removed511[i] * Scalar(2) * v511[m]
                            / (mu511_removal + Scalar(std::numeric_limits<double>::min()));
                    if (brems_term) r += brems_vacancies[i][m];
                    rate[i][m] = r;
                }
                rate[i][0] += impact_k;
            }
            // L vacancies from the shield's own K x-rays absorbed where they are
            // born (most are, in a high-Z shield): per depth, per K vacancy.
            std::array<std::array<Scalar, 3>, kDepth> k_cascade{};
            bool any_l_line = false;
            for (const auto& line : element.lines) any_l_line |= line.family == 2;
            if (any_l_line) {
                for (const auto& k_line : element.lines) {
                    if (k_line.family == 2 || !(k_line.per_k > 0.0)) continue;
                    const auto mu_k = cross_sections.evaluate(
                        k_line.energy_keV, effective_Z, hydrogen_mass_fraction);
                    const Scalar removal_k = mu_k[0] + mu_k[1] + mu_k[3];
                    const auto v_k = vacancies(k_line.energy_keV);
                    for (int i = 0; i < kDepth; ++i) {
                        const Scalar absorbed = Scalar(k_line.per_k)
                            * (Scalar(1) - escape_from(i, removal_k))
                            / (removal_k + Scalar(std::numeric_limits<double>::min()));
                        for (int s_l = 0; s_l < 3; ++s_l)
                            k_cascade[i][s_l] += absorbed * v_k[1 + s_l];
                    }
                }
            }
            for (const FluorescenceLineYield& line : element.lines) {
                const auto mu_line = cross_sections.evaluate(
                    line.energy_keV, effective_Z, hydrogen_mass_fraction);
                const Scalar removal_line = mu_line[0] + mu_line[1] + mu_line[3];
                Scalar base = Scalar(0);
                for (int i = 0; i < kDepth; ++i) {
                    Scalar line_rate = rate[i][0] * Scalar(line.per_k)
                        + rate[i][1] * Scalar(line.per_l[0])
                        + rate[i][2] * Scalar(line.per_l[1])
                        + rate[i][3] * Scalar(line.per_l[2]);
                    if (line.family == 2)
                        line_rate += rate[i][0]
                            * (k_cascade[i][0] * Scalar(line.per_l[0])
                               + k_cascade[i][1] * Scalar(line.per_l[1])
                               + k_cascade[i][2] * Scalar(line.per_l[2]));
                    base += Scalar(depth_weight[i]) * primary[i] * line_rate
                        * escape_from(i, removal_line);
                }
                result.fluorescence_lines[line.family].emplace_back(line.energy_keV, base);
                result.base[kFluorescenceKAlpha + line.family] += base;
            }
        }
    }
    return result;
}

/// CSDA range (g/cm^2) of a charged particle in the two-component shield: the
/// effective-Z blend of element ranges, combined with hydrogen reciprocally
/// (stopping powers add by mass fraction). The charged-particle evaluator and
/// its feature helper both take the range from here.
template <typename Scalar>
Scalar charged_blend_csda_range(
    double kinetic_energy_keV, const Scalar& effective_Z,
    const Scalar& hydrogen_mass_fraction, bool is_positron, double z_blend_width)
{
    Scalar weight_sum = Scalar(0), range_sum = Scalar(0);
    const double inverse_blend_width2 = 1.0 / (z_blend_width * z_blend_width);
    for (int Z = 1; Z <= 92; ++Z) {
        const Scalar delta_Z = effective_Z - Scalar(Z);
        const Scalar weight = scalar_exp(-delta_Z * delta_Z * inverse_blend_width2);
        weight_sum += weight;
        range_sum += weight * element_csda_range_g_cm2(Z, kinetic_energy_keV, is_positron);
    }
    const Scalar range = range_sum / (weight_sum + Scalar(std::numeric_limits<double>::min()));
    const double hydrogen_range = element_csda_range_g_cm2(1, kinetic_energy_keV, is_positron);
    return range
        / (Scalar(1) + hydrogen_mass_fraction * (range / Scalar(hydrogen_range) - Scalar(1)));
}

/// Shield fluorescence of a charged particle stopping in the shield, per
/// family and per line: `a`, from electron-impact ionization along its path and
/// its bremsstrahlung absorbed in the shield (both particles); `b`, from the two
/// annihilation photons absorbed, per unit positron containment. A line's
/// base is a + contained x b. The particle's path is taken uniformly along the
/// radius out to the stop radius (the smooth min(1, R/t) the annihilation atom
/// uses). Impact-ionization vacancies sit on that path. A photon born there
/// (bremsstrahlung, annihilation) is followed along its own ray, and its line
/// escapes from where it is absorbed -- in a thick shield, near the surface,
/// far from where it was born. The shield's own K x-rays absorbed where they
/// are born feed its L lines.
template <typename Scalar>
struct ChargedFluorescenceBases {
    std::array<Scalar, kFluorescenceFamilies> a{}, b{};
    struct Line {
        double energy_keV = 0.0;
        int family = 0;
        Scalar a = Scalar(0), b = Scalar(0);
    };
    std::vector<Line> lines;
};

template <typename Scalar>
ChargedFluorescenceBases<Scalar> charged_fluorescence_bases(
    double kinetic_energy_keV,
    const Scalar& areal_density_g_cm2,
    const Scalar& effective_Z,
    const Scalar& hydrogen_mass_fraction,
    bool is_positron,
    const SmoothElementCrossSection& cross_sections)
{
    constexpr double electron_mass_keV = 510.99895;
    ChargedFluorescenceBases<Scalar> result;
    const double width = cross_sections.blend_width_Z();
    const Scalar csda_range = charged_blend_csda_range(
        kinetic_energy_keV, effective_Z, hydrogen_mass_fraction, is_positron, width);
    const Scalar range_ratio =
        areal_density_g_cm2 / (csda_range + Scalar(std::numeric_limits<double>::min()));
    const Scalar path_fraction = Scalar(1) - scalar_exp(-range_ratio);
    const Scalar r2 = range_ratio * range_ratio, r4 = r2 * r2, r8 = r4 * r4;
    const Scalar stop_radius = scalar_exp(Scalar(-0.125) * scalar_log(Scalar(1) + r8));
    // Depth nodes u in [0, 1] refined toward u = 1 (at x = u x stop radius the
    // surface is reached when the particle crosses the whole shield), the
    // escape angle split at its kink.
    constexpr std::array<double, 7> kSegment{0.0, 0.9, 0.99, 0.999, 0.9999, 0.99999, 1.0};
    constexpr int kDepth = 48;
    std::array<double, kDepth> u_node{}, u_weight{};
    for (int segment = 0; segment < 6; ++segment) {
        const double lo = kSegment[segment], hi = kSegment[segment + 1];
        for (int n = 0; n < 8; ++n) {
            u_node[8 * segment + n] = lo + 0.5 * (hi - lo) * (kGaussLegendre8Nodes[n] + 1.0);
            u_weight[8 * segment + n] = 0.5 * (hi - lo) * kGaussLegendre8Weights[n];
        }
    }
    std::array<double, 16> cosine{}, angle_weight{};
    for (int n = 0; n < 8; ++n) {
        const double half = 0.5 * (kGaussLegendre8Nodes[n] + 1.0);
        cosine[n] = -half;
        cosine[8 + n] = half;
        angle_weight[n] = angle_weight[8 + n] = 0.25 * kGaussLegendre8Weights[n];
    }
    std::array<std::array<Scalar, 16>, kDepth> chord{};
    for (int i = 0; i < kDepth; ++i) {
        const Scalar x = stop_radius * Scalar(u_node[i]);
        for (int j = 0; j < 16; ++j)
            chord[i][j] = -x * cosine[j]
                + scalar_sqrt(Scalar(1) - x * x * (1.0 - cosine[j] * cosine[j]));
    }
    // Escape averaged over the vacancy depths and directions.
    const auto escape = [&](const Scalar& removal) {
        Scalar e = Scalar(0);
        for (int i = 0; i < kDepth; ++i) {
            Scalar at_depth = Scalar(0);
            for (int j = 0; j < 16; ++j)
                at_depth += Scalar(angle_weight[j])
                    * scalar_exp(-areal_density_g_cm2 * removal * chord[i][j]);
            e += Scalar(u_weight[i]) * at_depth;
        }
        return e;
    };
    const auto removal_at = [&](double energy_keV) {
        const auto mu = cross_sections.evaluate(energy_keV, effective_Z, hydrogen_mass_fraction);
        return mu[0] + mu[1] + mu[3];
    };
    // Photon rays: 8 birth sites along the path x 16 directions, each with its
    // length L to the surface and the cosine eta at which it meets it (both
    // in g/cm^2 for L and eta x L).
    constexpr int kSite = 8, kRay = kSite * 16;
    std::array<Scalar, kRay> ray_length{}, ray_cosine{}, ray_depth{};
    std::array<double, kRay> ray_weight{};
    for (int i = 0; i < kSite; ++i) {
        const Scalar x = stop_radius * Scalar(0.5 * (kGaussLegendre8Nodes[i] + 1.0));
        for (int j = 0; j < 16; ++j) {
            const Scalar exit_cosine = scalar_sqrt(Scalar(1) - x * x * (1.0 - cosine[j] * cosine[j]));
            const int r = 16 * i + j;
            ray_length[r] = areal_density_g_cm2 * (-x * cosine[j] + exit_cosine);
            ray_cosine[r] = exit_cosine;
            ray_depth[r] = exit_cosine * ray_length[r];
            ray_weight[r] = 0.5 * kGaussLegendre8Weights[i] * angle_weight[j];
        }
    }
    const Scalar diameter = Scalar(2) * areal_density_g_cm2;
    // A line of removal mu_x escaping a photon of removal mu absorbed along
    // the rays, per unit photoelectric mass coefficient: the integral over
    // each ray of exp(-mu s) times the line's escape from the point s, which
    // is (E2(mu_x d) + E2(mu_x (2 t - d))) / 2 at depth d = eta (L - s) below
    // the plane tangent at the exit, out through that plane or back through
    // the diameter. Exact in a thin shield (escape 1), the half-space in a
    // thick one. E2 is the angle integral over c in (0, 1], done by GL8, which
    // leaves each ray in closed form.
    struct LineEscape {
        std::array<Scalar, 8 * kRay> front{}, back{};  // exp(-b L), exp((mu_x/c)(eta L - 2t))
        std::array<Scalar, 8> through_diameter{};       // exp(-2 t mu_x / c)
        std::array<Scalar, 8> inverse_cosine_mu{};      // mu_x / c
    };
    const auto line_escape = [&](const Scalar& removal_line) {
        LineEscape out;
        for (int c = 0; c < 8; ++c) {
            const double cosine_c = 0.5 * (kGaussLegendre8Nodes[c] + 1.0);
            out.inverse_cosine_mu[c] = removal_line / cosine_c;
            out.through_diameter[c] = scalar_exp(-out.inverse_cosine_mu[c] * diameter);
            for (int r = 0; r < kRay; ++r) {
                out.front[8 * r + c] = scalar_exp(-out.inverse_cosine_mu[c] * ray_depth[r]);
                out.back[8 * r + c] =
                    scalar_exp(out.inverse_cosine_mu[c] * (ray_depth[r] - diameter));
            }
        }
        return out;
    };
    const auto absorbed_escape = [&](const Scalar& removal, const LineEscape& line) {
        Scalar total = Scalar(0);
        for (int r = 0; r < kRay; ++r) {
            const Scalar& length = ray_length[r];
            const Scalar transmitted = scalar_exp(-removal * length);
            Scalar along = Scalar(0);
            for (int c = 0; c < 8; ++c) {
                const Scalar b = line.inverse_cosine_mu[c] * ray_cosine[r];
                // (exp(-b L) - exp(-mu L)) / (mu - b), its limit L exp(-b L)
                // (1 - delta/2 + delta^2/6) at delta = (mu - b) L -> 0.
                const Scalar delta = (removal - b) * length;
                const Scalar front = delta < Scalar(1e-4) && delta > Scalar(-1e-4)
                    ? length * line.front[8 * r + c]
                        * (Scalar(1) - Scalar(0.5) * delta + delta * delta / Scalar(6))
                    : (line.front[8 * r + c] - transmitted) / (removal - b);
                const Scalar back = (line.back[8 * r + c] - transmitted * line.through_diameter[c])
                    / (removal + b);
                along += Scalar(0.25 * kGaussLegendre8Weights[c]) * (front + back);
            }
            total += Scalar(ray_weight[r]) * along;
        }
        return total;
    };
    // Bremsstrahlung photons per unit ln k: the thick-target (Kramers)
    // spectrum over the particle's effective path in the shield, with the
    // brems bases' 1/X0 normalization.
    const Scalar photons_per_lnk = csda_range * path_fraction
        * smooth_inverse_radiation_length_cm2_g(effective_Z, hydrogen_mass_fraction, width);
    const Scalar non_h = Scalar(1) - hydrogen_mass_fraction;
    const Scalar removal511 = removal_at(electron_mass_keV);
    for (const auto& [Z, weight] : fluorescence_element_weights(effective_Z, width)) {
        const ElementFluorescence& element = element_fluorescence(Z);
        if (element.lines.empty()) continue;
        const auto vacancies = [&](double energy_keV) {
            const double pe = element_mass_cross_sections(Z, energy_keV).pe;
            const auto share = element_photoelectric_vacancy_shares(Z, energy_keV);
            std::array<Scalar, 4> out{};
            for (int k = 0; k < 4; ++k) out[k] = non_h * weight * Scalar(pe * share[k]);
            return out;
        };
        // Impact ionization over the path the particle spends in the shield.
        std::array<Scalar, 4> impact{};
        if (element.omega_K > 0.0)
            impact[0] = non_h * weight * path_fraction
                * Scalar(element_impact_xray_yield(Z, kinetic_energy_keV) / element.omega_K);
        // Photons born on the path: its bremsstrahlung, GL8 in ln k between
        // the edges (a photon below the lowest one leaves no vacancy with a
        // retained line) up to the particle's energy, as {photons, removal,
        // vacancies per unit absorbed path}; for a positron, the two
        // annihilation photons on their own.
        struct BornPhoton {
            Scalar photons, removal;
            std::array<Scalar, 4> vacancies;
        };
        std::vector<BornPhoton> brems;
        std::vector<double> bounds;
        for (const double edge : {element.l_edge_keV[2], element.l_edge_keV[1],
                                  element.l_edge_keV[0], element.k_edge_keV})
            if (edge > 0.0) bounds.push_back(std::max(edge, kFluorescenceMinimumKeV));
        std::sort(bounds.begin(), bounds.end());
        while (!bounds.empty() && !(bounds.back() < kinetic_energy_keV)) bounds.pop_back();
        bounds.push_back(kinetic_energy_keV);
        for (std::size_t piece = 0; piece + 1 < bounds.size(); ++piece) {
            const double lo = bounds[piece], hi = bounds[piece + 1];
            if (!(hi > lo)) continue;
            const double span = std::log(hi / lo);
            for (int q = 0; q < 8; ++q) {
                const double k = lo * std::exp(0.5 * span * (kGaussLegendre8Nodes[q] + 1.0));
                brems.push_back({photons_per_lnk * Scalar(0.5 * span * kGaussLegendre8Weights[q]
                                                          * (1.0 - k / kinetic_energy_keV)),
                                 removal_at(k), vacancies(k)});
            }
        }
        const auto v511 = vacancies(electron_mass_keV);
        // L vacancies from the element's own K x-rays absorbed, per K vacancy.
        std::array<Scalar, 3> k_cascade{};
        for (const auto& k_line : element.lines) {
            if (k_line.family == 2 || !(k_line.per_k > 0.0)) continue;
            const Scalar removal_k = removal_at(k_line.energy_keV);
            const Scalar absorbed = Scalar(k_line.per_k) * (Scalar(1) - escape(removal_k))
                / (removal_k + Scalar(std::numeric_limits<double>::min()));
            const auto v_k = vacancies(k_line.energy_keV);
            for (int s_l = 0; s_l < 3; ++s_l) k_cascade[s_l] += absorbed * v_k[1 + s_l];
        }
        for (const FluorescenceLineYield& line : element.lines) {
            const auto line_source = [&](const std::array<Scalar, 4>& v) {
                Scalar r = v[0] * Scalar(line.per_k) + v[1] * Scalar(line.per_l[0])
                    + v[2] * Scalar(line.per_l[1]) + v[3] * Scalar(line.per_l[2]);
                if (line.family == 2)
                    r += v[0] * (k_cascade[0] * Scalar(line.per_l[0])
                                 + k_cascade[1] * Scalar(line.per_l[1])
                                 + k_cascade[2] * Scalar(line.per_l[2]));
                return r;
            };
            const Scalar removal_line = removal_at(line.energy_keV);
            const LineEscape along_rays = line_escape(removal_line);
            typename ChargedFluorescenceBases<Scalar>::Line out;
            out.energy_keV = line.energy_keV;
            out.family = line.family;
            out.a = line_source(impact) * escape(removal_line);
            for (const BornPhoton& photon : brems)
                out.a += photon.photons * line_source(photon.vacancies)
                    * absorbed_escape(photon.removal, along_rays);
            if (is_positron)
                out.b = Scalar(2) * line_source(v511) * absorbed_escape(removal511, along_rays);
            result.a[line.family] += out.a;
            result.b[line.family] += out.b;
            result.lines.push_back(out);
        }
    }
    return result;
}

/// Build the 7-D "process" feature vector for one secondary component.
///
/// This is the SINGLE implementation: the Python fitter does not recompute it,
/// it reads the finished vector from
/// `ceelo_source_cross_section_features --emit-process-features`, which calls
/// straight through here. A second copy of this arithmetic would be a lockstep
/// hazard of the same family as the basis normalization that once cost an 83%
/// shape-TV error.
///
/// Axes: log incident energy, Z/92, log total optical depth, log parent-process
/// optical depth, log descendant absorption depth, log thickness in radiation
/// lengths, hydrogen mass fraction. `tau_total` is the caller's own
/// areal_density * mu_total, so the first three axes equal the "raw"
/// coordinate exactly.
template <typename Scalar>
void source_process_features(
    double E0_keV,
    const Scalar& areal_density_g_cm2,
    const Scalar& effective_Z,
    const Scalar& hydrogen_mass_fraction,
    const Scalar& tau_total,
    const std::array<Scalar, 4>& mu0,
    const SmoothElementCrossSection& cross_sections,
    double z_blend_width,
    SourceSecondaryProcess process,
    std::array<Scalar, kSourceRbfMaxFeatures>& feature)
{
    constexpr double electron_mass_keV = 510.99895;
    const double scale = std::log(25.0);

    const int fluorescence_family =
        process == SourceSecondaryProcess::FluorescenceKAlpha ? 0
        : process == SourceSecondaryProcess::FluorescenceKBeta ? 1
        : process == SourceSecondaryProcess::FluorescenceL ? 2 : -1;
    Scalar mu_parent = mu0[3];
    if (process == SourceSecondaryProcess::PhotoelectronBrems
        || fluorescence_family >= 0)
        mu_parent = mu0[0];
    else if (process == SourceSecondaryProcess::ComptonElectronBrems)
        mu_parent = mu0[1];

    // Annihilation descendants are monoenergetic 511 keV photons; every
    // bremsstrahlung descendant is a broad spectrum, represented by the
    // eight-downscatter reference energy the feature helper already uses; a
    // fluorescence family is its lines' yield-weighted removal.
    Scalar child_tau;
    if (fluorescence_family >= 0) {
        child_tau = areal_density_g_cm2
            * fluorescence_family_removal(
                fluorescence_family, effective_Z, hydrogen_mass_fraction,
                cross_sections);
    } else if (process == SourceSecondaryProcess::Pair511Atom
        || process == SourceSecondaryProcess::PairAnnihilationContinuum) {
        const auto mu511 = cross_sections.evaluate(
            electron_mass_keV, effective_Z, hydrogen_mass_fraction);
        child_tau =
            areal_density_g_cm2 * (mu511[0] + mu511[1] + mu511[3]);
    } else {
        const double reference8_keV =
            E0_keV / (1.0 + 8.0 * E0_keV / electron_mass_keV);
        const auto mu8 = cross_sections.evaluate(
            reference8_keV, effective_Z, hydrogen_mass_fraction);
        child_tau = areal_density_g_cm2 * (mu8[0] + mu8[3]);
    }

    const Scalar inverse_x0 = smooth_inverse_radiation_length_cm2_g(
        effective_Z, hydrogen_mass_fraction, z_blend_width);

    feature[0] = Scalar(std::log(E0_keV / 10.0) / std::log(1200.0));
    feature[1] = effective_Z / Scalar(92);
    feature[2] = scalar_log(Scalar(1) + tau_total) / Scalar(scale);
    feature[3] =
        scalar_log(Scalar(1) + areal_density_g_cm2 * mu_parent)
        / Scalar(scale);
    feature[4] = scalar_log(Scalar(1) + child_tau) / Scalar(scale);
    feature[5] =
        scalar_log(Scalar(1) + areal_density_g_cm2 * inverse_x0)
        / Scalar(scale);
    feature[6] = hydrogen_mass_fraction;
}

}  // namespace detail

template <typename Scalar>
SpectralMeasure<Scalar> SourceGeomFluxModel::evaluate(
    double E0_keV,
    const Scalar& areal_density_g_cm2,
    const Scalar& effective_Z,
    const std::vector<double>& output_bins_keV) const
{
    return evaluate(E0_keV, areal_density_g_cm2, effective_Z, Scalar(0),
                    output_bins_keV);
}

template <typename Scalar>
SpectralMeasure<Scalar> SourceGeomFluxModel::evaluate(
    double E0_keV,
    const Scalar& areal_density_g_cm2,
    const Scalar& effective_Z,
    const Scalar& hydrogen_mass_fraction,
    const std::vector<double>& output_bins_keV) const
{
    validate_bin_edges(output_bins_keV);
    SpectralMeasure<Scalar> result;
    if (!(E0_keV > 0.0))
        return result;

    const auto mu = cross_sections_.evaluate(E0_keV, effective_Z, hydrogen_mass_fraction);
    const auto photonuclear_mu = [&](double energy_keV) {
        return payload_.enable_photonuclear
            ? cross_sections_.photonuclear_mass_attenuation(
                energy_keV, effective_Z, hydrogen_mass_fraction, payload_.photonuclear)
            : Scalar(0);
    };
    const Scalar mu_pn = photonuclear_mu(E0_keV);
    const Scalar tau_pe = areal_density_g_cm2 * mu[0];
    const Scalar tau_cs = areal_density_g_cm2 * mu[1];
    const Scalar tau_rs = areal_density_g_cm2 * mu[2];
    const Scalar tau_pp = areal_density_g_cm2 * mu[3];
    const Scalar tau_pn = areal_density_g_cm2 * mu_pn;
    const Scalar tau_abs = tau_pe + tau_pp + tau_pn;
    const Scalar tau_transport =
        tau_abs + tau_cs
        + (payload_.enable_rayleigh ? tau_rs : Scalar(0));
    const Scalar direct =
        detail::scalar_exp(-tau_transport);

    result.atoms.push_back(
        {E0_keV, direct, SpectralComponent::Uncollided});

    if (payload_.enable_rayleigh) {
        const Scalar elastic = payload_.elastic_escape_scale
            * detail::scalar_exp(-tau_abs)
            * (Scalar(1) - detail::scalar_exp(-tau_rs))
            * detail::scalar_exp(
                Scalar(-payload_.elastic_compton_damping) * tau_cs);
        result.atoms.push_back(
            {E0_keV, elastic, SpectralComponent::ElasticLine});
    }

    // Exactly one Compton interaction in a centered homogeneous sphere.
    // Integrate the normalized KN × S(x,Z) angular law, the bound-electron
    // shell-profile energy distribution, and the interaction depth. The
    // Doppler term uses the same subshell occupancies, binding caps, and
    // truncated profiles as CeeLo's event sampler. Disabling it recovers the
    // free-electron moving support [E0/(1+2E0/m_e), E0].
    // Aliases for the shared rule in `detail`;
    // `detail::source_component_bases` uses the same nodes.
    constexpr const double* quadrature_nodes = detail::kGaussLegendre8Nodes;
    constexpr const double* quadrature_weights = detail::kGaussLegendre8Weights;
    constexpr double electron_mass_keV = 510.99895;
    const double single_min_energy =
        E0_keV / (1.0 + 2.0 * E0_keV / electron_mass_keV);
    const Scalar mu0_total =
        mu[0] + mu[1] + mu[3]
        + mu_pn
        + (payload_.enable_rayleigh ? mu[2] : Scalar(0));
    for (std::size_t b = 0; b + 1 < output_bins_keV.size(); ++b) {
        const double energy_lo = std::max(
            output_bins_keV[b],
            payload_.enable_compton_doppler
                ? 0.0 : single_min_energy);
        const double energy_hi =
            std::min(output_bins_keV[b + 1], E0_keV);
        if (!(energy_hi > energy_lo))
            continue;
        const auto energy_to_cosine = [&](double energy) {
            if (!(energy > 0.0))
                return -1.0;
            return std::max(-1.0, std::min(
                1.0, 1.0 - electron_mass_keV
                    * (1.0 / energy - 1.0 / E0_keV)));
        };
        const double free_cosine_lo = energy_to_cosine(energy_lo);
        const double free_cosine_hi = energy_to_cosine(energy_hi);
        std::array<double, 6> angle_panel_bounds{};
        int angle_panel_count = 1;
        if (payload_.enable_compton_doppler) {
            // Center a panel on the free-electron preimage of this energy
            // bin, add narrow panels on both sides to resolve the Doppler
            // shoulder, and retain two broad tail panels. This avoids the
            // discrete-angle ripples that a single global quadrature produces
            // when the Compton profile is narrow (especially low Z).
            const double free_width =
                std::max(0.0, free_cosine_hi - free_cosine_lo);
            const double shoulder_width =
                std::max(0.04, 8.0 * free_width);
            angle_panel_bounds = {
                -1.0,
                std::max(-1.0, free_cosine_lo - shoulder_width),
                free_cosine_lo,
                free_cosine_hi,
                std::min(1.0, free_cosine_hi + shoulder_width),
                1.0};
            angle_panel_count = 5;
        } else {
            angle_panel_bounds[0] = free_cosine_lo;
            angle_panel_bounds[1] = free_cosine_hi;
        }
        const double representative_energy =
            std::max(1e-6, 0.5 * (energy_lo + energy_hi));
        Scalar bin_yield = Scalar(0);
        Scalar rayleigh_mixed_bin_yield = Scalar(0);
        for (int panel_index = 0;
             panel_index < angle_panel_count; ++panel_index) {
            const double cosine_lo = angle_panel_bounds[panel_index];
            const double cosine_hi = angle_panel_bounds[panel_index + 1];
            if (!(cosine_hi > cosine_lo))
                continue;
            const double cosine_mid =
                0.5 * (cosine_lo + cosine_hi);
            const double cosine_half =
                0.5 * (cosine_hi - cosine_lo);
            for (int angle_index = 0; angle_index < 8; ++angle_index) {
            const double cosine = cosine_mid
                + cosine_half * quadrature_nodes[angle_index];
            const double free_scattered_energy =
                detail::compton_scattered_energy_keV(E0_keV, cosine);
            const Scalar energy_bin_probability =
                payload_.enable_compton_doppler
                ? cross_sections_.compton_doppler_bin_probability(
                    E0_keV, effective_Z, hydrogen_mass_fraction, cosine,
                    energy_lo, energy_hi)
                : Scalar(1);
            const double scattered_energy =
                payload_.enable_compton_doppler
                ? representative_energy : free_scattered_energy;
            const auto scattered_mu =
                cross_sections_.evaluate(scattered_energy, effective_Z, hydrogen_mass_fraction);
            const Scalar scattered_mu_total =
                scattered_mu[0] + scattered_mu[1] + scattered_mu[3]
                + photonuclear_mu(scattered_energy)
                + (payload_.enable_rayleigh
                    ? scattered_mu[2] : Scalar(0));
            Scalar depth_integral = Scalar(0);
            Scalar rayleigh_mixed_depth_integral = Scalar(0);
            for (int depth_index = 0; depth_index < 8; ++depth_index) {
                const double depth_fraction =
                    0.5 * (quadrature_nodes[depth_index] + 1.0);
                const double transverse2 = depth_fraction * depth_fraction
                    * (1.0 - cosine * cosine);
                const double exit_path_fraction =
                    -depth_fraction * cosine
                    + std::sqrt(std::max(0.0, 1.0 - transverse2));
                const Scalar probability_density =
                    areal_density_g_cm2 * mu[1]
                    * detail::scalar_exp(
                        -areal_density_g_cm2
                        * (mu0_total * depth_fraction
                           + scattered_mu_total * exit_path_fraction));
                depth_integral +=
                    0.5 * quadrature_weights[depth_index]
                    * probability_density;

                // Forward-Rayleigh approximation for exactly one Compton and
                // one-or-more Rayleigh interactions. Rayleigh collisions are
                // an unrestricted Poisson process along the incoming and
                // outgoing legs; absorption and further Compton collisions
                // remain forbidden. This preserves the strong process/Z
                // dependence absent from the former empirical mixed term.
                if (payload_.enable_rayleigh) {
                    const Scalar nonrayleigh_attenuation =
                        detail::scalar_exp(
                            -areal_density_g_cm2
                            * ((mu0_total - mu[2]) * depth_fraction
                               + (scattered_mu_total - scattered_mu[2])
                                   * exit_path_fraction));
                    const Scalar rayleigh_mean =
                        areal_density_g_cm2
                        * (mu[2] * depth_fraction
                           + scattered_mu[2] * exit_path_fraction);
                    const Scalar rayleigh_probability =
                        Scalar(1) - detail::scalar_exp(-rayleigh_mean);
                    rayleigh_mixed_depth_integral +=
                        0.5 * quadrature_weights[depth_index]
                        * areal_density_g_cm2 * mu[1]
                        * nonrayleigh_attenuation * rayleigh_probability;
                }
            }
            const Scalar angular_density =
                cross_sections_.compton_angular_density(
                    E0_keV, effective_Z, hydrogen_mass_fraction, cosine);
            bin_yield += cosine_half * quadrature_weights[angle_index]
                * angular_density * energy_bin_probability * depth_integral;
            rayleigh_mixed_bin_yield +=
                cosine_half * quadrature_weights[angle_index]
                * angular_density * energy_bin_probability
                * rayleigh_mixed_depth_integral;
            }
        }
        result.continuous.push_back({
            output_bins_keV[b], output_bins_keV[b + 1],
            payload_.single_escape_scale * bin_yield,
            SpectralComponent::SingleCompton});
        if (payload_.enable_rayleigh) {
            result.continuous.push_back({
                output_bins_keV[b], output_bins_keV[b + 1],
                payload_.compton_rayleigh_scale
                    * rayleigh_mixed_bin_yield,
                SpectralComponent::SingleComptonRayleigh});
        }
    }

    // Exactly two Compton interactions, including the spatial/angular
    // correlation in the sphere. The first collision is at radius x along the
    // original ray. The second is at distance y along the once-scattered ray;
    // after the second scatter its azimuth determines the remaining chord.
    // Integrating both collision densities and final survival gives a
    // parameter-free double-Compton component with exact energy support.
    //
    // The first-scatter cosine integral must be mapped onto the kinematic
    // support of each requested output-energy bin. A single fixed quadrature
    // over [-1, 1] makes a node's finite weight disappear when that node's
    // first-scattered energy falls below the bin. That produces nonphysical
    // steps at the first-scattered energies of the quadrature nodes. The
    // moving support depends only on E0 and the fixed output bins, so this
    // support-aware integration remains free of branches on fitted Scalar
    // values.
    static constexpr double q4_nodes[4] = {
        -0.8611363115940525752, -0.3399810435848562648, 0.3399810435848562648,
        0.8611363115940525752};
    static constexpr double q4_weights[4] = {
        0.3478548451374538574, 0.6521451548625461426, 0.6521451548625461426,
        0.3478548451374538574};
    static constexpr double q3_nodes[3] = {-0.7745966692414833770, 0.0,
                                           0.7745966692414833770};
    static constexpr double q3_weights[3] = {
        0.5555555555555555556, 0.8888888888888888889, 0.5555555555555555556};
    for (std::size_t b = 0; b + 1 < output_bins_keV.size(); ++b) {
        Scalar bin_yield = Scalar(0);
        Scalar rayleigh_mixed_bin_yield = Scalar(0);
        const double output_energy_lo = output_bins_keV[b];
        const double output_energy_hi = output_bins_keV[b + 1];
        std::array<double, 6> first_cosine_panel_bounds{};
        std::size_t first_cosine_panel_bound_count = 0;
        const auto add_first_energy_boundary = [&](double first_energy_keV) {
            if (!(first_energy_keV > single_min_energy &&
                  first_energy_keV < E0_keV)) {
                return;
            }
            first_cosine_panel_bounds[first_cosine_panel_bound_count++] =
                std::max(-1.0, std::min(1.0, 1.0 - electron_mass_keV *
                                                       (1.0 / first_energy_keV -
                                                        1.0 / E0_keV)));
        };
        const auto add_second_minimum_boundary = [&](double second_energy_keV) {
            // Invert E2_min(E1) = E1 / (1 + 2 E1 / m_e).
            // No finite boundary exists at or above m_e / 2.
            if (!(second_energy_keV > 0.0 &&
                  second_energy_keV < 0.5 * electron_mass_keV)) {
                return;
            }
            const double first_energy_keV =
                second_energy_keV /
                (1.0 - 2.0 * second_energy_keV / electron_mass_keV);
            add_first_energy_boundary(first_energy_keV);
        };
        first_cosine_panel_bounds[first_cosine_panel_bound_count++] = -1.0;
        first_cosine_panel_bounds[first_cosine_panel_bound_count++] = 1.0;
        add_first_energy_boundary(output_energy_lo);
        add_first_energy_boundary(output_energy_hi);
        add_second_minimum_boundary(output_energy_lo);
        add_second_minimum_boundary(output_energy_hi);
        std::sort(first_cosine_panel_bounds.begin(),
                  first_cosine_panel_bounds.begin() +
                      first_cosine_panel_bound_count);
        const auto unique_end = std::unique(
            first_cosine_panel_bounds.begin(),
            first_cosine_panel_bounds.begin() + first_cosine_panel_bound_count,
            [](double left, double right) {
                return std::abs(left - right) <= 1e-14;
            });
        first_cosine_panel_bound_count = static_cast<std::size_t>(
            unique_end - first_cosine_panel_bounds.begin());

        for (std::size_t first_panel_index = 0;
             first_panel_index + 1 < first_cosine_panel_bound_count;
             ++first_panel_index) {
            const double first_cosine_lo =
                first_cosine_panel_bounds[first_panel_index];
            const double first_cosine_hi =
                first_cosine_panel_bounds[first_panel_index + 1];
            if (!(first_cosine_hi > first_cosine_lo))
                continue;
            const double first_cosine_mid =
                0.5 * (first_cosine_lo + first_cosine_hi);
            const double first_cosine_half =
                0.5 * (first_cosine_hi - first_cosine_lo);
            const double representative_energy1 =
                detail::compton_scattered_energy_keV(E0_keV, first_cosine_mid);
            const double representative_minimum2 =
                representative_energy1 /
                (1.0 + 2.0 * representative_energy1 / electron_mass_keV);
            if (!(std::min(output_energy_hi, representative_energy1) >
                  std::max(output_energy_lo, representative_minimum2))) {
                continue;
            }

            for (int first_angle_index = 0; first_angle_index < 4;
                 ++first_angle_index) {
                const double cosine1 =
                    first_cosine_mid +
                    first_cosine_half * q4_nodes[first_angle_index];
                const double sine1 =
                    std::sqrt(std::max(0.0, 1.0 - cosine1 * cosine1));
                const double energy1 =
                    detail::compton_scattered_energy_keV(E0_keV, cosine1);
                const double minimum2 =
                    energy1 / (1.0 + 2.0 * energy1 / electron_mass_keV);
                const double energy_lo = std::max(output_energy_lo, minimum2);
                const double energy_hi = std::min(output_energy_hi, energy1);
                if (!(energy_hi > energy_lo))
                    continue;
                const auto second_energy_to_cosine = [&](double energy) {
                    return std::max(-1.0,
                                    std::min(1.0, 1.0 - electron_mass_keV *
                                                            (1.0 / energy -
                                                             1.0 / energy1)));
                };
                const double cosine2_lo = second_energy_to_cosine(energy_lo);
                const double cosine2_hi = second_energy_to_cosine(energy_hi);
                const double cosine2_mid = 0.5 * (cosine2_lo + cosine2_hi);
                const double cosine2_half = 0.5 * (cosine2_hi - cosine2_lo);
                const auto mu1 = cross_sections_.evaluate(energy1, effective_Z, hydrogen_mass_fraction);
                const Scalar mu1_total =
                    mu1[0] + mu1[1] + mu1[3] + photonuclear_mu(energy1) +
                    (payload_.enable_rayleigh ? mu1[2] : Scalar(0));
                const Scalar angular1 = cross_sections_.compton_angular_density(
                    E0_keV, effective_Z, hydrogen_mass_fraction, cosine1);
                for (int second_angle_index = 0; second_angle_index < 4;
                     ++second_angle_index) {
                    const double cosine2 =
                        cosine2_mid +
                        cosine2_half * q4_nodes[second_angle_index];
                    const double sine2 =
                        std::sqrt(std::max(0.0, 1.0 - cosine2 * cosine2));
                    const double energy2 =
                        detail::compton_scattered_energy_keV(energy1, cosine2);
                    const auto mu2 =
                        cross_sections_.evaluate(energy2, effective_Z, hydrogen_mass_fraction);
                    const Scalar mu2_total =
                        mu2[0] + mu2[1] + mu2[3] + photonuclear_mu(energy2) +
                        (payload_.enable_rayleigh ? mu2[2] : Scalar(0));
                    const Scalar angular2 =
                        cross_sections_.compton_angular_density(
                            energy1, effective_Z, hydrogen_mass_fraction, cosine2);
                    Scalar geometry_integral = Scalar(0);
                    Scalar rayleigh_geometry_integral = Scalar(0);
                    for (int first_depth_index = 0; first_depth_index < 3;
                         ++first_depth_index) {
                        const double x =
                            0.5 * (q3_nodes[first_depth_index] + 1.0);
                        const double first_exit =
                            -x * cosine1 +
                            std::sqrt(
                                std::max(0.0, 1.0 - x * x * sine1 * sine1));
                        for (int second_depth_index = 0; second_depth_index < 3;
                             ++second_depth_index) {
                            const double y =
                                0.5 * first_exit *
                                (q3_nodes[second_depth_index] + 1.0);
                            const double radius2 =
                                x * x + y * y + 2.0 * x * y * cosine1;
                            for (int azimuth_index = 0; azimuth_index < 3;
                                 ++azimuth_index) {
                                const double phi =
                                    3.14159265358979323846 *
                                    (q3_nodes[azimuth_index] + 1.0);
                                const double radial_dot_direction =
                                    cosine2 * (y + x * cosine1) -
                                    x * sine1 * sine2 * std::cos(phi);
                                const double final_exit =
                                    -radial_dot_direction +
                                    std::sqrt(std::max(
                                        0.0, radial_dot_direction *
                                                     radial_dot_direction +
                                                 1.0 - radius2));
                                const Scalar attenuation = detail::scalar_exp(
                                    -areal_density_g_cm2 *
                                    (mu0_total * x + mu1_total * y +
                                     mu2_total * final_exit));
                                geometry_integral +=
                                    0.125 * q3_weights[first_depth_index] *
                                    first_exit *
                                    q3_weights[second_depth_index] *
                                    q3_weights[azimuth_index] * attenuation;

                                if (payload_.enable_rayleigh) {
                                    const Scalar nonrayleigh_attenuation =
                                        detail::scalar_exp(
                                            -areal_density_g_cm2 *
                                            ((mu0_total - mu[2]) * x +
                                             (mu1_total - mu1[2]) * y +
                                             (mu2_total - mu2[2]) *
                                                 final_exit));
                                    const Scalar rayleigh_mean =
                                        areal_density_g_cm2 *
                                        (mu[2] * x + mu1[2] * y +
                                         mu2[2] * final_exit);
                                    rayleigh_geometry_integral +=
                                        0.125 * q3_weights[first_depth_index] *
                                        first_exit *
                                        q3_weights[second_depth_index] *
                                        q3_weights[azimuth_index] *
                                        nonrayleigh_attenuation *
                                        (Scalar(1) -
                                         detail::scalar_exp(-rayleigh_mean));
                                }
                            }
                        }
                    }
                    bin_yield += first_cosine_half *
                                 q4_weights[first_angle_index] *
                                 cosine2_half * q4_weights[second_angle_index] *
                                 angular1 * angular2 * areal_density_g_cm2 *
                                 mu[1] * areal_density_g_cm2 * mu1[1] *
                                 geometry_integral;
                    rayleigh_mixed_bin_yield +=
                        first_cosine_half *
                        q4_weights[first_angle_index] * cosine2_half *
                        q4_weights[second_angle_index] * angular1 * angular2 *
                        areal_density_g_cm2 * mu[1] * areal_density_g_cm2 *
                        mu1[1] * rayleigh_geometry_integral;
                }
            }
        }
        result.continuous.push_back({output_bins_keV[b], output_bins_keV[b + 1],
                                     payload_.double_escape_scale * bin_yield,
                                     SpectralComponent::DoubleCompton});
        if (payload_.enable_rayleigh) {
            result.continuous.push_back(
                {output_bins_keV[b], output_bins_keV[b + 1],
                 payload_.double_rayleigh_scale * rayleigh_mixed_bin_yield,
                 SpectralComponent::DoubleComptonRayleigh});
        }
    }

    // Residual three-plus-Compton term. With Rayleigh enabled, the high-order
    // Compton+Rayleigh split is either merged into this aggregate
    // (merge_higher_rayleigh) or carried as its own component below.
    const Scalar poisson_three_plus =
        Scalar(1) - detail::scalar_exp(-tau_cs)
            * (Scalar(1) + tau_cs + Scalar(0.5) * tau_cs * tau_cs);
    const double higher_reference_energy =
        E0_keV / (1.0 + 3.0 * E0_keV / electron_mass_keV);
    const auto higher_mu =
        cross_sections_.evaluate(higher_reference_energy, effective_Z, hydrogen_mass_fraction);
    const Scalar higher_absorption_tau =
        areal_density_g_cm2
        * (higher_mu[0] + higher_mu[3]
           + photonuclear_mu(higher_reference_energy));
    const Scalar higher_scatter_tau =
        areal_density_g_cm2
        * (higher_mu[1]
           + (payload_.enable_rayleigh
              ? higher_mu[2] : Scalar(0)));
    // Scattering is not photon loss: in the zero-absorption limit every
    // photon ultimately escapes, even at very large optical depth. Scattering
    // enters only by lengthening the absorbing random walk.
    const Scalar higher_escape =
        Scalar(1)
        / (Scalar(1) + Scalar(0.5) * higher_absorption_tau
           + Scalar(0.15) * higher_absorption_tau
               * (higher_absorption_tau + higher_scatter_tau));
    const Scalar higher_base =
        detail::scalar_exp(-tau_abs) * poisson_three_plus * higher_escape;

    const auto evaluate_higher_yield =
        [&](const MultipleScatterYieldModel& yield_model) {
            const Scalar energy_feature =
                Scalar(std::log(E0_keV / 10.0)
                       / std::log(1200.0));
            const Scalar z_feature = effective_Z / Scalar(92);
            Scalar weight_sum = Scalar(0);
            std::array<Scalar, kMultipleScatterYieldBasisSize>
                coefficients{};
            for (const auto& node : yield_model.nodes) {
                const Scalar de =
                    (energy_feature - node.feature[0])
                    / yield_model.bandwidth[0];
                const Scalar dz =
                    (z_feature - node.feature[1])
                    / yield_model.bandwidth[1];
                const Scalar weight = detail::scalar_exp(
                    Scalar(-0.5) * (de * de + dz * dz));
                weight_sum += weight;
                for (std::size_t k = 0;
                     k < kMultipleScatterYieldBasisSize; ++k) {
                    coefficients[k] +=
                        weight * node.coefficients[k];
                }
            }
            const Scalar denominator =
                weight_sum
                + Scalar(std::numeric_limits<double>::min());
            for (auto& coefficient : coefficients)
                coefficient /= denominator;
            const Scalar tau_total =
                areal_density_g_cm2 * mu0_total;
            const double tau_log_lo =
                std::log1p(yield_model.minimum_tau);
            const double tau_log_hi =
                std::log1p(yield_model.maximum_tau);
            const Scalar raw_coordinate =
                Scalar(2)
                    * (detail::scalar_log(Scalar(1) + tau_total)
                       - Scalar(tau_log_lo))
                    / Scalar(tau_log_hi - tau_log_lo)
                - Scalar(1);
            // The Chebyshev polynomial is fitted on tau in
            // [minimum_tau, maximum_tau] <-> [-1, 1]. Above the fitted range
            // its cubic extrapolation can overflow exp() and poison the
            // measure with inf * 0 = NaN at extreme optical depth (near
            // tau ~ 150). Freeze the conditional yield at its tau_max value
            // instead: unchanged for tau <= maximum_tau (and below
            // minimum_tau); the tail envelope owns the deep-tau regime.
            const Scalar coordinate =
                raw_coordinate < Scalar(1) ? raw_coordinate : Scalar(1);
            Scalar reduced_log_yield = coefficients[0];
            if (yield_model.degree >= 1) {
                Scalar previous = Scalar(1);
                Scalar current = coordinate;
                reduced_log_yield += coefficients[1] * current;
                for (std::size_t k = 2;
                     k <= yield_model.degree
                     && k < kMultipleScatterYieldBasisSize; ++k) {
                    const Scalar next =
                        Scalar(2) * coordinate * current - previous;
                    reduced_log_yield += coefficients[k] * next;
                    previous = current;
                    current = next;
                }
            }
            return tau_total * tau_total * tau_total
                * detail::scalar_exp(reduced_log_yield);
        };

    // The "raw" 3-D coordinate. Every 3-D component (higher Compton, higher
    // Compton+Rayleigh, photonuclear) queries exactly these three numbers.
    std::array<Scalar, kSourceRbfMaxFeatures> raw_feature{};
    raw_feature[0] =
        Scalar(std::log(E0_keV / 10.0) / std::log(1200.0));
    raw_feature[1] = effective_Z / Scalar(92);
    raw_feature[2] =
        detail::scalar_log(Scalar(1)
            + areal_density_g_cm2 * mu0_total)
        / std::log(25.0);

    // 7-D "process" coordinate, built by the one shared implementation the
    // feature helper binary also calls (see detail::source_process_features).
    const auto process_feature_for =
        [&](detail::SourceSecondaryProcess process) {
            std::array<Scalar, kSourceRbfMaxFeatures> feature{};
            detail::source_process_features(
                E0_keV, areal_density_g_cm2, effective_Z,
                hydrogen_mass_fraction,
                areal_density_g_cm2 * mu0_total, mu,
                cross_sections_, payload_.z_blend_width,
                process, feature);
            return feature;
        };

    // The six mechanistic baselines the fitted conditional yields multiply.
    // Computed by the one shared implementation that also feeds the Python
    // fitter through `ceelo_source_cross_section_features
    // --emit-component-bases`, so the two sides cannot drift.
    const auto& fluorescence_surrogates = payload_.secondary_surrogate;
    const bool fluorescence_fitted = payload_.enable_fluorescence
        && (!fluorescence_surrogates.fluorescence_k_alpha_atom.nodes.empty()
            || !fluorescence_surrogates.fluorescence_k_beta_atom.nodes.empty()
            || !fluorescence_surrogates.fluorescence_l_atom.nodes.empty());
    const auto component_bases = detail::source_component_bases(
        E0_keV, areal_density_g_cm2, effective_Z, hydrogen_mass_fraction,
        mu, mu0_total, cross_sections_, payload_.secondary_surrogate,
        fluorescence_fitted ? detail::kFluorescenceAllTerms : 0);

    const auto interpolate_residual =
        [&](const MultipleScatterRbfComponent& surrogate,
            const Scalar& base_yield,
            bool apply_conditional_yield,
            const std::array<Scalar, kSourceRbfMaxFeatures>*
                process_feature = nullptr) {
            std::array<Scalar, kMultipleScatterBasisSize> coefficients{};
            const std::size_t dimension = std::min<std::size_t>(
                surrogate.feature_dimension, kSourceRbfMaxFeatures);
            const auto& query =
                (process_feature != nullptr && dimension > 3)
                    ? *process_feature : raw_feature;
            const double ridge = surrogate.local_linear_ridge;
            const std::size_t slope_index = std::min<std::size_t>(
                surrogate.local_linear_feature_index,
                dimension > 0 ? dimension - 1 : 0);

            Scalar yield_weight_sum = Scalar(0);
            Scalar shape_weight_sum = Scalar(0);
            Scalar log_yield_sum = Scalar(0);
            std::array<Scalar, kMultipleScatterBasisSize>
                weighted_shape_sum{};
            // Local-linear accumulators along `slope_index`, centred on the
            // query so the fitted intercept IS the prediction there.
            Scalar yield_slope_sum = Scalar(0);
            Scalar yield_slope2_sum = Scalar(0);
            Scalar log_yield_slope_sum = Scalar(0);
            Scalar shape_slope_sum = Scalar(0);
            Scalar shape_slope2_sum = Scalar(0);
            std::array<Scalar, kMultipleScatterBasisSize>
                shape_slope_weighted_sum{};
            for (const auto& node : surrogate.nodes) {
                Scalar yield_distance2 = Scalar(0);
                Scalar shape_distance2 = Scalar(0);
                Scalar yield_slope_delta = Scalar(0);
                Scalar shape_slope_delta = Scalar(0);
                for (std::size_t index = 0; index < dimension; ++index) {
                    const Scalar yield_delta =
                        (query[index] - node.feature[index])
                        / surrogate.yield_bandwidth[index];
                    const Scalar shape_delta =
                        (query[index] - node.feature[index])
                        / surrogate.shape_bandwidth[index];
                    yield_distance2 += yield_delta * yield_delta;
                    shape_distance2 += shape_delta * shape_delta;
                    if (index == slope_index) {
                        yield_slope_delta = yield_delta;
                        shape_slope_delta = shape_delta;
                    }
                }
                const Scalar yield_weight = detail::scalar_exp(
                    Scalar(-0.5) * yield_distance2);
                const Scalar shape_weight = detail::scalar_exp(
                    Scalar(-0.5) * shape_distance2);
                yield_weight_sum += yield_weight;
                shape_weight_sum += shape_weight;
                log_yield_sum +=
                    yield_weight * node.log_conditional_yield;
                if (ridge > 0.0) {
                    yield_slope_sum += yield_weight * yield_slope_delta;
                    yield_slope2_sum +=
                        yield_weight * yield_slope_delta
                        * yield_slope_delta;
                    log_yield_slope_sum +=
                        yield_weight * yield_slope_delta
                        * node.log_conditional_yield;
                    shape_slope_sum += shape_weight * shape_slope_delta;
                    shape_slope2_sum +=
                        shape_weight * shape_slope_delta
                        * shape_slope_delta;
                }
                for (std::size_t k = 0;
                     k < kMultipleScatterBasisSize; ++k) {
                    weighted_shape_sum[k] +=
                        shape_weight * node.shape_coefficients[k];
                    if (ridge > 0.0) {
                        shape_slope_weighted_sum[k] +=
                            shape_weight * shape_slope_delta
                            * node.shape_coefficients[k];
                    }
                }
            }
            const Scalar yield_denominator =
                yield_weight_sum
                + Scalar(std::numeric_limits<double>::min());
            const Scalar shape_weight_denominator =
                shape_weight_sum
                + Scalar(std::numeric_limits<double>::min());
            Scalar interpolated_log_yield =
                log_yield_sum / yield_denominator;
            if (ridge > 0.0) {
                const Scalar stabilized =
                    yield_slope2_sum
                    + Scalar(ridge) * yield_weight_sum;
                const Scalar determinant =
                    yield_weight_sum * stabilized
                    - yield_slope_sum * yield_slope_sum
                    + Scalar(std::numeric_limits<double>::min());
                interpolated_log_yield =
                    (log_yield_sum * stabilized
                     - log_yield_slope_sum * yield_slope_sum)
                    / determinant;
            }
            // Mechanistic tail decay beyond the deepest populated optical
            // depth. The kernel field would otherwise extrapolate flat off
            // the last node and keep spurious yield in fully opaque shields.
            // The excess uses a smooth softplus rather than max(0, .) so the
            // whole expression stays differentiable for ceres::Jet.
            if (surrogate.tail_decay_kappa > 0.0
                && surrogate.tail_decay_tau > 0.0) {
                const Scalar excess =
                    areal_density_g_cm2 * mu0_total
                    - Scalar(surrogate.tail_decay_tau);
                constexpr double softness = 1.0e-6;
                const Scalar positive_excess =
                    Scalar(0.5)
                    * (excess
                       + detail::scalar_sqrt(
                           excess * excess
                           + Scalar(softness * softness)));
                interpolated_log_yield -=
                    Scalar(surrogate.tail_decay_kappa) * positive_excess;
            }
            const Scalar conditional_yield =
                detail::scalar_exp(interpolated_log_yield);
            Scalar coefficient_sum = Scalar(0);
            for (std::size_t k = 0;
                 k < kMultipleScatterBasisSize; ++k) {
                Scalar coefficient =
                    weighted_shape_sum[k] / shape_weight_denominator;
                if (ridge > 0.0) {
                    const Scalar stabilized =
                        shape_slope2_sum
                        + Scalar(ridge) * shape_weight_sum;
                    const Scalar determinant =
                        shape_weight_sum * stabilized
                        - shape_slope_sum * shape_slope_sum
                        + Scalar(std::numeric_limits<double>::min());
                    coefficient =
                        (weighted_shape_sum[k] * stabilized
                         - shape_slope_weighted_sum[k]
                             * shape_slope_sum)
                        / determinant;
                    // Local-linear coefficients can go negative where a
                    // kernel mean cannot; a smooth softplus keeps the shape
                    // nonnegative without a Jet-hostile branch.
                    constexpr double positivity_width = 1.0e-12;
                    coefficient = Scalar(0.5) * (
                        coefficient
                        + detail::scalar_sqrt(
                            coefficient * coefficient
                            + Scalar(
                                positivity_width * positivity_width)));
                }
                coefficients[k] = coefficient;
                coefficient_sum += coefficient;
            }
            const Scalar shape_denominator =
                coefficient_sum
                + Scalar(std::numeric_limits<double>::min());
            for (std::size_t k = 0;
                 k < kMultipleScatterBasisSize; ++k) {
                coefficients[k] *=
                    base_yield
                    * (apply_conditional_yield
                           ? conditional_yield : Scalar(1))
                    / shape_denominator;
            }
            return coefficients;
        };

    const Scalar higher_compton_rayleigh_base =
        payload_.enable_rayleigh
        ? detail::scalar_exp(-tau_abs)
            * poisson_three_plus
            * (Scalar(1) - detail::scalar_exp(
                -areal_density_g_cm2 * higher_mu[2]))
            * higher_escape
        : Scalar(0);
    const bool merged_higher =
        payload_.enable_rayleigh
        && payload_.multiple_scatter_surrogate.merge_higher_rayleigh
        && !payload_.multiple_scatter_surrogate
                .higher_compton.nodes.empty();
    const bool use_tau_yield =
        merged_higher
        && !payload_.multiple_scatter_surrogate
                .higher_yield.nodes.empty();
    if (payload_.multiple_scatter_surrogate.higher_compton.nodes.empty()) {
        const Scalar higher =
            payload_.multiple_escape_scale * higher_base;
        detail::append_bernstein_density(
            result, output_bins_keV, E0_keV, higher,
            payload_.multiple_compton_shape,
            SpectralComponent::HigherCompton);
    } else {
        const auto coefficients = interpolate_residual(
            payload_.multiple_scatter_surrogate.higher_compton,
            use_tau_yield
                ? evaluate_higher_yield(
                    payload_.multiple_scatter_surrogate.higher_yield)
                : (merged_higher
                    ? higher_base + higher_compton_rayleigh_base
                    : higher_base),
            !use_tau_yield);
        detail::append_warped_bernstein_density(
            result, output_bins_keV, E0_keV,
            payload_.multiple_scatter_surrogate
                .inverse_energy_warp_scale,
            payload_.multiple_scatter_surrogate
                .absorption_envelope_strength,
            payload_.multiple_scatter_surrogate
                .absorption_envelope_power,
            payload_.multiple_scatter_surrogate
                .absorption_envelope_max_ratio,
            payload_.multiple_scatter_surrogate
                .use_cubic_bspline,
            areal_density_g_cm2 * mu[0],
            coefficients,
            SpectralComponent::HigherCompton);
    }

    if (!payload_.enable_rayleigh) {
        // The far-field v1 source kernel removes coherent scattering at the
        // transport level, so no Rayleigh-tagged measure is emitted.
    } else if (merged_higher) {
        // The selected residual representation models the aggregate
        // three-plus-Compton continuum. Lower interaction orders remain
        // separately tallied; only this statistically sparse high-order
        // Rayleigh split is merged.
    } else if (payload_.multiple_scatter_surrogate
            .higher_compton_rayleigh.nodes.empty()) {
        const Scalar higher_compton_rayleigh =
            payload_.multiple_rayleigh_scale
            * higher_compton_rayleigh_base;
        detail::append_bernstein_density(
            result, output_bins_keV, E0_keV, higher_compton_rayleigh,
            payload_.compton_rayleigh_shape,
            SpectralComponent::HigherComptonRayleigh);
    } else {
        const auto coefficients = interpolate_residual(
            payload_.multiple_scatter_surrogate.higher_compton_rayleigh,
            higher_compton_rayleigh_base, true);
        detail::append_warped_bernstein_density(
            result, output_bins_keV, E0_keV,
            payload_.multiple_scatter_surrogate
                .inverse_energy_warp_scale,
            payload_.multiple_scatter_surrogate
                .absorption_envelope_strength,
            payload_.multiple_scatter_surrogate
                .absorption_envelope_power,
            payload_.multiple_scatter_surrogate
                .absorption_envelope_max_ratio,
            payload_.multiple_scatter_surrogate
                .use_cubic_bspline,
            areal_density_g_cm2 * mu[0],
            coefficients,
            SpectralComponent::HigherComptonRayleigh);
    }

    if (payload_.enable_pair_annihilation && E0_keV > 1021.9979) {
        // Still needed here for the absorption-envelope depth below; the
        // bases themselves come from the shared helper.
        const auto mu511 = cross_sections_.evaluate(510.99895, effective_Z, hydrogen_mass_fraction);
        const Scalar& pair_atom_base =
            component_bases.base[detail::kPair511Atom];
        Scalar annihilation_yield = pair_atom_base;
        const auto& secondary = payload_.secondary_surrogate;
        if (!secondary.pair_511_atom.nodes.empty()) {
            const auto feature = process_feature_for(
                detail::SourceSecondaryProcess::Pair511Atom);
            const auto coefficients = interpolate_residual(
                secondary.pair_511_atom, pair_atom_base, true,
                &feature);
            annihilation_yield = Scalar(0);
            for (const auto& coefficient : coefficients)
                annihilation_yield += coefficient;
        }
        result.atoms.push_back(
            {510.99895, annihilation_yield,
             SpectralComponent::PairAnnihilation});

        if (!secondary.pair_annihilation_continuum.nodes.empty()) {
            const Scalar& pair_scatter_base =
                component_bases.base[detail::kPairAnnihilationContinuum];
            const auto feature = process_feature_for(
                detail::SourceSecondaryProcess
                    ::PairAnnihilationContinuum);
            const auto coefficients = interpolate_residual(
                secondary.pair_annihilation_continuum,
                pair_scatter_base, true, &feature);
            detail::append_warped_bernstein_density(
                result, output_bins_keV, 510.99895,
                secondary.inverse_energy_warp_scale,
                secondary.absorption_envelope_strength,
                secondary.absorption_envelope_power,
                secondary.absorption_envelope_max_ratio,
                secondary.use_cubic_bspline,
                areal_density_g_cm2 * mu511[0],
                coefficients,
                SpectralComponent::PairAnnihilationContinuum);
        }
        if (!secondary.pair_annihilation_in_flight.nodes.empty()) {
            const Scalar& in_flight_base =
                component_bases.base[detail::kPairAnnihilationInFlight];
            const auto feature = process_feature_for(
                detail::SourceSecondaryProcess
                    ::PairAnnihilationInFlight);
            const auto coefficients = interpolate_residual(
                secondary.pair_annihilation_in_flight,
                in_flight_base, true, &feature);
            detail::append_warped_bernstein_density(
                result, output_bins_keV, E0_keV,
                secondary.inverse_energy_warp_scale,
                secondary.absorption_envelope_strength,
                secondary.absorption_envelope_power,
                secondary.absorption_envelope_max_ratio,
                secondary.use_cubic_bspline,
                tau_pe, coefficients,
                SpectralComponent::PairAnnihilationInFlight);
        }
    }

    if (payload_.enable_fluorescence) {
        // One fitted conditional yield per family, spread over the family's
        // lines (all blended elements) in proportion to their bases.
        const auto& secondary = payload_.secondary_surrogate;
        const std::array<const MultipleScatterRbfComponent*, kFluorescenceFamilies>
            surrogates{{&secondary.fluorescence_k_alpha_atom,
                        &secondary.fluorescence_k_beta_atom,
                        &secondary.fluorescence_l_atom}};
        const std::array<detail::SourceSecondaryProcess, kFluorescenceFamilies>
            processes{{detail::SourceSecondaryProcess::FluorescenceKAlpha,
                       detail::SourceSecondaryProcess::FluorescenceKBeta,
                       detail::SourceSecondaryProcess::FluorescenceL}};
        const std::array<SpectralComponent, kFluorescenceFamilies>
            labels{{SpectralComponent::SourceFluorescenceKAlpha,
                    SpectralComponent::SourceFluorescenceKBeta,
                    SpectralComponent::SourceFluorescenceL}};
        for (int family = 0; family < kFluorescenceFamilies; ++family) {
            const Scalar& family_base =
                component_bases.base[detail::kFluorescenceKAlpha + family];
            if (surrogates[family]->nodes.empty() || !(family_base > Scalar(0)))
                continue;
            const auto feature = process_feature_for(processes[family]);
            const auto coefficients = interpolate_residual(
                *surrogates[family], family_base, true, &feature);
            Scalar family_yield = Scalar(0);
            for (const auto& coefficient : coefficients)
                family_yield += coefficient;
            for (const auto& [energy_keV, line_base] :
                 component_bases.fluorescence_lines[family]) {
                result.atoms.push_back(
                    {energy_keV, family_yield * line_base / family_base,
                     labels[family]});
            }
        }
    }

    if (payload_.enable_source_electron_brems) {
        const auto& secondary = payload_.secondary_surrogate;
        const bool has_process_surrogates =
            !secondary.photoelectron_brems.nodes.empty()
            || !secondary.compton_electron_brems.nodes.empty()
            || !secondary.pair_electron_brems.nodes.empty();
        if (has_process_surrogates) {
            const std::array<
                std::pair<const MultipleScatterRbfComponent*,
                          SpectralComponent>, 3> components{{
                {&secondary.photoelectron_brems,
                 SpectralComponent::PhotoelectronBrems},
                {&secondary.compton_electron_brems,
                 SpectralComponent::ComptonElectronBrems},
                {&secondary.pair_electron_brems,
                 SpectralComponent::PairElectronBrems},
            }};
            const std::array<Scalar, 3> bases{{
                component_bases.base[detail::kPhotoelectronBrems],
                component_bases.base[detail::kComptonElectronBrems],
                component_bases.base[detail::kPairElectronBrems],
            }};
            const std::array<detail::SourceSecondaryProcess, 3>
                processes{{
                    detail::SourceSecondaryProcess::PhotoelectronBrems,
                    detail::SourceSecondaryProcess::ComptonElectronBrems,
                    detail::SourceSecondaryProcess::PairElectronBrems,
                }};
            for (std::size_t index = 0; index < components.size(); ++index) {
                if (components[index].first->nodes.empty())
                    continue;
                const auto feature =
                    process_feature_for(processes[index]);
                const auto coefficients = interpolate_residual(
                    *components[index].first, bases[index], true,
                    &feature);
                detail::append_warped_bernstein_density(
                    result, output_bins_keV, E0_keV,
                    secondary.inverse_energy_warp_scale,
                    secondary.absorption_envelope_strength,
                    secondary.absorption_envelope_power,
                    secondary.absorption_envelope_max_ratio,
                    secondary.use_cubic_bspline,
                    tau_pe, coefficients,
                    components[index].second);
            }
        } else {
            const Scalar charged_creation =
                Scalar(1) - detail::scalar_exp(-(tau_pe + tau_pp));
            const double energy_turn_on =
                1.0 - std::exp(
                    -std::max(0.0, E0_keV - 50.0) / 650.0);
            const Scalar brems = payload_.brems_yield_scale * energy_turn_on
                * charged_creation
                * detail::scalar_exp(Scalar(-0.25) * tau_cs);
            detail::append_bernstein_density(
                result, output_bins_keV, E0_keV, brems,
                payload_.brems_shape,
                SpectralComponent::SourceElectronBrems);
        }
    }
    if (payload_.enable_photonuclear
        && !payload_.photonuclear.cross_sections.empty()) {
        // Mechanistic removal baseline shared by both representations. The
        // Geant4 attenuation table supplies tau_pn; CeeLo transport has no
        // photonuclear truth of its own.
        const Scalar photonuclear_base =
            Scalar(1) - detail::scalar_exp(-tau_pn);
        if (payload_.photonuclear.secondary.nodes.empty()) {
            detail::append_bernstein_density(
                result, output_bins_keV, E0_keV,
                payload_.photonuclear.escaped_yield_scale
                    * photonuclear_base,
                payload_.photonuclear.shape,
                SpectralComponent::PhotonuclearSecondary);
        } else {
            // The fitted conditional yield absorbs the v1 global scale. The
            // escape factor carries the steep tau falloff mechanistically so
            // the interpolated residual stays nearly flat in tau.
            const auto coefficients = interpolate_residual(
                payload_.photonuclear.secondary,
                photonuclear_base
                    * detail::photonuclear_escape_factor(
                        areal_density_g_cm2 * mu0_total,
                        payload_.photonuclear.escape_kappa),
                true);
            detail::append_warped_bernstein_density(
                result, output_bins_keV, E0_keV,
                payload_.photonuclear.inverse_energy_warp_scale,
                0.0, 0.0, 1.0,
                payload_.photonuclear.use_cubic_bspline,
                areal_density_g_cm2 * mu[0],
                coefficients,
                SpectralComponent::PhotonuclearSecondary);
        }
    }
    return result;
}

template <typename Scalar>
SpectralMeasure<Scalar> ChargedParticleSourceFluxModel::evaluate(
    double kinetic_energy_keV,
    const Scalar& areal_density_g_cm2,
    const Scalar& effective_Z,
    bool is_positron,
    const std::vector<double>& output_bins_keV) const
{
    return evaluate(kinetic_energy_keV, areal_density_g_cm2, effective_Z,
                    Scalar(0), is_positron, output_bins_keV);
}

template <typename Scalar>
SpectralMeasure<Scalar> ChargedParticleSourceFluxModel::evaluate(
    double kinetic_energy_keV,
    const Scalar& areal_density_g_cm2,
    const Scalar& effective_Z,
    const Scalar& hydrogen_mass_fraction,
    bool is_positron,
    const std::vector<double>& output_bins_keV) const
{
    validate_bin_edges(output_bins_keV);
    SpectralMeasure<Scalar> result;
    if (!(kinetic_energy_keV > 0.0)
        || (!is_positron
            && !(kinetic_energy_keV
                 > payload_.min_photon_energy_keV))) {
        return result;
    }

    Scalar weight_sum = Scalar(0);
    Scalar range_sum = Scalar(0);
    Scalar electron_range_sum = Scalar(0);
    Scalar inverse_radiation_length_sum = Scalar(0);
    const double inverse_blend_width2 =
        1.0 / (payload_.z_blend_width * payload_.z_blend_width);
    for (int Z = 1; Z <= 92; ++Z) {
        const Scalar delta_Z = effective_Z - Scalar(Z);
        const Scalar weight = detail::scalar_exp(
            -delta_Z * delta_Z * inverse_blend_width2);
        weight_sum += weight;
        range_sum += weight * element_csda_range_g_cm2(
            Z, kinetic_energy_keV, is_positron);
        electron_range_sum += weight * element_csda_range_g_cm2(
            Z, kinetic_energy_keV, false);
        inverse_radiation_length_sum +=
            weight * element_inverse_radiation_length_cm2_g(Z);
    }
    const Scalar denominator =
        weight_sum + Scalar(std::numeric_limits<double>::min());
    Scalar csda_range_g_cm2 = range_sum / denominator;
    Scalar electron_csda_range_g_cm2 =
        electron_range_sum / denominator;
    Scalar inverse_radiation_length_cm2_g =
        inverse_radiation_length_sum / denominator;
    {
        // Hydrogen mixture: stopping powers add by mass fraction so CSDA
        // ranges combine reciprocally (R = R_blend / (1 + f_H (R_blend/R_H
        // - 1))); radiation length combines linearly in X0^-1. Both forms
        // are exact at f_H = 0.
        const double hydrogen_range = element_csda_range_g_cm2(
            1, kinetic_energy_keV, is_positron);
        const double hydrogen_electron_range = element_csda_range_g_cm2(
            1, kinetic_energy_keV, false);
        const double hydrogen_inverse_x0 =
            element_inverse_radiation_length_cm2_g(1);
        csda_range_g_cm2 = csda_range_g_cm2
            / (Scalar(1)
               + hydrogen_mass_fraction
                     * (csda_range_g_cm2 / Scalar(hydrogen_range)
                        - Scalar(1)));
        electron_csda_range_g_cm2 = electron_csda_range_g_cm2
            / (Scalar(1)
               + hydrogen_mass_fraction
                     * (electron_csda_range_g_cm2
                            / Scalar(hydrogen_electron_range)
                        - Scalar(1)));
        inverse_radiation_length_cm2_g += hydrogen_mass_fraction
            * (Scalar(hydrogen_inverse_x0)
               - inverse_radiation_length_cm2_g);
    }
    const Scalar range_ratio =
        areal_density_g_cm2
        / (csda_range_g_cm2
           + Scalar(std::numeric_limits<double>::min()));

    const std::array<Scalar, 3> query_feature{{
        Scalar(
            std::log(kinetic_energy_keV / 10.0)
            / std::log(1200.0)),
        effective_Z / Scalar(92),
        detail::scalar_log(Scalar(1) + range_ratio)
            / std::log1p(payload_.range_feature_scale),
    }};

    // Tensor cubic B-spline log-yield field plus a rank-R softmax
    // archetype-shape mixture on freeze-clamped coordinates. The returned
    // coefficients sum to base_yield * conditional_yield, distributed by a
    // normalized nonnegative shape. See ChargedParticleSplineField.
    const auto spline_tensor_field =
        [](const ChargedParticleSplineField& field,
           const std::vector<std::uint32_t>& grid,
           const std::vector<double>& coefficients,
           const std::array<Scalar, 3>& interpolation_feature) {
            std::array<std::array<Scalar, kChargedSplineMaxBasis>, 3>
                basis{};
            for (std::size_t axis = 0; axis < 3; ++axis) {
                const double lo = field.domain[2 * axis];
                const double hi = field.domain[2 * axis + 1];
                Scalar u = (interpolation_feature[axis] - Scalar(lo))
                    / Scalar(hi - lo);
                // Freeze outside the trained domain (constant field; the
                // mechanistic base carries the physical trend).
                u = u < Scalar(0)
                    ? Scalar(0)
                    : (u > Scalar(1) ? Scalar(1) : u);
                detail::charged_cubic_bspline_values_runtime(
                    grid[axis], u, basis[axis].data());
            }
            Scalar answer = Scalar(0);
            const std::size_t n1 = grid[1];
            const std::size_t n2 = grid[2];
            for (std::size_t i = 0; i < grid[0]; ++i) {
                for (std::size_t j = 0; j < n1; ++j) {
                    const Scalar b01 = basis[0][i] * basis[1][j];
                    const std::size_t offset = (i * n1 + j) * n2;
                    for (std::size_t k = 0; k < n2; ++k) {
                        answer += b01 * basis[2][k]
                            * coefficients[offset + k];
                    }
                }
            }
            return answer;
        };
    const auto spline_component =
        [&](const ChargedParticleSplineField& field,
            const Scalar& base_yield,
            const std::array<Scalar, 3>& interpolation_feature) {
            std::array<Scalar, kMultipleScatterBasisSize> answer{};
            const Scalar conditional_yield = detail::scalar_exp(
                spline_tensor_field(
                    field, field.yield_grid, field.yield_coefficients,
                    interpolation_feature));
            const Scalar total = base_yield * conditional_yield;
            if (field.shape_rank == 0) {
                // Yield-only component: a uniform fill keeps the sum
                // contract.
                for (auto& coefficient : answer)
                    coefficient = total
                        / Scalar(double(kMultipleScatterBasisSize));
                return answer;
            }
            std::array<Scalar, 16> log_odds{};
            const std::size_t rank = field.shape_rank;
            Scalar highest = Scalar(0);
            for (std::size_t r = 0; r + 1 < rank; ++r) {
                log_odds[r] = spline_tensor_field(
                    field, field.weight_grid, field.weight_log_odds[r],
                    interpolation_feature);
                if (log_odds[r] > highest)
                    highest = log_odds[r];
            }
            // Max-shifted softmax weights: exactly invariant to the shift,
            // so the value-part branch above is Jet-safe. The softmax
            // denominator cancels against the shape renormalization below
            // (archetypes are sum-1), so it is never formed explicitly.
            std::array<Scalar, 16> weights{};
            for (std::size_t r = 0; r < rank; ++r) {
                weights[r] = detail::scalar_exp(
                    (r + 1 < rank ? log_odds[r] : Scalar(0)) - highest);
            }
            // shape_k = sum_r softmax_r * archetype_rk, defensively
            // renormalized to sum 1 (archetypes are stored sum-1, so the
            // normalization is weight_sum up to rounding).
            Scalar shape_normalization = Scalar(0);
            for (std::size_t k = 0; k < kMultipleScatterBasisSize; ++k) {
                Scalar coefficient = Scalar(0);
                for (std::size_t r = 0; r < rank; ++r)
                    coefficient += weights[r] * field.shape_basis[r][k];
                answer[k] = coefficient;
                shape_normalization += coefficient;
            }
            for (auto& coefficient : answer)
                coefficient = total * coefficient
                    / (shape_normalization
                       + Scalar(std::numeric_limits<double>::min()));
            return answer;
        };
    // A component without a fitted field contributes nothing.
    const auto component_coefficients =
        [&](const ChargedParticleSplineField& spline,
            const Scalar& base_yield,
            const std::array<Scalar, 3>& interpolation_feature) {
            return spline.enabled()
                ? spline_component(
                      spline, base_yield, interpolation_feature)
                : std::array<Scalar, kMultipleScatterBasisSize>{};
        };

    // The effective charged-particle path tends to t for a thin source and to
    // the CSDA range for a fully contained particle. Multiplying by 1/X0 and
    // log(T/k_min) is the leading 1/k bremsstrahlung photon-yield integral.
    if (kinetic_energy_keV > payload_.min_photon_energy_keV) {
        const Scalar effective_path_g_cm2 =
            csda_range_g_cm2
            * (Scalar(1) - detail::scalar_exp(-range_ratio));
        const Scalar brems_base =
            effective_path_g_cm2 * inverse_radiation_length_cm2_g
            * std::log(
                kinetic_energy_keV / payload_.min_photon_energy_keV);
        const auto& brems_spline = is_positron
            ? payload_.positron_brems_spline
            : payload_.electron_brems_spline;
        const auto brems_coefficients =
            component_coefficients(
                brems_spline, brems_base, query_feature);
        std::vector<Scalar> brems_escape_envelope(
            output_bins_keV.size() - 1, Scalar(1));
        for (std::size_t bin = 0;
             bin + 1 < output_bins_keV.size(); ++bin) {
            const double energy_lo = std::max(
                payload_.min_photon_energy_keV,
                output_bins_keV[bin]);
            const double energy_hi = std::min(
                kinetic_energy_keV, output_bins_keV[bin + 1]);
            if (!(energy_hi > energy_lo))
                continue;
            const double center_keV = std::sqrt(energy_lo * energy_hi);
            const auto photon_xs =
                cross_sections_.evaluate(center_keV, effective_Z, hydrogen_mass_fraction);
            const Scalar absorption_tau =
                areal_density_g_cm2 * (photon_xs[0] + photon_xs[3]);
            brems_escape_envelope[bin] = detail::scalar_exp(
                -payload_.photon_absorption_path_fraction
                * absorption_tau);
        }
        detail::append_attenuated_log_energy_bspline_density(
            result, output_bins_keV, payload_.min_photon_energy_keV,
            kinetic_energy_keV, brems_escape_envelope,
            brems_coefficients,
            is_positron
                ? SpectralComponent::PrimaryPositronBrems
                : SpectralComponent::PrimaryElectronBrems);
    }

    // Shield K and L x-rays: one yield-only field per family, spread over the
    // family's lines by their mechanistic bases (a + contained x b).
    const std::array<const ChargedParticleSplineField*, kFluorescenceFamilies>
        fluorescence_splines = is_positron
            ? std::array<const ChargedParticleSplineField*, kFluorescenceFamilies>{{
                  &payload_.positron_fluorescence_k_alpha_spline,
                  &payload_.positron_fluorescence_k_beta_spline,
                  &payload_.positron_fluorescence_l_spline}}
            : std::array<const ChargedParticleSplineField*, kFluorescenceFamilies>{{
                  &payload_.electron_fluorescence_k_alpha_spline,
                  &payload_.electron_fluorescence_k_beta_spline,
                  &payload_.electron_fluorescence_l_spline}};
    const bool any_fluorescence = fluorescence_splines[0]->enabled()
        || fluorescence_splines[1]->enabled() || fluorescence_splines[2]->enabled();
    detail::ChargedFluorescenceBases<Scalar> fluorescence;
    if (any_fluorescence)
        fluorescence = detail::charged_fluorescence_bases(
            kinetic_energy_keV, areal_density_g_cm2, effective_Z,
            hydrogen_mass_fraction, is_positron, cross_sections_);
    const auto emit_fluorescence = [&](const Scalar& contained) {
        static constexpr std::array<SpectralComponent, kFluorescenceFamilies> labels{{
            SpectralComponent::SourceFluorescenceKAlpha,
            SpectralComponent::SourceFluorescenceKBeta,
            SpectralComponent::SourceFluorescenceL}};
        for (int family = 0; family < kFluorescenceFamilies; ++family) {
            if (!fluorescence_splines[family]->enabled()) continue;
            const Scalar base = fluorescence.a[family] + contained * fluorescence.b[family];
            if (!(base > Scalar(0))) continue;
            Scalar total = Scalar(0);
            for (const auto& coefficient :
                 component_coefficients(*fluorescence_splines[family], base, query_feature))
                total += coefficient;
            for (const auto& line : fluorescence.lines)
                if (line.family == family)
                    result.atoms.push_back(
                        {line.energy_keV, total * (line.a + contained * line.b) / base,
                         labels[family]});
        }
    };

    if (!is_positron) {
        if (any_fluorescence) emit_fluorescence(Scalar(0));
        return result;
    }

    // Containment is learned separately from photon escape. CeeLo's validated
    // source-escape gate is organized by radial depth in its extrapolated
    // electron range, including a sharp lower detour clamp. That clamp creates
    // moving energy/Z landmarks in positron containment. Use the same physical
    // coordinate for the RBF's third axis so interpolation follows the landmark
    // instead of smoothing across it in raw (energy, Z, t/R_pos) space.
    //
    // The 0.002-wide smooth clamps retain useful Jet derivatives at the fitted
    // effective Z while being narrow relative to the empirical gate's physics
    // uncertainty.
    const double energy_MeV = kinetic_energy_keV * 1.0e-3;
    const double detour_energy_modulation =
        1.0 / (
            1.0 + std::pow(energy_MeV / 0.35, 1.6));
    const Scalar raw_detour =
        Scalar(1)
        - Scalar(0.0645 * detour_energy_modulation) * effective_Z;
    constexpr double detour_smoothing = 0.002;
    const Scalar lower_delta = raw_detour - Scalar(0.10);
    const Scalar lower_clamped_detour =
        Scalar(0.5)
        * (raw_detour + Scalar(0.10)
           + detail::scalar_sqrt(
               lower_delta * lower_delta
               + Scalar(detour_smoothing * detour_smoothing)));
    const Scalar upper_delta =
        lower_clamped_detour - Scalar(0.98);
    const Scalar gate_detour =
        Scalar(0.5)
        * (lower_clamped_detour + Scalar(0.98)
           - detail::scalar_sqrt(
               upper_delta * upper_delta
               + Scalar(detour_smoothing * detour_smoothing)));
    const Scalar gate_range_ratio =
        areal_density_g_cm2
        / (electron_csda_range_g_cm2 * gate_detour
           + Scalar(std::numeric_limits<double>::min()));
    auto containment_feature = query_feature;
    containment_feature[2] =
        detail::scalar_log(Scalar(1) + gate_range_ratio)
        / std::log1p(payload_.range_feature_scale);

    // The exponential factor enforces exactly zero containment at zero
    // thickness.
    const Scalar containment_turn_on =
        Scalar(1) - detail::scalar_exp(Scalar(-100) * range_ratio);
    const auto containment_coefficients = component_coefficients(
        payload_.positron_containment_spline, containment_turn_on,
        containment_feature);
    Scalar contained_probability = Scalar(0);
    for (const auto& coefficient : containment_coefficients)
        contained_probability += coefficient;

    const auto mu511 =
        cross_sections_.evaluate(510.99895, effective_Z, hydrogen_mass_fraction);
    const Scalar tau511_absorb =
        areal_density_g_cm2 * (mu511[0] + mu511[3]);
    const Scalar tau511_compton =
        areal_density_g_cm2 * mu511[1];

    // The positron stop radius approaches the source surface for t/R <= 1
    // and R/t for thick sources. A smooth p-norm approximation keeps the
    // geometry Jet-compatible. Average isotropic photon escape from that
    // radius rather than assuming every annihilation occurs at the center.
    const Scalar range_ratio2 = range_ratio * range_ratio;
    const Scalar range_ratio4 = range_ratio2 * range_ratio2;
    const Scalar range_ratio8 = range_ratio4 * range_ratio4;
    const Scalar stop_radius_fraction = detail::scalar_exp(
        Scalar(-0.125)
        * detail::scalar_log(Scalar(1) + range_ratio8));
    static constexpr double quadrature_nodes[8] = {
        -0.9602898564975363, -0.7966664774136267,
        -0.5255324099163290, -0.1834346424956498,
         0.1834346424956498,  0.5255324099163290,
         0.7966664774136267,  0.9602898564975363};
    static constexpr double quadrature_weights[8] = {
        0.1012285362903763, 0.2223810344533745,
        0.3137066458778873, 0.3626837833783620,
        0.3626837833783620, 0.3137066458778873,
        0.2223810344533745, 0.1012285362903763};
    Scalar atom_escape_average = Scalar(0);
    Scalar compton_escape_average = Scalar(0);
    for (int index = 0; index < 8; ++index) {
        const double cosine = quadrature_nodes[index];
        const Scalar transverse2 =
            stop_radius_fraction * stop_radius_fraction
            * (1.0 - cosine * cosine);
        const Scalar path_fraction =
            -stop_radius_fraction * cosine
            + detail::scalar_sqrt(Scalar(1) - transverse2);
        const Scalar survives_absorption =
            detail::scalar_exp(-tau511_absorb * path_fraction);
        const Scalar survives_without_compton =
            detail::scalar_exp(-tau511_compton * path_fraction);
        atom_escape_average +=
            Scalar(0.5 * quadrature_weights[index])
            * survives_absorption * survives_without_compton;
        compton_escape_average +=
            Scalar(0.5 * quadrature_weights[index])
            * survives_absorption
            * (Scalar(1) - survives_without_compton);
    }
    if (any_fluorescence) emit_fluorescence(contained_probability);

    const Scalar positron_atom_base =
        Scalar(2) * contained_probability
        * atom_escape_average;
    const auto atom_coefficients = component_coefficients(
        payload_.positron_511_atom_spline, positron_atom_base,
        query_feature);
    Scalar atom_yield = Scalar(0);
    for (const auto& coefficient : atom_coefficients)
        atom_yield += coefficient;
    if (payload_.positron_511_atom_spline.enabled()) {
        result.atoms.push_back(
            {510.99895, atom_yield,
             SpectralComponent::PrimaryPositronAnnihilation});
    }

    const Scalar positron_continuum_base =
        Scalar(2) * contained_probability
        * compton_escape_average;
    const auto continuum_coefficients = component_coefficients(
        payload_.positron_annihilation_continuum_spline,
        positron_continuum_base, query_feature);
    detail::append_log_energy_bspline_density(
        result, output_bins_keV, payload_.min_photon_energy_keV,
        510.99895, continuum_coefficients,
        SpectralComponent::PrimaryPositronAnnihilationContinuum);

    const double in_flight_energy_onset =
        1.0 - std::exp(-kinetic_energy_keV / 510.99895);
    const Scalar in_flight_base =
        Scalar(2.0 * in_flight_energy_onset)
        * (Scalar(1) - detail::scalar_exp(-range_ratio));
    const auto in_flight_direct_coefficients = component_coefficients(
        payload_.positron_annihilation_in_flight_direct_spline,
        in_flight_base, query_feature);
    detail::append_linear_energy_bspline_density(
        result, output_bins_keV, 510.99895,
        kinetic_energy_keV + 510.99895,
        in_flight_direct_coefficients,
        SpectralComponent::PrimaryPositronAnnihilationInFlight);
    const auto in_flight_continuum_coefficients = component_coefficients(
        payload_.positron_annihilation_in_flight_continuum_spline,
        in_flight_base, query_feature);
    detail::append_log_energy_bspline_density(
        result, output_bins_keV, payload_.min_photon_energy_keV,
        kinetic_energy_keV + 510.99895,
        in_flight_continuum_coefficients,
        SpectralComponent::
            PrimaryPositronAnnihilationInFlightContinuum);
    return result;
}

template <typename Scalar>
SpectralMeasure<Scalar> DetectorPulseResponseModel::evaluate(
    double Ein_keV,
    const Scalar& cos_theta,
    const std::vector<double>& deposit_bins_keV) const
{
    validate_bin_edges(deposit_bins_keV);
    SpectralMeasure<Scalar> result;
    if (!(Ein_keV > 0.0)) {
        result.atoms.push_back(
            {0.0, Scalar(1), SpectralComponent::DetectorZeroDeposit});
        return result;
    }

    if (payload_.use_response_surrogate) {
        struct InterpolatedComponent {
            Scalar yield = Scalar(0);
            std::array<Scalar, kDetectorResponseBasisSize> shape{};
            std::size_t basis_size = kDetectorResponseBasisSize;
            std::array<std::uint32_t, detail::kDetectorMaxBlocks>
                block_sizes{};

            bool uses_blocks() const { return block_sizes[0] != 0; }
        };
        std::array<Scalar, 3> feature;
        // Feature coordinates: (log energy, polar angle / pi, passive-layer
        // attenuation landmark). Pair-history components replace the energy
        // coordinate with one expanded from the pair threshold.
        feature[0] =
            Scalar(std::log(Ein_keV / 10.0) / std::log(1200.0));
        {
            constexpr double pi = 3.14159265358979323846;
            // The tiny contraction keeps acos and its Jet derivative finite
            // at exactly end-on incidence without materially changing angle.
            feature[1] = detail::scalar_acos(
                Scalar(0.99999999)
                * detail::clamped_cosine(cos_theta)) / pi;
            const double passive_ratio =
                Ein_keV / payload_.passive_landmark_keV;
            feature[2] = Scalar(
                1.0 / (
                    1.0 + passive_ratio * passive_ratio
                        * passive_ratio));
        }
        const bool pair_enabled = Ein_keV > 1021.9979;
        std::array<Scalar, 3> pair_feature = feature;
        if (pair_enabled) {
            constexpr double threshold_keV = 1021.9979;
            constexpr double threshold_scale_keV = 50.0;
            pair_feature[0] = Scalar(
                std::log1p(
                    (Ein_keV - threshold_keV)
                    / threshold_scale_keV)
                / std::log1p(
                    (payload_.max_energy_keV - threshold_keV)
                    / threshold_scale_keV));
        }
        const auto interpolate =
            [&](const DetectorResponseRbfComponent& component,
                const std::array<Scalar, 3>& query_feature,
                bool enabled = true,
                bool interpolate_shape = true,
                bool requires_blocks = false) {
                InterpolatedComponent answer;
                if (!enabled || component.nodes.empty())
                    return answer;
                answer.basis_size = component.basis_size;
                answer.block_sizes = component.block_sizes;
                if (answer.basis_size > kDetectorResponseBasisSize) {
                    throw std::runtime_error(
                        "invalid detector component basis size");
                }
                // A fitted continuum must carry its block partition; an
                // unpartitioned one would silently contribute nothing and
                // break the measure's normalization.
                if (requires_blocks && !component.uses_blocks()) {
                    throw std::runtime_error(
                        "detector continuum requires a block partition");
                }
                if (component.uses_blocks()) {
                    std::size_t total = 0;
                    for (std::uint32_t size : component.block_sizes) {
                        if (size == 0)
                            continue;
                        if (size < 4
                            || size > detail::kDetectorMaxBlockSize) {
                            throw std::runtime_error(
                                "invalid detector block size");
                        }
                        total += size;
                    }
                    if (total != answer.basis_size) {
                        throw std::runtime_error(
                            "detector block sizes must sum to basis size");
                    }
                } else if (interpolate_shape
                           && !detail::supported_detector_basis_size(
                               answer.basis_size)) {
                    throw std::runtime_error(
                        "invalid detector component basis size");
                }
                Scalar yield_weight_sum = Scalar(0);
                Scalar shape_weight_sum = Scalar(0);
                Scalar weighted_log_yield = Scalar(0);
                for (const auto& node : component.nodes) {
                    Scalar yield_distance2 = Scalar(0);
                    Scalar shape_distance2 = Scalar(0);
                    for (std::size_t index = 0;
                         index < query_feature.size();
                         ++index) {
                        const Scalar yield_delta =
                            (query_feature[index] - node.feature[index])
                            / component.yield_bandwidth[index];
                        const Scalar shape_delta =
                            (query_feature[index] - node.feature[index])
                            / component.shape_bandwidth[index];
                        yield_distance2 += yield_delta * yield_delta;
                        shape_distance2 += shape_delta * shape_delta;
                    }
                    const Scalar yield_weight = detail::scalar_exp(
                        Scalar(-0.5) * yield_distance2);
                    const Scalar shape_weight = detail::scalar_exp(
                        Scalar(-0.5) * shape_distance2);
                    yield_weight_sum += yield_weight;
                    if (interpolate_shape)
                        shape_weight_sum += shape_weight;
                    weighted_log_yield +=
                        yield_weight * node.log_conditional_yield;
                    if (interpolate_shape) {
                        for (std::size_t index = 0;
                             index < answer.basis_size; ++index) {
                            answer.shape[index] +=
                                shape_weight
                                * node.shape_coefficients[index];
                        }
                    }
                }
                const Scalar tiny =
                    Scalar(std::numeric_limits<double>::min());
                answer.yield = detail::scalar_exp(
                    weighted_log_yield / (yield_weight_sum + tiny));
                if (interpolate_shape) {
                    Scalar coefficient_sum = Scalar(0);
                    for (std::size_t index = 0;
                         index < answer.basis_size; ++index) {
                        answer.shape[index] /=
                            shape_weight_sum + tiny;
                        coefficient_sum += answer.shape[index];
                    }
                    for (std::size_t index = 0;
                         index < answer.basis_size; ++index) {
                        answer.shape[index] /= coefficient_sum + tiny;
                    }
                }
                return answer;
            };
        auto zero = interpolate(
            payload_.zero_deposit, feature,
            !payload_.detected_rate_surface.enabled, false);
        auto fep = interpolate(
            payload_.full_energy, feature, true, false);
        // Pair model: pair_fraction is the odds of pair history among
        // detected pulses; pair outcomes are normalized within that parent
        // population. A payload without pair nodes has no pair response.
        auto pair_fraction = interpolate(
            payload_.pair_fraction, pair_feature, pair_enabled, false);
        auto pair_fep = interpolate(
            payload_.pair_full_energy, pair_feature, pair_enabled, false);
        auto single_escape =
            interpolate(
                payload_.single_escape, pair_feature,
                pair_enabled, false);
        auto double_escape =
            interpolate(
                payload_.double_escape, pair_feature,
                pair_enabled, false);
        auto annihilation_511 =
            interpolate(
                payload_.annihilation_511, pair_feature,
                pair_enabled, false);
        const bool k_escape_enabled =
            payload_.k_escape_alpha_energy_keV > 0.0
            && payload_.k_escape_beta_energy_keV > 0.0
            && Ein_keV > payload_.k_escape_edge_keV;
        // Each component carries its own support; they differ by up to 2x.
        const auto within_support = [&](std::size_t index) {
            const double limit = payload_.k_escape_support_keV[index];
            return k_escape_enabled
                && (limit <= 0.0 || Ein_keV <= limit);
        };
        auto k_alpha_escape = interpolate(
            payload_.k_alpha_escape, feature, within_support(0), false);
        auto k_beta_escape = interpolate(
            payload_.k_beta_escape, feature, within_support(1), false);
        auto k_alpha_absorbed = interpolate(
            payload_.k_alpha_absorbed, feature, within_support(2), false);
        auto k_beta_absorbed = interpolate(
            payload_.k_beta_absorbed, feature, within_support(3), false);
        auto photoelectric_escape =
            interpolate(
                payload_.detector_photoelectric_escape, feature);
        auto single = interpolate(
            payload_.detector_single_compton, feature, true, true, true);
        auto multiple = interpolate(
            payload_.detector_multiple_compton, feature, true, true, true);
        auto pair_no_compton = interpolate(
            payload_.detector_pair_no_compton, pair_feature,
            pair_enabled, true, true);
        auto pair_single_compton = interpolate(
            payload_.detector_pair_single_compton, pair_feature,
            pair_enabled, true, true);
        auto pair_multiple_compton = interpolate(
            payload_.detector_pair_multiple_compton, pair_feature,
            pair_enabled, true, true);
        Scalar detected_yield;
        if (payload_.detected_rate_surface.enabled) {
            const double energy_u =
                std::log(Ein_keV / 10.0) / std::log(1200.0);
            const double landmark = payload_.passive_landmark_keV;
            const double passive_at_ten =
                1.0 / (1.0 + std::pow(10.0 / landmark, 3.0));
            const double passive = 1.0 / (
                1.0 + std::pow(Ein_keV / landmark, 3.0));
            const double alpha =
                payload_.detected_rate_surface.energy_warp_alpha;
            const double rate_energy_coordinate =
                (energy_u + alpha * (passive_at_ten - passive))
                / (1.0 + alpha * passive_at_ten);
            const std::size_t energy_size =
                payload_.detected_rate_surface.energy_size;
            const std::size_t angle_size =
                payload_.detected_rate_surface.angle_size;
            if (energy_size < 2
                || energy_size > kDetectorRateMaxEnergyBasisSize
                || angle_size < 2
                || angle_size > kDetectorRateMaxAngleBasisSize) {
                throw std::runtime_error(
                    "invalid detected-rate surface dimensions");
            }
            std::array<double, kDetectorRateMaxEnergyBasisSize>
                energy_basis{};
            detail::bernstein_values_runtime(
                energy_size, rate_energy_coordinate, energy_basis.data());
            constexpr double pi = 3.14159265358979323846;
            const Scalar polar = detail::scalar_acos(
                Scalar(0.99999999)
                * detail::clamped_cosine(cos_theta)) / pi;
            const double exponent =
                payload_.detected_rate_surface.angle_warp_exponent;
            const double epsilon =
                payload_.detected_rate_surface.angle_warp_epsilon;
            const Scalar centred = Scalar(2) * polar - Scalar(1);
            const double power = 0.5 * (exponent - 1.0);
            const Scalar shaped =
                centred
                * detail::scalar_pow(
                    centred * centred + Scalar(epsilon * epsilon), power)
                / std::pow(1.0 + epsilon * epsilon, power);
            const Scalar angle_coordinate =
                Scalar(0.5) * (Scalar(1) + shaped);
            std::array<Scalar, kDetectorRateMaxAngleBasisSize> angle_basis{};
            detail::bernstein_values_runtime(
                angle_size, angle_coordinate, angle_basis.data());
            Scalar log_odds = Scalar(0);
            for (std::size_t angle = 0; angle < angle_size; ++angle) {
                for (std::size_t energy = 0; energy < energy_size; ++energy) {
                    const std::size_t index = angle * energy_size + energy;
                    log_odds += angle_basis[angle]
                        * energy_basis[energy]
                        * payload_.detected_rate_surface
                              .log_odds_coefficients[index];
                }
            }
            log_odds += rate_surface_anchor_log_odds(
                payload_.detected_rate_surface.anchor, Ein_keV, cos_theta);
            detected_yield =
                Scalar(1)
                / (Scalar(1) + detail::scalar_exp(-log_odds));
        } else {
            const Scalar detected_odds = zero.yield;
            detected_yield =
                detected_odds / (Scalar(1) + detected_odds);
        }
        zero.yield = Scalar(1) - detected_yield;
        const Scalar tiny =
            Scalar(std::numeric_limits<double>::min());
        {
        const Scalar pair_probability =
            pair_enabled
            ? (
                pair_fraction.yield
                * Scalar(
                    (Ein_keV - 1021.9979) / 510.99895))
                / (
                    Scalar(1)
                    + pair_fraction.yield
                        * Scalar(
                            (Ein_keV - 1021.9979)
                            / 510.99895))
            : Scalar(0);
        const Scalar nonpair_probability =
            Scalar(1) - pair_probability;
        const Scalar nonpair_normalization =
            fep.yield + photoelectric_escape.yield
            + k_alpha_escape.yield + k_beta_escape.yield
            + k_alpha_absorbed.yield + k_beta_absorbed.yield
            + single.yield + multiple.yield + tiny;
        const Scalar pair_normalization =
            pair_fep.yield + single_escape.yield
            + double_escape.yield + annihilation_511.yield
            + pair_no_compton.yield
            + pair_single_compton.yield
            + pair_multiple_compton.yield + tiny;
        const Scalar nonpair_scale =
            detected_yield * nonpair_probability
            / nonpair_normalization;
        const Scalar pair_scale =
            detected_yield * pair_probability
            / pair_normalization;
        fep.yield *= nonpair_scale;
        photoelectric_escape.yield *= nonpair_scale;
        k_alpha_escape.yield *= nonpair_scale;
        k_beta_escape.yield *= nonpair_scale;
        k_alpha_absorbed.yield *= nonpair_scale;
        k_beta_absorbed.yield *= nonpair_scale;
        single.yield *= nonpair_scale;
        multiple.yield *= nonpair_scale;
        pair_fep.yield *= pair_scale;
        single_escape.yield *= pair_scale;
        double_escape.yield *= pair_scale;
        annihilation_511.yield *= pair_scale;
        pair_no_compton.yield *= pair_scale;
        pair_single_compton.yield *= pair_scale;
        pair_multiple_compton.yield *= pair_scale;
        }
        result.atoms.push_back(
            {0.0, zero.yield, SpectralComponent::DetectorZeroDeposit});
        result.atoms.push_back(
            {Ein_keV, fep.yield, SpectralComponent::DetectorFullEnergy});
        if (k_escape_enabled) {
            // Yields are already zero outside each component's own support,
            // so emitting the atoms unconditionally here is exact; the
            // per-component gate lives in within_support above.
            result.atoms.push_back(
                {Ein_keV - payload_.k_escape_alpha_energy_keV,
                 k_alpha_escape.yield,
                 SpectralComponent::DetectorKAlphaEscape});
            result.atoms.push_back(
                {Ein_keV - payload_.k_escape_beta_energy_keV,
                 k_beta_escape.yield,
                 SpectralComponent::DetectorKBetaEscape});
            // Absorbed dead-layer fluorescence deposits the line energy
            // itself, so these atoms sit at a FIXED energy (like the 511
            // annihilation atom) rather than tracking Ein.
            result.atoms.push_back(
                {payload_.k_escape_alpha_energy_keV,
                 k_alpha_absorbed.yield,
                 SpectralComponent::DetectorKAlphaAbsorbed});
            result.atoms.push_back(
                {payload_.k_escape_beta_energy_keV,
                 k_beta_absorbed.yield,
                 SpectralComponent::DetectorKBetaAbsorbed});
        }
        if (pair_enabled) {
            result.atoms.push_back(
                {Ein_keV, pair_fep.yield,
                 SpectralComponent::DetectorPairFullEnergy});
            result.atoms.push_back(
                {Ein_keV - 510.99895, single_escape.yield,
                 SpectralComponent::DetectorSingleEscape});
            result.atoms.push_back(
                {Ein_keV - 1021.9979, double_escape.yield,
                 SpectralComponent::DetectorDoubleEscape});
            result.atoms.push_back(
                {510.99895, annihilation_511.yield,
                 SpectralComponent::DetectorAnnihilation511});
        }
        auto scale_shape = [](InterpolatedComponent& component) {
            for (std::size_t index = 0;
                 index < component.basis_size; ++index) {
                component.shape[index] *= component.yield;
            }
        };
        scale_shape(single);
        scale_shape(multiple);
        scale_shape(pair_no_compton);
        scale_shape(pair_single_compton);
        scale_shape(pair_multiple_compton);
        scale_shape(photoelectric_escape);
        detail::append_detector_linear_density(
            result, deposit_bins_keV, 0.0, Ein_keV,
            photoelectric_escape.shape,
            photoelectric_escape.basis_size,
            SpectralComponent::DetectorPhotoelectricEscape);
        // Each Compton and pair continuum is expanded in independent
        // B-spline blocks between its moving landmarks, so a genuine step
        // (the free-electron Compton edge, the escape edges) is
        // representable. A single spline over one warped coordinate could
        // not represent those steps, and with nonnegative coefficients it
        // traded a deficit below the edge against a spike at it.
        auto append_compton = [&](const InterpolatedComponent& component,
                                  SpectralComponent label) {
            detail::append_detector_block_density(
                result, deposit_bins_keV,
                detail::detector_compton_block_layout(
                    Ein_keV, component.block_sizes[0],
                    component.block_sizes[1]),
                component.shape, label);
        };
        auto append_pair = [&](const InterpolatedComponent& component,
                               SpectralComponent label) {
            detail::append_detector_block_density(
                result, deposit_bins_keV,
                detail::detector_pair_block_layout(
                    Ein_keV,
                    {component.block_sizes[0], component.block_sizes[1],
                     component.block_sizes[2], component.block_sizes[3]}),
                component.shape, label);
        };
        append_compton(single, SpectralComponent::SingleCompton);
        append_compton(multiple, SpectralComponent::MultipleCompton);
        if (pair_enabled) {
            append_pair(pair_no_compton,
                        SpectralComponent::DetectorPairNoCompton);
            append_pair(pair_single_compton,
                        SpectralComponent::DetectorPairSingleCompton);
            append_pair(pair_multiple_compton,
                        SpectralComponent::DetectorPairMultipleCompton);
        }
        return result;
    }

    const PartialMassAttenuation xs = mass_cross_sections(Ein_keV);
    const double mu_total = xs.pe + xs.cs + xs.rs + xs.pp;
    const double diameter = 2.0 * payload_.radius_cm;
    const Scalar c2 = cos_theta * cos_theta;
    const Scalar inverse_chord = detail::scalar_sqrt(
        c2 / (payload_.length_cm * payload_.length_cm)
        + (Scalar(1) - c2) / (diameter * diameter)
        + Scalar(1e-18));
    const Scalar chord = Scalar(1) / inverse_chord;
    const Scalar back_factor =
        Scalar(1) + payload_.back_incidence_scale
        * Scalar(0.5) * (Scalar(1) - cos_theta);
    const Scalar tau =
        mu_total * payload_.density_g_cm3 * chord * back_factor;
    const Scalar zero = detail::scalar_exp(-tau);
    const Scalar interacted = Scalar(1) - zero;

    result.atoms.push_back(
        {0.0, zero, SpectralComponent::DetectorZeroDeposit});

    const double inv_mu = 1.0 / std::max(mu_total, 1e-300);
    const double f_pe = xs.pe * inv_mu;
    const double f_cs = xs.cs * inv_mu;
    const double f_pp = xs.pp * inv_mu;
    const Scalar fep_raw = payload_.fep_base
        + payload_.fep_photoelectric_scale * f_pe
        + payload_.fep_reinteraction_scale * interacted;
    const Scalar single_raw =
        Scalar(0.25 + 1.6 * f_cs);
    const Scalar multiple_raw =
        Scalar(0.18 + 0.9 * f_cs) * (Scalar(0.3) + interacted);
    const double single_escape_raw =
        Ein_keV > 1021.9979 ? payload_.single_escape_scale * f_pp : 0.0;
    const double double_escape_raw =
        Ein_keV > 1021.9979 ? payload_.double_escape_scale * f_pp : 0.0;
    const Scalar normalization = fep_raw + single_raw + multiple_raw
        + single_escape_raw + double_escape_raw;

    const Scalar fep = interacted * fep_raw / normalization;
    const Scalar single = interacted * single_raw / normalization;
    const Scalar multiple = interacted * multiple_raw / normalization;
    result.atoms.push_back(
        {Ein_keV, fep, SpectralComponent::DetectorFullEnergy});
    if (single_escape_raw > 0.0) {
        result.atoms.push_back(
            {Ein_keV - 510.99895,
             interacted * single_escape_raw / normalization,
             SpectralComponent::DetectorSingleEscape});
    }
    if (double_escape_raw > 0.0) {
        result.atoms.push_back(
            {Ein_keV - 1021.9979,
             interacted * double_escape_raw / normalization,
             SpectralComponent::DetectorDoubleEscape});
    }
    detail::append_bernstein_density(
        result, deposit_bins_keV, Ein_keV, single,
        payload_.single_compton_shape, SpectralComponent::SingleCompton);
    detail::append_bernstein_density(
        result, deposit_bins_keV, Ein_keV, multiple,
        payload_.multiple_compton_shape, SpectralComponent::MultipleCompton);
    return result;
}

template <typename Scalar>
std::vector<Scalar> FastSpectrumModel::predict(
    const std::vector<EmissionLine>& lines,
    const Scalar& activity,
    const ShieldParameters<Scalar>& shield,
    const ResolutionParameters<Scalar>& resolution,
    const std::vector<double>& channels_keV) const
{
    return predict(
        lines, std::vector<ChargedParticleLine>{},
        activity, shield, resolution, channels_keV);
}

template <typename Scalar>
std::vector<Scalar> FastSpectrumModel::predict(
    const std::vector<EmissionLine>& lines,
    const std::vector<ChargedParticleLine>& charged_lines,
    const Scalar& activity,
    const ShieldParameters<Scalar>& shield,
    const ResolutionParameters<Scalar>& resolution,
    const std::vector<double>& channels_keV) const
{
    validate_bin_edges(channels_keV);
    const std::size_t n_channels = channels_keV.size() - 1;
    std::vector<Scalar> raw(n_channels, Scalar(0));
    struct WeightedPoint {
        double energy_keV;
        Scalar weight;
    };
    std::vector<WeightedPoint> points;
    points.reserve(lines.size() * 8);

    // The detector response is per photon crossing the beam envelope the
    // detector campaign aims at, so a far source's acceptance is that
    // envelope's projected area over 4 pi d^2. A pilot payload, which has no
    // geometry, falls back to the crystal's face.
    constexpr double pi = 3.14159265358979323846;
    const auto& envelope = payload_.detector.detected_rate_surface.anchor;
    double projected_area_cm2 =
        pi * payload_.detector.radius_cm * payload_.detector.radius_cm;
    if (envelope.envelope_radius_cm > 0.0) {
        const double c = std::abs(payload_.incidence_cos_theta);
        const double s = std::sqrt(std::max(0.0, 1.0 - c * c));
        projected_area_cm2 = pi * envelope.envelope_radius_cm
                * envelope.envelope_radius_cm * c
            + 2.0 * envelope.envelope_radius_cm
                * envelope.envelope_length_cm * s;
    }
    const double geometric_acceptance = projected_area_cm2
        / (4.0 * pi * payload_.distance_cm * payload_.distance_cm);
    const Scalar cos_theta = Scalar(payload_.incidence_cos_theta);

    // The point source's finite distance (FiniteDistanceDepth): per incident
    // energy, the tracking atoms take the full-energy depth's factor and the
    // rest of the detected measure whatever keeps the detected total on the
    // detected depth's factor.
    const auto& depth = payload_.detector.finite_distance;
    auto accumulate_detector = [&](double incident_energy,
                                   const Scalar& incident_weight_far,
                                   double source_bin_width_keV) {
        const auto detector_measure =
            detector_.evaluate(incident_energy, cos_theta, channels_keV);
        Scalar tracking_weight = incident_weight_far;
        Scalar other_weight = incident_weight_far;
        if (depth.enabled()
            && payload_.incidence_cos_theta >= kFiniteDistanceMinCos) {
            const auto at = [&](const std::vector<double>& values) {
                return detail::finite_distance_depth(depth.energies_keV, values,
                                                     incident_energy);
            };
            const double full = detail::finite_distance_factor(
                payload_.distance_cm, at(depth.full_energy_cm),
                at(depth.full_energy_radius_cm));
            const double detected = detail::finite_distance_factor(
                payload_.distance_cm, at(depth.detected_cm),
                at(depth.detected_radius_cm));
            Scalar total = Scalar(0), tracking = Scalar(0);
            for (const auto& atom : detector_measure.atoms) {
                if (atom.component == SpectralComponent::DetectorZeroDeposit)
                    continue;
                total += atom.yield;
                if (detail::detector_atom_tracks_incident(atom.component))
                    tracking += atom.yield;
            }
            for (const auto& bin : detector_measure.continuous)
                total += bin.integral;
            const Scalar rest = total - tracking;
            Scalar other = Scalar(detected);
            if (rest > Scalar(1e-12) * total)
                other = (Scalar(detected) * total - Scalar(full) * tracking)
                    / rest;
            tracking_weight = incident_weight_far * Scalar(full);
            other_weight = incident_weight_far * other;
        }
        for (const auto& atom : detector_measure.atoms) {
            if (atom.component == SpectralComponent::DetectorZeroDeposit)
                continue;
            const Scalar& incident_weight =
                detail::detector_atom_tracks_incident(atom.component)
                    ? tracking_weight : other_weight;
            // An atom at a constant offset from the incident energy (full
            // energy, escapes, K escape) has, when the incident "energy" is
            // really a source-quadrature bin, a uniform box of the bin width as
            // its locus. Depositing the box instead of a point is what makes
            // the source quadrature converge: a point renders the smooth
            // escaping continuum as nq artificial full-energy/escape peaks that
            // relocate under nq refinement. Fixed-energy atoms (511 keV, the
            // absorbed K lines) stay points.
            const bool tracks_incident =
                detail::detector_atom_tracks_incident(atom.component);
            if (source_bin_width_keV > 0.0 && tracks_incident) {
                const double half = 0.5 * source_bin_width_keV;
                const double lo = atom.energy_keV - half;
                const double hi = atom.energy_keV + half;
                if (!(hi > channels_keV.front())
                    || !(lo < channels_keV.back()))
                    continue;
                const Scalar box_weight = incident_weight * atom.yield;
                const auto it = std::upper_bound(
                    channels_keV.begin(), channels_keV.end(), lo);
                std::size_t j = static_cast<std::size_t>(std::max<
                    std::ptrdiff_t>(0, (it - channels_keV.begin()) - 1));
                for (; j < n_channels && channels_keV[j] < hi; ++j) {
                    const double overlap =
                        std::min(hi, channels_keV[j + 1])
                        - std::max(lo, channels_keV[j]);
                    if (overlap > 0.0)
                        raw[j] += box_weight
                            * (overlap / source_bin_width_keV);
                }
                continue;
            }
            points.push_back(
                {atom.energy_keV, incident_weight * atom.yield});
        }
        for (const auto& bin : detector_measure.continuous) {
            // Detector evaluation appends one entry per populated requested bin
            // for each continuum component, so map by bin edge rather than
            // relying on append order.
            const auto it = std::lower_bound(
                channels_keV.begin(), channels_keV.end(), bin.lower_keV);
            if (it != channels_keV.end() && *it == bin.lower_keV) {
                const std::size_t index =
                    static_cast<std::size_t>(it - channels_keV.begin());
                if (index < n_channels)
                    raw[index] += other_weight * bin.integral;
            }
        }
    };

    const int nq = std::max(2, payload_.source_quadrature_bins);
    // Shield fluorescence lines sit at the same element energies for every
    // source line: they are merged across all lines and the detector is
    // evaluated once per distinct energy after the loops.
    std::vector<WeightedPoint> fluorescence_atoms;
    // Shared tail for both source kinds: merge the measure's discrete atoms
    // by energy and feed them at their exact energies, then reduce the
    // continuum to nq quadrature bins whose yields are spread over the bin
    // width inside accumulate_detector.
    auto process_source_measure = [&](double support_keV,
                                      const SpectralMeasure<Scalar>&
                                          source_measure,
                                      const Scalar& line_weight) {
        // Detector transport is linear in incident photon yield. Merge source
        // process labels that share an energy before evaluating the detector;
        // otherwise a 16-bin source measure with many diagnostic components
        // needlessly evaluates the same detector response dozens of times.
        std::vector<WeightedPoint> source_atoms;
        source_atoms.reserve(source_measure.atoms.size());
        for (const auto& atom : source_measure.atoms) {
            if (atom.component == SpectralComponent::SourceFluorescenceKAlpha
                || atom.component == SpectralComponent::SourceFluorescenceKBeta
                || atom.component == SpectralComponent::SourceFluorescenceL) {
                auto merged = std::find_if(
                    fluorescence_atoms.begin(), fluorescence_atoms.end(),
                    [&](const WeightedPoint& point) {
                        return point.energy_keV == atom.energy_keV;
                    });
                if (merged == fluorescence_atoms.end())
                    fluorescence_atoms.push_back(
                        {atom.energy_keV, line_weight * atom.yield});
                else
                    merged->weight += line_weight * atom.yield;
                continue;
            }
            auto existing = std::find_if(
                source_atoms.begin(), source_atoms.end(),
                [&](const WeightedPoint& point) {
                    return point.energy_keV == atom.energy_keV;
                });
            if (existing == source_atoms.end())
                source_atoms.push_back({atom.energy_keV, atom.yield});
            else
                existing->weight += atom.yield;
        }
        for (const auto& atom : source_atoms)
            accumulate_detector(
                atom.energy_keV, line_weight * atom.weight, 0.0);
        std::vector<Scalar> source_bin_yields(
            static_cast<std::size_t>(nq), Scalar(0));
        for (const auto& bin : source_measure.continuous) {
            const int index = std::clamp(
                static_cast<int>(std::llround(
                    bin.lower_keV / support_keV * nq)),
                0, nq - 1);
            source_bin_yields[static_cast<std::size_t>(index)]
                += bin.integral;
        }
        const double bin_width = support_keV / nq;
        for (int index = 0; index < nq; ++index) {
            const double energy =
                support_keV * (static_cast<double>(index) + 0.5) / nq;
            accumulate_detector(
                energy,
                line_weight
                    * source_bin_yields[static_cast<std::size_t>(index)],
                bin_width);
        }
    };

    for (const EmissionLine& line : lines) {
        if (!(line.energy_keV > 0.0) || !(line.photons_per_decay > 0.0))
            continue;
        std::vector<double> source_edges(static_cast<std::size_t>(nq) + 1);
        for (int i = 0; i <= nq; ++i)
            source_edges[static_cast<std::size_t>(i)] =
                line.energy_keV * static_cast<double>(i) / nq;
        const auto source_measure = source_.evaluate(
            line.energy_keV, shield.areal_density_g_cm2,
            shield.effective_Z, shield.hydrogen_mass_fraction,
            source_edges);
        const Scalar line_weight = activity * line.photons_per_decay
            * geometric_acceptance;
        process_source_measure(line.energy_keV, source_measure, line_weight);
    }

    for (const ChargedParticleLine& line : charged_lines) {
        if (!(line.kinetic_energy_keV > 0.0)
            || (!line.is_positron
                && !(line.kinetic_energy_keV
                     > charged_source_.payload().min_photon_energy_keV))
            || !(line.particles_per_decay > 0.0)) {
            continue;
        }
        const double source_support_keV =
            line.kinetic_energy_keV
            + (line.is_positron ? 510.99895 : 0.0);
        std::vector<double> source_edges(
            static_cast<std::size_t>(nq) + 1);
        for (int i = 0; i <= nq; ++i) {
            source_edges[static_cast<std::size_t>(i)] =
                source_support_keV * static_cast<double>(i) / nq;
        }
        const auto source_measure = charged_source_.evaluate(
            line.kinetic_energy_keV,
            shield.areal_density_g_cm2,
            shield.effective_Z,
            shield.hydrogen_mass_fraction,
            line.is_positron,
            source_edges);
        const Scalar line_weight =
            activity * line.particles_per_decay * geometric_acceptance;
        process_source_measure(
            source_support_keV, source_measure, line_weight);
    }
    for (const auto& atom : fluorescence_atoms)
        accumulate_detector(atom.energy_keV, atom.weight, 0.0);

    std::vector<Scalar> answer(n_channels, Scalar(0));
    const int half_window =
        std::max(1, payload_.detector.resolution_half_window_channels);
    constexpr double kInvSqrt2 = 0.7071067811865475244;
    constexpr double kFwhmToSigma2 = 1.0 / (8.0 * 0.6931471805599453094);

    auto smear_point = [&](double energy, const Scalar& weight) {
        if (!(energy >= 0.0
              && energy >= channels_keV.front()
              && energy <= channels_keV.back()))
            return;
        const Scalar fwhm2 = resolution.a0_keV2
            + resolution.a1_keV * energy
            + resolution.a2 * energy * energy;
        const Scalar sigma =
            detail::scalar_sqrt(fwhm2 * kFwhmToSigma2 + Scalar(1e-18));
        const auto center_it =
            std::upper_bound(channels_keV.begin(), channels_keV.end(), energy);
        const int center = std::max(
            0, static_cast<int>(center_it - channels_keV.begin()) - 1);
        const int first = std::max(0, center - half_window);
        const int last = std::min(
            static_cast<int>(n_channels) - 1, center + half_window);
        const Scalar window_zlo =
            (Scalar(channels_keV[static_cast<std::size_t>(first)]) - energy)
            * kInvSqrt2 / sigma;
        const Scalar window_zhi =
            (Scalar(channels_keV[static_cast<std::size_t>(last + 1)]) - energy)
            * kInvSqrt2 / sigma;
        const Scalar normalization =
            Scalar(0.5)
                * (detail::scalar_erf(window_zhi)
                   - detail::scalar_erf(window_zlo))
            + Scalar(std::numeric_limits<double>::min());
        for (int j = first; j <= last; ++j) {
            const Scalar zlo =
                (Scalar(channels_keV[static_cast<std::size_t>(j)]) - energy)
                * kInvSqrt2 / sigma;
            const Scalar zhi =
                (Scalar(channels_keV[static_cast<std::size_t>(j + 1)]) - energy)
                * kInvSqrt2 / sigma;
            const Scalar probability =
                Scalar(0.5) * (detail::scalar_erf(zhi)
                               - detail::scalar_erf(zlo));
            answer[static_cast<std::size_t>(j)] +=
                weight * probability / normalization;
        }
    };

    for (std::size_t i = 0; i < n_channels; ++i) {
        const double center = 0.5 * (channels_keV[i] + channels_keV[i + 1]);
        smear_point(center, raw[i]);
    }
    for (const auto& point : points)
        smear_point(point.energy_keV, point.weight);
    return answer;
}

} // namespace ceelo::full_spectrum
