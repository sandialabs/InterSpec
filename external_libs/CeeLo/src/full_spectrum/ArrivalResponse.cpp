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

#include "full_spectrum/ArrivalResponse.h"

#include "geometry/Geometry.h"

#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <thread>
#include <tuple>

namespace ceelo::full_spectrum {

namespace {

constexpr double kPi = 3.14159265358979323846;

/// Integrals over [a, b] of the piecewise-linear hats on ascending `nodes`,
/// extended as constants beyond the end nodes; nonzero entries only.
template <typename Add>
void hat_integrals(const std::vector<double>& nodes, double a, double b,
                   Add add)
{
    const std::size_t n = nodes.size();
    if (!(b > a) || n == 0)
        return;
    if (n == 1) {
        add(0, b - a);
        return;
    }
    if (a < nodes.front())
        add(0, std::min(b, nodes.front()) - a);
    if (b > nodes.back())
        add(n - 1, b - std::max(a, nodes.back()));
    // Segments [nodes[i], nodes[i + 1]] that overlap [a, b].
    std::size_t i = static_cast<std::size_t>(
        std::upper_bound(nodes.begin(), nodes.end(), a) - nodes.begin());
    i = i == 0 ? 0 : i - 1;
    for (; i + 1 < n && nodes[i] < b; ++i) {
        const double x0 = nodes[i], x1 = nodes[i + 1];
        const double s0 = std::max(a, x0), s1 = std::min(b, x1);
        if (!(s1 > s0))
            continue;
        const double h = x1 - x0;
        add(i, ((x1 - s0) * (x1 - s0) - (x1 - s1) * (x1 - s1)) / (2.0 * h));
        add(i + 1, ((s1 - x0) * (s1 - x0) - (s0 - x0) * (s0 - x0)) / (2.0 * h));
    }
}

// The atoms whose energy moves with the incident energy, in a fixed order.
constexpr SpectralComponent kTracked[] = {
    SpectralComponent::DetectorFullEnergy,
    SpectralComponent::DetectorPairFullEnergy,
    SpectralComponent::DetectorSingleEscape,
    SpectralComponent::DetectorDoubleEscape,
    SpectralComponent::DetectorKAlphaEscape,
    SpectralComponent::DetectorKBetaEscape,
};
constexpr std::size_t kTrackedCount = sizeof(kTracked) / sizeof(kTracked[0]);

std::size_t tracked_slot(SpectralComponent component)
{
    for (std::size_t t = 0; t < kTrackedCount; ++t)
        if (kTracked[t] == component)
            return t;
    return kTrackedCount;
}

/// One response on deposit bins [0, length) of `edges`, the last bin closed
/// at the incident energy: the continuum and the fixed-energy atoms into
/// `row`; the tracking atoms' yields into `yield[slot]` and their offsets
/// into `offset[slot]`.
void tabulate(const DetectorPulseResponseModel& model, double energy,
              double mu, const std::vector<double>& edges, std::size_t length,
              float* row, double* yield, double* offset)
{
    std::vector<double> local(edges.begin(),
                              edges.begin() + static_cast<long>(length + 1));
    local.back() = std::min(local.back(), energy);
    const SpectralMeasure<double> measure =
        model.evaluate(energy, mu, local);
    for (const auto& bin : measure.continuous) {
        const auto at = std::upper_bound(local.begin(), local.end(),
                                         bin.lower_keV + 1e-9 * energy);
        const std::size_t d = static_cast<std::size_t>(at - local.begin()) - 1;
        if (d < length)
            row[d] += static_cast<float>(bin.integral);
    }
    for (const auto& atom : measure.atoms) {
        if (atom.component == SpectralComponent::DetectorZeroDeposit)
            continue;
        const std::size_t slot = tracked_slot(atom.component);
        if (slot < kTrackedCount) {
            yield[slot] += atom.yield;
            offset[slot] = energy - atom.energy_keV;
            continue;
        }
        if (atom.energy_keV < local.front() || atom.energy_keV >= local.back())
            continue;
        const auto at = std::upper_bound(local.begin(), local.end(),
                                         atom.energy_keV);
        row[static_cast<std::size_t>(at - local.begin()) - 1] +=
            static_cast<float>(atom.yield);
    }
}

} // namespace

ArrivalBand arrival_band(const std::string& name)
{
    if (name == "iso")
        return {-1.0, 1.0};
    if (name == "front")
        return {0.0, 1.0};
    if (name == "back")
        return {-1.0, 0.0};
    if (name == "side")
        return {-0.3, 0.3};
    const std::string prefix = "mu-band:";
    if (name.compare(0, prefix.size(), prefix) == 0) {
        const std::string range = name.substr(prefix.size());
        const auto comma = range.find(',');
        ArrivalBand band;
        std::size_t used_lo = 0, used_hi = 0;
        try {
            if (comma == std::string::npos)
                throw std::invalid_argument("no comma");
            band.mu_lo = std::stod(range.substr(0, comma), &used_lo);
            band.mu_hi = std::stod(range.substr(comma + 1), &used_hi);
        } catch (const std::exception&) {
            throw std::invalid_argument(
                "arrival band '" + name + "' is not mu-band:LO,HI");
        }
        if (used_lo != comma || used_hi != range.size() - comma - 1
            || !(band.mu_lo >= -1.0 && band.mu_lo < band.mu_hi
                 && band.mu_hi <= 1.0)) {
            throw std::invalid_argument(
                "arrival band '" + name + "' needs -1 <= LO < HI <= 1");
        }
        return band;
    }
    throw std::invalid_argument(
        "arrival band '" + name
        + "' is not iso, front, back, side, or mu-band:LO,HI");
}

double envelope_projected_area_cm2(double envelope_radius_cm,
                                   double envelope_length_cm, double mu)
{
    const double sine = std::sqrt(std::max(0.0, 1.0 - mu * mu));
    return kPi * envelope_radius_cm * envelope_radius_cm * std::abs(mu)
        + 2.0 * envelope_radius_cm * envelope_length_cm * sine;
}

ArrivalQuadrature arrival_quadrature(const ArrivalBand& band,
                                     double envelope_radius_cm,
                                     double envelope_length_cm,
                                     std::size_t nodes_per_hemisphere)
{
    if (!(band.mu_lo >= -1.0 && band.mu_lo < band.mu_hi
          && band.mu_hi <= 1.0))
        throw std::invalid_argument("arrival band needs -1 <= LO < HI <= 1");
    if (!(envelope_radius_cm > 0.0 && envelope_length_cm > 0.0))
        throw std::invalid_argument("arrival quadrature needs an envelope");
    if (nodes_per_hemisphere < 1)
        throw std::invalid_argument("arrival quadrature needs nodes");

    std::vector<double> x(nodes_per_hemisphere), w(nodes_per_hemisphere);
    detail::gauss_legendre_nodes(nodes_per_hemisphere, x.data(), w.data());

    std::vector<std::pair<double, double>> pieces;
    if (band.mu_lo < 0.0 && band.mu_hi > 0.0) {
        pieces = {{band.mu_lo, 0.0}, {0.0, band.mu_hi}};
    } else {
        pieces = {{band.mu_lo, band.mu_hi}};
    }

    ArrivalQuadrature rule;
    double total = 0.0;
    for (const auto& [lo, hi] : pieces) {
        // theta runs from acos(hi) to acos(lo); d(mu) = sin(theta) d(theta).
        const double theta_0 = std::acos(hi), theta_1 = std::acos(lo);
        const double half = 0.5 * (theta_1 - theta_0);
        const double middle = 0.5 * (theta_0 + theta_1);
        for (std::size_t k = 0; k < nodes_per_hemisphere; ++k) {
            const double theta = middle + half * x[k];
            const double mu = std::cos(theta);
            const double d_mu = half * w[k] * std::sin(theta);
            const double weight = d_mu * envelope_projected_area_cm2(
                envelope_radius_cm, envelope_length_cm, mu);
            rule.mu.push_back(mu);
            rule.weight.push_back(weight);
            total += weight;
        }
    }
    for (double& weight : rule.weight)
        weight /= total;
    rule.mean_projected_area_cm2 = total / (band.mu_hi - band.mu_lo);
    return rule;
}

SpectralMeasure<double> evaluate_arrival(
    const DetectorPulseResponseModel& model, double energy_keV,
    const ArrivalBand& band, const std::vector<double>& deposit_bins_keV,
    std::size_t nodes_per_hemisphere)
{
    const auto& anchor = model.payload().detected_rate_surface.anchor;
    if (!(anchor.envelope_radius_cm > 0.0 && anchor.envelope_length_cm > 0.0))
        throw std::invalid_argument(
            "evaluate_arrival: the payload carries no detector envelope");
    const ArrivalQuadrature rule = arrival_quadrature(
        band, anchor.envelope_radius_cm, anchor.envelope_length_cm,
        nodes_per_hemisphere);

    SpectralMeasure<double> answer;
    std::map<std::tuple<int, double, double>, std::size_t> bin_index;
    for (std::size_t node = 0; node < rule.mu.size(); ++node) {
        const double weight = rule.weight[node];
        const SpectralMeasure<double> measure =
            model.evaluate(energy_keV, rule.mu[node], deposit_bins_keV);
        for (const auto& atom : measure.atoms) {
            auto found = std::find_if(
                answer.atoms.begin(), answer.atoms.end(),
                [&](const SpectralAtom<double>& existing) {
                    return existing.component == atom.component
                        && existing.energy_keV == atom.energy_keV;
                });
            if (found == answer.atoms.end()) {
                answer.atoms.push_back(atom);
                answer.atoms.back().yield = weight * atom.yield;
            } else {
                found->yield += weight * atom.yield;
            }
        }
        for (const auto& bin : measure.continuous) {
            const auto key = std::make_tuple(
                static_cast<int>(bin.component), bin.lower_keV, bin.upper_keV);
            const auto [where, inserted] =
                bin_index.emplace(key, answer.continuous.size());
            if (inserted) {
                answer.continuous.push_back(bin);
                answer.continuous.back().integral = weight * bin.integral;
            } else {
                answer.continuous[where->second].integral +=
                    weight * bin.integral;
            }
        }
    }
    return answer;
}

ArrivalBeam::ArrivalBeam(const Geometry& geometry, const ArrivalBand& band)
    : geometry_(geometry), band_(band)
{
    if (!(band.mu_lo >= -1.0 && band.mu_lo <= band.mu_hi
          && band.mu_hi <= 1.0))
        throw std::invalid_argument("arrival band needs -1 <= LO <= HI <= 1");
    const auto z_extent = geometry.outer_z_extent();
    center_ = Eigen::Vector3d(
        0.0, 0.0, 0.5 * (z_extent.first + z_extent.second));
    const double half_z = 0.5 * (z_extent.second - z_extent.first);
    bound_cm_ = 1.01 * std::sqrt(
        geometry.outer_bounding_radius() * geometry.outer_bounding_radius()
        + half_z * half_z);
}

void ArrivalBeam::sample(std::mt19937_64& rng, Eigen::Vector3d& origin,
                         Eigen::Vector3d& direction)
{
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    do {
        const double mu =
            band_.mu_lo + (band_.mu_hi - band_.mu_lo) * uniform(rng);
        const double sine = std::sqrt(std::max(0.0, 1.0 - mu * mu));
        const double azimuth = 2.0 * kPi * uniform(rng);
        direction = Eigen::Vector3d(
            sine * std::cos(azimuth), sine * std::sin(azimuth), mu);
        const Eigen::Vector3d trial = std::abs(direction.z()) < 0.8
            ? Eigen::Vector3d::UnitZ() : Eigen::Vector3d::UnitX();
        const Eigen::Vector3d e1 = direction.cross(trial).normalized();
        const Eigen::Vector3d e2 = direction.cross(e1).normalized();
        const double r = bound_cm_ * std::sqrt(uniform(rng));
        const double phi = 2.0 * kPi * uniform(rng);
        origin = center_ - direction * (2.0 * bound_cm_)
            + e1 * (r * std::cos(phi)) + e2 * (r * std::sin(phi));
        ++tries_;
    } while (!geometry_.ray_hits_outer_boundary(origin, direction));
    ++accepted_;
}

double ArrivalBeam::mean_projected_area_cm2() const
{
    if (tries_ == 0)
        return 0.0;
    return kPi * bound_cm_ * bound_cm_ * static_cast<double>(accepted_)
        / static_cast<double>(tries_);
}

double ArrivalBeam::mean_projected_area_sigma_cm2() const
{
    if (tries_ == 0)
        return 0.0;
    const double p =
        static_cast<double>(accepted_) / static_cast<double>(tries_);
    return kPi * bound_cm_ * bound_cm_
        * std::sqrt(p * (1.0 - p) / static_cast<double>(tries_));
}

ArrivalResponseTable::ArrivalResponseTable(
    const DetectorPulseResponseModel& model,
    std::vector<double> deposit_edges_keV, ArrivalTableOptions options)
    : deposit_edges_(std::move(deposit_edges_keV))
{
    validate_bin_edges(deposit_edges_);
    if (deposit_edges_.front() < 0.0)
        throw std::invalid_argument("arrival table: negative deposit edge");
    if (!(options.energy_min_keV > 0.0
          && options.energy_max_keV > options.energy_min_keV
          && options.energy_step > 0.0 && options.nodes_per_hemisphere >= 1))
        throw std::invalid_argument("arrival table: invalid node options");
    const auto& anchor = model.payload().detected_rate_surface.anchor;
    if (!(anchor.envelope_radius_cm > 0.0 && anchor.envelope_length_cm > 0.0))
        throw std::invalid_argument(
            "arrival table: the payload carries no detector envelope");

    // Energy nodes: log-spaced, both ends exact.
    const double ratio = options.energy_max_keV / options.energy_min_keV;
    const auto steps = std::max<std::size_t>(1, static_cast<std::size_t>(
        std::ceil(std::log(ratio) / std::log1p(options.energy_step) - 1e-9)));
    for (std::size_t i = 0; i <= steps; ++i)
        energies_.push_back(options.energy_min_keV
                            * std::pow(ratio, static_cast<double>(i) / steps));
    energies_.back() = options.energy_max_keV;

    // Mu nodes: Gauss-Legendre in theta, front hemisphere then back.
    const std::size_t per = options.nodes_per_hemisphere;
    std::vector<double> x(per), w(per);
    detail::gauss_legendre_nodes(per, x.data(), w.data());
    for (int h = 0; h < 2; ++h) {
        const double t0 = h == 0 ? 0.0 : 0.5 * kPi;
        const double t1 = h == 0 ? 0.5 * kPi : kPi;
        for (std::size_t i = 0; i < per; ++i)
            mu_.push_back(std::cos(0.5 * (t0 + t1) + 0.5 * (t1 - t0) * x[i]));
    }
    envelope_radius_cm_ = anchor.envelope_radius_cm;
    envelope_length_cm_ = anchor.envelope_length_cm;
    for (std::size_t k = 0; k < mu_.size(); ++k)
        (mu_[k] > 0.0 ? front_index_ : back_index_).push_back(k);
    for (auto* order : {&front_index_, &back_index_}) {
        std::sort(order->begin(), order->end(),
                  [&](std::size_t a, std::size_t b) { return mu_[a] < mu_[b]; });
    }
    for (std::size_t k : front_index_)
        front_mu_.push_back(mu_[k]);
    for (std::size_t k : back_index_)
        back_mu_.push_back(mu_[k]);

    // Each energy node's row covers the deposit bins below it.
    const std::size_t n_mu = mu_.size(), n_energy = energies_.size();
    const auto bins_below = [&](double energy) {
        return static_cast<std::size_t>(
            std::lower_bound(deposit_edges_.begin(), deposit_edges_.end() - 1,
                             energy)
            - deposit_edges_.begin());
    };
    for (double energy : energies_)
        length_.push_back(bins_below(energy));
    offset_.resize(n_mu * n_energy);
    std::size_t total = 0;
    for (std::size_t k = 0; k < n_mu; ++k) {
        for (std::size_t j = 0; j < n_energy; ++j) {
            offset_[k * n_energy + j] = total;
            total += length_[j];
        }
    }
    rows_.assign(total, 0.0f);

    std::vector<double> lines = options.exact_line_energies_keV;
    std::sort(lines.begin(), lines.end());
    lines.erase(std::unique(lines.begin(), lines.end()), lines.end());
    for (double energy : lines) {
        if (!(energy > 0.0))
            throw std::invalid_argument("arrival table: invalid line energy");
        ExactLine line{energy, bins_below(energy), {}, {}};
        line.rows.assign(line.length * n_mu, 0.0f);
        line.yield.assign(kTrackedCount * n_mu, 0.0);
        lines_.push_back(std::move(line));
    }

    // Tabulate: tasks are energy nodes, then exact lines; every task writes
    // its own rows, yields and observed offsets.
    std::vector<double> yield(kTrackedCount * n_mu * n_energy, 0.0);
    const std::size_t tasks = n_energy + lines_.size();
    std::vector<double> offsets(kTrackedCount * tasks,
                                std::numeric_limits<double>::quiet_NaN());
    std::atomic<std::size_t> next{0};
    auto worker = [&]() {
        std::vector<double> y(kTrackedCount);
        for (std::size_t task = next++; task < tasks; task = next++) {
            double* seen = offsets.data() + task * kTrackedCount;
            for (std::size_t k = 0; k < n_mu; ++k) {
                std::fill(y.begin(), y.end(), 0.0);
                if (task < n_energy) {
                    const std::size_t j = task;
                    if (length_[j] == 0)
                        continue;
                    tabulate(model, energies_[j], mu_[k], deposit_edges_,
                             length_[j], rows_.data() + offset_[k * n_energy + j],
                             y.data(), seen);
                    for (std::size_t t = 0; t < kTrackedCount; ++t)
                        yield[(t * n_mu + k) * n_energy + j] = y[t];
                } else {
                    ExactLine& line = lines_[task - n_energy];
                    if (line.length == 0)
                        continue;
                    tabulate(model, line.energy_keV, mu_[k], deposit_edges_,
                             line.length, line.rows.data() + k * line.length,
                             y.data(), seen);
                    for (std::size_t t = 0; t < kTrackedCount; ++t)
                        line.yield[t * n_mu + k] = y[t];
                }
            }
        }
    };
    unsigned threads = options.threads ? options.threads
                                       : std::thread::hardware_concurrency();
    threads = std::max(1u, std::min<unsigned>(threads,
                                              static_cast<unsigned>(tasks)));
    std::vector<std::thread> pool;
    for (unsigned i = 1; i < threads; ++i)
        pool.emplace_back(worker);
    worker();
    for (auto& thread : pool)
        thread.join();

    // One offset per tracking atom (constant: 0, 511, 1022 keV, the K-line
    // family means); components the payload never produced are dropped.
    std::vector<std::size_t> kept_slots;
    for (std::size_t t = 0; t < kTrackedCount; ++t) {
        double offset = std::numeric_limits<double>::quiet_NaN();
        for (std::size_t task = 0; task < tasks; ++task) {
            const double seen = offsets[task * kTrackedCount + t];
            if (std::isnan(seen))
                continue;
            if (std::isnan(offset))
                offset = seen;
            else if (std::abs(seen - offset) > 1e-6 * std::max(1.0, std::abs(offset)))
                throw std::runtime_error(
                    "arrival table: a tracking atom's offset varies with energy");
        }
        if (std::isnan(offset))
            continue;
        Tracked tracked{kTracked[t], offset, {}};
        tracked.yield.assign(yield.begin() + static_cast<long>(t * n_mu * n_energy),
                             yield.begin() + static_cast<long>((t + 1) * n_mu * n_energy));
        tracked_.push_back(std::move(tracked));
        kept_slots.push_back(t);
    }
    for (ExactLine& line : lines_) {
        std::vector<double> kept(kept_slots.size() * n_mu);
        for (std::size_t i = 0; i < kept_slots.size(); ++i)
            for (std::size_t k = 0; k < n_mu; ++k)
                kept[i * n_mu + k] = line.yield[kept_slots[i] * n_mu + k];
        line.yield = std::move(kept);
    }
}

std::size_t ArrivalResponseTable::table_bytes() const
{
    std::size_t bytes = rows_.size() * sizeof(float);
    for (const auto& tracked : tracked_)
        bytes += tracked.yield.size() * sizeof(double);
    for (const auto& line : lines_)
        bytes += line.rows.size() * sizeof(float)
            + line.yield.size() * sizeof(double);
    return bytes;
}

std::vector<ArrivalResponseTable::Weight>
ArrivalResponseTable::energy_cell_weights(double lo, double hi) const
{
    const double tolerance = 1e-9 * energies_.back();
    if (lo < energies_.front() - tolerance || hi > energies_.back() + tolerance)
        throw std::invalid_argument(
            "fold: a continuum cell lies outside the table's energy range");
    std::vector<Weight> out;
    hat_integrals(energies_, lo, hi, [&](std::size_t j, double integral) {
        if (integral > 0.0)
            out.push_back({j, integral / (hi - lo)});
    });
    return out;
}

std::vector<ArrivalResponseTable::Weight>
ArrivalResponseTable::energy_point_weights(double energy) const
{
    const double tolerance = 1e-9 * energies_.back();
    if (energy < energies_.front() - tolerance
        || energy > energies_.back() + tolerance)
        throw std::invalid_argument(
            "fold: a line lies outside the table's energy range");
    const std::size_t n = energies_.size();
    std::size_t j = static_cast<std::size_t>(
        std::upper_bound(energies_.begin(), energies_.end(), energy)
        - energies_.begin());
    j = j == 0 ? 0 : std::min(j - 1, n - 2);
    const double t = std::clamp(
        (energy - energies_[j]) / (energies_[j + 1] - energies_[j]), 0.0, 1.0);
    std::vector<Weight> out;
    if (t < 1.0)
        out.push_back({j, 1.0 - t});
    if (t > 0.0)
        out.push_back({j + 1, t});
    return out;
}

std::vector<ArrivalResponseTable::Weight>
ArrivalResponseTable::mu_bin_weights(double lo, double hi) const
{
    // (1 / dmu) int_bin A(mu) R(mu) dmu, with R the polynomial in theta
    // through each hemisphere's Gauss-Legendre nodes (barycentric Lagrange)
    // and A exact: node k weighs int A l_k. The nodes integrate a whole
    // hemisphere spectrally; a piecewise-linear R between them converges
    // only as h^2, too slowly where the response turns over at grazing
    // incidence below ~100 keV (+0.47% at 59.5 keV with 16 nodes).
    std::vector<Weight> out;
    const double width = hi - lo;
    constexpr std::size_t kPoints = 32;
    static const auto rule = [] {
        std::array<std::array<double, kPoints>, 2> r{};
        detail::gauss_legendre_nodes(kPoints, r[0].data(), r[1].data());
        return r;
    }();
    const auto add_part = [&](const std::vector<double>& nodes,
                              const std::vector<std::size_t>& index,
                              double a, double b) {
        const std::size_t n = nodes.size();
        if (!(b > a) || n == 0)
            return;
        std::vector<double> theta(n), barycentric(n, 1.0), integral(n, 0.0);
        for (std::size_t i = 0; i < n; ++i)
            theta[i] = std::acos(nodes[i]);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                if (j != i)
                    barycentric[i] /= theta[i] - theta[j];
        const double t0 = std::acos(std::min(1.0, b));
        const double t1 = std::acos(std::max(-1.0, a));
        const double half = 0.5 * (t1 - t0), middle = 0.5 * (t0 + t1);
        std::vector<double> basis(n);
        for (std::size_t q = 0; q < kPoints; ++q) {
            const double th = middle + half * rule[0][q];
            const double mu = std::cos(th);
            const double dw = half * rule[1][q] * std::sin(th)
                * envelope_projected_area_cm2(envelope_radius_cm_,
                                              envelope_length_cm_, mu);
            double denominator = 0.0;
            std::size_t on_node = n;
            for (std::size_t i = 0; i < n; ++i) {
                const double gap = th - theta[i];
                if (gap == 0.0) {
                    on_node = i;
                    break;
                }
                basis[i] = barycentric[i] / gap;
                denominator += basis[i];
            }
            if (on_node < n) {
                integral[on_node] += dw;
                continue;
            }
            for (std::size_t i = 0; i < n; ++i)
                integral[i] += dw * basis[i] / denominator;
        }
        for (std::size_t i = 0; i < n; ++i)
            if (integral[i] != 0.0)
                out.push_back({index[i], integral[i] / width});
    };
    if (lo < 0.0)
        add_part(back_mu_, back_index_, lo, std::min(hi, 0.0));
    if (hi > 0.0)
        add_part(front_mu_, front_index_, std::max(lo, 0.0), hi);
    return out;
}

int ArrivalResponseTable::exact_line(double energy) const
{
    for (std::size_t l = 0; l < lines_.size(); ++l)
        if (std::abs(lines_[l].energy_keV - energy)
            <= 1e-9 * std::max(1.0, energy))
            return static_cast<int>(l);
    return -1;
}

} // namespace ceelo::full_spectrum
