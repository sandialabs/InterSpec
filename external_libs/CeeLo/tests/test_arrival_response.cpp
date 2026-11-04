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

#define BOOST_TEST_MODULE ArrivalResponseTests
#include <boost/test/unit_test.hpp>

#include "full_spectrum/ArrivalResponse.h"
#include "geometry/Geometry.h"
#include "materials/Material.h"
#include "transport/PhotonTransport.h"

#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>

namespace fs = ceelo::full_spectrum;

namespace {

constexpr double kPi = 3.14159265358979323846;

// The envelope's mean projected area over [lo, hi] with mu uniform:
// pi R^2 int|mu| + 2 R L int sqrt(1 - mu^2), in closed form.
double analytic_mean_area(double R, double L, double lo, double hi)
{
    auto abs_integral = [](double a, double b) {
        auto f = [](double m) { return 0.5 * m * std::abs(m); };
        return f(b) - f(a);
    };
    auto root_integral = [](double a, double b) {
        auto f = [](double m) {
            return 0.5 * (m * std::sqrt(1.0 - m * m) + std::asin(m));
        };
        return f(b) - f(a);
    };
    return (kPi * R * R * abs_integral(lo, hi)
            + 2.0 * R * L * root_integral(lo, hi))
        / (hi - lo);
}

const char* const kHpgeCoaxGeometry =
    "CEELO_GEOMETRY 1\n"
    "shape cylinder\n"
    "dimensions 4 8\n"
    "bullet_radius 0\n"
    "crystal_material 0\n"
    "bore 0.5 5.5 0\n"
    "dead_layer 0.05 0.08 0.08\n"
    "layer 1 0.1 0.1 -0.1 8.1\n"
    "symmetry axial\n"
    "reference_point crystal_face\n"
    "material 0 5.323 1 32 1 HPGe\n"
    "material 1 2.7 1 13 1 Al\n"
    "end\n";

} // namespace

BOOST_AUTO_TEST_CASE(named_bands_and_explicit_bands) {
    const auto iso = fs::arrival_band("iso");
    BOOST_TEST(iso.mu_lo == -1.0);
    BOOST_TEST(iso.mu_hi == 1.0);
    BOOST_TEST(fs::arrival_band("front").mu_lo == 0.0);
    BOOST_TEST(fs::arrival_band("back").mu_hi == 0.0);
    BOOST_TEST(fs::arrival_band("side").mu_lo == -0.3);
    const auto band = fs::arrival_band("mu-band:-0.25,0.75");
    BOOST_TEST(band.mu_lo == -0.25);
    BOOST_TEST(band.mu_hi == 0.75);
    for (const char* bad : {"sideways", "mu-band:0.5", "mu-band:0.5,0.2",
                            "mu-band:-1.5,0", "mu-band:0,1x", ""})
        BOOST_CHECK_THROW(fs::arrival_band(bad), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(quadrature_weights_are_the_projected_area) {
    const double R = 4.18, L = 8.3;
    for (const char* name : {"iso", "front", "back", "side"}) {
        const auto band = fs::arrival_band(name);
        const auto rule = fs::arrival_quadrature(band, R, L, 16);
        double total = 0.0;
        for (double weight : rule.weight)
            total += weight;
        BOOST_TEST(total == 1.0, boost::test_tools::tolerance(1e-12));
        // A trigonometric polynomial in theta per hemisphere: GL is exact.
        BOOST_TEST(rule.mean_projected_area_cm2
                       == analytic_mean_area(R, L, band.mu_lo, band.mu_hi),
                   boost::test_tools::tolerance(1e-10));
        // A band crossing mu = 0 gets a rule on each side of it.
        BOOST_TEST(rule.mu.size()
                   == (band.mu_lo < 0.0 && band.mu_hi > 0.0 ? 32u : 16u));
    }
    // A convex body's mean projected area over all directions is S / 4.
    const auto iso = fs::arrival_quadrature(fs::arrival_band("iso"), R, L);
    BOOST_TEST(iso.mean_projected_area_cm2
                   == (2.0 * kPi * R * R + 2.0 * kPi * R * L) / 4.0,
               boost::test_tools::tolerance(1e-10));
}

BOOST_AUTO_TEST_CASE(evaluate_arrival_is_the_weighted_fixed_angle_sum) {
    auto payload = fs::coaxial_hpge_payload();
    fs::set_detector_geometry(
        payload, ceelo::GeometryDescriptor::from_text(kHpgeCoaxGeometry));
    const fs::DetectorPulseResponseModel model(payload);
    const auto& anchor = model.payload().detected_rate_surface.anchor;
    std::vector<double> edges(129);
    const double energy = 300.0;
    for (std::size_t b = 0; b < edges.size(); ++b)
        edges[b] = energy * static_cast<double>(b) / 128.0;

    const auto band = fs::arrival_band("iso");
    const auto averaged = fs::evaluate_arrival(model, energy, band, edges);
    const auto rule = fs::arrival_quadrature(
        band, anchor.envelope_radius_cm, anchor.envelope_length_cm);
    double expected_total = 0.0, expected_fep = 0.0;
    for (std::size_t k = 0; k < rule.mu.size(); ++k) {
        const auto measure = model.evaluate(energy, rule.mu[k], edges);
        expected_total += rule.weight[k] * measure.integrated_yield();
        expected_fep += rule.weight[k] * measure.component_yield(
            fs::SpectralComponent::DetectorFullEnergy);
    }
    BOOST_TEST(averaged.integrated_yield() == expected_total,
               boost::test_tools::tolerance(1e-12));
    BOOST_TEST(averaged.component_yield(
                   fs::SpectralComponent::DetectorFullEnergy)
                   == expected_fep,
               boost::test_tools::tolerance(1e-12));
    // Merged: one atom per (component, energy), one bin per (component, bin).
    for (std::size_t i = 0; i < averaged.atoms.size(); ++i)
        for (std::size_t j = i + 1; j < averaged.atoms.size(); ++j)
            BOOST_TEST(!(averaged.atoms[i].component
                             == averaged.atoms[j].component
                         && averaged.atoms[i].energy_keV
                             == averaged.atoms[j].energy_keV));

    // Against a brute-force fine midpoint rule in mu: the per-hemisphere
    // rule has converged.
    double fine_total = 0.0, fine_weight = 0.0, fine_detected = 0.0;
    const int points = 2000;
    for (int i = 0; i < points; ++i) {
        const double mu = -1.0 + 2.0 * (i + 0.5) / points;
        const double area = fs::envelope_projected_area_cm2(
            anchor.envelope_radius_cm, anchor.envelope_length_cm, mu);
        const auto measure = model.evaluate(energy, mu, edges);
        fine_weight += area;
        fine_total += area * measure.integrated_yield();
        fine_detected += area * (1.0 - measure.component_yield(
            fs::SpectralComponent::DetectorZeroDeposit));
    }
    BOOST_TEST(averaged.integrated_yield() == fine_total / fine_weight,
               boost::test_tools::tolerance(1e-5));
    BOOST_TEST(1.0 - averaged.component_yield(
                         fs::SpectralComponent::DetectorZeroDeposit)
                   == fine_detected / fine_weight,
               boost::test_tools::tolerance(2e-5));

    payload.detected_rate_surface.anchor.envelope_radius_cm = 0.0;
    BOOST_CHECK_THROW(
        fs::evaluate_arrival(fs::DetectorPulseResponseModel(payload),
                             energy, band, edges),
        std::invalid_argument);
}

// Monte Carlo: the arrival beam's accepted photons carry the projected-area
// weight, so a band's sampled truth is the area-weighted combination of
// fixed-angle beams.
BOOST_AUTO_TEST_CASE(arrival_beam_truth_is_the_area_weighted_fixed_angle_sum) {
    const double R = 3.0, L = 5.0, energy = 200.0;
    ceelo::Material germanium = ceelo::make_HPGe();
    ceelo::Geometry geometry;
    geometry.set_detector(&germanium, ceelo::CylinderDims{R, L});
    ceelo::TransportConfig config;
    std::mt19937_64 rng(20261025);

    const auto band = fs::arrival_band("front");
    fs::ArrivalBeam beam(geometry, band);
    const std::uint64_t sampled = 150000;
    std::uint64_t detected = 0, steep = 0;
    Eigen::Vector3d origin, direction;
    for (std::uint64_t i = 0; i < sampled; ++i) {
        beam.sample(rng, origin, direction);
        if (direction.z() > 0.5)
            ++steep;
        const auto result = ceelo::transport_photon(
            origin, direction, energy, geometry, config, rng);
        if (result.energy_deposited_scoring > 0.0)
            ++detected;
    }
    BOOST_TEST(beam.accepted() == sampled);

    // The boundary's mean projected area: bounding disk x acceptance.
    const double area = analytic_mean_area(R, L, 0.0, 1.0);
    const double area_z = (beam.mean_projected_area_cm2() - area)
        / beam.mean_projected_area_sigma_cm2();
    BOOST_TEST_MESSAGE("mean projected area " << beam.mean_projected_area_cm2()
                       << " +- " << beam.mean_projected_area_sigma_cm2()
                       << " vs " << area << " (z " << area_z << ")");
    BOOST_TEST(std::abs(area_z) < 4.0);

    // Accepted directions follow A(mu): the share with mu > 0.5.
    const double steep_expected =
        0.5 * analytic_mean_area(R, L, 0.5, 1.0) / area;
    const double steep_share = static_cast<double>(steep) / sampled;
    const double steep_sigma =
        std::sqrt(steep_expected * (1.0 - steep_expected) / sampled);
    BOOST_TEST(std::abs(steep_share - steep_expected) < 4.0 * steep_sigma);

    // Fixed-angle beams at the quadrature nodes (a zero-width band is a beam
    // at that cosine; the detector is axisymmetric), area-weighted.
    const auto rule = fs::arrival_quadrature(band, R, L, 8);
    const std::uint64_t per_node = 25000;
    double combined = 0.0, combined_variance = 0.0;
    for (std::size_t k = 0; k < rule.mu.size(); ++k) {
        fs::ArrivalBeam fixed(geometry, fs::ArrivalBand{rule.mu[k], rule.mu[k]});
        std::uint64_t hits = 0;
        for (std::uint64_t i = 0; i < per_node; ++i) {
            fixed.sample(rng, origin, direction);
            const auto result = ceelo::transport_photon(
                origin, direction, energy, geometry, config, rng);
            if (result.energy_deposited_scoring > 0.0)
                ++hits;
        }
        const double p = static_cast<double>(hits) / per_node;
        combined += rule.weight[k] * p;
        combined_variance +=
            rule.weight[k] * rule.weight[k] * p * (1.0 - p) / per_node;
    }
    const double p_sampled = static_cast<double>(detected) / sampled;
    const double sigma = std::sqrt(
        combined_variance + p_sampled * (1.0 - p_sampled) / sampled);
    BOOST_TEST_MESSAGE("detected: sampled " << p_sampled << ", area-weighted "
                       << combined << " (z " << (p_sampled - combined) / sigma
                       << ")");
    BOOST_TEST(std::abs(p_sampled - combined) < 4.0 * sigma);
}

namespace {

// A forward-mode dual number: what a Jet needs from the fold, and no more.
struct Dual {
    double v = 0.0, d = 0.0;
    Dual() = default;
    Dual(double value, double derivative = 0.0) : v(value), d(derivative) {}
    Dual& operator+=(const Dual& o) { v += o.v; d += o.d; return *this; }
    friend Dual operator*(const Dual& a, double b) { return {a.v * b, a.d * b}; }
};

fs::DetectorPulseResponseModel hpge_model()
{
    auto payload = fs::coaxial_hpge_payload();
    fs::set_detector_geometry(
        payload, ceelo::GeometryDescriptor::from_text(kHpgeCoaxGeometry));
    return fs::DetectorPulseResponseModel(payload);
}

std::vector<double> uniform_edges(double top, std::size_t bins)
{
    std::vector<double> edges(bins + 1);
    for (std::size_t b = 0; b <= bins; ++b)
        edges[b] = top * static_cast<double>(b) / static_cast<double>(bins);
    return edges;
}

// Counts per bin of a measure (zero deposits left out). With slice > 0 the
// atoms that track the incident energy spread uniformly over the incident
// slice [E - slice/2, E + slice/2] they stand for.
std::vector<double> binned(const fs::SpectralMeasure<double>& m,
                           const std::vector<double>& edges, double scale,
                           double slice = 0.0)
{
    std::vector<double> out(edges.size() - 1, 0.0);
    for (const auto& bin : m.continuous) {
        const auto at = std::upper_bound(edges.begin(), edges.end(),
                                         bin.lower_keV + 1e-9);
        out[static_cast<std::size_t>(at - edges.begin()) - 1] += scale * bin.integral;
    }
    for (const auto& atom : m.atoms) {
        if (atom.component == fs::SpectralComponent::DetectorZeroDeposit)
            continue;
        if (slice > 0.0 && fs::detail::detector_atom_tracks_incident(atom.component)) {
            const double lo = atom.energy_keV - 0.5 * slice;
            const double hi = atom.energy_keV + 0.5 * slice;
            for (std::size_t b = 0; b + 1 < edges.size(); ++b) {
                const double overlap = std::min(hi, edges[b + 1]) - std::max(lo, edges[b]);
                if (overlap > 0.0)
                    out[b] += scale * atom.yield * overlap / slice;
            }
            continue;
        }
        if (atom.energy_keV >= edges.back())
            continue;
        const auto at = std::upper_bound(edges.begin(), edges.end(), atom.energy_keV);
        out[static_cast<std::size_t>(at - edges.begin()) - 1] += scale * atom.yield;
    }

    return out;
}

double total_variation(const std::vector<double>& a, const std::vector<double>& b)
{
    double sa = 0.0, sb = 0.0, tv = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) { sa += a[i]; sb += b[i]; }
    for (std::size_t i = 0; i < a.size(); ++i)
        tv += std::abs(a[i] / sa - b[i] / sb);
    return 0.5 * tv;
}

double sum(const std::vector<double>& v)
{
    double s = 0.0;
    for (double x : v) s += x;
    return s;
}

} // namespace

BOOST_AUTO_TEST_CASE(arrival_table_folds_a_line_as_evaluate_arrival_does) {
    const auto model = hpge_model();
    const auto edges = uniform_edges(400.0, 400);
    fs::ArrivalTableOptions options;
    options.energy_min_keV = 50.0;
    options.energy_max_keV = 400.0;
    options.exact_line_energies_keV = {300.0};
    const fs::ArrivalResponseTable table(model, edges, options);

    const auto& anchor = model.payload().detected_rate_surface.anchor;
    for (double energy : {300.0, 301.7}) {
        // An isotropic line, 2 photons per cm^2: counts = 2 x mean area x the
        // area-weighted response.
        fs::IncidentField<double> field;
        field.lines.push_back({energy, {2.0}});
        const auto folded = table.fold(field);
        const auto band = fs::arrival_band("iso");
        const double area = fs::arrival_quadrature(
            band, anchor.envelope_radius_cm, anchor.envelope_length_cm)
            .mean_projected_area_cm2;
        const auto expected = binned(
            fs::evaluate_arrival(model, energy, band, edges), edges, 2.0 * area);
        BOOST_TEST_MESSAGE(energy << " keV: integral " << sum(folded) / sum(expected) - 1
                           << ", TV " << total_variation(folded, expected));
        // Exact rows: only the mu step (the response linear between the
        // Gauss-Legendre nodes, the area exact). Interpolated rows: the two
        // neighbouring energy nodes, 2% apart, too.
        const bool exact = energy == 300.0;
        BOOST_TEST(std::abs(sum(folded) / sum(expected) - 1) < (exact ? 1e-5 : 1e-4));
        BOOST_TEST(total_variation(folded, expected) < (exact ? 1e-6 : 1e-3));
    }
}

BOOST_AUTO_TEST_CASE(arrival_table_folds_a_continuum_like_a_fine_direct_fold) {
    const auto model = hpge_model();
    const auto edges = uniform_edges(300.0, 150);
    fs::ArrivalTableOptions options;
    options.energy_min_keV = 150.0;
    options.energy_max_keV = 300.0;
    const fs::ArrivalResponseTable table(model, edges, options);
    const auto& anchor = model.payload().detected_rate_surface.anchor;

    // Three 20 keV cells x four mu bins, uneven fluences.
    fs::IncidentField<double> field;
    field.mu_edges = {-1.0, -0.5, 0.0, 0.5, 1.0};
    field.energy_edges_keV = {200.0, 220.0, 240.0, 260.0};
    field.continuum = {1.0, 0.5, 2.0, 3.0,  0.2, 0.7, 1.5, 0.4,  2.5, 1.0, 0.3, 0.8};
    const auto folded = table.fold(field);

    // Direct: 24 energies per cell (midpoints; the atoms that track the
    // incident energy spread over each one's slice), 8 Gauss-Legendre points
    // in mu per bin, each photon weighted by the envelope's projected area.
    std::vector<double> direct(edges.size() - 1, 0.0);
    std::vector<double> x(8), w(8);
    fs::detail::gauss_legendre_nodes(8, x.data(), w.data());
    for (std::size_t c = 0; c < 3; ++c) {
        for (int s = 0; s < 24; ++s) {
            const double energy = field.energy_edges_keV[c] + (s + 0.5) * 20.0 / 24;
            for (std::size_t m = 0; m < 4; ++m) {
                const double lo = field.mu_edges[m], hi = field.mu_edges[m + 1];
                for (int q = 0; q < 8; ++q) {
                    const double mu = 0.5 * (lo + hi) + 0.5 * (hi - lo) * x[q];
                    const double weight = field.continuum[c * 4 + m] / 24.0 * 0.5 * w[q]
                        * fs::envelope_projected_area_cm2(anchor.envelope_radius_cm,
                                                          anchor.envelope_length_cm, mu);
                    const auto part = binned(model.evaluate(energy, mu, edges), edges,
                                             weight, 20.0 / 24);
                    for (std::size_t b = 0; b < direct.size(); ++b)
                        direct[b] += part[b];
                }
            }
        }
    }
    BOOST_TEST_MESSAGE("continuum: integral " << sum(folded) / sum(direct) - 1
                       << ", TV " << total_variation(folded, direct));
    // The shape error is the continuum interpolated between energy nodes 2%
    // apart (the Compton edge), on unbroadened 2 keV bins.
    BOOST_TEST(std::abs(sum(folded) / sum(direct) - 1) < 1e-3);
    BOOST_TEST(total_variation(folded, direct) < 5e-3);
}

BOOST_AUTO_TEST_CASE(arrival_table_fold_is_linear_and_differentiable_in_the_field) {
    const auto model = hpge_model();
    const auto edges = uniform_edges(700.0, 140);
    fs::ArrivalTableOptions options;
    options.energy_min_keV = 100.0;
    options.energy_max_keV = 700.0;
    options.exact_line_energies_keV = {661.657};
    const fs::ArrivalResponseTable table(model, edges, options);

    fs::IncidentField<Dual> field;
    field.mu_edges = {-1.0, 0.0, 1.0};
    field.energy_edges_keV = {300.0, 400.0, 500.0};
    field.continuum = {Dual(1.0), Dual(2.0, 1.0), Dual(0.5), Dual(0.25)};
    field.lines.push_back({661.657, {Dual(3.0), Dual(1.0)}});
    field.lines.push_back({480.0, {Dual(0.5), Dual(0.5)}});
    const auto folded = table.fold(field);

    // The derivative with respect to one cell's fluence is that cell's fold.
    fs::IncidentField<double> unit;
    unit.mu_edges = field.mu_edges;
    unit.energy_edges_keV = field.energy_edges_keV;
    unit.continuum = {0.0, 1.0, 0.0, 0.0};
    const auto response = table.fold(unit);
    for (std::size_t b = 0; b < folded.size(); ++b)
        BOOST_TEST(folded[b].d == response[b], boost::test_tools::tolerance(1e-12));
}

BOOST_AUTO_TEST_CASE(arrival_table_rejects_bad_fields) {
    const auto model = hpge_model();
    fs::ArrivalTableOptions options;
    options.energy_min_keV = 50.0;
    options.energy_max_keV = 400.0;
    const fs::ArrivalResponseTable table(model, uniform_edges(400.0, 100), options);
    fs::IncidentField<double> outside;
    outside.lines.push_back({500.0, {1.0}});
    BOOST_CHECK_THROW(table.fold(outside), std::invalid_argument);
    fs::IncidentField<double> bad_mu;
    bad_mu.mu_edges = {-1.5, 1.0};
    BOOST_CHECK_THROW(table.fold(bad_mu), std::invalid_argument);
    fs::IncidentField<double> bad_size;
    bad_size.energy_edges_keV = {100.0, 200.0};
    bad_size.continuum = {1.0, 2.0};
    BOOST_CHECK_THROW(table.fold(bad_size), std::invalid_argument);
}
