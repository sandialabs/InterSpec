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

#define BOOST_TEST_MODULE BetaSpectrumTests
#include <boost/test/unit_test.hpp>

#include "cascade/BetaSpectrum.h"
#include "geometry/SourceGeometry.h"
#include "materials/Material.h"

#include <Eigen/Core>

#include <cmath>
#include <cstdint>
#include <random>

BOOST_AUTO_TEST_CASE(beta_density_normalizes_and_samples_within_support)
{
    const ceelo::BetaBranch branch{
        2290.0, 1.0, 92, 234, false, ceelo::BetaShape::Allowed};
    const ceelo::BetaSpectrumSampler sampler(branch);

    constexpr int steps = 20000;
    double integral = 0.0;
    for (int index = 0; index <= steps; ++index) {
        const double energy = branch.endpoint_keV * index / steps;
        const double weight =
            (index == 0 || index == steps) ? 0.5 : 1.0;
        integral += weight * sampler.normalized_density_per_keV(energy);
    }
    integral *= branch.endpoint_keV / steps;
    BOOST_TEST(integral == 1.0, boost::test_tools::tolerance(2e-3));

    std::mt19937_64 rng(0x4245544154455354ULL);
    double sample_sum = 0.0;
    for (int index = 0; index < 5000; ++index) {
        const double energy = sampler.sample_keV(rng);
        BOOST_TEST(energy >= 0.0);
        BOOST_TEST(energy <= branch.endpoint_keV);
        sample_sum += energy;
    }
    const double mean_fraction = sample_sum / 5000.0 / branch.endpoint_keV;
    BOOST_TEST(mean_fraction > 0.15);
    BOOST_TEST(mean_fraction < 0.65);
}

BOOST_AUTO_TEST_CASE(primary_electron_vacuum_escape_is_separate_from_photons)
{
    ceelo::SourceGeometry source;
    source.configure_point(Eigen::Vector3d::Zero());
    std::mt19937_64 rng(12345);
    const auto result = source.transport_source_electron(
        Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ(),
        1000.0, rng);
    BOOST_TEST(result.escaped);
    BOOST_TEST(result.photons.empty());
    BOOST_TEST(
        result.escaped_electron.energy_keV == 1000.0,
        boost::test_tools::tolerance(1e-14));
}

BOOST_AUTO_TEST_CASE(primary_electron_in_uranium_produces_valid_brems)
{
    ceelo::Material uranium(
        "uranium", 1.0, {{static_cast<std::uint8_t>(92), 1.0}});
    ceelo::SourceGeometry source;
    source.configure_point(Eigen::Vector3d::Zero());
    source.add_shield(&uranium, 1.0);
    source.set_source_brems(true);
    source.set_record_truth_labels(true);

    std::mt19937_64 rng(0x4252454d53544553ULL);
    std::size_t photons = 0;
    std::size_t created = 0;
    for (int event = 0; event < 1000; ++event) {
        const auto result = source.transport_source_electron(
            Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ(),
            2290.0, rng);
        for (const auto& photon : result.created_brems_photons) {
            BOOST_TEST(std::isfinite(photon.energy_keV));
            BOOST_TEST(photon.energy_keV >= 10.0);
            BOOST_TEST(photon.energy_keV <= 2290.0);
            BOOST_TEST(photon.position.norm() < 1.00001);
            BOOST_TEST(
                photon.direction.norm() == 1.0,
                boost::test_tools::tolerance(1e-12));
            ++created;
        }
        for (const auto& photon : result.photons) {
            BOOST_TEST(std::isfinite(photon.energy_keV));
            BOOST_TEST(photon.energy_keV >= 10.0);
            BOOST_TEST(photon.energy_keV <= 2290.0);
            ++photons;
        }
    }
    BOOST_TEST(created > 0u);
    BOOST_TEST(photons > 0u);
}

BOOST_AUTO_TEST_CASE(low_energy_conversion_electron_uses_extended_sb_grid)
{
    ceelo::Material uranium(
        "uranium", 1.0, {{static_cast<std::uint8_t>(92), 1.0}});
    ceelo::SourceGeometry source;
    source.configure_point(Eigen::Vector3d::Zero());
    source.add_shield(&uranium, 1.0);
    source.set_source_brems(true);
    source.set_record_truth_labels(true);

    std::mt19937_64 rng(0x4c4f57454e455247ULL);
    std::size_t created = 0;
    for (int event = 0; event < 5000; ++event) {
        const auto result = source.transport_source_electron(
            Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ(),
            40.0, rng);
        for (const auto& photon : result.created_brems_photons) {
            BOOST_TEST(photon.energy_keV >= 10.0);
            BOOST_TEST(photon.energy_keV <= 40.0);
            ++created;
        }
    }
    BOOST_TEST(created > 0u);
}
