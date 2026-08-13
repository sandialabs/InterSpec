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

#include <cmath>
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
