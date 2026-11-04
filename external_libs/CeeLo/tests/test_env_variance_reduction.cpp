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

#define BOOST_TEST_MODULE EnvVarianceReductionTest
#include <boost/test/unit_test.hpp>

#include "geometry/EnvironmentGeometry.h"
#include "materials/Material.h"
#include "transport/NextEventEstimator.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <future>
#include <numeric>
#include <vector>

using namespace ceelo;

namespace {

double sum(const std::vector<double>& values)
{
    return std::accumulate(values.begin(), values.end(), 0.0);
}

void check_exact_closure(const FluenceTally& tally)
{
    const double total = tally.total_sum_w();
    const double scale = std::max(1.0, std::abs(total));
    const double common = sum(tally.continuum_sum_w()) +
                          sum(tally.line_total_sum_w());
    const double components = sum(tally.component_total_sum_w());
    const double line_components = sum(tally.line_component_sum_w());
    const double lines = sum(tally.line_total_sum_w());

    BOOST_TEST(common == total,
               boost::test_tools::tolerance(2e-13 * scale));
    BOOST_TEST(components == total,
               boost::test_tools::tolerance(2e-13 * scale));
    BOOST_TEST(line_components == lines,
               boost::test_tools::tolerance(2e-13 * scale));

    double origins = 0.0;
    for (int origin = 0; origin < static_cast<int>(FluxOrigin::Count);
         ++origin)
        origins += tally.total_mean_origin(static_cast<FluxOrigin>(origin));
    BOOST_TEST(origins == tally.total_mean(),
               boost::test_tools::tolerance(2e-13));
}

EnvScatterConfig test_config(double e0_keV)
{
    EnvScatterConfig cfg;
    cfg.e_floor_keV = 10.0;
    cfg.axes = make_registered_landmark_axes(e0_keV, 24, 8, 10.0);
    cfg.scoring_radius_cm = 4.0;
    cfg.max_vertices = 1;
    cfg.include_annihilation = false;
    cfg.include_brems = false;
    return cfg;
}

void check_consistent_means(const EnvScatterResult& standard,
                            const EnvScatterResult& candidate)
{
    const double a = standard.tally.total_mean();
    const double b = candidate.tally.total_mean();
    const double se = std::hypot(standard.tally.total_sigma(),
                                 candidate.tally.total_sigma());
    BOOST_REQUIRE_GT(a, 0.0);
    BOOST_REQUIRE_GT(b, 0.0);
    BOOST_REQUIRE_GT(se, 0.0);
    // This is a guard against missing 1/2 weights, not a performance
    // promotion gate.  The standalone benchmark performs paired-replica bias
    // tests and deliberately reports inconclusive when it lacks resolution.
    BOOST_TEST(std::abs(a - b) < 5.0 * se);
}

} // namespace

BOOST_AUTO_TEST_CASE(experimental_policies_are_off_by_default)
{
    const EnvScatterConfig cfg;
    BOOST_CHECK(cfg.experimental_primary_sampling ==
                ExperimentalPrimarySampling::Analog);
    BOOST_CHECK(cfg.experimental_cone_sampling ==
                ExperimentalConeSampling::Single);
}

BOOST_AUTO_TEST_CASE(explicit_off_policies_preserve_the_default_rng_stream)
{
    constexpr double energy = 122.0;
    const Material air = make_Air();
    const EnvironmentGeometry environment(air);
    const Eigen::Vector3d source(0.0, 0.0, -100.0);
    EnvScatterConfig implicit = test_config(energy);
    EnvScatterConfig explicit_off = implicit;
    explicit_off.experimental_primary_sampling =
        ExperimentalPrimarySampling::Analog;
    explicit_off.experimental_cone_sampling =
        ExperimentalConeSampling::Single;

    EnvScatterResult a(implicit.axes, energy);
    EnvScatterResult b(explicit_off.axes, energy);
    NextEventEstimator(environment, implicit, source)
        .run(energy, 128, 0xd3fa017ULL, a);
    NextEventEstimator(environment, explicit_off, source)
        .run(energy, 128, 0xd3fa017ULL, b);

    BOOST_TEST(a.tally.continuum_sum_w() == b.tally.continuum_sum_w(),
               boost::test_tools::per_element());
    BOOST_TEST(a.tally.continuum_sum_w2() == b.tally.continuum_sum_w2(),
               boost::test_tools::per_element());
    BOOST_TEST(a.tally.line_sum_w() == b.tally.line_sum_w(),
               boost::test_tools::per_element());
    BOOST_TEST(a.diagnostics.interaction_counts ==
                   b.diagnostics.interaction_counts,
               boost::test_tools::per_element());
}

BOOST_AUTO_TEST_CASE(antipodal_pair_is_one_half_weight_statistical_history)
{
    constexpr double energy = 122.0;
    constexpr std::uint64_t histories = 8192;
    const Material air = make_Air();
    const EnvironmentGeometry environment(air);
    const Eigen::Vector3d source(0.0, 0.0, -100.0);

    EnvScatterConfig standard_cfg = test_config(energy);
    EnvScatterConfig paired_cfg = standard_cfg;
    paired_cfg.experimental_primary_sampling =
        ExperimentalPrimarySampling::AntipodalPair;

    EnvScatterResult standard(standard_cfg.axes, energy);
    EnvScatterResult paired(paired_cfg.axes, energy);
    NextEventEstimator(environment, standard_cfg, source)
        .run(energy, histories, 0x51a7ULL, standard, 3);
    NextEventEstimator(environment, paired_cfg, source)
        .run(energy, histories, 0x51a7ULL, paired, 3);

    BOOST_TEST(paired.histories == histories);
    BOOST_TEST(paired.tally.num_histories() == histories);
    BOOST_TEST(paired.tally.replica_histories()[3] == histories);

    // With one vertex in unbounded air, every transported primary undergoes
    // exactly one analog event.  The paired policy therefore transports two
    // half-weight branches while closing only one statistical history.
    const auto interaction_total = [](const EnvScatterResult& result) {
        return std::accumulate(result.diagnostics.interaction_counts.begin(),
                               result.diagnostics.interaction_counts.end(),
                               std::uint64_t{0});
    };
    BOOST_TEST(interaction_total(standard) == histories);
    BOOST_TEST(interaction_total(paired) == 2 * histories);
    check_consistent_means(standard, paired);
    check_exact_closure(paired.tally);
}

BOOST_AUTO_TEST_CASE(detector_cone_antithetic_pair_is_unbiased_and_closes)
{
    constexpr double energy = 122.0;
    constexpr std::uint64_t histories = 8192;
    const Material air = make_Air();
    const Material soil = make_Soil();
    const Eigen::Vector3d source(0.0, 0.0, -40.0);
    const EnvironmentGeometry environment =
        EnvironmentGeometry::field(air, soil, source, 25.0, 25.0, 60.0);

    EnvScatterConfig standard_cfg = test_config(energy);
    EnvScatterConfig antithetic_cfg = standard_cfg;
    antithetic_cfg.experimental_cone_sampling =
        ExperimentalConeSampling::AntitheticPair;

    EnvScatterResult standard(standard_cfg.axes, energy);
    EnvScatterResult antithetic(antithetic_cfg.axes, energy);
    NextEventEstimator(environment, standard_cfg, source)
        .run(energy, histories, 0xc0aeULL, standard);
    NextEventEstimator(environment, antithetic_cfg, source)
        .run(energy, histories, 0xc0aeULL, antithetic);

    BOOST_TEST(antithetic.histories == histories);
    BOOST_TEST(antithetic.tally.num_histories() == histories);
    check_consistent_means(standard, antithetic);
    check_exact_closure(antithetic.tally);
}

BOOST_AUTO_TEST_CASE(experimental_streams_are_thread_schedule_independent)
{
    constexpr double energy = 122.0;
    constexpr std::uint64_t histories_per_stream = 32;
    const Material air = make_Air();
    const EnvironmentGeometry environment(air);
    const Eigen::Vector3d source(0.0, 0.0, -100.0);
    EnvScatterConfig cfg = test_config(energy);
    cfg.experimental_primary_sampling =
        ExperimentalPrimarySampling::AntipodalPair;
    cfg.experimental_cone_sampling =
        ExperimentalConeSampling::AntitheticPair;
    const std::vector<std::uint64_t> seeds{
        0x101ULL, 0x202ULL, 0x303ULL, 0x404ULL};

    auto run_stream = [&](std::uint64_t seed) {
        EnvScatterResult result(cfg.axes, energy);
        NextEventEstimator(environment, cfg, source)
            .run(energy, histories_per_stream, seed, result);
        return result;
    };

    std::vector<EnvScatterResult> serial;
    for (std::uint64_t seed : seeds) serial.push_back(run_stream(seed));

    std::vector<std::future<EnvScatterResult>> futures;
    for (std::uint64_t seed : seeds)
        futures.push_back(std::async(std::launch::async, run_stream, seed));
    std::vector<EnvScatterResult> parallel;
    for (auto& future : futures) parallel.push_back(future.get());

    FluenceTally serial_total(cfg.axes, env_line_channels(energy));
    FluenceTally parallel_total(cfg.axes, env_line_channels(energy));
    for (std::size_t i = 0; i < seeds.size(); ++i) {
        serial_total.merge(serial[i].tally);
        parallel_total.merge(parallel[i].tally);
    }

    BOOST_TEST(serial_total.num_histories() ==
               histories_per_stream * seeds.size());
    BOOST_TEST(serial_total.continuum_sum_w() ==
                   parallel_total.continuum_sum_w(),
               boost::test_tools::per_element());
    BOOST_TEST(serial_total.continuum_sum_w2() ==
                   parallel_total.continuum_sum_w2(),
               boost::test_tools::per_element());
    BOOST_TEST(serial_total.component_total_sum_w() ==
                   parallel_total.component_total_sum_w(),
               boost::test_tools::per_element());
    BOOST_TEST(serial_total.line_component_sum_w() ==
                   parallel_total.line_component_sum_w(),
               boost::test_tools::per_element());
    BOOST_TEST(serial_total.total_sum_w() == parallel_total.total_sum_w());
    BOOST_TEST(serial_total.total_sum_w2() == parallel_total.total_sum_w2());
    check_exact_closure(serial_total);
}

BOOST_AUTO_TEST_CASE(analog_truth_rejects_next_event_only_experiments)
{
    constexpr double energy = 122.0;
    const Material air = make_Air();
    const EnvironmentGeometry environment(air);
    const Eigen::Vector3d source(0.0, 0.0, -100.0);
    EnvScatterConfig cfg = test_config(energy);
    cfg.experimental_primary_sampling =
        ExperimentalPrimarySampling::AntipodalPair;
    EnvScatterResult result(cfg.axes, energy);
    BOOST_CHECK_THROW(
        NextEventEstimator(environment, cfg, source)
            .run_brute(energy, 1, 1, result),
        std::invalid_argument);
}
