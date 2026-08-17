/* InterSpec: focused fitted-age gamma-yield derivative tests. */

#include "InterSpec_config.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <future>
#include <limits>
#include <map>
#include <vector>

#define BOOST_TEST_MODULE RelActCalcAuto_AgeDerivative_suite
#include <boost/test/included/unit_test.hpp>

#include "ceres/jet.h"

#include "InterSpec/RelActCalcAuto_AgeDerivative.h"

using namespace RelActCalcAutoImp;

namespace
{

std::vector<AgeGammaYield> smooth_sample( const double age )
{
  // Deliberately change order on every side of 6 seconds.  The production derivative must match
  // stable energy identity rather than assuming that decay output position or transition metadata
  // is invariant with age.
  std::vector<AgeGammaYield> answer{
    { age_gamma_energy_identity(100.0), std::exp(-age/7.0) },
    { age_gamma_energy_identity(200.0), 0.5 + 0.2*age + 0.03*age*age },
    { age_gamma_energy_identity(300.0), std::cos(age/9.0) }
  };
  if( age > 6.0 )
    std::reverse( answer.begin(), answer.end() );
  return answer;
}

template<int NumLanes>
void check_jet_lanes_at_boundary( const bool upper_boundary )
{
  const double age = upper_boundary ? 10.0 : 0.0;
  const auto sample = []( const double t ) {
    return std::vector<AgeGammaYield>{
      { age_gamma_energy_identity(999.0), std::exp(-t/4.0) }
    };
  };
  const AgeDerivativeResult derivative = richardson_age_yield_derivative(
      age, 0.0, 10.0, 4.0, sample(age), sample );
  BOOST_REQUIRE( derivative.status == AgeDerivativeStatus::Converged );
  BOOST_REQUIRE_EQUAL( derivative.derivatives.size(), 1U );
  BOOST_CHECK( derivative.stencil == (upper_boundary ? AgeDerivativeStencil::Backward
                                                     : AgeDerivativeStencil::Forward) );

  ceres::Jet<double,NumLanes> jet_age(age);
  ceres::Jet<double,NumLanes> jet_yield(std::exp(-age/4.0));
  for( int lane = 0; lane < NumLanes; ++lane )
    jet_age.v[lane] = (lane % 2) ? -0.25*(lane + 1.0) : 0.5*(lane + 1.0);

  attach_age_yield_derivative( jet_yield, derivative.derivatives[0], jet_age );
  const double exact_derivative = -std::exp(-age/4.0)/4.0;
  BOOST_CHECK_CLOSE_FRACTION( derivative.derivatives[0], exact_derivative, 2.0e-7 );
  for( int lane = 0; lane < NumLanes; ++lane )
    BOOST_CHECK_EQUAL( jet_yield.v[lane], derivative.derivatives[0]*jet_age.v[lane] );
}

} // namespace


BOOST_AUTO_TEST_CASE( centered_richardson_converges_and_matches_by_identity )
{
  const double age = 6.0;
  // Reference order is neither of the sampler's orders, exercising output realignment too.
  const std::vector<AgeGammaYield> reference{
    { age_gamma_energy_identity(300.0), std::cos(age/9.0) },
    { age_gamma_energy_identity(100.0), std::exp(-age/7.0) },
    { age_gamma_energy_identity(200.0), 0.5 + 0.2*age + 0.03*age*age }
  };

  const AgeDerivativeResult result = richardson_age_yield_derivative(
      age, 0.0, 20.0, 7.0, reference, smooth_sample );

  BOOST_REQUIRE( result.status == AgeDerivativeStatus::Converged );
  BOOST_CHECK( result.stencil == AgeDerivativeStencil::Centered );
  BOOST_REQUIRE_EQUAL( result.derivatives.size(), 3U );
  BOOST_CHECK_CLOSE_FRACTION( result.derivatives[0], -std::sin(age/9.0)/9.0, 2.0e-8 );
  BOOST_CHECK_CLOSE_FRACTION( result.derivatives[1], -std::exp(-age/7.0)/7.0, 2.0e-8 );
  BOOST_CHECK_CLOSE_FRACTION( result.derivatives[2], 0.2 + 0.06*age, 2.0e-10 );
  BOOST_CHECK_GE( result.refinements, 3U );
}


BOOST_AUTO_TEST_CASE( physical_lower_bound_uses_second_order_forward_ladder )
{
  const double age = 0.0;
  const auto sample = []( const double t ) {
    return std::vector<AgeGammaYield>{
      { age_gamma_energy_identity(111.0), std::exp(-t/4.0) },
      { age_gamma_energy_identity(222.0), 2.0 + 3.0*t + 0.5*t*t }
    };
  };
  const std::vector<AgeGammaYield> reference = sample(age);

  const AgeDerivativeResult result = richardson_age_yield_derivative(
      age, 0.0, 12.0, 4.0, reference, sample );

  BOOST_REQUIRE( result.status == AgeDerivativeStatus::Converged );
  BOOST_CHECK( result.stencil == AgeDerivativeStencil::Forward );
  BOOST_CHECK_CLOSE_FRACTION( result.derivatives[0], -0.25, 5.0e-8 );
  BOOST_CHECK_CLOSE_FRACTION( result.derivatives[1], 3.0, 2.0e-11 );
  BOOST_CHECK_GT( result.initial_step, 0.0 );
}


BOOST_AUTO_TEST_CASE( input_upper_bound_uses_second_order_backward_ladder )
{
  const double age = 10.0;
  const auto sample = []( const double t ) {
    return std::vector<AgeGammaYield>{
      { age_gamma_energy_identity(444.0), std::exp(-t/3.0) },
      { age_gamma_energy_identity(333.0), 1.0 - 0.1*t + 0.02*t*t }
    };
  };

  const AgeDerivativeResult result = richardson_age_yield_derivative(
      age, 2.0, 10.0, 3.0, sample(age), sample );

  BOOST_REQUIRE( result.status == AgeDerivativeStatus::Converged );
  BOOST_CHECK( result.stencil == AgeDerivativeStencil::Backward );
  BOOST_CHECK_CLOSE_FRACTION( result.derivatives[0], -std::exp(-age/3.0)/3.0, 8.0e-8 );
  BOOST_CHECK_CLOSE_FRACTION( result.derivatives[1], -0.1 + 0.04*age, 2.0e-10 );
}


BOOST_AUTO_TEST_CASE( near_bound_prefers_resolvable_inward_step )
{
  const double age = 1.0e-10;
  double minimum_sampled_age = age;
  const auto sample = [&minimum_sampled_age]( const double t ) {
    minimum_sampled_age = (std::min)(minimum_sampled_age, t);
    return std::vector<AgeGammaYield>{
      { age_gamma_energy_identity(555.0), std::exp(-t/5.0) }
    };
  };
  const std::vector<AgeGammaYield> reference{
    { age_gamma_energy_identity(555.0), std::exp(-age/5.0) }
  };

  const AgeDerivativeResult result = richardson_age_yield_derivative(
      age, 0.0, 20.0, 5.0, reference, sample );

  BOOST_REQUIRE( result.status == AgeDerivativeStatus::Converged );
  BOOST_CHECK( result.stencil == AgeDerivativeStencil::Forward );
  BOOST_CHECK_GE( minimum_sampled_age, 0.0 );
  BOOST_CHECK_CLOSE_FRACTION( result.derivatives[0], -std::exp(-age/5.0)/5.0, 5.0e-8 );
}


BOOST_AUTO_TEST_CASE( scaled_parameter_one_ulp_bound_roundoff_is_accepted )
{
  const double input_lower_bound = 3.0;
  const double evaluated_age = std::nextafter(input_lower_bound,0.0);
  const auto sample = []( const double t ) {
    return std::vector<AgeGammaYield>{
      { age_gamma_energy_identity(556.0), std::exp(-t/5.0) }
    };
  };
  const AgeDerivativeResult result = richardson_age_yield_derivative(
      evaluated_age,input_lower_bound,20.0,5.0,sample(evaluated_age),sample );
  BOOST_REQUIRE( result.status == AgeDerivativeStatus::Converged );
  BOOST_CHECK( result.stencil == AgeDerivativeStencil::Forward );
  BOOST_CHECK_CLOSE_FRACTION( result.derivatives[0], -std::exp(-evaluated_age/5.0)/5.0,
                              2.0e-7 );
}


BOOST_AUTO_TEST_CASE( derivatives_do_not_jump_at_legacy_half_life_thresholds )
{
  const double half_life = 8.0;
  const auto sample = [half_life]( const double t ) {
    return std::vector<AgeGammaYield>{
      { age_gamma_energy_identity(666.0), std::exp(-t/half_life) }
    };
  };

  for( const double ratio : { 0.1, 2.5, 5.0 } )
  {
    const double epsilon = 1.0e-7*half_life;
    const double below = ratio*half_life - epsilon;
    const double above = ratio*half_life + epsilon;
    const AgeDerivativeResult left = richardson_age_yield_derivative(
        below, 0.0, 80.0, half_life, sample(below), sample );
    const AgeDerivativeResult right = richardson_age_yield_derivative(
        above, 0.0, 80.0, half_life, sample(above), sample );
    BOOST_REQUIRE( left.status == AgeDerivativeStatus::Converged );
    BOOST_REQUIRE( right.status == AgeDerivativeStatus::Converged );
    BOOST_CHECK_CLOSE_FRACTION( left.derivatives[0], -std::exp(-below/half_life)/half_life,
                                5.0e-8 );
    BOOST_CHECK_CLOSE_FRACTION( right.derivatives[0], -std::exp(-above/half_life)/half_life,
                                5.0e-8 );
    BOOST_CHECK_CLOSE_FRACTION( left.derivatives[0], right.derivatives[0], 3.0e-7 );
  }
}


BOOST_AUTO_TEST_CASE( nonconvergence_and_identity_failure_are_explicit )
{
  AgeDerivativeOptions strict;
  strict.relative_tolerance = 0.0;
  strict.absolute_scaled_tolerance = 0.0;
  strict.minimum_refinements = 3;
  strict.maximum_refinements = 4;

  const double age = 1.0;
  const std::vector<AgeGammaYield> discontinuous_reference{
    { age_gamma_energy_identity(777.0), 1.0 }
  };
  const auto discontinuous = []( const double t ) {
    return std::vector<AgeGammaYield>{
      { age_gamma_energy_identity(777.0), (t >= 1.0) ? 1.0 : 0.0 }
    };
  };
  const AgeDerivativeResult best_effort = richardson_age_yield_derivative(
      age, 0.0, 2.0, 1.0, discontinuous_reference, discontinuous, strict );
  BOOST_CHECK( best_effort.status == AgeDerivativeStatus::BestFiniteEstimate );
  BOOST_CHECK_EQUAL( best_effort.unconverged_lines, 1U );
  BOOST_CHECK( std::isfinite(best_effort.derivatives[0]) );
  BOOST_CHECK( best_effort.diagnostic.find("did not satisfy") != std::string::npos );

  const auto missing_identity = []( const double ) {
    return std::vector<AgeGammaYield>{};
  };
  const AgeDerivativeResult failed = richardson_age_yield_derivative(
      age, 0.0, 2.0, 1.0, discontinuous_reference, missing_identity );
  BOOST_CHECK( failed.status == AgeDerivativeStatus::Failed );
  BOOST_CHECK( failed.failure == AgeDerivativeFailure::IdentityMismatch );
  BOOST_CHECK( failed.diagnostic.find("identit") != std::string::npos );
}


BOOST_AUTO_TEST_CASE( duplicate_reference_identity_is_rejected )
{
  const std::vector<AgeGammaYield> reference{
    { age_gamma_energy_identity(888.0), 1.0 },
    { age_gamma_energy_identity(888.0), 2.0 }
  };
  const AgeDerivativeResult result = richardson_age_yield_derivative(
      1.0, 0.0, 2.0, 1.0, reference,
      []( const double ){ return std::vector<AgeGammaYield>{}; } );
  BOOST_CHECK( result.status == AgeDerivativeStatus::Failed );
  BOOST_CHECK( result.failure == AgeDerivativeFailure::DuplicateIdentity );
}


BOOST_AUTO_TEST_CASE( production_cache_key_isolates_age_exclusions_bounds_and_source )
{
  int first_source = 0;
  int second_source = 0;
  const double age = 10.0;
  const AgeDerivativeCacheKey base{ &first_source, age, {100.0,200.0,100.0},
                                    0.0, 20.0, 5.0 };
  const AgeDerivativeCacheKey equivalent{ &first_source, age, {200.0,100.0},
                                          0.0, 20.0, 5.0 };
  const AgeDerivativeCacheKey adjacent_age{
      &first_source, std::nextafter(age,std::numeric_limits<double>::infinity()), {100.0,200.0},
      0.0, 20.0, 5.0 };
  const AgeDerivativeCacheKey different_exclusion{ &first_source, age, {100.0},
                                                   0.0, 20.0, 5.0 };
  const AgeDerivativeCacheKey different_lower_bound{ &first_source, age, {100.0,200.0},
                                                     1.0, 20.0, 5.0 };
  const AgeDerivativeCacheKey different_upper_bound{ &first_source, age, {100.0,200.0},
                                                     0.0, 21.0, 5.0 };
  const AgeDerivativeCacheKey different_time_scale{ &first_source, age, {100.0,200.0},
                                                    0.0, 20.0, 6.0 };
  const AgeDerivativeCacheKey different_source{ &second_source, age, {100.0,200.0},
                                                0.0, 20.0, 5.0 };

  std::map<AgeDerivativeCacheKey,int,AgeDerivativeCacheKeyLess> cache;
  cache[base] = 1;
  cache[equivalent] = 2;
  BOOST_REQUIRE_EQUAL( cache.size(), 1U );
  BOOST_CHECK_EQUAL( cache.begin()->second, 2 );

  cache[adjacent_age] = 3;
  cache[different_exclusion] = 4;
  cache[different_lower_bound] = 5;
  cache[different_upper_bound] = 6;
  cache[different_time_scale] = 7;
  cache[different_source] = 8;
  BOOST_CHECK_EQUAL( cache.size(), 7U );
}


BOOST_AUTO_TEST_CASE( scalar_derivative_attaches_identically_across_dynamic_autodiff_boundaries )
{
  check_jet_lanes_at_boundary<31>(false);
  check_jet_lanes_at_boundary<32>(false);
  check_jet_lanes_at_boundary<33>(false);
  check_jet_lanes_at_boundary<63>(true);
  check_jet_lanes_at_boundary<64>(true);
  check_jet_lanes_at_boundary<65>(true);
}


BOOST_AUTO_TEST_CASE( richardson_results_are_thread_deterministic )
{
  const double age = 6.0;
  const std::vector<AgeGammaYield> reference = smooth_sample(age);
  const auto calculate = [age,reference]() {
    return richardson_age_yield_derivative(age,0.0,20.0,7.0,reference,smooth_sample);
  };

  std::vector<std::future<AgeDerivativeResult>> work;
  for( int i = 0; i < 8; ++i )
    work.push_back( std::async(std::launch::async,calculate) );
  const AgeDerivativeResult expected = work.front().get();
  BOOST_REQUIRE( expected.status == AgeDerivativeStatus::Converged );
  for( std::size_t i = 1; i < work.size(); ++i )
  {
    const AgeDerivativeResult actual = work[i].get();
    BOOST_CHECK( actual.status == expected.status );
    BOOST_CHECK( actual.stencil == expected.stencil );
    BOOST_CHECK_EQUAL( actual.refinements, expected.refinements );
    BOOST_CHECK_EQUAL_COLLECTIONS( actual.derivatives.begin(), actual.derivatives.end(),
                                   expected.derivatives.begin(), expected.derivatives.end() );
  }
}
