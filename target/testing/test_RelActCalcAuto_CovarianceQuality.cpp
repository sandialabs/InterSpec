/* Focused tests for quantity-specific RelActAuto covariance diagnostics. */

#include "InterSpec_config.h"

#include <vector>

#define BOOST_TEST_MODULE RelActCalcAuto_CovarianceQuality_suite
#include <boost/test/included/unit_test.hpp>

#include "InterSpec/RelActCalcAuto_CovarianceQuality.h"

using RelActCalcAutoImp::QuantityCovarianceDiagnostic;
using RelActCalcAutoImp::assess_quantity_covariance;

BOOST_AUTO_TEST_CASE( unrelated_dropped_direction_does_not_poison_quantity )
{
  const std::vector<std::vector<double>> basis{
    {1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,1.0}
  };
  const std::vector<std::vector<double>> dropped{{0.0,0.0,1.0}};

  const QuantityCovarianceDiagnostic unrelated
      = assess_quantity_covariance({2.0,-1.0,0.0},basis,dropped,{});
  BOOST_REQUIRE( unrelated.dimensions_valid );
  BOOST_REQUIRE( unrelated.finite );
  BOOST_CHECK( !unrelated.depends_on_dropped_direction );
  BOOST_CHECK_SMALL( unrelated.dropped_fraction, 1.0e-15 );

  const QuantityCovarianceDiagnostic dependent
      = assess_quantity_covariance({2.0,-1.0,0.2},basis,dropped,{});
  BOOST_REQUIRE( dependent.finite );
  BOOST_CHECK( dependent.depends_on_dropped_direction );
  BOOST_CHECK_GT( dependent.dropped_fraction, 0.005 );
}

BOOST_AUTO_TEST_CASE( active_bound_is_quantity_specific_and_scale_aware )
{
  // Parameter zero has a tiny raw derivative but a large natural step; parameter
  // one has the converse.  In equilibrated coordinates they contribute equally.
  const std::vector<std::vector<double>> basis{{100.0,0.0},{0.0,0.01}};
  const std::vector<char> first_at_bound{1,0};

  const QuantityCovarianceDiagnostic dependent
      = assess_quantity_covariance({1.0e-4,1.0},basis,{},first_at_bound);
  BOOST_REQUIRE( dependent.finite );
  BOOST_CHECK( dependent.depends_on_active_bound );
  BOOST_CHECK_CLOSE_FRACTION( dependent.active_bound_fraction, 0.5, 1.0e-12 );

  const QuantityCovarianceDiagnostic unrelated
      = assess_quantity_covariance({0.0,1.0},basis,{},first_at_bound);
  BOOST_REQUIRE( unrelated.finite );
  BOOST_CHECK( !unrelated.depends_on_active_bound );
  BOOST_CHECK_SMALL( unrelated.active_bound_fraction, 1.0e-15 );
}

BOOST_AUTO_TEST_CASE( numerical_leakage_is_not_material )
{
  const std::vector<std::vector<double>> basis{{1.0,0.0},{0.0,1.0}};
  const std::vector<std::vector<double>> dropped{{0.0,1.0}};
  const QuantityCovarianceDiagnostic diagnostic
      = assess_quantity_covariance({1.0,1.0e-4},basis,dropped,{});
  BOOST_REQUIRE( diagnostic.finite );
  BOOST_CHECK( !diagnostic.depends_on_dropped_direction );
  BOOST_CHECK_LT( diagnostic.dropped_fraction, 1.0e-6 );
}
