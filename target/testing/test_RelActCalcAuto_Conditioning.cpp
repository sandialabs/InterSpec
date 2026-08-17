/* Focused tests for RelActAuto empirical fixed-pivot/QR conditioning. */
#include "InterSpec_config.h"

#include <cmath>
#include <vector>

#define BOOST_TEST_MODULE RelActCalcAuto_Conditioning_suite
#include <boost/test/included/unit_test.hpp>

#include "ceres/jet.h"

#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/RelActCalcAuto_Conditioning.hpp"
#include "InterSpec/RelActCalc_imp.hpp"

using namespace boost::unit_test;
using RelActCalcAutoImp::EmpiricalBasisTransform;

namespace
{
std::vector<double> frozen_energies()
{
  std::vector<double> answer;
  for( size_t i = 0; i < 96; ++i )
    answer.push_back( 610.0 + (775.0 - 610.0)*(static_cast<double>(i) + 0.5)/96.0 );
  return answer;
}

double evaluate( const RelActCalc::RelEffEqnForm form,
                 const std::vector<double> &coefficients,
                 const double energy )
{
  return RelActCalc::eval_eqn_imp( energy, form, coefficients.data(), coefficients.size() );
}

void check_exact_seed_mapping( const RelActCalc::RelEffEqnForm form,
                               const std::vector<double> &legacy )
{
  const double pivot = std::sqrt(610.0*775.0);
  const EmpiricalBasisTransform transform( form, legacy.size(), pivot, frozen_energies() );
  const EmpiricalBasisTransform::Seed seed = transform.orthogonalize_legacy_seed( legacy );
  const std::vector<double> round_trip = transform.legacy_coefficients(
                                            seed.orthogonal_coefficients );

  BOOST_REQUIRE_EQUAL( round_trip.size(), legacy.size() );
  for( size_t i = 0; i < legacy.size(); ++i )
    BOOST_CHECK_CLOSE_FRACTION( round_trip[i], seed.normalized_legacy_coefficients[i], 2.0e-11 );

  BOOST_CHECK_CLOSE_FRACTION( evaluate(form,round_trip,pivot), 1.0, 2.0e-12 );
  for( const double energy : {610.0, 641.3, 700.0, 742.8, 775.0} )
  {
    const double predictor = transform.linear_predictor(
      energy, seed.orthogonal_coefficients.data(), seed.orthogonal_coefficients.size() );
    const double direct = (form == RelActCalc::RelEffEqnForm::LnX)
                            ? predictor : std::exp(predictor);
    BOOST_CHECK_CLOSE_FRACTION( direct, evaluate(form,round_trip,energy), 3.0e-12 );

    const double old_objective_amplitude = 3.75 * evaluate(form,legacy,energy);
    const double new_objective_amplitude
      = (3.75*seed.activity_scale) * evaluate(form,round_trip,energy);
    BOOST_CHECK_CLOSE_FRACTION( new_objective_amplitude, old_objective_amplitude, 4.0e-11 );
  }
}

void check_orthonormal_design( const RelActCalc::RelEffEqnForm form,
                               const size_t num_coefficients )
{
  const std::vector<double> energies = frozen_energies();
  const double pivot = std::sqrt(610.0*775.0);
  const EmpiricalBasisTransform transform( form, num_coefficients, pivot, energies );
  const size_t num_shape = num_coefficients - 1;

  std::vector<std::vector<double>> columns( num_shape,
                                            std::vector<double>(energies.size(),0.0) );
  for( size_t col = 0; col < num_shape; ++col )
  {
    std::vector<double> orthogonal( num_coefficients, 0.0 );
    orthogonal[col + 1] = 1.0;
    const std::vector<double> legacy = transform.legacy_coefficients( orthogonal );
    for( size_t row = 0; row < energies.size(); ++row )
    {
      // Derivative of the pre-exponential linear predictor (and of LnX itself).
      double value = 0.0;
      for( size_t j = 0; j < legacy.size(); ++j )
        value += legacy[j] * RelActCalcAutoImp::empirical_legacy_basis_value(
                                      form, j, energies[row] );
      if( form == RelActCalc::RelEffEqnForm::LnX )
        value -= 1.0; //remove the affine pivot-normalization baseline; retain d(f)/dq
      columns[col][row] = value;
    }
  }

  for( size_t lhs = 0; lhs < num_shape; ++lhs )
  {
    for( size_t rhs = 0; rhs < num_shape; ++rhs )
    {
      double dot = 0.0;
      for( size_t row = 0; row < energies.size(); ++row )
        dot += columns[lhs][row]*columns[rhs][row];
      // Mapping back through a high-order legacy log-power basis loses a few ulps to
      // cancellation; the columns remain orthonormal to substantially better than 1e-7.
      BOOST_CHECK_SMALL( dot - ((lhs == rhs) ? 1.0 : 0.0), 2.0e-8 );
    }
  }
}
}//namespace


BOOST_AUTO_TEST_CASE( lny_fixed_pivot_round_trip_preserves_objective )
{
  check_exact_seed_mapping( RelActCalc::RelEffEqnForm::LnY,
                            {-4.2, 0.0065, 310.0} );
  check_orthonormal_design( RelActCalc::RelEffEqnForm::LnY, 3 );
}


BOOST_AUTO_TEST_CASE( lnxlny_fixed_pivot_round_trip_preserves_objective )
{
  check_exact_seed_mapping( RelActCalc::RelEffEqnForm::LnXLnY,
                            {-2.0, 0.55, -0.035, 0.0021, -0.00008} );
  check_orthonormal_design( RelActCalc::RelEffEqnForm::LnXLnY, 5 );
}


BOOST_AUTO_TEST_CASE( lnx_affine_normalization_preserves_objective )
{
  check_exact_seed_mapping( RelActCalc::RelEffEqnForm::LnX,
                            {0.7, 0.11, -0.006} );
  check_orthonormal_design( RelActCalc::RelEffEqnForm::LnX, 3 );
}


BOOST_AUTO_TEST_CASE( orthogonal_jet_matches_scalar_finite_difference )
{
  constexpr int num_derivatives = 4;
  using Jet = ceres::Jet<double,num_derivatives>;

  const RelActCalc::RelEffEqnForm form = RelActCalc::RelEffEqnForm::LnXLnY;
  const EmpiricalBasisTransform transform( form, 4, std::sqrt(610.0*775.0),
                                           frozen_energies() );
  std::vector<double> q{0.0, 0.25, -0.18, 0.08};
  std::vector<Jet> q_jet( q.begin(), q.end() );
  for( int i = 0; i < num_derivatives; ++i )
    q_jet[static_cast<size_t>(i)].v[i] = 1.0;

  const Jet predictor_jet = transform.linear_predictor( 693.4,
                                            q_jet.data(), q_jet.size() );
  const Jet value_jet = RelActCalc::continued_exp( predictor_jet, 50.0 );

  const double step = 1.0e-5;
  for( size_t i = 0; i < q.size(); ++i )
  {
    std::vector<double> lower = q, upper = q;
    lower[i] -= step;
    upper[i] += step;
    const double lower_value = std::exp(transform.linear_predictor(
                                693.4,lower.data(),lower.size()));
    const double upper_value = std::exp(transform.linear_predictor(
                                693.4,upper.data(),upper.size()));
    const double numeric = (upper_value - lower_value)/(2.0*step);
    BOOST_CHECK_SMALL( value_jet.v[static_cast<int>(i)] - numeric,
                       2.0e-10*(1.0 + std::fabs(numeric)) );
  }
  BOOST_CHECK_SMALL( value_jet.v[0], 1.0e-14 ); //the fixed gauge slot is mathematically absent
}


BOOST_AUTO_TEST_CASE( covariance_maps_to_legacy_coefficients_exactly )
{
  const RelActCalc::RelEffEqnForm form = RelActCalc::RelEffEqnForm::LnY;
  const EmpiricalBasisTransform transform( form, 3, std::sqrt(610.0*775.0),
                                           frozen_energies() );
  const std::vector<std::vector<double>> q_cov{
    {0.0, 0.0, 0.0},
    {0.0, 0.09, 0.015},
    {0.0, 0.015, 0.04}
  };
  const std::vector<std::vector<double>> legacy_cov = transform.legacy_covariance(q_cov);
  const std::vector<std::vector<double>> jacobian = transform.legacy_jacobian();

  // Compare the variance of an arbitrary linear functional in both coordinates.
  const std::vector<double> legacy_gradient{1.3, -0.2, 4.7};
  std::vector<double> q_gradient( 3, 0.0 );
  for( size_t q = 0; q < 3; ++q )
    for( size_t legacy = 0; legacy < 3; ++legacy )
      q_gradient[q] += legacy_gradient[legacy]*jacobian[legacy][q];

  double legacy_variance = 0.0, q_variance = 0.0;
  for( size_t i = 0; i < 3; ++i )
  {
    for( size_t j = 0; j < 3; ++j )
    {
      legacy_variance += legacy_gradient[i]*legacy_cov[i][j]*legacy_gradient[j];
      q_variance += q_gradient[i]*q_cov[i][j]*q_gradient[j];
    }
  }
  BOOST_CHECK_CLOSE_FRACTION( legacy_variance, q_variance, 3.0e-12 );
}
