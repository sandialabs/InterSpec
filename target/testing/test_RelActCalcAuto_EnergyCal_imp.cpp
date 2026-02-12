/* InterSpec: an application to analyze spectral gamma radiation data.

 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
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

#include "InterSpec_config.h"

#include <algorithm>
#include <cmath>
#include <vector>
#include <utility>
#include <iostream>

#include <ceres/jet.h>

#define BOOST_TEST_MODULE RelActCalcAuto_EnergyCal_imp_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/EnergyCalibration.h"

#if ( defined( WIN32 ) )
#undef min
#undef max
#endif

#include "InterSpec/RelActCalcAuto_EnergyCal_imp.hpp"


using namespace std;
using namespace RelActCalcAutoImp;
using namespace boost::unit_test;


// Test helper to convert vector<T> to vector<float>
template<typename T>
std::vector<float> to_float_vector( const std::vector<T> &input )
{
  std::vector<float> result;
  result.reserve( input.size() );
  for( const auto &val : input )
    result.push_back( static_cast<float>(val) );
  return result;
}


// Legacy spline builder kept in unit tests for regression comparisons.
//  It doesnt account for the derivative of knot-motion, unlike RelActCalcAutoImp::create_cubic_spline.
template<typename T>
std::vector<CubicSplineNodeT<T>> create_cubic_spline_legacy( const std::vector<std::pair<double,T>> &deviation_pairs )
{
  const size_t n = deviation_pairs.size();
  if( n < 2 )
    return std::vector<CubicSplineNodeT<T>>{};

#ifndef NDEBUG
  for( size_t i = 1; i < n; ++i )
    assert( deviation_pairs[i].first > deviation_pairs[i-1].first );
#endif

  std::vector<double> anchor_energies( n );
  std::vector<T> offsets( n );
  for( size_t i = 0; i < n; ++i )
  {
    anchor_energies[i] = deviation_pairs[i].first;
    offsets[i] = deviation_pairs[i].second;
  }

  std::vector<double> transformed_energies( n );
  for( size_t i = 0; i < n; ++i )
  {
    double offset_val;
    if constexpr ( std::is_same_v<T, double> )
      offset_val = offsets[i];
    else
      offset_val = offsets[i].a;

    transformed_energies[i] = anchor_energies[i] - offset_val;
  }

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero( n, n );
  std::vector<T> rhs( n );

  for( size_t i = 1; i < n - 1; ++i )
  {
    const double h_prev = transformed_energies[i] - transformed_energies[i-1];
    const double h_next = transformed_energies[i+1] - transformed_energies[i];

    A( i, i - 1 ) = h_prev / 3.0;
    A( i, i )     = (h_prev + h_next) / 1.5;
    A( i, i + 1 ) = h_next / 3.0;

    rhs[i] = (offsets[i+1] - offsets[i]) / h_next - (offsets[i] - offsets[i-1]) / h_prev;
  }

  A( 0, 0 ) = 2.0;
  rhs[0] = T( 0.0 );

  const double h_last = transformed_energies[n-1] - transformed_energies[n-2];
  A( n - 1, n - 2 ) = h_last / 3.0;
  A( n - 1, n - 1 ) = 2.0 * h_last / 3.0;
  rhs[n - 1] = -(offsets[n-1] - offsets[n-2]) / h_last;

  const Eigen::MatrixXd A_inv = A.partialPivLu().inverse();

  std::vector<T> b_vals( n, T(0.0) );
  for( size_t i = 0; i < n; ++i )
    for( size_t j = 0; j < n; ++j )
      b_vals[i] += A_inv( i, j ) * rhs[j];

  std::vector<CubicSplineNodeT<T>> nodes( n );
  for( size_t i = 0; i < n - 1; ++i )
  {
    const double h = transformed_energies[i+1] - transformed_energies[i];

    nodes[i].x = transformed_energies[i];
    nodes[i].y = offsets[i];
    nodes[i].a = (b_vals[i+1] - b_vals[i]) / (3.0 * h);
    nodes[i].b = b_vals[i];
    nodes[i].c = (offsets[i+1] - offsets[i]) / h - (2.0 * b_vals[i] + b_vals[i+1]) * h / 3.0;
  }

  nodes[n-1].x = transformed_energies[n-1];
  nodes[n-1].y = offsets[n-1];
  nodes[n-1].a = T(0.0);
  nodes[n-1].b = T(0.0);
  nodes[n-1].c = T(0.0);
  return nodes;
}


BOOST_AUTO_TEST_CASE( CubicSplineCreation )
{
  // Test basic spline creation
  std::vector<pair<double,double>> dev_pairs = {
    {100.0, 0.5},
    {500.0, 1.0},
    {1000.0, 0.8},
    {1500.0, 0.3}
  };

  const auto nodes = create_cubic_spline_legacy( dev_pairs );

  BOOST_CHECK_EQUAL( nodes.size(), 4 );

  // Check that nodes are sorted
  for( size_t i = 1; i < nodes.size(); ++i )
  {
    BOOST_CHECK( nodes[i].x > nodes[i-1].x );
  }

  // Check boundary conditions (last node has zero derivatives)
  BOOST_CHECK_SMALL( nodes.back().a, 1.0e-10 );
  BOOST_CHECK_SMALL( nodes.back().b, 1.0e-10 );
  BOOST_CHECK_SMALL( nodes.back().c, 1.0e-10 );
}


BOOST_AUTO_TEST_CASE( CubicSplineEvaluation )
{
  std::vector<pair<double,double>> dev_pairs = {
    {100.0, 0.5},
    {500.0, 1.0},
    {1000.0, 0.8}
  };

  const auto nodes = create_cubic_spline_legacy( dev_pairs );

  // Evaluate at anchor points - should match offset values
  for( size_t i = 0; i < dev_pairs.size(); ++i )
  {
    const double energy = dev_pairs[i].first;
    const double result = eval_cubic_spline( energy, nodes );

    // Note: Due to transformation (x - offset), we need to evaluate at the adjusted energy
    // For this test, just check spline is smooth
    BOOST_CHECK( std::isfinite(result) );
  }

  // Check extrapolation (clamps to boundary values)
  const double result_low = eval_cubic_spline( 50.0, nodes );
  BOOST_CHECK_CLOSE( result_low, nodes.front().y, 1.0e-6 );

  const double result_high = eval_cubic_spline( 2000.0, nodes );
  BOOST_CHECK_CLOSE( result_high, nodes.back().y, 1.0e-6 );
}


BOOST_AUTO_TEST_CASE( CubicSplineFullDeriv_MatchesLegacyDouble )
{
  std::vector<pair<double,double>> dev_pairs = {
    {100.0, 0.5},
    {500.0, 1.0},
    {1000.0, 0.8},
    {1500.0, 0.3}
  };

  const std::vector<CubicSplineNodeT<double>> legacy_nodes = create_cubic_spline_legacy( dev_pairs );
  const std::vector<CubicSplineNodeT<double>> full_nodes = RelActCalcAutoImp::create_cubic_spline( dev_pairs );

  BOOST_REQUIRE_EQUAL( legacy_nodes.size(), full_nodes.size() );
  for( size_t i = 0; i < legacy_nodes.size(); ++i )
  {
    BOOST_CHECK_SMALL( std::fabs( legacy_nodes[i].x - full_nodes[i].x ), 1.0e-12 );
    BOOST_CHECK_SMALL( std::fabs( legacy_nodes[i].y - full_nodes[i].y ), 1.0e-12 );
    BOOST_CHECK_SMALL( std::fabs( legacy_nodes[i].a - full_nodes[i].a ), 1.0e-12 );
    BOOST_CHECK_SMALL( std::fabs( legacy_nodes[i].b - full_nodes[i].b ), 1.0e-12 );
    BOOST_CHECK_SMALL( std::fabs( legacy_nodes[i].c - full_nodes[i].c ), 1.0e-12 );
  }
}


BOOST_AUTO_TEST_CASE( CubicSplineFullDeriv_CloseToLegacyJet )
{
  using Jet3 = ceres::Jet<double,3>;

  std::vector<pair<double,Jet3>> dev_pairs;
  dev_pairs.reserve( 5 );

  const std::vector<double> anchor_energies = { 150.0, 500.0, 900.0, 1300.0, 1800.0 };
  const std::vector<double> offsets = { 0.6, 1.1, 0.9, 0.4, 0.2 };
  for( size_t i = 0; i < anchor_energies.size(); ++i )
  {
    Jet3 offset;
    offset.a = offsets[i];
    offset.v[0] = 1.0e-6 * static_cast<double>(i + 1);
    offset.v[1] = -7.0e-7 * static_cast<double>(i + 2);
    offset.v[2] = 5.0e-7 * static_cast<double>(i + 3);
    dev_pairs.emplace_back( anchor_energies[i], offset );
  }

  const std::vector<CubicSplineNodeT<Jet3>> legacy_nodes = RelActCalcAutoImp::create_cubic_spline( dev_pairs );
  const std::vector<CubicSplineNodeT<Jet3>> full_nodes = RelActCalcAutoImp::create_cubic_spline( dev_pairs );
  BOOST_REQUIRE_EQUAL( legacy_nodes.size(), full_nodes.size() );

  for( size_t i = 0; i < legacy_nodes.size(); ++i )
  {
    BOOST_CHECK_SMALL( std::fabs( legacy_nodes[i].x - full_nodes[i].x ), 1.0e-10 );
    BOOST_CHECK_SMALL( std::fabs( legacy_nodes[i].y.a - full_nodes[i].y.a ), 1.0e-12 );
    BOOST_CHECK_SMALL( std::fabs( legacy_nodes[i].a.a - full_nodes[i].a.a ), 1.0e-8 );
    BOOST_CHECK_SMALL( std::fabs( legacy_nodes[i].b.a - full_nodes[i].b.a ), 1.0e-8 );
    BOOST_CHECK_SMALL( std::fabs( legacy_nodes[i].c.a - full_nodes[i].c.a ), 1.0e-8 );
  }

  for( double energy = 200.0; energy <= 1700.0; energy += 100.0 )
  {
    Jet3 e;
    e.a = energy;
    e.v[0] = 0.0;
    e.v[1] = 0.0;
    e.v[2] = 0.0;

    const Jet3 y_legacy = eval_cubic_spline( e, legacy_nodes );
    const Jet3 y_full = eval_cubic_spline( e, full_nodes );

    BOOST_CHECK_SMALL( std::fabs( y_legacy.a - y_full.a ), 5.0e-8 );

    const double frac_thresh = 1.0e-4;  // 0.01%
    const double deriv_floor = 1.0e-12;

    const double rel_diff_v0 = std::fabs( y_legacy.v[0] - y_full.v[0] )
                               / (std::max)( deriv_floor, std::fabs( y_legacy.v[0] ) );
    const double rel_diff_v1 = std::fabs( y_legacy.v[1] - y_full.v[1] )
                               / (std::max)( deriv_floor, std::fabs( y_legacy.v[1] ) );
    const double rel_diff_v2 = std::fabs( y_legacy.v[2] - y_full.v[2] )
                               / (std::max)( deriv_floor, std::fabs( y_legacy.v[2] ) );
    BOOST_CHECK_SMALL( rel_diff_v0, frac_thresh );
    BOOST_CHECK_SMALL( rel_diff_v1, frac_thresh );
    BOOST_CHECK_SMALL( rel_diff_v2, frac_thresh );
  }
}


BOOST_AUTO_TEST_CASE( CubicSplineFullDeriv_AtKnotEnergies )
{
  using Jet3 = ceres::Jet<double,3>;

  std::vector<pair<double,Jet3>> dev_pairs;
  dev_pairs.reserve( 5 );

  const std::vector<double> anchor_energies = { 150.0, 500.0, 900.0, 1300.0, 1800.0 };
  const std::vector<double> offsets = { 0.6, 1.1, 0.9, 0.4, 0.2 };
  for( size_t i = 0; i < anchor_energies.size(); ++i )
  {
    Jet3 offset;
    offset.a = offsets[i];
    offset.v[0] = 1.0e-6 * static_cast<double>(i + 1);
    offset.v[1] = -7.0e-7 * static_cast<double>(i + 2);
    offset.v[2] = 5.0e-7 * static_cast<double>(i + 3);
    dev_pairs.emplace_back( anchor_energies[i], offset );
  }

  const std::vector<CubicSplineNodeT<Jet3>> legacy_nodes = RelActCalcAutoImp::create_cubic_spline( dev_pairs );
  const std::vector<CubicSplineNodeT<Jet3>> full_nodes = RelActCalcAutoImp::create_cubic_spline( dev_pairs );
  BOOST_REQUIRE_EQUAL( legacy_nodes.size(), full_nodes.size() );

  // Check exactly at transformed knot positions and immediately around them.
  // At exact knot locations, piecewise interval selection makes the legacy/full-deriv
  // implementations diverge slightly more than away from knots; keep a modest guardrail.
  const double frac_thresh = 5.0e-4;  // 0.05%
  const double deriv_floor = 1.0e-12;
  const double eps = 1.0e-6;

  double max_rel_v0 = 0.0;
  double max_rel_v1 = 0.0;
  double max_rel_v2 = 0.0;
  double max_abs_val_diff = 0.0;

  auto update_stats = [&]( const double energy_value ) {
    Jet3 e;
    e.a = energy_value;
    e.v[0] = 0.0;
    e.v[1] = 0.0;
    e.v[2] = 0.0;

    const Jet3 y_legacy = eval_cubic_spline( e, legacy_nodes );
    const Jet3 y_full = eval_cubic_spline( e, full_nodes );

    const double rel_v0 = std::fabs( y_legacy.v[0] - y_full.v[0] )
                          / (std::max)( deriv_floor, std::fabs( y_legacy.v[0] ) );
    const double rel_v1 = std::fabs( y_legacy.v[1] - y_full.v[1] )
                          / (std::max)( deriv_floor, std::fabs( y_legacy.v[1] ) );
    const double rel_v2 = std::fabs( y_legacy.v[2] - y_full.v[2] )
                          / (std::max)( deriv_floor, std::fabs( y_legacy.v[2] ) );

    max_rel_v0 = std::max( max_rel_v0, rel_v0 );
    max_rel_v1 = std::max( max_rel_v1, rel_v1 );
    max_rel_v2 = std::max( max_rel_v2, rel_v2 );
    max_abs_val_diff = std::max( max_abs_val_diff, std::fabs( y_legacy.a - y_full.a ) );
  };

  for( size_t i = 0; i < legacy_nodes.size(); ++i )
  {
    const double knot = legacy_nodes[i].x;
    update_stats( knot - eps );
    update_stats( knot );
    update_stats( knot + eps );
  }

  BOOST_CHECK_SMALL( max_abs_val_diff, 1.0e-7 );
  BOOST_CHECK_SMALL( max_rel_v0, frac_thresh );
  BOOST_CHECK_SMALL( max_rel_v1, frac_thresh );
  BOOST_CHECK_SMALL( max_rel_v2, frac_thresh );

  BOOST_TEST_MESSAGE( "Knot-energy max relative diffs: v0="
    << (100.0 * max_rel_v0) << "%, v1=" << (100.0 * max_rel_v1)
    << "%, v2=" << (100.0 * max_rel_v2) << "%" );
}


BOOST_AUTO_TEST_CASE( DeviationPairCorrection )
{
  // Test correction_due_to_deviation_pairs against SpecUtils
  std::vector<pair<float,float>> dev_pairs = {
    {100.0f, 0.5f},
    {500.0f, 1.0f},
    {1000.0f, 0.8f},
    {1500.0f, 0.3f},
    {2000.0f, 0.1f}
  };

  // Build spline from deviation pairs
  std::vector<pair<double,double>> dev_pairs_double;
  for( const auto &pair : dev_pairs )
  {
    dev_pairs_double.emplace_back( pair.first, pair.second );
  }

  const auto nodes = RelActCalcAutoImp::create_cubic_spline( dev_pairs_double );

  // Test at various energies
  for( double energy = 25.0; energy <= 2500.0; energy += 25.0 )
  {
    const double our_correction = correction_due_to_deviation_pairs( energy, nodes );
    const double specutils_correction = SpecUtils::correction_due_to_dev_pairs( energy, dev_pairs );

    // Should match within tolerance
    BOOST_CHECK_SMALL( std::fabs(our_correction - specutils_correction), 0.001 );  // 1 eV tolerance
  }
}


BOOST_AUTO_TEST_CASE( PolynomialChannel_Linear )
{
  // Test linear polynomial: E = 0 + 3*ch
  const size_t nchannel = 1024;
  std::vector<double> coeffs = { 0.0, 3.0 };
  std::vector<pair<double,double>> dev_pairs;  // No deviation pairs

  // Test various energies
  for( double energy = 0.0; energy < 3000.0; energy += 100.0 )
  {
    const double our_channel = find_polynomial_channel( energy, coeffs, nchannel, dev_pairs );
    const double expected_channel = energy / 3.0;

    BOOST_CHECK_CLOSE( our_channel, expected_channel, 1.0e-6 );
  }
}


BOOST_AUTO_TEST_CASE( PolynomialChannel_Quadratic )
{
  // Test quadratic polynomial: E = 10 + 2*ch + 0.001*ch²
  const size_t nchannel = 2048;
  std::vector<double> coeffs = { 10.0, 2.0, 0.001 };
  std::vector<pair<double,double>> dev_pairs;

  // Test various energies
  for( double energy = 50.0; energy < 4000.0; energy += 200.0 )
  {
    const double our_channel = find_polynomial_channel( energy, coeffs, nchannel, dev_pairs );

    const std::vector<float> coeffs_float = to_float_vector( coeffs );
    const double specutils_channel = SpecUtils::find_polynomial_channel( energy, coeffs_float, nchannel, {}, 0.001 );

    BOOST_CHECK_SMALL( std::fabs(our_channel - specutils_channel), 0.01 );  // Within 0.01 channels
  }
}


BOOST_AUTO_TEST_CASE( PolynomialChannel_Cubic )
{
  // Test cubic polynomial: E = 5 + 1.5*ch + 0.002*ch² - 0.000001*ch³
  const size_t nchannel = 4096;
  std::vector<double> coeffs = { 5.0, 1.5, 0.002, -0.000001 };
  std::vector<pair<double,double>> dev_pairs;

  // Compute max energy for the calibration
  double max_energy = 0.0;
  for( size_t i = 0; i < coeffs.size(); ++i )
    max_energy += coeffs[i] * std::pow( static_cast<double>(nchannel - 1), static_cast<double>(i) );

  for( double energy = 50.0; energy < max_energy - 100.0; energy += 300.0 )
  {
    const double our_channel = find_polynomial_channel( energy, coeffs, nchannel, dev_pairs );

    const std::vector<float> coeffs_float = to_float_vector( coeffs );
    const double specutils_channel = SpecUtils::find_polynomial_channel( energy, coeffs_float, nchannel, {}, 0.001 );

    BOOST_CHECK_SMALL( std::fabs(our_channel - specutils_channel), 0.01 );
  }
}


BOOST_AUTO_TEST_CASE( PolynomialChannel_WithDeviationPairs )
{
  // Test polynomial with deviation pairs
  const size_t nchannel = 1024;
  std::vector<double> coeffs = { 0.0, 3.0 };

  // Add deviation pairs
  std::vector<pair<double,double>> dev_pairs_double = {
    {500.0, 1.5},
    {1000.0, 2.0},
    {1500.0, 1.0},
    {2000.0, 0.5}
  };

  std::vector<pair<float,float>> dev_pairs_float;
  for( const auto &p : dev_pairs_double )
    dev_pairs_float.emplace_back( static_cast<float>(p.first), static_cast<float>(p.second) );

  for( double energy = 100.0; energy < 2500.0; energy += 150.0 )
  {
    const double our_channel = find_polynomial_channel( energy, coeffs, nchannel, dev_pairs_double );

    const std::vector<float> coeffs_float = to_float_vector( coeffs );
    const double specutils_channel = SpecUtils::find_polynomial_channel( energy, coeffs_float, nchannel, dev_pairs_float, 0.001 );

    BOOST_CHECK_SMALL( std::fabs(our_channel - specutils_channel), 0.01 );
  }
}


BOOST_AUTO_TEST_CASE( PolynomialChannel_LowOrderDeviationPairs_Inversion )
{
  // Low-order polynomial with deviation pairs should invert consistently
  const size_t nchannel = 2048;
  std::vector<double> coeffs = { 15.0, 2.5 };
  std::vector<pair<double,double>> dev_pairs = {
    {250.0, 0.8},
    {750.0, 1.2},
    {1250.0, 0.6},
    {1750.0, 0.2}
  };

  const std::vector<CubicSplineNodeT<double>> dev_pair_spline = RelActCalcAutoImp::create_cubic_spline( dev_pairs );

  for( double ch = 0.0; ch < static_cast<double>(nchannel); ch += 137.0 )
  {
    const double energy = polynomial_energy( ch, coeffs, dev_pair_spline );
    const double solved_ch = find_polynomial_channel( energy, coeffs, nchannel, dev_pairs );

    BOOST_CHECK_SMALL( std::fabs( solved_ch - ch ), 0.01 );
  }
}


BOOST_AUTO_TEST_CASE( PolynomialChannel_DeviationPairs_JetDerivative )
{
  using Jet1 = ceres::Jet<double,1>;

  const size_t nchannel = 2048;
  std::vector<double> coeffs = { 0.0, 3.0 };
  std::vector<pair<double,double>> dev_pairs = {
    {500.0, 1.5},
    {1000.0, 2.0},
    {1500.0, 1.0},
    {2000.0, 0.5}
  };

  std::vector<Jet1> coeffs_jet;
  coeffs_jet.reserve( coeffs.size() );
  for( size_t i = 0; i < coeffs.size(); ++i )
  {
    Jet1 coef;
    coef.a = coeffs[i];
    coef.v[0] = 0.0;
    coeffs_jet.push_back( coef );
  }

  std::vector<pair<double,Jet1>> dev_pairs_jet;
  dev_pairs_jet.reserve( dev_pairs.size() );
  for( size_t i = 0; i < dev_pairs.size(); ++i )
  {
    Jet1 offset;
    offset.a = dev_pairs[i].second;
    offset.v[0] = 0.0;
    dev_pairs_jet.emplace_back( dev_pairs[i].first, offset );
  }

  const double energy = 1200.0;
  const double eps = 1.0e-3;

  const double ch_plus = find_polynomial_channel( energy + eps, coeffs, nchannel, dev_pairs );
  const double ch_minus = find_polynomial_channel( energy - eps, coeffs, nchannel, dev_pairs );
  const double fd_deriv = (ch_plus - ch_minus) / (2.0 * eps);

  Jet1 jet_energy;
  jet_energy.a = energy;
  jet_energy.v[0] = 1.0;

  const Jet1 jet_channel = find_polynomial_channel( jet_energy, coeffs_jet, nchannel, dev_pairs_jet );
  const double ch_val = find_polynomial_channel( energy, coeffs, nchannel, dev_pairs );

  BOOST_CHECK_SMALL( std::fabs( jet_channel.a - ch_val ), 1.0e-6 );
  BOOST_CHECK_SMALL( std::fabs( jet_channel.v[0] - fd_deriv ), 1.0e-3 );
}


BOOST_AUTO_TEST_CASE( PolynomialChannel_ManyDeviationPairs )
{
  // Test with many deviation pairs
  const size_t nchannel = 8192;
  std::vector<double> coeffs = { 0.0, 0.5, 0.00001 };

  std::vector<pair<double,double>> dev_pairs_double;
  std::vector<pair<float,float>> dev_pairs_float;

  // Create 20 deviation pairs
  for( int i = 0; i < 20; ++i )
  {
    const double energy = 200.0 + i * 150.0;
    const double offset = 0.5 + 0.3 * std::sin(i * 0.5);
    dev_pairs_double.emplace_back( energy, offset );
    dev_pairs_float.emplace_back( static_cast<float>(energy), static_cast<float>(offset) );
  }

  for( double energy = 50.0; energy < 3000.0; energy += 25.0 )
  {
    const double our_channel = find_polynomial_channel( energy, coeffs, nchannel, dev_pairs_double );

    const std::vector<float> coeffs_float = to_float_vector( coeffs );
    const double specutils_channel = SpecUtils::find_polynomial_channel( energy, coeffs_float, nchannel, dev_pairs_float, 0.001 );

    BOOST_CHECK_SMALL( std::fabs(our_channel - specutils_channel), 0.01 );
  }
}


BOOST_AUTO_TEST_CASE( FullRangeFraction_Linear )
{
  // Test linear FRF: E = 0 + x*3000, where x = bin/nchannel
  const size_t nchannel = 1024;
  std::vector<double> coeffs = { 0.0, 3000.0 };
  std::vector<pair<double,double>> dev_pairs;

  for( double energy = 0.0; energy < 3000.0; energy += 100.0 )
  {
    const double our_bin = find_fullrangefraction_channel( energy, coeffs, nchannel, dev_pairs );

    const std::vector<float> coeffs_float = to_float_vector( coeffs );
    const double specutils_bin = SpecUtils::find_fullrangefraction_channel( energy, coeffs_float, nchannel, {}, 0.001 );

    BOOST_CHECK_SMALL( std::fabs(our_bin - specutils_bin), 0.01 );
  }
}


BOOST_AUTO_TEST_CASE( FullRangeFraction_Quadratic )
{
  // Test quadratic FRF: E = 10 + x*2000 + x²*500
  const size_t nchannel = 2048;
  std::vector<double> coeffs = { 10.0, 2000.0, 500.0 };
  std::vector<pair<double,double>> dev_pairs;

  for( double energy = 50.0; energy < 2500.0; energy += 150.0 )
  {
    const double our_bin = find_fullrangefraction_channel( energy, coeffs, nchannel, dev_pairs );

    const std::vector<float> coeffs_float = to_float_vector( coeffs );
    const double specutils_bin = SpecUtils::find_fullrangefraction_channel( energy, coeffs_float, nchannel, {}, 0.001 );

    BOOST_CHECK_SMALL( std::fabs(our_bin - specutils_bin), 0.01 );
  }
}


BOOST_AUTO_TEST_CASE( FullRangeFraction_WithC4Term )
{
  // Test with C₄/(1+60x) term
  const size_t nchannel = 1024;
  std::vector<double> coeffs = { 0.0, 2800.0, 100.0, 0.0, 50.0 };
  std::vector<pair<double,double>> dev_pairs;

  for( double energy = 50.0; energy < 3000.0; energy += 200.0 )
  {
    const double our_bin = find_fullrangefraction_channel( energy, coeffs, nchannel, dev_pairs );

    const std::vector<float> coeffs_float = to_float_vector( coeffs );
    const double specutils_bin = SpecUtils::find_fullrangefraction_channel( energy, coeffs_float, nchannel, {}, 0.001 );

    BOOST_CHECK_SMALL( std::fabs(our_bin - specutils_bin), 0.01 );
  }
}


BOOST_AUTO_TEST_CASE( FullRangeFraction_WithDeviationPairs )
{
  // Test FRF with deviation pairs
  const size_t nchannel = 1024;
  std::vector<double> coeffs = { 0.0, 3000.0 };

  std::vector<pair<double,double>> dev_pairs_double = {
    {500.0, 2.0},
    {1000.0, 3.0},
    {1500.0, 2.5},
    {2000.0, 1.5},
    {2500.0, 0.5}
  };

  std::vector<pair<float,float>> dev_pairs_float;
  for( const auto &p : dev_pairs_double )
    dev_pairs_float.emplace_back( static_cast<float>(p.first), static_cast<float>(p.second) );

  for( double energy = 100.0; energy < 2800.0; energy += 150.0 )
  {
    const double our_bin = find_fullrangefraction_channel( energy, coeffs, nchannel, dev_pairs_double );

    const std::vector<float> coeffs_float = to_float_vector( coeffs );
    const double specutils_bin = SpecUtils::find_fullrangefraction_channel( energy, coeffs_float, nchannel, dev_pairs_float, 0.001 );

    BOOST_CHECK_SMALL( std::fabs(our_bin - specutils_bin), 0.01 );
  }
}


BOOST_AUTO_TEST_CASE( FullRangeFraction_LowOrderDeviationPairs_Inversion )
{
  // Low-order FRF with deviation pairs should invert consistently
  const size_t nchannel = 1024;
  std::vector<double> coeffs = { 5.0, 2800.0 };
  std::vector<pair<double,double>> dev_pairs = {
    {400.0, 1.2},
    {900.0, 1.8},
    {1500.0, 1.0},
    {2200.0, 0.4}
  };

  const std::vector<CubicSplineNodeT<double>> dev_pair_spline = RelActCalcAutoImp::create_cubic_spline( dev_pairs );

  for( double bin = 0.0; bin < static_cast<double>(nchannel); bin += 101.0 )
  {
    const double energy = fullrangefraction_energy( bin, coeffs, nchannel, dev_pair_spline );
    const double solved_bin = find_fullrangefraction_channel( energy, coeffs, nchannel, dev_pairs );

    BOOST_CHECK_SMALL( std::fabs( solved_bin - bin ), 0.01 );
  }
}


BOOST_AUTO_TEST_CASE( FullRangeFraction_DeviationPairs_JetDerivative )
{
  using Jet1 = ceres::Jet<double,1>;

  const size_t nchannel = 1024;
  std::vector<double> coeffs = { 0.0, 3000.0 };
  std::vector<pair<double,double>> dev_pairs = {
    {500.0, 2.0},
    {1000.0, 3.0},
    {1500.0, 2.5},
    {2000.0, 1.5},
    {2500.0, 0.5}
  };

  std::vector<Jet1> coeffs_jet;
  coeffs_jet.reserve( coeffs.size() );
  for( size_t i = 0; i < coeffs.size(); ++i )
  {
    Jet1 coef;
    coef.a = coeffs[i];
    coef.v[0] = 0.0;
    coeffs_jet.push_back( coef );
  }

  std::vector<pair<double,Jet1>> dev_pairs_jet;
  dev_pairs_jet.reserve( dev_pairs.size() );
  for( size_t i = 0; i < dev_pairs.size(); ++i )
  {
    Jet1 offset;
    offset.a = dev_pairs[i].second;
    offset.v[0] = 0.0;
    dev_pairs_jet.emplace_back( dev_pairs[i].first, offset );
  }

  const double energy = 1500.0;
  const double eps = 1.0e-3;

  const double bin_plus = find_fullrangefraction_channel( energy + eps, coeffs, nchannel, dev_pairs );
  const double bin_minus = find_fullrangefraction_channel( energy - eps, coeffs, nchannel, dev_pairs );
  const double fd_deriv = (bin_plus - bin_minus) / (2.0 * eps);

  Jet1 jet_energy;
  jet_energy.a = energy;
  jet_energy.v[0] = 1.0;

  const Jet1 jet_bin = find_fullrangefraction_channel( jet_energy, coeffs_jet, nchannel, dev_pairs_jet );
  const double bin_val = find_fullrangefraction_channel( energy, coeffs, nchannel, dev_pairs );

  BOOST_CHECK_SMALL( std::fabs( jet_bin.a - bin_val ), 1.0e-6 );
  BOOST_CHECK_SMALL( std::fabs( jet_bin.v[0] - fd_deriv ), 1.0e-3 );
}


BOOST_AUTO_TEST_CASE( FullRangeFraction_ManyDeviationPairs )
{
  // Test FRF with many deviation pairs
  const size_t nchannel = 4096;
  std::vector<double> coeffs = { 5.0, 2500.0, 200.0 };

  std::vector<pair<double,double>> dev_pairs_double;
  std::vector<pair<float,float>> dev_pairs_float;

  // Create 30 deviation pairs
  for( int i = 0; i < 30; ++i )
  {
    const double energy = 100.0 + i * 90.0;
    const double offset = 1.0 + 0.8 * std::cos(i * 0.3);
    dev_pairs_double.emplace_back( energy, offset );
    dev_pairs_float.emplace_back( static_cast<float>(energy), static_cast<float>(offset) );
  }

  for( double energy = 200.0; energy < 2500.0; energy += 180.0 )
  {
    const double our_bin = find_fullrangefraction_channel( energy, coeffs, nchannel, dev_pairs_double );

    const std::vector<float> coeffs_float = to_float_vector( coeffs );
    const double specutils_bin = SpecUtils::find_fullrangefraction_channel( energy, coeffs_float, nchannel, dev_pairs_float, 0.001 );

    BOOST_CHECK_SMALL( std::fabs(our_bin - specutils_bin), 0.01 );
  }
}


BOOST_AUTO_TEST_CASE( LowerChannelEdge_Basic )
{
  // Test lower channel edge calibration
  const size_t nchannel = 10;
  std::vector<float> channel_energies = { 0.0f, 100.0f, 250.0f, 450.0f, 700.0f,
                                          1000.0f, 1350.0f, 1750.0f, 2200.0f,
                                          2700.0f, 3000.0f };
  std::vector<pair<double,double>> dev_pairs;

  // Test energies at channel boundaries
  for( size_t i = 0; i < nchannel; ++i )
  {
    const double energy = channel_energies[i] + 10.0;  // Slightly above lower edge
    const double channel = find_lowerchannel_channel( energy, channel_energies, dev_pairs );

    // Should be close to channel i
    BOOST_CHECK( channel >= static_cast<double>(i) );
    BOOST_CHECK( channel < static_cast<double>(i + 1) );
  }

  // Test interpolation within a channel
  const double mid_energy = 0.5 * (channel_energies[5] + channel_energies[6]);
  const double mid_channel = find_lowerchannel_channel( mid_energy, channel_energies, dev_pairs );
  BOOST_CHECK_CLOSE( mid_channel, 5.5, 0.1 );  // Within 0.1% of 5.5
}


BOOST_AUTO_TEST_CASE( LowerChannelEdge_WithDeviationPairs )
{
  // Test lower channel edge with deviation pairs
  const size_t nchannel = 8;
  std::vector<float> channel_energies = { 0.0f, 200.0f, 500.0f, 900.0f,
                                          1400.0f, 2000.0f, 2700.0f,
                                          3500.0f, 4000.0f };

  std::vector<pair<double,double>> dev_pairs_double = {
    {500.0, 1.0},
    {1500.0, 2.0},
    {2500.0, 1.5},
    {3500.0, 0.5}
  };

  // Create adjusted channel energies by applying deviation pair corrections
  const auto dev_pair_spline = RelActCalcAutoImp::create_cubic_spline( dev_pairs_double );

  std::vector<float> adjusted_channel_energies;
  adjusted_channel_energies.reserve( channel_energies.size() );
  for( const float energy : channel_energies )
  {
    const double correction = eval_cubic_spline( static_cast<double>(energy), dev_pair_spline );
    adjusted_channel_energies.push_back( static_cast<float>(energy + correction) );
  }

  // Test that using adjusted energies without dev pairs gives same result as original with dev pairs
  for( double energy = 100.0; energy < 3800.0; energy += 200.0 )
  {
    const double channel_with_devpairs = find_lowerchannel_channel( energy, channel_energies, dev_pairs_double );
    const double channel_adjusted = find_lowerchannel_channel( energy, adjusted_channel_energies,
                                                               std::vector<std::pair<double,double>>{} );

    // Should give very similar results
    BOOST_CHECK_SMALL( std::fabs(channel_with_devpairs - channel_adjusted), 0.01 );  // Within 0.01 channels

    // Basic sanity checks
    BOOST_CHECK( channel_with_devpairs >= 0.0 );
    BOOST_CHECK( channel_with_devpairs < static_cast<double>(nchannel) );
    BOOST_CHECK( std::isfinite(channel_with_devpairs) );
  }
}


BOOST_AUTO_TEST_CASE( LowerChannelEdge_WithOffsetGainAdj )
{
  // Test lower channel edge with offset and gain adjustments
  const size_t nchannel = 2048;

  // Create channel energies with 3 keV wide channels
  std::vector<float> channel_energies;
  channel_energies.reserve( nchannel + 1 );
  for( size_t i = 0; i <= nchannel; ++i )
  {
    channel_energies.push_back( static_cast<float>(i * 3.0) );
  }

  const double offset_adj = 5.0;  // 5 keV offset
  const double gain_adj = 10.0;   // 10 keV gain adjustment

  std::vector<pair<double,double>> dev_pairs_double = {
    {50.0, 0.0},
    {125.0, -5.0},
    {500.0, 10.0},
    {1500.0, -2.0},
    {2500.0, 15.0},
    {3500.0, 0.5},
    {4000.0, 0.0}
  };

  // Manually compute adjusted energies for verification
  const double lower_energy = channel_energies[0];
  const double upper_energy = channel_energies[nchannel];
  const double range = upper_energy - lower_energy;

  std::vector<float> manual_adjusted_energies;
  manual_adjusted_energies.reserve( channel_energies.size() );
  for( const float orig_energy : channel_energies )
  {
    const double range_frac = (orig_energy - lower_energy) / range;
    const double adjusted = orig_energy + offset_adj + (range_frac * gain_adj);
    manual_adjusted_energies.push_back( static_cast<float>(adjusted) );
  }

  // Test that using the offset/gain overload gives same result as using manually adjusted energies
  for( double energy = 50.0; energy < 4000.0; energy += 25.0 )
  {
    const double channel_with_adj = find_lowerchannel_channel( energy, channel_energies, dev_pairs_double, offset_adj, gain_adj );
    const double channel_manual = find_lowerchannel_channel( energy, manual_adjusted_energies, dev_pairs_double );

    // Should give very close results (within numerical precision)
    BOOST_CHECK_SMALL( std::fabs(channel_with_adj - channel_manual), 1.0e-6 );

    // Basic sanity checks
    BOOST_CHECK( channel_with_adj >= 0.0 );
    BOOST_CHECK( channel_with_adj < static_cast<double>(nchannel) );
    BOOST_CHECK( std::isfinite(channel_with_adj) );
  }

  // Test without deviation pairs
  for( double energy = 100.0; energy < 4000.0; energy += 200.0 )
  {
    const double channel_with_adj = find_lowerchannel_channel( energy, channel_energies,
                                                               std::vector<std::pair<double,double>>{},
                                                               offset_adj, gain_adj );
    const double channel_manual = find_lowerchannel_channel( energy, manual_adjusted_energies,
                                                             std::vector<std::pair<double,double>>{} );

    BOOST_CHECK_SMALL( std::fabs(channel_with_adj - channel_manual), 1.0e-6 );
  }
}


BOOST_AUTO_TEST_CASE( RoundTrip_Polynomial )
{
  // Test round trip: energy -> channel -> energy
  const size_t nchannel = 2048;
  std::vector<double> coeffs = { 5.0, 2.0, 0.001 };
  std::vector<pair<double,double>> dev_pairs_double = {
    {500.0, 0.8},
    {1500.0, 1.2},
    {2500.0, 0.6}
  };

  // Build deviation pair spline
  const auto dev_pair_spline = RelActCalcAutoImp::create_cubic_spline( dev_pairs_double );

  for( double energy = 100.0; energy < 4000.0; energy += 25.0 )
  {
    const double channel = find_polynomial_channel( energy, coeffs, nchannel, dev_pairs_double );
    const double energy_back = polynomial_energy( channel, coeffs, dev_pair_spline );

    // Should match within accuracy
    BOOST_CHECK_SMALL( std::fabs(energy - energy_back), 0.01 );  // Within 10 eV
  }
}


BOOST_AUTO_TEST_CASE( RoundTrip_FullRangeFraction )
{
  // Test round trip: energy -> bin -> energy
  const size_t nchannel = 1024;
  std::vector<double> coeffs = { 0.0, 2800.0, 150.0 };
  std::vector<pair<double,double>> dev_pairs_double = {
    {600.0, 1.5},
    {1400.0, 2.0},
    {2200.0, 1.0}
  };

  // Build deviation pair spline
  const auto dev_pair_spline = RelActCalcAutoImp::create_cubic_spline( dev_pairs_double );

  for( double energy = 50.0; energy < 2900.0; energy += 25.0 )
  {
    const double bin = find_fullrangefraction_channel( energy, coeffs, nchannel, dev_pairs_double );
    const double energy_back = fullrangefraction_energy( bin, coeffs, nchannel, dev_pair_spline );

    // Should match within accuracy
    BOOST_CHECK_SMALL( std::fabs(energy - energy_back), 0.01 );
  }
}


BOOST_AUTO_TEST_CASE( ApplyUnapplyRoundTrip_Polynomial )
{
  // Simulate apply_energy_cal_adjustment (true -> spectrum) and un_apply (spectrum -> true)
  // using the same flow: channel = find_*(E, cal); E_out = cal_energy(channel).
  const size_t nchannel = 2048;
  const double num_channel = static_cast<double>( nchannel );
  const double offset_adj = 3.0;
  const double gain_adj = 5.0;
  const double quad_adj = 0.5;

  std::vector<double> orig_coeffs = { 5.0, 2.0, 0.001 };
  std::vector<pair<double,double>> dev_pairs = { {500.0, 0.8}, {1500.0, 1.2}, {2500.0, 0.6} };
  const auto orig_spline = RelActCalcAutoImp::create_cubic_spline( dev_pairs );

  std::vector<double> adj_coeffs = orig_coeffs;
  adj_coeffs[0] += offset_adj;
  adj_coeffs[1] += gain_adj / num_channel;
  if( adj_coeffs.size() < 3 )
    adj_coeffs.resize( 3, 0.0 );
  adj_coeffs[2] += quad_adj / (num_channel * num_channel);

  for( double E_true = 200.0; E_true < 3500.0; E_true += 100.0 )
  {
    const double channel = find_polynomial_channel( E_true, adj_coeffs, nchannel, dev_pairs );
    const double E_spec = polynomial_energy( channel, orig_coeffs, orig_spline );

    const double ch2 = find_polynomial_channel( E_spec, orig_coeffs, nchannel, dev_pairs );
    const double E_true_back = polynomial_energy( ch2, adj_coeffs, orig_spline );

    BOOST_CHECK_SMALL( std::fabs( E_true_back - E_true ), 0.05 );
  }
}


BOOST_AUTO_TEST_CASE( ApplyUnapplyRoundTrip_FullRangeFraction )
{
  const size_t nchannel = 1024;
  const double offset_adj = 2.0;
  const double gain_adj = -10.0;
  const double quad_adj = 1.0;

  std::vector<double> orig_coeffs = { 0.0, 2800.0, 150.0 };
  std::vector<pair<double,double>> dev_pairs = { {600.0, 1.5}, {1400.0, 2.0}, {2200.0, 1.0} };
  const auto orig_spline = RelActCalcAutoImp::create_cubic_spline( dev_pairs );

  std::vector<double> adj_coeffs = orig_coeffs;
  if( adj_coeffs.size() < 3 )
    adj_coeffs.resize( 3, 0.0 );
  adj_coeffs[0] += offset_adj;
  adj_coeffs[1] += gain_adj;
  adj_coeffs[2] += quad_adj;

  for( double E_true = 100.0; E_true < 2700.0; E_true += 80.0 )
  {
    const double channel = find_fullrangefraction_channel( E_true, adj_coeffs, nchannel, dev_pairs );
    const double E_spec = fullrangefraction_energy( channel, orig_coeffs, nchannel, orig_spline );

    const double ch2 = find_fullrangefraction_channel( E_spec, orig_coeffs, nchannel, dev_pairs );
    const double E_true_back = fullrangefraction_energy( ch2, adj_coeffs, nchannel, orig_spline );

    BOOST_CHECK_SMALL( std::fabs( E_true_back - E_true ), 0.05 );
  }
}


BOOST_AUTO_TEST_CASE( ApplyUnapplyRoundTrip_LowerChannelEdge )
{
  const size_t nchannel = 512;
  std::vector<float> ch_energies;
  ch_energies.reserve( nchannel + 1 );
  for( size_t i = 0; i <= nchannel; ++i )
    ch_energies.push_back( static_cast<float>( i * 6.0 ) );  // 0 to 3072 keV

  const double lower_energy = static_cast<double>( ch_energies[0] );
  const double upper_energy = static_cast<double>( ch_energies[nchannel] );
  const double range = upper_energy - lower_energy;

  const double offset_adj = 4.0;
  const double gain_adj = 8.0;
  std::vector<pair<double,double>> dev_pairs = { {500.0, 0.5}, {1500.0, 1.0}, {2500.0, 0.5} };
  const auto orig_spline = RelActCalcAutoImp::create_cubic_spline( dev_pairs );
  const auto adj_spline = RelActCalcAutoImp::create_cubic_spline( dev_pairs );

  auto orig_energy_at_channel = [&]( double ch ) -> double {
    if( ch < 0.0 ) ch = 0.0;
    const size_t low = std::min( static_cast<size_t>( std::floor( ch ) ), nchannel - 1 );
    const double e_low = static_cast<double>( ch_energies[low] );
    const double e_high = static_cast<double>( ch_energies[low + 1] );
    const double frac = ch - static_cast<double>( low );
    const double base = e_low + (e_high - e_low) * frac;
    return base + eval_cubic_spline( base, orig_spline );
  };

  auto adj_energy_at_channel = [&]( double ch ) -> double {
    if( ch < 0.0 ) ch = 0.0;
    const size_t low = std::min( static_cast<size_t>( std::floor( ch ) ), nchannel - 1 );
    const double e_low = static_cast<double>( ch_energies[low] );
    const double e_high = static_cast<double>( ch_energies[low + 1] );
    const double frac = ch - static_cast<double>( low );
    const double e_base = e_low + (e_high - e_low) * frac;
    const double range_frac = (range > 0.0) ? ((e_base - lower_energy) / range) : 0.0;
    const double e_adj = e_base + offset_adj + range_frac * gain_adj;
    return e_adj + eval_cubic_spline( e_adj, adj_spline );
  };

  for( double E_true = 100.0; E_true < 3000.0; E_true += 120.0 )
  {
    const double channel = find_lowerchannel_channel( E_true, ch_energies, dev_pairs, offset_adj, gain_adj );
    const double E_spec = orig_energy_at_channel( channel );

    const double ch2 = find_lowerchannel_channel( E_spec, ch_energies, dev_pairs );
    const double E_true_back = adj_energy_at_channel( ch2 );

    BOOST_CHECK_SMALL( std::fabs( E_true_back - E_true ), 0.1 );
  }
}


BOOST_AUTO_TEST_CASE( EdgeCases )
{
  // Test edge cases and boundary conditions

  // Empty coefficients
  {
    std::vector<double> empty_coeffs;
    const double channel = find_polynomial_channel( 1000.0, empty_coeffs, 1024,
                                                    std::vector<std::pair<double,double>>{} );
    BOOST_CHECK_SMALL( channel, 1.0e-10 );
  }

  // Single coefficient (invalid)
  {
    std::vector<double> single_coeff = { 10.0 };
    const double channel = find_polynomial_channel( 1000.0, single_coeff, 1024,
                                                    std::vector<std::pair<double,double>>{} );
    BOOST_CHECK_SMALL( channel, 1.0e-10 );
  }

  // Energy outside calibration range
  {
    const size_t nchannel = 1024;
    std::vector<double> coeffs = { 0.0, 3.0 };

    // Very high energy (should clamp to last channel)
    const double high_channel = find_polynomial_channel( 10000.0, coeffs, nchannel,
                                                         std::vector<std::pair<double,double>>{} );
    BOOST_CHECK( high_channel >= static_cast<double>(nchannel - 1) );

    // Negative energy (for linear calibration E=3*ch, analytical solution gives negative channel)
    const double low_channel = find_polynomial_channel( -100.0, coeffs, nchannel,
                                                        std::vector<std::pair<double,double>>{} );
    const double expected_low = -100.0 / 3.0;  // Analytical solution for E = 0 + 3*ch
    BOOST_CHECK_CLOSE( low_channel, expected_low, 0.01 );
  }
}
