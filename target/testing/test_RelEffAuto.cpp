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

#include <string>
#include <iostream>
#include <atomic>
#include <cmath>
#include <limits>

#include <ceres/jet.h>

#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/RelActCalcAuto_imp.hpp"
#include "InterSpec/RelActCalc_imp.hpp"

#include "rapidxml/rapidxml.hpp"

//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE RelEffManual_suite
//#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/RelActCalcManual.h"
#include "InterSpec/DecayDataBaseServer.h"

#include "InterSpec/PeakFit_imp.hpp"


using namespace std;
using namespace boost::unit_test;

namespace RelActCalcAutoImp
{
size_t combine_nearby_peaks( std::vector<RelActCalcManual::GenericPeakInfo> &input_peaks,
                             double combine_threshold );
}

static_assert( USE_REL_ACT_TOOL, "Compile-time support for Rel Act tool is required." );


BOOST_AUTO_TEST_CASE( empirical_exponential_continuation_is_finite_monotone_and_c1 )
{
  using Jet = ceres::Jet<double,1>;
  const double values[] = { -100.0, -51.0, -50.001, -50.0, -49.999,
                              0.0,  49.999,  50.0,  50.001,  51.0, 100.0 };
  std::atomic<size_t> occupied{0};
  double previous = -1.0;
  for( const double value : values )
  {
    Jet input(value);
    input.v[0] = 1.0;
    const Jet output = RelActCalc::continued_exp( input, 50.0, &occupied );
    BOOST_CHECK( std::isfinite(output.a) );
    BOOST_CHECK( std::isfinite(output.v[0]) );
    BOOST_CHECK_GT( output.a, 0.0 );
    BOOST_CHECK_GT( output.v[0], 0.0 );
    if( previous >= 0.0 )
      BOOST_CHECK_GT( output.a, previous );
    previous = output.a;
  }
  BOOST_CHECK_EQUAL( occupied.load(), 6u );

  for( const double join : {-50.0, 50.0} )
  {
    Jet at_join(join);
    at_join.v[0] = 1.0;
    const Jet joined = RelActCalc::continued_exp( at_join, 50.0 );
    BOOST_CHECK_CLOSE_FRACTION( joined.a, std::exp(join), 2.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( joined.v[0], std::exp(join), 2.0e-15 );

    const double tail_side = std::nextafter(
        join, (join < 0.0) ? -std::numeric_limits<double>::infinity()
                           :  std::numeric_limits<double>::infinity() );
    Jet just_in_tail(tail_side);
    just_in_tail.v[0] = 1.0;
    const Jet continued = RelActCalc::continued_exp( just_in_tail, 50.0 );
    BOOST_CHECK_CLOSE_FRACTION( continued.a, joined.a, 3.0e-14 );
    BOOST_CHECK_CLOSE_FRACTION( continued.v[0], joined.v[0], 3.0e-14 );
  }
}


BOOST_AUTO_TEST_CASE( empirical_exponential_continuation_handles_extreme_finite_inputs )
{
  using Jet = ceres::Jet<double,1>;
  const double largest = std::numeric_limits<double>::max();
  const double inputs[] = {
    -largest, -0.5*largest, -1.0e300, -1.0e150, -100.0, -50.001,
    -50.0, 0.0, 50.0, 50.001, 100.0, 1.0e150, 1.0e300,
    0.5*largest, largest
  };

  double previous_value = -1.0;
  for( const double input_value : inputs )
  {
    Jet input(input_value);
    input.v[0] = 1.0;
    const Jet output = RelActCalc::continued_exp( input, 50.0 );
    const double scalar_output = RelActCalc::continued_exp( input_value, 50.0 );

    BOOST_CHECK( std::isfinite(output.a) );
    BOOST_CHECK( std::isfinite(output.v[0]) );
    BOOST_CHECK( std::isfinite(scalar_output) );
    BOOST_CHECK_GT( output.a, 0.0 );
    BOOST_CHECK_GT( scalar_output, 0.0 );
    BOOST_CHECK_GE( output.v[0], 0.0 );
    BOOST_CHECK_EQUAL( output.a, scalar_output );
    if( previous_value >= 0.0 )
      BOOST_CHECK_GE( output.a, previous_value );
    previous_value = output.a;
  }

  // The upper logarithmic overflow branch retains a useful recovery direction even at DBL_MAX.
  Jet upper(largest);
  upper.v[0] = 1.0;
  const Jet upper_output = RelActCalc::continued_exp( upper, 50.0 );
  BOOST_CHECK_GT( upper_output.v[0], 0.0 );

  // The reciprocal derivative remains positive throughout every representable, practically
  // relevant scale.  At still larger magnitudes it may underflow to zero, while the value remains
  // at its positive normal asymptote (a positive double slope over the entire DBL_MAX-wide lower
  // interval is numerically impossible without exhausting the available value range).
  Jet lower(-1.0e150);
  lower.v[0] = 1.0;
  const Jet lower_output = RelActCalc::continued_exp( lower, 50.0 );
  BOOST_CHECK_GT( lower_output.v[0], 0.0 );
  BOOST_CHECK_GE( RelActCalc::continued_exp(-largest, 50.0),
                  std::numeric_limits<double>::min() );

  // Values in the uncontinued region are exactly the native exponential, not merely close to it.
  for( const double input_value : {-49.0, -1.0, 0.0, 1.0, 49.0} )
    BOOST_CHECK_EQUAL( RelActCalc::continued_exp(input_value,50.0), std::exp(input_value) );
}


BOOST_AUTO_TEST_CASE( fwhm_continuation_is_positive_monotone_and_c1 )
{
  using Jet = ceres::Jet<double,1>;
  const double channel_width = 0.25;
  const double lower = 1.25*channel_width;
  const double join = 2.0*lower;
  double previous = -1.0;
  for( const double raw : {-100.0,-1.0,0.0,0.5*join,join-1.0e-8,join,join+0.1,5.0} )
  {
    Jet input(raw);
    input.v[0] = 1.0;
    const Jet output = RelActCalcAutoImp::resolvable_fwhm_continuation(input,channel_width);
    BOOST_CHECK( std::isfinite(output.a) );
    BOOST_CHECK( std::isfinite(output.v[0]) );
    BOOST_CHECK_GE( output.a, lower );
    BOOST_CHECK_GT( output.v[0], 0.0 );
    if( previous >= 0.0 )
      BOOST_CHECK_GT( output.a, previous );
    previous = output.a;
  }

  Jet at_join(join);
  at_join.v[0] = 1.0;
  const Jet joined = RelActCalcAutoImp::resolvable_fwhm_continuation(at_join,channel_width);
  BOOST_CHECK_CLOSE_FRACTION( joined.a,join,2.0e-15 );
  BOOST_CHECK_CLOSE_FRACTION( joined.v[0],1.0,2.0e-15 );
}


BOOST_AUTO_TEST_CASE( frozen_gamma_membership_covers_calibration_motion_bounds )
{
  using RelActCalcAutoImp::frozen_gamma_calibration_motion_guard_keV;
  using RelActCalcAutoImp::frozen_gamma_membership_energy_limits;

  BOOST_CHECK_EQUAL( frozen_gamma_calibration_motion_guard_keV(
                       false, false, false, 30.0, 12.5), 0.0 );

  const double linear_one_point
    = RelActCalcAuto::RelActAutoSolution::sm_energy_offset_range_keV
      + RelActCalcAuto::RelActAutoSolution::sm_energy_gain_range_keV
      + ((RelActCalcAuto::RelActAutoSolution::sm_num_energy_cal_pars > 2)
         ? RelActCalcAuto::RelActAutoSolution::sm_energy_quad_range_keV : 0.0);
  const double linear_guard = frozen_gamma_calibration_motion_guard_keV(
                                true, false, false, linear_one_point, 12.5 );
  BOOST_CHECK_EQUAL( linear_guard, linear_one_point );
  BOOST_CHECK_EQUAL( frozen_gamma_calibration_motion_guard_keV(
                       true, false, true, linear_one_point, 12.5 ),
                     2.0*linear_one_point );

  const double deviation_l1_budget = 7.5;
  const double nonlinear_guard = frozen_gamma_calibration_motion_guard_keV(
                                   true, true, false, linear_one_point,
                                   deviation_l1_budget );
  BOOST_CHECK_EQUAL( nonlinear_guard,
                     linear_one_point + deviation_l1_budget );

  // A line just beyond the peak-shape envelope is admitted through the complete feasible
  // calibration movement, while one beyond that closed guard is not.
  const double roi_lower = 610.0;
  const double roi_upper = 775.0;
  const double shape_margin = 15.0;
  const double admitted_line = roi_lower - shape_margin - nonlinear_guard + 0.5e-6;
  const double excluded_line = roi_lower - shape_margin - nonlinear_guard - 2.0e-6;
  const std::pair<double,double> frozen_limits = frozen_gamma_membership_energy_limits(
      roi_lower,roi_upper,shape_margin,shape_margin,nonlinear_guard );
  BOOST_CHECK_GE( admitted_line, frozen_limits.first );
  BOOST_CHECK_LT( excluded_line, frozen_limits.first );
  BOOST_CHECK_LE( roi_upper + shape_margin + nonlinear_guard, frozen_limits.second );
}

std::string g_test_file_dir;

// We need to set the static data directory, so the code knows where
//  like sandia.decay.xml is located.
void set_data_dir()
{
  // We only need to initialize things once
  static bool s_have_set = false;
  if( s_have_set )
    return;
  
  s_have_set = true;
  
  int argc = boost::unit_test::framework::master_test_suite().argc;
  char **argv = boost::unit_test::framework::master_test_suite().argv;
  
  string datadir;
  
  for( int i = 1; i < argc; ++i )
  {
    cout << "Arg " << i << ": '" << argv[i] << "'" << endl;
    const string arg = argv[i];
    if( SpecUtils::istarts_with( arg, "--datadir=" ) )
      datadir = arg.substr( 10 );
    
    if( SpecUtils::istarts_with( arg, "--testfiledir=" ) )
      g_test_file_dir = arg.substr( 14 );
  }//for( int arg = 1; arg < argc; ++ arg )
  
  SpecUtils::ireplace_all( datadir, "%20", " " );
  SpecUtils::ireplace_all( g_test_file_dir, "%20", " " );
  
  // Search around a little for the data directory, if it wasnt specified
  if( datadir.empty() )
  {
    for( const auto &d : { "data", "../data", "../../data", "../../../data", "/Users/wcjohns/rad_ana/InterSpec/data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) )
      {
        datadir = d;
        break;
      }
    }//for( loop over candidate dirs )
  }//if( datadir.empty() )
  
  const string required_data_file = "findCharacteristics/202204_example_problem_1.n42";
  if( g_test_file_dir.empty() )
  {
    for( const auto &d : { "test_data", "../test_data", "../../test_data", "/Users/wcjohns/rad_ana/InterSpec/target/testing/test_data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, required_data_file) ) )
      {
        g_test_file_dir = d;
        break;
      }
    }//for( loop over candidate dirs )
  }
  
  const string sandia_deacay_file = SpecUtils::append_path(datadir, "sandia.decay.xml");
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( sandia_deacay_file ), "sandia.decay.xml not at '" << sandia_deacay_file << "'" );
  
  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );
  
  // Make sure we can actually init the decay database
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );
  BOOST_REQUIRE_MESSAGE( db->nuclide("U238"), "SandiaDecayDataBase empty?" );
  
  const string required_data_path = SpecUtils::append_path(g_test_file_dir, required_data_file);
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( required_data_path ), "'" << required_data_file << "' not at '" << required_data_path << "'" );
}//void set_data_dir()



BOOST_AUTO_TEST_CASE( SkewBoundsAndNoSkewValue )
{
  // Review item A12: GaussExp/ExpGaussExp upper skew bound widened 3.25 -> 4.0 and CrystalBall/DSCB
  //  alpha 4.0 -> 5.0, plus a new PeakDef::skew_no_skew_value() giving the value to fix a skew
  //  coefficient at to remove skew.  The "no skew" end differs per model: GaussExp/CrystalBall reach
  //  it only at the (asymptotic) upper bound, Bortel/GaussPlusBortel at their lower bound.
  double lo = 0, hi = 0, start = 0, step = 0, nsv = 0;

  // --- Widened bounds ---
  BOOST_REQUIRE( PeakDef::skew_parameter_range( PeakDef::SkewType::GaussExp, PeakDef::CoefficientType::SkewPar0, lo, hi, start, step ) );
  BOOST_CHECK_CLOSE( hi, 4.0, 1.0E-9 );
  BOOST_CHECK_CLOSE( lo, 0.15, 1.0E-9 );

  BOOST_REQUIRE( PeakDef::skew_parameter_range( PeakDef::SkewType::ExpGaussExp, PeakDef::CoefficientType::SkewPar1, lo, hi, start, step ) );
  BOOST_CHECK_CLOSE( hi, 4.0, 1.0E-9 );

  BOOST_REQUIRE( PeakDef::skew_parameter_range( PeakDef::SkewType::CrystalBall, PeakDef::CoefficientType::SkewPar0, lo, hi, start, step ) );
  BOOST_CHECK_CLOSE( hi, 5.0, 1.0E-9 );  // alpha (left)
  BOOST_REQUIRE( PeakDef::skew_parameter_range( PeakDef::SkewType::DoubleSidedCrystalBall, PeakDef::CoefficientType::SkewPar2, lo, hi, start, step ) );
  BOOST_CHECK_CLOSE( hi, 5.0, 1.0E-9 );  // alpha (right)

  // --- skew_no_skew_value: GaussExp/ExpGaussExp/CrystalBall -> upper bound; Bortel/GaussPlusBortel -> 0 ---
  BOOST_CHECK( PeakDef::skew_no_skew_value( PeakDef::SkewType::GaussExp, PeakDef::CoefficientType::SkewPar0, nsv ) );
  BOOST_CHECK_CLOSE( nsv, 4.0, 1.0E-9 );

  BOOST_CHECK( PeakDef::skew_no_skew_value( PeakDef::SkewType::ExpGaussExp, PeakDef::CoefficientType::SkewPar1, nsv ) );
  BOOST_CHECK_CLOSE( nsv, 4.0, 1.0E-9 );

  BOOST_CHECK( PeakDef::skew_no_skew_value( PeakDef::SkewType::CrystalBall, PeakDef::CoefficientType::SkewPar0, nsv ) );
  BOOST_CHECK_CLOSE( nsv, 5.0, 1.0E-9 );  // alpha -> upper bound

  BOOST_CHECK( PeakDef::skew_no_skew_value( PeakDef::SkewType::Bortel, PeakDef::CoefficientType::SkewPar0, nsv ) );
  BOOST_CHECK_SMALL( nsv, 1.0E-12 );      // tau -> 0 is an exact pure Gaussian

  BOOST_CHECK( PeakDef::skew_no_skew_value( PeakDef::SkewType::GaussPlusBortel, PeakDef::CoefficientType::SkewPar0, nsv ) );
  BOOST_CHECK_SMALL( nsv, 1.0E-12 );      // R = 0 is a pure Gaussian

  // CrystalBall power-law `n` is a "don't care" once alpha is at no-skew: returns true at its neutral start.
  double n_lo = 0, n_hi = 0, n_start = 0, n_step = 0;
  BOOST_REQUIRE( PeakDef::skew_parameter_range( PeakDef::SkewType::CrystalBall, PeakDef::CoefficientType::SkewPar1, n_lo, n_hi, n_start, n_step ) );
  BOOST_CHECK( PeakDef::skew_no_skew_value( PeakDef::SkewType::CrystalBall, PeakDef::CoefficientType::SkewPar1, nsv ) );
  BOOST_CHECK_CLOSE( nsv, n_start, 1.0E-9 );

  // --- No removal offered: DoubleBortel (no pure-Gaussian limit), NoSkew, inapplicable coefficient ---
  BOOST_CHECK( !PeakDef::skew_no_skew_value( PeakDef::SkewType::DoubleBortel, PeakDef::CoefficientType::SkewPar0, nsv ) );
  BOOST_CHECK( !PeakDef::skew_no_skew_value( PeakDef::SkewType::NoSkew, PeakDef::CoefficientType::SkewPar0, nsv ) );
  BOOST_CHECK( !PeakDef::skew_no_skew_value( PeakDef::SkewType::GaussExp, PeakDef::CoefficientType::SkewPar1, nsv ) ); // GaussExp uses only SkewPar0
}//BOOST_AUTO_TEST_CASE( SkewBoundsAndNoSkewValue )


BOOST_AUTO_TEST_CASE( FitContinuum )
{
  //set_data_dir();
  
  
  vector<PeakDef> peaks;
  peaks.emplace_back( 103.5, 2.5, 150.0 );

  constexpr size_t nbin = 7;
  constexpr float energies[nbin+1] = {100.0f, 101.0f, 102.0f, 103.0f, 104.0f, 105.0f, 106.0f, 107.0f};
  constexpr float data[nbin] = {900.0f, 1090.0f, 990.0f, 1090.0f, 910.0f, 1090.0f, 950.0f};
  const PeakContinuum::OffsetType offset_type = PeakContinuum::OffsetType::Linear;
  const size_t num_polynomial_terms = PeakContinuum::num_parameters( offset_type );
  constexpr double ref_energy = energies[0];
  vector<double> continuum_coeffs(num_polynomial_terms, 0.0);
  double peak_counts[nbin];

  PeakFit::fit_continuum( energies, data, nullptr, nbin, offset_type,
                                  ref_energy, peaks, false, continuum_coeffs.data(), peak_counts );


  vector<double> dummy_amps, continuum_coeffs_old, dummy_amp_uncert, continuum_uncerts;
  fit_amp_and_offset( energies, data, nbin, offset_type,
                      ref_energy, {}, {}, peaks,
                      PeakDef::SkewType::NoSkew, nullptr, dummy_amps,
                      continuum_coeffs_old, dummy_amp_uncert, continuum_uncerts );

  BOOST_REQUIRE( continuum_coeffs_old.size() == num_polynomial_terms );
  for( size_t i = 0; i < num_polynomial_terms; ++i )
  {
    const double new_coef = continuum_coeffs[i];
    const double old_coef = continuum_coeffs_old[i];
    const double diff = fabs(new_coef - old_coef);
    BOOST_CHECK( (diff < 0.00001*(std::max)(fabs(new_coef), fabs(old_coef))) || (diff < 1.0E-12) );
  }
  
  double old_way_peak_counts[nbin] = { 0.0 };
  
  vector<const PeakDef *> pk_ptrs;
  for( const PeakDef &peak : peaks )
  {
    peak.gauss_integral( energies, old_way_peak_counts, nbin );
    pk_ptrs.push_back( &peak );
  }

  std::shared_ptr<PeakContinuum> continuum = peaks[0].continuum();
  BOOST_REQUIRE( !!continuum );
  continuum->setType( offset_type );
  continuum->setParameters( ref_energy, continuum_coeffs, {} );

  for( size_t bin = 0; bin < nbin; ++bin )
    {
      std::shared_ptr<const SpecUtils::Measurement> null_data;
      old_way_peak_counts[bin] += continuum->offset_integral( energies[bin], energies[bin+1], null_data, pk_ptrs.data(), pk_ptrs.size() );
    }

  for( size_t bin = 0; bin < nbin; ++bin )
  {
    const double new_val = peak_counts[bin];
    const double old_coef = old_way_peak_counts[bin];
    const double diff = fabs(new_val - old_coef);
    BOOST_CHECK( (diff < 0.00001*(std::max)(fabs(new_val), fabs(old_coef))) || (diff < 1.0E-12) );
  }

}//BOOST_AUTO_TEST_CASE( FitRelActManualToKnown )


BOOST_AUTO_TEST_CASE( test_pu242_by_correlation )
{
  // We will first roughly test PuCorrMethod::ByPu239Only to data given in paper.
  //
  // Fig 3 in Swinhoe 2010 gives Pu239 content vs Pu242 content; I manually
  //  extracted the following values from the fit line in the PDF.
  const vector<pair<double,double>> swinhoe_approx_fig_3_data = {
    {0.55496, 0.06998},
    {0.55894, 0.06821},
    {0.56438, 0.06619},
    {0.57073, 0.06399},
    {0.57829, 0.06154},
    {0.58453, 0.05952},
    {0.59068, 0.05756},
    {0.59471, 0.05634},
    {0.59844, 0.05524},
    {0.60146, 0.05444},
    {0.60579, 0.05322},
    {0.61526, 0.05077},
    {0.62101, 0.04906},
    {0.62605, 0.04802},
    {0.63088, 0.04691},
    {0.63542, 0.04594},
    {0.6402, 0.0449}
  };//swinhoe_approx_fig_3_data
  
  
  for( const auto x_y : swinhoe_approx_fig_3_data )
  {
    const double x = x_y.first;
    const double y = x_y.second;
    const double gamma_spec_pu239 = x / (1.0 - y);
    const double gamma_spec_pu_other = (1.0 - x - y)/(1.0 - y);
    
    // gamma_spec_pu239 plus gamma_spec_pu_other should sum to 1.0
    BOOST_CHECK( fabs(1.0 - (gamma_spec_pu239 + gamma_spec_pu_other)) < 0.001 );
    
    RelActCalc::Pu242ByCorrelationInput<double> input;
    input.pu_age = 0.0;
    input.pu238_rel_mass = gamma_spec_pu_other;
    input.pu239_rel_mass = gamma_spec_pu239;
    // Pu240, and Pu241/Am241 are irrelevant, all that
    
    RelActCalc::Pu242ByCorrelationOutput<double> output = RelActCalc::correct_pu_mass_fractions_for_pu242( input, RelActCalc::PuCorrMethod::ByPu239Only );
    
    //cout << "For Swinhoe [" << x << ", " << y << "]: Pu239: " << output.pu239_mass_frac
    //     << ", Pu242: " << output.pu242_mass_frac << " +- " << 100.0*output.pu242_uncert << "%\n";
    
    BOOST_CHECK( fabs(output.pu239_mass_frac - x) < 0.005 );
    BOOST_CHECK( fabs(output.pu242_mass_frac - y) < 0.0005 );
  }//for( const auto x_y : swinhoe_approx_fig_3_data )
  
  
  // For PuCorrMethod::Bignan95_BWR and PuCorrMethod::Bignan95_PWR, we dont have nearly as good
  //  of comparison data
  RelActCalc::Pu242ByCorrelationInput<double> input;
  input.pu_age = 0.0;
  input.pu238_rel_mass = 0.0120424;
  input.pu239_rel_mass = 0.6649628;
  input.pu240_rel_mass = 0.2327493;
  input.pu241_rel_mass = 0.0501864;
  //input.pu241_rel_mass = 0.0361259;
  RelActCalc::Pu242ByCorrelationOutput<double> output = RelActCalc::correct_pu_mass_fractions_for_pu242( input, RelActCalc::PuCorrMethod::Bignan95_BWR );
  //cout << "For Bignan95_BWR: Pu239: " << output.pu239_mass_frac
  //     << ", Pu242: " << output.pu242_mass_frac << " +- " << 100.0*output.pu242_uncert << "%\n";
  //For Bignan95_BWR: Pu239: 0.668767, Pu242: 0.0345679 +- 14%
  BOOST_CHECK( fabs(output.pu239_mass_frac - 0.668767) < 0.0001 );
  BOOST_CHECK( fabs(output.pu242_mass_frac - 0.0345679) < 0.0001 );
  BOOST_CHECK( fabs(output.pu242_uncert - 0.14) < 0.0001 );
  
  output = RelActCalc::correct_pu_mass_fractions_for_pu242( input, RelActCalc::PuCorrMethod::Bignan95_PWR );
  //cout << "For Bignan95_PWR: Pu239: " << output.pu239_mass_frac
  //     << ", Pu242: " << output.pu242_mass_frac << " +- " << 100.0*output.pu242_uncert << "%\n";
  //For Bignan95_PWR: Pu239: 0.664565, Pu242: 0.0406335 +- 6%
  BOOST_CHECK( fabs(output.pu239_mass_frac - 0.664565) < 0.0001 );
  BOOST_CHECK( fabs(output.pu242_mass_frac - 0.0406335) < 0.0001 );
  BOOST_CHECK( fabs(output.pu242_uncert - 0.06) < 0.0001 );
}//BOOST_AUTO_TEST_CASE( test_pu242_by_correlation )


BOOST_AUTO_TEST_CASE( profile_safe_pu242_correlation_is_finite_physical_and_c1 )
{
  using Jet = ceres::Jet<double,1>;
  using DoubleOutput = RelActCalc::Pu242ByCorrelationOutput<double>;

  const auto check_physical = []( const DoubleOutput &output ) {
    const double values[] = { output.pu238_mass_frac,output.pu239_mass_frac,
                              output.pu240_mass_frac,output.pu241_mass_frac,
                              output.pu242_mass_frac };
    double sum = 0.0;
    for( const double value : values )
    {
      BOOST_CHECK( std::isfinite(value) );
      BOOST_CHECK_GE( value,0.0 );
      BOOST_CHECK_LE( value,1.0 );
      sum += value;
    }
    BOOST_CHECK_SMALL( sum-1.0,2.0E-14 );
  };

  // The closed activity box contains an all-zero corner.  It is not a defined nominal
  // composition, but the optimizer-only singular-corner guard must give every correlation method
  // a finite physical residual there.
  for( const RelActCalc::PuCorrMethod method : {
       RelActCalc::PuCorrMethod::ByPu239Only,RelActCalc::PuCorrMethod::Bignan95_PWR,
       RelActCalc::PuCorrMethod::Bignan95_BWR } )
  {
    RelActCalc::Pu242ByCorrelationInput<double> zero;
    BOOST_CHECK_THROW(
        RelActCalc::correct_pu_mass_fractions_for_pu242(zero,method),std::domain_error );
    const DoubleOutput output
        = RelActCalc::correct_pu_mass_fractions_for_pu242(zero,method,true);
    check_physical(output);
    BOOST_CHECK_EQUAL( output.pu238_mass_frac,0.0 );
    BOOST_CHECK_EQUAL( output.pu239_mass_frac,0.0 );
    BOOST_CHECK_EQUAL( output.pu240_mass_frac,0.0 );
    BOOST_CHECK_EQUAL( output.pu241_mass_frac,0.0 );
    BOOST_CHECK_EQUAL( output.pu242_mass_frac,1.0 );
  }

  RelActCalc::Pu242ByCorrelationInput<double> no_pu239;
  no_pu239.pu238_rel_mass = 0.1;
  no_pu239.pu240_rel_mass = 0.9;
  BOOST_CHECK_THROW(
      RelActCalc::correct_pu_mass_fractions_for_pu242(
          no_pu239,RelActCalc::PuCorrMethod::ByPu239Only),std::domain_error );
  const DoubleOutput no_pu239_output
      = RelActCalc::correct_pu_mass_fractions_for_pu242(
          no_pu239,RelActCalc::PuCorrMethod::ByPu239Only,true);
  check_physical(no_pu239_output);
  BOOST_CHECK_EQUAL( no_pu239_output.pu242_mass_frac,1.0 );

  RelActCalc::Pu242ByCorrelationInput<double> unrepresentable;
  unrepresentable.pu239_rel_mass = 0.9;
  unrepresentable.other_pu_mass = 0.1;
  BOOST_CHECK_THROW(
      RelActCalc::correct_pu_mass_fractions_for_pu242(
          unrepresentable,RelActCalc::PuCorrMethod::ByPu239Only),std::domain_error );

  RelActCalc::Pu242ByCorrelationInput<double> nonfinite;
  nonfinite.pu239_rel_mass = std::numeric_limits<double>::infinity();
  BOOST_CHECK_THROW(
      RelActCalc::correct_pu_mass_fractions_for_pu242(
          nonfinite,RelActCalc::PuCorrMethod::ByPu239Only),std::domain_error );
  check_physical( RelActCalc::correct_pu_mass_fractions_for_pu242(
      nonfinite,RelActCalc::PuCorrMethod::ByPu239Only,true) );

  RelActCalc::Pu242ByCorrelationInput<Jet> zero_jet;
  const auto zero_jet_output = RelActCalc::correct_pu_mass_fractions_for_pu242(
      zero_jet,RelActCalc::PuCorrMethod::ByPu239Only,true);
  for( const Jet * const value : { &zero_jet_output.pu238_mass_frac,
                                  &zero_jet_output.pu239_mass_frac,
                                  &zero_jet_output.pu240_mass_frac,
                                  &zero_jet_output.pu241_mass_frac,
                                  &zero_jet_output.pu242_mass_frac } )
  {
    BOOST_CHECK( std::isfinite(value->a) );
    BOOST_CHECK( std::isfinite(value->v[0]) );
  }
  BOOST_CHECK_EQUAL( zero_jet_output.pu242_mass_frac.a,1.0 );
  BOOST_CHECK_EQUAL( zero_jet_output.pu242_mass_frac.v[0],0.0 );

  // Pu93-like ordinary input is well below the continuation seam.  Both scalar values and Jet
  // derivatives must take the literal historical branch, so nominal/profile coordinates agree
  // exactly and the production change cannot perturb the required numeric Pu-242 extrapolation.
  RelActCalc::Pu242ByCorrelationInput<double> pu93;
  pu93.pu238_rel_mass = 0.000091612433;
  pu93.pu239_rel_mass = 0.935920939673;
  pu93.pu240_rel_mass = 0.063100423177;
  pu93.pu241_rel_mass = 0.000490931012;
  const DoubleOutput pu93_legacy = RelActCalc::correct_pu_mass_fractions_for_pu242(
      pu93,RelActCalc::PuCorrMethod::ByPu239Only);
  const DoubleOutput pu93_safe = RelActCalc::correct_pu_mass_fractions_for_pu242(
      pu93,RelActCalc::PuCorrMethod::ByPu239Only,true);
  const double pu93_expected_pu242 = 9.66E-3*std::pow(
      pu93.pu239_rel_mass/(pu93.pu238_rel_mass+pu93.pu239_rel_mass
                            +pu93.pu240_rel_mass+pu93.pu241_rel_mass),-3.83);
  BOOST_CHECK_EQUAL( pu93_legacy.pu242_mass_frac,pu93_expected_pu242 );
  BOOST_CHECK_EQUAL( pu93_safe.pu238_mass_frac,pu93_legacy.pu238_mass_frac );
  BOOST_CHECK_EQUAL( pu93_safe.pu239_mass_frac,pu93_legacy.pu239_mass_frac );
  BOOST_CHECK_EQUAL( pu93_safe.pu240_mass_frac,pu93_legacy.pu240_mass_frac );
  BOOST_CHECK_EQUAL( pu93_safe.pu241_mass_frac,pu93_legacy.pu241_mass_frac );
  BOOST_CHECK_EQUAL( pu93_safe.pu242_mass_frac,pu93_legacy.pu242_mass_frac );
  BOOST_CHECK_EQUAL( pu93_safe.pu242_uncert,pu93_legacy.pu242_uncert );
  BOOST_CHECK_EQUAL( pu93_safe.is_within_range,pu93_legacy.is_within_range );

  RelActCalc::Pu242ByCorrelationInput<Jet> pu93_jet;
  pu93_jet.pu238_rel_mass = Jet(pu93.pu238_rel_mass);
  pu93_jet.pu239_rel_mass = Jet(pu93.pu239_rel_mass);
  pu93_jet.pu240_rel_mass = Jet(pu93.pu240_rel_mass);
  pu93_jet.pu241_rel_mass = Jet(pu93.pu241_rel_mass);
  pu93_jet.pu239_rel_mass.v[0] = 1.0;
  const auto pu93_legacy_jet = RelActCalc::correct_pu_mass_fractions_for_pu242(
      pu93_jet,RelActCalc::PuCorrMethod::ByPu239Only);
  const auto pu93_safe_jet = RelActCalc::correct_pu_mass_fractions_for_pu242(
      pu93_jet,RelActCalc::PuCorrMethod::ByPu239Only,true);
  BOOST_CHECK_EQUAL( pu93_safe_jet.pu239_mass_frac.a,pu93_legacy_jet.pu239_mass_frac.a );
  BOOST_CHECK_EQUAL( pu93_safe_jet.pu239_mass_frac.v[0],pu93_legacy_jet.pu239_mass_frac.v[0] );
  BOOST_CHECK_EQUAL( pu93_safe_jet.pu242_mass_frac.a,pu93_legacy_jet.pu242_mass_frac.a );
  BOOST_CHECK_EQUAL( pu93_safe_jet.pu242_mass_frac.v[0],pu93_legacy_jet.pu242_mass_frac.v[0] );

  // At q=0.9 the exact negative-power branch joins the bounded log-tail.  Vary Pu239 against Pu240
  // so total mass stays one and compare derivatives from both sides.
  static constexpr double A = 9.66E-3;
  static constexpr double C = -3.83;
  static constexpr double q_join = 0.9;
  const double x_join = std::pow(q_join/A,1.0/C);
  const auto at_pu239 = []( const double x, const bool profile_safe = false ) {
    RelActCalc::Pu242ByCorrelationInput<Jet> input;
    input.pu239_rel_mass = Jet(x);
    input.pu239_rel_mass.v[0] = 1.0;
    input.pu240_rel_mass = Jet(1.0-x);
    input.pu240_rel_mass.v[0] = -1.0;
    return RelActCalc::correct_pu_mass_fractions_for_pu242(
        input,RelActCalc::PuCorrMethod::ByPu239Only,profile_safe).pu242_mass_frac;
  };
  const Jet below = at_pu239(x_join*(1.0-1.0E-8));
  const Jet seam = at_pu239(x_join);
  const Jet above = at_pu239(x_join*(1.0+1.0E-8));
  const Jet profile_seam = at_pu239(x_join,true);
  const double expected_seam_derivative = C*q_join/x_join;
  BOOST_CHECK_CLOSE_FRACTION( seam.a,q_join,2.0E-14 );
  BOOST_CHECK_CLOSE_FRACTION( seam.v[0],expected_seam_derivative,2.0E-14 );
  BOOST_CHECK_CLOSE_FRACTION( below.v[0],seam.v[0],5.0E-7 );
  BOOST_CHECK_CLOSE_FRACTION( above.v[0],seam.v[0],5.0E-7 );
  BOOST_CHECK_EQUAL( profile_seam.a,seam.a );
  BOOST_CHECK_EQUAL( profile_seam.v[0],seam.v[0] );
  check_physical( RelActCalc::correct_pu_mass_fractions_for_pu242(
      RelActCalc::Pu242ByCorrelationInput<double>{0.0,x_join*0.5,
                                                  1.0-x_join*0.5,0.0,0.0,0.0},
      RelActCalc::PuCorrMethod::ByPu239Only) );

  // The same bounded tail is production behavior for finite positive Bignan inputs.  Exercise a
  // gross extrapolation and the q=0.9 seam without enabling the optimizer-only singular guards.
  RelActCalc::Pu242ByCorrelationInput<double> bignan_tail;
  bignan_tail.pu238_rel_mass = 0.20;
  bignan_tail.pu239_rel_mass = 0.01;
  bignan_tail.pu240_rel_mass = 0.79;
  check_physical( RelActCalc::correct_pu_mass_fractions_for_pu242(
      bignan_tail,RelActCalc::PuCorrMethod::Bignan95_PWR) );

  static constexpr double bignan_c0 = 1.313;
  static constexpr double bignan_pu238 = 0.01;
  static constexpr double bignan_pu240 = 0.20;
  const double bignan_scale = bignan_c0*std::pow(bignan_pu238,0.33)
                            * std::pow(bignan_pu240,1.7);
  const double bignan_x_join = std::pow(bignan_scale/q_join,1.0/1.03);
  const auto at_bignan_pu239 = []( const double x ) {
    RelActCalc::Pu242ByCorrelationInput<Jet> input;
    input.pu238_rel_mass = Jet(bignan_pu238);
    input.pu239_rel_mass = Jet(x);
    input.pu239_rel_mass.v[0] = 1.0;
    input.pu240_rel_mass = Jet(bignan_pu240);
    input.pu241_rel_mass = Jet(1.0-bignan_pu238-bignan_pu240-x);
    input.pu241_rel_mass.v[0] = -1.0;
    return RelActCalc::correct_pu_mass_fractions_for_pu242(
        input,RelActCalc::PuCorrMethod::Bignan95_PWR).pu242_mass_frac;
  };
  const Jet bignan_below = at_bignan_pu239(bignan_x_join*(1.0-1.0E-8));
  const Jet bignan_seam = at_bignan_pu239(bignan_x_join);
  const Jet bignan_above = at_bignan_pu239(bignan_x_join*(1.0+1.0E-8));
  const double bignan_expected_derivative = -1.03*q_join/bignan_x_join;
  BOOST_CHECK_CLOSE_FRACTION( bignan_seam.a,q_join,2.0E-14 );
  BOOST_CHECK_CLOSE_FRACTION( bignan_seam.v[0],bignan_expected_derivative,2.0E-13 );
  BOOST_CHECK_CLOSE_FRACTION( bignan_below.v[0],bignan_seam.v[0],5.0E-7 );
  BOOST_CHECK_CLOSE_FRACTION( bignan_above.v[0],bignan_seam.v[0],5.0E-7 );
}


BOOST_AUTO_TEST_CASE( pu242_correlation_range_is_classified_at_back_decayed_epoch )
{
  set_data_dir();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );
  const SandiaDecay::Nuclide * const pu238 = db->nuclide("Pu238");
  const SandiaDecay::Nuclide * const pu239 = db->nuclide("Pu239");
  const SandiaDecay::Nuclide * const pu240 = db->nuclide("Pu240");
  const SandiaDecay::Nuclide * const pu241 = db->nuclide("Pu241");
  BOOST_REQUIRE( pu238 && pu239 && pu240 && pu241 );

  const double age = 120.0 * 365.25 * 24.0 * 60.0 * 60.0;
  const auto aged_mass = [age]( const double at_epoch,
                                const SandiaDecay::Nuclide * const nuc ) {
    return at_epoch * std::exp(-age*nuc->decayConstant());
  };

  // At the correlation epoch Pu-239 is 50%, below the Pu239-only method's validated range.  After
  // 120 years of Pu-241 decay the returned acquisition-time Pu-239 fraction moves into 55--80%; it
  // must remain marked extrapolated and retain the out-of-range uncertainty bucket.
  RelActCalc::Pu242ByCorrelationInput<double> by_pu239;
  by_pu239.pu239_rel_mass = aged_mass(0.50,pu239);
  by_pu239.pu240_rel_mass = aged_mass(0.10,pu240);
  by_pu239.pu241_rel_mass = aged_mass(0.40,pu241);
  const double by_sum = by_pu239.pu239_rel_mass + by_pu239.pu240_rel_mass
                        + by_pu239.pu241_rel_mass;
  by_pu239.pu239_rel_mass /= by_sum;
  by_pu239.pu240_rel_mass /= by_sum;
  by_pu239.pu241_rel_mass /= by_sum;
  by_pu239.pu_age = age;
  const auto by_output = RelActCalc::correct_pu_mass_fractions_for_pu242(
      by_pu239,RelActCalc::PuCorrMethod::ByPu239Only );
  BOOST_CHECK_GE( by_output.pu239_mass_frac,0.55 );
  BOOST_CHECK_LE( by_output.pu239_mass_frac,0.80 );
  BOOST_CHECK( !by_output.is_within_range );
  BOOST_CHECK_CLOSE_FRACTION( by_output.pu242_uncert,0.10,1.0E-14 );

  // For Bignan, start just above the correlation-epoch Pu238/Pu239 ratio limit while keeping the
  // other two ratios valid.  Pu-238 decay moves all three acquisition-time ratios into their
  // published windows; provenance must still reflect the out-of-range epoch used by the formula.
  const double epoch_pu239 = 1.0/(1.0 + 0.04 + 0.35);
  RelActCalc::Pu242ByCorrelationInput<double> bignan;
  bignan.pu238_rel_mass = aged_mass(0.04*epoch_pu239,pu238);
  bignan.pu239_rel_mass = aged_mass(epoch_pu239,pu239);
  bignan.pu240_rel_mass = aged_mass(0.35*epoch_pu239,pu240);
  const double bignan_sum = bignan.pu238_rel_mass + bignan.pu239_rel_mass
                            + bignan.pu240_rel_mass;
  bignan.pu238_rel_mass /= bignan_sum;
  bignan.pu239_rel_mass /= bignan_sum;
  bignan.pu240_rel_mass /= bignan_sum;
  bignan.pu_age = age;
  const auto bignan_output = RelActCalc::correct_pu_mass_fractions_for_pu242(
      bignan,RelActCalc::PuCorrMethod::Bignan95_PWR );
  const double reported_238_239 = bignan_output.pu238_mass_frac/bignan_output.pu239_mass_frac;
  const double reported_240_239 = bignan_output.pu240_mass_frac/bignan_output.pu239_mass_frac;
  const double reported_242_239 = bignan_output.pu242_mass_frac/bignan_output.pu239_mass_frac;
  BOOST_CHECK_GE( reported_238_239,0.007851 );
  BOOST_CHECK_LE( reported_238_239,0.02952 );
  BOOST_CHECK_GE( reported_240_239,0.2688 );
  BOOST_CHECK_LE( reported_240_239,0.4586 );
  BOOST_CHECK_GE( reported_242_239,0.03323 );
  BOOST_CHECK_LE( reported_242_239,0.1152 );
  BOOST_CHECK( !bignan_output.is_within_range );
}


BOOST_AUTO_TEST_CASE( pu241_without_pu240_is_back_decayed_for_correlation )
{
  set_data_dir();

  // With no Pu-242 correction, back-decaying the input and then restoring it to the measurement
  // age must be an identity operation.  This specifically exercises Pu-241 with Pu-240 absent:
  // the old guard accidentally tested pu240_rel_mass and therefore aged Pu-241 only in the forward
  // direction a second time.
  RelActCalc::Pu242ByCorrelationInput<double> input;
  input.pu239_rel_mass = 0.9;
  input.pu240_rel_mass = 0.0;
  input.pu241_rel_mass = 0.1;
  input.pu_age = 20.0 * 365.25 * 24.0 * 60.0 * 60.0;

  const auto output = RelActCalc::correct_pu_mass_fractions_for_pu242(
      input, RelActCalc::PuCorrMethod::NotApplicable );

  BOOST_CHECK_SMALL( output.pu239_mass_frac - 0.9, 1.0E-12 );
  BOOST_CHECK_SMALL( output.pu240_mass_frac, 1.0E-15 );
  BOOST_CHECK_SMALL( output.pu241_mass_frac - 0.1, 1.0E-12 );
  BOOST_CHECK_SMALL( output.pu242_mass_frac, 1.0E-15 );
}


BOOST_AUTO_TEST_CASE( FrozenContinuumRankActiveSetScalarJetAndOrder )
{
  using Jet = ceres::Jet<double,1>;
  constexpr size_t nbin = 9;
  constexpr float energies[nbin+1]
      = {100.0f,101.0f,102.0f,103.0f,104.0f,105.0f,106.0f,107.0f,108.0f,109.0f};
  constexpr float data[nbin]
      = {65.0f,71.0f,78.0f,93.0f,118.0f,91.0f,79.0f,70.0f,66.0f};
  constexpr double ref_energy = energies[0];
  const PeakContinuum::OffsetType type = PeakContinuum::OffsetType::Quadratic;
  const size_t ncoef = PeakContinuum::num_parameters(type);

  const auto scalar_fit = [&]( const vector<PeakDef> &peaks,
                               const PeakFit::ContinuumFitPolicy *policy,
                               PeakFit::ContinuumFitPolicy *observed ) {
    vector<double> coefficients(ncoef,0.0), counts(nbin,0.0);
    PeakFit::fit_continuum( energies,data,nullptr,nbin,type,ref_energy,peaks,false,
                            coefficients.data(),counts.data(),policy,observed );
    return std::make_pair(coefficients,counts);
  };

  vector<PeakDef> peaks;
  peaks.emplace_back(103.7,0.82,54.0);
  peaks.emplace_back(105.4,1.05,31.0);

  PeakFit::ContinuumFitPolicy classified;
  const auto unfrozen = scalar_fit(peaks,nullptr,&classified);
  BOOST_REQUIRE( classified.initialized );
  BOOST_CHECK_EQUAL( classified.continuum_type,type );
  BOOST_CHECK_EQUAL( classified.num_channels,nbin );
  BOOST_CHECK_EQUAL( classified.num_polynomial_terms,ncoef );
  BOOST_CHECK_EQUAL( classified.num_lls_terms,ncoef );
  BOOST_CHECK_EQUAL( classified.numerical_rank,ncoef );
  BOOST_REQUIRE_EQUAL( classified.clamped_channels.size(),nbin );

  // Reversing semantic peak order must not change the classified/frozen mathematical problem.
  vector<PeakDef> reversed(peaks.rbegin(),peaks.rend());
  PeakFit::ContinuumFitPolicy reverse_classified;
  const auto reverse_unfrozen = scalar_fit(reversed,nullptr,&reverse_classified);
  BOOST_CHECK( PeakFit::same_continuum_fit_policy(classified,reverse_classified) );
  for( size_t bin = 0; bin < nbin; ++bin )
    BOOST_CHECK_SMALL( unfrozen.second[bin] - reverse_unfrozen.second[bin],1.0e-11 );

  // Force a frozen clamp decision even if the raw continuum later changes sign.  The selected
  // channel is then exactly peak-only, while every other channel follows the fixed-rank LLS.
  PeakFit::ContinuumFitPolicy frozen = classified;
  frozen.clamped_channels[0] = 1;
  const auto frozen_scalar = scalar_fit(peaks,&frozen,nullptr);
  double peak_only[nbin] = {0.0};
  for( const PeakDef &peak : peaks )
    peak.gauss_integral(energies,peak_only,nbin);
  BOOST_CHECK_SMALL( frozen_scalar.second[0] - peak_only[0],1.0e-12 );

  // Rank is a genuine frozen decision, not merely diagnostic metadata.
  PeakFit::ContinuumFitPolicy zero_rank = classified;
  zero_rank.numerical_rank = 0;
  std::fill(zero_rank.clamped_channels.begin(),zero_rank.clamped_channels.end(),0);
  const auto peak_only_fit = scalar_fit(peaks,&zero_rank,nullptr);
  for( size_t coefficient = 0; coefficient < ncoef; ++coefficient )
    BOOST_CHECK_SMALL( peak_only_fit.first[coefficient],1.0e-14 );
  for( size_t bin = 0; bin < nbin; ++bin )
    BOOST_CHECK_SMALL( peak_only_fit.second[bin] - peak_only[bin],1.0e-12 );

  // The Jet path must apply the same pseudo-inverse and active mask.  Check both values and the
  // derivative with respect to the first peak's amplitude against a fresh scalar finite difference.
  vector<RelActCalcAuto::PeakDefImp<Jet>> jet_peaks(2);
  for( size_t index = 0; index < peaks.size(); ++index )
  {
    jet_peaks[index].m_mean = Jet(peaks[index].mean());
    jet_peaks[index].m_sigma = Jet(peaks[index].sigma());
    jet_peaks[index].m_amplitude = Jet(peaks[index].amplitude());
  }
  jet_peaks[0].m_amplitude.v[0] = 1.0;
  vector<Jet> jet_coefficients(ncoef),jet_counts(nbin);
  PeakFit::fit_continuum( energies,data,nullptr,nbin,type,Jet(ref_energy),jet_peaks,false,
                          jet_coefficients.data(),jet_counts.data(),&frozen,nullptr );
  for( size_t bin = 0; bin < nbin; ++bin )
    BOOST_CHECK_SMALL( jet_counts[bin].a - frozen_scalar.second[bin],1.0e-10 );

  const double epsilon = 1.0e-4;
  vector<PeakDef> plus = peaks,minus = peaks;
  plus[0].setAmplitude(peaks[0].amplitude() + epsilon);
  minus[0].setAmplitude(peaks[0].amplitude() - epsilon);
  const auto plus_fit = scalar_fit(plus,&frozen,nullptr);
  const auto minus_fit = scalar_fit(minus,&frozen,nullptr);
  for( size_t bin = 0; bin < nbin; ++bin )
  {
    const double finite_difference
        = (plus_fit.second[bin] - minus_fit.second[bin])/(2.0*epsilon);
    BOOST_CHECK_SMALL( jet_counts[bin].v[0] - finite_difference,
                       2.0e-6*(1.0 + fabs(finite_difference)) );
  }

  PeakFit::ContinuumFitPolicy wrong_basis = frozen;
  wrong_basis.num_channels += 1;
  vector<double> coefficients(ncoef),counts(nbin);
  BOOST_CHECK_THROW(
    PeakFit::fit_continuum(energies,data,nullptr,nbin,type,ref_energy,peaks,false,
                           coefficients.data(),counts.data(),&wrong_basis,nullptr),
    std::logic_error );
}


BOOST_AUTO_TEST_CASE( combine_nearby_peaks_preserves_moments_and_uncertainty )
{
  RelActCalcManual::GenericPeakInfo first;
  first.m_energy = 100.0;
  first.m_mean = 100.1;
  first.m_fwhm = 2.0;
  first.m_counts = 100.0;
  first.m_counts_uncert = 10.0;

  RelActCalcManual::GenericPeakInfo second;
  second.m_energy = 100.5;
  second.m_mean = 100.7;
  second.m_fwhm = 4.0;
  second.m_counts = 300.0;
  second.m_counts_uncert = 20.0;

  RelActCalcManual::GenericPeakInfo far_peak;
  far_peak.m_energy = 110.0;
  far_peak.m_mean = 110.0;
  far_peak.m_fwhm = 2.0;
  far_peak.m_counts = 50.0;
  far_peak.m_counts_uncert = 7.0;

  vector<RelActCalcManual::GenericPeakInfo> peaks{ far_peak, second, first };
  const size_t combined = RelActCalcAutoImp::combine_nearby_peaks( peaks, 1.85 );

  BOOST_REQUIRE_EQUAL( combined, 1 );
  BOOST_REQUIRE_EQUAL( peaks.size(), 2 );

  const double expected_energy = (100.0*100.0 + 300.0*100.5) / 400.0;
  const double expected_mean = (100.0*100.1 + 300.0*100.7) / 400.0;
  const double sigma_1 = 2.0 / 2.35482;
  const double sigma_2 = 4.0 / 2.35482;
  const double expected_sigma = sqrt(
      (100.0*(sigma_1*sigma_1 + pow(100.0 - expected_energy, 2.0))
       + 300.0*(sigma_2*sigma_2 + pow(100.5 - expected_energy, 2.0))) / 400.0 );

  BOOST_CHECK_SMALL( peaks[0].m_counts - 400.0, 1.0E-12 );
  BOOST_CHECK_SMALL( peaks[0].m_counts_uncert - sqrt(500.0), 1.0E-12 );
  BOOST_CHECK_SMALL( peaks[0].m_energy - expected_energy, 1.0E-12 );
  BOOST_CHECK_SMALL( peaks[0].m_mean - expected_mean, 1.0E-12 );
  BOOST_CHECK_SMALL( peaks[0].m_fwhm - 2.35482*expected_sigma, 1.0E-12 );
  BOOST_CHECK_EQUAL( peaks[1].m_energy, far_peak.m_energy );
}


// ============================================================================
// Tests for deviation pair energy calibration corrections
// ============================================================================


BOOST_AUTO_TEST_CASE( test_deviation_pair_xml_serialization )
{
  set_data_dir();

  // Test that energy_cal_type is properly serialized and deserialized
  using namespace RelActCalcAuto;

  Options opts_orig;
  opts_orig.energy_cal_type = EnergyCalFitType::NonLinearFit;
  opts_orig.fwhm_form = FwhmForm::Polynomial_2;

  // Add a dummy ROI and nuclide to make XML valid
  RoiRange roi;
  roi.lower_energy = 100.0;
  roi.upper_energy = 200.0;
  opts_orig.rois.push_back( roi );

  RelEffCurveInput curve;
  curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::LnX;
  curve.rel_eff_eqn_order = 3;

  NucInputInfo nuc_info;
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  if( db )
  {
    nuc_info.source = db->nuclide("Co60");
    curve.nuclides.push_back( nuc_info );
    opts_orig.rel_eff_curves.push_back( curve );
  }

  // Serialize to XML
  rapidxml::xml_document<char> doc;
  rapidxml::xml_node<char> *root = doc.allocate_node( rapidxml::node_element, "Root" );
  doc.append_node( root );

  opts_orig.toXml( root );

  // Deserialize from XML
  Options opts_loaded;
  const rapidxml::xml_node<char> *opts_node = root->first_node( "Options" );
  BOOST_REQUIRE( opts_node );

  opts_loaded.fromXml( opts_node );

  BOOST_CHECK( opts_loaded.energy_cal_type == opts_orig.energy_cal_type );
}//BOOST_AUTO_TEST_CASE( test_deviation_pair_xml_serialization )



BOOST_AUTO_TEST_CASE( test_options_equal_enough_with_deviation_pairs )
{
#if( PERFORM_DEVELOPER_CHECKS )
  // Test that Options::equalEnough properly checks fit_deviation_pairs
  using namespace RelActCalcAuto;

  Options opts1, opts2;

  // Both NoFit - should be equal
  opts1.energy_cal_type = RelActCalcAuto::EnergyCalFitType::NoFit;
  opts2.energy_cal_type = RelActCalcAuto::EnergyCalFitType::NoFit;
  BOOST_CHECK_NO_THROW( Options::equalEnough( opts1, opts2 ) );

  // Both LinearFit - should be equal
  opts1.energy_cal_type = RelActCalcAuto::EnergyCalFitType::LinearFit;
  opts2.energy_cal_type = RelActCalcAuto::EnergyCalFitType::LinearFit;
  BOOST_CHECK_NO_THROW( Options::equalEnough( opts1, opts2 ) );

  // Both NonLinearFit - should be equal
  opts1.energy_cal_type = RelActCalcAuto::EnergyCalFitType::NonLinearFit;
  opts2.energy_cal_type = RelActCalcAuto::EnergyCalFitType::NonLinearFit;
  BOOST_CHECK_NO_THROW( Options::equalEnough( opts1, opts2 ) );

  // Different - should throw
  opts1.energy_cal_type = RelActCalcAuto::EnergyCalFitType::LinearFit;
  opts2.energy_cal_type = RelActCalcAuto::EnergyCalFitType::NoFit;
  BOOST_CHECK_THROW( Options::equalEnough( opts1, opts2 ), std::runtime_error );
#endif
}//BOOST_AUTO_TEST_CASE( test_options_equal_enough_with_deviation_pairs )


BOOST_AUTO_TEST_CASE( test_rel_act_auto_gui_state_equality )
{
  // `RelActAutoGuiState::operator==` is what `RelActAutoGui` uses to decide whether a GUI change
  //  warrants an undo/redo step, so every field that the user can change must take part in it -
  //  a field left out means that change silently gets no undo step.
  using namespace RelActCalcAuto;

  const RelActAutoGuiState default_state;
  BOOST_CHECK( default_state == RelActAutoGuiState{} );
  BOOST_CHECK( !(default_state != RelActAutoGuiState{}) );

  {
    RelActAutoGuiState state;
    state.note = "a user note";
    BOOST_CHECK( state != default_state );
  }

  {
    RelActAutoGuiState state;
    state.description = "a description";
    BOOST_CHECK( state != default_state );
  }

  {
    RelActAutoGuiState state;
    state.background_subtract = !default_state.background_subtract;
    BOOST_CHECK( state != default_state );
  }

  {
    RelActAutoGuiState state;
    state.show_ref_lines = !default_state.show_ref_lines;
    BOOST_CHECK( state != default_state );
  }

  {
    // The chart zoom is part of the state (it is saved/restored with the analysis), so a zoom
    //  must compare unequal - `RelActAutoGui` keeps its undo baseline in sync with the chart so
    //  that a zoom does not get folded into the next edits undo step.
    RelActAutoGuiState state;
    state.lower_display_energy = 100.0;
    state.upper_display_energy = 2000.0;
    BOOST_CHECK( state != default_state );
  }

  {
    RelActAutoGuiState state;
    state.options.energy_cal_type = EnergyCalFitType::LinearFit;
    BOOST_CHECK( state != default_state );
  }

  {
    RelActAutoGuiState state;
    state.options.skew_type = PeakDef::SkewType::CrystalBall;
    BOOST_CHECK( state != default_state );
  }

  {
    // ROI *order* is part of the compared state - `getRoiRanges()` walks the energy-range widgets
    //  in display order, so re-sorting them in the GUI is a real change that needs an undo step.
    RoiRange lower, upper;
    lower.lower_energy = 100.0;  lower.upper_energy = 200.0;
    upper.lower_energy = 300.0;  upper.upper_energy = 400.0;

    RelActAutoGuiState ascending, descending;
    ascending.options.rois = { lower, upper };
    descending.options.rois = { upper, lower };

    BOOST_CHECK( ascending != descending );
    BOOST_CHECK( ascending != default_state );
  }
}//BOOST_AUTO_TEST_CASE( test_rel_act_auto_gui_state_equality )


BOOST_AUTO_TEST_CASE( test_rel_act_auto_gui_state_xml_round_trip )
{
  using namespace RelActCalcAuto;

  RelActAutoGuiState state;
  state.note = "round trip note";
  state.description = "round trip description";
  state.background_subtract = true;
  state.show_ref_lines = true;
  state.lower_display_energy = 55.5;
  state.upper_display_energy = 2750.25;
  state.options.energy_cal_type = EnergyCalFitType::LinearFit;
  state.options.skew_type = PeakDef::SkewType::CrystalBall;

  rapidxml::xml_document<char> doc;
  BOOST_REQUIRE_NO_THROW( state.serialize( &doc ) );

  const rapidxml::xml_node<char> * const base_node = doc.first_node( "RelActCalcAuto" );
  BOOST_REQUIRE( base_node );

  RelActAutoGuiState loaded;
  BOOST_REQUIRE_NO_THROW( loaded.deSerialize( base_node ) );

  BOOST_CHECK( loaded == state );
}//BOOST_AUTO_TEST_CASE( test_rel_act_auto_gui_state_xml_round_trip )
