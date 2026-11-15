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

/* Tests of which gamma lines outside a RelActCalcAuto ROI that ROI models (RelActCalcAuto_LineWindow_imp.hpp):
   the weight functions, and that the objective stays continuous as the set of modelled lines changes. */

#include "InterSpec_config.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <limits>
#include <algorithm>

#include "ceres/jet.h"

#define BOOST_TEST_MODULE test_RelActCalcAuto_LineWindow
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/RelActCalcAuto_LineWindow_imp.hpp"

using namespace std;
using namespace boost::unit_test;
namespace LineWindow = RelActCalcAutoImp::LineWindow;

static_assert( USE_REL_ACT_TOOL, "Compile-time support for Rel Act tool is required." );

typedef ceres::Jet<double,1> Jet1;


BOOST_AUTO_TEST_CASE( smoothstep_is_c2_and_saturates )
{
  BOOST_CHECK_EQUAL( LineWindow::smoothstep( -1.0 ), 0.0 );
  BOOST_CHECK_EQUAL( LineWindow::smoothstep( 0.0 ), 0.0 );
  BOOST_CHECK_EQUAL( LineWindow::smoothstep( 1.0 ), 1.0 );
  BOOST_CHECK_EQUAL( LineWindow::smoothstep( 2.0 ), 1.0 );
  BOOST_CHECK_CLOSE( LineWindow::smoothstep( 0.5 ), 0.5, 1.0E-12 );

  double previous = 0.0;
  for( int i = 1; i <= 1000; ++i )
  {
    const double u = i / 1000.0;
    const double value = LineWindow::smoothstep( u );
    BOOST_CHECK( value >= previous );
    BOOST_CHECK_SMALL( value + LineWindow::smoothstep( 1.0 - u ) - 1.0, 1.0E-12 );
    previous = value;
  }

  // The first and second derivatives go to zero at both ends, so the weight joins its constant parts with
  //  a continuous gradient and curvature.
  for( const double u : { 1.0E-4, 1.0 - 1.0E-4 } )
  {
    const Jet1 s = LineWindow::smoothstep( Jet1( u, 0 ) );
    BOOST_CHECK_SMALL( s.v[0], 1.0E-6 );
  }
  const Jet1 mid = LineWindow::smoothstep( Jet1( 0.5, 0 ) );
  BOOST_CHECK_CLOSE( mid.v[0], 1.875, 1.0E-9 );
}//smoothstep_is_c2_and_saturates


BOOST_AUTO_TEST_CASE( strength_weight_follows_counts_over_noise )
{
  const double lo = LineWindow::sm_min_counts_in_noise, hi = LineWindow::sm_full_counts_in_noise;

  BOOST_CHECK_EQUAL( LineWindow::strength_weight( 0.0 ), 0.0 );
  BOOST_CHECK_EQUAL( LineWindow::strength_weight( lo ), 0.0 );
  BOOST_CHECK_EQUAL( LineWindow::strength_weight( hi ), 1.0 );
  BOOST_CHECK_EQUAL( LineWindow::strength_weight( 1.0E9 ), 1.0 );
  BOOST_CHECK_EQUAL( LineWindow::strength_weight( std::numeric_limits<double>::quiet_NaN() ), 0.0 );
  BOOST_CHECK_CLOSE( LineWindow::strength_weight( std::sqrt(lo*hi) ), 0.5, 1.0E-9 );

  // Monotone, continuous across the band, and the Jet derivative matches a central difference.
  double previous = 0.0;
  for( int i = 0; i <= 400; ++i )
  {
    const double g = lo * std::pow( hi/lo, -0.25 + 1.5*i/400.0 );
    const double w = LineWindow::strength_weight( g );
    BOOST_CHECK( (w >= previous) && (w >= 0.0) && (w <= 1.0) );
    previous = w;

    const double h = 1.0E-6 * g;
    const double numeric = (LineWindow::strength_weight( g + h ) - LineWindow::strength_weight( g - h )) / (2.0*h);
    const Jet1 jet = LineWindow::strength_weight( Jet1( g, 0 ) );
    BOOST_CHECK_SMALL( jet.a - w, 1.0E-12 );
    BOOST_CHECK_SMALL( jet.v[0] - numeric, 1.0E-5 * (1.0 + std::fabs(numeric)) );
  }
}//strength_weight_follows_counts_over_noise


BOOST_AUTO_TEST_CASE( fraction_weight_follows_own_area_in_roi )
{
  const double lo = LineWindow::sm_min_fraction_in_roi, hi = LineWindow::sm_full_fraction_in_roi;

  BOOST_CHECK_EQUAL( LineWindow::fraction_weight( 0.0 ), 0.0 );
  BOOST_CHECK_EQUAL( LineWindow::fraction_weight( lo ), 0.0 );
  BOOST_CHECK_EQUAL( LineWindow::fraction_weight( hi ), 1.0 );
  BOOST_CHECK_EQUAL( LineWindow::fraction_weight( 0.5 ), 1.0 );
  BOOST_CHECK_CLOSE( LineWindow::fraction_weight( std::sqrt(lo*hi) ), 0.5, 1.0E-9 );

  // The Crystal Ball scale moves the band up by ten.
  BOOST_CHECK_EQUAL( LineWindow::fraction_weight( hi, 10.0 ), 0.0 );
  BOOST_CHECK_EQUAL( LineWindow::fraction_weight( 10.0*hi, 10.0 ), 1.0 );
  BOOST_CHECK_CLOSE( LineWindow::fraction_weight( 10.0*std::sqrt(lo*hi), 10.0 ), 0.5, 1.0E-9 );

  // The general band is what both weights use.
  for( const double v : { 0.02, 0.03, 0.05, 0.07 } )
  {
    BOOST_CHECK_EQUAL( LineWindow::log_band_weight( v, LineWindow::sm_min_counts_in_noise,
                                                    LineWindow::sm_full_counts_in_noise ),
                       LineWindow::strength_weight( v ) );
  }
}//fraction_weight_follows_own_area_in_roi


BOOST_AUTO_TEST_CASE( cap_weight_tapers_beyond_the_cap )
{
  const double cap = 15.0, end = cap * (1.0 + LineWindow::sm_cap_taper_fraction);

  BOOST_CHECK_EQUAL( LineWindow::cap_weight( 0.0, cap ), 1.0 );
  BOOST_CHECK_EQUAL( LineWindow::cap_weight( cap, cap ), 1.0 );
  BOOST_CHECK_EQUAL( LineWindow::cap_weight( end, cap ), 0.0 );
  BOOST_CHECK_EQUAL( LineWindow::cap_weight( 100.0, cap ), 0.0 );
  BOOST_CHECK_EQUAL( LineWindow::cap_weight( 1.0E6, 0.0 ), 1.0 );   // no cap
  BOOST_CHECK_CLOSE( LineWindow::cap_weight( 0.5*(cap + end), cap ), 0.5, 1.0E-9 );

  double previous = 1.0;
  for( int i = 0; i <= 200; ++i )
  {
    const double n = cap - 1.0 + (end - cap + 2.0)*i/200.0;
    const double w = LineWindow::cap_weight( n, cap );
    BOOST_CHECK( (w <= previous) && (w >= 0.0) && (w <= 1.0) );
    previous = w;
  }

  for( const double n : { cap + 1.0E-4, end - 1.0E-4 } )
  {
    const Jet1 w = LineWindow::cap_weight( Jet1( n, 0 ), cap );
    BOOST_CHECK_SMALL( w.v[0], 1.0E-6 );
  }
}//cap_weight_tapers_beyond_the_cap


BOOST_AUTO_TEST_CASE( reach_beyond_edge_bounds_where_the_cap_weight_is_nonzero )
{
  const double cap = 15.0, k = cap * (1.0 + LineWindow::sm_cap_taper_fraction);

  // A detector-like resolution: sigma grows with energy.
  const auto sigma = []( const double energy ) -> double { return 0.2 + 0.0004*energy; };

  for( const double edge : { 60.0, 200.0, 1000.0 } )
  {
    const double above = LineWindow::reach_beyond_edge( sigma, edge, +1, cap, 1.0E4 );
    const double below = LineWindow::reach_beyond_edge( sigma, edge, -1, cap, edge );
    BOOST_CHECK( (above > 0.0) && (below > 0.0) );

    // Every line further out has a zero cap weight.
    for( int i = 0; i <= 100; ++i )
    {
      const double d_above = above * (1.0 + 0.01 + 3.0*i/100.0);
      BOOST_CHECK_EQUAL( LineWindow::cap_weight( d_above / sigma(edge + d_above), cap ), 0.0 );

      const double d_below = std::min( edge - 1.0E-3, below * (1.0 + 0.01 + 0.5*i/100.0) );
      if( d_below > below )
        BOOST_CHECK_EQUAL( LineWindow::cap_weight( d_below / sigma(edge - d_below), cap ), 0.0 );
    }

    // ... and the reach is not much more than needed.
    BOOST_CHECK_LE( above, 1.01 * k * sigma(edge + above) );
    BOOST_CHECK_LE( below, 1.01 * k * sigma(edge) );
  }//for( loop over edges )

  // The largest reach considered, and an unusable width, give the whole range.
  BOOST_CHECK_EQUAL( LineWindow::reach_beyond_edge( sigma, 200.0, +1, cap, 1.0 ), 1.0 );
  const auto bad_sigma = []( const double ) -> double { return std::numeric_limits<double>::quiet_NaN(); };
  BOOST_CHECK_EQUAL( LineWindow::reach_beyond_edge( bad_sigma, 200.0, +1, cap, 500.0 ), 500.0 );
  BOOST_CHECK_EQUAL( LineWindow::reach_beyond_edge( sigma, 200.0, +1, 0.0, 500.0 ), 500.0 );
  BOOST_CHECK_EQUAL( LineWindow::reach_beyond_edge( sigma, 200.0, -1, cap, -3.0 ), 0.0 );
}//reach_beyond_edge_bounds_where_the_cap_weight_is_nonzero


BOOST_AUTO_TEST_CASE( edge_noise_is_smallest_edge_uncertainty )
{
  const vector<float> uncerts{ 9.0f, 4.0f, 7.0f, 30.0f, 30.0f, 12.0f, 5.0f, 8.0f };
  BOOST_CHECK_EQUAL( LineWindow::edge_noise( uncerts.data(), uncerts.size(), true ), 4.0 );
  BOOST_CHECK_EQUAL( LineWindow::edge_noise( uncerts.data(), uncerts.size(), false ), 5.0 );

  const vector<float> tiny{ 0.0f, 0.5f };
  BOOST_CHECK_EQUAL( LineWindow::edge_noise( tiny.data(), tiny.size(), true ), 1.0 );
  BOOST_CHECK_EQUAL( LineWindow::edge_noise( tiny.data(), 0, true ), 1.0 );
}//edge_noise_is_smallest_edge_uncertainty


#if( PERFORM_DEVELOPER_CHECKS )
namespace
{
  string g_test_file_dir;

  void set_data_dir()
  {
    static bool s_have_set = false;
    if( s_have_set )
      return;
    s_have_set = true;

    const int argc = framework::master_test_suite().argc;
    char ** const argv = framework::master_test_suite().argv;

    string datadir;
    for( int i = 1; i < argc; ++i )
    {
      const string arg = argv[i];
      if( SpecUtils::istarts_with( arg, "--datadir=" ) )
        datadir = arg.substr( 10 );
      if( SpecUtils::istarts_with( arg, "--testfiledir=" ) )
        g_test_file_dir = arg.substr( 14 );
    }

    SpecUtils::ireplace_all( datadir, "%20", " " );
    SpecUtils::ireplace_all( g_test_file_dir, "%20", " " );

    if( datadir.empty() )
    {
      for( const char *d : { "data", "../data", "../../data", "../../../data" } )
      {
        if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) )
        {
          datadir = d;
          break;
        }
      }
    }
    if( g_test_file_dir.empty() )
    {
      for( const char *d : { "test_data", "../test_data", "../../test_data", "../../target/testing/test_data" } )
      {
        if( SpecUtils::is_directory( SpecUtils::append_path(d, "RelActAutoReport") ) )
        {
          g_test_file_dir = d;
          break;
        }
      }
    }

    const string sandia_decay_file = SpecUtils::append_path( datadir, "sandia.decay.xml" );
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( sandia_decay_file ),
                           "sandia.decay.xml not at '" << sandia_decay_file << "'" );

    BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );

    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    BOOST_REQUIRE_MESSAGE( db && db->nuclide("U238"), "Error initing SandiaDecayDataBase" );
  }//void set_data_dir()
}//namespace


/** Along each FWHM and skew coefficient, from the converged answer of the uranium HPGe fixture, the set of
 lines the ROIs model changes at many points (lines outside a ROI enter and leave it as their peaks widen or
 grow tails).  At every such change the objective must not jump.

 Before the lines outside a ROI were weighted smoothly, membership was a hard window from the peak-coverage
 limits at the ROI edges, and the objective jumped by the line's contribution at each change: here up to
 2.6E-5 (at chi2 138) over 70 changes, and up to tens of chi2 on the benchmark spectra. */
BOOST_AUTO_TEST_CASE( objective_is_continuous_as_line_membership_changes )
{
  set_data_dir();

  const string spec_path = SpecUtils::append_path(
      SpecUtils::append_path( g_test_file_dir, "RelActAutoReport" ), "U235_Unshielded_6000.n42" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( spec_path ), "Test spectrum not at '" << spec_path << "'" );

  auto meas = make_shared<SpecMeas>();
  BOOST_REQUIRE( meas->load_file( spec_path, SpecUtils::ParserType::Auto ) );

  shared_ptr<const SpecUtils::Measurement> foreground, background;
  for( const shared_ptr<const SpecUtils::Measurement> &m : meas->measurements() )
  {
    if( !m )
      continue;
    if( m->source_type() == SpecUtils::SourceType::Background )
      background = m;
    else if( !foreground )
      foreground = m;
  }
  BOOST_REQUIRE( foreground );
  BOOST_REQUIRE( meas->detector() );

  unique_ptr<RelActCalcAuto::RelActAutoGuiState> state = meas->getRelActAutoGuiState();
  BOOST_REQUIRE( state );
  RelActCalcAuto::Options options = state->options;
  options.energy_cal_type = RelActCalcAuto::EnergyCalFitType::LinearFit;
  options.auto_simplify_model = false;

  const PeakFitUtils::CoarseResolutionType det_type
      = PeakFitUtils::effective_det_type( meas->peakFitDetPrefs(), foreground, meas );

  RelActCalcAuto::RelActAutoSolution sol;
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, foreground, background, meas->detector(), {},
                                                       det_type, nullptr ) );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status( sol.m_status ),
                         "Solve failed: '" << sol.m_error_message << "'" );
  BOOST_REQUIRE_EQUAL( sol.m_parameter_names.size(), sol.m_final_parameters.size() );

  const vector<double> &x0 = sol.m_final_parameters;
  double chi2_0 = 0.0;
  BOOST_REQUIRE_NO_THROW( chi2_0 = sol.objective_chi2( x0 ) );
  BOOST_REQUIRE( std::isfinite(chi2_0) );

  // The FWHM and skew coefficients that are free (the Berstein FWHM forms carry their energy range as two
  //  fixed parameters, and an unused upper skew value is a -999.9 placeholder).
  vector<size_t> directions;
  for( size_t i = 0; i < x0.size(); ++i )
  {
    const string &name = sol.m_parameter_names[i];
    const bool fit = (i < sol.m_parameter_were_fit.size()) ? sol.m_parameter_were_fit[i] : true;
    if( fit && (SpecUtils::starts_with( name, "FWHM_" ) || SpecUtils::starts_with( name, "Skew_" ))
       && (x0[i] > -999.0) )
      directions.push_back( i );
  }
  BOOST_REQUIRE( !directions.empty() );

  size_t num_changes = 0;
  double max_jump = 0.0;
  string max_jump_desc;
  for( const size_t index : directions )
  {
    // A point the model cannot be evaluated at (e.g., a width too narrow to resolve) counts as NaN.
    const auto eval = [&]( const double value, size_t &num_peaks ) -> double {
      vector<double> x = x0;
      x[index] = value;
      num_peaks = 0;
      try
      {
        return sol.objective_chi2( x, &num_peaks );
      }catch( std::exception & )
      {
        return std::numeric_limits<double>::quiet_NaN();
      }
    };

    const double center = x0[index];
    const double half_range = 0.25 * std::max( std::fabs(center), 1.0E-3 );
    const size_t npoints = 81;

    double prev_value = center - half_range;
    size_t prev_count = 0;
    double prev_chi2 = eval( prev_value, prev_count );

    for( size_t i = 1; i < npoints; ++i )
    {
      const double value = center - half_range + 2.0*half_range*i/(npoints - 1.0);
      size_t count = 0;
      const double chi2 = eval( value, count );

      if( (count != prev_count) && std::isfinite(chi2) && std::isfinite(prev_chi2) )
      {
        // Narrow down to where the modelled lines change, then compare the objective just either side.
        double a = prev_value, b = value;
        const size_t count_a = prev_count;
        for( int iter = 0; iter < 60; ++iter )
        {
          const double m = 0.5*(a + b);
          if( (m <= a) || (m >= b) )
            break;
          size_t count_m = 0;
          if( !std::isfinite( eval( m, count_m ) ) )
            break;
          if( count_m == count_a )
            a = m;
          else
            b = m;
        }//for( bisection )

        size_t count_lo = 0, count_hi = 0;
        const double chi2_lo = eval( a, count_lo ), chi2_hi = eval( b, count_hi );
        if( !std::isfinite(chi2_lo) || !std::isfinite(chi2_hi) || (count_lo == count_hi) )
        {
          prev_value = value;
          prev_count = count;
          prev_chi2 = chi2;
          continue;
        }
        const double jump = std::fabs( chi2_hi - chi2_lo );
        num_changes += 1;
        if( jump > max_jump )
        {
          max_jump = jump;
          max_jump_desc = sol.m_parameter_names[index] + " at " + std::to_string( a )
                          + " (" + std::to_string(count_lo) + " -> " + std::to_string(count_hi) + " peaks, chi2 "
                          + std::to_string(chi2_lo) + ")";
        }

        // Round-off is ~1E-15 of chi2; before the smooth weights the jumps here were ~1E-7 to 1E-5 of chi2.
        BOOST_CHECK_MESSAGE( jump <= 1.0E-9 * std::max( 1.0, std::fabs(chi2_lo) ),
                             "Objective jumps by " << jump << " where the modelled lines change from "
                             << count_lo << " to " << count_hi << " peaks, along " << sol.m_parameter_names[index]
                             << " at " << a << " (chi2 " << chi2_lo << ")" );
      }//if( the modelled lines changed )

      prev_value = value;
      prev_count = count;
      prev_chi2 = chi2;
    }//for( loop over scan points )
  }//for( loop over directions )

  BOOST_TEST_MESSAGE( "Line-membership changes found: " << num_changes << "; largest objective change across one: "
                      << max_jump << (max_jump_desc.empty() ? string() : (" - " + max_jump_desc)) );
  BOOST_CHECK_MESSAGE( num_changes > 0, "Expected the scan to cross changes in which lines are modelled." );
}//objective_is_continuous_as_line_membership_changes
#endif //PERFORM_DEVELOPER_CHECKS
