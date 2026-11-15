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
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "InterSpec_config.h"

#include <string>
#include <vector>
#include <cstddef>
#include <utility>
#include <algorithm>

#define BOOST_TEST_MODULE RelActCalcAuto_PreFit_suite
#include <boost/test/included/unit_test.hpp>

#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/RelActCalcManual.h"
#include "InterSpec/RelActCalcAuto_PreFit.hpp"

using namespace std;
using RelActCalcAutoImp::PreFit::StageGroups;
using RelActCalcAutoImp::PreFit::peaks_mostly_sums;

// How RelActCalcAuto starts a fit that models sum peaks (see RelActCalcAuto_PreFit.hpp).  The fits these
//  protect - close or high-rate uranium spectra - are covered end to end by the isotopics benchmark, and by
//  `solver_close_source_with_sum_peaks` in test_RelActCalcCascade.cpp.

namespace
{
/** A peak of `fwhm` keV at `energy`, of `counts`, matched to the given lines {nuclide, yield}. */
RelActCalcManual::GenericPeakInfo peak( const double energy, const double counts,
                                        const vector<pair<string,double>> &lines, const double fwhm = 1.5 )
{
  RelActCalcManual::GenericPeakInfo answer;
  answer.m_energy = answer.m_mean = energy;
  answer.m_fwhm = fwhm;
  answer.m_counts = counts;
  answer.m_counts_uncert = std::sqrt( counts );
  for( const pair<string,double> &line : lines )
    answer.m_source_gammas.emplace_back( line.second, line.first );
  return answer;
}//peak(...)

const double sm_cluster_num_sigma = 1.5;  //as the solver matches lines to peaks
const double sm_max_share = 0.1;          //RelActAutoCostFcn::sm_initial_estimate_max_sum_share
}//namespace


// Close to the detector, U-235's 345.9 and 387.8 keV peaks are mostly coincidence and pile-up sums; they say
//  little about the lines, so the initial estimate does not use them.  Strong lines with small sums on them do.
BOOST_AUTO_TEST_CASE( peaks_that_are_mostly_sums_are_left_out )
{
  const vector<RelActCalcManual::GenericPeakInfo> peaks{
    peak( 185.715, 3.8E7, {{"U235", 0.572}} ),
    peak( 345.900, 3.1E4, {{"U235", 3.0E-4}} ),
    peak( 387.820, 1.3E4, {{"U235", 3.2E-5}, {"U238", 9.5E-6}} ),
    peak( 1000.99, 1.5E6, {{"U238", 8.4E-3}} )
  };

  const vector<pair<double,double>> sums{
    { 345.9, 2.5E4 },     //most of the peak
    { 387.8, 6.0E3 },     //two sums inside the window: together most of the peak
    { 387.4, 3.0E3 },
    { 1001.03, 1.0E3 },   //a tiny fraction of the peak
    { 186.4, 1.0E5 }      //a quarter of a percent of the 185.7 keV peak
  };

  const vector<bool> leave_out = peaks_mostly_sums( peaks, sums, sm_cluster_num_sigma, sm_max_share );
  BOOST_REQUIRE_EQUAL( leave_out.size(), peaks.size() );
  BOOST_CHECK( !leave_out[0] );
  BOOST_CHECK( leave_out[1] );
  BOOST_CHECK( leave_out[2] );
  BOOST_CHECK( !leave_out[3] );
}//peaks_that_are_mostly_sums_are_left_out


// Only sums within the line-matching window of a peak count against it, and only above the set share.
BOOST_AUTO_TEST_CASE( window_and_share )
{
  const double fwhm = 1.5, sigma = fwhm / PhysicalUnits::fwhm_nsigma;
  const vector<RelActCalcManual::GenericPeakInfo> peaks{
    peak( 185.715, 1.0E6, {{"U235", 0.572}}, fwhm ),
    peak( 300.0, 1.0E4, {{"U235", 1.0E-3}}, fwhm )
  };

  // Just outside the window of the 300 keV peak, however large: not counted.
  vector<pair<double,double>> sums{ { 300.0 + 1.01*sm_cluster_num_sigma*sigma, 1.0E6 } };
  BOOST_CHECK( !peaks_mostly_sums( peaks, sums, sm_cluster_num_sigma, sm_max_share )[1] );

  // Just inside it.
  sums = { { 300.0 + 0.99*sm_cluster_num_sigma*sigma, 1.0E6 } };
  BOOST_CHECK( peaks_mostly_sums( peaks, sums, sm_cluster_num_sigma, sm_max_share )[1] );

  // At the peak, either side of the share.
  sums = { { 300.0, 0.99*sm_max_share*1.0E4 } };
  BOOST_CHECK( !peaks_mostly_sums( peaks, sums, sm_cluster_num_sigma, sm_max_share )[1] );
  sums = { { 300.0, 1.01*sm_max_share*1.0E4 } };
  BOOST_CHECK( peaks_mostly_sums( peaks, sums, sm_cluster_num_sigma, sm_max_share )[1] );

  // No sums modeled: every peak is used.
  const vector<bool> none = peaks_mostly_sums( peaks, {}, sm_cluster_num_sigma, sm_max_share );
  BOOST_REQUIRE_EQUAL( none.size(), peaks.size() );
  BOOST_CHECK( std::none_of( begin(none), end(none), []( const bool b ){ return b; } ) );
}//window_and_share


// A source whose every peak may be mostly sums keeps the one with its largest yield, so the estimate still
//  has something for it.
BOOST_AUTO_TEST_CASE( a_source_keeps_a_peak )
{
  const vector<RelActCalcManual::GenericPeakInfo> peaks{
    peak( 185.715, 1.0E6, {{"U235", 0.572}} ),
    peak( 258.26, 2.0E3, {{"U238", 7.3E-4}} ),
    peak( 1000.99, 5.0E2, {{"U238", 8.4E-3}} ),
    peak( 345.9, 3.0E3, {{"U235", 3.0E-4}} )
  };
  const vector<pair<double,double>> sums{ { 258.0, 1.5E3 }, { 1001.0, 4.0E2 }, { 345.9, 2.0E3 } };

  const vector<bool> leave_out = peaks_mostly_sums( peaks, sums, sm_cluster_num_sigma, sm_max_share );
  BOOST_REQUIRE_EQUAL( leave_out.size(), peaks.size() );
  BOOST_CHECK( !leave_out[0] );
  BOOST_CHECK( leave_out[1] );    //U-238's weaker-yield line goes
  BOOST_CHECK( !leave_out[2] );   //its 1001 keV line stays: otherwise U-238 would have no peak
  BOOST_CHECK( leave_out[3] );    //U-235 still has 185.7 keV
}//a_source_keeps_a_peak


namespace
{
/** A parameter layout like a single Physical-model curve's, with the energy calibration fit (offset, gain
 and deviation pairs), a correction function, and - if `with_sums` - two sum-peak parameters.
 */
StageGroups example_groups( const bool with_sums, const bool fit_energy_cal = true, const bool with_correction = true )
{
  StageGroups g;
  size_t index = 0;
  auto take = [&index]( const size_t n ) -> vector<size_t> {
    vector<size_t> answer;
    for( size_t i = 0; i < n; ++i )
      answer.push_back( index++ );
    return answer;
  };

  if( fit_energy_cal )
  {
    g.linear_cal = take( 2 );
    g.dev_pairs = take( 5 );
    g.energy_cal = g.linear_cal;
    g.energy_cal.insert( end(g.energy_cal), begin(g.dev_pairs), end(g.dev_pairs) );
  }
  g.fwhm = take( 3 );
  g.self_atten = take( 1 );
  g.ext_atten = take( 1 );
  if( with_correction )
    g.correction = take( 2 );
  const vector<size_t> activities = take( 4 );
  g.skew = take( 2 );
  if( with_sums )
  {
    g.sums = take( 2 );
    g.skew.insert( end(g.skew), begin(g.sums), end(g.sums) );
  }

  for( size_t i = 0; i < index; ++i )
    g.free.push_back( i );
  (void)activities;
  return g;
}//example_groups(...)


bool holds_all( const vector<size_t> &held, const vector<size_t> &group )
{
  return std::all_of( begin(group), end(group), [&held]( const size_t i ){
    return std::find( begin(held), end(held), i ) != end(held);
  } );
}

bool holds_none( const vector<size_t> &held, const vector<size_t> &group )
{
  return std::none_of( begin(group), end(group), [&held]( const size_t i ){
    return std::find( begin(held), end(held), i ) != end(held);
  } );
}

vector<size_t> joined( const vector<vector<size_t>> &lists )
{
  vector<size_t> answer;
  for( const vector<size_t> &list : lists )
    for( const size_t i : list )
      if( std::find( begin(answer), end(answer), i ) == end(answer) )
        answer.push_back( i );
  return answer;
}

vector<size_t> sorted( vector<size_t> v )
{
  std::sort( begin(v), end(v) );
  return v;
}
}//namespace


// Without sum peaks the stages are what they were before sum peaks were modeled.
BOOST_AUTO_TEST_CASE( stages_without_sums )
{
  {
    const StageGroups g = example_groups( false );
    const vector<vector<size_t>> stages = RelActCalcAutoImp::PreFit::stages( g );
    const vector<vector<size_t>> expected{
      joined( {g.energy_cal, g.fwhm, g.correction, g.self_atten, g.ext_atten, g.other_rel_eff, g.skew} ),
      joined( {g.correction, g.skew} ),
      joined( {g.dev_pairs, g.correction, g.skew} ),
      joined( {g.correction, g.skew} )
    };
    BOOST_REQUIRE_EQUAL( stages.size(), expected.size() );
    for( size_t i = 0; i < stages.size(); ++i )
      BOOST_CHECK( sorted( stages[i] ) == sorted( expected[i] ) );
  }

  {
    // No energy calibration fit and no correction function: the widths get a stage of their own.
    const StageGroups g = example_groups( false, false, false );
    const vector<vector<size_t>> stages = RelActCalcAutoImp::PreFit::stages( g );
    BOOST_REQUIRE_EQUAL( stages.size(), 2u );
    BOOST_CHECK( sorted( stages[1] ) == sorted( g.skew ) );
  }
}//stages_without_sums


// With sum peaks: one stage fits the sum parameters by themselves, and until then the deviation pairs are
//  held (a deviation at a ROI of sum peaks would otherwise run to its limit before the sums are fit).
BOOST_AUTO_TEST_CASE( sums_fit_alone_before_deviation_pairs )
{
  for( const bool with_correction : { true, false } )
  {
    const StageGroups g = example_groups( true, true, with_correction );
    const vector<vector<size_t>> stages = RelActCalcAutoImp::PreFit::stages( g );

    vector<size_t> all_but_sums;
    for( const size_t i : g.free )
      if( std::find( begin(g.sums), end(g.sums), i ) == end(g.sums) )
        all_but_sums.push_back( i );

    size_t sums_stage = stages.size(), first_dev_pair_stage = stages.size();
    size_t num_sums_stages = 0;
    for( size_t i = 0; i < stages.size(); ++i )
    {
      if( sorted( stages[i] ) == sorted( all_but_sums ) )
      {
        sums_stage = i;
        ++num_sums_stages;
      }else
      {
        BOOST_CHECK_MESSAGE( holds_all( stages[i], g.sums ), "stage " << i << " frees a sum parameter" );
      }

      if( !holds_all( stages[i], g.dev_pairs ) && (first_dev_pair_stage == stages.size()) )
        first_dev_pair_stage = i;
    }//for( loop over stages )

    BOOST_CHECK_EQUAL( num_sums_stages, 1u );
    BOOST_REQUIRE_LT( sums_stage, stages.size() );
    BOOST_CHECK_MESSAGE( first_dev_pair_stage > sums_stage,
                         "the deviation pairs are freed (stage " << first_dev_pair_stage
                         << ") before the sums are fit (stage " << sums_stage << ")" );
    BOOST_CHECK_LT( first_dev_pair_stage, stages.size() );   //but they are freed in the end
    BOOST_CHECK( holds_none( stages[stages.size() - 1], g.dev_pairs ) );
  }//for( with and without a correction function )

  {
    // Without the energy calibration fit, the sums are fit last.
    const StageGroups g = example_groups( true, false, true );
    const vector<vector<size_t>> stages = RelActCalcAutoImp::PreFit::stages( g );
    BOOST_REQUIRE( !stages.empty() );
    BOOST_CHECK( holds_none( stages.back(), g.sums ) );
    BOOST_CHECK_EQUAL( stages.back().size(), g.free.size() - g.sums.size() );
  }
}//sums_fit_alone_before_deviation_pairs
