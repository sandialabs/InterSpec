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
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include "InterSpec_config.h"

#include <cmath>
#include <limits>

#define BOOST_TEST_MODULE RelActCalcAuto_ObsEffFilter_suite
#include <boost/test/included/unit_test.hpp>

#include "InterSpec/PeakDef.h"
#include "InterSpec/RelActCalcAuto.h"

using namespace std;

static_assert( USE_REL_ACT_TOOL, "Compile-time support for the Rel. Act. tool is required for this test." );

using ObsEff = RelActCalcAuto::RelActAutoSolution::ObsEff;
using ExclusionReason = ObsEff::ExclusionReason;

namespace
{
/** An `ObsEff` that `show_obs_eff_point(...)` should accept; individual tests then break one thing. */
ObsEff good_point()
{
  ObsEff eff;
  eff.energy = 185.7;
  eff.fit_peaks.push_back( PeakDef( 185.7, 0.5, 1000.0 ) );
  eff.observed_efficiency = 0.75;
  eff.curve_model_fraction = 1.0;
  eff.curve_fit_amplitude = 1000.0;
  eff.curve_num_sigma_significance = 30.0;
  eff.neighbor_leak_fraction = 0.1;
  eff.within_roi = true;

  return eff;
}//good_point()
}//namespace


BOOST_AUTO_TEST_CASE( obsEffDefaultInitialized )
{
  // Every `ObsEff` member needs a default initializer: `fit_free_peak_amplitudes()` used to leave
  //  `cluster_upper_energy` and `roi_lower_energy` unwritten, which is UB for any consumer.
  const ObsEff eff;

  BOOST_CHECK_EQUAL( eff.energy, 0.0 );
  BOOST_CHECK_EQUAL( eff.cluster_lower_energy, 0.0 );
  BOOST_CHECK_EQUAL( eff.cluster_upper_energy, 0.0 );
  BOOST_CHECK_EQUAL( eff.roi_lower_energy, 0.0 );
  BOOST_CHECK_EQUAL( eff.roi_upper_energy, 0.0 );
  BOOST_CHECK_EQUAL( eff.resolution_sigma, 0.0 );
  BOOST_CHECK_EQUAL( eff.neighbor_leak_fraction, 0.0 );
  BOOST_CHECK_EQUAL( eff.curve_model_fraction, 1.0 ); //a lone curve owns its clusters outright
  BOOST_CHECK_EQUAL( eff.within_roi, false );
  BOOST_CHECK( eff.exclusion_reason == ExclusionReason::NotExcluded );
}


BOOST_AUTO_TEST_CASE( showObsEffPointAcceptsGoodPoint )
{
  ObsEff eff = good_point();
  BOOST_CHECK( RelActCalcAuto::RelActAutoSolution::show_obs_eff_point( eff ) );
  BOOST_CHECK( eff.exclusion_reason == ExclusionReason::NotExcluded );
}


BOOST_AUTO_TEST_CASE( showObsEffPointExclusionReasons )
{
  // No peaks from this rel. eff. curve in the cluster.
  {
    ObsEff eff = good_point();
    eff.fit_peaks.clear();
    BOOST_CHECK( !RelActCalcAuto::RelActAutoSolution::show_obs_eff_point( eff ) );
    BOOST_CHECK( eff.exclusion_reason == ExclusionReason::NoCountsForCurve );
  }

  {
    ObsEff eff = good_point();
    eff.curve_model_fraction = 0.0;
    BOOST_CHECK( !RelActCalcAuto::RelActAutoSolution::show_obs_eff_point( eff ) );
    BOOST_CHECK( eff.exclusion_reason == ExclusionReason::NoCountsForCurve );
  }

  // The free-amplitude fit went negative.
  {
    ObsEff eff = good_point();
    eff.observed_efficiency = -0.01;
    BOOST_CHECK( !RelActCalcAuto::RelActAutoSolution::show_obs_eff_point( eff ) );
    BOOST_CHECK( eff.exclusion_reason == ExclusionReason::NonPositiveEff );
  }

  {
    ObsEff eff = good_point();
    eff.curve_num_sigma_significance = 1.0;
    BOOST_CHECK( !RelActCalcAuto::RelActAutoSolution::show_obs_eff_point( eff ) );
    BOOST_CHECK( eff.exclusion_reason == ExclusionReason::Insignificant );
  }

  // A minor gamma riding the tail of a much larger neighbor - the case the omission logic exists for.
  {
    ObsEff eff = good_point();
    eff.neighbor_leak_fraction = 0.9;
    BOOST_CHECK( !RelActCalcAuto::RelActAutoSolution::show_obs_eff_point( eff ) );
    BOOST_CHECK( eff.exclusion_reason == ExclusionReason::TailLeakage );
  }

  {
    ObsEff eff = good_point();
    eff.within_roi = false;
    BOOST_CHECK( !RelActCalcAuto::RelActAutoSolution::show_obs_eff_point( eff ) );
    BOOST_CHECK( eff.exclusion_reason == ExclusionReason::OutsideRoi );
  }
}


BOOST_AUTO_TEST_CASE( showObsEffPointRejectsNaNs )
{
  // A NaN must exclude the point, not sneak it through; the predicate is written with negated comparisons
  //  so that every NaN comparison is false.
  const double nan_val = std::numeric_limits<double>::quiet_NaN();

  {
    ObsEff eff = good_point();
    eff.observed_efficiency = nan_val;
    BOOST_CHECK( !RelActCalcAuto::RelActAutoSolution::show_obs_eff_point( eff ) );
    BOOST_CHECK( eff.exclusion_reason == ExclusionReason::NonPositiveEff );
  }

  {
    ObsEff eff = good_point();
    eff.curve_num_sigma_significance = nan_val;
    BOOST_CHECK( !RelActCalcAuto::RelActAutoSolution::show_obs_eff_point( eff ) );
    BOOST_CHECK( eff.exclusion_reason == ExclusionReason::Insignificant );
  }

  {
    ObsEff eff = good_point();
    eff.neighbor_leak_fraction = nan_val;
    BOOST_CHECK( !RelActCalcAuto::RelActAutoSolution::show_obs_eff_point( eff ) );
    BOOST_CHECK( eff.exclusion_reason == ExclusionReason::TailLeakage );
  }
}


BOOST_AUTO_TEST_CASE( sharedPointSignificanceIsStatisticalTimesShare )
{
  // For a cluster shared between rel. eff. curves, a point's significance is the statistically-measured signal
  //  _attributed_ to this curve: `curve_num_sigma_significance = curve_model_fraction * (area/area_uncert)`
  //  (see `fit_free_peak_amplitudes()`).  The unmeasurable split between curves is a model assumption, so it is
  //  deliberately _not_ folded into the significance - the chart conveys it by fading the point's opacity instead.
  //  Hence a modest share of a strong peak is still shown, while any share of a weak peak is dropped.
  const double strong_cluster_sigma = 50.0; //a 2% peak - by itself, 50 sigma
  const double weak_cluster_sigma   = 4.0;

  // Half of a strong peak is still 25 sigma, so it is shown (the pre-blend-decoupling code dropped this).
  ObsEff strong_half = good_point();
  strong_half.curve_model_fraction = 0.5;
  strong_half.curve_num_sigma_significance = 0.5 * strong_cluster_sigma;
  BOOST_CHECK( RelActCalcAuto::RelActAutoSolution::show_obs_eff_point( strong_half ) );

  // Half of a weak peak is only 2 sigma, so it does not pass the significance gate.
  ObsEff weak_half = good_point();
  weak_half.curve_model_fraction = 0.5;
  weak_half.curve_num_sigma_significance = 0.5 * weak_cluster_sigma;
  BOOST_CHECK( !RelActCalcAuto::RelActAutoSolution::show_obs_eff_point( weak_half ) );
  BOOST_CHECK( weak_half.exclusion_reason == ExclusionReason::Insignificant );

  // A tiny share (4%) of even a strong peak eventually falls below the gate (here 2 sigma).
  ObsEff sliver = good_point();
  sliver.curve_model_fraction = 0.04;
  sliver.curve_num_sigma_significance = 0.04 * strong_cluster_sigma;
  BOOST_CHECK( !RelActCalcAuto::RelActAutoSolution::show_obs_eff_point( sliver ) );
  BOOST_CHECK( sliver.exclusion_reason == ExclusionReason::Insignificant );
}
