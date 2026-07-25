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


BOOST_AUTO_TEST_CASE( blendedPointIsHeavilyDiscounted )
{
  // A cluster split 50/50 between two rel. eff. curves carries a ~50% uncertainty, because how the counts
  //  divide between curves is a model assumption rather than a measurement (co-located gammas are perfectly
  //  degenerate in the fit).  This mirrors the arithmetic in `fit_free_peak_amplitudes()`.
  const double area_rel_uncert = 0.02; //a 2% peak - by itself, 50 sigma
  const double curve_model_fraction = 0.5;
  const double blend_rel_uncert = 1.0 - curve_model_fraction;
  const double total_rel_uncert = sqrt( area_rel_uncert*area_rel_uncert + blend_rel_uncert*blend_rel_uncert );

  BOOST_CHECK_CLOSE( total_rel_uncert, 0.5004, 0.1 );

  // The point is now ~2 sigma, so it does not pass the significance gate, while the same peak on a curve
  //  that owns it outright would be 50 sigma and pass.
  ObsEff blended = good_point();
  blended.curve_model_fraction = curve_model_fraction;
  blended.curve_num_sigma_significance = 1.0/total_rel_uncert;
  BOOST_CHECK( !RelActCalcAuto::RelActAutoSolution::show_obs_eff_point( blended ) );
  BOOST_CHECK( blended.exclusion_reason == ExclusionReason::Insignificant );

  ObsEff unblended = good_point();
  unblended.curve_num_sigma_significance = 1.0/area_rel_uncert;
  BOOST_CHECK( RelActCalcAuto::RelActAutoSolution::show_obs_eff_point( unblended ) );

  // A lightly blended point (this curve owns 95%) is still shown.
  ObsEff light = good_point();
  light.curve_model_fraction = 0.95;
  const double light_uncert = sqrt( area_rel_uncert*area_rel_uncert + 0.05*0.05 );
  light.curve_num_sigma_significance = 1.0/light_uncert;
  BOOST_CHECK( RelActCalcAuto::RelActAutoSolution::show_obs_eff_point( light ) );
}
