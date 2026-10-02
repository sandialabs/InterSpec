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

/* Unit tests for resolving RelActCalcAuto ROI specifications into the ROIs a fit uses
 (`RelActCalcAutoRoi::resolve(...)`), and for the ROI XML serialization.  No data files needed.
 */

#include "InterSpec_config.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>

#define BOOST_TEST_MODULE RelActCalcAuto_Roi_suite
#include <boost/test/included/unit_test.hpp>

#include <boost/math/distributions/normal.hpp>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakDists.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/PeakFitDetPrefs.h"
#include "InterSpec/RelActCalcAuto_Roi.h"

using namespace std;

using RelActCalcAuto::RoiEdge;
using RelActCalcAuto::RoiRange;
using RelActCalcAutoRoi::Line;
using RelActCalcAutoRoi::PeakShape;
using RelActCalcAutoRoi::ResolvedRoi;

namespace
{
  /** 8192 channels over ~3 MeV, like a typical HPGe. */
  shared_ptr<const SpecUtils::EnergyCalibration> hpge_cal()
  {
    auto cal = make_shared<SpecUtils::EnergyCalibration>();
    cal->set_polynomial( 8192, { 0.0f, 0.366f }, {} );
    return cal;
  }


  /** A Gaussian peak shape with FWHM = 0.8 + 0.001*E keV. */
  PeakShape gaussian_shape( const double energy )
  {
    PeakShape shape;
    shape.fwhm = 0.8 + 0.001*energy;
    shape.skew_type = PeakDef::SkewType::NoSkew;
    return shape;
  }


  RoiRange line_anchored( const double lower, const double upper,
                          const PeakContinuum::OffsetType cont = PeakContinuum::OffsetType::Linear )
  {
    RoiRange roi;
    roi.lower_energy = lower;
    roi.upper_energy = upper;
    roi.continuum_type = cont;
    roi.range_limits_type = RoiRange::RangeLimitsType::LineAnchored;
    return roi;
  }


  RoiRange fixed_roi( const double lower, const double upper )
  {
    RoiRange roi;
    roi.lower_energy = lower;
    roi.upper_energy = upper;
    roi.continuum_type = PeakContinuum::OffsetType::Linear;
    roi.range_limits_type = RoiRange::RangeLimitsType::Fixed;
    return roi;
  }


  /** The nearest-channel rounding RelActCalcAuto integrates a ROI's energy bounds with. */
  size_t nearest_channel( const SpecUtils::EnergyCalibration &cal, const double energy )
  {
    return static_cast<size_t>( std::floor( cal.channel_for_energy( energy ) + 0.5 ) );
  }


  /** The default ROI edges to resolve with (by default, the `RoiEdge` defaults). */
  struct EdgeDefaults
  {
    RoiEdge lower_edge, upper_edge;
  };


  vector<ResolvedRoi> resolve( const vector<RoiRange> &rois, const vector<Line> &lines,
                               vector<string> &warnings,
                               const EdgeDefaults &settings = {},
                               const vector<RelActCalcAutoRoi::SidebandObstacle> &obstacles = {} )
  {
    return RelActCalcAutoRoi::resolve( rois, settings.lower_edge, settings.upper_edge, gaussian_shape,
                                       nullptr, hpge_cal(), lines, obstacles, warnings );
  }


  /** Checks the invariants every resolution must satisfy: sorted, fixed, and no two ROIs sharing a
   channel under the nearest-channel convention RelActCalcAuto uses; ROIs we resolved (i.e., not
   input fixed ROIs, whose bounds are used as given) must not share one under the floor convention
   either. */
  void check_invariants( const vector<ResolvedRoi> &rois )
  {
    const shared_ptr<const SpecUtils::EnergyCalibration> cal = hpge_cal();
    for( size_t i = 0; i < rois.size(); ++i )
    {
      BOOST_CHECK( rois[i].roi.range_limits_type == RoiRange::RangeLimitsType::Fixed );
      BOOST_CHECK( !rois[i].roi.auto_continuum );
      BOOST_CHECK_LT( rois[i].roi.lower_energy, rois[i].roi.upper_energy );

      if( i == 0 )
        continue;

      BOOST_CHECK_LT( rois[i-1].roi.upper_energy, rois[i].roi.lower_energy );
      BOOST_CHECK_LT( nearest_channel( *cal, rois[i-1].roi.upper_energy ),
                      nearest_channel( *cal, rois[i].roi.lower_energy ) );

      const bool both_resolved = !std::isnan( rois[i-1].info.lower_anchor_energy )
                                 && !std::isnan( rois[i].info.lower_anchor_energy );
      if( both_resolved )
        BOOST_CHECK_LT( std::floor( cal->channel_for_energy( rois[i-1].roi.upper_energy ) ),
                        std::floor( cal->channel_for_energy( rois[i].roi.lower_energy ) ) );
    }
  }//check_invariants(...)
}//namespace


BOOST_AUTO_TEST_CASE( gaussian_edge_distance_matches_quantile )
{
  const boost::math::normal_distribution<double> unit_normal( 0.0, 1.0 );

  for( const double tail : { 1.0E-2, 1.0E-3, 1.0E-4 } )
  {
    for( const double sideband : { 0.0, 1.0, 2.5 } )
    {
      const RoiEdge edge{ tail, sideband };
      const PeakShape shape = gaussian_shape( 500.0 );
      const double sigma = shape.fwhm / 2.35482;
      const double expected = sigma*boost::math::quantile( boost::math::complement(unit_normal, tail) )
                              + sideband*shape.fwhm;

      BOOST_CHECK_CLOSE( RelActCalcAutoRoi::edge_distance( 500.0, shape, edge, true ), expected, 1.0E-6 );
      BOOST_CHECK_CLOSE( RelActCalcAutoRoi::edge_distance( 500.0, shape, edge, false ), expected, 1.0E-6 );
    }
  }

  BOOST_CHECK_THROW( RelActCalcAutoRoi::edge_distance( 500.0, gaussian_shape(500.0), RoiEdge{0.0, 1.0}, true ),
                     std::exception );
  BOOST_CHECK_THROW( RelActCalcAutoRoi::edge_distance( 500.0, gaussian_shape(500.0), RoiEdge{0.6, 1.0}, true ),
                     std::exception );
  BOOST_CHECK_THROW( RelActCalcAutoRoi::edge_distance( 500.0, gaussian_shape(500.0), RoiEdge{1.0E-3, -1.0}, true ),
                     std::exception );
}//gaussian_edge_distance_matches_quantile


BOOST_AUTO_TEST_CASE( skewed_edge_distance_leaves_requested_tail )
{
  // GaussExp has an exponential low-energy tail, so the low side must extend further than the high side,
  //  and exactly `tail_fraction` of the area must lie beyond each edge's coverage point.
  const double energy = 185.7;
  PeakShape shape;
  shape.fwhm = 1.2;
  shape.skew_type = PeakDef::SkewType::GaussExp;
  shape.skew_pars[0] = 1.2;
  const double sigma = shape.fwhm / 2.35482;

  const double tail = 1.0E-3;
  const RoiEdge edge{ tail, 0.0 };
  const double lower_dist = RelActCalcAutoRoi::edge_distance( energy, shape, edge, true );
  const double upper_dist = RelActCalcAutoRoi::edge_distance( energy, shape, edge, false );
  BOOST_CHECK_GT( lower_dist, upper_dist );

  const double total = PeakDists::gauss_exp_integral( energy, sigma, shape.skew_pars[0], energy - 200.0*sigma, energy + 200.0*sigma );
  const double below = PeakDists::gauss_exp_integral( energy, sigma, shape.skew_pars[0], energy - 200.0*sigma, energy - lower_dist );
  const double above = PeakDists::gauss_exp_integral( energy, sigma, shape.skew_pars[0], energy + upper_dist, energy + 200.0*sigma );
  BOOST_CHECK_CLOSE( below/total, tail, 5.0 );
  BOOST_CHECK_CLOSE( above/total, tail, 5.0 );
}//skewed_edge_distance_leaves_requested_tail


BOOST_AUTO_TEST_CASE( merged_continuum_type_rule )
{
  using OT = PeakContinuum::OffsetType;
  const auto check = []( const OT a, const OT b, const OT expected ){
    BOOST_CHECK_EQUAL( static_cast<int>(RelActCalcAutoRoi::merged_continuum_type(a, b)), static_cast<int>(expected) );
    BOOST_CHECK_EQUAL( static_cast<int>(RelActCalcAutoRoi::merged_continuum_type(b, a)), static_cast<int>(expected) );
  };

  check( OT::Linear, OT::Linear, OT::Linear );
  check( OT::Linear, OT::Quadratic, OT::Quadratic );
  check( OT::Constant, OT::Cubic, OT::Cubic );
  check( OT::Linear, OT::FlatStep, OT::LinearStep );
  check( OT::FlatStep, OT::FlatStepCDF, OT::FlatStepCDF );
  check( OT::Linear, OT::FlatStepCDF, OT::LinearStepCDF );
  check( OT::Quadratic, OT::FlatStep, OT::BiLinearStep );
  check( OT::Cubic, OT::LinearStepCDF, OT::BiLinearStepCDF );
  check( OT::BiLinearStep, OT::Constant, OT::BiLinearStep );
  check( OT::External, OT::Quadratic, OT::Quadratic );
}//merged_continuum_type_rule


BOOST_AUTO_TEST_CASE( fixed_rois_pass_through_unchanged )
{
  const vector<RoiRange> input{ fixed_roi( 760.0, 770.0 ), fixed_roi( 180.7, 188.8 ), fixed_roi( 997.1, 1005.9 ) };
  vector<string> warnings;
  const vector<ResolvedRoi> resolved = resolve( input, {}, warnings );

  BOOST_REQUIRE_EQUAL( resolved.size(), 3 );
  BOOST_CHECK( warnings.empty() );
  BOOST_CHECK( resolved[0].roi == input[1] );
  BOOST_CHECK( resolved[1].roi == input[0] );
  BOOST_CHECK( resolved[2].roi == input[2] );
  BOOST_CHECK( resolved[0].info.input_roi_indices == vector<size_t>{1} );
  BOOST_CHECK( std::isnan( resolved[0].info.lower_anchor_energy ) );

  // Overlapping fixed ROIs are an error, as is a ROI with its bounds reversed.
  vector<string> w2;
  BOOST_CHECK_THROW( resolve( { fixed_roi(100.0, 120.0), fixed_roi(119.0, 130.0) }, {}, w2 ), std::exception );
  BOOST_CHECK_THROW( resolve( { fixed_roi(120.0, 100.0) }, {}, w2 ), std::exception );

  // A fixed ROI mostly outside the spectrum is left out, with a warning.
  vector<string> w3;
  const vector<ResolvedRoi> outside = resolve( { fixed_roi(2990.0, 3050.0), fixed_roi(500.0, 510.0) }, {}, w3 );
  BOOST_REQUIRE_EQUAL( outside.size(), 1 );
  BOOST_CHECK_EQUAL( w3.size(), 1 );
}//fixed_rois_pass_through_unchanged


BOOST_AUTO_TEST_CASE( line_anchored_roi_extends_by_peak_shape )
{
  const shared_ptr<const SpecUtils::EnergyCalibration> cal = hpge_cal();
  EdgeDefaults settings;
  settings.lower_edge = RoiEdge{ 1.0E-3, 1.0 };
  settings.upper_edge = RoiEdge{ 1.0E-4, 0.5 };

  // A single line, and a pair of lines.
  const vector<RoiRange> input{ line_anchored( 1000.99, 1000.99 ), line_anchored( 182.52, 185.71, PeakContinuum::OffsetType::FlatStep ) };
  vector<string> warnings;
  const vector<ResolvedRoi> resolved = resolve( input, {}, warnings, settings );
  BOOST_REQUIRE_EQUAL( resolved.size(), 2 );
  BOOST_CHECK( warnings.empty() );
  check_invariants( resolved );

  const ResolvedRoi &pair_roi = resolved[0], &single_roi = resolved[1];
  BOOST_CHECK( pair_roi.info.input_roi_indices == vector<size_t>{1} );
  BOOST_CHECK_EQUAL( pair_roi.info.lower_anchor_energy, 182.52 );
  BOOST_CHECK_EQUAL( pair_roi.info.upper_anchor_energy, 185.71 );
  BOOST_CHECK( pair_roi.roi.continuum_type == PeakContinuum::OffsetType::FlatStep );

  for( const ResolvedRoi *r : { &pair_roi, &single_roi } )
  {
    const double lower_line = r->info.lower_anchor_energy, upper_line = r->info.upper_anchor_energy;
    const double want_lower = lower_line - RelActCalcAutoRoi::edge_distance( lower_line, gaussian_shape(lower_line), settings.lower_edge, true );
    const double want_upper = upper_line + RelActCalcAutoRoi::edge_distance( upper_line, gaussian_shape(upper_line), settings.upper_edge, false );

    // Bounds are placed a quarter of the way into the channel containing the wanted edge.
    BOOST_CHECK_EQUAL( std::floor( cal->channel_for_energy( r->roi.lower_energy ) ),
                       std::floor( cal->channel_for_energy( want_lower ) ) );
    BOOST_CHECK_EQUAL( std::floor( cal->channel_for_energy( r->roi.upper_energy ) ),
                       std::floor( cal->channel_for_energy( want_upper ) ) );
    for( const double bound : { r->roi.lower_energy, r->roi.upper_energy } )
    {
      const double frac = cal->channel_for_energy( bound ) - std::floor( cal->channel_for_energy( bound ) );
      BOOST_CHECK_CLOSE( frac, 0.25, 1.0E-3 );
    }
  }

  // Above the line, the higher coverage (0.01% tail vs 0.1%) does not make up for the smaller sideband
  //  (0.5 vs 1 FWHM), so the ROI extends less far above the line than below it.
  BOOST_CHECK_LT( single_roi.roi.upper_energy - 1000.99, 1000.99 - single_roi.roi.lower_energy );
}//line_anchored_roi_extends_by_peak_shape


BOOST_AUTO_TEST_CASE( overlapping_rois_merge_and_fixed_rois_win )
{
  // Two line-anchored ROIs whose extents overlap become one ROI covering both (at ~1 keV FWHM, each
  //  extends ~2.3 keV past its line by default).
  vector<string> warnings;
  const vector<ResolvedRoi> merged = resolve( { line_anchored( 202.11, 205.31, PeakContinuum::OffsetType::FlatStep ),
                                                line_anchored( 200.0, 200.0, PeakContinuum::OffsetType::Linear ) },
                                              {}, warnings );
  BOOST_REQUIRE_EQUAL( merged.size(), 1 );
  check_invariants( merged );
  BOOST_CHECK( (merged[0].info.input_roi_indices == vector<size_t>{0, 1}) );
  BOOST_CHECK_EQUAL( merged[0].info.lower_anchor_energy, 200.0 );
  BOOST_CHECK_EQUAL( merged[0].info.upper_anchor_energy, 205.31 );
  BOOST_CHECK( merged[0].roi.continuum_type == PeakContinuum::OffsetType::LinearStep );

  // But lines further apart than the peak widths stay separate ROIs.
  vector<string> w1;
  const vector<ResolvedRoi> separate = resolve( { line_anchored( 202.11, 205.31 ), line_anchored( 194.94, 194.94 ) }, {}, w1 );
  BOOST_CHECK_EQUAL( separate.size(), 2 );
  check_invariants( separate );

  // A fixed ROI next to a line-anchored ROI keeps its bounds, and the line-anchored ROI is limited to
  //  not share a channel with it.
  vector<string> w2;
  const RoiRange fixed = fixed_roi( 187.0, 192.0 );
  const vector<ResolvedRoi> with_fixed = resolve( { line_anchored( 185.71, 185.71 ), fixed }, {}, w2 );
  BOOST_REQUIRE_EQUAL( with_fixed.size(), 2 );
  check_invariants( with_fixed );
  BOOST_CHECK( with_fixed[1].roi == fixed );
  BOOST_CHECK( w2.empty() );  // only the continuum sideband was cut, not a line

  // A fixed ROI over the line itself is worth a warning.
  vector<string> w3;
  const vector<ResolvedRoi> over_line = resolve( { line_anchored( 185.71, 185.71 ), fixed_roi( 184.0, 187.0 ) }, {}, w3 );
  check_invariants( over_line );
  BOOST_CHECK_GE( w3.size(), 1 );
}//overlapping_rois_merge_and_fixed_rois_win


BOOST_AUTO_TEST_CASE( overlapping_sidebands_are_shared )
{
  // At ~1 keV FWHM, the peaks at 200 and 204 keV each cover ~1.31 keV to either side (0.1% tails),
  //  leaving 1.38 keV between them - less than their combined 1 FWHM sidebands, but more than half of
  //  that; so they stay separate ROIs, which divide the continuum between them.
  const shared_ptr<const SpecUtils::EnergyCalibration> cal = hpge_cal();
  vector<string> warnings;
  const vector<ResolvedRoi> shared = resolve( { line_anchored( 200.0, 200.0 ), line_anchored( 204.0, 204.0 ) },
                                              {}, warnings );
  BOOST_REQUIRE_EQUAL( shared.size(), 2 );
  check_invariants( shared );
  BOOST_CHECK( warnings.empty() );

  const double lower_cover = 200.0 + RelActCalcAutoRoi::edge_distance( 200.0, gaussian_shape(200.0), RoiEdge{1.0E-3, 0.0}, false );
  const double upper_cover = 204.0 - RelActCalcAutoRoi::edge_distance( 204.0, gaussian_shape(204.0), RoiEdge{1.0E-3, 0.0}, true );
  const double boundary = 0.5*(lower_cover + upper_cover);
  const double channel_width = 0.366;
  BOOST_CHECK_GE( shared[0].roi.upper_energy, lower_cover );
  BOOST_CHECK_LE( shared[1].roi.lower_energy, upper_cover );
  BOOST_CHECK_SMALL( shared[0].roi.upper_energy - boundary, channel_width );
  BOOST_CHECK_SMALL( shared[1].roi.lower_energy - boundary, 1.5*channel_width );
  BOOST_CHECK_EQUAL( std::floor( cal->channel_for_energy( shared[0].roi.upper_energy ) ) + 1.0,
                     std::floor( cal->channel_for_energy( shared[1].roi.lower_energy ) ) );

  // With 2 FWHM sidebands the gap is less than half the combined sidebands, so they are merged.
  EdgeDefaults wide;
  wide.lower_edge.sideband_fwhm = wide.upper_edge.sideband_fwhm = 2.0;
  vector<string> w2;
  const vector<ResolvedRoi> merged = resolve( { line_anchored( 200.0, 200.0 ), line_anchored( 204.0, 204.0 ) },
                                              {}, w2, wide );
  BOOST_REQUIRE_EQUAL( merged.size(), 1 );
  check_invariants( merged );
  BOOST_CHECK( (merged[0].info.input_roi_indices == vector<size_t>{0, 1}) );
}//overlapping_sidebands_are_shared


BOOST_AUTO_TEST_CASE( sidebands_stop_short_of_other_peaks )
{
  // A peak found in the spectrum a few keV below a 120.9 keV line (like the uranium K-beta x-rays below
  //  the U234 line), where the ROI's 3 FWHM continuum sideband (to ~116.9 keV) would otherwise reach.
  EdgeDefaults settings;
  settings.lower_edge.sideband_fwhm = settings.upper_edge.sideband_fwhm = 3.0;
  const RoiRange roi = line_anchored( 120.9, 120.9 );
  const RelActCalcAutoRoi::SidebandObstacle xray{ 116.0, 118.0 };
  const RelActCalcAutoRoi::SidebandObstacle own_peak{ 120.9 - 1.1, 120.9 + 1.1 };

  vector<string> warnings;
  const vector<ResolvedRoi> free_roi = resolve( { roi }, {}, warnings, settings );
  const vector<ResolvedRoi> blocked = resolve( { roi }, {}, warnings, settings, { xray, own_peak } );
  BOOST_REQUIRE_EQUAL( free_roi.size(), 1 );
  BOOST_REQUIRE_EQUAL( blocked.size(), 1 );
  check_invariants( blocked );
  BOOST_CHECK( warnings.empty() );

  // The lower sideband stops at the x-ray peak; the upper side is not affected - including by the
  //  obstacle for the ROI's own peak.
  BOOST_CHECK_LT( free_roi[0].roi.lower_energy, xray.upper );
  BOOST_CHECK_SMALL( blocked[0].roi.lower_energy - xray.upper, 0.366 );
  BOOST_CHECK_EQUAL( blocked[0].roi.upper_energy, free_roi[0].roi.upper_energy );

  // An obstacle reaching into the region the peak covers only removes the sideband, not that region.
  const double coverage = RelActCalcAutoRoi::edge_distance( 120.9, gaussian_shape(120.9), RoiEdge{1.0E-3, 0.0}, true );
  const RelActCalcAutoRoi::SidebandObstacle close{ 120.9 - coverage - 2.0, 120.9 - coverage + 0.3 };
  const vector<ResolvedRoi> tight = resolve( { roi }, {}, warnings, settings, { close } );
  BOOST_REQUIRE_EQUAL( tight.size(), 1 );
  BOOST_CHECK_SMALL( tight[0].roi.lower_energy - (120.9 - coverage), 0.366 );
}//sidebands_stop_short_of_other_peaks


BOOST_AUTO_TEST_CASE( broken_up_roi_uses_significant_lines )
{
  RoiRange range;
  range.lower_energy = 100.0;
  range.upper_energy = 1000.0;
  range.continuum_type = PeakContinuum::OffsetType::Linear;
  range.range_limits_type = RoiRange::RangeLimitsType::CanBeBrokenUp;

  const vector<Line> lines{ { 143.76, true }, { 163.33, true }, { 185.71, true }, { 186.2, true },
                            { 400.0, false }, { 999.0, true }, { 1001.03, true }, { 50.0, true } };
  vector<string> warnings;
  const vector<ResolvedRoi> resolved = resolve( { range }, lines, warnings );
  check_invariants( resolved );

  // 143.76, 163.33, 185.71 + 186.2 (merged), and 999 keV - the insignificant 400 keV line and the lines
  //  outside the range get no ROI, and the window around 999 keV is limited to the range.
  BOOST_REQUIRE_EQUAL( resolved.size(), 4 );
  BOOST_CHECK_EQUAL( resolved[2].info.lower_anchor_energy, 185.71 );
  BOOST_CHECK_EQUAL( resolved[2].info.upper_anchor_energy, 186.2 );
  BOOST_CHECK_EQUAL( resolved[3].info.lower_anchor_energy, 999.0 );
  BOOST_CHECK_LE( resolved[3].roi.upper_energy, 1000.0 + 0.366 );

  // No significant lines at all gives no ROIs, and a warning.
  vector<string> w2;
  BOOST_CHECK( resolve( { range }, { { 400.0, false } }, w2 ).empty() );
  BOOST_CHECK_EQUAL( w2.size(), 1 );
}//broken_up_roi_uses_significant_lines


BOOST_AUTO_TEST_CASE( floating_peak_next_to_line_anchored_roi_is_covered )
{
  // A floating peak (e.g., for the Bi-214 768.36 keV background line) just above a line-anchored ROI for
  //  the 763.13 to 766.37 keV lines.  At 0.9 keV FWHM the ROI alone ends below 768.36 keV (it would
  //  not at the ~1.6 keV FWHM of `gaussian_shape`), so the ROI is extended to cover the floating peak
  //  like one of its own lines.
  const shared_ptr<const SpecUtils::EnergyCalibration> cal = hpge_cal();
  const auto sharp_shape = []( const double ) -> PeakShape {
    PeakShape shape;
    shape.fwhm = 0.9;
    return shape;
  };

  EdgeDefaults settings;
  settings.lower_edge = settings.upper_edge = RoiEdge{ 1.0E-3, 0.5 };
  const vector<RoiRange> input{ line_anchored( 763.13, 766.37 ) };

  const auto resolve_with = [&]( const vector<RoiRange> &rois, const vector<Line> &lines,
                                 const std::function<double(double)> &line_position = nullptr ) {
    vector<string> warnings;
    const vector<ResolvedRoi> answer = RelActCalcAutoRoi::resolve( rois, settings.lower_edge, settings.upper_edge,
                                                                   sharp_shape, line_position,
                                                                   cal, lines, {}, warnings );
    BOOST_CHECK( warnings.empty() );
    check_invariants( answer );
    return answer;
  };

  const Line floating{ 768.36, true, true, false };
  const vector<ResolvedRoi> alone = resolve_with( input, {} );
  const vector<ResolvedRoi> with = resolve_with( input, { floating } );
  BOOST_REQUIRE_EQUAL( alone.size(), 1 );
  BOOST_REQUIRE_EQUAL( with.size(), 1 );
  BOOST_CHECK_LT( alone[0].roi.upper_energy, floating.energy );

  const double want_upper = floating.energy + RelActCalcAutoRoi::edge_distance( floating.energy,
                                                   sharp_shape(floating.energy), settings.upper_edge, false );
  BOOST_CHECK_EQUAL( std::floor( cal->channel_for_energy( with[0].roi.upper_energy ) ),
                     std::floor( cal->channel_for_energy( want_upper ) ) );
  BOOST_CHECK_EQUAL( with[0].roi.lower_energy, alone[0].roi.lower_energy );
  BOOST_CHECK( with[0].info.input_roi_indices == vector<size_t>{0} );
  BOOST_CHECK_EQUAL( with[0].info.lower_anchor_energy, 763.13 );
  BOOST_CHECK_EQUAL( with[0].info.upper_anchor_energy, 768.36 );

  // A floating peak no ROI reaches gets no ROI (RelActCalcAuto reports it as not being in a ROI).
  const vector<ResolvedRoi> far = resolve_with( input, { Line{ 800.0, true, true, false } } );
  BOOST_REQUIRE_EQUAL( far.size(), 1 );
  BOOST_CHECK( far[0].roi == alone[0].roi );

  // One within a fixed ROI is left to it; the line-anchored ROI is not extended into the fixed one.
  const RoiRange fixed = fixed_roi( 768.1, 772.0 );
  const vector<ResolvedRoi> with_fixed = resolve_with( { input[0], fixed }, { Line{ 769.0, true, true, false } } );
  BOOST_REQUIRE_EQUAL( with_fixed.size(), 2 );
  BOOST_CHECK( with_fixed[0].roi == alone[0].roi );
  BOOST_CHECK( with_fixed[1].roi == fixed );

  // A floating peak with a known energy moves with the lines under a fitted energy calibration; one
  //  whose energy was read off the spectrum does not.
  const auto shifted = []( const double energy ) -> double { return energy + 1.5; };
  const vector<ResolvedRoi> known = resolve_with( input, { floating }, shifted );
  const vector<ResolvedRoi> observed = resolve_with( input, { Line{ 768.36, true, true, true } }, shifted );
  BOOST_REQUIRE_EQUAL( known.size(), 1 );
  BOOST_REQUIRE_EQUAL( observed.size(), 1 );
  BOOST_CHECK_EQUAL( std::floor( cal->channel_for_energy( known[0].roi.upper_energy ) ),
                     std::floor( cal->channel_for_energy( want_upper + 1.5 ) ) );
  BOOST_CHECK_EQUAL( std::floor( cal->channel_for_energy( observed[0].roi.upper_energy ) ),
                     std::floor( cal->channel_for_energy( want_upper ) ) );
}//floating_peak_next_to_line_anchored_roi_is_covered


BOOST_AUTO_TEST_CASE( roi_covers_line_when_energy_adjustment_went_astray )
{
  // A fitted energy-calibration adjustment that went astray (e.g., a deviation pair at its limit where
  //  few lines constrain it) puts the 867.38 keV line 5.4 keV low; the ROI must still cover the line
  //  itself, as well as where the adjustment puts it - so the next fit can see the peak, and correct.
  const shared_ptr<const SpecUtils::EnergyCalibration> cal = hpge_cal();
  const auto astray = []( const double energy ) -> double {
    return (std::fabs( energy - 867.38 ) < 1.0) ? (energy - 5.4) : energy;
  };

  vector<string> warnings;
  const vector<ResolvedRoi> nominal = RelActCalcAutoRoi::resolve( { line_anchored( 867.38, 867.38 ) }, {}, {},
                                                   gaussian_shape, nullptr, cal, {}, {}, warnings );
  const vector<ResolvedRoi> shifted = RelActCalcAutoRoi::resolve( { line_anchored( 867.38, 867.38 ) }, {}, {},
                                                   gaussian_shape, astray, cal, {}, {}, warnings );
  BOOST_REQUIRE_EQUAL( nominal.size(), 1 );
  BOOST_REQUIRE_EQUAL( shifted.size(), 1 );
  BOOST_CHECK( warnings.empty() );
  check_invariants( shifted );

  // Same upper edge as without the adjustment (the line's own position), and extended below to cover
  //  where the adjustment puts the line.
  BOOST_CHECK_EQUAL( shifted[0].roi.upper_energy, nominal[0].roi.upper_energy );
  const double coverage = RelActCalcAutoRoi::edge_distance( 867.38, gaussian_shape(867.38), RoiEdge{1.0E-3, 0.0}, true );
  const double sideband = gaussian_shape(867.38).fwhm;
  BOOST_CHECK_EQUAL( std::floor( cal->channel_for_energy( shifted[0].roi.lower_energy ) ),
                     std::floor( cal->channel_for_energy( 867.38 - 5.4 - coverage - sideband ) ) );

  // "Split by lines" windows behave the same way.
  RoiRange range;
  range.lower_energy = 800.0;
  range.upper_energy = 900.0;
  range.continuum_type = PeakContinuum::OffsetType::Linear;
  range.range_limits_type = RoiRange::RangeLimitsType::CanBeBrokenUp;
  const vector<ResolvedRoi> split = RelActCalcAutoRoi::resolve( { range }, {}, {}, gaussian_shape, astray, cal,
                                                                { Line{ 867.38, true } }, {}, warnings );
  BOOST_REQUIRE_EQUAL( split.size(), 1 );
  BOOST_CHECK_LT( split[0].roi.lower_energy, 867.38 - 5.4 );
  BOOST_CHECK_GT( split[0].roi.upper_energy, 867.38 );
}//roi_covers_line_when_energy_adjustment_went_astray


BOOST_AUTO_TEST_CASE( very_wide_line_anchored_roi_warns )
{
  // A "Lines + peak width" ROI spanning most of the spectrum is fit as one ROI, which its continuum can
  //  not follow - the user almost certainly wanted "Split by lines", so is told so; ordinary multi-line
  //  ROIs are not warned about.
  const auto split_advice = []( const vector<string> &warnings ) -> bool {
    for( const string &w : warnings )
    {
      if( w.find( "Split by lines" ) != string::npos )
        return true;
    }
    return false;
  };

  vector<string> wide_warnings;
  const vector<ResolvedRoi> wide = resolve( { line_anchored( 100.0, 2900.0 ) }, {}, wide_warnings );
  BOOST_REQUIRE_EQUAL( wide.size(), 1 );
  BOOST_CHECK( split_advice( wide_warnings ) );

  vector<string> normal_warnings;
  const vector<ResolvedRoi> normal = resolve( { line_anchored( 182.52, 185.71 ), line_anchored( 742.81, 786.25 ) },
                                              {}, normal_warnings );
  BOOST_REQUIRE_EQUAL( normal.size(), 2 );
  BOOST_CHECK( !split_advice( normal_warnings ) );
}//very_wide_line_anchored_roi_warns


BOOST_AUTO_TEST_CASE( lines_outside_spectrum_are_dropped )
{
  vector<string> warnings;
  const vector<ResolvedRoi> resolved = resolve( { line_anchored( 3500.0, 3500.0 ), line_anchored( 661.66, 661.66 ) }, {}, warnings );
  BOOST_REQUIRE_EQUAL( resolved.size(), 1 );
  BOOST_CHECK_EQUAL( resolved[0].info.lower_anchor_energy, 661.66 );
  BOOST_CHECK_EQUAL( warnings.size(), 1 );

  // Deterministic: resolving again gives exactly the same answer.
  vector<string> w2;
  const vector<ResolvedRoi> again = resolve( { line_anchored( 3500.0, 3500.0 ), line_anchored( 661.66, 661.66 ) }, {}, w2 );
  BOOST_REQUIRE_EQUAL( again.size(), 1 );
  BOOST_CHECK( again[0].roi == resolved[0].roi );
  BOOST_CHECK( again[0].info == resolved[0].info );
}//lines_outside_spectrum_are_dropped


BOOST_AUTO_TEST_CASE( lines_partly_outside_spectrum_are_clipped )
{
  // The `hpge_cal()` spectrum ends at 2998.3 keV; the 2990 keV line is in it, the 3010 keV one is not.
  const shared_ptr<const SpecUtils::EnergyCalibration> cal = hpge_cal();
  const double spectrum_upper = cal->energy_for_channel( static_cast<double>( cal->num_channels() ) );
  vector<string> warnings;
  const vector<ResolvedRoi> resolved = resolve( { line_anchored( 2990.0, 3010.0 ) }, {}, warnings );
  BOOST_REQUIRE_EQUAL( resolved.size(), 1 );
  BOOST_CHECK( warnings.empty() );
  check_invariants( resolved );
  BOOST_CHECK_LT( resolved[0].roi.lower_energy, 2990.0 );
  BOOST_CHECK_LE( resolved[0].roi.upper_energy, spectrum_upper );
  BOOST_CHECK_GT( resolved[0].roi.upper_energy, spectrum_upper - 1.0 );
}//lines_partly_outside_spectrum_are_clipped


BOOST_AUTO_TEST_CASE( unusual_energy_calibrations )
{
  vector<string> warnings;

  // Lower-channel-energy calibrations throw for energies above the spectrum, so a ROI reaching the top
  //  of the spectrum must not ask for one.
  {
    const size_t nchannel = 1024;
    vector<float> energies;
    for( size_t i = 0; i <= nchannel; ++i )
      energies.push_back( 3.0f*i );
    auto cal = make_shared<SpecUtils::EnergyCalibration>();
    cal->set_lower_channel_energy( nchannel, energies );

    vector<ResolvedRoi> resolved;
    BOOST_REQUIRE_NO_THROW( resolved = RelActCalcAutoRoi::resolve( { line_anchored( 3070.0, 3070.0 ) },
                                         RoiEdge{}, RoiEdge{}, gaussian_shape, nullptr, cal, {}, {}, warnings ) );
    BOOST_REQUIRE_EQUAL( resolved.size(), 1 );
    BOOST_CHECK_LE( resolved[0].roi.upper_energy, 3.0*nchannel );
    BOOST_CHECK_LT( resolved[0].roi.lower_energy, 3070.0 );
  }

  // With a negative energy offset, a ROI near the bottom of the spectrum is not taken below zero keV.
  {
    auto cal = make_shared<SpecUtils::EnergyCalibration>();
    cal->set_polynomial( 1024, { -8.0f, 3.0f }, {} );
    const auto wide_shape = []( const double ) -> PeakShape { PeakShape shape; shape.fwhm = 3.0; return shape; };

    vector<ResolvedRoi> resolved;
    BOOST_REQUIRE_NO_THROW( resolved = RelActCalcAutoRoi::resolve( { line_anchored( 6.0, 6.0 ) },
                                         RoiEdge{}, RoiEdge{}, wide_shape, nullptr, cal, {}, {}, warnings ) );
    BOOST_REQUIRE_EQUAL( resolved.size(), 1 );
    BOOST_CHECK_GE( resolved[0].roi.lower_energy, 0.0 );

    // ... so the resolved ROI is itself a valid (fixed) input ROI.
    vector<RoiRange> as_input{ resolved[0].roi };
    BOOST_CHECK_NO_THROW( RelActCalcAutoRoi::resolve( as_input, RoiEdge{}, RoiEdge{}, wide_shape,
                                                       nullptr, cal, {}, {}, warnings ) );
  }
}//unusual_energy_calibrations


BOOST_AUTO_TEST_CASE( windows_are_not_merged_through_fixed_roi )
{
  // Line-anchored ROIs on either side of a narrow fixed ROI, whose extents would overlap through it,
  //  stay separate ROIs - each with its own provenance and continuum.
  vector<string> warnings;
  const RoiRange fixed = fixed_roi( 186.0, 186.6 );
  const vector<ResolvedRoi> resolved = resolve( { line_anchored( 185.0, 185.0, PeakContinuum::OffsetType::Linear ),
                                                  fixed,
                                                  line_anchored( 187.6, 187.6, PeakContinuum::OffsetType::FlatStep ) },
                                                {}, warnings );
  BOOST_REQUIRE_EQUAL( resolved.size(), 3 );
  check_invariants( resolved );
  BOOST_CHECK( resolved[0].info.input_roi_indices == vector<size_t>{0} );
  BOOST_CHECK( resolved[0].roi.continuum_type == PeakContinuum::OffsetType::Linear );
  BOOST_CHECK_EQUAL( resolved[0].info.lower_anchor_energy, 185.0 );
  BOOST_CHECK( resolved[1].roi == fixed );
  BOOST_CHECK( resolved[2].info.input_roi_indices == vector<size_t>{2} );
  BOOST_CHECK( resolved[2].roi.continuum_type == PeakContinuum::OffsetType::FlatStep );
  BOOST_CHECK_EQUAL( resolved[2].info.upper_anchor_energy, 187.6 );

  // A "Split by lines" range overlapping a fixed ROI gives one warning, not one per line.
  RoiRange split = fixed_roi( 180.0, 200.0 );
  split.range_limits_type = RoiRange::RangeLimitsType::CanBeBrokenUp;
  vector<string> w2;
  const vector<ResolvedRoi> with_split = resolve( { split, fixed_roi( 184.0, 192.0 ) },
                                                  { Line{ 185.0 }, Line{ 186.0 }, Line{ 190.0 }, Line{ 196.0 } }, w2 );
  check_invariants( with_split );
  BOOST_CHECK_EQUAL( w2.size(), 1 );
}//windows_are_not_merged_through_fixed_roi


BOOST_AUTO_TEST_CASE( nearby_lines_are_always_covered )
{
  // However two lines' windows meet (merged, dividing a shared sideband, or just touching within a
  //  channel), both lines end up in a ROI, and no two ROIs share a channel.
  EdgeDefaults settings;
  settings.lower_edge = settings.upper_edge = RoiEdge{ 1.0E-3, 0.5 };
  for( double second = 202.0; second < 206.0; second += 0.01 )
  {
    vector<string> warnings;
    const vector<ResolvedRoi> resolved = resolve( { line_anchored( 200.0, 200.0 ), line_anchored( second, second ) },
                                                  {}, warnings, settings );
    check_invariants( resolved );
    BOOST_CHECK( (resolved.size() == 1) || (resolved.size() == 2) );
    for( const double line : { 200.0, second } )
    {
      const bool covered = std::any_of( begin(resolved), end(resolved), [line]( const ResolvedRoi &r ){
        return (line > r.roi.lower_energy) && (line < r.roi.upper_energy);
      } );
      BOOST_CHECK_MESSAGE( covered, "Line at " << line << " keV is not in a ROI (other line at "
                           << ((line == second) ? 200.0 : second) << " keV)" );
    }
  }//for( second line energies )
}//nearby_lines_are_always_covered


BOOST_AUTO_TEST_CASE( roi_xml_round_trip )
{
  // A ROI using only version-0 features is still written as version 0, with the legacy elements.
  {
    rapidxml::xml_document<char> doc;
    rapidxml::xml_node<char> *root = doc.allocate_node( rapidxml::node_element, "Root" );
    doc.append_node( root );
    fixed_roi( 180.7, 188.8 ).toXml( root );

    string xml;
    rapidxml::print( std::back_inserter(xml), doc, 0 );
    BOOST_CHECK( xml.find( "version=\"0\"" ) != string::npos );
    BOOST_CHECK( xml.find( "ForceFullRange" ) != string::npos );

    RoiRange back;
    back.fromXml( root->first_node( "RoiRange" ) );
    BOOST_CHECK( back == fixed_roi( 180.7, 188.8 ) );
  }

  // Version 1 features round trip.
  {
    RoiRange roi = line_anchored( 182.52, 185.71, PeakContinuum::OffsetType::FlatStepCDF );
    roi.auto_continuum = true;
    roi.lower_edge.tail_fraction = 2.0E-3;
    roi.upper_edge.sideband_fwhm = 1.75;

    rapidxml::xml_document<char> doc;
    rapidxml::xml_node<char> *root = doc.allocate_node( rapidxml::node_element, "Root" );
    doc.append_node( root );
    roi.toXml( root );

    string xml;
    rapidxml::print( std::back_inserter(xml), doc, 0 );
    BOOST_CHECK( xml.find( "version=\"1\"" ) != string::npos );
    BOOST_CHECK( xml.find( "ForceFullRange" ) == string::npos );

    RoiRange back;
    back.fromXml( root->first_node( "RoiRange" ) );
    BOOST_CHECK( back.range_limits_type == RoiRange::RangeLimitsType::LineAnchored );
    BOOST_CHECK( back.auto_continuum );
    BOOST_CHECK( back.continuum_type == PeakContinuum::OffsetType::FlatStepCDF );
    BOOST_REQUIRE( back.lower_edge.tail_fraction.has_value() );
    BOOST_CHECK_CLOSE( back.lower_edge.tail_fraction.value(), 2.0E-3, 1.0E-4 );
    BOOST_CHECK( !back.lower_edge.sideband_fwhm.has_value() );
    BOOST_REQUIRE( back.upper_edge.sideband_fwhm.has_value() );
    BOOST_CHECK_CLOSE( back.upper_edge.sideband_fwhm.value(), 1.75, 1.0E-4 );
#if( PERFORM_DEVELOPER_CHECKS )
    BOOST_CHECK_NO_THROW( RoiRange::equalEnough( roi, back ) );
#endif
  }

  // A legacy file that only has the old boolean elements still parses.
  {
    string xml = "<RoiRange version=\"0\"><LowerEnergy>100</LowerEnergy><UpperEnergy>200</UpperEnergy>"
                 "<ContinuumType>Linear</ContinuumType><ForceFullRange>false</ForceFullRange>"
                 "<AllowExpandForPeakWidth>true</AllowExpandForPeakWidth></RoiRange>";
    rapidxml::xml_document<char> doc;
    doc.parse<0>( &xml[0] );
    RoiRange back;
    back.fromXml( doc.first_node( "RoiRange" ) );
    BOOST_CHECK( back.range_limits_type == RoiRange::RangeLimitsType::CanExpandForFwhm );
    BOOST_CHECK( !back.auto_continuum );
  }

  BOOST_CHECK( RoiRange::range_limits_type_from_str( "LineAnchored" ) == RoiRange::RangeLimitsType::LineAnchored );
  BOOST_CHECK_EQUAL( string( RoiRange::to_str( RoiRange::RangeLimitsType::LineAnchored ) ), "LineAnchored" );
}//roi_xml_round_trip


BOOST_AUTO_TEST_CASE( skew_from_peak_fit_prefs )
{
  using SkewPrefsUsage = RelActCalcAuto::Options::SkewPrefsUsage;

  // GADRAS-style preferences: a 6-parameter skew, with every parameter given a (fixed) value.
  PeakFitDetPrefs gadras;
  gadras.m_peak_skew_type = PeakDef::SkewType::GadrasGeneric;
  gadras.m_roi_independent_skew = false;
  for( size_t i = 0; i < 6; ++i )
    gadras.m_lower_energy_skew[i] = 0.1*(i + 1);

  RelActCalcAuto::Options preset;
  preset.skew_type = PeakDef::SkewType::GaussExp;

  // Ignore (the default): the preset's skew is untouched.
  RelActCalcAuto::Options ignore = preset;
  ignore.apply_peak_fit_prefs( &gadras );
  BOOST_CHECK( ignore.skew_type == PeakDef::SkewType::GaussExp );
  BOOST_CHECK( !ignore.fixed_lower_skew[0].has_value() );

  // As the preferences specify: the type, and all six values are held fixed.
  RelActCalcAuto::Options as_prefs = preset;
  as_prefs.skew_prefs_usage = SkewPrefsUsage::AsPreferencesSpecify;
  as_prefs.apply_peak_fit_prefs( &gadras );
  BOOST_CHECK( as_prefs.skew_type == PeakDef::SkewType::GadrasGeneric );
  for( size_t i = 0; i < 6; ++i )
  {
    BOOST_REQUIRE( as_prefs.fixed_lower_skew[i].has_value() );
    BOOST_CHECK_CLOSE( as_prefs.fixed_lower_skew[i].value(), 0.1*(i + 1), 1.0E-9 );
    BOOST_CHECK( !as_prefs.start_lower_skew[i].has_value() );
  }

  // Starting values only: the same values, but as starting points, so all parameters are fit.
  RelActCalcAuto::Options start = preset;
  start.skew_prefs_usage = SkewPrefsUsage::StartingValuesOnly;
  start.apply_peak_fit_prefs( &gadras );
  BOOST_CHECK( start.skew_type == PeakDef::SkewType::GadrasGeneric );
  for( size_t i = 0; i < 6; ++i )
  {
    BOOST_CHECK( !start.fixed_lower_skew[i].has_value() );
    BOOST_REQUIRE( start.start_lower_skew[i].has_value() );
  }

  // A parameter without a value in the preferences is fit (not fixed), and an energy-dependent one
  //  keeps its upper-energy value.
  PeakFitDetPrefs gauss_exp;
  gauss_exp.m_peak_skew_type = PeakDef::SkewType::ExpGaussExp;
  gauss_exp.m_roi_independent_skew = false;
  gauss_exp.m_lower_energy_skew[0] = 1.5;
  gauss_exp.m_upper_energy_skew[0] = 2.5;
  RelActCalcAuto::Options partial = preset;
  partial.skew_prefs_usage = SkewPrefsUsage::AsPreferencesSpecify;
  partial.apply_peak_fit_prefs( &gauss_exp );
  BOOST_CHECK( partial.skew_type == PeakDef::SkewType::ExpGaussExp );
  BOOST_REQUIRE( partial.fixed_lower_skew[0].has_value() && partial.fixed_upper_skew[0].has_value() );
  BOOST_CHECK_CLOSE( partial.fixed_upper_skew[0].value(), 2.5, 1.0E-9 );
  BOOST_CHECK( !partial.fixed_lower_skew[1].has_value() );

  // ROI-independent preferences give only the type; preferences without a skew are not used.
  gauss_exp.m_roi_independent_skew = true;
  RelActCalcAuto::Options independent = preset;
  independent.skew_prefs_usage = SkewPrefsUsage::AsPreferencesSpecify;
  independent.apply_peak_fit_prefs( &gauss_exp );
  BOOST_CHECK( independent.skew_type == PeakDef::SkewType::ExpGaussExp );
  BOOST_CHECK( !independent.fixed_lower_skew[0].has_value() );

  PeakFitDetPrefs no_skew;
  RelActCalcAuto::Options unchanged = preset;
  unchanged.skew_prefs_usage = SkewPrefsUsage::AsPreferencesSpecify;
  unchanged.apply_peak_fit_prefs( &no_skew );
  BOOST_CHECK( unchanged.skew_type == PeakDef::SkewType::GaussExp );
  unchanged.apply_peak_fit_prefs( nullptr );
  BOOST_CHECK( unchanged.skew_type == PeakDef::SkewType::GaussExp );

  // The usage and starting values survive the XML round trip.
  rapidxml::xml_document<char> doc;
  rapidxml::xml_node<char> *root = doc.allocate_node( rapidxml::node_element, "Root" );
  doc.append_node( root );
  start.toXml( root );
  RelActCalcAuto::Options back;
  back.fromXml( root->first_node( "Options" ) );
  BOOST_CHECK( back.skew_prefs_usage == SkewPrefsUsage::StartingValuesOnly );
  for( size_t i = 0; i < 6; ++i )
  {
    BOOST_REQUIRE( back.start_lower_skew[i].has_value() );
    BOOST_CHECK_CLOSE( back.start_lower_skew[i].value(), 0.1*(i + 1), 1.0E-4 );
  }
  BOOST_CHECK( RelActCalcAuto::Options::skew_prefs_usage_from_str( "AsPreferencesSpecify" )
               == SkewPrefsUsage::AsPreferencesSpecify );
}//skew_from_peak_fit_prefs
