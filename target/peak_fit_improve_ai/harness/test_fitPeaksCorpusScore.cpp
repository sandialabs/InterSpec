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

#include <cmath>
#include <string>
#include <vector>
#include <memory>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/FitPeaksForNuclides.h"

#include "FitPeaksCorpusScore.h"

#define BOOST_TEST_MODULE FitPeaksCorpusScore_suite
#include <boost/test/included/unit_test.hpp>

using namespace std;
using namespace FitPeaksCorpus;

namespace
{
  /** A flat spectrum of `density` counts per 1 keV channel. */
  shared_ptr<SpecUtils::Measurement> make_flat_spectrum( const size_t nchannel = 2000, const float density = 100.0f )
  {
    auto meas = make_shared<SpecUtils::Measurement>();
    auto counts = make_shared<vector<float>>( nchannel, density );
    meas->set_gamma_counts( counts, 300.0f, 300.0f );
    auto cal = make_shared<SpecUtils::EnergyCalibration>();
    cal->set_polynomial( nchannel, { 0.0f, 1.0f }, {} );
    meas->set_energy_calibration( cal );
    return meas;
  }

  shared_ptr<PeakContinuum> make_continuum( const double lower, const double upper,
                                            const PeakContinuum::OffsetType type, const double density )
  {
    auto cont = make_shared<PeakContinuum>();
    cont->setType( type );
    cont->setRange( lower, upper );
    vector<double> pars( PeakContinuum::num_parameters( type ), 0.0 );
    if( !pars.empty() )
      pars[0] = density;
    cont->setParameters( lower, pars, {} );
    return cont;
  }

  shared_ptr<const PeakDef> make_peak( const double mean, const double sigma, const double amplitude,
                                       const shared_ptr<PeakContinuum> &cont )
  {
    auto peak = make_shared<PeakDef>( mean, sigma, amplitude );
    peak->setContinuum( cont );
    return peak;
  }

  const double sigma = 1.0;              // FWHM = 2.35482 keV on a 1 keV/channel spectrum
  const double fwhm = 2.35482 * sigma;
}//namespace


BOOST_AUTO_TEST_CASE( gaussian_fraction_and_detection_z )
{
  BOOST_CHECK_CLOSE( gaussian_fraction_within_num_fwhm( 1.0 ), 0.9815, 0.05 );
  BOOST_CHECK_CLOSE( gaussian_fraction_within_num_fwhm( 0.5 ), 0.7607, 0.05 );

  const auto data = make_flat_spectrum();
  auto cont = make_continuum( 490.0, 510.0, PeakContinuum::OffsetType::Linear, 100.0 );
  const auto peak = make_peak( 500.0, sigma, 10000.0, cont );
  double b = 0.0;
  const double z = detection_z( *peak, data, { peak }, &b );
  // B = 100 counts/keV over 2 FWHM = 470.96; S = 0.9815*10000
  BOOST_CHECK_CLOSE( b, 100.0 * 2.0 * fwhm, 0.1 );
  BOOST_CHECK_CLOSE( z, 9815.0 / std::sqrt( 9815.0 + 470.96 ), 0.1 );
}


BOOST_AUTO_TEST_CASE( identical_sets_score_zero )
{
  const auto data = make_flat_spectrum();
  auto cont_t = make_continuum( 490.0, 510.0, PeakContinuum::OffsetType::Linear, 100.0 );
  auto cont_f = make_continuum( 490.0, 510.0, PeakContinuum::OffsetType::Linear, 100.0 );
  PeakSet truth = make_peak_set( { make_peak( 500.0, sigma, 5000.0, cont_t ) }, data, 1.5 );
  PeakSet fitted = make_peak_set( { make_peak( 500.0, sigma, 5000.0, cont_f ) }, data, -1.0 );
  const ScoreWeights w;
  const ProblemScore s = score_problem( truth, fitted, w );
  BOOST_CHECK_EQUAL( s.n_truth_scored, 1u );
  BOOST_CHECK_EQUAL( s.n_matched, 1u );
  BOOST_CHECK_SMALL( s.raw_cost(), 1.0e-9 );
  BOOST_CHECK_EQUAL( truth.peaks[0].verdict, "matched" );
  BOOST_CHECK_EQUAL( fitted.peaks[0].verdict, "matched" );
}


BOOST_AUTO_TEST_CASE( missed_peaks_by_class )
{
  const auto data = make_flat_spectrum();
  // strong (z ~ 97), moderate (z ~ 4.5), weak (z ~ 2.2) reference peaks; fitted set is empty
  auto c1 = make_continuum( 190.0, 210.0, PeakContinuum::OffsetType::Linear, 100.0 );
  auto c2 = make_continuum( 390.0, 410.0, PeakContinuum::OffsetType::Linear, 100.0 );
  auto c3 = make_continuum( 590.0, 610.0, PeakContinuum::OffsetType::Linear, 100.0 );
  PeakSet truth = make_peak_set( { make_peak( 200.0, sigma, 10000.0, c1 ),
                                   make_peak( 400.0, sigma, 100.0, c2 ),
                                   make_peak( 600.0, sigma, 50.0, c3 ) }, data, 1.5 );
  PeakSet fitted = make_peak_set( {}, data, -1.0 );
  const ScoreWeights w;
  const ProblemScore s = score_problem( truth, fitted, w );
  BOOST_CHECK_EQUAL( s.missed_strong, 1u );
  BOOST_CHECK_EQUAL( s.missed_moderate, 1u );
  BOOST_CHECK_EQUAL( s.missed_weak, 1u );
  BOOST_CHECK_CLOSE( s.cost_missed, w.missed_strong + w.missed_moderate + w.missed_weak, 1.0e-9 );
  BOOST_CHECK_CLOSE( s.norm_cost(), s.raw_cost() / 3.0, 1.0e-9 );
}


BOOST_AUTO_TEST_CASE( extras_and_ghosts )
{
  const auto data = make_flat_spectrum();
  auto c1 = make_continuum( 190.0, 210.0, PeakContinuum::OffsetType::Linear, 100.0 );
  auto c2 = make_continuum( 390.0, 410.0, PeakContinuum::OffsetType::Linear, 100.0 );
  auto c3 = make_continuum( 590.0, 610.0, PeakContinuum::OffsetType::Linear, 100.0 );
  PeakSet truth = make_peak_set( {}, data, 1.5 );
  PeakSet fitted = make_peak_set( { make_peak( 200.0, sigma, 10000.0, c1 ),   // significant extra
                                    make_peak( 400.0, sigma, 40.0, c2 ),      // weak extra (z ~ 1.8)
                                    make_peak( 600.0, sigma, 2.0, c3 ) },     // ghost
                                  data, -1.0 );
  const ScoreWeights w;
  const ProblemScore s = score_problem( truth, fitted, w );
  BOOST_CHECK_EQUAL( s.extra_significant, 1u );
  BOOST_CHECK_EQUAL( s.extra_weak, 1u );
  BOOST_CHECK_EQUAL( s.extra_ghost, 1u );
  BOOST_CHECK_EQUAL( fitted.peaks[2].verdict, "ghost" );
  BOOST_CHECK_CLOSE( s.cost_extra, w.extra_significant + w.extra_weak + w.extra_ghost, 1.0e-9 );
}


BOOST_AUTO_TEST_CASE( share_separate_disagreement )
{
  const auto data = make_flat_spectrum();
  // Reference: two peaks 3 FWHM apart in separate ROIs; fitted: same peaks sharing one ROI
  const double e1 = 500.0, e2 = 500.0 + 3.0*fwhm;
  auto t1 = make_continuum( e1 - 6.0, e1 + 3.0, PeakContinuum::OffsetType::Linear, 100.0 );
  auto t2 = make_continuum( e2 - 3.0, e2 + 6.0, PeakContinuum::OffsetType::Linear, 100.0 );
  auto f = make_continuum( e1 - 6.0, e2 + 6.0, PeakContinuum::OffsetType::Linear, 100.0 );
  PeakSet truth = make_peak_set( { make_peak( e1, sigma, 3000.0, t1 ), make_peak( e2, sigma, 3000.0, t2 ) }, data, 1.5 );
  PeakSet fitted = make_peak_set( { make_peak( e1, sigma, 3000.0, f ), make_peak( e2, sigma, 3000.0, f ) }, data, -1.0 );
  const ScoreWeights w;
  const ProblemScore s = score_problem( truth, fitted, w );
  BOOST_CHECK_EQUAL( s.n_matched, 2u );
  BOOST_CHECK_EQUAL( s.pairs, 1u );
  BOOST_CHECK_EQUAL( s.share_disagree, 1u );
  BOOST_CHECK_CLOSE( s.cost_share, w.share_disagree, 1.0e-9 );
  BOOST_REQUIRE_EQUAL( s.pair_records.size(), 1u );
  BOOST_CHECK( !s.pair_records[0].truth_share );
  BOOST_CHECK( s.pair_records[0].fit_share );
  BOOST_CHECK_CLOSE( s.pair_records[0].separation_fwhm, 3.0, 1.0e-6 );

  // The reverse (reference shares, fit separates) is also a disagreement
  PeakSet truth2 = make_peak_set( { make_peak( e1, sigma, 3000.0, f ), make_peak( e2, sigma, 3000.0, f ) }, data, 1.5 );
  PeakSet fitted2 = make_peak_set( { make_peak( e1, sigma, 3000.0, t1 ), make_peak( e2, sigma, 3000.0, t2 ) }, data, -1.0 );
  const ProblemScore s2 = score_problem( truth2, fitted2, w );
  BOOST_CHECK_EQUAL( s2.share_disagree, 1u );

  // Far apart peaks do not form a pair
  const double e3 = 900.0;
  auto t3 = make_continuum( e3 - 6.0, e3 + 6.0, PeakContinuum::OffsetType::Linear, 100.0 );
  auto f3 = make_continuum( e3 - 6.0, e3 + 6.0, PeakContinuum::OffsetType::Linear, 100.0 );
  PeakSet truth3 = make_peak_set( { make_peak( e1, sigma, 3000.0, t1 ), make_peak( e3, sigma, 3000.0, t3 ) }, data, 1.5 );
  PeakSet fitted3 = make_peak_set( { make_peak( e1, sigma, 3000.0, f ), make_peak( e3, sigma, 3000.0, f3 ) }, data, -1.0 );
  const ProblemScore s3 = score_problem( truth3, fitted3, w );
  BOOST_CHECK_EQUAL( s3.pairs, 0u );
}


BOOST_AUTO_TEST_CASE( continuum_family_disagreement )
{
  const auto data = make_flat_spectrum();
  const ScoreWeights w;

  auto t_step = make_continuum( 490.0, 510.0, PeakContinuum::OffsetType::FlatStepCDF, 100.0 );
  auto f_lin = make_continuum( 490.0, 510.0, PeakContinuum::OffsetType::Linear, 100.0 );
  PeakSet truth = make_peak_set( { make_peak( 500.0, sigma, 20000.0, t_step ) }, data, 1.5 );
  PeakSet fitted = make_peak_set( { make_peak( 500.0, sigma, 20000.0, f_lin ) }, data, -1.0 );
  const ProblemScore s = score_problem( truth, fitted, w );
  BOOST_CHECK_EQUAL( s.family_disagree, 1u );
  BOOST_CHECK_CLOSE( s.cost_family, w.family_step_disagree, 1.0e-9 );

  auto t_quad = make_continuum( 490.0, 510.0, PeakContinuum::OffsetType::Quadratic, 100.0 );
  auto f_lin2 = make_continuum( 490.0, 510.0, PeakContinuum::OffsetType::Linear, 100.0 );
  PeakSet truth2 = make_peak_set( { make_peak( 500.0, sigma, 20000.0, t_quad ) }, data, 1.5 );
  PeakSet fitted2 = make_peak_set( { make_peak( 500.0, sigma, 20000.0, f_lin2 ) }, data, -1.0 );
  const ProblemScore s2 = score_problem( truth2, fitted2, w );
  BOOST_CHECK_CLOSE( s2.cost_family, w.family_poly_disagree, 1.0e-9 );

  // FlatStep (non-CDF) and FlatStepCDF are the same family
  auto t_step_noncdf = make_continuum( 490.0, 510.0, PeakContinuum::OffsetType::FlatStep, 100.0 );
  auto f_step_cdf = make_continuum( 490.0, 510.0, PeakContinuum::OffsetType::FlatStepCDF, 100.0 );
  PeakSet truth3 = make_peak_set( { make_peak( 500.0, sigma, 20000.0, t_step_noncdf ) }, data, 1.5 );
  PeakSet fitted3 = make_peak_set( { make_peak( 500.0, sigma, 20000.0, f_step_cdf ) }, data, -1.0 );
  const ProblemScore s3 = score_problem( truth3, fitted3, w );
  BOOST_CHECK_EQUAL( s3.family_disagree, 0u );
}


BOOST_AUTO_TEST_CASE( roi_extent_cost )
{
  const auto data = make_flat_spectrum();
  const ScoreWeights w;
  // Fitted ROI is 2 FWHM wider on each side than the reference ROI
  auto t = make_continuum( 500.0 - 3.0*fwhm, 500.0 + 3.0*fwhm, PeakContinuum::OffsetType::Linear, 100.0 );
  auto f = make_continuum( 500.0 - 5.0*fwhm, 500.0 + 5.0*fwhm, PeakContinuum::OffsetType::Linear, 100.0 );
  PeakSet truth = make_peak_set( { make_peak( 500.0, sigma, 5000.0, t ) }, data, 1.5 );
  PeakSet fitted = make_peak_set( { make_peak( 500.0, sigma, 5000.0, f ) }, data, -1.0 );
  const ProblemScore s = score_problem( truth, fitted, w );
  BOOST_CHECK_EQUAL( s.extent_sides, 2u );
  BOOST_CHECK_EQUAL( s.extent_gt1fwhm_sides, 2u );
  const double expected = 2.0 * w.extent_per_fwhm * (2.0 - w.extent_deadband_fwhm);
  BOOST_CHECK_CLOSE( s.cost_extent, expected, 1.0e-6 );
  BOOST_REQUIRE_EQUAL( s.roi_records.size(), 1u );
  BOOST_CHECK_CLOSE( s.roi_records[0].d_lower_fwhm, -2.0, 1.0e-6 );
  BOOST_CHECK_CLOSE( s.roi_records[0].d_upper_fwhm, 2.0, 1.0e-6 );
}


BOOST_AUTO_TEST_CASE( dont_care_reference_peaks )
{
  const auto data = make_flat_spectrum();
  const ScoreWeights w;
  auto t = make_continuum( 490.0, 510.0, PeakContinuum::OffsetType::Linear, 100.0 );
  auto f = make_continuum( 490.0, 510.0, PeakContinuum::OffsetType::Linear, 100.0 );
  // amplitude 3 on 470 continuum counts: z ~ 0.14 -> don't care
  PeakSet truth = make_peak_set( { make_peak( 500.0, sigma, 3.0, t ) }, data, 1.5 );
  BOOST_CHECK( truth.peaks[0].dont_care );

  PeakSet fitted_empty = make_peak_set( {}, data, -1.0 );
  const ProblemScore s = score_problem( truth, fitted_empty, w );
  BOOST_CHECK_EQUAL( s.n_truth_scored, 0u );
  BOOST_CHECK_EQUAL( s.n_truth_dontcare, 1u );
  BOOST_CHECK_SMALL( s.raw_cost(), 1.0e-9 );
  BOOST_CHECK_EQUAL( truth.peaks[0].verdict, "dontcare" );

  PeakSet fitted = make_peak_set( { make_peak( 500.0, sigma, 3.0, f ) }, data, -1.0 );
  const ProblemScore s2 = score_problem( truth, fitted, w );
  BOOST_CHECK_EQUAL( s2.neutral, 1u );
  BOOST_CHECK_EQUAL( fitted.peaks[0].verdict, "neutral" );
  BOOST_CHECK_SMALL( s2.raw_cost(), 1.0e-9 );
}


BOOST_AUTO_TEST_CASE( matching_window_and_parameter_costs )
{
  const auto data = make_flat_spectrum();
  const ScoreWeights w;
  auto t = make_continuum( 490.0, 510.0, PeakContinuum::OffsetType::Linear, 100.0 );
  auto f = make_continuum( 490.0, 510.0, PeakContinuum::OffsetType::Linear, 100.0 );
  PeakSet truth = make_peak_set( { make_peak( 500.0, sigma, 5000.0, t ) }, data, 1.5 );

  // Offset of 0.5 FWHM matches, with a mean-offset cost; area 20% off gives a capped pull cost
  PeakSet fitted = make_peak_set( { make_peak( 500.0 + 0.5*fwhm, sigma, 6000.0, f ) }, data, -1.0 );
  const ProblemScore s = score_problem( truth, fitted, w );
  BOOST_CHECK_EQUAL( s.n_matched, 1u );
  BOOST_CHECK_CLOSE( s.cost_mean, w.mean_offset_weight * (0.5 - w.mean_offset_deadband_fwhm), 1.0e-6 );
  const double sigma_area = std::sqrt( 5000.0 + std::pow( w.area_rel_floor*5000.0, 2.0 ) );
  const double pull = 1000.0 / sigma_area;
  BOOST_CHECK_CLOSE( s.cost_area, w.area_pull_weight * std::min( pull, w.area_pull_cap ), 1.0e-6 );

  // Offset of 1.5 FWHM does not match
  auto f2 = make_continuum( 490.0, 510.0, PeakContinuum::OffsetType::Linear, 100.0 );
  PeakSet fitted2 = make_peak_set( { make_peak( 500.0 + 1.5*fwhm, sigma, 5000.0, f2 ) }, data, -1.0 );
  const ProblemScore s2 = score_problem( truth, fitted2, w );
  BOOST_CHECK_EQUAL( s2.n_matched, 0u );
  BOOST_CHECK_EQUAL( s2.missed_strong, 1u );
  BOOST_CHECK_EQUAL( s2.extra_significant, 1u );
}


BOOST_AUTO_TEST_CASE( config_reflection_roundtrip )
{
  using FitPeaksForNuclides::PeakFitForNuclideConfig;
  PeakFitForNuclideConfig config = PeakFitForNuclideConfig::default_config( PeakFitUtils::CoarseResolutionType::High );
  const vector<string> names = PeakFitForNuclideConfig::field_names();
  BOOST_CHECK_GT( names.size(), 30u );

  for( const string &name : names )
  {
    const string value = config.get_field( name );
    BOOST_CHECK_MESSAGE( config.set_field( name, value ), "round-trip of " << name << "=" << value );
    BOOST_CHECK_EQUAL( config.get_field( name ), value );
  }

  BOOST_CHECK( config.set_field( "auto_keep_significance_z", "3.25" ) );
  BOOST_CHECK_CLOSE( config.auto_keep_significance_z, 3.25, 1.0e-12 );
  BOOST_CHECK( config.set_field( "fit_energy_cal", "false" ) );
  BOOST_CHECK( !config.fit_energy_cal );
  BOOST_CHECK( config.set_field( "rel_eff_eqn_order", "3" ) );
  BOOST_CHECK_EQUAL( config.rel_eff_eqn_order, 3u );
  BOOST_CHECK( config.set_field( "skew_type", "NoSkew" ) );
  BOOST_CHECK( config.skew_type == PeakDef::SkewType::NoSkew );
  BOOST_CHECK( !config.set_field( "no_such_field", "1" ) );
  BOOST_CHECK( !config.set_field( "auto_keep_significance_z", "abc" ) );
  BOOST_CHECK_THROW( config.get_field( "no_such_field" ), std::invalid_argument );
  BOOST_CHECK( config.to_string( ";" ).find( "auto_keep_significance_z=3.25" ) != string::npos );
}
