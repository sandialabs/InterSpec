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
 */

// Must be defined before Windows.h (or any header that includes it) is included
#ifdef _WIN32
  #define WIN32_LEAN_AND_MEAN
  #include <winsock2.h>
  #include <windows.h>
#endif

#include "InterSpec_config.h"

#include <cmath>
#include <deque>
#include <memory>
#include <string>
#include <vector>

#define BOOST_TEST_MODULE BatchPeakMda_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/BatchPeak.h"
#include "InterSpec/BatchInfoLog.h"
#include "InterSpec/DetectionLimitCalc.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace std;
using namespace boost::unit_test;


namespace
{
/** Sets the static data directory, so `sandia.decay.xml` and friends can be found. */
void set_data_dir()
{
  static bool s_have_set = false;
  if( s_have_set )
    return;
  s_have_set = true;

  const int argc = boost::unit_test::framework::master_test_suite().argc;
  char ** const argv = boost::unit_test::framework::master_test_suite().argv;

  string datadir;
  for( int i = 1; i < argc; ++i )
  {
    const string arg = argv[i];
    if( SpecUtils::istarts_with( arg, "--datadir=" ) )
      datadir = arg.substr( 10 );
  }

  SpecUtils::ireplace_all( datadir, "%20", " " );

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
  }//if( datadir.empty() )

  const string decay_file = SpecUtils::append_path( datadir, "sandia.decay.xml" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file(decay_file), "sandia.decay.xml not at '" << decay_file << "'" );
  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );
}//void set_data_dir()


const size_t sm_num_channels = 1024;
const float sm_channel_width = 2.0f;    //keV per channel
const float sm_live_time = 300.0f;      //seconds
const float sm_continuum_cps = 0.5f;    //counts per channel, per second
const double sm_fwhm = 8.0;             //keV; ~3 keV sigma

/** Returns the FWHM the synthetic spectrum peaks are made with. */
double peak_fwhm(){ return sm_fwhm; }


/** Builds a spectrum with a flat continuum, and Gaussian peaks of the requested areas. */
shared_ptr<SpecMeas> make_spec_meas( const vector<pair<double,double>> &energy_and_area )
{
  auto counts = make_shared<vector<float>>( sm_num_channels, sm_continuum_cps*sm_live_time );

  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( sm_num_channels, { 0.0f, sm_channel_width }, {} );

  const double sigma = sm_fwhm / 2.35482;

  for( const pair<double,double> &ea : energy_and_area )
  {
    for( size_t i = 0; i < sm_num_channels; ++i )
    {
      const double lower = cal->energy_for_channel( i );
      const double upper = cal->energy_for_channel( i + 1 );
      const double frac = 0.5*( erf( (upper - ea.first)/(sqrt(2.0)*sigma) )
                                - erf( (lower - ea.first)/(sqrt(2.0)*sigma) ) );
      (*counts)[i] += static_cast<float>( ea.second * frac );
    }
  }//for( const pair<double,double> &ea : energy_and_area )

  auto m = make_shared<SpecUtils::Measurement>();
  m->set_gamma_counts( counts, sm_live_time, sm_live_time );
  m->set_energy_calibration( cal );
  m->set_sample_number( 1 );
  m->set_detector_name( "Det1" );
  m->set_source_type( SpecUtils::SourceType::Foreground );

  auto answer = make_shared<SpecMeas>();
  answer->add_measurement( m, false );
  answer->cleanup_after_load();

  return answer;
}//shared_ptr<SpecMeas> make_spec_meas(...)


/** Returns a spectrum file with peaks at the requested energies set into it, so it can be used
 as an exemplar.  The peak amplitudes only matter for being non-zero.
 */
shared_ptr<SpecMeas> make_exemplar( const vector<double> &energies )
{
  vector<pair<double,double>> energy_and_area;
  for( const double energy : energies )
    energy_and_area.push_back( { energy, 5000.0 } );

  shared_ptr<SpecMeas> answer = make_spec_meas( energy_and_area );

  const double sigma = sm_fwhm / 2.35482;

  deque<shared_ptr<const PeakDef>> peaks;
  for( const double energy : energies )
  {
    auto peak = make_shared<PeakDef>( energy, sigma, 5000.0 );
    peak->continuum()->setType( PeakContinuum::OffsetType::Linear );
    peak->continuum()->setRange( energy - 3.0*sigma, energy + 3.0*sigma );
    peaks.push_back( peak );
  }

  answer->setPeaks( peaks, answer->sample_numbers() );

  return answer;
}//shared_ptr<SpecMeas> make_exemplar(...)


/** Splits one line of CSV, honoring double-quoted fields; enough to check field counts. */
vector<string> split_csv_line( const string &line )
{
  vector<string> fields;
  string current;
  bool in_quote = false;

  for( size_t i = 0; i < line.size(); ++i )
  {
    const char c = line[i];
    // RFC 4180: a quote only opens a quoted field when it immediately follows the delimiter.
    //  A quote preceded by anything else (e.g. a space) is a literal character, and the field is
    //  NOT quoted - which is exactly the mistake this helper exists to catch.
    if( (c == '"') && (in_quote || current.empty()) )
      in_quote = !in_quote;
    else if( (c == ',') && !in_quote )
    {
      fields.push_back( current );
      current.clear();
    }else
    {
      current += c;
    }
  }//for( loop over characters )

  fields.push_back( current );

  return fields;
}//vector<string> split_csv_line( const string & )


BatchPeak::BatchPeakFitOptions default_options()
{
  BatchPeak::BatchPeakFitOptions options;
  options.fit_all_peaks = false;
  options.to_stdout = false;
  options.refit_energy_cal = false;
  options.use_exemplar_energy_cal = false;
  options.write_n42_with_results = false;
  options.show_nonfit_peaks = false;
  options.overwrite_output_files = false;
  options.create_csv_output = false;
  options.create_json_output = false;
  options.concatenate_to_n42 = false;
  options.use_existing_background_peaks = false;
  options.use_exemplar_energy_cal_for_background = false;
  options.peak_stat_threshold = 2.0;
  options.peak_hypothesis_threshold = 1.0;
  options.template_include_dir = "none";

  return options;
}//BatchPeak::BatchPeakFitOptions default_options()


/** Runs a peak fit of `exemplar_energies` against a spectrum that only contains peaks at
 `present_energies`.
 */
BatchPeak::BatchPeakFitResult run_fit( const vector<double> &exemplar_energies,
                                       const vector<pair<double,double>> &present_energies,
                                       const BatchPeak::BatchPeakFitOptions &options )
{
  const shared_ptr<SpecMeas> exemplar = make_exemplar( exemplar_energies );
  const shared_ptr<SpecMeas> measurement = make_spec_meas( present_energies );

  return BatchPeak::fit_peaks_in_file( "exemplar.n42", {}, exemplar, "input.n42", measurement,
                                       {}, options );
}//BatchPeak::BatchPeakFitResult run_fit(...)
}//namespace


BOOST_AUTO_TEST_CASE( UnfitExemplarPeaksAreReported )
{
  set_data_dir();

  // The exemplar has peaks at 661 and 1173 keV, but only the 661 keV peak is in the spectrum.
  //  Before the fix to the shadowed `unused_exemplar_peaks`, `unfit_exemplar_peaks` was always
  //  empty on a successful fit, so `NotFitPeaks` never made it into any report.
  BatchPeak::BatchPeakFitOptions options = default_options();
  options.not_fit_peak_mda = BatchPeak::NotFitPeakMdaMethod::None;

  const BatchPeak::BatchPeakFitResult result = run_fit( {661.0, 1173.0}, {{661.0, 5000.0}}, options );

  BOOST_REQUIRE( result.success );
  BOOST_CHECK_EQUAL( result.exemplar_peaks.size(), 2 );

  BOOST_REQUIRE_EQUAL( result.fit_peaks.size(), 1 );
  BOOST_CHECK_CLOSE( result.fit_peaks.front()->mean(), 661.0, 1.0 );

  BOOST_REQUIRE_EQUAL( result.unfit_exemplar_peaks.size(), 1 );
  BOOST_CHECK_CLOSE( result.unfit_exemplar_peaks.front()->mean(), 1173.0, 1.0 );

  // Limits were turned off for this case
  BOOST_CHECK( result.not_fit_peak_mdas.empty() );
}//BOOST_AUTO_TEST_CASE( UnfitExemplarPeaksAreReported )


BOOST_AUTO_TEST_CASE( AllPeaksFitLeavesNoLimits )
{
  set_data_dir();

  BatchPeak::BatchPeakFitOptions options = default_options();

  const BatchPeak::BatchPeakFitResult result
        = run_fit( {661.0, 1173.0}, {{661.0, 5000.0}, {1173.0, 5000.0}}, options );

  BOOST_REQUIRE( result.success );
  BOOST_CHECK_EQUAL( result.fit_peaks.size(), 2 );
  BOOST_CHECK( result.unfit_exemplar_peaks.empty() );
  BOOST_CHECK( result.not_fit_peak_mdas.empty() );
}//BOOST_AUTO_TEST_CASE( AllPeaksFitLeavesNoLimits )


BOOST_AUTO_TEST_CASE( CurrieLimitMatchesDirectCalculation )
{
  set_data_dir();

  BatchPeak::BatchPeakFitOptions options = default_options();
  options.not_fit_peak_mda = BatchPeak::NotFitPeakMdaMethod::Currie;
  options.mda_confidence_level = 0.95;
  options.mda_num_side_channels = 4;

  const BatchPeak::BatchPeakFitResult result = run_fit( {661.0, 1173.0}, {{661.0, 5000.0}}, options );

  BOOST_REQUIRE( result.success );
  BOOST_REQUIRE_EQUAL( result.unfit_exemplar_peaks.size(), 1 );
  BOOST_REQUIRE_EQUAL( result.not_fit_peak_mdas.size(), 1 );

  const BatchPeak::NotFitPeakMda &mda = result.not_fit_peak_mdas.front();
  BOOST_CHECK( mda.exemplar_peak == result.unfit_exemplar_peaks.front() );
  BOOST_REQUIRE( mda.currie_computed );
  BOOST_CHECK( mda.result_type == BatchPeak::NotFitPeakMda::MdaResultType::NotDetected );

  // There should be no activity information for a plain peak fit
  BOOST_CHECK( !mda.has_activity );
  BOOST_CHECK( !mda.decon_attempted );

  // The limit should match a calculation done directly, with the same inputs
  DetectionLimitCalc::CurrieMdaInput input;
  input.spectrum = result.spectrum;
  input.gamma_energy = 1173.0f;
  input.roi_lower_energy = static_cast<float>( 1173.0 - 1.25*peak_fwhm() );
  input.roi_upper_energy = static_cast<float>( 1173.0 + 1.25*peak_fwhm() );
  input.num_lower_side_channels = 4;
  input.num_upper_side_channels = 4;
  input.detection_probability = 0.95;
  input.additional_uncertainty = 0.0f;

  const DetectionLimitCalc::CurrieMdaResult expected = DetectionLimitCalc::currie_mda_calc( input );

  BOOST_CHECK_CLOSE( mda.currie_result.decision_threshold, expected.decision_threshold, 0.01 );
  BOOST_CHECK_CLOSE( mda.currie_result.detection_limit, expected.detection_limit, 0.01 );
  BOOST_CHECK_CLOSE( mda.currie_result.upper_limit, expected.upper_limit, 0.01 );

  // A flat continuum of 150 counts/channel over ~10 channels; the limit should be a few hundred
  //  counts - a sanity check that we are not off by orders of magnitude.
  BOOST_CHECK( mda.currie_result.detection_limit > 10.0f );
  BOOST_CHECK( mda.currie_result.detection_limit < 1000.0f );

  // The description should be usable directly in a report
  BOOST_CHECK( mda.description.find("Not detected") != string::npos );
  BOOST_CHECK( mda.description.find("Lc") != string::npos );
  BOOST_CHECK( mda.description.find("95%") != string::npos );

  // ...and there should be a brief phrase, short enough for a table cell
  BOOST_CHECK_EQUAL( mda.short_description, string("Less than Lc") );
}//BOOST_AUTO_TEST_CASE( CurrieLimitMatchesDirectCalculation )


BOOST_AUTO_TEST_CASE( ConfidenceLevelAndSideChannelsAreUsed )
{
  set_data_dir();

  BatchPeak::BatchPeakFitOptions options = default_options();
  options.not_fit_peak_mda = BatchPeak::NotFitPeakMdaMethod::Currie;
  options.mda_confidence_level = 0.99;
  options.mda_num_side_channels = 8;

  const BatchPeak::BatchPeakFitResult result = run_fit( {661.0, 1173.0}, {{661.0, 5000.0}}, options );

  BOOST_REQUIRE_EQUAL( result.not_fit_peak_mdas.size(), 1 );
  const BatchPeak::NotFitPeakMda &mda = result.not_fit_peak_mdas.front();

  BOOST_REQUIRE( mda.currie_computed );
  BOOST_CHECK_CLOSE( mda.currie_result.input.detection_probability, 0.99, 1.0E-6 );
  BOOST_CHECK_EQUAL( mda.currie_result.input.num_lower_side_channels, 8 );
  BOOST_CHECK_EQUAL( mda.currie_result.input.num_upper_side_channels, 8 );
  BOOST_CHECK( mda.description.find("99%") != string::npos );

  // A higher confidence level must give a larger limit
  BatchPeak::BatchPeakFitOptions ninety_five = options;
  ninety_five.mda_confidence_level = 0.95;
  const BatchPeak::BatchPeakFitResult result_95
        = run_fit( {661.0, 1173.0}, {{661.0, 5000.0}}, ninety_five );

  BOOST_REQUIRE_EQUAL( result_95.not_fit_peak_mdas.size(), 1 );
  BOOST_CHECK( result_95.not_fit_peak_mdas.front().currie_result.detection_limit
               < mda.currie_result.detection_limit );

  BOOST_CHECK( BatchPeak::confidence_level_str(0.95) == "95%" );
  BOOST_CHECK( BatchPeak::confidence_level_str(0.99) == "99%" );
}//BOOST_AUTO_TEST_CASE( ConfidenceLevelAndSideChannelsAreUsed )


BOOST_AUTO_TEST_CASE( RoiWidthOption )
{
  set_data_dir();

  // The fraction of a Gaussian inside a region of the given width, in FWHM
  BOOST_CHECK_CLOSE( BatchPeak::gaussian_fraction_in_roi(2.5), 0.996755, 0.01 );
  BOOST_CHECK_CLOSE( BatchPeak::gaussian_fraction_in_roi(2.0), 0.981468, 0.01 );
  BOOST_CHECK_CLOSE( BatchPeak::gaussian_fraction_in_roi(1.0), 0.760968, 0.01 );
  BOOST_CHECK( BatchPeak::gaussian_fraction_in_roi(0.0) == 0.0 );

  BatchPeak::BatchPeakFitOptions wide = default_options();
  wide.not_fit_peak_mda = BatchPeak::NotFitPeakMdaMethod::Currie;
  wide.mda_roi_num_fwhm = 2.5;  // the default

  BatchPeak::BatchPeakFitOptions narrow = wide;
  narrow.mda_roi_num_fwhm = 1.0;

  const BatchPeak::BatchPeakFitResult wide_result
        = run_fit( {661.0, 1173.0}, {{661.0, 5000.0}}, wide );
  const BatchPeak::BatchPeakFitResult narrow_result
        = run_fit( {661.0, 1173.0}, {{661.0, 5000.0}}, narrow );

  BOOST_REQUIRE_EQUAL( wide_result.not_fit_peak_mdas.size(), 1 );
  BOOST_REQUIRE_EQUAL( narrow_result.not_fit_peak_mdas.size(), 1 );

  const BatchPeak::NotFitPeakMda &w = wide_result.not_fit_peak_mdas.front();
  const BatchPeak::NotFitPeakMda &n = narrow_result.not_fit_peak_mdas.front();
  BOOST_REQUIRE( w.currie_computed && n.currie_computed );

  // The peak region should be the requested number of FWHM wide, centered on the peak
  const double fwhm = peak_fwhm();
  BOOST_CHECK_CLOSE( w.currie_result.input.roi_upper_energy - w.currie_result.input.roi_lower_energy,
                    2.5*fwhm, 0.1 );
  BOOST_CHECK_CLOSE( n.currie_result.input.roi_upper_energy - n.currie_result.input.roi_lower_energy,
                    1.0*fwhm, 0.1 );

  // Less continuum in the region means a smaller limit...
  BOOST_CHECK( n.currie_result.detection_limit < w.currie_result.detection_limit );

  // ...but the narrow region only holds part of the peak, which must be reported
  BOOST_CHECK_CLOSE( w.signal_fraction_in_roi, 0.996755, 0.01 );
  BOOST_CHECK_CLOSE( n.signal_fraction_in_roi, 0.760968, 0.01 );
  BOOST_CHECK_MESSAGE( n.description.find("holds only") != string::npos,
                      "Narrow-region caveat missing from: " << n.description );
  BOOST_CHECK_MESSAGE( w.description.find("holds only") == string::npos,
                      "Caveat should not appear at the default width: " << w.description );

  // Caveats belong at the end of the paragraph, after the result and any activity information
  BOOST_CHECK( n.description == (n.result_summary + "  " + n.caveats) );
  BOOST_CHECK( n.short_description == "Less than Lc" );

  // The width used, and the fraction of the peak it holds, must reach templates
  nlohmann::json data;
  BatchInfoLog::add_peak_fit_results_to_json( data, narrow_result );
  BatchInfoLog::add_peak_fit_options_to_json( data, narrow );

#if( ALLOW_SPECIFY_MDA_ROI_WIDTH )
  BOOST_CHECK_CLOSE( data["PeakFitOptions"]["MdaRoiNumFwhm"].get<double>(), 1.0, 1.0E-6 );
#else
  // The width isnt a user option, so it isnt reported as one
  BOOST_CHECK( !data["PeakFitOptions"].contains("MdaRoiNumFwhm") );
#endif

  const nlohmann::json &mda_json = data["NotFitPeaks"]["Peaks"][0]["Mda"];
  BOOST_CHECK_CLOSE( mda_json["SignalFractionInRoi"].get<double>(), 0.760968, 0.01 );
  BOOST_CHECK_CLOSE( mda_json["RoiWidth_keV"].get<double>(), fwhm, 0.1 );
}//BOOST_AUTO_TEST_CASE( RoiWidthOption )


BOOST_AUTO_TEST_CASE( PeakOffSpectrumGivesErrorResult )
{
  set_data_dir();

  // The spectrum starts at 0 keV, so a peak at 2 keV has no room for the lower side channels
  BatchPeak::BatchPeakFitOptions options = default_options();
  options.not_fit_peak_mda = BatchPeak::NotFitPeakMdaMethod::Currie;

  const BatchPeak::BatchPeakFitResult result = run_fit( {2.0, 661.0}, {{661.0, 5000.0}}, options );

  BOOST_REQUIRE( result.success );
  BOOST_REQUIRE_EQUAL( result.not_fit_peak_mdas.size(), 1 );

  const BatchPeak::NotFitPeakMda &mda = result.not_fit_peak_mdas.front();
  BOOST_CHECK( !mda.currie_computed );
  BOOST_CHECK( mda.result_type == BatchPeak::NotFitPeakMda::MdaResultType::Error );
  BOOST_CHECK( !mda.currie_error.empty() );
  BOOST_CHECK( mda.description.find("could not be computed") != string::npos );
  BOOST_CHECK_EQUAL( mda.short_description, string("Not computed") );
}//BOOST_AUTO_TEST_CASE( PeakOffSpectrumGivesErrorResult )


BOOST_AUTO_TEST_CASE( DeconLimitInCounts )
{
  set_data_dir();

  BatchPeak::BatchPeakFitOptions options = default_options();
  options.not_fit_peak_mda = BatchPeak::NotFitPeakMdaMethod::CurrieAndDecon;

  const BatchPeak::BatchPeakFitResult result = run_fit( {661.0, 1173.0}, {{661.0, 5000.0}}, options );

  BOOST_REQUIRE_EQUAL( result.not_fit_peak_mdas.size(), 1 );
  const BatchPeak::NotFitPeakMda &mda = result.not_fit_peak_mdas.front();

  BOOST_CHECK( mda.decon_attempted );
  BOOST_CHECK( mda.decon_quantity_is_counts );
  BOOST_REQUIRE_MESSAGE( mda.decon_computed, "Deconvolution limit failed: " << mda.decon_error );
  BOOST_REQUIRE( mda.decon_result );
  BOOST_REQUIRE( mda.decon_result->foundUpperCl );

  // The two methods should agree to within a factor of a few
  const double currie_limit = mda.currie_result.upper_limit;
  const double decon_limit = mda.decon_result->upperLimit;
  BOOST_CHECK_MESSAGE( (decon_limit > 0.2*currie_limit) && (decon_limit < 5.0*currie_limit),
                      "Deconvolution limit of " << decon_limit << " counts is not comparable to the"
                      " Currie limit of " << currie_limit << " counts." );

  // This spectrum is an isolated Gaussian on a flat continuum - the case both methods model
  //  exactly - so the "these two disagree" diagnostic must not fire.  If it does here, the
  //  threshold is too tight to be useful on real data.
  BOOST_CHECK_CLOSE( mda.decon_over_currie_ratio, decon_limit/currie_limit, 0.01 );
  BOOST_CHECK_MESSAGE( !mda.methods_disagree,
                      "Methods flagged as disagreeing on a clean, isolated peak: decon/currie = "
                      << mda.decon_over_currie_ratio << " (" << mda.caveats << ")" );

  // Pin the absolute value.  The spectrum is deterministic and the peak width is exactly known,
  //  so this catches a change to the FWHM<->sigma conversion inside `decon_compute_peaks`
  //  (`sigma = fwhm/2.35482`); the loose ratio band above does not - reverting that constant to
  //  the old, wrong 2.634 moves this number by ~11% but keeps the ratio inside the band.
  BOOST_CHECK_CLOSE( decon_limit, 101.6, 3.0 );
}//BOOST_AUTO_TEST_CASE( DeconLimitInCounts )


BOOST_AUTO_TEST_CASE( DetectedButNotFitPeak )
{
  set_data_dir();

  // A peak that is present in the data, but which we make the fit reject by demanding an
  //  unreachable improvement to the chi2.
  BatchPeak::BatchPeakFitOptions options = default_options();
  options.not_fit_peak_mda = BatchPeak::NotFitPeakMdaMethod::Currie;
  options.peak_stat_threshold = 1000.0;
  options.peak_hypothesis_threshold = 1000.0;

  const BatchPeak::BatchPeakFitResult result
        = run_fit( {1173.0}, {{1173.0, 5000.0}}, options );

  BOOST_REQUIRE( result.success );
  BOOST_REQUIRE_EQUAL( result.unfit_exemplar_peaks.size(), 1 );
  BOOST_REQUIRE_EQUAL( result.not_fit_peak_mdas.size(), 1 );

  const BatchPeak::NotFitPeakMda &mda = result.not_fit_peak_mdas.front();
  BOOST_REQUIRE( mda.currie_computed );
  BOOST_CHECK( mda.result_type == BatchPeak::NotFitPeakMda::MdaResultType::Detected );
  BOOST_CHECK( mda.description.find("Signal present but peak not fit") != string::npos );
  BOOST_CHECK_EQUAL( mda.short_description, string("Greater than Lc") );
}//BOOST_AUTO_TEST_CASE( DetectedButNotFitPeak )


BOOST_AUTO_TEST_CASE( JsonHasMdaOnlyForNotFitPeaks )
{
  set_data_dir();

  BatchPeak::BatchPeakFitOptions options = default_options();
  options.not_fit_peak_mda = BatchPeak::NotFitPeakMdaMethod::Currie;

  const BatchPeak::BatchPeakFitResult result = run_fit( {661.0, 1173.0}, {{661.0, 5000.0}}, options );
  BOOST_REQUIRE_EQUAL( result.not_fit_peak_mdas.size(), 1 );

  nlohmann::json data;
  BOOST_REQUIRE_NO_THROW( BatchInfoLog::add_peak_fit_results_to_json( data, result ) );
  BOOST_REQUIRE_NO_THROW( BatchInfoLog::add_peak_fit_options_to_json( data, options ) );

  BOOST_REQUIRE( data.contains("NotFitPeaks") );
  BOOST_CHECK_EQUAL( data["FitAllPeaks"].get<bool>(), false );
  BOOST_CHECK_EQUAL( data["AnyNotFitPeakMda"].get<bool>(), true );
  BOOST_CHECK_EQUAL( data["NotFitPeaks"]["HasMdas"].get<bool>(), true );

  const nlohmann::json &peak = data["NotFitPeaks"]["Peaks"][0];
  BOOST_REQUIRE( peak.contains("Mda") );
  BOOST_CHECK_EQUAL( peak["HasMda"].get<bool>(), true );

  const nlohmann::json &mda = peak["Mda"];
  BOOST_CHECK_EQUAL( mda["ResultType"].get<string>(), string("NotDetected") );
  BOOST_CHECK_EQUAL( mda["ShortDescription"].get<string>(), string("Less than Lc") );
  BOOST_CHECK_EQUAL( mda["HasCaveats"].get<bool>(), false );
  BOOST_CHECK_EQUAL( mda["CurrieComputed"].get<bool>(), true );
  BOOST_CHECK_EQUAL( mda["HasMdaActivity"].get<bool>(), false );
  BOOST_CHECK_EQUAL( mda["DeconAttempted"].get<bool>(), false );
  BOOST_CHECK_EQUAL( mda["ConfidenceLevelStr"].get<string>(), string("95%") );
  BOOST_CHECK_EQUAL( mda["NumSideChannels"].get<int>(), 4 );
  BOOST_CHECK( mda["DetectionLimit_counts"].get<double>() > 0.0 );
  BOOST_CHECK( !mda["Description"].get<string>().empty() );

  // Peaks that were fit must not gain any new entries, so that reports written before detection
  //  limits existed see no change at all.
  BOOST_REQUIRE( data.contains("FitPeaks") );
  BOOST_CHECK( !data["FitPeaks"].contains("HasMdas") );
  BOOST_CHECK( !data["FitPeaks"]["Peaks"][0].contains("HasMda") );
  BOOST_CHECK( !data["FitPeaks"]["Peaks"][0].contains("Mda") );
  BOOST_CHECK( !data["ExemplarPeaks"]["Peaks"][0].contains("HasMda") );

  // The options used must be available to templates
  const nlohmann::json &opts = data["PeakFitOptions"];
  BOOST_CHECK_EQUAL( opts["NotFitPeakMdaMethod"].get<string>(), string("currie") );
  BOOST_CHECK_EQUAL( opts["MdaNumSideChannels"].get<int>(), 4 );
  BOOST_CHECK_CLOSE( opts["MdaConfidenceLevel"].get<double>(), 0.95, 1.0E-6 );
  BOOST_CHECK_CLOSE( opts["MdaConfidenceLevelPercent"].get<double>(), 95.0, 1.0E-6 );
}//BOOST_AUTO_TEST_CASE( JsonHasMdaOnlyForNotFitPeaks )


BOOST_AUTO_TEST_CASE( DefaultTemplatesRender )
{
  set_data_dir();

  // The default report templates must render both with, and without, the detection limit
  //  information - the latter being what a user re-rendering older results, or using the "none"
  //  method, will have.
  BatchPeak::BatchPeakFitOptions with_mda = default_options();
  with_mda.not_fit_peak_mda = BatchPeak::NotFitPeakMdaMethod::Currie;
  with_mda.template_include_dir = "default";

  BatchPeak::BatchPeakFitOptions without_mda = with_mda;
  without_mda.not_fit_peak_mda = BatchPeak::NotFitPeakMdaMethod::None;

  const string tmplt_dir = BatchInfoLog::template_include_dir( with_mda );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_directory(tmplt_dir),
                        "Report template directory not at '" << tmplt_dir << "'" );

  inja::Environment env = BatchInfoLog::get_default_inja_env( tmplt_dir );

  for( const BatchPeak::BatchPeakFitOptions &options : { with_mda, without_mda } )
  {
    const BatchPeak::BatchPeakFitResult result
          = run_fit( {661.0, 1173.0}, {{661.0, 5000.0}}, options );
    BOOST_REQUIRE( result.success );
    BOOST_REQUIRE_EQUAL( result.unfit_exemplar_peaks.size(), 1 );

    nlohmann::json data;
    BatchInfoLog::add_peak_fit_results_to_json( data, result );
    BatchInfoLog::add_peak_fit_options_to_json( data, options );
    BatchInfoLog::add_exe_info_to_json( data );
    data["Filename"] = "input.n42";
    data["ExemplarFile"] = "exemplar.n42";

    // The HTML reports embed the spectrum chart
    for( const pair<string,string> &key_val : BatchInfoLog::load_spectrum_chart_js_and_css() )
      data[key_val.first] = key_val.second;

    // The summary templates loop over a "Files" array of the per-file data
    nlohmann::json summary = data;
    summary["Files"] = nlohmann::json::array( { data } );

    const auto render = [&env,&options]( const char *type,
                                        const BatchInfoLog::TemplateRenderType render_type,
                                        const nlohmann::json &json ) -> string {
      try
      {
        return BatchInfoLog::render_template( type, env, render_type, options, json );
      }catch( std::exception &e )
      {
        BOOST_ERROR( "Failed to render '" << type << "' template: " << e.what() );
      }
      return string();
    };//render lambda

    for( const char *type : { "txt", "html" } )
      render( type, BatchInfoLog::TemplateRenderType::PeakFitIndividual, data );

    for( const char *type : { "csv", "html" } )
      render( type, BatchInfoLog::TemplateRenderType::PeakFitSummary, summary );

    // The activity/shielding templates also render the not-fit peaks, and are guarded separately.
    //  These templates need a little more scaffolding than the peak-fit ones; `BatchActivity`
    //  always provides it.
    nlohmann::json act_data = data;
    BatchInfoLog::add_not_fit_peaks_to_act_shield_json( act_data, result );
    act_data["HasFitResults"] = false;
    act_data["HasErrorMessage"] = false;
    act_data["Success"] = true;
    for( const pair<string,string> &key_val : BatchInfoLog::load_shielding_fit_plot_js_and_css() )
      act_data[key_val.first] = key_val.second;

    nlohmann::json act_summary = act_data;
    act_summary["Files"] = nlohmann::json::array( { act_data } );

    for( const char *type : { "txt", "html" } )
      render( type, BatchInfoLog::TemplateRenderType::ActShieldIndividual, act_data );

    for( const char *type : { "csv", "html" } )
      render( type, BatchInfoLog::TemplateRenderType::ActShieldSummary, act_summary );

    // The limit information should actually make it into the rendered reports
    if( options.not_fit_peak_mda != BatchPeak::NotFitPeakMdaMethod::None )
    {
      const string txt = render( "txt", BatchInfoLog::TemplateRenderType::PeakFitIndividual, data );
      BOOST_CHECK_MESSAGE( txt.find("Not detected") != string::npos,
                          "Detection limit description missing from txt report:\n" << txt );

      const string html = render( "html", BatchInfoLog::TemplateRenderType::PeakFitIndividual, data );
      BOOST_CHECK( html.find("Not detected") != string::npos );

      const string summary_csv = render( "csv", BatchInfoLog::TemplateRenderType::PeakFitSummary, summary );
      BOOST_CHECK_MESSAGE( summary_csv.find("could not be fit") != string::npos,
                          "Not-fit peak section missing from summary CSV:\n" << summary_csv );

      const string summary_html = render( "html", BatchInfoLog::TemplateRenderType::PeakFitSummary, summary );
      BOOST_CHECK( summary_html.find("Not detected") != string::npos );
    }//if( limits were computed )
  }//for( const BatchPeak::BatchPeakFitOptions &options : { with_mda, without_mda } )
}//BOOST_AUTO_TEST_CASE( DefaultTemplatesRender )


BOOST_AUTO_TEST_CASE( SummaryCsvIsWellFormed )
{
  set_data_dir();

  // Descriptions contain commas (the caveat sentences do), so the field must be quoted with the
  //  quote immediately after the delimiter - a space before it defeats CSV quoting, and every row
  //  whose description has a comma then gets extra fields.
  BatchPeak::BatchPeakFitOptions options = default_options();
  options.not_fit_peak_mda = BatchPeak::NotFitPeakMdaMethod::Currie;
  options.template_include_dir = "default";

  // Put two exemplar peaks close enough that one lands in the other's evaluated region, so the
  //  "a fit peak overlaps..." caveat - which contains a comma - is emitted.
  const BatchPeak::BatchPeakFitResult result
        = run_fit( {661.0, 673.0, 1173.0}, {{661.0, 5000.0}}, options );
  BOOST_REQUIRE( !result.not_fit_peak_mdas.empty() );

  bool any_comma_in_description = false;
  for( const BatchPeak::NotFitPeakMda &mda : result.not_fit_peak_mdas )
    any_comma_in_description |= (mda.description.find(',') != string::npos);
  BOOST_REQUIRE_MESSAGE( any_comma_in_description,
                        "Test needs a description containing a comma to be meaningful." );

  nlohmann::json data;
  BatchInfoLog::add_peak_fit_results_to_json( data, result );
  BatchInfoLog::add_peak_fit_options_to_json( data, options );
  BatchInfoLog::add_exe_info_to_json( data );
  data["Filename"] = "input.n42";

  nlohmann::json summary = data;
  summary["Files"] = nlohmann::json::array( { data } );
  summary["ExemplarFile"] = "exemplar.n42";

  inja::Environment env = BatchInfoLog::get_default_inja_env( BatchInfoLog::template_include_dir(options) );
  const string csv = BatchInfoLog::render_template( "csv", env,
                          BatchInfoLog::TemplateRenderType::PeakFitSummary, options, summary );

  vector<string> lines;
  SpecUtils::split( lines, csv, "\n" );

  // Find the not-fit section, then check every row has the same number of fields as its header
  size_t header_index = lines.size();
  for( size_t i = 0; i < lines.size(); ++i )
  {
    if( lines[i].find("could not be fit") != string::npos )
      header_index = i + 1;
  }

  BOOST_REQUIRE_MESSAGE( header_index < lines.size(), "Not-fit peak section missing from CSV" );

  const size_t num_columns = split_csv_line( lines[header_index] ).size();
  BOOST_CHECK( num_columns > 5 );

  size_t num_rows = 0;
  for( size_t i = header_index + 1; i < lines.size(); ++i )
  {
    const string &line = lines[i];
    if( line.empty() || (line.find(',') == string::npos) )
      break;

    num_rows += 1;
    BOOST_CHECK_MESSAGE( split_csv_line(line).size() == num_columns,
                        "CSV row " << num_rows << " has " << split_csv_line(line).size()
                        << " fields, but the header has " << num_columns << ": " << line );
  }//for( loop over rows )

  BOOST_CHECK_MESSAGE( num_rows == result.not_fit_peak_mdas.size(),
                      "Expected " << result.not_fit_peak_mdas.size() << " CSV rows, got " << num_rows
                      << " - rows may have run together onto one line." );
}//BOOST_AUTO_TEST_CASE( SummaryCsvIsWellFormed )


BOOST_AUTO_TEST_CASE( CsvSplitterMatchesRfc4180 )
{
  // The helper `SummaryCsvIsWellFormed` relies on must treat a field whose quote does not
  //  immediately follow the delimiter as UNQUOTED - that leading space is exactly the mistake the
  //  test exists to catch, and a lenient splitter would pass either way.
  BOOST_CHECK_EQUAL( split_csv_line("a,b,c").size(), 3 );
  BOOST_CHECK_EQUAL( split_csv_line("a,b,\"x, y\"").size(), 3 );          //correctly quoted
  BOOST_CHECK_MESSAGE( split_csv_line("a,b, \"x, y\"").size() == 4,       //space before the quote
                      "The splitter treats ', \"x, y\"' as a quoted field, so it cannot catch the"
                      " CSV-quoting bug it was written to guard against." );
}//BOOST_AUTO_TEST_CASE( CsvSplitterMatchesRfc4180 )


BOOST_AUTO_TEST_CASE( EmptyPeakRegion )
{
  set_data_dir();

  // A peak region holding no counts makes the Gaussian statistics degenerate: Lc and the upper
  //  limit both collapse to zero.  Left unguarded that reports "less than 0 counts", and a single
  //  stray count reads as a detection.
  BatchPeak::BatchPeakFitOptions options = default_options();
  options.not_fit_peak_mda = BatchPeak::NotFitPeakMdaMethod::Currie;

  // Counts below 1000 keV, nothing above it - so the 661 keV peak fits normally while the
  //  1173 keV region is empty.  (An entirely empty spectrum instead trips a pre-existing Minuit
  //  assertion in the peak fitter, which is a separate issue.)
  auto counts = make_shared<vector<float>>( sm_num_channels, 0.0f );
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( sm_num_channels, { 0.0f, sm_channel_width }, {} );

  const double sigma = sm_fwhm / 2.35482;
  for( size_t i = 0; i < sm_num_channels; ++i )
  {
    const double lower = cal->energy_for_channel( i );
    const double upper = cal->energy_for_channel( i + 1 );
    if( upper > 1000.0 )
      break;

    const double frac = 0.5*( erf( (upper - 661.0)/(sqrt(2.0)*sigma) )
                              - erf( (lower - 661.0)/(sqrt(2.0)*sigma) ) );
    (*counts)[i] = sm_continuum_cps*sm_live_time + static_cast<float>( 5000.0*frac );
  }

  auto m = make_shared<SpecUtils::Measurement>();
  m->set_gamma_counts( counts, sm_live_time, sm_live_time );
  m->set_energy_calibration( cal );
  m->set_sample_number( 1 );
  m->set_detector_name( "Det1" );
  m->set_source_type( SpecUtils::SourceType::Foreground );

  auto measurement = make_shared<SpecMeas>();
  measurement->add_measurement( m, false );
  measurement->cleanup_after_load();

  const shared_ptr<SpecMeas> exemplar = make_exemplar( {661.0, 1173.0} );
  const BatchPeak::BatchPeakFitResult result
        = BatchPeak::fit_peaks_in_file( "exemplar.n42", {}, exemplar, "empty.n42", measurement,
                                        {}, options );

  BOOST_REQUIRE( result.success );
  BOOST_REQUIRE( !result.not_fit_peak_mdas.empty() );

  for( const BatchPeak::NotFitPeakMda &mda : result.not_fit_peak_mdas )
  {
    BOOST_REQUIRE( mda.currie_computed );
    BOOST_CHECK_MESSAGE( mda.region_is_empty,
                        "Empty region not flagged at " << mda.exemplar_peak->mean() << " keV" );

    // Never a detection, whatever the classification arithmetic says
    BOOST_CHECK( mda.result_type == BatchPeak::NotFitPeakMda::MdaResultType::NotDetected );
    BOOST_CHECK_EQUAL( mda.short_description, string("No counts in region") );

    // Must not claim "less than 0 counts"; Ld is still meaningful and should be quoted
    BOOST_CHECK( mda.description.find("less than 0 counts") == string::npos );
    BOOST_CHECK( mda.description.find("Minimum reliably detectable") != string::npos );
    BOOST_CHECK( mda.currie_result.detection_limit > 0.0f );
  }//for( const BatchPeak::NotFitPeakMda &mda : result.not_fit_peak_mdas )
}//BOOST_AUTO_TEST_CASE( EmptyPeakRegion )


BOOST_AUTO_TEST_CASE( OptionsReportedEvenWhenLimitFails )
{
  set_data_dir();

  // A peak too close to the start of the spectrum for the side channels to fit; the options that
  //  were used must still be reported, rather than the zeros of a default-constructed result.
  BatchPeak::BatchPeakFitOptions options = default_options();
  options.not_fit_peak_mda = BatchPeak::NotFitPeakMdaMethod::Currie;
  options.mda_confidence_level = 0.99;
  options.mda_num_side_channels = 6;

  const BatchPeak::BatchPeakFitResult result = run_fit( {2.0, 661.0}, {{661.0, 5000.0}}, options );
  BOOST_REQUIRE_EQUAL( result.not_fit_peak_mdas.size(), 1 );
  BOOST_REQUIRE( !result.not_fit_peak_mdas.front().currie_computed );

  nlohmann::json data;
  BatchInfoLog::add_peak_fit_results_to_json( data, result );

  const nlohmann::json &mda = data["NotFitPeaks"]["Peaks"][0]["Mda"];
  BOOST_CHECK_EQUAL( mda["ResultType"].get<string>(), string("Error") );
  BOOST_CHECK_EQUAL( mda["ConfidenceLevelStr"].get<string>(), string("99%") );
  BOOST_CHECK_CLOSE( mda["ConfidenceLevel"].get<double>(), 0.99, 1.0E-6 );
  BOOST_CHECK_EQUAL( mda["NumSideChannels"].get<int>(), 6 );
}//BOOST_AUTO_TEST_CASE( OptionsReportedEvenWhenLimitFails )


BOOST_AUTO_TEST_CASE( MdaMethodStringRoundTrip )
{
  BOOST_CHECK_EQUAL( string(BatchPeak::to_str(BatchPeak::NotFitPeakMdaMethod::None)), string("none") );
  BOOST_CHECK_EQUAL( string(BatchPeak::to_str(BatchPeak::NotFitPeakMdaMethod::Currie)), string("currie") );
  BOOST_CHECK_EQUAL( string(BatchPeak::to_str(BatchPeak::NotFitPeakMdaMethod::CurrieAndDecon)),
                    string("currie-and-decon") );

  BOOST_CHECK( BatchPeak::not_fit_peak_mda_method_from_str("none") == BatchPeak::NotFitPeakMdaMethod::None );
  BOOST_CHECK( BatchPeak::not_fit_peak_mda_method_from_str("Currie") == BatchPeak::NotFitPeakMdaMethod::Currie );
  BOOST_CHECK( BatchPeak::not_fit_peak_mda_method_from_str(" currie-and-decon " )
              == BatchPeak::NotFitPeakMdaMethod::CurrieAndDecon );

  BOOST_CHECK_THROW( BatchPeak::not_fit_peak_mda_method_from_str("bogus"), std::exception );
}//BOOST_AUTO_TEST_CASE( MdaMethodStringRoundTrip )
