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

#include <map>
#include <set>
#include <cmath>
#include <chrono>
#include <fstream>
#include <deque>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <vector>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <algorithm>

#include <Wt/Utils>

#ifdef _WIN32
// For some reason, we need to include the following includes, before unit_test.hpp,
//  or we get a bunch of errors relating to winsock.k being included multiple times,
//  or something
#include "winsock2.h"
#include "Windows.h"
#endif

//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ShieldingSourceFitCalc_suite
//#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectionLimitCalc.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using namespace boost::unit_test;


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
}// void set_data_dir()

/** Builds a deterministic, isolated Gaussian on a flat continuum for exercising the profile scan.
 The fixed-geometry detector has unity efficiency, so activity is numerically equal to peak counts.
 */
shared_ptr<const DetectionLimitCalc::DeconComputeInput> make_decon_input(
  const double signal_counts,
  const double continuum_counts,
  const DetectionLimitCalc::DeconContinuumNorm continuum_norm = DetectionLimitCalc::DeconContinuumNorm::Floating )
{
  const size_t num_channels = 1024;
  const float channel_width = 2.0f;
  const float peak_energy = 1173.0f;
  const float fwhm = 8.0f;
  const double sigma = fwhm / PhysicalUnits::fwhm_nsigma;

  shared_ptr<vector<float>> counts = make_shared<vector<float>>( num_channels, continuum_counts );
  shared_ptr<SpecUtils::EnergyCalibration> calibration = make_shared<SpecUtils::EnergyCalibration>();
  calibration->set_polynomial( num_channels, { 0.0f, channel_width }, {} );

  for( size_t channel = 0; channel < num_channels; ++channel )
  {
    const double lower = calibration->energy_for_channel( channel );
    const double upper = calibration->energy_for_channel( channel + 1 );
    const double fraction = 0.5 * ( erf( ( upper - peak_energy ) / ( sqrt( 2.0 ) * sigma ) ) -
                                    erf( ( lower - peak_energy ) / ( sqrt( 2.0 ) * sigma ) ) );
    counts->at( channel ) += static_cast<float>( signal_counts * fraction );
  }

  shared_ptr<SpecUtils::Measurement> measurement = make_shared<SpecUtils::Measurement>();
  measurement->set_gamma_counts( counts, 1.0f, 1.0f );
  measurement->set_energy_calibration( calibration );

  shared_ptr<DetectorPeakResponse> drf = make_shared<DetectorPeakResponse>();
  drf->setIntrinsicEfficiencyFormula( "1.0",
                                      2.54f * PhysicalUnits::cm,
                                      PhysicalUnits::keV,
                                      0.0f,
                                      0.0f,
                                      DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct );
  drf->setFwhmCoefficients( vector<float>{ fwhm * fwhm, 0.0f },
                            DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial );

  DetectionLimitCalc::DeconRoiInfo roi;
  roi.roi_start = 1162.0f;
  roi.roi_end = 1184.0f;
  roi.continuum_type = PeakContinuum::OffsetType::Linear;
  roi.cont_norm_method = continuum_norm;
  roi.num_lower_side_channels = 4;
  roi.num_upper_side_channels = 4;

  DetectionLimitCalc::DeconRoiInfo::PeakInfo peak;
  peak.energy = peak_energy;
  peak.fwhm = fwhm;
  peak.counts_per_bq_into_4pi = 1.0;
  roi.peak_infos.push_back( peak );

  shared_ptr<DetectionLimitCalc::DeconComputeInput> input = make_shared<DetectionLimitCalc::DeconComputeInput>();
  input->activity = 0.0;
  input->distance = 0.0;
  input->include_air_attenuation = false;
  input->shielding_thickness = 0.0;
  input->drf = drf;
  input->measurement = measurement;
  input->roi_info.push_back( roi );

  return input;
}// make_decon_input(...)

/** Builds a background-reference input whose reference spectrum was counted for
 \p reference_seconds, and which predicts a measurement of \p sample_seconds.

 The channel counts and `counts_per_bq_into_4pi` both scale with the reference exposure, exactly as
 they would for a real spectrum, so "activity" keeps one meaning as the reference is lengthened and
 only its *statistical precision* changes.  That is what lets a test vary the reference/sample
 exposure ratio without also moving the quantity being measured.
 */
shared_ptr<const DetectionLimitCalc::DeconComputeInput>
make_background_reference_input( const double continuum_per_second,
                                 const double reference_seconds,
                                 const double sample_seconds )
{
  shared_ptr<DetectionLimitCalc::DeconComputeInput> input
      = make_shared<DetectionLimitCalc::DeconComputeInput>(
          *make_decon_input( 0.0, continuum_per_second*reference_seconds ) );

  shared_ptr<SpecUtils::Measurement> measurement = make_shared<SpecUtils::Measurement>();
  measurement->set_gamma_counts( make_shared<vector<float>>( *input->measurement->gamma_counts() ),
                                 static_cast<float>(reference_seconds),
                                 static_cast<float>(reference_seconds) );
  measurement->set_energy_calibration( input->measurement->energy_calibration() );
  input->measurement = measurement;

  // Counts per Bq are quoted against the live time of the spectrum they came from.
  input->roi_info.front().peak_infos.front().counts_per_bq_into_4pi = reference_seconds;

  input->measurement_model = DetectionLimitCalc::DeconMeasurementModel::BackgroundReference;
  input->sample_exposure = sample_seconds;

  return input;
}// make_background_reference_input(...)

DetectionLimitCalc::CurrieMdaResult
currie_for_decon_input( const shared_ptr<const DetectionLimitCalc::DeconComputeInput> &input )
{
  DetectionLimitCalc::CurrieMdaInput currie;
  currie.spectrum = input->measurement;
  currie.gamma_energy = input->roi_info.front().peak_infos.front().energy;
  currie.roi_lower_energy = input->roi_info.front().roi_start;
  currie.roi_upper_energy = input->roi_info.front().roi_end;
  currie.num_lower_side_channels = input->roi_info.front().num_lower_side_channels;
  currie.num_upper_side_channels = input->roi_info.front().num_upper_side_channels;
  currie.detection_probability = 0.95;
  currie.additional_uncertainty = 0.0;
  return DetectionLimitCalc::currie_mda_calc( currie );
}// currie_for_decon_input(...)

/** Suppresses the verbose profile-scan diagnostics while a grid or Monte Carlo study is running. */
class ScopedCoutSilence
{
public:
  ScopedCoutSilence() : m_previous( cout.rdbuf( m_sink.rdbuf() ) )
  {}

  ~ScopedCoutSilence()
  {
    cout.rdbuf( m_previous );
  }

private:
  ostringstream m_sink;
  streambuf *m_previous;
};// class ScopedCoutSilence

struct StudyProfileResult
{
  bool found_upper = false;
  bool found_lower = false;
  double best = 0.0;
  double upper = 0.0;
  /** The crossing below #best; only filled when a lower root was asked for and #found_lower. */
  double lower = 0.0;
};// struct StudyProfileResult

struct StudyPeakModel
{
  vector<double> observed;
  vector<double> variance;
  vector<double> relative_energy;
  vector<double> peak_fraction;
  vector<double> edge_expected;
  size_t num_score_bins = 0;

  /** Channel lower energies, `observed.size()+1` of them, and the channel-integrated polynomial
   basis over them - `basis0` is the channel width and `basis1` the integral of `(E - reference)`.
   Continuum coefficients in this model are therefore in exactly the same units production uses,
   which is what lets the Cash harness call the production optimizer directly.
   */
  vector<float> edges;
  vector<double> basis0;
  vector<double> basis1;
  double reference_energy = 0.0;
  double min_expected = 0.0;

  /** The channels either side of the ROI, exactly as production's `FixedByEdges` sees them.

   Kept as `PoissonChannel`s so the Cash harness can build the same sideband constraint production
   does by calling the same optimizer, rather than re-deriving a line by hand.  `edge_expected`
   above is the *old* frozen-line construction, retained only for the legacy Neyman baseline.
   */
  vector<DetectionLimitCalc::PoissonChannel> sidebands;
};// struct StudyPeakModel

pair<double, double> study_weighted_continuum( const StudyPeakModel &model, const double activity )
{
  double s00 = 0.0, s01 = 0.0, s11 = 0.0, t0 = 0.0, t1 = 0.0;
  for( size_t index = 0; index < model.observed.size(); ++index )
  {
    const double weight = 1.0 / model.variance[index];
    const double x = model.relative_energy[index];
    // `fit_amp_and_offset_imp(...)` clips the data remaining after fixed peaks are subtracted.
    // Reproduce that production behavior exactly here; the review discusses its statistical
    // consequence separately.
    const double residual = (std::max)( 0.0, model.observed[index] - activity * model.peak_fraction[index] );
    s00 += weight;
    s01 += weight * x;
    s11 += weight * x * x;
    t0 += weight * residual;
    t1 += weight * x * residual;
  }

  const double determinant = s00 * s11 - s01 * s01;
  if( fabs( determinant ) < 1.0E-12 )
    return { t0 / s00, 0.0 };

  return { ( t0 * s11 - t1 * s01 ) / determinant, ( t1 * s00 - t0 * s01 ) / determinant };
}// study_weighted_continuum(...)

/** The channels beside a scored region, clamped to the spectrum.

 Returns false for a side that has no channels available.  `first_channel - num_lower` on a
 `size_t` wraps to a huge number, after which `first <= last` is false and the range is silently
 empty - the caller then neither randomizes nor sums anything, and nothing says so.  Benign for the
 fixture this study has always used (the region starts around channel 581), and a live hazard for
 any relocated ROI, wider side-channel count, or second region.
 */
bool study_side_channel_range( const size_t first_channel,
                               const size_t last_channel,
                               const size_t num_lower,
                               const size_t num_upper,
                               const size_t num_channels,
                               pair<size_t,size_t> &lower_range,
                               pair<size_t,size_t> &upper_range )
{
  bool have_any = false;
  lower_range = { 1, 0 };   // an empty range, so an unchecked `for` loop does nothing
  upper_range = { 1, 0 };

  if( (num_lower > 0) && (first_channel > 0) )
  {
    lower_range.first = (first_channel > num_lower) ? (first_channel - num_lower) : 0;
    lower_range.second = first_channel - 1;
    have_any = true;
  }

  if( (num_upper > 0) && ((last_channel + 1) < num_channels) )
  {
    upper_range.first = last_channel + 1;
    upper_range.second = (std::min)( last_channel + num_upper, num_channels - 1 );
    have_any = true;
  }

  return have_any;
}// study_side_channel_range(...)

/** \param want_legacy Builds `edge_expected`, the frozen straight line the pre-Cash Neyman score
        needs.  Only `study_neyman_*` reads it, and it costs two side-channel sums plus a
        least-squares solve on every call - which the Monte Carlo would otherwise pay on each of
        its ~100,000 trials for a field it never looks at.
 */
StudyPeakModel make_study_peak_model( const DetectionLimitCalc::DeconComputeInput &input,
                                      const bool want_legacy = true )
{
  const DetectionLimitCalc::DeconRoiInfo &roi = input.roi_info.front();
  const DetectionLimitCalc::DeconRoiInfo::PeakInfo &peak = roi.peak_infos.front();
  const shared_ptr<const SpecUtils::EnergyCalibration> calibration = input.measurement->energy_calibration();

  // Production resolves the ROI with `round_roi_to_channels(...)` and uses that one range for both
  // the continuum fit and the score; they used to differ by a channel.
  const pair<size_t, size_t> roi_channels =
    DetectionLimitCalc::round_roi_to_channels( input.measurement, roi.roi_start, roi.roi_end );
  const size_t first_channel = roi_channels.first;
  const size_t fit_last_channel = roi_channels.second;
  const double sigma = peak.fwhm / PhysicalUnits::fwhm_nsigma;

  StudyPeakModel model;
  model.num_score_bins = 1 + fit_last_channel - first_channel;
  model.reference_energy = peak.energy;
  double observed_sum = 0.0;
  for( size_t channel = first_channel; channel <= fit_last_channel; ++channel )
  {
    const double lower = calibration->energy_for_channel( channel );
    const double upper = calibration->energy_for_channel( channel + 1 );
    const double center = 0.5 * ( lower + upper );
    const double fraction = 0.5 * ( erf( ( upper - peak.energy ) / ( sqrt( 2.0 ) * sigma ) ) -
                                    erf( ( lower - peak.energy ) / ( sqrt( 2.0 ) * sigma ) ) );
    const double observed = input.measurement->gamma_channel_content( channel );
    model.observed.push_back( observed );
    model.variance.push_back( (std::max)( observed, 1.0 ) );
    model.relative_energy.push_back( center - peak.energy );
    model.peak_fraction.push_back( fraction );

    const double x0 = lower - peak.energy;
    const double x1 = upper - peak.energy;
    model.edges.push_back( static_cast<float>( lower ) );
    model.basis0.push_back( x1 - x0 );
    model.basis1.push_back( 0.5 * ( x1 * x1 - x0 * x0 ) );
    observed_sum += observed;
  }
  model.edges.push_back( static_cast<float>( calibration->energy_for_channel( fit_last_channel + 1 ) ) );
  // Taken from production rather than copied: a study that floors expected counts differently is
  //  characterizing a different estimator than the one that ships.
  model.min_expected = DetectionLimitCalc::min_expected_channel_counts(
                          observed_sum / static_cast<double>( model.observed.size() ) );

  // Needed by both the legacy construction and the sideband list below.
  const size_t num_lower = roi.num_lower_side_channels;
  const size_t num_upper = roi.num_upper_side_channels;

  if( want_legacy )
  {
  double lower_cont_bound = calibration->channel_for_energy( roi.roi_start );
  double upper_cont_bound = calibration->channel_for_energy( roi.roi_end );
  double fractional = lower_cont_bound - floor( lower_cont_bound );
  lower_cont_bound = ( fractional > 0.9 ) ? round( lower_cont_bound ) : floor( lower_cont_bound );
  fractional = upper_cont_bound - floor( upper_cont_bound );
  upper_cont_bound = ( fractional < 0.1 ) ? round( upper_cont_bound ) : floor( upper_cont_bound );
  upper_cont_bound = (std::max)( upper_cont_bound, 1.0 );

  const size_t roi_first = static_cast<size_t>( lower_cont_bound );
  const size_t roi_last = (std::max)( roi_first, static_cast<size_t>( upper_cont_bound - 1.0 ) );
  // The legacy frozen-line construction keeps its own region rounding (it is the "before" baseline
  //  and must not move), but takes its channel ranges through the same guard as everything else -
  //  the subtraction below used to underflow to a huge `size_t` for a region near channel zero.
  pair<size_t,size_t> legacy_lower, legacy_upper;
  study_side_channel_range( roi_first, roi_last, num_lower, num_upper,
                           input.measurement->num_gamma_channels(), legacy_lower, legacy_upper );
  const size_t lower_first = legacy_lower.first;
  const size_t lower_last = legacy_lower.second;
  const size_t upper_first = legacy_upper.first;
  const size_t upper_last = legacy_upper.second;
  double lower_sum = 0.0, upper_sum = 0.0;
  for( size_t channel = lower_first; channel <= lower_last; ++channel )
    lower_sum += input.measurement->gamma_channel_content( channel );
  for( size_t channel = upper_first; channel <= upper_last; ++channel )
    upper_sum += input.measurement->gamma_channel_content( channel );

  const double lower_low = calibration->energy_for_channel( lower_first ) - peak.energy;
  const double lower_up = calibration->energy_for_channel( lower_last + 1 ) - peak.energy;
  const double upper_low = calibration->energy_for_channel( upper_first ) - peak.energy;
  const double upper_up = calibration->energy_for_channel( upper_last + 1 ) - peak.energy;
  const double lower_width = lower_up - lower_low;
  const double upper_width = upper_up - upper_low;
  const double lower_square_difference = lower_up * lower_up - lower_low * lower_low;
  const double upper_square_difference = upper_up * upper_up - upper_low * upper_low;
  const double denominator = upper_square_difference - upper_width * lower_square_difference / lower_width;
  const double slope = ( fabs( denominator ) < FLT_EPSILON )
                         ? 0.0
                         : 2.0 * ( upper_sum - upper_width * lower_sum / lower_width ) / denominator;
  const double intercept = ( fabs( denominator ) < FLT_EPSILON )
                             ? 0.5 * lower_sum / lower_width + 0.5 * upper_sum / upper_width
                             : ( lower_sum - 0.5 * slope * lower_square_difference ) / lower_width;

  for( size_t index = 0; index < model.num_score_bins; ++index )
  {
    const size_t channel = first_channel + index;
    const double lower = calibration->energy_for_channel( channel ) - peak.energy;
    const double upper = calibration->energy_for_channel( channel + 1 ) - peak.energy;
    model.edge_expected.push_back( intercept * ( upper - lower ) + 0.5 * slope * ( upper * upper - lower * lower ) );
  }
  }//if( want_legacy )

  // The same side channels, kept in the form production's sideband constraint uses.  Taken
  // relative to the scored ROI (`first_channel`..`fit_last_channel`, from
  // `round_roi_to_channels`), so the harness anchors on exactly the channels production does.
  const size_t num_measurement_channels = input.measurement->num_gamma_channels();
  const auto append_sideband = [&]( const size_t from, const size_t to ){
    for( size_t channel = from; channel <= to; ++channel )
    {
      DetectionLimitCalc::PoissonChannel side;
      side.lower_energy = calibration->energy_for_channel( channel );
      side.upper_energy = calibration->energy_for_channel( channel + 1 );
      side.observed = (std::max)( 0.0,
                        static_cast<double>( input.measurement->gamma_channel_content( channel ) ) );
      side.fixed_signal = 0.0;
      side.continuum_scale = 1.0;
      model.sidebands.push_back( side );
    }
  };

  pair<size_t,size_t> sideband_lower, sideband_upper;
  study_side_channel_range( first_channel, fit_last_channel, num_lower, num_upper,
                           num_measurement_channels, sideband_lower, sideband_upper );
  append_sideband( sideband_lower.first, sideband_lower.second );
  append_sideband( sideband_upper.first, sideband_upper.second );

  return model;
}// make_study_peak_model(...)

double study_neyman_score( const StudyPeakModel &model,
                           const DetectionLimitCalc::DeconContinuumNorm norm,
                           const double activity )
{
  pair<double, double> continuum;
  switch( norm )
  {
  case DetectionLimitCalc::DeconContinuumNorm::Floating:
    continuum = study_weighted_continuum( model, activity );
    break;

  case DetectionLimitCalc::DeconContinuumNorm::FixedByFullRange:
    continuum = study_weighted_continuum( model, 0.0 );
    break;

  case DetectionLimitCalc::DeconContinuumNorm::FixedByEdges:
    continuum = { 0.0, 0.0 };
    break;
  }

  double score = 0.0;
  for( size_t index = 0; index < model.num_score_bins; ++index )
  {
    const double expected_continuum = ( norm == DetectionLimitCalc::DeconContinuumNorm::FixedByEdges )
                                        ? model.edge_expected[index]
                                        : continuum.first + continuum.second * model.relative_energy[index];
    const double expected = expected_continuum + activity * model.peak_fraction[index];
    const double difference = model.observed[index] - expected;
    score += difference * difference / model.variance[index];
  }
  return score;
}// study_neyman_score(...)

StudyProfileResult study_neyman_limit( const DetectionLimitCalc::DeconComputeInput &input,
                                       const double initial_maximum )
{
  const StudyPeakModel model = make_study_peak_model( input );
  const DetectionLimitCalc::DeconContinuumNorm norm = input.roi_info.front().cont_norm_method;
  const auto score = [&model, norm]( const double activity ) { return study_neyman_score( model, norm, activity ); };

  if( norm == DetectionLimitCalc::DeconContinuumNorm::Floating )
  {
    double lower = 0.0;
    double upper = (std::max)( 1.0, initial_maximum );
    for( size_t iteration = 0; iteration < 24; ++iteration )
    {
      const double one_third = lower + ( upper - lower ) / 3.0;
      const double two_thirds = upper - ( upper - lower ) / 3.0;
      if( score( one_third ) < score( two_thirds ) )
        upper = two_thirds;
      else
        lower = one_third;
    }

    StudyProfileResult result;
    result.best = 0.5 * ( lower + upper );
    const double score0 = score( 0.0 );
    if( score0 <= score( result.best ) )
      result.best = 0.0;
    const double best_score = score( result.best );
    const double delta = DetectionLimitCalc::decon_limit_delta(
                            0.95, DetectionLimitCalc::DeconLimitType::OneSidedUpperLimit );
    double root_upper = (std::max)( initial_maximum, result.best + 1.0 );
    for( size_t expansion = 0; expansion < 20; ++expansion )
    {
      if( ( score( root_upper ) - best_score ) >= delta )
        break;
      root_upper *= 2.0;
    }
    if( ( score( root_upper ) - best_score ) < delta )
      return result;

    lower = result.best;
    for( size_t iteration = 0; iteration < 28; ++iteration )
    {
      const double midpoint = 0.5 * ( lower + root_upper );
      if( ( score( midpoint ) - best_score ) < delta )
        lower = midpoint;
      else
        root_upper = midpoint;
    }
    result.upper = 0.5 * ( lower + root_upper );
    result.found_upper = true;
    result.found_lower = ( result.best > 1.0E-8 ) && ( ( score0 - best_score ) >= delta );
    return result;
  }

  const double score0 = score( 0.0 );
  const double score1 = score( 1.0 );
  const double score2 = score( 2.0 );
  const double quadratic = 0.5 * ( score2 - 2.0 * score1 + score0 );
  const double linear = score1 - score0 - quadratic;

  StudyProfileResult result;
  if( quadratic <= 1.0E-12 )
    return result;

  result.best = (std::max)( 0.0, -linear / ( 2.0 * quadratic ) );
  const double best_score = quadratic * result.best * result.best + linear * result.best + score0;
  const double delta = DetectionLimitCalc::decon_limit_delta(
                          0.95, DetectionLimitCalc::DeconLimitType::OneSidedUpperLimit );
  const double discriminant = linear * linear - 4.0 * quadratic * ( score0 - best_score - delta );
  if( discriminant < 0.0 )
    return result;

  result.upper = ( -linear + sqrt( discriminant ) ) / ( 2.0 * quadratic );
  result.found_upper = ( result.upper >= result.best );
  result.found_lower = ( result.best > 1.0E-8 ) && ( ( score0 - best_score ) >= delta );
  return result;
}// study_neyman_limit(...)

/** The Poisson deviance of the model at the given continuum coefficients, in production's
 channel-integrated parameterization, using the production per-channel statistic.
 */
double study_cash_objective( const StudyPeakModel &model,
                             const double activity,
                             const pair<double, double> &continuum,
                             const size_t num_bins )
{
  double score = 0.0;
  for( size_t index = 0; index < num_bins; ++index )
  {
    const double expected =
      (std::max)( model.min_expected, continuum.first * model.basis0[index] +
                                        continuum.second * model.basis1[index] +
                                        activity * model.peak_fraction[index] );
    score += DetectionLimitCalc::poisson_deviance( (std::max)( model.observed[index], 0.0 ), expected );
  }
  return score;
}// study_cash_objective(...)

/** How the harness minimizes the joint objective for `FixedByEdges`.

 The two exist because the study wants opposite things in two places.
 */
enum class StudyOptimizer
{
  /** The harness's own coarse-to-fine coordinate search.

   Shares no code with the production continuum fit, which is the entire reason
   `DeconStudyProfileMatchesProduction` is evidence of anything: two independent implementations
   agreeing, rather than one checking itself.  About ten times slower than #Production, because it
   evaluates the joint objective some 640 times per activity.
   */
  Independent,

  /** Calls `fit_continuum_poisson` over the region-plus-sideband channel list, exactly as
   production does.

   What the Monte Carlo uses: the coordinate search costs 30-50 ms per trial limit, which the grid
   cannot afford tens of thousands of times.  This does narrow what the in-cell production audit
   proves - with this optimizer the audit compares two searches over the same fit, not two fits -
   so `DeconStudyProfileMatchesProduction` remains the test that covers the fit itself.
   */
  Production
};//enum class StudyOptimizer


/** Numerical health of one profiled limit, accumulated over the whole activity scan.

 Kept apart from the scan's own success/failure because they are different defects with different
 fixes: a continuum fit that does not converge is an optimizer problem, while a profile that never
 crosses the threshold is a search-range problem.  The study used to conflate them.
 */
struct StudyDiagnostics
{
  size_t score_evaluations = 0;
  size_t continuum_iterations = 0;
  size_t continuum_restarts = 0;
  size_t continuum_failures = 0;
  /** Coefficients of the LAST continuum fit this scan performed - a level, not a bias.

   The scan evaluates hundreds of activities, so which one this is depends on where the bisection
   finished; it moves with the probed activity and nothing is subtracted from it.  A sanity check
   that the continuum is in the right neighbourhood, and no more than that.
   */
  pair<double,double> best_continuum{ 0.0, 0.0 };
};//struct StudyDiagnostics


/** Profiles the continuum against the Poisson likelihood by calling the production optimizer, so
 the study cannot drift away from the code it is meant to characterize.
 */
pair<double, double> study_cash_continuum( const StudyPeakModel &model, const double activity,
                                           StudyDiagnostics * const diagnostics = nullptr )
{
  vector<double> signal( model.observed.size(), 0.0 );
  for( size_t index = 0; index < model.observed.size(); ++index )
    signal[index] = activity * model.peak_fraction[index];

  const DetectionLimitCalc::PoissonContinuumFit fit = DetectionLimitCalc::fit_continuum_poisson(
    model.edges.data(), model.observed.data(), signal.data(), model.observed.size(), 2,
    model.reference_energy, {} );

  if( diagnostics )
  {
    diagnostics->continuum_iterations += fit.num_iterations;
    diagnostics->continuum_restarts += fit.num_restarts;
    diagnostics->continuum_failures += (fit.converged ? 0u : 1u);
    if( fit.converged && (fit.coefficients.size() >= 2) )
      diagnostics->best_continuum = { fit.coefficients[0], fit.coefficients[1] };
  }

  // A non-converged fit used to be swallowed as a zero continuum, which is not a continuum at all -
  //  it silently produces a wrong limit that looks like a successful scan.  Still returned so the
  //  scan can proceed, but `diagnostics` now records that it happened.
  if( !fit.converged )
    return { 0.0, 0.0 };

  return { fit.coefficients[0], fit.coefficients[1] };
}// study_cash_continuum(...)

double study_cash_score( const StudyPeakModel &model,
                         const DetectionLimitCalc::DeconContinuumNorm norm,
                         const pair<double, double> &full_range_continuum,
                         const double activity,
                         const StudyOptimizer optimizer = StudyOptimizer::Independent,
                         StudyDiagnostics * const diagnostics = nullptr )
{
  if( diagnostics )
    diagnostics->score_evaluations += 1;

  if( (norm == DetectionLimitCalc::DeconContinuumNorm::FixedByEdges)
     && (optimizer == StudyOptimizer::Production) )
  {
    // The same joint channel list production builds - region channels carrying the trial signal,
    //  then the side channels carrying none - handed to the production optimizer.
    vector<DetectionLimitCalc::PoissonChannel> channels;
    channels.reserve( model.num_score_bins + model.sidebands.size() );

    for( size_t index = 0; index < model.num_score_bins; ++index )
    {
      DetectionLimitCalc::PoissonChannel channel;
      channel.lower_energy = model.edges[index];
      channel.upper_energy = model.edges[index + 1];
      channel.observed = (std::max)( model.observed[index], 0.0 );
      channel.fixed_signal = activity * model.peak_fraction[index];
      channel.continuum_scale = 1.0;
      channels.push_back( channel );
    }

    channels.insert( channels.end(), model.sidebands.begin(), model.sidebands.end() );

    const DetectionLimitCalc::PoissonContinuumFit fit = DetectionLimitCalc::fit_continuum_poisson(
      channels.data(), channels.size(), 2, model.reference_energy, {} );

    if( diagnostics )
    {
      diagnostics->continuum_iterations += fit.num_iterations;
      diagnostics->continuum_restarts += fit.num_restarts;
      diagnostics->continuum_failures += (fit.converged ? 0u : 1u);
    }

    const double c0 = fit.converged ? fit.coefficients[0] : 0.0;
    const double c1 = fit.converged ? fit.coefficients[1] : 0.0;

    double score = study_cash_objective( model, activity, { c0, c1 }, model.num_score_bins );
    for( const DetectionLimitCalc::PoissonChannel &side : model.sidebands )
    {
      const double x0 = side.lower_energy - model.reference_energy;
      const double x1 = side.upper_energy - model.reference_energy;
      const double expected = (std::max)( model.min_expected,
                                          c0*(x1 - x0) + c1*0.5*(x1*x1 - x0*x0) );
      score += DetectionLimitCalc::poisson_deviance( (std::max)( side.observed, 0.0 ), expected );
    }

    if( diagnostics )
      diagnostics->best_continuum = { c0, c1 };

    return score;
  }//if( FixedByEdges via the production optimizer )

  if( norm == DetectionLimitCalc::DeconContinuumNorm::FixedByEdges )
  {
    // Production folds the side channels into the same Poisson likelihood as the region, carrying
    // no signal, and profiles the continuum over the union.  The harness models the same thing -
    // but minimizes it with its OWN coarse-to-fine coordinate search rather than calling
    // production's optimizer, so `DeconStudyProfileMatchesProduction` stays a comparison of two
    // independent implementations instead of one implementation with itself.
    const auto joint_score = [&]( const double c0, const double c1 ) -> double {
      double score = study_cash_objective( model, activity, { c0, c1 }, model.num_score_bins );
      for( const DetectionLimitCalc::PoissonChannel &side : model.sidebands )
      {
        const double x0 = side.lower_energy - model.reference_energy;
        const double x1 = side.upper_energy - model.reference_energy;
        const double expected = (std::max)( model.min_expected,
                                            c0*(x1 - x0) + c1*0.5*(x1*x1 - x0*x0) );
        score += DetectionLimitCalc::poisson_deviance( (std::max)( side.observed, 0.0 ), expected );
      }
      return score;
    };

    pair<double, double> best = study_cash_continuum( model, activity, diagnostics );
    double best_score = joint_score( best.first, best.second );
    double step0 = (std::max)( 1.0, 0.25*fabs( best.first ) );
    double step1 = (std::max)( 1.0E-3, 0.25*fabs( best.second ) );

    for( size_t refine = 0; refine < 80; ++refine )
    {
      bool improved = false;
      for( const double d0 : { -step0, 0.0, step0 } )
      {
        for( const double d1 : { -step1, 0.0, step1 } )
        {
          if( ( d0 == 0.0 ) && ( d1 == 0.0 ) )
            continue;

          const double trial = joint_score( best.first + d0, best.second + d1 );
          if( trial < best_score )
          {
            best_score = trial;
            best = { best.first + d0, best.second + d1 };
            improved = true;
          }
        }
      }

      if( !improved )
      {
        step0 *= 0.5;
        step1 *= 0.5;
      }
    }

    if( diagnostics )
      diagnostics->best_continuum = best;

    return best_score;
  }


  pair<double, double> continuum;
  switch( norm )
  {
  case DetectionLimitCalc::DeconContinuumNorm::Floating:
    continuum = study_cash_continuum( model, activity, diagnostics );
    break;

  case DetectionLimitCalc::DeconContinuumNorm::FixedByFullRange:
    continuum = full_range_continuum;
    break;

  case DetectionLimitCalc::DeconContinuumNorm::FixedByEdges:
    continuum = { 0.0, 0.0 };
    break;
  }
  return study_cash_objective( model, activity, continuum, model.num_score_bins );
}// study_cash_score(...)

/** Profiles activity against the Cash statistic and returns the crossings of \p delta.

 \p delta comes from `DetectionLimitCalc::decon_limit_delta(...)` rather than a literal, so a change
 to how the threshold is defined cannot leave the study still measuring the old one - which is what
 three separate copies of 2.705543454095404 in this file used to allow.

 \p want_lower_root additionally locates the crossing below the best fit, which is what a central
 interval needs; the one-sided cases do not pay for it.
 */
StudyProfileResult study_cash_limit( const DetectionLimitCalc::DeconComputeInput &input,
                                     const double initial_maximum,
                                     const double delta,
                                     const bool want_lower_root = false,
                                     const StudyOptimizer optimizer = StudyOptimizer::Independent,
                                     StudyDiagnostics * const diagnostics = nullptr )
{
  // The legacy frozen-line construction is only needed by `study_neyman_*`; skipping it here saves
  //  building it on every one of the Monte Carlo's trials.
  const bool want_legacy = false;
  const StudyPeakModel model = make_study_peak_model( input, want_legacy );
  const DetectionLimitCalc::DeconContinuumNorm norm = input.roi_info.front().cont_norm_method;
  const pair<double, double> full_range_continuum = study_cash_continuum( model, 0.0, diagnostics );
  const auto score = [&]( const double activity )
  { return study_cash_score( model, norm, full_range_continuum, activity, optimizer, diagnostics ); };

  double lower = 0.0;
  double upper = (std::max)( 1.0, initial_maximum );
  for( size_t iteration = 0; iteration < 36; ++iteration )
  {
    const double one_third = lower + ( upper - lower ) / 3.0;
    const double two_thirds = upper - ( upper - lower ) / 3.0;
    if( score( one_third ) < score( two_thirds ) )
      upper = two_thirds;
    else
      lower = one_third;
  }

  StudyProfileResult result;
  result.best = 0.5 * ( lower + upper );
  if( score( 0.0 ) <= score( result.best ) )
    result.best = 0.0;
  const double best_score = score( result.best );

  double root_upper = (std::max)( initial_maximum, result.best + 1.0 );
  for( size_t expansion = 0; expansion < 20; ++expansion )
  {
    if( ( score( root_upper ) - best_score ) >= delta )
      break;
    root_upper *= 2.0;
  }

  if( ( score( root_upper ) - best_score ) >= delta )
  {
    double bracket_low = result.best;
    for( size_t iteration = 0; iteration < 36; ++iteration )
    {
      const double midpoint = 0.5 * ( bracket_low + root_upper );
      if( ( score( midpoint ) - best_score ) < delta )
        bracket_low = midpoint;
      else
        root_upper = midpoint;
    }
    result.upper = 0.5 * ( bracket_low + root_upper );
    result.found_upper = true;
  }//if( the upper threshold was bracketed )

  result.found_lower = ( result.best > 1.0E-8 ) && ( ( score( 0.0 ) - best_score ) >= delta );

  if( want_lower_root && result.found_lower )
  {
    // Between zero (above the threshold) and the best fit (below it), so the crossing is bracketed
    //  by construction and needs no expansion step.
    double bracket_low = 0.0, bracket_high = result.best;
    for( size_t iteration = 0; iteration < 36; ++iteration )
    {
      const double midpoint = 0.5 * ( bracket_low + bracket_high );
      if( ( score( midpoint ) - best_score ) >= delta )
        bracket_low = midpoint;
      else
        bracket_high = midpoint;
    }
    result.lower = 0.5 * ( bracket_low + bracket_high );
  }//if( want_lower_root && result.found_lower )

  return result;
}// study_cash_limit(...)


/** The one-sided 95% form every existing caller wants, so those call sites stay unchanged. */
StudyProfileResult study_cash_limit( const DetectionLimitCalc::DeconComputeInput &input,
                                     const double initial_maximum )
{
  const double delta = DetectionLimitCalc::decon_limit_delta(
                          0.95, DetectionLimitCalc::DeconLimitType::OneSidedUpperLimit );
  return study_cash_limit( input, initial_maximum, delta );
}// study_cash_limit(...)

/** A uniform on (0,1) taken straight from the engine's 32-bit output.

 `std::generate_canonical` and `std::uniform_real_distribution` are both free to consume a
 different number of words, so they do not give the same stream on every standard library; the
 engine's own output sequence is specified exactly.
 */
double study_uniform( mt19937 &generator )
{
  return ( static_cast<double>( generator() ) + 0.5 ) / 4294967296.0;
}// study_uniform(...)


/** A standard normal deviate, portable for the same reason `study_poisson` is.

 Box-Muller from two engine uniforms; no cached second variate, so a call consumes exactly two
 words of the stream whatever the standard library.
 */
double study_normal( mt19937 &generator )
{
  const double u1 = study_uniform( generator );
  const double u2 = study_uniform( generator );
  return std::sqrt( -2.0*std::log( u1 ) ) * std::cos( 2.0*3.14159265358979323846*u2 );
}// study_normal(...)


/** A Gamma(shape, 1) deviate, portable, for any shape > 0.

 Marsaglia-Tsang (2000) squeeze method, with their `shape < 1` boost.  Needed because
 `std::gamma_distribution` is as implementation-defined as `std::poisson_distribution`.
 */
double study_gamma( const double shape, mt19937 &generator )
{
  if( !(shape > 0.0) )
    return 0.0;

  if( shape < 1.0 )
  {
    // Gamma(a) == Gamma(a+1) * U^(1/a)
    const double g = study_gamma( shape + 1.0, generator );
    return g * std::pow( study_uniform( generator ), 1.0/shape );
  }

  const double d = shape - 1.0/3.0;
  const double c = 1.0 / std::sqrt( 9.0*d );

  for( size_t iteration = 0; iteration < 1000; ++iteration )
  {
    double x = 0.0, v = 0.0;
    do
    {
      x = study_normal( generator );
      v = 1.0 + c*x;
    }while( v <= 0.0 );

    v = v*v*v;
    const double u = study_uniform( generator );
    const double x2 = x*x;

    if( u < (1.0 - 0.0331*x2*x2) )
      return d*v;
    if( std::log(u) < (0.5*x2 + d*(1.0 - v + std::log(v))) )
      return d*v;
  }

  return d;   // unreachable in practice; keeps the function total
}// study_gamma(...)


/** A Poisson deviate with a fixed, portable algorithm.

 `std::poisson_distribution` may use any algorithm, and libstdc++, libc++ and MSVC do use different
 ones, so the same seed gives different spectra on different platforms.  The whole product of this
 study is a per-cell seed that lets someone reproduce a number on their own machine, which that
 makes impossible - hence our own sampler.

 Knuth's product method below `sm_knuth_mean_limit` (O(mean), fine for small means), and Hoermann's
 PTRS transformed-rejection above it (O(1), and what the small-mean method would be far too slow
 for at the hundreds of counts per channel the denser cells reach).  Both are exact.
 */
int study_poisson( const double mean, mt19937 &generator )
{
  const double sm_knuth_mean_limit = 10.0;

  if( !(mean > 0.0) )
    return 0;

  if( mean < sm_knuth_mean_limit )
  {
    const double limit = std::exp( -mean );
    double product = 1.0;
    int count = 0;
    for( ; count < 10000; ++count )
    {
      product *= study_uniform( generator );
      if( product <= limit )
        break;
    }
    return count;
  }//if( mean < sm_knuth_mean_limit )

  // Hoermann (1993), "The transformed rejection method for generating Poisson random variables".
  const double b = 0.931 + 2.53 * std::sqrt( mean );
  const double a = -0.059 + 0.02483 * b;
  const double inverse_alpha = 1.1239 + 1.1328 / ( b - 3.4 );
  const double v_r = 0.9277 - 3.6224 / ( b - 2.0 );

  for( size_t iteration = 0; iteration < 10000; ++iteration )
  {
    const double u = study_uniform( generator ) - 0.5;
    const double v = study_uniform( generator );
    const double us = 0.5 - std::fabs( u );
    const double k = std::floor( ( 2.0 * a / us + b ) * u + mean + 0.43 );

    if( (us >= 0.07) && (v <= v_r) )
      return static_cast<int>( k );

    if( (k < 0.0) || ((us < 0.013) && (v > us)) )
      continue;

    if( std::log( v * inverse_alpha / ( a / ( us * us ) + b ) )
       <= ( -mean + k * std::log( mean ) - std::lgamma( k + 1.0 ) ) )
    {
      return static_cast<int>( k );
    }
  }//for( bounded rejection loop )

  return static_cast<int>( mean );
}// study_poisson(...)




/** A Poisson realization of \p expected, over every region it describes.

 Randomizes each region's channels *and* its side channels: an extra region left at its exact
 expectation contributes a noiseless, perfectly-matching likelihood term, which inflates coverage
 for free.  Channels are marked first and sampled afterwards in ascending order, so a channel two
 regions share is drawn once - two draws would be two different realizations of one measurement -
 and so the draw order does not depend on the order the regions happen to be listed in.
 */
shared_ptr<const DetectionLimitCalc::DeconComputeInput>
make_poisson_trial( const shared_ptr<const DetectionLimitCalc::DeconComputeInput> &expected,
                    mt19937 &generator )
{
  const shared_ptr<const SpecUtils::Measurement> reference = expected->measurement;
  const size_t num_channels = reference->num_gamma_channels();
  shared_ptr<vector<float>> counts = make_shared<vector<float>>( *reference->gamma_counts() );

  vector<bool> randomize( num_channels, false );
  const auto mark = [&randomize,num_channels]( const size_t from, const size_t to ){
    for( size_t channel = from; (channel <= to) && (channel < num_channels); ++channel )
      randomize[channel] = true;
  };

  for( const DetectionLimitCalc::DeconRoiInfo &roi : expected->roi_info )
  {
    const size_t roi_first = reference->find_gamma_channel( roi.roi_start );
    const size_t roi_last = reference->find_gamma_channel( roi.roi_end - 0.0000001f );

    pair<size_t,size_t> lower_range, upper_range;
    study_side_channel_range( roi_first, roi_last, roi.num_lower_side_channels,
                             roi.num_upper_side_channels, num_channels, lower_range, upper_range );

    mark( roi_first, roi_last );
    mark( lower_range.first, lower_range.second );
    mark( upper_range.first, upper_range.second );
  }//for( loop over regions )

  for( size_t channel = 0; channel < num_channels; ++channel )
  {
    if( !randomize[channel] )
      continue;

    const double mean = (std::max)( 0.0,
                          static_cast<double>( reference->gamma_channel_content( channel ) ) );
    counts->at( channel ) = static_cast<float>( study_poisson( mean, generator ) );
  }

  shared_ptr<SpecUtils::Measurement> measurement = make_shared<SpecUtils::Measurement>();

  // The exposure is part of the input, not decoration: `decon_compute_peaks` forms
  //  `sample_exposure / measurement->live_time()`, so flattening it to one second silently
  //  reprojects a `BackgroundReference` trial to a different measurement than the one asked for.
  measurement->set_gamma_counts( counts, reference->live_time(), reference->real_time() );
  measurement->set_energy_calibration( reference->energy_calibration() );

  shared_ptr<DetectionLimitCalc::DeconComputeInput> trial =
    make_shared<DetectionLimitCalc::DeconComputeInput>( *expected );
  trial->measurement = measurement;
  return trial;
}// make_poisson_trial(...)

/** Formats a double for a cell descriptor, so the text is byte-identical on every platform. */
string study_number_token( const double value )
{
  char buffer[32] = { '\0' };
  snprintf( buffer, sizeof(buffer), "%.6g", value );
  return string( buffer );
}// study_number_token(...)


/** The seed for a study cell, derived from its descriptor.

 FNV-1a rather than `std::hash`, whose result is implementation-defined and differs between
 libstdc++ and libc++ - a "rerun this cell to reproduce the number" instruction that only works on
 the machine that produced it is not a reproduction instruction.  Unsigned wrap-around is
 well-defined, so the addition needs no guard.

 Deriving the seed from the same string the CSV prints as the cell key is what stops the two
 disagreeing about which cell a number belongs to.
 */
uint32_t study_cell_seed( const string &descriptor )
{
  const uint64_t offset_basis = 14695981039346656037ull;
  const uint64_t prime = 1099511628211ull;

  uint64_t hash = offset_basis;
  for( const char character : descriptor )
  {
    hash ^= static_cast<uint64_t>( static_cast<unsigned char>( character ) );
    hash *= prime;
  }

  return 20260807u + static_cast<uint32_t>( hash ^ ( hash >> 32 ) );
}// study_cell_seed(...)


pair<double, double> wilson_interval( const size_t successes, const size_t trials )
{
  const double z = 1.959963984540054;
  const double n = static_cast<double>( trials );
  const double fraction = static_cast<double>( successes ) / n;
  const double denominator = 1.0 + z * z / n;
  const double center = ( fraction + z * z / ( 2.0 * n ) ) / denominator;
  const double half_width = z * sqrt( ( fraction * ( 1.0 - fraction ) / n ) + z * z / ( 4.0 * n * n ) ) / denominator;
  return { center - half_width, center + half_width };
}// wilson_interval(...)

/** A small, deliberately naive description of a Poisson continuum-fitting problem, so the tests
 can build one without a `SpecUtils::Measurement`.
 */
struct ContinuumFitProblem
{
  vector<float> edges;         //!< nbin+1 channel lower energies
  vector<double> observed;     //!< nbin observed counts
  vector<double> signal;       //!< nbin fixed signal counts
  size_t num_coefficients = 2;
  double reference_energy = 0.0;

  size_t nbin() const { return observed.size(); }

  const double *signal_ptr() const { return signal.empty() ? nullptr : signal.data(); }
};// struct ContinuumFitProblem

/** Builds a problem with `nbin` uniform channels, a polynomial continuum given by `truth`, and an
 optional Gaussian signal.  Counts are the *expected* values (not Poisson sampled) unless a
 generator is supplied, which makes the tests deterministic by default.
 */
ContinuumFitProblem make_continuum_problem( const size_t nbin,
                                            const double channel_width,
                                            const double first_energy,
                                            const vector<double> &truth,
                                            const double signal_counts,
                                            const double signal_sigma,
                                            mt19937 *generator = nullptr )
{
  ContinuumFitProblem problem;
  problem.num_coefficients = truth.size();
  problem.reference_energy = first_energy + 0.5 * channel_width * static_cast<double>( nbin );

  for( size_t i = 0; i <= nbin; ++i )
    problem.edges.push_back( static_cast<float>( first_energy + channel_width * static_cast<double>( i ) ) );

  const double signal_mean = problem.reference_energy;

  for( size_t i = 0; i < nbin; ++i )
  {
    const double x0 = problem.edges[i] - problem.reference_energy;
    const double x1 = problem.edges[i + 1] - problem.reference_energy;

    double continuum = 0.0;
    double x0_power = x0, x1_power = x1;
    for( size_t k = 0; k < truth.size(); ++k )
    {
      continuum += truth[k] * ( x1_power - x0_power ) / static_cast<double>( k + 1 );
      x0_power *= x0;
      x1_power *= x1;
    }

    const double fraction =
      ( signal_counts > 0.0 )
        ? 0.5 * ( erf( ( x1 - ( signal_mean - problem.reference_energy ) ) / ( sqrt( 2.0 ) * signal_sigma ) ) -
                  erf( ( x0 - ( signal_mean - problem.reference_energy ) ) / ( sqrt( 2.0 ) * signal_sigma ) ) )
        : 0.0;
    const double signal = signal_counts * fraction;

    double counts = continuum + signal;
    if( generator )
    {
      poisson_distribution<int> distribution( (std::max)( counts, 0.0 ) );
      counts = distribution( *generator );
    }

    problem.observed.push_back( (std::max)( counts, 0.0 ) );
    problem.signal.push_back( signal );
  }

  return problem;
}// make_continuum_problem(...)

/** The Poisson deviance of `problem` at the given coefficients, computed independently of the
 production optimizer (though necessarily with the same basis definition).
 */
double reference_objective( const ContinuumFitProblem &problem, const vector<double> &coefs )
{
  double answer = 0.0;
  for( size_t i = 0; i < problem.nbin(); ++i )
  {
    const double x0 = problem.edges[i] - problem.reference_energy;
    const double x1 = problem.edges[i + 1] - problem.reference_energy;

    double expected = problem.signal.empty() ? 0.0 : problem.signal[i];
    double x0_power = x0, x1_power = x1;
    for( size_t k = 0; k < coefs.size(); ++k )
    {
      expected += coefs[k] * ( x1_power - x0_power ) / static_cast<double>( k + 1 );
      x0_power *= x0;
      x1_power *= x1;
    }

    if( expected <= 0.0 )
      return numeric_limits<double>::infinity();

    const double observed = problem.observed[i];
    answer += 2.0 * ( expected - observed + ( ( observed > 0.0 ) ? observed * std::log( observed / expected ) : 0.0 ) );
  }

  return answer;
}// reference_objective(...)

/** A deliberately slow, independent minimizer of `reference_objective`: cyclic coordinate descent
 with a golden-section line search and geometric step shrinking.  It shares no code with the
 production optimizer, so agreement between them is real evidence.
 */
vector<double> reference_minimize( const ContinuumFitProblem &problem, const vector<double> &start )
{
  vector<double> best = start;
  double best_value = reference_objective( problem, best );

  double step_scale = 1.0;
  for( size_t sweep = 0; sweep < 400; ++sweep )
  {
    bool improved = false;
    for( size_t k = 0; k < best.size(); ++k )
    {
      const double magnitude = (std::max)( fabs( best[k] ), 1.0E-6 );
      const double step = step_scale * 0.5 * magnitude;

      for( const double direction : { -1.0, 1.0 } )
      {
        vector<double> candidate = best;
        candidate[k] += direction * step;
        const double value = reference_objective( problem, candidate );
        if( value < best_value )
        {
          best = candidate;
          best_value = value;
          improved = true;
        }
      }
    }

    if( !improved )
    {
      step_scale *= 0.5;
      if( step_scale < 1.0E-12 )
        break;
    }
  }

  return best;
}// reference_minimize(...)

/** Builds a problem from explicit observed counts, so a test can construct the steeply-falling,
 partly-empty regions where the continuum's positivity constraint actually binds.
 */
ContinuumFitProblem make_explicit_problem( const vector<double> &observed,
                                           const double signal_counts,
                                           const double signal_mean_offset )
{
  const double channel_width = 2.0;
  const double first_energy = 1162.0;
  const size_t nbin = observed.size();

  ContinuumFitProblem problem;
  problem.num_coefficients = 2;
  problem.reference_energy = first_energy + 0.5 * channel_width * static_cast<double>( nbin );
  problem.observed = observed;

  for( size_t i = 0; i <= nbin; ++i )
    problem.edges.push_back( static_cast<float>( first_energy + channel_width * static_cast<double>( i ) ) );

  const double sigma = 8.0 / PhysicalUnits::fwhm_nsigma;
  const double mean = problem.reference_energy + signal_mean_offset;
  for( size_t i = 0; i < nbin; ++i )
  {
    const double lower = problem.edges[i];
    const double upper = problem.edges[i + 1];
    const double fraction = 0.5 * ( erf( ( upper - mean ) / ( sqrt( 2.0 ) * sigma ) ) -
                                    erf( ( lower - mean ) / ( sqrt( 2.0 ) * sigma ) ) );
    problem.signal.push_back( signal_counts * fraction );
  }

  return problem;
}// make_explicit_problem(...)

/** Searches for any nearby *feasible* point that improves on `coefs`.  Returns the largest
 improvement found, so zero means `coefs` really is a local minimum.

 This assumes nothing about what the right answer is, which is what makes it a sound test of "did
 the optimizer converge?" - unlike comparing against a second optimizer, which can stall in its own
 way on a constrained problem.
 */
double best_feasible_improvement( const ContinuumFitProblem &problem, const vector<double> &coefs )
{
  const double base = reference_objective( problem, coefs );
  if( !isfinite( base ) )
    return numeric_limits<double>::infinity();

  double best = 0.0;
  const double scale0 = (std::max)( fabs( coefs[0] ), 1.0E-3 );
  const double scale1 = (std::max)( fabs( coefs[1] ), 1.0E-6 );

  for( int e = -1; e >= -14; --e )
  {
    const double step = pow( 10.0, static_cast<double>( e ) );
    for( const double d0 : { -1.0, 0.0, 1.0 } )
    {
      for( const double d1 : { -1.0, 0.0, 1.0 } )
      {
        if( ( d0 == 0.0 ) && ( d1 == 0.0 ) )
          continue;
        const vector<double> candidate = { coefs[0] + d0 * step * scale0, coefs[1] + d1 * step * scale1 };
        const double value = reference_objective( problem, candidate );
        if( isfinite( value ) )
          best = (std::max)( best, base - value );
      }
    }
  }

  return best;
}// best_feasible_improvement(...)

/** Steeply-falling, partly-empty regions drive the continuum onto its positivity floor.  A Newton
 direction that leaves the feasible region does NOT mean no feasible descent direction remains -
 one can still exist along the active face - so the optimizer must not mistake "this step is
 blocked" for "converged".  Getting this wrong makes the reported limit too *small*, which is the
 dangerous direction for a detection limit.
 */
BOOST_AUTO_TEST_CASE( PoissonContinuumConvergesWhenPositivityBinds )
{
  const vector<double> falling = { 9, 8, 7, 6, 5, 4, 3, 2, 2, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  const vector<double> low_side = { 12, 10, 9, 7, 6, 5, 4, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  const vector<double> sparse = { 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

  struct Case
  {
    const char *name;
    const vector<double> *observed;
    double signal;
    double signal_offset;
  };

  const Case cases[] = {
    { "falling_no_signal", &falling, 0.0, 0.0 },
    { "falling_weak_signal", &falling, 30.0, 14.0 },
    { "low_side_strong_signal", &low_side, 200.0, 12.0 },
    { "sparse_strong_signal", &sparse, 80.0, 16.0 },
  };

  for( const Case &test : cases )
  {
    const ContinuumFitProblem problem = make_explicit_problem( *test.observed, test.signal, test.signal_offset );

    const DetectionLimitCalc::PoissonContinuumFit fit = DetectionLimitCalc::fit_continuum_poisson(
      problem.edges.data(), problem.observed.data(), problem.signal_ptr(), problem.nbin(),
      problem.num_coefficients, problem.reference_energy, {} );

    BOOST_REQUIRE_MESSAGE( fit.converged, test.name << ": " << fit.error );

    const double improvement = best_feasible_improvement( problem, fit.coefficients );

    // The scale that matters is the confidence delta itself, 2.705543: an error of that size in
    // the profiled statistic moves the reported limit by roughly its own uncertainty.  Demand far
    // better than that.
    BOOST_CHECK_MESSAGE( improvement < 0.01,
                         test.name << ": reported statistic " << fit.statistic
                                   << " can be improved by " << improvement
                                   << " with a nearby feasible point, so the fit did not converge"
                                      " (restarts=" << fit.num_restarts << ")" );
  }
}// BOOST_AUTO_TEST_CASE( PoissonContinuumConvergesWhenPositivityBinds )

/** The continuum fit runs a few hundred times per limit, so its iteration budget is part of
 whether the tool is usable.  A well-conditioned region should converge from the first starting
 estimate with plain Newton and never reach a fallback.
 */
BOOST_AUTO_TEST_CASE( PoissonContinuumIterationBudget )
{
  set_data_dir();

  cout << "DECON_FIT_BUDGET,norm,continuum_counts_per_channel,signal_counts,iterations,restarts\n";

  // Every continuum treatment is covered, not just Floating: FixedByEdges runs the optimizer a
  // second time per trial to fit its side channels, and BackgroundReference runs it a second time
  // to solve the null the expected-counts sample is built from.  Both of those are per-activity costs in a
  // scan that evaluates a few hundred activities, so they are worth watching.
  const pair<const char *, DetectionLimitCalc::DeconContinuumNorm> norms[] = {
    { "Floating", DetectionLimitCalc::DeconContinuumNorm::Floating },
    { "FixedByEdges", DetectionLimitCalc::DeconContinuumNorm::FixedByEdges },
    { "FixedByFullRange", DetectionLimitCalc::DeconContinuumNorm::FixedByFullRange },
  };

  for( const auto &norm : norms )
  {
    for( const double continuum : { 0.1, 1.0, 10.0, 100.0 } )
    {
      for( const double signal : { 0.0, 100.0 } )
      {
        DetectionLimitCalc::DeconComputeInput input = *make_decon_input( signal, continuum, norm.second );
        input.activity = signal;

        const DetectionLimitCalc::DeconComputeResults results = DetectionLimitCalc::decon_compute_peaks( input );

        cout << "DECON_FIT_BUDGET," << norm.first << ',' << continuum << ',' << signal << ','
             << results.num_continuum_iterations << ',' << results.num_continuum_restarts << '\n';

        BOOST_CHECK_MESSAGE( results.num_continuum_restarts == 0,
                             norm.first << " continuum=" << continuum << " signal=" << signal
                                          << ": needed " << results.num_continuum_restarts
                                          << " fallback attempts on a well-behaved region" );
        BOOST_CHECK_MESSAGE( results.num_continuum_iterations < 60,
                             norm.first << " continuum=" << continuum << " signal=" << signal
                                          << ": took " << results.num_continuum_iterations
                                          << " iterations" );
      }
    }
  }

  // The background-reference model scores twice as many channels and solves an extra null fit, so
  // give it its own row rather than folding it in above.
  for( const double continuum : { 1.0, 100.0 } )
  {
    DetectionLimitCalc::DeconComputeInput input
        = *make_background_reference_input( continuum, 1.0, 1.0 );
    input.activity = 0.0;

    const DetectionLimitCalc::DeconComputeResults results = DetectionLimitCalc::decon_compute_peaks( input );
    cout << "DECON_FIT_BUDGET,BackgroundReference," << continuum << ",0,"
         << results.num_continuum_iterations << ',' << results.num_continuum_restarts << '\n';

    BOOST_CHECK_MESSAGE( results.num_continuum_restarts == 0,
                         "BackgroundReference continuum=" << continuum << ": needed "
                         << results.num_continuum_restarts << " fallback attempts" );
    BOOST_CHECK_MESSAGE( results.num_continuum_iterations < 60,
                         "BackgroundReference continuum=" << continuum << ": took "
                         << results.num_continuum_iterations << " iterations" );
  }
}// BOOST_AUTO_TEST_CASE( PoissonContinuumIterationBudget )


BOOST_AUTO_TEST_CASE( StudyPoissonSamplerIsPoisson )
{
  // Every coverage number in `DeconCoverageStudy` rests on this sampler, and a wrong one would not
  //  announce itself - it would just quietly shift coverage.  Both branches are exercised: Knuth
  //  below mean 10, Hoermann's transformed rejection above it.
  const double means[] = { 0.05, 0.5, 3.0, 9.5, 10.5, 40.0, 250.0 };

  for( const double mean : means )
  {
    const size_t num_draws = 200000;
    mt19937 generator( 20260809u );

    double sum = 0.0, sum_squares = 0.0;
    map<int,size_t> histogram;
    for( size_t i = 0; i < num_draws; ++i )
    {
      const int draw = study_poisson( mean, generator );
      BOOST_REQUIRE_MESSAGE( draw >= 0, "study_poisson(" << mean << ") returned " << draw );
      sum += draw;
      sum_squares += static_cast<double>(draw) * draw;
      histogram[draw] += 1;
    }

    const double n = static_cast<double>( num_draws );
    const double sample_mean = sum / n;
    const double sample_variance = sum_squares/n - sample_mean*sample_mean;

    // A Poisson has mean == variance.  The standard error of the mean is sqrt(mean/n), so five of
    //  them is a wide berth for a correct sampler and still far tighter than any plausible bug.
    const double mean_tolerance = 5.0 * sqrt( (std::max)( mean, 1.0E-3 ) / n );
    BOOST_CHECK_MESSAGE( fabs( sample_mean - mean ) < (mean_tolerance + 1.0E-3),
                        "mean " << mean << ": sampled mean " << sample_mean
                        << " is off by more than " << mean_tolerance );

    BOOST_CHECK_MESSAGE( fabs( sample_variance - mean ) < (0.05*mean + 0.01),
                        "mean " << mean << ": sampled variance " << sample_variance
                        << " should equal the mean for a Poisson" );

    // Pearson chi-square against the exact pmf, pooling neighbouring values until each cell has a
    //  usefully large expectation.  3-sigma on a chi-square with `dof` degrees of freedom is
    //  dof + 3*sqrt(2*dof).
    double chi2 = 0.0;
    size_t num_cells = 0;
    double pmf = exp( -mean );           // p(0)
    double pooled_expected = 0.0;
    size_t pooled_observed = 0;
    for( int k = 0; k < 10000; ++k )
    {
      if( k > 0 )
        pmf *= mean / k;

      const double expected = n * pmf;
      const map<int,size_t>::const_iterator pos = histogram.find( k );
      const size_t observed = (pos == histogram.end()) ? 0 : pos->second;

      pooled_expected += expected;
      pooled_observed += observed;

      if( pooled_expected >= 10.0 )
      {
        const double residual = static_cast<double>(pooled_observed) - pooled_expected;
        chi2 += residual*residual/pooled_expected;
        ++num_cells;
        pooled_expected = 0.0;
        pooled_observed = 0;
      }

      if( (k > mean) && (expected < 1.0E-6) )
        break;
    }//for( loop over the pmf )

    BOOST_REQUIRE_MESSAGE( num_cells > 2,
                          "mean " << mean << ": only " << num_cells << " usable cells" );

    const double dof = static_cast<double>( num_cells - 1 );
    const double limit = dof + 3.0*sqrt( 2.0*dof );
    BOOST_CHECK_MESSAGE( chi2 < limit,
                        "mean " << mean << ": goodness-of-fit chi2 " << chi2 << " over "
                        << num_cells << " cells exceeds " << limit
                        << " - the sampler is not Poisson" );
  }//for( loop over means )
}// BOOST_AUTO_TEST_CASE( StudyPoissonSamplerIsPoisson )


BOOST_AUTO_TEST_CASE( CashStatisticHandCalculated )
{
  // observed == 0 leaves only the `2*expected` term, since x*log(x) -> 0.
  BOOST_CHECK_CLOSE( DetectionLimitCalc::poisson_deviance( 0.0, 1.0 ), 2.0, 1.0E-9 );
  BOOST_CHECK_CLOSE( DetectionLimitCalc::poisson_deviance( 0.0, 7.5 ), 15.0, 1.0E-9 );

  // A perfectly predicted channel contributes exactly nothing.
  BOOST_CHECK_SMALL( DetectionLimitCalc::poisson_deviance( 4.0, 4.0 ), 1.0E-12 );
  BOOST_CHECK_SMALL( DetectionLimitCalc::poisson_deviance( 1.0E6, 1.0E6 ), 1.0E-6 );

  // 2*(2 - 1 + 1*ln(1/2)) = 2 - 2*ln(2)
  BOOST_CHECK_CLOSE( DetectionLimitCalc::poisson_deviance( 1.0, 2.0 ), 2.0 - 2.0 * std::log( 2.0 ), 1.0E-9 );
  // 2*(1 - 2 + 2*ln(2/1)) = -2 + 4*ln(2)
  BOOST_CHECK_CLOSE( DetectionLimitCalc::poisson_deviance( 2.0, 1.0 ), -2.0 + 4.0 * std::log( 2.0 ), 1.0E-9 );

  // The deviance is non-negative everywhere, and at high counts approaches the Gaussian
  //  (observed-expected)^2/observed that the old modified-Neyman score used.
  const double observed = 10000.0, expected = 10100.0;
  BOOST_CHECK_CLOSE( DetectionLimitCalc::poisson_deviance( observed, expected ),
                     ( expected - observed ) * ( expected - observed ) / observed, 1.0 );

  BOOST_CHECK_THROW( DetectionLimitCalc::poisson_deviance( 1.0, 0.0 ), runtime_error );
  BOOST_CHECK_THROW( DetectionLimitCalc::poisson_deviance( -1.0, 1.0 ), runtime_error );
  BOOST_CHECK_THROW( DetectionLimitCalc::poisson_deviance( numeric_limits<double>::quiet_NaN(), 1.0 ),
                     runtime_error );
}// BOOST_AUTO_TEST_CASE( CashStatisticHandCalculated )

BOOST_AUTO_TEST_CASE( ConfidenceLevelPctStrDistinguishesSigmaChoices )
{
  using namespace DetectionLimitCalc;

  // The values the confidence-level selectors offer.
  BOOST_CHECK_EQUAL( confidence_level_pct_str( 0.95 ), "95%" );
  BOOST_CHECK_EQUAL( confidence_level_pct_str( 0.99 ), "99%" );
  BOOST_CHECK_EQUAL( confidence_level_pct_str( 0.682689492137086 ), "68.27%" );
  BOOST_CHECK_EQUAL( confidence_level_pct_str( 0.954499736103642 ), "95.45%" );
  BOOST_CHECK_EQUAL( confidence_level_pct_str( 0.997300203936740 ), "99.73%" );

  // The whole reason this function exists: the old "%.1f%%" formatting rendered both of these as
  //  "100.0%", making two different confidence levels indistinguishable on screen.
  BOOST_CHECK_EQUAL( confidence_level_pct_str( 0.999936657516334 ), "99.9937%" );
  BOOST_CHECK_EQUAL( confidence_level_pct_str( 0.999999426696856 ), "99.999943%" );
  BOOST_CHECK( confidence_level_pct_str( 0.999936657516334 )
               != confidence_level_pct_str( 0.999999426696856 ) );

  // Out-of-domain input must not produce a confident-looking percentage.
  BOOST_CHECK_EQUAL( confidence_level_pct_str( 0.0 ), "?" );
  BOOST_CHECK_EQUAL( confidence_level_pct_str( 1.0 ), "?" );
  BOOST_CHECK_EQUAL( confidence_level_pct_str( -0.5 ), "?" );
  BOOST_CHECK_EQUAL( confidence_level_pct_str( numeric_limits<double>::quiet_NaN() ), "?" );

  // `confidence_level_str` keeps its complement form - it is a documented batch-report contract.
  BOOST_CHECK_EQUAL( confidence_level_str( 0.95 ), "95%" );
  BOOST_CHECK_EQUAL( confidence_level_str( 0.999999426696856 ), "1-5.7E-07" );
}// BOOST_AUTO_TEST_CASE( ConfidenceLevelPctStrDistinguishesSigmaChoices )


BOOST_AUTO_TEST_CASE( LimitTextKindMatrix )
{
  using namespace DetectionLimitCalc;

  const DeconLimitType upper = DeconLimitType::OneSidedUpperLimit;
  const DeconLimitType central = DeconLimitType::CentralInterval;

  // ---- Activity scans: reported from the UPPER crossing. ----
  // Upper only, ordinary measurement model.
  BOOST_CHECK( decon_limit_text_kind( false, false, true, false, upper )
               == DeconLimitTextKind::OneSidedUpper );

  // Upper only, background-reference model: a prediction, never an upper limit on this spectrum.
  BOOST_CHECK( decon_limit_text_kind( false, false, true, true, upper )
               == DeconLimitTextKind::PredictedSensitivity );

  // Both crossings at a ONE-sided threshold are two one-sided bounds, not an interval.  This is the
  //  case that used to be worded "Between L and U at 95% CL", overstating the coverage.
  BOOST_CHECK( decon_limit_text_kind( false, true, true, false, upper )
               == DeconLimitTextKind::TwoOneSidedBounds );

  // Both crossings at a central threshold really are an interval.
  BOOST_CHECK( decon_limit_text_kind( false, true, true, false, central )
               == DeconLimitTextKind::CentralInterval );

  // No upper crossing bounds nothing, whether or not a lower one was found.
  BOOST_CHECK( decon_limit_text_kind( false, false, false, false, upper ) == DeconLimitTextKind::None );
  BOOST_CHECK( decon_limit_text_kind( false, true, false, false, upper ) == DeconLimitTextKind::None );

  // ---- Distance scans: reported from the LOWER crossing, so the asymmetry is inverted. ----
  BOOST_CHECK( decon_limit_text_kind( true, true, false, false, upper )
               == DeconLimitTextKind::DistanceLowerBound );
  BOOST_CHECK( decon_limit_text_kind( true, true, true, false, upper )
               == DeconLimitTextKind::TwoOneSidedBounds );
  BOOST_CHECK( decon_limit_text_kind( true, true, true, false, central )
               == DeconLimitTextKind::CentralInterval );

  // An upper crossing alone states nothing about a distance.
  BOOST_CHECK( decon_limit_text_kind( true, false, true, false, upper ) == DeconLimitTextKind::None );
  BOOST_CHECK( decon_limit_text_kind( true, false, false, false, upper ) == DeconLimitTextKind::None );

  // A predicted sensitivity forces `found_lower_cl` false, so a distance limit under that model can
  //  never be reported - which is why the UI does not offer the combination.
  BOOST_CHECK( decon_limit_text_kind( true, false, true, true, upper ) == DeconLimitTextKind::None );

  // ---- A central interval that only found its upper crossing. ----
  // This is the ORDINARY outcome of asking for a central interval on a spectrum with no signal:
  //  the lower endpoint runs into the physical boundary at zero.  It must NOT be worded as a
  //  one-sided bound at this CL - the root sits at quantile(chi2(1), CL), which is a (1+CL)/2
  //  one-sided bound, so calling it one-sided at CL understates its coverage.
  BOOST_CHECK( decon_limit_text_kind( false, false, true, false, central )
               == DeconLimitTextKind::CentralIntervalUpperBound );
  // ... and the same scan asked for as a one-sided limit is still a one-sided limit.
  BOOST_CHECK( decon_limit_text_kind( false, false, true, false, upper )
               == DeconLimitTextKind::OneSidedUpper );

  // A predicted sensitivity wins over the interval wordings if both crossings are somehow present.
  //  `foundLowerCl` is forced false upstream so this should be unreachable, but that is enforced by
  //  an assert() which is compiled out in Release - so pin the classifier's own behaviour.
  BOOST_CHECK( decon_limit_text_kind( false, true, true, true, upper )
               == DeconLimitTextKind::PredictedSensitivity );
  BOOST_CHECK( decon_limit_text_kind( false, true, true, true, central )
               == DeconLimitTextKind::PredictedSensitivity );
}// BOOST_AUTO_TEST_CASE( LimitTextKindMatrix )


BOOST_AUTO_TEST_CASE( PlannedMeasurementMatchesScaledSpectrum )
{
  using namespace DetectionLimitCalc;

  // Deliberately a non-unity dead-time fraction: live/real = 0.8.  With live == real the bug this
  //  test exists to catch - treating the planned REAL time as `sample_exposure`, which is compared
  //  against a LIVE time - is invisible.
  const float live_time = 80.0f, real_time = 100.0f;
  auto reference = make_shared<SpecUtils::Measurement>();
  const size_t num_channels = 256;
  reference->set_gamma_counts( make_shared<vector<float>>( num_channels, 10.0f ), live_time, real_time );
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( num_channels, { 0.0f, 3.0f }, {} );
  reference->set_energy_calibration( cal );

  const double planned_real_s = 250.0;

  // The invariant that makes the Currie and deconvolution paths describe the SAME measurement.
  const PlannedMeasurement backref = plan_measurement( reference, planned_real_s,
                                                       DeconMeasurementModel::BackgroundReference );
  const shared_ptr<const SpecUtils::Measurement> scaled
                  = scale_spectrum_for_dwell( reference, static_cast<float>(planned_real_s) );
  BOOST_CHECK_CLOSE( backref.sample_exposure, static_cast<double>(scaled->live_time()), 1.0E-4 );

  // ...which is a LIVE time, i.e. 250 * (80/100) == 200, NOT the 250 the control was given.
  BOOST_CHECK_CLOSE( backref.sample_exposure, 200.0, 1.0E-6 );
  BOOST_CHECK_CLOSE( backref.exposure_ratio, 200.0/80.0, 1.0E-6 );

  // The reference counts must stay real counts at their own exposure, or their counting statistics
  //  cannot propagate; the Currie path gets the projected spectrum.
  BOOST_CHECK( backref.decon == reference );
  BOOST_CHECK( backref.currie != reference );
  BOOST_CHECK_CLOSE( static_cast<double>(backref.currie->live_time()), 200.0, 1.0E-4 );

  // Current-spectrum model: one projected spectrum for both, and no separate exposure.
  const PlannedMeasurement current = plan_measurement( reference, planned_real_s,
                                                       DeconMeasurementModel::CurrentSpectrum );
  BOOST_CHECK( current.currie == current.decon );
  BOOST_CHECK( current.decon != reference );
  BOOST_CHECK_SMALL( current.sample_exposure, 1.0E-12 );
  BOOST_CHECK_CLOSE( current.exposure_ratio, 1.0, 1.0E-9 );

  // No planned time requested: the reference is handed back untouched by either model, so a
  //  cleared control cannot silently change a number.
  for( const DeconMeasurementModel model : { DeconMeasurementModel::CurrentSpectrum,
                                             DeconMeasurementModel::BackgroundReference } )
  {
    const PlannedMeasurement none = plan_measurement( reference, 0.0, model );
    BOOST_CHECK( none.currie == reference );
    BOOST_CHECK( none.decon == reference );
    BOOST_CHECK_SMALL( none.sample_exposure, 1.0E-12 );
    BOOST_CHECK_CLOSE( none.exposure_ratio, 1.0, 1.0E-9 );
  }

  // A null reference must not throw; the tools call this before a spectrum is loaded.
  const PlannedMeasurement empty = plan_measurement( nullptr, planned_real_s,
                                                     DeconMeasurementModel::BackgroundReference );
  BOOST_CHECK( !empty.currie );
  BOOST_CHECK( !empty.decon );
}// BOOST_AUTO_TEST_CASE( PlannedMeasurementMatchesScaledSpectrum )


BOOST_AUTO_TEST_CASE( PoissonContinuumGradientsMatchFiniteDifference )
{
  // The production optimizer's gradient/Hessian are not public, so check them the way that
  //  actually matters: at the reported solution the objective must be stationary, and a
  //  finite-difference Newton step from there must not improve it.
  struct Case
  {
    const char *name;
    vector<double> truth;
    double signal;
  };

  const Case cases[] = {
    { "flat", { 20.0 }, 0.0 },
    { "flat_with_signal", { 20.0 }, 300.0 },
    { "sloped", { 20.0, 0.35 }, 0.0 },
    { "sloped_with_signal", { 20.0, 0.35 }, 300.0 },
    { "quadratic", { 20.0, 0.35, 0.01 }, 150.0 },
    { "sparse", { 0.05 }, 0.0 },
    { "sparse_sloped", { 0.05, 0.001 }, 2.0 },
  };

  for( const Case &test : cases )
  {
    const ContinuumFitProblem problem =
      make_continuum_problem( 24, 2.0, 1000.0, test.truth, test.signal, 3.4 );

    const DetectionLimitCalc::PoissonContinuumFit fit = DetectionLimitCalc::fit_continuum_poisson(
      problem.edges.data(), problem.observed.data(), problem.signal_ptr(), problem.nbin(),
      problem.num_coefficients, problem.reference_energy, {} );

    BOOST_REQUIRE_MESSAGE( fit.converged, test.name << ": " << fit.error );

    const double best = reference_objective( problem, fit.coefficients );
    BOOST_CHECK_MESSAGE( fabs( best - fit.statistic ) < 1.0E-6 * (std::max)( 1.0, fabs( best ) ),
                         test.name << ": reported statistic " << fit.statistic
                                   << " does not match an independent evaluation " << best );

    // Central-difference gradient must vanish at the solution, relative to the curvature scale.
    for( size_t k = 0; k < fit.coefficients.size(); ++k )
    {
      const double magnitude = (std::max)( fabs( fit.coefficients[k] ), 1.0E-8 );
      const double delta = 1.0E-5 * magnitude;

      vector<double> up = fit.coefficients, down = fit.coefficients;
      up[k] += delta;
      down[k] -= delta;

      const double value_up = reference_objective( problem, up );
      const double value_down = reference_objective( problem, down );
      BOOST_REQUIRE( isfinite( value_up ) && isfinite( value_down ) );

      const double gradient = ( value_up - value_down ) / ( 2.0 * delta );
      const double curvature = ( value_up - 2.0 * best + value_down ) / ( delta * delta );
      BOOST_REQUIRE_MESSAGE( curvature > 0.0, test.name << ": non-positive curvature" );

      // |g|/sqrt(H) is the distance to the stationary point in units of the parameter's own
      //  standard deviation; anything under a thousandth of a sigma is converged.
      const double sigmas = fabs( gradient ) / sqrt( curvature );
      BOOST_CHECK_MESSAGE( sigmas < 1.0E-3,
                           test.name << ": coefficient " << k << " is " << sigmas
                                     << " sigma from stationary" );
    }
  }
}// BOOST_AUTO_TEST_CASE( PoissonContinuumGradientsMatchFiniteDifference )

BOOST_AUTO_TEST_CASE( PoissonContinuumMatchesReferenceOptimizer )
{
  struct Case
  {
    const char *name;
    size_t nbin;
    vector<double> truth;
    double signal;
    bool sample;
  };

  // Flat and sloped continua; very low and very high counts; no/weak/strong signal; and, through
  //  `sample`, Poisson realizations that produce genuinely empty channels.
  const Case cases[] = {
    { "flat_no_signal", 24, { 20.0 }, 0.0, false },
    { "flat_weak_signal", 24, { 20.0 }, 30.0, false },
    { "flat_strong_signal", 24, { 20.0 }, 3000.0, false },
    { "sloped_no_signal", 24, { 20.0, 0.4 }, 0.0, false },
    { "sloped_strong_signal", 24, { 20.0, 0.4 }, 3000.0, false },
    { "quadratic", 32, { 20.0, 0.4, 0.02 }, 400.0, false },
    { "high_counts", 24, { 5000.0, 3.0 }, 20000.0, false },
    { "low_counts_zero_channels", 24, { 0.05 }, 0.0, true },
    { "low_counts_sloped_zero_channels", 24, { 0.05, 0.002 }, 1.0, true },
    // The continuum density falls to ~9% of its central value at the low edge of the ROI, so the
    //  fit sits close to - but not across - the positivity boundary.
    { "near_positivity_boundary", 24, { 0.5, 0.019 }, 0.0, false },
  };

  mt19937 generator( 20260808u );

  for( const Case &test : cases )
  {
    const ContinuumFitProblem problem = make_continuum_problem(
      test.nbin, 2.0, 1000.0, test.truth, test.signal, 3.4, test.sample ? &generator : nullptr );

    const DetectionLimitCalc::PoissonContinuumFit fit = DetectionLimitCalc::fit_continuum_poisson(
      problem.edges.data(), problem.observed.data(), problem.signal_ptr(), problem.nbin(),
      problem.num_coefficients, problem.reference_energy, {} );

    BOOST_REQUIRE_MESSAGE( fit.converged, test.name << ": " << fit.error );
    BOOST_CHECK_MESSAGE( fit.covariance.size() == problem.num_coefficients * problem.num_coefficients,
                         test.name << ": no covariance returned" );

    // Start the independent minimizer somewhere deliberately unhelpful, so that agreement is not
    //  an artifact of a shared starting estimate.
    vector<double> start( problem.num_coefficients, 0.0 );
    double observed_sum = 0.0, width_sum = 0.0, signal_sum = 0.0;
    for( size_t i = 0; i < problem.nbin(); ++i )
    {
      observed_sum += problem.observed[i];
      signal_sum += problem.signal[i];
      width_sum += problem.edges[i + 1] - problem.edges[i];
    }
    start[0] = (std::max)( 1.0E-6, ( observed_sum - signal_sum ) / width_sum );

    const vector<double> reference = reference_minimize( problem, start );
    const double reference_value = reference_objective( problem, reference );

    // The production optimizer must do at least as well as the slow reference, up to the
    //  reference's own coarse convergence.
    BOOST_CHECK_MESSAGE( fit.statistic <= reference_value + 1.0E-4 * (std::max)( 1.0, fabs( reference_value ) ),
                         test.name << ": production statistic " << fit.statistic
                                   << " is worse than the reference " << reference_value );

    // The statistic the optimizer reports must match what its own coefficients actually give when
    //  evaluated independently.  This is what catches a coefficient transform bug: the internal
    //  objective can be perfect while the coefficients handed back are wrong.
    const double reference_at_production = reference_objective( problem, fit.coefficients );
    BOOST_CHECK_MESSAGE( fabs( reference_at_production - fit.statistic ) <
                           1.0E-6 * (std::max)( 1.0, fabs( fit.statistic ) ),
                         test.name << ": reported statistic " << fit.statistic
                                   << " does not match an independent evaluation of the returned"
                                      " coefficients, " << reference_at_production );

    // Comparing the two parameter vectors directly is only meaningful where the likelihood
    //  actually determines them.  At a fraction of a count per channel it does not: the
    //  coefficients trade off along a nearly flat valley, and two correct optimizers stop at
    //  different points on it with indistinguishable objectives.  Keep this only as a loose
    //  sanity bound - the inverted-equilibration bug it first caught showed 20-350 sigma.
    for( size_t k = 0; k < problem.num_coefficients; ++k )
    {
      const double variance = fit.covariance.empty()
                                ? 0.0
                                : fit.covariance[k * problem.num_coefficients + k];
      const double sigma = ( variance > 0.0 ) ? sqrt( variance ) : 0.0;
      if( sigma <= 0.0 )
        continue;

      const double difference = fabs( fit.coefficients[k] - reference[k] );
      BOOST_CHECK_MESSAGE( difference < 0.5 * sigma,
                           test.name << ": coefficient " << k << " differs from the reference by "
                                     << ( difference / sigma ) << " sigma" );
    }
  }
}// BOOST_AUTO_TEST_CASE( PoissonContinuumMatchesReferenceOptimizer )

BOOST_AUTO_TEST_CASE( PoissonContinuumConstraintAndInvalidInput )
{
  const ContinuumFitProblem problem = make_continuum_problem( 24, 2.0, 1000.0, { 20.0, 0.4 }, 0.0, 3.4 );

  const DetectionLimitCalc::PoissonContinuumFit unconstrained =
    DetectionLimitCalc::fit_continuum_poisson( problem.edges.data(), problem.observed.data(),
                                               problem.signal_ptr(), problem.nbin(),
                                               problem.num_coefficients, problem.reference_energy, {} );
  BOOST_REQUIRE( unconstrained.converged );
  BOOST_REQUIRE_EQUAL( unconstrained.covariance.size(), size_t( 4 ) );

  // A constraint centered far from the data, and tight, must drag the solution toward its center.
  const vector<double> center = { 2.0 * unconstrained.coefficients[0], 0.0 };
  const double offset_variance = 1.0E-6 * center[0] * center[0];
  const vector<double> precision = { 1.0 / offset_variance, 0.0, 0.0, 1.0E12 };

  const DetectionLimitCalc::PoissonContinuumFit constrained =
    DetectionLimitCalc::fit_continuum_poisson( problem.edges.data(), problem.observed.data(),
                                               problem.signal_ptr(), problem.nbin(),
                                               problem.num_coefficients, problem.reference_energy, {},
                                               center, precision );
  BOOST_REQUIRE_MESSAGE( constrained.converged, constrained.error );
  BOOST_CHECK_MESSAGE( constrained.coefficients[0] > 1.5 * unconstrained.coefficients[0],
                       "a tight constraint did not pull the offset toward its center: "
                         << constrained.coefficients[0] << " vs " << center[0] );
  BOOST_CHECK_MESSAGE( constrained.covariance[0] < unconstrained.covariance[0],
                       "a constraint should reduce the reported coefficient variance" );

  // Structurally invalid input throws; merely difficult input does not.
  BOOST_CHECK_THROW( DetectionLimitCalc::fit_continuum_poisson( nullptr, problem.observed.data(), nullptr,
                                                                problem.nbin(), 2, problem.reference_energy, {} ),
                     runtime_error );
  BOOST_CHECK_THROW( DetectionLimitCalc::fit_continuum_poisson( problem.edges.data(), problem.observed.data(),
                                                                nullptr, problem.nbin(), 0,
                                                                problem.reference_energy, {} ),
                     runtime_error );
  BOOST_CHECK_THROW( DetectionLimitCalc::fit_continuum_poisson( problem.edges.data(), problem.observed.data(),
                                                                nullptr, problem.nbin(), 2,
                                                                problem.reference_energy, { 1.0 } ),
                     runtime_error );
  BOOST_CHECK_THROW( DetectionLimitCalc::fit_continuum_poisson( problem.edges.data(), problem.observed.data(),
                                                                nullptr, problem.nbin(), 2,
                                                                problem.reference_energy, {}, center, { 1.0 } ),
                     runtime_error );

  // An all-empty ROI is degenerate but must still return a usable, positive continuum rather than
  //  aborting - this is the low-count regime the whole exercise is about.
  ContinuumFitProblem empty = problem;
  std::fill( begin( empty.observed ), end( empty.observed ), 0.0 );
  const DetectionLimitCalc::PoissonContinuumFit empty_fit =
    DetectionLimitCalc::fit_continuum_poisson( empty.edges.data(), empty.observed.data(), nullptr,
                                               empty.nbin(), empty.num_coefficients,
                                               empty.reference_energy, {} );
  BOOST_CHECK_MESSAGE( empty_fit.converged, "all-empty ROI: " << empty_fit.error );
  BOOST_CHECK_MESSAGE( empty_fit.statistic >= 0.0, "all-empty ROI gave a negative deviance" );
}// BOOST_AUTO_TEST_CASE( PoissonContinuumConstraintAndInvalidInput )

BOOST_AUTO_TEST_CASE( PoissonChannelOverloadMatchesContiguous )
{
  set_data_dir();

  // The two entry points must describe the same problem, or the sideband and background-reference
  // paths would silently be solving something different from the ordinary one.
  const size_t nbin = 24;
  const double reference_energy = 1173.0;
  vector<float> edges( nbin + 1 );
  for( size_t i = 0; i <= nbin; ++i )
    edges[i] = static_cast<float>( 1150.0 + 2.0*i );

  mt19937 generator( 20260808u );
  vector<double> observed( nbin, 0.0 ), signal( nbin, 0.0 );
  for( size_t i = 0; i < nbin; ++i )
  {
    // Steeply falling and partly empty: the regime where positivity actually binds.  A noiseless
    // fixture would let the deviance self-barrier, and the constrained path would never run.
    const double mean = 40.0*std::exp( -0.35*static_cast<double>(i) );
    observed[i] = poisson_distribution<int>( mean )( generator );
    signal[i] = 0.05*mean;
  }

  const DetectionLimitCalc::PoissonContinuumFit array_fit
      = DetectionLimitCalc::fit_continuum_poisson( edges.data(), observed.data(), signal.data(),
                                                  nbin, 2, reference_energy, {} );

  vector<DetectionLimitCalc::PoissonChannel> channels( nbin );
  for( size_t i = 0; i < nbin; ++i )
  {
    channels[i].lower_energy = edges[i];
    channels[i].upper_energy = edges[i+1];
    channels[i].observed = observed[i];
    channels[i].fixed_signal = signal[i];
    channels[i].continuum_scale = 1.0;
  }

  const DetectionLimitCalc::PoissonContinuumFit channel_fit
      = DetectionLimitCalc::fit_continuum_poisson( channels.data(), nbin, 2, reference_energy, {} );

  BOOST_REQUIRE_MESSAGE( array_fit.converged, array_fit.error );
  BOOST_REQUIRE_MESSAGE( channel_fit.converged, channel_fit.error );
  BOOST_REQUIRE_EQUAL( array_fit.coefficients.size(), channel_fit.coefficients.size() );
  for( size_t k = 0; k < array_fit.coefficients.size(); ++k )
    BOOST_CHECK_SMALL( fabs( array_fit.coefficients[k] - channel_fit.coefficients[k] ), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs( array_fit.statistic - channel_fit.statistic ), 1.0E-12 );

  // The reported information must really be the inverse of the reported covariance, since
  // FixedByEdges hands it straight back in as a precision matrix.
  BOOST_REQUIRE_EQUAL( channel_fit.information.size(), 4u );
  BOOST_REQUIRE_EQUAL( channel_fit.covariance.size(), 4u );
  for( size_t r = 0; r < 2; ++r )
  {
    for( size_t c = 0; c < 2; ++c )
    {
      double product = 0.0;
      for( size_t k = 0; k < 2; ++k )
        product += channel_fit.information[r*2 + k] * channel_fit.covariance[k*2 + c];
      BOOST_CHECK_SMALL( fabs( product - ((r == c) ? 1.0 : 0.0) ), 1.0E-6 );
    }
  }

  // Disjoint channels: dropping the middle of a noiseless linear continuum must still recover it.
  // A contiguous-only optimizer cannot even be asked this question, which is why the overload
  // exists - the sidebands of an ROI are two blocks with a gap between them.
  vector<DetectionLimitCalc::PoissonChannel> disjoint;
  for( size_t i = 0; i < nbin; ++i )
  {
    if( (i >= 8) && (i < 16) )
      continue;

    DetectionLimitCalc::PoissonChannel channel;
    channel.lower_energy = edges[i];
    channel.upper_energy = edges[i+1];
    channel.continuum_scale = 1.0;
    channel.fixed_signal = 0.0;

    // Exactly the integral of (30 - 0.5*(E - ref)) over the channel.
    const double x0 = channel.lower_energy - reference_energy;
    const double x1 = channel.upper_energy - reference_energy;
    channel.observed = 30.0*(x1 - x0) - 0.25*(x1*x1 - x0*x0);
    disjoint.push_back( channel );
  }

  const DetectionLimitCalc::PoissonContinuumFit disjoint_fit
      = DetectionLimitCalc::fit_continuum_poisson( disjoint.data(), disjoint.size(), 2,
                                                  reference_energy, {} );
  BOOST_REQUIRE_MESSAGE( disjoint_fit.converged, disjoint_fit.error );
  BOOST_CHECK_CLOSE( disjoint_fit.coefficients[0], 30.0, 0.01 );
  BOOST_CHECK_CLOSE( disjoint_fit.coefficients[1], -0.5, 0.01 );
}// BOOST_AUTO_TEST_CASE( PoissonChannelOverloadMatchesContiguous )


BOOST_AUTO_TEST_CASE( PoissonInformationAndConstraintMagnitude )
{
  set_data_dir();

  // `information` is a reported diagnostic and the natural input to any future constrained
  // nuisance term, so its MAGNITUDE needs pinning.  Note that checking
  // `information * covariance == I` would be vacuous: `covariance` is computed as
  // `pinv(information)`, so that identity holds for ANY rescaling.  These two checks pin it
  // against arithmetic done by hand instead.

  const double reference_energy = 1000.0;
  const double channel_width = 2.0;
  const double counts_per_channel = 100.0;
  const size_t nbin = 10;

  vector<DetectionLimitCalc::PoissonChannel> channels( nbin );
  for( size_t i = 0; i < nbin; ++i )
  {
    channels[i].lower_energy = reference_energy + channel_width*static_cast<double>(i);
    channels[i].upper_energy = channels[i].lower_energy + channel_width;
    channels[i].observed = counts_per_channel;
    channels[i].fixed_signal = 0.0;
    channels[i].continuum_scale = 1.0;
  }

  const DetectionLimitCalc::PoissonContinuumFit fit
      = DetectionLimitCalc::fit_continuum_poisson( channels.data(), nbin, 2, reference_energy, {} );
  BOOST_REQUIRE_MESSAGE( fit.converged, fit.error );
  BOOST_REQUIRE_EQUAL( fit.information.size(), 4u );

  // Flat noiseless data, so the fit recovers `counts_per_channel` per channel exactly and every
  // expected count equals it.  The Fisher information is `F_kl = sum_i I_ik*I_il/E_i` with
  // `I_i0` the channel width, so  F_00 = nbin * width^2 / counts.
  const double expected_f00 = nbin * channel_width * channel_width / counts_per_channel;
  BOOST_CHECK_MESSAGE( fabs( fit.information[0] - expected_f00 ) < 0.01*expected_f00,
                      "Fisher information F_00 is " << fit.information[0] << ", hand calculation"
                      " gives " << expected_f00 << ".  A wrong scale here silently mis-calibrates"
                      " every FixedByEdges limit." );

  // Now pin how the optimizer USES a supplied precision matrix.  At the constrained solution the
  // gradient of `deviance(c) + (c-mu)^T P (c-mu)` vanishes, so the channel-only gradient must
  // equal `-2*P*(c_hat - mu)`.  This catches a factor of two, a sign error, or the constraint
  // being ignored - independently of the optimizer's internals, since the test differentiates the
  // channel deviance itself.
  const vector<double> center{ 0.9*counts_per_channel/channel_width, 0.0 };
  const vector<double> precision{ 5.0, 0.0, 0.0, 500.0 };

  const DetectionLimitCalc::PoissonContinuumFit constrained
      = DetectionLimitCalc::fit_continuum_poisson( channels.data(), nbin, 2, reference_energy, {},
                                                  center, precision );
  BOOST_REQUIRE_MESSAGE( constrained.converged, constrained.error );

  const auto channel_deviance = [&]( const vector<double> &coefs ) -> double {
    double total = 0.0;
    for( size_t i = 0; i < nbin; ++i )
    {
      const double x0 = channels[i].lower_energy - reference_energy;
      const double x1 = channels[i].upper_energy - reference_energy;
      const double expected = coefs[0]*(x1 - x0) + coefs[1]*0.5*(x1*x1 - x0*x0);
      total += DetectionLimitCalc::poisson_deviance( channels[i].observed, expected );
    }
    return total;
  };

  for( size_t k = 0; k < 2; ++k )
  {
    vector<double> up = constrained.coefficients, down = constrained.coefficients;
    const double step = 1.0E-5 * (std::max)( 1.0, fabs(constrained.coefficients[k]) );
    up[k] += step;
    down[k] -= step;
    const double numeric_gradient = (channel_deviance(up) - channel_deviance(down)) / (2.0*step);

    double penalty_gradient = 0.0;
    for( size_t c = 0; c < 2; ++c )
      penalty_gradient += 2.0 * precision[k*2 + c] * (constrained.coefficients[c] - center[c]);

    const double scale = (std::max)( 1.0, fabs(numeric_gradient) );
    BOOST_CHECK_MESSAGE( fabs( numeric_gradient + penalty_gradient ) < 1.0E-3*scale,
                        "coefficient " << k << ": channel-deviance gradient " << numeric_gradient
                        << " does not cancel the constraint gradient " << penalty_gradient
                        << "; the penalty is not entering the objective as (c-mu)^T P (c-mu)" );
  }
}// BOOST_AUTO_TEST_CASE( PoissonInformationAndConstraintMagnitude )


BOOST_AUTO_TEST_CASE( PoissonChannelScaleIsApplied )
{
  set_data_dir();

  // `continuum_scale` is the reason the channel overload exists, and every other test leaves it at
  // 1.  Two blocks over the same energies, one counted `scale` times longer, must recover the same
  // per-unit-exposure continuum - exactly what the background-reference model relies on.
  const double reference_energy = 1000.0;
  const double scale = 25.0;
  const size_t nbin = 12;

  vector<DetectionLimitCalc::PoissonChannel> channels;
  for( size_t block = 0; block < 2; ++block )
  {
    const double block_scale = block ? scale : 1.0;
    for( size_t i = 0; i < nbin; ++i )
    {
      DetectionLimitCalc::PoissonChannel channel;
      channel.lower_energy = reference_energy + 2.0*static_cast<double>(i);
      channel.upper_energy = channel.lower_energy + 2.0;
      channel.fixed_signal = 0.0;
      channel.continuum_scale = block_scale;

      // Exactly `block_scale` times the integral of (40 - 0.3*(E - ref)) over the channel.
      const double x0 = channel.lower_energy - reference_energy;
      const double x1 = channel.upper_energy - reference_energy;
      channel.observed = block_scale * ( 40.0*(x1 - x0) - 0.15*(x1*x1 - x0*x0) );
      channels.push_back( channel );
    }
  }

  const DetectionLimitCalc::PoissonContinuumFit fit
      = DetectionLimitCalc::fit_continuum_poisson( channels.data(), channels.size(), 2,
                                                  reference_energy, {} );
  BOOST_REQUIRE_MESSAGE( fit.converged, fit.error );

  // Ignoring `continuum_scale` would fit a compromise between the two blocks rather than
  // recovering the generating continuum.
  BOOST_CHECK_CLOSE( fit.coefficients[0], 40.0, 0.1 );
  BOOST_CHECK_CLOSE( fit.coefficients[1], -0.3, 0.1 );
}// BOOST_AUTO_TEST_CASE( PoissonChannelScaleIsApplied )


BOOST_AUTO_TEST_CASE( SidebandConstraintUsesFullCovariance )
{
  set_data_dir();

  // `fit_continuum_poisson`'s Gaussian-constraint facility is no longer used by FixedByEdges -
  // that now folds the side channels into the likelihood directly - but the facility itself
  // remains, and a constrained nuisance term is what a systematic-uncertainty model will need.
  // Check it is genuinely multivariate: zeroing the offset/slope off-diagonal must move the
  // answer, or only two independent scalar penalties are reaching the objective.
  const size_t nbin = 20;
  const double reference_energy = 1173.0;
  mt19937 generator( 20260808u );

  vector<DetectionLimitCalc::PoissonChannel> channels( nbin );
  for( size_t i = 0; i < nbin; ++i )
  {
    channels[i].lower_energy = 1150.0 + 2.0*i;
    channels[i].upper_energy = 1152.0 + 2.0*i;
    channels[i].observed = poisson_distribution<int>( 60.0 )( generator );
    channels[i].fixed_signal = 0.0;
    channels[i].continuum_scale = 1.0;
  }

  const vector<double> center{ 30.0, -0.4 };
  const vector<double> full{ 4.0, 1.8, 1.8, 4.0 };      // deliberately strongly correlated
  const vector<double> diagonal{ 4.0, 0.0, 0.0, 4.0 };

  const DetectionLimitCalc::PoissonContinuumFit full_fit
      = DetectionLimitCalc::fit_continuum_poisson( channels.data(), nbin, 2, reference_energy, {},
                                                  center, full );
  const DetectionLimitCalc::PoissonContinuumFit diagonal_fit
      = DetectionLimitCalc::fit_continuum_poisson( channels.data(), nbin, 2, reference_energy, {},
                                                  center, diagonal );

  BOOST_REQUIRE_MESSAGE( full_fit.converged, full_fit.error );
  BOOST_REQUIRE_MESSAGE( diagonal_fit.converged, diagonal_fit.error );

  const double coefficient_shift = fabs( full_fit.coefficients[0] - diagonal_fit.coefficients[0] )
                                   + fabs( full_fit.coefficients[1] - diagonal_fit.coefficients[1] );
  BOOST_CHECK_MESSAGE( coefficient_shift > 1.0E-6,
                      "zeroing the offset/slope correlation left the fit unchanged, so the"
                      " off-diagonal is not reaching the objective" );
}// BOOST_AUTO_TEST_CASE( SidebandConstraintUsesFullCovariance )


BOOST_AUTO_TEST_CASE( SidebandConstraintWidensLimit )
{
  set_data_dir();

  // FixedByEdges scores the side channels as ordinary Poisson observations alongside the region,
  // so fewer of them means a weaker anchor, the continuum is freer to absorb signal, and the limit
  // must widen.  The original frozen-line implementation propagated none of this and was therefore
  // insensitive to the side-channel count - exactly why it under-covered.
  const auto limit_for = []( const size_t num_side ) -> double {
    shared_ptr<DetectionLimitCalc::DeconComputeInput> input
        = make_shared<DetectionLimitCalc::DeconComputeInput>(
            *make_decon_input( 0.0, 100.0, DetectionLimitCalc::DeconContinuumNorm::FixedByEdges ) );
    input->roi_info.front().num_lower_side_channels = num_side;
    input->roi_info.front().num_upper_side_channels = num_side;

    ScopedCoutSilence silence;
    const DetectionLimitCalc::DeconActivityOrDistanceLimitResult result
        = DetectionLimitCalc::get_activity_or_distance_limits( 0.95f, input, false, 0.0, 800.0, true );
    BOOST_REQUIRE_MESSAGE( result.foundUpperCl, "no limit for " << num_side << " side channels: "
                          << result.errorMessage );
    return result.upperLimit;
  };

  const double few = limit_for( 2 );
  const double some = limit_for( 8 );
  const double many = limit_for( 64 );

  BOOST_CHECK_MESSAGE( few > some,
                      "2 side channels gave " << few << ", 8 gave " << some
                      << "; a weaker anchor must give a wider limit" );
  BOOST_CHECK_MESSAGE( some > many,
                      "8 side channels gave " << some << ", 64 gave " << many );
}// BOOST_AUTO_TEST_CASE( SidebandConstraintWidensLimit )


/** The production projection Monte Carlo: the properties a caller relies on.

 The statistical validation - that the median tracks and the band covers - is
 `ProjectedMeasurementMonteCarlo`, which is opt-in because it costs minutes.  This pins the
 contract that has to hold on every run.
 */
BOOST_AUTO_TEST_CASE( ProjectedLimitApi )
{
  using namespace DetectionLimitCalc;

  set_data_dir();

  const size_t num_channels = 1024;
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( num_channels, { 0.0f, 2.0f }, {} );

  auto counts = make_shared<vector<float>>( num_channels, 20.0f );
  auto reference = make_shared<SpecUtils::Measurement>();
  reference->set_gamma_counts( counts, 100.0f, 125.0f );   // real dead time, live != real
  reference->set_energy_calibration( cal );

  const size_t first_channel = 480, last_channel = 520;

  // --- draw_projected_measurement -------------------------------------------------------------
  {
    mt19937 generator( 12345 );
    const shared_ptr<const SpecUtils::Measurement> drawn
        = draw_projected_measurement( reference, 400.0f, first_channel, last_channel, generator );

    BOOST_REQUIRE( drawn );

    // The dwell is entered as a REAL time and the counts scale by the ratio of REAL times, the same
    //  convention `scale_spectrum_for_dwell` uses - not the ratio of live times.  The reference here
    //  has genuine dead time (100 s live in 125 s real) precisely so the two cannot be confused.
    const double ratio = 400.0/125.0;
    BOOST_CHECK_CLOSE_FRACTION( drawn->real_time(), 400.0f, 1.0E-5 );
    BOOST_CHECK_CLOSE_FRACTION( drawn->live_time(), static_cast<float>(100.0*ratio), 1.0E-5 );
    BOOST_CHECK_EQUAL( drawn->num_gamma_channels(), num_channels );
    BOOST_REQUIRE( drawn->energy_calibration() );
    BOOST_CHECK( (*drawn->energy_calibration()) == (*cal) );

    // Outside the window the expectation is written straight in.
    const double expected = ratio*20.0;
    BOOST_CHECK_CLOSE_FRACTION( drawn->gamma_channel_content( 10 ),
                               static_cast<float>(expected), 1.0E-4 );
    BOOST_CHECK_CLOSE_FRACTION( drawn->gamma_channel_content( 900 ),
                               static_cast<float>(expected), 1.0E-4 );

    // Inside it, the counts are drawn - so they must not all sit exactly on the expectation.
    size_t num_differing = 0;
    for( size_t i = first_channel; i <= last_channel; ++i )
      num_differing += ( fabs( drawn->gamma_channel_content(i) - expected ) > 1.0E-3 );
    BOOST_CHECK_MESSAGE( num_differing > (last_channel - first_channel)/2,
                        "only " << num_differing << " channels in the window were actually drawn" );

    // A channel that recorded zero must be able to come out non-zero - the whole reason the rate
    //  is drawn from a Gamma rather than the count from a Poisson.  \sa projected_limit_prior_strength
    auto zero_counts = make_shared<vector<float>>( num_channels, 0.0f );
    auto zero_ref = make_shared<SpecUtils::Measurement>();
    zero_ref->set_gamma_counts( zero_counts, 100.0f, 100.0f );
    zero_ref->set_energy_calibration( cal );

    size_t num_nonzero = 0;
    for( size_t trial = 0; trial < 40; ++trial )
    {
      const shared_ptr<const SpecUtils::Measurement> z
          = draw_projected_measurement( zero_ref, 1000.0f, first_channel, last_channel, generator );
      for( size_t i = first_channel; i <= last_channel; ++i )
        num_nonzero += ( z->gamma_channel_content(i) > 0.0f );
    }
    BOOST_CHECK_MESSAGE( num_nonzero > 0,
                        "an all-zero reference produced no counts in any of 40 realisations,"
                        " which means the rate draw is locked at zero" );

    BOOST_CHECK_THROW( draw_projected_measurement( reference, -1.0f, first_channel, last_channel,
                                                  generator ), std::exception );
    BOOST_CHECK_THROW( draw_projected_measurement( reference, 400.0f, 10, num_channels + 5,
                                                  generator ), std::exception );
    BOOST_CHECK_THROW( draw_projected_measurement( nullptr, 400.0f, 0, 10, generator ),
                      std::exception );
  }

  // --- currie_projected_limit -----------------------------------------------------------------
  CurrieMdaInput input;
  input.spectrum = reference;
  input.gamma_energy = 1000.0f;
  input.roi_lower_energy = 988.0f;
  input.roi_upper_energy = 1012.0f;
  input.num_lower_side_channels = 4;
  input.num_upper_side_channels = 4;
  input.detection_probability = 0.95;
  input.additional_uncertainty = 0.0f;

  const ProjectedLimit projected = currie_projected_limit( input, 400.0, 400 );
  BOOST_REQUIRE_MESSAGE( projected.valid, "no projected Currie limit" );
  BOOST_CHECK_EQUAL( projected.num_attempted, 400 );
  BOOST_CHECK_GT( projected.num_used, 350 );

  // A band, in the right order, and actually a band rather than a point.
  BOOST_CHECK_LE( projected.lower, projected.median );
  BOOST_CHECK_LE( projected.median, projected.upper );
  BOOST_CHECK_MESSAGE( (projected.upper - projected.lower) > 0.02*projected.median,
                      "the 68% band is " << (projected.upper - projected.lower)
                      << " wide about a median of " << projected.median );

  // Deterministic: the same input must give the same number, or a displayed MDA jitters for no
  //  reason the user can see.
  const ProjectedLimit repeat = currie_projected_limit( input, 400.0, 400 );
  BOOST_CHECK_EQUAL( repeat.median, projected.median );
  BOOST_CHECK_EQUAL( repeat.lower, projected.lower );
  BOOST_CHECK_EQUAL( repeat.upper, projected.upper );

  // Four times the dwell must give a smaller limit *rate*; in counts it grows about as its root.
  const ProjectedLimit longer = currie_projected_limit( input, 1600.0, 400 );
  BOOST_REQUIRE( longer.valid );
  BOOST_CHECK_MESSAGE( longer.median > projected.median,
                      "a longer dwell gave fewer limit counts: " << longer.median
                      << " against " << projected.median );
  BOOST_CHECK_MESSAGE( (longer.median/projected.median) < 4.0,
                      "limit counts grew by " << (longer.median/projected.median)
                      << "x for a 4x dwell, which is faster than counting statistics allow" );

  // What the band has to do as the dwell grows, and it is not "get wider".
  //
  //  The detection limit goes about as the root of the region counts B, so its relative spread is
  //  about half B's.  Under the predictive draw `Var(B) = k(n+alpha)(1+k)` about a mean of `k*n`,
  //  so B's relative spread is `sqrt((1+k)/(k*n))` and the band goes as `sqrt((1+k)/k)`.  That
  //  *falls* with k, towards the floor `sqrt(1/n)` set by how well the reference pins the rate:
  //  a longer projection is built on more counts and so is better determined, relatively.
  //
  //  A plain scaling instead claims `Var(B) = k*n`, giving `sqrt(1/(k*n))` - the same floor-free
  //  `1/sqrt(k)` decay all the way down.  So the test is that the band decays *slower* than that:
  //  quadrupling the dwell must not come close to halving it.  Here k goes 3.2 -> 12.8, for a
  //  predicted ratio of sqrt(1.078/1.3125) = 0.91 against the 0.50 a plain scaling implies.
  const double near_band = (projected.upper - projected.lower)/projected.median;
  const double far_band = (longer.upper - longer.lower)/longer.median;
  BOOST_CHECK_MESSAGE( far_band < near_band,
                      "relative band did not tighten with dwell: " << near_band << " -> "
                      << far_band );
  BOOST_CHECK_MESSAGE( far_band > 0.75*near_band,
                      "relative band fell from " << near_band << " to " << far_band << " for a 4x"
                      " dwell - close to the 0.5 of a plain scaling, so the reference's own rate"
                      " uncertainty is not reaching the band" );

  BOOST_CHECK( !currie_projected_limit( input, -1.0, 400 ).valid );
  BOOST_CHECK( !currie_projected_limit( input, 400.0, 2 ).valid );
}// BOOST_AUTO_TEST_CASE( ProjectedLimitApi )


/** The threaded deconvolution projection: the answer must not depend on the thread count.

 The realisations are drawn up front, single threaded, and written back by index precisely so that
 a machine with four cores and a machine with one produce the same limit.  A Monte Carlo whose
 result depended on how the work happened to be split would be indefensible in a report.
 */
BOOST_AUTO_TEST_CASE( DeconProjectedLimitThreading )
{
  using namespace DetectionLimitCalc;

  set_data_dir();

  // 100 counts/channel of flat continuum over a 1 s reference, and no signal.
  const shared_ptr<const DeconComputeInput> base = make_decon_input( 0.0, 100.0 );
  BOOST_REQUIRE( base && base->measurement );

  const size_t trials = 48;

  ProjectedLimit one, four;
  {
    ScopedCoutSilence silence;
    one = decon_projected_limit( base, 0.95, 4.0, 5000.0, true,
                                DeconLimitType::OneSidedUpperLimit,
                                ProjectedLimitScoring::SampleOnly, trials, 1 );
    four = decon_projected_limit( base, 0.95, 4.0, 5000.0, true,
                                 DeconLimitType::OneSidedUpperLimit,
                                 ProjectedLimitScoring::SampleOnly, trials, 4 );
  }

  BOOST_REQUIRE_MESSAGE( one.valid, "single-threaded projection gave no limit" );
  BOOST_REQUIRE_MESSAGE( four.valid, "four-threaded projection gave no limit" );

  // The whole point: identical, not merely close.
  BOOST_CHECK_EQUAL( one.median, four.median );
  BOOST_CHECK_EQUAL( one.lower, four.lower );
  BOOST_CHECK_EQUAL( one.upper, four.upper );
  BOOST_CHECK_EQUAL( one.num_used, four.num_used );

  BOOST_CHECK_LE( one.lower, one.median );
  BOOST_CHECK_LE( one.median, one.upper );
  BOOST_CHECK_MESSAGE( (one.upper - one.lower) > 0.0,
                      "the projected deconvolution limit has no spread at all" );
  BOOST_CHECK_GT( one.num_used, trials/2 );

  BOOST_CHECK( !decon_projected_limit( base, 0.95, -1.0, 5000.0, true,
                                      DeconLimitType::OneSidedUpperLimit,
                                      ProjectedLimitScoring::SampleOnly, trials, 4 ).valid );
  BOOST_CHECK( !decon_projected_limit( nullptr, 0.95, 4.0, 5000.0, true,
                                      DeconLimitType::OneSidedUpperLimit,
                                      ProjectedLimitScoring::SampleOnly, trials, 4 ).valid );

  // Joint scoring only means something against an actual background reference, so asking for it on
  //  a current-spectrum input is refused rather than silently answered with the wrong estimand.
  BOOST_CHECK( !decon_projected_limit( base, 0.95, 4.0, 5000.0, true,
                                      DeconLimitType::OneSidedUpperLimit,
                                      ProjectedLimitScoring::JointWithReference, trials, 4 ).valid );

  // --- joint scoring, against a real background reference ------------------------------------
  const shared_ptr<const DeconComputeInput> backref
      = make_background_reference_input( 100.0, 4.0, 4.0 );
  BOOST_REQUIRE( backref );

  ProjectedLimit joint_one, joint_four;
  {
    ScopedCoutSilence silence;
    joint_one = decon_projected_limit( backref, 0.95, 8.0, 5000.0, true,
                                      DeconLimitType::OneSidedUpperLimit,
                                      ProjectedLimitScoring::JointWithReference, trials, 1 );
    joint_four = decon_projected_limit( backref, 0.95, 8.0, 5000.0, true,
                                       DeconLimitType::OneSidedUpperLimit,
                                       ProjectedLimitScoring::JointWithReference, trials, 4 );
  }

  BOOST_REQUIRE_MESSAGE( joint_one.valid, "joint-scored projection gave no limit" );
  BOOST_CHECK_EQUAL( joint_one.median, joint_four.median );
  BOOST_CHECK_EQUAL( joint_one.lower, joint_four.lower );
  BOOST_CHECK_EQUAL( joint_one.upper, joint_four.upper );

  BOOST_CHECK_LE( joint_one.lower, joint_one.median );
  BOOST_CHECK_LE( joint_one.median, joint_one.upper );
  BOOST_CHECK_MESSAGE( (joint_one.upper - joint_one.lower) > 0.0,
                      "the joint-scored projection has no spread at all" );

  // The point of the mode: a spread where `BackgroundReference` alone reports a single number.
  ProjectedLimit sample_only;
  {
    ScopedCoutSilence silence;
    sample_only = decon_projected_limit( backref, 0.95, 8.0, 5000.0, true,
                                        DeconLimitType::OneSidedUpperLimit,
                                        ProjectedLimitScoring::SampleOnly, trials, 4 );
  }
  BOOST_REQUIRE( sample_only.valid );

  // Scoring with the reference can only add information, so the joint median must not be the
  //  looser of the two by any margin worth speaking of.
  BOOST_CHECK_MESSAGE( joint_one.median < 1.15*sample_only.median,
                      "joint median " << joint_one.median << " against sample-only "
                      << sample_only.median << "; more data cannot loosen the limit" );
}// BOOST_AUTO_TEST_CASE( DeconProjectedLimitThreading )


/** Does the *joint* projection cover?  The one that would give `BackgroundReference` a spread.

 `ProjectedMeasurementMonteCarlo` validates `ProjectedLimitScoring::SampleOnly` - each realisation
 analysed on its own.  `BackgroundReference` reports something else: the reference and the sample
 scored together, sharing one continuum, which is legitimately tighter.  Its coverage does not
 follow from the sample-only result and has to be measured separately, which is what this does.

 Unlike the sweep, this calls the production `decon_projected_limit` rather than re-drawing here, so
 what is validated is the function a caller would actually use.
 */
BOOST_AUTO_TEST_CASE( JointProjectedLimitCoverage )
{
  using namespace DetectionLimitCalc;

  set_data_dir();

  const char * const enabled = getenv( "INTERSPEC_RUN_DECON_MC_PROJECTION" );
  if( !enabled || ( string( enabled ) != "1" ) )
  {
    BOOST_TEST_MESSAGE( "Set INTERSPEC_RUN_DECON_MC_PROJECTION=1 to run the joint projection"
                       " coverage check." );
    return;
  }

  const string path = SpecUtils::append_path( g_test_file_dir,
                                             "FitPeaksForSource/trinitite_sample_b_background.n42" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( path ), "Background spectrum not at '" << path << "'" );

  SpecMeas file;
  BOOST_REQUIRE( file.load_N42_file( path ) );

  shared_ptr<const SpecUtils::Measurement> truth_rate;
  for( const shared_ptr<const SpecUtils::Measurement> &candidate : file.measurements() )
  {
    if( candidate && (candidate->num_gamma_channels() > 1024) && (candidate->live_time() > 0.0f) )
    {
      truth_rate = candidate;
      break;
    }
  }
  BOOST_REQUIRE( truth_rate );

  const double source_live = truth_rate->live_time();
  const double live_over_real = truth_rate->live_time() / truth_rate->real_time();

  auto drf = make_shared<DetectorPeakResponse>();
  drf->setIntrinsicEfficiencyFormula( "1.0", 2.54f*PhysicalUnits::cm, PhysicalUnits::keV,
                                     0.0f, 0.0f,
                                     DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct );
  drf->setFwhmCoefficients( vector<float>{ 4.0f, 0.02f },
                           DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial );

  const float peak_energy = 1000.0f;
  const float fwhm = static_cast<float>( drf->peakResolutionFWHM( peak_energy ) );

  const auto draw_truth = [&]( const double live_seconds, mt19937 &generator )
        -> shared_ptr<const SpecUtils::Measurement> {
    const shared_ptr<const vector<float>> rate = truth_rate->gamma_counts();
    auto counts = make_shared<vector<float>>( rate->size(), 0.0f );
    const double scale = live_seconds / source_live;
    for( size_t i = 0; i < rate->size(); ++i )
      counts->at(i) = static_cast<float>( study_poisson( (std::max)(0.0, scale*rate->at(i)),
                                                        generator ) );
    auto m = make_shared<SpecUtils::Measurement>();
    m->set_gamma_counts( counts, static_cast<float>(live_seconds),
                        static_cast<float>(live_seconds/live_over_real) );
    m->set_energy_calibration( truth_rate->energy_calibration() );
    return m;
  };

  const auto make_backref_input = [&]( const shared_ptr<const SpecUtils::Measurement> &reference,
                                       const double sample_live )
        -> shared_ptr<DeconComputeInput> {
    DeconRoiInfo roi;
    roi.roi_start = peak_energy - 1.25f*fwhm;
    roi.roi_end = peak_energy + 1.25f*fwhm;
    roi.continuum_type = PeakContinuum::OffsetType::Linear;
    roi.cont_norm_method = DeconContinuumNorm::Floating;
    roi.num_lower_side_channels = 4;
    roi.num_upper_side_channels = 4;

    DeconRoiInfo::PeakInfo peak;
    peak.energy = peak_energy;
    peak.fwhm = fwhm;
    peak.counts_per_bq_into_4pi = reference->live_time();
    roi.peak_infos.push_back( peak );

    auto input = make_shared<DeconComputeInput>();
    input->activity = 0.0;
    input->distance = 0.0;
    input->include_air_attenuation = false;
    input->shielding_thickness = 0.0;
    input->drf = drf;
    input->measurement = reference;
    input->roi_info.push_back( roi );
    input->measurement_model = DeconMeasurementModel::BackgroundReference;
    input->sample_exposure = sample_live;
    return input;
  };

  cout << "DECON_JOINT_MC_COVERAGE,k,pairs,mc_trials,median_pred_over_actual,cover68\n";

  const double reference_live = 600.0;
  const size_t pairs = 40;
  const size_t mc_trials = 60;

  struct Row { double k, ratio, cover68; };
  vector<Row> rows;

  for( const double k : { 0.25, 1.0, 4.0, 16.0 } )
  {
    const double sample_live = k*reference_live;
    mt19937 generator( study_cell_seed( "jointmc|k=" + study_number_token( k ) ) );

    vector<double> ratios;
    size_t inside68 = 0, counted = 0;

    for( size_t pair = 0; pair < pairs; ++pair )
    {
      const shared_ptr<const SpecUtils::Measurement> reference
            = draw_truth( reference_live, generator );
      const shared_ptr<const SpecUtils::Measurement> future = draw_truth( sample_live, generator );

      // What the joint analysis actually returns once the future measurement is in hand.
      double actual = numeric_limits<double>::quiet_NaN();
      try
      {
        ScopedCoutSilence silence;
        const shared_ptr<DeconComputeInput> observed = make_backref_input( reference, sample_live );
        observed->observed_sample = future;
        const DeconActivityOrDistanceLimitResult r =
            get_activity_or_distance_limits( 0.95f, observed, false, 0.0, 5.0, true );
        if( r.foundUpperCl )
          actual = r.upperLimit;
      }catch( std::exception & ){}

      if( IsNan(actual) )
        continue;

      // What the projection promised, from the reference alone.
      ProjectedLimit predicted;
      {
        ScopedCoutSilence silence;
        predicted = decon_projected_limit( make_backref_input( reference, sample_live ), 0.95,
                                          sample_live/live_over_real, 5.0, true,
                                          DeconLimitType::OneSidedUpperLimit,
                                          ProjectedLimitScoring::JointWithReference,
                                          mc_trials, 4 );
      }

      if( !predicted.valid || !(predicted.median > 0.0) )
        continue;

      ++counted;
      ratios.push_back( predicted.median/actual );
      // `ProjectedLimit` carries only the 16th/84th percentiles, so 68% is the one band there is
      //  anything to check.  A wider band would need more percentiles on the struct.
      if( (actual >= predicted.lower) && (actual <= predicted.upper) ) ++inside68;
    }//for( loop over pairs )

    BOOST_REQUIRE_MESSAGE( counted >= pairs/2,
                          "k=" << k << ": only " << counted << " of " << pairs << " pairs usable" );

    sort( begin(ratios), end(ratios) );
    const double median_ratio = ratios[ratios.size()/2];
    const double c68 = static_cast<double>(inside68)/static_cast<double>(counted);

    cout << "DECON_JOINT_MC_COVERAGE," << k << ',' << counted << ',' << mc_trials << ','
         << median_ratio << ',' << c68 << '\n';

    rows.push_back( { k, median_ratio, c68 } );
  }//for( loop over dwell ratios )

  for( const Row &row : rows )
  {
    BOOST_CHECK_MESSAGE( (row.ratio > 0.75) && (row.ratio < 1.35),
                        "k=" << row.k << ": predicted joint median is " << row.ratio
                        << "x what the joint analysis delivered" );
    BOOST_CHECK_MESSAGE( (row.cover68 > 0.45) && (row.cover68 < 0.90),
                        "k=" << row.k << ": 68% band covered " << row.cover68 );
  }
}// BOOST_AUTO_TEST_CASE( JointProjectedLimitCoverage )


/** Does a two-stage Monte Carlo predict a future measurement's limit, and does its band cover?

 Projecting a dwell today scales the counts by `k` and treats the result as Poisson with variance
 equal to the scaled count, which asserts the loaded spectrum's rate is known exactly.  D-21 measures
 what that costs: the median is right but the spread is too narrow by ~sqrt(1+k).  The deconvolution
 route instead invents the future measurement as the *fitted continuum*, which is smooth by
 construction and so cannot show the region's real shape at all (D-23).

 The alternative tested here draws the future measurement properly, per channel:

     mu ~ Gamma(n + 1, 1)        // the rate, given n counts observed - flat prior posterior
     s  ~ Poisson(k * mu)        // what a dwell of k times the reference would then record

 which gives `E[s|n] = k(n+1)` and `Var[s|n] = k(n+1)(1+k)` - the predictive variance, rather than
 the `k*n` a plain scaling asserts.  Two properties matter and both are measured below:

   - drawing the *rate* rather than a count means a channel that happened to record zero is not
     locked at zero for every trial, which matters because that is most channels in the low-count
     regime a detection limit lives in;
   - each trial is built from the observed counts channel by channel, so curvature and interfering
     lines survive into the trials instead of being replaced by a straight line.

 A median over trials is then a real median prediction, and the percentiles are a real band.  The
 test is whether the band covers: over many independent (reference, future) pairs, the actual future
 limit should fall inside the predicted 68% band about 68% of the time, and inside the 90% band about
 90% of the time.  A median that tracks but a band that does not cover would mean the spread is still
 wrong, which is the whole point of the exercise.
 */
BOOST_AUTO_TEST_CASE( ProjectedMeasurementMonteCarlo )
{
  using namespace DetectionLimitCalc;

  set_data_dir();

  // ~25 minutes of CPU in Release: 400 pairs x 200 trials for Currie and 171 x 198 for the
  //  deconvolution, over three prior strengths and four dwell ratios.  This is the evidence behind
  //  the Jeffreys choice in `projected_limit_prior_strength()`, kept reproducible and opt-in.
  //  Most of that is the deconvolution arm - see `decon_scale` below for why it is not cheaper.
  //  The measured result, at alpha = 0.5:
  //
  //    method | k    | predicted centre / actual | 68% band covers | 90% band covers
  //    Currie | 0.25 | 1.00                      | 0.74            | 0.93
  //    Currie | 1    | 1.04                      | 0.66            | 0.89
  //    Currie | 4    | 1.03                      | 0.66            | 0.88
  //    Currie | 16   | 1.03                      | 0.66            | 0.89
  //    decon  | 0.25 | 1.01                      | 0.67            | 0.89
  //    decon  | 1    | 1.06                      | 0.65            | 0.94
  //    decon  | 4    | 1.13                      | 0.70            | 0.91
  //    decon  | 16   | 1.00                      | 0.70            | 0.89
  //
  //  Currie at 400 pairs x 200 trials, deconvolution at 171 x 198.  The Currie arm is scored on
  //  `detection_limit` and the deconvolution arm on the profile upper limit, because those are the
  //  quantities the two tools display the band beside - see `currie_projected_limit`.
  const char * const enabled = getenv( "INTERSPEC_RUN_DECON_MC_PROJECTION" );
  if( !enabled || ( string( enabled ) != "1" ) )
  {
    BOOST_TEST_MESSAGE( "Set INTERSPEC_RUN_DECON_MC_PROJECTION=1 to run the projection Monte Carlo"
                       " validation." );
    return;
  }

  const string path = SpecUtils::append_path( g_test_file_dir,
                                             "FitPeaksForSource/trinitite_sample_b_background.n42" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( path ), "Background spectrum not at '" << path << "'" );

  SpecMeas file;
  BOOST_REQUIRE_MESSAGE( file.load_N42_file( path ), "Failed loading '" << path << "'" );

  shared_ptr<const SpecUtils::Measurement> truth_rate;
  for( const shared_ptr<const SpecUtils::Measurement> &candidate : file.measurements() )
  {
    if( candidate && (candidate->num_gamma_channels() > 1024) && (candidate->live_time() > 0.0f) )
    {
      truth_rate = candidate;
      break;
    }
  }
  BOOST_REQUIRE( truth_rate );

  const double source_live = truth_rate->live_time();
  const double live_over_real = truth_rate->live_time() / truth_rate->real_time();

  auto drf = make_shared<DetectorPeakResponse>();
  drf->setIntrinsicEfficiencyFormula( "1.0", 2.54f*PhysicalUnits::cm, PhysicalUnits::keV,
                                     0.0f, 0.0f,
                                     DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct );
  drf->setFwhmCoefficients( vector<float>{ 4.0f, 0.02f },
                           DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial );

  const float peak_energy = 1000.0f;
  const float fwhm = static_cast<float>( drf->peakResolutionFWHM( peak_energy ) );
  BOOST_REQUIRE( fwhm > 0.0f );

  // A real, independent measurement of `live_seconds` drawn from the source spectrum's rate.
  const auto draw_truth = [&]( const double live_seconds, mt19937 &generator )
        -> shared_ptr<const SpecUtils::Measurement> {
    const shared_ptr<const vector<float>> rate = truth_rate->gamma_counts();
    auto counts = make_shared<vector<float>>( rate->size(), 0.0f );
    const double scale = live_seconds / source_live;
    for( size_t i = 0; i < rate->size(); ++i )
      counts->at(i) = static_cast<float>( study_poisson( (std::max)(0.0, scale*rate->at(i)),
                                                        generator ) );

    auto m = make_shared<SpecUtils::Measurement>();
    m->set_gamma_counts( counts, static_cast<float>(live_seconds),
                        static_cast<float>(live_seconds/live_over_real) );
    m->set_energy_calibration( truth_rate->energy_calibration() );
    return m;
  };

  // The proposal under test: one Monte Carlo realisation of what a dwell of `k` times the
  //  reference's would record, given only the reference.
  //  \p alpha is the prior strength added to each channel's observed count before the rate is
  //  drawn.  It is the whole design question: alpha=0 reproduces the plain two-stage
  //  Poisson-of-Poisson, which has exactly the right mean but locks a channel that recorded zero at
  //  zero forever; alpha=1 is the flat-prior posterior, which never locks but adds one count to
  //  *every* channel - and a region holding a couple of counts per channel is then inflated by
  //  tens of percent, which the limit inherits as its square root.  alpha=0.5 is Jeffreys.
  //  Only the channels the limit actually reads are drawn - the region, its side channels and a
  //  generous margin.  Everywhere else the expectation is written straight in, since no arm of this
  //  study looks there.  Drawing all 16k channels made the Monte Carlo cost 200x what the limit
  //  did, which capped the pair count low enough that the medians were pure noise.
  const shared_ptr<const SpecUtils::EnergyCalibration> cal = truth_rate->energy_calibration();
  const size_t nchannel = cal->num_channels();
  const size_t window_lower = static_cast<size_t>( (std::max)( 0.0,
                                 cal->channel_for_energy( peak_energy - 5.0*fwhm ) ) );
  const size_t window_upper = (std::min)( nchannel - 1, static_cast<size_t>(
                                 (std::max)( 0.0, cal->channel_for_energy( peak_energy + 5.0*fwhm ) ) ) );

  const auto draw_projection = [&]( const shared_ptr<const SpecUtils::Measurement> &reference,
                                    const double k, const double alpha, mt19937 &generator )
        -> shared_ptr<const SpecUtils::Measurement> {
    const shared_ptr<const vector<float>> observed = reference->gamma_counts();
    auto counts = make_shared<vector<float>>( observed->size(), 0.0f );
    for( size_t i = 0; i < observed->size(); ++i )
      counts->at(i) = static_cast<float>( k*(std::max)( 0.0, static_cast<double>(observed->at(i)) ) );

    for( size_t i = window_lower; i <= window_upper; ++i )
    {
      const double n = (std::max)( 0.0, static_cast<double>( observed->at(i) ) );
      const double rate = study_gamma( n + alpha, generator );
      counts->at(i) = static_cast<float>( study_poisson( k*rate, generator ) );
    }

    auto m = make_shared<SpecUtils::Measurement>();
    m->set_gamma_counts( counts, static_cast<float>(k*reference->live_time()),
                        static_cast<float>(k*reference->real_time()) );
    m->set_energy_calibration( reference->energy_calibration() );
    return m;
  };

  const auto make_input = [&]( const shared_ptr<const SpecUtils::Measurement> &m )
        -> shared_ptr<DeconComputeInput> {
    DeconRoiInfo roi;
    roi.roi_start = peak_energy - 1.25f*fwhm;
    roi.roi_end = peak_energy + 1.25f*fwhm;
    roi.continuum_type = PeakContinuum::OffsetType::Linear;
    roi.cont_norm_method = DeconContinuumNorm::Floating;
    roi.num_lower_side_channels = 4;
    roi.num_upper_side_channels = 4;

    DeconRoiInfo::PeakInfo peak;
    peak.energy = peak_energy;
    peak.fwhm = fwhm;
    peak.counts_per_bq_into_4pi = m->live_time();
    roi.peak_infos.push_back( peak );

    auto input = make_shared<DeconComputeInput>();
    input->activity = 0.0;
    input->distance = 0.0;
    input->include_air_attenuation = false;
    input->shielding_thickness = 0.0;
    input->drf = drf;
    input->measurement = m;
    input->roi_info.push_back( roi );
    input->measurement_model = DeconMeasurementModel::CurrentSpectrum;
    return input;
  };

  // Both limits as a RATE, so measurements of different lengths compare directly.
  //  `detection_limit` for the Currie arm, because that is what `currie_projected_limit` scores and
  //  what the tools quote - a coverage measurement of `upper_limit` would validate a band nothing
  //  ships.  The deconvolution arm keeps the profile upper limit, which is what its band does score.
  const auto currie_rate = [&]( const shared_ptr<const SpecUtils::Measurement> &m ) -> double {
    try
    {
      ScopedCoutSilence silence;
      const CurrieMdaResult r = currie_for_decon_input( make_input( m ) );
      return (m->live_time() > 0.0f) ? (r.detection_limit / m->live_time())
                                     : numeric_limits<double>::quiet_NaN();
    }catch( std::exception & ){ return numeric_limits<double>::quiet_NaN(); }
  };

  const auto decon_rate = [&]( const shared_ptr<const SpecUtils::Measurement> &m ) -> double {
    try
    {
      ScopedCoutSilence silence;
      const DeconActivityOrDistanceLimitResult r =
          get_activity_or_distance_limits( 0.95f, make_input( m ), false, 0.0, 5.0, true );
      return r.foundUpperCl ? r.upperLimit : numeric_limits<double>::quiet_NaN();
    }catch( std::exception & ){ return numeric_limits<double>::quiet_NaN(); }
  };

  const auto quantile_of = []( vector<double> v, const double q ) -> double {
    if( v.empty() )
      return numeric_limits<double>::quiet_NaN();
    sort( begin(v), end(v) );
    const size_t index = (std::min)( v.size() - 1,
                          static_cast<size_t>( q*static_cast<double>(v.size()) ) );
    return v[index];
  };

  cout << "DECON_MC_PROJECTION,method,alpha,k,pairs,mc_trials,median_of_ratios,"
          "ratio_of_centres,cover68,cover90,median_band_width_over_median\n";

  const double reference_live = 600.0;

  const size_t pairs = 400;      // independent (reference, future) draws per cell
  const size_t mc_trials = 200;  // Monte Carlo realisations inside each prediction

  // The deconvolution arm runs a fraction of the pairs and trials, because each of its limits is a
  //  full profile scan.  The base fraction - a seventh of the pairs, a third of the trials - turned
  //  out to be too coarse to gate: at 57 pairs the k = 16 cell read 1.50x, and at 171 pairs the same
  //  cell reads 1.00x, so the deviation was sampling noise in the summary and not a bias in the
  //  prediction.  Hence a default of 3, which is what the recorded table above was measured at, and
  //  INTERSPEC_DECON_MC_SCALE to go finer still.  It scales only this arm: the Currie arm is
  //  already at 400 pairs, and scaling both would cost the square of the factor for nothing.
  double decon_scale = 3.0;
  if( const char * const scale_env = getenv( "INTERSPEC_DECON_MC_SCALE" ) )
    decon_scale = (std::max)( 1.0, atof( scale_env ) );

  struct McRow { string method; double alpha, k, ratio, cover68, cover90; };
  vector<McRow> rows;

  for( const bool currie : { true, false } )
  {
   for( const double alpha : { 0.0, 0.5, 1.0 } )
   {
    for( const double k : { 0.25, 1.0, 4.0, 16.0 } )
    {
      // The deconvolution arm costs a full profile scan per Monte Carlo trial, so it runs a
      //  coarser grid than the Currie arm, which is closed-form arithmetic.
      // A deconvolution trial is a full profile scan where a Currie trial is arithmetic, so it
      //  runs ~7x fewer pairs and ~3x fewer trials for a comparable wall clock.
      const size_t trials = currie ? mc_trials
                                   : static_cast<size_t>( (mc_trials/3)*decon_scale );
      const size_t num_pairs = currie ? pairs
                                      : static_cast<size_t>( (pairs/7)*decon_scale );

      mt19937 generator( study_cell_seed( string("mcproj|") + (currie ? "currie" : "decon")
                                         + "|a=" + study_number_token( alpha )
                                         + "|k=" + study_number_token( k ) ) );

      vector<double> ratio_of_medians;
      size_t inside68 = 0, inside90 = 0, counted = 0;
      vector<double> relative_band;

      // Two summaries of the same question, kept side by side as a check on each other.  The
      //  per-pair ratio p50/actual puts a broad random variable in a denominator, so it could in
      //  principle drift from the ratio of the two distributions' centres; measured, it does not -
      //  the two agree to a few percent in every cell, including the wide ones.  That agreement is
      //  worth keeping visible, because it rules out the summary statistic as an explanation when
      //  a cell comes out far from one.
      vector<double> predicted_centres, actual_centres;

      for( size_t pair = 0; pair < num_pairs; ++pair )
      {
        const shared_ptr<const SpecUtils::Measurement> reference
              = draw_truth( reference_live, generator );
        const shared_ptr<const SpecUtils::Measurement> future
              = draw_truth( k*reference_live, generator );

        const double actual = currie ? currie_rate( future ) : decon_rate( future );
        if( IsNan(actual) )
          continue;

        vector<double> predicted;
        for( size_t trial = 0; trial < trials; ++trial )
        {
          const shared_ptr<const SpecUtils::Measurement> projected
                = draw_projection( reference, k, alpha, generator );
          const double value = currie ? currie_rate( projected ) : decon_rate( projected );
          if( !IsNan(value) )
            predicted.push_back( value );
        }

        if( predicted.size() < (trials/2) )
          continue;

        const double p16 = quantile_of( predicted, 0.16 );
        const double p50 = quantile_of( predicted, 0.50 );
        const double p84 = quantile_of( predicted, 0.84 );
        const double p05 = quantile_of( predicted, 0.05 );
        const double p95 = quantile_of( predicted, 0.95 );

        ++counted;
        if( (actual >= p16) && (actual <= p84) ) ++inside68;
        if( (actual >= p05) && (actual <= p95) ) ++inside90;
        if( p50 > 0.0 )
        {
          ratio_of_medians.push_back( p50/actual );
          predicted_centres.push_back( p50 );
          actual_centres.push_back( actual );
          relative_band.push_back( (p84 - p16)/p50 );
        }
      }//for( loop over independent pairs )

      BOOST_REQUIRE_MESSAGE( counted >= (num_pairs/2),
                            (currie ? "currie" : "decon") << " k=" << k << ": only " << counted
                            << " of " << num_pairs << " pairs usable" );

      const double median_ratio = quantile_of( ratio_of_medians, 0.50 );
      const double centre_ratio = quantile_of( predicted_centres, 0.50 )
                                    / quantile_of( actual_centres, 0.50 );
      const double c68 = static_cast<double>(inside68)/static_cast<double>(counted);
      const double c90 = static_cast<double>(inside90)/static_cast<double>(counted);

      cout << "DECON_MC_PROJECTION," << (currie ? "currie" : "decon") << ',' << alpha << ','
           << k << ',' << counted << ',' << trials << ',' << median_ratio << ',' << centre_ratio
           << ',' << c68 << ',' << c90 << ',' << quantile_of( relative_band, 0.50 ) << '\n';

      rows.push_back( { currie ? "currie" : "decon", alpha, k, centre_ratio, c68, c90 } );
    }//for( loop over dwell ratios )
   }//for( loop over prior strength )
  }//for( Currie / deconvolution )

  // What makes this usable at all: the prediction has to sit on the answer, and its band has to
  //  mean what it says.  Deliberately loose - these are 30-60 pairs, so the coverage estimates
  //  carry several percent of noise; the point is that they are near nominal rather than far off.
  for( const McRow &row : rows )
  {
    // Gate the configuration that actually ships: `projected_limit_prior_strength()` is Jeffreys.
    //  This used to gate alpha = 0, which is the variant the sweep exists to rule out, so the
    //  shipped setting was measured and printed but never asserted.
    if( row.alpha != DetectionLimitCalc::projected_limit_prior_strength() )
      continue;   // the other prior strengths are measured and printed, not gated

    // On the ratio of centres, not the median of per-pair ratios - see where they are computed.
    BOOST_CHECK_MESSAGE( (row.ratio > 0.75) && (row.ratio < 1.35),
                        row.method << " k=" << row.k << ": predicted centre is " << row.ratio
                        << "x where the dwell's outcomes centre" );
    BOOST_CHECK_MESSAGE( (row.cover68 > 0.45) && (row.cover68 < 0.90),
                        row.method << " k=" << row.k << ": 68% band covered " << row.cover68 );
    BOOST_CHECK_MESSAGE( (row.cover90 > 0.70) && (row.cover90 < 1.0001),
                        row.method << " k=" << row.k << ": 90% band covered " << row.cover90 );
  }
}// BOOST_AUTO_TEST_CASE( ProjectedMeasurementMonteCarlo )


/** Is the `BackgroundReference` prediction bias a function of the dwell ratio alone?

 `RealSpectrumDwellProjectionAccuracy` measures it at one continuum level, one region and one
 reference length, and finds it unbiased to ~4% for `k = T_sample/T_reference` up to 1 and about a
 third optimistic at 4 and 16.  If that curve depends only on `k`, it can be calibrated once and
 applied for free; if it also moves with how many counts are under the peak, a correction would need
 a Monte Carlo per limit, which costs ~0.5 s of CPU per region and is not viable interactively.

 So this sweeps `k` against the two things that set the counts under the peak - the continuum level
 and the peak width - independently, which lets the two hypotheses be told apart:

   - bias is a function of `k` alone: every cell at a given `k` agrees within noise;
   - bias is a function of `(k, B)`: cells at a given `k` spread systematically with `B`, the
     region's background counts, and continuum and FWHM enter only through their product.

 `B` is varied over four decades here (0.5 to 4000 counts) by moving continuum and FWHM in ways that
 sometimes agree and sometimes cancel, so the two are separable rather than confounded.
 */
BOOST_AUTO_TEST_CASE( PredictedSampleBiasParameterSweep )
{
  using namespace DetectionLimitCalc;

  set_data_dir();

  // 64 cells of 150 realisations, two profile scans each: ~5.5 minutes of CPU in Release and far
  //  more in Debug, so it is opt-in like the coverage study rather than part of every run.  Its
  //  conclusion is recorded in the register (D-17) and the always-on
  //  `ContinuumMisspecificationCollapsesLimit` covers the finding that came out of it.
  const char * const enabled = getenv( "INTERSPEC_RUN_DECON_PREDICTION_SWEEP" );
  if( !enabled || ( string( enabled ) != "1" ) )
  {
    BOOST_TEST_MESSAGE( "Set INTERSPEC_RUN_DECON_PREDICTION_SWEEP=1 to run the expected-counts-bias sweep." );
    return;
  }

  const size_t num_channels = 1024;
  const float channel_width = 2.0f;
  const float peak_energy = 1173.0f;

  auto calibration = make_shared<SpecUtils::EnergyCalibration>();
  calibration->set_polynomial( num_channels, { 0.0f, channel_width }, {} );

  auto drf = make_shared<DetectorPeakResponse>();
  drf->setIntrinsicEfficiencyFormula( "1.0", 2.54f*PhysicalUnits::cm, PhysicalUnits::keV,
                                     0.0f, 0.0f,
                                     DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct );
  drf->setFwhmCoefficients( vector<float>{ 64.0f, 0.0f },
                           DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial );

  // Draws a flat-continuum spectrum of mean `counts_per_channel`, at exposure `live`.
  const auto draw = [&]( const double counts_per_channel, const double live, mt19937 &generator )
        -> shared_ptr<const SpecUtils::Measurement> {
    auto counts = make_shared<vector<float>>( num_channels, 0.0f );
    for( size_t c = 0; c < num_channels; ++c )
      counts->at(c) = static_cast<float>( study_poisson( counts_per_channel, generator ) );

    auto m = make_shared<SpecUtils::Measurement>();
    m->set_gamma_counts( counts, static_cast<float>(live), static_cast<float>(live) );
    m->set_energy_calibration( calibration );
    return m;
  };

  const auto make_input = [&]( const shared_ptr<const SpecUtils::Measurement> &reference,
                               const double fwhm ) -> shared_ptr<DeconComputeInput> {
    DeconRoiInfo roi;
    roi.roi_start = static_cast<float>( peak_energy - 1.25*fwhm );
    roi.roi_end = static_cast<float>( peak_energy + 1.25*fwhm );
    roi.continuum_type = PeakContinuum::OffsetType::Linear;
    roi.cont_norm_method = DeconContinuumNorm::Floating;
    roi.num_lower_side_channels = 4;
    roi.num_upper_side_channels = 4;

    DeconRoiInfo::PeakInfo peak;
    peak.energy = peak_energy;
    peak.fwhm = static_cast<float>( fwhm );
    // Reference live time, so "activity" reads as counts at the reference exposure and the ratio
    //  this case reports is independent of the normalisation.
    peak.counts_per_bq_into_4pi = reference->live_time();
    roi.peak_infos.push_back( peak );

    auto input = make_shared<DeconComputeInput>();
    input->activity = 0.0;
    input->distance = 0.0;
    input->include_air_attenuation = false;
    input->shielding_thickness = 0.0;
    input->drf = drf;
    input->measurement = reference;
    input->roi_info.push_back( roi );
    input->measurement_model = DeconMeasurementModel::BackgroundReference;
    return input;
  };

  const auto limit_of = [&]( const shared_ptr<const DeconComputeInput> &input,
                             const double max_search ) -> double {
    try
    {
      ScopedCoutSilence silence;
      const DeconActivityOrDistanceLimitResult r =
          get_activity_or_distance_limits( 0.95f, input, false, 0.0, max_search, true );
      return r.foundUpperCl ? r.upperLimit : numeric_limits<double>::quiet_NaN();
    }catch( std::exception & )
    {
      return numeric_limits<double>::quiet_NaN();
    }
  };

  const auto median_of = []( vector<double> v ) -> double {
    if( v.empty() )
      return numeric_limits<double>::quiet_NaN();
    const size_t mid = v.size()/2;
    nth_element( v.begin(), v.begin() + mid, v.end() );
    return v[mid];
  };

  cout << "DECON_PREDICTION_SWEEP,k,continuum_per_channel,fwhm_keV,roi_channels,"
          "roi_background_counts,reps,attempted,predicted,joint_actual,predicted_over_joint\n";

  const double reference_live = 1.0;
  const size_t reps = 150;
  // Wall-clock guard: the low-continuum cells are the slow ones (D-18), and a cell truncated by its
  //  slice must report how many reps it actually ran rather than assert on noise.
  const double seconds_per_cell = 12.0;

  struct SweepRow
  {
    double k, continuum, fwhm, background, ratio;
    size_t reps;
  };
  vector<SweepRow> rows;

  for( const double k : { 0.25, 1.0, 4.0, 16.0 } )
  {
    for( const double continuum : { 0.1, 1.0, 10.0, 100.0 } )
    {
      for( const double fwhm : { 4.0, 8.0, 16.0, 32.0 } )
      {
        // The "asimov" token is the old name for the expected-counts step, kept because this string
        //  *is* the seed: changing a character changes every draw, and so every number this study
        //  published.  It is an opaque token, not a description.
        const string cell = "asimov|k=" + study_number_token( k )
                            + "|c=" + study_number_token( continuum )
                            + "|f=" + study_number_token( fwhm );
        mt19937 generator( study_cell_seed( cell ) );

        const size_t roi_channels
              = static_cast<size_t>( std::round( 2.5*fwhm / channel_width ) );
        const double background = continuum * static_cast<double>( roi_channels );

        // Comfortably above the limit's own scale, so the scan brackets it without wasting grid
        //  points: a 95% bound sits near 3.3*sqrt(B_sample) counts, which is
        //  3.3*sqrt(B/k) once expressed at the reference exposure.
        const double max_search = 50.0 * ( 3.3*sqrt( (std::max)(background,1.0)/k ) + 2.0 );

        vector<double> prediction, joint;
        size_t attempted = 0;
        const std::chrono::steady_clock::time_point started = std::chrono::steady_clock::now();

        for( size_t rep = 0; rep < reps; ++rep )
        {
          const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - started;
          if( (rep > 8) && (elapsed.count() > seconds_per_cell) )
            break;
          ++attempted;

          const shared_ptr<const SpecUtils::Measurement> reference
                = draw( continuum, reference_live, generator );
          const shared_ptr<const SpecUtils::Measurement> sample
                = draw( continuum*k, reference_live*k, generator );

          const shared_ptr<DeconComputeInput> predictive = make_input( reference, fwhm );
          predictive->sample_exposure = reference_live*k;
          const double p = limit_of( predictive, max_search );

          const shared_ptr<DeconComputeInput> observed = make_input( reference, fwhm );
          observed->observed_sample = sample;
          const double j = limit_of( observed, max_search );

          if( !IsNan(p) && !IsNan(j) )
          {
            prediction.push_back( p );
            joint.push_back( j );
          }
        }//for( loop over repetitions )

        const double mp = median_of( prediction );
        const double mj = median_of( joint );
        const double ratio = (mj > 0.0) ? (mp/mj) : numeric_limits<double>::quiet_NaN();

        cout << "DECON_PREDICTION_SWEEP," << k << ',' << continuum << ',' << fwhm << ','
             << roi_channels << ',' << background << ',' << prediction.size() << ','
             << attempted << ',' << mp << ',' << mj << ',' << ratio << '\n';

        // Only cells with enough successful pairs to have a meaningful median are carried into the
        //  hypothesis test below; the row is still printed either way.
        if( (prediction.size() >= 40) && !IsNan(ratio) )
          rows.push_back( { k, continuum, fwhm, background, ratio, prediction.size() } );
      }//for( loop over FWHM )
    }//for( loop over continuum )
  }//for( loop over dwell ratio )

  BOOST_REQUIRE_MESSAGE( rows.size() >= 40,
                        "only " << rows.size() << " of 64 cells produced a usable median" );

  // The question this case exists to answer: at fixed k, does the bias hold still as the counts
  //  under the peak move over four decades?  Reported as the spread of the ratio within each k.
  cout << "DECON_PREDICTION_SWEEP_BY_K,k,cells,min_ratio,median_ratio,max_ratio,spread,"
          "min_background,max_background\n";

  for( const double k : { 0.25, 1.0, 4.0, 16.0 } )
  {
    vector<double> ratios;
    double min_b = numeric_limits<double>::max(), max_b = 0.0;
    for( const SweepRow &row : rows )
    {
      if( row.k != k )
        continue;
      ratios.push_back( row.ratio );
      min_b = (std::min)( min_b, row.background );
      max_b = (std::max)( max_b, row.background );
    }

    if( ratios.size() < 4 )
      continue;

    sort( begin(ratios), end(ratios) );
    const double lo = ratios.front(), hi = ratios.back();
    const double mid = median_of( ratios );

    cout << "DECON_PREDICTION_SWEEP_BY_K," << k << ',' << ratios.size() << ',' << lo << ',' << mid
         << ',' << hi << ',' << (hi - lo) << ',' << min_b << ',' << max_b << '\n';
  }//for( loop over dwell ratio )

  // The measured answer, pinned so it cannot quietly drift: on a continuum the model can represent,
  //  the expected-counts step is close to unbiased at every dwell ratio and over four decades of counts.
  //  What bias there is on real spectra is the continuum *model*, not the noise-free idealisation - see
  //  `ContinuumMisspecificationCollapsesLimit`.
  //
  //  Cells whose reference region holds fewer than ~5 counts are excluded rather than gated: below
  //  that the median of a limit is not a stable quantity, and the ratio there ranges 0.21 to 1.51
  //  in this very sweep.  They are still printed, and the register records them.
  size_t gated = 0;
  for( const SweepRow &row : rows )
  {
    if( row.background < 5.0 )
      continue;

    ++gated;
    BOOST_CHECK_MESSAGE( (row.ratio > 0.70) && (row.ratio < 1.35),
                        "k=" << row.k << " continuum=" << row.continuum << " fwhm=" << row.fwhm
                        << " (B=" << row.background << "): expected-counts/joint is " << row.ratio
                        << ", outside the band a correctly specified continuum gives" );
  }

  BOOST_CHECK_MESSAGE( gated >= 40, "only " << gated << " cells had enough counts to gate" );
}// BOOST_AUTO_TEST_CASE( PredictedSampleBiasParameterSweep )


/** Is the `BackgroundReference` bias the expected-counts step, or the continuum model being wrong?

 `PredictedSampleBiasParameterSweep` finds essentially no bias on a *flat* continuum at any dwell ratio, once
 the region holds more than a few counts — yet `RealSpectrumDwellProjectionAccuracy` finds 0.70 and
 0.63 at `k` = 4 and 16 on a real background whose region holds ~140 counts. Something other than
 the noise-free idealisation is doing the work.

 The candidate is misspecification. The expected-counts sample block is built from the continuum *fitted to
 the reference*, so it matches the fitted model exactly, by construction — it cannot show curvature
 or structure the model lacks. A real sample can, and does. That mismatch contributes counts the fit
 must absorb or attribute to signal, and it grows **linearly** with exposure while Poisson noise
 grows as its square root, so it is invisible at short dwells and dominates at long ones. Which is
 the shape the real spectrum shows.

 Four truths, all analysed with the same `Linear` continuum:
   - flat, and sloped: both exactly representable, so the model is *correct* — the control;
   - curved: a quadratic the line cannot follow;
   - bump: a flat continuum with a small unmodelled Gaussian off to one side of the region, which is
     the interference case a real background presents.

 If misspecification is the mechanism, the two correct-model truths stay near 1 at every `k` and the
 two wrong-model truths fall away above `k = 1`.
 */
BOOST_AUTO_TEST_CASE( ContinuumMisspecificationCollapsesLimit )
{
  using namespace DetectionLimitCalc;

  set_data_dir();

  const size_t num_channels = 1024;
  const float channel_width = 2.0f;
  const float peak_energy = 1173.0f;
  const double fwhm = 32.0;                 // 40 channels across the region
  const double continuum_level = 10.0;      // counts/channel at the reference exposure
  const double roi_half_width = 1.25*fwhm;

  auto calibration = make_shared<SpecUtils::EnergyCalibration>();
  calibration->set_polynomial( num_channels, { 0.0f, channel_width }, {} );

  auto drf = make_shared<DetectorPeakResponse>();
  drf->setIntrinsicEfficiencyFormula( "1.0", 2.54f*PhysicalUnits::cm, PhysicalUnits::keV,
                                     0.0f, 0.0f,
                                     DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct );
  drf->setFwhmCoefficients( vector<float>{ 64.0f, 0.0f },
                           DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial );

  struct Shape
  {
    const char *name;
    // Mean counts in a channel whose centre sits `dE` from the peak, at unit exposure.
    std::function<double(double)> mean;
    bool model_is_correct;
  };

  const double bump_offset = 0.75*roi_half_width;   // inside the region, off to one side
  const double bump_sigma = fwhm / PhysicalUnits::fwhm_nsigma;

  const vector<Shape> shapes = {
    { "flat",   []( const double ){ return 10.0; }, true },
    { "sloped", [&]( const double dE ){ return continuum_level*(1.0 + 0.30*dE/roi_half_width); }, true },
    // ~5% peak-to-line deviation across the region; a line fitted to a parabola is high at the ends
    //  and low in the middle, and no linear continuum can remove it.
    { "curved", [&]( const double dE ){
        const double u = dE/roi_half_width;
        return continuum_level*(1.0 + 0.15*u*u);
      }, false },
    // A 4% unmodelled line sitting off-centre in the region.
    { "bump",   [&]( const double dE ){
        const double z = (dE - bump_offset)/bump_sigma;
        return continuum_level*(1.0 + 0.60*exp( -0.5*z*z ));
      }, false },
  };

  const auto draw = [&]( const Shape &shape, const double exposure, mt19937 &generator )
        -> shared_ptr<const SpecUtils::Measurement> {
    auto counts = make_shared<vector<float>>( num_channels, 0.0f );
    for( size_t c = 0; c < num_channels; ++c )
    {
      const double centre = 0.5*( calibration->energy_for_channel(c)
                                 + calibration->energy_for_channel(c+1) );
      const double dE = centre - peak_energy;
      // Only shape the region and its surroundings; far away it is flat, which keeps the fixture's
      //  total counts from depending on the shape.
      const double mean = (fabs(dE) <= 2.0*roi_half_width) ? shape.mean( dE ) : continuum_level;
      counts->at(c) = static_cast<float>( study_poisson( exposure*mean, generator ) );
    }

    auto m = make_shared<SpecUtils::Measurement>();
    m->set_gamma_counts( counts, static_cast<float>(exposure), static_cast<float>(exposure) );
    m->set_energy_calibration( calibration );
    return m;
  };

  const auto make_input = [&]( const shared_ptr<const SpecUtils::Measurement> &reference )
        -> shared_ptr<DeconComputeInput> {
    DeconRoiInfo roi;
    roi.roi_start = static_cast<float>( peak_energy - roi_half_width );
    roi.roi_end = static_cast<float>( peak_energy + roi_half_width );
    roi.continuum_type = PeakContinuum::OffsetType::Linear;
    roi.cont_norm_method = DeconContinuumNorm::Floating;
    roi.num_lower_side_channels = 4;
    roi.num_upper_side_channels = 4;

    DeconRoiInfo::PeakInfo peak;
    peak.energy = peak_energy;
    peak.fwhm = static_cast<float>( fwhm );
    peak.counts_per_bq_into_4pi = reference->live_time();
    roi.peak_infos.push_back( peak );

    auto input = make_shared<DeconComputeInput>();
    input->activity = 0.0;
    input->distance = 0.0;
    input->include_air_attenuation = false;
    input->shielding_thickness = 0.0;
    input->drf = drf;
    input->measurement = reference;
    input->roi_info.push_back( roi );
    input->measurement_model = DeconMeasurementModel::BackgroundReference;
    return input;
  };

  const auto limit_of = [&]( const shared_ptr<const DeconComputeInput> &input ) -> double {
    try
    {
      ScopedCoutSilence silence;
      const DeconActivityOrDistanceLimitResult r =
          get_activity_or_distance_limits( 0.95f, input, false, 0.0, 20000.0, true );
      return r.foundUpperCl ? r.upperLimit : numeric_limits<double>::quiet_NaN();
    }catch( std::exception & )
    {
      return numeric_limits<double>::quiet_NaN();
    }
  };

  const auto median_of = []( vector<double> v ) -> double {
    if( v.empty() )
      return numeric_limits<double>::quiet_NaN();
    const size_t mid = v.size()/2;
    nth_element( v.begin(), v.begin() + mid, v.end() );
    return v[mid];
  };

  cout << "DECON_SHAPE_SENSITIVITY,shape,model_is_correct,k,reps,rms_misspecification_frac,"
          "predicted,joint_actual,predicted_over_joint,sample_only,sample_only_over_flat_expectation\n";

  const double reference_live = 1.0;
  // Enough to resolve the effect - which is a factor of several, not a few percent - while keeping
  //  this always-on rather than opt-in.  The ratio of two medians carries ~5% noise at this count,
  //  which is why the correct-model gate below is a 20% band.
  const size_t reps = 60;

  struct ShapeRow { string shape; bool correct; double k, ratio; };
  vector<ShapeRow> rows;

  for( const Shape &shape : shapes )
  {
    // How far the truth sits from the best straight line through it, over the region, as a fraction
    //  of the continuum level.  This is the quantity the mechanism says should predict the bias.
    double sum_w = 0.0, sum_x = 0.0, sum_y = 0.0, sum_xx = 0.0, sum_xy = 0.0;
    vector<pair<double,double>> region;
    for( size_t c = 0; c < num_channels; ++c )
    {
      const double centre = 0.5*( calibration->energy_for_channel(c)
                                 + calibration->energy_for_channel(c+1) );
      const double dE = centre - peak_energy;
      if( fabs(dE) > roi_half_width )
        continue;
      const double y = shape.mean( dE );
      region.emplace_back( dE, y );
      sum_w += 1.0; sum_x += dE; sum_y += y; sum_xx += dE*dE; sum_xy += dE*y;
    }
    const double denominator = sum_w*sum_xx - sum_x*sum_x;
    const double slope = (fabs(denominator) > 0.0) ? ((sum_w*sum_xy - sum_x*sum_y)/denominator) : 0.0;
    const double intercept = (sum_y - slope*sum_x)/sum_w;
    double sum_sq = 0.0;
    for( const pair<double,double> &p : region )
    {
      const double residual = p.second - (intercept + slope*p.first);
      sum_sq += residual*residual;
    }
    const double rms_frac = sqrt( sum_sq/region.size() ) / continuum_level;

    for( const double k : { 0.25, 1.0, 4.0, 16.0 } )
    {
      mt19937 generator( study_cell_seed( string("shape|") + shape.name
                                         + "|k=" + study_number_token( k ) ) );

      vector<double> prediction, joint, sample_only;
      for( size_t rep = 0; rep < reps; ++rep )
      {
        const shared_ptr<const SpecUtils::Measurement> reference
              = draw( shape, reference_live, generator );
        const shared_ptr<const SpecUtils::Measurement> sample
              = draw( shape, reference_live*k, generator );

        const shared_ptr<DeconComputeInput> predictive = make_input( reference );
        predictive->sample_exposure = reference_live*k;
        const double p = limit_of( predictive );

        const shared_ptr<DeconComputeInput> observed = make_input( reference );
        observed->observed_sample = sample;
        const double j = limit_of( observed );

        // The control that decides whether any collapse below belongs to the joint path this
        //  increment added, or to the ordinary single-spectrum path the tools have always shipped.
        //  Same misspecified truth, no reference block, no expected-counts step - just `CurrentSpectrum`.
        const shared_ptr<DeconComputeInput> alone = make_input( sample );
        alone->measurement_model = DeconMeasurementModel::CurrentSpectrum;
        alone->sample_exposure = 0.0;
        const double s = limit_of( alone );

        if( !IsNan(p) && !IsNan(j) )
        {
          prediction.push_back( p );
          joint.push_back( j );
        }
        if( !IsNan(s) )
          sample_only.push_back( s );
      }//for( loop over repetitions )

      BOOST_REQUIRE_MESSAGE( prediction.size() > reps/2,
                            shape.name << " k=" << k << ": only " << prediction.size()
                            << " of " << reps << " pairs gave both limits" );

      const double mp = median_of( prediction );
      const double mj = median_of( joint );
      const double ms = median_of( sample_only );
      const double ratio = mp/mj;

      cout << "DECON_SHAPE_SENSITIVITY," << shape.name << ',' << (shape.model_is_correct ? 1 : 0)
           << ',' << k << ',' << prediction.size() << ',' << rms_frac << ',' << mp << ',' << mj
           << ',' << ratio << ',' << ms << ',' << (mp > 0.0 ? ms/mp : 0.0) << '\n';

      rows.push_back( { shape.name, shape.model_is_correct, k, ratio } );
    }//for( loop over dwell ratio )
  }//for( loop over continuum shapes )

  // The control has to hold: where the continuum model can represent the truth, nothing here moves
  //  by a factor, at any dwell ratio.  The band is the one `PredictedSampleBiasParameterSweep` measured over
  //  64 flat-continuum cells rather than a tighter one picked here, so the two cases agree on what
  //  "no effect" means.
  //
  //  Note `sloped` at k=16 sits near the top of it - 1.16 at 150 realisations, 1.22 at 60 - where
  //  `flat` at the same k is 1.02.  That is a real residual and it is not explained: a slope is
  //  exactly representable by the Linear continuum, so misspecification cannot be the cause.  It is
  //  an order of magnitude below the effect this case exists to show, and is recorded rather than
  //  chased.
  for( const ShapeRow &row : rows )
  {
    if( row.correct )
      BOOST_CHECK_MESSAGE( (row.ratio > 0.70) && (row.ratio < 1.35),
                          row.shape << " k=" << row.k << ": ratio " << row.ratio
                          << ", outside the band a correctly specified continuum gives" );
  }

  // The finding this case exists to surface, stated rather than asserted - asserting it would turn
  //  a future fix into a test failure.  Where the `Linear` continuum cannot represent the truth,
  //  the limit collapses, and it does so on the *ordinary* single-spectrum path as much as on the
  //  prediction path, so this is not a property of the background-reference model.  Every coverage
  //  cell in this suite runs on a flat continuum and so cannot see it.
  for( const ShapeRow &row : rows )
  {
    if( !row.correct && (row.k >= 4.0) && (row.ratio > 2.0) )
      BOOST_TEST_MESSAGE( "Continuum misspecification (" << row.shape << ", k=" << row.k
                         << "): expected-counts/joint = " << row.ratio
                         << " - driven by the joint and single-spectrum limits collapsing, not by"
                            " the expected-counts prediction, which barely moves.  See DECON_SHAPE_SENSITIVITY." );
  }
}// BOOST_AUTO_TEST_CASE( ContinuumMisspecificationCollapsesLimit )


/** `DeconComputeInput::observed_sample` turns the background-reference prediction into an ordinary
 joint analysis of a sample that has actually been taken.

 Pins the three things that must change with it, and the one that must not:
   - unset, every number is bit-identical to what the prediction always gave;
   - set, the result is *not* a predicted sensitivity and must not be worded as one;
   - set, the sample's own live time supersedes `sample_exposure`;
   - set, the answer bounds the sample rather than predicting - so with a sample that really does
     hold signal, zero can be excluded, which an null expectation can never do.
 */
BOOST_AUTO_TEST_CASE( ObservedSampleMakesBackgroundReferenceJoint )
{
  set_data_dir();

  const double continuum_per_second = 100.0;
  const double reference_seconds = 4.0;
  const double sample_seconds = 4.0;

  const shared_ptr<const DetectionLimitCalc::DeconComputeInput> predictive
      = make_background_reference_input( continuum_per_second, reference_seconds, sample_seconds );

  const auto limits_of = []( const shared_ptr<const DetectionLimitCalc::DeconComputeInput> &in ){
    ScopedCoutSilence silence;
    return DetectionLimitCalc::get_activity_or_distance_limits( 0.95f, in, false, 0.0, 5000.0, true );
  };

  const DetectionLimitCalc::DeconActivityOrDistanceLimitResult predicted = limits_of( predictive );
  BOOST_REQUIRE( predicted.foundUpperCl );
  BOOST_CHECK( predicted.is_predicted_sensitivity );

  // A sample drawn from the same flat continuum, with a deliberate excess in the peak region so
  //  that the joint fit has something to find.
  const size_t nchannel = predictive->measurement->num_gamma_channels();
  auto sample_counts = make_shared<vector<float>>( *predictive->measurement->gamma_counts() );
  BOOST_REQUIRE_EQUAL( sample_counts->size(), nchannel );

  const float peak_energy = predictive->roi_info.front().peak_infos.front().energy;
  const size_t peak_channel = predictive->measurement->find_gamma_channel( peak_energy );
  for( size_t i = (peak_channel > 3 ? peak_channel - 3 : 0);
      (i <= peak_channel + 3) && (i < nchannel); ++i )
  {
    sample_counts->at(i) += 400.0f;
  }

  auto sample = make_shared<SpecUtils::Measurement>();
  sample->set_gamma_counts( sample_counts, static_cast<float>(sample_seconds),
                           static_cast<float>(sample_seconds) );
  sample->set_energy_calibration( predictive->measurement->energy_calibration() );

  auto joint = make_shared<DetectionLimitCalc::DeconComputeInput>( *predictive );
  joint->observed_sample = sample;
  // Deliberately wrong, to prove the sample's own live time is what gets used.
  joint->sample_exposure = 1000.0*sample_seconds;

  const DetectionLimitCalc::DeconActivityOrDistanceLimitResult observed = limits_of( joint );
  BOOST_REQUIRE_MESSAGE( observed.foundUpperCl, "no joint limit: " << observed.errorMessage );

  BOOST_CHECK_MESSAGE( !observed.is_predicted_sensitivity,
                      "a joint analysis of a measured sample is a bound, not a prediction" );
  BOOST_CHECK_CLOSE_FRACTION( observed.sampleExposure, sample_seconds, 1.0E-6 );
  BOOST_CHECK_MESSAGE( observed.foundLowerCl,
                      "a sample holding a real excess must be able to exclude zero, which the"
                      " null expectation never can" );

  // Leaving `observed_sample` null must change nothing at all about the prediction.
  auto unchanged = make_shared<DetectionLimitCalc::DeconComputeInput>( *predictive );
  unchanged->observed_sample = nullptr;
  const DetectionLimitCalc::DeconActivityOrDistanceLimitResult again = limits_of( unchanged );
  BOOST_REQUIRE( again.foundUpperCl );
  BOOST_CHECK_EQUAL( again.upperLimit, predicted.upperLimit );
  BOOST_CHECK_EQUAL( again.is_predicted_sensitivity, predicted.is_predicted_sensitivity );
  BOOST_CHECK_EQUAL( again.sampleExposure, predicted.sampleExposure );

  // A sample that cannot be paired channel-for-channel with the reference is rejected rather than
  //  silently mis-indexed.
  auto mismatched = make_shared<SpecUtils::Measurement>();
  mismatched->set_gamma_counts( make_shared<vector<float>>( nchannel/2, 1.0f ),
                               static_cast<float>(sample_seconds),
                               static_cast<float>(sample_seconds) );
  auto bad = make_shared<DetectionLimitCalc::DeconComputeInput>( *predictive );
  bad->observed_sample = mismatched;
  BOOST_CHECK_THROW( DetectionLimitCalc::decon_compute_peaks( *bad ), std::exception );
}// BOOST_AUTO_TEST_CASE( ObservedSampleMakesBackgroundReferenceJoint )


BOOST_AUTO_TEST_CASE( BackgroundReferenceApproachesIdeal )
{
  set_data_dir();

  // A finite reference background is itself uncertain, and that uncertainty belongs in the
  // predicted sensitivity.  Lengthening the reference while holding the predicted measurement
  // fixed must therefore tighten the answer monotonically, converging on the ideal
  // known-background case - rather than being independent of the reference all along, which is
  // what the old zero-signal-continuum treatment effectively assumed.
  const double continuum_per_second = 100.0;
  const double sample_seconds = 1.0;

  vector<double> limits;
  for( const double reference_seconds : { 0.25, 1.0, 4.0, 20.0, 1000.0 } )
  {
    const shared_ptr<const DetectionLimitCalc::DeconComputeInput> input
        = make_background_reference_input( continuum_per_second, reference_seconds, sample_seconds );

    DetectionLimitCalc::DeconActivityOrDistanceLimitResult result;
    {
      ScopedCoutSilence silence;
      result = DetectionLimitCalc::get_activity_or_distance_limits( 0.95f, input, false,
                                                                   0.0, 5000.0, true );
    }

    BOOST_REQUIRE_MESSAGE( result.foundUpperCl,
                          "no sensitivity at reference " << reference_seconds << " s: "
                          << result.errorMessage );
    BOOST_CHECK_MESSAGE( result.is_predicted_sensitivity,
                        "a background-reference result must be flagged as a prediction" );
    BOOST_CHECK_MESSAGE( !result.foundLowerCl,
                        "the null expectation holds no excess, so zero can never be excluded" );
    limits.push_back( result.upperLimit );
  }

  for( size_t i = 1; i < limits.size(); ++i )
    BOOST_CHECK_MESSAGE( limits[i] <= limits[i-1]*1.01,
                        "limit rose from " << limits[i-1] << " to " << limits[i]
                        << " as the reference got longer" );

  // Monotonicity alone is far too weak: an implementation propagating the reference's statistics
  // at a tenth of the correct strength is still monotone and still converges.  The SIZE of the
  // degradation has to be pinned.
  const double reference_seconds[] = { 0.25, 1.0, 4.0, 20.0, 1000.0 };

  cout << "DECON_BACKREF_LADDER,reference_seconds,sample_seconds,limit,limit_over_converged\n";
  for( size_t i = 0; i < limits.size(); ++i )
    cout << "DECON_BACKREF_LADDER," << reference_seconds[i] << ',' << sample_seconds << ','
         << limits[i] << ',' << (limits[i]/limits[4]) << '\n';

  // Theoretical bracket, valid whatever the region width.  Subtracting a reference of exposure
  // T_ref from a sample of T_s inflates the background variance by (1 + T_s/T_ref), so
  // sqrt(1 + T_s/T_ref) is the degradation when the reference is the *only* thing determining the
  // continuum.  Here it is not: the region's own signal-free channels also constrain it, so the
  // real degradation is strictly smaller - but still strictly greater than 1, which is what
  // "the reference's statistics reach the limit at all" means.
  for( size_t i = 0; i + 1 < limits.size(); ++i )
  {
    const double ceiling = sqrt( 1.0 + sample_seconds/reference_seconds[i] );
    const double observed = limits[i] / limits[4];
    BOOST_CHECK_MESSAGE( (observed > 1.0) && (observed < 1.05*ceiling),
                        "reference " << reference_seconds[i] << " s: limit/converged is "
                        << observed << ", outside the bracket (1, " << ceiling << "]" );
  }

  // ...and pin the two informative cells to what the construction actually produces, tightly
  // enough to reject a materially weakened propagation.  An implementation carrying the reference
  // statistics at a tenth strength would give roughly 1.03 and 1.02 here.
  BOOST_CHECK_CLOSE( limits[0]/limits[4], 1.3218, 8.0 );
  BOOST_CHECK_CLOSE( limits[1]/limits[4], 1.1645, 8.0 );
}// BOOST_AUTO_TEST_CASE( BackgroundReferenceApproachesIdeal )


BOOST_AUTO_TEST_CASE( BackgroundReferenceScalesWithSampleExposure )
{
  set_data_dir();

  // In the background-limited regime the detectable *counts* grow like sqrt(T_s) while the counts
  // per unit activity grow like T_s, so the detectable activity improves like 1/sqrt(T_s).  A long
  // reference keeps its own statistics negligible so the scaling is clean.
  const double continuum_per_second = 100.0;
  const double reference_seconds = 2000.0;

  const auto limit_for = [&]( const double sample_seconds ) -> double {
    const shared_ptr<const DetectionLimitCalc::DeconComputeInput> input
        = make_background_reference_input( continuum_per_second, reference_seconds, sample_seconds );

    DetectionLimitCalc::DeconActivityOrDistanceLimitResult result;
    {
      ScopedCoutSilence silence;
      result = DetectionLimitCalc::get_activity_or_distance_limits( 0.95f, input, false,
                                                                   0.0, 50000.0, true );
    }
    BOOST_REQUIRE_MESSAGE( result.foundUpperCl,
                          "no sensitivity for a " << sample_seconds << " s measurement: "
                          << result.errorMessage );
    BOOST_CHECK_CLOSE( result.sampleExposure, sample_seconds, 1.0E-3 );
    return result.upperLimit;
  };

  const double one = limit_for( 1.0 );
  const double four = limit_for( 4.0 );
  const double sixteen = limit_for( 16.0 );

  BOOST_CHECK_MESSAGE( four < one, "a longer measurement must be more sensitive" );
  BOOST_CHECK_MESSAGE( sixteen < four, "a longer measurement must be more sensitive" );

  // Quadrupling the dwell should halve the detectable activity.  A generous band, since at these
  // count rates the Poisson limit is not exactly its asymptotic square-root form.
  BOOST_CHECK_CLOSE( one/four, 2.0, 20.0 );
  BOOST_CHECK_CLOSE( four/sixteen, 2.0, 20.0 );
}// BOOST_AUTO_TEST_CASE( BackgroundReferenceScalesWithSampleExposure )


BOOST_AUTO_TEST_CASE( BackgroundReferencePeaksMatchTheSpectrumTheyAreDrawnOver )
{
  using namespace DetectionLimitCalc;

  set_data_dir();

  // `fit_peaks` are drawn over `DeconComputeInput::measurement`, which under a background reference
  //  is the *reference* spectrum.  The Gaussian used to be handed back projected to the future
  //  measurement while the continuum on the very same peak stayed at reference exposure, so the
  //  drawn peak stood too tall over both, by exactly the exposure ratio.
  //
  //  Lengthening only the predicted measurement must therefore leave the returned peak alone: the
  //  reference data, and so the continuum solved from it, have not changed.  What changes is the
  //  activity the limit reports, which other cases cover.
  const double continuum_per_second = 10.0;
  const double reference_seconds = 4.0;
  const double activity = 25.0;

  vector<double> amplitudes, continuum_areas;
  for( const double sample_seconds : { 1.0, 4.0, 16.0 } )
  {
    shared_ptr<DeconComputeInput> input = make_shared<DeconComputeInput>(
        *make_background_reference_input( continuum_per_second, reference_seconds,
                                          sample_seconds ) );
    input->activity = activity;

    const DeconComputeResults results = decon_compute_peaks( *input );

    BOOST_REQUIRE_EQUAL( results.fit_peaks.size(), 1 );
    const PeakDef &peak = results.fit_peaks.front();
    BOOST_REQUIRE( peak.continuum() );

    const double lower = peak.continuum()->lowerEnergy();
    const double upper = peak.continuum()->upperEnergy();
    BOOST_REQUIRE( upper > lower );

    amplitudes.push_back( peak.amplitude() );
    // Linear continuum, so the peak list only matters for stepped types.
    const vector<shared_ptr<const PeakDef>> no_roi_peaks;
    continuum_areas.push_back( peak.continuum()->offset_integral( lower, upper, input->measurement,
                                                                  no_roi_peaks ) );

    BOOST_CHECK_CLOSE( results.exposure_ratio, sample_seconds/reference_seconds, 1.0E-6 );
  }//for( loop over predicted measurement lengths )

  // The Gaussian describes the reference spectrum, which did not change, so it must not move at
  //  all.  This is the assertion the defect fails: it used to scale with the predicted length.
  for( size_t i = 1; i < amplitudes.size(); ++i )
  {
    BOOST_CHECK_MESSAGE( fabs( amplitudes[i] - amplitudes[0] ) <= 1.0E-6*fabs( amplitudes[0] ),
                        "Returned peak amplitude moved with the *predicted* measurement length: "
                        << amplitudes[0] << " -> " << amplitudes[i]
                        << ".  It describes the reference spectrum, which did not change." );
  }

  // The continuum legitimately moves a little, and it is worth being clear why rather than
  //  asserting it constant: the coefficients are profiled over the reference block and the
  //  projected sample block *together*, so as the predicted measurement lengthens the sample block
  //  carries more weight and its fixed trial signal pulls the shared continuum down.  That is the
  //  joint model working, not a units error.  Sixteen-fold in predicted length moves it by ~12%.
  for( size_t i = 1; i < continuum_areas.size(); ++i )
  {
    BOOST_CHECK_MESSAGE( fabs( continuum_areas[i] - continuum_areas[0] )
                          <= 0.25*fabs( continuum_areas[0] ),
                        "Continuum under the returned peak moved by more than the joint fit can "
                        "account for: " << continuum_areas[0] << " -> " << continuum_areas[i] );
  }

  // The invariant that actually matters, stated directly: the Gaussian and the continuum it sits on
  //  are in the same exposure.  Their ratio may drift with the joint fit, but it must not *scale*
  //  with the exposure ratio - under the defect it grew by the full 16x across these three cases.
  const double first_fraction = amplitudes.front() / continuum_areas.front();
  const double last_fraction = amplitudes.back() / continuum_areas.back();
  BOOST_CHECK_MESSAGE( (last_fraction / first_fraction) < 1.5,
                      "Peak-to-continuum ratio grew from " << first_fraction << " to "
                      << last_fraction << " over a 16x change in predicted measurement length - "
                      "the Gaussian and its continuum are in different exposures" );
}// BOOST_AUTO_TEST_CASE( BackgroundReferencePeaksMatchTheSpectrumTheyAreDrawnOver )


BOOST_AUTO_TEST_CASE( DeconProfileThreshold )
{
  set_data_dir();

  const shared_ptr<const DetectionLimitCalc::DeconComputeInput> input = make_decon_input( 250.0, 100.0 );

  struct TestCase
  {
    float confidence_level;
    double expected_delta_chi2;
  };

  // z_0.95^2 and z_0.99^2: the one-sided normal thresholds expressed as chi-square(1).
  const TestCase cases[] = { { 0.95f, 2.705543454095404 }, { 0.99f, 5.411894431054342 } };

  for( const TestCase &test : cases )
  {
    const DetectionLimitCalc::DeconActivityOrDistanceLimitResult result =
      DetectionLimitCalc::get_activity_or_distance_limits( test.confidence_level, input, false, 0.0, 600.0, true );
    BOOST_REQUIRE( result.foundLowerCl );
    BOOST_REQUIRE( result.foundUpperCl );
    BOOST_CHECK_SMALL( fabs( ( result.lowerLimitChi2 - result.overallBestChi2 ) - test.expected_delta_chi2 ), 0.02 );
    BOOST_CHECK_SMALL( fabs( ( result.upperLimitChi2 - result.overallBestChi2 ) - test.expected_delta_chi2 ), 0.02 );
  }

  // A central interval is a different product from a one-sided bound and uses a different
  // threshold: quantile(chi2[1], CL) rather than quantile(chi2[1], 2*CL-1).  Taking the two roots
  // of the one-sided threshold would give roughly 2*CL-1 central coverage, i.e. about 90% for a
  // requested 95%, which is why the two are kept explicitly distinct.
  // Tolerances are in percent.  Boost's inverse incomplete gamma is good to roughly 1E-7 relative,
  // which is many orders of magnitude tighter than the 0.025 tolerance the root search itself uses.
  BOOST_CHECK_CLOSE( DetectionLimitCalc::decon_limit_delta( 0.95, DetectionLimitCalc::DeconLimitType::OneSidedUpperLimit ),
                     2.705543454095404, 1.0E-4 );
  BOOST_CHECK_CLOSE( DetectionLimitCalc::decon_limit_delta( 0.95, DetectionLimitCalc::DeconLimitType::CentralInterval ),
                     3.841458820694124, 1.0E-4 );
  BOOST_CHECK_THROW( DetectionLimitCalc::decon_limit_delta( 0.5, DetectionLimitCalc::DeconLimitType::CentralInterval ),
                     runtime_error );

  {
    DetectionLimitCalc::DeconActivityOrDistanceLimitResult central;
    {
      ScopedCoutSilence silence;
      central = DetectionLimitCalc::get_activity_or_distance_limits(
        0.95f, input, false, 0.0, 600.0, true, DetectionLimitCalc::DeconLimitType::CentralInterval );
    }
    BOOST_REQUIRE( central.foundLowerCl );
    BOOST_REQUIRE( central.foundUpperCl );
    BOOST_CHECK_CLOSE( central.limitDelta, 3.841458820694124, 1.0E-4 );
    BOOST_CHECK_SMALL( fabs( ( central.lowerLimitChi2 - central.overallBestChi2 ) - 3.841458820694124 ), 0.02 );
    BOOST_CHECK_SMALL( fabs( ( central.upperLimitChi2 - central.overallBestChi2 ) - 3.841458820694124 ), 0.02 );

    // A central interval is strictly wider than the pair of one-sided roots at the same CL.
    DetectionLimitCalc::DeconActivityOrDistanceLimitResult one_sided;
    {
      ScopedCoutSilence silence;
      one_sided = DetectionLimitCalc::get_activity_or_distance_limits( 0.95f, input, false, 0.0, 600.0, true );
    }
    BOOST_REQUIRE( one_sided.foundUpperCl );
    BOOST_CHECK( central.upperLimit > one_sided.upperLimit );
    BOOST_CHECK( central.lowerLimit < one_sided.lowerLimit );
  }

  BOOST_CHECK_THROW( DetectionLimitCalc::get_activity_or_distance_limits( 0.5f, input, false, 0.0, 600.0, true ),
                     runtime_error );
  BOOST_CHECK_THROW( DetectionLimitCalc::get_activity_or_distance_limits( 1.0f, input, false, 0.0, 600.0, true ),
                     runtime_error );
  BOOST_CHECK_THROW( DetectionLimitCalc::get_activity_or_distance_limits(
                       numeric_limits<float>::quiet_NaN(), input, false, 0.0, 600.0, true ),
                     runtime_error );
}// BOOST_AUTO_TEST_CASE( DeconProfileThreshold )

BOOST_AUTO_TEST_CASE( DeconUnbracketedScanIsAResultNotAnAbort )
{
  set_data_dir();

  // A search range that is far too small no longer loses the limit: the activity scan grows its
  // upper bound geometrically until the profile clears the threshold.  Callers seed this range
  // from a Currie estimate times an arbitrary factor of 50, so it cannot be relied on.
  const shared_ptr<const DetectionLimitCalc::DeconComputeInput> detected = make_decon_input( 250.0, 100.0 );
  DetectionLimitCalc::DeconActivityOrDistanceLimitResult result;
  {
    ScopedCoutSilence silence;
    BOOST_REQUIRE_NO_THROW(
      result = DetectionLimitCalc::get_activity_or_distance_limits( 0.95f, detected, false, 0.0, 260.0, true ) );
  }
  BOOST_CHECK( result.foundLowerCl );
  BOOST_CHECK_MESSAGE( result.foundUpperCl,
                       "the scan did not expand its range to bracket the limit: " << result.errorMessage );
  BOOST_CHECK_MESSAGE( result.upperLimit > 260.0,
                       "expected the limit to lie beyond the requested maximum, got " << result.upperLimit );
  BOOST_CHECK( result.errorMessage.empty() );

  // Answering outside the range that was asked for has to be visible, and `maxSearchValue` must
  // describe the range actually searched rather than the one requested.
  BOOST_CHECK_MESSAGE( !result.warnings.empty(), "extending the search range was not reported" );
  BOOST_CHECK_MESSAGE( result.maxSearchValue > 260.0,
                       "maxSearchValue still reports the un-extended range: " << result.maxSearchValue );
  BOOST_CHECK( result.maxSearchValue >= result.upperLimit );

  const shared_ptr<const DetectionLimitCalc::DeconComputeInput> boundary = make_decon_input( 0.0, 100.0 );
  {
    ScopedCoutSilence silence;
    BOOST_REQUIRE_NO_THROW(
      result = DetectionLimitCalc::get_activity_or_distance_limits( 0.95f, boundary, false, 0.0, 1.0, true ) );
  }
  BOOST_CHECK( !result.foundLowerCl );
  BOOST_CHECK_MESSAGE( result.foundUpperCl, "zero-signal upper limit: " << result.errorMessage );

  // A profile that genuinely carries no information about the activity - here because the source
  // produces essentially no counts - still cannot be bracketed at any range, and must come back as
  // an explicit error result rather than an assertion or an invented number.
  DetectionLimitCalc::DeconComputeInput insensitive = *boundary;
  insensitive.roi_info.front().peak_infos.front().counts_per_bq_into_4pi = 1.0E-30;
  {
    ScopedCoutSilence silence;
    BOOST_REQUIRE_NO_THROW( result = DetectionLimitCalc::get_activity_or_distance_limits(
                              0.95f, make_shared<const DetectionLimitCalc::DeconComputeInput>( insensitive ), false,
                              0.0, 1.0, true ) );
  }
  BOOST_CHECK( !result.foundUpperCl );
  BOOST_CHECK( !result.errorMessage.empty() );
}// BOOST_AUTO_TEST_CASE( DeconUnbracketedScanIsAResultNotAnAbort )

BOOST_AUTO_TEST_CASE( DeconDofBookkeeping )
{
  set_data_dir();

  shared_ptr<const DetectionLimitCalc::DeconComputeInput> floating =
    make_decon_input( 50.0, 100.0, DetectionLimitCalc::DeconContinuumNorm::Floating );
  shared_ptr<const DetectionLimitCalc::DeconComputeInput> full_range =
    make_decon_input( 50.0, 100.0, DetectionLimitCalc::DeconContinuumNorm::FixedByFullRange );
  shared_ptr<const DetectionLimitCalc::DeconComputeInput> edges =
    make_decon_input( 50.0, 100.0, DetectionLimitCalc::DeconContinuumNorm::FixedByEdges );

  DetectionLimitCalc::DeconComputeInput floating_input = *floating;
  DetectionLimitCalc::DeconComputeInput full_range_input = *full_range;
  DetectionLimitCalc::DeconComputeInput edges_input = *edges;
  floating_input.activity = full_range_input.activity = edges_input.activity = 50.0;

  const DetectionLimitCalc::DeconComputeResults floating_result =
    DetectionLimitCalc::decon_compute_peaks( floating_input );
  const DetectionLimitCalc::DeconComputeResults full_range_result =
    DetectionLimitCalc::decon_compute_peaks( full_range_input );
  const DetectionLimitCalc::DeconComputeResults edges_result = DetectionLimitCalc::decon_compute_peaks( edges_input );

  // Floating and FixedByFullRange score the same channels and both spend the two continuum
  // parameters, so they agree.  FixedByEdges scores those channels *plus* its side channels, which
  // are ordinary Poisson observations in the joint likelihood - 4 below and 4 above here - so it
  // reports eight more.  (Under the earlier Gaussian-constraint construction the sidebands were
  // compressed into a two-term penalty and it reported only two more.)
  BOOST_CHECK_EQUAL( floating_result.num_degree_of_freedom, full_range_result.num_degree_of_freedom );
  BOOST_CHECK_EQUAL( edges_result.num_degree_of_freedom, floating_result.num_degree_of_freedom + 8 );

  DetectionLimitCalc::DeconComputeInput second_roi_input = edges_input;
  DetectionLimitCalc::DeconRoiInfo &second_roi = second_roi_input.roi_info.front();
  second_roi.roi_start = 890.0f;
  second_roi.roi_end = 912.0f;
  second_roi.peak_infos.front().energy = 901.0f;

  const DetectionLimitCalc::DeconComputeResults second_result =
    DetectionLimitCalc::decon_compute_peaks( second_roi_input );
  DetectionLimitCalc::DeconComputeInput combined_input = edges_input;
  combined_input.roi_info.push_back( second_roi );
  const DetectionLimitCalc::DeconComputeResults combined_result =
    DetectionLimitCalc::decon_compute_peaks( combined_input );

  BOOST_REQUIRE_EQUAL( combined_result.fit_peaks.size(), 2 );
  BOOST_REQUIRE_EQUAL( edges_result.fit_peaks.size(), 1 );
  BOOST_REQUIRE_EQUAL( second_result.fit_peaks.size(), 1 );
  BOOST_CHECK( combined_result.fit_peaks[0].continuum() != combined_result.fit_peaks[1].continuum() );

  // Regions of interest come back ordered by energy, because they are resolved onto channels and
  // sorted before overlapping ones are combined.  Match on the peak mean rather than on position,
  // so this stays meaningful regardless of the order they were supplied in.
  const auto peak_at = []( const vector<PeakDef> &peaks, const double energy ) -> const PeakDef & {
    for( const PeakDef &p : peaks )
    {
      if( fabs( p.mean() - energy ) < 0.5 )
        return p;
    }
    BOOST_REQUIRE_MESSAGE( false, "no fit peak near " << energy << " keV" );
    return peaks.front();
  };

  BOOST_CHECK_SMALL( fabs( peak_at( combined_result.fit_peaks, 1173.0 ).chi2dof() -
                           peak_at( edges_result.fit_peaks, 1173.0 ).chi2dof() ),
                     1.0E-8 );
  BOOST_CHECK_SMALL( fabs( peak_at( combined_result.fit_peaks, 901.0 ).chi2dof() -
                           peak_at( second_result.fit_peaks, 901.0 ).chi2dof() ),
                     1.0E-8 );
}// BOOST_AUTO_TEST_CASE( DeconDofBookkeeping )

BOOST_AUTO_TEST_CASE( DeconMultiRoiStudy )
{
  set_data_dir();

  const shared_ptr<const DetectionLimitCalc::DeconComputeInput> input = make_decon_input( 0.0, 100.0 );
  DetectionLimitCalc::DeconComputeInput independent = *input;
  DetectionLimitCalc::DeconRoiInfo second = independent.roi_info.front();
  second.roi_start = 890.0f;
  second.roi_end = 912.0f;
  second.peak_infos.front().energy = 901.0f;
  independent.roi_info.push_back( second );

  // Split one line's yield across two identical, exactly overlapping ROIs.  The *signal model* is
  // then unchanged from the single-ROI case, so any difference in the limit could only come from
  // counting the shared channels' data twice.  Overlaps are combined into one joint region, so
  // there must be no difference at all.
  DetectionLimitCalc::DeconComputeInput overlapping = *input;
  overlapping.roi_info.front().peak_infos.front().counts_per_bq_into_4pi *= 0.5;
  overlapping.roi_info.push_back( overlapping.roi_info.front() );

  // A realistic partial overlap: a second line close enough that the two ISO-recommended
  // +-1.25 FWHM regions share channels, which is what the GUI produces for any nearby doublet.
  DetectionLimitCalc::DeconComputeInput partial = *input;
  DetectionLimitCalc::DeconRoiInfo neighbor = partial.roi_info.front();
  neighbor.roi_start = 1176.0f;
  neighbor.roi_end = 1198.0f;
  neighbor.peak_infos.front().energy = 1187.0f;
  partial.roi_info.push_back( neighbor );

  DetectionLimitCalc::DeconActivityOrDistanceLimitResult single_result;
  DetectionLimitCalc::DeconActivityOrDistanceLimitResult independent_result;
  DetectionLimitCalc::DeconActivityOrDistanceLimitResult overlapping_result;
  {
    ScopedCoutSilence silence;
    single_result = DetectionLimitCalc::get_activity_or_distance_limits( 0.95f, input, false, 0.0, 500.0, true );
    independent_result = DetectionLimitCalc::get_activity_or_distance_limits(
      0.95f, make_shared<const DetectionLimitCalc::DeconComputeInput>( independent ), false, 0.0, 500.0, true );
    overlapping_result = DetectionLimitCalc::get_activity_or_distance_limits(
      0.95f, make_shared<const DetectionLimitCalc::DeconComputeInput>( overlapping ), false, 0.0, 500.0, true );
  }

  BOOST_REQUIRE( single_result.foundUpperCl );
  BOOST_REQUIRE( independent_result.foundUpperCl );
  BOOST_REQUIRE( overlapping_result.foundUpperCl );
  BOOST_CHECK_CLOSE( independent_result.upperLimit / single_result.upperLimit, 1.0 / sqrt( 2.0 ), 4.0 );

  // No manufactured information: splitting one line across two coincident regions must give
  // exactly the single-region answer.  Before overlapping regions were combined, this produced the
  // full sqrt(2) improvement purely from double-counted channels.
  BOOST_CHECK_CLOSE( overlapping_result.upperLimit / single_result.upperLimit, 1.0, 0.5 );

  // The merge is reported rather than done silently, and the two regions become one continuum.
  const DetectionLimitCalc::DeconComputeResults partial_results =
    DetectionLimitCalc::decon_compute_peaks( partial );
  BOOST_CHECK_MESSAGE( !partial_results.warnings.empty(),
                       "combining overlapping regions of interest was not reported" );
  BOOST_REQUIRE_EQUAL( partial_results.fit_peaks.size(), 2 );
  BOOST_CHECK( partial_results.fit_peaks[0].continuum() == partial_results.fit_peaks[1].continuum() );

  cout << "DECON_MULTI_ROI,single_upper,independent_two_roi_upper,overlapping_two_roi_upper\n"
       << "DECON_MULTI_ROI," << single_result.upperLimit << ',' << independent_result.upperLimit << ','
       << overlapping_result.upperLimit << '\n';

  const shared_ptr<const DetectionLimitCalc::DeconComputeInput> one_observed_line = make_decon_input( 100.0, 100.0 );
  DetectionLimitCalc::DeconComputeInput wrong_yield = *one_observed_line;
  wrong_yield.roi_info.push_back( second );
  DetectionLimitCalc::DeconActivityOrDistanceLimitResult one_line_result;
  DetectionLimitCalc::DeconActivityOrDistanceLimitResult wrong_yield_result;
  {
    ScopedCoutSilence silence;
    one_line_result =
      DetectionLimitCalc::get_activity_or_distance_limits( 0.95f, one_observed_line, false, 0.0, 500.0, true );
    wrong_yield_result = DetectionLimitCalc::get_activity_or_distance_limits(
      0.95f, make_shared<const DetectionLimitCalc::DeconComputeInput>( wrong_yield ), false, 0.0, 500.0, true );
  }
  BOOST_REQUIRE( one_line_result.foundUpperCl );
  BOOST_REQUIRE( wrong_yield_result.foundUpperCl );
  BOOST_CHECK( wrong_yield_result.overallBestQuantity < one_line_result.overallBestQuantity );
  BOOST_CHECK( wrong_yield_result.upperLimit < one_line_result.upperLimit );
  cout << "DECON_BAD_LINE,one_line_best,one_line_upper,missing_second_line_best,"
          "missing_second_line_upper\n"
       << "DECON_BAD_LINE," << one_line_result.overallBestQuantity << ',' << one_line_result.upperLimit << ','
       << wrong_yield_result.overallBestQuantity << ',' << wrong_yield_result.upperLimit << '\n';
}// BOOST_AUTO_TEST_CASE( DeconMultiRoiStudy )

BOOST_AUTO_TEST_CASE( DeconCurrieCleanSpectrumGrid )
{
  set_data_dir();

  struct RatioResult
  {
    double continuum;
    double signal_sigma;
    double currie_limit;
    double decon_limit;
    double ratio;
  };

  vector<RatioResult> results;
  const double continua[] = { 0.1, 1.0, 10.0, 100.0 };
  const double signal_sigmas[] = { 0.0, 1.0, 3.0 };

  {
    ScopedCoutSilence silence;
    for( const double continuum : continua )
    {
      const shared_ptr<const DetectionLimitCalc::DeconComputeInput> background = make_decon_input( 0.0, continuum );
      const pair<size_t, size_t> roi_channels = DetectionLimitCalc::round_roi_to_channels(
        background->measurement, background->roi_info.front().roi_start, background->roi_info.front().roi_end );
      const double roi_background = continuum * ( 1 + roi_channels.second - roi_channels.first );

      for( const double signal_sigma : signal_sigmas )
      {
        const double signal_counts = signal_sigma * sqrt( roi_background );
        const shared_ptr<const DetectionLimitCalc::DeconComputeInput> input =
          make_decon_input( signal_counts, continuum );
        const DetectionLimitCalc::CurrieMdaResult currie = currie_for_decon_input( input );
        const double max_quantity =
          (std::max)( 25.0, signal_counts + 50.0 * sqrt( roi_background + signal_counts + 1.0 ) );
        const DetectionLimitCalc::DeconActivityOrDistanceLimitResult decon =
          DetectionLimitCalc::get_activity_or_distance_limits( 0.95f, input, false, 0.0, max_quantity, true );
        BOOST_REQUIRE( decon.foundUpperCl );
        BOOST_REQUIRE_GT( currie.upper_limit, 0.0 );
        const double ratio = decon.upperLimit / currie.upper_limit;
        results.push_back( { continuum, signal_sigma, currie.upper_limit, decon.upperLimit, ratio } );
        BOOST_CHECK_MESSAGE( ( ratio > 0.25 ) && ( ratio < 4.0 ), "Implausible clean-spectrum ratio " << ratio );
      }
    }
  }

  cout << "DECON_CURRIE_GRID,continuum_counts_per_channel,signal_over_sqrt_roi_background,"
          "currie_upper_counts,decon_upper_counts,decon_over_currie\n";
  for( const RatioResult &result : results )
  {
    cout << "DECON_CURRIE_GRID," << result.continuum << ',' << result.signal_sigma << ',' << result.currie_limit << ','
         << result.decon_limit << ',' << result.ratio << '\n';
  }
}// BOOST_AUTO_TEST_CASE( DeconCurrieCleanSpectrumGrid )

/** `currie_check_for_peak` on a spectrum whose channels are not uniformly wide.

 This is the entry point `ShieldingSourceFitCalc::fit_model` calls for every peak, so its numbers
 reach the activity/shielding fit report and, through it, the batch output - see
 `BatchActivity.cpp` and `BatchInfoLog.cpp`.  It had no test of its own, and the batch tests that do
 exist build a uniform 2 keV/channel polynomial with no deviation pairs, which is exactly the case
 where the width-weighted continuum fit reduces to the mean of the two side-band densities it
 replaced.  So nothing exercised this path on the binning the fix actually changes.

 The Ba-133 spectrum's deviation pairs make the two side bands span unequal energy, which is what
 used to bias the continuum - by -0.9 sigma at the 81 keV line, in the direction that manufactures
 signal that is not there.  `Detected` is decided by `source_counts > decision_threshold`, so that
 bias moved detection decisions and not just a displayed number.  \sa CurrieContinuumNonUniformBinning
 */
BOOST_AUTO_TEST_CASE( CurrieCheckForPeakNonUniformBinning )
{
  set_data_dir();

  const string path = SpecUtils::append_path( g_test_file_dir, "PeakFitLM/Ba133_Unshielded.n42" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( path ), "Ba-133 test spectrum not at '" << path << "'" );

  SpecMeas file;
  BOOST_REQUIRE_MESSAGE( file.load_N42_file( path ), "Failed loading '" << path << "'" );
  const set<set<int>> peak_sample_sets = file.sampleNumsWithAutomatedSearchPeaks();
  BOOST_REQUIRE( !peak_sample_sets.empty() );
  const set<int> &peak_samples = *begin( peak_sample_sets );
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks
                                                    = file.automatedSearchPeaks( peak_samples );
  BOOST_REQUIRE( peaks && !peaks->empty() );

  shared_ptr<const SpecUtils::Measurement> measurement;
  for( const shared_ptr<const SpecUtils::Measurement> &candidate : file.measurements() )
  {
    if( candidate && (candidate->num_gamma_channels() > 16)
       && peak_samples.count( candidate->sample_number() ) )
    {
      measurement = candidate;
      break;
    }
  }
  BOOST_REQUIRE( measurement );

  // The premise of the case: without deviation pairs this spectrum would not distinguish the two
  //  continuum estimators at all, and every assertion below would pass vacuously.
  const shared_ptr<const SpecUtils::EnergyCalibration> cal = measurement->energy_calibration();
  BOOST_REQUIRE( cal && cal->valid() );
  BOOST_REQUIRE_MESSAGE( !cal->deviation_pairs().empty(),
                        "this spectrum no longer has deviation pairs, so it does not test"
                        " non-uniform binning and this case is vacuous" );

  DetectionLimitCalc::PeakCurrieCheckOptions options;
  options.confidence_level = 0.95;
  options.num_side_channels = 4;
  options.roi_num_fwhm = 2.5;

  size_t checked = 0;
  cout << "CURRIE_CHECK_PEAK,energy_keV,result_type,continuum_counts,continuum_from_eqn,"
          "rel_diff,source_counts,decision_threshold\n";

  for( const shared_ptr<const PeakDef> &peak : *peaks )
  {
    if( !peak || !peak->gausPeak() || (peak->fwhm() <= 0.0) )
      continue;

    const DetectionLimitCalc::PeakCurrieCheck check
                    = DetectionLimitCalc::currie_check_for_peak( *peak, measurement, options, true );

    BOOST_REQUIRE_MESSAGE( check.computed,
                          "no Currie check at " << peak->mean() << " keV: " << check.error_message );

    const DetectionLimitCalc::CurrieMdaResult &res = check.result;

    // Integrate the continuum line this check reports over the peak region it reports.  Under the
    //  old estimator these two disagreed on this spectrum - the drawn continuum and the counted
    //  continuum were different numbers - and the agreement is now an algebraic identity, so it
    //  holds on any binning rather than only on uniform binning.
    //
    //  Over the region's *channel edges*, not `input.roi_lower_energy`/`roi_upper_energy`: the
    //  requested ROI is snapped outwards to whole channels, and on this spectrum's widening
    //  channels the two differ by enough to fail this check by up to 14% if you use the request.
    const double lower = measurement->gamma_channel_lower( res.first_peak_region_channel );
    const double upper = measurement->gamma_channel_upper( res.last_peak_region_channel );
    const double x0 = lower - res.input.gamma_energy;
    const double x1 = upper - res.input.gamma_energy;
    const double from_eqn = res.continuum_eqn[0]*(x1 - x0)
                              + 0.5*res.continuum_eqn[1]*(x1*x1 - x0*x0);

    const double scale = (std::max)( 1.0, static_cast<double>(res.estimated_peak_continuum_counts) );
    const double rel_diff = fabs( from_eqn - res.estimated_peak_continuum_counts ) / scale;

    cout << "CURRIE_CHECK_PEAK," << peak->mean() << ',' << static_cast<int>(check.result_type) << ','
         << res.estimated_peak_continuum_counts << ',' << from_eqn << ',' << rel_diff << ','
         << res.source_counts << ',' << res.decision_threshold << '\n';

    BOOST_CHECK_MESSAGE( rel_diff < 1.0E-5,
                        "at " << peak->mean() << " keV the reported continuum counts ("
                        << res.estimated_peak_continuum_counts << ") disagree with the integral of"
                        " the reported continuum equation (" << from_eqn << ")" );

    // A continuum estimate cannot be negative, and the detection decision has to be the one the
    //  reported numbers imply - not an independently computed flag that can drift from them.
    BOOST_CHECK_GE( res.estimated_peak_continuum_counts, 0.0f );

    if( check.result_type == DetectionLimitCalc::PeakCurrieCheck::ResultType::Detected )
      BOOST_CHECK_MESSAGE( res.source_counts > res.decision_threshold,
                          "at " << peak->mean() << " keV the check says Detected, but its own"
                          " source counts " << res.source_counts << " do not exceed the decision"
                          " threshold " << res.decision_threshold );
    else if( check.result_type == DetectionLimitCalc::PeakCurrieCheck::ResultType::NotDetected )
      BOOST_CHECK_MESSAGE( res.source_counts <= res.decision_threshold,
                          "at " << peak->mean() << " keV the check says NotDetected, but its own"
                          " source counts " << res.source_counts << " exceed the decision"
                          " threshold " << res.decision_threshold );

    ++checked;
  }//for( loop over the fit peaks )

  BOOST_REQUIRE_MESSAGE( checked >= 4,
                        "only " << checked << " peaks were checked; this spectrum should offer"
                        " several, so the case is not covering what it claims" );
}// BOOST_AUTO_TEST_CASE( CurrieCheckForPeakNonUniformBinning )


BOOST_AUTO_TEST_CASE( DeconCurrieBundledSpectrumGrid )
{
  set_data_dir();

  // This used to bail out of Debug builds: `currie_mda_calc(...)` estimated the peak-region
  //  continuum as the mean of the two side-band densities, which its own developer assertion
  //  contradicts whenever the bands span unequal energy - and the Ba-133 spectrum's deviation pairs
  //  make them unequal.  The estimate is now a width-weighted least-squares fit, so the assertion
  //  is an identity and the suite's only real spectrum runs in both configurations.
  //  \sa CurrieContinuumNonUniformBinning

  const string path = SpecUtils::append_path( g_test_file_dir, "PeakFitLM/Ba133_Unshielded.n42" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( path ), "Ba-133 test spectrum not at '" << path << "'" );

  SpecMeas file;
  BOOST_REQUIRE_MESSAGE( file.load_N42_file( path ), "Failed loading '" << path << "'" );
  const set<set<int>> peak_sample_sets = file.sampleNumsWithAutomatedSearchPeaks();
  BOOST_REQUIRE( !peak_sample_sets.empty() );
  const set<int> &peak_samples = *begin( peak_sample_sets );
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks = file.automatedSearchPeaks( peak_samples );
  BOOST_REQUIRE( peaks && !peaks->empty() );

  shared_ptr<const SpecUtils::Measurement> measurement;
  for( const shared_ptr<const SpecUtils::Measurement> &candidate : file.measurements() )
  {
    if( candidate && ( candidate->num_gamma_channels() > 16 ) && peak_samples.count( candidate->sample_number() ) )
    {
      measurement = candidate;
      break;
    }
  }
  BOOST_REQUIRE( measurement );

  shared_ptr<DetectorPeakResponse> drf = make_shared<DetectorPeakResponse>();
  drf->setIntrinsicEfficiencyFormula( "1.0",
                                      2.54f * PhysicalUnits::cm,
                                      PhysicalUnits::keV,
                                      0.0f,
                                      0.0f,
                                      DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct );
  drf->setFwhmCoefficients( vector<float>{ 64.0f, 0.0f }, DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial );

  size_t num_ratios = 0;
  cout << "DECON_CURRIE_BUNDLED,energy_keV,continuum_counts_per_channel,"
          "signal_over_sqrt_roi_background,currie_upper_counts,decon_upper_counts,"
          "decon_over_currie\n";
  for( const shared_ptr<const PeakDef> &peak : *peaks )
  {
    if( !peak || !peak->gausPeak() || ( peak->skewType() != PeakDef::SkewType::NoSkew ) || ( peak->fwhm() <= 0.0 ) ||
        !peak->continuum() || ( peak->continuum()->type() != PeakContinuum::OffsetType::Linear ) )
      continue;

    bool isolated = true;
    for( const shared_ptr<const PeakDef> &other : *peaks )
    {
      if( other && ( other != peak ) &&
          ( fabs( other->mean() - peak->mean() ) < 3.0 * (std::max)( peak->fwhm(), other->fwhm() ) ) )
      {
        isolated = false;
        break;
      }
    }
    if( !isolated )
      continue;

    DetectionLimitCalc::DeconRoiInfo roi;
    roi.roi_start = static_cast<float>( peak->mean() - 1.25 * peak->fwhm() );
    roi.roi_end = static_cast<float>( peak->mean() + 1.25 * peak->fwhm() );
    roi.continuum_type = PeakContinuum::OffsetType::Linear;
    roi.cont_norm_method = DetectionLimitCalc::DeconContinuumNorm::Floating;
    roi.num_lower_side_channels = 4;
    roi.num_upper_side_channels = 4;

    DetectionLimitCalc::DeconRoiInfo::PeakInfo peak_info;
    peak_info.energy = static_cast<float>( peak->mean() );
    peak_info.fwhm = static_cast<float>( peak->fwhm() );
    peak_info.counts_per_bq_into_4pi = 1.0;
    roi.peak_infos.push_back( peak_info );

    shared_ptr<DetectionLimitCalc::DeconComputeInput> input = make_shared<DetectionLimitCalc::DeconComputeInput>();
    input->activity = 0.0;
    input->distance = 0.0;
    input->include_air_attenuation = false;
    input->shielding_thickness = 0.0;
    input->drf = drf;
    input->measurement = measurement;
    input->roi_info.push_back( roi );

    DetectionLimitCalc::CurrieMdaResult currie;
    try
    {
      currie = currie_for_decon_input( input );
    } catch( std::exception &error )
    {
      cout << "DECON_CURRIE_BUNDLED_SKIP," << peak->mean() << ",currie_error," << error.what() << '\n';
      continue;
    }
    if( currie.upper_limit <= 0.0f )
    {
      cout << "DECON_CURRIE_BUNDLED_SKIP," << peak->mean() << ",invalid_currie_upper\n";
      continue;
    }

    const double multiple = 50.0;
    double minimum = 0.0;
    double maximum = 0.0;
    if( currie.source_counts > currie.decision_threshold )
    {
      const double lower_difference = fabs( currie.source_counts - currie.lower_limit );
      const double upper_difference = fabs( currie.upper_limit - currie.source_counts );
      minimum = (std::max)( 0.0, currie.source_counts - multiple * (std::max)( lower_difference, upper_difference ) );
      maximum = currie.source_counts + multiple * (std::max)( lower_difference, upper_difference );
    } else
    {
      maximum = multiple * currie.upper_limit;
    }
    maximum = (std::max)( 1.0, maximum );

    DetectionLimitCalc::DeconActivityOrDistanceLimitResult decon;
    {
      ScopedCoutSilence silence;
      decon = DetectionLimitCalc::get_activity_or_distance_limits( 0.95f, input, false, minimum, maximum, true );
    }
    if( !decon.foundUpperCl )
    {
      cout << "DECON_CURRIE_BUNDLED_SKIP," << peak->mean() << ",unbracketed_decon_upper\n";
      continue;
    }

    const pair<size_t, size_t> roi_channels =
      DetectionLimitCalc::round_roi_to_channels( measurement, roi.roi_start, roi.roi_end );
    double roi_counts = 0.0;
    for( size_t channel = roi_channels.first; channel <= roi_channels.second; ++channel )
      roi_counts += measurement->gamma_channel_content( channel );
    const size_t num_roi_channels = 1 + roi_channels.second - roi_channels.first;
    const double continuum =
      peak->continuum()->offset_integral( roi.roi_start, roi.roi_end, measurement, nullptr, 0 ) / num_roi_channels;
    const double signal_sigma =
      ( roi_counts - continuum * num_roi_channels ) / sqrt( (std::max)( 1.0, continuum * num_roi_channels ) );
    const double ratio = decon.upperLimit / currie.upper_limit;
    cout << "DECON_CURRIE_BUNDLED," << peak->mean() << ',' << continuum << ',' << signal_sigma << ','
         << currie.upper_limit << ',' << decon.upperLimit << ',' << ratio << '\n';
    ++num_ratios;
  }

  BOOST_REQUIRE_MESSAGE( num_ratios >= 2, "Expected at least two isolated, unskewed peaks in the bundled spectrum" );
}// BOOST_AUTO_TEST_CASE( DeconCurrieBundledSpectrumGrid )

/** Does a limit predicted for a *different* dwell match what that dwell actually delivers?

 Every other Monte Carlo here runs on a synthetic 1024-channel, flat-continuum, single-Gaussian
 fixture.  This one takes its rate model, energy calibration and dead-time fraction from a real
 spectrum, then draws two *independent* measurements from it - a reference of `T_ref` and an actual
 future measurement of `T_s` - and asks what each projection would have told the user beforehand,
 against what the future measurement then gave.

 Two projections are compared, because the tools offer both:
   - Currie: `scale_spectrum_for_dwell` multiplies the reference's counts up to `T_s` and takes an
     ordinary limit.  The scaled counts are then treated as Poisson with variance equal to the
     scaled count, which asserts the reference's rate is known exactly.
   - Deconvolution `BackgroundReference`: the reference and the projected sample share one
     continuum in a joint likelihood, so the reference's own counting error propagates.

 Taking the real spectrum's observed counts as the rate model makes the "truth" only as good as
 that spectrum's own statistics, but that error is common to all arms, so the comparison between
 them is unaffected.

 Two things beyond the medians are measured here, both of which the first revision of this case
 could not see:

 **D-17 / D-20 - what the background-reference optimism is made of.**  Comparing the prediction to
 a limit taken on the future measurement *alone* mixes two effects: the expected-counts step evaluating at
 unfluctuated data, and the joint model legitimately knowing more than a sample-only analysis.
 `DeconComputeInput::observed_sample` re-scores the same two-block likelihood with the future
 measurement's real counts, which splits the ratio into those two factors:

     prediction/sample-only  ==  (prediction/joint) * (joint/sample-only)

 the first factor being the prediction calibration D-17 asks for, the second the information the
 reference genuinely adds.

 **D-21 - the spread of the Currie projection, not just its median.**  Projecting `k = T_s/T_ref`
 makes the prediction inherit the *reference's* counting noise: the projected limit's relative
 scatter goes as `1/sqrt(counts at T_ref)` where the future limit's goes as `1/sqrt(counts at T_s)`,
 so their ratio should be about `sqrt(k)` and the honest predictive spread is understated by about
 `sqrt(1 + k)`.  Both are measured against `sqrt(k)` below rather than asserted.
 */
BOOST_AUTO_TEST_CASE( RealSpectrumDwellProjectionAccuracy )
{
  using namespace DetectionLimitCalc;

  set_data_dir();

  // A real 11.5-hour background: non-trivial dead time (41610.3 s live against 42101.5 s real), a
  //  real energy calibration, and a real continuum shape.
  const string path = SpecUtils::append_path( g_test_file_dir,
                                             "FitPeaksForSource/trinitite_sample_b_background.n42" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( path ), "Background spectrum not at '" << path << "'" );

  SpecMeas file;
  BOOST_REQUIRE_MESSAGE( file.load_N42_file( path ), "Failed loading '" << path << "'" );

  shared_ptr<const SpecUtils::Measurement> truth_rate;
  for( const shared_ptr<const SpecUtils::Measurement> &candidate : file.measurements() )
  {
    if( candidate && (candidate->num_gamma_channels() > 1024) && (candidate->live_time() > 0.0f) )
    {
      truth_rate = candidate;
      break;
    }
  }
  BOOST_REQUIRE( truth_rate );
  BOOST_REQUIRE( truth_rate->real_time() > truth_rate->live_time() );  // it really does have dead time

  const double source_live = truth_rate->live_time();
  const double live_over_real = truth_rate->live_time() / truth_rate->real_time();

  shared_ptr<DetectorPeakResponse> drf = make_shared<DetectorPeakResponse>();
  drf->setIntrinsicEfficiencyFormula( "1.0", 2.54f*PhysicalUnits::cm, PhysicalUnits::keV,
                                      0.0f, 0.0f,
                                      DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct );
  drf->setFwhmCoefficients( vector<float>{ 4.0f, 0.02f },
                            DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial );

  // A clean stretch of continuum away from the strong natural lines - the situation a background
  //  reference is actually used in.
  const float peak_energy = 1000.0f;
  const float fwhm = static_cast<float>( drf->peakResolutionFWHM( peak_energy ) );
  BOOST_REQUIRE( fwhm > 0.0f );

  // Draws an independent measurement of `live_seconds` from the real spectrum's rate, keeping the
  //  source's dead-time fraction so the live/real distinction is exercised end to end.
  const auto draw = [&]( const double live_seconds, mt19937 &generator )
        -> shared_ptr<const SpecUtils::Measurement> {
    const shared_ptr<const vector<float>> rate = truth_rate->gamma_counts();
    shared_ptr<vector<float>> counts = make_shared<vector<float>>( rate->size(), 0.0f );
    const double scale = live_seconds / source_live;
    for( size_t i = 0; i < rate->size(); ++i )
      counts->at(i) = static_cast<float>( study_poisson( (std::max)(0.0, scale*rate->at(i)),
                                                         generator ) );

    shared_ptr<SpecUtils::Measurement> m = make_shared<SpecUtils::Measurement>();
    m->set_gamma_counts( counts, static_cast<float>(live_seconds),
                         static_cast<float>(live_seconds/live_over_real) );
    m->set_energy_calibration( truth_rate->energy_calibration() );
    return m;
  };

  const auto make_input = [&]( const shared_ptr<const SpecUtils::Measurement> &m,
                               const DeconMeasurementModel model, const double sample_live )
        -> shared_ptr<DeconComputeInput> {
    DeconRoiInfo roi;
    roi.roi_start = peak_energy - 1.25f*fwhm;
    roi.roi_end = peak_energy + 1.25f*fwhm;
    roi.continuum_type = PeakContinuum::OffsetType::Linear;
    roi.cont_norm_method = DeconContinuumNorm::Floating;
    roi.num_lower_side_channels = 4;
    roi.num_upper_side_channels = 4;

    DeconRoiInfo::PeakInfo peak;
    peak.energy = peak_energy;
    peak.fwhm = fwhm;
    // Counts per Bq carries the spectrum's own live time, so "activity" reads as counts per second
    //  and limits from different dwells are directly comparable.
    peak.counts_per_bq_into_4pi = m->live_time();
    roi.peak_infos.push_back( peak );

    shared_ptr<DeconComputeInput> input = make_shared<DeconComputeInput>();
    input->activity = 0.0;
    input->distance = 0.0;
    input->include_air_attenuation = false;
    input->shielding_thickness = 0.0;
    input->drf = drf;
    input->measurement = m;
    input->roi_info.push_back( roi );
    input->measurement_model = model;
    input->sample_exposure = (model == DeconMeasurementModel::BackgroundReference)
                               ? sample_live : 0.0;
    return input;
  };

  // The same two-block joint likelihood the background-reference prediction uses, but with the
  //  future measurement's real counts in place of the expected-counts sample block.  This is the estimand
  //  the prediction is supposed to be the median of, and nothing measured it before.
  const auto make_joint_input = [&]( const shared_ptr<const SpecUtils::Measurement> &reference,
                                     const shared_ptr<const SpecUtils::Measurement> &sample )
        -> shared_ptr<DeconComputeInput> {
    const shared_ptr<DeconComputeInput> input
        = make_input( reference, DeconMeasurementModel::BackgroundReference, sample->live_time() );
    input->observed_sample = sample;
    return input;
  };

  const auto limit_of = [&]( const shared_ptr<const DeconComputeInput> &input ) -> double {
    try
    {
      ScopedCoutSilence silence;
      const DeconActivityOrDistanceLimitResult r =
          get_activity_or_distance_limits( 0.95f, input, false, 0.0, 5.0, true );
      return r.foundUpperCl ? r.upperLimit : numeric_limits<double>::quiet_NaN();
    }catch( std::exception & )
    {
      return numeric_limits<double>::quiet_NaN();
    }
  };

  // A Currie limit on `m`, as a rate, so it is directly comparable with the profile limits above.
  //  `detection_limit` and not `upper_limit`, because the detection limit is what the tools quote
  //  when a measurement time is entered, and so is what the projection's spread has to describe.
  //  The two scatter differently: `upper_limit` also carries the observed peak-region excess, while
  //  `detection_limit` is a function of the continuum estimate alone.
  const auto currie_rate_of = [&]( const shared_ptr<const SpecUtils::Measurement> &m ) -> double {
    try
    {
      ScopedCoutSilence silence;
      const CurrieMdaResult r = currie_for_decon_input(
                                  make_input( m, DeconMeasurementModel::CurrentSpectrum, 0.0 ) );
      return (m->live_time() > 0.0f) ? (r.detection_limit / m->live_time())
                                     : numeric_limits<double>::quiet_NaN();
    }catch( std::exception & )
    {
      return numeric_limits<double>::quiet_NaN();
    }
  };

  cout << "DECON_DWELL_PROJECTION,reference_live_s,sample_live_s,ratio,reps,"
          "decon_actual,decon_backref_pred,decon_scaled_pred,currie_actual,currie_scaled_pred,"
          "backref_over_actual,scaled_decon_over_actual,currie_scaled_over_currie_actual\n";

  // D-17: does the expected-counts prediction sit at the median of the *same* joint procedure run on real
  //  sample counts?  And how much of the prediction-vs-sample-only ratio is that, against the
  //  information the reference legitimately adds?
  cout << "DECON_PREDICTION_CALIBRATION,ratio,pairs,predicted,joint_actual,sample_only_actual,"
          "predicted_over_joint,joint_over_sample_only,predicted_over_sample_only,product_check\n";

  // D-21: the projected Currie limit inherits the reference's counting noise, so its scatter should
  //  exceed the future limit's by about sqrt(k) - which makes the honest predictive spread about
  //  sqrt(1+k) wider than the projection implies.
  cout << "DECON_PROJECTION_SPREAD,ratio,pairs,currie_actual_rel_iqr,currie_projected_rel_iqr,"
          "spread_ratio,sqrt_k,implied_spread_understatement\n";

  const auto median_of = []( vector<double> v ) -> double {
    if( v.empty() )
      return numeric_limits<double>::quiet_NaN();
    const size_t mid = v.size()/2;
    nth_element( v.begin(), v.begin() + mid, v.end() );
    return v[mid];
  };

  // Interquartile range over the median: a scatter measure that a handful of failed scans or a
  //  fat tail cannot dominate, which matters at 80 realisations.
  const auto relative_iqr = []( vector<double> v ) -> double {
    if( v.size() < 8 )
      return numeric_limits<double>::quiet_NaN();
    sort( begin(v), end(v) );
    const double q1 = v[v.size()/4];
    const double q2 = v[v.size()/2];
    const double q3 = v[(3*v.size())/4];
    return (q2 > 0.0) ? ((q3 - q1)/q2) : numeric_limits<double>::quiet_NaN();
  };

  const double reference_live = 600.0;
  const size_t reps = 80;

  for( const double ratio : { 0.25, 1.0, 4.0, 16.0 } )
  {
    const double sample_live = ratio * reference_live;
    mt19937 generator( study_cell_seed( "dwell|ratio=" + study_number_token( ratio ) ) );

    vector<double> actual, backref, scaled_decon, currie_actual, currie_scaled, joint;

    // The three-way decomposition has to be taken over the realisations where *all three* arms
    //  gave a limit, or the ratio of medians would be a ratio over different sub-samples.
    vector<double> paired_actual, paired_backref, paired_joint;

    for( size_t rep = 0; rep < reps; ++rep )
    {
      const shared_ptr<const SpecUtils::Measurement> reference = draw( reference_live, generator );
      const shared_ptr<const SpecUtils::Measurement> future = draw( sample_live, generator );

      // What the future measurement actually delivered.
      const double a = limit_of( make_input( future, DeconMeasurementModel::CurrentSpectrum, 0.0 ) );

      // What a background reference predicted for it, beforehand.
      const double b = limit_of( make_input( reference, DeconMeasurementModel::BackgroundReference,
                                             sample_live ) );

      // And what the *same joint model* gives once the future measurement is in hand - the estimand
      //  the expected-counts prediction `b` is meant to be the median of.  \sa D-17
      const double f = limit_of( make_joint_input( reference, future ) );

      // What scaling the reference up to the future dwell predicted, beforehand - by both methods.
      //  The entered time is a REAL time, hence the division by the live/real fraction.
      double c = numeric_limits<double>::quiet_NaN(), e = numeric_limits<double>::quiet_NaN();
      try
      {
        const shared_ptr<const SpecUtils::Measurement> scaled
            = scale_spectrum_for_dwell( reference,
                                        static_cast<float>(sample_live/live_over_real) );
        c = limit_of( make_input( scaled, DeconMeasurementModel::CurrentSpectrum, 0.0 ) );
        e = currie_rate_of( scaled );
      }catch( std::exception & )
      {
      }

      // And the Currie method on the future measurement itself, so the Currie projection has a
      //  like-for-like truth to be judged against rather than being compared to a profile limit.
      const double d = currie_rate_of( future );

      if( !IsNan(a) ) actual.push_back( a );
      if( !IsNan(b) ) backref.push_back( b );
      if( !IsNan(c) ) scaled_decon.push_back( c );
      if( !IsNan(d) ) currie_actual.push_back( d );
      if( !IsNan(e) ) currie_scaled.push_back( e );
      if( !IsNan(f) ) joint.push_back( f );

      if( !IsNan(a) && !IsNan(b) && !IsNan(f) )
      {
        paired_actual.push_back( a );
        paired_backref.push_back( b );
        paired_joint.push_back( f );
      }
    }//for( loop over repetitions )

    BOOST_REQUIRE_MESSAGE( actual.size() > reps/2,
                          "ratio " << ratio << ": only " << actual.size() << " of " << reps
                          << " future measurements gave a limit" );

    const double ma = median_of( actual );
    const double mb = median_of( backref );
    const double mc = median_of( scaled_decon );
    const double md = median_of( currie_actual );
    const double me = median_of( currie_scaled );

    cout << "DECON_DWELL_PROJECTION," << reference_live << ',' << sample_live << ',' << ratio << ','
         << actual.size() << ',' << ma << ',' << mb << ',' << mc << ',' << md << ',' << me << ','
         << (ma > 0.0 ? mb/ma : 0.0) << ',' << (ma > 0.0 ? mc/ma : 0.0) << ','
         << (md > 0.0 ? me/md : 0.0) << '\n';

    // D-17 / D-20: split the prediction-vs-sample-only ratio into the expected-counts step and the joint
    //  model's genuine information advantage.  All three medians over the same realisations.
    BOOST_REQUIRE_MESSAGE( paired_joint.size() > reps/2,
                          "ratio " << ratio << ": only " << paired_joint.size() << " of " << reps
                          << " realisations gave all three of prediction, joint and sample-only" );

    const double pa = median_of( paired_actual );
    const double pb = median_of( paired_backref );
    const double pf = median_of( paired_joint );

    const double predicted_over_joint = (pf > 0.0) ? (pb/pf) : 0.0;
    const double joint_over_sample = (pa > 0.0) ? (pf/pa) : 0.0;
    const double predicted_over_sample = (pa > 0.0) ? (pb/pa) : 0.0;

    cout << "DECON_PREDICTION_CALIBRATION," << ratio << ',' << paired_joint.size() << ','
         << pb << ',' << pf << ',' << pa << ',' << predicted_over_joint << ',' << joint_over_sample
         << ',' << predicted_over_sample << ','
         << (predicted_over_joint*joint_over_sample) << '\n';

    // The decomposition is arithmetic, not a new estimand: the two factors must multiply back to
    //  the ratio already published.  A failure here means the medians were taken over different
    //  sub-samples, not that the physics changed.
    BOOST_CHECK_CLOSE_FRACTION( predicted_over_joint*joint_over_sample, predicted_over_sample, 1.0E-9 );

    // The joint model sees the sample *and* the reference, so it cannot be the looser of the two by
    //  any margin worth speaking of.  Deliberately loose - this catches the joint arm being wired
    //  up wrongly, not a subtle bias; the printed factors are the actual result.
    BOOST_CHECK_MESSAGE( joint_over_sample < 1.1,
                        "ratio " << ratio << ": joint limit " << pf << " is looser than the"
                        " sample-only limit " << pa << ", which more data cannot make it" );

    // D-21: how much wider the projected Currie limit's own scatter is than the future limit's.
    //  Expected about sqrt(k), which implies the predictive spread is understated by sqrt(1 + k).
    const double actual_spread = relative_iqr( currie_actual );
    const double projected_spread = relative_iqr( currie_scaled );
    const double spread_ratio = (actual_spread > 0.0) ? (projected_spread/actual_spread) : 0.0;

    cout << "DECON_PROJECTION_SPREAD," << ratio << ',' << currie_scaled.size() << ','
         << actual_spread << ',' << projected_spread << ',' << spread_ratio << ','
         << sqrt(ratio) << ',' << sqrt( 1.0 + spread_ratio*spread_ratio ) << '\n';

    // A prediction used for planning has to be in the right ballpark of what the dwell delivers;
    //  a factor of two either way would make it useless.  Deliberately loose - the point of the
    //  case is the printed ratios, and what they say about each method's bias.
    BOOST_CHECK_MESSAGE( (mb > 0.5*ma) && (mb < 2.0*ma),
                        "ratio " << ratio << ": background-reference prediction " << mb
                        << " against an actual median of " << ma );
    BOOST_CHECK_MESSAGE( (mc > 0.5*ma) && (mc < 2.0*ma),
                        "ratio " << ratio << ": scaled-spectrum profile projection " << mc
                        << " against an actual median of " << ma );
    BOOST_CHECK_MESSAGE( (me > 0.5*md) && (me < 2.0*md),
                        "ratio " << ratio << ": Currie projection " << me
                        << " against an actual Currie median of " << md );
  }//for( loop over dwell ratios )
}// BOOST_AUTO_TEST_CASE( RealSpectrumDwellProjectionAccuracy )
/** The coverage study runs far more trials than the real profile scan could, so it uses its own
 lightweight bisection.  That is only legitimate while it lands on the same answer production does,
 which is what this checks.  The harness computes its statistic with the production
 `poisson_deviance` / `fit_continuum_poisson` entry points, so the only thing that can differ is
 the search itself.
 */
BOOST_AUTO_TEST_CASE( DeconStudyProfileMatchesProduction )
{
  set_data_dir();

  const DetectionLimitCalc::DeconContinuumNorm norms[] = { DetectionLimitCalc::DeconContinuumNorm::Floating,
                                                           DetectionLimitCalc::DeconContinuumNorm::FixedByEdges,
                                                           DetectionLimitCalc::DeconContinuumNorm::FixedByFullRange };

  for( const DetectionLimitCalc::DeconContinuumNorm norm : norms )
  {
    for( const double continuum : { 0.1, 1.0, 100.0 } )
    {
      for( const double signal : { 0.0, 25.0 } )
      {
        const shared_ptr<const DetectionLimitCalc::DeconComputeInput> input =
          make_decon_input( signal, continuum, norm );
        const double maximum = (std::max)( 50.0, signal + 100.0 * sqrt( continuum + 1.0 ) );
        const StudyProfileResult study = study_cash_limit( *input, maximum );
        DetectionLimitCalc::DeconActivityOrDistanceLimitResult production;
        {
          ScopedCoutSilence silence;
          production = DetectionLimitCalc::get_activity_or_distance_limits( 0.95f, input, false, 0.0, maximum, true );
        }
        // Report the cell and carry on rather than aborting the whole case: a bare fatal check
        // here says only "something did not bracket", which is the least useful thing it could
        // say when the point of the test is to compare every cell.
        if( !study.found_upper || !production.foundUpperCl )
        {
          BOOST_ERROR( "no upper limit ("
                      << (study.found_upper ? "production" : "study") << ") for continuum="
                      << continuum << ", signal=" << signal << ", norm=" << static_cast<int>( norm )
                      << ", max=" << maximum << "; " << production.errorMessage );
          continue;
        }

        BOOST_CHECK_MESSAGE( fabs( study.upper - production.upperLimit ) < 0.03 * ( 1.0 + production.upperLimit ),
                             "Study profile " << study.upper << " differs from production " << production.upperLimit
                                              << " for continuum=" << continuum << ", signal=" << signal
                                              << ", norm=" << static_cast<int>( norm ) );
      }
    }
  }
}// BOOST_AUTO_TEST_CASE( DeconStudyProfileMatchesProduction )


namespace
{
/** One cell of the coverage study: a fully specified question, and how many trials answer it. */
struct StudyCell
{
  string block = "core";
  string mode = "Floating";
  string perturbation = "baseline";
  DetectionLimitCalc::DeconContinuumNorm norm = DetectionLimitCalc::DeconContinuumNorm::Floating;
  DetectionLimitCalc::DeconMeasurementModel model
        = DetectionLimitCalc::DeconMeasurementModel::CurrentSpectrum;
  DetectionLimitCalc::DeconLimitType limit_type
        = DetectionLimitCalc::DeconLimitType::OneSidedUpperLimit;
  double continuum = 10.0;             //!< counts per channel
  double signal_sigma = 1.0;           //!< truth signal, in units of sqrt(region background)
  double reference_over_sample = 0.0;  //!< BackgroundReference only; zero elsewhere
  size_t trials = 2000;

  /** Wall-clock a cell may use before it stops early, in seconds; zero means no cap.

   The per-trial cost varies by more than an order of magnitude across the grid - the continuum
   optimizer needs some 10,900 iterations per trial at one count per channel against 423 at a
   hundred - so a fixed trial count either makes the low-count cells unaffordable or makes the
   dense ones pointlessly precise.  A cell that stops early reports the trials it actually ran and
   emits a `DECON_COVERAGE_CUT` row; its Wilson interval widens accordingly, which is the honest
   consequence rather than a hidden one.
   */
  double max_seconds = 0.0;
};//struct StudyCell


/** The cell's identity as one string.  The CSV key and the seed both come from it, so the two
 cannot disagree about which cell a number belongs to.
 */
string study_cell_key( const StudyCell &cell )
{
  string key = cell.block + "|" + cell.mode + "|" + cell.perturbation
               + "|cont=" + study_number_token( cell.continuum )
               + "|sig=" + study_number_token( cell.signal_sigma );

  if( cell.limit_type == DetectionLimitCalc::DeconLimitType::CentralInterval )
    key += "|central";

  if( cell.model == DetectionLimitCalc::DeconMeasurementModel::BackgroundReference )
    key += "|ref=" + study_number_token( cell.reference_over_sample );

  return key;
}// study_cell_key(...)


/** The domain this study claims to have validated.

 Emitted as a `DECON_COVERAGE_DOMAIN` row so the report carries it rather than relying on someone
 remembering it.  Inside the domain the coverage assertions are `BOOST_CHECK`; outside they stay
 `BOOST_WARN` and the measured value is still recorded - the asymptotic chi-square threshold that
 defines the limit is not expected to hold at fractions of a count per channel, and a deliberately
 mis-specified model has no coverage guarantee to violate.
 */
bool in_supported_domain( const StudyCell &cell )
{
  // The always-on gate runs the same cells at 300 trials instead of 2000, where the study's
  //  thresholds are only one or two standard errors from the truth - a 0.075 false-detection
  //  ceiling against a true 0.05 is ~2 sigma at n=300, so an always-on test would go red a few
  //  percent of runs for no reason.  The gate asserts its own, deliberately looser bounds instead;
  //  excluding it here is what keeps those from being shadowed by these.
  if( cell.block == "gate" )
    return false;

  if( cell.limit_type != DetectionLimitCalc::DeconLimitType::OneSidedUpperLimit )
    return false;

  // Measured, not assumed: at 5*sqrt(B) coverage falls to about 0.93 even for `Floating` at one
  //  count per channel, so the claim stops at 3*sqrt(B).  Above that the cells are still run and
  //  reported, they are just not asserted against.
  if( cell.signal_sigma > 3.0 )
    return false;

  if( cell.model != DetectionLimitCalc::DeconMeasurementModel::CurrentSpectrum )
    return false;

  if( cell.perturbation != "baseline" )
    return false;

  if( cell.norm == DetectionLimitCalc::DeconContinuumNorm::Floating )
    return ( cell.continuum >= 1.0 );

  // `FixedByEdges` reads its background out of a handful of side channels; below roughly three
  //  counts per channel those channels no longer pin an offset and a slope well enough to trust.
  if( cell.norm == DetectionLimitCalc::DeconContinuumNorm::FixedByEdges )
    return ( cell.continuum >= 3.0 );

  return false;   // FixedByFullRange is retired: measured, never gated.
}// in_supported_domain(...)


/** Everything one cell reports. */
struct StudyResult
{
  StudyCell cell;
  uint32_t seed = 0;
  double truth = 0.0;
  double roi_background = 0.0;

  size_t trials = 0;
  size_t covered = 0;              //!< upper >= truth, over ALL trials (see coverage_strict)
  size_t scan_failures = 0;        //!< the profile never crossed the threshold
  size_t optimizer_failures = 0;   //!< at least one continuum fit reported non-convergence
  size_t found_lower = 0;          //!< claimed a detection; at truth zero this is the false-detection rate
  size_t interval_contains = 0;    //!< central cells only

  double estimator_bias_sum = 0.0;
  double upper_bias_sum = 0.0;
  double continuum_offset_sum = 0.0;
  double continuum_slope_sum = 0.0;
  size_t continuum_samples = 0;
  size_t iteration_sum = 0;

  vector<double> cpu_ms;
  vector<double> upper_limits;

  bool have_legacy = false;
  size_t legacy_covered = 0;
  size_t legacy_failures = 0;

  bool have_backref = false;
  size_t prediction_above = 0;
  size_t prediction_pairs = 0;
  double backref_ratio_median = 0.0;
  double backref_naive_ratio = 0.0;

  size_t production_checks = 0;
  size_t production_failures = 0;
  double production_max_relative_difference = 0.0;
};//struct StudyResult


double study_quantile( vector<double> values, const double fraction )
{
  if( values.empty() )
    return numeric_limits<double>::quiet_NaN();

  const size_t index = (std::min)( values.size() - 1,
                                   static_cast<size_t>( fraction * values.size() ) );
  nth_element( values.begin(), values.begin() + index, values.end() );
  return values[index];
}// study_quantile(...)


/** Writes one CSV row to `cout`, and to the file named by INTERSPEC_DECON_COVERAGE_CSV when set.

 Both, deliberately: a run of this length has to survive a crash near the end and be openable in a
 spreadsheet, while still leaving the numbers visible in the ordinary test log.
 */
void study_emit( const string &row )
{
  static ofstream *s_file = nullptr;
  static bool s_tried = false;

  cout << row << flush;

  if( !s_tried )
  {
    s_tried = true;
    const char * const path = getenv( "INTERSPEC_DECON_COVERAGE_CSV" );
    if( path && path[0] )
    {
      s_file = new ofstream( path, ios::out | ios::trunc );
      if( !s_file->is_open() )
      {
        BOOST_ERROR( "Could not open INTERSPEC_DECON_COVERAGE_CSV='" << path
                    << "'; continuing on stdout only" );
        delete s_file;
        s_file = nullptr;
      }
    }
  }

  if( s_file )
    (*s_file) << row << flush;
}// study_emit(...)


const char * const sm_study_csv_header =
  "DECON_COVERAGE,cell,seed,block,mode,measurement_model,limit_type,perturbation,"
  "continuum_counts_per_channel,signal_over_sqrt_roi_background,reference_over_sample_exposure,"
  "truth_counts,roi_background_counts,trials,in_supported_domain,statistic_kind,statistic,"
  "statistic_target,wilson_low,wilson_high,coverage_strict,coverage_conditional,estimator_bias,"
  "upper_bias,median_upper,interval_containment,detection_rate,scan_failure_rate,"
  "optimizer_failure_rate,last_continuum_offset,last_continuum_slope,mean_iterations,"
  "median_cpu_ms,p95_cpu_ms,backref_ratio_median,backref_naive_ratio,legacy_coverage,"
  "production_checks,production_failure_rate,production_max_relative_difference_vs_cash\n";


/** A search range generous enough that the profile scan is not the thing under test.

 Seeded from the Currie result exactly as the pre-restructure study did, so scan-failure rates stay
 comparable across the change.
 */
double study_seeded_maximum( const DetectionLimitCalc::CurrieMdaResult &currie )
{
  const double diff_multiple = 50.0;

  if( currie.source_counts > currie.decision_threshold )
  {
    const double lower_difference = fabs( currie.source_counts - currie.lower_limit );
    const double upper_difference = fabs( currie.upper_limit - currie.source_counts );
    return (std::max)( 1.0, currie.source_counts
                            + diff_multiple * (std::max)( lower_difference, upper_difference ) );
  }

  if( currie.upper_limit < 0.0f )
    return (std::max)( 1.0,
                      diff_multiple * sqrt( (std::max)( 0.0f, currie.peak_region_counts_sum ) ) );

  return (std::max)( 1.0, diff_multiple * currie.upper_limit );
}// study_seeded_maximum(...)


/** Builds the truth spectrum a cell describes.

 Perturbations that only re-describe the same data - a different region width, a different number
 of side channels, a mis-stated FWHM or branching ratio - are applied to the input rather than the
 counts, which is exactly what makes them model mis-specifications rather than different problems.
 */
shared_ptr<const DetectionLimitCalc::DeconComputeInput>
make_study_input( const StudyCell &cell, const double truth )
{
  shared_ptr<DetectionLimitCalc::DeconComputeInput> input
      = make_shared<DetectionLimitCalc::DeconComputeInput>(
          *make_decon_input( truth, cell.continuum, cell.norm ) );

  DetectionLimitCalc::DeconRoiInfo &roi = input->roi_info.front();
  DetectionLimitCalc::DeconRoiInfo::PeakInfo &peak = roi.peak_infos.front();
  const float centre = peak.energy;
  const float fwhm = peak.fwhm;

  if( cell.perturbation == "roi_narrow" )
  {
    roi.roi_start = centre - 0.75f*fwhm;
    roi.roi_end = centre + 0.75f*fwhm;
  } else if( cell.perturbation == "roi_wide" )
  {
    roi.roi_start = centre - 2.5f*fwhm;
    roi.roi_end = centre + 2.5f*fwhm;
  } else if( cell.perturbation == "side_2" )
  {
    roi.num_lower_side_channels = 2;
    roi.num_upper_side_channels = 2;
  } else if( cell.perturbation == "side_16" )
  {
    roi.num_lower_side_channels = 16;
    roi.num_upper_side_channels = 16;
  } else if( cell.perturbation == "fwhm_low" )
  {
    peak.fwhm = 0.8f*fwhm;      // the model is narrower than the data
  } else if( cell.perturbation == "fwhm_high" )
  {
    peak.fwhm = 1.25f*fwhm;     // and here wider
  } else if( cell.perturbation == "wrong_branching_ratio" )
  {
    peak.counts_per_bq_into_4pi = 0.5;
  }

  return input;
}// make_study_input(...)


/** Runs every trial of one `CurrentSpectrum` cell.

 \param audit_production How many trials additionally run the real
        `get_activity_or_distance_limits` scan, so the study's own search is checked against the
        one that ships.  With `StudyOptimizer::Production` this compares two searches over the
        same continuum fit; `DeconStudyProfileMatchesProduction` is what still compares two
        independent fits.
 \param want_legacy Also run the pre-Cash Neyman scan on the same trial spectra, for the report's
        before/after table.  Off everywhere but the `legacy` block, since it doubles the cost.
 */
void run_coverage_cell( StudyResult &out, const size_t audit_production, const bool want_legacy )
{
  const StudyCell &cell = out.cell;
  mt19937 generator( out.seed );

  const double delta = DetectionLimitCalc::decon_limit_delta( 0.95, cell.limit_type );
  const bool want_lower_root
      = ( cell.limit_type == DetectionLimitCalc::DeconLimitType::CentralInterval );

  const shared_ptr<const DetectionLimitCalc::DeconComputeInput> expected
      = make_study_input( cell, out.truth );

  const pair<size_t,size_t> roi_channels = DetectionLimitCalc::round_roi_to_channels(
      expected->measurement, expected->roi_info.front().roi_start,
      expected->roi_info.front().roi_end );
  out.roi_background = cell.continuum * ( 1 + roi_channels.second - roi_channels.first );

  const chrono::steady_clock::time_point cell_started = chrono::steady_clock::now();

  for( size_t trial_index = 0; trial_index < cell.trials; ++trial_index )
  {
    if( (cell.max_seconds > 0.0) && (trial_index > 0) )
    {
      const double elapsed = 1.0E-9 * chrono::duration_cast<chrono::nanoseconds>(
                               chrono::steady_clock::now() - cell_started ).count();
      if( elapsed > cell.max_seconds )
      {
        study_emit( "DECON_COVERAGE_CUT,trials," + study_cell_key( cell ) + ","
                    + std::to_string( trial_index ) + " of "
                    + std::to_string( cell.trials )
                    + " trials run before the per-cell time budget was reached\n" );
        break;
      }
    }

    const shared_ptr<const DetectionLimitCalc::DeconComputeInput> trial
        = make_poisson_trial( expected, generator );

    double seeded_maximum = 1.0;
    try
    {
      seeded_maximum = study_seeded_maximum( currie_for_decon_input( trial ) );
    }catch( std::exception & )
    {
      seeded_maximum = (std::max)( 50.0, out.truth + 100.0*sqrt( out.roi_background + 1.0 ) );
    }

    StudyDiagnostics diagnostics;
    const chrono::steady_clock::time_point started = chrono::steady_clock::now();
    const StudyProfileResult result = study_cash_limit( *trial, seeded_maximum, delta,
                                                       want_lower_root,
                                                       StudyOptimizer::Production, &diagnostics );
    const chrono::steady_clock::time_point finished = chrono::steady_clock::now();

    ++out.trials;
    out.cpu_ms.push_back( 1.0E-6 * chrono::duration_cast<chrono::nanoseconds>(
                                     finished - started ).count() );
    out.iteration_sum += diagnostics.continuum_iterations;
    out.optimizer_failures += ( diagnostics.continuum_failures > 0 ) ? 1 : 0;

    if( diagnostics.best_continuum.first != 0.0 )
    {
      out.continuum_offset_sum += diagnostics.best_continuum.first;
      out.continuum_slope_sum += diagnostics.best_continuum.second;
      ++out.continuum_samples;
    }

    // A scan that never crossed the threshold is NOT covered.  Dropping such trials from the
    //  denominator - which the study used to do - inflates coverage exactly in the cells where a
    //  mode fails to bracket, which are the cells most likely to be under-covering.
    if( !result.found_upper )
    {
      ++out.scan_failures;
    } else
    {
      out.covered += ( result.upper >= out.truth ) ? 1 : 0;
      out.estimator_bias_sum += result.best - out.truth;
      out.upper_bias_sum += result.upper - out.truth;
      out.upper_limits.push_back( result.upper );
      out.found_lower += result.found_lower ? 1 : 0;

      if( want_lower_root && result.found_lower )
        out.interval_contains += ( (result.lower <= out.truth) && (out.truth <= result.upper) )
                                   ? 1 : 0;
      else if( want_lower_root )
        out.interval_contains += ( out.truth <= result.upper ) ? 1 : 0;
    }

    if( out.production_checks < audit_production )
    {
      DetectionLimitCalc::DeconActivityOrDistanceLimitResult production;
      bool production_threw = false;
      try
      {
        ScopedCoutSilence silence;
        // With the cell's own limit type.  Defaulting to one-sided while the study profiled a
        //  central interval compared Delta=3.8415 against Delta=2.7055, which left every central
        //  cell unaudited and showed the mismatch as a scan disagreement.
        production = DetectionLimitCalc::get_activity_or_distance_limits(
                        0.95f, trial, false, 0.0, seeded_maximum, true, cell.limit_type );
      }catch( std::exception & )
      {
        production_threw = true;
      }

      ++out.production_checks;
      if( production_threw || !production.foundUpperCl || !result.found_upper )
      {
        ++out.production_failures;
      } else
      {
        const double relative_difference = fabs( production.upperLimit - result.upper )
                                           / ( 1.0 + fabs( production.upperLimit ) );
        out.production_max_relative_difference =
            (std::max)( out.production_max_relative_difference, relative_difference );
      }
    }//if( auditing production on this trial )

    if( want_legacy )
    {
      const StudyProfileResult legacy = study_neyman_limit( *trial, seeded_maximum );
      out.have_legacy = true;
      if( !legacy.found_upper || (legacy.upper > seeded_maximum) )
        ++out.legacy_failures;
      else
        out.legacy_covered += ( legacy.upper >= out.truth ) ? 1 : 0;
    }
  }//for( loop over trials )
}// run_coverage_cell(...)


/** What a `BackgroundReference` prediction looks like against an independently drawn future
 measurement.

 There is no "coverage" to measure here: the mode reports a *predicted median sensitivity* built
 from an expected-counts sample, not an interval over observed data.

 What this block measures is the coarser of the two available questions - whether the prediction
 still describes what the analyst gets if they take the future measurement and analyse it **on its
 own**, discarding the reference.  Each pair draws a reference `b`, takes production's prediction
 `S(b)` from it, draws an independent future sample `n` at truth zero, and takes an ordinary
 `Floating` limit on `n` alone.

 That is deliberately **not gated**, because the two are different estimands and the difference is
 not a defect: the joint model behind `S` also sees the reference block, so it is legitimately
 tighter, by an amount that grows as the reference shortens.  The measured ratio below (roughly
 0.8 at moderate count rates) is that effect, and reading it as an error would be a mistake.

 The sharper question - whether the expected-counts prediction sits at the *median* of what the same joint
 procedure returns when the sample is really observed - needs the sample block re-scored with `n`
 in place of its expectation, over the identical two-block channel list.  There is no production
 entry point that does that and this harness does not build one yet, so that calibration is
 **not measured**.  It is the first thing the next increment should add; until then nothing here
 licenses a claim that the expected-counts step is unbiased.
 */
void run_background_reference_cell( StudyResult &out )
{
  const StudyCell &cell = out.cell;
  mt19937 generator( out.seed );
  const double delta = DetectionLimitCalc::decon_limit_delta(
                          0.95, DetectionLimitCalc::DeconLimitType::OneSidedUpperLimit );

  const double sample_seconds = 1.0;
  const double reference_seconds = cell.reference_over_sample * sample_seconds;

  const shared_ptr<const DetectionLimitCalc::DeconComputeInput> reference_expectation
      = make_background_reference_input( cell.continuum, reference_seconds, sample_seconds );
  const shared_ptr<const DetectionLimitCalc::DeconComputeInput> sample_expectation
      = make_decon_input( 0.0, cell.continuum * sample_seconds,
                          DetectionLimitCalc::DeconContinuumNorm::Floating );

  vector<double> predicted, naive;
  out.have_backref = true;

  const chrono::steady_clock::time_point cell_started = chrono::steady_clock::now();

  for( size_t pair_index = 0; pair_index < cell.trials; ++pair_index )
  {
    // Same budget the coverage cells obey; without it these ran well past their slice and the
    //  run's total stopped matching what the budget arithmetic promised.
    if( (cell.max_seconds > 0.0) && (pair_index > 0) )
    {
      const double elapsed = 1.0E-9 * chrono::duration_cast<chrono::nanoseconds>(
                               chrono::steady_clock::now() - cell_started ).count();
      if( elapsed > cell.max_seconds )
      {
        study_emit( "DECON_COVERAGE_CUT,trials," + study_cell_key( cell ) + ","
                    + std::to_string( pair_index ) + " of " + std::to_string( cell.trials )
                    + " pairs run before the per-cell time budget was reached\n" );
        break;
      }
    }

    const shared_ptr<const DetectionLimitCalc::DeconComputeInput> reference
        = make_poisson_trial( reference_expectation, generator );

    const chrono::steady_clock::time_point started = chrono::steady_clock::now();

    double prediction = numeric_limits<double>::quiet_NaN();
    try
    {
      ScopedCoutSilence silence;
      const DetectionLimitCalc::DeconActivityOrDistanceLimitResult production =
          DetectionLimitCalc::get_activity_or_distance_limits(
              0.95f, reference, false, 0.0,
              (std::max)( 50.0, 200.0*sqrt( cell.continuum*sample_seconds + 1.0 ) ), true );
      if( production.foundUpperCl )
        prediction = production.upperLimit;
    }catch( std::exception & )
    {
    }

    // An independent future measurement of the same background, at truth zero.
    const shared_ptr<const DetectionLimitCalc::DeconComputeInput> sample
        = make_poisson_trial( sample_expectation, generator );

    const StudyProfileResult sample_limit = study_cash_limit(
        *sample, (std::max)( 50.0, 200.0*sqrt( cell.continuum*sample_seconds + 1.0 ) ),
        delta, false, StudyOptimizer::Production, nullptr );

    const chrono::steady_clock::time_point finished = chrono::steady_clock::now();
    out.cpu_ms.push_back( 1.0E-6 * chrono::duration_cast<chrono::nanoseconds>(
                                     finished - started ).count() );
    ++out.trials;

    if( IsNan(prediction) || !sample_limit.found_upper )
    {
      ++out.scan_failures;
      continue;
    }

    ++out.prediction_pairs;
    out.prediction_above += ( sample_limit.upper > prediction ) ? 1 : 0;
    predicted.push_back( prediction );
    naive.push_back( sample_limit.upper );
  }//for( loop over paired draws )

  const double median_predicted = study_quantile( predicted, 0.5 );
  const double median_naive = study_quantile( naive, 0.5 );

  // Left NaN on purpose: the joint re-scoring this would be the ratio of is not built (see above),
  //  and emitting the naive ratio under this name would quietly answer a different question.
  out.backref_ratio_median = numeric_limits<double>::quiet_NaN();
  out.backref_naive_ratio = ( median_naive > 0.0 ) ? (median_predicted/median_naive)
                                                   : numeric_limits<double>::quiet_NaN();
  out.upper_limits = naive;
}// run_background_reference_cell(...)


/** Formats one row, and makes the assertions the cell's position in the supported domain calls
 for.  One place, so a block cannot accidentally report a number it does not also gate.
 */
void study_report_cell( const StudyResult &out )
{
  const StudyCell &cell = out.cell;
  const bool supported = in_supported_domain( cell );
  const size_t bracketed = out.trials - out.scan_failures;

  const bool is_backref = out.have_backref;
  // The fraction of trials that claimed a detection.  At truth zero that is the false-detection
  //  rate; above it, the power.  One measurement with two readings, so the column is named for
  //  what it counts rather than for one of its interpretations.
  const double detection_rate = bracketed
        ? static_cast<double>( out.found_lower ) / static_cast<double>( bracketed )
        : numeric_limits<double>::quiet_NaN();

  const double coverage_strict = out.trials
        ? static_cast<double>( out.covered ) / static_cast<double>( out.trials )
        : numeric_limits<double>::quiet_NaN();
  const double coverage_conditional = bracketed
        ? static_cast<double>( out.covered ) / static_cast<double>( bracketed )
        : numeric_limits<double>::quiet_NaN();

  // At truth zero the quantity with a 0.05 target is the *detection* rate, not coverage - which is
  //  trivially 1 there, since any upper limit covers a truth of zero.  Reporting coverage against a
  //  0.05 target left three adjacent self-describing columns reading `1, 0.05, [0.998, 1]`, which
  //  looks like a catastrophic miss and is simply the wrong pairing.
  const bool zero_truth = !is_backref && !(out.truth > 0.0);

  const double statistic = is_backref
        ? ( out.prediction_pairs ? static_cast<double>( out.prediction_above )
                                 / static_cast<double>( out.prediction_pairs )
                             : numeric_limits<double>::quiet_NaN() )
        : ( zero_truth ? detection_rate : coverage_strict );
  const double statistic_target = is_backref ? 0.50 : ( zero_truth ? 0.05 : 0.95 );
  const pair<double,double> interval = is_backref
        ? wilson_interval( out.prediction_above, (std::max)( size_t(1), out.prediction_pairs ) )
        : ( zero_truth ? wilson_interval( out.found_lower, (std::max)( size_t(1), bracketed ) )
                       : wilson_interval( out.covered, (std::max)( size_t(1), out.trials ) ) );

  ostringstream row;
  row.precision( 6 );
  row << "DECON_COVERAGE," << study_cell_key( cell ) << ',' << out.seed << ','
      << cell.block << ',' << cell.mode << ','
      << ((cell.model == DetectionLimitCalc::DeconMeasurementModel::BackgroundReference)
            ? "BackgroundReference" : "CurrentSpectrum") << ','
      << ((cell.limit_type == DetectionLimitCalc::DeconLimitType::CentralInterval)
            ? "Central" : "OneSidedUpper") << ','
      << cell.perturbation << ',' << cell.continuum << ',' << cell.signal_sigma << ','
      << cell.reference_over_sample << ',' << out.truth << ',' << out.roi_background << ','
      << out.trials << ',' << (supported ? 1 : 0) << ','
      << (is_backref ? "PredictionVsSampleOnly" : "Coverage") << ',' << statistic << ','
      << statistic_target << ',' << interval.first << ',' << interval.second << ',';

  if( is_backref )
    row << ",,";                                   // coverage is undefined for a prediction
  else
    row << coverage_strict << ',' << coverage_conditional << ',';

  row << ( bracketed ? out.estimator_bias_sum/bracketed : numeric_limits<double>::quiet_NaN() )
      << ',' << ( bracketed ? out.upper_bias_sum/bracketed : numeric_limits<double>::quiet_NaN() )
      << ',' << study_quantile( out.upper_limits, 0.5 ) << ',';

  if( cell.limit_type == DetectionLimitCalc::DeconLimitType::CentralInterval )
    row << ( out.trials ? static_cast<double>(out.interval_contains)/out.trials : 0.0 );
  row << ',';

  if( is_backref )
    row << ',';                                    // no false-detection rate for a prediction
  else
    row << detection_rate << ',';

  row << ( out.trials ? static_cast<double>(out.scan_failures)/out.trials : 0.0 ) << ','
      << ( out.trials ? static_cast<double>(out.optimizer_failures)/out.trials : 0.0 ) << ','
      << ( out.continuum_samples ? out.continuum_offset_sum/out.continuum_samples : 0.0 ) << ','
      << ( out.continuum_samples ? out.continuum_slope_sum/out.continuum_samples : 0.0 ) << ','
      << ( out.trials ? static_cast<double>(out.iteration_sum)/out.trials : 0.0 ) << ','
      << study_quantile( out.cpu_ms, 0.5 ) << ',' << study_quantile( out.cpu_ms, 0.95 ) << ',';

  if( is_backref )
    row << out.backref_ratio_median << ',' << out.backref_naive_ratio << ',';
  else
    row << ",,";

  if( out.have_legacy )
  {
    const size_t legacy_valid = out.trials - out.legacy_failures;
    row << ( legacy_valid ? static_cast<double>(out.legacy_covered)/legacy_valid : 0.0 );
  }
  row << ',';

  row << out.production_checks << ','
      << ( out.production_checks
             ? static_cast<double>(out.production_failures)/out.production_checks : 0.0 )
      << ',' << out.production_max_relative_difference << '\n';

  study_emit( row.str() );

  const string what = study_cell_key( cell );

  // The production audit has to hold everywhere; it is a statement about the search, not about
  //  whether the statistic covers.
  BOOST_WARN_MESSAGE( out.production_failures == 0,
                     what << ": production scan audit failed " << out.production_failures
                     << " of " << out.production_checks << " times" );
  BOOST_WARN_MESSAGE( out.production_max_relative_difference < 0.03,
                     what << ": study and production scans differ by "
                     << out.production_max_relative_difference );

  if( is_backref )
  {
    // Reported, never gated.  This compares two different estimands - the joint prediction against
    //  a sample-only limit - so a departure from 0.5 is the information the joint model carries,
    //  not an error.  Gating it would assert that the reference contributes nothing.
    BOOST_TEST_MESSAGE( what << ": prediction exceeds a sample-only limit in "
                       << (1.0 - statistic) << " of pairs; predicted/sample-only median ratio "
                       << out.backref_naive_ratio );

    // What *would* be a defect is failing to produce a prediction at all.
    BOOST_CHECK_MESSAGE( out.prediction_pairs > (out.trials/2),
                        what << ": only " << out.prediction_pairs << " of " << out.trials
                        << " pairs produced both a prediction and a sample limit" );
    return;
  }//if( is_backref )

  // The sharp statement - "the Wilson interval contains the nominal rate" - stays a WARN forever.
  //  It is false by construction for 5% of cells, so across the grid at least one would go red
  //  most runs, which is how a test gets ignored.  What is CHECKed is an absolute floor that
  //  Monte Carlo noise cannot reach: at 2000 trials the Wilson low for a true 0.93 is 0.918.
  if( out.truth > 0.0 )
  {
    BOOST_WARN_MESSAGE( (interval.first <= 0.95) && (interval.second >= 0.95),
                       what << ": coverage " << coverage_strict << " ["
                       << interval.first << ", " << interval.second << "] excludes 95%" );

    if( supported )
    {
      BOOST_CHECK_MESSAGE( coverage_strict >= 0.93,
                          what << ": coverage " << coverage_strict << " is below the 0.93 floor" );

      // Only the floor is asserted.  Over-coverage is *expected* for a parameter bounded at zero:
      //  when the truth sits near the boundary the upper limit almost never falls below it, so a
      //  one-sided bound at 0.5*sqrt(B) legitimately covers ~100% of the time.  Asserting an upper
      //  bound on coverage would be asserting that the method is not conservative near zero, which
      //  is not a property anyone wants.  It is reported, and warned on, but not gated.
      BOOST_WARN_MESSAGE( coverage_strict <= 0.99,
                         what << ": coverage " << coverage_strict << " exceeds 0.99 - conservative "
                         "here, which is expected close to the boundary at zero" );
    } else
    {
      BOOST_WARN_MESSAGE( coverage_strict >= 0.93,
                         what << ": coverage " << coverage_strict
                         << " below 0.93 (outside the supported domain)" );
    }
  } else
  {
    const pair<double,double> false_interval
        = wilson_interval( out.found_lower, (std::max)( size_t(1), bracketed ) );
    BOOST_WARN_MESSAGE( (false_interval.first <= 0.05) && (false_interval.second >= 0.05),
                       what << ": false-detection rate " << detection_rate << " ["
                       << false_interval.first << ", " << false_interval.second
                       << "] excludes 5%" );

    if( supported )
      BOOST_CHECK_MESSAGE( detection_rate <= 0.075,
                          what << ": false-detection rate " << detection_rate
                          << " exceeds the 0.075 ceiling" );
    else
      BOOST_WARN_MESSAGE( detection_rate <= 0.075,
                         what << ": false-detection rate " << detection_rate
                         << " exceeds 0.075 (outside the supported domain)" );
  }

  if( supported )
  {
    BOOST_CHECK_MESSAGE( out.scan_failures <= (out.trials/200 + 1),
                        what << ": " << out.scan_failures << " of " << out.trials
                        << " scans failed to bracket" );
    BOOST_CHECK_MESSAGE( out.optimizer_failures == 0,
                        what << ": " << out.optimizer_failures
                        << " trials had a non-converged continuum fit" );
  } else
  {
    BOOST_WARN_MESSAGE( out.optimizer_failures == 0,
                       what << ": " << out.optimizer_failures
                       << " trials had a non-converged continuum fit" );
  }
}// study_report_cell(...)


/** Runs a cell end to end: seed, trials, report, assertions. */
void study_run_and_report( StudyCell cell, const size_t audit_production = 10,
                           const bool want_legacy = false )
{
  StudyResult out;
  out.cell = cell;
  out.seed = study_cell_seed( study_cell_key( cell ) );

  if( cell.model == DetectionLimitCalc::DeconMeasurementModel::BackgroundReference )
  {
    out.truth = 0.0;
    run_background_reference_cell( out );
  } else
  {
    // The truth signal is quoted against the region background, so it has to be computed from the
    //  same channel rounding production uses.
    const shared_ptr<const DetectionLimitCalc::DeconComputeInput> background
        = make_study_input( cell, 0.0 );
    const pair<size_t,size_t> roi_channels = DetectionLimitCalc::round_roi_to_channels(
        background->measurement, background->roi_info.front().roi_start,
        background->roi_info.front().roi_end );
    const double roi_background = cell.continuum * ( 1 + roi_channels.second - roi_channels.first );
    out.truth = cell.signal_sigma * sqrt( roi_background );
    run_coverage_cell( out, audit_production, want_legacy );
  }

  study_report_cell( out );
}// study_run_and_report(...)
}//namespace


BOOST_AUTO_TEST_CASE( DeconCoverageGate )
{
  set_data_dir();

  // `DeconCoverageStudy` is a report: opt-in, tens of minutes, and nobody runs it on a pull
  //  request.  A statistic that nothing checks is one that silently rots, so a handful of its
  //  cells run on every build - with thresholds loose enough that Monte Carlo noise cannot trip
  //  them, and tight enough that the kind of regression that took `FixedByFullRange` to 40%
  //  cannot pass.
  study_emit( sm_study_csv_header );

  // Deliberately 10 and 100 counts per channel, not 1: a trial at one count per channel costs
  //  about eighteen times as much (the continuum optimizer needs ~10,900 iterations against ~423),
  //  which would put this test into minutes on every build.  The low-count regime is inside the
  //  supported domain and is covered by the opt-in study, not here - said out loud rather than
  //  left to be inferred from the cell list.
  study_emit( "DECON_COVERAGE_CUT,gate,continuum=0.1|1,"
              "low-count cells cost ~18x per trial; covered by the opt-in study instead\n" );

  const double continua[] = { 10.0, 100.0 };
  const double signals[] = { 0.0, 1.0 };

  for( const double continuum : continua )
  {
    for( const double signal : signals )
    {
      StudyCell cell;
      cell.block = "gate";
      cell.mode = "Floating";
      cell.norm = DetectionLimitCalc::DeconContinuumNorm::Floating;
      cell.continuum = continuum;
      cell.signal_sigma = signal;
      cell.trials = 300;
      cell.max_seconds = 25.0;

      StudyResult out;
      out.cell = cell;
      out.seed = study_cell_seed( study_cell_key( cell ) );

      const shared_ptr<const DetectionLimitCalc::DeconComputeInput> background
          = make_study_input( cell, 0.0 );
      const pair<size_t,size_t> roi_channels = DetectionLimitCalc::round_roi_to_channels(
          background->measurement, background->roi_info.front().roi_start,
          background->roi_info.front().roi_end );
      const double roi_background = continuum * ( 1 + roi_channels.second - roi_channels.first );
      out.truth = signal * sqrt( roi_background );

      run_coverage_cell( out, 5, false );

      // Deliberately looser than the study's own gates: at 300 trials the Wilson low for a true
      //  0.95 is about 0.918, so a 0.88 floor is roughly three half-widths clear of noise while
      //  still catching a drop to 0.90.
      const double coverage = static_cast<double>( out.covered ) / out.trials;
      const size_t bracketed = out.trials - out.scan_failures;
      const double false_rate = bracketed
          ? static_cast<double>( out.found_lower ) / bracketed : 0.0;

      study_report_cell( out );

      // A cell that ran out of its time slice on a loaded machine must not then assert on a
      //  handful of trials, where these bounds are coin flips.
      BOOST_REQUIRE_MESSAGE( out.trials >= 200,
                            "gate continuum=" << continuum << ": only " << out.trials
                            << " trials completed; the thresholds below are not meaningful" );

      if( out.truth > 0.0 )
        BOOST_CHECK_MESSAGE( coverage >= 0.88,
                            "gate continuum=" << continuum << ": coverage " << coverage );
      else
        BOOST_CHECK_MESSAGE( false_rate <= 0.12,
                            "gate continuum=" << continuum << ": false-detection rate "
                            << false_rate );

      BOOST_CHECK_MESSAGE( out.scan_failures == 0,
                          "gate continuum=" << continuum << ": " << out.scan_failures
                          << " scans failed to bracket" );
      BOOST_CHECK_MESSAGE( out.optimizer_failures == 0,
                          "gate continuum=" << continuum << ": " << out.optimizer_failures
                          << " non-converged continuum fits" );
      BOOST_CHECK_MESSAGE( out.production_max_relative_difference < 0.03,
                          "gate continuum=" << continuum << ": study and production differ by "
                          << out.production_max_relative_difference );
    }//for( loop over signal levels )
  }//for( loop over continuum levels )
}// BOOST_AUTO_TEST_CASE( DeconCoverageGate )


#if ( PERFORM_DEVELOPER_CHECKS )
BOOST_AUTO_TEST_CASE( DeconCoverageStudy )
{
  const char * const enabled = getenv( "INTERSPEC_RUN_DECON_COVERAGE" );
  if( !enabled || ( string( enabled ) != "1" ) )
  {
    BOOST_TEST_MESSAGE( "Set INTERSPEC_RUN_DECON_COVERAGE=1 to run the deconvolution coverage study." );
    return;
  }

  set_data_dir();

  // The whole run is held to a wall-clock budget, divided evenly over the cells.  Without it the
  //  low-count cells alone run for many hours, and a study nobody can afford to run is a study
  //  nobody runs.
  double budget_hours = 1.0;
  if( const char * const budget = getenv( "INTERSPEC_DECON_COVERAGE_BUDGET_HOURS" ) )
    budget_hours = (std::max)( 0.05, atof( budget ) );

  const size_t planned_cells = 48 + 12 + 16 + 14 + 12;
  const double seconds_per_cell = 3600.0*budget_hours/static_cast<double>( planned_cells );
  study_emit( "DECON_COVERAGE_CALIBRATION,budget_hours," + study_number_token( budget_hours )
              + ",cells," + std::to_string( planned_cells ) + ",seconds_per_cell,"
              + study_number_token( seconds_per_cell ) + "\n" );

  // Must state exactly what `in_supported_domain()` gates, and nothing broader: this row is the
  //  study's central claim and travels with the data.
  study_emit( "DECON_COVERAGE_DOMAIN,One-sided upper limit; CurrentSpectrum; a single unskewed "
              "line whose FWHM is correctly modelled; linear continuum; the baseline region of "
              "+/-1.375 FWHM with 4 side channels each side; truth 0 to 3*sqrt(B); "
              "Floating at >=1 count/channel; FixedByEdges at >=3 counts/channel\n" );

  study_emit( sm_study_csv_header );

  const double continua[] = { 0.1, 1.0, 10.0, 100.0 };
  const double signals[] = { 0.0, 0.5, 1.0, 2.0, 3.0, 5.0 };

  // ---- core: the two selectable continuum treatments, on the spectrum in hand ----------------
  const pair<const char *, DetectionLimitCalc::DeconContinuumNorm> core_modes[] = {
    { "Floating", DetectionLimitCalc::DeconContinuumNorm::Floating },
    { "FixedByEdges", DetectionLimitCalc::DeconContinuumNorm::FixedByEdges }
  };

  for( const auto &mode : core_modes )
  {
    for( const double continuum : continua )
    {
      for( const double signal : signals )
      {
        StudyCell cell;
        cell.block = "core";
        cell.mode = mode.first;
        cell.norm = mode.second;
        cell.continuum = continuum;
        cell.signal_sigma = signal;
        cell.trials = 2000;
        cell.max_seconds = seconds_per_cell;
        study_run_and_report( cell );
      }
    }
  }

  // ---- central: the two-sided interval, on the validated Floating subset ----------------------
  for( const double continuum : continua )
  {
    for( const double signal : { 2.0, 3.0, 5.0 } )
    {
      StudyCell cell;
      cell.block = "central";
      cell.mode = "Floating";
      cell.norm = DetectionLimitCalc::DeconContinuumNorm::Floating;
      cell.limit_type = DetectionLimitCalc::DeconLimitType::CentralInterval;
      cell.continuum = continuum;
      cell.signal_sigma = signal;
      cell.trials = 2000;
        cell.max_seconds = seconds_per_cell;
      study_run_and_report( cell );
    }
  }

  // ---- backref: is the expected-counts prediction really the median? ------------------------------------
  for( const double continuum : continua )
  {
    for( const double ratio : { 0.25, 1.0, 4.0, 20.0 } )
    {
      StudyCell cell;
      cell.block = "backref";
      cell.mode = "BackgroundReference";
      cell.norm = DetectionLimitCalc::DeconContinuumNorm::Floating;
      cell.model = DetectionLimitCalc::DeconMeasurementModel::BackgroundReference;
      cell.continuum = continuum;
      cell.signal_sigma = 0.0;
      cell.reference_over_sample = ratio;
      cell.trials = 1000;
        cell.max_seconds = seconds_per_cell;
      study_run_and_report( cell, 0 );
    }
  }

  // ---- perturbations: one factor at a time off the documented baseline --------------------------
  const char * const perturbations[] = { "roi_narrow", "roi_wide", "side_2", "side_16",
                                         "fwhm_low", "fwhm_high" };
  for( const char * const perturbation : perturbations )
  {
    for( const double signal : { 0.0, 1.0 } )
    {
      StudyCell cell;
      cell.block = "perturbation";
      cell.mode = "Floating";
      cell.perturbation = perturbation;
      cell.norm = DetectionLimitCalc::DeconContinuumNorm::Floating;
      cell.continuum = 10.0;
      cell.signal_sigma = signal;
      cell.trials = 1000;
        cell.max_seconds = seconds_per_cell;
      study_run_and_report( cell );
    }
  }

  // Not silently absent.  These need spectra this fixture cannot synthesize yet - a sloped or
  //  curved continuum, a second line inside or beside the region, and overlapping regions that
  //  production must merge.  They are the next thing to add, and until then the supported domain
  //  above is not a claim about them.
  const char * const not_covered[] = { "sloped", "curved", "second_line_strong",
                                       "second_line_weak", "separated_rois", "overlapping_rois",
                                       "interference", "skew_mismatch" };
  for( const char * const perturbation : not_covered )
    study_emit( string("DECON_COVERAGE_CUT,perturbation,") + perturbation
                + ",needs a spectrum builder this fixture does not have yet\n" );

  // Not merely unimplemented - it cannot be measured by this harness at all.  `make_study_peak_model`
  //  builds `peak_fraction` from the line shape alone and never reads `counts_per_bq_into_4pi`, so
  //  the study's "activity" is counts where production's is becquerel.  Mis-stating the branching
  //  ratio therefore cancels inside the harness and the cell silently reproduces the baseline - it
  //  did, at 0.976 against the baseline's 0.978 - while the production audit on the same cells
  //  diverges by the full factor.  Measuring it needs the harness to carry efficiency.
  study_emit( "DECON_COVERAGE_CUT,perturbation,wrong_branching_ratio,"
              "the study harness works in counts and ignores counts_per_bq_into_4pi, so this "
              "cancels internally and would report the baseline\n" );

  // ---- legacy: the pre-Cash Neyman scan, for the report's before/after table --------------------
  const pair<const char *, DetectionLimitCalc::DeconContinuumNorm> legacy_modes[] = {
    { "Floating", DetectionLimitCalc::DeconContinuumNorm::Floating },
    { "FixedByEdges", DetectionLimitCalc::DeconContinuumNorm::FixedByEdges },
    { "FixedByFullRange", DetectionLimitCalc::DeconContinuumNorm::FixedByFullRange }
  };

  for( const auto &mode : legacy_modes )
  {
    for( const double continuum : { 1.0, 100.0 } )
    {
      for( const double signal : { 0.0, 1.0 } )
      {
        StudyCell cell;
        cell.block = "legacy";
        cell.mode = mode.first;
        cell.norm = mode.second;
        cell.continuum = continuum;
        cell.signal_sigma = signal;
        cell.trials = 1000;
        cell.max_seconds = seconds_per_cell;
        study_run_and_report( cell, 10, true );
      }
    }
  }
}// BOOST_AUTO_TEST_CASE( DeconCoverageStudy )
#endif


BOOST_AUTO_TEST_CASE( DeconCharacteristicLimitsFisherRootAndSignedEstimate )
{
  using namespace DetectionLimitCalc;

  DeconCharacteristicLimitInput characteristic;
  characteristic.decon_input = *make_decon_input( 0.0, 100.0 );
  characteristic.alpha = 0.05;
  characteristic.beta = 0.05;

  const DeconCharacteristicLimitResult result = decon_characteristic_limits( characteristic );
  BOOST_REQUIRE_MESSAGE( result.status == DeconCharacteristicLimitStatus::Success,
                         result.error_message );
  BOOST_CHECK_SMALL( result.primary_activity, 2.0E-3 );
  BOOST_CHECK_GT( result.decision_threshold, 0.0 );
  BOOST_CHECK_GT( result.detection_limit, result.decision_threshold );

  // Hand calculation of the 3x3 expected Fisher matrix for activity plus a linear continuum.
  // The ROI is 1162--1184 keV in 2 keV channels, centred on the 1173 keV line.  Counts are exactly
  // 100 per channel at A=0.  Schur-complementing the two continuum coefficients must give the same
  // complete activity covariance as the artificial-spectrum calculation.
  const double peak_energy = 1173.0;
  const double sigma = 8.0 / PhysicalUnits::fwhm_nsigma;
  const double reference_energy = peak_energy;
  double information_aa = 0.0;
  double information_ac0 = 0.0, information_ac1 = 0.0;
  double information_c00 = 0.0, information_c01 = 0.0, information_c11 = 0.0;

  for( size_t channel = 581; channel <= 591; ++channel )
  {
    const double lower = 2.0*static_cast<double>(channel);
    const double upper = lower + 2.0;
    const double signal = 0.5*( std::erf((upper - peak_energy)/(std::sqrt(2.0)*sigma))
                                - std::erf((lower - peak_energy)/(std::sqrt(2.0)*sigma)) );
    const double basis0 = upper - lower;
    const double x0 = lower - reference_energy;
    const double x1 = upper - reference_energy;
    const double basis1 = 0.5*(x1*x1 - x0*x0);
    const double inverse_expected = 1.0 / 100.0;

    information_aa += signal*signal*inverse_expected;
    information_ac0 += signal*basis0*inverse_expected;
    information_ac1 += signal*basis1*inverse_expected;
    information_c00 += basis0*basis0*inverse_expected;
    information_c01 += basis0*basis1*inverse_expected;
    information_c11 += basis1*basis1*inverse_expected;
  }

  const double continuum_determinant
          = information_c00*information_c11 - information_c01*information_c01;
  const double nuisance_projection
          = ( information_ac0*information_ac0*information_c11
              - 2.0*information_ac0*information_ac1*information_c01
              + information_ac1*information_ac1*information_c00 ) / continuum_determinant;
  const double hand_uncertainty = 1.0 / std::sqrt( information_aa - nuisance_projection );
  BOOST_CHECK_CLOSE_FRACTION( result.uncertainty_at_zero, hand_uncertainty, 2.0E-5 );

  const double k = 1.6448536269514722;
  BOOST_CHECK_CLOSE_FRACTION( result.decision_threshold,
                              k*result.uncertainty_at_zero, 2.0E-8 );
  BOOST_CHECK_CLOSE_FRACTION( result.detection_limit,
                              result.decision_threshold
                                + k*result.uncertainty_at_detection_limit, 2.0E-8 );

  characteristic.decon_input = *make_decon_input( -12.0, 100.0 );
  const DeconCharacteristicLimitResult negative = decon_characteristic_limits( characteristic );
  BOOST_REQUIRE_MESSAGE( negative.status == DeconCharacteristicLimitStatus::Success,
                         negative.error_message );
  BOOST_CHECK_LT( negative.primary_activity, 0.0 );
  BOOST_CHECK_CLOSE_FRACTION( negative.primary_activity, -12.0, 2.0E-3 );
}//BOOST_AUTO_TEST_CASE( DeconCharacteristicLimitsFisherRootAndSignedEstimate )


BOOST_AUTO_TEST_CASE( DeconCharacteristicLimitsIndependentProbabilitiesAndPositivity )
{
  using namespace DetectionLimitCalc;

  DeconCharacteristicLimitInput nominal;
  nominal.decon_input = *make_decon_input( 0.0, 10.0 );
  nominal.alpha = 0.05;
  nominal.beta = 0.05;
  const DeconCharacteristicLimitResult base = decon_characteristic_limits( nominal );
  BOOST_REQUIRE_MESSAGE( base.status == DeconCharacteristicLimitStatus::Success,
                         base.error_message );

  DeconCharacteristicLimitInput lower_alpha = nominal;
  lower_alpha.alpha = 0.01;
  const DeconCharacteristicLimitResult alpha_result
                                            = decon_characteristic_limits( lower_alpha );
  BOOST_REQUIRE_MESSAGE( alpha_result.status == DeconCharacteristicLimitStatus::Success,
                         alpha_result.error_message );
  BOOST_CHECK_GT( alpha_result.decision_threshold, base.decision_threshold );
  BOOST_CHECK_GT( alpha_result.detection_limit, base.detection_limit );

  DeconCharacteristicLimitInput lower_beta = nominal;
  lower_beta.beta = 0.01;
  const DeconCharacteristicLimitResult beta_result = decon_characteristic_limits( lower_beta );
  BOOST_REQUIRE_MESSAGE( beta_result.status == DeconCharacteristicLimitStatus::Success,
                         beta_result.error_message );
  BOOST_CHECK_CLOSE_FRACTION( beta_result.decision_threshold, base.decision_threshold, 1.0E-12 );
  BOOST_CHECK_GT( beta_result.detection_limit, base.detection_limit );

  // A signed fixed component is legal.  The fitted constant continuum must offset it enough that
  // every *total* expectation remains positive; the activity component itself is not clipped.
  std::vector<PoissonChannel> channels( 3 );
  for( size_t i = 0; i < channels.size(); ++i )
  {
    channels[i].lower_energy = static_cast<double>(i);
    channels[i].upper_energy = static_cast<double>(i + 1);
    channels[i].observed = 1.0;
    channels[i].fixed_signal = (i == 1) ? -20.0 : -2.0;
  }
  const PoissonContinuumFit fit
                    = fit_continuum_poisson( channels.data(), channels.size(), 1, 0.0, {} );
  BOOST_REQUIRE_MESSAGE( fit.converged, fit.error );
  for( const PoissonChannel &channel : channels )
    BOOST_CHECK_GT( fit.coefficients[0]*(channel.upper_energy - channel.lower_energy)
                      + channel.fixed_signal, 0.0 );

  DeconComputeInput signed_trial = *make_decon_input( 0.0, 10.0 );
  signed_trial.activity = -1.0;
  BOOST_CHECK_THROW( decon_compute_peaks(signed_trial), std::runtime_error );
  signed_trial.allow_negative_activity = true;
  BOOST_CHECK_NO_THROW( decon_compute_peaks(signed_trial) );
}//BOOST_AUTO_TEST_CASE( DeconCharacteristicLimitsIndependentProbabilitiesAndPositivity )


BOOST_AUTO_TEST_CASE( DeconCharacteristicLimitsModelsAndTypedFailures )
{
  using namespace DetectionLimitCalc;

  // Two channels and two continuum coefficients leave at most rank two for three joint
  // parameters (activity, offset, slope), so this model has a known singular information matrix.
  DeconCharacteristicLimitInput singular;
  singular.decon_input = *make_decon_input( 0.0, 10.0 );
  singular.decon_input.roi_info[0].roi_start = 1172.0f;
  singular.decon_input.roi_info[0].roi_end = 1176.0f;
  const DeconCharacteristicLimitResult singular_result
                                              = decon_characteristic_limits( singular );
  BOOST_CHECK_EQUAL( static_cast<int>(singular_result.status),
                     static_cast<int>(DeconCharacteristicLimitStatus::SingularInformation) );
  BOOST_CHECK( !std::isfinite(singular_result.detection_limit) );

  DeconCharacteristicLimitInput unsupported;
  unsupported.decon_input = *make_decon_input( 0.0, 10.0 );
  unsupported.decon_input.measurement_model = DeconMeasurementModel::BackgroundReference;
  const DeconCharacteristicLimitResult unsupported_result
                                            = decon_characteristic_limits( unsupported );
  BOOST_CHECK_EQUAL( static_cast<int>(unsupported_result.status),
                     static_cast<int>(DeconCharacteristicLimitStatus::UnsupportedModel) );

  DeconCharacteristicLimitInput invalid = unsupported;
  invalid.decon_input.measurement_model = DeconMeasurementModel::CurrentSpectrum;
  invalid.decon_input.roi_info[0].peak_infos.clear();
  const DeconCharacteristicLimitResult invalid_result
                                            = decon_characteristic_limits( invalid );
  BOOST_CHECK_EQUAL( static_cast<int>(invalid_result.status),
                     static_cast<int>(DeconCharacteristicLimitStatus::InvalidInput) );

  // Exercise overlapping ROIs, multiple lines, nonuniform channels, and a quadratic FixedByEdges
  // continuum in one identifiable model.  The sidebands remain in the joint likelihood after the
  // overlap merge.
  DeconCharacteristicLimitInput combined;
  combined.decon_input
        = *make_decon_input( 20.0, 100.0, DeconContinuumNorm::FixedByEdges );
  combined.decon_input.roi_info[0].continuum_type = PeakContinuum::OffsetType::Quadratic;

  DeconRoiInfo second = combined.decon_input.roi_info[0];
  second.roi_start = 1170.0f;
  second.roi_end = 1190.0f;
  second.peak_infos.clear();
  DeconRoiInfo::PeakInfo second_line;
  second_line.energy = 1179.0f;
  second_line.fwhm = 8.0f;
  second_line.counts_per_bq_into_4pi = 0.4;
  second.peak_infos.push_back( second_line );
  combined.decon_input.roi_info.push_back( second );

  const size_t num_channels = combined.decon_input.measurement->num_gamma_channels();
  std::vector<float> nonuniform_edges( num_channels + 1, 0.0f );
  for( size_t channel = 0; channel <= num_channels; ++channel )
  {
    const double x = static_cast<double>(channel);
    nonuniform_edges[channel] = static_cast<float>( 2.0*x + 1.0E-6*x*x );
  }
  std::shared_ptr<SpecUtils::EnergyCalibration> calibration
                                          = std::make_shared<SpecUtils::EnergyCalibration>();
  calibration->set_lower_channel_energy( num_channels, nonuniform_edges );
  std::shared_ptr<SpecUtils::Measurement> measurement
                                          = std::make_shared<SpecUtils::Measurement>();
  measurement->set_gamma_counts(
                std::make_shared<std::vector<float>>(*combined.decon_input.measurement->gamma_counts()),
                1.0f, 1.0f );
  measurement->set_energy_calibration( calibration );
  combined.decon_input.measurement = measurement;

  const DeconCharacteristicLimitResult combined_result
                                              = decon_characteristic_limits( combined );
  BOOST_REQUIRE_MESSAGE( combined_result.status == DeconCharacteristicLimitStatus::Success,
                         combined_result.error_message );
  BOOST_CHECK( std::isfinite(combined_result.detection_limit) );
  BOOST_CHECK( !combined_result.warnings.empty() );

  // Two separated regions retain independent nuisance coefficients but share the one activity.
  // The second line must add information rather than being treated as a separate activity fit.
  DeconCharacteristicLimitInput separated;
  separated.decon_input = *make_decon_input( 0.0, 100.0, DeconContinuumNorm::Floating );
  const DeconCharacteristicLimitResult one_region
                                            = decon_characteristic_limits( separated );
  BOOST_REQUIRE_MESSAGE( one_region.status == DeconCharacteristicLimitStatus::Success,
                         one_region.error_message );
  DeconRoiInfo separated_roi = separated.decon_input.roi_info[0];
  separated_roi.roi_start = 1388.0f;
  separated_roi.roi_end = 1412.0f;
  separated_roi.peak_infos[0].energy = 1400.0f;
  separated.decon_input.roi_info.push_back( separated_roi );
  const DeconCharacteristicLimitResult separated_result
                                            = decon_characteristic_limits( separated );
  BOOST_REQUIRE_MESSAGE( separated_result.status == DeconCharacteristicLimitStatus::Success,
                         separated_result.error_message );
  BOOST_CHECK_LT( separated_result.uncertainty_at_zero, one_region.uncertainty_at_zero );

  // FixedByEdges at 0.1 counts/channel still produces a limit, but it is the cell the Monte Carlo
  // excludes as outside the checked range - it measured a false-positive probability of 0.083.  A
  // bare `Success` would present an unchecked answer as a checked one, so the caller must be
  // warned.  A 1-count control region is inside the checked range and must not be warned about,
  // or the warning would be noise everywhere.
  DeconCharacteristicLimitInput sparse_edges;
  sparse_edges.decon_input
        = *make_decon_input( 0.0, 0.1, DeconContinuumNorm::FixedByEdges );
  const DeconCharacteristicLimitResult sparse_result
                                            = decon_characteristic_limits( sparse_edges );
  BOOST_REQUIRE_MESSAGE( sparse_result.status == DeconCharacteristicLimitStatus::Success,
                         sparse_result.error_message );
  const auto has_control_warning = []( const DeconCharacteristicLimitResult &result ) -> bool {
    return std::any_of( begin(result.warnings), end(result.warnings),
                       []( const std::string &warning ){
                         return warning.find("anchor its continuum") != std::string::npos;
                       } );
  };
  BOOST_CHECK( has_control_warning(sparse_result) );

  DeconCharacteristicLimitInput checked_edges;
  checked_edges.decon_input
        = *make_decon_input( 0.0, 1.0, DeconContinuumNorm::FixedByEdges );
  const DeconCharacteristicLimitResult checked_result
                                            = decon_characteristic_limits( checked_edges );
  BOOST_REQUIRE_MESSAGE( checked_result.status == DeconCharacteristicLimitStatus::Success,
                         checked_result.error_message );
  BOOST_CHECK( !has_control_warning(checked_result) );

  // A floating continuum is determined by the whole region rather than by a few edge channels, so
  // it stayed calibrated at 0.1 counts/channel and must not raise the control-channel warning.
  DeconCharacteristicLimitInput sparse_floating;
  sparse_floating.decon_input
        = *make_decon_input( 0.0, 0.1, DeconContinuumNorm::Floating );
  const DeconCharacteristicLimitResult sparse_floating_result
                                            = decon_characteristic_limits( sparse_floating );
  BOOST_REQUIRE_MESSAGE( sparse_floating_result.status == DeconCharacteristicLimitStatus::Success,
                         sparse_floating_result.error_message );
  BOOST_CHECK( !has_control_warning(sparse_floating_result) );

  // The user picks the sideband geometry, and the Monte Carlo only varies the counts, so the
  // warning must key on the total control counts as well as the per-channel mean.  One channel a
  // side at one count each is two control counts in total - fewer than the 2.4 the failing 4+4
  // cell had - and must be warned about even though its per-channel mean is a healthy 1.0.
  DeconCharacteristicLimitInput thin_edges;
  thin_edges.decon_input = *make_decon_input( 0.0, 1.0, DeconContinuumNorm::FixedByEdges );
  thin_edges.decon_input.roi_info[0].num_lower_side_channels = 1;
  thin_edges.decon_input.roi_info[0].num_upper_side_channels = 1;
  const DeconCharacteristicLimitResult thin_result = decon_characteristic_limits( thin_edges );
  BOOST_REQUIRE_MESSAGE( thin_result.status == DeconCharacteristicLimitStatus::Success,
                         thin_result.error_message );
  BOOST_CHECK( has_control_warning(thin_result) );

  // A region whose lines land outside it carries no activity information.  It must be dropped with
  // a warning rather than failing the request, so the regions that *do* carry information still
  // produce a limit.
  DeconCharacteristicLimitInput signal_free;
  signal_free.decon_input = *make_decon_input( 0.0, 100.0 );
  DeconRoiInfo empty_roi = signal_free.decon_input.roi_info[0];
  empty_roi.roi_start = 1388.0f;
  empty_roi.roi_end = 1412.0f;   // its 1173 keV line is nowhere near these channels
  signal_free.decon_input.roi_info.push_back( empty_roi );
  const DeconCharacteristicLimitResult signal_free_result
                                            = decon_characteristic_limits( signal_free );
  BOOST_REQUIRE_MESSAGE( signal_free_result.status == DeconCharacteristicLimitStatus::Success,
                         signal_free_result.error_message );
  BOOST_CHECK( std::isfinite(signal_free_result.detection_limit) );
  BOOST_CHECK( std::any_of( begin(signal_free_result.warnings), end(signal_free_result.warnings),
                           []( const std::string &warning ){
                             return warning.find("no counts inside it") != std::string::npos;
                           } ) );

  // A peak-like observation over an otherwise empty Floating region drives the signed primary
  // fit's continuum onto its positivity boundary.  The local-normal 1/E information must be
  // rejected explicitly; using the optimizer's tiny numerical floor as physical information would
  // fabricate an extremely small uncertainty and threshold.
  DeconCharacteristicLimitInput boundary;
  boundary.decon_input = *make_decon_input( 0.0, 0.0, DeconContinuumNorm::Floating );
  std::shared_ptr<std::vector<float>> boundary_counts
      = std::make_shared<std::vector<float>>(
                                  boundary.decon_input.measurement->num_gamma_channels(), 0.0f );
  boundary_counts->at(586) = 5.0f;
  std::shared_ptr<SpecUtils::Measurement> boundary_measurement
                                            = std::make_shared<SpecUtils::Measurement>();
  boundary_measurement->set_gamma_counts( boundary_counts, 1.0f, 1.0f );
  boundary_measurement->set_energy_calibration(
                                  boundary.decon_input.measurement->energy_calibration() );
  boundary.decon_input.measurement = boundary_measurement;
  const DeconCharacteristicLimitResult boundary_result
                                            = decon_characteristic_limits( boundary );
  BOOST_CHECK_EQUAL( static_cast<int>(boundary_result.status),
                     static_cast<int>(DeconCharacteristicLimitStatus::ApproximationInvalid) );
  BOOST_CHECK( std::isfinite(boundary_result.primary_activity) );
  BOOST_CHECK( !std::isfinite(boundary_result.decision_threshold) );
}//BOOST_AUTO_TEST_CASE( DeconCharacteristicLimitsModelsAndTypedFailures )


BOOST_AUTO_TEST_CASE( DeconCharacteristicLimitsMonteCarlo )
{
  using namespace DetectionLimitCalc;

  const char * const enabled = std::getenv( "INTERSPEC_RUN_DECON_CHARACTERISTIC_MC" );
  if( !enabled || (std::string(enabled) != "1") )
  {
    BOOST_TEST_MESSAGE( "Set INTERSPEC_RUN_DECON_CHARACTERISTIC_MC=1 to run the Lc/Ld Monte Carlo." );
    return;
  }

  size_t num_trials = 2000;
  const char * const requested_trials
                        = std::getenv( "INTERSPEC_DECON_CHARACTERISTIC_MC_TRIALS" );
  if( requested_trials )
    num_trials = static_cast<size_t>( std::stoul(requested_trials) );
  BOOST_REQUIRE_GE( num_trials, size_t(100) );

  std::mt19937_64 random( 0x11929D1CULL );
  const std::pair<const char *,DeconContinuumNorm> modes[] = {
    { "Floating", DeconContinuumNorm::Floating },
    { "FixedByEdges", DeconContinuumNorm::FixedByEdges }
  };

  for( const std::pair<const char *,DeconContinuumNorm> &mode : modes )
  {
    for( const double continuum : { 0.1, 1.0, 10.0, 100.0 } )
    {
      // The supplied sideband geometry has only a handful of control channels.  At 0.1 expected
      // count/channel, FixedByEdges failed the fixed-seed false-positive gate (0.083 versus 0.05),
      // so this local-normal Lc/Ld construction has no validated claim in that cell.  Retain the
      // 0.1-count Floating cell and validate FixedByEdges from 1 count/channel upward; callers see
      // ApproximationInvalid when an individual fitted artificial model reaches its positivity
      // boundary.
      if( (mode.second == DeconContinuumNorm::FixedByEdges) && (continuum == 0.1) )
      {
        BOOST_TEST_MESSAGE( "FixedByEdges, 0.1 count/channel excluded: below validated "
                            "local-normal applicability range." );
        continue;
      }

      DeconCharacteristicLimitInput design_input;
      design_input.decon_input = *make_decon_input( 0.0, continuum, mode.second );
      design_input.alpha = 0.05;
      design_input.beta = 0.05;
      const DeconCharacteristicLimitResult design
                                            = decon_characteristic_limits( design_input );
      BOOST_REQUIRE_MESSAGE( design.status == DeconCharacteristicLimitStatus::Success,
                             mode.first << " at " << continuum << ": " << design.error_message );

      size_t false_positives = 0, detections_at_ld = 0, failures = 0;
      std::map<DeconCharacteristicLimitStatus,size_t> failure_statuses;
      std::string first_failure;
      for( size_t truth_index = 0; truth_index < 2; ++truth_index )
      {
        const double truth = truth_index ? design.detection_limit : 0.0;
        DeconComputeInput expectation = *make_decon_input( truth, continuum, mode.second );
        const std::vector<float> &expected = *expectation.measurement->gamma_counts();

        for( size_t trial = 0; trial < num_trials; ++trial )
        {
          std::shared_ptr<std::vector<float>> observed
                              = std::make_shared<std::vector<float>>( expected.size(), 0.0f );
          for( size_t channel = 0; channel < expected.size(); ++channel )
          {
            std::poisson_distribution<unsigned int> poisson(
                                                (std::max)(0.0, double(expected[channel])) );
            observed->at(channel) = static_cast<float>( poisson(random) );
          }

          // Every channel is redrawn, including both FixedByEdges sidebands; no control counts are
          // reused between realizations.
          std::shared_ptr<SpecUtils::Measurement> spectrum
                                              = std::make_shared<SpecUtils::Measurement>();
          spectrum->set_gamma_counts( observed, 1.0f, 1.0f );
          spectrum->set_energy_calibration( expectation.measurement->energy_calibration() );

          DeconCharacteristicLimitInput realization;
          realization.decon_input = expectation;
          realization.decon_input.measurement = spectrum;
          realization.alpha = design_input.alpha;
          realization.beta = design_input.beta;
          const DeconCharacteristicLimitResult result
                                                = decon_characteristic_limits( realization );
          // Lc and Ld define one measurement procedure and are fixed by the artificial design
          // spectrum above.  Recomputing Lc for every random realization would measure the random
          // threshold A_hat > Lc(X), not the required probability P(A_hat > Lc); measured that
          // way the false-positive probability inflates to as much as 0.15, which is why
          // `decon_characteristic_limits` returns no detection boolean of its own.  Only the
          // signed primary estimate is therefore taken from each realization.  A later
          // artificial-spectrum failure (possible for sparse spectra whose fitted continuum
          // reaches its positivity boundary once A_hat is replaced by zero) does not invalidate
          // that already-completed primary fit.
          if( !std::isfinite(result.primary_activity) )
          {
            ++failures;
            ++failure_statuses[result.status];
            if( first_failure.empty() )
              first_failure = result.error_message;
            continue;
          }

          const bool detected = (result.primary_activity > design.decision_threshold);
          if( truth_index )
            detections_at_ld += detected ? 1 : 0;
          else
            false_positives += detected ? 1 : 0;
        }//for( each trial )
      }//for( null and Ld truth )

      BOOST_REQUIRE_MESSAGE( failures == 0,
                             mode.first << " at " << continuum << " had " << failures
                                        << " failed signed-primary fits; first: "
                                        << first_failure << "; status counts: "
                                        << [&failure_statuses]() {
                                             std::ostringstream out;
                                             for( const auto &entry : failure_statuses )
                                               out << static_cast<int>(entry.first) << '='
                                                   << entry.second << ' ';
                                             return out.str();
                                           }() );
      const double false_positive_probability
                          = static_cast<double>(false_positives) / static_cast<double>(num_trials);
      const double detection_power
                          = static_cast<double>(detections_at_ld) / static_cast<double>(num_trials);
      const double false_positive_se
                          = std::sqrt(0.05*0.95/static_cast<double>(num_trials));
      const double power_se = false_positive_se;
      const double false_positive_tolerance = (std::max)( 0.005, 3.0*false_positive_se );
      const double power_tolerance = (std::max)( 0.01, 3.0*power_se );

      BOOST_TEST_MESSAGE( mode.first << ", " << continuum << " counts/channel: P(false+)="
                          << false_positive_probability << ", power(Ld)=" << detection_power );
      BOOST_CHECK_LE( std::fabs(false_positive_probability - 0.05),
                      false_positive_tolerance );
      BOOST_CHECK_LE( std::fabs(detection_power - 0.95), power_tolerance );
    }//for( continuum level )
  }//for( continuum method )
}//BOOST_AUTO_TEST_CASE( DeconCharacteristicLimitsMonteCarlo )


BOOST_AUTO_TEST_CASE( BasicCurrieMda )
{
  using namespace DetectionLimitCalc;

  set_data_dir();

  auto spectrum = make_shared<SpecUtils::Measurement>();
  float liveTime = 1.0, realTime = 1.0;
  size_t num_channels = 1024;
  auto counts = make_shared<vector<float>>( num_channels, 10.0f );
  spectrum->set_gamma_counts( counts, liveTime, realTime );

  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( num_channels, { 0.0f, 3.0f }, {} );
  spectrum->set_energy_calibration( cal );

  CurrieMdaInput input;
  input.spectrum = spectrum;
  input.gamma_energy = 512.0f;    // channel: 170.667
  input.roi_lower_energy = 500.0f;// channel: 166.667
  input.roi_upper_energy = 524.0f;// channel: 174.667
  input.num_lower_side_channels = 4;
  input.num_upper_side_channels = 4;
  input.detection_probability = 0.95;
  input.additional_uncertainty = 0.0;

  CurrieMdaResult result;

  BOOST_REQUIRE_NO_THROW( result = currie_mda_calc( input ) );

  BOOST_CHECK_EQUAL( result.first_lower_continuum_channel, 163 );
  BOOST_CHECK_EQUAL( result.last_lower_continuum_channel, 166 );
  BOOST_CHECK_EQUAL( result.lower_continuum_counts_sum, 40.0f );
  BOOST_CHECK_EQUAL( result.first_upper_continuum_channel, 175 );
  BOOST_CHECK_EQUAL( result.last_upper_continuum_channel, 178 );
  BOOST_CHECK_EQUAL( result.upper_continuum_counts_sum, 40.0f );
  BOOST_CHECK_EQUAL( result.first_peak_region_channel, 167 );
  BOOST_CHECK_EQUAL( result.last_peak_region_channel, 174 );
  BOOST_CHECK_EQUAL( result.peak_region_counts_sum, 80.0f );
  BOOST_CHECK_CLOSE_FRACTION( result.continuum_eqn[0], 3.3333333333333335, 1.0E-5 );
  BOOST_CHECK_LT( fabs(result.continuum_eqn[1]), 1.0E-5 );
  BOOST_CHECK_EQUAL( result.estimated_peak_continuum_counts, 80.0f );
  BOOST_CHECK_CLOSE_FRACTION( result.estimated_peak_continuum_uncert, sqrt(80.0f), 1.0E-5 ); //4 side channels each side
  
  const double k = 1.64485;
  const double region_sigma = sqrt( 80 + 80 );
  BOOST_CHECK_CLOSE_FRACTION( result.decision_threshold, k*region_sigma, 1.0E-5 ); //continuum uncert^2 + expected amount uncert^2
  BOOST_CHECK_CLOSE_FRACTION( result.detection_limit, (2.0*k*region_sigma + k*k), 1.0E-5 );
  BOOST_CHECK_EQUAL( result.source_counts, 0.0f );
  BOOST_CHECK_CLOSE_FRACTION( result.lower_limit, -k*region_sigma, 1.0E-5 );
  BOOST_CHECK_CLOSE_FRACTION( result.upper_limit, k*region_sigma, 1.0E-5 );

  // The denominator is 1-k^2*u^2.  A relative uncertainty of 0.5 is valid at 95% even though
  // the old guard (which omitted one factor of u) rejected it; 0.7 makes the denominator invalid.
  input.additional_uncertainty = 0.5;
  BOOST_REQUIRE_NO_THROW( result = currie_mda_calc( input ) );
  BOOST_CHECK_CLOSE_FRACTION( result.detection_limit,
          (2.0*result.decision_threshold + k*k) / (1.0 - k*k*0.5*0.5), 1.0E-5 );

  input.additional_uncertainty = 0.7;
  BOOST_REQUIRE_NO_THROW( result = currie_mda_calc( input ) );
  BOOST_CHECK_EQUAL( result.detection_limit, -999.0f );
}// BOOST_AUTO_TEST_CASE( BasicCurrieMda )


/** The Currie peak-region continuum, on binning whose channel widths change across the region.

 `currie_mda_calc` used to estimate the continuum as the mean of the two side-band *densities*.
 That equals the fitted continuum at the region centre only when the two bands span equal energy,
 and the bands hold an equal number of *channels* - so any deviation-pair or lower-channel-edge
 calibration biased it by `slope*(w_upper - w_lower)/4*w_peak`.  Measured at -0.9 sigma of the
 continuum's own Poisson noise on the 81 keV line of the bundled Ba-133 NaI spectrum, whose
 deviation pairs make the two bands 11.5 and 12.9 keV wide.

 Fed the *exact* integrals of a known linear continuum, so there is a right answer to hit rather
 than a tolerance to negotiate: a width-weighted least-squares line through the side bands
 reproduces that continuum exactly, whatever the binning.
 */
BOOST_AUTO_TEST_CASE( CurrieContinuumNonUniformBinning )
{
  using namespace DetectionLimitCalc;

  set_data_dir();

  const size_t num_channels = 256;
  const size_t first_peak_channel = 108, last_peak_channel = 127;
  const size_t num_side = 4;

  // The continuum the channels are filled from, as a density in counts/keV about `reference_energy`.
  const double reference_energy = 100.0;
  const double density_at_ref = 100.0;
  const double density_slope = -5.0;
  const auto density = [&]( const double energy ) -> double {
    return density_at_ref + density_slope*(energy - reference_energy);
  };

  // Two calibrations over the same channels: one whose widths grow across the region (the case
  //  under test), and a uniform one that has to keep giving bit-for-bit what it always did.
  const auto run = [&]( const bool non_uniform ) -> CurrieMdaResult {
    vector<float> edges( num_channels + 1, 0.0f );
    for( size_t c = 0; c <= num_channels; ++c )
    {
      const double ch = static_cast<double>( c );
      // 0.5 keV/channel, with a quadratic term that widens channels by ~13% across the region.
      edges[c] = static_cast<float>( non_uniform ? (0.5*ch + 0.002*ch*ch) : (0.5*ch) );
    }

    auto cal = make_shared<SpecUtils::EnergyCalibration>();
    cal->set_lower_channel_energy( num_channels, edges );

    // Exactly the integral of `density` over each channel; for a linear density that is its value
    //  at the channel midpoint times the channel width, with no quadrature error at all.
    auto counts = make_shared<vector<float>>( num_channels, 0.0f );
    for( size_t c = 0; c < num_channels; ++c )
    {
      const double lower = edges[c], upper = edges[c+1];
      counts->at(c) = static_cast<float>( density( 0.5*(lower + upper) ) * (upper - lower) );
    }

    auto spectrum = make_shared<SpecUtils::Measurement>();
    spectrum->set_gamma_counts( counts, 1.0f, 1.0f );
    spectrum->set_energy_calibration( cal );

    CurrieMdaInput input;
    input.spectrum = spectrum;
    input.roi_lower_energy = edges[first_peak_channel];
    input.roi_upper_energy = edges[last_peak_channel + 1];
    // Deliberately off the region's centre, so the shift of `continuum_eqn` off the band centroid
    //  is exercised rather than cancelling.
    input.gamma_energy = static_cast<float>( input.roi_lower_energy
                            + 0.4*(input.roi_upper_energy - input.roi_lower_energy) );
    input.num_lower_side_channels = num_side;
    input.num_upper_side_channels = num_side;
    input.detection_probability = 0.95;
    input.additional_uncertainty = 0.0;

    CurrieMdaResult result;
    BOOST_REQUIRE_NO_THROW( result = currie_mda_calc( input ) );
    BOOST_REQUIRE_EQUAL( result.first_peak_region_channel, first_peak_channel );
    BOOST_REQUIRE_EQUAL( result.last_peak_region_channel, last_peak_channel );

    const double peak_lower = spectrum->gamma_channel_lower( first_peak_channel );
    const double peak_upper = spectrum->gamma_channel_upper( last_peak_channel );
    const double peak_width = peak_upper - peak_lower;

    // What the continuum actually puts in the peak region: a line integrated over a range is its
    //  value at the range's centre times the width.
    const double analytic = density( 0.5*(peak_lower + peak_upper) ) * peak_width;

    // What averaging the two side-band densities would have given.
    const double lower_width = spectrum->gamma_channel_upper( result.last_lower_continuum_channel )
                               - spectrum->gamma_channel_lower( result.first_lower_continuum_channel );
    const double upper_width = spectrum->gamma_channel_upper( result.last_upper_continuum_channel )
                               - spectrum->gamma_channel_lower( result.first_upper_continuum_channel );
    const double mean_of_densities = 0.5*( (result.lower_continuum_counts_sum / lower_width)
                                          + (result.upper_continuum_counts_sum / upper_width) )
                                     * peak_width;

    BOOST_TEST_MESSAGE( (non_uniform ? "non-uniform" : "uniform")
        << " binning: lower band " << lower_width << " keV, upper band " << upper_width << " keV;"
        << " analytic " << analytic << ", fit " << result.estimated_peak_continuum_counts
        << ", mean-of-densities " << mean_of_densities
        << " (" << 100.0*(mean_of_densities - analytic)/analytic << "%)" );

    // The fit has to land on the continuum it was given, whatever the binning.
    BOOST_CHECK_CLOSE_FRACTION( static_cast<double>(result.estimated_peak_continuum_counts),
                               analytic, 1.0E-4 );

    // And `continuum_eqn` has to be the same line - it is what the MDA chart draws under the peak,
    //  and the two used to disagree with each other on exactly this binning.
    const double eqn_at_peak_mid = result.continuum_eqn[0]
        + result.continuum_eqn[1]*(0.5*(peak_lower + peak_upper) - input.gamma_energy);
    BOOST_CHECK_CLOSE_FRACTION( eqn_at_peak_mid*peak_width, analytic, 1.0E-4 );

    if( non_uniform )
    {
      // States the size of the defect this case pins: the bands differ in width by more than 10%,
      //  and the old estimator was wrong by more than 0.2% of the continuum because of it.
      BOOST_CHECK_GT( upper_width/lower_width, 1.10 );
      BOOST_CHECK_GT( fabs(mean_of_densities - analytic)/analytic, 2.0E-3 );
    }else
    {
      // The reduction property: with uniform channels and an equal number of side channels each
      //  side, the fit *is* the mean of the two densities.  Every production caller is in this
      //  case, so nothing about the change may move here.
      BOOST_CHECK_CLOSE_FRACTION( static_cast<double>(result.estimated_peak_continuum_counts),
                                 mean_of_densities, 1.0E-6 );
    }

    return result;
  };//run(...)

  run( true );
  run( false );
}// BOOST_AUTO_TEST_CASE( CurrieContinuumNonUniformBinning )


BOOST_AUTO_TEST_CASE( Table16OfAQ48 )
{
  //Example given on page 68 of:
  //  INTERNATIONAL ATOMIC ENERGY AGENCY, Determination and Interpretation of
  //  Characteristic Limits for Radioactivity Measurements, IAEA Analytical Quality
  //  in Nuclear Applications Series No. 48, IAEA, Vienna (2017)
  using namespace DetectionLimitCalc;
  
  set_data_dir();
  
  auto spectrum = make_shared<SpecUtils::Measurement>();
  float liveTime = 1.0, realTime = 1.0;
  size_t num_channels = 1024;
  auto counts = make_shared<vector<float>>(num_channels, 0.0f);
  (*counts)[992]  = 59;
  (*counts)[993]  = 56;
  (*counts)[994]  = 74;
  (*counts)[995]  = 72;
  (*counts)[996]  = 64;
  (*counts)[997]  = 46; //First channel of ROI
  (*counts)[998]  = 64;
  (*counts)[999]  = 67;
  (*counts)[1000] = 82;
  (*counts)[1001] = 95;
  (*counts)[1002] = 157;
  (*counts)[1003] = 398;
  (*counts)[1004] = 807;
  (*counts)[1005] = 1480;
  (*counts)[1006] = 1814;
  (*counts)[1007] = 1936;
  (*counts)[1008] = 1575;
  (*counts)[1009] = 940;
  (*counts)[1010] = 457;
  (*counts)[1011] = 207;
  (*counts)[1012] = 82;
  (*counts)[1013] = 49;
  (*counts)[1014] = 50;
  (*counts)[1015] = 45;
  (*counts)[1016] = 43; //Last channel of ROI
  (*counts)[1017] = 45;
  (*counts)[1018] = 51;
  (*counts)[1019] = 35;
  (*counts)[1020] = 45;
  (*counts)[1021] = 50;
  
  spectrum->set_gamma_counts( counts, liveTime, realTime );
  
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( num_channels, {0.0f, 1.0f}, {} );
  spectrum->set_energy_calibration( cal );
  
  CurrieMdaInput input;
  input.spectrum = spectrum;
  input.gamma_energy = 1007.0f;
  input.roi_lower_energy = 997.0f;
  input.roi_upper_energy = 1017.0f;
  input.num_lower_side_channels = 5;
  input.num_upper_side_channels = 5;
  input.detection_probability = 0.95;
  input.additional_uncertainty = 0.0;
  
  CurrieMdaResult result;
  
  BOOST_REQUIRE_NO_THROW( result = currie_mda_calc( input ) );
  
  BOOST_CHECK_EQUAL( result.peak_region_counts_sum, 10394 );
  BOOST_CHECK_EQUAL( result.lower_continuum_counts_sum, 325 );  //n_1
  BOOST_CHECK_EQUAL( result.upper_continuum_counts_sum, 226 );  //n_2
  BOOST_CHECK_EQUAL( result.lower_continuum_counts_sum + result.upper_continuum_counts_sum, 551 );  //n_0
  BOOST_CHECK_EQUAL( 1 + result.last_lower_continuum_channel - result.first_lower_continuum_channel, 5 );
  BOOST_CHECK_EQUAL( 1 + result.last_upper_continuum_channel - result.first_upper_continuum_channel, 5 );
  BOOST_CHECK_EQUAL( 1 + result.last_peak_region_channel - result.first_peak_region_channel, 20 ); //n_g - num channels in peak reagion
  BOOST_CHECK_EQUAL( result.estimated_peak_continuum_counts, 1102.0f );
  const double w = 0.002458; //weight: counts to bq
  //const double w_uncert_rel = 0.0478; //We dont take this into account, yet
  const double k = 1.64485;
  const double a_star = w * result.decision_threshold;
  BOOST_CHECK_CLOSE_FRACTION( a_star, 0.232506, 1.0E-3 ); //We get: 0.232467
  
  const double a_pound = w * result.detection_limit;
  BOOST_CHECK_CLOSE_FRACTION( a_pound, 0.471705, 1.0E-3 ); //We get: 0.471583
}//BOOST_AUTO_TEST_CASE( Table16OfAQ48 )
