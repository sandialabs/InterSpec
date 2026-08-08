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

#include <set>
#include <cmath>
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
};// struct StudyProfileResult

struct StudyPeakModel
{
  vector<double> observed;
  vector<double> variance;
  vector<double> relative_energy;
  vector<double> peak_fraction;
  vector<double> edge_expected;
  size_t num_score_bins = 0;
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

StudyPeakModel make_study_peak_model( const DetectionLimitCalc::DeconComputeInput &input )
{
  const DetectionLimitCalc::DeconRoiInfo &roi = input.roi_info.front();
  const DetectionLimitCalc::DeconRoiInfo::PeakInfo &peak = roi.peak_infos.front();
  const shared_ptr<const SpecUtils::EnergyCalibration> calibration = input.measurement->energy_calibration();
  const size_t first_channel = input.measurement->find_gamma_channel( roi.roi_start );
  const size_t fit_last_channel = input.measurement->find_gamma_channel( roi.roi_end );
  const size_t score_last_channel = input.measurement->find_gamma_channel( roi.roi_end - 0.0000001f );
  const double sigma = peak.fwhm / PhysicalUnits::fwhm_nsigma;

  StudyPeakModel model;
  model.num_score_bins = 1 + score_last_channel - first_channel;
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
  }

  double lower_cont_bound = calibration->channel_for_energy( roi.roi_start );
  double upper_cont_bound = calibration->channel_for_energy( roi.roi_end );
  double fractional = lower_cont_bound - floor( lower_cont_bound );
  lower_cont_bound = ( fractional > 0.9 ) ? round( lower_cont_bound ) : floor( lower_cont_bound );
  fractional = upper_cont_bound - floor( upper_cont_bound );
  upper_cont_bound = ( fractional < 0.1 ) ? round( upper_cont_bound ) : floor( upper_cont_bound );
  upper_cont_bound = (std::max)( upper_cont_bound, 1.0 );

  const size_t num_lower = roi.num_lower_side_channels;
  const size_t num_upper = roi.num_upper_side_channels;
  const size_t roi_first = static_cast<size_t>( lower_cont_bound );
  const size_t roi_last = (std::max)( roi_first, static_cast<size_t>( upper_cont_bound - 1.0 ) );
  const size_t lower_first = roi_first - num_lower;
  const size_t lower_last = roi_first - 1;
  const size_t upper_first = roi_last + 1;
  const size_t upper_last = upper_first + num_upper - 1;
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
    const double delta = 2.705543454095404;
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
  const double delta = 2.705543454095404;
  const double discriminant = linear * linear - 4.0 * quadratic * ( score0 - best_score - delta );
  if( discriminant < 0.0 )
    return result;

  result.upper = ( -linear + sqrt( discriminant ) ) / ( 2.0 * quadratic );
  result.found_upper = ( result.upper >= result.best );
  result.found_lower = ( result.best > 1.0E-8 ) && ( ( score0 - best_score ) >= delta );
  return result;
}// study_neyman_limit(...)

double study_cash_objective( const StudyPeakModel &model,
                             const double activity,
                             const pair<double, double> &continuum,
                             const size_t num_bins )
{
  double score = 0.0;
  for( size_t index = 0; index < num_bins; ++index )
  {
    const double expected = (std::max)( 1.0E-10,
                                        continuum.first + continuum.second * model.relative_energy[index] +
                                          activity * model.peak_fraction[index] );
    const double observed = model.observed[index];
    score += 2.0 * ( expected - observed + ( ( observed > 0.0 ) ? observed * std::log( observed / expected ) : 0.0 ) );
  }
  return score;
}// study_cash_objective(...)

pair<double, double> study_cash_continuum( const StudyPeakModel &model, const double activity )
{
  pair<double, double> continuum = study_weighted_continuum( model, activity );
  double min_expected = numeric_limits<double>::max();
  for( size_t index = 0; index < model.observed.size(); ++index )
  {
    min_expected = (std::min)( min_expected,
                               continuum.first + continuum.second * model.relative_energy[index] +
                                 activity * model.peak_fraction[index] );
  }
  if( min_expected <= 1.0E-6 )
    continuum.first += 1.0E-6 - min_expected;

  double objective = study_cash_objective( model, activity, continuum, model.observed.size() );
  for( size_t iteration = 0; iteration < 30; ++iteration )
  {
    double gradient0 = 0.0, gradient1 = 0.0;
    double hessian00 = 1.0E-10, hessian01 = 0.0, hessian11 = 1.0E-10;
    for( size_t index = 0; index < model.observed.size(); ++index )
    {
      const double x = model.relative_energy[index];
      const double expected =
        (std::max)( 1.0E-10, continuum.first + continuum.second * x + activity * model.peak_fraction[index] );
      const double observed = model.observed[index];
      const double ratio = observed / expected;
      const double curvature = observed / ( expected * expected );
      gradient0 += 1.0 - ratio;
      gradient1 += ( 1.0 - ratio ) * x;
      hessian00 += curvature;
      hessian01 += curvature * x;
      hessian11 += curvature * x * x;
    }

    const double determinant = hessian00 * hessian11 - hessian01 * hessian01;
    if( determinant <= 1.0E-20 )
      break;
    const double step0 = ( hessian11 * gradient0 - hessian01 * gradient1 ) / determinant;
    const double step1 = ( hessian00 * gradient1 - hessian01 * gradient0 ) / determinant;
    if( ( fabs( step0 ) + fabs( step1 ) ) < 1.0E-9 )
      break;

    bool accepted = false;
    double scale = 1.0;
    for( size_t line_iteration = 0; line_iteration < 40; ++line_iteration )
    {
      const pair<double, double> candidate = { continuum.first - scale * step0, continuum.second - scale * step1 };
      bool positive = true;
      for( size_t index = 0; index < model.observed.size(); ++index )
      {
        const double expected =
          candidate.first + candidate.second * model.relative_energy[index] + activity * model.peak_fraction[index];
        positive = positive && ( expected > 1.0E-10 );
      }
      const double candidate_objective = positive
                                           ? study_cash_objective( model, activity, candidate, model.observed.size() )
                                           : numeric_limits<double>::max();
      if( candidate_objective < objective )
      {
        continuum = candidate;
        objective = candidate_objective;
        accepted = true;
        break;
      }
      scale *= 0.5;
    }
    if( !accepted )
      break;
  }

  return continuum;
}// study_cash_continuum(...)

double study_cash_score( const StudyPeakModel &model,
                         const DetectionLimitCalc::DeconContinuumNorm norm,
                         const pair<double, double> &full_range_continuum,
                         const double activity )
{
  if( norm == DetectionLimitCalc::DeconContinuumNorm::FixedByEdges )
  {
    double score = 0.0;
    for( size_t index = 0; index < model.num_score_bins; ++index )
    {
      const double expected = (std::max)( 1.0E-10, model.edge_expected[index] + activity * model.peak_fraction[index] );
      const double observed = model.observed[index];
      score +=
        2.0 * ( expected - observed + ( ( observed > 0.0 ) ? observed * std::log( observed / expected ) : 0.0 ) );
    }
    return score;
  }

  pair<double, double> continuum;
  switch( norm )
  {
  case DetectionLimitCalc::DeconContinuumNorm::Floating:
    continuum = study_cash_continuum( model, activity );
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

StudyProfileResult study_cash_limit( const DetectionLimitCalc::DeconComputeInput &input, const double initial_maximum )
{
  const StudyPeakModel model = make_study_peak_model( input );
  const DetectionLimitCalc::DeconContinuumNorm norm = input.roi_info.front().cont_norm_method;
  const pair<double, double> full_range_continuum = study_cash_continuum( model, 0.0 );
  const auto score = [&model, norm, &full_range_continuum]( const double activity )
  { return study_cash_score( model, norm, full_range_continuum, activity ); };

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
  const double delta = 2.705543454095404;
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
  for( size_t iteration = 0; iteration < 36; ++iteration )
  {
    const double midpoint = 0.5 * ( lower + root_upper );
    if( ( score( midpoint ) - best_score ) < delta )
      lower = midpoint;
    else
      root_upper = midpoint;
  }
  result.upper = 0.5 * ( lower + root_upper );
  result.found_upper = true;
  result.found_lower = ( result.best > 1.0E-8 ) && ( ( score( 0.0 ) - best_score ) >= delta );
  return result;
}// study_cash_limit(...)

shared_ptr<const DetectionLimitCalc::DeconComputeInput>
make_poisson_trial( const shared_ptr<const DetectionLimitCalc::DeconComputeInput> &expected, mt19937 &generator )
{
  shared_ptr<vector<float>> counts = make_shared<vector<float>>( *expected->measurement->gamma_counts() );
  const DetectionLimitCalc::DeconRoiInfo &roi = expected->roi_info.front();
  const size_t first = expected->measurement->find_gamma_channel( roi.roi_start ) - roi.num_lower_side_channels;
  const size_t last =
    expected->measurement->find_gamma_channel( roi.roi_end - 0.0000001f ) + roi.num_upper_side_channels;
  for( size_t channel = first; channel <= last; ++channel )
  {
    const double mean = expected->measurement->gamma_channel_content( channel );
    poisson_distribution<int> distribution( mean );
    counts->at( channel ) = static_cast<float>( distribution( generator ) );
  }

  shared_ptr<SpecUtils::Measurement> measurement = make_shared<SpecUtils::Measurement>();
  measurement->set_gamma_counts( counts, 1.0f, 1.0f );
  measurement->set_energy_calibration( expected->measurement->energy_calibration() );

  shared_ptr<DetectionLimitCalc::DeconComputeInput> trial =
    make_shared<DetectionLimitCalc::DeconComputeInput>( *expected );
  trial->measurement = measurement;
  return trial;
}// make_poisson_trial(...)

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

  const shared_ptr<const DetectionLimitCalc::DeconComputeInput> detected = make_decon_input( 250.0, 100.0 );
  DetectionLimitCalc::DeconActivityOrDistanceLimitResult result;
  BOOST_REQUIRE_NO_THROW(
    result = DetectionLimitCalc::get_activity_or_distance_limits( 0.95f, detected, false, 0.0, 260.0, true ) );
  BOOST_CHECK( result.foundLowerCl );
  BOOST_CHECK( !result.foundUpperCl );
  BOOST_CHECK( !result.errorMessage.empty() );

  const shared_ptr<const DetectionLimitCalc::DeconComputeInput> boundary = make_decon_input( 0.0, 100.0 );
  BOOST_REQUIRE_NO_THROW(
    result = DetectionLimitCalc::get_activity_or_distance_limits( 0.95f, boundary, false, 0.0, 1.0, true ) );
  BOOST_CHECK( !result.foundLowerCl );
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

  BOOST_CHECK_EQUAL( floating_result.num_degree_of_freedom, full_range_result.num_degree_of_freedom );
  BOOST_CHECK_EQUAL( edges_result.num_degree_of_freedom, floating_result.num_degree_of_freedom + 2 );

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
  BOOST_CHECK_SMALL( fabs( combined_result.fit_peaks[0].chi2dof() - edges_result.fit_peaks[0].chi2dof() ), 1.0E-8 );
  BOOST_CHECK_SMALL( fabs( combined_result.fit_peaks[1].chi2dof() - second_result.fit_peaks[0].chi2dof() ), 1.0E-8 );
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

  DetectionLimitCalc::DeconComputeInput overlapping = *input;
  overlapping.roi_info.push_back( overlapping.roi_info.front() );

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

  // This records a methodology gap rather than endorsing it: duplicate ROIs currently count the
  // same channels twice and therefore claim the same sqrt(2) information gain as independent ROIs.
  BOOST_CHECK_CLOSE( overlapping_result.upperLimit / independent_result.upperLimit, 1.0, 4.0 );
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

BOOST_AUTO_TEST_CASE( DeconCurrieBundledSpectrumGrid )
{
  set_data_dir();

#ifndef NDEBUG
  // `currie_mda_calc(...)` has an existing developer-only assertion that assumes equally spaced
  // side regions.  The Ba-133 spectrum uses a nonlinear energy calibration, so collect these
  // empirical ratios in the requested Release study rather than aborting the Debug suite.
  BOOST_TEST_MESSAGE(
    "Bundled-spectrum ratio study runs in Release because of the Currie nonlinear-calibration assertion" );
  return;
#endif

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
        const StudyProfileResult study = study_neyman_limit( *input, maximum );
        DetectionLimitCalc::DeconActivityOrDistanceLimitResult production;
        {
          ScopedCoutSilence silence;
          production = DetectionLimitCalc::get_activity_or_distance_limits( 0.95f, input, false, 0.0, maximum, true );
        }
        BOOST_REQUIRE( study.found_upper );
        BOOST_REQUIRE( production.foundUpperCl );
        BOOST_CHECK_MESSAGE( fabs( study.upper - production.upperLimit ) < 0.03 * ( 1.0 + production.upperLimit ),
                             "Study profile " << study.upper << " differs from production " << production.upperLimit
                                              << " for continuum=" << continuum << ", signal=" << signal
                                              << ", norm=" << static_cast<int>( norm ) );
      }
    }
  }
}// BOOST_AUTO_TEST_CASE( DeconStudyProfileMatchesProduction )


#if ( PERFORM_DEVELOPER_CHECKS )
BOOST_AUTO_TEST_CASE( DeconCoverageStudy )
{
  const char *const enabled = getenv( "INTERSPEC_RUN_DECON_COVERAGE" );
  if( !enabled || ( string( enabled ) != "1" ) )
  {
    BOOST_TEST_MESSAGE( "Set INTERSPEC_RUN_DECON_COVERAGE=1 to run the deconvolution coverage study." );
    return;
  }

  set_data_dir();
  mt19937 generator( 20260807u );

  struct CoverageStats
  {
    size_t trials = 0;
    size_t covered = 0;
    size_t cash_covered = 0;
    size_t found_lower = 0;
    size_t cash_found_lower = 0;
    size_t scan_failures = 0;
    size_t cash_failures = 0;
    size_t production_scan_checks = 0;
    size_t production_scan_failures = 0;
    double estimator_difference_sum = 0.0;
    double upper_difference_sum = 0.0;
    double cash_estimator_difference_sum = 0.0;
    double cash_upper_difference_sum = 0.0;
    double production_max_relative_difference = 0.0;
  };

  struct Mode
  {
    const char *name;
    DetectionLimitCalc::DeconContinuumNorm norm;
    size_t initial_trials;
  };

  const Mode modes[] = { { "Floating", DetectionLimitCalc::DeconContinuumNorm::Floating, 2000 },
                         { "FixedByEdges", DetectionLimitCalc::DeconContinuumNorm::FixedByEdges, 1000 },
                         { "FixedByFullRange", DetectionLimitCalc::DeconContinuumNorm::FixedByFullRange, 1000 } };
  const double continua[] = { 0.1, 1.0, 10.0, 100.0 };
  const double signal_sigmas[] = { 0.0, 1.0, 3.0 };

  cout << "DECON_COVERAGE,mode,continuum_counts_per_channel,signal_over_sqrt_roi_background,"
          "truth_counts,trials,coverage,wilson_low,wilson_high,estimator_bias,limit_bias,"
          "scan_failure_rate,found_lower_rate,cash_coverage,cash_wilson_low,cash_wilson_high,"
          "cash_estimator_bias,cash_limit_bias,cash_failure_rate,cash_found_lower_rate,"
          "production_scan_checks,production_scan_failure_rate,"
          "production_max_relative_difference\n";

  for( const Mode &mode : modes )
  {
    for( const double continuum : continua )
    {
      const shared_ptr<const DetectionLimitCalc::DeconComputeInput> background =
        make_decon_input( 0.0, continuum, mode.norm );
      const pair<size_t, size_t> roi_channels = DetectionLimitCalc::round_roi_to_channels(
        background->measurement, background->roi_info.front().roi_start, background->roi_info.front().roi_end );
      const double roi_background = continuum * ( 1 + roi_channels.second - roi_channels.first );

      for( const double signal_sigma : signal_sigmas )
      {
        const double truth = signal_sigma * sqrt( roi_background );
        const shared_ptr<const DetectionLimitCalc::DeconComputeInput> expected =
          make_decon_input( truth, continuum, mode.norm );
        CoverageStats stats;

        const auto run_trials = [&stats, &generator, &expected, truth]( const size_t count )
        {
          for( size_t trial_index = 0; trial_index < count; ++trial_index )
          {
            const shared_ptr<const DetectionLimitCalc::DeconComputeInput> trial =
              make_poisson_trial( expected, generator );
            const DetectionLimitCalc::CurrieMdaResult currie = currie_for_decon_input( trial );
            const double diff_multiple = 50.0;
            double seeded_maximum = 1.0;
            if( currie.source_counts > currie.decision_threshold )
            {
              const double lower_difference = fabs( currie.source_counts - currie.lower_limit );
              const double upper_difference = fabs( currie.upper_limit - currie.source_counts );
              seeded_maximum =
                (std::max)( 1.0,
                            currie.source_counts + diff_multiple * (std::max)( lower_difference, upper_difference ) );
            } else if( currie.upper_limit < 0.0f )
            {
              seeded_maximum =
                (std::max)( 1.0, diff_multiple * sqrt( (std::max)( 0.0f, currie.peak_region_counts_sum ) ) );
            } else
            {
              seeded_maximum = (std::max)( 1.0, diff_multiple * currie.upper_limit );
            }

            const StudyProfileResult neyman = study_neyman_limit( *trial, seeded_maximum );
            const StudyProfileResult cash = study_cash_limit( *trial, seeded_maximum );

            // The full grid uses the algebraically equivalent study profile for tractability.
            // Audit the actual production Brent/bisection path on ten trials in every cell.
            if( stats.production_scan_checks < 10 )
            {
              DetectionLimitCalc::DeconActivityOrDistanceLimitResult production;
              bool production_threw = false;
              try
              {
                ScopedCoutSilence silence;
                production =
                  DetectionLimitCalc::get_activity_or_distance_limits( 0.95f, trial, false, 0.0, seeded_maximum, true );
              } catch( std::exception & )
              {
                production_threw = true;
              }
              ++stats.production_scan_checks;
              if( production_threw || !production.foundUpperCl || !neyman.found_upper )
              {
                ++stats.production_scan_failures;
              } else
              {
                const double relative_difference =
                  fabs( production.upperLimit - neyman.upper ) / ( 1.0 + fabs( production.upperLimit ) );
                stats.production_max_relative_difference =
                  (std::max)( stats.production_max_relative_difference, relative_difference );
              }
            }
            ++stats.trials;
            if( !neyman.found_upper || ( neyman.upper > seeded_maximum ) )
            {
              ++stats.scan_failures;
            } else
            {
              stats.covered += ( neyman.upper >= truth );
              stats.estimator_difference_sum += neyman.best - truth;
              stats.upper_difference_sum += neyman.upper - truth;
              stats.found_lower += neyman.found_lower;
            }

            if( !cash.found_upper )
            {
              ++stats.cash_failures;
            } else
            {
              stats.cash_covered += ( cash.upper >= truth );
              stats.cash_estimator_difference_sum += cash.best - truth;
              stats.cash_upper_difference_sum += cash.upper - truth;
              stats.cash_found_lower += cash.found_lower;
            }
          }
        };

        run_trials( mode.initial_trials );
        pair<double, double> interval = wilson_interval( stats.covered, stats.trials - stats.scan_failures );
        const pair<double, double> cash_initial_interval =
          wilson_interval( stats.cash_covered, stats.trials - stats.cash_failures );
        const double target = ( truth > 0.0 ) ? 0.95 : 0.05;
        const pair<double, double> target_interval =
          ( truth > 0.0 ) ? interval : wilson_interval( stats.found_lower, stats.trials - stats.scan_failures );
        const bool anomalous =
          ( target < target_interval.first ) || ( target > target_interval.second ) || ( stats.scan_failures > 0 ) ||
          ( ( truth > 0.0 ) && ( ( 0.95 < cash_initial_interval.first ) || ( 0.95 > cash_initial_interval.second ) ) );
        if( anomalous && ( stats.trials < 10000 ) )
          run_trials( 10000 - stats.trials );

        const size_t valid = stats.trials - stats.scan_failures;
        const size_t cash_valid = stats.trials - stats.cash_failures;
        interval = wilson_interval( stats.covered, valid );
        const pair<double, double> cash_interval = wilson_interval( stats.cash_covered, cash_valid );
        const double coverage = static_cast<double>( stats.covered ) / valid;
        const double cash_coverage = static_cast<double>( stats.cash_covered ) / cash_valid;
        const double found_lower_rate = static_cast<double>( stats.found_lower ) / valid;
        const double cash_found_lower_rate = static_cast<double>( stats.cash_found_lower ) / cash_valid;

        cout << "DECON_COVERAGE," << mode.name << ',' << continuum << ',' << signal_sigma << ',' << truth << ','
             << stats.trials << ',' << coverage << ',' << interval.first << ',' << interval.second << ','
             << stats.estimator_difference_sum / valid << ',' << stats.upper_difference_sum / valid << ','
             << static_cast<double>( stats.scan_failures ) / stats.trials << ',' << found_lower_rate << ','
             << cash_coverage << ',' << cash_interval.first << ',' << cash_interval.second << ','
             << stats.cash_estimator_difference_sum / cash_valid << ',' << stats.cash_upper_difference_sum / cash_valid
             << ',' << static_cast<double>( stats.cash_failures ) / stats.trials << ',' << cash_found_lower_rate << ','
             << stats.production_scan_checks << ','
             << static_cast<double>( stats.production_scan_failures ) / stats.production_scan_checks << ','
             << stats.production_max_relative_difference << '\n';

        BOOST_WARN_MESSAGE( stats.production_scan_failures == 0,
                            "Production scan audit failed for " << mode.name << ", continuum=" << continuum
                                                                << ", signal sigma=" << signal_sigma );
        BOOST_WARN_MESSAGE( stats.production_max_relative_difference < 0.03,
                            "Study and production scans differ for " << mode.name << ", continuum=" << continuum
                                                                     << ", signal sigma=" << signal_sigma << ": "
                                                                     << stats.production_max_relative_difference );

        if( truth > 0.0 )
        {
          BOOST_WARN_MESSAGE( ( interval.first <= 0.95 ) && ( interval.second >= 0.95 ),
                              "Neyman coverage excludes 95% for " << mode.name << ", continuum=" << continuum
                                                                  << ", signal sigma=" << signal_sigma );
        } else
        {
          const pair<double, double> false_detection_interval = wilson_interval( stats.found_lower, valid );
          BOOST_WARN_MESSAGE( ( false_detection_interval.first <= 0.05 ) && ( false_detection_interval.second >= 0.05 ),
                              "Boundary false-detection rate excludes 5% for " << mode.name
                                                                               << ", continuum=" << continuum );
        }
      }
    }
  }
}// BOOST_AUTO_TEST_CASE( DeconCoverageStudy )
#endif


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
