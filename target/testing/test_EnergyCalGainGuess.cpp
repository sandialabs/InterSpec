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
#include <random>
#include <iostream>
#include <algorithm>

#if( defined(WIN32) )
#undef min
#undef max
#endif

#define BOOST_TEST_MODULE EnergyCalGainGuess_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/EnergyCalibration.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/EnergyCalGainGuess.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace std;
using namespace boost::unit_test;

using SpecUtils::Measurement;
using SpecUtils::EnergyCalibration;

using namespace GainGuessCalc;

string g_data_dir;
string g_test_file_dir;


void set_data_dir()
{
  static bool s_have_set = false;
  if( s_have_set )
    return;

  s_have_set = true;

  int argc = boost::unit_test::framework::master_test_suite().argc;
  char **argv = boost::unit_test::framework::master_test_suite().argv;

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
    for( const auto &d : { "data", "../data", "../../data", "../../../data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path( d, "sandia.decay.xml" ) ) )
      {
        datadir = d;
        break;
      }
    }
  }

  const string sandia_decay_file = SpecUtils::append_path( datadir, "sandia.decay.xml" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( sandia_decay_file ),
    "Could not find 'sandia.decay.xml' at '" << sandia_decay_file << "'" );

  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );

  g_data_dir = datadir;
}//void set_data_dir()


namespace
{
  struct SyntheticPeak
  {
    double energy;    //keV
    double area;      //counts
    double fwhm;      //keV
  };


  /** A nominal HPGe FWHM (keV): ~1.05 at 122 keV, ~1.9 at 1332 keV. */
  double hpge_fwhm( const double energy )
  {
    return sqrt( 0.9 + 0.002*energy );
  }


  /** Builds a synthetic spectrum with the given true linear calibration, Gaussian peaks on an
   exponential continuum, and Poisson noise (fixed seed, so tests are deterministic).
   */
  shared_ptr<Measurement> make_synthetic_spectrum( const size_t nchan,
                                                   const double offset, const double gain,
                                                   const vector<SyntheticPeak> &peaks,
                                                   const double continuum_amp = 150.0,
                                                   const double continuum_decay_kev = 500.0 )
  {
    std::mt19937 rng( 987654321u );

    auto counts = make_shared<vector<float>>( nchan, 0.0f );
    for( size_t i = 0; i < nchan; ++i )
    {
      const double energy = offset + gain*(static_cast<double>(i) + 0.5);
      double mu = continuum_amp*exp( -energy/continuum_decay_kev ) + 4.0;

      for( const SyntheticPeak &peak : peaks )
      {
        const double sigma = peak.fwhm / 2.3548;
        const double z = (energy - peak.energy) / sigma;
        if( fabs(z) < 8.0 )
          mu += gain * peak.area * exp( -0.5*z*z ) / ( sigma*sqrt(2.0*M_PI) );
      }

      std::poisson_distribution<long> pois( mu );
      (*counts)[i] = static_cast<float>( pois(rng) );
    }//for( channels )

    auto meas = make_shared<Measurement>();
    meas->set_gamma_counts( counts, 600.0f, 600.0f );

    // The spectrum "arrives" with a default/nonsense calibration - a 0 to 3 MeV linear guess -
    //  like a file with no calibration info would.
    auto cal = make_shared<EnergyCalibration>();
    cal->set_default_polynomial( nchan, { 0.0f, 3000.0f/static_cast<float>(nchan) }, {} );
    meas->set_energy_calibration( cal );

    return meas;
  }//make_synthetic_spectrum(...)


  /** The energy the guessed cal predicts for a channel vs what the truth predicts. */
  void check_cal_close( const GuessedCal &guess, const double true_offset,
                        const double true_gain, const size_t nchan,
                        const double frac_tol, const double offset_tol )
  {
    BOOST_CHECK_MESSAGE( fabs( guess.gain - true_gain ) < frac_tol*true_gain,
        "gain " << guess.gain << " not within " << 100*frac_tol << "% of " << true_gain
        << " (provenance: " << guess.provenance << ")" );
    BOOST_CHECK_MESSAGE( fabs( guess.offset - true_offset ) < offset_tol,
        "offset " << guess.offset << " not within " << offset_tol << " keV of " << true_offset
        << " (provenance: " << guess.provenance << ")" );

    const double e_true = true_offset + true_gain*static_cast<double>(nchan);
    const double e_guess = guess.offset + guess.gain*static_cast<double>(nchan);
    BOOST_CHECK_MESSAGE( fabs(e_guess - e_true) < frac_tol*e_true + offset_tol,
        "top-of-spectrum energy " << e_guess << " vs true " << e_true );
  }//check_cal_close(...)


  bool has_nuclide( const GuessedCal &guess, const string &symbol )
  {
    for( const SandiaDecay::Nuclide * const nuc : guess.nuclides )
      if( nuc && (nuc->symbol == symbol) )
        return true;
    return false;
  }


  /** Cs-137 + K-40 + the major Th-232 chain lines (Pb-212, Tl-208, Ac-228) - a physically
   plausible "Cs-137 plus NORM background" HPGe spectrum.
   */
  vector<SyntheticPeak> cs137_norm_peaks()
  {
    return vector<SyntheticPeak>{
      { 238.6,   9000.0, hpge_fwhm(238.6) },   //Pb-212
      { 338.3,   1300.0, hpge_fwhm(338.3) },   //Ac-228
      { 583.2,   3000.0, hpge_fwhm(583.2) },   //Tl-208
      { 661.7,  20000.0, hpge_fwhm(661.7) },   //Cs-137
      { 911.2,   2600.0, hpge_fwhm(911.2) },   //Ac-228
      { 969.0,   1500.0, hpge_fwhm(969.0) },   //Ac-228
      { 1460.8,  6000.0, hpge_fwhm(1460.8) },  //K-40
      { 2614.5,  2500.0, hpge_fwhm(2614.5) },  //Tl-208
    };
  }
}//namespace


BOOST_AUTO_TEST_CASE( SyntheticHpgeRecovery )
{
  set_data_dir();

  const size_t nchan = 8192;
  const double true_offset = 3.0, true_gain = 0.35;

  // Cs-137 + NORM (K-40 and the Th-232 chain): the bread-and-butter case.
  const shared_ptr<Measurement> meas
              = make_synthetic_spectrum( nchan, true_offset, true_gain, cs137_norm_peaks() );

  GuessOptions options;
  const vector<GuessedCal> results = guess_energy_cal( meas, options );

  BOOST_REQUIRE_MESSAGE( !results.empty(), "No calibration guesses returned" );

  const GuessedCal &best = results.front();
  check_cal_close( best, true_offset, true_gain, nchan, 0.005, 3.0 );

  BOOST_CHECK_MESSAGE( has_nuclide( best, "Cs137" ), "Cs137 not in implied nuclides" );
  BOOST_CHECK_MESSAGE( has_nuclide( best, "K40" ), "K40 not in implied nuclides" );

  BOOST_CHECK_GT( best.explained_frac, sm_min_confident_explained_frac );
}//BOOST_AUTO_TEST_CASE( SyntheticHpgeRecovery )


BOOST_AUTO_TEST_CASE( Wide511Strategy )
{
  set_data_dir();

  const size_t nchan = 8192;
  const double true_offset = -1.0, true_gain = 0.25;

  // Na-22: a Doppler-broadened 511 plus the 1274.5 line.
  const vector<SyntheticPeak> peaks{
    { 511.0,  30000.0, 2.2*hpge_fwhm(511.0) },  //annihilation: much wider than the trend
    { 1274.5,  9000.0, hpge_fwhm(1274.5) },
  };

  const shared_ptr<Measurement> meas
              = make_synthetic_spectrum( nchan, true_offset, true_gain, peaks );

  GuessOptions options;
  const vector<GuessedCal> results = guess_energy_cal( meas, options );

  BOOST_REQUIRE_MESSAGE( !results.empty(), "No calibration guesses returned" );

  const GuessedCal &best = results.front();
  check_cal_close( best, true_offset, true_gain, nchan, 0.005, 3.0 );

  bool have_511 = false;
  for( const PeakAssignment &assign : best.assignments )
    have_511 = ( have_511 || (assign.line.type == LineSource::Type::Annihilation) );
  BOOST_CHECK_MESSAGE( have_511, "511 keV assignment missing from best guess" );
}//BOOST_AUTO_TEST_CASE( Wide511Strategy )


BOOST_AUTO_TEST_CASE( FluorescenceXrayPattern )
{
  set_data_dir();

  const size_t nchan = 16384;
  const double true_offset = 0.5, true_gain = 0.125;

  // Tungsten K x-rays (e.g. a W collimator) plus Cs-137.
  const vector<SyntheticPeak> peaks{
    { 57.98,   8000.0, hpge_fwhm(57.98) },   //W Ka2
    { 59.32,  15000.0, hpge_fwhm(59.32) },   //W Ka1
    { 67.24,   5000.0, hpge_fwhm(67.24) },   //W Kb1
    { 661.7,  12000.0, hpge_fwhm(661.7) },
  };

  const shared_ptr<Measurement> meas
              = make_synthetic_spectrum( nchan, true_offset, true_gain, peaks );

  GuessOptions options;
  const vector<GuessedCal> results = guess_energy_cal( meas, options );

  BOOST_REQUIRE_MESSAGE( !results.empty(), "No calibration guesses returned" );

  const GuessedCal &best = results.front();
  check_cal_close( best, true_offset, true_gain, nchan, 0.005, 3.0 );
}//BOOST_AUTO_TEST_CASE( FluorescenceXrayPattern )


BOOST_AUTO_TEST_CASE( DegenerateGainRejected )
{
  set_data_dir();

  const size_t nchan = 8192;
  const double true_offset = 3.0, true_gain = 0.35;

  const shared_ptr<Measurement> meas
              = make_synthetic_spectrum( nchan, true_offset, true_gain, cs137_norm_peaks() );

  bool is_high_res = false;
  const vector<CandidatePeak> candidates = find_candidate_peaks_channelspace( meas, is_high_res );
  BOOST_REQUIRE( !candidates.empty() );
  BOOST_CHECK( is_high_res );

  const LineLibrary library = build_line_library( {} );
  BOOST_REQUIRE( !library.lines.empty() );

  // A tiny gain packs the whole line library below every peak; without the accidental-match
  //  density discount it would "explain" everything.  The truth must win.
  const double truth_score = score_hypothesis( true_offset, true_gain, candidates, nchan,
                                               library, is_high_res );
  const double tiny_score = score_hypothesis( 0.0, 150.0/static_cast<double>(nchan), candidates,
                                              nchan, library, is_high_res );

  BOOST_CHECK_MESSAGE( truth_score > tiny_score,
      "Truth score " << truth_score << " did not beat degenerate tiny-gain score " << tiny_score );
  BOOST_CHECK_GT( truth_score, 0.0 );
}//BOOST_AUTO_TEST_CASE( DegenerateGainRejected )


BOOST_AUTO_TEST_CASE( RealUraniumFileRecovery )
{
  set_data_dir();

  string path = SpecUtils::append_path( g_test_file_dir.empty() ? string("test_data")
                                        : g_test_file_dir, "uranium_40pct_selfatten_HPGe_15cm.n42" );
  if( !SpecUtils::is_file( path ) )
    path = SpecUtils::append_path( "../test_data", "uranium_40pct_selfatten_HPGe_15cm.n42" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( path ), "Test file not found: " << path );

  SpecUtils::SpecFile specfile;
  BOOST_REQUIRE_MESSAGE( specfile.load_file( path, SpecUtils::ParserType::Auto ),
                         "Failed to parse " << path );

  shared_ptr<const Measurement> orig;
  for( const shared_ptr<const Measurement> &m : specfile.measurements() )
  {
    if( m && (m->num_gamma_channels() >= 4096) && m->energy_calibration()
        && m->energy_calibration()->valid() )
    {
      orig = m;
      break;
    }
  }
  BOOST_REQUIRE_MESSAGE( orig, "No suitable measurement in test file" );

  const size_t nchan = orig->num_gamma_channels();
  const shared_ptr<const EnergyCalibration> true_cal = orig->energy_calibration();

  // Strip the calibration: pretend the file arrived with none (a default 0-3 MeV guess).
  auto blind = make_shared<Measurement>( *orig );
  auto default_cal = make_shared<EnergyCalibration>();
  default_cal->set_default_polynomial( nchan, { 0.0f, 3000.0f/static_cast<float>(nchan) }, {} );
  blind->set_energy_calibration( default_cal );

  GuessOptions options;
  const vector<GuessedCal> results = guess_energy_cal( blind, options );

  BOOST_REQUIRE_MESSAGE( !results.empty(), "No calibration guesses returned for uranium file" );

  // Compare implied energies at the middle and top of the spectrum against the true cal.
  const GuessedCal &best = results.front();
  for( const double channel : { 0.5*nchan, 0.95*nchan } )
  {
    const double e_true = true_cal->energy_for_channel( channel );
    const double e_guess = best.offset + best.gain*channel;
    BOOST_CHECK_MESSAGE( fabs(e_guess - e_true) < 0.01*e_true + 2.0,
        "channel " << channel << ": guessed " << e_guess << " keV vs true " << e_true
        << " keV (provenance: " << best.provenance << ")" );
  }

  BOOST_CHECK_MESSAGE( has_nuclide( best, "U235" ) || has_nuclide( best, "U238" ),
                       "Uranium not among implied nuclides" );
}//BOOST_AUTO_TEST_CASE( RealUraniumFileRecovery )


BOOST_AUTO_TEST_CASE( FeaturelessSpectrumAndCancel )
{
  set_data_dir();

  // Flat Poisson noise - no peaks to find; must not crash, and must not confidently "find" a
  //  calibration.
  const size_t nchan = 4096;
  std::mt19937 rng( 13579u );
  auto counts = make_shared<vector<float>>( nchan );
  std::poisson_distribution<long> pois( 100.0 );
  for( size_t i = 0; i < nchan; ++i )
    (*counts)[i] = static_cast<float>( pois(rng) );

  auto meas = make_shared<Measurement>();
  meas->set_gamma_counts( counts, 600.0f, 600.0f );
  auto cal = make_shared<EnergyCalibration>();
  cal->set_default_polynomial( nchan, { 0.0f, 3000.0f/static_cast<float>(nchan) }, {} );
  meas->set_energy_calibration( cal );

  GuessOptions options;
  const vector<GuessedCal> results = guess_energy_cal( meas, options );

  for( const GuessedCal &guess : results )
    BOOST_CHECK_LT( guess.explained_frac, 0.95 );

  // Cancellation: a pre-set cancel flag must give an empty result (and not hang/crash).
  const vector<SyntheticPeak> peaks{ { 661.7, 20000.0, hpge_fwhm(661.7) } };
  const shared_ptr<Measurement> real_meas = make_synthetic_spectrum( 8192, 0.0, 0.35, peaks );

  auto cancel = make_shared<atomic<bool>>( true );
  GuessOptions cancel_options;
  cancel_options.cancel = cancel;
  const vector<GuessedCal> cancelled = guess_energy_cal( real_meas, cancel_options );
  BOOST_CHECK( cancelled.empty() );

  // Null / tiny spectra: never throws.
  BOOST_CHECK_NO_THROW( guess_energy_cal( nullptr, GuessOptions{} ) );
  auto tiny = make_shared<Measurement>();
  BOOST_CHECK_NO_THROW( guess_energy_cal( tiny, GuessOptions{} ) );
}//BOOST_AUTO_TEST_CASE( FeaturelessSpectrumAndCancel )


BOOST_AUTO_TEST_CASE( UserNuclideHint )
{
  set_data_dir();

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );

  // A nuclide NOT in the built-in common list: Se-75 (264.7, 136.0, 279.5, 400.7 keV).
  const SandiaDecay::Nuclide * const se75 = db->nuclide( "Se75" );
  BOOST_REQUIRE( se75 );

  const size_t nchan = 8192;
  const double true_offset = 0.0, true_gain = 0.2;

  const vector<SyntheticPeak> peaks{
    { 136.0,   8000.0, hpge_fwhm(136.0) },
    { 264.7,  20000.0, hpge_fwhm(264.7) },
    { 279.5,   9000.0, hpge_fwhm(279.5) },
    { 400.7,   4000.0, hpge_fwhm(400.7) },
  };

  const shared_ptr<Measurement> meas
              = make_synthetic_spectrum( nchan, true_offset, true_gain, peaks );

  GuessOptions options;
  options.user_nuclides.push_back( se75 );
  const vector<GuessedCal> results = guess_energy_cal( meas, options );

  BOOST_REQUIRE_MESSAGE( !results.empty(), "No guesses returned with Se75 hint" );

  const GuessedCal &best = results.front();
  check_cal_close( best, true_offset, true_gain, nchan, 0.005, 3.0 );
  BOOST_CHECK_MESSAGE( has_nuclide( best, "Se75" ), "Se75 not in implied nuclides" );
}//BOOST_AUTO_TEST_CASE( UserNuclideHint )
