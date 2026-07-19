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

#define BOOST_TEST_MODULE EnergyCalGainMatch_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/EnergyCalGainMatch.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace std;
using namespace boost::unit_test;

using SpecUtils::Measurement;
using SpecUtils::EnergyCalibration;

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
  /** Path to the 8-detector (A1-A4, B1-B4) Ba-133 portal example spectrum. */
  string ba133_example_path()
  {
    set_data_dir();
    const string path = SpecUtils::append_path( g_data_dir,
                                        "../example_spectra/ba133_source_640s_20100317.n42" );
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( path ), "Example spectrum not found: " << path );
    return path;
  }//ba133_example_path()


  /** Loads the 8-detector Ba-133 portal example spectrum as a plain SpecFile. */
  shared_ptr<SpecUtils::SpecFile> load_ba133_example()
  {
    auto spec = make_shared<SpecUtils::SpecFile>();
    BOOST_REQUIRE_MESSAGE( spec->load_file( ba133_example_path(), SpecUtils::ParserType::Auto ),
                           "Failed to load ba133 example" );
    return spec;
  }//load_ba133_example()


  /** Loads the Ba-133 portal example as a SpecMeas (needed for the auto-match analysis). */
  shared_ptr<SpecMeas> load_ba133_specmeas()
  {
    auto spec = make_shared<SpecMeas>();
    BOOST_REQUIRE_MESSAGE( spec->load_file( ba133_example_path(), SpecUtils::ParserType::Auto ),
                           "Failed to load ba133 example as SpecMeas" );
    return spec;
  }//load_ba133_specmeas()


  /** Loads a single-record GADRAS-simulated reference spectrum (RADDATA URI .txt). */
  shared_ptr<const Measurement> load_reference_spectrum( const string &det_dir,
                                                         const string &filename )
  {
    set_data_dir();

    const string path = SpecUtils::append_path(
        SpecUtils::append_path(
          SpecUtils::append_path( g_data_dir, "reference_spectra/Common_Field_Nuclides" ),
          det_dir ),
        filename );
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( path ), "Reference spectrum not found: " << path );

    SpecUtils::SpecFile specfile;
    BOOST_REQUIRE_MESSAGE( specfile.load_file( path, SpecUtils::ParserType::Auto ),
                           "Failed to load: " << path );

    shared_ptr<const Measurement> meas = specfile.measurement_at_index( 0 );
    BOOST_REQUIRE( meas && (meas->num_gamma_channels() > 4)
                   && meas->energy_calibration() && meas->energy_calibration()->valid() );
    return meas;
  }//load_reference_spectrum()


  /** First valid gamma energy calibration for the named detector. */
  shared_ptr<const EnergyCalibration> detector_cal( const SpecUtils::SpecFile &spec,
                                                    const string &det )
  {
    for( const shared_ptr<const Measurement> &m : spec.measurements() )
    {
      if( (m->detector_name() == det) && (m->num_gamma_channels() > 4)
          && m->energy_calibration() && m->energy_calibration()->valid() )
        return m->energy_calibration();
    }
    return nullptr;
  }//detector_cal()


  /** One GainMatchCalc input per gamma detector, each summed over all samples. */
  vector<GainMatchCalc::SpectrumInput> make_detector_inputs( const SpecUtils::SpecFile &spec )
  {
    vector<GainMatchCalc::SpectrumInput> inputs;
    const set<int> &samples = spec.sample_numbers();

    for( const string &det : spec.gamma_detector_names() )
    {
      const shared_ptr<const EnergyCalibration> cal = detector_cal( spec, det );
      if( !cal )
        continue;

      GainMatchCalc::SpectrumInput input;
      input.name = det;
      input.spectrum = spec.sum_measurements( samples, {det}, cal );
      if( input.spectrum )
        inputs.push_back( input );
    }

    return inputs;
  }//make_detector_inputs()


  /** Applies E' = gain*E + offset to every gamma measurement of the named detector. */
  void perturb_detector( SpecUtils::SpecFile &spec, const string &det,
                         const double gain, const double offset )
  {
    for( const shared_ptr<const Measurement> &m : spec.measurements() )
    {
      if( (m->detector_name() == det) && (m->num_gamma_channels() > 4)
          && m->energy_calibration() && m->energy_calibration()->valid() )
      {
        const auto newcal = GainMatchCalc::transform_calibration( m->energy_calibration(),
                                                                  gain, offset );
        spec.set_energy_calibration( newcal, m );
      }
    }
  }//perturb_detector()


  /** Adds a quadratic term to the polynomial cal of the named detector, so its energy scale is
   distorted in a way a purely linear (gain/offset) match cannot correct.
   */
  void perturb_detector_quadratic( SpecUtils::SpecFile &spec, const string &det,
                                   const double c2_add )
  {
    for( const shared_ptr<const Measurement> &m : spec.measurements() )
    {
      const shared_ptr<const EnergyCalibration> cal = m->energy_calibration();
      if( (m->detector_name() == det) && (m->num_gamma_channels() > 4) && cal && cal->valid() )
      {
        vector<float> coefs = cal->coefficients();
        if( coefs.size() < 3 )
          coefs.resize( 3, 0.0f );
        coefs[2] += static_cast<float>( c2_add );
        auto newcal = make_shared<EnergyCalibration>();
        newcal->set_polynomial( cal->num_channels(), coefs, cal->deviation_pairs() );
        spec.set_energy_calibration( newcal, m );
      }
    }
  }//perturb_detector_quadratic()


  /** A Poisson-resampled clone of a spectrum, at reduced intensity, with its energy
   calibration perturbed by E' = gain*E + offset.
   */
  shared_ptr<Measurement> poisson_clone( const Measurement &orig,
                                         const double intensity_frac,
                                         const double gain, const double offset,
                                         const string &det_name,
                                         std::mt19937 &rng )
  {
    const shared_ptr<const vector<float>> &counts = orig.gamma_counts();
    BOOST_REQUIRE( counts && !counts->empty() );

    auto newcounts = make_shared<vector<float>>( counts->size(), 0.0f );
    for( size_t i = 0; i < counts->size(); ++i )
    {
      const double mu = intensity_frac * std::max( 0.0f, (*counts)[i] );
      if( mu > 0.0 )
      {
        std::poisson_distribution<unsigned int> dist( mu );
        (*newcounts)[i] = static_cast<float>( dist(rng) );
      }
    }

    auto meas = make_shared<Measurement>( orig );
    meas->set_gamma_counts( newcounts,
                            static_cast<float>( intensity_frac * orig.live_time() ),
                            static_cast<float>( intensity_frac * orig.real_time() ) );
    meas->set_detector_name( det_name );
    meas->set_energy_calibration( GainMatchCalc::transform_calibration(
                                              orig.energy_calibration(), gain, offset ) );
    return meas;
  }//poisson_clone()
}//namespace


/** Matching an unperturbed file, then the same file with known per-detector gain
 perturbations, must recover exactly the injected perturbations - independent of any
 real-world misalignment already present between the detectors.
 */
BOOST_AUTO_TEST_CASE( Ba133EightDetectorGainRecovery )
{
  const vector<double> perturbs = { 1.0, 0.97, 1.03, 1.05, 0.99, 1.5, 0.98, 1.02 };

  GainMatchCalc::MatchOptions options;
  options.lower_energy = 75.0f;  //Ba-133 emissions top out at 384 keV
  options.upper_energy = 0.0f;

  // Baseline: match the file as-is, reference fixed to detector 0
  const shared_ptr<SpecUtils::SpecFile> baseline_file = load_ba133_example();
  const vector<GainMatchCalc::SpectrumInput> baseline_inputs = make_detector_inputs( *baseline_file );
  BOOST_REQUIRE_EQUAL( baseline_inputs.size(), 8 );

  const GainMatchCalc::MatchResults baseline = GainMatchCalc::match( baseline_inputs, 0, options );

  // Perturbed: same file with known gain perturbations applied to detectors 1-7
  const shared_ptr<SpecUtils::SpecFile> perturbed_file = load_ba133_example();
  const vector<string> dets = perturbed_file->gamma_detector_names();
  BOOST_REQUIRE_EQUAL( dets.size(), 8 );
  for( size_t i = 1; i < dets.size(); ++i )
    perturb_detector( *perturbed_file, dets[i], perturbs[i], 0.0 );

  const vector<GainMatchCalc::SpectrumInput> inputs = make_detector_inputs( *perturbed_file );
  BOOST_REQUIRE_EQUAL( inputs.size(), 8 );

  const GainMatchCalc::MatchResults matched = GainMatchCalc::match( inputs, 0, options );
  BOOST_REQUIRE_EQUAL( matched.reference_index, 0 );
  BOOST_REQUIRE_EQUAL( matched.results.size(), 8 );

  // Auto-reference selection must pick a usable spectrum
  const GainMatchCalc::MatchResults auto_ref = GainMatchCalc::match( inputs, -1, options );
  BOOST_CHECK( auto_ref.results.at(auto_ref.reference_index).used );

  for( size_t i = 1; i < 8; ++i )
  {
    const GainMatchCalc::SpectrumResult &res = matched.results[i];
    const GainMatchCalc::SpectrumResult &base = baseline.results[i];

    BOOST_REQUIRE_MESSAGE( res.used, "Detector " << dets[i] << " not matched" );
    BOOST_REQUIRE_MESSAGE( base.used, "Detector " << dets[i] << " not matched in baseline" );
    BOOST_CHECK_MESSAGE( res.correlation > 0.5,
        "Detector " << dets[i] << " correlation only " << res.correlation );

    // Real portal detectors are already roughly matched, so the baseline (unperturbed)
    //  match must be close to identity - guards against degenerate solutions
    BOOST_CHECK_MESSAGE( fabs(base.gain - 1.0) < 0.1,
        "Detector " << dets[i] << " baseline gain " << base.gain << " is far from 1" );

    // Matching the perturbed data must compose with the perturbation to give the
    //  baseline transform: gain_matched * perturb == gain_baseline
    const double recovered = res.gain * perturbs[i];
    BOOST_CHECK_MESSAGE( fabs(recovered/base.gain - 1.0) < 2.0E-3,
        "Detector " << dets[i] << ": perturb " << perturbs[i] << ", matched gain " << res.gain
        << ", baseline gain " << base.gain << " (ratio-1 = " << (recovered/base.gain - 1.0) << ")" );

    BOOST_REQUIRE( res.updated_cal && res.updated_cal->valid() );
  }//for( each non-reference detector )
}//BOOST_AUTO_TEST_CASE( Ba133EightDetectorGainRecovery )


/** Gain plus offset perturbation, with offset fitting enabled. */
BOOST_AUTO_TEST_CASE( Ba133GainAndOffsetRecovery )
{
  const double gain_pert = 1.02, offset_pert = 15.0;

  GainMatchCalc::MatchOptions options;
  options.lower_energy = 75.0f;
  options.fit_offset = true;

  const shared_ptr<SpecUtils::SpecFile> baseline_file = load_ba133_example();
  const vector<GainMatchCalc::SpectrumInput> baseline_inputs = make_detector_inputs( *baseline_file );
  BOOST_REQUIRE( baseline_inputs.size() == 8 );
  const GainMatchCalc::MatchResults baseline = GainMatchCalc::match( baseline_inputs, 0, options );

  const shared_ptr<SpecUtils::SpecFile> perturbed_file = load_ba133_example();
  const vector<string> dets = perturbed_file->gamma_detector_names();
  perturb_detector( *perturbed_file, dets[1], gain_pert, offset_pert );

  const vector<GainMatchCalc::SpectrumInput> inputs = make_detector_inputs( *perturbed_file );
  const GainMatchCalc::MatchResults matched = GainMatchCalc::match( inputs, 0, options );

  const GainMatchCalc::SpectrumResult &res = matched.results.at(1);
  const GainMatchCalc::SpectrumResult &base = baseline.results.at(1);
  BOOST_REQUIRE( res.used && base.used );

  // Gain and offset trade off along a shallow chi2 valley, so compare the composed energy
  //  maps (which is what actually matters) rather than the raw parameters:
  //  gain*(pert_g*E + pert_b) + offset  must equal  base_gain*E + base_offset
  //  to within a fraction of the detector resolution, across the fitted window.
  for( const double energy : { 100.0, 200.0, 400.0, 800.0, 1600.0, 2800.0 } )
  {
    const double matched_energy = res.gain*(gain_pert*energy + offset_pert) + res.offset;
    const double baseline_energy = base.gain*energy + base.offset;
    // The chi2 only pins the energy map where the spectrum has features, so allow a tenth
    //  of the (energy-dependent) NaI FWHM, with a 3 keV floor
    const double tolerance = std::max( 3.0,
                          0.1*PeakFitUtils::nai_fwhm_fcn( static_cast<float>(energy) ) );
    BOOST_CHECK_MESSAGE( fabs(matched_energy - baseline_energy) < tolerance,
        "at " << energy << " keV: matched maps to " << matched_energy
        << ", baseline to " << baseline_energy
        << " (matched g=" << res.gain << " b=" << res.offset
        << "; baseline g=" << base.gain << " b=" << base.offset << ")" );
  }
}//BOOST_AUTO_TEST_CASE( Ba133GainAndOffsetRecovery )


/** Synthetic multi-detector NaI files: Poisson-resampled clones of a single GADRAS-simulated
 spectrum with known gain perturbations must be matched back to the unperturbed clone.
 */
BOOST_AUTO_TEST_CASE( SyntheticNaIGainRecovery )
{
  const vector<double> perturbs = { 1.0, 0.95, 1.07, 1.25 };

  for( const string filename : { "Cs137_Unshielded.txt", "K40_Unshielded.txt" } )
  {
    const shared_ptr<const Measurement> orig = load_reference_spectrum( "2x4x16-NaI", filename );

    std::mt19937 rng( 20260716u );
    vector<GainMatchCalc::SpectrumInput> inputs;
    for( size_t i = 0; i < perturbs.size(); ++i )
    {
      GainMatchCalc::SpectrumInput input;
      input.name = "det" + std::to_string(i);
      input.spectrum = poisson_clone( *orig, 0.25, perturbs[i], 0.0, input.name, rng );
      inputs.push_back( input );
    }

    GainMatchCalc::MatchOptions options;
    options.lower_energy = 100.0f;

    const GainMatchCalc::MatchResults matched = GainMatchCalc::match( inputs, 0, options );

    for( size_t i = 1; i < perturbs.size(); ++i )
    {
      const GainMatchCalc::SpectrumResult &res = matched.results.at(i);
      BOOST_REQUIRE_MESSAGE( res.used, filename << " det" << i << " not matched" );

      const double recovered = res.gain * perturbs[i];
      BOOST_CHECK_MESSAGE( fabs(recovered - 1.0) < 5.0E-3,
          filename << " det" << i << ": perturb " << perturbs[i] << ", matched gain "
          << res.gain << " (product-1 = " << (recovered - 1.0) << ")" );
    }
  }//for( loop over source spectra )
}//BOOST_AUTO_TEST_CASE( SyntheticNaIGainRecovery )


/** Same as SyntheticNaIGainRecovery, but for a HPGe detector, where the sharp peaks should
 pin the recovered gains much tighter.
 */
BOOST_AUTO_TEST_CASE( SyntheticHPGeGainRecovery )
{
  const vector<double> perturbs = { 1.0, 0.99, 1.02, 1.2 };

  const shared_ptr<const Measurement> orig
                          = load_reference_spectrum( "Detective X", "Th232_Unshielded.txt" );

  std::mt19937 rng( 20260716u );
  vector<GainMatchCalc::SpectrumInput> inputs;
  for( size_t i = 0; i < perturbs.size(); ++i )
  {
    GainMatchCalc::SpectrumInput input;
    input.name = "det" + std::to_string(i);
    input.spectrum = poisson_clone( *orig, 0.25, perturbs[i], 0.0, input.name, rng );
    inputs.push_back( input );
  }

  GainMatchCalc::MatchOptions options;
  options.lower_energy = 100.0f;

  const GainMatchCalc::MatchResults matched = GainMatchCalc::match( inputs, 0, options );

  for( size_t i = 1; i < perturbs.size(); ++i )
  {
    const GainMatchCalc::SpectrumResult &res = matched.results.at(i);
    BOOST_REQUIRE_MESSAGE( res.used, "det" << i << " not matched" );

    const double recovered = res.gain * perturbs[i];
    BOOST_CHECK_MESSAGE( fabs(recovered - 1.0) < 5.0E-4,
        "det" << i << ": perturb " << perturbs[i] << ", matched gain " << res.gain
        << " (product-1 = " << (recovered - 1.0) << ")" );
  }
}//BOOST_AUTO_TEST_CASE( SyntheticHPGeGainRecovery )


/** Peak refinement: after matching with a refinement peak, the corrected spectra must have
 that peak within a small fraction of a FWHM of the references fitted mean.
 */
BOOST_AUTO_TEST_CASE( PeakRefinement )
{
  struct RefineCase
  {
    string det_dir, filename;
    double peak_energy;
    float fwhm;
  };

  const RefineCase cases[] = {
    { "2x4x16-NaI", "Cs137_Unshielded.txt", 661.657, PeakFitUtils::nai_fwhm_fcn(661.657f) },
    { "Detective X", "Th232_Unshielded.txt", 2614.511, PeakFitUtils::hpge_fwhm_fcn(2614.511f) },
  };

  for( const RefineCase &test_case : cases )
  {
    const shared_ptr<const Measurement> orig
                     = load_reference_spectrum( test_case.det_dir, test_case.filename );

    const vector<double> perturbs = { 1.0, 0.95, 1.15 };
    std::mt19937 rng( 20260716u );
    vector<GainMatchCalc::SpectrumInput> inputs;
    for( size_t i = 0; i < perturbs.size(); ++i )
    {
      GainMatchCalc::SpectrumInput input;
      input.name = "det" + std::to_string(i);
      input.spectrum = poisson_clone( *orig, 0.25, perturbs[i], 0.0, input.name, rng );
      inputs.push_back( input );
    }

    GainMatchCalc::MatchOptions options;
    options.lower_energy = 100.0f;
    options.ref_peak_energy = test_case.peak_energy;

    const GainMatchCalc::MatchResults matched = GainMatchCalc::match( inputs, 0, options );

    BOOST_REQUIRE_MESSAGE( matched.ref_peak_mean > 0.0,
        test_case.filename << ": reference peak fit failed" );

    for( size_t i = 1; i < perturbs.size(); ++i )
    {
      const GainMatchCalc::SpectrumResult &res = matched.results.at(i);
      BOOST_REQUIRE_MESSAGE( res.used, test_case.filename << " det" << i << " not matched" );
      BOOST_CHECK_MESSAGE( res.fit_peak_mean > 0.0,
          test_case.filename << " det" << i << ": peak refinement did not run" );

      // Re-fit the peak on the corrected spectrum; it must land on the references mean
      auto corrected = make_shared<Measurement>( *inputs[i].spectrum );
      corrected->set_energy_calibration( res.updated_cal );

      double corrected_mean = 0.0;
      BOOST_REQUIRE_NO_THROW( corrected_mean = GainMatchCalc::fit_peak_mean(
                                          corrected, matched.ref_peak_mean, nullptr ) );
      BOOST_CHECK_MESSAGE( fabs(corrected_mean - matched.ref_peak_mean) < 0.1*test_case.fwhm,
          test_case.filename << " det" << i << ": corrected peak at " << corrected_mean
          << " keV vs reference " << matched.ref_peak_mean << " keV (FWHM "
          << test_case.fwhm << ")" );
    }
  }//for( loop over test cases )
}//BOOST_AUTO_TEST_CASE( PeakRefinement )


/** transform_calibration must be an exact linear map for all calibration types, and
 composing with the inverse transform must give back the original energies.
 */
BOOST_AUTO_TEST_CASE( TransformCalibrationRoundTrip )
{
  const double gain = 1.35, offset = 7.5;
  const size_t nchannel = 1024;
  const vector<pair<float,float>> dev_pairs = { {59.5f, -1.2f}, {661.7f, 2.1f}, {2614.5f, 0.0f} };

  auto poly = make_shared<EnergyCalibration>();
  poly->set_polynomial( nchannel, {0.0f, 3.0f}, dev_pairs );

  auto frf = make_shared<EnergyCalibration>();
  frf->set_full_range_fraction( nchannel, {0.0f, 3072.0f, 12.0f}, dev_pairs );

  auto lower_chan = make_shared<EnergyCalibration>();
  lower_chan->set_lower_channel_energy( nchannel, vector<float>(*poly->channel_energies()) );

  for( const shared_ptr<const EnergyCalibration> cal :
                  vector<shared_ptr<const EnergyCalibration>>{ poly, frf, lower_chan } )
  {
    const auto transformed = GainMatchCalc::transform_calibration( cal, gain, offset );
    BOOST_REQUIRE( transformed && transformed->valid() );
    BOOST_REQUIRE_EQUAL( static_cast<int>(transformed->type()), static_cast<int>(cal->type()) );

    // Forward: E' = gain*E + offset (deviation pairs are identical in both, so they cancel
    //  when comparing at the same channel for polynomial/FRF; exact for lower channel edge)
    const auto roundtrip = GainMatchCalc::transform_calibration( transformed,
                                                                 1.0/gain, -offset/gain );

    for( const double channel : { 0.0, 100.5, 512.0, 1000.25, 1024.0 } )
    {
      const double orig_energy = cal->energy_for_channel( channel );
      const double rt_energy = roundtrip->energy_for_channel( channel );
      BOOST_CHECK_MESSAGE( fabs(rt_energy - orig_energy) < 0.01,
          "Round trip at channel " << channel << ": " << rt_energy << " vs " << orig_energy );
    }
  }//for( loop over calibration types )

  // Invalid inputs must throw
  BOOST_CHECK_THROW( GainMatchCalc::transform_calibration( poly, 0.0, 0.0 ), std::exception );
  BOOST_CHECK_THROW( GainMatchCalc::transform_calibration( poly, -1.0, 0.0 ), std::exception );
  BOOST_CHECK_THROW( GainMatchCalc::transform_calibration( nullptr, 1.0, 0.0 ), std::exception );
}//BOOST_AUTO_TEST_CASE( TransformCalibrationRoundTrip )


/** Edge cases: empty spectra excluded, too-few-usable throws, and perturbations beyond the
 gain search window flagged rather than crashing.
 */
BOOST_AUTO_TEST_CASE( EdgeCases )
{
  const shared_ptr<const Measurement> orig
                       = load_reference_spectrum( "2x4x16-NaI", "Cs137_Unshielded.txt" );

  std::mt19937 rng( 20260716u );

  GainMatchCalc::MatchOptions options;
  options.lower_energy = 100.0f;

  // A zero-count spectrum gets excluded with a LowStatistics warning
  {
    auto empty_meas = make_shared<Measurement>( *orig );
    auto zero_counts = make_shared<vector<float>>( orig->gamma_counts()->size(), 0.0f );
    empty_meas->set_gamma_counts( zero_counts, orig->live_time(), orig->real_time() );

    vector<GainMatchCalc::SpectrumInput> inputs( 3 );
    inputs[0].name = "det0";
    inputs[0].spectrum = poisson_clone( *orig, 0.25, 1.0, 0.0, "det0", rng );
    inputs[1].name = "det1";
    inputs[1].spectrum = poisson_clone( *orig, 0.25, 1.03, 0.0, "det1", rng );
    inputs[2].name = "empty";
    inputs[2].spectrum = empty_meas;

    const GainMatchCalc::MatchResults matched = GainMatchCalc::match( inputs, 0, options );
    BOOST_CHECK( !matched.results.at(2).used );
    const vector<GainMatchCalc::WarningType> &warns = matched.results.at(2).warnings;
    BOOST_CHECK( std::count( warns.begin(), warns.end(),
                             GainMatchCalc::WarningType::LowStatistics ) );
    BOOST_CHECK( matched.results.at(1).used );
  }

  // Fewer than two usable spectra throws
  {
    vector<GainMatchCalc::SpectrumInput> inputs( 2 );
    inputs[0].name = "det0";
    inputs[0].spectrum = poisson_clone( *orig, 0.25, 1.0, 0.0, "det0", rng );
    inputs[1].name = "null";

    BOOST_CHECK_THROW( GainMatchCalc::match( inputs, 0, options ), std::exception );
  }

  // Choosing an unusable reference throws
  {
    vector<GainMatchCalc::SpectrumInput> inputs( 3 );
    inputs[0].name = "det0";
    inputs[0].spectrum = poisson_clone( *orig, 0.25, 1.0, 0.0, "det0", rng );
    inputs[1].name = "det1";
    inputs[1].spectrum = poisson_clone( *orig, 0.25, 1.0, 0.0, "det1", rng );
    inputs[2].name = "null";

    BOOST_CHECK_THROW( GainMatchCalc::match( inputs, 2, options ), std::exception );
  }

  // An explicit upper energy limit works
  {
    vector<GainMatchCalc::SpectrumInput> inputs( 2 );
    inputs[0].name = "det0";
    inputs[0].spectrum = poisson_clone( *orig, 0.25, 1.0, 0.0, "det0", rng );
    inputs[1].name = "det1";
    inputs[1].spectrum = poisson_clone( *orig, 0.25, 1.04, 0.0, "det1", rng );

    options.upper_energy = 800.0f;
    const GainMatchCalc::MatchResults matched = GainMatchCalc::match( inputs, 0, options );
    options.upper_energy = 0.0f;

    BOOST_REQUIRE( matched.results.at(1).used );
    BOOST_CHECK_MESSAGE( fabs(matched.results.at(1).gain * 1.04 - 1.0) < 5.0E-3,
        "restricted-range matched gain " << matched.results.at(1).gain );
  }

  // A perturbation beyond the +-8x search window must not crash; it either fails the row
  //  or flags the correlation as being at the window edge
  {
    vector<GainMatchCalc::SpectrumInput> inputs( 2 );
    inputs[0].name = "det0";
    inputs[0].spectrum = poisson_clone( *orig, 0.25, 1.0, 0.0, "det0", rng );
    inputs[1].name = "way_off";
    inputs[1].spectrum = poisson_clone( *orig, 0.25, 12.0, 0.0, "way_off", rng );

    GainMatchCalc::MatchResults matched;
    BOOST_REQUIRE_NO_THROW( matched = GainMatchCalc::match( inputs, 0, options ) );

    const GainMatchCalc::SpectrumResult &res = matched.results.at(1);
    BOOST_TEST_MESSAGE( "12x perturbation: used=" << res.used << " gain=" << res.gain
                        << " correlation=" << res.correlation
                        << " nwarnings=" << res.warnings.size() );
    const bool flagged = !res.used
      || std::count( res.warnings.begin(), res.warnings.end(),
                     GainMatchCalc::WarningType::CorrelationAtEdge )
      || std::count( res.warnings.begin(), res.warnings.end(),
                     GainMatchCalc::WarningType::WeakCorrelation )
      || std::count( res.warnings.begin(), res.warnings.end(),
                     GainMatchCalc::WarningType::MatchFailed );
    BOOST_CHECK_MESSAGE( flagged, "12x perturbation was not flagged in any way" );
  }
}//BOOST_AUTO_TEST_CASE( EdgeCases )


/** The on-load auto-match analysis compares the implied inter-detector shift to the intrinsic
 peak width, so it stays quiet when a low-resolution detector's broad peaks tolerate the
 misalignment, but flags a grossly mis-calibrated detector and returns a correcting gain.
 */
BOOST_AUTO_TEST_CASE( AutoMatchDetection )
{
  // The real Ba-133 portal detectors differ in gain by ~1% (~4 keV near 356 keV), but the NaI
  //  peaks there are ~30 keV wide, so matching would not visibly sharpen them - no suggestion.
  const shared_ptr<SpecMeas> mild = load_ba133_specmeas();
  const set<int> samples = mild->sample_numbers();

  const shared_ptr<GainMatchCalc::DetectorMatchSuggestion> none
      = GainMatchCalc::analyzeForAutoMatch( SpecUtils::SpectrumType::Foreground, mild, samples );
  BOOST_CHECK_MESSAGE( !none, "Auto-match suggested a change for a low-res file within resolution"
                       << (none ? " (shift " + std::to_string(none->max_shift_kev) + " keV)" : "") );

  // Grossly mis-calibrate one detector (+20% gain -> tens of keV at 356 keV, far exceeding the
  //  NaI peak width); now the analysis should offer to fix it.
  const shared_ptr<SpecMeas> perturbed = load_ba133_specmeas();
  const vector<string> dets = perturbed->gamma_detector_names();
  BOOST_REQUIRE( dets.size() == 8 );
  perturb_detector( *perturbed, dets[3], 1.20, 0.0 );

  const shared_ptr<GainMatchCalc::DetectorMatchSuggestion> suggestion
      = GainMatchCalc::analyzeForAutoMatch( SpecUtils::SpectrumType::Foreground, perturbed, samples );
  BOOST_REQUIRE_MESSAGE( suggestion, "Auto-match did not flag a grossly mis-calibrated detector" );
  BOOST_CHECK( suggestion->max_shift_kev > 3.0 );
  BOOST_CHECK( !suggestion->per_detector.empty() );

  // The perturbed detector should be among the suggested corrections
  bool found = false;
  for( const auto &det_gain_offset : suggestion->per_detector )
  {
    if( std::get<0>(det_gain_offset) == dets[3] )
      found = true;
  }
  BOOST_CHECK_MESSAGE( found, "The mis-calibrated detector was not in the suggestion" );
}//BOOST_AUTO_TEST_CASE( AutoMatchDetection )


/** Multi-peak refinement must correct a quadratic distortion that a linear (gain/offset) match
 cannot: after fitting an order-2 polynomial to the shared Ba-133 peaks, the distorted detector's
 peaks should land on the reference detector's peak energies far more tightly than linear-only.
 */
BOOST_AUTO_TEST_CASE( MultiPeakRefinement )
{
  const shared_ptr<SpecMeas> file = load_ba133_specmeas();
  const set<int> samples = file->sample_numbers();
  const vector<string> dets = file->gamma_detector_names();
  BOOST_REQUIRE( dets.size() == 8 );

  // Distort detector index 2 with a quadratic term (a few keV of curvature across the range),
  //  which a linear match structurally cannot remove.
  perturb_detector_quadratic( *file, dets[2], 1.0e-6 );

  const vector<GainMatchCalc::SpectrumInput> inputs = GainMatchCalc::buildDetectorInputs( file, samples );
  BOOST_REQUIRE_EQUAL( inputs.size(), 8 );

  size_t di = inputs.size();
  for( size_t i = 0; i < inputs.size(); ++i )
    if( inputs[i].name == dets[2] )
      di = i;
  BOOST_REQUIRE( di < inputs.size() );

  GainMatchCalc::MatchOptions options;
  options.lower_energy = 75.0f;   //Ba-133 lines span ~81-384 keV

  // Fix the reference to an un-distorted detector so the distorted one is a test spectrum.
  const size_t ref = (di == 0) ? 1 : 0;
  const GainMatchCalc::MatchResults stage2 = GainMatchCalc::match( inputs, static_cast<int>(ref), options );
  BOOST_REQUIRE_EQUAL( stage2.reference_index, ref );

  const deque<shared_ptr<const PeakDef>> no_id;
  const vector<GainMatchCalc::SharedPeak> shared
      = GainMatchCalc::findSharedPeaks( inputs, stage2, no_id, 75.0f, 0.0f );
  BOOST_TEST_MESSAGE( "found " << shared.size() << " shared Ba-133 peaks" );
  BOOST_REQUIRE_MESSAGE( shared.size() >= 3, "need >= 3 shared peaks for a quadratic fit; got "
                         << shared.size() );

  // Residual of the distorted detector BEFORE refinement (linear stage-2 only): its peaks,
  //  mapped through the stage-2 gain/offset, vs the reference's peak energies.
  const shared_ptr<const EnergyCalibration> orig_cal = inputs[di].spectrum->energy_calibration();
  const shared_ptr<const EnergyCalibration> lin_cal
      = GainMatchCalc::transform_calibration( orig_cal, stage2.results[di].gain,
                                              stage2.results[di].offset );
  double linear_resid = 0.0;
  for( const GainMatchCalc::SharedPeak &pk : shared )
  {
    if( pk.detector_channels[di] >= 0.0 )
      linear_resid = std::max( linear_resid,
          std::fabs( lin_cal->energy_for_channel(pk.detector_channels[di]) - pk.target_energy ) );
  }

  // Refine with an order-2 (quadratic) polynomial from all shared peaks.
  const vector<bool> use( shared.size(), true );
  const vector<shared_ptr<const EnergyCalibration>> cals
      = GainMatchCalc::refineWithSharedPeaks( inputs, stage2, shared, use, 2 );
  BOOST_REQUIRE_EQUAL( cals.size(), 8 );
  BOOST_REQUIRE_MESSAGE( cals[di], "distorted detector was not refined" );
  BOOST_CHECK_MESSAGE( !cals[ref], "reference detector should be untouched in relative mode" );

  double refined_resid = 0.0;
  for( const GainMatchCalc::SharedPeak &pk : shared )
  {
    if( pk.detector_channels[di] >= 0.0 )
      refined_resid = std::max( refined_resid,
          std::fabs( cals[di]->energy_for_channel(pk.detector_channels[di]) - pk.target_energy ) );
  }

  BOOST_TEST_MESSAGE( "linear residual " << linear_resid << " keV, refined residual "
                      << refined_resid << " keV" );
  BOOST_CHECK_MESSAGE( refined_resid < 0.5,
      "quadratic refinement residual too large: " << refined_resid << " keV" );
  BOOST_CHECK_MESSAGE( refined_resid < 0.25*linear_resid,
      "refinement did not substantially beat linear (" << refined_resid << " vs " << linear_resid << ")" );
}//BOOST_AUTO_TEST_CASE( MultiPeakRefinement )
