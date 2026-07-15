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
#include <deque>
#include <cmath>
#include <string>
#include <memory>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

#include <Wt/Utils>
#include <Wt/WApplication>
#include <Wt/Test/WTestEnvironment>

#ifdef _WIN32
#include "winsock2.h"
#include "Windows.h"
#endif

#define BOOST_TEST_MODULE ExportSpecFile_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/ExportSpecFile.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace std;
using namespace boost::unit_test;

std::string g_test_file_dir;

// Set the static data directory
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

  // Search for the data directory if not specified
  if( datadir.empty() )
  {
    for( const auto &d : { "data", "../data", "../../data", "../../../data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) )
      {
        datadir = d;
        break;
      }
    }
  }//if( datadir.empty() )

  // Search for test file directory if not specified
  if( g_test_file_dir.empty() )
  {
    const vector<string> candidate_test_dirs{
      "test_data",
      "../test_data",
      "target/testing/test_data",
      "../target/testing/test_data",
      "../../target/testing/test_data",
      "../../../target/testing/test_data"
    };

    for( const auto &d : candidate_test_dirs )
    {
      if( SpecUtils::is_directory(d) )
      {
        g_test_file_dir = d;
        break;
      }
    }
  }//if( g_test_file_dir.empty() )

  if( datadir.empty() )
    throw runtime_error( "Could not find data directory" );

  if( g_test_file_dir.empty() )
    throw runtime_error( "Could not find test file directory" );

  // Initialize the decay database
  DecayDataBaseServer::setDecayXmlFile( SpecUtils::append_path(datadir, "sandia.decay.xml") );

  // Set other data directories
  InterSpec::setStaticDataDirectory( datadir );
}//void set_data_dir()


namespace
{
/** Returns a deque holding a single Gaussian peak at 661.7 keV. */
deque<shared_ptr<const PeakDef>> make_test_peaks()
{
  auto cont = make_shared<PeakContinuum>();
  cont->setType( PeakContinuum::OffsetType::Linear );
  cont->setRange( 630.0, 690.0 );

  auto peak = make_shared<PeakDef>( 661.7, 6.0, 1000.0 );
  peak->setContinuum( cont );

  return deque<shared_ptr<const PeakDef>>{ peak };
}//make_test_peaks()


/** One sample of a synthetic file: its sample number, the detectors present for it,
 and its source type.
 */
struct SynthSample
{
  int sample_number;
  vector<string> detector_names;
  SpecUtils::SourceType source_type;
};//struct SynthSample


/** Creates a SpecMeas from the given sample descriptions: 1024 channels at 3 keV/channel,
 flat continuum plus a Gaussian bump at 661.7 keV whose area scales with sample number and
 detector index (so exported records can be told apart when manually checking in the GUI).
 */
shared_ptr<SpecMeas> create_synthetic_file( const vector<SynthSample> &samples )
{
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( 1024, {0.0f, 3.0f}, {} );

  auto meas = make_shared<SpecMeas>();
  for( const SynthSample &sample : samples )
  {
    for( size_t det_index = 0; det_index < sample.detector_names.size(); ++det_index )
    {
      const double amp = 100.0 * sample.sample_number * (1.0 + det_index);
      auto counts = make_shared<vector<float>>( 1024, 10.0f );
      // Gaussian at 661.7 keV (channel ~220.57), sigma 2 channels
      for( size_t ch = 200; ch < 240; ++ch )
      {
        const double x = (ch + 0.5) - 220.57;
        (*counts)[ch] += static_cast<float>( amp * std::exp( -0.5*x*x/(2.0*2.0) )
                                            / (2.0*std::sqrt(2.0*3.14159265358979)) );
      }

      auto m = make_shared<SpecUtils::Measurement>();
      m->set_gamma_counts( counts, 300.0f, 315.0f );
      m->set_energy_calibration( cal );
      m->set_sample_number( sample.sample_number );
      m->set_detector_name( sample.detector_names[det_index] );
      m->set_title( "Sample " + std::to_string(sample.sample_number) );
      m->set_source_type( sample.source_type );
      meas->add_measurement( m, false );
    }//for( loop over detectors of this sample )
  }//for( const SynthSample &sample : samples )

  meas->cleanup_after_load( SpecUtils::SpecFile::DontChangeOrReorderSamples );

  return meas;
}//create_synthetic_file(...)


/** Synthesizes the N42 test-file variants into `<testfiledir>/ExportSpecFile/`, overwriting
 on each run.  Files are left on disk (not temp files) so they can be manually loaded in the
 InterSpec GUI to reproduce/inspect export behavior.  Returns the directory path.
 */
string ensure_synthetic_files()
{
  static string s_dir;
  if( !s_dir.empty() )
    return s_dir;

  set_data_dir();

  const string dir = SpecUtils::append_path( g_test_file_dir, "ExportSpecFile" );
  if( !SpecUtils::is_directory(dir) )
    BOOST_REQUIRE_EQUAL( SpecUtils::create_directory(dir), 1 );

  const SpecUtils::SourceType fore_type = SpecUtils::SourceType::Foreground;
  const SpecUtils::SourceType back_type = SpecUtils::SourceType::Background;

  const auto write_file = [&dir]( const shared_ptr<SpecMeas> &meas, const string &name ){
    const string path = SpecUtils::append_path( dir, name );
    ofstream out( path.c_str(), ios::binary );
    BOOST_REQUIRE( out.is_open() );
    BOOST_REQUIRE( meas->write_2012_N42( out ) );
  };

  // The reported-bug scenario: 13 samples, one detector, user peaks on sample 8 only
  {
    vector<SynthSample> samples;
    for( int i = 1; i <= 13; ++i )
      samples.push_back( SynthSample{i, {"Aa1"}, fore_type} );
    const shared_ptr<SpecMeas> meas = create_synthetic_file( samples );
    meas->setPeaks( make_test_peaks(), {8} );
    write_file( meas, "13samples_1det_peaks_on_8.n42" );
  }

  // Multiple detectors per sample: 5 samples, each with two detectors, peaks on sample 3
  {
    vector<SynthSample> samples;
    for( int i = 1; i <= 5; ++i )
      samples.push_back( SynthSample{i, {"Aa1", "Ba2"}, fore_type} );
    const shared_ptr<SpecMeas> meas = create_synthetic_file( samples );
    meas->setPeaks( make_test_peaks(), {3} );
    write_file( meas, "5samples_2dets_peaks_on_3.n42" );
  }

  // Non-uniform detectors: "Aa1" on every sample, "Ba2" only on odd samples
  {
    vector<SynthSample> samples;
    for( int i = 1; i <= 6; ++i )
    {
      const vector<string> dets = (i % 2) ? vector<string>{"Aa1", "Ba2"} : vector<string>{"Aa1"};
      samples.push_back( SynthSample{i, dets, fore_type} );
    }
    const shared_ptr<SpecMeas> meas = create_synthetic_file( samples );
    meas->setPeaks( make_test_peaks(), {3} );
    write_file( meas, "6samples_nonuniform_dets_peaks_on_3.n42" );
  }

  // Non-monotonic / gapped sample numbers, peaks on sample 9
  {
    vector<SynthSample> samples;
    for( const int i : {2, 3, 7, 9, 15} )
      samples.push_back( SynthSample{i, {"Aa1"}, fore_type} );
    const shared_ptr<SpecMeas> meas = create_synthetic_file( samples );
    meas->setPeaks( make_test_peaks(), {9} );
    write_file( meas, "gapped_samples_peaks_on_9.n42" );
  }

  // Foreground samples plus a background sample, peaks on the foreground samples
  {
    vector<SynthSample> samples;
    samples.push_back( SynthSample{1, {"Aa1"}, fore_type} );
    samples.push_back( SynthSample{2, {"Aa1"}, fore_type} );
    samples.push_back( SynthSample{3, {"Aa1"}, back_type} );
    const shared_ptr<SpecMeas> meas = create_synthetic_file( samples );
    meas->setPeaks( make_test_peaks(), {1, 2} );
    write_file( meas, "fore_plus_back_peaks_on_1_2.n42" );
  }

  s_dir = dir;
  return s_dir;
}//ensure_synthetic_files()
}//namespace


// Helper class to manage InterSpec instance for tests
class InterSpecTestFixture
{
public:
  std::unique_ptr<Wt::Test::WTestEnvironment> m_env;
  InterSpecApp *m_app;
  std::unique_ptr<Wt::WApplication::UpdateLock> m_update_lock;
  InterSpec *m_interspec;

  InterSpecTestFixture()
  : m_env( nullptr ),
    m_app( nullptr ),
    m_interspec( nullptr )
  {
    set_data_dir();
    ensure_synthetic_files();

    string wt_app_root = SpecUtils::append_path( InterSpec::staticDataDirectory(), "..");
    wt_app_root = SpecUtils::lexically_normalize_path(wt_app_root);

    // Create a test environment
    const std::string applicationPath = "";
    const std::string configurationFile = "";
    m_env.reset( new Wt::Test::WTestEnvironment( applicationPath, configurationFile, Wt::Application ) );
    m_env->setAppRoot( wt_app_root );

    // Create the app
    m_app = new InterSpecApp( *m_env );

    m_update_lock.reset( new Wt::WApplication::UpdateLock(m_app) );

    // Get the InterSpec viewer instance
    m_interspec = m_app->viewer();
    BOOST_REQUIRE( m_interspec );
  }

  /** Loads one of the synthesized N42 variants as the foreground. */
  shared_ptr<SpecMeas> loadForeground( const string &filename )
  {
    const string path = SpecUtils::append_path( ensure_synthetic_files(), filename );
    BOOST_REQUIRE( SpecUtils::is_file(path) );

    m_interspec->userOpenFileFromFilesystem( path, path );
    Wt::WApplication::instance()->processEvents();

    shared_ptr<SpecMeas> fore = m_interspec->measurment( SpecUtils::SpectrumType::Foreground );
    BOOST_REQUIRE( fore );

    return fore;
  }//loadForeground(...)

  /** Returns the (single) sample-number set the loaded file has peaks for.
   Using the loaded key - rather than the key the file was written with - keeps the tests
   robust to any sample-number normalization the N42 load path performs.
   */
  set<int> loadedPeakKey( const shared_ptr<SpecMeas> &meas )
  {
    const set<set<int>> keys = meas->sampleNumsWithPeaks();
    BOOST_REQUIRE_EQUAL( static_cast<int>(keys.size()), 1 );
    return *begin(keys);
  }//loadedPeakKey(...)

  ~InterSpecTestFixture()
  {
    m_update_lock.reset();
    m_env.reset();
    m_interspec = nullptr;
  }
};//class InterSpecTestFixture


/** SpecMeas-level guards for the peak-map access semantics the export code relies on. */
BOOST_AUTO_TEST_CASE( specMeasPeakMapGuards )
{
  set_data_dir();

  vector<SynthSample> sample_defs;
  for( int i = 1; i <= 13; ++i )
    sample_defs.push_back( SynthSample{i, {"Aa1"}, SpecUtils::SourceType::Foreground} );

  shared_ptr<SpecMeas> meas = create_synthetic_file( sample_defs );
  meas->setPeaks( make_test_peaks(), {8} );
  const SpecMeas &const_meas = *meas;

  // The const `peaks(...)` overload must not insert an entry on a lookup miss
  BOOST_CHECK( !const_meas.peaks( set<int>{1,2,3} ) );
  BOOST_CHECK( !const_meas.peaks( set<int>{1,2,3} ) );
  BOOST_CHECK_EQUAL( static_cast<int>(meas->sampleNumsWithPeaks().size()), 1 );

  // `removePeaks(...)` must not touch a previously obtained deque - only the map entry
  const shared_ptr<deque<shared_ptr<const PeakDef>>> live_deque = meas->peaks( set<int>{8} );
  BOOST_REQUIRE( live_deque );
  BOOST_CHECK_EQUAL( static_cast<int>(live_deque->size()), 1 );
  meas->removePeaks( set<int>{8} );
  BOOST_CHECK_EQUAL( static_cast<int>(live_deque->size()), 1 );
  BOOST_CHECK( !const_meas.peaks( set<int>{8} ) );

  // `setPeaks(...)` over an existing entry clears/re-fills the stored deque IN-PLACE;
  //  PeakModel relies on this aliasing, so pin the behavior (see notes in SpecMeas.h).
  meas->setPeaks( *live_deque, {8} );
  const shared_ptr<deque<shared_ptr<const PeakDef>>> alias = meas->peaks( set<int>{8} );
  BOOST_REQUIRE( alias );
  BOOST_CHECK_EQUAL( static_cast<int>(alias->size()), 1 );
  meas->setPeaks( {}, {8} );
  BOOST_CHECK( alias->empty() );

  // `change_sample_numbers(...)` must remap the peak-map keys along with the measurements
  meas->removePeaks( set<int>{8} );
  meas->setPeaks( make_test_peaks(), {8} );
  vector<pair<int,int>> renumber;
  for( int sample = 1; sample <= 13; ++sample )
    renumber.emplace_back( sample, sample + 100 );
  meas->change_sample_numbers( renumber );
  BOOST_CHECK( !const_meas.peaks( set<int>{8} ) );
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> remapped = const_meas.peaks( set<int>{108} );
  BOOST_REQUIRE( remapped );
  BOOST_CHECK_EQUAL( static_cast<int>(remapped->size()), 1 );

  // `write_iaea_spe` includes the peak CSV section when (and only when) the sample
  //  numbers passed exactly match the peak-map key
  {
    stringstream matching;
    BOOST_REQUIRE( meas->write_iaea_spe( matching, set<int>{108}, set<int>{} ) );
    BOOST_CHECK( matching.str().find( "$PEAK_INFO_CSV:" ) != string::npos );

    stringstream non_matching;
    BOOST_REQUIRE( meas->write_iaea_spe( non_matching, meas->sample_numbers(), set<int>{} ) );
    BOOST_CHECK( non_matching.str().find( "$PEAK_INFO_CSV:" ) == string::npos );
  }
}//BOOST_AUTO_TEST_CASE( specMeasPeakMapGuards )


/** The reported bug: 13-sample file, displayed foreground = sample 8 only, peaks fit on it,
 exported to IAEA SPE.  The output must be a faithful copy of the sample-8 record (not a
 sum artifact), renumbered to sample 1, with the peaks preserved and written to the SPE file.
 */
BOOST_FIXTURE_TEST_CASE( exportDisplayedForegroundToSpeKeepsPeaks, InterSpecTestFixture )
{
  const shared_ptr<SpecMeas> fore = loadForeground( "13samples_1det_peaks_on_8.n42" );
  BOOST_REQUIRE_EQUAL( static_cast<int>(fore->sample_numbers().size()), 13 );

  const set<int> peak_key = loadedPeakKey( fore );
  BOOST_REQUIRE_EQUAL( static_cast<int>(peak_key.size()), 1 );

  m_interspec->changeDisplayedSampleNums( peak_key, SpecUtils::SpectrumType::Foreground );
  Wt::WApplication::instance()->processEvents();

  // The tool must be attached to the widget tree - `generateFileToSave()` checks
  //  `isVisible()` on the option checkboxes, and orphan widgets report not-visible.
  ExportSpecFileTool *tool = new ExportSpecFileTool( m_interspec, m_app->root() );
  tool->handleAppUrl( "V=1&FORE=1&FORMAT=IAEA-SPE&DISPFORE=1" );

  const shared_ptr<const SpecMeas> result = tool->generateFileToSave();
  BOOST_REQUIRE( result );
  BOOST_REQUIRE_EQUAL( static_cast<int>(result->num_measurements()), 1 );

  const string det_name = fore->detector_names().empty() ? string() : fore->detector_names()[0];
  const shared_ptr<const SpecUtils::Measurement> orig = fore->measurement( *begin(peak_key), det_name );
  const shared_ptr<const SpecUtils::Measurement> rec = result->measurements()[0];
  BOOST_REQUIRE( orig );
  BOOST_REQUIRE( rec );

  // Output sample numbers should start at 1
  BOOST_CHECK_EQUAL( rec->sample_number(), 1 );

  // The record must be a faithful copy of the displayed sample - not a sum artifact
  BOOST_CHECK_EQUAL( rec->title(), orig->title() );
  BOOST_CHECK_EQUAL( rec->live_time(), orig->live_time() );
  BOOST_CHECK_EQUAL( rec->real_time(), orig->real_time() );
  BOOST_CHECK( rec->source_type() == SpecUtils::SourceType::Foreground );
  BOOST_REQUIRE( rec->gamma_counts() );
  BOOST_REQUIRE( orig->gamma_counts() );
  BOOST_CHECK( (*rec->gamma_counts()) == (*orig->gamma_counts()) );

  // The user peaks must have survived, re-keyed to sample 1
  const SpecMeas &const_result = *result;
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> out_peaks = const_result.peaks( set<int>{1} );
  BOOST_REQUIRE( out_peaks );
  BOOST_CHECK_EQUAL( static_cast<int>(out_peaks->size()), 1 );

  // And they must make it into the IAEA SPE output
  stringstream spe;
  BOOST_REQUIRE( result->write_iaea_spe( spe, result->sample_numbers(), set<int>{} ) );
  BOOST_CHECK( spe.str().find( "$PEAK_INFO_CSV:" ) != string::npos );
}//BOOST_FIXTURE_TEST_CASE( exportDisplayedForegroundToSpeKeepsPeaks, ... )


/** Summing all samples to a single record: peaks keyed on the full sample set must survive
 (re-keyed to sample 1); peaks keyed on a subset (e.g. just sample 8) are - by design - not
 carried over to the sum (exact-match peak-key semantics).
 */
BOOST_FIXTURE_TEST_CASE( exportAllSamplesSummedPeakSemantics, InterSpecTestFixture )
{
  const shared_ptr<SpecMeas> fore = loadForeground( "13samples_1det_peaks_on_8.n42" );
  const int nsamples = static_cast<int>( fore->sample_numbers().size() );
  BOOST_REQUIRE_EQUAL( nsamples, 13 );

  ExportSpecFileTool *tool = new ExportSpecFileTool( m_interspec, m_app->root() );

  // Case 1: peaks keyed on a single-sample subset, exporting all samples summed:
  //  exact-match semantics mean the summed record carries no peaks.
  tool->handleAppUrl( "V=1&FORE=1&FORMAT=IAEA-SPE&ALLSAMPLES=1" );
  shared_ptr<const SpecMeas> result = tool->generateFileToSave();
  BOOST_REQUIRE( result );
  BOOST_REQUIRE_EQUAL( static_cast<int>(result->num_measurements()), 1 );
  BOOST_CHECK_EQUAL( result->measurements()[0]->sample_number(), 1 );

  {
    const SpecMeas &const_result = *result;
    const shared_ptr<const deque<shared_ptr<const PeakDef>>> out_peaks = const_result.peaks( set<int>{1} );
    BOOST_CHECK( !out_peaks || out_peaks->empty() );

    stringstream spe;
    BOOST_REQUIRE( result->write_iaea_spe( spe, result->sample_numbers(), set<int>{} ) );
    BOOST_CHECK( spe.str().find( "$PEAK_INFO_CSV:" ) == string::npos );
  }

  // Case 2: peaks keyed on the full sample set must survive the sum, re-keyed to sample 1.
  fore->setPeaks( make_test_peaks(), fore->sample_numbers() );

  tool->handleAppUrl( "V=1&FORE=1&FORMAT=IAEA-SPE&ALLSAMPLES=1" );
  result = tool->generateFileToSave();
  BOOST_REQUIRE( result );
  BOOST_REQUIRE_EQUAL( static_cast<int>(result->num_measurements()), 1 );

  {
    const SpecMeas &const_result = *result;
    const shared_ptr<const deque<shared_ptr<const PeakDef>>> out_peaks = const_result.peaks( set<int>{1} );
    BOOST_REQUIRE( out_peaks );
    BOOST_CHECK_EQUAL( static_cast<int>(out_peaks->size()), 1 );

    // The summed record should cover all 13 samples worth of live time
    const shared_ptr<const SpecUtils::Measurement> rec = result->measurements()[0];
    BOOST_CHECK_CLOSE( rec->live_time(), nsamples*300.0f, 1.0e-3 );

    stringstream spe;
    BOOST_REQUIRE( result->write_iaea_spe( spe, result->sample_numbers(), set<int>{} ) );
    BOOST_CHECK( spe.str().find( "$PEAK_INFO_CSV:" ) != string::npos );
  }
}//BOOST_FIXTURE_TEST_CASE( exportAllSamplesSummedPeakSemantics, ... )


/** Multiple detectors per sample: exporting the displayed (single) sample to SPE sums the
 two detectors into one record; the peaks for that sample must survive.
 */
BOOST_FIXTURE_TEST_CASE( exportMultiDetectorSampleToSpeKeepsPeaks, InterSpecTestFixture )
{
  const shared_ptr<SpecMeas> fore = loadForeground( "5samples_2dets_peaks_on_3.n42" );
  BOOST_REQUIRE_EQUAL( static_cast<int>(fore->detector_names().size()), 2 );

  const set<int> peak_key = loadedPeakKey( fore );
  BOOST_REQUIRE_EQUAL( static_cast<int>(peak_key.size()), 1 );

  m_interspec->changeDisplayedSampleNums( peak_key, SpecUtils::SpectrumType::Foreground );
  Wt::WApplication::instance()->processEvents();

  ExportSpecFileTool *tool = new ExportSpecFileTool( m_interspec, m_app->root() );
  tool->handleAppUrl( "V=1&FORE=1&FORMAT=IAEA-SPE&DISPFORE=1" );

  const shared_ptr<const SpecMeas> result = tool->generateFileToSave();
  BOOST_REQUIRE( result );
  BOOST_REQUIRE_EQUAL( static_cast<int>(result->num_measurements()), 1 );

  const shared_ptr<const SpecUtils::Measurement> rec = result->measurements()[0];
  BOOST_REQUIRE( rec );
  BOOST_CHECK_EQUAL( rec->sample_number(), 1 );

  // Both detectors were summed for the displayed sample
  BOOST_CHECK_CLOSE( rec->live_time(), 2.0f*300.0f, 1.0e-3 );

  const SpecMeas &const_result = *result;
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> out_peaks = const_result.peaks( set<int>{1} );
  BOOST_REQUIRE( out_peaks );
  BOOST_CHECK_EQUAL( static_cast<int>(out_peaks->size()), 1 );

  stringstream spe;
  BOOST_REQUIRE( result->write_iaea_spe( spe, result->sample_numbers(), set<int>{} ) );
  BOOST_CHECK( spe.str().find( "$PEAK_INFO_CSV:" ) != string::npos );
}//BOOST_FIXTURE_TEST_CASE( exportMultiDetectorSampleToSpeKeepsPeaks, ... )


/** Non-uniform detectors per sample: peaks fit over displayed samples where one detector is
 missing from some samples must still survive an SPE export of the displayed samples.
 */
BOOST_FIXTURE_TEST_CASE( exportNonUniformDetectorsToSpeKeepsPeaks, InterSpecTestFixture )
{
  const shared_ptr<SpecMeas> fore = loadForeground( "6samples_nonuniform_dets_peaks_on_3.n42" );
  BOOST_REQUIRE_EQUAL( static_cast<int>(fore->detector_names().size()), 2 );

  // Display samples 1 and 2: sample 1 has both detectors, sample 2 only "Aa1"
  const set<int> disp_samples{1, 2};
  m_interspec->changeDisplayedSampleNums( disp_samples, SpecUtils::SpectrumType::Foreground );
  Wt::WApplication::instance()->processEvents();

  fore->setPeaks( make_test_peaks(), disp_samples );

  ExportSpecFileTool *tool = new ExportSpecFileTool( m_interspec, m_app->root() );
  tool->handleAppUrl( "V=1&FORE=1&FORMAT=IAEA-SPE&DISPFORE=1" );

  const shared_ptr<const SpecMeas> result = tool->generateFileToSave();
  BOOST_REQUIRE( result );
  BOOST_REQUIRE_EQUAL( static_cast<int>(result->num_measurements()), 1 );

  const shared_ptr<const SpecUtils::Measurement> rec = result->measurements()[0];
  BOOST_REQUIRE( rec );
  BOOST_CHECK_EQUAL( rec->sample_number(), 1 );

  // Three records summed: sample 1 (2 dets) + sample 2 (1 det)
  BOOST_CHECK_CLOSE( rec->live_time(), 3.0f*300.0f, 1.0e-3 );

  const SpecMeas &const_result = *result;
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> out_peaks = const_result.peaks( set<int>{1} );
  BOOST_REQUIRE( out_peaks );
  BOOST_CHECK_EQUAL( static_cast<int>(out_peaks->size()), 1 );

  stringstream spe;
  BOOST_REQUIRE( result->write_iaea_spe( spe, result->sample_numbers(), set<int>{} ) );
  BOOST_CHECK( spe.str().find( "$PEAK_INFO_CSV:" ) != string::npos );
}//BOOST_FIXTURE_TEST_CASE( exportNonUniformDetectorsToSpeKeepsPeaks, ... )


/** Non-monotonic / gapped sample numbers (2, 3, 7, 9, 15): the SPE export of the displayed
 sample must produce a faithful copy renumbered to sample 1; a whole-file N42 export must
 renumber samples to 1..N with the peak-map key remapped to match.
 */
BOOST_FIXTURE_TEST_CASE( exportGappedSampleNumbers, InterSpecTestFixture )
{
  const shared_ptr<SpecMeas> fore = loadForeground( "gapped_samples_peaks_on_9.n42" );
  BOOST_REQUIRE_EQUAL( static_cast<int>(fore->sample_numbers().size()), 5 );

  // InterSpec-written N42 files must keep their (gapped) sample numbers on load,
  //  so the peak-map key still refers to valid samples
  BOOST_CHECK( fore->sample_numbers().count(9) );

  const set<int> peak_key = loadedPeakKey( fore );
  BOOST_REQUIRE_EQUAL( static_cast<int>(peak_key.size()), 1 );
  const int peak_sample = *begin(peak_key);

  // Zero-based index of the peak sample within the file's ordered sample numbers
  const set<int> &all_samples = fore->sample_numbers();
  const int peak_sample_index = static_cast<int>( std::distance( begin(all_samples),
                                                            all_samples.find(peak_sample) ) );

  m_interspec->changeDisplayedSampleNums( peak_key, SpecUtils::SpectrumType::Foreground );
  Wt::WApplication::instance()->processEvents();

  ExportSpecFileTool *tool = new ExportSpecFileTool( m_interspec, m_app->root() );

  // SPE of the displayed sample: faithful copy, renumbered to 1, peaks intact
  {
    tool->handleAppUrl( "V=1&FORE=1&FORMAT=IAEA-SPE&DISPFORE=1" );
    const shared_ptr<const SpecMeas> result = tool->generateFileToSave();
    BOOST_REQUIRE( result );
    BOOST_REQUIRE_EQUAL( static_cast<int>(result->num_measurements()), 1 );

    const shared_ptr<const SpecUtils::Measurement> rec = result->measurements()[0];
    BOOST_CHECK_EQUAL( rec->sample_number(), 1 );

    const string det_name = fore->detector_names().empty() ? string() : fore->detector_names()[0];
    const shared_ptr<const SpecUtils::Measurement> orig = fore->measurement( peak_sample, det_name );
    BOOST_REQUIRE( orig );
    BOOST_REQUIRE( rec->gamma_counts() && orig->gamma_counts() );
    BOOST_CHECK( (*rec->gamma_counts()) == (*orig->gamma_counts()) );

    const SpecMeas &const_result = *result;
    const shared_ptr<const deque<shared_ptr<const PeakDef>>> out_peaks = const_result.peaks( set<int>{1} );
    BOOST_REQUIRE( out_peaks );
    BOOST_CHECK_EQUAL( static_cast<int>(out_peaks->size()), 1 );

    stringstream spe;
    BOOST_REQUIRE( result->write_iaea_spe( spe, result->sample_numbers(), set<int>{} ) );
    BOOST_CHECK( spe.str().find( "$PEAK_INFO_CSV:" ) != string::npos );
  }

  // Whole file to N42: all samples renumbered 1..N, peak key remapped along with its sample
  {
    tool->handleAppUrl( "V=1&FORE=1&FORMAT=N42-2012&ALLSAMPLES=1" );
    const shared_ptr<const SpecMeas> result = tool->generateFileToSave();
    BOOST_REQUIRE( result );
    BOOST_REQUIRE_EQUAL( static_cast<int>(result->num_measurements()), 5 );

    set<int> expected_samples;
    for( int i = 1; i <= 5; ++i )
      expected_samples.insert( i );
    BOOST_CHECK( result->sample_numbers() == expected_samples );

    const SpecMeas &const_result = *result;
    const int expected_peak_sample = peak_sample_index + 1;
    const shared_ptr<const deque<shared_ptr<const PeakDef>>> out_peaks
                                    = const_result.peaks( set<int>{expected_peak_sample} );
    BOOST_REQUIRE( out_peaks );
    BOOST_CHECK_EQUAL( static_cast<int>(out_peaks->size()), 1 );
  }
}//BOOST_FIXTURE_TEST_CASE( exportGappedSampleNumbers, ... )
