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
#include <deque>
#include <cmath>
#include <chrono>
#include <string>
#include <memory>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

#define BOOST_TEST_MODULE SpecMeasN42RoundTrip_suite
#include <boost/test/included/unit_test.hpp>

#include "rapidxml/rapidxml.hpp"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/RapidXmlUtils.hpp"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/DecayDataBaseServer.h"

#if( USE_LLM_INTERFACE )
#include "InterSpec/LlmConversationHistory.h"
#endif

/* Checks that the sample numbers of a file InterSpec writes to N42-2012 survive being read back in,
 along with everything InterSpec keys by sample number (peaks, auto-search peaks, displayed samples,
 LLM conversations).  Search-mode / portal files whose first sample is a long background are the
 interesting case: SpecUtils writes their <RadMeasurement> ids positionally ("Background",
 "Survey<i>"), which used to shift every sample number down by one on reload.
 */

using namespace std;
using namespace boost::unit_test;

std::string g_test_file_dir;
std::string g_data_dir;


void set_data_dir()
{
  static bool s_have_set = false;
  if( s_have_set )
    return;

  s_have_set = true;

  const int argc = boost::unit_test::framework::master_test_suite().argc;
  char ** const argv = boost::unit_test::framework::master_test_suite().argv;

  for( int i = 1; i < argc; ++i )
  {
    const string arg = argv[i];
    if( SpecUtils::istarts_with( arg, "--datadir=" ) )
      g_data_dir = arg.substr( 10 );

    if( SpecUtils::istarts_with( arg, "--testfiledir=" ) )
      g_test_file_dir = arg.substr( 14 );
  }

  SpecUtils::ireplace_all( g_data_dir, "%20", " " );
  SpecUtils::ireplace_all( g_test_file_dir, "%20", " " );

  if( g_data_dir.empty() )
  {
    for( const char *d : { "data", "../data", "../../data", "../../../data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) )
      {
        g_data_dir = d;
        break;
      }
    }
  }//if( g_data_dir.empty() )

  if( g_test_file_dir.empty() )
  {
    for( const char *d : { "test_data", "../test_data", "target/testing/test_data",
                           "../target/testing/test_data", "../../target/testing/test_data" } )
    {
      if( SpecUtils::is_directory(d) )
      {
        g_test_file_dir = d;
        break;
      }
    }
  }//if( g_test_file_dir.empty() )

  if( g_data_dir.empty() )
    throw runtime_error( "Could not find data directory" );

  if( g_test_file_dir.empty() )
    throw runtime_error( "Could not find test file directory" );

  DecayDataBaseServer::setDecayXmlFile( SpecUtils::append_path(g_data_dir, "sandia.decay.xml") );
  InterSpec::setStaticDataDirectory( g_data_dir );
}//void set_data_dir()


namespace
{
const SpecUtils::SourceType sm_fore = SpecUtils::SourceType::Foreground;
const SpecUtils::SourceType sm_back = SpecUtils::SourceType::Background;


/** One sample of a synthetic file. */
struct SynthSample
{
  int sample_number;
  SpecUtils::SourceType source_type;
  float real_time;
};//struct SynthSample


/** Search-mode layout: a 300 s background, then `num_fore` 4 s foreground samples, numbered by
 `sample_numbers` (background first).
 */
vector<SynthSample> search_mode_layout( const vector<int> &sample_numbers, const float fore_real_time = 4.0f )
{
  vector<SynthSample> samples;
  for( size_t i = 0; i < sample_numbers.size(); ++i )
    samples.push_back( SynthSample{ sample_numbers[i], (i ? sm_fore : sm_back), (i ? fore_real_time : 300.0f) } );
  return samples;
}


vector<int> range_inclusive( const int first, const int last )
{
  vector<int> answer;
  for( int i = first; i <= last; ++i )
    answer.push_back( i );
  return answer;
}


/** Creates a SpecMeas from the sample descriptions: 512 channels at 3 keV/channel, a flat
 continuum whose level is unique to each sample (so each sample can be told apart by its gamma
 sum), and a Cs-137-like bump at 661.7 keV on the sixth record only.  The second detector, if
 any, also has neutron counts, so neutron <GrossCounts> ids get written.
 */
shared_ptr<SpecMeas> create_file( const vector<SynthSample> &samples, const size_t num_dets = 1 )
{
  const size_t nchannel = 512;
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( nchannel, {0.0f, 3.0f}, {} );

  const SpecUtils::time_point_t start = chrono::time_point_cast<chrono::microseconds>(
                                         chrono::system_clock::from_time_t( 1760000000 ) );
  double elapsed = 0.0;

  auto meas = make_shared<SpecMeas>();
  for( size_t index = 0; index < samples.size(); ++index )
  {
    const SynthSample &sample = samples[index];
    for( size_t det = 0; det < num_dets; ++det )
    {
      auto counts = make_shared<vector<float>>( nchannel, static_cast<float>(1 + index + 20*det) );
      if( index == 5 )
      {
        for( size_t ch = 210; ch < 232; ++ch )
        {
          const double x = (ch + 0.5) - 220.57;
          (*counts)[ch] += static_cast<float>( std::round( 2000.0 * std::exp( -0.5*x*x/4.0 ) ) );
        }
      }//if( the "Cs-137" record )

      auto m = make_shared<SpecUtils::Measurement>();
      m->set_gamma_counts( counts, 0.98f*sample.real_time, sample.real_time );
      m->set_energy_calibration( cal );
      m->set_sample_number( sample.sample_number );
      m->set_detector_name( "Aa" + std::to_string(det + 1) );
      m->set_source_type( sample.source_type );
      m->set_start_time( start + chrono::microseconds( static_cast<int64_t>(1.0E6*elapsed) ) );
      if( det == 1 )
        m->set_neutron_counts( { static_cast<float>(3 + index) }, 0.98f*sample.real_time );
      meas->add_measurement( m, false );
    }//for( loop over detectors )

    elapsed += sample.real_time;
  }//for( loop over samples )

  meas->cleanup_after_load( SpecUtils::SpecFile::DontChangeOrReorderSamples );

  return meas;
}//create_file(...)


deque<shared_ptr<const PeakDef>> make_peaks( const double mean, const double amplitude )
{
  auto cont = make_shared<PeakContinuum>();
  cont->setType( PeakContinuum::OffsetType::Linear );
  cont->setRange( mean - 30.0, mean + 30.0 );

  auto peak = make_shared<PeakDef>( mean, 6.0, amplitude );
  peak->setContinuum( cont );

  return deque<shared_ptr<const PeakDef>>{ peak };
}//make_peaks(...)


/** Attaches state to `meas`: user peaks on `single` and on `summed`, auto-search peaks on
 `single`, `summed` as the displayed foreground, and an LLM conversation on `single`.  The peak
 amplitudes identify which set they were put on.
 */
void add_state( const shared_ptr<SpecMeas> &meas, const set<int> &single, const set<int> &summed )
{
  meas->setPeaks( make_peaks( 661.7, 1000.0 ), single );
  meas->setPeaks( make_peaks( 661.7, 2000.0 ), summed );
  meas->setAutomatedSearchPeaks( single, make_shared<deque<shared_ptr<const PeakDef>>>( make_peaks( 661.7, 3000.0 ) ) );
  meas->displayedSpectrumChangedCallback( SpecUtils::SpectrumType::Foreground, meas, summed, meas->detector_names() );

#if( USE_LLM_INTERFACE )
  LlmConversationHistory history;
  history.addUserMessageToMainConversation( "What is in sample " + SpecUtils::sequencesToBriefString(single) + "?" );
  meas->setLlmConversationHistory( single,
                     make_shared<vector<shared_ptr<LlmInteraction>>>( history.getConversations() ) );
#endif
}//add_state(...)


/** Sample number -> (gamma sum, neutron sum, source type) - what identifies each time slice. */
struct SampleFingerprint
{
  double gamma_sum = 0.0;
  double neutron_sum = 0.0;
  int source_type = -1;

  bool operator==( const SampleFingerprint &rhs ) const
  {
    return (fabs(gamma_sum - rhs.gamma_sum) < 1.0E-6*std::max(1.0, gamma_sum))
           && (fabs(neutron_sum - rhs.neutron_sum) < 1.0E-6*std::max(1.0, neutron_sum))
           && (source_type == rhs.source_type);
  }
};//struct SampleFingerprint


std::ostream &operator<<( std::ostream &os, const SampleFingerprint &f )
{
  return os << "{gamma=" << f.gamma_sum << ", neut=" << f.neutron_sum << ", type=" << f.source_type << "}";
}


map<int,SampleFingerprint> fingerprints( const SpecUtils::SpecFile &spec )
{
  map<int,SampleFingerprint> answer;
  for( const shared_ptr<const SpecUtils::Measurement> &m : spec.measurements() )
  {
    SampleFingerprint &f = answer[m->sample_number()];
    f.gamma_sum += m->gamma_count_sum();
    f.neutron_sum += m->neutron_counts_sum();
    f.source_type = std::max( f.source_type, static_cast<int>(m->source_type()) );
  }
  return answer;
}//fingerprints(...)


void check_same_samples( const SpecUtils::SpecFile &orig, const SpecUtils::SpecFile &reloaded, const string &label )
{
  BOOST_CHECK_MESSAGE( orig.sample_numbers() == reloaded.sample_numbers(),
                       label << ": sample numbers changed from {"
                       << SpecUtils::sequencesToBriefString(orig.sample_numbers()) << "} to {"
                       << SpecUtils::sequencesToBriefString(reloaded.sample_numbers()) << "}" );

  const map<int,SampleFingerprint> before = fingerprints( orig ), after = fingerprints( reloaded );
  for( const auto &sample_fp : before )
  {
    const auto pos = after.find( sample_fp.first );
    BOOST_CHECK_MESSAGE( (pos != end(after)) && (pos->second == sample_fp.second),
                         label << ": sample " << sample_fp.first << " was " << sample_fp.second
                         << " but after reload is "
                         << ((pos == end(after)) ? string("missing") : [&pos](){
                              stringstream s; s << pos->second; return s.str(); }()) );
  }
}//check_same_samples(...)


/** Amplitudes of the peaks stored for each sample-number set. */
map<set<int>,vector<double>> peak_amplitudes( const SpecMeas &meas, const bool auto_search )
{
  map<set<int>,vector<double>> answer;
  const set<set<int>> sample_sets = auto_search ? meas.sampleNumsWithAutomatedSearchPeaks()
                                                : meas.sampleNumsWithPeaks();
  for( const set<int> &samples : sample_sets )
  {
    const shared_ptr<const SpecMeas::PeakDeque> peaks = auto_search ? meas.automatedSearchPeaks( samples )
                                                                    : meas.peaks( samples );
    vector<double> &amps = answer[samples];
    if( peaks )
    {
      for( const shared_ptr<const PeakDef> &p : *peaks )
        amps.push_back( p->amplitude() );
    }
  }//for( loop over sample sets )

  return answer;
}//peak_amplitudes(...)


void check_same_state( const SpecMeas &orig, const SpecMeas &reloaded, const string &label )
{
  BOOST_CHECK_MESSAGE( peak_amplitudes( orig, false ) == peak_amplitudes( reloaded, false ),
                       label << ": user peaks are not on the same sample numbers after reload" );
  BOOST_CHECK_MESSAGE( peak_amplitudes( orig, true ) == peak_amplitudes( reloaded, true ),
                       label << ": auto-search peaks are not on the same sample numbers after reload" );
  BOOST_CHECK_MESSAGE( orig.displayedSampleNumbers() == reloaded.displayedSampleNumbers(),
                       label << ": displayed samples were {"
                       << SpecUtils::sequencesToBriefString(orig.displayedSampleNumbers()) << "}, after reload {"
                       << SpecUtils::sequencesToBriefString(reloaded.displayedSampleNumbers()) << "}" );
#if( USE_LLM_INTERFACE )
  BOOST_CHECK_MESSAGE( orig.sampleNumsWithLlmHistory() == reloaded.sampleNumsWithLlmHistory(),
                       label << ": LLM history is not on the same sample numbers after reload" );
#endif
}//check_same_state(...)


string write_n42( const SpecMeas &meas )
{
  stringstream strm;
  BOOST_REQUIRE( meas.write_2012_N42( strm ) );
  return strm.str();
}


shared_ptr<SpecMeas> load_n42( string data )
{
  auto meas = make_shared<SpecMeas>();
  BOOST_REQUIRE( !data.empty() );
  BOOST_REQUIRE( meas->load_N42_from_data( &data[0], &data[0] + data.size() ) );
  return meas;
}


/** The "id" attributes of the given element, in document order, anywhere in the document. */
vector<string> ids_of( const string &xml, const string &element )
{
  vector<string> answer;
  const string tag = "<" + element + " ";
  for( size_t pos = xml.find( tag ); pos != string::npos; pos = xml.find( tag, pos + 1 ) )
  {
    const size_t end_pos = xml.find( '>', pos );
    const size_t id_pos = xml.find( " id=\"", pos );
    if( (id_pos == string::npos) || (id_pos > end_pos) )
      continue;
    const size_t value_start = id_pos + 5;
    answer.push_back( xml.substr( value_start, xml.find( '"', value_start ) - value_start ) );
  }
  return answer;
}//ids_of(...)


/** Writes, reloads and compares a file and its state; returns the written XML. */
string check_round_trip( const shared_ptr<SpecMeas> &meas, const string &label )
{
  const string xml = write_n42( *meas );
  const shared_ptr<SpecMeas> reloaded = load_n42( xml );
  check_same_samples( *meas, *reloaded, label );
  check_same_state( *meas, *reloaded, label );
  return xml;
}//check_round_trip(...)


/** The legacy fixture: written by the N42 writer from before sample numbers were preserved for
 search-mode files (positional Spectrum ids, no <AllSampleNumbers>).  Created by this test if
 missing, so it must be (re)created by running this test against the old code; the test checks it
 really has the old form.
 */
string legacy_fixture_path()
{
  set_data_dir();
  const string dir = SpecUtils::append_path( g_test_file_dir, "SpecMeasN42RoundTrip" );
  if( !SpecUtils::is_directory(dir) )
    BOOST_REQUIRE_EQUAL( SpecUtils::create_directory(dir), 1 );
  return SpecUtils::append_path( dir, "legacy_v2_search_bg_peaks.n42" );
}


shared_ptr<SpecMeas> legacy_fixture_source()
{
  shared_ptr<SpecMeas> meas = create_file( search_mode_layout( range_inclusive( 1, 11 ) ) );
  add_state( meas, {6}, {2,3,4,5,6,7,8,9,10,11} );
  return meas;
}


string read_file( const string &path )
{
  ifstream input( path.c_str(), ios::binary );
  BOOST_REQUIRE_MESSAGE( input.is_open(), "Could not open " << path );
  stringstream strm;
  strm << input.rdbuf();
  return strm.str();
}
}//namespace


BOOST_AUTO_TEST_CASE( SearchModeLeadingBackgroundFromOne )
{
  set_data_dir();

  // Make sure the legacy fixture exists before anything else (see `legacy_fixture_path()`).
  const string legacy_path = legacy_fixture_path();
  if( !SpecUtils::is_file( legacy_path ) )
  {
    ofstream out( legacy_path.c_str(), ios::binary );
    BOOST_REQUIRE( out.is_open() );
    BOOST_REQUIRE( legacy_fixture_source()->write_2012_N42( out ) );
  }

  const shared_ptr<SpecMeas> meas = create_file( search_mode_layout( range_inclusive( 1, 11 ) ) );
  BOOST_REQUIRE( meas->passthrough() );
  add_state( meas, {6}, {2,3,4,5,6,7,8,9,10,11} );

  const string xml = check_round_trip( meas, "search-mode 1..11" );

  // The RadMeasurement ids GADRAS relies on are unchanged...
  vector<string> expected_meas_ids{ "Background" };
  for( int i = 1; i <= 10; ++i )
    expected_meas_ids.push_back( "Survey" + std::to_string(i) );
  const vector<string> meas_ids = ids_of( xml, "RadMeasurement" );
  BOOST_CHECK_EQUAL_COLLECTIONS( begin(expected_meas_ids), end(expected_meas_ids),
                                 begin(meas_ids), end(meas_ids) );

  // ... but the Spectrum ids carry the real sample numbers.
  vector<string> expected_spec_ids;
  for( int i = 1; i <= 11; ++i )
    expected_spec_ids.push_back( "Sample" + std::to_string(i) + "Det0Spectrum" );
  const vector<string> spec_ids = ids_of( xml, "Spectrum" );
  BOOST_CHECK_EQUAL_COLLECTIONS( begin(expected_spec_ids), end(expected_spec_ids),
                                 begin(spec_ids), end(spec_ids) );
}//BOOST_AUTO_TEST_CASE( SearchModeLeadingBackgroundFromOne )


BOOST_AUTO_TEST_CASE( SearchModeOtherNumberings )
{
  set_data_dir();

  {// Numbered from 5
    const shared_ptr<SpecMeas> meas = create_file( search_mode_layout( range_inclusive( 5, 15 ) ) );
    BOOST_REQUIRE( meas->passthrough() );
    add_state( meas, {10}, {6,7,8,9,10,11,12,13,14,15} );
    check_round_trip( meas, "search-mode 5..15" );
  }

  {// Gaps in the numbering
    const shared_ptr<SpecMeas> meas = create_file( search_mode_layout( {1, 2, 3, 4, 7, 8, 9, 10, 12, 13, 14} ) );
    BOOST_REQUIRE( meas->passthrough() );
    add_state( meas, {8}, {2,3,4,7,8,9,10,12,13,14} );
    check_round_trip( meas, "search-mode with gaps" );
  }

  {// Two detectors, the second with neutrons (exercises the neutron <GrossCounts> ids)
    const shared_ptr<SpecMeas> meas = create_file( search_mode_layout( range_inclusive( 1, 11 ) ), 2 );
    BOOST_REQUIRE( meas->passthrough() );
    add_state( meas, {6}, {2,3,4,5,6,7,8,9,10,11} );
    check_round_trip( meas, "search-mode 2 detectors with neutrons" );
  }
}//BOOST_AUTO_TEST_CASE( SearchModeOtherNumberings )


BOOST_AUTO_TEST_CASE( NotSearchModeControl )
{
  set_data_dir();

  // Same layout with 60 s samples, so not search-mode: always round-tripped correctly.
  const shared_ptr<SpecMeas> meas = create_file( search_mode_layout( range_inclusive( 1, 11 ), 60.0f ) );
  BOOST_REQUIRE( !meas->passthrough() );
  add_state( meas, {6}, {2,3,4,5,6,7,8,9,10,11} );
  const string xml = check_round_trip( meas, "60 s samples" );

  const vector<string> meas_ids = ids_of( xml, "RadMeasurement" );
  BOOST_REQUIRE( !meas_ids.empty() );
  BOOST_CHECK_EQUAL( meas_ids.front(), "Sample1" );
}//BOOST_AUTO_TEST_CASE( NotSearchModeControl )


BOOST_AUTO_TEST_CASE( LegacyFileFromBeforeFix )
{
  set_data_dir();

  const string legacy_path = legacy_fixture_path();
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( legacy_path ),
                         "Missing " << legacy_path << " - run SearchModeLeadingBackgroundFromOne first" );

  const string xml = read_file( legacy_path );

  // Make sure the fixture really was written by the old writer.
  const vector<string> spec_ids = ids_of( xml, "Spectrum" );
  BOOST_REQUIRE( !spec_ids.empty() );
  BOOST_REQUIRE_MESSAGE( spec_ids.front() == "BackgroundDet0Spectrum",
                         "Fixture does not have the old positional Spectrum ids; delete it and"
                         " re-create it by running this test against the old code" );
  BOOST_REQUIRE( xml.find( "<AllSampleNumbers>" ) == string::npos );

  const shared_ptr<SpecMeas> orig = legacy_fixture_source();
  const shared_ptr<SpecMeas> reloaded = load_n42( xml );
  check_same_samples( *orig, *reloaded, "legacy fixture" );
  check_same_state( *orig, *reloaded, "legacy fixture" );
}//BOOST_AUTO_TEST_CASE( LegacyFileFromBeforeFix )


BOOST_AUTO_TEST_CASE( StoredSampleNumbersReconcile )
{
  set_data_dir();

  const shared_ptr<SpecMeas> meas = create_file( search_mode_layout( range_inclusive( 1, 11 ) ) );
  add_state( meas, {6}, {2,3,4,5,6,7,8,9,10,11} );
  const string xml = write_n42( *meas );
  BOOST_REQUIRE( xml.find( "<AllSampleNumbers>1 2 3 4 5 6 7 8 9 10 11</AllSampleNumbers>" ) != string::npos );

  // Put the old positional Spectrum ids back, as if SpecUtils changed its naming: the
  //  <AllSampleNumbers> InterSpec stores must still restore the original numbers.
  {
    string positional = xml;
    SpecUtils::ireplace_all( positional, "id=\"Sample1Det0Spectrum\"", "id=\"BackgroundDet0Spectrum\"" );
    for( int i = 2; i <= 11; ++i )
      SpecUtils::ireplace_all( positional, ("id=\"Sample" + std::to_string(i) + "Det0Spectrum\"").c_str(),
                               ("id=\"Survey" + std::to_string(i - 1) + "Det0Spectrum\"").c_str() );
    BOOST_REQUIRE( positional != xml );

    const shared_ptr<SpecMeas> reloaded = load_n42( positional );
    check_same_samples( *meas, *reloaded, "positional ids + AllSampleNumbers" );
    check_same_state( *meas, *reloaded, "positional ids + AllSampleNumbers" );
  }

  // Same, but also drop <AllSampleNumbers> and mark it as the previous version: the legacy
  //  handling must kick in.
  {
    string legacy = xml;
    SpecUtils::ireplace_all( legacy, "id=\"Sample1Det0Spectrum\"", "id=\"BackgroundDet0Spectrum\"" );
    for( int i = 2; i <= 11; ++i )
      SpecUtils::ireplace_all( legacy, ("id=\"Sample" + std::to_string(i) + "Det0Spectrum\"").c_str(),
                               ("id=\"Survey" + std::to_string(i - 1) + "Det0Spectrum\"").c_str() );
    SpecUtils::ireplace_all( legacy, "<AllSampleNumbers>1 2 3 4 5 6 7 8 9 10 11</AllSampleNumbers>", "" );
    BOOST_REQUIRE( legacy.find( "<AllSampleNumbers>" ) == string::npos );

    const shared_ptr<SpecMeas> reloaded = load_n42( legacy );
    check_same_samples( *meas, *reloaded, "synthesized legacy" );
    check_same_state( *meas, *reloaded, "synthesized legacy" );
  }

  // A stored list that can't be matched up: keep what was parsed, and warn.
  {
    string mismatched = xml;
    SpecUtils::ireplace_all( mismatched, "<AllSampleNumbers>1 2 3", "<AllSampleNumbers>0 1 2 3" );
    const shared_ptr<SpecMeas> reloaded = load_n42( mismatched );
    BOOST_CHECK( reloaded->sample_numbers() == meas->sample_numbers() );

    bool warned = false;
    for( const string &warning : reloaded->parse_warnings() )
      warned |= SpecUtils::icontains( warning, "sample numbers" );
    BOOST_CHECK_MESSAGE( warned, "No parse warning for mismatched <AllSampleNumbers>" );
  }
}//BOOST_AUTO_TEST_CASE( StoredSampleNumbersReconcile )


BOOST_AUTO_TEST_CASE( ChangeSampleNumbersRemapsState )
{
  set_data_dir();

  const shared_ptr<SpecMeas> meas = create_file( search_mode_layout( range_inclusive( 1, 11 ) ) );
  add_state( meas, {6}, {2,3,4,5,6,7,8,9,10,11} );

  // Shift every sample up by 10 (overlapping from/to ranges are fine).
  vector<pair<int,int>> from_to;
  for( const int sample : meas->sample_numbers() )
    from_to.emplace_back( sample, sample + 10 );
  meas->change_sample_numbers( from_to );

  BOOST_CHECK( meas->sample_numbers() == set<int>( {11,12,13,14,15,16,17,18,19,20,21} ) );
  BOOST_CHECK( meas->sampleNumsWithPeaks() == set<set<int>>( { {16}, {12,13,14,15,16,17,18,19,20,21} } ) );
  BOOST_CHECK( meas->sampleNumsWithAutomatedSearchPeaks() == set<set<int>>( { {16} } ) );
  BOOST_CHECK( meas->displayedSampleNumbers() == set<int>( {12,13,14,15,16,17,18,19,20,21} ) );
#if( USE_LLM_INTERFACE )
  BOOST_CHECK( meas->sampleNumsWithLlmHistory() == set<set<int>>( { {16} } ) );
#endif
}//BOOST_AUTO_TEST_CASE( ChangeSampleNumbersRemapsState )


BOOST_AUTO_TEST_CASE( RealFilesRoundTrip )
{
  set_data_dir();

  string data_dir = g_data_dir;
  BOOST_REQUIRE( SpecUtils::make_canonical_path( data_dir ) );
  const string repo_root = SpecUtils::parent_path( data_dir );
  vector<string> files{ SpecUtils::append_path( repo_root, "example_spectra/passthrough.n42" ) };
  const string time_history_dir = SpecUtils::append_path( g_test_file_dir, "llm_time_history" );
  for( const string &f : SpecUtils::recursive_ls( time_history_dir ) )
    files.push_back( f );

  BOOST_REQUIRE( files.size() > 1 );

  for( const string &path : files )
  {
    BOOST_TEST_MESSAGE( "Round-tripping " << path );
    auto meas = make_shared<SpecMeas>();
    BOOST_REQUIRE_MESSAGE( meas->load_file( path, SpecUtils::ParserType::Auto, path ), "Failed to load " << path );

    const set<int> &samples = meas->sample_numbers();
    BOOST_REQUIRE( !samples.empty() );
    const int middle = *std::next( begin(samples), samples.size()/2 );
    add_state( meas, {middle}, samples );

    check_round_trip( meas, SpecUtils::filename(path) );
  }//for( const string &path : files )
}//BOOST_AUTO_TEST_CASE( RealFilesRoundTrip )
