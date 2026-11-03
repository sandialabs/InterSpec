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

#if defined(_WIN32)
  #define WIN32_LEAN_AND_MEAN
  #include <winsock2.h>
  #include <windows.h>
#endif

#include "InterSpec_config.h"

#include <set>
#include <deque>
#include <string>
#include <memory>
#include <vector>

#define BOOST_TEST_MODULE TestSpecMeasManager
#include <boost/test/included/unit_test.hpp>

#include <Wt/WApplication.h>
#include <Wt/WModelIndex.h>
#include <Wt/Test/WTestEnvironment.h>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/RowStretchTreeView.h"

using namespace std;

/** Tests of how the "Spectrum Manager" (SpecMeasManager) and its file model track the displayed
 spectrum files, and of its sum / selection handling.

 The "Spectrum Files" tab finds the displayed files by looking them up in the file model, so the
 central invariant checked here is that every displayed SpecMeas is in the model.
 */

namespace
{
void set_test_data_dir()
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
    for( const char * const d : { "data", "../data", "../../data", "../../../data" } )
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
}//set_test_data_dir()


/** A minimal InterSpec session. */
class SessionFixture
{
public:
  std::unique_ptr<Wt::Test::WTestEnvironment> m_env;
  InterSpecApp *m_app = nullptr;
  std::unique_ptr<Wt::WApplication::UpdateLock> m_lock;
  InterSpec *m_interspec = nullptr;
  SpecMeasManager *m_manager = nullptr;

  SessionFixture()
  {
    set_test_data_dir();

    string app_root = SpecUtils::append_path( InterSpec::staticDataDirectory(), ".." );
    app_root = SpecUtils::lexically_normalize_path( app_root );

    m_env.reset( new Wt::Test::WTestEnvironment( "", "", Wt::EntryPointType::Application ) );
    m_env->setAppRoot( app_root );
    m_app = new InterSpecApp( *m_env );
    m_lock.reset( new Wt::WApplication::UpdateLock(m_app) );
    m_interspec = m_app->viewer();
    BOOST_REQUIRE( m_interspec );
    m_manager = m_interspec->fileManager();
    BOOST_REQUIRE( m_manager && m_manager->model() && m_manager->treeView() );
  }

  ~SessionFixture()
  {
    // `m_app` is owned by the environment
    m_lock.reset();
    m_env.reset();
  }

  /** A file of three records (sample numbers 1, 2, 3; the third marked as a background), each a
   copy of a reference spectrum, so tests dont depend on a particular multi-record test file.
   */
  static shared_ptr<SpecMeas> make_three_record_file( const string &filename )
  {
    const string spec = SpecUtils::append_path( InterSpec::staticDataDirectory(),
                    "reference_spectra/Common_Field_Nuclides/Detective X/Br82_Unshielded.txt" );
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file(spec), "test spectrum not at '" << spec << "'" );

    SpecMeas src;
    BOOST_REQUIRE( src.load_file( spec, SpecUtils::ParserType::Auto, spec ) );
    BOOST_REQUIRE( src.num_measurements() > 0 );
    const shared_ptr<const SpecUtils::Measurement> base = src.measurement_at_index( 0 );
    BOOST_REQUIRE( base );

    auto meas = make_shared<SpecMeas>();
    for( int sample = 1; sample <= 3; ++sample )
    {
      auto m = make_shared<SpecUtils::Measurement>( *base );
      m->set_sample_number( sample );
      m->set_source_type( (sample == 3) ? SpecUtils::SourceType::Background
                                        : SpecUtils::SourceType::Foreground );
      meas->add_measurement( m, false );
    }

    meas->cleanup_after_load( SpecUtils::SpecFile::CleanupAfterLoadFlags::StandardCleanup
                             | SpecUtils::SpecFile::CleanupAfterLoadFlags::DontChangeOrReorderSamples );
    meas->set_filename( filename );
    BOOST_REQUIRE( (meas->sample_numbers() == set<int>{1, 2, 3}) );

    return meas;
  }//make_three_record_file(...)
};//class SessionFixture
}//namespace


// `setSpectrum` should tell listeners the samples it actually displays (here, the files valid
//  foreground records), not the empty set it was asked for; the "Spectrum Files" tab and undo/redo
//  both use what it emits.  And a displayed file that isnt in the file model gets added to it.
BOOST_FIXTURE_TEST_CASE( SetSpectrumEmitsDisplayedSamplesAndAddsToModel, SessionFixture )
{
  const shared_ptr<SpecMeas> meas = make_three_record_file( "three_records.txt" );
  BOOST_REQUIRE( !m_manager->model()->index( meas ).isValid() );

  set<int> emitted_samples;
  shared_ptr<SpecMeas> emitted_meas;
  Wt::Signals::connection conn = m_interspec->displayedSpectrumChanged().connect(
      [&]( SpecUtils::SpectrumType type, shared_ptr<SpecMeas> m, set<int> samples, vector<string> ){
        if( type == SpecUtils::SpectrumType::Foreground )
        {
          emitted_meas = m;
          emitted_samples = samples;
        }
      } );

  m_interspec->setSpectrum( meas, {}, SpecUtils::SpectrumType::Foreground, {} );
  conn.disconnect();

  const set<int> &displayed = m_interspec->displayedSamples( SpecUtils::SpectrumType::Foreground );
  BOOST_CHECK( emitted_meas == meas );
  BOOST_CHECK( !displayed.empty() );
  BOOST_CHECK( (emitted_samples == displayed) );
  BOOST_CHECK( (displayed == set<int>{1, 2}) ); // the background record isnt a valid foreground

  BOOST_CHECK_MESSAGE( m_manager->model()->index( meas ).isValid(),
                       "The displayed file was not added to the file model" );
}//SetSpectrumEmitsDisplayedSamplesAndAddsToModel


// Removing a file from the Spectrum Manager should un-display it from every spectrum type.
BOOST_FIXTURE_TEST_CASE( RemovingAFileUndisplaysItEverywhere, SessionFixture )
{
  const shared_ptr<SpecMeas> meas = make_three_record_file( "fore_and_back.txt" );
  m_manager->addFile( meas->filename(), meas );

  m_interspec->setSpectrum( meas, {1}, SpecUtils::SpectrumType::Foreground, {} );
  m_interspec->setSpectrum( meas, {3}, SpecUtils::SpectrumType::Background, {} );
  BOOST_REQUIRE( m_interspec->measurment(SpecUtils::SpectrumType::Background) == meas );

  m_manager->removeSpecMeas( meas, true );

  BOOST_CHECK( !m_interspec->measurment(SpecUtils::SpectrumType::Foreground) );
  BOOST_CHECK( !m_interspec->measurment(SpecUtils::SpectrumType::Background) );
  BOOST_CHECK( !m_manager->model()->index( meas ).isValid() );

  // And `addFile` of a file already in the model returns its existing row, rather than adding one.
  const shared_ptr<SpecMeas> other = make_three_record_file( "other.txt" );
  const shared_ptr<SpectraFileHeader> header = m_manager->addFile( other->filename(), other );
  const int nrows = m_manager->model()->rowCount();
  BOOST_CHECK( m_manager->addFile( other->filename(), other ) == header );
  BOOST_CHECK_EQUAL( m_manager->model()->rowCount(), nrows );
}//RemovingAFileUndisplaysItEverywhere


// A selected file row means all its records, unless some of its records are also selected (e.g.,
//  the user selected the file, then de-selected a record); in which case only those records.
BOOST_FIXTURE_TEST_CASE( SelectionHonorsSelectedRecords, SessionFixture )
{
  const shared_ptr<SpecMeas> meas = make_three_record_file( "selection.txt" );
  const shared_ptr<SpectraFileHeader> header = m_manager->addFile( meas->filename(), meas );
  const Wt::WModelIndex file_index = m_manager->model()->index( header );
  BOOST_REQUIRE( file_index.isValid() );
  BOOST_REQUIRE_EQUAL( m_manager->model()->rowCount( file_index ), 3 );

  RowStretchTreeView * const table = m_manager->treeView();

  table->setSelectedIndexes( { file_index } );
  BOOST_CHECK( (m_manager->selectedSampleNumbers() == set<int>{1, 2, 3}) );
  BOOST_CHECK( m_manager->selectedFile() == header );

  table->setSelectedIndexes( { file_index, file_index.child(0, 0), file_index.child(2, 0) } );
  BOOST_CHECK( (m_manager->selectedSampleNumbers() == set<int>{1, 3}) );

  // Just records, without the file row
  table->setSelectedIndexes( { file_index.child(1, 0) } );
  BOOST_CHECK( (m_manager->selectedSampleNumbers() == set<int>{2}) );
  BOOST_CHECK( m_manager->selectedFile() == header );

  // Records from two files is not a single-file selection
  const shared_ptr<SpecMeas> other = make_three_record_file( "selection_other.txt" );
  const shared_ptr<SpectraFileHeader> other_header = m_manager->addFile( other->filename(), other );
  const Wt::WModelIndex other_index = m_manager->model()->index( other_header );
  table->setSelectedIndexes( { file_index, other_index } );
  BOOST_CHECK( !m_manager->selectedFile() );
  BOOST_CHECK( m_manager->selectedSampleNumbers().empty() );

  // ...and asking to display such a selection shouldnt un-load whats currently displayed.
  m_interspec->setSpectrum( meas, {1}, SpecUtils::SpectrumType::Foreground, {} );
  m_manager->loadSelected( SpecUtils::SpectrumType::Foreground, false );
  BOOST_CHECK( m_interspec->measurment(SpecUtils::SpectrumType::Foreground) == meas );
}//SelectionHonorsSelectedRecords


// Summing a whole file gives a new one-record file, which must not keep the source files peaks,
//  displayed sample numbers, etc. (which are keyed by the source files sample numbers).
BOOST_FIXTURE_TEST_CASE( SumOfWholeFileDropsPerSampleState, SessionFixture )
{
  const shared_ptr<SpecMeas> meas = make_three_record_file( "to_sum.txt" );
  const shared_ptr<SpectraFileHeader> header = m_manager->addFile( meas->filename(), meas );

  // Displaying record 2 records {2} as the files displayed samples.
  m_interspec->setSpectrum( meas, {2}, SpecUtils::SpectrumType::Foreground, {} );
  BOOST_REQUIRE( (meas->displayedSampleNumbers() == set<int>{2}) );

  // A peak fit to record 1 only
  const deque<shared_ptr<const PeakDef>> peaks{ make_shared<const PeakDef>( 776.5, 1.5, 1000.0 ) };
  meas->setPeaks( peaks, {1} );

  const float single_live_time = meas->measurement( 1, meas->detector_names().front() )->live_time();

  const int nrows_before = m_manager->model()->rowCount();
  m_manager->treeView()->setSelectedIndexes( { m_manager->model()->index( header ) } );
  m_manager->sumSelectedSpectra();
  BOOST_REQUIRE_EQUAL( m_manager->model()->rowCount(), nrows_before + 1 );

  const shared_ptr<SpectraFileHeader> sum_header = m_manager->model()->fileHeader( nrows_before );
  BOOST_REQUIRE( sum_header );
  const shared_ptr<const SpecMeas> sum = sum_header->parseFile();
  BOOST_REQUIRE( sum && (sum != meas) );

  BOOST_CHECK_EQUAL( sum->num_measurements(), 1 );
  BOOST_REQUIRE_EQUAL( sum->sample_numbers().size(), 1 );
  BOOST_CHECK_CLOSE( sum->gamma_live_time(), 3.0f*single_live_time, 1.0E-3 );

  const set<int> sum_samples = sum->sample_numbers();
  const auto sum_peaks = sum->peaks( sum_samples );
  BOOST_CHECK_MESSAGE( !sum_peaks || sum_peaks->empty(),
                       "The sum picked up the peaks of a record of the source file" );
  BOOST_CHECK_MESSAGE( sum->displayedSampleNumbers().empty(),
                       "The sum kept the source files displayed sample numbers" );

  // The source file is untouched
  const auto src_peaks = std::const_pointer_cast<const SpecMeas>( meas )->peaks( {1} );
  BOOST_CHECK( src_peaks && (src_peaks->size() == 1) );
}//SumOfWholeFileDropsPerSampleState


// `displayFile` shows a files first foreground record as the foreground, and loads its single
//  background record as the background.
BOOST_FIXTURE_TEST_CASE( DisplayFilePicksForegroundAndBackgroundRecords, SessionFixture )
{
  const shared_ptr<SpecMeas> meas = make_three_record_file( "display_file.txt" );
  const shared_ptr<SpectraFileHeader> header = m_manager->addFile( meas->filename(), meas );
  const int row = m_manager->model()->index( header ).row();
  BOOST_REQUIRE( row >= 0 );

  m_manager->displayFile( row, meas, SpecUtils::SpectrumType::Foreground, false, false,
                          SpecMeasManager::VariantChecksToDo::None );

  BOOST_CHECK( m_interspec->measurment(SpecUtils::SpectrumType::Foreground) == meas );
  BOOST_CHECK( (m_interspec->displayedSamples(SpecUtils::SpectrumType::Foreground) == set<int>{1}) );
  BOOST_CHECK( m_interspec->measurment(SpecUtils::SpectrumType::Background) == meas );
  BOOST_CHECK( (m_interspec->displayedSamples(SpecUtils::SpectrumType::Background) == set<int>{3}) );
}//DisplayFilePicksForegroundAndBackgroundRecords


// When a file is replaced by another version of itself (e.g., "resume previous work"), the new
//  version takes its place everywhere it is displayed, so removing the old one un-displays nothing.
BOOST_FIXTURE_TEST_CASE( ReplacedFileTakesItsPlaceEverywhere, SessionFixture )
{
  const shared_ptr<SpecMeas> original = make_three_record_file( "original.txt" );
  const shared_ptr<SpecMeas> resumed = make_three_record_file( "resumed.txt" );
  m_manager->addFile( original->filename(), original );
  m_manager->addFile( resumed->filename(), resumed );

  m_interspec->setSpectrum( original, {2}, SpecUtils::SpectrumType::Foreground, {} );
  m_interspec->setSpectrum( original, {3}, SpecUtils::SpectrumType::Background, {} );

  m_manager->displayInPlaceOf( original, resumed );
  m_manager->removeSpecMeas( original, true );

  BOOST_CHECK( m_interspec->measurment(SpecUtils::SpectrumType::Foreground) == resumed );
  BOOST_CHECK( (m_interspec->displayedSamples(SpecUtils::SpectrumType::Foreground) == set<int>{2}) );
  BOOST_CHECK( m_interspec->measurment(SpecUtils::SpectrumType::Background) == resumed );
  BOOST_CHECK( (m_interspec->displayedSamples(SpecUtils::SpectrumType::Background) == set<int>{3}) );
  BOOST_CHECK( !m_manager->model()->index( original ).isValid() );
  BOOST_CHECK( m_manager->model()->index( resumed ).isValid() );
}//ReplacedFileTakesItsPlaceEverywhere
