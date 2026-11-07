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

#include <map>
#include <set>
#include <memory>
#include <string>
#include <vector>

#define BOOST_TEST_MODULE TestAppStateDb
#include <boost/test/included/unit_test.hpp>

#include <Wt/WApplication.h>
#include <Wt/Dbo/Dbo.h>
#include <Wt/Test/WTestEnvironment.h>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/SpectraFileModel.h"

using namespace std;

/** Tests of saving and restoring the application state to/from the database (the "Store",
 "Store As...", and "Tag..." menu items, auto-saved work, and end-of-session states).

 The central invariants (checked after every test, see `check_invariants`) are that every state
 exclusively owns the spectrum file rows it references, that no state-owned rows are orphaned, and
 that only explicitly storing a state ever changes that states rows.
 */

#if( USE_DB_TO_STORE_SPECTRA )
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


/** What the database holds for a `UserState`, read independently of the apps database session. */
struct StateRow
{
  long long int id = -1;
  int type = -1;
  long long int parent = -1;
  long long int user = -1;
  int fg = -1, bg = -1, sec = -1;
  string fgSamples, bgSamples, secSamples;
  int bgSubMode = -1;
};//struct StateRow


/** What the database holds for a `UserFileInDb` (and its data). */
struct FileRow
{
  long long int id = -1;
  bool partOfState = false;
  bool hasSnapshotParent = false;
  string uuid;
  vector<unsigned char> data;
};//struct FileRow


struct DbSnapshot
{
  map<long long int,StateRow> states;
  map<long long int,FileRow> files;

  /** The (non-`kEndOfSessionTemp`) top-level user-saved states. */
  vector<StateRow> topLevelUserStates() const
  {
    vector<StateRow> answer;
    for( const auto &id_state : states )
    {
      if( (id_state.second.type == UserState::kUserSaved) && (id_state.second.parent < 0) )
        answer.push_back( id_state.second );
    }
    return answer;
  }

  vector<StateRow> children( const long long int parent, const int type ) const
  {
    vector<StateRow> answer;
    for( const auto &id_state : states )
    {
      if( (id_state.second.parent == parent) && (id_state.second.type == type) )
        answer.push_back( id_state.second );
    }
    return answer;
  }

  size_t numStandaloneFilesWithUuid( const string &uuid ) const
  {
    size_t n = 0;
    for( const auto &id_file : files )
      n += (!id_file.second.partOfState && (id_file.second.uuid == uuid));
    return n;
  }
};//struct DbSnapshot


/** Reads the database through its own session, so nothing cached by the app is seen.
 Must not be called while the app holds a transaction (there is a single in-memory connection).
 */
DbSnapshot snapshot()
{
  DbSnapshot answer;

  DataBaseUtils::DbSession probe;
  DataBaseUtils::DbTransaction transaction( probe );
  Wt::Dbo::Session &session = *probe.session();

  Wt::Dbo::collection<Wt::Dbo::ptr<UserState>> states = session.find<UserState>();
  for( const Wt::Dbo::ptr<UserState> &st : states )
  {
    StateRow row;
    row.id = st.id();
    row.type = static_cast<int>( st->stateType );
    row.parent = st->snapshotTagParent ? st->snapshotTagParent.id() : -1;
    row.user = st->user ? st->user.id() : -1;
    row.fg = st->foregroundId;
    row.bg = st->backgroundId;
    row.sec = st->secondForegroundId;
    row.fgSamples = st->foregroundSampleNumsCsvIds;
    row.bgSamples = st->backgroundSampleNumsCsvIds;
    row.secSamples = st->secondForegroundSampleNumsCsvIds;
    row.bgSubMode = static_cast<int>( st->backgroundSubMode );
    answer.states[row.id] = row;
  }//for( loop over states )

  Wt::Dbo::collection<Wt::Dbo::ptr<UserFileInDb>> files = session.find<UserFileInDb>();
  for( const Wt::Dbo::ptr<UserFileInDb> &f : files )
  {
    FileRow row;
    row.id = f.id();
    row.partOfState = f->isPartOfSaveState;
    row.hasSnapshotParent = !!f->snapshotParent;
    row.uuid = f->uuid;
    answer.files[row.id] = row;
  }//for( loop over files )

  // The spectrum data is read with its own query (rather than through each files collection)
  Wt::Dbo::collection<Wt::Dbo::ptr<UserFileInDbData>> datas = session.find<UserFileInDbData>();
  for( const Wt::Dbo::ptr<UserFileInDbData> &d : datas )
  {
    const long long int fileid = d->fileInfo.id();
    const auto pos = answer.files.find( fileid );
    BOOST_CHECK_MESSAGE( pos != end(answer.files), "Spectrum data for missing file " << fileid );
    if( pos != end(answer.files) )
      pos->second.data.assign( begin(d->fileData), end(d->fileData) );
  }//for( loop over file data )

  transaction.commit();

  return answer;
}//DbSnapshot snapshot()


/** Checks the database invariants the state saving/loading code must maintain. */
void check_invariants( const DbSnapshot &snap, const bool allowMultipleAutoSaved = false )
{
  map<long long int,set<long long int>> statesUsingFile;

  for( const auto &id_state : snap.states )
  {
    const StateRow &st = id_state.second;

    for( const int fileid : { st.fg, st.bg, st.sec } )
    {
      if( fileid < 0 )
        continue;

      statesUsingFile[fileid].insert( st.id );

      const auto pos = snap.files.find( fileid );
      BOOST_CHECK_MESSAGE( pos != end(snap.files),
                           "State " << st.id << " references missing file " << fileid );
      if( pos == end(snap.files) )
        continue;

      BOOST_CHECK_MESSAGE( !pos->second.data.empty(), "State " << st.id << " file " << fileid << " has no data" );
      BOOST_CHECK_MESSAGE( pos->second.partOfState,
                           "State " << st.id << " file " << fileid << " not marked as part of a state" );
    }//for( loop over file ids )

    // Every child states top-level ancestor is a user-saved state
    if( st.parent >= 0 )
    {
      long long int top = st.parent;
      while( snap.states.count(top) && (snap.states.at(top).parent >= 0) )
        top = snap.states.at(top).parent;
      BOOST_CHECK_MESSAGE( snap.states.count(top) && (snap.states.at(top).type == UserState::kUserSaved),
                           "State " << st.id << " top-level ancestor " << top << " is not a user-saved state" );
    }//if( st.parent >= 0 )

    if( !allowMultipleAutoSaved && (st.type == UserState::kUserSaved) && (st.parent < 0) )
    {
      const size_t nauto = snap.children( st.id, UserState::kUserStateAutoSavedWork ).size();
      BOOST_CHECK_MESSAGE( nauto <= 1, "State " << st.id << " has " << nauto << " auto-saved states" );
    }
  }//for( loop over states )

  for( const auto &file_states : statesUsingFile )
  {
    BOOST_CHECK_MESSAGE( file_states.second.size() == 1,
                         "File " << file_states.first << " is shared by " << file_states.second.size() << " states" );
  }

  for( const auto &id_file : snap.files )
  {
    if( id_file.second.partOfState )
    {
      BOOST_CHECK_MESSAGE( statesUsingFile.count( id_file.first ),
                           "File " << id_file.first << " is part of a state, but no state references it" );
    }
  }//for( loop over files )
}//void check_invariants(...)


/** Checks the states, and the file rows they own, are identical in two snapshots. */
void check_same_states( const DbSnapshot &before, const DbSnapshot &after )
{
  BOOST_CHECK_EQUAL( before.states.size(), after.states.size() );
  for( const auto &id_state : before.states )
  {
    const auto pos = after.states.find( id_state.first );
    BOOST_CHECK_MESSAGE( pos != end(after.states), "State " << id_state.first << " removed" );
    if( pos == end(after.states) )
      continue;
    
    const StateRow &a = id_state.second, &b = pos->second;
    BOOST_CHECK_MESSAGE( (a.type == b.type) && (a.parent == b.parent) && (a.fg == b.fg) && (a.bg == b.bg)
                         && (a.sec == b.sec) && (a.fgSamples == b.fgSamples) && (a.bgSamples == b.bgSamples)
                         && (a.secSamples == b.secSamples) && (a.bgSubMode == b.bgSubMode),
                         "State " << a.id << " was changed" );
  }//for( loop over states )
  
  for( const auto &id_file : before.files )
  {
    if( !id_file.second.partOfState )
      continue;
    const auto pos = after.files.find( id_file.first );
    BOOST_CHECK_MESSAGE( (pos != end(after.files)) && (pos->second.data == id_file.second.data),
                         "State file " << id_file.first << " was removed or changed" );
  }//for( loop over files )
}//void check_same_states(...)


/** A minimal InterSpec session, with helpers to drive the state saving. */
class StateFixture
{
public:
  std::unique_ptr<Wt::Test::WTestEnvironment> m_env;
  InterSpecApp *m_app = nullptr;
  std::unique_ptr<Wt::WApplication::UpdateLock> m_lock;
  InterSpec *m_interspec = nullptr;
  SpecMeasManager *m_manager = nullptr;

  StateFixture()
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
    BOOST_REQUIRE( m_manager && m_manager->model() );

    UserPreferences::setPreferenceValue<bool>( "AutoSaveSpectraToDb", true, m_interspec );
  }

  ~StateFixture()
  {
    // `m_app` is owned by the environment
    m_lock.reset();
    m_env.reset();

    // The in-memory database is shared by all the tests, so clean up our states and files (otherwise
    //  the next tests session would pick up our end-of-session state, and counts would be off).
    try
    {
      DataBaseUtils::DbSession cleanup;
      DataBaseUtils::DbTransaction transaction( cleanup );
      cleanup.session()->execute( "DELETE FROM UserState" );
      cleanup.session()->execute( "DELETE FROM UserFileInDb" );
      transaction.commit();
    }catch( std::exception &e )
    {
      BOOST_ERROR( "Failed to clean up database: " << e.what() );
    }
  }//~StateFixture()

  /** A file of three foreground records (sample numbers 1, 2, 3), each a copy of a reference spectrum. */
  static shared_ptr<SpecMeas> make_file( const string &filename )
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
      m->set_source_type( SpecUtils::SourceType::Foreground );
      meas->add_measurement( m, false );
    }

    meas->cleanup_after_load( SpecUtils::SpecFile::CleanupAfterLoadFlags::StandardCleanup
                             | SpecUtils::SpecFile::CleanupAfterLoadFlags::DontChangeOrReorderSamples );
    meas->set_filename( filename );
    meas->set_uuid( filename + "-uuid" );
    BOOST_REQUIRE( (meas->sample_numbers() == set<int>{1, 2, 3}) );

    return meas;
  }//make_file(...)

  /** A file too large to be stored in the database (its serialization exceeds
   `UserFileInDb::sm_maxFileSizeBytes`), with the same energy binning as `make_file`.
   */
  static shared_ptr<SpecMeas> make_large_file( const string &filename )
  {
    const shared_ptr<SpecMeas> small = make_file( filename );
    const shared_ptr<const SpecUtils::Measurement> base = small->measurement_at_index( 0 );
    BOOST_REQUIRE( base && (base->num_gamma_channels() > 16) );
    
    const size_t nchannel = base->num_gamma_channels();
    const auto counts = make_shared<const vector<float>>( nchannel, 12345.678f );
    const size_t nsample = 1 + (2 * UserFileInDb::sm_maxFileSizeBytes) / (8 * nchannel);
    
    auto meas = make_shared<SpecMeas>();
    for( size_t sample = 1; sample <= nsample; ++sample )
    {
      auto m = make_shared<SpecUtils::Measurement>( *base );
      m->set_sample_number( static_cast<int>(sample) );
      m->set_gamma_counts( counts, base->live_time(), base->real_time() );
      meas->add_measurement( m, false );
    }
    
    meas->cleanup_after_load( SpecUtils::SpecFile::CleanupAfterLoadFlags::StandardCleanup
                             | SpecUtils::SpecFile::CleanupAfterLoadFlags::DontChangeOrReorderSamples );
    meas->set_filename( filename );
    meas->set_uuid( filename + "-uuid" );
    
    return meas;
  }//make_large_file(...)

  /** Uses the apps database session (which flushes any changes it has pending), like the app does all the time. */
  void touchAppSession()
  {
    std::shared_ptr<DataBaseUtils::DbSession> sql = m_interspec->sql();
    DataBaseUtils::DbTransaction transaction( *sql );
    const size_t nstates = sql->session()->find<UserState>().resultList().size();
    BOOST_TEST_MESSAGE( "App session sees " << nstates << " states" );
    transaction.commit();
  }

  /** Adds a file to the file manager, and displays it as `type`. */
  shared_ptr<SpecMeas> show( const string &filename, const set<int> &samples,
                             const SpecUtils::SpectrumType type = SpecUtils::SpectrumType::Foreground )
  {
    const shared_ptr<SpecMeas> meas = make_file( filename );
    m_manager->addFile( meas->filename(), meas );
    m_interspec->setSpectrum( meas, samples, type, {} );
    BOOST_REQUIRE( m_interspec->measurment( type ) == meas );
    return meas;
  }//show(...)

  shared_ptr<SpecMeas> foreground() const
  {
    return m_interspec->measurment( SpecUtils::SpectrumType::Foreground );
  }

  /** Adds a peak to the foreground (which marks it as modified). */
  void addPeak( const double energy )
  {
    PeakDef peak( energy, 1.0, 1000.0 );
    peak.continuum()->setRange( energy - 6.0, energy + 6.0 );
    m_interspec->peakModel()->addNewPeak( peak );
  }

  size_t npeaks() const
  {
    return m_interspec->peakModel()->npeaks();
  }

  long long int link() const
  {
    return m_interspec->currentAppStateDbId();
  }

  /** Like the "Store As..." menu item; returns the new states id. */
  long long int storeAs( const string &name )
  {
    m_interspec->finishStoreStateInDb( name, "", Wt::Dbo::ptr<UserState>() );
    const long long int id = link();
    BOOST_REQUIRE_MESSAGE( id >= 0, "Store As '" << name << "' did not connect to a state" );
    return id;
  }

  /** Like the "Tag..." menu item. */
  void tag( const string &name )
  {
    const long long int parent = link();
    BOOST_REQUIRE( parent >= 0 );
    m_interspec->finishStoreStateInDb( name, "", stateById( parent ) );
    BOOST_CHECK_EQUAL( link(), parent ); // We stay connected to the parent
  }

  /** The state, from the apps database session. */
  Wt::Dbo::ptr<UserState> stateById( const long long int id )
  {
    std::shared_ptr<DataBaseUtils::DbSession> sql = m_interspec->sql();
    DataBaseUtils::DbTransaction transaction( *sql );
    Wt::Dbo::ptr<UserState> state = sql->session()->find<UserState>().where( "id = ?" ).bind( id ).resultValue();
    transaction.commit();
    BOOST_REQUIRE_MESSAGE( state, "No state with id " << id );
    return state;
  }

  void load( const long long int id )
  {
    BOOST_REQUIRE_NO_THROW( m_interspec->loadStateFromDb( stateById( id ) ) );
  }

  /** Displays the file as foreground the same way the "Spectrum Files" tab does
   (see `CompactFileManager::handleFileChangeRequest`).
   */
  void displayViaFileManager( const shared_ptr<SpectraFileHeader> &header )
  {
    BOOST_REQUIRE( header );
    const Wt::WModelIndex index = m_manager->model()->index( header );
    BOOST_REQUIRE( index.isValid() );
    const shared_ptr<SpecMeas> meas = header->parseFile();
    m_manager->displayFile( index.row(), meas, SpecUtils::SpectrumType::Foreground, false, true,
                            SpecMeasManager::VariantChecksToDo::None );
  }
};//class StateFixture
}//namespace


// A state with summed foreground and background samples should restore those samples (the sample
//  numbers are stored as comma-separated values, which used to be parsed as just the first number).
BOOST_FIXTURE_TEST_CASE( SummedSamplesRoundTrip, StateFixture )
{
  show( "fore.n42", {1, 2} );
  show( "back.n42", {1, 2}, SpecUtils::SpectrumType::Background );
  show( "second.n42", {2, 3}, SpecUtils::SpectrumType::SecondForeground );
  addPeak( 554.3 );
  addPeak( 776.5 );
  BOOST_REQUIRE_EQUAL( npeaks(), 2 );

  const long long int state = storeAs( "Summed" );
  BOOST_CHECK_EQUAL( snapshot().states.at(state).bgSamples, "1,2" );
  BOOST_CHECK_EQUAL( snapshot().states.at(state).secSamples, "2,3" );

  show( "other.n42", {1} );
  BOOST_CHECK_EQUAL( link(), -1 );

  const DbSnapshot beforeLoad = snapshot();
  load( state );
  check_same_states( beforeLoad, snapshot() ); // Loading a state doesnt change the database
  BOOST_CHECK( (m_interspec->displayedSamples( SpecUtils::SpectrumType::Foreground ) == set<int>{1, 2}) );
  BOOST_CHECK( (m_interspec->displayedSamples( SpecUtils::SpectrumType::Background ) == set<int>{1, 2}) );
  BOOST_CHECK( (m_interspec->displayedSamples( SpecUtils::SpectrumType::SecondForeground ) == set<int>{2, 3}) );
  BOOST_CHECK_EQUAL( npeaks(), 2 );
  BOOST_CHECK_EQUAL( link(), state );

  check_invariants( snapshot() );
}//SummedSamplesRoundTrip


// Foreground and background from the same file, with different samples.
BOOST_FIXTURE_TEST_CASE( SameFileForeAndBackRoundTrip, StateFixture )
{
  const shared_ptr<SpecMeas> meas = show( "fore_and_back.n42", {1, 2} );
  m_interspec->setSpectrum( meas, {3}, SpecUtils::SpectrumType::Background, {} );

  const long long int state = storeAs( "SameFile" );
  const DbSnapshot snap = snapshot();
  BOOST_CHECK_EQUAL( snap.states.at(state).fg, snap.states.at(state).bg );

  show( "other.n42", {1} );
  load( state );
  BOOST_CHECK( (m_interspec->displayedSamples( SpecUtils::SpectrumType::Foreground ) == set<int>{1, 2}) );
  BOOST_CHECK( (m_interspec->displayedSamples( SpecUtils::SpectrumType::Background ) == set<int>{3}) );
  BOOST_CHECK( m_interspec->measurment( SpecUtils::SpectrumType::Background ) == foreground() );

  check_invariants( snapshot() );
}//SameFileForeAndBackRoundTrip


// Changing which samples are displayed must not disconnect from the state; "Store" then updates
//  the same database entry, recording the newly displayed samples.
BOOST_FIXTURE_TEST_CASE( ChangeSamplesThenStoreUpdatesSameEntry, StateFixture )
{
  show( "multi.n42", {1} );
  const long long int state = storeAs( "Multi" );
  const int fgRow = snapshot().states.at(state).fg;

  m_interspec->changeDisplayedSampleNums( {2}, SpecUtils::SpectrumType::Foreground );
  BOOST_CHECK_EQUAL( link(), state );

  m_interspec->stateSave();

  const DbSnapshot snap = snapshot();
  BOOST_CHECK_EQUAL( snap.topLevelUserStates().size(), 1 );
  BOOST_CHECK_EQUAL( snap.states.at(state).fg, fgRow );
  BOOST_CHECK_EQUAL( snap.states.at(state).fgSamples, "2" );
  check_invariants( snap );

  show( "other.n42", {1} );
  load( state );
  BOOST_CHECK( (m_interspec->displayedSamples( SpecUtils::SpectrumType::Foreground ) == set<int>{2}) );
}//ChangeSamplesThenStoreUpdatesSameEntry


// "Store As..." then "Store" (after changing samples back) must update the new state, and leave
//  the original state, and its rows, alone.
BOOST_FIXTURE_TEST_CASE( StoreAsThenStore, StateFixture )
{
  show( "multi.n42", {1} );
  const long long int stateA = storeAs( "A" );
  const StateRow beforeA = snapshot().states.at(stateA);
  const vector<unsigned char> beforeAData = snapshot().files.at(beforeA.fg).data;

  m_interspec->changeDisplayedSampleNums( {2}, SpecUtils::SpectrumType::Foreground );
  const long long int stateB = storeAs( "B" );
  BOOST_CHECK( stateB != stateA );
  BOOST_CHECK_EQUAL( link(), stateB );

  m_interspec->changeDisplayedSampleNums( {1}, SpecUtils::SpectrumType::Foreground );
  addPeak( 554.3 );
  m_interspec->stateSave();

  const DbSnapshot snap = snapshot();
  const StateRow &a = snap.states.at(stateA);
  const StateRow &b = snap.states.at(stateB);
  BOOST_CHECK_EQUAL( b.fgSamples, "1" );
  BOOST_CHECK_EQUAL( a.fgSamples, beforeA.fgSamples );
  BOOST_CHECK_EQUAL( a.fg, beforeA.fg );
  BOOST_CHECK( a.fg != b.fg );
  BOOST_CHECK( snap.files.at(a.fg).data == beforeAData );
  check_invariants( snap );
}//StoreAsThenStore


// Loading a tag must not change its parent, and storing afterwards updates the parent, not the tag.
BOOST_FIXTURE_TEST_CASE( TagIsImmutable, StateFixture )
{
  show( "fore.n42", {1} );
  const long long int state = storeAs( "Parent" );
  BOOST_REQUIRE_EQUAL( snapshot().states.at(state).bg, -1 );

  show( "back.n42", {1}, SpecUtils::SpectrumType::Background );
  addPeak( 554.3 );
  tag( "Tag" );

  DbSnapshot snap = snapshot();
  const vector<StateRow> tags = snap.children( state, UserState::kUserSaved );
  BOOST_REQUIRE_EQUAL( tags.size(), 1 );
  const StateRow tagRow = tags.front();
  BOOST_REQUIRE( tagRow.bg >= 0 );
  const vector<unsigned char> tagFgData = snap.files.at(tagRow.fg).data;
  const vector<unsigned char> tagBgData = snap.files.at(tagRow.bg).data;
  BOOST_CHECK_EQUAL( snap.states.at(state).bg, -1 ); // The parent wasnt changed by tagging
  check_invariants( snap );

  show( "other.n42", {1} );
  load( tagRow.id );
  BOOST_CHECK_EQUAL( link(), state );
  BOOST_CHECK_EQUAL( npeaks(), 1 );
  BOOST_CHECK_EQUAL( snapshot().states.at(state).bg, -1 ); // Loading doesnt change the parent

  addPeak( 776.5 );
  m_interspec->stateSave();

  snap = snapshot();
  const StateRow &parent = snap.states.at(state);
  const StateRow &tagAfter = snap.states.at(tagRow.id);
  BOOST_CHECK( parent.bg >= 0 );
  BOOST_CHECK( parent.bg != tagAfter.bg );
  BOOST_CHECK( parent.fg != tagAfter.fg );
  BOOST_CHECK_EQUAL( tagAfter.fg, tagRow.fg );
  BOOST_CHECK_EQUAL( tagAfter.bg, tagRow.bg );
  BOOST_CHECK( snap.files.at(tagAfter.fg).data == tagFgData );
  BOOST_CHECK( snap.files.at(tagAfter.bg).data == tagBgData );
  check_invariants( snap );
}//TagIsImmutable


// Auto-saving work when the foreground changes: when connected, a single auto-saved state is kept,
//  and the stored state is untouched; when not connected (but the files came from a state), a
//  standalone spectrum entry is written, and the state is again untouched.
BOOST_FIXTURE_TEST_CASE( AutoSaveOnForegroundChange, StateFixture )
{
  show( "fore.n42", {1} );
  const long long int state = storeAs( "Stored" );
  const DbSnapshot stored = snapshot();
  const vector<unsigned char> storedData = stored.files.at( stored.states.at(state).fg ).data;

  addPeak( 554.3 );
  m_interspec->saveStateAtForegroundChange( false );

  DbSnapshot snap = snapshot();
  vector<StateRow> autoSaved = snap.children( state, UserState::kUserStateAutoSavedWork );
  BOOST_REQUIRE_EQUAL( autoSaved.size(), 1 );
  const long long int firstAutoSaved = autoSaved.front().id;
  BOOST_CHECK( autoSaved.front().fg != snap.states.at(state).fg );
  BOOST_CHECK( snap.files.at( snap.states.at(state).fg ).data == storedData );
  check_invariants( snap );

  addPeak( 776.5 );
  m_interspec->saveStateAtForegroundChange( false );

  snap = snapshot();
  autoSaved = snap.children( state, UserState::kUserStateAutoSavedWork );
  BOOST_REQUIRE_EQUAL( autoSaved.size(), 1 );
  BOOST_CHECK( autoSaved.front().id != firstAutoSaved );
  BOOST_CHECK( snap.files.at( snap.states.at(state).fg ).data == storedData );
  check_invariants( snap );

  // Now load the state (so the file headers refer to its rows), and disconnect from it.
  show( "other.n42", {1} );
  load( state );
  m_interspec->resetCurrentAppStateDbId();
  const string uuid = foreground()->uuid();
  const size_t nstandalone = snapshot().numStandaloneFilesWithUuid( uuid );

  addPeak( 554.3 );
  m_interspec->saveStateAtForegroundChange( false );

  snap = snapshot();
  BOOST_CHECK( snap.files.at( snap.states.at(state).fg ).data == storedData );
  BOOST_CHECK_EQUAL( snap.numStandaloneFilesWithUuid( uuid ), nstandalone + 1 );
  check_invariants( snap );
}//AutoSaveOnForegroundChange


// The generic spectrum saving (used when files are removed from memory, etc.) must never write
//  to, or delete, a database entry that belongs to a state.
BOOST_FIXTURE_TEST_CASE( HeaderSaveNeverAltersState, StateFixture )
{
  show( "fore.n42", {1} );
  const long long int state = storeAs( "Stored" );
  show( "other.n42", {1} );
  load( state );

  const shared_ptr<SpecMeas> fg = foreground();
  const shared_ptr<SpectraFileHeader> header = m_manager->model()->fileHeader( fg );
  BOOST_REQUIRE( header );

  DbSnapshot snap = snapshot();
  const int stateRow = snap.states.at(state).fg;
  const vector<unsigned char> storedData = snap.files.at(stateRow).data;
  BOOST_REQUIRE( header->dbEntry() );
  BOOST_CHECK_EQUAL( header->dbEntry().id(), stateRow );

  // No changes since loading: nothing to write
  const size_t nfiles = snap.files.size();
  header->saveToDatabase( fg );
  BOOST_CHECK_EQUAL( snapshot().files.size(), nfiles );

  addPeak( 554.3 );
  header->saveToDatabase( fg );

  snap = snapshot();
  BOOST_CHECK( snap.files.at(stateRow).data == storedData );
  BOOST_CHECK_EQUAL( snap.files.size(), nfiles + 1 );
  BOOST_REQUIRE( header->dbEntry() );
  const long long int standalone = header->dbEntry().id();
  BOOST_CHECK( standalone != stateRow );
  BOOST_CHECK( !snap.files.at(standalone).partOfState );

  // A subsequent save updates the same standalone entry
  addPeak( 776.5 );
  header->saveToDatabase( fg );
  snap = snapshot();
  BOOST_CHECK_EQUAL( snap.files.size(), nfiles + 1 );
  BOOST_CHECK_EQUAL( header->dbEntry().id(), standalone );
  BOOST_CHECK( snap.files.at(stateRow).data == storedData );

  // Re-pointing a header that refers to a state row must not delete that row
  load( state );
  const shared_ptr<SpectraFileHeader> header2 = m_manager->model()->fileHeader( foreground() );
  BOOST_REQUIRE( header2 && header2->dbEntry() );
  BOOST_REQUIRE_EQUAL( header2->dbEntry().id(), stateRow );
  {
    std::shared_ptr<DataBaseUtils::DbSession> sql = m_interspec->sql();
    Wt::Dbo::ptr<UserFileInDb> other;
    {
      DataBaseUtils::DbTransaction transaction( *sql );
      other = sql->session()->find<UserFileInDb>().where( "id = ?" ).bind( standalone ).resultValue();
      transaction.commit();
    }
    BOOST_REQUIRE( other );
    header2->setDbEntry( other );
  }
  snap = snapshot();
  BOOST_CHECK( snap.files.count(stateRow) );
  check_invariants( snap );
}//HeaderSaveNeverAltersState


// After an explicit "Store", there is no unsaved work, so changing foreground shouldnt create an
//  auto-saved state; but changing samples is a change worth auto-saving.
BOOST_FIXTURE_TEST_CASE( StoreThenForegroundChangeMakesNoAutoSave, StateFixture )
{
  show( "fore.n42", {1} );
  const long long int state = storeAs( "Stored" );
  addPeak( 554.3 );
  BOOST_CHECK( foreground()->modified() );
  m_interspec->stateSave();
  BOOST_CHECK( !foreground()->modified() );

  m_interspec->saveStateAtForegroundChange( false );
  BOOST_CHECK( snapshot().children( state, UserState::kUserStateAutoSavedWork ).empty() );

  m_interspec->changeDisplayedSampleNums( {2}, SpecUtils::SpectrumType::Foreground );
  m_interspec->saveStateAtForegroundChange( false );
  BOOST_CHECK_EQUAL( snapshot().children( state, UserState::kUserStateAutoSavedWork ).size(), 1 );

  check_invariants( snapshot() );
}//StoreThenForegroundChangeMakesNoAutoSave


// Loading the auto-saved work of the connected state, while there are further unsaved changes,
//  must not delete the state being loaded.
BOOST_FIXTURE_TEST_CASE( LoadAutoSavedWorkWhileModified, StateFixture )
{
  show( "fore.n42", {1} );
  addPeak( 554.3 );
  const long long int state = storeAs( "Stored" );
  addPeak( 776.5 );
  m_interspec->saveStateAtForegroundChange( false );

  const vector<StateRow> autoSaved = snapshot().children( state, UserState::kUserStateAutoSavedWork );
  BOOST_REQUIRE_EQUAL( autoSaved.size(), 1 );

  addPeak( 619.1 );
  BOOST_REQUIRE_EQUAL( npeaks(), 3 );
  load( autoSaved.front().id );
  BOOST_CHECK_EQUAL( npeaks(), 2 );
  BOOST_CHECK_EQUAL( link(), state );

  // The just-discarded work (3 peaks) was itself auto-saved, alongside the one we loaded
  const DbSnapshot snap = snapshot();
  BOOST_CHECK_EQUAL( snap.children( state, UserState::kUserStateAutoSavedWork ).size(), 2 );
  check_invariants( snap, true );
}//LoadAutoSavedWorkWhileModified


// End-of-session states: without a connected state, loading it leaves us unconnected; with one,
//  it is a child of that state, and loading it re-connects to the state.
BOOST_FIXTURE_TEST_CASE( EndOfSessionStates, StateFixture )
{
  show( "unconnected.n42", {1} );
  addPeak( 554.3 );
  m_interspec->saveStateForEndOfSession();

  DbSnapshot snap = snapshot();
  vector<StateRow> eos = snap.children( -1, UserState::kEndOfSessionTemp );
  BOOST_REQUIRE_EQUAL( eos.size(), 1 );
  check_invariants( snap );

  show( "other.n42", {1} );
  load( eos.front().id );
  BOOST_CHECK_EQUAL( link(), -1 );
  BOOST_CHECK_EQUAL( npeaks(), 1 );

  show( "fore.n42", {1} );
  const long long int state = storeAs( "Stored" );
  const vector<unsigned char> storedData = snapshot().files.at( snapshot().states.at(state).fg ).data;
  addPeak( 554.3 );
  addPeak( 776.5 );
  m_interspec->saveStateForEndOfSession();

  snap = snapshot();
  BOOST_CHECK( snap.children( -1, UserState::kEndOfSessionTemp ).empty() ); // previous one removed
  eos = snap.children( state, UserState::kEndOfSessionTemp );
  BOOST_REQUIRE_EQUAL( eos.size(), 1 );
  check_invariants( snap );

  show( "other.n42", {1} );
  load( eos.front().id );
  BOOST_CHECK_EQUAL( link(), state );
  BOOST_CHECK_EQUAL( npeaks(), 2 );
  BOOST_CHECK( snapshot().files.at( snapshot().states.at(state).fg ).data == storedData );

  m_interspec->saveStateForEndOfSession();
  snap = snapshot();
  BOOST_CHECK_EQUAL( snap.children( state, UserState::kEndOfSessionTemp ).size(), 1 );
  check_invariants( snap );
}//EndOfSessionStates


// If the connected state is deleted (e.g., from another window), saving work must still work, and
//  "Store As..." can create a new state.
BOOST_FIXTURE_TEST_CASE( DeletedConnectedState, StateFixture )
{
  show( "fore.n42", {1} );
  const long long int state = storeAs( "Stored" );
  tag( "Tag" );
  show( "other.n42", {1} );
  load( state );

  {
    std::shared_ptr<DataBaseUtils::DbSession> sql = m_interspec->sql();
    BOOST_REQUIRE_NO_THROW( UserState::removeFromDatabase( stateById( state ), *sql ) );
  }
  BOOST_CHECK( snapshot().states.empty() );

  const string uuid = foreground()->uuid();
  addPeak( 554.3 );
  BOOST_CHECK_NO_THROW( m_interspec->saveStateAtForegroundChange( false ) );

  DbSnapshot snap = snapshot();
  BOOST_CHECK( snap.states.empty() );
  BOOST_CHECK_EQUAL( snap.numStandaloneFilesWithUuid( uuid ), 1 );

  m_interspec->stateSave(); // Not connected to anything, so should do nothing besides show a dialog
  BOOST_CHECK( snapshot().states.empty() );

  storeAs( "New" );
  snap = snapshot();
  BOOST_CHECK_EQUAL( snap.topLevelUserStates().size(), 1 );
  check_invariants( snap );
}//DeletedConnectedState


// Historical bugs could leave two states sharing a spectrum row; storing one of them must not
//  modify the other.
BOOST_FIXTURE_TEST_CASE( LegacySharedRowIsSplitOnStore, StateFixture )
{
  show( "fore.n42", {1} );
  const long long int state1 = storeAs( "One" );
  m_interspec->changeDisplayedSampleNums( {2}, SpecUtils::SpectrumType::Foreground );
  const long long int state2 = storeAs( "Two" );

  DbSnapshot snap = snapshot();
  const int sharedRow = snap.states.at(state1).fg;
  const int ownRow2 = snap.states.at(state2).fg;
  const vector<unsigned char> data1 = snap.files.at(sharedRow).data;

  {// Simulate the legacy corruption, through the apps session (states use optimistic locking)
    std::shared_ptr<DataBaseUtils::DbSession> sql = m_interspec->sql();
    DataBaseUtils::DbTransaction transaction( *sql );
    Wt::Dbo::ptr<UserState> st2 = sql->session()->find<UserState>().where( "id = ?" ).bind( state2 ).resultValue();
    st2.modify()->foregroundId = sharedRow;
    Wt::Dbo::ptr<UserFileInDb> own = sql->session()->find<UserFileInDb>().where( "id = ?" ).bind( ownRow2 ).resultValue();
    own.remove();
    transaction.commit();
  }

  BOOST_REQUIRE_EQUAL( link(), state2 );
  addPeak( 554.3 );
  m_interspec->stateSave();

  snap = snapshot();
  BOOST_CHECK_EQUAL( snap.states.at(state1).fg, sharedRow );
  BOOST_CHECK( snap.states.at(state2).fg != sharedRow );
  BOOST_CHECK( snap.files.at(sharedRow).data == data1 );
  check_invariants( snap );
}//LegacySharedRowIsSplitOnStore


// Switching the foreground to another file and back reconnects to the state; and the connection
//  survives the file being removed from memory and re-read.
BOOST_FIXTURE_TEST_CASE( SwitchAwayAndBackReconnects, StateFixture )
{
  const shared_ptr<SpecMeas> fore = show( "fore.n42", {1} );
  const long long int state = storeAs( "Stored" );

  show( "other.n42", {1} );
  BOOST_CHECK_EQUAL( link(), -1 );

  m_interspec->setSpectrum( fore, {2}, SpecUtils::SpectrumType::Foreground, {} );
  BOOST_CHECK_EQUAL( link(), state );

  // A file written to disk when removed from memory, then parsed back in, keeps its connection
  const shared_ptr<SpectraFileHeader> header = make_shared<SpectraFileHeader>( false, m_interspec );
  shared_ptr<SpecMeas> meas = make_file( "evicted.n42" );
  header->setFile( meas->filename(), meas );
  meas->setDbStateId( 42 );
  const weak_ptr<SpecMeas> weak = meas;
  meas.reset();
  BOOST_REQUIRE( weak.expired() );
  const shared_ptr<SpecMeas> reread = header->parseFile();
  BOOST_REQUIRE( reread );
  BOOST_CHECK_EQUAL( reread->dbStateId(), 42 );

  // But a copy of a file is a different file, so isnt connected
  SpecMeas copy;
  copy.uniqueCopyContents( *reread );
  BOOST_CHECK_EQUAL( copy.dbStateId(), -1 );
}//SwitchAwayAndBackReconnects


// Open a state from the database, switch the foreground to a different file, then back using the
//  "Spectrum Files" tab: we should be re-connected to the state (with the samples we were last
//  looking at), and "Store" should update it.  Also when the file was removed from memory meanwhile.
BOOST_FIXTURE_TEST_CASE( ReconnectAfterSwitchingViaFileManager, StateFixture )
{
  show( "fore.n42", {1, 2} );
  const long long int state = storeAs( "Stored" );
  show( "temp.n42", {1} );
  load( state ); // the foreground now comes from the database
  BOOST_REQUIRE_EQUAL( link(), state );

  shared_ptr<SpecMeas> fromState = foreground();
  const shared_ptr<SpectraFileHeader> stateHeader = m_manager->model()->fileHeader( fromState );
  BOOST_REQUIRE( stateHeader );

  // Some unsaved work
  m_interspec->changeDisplayedSampleNums( {2}, SpecUtils::SpectrumType::Foreground );
  addPeak( 554.3 );

  // Switch to a different file.  (In the app, `setSpectrum` first auto-saves the outgoing work; the
  //  unit-test build skips that, so we do it here.)
  const shared_ptr<SpecMeas> other = make_file( "other.n42" );
  const shared_ptr<SpectraFileHeader> otherHeader = m_manager->addFile( other->filename(), other );
  m_interspec->saveStateAtForegroundChange( false );
  displayViaFileManager( otherHeader );
  BOOST_REQUIRE( foreground() == other );
  BOOST_CHECK_EQUAL( link(), -1 );
  BOOST_CHECK_EQUAL( snapshot().children( state, UserState::kUserStateAutoSavedWork ).size(), 1 );

  // And back
  displayViaFileManager( stateHeader );
  BOOST_REQUIRE( foreground() == fromState );
  BOOST_CHECK_EQUAL( link(), state );
  BOOST_CHECK( (m_interspec->displayedSamples( SpecUtils::SpectrumType::Foreground ) == set<int>{2}) );
  BOOST_CHECK_EQUAL( npeaks(), 1 );

  m_interspec->stateSave();
  DbSnapshot snap = snapshot();
  BOOST_CHECK_EQUAL( snap.topLevelUserStates().size(), 1 );
  BOOST_CHECK_EQUAL( snap.states.at(state).fgSamples, "2" );
  BOOST_CHECK( snap.children( state, UserState::kUserStateAutoSavedWork ).empty() );
  check_invariants( snap );

  // Now switch away again, and have the file removed from memory (written to a temporary file),
  //  before switching back.
  displayViaFileManager( otherHeader );
  BOOST_REQUIRE( foreground() == other );
  m_manager->removeFromSpectrumInfoCache( fromState, false );
  const weak_ptr<SpecMeas> weakFromState = fromState;
  fromState.reset();
  BOOST_REQUIRE_MESSAGE( weakFromState.expired(), "File from state is still referenced elsewhere" );

  displayViaFileManager( stateHeader );
  BOOST_CHECK( foreground() && (foreground() != other) );
  BOOST_CHECK_EQUAL( link(), state );
  BOOST_CHECK( (m_interspec->displayedSamples( SpecUtils::SpectrumType::Foreground ) == set<int>{2}) );
  BOOST_CHECK_EQUAL( npeaks(), 1 );
}//ReconnectAfterSwitchingViaFileManager


// A failed save (here, a file too large to store in the database) must not leave anything behind -
//  not even once the apps database session is used again - nor change the stored state.
BOOST_FIXTURE_TEST_CASE( FailedSavesLeaveNothingBehind, StateFixture )
{
  const shared_ptr<SpecMeas> big = make_large_file( "big.n42" );
  m_manager->addFile( big->filename(), big );
  m_interspec->setSpectrum( big, {1}, SpecUtils::SpectrumType::Foreground, {} );
  BOOST_REQUIRE( foreground() == big );

  m_interspec->finishStoreStateInDb( "TooBig", "", Wt::Dbo::ptr<UserState>() );
  BOOST_CHECK_EQUAL( link(), -1 );
  touchAppSession();
  DbSnapshot snap = snapshot();
  BOOST_CHECK( snap.states.empty() );
  check_invariants( snap );

  // Now a stored state, whose auto-saving (and "Store") fails because of a too-large background
  show( "fore.n42", {1} );
  const long long int state = storeAs( "Small" );
  const DbSnapshot stored = snapshot();
  m_interspec->setSpectrum( big, {1}, SpecUtils::SpectrumType::Background, {} );
  BOOST_REQUIRE( m_interspec->measurment( SpecUtils::SpectrumType::Background ) == big );
  addPeak( 554.3 );

  m_interspec->saveStateAtForegroundChange( false );
  touchAppSession();
  snap = snapshot();
  BOOST_CHECK( snap.children( state, UserState::kUserStateAutoSavedWork ).empty() );
  BOOST_CHECK( foreground()->modified() ); // The work is still unsaved
  check_same_states( stored, snap );

  m_interspec->stateSave();
  touchAppSession();
  snap = snapshot();
  check_same_states( stored, snap );
  BOOST_CHECK_EQUAL( link(), state );
  check_invariants( snap );
}//FailedSavesLeaveNothingBehind


// Each time an (unconnected) end-of-session state is resumed and worked on, the files work is saved
//  to a standalone database entry; this should reuse the files entry, not add a new one each session.
BOOST_FIXTURE_TEST_CASE( ResumingEndOfSessionDoesntAccumulateRows, StateFixture )
{
  show( "resumed.n42", {1} );
  const string uuid = foreground()->uuid();

  for( size_t session = 0; session < 3; ++session )
  {
    addPeak( 500.0 + 50.0*session );

    // Like the end of a session: auto-save the work, then save the end-of-session state
    m_interspec->saveStateAtForegroundChange( false );
    m_interspec->saveStateForEndOfSession();

    // Like starting the next session: resume the end-of-session state
    const vector<StateRow> eos = snapshot().children( -1, UserState::kEndOfSessionTemp );
    BOOST_REQUIRE_EQUAL( eos.size(), 1 );
    load( eos.front().id );
    BOOST_CHECK_EQUAL( link(), -1 );
    BOOST_CHECK_EQUAL( npeaks(), session + 1 );
  }//for( size_t session = 0; session < 3; ++session )

  const DbSnapshot snap = snapshot();
  BOOST_CHECK_EQUAL( snap.numStandaloneFilesWithUuid( uuid ), 1 );
  check_invariants( snap );
}//ResumingEndOfSessionDoesntAccumulateRows


// The generic file saving (e.g., when a file is removed from memory) must not hide unsaved work
//  from the connected state, so it still gets auto-saved.
BOOST_FIXTURE_TEST_CASE( GenericSaveKeepsWorkUnsavedForState, StateFixture )
{
  show( "fore.n42", {1} );
  const long long int state = storeAs( "Stored" );
  addPeak( 554.3 );

  const shared_ptr<SpectraFileHeader> header = m_manager->model()->fileHeader( foreground() );
  BOOST_REQUIRE( header );
  header->saveToDatabase( foreground() );
  BOOST_CHECK( foreground()->modified() );

  m_interspec->saveStateAtForegroundChange( false );
  BOOST_CHECK_EQUAL( snapshot().children( state, UserState::kUserStateAutoSavedWork ).size(), 1 );
  check_invariants( snapshot() );
}//GenericSaveKeepsWorkUnsavedForState


// Legacy databases may have rows whose `SnapshotParent` is a files standalone entry; deleting that
//  entry would cascade-delete them, so it must be treated as owned.
BOOST_FIXTURE_TEST_CASE( LegacySnapshotParentRowsAreProtected, StateFixture )
{
  const shared_ptr<SpecMeas> fore = show( "fore.n42", {1} );
  const shared_ptr<SpectraFileHeader> header = m_manager->model()->fileHeader( fore );
  BOOST_REQUIRE( header );
  header->saveToDatabase( fore ); // creates the files standalone entry
  BOOST_REQUIRE( header->dbEntry() );
  const long long int parentId = header->dbEntry().id();

  std::shared_ptr<DataBaseUtils::DbSession> sql = m_interspec->sql();
  long long int childId = -1;
  {
    DataBaseUtils::DbTransaction transaction( *sql );
    std::unique_ptr<UserFileInDb> child = std::make_unique<UserFileInDb>();
    child->user = m_interspec->user();
    child->uuid = "legacy-child";
    child->snapshotParent = header->dbEntry();
    Wt::Dbo::ptr<UserFileInDb> childPtr = sql->session()->add( std::move(child) );
    sql->session()->flush();
    childId = childPtr.id();
    BOOST_CHECK( UserFileInDb::isOwnedByState( header->dbEntry(), *sql->session() ) );
    transaction.commit();
  }

  // Re-pointing the header at a different entry must not delete the parent row
  std::unique_ptr<UserFileInDb> replacement = std::make_unique<UserFileInDb>();
  replacement->user = m_interspec->user();
  header->setDbEntry( Wt::Dbo::ptr<UserFileInDb>( std::move(replacement) ) );

  const DbSnapshot snap = snapshot();
  BOOST_CHECK( snap.files.count( parentId ) );
  BOOST_CHECK( snap.files.count( childId ) );
}//LegacySnapshotParentRowsAreProtected


// After the auto-saved state we loaded from is removed (by storing), the file header that referred to
//  its rows must not hand out the deleted row.
BOOST_FIXTURE_TEST_CASE( DeletedRowsArentHandedOut, StateFixture )
{
  show( "fore.n42", {1} );
  const long long int state = storeAs( "Stored" );
  addPeak( 554.3 );
  m_interspec->saveStateAtForegroundChange( false );
  show( "other.n42", {1} );
  load( state );
  m_interspec->loadNewestAutoSavedWork( state );
  BOOST_REQUIRE_EQUAL( npeaks(), 1 );

  m_interspec->stateSave(); // removes the auto-saved state, and its rows
  BOOST_CHECK( snapshot().children( state, UserState::kUserStateAutoSavedWork ).empty() );

  const Wt::Dbo::ptr<UserFileInDb> fromDb
                 = m_interspec->measurementFromDb( SpecUtils::SpectrumType::Foreground, false );
  BOOST_CHECK( !fromDb || !fromDb.isTransient() );
  check_invariants( snapshot() );
}//DeletedRowsArentHandedOut


// The background-subtract display option is restored.
BOOST_FIXTURE_TEST_CASE( BackgroundSubtractRoundTrip, StateFixture )
{
  show( "fore.n42", {1} );
  show( "back.n42", {1}, SpecUtils::SpectrumType::Background );
  m_interspec->setBackgroundSub( true );
  const long long int state1 = storeAs( "WithSub" );
  BOOST_CHECK_EQUAL( snapshot().states.at(state1).bgSubMode, UserState::kBackgorundSubtract );

  m_interspec->setBackgroundSub( false );
  show( "other.n42", {1} );
  load( state1 );
  const long long int state2 = storeAs( "Resaved" );
  BOOST_CHECK_EQUAL( snapshot().states.at(state2).bgSubMode, UserState::kBackgorundSubtract );

  check_invariants( snapshot() );
}//BackgroundSubtractRoundTrip


// Loading a stored state that has newer auto-saved work: the link in the notification loads the work.
BOOST_FIXTURE_TEST_CASE( LoadNewestAutoSavedWork, StateFixture )
{
  show( "fore.n42", {1} );
  const long long int state = storeAs( "Stored" );
  addPeak( 554.3 );
  m_interspec->saveStateAtForegroundChange( false );
  BOOST_REQUIRE_EQUAL( snapshot().children( state, UserState::kUserStateAutoSavedWork ).size(), 1 );

  show( "other.n42", {1} );
  load( state );
  BOOST_CHECK_EQUAL( npeaks(), 0 );

  m_interspec->loadNewestAutoSavedWork( state );
  BOOST_CHECK_EQUAL( npeaks(), 1 );
  BOOST_CHECK_EQUAL( link(), state );

  check_invariants( snapshot() );
}//LoadNewestAutoSavedWork

#else

BOOST_AUTO_TEST_CASE( NoDatabase )
{
  BOOST_TEST_MESSAGE( "USE_DB_TO_STORE_SPECTRA is off; nothing to test." );
}

#endif //#if( USE_DB_TO_STORE_SPECTRA )
