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
#include <deque>
#include <string>
#include <memory>
#include <vector>

#include <Wt/WApplication.h>
#include <Wt/Test/WTestEnvironment.h>

#ifdef _WIN32
#include "winsock2.h"
#include "Windows.h"
#endif

#define BOOST_TEST_MODULE EnergyCalTool_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/EnergyCalTool.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace std;
using namespace boost::unit_test;

typedef EnergyCalTool::ApplyToCbIndex CbIndex;
typedef EnergyCalTool::ApplyToSelection Selection;
typedef EnergyCalTool::ApplyToDisplayState DisplayState;


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
  }
  SpecUtils::ireplace_all( datadir, "%20", " " );

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

  if( datadir.empty() )
    throw runtime_error( "Could not find data directory" );

  DecayDataBaseServer::setDecayXmlFile( SpecUtils::append_path(datadir, "sandia.decay.xml") );
  InterSpec::setStaticDataDirectory( datadir );
}//void set_data_dir()


namespace
{
const float sm_peak_energy = 661.7f;

/** Creates a SpecMeas with the given samples and detectors (1024 channels).  If `per_det_cals`, each
 detector gets its own calibration (shared by all its samples), with detector `i` having gain
 `3 + 0.01*i` keV/channel and offset `0.5*i` keV - so the last detector has the largest energy range,
 and so defines the calibration of the summed (displayed) spectrum; otherwise every measurement
 shares a single {0, 3 keV/chan} calibration.  Measurements listed in `missing` (sample, detector)
 are not created.
 */
shared_ptr<SpecMeas> make_file( const vector<int> &samples, const vector<string> &dets,
                                const bool per_det_cals,
                                const set<pair<int,string>> &missing = {} )
{
  vector<shared_ptr<const SpecUtils::EnergyCalibration>> cals;
  for( size_t i = 0; i < dets.size(); ++i )
  {
    if( !per_det_cals && !cals.empty() )
    {
      cals.push_back( cals.front() );
      continue;
    }

    auto cal = make_shared<SpecUtils::EnergyCalibration>();
    const float offset = per_det_cals ? 0.5f*i : 0.0f;
    const float gain = per_det_cals ? (3.0f + 0.01f*i) : 3.0f;
    cal->set_polynomial( 1024, {offset, gain}, {} );
    cals.push_back( cal );
  }//for( size_t i = 0; i < dets.size(); ++i )

  auto meas = make_shared<SpecMeas>();
  for( const int sample : samples )
  {
    for( size_t i = 0; i < dets.size(); ++i )
    {
      if( missing.count( {sample, dets[i]} ) )
        continue;

      // Flat continuum plus a Gaussian at sm_peak_energy, sigma 2 channels
      auto counts = make_shared<vector<float>>( 1024, 10.0f );
      const double peak_chan = cals[i]->channel_for_energy( sm_peak_energy );
      for( int ch = 0; ch < 1024; ++ch )
      {
        const double x = (ch + 0.5) - peak_chan;
        if( fabs(x) < 20 )
          (*counts)[ch] += static_cast<float>( 1000.0*exp( -0.5*x*x/4.0 ) );
      }

      auto m = make_shared<SpecUtils::Measurement>();
      m->set_gamma_counts( counts, 300.0f, 315.0f );
      m->set_energy_calibration( cals[i] );
      m->set_sample_number( sample );
      m->set_detector_name( dets[i] );
      m->set_source_type( SpecUtils::SourceType::Foreground );
      meas->add_measurement( m, false );
    }//for( loop over detectors )
  }//for( const int sample : samples )

  meas->cleanup_after_load( SpecUtils::SpecFile::DontChangeOrReorderSamples );
  return meas;
}//make_file(...)


/** A single-peak deque, with the peak at `mean`. */
deque<shared_ptr<const PeakDef>> make_peaks( const double mean )
{
  auto cont = make_shared<PeakContinuum>();
  cont->setType( PeakContinuum::OffsetType::Linear );
  cont->setRange( mean - 30.0, mean + 30.0 );

  auto peak = make_shared<PeakDef>( mean, 6.0, 1000.0 );
  peak->setContinuum( cont );

  return deque<shared_ptr<const PeakDef>>{ peak };
}//make_peaks(...)


DisplayState display_state( const shared_ptr<SpecMeas> &fore, const set<int> &fore_samples,
                            const shared_ptr<SpecMeas> &back = nullptr, const set<int> &back_samples = {},
                            const shared_ptr<SpecMeas> &sec = nullptr, const set<int> &sec_samples = {} )
{
  DisplayState state;
  const shared_ptr<SpecMeas> files[3] = { fore, back, sec };
  const set<int> samples[3] = { fore_samples, back_samples, sec_samples };
  for( int i = 0; i < 3; ++i )
  {
    state.meas[i] = files[i];
    if( files[i] )
    {
      state.samples[i] = samples[i];
      state.detectors[i] = files[i]->detector_names();
    }
  }
  return state;
}//display_state(...)


const MeasToApplyCoefChangeTo *entry_for( const vector<MeasToApplyCoefChangeTo> &scope,
                                          const shared_ptr<SpecMeas> &meas )
{
  for( const MeasToApplyCoefChangeTo &c : scope )
  {
    if( c.meas == meas )
      return &c;
  }
  return nullptr;
}//entry_for(...)


set<string> dets( std::initializer_list<string> names )
{
  return set<string>( names );
}

set<int> samps( std::initializer_list<int> nums )
{
  return set<int>( nums );
}


shared_ptr<const SpecUtils::EnergyCalibration> cal_of( const shared_ptr<SpecMeas> &meas,
                                                        const int sample, const string &det )
{
  const shared_ptr<const SpecUtils::Measurement> m = meas->measurement( sample, det );
  BOOST_REQUIRE( m );
  return m->energy_calibration();
}//cal_of(...)


/** A copy of `cal` with its gain scaled by `scale`. */
shared_ptr<const SpecUtils::EnergyCalibration> scaled_gain(
                                     const shared_ptr<const SpecUtils::EnergyCalibration> &cal,
                                     const float scale )
{
  vector<float> coefs = cal->coefficients();
  coefs[1] *= scale;
  auto answer = make_shared<SpecUtils::EnergyCalibration>();
  answer->set_polynomial( cal->num_channels(), coefs, cal->deviation_pairs() );
  return answer;
}//scaled_gain(...)


double peak_mean( const shared_ptr<SpecMeas> &meas, const set<int> &samples )
{
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks
                                          = static_cast<const SpecMeas &>(*meas).peaks( samples );
  BOOST_REQUIRE( peaks && (peaks->size() == 1) );
  return peaks->front()->mean();
}//peak_mean(...)
}//namespace


/** Which "Apply Changes To" rows are applicable, over a number of display configurations. */
BOOST_AUTO_TEST_CASE( applicableRows )
{
  set_data_dir();

  // No foreground: nothing applicable
  {
    const array<bool,CbIndex::NumApplyToCbIndex> rows = EnergyCalTool::applicable_apply_to_rows( DisplayState{} );
    for( const bool r : rows )
      BOOST_CHECK( !r );
  }

  // Single file, detector, and sample: nothing applicable
  {
    const shared_ptr<SpecMeas> f = make_file( {1}, {"Aa1"}, false );
    const auto rows = EnergyCalTool::applicable_apply_to_rows( display_state( f, {1} ) );
    for( const bool r : rows )
      BOOST_CHECK( !r );
  }

  // Multiple detectors and samples, all detectors displayed, only some samples
  {
    const shared_ptr<SpecMeas> f = make_file( {1,2,3}, {"Aa1","Aa2","Aa3"}, true );
    const auto rows = EnergyCalTool::applicable_apply_to_rows( display_state( f, {2} ) );
    BOOST_CHECK( !rows[CbIndex::ApplyToForeground] );
    BOOST_CHECK( !rows[CbIndex::ApplyToBackground] );
    BOOST_CHECK( !rows[CbIndex::ApplyToDisplayedDetectors] );
    BOOST_CHECK( !rows[CbIndex::ApplyToAllDetectors] );
    BOOST_CHECK( rows[CbIndex::ApplyToSelectedDetOnly] );
    BOOST_CHECK( rows[CbIndex::ApplyToDisplayedSamples] );
    BOOST_CHECK( rows[CbIndex::ApplyToAllSamples] );
  }

  // Only a subset of detectors displayed
  {
    const shared_ptr<SpecMeas> f = make_file( {1}, {"Aa1","Aa2","Aa3"}, true );
    DisplayState state = display_state( f, {1} );
    state.detectors[0] = { "Aa1", "Aa2" };
    const auto rows = EnergyCalTool::applicable_apply_to_rows( state );
    BOOST_CHECK( rows[CbIndex::ApplyToDisplayedDetectors] );
    BOOST_CHECK( rows[CbIndex::ApplyToAllDetectors] );
    BOOST_CHECK( rows[CbIndex::ApplyToSelectedDetOnly] );
    BOOST_CHECK( !rows[CbIndex::ApplyToDisplayedSamples] );
  }

  // Single-detector foreground and background: "Indep. Dets." not applicable
  {
    const shared_ptr<SpecMeas> f = make_file( {1}, {"Aa1"}, false );
    const shared_ptr<SpecMeas> b = make_file( {1}, {"Bb1"}, false );
    const auto rows = EnergyCalTool::applicable_apply_to_rows( display_state( f, {1}, b, {1} ) );
    BOOST_CHECK( rows[CbIndex::ApplyToForeground] );
    BOOST_CHECK( rows[CbIndex::ApplyToBackground] );
    BOOST_CHECK( !rows[CbIndex::ApplyToSecondary] );
    BOOST_CHECK( !rows[CbIndex::ApplyToSelectedDetOnly] );
  }

  // Same file as foreground and background, each displaying some of the samples - all samples
  //  of the file are displayed between the two, so the sample rows dont apply
  {
    const shared_ptr<SpecMeas> f = make_file( {1,2}, {"Aa1"}, false );
    const auto rows = EnergyCalTool::applicable_apply_to_rows( display_state( f, {1}, f, {2} ) );
    BOOST_CHECK( !rows[CbIndex::ApplyToDisplayedSamples] );
    BOOST_CHECK( !rows[CbIndex::ApplyToAllSamples] );
  }
}//BOOST_AUTO_TEST_CASE( applicableRows )


/** The measurements the "Apply Changes To" selection maps to, over many combinations. */
BOOST_AUTO_TEST_CASE( applyToScopeMatrix )
{
  set_data_dir();

  const Selection defaults;

  // Single file/detector/sample: everything
  {
    const shared_ptr<SpecMeas> f = make_file( {1}, {"Aa1"}, false );
    const vector<MeasToApplyCoefChangeTo> scope = EnergyCalTool::apply_to_scope( defaults, display_state( f, {1} ) );
    BOOST_REQUIRE_EQUAL( scope.size(), 1 );
    BOOST_CHECK( scope[0].detectors == dets({"Aa1"}) );
    BOOST_CHECK( scope[0].sample_numbers == samps({1}) );
  }

  const shared_ptr<SpecMeas> multi = make_file( {1,2,3}, {"Aa1","Aa2","Aa3"}, true );

  // Default: all detectors, all samples
  {
    const auto scope = EnergyCalTool::apply_to_scope( defaults, display_state( multi, {2} ) );
    BOOST_REQUIRE_EQUAL( scope.size(), 1 );
    BOOST_CHECK( scope[0].detectors == dets({"Aa1","Aa2","Aa3"}) );
    BOOST_CHECK( scope[0].sample_numbers == samps({1,2,3}) );
  }

  // Displayed samples only
  {
    Selection sel;
    sel.checked[CbIndex::ApplyToDisplayedSamples] = true;
    sel.checked[CbIndex::ApplyToAllSamples] = false;
    const auto scope = EnergyCalTool::apply_to_scope( sel, display_state( multi, {2} ) );
    BOOST_REQUIRE_EQUAL( scope.size(), 1 );
    BOOST_CHECK( scope[0].sample_numbers == samps({2}) );
  }

  // Selected detector only
  {
    Selection sel;
    sel.checked[CbIndex::ApplyToSelectedDetOnly] = true;
    sel.selected_detector = "Aa2";
    const auto scope = EnergyCalTool::apply_to_scope( sel, display_state( multi, {2} ) );
    BOOST_REQUIRE_EQUAL( scope.size(), 1 );
    BOOST_CHECK( scope[0].detectors == dets({"Aa2"}) );
    BOOST_CHECK( scope[0].sample_numbers == samps({1,2,3}) );
  }

  // Subset of detectors displayed: "Displayed Detectors" vs "All Detectors", and "Selected Det.
  //  Only" applies to the selected detector even if it isnt displayed
  {
    DisplayState state = display_state( multi, {1,2,3} );
    state.detectors[0] = { "Aa1", "Aa2" };

    Selection sel;
    sel.checked[CbIndex::ApplyToDisplayedDetectors] = true;
    sel.checked[CbIndex::ApplyToAllDetectors] = false;
    auto scope = EnergyCalTool::apply_to_scope( sel, state );
    BOOST_REQUIRE_EQUAL( scope.size(), 1 );
    BOOST_CHECK( scope[0].detectors == dets({"Aa1","Aa2"}) );

    scope = EnergyCalTool::apply_to_scope( defaults, state );
    BOOST_REQUIRE_EQUAL( scope.size(), 1 );
    BOOST_CHECK( scope[0].detectors == dets({"Aa1","Aa2","Aa3"}) );

    sel.checked[CbIndex::ApplyToSelectedDetOnly] = true;
    sel.selected_detector = "Aa3";
    scope = EnergyCalTool::apply_to_scope( sel, state );
    BOOST_REQUIRE_EQUAL( scope.size(), 1 );
    BOOST_CHECK( scope[0].detectors == dets({"Aa3"}) );

    // An inconsistent (both unchecked) detector choice falls back to all detectors
    Selection both_off;
    both_off.checked[CbIndex::ApplyToDisplayedDetectors] = false;
    both_off.checked[CbIndex::ApplyToAllDetectors] = false;
    scope = EnergyCalTool::apply_to_scope( both_off, state );
    BOOST_REQUIRE_EQUAL( scope.size(), 1 );
    BOOST_CHECK( scope[0].detectors == dets({"Aa1","Aa2","Aa3"}) );
  }

  // A detector that only has data in samples that arent displayed is still included by
  //  "All Detectors" + "All Samples" (previously only the detectors of displayed samples were used)
  {
    const shared_ptr<SpecMeas> f = make_file( {1,2,3,4}, {"Aa1","Ba2"}, true, {{1,"Ba2"},{2,"Ba2"}} );
    const auto scope = EnergyCalTool::apply_to_scope( defaults, display_state( f, {1} ) );
    BOOST_REQUIRE_EQUAL( scope.size(), 1 );
    BOOST_CHECK( scope[0].detectors == dets({"Aa1","Ba2"}) );
    BOOST_CHECK( scope[0].sample_numbers == samps({1,2,3,4}) );
  }

  // Foreground and background different files, with the same detector names
  {
    const shared_ptr<SpecMeas> f = make_file( {1}, {"Aa1","Aa2"}, true );
    const shared_ptr<SpecMeas> b = make_file( {1}, {"Aa1","Aa2"}, true );

    auto scope = EnergyCalTool::apply_to_scope( defaults, display_state( f, {1}, b, {1} ) );
    BOOST_CHECK_EQUAL( scope.size(), 2 );
    BOOST_CHECK( entry_for( scope, f ) && entry_for( scope, b ) );

    Selection no_back;
    no_back.checked[CbIndex::ApplyToBackground] = false;
    scope = EnergyCalTool::apply_to_scope( no_back, display_state( f, {1}, b, {1} ) );
    BOOST_REQUIRE_EQUAL( scope.size(), 1 );
    BOOST_CHECK( scope[0].meas == f );

    Selection no_fore;
    no_fore.checked[CbIndex::ApplyToForeground] = false;
    scope = EnergyCalTool::apply_to_scope( no_fore, display_state( f, {1}, b, {1} ) );
    BOOST_REQUIRE_EQUAL( scope.size(), 1 );
    BOOST_CHECK( scope[0].meas == b );

    // "Indep. Dets." matches the detector by name in both files
    Selection sel;
    sel.checked[CbIndex::ApplyToSelectedDetOnly] = true;
    sel.selected_detector = "Aa2";
    scope = EnergyCalTool::apply_to_scope( sel, display_state( f, {1}, b, {1} ) );
    BOOST_REQUIRE_EQUAL( scope.size(), 2 );
    BOOST_CHECK( entry_for( scope, f )->detectors == dets({"Aa2"}) );
    BOOST_CHECK( entry_for( scope, b )->detectors == dets({"Aa2"}) );
  }

  // Background with different detector names: "Indep. Dets." leaves it alone
  {
    const shared_ptr<SpecMeas> f = make_file( {1}, {"Aa1","Aa2"}, true );
    const shared_ptr<SpecMeas> b = make_file( {1}, {"Bb1","Bb2"}, true );
    Selection sel;
    sel.checked[CbIndex::ApplyToSelectedDetOnly] = true;
    sel.selected_detector = "Aa1";
    const auto scope = EnergyCalTool::apply_to_scope( sel, display_state( f, {1}, b, {1} ) );
    BOOST_REQUIRE_EQUAL( scope.size(), 1 );
    BOOST_CHECK( scope[0].meas == f );
    BOOST_CHECK( scope[0].detectors == dets({"Aa1"}) );
  }

  // Same file as foreground (sample 1) and background (sample 2), of a 4 sample file
  {
    const shared_ptr<SpecMeas> f = make_file( {1,2,3,4}, {"Aa1"}, false );
    const DisplayState state = display_state( f, {1}, f, {2} );

    auto scope = EnergyCalTool::apply_to_scope( defaults, state );
    BOOST_REQUIRE_EQUAL( scope.size(), 1 );  //never duplicate entries for a file
    BOOST_CHECK( scope[0].sample_numbers == samps({1,2,3,4}) );

    Selection no_back;
    no_back.checked[CbIndex::ApplyToBackground] = false;
    scope = EnergyCalTool::apply_to_scope( no_back, state );
    BOOST_REQUIRE_EQUAL( scope.size(), 1 );
    BOOST_CHECK( scope[0].sample_numbers == samps({1,3,4}) );

    Selection no_fore;
    no_fore.checked[CbIndex::ApplyToForeground] = false;
    scope = EnergyCalTool::apply_to_scope( no_fore, state );
    BOOST_REQUIRE_EQUAL( scope.size(), 1 );
    BOOST_CHECK( scope[0].sample_numbers == samps({2,3,4}) );

    Selection disp_only = no_back;
    disp_only.checked[CbIndex::ApplyToDisplayedSamples] = true;
    disp_only.checked[CbIndex::ApplyToAllSamples] = false;
    scope = EnergyCalTool::apply_to_scope( disp_only, state );
    BOOST_REQUIRE_EQUAL( scope.size(), 1 );
    BOOST_CHECK( scope[0].sample_numbers == samps({1}) );
  }

  // Secondary spectrum
  {
    const shared_ptr<SpecMeas> f = make_file( {1}, {"Aa1"}, false );
    const shared_ptr<SpecMeas> s = make_file( {1}, {"Aa1"}, false );
    Selection no_sec;
    no_sec.checked[CbIndex::ApplyToSecondary] = false;
    const auto scope = EnergyCalTool::apply_to_scope( no_sec, display_state( f, {1}, nullptr, {}, s, {1} ) );
    BOOST_REQUIRE_EQUAL( scope.size(), 1 );
    BOOST_CHECK( scope[0].meas == f );
  }

  // A checkbox that isnt applicable (so isnt shown) is ignored: "Foreground" unchecked while there
  //  is no background/secondary (e.g., the background was unloaded after unchecking it), or
  //  "Indep. Dets." checked while there is only a single detector.
  {
    const shared_ptr<SpecMeas> f = make_file( {1}, {"Aa1"}, false );
    Selection sel;
    sel.checked[CbIndex::ApplyToForeground] = false;
    sel.checked[CbIndex::ApplyToSelectedDetOnly] = true;
    sel.selected_detector = "NotADet";
    const auto scope = EnergyCalTool::apply_to_scope( sel, display_state( f, {1} ) );
    BOOST_REQUIRE_EQUAL( scope.size(), 1 );
    BOOST_CHECK( scope[0].detectors == dets({"Aa1"}) );
  }

  // "Indep. Dets." with a detector not in any file: nothing to change
  {
    Selection sel;
    sel.checked[CbIndex::ApplyToSelectedDetOnly] = true;
    sel.selected_detector = "Zz9";
    const auto scope = EnergyCalTool::apply_to_scope( sel, display_state( multi, {1} ) );
    BOOST_CHECK( scope.empty() );
  }
}//BOOST_AUTO_TEST_CASE( applyToScopeMatrix )


namespace
{
/** A headless InterSpec session, to exercise applying the calibration changes. */
class EnergyCalSessionFixture
{
public:
  std::unique_ptr<Wt::Test::WTestEnvironment> m_env;
  InterSpecApp *m_app = nullptr;
  std::unique_ptr<Wt::WApplication::UpdateLock> m_lock;
  InterSpec *m_interspec = nullptr;
  EnergyCalTool *m_tool = nullptr;

  EnergyCalSessionFixture()
  {
    set_data_dir();

    string app_root = SpecUtils::append_path( InterSpec::staticDataDirectory(), ".." );
    app_root = SpecUtils::lexically_normalize_path( app_root );

    m_env.reset( new Wt::Test::WTestEnvironment( "", "", Wt::EntryPointType::Application ) );
    m_env->setAppRoot( app_root );
    m_app = new InterSpecApp( *m_env );
    m_lock.reset( new Wt::WApplication::UpdateLock(m_app) );
    m_interspec = m_app->viewer();
    BOOST_REQUIRE( m_interspec );
    m_tool = m_interspec->energyCalTool();
    BOOST_REQUIRE( m_tool );
  }

  void show( const shared_ptr<SpecMeas> &meas, const set<int> &samples, const SpecUtils::SpectrumType type )
  {
    m_interspec->setSpectrum( meas, samples, type, {} );
    BOOST_REQUIRE( m_interspec->measurment(type) == meas );
  }

  ~EnergyCalSessionFixture()
  {
    // `m_app` is owned by the environment - deleting it here as well crashes during teardown
    m_lock.reset();
    m_env.reset();
    m_interspec = nullptr;
    m_app = nullptr;
  }
};//class EnergyCalSessionFixture
}//namespace


/** With "Displayed Samples", only the displayed samples get the change, and only their peaks move -
 for coefficient changes, adding a deviation pair, and editing deviation pairs (previously the
 deviation pair paths moved the peaks of every sample sharing the calibration).
 */
BOOST_AUTO_TEST_CASE( displayedSamplesOnlyMovesTheirPeaks )
{
  enum class ChangeType { Coefficient, AddDevPair, EditDevPairs };

  for( const ChangeType change_type : { ChangeType::Coefficient, ChangeType::AddDevPair, ChangeType::EditDevPairs } )
  {
    EnergyCalSessionFixture session;

    const shared_ptr<SpecMeas> f = make_file( {1,2,3,4}, {"Aa1","Aa2"}, false );
    const shared_ptr<const SpecUtils::EnergyCalibration> orig_cal = cal_of( f, 1, "Aa1" );
    f->setPeaks( make_peaks( sm_peak_energy ), {2} );
    f->setPeaks( make_peaks( sm_peak_energy ), {3} );
    session.show( f, {2}, SpecUtils::SpectrumType::Foreground );

    Selection sel;
    sel.checked[CbIndex::ApplyToDisplayedSamples] = true;
    sel.checked[CbIndex::ApplyToAllSamples] = false;
    const vector<MeasToApplyCoefChangeTo> scope
                    = EnergyCalTool::apply_to_scope( sel, session.m_tool->currentApplyToDisplayState() );
    BOOST_REQUIRE_EQUAL( scope.size(), 1 );
    BOOST_REQUIRE( scope[0].sample_numbers == samps({2}) );

    switch( change_type )
    {
      case ChangeType::Coefficient:
        session.m_tool->applyCalChange( orig_cal, scaled_gain(orig_cal, 1.01f), scope, false );
        break;
      case ChangeType::AddDevPair:
        session.m_tool->addDeviationPair( {600.0f, 5.0f}, scope );
        break;
      case ChangeType::EditDevPairs:
        session.m_tool->applyDeviationPairEdit( {{0.0f, 0.0f}, {600.0f, 5.0f}, {1500.0f, 0.0f}}, scope );
        break;
    }//switch( change_type )

    for( const string det : {"Aa1", "Aa2"} )
    {
      BOOST_CHECK( cal_of( f, 2, det ) != orig_cal );
      for( const int sample : {1, 3, 4} )
        BOOST_CHECK( cal_of( f, sample, det ) == orig_cal );
    }

    // Detectors of a sample that shared a calibration, still share one
    BOOST_CHECK( cal_of( f, 2, "Aa1" ) == cal_of( f, 2, "Aa2" ) );

    BOOST_CHECK( fabs( peak_mean( f, {2} ) - sm_peak_energy ) > 0.5 );
    BOOST_CHECK_EQUAL( peak_mean( f, {3} ), sm_peak_energy );
  }//for( loop over change types )
}//BOOST_AUTO_TEST_CASE( displayedSamplesOnlyMovesTheirPeaks )


/** "Indep. Dets.": only the selected detector changes, matched by name in the background.
 The peaks were fit to the sum of all three detectors, and only one of them changes, so the peaks
 are left alone - whether or not the selected detector defines the summed calibration.
 */
BOOST_AUTO_TEST_CASE( selectedDetectorOnly )
{
  for( const string selected : { string("Aa3"), string("Aa1") } )
  {
    EnergyCalSessionFixture session;

    const shared_ptr<SpecMeas> f = make_file( {1,2}, {"Aa1","Aa2","Aa3"}, true );
    const shared_ptr<SpecMeas> b = make_file( {1,2}, {"Aa1","Aa2","Aa3"}, true );
    f->setPeaks( make_peaks( sm_peak_energy ), {1,2} );
    session.show( f, {1,2}, SpecUtils::SpectrumType::Foreground );
    session.show( b, {1,2}, SpecUtils::SpectrumType::Background );

    // The detector with the largest energy range (Aa3) defines the summed spectrums calibration
    const shared_ptr<const SpecUtils::Measurement> disp_fore
                     = session.m_interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
    BOOST_REQUIRE( disp_fore );
    BOOST_REQUIRE( disp_fore->energy_calibration() == cal_of( f, 1, "Aa3" ) );

    map<pair<shared_ptr<SpecMeas>,string>,shared_ptr<const SpecUtils::EnergyCalibration>> orig;
    for( const shared_ptr<SpecMeas> &meas : {f, b} )
      for( const string det : {"Aa1","Aa2","Aa3"} )
        orig[make_pair(meas,det)] = cal_of( meas, 1, det );

    Selection sel;
    sel.checked[CbIndex::ApplyToSelectedDetOnly] = true;
    sel.selected_detector = selected;
    const vector<MeasToApplyCoefChangeTo> scope
                  = EnergyCalTool::apply_to_scope( sel, session.m_tool->currentApplyToDisplayState() );
    BOOST_REQUIRE_EQUAL( scope.size(), 2 );

    const shared_ptr<const SpecUtils::EnergyCalibration> prev = orig[make_pair(f,selected)];
    const shared_ptr<const SpecUtils::EnergyCalibration> updated = scaled_gain( prev, 1.01f );
    BOOST_REQUIRE_NO_THROW( session.m_tool->applyCalChange( prev, updated, scope, false ) );

    for( const shared_ptr<SpecMeas> &meas : {f, b} )
    {
      for( const string det : {"Aa1","Aa2","Aa3"} )
      {
        for( const int sample : {1, 2} )
        {
          if( det == selected )
            BOOST_CHECK( cal_of( meas, sample, det ) != orig[make_pair(meas,det)] );
          else
            BOOST_CHECK( cal_of( meas, sample, det ) == orig[make_pair(meas,det)] );
        }
      }
    }//for( const shared_ptr<SpecMeas> &meas : {f, b} )

    BOOST_CHECK( cal_of( f, 1, selected ) == updated );
    // The background detector got the equivalent (propagated) change, not the foregrounds cal
    BOOST_CHECK( cal_of( b, 1, selected ) != updated );
    BOOST_CHECK_CLOSE( cal_of( b, 1, selected )->coefficients()[1],
                       1.01f * orig[make_pair(b,selected)]->coefficients()[1], 1.0E-3 );

    BOOST_CHECK_EQUAL( peak_mean( f, {1,2} ), sm_peak_energy );
  }//for( loop over selected detector )
}//BOOST_AUTO_TEST_CASE( selectedDetectorOnly )


/** With detectors sharing a single calibration object (common binning): changing only one of them
 ("Indep. Dets.", or a deviation pair edit) doesnt move the peaks of samples where the other
 detector is summed in, but does for samples that only have the changed detector.
 */
BOOST_AUTO_TEST_CASE( partialChangeOfSharedCalibration )
{
  enum class ChangeType { SelectedDetCoefficient, DevPairEdit };
  for( const ChangeType change_type : { ChangeType::SelectedDetCoefficient, ChangeType::DevPairEdit } )
  {
    EnergyCalSessionFixture session;
    
    // "Aa1" has no data for samples 3 and 4
    const shared_ptr<SpecMeas> f = make_file( {1,2,3,4}, {"Aa1","Ba2"}, false, {{3,"Aa1"},{4,"Aa1"}} );
    const shared_ptr<const SpecUtils::EnergyCalibration> orig_cal = cal_of( f, 1, "Aa1" );
    BOOST_REQUIRE( cal_of( f, 1, "Ba2" ) == orig_cal );
    f->setPeaks( make_peaks( sm_peak_energy ), {1,2} );
    f->setPeaks( make_peaks( sm_peak_energy ), {3,4} );
    session.show( f, {1,2,3,4}, SpecUtils::SpectrumType::Foreground );
    
    const DisplayState state = session.m_tool->currentApplyToDisplayState();
    
    switch( change_type )
    {
      case ChangeType::SelectedDetCoefficient:
      {
        Selection sel;
        sel.checked[CbIndex::ApplyToSelectedDetOnly] = true;
        sel.selected_detector = "Ba2";
        const vector<MeasToApplyCoefChangeTo> scope = EnergyCalTool::apply_to_scope( sel, state );
        BOOST_REQUIRE_EQUAL( scope.size(), 1 );
        BOOST_REQUIRE( scope[0].detectors == dets({"Ba2"}) );
        session.m_tool->applyCalChange( orig_cal, scaled_gain(orig_cal, 0.99f), scope, false );
        break;
      }//case ChangeType::SelectedDetCoefficient:
        
      case ChangeType::DevPairEdit:
      {
        const vector<MeasToApplyCoefChangeTo> scope
                       = EnergyCalTool::deviation_pair_edit_scope( Selection{}, state, f, "Ba2" );
        BOOST_REQUIRE_EQUAL( scope.size(), 1 );
        BOOST_REQUIRE( scope[0].detectors == dets({"Ba2"}) );
        session.m_tool->applyDeviationPairEdit( {{0.0f, 0.0f}, {600.0f, 5.0f}, {1500.0f, 0.0f}}, scope );
        break;
      }//case ChangeType::DevPairEdit:
    }//switch( change_type )
    
    for( const int sample : {1, 2} )
      BOOST_CHECK( cal_of( f, sample, "Aa1" ) == orig_cal );
    for( const int sample : {1, 2, 3, 4} )
      BOOST_CHECK( cal_of( f, sample, "Ba2" ) != orig_cal );
    
    BOOST_CHECK_EQUAL( peak_mean( f, {1,2} ), sm_peak_energy );   //Aa1 summed in, and didnt change
    BOOST_CHECK( fabs( peak_mean( f, {3,4} ) - sm_peak_energy ) > 0.5 ); //only Ba2, which changed
  }//for( loop over change types )
}//BOOST_AUTO_TEST_CASE( partialChangeOfSharedCalibration )


/** Deviation pair edits only apply to the edited detector, matched by name - except when files have
 only a single detector, in which case that detector is used regardless of its name.
 */
BOOST_AUTO_TEST_CASE( deviationPairEditScope )
{
  set_data_dir();
  const Selection defaults;
  
  // Same detector names: just the edited detector, in both files, even with "Displayed Detectors"
  //  checked and the edited detector not displayed
  {
    const shared_ptr<SpecMeas> f = make_file( {1}, {"Aa1","Aa2"}, true );
    const shared_ptr<SpecMeas> b = make_file( {1}, {"Aa1","Aa2"}, true );
    DisplayState state = display_state( f, {1}, b, {1} );
    state.detectors[0] = state.detectors[1] = { "Aa1" };
    
    Selection sel;
    sel.checked[CbIndex::ApplyToDisplayedDetectors] = true;
    sel.checked[CbIndex::ApplyToAllDetectors] = false;
    const auto scope = EnergyCalTool::deviation_pair_edit_scope( sel, state, f, "Aa2" );
    BOOST_REQUIRE_EQUAL( scope.size(), 2 );
    BOOST_CHECK( entry_for( scope, f )->detectors == dets({"Aa2"}) );
    BOOST_CHECK( entry_for( scope, b )->detectors == dets({"Aa2"}) );
    
    Selection no_back;
    no_back.checked[CbIndex::ApplyToBackground] = false;
    BOOST_CHECK_EQUAL( EnergyCalTool::deviation_pair_edit_scope( no_back, state, f, "Aa2" ).size(), 1 );
  }
  
  // Single detector files, with different names: both
  {
    const shared_ptr<SpecMeas> f = make_file( {1}, {""}, false );
    const shared_ptr<SpecMeas> b = make_file( {1}, {"Aa1"}, false );
    const auto scope = EnergyCalTool::deviation_pair_edit_scope( defaults, display_state( f, {1}, b, {1} ), f, "" );
    BOOST_REQUIRE_EQUAL( scope.size(), 2 );
    BOOST_CHECK( entry_for( scope, b )->detectors == dets({"Aa1"}) );
  }
  
  // Multi-detector foreground, single (differently named) detector background: not the background
  {
    const shared_ptr<SpecMeas> f = make_file( {1}, {"Aa1","Aa2"}, true );
    const shared_ptr<SpecMeas> b = make_file( {1}, {"Bb1"}, false );
    const auto scope = EnergyCalTool::deviation_pair_edit_scope( defaults, display_state( f, {1}, b, {1} ), f, "Aa1" );
    BOOST_REQUIRE_EQUAL( scope.size(), 1 );
    BOOST_CHECK( scope[0].meas == f );
    BOOST_CHECK( scope[0].detectors == dets({"Aa1"}) );
  }
}//BOOST_AUTO_TEST_CASE( deviationPairEditScope )


/** Whether a change made relative to the displayed foreground (e.g., dragging the spectrum) would be
 applied to it.
 */
BOOST_AUTO_TEST_CASE( scopeIncludesDisplayedForeground )
{
  set_data_dir();
  
  const shared_ptr<SpecMeas> f = make_file( {1,2}, {"A1"}, false );
  const shared_ptr<SpecMeas> b = make_file( {1}, {"B1","B2"}, true );
  const DisplayState state = display_state( f, {1}, b, {1} );
  
  BOOST_CHECK( EnergyCalTool::scope_includes_displayed_foreground(
                                       EnergyCalTool::apply_to_scope( Selection{}, state ), state ) );
  
  Selection no_fore;
  no_fore.checked[CbIndex::ApplyToForeground] = false;
  BOOST_CHECK( !EnergyCalTool::scope_includes_displayed_foreground(
                                       EnergyCalTool::apply_to_scope( no_fore, state ), state ) );
  
  // "Foreground" is checked, but the selected detector is only in the background
  Selection sel;
  sel.checked[CbIndex::ApplyToSelectedDetOnly] = true;
  sel.selected_detector = "B1";
  BOOST_CHECK( !EnergyCalTool::scope_includes_displayed_foreground(
                                       EnergyCalTool::apply_to_scope( sel, state ), state ) );
  
  // Same file as foreground (sample 1) and background (sample 2), with only the background selected
  const DisplayState same_file = display_state( f, {1}, f, {2} );
  BOOST_CHECK( !EnergyCalTool::scope_includes_displayed_foreground(
                                       EnergyCalTool::apply_to_scope( no_fore, same_file ), same_file ) );
}//BOOST_AUTO_TEST_CASE( scopeIncludesDisplayedForeground )


/** Duplicate entries for the same file are tolerated (e.g., the Isotopics tool passing both the
 foreground and background, when they are the same file), and dont double-apply the change.
 */
BOOST_AUTO_TEST_CASE( duplicateEntriesApplyOnce )
{
  EnergyCalSessionFixture session;

  const shared_ptr<SpecMeas> f = make_file( {1,2}, {"Aa1"}, false );
  f->setPeaks( make_peaks( sm_peak_energy ), {1} );
  session.show( f, {1}, SpecUtils::SpectrumType::Foreground );
  session.show( f, {2}, SpecUtils::SpectrumType::Background );

  const shared_ptr<const SpecUtils::EnergyCalibration> orig_cal = cal_of( f, 1, "Aa1" );
  const shared_ptr<const SpecUtils::EnergyCalibration> updated = scaled_gain( orig_cal, 1.01f );

  MeasToApplyCoefChangeTo fore, back;
  fore.meas = back.meas = f;
  fore.detectors = back.detectors = dets({"Aa1"});
  fore.sample_numbers = samps({1});
  back.sample_numbers = samps({1,2});

  BOOST_REQUIRE_NO_THROW( session.m_tool->applyCalChange( orig_cal, updated, {fore, back}, false ) );
  BOOST_CHECK( cal_of( f, 1, "Aa1" ) == updated );
  BOOST_CHECK( cal_of( f, 2, "Aa1" ) == updated );
  BOOST_CHECK_CLOSE( peak_mean( f, {1} ), sm_peak_energy * 1.01, 0.01 );
}//BOOST_AUTO_TEST_CASE( duplicateEntriesApplyOnce )


/** A multi-detector CALp file, applied to detectors sharing one calibration: each detector gets
 its own calibration, and the peaks are translated exactly once.
 */
BOOST_AUTO_TEST_CASE( calpPerDetectorTranslatesPeaksOnce )
{
  EnergyCalSessionFixture session;

  const shared_ptr<SpecMeas> f = make_file( {1,2}, {"Aa1","Aa2"}, false );
  f->setPeaks( make_peaks( sm_peak_energy ), {1,2} );
  session.show( f, {1,2}, SpecUtils::SpectrumType::Foreground );

  auto x1 = make_shared<SpecUtils::EnergyCalibration>();
  x1->set_polynomial( 1024, {1.0f, 3.0f}, {} );
  auto x2 = make_shared<SpecUtils::EnergyCalibration>();
  x2->set_polynomial( 1024, {2.0f, 3.0f}, {} );

  map<string,shared_ptr<const SpecUtils::EnergyCalibration>> det_to_cal;
  det_to_cal["Aa1"] = x1;
  det_to_cal["Aa2"] = x2;

  BOOST_REQUIRE_NO_THROW( session.m_tool->applyCALpEnergyCal( det_to_cal,
                                              SpecUtils::SpectrumType::Foreground, true, true ) );

  for( const int sample : {1, 2} )
  {
    BOOST_CHECK( cal_of( f, sample, "Aa1" ) == x1 );
    BOOST_CHECK( cal_of( f, sample, "Aa2" ) == x2 );
  }

  // Translated once, according to the new calibration of the first displayed measurement that had
  //  the summed calibration (Aa1: +1 keV); with the detectors now differing, there is no exact answer
  BOOST_CHECK_CLOSE( peak_mean( f, {1,2} ), sm_peak_energy + 1.0, 0.01 );
}//BOOST_AUTO_TEST_CASE( calpPerDetectorTranslatesPeaksOnce )


/** A change that would change nothing throws, and leaves everything alone. */
BOOST_AUTO_TEST_CASE( nothingToChangeThrows )
{
  EnergyCalSessionFixture session;

  const shared_ptr<SpecMeas> f = make_file( {1}, {"Aa1","Aa2"}, false );
  session.show( f, {1}, SpecUtils::SpectrumType::Foreground );
  const shared_ptr<const SpecUtils::EnergyCalibration> orig_cal = cal_of( f, 1, "Aa1" );

  Selection sel;
  sel.checked[CbIndex::ApplyToSelectedDetOnly] = true;
  sel.selected_detector = "Zz9";
  const vector<MeasToApplyCoefChangeTo> scope
                  = EnergyCalTool::apply_to_scope( sel, session.m_tool->currentApplyToDisplayState() );
  BOOST_CHECK( scope.empty() );

  BOOST_CHECK_THROW( session.m_tool->applyCalChange( orig_cal, scaled_gain(orig_cal, 1.01f), scope, false ),
                     std::exception );
  BOOST_CHECK( cal_of( f, 1, "Aa1" ) == orig_cal );
  BOOST_CHECK( cal_of( f, 1, "Aa2" ) == orig_cal );
}//BOOST_AUTO_TEST_CASE( nothingToChangeThrows )


/** The widget-default selection, through the tool, applies to everything; and a single undo
 restores both the calibrations and the peaks.
 */
BOOST_AUTO_TEST_CASE( defaultScopeAndUndo )
{
  EnergyCalSessionFixture session;

  const shared_ptr<SpecMeas> f = make_file( {1,2}, {"Aa1","Aa2"}, true );
  f->setPeaks( make_peaks( sm_peak_energy ), {1,2} );
  session.show( f, {1,2}, SpecUtils::SpectrumType::Foreground );

  const vector<MeasToApplyCoefChangeTo> scope = session.m_tool->measurementsToApplyCoeffChangeTo();
  BOOST_REQUIRE_EQUAL( scope.size(), 1 );
  BOOST_CHECK( scope[0].detectors == dets({"Aa1","Aa2"}) );
  BOOST_CHECK( scope[0].sample_numbers == samps({1,2}) );

  const shared_ptr<const SpecUtils::EnergyCalibration> a1 = cal_of( f, 1, "Aa1" ), a2 = cal_of( f, 1, "Aa2" );
  const shared_ptr<const SpecUtils::EnergyCalibration> updated = scaled_gain( a2, 1.01f );
  session.m_tool->applyCalChange( a2, updated, scope, false );

  BOOST_CHECK( cal_of( f, 2, "Aa2" ) == updated );
  BOOST_CHECK( cal_of( f, 1, "Aa1" ) != a1 );
  BOOST_CHECK( cal_of( f, 1, "Aa1" ) == cal_of( f, 2, "Aa1" ) );  //samples still share the cal
  BOOST_CHECK_CLOSE( peak_mean( f, {1,2} ), sm_peak_energy * 1.01, 0.01 );

  UndoRedoManager * const undo = session.m_interspec->undoRedoManager();
  if( !undo || !undo->canUndo() )
  {
    BOOST_TEST_MESSAGE( "Undo/redo not available in this test environment - not testing undo." );
    return;
  }

  undo->executeUndo();
  BOOST_CHECK( cal_of( f, 1, "Aa1" ) == a1 );
  BOOST_CHECK( cal_of( f, 2, "Aa2" ) == a2 );
  BOOST_CHECK_CLOSE( peak_mean( f, {1,2} ), sm_peak_energy, 1.0E-6 );
}//BOOST_AUTO_TEST_CASE( defaultScopeAndUndo )
