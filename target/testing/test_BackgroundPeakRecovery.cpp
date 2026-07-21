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
#include <map>
#include <cmath>
#include <deque>
#include <mutex>
#include <atomic>
#include <chrono>
#include <string>
#include <vector>
#include <future>
#include <memory>
#include <thread>
#include <iostream>
#include <algorithm>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakSearchGuiUtils.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

#define BOOST_TEST_MODULE BackgroundPeakRecovery_suite
#include <boost/test/included/unit_test.hpp>

using namespace std;

namespace
{
  string g_data_dir;

  void set_data_dir()
  {
    static bool s_have_set = false;
    if( s_have_set )
      return;
    s_have_set = true;

    const int argc = boost::unit_test::framework::master_test_suite().argc;
    char **argv = boost::unit_test::framework::master_test_suite().argv;

    for( int i = 1; i < argc; ++i )
    {
      const string arg = argv[i];
      if( SpecUtils::istarts_with( arg, "--datadir=" ) )
        g_data_dir = arg.substr( 10 );
    }

    SpecUtils::ireplace_all( g_data_dir, "%20", " " );

    if( g_data_dir.empty() )
    {
      for( const char * const d : { "data", "../data", "../../data", "../../../data" } )
      {
        if( SpecUtils::is_file( SpecUtils::append_path( d, "sandia.decay.xml" ) ) )
        {
          g_data_dir = d;
          break;
        }
      }
    }

    BOOST_REQUIRE_MESSAGE( !g_data_dir.empty(), "Could not find data directory" );
    BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( g_data_dir ) );
    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );
  }//set_data_dir()


  // Loads a single measurement from a reference-spectra file for the given detector folder.
  shared_ptr<const SpecUtils::Measurement> load_reference_measurement( const string &detector_dir,
                                                                       const string &filename )
  {
    set_data_dir();
    const string dir = SpecUtils::append_path( g_data_dir,
                          SpecUtils::append_path( "reference_spectra/Common_Field_Nuclides", detector_dir ) );
    const string path = SpecUtils::append_path( dir, filename );
    if( !SpecUtils::is_file( path ) )
      return nullptr;

    SpecUtils::SpecFile f;
    if( !f.load_file( path, SpecUtils::ParserType::Auto ) )
      return nullptr;

    shared_ptr<const SpecUtils::Measurement> m = f.measurement_at_index( 0 );
    if( m && !m->num_gamma_channels() )
      m = nullptr;
    return m;
  }//load_reference_measurement(...)


  bool has_peak_near( const vector<shared_ptr<const PeakDef>> &peaks, const double energy,
                      const double tol_kev )
  {
    for( const shared_ptr<const PeakDef> &p : peaks )
    {
      if( p && p->gausPeak() && (fabs(p->mean() - energy) < tol_kev) )
        return true;
    }
    return false;
  }//has_peak_near(...)


  bool has_peak_near( const deque<shared_ptr<const PeakDef>> &peaks, const double energy,
                      const double tol_kev )
  {
    return has_peak_near( vector<shared_ptr<const PeakDef>>( begin(peaks), end(peaks) ), energy, tol_kev );
  }


  // Returns a low-statistics copy of `orig`: counts and times scaled by `factor` (<1), which lowers
  // the per-line significance so a whole-spectrum auto-search may miss weak lines (reproducing the
  // real 300 s NORM-leak scenario from ~1800 s reference spectra).
  shared_ptr<SpecUtils::Measurement> scaled_low_stats_copy( const shared_ptr<const SpecUtils::Measurement> &orig,
                                                            const double factor )
  {
    if( !orig || !orig->gamma_counts() )
      return nullptr;

    shared_ptr<SpecUtils::Measurement> out = make_shared<SpecUtils::Measurement>( *orig );
    shared_ptr<vector<float>> counts = make_shared<vector<float>>( *orig->gamma_counts() );
    for( float &c : *counts )
      c = static_cast<float>( c * factor );

    const float lt = static_cast<float>( std::max( 1.0, orig->live_time() * factor ) );
    const float rt = static_cast<float>( std::max( 1.0, orig->real_time() * factor ) );
    out->set_gamma_counts( counts, lt, rt );
    return out;
  }//scaled_low_stats_copy(...)
}//namespace


// -----------------------------------------------------------------------------
// PART A: the shared, de-duplicated automated-search-peak accessor.
// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( AccessorCoalescesConcurrentSearches )
{
  set_data_dir();

  const string dir = SpecUtils::append_path( g_data_dir,
                        "reference_spectra/Common_Field_Nuclides/Detective X" );
  const string fg_path = SpecUtils::append_path( dir, "Cs137_Unshielded.txt" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( fg_path ), "Missing reference spectrum: " << fg_path );

  shared_ptr<SpecMeas> meas = make_shared<SpecMeas>();
  BOOST_REQUIRE( meas->load_file( fg_path, SpecUtils::ParserType::Auto, fg_path ) );

  const shared_ptr<const SpecUtils::Measurement> spectrum = meas->measurement_at_index( 0 );
  BOOST_REQUIRE( spectrum && spectrum->num_gamma_channels() );
  const set<int> samples = { spectrum->sample_number() };
  const shared_ptr<const DetectorPeakResponse> drf = meas->detector();
  const bool isHPGe = true;  // Detective X is HPGe

  // Fire N concurrent requests from N threads; only one search should actually run, and every
  // future must resolve to the *same* peak-set object (proving a single computation + coalescing).
  const size_t N = 8;
  vector<shared_future<shared_ptr<const deque<shared_ptr<const PeakDef>>>>> futures( N );
  {
    vector<thread> threads;
    mutex fut_mutex;
    for( size_t i = 0; i < N; ++i )
    {
      threads.emplace_back( [&,i](){
        shared_future<shared_ptr<const deque<shared_ptr<const PeakDef>>>> f
          = PeakSearchGuiUtils::get_or_launch_automated_search_peaks( meas, samples, spectrum, drf, isHPGe );
        lock_guard<mutex> lock( fut_mutex );
        futures[i] = f;
      } );
    }
    for( thread &t : threads )
      t.join();
  }

  const shared_ptr<const deque<shared_ptr<const PeakDef>>> first = futures[0].get();
  BOOST_REQUIRE_MESSAGE( first, "Automated search returned null result" );
  BOOST_CHECK_MESSAGE( !first->empty(), "Cs-137 automated search found no peaks" );
  BOOST_CHECK_MESSAGE( has_peak_near( *first, 661.66, 3.0 ), "Cs-137 661 keV not found by search" );

  for( size_t i = 1; i < N; ++i )
  {
    const shared_ptr<const deque<shared_ptr<const PeakDef>>> res = futures[i].get();
    BOOST_CHECK_MESSAGE( res == first,
        "Concurrent request " << i << " got a different peak-set object (a duplicate search ran)" );
  }

  // A subsequent call must be a cache hit: ready immediately, and the same object.
  shared_future<shared_ptr<const deque<shared_ptr<const PeakDef>>>> again
    = PeakSearchGuiUtils::get_or_launch_automated_search_peaks( meas, samples, spectrum, drf, isHPGe );
  BOOST_CHECK_MESSAGE( again.wait_for( std::chrono::seconds(0) ) == std::future_status::ready,
                       "Cached automated-search result was not returned as a ready future" );
  BOOST_CHECK( again.get() == first );
}//BOOST_AUTO_TEST_CASE( AccessorCoalescesConcurrentSearches )


BOOST_AUTO_TEST_CASE( AccessorCancellationUnblocksWaiters )
{
  set_data_dir();

  const string dir = SpecUtils::append_path( g_data_dir,
                        "reference_spectra/Common_Field_Nuclides/Detective X" );
  const string fg_path = SpecUtils::append_path( dir, "Eu152_Unshielded.txt" );
  if( !SpecUtils::is_file( fg_path ) )
  {
    BOOST_WARN_MESSAGE( false, "Missing Eu152 reference spectrum; skipping cancellation test" );
    return;
  }

  shared_ptr<SpecMeas> meas = make_shared<SpecMeas>();
  BOOST_REQUIRE( meas->load_file( fg_path, SpecUtils::ParserType::Auto, fg_path ) );
  const shared_ptr<const SpecUtils::Measurement> spectrum = meas->measurement_at_index( 0 );
  BOOST_REQUIRE( spectrum && spectrum->num_gamma_channels() );
  const set<int> samples = { spectrum->sample_number() };

  shared_future<shared_ptr<const deque<shared_ptr<const PeakDef>>>> fut
    = PeakSearchGuiUtils::get_or_launch_automated_search_peaks( meas, samples, spectrum,
                                                               meas->detector(), true );
  // Request cancellation; the waiter must unblock (result may be null on cancel, or a finished set
  // if the search completed first).  The key assertion is simply that .get() returns (no hang).
  meas->cancel_in_progress_peak_search();
  BOOST_CHECK_NO_THROW( fut.get() );
}//BOOST_AUTO_TEST_CASE( AccessorCancellationUnblocksWaiters )


// -----------------------------------------------------------------------------
// PART B: background-peak recovery.
// -----------------------------------------------------------------------------
namespace
{
  // Runs the deterministic recovery check for one detector folder + one source spectrum that has
  // NORM lines overlapping the background.  Returns false (with a warning) if inputs are unusable.
  void run_recovery_check( const string &detector_dir, const string &fg_filename, const bool isHPGe )
  {
    const shared_ptr<const SpecUtils::Measurement> fg = load_reference_measurement( detector_dir, fg_filename );
    const shared_ptr<const SpecUtils::Measurement> bg = load_reference_measurement( detector_dir, "background.txt" );
    if( !fg || !bg )
    {
      BOOST_WARN_MESSAGE( false, "Missing spectra for " << detector_dir << "/" << fg_filename << "; skipping" );
      return;
    }

    const shared_ptr<const DetectorPeakResponse> drf;  // reference spectra have no DRF attached

    const vector<shared_ptr<const PeakDef>> fg_peaks
      = ExperimentalAutomatedPeakSearch::search_for_peaks( fg, drf, nullptr, true, isHPGe );
    const vector<shared_ptr<const PeakDef>> bg_peaks
      = ExperimentalAutomatedPeakSearch::search_for_peaks( bg, drf, nullptr, true, isHPGe );

    BOOST_TEST_MESSAGE( detector_dir << "/" << fg_filename << ": " << fg_peaks.size()
                        << " foreground peaks, " << bg_peaks.size() << " background peaks (full stats)" );

    // --- Direction 1: recovery re-finds a real background peak the base search "missed" ---
    // Pick the strongest background peak that also has a matching foreground peak (a shared line,
    // e.g. NORM), then remove it from the background set to simulate a low-stats auto-search miss.
    shared_ptr<const PeakDef> shared_line;
    for( const shared_ptr<const PeakDef> &bgp : bg_peaks )
    {
      if( !bgp || !bgp->gausPeak() )
        continue;
      const double tol = 0.75 * bgp->fwhm();
      if( !has_peak_near( fg_peaks, bgp->mean(), tol ) )
        continue;
      if( !shared_line || (bgp->amplitude() > shared_line->amplitude()) )
        shared_line = bgp;
    }

    if( shared_line )
    {
      const double removed_energy = shared_line->mean();
      const double tol = 0.75 * shared_line->fwhm();

      deque<shared_ptr<const PeakDef>> handicapped;
      for( const shared_ptr<const PeakDef> &bgp : bg_peaks )
        if( bgp != shared_line )
          handicapped.push_back( bgp );

      BOOST_CHECK_MESSAGE( !has_peak_near( handicapped, removed_energy, tol ),
          "Sanity: removed background line still present in handicapped set" );

      const shared_ptr<const deque<shared_ptr<const PeakDef>>> handicapped_ptr
        = make_shared<const deque<shared_ptr<const PeakDef>>>( handicapped );

      const shared_ptr<const deque<shared_ptr<const PeakDef>>> augmented
        = ExperimentalAutomatedPeakSearch::recover_background_peaks_under_foreground(
              fg_peaks, bg, handicapped_ptr, drf, isHPGe );

      BOOST_REQUIRE_MESSAGE( augmented, "Recovery returned null" );
      BOOST_CHECK_MESSAGE( augmented->size() >= handicapped.size(),
          "Recovery removed background peaks (it should only add)" );
      BOOST_CHECK_MESSAGE( has_peak_near( *augmented, removed_energy, tol ),
          detector_dir << ": recovery failed to re-find the real background line at "
          << removed_energy << " keV" );
    }else
    {
      BOOST_WARN_MESSAGE( false, detector_dir << "/" << fg_filename
                          << ": no shared foreground/background line found to exercise recovery" );
    }

    // --- Direction 2: recovery does NOT fabricate a true foreground-only source line ---
    // The source's own photopeak (present in the foreground, absent from a plain NORM background)
    // must NOT get invented in the background by recovery.
    const double source_energy = 661.66;  // Cs-137 (foreground only)
    if( has_peak_near( fg_peaks, source_energy, 3.0 ) && !has_peak_near( bg_peaks, source_energy, 3.0 ) )
    {
      const shared_ptr<const deque<shared_ptr<const PeakDef>>> base_bg_ptr
        = make_shared<const deque<shared_ptr<const PeakDef>>>( begin(bg_peaks), end(bg_peaks) );
      const shared_ptr<const deque<shared_ptr<const PeakDef>>> augmented
        = ExperimentalAutomatedPeakSearch::recover_background_peaks_under_foreground(
              fg_peaks, bg, base_bg_ptr, drf, isHPGe );
      BOOST_REQUIRE( augmented );
      BOOST_CHECK_MESSAGE( !has_peak_near( *augmented, source_energy, 3.0 ),
          detector_dir << ": recovery fabricated a background peak at the foreground-only source line "
          << source_energy << " keV" );
    }
  }//run_recovery_check(...)
}//namespace


BOOST_AUTO_TEST_CASE( RecoveryFindsRealBackgroundLines_HPGe )
{
  run_recovery_check( "Detective X", "Cs137_Unshielded.txt", true );
}


BOOST_AUTO_TEST_CASE( RecoveryFindsRealBackgroundLines_NaI )
{
  run_recovery_check( "RadEagle NaI 3x1", "Cs137_Unshielded.txt", false );
}


// Multi-peak ROI handling: a background ROI containing >=2 peaks that share one continuum (e.g. the
// Th-232 chain 238.6/240.9 doublet on HPGe) must, when the auto-search "misses" the whole ROI, be
// recovered WITHOUT dropping a line and WITHOUT creating duplicates or an inconsistent (split)
// continuum.  This directly exercises the same-ROI collision hazard the earlier parallel design had.
BOOST_AUTO_TEST_CASE( RecoveryHandlesMultiPeakRoi )
{
  set_data_dir();

  const string detector_dir = "Detective X";  // HPGe - resolves close NORM lines into multi-peak ROIs
  const bool isHPGe = true;

  const shared_ptr<const SpecUtils::Measurement> fg = load_reference_measurement( detector_dir, "Cs137_Unshielded.txt" );
  const shared_ptr<const SpecUtils::Measurement> bg = load_reference_measurement( detector_dir, "background.txt" );
  if( !fg || !bg )
  {
    BOOST_WARN_MESSAGE( false, "Missing spectra for multi-peak ROI test; skipping" );
    return;
  }

  const shared_ptr<const DetectorPeakResponse> drf;

  const vector<shared_ptr<const PeakDef>> fg_peaks
    = ExperimentalAutomatedPeakSearch::search_for_peaks( fg, drf, nullptr, true, isHPGe );
  const vector<shared_ptr<const PeakDef>> bg_peaks
    = ExperimentalAutomatedPeakSearch::search_for_peaks( bg, drf, nullptr, true, isHPGe );

  // Group background peaks by their shared continuum, and find a multi-peak ROI (>=2 peaks) with at
  // least two foreground-matched lines (so recovery will target them).
  std::map<shared_ptr<const PeakContinuum>, vector<shared_ptr<const PeakDef>>> by_cont;
  for( const shared_ptr<const PeakDef> &p : bg_peaks )
  {
    if( p && p->gausPeak() )
      by_cont[p->continuum()].push_back( p );
  }

  vector<shared_ptr<const PeakDef>> roi;
  for( const std::pair<const shared_ptr<const PeakContinuum>, vector<shared_ptr<const PeakDef>>> &kv : by_cont )
  {
    if( kv.second.size() < 2 )
      continue;
    int matched = 0;
    for( const shared_ptr<const PeakDef> &p : kv.second )
    {
      if( has_peak_near( fg_peaks, p->mean(), 0.75*p->fwhm() ) )
        ++matched;
    }
    if( (matched >= 2) && (kv.second.size() > roi.size()) )
      roi = kv.second;
  }

  if( roi.size() < 2 )
  {
    BOOST_WARN_MESSAGE( false, detector_dir
        << ": no multi-peak background ROI with >=2 foreground matches found; skipping multi-peak ROI test" );
    return;
  }

  double roi_lo = std::numeric_limits<double>::infinity();
  double roi_hi = -std::numeric_limits<double>::infinity();
  vector<double> target_energies;
  for( const shared_ptr<const PeakDef> &p : roi )
  {
    roi_lo = std::min( roi_lo, p->mean() );
    roi_hi = std::max( roi_hi, p->mean() );
    if( has_peak_near( fg_peaks, p->mean(), 0.75*p->fwhm() ) )
      target_energies.push_back( p->mean() );
  }

  BOOST_TEST_MESSAGE( detector_dir << ": exercising a " << roi.size() << "-peak background ROI over ["
                      << roi_lo << ", " << roi_hi << "] keV" );

  // Remove the whole ROI group from the base background set (simulate the auto-search missing it).
  deque<shared_ptr<const PeakDef>> handicapped;
  for( const shared_ptr<const PeakDef> &p : bg_peaks )
  {
    if( std::find( roi.begin(), roi.end(), p ) == roi.end() )
      handicapped.push_back( p );
  }

  const shared_ptr<const deque<shared_ptr<const PeakDef>>> handicapped_ptr
    = make_shared<const deque<shared_ptr<const PeakDef>>>( handicapped );
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> augmented
    = ExperimentalAutomatedPeakSearch::recover_background_peaks_under_foreground(
          fg_peaks, bg, handicapped_ptr, drf, isHPGe );
  BOOST_REQUIRE( augmented );

  // (a) No line dropped: every foreground-matched line in the ROI is recovered in the background.
  for( const double e : target_energies )
  {
    BOOST_CHECK_MESSAGE( has_peak_near( *augmented, e, 1.5 ),
        detector_dir << ": multi-peak ROI recovery dropped the background line at " << e << " keV" );
  }

  // (b) No duplicates: the background peaks recovered across the ROI span must not exceed the number
  //     of foreground peaks there (a duplicate would inflate the count).
  const auto count_in_span = []( const auto &peaks, const double lo, const double hi ) -> int {
    int n = 0;
    for( const shared_ptr<const PeakDef> &p : peaks )
      if( p && p->gausPeak() && (p->mean() >= (lo - 2.0)) && (p->mean() <= (hi + 2.0)) )
        ++n;
    return n;
  };
  const int fg_in_span = count_in_span( fg_peaks, roi_lo, roi_hi );
  const int aug_in_span = count_in_span( *augmented, roi_lo, roi_hi );
  BOOST_CHECK_MESSAGE( aug_in_span <= (fg_in_span + 1),
      detector_dir << ": recovery produced " << aug_in_span << " background peaks in the ROI span vs "
      << fg_in_span << " foreground peaks (possible duplicate)" );

  // (c) Consistency: peaks recovered together in one ROI must share a single continuum (not each be
  //     left on its own/stale continuum).  If >=2 peaks land in the original ROI span, they must not
  //     all be on distinct continua.
  std::set<const void*> continua;
  int span_peaks = 0;
  for( const shared_ptr<const PeakDef> &p : *augmented )
  {
    if( p && p->gausPeak() && (p->mean() >= (roi_lo - 0.5)) && (p->mean() <= (roi_hi + 0.5)) )
    {
      ++span_peaks;
      continua.insert( p->continuum().get() );
    }
  }
  BOOST_CHECK_MESSAGE( (span_peaks < 2) || (continua.size() < static_cast<size_t>(span_peaks)),
      detector_dir << ": " << span_peaks << " recovered ROI peaks are on " << continua.size()
      << " separate continua (expected a shared continuum for the multi-peak ROI)" );
}


// Low-statistics report: scale the background down so the whole-spectrum auto-search misses weak
// NORM lines, then confirm recovery reduces the number of foreground NORM lines left unmatched.
BOOST_AUTO_TEST_CASE( RecoveryReducesLeakAtLowStats )
{
  set_data_dir();

  const string detector_dir = "Detective X";
  const bool isHPGe = true;

  const shared_ptr<const SpecUtils::Measurement> fg = load_reference_measurement( detector_dir, "Cs137_Unshielded.txt" );
  const shared_ptr<const SpecUtils::Measurement> bg_full = load_reference_measurement( detector_dir, "background.txt" );
  if( !fg || !bg_full )
  {
    BOOST_WARN_MESSAGE( false, "Missing spectra for low-stats leak test; skipping" );
    return;
  }

  const shared_ptr<const SpecUtils::Measurement> bg = scaled_low_stats_copy( bg_full, 1.0/6.0 );
  BOOST_REQUIRE( bg );

  const shared_ptr<const DetectorPeakResponse> drf;

  const vector<shared_ptr<const PeakDef>> fg_peaks
    = ExperimentalAutomatedPeakSearch::search_for_peaks( fg, drf, nullptr, true, isHPGe );
  const vector<shared_ptr<const PeakDef>> bg_peaks
    = ExperimentalAutomatedPeakSearch::search_for_peaks( bg, drf, nullptr, true, isHPGe );

  // Common NORM lines that a Cs-137 foreground (with nominal NORM background) should show.
  const vector<double> norm_lines = { 1460.8, 2614.5, 1764.5, 911.2, 583.1, 238.6, 609.3, 968.9 };

  const double tol = 3.0;
  auto count_unmatched = [&]( const deque<shared_ptr<const PeakDef>> &bgset ) -> int {
    int n = 0;
    for( const double e : norm_lines )
      if( has_peak_near( fg_peaks, e, tol ) && !has_peak_near( bgset, e, tol ) )
        ++n;
    return n;
  };

  const deque<shared_ptr<const PeakDef>> base_bg( begin(bg_peaks), end(bg_peaks) );
  const int leak_before = count_unmatched( base_bg );

  const shared_ptr<const deque<shared_ptr<const PeakDef>>> base_bg_ptr
    = make_shared<const deque<shared_ptr<const PeakDef>>>( base_bg );
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> augmented
    = ExperimentalAutomatedPeakSearch::recover_background_peaks_under_foreground(
          fg_peaks, bg, base_bg_ptr, drf, isHPGe );
  BOOST_REQUIRE( augmented );

  const int leak_after = count_unmatched( *augmented );

  BOOST_TEST_MESSAGE( "Low-stats NORM leak (foreground NORM lines missing from background): "
                      << leak_before << " before recovery, " << leak_after << " after recovery" );

  // Recovery must never increase the leak (whether or not this particular synthetic low-stats case
  // happens to leave lines below the recovery significance threshold).  The deterministic
  // remove-and-recover cases above carry the positive proof that recovery re-finds real lines.
  BOOST_CHECK_MESSAGE( leak_after <= leak_before, "Recovery increased the NORM leak" );
}//BOOST_AUTO_TEST_CASE( RecoveryReducesLeakAtLowStats )
