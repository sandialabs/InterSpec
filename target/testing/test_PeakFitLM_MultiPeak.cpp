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
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

#define BOOST_TEST_MODULE PeakFitLM_MultiPeak_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/PeakFitUtils.h"

using namespace std;
using namespace boost::unit_test;

// Regression test for PeakFitLM review issue #1 (residual-index off-by-one).
//
// `PeakFitLM::PeakFitDiffCostFunction` lays out its Ceres residuals as:
//   data residuals          -> [0, nchannel)
//   peak-proximity residuals -> [nchannel, nchannel + npeaks - 1)        (npeaks-1 of them)
//   optional stat-insig slot -> nchannel + npeaks - 1
// The proximity-punishment residual (written for every ROI with >=2 peaks, since punishing for
// peaks being too close is on by default) was indexed off `end_data_residual`, which was set to
// `m_upper_channel - m_lower_channel` (== nchannel - 1) instead of `nchannel`.  Consequently the
// first proximity write landed on `residuals[nchannel-1]` -- the last data channel's residual --
// clobbering it (and leaving the real last proximity slot unwritten), so the highest-energy
// channel was silently dropped from chi^2 and the Jacobian on every multi-peak ROI.
//
// The fix ties `end_data_residual` to `nchannel` and adds the invariant asserts
//   assert( prox_resid_index >= nchannel );           // never clobber a data residual
//   assert( prox_resid_index < number_residuals() );  // stay within the residual buffer
// in the proximity loop.  This test refits the *overlapping* (shared-continuum) ROIs of a real
// Ba-133 spectrum, which drives the multi-peak proximity path end-to-end.  The test build forces
// PERFORM_DEVELOPER_CHECKS=ON (CMakeLists.txt), so those asserts run on every Ceres iteration and
// abort the test if the off-by-one is ever reintroduced -- that assert IS the regression guard
// (a single dropped edge channel is too small to detect from the converged parameters alone,
// and the continuum/amplitude least-squares solve already uses every channel regardless).

namespace
{
string g_data_dir;
string g_test_file_dir;

void set_data_dir()
{
  static bool s_set = false;
  if( s_set )
    return;
  s_set = true;

  const int argc = framework::master_test_suite().argc;
  char ** const argv = framework::master_test_suite().argv;
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

  if( g_data_dir.empty()
     || !SpecUtils::is_file( SpecUtils::append_path(g_data_dir, "sandia.decay.xml") ) )
  {
    for( const char * const d : { "data", "../data", "../../data", "../../../data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) )
      {
        g_data_dir = d;
        break;
      }
    }
  }

  const char * const ba133_rel = "PeakFitLM/Ba133_Unshielded.n42";
  if( g_test_file_dir.empty()
     || !SpecUtils::is_file( SpecUtils::append_path(g_test_file_dir, ba133_rel) ) )
  {
    for( const char * const d : { "test_data", "../test_data", "../../test_data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, ba133_rel) ) )
      {
        g_test_file_dir = d;
        break;
      }
    }
  }

  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( SpecUtils::append_path(g_data_dir, "sandia.decay.xml") ),
                         "Could not locate the data directory (pass --datadir=...)" );
  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( g_data_dir ) );
}//set_data_dir()


struct LoadedSpectrum
{
  shared_ptr<const SpecUtils::Measurement> data;
  vector<shared_ptr<const PeakDef>> peaks;  // copied out, so the SpecMeas need not outlive us
};


LoadedSpectrum load_ba133()
{
  set_data_dir();

  const string path = SpecUtils::append_path( g_test_file_dir, "PeakFitLM/Ba133_Unshielded.n42" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file(path), "Ba133 spectrum not at '" << path << "'" );

  SpecMeas meas;
  BOOST_REQUIRE_MESSAGE( meas.load_N42_file(path), "Failed to load Ba133 spectrum '" << path << "'" );

  // The saved InterSpec peaks are keyed by the sample-number set they were fit on.
  const set<set<int>> peak_sample_sets = meas.sampleNumsWithPeaks();
  BOOST_REQUIRE_MESSAGE( !peak_sample_sets.empty(), "Ba133 spectrum has no saved peaks" );
  const set<int> &peak_samples = *peak_sample_sets.begin();

  const shared_ptr<deque<shared_ptr<const PeakDef>>> peaks = meas.peaks( peak_samples );
  BOOST_REQUIRE_MESSAGE( peaks && !peaks->empty(), "Ba133 spectrum has no peaks for its peak sample-set" );

  // The data measurement that owns those peaks' sample number.
  shared_ptr<const SpecUtils::Measurement> data;
  for( const shared_ptr<const SpecUtils::Measurement> &m : meas.measurements() )
  {
    if( m && (m->num_gamma_channels() > 16) && peak_samples.count(m->sample_number()) )
    {
      data = m;
      break;
    }
  }
  BOOST_REQUIRE_MESSAGE( data, "Could not find the Ba133 data measurement holding the saved peaks" );
  BOOST_REQUIRE( data->energy_calibration() && data->energy_calibration()->valid() );

  LoadedSpectrum out;
  out.data = data;
  out.peaks.assign( begin(*peaks), end(*peaks) );
  return out;
}//load_ba133()


// Group peaks into ROIs by their shared continuum (peaks in one ROI share a PeakContinuum).
vector<vector<shared_ptr<const PeakDef>>>
group_by_continuum( const vector<shared_ptr<const PeakDef>> &peaks )
{
  vector<vector<shared_ptr<const PeakDef>>> groups;
  map<shared_ptr<const PeakContinuum>,size_t> cont_to_group;
  for( const shared_ptr<const PeakDef> &p : peaks )
  {
    const shared_ptr<const PeakContinuum> cont = p->continuum();
    const auto pos = cont_to_group.find( cont );
    if( pos == end(cont_to_group) )
    {
      cont_to_group[cont] = groups.size();
      groups.push_back( { p } );
    }else
    {
      groups[pos->second].push_back( p );
    }
  }//for( loop over peaks )

  return groups;
}//group_by_continuum(...)

}//namespace


// Refit every overlapping (>=2 peak) ROI of a real Ba-133 spectrum through the default
//  `fit_peaks_in_roi_LM` path (peak-proximity punishment on).  Exercises the multi-peak residual
//  layout that issue #1 corrupted; the in-fitter `prox_resid_index` asserts validate the indexing.
BOOST_AUTO_TEST_CASE( refitOverlappingROIs )
{
  const LoadedSpectrum spec = load_ba133();
  const bool isHPGe = PeakFitUtils::is_high_res( spec.data );

  const vector<vector<shared_ptr<const PeakDef>>> rois = group_by_continuum( spec.peaks );

  size_t num_multi = 0, max_peaks = 0;
  for( const vector<shared_ptr<const PeakDef>> &roi : rois )
  {
    num_multi += (roi.size() >= 2);
    max_peaks = std::max( max_peaks, roi.size() );
  }
  cout << "Ba133: " << spec.peaks.size() << " peaks in " << rois.size() << " ROIs; "
       << num_multi << " overlapping (>=2 peak) ROIs, largest has " << max_peaks << " peaks." << endl;

  // The whole point of the test is the multi-peak path; make sure the fixture really has it.
  BOOST_REQUIRE_MESSAGE( num_multi >= 1,
                         "Expected at least one overlapping ROI in the Ba133 spectrum, found " << num_multi );

  size_t num_tested = 0;
  for( const vector<shared_ptr<const PeakDef>> &roi_peaks : rois )
  {
    if( roi_peaks.size() < 2 )
      continue;

    ++num_tested;

    // Default fit_options (== 0) keeps "punish for peaks being too close" enabled -> the proximity
    //  residuals are written, which is exactly the code path issue #1 broke.  Under
    //  PERFORM_DEVELOPER_CHECKS the prox_resid_index asserts fire here on any index regression.
    vector<shared_ptr<const PeakDef>> fit;
    BOOST_REQUIRE_NO_THROW( fit = PeakFitLM::fit_peaks_in_roi_LM( roi_peaks, spec.data, isHPGe ) );

    BOOST_CHECK_MESSAGE( fit.size() == roi_peaks.size(),
                         "Refit returned " << fit.size() << " peaks for a " << roi_peaks.size()
                         << "-peak ROI" );
    if( fit.size() != roi_peaks.size() )
      continue;

    vector<shared_ptr<const PeakDef>> orig = roi_peaks;
    std::sort( begin(orig), end(orig), &PeakDef::lessThanByMeanShrdPtr );
    std::sort( begin(fit),  end(fit),  &PeakDef::lessThanByMeanShrdPtr );

    for( size_t i = 0; i < fit.size(); ++i )
    {
      const double orig_sigma = orig[i]->gausPeak() ? orig[i]->sigma() : 0.0;
      const double mean_tol = std::max( 2.0, 0.5*orig_sigma );
      BOOST_CHECK_MESSAGE( std::fabs(fit[i]->mean() - orig[i]->mean()) < mean_tol,
                           "Refit peak mean " << fit[i]->mean() << " keV drifted > " << mean_tol
                           << " keV from the saved mean " << orig[i]->mean() << " keV" );

      BOOST_CHECK_MESSAGE( std::isfinite(fit[i]->amplitude()) && (fit[i]->amplitude() > 0.0),
                           "Refit peak (mean " << fit[i]->mean() << ") has non-positive/non-finite amplitude "
                           << fit[i]->amplitude() );

      if( fit[i]->gausPeak() )
        BOOST_CHECK_MESSAGE( std::isfinite(fit[i]->sigma()) && (fit[i]->sigma() > 0.0),
                             "Refit peak (mean " << fit[i]->mean() << ") has non-positive/non-finite sigma "
                             << fit[i]->sigma() );

      BOOST_CHECK_MESSAGE( std::isfinite(fit[i]->chi2dof()) && (fit[i]->chi2dof() > 0.0),
                           "Refit peak (mean " << fit[i]->mean() << ") has non-positive/non-finite chi2dof "
                           << fit[i]->chi2dof() );
    }//for( loop over fitted peaks )
  }//for( loop over ROIs )

  BOOST_CHECK_EQUAL( num_tested, num_multi );
}//BOOST_AUTO_TEST_CASE( refitOverlappingROIs )


// Issue #2: the peak-proximity punishment (PeakFitLM::peaks_too_close_punishment) must be
//  continuous (no jump at the old reldist==1.25 hard cutoff) and have compact support (identically
//  zero at/beyond the 1.75 sigma cutoff), so Ceres' L-M trust region stays well-behaved.  The dense
//  scan below is formula-agnostic, so any reintroduced discontinuity or C1 kink trips it.
BOOST_AUTO_TEST_CASE( proximityPunishmentContinuity )
{
  const double pf = 7.0;        // arbitrary positive punishment_factor
  const double cutoff = 1.75;   // must match peaks_too_close_punishment()

  auto P = [pf]( const double reldist ){ return PeakFitLM::peaks_too_close_punishment( reldist, pf ); };

  // Compact support: identically zero at and beyond the cutoff.
  for( const double r : { 1.75, 1.8, 2.0, 2.5, 3.0, 5.0 } )
    BOOST_CHECK_MESSAGE( P(r) == 0.0, "punishment at reldist=" << r << " must be 0, got " << P(r) );

  // Positive and strictly decreasing on (0, cutoff).
  double prev = std::numeric_limits<double>::infinity();
  for( double r = 0.05; r < (cutoff - 1.0e-9); r += 0.005 )
  {
    const double v = P(r);
    BOOST_CHECK_MESSAGE( v > 0.0, "punishment at reldist=" << r << " must be > 0, got " << v );
    BOOST_CHECK_MESSAGE( v < prev, "punishment must strictly decrease; reldist=" << r << " gave "
                         << v << " (>= previous " << prev << ")" );
    prev = v;
  }

  // C0 (no jumps) and C1 (no slope kinks) across [0.5, 3.0] - in particular catches the old hard
  //  cutoff at 1.25 (a ~0.8*pf jump) and any kink at the 1.75 cutoff.  On [0.5,3] the true function
  //  has max|dP/dr| ~ 5*pf and max|d2P/dr2| ~ 16*pf (near 0.5); the loose limits below pass it by a
  //  wide margin, while a real jump (~0.8*pf, i.e. ~8000*pf over one step) or kink blows past them.
  const double h = 1.0e-4;
  double prev_slope = 0.0;
  bool have_prev_slope = false;
  for( double r = 0.5; r <= 3.0; r += h )
  {
    const double v0 = P(r), v1 = P(r + h);
    BOOST_CHECK_MESSAGE( std::fabs(v1 - v0) < (20.0*pf*h + 1.0e-9),
                         "discontinuity near reldist=" << r << ": P jumped " << (v1 - v0)
                         << " over " << h );
    const double slope = (v1 - v0) / h;
    if( have_prev_slope )
      BOOST_CHECK_MESSAGE( std::fabs(slope - prev_slope) < (0.02*pf),
                           "slope kink near reldist=" << r << ": slope " << prev_slope
                           << " -> " << slope );
    prev_slope = slope;
    have_prev_slope = true;
  }

  // Pin the specific c=1.75 form (not merely "some continuous function").
  auto expected = [pf,cutoff]( const double r ){
    const double frac = r/cutoff, d = 1.0 - frac*frac;
    return (pf*d*d)/r;
  };
  for( const double r : { 0.5, 1.0, 1.25, 1.5, 1.7 } )
    BOOST_CHECK_CLOSE( P(r), expected(r), 1.0e-6 );

  // Hand-computed sanity values (per the design table).
  BOOST_CHECK_CLOSE( P(1.0),  0.453561*pf, 1.0e-2 );
  BOOST_CHECK_CLOSE( P(1.25), 0.191920*pf, 1.0e-2 );
}//BOOST_AUTO_TEST_CASE( proximityPunishmentContinuity )
