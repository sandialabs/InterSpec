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

/** Tests MakeDrfCalc::samplePoints - how "Make Detector Response" turns each peak of a source
 spectrum into an absolute-efficiency point: which gamma lines count toward which peak, the decay
 during the measurement, the source shielding and air, background subtraction, the guard against
 peaks an unusable source would add to, and which peaks the FWHM fit uses.

 Expected values come from SandiaDecay directly (aged mixtures, line by line), not from the code
 under test.
 */

#include "InterSpec_config.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <iostream>

#define BOOST_TEST_MODULE test_MakeDrfSamplePoints_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/EnergyCalibration.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/MakeDrfCalc.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/GammaInteractionCalc.h"

using namespace std;
using MakeDrfCalc::CalSource;
using MakeDrfCalc::PeakPoint;
using MakeDrfCalc::SamplePeak;
using MakeDrfCalc::SampleInput;
using MakeDrfCalc::SampleResult;
using MakeDrfCalc::BackgroundPeak;

namespace
{
void set_data_dir()
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
      if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.reactiongamma.xml") ) )
      {
        datadir = d;
        break;
      }
    }
  }//if( datadir.empty() )

  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( SpecUtils::append_path(datadir, "sandia.decay.xml") ),
                         "sandia.decay.xml not in '" << datadir << "'; pass --datadir=" );
  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );
  DecayDataBaseServer::setDecayXmlFile( SpecUtils::append_path( datadir, "sandia.decay.xml" ) );
}//set_data_dir()


const SandiaDecay::Nuclide *nuclide( const string &name )
{
  set_data_dir();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );
  const SandiaDecay::Nuclide * const nuc = db->nuclide( name );
  BOOST_REQUIRE_MESSAGE( nuc, "No nuclide " << name );
  return nuc;
}//nuclide(...)


/** A Gaussian peak; with a nuclide, assigned to its line nearest `energy` (and centered on it). */
shared_ptr<const PeakDef> make_peak( const SandiaDecay::Nuclide *nuc, const double energy,
                                     const double sigma, const double area, const double area_uncert )
{
  auto peak = make_shared<PeakDef>( energy, sigma, area );
  peak->setAmplitudeUncert( area_uncert );
  peak->setSigmaUncert( 0.01*sigma );

  if( nuc )
  {
    const SandiaDecay::Transition *transition = nullptr;
    size_t index = 0;
    PeakDef::SourceGammaType type = PeakDef::NormalGamma;
    PeakDef::findNearestPhotopeak( nuc, energy, 0.5, false, -1.0, transition, index, type );
    BOOST_REQUIRE_MESSAGE( transition || (type == PeakDef::AnnihilationGamma),
                           "No " << nuc->symbol << " line near " << energy << " keV" );
    peak->setNuclearTransition( nuc, transition, static_cast<int>(index), type );
    peak->setMean( peak->gammaParticleEnergy() );
  }//if( nuc )

  return peak;
}//make_peak(...)


/** Gammas per second `nuc` (activity `act` at age `age`) puts into peak `which` of `peaks`: the
 lines whose nearest peak, in units of the peaks' sigma, is `which`, and within 1.25 sigma of it.
 Optionally averaged over a measurement of `real_time` (for a nuclide with no progeny). */
double expected_rate( const SandiaDecay::Nuclide *nuc, const double act, const double age,
                      const vector<shared_ptr<const PeakDef>> &peaks, const size_t which,
                      const double real_time = 0.0 )
{
  SandiaDecay::NuclideMixture mix;
  mix.addAgedNuclideByActivity( nuc, act, age );

  double rate = 0.0;
  for( const SandiaDecay::EnergyRatePair &line : mix.photons( 0.0 ) )
  {
    double nearest = 1.0E9;
    size_t nearest_index = 0;
    for( size_t i = 0; i < peaks.size(); ++i )
    {
      const double nsigma = fabs( line.energy - peaks[i]->gammaParticleEnergy() ) / peaks[i]->sigma();
      if( nsigma < nearest )
      {
        nearest = nsigma;
        nearest_index = i;
      }
    }

    if( (nearest < MakeDrfCalc::sm_photopeakClusterSigma) && (nearest_index == which) )
      rate += line.numPerSecond;
  }//for( loop over lines )

  if( real_time > 0.0 )
  {
    const double lambda = log(2.0) / nuc->halfLife;
    rate *= (1.0 - exp(-lambda*real_time)) / (lambda*real_time);
  }

  return rate;
}//expected_rate(...)


/** Finds the result for peak `index`; fails the test if there is none. */
const PeakPoint &point_for( const SampleResult &result, const size_t index )
{
  for( const PeakPoint &p : result.points )
  {
    if( p.peak == index )
      return p;
  }
  BOOST_REQUIRE_MESSAGE( false, "No result for peak " << index );
  return result.points.front();
}//point_for(...)


bool has_point_for( const SampleResult &result, const size_t index )
{
  for( const PeakPoint &p : result.points )
  {
    if( p.peak == index )
      return true;
  }
  return false;
}//has_point_for(...)

}//namespace


BOOST_AUTO_TEST_CASE( single_source_efficiency )
{
  const SandiaDecay::Nuclide * const co60 = nuclide( "Co60" );

  CalSource src;
  src.nuclide = co60;
  src.activity = 1.0E4 * PhysicalUnits::bq;
  src.fracActivityUncert = 0.03;
  src.age = 1.0 * PhysicalUnits::year;
  src.distance = 25.0 * PhysicalUnits::cm;
  src.distanceUncert = 0.1 * PhysicalUnits::cm;

  SampleInput sample;
  sample.liveTime = 300.0;
  sample.realTime = 301.0;  //far too short for Co60 to decay during
  sample.sources.push_back( src );

  const vector<shared_ptr<const PeakDef>> peaks = { make_peak( co60, 1173.2, 1.0, 5000.0, 80.0 ),
                                                    make_peak( co60, 1332.5, 1.1, 4000.0, 70.0 ) };
  for( const auto &p : peaks )
  {
    SamplePeak sp;
    sp.peak = p;
    sp.source = 0;
    sample.peaks.push_back( sp );
  }

  const SampleResult result = MakeDrfCalc::samplePoints( sample, {}, false );
  BOOST_REQUIRE_EQUAL( result.points.size(), 2 );
  BOOST_CHECK( result.sourceErrors.empty() );

  for( size_t i = 0; i < peaks.size(); ++i )
  {
    const PeakPoint &pt = point_for( result, i );
    BOOST_REQUIRE( pt.status == PeakPoint::Status::Used );

    const double rate = expected_rate( co60, src.activity, src.age, peaks, i );
    BOOST_CHECK_GT( rate, 0.99*src.activity );  //both lines are ~100%
    BOOST_CHECK_CLOSE( pt.sourceRate, rate, 1.0E-6 );
    BOOST_CHECK_CLOSE( pt.fracSourceUncert, 0.03, 1.0E-6 );
    BOOST_CHECK_CLOSE( pt.energy, peaks[i]->gammaParticleEnergy(), 1.0E-6 );
    BOOST_CHECK( pt.useForFwhm );

    const MeasuredEffPoint &m = pt.measured;
    BOOST_CHECK_CLOSE( m.efficiency, peaks[i]->peakArea() / (rate * sample.liveTime), 1.0E-4 );
    BOOST_CHECK_CLOSE( m.fracStatUncert, peaks[i]->peakAreaUncert() / peaks[i]->peakArea(), 1.0E-4 );
    BOOST_CHECK_CLOSE( m.fracCertUncert, 0.03, 1.0E-4 );
    BOOST_CHECK_CLOSE( m.distance, src.distance, 1.0E-4 );
    BOOST_CHECK_CLOSE( m.distanceUncert, src.distanceUncert, 1.0E-4 );
    BOOST_CHECK_CLOSE( m.liveTime, sample.liveTime, 1.0E-4 );
    BOOST_CHECK_CLOSE( m.peakArea, peaks[i]->peakArea(), 1.0E-4 );
    BOOST_CHECK_LT( m.bkgPeakArea, 0.0f );  //no background subtracted
  }//for( size_t i = 0; i < peaks.size(); ++i )

  // A fixed-geometry source has no distance; the point carries none
  sample.sources[0].distance = -1.0;
  sample.sources[0].distanceUncert = 0.0;
  const SampleResult fixed = MakeDrfCalc::samplePoints( sample, {}, true );
  BOOST_REQUIRE_EQUAL( fixed.points.size(), 2 );
  BOOST_CHECK_LT( fixed.points[0].measured.distance, 0.0f );
  BOOST_CHECK_EQUAL( fixed.points[0].measured.distanceUncert, 0.0f );
  BOOST_CHECK_CLOSE( fixed.points[0].sourceRate, result.points[0].sourceRate, 1.0E-6 );  //no air path
}//BOOST_AUTO_TEST_CASE( single_source_efficiency )


BOOST_AUTO_TEST_CASE( lines_go_to_nearest_peak )
{
  // Ba133's 79.61 and 81.00 keV lines are 1.38 keV apart
  const SandiaDecay::Nuclide * const ba133 = nuclide( "Ba133" );

  CalSource src;
  src.nuclide = ba133;
  src.activity = 5.0E4 * PhysicalUnits::bq;
  src.age = 2.0 * PhysicalUnits::year;
  src.distance = 25.0 * PhysicalUnits::cm;

  SampleInput sample;
  sample.liveTime = sample.realTime = 600.0;
  sample.sources.push_back( src );

  // Resolved: each peak gets only its own line, though each line is within 1.25 sigma of both
  //  peaks - a line between two peaks counts once, toward the nearer.
  {
    const vector<shared_ptr<const PeakDef>> peaks = { make_peak( ba133, 79.6, 1.2, 1000.0, 40.0 ),
                                                      make_peak( ba133, 81.0, 1.2, 9000.0, 100.0 ) };
    sample.peaks.clear();
    for( const auto &p : peaks )
      sample.peaks.push_back( SamplePeak{ p, 0, -1.0 } );

    const SampleResult result = MakeDrfCalc::samplePoints( sample, {}, false );
    BOOST_REQUIRE_EQUAL( result.points.size(), 2 );

    const double rate_79 = expected_rate( ba133, src.activity, src.age, peaks, 0 );
    const double rate_81 = expected_rate( ba133, src.activity, src.age, peaks, 1 );
    BOOST_CHECK_CLOSE( point_for( result, 0 ).sourceRate, rate_79, 1.0E-6 );
    BOOST_CHECK_CLOSE( point_for( result, 1 ).sourceRate, rate_81, 1.0E-6 );

    // and the 81 keV line alone is ~12x the 79.6 keV one
    BOOST_CHECK_GT( rate_81, 10.0*rate_79 );

    // The sum is all of the source's gammas in the two peaks: nothing counted twice
    const vector<shared_ptr<const PeakDef>> merged = { make_peak( ba133, 81.0, 1.2, 0.0, 0.0 ) };
    BOOST_CHECK_CLOSE( rate_79 + rate_81, expected_rate( ba133, src.activity, src.age, merged, 0 ), 1.0E-6 );
  }

  // Only the 81 keV peak is fit, but the 79.6 keV line is within 1.25 sigma of it: both count.
  //  A peak that is fit but not used for the efficiency still takes its lines.
  {
    const vector<shared_ptr<const PeakDef>> peaks = { make_peak( ba133, 81.0, 1.2, 9000.0, 100.0 ) };
    sample.peaks.clear();
    sample.peaks.push_back( SamplePeak{ peaks[0], 0, -1.0 } );

    const SampleResult result = MakeDrfCalc::samplePoints( sample, {}, false );
    BOOST_REQUIRE_EQUAL( result.points.size(), 1 );
    const double rate = expected_rate( ba133, src.activity, src.age, peaks, 0 );
    BOOST_CHECK_CLOSE( result.points[0].sourceRate, rate, 1.0E-6 );

    sample.peaks.push_back( SamplePeak{ make_peak( ba133, 79.6, 1.2, 1000.0, 40.0 ), -1, -1.0 } );
    const SampleResult with_unused = MakeDrfCalc::samplePoints( sample, {}, false );
    BOOST_REQUIRE_EQUAL( with_unused.points.size(), 1 );
    BOOST_CHECK_LT( with_unused.points[0].sourceRate, 0.95*rate );
  }
}//BOOST_AUTO_TEST_CASE( lines_go_to_nearest_peak )


BOOST_AUTO_TEST_CASE( decay_during_measurement )
{
  // Na24 (15 h) has no progeny gammas, so the average rate over the measurement is analytic
  const SandiaDecay::Nuclide * const na24 = nuclide( "Na24" );

  CalSource src;
  src.nuclide = na24;
  src.activity = 1.0E5 * PhysicalUnits::bq;
  src.age = 0.0;
  src.distance = 50.0 * PhysicalUnits::cm;

  const vector<shared_ptr<const PeakDef>> peaks = { make_peak( na24, 1368.6, 1.0, 2.0E5, 500.0 ),
                                                    make_peak( na24, 2754.0, 1.5, 1.0E5, 350.0 ) };
  SampleInput sample;
  sample.sources.push_back( src );
  for( const auto &p : peaks )
    sample.peaks.push_back( SamplePeak{ p, 0, -1.0 } );

  // A fifth of a half-life: corrected, the average rate is ~6.7% below the starting one
  sample.realTime = 0.2 * na24->halfLife;
  sample.liveTime = 0.98 * sample.realTime;
  const SampleResult long_meas = MakeDrfCalc::samplePoints( sample, {}, false );
  BOOST_REQUIRE_EQUAL( long_meas.points.size(), 2 );
  for( size_t i = 0; i < peaks.size(); ++i )
  {
    const double start_rate = expected_rate( na24, src.activity, src.age, peaks, i );
    const double avrg_rate = expected_rate( na24, src.activity, src.age, peaks, i, sample.realTime );
    BOOST_CHECK_LT( avrg_rate, 0.94*start_rate );
    BOOST_CHECK_CLOSE( point_for( long_meas, i ).sourceRate, avrg_rate, 1.0E-3 );
    BOOST_CHECK_CLOSE( point_for( long_meas, i ).measured.efficiency,
                       peaks[i]->peakArea() / (avrg_rate * sample.liveTime), 1.0E-3 );
  }

  // Below 0.5% of the half-life it is not corrected (the Activity/Shielding fit's default)
  sample.realTime = 0.004 * na24->halfLife;
  sample.liveTime = sample.realTime;
  const SampleResult short_meas = MakeDrfCalc::samplePoints( sample, {}, false );
  BOOST_REQUIRE_EQUAL( short_meas.points.size(), 2 );
  for( size_t i = 0; i < peaks.size(); ++i )
  {
    const double start_rate = expected_rate( na24, src.activity, src.age, peaks, i );
    BOOST_CHECK_CLOSE( point_for( short_meas, i ).sourceRate, start_rate, 1.0E-6 );
  }
}//BOOST_AUTO_TEST_CASE( decay_during_measurement )


BOOST_AUTO_TEST_CASE( shielding_and_air )
{
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );
  const shared_ptr<const Material> iron = matdb->material( "Fe (iron)" );
  BOOST_REQUIRE( iron );

  const SandiaDecay::Nuclide * const cs137 = nuclide( "Cs137" );

  CalSource bare;
  bare.nuclide = cs137;
  bare.activity = 3.0E4 * PhysicalUnits::bq;
  bare.age = 5.0 * PhysicalUnits::year;
  bare.distance = 200.0 * PhysicalUnits::cm;

  const vector<shared_ptr<const PeakDef>> peaks = { make_peak( cs137, 661.7, 0.9, 20000.0, 150.0 ) };
  SampleInput sample;
  sample.liveTime = sample.realTime = 300.0;
  sample.peaks.push_back( SamplePeak{ peaks[0], 0, -1.0 } );

  auto rate_for = [&sample]( const CalSource &src, const bool air ) -> double {
    sample.sources = { src };
    const SampleResult result = MakeDrfCalc::samplePoints( sample, {}, air );
    BOOST_REQUIRE_EQUAL( result.points.size(), 1 );
    BOOST_REQUIRE( result.points[0].status == PeakPoint::Status::Used );
    return result.points[0].sourceRate;
  };

  const double energy = peaks[0]->gammaParticleEnergy();
  const double emitted = expected_rate( cs137, bare.activity, bare.age, peaks, 0 );
  BOOST_CHECK_CLOSE( rate_for( bare, false ), emitted, 1.0E-6 );

  // Two meters of air: ~2% at 662 keV
  const double air_trans = exp( -GammaInteractionCalc::transmission_length_coefficient_air( energy ) * bare.distance );
  BOOST_CHECK( (air_trans > 0.97) && (air_trans < 0.99) );
  BOOST_CHECK_CLOSE( rate_for( bare, true ), emitted * air_trans, 1.0E-6 );

  // Generic shielding: 7.874 g/cm2 of Z=26
  CalSource generic = bare;
  generic.shieldAtomicNumber = 26.0;
  generic.shieldArealDensity = 7.874 * PhysicalUnits::g / PhysicalUnits::cm2;
  const double generic_trans = exp( -MassAttenuation::massAttenuationCoefficientFracAN( 26.0, energy )
                                    * generic.shieldArealDensity );
  BOOST_CHECK( (generic_trans > 0.5) && (generic_trans < 0.6) );  //~0.56 at 662 keV
  BOOST_CHECK_CLOSE( rate_for( generic, false ), emitted * generic_trans, 1.0E-6 );

  // 1 cm of iron, by its composition: the same areal density of (almost only) Z=26
  CalSource named = bare;
  named.shieldMaterial = iron;
  named.shieldThickness = 1.0 * PhysicalUnits::cm;
  BOOST_CHECK_CLOSE( rate_for( named, false ), emitted * generic_trans, 1.0 );

  // An atomic number below 1, or no areal density, is no shielding
  CalSource none = bare;
  none.shieldAtomicNumber = 0.0;
  none.shieldArealDensity = 5.0 * PhysicalUnits::g / PhysicalUnits::cm2;
  BOOST_CHECK_CLOSE( rate_for( none, false ), emitted, 1.0E-6 );
}//BOOST_AUTO_TEST_CASE( shielding_and_air )


BOOST_AUTO_TEST_CASE( other_sources_and_unusable_ones )
{
  // Co57's 122.06 keV line and Eu152's 121.78 keV line land in one peak
  const SandiaDecay::Nuclide * const co57 = nuclide( "Co57" );
  const SandiaDecay::Nuclide * const eu152 = nuclide( "Eu152" );

  CalSource co57_src;
  co57_src.nuclide = co57;
  co57_src.activity = 2.0E4 * PhysicalUnits::bq;
  co57_src.fracActivityUncert = 0.02;
  co57_src.age = 0.5 * PhysicalUnits::year;
  co57_src.distance = 25.0 * PhysicalUnits::cm;

  CalSource eu152_src = co57_src;
  eu152_src.nuclide = eu152;
  eu152_src.activity = 1.0E4 * PhysicalUnits::bq;
  eu152_src.fracActivityUncert = 0.05;
  eu152_src.age = 3.0 * PhysicalUnits::year;

  const vector<shared_ptr<const PeakDef>> peaks = { make_peak( co57, 122.06, 0.5, 30000.0, 200.0 ) };

  SampleInput sample;
  sample.liveTime = sample.realTime = 600.0;
  sample.sources = { co57_src, eu152_src };  //Eu152 has no peaks used, but its lines still count
  sample.peaks.push_back( SamplePeak{ peaks[0], 0, -1.0 } );

  const double co57_rate = expected_rate( co57, co57_src.activity, co57_src.age, peaks, 0 );
  const double eu152_rate = expected_rate( eu152, eu152_src.activity, eu152_src.age, peaks, 0 );
  BOOST_REQUIRE_GT( eu152_rate, 0.1*co57_rate );

  // Both usable: the peak's gammas are from both, and their activity uncertainties add
  {
    const SampleResult result = MakeDrfCalc::samplePoints( sample, {}, false );
    BOOST_REQUIRE_EQUAL( result.points.size(), 1 );
    const PeakPoint &pt = result.points[0];
    BOOST_REQUIRE( pt.status == PeakPoint::Status::Used );
    BOOST_CHECK_CLOSE( pt.sourceRate, co57_rate + eu152_rate, 1.0E-6 );
    const double uncert = sqrt( pow(0.02*co57_rate, 2.0) + pow(0.05*eu152_rate, 2.0) );
    BOOST_CHECK_CLOSE( pt.fracSourceUncert, uncert / (co57_rate + eu152_rate), 1.0E-6 );
    BOOST_CHECK_CLOSE( pt.measured.efficiency,
                       peaks[0]->peakArea() / ((co57_rate + eu152_rate) * sample.liveTime), 1.0E-4 );
  }

  // Eu152 unusable: its contribution can't be known, so the peak is left out, saying why
  {
    sample.sources[1].usable = false;
    const SampleResult result = MakeDrfCalc::samplePoints( sample, {}, false );
    BOOST_REQUIRE_EQUAL( result.points.size(), 1 );
    BOOST_CHECK( result.points[0].status == PeakPoint::Status::NeedsUnusableSource );
    BOOST_CHECK_EQUAL( result.points[0].interferer, 1 );
  }

  // ... unless the peak is narrow enough that the Eu152 line is outside it
  {
    const vector<shared_ptr<const PeakDef>> narrow = { make_peak( co57, 122.06, 0.1, 30000.0, 200.0 ) };
    sample.peaks[0].peak = narrow[0];
    const SampleResult result = MakeDrfCalc::samplePoints( sample, {}, false );
    BOOST_REQUIRE_EQUAL( result.points.size(), 1 );
    BOOST_REQUIRE( result.points[0].status == PeakPoint::Status::Used );
    BOOST_CHECK_CLOSE( result.points[0].sourceRate,
                       expected_rate( co57, co57_src.activity, co57_src.age, narrow, 0 ), 1.0E-6 );
  }

  // An unusable source's own peaks give no result at all (its problems say why)
  {
    sample.sources[0].usable = false;
    sample.sources[1].usable = true;
    const SampleResult result = MakeDrfCalc::samplePoints( sample, {}, false );
    BOOST_CHECK( result.points.empty() );
  }
}//BOOST_AUTO_TEST_CASE( other_sources_and_unusable_ones )


BOOST_AUTO_TEST_CASE( background_subtraction )
{
  const SandiaDecay::Nuclide * const cs137 = nuclide( "Cs137" );

  CalSource src;
  src.nuclide = cs137;
  src.activity = 1.0E4 * PhysicalUnits::bq;
  src.age = 5.0 * PhysicalUnits::year;
  src.distance = 25.0 * PhysicalUnits::cm;

  const vector<shared_ptr<const PeakDef>> peaks = { make_peak( cs137, 661.7, 1.0, 10000.0, 100.0 ) };
  SampleInput sample;
  sample.liveTime = sample.realTime = 300.0;
  sample.sources = { src };
  sample.peaks.push_back( SamplePeak{ peaks[0], 0, -1.0 } );

  // Background spectrum 0 (1000 s) has two peaks that match - its live time counts once; spectrum
  //  1 (500 s) has one; a far-away peak doesn't match.
  vector<BackgroundPeak> background = {
    BackgroundPeak{ make_peak( nullptr, 661.0, 1.0, 100.0, 10.0 ), 1000.0, 0 },
    BackgroundPeak{ make_peak( nullptr, 662.0, 1.0, 200.0, 15.0 ), 1000.0, 0 },
    BackgroundPeak{ make_peak( nullptr, 1460.8, 1.5, 300.0, 20.0 ), 1000.0, 0 },
    BackgroundPeak{ make_peak( nullptr, 661.7, 1.0, 150.0, 12.0 ), 500.0, 1 }
  };

  const SampleResult result = MakeDrfCalc::samplePoints( sample, background, false );
  BOOST_REQUIRE_EQUAL( result.points.size(), 1 );
  const PeakPoint &pt = result.points[0];
  BOOST_REQUIRE( pt.status == PeakPoint::Status::Used );

  BOOST_CHECK_CLOSE( pt.bkgArea, 450.0, 1.0E-6 );
  BOOST_CHECK_CLOSE( pt.bkgLiveTime, 1500.0, 1.0E-6 );

  const double scaled_back = 300.0 * 450.0 / 1500.0;  //90 counts
  const double scaled_back_uncert = scaled_back * sqrt( 10.0*10.0 + 15.0*15.0 + 12.0*12.0 ) / 450.0;
  BOOST_CHECK_CLOSE( pt.peakArea, 10000.0 - scaled_back, 1.0E-6 );
  BOOST_CHECK_CLOSE( pt.peakAreaUncert, sqrt( 100.0*100.0 + scaled_back_uncert*scaled_back_uncert ), 1.0E-6 );
  BOOST_CHECK_CLOSE( pt.measured.bkgPeakArea, scaled_back, 1.0E-4 );
  BOOST_CHECK_CLOSE( pt.measured.bkgPeakAreaUncert, scaled_back_uncert, 1.0E-4 );
  BOOST_CHECK_CLOSE( pt.measured.efficiency, (10000.0 - scaled_back) / (pt.sourceRate * 300.0), 1.0E-4 );

  // Background bigger than the peak: no efficiency, and it says so (not "no source gammas")
  sample.peaks[0].peak = make_peak( cs137, 661.7, 1.0, 50.0, 10.0 );
  const SampleResult swamped = MakeDrfCalc::samplePoints( sample, background, false );
  BOOST_REQUIRE_EQUAL( swamped.points.size(), 1 );
  BOOST_CHECK( swamped.points[0].status == PeakPoint::Status::NoNetArea );
}//BOOST_AUTO_TEST_CASE( background_subtraction )


BOOST_AUTO_TEST_CASE( user_branching_ratio_and_no_gammas )
{
  const SandiaDecay::Nuclide * const cs137 = nuclide( "Cs137" );

  // A source with no nuclide: the gammas are the user's branching ratio times the activity
  CalSource user_src;
  user_src.activity = 2.0E4 * PhysicalUnits::bq;
  user_src.distance = 25.0 * PhysicalUnits::cm;

  SampleInput sample;
  sample.liveTime = sample.realTime = 100.0;
  sample.sources = { user_src };
  sample.peaks.push_back( SamplePeak{ make_peak( nullptr, 1000.0, 1.2, 5000.0, 80.0 ), 0, 0.5 } );

  {
    const SampleResult result = MakeDrfCalc::samplePoints( sample, {}, false );
    BOOST_REQUIRE_EQUAL( result.points.size(), 1 );
    BOOST_REQUIRE( result.points[0].status == PeakPoint::Status::Used );
    BOOST_CHECK_CLOSE( result.points[0].sourceRate, 1.0E4, 1.0E-6 );
    BOOST_CHECK_CLOSE( result.points[0].energy, 1000.0, 1.0E-6 );
  }

  // No (valid) branching ratio entered
  sample.peaks[0].userBr = -1.0;
  {
    const SampleResult result = MakeDrfCalc::samplePoints( sample, {}, false );
    BOOST_REQUIRE_EQUAL( result.points.size(), 1 );
    BOOST_CHECK( result.points[0].status == PeakPoint::Status::NoUserBranchRatio );
  }

  // A nuclide with no line near the peak, or no activity, puts no gammas in it
  CalSource cs_src;
  cs_src.nuclide = cs137;
  cs_src.activity = 1.0E4 * PhysicalUnits::bq;
  cs_src.age = 5.0 * PhysicalUnits::year;
  sample.sources = { cs_src };
  sample.peaks[0] = SamplePeak{ make_peak( cs137, 661.7, 0.9, 5000.0, 80.0 ), 0, -1.0 };
  sample.peaks.push_back( SamplePeak{ make_peak( nullptr, 150.0, 0.7, 500.0, 30.0 ), 0, -1.0 } );
  {
    const SampleResult result = MakeDrfCalc::samplePoints( sample, {}, false );
    BOOST_REQUIRE_EQUAL( result.points.size(), 2 );
    BOOST_CHECK( point_for( result, 0 ).status == PeakPoint::Status::Used );
    BOOST_CHECK( point_for( result, 1 ).status == PeakPoint::Status::NoSourceGammas );
    BOOST_CHECK( point_for( result, 1 ).useForFwhm );  //its width is still good
  }

  sample.sources[0].activity = 0.0;
  {
    const SampleResult result = MakeDrfCalc::samplePoints( sample, {}, false );
    BOOST_REQUIRE_EQUAL( result.points.size(), 2 );
    BOOST_CHECK( point_for( result, 0 ).status == PeakPoint::Status::NoSourceGammas );
  }
}//BOOST_AUTO_TEST_CASE( user_branching_ratio_and_no_gammas )


BOOST_AUTO_TEST_CASE( annihilation_left_out_of_fwhm )
{
  const SandiaDecay::Nuclide * const na22 = nuclide( "Na22" );

  CalSource src;
  src.nuclide = na22;
  src.activity = 1.0E4 * PhysicalUnits::bq;
  src.age = 1.0 * PhysicalUnits::year;
  src.distance = 25.0 * PhysicalUnits::cm;

  const vector<shared_ptr<const PeakDef>> peaks = { make_peak( na22, 511.0, 1.6, 8000.0, 100.0 ),
                                                    make_peak( na22, 1274.5, 1.1, 3000.0, 60.0 ) };
  BOOST_CHECK( peaks[0]->sourceGammaType() == PeakDef::AnnihilationGamma );

  SampleInput sample;
  sample.liveTime = sample.realTime = 600.0;
  sample.sources = { src };
  for( const auto &p : peaks )
    sample.peaks.push_back( SamplePeak{ p, 0, -1.0 } );

  const SampleResult result = MakeDrfCalc::samplePoints( sample, {}, false );
  BOOST_REQUIRE_EQUAL( result.points.size(), 2 );

  // The 511 keV peak still gives an efficiency (two photons per positron), just not a width
  const PeakPoint &annih = point_for( result, 0 );
  BOOST_CHECK( annih.status == PeakPoint::Status::Used );
  BOOST_CHECK( !annih.useForFwhm );
  BOOST_CHECK_GT( annih.sourceRate, 1.7*src.activity );
  BOOST_CHECK_CLOSE( annih.sourceRate, expected_rate( na22, src.activity, src.age, peaks, 0 ), 1.0E-6 );

  BOOST_CHECK( point_for( result, 1 ).useForFwhm );

  // A peak near 511 keV with no source assigned is treated as annihilation too
  sample.sources = { CalSource{} };
  sample.sources[0].activity = 1.0E4 * PhysicalUnits::bq;
  sample.peaks = { SamplePeak{ make_peak( nullptr, 511.2, 1.6, 8000.0, 100.0 ), 0, 1.8 } };
  const SampleResult user = MakeDrfCalc::samplePoints( sample, {}, false );
  BOOST_REQUIRE_EQUAL( user.points.size(), 1 );
  BOOST_CHECK( !user.points[0].useForFwhm );
}//BOOST_AUTO_TEST_CASE( annihilation_left_out_of_fwhm )


BOOST_AUTO_TEST_CASE( data_defined_peak )
{
  // A data-defined peak has no sigma (PeakDef::sigma() throws); a quarter of its ROI stands in for
  //  it - including when matching background peaks.
  const SandiaDecay::Nuclide * const cs137 = nuclide( "Cs137" );

  // 3 keV channels, flat at 10 counts, with ~6000 counts over 657-669 keV
  auto counts = make_shared<vector<float>>( 1024, 10.0f );
  for( size_t i = 219; i <= 222; ++i )
    (*counts)[i] += 1500.0f;
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( counts->size(), { 0.0f, 3.0f }, {} );
  auto meas = make_shared<SpecUtils::Measurement>();
  meas->set_gamma_counts( counts, 300.0f, 300.0f );
  meas->set_energy_calibration( cal );

  auto peak = make_shared<PeakDef>( 649.0, 673.0, 661.657, meas, nullptr );
  BOOST_REQUIRE( !peak->gausPeak() );
  {
    // Assign it to Cs137's 661.657 keV line
    const SandiaDecay::Transition *transition = nullptr;
    size_t index = 0;
    PeakDef::SourceGammaType type = PeakDef::NormalGamma;
    PeakDef::findNearestPhotopeak( cs137, 661.657, 0.5, false, -1.0, transition, index, type );
    BOOST_REQUIRE( transition );
    peak->setNuclearTransition( cs137, transition, static_cast<int>(index), type );
  }
  BOOST_REQUIRE_GT( peak->peakArea(), 1000.0 );

  CalSource src;
  src.nuclide = cs137;
  src.activity = 1.0E4 * PhysicalUnits::bq;
  src.age = 5.0 * PhysicalUnits::year;
  src.distance = 25.0 * PhysicalUnits::cm;

  SampleInput sample;
  sample.liveTime = sample.realTime = 300.0;
  sample.sources = { src };
  sample.peaks.push_back( SamplePeak{ peak, 0, -1.0 } );

  // Within 1.5 x (24 keV ROI / 4) of the peak mean
  const vector<BackgroundPeak> background = { BackgroundPeak{ make_peak( nullptr, 664.0, 1.0, 100.0, 10.0 ), 1000.0, 0 } };

  SampleResult result;
  BOOST_REQUIRE_NO_THROW( result = MakeDrfCalc::samplePoints( sample, background, false ) );
  BOOST_REQUIRE_EQUAL( result.points.size(), 1 );
  const PeakPoint &pt = result.points[0];
  BOOST_CHECK( pt.status == PeakPoint::Status::Used );
  BOOST_CHECK( !pt.useForFwhm );
  BOOST_CHECK_CLOSE( pt.bkgArea, 100.0, 1.0E-6 );
  BOOST_CHECK_CLOSE( pt.peakArea, peak->peakArea() - 30.0, 1.0E-6 );

  const vector<shared_ptr<const PeakDef>> as_gaus = { make_peak( cs137, 661.7, 0.25*peak->roiWidth(), 1.0, 1.0 ) };
  BOOST_CHECK_CLOSE( pt.sourceRate, expected_rate( cs137, src.activity, src.age, as_gaus, 0 ), 1.0E-6 );
}//BOOST_AUTO_TEST_CASE( data_defined_peak )


BOOST_AUTO_TEST_CASE( default_equation_terms )
{
  // Two points more than coefficients, so chi2/dof means something; one more with three points
  BOOST_CHECK_EQUAL( MakeDrfCalc::defaultEquationTerms( 0, 7 ), 0 );
  BOOST_CHECK_EQUAL( MakeDrfCalc::defaultEquationTerms( 1, 7 ), 1 );
  BOOST_CHECK_EQUAL( MakeDrfCalc::defaultEquationTerms( 2, 7 ), 2 );
  BOOST_CHECK_EQUAL( MakeDrfCalc::defaultEquationTerms( 3, 7 ), 2 );
  BOOST_CHECK_EQUAL( MakeDrfCalc::defaultEquationTerms( 4, 7 ), 2 );
  BOOST_CHECK_EQUAL( MakeDrfCalc::defaultEquationTerms( 5, 7 ), 3 );
  BOOST_CHECK_EQUAL( MakeDrfCalc::defaultEquationTerms( 9, 7 ), 7 );
  BOOST_CHECK_EQUAL( MakeDrfCalc::defaultEquationTerms( 15, 7 ), 7 );
  BOOST_CHECK_EQUAL( MakeDrfCalc::defaultEquationTerms( 6, 5 ), 4 );
  BOOST_CHECK_EQUAL( MakeDrfCalc::defaultEquationTerms( 12, 5 ), 5 );
}//BOOST_AUTO_TEST_CASE( default_equation_terms )
