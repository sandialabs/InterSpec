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

/* Tests of the Flux Tool's per-peak math (FluxToolCalc), and the background-peak helpers it shares
 with the Activity/Shielding fit.  The tool has had two serious bugs - flux per mm2 instead of per
 cm2, and dividing by the live time of the whole file rather than of the displayed spectrum - so
 these pin the numbers against independent, by-hand, expectations.
 */

#define _USE_MATH_DEFINES

#include "InterSpec_config.h"

#define BOOST_TEST_MODULE test_FluxTool_suite
#include <boost/test/included/unit_test.hpp>

#if ( defined( WIN32 ) )
#undef min
#undef max
#endif

#include <array>
#include <cmath>
#include <deque>
#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include <iostream>

#include "io/DetectorResponse.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/FluxToolCalc.h"
#include "InterSpec/DrfModifyCalc.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/PeakSearchGuiUtils.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;

namespace
{
  string g_data_dir;
  string g_test_data_dir;

  typedef DetectorPeakResponse::EffGeometryType GeomType;


  /** A displayed-histogram stand-in with the given live time. */
  shared_ptr<SpecUtils::Measurement> make_measurement( const float live_time )
  {
    auto meas = make_shared<SpecUtils::Measurement>();
    auto counts = make_shared<vector<float>>( 1024, 10.0f );
    meas->set_gamma_counts( counts, live_time, 1.1f*live_time );
    auto cal = make_shared<SpecUtils::EnergyCalibration>();
    cal->set_polynomial( 1024, { 0.0f, 3.0f }, {} );
    meas->set_energy_calibration( cal );
    return meas;
  }//make_measurement(...)


  /** A Gaussian peak with the given area, and Poisson area uncertainty. */
  PeakDef make_peak( const double mean, const double sigma, const double area )
  {
    PeakDef peak( mean, sigma, area );
    peak.setPeakAreaUncert( std::sqrt( area ) );
    return peak;
  }


  /** The fraction of 4pi a disk of diameter `diameter` subtends from a point `dist` away on its axis -
   written out here, rather than calling DetectorPeakResponse::fractionalSolidAngle, so the test does
   not lean on the code under test.
   */
  double disk_fraction( const double diameter, const double dist )
  {
    const double r = 0.5*diameter;
    return 0.5*(1.0 - dist/std::sqrt(dist*dist + r*r));
  }


  /** A far-field DRF with a constant intrinsic efficiency, a known diameter, and a setback. */
  shared_ptr<DetectorPeakResponse> make_far_field_drf( const double intrinsic, const double diameter,
                                                       const double setback, const float lowerE = 0.0f,
                                                       const float upperE = 0.0f )
  {
    auto drf = make_shared<DetectorPeakResponse>( "FarFieldTest", "constant intrinsic efficiency" );
    drf->setIntrinsicEfficiencyFormula( std::to_string(intrinsic), static_cast<float>(diameter),
                                        PhysicalUnits::keV, lowerE, upperE, GeomType::FarFieldIntrinsic );
    // Must come after the efficiency, which resets it
    if( setback > 0.0 )
      drf->setDetectorSetback( setback );
    drf->setFwhmCoefficients( { 1.0f, 0.03f, 0.0f }, DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn );
    BOOST_REQUIRE( drf->isValid() );
    return drf;
  }//make_far_field_drf(...)


  /** A far-field DRF that states its own efficiency uncertainty (per-point), so the systematic parts
   are non-zero.
   */
  shared_ptr<DetectorPeakResponse> make_drf_with_uncert()
  {
    vector<DetectorPeakResponse::EnergyEffPoint> points;
    for( const array<float,3> &eeu : vector<array<float,3>>{ {40.0f,0.55f,0.03f}, {122.0f,0.62f,0.02f},
                         {344.0f,0.40f,0.02f}, {662.0f,0.26f,0.015f}, {1332.0f,0.16f,0.012f}, {3000.0f,0.08f,0.01f} } )
    {
      DetectorPeakResponse::EnergyEffPoint p;
      p.energy = eeu[0];
      p.efficiency = eeu[1];
      p.efficiencyUncert = eeu[1] * eeu[2];
      points.push_back( p );
    }

    auto drf = make_shared<DetectorPeakResponse>( "WithUncert", "energy-efficiency pairs, with uncertainties" );
    drf->setEfficiencyPoints( points, 7.62f*PhysicalUnits::cm, 0.0, GeomType::FarFieldIntrinsic );
    drf->setDetectorSetback( 1.0*PhysicalUnits::cm );
    BOOST_REQUIRE( drf->isValid() );
    BOOST_REQUIRE( drf->statesOwnEfficiencyUncert() );
    return drf;
  }//make_drf_with_uncert()


  /** The Flux Tool's per-peak math exactly as it was written in FluxToolWidget::refreshPeakTable() before
   it was moved to FluxToolCalc::peak_flux - kept here so we notice if a refactor changes any number,
   even in the last bit.
   */
  struct LegacyFlux
  {
    double cps, cpsUncert, geomEff, intrinsic, intrinsicSyst;
    double fluxOnDet, fluxOnDetStat, fluxOnDetSyst;
    double flux, fluxStat, fluxSyst;
    double gammaInto4pi, gammaInto4piStat, gammaInto4piSyst;
  };

  LegacyFlux legacy_flux( const PeakDef &peak, const shared_ptr<const DetectorPeakResponse> &det,
                          const float live_time, const double distance )
  {
    LegacyFlux answer{};
    const bool fixed_geom = det->isFixedGeometry();
    const double energy = peak.mean();

    const double amp = peak.peakArea();
    const double ampUncert = peak.peakAreaUncert();
    const double cps = amp / live_time;
    const double cpsUncert = ampUncert / live_time;
    const double geomEff = fixed_geom ? 1.0 : det->fractionalSolidAngle( det->detectorDiameter(), distance + det->detectorSetback() );
    const double totaleff = det->efficiency( energy, distance );  //distance ignored for fixed geometry
    const double intrinsic = (geomEff > 0.0) ? (totaleff / geomEff) : 0.0;
    const DrfModifyCalc::UncertSummary drf_uncert = DrfModifyCalc::uncertSummary( *det, energy, distance );
    const double eff_uncert = drf_uncert.usedFrac( false );

    answer.cps = cps;
    answer.cpsUncert = cpsUncert;
    answer.geomEff = geomEff;
    answer.intrinsic = intrinsic;
    answer.intrinsicSyst = intrinsic * eff_uncert;

    if( totaleff <= 0.0 || intrinsic <= 0.0 )
    {
      answer.fluxOnDet = answer.flux = answer.gammaInto4pi = std::numeric_limits<double>::infinity();
    }else
    {
      const double fluxOnDet = cps / intrinsic;
      const double gammaInto4pi = cps / totaleff;
      const double distance_cm = distance / PhysicalUnits::cm;
      const double area_4pi = 4*M_PI*distance_cm*distance_cm;
      const double flux = gammaInto4pi / area_4pi;

      answer.fluxOnDet = fluxOnDet;
      answer.fluxOnDetStat = cpsUncert / intrinsic;
      answer.fluxOnDetSyst = fluxOnDet * eff_uncert;
      answer.flux = flux;
      answer.fluxStat = cpsUncert / totaleff / area_4pi;
      answer.fluxSyst = flux * eff_uncert;
      answer.gammaInto4pi = gammaInto4pi;
      answer.gammaInto4piStat = cpsUncert / totaleff;
      answer.gammaInto4piSyst = gammaInto4pi * eff_uncert;
    }

    return answer;
  }//legacy_flux(...)


  string read_file( const string &path )
  {
    ifstream input( path.c_str(), ios::in | ios::binary );
    BOOST_REQUIRE_MESSAGE( input.is_open(), "Could not open '" + path + "'" );
    stringstream strm;
    strm << input.rdbuf();
    return strm.str();
  }
}//namespace


struct GlobalFixture
{
  GlobalFixture()
  {
    const int argc = boost::unit_test::framework::master_test_suite().argc;
    char **argv = boost::unit_test::framework::master_test_suite().argv;

    for( int i = 1; i < argc; ++i )
    {
      const string arg = argv[i];
      if( arg.find("--datadir=") == 0 )
        g_data_dir = arg.substr( 10 );
      else if( arg.find("--testfiledir=") == 0 )
        g_test_data_dir = arg.substr( 14 );
    }

    // Search some likely places, if not given
    for( const char *dir : { "data", "../data", "../../data", "../../../data" } )
    {
      if( g_data_dir.empty() && SpecUtils::is_file( SpecUtils::append_path(dir, "sandia.decay.xml") ) )
        g_data_dir = dir;
    }

    if( !g_data_dir.empty() )
    {
      try
      {
        InterSpec::setStaticDataDirectory( g_data_dir );
      }catch( std::exception &e )
      {
        cerr << "Failed to set static data directory: " << e.what() << endl;
      }
    }
  }//GlobalFixture()
};//struct GlobalFixture

BOOST_GLOBAL_FIXTURE( GlobalFixture );


// A point source of known emission rate, at a known distance from a detector of known diameter,
//  setback, and intrinsic efficiency: the tool must recover the emission rate and fluxes.
BOOST_AUTO_TEST_CASE( far_field_round_trip )
{
  const double diameter = 5.0*PhysicalUnits::cm;
  const double setback = 0.5*PhysicalUnits::cm;
  const double intrinsic = 0.25;
  const shared_ptr<DetectorPeakResponse> drf = make_far_field_drf( intrinsic, diameter, setback );

  const double emission_rate = 1.0E5;  //gammas/s, into 4pi
  const float live_time = 300.0f;

  for( const double distance : { 1.0*PhysicalUnits::cm, 50.0*PhysicalUnits::cm, 100.0*PhysicalUnits::cm, 10.0*PhysicalUnits::m } )
  {
    // The crystal is `setback` behind the face the distance is measured to
    const double omega = disk_fraction( diameter, distance + setback );
    const double cps = emission_rate * omega * intrinsic;
    const PeakDef peak = make_peak( 661.66, 1.5, cps*live_time );

    FluxToolCalc::Options options;
    options.distance = distance;
    options.correct_for_air = false;

    const shared_ptr<SpecUtils::Measurement> meas = make_measurement( live_time );
    const FluxToolCalc::PeakFlux flux = FluxToolCalc::peak_flux( peak, *drf, *meas, options );

    const double tol = 1.0E-8;  //percent
    BOOST_CHECK_CLOSE( flux.energy, 661.66, tol );
    BOOST_CHECK_CLOSE( flux.cps, cps, tol );
    BOOST_CHECK_CLOSE( flux.cps_stat, std::sqrt(cps*live_time)/live_time, tol );
    BOOST_CHECK_CLOSE( flux.geom_eff, omega, tol );
    BOOST_CHECK_CLOSE( flux.intrinsic_eff, intrinsic, tol );
    BOOST_CHECK_CLOSE( flux.gammas_4pi, emission_rate, tol );
    BOOST_CHECK_CLOSE( flux.flux_on_det, emission_rate*omega, tol );

    // Flux per cm2 at the detector face - the distance entered, not on to the crystal - and per *cm2*
    //  (the tool once gave it per mm2, 100x too small).
    const double distance_cm = distance / PhysicalUnits::cm;
    const double expected_flux_per_cm2 = emission_rate / (4.0*M_PI*distance_cm*distance_cm);
    BOOST_CHECK_CLOSE( flux.flux_per_cm2, expected_flux_per_cm2, tol );

    BOOST_CHECK_EQUAL( flux.air_transmission, 1.0 );
    BOOST_CHECK( !flux.background_subtracted );
    BOOST_CHECK_EQUAL( flux.background_cps, 0.0 );
    BOOST_CHECK( !flux.outside_drf_energy_range );
    BOOST_CHECK_EQUAL( flux.eff_uncert_frac, 0.0 );  //this DRF states no uncertainty

    // The CPS comes from the live time of the histogram passed in (i.e., the displayed spectrum).
    const shared_ptr<SpecUtils::Measurement> longer = make_measurement( 2.0f*live_time );
    const FluxToolCalc::PeakFlux half = FluxToolCalc::peak_flux( peak, *drf, *longer, options );
    BOOST_CHECK_CLOSE( half.cps, 0.5*cps, tol );
    BOOST_CHECK_CLOSE( half.gammas_4pi, 0.5*emission_rate, tol );
  }//for( loop over distances )

  // An invalid live time or distance is an error, not a table of nonsense
  FluxToolCalc::Options options;
  options.distance = 1.0*PhysicalUnits::m;
  const PeakDef peak = make_peak( 661.66, 1.5, 1000.0 );
  BOOST_CHECK_THROW( FluxToolCalc::peak_flux( peak, *drf, *make_measurement( 0.0f ), options ), std::exception );
  options.distance = -1.0;
  BOOST_CHECK_THROW( FluxToolCalc::peak_flux( peak, *drf, *make_measurement( 10.0f ), options ), std::exception );
}//far_field_round_trip


// With no options on, the refactored math must give exactly - to the last bit - what the tool gave
//  before it was moved out of the widget.
BOOST_AUTO_TEST_CASE( matches_previous_math_exactly )
{
  vector<shared_ptr<DetectorPeakResponse>> drfs;
  drfs.push_back( make_far_field_drf( 0.25, 5.0*PhysicalUnits::cm, 0.5*PhysicalUnits::cm ) );
  drfs.push_back( make_far_field_drf( 0.6, 7.62*PhysicalUnits::cm, 0.0 ) );
  drfs.push_back( make_drf_with_uncert() );

  // An absolute efficiency (at 25 cm), which the DRF converts to intrinsic itself
  {
    vector<DetectorPeakResponse::EnergyEffPoint> points;
    for( const pair<float,float> &ee : vector<pair<float,float>>{ {60.0f,0.004f}, {122.0f,0.006f}, {344.0f,0.003f},
                                                                 {662.0f,0.0016f}, {1332.0f,0.0009f}, {2614.0f,0.0005f} } )
    {
      DetectorPeakResponse::EnergyEffPoint p;
      p.energy = ee.first;
      p.efficiency = ee.second;
      points.push_back( p );
    }
    auto abs_drf = make_shared<DetectorPeakResponse>( "Absolute", "absolute efficiency at 25 cm" );
    abs_drf->setEfficiencyPoints( points, 5.08f*PhysicalUnits::cm, 25.0*PhysicalUnits::cm, GeomType::FarFieldAbsolute );
    abs_drf->setDetectorSetback( 0.3*PhysicalUnits::cm );
    BOOST_REQUIRE( abs_drf->isValid() );
    drfs.push_back( abs_drf );
  }

  const vector<double> energies{ 59.54, 122.06, 661.66, 1332.49, 2614.5 };
  const vector<double> distances{ 1.0*PhysicalUnits::cm, 25.0*PhysicalUnits::cm, 100.0*PhysicalUnits::cm,
                                  3.0*PhysicalUnits::m, 50.0*PhysicalUnits::m };
  const float live_time = 123.4f;
  const shared_ptr<SpecUtils::Measurement> meas = make_measurement( live_time );

  for( const shared_ptr<DetectorPeakResponse> &drf : drfs )
  {
    for( const double energy : energies )
    {
      for( const double distance : distances )
      {
        const PeakDef peak = make_peak( energy, 1.2, 12345.6 );

        FluxToolCalc::Options options;
        options.distance = distance;
        options.correct_for_air = false;

        const FluxToolCalc::PeakFlux flux = FluxToolCalc::peak_flux( peak, *drf, *meas, options );
        const LegacyFlux expected = legacy_flux( peak, drf, live_time, distance );

        BOOST_CHECK_EQUAL( flux.cps, expected.cps );
        BOOST_CHECK_EQUAL( flux.cps_stat, expected.cpsUncert );
        BOOST_CHECK_EQUAL( flux.geom_eff, expected.geomEff );
        BOOST_CHECK_EQUAL( flux.intrinsic_eff, expected.intrinsic );
        BOOST_CHECK_EQUAL( flux.intrinsic_syst, expected.intrinsicSyst );
        BOOST_CHECK_EQUAL( flux.flux_on_det, expected.fluxOnDet );
        BOOST_CHECK_EQUAL( flux.flux_on_det_stat, expected.fluxOnDetStat );
        BOOST_CHECK_EQUAL( flux.flux_on_det_syst, expected.fluxOnDetSyst );
        BOOST_CHECK_EQUAL( flux.flux_per_cm2, expected.flux );
        BOOST_CHECK_EQUAL( flux.flux_per_cm2_stat, expected.fluxStat );
        BOOST_CHECK_EQUAL( flux.flux_per_cm2_syst, expected.fluxSyst );
        BOOST_CHECK_EQUAL( flux.gammas_4pi, expected.gammaInto4pi );
        BOOST_CHECK_EQUAL( flux.gammas_4pi_stat, expected.gammaInto4piStat );
        BOOST_CHECK_EQUAL( flux.gammas_4pi_syst, expected.gammaInto4piSyst );
      }//for( const double distance : distances )
    }//for( const double energy : energies )
  }//for( const shared_ptr<DetectorPeakResponse> &drf : drfs )

  // Make sure the DRF with a stated uncertainty actually exercised the systematic parts
  {
    const FluxToolCalc::Options options{ 1.0*PhysicalUnits::m, false, nullptr };
    const FluxToolCalc::PeakFlux flux = FluxToolCalc::peak_flux( make_peak( 661.66, 1.2, 1.0E4 ),
                                                                 *drfs[2], *meas, options );
    BOOST_CHECK( flux.eff_uncert_frac > 0.0 );
    BOOST_CHECK( flux.gammas_4pi_syst > 0.0 );
  }
}//matches_previous_math_exactly


// A fixed-geometry efficiency gives the sources emission rate directly - per cm2, here - with no
//  distance, flux at the detector, or air correction.
BOOST_AUTO_TEST_CASE( fixed_geometry )
{
  auto drf = make_shared<DetectorPeakResponse>( "PerCm2", "fixed geometry, per cm2" );
  drf->setIntrinsicEfficiencyFormula( "1.0E-4", 2.54f*PhysicalUnits::cm, PhysicalUnits::keV, 0.0f, 0.0f,
                                      GeomType::FixedGeomActPerCm2 );
  BOOST_REQUIRE( drf->isValid() && drf->isFixedGeometry() );
  BOOST_CHECK_EQUAL( DetectorPeakResponse::det_eff_geom_type_postfix( drf->geometryType() ), "/cm2" );

  const float live_time = 100.0f;
  const double cps = 10.0;
  const PeakDef peak = make_peak( 661.66, 1.5, cps*live_time );
  const shared_ptr<SpecUtils::Measurement> meas = make_measurement( live_time );

  FluxToolCalc::Options options;
  options.correct_for_air = true;  //must be ignored

  options.distance = 1.0*PhysicalUnits::m;
  const FluxToolCalc::PeakFlux near_flux = FluxToolCalc::peak_flux( peak, *drf, *meas, options );

  options.distance = 3.0*PhysicalUnits::m;
  const FluxToolCalc::PeakFlux far_flux = FluxToolCalc::peak_flux( peak, *drf, *meas, options );

  options.distance = -1.0;  //not even looked at
  const FluxToolCalc::PeakFlux no_dist_flux = FluxToolCalc::peak_flux( peak, *drf, *meas, options );

  // 10 cps / (1E-4 counts per gamma/s/cm2) = 1E5 gammas/s/cm2, wherever the detector is
  BOOST_CHECK_CLOSE( near_flux.gammas_4pi, 1.0E5, 1.0E-4 );
  BOOST_CHECK_EQUAL( near_flux.gammas_4pi, far_flux.gammas_4pi );
  BOOST_CHECK_EQUAL( near_flux.gammas_4pi, no_dist_flux.gammas_4pi );
  BOOST_CHECK_EQUAL( near_flux.geom_eff, 1.0 );
  BOOST_CHECK_EQUAL( near_flux.air_transmission, 1.0 );
  BOOST_CHECK_EQUAL( near_flux.flux_on_det, 0.0 );
  BOOST_CHECK_EQUAL( near_flux.flux_per_cm2, 0.0 );
  BOOST_CHECK_EQUAL( drf->geometricEfficiency( 1.0*PhysicalUnits::m ), 1.0 );
  BOOST_CHECK_EQUAL( drf->modeledCrystalSetback(), 0.0 );
}//fixed_geometry


// The air correction only changes gammas-into-4pi; the fluxes are what reach the detector.
BOOST_AUTO_TEST_CASE( air_attenuation )
{
  const shared_ptr<DetectorPeakResponse> drf = make_far_field_drf( 0.3, 5.0*PhysicalUnits::cm, 0.5*PhysicalUnits::cm );
  const float live_time = 600.0f;
  const shared_ptr<SpecUtils::Measurement> meas = make_measurement( live_time );
  const double distance = 10.0*PhysicalUnits::m;

  FluxToolCalc::Options without_air;
  without_air.distance = distance;
  without_air.correct_for_air = false;

  FluxToolCalc::Options with_air = without_air;
  with_air.correct_for_air = true;

  for( const double energy : { 59.54, 661.66 } )
  {
    const PeakDef peak = make_peak( energy, 1.2, 5.0E4 );
    const FluxToolCalc::PeakFlux plain = FluxToolCalc::peak_flux( peak, *drf, *meas, without_air );
    const FluxToolCalc::PeakFlux corrected = FluxToolCalc::peak_flux( peak, *drf, *meas, with_air );

    // The air path is to the detector face, not on to the crystal
    const double transmission = std::exp( -GammaInteractionCalc::transmission_coefficient_air( energy, distance ) );
    BOOST_TEST_MESSAGE( "Air attenuation at " << energy << " keV: "
                        << 100.0*(1.0 - std::exp( -GammaInteractionCalc::transmission_coefficient_air( energy, 1.0*PhysicalUnits::m ) ))
                        << "% at 1 m, " << 100.0*(1.0 - transmission) << "% at 10 m" );
    BOOST_CHECK_CLOSE( corrected.air_transmission, transmission, 1.0E-10 );
    BOOST_CHECK_CLOSE( corrected.gammas_4pi, plain.gammas_4pi / transmission, 1.0E-10 );
    BOOST_CHECK_CLOSE( corrected.gammas_4pi_stat, plain.gammas_4pi_stat / transmission, 1.0E-10 );

    // Roughly 9% lost through 10 m of air at 662 keV, and 18% at 60 keV
    if( energy > 600.0 )
      BOOST_CHECK_MESSAGE( (1.0 - transmission > 0.05) && (1.0 - transmission < 0.15),
                           "Unexpected air attenuation at 662 keV, 10 m: " << (1.0 - transmission) );
    else
      BOOST_CHECK_MESSAGE( (1.0 - transmission > 0.10) && (1.0 - transmission < 0.30),
                           "Unexpected air attenuation at 60 keV, 10 m: " << (1.0 - transmission) );

    BOOST_CHECK_EQUAL( corrected.flux_per_cm2, plain.flux_per_cm2 );
    BOOST_CHECK_EQUAL( corrected.flux_on_det, plain.flux_on_det );
    BOOST_CHECK_EQUAL( corrected.cps, plain.cps );
    BOOST_CHECK_EQUAL( plain.air_transmission, 1.0 );
  }//for( const double energy : { 59.54, 661.66 } )
}//air_attenuation


// Background-peak subtraction, done the way the Activity/Shielding fit does it.
BOOST_AUTO_TEST_CASE( background_subtraction )
{
  // The user scale factor wins when set; otherwise the live-time ratio.
  BOOST_CHECK_EQUAL( GammaInteractionCalc::background_peak_scale( 0.25, 300.0, 600.0 ), 0.25 );
  BOOST_CHECK_EQUAL( GammaInteractionCalc::background_peak_scale( 0.0, 300.0, 600.0 ), 0.5 );
  BOOST_CHECK_EQUAL( GammaInteractionCalc::background_peak_scale( -1.0, 300.0, 600.0 ), 0.5 );
  BOOST_CHECK_EQUAL( GammaInteractionCalc::background_peak_scale( std::numeric_limits<double>::quiet_NaN(), 300.0, 600.0 ), 0.5 );
  BOOST_CHECK_EQUAL( GammaInteractionCalc::background_peak_scale( std::numeric_limits<double>::infinity(), 300.0, 600.0 ), 0.5 );

  // Only Gaussian background peaks are scaled and kept
  const shared_ptr<SpecUtils::Measurement> back_hist = make_measurement( 600.0f );
  const PeakDef data_defined( 1170.0, 1176.0, 1173.2, back_hist, nullptr );
  BOOST_REQUIRE( !data_defined.gausPeak() );

  const vector<PeakDef> raw_back{ make_peak( 609.4, 1.0, 800.0 ),     //within 1 sigma of 609.3
                                  make_peak( 1461.3, 1.3, 1200.0 ),   //within 1 sigma of 1460.8
                                  make_peak( 1460.2, 1.3, 400.0 ),    //also within 1 sigma of 1460.8
                                  make_peak( 1462.5, 1.3, 9999.0 ),   //1.7 keV away - outside 1 sigma
                                  data_defined };
  vector<double> skipped;
  const vector<PeakDef> scaled = GammaInteractionCalc::scale_background_peaks( raw_back, 0.5, &skipped );
  BOOST_REQUIRE_EQUAL( scaled.size(), 4u );
  BOOST_REQUIRE_EQUAL( skipped.size(), 1u );
  BOOST_CHECK_CLOSE( skipped[0], 1173.2, 1.0E-6 );
  BOOST_CHECK_CLOSE( scaled[0].peakArea(), 400.0, 1.0E-10 );
  BOOST_CHECK_CLOSE( scaled[0].peakAreaUncert(), 0.5*std::sqrt(800.0), 1.0E-10 );

  // Matching: within 1 sigma of the foreground peak is subtracted; several matches are summed, and
  //  uncertainties added in quadrature.
  const PeakDef fore_1460 = make_peak( 1460.8, 1.2, 9000.0 );
  const GammaInteractionCalc::BackgroundOverlap overlap = GammaInteractionCalc::overlapping_background( fore_1460, scaled );
  BOOST_CHECK_CLOSE( overlap.counts, 0.5*1200.0 + 0.5*400.0, 1.0E-10 );
  BOOST_CHECK_CLOSE( overlap.uncert2, 0.25*1200.0 + 0.25*400.0, 1.0E-10 );

  const PeakDef fore_1000 = make_peak( 1000.0, 1.2, 9000.0 );
  BOOST_CHECK_EQUAL( GammaInteractionCalc::overlapping_background( fore_1000, scaled ).counts, 0.0 );

  // Through the Flux Tool math
  const shared_ptr<DetectorPeakResponse> drf = make_far_field_drf( 0.2, 7.62*PhysicalUnits::cm, 0.0 );
  const float live_time = 300.0f;
  const shared_ptr<SpecUtils::Measurement> meas = make_measurement( live_time );

  FluxToolCalc::Options options;
  options.distance = 1.0*PhysicalUnits::m;
  options.correct_for_air = false;

  const FluxToolCalc::PeakFlux gross = FluxToolCalc::peak_flux( fore_1460, *drf, *meas, options );
  options.scaled_background_peaks = &scaled;
  const FluxToolCalc::PeakFlux net = FluxToolCalc::peak_flux( fore_1460, *drf, *meas, options );

  BOOST_CHECK( !gross.background_subtracted );
  BOOST_CHECK( net.background_subtracted );
  BOOST_CHECK_CLOSE( net.background_cps, 800.0 / live_time, 1.0E-10 );
  BOOST_CHECK_CLOSE( net.background_cps_stat, std::sqrt(400.0) / live_time, 1.0E-10 );
  BOOST_CHECK_CLOSE( net.cps, (9000.0 - 800.0) / live_time, 1.0E-10 );
  BOOST_CHECK_CLOSE( net.cps_stat, std::sqrt(9000.0 + 400.0) / live_time, 1.0E-10 );
  BOOST_CHECK_CLOSE( net.gammas_4pi, gross.gammas_4pi * (9000.0 - 800.0) / 9000.0, 1.0E-8 );
  BOOST_CHECK_CLOSE( net.flux_per_cm2, gross.flux_per_cm2 * (9000.0 - 800.0) / 9000.0, 1.0E-8 );

  // A peak with no background under it is untouched
  const FluxToolCalc::PeakFlux untouched = FluxToolCalc::peak_flux( fore_1000, *drf, *meas, options );
  BOOST_CHECK( !untouched.background_subtracted );
  BOOST_CHECK_EQUAL( untouched.cps, 9000.0 / live_time );

  // More background than foreground (e.g., K-40 in a shielded measurement) gives negative net values,
  //  whose uncertainties must still be positive - in particular the systematic (efficiency) parts.
  {
    const shared_ptr<DetectorPeakResponse> drf_uncert = make_drf_with_uncert();
    const PeakDef small_fore = make_peak( 1460.8, 1.2, 500.0 );
    const FluxToolCalc::PeakFlux negative = FluxToolCalc::peak_flux( small_fore, *drf_uncert, *meas, options );
    BOOST_CHECK( negative.background_subtracted );
    BOOST_CHECK_CLOSE( negative.cps, (500.0 - 800.0) / live_time, 1.0E-10 );
    BOOST_REQUIRE( negative.eff_uncert_frac > 0.0 );
    BOOST_CHECK( negative.gammas_4pi < 0.0 );
    BOOST_CHECK( negative.flux_per_cm2 < 0.0 );
    BOOST_CHECK( negative.flux_on_det < 0.0 );
    BOOST_CHECK( negative.cps_stat > 0.0 );
    BOOST_CHECK( negative.gammas_4pi_stat > 0.0 );
    BOOST_CHECK_CLOSE( negative.gammas_4pi_syst, -negative.gammas_4pi * negative.eff_uncert_frac, 1.0E-10 );
    BOOST_CHECK_CLOSE( negative.flux_per_cm2_syst, -negative.flux_per_cm2 * negative.eff_uncert_frac, 1.0E-10 );
    BOOST_CHECK_CLOSE( negative.flux_on_det_syst, -negative.flux_on_det * negative.eff_uncert_frac, 1.0E-10 );
  }
}//background_subtraction


// Which unfit background peaks the "missing background peaks" prompt offers (shared with the
//  Activity/Shielding fit).
BOOST_AUTO_TEST_CASE( missing_background_peak_candidates )
{
  const auto peak_ptr = []( const double mean, const double sigma ){
    return make_shared<const PeakDef>( make_peak( mean, sigma, 1000.0 ) );
  };

  const vector<shared_ptr<const PeakDef>> fore{ peak_ptr( 609.3, 1.0 ), peak_ptr( 1460.8, 1.2 ),
                                                peak_ptr( 1461.2, 1.2 ), peak_ptr( 2614.5, 1.5 ) };

  // 609 is covered (within 1 sigma); the fitted 1465.0 is neither within 1 sigma of 1460.8 or 1461.2,
  //  nor within a FWHM of the 1461.0 candidate, so doesnt stop it being offered.
  const deque<shared_ptr<const PeakDef>> user_back{ peak_ptr( 609.5, 1.0 ), peak_ptr( 1465.0, 1.2 ) };

  // 1461.0 is the closest auto-search peak to both 1460.8 and 1461.2, so should be offered once;
  //  nothing is within a FWHM of 2614.5.
  const deque<shared_ptr<const PeakDef>> auto_back{ peak_ptr( 609.4, 1.0 ), peak_ptr( 1461.0, 1.2 ),
                                                    peak_ptr( 1463.0, 1.2 ), peak_ptr( 2620.0, 1.5 ) };

  const vector<shared_ptr<const PeakDef>> candidates
                     = PeakSearchGuiUtils::missing_background_peak_candidates( fore, user_back, auto_back );
  BOOST_REQUIRE_EQUAL( candidates.size(), 1u );
  BOOST_CHECK_CLOSE( candidates[0]->mean(), 1461.0, 1.0E-8 );

  // A background peak the user already fit, though too far off (1.5 sigma) to count as under the
  //  foreground peak, is not offered again - the fit-and-add step would skip it, and the prompt would
  //  just come back every time.
  {
    const vector<shared_ptr<const PeakDef>> fore_662{ peak_ptr( 661.7, 0.6 ) };
    const deque<shared_ptr<const PeakDef>> user_662{ peak_ptr( 662.6, 0.6 ) };
    const deque<shared_ptr<const PeakDef>> auto_662{ peak_ptr( 662.4, 0.6 ) };
    BOOST_CHECK( PeakSearchGuiUtils::missing_background_peak_candidates( fore_662, user_662, auto_662 ).empty() );
    BOOST_CHECK_EQUAL( PeakSearchGuiUtils::missing_background_peak_candidates( fore_662, {}, auto_662 ).size(), 1u );
  }

  // A non-Gaussian background peak the user fit doesnt count as covering the foreground peak (it
  //  wouldnt be subtracted), but an automated-search peak on top of it isnt offered either - fitting
  //  and adding it would skip it as already there.  And a non-Gaussian automated-search peak is never
  //  offered.
  {
    const shared_ptr<SpecUtils::Measurement> hist = make_measurement( 600.0f );
    const auto data_defined = make_shared<const PeakDef>( PeakDef( 1456.0, 1466.0, 1460.9, hist, nullptr ) );
    BOOST_REQUIRE( !data_defined->gausPeak() );
    const vector<shared_ptr<const PeakDef>> fore_1460{ peak_ptr( 1460.8, 1.2 ) };
    const deque<shared_ptr<const PeakDef>> user_dd{ data_defined };
    const deque<shared_ptr<const PeakDef>> auto_1461{ peak_ptr( 1461.0, 1.2 ) };
    BOOST_CHECK( PeakSearchGuiUtils::missing_background_peak_candidates( fore_1460, user_dd, auto_1461 ).empty() );
    const deque<shared_ptr<const PeakDef>> auto_dd{ data_defined };
    BOOST_CHECK( PeakSearchGuiUtils::missing_background_peak_candidates( fore_1460, {}, auto_dd ).empty() );
  }

  // With no user background peaks, 609 is offered too (and results are sorted by energy)
  const vector<shared_ptr<const PeakDef>> all
                     = PeakSearchGuiUtils::missing_background_peak_candidates( fore, {}, auto_back );
  BOOST_REQUIRE_EQUAL( all.size(), 2u );
  BOOST_CHECK_CLOSE( all[0]->mean(), 609.4, 1.0E-8 );
  BOOST_CHECK_CLOSE( all[1]->mean(), 1461.0, 1.0E-8 );
}//missing_background_peak_candidates


// Peaks outside the energy range the DRF states it is valid for are flagged; DRFs that dont state one
//  never flag anything.
BOOST_AUTO_TEST_CASE( energy_range )
{
  const shared_ptr<DetectorPeakResponse> ranged = make_far_field_drf( 0.3, 5.0*PhysicalUnits::cm, 0.0, 50.0f, 3000.0f );
  BOOST_CHECK( FluxToolCalc::outside_energy_range( *ranged, 30.0 ) );
  BOOST_CHECK( !FluxToolCalc::outside_energy_range( *ranged, 661.66 ) );
  BOOST_CHECK( FluxToolCalc::outside_energy_range( *ranged, 3500.0 ) );

  const shared_ptr<DetectorPeakResponse> unranged = make_far_field_drf( 0.3, 5.0*PhysicalUnits::cm, 0.0 );
  BOOST_CHECK( !FluxToolCalc::outside_energy_range( *unranged, 30.0 ) );
  BOOST_CHECK( !FluxToolCalc::outside_energy_range( *unranged, 3500.0 ) );

  // Some GADRAS detectors give only a lower energy, which isnt a range
  shared_ptr<DetectorPeakResponse> half = make_far_field_drf( 0.3, 5.0*PhysicalUnits::cm, 0.0 );
  half->setEnergyRange( 25.0, 0.0 );
  BOOST_CHECK( !FluxToolCalc::outside_energy_range( *half, 10.0 ) );
  BOOST_CHECK( !FluxToolCalc::outside_energy_range( *half, 661.66 ) );

  FluxToolCalc::Options options;
  options.distance = 1.0*PhysicalUnits::m;
  const shared_ptr<SpecUtils::Measurement> meas = make_measurement( 100.0f );
  BOOST_CHECK( FluxToolCalc::peak_flux( make_peak( 30.85, 0.5, 1000.0 ), *ranged, *meas, options ).outside_drf_energy_range );
  BOOST_CHECK( !FluxToolCalc::peak_flux( make_peak( 356.0, 0.8, 1000.0 ), *ranged, *meas, options ).outside_drf_energy_range );
}//energy_range


// A DRF with a geometry-modeled (CeeLo) response: the efficiency, and so gammas-into-4pi, come from
//  the full response; the geometric efficiency is the responses crystal, recessed behind the face by
//  its endcap - whatever the DRFs own stated diameter and setback say.
BOOST_AUTO_TEST_CASE( geometry_modeled_response )
{
  BOOST_REQUIRE_MESSAGE( !g_test_data_dir.empty(), "--testfiledir not given" );
  const string path = SpecUtils::append_path( SpecUtils::append_path( g_test_data_dir, "ceelo_drf" ),
                                               "nai3x3_response.xml" );

  auto drf = make_shared<DetectorPeakResponse>( "nai3x3", "golden fixture" );
  drf->setIntrinsicEfficiencyFormula( "1.0", 5.0*PhysicalUnits::cm, PhysicalUnits::keV, 0.0f, 0.0f,
                                      GeomType::FarFieldIntrinsic );
  drf->setDetectorSetback( 3.0*PhysicalUnits::cm );  //deliberately not the responses
  drf->setCeeloResponse( ceelo::DetectorResponse::from_xml_string( read_file( path ) ) );
  BOOST_REQUIRE( drf->ceeloResponse() );

  const ceelo::GeometryDescriptor &gd = drf->ceeloResponse()->descriptor;
  const double crystal_diam = 2.0 * CeeLoUtils::crystalHalfExtent( gd ) * PhysicalUnits::cm;
  const double crystal_setback = gd.endcap_front_offset_cm() * PhysicalUnits::cm;
  BOOST_CHECK_CLOSE( drf->modeledCrystalDiameter(), crystal_diam, 1.0E-10 );
  BOOST_CHECK_CLOSE( drf->modeledCrystalSetback(), crystal_setback, 1.0E-10 );

  const double distance = 25.0*PhysicalUnits::cm;
  const float live_time = 500.0f;
  const shared_ptr<SpecUtils::Measurement> meas = make_measurement( live_time );
  const double emission_rate = 2.0E5;

  for( const double energy : { 121.78, 661.66, 1332.49 } )
  {
    const double eff = drf->efficiency( energy, distance );
    BOOST_REQUIRE( eff > 0.0 );
    const PeakDef peak = make_peak( energy, 3.0, emission_rate*eff*live_time );

    FluxToolCalc::Options options;
    options.distance = distance;
    options.correct_for_air = false;
    const FluxToolCalc::PeakFlux flux = FluxToolCalc::peak_flux( peak, *drf, *meas, options );

    BOOST_CHECK_CLOSE( flux.gammas_4pi, emission_rate, 1.0E-6 );
    BOOST_CHECK_CLOSE( flux.geom_eff, disk_fraction( crystal_diam, distance + crystal_setback ), 1.0E-8 );
    BOOST_CHECK_CLOSE( flux.geom_eff * flux.intrinsic_eff, eff, 1.0E-8 );
    BOOST_CHECK_CLOSE( flux.flux_on_det, emission_rate * flux.geom_eff, 1.0E-6 );

    // Responses are computed in vacuum, so the air correction applies to them as to any other
    options.correct_for_air = true;
    const FluxToolCalc::PeakFlux with_air = FluxToolCalc::peak_flux( peak, *drf, *meas, options );
    BOOST_CHECK( with_air.air_transmission < 1.0 );
    BOOST_CHECK_CLOSE( with_air.gammas_4pi, flux.gammas_4pi / with_air.air_transmission, 1.0E-10 );
  }//for( loop over energies )
}//geometry_modeled_response
