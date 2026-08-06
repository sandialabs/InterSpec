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
#include <cstring>
#include <string>
#include <memory>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

#ifdef _WIN32
#include "winsock2.h"
#include "Windows.h"
#endif

#define BOOST_TEST_MODULE ExportCAM_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/CAMIO.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ExportSpecFileCAM.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using namespace boost::unit_test;

namespace
{
  string g_test_file_dir;

  /** A nominal HPGe `FWHM = A0 + A1*sqrt(E)` shape calibration, used as the line-clustering
   resolution in the `build_genie_library(...)` tests.
   */
  const pair<float,float> s_test_hpge_fwhm{ 1.0f, 0.035f };

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
        if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) )
        {
          datadir = d;
          break;
        }
      }
    }

    if( g_test_file_dir.empty() )
    {
      const vector<string> candidate_test_dirs{
        "test_data", "../test_data", "target/testing/test_data",
        "../target/testing/test_data", "../../target/testing/test_data",
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
    }

    if( datadir.empty() )
      throw runtime_error( "Could not find data directory" );
    if( g_test_file_dir.empty() )
      throw runtime_error( "Could not find test file directory" );

    DecayDataBaseServer::setDecayXmlFile( SpecUtils::append_path(datadir, "sandia.decay.xml") );
    InterSpec::setStaticDataDirectory( datadir );
  }//void set_data_dir()


  /** Finds a gamma transition anywhere in `nuc`'s decay chain (not just its direct decays -
   e.g., Cs137's well-known 661.7 keV gamma actually comes from its daughter Ba137m) with
   energy near `energy_kev`.
   */
  pair<const SandiaDecay::Transition *, int> find_gamma_transition( const SandiaDecay::Nuclide *nuc,
                                                                    const float energy_kev,
                                                                    const float tol_kev = 1.0f )
  {
    vector<const SandiaDecay::Nuclide *> to_visit{ nuc };
    set<const SandiaDecay::Nuclide *> visited;

    while( !to_visit.empty() )
    {
      const SandiaDecay::Nuclide * const cur = to_visit.back();
      to_visit.pop_back();
      if( !cur || !visited.insert(cur).second )
        continue;

      for( const SandiaDecay::Transition *trans : cur->decaysToChildren )
      {
        if( !trans )
          continue;

        for( size_t i = 0; i < trans->products.size(); ++i )
        {
          const SandiaDecay::RadParticle &p = trans->products[i];
          if( (p.type == SandiaDecay::ProductType::GammaParticle)
             && (std::fabs(p.energy - energy_kev) < tol_kev) )
          {
            return { trans, static_cast<int>(i) };
          }
        }

        if( trans->child )
          to_visit.push_back( trans->child );
      }//for( loop over cur's decay transitions )
    }//while( !to_visit.empty() )

    return { nullptr, -1 };
  }//find_gamma_transition(...)


  /** Builds a simple Gaussian peak at `energy_kev`, optionally assigned to a nuclide gamma. */
  shared_ptr<const PeakDef> make_peak( const float energy_kev, const SandiaDecay::Nuclide *nuc = nullptr )
  {
    auto cont = make_shared<PeakContinuum>();
    cont->setType( PeakContinuum::OffsetType::Linear );
    cont->setRange( energy_kev - 15.0, energy_kev + 15.0 );

    auto peak = make_shared<PeakDef>( energy_kev, 1.5, 1000.0 );
    peak->setContinuum( cont );

    if( nuc )
    {
      const pair<const SandiaDecay::Transition *, int> trans = find_gamma_transition( nuc, energy_kev, 2.0f );
      BOOST_REQUIRE_MESSAGE( trans.first, "Could not find a gamma transition for "
                             << nuc->symbol << " near " << energy_kev << " keV" );
      peak->setNuclearTransition( nuc, trans.first, trans.second, PeakDef::SourceGammaType::NormalGamma );
    }

    return peak;
  }//make_peak(...)
}//namespace


BOOST_AUTO_TEST_CASE( readExampleMultiNuclideLibrary )
{
  set_data_dir();

  const string path = SpecUtils::append_path( g_test_file_dir, "CamLibrary/ExampleMultiNuclide.nlb" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file(path), "Missing test file: " << path );

  ifstream infile( path.c_str(), ios::in | ios::binary );
  BOOST_REQUIRE( infile.is_open() );
  vector<unsigned char> filedata( (istreambuf_iterator<char>(infile)), istreambuf_iterator<char>() );
  BOOST_REQUIRE( !filedata.empty() );

  CAMInputOutput::CAMIO reader;
  BOOST_REQUIRE_NO_THROW( reader.ReadFile(filedata) );

  const vector<CAMInputOutput::Line> &lines = reader.GetLines();
  const vector<CAMInputOutput::Nuclide> &nuclides = reader.GetNuclides();

  BOOST_REQUIRE_MESSAGE( !nuclides.empty(), "Should have read at least one nuclide" );
  BOOST_CHECK_MESSAGE( !lines.empty(), "Should have read at least one line" );

  const set<string> expected_substrs = {
    "47", "52", "57", "76", "89", "111", "127", "196", "198", "237", "239"
  };
  set<string> found_substrs;
  for( const CAMInputOutput::Nuclide &nuc : nuclides )
  {
    for( const string &s : expected_substrs )
    {
      if( nuc.Name.find(s) != string::npos )
        found_substrs.insert( s );
    }
  }

  BOOST_CHECK_MESSAGE( nuclides.size() >= 11, "Expected at least 11 nuclides, got " << nuclides.size() );
  BOOST_CHECK_MESSAGE( found_substrs.size() == expected_substrs.size(),
                       "Only matched " << found_substrs.size() << " of " << expected_substrs.size()
                       << " expected nuclide mass-number substrings" );
}//BOOST_AUTO_TEST_CASE( readExampleMultiNuclideLibrary )


BOOST_AUTO_TEST_CASE( buildLibraryPeakLinesOnly )
{
  set_data_dir();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );

  const SandiaDecay::Nuclide * const co60 = db->nuclide( "Co60" );
  BOOST_REQUIRE( co60 );

  deque<shared_ptr<const PeakDef>> peaks;
  peaks.push_back( make_peak( 1173.2f, co60 ) );
  peaks.push_back( make_peak( 1332.5f, co60 ) );

  vector<string> warnings;
  const vector<ExportSpecFileCAM::GenieLibrarySource> sources = ExportSpecFileCAM::build_genie_library(
                                peaks, ExportSpecFileCAM::GenieLibraryLineMode::PeakLinesOnly,
                                1.0, true, s_test_hpge_fwhm, {}, &warnings );

  BOOST_REQUIRE_EQUAL( sources.size(), 1u );
  const ExportSpecFileCAM::GenieLibrarySource &src = sources[0];
  BOOST_CHECK_EQUAL( src.name, "Co60" );
  BOOST_CHECK( src.nuclide == co60 );
  BOOST_CHECK_CLOSE( src.half_life_seconds, static_cast<float>(co60->halfLife / PhysicalUnits::second), 0.1 );
  BOOST_REQUIRE_EQUAL( src.lines.size(), 2u );

  bool have_key_line = false;
  for( const ExportSpecFileCAM::GenieLibraryLine &line : src.lines )
  {
    BOOST_CHECK( line.yield > 0.0f );
    BOOST_CHECK( line.yield <= 1.0f );
    BOOST_CHECK( !line.is_xray );
    have_key_line = have_key_line || line.is_key_line;
  }
  BOOST_CHECK_MESSAGE( have_key_line, "Exactly one line should have been marked the key line" );

  // Only the 1332 keV line should be marked key (higher score: higher energy, similar yield)
  for( const ExportSpecFileCAM::GenieLibraryLine &line : src.lines )
  {
    if( line.is_key_line )
      BOOST_CHECK_CLOSE( line.energy, 1332.5f, 0.1 );
  }
}//BOOST_AUTO_TEST_CASE( buildLibraryPeakLinesOnly )


BOOST_AUTO_TEST_CASE( buildLibraryAllLinesAboveThreshold )
{
  set_data_dir();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );

  const SandiaDecay::Nuclide * const co60 = db->nuclide( "Co60" );
  BOOST_REQUIRE( co60 );

  deque<shared_ptr<const PeakDef>> peaks;
  peaks.push_back( make_peak( 1332.5f, co60 ) );

  // Very low threshold - should end up with (at least) the two well known Co60 lines.
  vector<ExportSpecFileCAM::GenieLibrarySource> sources = ExportSpecFileCAM::build_genie_library(
                                peaks, ExportSpecFileCAM::GenieLibraryLineMode::AllLinesAboveThreshold,
                                1.0, false, s_test_hpge_fwhm, {} );
  BOOST_REQUIRE_EQUAL( sources.size(), 1u );
  BOOST_CHECK_MESSAGE( sources[0].lines.size() >= 2, "Expected at least the 1173/1332 keV lines" );

  // A very high threshold should exclude everything but the most intense line(s).
  sources = ExportSpecFileCAM::build_genie_library(
                                peaks, ExportSpecFileCAM::GenieLibraryLineMode::AllLinesAboveThreshold,
                                99.0, false, s_test_hpge_fwhm, {} );
  BOOST_REQUIRE_EQUAL( sources.size(), 1u );
  BOOST_CHECK_MESSAGE( sources[0].lines.size() <= 2, "A 99% threshold should leave very few lines" );
}//BOOST_AUTO_TEST_CASE( buildLibraryAllLinesAboveThreshold )


BOOST_AUTO_TEST_CASE( combineUnresolvableLines )
{
  set_data_dir();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );

  // Two lines 0.1 keV apart should merge into one when combining is enabled, using synthetic
  //  peaks assigned to the same nuclide (Co60) but hand-set energies close together via the
  //  pure GenieLibraryLine merge logic - exercised indirectly through build_genie_library isn't
  //  practical for exact energy control (real line energies aren't adjustable), so we go through
  //  the AllLinesAboveThreshold path on a nuclide with closely-spaced lines (Np-237 or similar
  //  isn't guaranteed available), and instead verify the general property: combining never
  //  produces *more* lines than not combining.
  const SandiaDecay::Nuclide * const co60 = db->nuclide( "Co60" );
  BOOST_REQUIRE( co60 );

  deque<shared_ptr<const PeakDef>> peaks;
  peaks.push_back( make_peak( 1332.5f, co60 ) );

  const vector<ExportSpecFileCAM::GenieLibrarySource> combined = ExportSpecFileCAM::build_genie_library(
                                peaks, ExportSpecFileCAM::GenieLibraryLineMode::AllLinesAboveThreshold,
                                0.01, true, s_test_hpge_fwhm, {} );
  const vector<ExportSpecFileCAM::GenieLibrarySource> uncombined = ExportSpecFileCAM::build_genie_library(
                                peaks, ExportSpecFileCAM::GenieLibraryLineMode::AllLinesAboveThreshold,
                                0.01, false, s_test_hpge_fwhm, {} );

  BOOST_REQUIRE_EQUAL( combined.size(), 1u );
  BOOST_REQUIRE_EQUAL( uncombined.size(), 1u );
  BOOST_CHECK_LE( combined[0].lines.size(), uncombined[0].lines.size() );
}//BOOST_AUTO_TEST_CASE( combineUnresolvableLines )


BOOST_AUTO_TEST_CASE( toLibraryLinesRespectsIncluded )
{
  set_data_dir();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );

  const SandiaDecay::Nuclide * const co60 = db->nuclide( "Co60" );
  BOOST_REQUIRE( co60 );

  deque<shared_ptr<const PeakDef>> peaks;
  peaks.push_back( make_peak( 1173.2f, co60 ) );
  peaks.push_back( make_peak( 1332.5f, co60 ) );

  vector<ExportSpecFileCAM::GenieLibrarySource> sources = ExportSpecFileCAM::build_genie_library(
                                peaks, ExportSpecFileCAM::GenieLibraryLineMode::PeakLinesOnly,
                                1.0, true, s_test_hpge_fwhm, {} );
  BOOST_REQUIRE_EQUAL( sources.size(), 1u );
  BOOST_REQUIRE_EQUAL( sources[0].lines.size(), 2u );

  // All included by default
  vector<CAMInputOutput::CnfGenieExtras::LibraryLine> ll = ExportSpecFileCAM::to_library_lines( sources );
  BOOST_CHECK_EQUAL( ll.size(), 2u );

  // Un-checking one line excludes just it
  sources[0].lines[0].included = false;
  ll = ExportSpecFileCAM::to_library_lines( sources );
  BOOST_CHECK_EQUAL( ll.size(), 1u );

  // Un-checking the whole source excludes everything, regardless of line-level checks
  sources[0].lines[0].included = true;
  sources[0].included = false;
  ll = ExportSpecFileCAM::to_library_lines( sources );
  BOOST_CHECK_EQUAL( ll.size(), 0u );
}//BOOST_AUTO_TEST_CASE( toLibraryLinesRespectsIncluded )


BOOST_AUTO_TEST_CASE( writeCnfWithLibraryRoundTrips )
{
  set_data_dir();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );

  const SandiaDecay::Nuclide * const co60 = db->nuclide( "Co60" );
  const SandiaDecay::Nuclide * const cs137 = db->nuclide( "Cs137" );
  BOOST_REQUIRE( co60 );
  BOOST_REQUIRE( cs137 );

  deque<shared_ptr<const PeakDef>> peaks;
  peaks.push_back( make_peak( 1173.2f, co60 ) );
  peaks.push_back( make_peak( 1332.5f, co60 ) );
  peaks.push_back( make_peak( 661.66f, cs137 ) );

  const vector<ExportSpecFileCAM::GenieLibrarySource> sources = ExportSpecFileCAM::build_genie_library(
                                peaks, ExportSpecFileCAM::GenieLibraryLineMode::PeakLinesOnly,
                                1.0, true, s_test_hpge_fwhm, {} );
  BOOST_REQUIRE_EQUAL( sources.size(), 2u );

  // Build a minimal synthetic spectrum to write, along with the library.
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( 512, {0.0f, 3.0f}, {} );

  auto m = make_shared<SpecUtils::Measurement>();
  auto counts = make_shared<vector<float>>( 512, 10.0f );
  m->set_gamma_counts( counts, 300.0f, 300.0f );
  m->set_energy_calibration( cal );
  m->set_sample_number( 1 );

  SpecMeas meas;
  meas.add_measurement( m, false );

  CAMInputOutput::CnfGenieExtras extras;
  extras.library_lines = ExportSpecFileCAM::to_library_lines( sources );
  extras.shape_cal = make_pair( 1.0f, 0.03f );

  ostringstream strm;
  const bool success = meas.write_cnf( strm, {}, {}, &extras );
  BOOST_REQUIRE( success );

  const string filedata_str = strm.str();
  BOOST_REQUIRE( !filedata_str.empty() );
  vector<unsigned char> filedata( filedata_str.begin(), filedata_str.end() );

  CAMInputOutput::CAMIO reader;
  BOOST_REQUIRE_NO_THROW( reader.ReadFile(filedata) );

  // GetNuclides() associates lines to nuclides via NuclideIndex, so GetLines() must be called
  //  first; also, neither accessor clears its internal cache, so each must only be called once.
  const vector<CAMInputOutput::Line> &read_lines = reader.GetLines();
  const vector<CAMInputOutput::Nuclide> &read_nucs = reader.GetNuclides();

  BOOST_REQUIRE_EQUAL( read_nucs.size(), 2u );
  BOOST_REQUIRE_EQUAL( read_lines.size(), 3u );

  // Names are written padded to a fixed width (per the CAM format), so trim before comparing.
  set<string> read_names;
  for( const CAMInputOutput::Nuclide &n : read_nucs )
    read_names.insert( SpecUtils::trim_copy(n.Name) );
  BOOST_CHECK( read_names.count("Co60") );
  BOOST_CHECK( read_names.count("Cs137") );

  vector<float> shape_cal = reader.GetShapeCalibration();
  while( !shape_cal.empty() && (shape_cal.back() == 0.0f) )
    shape_cal.pop_back();
  BOOST_REQUIRE_EQUAL( shape_cal.size(), 2u );
  BOOST_CHECK_CLOSE( shape_cal[0], 1.0f, 1.0 );
  BOOST_CHECK_CLOSE( shape_cal[1], 0.03f, 1.0 );
}//BOOST_AUTO_TEST_CASE( writeCnfWithLibraryRoundTrips )


BOOST_AUTO_TEST_CASE( fitGenieFwhmFromDrf )
{
  set_data_dir();

  auto drf = make_shared<DetectorPeakResponse>( "TestDRF", "Test" );
  const vector<float> coefs{ 1.0f, 0.03f };
  drf->setFwhmCoefficients( coefs, DetectorPeakResponse::ResolutionFnctForm::kConstantPlusSqrtEnergy );
  drf->setEnergyRange( 59.0, 2614.0 );

  const pair<float,float> fit = ExportSpecFileCAM::fit_genie_fwhm_from_drf( *drf );

  // Since the source is already exactly the Genie form, the fit should reproduce it closely.
  BOOST_CHECK_CLOSE( fit.first, coefs[0], 1.0 );
  BOOST_CHECK_CLOSE( fit.second, coefs[1], 1.0 );

  // Round-trip through CAMIO
  CAMInputOutput::CAMIO writer;
  writer.AddShapeCalibration( fit.first, fit.second );
  writer.AddEnergyCalibration( {0.0f, 3.0f} );
  writer.AddSpectrum( vector<uint32_t>( 512, 10 ) );
  vector<unsigned char> &filedata = writer.CreateFile();
  BOOST_REQUIRE( !filedata.empty() );

  CAMInputOutput::CAMIO reader;
  BOOST_REQUIRE_NO_THROW( reader.ReadFile(filedata) );
  vector<float> read_shape = reader.GetShapeCalibration();
  while( !read_shape.empty() && (read_shape.back() == 0.0f) )
    read_shape.pop_back();
  BOOST_REQUIRE_EQUAL( read_shape.size(), 2u );
  BOOST_CHECK_CLOSE( read_shape[0], fit.first, 1.0 );
  BOOST_CHECK_CLOSE( read_shape[1], fit.second, 1.0 );
}//BOOST_AUTO_TEST_CASE( fitGenieFwhmFromDrf )


BOOST_AUTO_TEST_CASE( defaultGenieFwhm )
{
  const pair<float,float> hpge = ExportSpecFileCAM::default_genie_fwhm( true );
  const pair<float,float> nai = ExportSpecFileCAM::default_genie_fwhm( false );

  // HPGe has to match what `CAMIO::AddDetectorType(...)` writes for a Ge detector, so an export
  //  with the FWHM option off and one left at the HPGe default agree.
  BOOST_CHECK_CLOSE( hpge.first, 1.0f, 0.1 );
  BOOST_CHECK_CLOSE( hpge.second, 0.035f, 0.1 );

  // NaI is a fit of the Genie equation form to InterSpec's nominal "NaI 3x3" resolution, so check
  //  the physics rather than the exact coefficients: ~6-8% FWHM at 661.7 keV, monotonically
  //  increasing, and never negative over the range NaI spectra are used over.
  const auto genie_fwhm = []( const pair<float,float> &c, const float energy ) -> float {
    return c.first + c.second*std::sqrt(energy);
  };

  const float nai_at_662 = genie_fwhm( nai, 661.657f );
  BOOST_CHECK_MESSAGE( (nai_at_662 > 0.05f*661.657f) && (nai_at_662 < 0.09f*661.657f),
      "Default NaI FWHM at 661.7 keV was " << nai_at_662 << " keV ("
      << 100.0f*nai_at_662/661.657f << "%), expected roughly 5-9%" );

  float prev = -1.0f;
  for( float energy = 60.0f; energy <= 2614.0f; energy += 50.0f )
  {
    const float fwhm = genie_fwhm( nai, energy );
    BOOST_CHECK_MESSAGE( fwhm > 0.0f, "Default NaI FWHM went non-positive (" << fwhm
                        << " keV) at " << energy << " keV" );
    BOOST_CHECK_MESSAGE( fwhm > prev, "Default NaI FWHM was not monotonic at " << energy << " keV" );
    prev = fwhm;
  }

  // ...and it should be far wider than the HPGe default, or the two got swapped.
  BOOST_CHECK( nai_at_662 > 10.0f*genie_fwhm(hpge, 661.657f) );
}//BOOST_AUTO_TEST_CASE( defaultGenieFwhm )


BOOST_AUTO_TEST_CASE( efficiencyPairsConvertToSplineAndRoundTrip )
{
  set_data_dir();

  auto drf = make_shared<DetectorPeakResponse>( "PairsDRF", "Pairs test" );
  vector<DetectorPeakResponse::EnergyEffPoint> eff_pairs;
  const vector<pair<float,float>> points = { {59.0f, 0.1f}, {200.0f, 0.3f}, {661.0f, 0.5f}, {1332.0f, 0.2f} };
  for( const auto &pt : points )
    eff_pairs.push_back( {pt.first, pt.second, {}} );
  drf->setEfficiencyPoints( eff_pairs, 5.0f, 0.0, DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );

  // Genie wants absolute efficiency, so a far-field DRF needs a distance...
  BOOST_CHECK_THROW( ExportSpecFileCAM::convert_efficiency_to_genie( *drf, -1.0 ), std::exception );

  const double distance = 100.0*PhysicalUnits::cm;
  const ExportSpecFileCAM::GenieEfficiencyResult result
                      = ExportSpecFileCAM::convert_efficiency_to_genie( *drf, distance );
  BOOST_CHECK( result.model == CAMInputOutput::CAMIO::EfficiencyModel::SPLINE );
  BOOST_REQUIRE_EQUAL( result.points.size(), points.size() );
  for( size_t i = 0; i < points.size(); ++i )
  {
    BOOST_CHECK_CLOSE( result.points[i].Energy, points[i].first, 0.1 );

    // ...which is the intrinsic efficiency times the fractional solid angle, not the intrinsic
    //  efficiency itself - writing the latter would over-state the efficiency by orders of
    //  magnitude, and Genie would report correspondingly wrong activities.
    const double expected = drf->efficiency( points[i].first, distance );
    BOOST_CHECK_CLOSE( result.points[i].Efficiency, static_cast<float>(expected), 0.1 );
    BOOST_CHECK_MESSAGE( result.points[i].Efficiency < points[i].second,
        "Absolute efficiency at 1 m should be well below the intrinsic efficiency" );
  }

  // Round-trip the GEOM block through CAMIO itself.
  CAMInputOutput::CAMIO writer;
  writer.AddEfficiencyModel( result.model );
  writer.AddEfficiencyPoints( result.points );
  writer.AddEnergyCalibration( {0.0f, 3.0f} );
  writer.AddSpectrum( vector<uint32_t>( 512, 10 ) );
  vector<unsigned char> &filedata = writer.CreateFile();
  BOOST_REQUIRE( !filedata.empty() );

  CAMInputOutput::CAMIO reader;
  BOOST_REQUIRE_NO_THROW( reader.ReadFile(filedata) );
  const vector<CAMInputOutput::EfficiencyPoint> &read_points = reader.GetEfficiencyPoints();
  BOOST_CHECK( reader.GetEfficiencyModel() == CAMInputOutput::CAMIO::EfficiencyModel::SPLINE );
  BOOST_REQUIRE_EQUAL( read_points.size(), points.size() );

  for( size_t i = 0; i < points.size(); ++i )
  {
    BOOST_CHECK_CLOSE( read_points[i].Energy, result.points[i].Energy, 1.0 );
    BOOST_CHECK_CLOSE( read_points[i].Efficiency, result.points[i].Efficiency, 1.0 );
  }
}//BOOST_AUTO_TEST_CASE( efficiencyPairsConvertToSplineAndRoundTrip )


BOOST_AUTO_TEST_CASE( efficiencyExpOfLogsConvertsToDualExactly )
{
  set_data_dir();

  auto drf = make_shared<DetectorPeakResponse>( "ExpLogDRF", "Exp-log test" );
  const vector<float> coefs{ -1.0f, 0.5f, -0.1f };
  drf->fromExpOfLogPowerSeries( coefs, {}, 0.0, 5.0f, 1000.0f /*MeV*/, 59.0f, 3000.0f,
                                DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );

  const double distance = 50.0*PhysicalUnits::cm;
  const ExportSpecFileCAM::GenieEfficiencyResult result
                      = ExportSpecFileCAM::convert_efficiency_to_genie( *drf, distance );
  BOOST_CHECK( result.model == CAMInputOutput::CAMIO::EfficiencyModel::DUAL );
  BOOST_REQUIRE( !result.points.empty() );

  // The sampled points should exactly match the DRF's own absolute efficiency evaluation,
  //  since no fitting is involved for this (exact, analytic) source form.
  for( const CAMInputOutput::EfficiencyPoint &pt : result.points )
  {
    const float expected = static_cast<float>( drf->efficiency( pt.Energy, distance ) );
    BOOST_CHECK_CLOSE( pt.Efficiency, expected, 0.1 );
  }
}//BOOST_AUTO_TEST_CASE( efficiencyExpOfLogsConvertsToDualExactly )


BOOST_AUTO_TEST_CASE( efficiencyFixedGeometryIgnoresDistance )
{
  set_data_dir();

  // A fixed-geometry DRF's "intrinsic" efficiency already *is* the absolute efficiency, so no
  //  solid-angle correction should be applied, and no distance should be required.
  auto drf = make_shared<DetectorPeakResponse>( "FixedGeomDRF", "Fixed geometry test" );
  const vector<float> coefs{ -1.0f, 0.5f, -0.1f };
  drf->fromExpOfLogPowerSeries( coefs, {}, 0.0, 5.0f, 1000.0f /*MeV*/, 59.0f, 3000.0f,
                                DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct );
  BOOST_REQUIRE( drf->isFixedGeometry() );

  const ExportSpecFileCAM::GenieEfficiencyResult result
                      = ExportSpecFileCAM::convert_efficiency_to_genie( *drf, -1.0 );
  BOOST_REQUIRE( !result.points.empty() );

  for( const CAMInputOutput::EfficiencyPoint &pt : result.points )
    BOOST_CHECK_CLOSE( pt.Efficiency, drf->intrinsicEfficiency( pt.Energy ), 0.1 );
}//BOOST_AUTO_TEST_CASE( efficiencyFixedGeometryIgnoresDistance )


BOOST_AUTO_TEST_CASE( omitEnergyCalibration )
{
  set_data_dir();

  SpecUtils::SpecFile spec;
  auto m = make_shared<SpecUtils::Measurement>();
  m->set_gamma_counts( make_shared<vector<float>>( 1024, 10.0f ), 100.0f, 110.0f );
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( 1024, {0.0f, 3.0f}, {} );
  m->set_energy_calibration( cal );
  spec.add_measurement( m, true );

  const auto read_back_cal = []( const string &filedata ) -> vector<float> {
    const vector<uint8_t> data( begin(filedata), end(filedata) );
    CAMInputOutput::CAMIO cam;
    cam.ReadFile( data );
    return cam.GetEnergyCalibration();
  };

  CAMInputOutput::CnfGenieExtras extras;
  extras.omit_energy_calibration = false;
  stringstream with_cal;
  BOOST_REQUIRE( spec.write_cnf( with_cal, {}, {}, &extras ) );
  const vector<float> written = read_back_cal( with_cal.str() );
  BOOST_REQUIRE( written.size() >= 2 );
  BOOST_CHECK_CLOSE( written[1], 3.0f, 0.1 );

  extras.omit_energy_calibration = true;
  stringstream without_cal;
  BOOST_REQUIRE( spec.write_cnf( without_cal, {}, {}, &extras ) );
  const vector<float> omitted = read_back_cal( without_cal.str() );
  for( const float coef : omitted )
    BOOST_CHECK_SMALL( coef, 1.0E-6f );
}//BOOST_AUTO_TEST_CASE( omitEnergyCalibration )


BOOST_AUTO_TEST_CASE( multiBlockLibraryRoundTrips )
{
  set_data_dir();

  // NUCL/NLINES records spill into multiple chained blocks past 29 nuclides / 125 lines, and the
  //  chain's "next block" pointers are indices into the file's block directory - so how many
  //  optional blocks (SAMP/SPEC/GEOM) precede the chain must not shift them.  Exercise both with
  //  and without a GEOM (efficiency) block.
  for( const bool with_efficiency : { false, true } )
  {
    for( const int num_nuc : { 3, 29, 30, 70 } )
    {
      SpecUtils::SpecFile spec;
      auto m = make_shared<SpecUtils::Measurement>();
      m->set_gamma_counts( make_shared<vector<float>>( 1024, 10.0f ), 100.0f, 110.0f );
      auto cal = make_shared<SpecUtils::EnergyCalibration>();
      cal->set_polynomial( 1024, {0.0f, 3.0f}, {} );
      m->set_energy_calibration( cal );
      spec.add_measurement( m, true );

      const int lines_per_nuc = 5;

      CAMInputOutput::CnfGenieExtras extras;
      for( int n = 0; n < num_nuc; ++n )
      {
        for( int l = 0; l < lines_per_nuc; ++l )
        {
          CAMInputOutput::CnfGenieExtras::LibraryLine ll;
          ll.nuclide_name = "Nuc" + std::to_string(n);
          ll.half_life_seconds = 1000.0f*(n + 1);
          ll.energy = 50.0f + 7.0f*n + 0.9f*l;
          ll.yield = 0.1f*(l + 1);
          extras.library_lines.push_back( ll );
        }
      }

      if( with_efficiency )
      {
        extras.eff_model = CAMInputOutput::CAMIO::EfficiencyModel::DUAL;
        for( int i = 0; i < 5; ++i )
        {
          CAMInputOutput::EfficiencyPoint pt;
          pt.Index = i;
          pt.Energy = 100.0f*(i + 1);
          pt.Efficiency = 0.01f*(i + 1);
          pt.EfficiencyUncertainty = 0.0f;
          extras.eff_points.push_back( pt );
        }
      }//if( with_efficiency )

      stringstream ss;
      BOOST_REQUIRE( spec.write_cnf( ss, {}, {}, &extras ) );

      const string filedata = ss.str();
      const vector<uint8_t> data( begin(filedata), end(filedata) );

      CAMInputOutput::CAMIO cam;
      BOOST_REQUIRE_NO_THROW( cam.ReadFile( data ) );

      // Note: `GetLines()` appends to its cache, and `GetNuclides()` populates that cache as a
      //  side effect - so read the lines first, and only once (see the gotcha on `GetLines()`).
      const size_t num_lines = cam.GetLines().size();
      const size_t num_nuclides = cam.GetNuclides().size();

      BOOST_CHECK_MESSAGE( num_nuclides == static_cast<size_t>(num_nuc),
          "Read " << num_nuclides << " nuclides, expected " << num_nuc
          << " (with_efficiency=" << with_efficiency << ")" );
      BOOST_CHECK_MESSAGE( num_lines == static_cast<size_t>(num_nuc*lines_per_nuc),
          "Read " << num_lines << " lines, expected " << (num_nuc*lines_per_nuc)
          << " (with_efficiency=" << with_efficiency << ")" );

      // Every chained NUCL/NLINES block must point at the block that actually follows it (or at
      //  nothing, if it ends the chain); the pointer is `0x2800 + <block directory index>`.
      const auto read_u32 = [&data]( const size_t offset ) -> uint32_t {
        uint32_t value = 0;
        memcpy( &value, &data[offset], sizeof(uint32_t) );
        return value;
      };
      const auto read_u16 = [&data]( const size_t offset ) -> uint16_t {
        uint16_t value = 0;
        memcpy( &value, &data[offset], sizeof(uint16_t) );
        return value;
      };

      const uint32_t nucl_id = 0x00012007, nlines_id = 0x00012008;
      vector<uint32_t> block_ids;
      for( size_t i = 0; i < 64; ++i )
      {
        const size_t offset = 0x70 + i*0x30;
        if( (offset + 0x30) > data.size() )
          break;
        const uint32_t id = read_u32( offset );
        if( id == 0 )
          break;
        block_ids.push_back( id );
      }

      for( size_t i = 0; i < block_ids.size(); ++i )
      {
        if( (block_ids[i] != nucl_id) && (block_ids[i] != nlines_id) )
          continue;

        const uint16_t link = read_u16( 0x70 + i*0x30 + 0x0E );
        BOOST_REQUIRE_MESSAGE( (link & 0xFF00) == 0x2800,
            "Unexpected chain link 0x" << std::hex << link << std::dec << " in block " << i );

        const size_t next = (link & 0x00FF);
        if( next == 0 )
          continue;  //end of the chain

        BOOST_CHECK_MESSAGE( next == (i + 1),
            "Block " << i << " points at block " << next << ", but the next block is " << (i + 1)
            << " (with_efficiency=" << with_efficiency << ", num_nuc=" << num_nuc << ")" );
        BOOST_CHECK_MESSAGE( (next < block_ids.size()) && (block_ids[next] == block_ids[i]),
            "Block " << i << " chains to a block of a different type" );
      }
    }//for( num_nuc )
  }//for( with_efficiency )
}//BOOST_AUTO_TEST_CASE( multiBlockLibraryRoundTrips )
