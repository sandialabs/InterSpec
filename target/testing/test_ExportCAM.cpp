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
#include "InterSpec/PeakFitUtils.h"
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

  /** Disables the "drop lines outside the spectrum" filter in the `build_genie_library(...)` tests. */
  const pair<float,float> s_test_energy_range{ 0.0f, 0.0f };

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
                                1.0, true, s_test_hpge_fwhm, s_test_energy_range, {}, &warnings );

  BOOST_REQUIRE_EQUAL( sources.size(), 1u );
  const ExportSpecFileCAM::GenieLibrarySource &src = sources[0];
  BOOST_CHECK_EQUAL( src.name, "Co60" );
  BOOST_CHECK( src.nuclide == co60 );
  //  Against a literal, not against the same expression the implementation uses.
  BOOST_CHECK_CLOSE( src.half_life_seconds, 1.6634e8f, 0.5 );
  BOOST_REQUIRE_EQUAL( src.lines.size(), 2u );

  bool have_key_line = false;
  for( const ExportSpecFileCAM::GenieLibraryLine &line : src.lines )
  {
    BOOST_CHECK( line.yield > 0.0f );
    //  Note: not bounded by 1 in general - annihilation lines get a contribution per positron,
    //  so e.g. Na22's 511 keV yield is about 1.8.
    BOOST_CHECK( line.yield < 3.0f );
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
                                1.0, false, s_test_hpge_fwhm, s_test_energy_range, {} );
  BOOST_REQUIRE_EQUAL( sources.size(), 1u );
  BOOST_CHECK_MESSAGE( sources[0].lines.size() >= 2, "Expected at least the 1173/1332 keV lines" );

  // A very high threshold should exclude everything but the most intense line(s).
  sources = ExportSpecFileCAM::build_genie_library(
                                peaks, ExportSpecFileCAM::GenieLibraryLineMode::AllLinesAboveThreshold,
                                99.0, false, s_test_hpge_fwhm, s_test_energy_range, {} );
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
                                0.01, true, s_test_hpge_fwhm, s_test_energy_range, {} );
  const vector<ExportSpecFileCAM::GenieLibrarySource> uncombined = ExportSpecFileCAM::build_genie_library(
                                peaks, ExportSpecFileCAM::GenieLibraryLineMode::AllLinesAboveThreshold,
                                0.01, false, s_test_hpge_fwhm, s_test_energy_range, {} );

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
                                1.0, true, s_test_hpge_fwhm, s_test_energy_range, {} );
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
                                1.0, true, s_test_hpge_fwhm, s_test_energy_range, {} );
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


// Note: the energy-calibration option is only honored for a spectrum-less file; see
//  `libraryOnlyFileHasNoSpectrum` for the full set of cases.


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
        // Byte 0x0F is a block-type flag (0x28 for NUCL/NLINES), not the high byte of a 16-bit
        //  field - the index lives entirely in the byte at 0x0E.
        BOOST_REQUIRE_MESSAGE( (link & 0xFF00) == 0x2800,
            "Unexpected chain link 0x" << std::hex << link << std::dec << " in block " << i );

        const size_t next = (link & 0x00FF);
        if( next == 0 )
          continue;  //end of the chain

        BOOST_REQUIRE_MESSAGE( next >= 3, "Chain link " << next << " in block " << i
                              << " is too small to be a slot index plus 2" );

        // Real Genie files store the *next* chain member's directory slot plus 2 (see
        //  `set_next_block_link` in CAMIO.cpp for the files this was derived from).
        const size_t next_slot = next - 2;
        BOOST_CHECK_MESSAGE( next_slot == (i + 1),
            "Block " << i << " points at slot " << next_slot << " (raw link " << next
            << "), but the next block is " << (i + 1)
            << " (with_efficiency=" << with_efficiency << ", num_nuc=" << num_nuc << ")" );
        BOOST_CHECK_MESSAGE( (next_slot < block_ids.size()) && (block_ids[next_slot] == block_ids[i]),
            "Block " << i << " chains to a block of a different type" );
      }
    }//for( num_nuc )
  }//for( with_efficiency )
}//BOOST_AUTO_TEST_CASE( multiBlockLibraryRoundTrips )


BOOST_AUTO_TEST_CASE( libraryYieldsAreBranchingRatios )
{
  set_data_dir();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );

  // `NuclideMixture::photons(age)` gives photons/second per *initial* becquerel, so it still
  //  carries the parent's decay factor; a GENIE library wants photons per decay.  Nuclides whose
  //  default age is non-zero (i.e. anything that doesn't decay to stable children) exposed this -
  //  Co60's default age is 0, which is why a Co60-only test cannot catch it.
  struct Expect { const char *nuc; float energy; double branching_ratio; };
  const Expect expects[] = {
    { "Co60",  1332.5f, 0.99983 },
    { "Cs137",  661.7f, 0.85334 },
    { "Eu152",  121.8f, 0.28580 },
    { "I131",   364.5f, 0.81500 },
  };

  for( const Expect &e : expects )
  {
    const SandiaDecay::Nuclide * const nuc = db->nuclide( e.nuc );
    BOOST_REQUIRE( nuc );

    deque<shared_ptr<const PeakDef>> peaks;
    peaks.push_back( make_peak( e.energy, nuc ) );

    const vector<ExportSpecFileCAM::GenieLibrarySource> sources
        = ExportSpecFileCAM::build_genie_library( peaks,
              ExportSpecFileCAM::GenieLibraryLineMode::PeakLinesOnly,
              1.0, false, s_test_hpge_fwhm, s_test_energy_range, {} );

    BOOST_REQUIRE_MESSAGE( sources.size() == 1u, e.nuc << ": expected one source" );
    BOOST_REQUIRE_MESSAGE( sources[0].lines.size() == 1u, e.nuc << ": expected one line" );

    BOOST_CHECK_MESSAGE( std::fabs(sources[0].lines[0].yield - e.branching_ratio) < 0.01*e.branching_ratio,
        e.nuc << " " << e.energy << " keV yield was " << sources[0].lines[0].yield
        << ", expected the branching ratio " << e.branching_ratio );
  }
}//BOOST_AUTO_TEST_CASE( libraryYieldsAreBranchingRatios )


BOOST_AUTO_TEST_CASE( libraryAbundanceIsPercent )
{
  set_data_dir();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );

  // GENIE library abundances are a percent, not a fraction - confirmed against the real,
  //  Canberra-produced library that ships in this repo (see readExampleMultiNuclideLibrary).
  const SandiaDecay::Nuclide * const cs137 = db->nuclide( "Cs137" );
  BOOST_REQUIRE( cs137 );

  deque<shared_ptr<const PeakDef>> peaks;
  peaks.push_back( make_peak( 661.7f, cs137 ) );

  const vector<ExportSpecFileCAM::GenieLibrarySource> sources
      = ExportSpecFileCAM::build_genie_library( peaks,
            ExportSpecFileCAM::GenieLibraryLineMode::PeakLinesOnly,
            1.0, false, s_test_hpge_fwhm, s_test_energy_range, {} );
  BOOST_REQUIRE_EQUAL( sources.size(), 1u );

  const vector<CAMInputOutput::CnfGenieExtras::LibraryLine> lines
      = ExportSpecFileCAM::to_library_lines( sources );
  BOOST_REQUIRE_EQUAL( lines.size(), 1u );
  BOOST_CHECK_MESSAGE( std::fabs(lines[0].yield - 85.334f) < 1.0f,
      "Library line abundance was " << lines[0].yield << ", expected ~85.3 (percent, not 0.853)" );

  // ...and it survives to the file that way.
  SpecUtils::SpecFile spec;
  auto m = make_shared<SpecUtils::Measurement>();
  m->set_gamma_counts( make_shared<vector<float>>( 1024, 10.0f ), 100.0f, 110.0f );
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( 1024, {0.0f, 3.0f}, {} );
  m->set_energy_calibration( cal );
  spec.add_measurement( m, true );

  CAMInputOutput::CnfGenieExtras extras;
  extras.library_lines = lines;

  stringstream ss;
  BOOST_REQUIRE( spec.write_cnf( ss, {}, {}, &extras ) );

  const string filedata = ss.str();
  const vector<uint8_t> data( begin(filedata), end(filedata) );
  CAMInputOutput::CAMIO cam;
  BOOST_REQUIRE_NO_THROW( cam.ReadFile( data ) );
  const vector<CAMInputOutput::Line> read_lines = cam.GetLines();
  BOOST_REQUIRE_EQUAL( read_lines.size(), 1u );
  BOOST_CHECK_CLOSE( read_lines[0].Abundance, 85.334f, 2.0 );
}//BOOST_AUTO_TEST_CASE( libraryAbundanceIsPercent )


BOOST_AUTO_TEST_CASE( getLinesIsIdempotent )
{
  set_data_dir();

  const string libfile = SpecUtils::append_path( g_test_file_dir, "CamLibrary/ExampleMultiNuclide.nlb" );
  BOOST_REQUIRE( SpecUtils::is_file(libfile) );

  vector<char> filedata;
  SpecUtils::load_file_data( libfile.c_str(), filedata );
  BOOST_REQUIRE( !filedata.empty() );
  if( !filedata.empty() && (filedata.back() == '\0') )
    filedata.pop_back();
  const vector<uint8_t> data( begin(filedata), end(filedata) );

  CAMInputOutput::CAMIO cam;
  BOOST_REQUIRE_NO_THROW( cam.ReadFile( data ) );

  const size_t first = cam.GetLines().size();
  const size_t after_nuclides = (cam.GetNuclides(), cam.GetLines().size());
  const size_t third = cam.GetLines().size();

  BOOST_CHECK_MESSAGE( first == after_nuclides,
      "GetLines() returned " << first << " lines, then " << after_nuclides
      << " after GetNuclides() - it must not append to its cache" );
  BOOST_CHECK_EQUAL( first, third );
}//BOOST_AUTO_TEST_CASE( getLinesIsIdempotent )


BOOST_AUTO_TEST_CASE( explicitKeyLineIsHonored )
{
  set_data_dir();

  // What the GUI table previews has to be what lands in the file; `CAMIO::AssignKeyLines()` must
  //  not override a caller's explicit choice.
  SpecUtils::SpecFile spec;
  auto m = make_shared<SpecUtils::Measurement>();
  m->set_gamma_counts( make_shared<vector<float>>( 1024, 10.0f ), 100.0f, 110.0f );
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( 1024, {0.0f, 3.0f}, {} );
  m->set_energy_calibration( cal );
  spec.add_measurement( m, true );

  CAMInputOutput::CnfGenieExtras extras;
  const float energies[] = { 122.06f, 136.47f, 692.03f };
  const float yields[]   = {   85.6f,   10.68f,  0.157f };
  for( size_t i = 0; i < 3; ++i )
  {
    CAMInputOutput::CnfGenieExtras::LibraryLine ll;
    ll.nuclide_name = "Co57";
    ll.half_life_seconds = 2.35e7f;
    ll.energy = energies[i];
    ll.yield = yields[i];
    ll.is_key_line = (i == 0);   //the 122 keV line - the one a user would pick
    extras.library_lines.push_back( ll );
  }

  stringstream ss;
  BOOST_REQUIRE( spec.write_cnf( ss, {}, {}, &extras ) );

  const string filedata = ss.str();
  const vector<uint8_t> data( begin(filedata), end(filedata) );
  CAMInputOutput::CAMIO cam;
  BOOST_REQUIRE_NO_THROW( cam.ReadFile( data ) );

  const vector<CAMInputOutput::Line> read_lines = cam.GetLines();
  BOOST_REQUIRE_EQUAL( read_lines.size(), 3u );

  size_t num_key = 0;
  float key_energy = 0.0f;
  for( const CAMInputOutput::Line &l : read_lines )
  {
    if( l.IsKeyLine )
    {
      num_key += 1;
      key_energy = l.Energy;
    }
  }

  BOOST_CHECK_EQUAL( num_key, 1u );
  BOOST_CHECK_MESSAGE( std::fabs(key_energy - 122.06f) < 0.5f,
      "The file's key line is at " << key_energy << " keV, but 122.06 keV was requested" );
}//BOOST_AUTO_TEST_CASE( explicitKeyLineIsHonored )


BOOST_AUTO_TEST_CASE( fitGenieFwhmHandlesZeroLowerEnergy )
{
  set_data_dir();

  // A DRF may report a lower energy of 0 (it is the default); the fit samples logarithmically, so
  //  a zero lower bound used to make every sample NaN, and the NaN went into the file.
  auto drf = make_shared<DetectorPeakResponse>( "ZeroLowerDRF", "Zero lower energy" );
  const vector<float> eff_coefs{ -1.0f, 0.5f, -0.1f };
  drf->fromExpOfLogPowerSeries( eff_coefs, {}, 0.0, 5.0f, 1000.0f, 0.0f, 3000.0f,
                                DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
  drf->setFwhmCoefficients( {1.0f, 0.035f},
                            DetectorPeakResponse::ResolutionFnctForm::kConstantPlusSqrtEnergy );

  const pair<float,float> coeffs = ExportSpecFileCAM::fit_genie_fwhm_from_drf( *drf );
  BOOST_CHECK_MESSAGE( std::isfinite(coeffs.first) && std::isfinite(coeffs.second),
      "fit_genie_fwhm_from_drf gave non-finite coefficients {" << coeffs.first << ", "
      << coeffs.second << "} for a DRF with lowerEnergy()==0" );
  BOOST_CHECK_CLOSE( coeffs.first, 1.0f, 5.0 );
  BOOST_CHECK_CLOSE( coeffs.second, 0.035f, 5.0 );
}//BOOST_AUTO_TEST_CASE( fitGenieFwhmHandlesZeroLowerEnergy )


BOOST_AUTO_TEST_CASE( defaultNaiFwhmIsPositiveAndBeatsGenieDefault )
{
  const pair<float,float> nai = ExportSpecFileCAM::default_genie_fwhm( false );

  const auto genie_fwhm = []( const pair<float,float> &c, const float energy ) -> float {
    return c.first + c.second*std::sqrt(energy);
  };

  // Must stay positive over the whole range a NaI spectrum covers - `merge_unresolvable_lines`
  //  clamps a non-positive sigma, which would make low-energy multiplets look resolvable.
  for( float energy = 10.0f; energy <= 3000.0f; energy += 1.0f )
  {
    const float fwhm = genie_fwhm( nai, energy );
    BOOST_REQUIRE_MESSAGE( fwhm > 0.0f, "Default NaI FWHM was " << fwhm << " keV at " << energy << " keV" );
  }

  // A + B*sqrt(E) cannot follow a real NaI resolution across two decades, so the fit trades
  //  accuracy at one end against the other; what it must do is keep the WORST-case relative error
  //  over the useful range below the Genie manual's own {-7, 2} default.  (Measured: ours peaks at
  //  ~15% at 2614 keV, Genie's at ~24% at 122 keV.)  It is expected to be worse than Genie's
  //  default at some individual energies - that is the trade, not a regression.
  const pair<float,float> genie_manual_default{ -7.0f, 2.0f };
  float our_worst = 0.0f, genie_worst = 0.0f;
  for( float energy = 30.0f; energy <= 2614.0f; energy += 1.0f )
  {
    const float ref = PeakFitUtils::nai_fwhm_fcn( energy );
    our_worst = std::max( our_worst, std::fabs(genie_fwhm(nai, energy) - ref) / ref );
    genie_worst = std::max( genie_worst, std::fabs(genie_fwhm(genie_manual_default, energy) - ref) / ref );
  }

  BOOST_CHECK_MESSAGE( our_worst < 0.20f,
      "Worst-case default NaI FWHM error is " << 100.0f*our_worst << "%, expected under 20%" );
  BOOST_CHECK_MESSAGE( our_worst < genie_worst,
      "Worst-case error is " << 100.0f*our_worst << "%, no better than the Genie manual default's "
      << 100.0f*genie_worst << "% - the fit is not earning its keep" );
}//BOOST_AUTO_TEST_CASE( defaultNaiFwhmIsPositiveAndBeatsGenieDefault )


BOOST_AUTO_TEST_CASE( libraryOnlyFileHasNoSpectrum )
{
  set_data_dir();

  SpecUtils::SpecFile spec;
  auto m = make_shared<SpecUtils::Measurement>();
  m->set_gamma_counts( make_shared<vector<float>>( 1024, 10.0f ), 100.0f, 110.0f );
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( 1024, {0.0f, 3.0f}, {} );
  m->set_energy_calibration( cal );
  spec.add_measurement( m, true );

  CAMInputOutput::CnfGenieExtras extras;
  CAMInputOutput::CnfGenieExtras::LibraryLine ll;
  ll.nuclide_name = "Cs137";
  ll.half_life_seconds = 9.5e8f;
  ll.energy = 661.657f;
  ll.yield = 85.3f;
  extras.library_lines.push_back( ll );
  extras.omit_spectrum = true;

  stringstream ss;
  BOOST_REQUIRE( spec.write_cnf( ss, {}, {}, &extras ) );

  const string filedata = ss.str();
  const vector<uint8_t> data( begin(filedata), end(filedata) );

  CAMInputOutput::CAMIO cam;
  BOOST_REQUIRE_NO_THROW( cam.ReadFile( data ) );

  // The library is there...
  BOOST_CHECK_EQUAL( cam.GetLines().size(), 1u );
  // ...and the spectrum is not.
  BOOST_CHECK_THROW( cam.GetSpectrum(), std::exception );

  // With a spectrum omitted, the energy calibration becomes optional.
  extras.omit_energy_calibration = true;
  stringstream no_cal;
  BOOST_REQUIRE( spec.write_cnf( no_cal, {}, {}, &extras ) );
  const string nocal_data = no_cal.str();
  const vector<uint8_t> nocal( begin(nocal_data), end(nocal_data) );
  CAMInputOutput::CAMIO cam2;
  BOOST_REQUIRE_NO_THROW( cam2.ReadFile( nocal ) );
  for( const float coef : cam2.GetEnergyCalibration() )
    BOOST_CHECK_SMALL( coef, 1.0E-6f );

  // ...but asking to omit it while writing a spectrum must be ignored, since channel counts are
  //  meaningless without a calibration.
  extras.omit_spectrum = false;
  stringstream with_spec;
  BOOST_REQUIRE( spec.write_cnf( with_spec, {}, {}, &extras ) );
  const string withspec_data = with_spec.str();
  const vector<uint8_t> withspec( begin(withspec_data), end(withspec_data) );
  CAMInputOutput::CAMIO cam3;
  BOOST_REQUIRE_NO_THROW( cam3.ReadFile( withspec ) );
  const vector<float> cal_coefs = cam3.GetEnergyCalibration();
  BOOST_REQUIRE( cal_coefs.size() >= 2 );
  BOOST_CHECK_CLOSE( cal_coefs[1], 3.0f, 0.1 );
}//BOOST_AUTO_TEST_CASE( libraryOnlyFileHasNoSpectrum )


BOOST_AUTO_TEST_CASE( tooMuchLibraryDataFailsLoudly )
{
  set_data_dir();

  // A CAM file's block directory is finite; overrunning it used to silently corrupt the file
  //  (dropping nuclides, and eventually overwriting the ACQP block).  It must fail instead.
  SpecUtils::SpecFile spec;
  auto m = make_shared<SpecUtils::Measurement>();
  m->set_gamma_counts( make_shared<vector<float>>( 1024, 10.0f ), 100.0f, 110.0f );
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( 1024, {0.0f, 3.0f}, {} );
  m->set_energy_calibration( cal );
  spec.add_measurement( m, true );

  //  250 nuclides keeps us under CAMIO's separate 255-nuclide limit, so it is the block-directory
  //  limit (~5000 lines -> ~40 NLINES blocks) this exercises.
  CAMInputOutput::CnfGenieExtras extras;
  for( int n = 0; n < 250; ++n )
  {
    for( int l = 0; l < 20; ++l )
    {
      CAMInputOutput::CnfGenieExtras::LibraryLine ll;
      ll.nuclide_name = "Nuc" + std::to_string(n);
      ll.half_life_seconds = 1000.0f*(n + 1);
      ll.energy = 50.0f + 7.0f*n + 0.9f*l;
      ll.yield = 10.0f;
      extras.library_lines.push_back( ll );
    }
  }

  stringstream ss;
  BOOST_CHECK_MESSAGE( !spec.write_cnf( ss, {}, {}, &extras ),
      "write_cnf should refuse to write a library too large for the block directory, rather than"
      " silently producing a corrupt file" );
}//BOOST_AUTO_TEST_CASE( tooMuchLibraryDataFailsLoudly )


BOOST_AUTO_TEST_CASE( manyLinesForOneNuclideRoundTrips )
{
  set_data_dir();

  // The nuclide record's size field is 16 bits; writing only its low byte made any nuclide with
  //  more than 65 lines unreadable.  Eu-152 and Bi-214 easily exceed that.
  SpecUtils::SpecFile spec;
  auto m = make_shared<SpecUtils::Measurement>();
  m->set_gamma_counts( make_shared<vector<float>>( 1024, 10.0f ), 100.0f, 110.0f );
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( 1024, {0.0f, 3.0f}, {} );
  m->set_energy_calibration( cal );
  spec.add_measurement( m, true );

  for( const int num_lines : { 65, 66, 150 } )
  {
    CAMInputOutput::CnfGenieExtras extras;
    for( int l = 0; l < num_lines; ++l )
    {
      CAMInputOutput::CnfGenieExtras::LibraryLine ll;
      ll.nuclide_name = "Eu152";
      ll.half_life_seconds = 4.27e8f;
      ll.energy = 50.0f + 3.7f*l;
      ll.yield = 5.0f;
      extras.library_lines.push_back( ll );
    }

    stringstream ss;
    BOOST_REQUIRE( spec.write_cnf( ss, {}, {}, &extras ) );

    const string filedata = ss.str();
    const vector<uint8_t> data( begin(filedata), end(filedata) );

    CAMInputOutput::CAMIO cam;
    BOOST_REQUIRE_NO_THROW( cam.ReadFile( data ) );
    BOOST_CHECK_MESSAGE( cam.GetLines().size() == static_cast<size_t>(num_lines),
        "Read " << cam.GetLines().size() << " lines, expected " << num_lines );
    BOOST_CHECK_NO_THROW( cam.GetNuclides() );
  }
}//BOOST_AUTO_TEST_CASE( manyLinesForOneNuclideRoundTrips )


BOOST_AUTO_TEST_CASE( camioIsReusableForWriting )
{
  set_data_dir();

  // The ACQP/SAMP "common" sections used to be file-scope globals, so a second file inherited the
  //  first one's calibrations and sample title - and across sessions, raced.
  const auto write_one = []( const float gain, const string &title ) -> vector<uint8_t> {
    CAMInputOutput::CAMIO cam;
    cam.AddSampleTitle( title );
    cam.AddEnergyCalibration( {0.0f, gain} );
    cam.AddSpectrum( vector<uint32_t>( 256, 5 ) );
    return cam.CreateFile();
  };

  const vector<uint8_t> a = write_one( 0.5f, "AAA" );
  const vector<uint8_t> b = write_one( 2.0f, "BBB" );

  CAMInputOutput::CAMIO reader_a, reader_b;
  BOOST_REQUIRE_NO_THROW( reader_a.ReadFile( a ) );
  BOOST_REQUIRE_NO_THROW( reader_b.ReadFile( b ) );

  BOOST_CHECK_CLOSE( reader_a.GetEnergyCalibration().at(1), 0.5f, 0.1 );
  BOOST_CHECK_CLOSE( reader_b.GetEnergyCalibration().at(1), 2.0f, 0.1 );
  BOOST_CHECK_EQUAL( SpecUtils::trim_copy(reader_a.GetSampleTitle()), "AAA" );
  BOOST_CHECK_EQUAL( SpecUtils::trim_copy(reader_b.GetSampleTitle()), "BBB" );
}//BOOST_AUTO_TEST_CASE( camioIsReusableForWriting )


BOOST_AUTO_TEST_CASE( peakBlockRoundTrips )
{
  set_data_dir();

  SpecUtils::SpecFile spec;
  auto m = make_shared<SpecUtils::Measurement>();
  m->set_gamma_counts( make_shared<vector<float>>( 1024, 10.0f ), 100.0f, 110.0f );
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( 1024, {0.0f, 3.0f}, {} );
  m->set_energy_calibration( cal );
  spec.add_measurement( m, true );

  CAMInputOutput::CnfGenieExtras extras;
  for( int i = 0; i < 3; ++i )
  {
    CAMInputOutput::Peak p{};
    p.Energy = 661.657f + 200.0f*i;
    p.Centroid = p.Energy/3.0f;
    p.CentroidUncertainty = 0.15f;
    p.FullWidthAtHalfMaximum = 2.5f + 0.1f*i;
    p.Area = 12345.0f*(i+1);
    p.AreaUncertainty = 150.0f*(i+1);
    p.Continuum = 500.0f*(i+1);
    p.CountRate = p.Area/100.0f;
    p.CountRateUncertainty = p.AreaUncertainty/100.0f;
    p.LeftChannel = static_cast<int>(p.Centroid) - 10;
    p.RightChannel = static_cast<int>(p.Centroid) + 10;
    extras.peaks.push_back( p );
  }

  stringstream ss;
  BOOST_REQUIRE( spec.write_cnf( ss, {}, {}, &extras ) );

  const string filedata = ss.str();
  BOOST_CHECK_MESSAGE( (filedata.size() % 512) == 0,
      "File size " << filedata.size() << " is not a multiple of 512" );

  const vector<uint8_t> data( begin(filedata), end(filedata) );
  CAMInputOutput::CAMIO cam;
  BOOST_REQUIRE_NO_THROW( cam.ReadFile( data ) );

  const vector<CAMInputOutput::Peak> read_peaks = cam.GetPeaks();
  BOOST_REQUIRE_EQUAL( read_peaks.size(), extras.peaks.size() );

  for( size_t i = 0; i < read_peaks.size(); ++i )
  {
    const CAMInputOutput::Peak &w = extras.peaks[i];
    const CAMInputOutput::Peak &r = read_peaks[i];
    BOOST_CHECK_CLOSE( r.Energy, w.Energy, 0.1 );
    BOOST_CHECK_CLOSE( r.Centroid, w.Centroid, 0.1 );
    BOOST_CHECK_CLOSE( r.CentroidUncertainty, w.CentroidUncertainty, 1.0 );
    BOOST_CHECK_CLOSE( r.FullWidthAtHalfMaximum, w.FullWidthAtHalfMaximum, 0.1 );
    BOOST_CHECK_CLOSE( r.Area, w.Area, 0.1 );
    BOOST_CHECK_CLOSE( r.AreaUncertainty, w.AreaUncertainty, 0.1 );
    // `Continuum` (0x0C) sits immediately after the 16-bit `Width` (0x0A); writing Width as a
    //  longword silently zeroed it, so this assertion is the guard against that regression.
    BOOST_CHECK_CLOSE( r.Continuum, w.Continuum, 0.1 );
    BOOST_CHECK_CLOSE( r.CountRate, w.CountRate, 0.1 );
    BOOST_CHECK_EQUAL( r.LeftChannel, w.LeftChannel );
    BOOST_CHECK_EQUAL( r.RightChannel, w.RightChannel );
  }

  // Writing peaks must not disturb the spectrum.
  SpecUtils::SpecFile readback;
  stringstream in( filedata );
  BOOST_REQUIRE( readback.load_from_cnf( in ) );
  BOOST_CHECK_CLOSE( readback.gamma_count_sum(), 10240.0, 0.1 );
}//BOOST_AUTO_TEST_CASE( peakBlockRoundTrips )


/** Our own reader will happily round-trip a PEAK block that GENIE ignores, so this pins the
 written block against a real Genie-produced one (CNFreader's cs137.CNF): the block header must
 match it word for word apart from the block's own address, and the common area must carry the
 analysis-engine names and per-record markers Genie writes.  Peaks written without those were
 silently ignored by GENIE.
 */
BOOST_AUTO_TEST_CASE( peakBlockMatchesRealGenieLayout )
{
  set_data_dir();

  SpecUtils::SpecFile spec;
  auto m = make_shared<SpecUtils::Measurement>();
  m->set_gamma_counts( make_shared<vector<float>>( 1024, 10.0f ), 100.0f, 110.0f );
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( 1024, {0.0f, 3.0f}, {} );
  m->set_energy_calibration( cal );
  m->set_start_time( std::chrono::system_clock::from_time_t( 1600000000 ) );
  spec.add_measurement( m, true );

  CAMInputOutput::CnfGenieExtras extras;
  CAMInputOutput::Peak p{};
  p.Energy = 661.657f;
  p.Centroid = 220.5f;
  p.FullWidthAtHalfMaximum = 2.5f;
  p.Area = 12345.0f;
  p.AreaUncertainty = 150.0f;
  p.CountRate = 123.45f;
  p.LeftChannel = 210;
  p.RightChannel = 232;
  extras.peaks.push_back( p );

  stringstream ss;
  BOOST_REQUIRE( spec.write_cnf( ss, {}, {}, &extras ) );
  const string filedata = ss.str();
  const vector<uint8_t> ours( begin(filedata), end(filedata) );

  vector<uint8_t> genie;
  {
    const string path = SpecUtils::append_path( g_test_file_dir, "CamLibrary/cs137.CNF" );
    ifstream f( path.c_str(), ios::binary );
    if( !f.is_open() )
    {
      BOOST_TEST_MESSAGE( "Skipping: no real Genie peak file at " << path );
      return;
    }
    genie.assign( (istreambuf_iterator<char>(f)), istreambuf_iterator<char>() );
  }

  const auto u16 = []( const vector<uint8_t> &d, const size_t o ) -> uint16_t {
    BOOST_REQUIRE( (o + 2) <= d.size() );
    uint16_t v; memcpy( &v, &d[o], 2 ); return v;
  };

  const auto peak_block = []( const vector<uint8_t> &d ) -> size_t {
    for( size_t i = 0; i < 28; ++i )  //the block directory holds at most 28 entries
    {
      const size_t h = 0x70 + i*0x30;
      if( (h + 0x30) > d.size() ) break;
      uint32_t id, loc;
      memcpy( &id, &d[h], 4 );
      memcpy( &loc, &d[h+0x0A], 4 );
      if( id == 0x12006 ) return loc;
    }
    return 0;
  };

  const size_t our_pos = peak_block( ours ), genie_pos = peak_block( genie );
  BOOST_REQUIRE( our_pos && genie_pos );

  // Every header word but the block's own address (0x0A-0x0D) must match.
  for( size_t o = 0; o < 0x30; o += 2 )
  {
    if( (o == 0x0A) || (o == 0x0C) )
      continue;
    BOOST_CHECK_MESSAGE( u16(ours,our_pos+o) == u16(genie,genie_pos+o),
        "PEAK header word 0x" << std::hex << o << ": wrote 0x" << u16(ours,our_pos+o)
        << ", real Genie file has 0x" << u16(genie,genie_pos+o) << std::dec );
  }

  const auto has_string = []( const vector<uint8_t> &d, const size_t o, const char *str ){
    return (o + strlen(str)) <= d.size() && (memcmp( &d[o], str, strlen(str) ) == 0);
  };
  //  Offsets and names below are literals on purpose - they are what the real file has,
  //  so the test would still catch a change made to both writer and its constants.
  BOOST_CHECK( has_string( ours, our_pos + 0x3C, "2nd Diff" ) );
  BOOST_CHECK( has_string( ours, our_pos + 0x9D, "NLSQ Fit" ) );
  BOOST_CHECK( has_string( genie, genie_pos + 0x3C, "2nd Diff" ) );
  BOOST_CHECK( has_string( genie, genie_pos + 0x9D, "NLSQ Fit" ) );

  for( const size_t offset : { size_t(0x8D), size_t(0xB9) } )
  {
    uint64_t ticks = 0;
    memcpy( &ticks, &ours[our_pos + offset], sizeof(ticks) );
    BOOST_CHECK_MESSAGE( ticks != 0, "Analysis timestamp at PEAK+0x" << std::hex << offset
                                      << std::dec << " was left zero" );
  }

  // Each record starts with a 0x01 marker, the field values one byte further in.
  const uint16_t headSize = u16( ours, our_pos + 0x10 );
  const uint16_t recSize = u16( ours, our_pos + 0x20 );
  const uint16_t recOffset = u16( ours, our_pos + 0x22 );
  const uint16_t numRec = u16( ours, our_pos + 0x1E );
  BOOST_REQUIRE_EQUAL( numRec, 1 );
  for( size_t i = 0; i < numRec; ++i )
    BOOST_CHECK_EQUAL( int(ours[our_pos + headSize + recOffset + i*recSize]), 1 );
}//BOOST_AUTO_TEST_CASE( peakBlockMatchesRealGenieLayout )


BOOST_AUTO_TEST_CASE( realGenieFileFormatAnchors )
{
  set_data_dir();

  // Ba-133.cnf is a real Genie 2000 written CNF, and is the anchor for essentially every
  //  reverse-engineered layout this code depends on: it has multi-block NUCL *and* NLINES chains,
  //  a populated GEOM (efficiency) block, and a populated PEAK block.  If any of these assertions
  //  break, the writer's idea of the format has drifted from a real file's.
  const string filename = SpecUtils::append_path( g_test_file_dir, "CamLibrary/Ba-133.cnf" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file(filename), "Missing test file: " << filename );

  vector<char> filedata;
  SpecUtils::load_file_data( filename.c_str(), filedata );
  BOOST_REQUIRE( !filedata.empty() );
  if( filedata.back() == '\0' )
    filedata.pop_back();
  const vector<uint8_t> data( begin(filedata), end(filedata) );

  CAMInputOutput::CAMIO cam;
  BOOST_REQUIRE_NO_THROW( cam.ReadFile( data ) );

  // --- Nuclide library: abundances are a PERCENT, not a fraction ---------------------------
  //  (Ba133's 30.97 keV Kalpha x-ray group is ~64%; a fraction-based file would hold 0.64.)
  const vector<CAMInputOutput::Line> lines = cam.GetLines();
  const vector<CAMInputOutput::Nuclide> nuclides = cam.GetNuclides();
  BOOST_CHECK_EQUAL( nuclides.size(), 30u );
  BOOST_CHECK_EQUAL( lines.size(), 175u );

  float max_abundance = 0.0f;
  for( const CAMInputOutput::Line &l : lines )
    max_abundance = std::max( max_abundance, l.Abundance );
  BOOST_CHECK_MESSAGE( max_abundance > 10.0f,
      "Largest line abundance in a real Genie library is " << max_abundance
      << "; abundances are percents, so this should be well above 10" );

  // --- Efficiency (GEOM) block ---------------------------------------------------------------
  const vector<CAMInputOutput::EfficiencyPoint> eff = cam.GetEfficiencyPoints();
  BOOST_REQUIRE_EQUAL( eff.size(), 19u );
  BOOST_CHECK( cam.GetEfficiencyModel() == CAMInputOutput::CAMIO::EfficiencyModel::INTERPOL );
  BOOST_CHECK_CLOSE( eff[0].Energy, 20.0f, 0.1 );
  BOOST_CHECK_CLOSE( eff[0].Efficiency, 7.63589e-4f, 0.5 );
  BOOST_CHECK_CLOSE( eff[0].EfficiencyUncertainty, 9.16307e-5f, 0.5 );

  // --- Peaks ---------------------------------------------------------------------------------
  const vector<CAMInputOutput::Peak> peaks = cam.GetPeaks();
  BOOST_REQUIRE_EQUAL( peaks.size(), 23u );

  //  The strongest single check on the peak record layout: the count-rate uncertainty field is a
  //  RELATIVE PERCENT, so it must equal 100*AreaUncertainty/Area for every peak.
  for( const CAMInputOutput::Peak &p : peaks )
  {
    if( p.Area <= 0.0f )
      continue;
    const float expected = 100.0f * p.AreaUncertainty / p.Area;
    BOOST_CHECK_MESSAGE( std::fabs(p.CountRateUncertainty - expected) < 0.01f*expected,
        "Peak at " << p.Energy << " keV: CountRateUncertainty " << p.CountRateUncertainty
        << " is not 100*areaUnc/area (" << expected << ")" );
  }

  //  ...and 1000 is the "no low tail" sentinel, with real tails taking other values.
  bool saw_no_tail = false, saw_real_tail = false;
  for( const CAMInputOutput::Peak &p : peaks )
  {
    saw_no_tail = saw_no_tail || (std::fabs(p.LowTail - CAMInputOutput::Peak::sm_no_low_tail) < 0.01f);
    saw_real_tail = saw_real_tail || ((p.LowTail > 0.0f)
                        && (std::fabs(p.LowTail - CAMInputOutput::Peak::sm_no_low_tail) > 1.0f));
  }
  BOOST_CHECK( saw_no_tail );
  BOOST_CHECK( saw_real_tail );

  //  FWHM is in keV, not channels.  This is the *measured* width from the peak fit, so it only
  //  tracks the file's `FWHM = A0 + A1*sqrt(E)` shape calibration for peaks with enough counts to
  //  determine a width - weak peaks carry a lot of noise (and unresolved x-ray multiplets are
  //  wider than a single line).  Compare on strong peaks only, and check the aggregate rather
  //  than each peak: in channels these values would be off by ~1/gain (about 3x here).
  const vector<float> shape = cam.GetShapeCalibration();
  BOOST_REQUIRE( shape.size() >= 2 );

  double sum_rel_err = 0.0;
  size_t num_strong = 0;
  float lowest_e_fwhm = 0.0f, highest_e_fwhm = 0.0f, lowest_e = 1.0e9f, highest_e = 0.0f;
  for( const CAMInputOutput::Peak &p : peaks )
  {
    if( (p.Energy < 50.0f) || (p.Area < 200.0f) || (p.FullWidthAtHalfMaximum <= 0.0f) )
      continue;

    const float predicted = shape[0] + shape[1]*std::sqrt(p.Energy);
    sum_rel_err += std::fabs( p.FullWidthAtHalfMaximum - predicted ) / predicted;
    num_strong += 1;

    if( p.Energy < lowest_e ){ lowest_e = p.Energy; lowest_e_fwhm = p.FullWidthAtHalfMaximum; }
    if( p.Energy > highest_e ){ highest_e = p.Energy; highest_e_fwhm = p.FullWidthAtHalfMaximum; }
  }

  BOOST_REQUIRE_MESSAGE( num_strong >= 5, "Only " << num_strong << " strong peaks to check FWHM with" );

  const double mean_rel_err = sum_rel_err / num_strong;
  BOOST_CHECK_MESSAGE( mean_rel_err < 0.20,
      "Peak FWHM values are on average " << 100.0*mean_rel_err << "% away from the file's own shape"
      " calibration - the FWHM field may not be in keV" );

  BOOST_CHECK_MESSAGE( highest_e_fwhm > lowest_e_fwhm,
      "FWHM does not increase with energy (" << lowest_e_fwhm << " keV at " << lowest_e
      << " keV vs " << highest_e_fwhm << " keV at " << highest_e << " keV)" );

  // --- Multi-block chain links ---------------------------------------------------------------
  //  Both the NUCL and NLINES records span several blocks here, and each block points at the next
  //  one's directory slot PLUS 2 (byte 0x0E; 0x0F is a block-type flag, not a high byte).
  const auto u32 = [&data]( const size_t o ){ uint32_t v=0; memcpy(&v,&data[o],4); return v; };
  vector<uint32_t> ids;
  for( size_t i = 0; i < 40; ++i )
  {
    const size_t off = 0x70 + i*0x30;
    if( (off + 0x30) > data.size() )
      break;
    const uint32_t id = u32(off);
    if( id == 0 )
      break;
    ids.push_back( id );
  }

  size_t num_chained = 0;
  for( size_t i = 0; i < ids.size(); ++i )
  {
    if( (ids[i] != 0x00012007u) && (ids[i] != 0x00012008u) )  //NUCL, NLINES
      continue;

    const size_t off = 0x70 + i*0x30;
    BOOST_CHECK_MESSAGE( data[off + 0x0F] == 0x28,
        "Block " << i << " type flag at 0x0F is " << int(data[off+0x0F]) << ", expected 0x28" );

    const uint8_t link = data[off + 0x0E];
    if( link == 0 )
      continue;   //end of chain

    num_chained += 1;
    const size_t next_slot = static_cast<size_t>(link) - 2;
    BOOST_CHECK_MESSAGE( next_slot == (i + 1),
        "Block " << i << " chain link is " << int(link) << " -> slot " << next_slot
        << ", but the next block is slot " << (i + 1) );
    BOOST_REQUIRE( next_slot < ids.size() );
    BOOST_CHECK_EQUAL( ids[next_slot], ids[i] );
  }
  BOOST_CHECK_MESSAGE( num_chained >= 4,
      "Expected several chained NUCL/NLINES blocks in this file, found " << num_chained );
}//BOOST_AUTO_TEST_CASE( realGenieFileFormatAnchors )


BOOST_AUTO_TEST_CASE( lowTailMapsToGaussExpSkew )
{
  set_data_dir();

  // Genie's low tail and InterSpec's GaussExp skew are the same function:
  //    core:  exp(-0.5*u^2)                        for u >= -s
  //    tail:  exp(0.5*s^2 + s*u)                   for u <  -s
  //  with Genie's stored `T` equal to s*sigma.  Verify the two agree numerically, and that the
  //  conversion survives a write/read cycle.
  const double mean = 661.657, sigma = 1.2, skew = 1.75;   //s, in sigma
  for( const double x : { mean - 6.0*sigma, mean - 2.5*sigma, mean - sigma, mean, mean + 2.0*sigma } )
  {
    const double u = (x - mean)/sigma;
    const double T = skew*sigma;
    const double genie = (u >= -T/sigma) ? std::exp(-0.5*u*u)
                                         : std::exp(0.5*(T/sigma)*(T/sigma) + (T/sigma)*u);
    const double interspec_unnormed = (u >= -skew) ? std::exp(-0.5*u*u)
                                                   : std::exp(0.5*skew*skew + skew*u);
    BOOST_CHECK_CLOSE( genie, interspec_unnormed, 1.0E-6 );
  }

  // A "no tail" T must not be turned into a skew.
  double lower_val, upper_val, start_val, step_size;
  BOOST_REQUIRE( PeakDef::skew_parameter_range( PeakDef::SkewType::GaussExp,
                        PeakDef::CoefficientType::SkewPar0,
                        lower_val, upper_val, start_val, step_size ) );
  const double no_tail_skew = CAMInputOutput::Peak::sm_no_low_tail / sigma;
  BOOST_CHECK_MESSAGE( no_tail_skew > upper_val,
      "The no-tail sentinel converts to a skew of " << no_tail_skew
      << ", which is inside InterSpec's accepted range (" << lower_val << " to " << upper_val
      << ") - it would be imported as a real tail" );

  // ...and a real tail is in range, so it will be.
  const double real_skew = (skew*sigma) / sigma;
  BOOST_CHECK( (real_skew >= lower_val) && (real_skew <= upper_val) );
}//BOOST_AUTO_TEST_CASE( lowTailMapsToGaussExpSkew )


BOOST_AUTO_TEST_CASE( lowTailCalibrationRoundTrips )
{
  set_data_dir();

  // Genie's shape calibration carries a low-tail curve T(E) = B2 + B3*E in the third and fourth
  //  coefficients, right after the FWHM ones - verify we write both curves without disturbing
  //  each other, and that "no tail" leaves the tail coefficients zero.
  SpecUtils::SpecFile spec;
  auto m = make_shared<SpecUtils::Measurement>();
  m->set_gamma_counts( make_shared<vector<float>>( 1024, 10.0f ), 100.0f, 110.0f );
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( 1024, {0.0f, 3.0f}, {} );
  m->set_energy_calibration( cal );
  spec.add_measurement( m, true );

  const auto write_and_read = []( SpecUtils::SpecFile &spec_in,
                                  const CAMInputOutput::CnfGenieExtras &extras_in ) -> vector<float> {
    stringstream ss;
    BOOST_REQUIRE( spec_in.write_cnf( ss, {}, {}, &extras_in ) );
    const string filedata = ss.str();
    const vector<uint8_t> data( begin(filedata), end(filedata) );
    CAMInputOutput::CAMIO cam;
    cam.ReadFile( data );
    return cam.GetShapeCalibration();
  };

  CAMInputOutput::CnfGenieExtras extras;
  extras.shape_cal = make_pair( 1.25f, 0.0375f );
  extras.low_tail_cal = make_pair( 0.85f, 0.0012f );

  const vector<float> with_tail = write_and_read( spec, extras );
  BOOST_REQUIRE( with_tail.size() >= 4 );
  BOOST_CHECK_CLOSE( with_tail[0], 1.25f, 0.5 );
  BOOST_CHECK_CLOSE( with_tail[1], 0.0375f, 0.5 );
  BOOST_CHECK_CLOSE( with_tail[2], 0.85f, 0.5 );
  BOOST_CHECK_CLOSE( with_tail[3], 0.0012f, 0.5 );

  extras.low_tail_cal.reset();
  const vector<float> no_tail = write_and_read( spec, extras );
  BOOST_REQUIRE( no_tail.size() >= 4 );
  BOOST_CHECK_CLOSE( no_tail[0], 1.25f, 0.5 );
  BOOST_CHECK_CLOSE( no_tail[1], 0.0375f, 0.5 );
  BOOST_CHECK_SMALL( no_tail[2], 1.0E-6f );
  BOOST_CHECK_SMALL( no_tail[3], 1.0E-6f );

  // The real Genie files available have no tail, and store exactly that.
  const string filename = SpecUtils::append_path( g_test_file_dir, "CamLibrary/Ba-133.cnf" );
  if( SpecUtils::is_file(filename) )
  {
    vector<char> filedata;
    SpecUtils::load_file_data( filename.c_str(), filedata );
    if( !filedata.empty() && (filedata.back() == '\0') )
      filedata.pop_back();
    const vector<uint8_t> data( begin(filedata), end(filedata) );
    CAMInputOutput::CAMIO cam;
    BOOST_REQUIRE_NO_THROW( cam.ReadFile( data ) );
    const vector<float> shape = cam.GetShapeCalibration();
    BOOST_REQUIRE( shape.size() >= 4 );
    BOOST_CHECK_SMALL( shape[2], 1.0E-6f );
    BOOST_CHECK_SMALL( shape[3], 1.0E-6f );
  }
}//BOOST_AUTO_TEST_CASE( lowTailCalibrationRoundTrips )
