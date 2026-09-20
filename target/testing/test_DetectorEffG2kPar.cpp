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

/*
 =============================================================================
  Validation test for the detector-characterization .PAR + DETECTOR.txt reader
 =============================================================================

 This test is deliberately NOT run by default: the files it validates against
 (binary .PAR grids, DETECTOR.txt geometry records, and the reference-tool
 .ecc/.gis efficiency outputs) are user-supplied and are NOT committed to the
 repository (vendor IP).  The executable is always BUILT (so the reader stays
 compile-checked), but registers a ctest only when INTERSPEC_ENABLE_G2K_PAR_TEST
 is ON, and even then graceful-skips every case if `--pardir=` is not given.

 Point it at a directory laid out one subdirectory per detector:

   <pardir-root>/
     <detector-A>/
       <something>.PAR            (required)
       DETECTOR.txt               (required; single- or multi-record)
       <class>/ *.ecc *.gis       (optional reference runs, grouped or flat)
     <detector-B>/
       ....PAR   DETECTOR.txt   *.ecc *.gis

 What is asserted, per discovered detector:
   1. the .PAR + DETECTOR.txt parse, and a DRF assembles from them;
   2. every eta_fep / near_field node round-trips: the assembled CeeLo response
      reproduces the grid evaluator at the node to <1e-6 relative (the algebraic
      node identity - independent of geometry precision);
   3. (optional) if python3 + the prototype decoder are reachable, the C++ grid
      evaluator matches `par_decode.py --query` at energy nodes to ~1e-6;
   4. for each reference .ecc run whose .gis geometry we can map to a point/angle
      (point, SPHERE off-axis, disk), the evaluator AND the assembled DRF match
      the .ecc within a per-geometry tolerance (tight for vacuum points/off-axis,
      looser for air and volumetric; Marinelli/theta>90 is informational only);
   5. the DRF serializes and round-trips (toXml -> fromXml), preserving the
      efficiency and the DrfSource provenance.
 =============================================================================
*/

// Must be defined before Windows.h (or any header that includes it) is included
#ifdef _WIN32
  #include <winsock2.h>
  #include <windows.h>
#endif

#define BOOST_TEST_MODULE test_DetectorEffG2kPar_suite
#include <boost/test/included/unit_test.hpp>

#include <cmath>
#include <array>
#include <string>
#include <vector>
#include <memory>
#include <limits>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <utility>
#include <iostream>
#include <algorithm>

#include <rapidxml/rapidxml.hpp>

#include <Eigen/Core>

#include "io/ResponseKernel.h"
#include "io/DetectorResponse.h"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectorEffG2kPar.h"
#include "InterSpec/CascadeSummingCalc.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;

namespace
{
  string g_data_dir;       ///< --datadir=   (nuclear data; for InterSpec statics)
  string g_test_data_dir;  ///< --testfiledir=
  string g_par_dir;        ///< --pardir=    (root of the uncommitted corpus)

  const double pi = 3.14159265358979323846;


  //===========================================================================
  //  A tiny .gis reader (source geometry definition emitted by the reference tool)
  //===========================================================================
  //  We only need: the template class (~Geometry), the radial range from the
  //  endcap face (~sd1, mm), the lateral offset that sets the polar angle for
  //  off-axis SPHERE (~sd4, mm), and whether air is present (~APressure > 0).  Format: "~key=value", one or more per line.
  struct GisGeometry
  {
    string templ;             ///< e.g. "SPHERE", "CIRCULAR_PLANE", "WELL_or_MARINELLI_BEAKER"
    double sd1_mm = numeric_limits<double>::quiet_NaN();
    double sd4_mm = numeric_limits<double>::quiet_NaN();
    double pressure = 0.0;    ///< ~APressure (0 => vacuum run)
    double disk_diam_mm = numeric_limits<double>::quiet_NaN();  ///< ~d1.2 (disk source diameter)
    bool parsed = false;
  };//struct GisGeometry


  GisGeometry parseGis( const string &path )
  {
    GisGeometry g;
    ifstream in( path.c_str(), ios::in | ios::binary );
    if( !in.is_open() )
      return g;

    string line;
    while( std::getline( in, line ) )
    {
      // Tokens are whitespace-separated "~key=value" pairs; comments start '#'.
      string trimmed = line;
      SpecUtils::trim( trimmed );
      if( trimmed.empty() || trimmed[0] == '#' )
        continue;

      vector<string> toks;
      SpecUtils::split( toks, trimmed, " \t" );
      for( const string &tok : toks )
      {
        const size_t eq = tok.find( '=' );
        if( (tok.empty()) || (tok[0] != '~') || (eq == string::npos) )
          continue;
        const string key = tok.substr( 1, eq - 1 );
        const string val = tok.substr( eq + 1 );
        if( val.empty() )
          continue;

        if( SpecUtils::iequals_ascii( key, "Geometry" ) )
          g.templ = val;
        else if( SpecUtils::iequals_ascii( key, "sd1" ) )
          SpecUtils::parse_double( val.c_str(), val.size(), g.sd1_mm );
        else if( SpecUtils::iequals_ascii( key, "sd4" ) )
          SpecUtils::parse_double( val.c_str(), val.size(), g.sd4_mm );
        else if( SpecUtils::iequals_ascii( key, "APressure" ) )
          SpecUtils::parse_double( val.c_str(), val.size(), g.pressure );
        else if( SpecUtils::iequals_ascii( key, "d1.2" ) )
          SpecUtils::parse_double( val.c_str(), val.size(), g.disk_diam_mm );
      }//for( each token )
    }//while( getline )

    g.parsed = !g.templ.empty() && !std::isnan( g.sd1_mm );
    return g;
  }//parseGis(...)


  //===========================================================================
  //  Corpus discovery
  //===========================================================================
  struct EccRun
  {
    string eccPath;
    string gisPath;   ///< may be empty if no like-named .gis sits beside it
  };//struct EccRun

  struct DetectorCase
  {
    string dir;
    string parPath;
    string detectorTxtPath;
    vector<EccRun> runs;
  };//struct DetectorCase


  /** Finds the .par and DETECTOR.txt at the top of a detector directory, then
   recurses within it for .ecc runs (each paired to its like-named .gis). */
  bool buildDetectorCase( const string &detDir, DetectorCase &out )
  {
    out.dir = detDir;

    // .par / .PAR at the top of the detector directory.
    const vector<string> top = SpecUtils::ls_files_in_directory( detDir );
    for( const string &f : top )
    {
      if( SpecUtils::iends_with( f, ".par" ) )
      {
        out.parPath = f;
        break;
      }
    }
    if( out.parPath.empty() )
      return false;

    // DETECTOR.txt (or a lone detector.txt) at the top.
    for( const string &f : top )
    {
      const string base = SpecUtils::to_lower_ascii_copy( SpecUtils::filename( f ) );
      if( base == "detector.txt" )
      {
        out.detectorTxtPath = f;
        break;
      }
    }
    // Fallback: any *.txt at the top mentioning a ".par" token would also work,
    //  but keep it strict - a detector dir is expected to name it DETECTOR.txt.
    if( out.detectorTxtPath.empty() )
      return false;

    // Reference runs: recurse for .ecc, pair with a like-named .gis in the same
    //  directory.  Matching is case-insensitive on the stem, because the sample
    //  corpus has case-mismatched pairs (e.g. Lab06_*.ecc vs lab06_*.gis) that a
    //  literal path build would miss on a case-sensitive filesystem.
    auto stem_lower = []( const string &path ) -> string
    {
      string base = SpecUtils::filename( path );
      const size_t dot = base.find_last_of( '.' );
      if( dot != string::npos )
        base = base.substr( 0, dot );
      return SpecUtils::to_lower_ascii_copy( base );
    };

    const vector<string> giss = SpecUtils::recursive_ls( detDir, ".gis" );
    const vector<string> eccs = SpecUtils::recursive_ls( detDir, ".ecc" );
    for( const string &ecc : eccs )
    {
      EccRun run;
      run.eccPath = ecc;
      const string ecc_dir = SpecUtils::parent_path( ecc );
      const string ecc_stem = stem_lower( ecc );
      for( const string &gis : giss )
      {
        if( SpecUtils::parent_path( gis ) == ecc_dir && stem_lower( gis ) == ecc_stem )
        {
          run.gisPath = gis;
          break;
        }
      }
      out.runs.push_back( run );
    }
    return true;
  }//buildDetectorCase(...)


  vector<DetectorCase> discoverCases()
  {
    vector<DetectorCase> cases;
    if( g_par_dir.empty() || !SpecUtils::is_directory( g_par_dir ) )
      return cases;

    // Canonical layout: one subdirectory per detector.
    const vector<string> subdirs = SpecUtils::ls_directories_in_directory( g_par_dir );
    for( const string &sub : subdirs )
    {
      DetectorCase c;
      if( buildDetectorCase( sub, c ) )
        cases.push_back( std::move(c) );
    }

    // Fallback for "drop files anywhere": if no per-detector subdir matched but
    //  the root itself holds a .par + DETECTOR.txt, treat the root as one case.
    if( cases.empty() )
    {
      DetectorCase c;
      if( buildDetectorCase( g_par_dir, c ) )
        cases.push_back( std::move(c) );
    }

    return cases;
  }//discoverCases(...)


  //===========================================================================
  //  Helpers
  //===========================================================================

  /** Reads .ecc energy/efficiency pairs via the production parser. */
  vector<pair<double,double>> readEcc( const string &path )
  {
    vector<pair<double,double>> pairs;
    ifstream in( path.c_str(), ios::in | ios::binary );
    if( !in.is_open() )
      return pairs;

    try
    {
      const DetectorPeakResponse::EccParseResult res
                                   = DetectorPeakResponse::parseEccFile( in );
      if( !res.drf )
        return pairs;
      // Sample the parsed legacy curve back at its node energies via the DRF.
      //  (parseEccFile stored the pairs into an efficiency curve; re-read them
      //  at the uncertainty node energies, which are the .ecc energies.)
      for( const float e : res.uncertEnergies )
      {
        const DetectorPeakResponse::EffEval ev = res.drf->intrinsicEfficiencyEval( e );
        pairs.emplace_back( static_cast<double>(e), ev.value );
      }
    }catch( const std::exception &e )
    {
      BOOST_TEST_MESSAGE( string("readEcc: parseEccFile threw: ") + e.what() );
    }
    return pairs;
  }//readEcc(...)


  /** Optional grid-vs-Python parity: shell out to the prototype decoder if it
   and python3 are reachable.  Returns false (skip) if not runnable. */
  bool pythonQuery( const string &parDecodePy, const string &parPath,
                    double e_keV, double d_mm, double theta_deg, double &eff_out )
  {
    char cmd[2048];
    std::snprintf( cmd, sizeof(cmd),
                   "python3 '%s' --query '%s' %.6f %.6f %.6f 2>/dev/null",
                   parDecodePy.c_str(), parPath.c_str(), e_keV, d_mm, theta_deg );
    FILE *pipe = popen( cmd, "r" );
    if( !pipe )
      return false;

    string output;
    char buf[512];
    while( fgets( buf, sizeof(buf), pipe ) )
      output += buf;
    const int rc = pclose( pipe );
    if( rc != 0 )
      return false;

    // Parse "eff_vac=1.234567e-04" out of the output.
    const size_t pos = output.find( "eff_vac=" );
    if( pos == string::npos )
      return false;
    const string num = output.substr( pos + 8 );
    return SpecUtils::parse_double( num.c_str(), num.size(), eff_out );
  }//pythonQuery(...)


  struct TestFixture
  {
    TestFixture()
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
        else if( arg.find("--pardir=") == 0 )
          g_par_dir = arg.substr( 9 );
      }

      if( !g_data_dir.empty() )
      {
        try
        {
          InterSpec::setStaticDataDirectory( g_data_dir );
        }catch( std::exception &e )
        {
          cerr << "Warning: failed to set static data dir: " << e.what() << endl;
        }
      }
    }
  };//struct TestFixture
}//anonymous namespace

BOOST_GLOBAL_FIXTURE( TestFixture );


//=============================================================================
//  Tests
//=============================================================================

BOOST_AUTO_TEST_CASE( DiscoveryAndParse )
{
  if( g_par_dir.empty() )
  {
    BOOST_TEST_MESSAGE( "test_DetectorEffG2kPar: skipped (pass --pardir=<root> to run; "
                        "the .PAR/DETECTOR.txt corpus is not committed)." );
    return;
  }

  const vector<DetectorCase> cases = discoverCases();
  BOOST_REQUIRE_MESSAGE( !cases.empty(),
      "No detector cases discovered under --pardir=" + g_par_dir
      + " (expected <root>/<detector>/{*.PAR, DETECTOR.txt})." );

  for( const DetectorCase &c : cases )
  {
    BOOST_TEST_MESSAGE( "Detector case: " + c.dir + "  (" + std::to_string(c.runs.size())
                        + " reference runs)" );

    // Parse both files and assemble a DRF.
    shared_ptr<DetectorPeakResponse> drf;
    BOOST_REQUIRE_NO_THROW( drf = DetEffG2kPar::makeDrfFromFiles( c.parPath, c.detectorTxtPath ) );
    BOOST_REQUIRE( !!drf );
    BOOST_CHECK( drf->isValid() );
    BOOST_CHECK( drf->drfSource() == DetectorPeakResponse::DrfSource::CharacterizationParFile );
    BOOST_CHECK( !!drf->ceeloResponse() );

    // These files characterize the FULL-ENERGY PEAK only, so nothing downstream
    //  may believe a total efficiency is available.  The bare crystal kernel is
    //  NOT an acceptable stand-in: it omits the housing and every peak-to-total
    //  effect, and at low energy it comes out BELOW the response's own eps_fep,
    //  which is impossible - a host gating cascade summing on "is there a
    //  total?" would get a confident yes and a wrong correction.
    BOOST_CHECK( !drf->ceeloResponse()->tot_eff.characterized() );
    BOOST_CHECK( !drf->hasAnyTotalEfficiencyInfo() );
    BOOST_CHECK( !GammaInteractionCalc::CascadeSummingCalc::drfHasNeededInfo( drf ) );
    for( const float E : { 60.0f, 300.0f, 1332.0f } )
    {
      BOOST_CHECK_EQUAL( drf->totalIntrinsicEfficiencyAny( E ), 0.0f );
      BOOST_CHECK( drf->farFieldIntrinsicEfficiency( E ) > 0.0f );   //FEP is unaffected
    }
    BOOST_CHECK_EQUAL(
      GammaInteractionCalc::CascadeSummingCalc::estimateMaxSummingMagnitude(
                                          drf, 25.0*PhysicalUnits::cm ), 0.0 );
  }//for( each case )
}//BOOST_AUTO_TEST_CASE( DiscoveryAndParse )


/** The layer-material tokens are elemental symbols, and ANY element must reach
 the geometry - not just the handful a hand-written table happened to list.  A
 table like that dropped a whole 1.5 mm magnesium endcap on one detector in the
 sample corpus, silently: the descriptor came out with the vacuum gap and the
 copper liner and no endcap at all.

 Needs no vendor files: a DETECTOR.txt record is plain ASCII, and a synthetic
 grid drives makeDrf.  `--datadir=` is required, since the element lookup goes
 through SandiaDecay.
 */
BOOST_AUTO_TEST_CASE( LayerMaterialsResolveAnyElement )
{
  if( g_data_dir.empty() )
  {
    BOOST_TEST_MESSAGE( "LayerMaterialsResolveAnyElement: skipped (needs --datadir=)." );
    return;
  }

  // Slots are front/side/back triplets: 0-2 front (Ge dead layer, window,
  //  endcap), 3-5 side, 6-8 back.  Mg/Be/Ti/Zn are all outside the shielding
  //  elements CeeLo ships a builtin for, so each one exercises the SandiaDecay
  //  path with the file's own density.
  const string detector_txt =
    "# 12345 - SYNTHETIC - S/N-TEST-1\n"
    "SynthDet,70.0,60.0,0,80.0,150.0,4.0,4.0,26,synth.par,4, #\n"
    "ge,0.7,5.35, #\n"
    "be,0.5,1.848, #\n"
    "mg,1.5,1.74, #\n"
    "ge,0.7,5.35, #\n"
    "ti,1.0,4.507, #\n"
    "zn,2.0,7.14, #\n"
    "ge,0.5,5.35, #\n"
    "cu,3.0,8.96, #\n"
    "mg,5.0,1.74\n";

  istringstream txt( detector_txt );
  const vector<DetEffG2kPar::DetectorDef> defs = DetEffG2kPar::parseDetectorTxt( txt );
  BOOST_REQUIRE_EQUAL( defs.size(), 1 );
  const DetEffG2kPar::DetectorDef &def = defs.front();
  BOOST_CHECK_EQUAL( def.layers[1].material, "be" );
  BOOST_CHECK_EQUAL( def.layers[2].material, "mg" );
  BOOST_CHECK_EQUAL( def.layers[4].material, "ti" );
  BOOST_CHECK_EQUAL( def.layers[5].material, "zn" );

  // A small but well-formed grid: two energies, monotone in distance/angle.
  DetEffG2kPar::ParFile par;
  par.emin_keV = 60.0;
  par.emax_keV = 1332.0;
  par.energies_keV = { 60.0, 1332.0 };
  for( size_t e = 0; e < par.energies_keV.size(); ++e )
  {
    DetEffG2kPar::ParGrid g;
    g.ncols = 19;                    // 0..180 deg in 10 deg steps
    g.nrows = 40;
    g.theta_step_rad = 10.0 * 3.14159265358979323846 / 180.0;
    g.r_step = 0.2;                  // ln(mm)
    g.V.resize( static_cast<size_t>(g.nrows) * g.ncols );
    for( int r = 0; r < g.nrows; ++r )
    {
      for( int c = 0; c < g.ncols; ++c )
      {
        // eff falls with distance and angle, and with energy: V rises.
        const double V = 2000.0 + 900.0*e + 40.0*r + 15.0*c;
        g.V[r*g.ncols + c] = static_cast<uint16_t>( V );
      }
    }
    par.grids.push_back( g );
  }

  shared_ptr<DetectorPeakResponse> drf;
  BOOST_REQUIRE_NO_THROW( drf = DetEffG2kPar::makeDrf( par, def ) );
  BOOST_REQUIRE( !!drf );
  BOOST_REQUIRE( !!drf->ceeloResponse() );

  // Every modeled token must appear in the descriptor's material table.  Only
  //  slots 1/2 (front) and 4/5 (side) become CeeLo layers - the back path is
  //  not modeled - so Mg, Be, Ti and Zn are the ones to look for.
  const ceelo::GeometryDescriptor &gd = drf->ceeloResponse()->descriptor;
  auto has_material = [&gd]( const string &name ) -> bool
  {
    for( const ceelo::MaterialSpec &m : gd.materials )
    {
      if( m.name == name )
        return true;
    }
    return false;
  };
  for( const char * const name : { "Be", "Mg", "Ti", "Zn" } )
    BOOST_CHECK_MESSAGE( has_material( name ),
                         string("layer material ") + name + " never reached the descriptor" );

  // And the germanium crystal must NOT be duplicated by the "ge" dead-layer
  //  tokens: the material table dedupes BY NAME, so a token resolved to a
  //  differently-named germanium would silently add a second one.
  size_t n_ge = 0;
  for( const ceelo::MaterialSpec &m : gd.materials )
  {
    if( !m.composition.empty() && (m.composition.size() == 1)
        && (m.composition.front().Z == 32) )
      ++n_ge;
  }
  BOOST_CHECK_EQUAL( n_ge, 1 );

  // No layer was skipped, so the description carries no omission note.
  BOOST_CHECK( drf->description().find( "not modeled" ) == string::npos );

  // A token that is NOT an element cannot be modeled: the layer is dropped, but
  //  the DRF still builds and says so, rather than losing it silently.
  {
    DetEffG2kPar::DetectorDef bad = def;
    bad.layers[2].material = "unobtainium";
    shared_ptr<DetectorPeakResponse> bad_drf;
    BOOST_REQUIRE_NO_THROW( bad_drf = DetEffG2kPar::makeDrf( par, bad ) );
    BOOST_REQUIRE( !!bad_drf );
    BOOST_CHECK( bad_drf->description().find( "not modeled" ) != string::npos );
    BOOST_CHECK( bad_drf->description().find( "unobtainium" ) != string::npos );
  }
}//BOOST_AUTO_TEST_CASE( LayerMaterialsResolveAnyElement )


BOOST_AUTO_TEST_CASE( NodeRoundTripIdentity )
{
  if( g_par_dir.empty() )
    return;

  const vector<DetectorCase> cases = discoverCases();
  for( const DetectorCase &c : cases )
  {
    DetEffG2kPar::ParFile par;
    BOOST_REQUIRE_NO_THROW( par = DetEffG2kPar::parseParFile( c.parPath ) );

    ifstream txt( c.detectorTxtPath.c_str(), ios::in | ios::binary );
    BOOST_REQUIRE( txt.is_open() );
    const vector<DetEffG2kPar::DetectorDef> defs = DetEffG2kPar::parseDetectorTxt( txt );
    const DetEffG2kPar::DetectorDef def
                = DetEffG2kPar::selectDetectorDef( defs, SpecUtils::filename( c.parPath ) );

    shared_ptr<DetectorPeakResponse> drf;
    BOOST_REQUIRE_NO_THROW( drf = DetEffG2kPar::makeDrf( par, def ) );
    const shared_ptr<const ceelo::DetectorResponse> resp = drf->ceeloResponse();
    BOOST_REQUIRE( !!resp );

    const DetEffG2kPar::ParEfficiency parEff( par );

    // For every eta_fep node, and every near_field distance node, the response
    //  must reproduce the grid value at that node to within a hair - this is
    //  the algebraic node identity, independent of any geometry approximation.
    const ceelo::EtaTable &eta = resp->eta_fep;
    const Eigen::Vector3d face = CeeLoUtils::detectorFacePosition( resp->descriptor );

    size_t nchecked = 0, nfailed = 0, nskipped = 0, nbadsigma = 0;
    double worst_rel = 0.0;

    auto check_node = [&]( double E, const Eigen::Vector3d &src )
    {
      // Mirror makeDrf's face_sample: the grid row is the RADIAL range |v| from
      //  the face centre, not the axial height (-v.z()) - the spherical-grid
      //  convention; see the DetectorEffG2kPar.h header.
      const Eigen::Vector3d v = src - face;
      const double n = v.norm();
      const double d_face_mm = n * 10.0;
      const double theta_face = (n > 0.0)
                    ? std::acos( std::max(-1.0, std::min(1.0, -v.z()/n)) ) : 0.0;

      // Rungs whose face-frame position lands in the grid's no-data region
      //  (V == 0, which would decode to a 100% efficiency) carry no file value to
      //  be identical TO: makeDrf truncates each angular column's ladder there.
      //  Exempt them from the exactness check, but pin the truncation contract -
      //  they must be marked fully uncertain, so a query cannot land on one and
      //  come back looking authoritative.
      bool no_data = false;
      parEff.efficiency( E, d_face_mm, theta_face, &no_data );
      if( no_data )
      {
        ++nskipped;
        const double cos_theta = (n > 0.0) ? (-src.z() / src.norm()) : 1.0;
        const double frac_sigma = resp->near_field.node_frac_sigma(
                                            E, cos_theta, src.norm() );
        if( !(frac_sigma >= 0.99) )
          ++nbadsigma;
        return;
      }

      const double expect = parEff.efficiency( E, d_face_mm, theta_face );
      const ceelo::EffResult got = resp->eps_fep_at( E, src );
      const double rel = std::fabs( got.value - expect ) / std::max( expect, 1e-300 );
      worst_rel = std::max( worst_rel, rel );
      ++nchecked;
      if( rel > 1e-6 )
        ++nfailed;
    };

    // eta_fep nodes (far reference distance).
    const double d_ref_cm = resp->near_field.dists_cm.empty()
                              ? 300.0 : resp->near_field.dists_cm.back();
    for( size_t ci = 0; ci < eta.cos_thetas.size(); ++ci )
    {
      const Eigen::Vector3d src = ceelo::source_position( d_ref_cm, eta.cos_thetas[ci], 0.0 );
      for( const double E : eta.energies_keV )
        check_node( E, src );
    }

    // near_field distance nodes.
    const ceelo::NearFieldModel &nf = resp->near_field;
    for( size_t ci = 0; ci < nf.cos_thetas.size(); ++ci )
    {
      for( const double d_cm : nf.dists_cm )
      {
        const Eigen::Vector3d src = ceelo::source_position( d_cm, nf.cos_thetas[ci], 0.0 );
        for( const double E : nf.energies_keV )
          check_node( E, src );
      }
    }

    BOOST_TEST_MESSAGE( "  node round-trip: " + std::to_string(nchecked) + " nodes, worst rel="
                        + std::to_string(worst_rel) + ", " + std::to_string(nskipped)
                        + " truncated rungs skipped" );
    BOOST_CHECK_MESSAGE( nfailed == 0,
        std::to_string(nfailed) + "/" + std::to_string(nchecked)
        + " nodes exceeded 1e-6 relative (worst " + std::to_string(worst_rel) + ")" );

    // Truncated rungs must be flagged as fully uncertain, not served silently.
    BOOST_CHECK_MESSAGE( nbadsigma == 0,
        std::to_string(nbadsigma) + "/" + std::to_string(nskipped)
        + " truncated rungs lacked the expected ~100% fractional sigma" );

    // The exemption must stay narrow: the no-data region is behind the endcap at
    //  close range only, so the vast majority of rungs are still checked exactly.
    BOOST_CHECK_MESSAGE( nskipped * 4 < nchecked,
        "too many rungs exempted (" + std::to_string(nskipped) + " skipped vs "
        + std::to_string(nchecked) + " checked) - the no-data region should be a"
        " small grazing-angle corner, so this suggests a frame or geometry error" );
  }//for( each case )
}//BOOST_AUTO_TEST_CASE( NodeRoundTripIdentity )


BOOST_AUTO_TEST_CASE( GridVsPython )
{
  if( g_par_dir.empty() )
    return;

  // Locate the prototype decoder, which lives in the uncommitted scratch/ tree.
  //  Checked relative to the CWD and to --datadir first so this is not tied to
  //  one machine; PAR_DECODE_PY overrides.  If it (or python3) is not reachable
  //  this test just skips - it is a cross-check, not a gate.
  string parDecodePy;
  {
    const char * const env = getenv( "PAR_DECODE_PY" );
    const string rel = "scratch/genie_par_dev/tools/par_decode.py";
    vector<string> candidates;
    if( env && env[0] )
      candidates.push_back( env );
    candidates.push_back( rel );
    candidates.push_back( "../" + rel );
    candidates.push_back( "../../" + rel );
    if( !g_data_dir.empty() )
      candidates.push_back( SpecUtils::append_path( SpecUtils::parent_path(g_data_dir), rel ) );

    for( const string &c : candidates )
    {
      if( SpecUtils::is_file( c ) )
      {
        parDecodePy = c;
        break;
      }
    }//for( const string &c : candidates )
  }

  if( parDecodePy.empty() )
  {
    BOOST_TEST_MESSAGE( "GridVsPython: skipped (prototype decoder not present)." );
    return;
  }

  const vector<DetectorCase> cases = discoverCases();
  for( const DetectorCase &c : cases )
  {
    DetEffG2kPar::ParFile par;
    BOOST_REQUIRE_NO_THROW( par = DetEffG2kPar::parseParFile( c.parPath ) );
    const DetEffG2kPar::ParEfficiency parEff( par );

    // Probe a spread of (energy node, distance, angle) points; compare to the
    //  Python at energy NODES (where PCHIP and the Python's linear-in-lnE agree).
    const array<double,3> dists_mm = { 25.0, 150.0, 1000.0 };
    const array<double,3> thetas_deg = { 0.0, 20.0, 60.0 };

    size_t nchecked = 0, nfailed = 0;
    bool python_ran = false;
    for( const double E : par.energies_keV )
    {
      for( const double d_mm : dists_mm )
      {
        for( const double th_deg : thetas_deg )
        {
          double py_eff = 0.0;
          if( !pythonQuery( parDecodePy, c.parPath, E, d_mm, th_deg, py_eff ) )
            continue;
          python_ran = true;
          const double cpp_eff = parEff.efficiency( E, d_mm, th_deg * pi / 180.0 );
          const double rel = std::fabs( cpp_eff - py_eff ) / std::max( py_eff, 1e-300 );
          ++nchecked;
          if( rel > 1e-6 )
          {
            ++nfailed;
            BOOST_TEST_MESSAGE( "  mismatch E=" + std::to_string(E) + " d=" + std::to_string(d_mm)
                                + "mm th=" + std::to_string(th_deg) + "deg  cpp=" + std::to_string(cpp_eff)
                                + " py=" + std::to_string(py_eff) + " rel=" + std::to_string(rel) );
          }
        }
      }
    }//for( energies )

    if( !python_ran )
    {
      BOOST_TEST_MESSAGE( "GridVsPython: python3 not runnable; skipped." );
      return;
    }
    BOOST_CHECK_MESSAGE( nfailed == 0, std::to_string(nfailed) + "/" + std::to_string(nchecked)
                         + " grid-vs-Python node queries exceeded 1e-6 for " + c.parPath );
  }//for( each case )
}//BOOST_AUTO_TEST_CASE( GridVsPython )


BOOST_AUTO_TEST_CASE( EccMatch )
{
  if( g_par_dir.empty() )
    return;

  const vector<DetectorCase> cases = discoverCases();
  for( const DetectorCase &c : cases )
  {
    DetEffG2kPar::ParFile par;
    BOOST_REQUIRE_NO_THROW( par = DetEffG2kPar::parseParFile( c.parPath ) );

    ifstream txt( c.detectorTxtPath.c_str(), ios::in | ios::binary );
    BOOST_REQUIRE( txt.is_open() );
    const vector<DetEffG2kPar::DetectorDef> defs = DetEffG2kPar::parseDetectorTxt( txt );
    const DetEffG2kPar::DetectorDef def
                = DetEffG2kPar::selectDetectorDef( defs, SpecUtils::filename( c.parPath ) );

    const DetEffG2kPar::ParEfficiency parEff( par );
    shared_ptr<DetectorPeakResponse> drf;
    BOOST_REQUIRE_NO_THROW( drf = DetEffG2kPar::makeDrf( par, def ) );

    for( const EccRun &run : c.runs )
    {
      if( run.gisPath.empty() )
        continue;

      const GisGeometry gis = parseGis( run.gisPath );
      if( !gis.parsed )
      {
        BOOST_TEST_MESSAGE( "  skip (unparsed .gis): " + run.gisPath );
        continue;
      }

      const vector<pair<double,double>> ecc = readEcc( run.eccPath );
      if( ecc.empty() )
        continue;

      const string &templ = gis.templ;  // icontains is case-insensitive
      const bool is_point = SpecUtils::icontains( templ, "SPHERE" )
                          || SpecUtils::icontains( templ, "POINT" )
                          || SpecUtils::icontains( templ, "CIRCULAR_PLANE" );
      const bool is_marinelli = SpecUtils::icontains( templ, "MARINELLI" )
                             || SpecUtils::icontains( templ, "WELL" );

      // Geometry -> (distance-from-face mm, theta from axis).  The PAR grid is
      //  SPHERICAL: its distance axis is the RADIAL range from the endcap-face
      //  centre = ~sd1, and ~sd4 only sets the direction, theta = atan2(sd4, sd1).
      //  The off-axis runs say so themselves (~Geometry=SPHERE); see the
      //  DetectorEffG2kPar.h header for the full evidence.  With this reading
      //  every off-axis reference run reproduces to <0.08% at node energies
      //  (worst over the corpus: 0.0785%, LAB06 at sd1=200/sd4=350 mm); the
      //  axial reading misses by tens of percent to 84x.
      const double sd1 = gis.sd1_mm;
      const double sd4 = std::isnan( gis.sd4_mm ) ? 0.0 : gis.sd4_mm;
      const double dist_mm = sd1;                     // radial range = grid row
      const double theta = std::atan2( sd4, sd1 );    // polar angle from the axis

      // A disk source (CIRCULAR_PLANE with a large ~d1.2 diameter) is extended;
      //  point/off-axis are the clean cases.  Marinelli reaches theta>90.
      const bool is_disk = SpecUtils::icontains( templ, "CIRCULAR_PLANE" )
                        && !std::isnan( gis.disk_diam_mm ) && (gis.disk_diam_mm > 5.0);

      // Tolerances differ fundamentally between ENERGY-NODE and BETWEEN-NODE
      //  points.  The .PAR grid stores efficiency only at ~15-20 energy nodes; at
      //  a node the reader reproduces the reference tool very closely: over the
      //  corpus's vacuum point/off-axis runs out to theta=84deg, both detectors,
      //  14 of 16 runs are within 0.003% and the two worst are 0.044% and
      //  0.079%.  (Air runs add the analytic attenuation factor and sit near
      //  0.03%, and are gated looser at 1% since the attenuation factor is our
      //  own analytic model rather than a grid lookup.)  The vacuum gate is
      //  0.1%, which those two worst runs genuinely need.
      //  BETWEEN nodes both tools interpolate the SAME sparse grid, but with
      //  different schemes (the reference behaves like a local forward-quadratic
      //  in (V, ln E); we use PCHIP) - so a ~1-2% gap is expected and is NOT our
      //  error against truth; an independent Monte-Carlo arbiter says PCHIP is the
      //  more physical of the two there.  See the DetectorEffG2kPar.h header for
      //  the measured comparison and the ~0.23% quantization floor.  We gate the
      //  two classes separately: tight at nodes (proves framing/grid/bilinear/angle
      //  mapping), loose between nodes (documents the grid's energy-resolution
      //  floor).  Air adds the exp(-mu_air*d) term, matched to NIST mu_total to
      //  ~1-3% (looser at 45 keV where coherent scatter forward-peaks).
      auto is_energy_node = [&par]( const double E ) -> bool
      {
        for( const double n : par.energies_keV )
          if( std::fabs( E - n ) <= 1e-3 * std::max(n,1.0) )
            return true;
        return false;
      };

      const double node_tol = (gis.pressure > 0.0) ? 0.01  : 0.001;
      const double btwn_tol = (gis.pressure > 0.0) ? 0.03  : 0.025;

      size_t n_node = 0, n_node_fail = 0, n_btwn = 0, n_btwn_fail = 0;
      double worst_node = 0.0, worst_btwn = 0.0;
      for( const pair<double,double> &pt : ecc )
      {
        const double E = pt.first;
        const double ecc_eff = pt.second;
        if( !(ecc_eff > 0.0) )
          continue;

        double eff = parEff.efficiency( E, dist_mm, theta );
        if( gis.pressure > 0.0 )
        {
          // The source sits at radial range sd1, so the air path is that same
          //  range - the grid row and the air path coincide here.
          eff *= parEff.airTransmission( E, dist_mm, gis.pressure / 760.0 );
        }

        const double rel = std::fabs( eff - ecc_eff ) / ecc_eff;
        if( is_energy_node( E ) )
        {
          worst_node = std::max( worst_node, rel );
          ++n_node;
          if( rel > node_tol )
            ++n_node_fail;
        }else
        {
          worst_btwn = std::max( worst_btwn, rel );
          ++n_btwn;
          if( rel > btwn_tol )
            ++n_btwn_fail;
        }
      }//for( each ecc point )

      const string tag = SpecUtils::filename( run.eccPath ) + " [" + templ
                       + (gis.pressure > 0.0 ? ", air" : ", vac") + "]";

      if( is_marinelli || is_disk || !is_point )
      {
        // Extended (disk), volumetric (Marinelli), or behind-plane source: the
        //  point model evaluated at the source centre cannot represent the
        //  reference tool's per-element integral; report only, do not gate.
        BOOST_TEST_MESSAGE( "  (informational) " + tag + "  worst node="
                            + std::to_string(worst_node) + " btwn=" + std::to_string(worst_btwn)
                            + " over " + std::to_string(n_node+n_btwn) + " pts" );
        continue;
      }

      BOOST_TEST_MESSAGE( "  " + tag + "  nodes " + std::to_string(n_node)
                          + " (worst " + std::to_string(worst_node) + "), between "
                          + std::to_string(n_btwn) + " (worst " + std::to_string(worst_btwn) + ")" );
      BOOST_CHECK_MESSAGE( n_node_fail == 0,
          tag + ": " + std::to_string(n_node_fail) + "/" + std::to_string(n_node)
          + " ENERGY-NODE points exceeded " + std::to_string(node_tol)
          + " (worst " + std::to_string(worst_node) + ")" );
      BOOST_CHECK_MESSAGE( n_btwn_fail == 0,
          tag + ": " + std::to_string(n_btwn_fail) + "/" + std::to_string(n_btwn)
          + " between-node points exceeded " + std::to_string(btwn_tol)
          + " (worst " + std::to_string(worst_btwn) + ")" );

      // Also cross-check the assembled DRF matches the grid evaluator at the SAME
      //  physical point, through the InterSpec dispatch.  fepEfficiencyEval's
      //  `distance` is the range from the face centre, which for this spherical
      //  grid IS the grid row, so we pass dist_mm (=sd1) directly.
      //
      // Do NOT pass hypot(sd1,sd4) here.  That was the earlier form, and it made
      //  this check blind to the very convention it is meant to police: because
      //  hypot(sd1,sd4)*cos(atan2(sd4,sd1)) == sd1 identically, an axial-sampling
      //  bug in makeDrf and a slant-range query here are exact inverses, so the
      //  check passed under either convention while assembled DRFs were in fact
      //  reading the wrong grid row off-axis (92% error at 45deg, 8427% at 84deg).
      //  Querying the TRUE radial range is what makes this check discriminating.
      //
      // This is a coarser check than the node round-trip: the assembled CeeLo
      //  tables subsample the grid (an ~18-node near-field distance ladder and
      //  the ~2.5 deg cos-theta nodes) and interpolate on their own manifold, so
      //  a query that is NOT exactly a fill node differs from the raw-grid
      //  bilinear by the table-vs-grid interpolation gap (<~0.6% at the mid-index
      //  energy probed below; that gap is much larger - tens of percent - at the
      //  10-16 keV low end and at extreme angles, which is why this check probes
      //  a mid energy rather than sweeping).
      //  Its job is to catch a gross frame/convention error in the dispatch path
      //  (which would show up as many-percent off-axis, as the pre-fix slant
      //  sampling did), not to re-prove the node identity.  Probe at an energy
      //  node so only the spatial dimension is between-node.
      if( gis.pressure <= 0.0 && !par.energies_keV.empty() )
      {
        const double E_node = par.energies_keV[ par.energies_keV.size()/2 ];
        const double eff_eval = parEff.efficiency( E_node, dist_mm, theta );
        const DetectorPeakResponse::EffEval ev = drf->fepEfficiencyEval(
                    static_cast<float>(E_node), theta, 0.0, (dist_mm/10.0)*PhysicalUnits::cm );
        BOOST_CHECK_MESSAGE(
            std::fabs(ev.value - eff_eval) <= 0.015 * std::max(eff_eval,1e-300) + 1e-30,
            tag + ": DRF dispatch " + std::to_string(ev.value) + " != evaluator "
            + std::to_string(eff_eval) + " at node " + std::to_string(E_node) + " keV" );
      }
    }//for( each ecc run )
  }//for( each case )
}//BOOST_AUTO_TEST_CASE( EccMatch )


/** Validates the AIR TERM IN ISOLATION, which the absolute-efficiency checks in
 EccMatch cannot do at their 1%/3% tolerances.

 The corpus pairs each air run with a vacuum run at the SAME position, so the
 ratio `ecc_air / ecc_vac` is the reference tool's own air factor with the grid,
 the spatial interpolation and the energy interpolation all divided out.  That
 makes this a ~0.03% probe of `airTransmission` alone, over five standoffs
 (15..999.9 mm) and 45-2000 keV, and it is what pins the removal coefficient to
 mu_total - mu_Rayleigh: with Rayleigh wrongly included the worst error here is
 0.296% and GROWS with path length (the signature of a wrong mu), versus 0.028%
 with it excluded.  Point-like sources only - an extended source needs a
 per-element path integral (see the airTransmission docs).
 */
BOOST_AUTO_TEST_CASE( AirTransmissionOnly )
{
  if( g_par_dir.empty() )
    return;

  // Slack over the 0.028% observed worst case: enough that print precision and
  //  a differing air composition/density do not make this flaky, tight enough
  //  that re-including Rayleigh (0.296%) fails loudly.
  const double air_tol = 0.0005;   // 0.05%

  size_t n_pairs = 0, n_pts = 0, n_fail = 0;
  const vector<DetectorCase> cases = discoverCases();
  for( const DetectorCase &c : cases )
  {
    DetEffG2kPar::ParFile par;
    BOOST_REQUIRE_NO_THROW( par = DetEffG2kPar::parseParFile( c.parPath ) );
    const DetEffG2kPar::ParEfficiency parEff( par );

    for( const EccRun &air_run : c.runs )
    {
      if( air_run.gisPath.empty() )
        continue;

      const GisGeometry air_gis = parseGis( air_run.gisPath );
      if( !air_gis.parsed || (air_gis.pressure <= 0.0) )
        continue;

      // Only point-like sources: a disk/Marinelli air path is not one scalar.
      const string &templ = air_gis.templ;
      const bool is_disk = SpecUtils::icontains( templ, "CIRCULAR_PLANE" )
                        && !std::isnan( air_gis.disk_diam_mm )
                        && (air_gis.disk_diam_mm > 5.0);
      const bool is_point = ( SpecUtils::icontains( templ, "SPHERE" )
                              || SpecUtils::icontains( templ, "POINT" )
                              || SpecUtils::icontains( templ, "CIRCULAR_PLANE" ) )
                            && !is_disk;
      if( !is_point )
        continue;

      // Find the like-named "*_vacuum" run at the same position.
      const string air_stem = SpecUtils::to_lower_ascii_copy(
                                  SpecUtils::filename( air_run.eccPath ) );
      const EccRun *vac_run = nullptr;
      for( const EccRun &r : c.runs )
      {
        const string s = SpecUtils::to_lower_ascii_copy( SpecUtils::filename( r.eccPath ) );
        if( (s.size() > air_stem.size()) && !r.gisPath.empty()
            && SpecUtils::istarts_with( s, air_stem.substr( 0, air_stem.find_last_of('.') ) )
            && SpecUtils::icontains( s, "vacuum" ) )
        {
          vac_run = &r;
          break;
        }
      }
      if( !vac_run )
        continue;

      const GisGeometry vac_gis = parseGis( vac_run->gisPath );
      if( !vac_gis.parsed || (vac_gis.pressure > 0.0) )
        continue;

      // Same position, or the ratio is not an air-only ratio.
      const double vac_sd4 = std::isnan( vac_gis.sd4_mm ) ? 0.0 : vac_gis.sd4_mm;
      const double air_sd4 = std::isnan( air_gis.sd4_mm ) ? 0.0 : air_gis.sd4_mm;
      if( (std::fabs( vac_gis.sd1_mm - air_gis.sd1_mm ) > 1.0e-6)
          || (std::fabs( vac_sd4 - air_sd4 ) > 1.0e-6) )
        continue;

      const vector<pair<double,double>> ecc_air = readEcc( air_run.eccPath );
      const vector<pair<double,double>> ecc_vac = readEcc( vac_run->eccPath );
      if( ecc_air.empty() || ecc_vac.empty() )
        continue;

      // The air path is the source's radial range from the face, which is ~sd1
      //  (the same quantity the grid row uses).
      const double air_path_mm = air_gis.sd1_mm;

      ++n_pairs;
      size_t n_here = 0, n_fail_here = 0;
      double worst = 0.0;
      for( const pair<double,double> &a : ecc_air )
      {
        double vac_eff = -1.0;
        for( const pair<double,double> &v : ecc_vac )
        {
          if( std::fabs( v.first - a.first ) <= 1.0e-3 * std::max( a.first, 1.0 ) )
          {
            vac_eff = v.second;
            break;
          }
        }
        if( !(vac_eff > 0.0) || !(a.second > 0.0) )
          continue;

        const double ref_T = a.second / vac_eff;
        const double our_T = parEff.airTransmission( a.first, air_path_mm,
                                                     air_gis.pressure / 760.0 );
        const double rel = std::fabs( our_T - ref_T ) / ref_T;
        worst = std::max( worst, rel );
        ++n_here;
        if( rel > air_tol )
          ++n_fail_here;
      }//for( each air/vacuum energy pair )

      n_pts += n_here;
      n_fail += n_fail_here;
      BOOST_CHECK_MESSAGE( n_fail_here == 0,
          SpecUtils::filename( air_run.eccPath ) + ": " + std::to_string(n_fail_here)
          + "/" + std::to_string(n_here) + " air-factor points exceeded "
          + std::to_string(air_tol) + " (worst " + std::to_string(worst)
          + ") at air path " + std::to_string(air_path_mm) + " mm" );
      BOOST_TEST_MESSAGE( "  air-only " + SpecUtils::filename( air_run.eccPath )
                          + "  air path " + std::to_string(air_path_mm) + " mm, "
                          + std::to_string(n_here) + " pts, worst "
                          + std::to_string(worst) );
    }//for( each air run )
  }//for( each case )

  BOOST_TEST_MESSAGE( "AirTransmissionOnly: " + std::to_string(n_pairs)
                      + " air/vacuum pairs, " + std::to_string(n_pts) + " points, "
                      + std::to_string(n_fail) + " failures" );
}//BOOST_AUTO_TEST_CASE( AirTransmissionOnly )


BOOST_AUTO_TEST_CASE( SerializationRoundTrip )
{
  if( g_par_dir.empty() )
    return;

  const vector<DetectorCase> cases = discoverCases();
  for( const DetectorCase &c : cases )
  {
    shared_ptr<DetectorPeakResponse> drf;
    BOOST_REQUIRE_NO_THROW( drf = DetEffG2kPar::makeDrfFromFiles( c.parPath, c.detectorTxtPath ) );
    BOOST_REQUIRE( !!drf );

    rapidxml::xml_document<char> doc;
    rapidxml::xml_node<char> *root = doc.allocate_node( rapidxml::node_element, "root" );
    doc.append_node( root );
    BOOST_REQUIRE_NO_THROW( drf->toXml( root, &doc ) );

    rapidxml::xml_node<char> *drf_node = root->first_node( "DetectorPeakResponse" );
    BOOST_REQUIRE( drf_node != nullptr );

    auto restored = make_shared<DetectorPeakResponse>();
    BOOST_REQUIRE_NO_THROW( restored->fromXml( drf_node ) );

    // The DrfSource provenance persists.
    BOOST_CHECK( restored->drfSource() == DetectorPeakResponse::DrfSource::CharacterizationParFile );

    // The CeeLo response survives the round-trip, and reproduces the same
    //  off-axis point efficiency before and after.
    BOOST_REQUIRE( !!restored->ceeloResponse() );

    // So does "no total efficiency" - a reload that quietly reverted to the
    //  bare-kernel tier would re-enable cascade summing on an FEP-only import.
    BOOST_CHECK( !restored->ceeloResponse()->tot_eff.characterized() );
    BOOST_CHECK( !restored->hasAnyTotalEfficiencyInfo() );

    const double d_cm = 25.0, theta = 0.3;
    for( const double E : { 60.0, 200.0, 662.0, 1332.0 } )
    {
      const DetectorPeakResponse::EffEval a
              = drf->fepEfficiencyEval( static_cast<float>(E), theta, 0.0, d_cm*PhysicalUnits::cm );
      const DetectorPeakResponse::EffEval b
              = restored->fepEfficiencyEval( static_cast<float>(E), theta, 0.0, d_cm*PhysicalUnits::cm );
      BOOST_CHECK_MESSAGE( std::fabs(a.value - b.value) <= 1e-9 * std::max(a.value,1e-300) + 1e-30,
          "round-trip eff mismatch at " + std::to_string(E) + " keV: "
          + std::to_string(a.value) + " vs " + std::to_string(b.value) );
    }
  }//for( each case )
}//BOOST_AUTO_TEST_CASE( SerializationRoundTrip )
