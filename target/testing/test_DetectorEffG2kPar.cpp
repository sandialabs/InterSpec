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
      efficiency and the DrfSource provenance;
   6. `GridReproductionByBand`: how well the ASSEMBLED response reproduces the
      grid, per energy band, gated on the file's own uint16 lattice (see that
      case for why that is the only unambiguous truth the file carries).

 Two cases need no vendor files and so run in CI on `--datadir=` alone:
 `LayerMaterialsResolveAnyElement` and `KEdgeSegmentsHaveBothFlanks`.
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


/** No K-edge segment of the assembled eta table may hold a single node, and the
 assembled efficiency must be CONTINUOUS across every edge.

 `EtaTable` interpolates ln eta in ln E segmented at the crystal K-edges, and its
 contract is that both flanks of each edge are nodes.  A segment left with one
 node can only be represented as a constant stub, which freezes ln eta across the
 whole sub-edge interval while `K(E)` keeps moving - fabricating a discontinuity
 that exists nowhere in the file.  That is not hypothetical: with the file's own
 20 energies as the axis, the Ge K-edge at 11.107 keV falls between the 10 and
 12 keV nodes, leaving 10 keV alone in its segment, and the assembled efficiency
 jumped by 3.8e8 across the edge - where the grid's 10 and 12 keV pages are
 bit-identical, so the truth ratio is exactly 1.

 Needs no vendor files (a synthetic grid straddling the edge drives makeDrf), and
 the edge at issue is a property of the crystal, so this runs in CI.  `--datadir=`
 is required for the element/cross-section lookups.
 */
BOOST_AUTO_TEST_CASE( KEdgeSegmentsHaveBothFlanks )
{
  if( g_data_dir.empty() )
  {
    BOOST_TEST_MESSAGE( "KEdgeSegmentsHaveBothFlanks: skipped (needs --datadir=)." );
    return;
  }

  // A germanium detector with a dead layer, so the Ge K-edge at 11.107 keV is in
  //  play, and a grid whose range STRADDLES it: crystal_k_edges keeps an edge only
  //  when it is inside the range by its 1.02/0.98 margins, which 11.107 is for a
  //  10 keV first node.  This is the geometry LAB06 does not have (its grid starts
  //  at 45 keV, so no edge is retained and it never showed the bug).
  const string detector_txt =
    "# 99999 - SYNTHETIC-KEDGE - S/N-TEST-2\n"
    "KEdgeDet,70.0,60.0,0,80.0,150.0,4.0,4.0,26,kedge.par,4, #\n"
    "ge,0.23,5.35, #\n"
    "be,0.5,1.848, #\n"
    "al,1.5,2.70, #\n"
    "ge,0.23,5.35, #\n"
    "al,1.0,2.70, #\n"
    "al,2.0,2.70, #\n"
    "ge,0.23,5.35, #\n"
    "cu,3.0,8.96, #\n"
    "al,5.0,2.70\n";

  istringstream txt( detector_txt );
  const vector<DetEffG2kPar::DetectorDef> defs = DetEffG2kPar::parseDetectorTxt( txt );
  BOOST_REQUIRE_EQUAL( defs.size(), 1 );
  const DetEffG2kPar::DetectorDef &def = defs.front();

  // The file's own energies bracket the edge without flanking it - exactly the
  //  vendor layout that produced the bug.
  DetEffG2kPar::ParFile par;
  par.emin_keV = 10.0;
  par.emax_keV = 1332.0;
  par.energies_keV = { 10.0, 12.0, 16.0, 22.0, 45.0, 100.0, 1332.0 };
  for( size_t e = 0; e < par.energies_keV.size(); ++e )
  {
    DetEffG2kPar::ParGrid g;
    g.ncols = 19;                    // 0..180 deg in 10 deg steps
    g.nrows = 40;
    g.theta_step_rad = 10.0 * pi / 180.0;
    g.r_step = 0.2;                  // ln(mm)
    g.V.resize( static_cast<size_t>(g.nrows) * g.ncols );
    for( int r = 0; r < g.nrows; ++r )
    {
      for( int c = 0; c < g.ncols; ++c )
      {
        // The 10 and 12 keV pages are deliberately IDENTICAL, mirroring the
        //  corpus file: the truth is then exactly flat across the edge, so any
        //  jump the assembled response shows is entirely fabricated.
        const double V = 2000.0 + 40.0*r + 15.0*c
                         + ((e < 2) ? 0.0 : 300.0*static_cast<double>(e-1));
        g.V[r*g.ncols + c] = static_cast<uint16_t>( V );
      }
    }
    par.grids.push_back( g );
  }

  shared_ptr<DetectorPeakResponse> drf;
  BOOST_REQUIRE_NO_THROW( drf = DetEffG2kPar::makeDrf( par, def ) );
  BOOST_REQUIRE( !!drf );
  const shared_ptr<const ceelo::DetectorResponse> resp = drf->ceeloResponse();
  BOOST_REQUIRE( !!resp );

  const ceelo::EtaTable &eta = resp->eta_fep;
  BOOST_REQUIRE( !eta.energies_keV.empty() );

  // The fixture must actually retain an edge, or everything below is vacuous.
  BOOST_REQUIRE_MESSAGE( !eta.edges_keV.empty(),
      "fixture retained no K-edge, so this case would pass without testing anything"
      " - crystal_k_edges keeps an edge only inside the range by its 1.02/0.98"
      " margins; check the crystal material and the grid's energy range" );

  BOOST_TEST_MESSAGE( "  K-edge fixture: " + std::to_string(par.energies_keV.size())
                      + " file nodes -> " + std::to_string(eta.energies_keV.size())
                      + " eta nodes, " + std::to_string(eta.edges_keV.size())
                      + " retained edge(s)" );

  // (1) Structural, tolerance-free: >= 2 nodes strictly each side of every edge.
  //     This is the assertion that cannot silently pass.
  for( const double edge : eta.edges_keV )
  {
    size_t below = 0, above = 0;
    for( const double E : eta.energies_keV )
    {
      if( E < edge )
        ++below;
      else if( E > edge )
        ++above;
    }
    BOOST_TEST_MESSAGE( "  edge " + std::to_string(edge) + " keV: "
                        + std::to_string(below) + " nodes below, "
                        + std::to_string(above) + " above" );
    BOOST_CHECK_MESSAGE( below >= 2 && above >= 2,
        "K-edge " + std::to_string(edge) + " keV has " + std::to_string(below)
        + " nodes below and " + std::to_string(above) + " above; a segment with"
        " fewer than 2 interpolates as a constant stub, freezing ln eta on that"
        " side of the edge and fabricating a discontinuity" );
  }

  // (2) The assembled efficiency is continuous across each edge.  The grid is
  //     flat across 10-12 keV here, so the truth ratio over +-1% is the ratio the
  //     KERNEL contributes - itself modest.  Pre-fix this returned ~3.8e8.
  const Eigen::Vector3d src
      = CeeLoUtils::sourcePositionFromFace( resp->descriptor, 0.0, 0.0, 25.0 );
  for( const double edge : eta.edges_keV )
  {
    const ceelo::EffResult lo = resp->eps_fep_at( edge*0.99, src );
    const ceelo::EffResult hi = resp->eps_fep_at( edge*1.01, src );
    BOOST_REQUIRE_MESSAGE( lo.value > 0.0 && hi.value > 0.0,
        "assembled efficiency vanished at the " + std::to_string(edge) + " keV edge" );

    const DetEffG2kPar::ParEfficiency parEff( par );
    const double t_lo = parEff.efficiency( edge*0.99, 250.0, 0.0 );
    const double t_hi = parEff.efficiency( edge*1.01, 250.0, 0.0 );
    BOOST_REQUIRE( t_lo > 0.0 && t_hi > 0.0 );

    // Compare the JUMP, not the absolute value: the ratio across the edge must
    //  track the grid's own ratio.  A stub segment breaks this by many decades
    //  while leaving each one-sided value superficially plausible.
    const double ours = hi.value / lo.value;
    const double truth = t_hi / t_lo;
    BOOST_TEST_MESSAGE( "  across edge " + std::to_string(edge) + " keV: response steps "
                        + std::to_string(ours) + "x, grid steps " + std::to_string(truth)
                        + "x (pre-fix the response stepped ~3.8e8x)" );
    BOOST_CHECK_MESSAGE( std::fabs( std::log( ours / truth ) ) < 0.05,
        "across the " + std::to_string(edge) + " keV edge the response steps by "
        + std::to_string(ours) + "x but the grid steps by " + std::to_string(truth)
        + "x - a lone-node segment freezes ln eta on one side of the edge" );
  }

  {
    // (2b) The edge-ratio check above is local to one edge; this sweeps the whole
    //  sub-knee range for a jump ANYWHERE.  The grid here is flat below the knee,
    //  so any step the response takes that the grid does not is ours - and a
    //  lone-node segment produces exactly that, at whatever energy it occurs,
    //  which is what this catches that (2) cannot.  Deliberately OUTSIDE the
    //  per-edge loop above: the sweep already tests every edge through
    //  `near_edge`, so running it once per edge would only reset its counters and
    //  repeat the work, reporting a per-iteration figure as if it were global.
    const DetEffG2kPar::ParEfficiency parEff( par );
    //
    //  A +-2% window around each edge is EXCLUDED, and the reason is physics, not
    //  tolerance shopping.  Inside it `K` is mid-cliff: it falls ~9 decades
    //  (2.9e-09 -> 4.5e-18 over 11.11-11.20 keV on this fixture) as the dead
    //  layer's tau jumps through the edge, and `eta` has to cancel that to
    //  reproduce a flat grid.  Cancellation is exact only where the two agree on
    //  an interpolation basis, and they do not: eta is a cubic in ln E while
    //  `K`'s attenuation is exp(-tau) off the mu table's own log-log grid.  The
    //  measured residue is a ~10% bump over ~0.09 keV that recovers to 1.0002 the
    //  moment `K` clears the transition.  On a REAL detector that window is deep
    //  in the sub-16 keV region whose peak efficiency is ~2e-09 (see the band
    //  table in "DetectorEffG2kPar.h"), which is why that band is gated on
    //  ABSOLUTE error and passes; on this synthetic fixture the efficiency there
    //  is much larger (~8e-04) because the fixture's grid is not a real detector's,
    //  so do not read a magnitude off this case.  Fixing the residue means putting
    //  eta on the kernel's basis, i.e. changing the interpolator - out of scope.
    //  The in-window figure is reported rather than hidden, and (2) pins the ratio
    //  across the edge.
    const auto near_edge = [&eta]( const double E ) -> bool
    {
      for( const double ed : eta.edges_keV )
        if( E > ed*0.98 && E < ed*1.02 )
          return true;
      return false;
    };

    double prev_ours = -1.0, prev_truth = -1.0, prev_E = 0.0;
    double worst_drop = 0.0, at_E = 0.0, worst_in_win = 0.0;
    for( double E = 10.0; E <= 45.0 + 1e-9; E *= 1.01 )
    {
      const double v = resp->eps_fep_at( E, src ).value;
      const double t = parEff.efficiency( E, 250.0, 0.0 );
      if( !(v > 0.0) || !(t > 0.0) )
      {
        prev_ours = prev_truth = -1.0;
        prev_E = E;   // reset too, else `near_edge(prev_E)` reads a stale energy
        continue;
      }
      if( prev_ours > 0.0 && t >= prev_truth )   // grid rising: we must rise too
      {
        const double drop = prev_ours / v - 1.0;
        if( near_edge( E ) || near_edge( prev_E ) )
          worst_in_win = std::max( worst_in_win, drop );
        else if( drop > worst_drop ){ worst_drop = drop; at_E = E; }
      }
      prev_ours = v;
      prev_truth = t;
      prev_E = E;
    }
    BOOST_TEST_MESSAGE( "  10-45 keV monotonicity: worst relative DROP where the grid"
                        " rises is " + std::to_string(worst_drop)
                        + " (at " + std::to_string(at_E) + " keV); inside the +-2%"
                        " edge windows, where the K cliff makes exact cancellation"
                        " unachievable, " + std::to_string(worst_in_win) );
    BOOST_CHECK_MESSAGE( worst_drop < 0.005,
        "the assembled response FALLS by " + std::to_string(worst_drop)
        + " across " + std::to_string(at_E) + " keV while the grid rises - a"
        " fabricated discontinuity in the sub-knee region" );

    // The in-window figure is GATED too, just loosely.  It must not be merely
    //  reported: the `reach <= kStubReach` branch in `EtaTable::finalize`
    //  deliberately keeps a stub for a lone FLANK node, and the artifact such a
    //  stub produces lands exactly inside this excluded window - where check (2),
    //  which samples only edge*0.99 and edge*1.01, would step straight over it.
    //  The bound is set by what the two failure modes look like, not by the
    //  measured value: the K-cliff basis mismatch is a ~10% bump, while a stub
    //  freezing ln eta across a segment is the ~95x (9400%) jump this whole change
    //  exists to prevent.  Anything between those is also a defect.
    BOOST_CHECK_MESSAGE( worst_in_win < 0.5,
        "inside the +-2% edge window the response falls by "
        + std::to_string(worst_in_win) + " where the grid rises.  The kernel-basis"
        " residue is ~10%; a drop this large is a frozen segment, not that" );
  }

  // (3) Direct library guard: hand-build an EtaTable whose edge list leaves the
  //     first node alone in its segment, and require eval_ln in the gap to lie
  //     BETWEEN the neighbouring node values rather than frozen at one of them.
  //     Reachable from any producer or any already-serialized response, so the
  //     library must not be able to produce the stub when a neighbour exists.
  {
    ceelo::EtaTable t;
    t.energies_keV = { 10.0, 12.0, 16.0, 22.0 };
    t.cos_thetas = { 1.0 };
    t.edges_keV = { 11.107 };              // between node 0 and node 1: lone node
    t.ln_eta = { -4.0, 0.0, 0.5, 0.6 };    // a big step 10 -> 12 keV
    t.frac_sigma.assign( t.ln_eta.size(), 0.0 );
    BOOST_REQUIRE_NO_THROW( t.finalize() );

    bool clamped = false;
    const double mid = t.eval_ln( 11.0, 1.0, 0.0, clamped );
    BOOST_CHECK_MESSAGE( mid > -4.0 + 1e-9 && mid < 0.0 - 1e-9,
        "EtaTable::eval_ln at 11.0 keV returned " + std::to_string(mid)
        + ", not strictly between the bracketing node values -4.0 and 0.0:"
        " the lone first segment is still being served as a constant stub" );
    BOOST_CHECK_MESSAGE( !clamped,
        "EtaTable::eval_ln clamped at 11.0 keV, which is interior to the node"
        " range - the query fell outside a stub segment's degenerate span" );

    // A genuinely single-node table must still evaluate, as a constant.
    ceelo::EtaTable one;
    one.energies_keV = { 661.657 };
    one.cos_thetas = { 1.0 };
    one.ln_eta = { -0.25 };
    one.frac_sigma.assign( 1, 0.0 );
    BOOST_REQUIRE_NO_THROW( one.finalize() );
    bool one_clamped = false;
    BOOST_CHECK_CLOSE( one.eval_ln( 661.657, 1.0, 0.0, one_clamped ), -0.25, 1e-9 );
  }
}//BOOST_AUTO_TEST_CASE( KEdgeSegmentsHaveBothFlanks )


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


/** How well the assembled response reproduces the grid, measured per energy band -
 the measurement the header's validity table should come from, rather than an
 ad-hoc script.

 THE METRIC, and why this one.  There are three distinct things one could call
 "error against the file", and conflating them is how the previous band table came
 to read as though the low energies were unusable:

   (a) ON the file's own lattice - E a file node, d an exact grid row, theta an
       exact grid column.  `ParEfficiency::efficiency` then returns the stored
       uint16 with NO interpolation in ANY dimension, so this is the file's literal
       content and the only unambiguous truth it has.  THIS is what we gate.
   (b) OFF-node in ENERGY.  Truth is `ParEfficiency`'s log-log PCHIP - our own
       curve.  An in-tree Monte-Carlo arbiter already ruled it more physical than
       the reference tool's interpolation of the same grid (see the
       DetectorEffG2kPar.h header), so it is worth matching, and it is reported.
   (c) OFF-node in d or theta.  Truth then includes `interp_grid_V`'s bilinear
       interpolation in RAW V, which puts a kink on every one of the 440 radial
       rows and 73 angular columns.  That is the vendor's storage scheme, not
       physics - the same category as the `.ecc` between-node gap, which the header
       explains we deliberately do not chase.  Reported, never gated.

 WHAT SETS THE FLOOR.  The file stores eff = 10^(-V/1000) with V an integer, so
 one raw-V step is 0.2305% and +-1 LSB is 0.1153%.  A smooth interpolant through
 lattice-quantized nodes cannot beat that scatter, so an on-lattice result at or
 below ~0.1% is representation-limited, not method-limited, and no tighter claim
 against this file is definable.  The gates below are set just above the measured
 values in that spirit; they are change-detectors at the floor, not headroom.
 */
BOOST_AUTO_TEST_CASE( GridReproductionByBand )
{
  if( g_par_dir.empty() )
    return;

  // Bands chosen to isolate the absorption knee, where the vendor's 20 energy
  //  nodes are genuinely sparse (the truth's own curvature in ln E across them is
  //  ~4.5 at 22-32 keV vs ~0.08 at 200-1000 keV).
  // `gate`/`corner_gate` are RELATIVE bounds on the on-lattice error; `abs_gate`
  //  is an ABSOLUTE bound in efficiency units.  Below 22 keV the absolute bound is
  //  the one that carries the argument, but BOTH are applied there: an absolute
  //  gate alone leaves the relative error completely unpoliced, so a regression
  //  that multiplied it while staying under a floor of ~1e-9 would pass silently.
  //  The sub-22 relative gates are therefore set loosely (~1.5x measured), to
  //  catch a change in kind rather than to assert a tolerance.
  //
  // Why the two lowest bands are gated absolutely.  A relative bound is only
  //  meaningful next to the efficiency it sits on, and below 22 keV this crystal's
  //  stored efficiency is 1e-11 to 1e-9 - seven to ten decades below its ~3e-1
  //  peak.  The worst 10-16 keV locus is a 49% error on eps = 2.2e-09, i.e. an
  //  absolute error of 6e-10, which cannot move any spectrum.  That residual is
  //  also not reachable by node density: it is the crystal-frame lattice reading a
  //  face-frame grid, whose origins differ by `endcap_front_offset_cm`, so the
  //  angular skew grows from +0.08 deg at 300 cm to +3 deg at 8 cm against the
  //  grid's 2.5 deg columns, shearing the file's sharp off-axis shoulder (a factor
  //  of 44 over 57 deg at 16 keV / 9 cm) into the radial direction.  A separable
  //  product lattice cannot align with a sheared kink; refining the ladder 4x
  //  (nd 36 -> 144) left the worst 10-16 locus at 0.067, and refining the angular
  //  axis alone made it slightly worse.  So the absolute gate is what states the
  //  guarantee that matters, instead of a 50% relative tolerance that would look
  //  like a measurement while policing nothing.  The loose relative gate is kept
  //  alongside it so the error cannot change by an order of magnitude unnoticed.
  //
  // Why the >= 22 keV relative gates carry headroom over the achieved figures.
  //  The achieved worst is not monotone in ladder density: rung phase against the
  //  file's fixed 3.03%-per-row grid matters as much as rung spacing, and no
  //  crystal-frame ladder can be row-aligned in the face frame anyway.  Measured
  //  across subdivisions 2/3/4, one band swung 0.031 / 0.0096 / 0.021.  Gates sit
  //  ~1.5x above the achieved value so a benign re-phasing is not a red build,
  //  while a real regression still trips them.
  struct Band { double lo, hi; const char *name; double gate, corner_gate,
                abs_gate, corner_abs_gate; };
  const Band bands[] = {
    {   10.0,    16.0, "10-16",  0.75,  1.45, 1.8e-09, 1.3e-08 },
    {   16.0,    22.0, "16-22",  0.065, 0.92, 4.6e-10, 6.0e-08 },
    {   22.0,    32.0, "22-32",  0.040, 0.40, 0.0,     0.0     },
    {   32.0,    45.0, "32-45",  0.025, 0.20, 0.0,     0.0     },
    {   45.0,    60.0, "45-60",  0.008, 0.15, 0.0,     0.0     },
    {   60.0,   200.0, "60-200", 0.006, 0.10, 0.0,     0.0     },
    {  200.0,  1000.0, "200-1k", 0.009, 0.06, 0.0,     0.0     },
    { 1000.0,  7100.0, "1k-7k",  0.015, 0.12, 0.0,     0.0     },
  };
  const size_t nbands = sizeof(bands)/sizeof(bands[0]);

  // The uint16 quantization figures, so the numbers below are interpretable.
  const double lsb_rel = std::pow( 10.0, 0.5/1000.0 ) - 1.0;   // +-1 LSB, ~0.1153%
  BOOST_TEST_MESSAGE( "  uint16 storage floor: one raw-V step "
                      + std::to_string( std::pow(10.0,1.0/1000.0) - 1.0 )
                      + ", +-1 LSB " + std::to_string(lsb_rel) );

  // Per REGIME as well as per band: counting only MAIN would leave every CORNER
  //  gate provable by an empty bucket.
  size_t global_count[2][16] = {{0}};
  bool any_case = false;

  const vector<DetectorCase> cases = discoverCases();
  for( const DetectorCase &c : cases )
  {
    DetEffG2kPar::ParFile par;
    BOOST_REQUIRE_NO_THROW( par = DetEffG2kPar::parseParFile( c.parPath ) );
    BOOST_REQUIRE( !par.grids.empty() && !par.energies_keV.empty() );

    ifstream txt( c.detectorTxtPath.c_str(), ios::in | ios::binary );
    BOOST_REQUIRE( txt.is_open() );
    const vector<DetEffG2kPar::DetectorDef> defs = DetEffG2kPar::parseDetectorTxt( txt );
    const DetEffG2kPar::DetectorDef def
                = DetEffG2kPar::selectDetectorDef( defs, SpecUtils::filename( c.parPath ) );

    shared_ptr<DetectorPeakResponse> drf;
    BOOST_REQUIRE_NO_THROW( drf = DetEffG2kPar::makeDrf( par, def ) );
    const shared_ptr<const ceelo::DetectorResponse> resp = drf->ceeloResponse();
    BOOST_REQUIRE( !!resp );
    any_case = true;

    const DetEffG2kPar::ParEfficiency parEff( par );
    const ceelo::EtaTable &eta = resp->eta_fep;
    const ceelo::NearFieldModel &nf = resp->near_field;
    const double r_step = par.grids.front().r_step;
    const double t_step = par.grids.front().theta_step_rad;

    BOOST_TEST_MESSAGE( "  " + SpecUtils::filename(c.dir) + ": "
        + std::to_string(par.energies_keV.size()) + " file energies -> eta "
        + std::to_string(eta.energies_keV.size()) + " x "
        + std::to_string(eta.cos_thetas.size()) + ", near field "
        + std::to_string(nf.energies_keV.size()) + " x "
        + std::to_string(nf.cos_thetas.size()) + " x "
        + std::to_string(nf.dists_cm.size()) + ", serialized "
        + std::to_string(resp->to_xml_string().size()) + " bytes" );

    // Probe positions: exact grid rows over the near field's span, and exact grid
    //  columns.  `sourcePositionFromFace` is the inverse of makeDrf's face_sample,
    //  so |src - face| == d_cm and the polar angle == theta identically; the
    //  truth-side arguments are then just 10*d_cm and theta, with no round trip
    //  and no chance of the axial-vs-radial inversion EccMatch warns about.  Do
    //  NOT use hypot here for the same reason.
    vector<double> d_rows;                       // exact grid rows, cm
    {
      const double d_lo_cm = nf.dists_cm.empty() ? 1.0 : nf.dists_cm.front();
      const double d_hi_cm = nf.dists_cm.empty() ? 300.0 : nf.dists_cm.back();
      const int i_lo = static_cast<int>( std::ceil( std::log( d_lo_cm*10.0 ) / r_step ) );
      const int i_hi = static_cast<int>( std::floor( std::log( d_hi_cm*10.0 ) / r_step ) );
      for( int i = i_lo; i <= i_hi; i += 4 )     // every 4th row: ~12% apart
        d_rows.push_back( std::exp( static_cast<double>(i) * r_step ) / 10.0 );
    }
    vector<double> thetas;                       // exact grid columns on [0, 90]
    for( int j = 0; j*t_step <= 0.5*pi + 1e-9; j += 2 )
      thetas.push_back( static_cast<double>(j) * t_step );

    BOOST_REQUIRE_MESSAGE( !d_rows.empty() && !thetas.empty(),
        "no on-lattice probe positions were generated - the grid step or the near"
        " field's distance span is not what this test assumes" );

    // Per-band accumulators: (a) on-lattice (gated), (b) off-node in energy only
    //  (reported), (c) off-node in d and theta too (reported).
    //
    // Each is split by REGIME, because two of the three regimes are limited by
    //  different things and folding them together hides both:
    //    MAIN   d >= 5 cm and theta <= 75 deg - the regime a user query lands in,
    //           limited only by our own table resolution.  Gated per band.
    //    CORNER near contact (d < 5 cm) or grazing (theta > 75 deg) - limited by
    //           the distance ladder's rung spacing against the grid's 440 radial
    //           rows and by the face-frame geometry at grazing incidence, which is
    //           entangled with the deferred theta>90 work.  Measured and reported
    //           with a loose change-detector gate, never folded into MAIN: a
    //           known-deferred regime must not mask an on-axis regression, and an
    //           on-axis number must not be flattered by being averaged with it.
    enum Regime { kMain = 0, kCorner = 1, kNRegimes = 2 };
    double worst_a[kNRegimes][16] = {{0}}, worst_b[kNRegimes][16] = {{0}};
    double worst_c[kNRegimes][16] = {{0}}, sum2_a[kNRegimes][16] = {{0}};
    size_t n_a[kNRegimes][16] = {{0}}, n_b[kNRegimes][16] = {{0}};
    // On-lattice ABSOLUTE error, and the efficiency the worst RELATIVE error sits
    //  on, so every relative figure below can be read next to its magnitude.
    double absw_a[kNRegimes][16] = {{0}}, eps_at_worst[kNRegimes][16] = {{0}};
    double eps_max[kNRegimes][16] = {{0}};
    size_t n_c[kNRegimes][16] = {{0}};

    auto band_of = [&]( const double E ) -> size_t
    {
      for( size_t b = 0; b < nbands; ++b )
      {
        if( E >= bands[b].lo && E < bands[b].hi )
          return b;
      }
      return nbands;
    };

    // Off-node energies: geometric midpoints of adjacent eta nodes are the
    //  worst-case location for cubic-Hermite error.  Assert they really are
    //  off-node, so a coincidence cannot turn this into a second node check.
    auto is_eta_node = [&eta]( const double E ) -> bool
    {
      for( const double n : eta.energies_keV )
        if( std::fabs( std::log( E / n ) ) < 1.0e-4 )
          return true;
      return false;
    };
    vector<double> e_offnode;
    for( size_t i = 0; (i + 1) < eta.energies_keV.size(); i += 17 )
    {
      const double E = std::sqrt( eta.energies_keV[i] * eta.energies_keV[i+1] );
      if( !is_eta_node( E ) )
        e_offnode.push_back( E );
    }

    for( const double theta : thetas )
    {
      for( const double d_cm : d_rows )
      {
        const Eigen::Vector3d src
            = CeeLoUtils::sourcePositionFromFace( resp->descriptor, theta, 0.0, d_cm );
        const double d_face_mm = d_cm * 10.0;
        const size_t rg = ((d_cm < 5.0) || (theta > 75.0*pi/180.0)) ? kCorner : kMain;

        // One quadrature per position, reused across every energy: the header
        //  documents this overload as ~30x faster and bit-identical.
        const shared_ptr<const DetectorPeakResponse::PositionedQuadrature> quad
            = drf->apertureQuadrature( theta, 0.0, d_cm*PhysicalUnits::cm );

        // Returns the signed-magnitude relative error, or -1 where the file has
        //  nothing to be compared TO (the V==0 no-data corner behind the endcap).
        //  `q` must be the quadrature traced at (th_truth, d_truth_mm) - a
        //  mismatch silently answers for the wrong position (and trips the
        //  developer-check assert, which is how this was caught).
        //  `absw`/`eps_w`/`eps_hi`, when given, additionally track the worst
        //  ABSOLUTE error, the efficiency the worst RELATIVE error sits on, and the
        //  band's peak efficiency - the magnitude context that makes a relative
        //  figure interpretable (and, below 22 keV, what is gated instead).
        auto compare = [&]( const double E, double &worst, double *sum2,
                            size_t &n, const double d_truth_mm, const double th_truth,
                  const shared_ptr<const DetectorPeakResponse::PositionedQuadrature> &q,
                            double *absw = nullptr, double *eps_w = nullptr,
                            double *eps_hi = nullptr )
                            -> double
        {
          bool no_data = false;
          const double truth = parEff.efficiency( E, d_truth_mm, th_truth, &no_data );
          if( no_data || !(truth > 1.0e-30) )
            return -1.0;
          const DetectorPeakResponse::EffEval ev = drf->fepEfficiencyEval(
                          static_cast<float>(E), th_truth, 0.0,
                          (d_truth_mm/10.0)*PhysicalUnits::cm, q );
          if( !(ev.value > 0.0) )
            return -1.0;
          const double rel = std::fabs( ev.value/truth - 1.0 );
          if( eps_w && (rel > worst) )
            *eps_w = truth;
          worst = std::max( worst, rel );
          if( absw )
            *absw = std::max( *absw, std::fabs( ev.value - truth ) );
          if( eps_hi )
            *eps_hi = std::max( *eps_hi, truth );
          if( sum2 )
            *sum2 += rel*rel;
          ++n;
          return rel;
        };

        // (a) ON the lattice: file-node energies at this exact row and column.
        for( const double E : par.energies_keV )
        {
          const size_t b = band_of( E );
          if( b >= nbands )
            continue;
          const double rel = compare( E, worst_a[rg][b], &sum2_a[rg][b], n_a[rg][b],
                                      d_face_mm, theta, quad, &absw_a[rg][b],
                                      &eps_at_worst[rg][b], &eps_max[rg][b] );
          if( rel >= 0.0 )
            ++global_count[rg][b];
        }

        // (b) OFF node in ENERGY only - still an exact row and column, so the
        //     same quadrature applies.
        for( const double E : e_offnode )
        {
          const size_t b = band_of( E );
          if( b < nbands )
            compare( E, worst_b[rg][b], nullptr, n_b[rg][b], d_face_mm, theta, quad );
        }

        // (c) OFF node in d and theta as well: half a grid row out and half a
        //     column over, where `interp_grid_V`'s bilinear-in-raw-V kinks live.
        //     A different position, so it needs its own traced quadrature.
        const double d_mid_mm = d_face_mm * std::exp( 0.5*r_step );
        const double th_mid = std::min( theta + 0.5*t_step, 0.5*pi );
        const shared_ptr<const DetectorPeakResponse::PositionedQuadrature> quad_mid
            = drf->apertureQuadrature( th_mid, 0.0, (d_mid_mm/10.0)*PhysicalUnits::cm );
        for( const double E : e_offnode )
        {
          const size_t b = band_of( E );
          if( b < nbands )
            compare( E, worst_c[rg][b], nullptr, n_c[rg][b], d_mid_mm, th_mid, quad_mid );
        }
      }//for( each grid row )
    }//for( each grid column )

    // Report every metric for both regimes; gate MAIN per band, CORNER loosely.
    for( size_t rg = 0; rg < kNRegimes; ++rg )
    {
      BOOST_TEST_MESSAGE( string("    ") + ((rg == kMain)
              ? "MAIN (d >= 5 cm, theta <= 75 deg) - gated per band"
              : "CORNER (d < 5 cm or theta > 75 deg) - ladder/frame limited, loose gate" ) );
      BOOST_TEST_MESSAGE( "      band      n     ON-lattice  (rms)      off-node-E"
                          "  off-node-d,theta   abs@worst   eps@worst    eps_max" );

      for( size_t b = 0; b < nbands; ++b )
      {
        if( !n_a[rg][b] && !n_b[rg][b] )
          continue;   // LAB06 starts at 45 keV: its four low bands are legitimately empty
        const double rms_a = n_a[rg][b]
            ? std::sqrt( sum2_a[rg][b]/static_cast<double>(n_a[rg][b]) ) : 0.0;
        char line[320];
        std::snprintf( line, sizeof(line),
            "      %-8s %6zu   %9.5f (%7.5f)   %9.5f    %9.5f   %9.2e  %9.2e  %9.2e",
            bands[b].name, n_a[rg][b], worst_a[rg][b], rms_a,
            n_b[rg][b] ? worst_b[rg][b] : 0.0, n_c[rg][b] ? worst_c[rg][b] : 0.0,
            absw_a[rg][b], eps_at_worst[rg][b], eps_max[rg][b] );
        BOOST_TEST_MESSAGE( line );

        if( !n_a[rg][b] )
          continue;

        // Absolute gate where one is set (below 22 keV); the relative gate below
        //  then applies IN ADDITION, not instead.
        const double abs_gate = (rg == kMain) ? bands[b].abs_gate
                                             : bands[b].corner_abs_gate;
        if( abs_gate > 0.0 )
        {
          char msg[512];
          std::snprintf( msg, sizeof(msg),
              "%s keV, %s: worst ON-LATTICE ABSOLUTE error %.3e over %zu probes"
              " exceeds %.3e.  This band is gated absolutely because its stored"
              " efficiency (peak %.2e here) is many decades below the crystal's"
              " ~3e-1 peak, where a relative bound polices nothing; the worst"
              " relative error was %.5f on eps %.2e.  The limit is the crystal-frame"
              " lattice reading a face-frame grid, not node density",
              bands[b].name, (rg == kMain) ? "MAIN" : "CORNER", absw_a[rg][b],
              n_a[rg][b], abs_gate, eps_max[rg][b], worst_a[rg][b],
              eps_at_worst[rg][b] );
          BOOST_CHECK_MESSAGE( absw_a[rg][b] <= abs_gate, msg );
          // no `continue` - the relative gate below is a second, looser guard
        }

        const double gate = (rg == kMain) ? bands[b].gate : bands[b].corner_gate;
        BOOST_CHECK_MESSAGE( worst_a[rg][b] <= gate,
            string(bands[b].name) + " keV, "
            + ((rg == kMain) ? "MAIN" : "CORNER") + ": worst ON-LATTICE error "
            + std::to_string(worst_a[rg][b]) + " over " + std::to_string(n_a[rg][b])
            + " probes exceeds " + std::to_string(gate)
            + ".  On the lattice the file returns its stored uint16 with no"
            " interpolation anywhere, so this is our representation error against"
            " the file's literal content" );
      }
    }

    // A band gate that passes because its bucket is empty proves nothing: any
    //  detector whose grid reaches 12 keV must populate the lowest band.
    if( par.energies_keV.front() <= 12.0 )
      BOOST_CHECK_MESSAGE( n_a[kMain][0] > 0,
          "the 10-16 keV band has no on-lattice MAIN probes even though the grid starts at "
          + std::to_string(par.energies_keV.front()) + " keV - the probe list is broken,"
          " not the band empty" );
  }//for( each case )

  // Across all detectors every band must have been exercised somewhere, or the
  //  sweep silently stopped covering part of the range.
  if( any_case )
  {
    for( size_t b = 0; b < nbands; ++b )
    {
      // `Regime` is scoped to the per-case block above; 0 is kMain, 1 is kCorner.
      for( size_t rg = 0; rg < 2; ++rg )
        BOOST_CHECK_MESSAGE( global_count[rg][b] > 0,
            string("band ") + bands[b].name + " keV, "
            + ((rg == 0) ? "MAIN" : "CORNER") + " was never probed by ANY"
            " detector - a per-band gate over an empty bucket passes without"
            " measuring" );
    }
  }
}//BOOST_AUTO_TEST_CASE( GridReproductionByBand )


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
      //  tables subsample the grid (a 60-rung near-field distance ladder against
      //  the grid's 440 radial rows) and interpolate on their own manifold, so a
      //  query that is NOT exactly a fill node differs from the raw-grid
      //  bilinear by the table-vs-grid gap.  Note the cos-theta nodes being the
      //  grid's own 37 columns does NOT make the angular axis node-exact, as an
      //  earlier version of this comment claimed: those nodes are CRYSTAL-frame
      //  cosines while the grid is indexed in the ENDCAP-FACE frame, and the
      //  origins differ by `endcap_front_offset_cm` (a +0.28 deg skew at 83 cm
      //  against 2.5 deg columns, +3.0 at 8 cm).  Unlike the vendor's `.ecc`
      //  between-node values, that gap is OURS, and it is what
      //  `GridReproductionByBand` measures per band.
      //  This check's job is to catch a gross frame/convention error in the
      //  dispatch path (which would show up as many-percent off-axis, as the
      //  pre-fix slant sampling did).  Probe at energy nodes so only the spatial
      //  dimension is between-node - including a LOW one, which is what turns the
      //  old "probes a mid energy rather than sweeping" caveat into a real check:
      //  the 10-16 keV end was where the fabricated K-edge discontinuity lived.
      if( gis.pressure <= 0.0 && !par.energies_keV.empty() )
      {
        const double probes[] = { par.energies_keV.front(),
                                  par.energies_keV[ par.energies_keV.size()/2 ] };
        for( const double E_node : probes )
        {
          const double eff_eval = parEff.efficiency( E_node, dist_mm, theta );
          if( !(eff_eval > 0.0) )
            continue;
          const DetectorPeakResponse::EffEval ev = drf->fepEfficiencyEval(
                      static_cast<float>(E_node), theta, 0.0, (dist_mm/10.0)*PhysicalUnits::cm );

          // Relative OR absolute, whichever is looser.  A purely relative bound
          //  is not meaningful where the file's own efficiency is ~1e-9 (10 keV at
          //  a grazing angle is eight decades below this crystal's peak), and that
          //  is exactly the regime `GridReproductionByBand` shows is limited by
          //  the crystal/face frame skew rather than by anything this check can
          //  police.  The absolute floor is still far below the defects this check
          //  exists to catch: the pre-fix slant-sampling bug was 92% at 45 deg and
          //  8427% at 84 deg, and the fabricated K-edge discontinuity - the reason
          //  a LOW probe is here at all - was ~95x, i.e. ~1e-7 absolute at this
          //  locus.
          //
          //  The floor is 1e-10 rather than a rounder 1e-8 because of how little
          //  the low probe is worth otherwise: at 10 keV and a grazing angle
          //  `eff_eval` is itself ~1e-9, so a 1e-8 floor would tolerate an absolute
          //  error TEN TIMES the entire efficiency at that point, and a recurrence
          //  of the K-edge defect merely 5x smaller than the original would pass.
          //  1e-10 keeps the floor below the value being checked while staying far
          //  above double-rounding on this path.
          const double rel_tol = 0.015 * std::max(eff_eval,1e-300);
          const double abs_tol = 1.0e-10;
          const double diff = std::fabs(ev.value - eff_eval);
          char buf[256];
          snprintf( buf, sizeof(buf), "%s: DRF dispatch %.6g != evaluator %.6g"
                    " (rel %.4f, abs %.3g) at node %.1f keV, %.1f mm, %.1f deg",
                    tag.c_str(), ev.value, eff_eval,
                    diff/std::max(eff_eval,1e-300), diff, E_node, dist_mm,
                    theta*180.0/pi );
          BOOST_CHECK_MESSAGE( diff <= std::max(rel_tol,abs_tol) + 1e-30, buf );
        }
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
