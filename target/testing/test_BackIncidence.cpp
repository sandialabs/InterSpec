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

/* Detector responses for sources BESIDE and BEHIND the detector - behind the plane of its face.

 The truth is Monte Carlo through the same detector description the responses are built from, at
 points placed in CLEARANCE coordinates: along a direction (theta, phi) from the crystal-face
 origin, a clearance s past the detector's outer surface
 (ceelo::GeometryDescriptor::Envelope::exit_distance_cm).  So no point ever lands inside the
 detector, and "near" means near the detector's outside whatever side it is on.

 Each detector carries an APPROXIMATE back structure - two layers behind the crystal, as much as a
 .par DETECTOR.txt record describes - and a DETAILED variant with more of what really sits there,
 used only to measure how much the approximation costs.

 Cases:
   GenerateBackPointTruth  - opt-in (--back-truth-out=<dir>): writes the Monte-Carlo truth banks and
                             the approximate-vs-detailed back comparison per detector.  CPU-hours;
                             results are committed to test_data/ceelo_drf/back_truth/.
   GenerateBackCoveredResponses - opt-in (--back-truth-out=<dir>): a full-MC response WITH Coarse
                             side & back coverage per detector - <name>_fullmc_coarse_response.xml.
                             The uncovered response the measurements also compare is that one with
                             its back columns dropped (uncovered_response): coverage changes no front
                             node, value or seed.
   GenerateVolumetricSceneTruth - opt-in (--back-truth-out=<dir>): Monte Carlo of uniform water
                             sources beside and behind gem35_70 and nai3x3 (the scenes of
                             test_VolumetricLinePath's ExtendedLinePathVsReference) -
                             <name>_volumetric.csv.
   VolumetricSceneMeasure  - InterSpec's line path, with each response, against that truth.
   BackEnvelopeMeasure     - compares the responses against the committed truth, by region
                             (beside / behind the crystal's back plane; near / far), and reports.

 Arguments (after `--`): --datadir= --testfiledir= [--back-truth-out=<dir>] [--threads=N]
                         [--only=<detector-name-substring>]
                         [--measure-out=<dir>]  (BackEnvelopeMeasure: per-row CSV of every FEP comparison)
 */

#include "InterSpec_config.h"

#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <algorithm>

#define BOOST_TEST_MODULE BackIncidence_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/AngleOutxImport.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/GammaInteractionCalc_imp.hpp"

#include "io/ResponseKernel.h"
#include "io/DetectorResponse.h"
#include "io/ResponseGenerator.h"
#include "materials/Material.h"
#include "efficiency/EfficiencyCalculator.h"


using namespace std;
using namespace boost::unit_test;

namespace
{
  string g_data_dir, g_test_data_dir, g_truth_out, g_only, g_measure_out;
  unsigned g_threads = 3;
  // `--back-near-angles=93,115,155,180`: also measure the covered responses with their side & back
  //  near-field tables cut to these angles (no new MC) - the angle-subset A/B.
  vector<double> g_back_near_angles;

  /** The detector body InterSpec's MC tool assumes behind the crystal face (GenerationOptions::
   back_body_depth_cm); the gates below do not judge sources inside it. */
  const double sm_body_depth_cm = 15.0;

  struct Fixture
  {
    Fixture()
    {
      const int argc = framework::master_test_suite().argc;
      char ** const argv = framework::master_test_suite().argv;
      for( int i = 1; i < argc; ++i )
      {
        const string arg = argv[i];
        if( SpecUtils::istarts_with( arg, "--datadir=" ) )
          g_data_dir = arg.substr( 10 );
        else if( SpecUtils::istarts_with( arg, "--testfiledir=" ) )
          g_test_data_dir = arg.substr( 14 );
        else if( SpecUtils::istarts_with( arg, "--back-truth-out=" ) )
          g_truth_out = arg.substr( 17 );
        else if( SpecUtils::istarts_with( arg, "--threads=" ) )
          g_threads = static_cast<unsigned>( std::max( 1, std::stoi( arg.substr( 10 ) ) ) );
        else if( SpecUtils::istarts_with( arg, "--only=" ) )
          g_only = arg.substr( 7 );
        else if( SpecUtils::istarts_with( arg, "--measure-out=" ) )
          g_measure_out = arg.substr( 14 );
        else if( SpecUtils::istarts_with( arg, "--back-near-angles=" ) )
        {
          vector<string> fields;
          SpecUtils::split( fields, arg.substr( 19 ), "," );
          for( const string &f : fields )
            g_back_near_angles.push_back( std::stod( f ) );
        }
      }
      if( !g_data_dir.empty() )
        InterSpec::setStaticDataDirectory( g_data_dir );
    }
  };//struct Fixture

  string truth_dir()
  {
    return SpecUtils::append_path( g_test_data_dir, "ceelo_drf/back_truth" );
  }

  string read_file( const string &path )
  {
    ifstream f( path.c_str(), ios::in | ios::binary );
    if( !f.is_open() )
      return string();
    return string( (istreambuf_iterator<char>(f)), istreambuf_iterator<char>() );
  }


  /** The response a "None" run gives: side & back coverage only adds the back (cos theta < 0)
   columns of the eta tables - every front node, value and seed is the same - so dropping them
   recovers it.  Saves committing an uncovered fixture next to each covered one. */
  shared_ptr<ceelo::DetectorResponse> uncovered_response( const ceelo::DetectorResponse &covered )
  {
    const shared_ptr<ceelo::DetectorResponse> r = ceelo::DetectorResponse::from_xml_string( covered.to_xml_string() );
    for( ceelo::EtaTable *t : { &r->eta_fep, &r->tot_eff.eta_tot } )
    {
      size_t nb = 0;
      while( (nb < t->cos_thetas.size()) && (t->cos_thetas[nb] < 0.0) )
        ++nb;
      if( !nb )
        continue;

      const size_t nc = t->cos_thetas.size(), np = t->phis_deg.empty() ? 1 : t->phis_deg.size();
      vector<double> ln_eta, frac_sigma;
      for( size_t e = 0; e < t->energies_keV.size(); ++e )
        for( size_t c = nb; c < nc; ++c )
          for( size_t p = 0; p < np; ++p )
          {
            ln_eta.push_back( t->ln_eta[t->index( e, c, p )] );
            frac_sigma.push_back( t->frac_sigma[t->index( e, c, p )] );
          }
      t->cos_thetas.erase( begin(t->cos_thetas), begin(t->cos_thetas) + nb );
      t->ln_eta = std::move( ln_eta );
      t->frac_sigma = std::move( frac_sigma );
    }//for( both eta tables )

    if( r->provenance.back_coverage )
      r->provenance.back_coverage = ceelo::BackCoverage::None;
    // ...and the side & back near-field tables, which only exist with coverage.
    r->back_near_field = ceelo::NearFieldModel{};
    r->tot_eff.back_near_field = ceelo::NearFieldModel{};
    r->provenance.back_body_depth_cm = 0.0;
    r->finalize();
    return r;
  }//uncovered_response(...)


  /** `covered` without its side & back near-field tables: every other table is the same nodes and
   values (the tables are their own generation stage), so this is the clean A/B of what they do. */
  shared_ptr<ceelo::DetectorResponse> without_back_near( const ceelo::DetectorResponse &covered )
  {
    const shared_ptr<ceelo::DetectorResponse> r = ceelo::DetectorResponse::from_xml_string( covered.to_xml_string() );
    r->back_near_field = ceelo::NearFieldModel{};
    r->tot_eff.back_near_field = ceelo::NearFieldModel{};
    r->provenance.back_body_depth_cm = 0.0;
    r->finalize();
    return r;
  }//without_back_near(...)


  /** `covered` with its side & back near-field tables cut to the angle columns nearest `thetas_deg`
   - what generating with those angles alone would have given, without new MC. */
  shared_ptr<ceelo::DetectorResponse> back_near_subset( const ceelo::DetectorResponse &covered,
                                                        const vector<double> &thetas_deg )
  {
    const shared_ptr<ceelo::DetectorResponse> r = ceelo::DetectorResponse::from_xml_string( covered.to_xml_string() );
    for( ceelo::NearFieldModel *m : { &r->back_near_field, &r->tot_eff.back_near_field } )
    {
      if( m->empty() )
        continue;
      vector<size_t> keep;
      for( size_t c = 0; c < m->cos_thetas.size(); ++c )
        for( const double th : thetas_deg )
          if( std::fabs( std::acos( m->cos_thetas[c] ) * 180.0 / M_PI - th ) < 0.5 )
            keep.push_back( c );
      std::sort( begin(keep), end(keep) );
      keep.erase( std::unique( begin(keep), end(keep) ), end(keep) );
      const size_t ne = m->energies_keV.size(), nd = m->dists_cm.size();
      vector<double> cts, ln_n, frac_sigma;
      for( const size_t c : keep )
        cts.push_back( m->cos_thetas[c] );
      for( size_t e = 0; e < ne; ++e )
        for( const size_t c : keep )
          for( size_t d = 0; d < nd; ++d )
          {
            ln_n.push_back( m->ln_n[m->index( e, c, d )] );
            frac_sigma.push_back( m->frac_sigma[m->index( e, c, d )] );
          }
      m->cos_thetas = cts;
      m->ln_n = ln_n;
      m->frac_sigma = frac_sigma;
      m->break_cos_thetas = { cts.front(), cts.back() };
    }
    r->finalize();
    return r;
  }//back_near_subset(...)


  /** The committed covered full-MC response of detector `name`, and what it implies: the uncovered
   one ("fullmc_none"), and - when it carries side & back near-field tables - itself without them
   ("fullmc_coarse_nonear") and, with --back-near-angles, cut to those angles ("fullmc_coarse_subset");
   the response itself is "fullmc_coarse".  Nothing when `dir` has no fixture for it. */
  vector<pair<const char *,shared_ptr<const ceelo::DetectorResponse>>> full_mc_responses( const string &dir,
                                                                                        const string &name )
  {
    const string xml = read_file( SpecUtils::append_path( dir, name + "_fullmc_coarse_response.xml" ) );
    if( xml.empty() )
      return {};
    const shared_ptr<const ceelo::DetectorResponse> covered = ceelo::DetectorResponse::from_xml_string( xml );
    vector<pair<const char *,shared_ptr<const ceelo::DetectorResponse>>> out{ { "fullmc_none", uncovered_response( *covered ) } };
    if( !covered->back_near_field.empty() )
    {
      out.emplace_back( "fullmc_coarse_nonear", without_back_near( *covered ) );
      if( !g_back_near_angles.empty() )
        out.emplace_back( "fullmc_coarse_subset", back_near_subset( *covered, g_back_near_angles ) );
    }
    out.emplace_back( "fullmc_coarse", covered );
    return out;
  }//full_mc_responses(...)


  /** Whether `src` is inside the detector body InterSpec assumes (sm_body_depth_cm) - outside the
   described detector, but where a real one keeps its cryostat or electronics. */
  bool inside_assumed_body( const ceelo::GeometryDescriptor &gd, const Eigen::Vector3d &src )
  {
    ceelo::GeometryDescriptor::Envelope body = gd.outer_envelope();
    body.z_max_cm = std::max( body.z_max_cm, sm_body_depth_cm );
    const ceelo::ConvexSolid s = body.solid();
    return !((src - s.project( src )).norm() > 0.0);
  }

  /** The detectors the side & back near-field gates are judged on.  nai_log is reported only: seen
   from its face its back occupies the last few degrees, its far-field back nodes sit within ~2a, and
   it is a box whose phi = 45 deg reuses the phi = 0 tables. */
  bool is_gated_detector( const string &name )
  {
    return (name == "gem35_70") || (name == "detective_x") || (name == "nai3x3");
  }


  // ---------------------------------------------------------------------------------------------
  //  The detectors
  // ---------------------------------------------------------------------------------------------

  int add_material( ceelo::GeometryDescriptor &gd, const ceelo::Material &m )
  {
    gd.materials.push_back( ceelo::MaterialSpec::from( m ) );
    return static_cast<int>( gd.materials.size() - 1 );
  }

  int add_vacuum( ceelo::GeometryDescriptor &gd )
  {
    for( size_t i = 0; i < gd.materials.size(); ++i )
    {
      if( SpecUtils::icontains( gd.materials[i].name, "vacuum" ) )
        return static_cast<int>( i );
    }
    gd.materials.push_back( CeeLoUtils::vacuumMaterialSpec() );
    return static_cast<int>( gd.materials.size() - 1 );
  }

  /** Inserts, at layer position `k`, a layer with neither front nor side wall whose back plate is
   `t_cm` - a disk behind everything inside it, as wide as that.  Layers outside it keep their own
   back plates (their z_end moves out with it). */
  void insert_back_layer( ceelo::GeometryDescriptor &gd, const size_t k, const int material_index,
                          const double t_cm )
  {
    ceelo::LayerSpec l;
    l.material_index = material_index;
    l.front_thickness_cm = 0.0;
    l.side_thickness_cm = 0.0;
    l.z_start_cm = 0.0;
    l.z_end_cm = gd.dimensions_cm.back();   // placeholder; sized below
    gd.layers.insert( gd.layers.begin() + static_cast<std::ptrdiff_t>(k), l );
    gd.set_layer_back_thickness_cm( k, t_cm );
  }

  /** Appends an outermost back-only layer. */
  void append_back_layer( ceelo::GeometryDescriptor &gd, const int material_index, const double t_cm )
  {
    insert_back_layer( gd, gd.layers.size(), material_index, t_cm );
  }


  ceelo::GeometryDescriptor gem35_70_base()
  {
    const string path = SpecUtils::append_path( SpecUtils::append_path( g_test_data_dir, "det_eff" ),
                                                "Angle-example-efficiency.outx" );
    ifstream in( path.c_str(), ios::in | ios::binary );
    BOOST_REQUIRE_MESSAGE( in.is_open(), "Failed to open " << path );
    const AngleOutxContents contents = DetectorPeakResponse::parseAngleOutxFileFull( in );
    BOOST_REQUIRE( contents.hasGeometry );
    vector<string> warnings;
    return CeeLoUtils::buildAngleGeometry( contents, warnings );
  }

  ceelo::GeometryDescriptor golden_descriptor( const string &preset )
  {
    const string xml = read_file( SpecUtils::append_path( g_test_data_dir,
                                                          "ceelo_drf/" + preset + "_response.xml" ) );
    BOOST_REQUIRE_MESSAGE( !xml.empty(), "missing golden " << preset );
    return ceelo::DetectorResponse::from_xml_string( xml )->descriptor;
  }

  ceelo::GeometryDescriptor nai_log_base()
  {
    // 2" x 4" x 16" NaI log, its long axis along the detector axis: side incidence sees the large
    //  faces, and phi matters.
    ceelo::GeometryDescriptor gd;
    gd.set_dimensions( ceelo::BoxDims{ 2.54, 5.08, 40.64 } );
    gd.crystal_material_index = add_material( gd, ceelo::make_NaI() );
    ceelo::LayerSpec can;
    can.material_index = add_material( gd, ceelo::make_Aluminum() );
    can.front_thickness_cm = 0.1;
    can.side_thickness_cm = 0.1;
    can.z_end_cm = 40.64;
    gd.layers.push_back( can );
    gd.symmetry = ceelo::ResponseSymmetry::Quadrant;
    gd.reference_point = ceelo::ReferencePoint::EndcapFront;
    return gd;
  }


  /** The innermost vacuum layer (an HPGe's gap), or none. */
  int innermost_vacuum_layer( const ceelo::GeometryDescriptor &gd )
  {
    for( size_t i = 0; i < gd.layers.size(); ++i )
    {
      const int m = gd.layers[i].material_index;
      if( (m >= 0) && (static_cast<size_t>(m) < gd.materials.size())
          && SpecUtils::icontains( gd.materials[m].name, "vacuum" ) )
        return static_cast<int>( i );
    }
    return -1;
  }


  /** HPGe: a vacuum gap behind the crystal, the endcap's back plate, and a copper disk standing in
   for the cold finger / mount - about what a .par record's back slots describe.  The detailed
   variant adds a mount cup, a thinner copper disk and a steel cryostat/preamp block. */
  ceelo::GeometryDescriptor hpge_with_back( ceelo::GeometryDescriptor gd, const bool detailed )
  {
    int gap = innermost_vacuum_layer( gd );
    if( gap < 0 )
    {
      insert_back_layer( gd, 0, add_vacuum( gd ), 1.0 );
      gap = 0;
    }else
    {
      gd.set_layer_back_thickness_cm( static_cast<size_t>(gap), 1.0 );
    }

    gd.set_layer_back_thickness_cm( gd.layers.size() - 1, 0.15 );   // endcap back plate
    const int cu = add_material( gd, ceelo::make_Copper() );
    if( !detailed )
    {
      append_back_layer( gd, cu, 0.5 );
    }else
    {
      // An aluminium holder cup just inside the gap, a 2 mm copper plate, and 2 cm of steel.
      ceelo::LayerSpec cup;
      cup.material_index = add_material( gd, ceelo::make_Aluminum() );
      cup.front_thickness_cm = 0.0;
      cup.side_thickness_cm = 0.08;
      cup.z_start_cm = 0.0;
      cup.z_end_cm = gd.dimensions_cm.back();
      gd.layers.insert( gd.layers.begin() + gap, cup );
      gd.set_layer_back_thickness_cm( static_cast<size_t>(gap), 0.3 );
      append_back_layer( gd, cu, 0.2 );
      append_back_layer( gd, add_material( gd, ceelo::make_StainlessSteel304() ), 2.0 );
    }
    return gd;
  }


  /** NaI: a glass window behind the crystal and a steel PMT shell.  The detailed variant has the
   PMT's evacuated interior, its back glass, a mu-metal shield and an aluminium housing. */
  ceelo::GeometryDescriptor nai_with_back( ceelo::GeometryDescriptor gd, const bool detailed )
  {
    const int glass = add_material( gd, ceelo::make_Pyrex() );
    insert_back_layer( gd, 0, glass, 0.6 );   // optical window, inside the can
    const int steel = add_material( gd, ceelo::make_StainlessSteel304() );
    if( !detailed )
    {
      append_back_layer( gd, steel, 0.3 );
    }else
    {
      append_back_layer( gd, add_vacuum( gd ), 4.0 );
      append_back_layer( gd, glass, 0.2 );
      append_back_layer( gd, steel, 0.1 );
      append_back_layer( gd, add_material( gd, ceelo::make_Aluminum() ), 0.15 );
    }
    return gd;
  }


  struct BackDet
  {
    string name;
    ceelo::GeometryDescriptor approx;     ///< what the responses are built from
    ceelo::GeometryDescriptor detailed;   ///< only to measure what the approximation costs
    vector<double> phis_deg;
  };


  vector<BackDet> back_detectors()
  {
    vector<BackDet> dets;
    {
      const ceelo::GeometryDescriptor base = gem35_70_base();
      dets.push_back( { "gem35_70", hpge_with_back( base, false ), hpge_with_back( base, true ), { 0.0 } } );
    }
    {
      const ceelo::GeometryDescriptor base = golden_descriptor( "detective_x" );
      dets.push_back( { "detective_x", hpge_with_back( base, false ), hpge_with_back( base, true ), { 0.0 } } );
    }
    {
      const ceelo::GeometryDescriptor base = golden_descriptor( "nai3x3" );
      dets.push_back( { "nai3x3", nai_with_back( base, false ), nai_with_back( base, true ), { 0.0 } } );
    }
    {
      const ceelo::GeometryDescriptor base = nai_log_base();
      dets.push_back( { "nai_log", nai_with_back( base, false ), nai_with_back( base, true ), { 0.0, 45.0 } } );
    }

    for( BackDet &d : dets )
    {
      BOOST_REQUIRE_MESSAGE( d.approx.problems().empty(), d.name << " approx descriptor has problems" );
      BOOST_REQUIRE_MESSAGE( d.detailed.problems().empty(), d.name << " detailed descriptor has problems" );
      BOOST_REQUIRE_CLOSE( d.approx.endcap_front_offset_cm(), d.detailed.endcap_front_offset_cm(), 1e-9 );
    }
    return dets;
  }


  // ---------------------------------------------------------------------------------------------
  //  Points and Monte Carlo
  // ---------------------------------------------------------------------------------------------

  const vector<double> sm_energies = { 40.0, 60.0, 122.0, 344.0, 662.0, 1332.0 };
  const vector<double> sm_thetas = { 92.0, 100.0, 110.0, 125.0, 145.0, 165.0, 180.0 };
  const vector<double> sm_clearances = { 0.05, 0.3, 1.0, 2.0, 4.0, 9.0 };  // in units of a

  /** Crystal-frame position at (theta, phi) from the crystal-face origin, `s_over_a` past the
   detector's outside.  cos(theta) = -z/d, as CeeLo measures it. */
  Eigen::Vector3d clearance_position( const ceelo::GeometryDescriptor &gd, const double theta_deg,
                                      const double phi_deg, const double s_over_a )
  {
    const double th = theta_deg * M_PI / 180.0, ph = phi_deg * M_PI / 180.0;
    const Eigen::Vector3d u( std::sin(th)*std::cos(ph), std::sin(th)*std::sin(ph), -std::cos(th) );
    const double d = gd.outer_envelope().exit_distance_cm( u ) + s_over_a * gd.transverse_half_extent();
    return d * u;
  }


  struct McPoint
  {
    double fep = 0.0, fep_sig = 0.0, tot = 0.0, tot_sig = 0.0;
    uint64_t events = 0;
    double cpu_s = 0.0;

    McPoint() = default;
    explicit McPoint( const ceelo::EfficiencyResult &r )
      : fep( r.full_energy_peak_efficiency ), fep_sig( r.fep_uncertainty ), tot( r.total_efficiency ),
        tot_sig( r.total_uncertainty ), events( r.num_events_simulated ), cpu_s( r.cpu_time_seconds )
    {
    }
  };

  McPoint run_mc( const ceelo::GeometryDescriptor &gd, const Eigen::Vector3d &src, const double energy,
                  const uint64_t seed )
  {
    vector<unique_ptr<ceelo::Material>> owned;
    ceelo::EfficiencyCalculator calc;
    ceelo::ResponseGenerator::configure_calculator( calc, gd, owned );
    calc.enable_general_cone_sampling( true );
    calc.set_point_source( src );

    ceelo::SimulationConfig cfg;
    cfg.energy_keV = energy;
    cfg.termination.target_fep_rel_precision = 0.003;
    cfg.termination.max_events = 40000000;
    cfg.termination.min_events = 20000;
    cfg.termination.max_wall_seconds = 120.0;
    cfg.num_threads = g_threads;
    cfg.batch_size = 20000;
    cfg.seed = seed;
    return McPoint( calc.compute( cfg ) );
  }


  /** A uniform water source beside or behind the detector, as an InterSpec assembly: geometry,
   outer dims (sphere r; cylinder r, half-length; box half width, height, depth), the centre's
   distance to the face plane along the axis, and a lateral offset. */
  struct VolScene
  {
    const char *name;
    GammaInteractionCalc::GeometryType geom;
    std::array<double,3> dims_cm;
    double distance_cm, offset_cm;
  };

  const vector<VolScene> &vol_scenes()
  {
    using GT = GammaInteractionCalc::GeometryType;
    static const vector<VolScene> scenes = {
      { "endon_straddle", GT::CylinderEndOn,  { 2.0, 5.0, 0.0 },   0.0,  9.0 },
      { "sphere_beside",  GT::Spherical,      { 3.0, 0.0, 0.0 },  -4.0, 10.0 },
      { "box_beside",     GT::Rectangular,    { 3.0, 3.0, 4.0 },  -2.0,  9.0 },
      { "sideon_beside",  GT::CylinderSideOn, { 2.0, 5.0, 0.0 },   0.0,  9.0 },
      { "sphere_behind",  GT::Spherical,      { 2.0, 0.0, 0.0 }, -20.0,  0.0 },
    };
    return scenes;
  }

  const vector<double> sm_vol_energies = { 60.0, 122.0, 661.7, 1332.5 };

  bool is_vol_detector( const string &name )
  {
    return (name == "gem35_70") || (name == "nai3x3");
  }

  double vol_scene_volume_cm3( const VolScene &sc )
  {
    const double pi = 3.14159265358979323846;
    switch( sc.geom )
    {
      case GammaInteractionCalc::GeometryType::Spherical:
        return 4.0/3.0 * pi * std::pow( sc.dims_cm[0], 3.0 );
      case GammaInteractionCalc::GeometryType::CylinderEndOn:
      case GammaInteractionCalc::GeometryType::CylinderSideOn:
        return pi * sc.dims_cm[0] * sc.dims_cm[0] * 2.0 * sc.dims_cm[1];
      case GammaInteractionCalc::GeometryType::Rectangular:
        return 8.0 * sc.dims_cm[0] * sc.dims_cm[1] * sc.dims_cm[2];
      case GammaInteractionCalc::GeometryType::NumGeometryType:
        break;
    }
    return 0.0;
  }

  shared_ptr<const Material> water_material()
  {
    BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
    const shared_ptr<const Material> water = MaterialDB::instance()->material( "water" );
    BOOST_REQUIRE( water );
    return water;
  }

  /** CeeLo MC of scene `sc` around `gd`: the mirror of detector_geom_from_config, which puts the
   detector at (-dx, 0, D) from the assembly (side-on: at (D, -dx, 0), its axis along +x); here the
   source centre goes to (dx, 0, -(off + D)) from the crystal face (side-on: (0, dx, ...), with the
   cylinder axis along x - perpendicular to both the detector axis and the offset, as there). */
  McPoint run_volumetric_mc( const ceelo::GeometryDescriptor &gd, const VolScene &sc, const double energy,
                             const uint64_t seed )
  {
    using GT = GammaInteractionCalc::GeometryType;
    vector<unique_ptr<ceelo::Material>> owned;
    ceelo::EfficiencyCalculator calc;
    ceelo::ResponseGenerator::configure_calculator( calc, gd, owned );

    const double z = -(gd.endcap_front_offset_cm() + sc.distance_cm);
    switch( sc.geom )
    {
      case GT::Spherical:
        calc.set_spherical_source( Eigen::Vector3d( sc.offset_cm, 0.0, z ), sc.dims_cm[0] );
        break;
      case GT::CylinderEndOn:
        calc.set_cylindrical_source( Eigen::Vector3d( sc.offset_cm, 0.0, z ), sc.dims_cm[0], sc.dims_cm[1] );
        break;
      case GT::CylinderSideOn:
        calc.set_cylindrical_source( Eigen::Vector3d( 0.0, sc.offset_cm, z ), sc.dims_cm[0], sc.dims_cm[1],
                                     Eigen::Matrix3d( Eigen::AngleAxisd( -0.5*3.14159265358979323846,
                                                                         Eigen::Vector3d::UnitY() ) ) );
        break;
      case GT::Rectangular:
        calc.set_rectangular_source( Eigen::Vector3d( sc.offset_cm, 0.0, z ),
                                     Eigen::Vector3d( sc.dims_cm[0], sc.dims_cm[1], sc.dims_cm[2] ) );
        break;
      case GT::NumGeometryType:
        BOOST_REQUIRE( false );
    }

    const ceelo::Material matrix = CeeLoUtils::to_ceelo_material( *water_material() ).to_material();
    calc.set_source_material( &matrix );

    ceelo::SimulationConfig cfg;
    cfg.energy_keV = energy;
    cfg.termination.target_fep_rel_precision = 0.005;
    cfg.termination.max_events = 2000000000ULL;
    cfg.termination.min_events = 100000;
    cfg.termination.max_wall_seconds = 300.0;
    cfg.num_threads = g_threads;
    cfg.seed = seed;
    return McPoint( calc.compute( cfg ) );
  }//run_volumetric_mc(...)


  /** The water matrix's removal coefficient as the MC scores the peak: mu_total less the Compton
   scatters that stay inside CeeLo's FEP window (VolumetricNearFieldHarness's fep_removal_coefficient,
   for a source matrix). */
  double water_fep_removal( const Material &mat, const double energy_keV )
  {
    const double mu_total = GammaInteractionCalc::transmition_length_coefficient( &mat, energy_keV );
    double mu_compton = 0.0;
    for( const Material::ElementFractionPair &p : mat.elements )
      mu_compton += p.second * mat.density
                    * MassAttenuation::massAttenuationCoefficientElement( p.first->atomicNumber, energy_keV,
                                                       MassAttenuation::GammaEmProcces::ComptonScatter );
    const ceelo::Material cmat = CeeLoUtils::to_ceelo_material( mat ).to_material();
    return mu_total - ceelo::kn_in_window_fraction( energy_keV, ceelo::kDefaultFepWindowKeV, cmat ) * mu_compton;
  }//water_fep_removal(...)


  /** InterSpec's line path for scene `sc` with response `resp`: FEP efficiency per emitted photon,
   and the worst response flag the source's chords met. */
  double line_path_eff( const shared_ptr<const ceelo::DetectorResponse> &resp, const VolScene &sc,
                        const double energy, ceelo::ResponseFlag &worst_flag )
  {
    using namespace GammaInteractionCalc;
    const double cm = PhysicalUnits::cm;
    const shared_ptr<const Material> water = water_material();

    DistributedSrcCalcT<double> calc;
    calc.m_geometry = sc.geom;
    calc.m_materialIndex = 0;
    calc.m_attenuateForAir = false;
    calc.m_airTransLenCoef = 0.0;
    calc.m_isInSituExponential = false;
    calc.m_inSituRelaxationLength = -1.0;
    calc.m_srcVolumetricActivity = 1.0;
    calc.m_normalizeByVolume = false;
    calc.m_energy = energy;
    calc.m_nuclide = nullptr;
    calc.integral = 0.0;
    calc.m_effResponse = resp;
    calc.m_effMethod = ShieldingSourceFitCalc::VolumetricEffMethod::MCTransfer;
    calc.m_detector = detector_geom_from_config<double>( sc.geom, sc.distance_cm*cm,
                                                         resp->descriptor.transverse_half_extent()*cm, 0.0,
                                                         sc.offset_cm*cm, 0.0 );
    DistributedSrcCalcT<double>::ShellInfo shell;
    shell.dims = { sc.dims_cm[0]*cm, sc.dims_cm[1]*cm, sc.dims_cm[2]*cm };
    shell.trans_len_coef = water_fep_removal( *water, energy );
    shell.type = ShellType::Material;
    calc.m_shells.push_back( shell );

    const double outer[3] = { shell.dims[0], shell.dims[1], shell.dims[2] };
    BOOST_REQUIRE_GT( detector_enclosure_clearance( resp->descriptor, sc.geom, outer, sc.distance_cm*cm,
                                                    sc.offset_cm*cm, 0.0 ), 0.0 );

    const std::array<double,3> det_pos = { calc.m_detector.position[0], calc.m_detector.position[1], calc.m_detector.position[2] };
    const std::array<double,3> det_axis = { calc.m_detector.axis[0], calc.m_detector.axis[1], calc.m_detector.axis[2] };
    calc.m_lineCache = build_volumetric_line_cache( resp, sc.geom, 0, shell.dims, det_pos, det_axis, 0.0, 1 << 17,
                                                    ShieldingSourceChi2Fcn::sm_default_volumetric_line_pad,
                                                    sm_default_volumetric_line_surface_frac, LineSampleParams(),
                                                    sm_volumetric_line_hemi_frac, true );
    std::vector<std::unique_ptr<DistributedSrcCalcT<double>>> calcs;
    calcs.push_back( std::make_unique<DistributedSrcCalcT<double>>( calc ) );
    {
      const ScopedVolumetricIntegratorOverride force( VolumetricIntegrator::Line );
      integrate_volumetric_calculators<double>( calcs, true );
    }
    worst_flag = calcs.front()->m_lineCache->worst_flag( energy );
    return calcs.front()->integral / (vol_scene_volume_cm3( sc ) * cm*cm*cm);
  }//line_path_eff(...)


  struct TruthRow
  {
    string region;      ///< "anchor", or "point"
    double energy = 0.0, theta = 0.0, phi = 0.0, s_over_a = 0.0;
    Eigen::Vector3d src = Eigen::Vector3d::Zero();
    McPoint mc;
  };

  const char * const sm_truth_header =
      "region,E_keV,theta_deg,phi_deg,s_over_a,x_cm,y_cm,z_cm,fep,fep_sig,tot,tot_sig,events,cpu_s";

  vector<TruthRow> read_truth( const string &path )
  {
    vector<TruthRow> rows;
    ifstream in( path.c_str() );
    string line;
    while( std::getline( in, line ) )
    {
      if( line.empty() || (line[0] == '#') || SpecUtils::istarts_with( line, "region" ) )
        continue;
      vector<string> f;
      SpecUtils::split( f, line, "," );
      if( f.size() < 14 )
        continue;
      TruthRow r;
      r.region = f[0];
      r.energy = stod( f[1] );
      r.theta = stod( f[2] );
      r.phi = stod( f[3] );
      r.s_over_a = stod( f[4] );
      r.src = Eigen::Vector3d( stod( f[5] ), stod( f[6] ), stod( f[7] ) );
      r.mc.fep = stod( f[8] );
      r.mc.fep_sig = stod( f[9] );
      r.mc.tot = stod( f[10] );
      r.mc.tot_sig = stod( f[11] );
      r.mc.events = static_cast<uint64_t>( stod( f[12] ) );
      r.mc.cpu_s = stod( f[13] );
      rows.push_back( r );
    }
    return rows;
  }


  /** An EFFTRAN curve transfer anchored on the truth bank's own on-axis anchor rows: the same
   physics as the truth, so the comparison isolates the angular model. */
  shared_ptr<ceelo::DetectorResponse> transfer_from_anchor( const ceelo::GeometryDescriptor &gd,
                                                            const vector<TruthRow> &rows,
                                                            const string &name )
  {
    CeeLoUtils::TransferAnchor anchor;
    anchor.curve_derived = false;
    double ref_z = 0.0;
    for( const TruthRow &r : rows )
    {
      if( r.region != "anchor" )
        continue;
      ref_z = r.src.z();
      anchor.curve.energies_keV.push_back( r.energy );
      anchor.curve.eff.push_back( r.mc.fep );
      anchor.curve.frac_sigma.push_back( r.mc.fep_sig / r.mc.fep );
    }
    if( anchor.curve.energies_keV.size() < 2 )
      return nullptr;
    anchor.ref_distance_cm = -ref_z - gd.endcap_front_offset_cm();   // crystal frame -> face
    return CeeLoUtils::makeTransferResponse( gd, anchor, ceelo::AnchorCurve{}, name , 1.0);
  }
}//namespace

BOOST_GLOBAL_FIXTURE( Fixture );


BOOST_AUTO_TEST_CASE( BackStructureDescriptors )
{
  // The descriptors themselves: back structure present, front untouched.
  for( const BackDet &d : back_detectors() )
  {
    BOOST_CHECK_GT( d.approx.back_structure_extent_cm(), 0.3 );
    BOOST_CHECK_GT( d.detailed.back_structure_extent_cm(), d.approx.back_structure_extent_cm() );
    const ceelo::GeometryDescriptor::Envelope env = d.approx.outer_envelope();
    BOOST_CHECK_GT( env.z_max_cm, d.approx.dimensions_cm.back() );
    BOOST_TEST_MESSAGE( d.name << ": back extent " << d.approx.back_structure_extent_cm()
                        << " cm (detailed " << d.detailed.back_structure_extent_cm() << " cm), envelope z ["
                        << env.z_min_cm << ", " << env.z_max_cm << "]" );
  }
}


BOOST_AUTO_TEST_CASE( GenerateBackPointTruth )
{
  if( g_truth_out.empty() )
  {
    BOOST_TEST_MESSAGE( "GenerateBackPointTruth: skipped (pass --back-truth-out=<dir>)" );
    return;
  }
  if( !SpecUtils::is_directory( g_truth_out ) )
    SpecUtils::create_directory( g_truth_out );

  for( const BackDet &det : back_detectors() )
  {
    if( !g_only.empty() && !SpecUtils::icontains( det.name, g_only ) )
      continue;

    const double a = det.approx.transverse_half_extent();

    // 1. The point bank: anchors on axis in front, then beside/behind in clearance coordinates.
    {
      const string path = SpecUtils::append_path( g_truth_out, det.name + "_points.csv" );
      ofstream out( path.c_str() );
      BOOST_REQUIRE( out.is_open() );
      out << std::setprecision( 10 );
      out << "# Monte-Carlo truth for " << det.name << " (approximate back structure), CeeLo point"
             " sources, general cone sampling, FEP target 0.3%\n";
      out << sm_truth_header << "\n";

      uint64_t seed = 1000;
      const double far_z = -std::max( 10.0*a, 10.0 );
      for( const double e : sm_energies )
      {
        const Eigen::Vector3d src( 0.0, 0.0, far_z );
        const McPoint p = run_mc( det.approx, src, e, ++seed );
        out << "anchor," << e << ",0,0," << (-far_z)/a << ",0,0," << far_z << ","
            << p.fep << "," << p.fep_sig << "," << p.tot << "," << p.tot_sig << ","
            << p.events << "," << p.cpu_s << "\n";
      }
      out.flush();

      for( const double phi : det.phis_deg )
      for( const double theta : sm_thetas )
      for( const double s : sm_clearances )
      {
        const Eigen::Vector3d src = clearance_position( det.approx, theta, phi, s );
        for( const double e : sm_energies )
        {
          const McPoint p = run_mc( det.approx, src, e, ++seed );
          out << "point," << e << "," << theta << "," << phi << "," << s << ","
              << src.x() << "," << src.y() << "," << src.z() << ","
              << p.fep << "," << p.fep_sig << "," << p.tot << "," << p.tot_sig << ","
              << p.events << "," << p.cpu_s << "\n";
        }
        out.flush();
        cout << det.name << ": theta " << theta << " phi " << phi << " s/a " << s << " done" << endl;
      }
      out << "#complete\n";
    }

    // 2. What the approximate back costs: the same points with the detailed back, behind the
    //    crystal's back plane where the difference lives.
    {
      const string path = SpecUtils::append_path( g_truth_out, det.name + "_back_approx_vs_detailed.csv" );
      ofstream out( path.c_str() );
      BOOST_REQUIRE( out.is_open() );
      out << std::setprecision( 10 );
      out << "E_keV,theta_deg,phi_deg,s_over_a,approx_fep,approx_sig,detailed_fep,detailed_sig,"
             "approx_tot,detailed_tot\n";
      uint64_t seed = 900000;
      for( const double theta : { 125.0, 145.0, 165.0, 180.0 } )
      for( const double s : { 0.3, 1.0, 4.0 } )
      {
        // Clearance is measured on the APPROXIMATE outside, so both see the same source point.
        const Eigen::Vector3d src = clearance_position( det.approx, theta, det.phis_deg.front(), s );
        // The detailed back may reach further; skip points that would sit inside it.
        const ceelo::GeometryDescriptor::Envelope denv = det.detailed.outer_envelope();
        const double dz = denv.exit_distance_cm( src.normalized() );
        if( src.norm() <= dz )
          continue;
        for( const double e : sm_energies )
        {
          const McPoint pa = run_mc( det.approx, src, e, ++seed );
          const McPoint pd = run_mc( det.detailed, src, e, ++seed );
          out << e << "," << theta << "," << det.phis_deg.front() << "," << s << ","
              << pa.fep << "," << pa.fep_sig << "," << pd.fep << "," << pd.fep_sig << ","
              << pa.tot << "," << pd.tot << "\n";
        }
        out.flush();
      }
      out << "#complete\n";
    }
  }//for( const BackDet &det : back_detectors() )
}//GenerateBackPointTruth


BOOST_AUTO_TEST_CASE( GenerateBackCoveredResponses )
{
  if( g_truth_out.empty() )
  {
    BOOST_TEST_MESSAGE( "GenerateBackCoveredResponses: skipped (pass --back-truth-out=<dir>)." );
    return;
  }

  for( const BackDet &det : back_detectors() )
  {
    if( !g_only.empty() && !SpecUtils::icontains( det.name, g_only ) )
      continue;
    // The full-MC response the measurements compare: General profile, 0.5% nodes, fixed seed.
    ceelo::GenerationOptions opts;
    opts.detector_name = det.name + " (approximate back)";
    opts.profile = ceelo::ResponseProfile::General;
    opts.node_fep_precision = 0.005;
    opts.num_threads = g_threads;
    opts.base_seed = 4242;
    opts.back_coverage = ceelo::BackCoverage::Coarse;
    opts.back_body_depth_cm = sm_body_depth_cm;   // as InterSpec's MC tool
    ceelo::GenerationStats stats;
    opts.stats_out = &stats;
    const ceelo::ResponseGenerator::NodePlan plan = ceelo::ResponseGenerator::plan_nodes( det.approx, opts );
    const shared_ptr<ceelo::DetectorResponse> resp = ceelo::ResponseGenerator::generate( det.approx, opts );
    BOOST_REQUIRE( resp );
    ofstream out( SpecUtils::append_path( g_truth_out, det.name + "_fullmc_coarse_response.xml" ).c_str() );
    BOOST_REQUIRE( out.is_open() );
    out << resp->to_xml_string();

    // What each stage cost, CPU seconds and MC nodes run (stage 5 skips refuse-grade nodes).
    size_t n_by_stage[8] = {};
    for( const ceelo::NodeStat &ns : stats.nodes )
      ++n_by_stage[std::min( ns.stage, 7u )];
    cout << det.name << ": covered response written; CPU s by stage:";
    for( unsigned st = 1; st <= 5; ++st )
      cout << " " << st << "=" << std::fixed << std::setprecision( 0 ) << stats.cpu_s_by_stage[st]
           << " (" << n_by_stage[st] << " nodes)";
    cout << std::defaultfloat << std::setprecision( 6 ) << "; stage 5 planned " << plan.n_back_near() << ", body depth "
         << resp->provenance.back_body_depth_cm << " cm; total " << stats.total_cpu_s << " s" << endl;
  }
}//GenerateBackCoveredResponses


BOOST_AUTO_TEST_CASE( GenerateVolumetricSceneTruth )
{
  if( g_truth_out.empty() )
  {
    BOOST_TEST_MESSAGE( "GenerateVolumetricSceneTruth: skipped (pass --back-truth-out=<dir>)." );
    return;
  }

  for( const BackDet &det : back_detectors() )
  {
    if( !is_vol_detector( det.name ) || (!g_only.empty() && !SpecUtils::icontains( det.name, g_only )) )
      continue;
    const string path = SpecUtils::append_path( g_truth_out, det.name + "_volumetric.csv" );
    ofstream out( path.c_str() );
    BOOST_REQUIRE( out.is_open() );
    out << "# Monte-Carlo truth for " << det.name << " (approximate back): uniform water sources beside"
           " and behind it; CeeLo FEP per emitted photon\n";
    out << "scene,E_keV,fep,fep_sig,events,cpu_s\n" << std::setprecision( 10 );
    uint64_t seed = 700000;
    for( const VolScene &sc : vol_scenes() )
    {
      for( const double e : sm_vol_energies )
      {
        const McPoint p = run_volumetric_mc( det.approx, sc, e, ++seed );
        out << sc.name << "," << e << "," << p.fep << "," << p.fep_sig << "," << p.events << "," << p.cpu_s << "\n";
        out.flush();
        cout << det.name << ": " << sc.name << " @ " << e << " keV: " << p.fep << " +- "
             << 100.0*p.fep_sig/std::max( p.fep, 1e-300 ) << "%" << endl;
      }
    }
    out << "#complete\n";
  }//for( detectors )
}//GenerateVolumetricSceneTruth


BOOST_AUTO_TEST_CASE( VolumetricSceneMeasure )
{
  const string dir = g_truth_out.empty() ? truth_dir() : g_truth_out;
  for( const BackDet &det : back_detectors() )
  {
    if( !is_vol_detector( det.name ) || (!g_only.empty() && !SpecUtils::icontains( det.name, g_only )) )
      continue;
    ifstream in( SpecUtils::append_path( dir, det.name + "_volumetric.csv" ).c_str() );
    if( !in.is_open() )
      continue;

    map<pair<string,double>,McPoint> truth;
    string line;
    while( std::getline( in, line ) )
    {
      if( line.empty() || (line[0] == '#') || SpecUtils::starts_with( line, "scene," ) )
        continue;
      vector<string> f;
      SpecUtils::split( f, line, "," );
      if( f.size() < 4 )
        continue;
      McPoint p;
      p.fep = std::stod( f[2] );
      p.fep_sig = std::stod( f[3] );
      truth[{ f[0], std::stod( f[1] ) }] = p;
    }

    for( const auto &kind_resp : full_mc_responses( dir, det.name ) )
    {
      const char * const kind = kind_resp.first;
      const shared_ptr<const ceelo::DetectorResponse> &resp = kind_resp.second;
      for( const VolScene &sc : vol_scenes() )
      {
        for( const double e : sm_vol_energies )
        {
          const auto pos = truth.find( { sc.name, e } );
          if( (pos == end(truth)) || !(pos->second.fep > 0.0) )
            continue;
          ceelo::ResponseFlag flag = ceelo::ResponseFlag::Ok;
          const double model = line_path_eff( resp, sc, e, flag );
          const McPoint &mc = pos->second;
          const double rel = model/mc.fep - 1.0;
          cout << det.name << " | " << kind << " | " << sc.name << " @ " << e << " keV: line/MC - 1 = "
               << std::showpos << std::fixed << std::setprecision( 2 ) << 100.0*rel
               << "%" << std::noshowpos << " (MC +- " << 100.0*mc.fep_sig/mc.fep << "%) "
               << ceelo::to_string( flag ) << std::defaultfloat << endl;

          // Gates (2026-11-05) for a covered response with side & back near-field tables: the
          //  scene behind within max(1.5%, 3 sigma of the MC), every scene within 2.7% - where the
          //  response gives an answer at all.  A refuse-grade one (NeedsMc: the direct path through
          //  the back nearly blocked, 60 keV behind) is left out, as in BackEnvelopeMeasure.
          if( (string( kind ) == "fullmc_coarse") && !resp->back_near_field.empty()
              && (flag != ceelo::ResponseFlag::NeedsMc) )
          {
            const double tol = (sc.name == "sphere_behind") ? std::max( 0.015, 3.0*mc.fep_sig/mc.fep ) : 0.027;
            BOOST_CHECK_MESSAGE( std::fabs( rel ) <= tol, det.name << " " << sc.name << " @ " << e << " keV: line/MC - 1 = "
                                 << 100.0*rel << "% exceeds " << 100.0*tol << "%" );
          }
        }
      }
    }//for( response kinds )
  }//for( detectors )
}//VolumetricSceneMeasure


BOOST_AUTO_TEST_CASE( PrefactorGridBehindResolution )
{
  // The line path reads the response prefactor off a grid in (ln d, cos theta) from the crystal-face
  //  origin (12 nodes per octave); close behind a detector a ln d step there is a larger step in the
  //  side & back near-field table's own distance (clearance past the housing).  Gate (2026-11-05):
  //  grid vs direct evaluation within 0.3%, for water spheres close beside and behind each detector.
  using GT = GammaInteractionCalc::GeometryType;
  const string dir = g_truth_out.empty() ? truth_dir() : g_truth_out;
  for( const BackDet &det : back_detectors() )
  {
    if( !g_only.empty() && !SpecUtils::icontains( det.name, g_only ) )
      continue;
    const string xml = read_file( SpecUtils::append_path( dir, det.name + "_fullmc_coarse_response.xml" ) );
    if( xml.empty() )
      continue;
    const shared_ptr<const ceelo::DetectorResponse> resp = ceelo::DetectorResponse::from_xml_string( xml );
    if( resp->back_near_field.empty() || (resp->descriptor.shape != ceelo::DetectorShape::Cylinder) )
      continue;

    const ceelo::GeometryDescriptor::Envelope env = resp->descriptor.outer_envelope();
    const double r = 1.5, gap = 0.5, side = env.radius_cm + r + gap;
    const vector<VolScene> scenes = {
      { "grid_beside_front", GT::Spherical, { r, 0.0, 0.0 }, -1.0, side },
      { "grid_beside_back",  GT::Spherical, { r, 0.0, 0.0 }, -0.7*env.z_max_cm, side },
      { "grid_behind_near",  GT::Spherical, { r, 0.0, 0.0 }, -(env.z_max_cm + r + gap), 0.0 },
      { "grid_behind_body",  GT::Spherical, { r, 0.0, 0.0 }, -(std::max( env.z_max_cm, sm_body_depth_cm ) + r + gap), 0.0 }
    };
    for( const VolScene &sc : scenes )
    {
      for( const double e : { 60.0, 661.7 } )
      {
        ceelo::ResponseFlag flag = ceelo::ResponseFlag::Ok;
        const double grid = line_path_eff( resp, sc, e, flag );
        GammaInteractionCalc::sm_prefactor_direct_eval = true;
        double direct = 0.0;
        try
        {
          direct = line_path_eff( resp, sc, e, flag );
        }catch( ... )
        {
          GammaInteractionCalc::sm_prefactor_direct_eval = false;
          throw;
        }
        GammaInteractionCalc::sm_prefactor_direct_eval = false;
        const double rel = grid / direct - 1.0;
        cout << det.name << " | " << sc.name << " @ " << e << " keV: grid/direct - 1 = " << std::showpos
             << 100.0*rel << "%" << std::noshowpos << endl;
        BOOST_CHECK_MESSAGE( std::fabs( rel ) <= 0.003, det.name << " " << sc.name << " @ " << e
                             << " keV: prefactor grid off by " << 100.0*rel << "%" );
      }
    }
  }//for( detectors )
}//PrefactorGridBehindResolution


BOOST_AUTO_TEST_CASE( BackEnvelopeMeasure )
{
  const string dir = g_truth_out.empty() ? truth_dir() : g_truth_out;
  size_t n_detectors = 0;

  for( const BackDet &det : back_detectors() )
  {
    if( !g_only.empty() && !SpecUtils::icontains( det.name, g_only ) )
      continue;
    const vector<TruthRow> rows = read_truth( SpecUtils::append_path( dir, det.name + "_points.csv" ) );
    if( rows.empty() )
      continue;
    ++n_detectors;

    vector<pair<string,shared_ptr<const ceelo::DetectorResponse>>> responses;
    if( const shared_ptr<ceelo::DetectorResponse> t = transfer_from_anchor( det.approx, rows, det.name + " transfer" ) )
      responses.emplace_back( "transfer", t );
    for( const auto &kind_resp : full_mc_responses( dir, det.name ) )
      responses.emplace_back( kind_resp.first, kind_resp.second );

    const double crystal_back = det.approx.dimensions_cm.back();
    for( const auto &named : responses )
    {
      // region -> (|rel err|, pull) for FEP and total
      map<string,vector<pair<double,double>>> fep, tot;
      map<string,size_t> flagged;
      // region -> the worst FEP rows, for a look at what drives the tail
      map<string,vector<pair<double,string>>> worst_rows;
      // Every FEP comparison, with the kernel's direct transmission s (transmitted / geometric solid
      //  angle of the active crystal) - what the error behind the detector is expected to track.
      ofstream csv;
      if( !g_measure_out.empty() )
      {
        csv.open( SpecUtils::append_path( g_measure_out, det.name + "_" + named.first + "_measure.csv" ).c_str() );
        csv << "region,E_keV,theta_deg,phi_deg,s_over_a,mc_fep,mc_rel_sig,model_fep,model_rel_sig,rel_err,"
               "transmission,flag,in_body\n";
      }

      // The pre-registered gates (2026-11-05), for a covered response with side & back near-field
      //  tables on a gated detector, outside the assumed body, the direct path mostly open
      //  (s >= 0.3) and not refuse-grade.  Pulls leave out back_structure_approx: it covers a back
      //  described differently from the real one, which MC of the same description cannot see.
      const shared_ptr<const ceelo::DetectorResponse> &resp = named.second;
      const bool gated = is_gated_detector( det.name ) && (named.first == "fullmc_coarse")
                         && !resp->back_near_field.empty();
      const shared_ptr<const ceelo::DetectorResponse> stripped = gated ? without_back_near( *resp ) : nullptr;
      map<string,vector<double>> gate_err, gate_contact_err, gate_pull;   // side -> values
      size_t n_in_body = 0, n_in_body_moved = 0;
      for( const TruthRow &r : rows )
      {
        if( (r.region != "point") || !(r.mc.fep > 0.0) )
          continue;
        const string side = (r.src.z() > crystal_back) ? "behind" : "beside";
        const string where = side + ((r.s_over_a < 1.0) ? "_near" : "_far");
        const ceelo::EffResult m = named.second->eps_fep_at( r.energy, r.src );
        const bool in_body = inside_assumed_body( det.approx, r.src );
        if( gated )
        {
          const ceelo::ApertureQuadrature q = resp->make_quadrature( r.src );
          const double s = (q.omega_frac_active > 0.0) ? resp->kernel_transmitted( r.energy, q ) / q.omega_frac_active : 0.0;
          if( in_body )
          {
            // Inside the body the tables do not apply: the answer is the one without them.
            const ceelo::EffResult o = stripped->eps_fep_at( r.energy, r.src );
            ++n_in_body;
            n_in_body_moved += ((o.value != m.value) || (o.sigma != m.sigma) || (o.flag != m.flag));
          }else if( (s >= 0.3) && (m.flag != ceelo::ResponseFlag::NeedsMc) && (m.value > 0.0) )
          {
            const double rel = m.value / r.mc.fep - 1.0;
            double sig2 = m.sigma*m.sigma;
            if( (r.src.z() > crystal_back) && (det.approx.back_structure_extent_cm() > 0.0) )
              sig2 -= std::pow( ceelo::model_sigma::back_structure_approx * m.value, 2.0 );
            gate_pull[where].push_back( rel * r.mc.fep / std::sqrt( std::max( 0.0, sig2 ) + r.mc.fep_sig*r.mc.fep_sig ) );
            if( r.s_over_a >= 0.3 - 1e-9 )
              gate_err[side].push_back( std::fabs( rel ) );
            else
              gate_contact_err[side].push_back( std::fabs( rel ) );
          }
        }//if( gated )
        if( m.flag == ceelo::ResponseFlag::NeedsMc )
          ++flagged[where];
        if( m.value > 0.0 )
        {
          const double rel = m.value / r.mc.fep - 1.0;
          const double sig = std::sqrt( m.sigma*m.sigma + r.mc.fep_sig*r.mc.fep_sig ) / r.mc.fep;
          // Refuse-grade (NeedsMc) answers are reported by count, not pooled into the pull.
          if( m.flag != ceelo::ResponseFlag::NeedsMc )
            fep[where].emplace_back( std::fabs(rel), rel / sig );
          ostringstream row;
          row << "E " << r.energy << " theta " << r.theta << " phi " << r.phi << " s/a " << r.s_over_a
              << ": MC " << r.mc.fep << " +- " << 100.0*r.mc.fep_sig/r.mc.fep << "%, model " << m.value
              << " (" << std::showpos << 100.0*rel << std::noshowpos << "%, sigma " << 100.0*m.sigma/m.value
              << "%) " << ceelo::to_string( m.flag );
          worst_rows[where].emplace_back( std::fabs(rel), row.str() );
          if( csv.is_open() )
          {
            const ceelo::ApertureQuadrature q = named.second->make_quadrature( r.src );
            const double s = (q.omega_frac_active > 0.0)
                               ? named.second->kernel_transmitted( r.energy, q ) / q.omega_frac_active : 0.0;
            csv << where << "," << r.energy << "," << r.theta << "," << r.phi << "," << r.s_over_a << ","
                << r.mc.fep << "," << r.mc.fep_sig/r.mc.fep << "," << m.value << "," << m.sigma/m.value << ","
                << rel << "," << s << "," << ceelo::to_string( m.flag ) << "," << in_body << "\n";
          }
        }
        try
        {
          const ceelo::EffResult t = named.second->eps_total_at( r.energy, r.src );
          if( (t.value > 0.0) && (r.mc.tot > 0.0) )
          {
            const double rel = t.value / r.mc.tot - 1.0;
            const double sig = std::sqrt( t.sigma*t.sigma + r.mc.tot_sig*r.mc.tot_sig ) / r.mc.tot;
            tot[where].emplace_back( std::fabs(rel), rel / sig );
          }
        }catch( std::exception & )
        {
        }
      }//for( rows )

      const auto report = [&]( const char *what, map<string,vector<pair<double,double>>> &m ){
        for( auto &kv : m )
        {
          vector<double> errs;
          double pull2 = 0.0;
          for( const auto &ep : kv.second )
          {
            errs.push_back( ep.first );
            pull2 += ep.second * ep.second;
          }
          std::sort( begin(errs), end(errs) );
          const size_t n = errs.size();
          const double median = errs[n/2];
          const double p95 = errs[std::min( n - 1, static_cast<size_t>( 0.95 * n ) )];
          const double rms_pull = std::sqrt( pull2 / n );
          cout << det.name << " | " << named.first << " | " << what << " | " << kv.first
               << ": n=" << n << " median " << 100*median << "% p95 " << 100*p95
               << "% max " << 100*errs.back() << "% rms pull " << rms_pull
               << " NeedsMc " << flagged[kv.first] << endl;
        }
      };
      report( "FEP (excl. NeedsMc)", fep );
      report( "total", tot );

      if( gated )
      {
        const auto pct = []( vector<double> v, const double f ) {
          std::sort( begin(v), end(v) );
          return v.empty() ? 0.0 : v[std::min( v.size() - 1, static_cast<size_t>( f * v.size() ) )];
        };
        for( const auto &kv : gate_err )
        {
          const double median = pct( kv.second, 0.5 ), p95 = pct( kv.second, 0.95 );
          cout << "GATE " << det.name << " " << kv.first << " s/a>=0.3: n=" << kv.second.size() << " median "
               << 100*median << "% p95 " << 100*p95 << "%" << endl;
          BOOST_CHECK_MESSAGE( median <= 0.01, det.name << " " << kv.first << ": median " << 100*median << "% > 1%" );
          BOOST_CHECK_MESSAGE( p95 <= 0.03, det.name << " " << kv.first << ": p95 " << 100*p95 << "% > 3%" );
        }
        for( const auto &kv : gate_contact_err )
        {
          const double median = pct( kv.second, 0.5 );
          cout << "GATE " << det.name << " " << kv.first << " s/a=0.05: n=" << kv.second.size() << " median "
               << 100*median << "% p95 " << 100*pct( kv.second, 0.95 ) << "%" << endl;
          BOOST_CHECK_MESSAGE( median <= 0.02, det.name << " " << kv.first << " at contact: median " << 100*median << "% > 2%" );
        }
        for( const auto &kv : gate_pull )
        {
          double pull2 = 0.0;
          for( const double p : kv.second )
            pull2 += p*p;
          const double rms = std::sqrt( pull2 / std::max<size_t>( 1, kv.second.size() ) );
          cout << "GATE " << det.name << " " << kv.first << ": rms pull " << rms << " (n=" << kv.second.size() << ")" << endl;
          BOOST_CHECK_MESSAGE( (rms >= 0.5) && (rms <= 1.3), det.name << " " << kv.first << ": rms pull " << rms << " outside [0.5, 1.3]" );
        }
        cout << "GATE " << det.name << " inside the assumed body: " << n_in_body << " rows, "
             << n_in_body_moved << " differ from the response without the tables" << endl;
        BOOST_CHECK_EQUAL( n_in_body_moved, 0u );
      }//if( gated )
      for( auto &kv : worst_rows )
      {
        std::sort( begin(kv.second), end(kv.second),
                   []( const pair<double,string> &a, const pair<double,string> &b ){ return a.first > b.first; } );
        for( size_t i = 0; i < std::min( kv.second.size(), size_t(5) ); ++i )
          cout << "    worst " << named.first << " " << kv.first << ": " << kv.second[i].second << endl;
      }
    }//for( responses )
  }//for( detectors )

  if( !n_detectors )
    BOOST_TEST_MESSAGE( "BackEnvelopeMeasure: no truth banks in " << dir );
}//BackEnvelopeMeasure

