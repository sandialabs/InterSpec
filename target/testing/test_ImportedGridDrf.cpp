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

/*
 A DRF imported from a .par efficiency grid + DETECTOR.txt, as a user carries it on from the import:
 keeping the grid's full-energy-peak (FEP) efficiency while giving it a Monte-Carlo total efficiency
 (DetEffG2kPar::attachTotalEfficiency), editing its geometry away from the imported one
 (DetectorPeakResponse::monteCarloGeometry), and the consumers that must then treat it consistently.

 Driven by the mock pair test_data/det_eff/mock_DETECTOR.txt + mock.par - NOT vendor data; see
 WriteMockParFiles in test_DetectorEffG2kPar.cpp for how, and from what, they were made.  Totals here
 come from synthetic donors rather than a Monte Carlo, which would take minutes per run; what is
 tested is that the grid carries a donor's total faithfully, not the Monte Carlo itself.
 */

// Must be defined before Windows.h (or any header that includes it) is included
#ifdef _WIN32
  #define WIN32_LEAN_AND_MEAN
  #include <winsock2.h>
  #include <windows.h>
#endif

#define BOOST_TEST_MODULE test_ImportedGridDrf
#include <boost/test/included/unit_test.hpp>

#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

#include <rapidxml/rapidxml.hpp>
#include <rapidxml/rapidxml_print.hpp>

#include <Eigen/Core>

#include "io/ResponseKernel.h"
#include "io/DetectorResponse.h"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectorEffG2kPar.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

#include "ParTestUtils.h"

using namespace std;

namespace
{
  string g_data_dir;
  string g_test_file_dir;

  const double pi = 3.14159265358979323846;

  void set_dirs()
  {
    static bool s_done = false;
    if( s_done )
      return;

    const int argc = boost::unit_test::framework::master_test_suite().argc;
    char ** const argv = boost::unit_test::framework::master_test_suite().argv;
    for( int i = 1; i < argc; ++i )
    {
      const string arg = argv[i];
      if( SpecUtils::starts_with( arg, "--datadir=" ) )
        g_data_dir = arg.substr( 10 );
      else if( SpecUtils::starts_with( arg, "--testfiledir=" ) )
        g_test_file_dir = arg.substr( 14 );
    }

    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( SpecUtils::append_path( g_data_dir, "sandia.decay.xml" ) ),
                           "sandia.decay.xml not in '" << g_data_dir << "'; pass --datadir=" );
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( SpecUtils::append_path( g_test_file_dir, "det_eff/mock.par" ) ),
                           "No det_eff/mock.par in '" << g_test_file_dir << "'; pass --testfiledir=" );

    BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( g_data_dir ) );
    DecayDataBaseServer::setDecayXmlFile( SpecUtils::append_path( g_data_dir, "sandia.decay.xml" ) );
    s_done = true;
  }//set_dirs()


  string det_eff_path( const string &name )
  {
    return SpecUtils::append_path( SpecUtils::append_path( g_test_file_dir, "det_eff" ), name );
  }


  DetEffG2kPar::DetectorDef mock_def()
  {
    ifstream txt( det_eff_path( "mock_DETECTOR.txt" ).c_str(), ios::in | ios::binary );
    const vector<DetEffG2kPar::DetectorDef> defs = DetEffG2kPar::parseDetectorTxt( txt );
    BOOST_REQUIRE_EQUAL( defs.size(), 1 );
    return defs.front();
  }


  /** The geometry the mock grid was generated from: the imported one, with its guessed bore and
   bulletizing corrected - what a user edits the geometry to. */
  ceelo::GeometryDescriptor corrected_geometry()
  {
    return ParTestUtils::mock_truth_geometry( mock_def(), det_eff_path( "Angle-detector-only.detx" ) );
  }


  /** The imported mock DRF - built once, since a build takes tens of seconds in a Debug build. */
  shared_ptr<const DetectorPeakResponse> mock_drf()
  {
    static shared_ptr<const DetectorPeakResponse> s_drf;
    if( !s_drf )
    {
      set_dirs();
      BOOST_REQUIRE_NO_THROW( s_drf = DetEffG2kPar::makeDrfFromFiles( det_eff_path( "mock.par" ),
                                                                      det_eff_path( "mock_DETECTOR.txt" ) ) );
      BOOST_REQUIRE( s_drf && s_drf->ceeloResponse() );
    }
    return s_drf;
  }


  /** A response for `gd` that carries only a (synthetic, smooth) EtaTotTable total on the grid's own
   energy x angle nodes - the stand-in for a Monte Carlo of `gd` - and, optionally, a synthetic
   near-field total table like a General-profile run's.  Its FEP table is the grid's, only so it is
   a complete response; nothing here reads it. */
  shared_ptr<ceelo::DetectorResponse> total_donor( const ceelo::DetectorResponse &grid,
                                                   const ceelo::GeometryDescriptor &gd,
                                                   const bool with_near_table = false )
  {
    auto donor = make_shared<ceelo::DetectorResponse>();
    donor->descriptor = gd;
    for( size_t i = 0; i < gd.materials.size(); ++i )
      donor->mu_tables.push_back( ceelo::MuTable::sample( gd.materials[i].to_material(), static_cast<int>(i) ) );
    donor->eta_fep = grid.eta_fep;
    donor->provenance = grid.provenance;
    donor->provenance.method = ceelo::ProductionMethod::FullMc;

    ceelo::EtaTable &eta = donor->tot_eff.eta_tot;
    donor->tot_eff.tier = ceelo::TotEffTier::EtaTotTable;
    eta.energies_keV = grid.eta_fep.energies_keV;
    eta.cos_thetas = grid.eta_fep.cos_thetas;
    eta.edges_keV = grid.eta_fep.edges_keV;
    for( const double energy : eta.energies_keV )
    {
      for( const double cos_theta : eta.cos_thetas )
      {
        // Total above the kernel by a factor rising with energy (Compton deposits), a little more
        //  off axis - the shape a real eta_tot has, which is all this needs.
        eta.ln_eta.push_back( 0.2 + 0.15*std::log( energy / 100.0 ) + 0.05*(1.0 - cos_theta) );
        eta.frac_sigma.push_back( 0.01 );
      }
    }

    if( with_near_table )
    {
      // A few-percent close-in correction, growing toward grazing and contact, fading to 1 at its
      //  anchor - the shape a measured one has.
      ceelo::NearFieldModel &tn = donor->tot_eff.near_field;
      tn.energies_keV = { 60.0, 300.0, 1500.0 };
      tn.cos_thetas = { 0.02, 0.3, 0.7, 1.0 };
      tn.dists_cm = { 1.0, 3.0, 8.0, 15.0, 28.0 };
      for( const double energy : tn.energies_keV )
        for( const double cos_theta : tn.cos_thetas )
          for( const double d : tn.dists_cm )
          {
            const bool anchor = (d == tn.dists_cm.back());
            tn.ln_n.push_back( anchor ? 0.0 : (0.06*(1.0 - cos_theta) + 0.02)*std::exp( -d / 6.0 )
                                              * (1.0 + 0.1*std::log( energy / 300.0 )) );
            tn.frac_sigma.push_back( 0.004 );
          }
      tn.break_cos_thetas = { 0.02, 1.0 };
      tn.break_d_cm.assign( tn.energies_keV.size() * 2, tn.dists_cm.back() );
    }//if( with_near_table )

    donor->finalize();
    return donor;
  }//total_donor(...)


  double rel_diff( const double a, const double b )
  {
    return std::fabs( a - b ) / (std::max)( std::fabs(b), 1.0e-300 );
  }
}//namespace


/** The mock pair imports as a grid: FEP only, the importer's guessed bore, the crystal's recess as
 the setback, and the grid's own values on its lattice. */
BOOST_AUTO_TEST_CASE( MockPairImports )
{
  const shared_ptr<const DetectorPeakResponse> drf = mock_drf();
  BOOST_CHECK( drf->hasImportedGrid() );
  BOOST_CHECK( drf->drfSource() == DetectorPeakResponse::DrfSource::CharacterizationParFile );
  BOOST_CHECK( !drf->hasAnyTotalEfficiencyInfo() );
  BOOST_CHECK( !drf->geometryModifiedFromImport() );

  // DETECTOR.txt cannot state the core: the importer guessed one, which is not the detector's.
  const ceelo::GeometryDescriptor &imported = drf->ceeloResponse()->descriptor;
  const ceelo::GeometryDescriptor corrected = corrected_geometry();
  BOOST_REQUIRE( imported.bore && corrected.bore );
  BOOST_CHECK_CLOSE( imported.bore->radius, 0.4, 1.0e-6 );
  BOOST_CHECK_CLOSE( corrected.bore->radius, 0.5, 1.0e-6 );
  BOOST_CHECK_CLOSE( corrected.bore->depth, 3.5, 1.0e-6 );
  BOOST_CHECK( imported.to_xml_string() != corrected.to_xml_string() );

  // The crystal sits the endcap-front offset behind the face (the GADRAS setback convention).
  BOOST_CHECK_CLOSE( drf->detectorSetback(), imported.endcap_front_offset_cm() * PhysicalUnits::cm, 1.0e-4 );

  // On the grid's own lattice - file energies, rows and columns - the DRF answers the file.
  const DetEffG2kPar::ParEfficiency par_eff( DetEffG2kPar::parseParFile( det_eff_path( "mock.par" ) ) );
  const DetEffG2kPar::ParGrid &grid0 = par_eff.parFile().grids.front();
  for( const double energy : { 60.0, 122.0, 662.0, 1400.0 } )
  {
    for( const int row : { 48, 55, 69 } )      // ~12, 25 and 99 cm
    {
      for( const int col : { 0, 6, 9 } )       // 0, 30 and 45 degrees
      {
        const double d_mm = std::exp( row * grid0.r_step );
        const double theta = col * grid0.theta_step_rad;
        const double expected = par_eff.efficiency( energy, d_mm, theta );
        const double value = drf->fepEfficiencyEval( static_cast<float>(energy), theta, 0.0,
                                                     0.1 * d_mm * PhysicalUnits::cm ).value;
        BOOST_CHECK_MESSAGE( rel_diff( value, expected ) < 0.01,
                             energy << " keV, " << 0.1*d_mm << " cm, " << col*5 << " deg: DRF "
                             << value << " vs grid " << expected );
      }
    }
  }
}//BOOST_AUTO_TEST_CASE( MockPairImports )


/** A total from a donor on the grid's own geometry is carried exactly - same kernel, same nodes -
 everywhere, while the FEP stays the grid's to the bit; stripping it restores the import. */
BOOST_AUTO_TEST_CASE( TotalOnImportedGeometry )
{
  const shared_ptr<const ceelo::DetectorResponse> grid = mock_drf()->ceeloResponse();
  const shared_ptr<ceelo::DetectorResponse> donor = total_donor( *grid, grid->descriptor );

  shared_ptr<ceelo::DetectorResponse> with_total;
  BOOST_REQUIRE_NO_THROW( with_total = DetEffG2kPar::attachTotalEfficiency( *grid, donor.get() ) );
  BOOST_REQUIRE( with_total );
  BOOST_CHECK( DetEffG2kPar::isGridResponse( with_total ) );
  BOOST_CHECK( with_total->tot_eff.characterized() );

  double worst = 0.0;
  for( const double d_cm : { 1.0, 5.0, 25.0, 50.0, 200.0 } )
  {
    for( const double theta_deg : { 0.0, 30.0, 60.0, 85.0 } )
    {
      const Eigen::Vector3d pos = CeeLoUtils::sourcePositionFromFace( grid->descriptor, theta_deg*pi/180.0, 0.0, d_cm );
      for( const double energy : { 45.0, 100.0, 300.0, 1000.0, 2500.0 } )
      {
        worst = (std::max)( worst, rel_diff( with_total->eps_total_at( energy, pos ).value,
                                             donor->eps_total_at( energy, pos ).value ) );
        BOOST_CHECK_EQUAL( with_total->eps_fep_at( energy, pos ).value, grid->eps_fep_at( energy, pos ).value );
      }
    }
  }
  BOOST_TEST_MESSAGE( "Total on the imported geometry: worst |graft/donor - 1| = " << worst );
  BOOST_CHECK_LT( worst, 1.0e-9 );

  shared_ptr<ceelo::DetectorResponse> stripped;
  BOOST_REQUIRE_NO_THROW( stripped = DetEffG2kPar::attachTotalEfficiency( *with_total, nullptr ) );
  BOOST_CHECK( !stripped->tot_eff.characterized() );
  BOOST_CHECK_EQUAL( stripped->content_hash(), grid->content_hash() );

  // A DRF carrying it can do cascade summing; the donor must have a total to give.
  auto drf = make_shared<DetectorPeakResponse>( *mock_drf() );
  drf->setCeeloResponse( with_total );
  BOOST_CHECK( drf->hasImportedGrid() );
  BOOST_CHECK( drf->hasAnyTotalEfficiencyInfo() );
  BOOST_CHECK_THROW( DetEffG2kPar::attachTotalEfficiency( *grid, grid.get() ), std::exception );
}//BOOST_AUTO_TEST_CASE( TotalOnImportedGeometry )


/** A donor on an EDITED geometry - the corrected core and bulletizing - with its own near-field
 total, as a General-profile Monte Carlo of it would have.  The grid reproduces it exactly at the
 far-field pin and at every node of its own distance table, and between nodes to what interpolation
 allows; the gates sit a little above the printed figures. */
BOOST_AUTO_TEST_CASE( TotalOnEditedGeometry )
{
  const shared_ptr<const ceelo::DetectorResponse> grid = mock_drf()->ceeloResponse();
  const ceelo::GeometryDescriptor corrected = corrected_geometry();
  const shared_ptr<ceelo::DetectorResponse> donor = total_donor( *grid, corrected, true );

  shared_ptr<ceelo::DetectorResponse> with_total;
  BOOST_REQUIRE_NO_THROW( with_total = DetEffG2kPar::attachTotalEfficiency( *grid, donor.get() ) );

  // The same physical point for each response: the grid's crystal-frame position, shifted by any
  //  difference in where the two put the endcap face.
  const double shift = corrected.endcap_front_offset_cm() - grid->descriptor.endcap_front_offset_cm();
  const auto donor_pos = [shift]( Eigen::Vector3d pos ) -> Eigen::Vector3d { pos.z() -= shift; return pos; };

  // At the far-field pin (from the crystal face), on axis, at the nodes.
  const double d_pin = (std::max)( 1000.0 * grid->descriptor.transverse_half_extent(), 100.0 );
  for( size_t i = 0; i < grid->eta_fep.energies_keV.size(); i += 7 )
  {
    const Eigen::Vector3d pos = ceelo::source_position( d_pin, 1.0 );
    const double energy = grid->eta_fep.energies_keV[i];
    BOOST_CHECK_LT( rel_diff( with_total->eps_total_at( energy, pos ).value,
                              donor->eps_total_at( energy, donor_pos( pos ) ).value ), 1.0e-6 );
  }

  // At every node of the grid's own near-field total table.
  const ceelo::NearFieldModel &tn = with_total->tot_eff.near_field;
  BOOST_REQUIRE( !tn.empty() );
  double worst_node = 0.0;
  for( size_t c = 0; c < tn.cos_thetas.size(); ++c )
    for( size_t d = 0; (d + 1) < tn.dists_cm.size(); ++d )
    {
      const Eigen::Vector3d pos = ceelo::source_position( tn.dists_cm[d], tn.cos_thetas[c] );
      for( size_t e = 0; e < tn.energies_keV.size(); e += 3 )
      {
        if( tn.frac_sigma[tn.index( e, c, d )] >= 1.0 )
          continue;   //held: the donor had nothing there
        worst_node = (std::max)( worst_node, rel_diff( with_total->eps_total_at( tn.energies_keV[e], pos ).value,
                                     donor->eps_total_at( tn.energies_keV[e], donor_pos( pos ) ).value ) );
      }
    }
  BOOST_TEST_MESSAGE( "Total on the edited geometry, at the near-field nodes: worst |graft/donor - 1| = " << worst_node );
  BOOST_CHECK_LT( worst_node, 1.0e-6 );

  // Between nodes: at the same place relative to the endcap face.
  //  Measured: 0.24, 0.87, 0.41, 0.13, 0.07, 0.02, 0.03 % - interpolation between table nodes.
  const vector<pair<double,double>> gates{ {1.0, 0.005}, {2.0, 0.012}, {5.0, 0.006}, {10.0, 0.003},
                                           {25.0, 0.002}, {100.0, 0.001}, {300.0, 0.001} };
  for( const pair<double,double> &gate : gates )
  {
    const double d_cm = gate.first;
    double worst_d = 0.0;
    for( const double theta_deg : { 0.0, 45.0, 75.0 } )
    {
      const double theta = theta_deg * pi / 180.0;
      const Eigen::Vector3d grid_pos = CeeLoUtils::sourcePositionFromFace( grid->descriptor, theta, 0.0, d_cm );
      const Eigen::Vector3d other_pos = CeeLoUtils::sourcePositionFromFace( corrected, theta, 0.0, d_cm );
      for( const double energy : { 60.0, 122.0, 344.0, 662.0, 1332.0, 2614.0 } )
        worst_d = (std::max)( worst_d, rel_diff( with_total->eps_total_at( energy, grid_pos ).value,
                                                 donor->eps_total_at( energy, other_pos ).value ) );
    }
    BOOST_TEST_MESSAGE( "Total on the edited geometry, " << d_cm << " cm (0-75 deg): worst |graft/donor - 1| = "
                        << worst_d );
    BOOST_CHECK_MESSAGE( worst_d < gate.second, "at " << d_cm << " cm: " << worst_d << " vs gate " << gate.second );
  }
}//BOOST_AUTO_TEST_CASE( TotalOnEditedGeometry )


/** Grounding points and transfer anchors for a grid come from the grid itself, at 50 cm - not from
 its legacy curve, which is sampled tens of metres out. */
BOOST_AUTO_TEST_CASE( GridAnchorsAndGrounding )
{
  const shared_ptr<const DetectorPeakResponse> drf = mock_drf();
  const shared_ptr<const ceelo::DetectorResponse> grid = drf->ceeloResponse();
  const Eigen::Vector3d pos = CeeLoUtils::sourcePositionFromFace( grid->descriptor, 0.0, 0.0, 50.0 );

  const vector<ceelo::GroundingPoint> points = DetEffG2kPar::groundingPoints( *grid );
  BOOST_REQUIRE_GE( points.size(), 24 );
  for( const ceelo::GroundingPoint &p : points )
  {
    BOOST_CHECK_CLOSE( p.distance_cm, 50.0, 1.0e-9 );
    BOOST_CHECK_LT( rel_diff( p.measured_eff, grid->eps_fep_at( p.energy_keV, pos ).value ), 1.0e-12 );
  }

  for( const bool with_cov : { false, true } )
  {
    const CeeLoUtils::TransferAnchor anchor = with_cov
          ? CeeLoUtils::curveAnchorWithCovarianceForDrf( drf, grid->descriptor, -1.0 )
          : CeeLoUtils::transferAnchorForDrf( drf, grid->descriptor, 0.0 );
    BOOST_CHECK_CLOSE( anchor.ref_distance_cm, 50.0, 1.0e-9 );
    BOOST_REQUIRE_GE( anchor.curve.energies_keV.size(), 2 );
    for( size_t i = 0; i < anchor.curve.energies_keV.size(); ++i )
      BOOST_CHECK_LT( rel_diff( anchor.curve.eff[i],
                                grid->eps_fep_at( anchor.curve.energies_keV[i], pos ).value ), 1.0e-12 );
  }

  // An explicit distance is honoured.
  BOOST_CHECK_CLOSE( CeeLoUtils::transferAnchorForDrf( drf, grid->descriptor, 25.0 ).ref_distance_cm, 25.0, 1.0e-9 );

  // A transfer built from a grid DRF with a total keeps it.
  const CeeLoUtils::TransferAnchor anchor = CeeLoUtils::transferAnchorForDrf( drf, grid->descriptor, 0.0 );
  BOOST_CHECK( CeeLoUtils::totalTransferAnchorForDrf( drf, anchor ).energies_keV.empty() );

  auto with_total = make_shared<DetectorPeakResponse>( *drf );
  with_total->setCeeloResponse( DetEffG2kPar::attachTotalEfficiency( *grid, total_donor( *grid, grid->descriptor ).get() ) );
  const ceelo::AnchorCurve tot = CeeLoUtils::totalTransferAnchorForDrf( with_total, anchor );
  BOOST_REQUIRE_EQUAL( tot.energies_keV.size(), anchor.curve.energies_keV.size() );
  for( size_t i = 0; i < tot.energies_keV.size(); ++i )  //sampled through a float-energy DRF call
    BOOST_CHECK_LT( rel_diff( tot.eff[i], with_total->ceeloResponse()->eps_total_at( tot.energies_keV[i], pos ).value ), 1.0e-6 );
}//BOOST_AUTO_TEST_CASE( GridAnchorsAndGrounding )


/** The grid-vs-Monte-Carlo FEP figure reads the donor's own near-field nodes: nothing against
 itself, and an injected offset back exactly. */
BOOST_AUTO_TEST_CASE( FepConsistency )
{
  const shared_ptr<const ceelo::DetectorResponse> grid = mock_drf()->ceeloResponse();

  const DetEffG2kPar::GridFepConsistency self = DetEffG2kPar::gridFepConsistency( *grid, *grid );
  BOOST_REQUIRE_GT( self.num_nodes, 50 );
  BOOST_CHECK_SMALL( self.rms_ln, 1.0e-12 );

  // A stand-in Monte Carlo 3% high at every near-field node, and 20% low at one on axis.
  const shared_ptr<ceelo::DetectorResponse> donor = ceelo::DetectorResponse::from_xml_string( grid->to_xml_string() );
  ceelo::NearFieldModel &nf = donor->near_field;
  for( double &v : nf.ln_n )
    v += 0.03;

  const double off = grid->descriptor.endcap_front_offset_cm();
  size_t d = 0;
  while( nf.dists_cm[d] < (off + 1.0) )
    ++d;
  BOOST_REQUIRE_LT( d + 1, nf.dists_cm.size() );
  const size_t c = nf.cos_thetas.size() - 1, e = nf.energies_keV.size() / 2;
  BOOST_REQUIRE_EQUAL( nf.cos_thetas[c], 1.0 );
  nf.ln_n[nf.index( e, c, d )] -= 0.23;
  nf.finalize();

  const DetEffG2kPar::GridFepConsistency r = DetEffG2kPar::gridFepConsistency( *grid, *donor );
  const double n = static_cast<double>( self.num_nodes );
  BOOST_REQUIRE_EQUAL( r.num_nodes, self.num_nodes );
  BOOST_CHECK_CLOSE( r.worst_ln, -0.20, 1.0e-6 );
  BOOST_CHECK_CLOSE( r.mean_ln, (0.03*(n - 1.0) - 0.20) / n, 1.0e-6 );
  BOOST_CHECK_CLOSE( r.rms_ln, std::sqrt( (0.0009*(n - 1.0) + 0.04) / n ), 1.0e-6 );
  BOOST_CHECK_CLOSE( r.worst_energy_keV, nf.energies_keV[e], 1.0e-9 );
  BOOST_CHECK_CLOSE( r.worst_dist_cm, nf.dists_cm[d] - off, 1.0e-9 );
  BOOST_CHECK_CLOSE( r.worst_cos_theta, 1.0, 1.0e-9 );
}//BOOST_AUTO_TEST_CASE( FepConsistency )


/** The edited geometry is kept beside the imported one: it is the DRF's Monte-Carlo geometry, part
 of its identity, and survives the file, the database column and a URL (which cannot carry the grid,
 and says so). */
BOOST_AUTO_TEST_CASE( EditedGeometryPersists )
{
  const shared_ptr<const DetectorPeakResponse> imported = mock_drf();
  const shared_ptr<const ceelo::DetectorResponse> grid = imported->ceeloResponse();
  const ceelo::GeometryDescriptor corrected = corrected_geometry();

  auto drf = make_shared<DetectorPeakResponse>( *imported );
  drf->setCeeloResponse( DetEffG2kPar::attachTotalEfficiency( *grid, total_donor( *grid, corrected ).get() ) );
  const uint64_t unedited_hash = drf->hashValue();

  // The imported geometry stated again is not an edit.
  drf->setGeometry( make_shared<const ceelo::GeometryDescriptor>( grid->descriptor ) );
  BOOST_CHECK( !drf->geometryModifiedFromImport() );
  BOOST_CHECK_EQUAL( drf->hashValue(), unedited_hash );

  drf->setGeometry( make_shared<const ceelo::GeometryDescriptor>( corrected ) );
  BOOST_CHECK( drf->geometryModifiedFromImport() );
  BOOST_CHECK( drf->hasImportedGrid() );
  BOOST_CHECK( drf->hashValue() != unedited_hash );
  BOOST_REQUIRE( drf->monteCarloGeometry() );
  BOOST_CHECK_EQUAL( drf->monteCarloGeometry()->to_xml_string(), corrected.to_xml_string() );
  BOOST_CHECK_EQUAL( drf->storedGeometry()->to_xml_string(), grid->descriptor.to_xml_string() );  //the grid's stays

  // The XML file.
  {
    rapidxml::xml_document<char> doc;
    rapidxml::xml_node<char> *root = doc.allocate_node( rapidxml::node_element, "root" );
    doc.append_node( root );
    drf->toXml( root, &doc );
    string xml;
    rapidxml::print( std::back_inserter(xml), doc, 0 );
    vector<char> buf( xml.begin(), xml.end() );
    buf.push_back( '\0' );
    rapidxml::xml_document<char> doc2;
    doc2.parse<0>( buf.data() );

    auto reread = make_shared<DetectorPeakResponse>();
    BOOST_REQUIRE_NO_THROW( reread->fromXml( doc2.first_node("root")->first_node("DetectorPeakResponse") ) );
    BOOST_CHECK( reread->hasImportedGrid() );
    BOOST_CHECK( reread->hasAnyTotalEfficiencyInfo() );
    BOOST_CHECK( reread->geometryModifiedFromImport() );
    BOOST_CHECK_EQUAL( reread->monteCarloGeometry()->to_xml_string(), corrected.to_xml_string() );
    BOOST_CHECK_EQUAL( reread->ceeloResponse()->content_hash(), drf->ceeloResponse()->content_hash() );
    BOOST_CHECK_EQUAL( reread->hashValue(), drf->hashValue() );
#if( PERFORM_DEVELOPER_CHECKS )
    BOOST_CHECK_NO_THROW( DetectorPeakResponse::equalEnough( *drf, *reread ) );
#endif
  }

  // The database's extra column.
  {
    const string extra = drf->drfExtraToXmlString();
    BOOST_REQUIRE( !extra.empty() );
    auto reread = make_shared<DetectorPeakResponse>( *imported );
    reread->setDrfExtraFromXmlString( extra );
    BOOST_CHECK( reread->geometryModifiedFromImport() );
    BOOST_CHECK( reread->hasAnyTotalEfficiencyInfo() );
    BOOST_CHECK_EQUAL( reread->monteCarloGeometry()->to_xml_string(), corrected.to_xml_string() );
  }

  // A URL carries no response, so it sends the edited geometry and says the grid is not included.
  {
    DetectorPeakResponse received;
    BOOST_REQUIRE_NO_THROW( received.fromAppUrl( drf->toAppUrl() ) );
    BOOST_CHECK( !received.ceeloResponse() );
    BOOST_CHECK( received.description().find( "efficiency grid not included" ) != string::npos );
    BOOST_REQUIRE( received.storedGeometry() );
    BOOST_CHECK_EQUAL( received.storedGeometry()->to_xml_string(), corrected.to_xml_string() );
  }

  // Detaching the grid (Flat Disk) keeps the edited geometry, as Modify's buildWorkingDrf does.
  {
    auto detached = make_shared<DetectorPeakResponse>( *drf );
    const shared_ptr<const ceelo::GeometryDescriptor> keep
                 = make_shared<const ceelo::GeometryDescriptor>( *detached->monteCarloGeometry() );
    detached->setCeeloResponse( nullptr );
    detached->setGeometry( keep );
    BOOST_CHECK( !detached->hasImportedGrid() );
    BOOST_CHECK( !detached->geometryModifiedFromImport() );
    BOOST_CHECK_EQUAL( detached->storedGeometry()->to_xml_string(), corrected.to_xml_string() );
  }
}//BOOST_AUTO_TEST_CASE( EditedGeometryPersists )
