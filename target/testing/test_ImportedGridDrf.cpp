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
#include <ctime>
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
#include "materials/Material.h"
#include "io/DetectorResponse.h"
#include "io/ResponseGenerator.h"

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


  /** A response for `gd` that carries a (synthetic, smooth) EtaTotTable total on the grid's own
   energy x angle nodes - the stand-in for a Monte Carlo of `gd` - and, optionally, a synthetic
   near-field total table like a General-profile run's.  Its FEP table is the grid's (without the
   grid's near-field table), so its total/FEP - all attachTotalEfficiency takes from it - is
   exp(eta_tot - eta_fep), times the near-field total table close in. */
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


  /** As #total_donor, but with a compact total tier - BCurve, eps_tot = b(E) * K_NoRayleigh, as CeeLo
   picks for scintillators.  Grounding scales such a response's FEP but not its total. */
  shared_ptr<ceelo::DetectorResponse> bcurve_donor( const ceelo::DetectorResponse &grid,
                                                    const ceelo::GeometryDescriptor &gd )
  {
    const shared_ptr<ceelo::DetectorResponse> donor = total_donor( grid, gd );
    donor->tot_eff = ceelo::TotEffPayload{};
    donor->tot_eff.tier = ceelo::TotEffTier::BCurve;
    for( const double energy : { 40.0, 150.0, 600.0, 3000.0 } )
    {
      donor->tot_eff.b_energies_keV.push_back( energy );
      donor->tot_eff.ln_b.push_back( 0.3 + 0.2*std::log( energy / 100.0 ) );
    }
    donor->finalize();
    return donor;
  }//bcurve_donor(...)


  double rel_diff( const double a, const double b )
  {
    return std::fabs( a - b ) / (std::max)( std::fabs(b), 1.0e-300 );
  }


  /** The total attachTotalEfficiency promises: the grid's FEP times the donor's total/FEP at the same
   physical point (`grid_pos` in the grid's frame, `donor_pos` in the donor's). */
  double grid_fep_times_donor_ratio( const ceelo::DetectorResponse &grid, const Eigen::Vector3d &grid_pos,
                                     const ceelo::DetectorResponse &donor, const Eigen::Vector3d &donor_pos,
                                     const double energy )
  {
    return grid.eps_fep_at( energy, grid_pos ).value * donor.eps_total_at( energy, donor_pos ).value
           / donor.eps_fep_at( energy, donor_pos ).value;
  }


  /** Whether node (c, d) of the total's close-in table `tn` is one of the grid's rungs inside the
   endcap - which makeDrf, and so attachTotalEfficiency, reads at the endcap face in front of it. */
  bool endcap_rung( const ceelo::DetectorResponse &grid, const ceelo::NearFieldModel &tn, const size_t c, const size_t d )
  {
    const ceelo::NearFieldModel &gn = grid.near_field;
    return (d < gn.dists_cm.size()) && (gn.dists_cm[d] == tn.dists_cm[d])
           && (gn.frac_sigma[gn.index( 0, c, d )] >= 1.0);
  }


  /** The endcap-face point in front of `pos`, both in `gd`'s crystal frame - where an endcap rung is read. */
  Eigen::Vector3d face_point_in_front( const ceelo::GeometryDescriptor &gd, const Eigen::Vector3d &pos )
  {
    const Eigen::Vector3d face = CeeLoUtils::detectorFacePosition( gd );
    return CeeLoUtils::sourcePositionFromFace( gd, 0.5 * 3.14159265358979323846, 0.0,
                                               std::hypot( pos.x() - face.x(), pos.y() - face.y() ) );
  }


  /** As above, with the caller's quadratures - built once per position, for loops over a lattice.
   A quadrature only serves the response instance that made it (it keys that instance's materials). */
  double grid_fep_times_donor_ratio( const ceelo::DetectorResponse &grid, const Eigen::Vector3d &grid_pos,
                                     const ceelo::ApertureQuadrature &q_grid,
                                     const ceelo::DetectorResponse &donor, const Eigen::Vector3d &donor_pos,
                                     const ceelo::ApertureQuadrature &q_donor, const double energy )
  {
    return grid.eps_fep_at( energy, grid_pos, q_grid ).value * donor.eps_total_at( energy, donor_pos, q_donor ).value
           / donor.eps_fep_at( energy, donor_pos, q_donor ).value;
  }
}//namespace


/** `t` without its side & back nodes (cos theta < 0) - the table a front-only import makes. */
ceelo::EtaTable strip_back_nodes( const ceelo::EtaTable &t )
{
  size_t nb = 0;
  while( (nb < t.cos_thetas.size()) && (t.cos_thetas[nb] < 0.0) )
    ++nb;

  ceelo::EtaTable u = t;
  u.cos_thetas.assign( t.cos_thetas.begin() + nb, t.cos_thetas.end() );
  u.ln_eta.assign( t.energies_keV.size() * u.cos_thetas.size(), 0.0 );
  u.frac_sigma.assign( u.ln_eta.size(), 0.0 );
  for( size_t e = 0; e < t.energies_keV.size(); ++e )
  {
    for( size_t c = 0; c < u.cos_thetas.size(); ++c )
    {
      u.ln_eta[u.index( e, c, 0 )] = t.ln_eta[t.index( e, c + nb, 0 )];
      u.frac_sigma[u.index( e, c, 0 )] = t.frac_sigma[t.index( e, c + nb, 0 )];
    }
  }
  u.finalize();
  return u;
}//strip_back_nodes(...)


/** The grid characterizes the whole sphere, so the import carries its columns past 90 degrees as
 side & back nodes, and the DETECTOR.txt record's endcap length as the space behind the crystal:
 the far-field grid is reproduced behind the detector, the front is unchanged by it, and a total
 attached from a front-only donor is charged for extrapolating there. */
BOOST_AUTO_TEST_CASE( SideAndBackFromTheGrid )
{
  const shared_ptr<const DetectorPeakResponse> drf = mock_drf();
  const shared_ptr<const ceelo::DetectorResponse> grid = drf->ceeloResponse();
  const ceelo::GeometryDescriptor &gd = grid->descriptor;

  // mock_DETECTOR.txt: a 50 mm long crystal, 4 mm front gap + 1.5 mm Al endcap, a 100 mm long can
  //  with nothing in the back slots - so 44.5 mm of vacuum behind the crystal, and side walls that
  //  run the can's length.
  BOOST_CHECK_CLOSE( gd.back_structure_extent_cm(), 10.0 - 0.55 - 5.0, 1.0e-6 );
  BOOST_CHECK_CLOSE( gd.outer_envelope().z_max_cm, 10.0 - 0.55, 1.0e-6 );
  BOOST_CHECK_CLOSE( gd.endcap_front_offset_cm(), 0.55, 1.0e-6 );

  // Every grid column past 90 degrees is a node, 180 degrees exactly; the file is v4.
  const ceelo::EtaTable &eta = grid->eta_fep;
  BOOST_REQUIRE( !eta.cos_thetas.empty() );
  BOOST_CHECK_EQUAL( eta.cos_thetas.front(), -1.0 );
  size_t n_back = 0;
  for( const double ct : eta.cos_thetas )
    n_back += (ct < 0.0) ? 1 : 0;
  BOOST_CHECK_EQUAL( n_back, 18 );    // 95-180 deg in 5 deg columns
  BOOST_CHECK( grid->to_xml_string().find( "version=\"4\"" ) != string::npos );

  // The far field behind the detector answers the grid on its lattice (row 80 is ~3 m, where the
  //  nodes are); closer in it is the kernel, through the descriptor's back, that carries it.
  const DetEffG2kPar::ParEfficiency par_eff( DetEffG2kPar::parseParFile( det_eff_path( "mock.par" ) ) );
  const DetEffG2kPar::ParGrid &grid0 = par_eff.parFile().grids.front();
  for( const double energy : { 60.0, 122.0, 662.0, 1400.0 } )
  {
    for( const int col : { 20, 27, 36 } )   // 100, 135 and 180 degrees
    {
      const double theta = col * grid0.theta_step_rad;
      for( const int row : { 69, 80 } )      // ~1 m and ~3 m
      {
        const double d_mm = std::exp( row * grid0.r_step );
        const double expected = par_eff.efficiency( energy, d_mm, theta );
        const DetectorPeakResponse::EffEval value
                  = drf->fepEfficiencyEval( static_cast<float>(energy), theta, 0.0, 0.1 * d_mm * PhysicalUnits::cm );
        BOOST_TEST_MESSAGE( energy << " keV, " << 0.1*d_mm << " cm, " << col*5 << " deg: DRF "
                            << value.value << " vs grid " << expected );
        if( row == 80 )
          BOOST_CHECK_MESSAGE( rel_diff( value.value, expected ) < 0.01,
                               energy << " keV, " << 0.1*d_mm << " cm, " << col*5 << " deg: DRF "
                               << value.value << " vs grid " << expected );
        const Eigen::Vector3d pos = CeeLoUtils::sourcePositionFromFace( gd, theta, 0.0, 0.1 * d_mm );
        const ceelo::EffResult r = grid->eps_fep_at( energy, pos );
        BOOST_CHECK( r.flag != ceelo::ResponseFlag::NeedsMc );
        BOOST_CHECK( r.flag != ceelo::ResponseFlag::OutOfRangeClamped );
      }
    }
  }

  // The front evaluates exactly as without the back nodes.
  const shared_ptr<ceelo::DetectorResponse> front_only = ceelo::DetectorResponse::from_xml_string( grid->to_xml_string() );
  front_only->eta_fep = strip_back_nodes( front_only->eta_fep );
  for( const double d_cm : { 1.0, 5.0, 25.0, 200.0 } )
  {
    for( const double theta_deg : { 0.0, 45.0, 85.0, 90.0 } )
    {
      const Eigen::Vector3d pos = CeeLoUtils::sourcePositionFromFace( gd, theta_deg*pi/180.0, 0.0, d_cm );
      for( const double energy : { 45.0, 300.0, 2500.0 } )
        BOOST_CHECK_EQUAL( front_only->eps_fep_at( energy, pos ).value, grid->eps_fep_at( energy, pos ).value );
    }
  }

  // A total from a donor characterized beside the detector keeps the back nodes; from a front-only
  //  donor it does not, and is charged for extrapolating there.
  const Eigen::Vector3d beside = CeeLoUtils::sourcePositionFromFace( gd, 95.0*pi/180.0, 0.0, 40.0 );
  BOOST_REQUIRE( beside.z() < gd.dimensions_cm.back() );   //not behind the crystal's back plane

  const shared_ptr<ceelo::DetectorResponse> donor = total_donor( *grid, gd );
  const shared_ptr<ceelo::DetectorResponse> covered = DetEffG2kPar::attachTotalEfficiency( *grid, donor.get() );
  BOOST_CHECK_EQUAL( covered->tot_eff.eta_tot.cos_thetas.front(), -1.0 );
  const ceelo::EffResult tot_covered = covered->eps_total_at( 662.0, beside );
  BOOST_CHECK_LT( rel_diff( tot_covered.value, grid_fep_times_donor_ratio( *grid, beside, *donor, beside, 662.0 ) ), 1.0e-6 );

  const shared_ptr<ceelo::DetectorResponse> front_donor = total_donor( *grid, gd );
  front_donor->tot_eff.eta_tot = strip_back_nodes( front_donor->tot_eff.eta_tot );
  const shared_ptr<ceelo::DetectorResponse> uncovered = DetEffG2kPar::attachTotalEfficiency( *grid, front_donor.get() );
  BOOST_CHECK( !(uncovered->tot_eff.eta_tot.cos_thetas.front() < 0.0) );
  const ceelo::EffResult tot_uncovered = uncovered->eps_total_at( 662.0, beside );
  BOOST_CHECK( tot_uncovered.value > 0.0 );
  BOOST_CHECK_GE( tot_uncovered.sigma_model,
                  0.999 * ceelo::model_sigma::behind_extrapolation * tot_uncovered.value );
  BOOST_CHECK_LT( tot_covered.sigma_model / tot_covered.value, tot_uncovered.sigma_model / tot_uncovered.value );
}//BOOST_AUTO_TEST_CASE( SideAndBackFromTheGrid )


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


/** On the grid's own geometry the total is the grid's FEP times the donor's total/FEP: exactly at
 the far-field pin's nodes (same kernel, same nodes), and - its close-in table being on the grid's
 own lattice - between them too, wherever the donor's ratio does not vary with distance (this
 donor's does not), while the FEP stays the grid's to the bit; stripping it restores the import.
 The gates sit a little above the printed figures. */
BOOST_AUTO_TEST_CASE( TotalOnImportedGeometry )
{
  const shared_ptr<const ceelo::DetectorResponse> grid = mock_drf()->ceeloResponse();
  const shared_ptr<ceelo::DetectorResponse> donor = total_donor( *grid, grid->descriptor );

  const std::clock_t cpu_start = std::clock();
  shared_ptr<ceelo::DetectorResponse> with_total;
  BOOST_REQUIRE_NO_THROW( with_total = DetEffG2kPar::attachTotalEfficiency( *grid, donor.get() ) );
  const double cpu_s = static_cast<double>( std::clock() - cpu_start ) / CLOCKS_PER_SEC;
  BOOST_REQUIRE( with_total );
  BOOST_CHECK( DetEffG2kPar::isGridResponse( with_total ) );
  BOOST_CHECK( with_total->tot_eff.characterized() );

  // The total's close-in table is on the grid's own lattice.
  const ceelo::NearFieldModel &gn = grid->near_field, &tn = with_total->tot_eff.near_field;
  BOOST_CHECK( tn.cos_thetas == gn.cos_thetas );
  for( const double energy : gn.energies_keV )
    BOOST_CHECK( std::find( begin(tn.energies_keV), end(tn.energies_keV), energy ) != end(tn.energies_keV) );
  BOOST_TEST_MESSAGE( "attachTotalEfficiency: " << cpu_s << " s CPU; total near-field table "
                      << tn.energies_keV.size() << " x " << tn.cos_thetas.size() << " x " << tn.dists_cm.size()
                      << " (the grid's " << gn.energies_keV.size() << " x " << gn.cos_thetas.size() << " x "
                      << gn.dists_cm.size() << "); serialized " << grid->to_xml_string().size() << " -> "
                      << with_total->to_xml_string().size() << " bytes" );

  const double d_pin = (std::max)( 1000.0 * grid->descriptor.transverse_half_extent(), 100.0 );
  for( size_t i = 0; i < grid->eta_fep.energies_keV.size(); i += 5 )
  {
    const Eigen::Vector3d pos = ceelo::source_position( d_pin, 1.0 );
    const double energy = grid->eta_fep.energies_keV[i];
    BOOST_CHECK_LT( rel_diff( with_total->eps_total_at( energy, pos ).value,
                              grid_fep_times_donor_ratio( *grid, pos, *donor, pos, energy ) ), 1.0e-9 );
  }

  // Past the pin, between the eta nodes: both tables' energy and angle interpolation.
  {
    double worst_far = 0.0;
    const vector<double> &energies = grid->eta_fep.energies_keV;
    for( const double cos_theta : { 1.0, 0.83, 0.5, 0.1 } )
    {
      const Eigen::Vector3d pos = ceelo::source_position( 3.0 * d_pin, cos_theta );
      for( size_t i = 0; (i + 1) < energies.size(); i += 3 )
      {
        const double energy = std::sqrt( energies[i] * energies[i + 1] );
        worst_far = (std::max)( worst_far, rel_diff( with_total->eps_total_at( energy, pos ).value,
                                                     grid_fep_times_donor_ratio( *grid, pos, *donor, pos, energy ) ) );
      }
    }
    BOOST_TEST_MESSAGE( "Total on the imported geometry, past the pin: worst |graft/(grid FEP x donor ratio) - 1| = "
                        << worst_far );
    BOOST_CHECK_LT( worst_far, 1.0e-6 );   //measured ~1e-16
  }

  //  Measured: 0.031 % at 1 cm - where the grazing rungs inside the endcap take the ratio at the
  //  endcap face in front of them - and ~1e-16 from 5 cm out.  (On a lattice of its own, which
  //  re-sampled the grid's close-in FEP: 0.27, 0.74, 0.16, 0.16, 0.13 %.)
  const vector<pair<double,double>> gates{ {1.0, 0.0005}, {5.0, 1.0e-6}, {25.0, 1.0e-6}, {50.0, 1.0e-6}, {200.0, 1.0e-6} };
  for( const pair<double,double> &gate : gates )
  {
    const double d_cm = gate.first;
    double worst_d = 0.0;
    for( const double theta_deg : { 0.0, 30.0, 60.0, 85.0 } )
    {
      const Eigen::Vector3d pos = CeeLoUtils::sourcePositionFromFace( grid->descriptor, theta_deg*pi/180.0, 0.0, d_cm );
      for( const double energy : { 45.0, 100.0, 300.0, 1000.0, 2500.0 } )
      {
        worst_d = (std::max)( worst_d, rel_diff( with_total->eps_total_at( energy, pos ).value,
                                                 grid_fep_times_donor_ratio( *grid, pos, *donor, pos, energy ) ) );
        BOOST_CHECK_EQUAL( with_total->eps_fep_at( energy, pos ).value, grid->eps_fep_at( energy, pos ).value );
      }
    }
    BOOST_TEST_MESSAGE( "Total on the imported geometry, " << d_cm << " cm: worst |graft/(grid FEP x donor ratio) - 1| = "
                        << worst_d );
    BOOST_CHECK_MESSAGE( worst_d < gate.second, "at " << d_cm << " cm: " << worst_d << " vs gate " << gate.second );
  }

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


/** Where nothing reaches the donor's crystal - 8 cm of lead around it underflows its kernel to 0 at
 the grid's lowest energy, at every angle - the donor's total/FEP comes from its kernel-free
 prefactors rather than the energy being refused; where the lead passes something, from its
 efficiencies, however small. */
BOOST_AUTO_TEST_CASE( TotalWhereKernelUnderflows )
{
  const shared_ptr<const ceelo::DetectorResponse> grid = mock_drf()->ceeloResponse();

  ceelo::GeometryDescriptor shielded = grid->descriptor;
  BOOST_REQUIRE( !shielded.layers.empty() );
  shielded.materials.push_back( ceelo::MaterialSpec::from( ceelo::make_Lead() ) );
  ceelo::LayerSpec lead = shielded.layers.back();
  lead.material_index = static_cast<int>( shielded.materials.size() ) - 1;
  lead.front_thickness_cm = lead.side_thickness_cm = 8.0;
  shielded.layers.push_back( lead );
  shielded.set_layer_back_thickness_cm( shielded.layers.size() - 1, 8.0 );  //behind it too: the grid has side & back nodes
  const shared_ptr<ceelo::DetectorResponse> donor = total_donor( *grid, shielded );

  // The same physical point for each response, as attachTotalEfficiency maps them.
  const double shift = shielded.endcap_front_offset_cm() - grid->descriptor.endcap_front_offset_cm();
  const auto donor_pos = [shift]( Eigen::Vector3d pos ) -> Eigen::Vector3d { pos.z() -= shift; return pos; };
  const double d_pin = (std::max)( 1000.0 * grid->descriptor.transverse_half_extent(), 100.0 );

  // The premise: at the lowest grid energy the donor has no total at any of the grid's angles.
  const ceelo::EtaTable &nodes = grid->eta_fep;
  const double e_lo = nodes.energies_keV.front();
  for( const double cos_theta : nodes.cos_thetas )
    BOOST_REQUIRE_EQUAL( donor->eps_total_at( e_lo, donor_pos( ceelo::source_position( d_pin, cos_theta ) ) ).value, 0.0 );

  shared_ptr<ceelo::DetectorResponse> with_total;
  BOOST_REQUIRE_NO_THROW( with_total = DetEffG2kPar::attachTotalEfficiency( *grid, donor.get() ) );

  // There, at the pin: the grid's FEP times the donor's own (kernel-free) total/FEP.
  for( const double cos_theta : { 1.0, 0.5 } )
  {
    const Eigen::Vector3d pos = ceelo::source_position( d_pin, cos_theta );
    const Eigen::Vector3d dpos = donor_pos( pos );
    const ceelo::ApertureQuadrature q_donor = donor->make_quadrature( dpos );
    const double ratio = donor->total_prefactor( e_lo, dpos, q_donor ).value
                         / donor->fep_prefactor( e_lo, dpos, q_donor ).value;
    const double fep = grid->eps_fep_at( e_lo, pos ).value;
    BOOST_REQUIRE_GT( fep, 0.0 );
    BOOST_CHECK_LT( rel_diff( with_total->eps_total_at( e_lo, pos ).value / fep, ratio ), 1.0e-9 );
  }

  // Close in, too, a finite positive total (the distance table holds where the donor has none).
  const double tot_close = with_total->eps_total_at( e_lo, CeeLoUtils::sourcePositionFromFace( grid->descriptor, 0.0, 0.0, 5.0 ) ).value;
  BOOST_CHECK( std::isfinite( tot_close ) && (tot_close > 0.0) );

  // Where the lead passes photons, the same product from the donor's efficiencies, at the nodes -
  //  however many decades the lead takes off them; the FEP is the grid's.
  for( const size_t e : { nodes.energies_keV.size() / 2, nodes.energies_keV.size() - 1 } )
  {
    const double energy = nodes.energies_keV[e];
    const Eigen::Vector3d pos = ceelo::source_position( d_pin, 1.0 );
    BOOST_REQUIRE_GT( donor->eps_total_at( energy, donor_pos( pos ) ).value, 0.0 );
    const double expected = grid_fep_times_donor_ratio( *grid, pos, *donor, donor_pos( pos ), energy );
    BOOST_CHECK_LT( rel_diff( with_total->eps_total_at( energy, pos ).value, expected ), 1.0e-6 );
    BOOST_CHECK_EQUAL( with_total->eps_fep_at( energy, pos ).value, grid->eps_fep_at( energy, pos ).value );
  }
}//BOOST_AUTO_TEST_CASE( TotalWhereKernelUnderflows )


/** A donor with a compact total tier, grounded as the app's worker grounds every Monte Carlo: the
 grounding k(E) moves its FEP but not its total, so the graft must take the total/FEP of the
 UNGROUNDED donor - at the far-field pin, at the close-in table's nodes, and, from the kernel-free
 prefactors, where 8 cm of lead underflows the donor's kernel. */
BOOST_AUTO_TEST_CASE( TotalFromGroundedCompactTierDonor )
{
  const shared_ptr<const ceelo::DetectorResponse> grid = mock_drf()->ceeloResponse();
  const double d_pin = (std::max)( 1000.0 * grid->descriptor.transverse_half_extent(), 100.0 );
  const ceelo::EtaTable &nodes = grid->eta_fep;

  for( const bool shielded : { false, true } )
  {
    ceelo::GeometryDescriptor gd = grid->descriptor;
    if( shielded )
    {
      gd.materials.push_back( ceelo::MaterialSpec::from( ceelo::make_Lead() ) );
      ceelo::LayerSpec lead = gd.layers.back();
      lead.material_index = static_cast<int>( gd.materials.size() ) - 1;
      lead.front_thickness_cm = lead.side_thickness_cm = 8.0;
      gd.layers.push_back( lead );
    }

    const shared_ptr<const ceelo::DetectorResponse> ungrounded = bcurve_donor( *grid, gd );
    const shared_ptr<ceelo::DetectorResponse> donor = ceelo::DetectorResponse::from_xml_string( ungrounded->to_xml_string() );

    // Grounded 13-20% high, on axis at 50 cm - where even the lead passes something.
    const Eigen::Vector3d ref_pos = CeeLoUtils::sourcePositionFromFace( gd, 0.0, 0.0, 50.0 );
    vector<ceelo::GroundingPoint> pts;
    for( const double energy : { 100.0, 300.0, 1000.0, 2500.0 } )
    {
      ceelo::GroundingPoint p;
      p.energy_keV = energy;
      p.model_eff = donor->eps_fep_at( energy, ref_pos ).value;
      p.measured_eff = p.model_eff * (1.2 - 0.05*std::log10( energy / 100.0 ));
      p.frac_stat_sigma = 0.01;
      p.distance_cm = 50.0;
      p.cos_theta = 1.0;
      BOOST_REQUIRE_GT( p.model_eff, 0.0 );
      pts.push_back( p );
    }
    BOOST_REQUIRE_NO_THROW( ceelo::ResponseGenerator::ground_to_points( *donor, pts, true ) );

    // The premise: grounding moved the FEP, and not the total.
    BOOST_REQUIRE_GT( donor->eps_fep_at( 300.0, ref_pos ).value / ungrounded->eps_fep_at( 300.0, ref_pos ).value, 1.1 );
    BOOST_REQUIRE_LT( rel_diff( donor->eps_total_at( 300.0, ref_pos ).value, ungrounded->eps_total_at( 300.0, ref_pos ).value ), 1.0e-12 );

    shared_ptr<ceelo::DetectorResponse> with_total;
    BOOST_REQUIRE_NO_THROW( with_total = DetEffG2kPar::attachTotalEfficiency( *grid, donor.get() ) );

    const double shift = gd.endcap_front_offset_cm() - grid->descriptor.endcap_front_offset_cm();
    const auto donor_pos = [shift]( Eigen::Vector3d pos ) -> Eigen::Vector3d { pos.z() -= shift; return pos; };

    // At the pin's nodes, where the donor's efficiencies are nonzero: the ungrounded donor's ratio,
    //  not the grounded one's (which grounding lowers by k).
    for( size_t e = 0; e < nodes.energies_keV.size(); e += 6 )
    {
      const double energy = nodes.energies_keV[e];
      for( const double cos_theta : { 1.0, 0.5 } )
      {
        const Eigen::Vector3d pos = ceelo::source_position( d_pin, cos_theta );
        if( !(ungrounded->eps_total_at( energy, donor_pos( pos ) ).value > 0.0) )
          continue;   //under the lead: checked from the prefactors below
        const double ours = with_total->eps_total_at( energy, pos ).value;
        BOOST_CHECK_LT( rel_diff( ours, grid_fep_times_donor_ratio( *grid, pos, *ungrounded, donor_pos( pos ), energy ) ), 1.0e-6 );
        BOOST_CHECK_GT( rel_diff( ours, grid_fep_times_donor_ratio( *grid, pos, *donor, donor_pos( pos ), energy ) ), 0.1 );
      }
    }

    if( !shielded )
    {
      // At every node of the close-in table too.
      const ceelo::NearFieldModel &tn = with_total->tot_eff.near_field;
      BOOST_REQUIRE( !tn.empty() );
      double worst_node = 0.0;
      for( size_t c = 0; c < tn.cos_thetas.size(); ++c )
        for( size_t d = 0; (d + 1) < tn.dists_cm.size(); ++d )
        {
          if( endcap_rung( *grid, tn, c, d ) )
            continue;   //read at the endcap face: checked in TotalOnEditedGeometry
          const Eigen::Vector3d pos = ceelo::source_position( tn.dists_cm[d], tn.cos_thetas[c] );
          const ceelo::ApertureQuadrature q_ours = with_total->make_quadrature( pos );
          const ceelo::ApertureQuadrature q_grid = grid->make_quadrature( pos );
          const ceelo::ApertureQuadrature q_donor = ungrounded->make_quadrature( pos );
          for( size_t e = 0; e < tn.energies_keV.size(); e += 3 )
          {
            if( tn.frac_sigma[tn.index( e, c, d )] >= 1.0 )
              continue;   //held: the donor had nothing there
            const double energy = tn.energies_keV[e];
            worst_node = (std::max)( worst_node, rel_diff( with_total->eps_total_at( energy, pos, q_ours ).value,
                                         grid_fep_times_donor_ratio( *grid, pos, q_grid, *ungrounded, pos, q_donor, energy ) ) );
          }
        }
      BOOST_TEST_MESSAGE( "Grounded BCurve donor, at the near-field nodes: worst |graft/(grid FEP x ungrounded ratio) - 1| = "
                          << worst_node );
      BOOST_CHECK_LT( worst_node, 1.0e-6 );
    }else
    {
      // Under the lead, at the grid's lowest energy: from the ungrounded donor's prefactors.
      const double e_lo = nodes.energies_keV.front();
      const Eigen::Vector3d pos = ceelo::source_position( d_pin, 1.0 );
      const Eigen::Vector3d dpos = donor_pos( pos );
      BOOST_REQUIRE_EQUAL( donor->eps_total_at( e_lo, dpos ).value, 0.0 );
      const ceelo::ApertureQuadrature q_donor = ungrounded->make_quadrature( dpos );
      const double ratio = ungrounded->total_prefactor( e_lo, dpos, q_donor ).value
                           / ungrounded->fep_prefactor( e_lo, dpos, q_donor ).value;
      const double fep = grid->eps_fep_at( e_lo, pos ).value;
      BOOST_REQUIRE_GT( fep, 0.0 );
      BOOST_CHECK_LT( rel_diff( with_total->eps_total_at( e_lo, pos ).value / fep, ratio ), 1.0e-9 );
    }//if( !shielded ) / else
  }//for( const bool shielded : { false, true } )
}//BOOST_AUTO_TEST_CASE( TotalFromGroundedCompactTierDonor )


/** A donor on an EDITED geometry - the corrected core and bulletizing - with its own near-field
 total, as a General-profile Monte Carlo of it would have.  The grid's total is its own FEP times the
 donor's total/FEP: exactly at the far-field pin and at every node of its own distance table, and
 between nodes to what interpolating the donor's ratio allows; the gates sit a little above the
 printed figures. */
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
                              grid_fep_times_donor_ratio( *grid, pos, *donor, donor_pos( pos ), energy ) ), 1.0e-6 );
  }

  // At every node of the grid's own near-field total table - which is the grid's lattice.  A rung
  //  inside the endcap takes the donor's total/FEP at the endcap face in front of it, where makeDrf
  //  read the grid's FEP for it, rather than at a point inside the can.
  const ceelo::NearFieldModel &tn = with_total->tot_eff.near_field;
  BOOST_REQUIRE( !tn.empty() );
  BOOST_CHECK( tn.cos_thetas == grid->near_field.cos_thetas );
  double worst_node = 0.0;
  size_t n_endcap_rungs = 0;
  for( size_t c = 0; c < tn.cos_thetas.size(); ++c )
    for( size_t d = 0; (d + 1) < tn.dists_cm.size(); ++d )
    {
      const Eigen::Vector3d pos = ceelo::source_position( tn.dists_cm[d], tn.cos_thetas[c] );
      const bool in_endcap = endcap_rung( *grid, tn, c, d );
      n_endcap_rungs += in_endcap ? 1 : 0;
      const Eigen::Vector3d ratio_pos = donor_pos( in_endcap ? face_point_in_front( grid->descriptor, pos ) : pos );
      const ceelo::ApertureQuadrature q_ours = with_total->make_quadrature( pos );
      const ceelo::ApertureQuadrature q_grid = grid->make_quadrature( pos );
      const ceelo::ApertureQuadrature q_donor = donor->make_quadrature( ratio_pos );
      for( size_t e = 0; e < tn.energies_keV.size(); e += 3 )
      {
        if( (tn.frac_sigma[tn.index( e, c, d )] >= 1.0) && !in_endcap )
          continue;   //held: the donor had nothing there
        const double energy = tn.energies_keV[e];
        const double tot = donor->eps_total_at( energy, ratio_pos, q_donor ).value;
        const double fep = donor->eps_fep_at( energy, ratio_pos, q_donor ).value;
        BOOST_REQUIRE( (tot > 0.0) && (fep > 0.0) );
        worst_node = (std::max)( worst_node, rel_diff( with_total->eps_total_at( energy, pos, q_ours ).value,
                                                       grid->eps_fep_at( energy, pos, q_grid ).value * tot / fep ) );
      }
    }
  BOOST_TEST_MESSAGE( "Total on the edited geometry, at the near-field nodes (" << n_endcap_rungs
                      << " rungs inside the endcap): worst |graft/(grid FEP x donor ratio) - 1| = " << worst_node );
  BOOST_CHECK_GT( n_endcap_rungs, 0 );
  BOOST_CHECK_LT( worst_node, 1.0e-6 );

  // Between nodes: at the same place relative to the endcap face.
  //  Measured: 0.015, 0.047, 0.018, 0.012, 0.004 %, then ~1e-16 past the donor's own table (28 cm),
  //  where its ratio no longer varies with distance.  (On a lattice of its own, which re-sampled the
  //  grid's close-in FEP: 0.34, 1.13, 0.37, 0.25, 0.27, 0.19, 0.05 %.)
  const vector<pair<double,double>> gates{ {1.0, 0.0003}, {2.0, 0.0007}, {5.0, 0.0003}, {10.0, 0.0002},
                                           {25.0, 0.0001}, {100.0, 1.0e-6}, {300.0, 1.0e-6} };
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
                                                 grid_fep_times_donor_ratio( *grid, grid_pos, *donor, other_pos, energy ) ) );
    }
    BOOST_TEST_MESSAGE( "Total on the edited geometry, " << d_cm << " cm (0-75 deg): worst |graft/(grid FEP x donor ratio) - 1| = "
                        << worst_d );
    BOOST_CHECK_MESSAGE( worst_d < gate.second, "at " << d_cm << " cm: " << worst_d << " vs gate " << gate.second );
  }
}//BOOST_AUTO_TEST_CASE( TotalOnEditedGeometry )


/** Between the close-in table's nodes - densely, from contact out, grazing to on-axis, at the file's
 energies and between them - the graft's total/FEP is the donor's ratio: the grid's own close-in FEP
 structure, interpolated the same way in both tables, cancels.  So the total cannot fall below the
 FEP where a (physical) donor's ratio is at least 1.  On the edited geometry, with the donor's own
 near-field total; the gate sits a little above the printed figure. */
BOOST_AUTO_TEST_CASE( TotalToPeakBetweenNodes )
{
  const shared_ptr<const ceelo::DetectorResponse> grid = mock_drf()->ceeloResponse();
  const ceelo::GeometryDescriptor corrected = corrected_geometry();
  const shared_ptr<ceelo::DetectorResponse> donor = total_donor( *grid, corrected, true );

  shared_ptr<ceelo::DetectorResponse> with_total;
  BOOST_REQUIRE_NO_THROW( with_total = DetEffG2kPar::attachTotalEfficiency( *grid, donor.get() ) );

  // The file's energies (every one is a node of the grid's close-in table) and the midpoints between them.
  const DetEffG2kPar::ParFile par = DetEffG2kPar::parseParFile( det_eff_path( "mock.par" ) );
  vector<double> energies;
  for( size_t i = 0; i < par.energies_keV.size(); ++i )
  {
    energies.push_back( par.energies_keV[i] );
    if( (i + 1) < par.energies_keV.size() )
      energies.push_back( std::sqrt( par.energies_keV[i] * par.energies_keV[i + 1] ) );
  }

  double worst = 0.0, worst_e = 0.0, worst_d = 0.0, worst_theta = 0.0;
  for( const double theta_deg : { 0.0, 12.0, 27.0, 41.0, 55.0, 66.0, 75.0, 82.0, 87.0 } )
  {
    const double theta = theta_deg * pi / 180.0;
    for( const double d_cm : { 0.3, 0.7, 1.3, 2.2, 3.7, 6.0, 9.5, 16.0, 27.0, 45.0, 80.0, 150.0 } )
    {
      const Eigen::Vector3d grid_pos = CeeLoUtils::sourcePositionFromFace( grid->descriptor, theta, 0.0, d_cm );
      const Eigen::Vector3d other_pos = CeeLoUtils::sourcePositionFromFace( corrected, theta, 0.0, d_cm );
      const ceelo::ApertureQuadrature q_ours = with_total->make_quadrature( grid_pos );
      const ceelo::ApertureQuadrature q_grid = grid->make_quadrature( grid_pos );
      const ceelo::ApertureQuadrature q_donor = donor->make_quadrature( other_pos );
      for( const double energy : energies )
      {
        const double ours = with_total->eps_total_at( energy, grid_pos, q_ours ).value;
        const double expected = grid_fep_times_donor_ratio( *grid, grid_pos, q_grid, *donor, other_pos, q_donor, energy );
        if( !(expected > 0.0) )
          continue;
        BOOST_REQUIRE( std::isfinite( ours ) );
        const double diff = rel_diff( ours, expected );
        if( diff > worst )
        {
          worst = diff;
          worst_e = energy;
          worst_d = d_cm;
          worst_theta = theta_deg;
        }
      }
    }
  }

  //  Measured: 0.13 %, at 40 keV, 0.7 cm and 87 deg - the grazing near-contact corner, where rungs
  //  inside the endcap take the ratio at the face in front of them.
  BOOST_TEST_MESSAGE( "Total/FEP between nodes, edited geometry: worst |graft/(grid FEP x donor ratio) - 1| = "
                      << worst << " (" << worst_e << " keV, " << worst_d << " cm, " << worst_theta << " deg)" );
  BOOST_CHECK_LT( worst, 2.0e-3 );
}//BOOST_AUTO_TEST_CASE( TotalToPeakBetweenNodes )


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
