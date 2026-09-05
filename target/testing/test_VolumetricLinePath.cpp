/* InterSpec: an application to analyze spectral gamma radiation data.

 Copyright 2026 National Technology & Engineering Solutions of Sandia, LLC
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


/** The detector-side LINE integration path (InterSpec/VolumetricLineIntegration_imp.hpp) against
 the per-element path on everything the line path provides BEYOND the plain volume integral (that
 integral itself is cross-checked by LineVsElementScenarioMatrix in test_VolumetricLadder.cpp):

   - the per-point shell walk (#shell_path_to_point_imp) that hands the cascade correction and the
     effective-shielding report the same per-shell path lengths the element walkers accumulate on
     their centre ray - pinned against eval_cylinder / eval_rect / eval_spherical directly, since the
     two walkers deliberately share no code;
   - the effective-shielding components accumulated per line (#integrate_effective_shielding_all)
     against the element path's #integrate_effective_shielding.

 No Monte Carlo; the whole file runs in seconds.
 */

#include <array>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>
#include <iostream>

#define BOOST_TEST_MODULE VolumetricLinePath_suite
#include <boost/test/included/unit_test.hpp>

#include "ceres/jet.h"

#include "VolumetricNearFieldHarness.h"


namespace
{
/** Snapshot of the per-ray scratch a shell walk leaves on a calculator. */
struct WalkScratch
{
  double ad = 0.0, an_ad = 0.0, ad_h = 0.0, mu_d = 0.0, an_mu_d = 0.0;
  std::vector<double> partner_mu_d;
  double cascade_ad = 0.0, cascade_an_ad = 0.0, air = 0.0;

  static WalkScratch take( const GammaInteractionCalc::DistributedSrcCalcT<double> &calc )
  {
    WalkScratch w;
    w.ad = calc.m_ray_ad;
    w.an_ad = calc.m_ray_an_ad;
    w.ad_h = calc.m_ray_ad_h;
    w.mu_d = calc.m_ray_mu_d;
    w.an_mu_d = calc.m_ray_an_mu_d;
    w.partner_mu_d = calc.m_ray_partner_mu_d;
    w.cascade_ad = calc.m_ray_cascade_ad;
    w.cascade_an_ad = calc.m_ray_cascade_an_ad;
    w.air = calc.m_ray_air_dist;
    return w;
  }
};//struct WalkScratch


void check_close( const char *what, const double a, const double b, const double rel, const std::string &where )
{
  const double scale = std::max( std::fabs(a), std::fabs(b) );
  BOOST_CHECK_MESSAGE( std::fabs(a - b) <= rel*scale + 1.0e-14,
                       where << ": " << what << " element=" << a << " walk=" << b );
}
}//namespace


/** The line-side shell walk against the element walkers' centre ray.

 Five nested shells - a steel core, a generic layer on it, the water SOURCE, a generic layer on
 that, and a steel jacket - in every geometry, so every row of the convention table in
 #ShellPathT is exercised: inner Generic counted only when crossed, outer Generic always, outer
 Material as its near segment, inner Material as its full nested chord, the source as its own near
 piece, and the air gap to the detector.  Each Material shell carries distinct per-partner cascade
 coefficients, so a shell attributed to the wrong neighbour shows up in the partner sums even when
 the total path length happens to agree.
 */
BOOST_AUTO_TEST_CASE( ShellWalkMatchesElementCentreRay )
{
  using namespace GammaInteractionCalc;
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const double cm = PhysicalUnits::cm;

  const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );
  const shared_ptr<const Material> steel = matdb->material( "Stainless steel SS-304" );
  const shared_ptr<const Material> water = matdb->material( "Water" );
  BOOST_REQUIRE( steel && water );

  const std::vector<double> partners = { 100.0, 500.0 };
  const double energy = 661.7;

  struct Stack { const char *name; GeometryType geom; std::array<double,3> core, src, outer; double dist; };
  const std::vector<Stack> stacks = {
    { "cylEnd",  GeometryType::CylinderEndOn,  {1.0,0.8,0.0},   {3.0,2.0,0.0},   {3.4,2.5,0.0},   4.0 },
    { "cylSide", GeometryType::CylinderSideOn, {1.0,1.5,0.0},   {2.5,3.0,0.0},   {2.9,3.4,0.0},   5.0 },
    { "rect",    GeometryType::Rectangular,    {1.0,0.8,0.6},   {2.5,2.0,1.5},   {2.9,2.4,1.9},   5.0 },
    { "sphere",  GeometryType::Spherical,      {1.0,0.0,0.0},   {2.5,0.0,0.0},   {2.9,0.0,0.0},   5.0 },
  };

  const auto material_shell = [&]( const shared_ptr<const Material> &mat, const std::array<double,3> &dims ) {
    DistributedSrcCalcT<double>::ShellInfo sh;
    for( int i = 0; i < 3; ++i )
      sh.dims[i] = dims[i]*cm;
    sh.trans_len_coef = transmition_length_coefficient( mat.get(), static_cast<float>(energy) );
    sh.type = ShellType::Material;
    sh.density = mat->density;
    sh.effective_an = material_mass_weighted_atomic_number( *mat );
    sh.hydrogen_mass_frac = material_hydrogen_mass_fraction( *mat );
    for( const double pe : partners )
      sh.cascade_mu.push_back( transmition_length_coefficient( mat.get(), static_cast<float>(pe) ) );
    return sh;
  };
  const auto generic_shell = [&]( const std::array<double,3> &dims, const double ad_gcm2, const double an ) {
    DistributedSrcCalcT<double>::ShellInfo sh;
    for( int i = 0; i < 3; ++i )
      sh.dims[i] = dims[i]*cm;
    sh.trans_len_coef = 0.05*ad_gcm2;   //dimensionless total attenuation of the layer (arbitrary)
    sh.type = ShellType::Generic;
    sh.areal_density = ad_gcm2 * PhysicalUnits::g / PhysicalUnits::cm2;
    sh.effective_an = an;
    sh.hydrogen_mass_frac = 0.0;
    sh.cascade_mu = { 0.03*ad_gcm2, 0.01*ad_gcm2 };
    return sh;
  };

  size_t num_points_checked = 0;
  for( const Stack &st : stacks )
  {
    DistributedSrcCalcT<double> calc;
    calc.m_geometry = st.geom;
    calc.m_materialIndex = 2;
    calc.m_attenuateForAir = true;      //so the air distance is recorded ...
    calc.m_airTransLenCoef = 0.0;       // ... without attenuating anything
    calc.m_isInSituExponential = false;
    calc.m_inSituRelaxationLength = -1.0;
    calc.m_srcVolumetricActivity = 1.0;
    calc.m_normalizeByVolume = false;
    calc.m_energy = energy;
    calc.m_effResponse = nullptr;       //flat-disk response: the walk is what is under test
    calc.m_effMethod = ShieldingSourceFitCalc::VolumetricEffMethod::FlatDisk;
    calc.m_detector = detector_geom_from_config<double>( st.geom, st.dist*cm, 3.0*cm, 0.0 );
    calc.m_accumulateEffectiveAnAd = true;

    // A cascade block with no engine: the walk fills the partner scratch, nothing consumes it.
    auto block = std::make_shared<DistributedSrcCalcT<double>::CascadeBlock>();
    block->calc = nullptr;
    block->partner_energies = partners;
    block->partner_fep_int = { 0.1, 0.1 };
    block->partner_tot_int = { 0.2, 0.2 };
    block->air_mu = { 0.0, 0.0 };
    calc.m_cascade = block;

    calc.m_shells.push_back( material_shell( steel, st.core ) );
    calc.m_shells.push_back( generic_shell( st.core, 1.3, 13.0 ) );
    calc.m_shells.push_back( material_shell( water, st.src ) );
    calc.m_shells.push_back( generic_shell( st.src, 0.7, 29.0 ) );
    calc.m_shells.push_back( material_shell( steel, st.outer ) );
    calc.finalize_shell_coefficients();

    const std::array<double,3> &Do = st.src, &Di = st.core;
    size_t checked = 0;
    for( int i = 1; i <= 40; ++i )
    {
      const double xx[3] = { ceelo::halton( i, 2 ), ceelo::halton( i, 3 ), ceelo::halton( i, 5 ) };

      // The emission point the element evaluator maps these unit-cube coordinates to, in the
      //  assembly frame (the on-axis detector keeps eval_spherical's rotated frame equal to it).
      double point[3] = { 0.0, 0.0, 0.0 };
      bool in_core = false;
      switch( st.geom )
      {
        case GeometryType::CylinderEndOn:
        case GeometryType::CylinderSideOn:
        {
          const double r = xx[0]*Do[0], theta = xx[1]*2.0*PhysicalUnits::pi, z = 2.0*Do[1]*(xx[2] - 0.5);
          point[0] = r*std::cos(theta);
          point[1] = r*std::sin(theta);
          point[2] = z;
          in_core = (r < Di[0]) && (std::fabs(z) < Di[1]);
          break;
        }
        case GeometryType::Rectangular:
          for( int k = 0; k < 3; ++k )
            point[k] = (xx[k] - 0.5)*2.0*Do[k];
          in_core = (std::fabs(point[0]) < Di[0]) && (std::fabs(point[1]) < Di[1]) && (std::fabs(point[2]) < Di[2]);
          break;
        case GeometryType::Spherical:
        {
          const double r = Di[0] + xx[0]*(Do[0] - Di[0]);
          const double theta = xx[1]*PhysicalUnits::pi, phi = xx[2]*2.0*PhysicalUnits::pi;
          point[0] = r*std::sin(theta)*std::cos(phi);
          point[1] = r*std::sin(theta)*std::sin(phi);
          point[2] = r*std::cos(theta);
          break;
        }
        case GeometryType::NumGeometryType:
          BOOST_REQUIRE( false );
      }
      if( in_core )
        continue;   //the element evaluator returns before walking anything
      for( int k = 0; k < 3; ++k )
        point[k] *= cm;

      double value = 0.0;
      switch( st.geom )
      {
        case GeometryType::CylinderEndOn:
        case GeometryType::CylinderSideOn: value = calc.eval_cylinder( xx, 3 ); break;
        case GeometryType::Rectangular:    value = calc.eval_rect( xx, 3 ); break;
        case GeometryType::Spherical:      value = calc.eval_spherical( xx, 3 ); break;
        case GeometryType::NumGeometryType: break;
      }
      BOOST_REQUIRE( std::isfinite(value) );
      const WalkScratch element = WalkScratch::take( calc );

      ShellPathT<double> path;
      shell_path_from_point_imp( calc, point, path );
      record_shell_path_imp( calc, path );
      const WalkScratch walk = WalkScratch::take( calc );

      std::ostringstream where;
      where << st.name << " point " << i << " (" << point[0]/cm << "," << point[1]/cm << "," << point[2]/cm << ") cm";
      check_close( "AD", element.ad, walk.ad, 1.0e-10, where.str() );
      check_close( "AN*AD", element.an_ad, walk.an_ad, 1.0e-10, where.str() );
      check_close( "AD_H", element.ad_h, walk.ad_h, 1.0e-10, where.str() );
      check_close( "mu*d", element.mu_d, walk.mu_d, 1.0e-10, where.str() );
      check_close( "AN*mu*d", element.an_mu_d, walk.an_mu_d, 1.0e-10, where.str() );
      BOOST_REQUIRE( element.partner_mu_d.size() == partners.size() );
      BOOST_REQUIRE( walk.partner_mu_d.size() == partners.size() );
      for( size_t j = 0; j < partners.size(); ++j )
        check_close( "partner mu*d", element.partner_mu_d[j], walk.partner_mu_d[j], 1.0e-10, where.str() );
      check_close( "cascade AD", element.cascade_ad, walk.cascade_ad, 1.0e-10, where.str() );
      check_close( "cascade AN*AD", element.cascade_an_ad, walk.cascade_an_ad, 1.0e-10, where.str() );
      check_close( "air", element.air, walk.air, 1.0e-10, where.str() );
      ++checked;
    }//for( points )

    BOOST_CHECK_MESSAGE( checked >= 20, st.name << ": only " << checked << " points fell outside the core" );
    BOOST_TEST_MESSAGE( "  " << st.name << ": " << checked << " emission points, walk == element centre ray" );
    num_points_checked += checked;
  }//for( stacks )

  BOOST_CHECK( num_points_checked >= 80 );
}//BOOST_AUTO_TEST_CASE( ShellWalkMatchesElementCentreRay )


/** The effective-shielding report on the line path against the element path.

 The two weight the per-shell quantities differently, by construction: the element path weights each
 element's CENTRE-RAY areal density by that element's whole (aperture-averaged) contribution, while
 the line path weights every line's OWN path by that line's contribution - so lines through less
 material, which contribute more, pull the line path's <AD> lower.  The gap is the within-aperture
 covariance of transmission and areal density.  MEASURED (2026-09-03, ANGLE GEM35-70 transfer):
 0.2% or less far from the detector, where the two must and do coincide (and where the accumulation
 arithmetic is therefore pinned), but 9-42% for an optically thick source at CONTACT - the centre
 ray of an element near the rim runs far more steel than the lines that actually carry that element's
 counts.  The line path's number is the attenuation-weighted path of the detected photons, which is
 what "effective shielding" means; the element's is the same quantity for a fictitious centre ray.
 The AN ratios are immune (the mix of materials along a path barely depends on the path).  Gates:
 far rows 1% on every ratio; contact rows 1% on the AN ratios, and <AD> must be LOWER on the line
 path (Jensen's inequality: the short paths dominate exp(-tau)) by at most a factor of two.  The
 c[0] component must equal the plain line integral exactly, since it is the same sum.
 */
BOOST_AUTO_TEST_CASE( EffectiveShieldingLineVsElement )
{
  using namespace GammaInteractionCalc;
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );
  const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );

  struct Case { const char *name; bool contact; };
  const std::vector<Case> cases = {
    { "small-far-light", false },
    { "large-far-dense", false },
    { "large-near-dense", true },
    { "shielded-near-dense", true },
    { "box-shielded-near-dense", true },
  };
  const int num_lines = 1 << 16;

  for( const Case &c : cases )
  {
    const Scenario s = find_scenario( c.name );
    const shared_ptr<const Material> matrix = matdb->material( scenario_matrix_material( s.dense ) );
    const shared_ptr<const Material> iron = matdb->material( scenario_shield_material() );
    BOOST_REQUIRE( matrix && iron );

    for( const double e : { 60.0, 661.7 } )
    {
      DistributedSrcCalcT<double> calc = build_scenario_calc( det, s, e, det.mc_transfer );
      // The harness builds shells without the per-shell metadata the report needs.
      BOOST_REQUIRE( (calc.m_shells.size() == 1) || (calc.m_shells.size() == 2) );
      for( size_t i = 0; i < calc.m_shells.size(); ++i )
      {
        const shared_ptr<const Material> &mat = (i == 0) ? matrix : iron;
        calc.m_shells[i].density = mat->density;
        calc.m_shells[i].effective_an = material_mass_weighted_atomic_number( *mat );
        calc.m_shells[i].hydrogen_mass_frac = material_hydrogen_mass_fraction( *mat );
      }

      // Element path: the reference report.
      EffShieldComponents elem;
      {
        DistributedSrcCalcT<double> ref = calc;
        const ScopedVolumetricIntegratorOverride force( VolumetricIntegrator::Element );
        elem = integrate_effective_shielding( ref );
      }

      // Line path: through the production entry point, plus the plain integral for the identity.
      EffShieldComponents line;
      double line_integral = 0.0;
      {
        DistributedSrcCalcT<double> lc = calc;
        attach_line_cache( lc, num_lines );
        std::vector<std::unique_ptr<DistributedSrcCalcT<double>>> v;
        v.push_back( std::make_unique<DistributedSrcCalcT<double>>( lc ) );
        const ScopedVolumetricIntegratorOverride force( VolumetricIntegrator::Line );
        const std::vector<EffShieldComponents> comps = integrate_effective_shielding_all( v, true );
        BOOST_REQUIRE( comps.size() == 1 );
        line = comps[0];
        integrate_volumetric_calculators<double>( v, true );
        line_integral = v[0]->integral;
      }

      BOOST_REQUIRE( (elem.c[0] > 0.0) && (elem.c[1] > 0.0) && (elem.c[4] > 0.0) );
      BOOST_REQUIRE( (line.c[0] > 0.0) && (line.c[1] > 0.0) && (line.c[4] > 0.0) );

      const double ad_e = elem.c[1]/elem.c[0], ad_l = line.c[1]/line.c[0];
      const double an_e = elem.c[2]/elem.c[1], an_l = line.c[2]/line.c[1];
      const double h_e = elem.c[3]/elem.c[1], h_l = line.c[3]/line.c[1];
      const double anxs_e = elem.c[5]/elem.c[4], anxs_l = line.c[5]/line.c[4];

      std::ostringstream row;
      row << "  " << std::left << std::setw(24) << c.name << std::right << " @ " << std::setw(6)
          << std::fixed << std::setprecision(1) << e << " keV:  c0 " << std::showpos
          << std::setprecision(3) << 100.0*(line.c[0]/elem.c[0] - 1.0) << "%  <AD> "
          << 100.0*(ad_l/ad_e - 1.0) << "%  <AN> " << 100.0*(an_l/an_e - 1.0) << "%  <AN_xs> "
          << 100.0*(anxs_l/anxs_e - 1.0) << "%" << std::noshowpos;
      if( h_e > 0.0 )
        row << "  <fracH> " << std::showpos << 100.0*(h_l/h_e - 1.0) << "%" << std::noshowpos;
      row << "   (<AD> element " << std::setprecision(4) << ad_e/(PhysicalUnits::g/PhysicalUnits::cm2) << " g/cm2)";
      BOOST_TEST_MESSAGE( row.str() );

      const std::string where = std::string(c.name) + " @ " + std::to_string(e) + " keV";
      BOOST_CHECK_MESSAGE( std::fabs(line.c[0] - line_integral) <= 1.0e-12*line_integral,
                           where << ": c[0] " << line.c[0] << " != the plain line integral " << line_integral );
      BOOST_CHECK_MESSAGE( std::fabs(line.c[0]/elem.c[0] - 1.0) < 0.0075,
                           where << ": c[0] differs from the element path by " << 100.0*(line.c[0]/elem.c[0] - 1.0) << "%" );
      if( c.contact )
        BOOST_CHECK_MESSAGE( (ad_l < ad_e) && (ad_l > 0.5*ad_e),
                             where << ": contact <AD> on the line path (" << ad_l << ") should be below the"
                             " element path's centre-ray value (" << ad_e << ") by less than a factor of two" );
      else
        BOOST_CHECK_MESSAGE( std::fabs(ad_l/ad_e - 1.0) < 0.01,
                             where << ": <AD> differs by " << 100.0*(ad_l/ad_e - 1.0) << "%" );
      BOOST_CHECK_MESSAGE( std::fabs(an_l/an_e - 1.0) < 0.01,
                           where << ": <AN> differs by " << 100.0*(an_l/an_e - 1.0) << "%" );
      BOOST_CHECK_MESSAGE( std::fabs(anxs_l/anxs_e - 1.0) < 0.01,
                           where << ": <AN_xs> differs by " << 100.0*(anxs_l/anxs_e - 1.0) << "%" );
      if( h_e > 0.0 )
        BOOST_CHECK_MESSAGE( std::fabs(h_l/h_e - 1.0) < (c.contact ? 0.05 : 0.01),
                             where << ": <fracH> differs by " << 100.0*(h_l/h_e - 1.0) << "%" );
    }//for( energies )
  }//for( cases )
}//BOOST_AUTO_TEST_CASE( EffectiveShieldingLineVsElement )


/** A COLLIMATED response on the line path.

 The collimator is part of the geometry every line is traced through, so the line kernel carries
 its shadow in the VALUE already; what the line path lacked was the shadow gate in
 `DetectorResponse::common_eval`, which only sets the flag and the sigma and which the prefactor grid
 used to skip by handing it an empty quadrature.  The line cache now builds a coarse grid of gate
 quadratures (#CollimatorGateGrid) for a collimated response, and the fit's per-peak flags take the
 worst flag along the source's chords into account.

 Synthetic, MC-free transfer response (as VolumetricLinePathZeroThicknessLimit builds): a 3"x3" NaI
 in an Al can with a 1 cm lead collimator tube reaching 5 cm in front of the crystal.  A small source
 on the axis sees the whole crystal - line and element agree, flag Ok; a wide disc centred on the same
 axis reaches out into the shadow - line and element still agree (both trace the lead), the point
 query at its centre still says Ok, and the line cache's chord-range flag says Shadowed or worse,
 which is the case the point query alone cannot see.
 */
BOOST_AUTO_TEST_CASE( CollimatedResponseLinePath )
{
  using namespace GammaInteractionCalc;
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const double cm = PhysicalUnits::cm;
  const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );
  const shared_ptr<const Material> water = matdb->material( "Water" );
  BOOST_REQUIRE( water );

  ceelo::GeometryDescriptor gd;
  gd.shape = ceelo::DetectorShape::Cylinder;
  gd.dimensions_cm = { 3.81, 7.62 };
  gd.materials = { ceelo::MaterialSpec::from( ceelo::make_NaI() ),
                   ceelo::MaterialSpec::from( ceelo::make_Aluminum() ),
                   ceelo::MaterialSpec::from( ceelo::make_Lead() ) };
  gd.crystal_material_index = 0;
  ceelo::LayerSpec can;
  can.material_index = 1;
  can.front_thickness_cm = 0.05;
  can.side_thickness_cm = 0.05;
  can.z_end_cm = 7.62;
  gd.layers.push_back( can );
  ceelo::CollimatorSpec collimator;
  collimator.material_index = 2;
  collimator.side_thickness_cm = 1.0;
  collimator.z_start_cm = -5.0;
  collimator.z_end_cm = 7.62;
  gd.collimator = collimator;

  ceelo::AnchorCurve anchor;
  anchor.energies_keV = { 60.0, 100.0, 300.0, 662.0, 1000.0 };
  anchor.eff = { 2.0e-2, 1.3e-2, 5.0e-3, 3.0e-3, 2.2e-3 };
  anchor.frac_sigma = { 0.003, 0.003, 0.003, 0.003, 0.003 };

  std::shared_ptr<const ceelo::DetectorResponse> response;
  BOOST_REQUIRE_NO_THROW( response = ceelo::make_transfer_response( gd, anchor, Eigen::Vector3d( 0.0, 0.0, -30.0 ) ) );
  BOOST_REQUIRE( response );
  BOOST_REQUIRE( response->descriptor.collimator );

  const double energy = 661.7;
  struct Case { const char *name; double radius_cm, half_len_cm, dist_cm; bool expect_shadow; };
  const std::vector<Case> cases = {
    { "small on-axis source at 30 cm",   1.0, 0.5, 30.0, false },
    { "wide disc reaching into the shadow", 12.0, 0.5, 20.0, true },
  };

  for( const Case &c : cases )
  {
    DistributedSrcCalcT<double> calc;
    calc.m_geometry = GeometryType::CylinderEndOn;
    calc.m_materialIndex = 0;
    calc.m_attenuateForAir = false;
    calc.m_isInSituExponential = false;
    calc.m_inSituRelaxationLength = -1.0;
    calc.m_srcVolumetricActivity = 1.0;
    calc.m_normalizeByVolume = false;
    calc.m_energy = energy;
    calc.m_effResponse = response;
    calc.m_effMethod = ShieldingSourceFitCalc::VolumetricEffMethod::MCTransfer;
    calc.m_detector = detector_geom_from_config<double>( GeometryType::CylinderEndOn, c.dist_cm*cm,
                                                         response->transverse_half_extent()*cm, 0.0 );
    DistributedSrcCalcT<double>::ShellInfo shell;
    shell.dims = { c.radius_cm*cm, c.half_len_cm*cm, 0.0 };
    shell.trans_len_coef = transmition_length_coefficient( water.get(), static_cast<float>(energy) );
    shell.type = ShellType::Material;
    calc.m_shells.push_back( shell );

    DistributedSrcCalcT<double> elem = calc, line = calc;
    integrate_on_path( elem, VolumetricIntegrator::Element, -1 );
    integrate_on_path( line, VolumetricIntegrator::Line, 1 << 15 );
    BOOST_REQUIRE( line.m_lineCache );
    BOOST_CHECK( line_path_applicable( line ) );
    BOOST_REQUIRE( (elem.integral > 0.0) && (line.integral > 0.0) );
    const double rel = 100.0*(line.integral/elem.integral - 1.0);

    const ceelo::EffResult centre = response->eps_fep( energy, 0.0, 0.0, c.dist_cm );
    const ceelo::ResponseFlag chord_flag = line.m_lineCache->worst_flag( energy );

    BOOST_TEST_MESSAGE( "  " << c.name << ": line/element - 1 = " << std::fixed << std::showpos
                        << std::setprecision(3) << rel << "%" << std::noshowpos << "; point-query flag "
                        << ceelo::to_string( centre.flag ) << ", chord-range flag "
                        << ceelo::to_string( chord_flag ) );

    BOOST_CHECK_MESSAGE( std::fabs(rel) < 1.0, c.name << ": line and element disagree by " << rel << "%" );
    BOOST_CHECK_MESSAGE( centre.flag == ceelo::ResponseFlag::Ok,
                         c.name << ": the on-axis point query should be Ok, got " << ceelo::to_string( centre.flag ) );
    if( c.expect_shadow )
      BOOST_CHECK_MESSAGE( static_cast<int>(chord_flag) >= static_cast<int>(ceelo::ResponseFlag::Shadowed),
                           c.name << ": the chord-range flag should report the shadow, got "
                                  << ceelo::to_string( chord_flag ) );
    else
      BOOST_CHECK_MESSAGE( chord_flag == ceelo::ResponseFlag::Ok,
                           c.name << ": the chord-range flag should be Ok, got " << ceelo::to_string( chord_flag ) );
  }//for( cases )
}//BOOST_AUTO_TEST_CASE( CollimatedResponseLinePath )


/** SPHERICAL sources: the line path against the element path.

 The scenario matrix has no spheres (it is cylinders and boxes), so this is the only place the two
 quadratures are compared on the geometry that is InterSpec's DEFAULT for a shielding stack.  Solid
 and hollow, bare and shielded, at contact and far, so the source chord, the inner-core split and the
 outer-shell walk are each exercised.

 Both paths must apply the response through the aperture the source element actually sees.  Until
 2026-09-03 `eval_spherical` did not: it used the flat-disk solid angle scaled by a single
 centre-ray response, which is the same class of error the rectangles carried before their per-ray
 kernel landed.  Measured here before the fix, the element path sat 52% below the line path on a
 solid steel sphere at contact (60 keV), 19% below at 662 keV, 37%/13% on a hollow one and 6.6%/3.3%
 on water - while the FAR rows already agreed to 0.05%, which is the signature of a purely
 near-field aperture error.

 Both quadratures are models.  Rung8_SphericalSourceTruth in test_VolumetricLadder.cpp is the
 Monte-Carlo leg (developer-only, real MC); this case is the cheap consistency check that runs every
 time.
 */
BOOST_AUTO_TEST_CASE( LineVsElementSphericalSource )
{
  using namespace GammaInteractionCalc;
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );
  const double cm = PhysicalUnits::cm;

  const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );
  const shared_ptr<const Material> steel = matdb->material( "Stainless steel SS-304" );
  const shared_ptr<const Material> water = matdb->material( "Water" );
  BOOST_REQUIRE( steel && water );

  // radii[i] is shell i's OUTER radius (cm); `src` is which of them emits.
  struct Case
  {
    const char *name;
    std::vector<std::pair<shared_ptr<const Material>,double>> shells;
    size_t src;
    double dist_cm;
  };
  const std::vector<Case> cases = {
    { "solid water, contact",        { {water,3.0} },                        0, 4.0 },
    { "solid water, far",            { {water,3.0} },                        0, 50.0 },
    { "solid steel, contact",        { {steel,3.0} },                        0, 4.0 },
    { "hollow water on steel core",  { {steel,1.0}, {water,3.0} },           1, 4.0 },
    { "hollow steel on steel core",  { {steel,1.0}, {steel,2.5} },           1, 4.0 },
    { "solid water in steel shield", { {water,2.5}, {steel,2.9} },           0, 4.0 },
    { "hollow water, shielded",      { {steel,1.0}, {water,2.5}, {steel,2.9} }, 1, 4.0 },
  };

  const int num_lines = 1 << 16;
  double worst = 0.0;
  string worst_where;

  for( const Case &c : cases )
  {
    for( const double e : { 60.0, 661.7 } )
    {
      DistributedSrcCalcT<double> calc;
      calc.m_geometry = GeometryType::Spherical;
      calc.m_materialIndex = c.src;
      calc.m_attenuateForAir = false;
      calc.m_isInSituExponential = false;
      calc.m_inSituRelaxationLength = -1.0;
      calc.m_srcVolumetricActivity = 1.0;
      calc.m_normalizeByVolume = false;
      calc.m_energy = e;
      calc.m_effResponse = det.mc_transfer;
      calc.m_effMethod = ShieldingSourceFitCalc::VolumetricEffMethod::MCTransfer;
      calc.m_detector = detector_geom_from_config<double>( GeometryType::Spherical, c.dist_cm*cm,
                                              det.gd.transverse_half_extent()*cm, 0.0 );
      for( const std::pair<shared_ptr<const Material>,double> &sh : c.shells )
      {
        DistributedSrcCalcT<double>::ShellInfo info;
        info.dims = { sh.second*cm, 0.0, 0.0 };
        info.trans_len_coef = transmition_length_coefficient( sh.first.get(), static_cast<float>(e) );
        info.type = ShellType::Material;
        calc.m_shells.push_back( info );
      }

      DistributedSrcCalcT<double> elem = calc, line = calc;
      integrate_on_path( elem, VolumetricIntegrator::Element, -1 );
      integrate_on_path( line, VolumetricIntegrator::Line, num_lines );
      BOOST_REQUIRE( (elem.integral > 0.0) && (line.integral > 0.0) );

      const double rel = 100.0*(line.integral/elem.integral - 1.0);
      std::ostringstream row;
      row << "    " << std::left << std::setw(30) << c.name << std::right << " @ " << std::setw(6)
          << std::fixed << std::setprecision(1) << e << " keV:  element " << std::scientific
          << std::setprecision(4) << elem.integral << "  line " << line.integral << "  ("
          << std::fixed << std::showpos << std::setprecision(2) << rel << "%" << std::noshowpos << ")";
      BOOST_TEST_MESSAGE( row.str() );
      if( std::fabs(rel) > worst )
      {
        worst = std::fabs( rel );
        worst_where = string(c.name) + " @ " + std::to_string(e) + " keV";
      }
    }//for( energies )
  }//for( cases )

  BOOST_TEST_MESSAGE( "  worst |line/element - 1|: " << std::fixed << std::setprecision(3) << worst
                      << "% (" << worst_where << ")" );
  // The same gate the other geometries carry (LineVsElementScenarioMatrix uses 0.75%), loosened to
  //  1% because a sphere's outer quadrature is 2-D and its adaptive refinement is coarser.
  BOOST_CHECK_MESSAGE( worst < 1.0,
                       "line and element disagree by " << worst << "% at " << worst_where
                       << " - spherical sources are not being integrated the same way" );
}//BOOST_AUTO_TEST_CASE( LineVsElementSphericalSource )


/** A vanishingly small volumetric source must reproduce the POINT query at its centre.

 The point-source path evaluates a ray fan from the point (`eps_fep`); the volumetric path
 integrates the detector-side line set over the source volume.  They are two quadratures of the same
 CeeLo kernel, and for a transparent sphere of radius 1e-3 (and 1e-5) of the standoff the volume
 average of the efficiency equals its value at the centre to O(1e-6) - so all three must agree:
 the point query, the element path (adaptive, within 1e-3) and the line path (65536 lines, within its
 own ~0.2% quadrature noise).

 The descriptor's `reference_point` is then flipped and the response rebuilt around the SAME
 physical anchor position: InterSpec measures every distance from the detector face and forms the
 query positions itself, so nothing on the InterSpec side may change - the field is an internal
 CeeLo convention InterSpec never consults.
 */
BOOST_AUTO_TEST_CASE( TinySourceMatchesPointQuery )
{
  using namespace GammaInteractionCalc;
  const double cm = PhysicalUnits::cm;

  // A 3"x3" NaI with a thin Al can, and a plausible decreasing anchor curve (no data files).
  ceelo::GeometryDescriptor gd;
  gd.shape = ceelo::DetectorShape::Cylinder;
  gd.dimensions_cm = { 3.81, 7.62 };
  gd.materials = { ceelo::MaterialSpec::from( ceelo::make_NaI() ),
                   ceelo::MaterialSpec::from( ceelo::make_Aluminum() ) };
  gd.crystal_material_index = 0;
  ceelo::LayerSpec can;
  can.material_index = 1;
  can.front_thickness_cm = 0.05;
  can.side_thickness_cm = 0.05;
  can.z_end_cm = 7.62;
  gd.layers.push_back( can );
  gd.reference_point = ceelo::ReferencePoint::EndcapFront;

  ceelo::AnchorCurve anchor;
  anchor.energies_keV = { 60.0, 100.0, 300.0, 662.0, 1000.0 };
  anchor.eff = { 2.0e-2, 1.3e-2, 5.0e-3, 3.0e-3, 2.2e-3 };
  anchor.frac_sigma = { 0.003, 0.003, 0.003, 0.003, 0.003 };
  const Eigen::Vector3d anchor_pos( 0.0, 0.0, -25.0 );   //crystal-face frame: a POSITION, not a distance

  ceelo::GeometryDescriptor gd_crystal = gd;
  gd_crystal.reference_point = ceelo::ReferencePoint::CrystalFace;

  const std::shared_ptr<const ceelo::DetectorResponse> resp_endcap
        = ceelo::make_transfer_response( gd, anchor, anchor_pos );
  const std::shared_ptr<const ceelo::DetectorResponse> resp_crystal
        = ceelo::make_transfer_response( gd_crystal, anchor, anchor_pos );
  BOOST_REQUIRE( resp_endcap && resp_crystal );

  const double dist_cm = 20.0;
  const int num_lines = 1 << 16;

  // Face-referenced point query: the distance InterSpec means, through the EndcapFront descriptor.
  const auto point_query = [&]( const double energy ) -> double {
    return resp_endcap->eps_fep( energy, 0.0, 0.0, dist_cm ).value;
  };

  const auto tiny_calc = [&]( const std::shared_ptr<const ceelo::DetectorResponse> &resp,
                              const double radius_cm, const double energy ) -> DistributedSrcCalcT<double>
  {
    DistributedSrcCalcT<double> calc;
    calc.m_geometry = GeometryType::Spherical;
    calc.m_materialIndex = 0;
    calc.m_attenuateForAir = false;
    calc.m_isInSituExponential = false;
    calc.m_inSituRelaxationLength = -1.0;
    calc.m_srcVolumetricActivity = 1.0;
    calc.m_normalizeByVolume = true;      //integral = volume-average efficiency
    calc.m_energy = energy;
    calc.m_effResponse = resp;
    calc.m_effMethod = ShieldingSourceFitCalc::VolumetricEffMethod::MCTransfer;
    calc.m_detector = detector_geom_from_config<double>( GeometryType::Spherical, dist_cm*cm,
                                                         gd.transverse_half_extent()*cm, 0.0 );
    DistributedSrcCalcT<double>::ShellInfo info;
    info.dims = { radius_cm*cm, 0.0, 0.0 };
    info.trans_len_coef = 0.0;            //transparent
    info.type = ShellType::Material;
    calc.m_shells.push_back( info );
    return calc;
  };

  for( const double ratio : { 1.0e-3, 1.0e-5 } )
  {
    for( const double e : { 60.0, 661.7 } )
    {
      const double point = point_query( e );
      BOOST_REQUIRE( point > 0.0 );

      DistributedSrcCalcT<double> elem = tiny_calc( resp_endcap, ratio*dist_cm, e );
      DistributedSrcCalcT<double> line = tiny_calc( resp_endcap, ratio*dist_cm, e );
      integrate_on_path( elem, VolumetricIntegrator::Element, -1 );
      integrate_on_path( line, VolumetricIntegrator::Line, num_lines );

      const double elem_rel = 100.0*(elem.integral/point - 1.0);
      const double line_rel = 100.0*(line.integral/point - 1.0);
      std::ostringstream row;
      row << "    r/d=" << std::scientific << std::setprecision(0) << ratio << " @ " << std::fixed
          << std::setprecision(1) << e << " keV: point " << std::scientific << std::setprecision(5)
          << point << "  element " << elem.integral << " (" << std::fixed << std::showpos
          << std::setprecision(3) << elem_rel << "%)  line " << std::scientific << std::setprecision(5)
          << line.integral << " (" << std::fixed << std::showpos << std::setprecision(3) << line_rel
          << "%)" << std::noshowpos;
      BOOST_TEST_MESSAGE( row.str() );
      BOOST_CHECK_MESSAGE( std::fabs(elem_rel) < 0.1,
                           "element path on a tiny source differs from the point query by " << elem_rel << "%" );
      BOOST_CHECK_MESSAGE( std::fabs(line_rel) < 0.5,
                           "line path on a tiny source differs from the point query by " << line_rel << "%" );

      // The reference-point field must be inert on the InterSpec side.
      DistributedSrcCalcT<double> elem_c = tiny_calc( resp_crystal, ratio*dist_cm, e );
      DistributedSrcCalcT<double> line_c = tiny_calc( resp_crystal, ratio*dist_cm, e );
      integrate_on_path( elem_c, VolumetricIntegrator::Element, -1 );
      integrate_on_path( line_c, VolumetricIntegrator::Line, num_lines );
      BOOST_CHECK_MESSAGE( std::fabs(elem_c.integral/elem.integral - 1.0) < 1.0e-9,
                           "element path depends on the descriptor's reference_point: "
                           << elem_c.integral << " vs " << elem.integral );
      BOOST_CHECK_MESSAGE( std::fabs(line_c.integral/line.integral - 1.0) < 1.0e-9,
                           "line path depends on the descriptor's reference_point: "
                           << line_c.integral << " vs " << line.integral );
    }//for( energies )
  }//for( ratio )
}//BOOST_AUTO_TEST_CASE( TinySourceMatchesPointQuery )


/** DEVELOPER PROBE (2026-09-05): a deeply opaque self-attenuating sphere on the FLAT-DISK path.

 Motivated by a real fit that moved between app builds: 18 cm of enriched uranium at 1 m, where the
 escaping signal comes from a skin of 1/mu.  At 1001 keV that skin is ~7 mm (4% of the radius) and
 the volume quadrature resolves it; at 121 keV it is ~0.1 mm (7e-6 of the radius) and no globally
 adaptive rule can find it within its evaluation budget.  Prints, per energy, the coefficient the
 model uses, what the quadrature returns, and the analytic deep-opacity limit 3/(4 mu R) that the
 volume-averaged escape probability must approach - so an unconverged row is visible as a departure
 from that limit rather than as a plausible-looking number.
 */
BOOST_AUTO_TEST_CASE( OpaqueSphereSelfAttenConvergence, * boost::unit_test::disabled() )
{
  using namespace GammaInteractionCalc;
  const double cm = PhysicalUnits::cm;
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );
  const std::shared_ptr<const Material> mat = matdb->material( "Enriched uranium alloy" );
  BOOST_REQUIRE( mat );

  struct Row { double energy, fwhm; };
  const std::vector<Row> rows = {
    {120.91,0.626},{143.76,0.661},{163.33,0.698},{185.68,0.726},{205.27,0.755},{238.63,0.844},
    {258.18,0.854},{583.21,1.087},{742.89,1.396},{766.48,1.445},{1001.09,1.713} };

  const AngleDetector det = load_angle_detector();
  const std::shared_ptr<const ceelo::DetectorResponse> resp = det.mc_transfer;
  BOOST_REQUIRE( resp );

  const double dist = 100.0*cm, det_rad = 3.2*cm;
  const double omega = DetectorPeakResponse::fractionalSolidAngle( 2.0*det_rad, dist );

  std::ostringstream hdr;
  hdr << "  density " << mat->density/(PhysicalUnits::g/PhysicalUnits::cm3) << " g/cm3, centre solid angle "
      << std::scientific << std::setprecision(5) << omega;
  BOOST_TEST_MESSAGE( hdr.str() );
  BOOST_TEST_MESSAGE( "   E(keV)   mu/rho   mu_fep/mu |   R(cm)   quadrature      analytic 3/(4muR)   quad/analytic" );

  for( const double radius_cm : { 18.01006 } )
  {
    for( const Row &r : rows )
    {
      const double mu = transmition_length_coefficient( mat.get(), static_cast<float>(r.energy) );
      const double mu_fep = fep_survival_removal_coefficient( mat.get(), static_cast<float>(r.energy),
                                                              0.5*r.fwhm, 0.0 );
      DistributedSrcCalcT<double> calc;
      calc.m_geometry = GeometryType::Spherical;
      calc.m_materialIndex = 0;
      calc.m_attenuateForAir = false;
      calc.m_isInSituExponential = false;
      calc.m_inSituRelaxationLength = -1.0;
      calc.m_srcVolumetricActivity = 1.0;
      calc.m_normalizeByVolume = false;
      calc.m_energy = r.energy;
      calc.m_detector = detector_geom_from_config<double>( GeometryType::Spherical, dist, det_rad, 0.0 );
      DistributedSrcCalcT<double>::ShellInfo info;
      info.dims = { radius_cm*cm, 0.0, 0.0 };
      info.trans_len_coef = mu;
      info.fep_trans_len_coef = mu_fep;
      info.type = ShellType::Material;
      info.density = mat->density;
      calc.m_shells.push_back( info );

      integrate_on_path( calc, VolumetricIntegrator::Element, -1 );

      // The same source through the LINE path, whose chord integral of exp(-mu s) is analytic and
      //  therefore cannot miss a thin skin.  It needs a response, so this leg is a RATIO test of
      //  the two quadratures on the identical integrand, not of absolute values.
      double line_over_elem = -1.0;
      {
        DistributedSrcCalcT<double> e2 = calc, l2 = calc;
        e2.m_effResponse = resp;  l2.m_effResponse = resp;
        e2.m_effMethod = ShieldingSourceFitCalc::VolumetricEffMethod::MCTransfer;
        l2.m_effMethod = ShieldingSourceFitCalc::VolumetricEffMethod::MCTransfer;
        integrate_on_path( e2, VolumetricIntegrator::Element, -1 );
        integrate_on_path( l2, VolumetricIntegrator::Line, 1 << 16 );
        if( e2.integral > 0.0 )
          line_over_elem = l2.integral / e2.integral;
      }

      const double volume = (4.0/3.0)*M_PI*std::pow( radius_cm*cm, 3.0 );
      const double quad = (calc.integral/volume)/omega;              //the reports "Shield Atten. Factor"
      const double analytic = 3.0/(4.0*mu_fep*radius_cm*cm);         //deep-opacity limit of <exp(-mu t)>

      std::ostringstream o;
      o << "  " << std::fixed << std::setw(8) << std::setprecision(2) << r.energy
        << std::setw(9) << std::setprecision(4) << mu/(mat->density*PhysicalUnits::cm2/PhysicalUnits::g)
        << std::setw(11) << std::setprecision(5) << (mu_fep/mu)
        << " | " << std::setw(8) << std::setprecision(3) << radius_cm
        << "  " << std::scientific << std::setprecision(4) << quad
        << "      " << analytic
        << "      " << std::fixed << std::setprecision(3) << (quad/analytic)
        << "   line/elem " << std::setprecision(3) << line_over_elem;
      BOOST_TEST_MESSAGE( o.str() );
    }
    BOOST_TEST_MESSAGE( "" );
  }
}//BOOST_AUTO_TEST_CASE( OpaqueSphereSelfAttenConvergence )
