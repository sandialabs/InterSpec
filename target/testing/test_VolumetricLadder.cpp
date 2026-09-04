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


/** The incremental validation LADDER for InterSpec's volumetric-source per-ray efficiency kernel.

 test_VolumetricNearField.cpp compares whole volume integrals against a Monte-Carlo truth bank and
 has accumulated many one-off diagnostics.  This file climbs from the smallest verifiable piece to the
 whole, one rung at a time, so that WHERE an error first appears is a measurement rather than an
 inference:

   rung 0  no Monte Carlo: fan geometry audit, per-ray optical-depth audit, and a zero-MC
           before/after of the aperture reference-point fix over the existing truth bank;
   rung 1  the transfer at bare points on the emitting skin (MC, 0.25%);
   rung 2  transparent whole volumes (MC);
   rung 3  ONE attenuated point at depth in a slab (MC) - the rung between points and volumes that
           the older suite never had;
   rung 4  slabs of increasing optical depth, then the scenario matrix (MC).

 MC-running rungs are developer-only (disabled) and go through the on-disk cache in
 VolumetricNearFieldHarness.h (`--cachedir=`), so the InterSpec side can be re-evaluated for free
 after any code change.  Findings are recorded in scratch/20260902_volumetric_ladder/FINDINGS.md.
 */

#include <array>
#include <cmath>
#include <string>
#include <map>
#include <vector>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <iostream>
#include <set>
#include <algorithm>

#define BOOST_TEST_MODULE VolumetricLadder_suite
#include <boost/test/included/unit_test.hpp>

#include "ceres/jet.h"

#include "VolumetricNearFieldHarness.h"


namespace
{
/** Assembly-frame description of the crystal solid, for the ray-hits-crystal audit: a cylinder of
 radius `radius_cm` from `z_face` to `z_face + length_cm` along +z, centred on the detector axis.
 Uses the CRYSTAL SOLID (dimensions_cm), not the dead-layer-inflated transverse_half_extent(). */
struct CrystalInAssembly
{
  double radius_cm = 0.0;
  double length_cm = 0.0;
  double z_face = 0.0;     ///< crystal-face plane, PhysicalUnits-free (cm), assembly frame
  double x0 = 0.0, y0 = 0.0;   ///< axis position (detector lateral offset)
};

CrystalInAssembly crystal_in_assembly( const AngleDetector &det,
                                        const GammaInteractionCalc::DetectorGeomT<double> &dg )
{
  CrystalInAssembly c;
  BOOST_REQUIRE( det.gd.shape == ceelo::DetectorShape::Cylinder );
  BOOST_REQUIRE( det.gd.dimensions_cm.size() >= 2 );
  c.radius_cm = det.gd.dimensions_cm[0];
  c.length_cm = det.gd.dimensions_cm[1];
  // det.axis points detector -> assembly (0,0,-1); the crystal face is `off` BEHIND the endcap front,
  //  i.e. further along -axis.
  const double cm = PhysicalUnits::cm;
  c.z_face = dg.position[2]/cm - dg.axis[2]*det.endcap_front_offset_cm;
  c.x0 = dg.position[0]/cm;
  c.y0 = dg.position[1]/cm;
  return c;
}

/** Does the ray from `e` (cm) along unit `d` enter the crystal solid? */
bool ray_hits_crystal( const CrystalInAssembly &c, const double e[3], const std::array<double,3> &d )
{
  if( !(d[2] > 0.0) )
    return false;
  const double ta = (c.z_face - e[2]) / d[2];
  const double tb = (c.z_face + c.length_cm - e[2]) / d[2];
  if( tb < 0.0 )
    return false;

  // rho^2(t) along the ray, minimised over [max(ta,0), tb]
  const double ex = e[0] - c.x0, ey = e[1] - c.y0;
  const auto rho2 = [&]( const double t ) {
    const double x = ex + t*d[0], y = ey + t*d[1];
    return x*x + y*y;
  };
  const double t0 = std::max( ta, 0.0 );
  double best = std::min( rho2(t0), rho2(tb) );
  const double dd = d[0]*d[0] + d[1]*d[1];
  if( dd > 0.0 )
  {
    const double tv = -(ex*d[0] + ey*d[1]) / dd;
    if( (tv > t0) && (tv < tb) )
      best = std::min( best, rho2(tv) );
  }
  return best <= c.radius_cm*c.radius_cm;
}

/** One element of the fan audit. */
struct FanProbe
{
  string label;
  double e_cm[3];   ///< element position, assembly frame, cm
};

/** Runs build_element_aperture for one element exactly as eval_rect does, under the current
 #sm_aperture_frame_legacy_origin, and returns (tilt_deg, miss_fraction, weighted <d.z>). */
struct FanMeasure
{
  double tilt_deg = 0.0;
  double miss_frac = 0.0;
  double mean_dz = 0.0;
  double predicted_tilt_deg = 0.0;   ///< theta - theta' from the endcap-offset parallax
  size_t n_rays = 0;
};

FanMeasure measure_fan( const AngleDetector &det,
                        const GammaInteractionCalc::DetectorGeomT<double> &dg,
                        const FanProbe &p, const double energy_keV, const int n_rays )
{
  using namespace GammaInteractionCalc;
  const double cm = PhysicalUnits::cm;

  const double e[3] = { p.e_cm[0]*cm, p.e_cm[1]*cm, p.e_cm[2]*cm };
  double to_det[3] = { dg.position[0]-e[0], dg.position[1]-e[1], dg.position[2]-e[2] };
  const double dist = sqrt( to_det[0]*to_det[0] + to_det[1]*to_det[1] + to_det[2]*to_det[2] );
  BOOST_REQUIRE( dist > 0.0 );
  for( int k = 0; k < 3; ++k )
    to_det[k] /= dist;
  double cos_theta = -(to_det[0]*dg.axis[0] + to_det[1]*dg.axis[1] + to_det[2]*dg.axis[2]);
  cos_theta = std::max( -1.0, std::min( 1.0, cos_theta ) );

  const ElementAperture ap = build_element_aperture( *det.mc_transfer, energy_keV, dist, cos_theta,
                                                     to_det, dg.axis, n_rays );

  FanMeasure m;
  m.n_rays = ap.dirs.size();

  // Where did CeeLo's detector axis a_c = (0,0,-1) land?
  const double a_c[3] = { 0.0, 0.0, -1.0 };
  double ra[3] = { 0.0, 0.0, 0.0 };
  for( int i = 0; i < 3; ++i )
    for( int j = 0; j < 3; ++j )
      ra[i] += ap.rotation[i][j] * a_c[j];
  const double dot = std::max( -1.0, std::min( 1.0, ra[0]*dg.axis[0] + ra[1]*dg.axis[1] + ra[2]*dg.axis[2] ) );
  m.tilt_deg = acos( dot ) * 180.0 / M_PI;

  const CrystalInAssembly cr = crystal_in_assembly( det, dg );
  double sum_w = 0.0, sum_wdz = 0.0, sum_w_miss = 0.0;
  for( size_t i = 0; i < ap.dirs.size(); ++i )
  {
    const double w = ap.weights[i];
    sum_w += w;
    sum_wdz += w * ap.dirs[i][2];
    if( !ray_hits_crystal( cr, p.e_cm, ap.dirs[i] ) )
      sum_w_miss += w;
  }
  m.miss_frac = (sum_w > 0.0) ? (sum_w_miss/sum_w) : 0.0;
  m.mean_dz = (sum_w > 0.0) ? (sum_wdz/sum_w) : 0.0;

  // The parallax prediction: polar angle to the endcap-front point vs to the crystal-face centre.
  const double r = sqrt( pow(p.e_cm[0]-cr.x0,2) + pow(p.e_cm[1]-cr.y0,2) );
  const double z_e = dg.position[2]/cm - p.e_cm[2];
  m.predicted_tilt_deg = (atan2( r, z_e ) - atan2( r, z_e + det.endcap_front_offset_cm )) * 180.0/M_PI;

  return m;
}

}//namespace


/** RUNG 0a - fan geometry audit.

 For elements at contact and in the far field, on and off axis, front face and side face, measure
 what build_element_aperture actually produced in the ASSEMBLY frame:

   tilt      angle between the image of CeeLo's detector axis and InterSpec's detector axis.  A
             correct frame mapping is a pure azimuthal rotation, so this must be 0.
   miss      weight fraction of aperture rays that do NOT enter the crystal solid.  CeeLo only
             returns rays with a positive active chord, so every one of them must hit.
   <d.z>     weighted mean direction cosine to the detector axis - for a front-face skin element
             this is the quantity a thick-source integral is proportional to.

 Measured under BOTH settings of #sm_aperture_frame_legacy_origin; the fixed setting is gated.  The
 legacy setting is predicted to show tilt == theta - theta' (the endcap-offset parallax), which is
 also printed.
 */
BOOST_AUTO_TEST_CASE( FanGeometryAudit )
{
  using namespace GammaInteractionCalc;
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );
  const double cm = PhysicalUnits::cm;

  BOOST_TEST_MESSAGE( "endcap_front_offset = " << det.endcap_front_offset_cm << " cm; crystal r="
                      << det.gd.dimensions_cm[0] << " len=" << det.gd.dimensions_cm[1] << " cm" );

  // Probes are laid out for the large-near-* geometry (r=4, hl=2, standoff 1): front-face skin,
  //  side-face skin, and the box corner; then the same at 50 cm standoff.
  const double hl = 2.0, R = 4.0;
  const vector<FanProbe> probes = {
    { "front r=0",      { 0.0,   0.0, hl-0.05 } },
    { "front r=1e-6",   { 1e-6,  0.0, hl-0.05 } },
    { "front r=1",      { 1.0,   0.0, hl-0.05 } },
    { "front r=2",      { 2.0,   0.0, hl-0.05 } },
    { "front r=3",      { 3.0,   0.0, hl-0.05 } },
    { "front r=4",      { 4.0,   0.0, hl-0.05 } },
    { "front corner",   { 4.0,   4.0, hl-0.05 } },
    { "side z=+1.5",    { R-0.05, 0.0,  1.5 } },
    { "side z=0",       { R-0.05, 0.0,  0.0 } },
    { "side z=-1.5",    { R-0.05, 0.0, -1.5 } },
  };

  const bool saved = sm_aperture_frame_legacy_origin;
  double worst_tilt_fixed = 0.0, worst_miss_fixed = 0.0;

  for( const double standoff : { 1.0, 50.0 } )
  {
    const DetectorGeomT<double> dg = detector_geom_from_config<double>( GeometryType::CylinderEndOn,
                (standoff + hl)*cm, det.gd.transverse_half_extent()*cm, 0.0 );

    for( const double energy : { 60.0, 1332.5 } )
    {
      ostringstream hdr;
      hdr << "standoff " << standoff << " cm, " << energy << " keV  (512 rays)\n"
          << "  element         | legacy: tilt   miss   <dz>  | fixed: tilt   miss   <dz>  | predicted tilt";
      BOOST_TEST_MESSAGE( hdr.str() );

      for( const FanProbe &p : probes )
      {
        sm_aperture_frame_legacy_origin = true;
        const FanMeasure leg = measure_fan( det, dg, p, energy, 512 );
        sm_aperture_frame_legacy_origin = false;
        const FanMeasure fix = measure_fan( det, dg, p, energy, 512 );

        worst_tilt_fixed = std::max( worst_tilt_fixed, fix.tilt_deg );
        worst_miss_fixed = std::max( worst_miss_fixed, fix.miss_frac );

        ostringstream o;
        o << "  " << left << setw(16) << p.label << "| " << fixed
          << setprecision(2) << setw(6) << leg.tilt_deg << "  " << setprecision(3) << setw(5) << leg.miss_frac
          << "  " << setprecision(3) << setw(5) << leg.mean_dz << " | "
          << setprecision(2) << setw(6) << fix.tilt_deg << "  " << setprecision(3) << setw(5) << fix.miss_frac
          << "  " << setprecision(3) << setw(5) << fix.mean_dz << " | "
          << setprecision(2) << leg.predicted_tilt_deg << " deg";
        BOOST_TEST_MESSAGE( o.str() );
      }//for( probes )
    }//for( energies )
  }//for( standoff )

  sm_aperture_frame_legacy_origin = saved;

  BOOST_TEST_MESSAGE( "fixed frame: worst tilt " << worst_tilt_fixed << " deg, worst miss fraction "
                      << worst_miss_fixed );
  BOOST_CHECK_MESSAGE( worst_tilt_fixed < 0.05,
                       "Aperture frame still tilts CeeLo's detector axis by " << worst_tilt_fixed << " deg" );
  BOOST_CHECK_MESSAGE( worst_miss_fixed < 1.0e-9,
                       "Aperture rays miss the crystal (weight fraction " << worst_miss_fixed << ")" );
}//BOOST_AUTO_TEST_CASE( FanGeometryAudit )


/** RUNG 0d - zero-Monte-Carlo before/after of the frame fix over the EXISTING truth bank.

 Deterministic: the same integrand is run with the legacy and the fixed aperture frame against the
 recorded truth rows.  Prediction if the tilt is the near-field high-tau deficit: the fixed/legacy
 ratio is +several % on every optically thick contact row at every energy, ~0 on far-field and
 low-tau rows, and the fixed residual collapses to the transfer's own point bias.
 Minutes (box rows dominate); developer-only.
 */
BOOST_AUTO_TEST_CASE( LegacyVsFixedFrameOverTruthBank, * boost::unit_test::disabled() )
{
  using namespace GammaInteractionCalc;
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );
  BOOST_REQUIRE( !sm_truth.empty() );

  const bool saved = sm_aperture_frame_legacy_origin;

  BOOST_TEST_MESSAGE( "scenario                 E(keV)    tau   legacy/truth-1  fixed/truth-1  fixed/legacy-1" );
  for( const Scenario &s : scenarios() )
  {
    for( const TruthRow &row : sm_truth )
    {
      if( (row.scenario != s.name) || !(row.fep_eff > 0.0) )
        continue;

      sm_aperture_frame_legacy_origin = true;
      const double legacy = interspec_volumetric_eff( det, s, row.energy_keV, det.mc_transfer );
      sm_aperture_frame_legacy_origin = false;
      const double fixed_v = interspec_volumetric_eff( det, s, row.energy_keV, det.mc_transfer );

      ostringstream o;
      o << left << setw(25) << s.name << right << fixed << setprecision(1) << setw(7) << row.energy_keV
        << setprecision(2) << setw(7) << scenario_optical_depth( s, row.energy_keV )
        << setprecision(2) << setw(13) << 100.0*(legacy/row.fep_eff - 1.0) << "%"
        << setw(14) << 100.0*(fixed_v/row.fep_eff - 1.0) << "%"
        << setw(15) << 100.0*(fixed_v/legacy - 1.0) << "%"
        << "   (mc sigma " << setprecision(2) << 100.0*row.fep_uncert/row.fep_eff << "%)";
      BOOST_TEST_MESSAGE( o.str() );
    }
  }

  sm_aperture_frame_legacy_origin = saved;
}//BOOST_AUTO_TEST_CASE( LegacyVsFixedFrameOverTruthBank )


namespace
{
/** "A point at depth": a tiny rectangular source (half-extent `sm_point_half_cm`) at the centre of a
 slab of the same material, `depth_cm` below the slab's near face, laterally `lateral_cm` off the
 detector axis, the near face `standoff_cm` from the endcap front.  Built identically on both sides:
   InterSpec: Rectangular tiny source shell + material shell of the same outer dims;
   CeeLo:     set_rectangular_source(tiny) + set_source_material + add_source_shield(hx, hy, depth).
 The slab is laterally wide (`hx_cm`) so lateral exits never matter.
 */
constexpr double sm_point_half_cm = 1.0e-3;

struct PointAtDepth
{
  double depth_cm = 0.1;
  double lateral_cm = 0.0;
  double standoff_cm = 1.0;
  bool dense = true;
  double hx_cm = 10.0;
};

double point_center_distance_cm( const PointAtDepth &p )
{
  return p.standoff_cm + p.depth_cm + sm_point_half_cm;
}

string point_at_depth_key( const PointAtDepth &p )
{
  ostringstream o;
  o << setprecision(10) << "pointInSlab(depth=" << p.depth_cm << ",lateral=" << p.lateral_cm
    << ",standoff=" << p.standoff_cm << ",hx=" << p.hx_cm << ",eps=" << sm_point_half_cm
    << ";matrix=" << scenario_matrix_material( p.dense ) << ")";
  return o.str();
}

GammaInteractionCalc::DistributedSrcCalcT<double>
build_point_at_depth_calc( const AngleDetector &det, const PointAtDepth &p, const double energy_keV,
                           const shared_ptr<const ceelo::DetectorResponse> &resp,
                           const int n_rays = -1, const bool transparent = false,
                           const double fep_window_keV = -1.0 )
{
  using namespace GammaInteractionCalc;
  const double cm = PhysicalUnits::cm;
  const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );
  const shared_ptr<const Material> mat = matdb->material( scenario_matrix_material( p.dense ) );
  BOOST_REQUIRE( mat );
  // `fep_window_keV > 0` uses the SURVIVAL-removal coefficient (mu_total minus the forward-Compton
  //  scatters whose energy loss stays inside the +-window), which is what the peak actually sees.
  const double mu = transparent ? 0.0
        : ( (fep_window_keV > 0.0)
              ? fep_removal_coefficient( *mat, energy_keV, fep_window_keV )
              : transmition_length_coefficient( mat.get(), static_cast<float>(energy_keV) ) );

  DistributedSrcCalcT<double> calc;
  calc.m_geometry = GeometryType::Rectangular;
  calc.m_materialIndex = 0;
  calc.m_attenuateForAir = false;
  calc.m_airTransLenCoef = 0.0;
  calc.m_isInSituExponential = false;
  calc.m_inSituRelaxationLength = -1.0;
  calc.m_srcVolumetricActivity = 1.0;
  calc.m_normalizeByVolume = false;
  calc.m_energy = energy_keV;
  calc.m_nuclide = nullptr;
  calc.integral = 0.0;
  calc.m_effResponse = resp;
  calc.m_effMethod = resp ? ShieldingSourceFitCalc::VolumetricEffMethod::MCTransfer
                          : ShieldingSourceFitCalc::VolumetricEffMethod::FlatDisk;
  if( n_rays > 0 )
    calc.m_effNumRays = n_rays;

  calc.m_detector = detector_geom_from_config<double>( GeometryType::Rectangular,
                        point_center_distance_cm( p )*cm, det.gd.transverse_half_extent()*cm, 0.0,
                        p.lateral_cm*cm, 0.0 );

  DistributedSrcCalcT<double>::ShellInfo src;
  src.dims = { sm_point_half_cm*cm, sm_point_half_cm*cm, sm_point_half_cm*cm };
  src.trans_len_coef = mu;
  src.type = ShellType::Material;
  calc.m_shells.push_back( src );

  DistributedSrcCalcT<double>::ShellInfo slab;
  slab.dims = { (p.hx_cm + sm_point_half_cm)*cm, (p.hx_cm + sm_point_half_cm)*cm,
                (p.depth_cm + sm_point_half_cm)*cm };
  slab.trans_len_coef = mu;
  slab.type = ShellType::Material;
  calc.m_shells.push_back( slab );

  return calc;
}//build_point_at_depth_calc(...)

void configure_point_at_depth_mc( ceelo::EfficiencyCalculator &calc, const AngleDetector &det,
                                  const PointAtDepth &p, ScenarioMcMaterials &mats )
{
  const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );
  ceelo::ResponseGenerator::configure_calculator( calc, det.gd, mats.owned );

  // CeeLo crystal-face frame: the detector's lateral offset in InterSpec is the source's lateral
  //  offset here (mirror-symmetric, sign irrelevant).
  const Eigen::Vector3d centre( p.lateral_cm, 0.0,
                                -(det.endcap_front_offset_cm + point_center_distance_cm( p )) );
  calc.set_rectangular_source( centre, Eigen::Vector3d( sm_point_half_cm, sm_point_half_cm, sm_point_half_cm ) );

  const shared_ptr<const Material> mat = matdb->material( scenario_matrix_material( p.dense ) );
  BOOST_REQUIRE( mat );
  mats.matrix.reset( new ceelo::Material( CeeLoUtils::to_ceelo_material( *mat ).to_material() ) );
  calc.set_source_material( mats.matrix.get() );
  calc.add_source_shield( mats.matrix.get(), p.hx_cm, p.hx_cm, p.depth_cm );
}//configure_point_at_depth_mc(...)


/** Independent optical depth from `e` along unit `d` through the calculator's shells, by a brute
 force march with point-in-shell material lookup (innermost containing shell wins).  Uses the shells'
 FEP coefficients, so call finalize_shell_coefficients() first. */
double march_optical_depth( const GammaInteractionCalc::DistributedSrcCalcT<double> &calc,
                            const double e[3], const std::array<double,3> &d, const double step )
{
  using namespace GammaInteractionCalc;
  const auto inside = [&]( const size_t i, const double x, const double y, const double z ) -> bool {
    const std::array<double,3> &dm = calc.m_shells[i].dims;
    if( calc.m_geometry == GeometryType::Rectangular )
      return (fabs(x) <= dm[0]) && (fabs(y) <= dm[1]) && (fabs(z) <= dm[2]);
    BOOST_REQUIRE( calc.m_geometry == GeometryType::CylinderEndOn );
    return ((x*x + y*y) <= dm[0]*dm[0]) && (fabs(z) <= dm[1]);
  };

  // Farthest any shell reaches, so the march knows when to stop.
  double extent = 0.0;
  for( const auto &sh : calc.m_shells )
    extent = std::max( extent, sqrt( sh.dims[0]*sh.dims[0] + sh.dims[1]*sh.dims[1] + sh.dims[2]*sh.dims[2] ) );
  const double s_max = 2.0*extent + sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );

  double tau = 0.0;
  for( double s0 = 0.0; s0 < s_max; s0 += step )
  {
    const double sm = s0 + 0.5*step;
    const double x = e[0] + sm*d[0], y = e[1] + sm*d[1], z = e[2] + sm*d[2];
    for( size_t i = 0; i < calc.m_shells.size(); ++i )
    {
      if( inside( i, x, y, z ) )
      {
        tau += calc.m_shells[i].fep_trans_len_coef * step;
        break;
      }
    }
  }
  return tau;
}//march_optical_depth(...)


/** The test-side per-point kernel: prefactor * sum_i w_i * exp(-tau_march,i), with the aperture from
 build_element_aperture (as production uses it) but the optical depths from #march_optical_depth. */
double marched_point_kernel( const AngleDetector &det,
                             const GammaInteractionCalc::DistributedSrcCalcT<double> &calc,
                             const double e[3], const double step )
{
  using namespace GammaInteractionCalc;
  const DetectorGeomT<double> &dg = calc.m_detector;
  double to_det[3] = { dg.position[0]-e[0], dg.position[1]-e[1], dg.position[2]-e[2] };
  const double dist = sqrt( to_det[0]*to_det[0] + to_det[1]*to_det[1] + to_det[2]*to_det[2] );
  for( int k = 0; k < 3; ++k )
    to_det[k] /= dist;
  const double cos_theta = std::max( -1.0, std::min( 1.0,
        -(to_det[0]*dg.axis[0] + to_det[1]*dg.axis[1] + to_det[2]*dg.axis[2]) ) );

  const ElementAperture ap = build_element_aperture( *calc.m_effResponse, calc.m_energy, dist,
                                                     cos_theta, to_det, dg.axis, calc.m_effNumRays );
  double sum = 0.0;
  for( size_t i = 0; i < ap.dirs.size(); ++i )
    sum += ap.weights[i] * exp( -march_optical_depth( calc, e, ap.dirs[i], step ) );
  return ap.prefactor * sum;
}//marched_point_kernel(...)


/** An EFFTRAN transfer anchored on Monte-Carlo point sources at `distance_cm` on axis, through the
 cache, at the ladder precision.  `distance_cm` is endcap-front referenced. */
shared_ptr<ceelo::DetectorResponse> ladder_anchor_response( const AngleDetector &det, McCache &cache,
                                                            const double distance_cm,
                                                            const vector<double> &energies,
                                                            const double target_rel = 0.0025 )
{
  vector<unique_ptr<ceelo::Material>> owned;
  ceelo::EfficiencyCalculator pt;
  ceelo::ResponseGenerator::configure_calculator( pt, det.gd, owned );
  pt.set_point_source( Eigen::Vector3d( 0.0, 0.0, -(det.endcap_front_offset_cm + distance_cm) ) );

  ostringstream k;
  k << "point(r=0,ze=" << setprecision(10) << distance_cm << ");bare";

  CeeLoUtils::TransferAnchor anchor;
  anchor.ref_distance_cm = distance_cm;
  anchor.curve_derived = false;
  for( const double e : energies )
  {
    const McResult r = run_mc( pt, cache, k.str(), e, target_rel );
    anchor.curve.energies_keV.push_back( e );
    anchor.curve.eff.push_back( r.eff );
    anchor.curve.frac_sigma.push_back( r.frac_sigma );
  }

  ostringstream nm;
  nm << "ladder anchor d=" << fixed << setprecision(2) << distance_cm << " cm";
  const shared_ptr<ceelo::DetectorResponse> resp =
        CeeLoUtils::makeTransferResponse( det.gd, anchor, ceelo::AnchorCurve{}, nm.str() );
  BOOST_REQUIRE( resp );
  return resp;
}//ladder_anchor_response(...)


string pct( const double x )
{
  ostringstream o;
  o << fixed << setprecision(2) << setw(7) << 100.0*x << "%";
  return o.str();
}

}//namespace


/** RUNG 0b - per-ray kernel identity against an independent ray march (no MC).

 For a point at depth in a slab (rectangular walk) and a tiny cylinder inside a cylinder of the same
 material (cylindrical walk), the production per-element kernel must equal
 prefactor * sum_i w_i * exp(-tau_i) with tau_i from a brute-force march through the same shells.
 Checks the shell walk (chords, exit faces, multi-shell bookkeeping) and the per-ray composition
 together, at the exact element positions the point-at-depth rung uses.
 */
BOOST_AUTO_TEST_CASE( PerRayKernelIdentityVsRayMarch )
{
  using namespace GammaInteractionCalc;
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );
  const double cm = PhysicalUnits::cm;
  const double step = 2.0e-4*cm;

  double worst = 0.0;
  string worst_where;

  // Rectangular walk: the point-at-depth configurations.
  for( const bool dense : { true, false } )
  {
    for( const double depth : { 0.05, 0.5 } )
    {
      for( const double lateral : { 0.0, 2.0, 4.0 } )
      {
        for( const double standoff : { 1.0, 50.0 } )
        {
          PointAtDepth p;
          p.depth_cm = depth;  p.lateral_cm = lateral;  p.standoff_cm = standoff;  p.dense = dense;
          const double energy = 60.0;
          DistributedSrcCalcT<double> calc = build_point_at_depth_calc( det, p, energy, det.mc_transfer, 128 );
          const double prod = point_kernel_eff( calc );   //finalizes coefficients
          const double e[3] = { 0.0, 0.0, 0.0 };
          const double ref = marched_point_kernel( det, calc, e, step );
          const double rel = fabs( prod/ref - 1.0 );
          ostringstream o;
          o << "  rect " << (dense ? "steel" : "water") << " depth=" << depth << " lateral=" << lateral
            << " standoff=" << standoff << ": production " << setprecision(8) << prod
            << " march " << ref << "  rel " << setprecision(3) << rel;
          BOOST_TEST_MESSAGE( o.str() );
          if( rel > worst ) { worst = rel; worst_where = o.str(); }
        }
      }
    }
  }

  // Cylindrical walk: tiny cylinder source at the centre of a cylinder of the same material, with the
  //  detector on and off axis; integrated (the tiny volume makes the outer rule trivial).
  for( const bool dense : { true, false } )
  {
    for( const double lateral : { 0.0, 2.0 } )
    {
      const double energy = 60.0;
      const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
      const shared_ptr<const Material> mat = matdb->material( scenario_matrix_material( dense ) );
      const double mu = transmition_length_coefficient( mat.get(), static_cast<float>(energy) );

      DistributedSrcCalcT<double> calc;
      calc.m_geometry = GeometryType::CylinderEndOn;
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
      calc.m_effResponse = det.mc_transfer;
      calc.m_effMethod = ShieldingSourceFitCalc::VolumetricEffMethod::MCTransfer;
      calc.m_effNumRays = 128;
      const double depth = 0.2;
      calc.m_detector = detector_geom_from_config<double>( GeometryType::CylinderEndOn,
                              (1.0 + depth + sm_point_half_cm)*cm, det.gd.transverse_half_extent()*cm,
                              0.0, lateral*cm, 0.0 );
      DistributedSrcCalcT<double>::ShellInfo src;
      src.dims = { sm_point_half_cm*cm, sm_point_half_cm*cm, 0.0 };
      src.trans_len_coef = mu;
      src.type = ShellType::Material;
      calc.m_shells.push_back( src );
      DistributedSrcCalcT<double>::ShellInfo outer;
      outer.dims = { (4.0 + sm_point_half_cm)*cm, (depth + sm_point_half_cm)*cm, 0.0 };
      outer.trans_len_coef = mu;
      outer.type = ShellType::Material;
      calc.m_shells.push_back( outer );

      DistributedSrcCalcT<double> ref_calc = calc;
      ref_calc.finalize_shell_coefficients();
      const double e[3] = { 0.0, 0.0, 0.0 };
      const double ref = marched_point_kernel( det, ref_calc, e, step );

      self_shielding_integration_imp<double>( calc );
      const double vol = M_PI * src.dims[0]*src.dims[0] * 2.0*src.dims[1];
      const double prod = calc.integral / vol;
      const double rel = fabs( prod/ref - 1.0 );
      ostringstream o;
      o << "  cyl " << (dense ? "steel" : "water") << " depth=" << depth << " lateral=" << lateral
        << ": production " << setprecision(8) << prod << " march " << ref << "  rel " << setprecision(3) << rel;
      BOOST_TEST_MESSAGE( o.str() );
      if( rel > worst ) { worst = rel; worst_where = o.str(); }
    }
  }

  BOOST_TEST_MESSAGE( "worst production-vs-march relative difference " << worst << " at " << worst_where );
  // The march's own error is O(step*mu) ~ 2e-4*9.4 ~ 2e-3 per e-folding at worst; the cylinder rows also
  //  average over the (tiny) source volume.  1e-2 catches any real walk/composition defect.
  BOOST_CHECK_MESSAGE( worst < 1.0e-2, "per-ray kernel differs from the ray-march reference: " << worst_where );
}//BOOST_AUTO_TEST_CASE( PerRayKernelIdentityVsRayMarch )


/** RUNG 1 - the transfer at bare points on the emitting skin (MC, cached).

 Anchored on MC point curves at the source-centre distances (2 and 3 cm) and at 25 cm, all at 0.25%;
 compared against MC at points on the front-face and side-face skin of the large-near geometry and
 at the box/slab corners.  The residual here is the floor every attenuating rung inherits.
 */
BOOST_AUTO_TEST_CASE( Rung1_TransferAtSkinPoints, * boost::unit_test::disabled() )
{
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  McCache cache;

  const vector<double> energies = { 60.0, 344.0, 1332.5 };
  map<double,shared_ptr<ceelo::DetectorResponse>> anchors;
  for( const double d : { 2.0, 3.0, 25.0 } )
    anchors[d] = ladder_anchor_response( det, cache, d, scenario_energies() );

  // (r, z_e) in cm, endcap-front referenced.
  struct Pt { string label; double r, ze; };
  const vector<Pt> pts = {
    { "front r=0",    0.0, 1.05 }, { "front r=1", 1.0, 1.05 }, { "front r=2", 2.0, 1.05 },
    { "front r=3",    3.0, 1.05 }, { "front r=4", 4.0, 1.05 },
    { "side z=1.5",   3.95, 1.5 }, { "side z=3", 3.95, 3.0 }, { "side z=5", 3.95, 5.0 },
    { "box corner",   5.657, 1.05 }, { "slab corner", 4.123, 1.05 },
    { "far r=0",      0.0, 51.0 }, { "far r=4", 4.0, 51.0 },
  };

  vector<unique_ptr<ceelo::Material>> owned;
  ceelo::EfficiencyCalculator pt;
  ceelo::ResponseGenerator::configure_calculator( pt, det.gd, owned );

  BOOST_TEST_MESSAGE( "point            E(keV)   theta   MC eff      sigma  | anchor 2cm | anchor 3cm | anchor 25cm   (transfer/MC - 1)" );
  for( const Pt &p : pts )
  {
    const double theta = atan2( p.r, p.ze );
    const double dist = sqrt( p.r*p.r + p.ze*p.ze );
    const Eigen::Vector3d pos = det.mc_transfer->query_position( theta, 0.0, dist );
    pt.set_point_source( pos );
    ostringstream k;
    k << "point(r=" << setprecision(10) << p.r << ",ze=" << p.ze << ");bare";
    for( const double e : energies )
    {
      const McResult mc = run_mc( pt, cache, k.str(), e );
      ostringstream o;
      o << left << setw(16) << p.label << right << fixed << setprecision(1) << setw(7) << e
        << setprecision(1) << setw(7) << theta*180.0/M_PI << "  " << scientific << setprecision(4) << mc.eff
        << fixed << "  " << setprecision(2) << 100.0*mc.frac_sigma << "%" << (mc.hit_cap() ? "CAP" : "   ");
      for( const double d : { 2.0, 3.0, 25.0 } )
      {
        const double xfer = anchors[d]->eps_fep( e, theta, 0.0, dist ).value;
        o << " | " << pct( xfer/mc.eff - 1.0 );
      }
      BOOST_TEST_MESSAGE( o.str() );
    }
  }
}//BOOST_AUTO_TEST_CASE( Rung1_TransferAtSkinPoints )


/** RUNG 2 - transparent whole volumes (MC, cached): the model with every coefficient zeroed against
 an emitting volume in vacuum, 0.25% both sides. */
BOOST_AUTO_TEST_CASE( Rung2_TransparentVolumes, * boost::unit_test::disabled() )
{
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  McCache cache;
  const shared_ptr<ceelo::DetectorResponse> far_anchor = ladder_anchor_response( det, cache, 25.0, anchor_energies() );

  BOOST_TEST_MESSAGE( "scenario (transparent)    E(keV)   MC eff      sigma  | model(25cm anchor)/MC-1 | model(centre anchor)/MC-1" );
  for( const char *name : { "small-near-light", "large-near-light", "box-large-near-light",
                            "box-slab-near-light", "shielded-near-light", "box-shielded-near-light" } )
  {
    const Scenario s = find_scenario( name );
    const shared_ptr<ceelo::DetectorResponse> near_anchor =
          ladder_anchor_response( det, cache, scenario_center_distance_cm( s ), scenario_energies() );
    ScenarioMcMaterials mats;
    ceelo::EfficiencyCalculator calc;
    configure_scenario_mc( calc, det, s, mats, true );
    const string key = scenario_mc_key( s, true );
    for( const double e : { 60.0, 344.0, 1332.5 } )
    {
      const McResult mc = run_mc( calc, cache, key, e );
      const double m_far = interspec_volumetric_eff( det, s, e, far_anchor, -1, -1.0, -1.0, true );
      const double m_near = interspec_volumetric_eff( det, s, e, near_anchor, -1, -1.0, -1.0, true );
      ostringstream o;
      o << left << setw(26) << s.name << right << fixed << setprecision(1) << setw(7) << e
        << "  " << scientific << setprecision(4) << mc.eff << fixed << "  " << setprecision(2)
        << 100.0*mc.frac_sigma << "%" << (mc.hit_cap() ? "CAP" : "   ")
        << " | " << pct( m_far/mc.eff - 1.0 ) << " | " << pct( m_near/mc.eff - 1.0 );
      BOOST_TEST_MESSAGE( o.str() );
    }
  }
}//BOOST_AUTO_TEST_CASE( Rung2_TransparentVolumes )


/** RUNG 3 - one attenuated point at depth in a slab (MC, cached).  See #PointAtDepth.

 Both the legacy and the fixed aperture frame are evaluated on the InterSpec side, against the same
 MC, with the transfer anchored at the point's own distance (removing the transfer's near-field point
 bias) and at 25 cm.
 */
BOOST_AUTO_TEST_CASE( Rung3_PointAtDepth, * boost::unit_test::disabled() )
{
  using namespace GammaInteractionCalc;
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  McCache cache;
  const bool saved = sm_aperture_frame_legacy_origin;

  const shared_ptr<ceelo::DetectorResponse> far_anchor = ladder_anchor_response( det, cache, 25.0, anchor_energies() );
  map<long long,shared_ptr<ceelo::DetectorResponse>> near_anchors;
  const auto anchor_for = [&]( const double d ) {
    const long long k = llround( d*1.0e6 );
    if( !near_anchors.count(k) )
      near_anchors[k] = ladder_anchor_response( det, cache, d, scenario_energies() );
    return near_anchors[k];
  };

  // The FEP window the truth was generated with: a photon scattering by less than this stays in the
  //  peak, so the model must not attenuate it away.  CeeLo's convention is a +-half-window.
  const double fep_win = ceelo::kDefaultFepWindowKeV;

  BOOST_TEST_MESSAGE( "model/MC - 1, near = transfer anchored at the point's own distance, 25cm = far anchor;"
                      " +win = survival-removal mu at a " + to_string(fep_win) + " keV FEP half-window" );
  BOOST_TEST_MESSAGE( "material depth lateral standoff  E(keV)   tau    MC eff      sigma   | legacy(near)  fixed(near)  fixed+win(near) | fixed(25cm)  fixed+win(25cm)" );

  for( const bool dense : { true, false } )
  {
    const vector<double> depths = dense ? vector<double>{ 0.02, 0.05, 0.1, 0.2, 0.5 }
                                        : vector<double>{ 0.5, 1.0, 2.0, 5.0 };
    const vector<double> energies = dense ? vector<double>{ 60.0, 344.0, 1332.5 } : vector<double>{ 60.0 };
    for( const double standoff : { 1.0, 5.0, 50.0 } )
    {
      for( const double lateral : { 0.0, 2.0, 4.0 } )
      {
        for( const double depth : depths )
        {
          PointAtDepth p;
          p.depth_cm = depth;  p.lateral_cm = lateral;  p.standoff_cm = standoff;  p.dense = dense;

          for( const double e : energies )
          {
            DistributedSrcCalcT<double> probe = build_point_at_depth_calc( det, p, e, far_anchor );
            const double tau = probe.m_shells[1].trans_len_coef * depth * PhysicalUnits::cm;
            // Far rows only exist as a null check; deep ones are infeasible at 0.25%.
            if( (standoff > 20.0) && (tau > 1.6) )
              continue;
            const double target = (standoff > 20.0) ? 0.01 : 0.0025;

            ScenarioMcMaterials mats;
            ceelo::EfficiencyCalculator mc_calc;
            configure_point_at_depth_mc( mc_calc, det, p, mats );
            const McResult mc = run_mc( mc_calc, cache, point_at_depth_key( p ), e, target );

            const shared_ptr<ceelo::DetectorResponse> near_anchor = anchor_for( point_center_distance_cm( p ) );
            // legacy(near), fixed(near), fixed+win(near), fixed(25cm), fixed+win(25cm)
            double v[5];
            {
              sm_aperture_frame_legacy_origin = true;
              DistributedSrcCalcT<double> c0 = build_point_at_depth_calc( det, p, e, near_anchor );
              v[0] = point_kernel_eff( c0 );
              sm_aperture_frame_legacy_origin = false;
              DistributedSrcCalcT<double> c1 = build_point_at_depth_calc( det, p, e, near_anchor );
              v[1] = point_kernel_eff( c1 );
              DistributedSrcCalcT<double> c2 = build_point_at_depth_calc( det, p, e, near_anchor, -1, false, fep_win );
              v[2] = point_kernel_eff( c2 );
              DistributedSrcCalcT<double> c3 = build_point_at_depth_calc( det, p, e, far_anchor );
              v[3] = point_kernel_eff( c3 );
              DistributedSrcCalcT<double> c4 = build_point_at_depth_calc( det, p, e, far_anchor, -1, false, fep_win );
              v[4] = point_kernel_eff( c4 );
            }
            sm_aperture_frame_legacy_origin = saved;

            ostringstream o;
            o << left << setw(8) << (dense ? "steel" : "water") << right << fixed
              << setprecision(2) << setw(6) << depth << setw(8) << lateral << setw(9) << standoff
              << setprecision(1) << setw(8) << e << setprecision(2) << setw(7) << tau
              << "  " << scientific << setprecision(4) << mc.eff << fixed << "  " << setprecision(2)
              << 100.0*mc.frac_sigma << "%" << (mc.hit_cap() ? "CAP" : "   ")
              << " | " << pct( v[0]/mc.eff - 1.0 ) << " " << pct( v[1]/mc.eff - 1.0 )
              << "  " << pct( v[2]/mc.eff - 1.0 )
              << " | " << pct( v[3]/mc.eff - 1.0 ) << "  " << pct( v[4]/mc.eff - 1.0 );
            BOOST_TEST_MESSAGE( o.str() );
          }//for( energies )
        }//for( depths )
      }//for( lateral )
    }//for( standoff )
  }//for( dense )

  sm_aperture_frame_legacy_origin = saved;
}//BOOST_AUTO_TEST_CASE( Rung3_PointAtDepth )


/** RUNG 4a - slabs of increasing optical depth (MC, cached): cylinders r=4 cm of thickness t at
 standoff 1, 15 and 50 cm, steel and water, legacy vs fixed frame, near-centre and 25 cm anchors. */
BOOST_AUTO_TEST_CASE( Rung4_SlabThicknessSweep, * boost::unit_test::disabled() )
{
  using namespace GammaInteractionCalc;
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  McCache cache;
  const bool saved = sm_aperture_frame_legacy_origin;
  const shared_ptr<ceelo::DetectorResponse> far_anchor = ladder_anchor_response( det, cache, 25.0, anchor_energies() );
  map<long long,shared_ptr<ceelo::DetectorResponse>> near_anchors;
  const auto anchor_for = [&]( const double d ) {
    const long long k = llround( d*1.0e6 );
    if( !near_anchors.count(k) )
      near_anchors[k] = ladder_anchor_response( det, cache, d, scenario_energies() );
    return near_anchors[k];
  };

  const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();

  BOOST_TEST_MESSAGE( "material  t(cm) standoff  E(keV)   tau    MC eff      sigma   | legacy(near)  fixed(near) | legacy(25cm)  fixed(25cm)   (model/MC-1)" );
  for( const bool dense : { true, false } )
  {
    const shared_ptr<const Material> mat = matdb->material( scenario_matrix_material( dense ) );
    for( const double standoff : { 1.0, 15.0, 50.0 } )
    {
      for( const double e : { 60.0, 344.0, 1332.5 } )
      {
        // transmition_length_coefficient is per PhysicalUnits length, so MULTIPLYING by `cm` gives
        //  the per-centimetre coefficient (the same way scenario_optical_depth forms its tau).
        const double mu = transmition_length_coefficient( mat.get(), static_cast<float>(e) ) * PhysicalUnits::cm;
        for( const double tau : { 0.1, 0.3, 1.0, 3.0, 10.0 } )
        {
          const double t = tau / mu;   //cm
          if( t > 20.0 )
            continue;   //water at high energy: absurdly long cylinders

          Scenario s;
          ostringstream nm;
          nm << "slab-" << (dense ? "steel" : "water") << "-t" << fixed << setprecision(4) << t << "-d" << standoff;
          s.name = nm.str();
          s.radius_cm = 4.0;
          s.half_length_cm = 0.5*t;
          s.standoff_cm = standoff;
          s.dense = dense;
          s.shield_cm = 0.0;

          if( (standoff > 20.0) && (tau > 3.5) )
            continue;
          const double target = (standoff > 20.0) ? 0.01 : 0.0025;

          ScenarioMcMaterials mats;
          ceelo::EfficiencyCalculator mc_calc;
          configure_scenario_mc( mc_calc, det, s, mats );
          const McResult mc = run_mc( mc_calc, cache, scenario_mc_key( s ), e, target );

          const shared_ptr<ceelo::DetectorResponse> near_anchor = anchor_for( scenario_center_distance_cm( s ) );
          double v[4];
          int idx = 0;
          for( const shared_ptr<ceelo::DetectorResponse> &resp : { near_anchor, far_anchor } )
          {
            for( const bool legacy : { true, false } )
            {
              sm_aperture_frame_legacy_origin = legacy;
              v[idx++] = interspec_volumetric_eff( det, s, e, resp );
            }
          }
          sm_aperture_frame_legacy_origin = saved;

          ostringstream o;
          o << left << setw(8) << (dense ? "steel" : "water") << right << fixed
            << setprecision(3) << setw(7) << t << setprecision(1) << setw(8) << standoff << setw(8) << e
            << setprecision(2) << setw(7) << tau
            << "  " << scientific << setprecision(4) << mc.eff << fixed << "  " << setprecision(2)
            << 100.0*mc.frac_sigma << "%" << (mc.hit_cap() ? "CAP" : "   ")
            << " | " << pct( v[0]/mc.eff - 1.0 ) << " " << pct( v[1]/mc.eff - 1.0 )
            << " | " << pct( v[2]/mc.eff - 1.0 ) << " " << pct( v[3]/mc.eff - 1.0 );
          BOOST_TEST_MESSAGE( o.str() );
        }//for( tau )
      }//for( energies )
    }//for( standoff )
  }//for( dense )
  sm_aperture_frame_legacy_origin = saved;
}//BOOST_AUTO_TEST_CASE( Rung4_SlabThicknessSweep )


/** RUNG 4b - the full scenario matrix at 0.25% (MC, cached): the new truth bank.  Prints the
 model/MC comparison (legacy and fixed, 25 cm anchor - the bank's own convention) and a paste-ready
 TruthRow table plus the 25 cm anchor curve for VolumetricNearFieldTruth.h. */
BOOST_AUTO_TEST_CASE( Rung4_ScenarioMatrixTruth, * boost::unit_test::disabled() )
{
  using namespace GammaInteractionCalc;
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  McCache cache;
  const bool saved = sm_aperture_frame_legacy_origin;

  // Anchor curve first (it is what the model column is grounded on).
  vector<unique_ptr<ceelo::Material>> owned;
  ceelo::EfficiencyCalculator pt;
  ceelo::ResponseGenerator::configure_calculator( pt, det.gd, owned );
  pt.set_point_source( Eigen::Vector3d( 0.0, 0.0, -(det.endcap_front_offset_cm + kMcAnchorDistanceCm) ) );
  ostringstream ak;
  ak << "point(r=0,ze=" << setprecision(10) << kMcAnchorDistanceCm << ");bare";
  cout << "// Regenerated by Rung4_ScenarioMatrixTruth at 0.25%.  FEP window: "
       << ceelo::kDefaultFepWindowKeV << " keV.\nstatic const std::vector<AnchorRow> sm_mc_anchor = {\n";
  for( const double e : anchor_energies() )
  {
    const McResult r = run_mc( pt, cache, ak.str(), e );
    cout << "  { " << setprecision(12) << e << ", " << setprecision(10) << r.eff << ", "
         << setprecision(4) << r.frac_sigma << " },\n";
  }
  cout << "};\n" << flush;
  const shared_ptr<ceelo::DetectorResponse> far_anchor = ladder_anchor_response( det, cache, 25.0, anchor_energies() );

  cout << "static const std::vector<TruthRow> sm_truth = {\n" << flush;
  BOOST_TEST_MESSAGE( "scenario                 E(keV)    tau    MC eff      sigma   | legacy/MC-1  fixed/MC-1" );
  for( const Scenario &s : scenarios() )
  {
    ScenarioMcMaterials mats;
    ceelo::EfficiencyCalculator mc_calc;
    configure_scenario_mc( mc_calc, det, s, mats );
    for( const double e : scenario_energies() )
    {
      const McResult mc = run_mc( mc_calc, cache, scenario_mc_key( s ), e );
      cout << "  { \"" << s.name << "\", " << setprecision(12) << e << ", "
           << setprecision(10) << mc.eff << ", " << setprecision(4) << mc.eff*mc.frac_sigma << " },"
           << (mc.hit_cap() ? "  // stopped on the event cap" : "") << "\n" << flush;

      sm_aperture_frame_legacy_origin = true;
      const double legacy = interspec_volumetric_eff( det, s, e, far_anchor );
      sm_aperture_frame_legacy_origin = false;
      const double fixed_v = interspec_volumetric_eff( det, s, e, far_anchor );
      sm_aperture_frame_legacy_origin = saved;

      ostringstream o;
      o << left << setw(25) << s.name << right << fixed << setprecision(1) << setw(7) << e
        << setprecision(2) << setw(7) << scenario_optical_depth( s, e )
        << "  " << scientific << setprecision(4) << mc.eff << fixed << "  " << setprecision(2)
        << 100.0*mc.frac_sigma << "%" << (mc.hit_cap() ? "CAP" : "   ")
        << " | " << pct( legacy/mc.eff - 1.0 ) << " " << pct( fixed_v/mc.eff - 1.0 );
      BOOST_TEST_MESSAGE( o.str() );
    }
  }
  cout << "};\n" << flush;
  sm_aperture_frame_legacy_origin = saved;
}//BOOST_AUTO_TEST_CASE( Rung4_ScenarioMatrixTruth )


namespace
{
/** The mu the FEP leg should use, as a function of the window half-width and of how deep the photon
 already is, so the credit can be studied as a one-line policy rather than a hard-wired choice.

   `win_keV <= 0`     : mu_total.  No credit - the choice every model/truth number before 2026-09-02
                        was silently making.
   `win_keV > 0`      : mu_total - f_win*mu_Compton (`fep_removal_coefficient`), the flat credit.
   `tau_c > 0`        : the same, with the Compton term scaled by exp(-tau_src/tau_c) - CeeLo's
                        depth-aware form (`fep_depth_survival_credit`, kFepDepthTauC = 5 mfp).  A
                        photon that has already crossed many mean free paths is unlikely to have
                        stayed in the peak while doing so, so the in-window credit must fade.

 `tau_src` is the mu_total optical depth the photon has to cross to leave, i.e. what the flat credit
 implicitly assumes is zero.  The shell walk carries ONE coefficient per shell, so a genuinely
 per-ray depth-aware mu is a kernel change; this evaluates the policy at the geometry's own tau
 first, which is what decides whether that change is worth making.
 */
double fep_leg_mu( const Material &mat, const double energy_keV, const double win_keV,
                   const double tau_src, const double tau_c )
{
  const double mu_total = GammaInteractionCalc::transmition_length_coefficient(
                                                  &mat, static_cast<float>(energy_keV) );
  if( !(win_keV > 0.0) )
    return mu_total;

  const double flat = fep_removal_coefficient( mat, energy_keV, win_keV );
  if( !(tau_c > 0.0) )
    return flat;

  // flat == mu_total - f_win*mu_cs, so (mu_total - flat) IS the credit; fade it with depth.
  return mu_total - (mu_total - flat)*std::exp( -std::max(0.0,tau_src)/tau_c );
}//fep_leg_mu(...)

}//namespace


/** RUNG 5 - the full-energy-peak window credit, studied against the CACHED Monte Carlo.

 Every model-vs-truth number in this project before 2026-09-02 used `mu_rem = mu_total`, which is not
 the neutral choice it looks like: it attenuates away the forward Compton scatters that lose less than
 the peak window and therefore stay in the peak.  Rung 3 showed the credit closing water at 60 keV from
 -6.4% to +0.5%, and over-crediting deep steel.  This case measures the policy across the whole cached
 grid so the choice is made on data:

   mu_total, then the flat credit at 0.375 / 0.75 / 1.5 keV half-windows (the truth bank was generated
   at CeeLo's 0.75 keV default, so that column is the like-for-like one and the other two are the
   sensitivity), then the depth-aware credit at 0.75 keV with tau_c = 5 mfp.

 Runs entirely off the Monte-Carlo cache: no new MC, so it is cheap to re-run after any model change.
 */
BOOST_AUTO_TEST_CASE( Rung5_FepWindowCredit, * boost::unit_test::disabled() )
{
  using namespace GammaInteractionCalc;
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  McCache cache;
  const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();

  const shared_ptr<ceelo::DetectorResponse> far_anchor
        = ladder_anchor_response( det, cache, 25.0, anchor_energies() );
  map<long long,shared_ptr<ceelo::DetectorResponse>> near_anchors;
  const auto anchor_for = [&]( const double d ) {
    const long long k = llround( d*1.0e6 );
    if( !near_anchors.count(k) )
      near_anchors[k] = ladder_anchor_response( det, cache, d, scenario_energies() );
    return near_anchors[k];
  };

  BOOST_TEST_MESSAGE( "model/MC - 1, transfer anchored at the source's own distance.  win = FEP"
                      " half-window (keV); depth = the same 0.75 keV credit faded by"
                      " exp(-tau/5mfp)." );
  BOOST_TEST_MESSAGE( "material depth lateral standoff  E(keV)   tau   | mu_total   win0.375   win0.75"
                      "   win1.50    depth0.75" );

  // (a) The point-at-depth grid: one eval_rect per row, no outer quadrature at all.
  for( const bool dense : { true, false } )
  {
    const shared_ptr<const Material> mat = matdb->material( scenario_matrix_material( dense ) );
    BOOST_REQUIRE( mat );
    const vector<double> depths = dense ? vector<double>{ 0.02, 0.05, 0.1, 0.2, 0.5 }
                                        : vector<double>{ 0.5, 1.0, 2.0, 5.0 };
    const vector<double> energies = dense ? vector<double>{ 60.0, 344.0, 1332.5 }
                                          : vector<double>{ 60.0 };
    for( const double standoff : { 1.0, 5.0, 50.0 } )
    {
      for( const double lateral : { 0.0, 2.0, 4.0 } )
      {
        for( const double depth : depths )
        {
          PointAtDepth p;
          p.depth_cm = depth;  p.lateral_cm = lateral;  p.standoff_cm = standoff;  p.dense = dense;

          for( const double e : energies )
          {
            const double mu_tot = transmition_length_coefficient( mat.get(), static_cast<float>(e) );
            const double tau = mu_tot * depth * PhysicalUnits::cm;
            if( (standoff > 20.0) && (tau > 1.6) )
              continue;   //the row rung 3 skipped; its MC is not in the cache
            const double target = (standoff > 20.0) ? 0.01 : 0.0025;

            ScenarioMcMaterials mats;
            ceelo::EfficiencyCalculator mc_calc;
            configure_point_at_depth_mc( mc_calc, det, p, mats );
            const McResult mc = run_mc( mc_calc, cache, point_at_depth_key( p ), e, target );
            BOOST_CHECK_MESSAGE( mc.from_cache, "Rung5 ran Monte Carlo (cache miss) - it is meant to"
                                 " re-score rung 3/4a results, so a miss means the key drifted" );

            const shared_ptr<ceelo::DetectorResponse> resp = anchor_for( point_center_distance_cm( p ) );

            ostringstream o;
            o << left << setw(8) << (dense ? "steel" : "water") << right << fixed
              << setprecision(2) << setw(6) << depth << setw(8) << lateral << setw(9) << standoff
              << setprecision(1) << setw(8) << e << setprecision(2) << setw(7) << tau << "  |";

            // mu_total, three flat windows, then the depth-aware one.
            const double wins[4] = { -1.0, 0.375, 0.75, 1.5 };
            for( int i = 0; i < 5; ++i )
            {
              const double win = (i < 4) ? wins[i] : 0.75;
              const double tau_c = (i < 4) ? -1.0 : ceelo::kFepDepthTauC;
              DistributedSrcCalcT<double> c = build_point_at_depth_calc( det, p, e, resp );
              const double mu = fep_leg_mu( *mat, e, win, tau, tau_c );
              for( DistributedSrcCalcT<double>::ShellInfo &sh : c.m_shells )
                sh.trans_len_coef = mu;
              o << " " << pct( point_kernel_eff( c )/mc.eff - 1.0 );
            }
            BOOST_TEST_MESSAGE( o.str() );
          }//for( energies )
        }//for( depths )
      }//for( lateral )
    }//for( standoff )
  }//for( dense )

  // (b) The slab volume grid, so the same policy is read on a whole integral rather than one point.
  BOOST_TEST_MESSAGE( "" );
  BOOST_TEST_MESSAGE( "volume slabs (r=4 cm cylinders):" );
  BOOST_TEST_MESSAGE( "material  t(cm) standoff  E(keV)   tau   | mu_total   win0.375   win0.75"
                      "   win1.50    depth0.75" );
  for( const bool dense : { true, false } )
  {
    const shared_ptr<const Material> mat = matdb->material( scenario_matrix_material( dense ) );
    for( const double standoff : { 1.0, 15.0, 50.0 } )
    {
      for( const double e : { 60.0, 344.0, 1332.5 } )
      {
        const double mu_per_cm = transmition_length_coefficient( mat.get(), static_cast<float>(e) )
                                 * PhysicalUnits::cm;
        for( const double tau : { 0.1, 0.3, 1.0, 3.0, 10.0 } )
        {
          const double t = tau / mu_per_cm;
          if( t > 20.0 )
            continue;
          if( (standoff > 20.0) && (tau > 3.5) )
            continue;
          const double target = (standoff > 20.0) ? 0.01 : 0.0025;

          Scenario s;
          ostringstream nm;
          nm << "slab-" << (dense ? "steel" : "water") << "-t" << fixed << setprecision(4) << t
             << "-d" << standoff;
          s.name = nm.str();
          s.radius_cm = 4.0;
          s.half_length_cm = 0.5*t;
          s.standoff_cm = standoff;
          s.dense = dense;
          s.shield_cm = 0.0;

          ScenarioMcMaterials mats;
          ceelo::EfficiencyCalculator mc_calc;
          configure_scenario_mc( mc_calc, det, s, mats );
          const McResult mc = run_mc( mc_calc, cache, scenario_mc_key( s ), e, target );
          BOOST_CHECK_MESSAGE( mc.from_cache, "Rung5 ran Monte Carlo (cache miss)" );

          const shared_ptr<ceelo::DetectorResponse> resp
                = anchor_for( scenario_center_distance_cm( s ) );

          ostringstream o;
          o << left << setw(8) << (dense ? "steel" : "water") << right << fixed
            << setprecision(3) << setw(7) << t << setprecision(1) << setw(8) << standoff
            << setw(8) << e << setprecision(2) << setw(7) << tau << "  |";

          const double wins[4] = { -1.0, 0.375, 0.75, 1.5 };
          for( int i = 0; i < 5; ++i )
          {
            const double win = (i < 4) ? wins[i] : 0.75;
            const double tau_c = (i < 4) ? -1.0 : ceelo::kFepDepthTauC;
            // interspec_volumetric_eff builds its own shells, so drive it through fep_window_keV for
            //  the flat columns and fall back to a hand-built calculator for the depth-aware one.
            double v = 0.0;
            if( tau_c <= 0.0 )
            {
              v = interspec_volumetric_eff( det, s, e, resp, -1, win );
            }else
            {
              DistributedSrcCalcT<double> c = build_scenario_calc( det, s, e, resp );
              const double mu = fep_leg_mu( *mat, e, win, tau, tau_c );
              for( DistributedSrcCalcT<double>::ShellInfo &sh : c.m_shells )
                sh.trans_len_coef = mu;
              self_shielding_integration_imp<double>( c );
              v = c.integral / (scenario_volume_cm3(s) * pow(PhysicalUnits::cm,3.0));
            }
            o << " " << pct( v/mc.eff - 1.0 );
          }
          BOOST_TEST_MESSAGE( o.str() );
        }//for( tau )
      }//for( energies )
    }//for( standoff )
  }//for( dense )
}//BOOST_AUTO_TEST_CASE( Rung5_FepWindowCredit )


/** COST: what one energy of the volumetric model actually costs, split into its pieces.

 The per-element aperture is the unit of work: each element asks the response for a fan of
 `m_effNumRays` ray directions (a CeeLo ray trace apiece) and then ray-traces every one of them
 through InterSpec's shells.  So the cost of a volume integral is
 (number of elements the adaptive rule uses) x (rays) x (one CeeLo trace + one InterSpec walk), and
 the useful numbers are the per-element cost, the element count, and the total.

 Reported as CPU seconds (the model is single-threaded, so CPU == wall here).
 */
BOOST_AUTO_TEST_CASE( CostPerEnergy )
{
  using namespace GammaInteractionCalc;
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );
  const double cm = PhysicalUnits::cm;

  const auto cpu_now = []() -> double {
    return static_cast<double>( std::clock() ) / CLOCKS_PER_SEC;
  };

  // (1) One bare transfer query: eps_fep at a point.  This is what the FLAT-disk-plus-response model
  //     costs per element, and what a point-source efficiency costs outright.
  {
    const int n = 2000;
    const double t0 = cpu_now();
    double sink = 0.0;
    for( int i = 0; i < n; ++i )
      sink += det.mc_transfer->eps_fep( 60.0 + 0.001*i, 0.3, 0.0, 3.0 ).value;
    const double dt = (cpu_now() - t0) / n;
    ostringstream o;
    o << "  one eps_fep transfer query: " << scientific << setprecision(3) << dt << " s"
      << " (2048-ray internal quadrature)   [sink " << sink << "]";
    BOOST_TEST_MESSAGE( o.str() );
  }

  // (2) One per-element aperture at the shipped 128 rays (was 512 until 2026-09-03), 512 and 2048, so the scaling is visible.
  for( const int rays : { 128, 512, 2048 } )
  {
    const int n = 200;
    const double to_det[3] = { 0.6, 0.0, 0.8 };
    const double axis[3] = { 0.0, 0.0, -1.0 };
    const double t0 = cpu_now();
    size_t sink = 0;
    for( int i = 0; i < n; ++i )
      sink += build_element_aperture( *det.mc_transfer, 60.0 + 0.001*i, 3.0*cm, 0.8, to_det, axis,
                                      rays ).weights.size();
    const double dt = (cpu_now() - t0) / n;
    ostringstream o;
    o << "  one element aperture, " << setw(4) << rays << " rays: " << scientific << setprecision(3)
      << dt << " s   [" << (sink/n) << " rays kept]";
    BOOST_TEST_MESSAGE( o.str() );
  }

  // (3) A whole volume integral per energy, for representative scenarios: total, element count, and
  //     the cost per element (which should be the aperture cost plus the shell walk).
  BOOST_TEST_MESSAGE( "" );
  BOOST_TEST_MESSAGE( "  scenario                   E(keV)   CPU s   elements   s/element  ndim" );
  for( const char *name : { "small-near-light", "large-near-dense", "shielded-near-dense",
                            "box-large-near-dense", "box-slab-near-dense", "large-far-dense" } )
  {
    const Scenario s = find_scenario( name );
    for( const double e : { 60.0, 1332.5 } )
    {
      size_t evals = 0;
      double est_err = 0.0;
      const double t0 = cpu_now();
      const double eff = interspec_volumetric_eff( det, s, e, det.mc_transfer, -1, -1.0, -1.0, false,
                                                   &evals, &est_err );
      const double dt = cpu_now() - t0;
      ostringstream o;
      o << "  " << left << setw(26) << s.name << right << fixed << setprecision(1) << setw(7) << e
        << setprecision(2) << setw(9) << dt << setw(11) << evals << "  " << scientific
        << setprecision(2) << (evals ? dt/evals : 0.0) << "   "
        << ((s.shape == Shape::Box) ? 3 : 2) << "   [eff " << setprecision(3) << eff << "]";
      BOOST_TEST_MESSAGE( o.str() );
    }
  }

  // (3b) The LINE path: what one cost-function evaluation of a whole fit costs.  This is the number
  //      that matters, and it is not (3) divided by anything: the detector-side line set is built
  //      ONCE per fit, every energy of a source shares it, and the fit runs the whole thing in
  //      ceres::Jet<double,16> for the Jacobian.  So: the one-off build, then one evaluation of
  //      `num_energies` energies together, in double and in Jet.
  {
    using Jet16 = ceres::Jet<double,16>;
    BOOST_TEST_MESSAGE( "" );
    BOOST_TEST_MESSAGE( "  LINE path, per cost-function evaluation (one shared line set, all energies):" );
    const std::vector<double> energies = { 60.0, 88.0, 122.0, 344.0, 411.0, 661.7, 778.0, 964.0,
                                           1112.0, 1332.5 };
    for( const char *name : { "large-near-dense", "box-large-near-dense" } )
    {
      const Scenario s = find_scenario( name );
      for( const int num_lines : { 1 << 16, 1 << 17 } )
      {
        DistributedSrcCalcT<double> proto = build_scenario_calc( det, s, energies.front(), det.mc_transfer );
        const double tb0 = cpu_now();
        attach_line_cache( proto, num_lines );
        const double build_s = cpu_now() - tb0;

        // One evaluation = every energy of this source, integrated together.
        std::vector<std::unique_ptr<DistributedSrcCalcT<double>>> dcalcs;
        for( const double e : energies )
        {
          auto c = std::make_unique<DistributedSrcCalcT<double>>( proto );
          c->m_energy = e;
          dcalcs.push_back( std::move(c) );
        }
        const ScopedVolumetricIntegratorOverride force( VolumetricIntegrator::Line );
        const double t0 = cpu_now();
        integrate_volumetric_calculators<double>( dcalcs, true );
        const double dbl_s = cpu_now() - t0;

        // The same, with the source dimensions carrying autodiff seeds (as a dimension fit does).
        std::vector<std::unique_ptr<DistributedSrcCalcT<Jet16>>> jcalcs;
        for( const double e : energies )
        {
          auto c = std::make_unique<DistributedSrcCalcT<Jet16>>();
          c->m_geometry = proto.m_geometry;
          c->m_materialIndex = proto.m_materialIndex;
          c->m_attenuateForAir = proto.m_attenuateForAir;
          c->m_airTransLenCoef = proto.m_airTransLenCoef;
          c->m_isInSituExponential = false;
          c->m_inSituRelaxationLength = -1.0;
          c->m_srcVolumetricActivity = Jet16( 1.0 );
          c->m_normalizeByVolume = false;
          c->m_energy = e;
          c->m_effResponse = proto.m_effResponse;
          c->m_effMethod = proto.m_effMethod;
          c->m_lineCache = proto.m_lineCache;
          for( int i = 0; i < 3; ++i )
          {
            c->m_detector.position[i] = Jet16( proto.m_detector.position[i] );
            c->m_detector.axis[i] = proto.m_detector.axis[i];
          }
          c->m_detector.radius = proto.m_detector.radius;
          c->m_detector.setback = proto.m_detector.setback;
          for( size_t sh = 0; sh < proto.m_shells.size(); ++sh )
          {
            DistributedSrcCalcT<Jet16>::ShellInfo info;
            for( int i = 0; i < 3; ++i )
            {
              info.dims[i] = Jet16( proto.m_shells[sh].dims[i] );
              if( sh == proto.m_materialIndex )
                info.dims[i].v[i] = 1.0;      //seed d/d(source dims)
            }
            info.trans_len_coef = Jet16( proto.m_shells[sh].trans_len_coef );
            info.type = proto.m_shells[sh].type;
            c->m_shells.push_back( info );
          }
          jcalcs.push_back( std::move(c) );
        }
        const double t1 = cpu_now();
        integrate_volumetric_calculators<Jet16>( jcalcs, true );
        const double jet_s = cpu_now() - t1;

        ostringstream o;
        o << "  " << left << setw(22) << s.name << right << "  " << setw(7) << num_lines
          << " lines:  build " << fixed << setprecision(3) << build_s << " s (once)"
          << ",  eval " << setprecision(3) << dbl_s << " s double"
          << ",  " << setprecision(3) << jet_s << " s Jet<16>"
          << "   [" << energies.size() << " energies]";
        BOOST_TEST_MESSAGE( o.str() );
      }
    }
  }

  // (4) The flat-disk model on the same scenarios, as the baseline the per-ray kernel is paid against.
  BOOST_TEST_MESSAGE( "" );
  BOOST_TEST_MESSAGE( "  flat-disk (no response, no per-ray kernel), same scenarios:" );
  for( const char *name : { "large-near-dense", "box-large-near-dense" } )
  {
    const Scenario s = find_scenario( name );
    for( const double e : { 60.0, 1332.5 } )
    {
      size_t evals = 0;
      double est_err = 0.0;
      const double t0 = cpu_now();
      interspec_volumetric_eff( det, s, e, nullptr, -1, -1.0, -1.0, false, &evals, &est_err );
      const double dt = cpu_now() - t0;
      ostringstream o;
      o << "  " << left << setw(26) << s.name << right << fixed << setprecision(1) << setw(7) << e
        << setprecision(4) << setw(11) << dt << setw(11) << evals;
      BOOST_TEST_MESSAGE( o.str() );
    }
  }
}//BOOST_AUTO_TEST_CASE( CostPerEnergy )


/** LINE vs ELEMENT - the detector-side line integration against the per-element aperture path.

 Both quadratures evaluate the SAME integral (VolumetricLineIntegration_imp.hpp derives the
 reversal), on the same transfer response, the same shells and the same attenuation
 coefficients; only the discretisation differs.  The element path's own aperture discretisation
 is ~0.1% (ApertureRayConvergence), so the two must agree to a few tenths of a percent across the
 whole scenario matrix - cylinders, boxes, off-axis and side-on, shielded, transparent - or one of
 them is wrong.  Reports the cost of each per energy as well.
 */
BOOST_AUTO_TEST_CASE( LineVsElementScenarioMatrix )
{
  using namespace GammaInteractionCalc;
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );

  const auto cpu_now = []() -> double {
    return static_cast<double>( std::clock() ) / CLOCKS_PER_SEC;
  };

  const int num_lines = 1 << 16;
  const std::vector<double> energies = { 60.0, 122.0, 661.7, 1332.5 };

  double worst = 0.0;
  string worst_where;
  double elem_cpu = 0.0, line_cpu = 0.0;

  // The element path is the VALIDATION REFERENCE now, and it is what this case spends its time
  //  on (the line side is 0.1 s per row throughout): the whole matrix - 28 scenarios x
  //  {attenuating, transparent} x 4 energies - costs ~25 minutes, essentially all of it element
  //  time, with the boxes and the off-axis cylinders at 5-120 s per energy.  Every run therefore
  //  does a QUICK SUBSET chosen to cover each convention once, ~50 s of element time in all:
  //    small-near-light         cylinder end-on at contact, well inside the crystal radius
  //    large-near-dense         optically thick, wider than the crystal
  //    shielded-near-dense      an outer shield in the walk
  //    sideon-tall-near-light   the side-on axis convention
  //    box-shielded-near-dense  the only 3-D geometry, and the cheapest box (1.8-3.5 s/energy)
  //  INTERSPEC_LINE_AB_SCENARIOS=all runs the full matrix; =name[,name...] runs those scenarios.
  //  Full-matrix element cost per energy, measured 2026-09-03: cheap (0.1-0.9 s) large-far-*,
  //  small-far-*, small-near-*, shielded-near-dense, wide-angle-far-light; medium (2-4 s)
  //  sideon-tall-near-light, box-shielded-near-dense; slow box-large-near-* 5-10, sideon-squat-near
  //  6-20, offaxis-small-near-* 9-29, box-slab-near-* 18-46, box-shielded-near-light 3-62,
  //  offaxis-large-near-light 62-120.  The full matrix measured a worst row of 0.44% (transparent
  //  box at 122 keV); the quick subset 0.29% (large-near-dense at 60 keV).
  std::set<string> only = { "small-near-light", "large-near-dense", "shielded-near-dense",
                            "sideon-tall-near-light", "box-shielded-near-dense" };
  if( const char *env = std::getenv("INTERSPEC_LINE_AB_SCENARIOS") )
  {
    only.clear();
    string spec = env;
    size_t pos = 0;
    while( pos <= spec.size() )
    {
      const size_t comma = spec.find( ',', pos );
      const string name = spec.substr( pos, (comma == string::npos) ? string::npos : (comma - pos) );
      if( !name.empty() )
        only.insert( name );
      if( comma == string::npos )
        break;
      pos = comma + 1;
    }
    if( only.count( "all" ) )
    {
      only.clear();
      BOOST_TEST_MESSAGE( "  (INTERSPEC_LINE_AB_SCENARIOS=all: the whole scenario matrix)" );
    }else
    {
      BOOST_TEST_MESSAGE( "  (restricted to " << only.size() << " scenario(s) by INTERSPEC_LINE_AB_SCENARIOS)" );
    }
  }else
  {
    BOOST_TEST_MESSAGE( "  (quick subset of " << only.size() << " scenarios; INTERSPEC_LINE_AB_SCENARIOS=all for the matrix)" );
  }

  // A misspelled name would otherwise run nothing and pass vacuously.
  for( const string &name : only )
  {
    bool known = false;
    for( const Scenario &s : scenarios() )
      known = known || (s.name == name);
    BOOST_REQUIRE_MESSAGE( known, "INTERSPEC_LINE_AB_SCENARIOS names an unknown scenario '" << name << "'" );
  }

  BOOST_TEST_MESSAGE( "  line/element - 1 (%), and CPU s per energy (element | line):" );
  for( const Scenario &s : scenarios() )
  {
    if( !only.empty() && !only.count(s.name) )
      continue;
    for( const bool transparent : { false, true } )
    {
      ostringstream row;
      row << "  " << left << setw(28) << (s.name + (transparent ? " (transp)" : "")) << right;
      for( const double e : energies )
      {
        size_t evals_e = 0, evals_l = 0;
        double err_e = 0.0, err_l = 0.0;
        const double t0 = cpu_now();
        const double elem = interspec_volumetric_eff( det, s, e, det.mc_transfer, -1, -1.0, -1.0,
                                                      transparent, &evals_e, &err_e );
        const double t1 = cpu_now();
        const double line = interspec_volumetric_eff( det, s, e, det.mc_transfer, -1, -1.0, -1.0,
                                                      transparent, &evals_l, &err_l, num_lines );
        const double t2 = cpu_now();
        elem_cpu += (t1 - t0);
        line_cpu += (t2 - t1);
        BOOST_REQUIRE( (elem > 0.0) && (line > 0.0) );
        const double rel = 100.0*(line/elem - 1.0);
        row << "  " << fixed << setprecision(0) << setw(5) << e << ":" << showpos << setprecision(2)
            << setw(6) << rel << "%" << noshowpos << " (" << setprecision(2) << (t1-t0) << "|"
            << (t2-t1) << "s)";
        if( fabs(rel) > worst )
        {
          worst = fabs( rel );
          worst_where = s.name + (transparent ? " (transp)" : "") + " @ " + to_string( e ) + " keV";
        }
      }//for( energies )
      BOOST_TEST_MESSAGE( row.str() );
    }//for( transparent )
  }//for( scenarios )

  ostringstream tail;
  tail << "  worst |line/element - 1|: " << fixed << setprecision(3) << worst << "% (" << worst_where
       << ");  total CPU element " << setprecision(1) << elem_cpu << " s, line " << line_cpu
       << " s (" << num_lines << " lines; the line set is rebuilt per call here - a fit builds it once)";
  BOOST_TEST_MESSAGE( tail.str() );

  // The gate is the two quadratures' COMBINED discretisation, measured rather than hoped for:
  //  the element path sits within 0.13% of its own 8192-ray reference (ApertureRayConvergence) and
  //  the line path within 0.21% of its 2^18-line one (LineCountConvergence), and the worst row
  //  measured 0.44% over the full matrix (transparent box at 122 keV, 2026-09-03; 0.29% over the
  //  quick subset).  0.75% leaves about one part in three of headroom over that - enough for
  //  platform floating-point drift, not enough to hide a real disagreement, which for two
  //  quadratures of the same integral would show up as a systematic trend across a scenario family
  //  rather than one loose row.
  BOOST_CHECK_MESSAGE( worst < 0.75,
                       "line and element quadratures disagree by " << worst << "% at " << worst_where );
}//BOOST_AUTO_TEST_CASE( LineVsElementScenarioMatrix )


/** LINE COUNT - how many lines the line path needs: value against a 2^18-line reference, on the
 contact rows (the hardest: short chords, steep prefactor), so the production default can be set
 where the change drops below 0.1%. */
BOOST_AUTO_TEST_CASE( LineCountConvergence )
{
  using namespace GammaInteractionCalc;
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );

  const std::vector<int> counts = { 1 << 12, 1 << 13, 1 << 14, 1 << 15, 1 << 16, 1 << 17 };
  const int reference = 1 << 18;

  double worst_64k = 0.0;
  string worst_where;
  for( const char *name : { "small-near-dense", "large-near-dense", "box-large-near-dense",
                            "shielded-near-dense", "large-far-dense" } )
  {
    const Scenario s = find_scenario( name );
    for( const double e : { 60.0, 661.7 } )
    {
      const double ref = interspec_volumetric_eff( det, s, e, det.mc_transfer, -1, -1.0, -1.0, false,
                                                   nullptr, nullptr, reference );
      ostringstream row;
      row << "  " << left << setw(22) << s.name << right << " @ " << setw(6) << fixed
          << setprecision(1) << e << " keV:";
      for( const int n : counts )
      {
        const double v = interspec_volumetric_eff( det, s, e, det.mc_transfer, -1, -1.0, -1.0, false,
                                                   nullptr, nullptr, n );
        const double rel = 100.0*(v/ref - 1.0);
        row << "  n=" << n << ":" << showpos << setprecision(3) << rel << "%" << noshowpos;
        if( (n == (1 << 16)) && (fabs(rel) > worst_64k) )
        {
          worst_64k = fabs( rel );
          worst_where = s.name + " @ " + to_string( e ) + " keV";
        }
      }
      BOOST_TEST_MESSAGE( row.str() );
    }
  }
  BOOST_TEST_MESSAGE( "  worst deviation of 65536 lines from the " << reference << "-line reference: "
                      << fixed << setprecision(3) << worst_64k << "% (" << worst_where << ")" );
  BOOST_CHECK_MESSAGE( worst_64k < 0.3, "65536 lines are not converged: " << worst_64k << "% at " << worst_where );
}//BOOST_AUTO_TEST_CASE( LineCountConvergence )


/** NESTED and MULTI-SHELL sources: line vs element on configurations the scenario matrix has none of.

 `LineVsElementScenarioMatrix` only ever builds ONE source shell at index 0, so two real
 Activity/Shielding configurations go unexercised by it:
   - a HOLLOW self-attenuating source (a non-source core inside the emitting shell), where the line
     path has to split each chord into the two pieces either side of the core and attenuate the far
     piece through it, and the element path takes a completely different route - it tiles the shell
     into smooth sub-domains (`split_source_subdomains`) to keep the inner-void indicator from
     wrecking the adaptive rule;
   - TWO self-attenuating materials in one problem, which produces two independent calculator groups
     with different `m_materialIndex`, hence two different line caches that must not be confused
     with each other.

 Both are integrated here through the production dispatcher, so the grouping and the per-material
 cache keying are under test along with the arithmetic.
 */
BOOST_AUTO_TEST_CASE( LineVsElementNestedAndMultiShell )
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

  // Builds a calculator whose source is shell `src_index` of a nested stack.
  const auto make_calc = [&]( const GeometryType geom,
                              const vector<pair<shared_ptr<const Material>,array<double,3>>> &shells,
                              const size_t src_index, const double energy,
                              const double standoff_cm ) -> DistributedSrcCalcT<double>
  {
    DistributedSrcCalcT<double> calc;
    calc.m_geometry = geom;
    calc.m_materialIndex = src_index;
    calc.m_attenuateForAir = false;
    calc.m_airTransLenCoef = 0.0;
    calc.m_isInSituExponential = false;
    calc.m_inSituRelaxationLength = -1.0;
    calc.m_srcVolumetricActivity = 1.0;
    calc.m_normalizeByVolume = false;
    calc.m_energy = energy;
    calc.m_effResponse = det.mc_transfer;
    calc.m_effMethod = ShieldingSourceFitCalc::VolumetricEffMethod::MCTransfer;

    // Detector far enough out that the outermost shell is comfortably inside it.
    const double det_radius = det.gd.transverse_half_extent() * cm;
    calc.m_detector = detector_geom_from_config<double>( geom, standoff_cm*cm, det_radius, 0.0 );

    for( const auto &sh : shells )
    {
      DistributedSrcCalcT<double>::ShellInfo info;
      for( int i = 0; i < 3; ++i )
        info.dims[i] = sh.second[i]*cm;
      info.trans_len_coef = transmition_length_coefficient( sh.first.get(),
                                                            static_cast<float>(energy) );
      info.type = ShellType::Material;
      calc.m_shells.push_back( info );
    }
    return calc;
  };

  const auto run = [&]( DistributedSrcCalcT<double> calc, const VolumetricIntegrator path,
                        const int num_lines ) -> double
  {
    integrate_on_path( calc, path, num_lines );
    return calc.integral;
  };

  double worst = 0.0;
  string worst_where;
  const int num_lines = 1 << 16;

  // ---- (a) HOLLOW source: steel core, emitting water shell, steel outer shield ----
  BOOST_TEST_MESSAGE( "  hollow source (non-source core inside the emitting shell):" );
  {
    struct Case { const char *name; GeometryType geom; array<double,3> in, mid, out; double standoff; };
    const vector<Case> cases = {
      { "cylEnd hollow", GeometryType::CylinderEndOn,  {1.0,0.8,0}, {3.0,2.0,0}, {3.3,2.3,0}, 3.0 },
      { "cylSide hollow", GeometryType::CylinderSideOn, {1.0,1.5,0}, {2.5,3.0,0}, {2.8,3.3,0}, 4.0 },
      { "rect hollow",   GeometryType::Rectangular,    {1.0,0.8,0.6}, {2.5,2.0,1.5}, {2.8,2.3,1.8}, 4.0 },
      { "sphere hollow", GeometryType::Spherical,      {1.0,0,0}, {2.5,0,0}, {2.8,0,0}, 4.0 },
    };
    for( const Case &c : cases )
    {
      for( const double e : { 60.0, 661.7 } )
      {
        const DistributedSrcCalcT<double> calc = make_calc( c.geom,
              { {steel,c.in}, {water,c.mid}, {steel,c.out} }, 1 /*source is the middle shell*/,
              e, c.standoff );
        const double elem = run( calc, VolumetricIntegrator::Element, -1 );
        const double line = run( calc, VolumetricIntegrator::Line, num_lines );
        BOOST_REQUIRE( elem > 0.0 );
        BOOST_REQUIRE( line > 0.0 );
        const double rel = 100.0*(line/elem - 1.0);
        ostringstream o;
        o << "    " << left << setw(16) << c.name << right << " @ " << setw(6) << fixed
          << setprecision(1) << e << " keV:  element " << scientific << setprecision(4) << elem
          << "  line " << line << "  (" << fixed << showpos << setprecision(2) << rel << "%"
          << noshowpos << ")";
        BOOST_TEST_MESSAGE( o.str() );
        // Spheres are gated like every other geometry as of 2026-09-03: `eval_spherical` used to
        //  apply the response through a single CENTRE RAY (the per-ray kernel had only ever been
        //  written for cylinders and rectangles), so this row was reported rather than gated
        //  because agreement would have pinned that defect.  It now runs the same aperture fan the
        //  other geometries do - see LineVsElementSphericalSource in test_VolumetricLinePath.cpp,
        //  which is where spheres are covered systematically.
        if( fabs(rel) > worst ){ worst = fabs(rel); worst_where = string(c.name) + " @ " + to_string(e); }
      }
    }
  }

  // ---- (b) TWO self-attenuating materials, integrated together through the dispatcher ----
  BOOST_TEST_MESSAGE( "  two self-attenuating shells in one problem (both source at once):" );
  {
    const vector<pair<shared_ptr<const Material>,array<double,3>>> stack =
        { {water,{1.5,1.2,0}}, {steel,{3.0,2.5,0}} };

    for( const double e : { 60.0, 661.7 } )
    {
      // One calculator per source shell, exactly as build_volumetric_calculators would emit.
      vector<unique_ptr<DistributedSrcCalcT<double>>> both;
      for( size_t src = 0; src < 2; ++src )
      {
        DistributedSrcCalcT<double> c = make_calc( GeometryType::CylinderEndOn, stack, src, e, 4.0 );
        attach_line_cache( c, num_lines );
        both.push_back( std::make_unique<DistributedSrcCalcT<double>>( c ) );
      }
      // The two source shells must have been given DIFFERENT line caches.
      BOOST_CHECK( both[0]->m_lineCache.get() != both[1]->m_lineCache.get() );

      {
        const ScopedVolumetricIntegratorOverride force( VolumetricIntegrator::Line );
        integrate_volumetric_calculators<double>( both, true );
      }

      for( size_t src = 0; src < 2; ++src )
      {
        DistributedSrcCalcT<double> ref = make_calc( GeometryType::CylinderEndOn, stack, src, e, 4.0 );
        const double elem = run( ref, VolumetricIntegrator::Element, -1 );
        const double line = both[src]->integral;
        BOOST_REQUIRE( elem > 0.0 );
        BOOST_REQUIRE( line > 0.0 );
        const double rel = 100.0*(line/elem - 1.0);
        ostringstream o;
        o << "    source shell " << src << " @ " << setw(6) << fixed << setprecision(1) << e
          << " keV:  element " << scientific << setprecision(4) << elem << "  line " << line
          << "  (" << fixed << showpos << setprecision(2) << rel << "%" << noshowpos << ")";
        BOOST_TEST_MESSAGE( o.str() );
        if( fabs(rel) > worst ){ worst = fabs(rel); worst_where = "multi-shell src " + to_string(src); }
      }
    }
  }

  BOOST_TEST_MESSAGE( "  worst |line/element - 1|: " << fixed << setprecision(3) << worst
                      << "% (" << worst_where << ")" );
  BOOST_CHECK_MESSAGE( worst < 1.0,
    "line and element disagree by " << worst << "% at " << worst_where
    << " - nested/multi-shell sources are not being integrated the same way" );
}//BOOST_AUTO_TEST_CASE( LineVsElementNestedAndMultiShell )


/** Is the objective CONTINUOUS in a source dimension across a line-set rebuild?

 The line set is a fixed importance-sampling proposal aimed at the (padded) source, and
 `VolumetricLineCache::matches` holds one set across a +-20% dimension window before rebuilding.
 A rebuild is the thing to be suspicious of: if it produced an INDEPENDENT quadrature the objective
 would step by that quadrature's own discretisation error (~0.2%) as a fitted dimension crossed the
 window edge, and Levenberg-Marquardt would be comparing a predicted reduction against an actual one
 that contains a jump it cannot model.

 It should NOT be independent - the Halton indices are the same, the hull points depend on the
 detector geometry rather than the source, and each line's direction is aimed at a point that moves
 CONTINUOUSLY with the proposal dimensions - so a rebuilt set is a small deformation of the old one
 and the estimate should carry across.  This measures that rather than assuming it: sweep a source
 radius finely across the rebuild boundary, driving the cache exactly as a fit does, and compare
 each step against its neighbour and against the (rebuild-free) element path.

 KNOWN-FAILING, and marked as such so a REGRESSION here is still visible: the answer is that a
 rebuild does NOT carry across.  Measured 2026-09-03 and unchanged since: a second difference of
 1.6e-3 in the line/element ratio at the window edge, i.e. the objective steps by ~0.16% when a
 fitted dimension crosses it, which Levenberg-Marquardt cannot model.  The proposal is re-aimed by
 the whole width of the window at once, so the sets either side really are independent quadratures.
 The fix is not a wider window - a fit crosses it eventually - but to aim the proposal at the
 dimension parameters' upper BOUNDS once per fit, so it spans the search domain and never rebuilds
 (equivalently, a source-SCALED proposal whose aim points deform continuously with the dimensions).
 Written up in TODO.md; only fits with FREE source dimensions can meet it.
 */
BOOST_AUTO_TEST_CASE( LineCacheRebuildContinuity, * boost::unit_test::expected_failures(1) )
{
  using namespace GammaInteractionCalc;
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );
  const double cm = PhysicalUnits::cm;

  const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  const shared_ptr<const Material> water = matdb->material( "Water" );
  BOOST_REQUIRE( water );

  const double energy = 661.7;
  const int num_lines = 1 << 16;
  const double half_len = 2.0;

  // Stand in for ShieldingSourceChi2Fcn::volumetricLineCache: hold one cache and rebuild it only
  //  when `matches` says the current dimensions have left its window.
  std::shared_ptr<const VolumetricLineCache> cache;
  int rebuilds = 0;

  const auto eff_at = [&]( const double radius_cm ) -> double
  {
    DistributedSrcCalcT<double> calc;
    calc.m_geometry = GeometryType::CylinderEndOn;
    calc.m_materialIndex = 0;
    calc.m_attenuateForAir = false;
    calc.m_isInSituExponential = false;
    calc.m_inSituRelaxationLength = -1.0;
    calc.m_srcVolumetricActivity = 1.0;
    calc.m_normalizeByVolume = true;      //fit-like: divide by the (changing) source volume
    calc.m_energy = energy;
    calc.m_effResponse = det.mc_transfer;
    calc.m_effMethod = ShieldingSourceFitCalc::VolumetricEffMethod::MCTransfer;
    const double det_radius = det.gd.transverse_half_extent() * cm;
    calc.m_detector = detector_geom_from_config<double>( GeometryType::CylinderEndOn, 4.0*cm,
                                                         det_radius, 0.0 );
    DistributedSrcCalcT<double>::ShellInfo info;
    info.dims = { radius_cm*cm, half_len*cm, 0.0 };
    info.trans_len_coef = transmition_length_coefficient( water.get(),
                                                          static_cast<float>(energy) );
    info.type = ShellType::Material;
    calc.m_shells.push_back( info );

    const std::array<double,3> dims = { radius_cm*cm, half_len*cm, 0.0 };
    const std::array<double,3> det_pos = { calc.m_detector.position[0], calc.m_detector.position[1],
                                           calc.m_detector.position[2] };
    const std::array<double,3> det_axis = { calc.m_detector.axis[0], calc.m_detector.axis[1],
                                            calc.m_detector.axis[2] };
    if( !cache || !cache->matches( det.mc_transfer.get(), GeometryType::CylinderEndOn, 0, dims,
                                   det_pos, det_axis, 0.0, num_lines ) )
    {
      cache = build_volumetric_line_cache( det.mc_transfer, GeometryType::CylinderEndOn, 0, dims,
                                           det_pos, det_axis, 0.0, num_lines );
      ++rebuilds;
    }
    calc.m_lineCache = cache;

    std::vector<std::unique_ptr<DistributedSrcCalcT<double>>> v;
    v.push_back( std::make_unique<DistributedSrcCalcT<double>>( calc ) );
    {
      const ScopedVolumetricIntegratorOverride force( VolumetricIntegrator::Line );
      integrate_volumetric_calculators<double>( v, true );
    }
    return v.front()->integral;
  };

  // Start at 2.0 cm; the window is +-20%, so 2.4 cm and 1.6 cm are the edges.  Sweep well past both.
  const double r0 = 2.0;
  vector<double> radii;
  for( double r = 1.50; r < 3.001; r += 0.025 )
    radii.push_back( r );

  cache.reset();
  (void)eff_at( r0 );          //seed the cache at the centre of the window, as a fit would
  const int rebuilds_before = rebuilds;

  vector<double> line_vals, elem_vals;
  for( const double r : radii )
    line_vals.push_back( eff_at( r ) );

  BOOST_TEST_MESSAGE( "  swept radius " << radii.front() << " -> " << radii.back() << " cm in "
                      << radii.size() << " steps; line-set rebuilds: "
                      << (rebuilds - rebuilds_before) );
  BOOST_CHECK_MESSAGE( (rebuilds - rebuilds_before) > 0,
                       "the sweep never crossed a rebuild boundary - it is not testing anything" );

  // The element path for the same sweep: no proposal, so any structure here is physical.
  {
    std::shared_ptr<const VolumetricLineCache> keep = cache;
    cache.reset();
    for( const double r : radii )
    {
      DistributedSrcCalcT<double> calc;
      calc.m_geometry = GeometryType::CylinderEndOn;
      calc.m_materialIndex = 0;
      calc.m_attenuateForAir = false;
      calc.m_isInSituExponential = false;
      calc.m_inSituRelaxationLength = -1.0;
      calc.m_srcVolumetricActivity = 1.0;
      calc.m_normalizeByVolume = true;
      calc.m_energy = energy;
      calc.m_effResponse = det.mc_transfer;
      calc.m_effMethod = ShieldingSourceFitCalc::VolumetricEffMethod::MCTransfer;
      const double det_radius = det.gd.transverse_half_extent() * cm;
      calc.m_detector = detector_geom_from_config<double>( GeometryType::CylinderEndOn, 4.0*cm,
                                                           det_radius, 0.0 );
      DistributedSrcCalcT<double>::ShellInfo info;
      info.dims = { r*cm, half_len*cm, 0.0 };
      info.trans_len_coef = transmition_length_coefficient( water.get(),
                                                            static_cast<float>(energy) );
      info.type = ShellType::Material;
      calc.m_shells.push_back( info );
      integrate_on_path( calc, VolumetricIntegrator::Element, -1 );
      elem_vals.push_back( calc.integral );
    }
    cache = keep;
  }

  // Second difference of the RATIO line/element: the physical curvature cancels, so what is left is
  //  the proposal's own contribution.  A rebuild that reset the quadrature would show up here as a
  //  single large spike at the window edge.
  double worst_jump = 0.0;
  size_t worst_i = 0;
  for( size_t i = 1; i + 1 < radii.size(); ++i )
  {
    const double a = line_vals[i-1]/elem_vals[i-1];
    const double b = line_vals[i]  /elem_vals[i];
    const double c = line_vals[i+1]/elem_vals[i+1];
    const double second = fabs( c - 2.0*b + a );
    if( second > worst_jump ){ worst_jump = second; worst_i = i; }
  }

  ostringstream o;
  o << "  worst second difference of line/element over the sweep: " << scientific
    << setprecision(3) << worst_jump << " at r = " << fixed << setprecision(3) << radii[worst_i]
    << " cm (ratio there " << setprecision(5) << line_vals[worst_i]/elem_vals[worst_i] << ")";
  BOOST_TEST_MESSAGE( o.str() );

  // An independent re-draw of the proposal would put a ~2e-3 step in the ratio, hence a second
  //  difference of the same order.  Anything at or below 1e-3 means the rebuild carried across.
  BOOST_CHECK_MESSAGE( worst_jump < 1.0e-3,
    "the line/element ratio jumps by " << worst_jump << " at r = " << radii[worst_i]
    << " cm - a line-set rebuild is resetting the quadrature instead of deforming it, which puts a"
    " step in the objective that Levenberg-Marquardt cannot model" );
}//BOOST_AUTO_TEST_CASE( LineCacheRebuildContinuity )


/** DIAGNOSTIC (developer-only): which side of the hollow-rectangle line-vs-element gap is converged?

 WAS: `LineVsElementNestedAndMultiShell` found the two paths ~2.4% apart on a HOLLOW box at 60 keV.
 RESOLVED 2026-09-03 - it was neither side's convergence but a real bug in `rectangle_intersections_imp`,
 which reported a MISS for any ray clipping a corner of the inner core, so the element path walked
 those rays as if the core were not there.  The rows now agree to 0.27%.  Kept because refining each
 side against itself - the line path in its number of lines, the element path in its requested
 relative error - is the right first move whenever a new gap appears: whichever one moves is the one
 that was not converged.
 */
BOOST_AUTO_TEST_CASE( NestedRectConvergenceProbe, * boost::unit_test::disabled() )
{
  using namespace GammaInteractionCalc;
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  BOOST_REQUIRE( det.mc_transfer );
  const double cm = PhysicalUnits::cm;

  const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  const shared_ptr<const Material> steel = matdb->material( "Stainless steel SS-304" );
  const shared_ptr<const Material> water = matdb->material( "Water" );
  BOOST_REQUIRE( steel && water );

  const auto make_calc = [&]( const double energy ) -> DistributedSrcCalcT<double>
  {
    DistributedSrcCalcT<double> calc;
    calc.m_geometry = GeometryType::Rectangular;
    calc.m_materialIndex = 1;
    calc.m_attenuateForAir = false;
    calc.m_isInSituExponential = false;
    calc.m_inSituRelaxationLength = -1.0;
    calc.m_srcVolumetricActivity = 1.0;
    calc.m_normalizeByVolume = false;
    calc.m_energy = energy;
    calc.m_effResponse = det.mc_transfer;
    calc.m_effMethod = ShieldingSourceFitCalc::VolumetricEffMethod::MCTransfer;
    const double det_radius = det.gd.transverse_half_extent() * cm;
    calc.m_detector = detector_geom_from_config<double>( GeometryType::Rectangular, 4.0*cm,
                                                         det_radius, 0.0 );
    const array<array<double,3>,3> dims = { array<double,3>{1.0,0.8,0.6},
                                            array<double,3>{2.5,2.0,1.5},
                                            array<double,3>{2.8,2.3,1.8} };
    const array<shared_ptr<const Material>,3> mats = { steel, water, steel };
    for( size_t i = 0; i < 3; ++i )
    {
      DistributedSrcCalcT<double>::ShellInfo info;
      for( int k = 0; k < 3; ++k )
        info.dims[k] = dims[i][k]*cm;
      info.trans_len_coef = transmition_length_coefficient( mats[i].get(),
                                                            static_cast<float>(energy) );
      info.type = ShellType::Material;
      calc.m_shells.push_back( info );
    }
    return calc;
  };

  for( const double energy : { 60.0, 661.7 } )
  {
    BOOST_TEST_MESSAGE( "  hollow rect @ " << energy << " keV" );

    ostringstream lrow;
    lrow << "    line ";
    for( const int n : { 1<<14, 1<<16, 1<<18, 1<<20 } )
    {
      DistributedSrcCalcT<double> calc = make_calc( energy );
      integrate_on_path( calc, VolumetricIntegrator::Line, n );
      lrow << "  n=" << n << ":" << scientific << setprecision(6) << calc.integral;
    }
    BOOST_TEST_MESSAGE( lrow.str() );

    ostringstream erow;
    erow << "    elem ";
    for( const double epsrel : { 1.0e-3, 1.0e-4, 1.0e-5 } )
    {
      DistributedSrcCalcT<double> calc = make_calc( energy );
      self_shielding_integration_imp<double>( calc, epsrel, 400000000 );
      erow << "  eps=" << scientific << setprecision(0) << epsrel << ":" << setprecision(6)
           << calc.integral << " (" << calc.m_num_evals << " evals, est "
           << setprecision(1) << calc.m_est_rel_error << ")";
    }
    BOOST_TEST_MESSAGE( erow.str() );
  }
}//BOOST_AUTO_TEST_CASE( NestedRectConvergenceProbe )


/** RUNG 6 - OFF-AXIS cylinders against Monte Carlo.

 Every other scenario in the matrix puts the detector on the source axis, where an end-on cylinder is
 azimuthally symmetric and the integration collapses to 2D: all its elements live in one half-plane.
 A radial offset forces the full 3D integral and puts elements at every azimuth, which is the only
 cylinder configuration where the aperture fan's AZIMUTHAL placement enters rather than just its
 polar orientation.  A SIDE-ON cylinder adds the other thing the matrix never had: the detector on
 +x rather than +z, looking at the curved surface - a first-class production geometry, with its own
 detector placement and its own in-situ depth convention, that no Monte-Carlo-backed test had ever
 exercised.  Until this case existed, off-axis cylinders were covered only by a
 self-consistency invariant (`OffAxisEndOnAzimuthalSymmetry` in test_VolumetricIntegration, which
 checks that an x offset equals a y offset - and runs the FLAT-DISK path, so it never touches the
 per-ray kernel at all) and by the ray-march identity in `PerRayKernelIdentityVsRayMarch` (which
 shares its aperture with the thing under test, so it cannot see an orientation error).

 Prints paste-ready TruthRow initializers for the off-axis scenarios plus the legacy/fixed
 comparison, so the frame fix is measured in the one geometry that isolates the azimuth.
 */
BOOST_AUTO_TEST_CASE( Rung6_OffAxisAndSideOnCylinderTruth, * boost::unit_test::disabled() )
{
  using namespace GammaInteractionCalc;
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  McCache cache;
  const bool saved = sm_aperture_frame_legacy_origin;

  const shared_ptr<ceelo::DetectorResponse> far_anchor
        = ladder_anchor_response( det, cache, 25.0, anchor_energies() );

  // Both orientations the on-axis end-on matrix never reached: a radial offset (elements at general
  //  azimuth) and a side-on cylinder (the detector on +x rather than +z, looking at the curved
  //  surface).  Both are always 3D integrations.
  vector<Scenario> offaxis;
  for( const Scenario &s : scenarios() )
  {
    if( (s.offset_cm != 0.0) || (s.shape == Shape::CylinderSideOn) )
      offaxis.push_back( s );
  }
  BOOST_REQUIRE( !offaxis.empty() );

  cout << "  // Off-axis and side-on cylinders - the orientations the on-axis end-on matrix,\n"
          "  //  which collapses to a 2D integration, could never reach.\n" << flush;
  BOOST_TEST_MESSAGE( "scenario                       E(keV)   tau    MC eff      sigma   |"
                      " legacy/MC-1  fixed/MC-1   CPU s" );
  for( const Scenario &s : offaxis )
  {
    ScenarioMcMaterials mats;
    ceelo::EfficiencyCalculator mc_calc;
    configure_scenario_mc( mc_calc, det, s, mats );

    for( const double e : scenario_energies() )
    {
      const McResult mc = run_mc( mc_calc, cache, scenario_mc_key( s ), e );
      cout << "  { \"" << s.name << "\", " << setprecision(12) << e << ", "
           << setprecision(10) << mc.eff << ", " << setprecision(4) << mc.eff*mc.frac_sigma << " },"
           << (mc.hit_cap() ? "  // stopped on a cap" : "") << "\n" << flush;

      const clock_t t0 = std::clock();
      sm_aperture_frame_legacy_origin = true;
      const double legacy = interspec_volumetric_eff( det, s, e, far_anchor );
      sm_aperture_frame_legacy_origin = false;
      const double fixed_v = interspec_volumetric_eff( det, s, e, far_anchor );
      sm_aperture_frame_legacy_origin = saved;
      const double cpu = 0.5*static_cast<double>(std::clock() - t0)/CLOCKS_PER_SEC;

      ostringstream o;
      o << left << setw(31) << s.name << right << fixed << setprecision(1) << setw(7) << e
        << setprecision(2) << setw(7) << scenario_optical_depth( s, e )
        << "  " << scientific << setprecision(4) << mc.eff << fixed << "  " << setprecision(2)
        << 100.0*mc.frac_sigma << "%" << (mc.hit_cap() ? "CAP" : "   ")
        << " | " << pct( legacy/mc.eff - 1.0 ) << " " << pct( fixed_v/mc.eff - 1.0 )
        << "  " << setprecision(1) << setw(7) << cpu;
      BOOST_TEST_MESSAGE( o.str() );
    }//for( energies )
  }//for( off-axis scenarios )

  sm_aperture_frame_legacy_origin = saved;
}//BOOST_AUTO_TEST_CASE( Rung6_OffAxisAndSideOnCylinderTruth )


/** RUNG 7 - record the CENTRE-ANCHORED transfer curves the committed comparison needs.

 `VolumetricNearFieldVsTruth` grounds its transfer on one recorded Monte-Carlo point-source curve at
 25 cm, because the committed test runs no Monte Carlo of its own.  That single choice is where most
 of its tolerance goes: the contact scenarios sit at 1.4-3 cm, so the transfer extrapolates inward by
 a factor of ten, and rung 2 measured the consequence with NO attenuating material anywhere - the same
 scenarios read within +-0.7% at every energy against a centre anchor but +4.6 to +6.4% at 344-1332
 keV against the 25 cm one.  The allowance was therefore paying for the anchor, not for the volume
 integral it is supposed to be gating.

 So record a curve per distinct source-centre distance as well.  A centre anchor is a bare point in
 vacuum on axis, so it depends only on that distance and the energy grid - not on the shape, matrix or
 shield of the scenario whose centre happens to be there - which is why one curve serves every
 scenario at that distance (a box and its cylinder twin share one exactly).

 Cheap: bare points at contact are high-efficiency, so each is seconds of Monte Carlo, and they all
 come from the ladder's cache.
 */
BOOST_AUTO_TEST_CASE( Rung7_CentreAnchorTruth, * boost::unit_test::disabled() )
{
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  McCache cache;

  // Every distinct centre distance in the matrix, rounded to the micron so a box and its cylinder
  //  twin collapse onto one entry.
  std::set<long long> keys;
  vector<double> dists;
  for( const Scenario &s : scenarios() )
  {
    const double d = scenario_center_distance_cm( s );
    if( keys.insert( llround( d*1.0e6 ) ).second )
      dists.push_back( d );
  }
  std::sort( begin(dists), end(dists) );
  BOOST_TEST_MESSAGE( "distinct source-centre distances: " << dists.size() );

  vector<unique_ptr<ceelo::Material>> owned;
  ceelo::EfficiencyCalculator pt;
  ceelo::ResponseGenerator::configure_calculator( pt, det.gd, owned );

  cout << "static const std::vector<CentreAnchorRow> sm_centre_anchor = {\n" << flush;
  for( const double d : dists )
  {
    pt.set_point_source( Eigen::Vector3d( 0.0, 0.0, -(det.endcap_front_offset_cm + d) ) );
    ostringstream k;
    k << "point(r=0,ze=" << setprecision(10) << d << ");bare";
    for( const double e : scenario_energies() )
    {
      const McResult r = run_mc( pt, cache, k.str(), e );
      cout << "  { " << setprecision(10) << d << ", " << setprecision(12) << e << ", "
           << setprecision(10) << r.eff << ", " << setprecision(4) << r.frac_sigma << " },"
           << (r.hit_cap() ? "  // stopped on a cap" : "") << "\n" << flush;
    }
  }
  cout << "};\n" << flush;
}//BOOST_AUTO_TEST_CASE( Rung7_CentreAnchorTruth )


/** RUNG 8 - SPHERICAL sources against Monte Carlo.

 The scenario matrix and its truth bank are cylinders and boxes, so until now NOTHING had ever
 compared a spherical volumetric source against Monte Carlo - the sphere claims in
 `LineVsElementSphericalSource` are the two model quadratures agreeing with each other, which says
 they are consistent, not that they are right.  This rung supplies the missing leg.

 Each row is grounded on its own CENTRE-ANCHORED transfer, built here from a bare point-source MC at
 the sphere's centre distance (the argument is spelled out at Rung7): anchoring at 25 cm and
 comparing a source at 4 cm would charge the volume integral for the transfer's inward extrapolation,
 which rung 2 measured at +4.6 to +6.4%.  With a centre anchor, what is left is the volumetric model.

 HOLLOW ROW - read the geometry.  CeeLo's spherical source takes an inner radius whose interior is a
 non-attenuating VOID, so the "steel core" of the model-vs-model test is not expressible here; the
 row below therefore puts a genuine void inside the emitting shell and the InterSpec side matches it
 with a zero attenuation coefficient on the core shell.  It still exercises the two-piece chord
 split, which is the part of the code the hollow case is there to test.

 Developer-only (real Monte Carlo).  Run with --cachedir= so the rows are paid for once.
 */
BOOST_AUTO_TEST_CASE( Rung8_SphericalSourceTruth, * boost::unit_test::disabled() )
{
  using namespace GammaInteractionCalc;
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  McCache cache;
  const double cm = PhysicalUnits::cm;

  const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );
  const shared_ptr<const Material> steel = matdb->material( "Stainless steel SS-304" );
  const shared_ptr<const Material> water = matdb->material( "Water" );
  BOOST_REQUIRE( steel && water );

  struct Case
  {
    const char *name;
    shared_ptr<const Material> mat;
    double outer_cm;
    double inner_cm;     ///< 0 = solid; > 0 = emitting shell around a NON-attenuating void
    double centre_cm;    ///< endcap front to the sphere centre
  };
  const vector<Case> cases = {
    { "solid steel, contact",     steel, 3.0, 0.0, 4.0 },
    { "solid water, contact",     water, 3.0, 0.0, 4.0 },
    { "solid water, far",         water, 3.0, 0.0, 50.0 },
    { "hollow steel on a void",   steel, 3.0, 1.0, 4.0 },
  };
  const vector<double> energies = { 60.0, 661.7 };
  const vector<double> anchor_energies = { 60.0, 88.0, 122.0, 344.0, 661.7, 1332.5 };
  const int num_lines = 1 << 16;

  // The source attenuation coefficient carries the FEP-WINDOW credit, as the truth bank does: a
  //  photon Compton-scattering by less than the peak's half-width is still in the full-energy peak,
  //  so charging it mu_TOTAL under-counts the efficiency.  The window is not free - CeeLo scores on
  //  a +-half-width, so it is FWHM/2 of the detector the transfer is anchored on, and for this
  //  GEM35-70 that is FWHM = sqrt(0.359 + 0.00230*E).  Run with mu_total instead and the light
  //  matrix reads ~2% low at 60 keV, where the Compton fraction is largest (measured 2026-09-04);
  //  a dense matrix at 60 keV is photoelectric-dominated and barely moves.
  const auto fep_window_keV = []( const double e ) -> double {
    return 0.5*std::sqrt( 0.359 + 0.00230*e );
  };

  // A bare point in vacuum at `d_cm`, as the transfer's anchor curve; cached by distance.
  map<long long,shared_ptr<const ceelo::DetectorResponse>> anchored;
  const auto centre_anchor = [&]( const double d_cm ) -> shared_ptr<const ceelo::DetectorResponse>
  {
    const long long key = llround( d_cm*1.0e6 );
    const map<long long,shared_ptr<const ceelo::DetectorResponse>>::const_iterator it = anchored.find( key );
    if( it != end(anchored) )
      return it->second;

    CeeLoUtils::TransferAnchor anchor;
    anchor.ref_distance_cm = d_cm;
    anchor.curve_derived = false;
    for( const double e : anchor_energies )
    {
      vector<unique_ptr<ceelo::Material>> owned;
      ceelo::EfficiencyCalculator pt;
      ceelo::ResponseGenerator::configure_calculator( pt, det.gd, owned );
      pt.set_point_source( Eigen::Vector3d( 0.0, 0.0, -(det.endcap_front_offset_cm + d_cm) ) );
      ostringstream k;
      k << "sph-anchor|point|d=" << setprecision(10) << d_cm;
      // 0.1%, not the 0.5% a bare anchor would normally get: the model prediction scales DIRECTLY
      //  with the anchor value at that energy, so the anchor's own statistics land in the quoted
      //  model/MC difference one-for-one.  At 0.5% the floor was 0.47% (anchor 0.40% + sphere MC
      //  0.25% in quadrature), which is bigger than most of the differences being reported; at 0.1%
      //  the floor is 0.27% and is dominated by the sphere MC, where it belongs.  A point source is
      //  cheap, so this costs seconds.
      const McResult r = run_mc( pt, cache, k.str(), e, 0.001, 4000000000ULL, 20260904, 900.0 );
      anchor.curve.energies_keV.push_back( e );
      anchor.curve.eff.push_back( r.eff );
      anchor.curve.frac_sigma.push_back( r.frac_sigma );
    }

    ostringstream nm;
    nm << "sphere centre anchor (d=" << fixed << setprecision(2) << d_cm << " cm)";
    shared_ptr<const ceelo::DetectorResponse> resp
          = CeeLoUtils::makeTransferResponse( det.gd, anchor, ceelo::AnchorCurve{}, nm.str() );
    BOOST_REQUIRE( resp );
    anchored[key] = resp;
    return resp;
  };

  BOOST_TEST_MESSAGE( "  sphere vs MC (centre-anchored transfer); model/MC - 1:" );
  BOOST_TEST_MESSAGE( "  case                            E(keV)      MC eff   MCsig |  element      line" );

  double worst = 0.0;
  string worst_where;
  for( const Case &c : cases )
  {
    const shared_ptr<const ceelo::DetectorResponse> resp = centre_anchor( c.centre_cm );

    for( const double e : energies )
    {
      // ---- Monte Carlo truth ----
      ScenarioMcMaterials mats;
      ceelo::EfficiencyCalculator mc_calc;
      ceelo::ResponseGenerator::configure_calculator( mc_calc, det.gd, mats.owned );
      mc_calc.set_spherical_source( Eigen::Vector3d( 0.0, 0.0, -(det.endcap_front_offset_cm + c.centre_cm) ),
                                    c.outer_cm, Eigen::Matrix3d::Identity(), c.inner_cm );
      mats.matrix.reset( new ceelo::Material( CeeLoUtils::to_ceelo_material( *c.mat ).to_material() ) );
      mc_calc.set_source_material( mats.matrix.get() );

      ostringstream key;
      key << "sphere(r=" << setprecision(10) << c.outer_cm << ",ri=" << c.inner_cm
          << ",d=" << c.centre_cm << ",mat=" << c.mat->name << ")";
      const McResult mc = run_mc( mc_calc, cache, key.str(), e, 0.0025, 16000000000ULL, 20260904, 1800.0 );
      BOOST_REQUIRE( mc.eff > 0.0 );

      // ---- InterSpec model, both quadratures ----
      DistributedSrcCalcT<double> calc;
      calc.m_geometry = GeometryType::Spherical;
      calc.m_materialIndex = (c.inner_cm > 0.0) ? 1 : 0;
      calc.m_attenuateForAir = false;
      calc.m_isInSituExponential = false;
      calc.m_inSituRelaxationLength = -1.0;
      calc.m_srcVolumetricActivity = 1.0;
      calc.m_normalizeByVolume = false;
      calc.m_energy = e;
      calc.m_effResponse = resp;
      calc.m_effMethod = ShieldingSourceFitCalc::VolumetricEffMethod::MCTransfer;
      calc.m_detector = detector_geom_from_config<double>( GeometryType::Spherical, c.centre_cm*cm,
                                              det.gd.transverse_half_extent()*cm, 0.0 );
      if( c.inner_cm > 0.0 )
      {
        DistributedSrcCalcT<double>::ShellInfo core;
        core.dims = { c.inner_cm*cm, 0.0, 0.0 };
        core.trans_len_coef = 0.0;    //a VOID, matching CeeLo's non-attenuating inner region
        core.type = ShellType::Material;
        calc.m_shells.push_back( core );
      }
      DistributedSrcCalcT<double>::ShellInfo src;
      src.dims = { c.outer_cm*cm, 0.0, 0.0 };
      src.trans_len_coef = fep_removal_coefficient( *c.mat, e, fep_window_keV(e) );
      src.type = ShellType::Material;
      calc.m_shells.push_back( src );

      const double vol_cm3 = (4.0/3.0)*PhysicalUnits::pi
                             * (c.outer_cm*c.outer_cm*c.outer_cm - c.inner_cm*c.inner_cm*c.inner_cm);

      DistributedSrcCalcT<double> elem = calc, line = calc;
      integrate_on_path( elem, VolumetricIntegrator::Element, -1 );
      integrate_on_path( line, VolumetricIntegrator::Line, num_lines );
      const double eff_elem = elem.integral / (vol_cm3*cm*cm*cm);
      const double eff_line = line.integral / (vol_cm3*cm*cm*cm);

      const double d_elem = 100.0*(eff_elem/mc.eff - 1.0);
      const double d_line = 100.0*(eff_line/mc.eff - 1.0);

      ostringstream o;
      o << "  " << left << setw(28) << c.name << right << fixed << setprecision(1) << setw(8) << e
        << "  " << scientific << setprecision(4) << mc.eff << fixed << "  " << setprecision(2)
        << setw(5) << 100.0*mc.frac_sigma << "%" << (mc.hit_cap() ? "CAP" : "   ")
        << " | " << showpos << setw(7) << setprecision(2) << d_elem << "%"
        << setw(9) << d_line << "%" << noshowpos;
      BOOST_TEST_MESSAGE( o.str() );

      // 3 sigma of the MC on top of the model allowance, as the truth bank does.  The ANCHOR's
      //  statistics are in here too - see centre_anchor - but at 0.1% they are well inside this.
      const double allow = 3.0 + 3.0*100.0*mc.frac_sigma;
      for( const pair<const char *,double> &m : { make_pair("element",d_elem), make_pair("line",d_line) } )
      {
        if( fabs(m.second) > worst ){ worst = fabs(m.second); worst_where = string(c.name) + " " + m.first; }
        BOOST_CHECK_MESSAGE( fabs(m.second) < allow,
          c.name << " @ " << e << " keV: the " << m.first << " path is " << m.second
                 << "% from MC (allowance " << allow << "%)" );
      }
    }//for( energies )
  }//for( cases )

  BOOST_TEST_MESSAGE( "  worst |model/MC - 1|: " << fixed << setprecision(2) << worst
                      << "% (" << worst_where << ")" );
}//BOOST_AUTO_TEST_CASE( Rung8_SphericalSourceTruth )


/** RUNG 9 - is the light-matrix residual really the FEP-WINDOW CREDIT, and if so WHY?

 Rung 8 leaves both quadratures reading ~1.1-1.5% LOW on a water sphere at 60 keV, at contact and at
 50 cm alike, while steel at the same energy is within 0.3%.  "Changing the credit changes the
 answer" is not evidence for the credit being the cause - it is nearly a tautology - so this rung
 tries to FALSIFY it, and then to pin the mechanism.

 The lever is that `f_win` is a property of (energy, window, material) ONLY.  It cannot depend on how
 far away the detector is or how big the source is.  So SOLVE for the credit each configuration
 demands - the f that makes model == MC - and see whether one value explains them all:

   - GEOMETRY: the same water sphere at 4, 6 and 50 cm.  Same material, same optical depth, very
     different solid angle and near-field regime.
   - OPTICAL DEPTH: water spheres of radius 1, 3 and 5 cm at one distance.  The credit acts through
     exp(-f*mu_Compton*chord), so an error in f shows up in proportion to tau; a geometric error
     does not.
   - MATERIAL: steel at 60 keV, photoelectric-dominated, so mu_Compton is a small share of mu_total.
   - TRANSPARENT TWIN for every geometry: mu_src = 0, so the credit is inert by construction and
     what is left is the geometry-and-response floor.  A transparent source is the same object
     whatever material it nominally is, so one twin serves every probe at that (radius, distance).
     The demanded f is solved on the residual AFTER that floor is removed - without this the contact
     rows are charged for a near-field response error that has nothing to do with attenuation.

 THE MECHANISM ON TRIAL.  The window is not a free parameter here: CeeLo's MC scores a full-energy
 event as |E_dep - E_src| < kDefaultFepWindowKeV, a HALF-width of 0.75 keV.  When the model is
 compared against THAT truth it must credit back the Compton scatters that same window keeps, i.e.
 it must use 0.75 keV - whereas the truth-bank convention (and rung 8) uses FWHM/2 of the anchoring
 detector, 0.35 keV at 60 keV.  Too narrow a window credits back too little, over-attenuates, and
 reads LOW in proportion to tau, which is the observed sign and the observed scaling.

 PREDICTIONS, stated before the numbers:
   1. the demanded f agrees across all six configurations once the floor is removed;
   2. it lands near kn_in_window_fraction(60 keV, 0.75 keV, water), NOT near the FWHM/2 value;
   3. re-running the model at the MC's own 0.75 keV window collapses the tau-proportional residual,
      leaving each row at its transparent floor.
 The credit explanation FAILS if the demanded f values disagree beyond their uncertainties, or if
 any exceeds 1 (f is a fraction of Compton scatters, so no physical credit could close the gap).

 Developer-only (real Monte Carlo); run with --cachedir=.
 */
BOOST_AUTO_TEST_CASE( Rung9_FepWindowCreditProbe, * boost::unit_test::disabled() )
{
  using namespace GammaInteractionCalc;
  set_data_dir();
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  const AngleDetector det = load_angle_detector();
  McCache cache;
  const double cm = PhysicalUnits::cm;

  const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );
  const shared_ptr<const Material> steel = matdb->material( "Stainless steel SS-304" );
  const shared_ptr<const Material> water = matdb->material( "Water" );
  BOOST_REQUIRE( steel && water );

  const double energy = 60.0;
  const int num_lines = 1 << 16;
  const vector<double> anchor_energies = { 60.0, 88.0, 122.0, 344.0, 661.7, 1332.5 };
  const double win_fwhm_half = 0.5*std::sqrt( 0.359 + 0.00230*energy );   //the truth-bank convention
  const double win_mc = ceelo::kDefaultFepWindowKeV;                      //what the MC actually scores

  const auto mu_compton_of = []( const Material &mat, const double e ) -> double {
    double mu = 0.0;
    for( const Material::ElementFractionPair &p : mat.elements )
      mu += p.second * mat.density
            * MassAttenuation::massAttenuationCoefficientElement( p.first->atomicNumber,
                  static_cast<float>(e), MassAttenuation::GammaEmProcces::ComptonScatter );
    for( const Material::NuclideFractionPair &p : mat.nuclides )
      mu += p.second * mat.density
            * MassAttenuation::massAttenuationCoefficientElement( p.first->atomicNumber,
                  static_cast<float>(e), MassAttenuation::GammaEmProcces::ComptonScatter );
    return mu;
  };

  map<long long,shared_ptr<const ceelo::DetectorResponse>> anchored;
  const auto centre_anchor = [&]( const double d_cm ) -> shared_ptr<const ceelo::DetectorResponse>
  {
    const long long key = llround( d_cm*1.0e6 );
    const map<long long,shared_ptr<const ceelo::DetectorResponse>>::const_iterator it = anchored.find( key );
    if( it != end(anchored) )
      return it->second;
    CeeLoUtils::TransferAnchor anchor;
    anchor.ref_distance_cm = d_cm;
    anchor.curve_derived = false;
    for( const double e : anchor_energies )
    {
      vector<unique_ptr<ceelo::Material>> owned;
      ceelo::EfficiencyCalculator pt;
      ceelo::ResponseGenerator::configure_calculator( pt, det.gd, owned );
      pt.set_point_source( Eigen::Vector3d( 0.0, 0.0, -(det.endcap_front_offset_cm + d_cm) ) );
      ostringstream k;
      k << "sph-anchor|point|d=" << setprecision(10) << d_cm;
      const McResult r = run_mc( pt, cache, k.str(), e, 0.001, 4000000000ULL, 20260904, 900.0 );
      anchor.curve.energies_keV.push_back( e );
      anchor.curve.eff.push_back( r.eff );
      anchor.curve.frac_sigma.push_back( r.frac_sigma );
    }
    ostringstream nm;
    nm << "sphere centre anchor (d=" << fixed << setprecision(2) << d_cm << " cm)";
    shared_ptr<const ceelo::DetectorResponse> resp
          = CeeLoUtils::makeTransferResponse( det.gd, anchor, ceelo::AnchorCurve{}, nm.str() );
    BOOST_REQUIRE( resp );
    anchored[key] = resp;
    return resp;
  };

  // One MC row: a sphere of radius r at centre distance d, of `mat` (null => transparent).
  const auto sphere_mc = [&]( const shared_ptr<const Material> &mat, const double r_cm,
                              const double d_cm ) -> McResult
  {
    ScenarioMcMaterials mats;
    ceelo::EfficiencyCalculator c;
    ceelo::ResponseGenerator::configure_calculator( c, det.gd, mats.owned );
    c.set_spherical_source( Eigen::Vector3d( 0.0, 0.0, -(det.endcap_front_offset_cm + d_cm) ),
                            r_cm, Eigen::Matrix3d::Identity(), 0.0 );
    ostringstream key;
    key << "sphere(r=" << setprecision(10) << r_cm << ",ri=0,d=" << d_cm
        << ",mat=" << (mat ? mat->name : string("TRANSPARENT")) << ")";
    if( mat )
    {
      mats.matrix.reset( new ceelo::Material( CeeLoUtils::to_ceelo_material( *mat ).to_material() ) );
      c.set_source_material( mats.matrix.get() );
    }
    return run_mc( c, cache, key.str(), energy, 0.0025, 16000000000ULL, 20260904, 1800.0 );
  };

  // The model at a given source removal coefficient.
  const auto model_eff = [&]( const shared_ptr<const ceelo::DetectorResponse> &resp,
                              const double r_cm, const double d_cm, const double mu_rem ) -> double
  {
    DistributedSrcCalcT<double> calc;
    calc.m_geometry = GeometryType::Spherical;
    calc.m_materialIndex = 0;
    calc.m_attenuateForAir = false;
    calc.m_isInSituExponential = false;
    calc.m_inSituRelaxationLength = -1.0;
    calc.m_srcVolumetricActivity = 1.0;
    calc.m_normalizeByVolume = false;
    calc.m_energy = energy;
    calc.m_effResponse = resp;
    calc.m_effMethod = ShieldingSourceFitCalc::VolumetricEffMethod::MCTransfer;
    calc.m_detector = detector_geom_from_config<double>( GeometryType::Spherical, d_cm*cm,
                                            det.gd.transverse_half_extent()*cm, 0.0 );
    DistributedSrcCalcT<double>::ShellInfo src;
    src.dims = { r_cm*cm, 0.0, 0.0 };
    src.trans_len_coef = mu_rem;
    src.type = ShellType::Material;
    calc.m_shells.push_back( src );
    integrate_on_path( calc, VolumetricIntegrator::Line, num_lines );
    return calc.integral / ((4.0/3.0)*PhysicalUnits::pi*r_cm*r_cm*r_cm*cm*cm*cm);
  };

  {
    const ceelo::Material cw = CeeLoUtils::to_ceelo_material( *water ).to_material();
    ostringstream o;
    o << "  60 keV.  window conventions: FWHM/2 = " << fixed << setprecision(3) << win_fwhm_half
      << " keV -> f_win(water) = " << setprecision(4)
      << ceelo::kn_in_window_fraction( energy, win_fwhm_half, cw )
      << ";  the MC's own scoring window = " << setprecision(3) << win_mc << " keV -> f_win(water) = "
      << setprecision(4) << ceelo::kn_in_window_fraction( energy, win_mc, cw );
    BOOST_TEST_MESSAGE( o.str() );
  }
  BOOST_TEST_MESSAGE( "  probe                      tau_d  transparent |  m/MC-1 @FWHM/2  @MC window |  f demanded" );

  struct Probe { const char *name; shared_ptr<const Material> mat; double r_cm; double centre_cm; };
  const vector<Probe> probes = {
    { "water R=3, d=4 cm",   water, 3.0,  4.0 },
    { "water R=3, d=6 cm",   water, 3.0,  6.0 },
    { "water R=3, d=50 cm",  water, 3.0, 50.0 },
    { "water R=1, d=6 cm",   water, 1.0,  6.0 },
    { "water R=5, d=6 cm",   water, 5.0,  6.0 },
    { "steel R=3, d=4 cm",   steel, 3.0,  4.0 },
  };

  vector<double> demanded_water;
  double worst_at_mc_window = 0.0;
  string worst_where;

  for( const Probe &p : probes )
  {
    const shared_ptr<const ceelo::DetectorResponse> resp = centre_anchor( p.centre_cm );
    const double mu_tot = transmition_length_coefficient( p.mat.get(), static_cast<float>(energy) );
    const double mu_c = mu_compton_of( *p.mat, energy );
    const ceelo::Material cmat = CeeLoUtils::to_ceelo_material( *p.mat ).to_material();
    const double f_fwhm = ceelo::kn_in_window_fraction( energy, win_fwhm_half, cmat );
    const double f_mcwin = ceelo::kn_in_window_fraction( energy, win_mc, cmat );

    const McResult mc = sphere_mc( p.mat, p.r_cm, p.centre_cm );
    const McResult mc_tr = sphere_mc( nullptr, p.r_cm, p.centre_cm );
    BOOST_REQUIRE( (mc.eff > 0.0) && (mc_tr.eff > 0.0) );

    // The geometry+response floor at THIS geometry, from the transparent twin.
    const double floor_ratio = model_eff( resp, p.r_cm, p.centre_cm, 0.0 ) / mc_tr.eff;

    const double at_fwhm = model_eff( resp, p.r_cm, p.centre_cm, mu_tot - f_fwhm*mu_c ) / mc.eff;
    const double at_mcwin = model_eff( resp, p.r_cm, p.centre_cm, mu_tot - f_mcwin*mu_c ) / mc.eff;

    // Solve for the f that closes the gap AFTER removing the floor: the attenuation model alone.
    double f_star = std::numeric_limits<double>::quiet_NaN();
    {
      double lo = 0.0, hi = 2.0;
      const auto resid = [&]( const double f ){
        return model_eff( resp, p.r_cm, p.centre_cm, mu_tot - f*mu_c )/mc.eff - floor_ratio;
      };
      if( (resid(lo)*resid(hi)) < 0.0 )
      {
        for( int it = 0; it < 40; ++it )
        {
          const double mid = 0.5*(lo + hi);
          if( (resid(lo)*resid(mid)) <= 0.0 ) hi = mid; else lo = mid;
        }
        f_star = 0.5*(lo + hi);
      }
    }

    ostringstream o;
    o << "  " << left << setw(24) << p.name << right << fixed << setprecision(2) << setw(7)
      << (mu_tot*2.0*p.r_cm*cm) << "  " << showpos << setprecision(2) << setw(8)
      << 100.0*(floor_ratio - 1.0) << "%" << " | " << setw(11) << 100.0*(at_fwhm - 1.0) << "%"
      << setw(11) << 100.0*(at_mcwin - 1.0) << "%" << noshowpos << " | ";
    if( std::isfinite(f_star) ) o << setprecision(4) << setw(8) << f_star; else o << "     n/a";
    o << "   (f_win " << setprecision(4) << f_fwhm << " -> " << f_mcwin << ")";
    BOOST_TEST_MESSAGE( o.str() );

    if( std::isfinite(f_star) && (p.mat == water) )
      demanded_water.push_back( f_star );

    // Prediction 3: at the MC's own window each row should sit at its transparent floor.  Judged on
    //  the rows that can actually see the credit - a steel sphere at tau = 50 is opaque, its integral
    //  is a surface skin, and changing mu barely moves it (which is itself the material lever's
    //  prediction, and is why its "demanded f" below is meaningless rather than discrepant).
    const double left_over = 100.0*(at_mcwin - floor_ratio);
    if( (p.mat == water) && (fabs(left_over) > worst_at_mc_window) )
    {
      worst_at_mc_window = fabs(left_over);
      worst_where = p.name;
    }
  }//for( probes )

  double f_min = 1.0e9, f_max = -1.0e9, f_sum = 0.0;
  for( const double f : demanded_water ){ f_min = std::min(f_min,f); f_max = std::max(f_max,f); f_sum += f; }
  const double f_mean = f_sum/std::max<size_t>(1,demanded_water.size());
  const ceelo::Material cw = CeeLoUtils::to_ceelo_material( *water ).to_material();
  const double f_mcwin_water = ceelo::kn_in_window_fraction( energy, win_mc, cw );

  BOOST_TEST_MESSAGE( "  water rows demand f = " << fixed << setprecision(4) << f_min << " .. " << f_max
                      << " (mean " << f_mean << ") over " << demanded_water.size()
                      << " configurations; f_win at the MC's own window = " << f_mcwin_water );
  BOOST_TEST_MESSAGE( "  after switching to the MC's window, the worst SENSITIVE row sits "
                      << setprecision(2) << worst_at_mc_window << "% from its transparent floor ("
                      << worst_where << "); steel is opaque at this energy and is reported as the"
                      " insensitivity check, not as a discrepancy" );

  BOOST_CHECK_MESSAGE( f_max <= 1.0,
    "a credit fraction above 1 is unphysical - no FEP-window credit can close this gap (max "
    << f_max << ")" );
  BOOST_CHECK_MESSAGE( (f_max - f_min) < 0.5*f_mean,
    "the demanded credit fraction is not consistent across geometry/optical depth (" << f_min
    << " .. " << f_max << "), so the residual is not the window credit" );
  BOOST_CHECK_MESSAGE( fabs(f_mean/f_mcwin_water - 1.0) < 0.35,
    "the demanded credit (" << f_mean << ") does not match the fraction the MC's own 0.75 keV"
    " scoring window implies (" << f_mcwin_water << ")" );
  BOOST_CHECK_MESSAGE( worst_at_mc_window < 0.75,
    "using the MC's own window does not collapse the residual onto the transparent floor; worst "
    << worst_at_mc_window << "% at " << worst_where );
}//BOOST_AUTO_TEST_CASE( Rung9_FepWindowCreditProbe )
