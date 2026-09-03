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
        const VolumetricIntegrator prev = sm_volumetric_integrator_override;
        sm_volumetric_integrator_override = VolumetricIntegrator::Line;
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
        sm_volumetric_integrator_override = prev;

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

  // The element path costs 20-40 s per energy on a box, so the whole matrix is an hours-long run.
  //  INTERSPEC_LINE_AB_SCENARIOS=name[,name...] restricts it while iterating; unset runs everything.
  std::set<string> only;
  if( const char *env = std::getenv("INTERSPEC_LINE_AB_SCENARIOS") )
  {
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
    BOOST_TEST_MESSAGE( "  (restricted to " << only.size() << " scenario(s) by INTERSPEC_LINE_AB_SCENARIOS)" );
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
  //  the line path within 0.21% of its 2^18-line one (LineCountConvergence), and the worst row here
  //  measured 0.44% (transparent box at 122 keV, 2026-09-03).  0.75% leaves about one part in three
  //  of headroom over that - enough for platform floating-point drift, not enough to hide a real
  //  disagreement, which for two quadratures of the same integral would show up as a systematic
  //  trend across a scenario family rather than one loose row.
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
