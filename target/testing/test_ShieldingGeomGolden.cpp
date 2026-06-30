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

/** Golden-value regression tests for the shielding geometry calculations
 (ray/volume intersections and the volumetric-source integrals), capturing both
 the scalar value AND the ceres::Jet derivative lanes in numerically healthy
 (non-degenerate) regimes.

 Purpose: the zero-thickness reformulation (Stage 2 - folding the 1/volume
 normalization into the integrand) is an exact algebraic identity away from the
 zero boundary, so it MUST reproduce these recorded values and derivatives bit-
 for-bit-to-tolerance.  These tests pin that contract down across a spread of
 dimensions and source types, so a refactor that is subtly wrong in the normal
 regime is caught immediately.

 To (re)record the goldens after an intended change: build with
 -DRECORD_GEOM_GOLDENS=1, run, and paste the printed map bodies directly into the
 g_chord / g_inter / g_integral initializers below (they are inlined here, not in
 separate .inc files).
 */

#include "InterSpec_config.h"

#include <map>
#include <array>
#include <cmath>
#include <string>
#include <vector>
#include <cstdio>
#include <iomanip>
#include <iostream>

#ifdef _WIN32
#include "winsock2.h"
#include "Windows.h"
#endif

#define BOOST_TEST_MODULE ShieldingGeomGolden_suite
#include <boost/test/included/unit_test.hpp>

#include "ceres/jet.h"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/GammaInteractionCalc_imp.hpp"

#include "VolumetricIntegrationFixtures.h"

#ifndef RECORD_GEOM_GOLDENS
#define RECORD_GEOM_GOLDENS 0
#endif

using namespace std;
using namespace boost::unit_test;

namespace
{
using GammaInteractionCalc::GeometryType;
using GammaInteractionCalc::DetectorGeomT;
using GammaInteractionCalc::detector_geom_from_config;
using GammaInteractionCalc::center_ray_exit_distance;
using GammaInteractionCalc::rectangle_exit_location_imp;
using GammaInteractionCalc::cylinder_line_intersection_imp;
using GammaInteractionCalc::exit_point_of_sphere_z_imp;
using GammaInteractionCalc::CylExitDir;

const double cm = PhysicalUnits::cm;
const double keV = PhysicalUnits::keV;

/** A recorded result: scalar value plus up to 3 derivative lanes. */
struct Golden
{
  double val;
  std::array<double,3> d;
  int nd;   // number of meaningful derivative lanes
};

// Relative+absolute tolerances: the reformulation is an exact identity in the
//  interior, so values/derivatives should match to a few ULP of the integrator
//  noise; we allow a small slack.
const double k_val_rel = 1.0e-7, k_val_abs = 1.0e-10;
const double k_drv_rel = 1.0e-6, k_drv_abs = 1.0e-10;

bool close( double a, double b, double rel, double abs_ )
{
  return std::fabs(a - b) <= (abs_ + rel*std::fabs(b));
}

#if RECORD_GEOM_GOLDENS
void emit( const string &name, const Golden &g )
{
  std::printf("  { \"%s\", { %.16g, { %.16g, %.16g, %.16g }, %d } },\n",
              name.c_str(), g.val, g.d[0], g.d[1], g.d[2], g.nd);
}
#endif

void check( const std::map<string,Golden> &goldens, const string &name, const Golden &got )
{
  BOOST_REQUIRE_MESSAGE( std::isfinite(got.val), name << ": non-finite value " << got.val );
  for( int i = 0; i < got.nd; ++i )
    BOOST_REQUIRE_MESSAGE( std::isfinite(got.d[i]), name << ": non-finite deriv[" << i << "] " << got.d[i] );

  const auto it = goldens.find( name );
  BOOST_REQUIRE_MESSAGE( it != goldens.end(), "No golden recorded for '" << name << "'" );
  const Golden &exp = it->second;

  BOOST_CHECK_MESSAGE( close(got.val, exp.val, k_val_rel, k_val_abs),
    name << ": value " << std::setprecision(12) << got.val << " != golden " << exp.val );
  BOOST_REQUIRE_EQUAL( got.nd, exp.nd );
  for( int i = 0; i < got.nd; ++i )
    BOOST_CHECK_MESSAGE( close(got.d[i], exp.d[i], k_drv_rel, k_drv_abs),
      name << ": deriv[" << i << "] " << std::setprecision(12) << got.d[i] << " != golden " << exp.d[i] );
}

// ----- value+derivative extraction from a ceres::Jet<double,N> -----
template<int N>
Golden from_jet( const ceres::Jet<double,N> &j, int nd )
{
  Golden g; g.val = j.a; g.nd = nd; g.d = {0,0,0};
  for( int i = 0; i < nd && i < N; ++i ) g.d[i] = j.v[i];
  return g;
}


// ============================ A) CHORD ============================
//  center_ray_exit_distance with all three outer dimensions seeded.
using Jet3 = ceres::Jet<double,3>;

Golden chord_golden( GeometryType g, const std::array<double,3> &dims,
                     double distance, double offx, double offy )
{
  DetectorGeomT<Jet3> det = detector_geom_from_config<Jet3>( g, Jet3(distance), 2.54*cm, 0.0, offx, offy );
  std::array<Jet3,3> d = { Jet3(dims[0],0), Jet3(dims[1],1), Jet3(dims[2],2) };
  return from_jet( center_ray_exit_distance<Jet3>(g,d,det), 3 );
}

struct ChordCase{ string name; GeometryType g; std::array<double,3> dims; double dist, offx, offy; };


// ====================== B) RAW INTERSECTIONS ======================
//  The ray/box, ray/cylinder, ray/sphere primitives the integrand uses, with an
//  interior source point and the dimensions seeded.
using Jet1 = ceres::Jet<double,1>;

Golden rect_inter_golden( double hw, double hh, double hd, int seed,
                          const std::array<double,3> &src, const std::array<double,3> &det )
{
  std::array<Jet1,3> half = { Jet1(hw), Jet1(hh), Jet1(hd) };
  half[seed].v[0] = 1.0;
  const Jet1 s[3] = { Jet1(src[0]), Jet1(src[1]), Jet1(src[2]) };
  const Jet1 dt[3] = { Jet1(det[0]), Jet1(det[1]), Jet1(det[2]) };
  Jet1 ep[3];
  return from_jet( rectangle_exit_location_imp<Jet1>(half[0],half[1],half[2],s,dt,ep), 1 );
}

Golden cyl_inter_golden( double rad, double hlen, int seed,
                         const std::array<double,3> &src, const std::array<double,3> &det )
{
  std::array<Jet1,2> dim = { Jet1(rad), Jet1(hlen) };
  dim[seed].v[0] = 1.0;
  const Jet1 s[3] = { Jet1(src[0]), Jet1(src[1]), Jet1(src[2]) };
  const Jet1 dt[3] = { Jet1(det[0]), Jet1(det[1]), Jet1(det[2]) };
  Jet1 ep[3];
  return from_jet( cylinder_line_intersection_imp<Jet1>(dim[0],dim[1],s,dt,CylExitDir::TowardDetector,ep), 1 );
}

Golden sphere_inter_golden( double rad, const std::array<double,3> &src, double obs )
{
  Jet1 r(rad); r.v[0] = 1.0;
  const Jet1 s[3] = { Jet1(src[0]), Jet1(src[1]), Jet1(src[2]) };
  Jet1 ep[3];
  return from_jet( exit_point_of_sphere_z_imp<Jet1>(s,ep,r,Jet1(obs)), 1 );
}


// ====================== C) VOLUME INTEGRALS =======================
//  self_shielding_integration_imp over a fixture, seeding the source-shell's
//  outer dimension `dim` (cumulative, so enclosing shells move with it).
Golden integral_golden( const VolumetricFixture::VolumetricSrcSpec &spec,
                        const MaterialDB &matdb, int dim, double epsrel, size_t max_evals )
{
  const GammaInteractionCalc::DistributedSrcCalcT<double> dbl
                              = VolumetricFixture::make_distributed_src_calc_t( spec, matdb );

  GammaInteractionCalc::DistributedSrcCalcT<Jet1> jet;
  jet.m_geometry = dbl.m_geometry;
  jet.m_materialIndex = dbl.m_materialIndex;
  jet.m_detector = detector_geom_from_config<Jet1>( spec.geometry, Jet1(spec.distance),
                          spec.detector_radius, spec.detector_setback, spec.offset_x, spec.offset_y );
  jet.m_attenuateForAir = dbl.m_attenuateForAir;
  jet.m_airTransLenCoef = dbl.m_airTransLenCoef;
  jet.m_isInSituExponential = dbl.m_isInSituExponential;
  jet.m_inSituRelaxationLength = dbl.m_inSituRelaxationLength;
  jet.m_srcVolumetricActivity = Jet1( dbl.m_srcVolumetricActivity );
  jet.m_energy = dbl.m_energy;

  for( size_t s = 0; s < dbl.m_shells.size(); ++s )
  {
    GammaInteractionCalc::DistributedSrcCalcT<Jet1>::ShellInfo info;
    for( size_t i = 0; i < 3; ++i ) info.dims[i] = Jet1( dbl.m_shells[s].dims[i] );
    info.trans_len_coef = Jet1( dbl.m_shells[s].trans_len_coef );
    info.type = dbl.m_shells[s].type;
    if( s >= spec.source_shell_index )
      info.dims[dim].v[0] = 1.0;   // cumulative outer dim
    jet.m_shells.push_back( info );
  }

  GammaInteractionCalc::self_shielding_integration_imp( jet, epsrel, max_evals );
  return from_jet( jet.integral, 1 );
}


// ============== D) TotalActivity CONTRIBUTION (act/vol)*I ==============
//  The production combination contrib = integral * m_srcVolumetricActivity, with
//  m_srcVolumetricActivity = act/vol the way TotalActivity trace sources set it.
//  This is the exact quantity Stage 2 reformulates (carry act, fold 1/vol into the
//  integrand weight); it MUST reproduce these golden value+derivative numbers.
//  Single source shell from the origin, so vol is the geometry's closed form
//  (computed inline as a Jet in the variant-0/1 branch below).
Golden total_activity_contrib_golden( const VolumetricFixture::VolumetricSrcSpec &spec,
                                      const MaterialDB &matdb, int dim, double act,
                                      double epsrel, size_t max_evals )
{
  // Build the integral as a Jet (dim seeded), exactly as integral_golden does.
  const GammaInteractionCalc::DistributedSrcCalcT<double> dbl
                              = VolumetricFixture::make_distributed_src_calc_t( spec, matdb );
  GammaInteractionCalc::DistributedSrcCalcT<Jet1> jet;
  jet.m_geometry = dbl.m_geometry;
  jet.m_materialIndex = dbl.m_materialIndex;
  jet.m_detector = detector_geom_from_config<Jet1>( spec.geometry, Jet1(spec.distance),
                          spec.detector_radius, spec.detector_setback, spec.offset_x, spec.offset_y );
  jet.m_attenuateForAir = dbl.m_attenuateForAir;
  jet.m_airTransLenCoef = dbl.m_airTransLenCoef;
  jet.m_isInSituExponential = dbl.m_isInSituExponential;
  jet.m_inSituRelaxationLength = dbl.m_inSituRelaxationLength;
  jet.m_energy = dbl.m_energy;

  std::array<Jet1,3> outer = { Jet1(0), Jet1(0), Jet1(0) };
  for( size_t s = 0; s < dbl.m_shells.size(); ++s )
  {
    GammaInteractionCalc::DistributedSrcCalcT<Jet1>::ShellInfo info;
    for( size_t i = 0; i < 3; ++i ) info.dims[i] = Jet1( dbl.m_shells[s].dims[i] );
    info.trans_len_coef = Jet1( dbl.m_shells[s].trans_len_coef );
    info.type = dbl.m_shells[s].type;
    if( s >= spec.source_shell_index )
      info.dims[dim].v[0] = 1.0;
    jet.m_shells.push_back( info );
    if( s == spec.source_shell_index )
      outer = info.dims;
  }
#if( VOL_CALC_VARIANT == 2 )
  // Production variant-2 path: carry act and let the integrand fold in 1/vol
  //  (m_normalizeByVolume); contrib = integral * act must reproduce the act/vol goldens.
  jet.m_srcVolumetricActivity = Jet1(act);
  jet.m_normalizeByVolume = true;
  (void)outer;

  GammaInteractionCalc::self_shielding_integration_imp( jet, epsrel, max_evals );

  const Jet1 contrib = jet.integral * jet.m_srcVolumetricActivity;
#else
  // Variants 0/1: integral is unnormalized (trans*dV); apply act/vol here, pinning the formula.
  jet.m_srcVolumetricActivity = Jet1(1.0);   // we apply act/vol ourselves below

  GammaInteractionCalc::self_shielding_integration_imp( jet, epsrel, max_evals );

  // vol as a Jet (single source shell from origin) carries the dim derivative.
  Jet1 vol;
  const double pi = PhysicalUnits::pi;
  switch( spec.geometry )
  {
    case GeometryType::Spherical:      vol = (4.0/3.0)*pi*outer[0]*outer[0]*outer[0]; break;
    case GeometryType::CylinderEndOn:
    case GeometryType::CylinderSideOn: vol = pi*outer[0]*outer[0]*2.0*outer[1]; break;
    case GeometryType::Rectangular:    vol = 8.0*outer[0]*outer[1]*outer[2]; break;
    default: break;
  }

  const Jet1 contrib = jet.integral * ( Jet1(act) / vol );
#endif

  return from_jet( contrib, 1 );
}


void set_data_dir()
{
  static bool s_have = false;
  if( s_have ) return;
  s_have = true;
  int argc = framework::master_test_suite().argc;
  char **argv = framework::master_test_suite().argv;
  string datadir;
  for( int i = 1; i < argc; ++i )
  {
    const string a = argv[i];
    if( SpecUtils::istarts_with(a,"--datadir=") )
      datadir = a.substr(10);
  }
  SpecUtils::ireplace_all( datadir, "%20", " " );
  if( datadir.empty() )
  {
    for( const auto &d : { "data", "../data", "../../data", "../../../data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d,"sandia.decay.xml") ) ){ datadir = d; break; }
    }
  }
  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db && db->nuclide("U238"), "Error initing SandiaDecayDataBase" );
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  BOOST_REQUIRE( MaterialDB::initialized() );
}

}//namespace


// ---- the recorded goldens (pasted directly below from a RECORD_GEOM_GOLDENS=1 run) ----
static const std::map<string,Golden> g_chord = {
  { "sph-r2", { 20, { 1, 0, 0 }, 3 } },
  { "sph-r0p5", { 5, { 1, 0, 0 }, 3 } },
  { "sph-r12", { 120, { 1, 0, 0 }, 3 } },
  { "rect-onaxis", { 30, { 0, 0, 1 }, 3 } },
  { "rect-wide", { 70, { 0, 0, 1 }, 3 } },
  { "rect-offaxis", { 30.0074990627343, { 0, 0, 1.00024996875781 }, 3 } },
  { "cylend-onaxis", { 80, { 0, 1, 0 }, 3 } },
  { "cylend-thin", { 120, { 0, 1, 0 }, 3 } },
  { "cylend-offaxis", { 80.00899949380694, { 0, 1.000112493672587, 0 }, 3 } },
  { "cylside-onaxis", { 60, { 1, 0, 0 }, 3 } },
  { "cylside-offaxis", { 60.01199880023994, { 1.000199980003999, 0, 0 }, 3 } },
};
static const std::map<string,Golden> g_inter = {
  { "rect-d-src", { 20.00127533940099, { 1.00006376697005, 0, 0 }, 1 } },
  { "rect-w-src", { 20.00127533940099, { 0, 0, 0 }, 1 } },
  { "rect-h-off", { 35.01590507529708, { 0, 0, 0 }, 1 } },
  { "cyl-r-src", { 70.00446368790348, { 0, 0, 0 }, 1 } },
  { "cyl-l-src", { 70.00446368790348, { 1.00006376697005, 0, 0 }, 1 } },
  { "cyl-r-side", { 60.00299992500375, { 1.000049998750062, 0, 0 }, 1 } },
  { "sph-r-src", { 19.18814717904376, { 1.029604974856527, 0, 0 }, 1 } },
  { "sph-r-deep", { 67.83144536704604, { 1.027663145877849, 0, 0 }, 1 } },
};
static const std::map<string,Golden> g_integral = {
  { "I:sph-U-self-atten-185keV:d0", { 0.07311230543633152, { 0.007313102406223047, 0, 0 }, 1 } },
  { "C:sph-U-self-atten-185keV:d0", { 2.181784652069365e-06, { -1.090333733189737e-07, 0, 0 }, 1 } },
  { "I:sph-U-self-atten-1001keV:d0", { 1.337520695702269, { 0.1418685626814017, 0, 0 }, 1 } },
  { "C:sph-U-self-atten-1001keV:d0", { 3.991369316413629e-05, { -1.753475908354042e-06, 0, 0 }, 1 } },
  { "I:sph-U-plus-Fe-661keV:d0", { 0.4158099388105139, { 0.0408480993052411, 0, 0 }, 1 } },
  { "I:cyl-end-Fe-661keV:d0", { 23.65254821587288, { 0.9442574301613307, 0, 0 }, 1 } },
  { "C:cyl-end-Fe-661keV:d0", { 3.011535972220448e-05, { -2.3484883967715e-09, 0, 0 }, 1 } },
  { "I:cyl-end-Fe-661keV:d1", { 23.65254821587288, { 0.05650853439299153, 0, 0 }, 1 } },
  { "I:cyl-end-U-plus-Fe-661keV:d1", { 0.5134100487969468, { 0.001064294353760296, 0, 0 }, 1 } },
  { "I:cyl-side-Fe-661keV:d0", { 55.7527771347177, { 1.24444762991284, 0, 0 }, 1 } },
  { "C:cyl-side-Fe-661keV:d0", { 3.549332028836448e-05, { -6.274928446560906e-07, 0, 0 }, 1 } },
  { "I:cyl-side-U-plus-Fe-661keV:d0", { 1.436796290616525, { 0.07138608012417172, 0, 0 }, 1 } },
  { "I:rect-Fe-661keV:d2", { 30.09581405107642, { 0.07185445176650262, 0, 0 }, 1 } },
  { "C:rect-Fe-661keV:d2", { 3.009581405107642e-05, { -5.300618292550259e-07, 0, 0 }, 1 } },
  { "I:rect-Fe-661keV:d0", { 30.09581405107642, { 0.6003532942618752, 0, 0 }, 1 } },
  { "I:rect-U-plus-Fe-661keV:d2", { 0.6536141650999343, { 0.001354767414810882, 0, 0 }, 1 } },
};


BOOST_AUTO_TEST_CASE( GeometryGoldens )
{
  set_data_dir();
  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );

  const double dist = 100.0*cm;

  // ---- A) Chord cases: a spread of dimensions, on- and off-axis ----
  const vector<ChordCase> chord = {
    { "sph-r2",        GeometryType::Spherical,    {2.0*cm,0,0},          dist, 0, 0 },
    { "sph-r0p5",      GeometryType::Spherical,    {0.5*cm,0,0},          dist, 0, 0 },
    { "sph-r12",       GeometryType::Spherical,    {12.0*cm,0,0},         dist, 0, 0 },
    { "rect-onaxis",   GeometryType::Rectangular,  {5.0*cm,4.0*cm,3.0*cm},dist, 0, 0 },
    { "rect-wide",     GeometryType::Rectangular,  {10.0*cm,1.0*cm,7.0*cm},dist,0, 0 },
    { "rect-offaxis",  GeometryType::Rectangular,  {5.0*cm,4.0*cm,3.0*cm},dist, 2.0*cm, 1.0*cm },
    { "cylend-onaxis", GeometryType::CylinderEndOn,{5.0*cm,8.0*cm,0},     dist, 0, 0 },
    { "cylend-thin",   GeometryType::CylinderEndOn,{1.0*cm,12.0*cm,0},    dist, 0, 0 },
    { "cylend-offaxis",GeometryType::CylinderEndOn,{5.0*cm,8.0*cm,0},     dist, 1.5*cm, 0 },
    { "cylside-onaxis",GeometryType::CylinderSideOn,{6.0*cm,10.0*cm,0},   dist, 0, 0 },
    { "cylside-offaxis",GeometryType::CylinderSideOn,{6.0*cm,10.0*cm,0},  dist, 0, 2.0*cm },
  };
#if RECORD_GEOM_GOLDENS
  std::printf("\n// ===== chord goldens =====\n");
#endif
  for( const ChordCase &c : chord )
  {
    const Golden g = chord_golden( c.g, c.dims, c.dist, c.offx, c.offy );
#if RECORD_GEOM_GOLDENS
    emit( c.name, g );
#else
    check( g_chord, c.name, g );
#endif
  }

  // ---- B) Raw intersection primitives, interior source points ----
#if RECORD_GEOM_GOLDENS
  std::printf("\n// ===== raw-intersection goldens =====\n");
#endif
  {
    const std::array<double,3> det_z = { 0, 0, dist };
    const std::array<double,3> det_x = { dist, 0, 0 };
    struct R{ string n; Golden g; };
    const vector<R> rs = {
      { "rect-d-src", rect_inter_golden(5*cm,4*cm,3*cm,2,{1*cm,0.5*cm,1*cm},det_z) },
      { "rect-w-src", rect_inter_golden(5*cm,4*cm,3*cm,0,{1*cm,0.5*cm,1*cm},det_z) },
      { "rect-h-off", rect_inter_golden(6*cm,5*cm,4*cm,1,{-1*cm,1*cm,0.5*cm},{2*cm,1*cm,dist}) },
      { "cyl-r-src",  cyl_inter_golden(5*cm,8*cm,0,{1*cm,0.5*cm,1*cm},det_z) },
      { "cyl-l-src",  cyl_inter_golden(5*cm,8*cm,1,{1*cm,0.5*cm,1*cm},det_z) },
      { "cyl-r-side", cyl_inter_golden(6*cm,10*cm,0,{0,0,1*cm},det_x) },
      { "sph-r-src",  sphere_inter_golden(3*cm,{0.5*cm,0.5*cm,1*cm},dist) },
      { "sph-r-deep", sphere_inter_golden(10*cm,{2*cm,-1*cm,3*cm},dist) },
    };
    for( const R &r : rs )
    {
#if RECORD_GEOM_GOLDENS
      emit( r.n, r.g );
#else
      check( g_inter, r.n, r.g );
#endif
    }
  }

  // ---- C) Volume integrals + D) TotalActivity contribution ----
  //  Only numerically-healthy (well-converged, single source-shell or simple
  //  point-source-plus-shielding) cases - the regression contract for Stage 2.
  const double epsrel = 1.0e-6;
  const size_t max_evals = 20000000;
  const vector<VolumetricFixture::VolumetricSrcSpec> fixtures
                                = VolumetricFixture::standard_volumetric_fixtures();
  struct IntCase{ string fixture; int dim; bool contrib; };
  const vector<IntCase> int_cases = {
    { "sph-U-self-atten-185keV", 0, true },
    { "sph-U-self-atten-1001keV", 0, true },
    { "sph-U-plus-Fe-661keV", 0, false },
    { "cyl-end-Fe-661keV", 0, true },
    { "cyl-end-Fe-661keV", 1, false },
    { "cyl-end-U-plus-Fe-661keV", 1, false },
    { "cyl-side-Fe-661keV", 0, true },
    { "cyl-side-U-plus-Fe-661keV", 0, false },
    { "rect-Fe-661keV", 2, true },
    { "rect-Fe-661keV", 0, false },
    { "rect-U-plus-Fe-661keV", 2, false },
  };
#if RECORD_GEOM_GOLDENS
  std::printf("\n// ===== integral (+contrib) goldens =====\n");
#endif
  for( const IntCase &ic : int_cases )
  {
    const VolumetricFixture::VolumetricSrcSpec *spec = nullptr;
    for( const auto &f : fixtures ) if( f.name == ic.fixture ) spec = &f;
    BOOST_REQUIRE_MESSAGE( spec, "missing fixture " << ic.fixture );

    const string key = ic.fixture + ":d" + std::to_string(ic.dim);
    const Golden gi = integral_golden( *spec, *matdb, ic.dim, epsrel, max_evals );
#if RECORD_GEOM_GOLDENS
    emit( "I:" + key, gi );
#else
    check( g_integral, "I:" + key, gi );
#endif

    // Contribution golden only for single-shell (origin) fixtures where vol is closed-form.
    if( ic.contrib )
    {
      const Golden gc = total_activity_contrib_golden( *spec, *matdb, ic.dim, 1.0,
                                                       epsrel, max_evals );
#if RECORD_GEOM_GOLDENS
      emit( "C:" + key, gc );
#else
      check( g_integral, "C:" + key, gc );
#endif
    }
  }
}//GeometryGoldens
