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

#include <array>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>

#ifdef _WIN32
#include "winsock2.h"
#include "Windows.h"
#endif

#define BOOST_TEST_MODULE ShieldingDimLimit_suite
#include <boost/test/included/unit_test.hpp>

#include "ceres/jet.h"

#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/GammaInteractionCalc_imp.hpp"

using namespace std;
using namespace boost::unit_test;

namespace
{
using GammaInteractionCalc::GeometryType;
using GammaInteractionCalc::DetectorGeomT;
using GammaInteractionCalc::detector_geom_from_config;
using GammaInteractionCalc::center_ray_exit_distance;

using Jet1 = ceres::Jet<double,1>;

const double cm = PhysicalUnits::cm;

// The detector face radius / setback are irrelevant to the geometric chord, but
//  must be something sane.
const double det_radius = 2.54 * cm;
const double det_setback = 0.0;

/** The point-source center-ray chord for a cumulative box/cylinder/sphere of the
 given outer dimensions, evaluated with plain doubles. */
double exit_double( GeometryType g, const std::array<double,3> &dims, double distance )
{
  const DetectorGeomT<double> det
        = detector_geom_from_config<double>( g, distance, det_radius, det_setback );
  return center_ray_exit_distance<double>( g, dims, det );
}

/** Same chord, with `swept` dimension carrying a unit derivative seed (lane 0), so
 the returned Jet's .v[0] is d(chord)/d(that dimension). */
Jet1 exit_jet( GeometryType g, const std::array<double,3> &dimvals, int swept, double distance )
{
  const DetectorGeomT<Jet1> det
        = detector_geom_from_config<Jet1>( g, Jet1(distance), det_radius, det_setback );
  std::array<Jet1,3> dims = { Jet1(dimvals[0]), Jet1(dimvals[1]), Jet1(dimvals[2]) };
  dims[swept].v[0] = 1.0;
  return center_ray_exit_distance<Jet1>( g, dims, det );
}

bool is_finite( double v ){ return std::isfinite(v); }

/** One chord-limit scenario: as the `swept` dimension of a single shell shrinks to
 zero (the other two held fixed), assert the chord value and its autodiff derivative
 stay finite, continuous, and converge to the documented zero-thickness limit.

 `expect_value(t)` is the analytic chord as a function of the swept dimension t.
 `expect_dval`    is the analytic d(chord)/d(swept) (constant for these cases: 1 for
                  the along-ray dimension, 0 for a perpendicular one). */
struct ChordCase
{
  string name;
  GeometryType geom;
  int swept;                       // index 0/1/2 of the dimension being shrunk
  std::array<double,3> fixed;      // values of all dims; fixed[swept] is overwritten by the sweep
  double distance;
  std::function<double(double)> expect_value;
  double expect_dval;
};

void run_chord_case( const ChordCase &c )
{
  // A geometric sweep of the swept dimension toward 0, then exactly 0, then a
  //  (non-physical) tiny negative to prove the branch is smooth across the bound.
  const vector<double> sweep = { 5.0*cm, 1.0*cm, 1.0e-2*cm, 1.0e-4*cm, 1.0e-6*cm,
                                 1.0e-9*cm, 1.0e-12*cm, 0.0, -1.0e-9*cm };

  double prev_deriv = std::numeric_limits<double>::quiet_NaN();
  for( const double t : sweep )
  {
    std::array<double,3> dims = c.fixed;
    dims[c.swept] = t;

    const Jet1 jet = exit_jet( c.geom, dims, c.swept, c.distance );
    const double val = jet.a;
    const double deriv = jet.v[0];

    const double exp_val = c.expect_value( t );

    BOOST_CHECK_MESSAGE( is_finite(val) && is_finite(deriv),
      c.name << " @ t=" << t/cm << "cm: non-finite chord (val=" << val << ", d=" << deriv << ")" );

    // Value: matches the analytic chord (absolute tol scaled to cm; chords are O(cm)).
    BOOST_CHECK_MESSAGE( fabs(val - exp_val) <= 1.0e-9*cm + 1.0e-7*fabs(exp_val),
      c.name << " @ t=" << t/cm << "cm: chord value " << val/cm << "cm != expected " << exp_val/cm << "cm" );

    // Derivative: matches the analytic limit derivative at every thickness, including 0.
    BOOST_CHECK_MESSAGE( fabs(deriv - c.expect_dval) <= 1.0e-6,
      c.name << " @ t=" << t/cm << "cm: d(chord)/d(dim) " << deriv << " != expected " << c.expect_dval );

    // Continuity: the derivative at successive (shrinking) thicknesses agrees - the
    //  "derivative at epsilon ~= derivative at zero" property the fit relies on.
    if( is_finite(prev_deriv) )
      BOOST_CHECK_MESSAGE( fabs(deriv - prev_deriv) <= 1.0e-6,
        c.name << " @ t=" << t/cm << "cm: derivative jumped " << prev_deriv << " -> " << deriv );
    prev_deriv = deriv;

    // Cross-check the autodiff derivative against a central finite difference of the
    //  double path, for strictly-positive (and not too tiny) thicknesses.
    if( t > 1.0e-5*cm )
    {
      const double h = 1.0e-3 * t;
      std::array<double,3> dp = dims, dm = dims;
      dp[c.swept] = t + h;
      dm[c.swept] = t - h;
      const double fd = ( exit_double(c.geom,dp,c.distance) - exit_double(c.geom,dm,c.distance) ) / (2.0*h);
      BOOST_CHECK_MESSAGE( fabs(deriv - fd) <= 1.0e-4 + 1.0e-4*fabs(fd),
        c.name << " @ t=" << t/cm << "cm: Jet deriv " << deriv << " != finite-diff " << fd );
    }
  }//for( sweep )
}//run_chord_case


using GammaInteractionCalc::DistributedSrcCalcT;
using GammaInteractionCalc::self_shielding_integration_imp;

/** Builds a single-material-shell volumetric calculator (the source is the only, innermost
 shell), with the `swept` outer dimension carrying the lane-0 derivative seed.  A fixed,
 finite attenuation coefficient is enough - the geometry, not the material, is under test. */
DistributedSrcCalcT<Jet1> make_vol_calc( GeometryType g, const std::array<double,3> &dims,
                                         int swept, double distance, bool exponential, double L )
{
  DistributedSrcCalcT<Jet1> calc;
  calc.m_geometry = g;
  calc.m_materialIndex = 0;
  calc.m_detector = detector_geom_from_config<Jet1>( g, Jet1(distance), det_radius, det_setback );
  calc.m_attenuateForAir = false;
  calc.m_airTransLenCoef = 0.0;
  calc.m_isInSituExponential = exponential;
  calc.m_inSituRelaxationLength = exponential ? L : -1.0;

  DistributedSrcCalcT<Jet1>::ShellInfo info;
  for( int i = 0; i < 3; ++i ) info.dims[i] = Jet1( dims[i] );
  info.dims[swept].v[0] = 1.0;
  info.trans_len_coef = Jet1( 0.10 / cm );   // a sane, finite per-length attenuation
  info.type = DistributedSrcCalcT<Jet1>::ShellType::Material;
  calc.m_shells.push_back( info );

  return calc;
}

/** Builds a nested two-shell calculator: a fixed non-source inner shell (index 0) surrounded by the
 source shell (index 1, m_materialIndex=1).  The source shell's `swept` cumulative *outer* dimension
 carries the lane-0 derivative seed; sweeping its thickness toward 0 collapses the source shell onto
 the inner core (its annular volume -> 0). */
DistributedSrcCalcT<Jet1> make_nested_vol_calc( GeometryType g, const std::array<double,3> &inner_dims,
                                                const std::array<double,3> &outer_dims, int swept,
                                                double distance )
{
  DistributedSrcCalcT<Jet1> calc;
  calc.m_geometry = g;
  calc.m_materialIndex = 1;   // source is the OUTER shell
  calc.m_detector = detector_geom_from_config<Jet1>( g, Jet1(distance), det_radius, det_setback );
  calc.m_attenuateForAir = false;
  calc.m_airTransLenCoef = 0.0;
  calc.m_isInSituExponential = false;
  calc.m_inSituRelaxationLength = -1.0;

  DistributedSrcCalcT<Jet1>::ShellInfo inner;
  for( int i = 0; i < 3; ++i ) inner.dims[i] = Jet1( inner_dims[i] );
  inner.trans_len_coef = Jet1( 0.05 / cm );
  inner.type = DistributedSrcCalcT<Jet1>::ShellType::Material;
  calc.m_shells.push_back( inner );

  DistributedSrcCalcT<Jet1>::ShellInfo outer;
  for( int i = 0; i < 3; ++i ) outer.dims[i] = Jet1( outer_dims[i] );
  outer.dims[swept].v[0] = 1.0;   // seed d/d(source-shell outer dim)
  outer.trans_len_coef = Jet1( 0.10 / cm );
  outer.type = DistributedSrcCalcT<Jet1>::ShellType::Material;
  calc.m_shells.push_back( outer );

  return calc;
}

/** Closed-form single-shell (origin) volume, as a Jet (carries the seeded dim derivative).
 Used only by the variant 0/1 reference normalization. */
[[maybe_unused]] Jet1 single_shell_vol_jet( GeometryType g, const std::array<Jet1,3> &o )
{
  const double pi = PhysicalUnits::pi;
  switch( g )
  {
    case GeometryType::Spherical:      return (4.0/3.0)*pi*o[0]*o[0]*o[0];
    case GeometryType::CylinderEndOn:
    case GeometryType::CylinderSideOn: return pi*o[0]*o[0]*2.0*o[1];
    case GeometryType::Rectangular:    return 8.0*o[0]*o[1]*o[2];
    default: return Jet1(0.0);
  }
}

/** Closed-form single-shell in-situ exponential depth-integral norm, as a Jet (mirrors the
 build switch).  Used only by the variant 0/1 reference normalization. */
[[maybe_unused]] Jet1 single_shell_exp_norm_jet( GeometryType g, const std::array<Jet1,3> &o, double L )
{
  using ceres::exp;
  const double pi = PhysicalUnits::pi;
  switch( g )
  {
    case GeometryType::Spherical:
    { const Jet1 R = o[0];        return 4.0*pi*L*(L*L*(2.0-2.0*exp(-R/L)) - 2.0*L*R + R*R); }
    case GeometryType::CylinderEndOn:
    { const Jet1 R = 2.0*o[1];    return L*(1.0 - exp(-R/L)); }
    case GeometryType::Rectangular:
    { const Jet1 R = 2.0*o[2];    return L*(1.0 - exp(-R/L)); }
    case GeometryType::CylinderSideOn:
    { const Jet1 R = o[0];        return 2.0*L*pi*(L*(exp(-R/L)-1.0)+R); }
    default: return Jet1(0.0);
  }
}

struct VolCase
{
  string name;
  GeometryType geom;
  int swept;                    // source dimension shrunk toward 0
  std::array<double,3> fixed;   // values of all dims; fixed[swept] overwritten by the sweep
  bool exponential;             // false = TotalActivity, true = in-situ ExponentialDistribution
};

/** Sweeps one source dimension of a single-shell volumetric source toward (and through) zero,
 building the contribution the way the compiled VOL_CALC_VARIANT does, and - at variant 2 -
 asserting the contribution value and its autodiff derivative stay finite and continuous
 through zero.  (At variants 0/1 the near-zero blow-up / floor is only printed, not asserted,
 so A/B builds stay green.) */
void run_vol_case( const VolCase &c )
{
  const double dist = 100.0*cm;
  const double L = 1.0*cm;       // relaxation length (exponential only)
  const double act = 3.7e10;     // arbitrary positive activity
  const double m2 = PhysicalUnits::m2;
  // Positive thicknesses shrinking to (and reaching) zero - the bounded fit's operating regime.
  //  (A negative thickness is non-physical and out of the source ray-tracing's domain, so it is
  //   not swept here; the bound at 0 is the relevant edge.)
  const vector<double> sweep = { 1.0*cm, 1.0e-2*cm, 1.0e-4*cm, 1.0e-6*cm,
                                 1.0e-9*cm, 1.0e-12*cm, 0.0 };

  vector<double> vals, derivs;
  for( const double t : sweep )
  {
    std::array<double,3> dims = c.fixed;
    dims[c.swept] = t;

    DistributedSrcCalcT<Jet1> calc = make_vol_calc( c.geom, dims, c.swept, dist, c.exponential, L );
    const std::array<Jet1,3> outer = calc.m_shells[0].dims;
    const double numerator = c.exponential ? (act / m2) : act;

#if( VOL_CALC_VARIANT == 2 )
    calc.m_srcVolumetricActivity = Jet1( numerator );
    calc.m_normalizeByVolume = true;
#elif( VOL_CALC_VARIANT == 0 )
    const Jet1 N0 = c.exponential ? single_shell_exp_norm_jet(c.geom,outer,L)
                                  : single_shell_vol_jet(c.geom,outer);
    calc.m_srcVolumetricActivity = Jet1(numerator) / N0;   // unguarded; -> Inf/NaN as shell -> 0
#else
    #error "Unknown VOL_CALC_VARIANT"
#endif

    self_shielding_integration_imp( calc, 1.0e-6, 2000000 );
    const Jet1 contrib = calc.integral * calc.m_srcVolumetricActivity;

    vals.push_back( contrib.a );
    derivs.push_back( contrib.v[0] );
    BOOST_TEST_MESSAGE( "    " << c.name << " @ t=" << t/cm << "cm: contrib=" << contrib.a
                        << " d/dt=" << contrib.v[0] );
  }//for( sweep )

#if( VOL_CALC_VARIANT == 2 )
  // (a) Finite VALUE at every thickness, including exactly 0.  This is the core win of variant 2:
  //     variant 0's act/vol (or act/m^2/norm) -> Inf/NaN as the source dimension -> 0, while
  //     variant 2 folds the normalization into the integrand and stays finite throughout.
  for( size_t i = 0; i < sweep.size(); ++i )
    BOOST_CHECK_MESSAGE( std::isfinite(vals[i]),
      c.name << " @ t=" << sweep[i]/cm << "cm: non-finite contribution value " << vals[i] );

  // (b) Finite, well-conditioned DERIVATIVE at every strictly-positive thickness (down to 1e-12 cm)
  //     - the finite Jacobian Ceres needs as a bounded dimension approaches its 0 lower bound.
  //     NOTE: at *exactly* 0 the off-center sphere/cylinder source ray-trace
  //     (exit_point_of_sphere_z_imp / cylinder_line_intersection_imp) has a pre-existing
  //     sqrt(arg)-with-arg->0 Jet singularity (and a value branch-switch for a zero-half-length
  //     cylinder) that is independent of this normalization reformulation - and outside the
  //     bounded fit's operating domain - so the exact-zero derivative is intentionally not asserted
  //     (see TODO.md).  rect-along-depth, whose ray-trace is the Stage-1 rectangle_exit_location,
  //     is in fact continuous through 0.
  const size_t last_pos = sweep.size() - 2;   // smallest strictly-positive thickness (1e-12 cm)
  for( size_t i = 0; i <= last_pos; ++i )
    BOOST_CHECK_MESSAGE( std::isfinite(derivs[i]),
      c.name << " @ t=" << sweep[i]/cm << "cm: non-finite contribution derivative " << derivs[i] );

  // (c) The value and derivative settle to a finite zero-thickness limit: the smallest positive
  //     thicknesses agree - the "derivative at epsilon ~= derivative at 0" property the fit relies
  //     on (here epsilon down to 1e-12 cm stands in for the 0 the ray-trace can't differentiate at).
  const double vref = vals[last_pos], dref = derivs[last_pos];
  for( size_t i = (last_pos >= 2 ? last_pos - 2 : 0); i <= last_pos; ++i )
  {
    BOOST_CHECK_MESSAGE( fabs(vals[i]-vref) <= 1.0e-4*fabs(vref) + 1.0e-9,
      c.name << ": contribution value not settling toward zero thickness: " << vals[i] << " vs " << vref );
    BOOST_CHECK_MESSAGE( fabs(derivs[i]-dref) <= 1.0e-3*fabs(dref) + 1.0e-3,
      c.name << ": contribution derivative not settling toward zero thickness: " << derivs[i] << " vs " << dref );
  }

  // (d) Geometries whose off-center source ray-trace is derivative-continuous through *exactly* zero
  //     thickness additionally require the value and derivative to be finite and continuous AT 0 (no
  //     NaN, no value jump).  Rectangular (Stage-1 rectangle_exit_location_imp) qualifies; sphere and
  //     cylinder are added here as their exit ray-traces (exit_point_of_sphere_z_imp /
  //     cylinder_line_intersection_imp) are reformulated to factor the vanishing radius out of the
  //     intersection discriminant (so sqrt() is never of a quadratically-vanishing argument).
  const bool smooth_through_zero =
        (c.geom == GeometryType::Rectangular)
      || (c.geom == GeometryType::Spherical)
      || (c.geom == GeometryType::CylinderEndOn)
      || (c.geom == GeometryType::CylinderSideOn);
  if( smooth_through_zero )
  {
    const size_t i0 = sweep.size() - 1;   // exactly 0
    BOOST_CHECK_MESSAGE( std::isfinite(derivs[i0]),
      c.name << ": non-finite contribution derivative at exactly 0: " << derivs[i0] );
    BOOST_CHECK_MESSAGE( fabs(vals[i0]-vref) <= 1.0e-4*fabs(vref) + 1.0e-9,
      c.name << ": contribution value discontinuous at exactly 0: " << vals[i0] << " vs " << vref );
    BOOST_CHECK_MESSAGE( fabs(derivs[i0]-dref) <= 1.0e-3*fabs(dref) + 1.0e-3,
      c.name << ": contribution derivative discontinuous at exactly 0: " << derivs[i0] << " vs " << dref );
  }
#else
  // Variants 0/1: the reformulation is inactive; only sanity-check the healthy regime so this
  //  case stays green when A/B-testing.  The near-zero behavior is printed above for the record.
  BOOST_CHECK_MESSAGE( std::isfinite(vals.front()) && std::isfinite(derivs.front()),
    c.name << ": non-finite contribution in the healthy regime (t=1cm)" );
#endif
}//run_vol_case

}//namespace


// The point-source center-ray chord through a single shell, as one of its
//  dimensions shrinks to zero.  On-axis, the chord equals the along-ray
//  dimension and is independent of the transverse dimensions; the derivative
//  w.r.t. the along-ray dimension is 1, and 0 w.r.t. a transverse one - and all
//  of this must hold continuously through (and at) zero thickness.
BOOST_AUTO_TEST_CASE( ChordZeroThicknessLimit_OnAxis )
{
  const double dist = 100.0*cm;
  vector<ChordCase> cases;

  // Spherical: chord = thickness (trivial baseline).
  cases.push_back( { "sphere-thickness", GeometryType::Spherical, 0, {1.0*cm,0,0}, dist,
                     [](double t){ return t; }, 1.0 } );

  // Rectangular: along-ray is the depth (index 2); width/height are transverse.
  cases.push_back( { "rect-depth(along)", GeometryType::Rectangular, 2, {5.0*cm,5.0*cm,1.0*cm}, dist,
                     [](double t){ return t; }, 1.0 } );
  cases.push_back( { "rect-width(perp)", GeometryType::Rectangular, 0, {1.0*cm,5.0*cm,5.0*cm}, dist,
                     [](double){ return 5.0*cm; }, 0.0 } );
  cases.push_back( { "rect-height(perp)", GeometryType::Rectangular, 1, {5.0*cm,1.0*cm,5.0*cm}, dist,
                     [](double){ return 5.0*cm; }, 0.0 } );

  // CylinderEndOn: along-ray is the half-length (index 1); radius is transverse.
  cases.push_back( { "cylEnd-halflen(along)", GeometryType::CylinderEndOn, 1, {5.0*cm,1.0*cm,0}, dist,
                     [](double t){ return t; }, 1.0 } );
  cases.push_back( { "cylEnd-radius(perp)", GeometryType::CylinderEndOn, 0, {1.0*cm,5.0*cm,0}, dist,
                     [](double){ return 5.0*cm; }, 0.0 } );

  // CylinderSideOn: along-ray is the radius (index 0); half-length is transverse.
  cases.push_back( { "cylSide-radius(along)", GeometryType::CylinderSideOn, 0, {1.0*cm,5.0*cm,0}, dist,
                     [](double t){ return t; }, 1.0 } );
  cases.push_back( { "cylSide-halflen(perp)", GeometryType::CylinderSideOn, 1, {5.0*cm,1.0*cm,0}, dist,
                     [](double){ return 5.0*cm; }, 0.0 } );

  for( const ChordCase &c : cases )
  {
    BOOST_TEST_MESSAGE( "  chord case: " << c.name );
    run_chord_case( c );
  }
}//ChordZeroThicknessLimit_OnAxis


// The volumetric-source contribution (integral * m_srcVolumetricActivity) as a single source
//  shell's dimension shrinks to zero.  At VOL_CALC_VARIANT 2 the reformulation keeps the value
//  and its autodiff derivative finite and continuous through zero (variant 0 -> Inf/NaN).  Covers
//  TotalActivity and in-situ ExponentialDistribution across all geometries.
BOOST_AUTO_TEST_CASE( VolumetricContributionZeroThicknessLimit )
{
  vector<VolCase> cases;

  // TotalActivity: any source dimension -> 0 collapses the volume (act/vol blows up at variant 0).
  cases.push_back( { "sph-TA-radius",      GeometryType::Spherical,     0, {1.0*cm,0,0},           false } );
  cases.push_back( { "cylEnd-TA-radius",   GeometryType::CylinderEndOn, 0, {1.0*cm,2.0*cm,0},      false } );
  cases.push_back( { "cylEnd-TA-halflen",  GeometryType::CylinderEndOn, 1, {2.0*cm,1.0*cm,0},      false } );
  cases.push_back( { "cylSide-TA-radius",  GeometryType::CylinderSideOn,0, {1.0*cm,2.0*cm,0},      false } );
  cases.push_back( { "cylSide-TA-halflen", GeometryType::CylinderSideOn,1, {2.0*cm,1.0*cm,0},      false } );
  cases.push_back( { "rect-TA-depth",      GeometryType::Rectangular,   2, {2.0*cm,2.0*cm,1.0*cm}, false } );
  cases.push_back( { "rect-TA-width",      GeometryType::Rectangular,   0, {1.0*cm,2.0*cm,2.0*cm}, false } );

  // In-situ ExponentialDistribution: the depth dimension -> 0 collapses the depth-integral norm.
  cases.push_back( { "sph-Exp-radius",     GeometryType::Spherical,     0, {1.0*cm,0,0},           true } );
  cases.push_back( { "cylEnd-Exp-halflen", GeometryType::CylinderEndOn, 1, {2.0*cm,1.0*cm,0},      true } );
  cases.push_back( { "cylSide-Exp-radius", GeometryType::CylinderSideOn,0, {1.0*cm,2.0*cm,0},      true } );
  cases.push_back( { "rect-Exp-depth",     GeometryType::Rectangular,   2, {2.0*cm,2.0*cm,1.0*cm}, true } );

  for( const VolCase &c : cases )
  {
    BOOST_TEST_MESSAGE( "  vol-contrib case: " << c.name );
    run_vol_case( c );
  }
}//VolumetricContributionZeroThicknessLimit


// Nested multi-shell probe: a non-source inner shell surrounded by the source shell, sweeping the
//  source shell's outer dimension so its (annular) volume -> 0.  This is the multi-shell
//  (m_materialIndex>0) path - NOT what the single-shell reformulation above targets (sphere uses the
//  cancelled r^2/(Ro^2+RoRi+Ri^2) weight; cylinder/rect use raw dV/vol_annular).  For POSITIVE
//  thickness the contribution value+derivative must stay finite - asserted, so this protects the
//  nested ray-trace/weight fixes against a regression to NaN.  The exactly-zero full-collapse point
//  is observe-only: it is a known degenerate corner (the sphere derivative is NaN there; cylinder/rect
//  collapse the value to 0), guarded in production by the trace-source dimension lower bound rather
//  than by the integrand (see the rect-innervoid follow-up).
BOOST_AUTO_TEST_CASE( NestedShellZeroThicknessProbe )
{
  const double dist = 100.0*cm;
  const double act = 3.7e10;

  struct NC { string name; GeometryType geom; std::array<double,3> inner; int swept; };
  const vector<NC> ncs = {
    { "sph-nested-radius",      GeometryType::Spherical,      {1.0*cm, 0,      0},      0 },
    { "cylSide-nested-radius",  GeometryType::CylinderSideOn, {1.0*cm, 1.0*cm, 0},      0 },
    { "rect-nested-width",      GeometryType::Rectangular,    {1.0*cm, 1.0*cm, 1.0*cm}, 0 },
  };

  const vector<double> sweep = { 1.0*cm, 1.0e-2*cm, 1.0e-4*cm, 1.0e-6*cm, 1.0e-9*cm, 0.0 };

  for( const NC &c : ncs )
  {
    BOOST_TEST_MESSAGE( "  nested-shell probe: " << c.name );
    for( const double t : sweep )
    {
      std::array<double,3> outer = c.inner;
      outer[c.swept] = c.inner[c.swept] + t;   // cumulative outer dim = inner + thickness t

      DistributedSrcCalcT<Jet1> calc = make_nested_vol_calc( c.geom, c.inner, outer, c.swept, dist );
#if( VOL_CALC_VARIANT == 2 )
      calc.m_srcVolumetricActivity = Jet1( act );
      calc.m_normalizeByVolume = true;
#else
      calc.m_srcVolumetricActivity = Jet1( 1.0 );
#endif
      try
      {
        // Modest eval cap: the inner-core silhouette near-discontinuity under-converges anyway (the
        //  point of the probe), and we only assert finiteness (robust to convergence), so a low cap
        //  keeps the probe fast.
        self_shielding_integration_imp( calc, 1.0e-6, 100000 );
      }catch( std::exception &e )
      {
        BOOST_TEST_MESSAGE( "    " << c.name << " @ t=" << t/cm << "cm: threw " << e.what() );
        BOOST_CHECK_MESSAGE( t <= 0.0, c.name << " @ t=" << t/cm << "cm: threw at positive thickness" );
        continue;
      }

      const Jet1 contrib = calc.integral * calc.m_srcVolumetricActivity;
      const bool ok = std::isfinite(contrib.a) && std::isfinite(contrib.v[0]);
      BOOST_TEST_MESSAGE( "    " << c.name << " @ t=" << t/cm << "cm: contrib=" << contrib.a
                          << " d/dt=" << contrib.v[0] << (ok ? "" : "   <-- NON-FINITE")
                          << "  (evals=" << calc.m_num_evals << ")" );
      if( t > 0.0 )
        BOOST_CHECK_MESSAGE( ok, c.name << " @ t=" << t/cm
                            << "cm: non-finite nested contribution at positive thickness" );
    }//for( sweep )
  }//for( nested cases )
}//NestedShellZeroThicknessProbe


// Configs where one thickness of the source shell is exactly 0 but the shell still has NON-zero
//  volume (so this is NOT the full-collapse case): a radial-only annulus (same length as the inner
//  cylinder), an endcaps-only shell (same radius), and a rectangular slab pair (zero thickness in
//  two of three dims).  Sub-domains with a zero extent are skipped, so the remaining region has a
//  real volume - value and derivative should be finite.  The NON-zero (would-be-fit) dimension is
//  seeded.
BOOST_AUTO_TEST_CASE( OneDimZeroThicknessProbe )
{
  const double dist = 100.0*cm, act = 3.7e10;

  struct OC {
    string name; GeometryType geom;
    std::array<double,3> inner, outer; int seed;   // seed = the non-zero (fit) thickness dimension
  };
  const vector<OC> cases = {
    // Same-length cylinder: outer adds radial thickness (R 1->1.5), endcap thickness exactly 0 (L=1=L).
    { "cyl-same-length(radial only)", GeometryType::CylinderSideOn, {1.0*cm,1.0*cm,0}, {1.5*cm,1.0*cm,0}, 0 },
    // Endcaps-only cylinder: radial thickness exactly 0 (R=1=1), outer adds length (L 1->1.5).
    { "cyl-endcaps-only",             GeometryType::CylinderSideOn, {1.0*cm,1.0*cm,0}, {1.0*cm,1.5*cm,0}, 1 },
    // Rect: width thickness 0.5, height & depth thickness exactly 0 (a +/- x slab pair).
    { "rect-width-slab-only",         GeometryType::Rectangular,    {1.0*cm,1.0*cm,1.0*cm}, {1.5*cm,1.0*cm,1.0*cm}, 0 },
  };

  for( const OC &c : cases )
  {
    DistributedSrcCalcT<Jet1> calc = make_nested_vol_calc( c.geom, c.inner, c.outer, c.seed, dist );
    calc.m_srcVolumetricActivity = Jet1( act );
    calc.m_normalizeByVolume = true;
    self_shielding_integration_imp( calc, 1.0e-6, 500000 );
    const Jet1 contrib = calc.integral * calc.m_srcVolumetricActivity;
    const bool ok = std::isfinite(contrib.a) && std::isfinite(contrib.v[0]);
    BOOST_TEST_MESSAGE( "  one-dim-zero: " << c.name << ": contrib=" << contrib.a
                        << " d/d(fit dim)=" << contrib.v[0] << (ok ? "" : "   <-- NON-FINITE")
                        << "  (evals=" << calc.m_num_evals << ")" );
    BOOST_CHECK_MESSAGE( ok, c.name << ": non-finite contribution (one thickness 0, volume non-zero)" );
  }
}//OneDimZeroThicknessProbe
