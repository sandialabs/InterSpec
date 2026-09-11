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
#ifndef VolumetricLineIntegration_imp_hpp
#define VolumetricLineIntegration_imp_hpp

/** Detector-side LINE integration of volumetric sources.

 This file is an implementation fragment of GammaInteractionCalc_imp.hpp - it is included from the
 end of that header (after `namespace GammaInteractionCalc` closes) and must not be included on its
 own.  It needs #DistributedSrcCalcT, the shell ray-trace helpers and #self_shielding_integration_imp
 from there, and is forward-declared before #ShieldingSourceChi2Fcn::expected_peak_counts_imp.

 THE IDEA.  The per-element extended-source kernel (eval_cylinder / eval_rect) integrates over
 SOURCE points and, at each, over a fan of ~500 rays to the crystal:

     eps = Int_V dV rho(r) P(r) (1/4pi) Int dOmega k(r,w) T(r,w)

 with P the response's prefactor (eta * near-field N * grounding k), k the per-ray crystal kernel,
 T the transmission through the source and shields.  Reversing the order of integration turns it
 into an integral over LINES that hit the active crystal, parameterised on the DETECTOR side by a
 hull point x and a direction w with the etendue measure dA |w.n| dOmega:

     eps = (1/4pi) Int_S dA Int dOmega |w.n| k(line) Int_chord ds rho(s) P(s) T(s)

 which is the same number (dV dOmega_from_r == dA |w.n| dOmega ds).  Everything on the detector
 side of a line - hull point, direction, the material segments through endcap/dead layer/crystal -
 depends on neither the source nor the energy, so a fixed line set is built ONCE per fit
 (#VolumetricLineCache) and reused by every energy and every cost-function evaluation; k(line;E)
 is memoized per energy.  Per evaluation the only new work is intersecting each line with the
 current (fit-parameter, T-valued) shells, and per energy one exponential per source piece: the
 chord integral of exp(-mu s) is analytic, and the smooth remainder (P, an in-situ profile) is a
 2-4 point Gauss-Legendre average in y = exp(-mu (s - s0)), s0 the near end of the piece.

 What this buys over the element path (measured there: 5.8e-4 s per element at 512 rays, boxes
 needing 11k-67k elements PER ENERGY, nothing shared across energies or evaluations): the cost is
 (lines) x (energies) exponentials, with the detector-side geometry of a line and its per-energy
 kernel hoisted out of the per-energy work, and the chords exact in T, so d(integral)/d(dims) does
 not carry the frozen-aperture staircase error the element path documents at eval_cylinder.  The
 lines themselves are re-traced through the crystal once per distinct set of scalar source
 dimensions (#VolumetricLineCache::traced), which is once or twice per optimizer step.

 DIRECTION PROPOSAL.  Directions are aimed at points OF THE SOURCE rather than into a cone around
 its bounding sphere: every line then crosses the source by construction, whatever its aspect ratio -
 where a cone around the bounding SPHERE would waste most of its lines on empty space for a needle
 or an edge-on sheet.  The aim points are frozen in NORMALISED coordinates (the unit solid, or a
 unit face) and scaled by the CURRENT source dimensions at every evaluation (#line_direction_imp),
 so the set deforms continuously with the source instead of the source sliding through a fixed set:
 chord, weight and direction are then smooth functions of the fitted dimensions, and the set never
 has to be re-drawn.  That is what a FIT needs - see below for why.

 The proposal is a MIXTURE of a volume component (a point uniform in the source padded by `pad`,
 direction density (s1^3 - s0^3)/(3 V_p) over the line's chord through it) and a surface component
 (a point uniform on the source's own boundary, sampled per face with fixed probabilities; density
 sum over the line's crossings P of p_face |P - x|^2 / (A_face |n.w|)).  #line_proposal_density_imp
 sums them.  A hollow source gets NO component on its inner boundary - see
 #VolumetricLineCache::frac_outer for the measurement that rules one out.

 WHY THE SURFACE COMPONENT, which is the whole reason the mixture exists.  Differentiating a volume
 integral whose boundary moves with the parameter gives, by the Reynolds transport theorem,
 d/dR Int_V f dV = Int_{dV} f (v.n) dA + Int_V d_R f dV.  A line estimator carries the transverse
 measure dA_perp dOmega and dA_perp = (n.w) dA_S, so a line's share of that BOUNDARY term comes with
 a 1/(n.w) at its exit point - which is exactly the R/sqrt(R^2 - b^2) that a chord's derivative has
 at the limb.  Under a volume-only proposal the per-line derivative therefore has DIVERGENT variance
 (integrable, so unbiased, but heavy-tailed): measured against a converged element-path finite
 difference the value was right to 0.1% at 65536 lines while the gradient wobbled +-5-14% and
 flipped sign between adjacent radii - an objective that is smooth in value but rough in slope,
 which is worse for Levenberg-Marquardt than a constant offset, since the roughness is largest
 exactly where the true gradient is smallest and it corrupts the curvature the dimension
 UNCERTAINTY is read from.  The surface component's density diverges at that same limb, so the
 mixture weight 1/p vanishes there and every line's contribution - value AND derivative - is
 bounded.  `sm_default_volumetric_line_surface_frac` is the mixture weight.

 A third, HEMISPHERE component (directions uniform above the hull normal, density 1/2pi, no
 dependence on the dimensions at all) is added for WIDE sources - a 20 m in-situ disk seen from 1 m
 - where the aim-point components' s^2 direction Jacobian leaves the weights spanning the source's
 (far/near)^2; it puts a floor under the mixture density and cut the replica scatter 7-16x there
 (#VolumetricLineCache::frac_hemi, #sm_volumetric_line_hemi_ratio).

 The mixture is unbiased whatever the fractions are (each component covers the source), so they are
 a variance knob only.  Aiming with the CURRENT dimensions is what makes the estimator's derivative
 the derivative of what it estimates: the crystal kernel k moves with the line, and its direction
 gradient is carried by a forward difference of two extra traces (#VolumetricLineCache::TracedLines,
 only when a dimension is being differentiated).  That term is NOT zero in the continuum limit - it
 is of order the relative variation of k across the set - which is why it is carried rather than
 dropped.

 LIMIT.  A source extent that is a vanishing fraction of the standoff makes chord/volume a 0/0 here.
 The intervals are computed from each line's closest approach to the assembly origin
 (#line_shell_intervals_imp), which keeps the chords accurate down to an extent/distance ratio of
 ~1e-10; below that the dispatcher (#integrate_volumetric_calculators) integrates a copy of the
 source whose vanishing extent is floored at that ratio (#sm_line_path_extent_ratio_floor), so the
 value and derivative stay finite and continuous all the way to exactly zero and the line path owns
 the whole domain (test_ShieldingDimLimit pins this).

 SEQUENCE AND ERROR.  The per-line unit coordinates come from a host-side stream
 (#LineSampleStream; Sobol' with a random digital shift by default, Halton kept for reproducing the
 original construction), so a REPLICA of the set - an independent randomisation - is one parameter
 (#LineSampleParams).  Replicas are how the quadrature's precision is measured; within one set the
 lines are summed in fixed contiguous blocks whose two-scale scatter gives `m_est_rel_error`,
 calibrated against replicas (LineErrorEstimateCalibration).  Against two independent per-voxel
 references (VolumetricReferenceIntegrator.h) the line path has no measurable bias
 (LinePathVsReference); its precision at the shipped line count is what the mixture fractions set.

 UNITS.  CeeLo works in cm; everything here is in PhysicalUnits.  A line's weight is its hull
 point's area share (cm^2) times the hull cosine, over the proposal density (per steradian) and
 4 pi and the line count, converted to PhysicalUnits; multiplied by an emission density (1/volume)
 and a chord (length) it is dimensionless.
 */

#include <map>
#include <array>
#include <cmath>
#include <deque>
#include <mutex>
#include <tuple>
#include <atomic>
#include <memory>
#include <random>
#include <vector>
#include <cassert>
#include <algorithm>
#include <type_traits>

#include <boost/random/sobol.hpp>

#include <Eigen/Geometry>

#include "io/DetectorEtendue.h"
#include "io/LowDiscrepancy.h"
#include "io/DetectorResponse.h"

#include "InterSpec/CeeLoUtils.h"

namespace GammaInteractionCalc
{

/** Which quadrature a volumetric calculator is integrated with. */
enum class VolumetricIntegrator : int
{
  /** Production choice: the line path whenever a response and a line set are attached (cascade
   summing, the effective-shielding report and collimated responses are all served by it), else -
   flat-disk, i.e. no response - the element quadrature underneath eval_*. */
  Auto,
  /** Force the per-element aperture path (the reference implementation). */
  Element,
  /** Force the line path (throws for calculators it cannot serve). */
  Line
};//enum class VolumetricIntegrator

/** TEST HOOK - overrides the integrator choice for every calculator in the process.  Leave at
 Auto in production; the tests use it to A/B the two quadratures on identical calculators. */
inline VolumetricIntegrator sm_volumetric_integrator_override = VolumetricIntegrator::Auto;

/** Sets #sm_volumetric_integrator_override for a scope and restores it on exit - a throw inside
 the scope must not leave the process forcing one path for every later test case. */
struct ScopedVolumetricIntegratorOverride
{
  const VolumetricIntegrator previous;
  explicit ScopedVolumetricIntegratorOverride( const VolumetricIntegrator path )
    : previous( sm_volumetric_integrator_override )
  {
    sm_volumetric_integrator_override = path;
  }
  ~ScopedVolumetricIntegratorOverride()
  {
    sm_volumetric_integrator_override = previous;
  }
  ScopedVolumetricIntegratorOverride( const ScopedVolumetricIntegratorOverride & ) = delete;
  ScopedVolumetricIntegratorOverride &operator=( const ScopedVolumetricIntegratorOverride & ) = delete;
};//struct ScopedVolumetricIntegratorOverride

/** Gauss-Legendre points on a source chord, in the y = exp(-mu_eff (s - s0)) substitution (s0 the
 near end of the piece; 2 = exact for a remainder linear in y; radial in-situ profiles use 4 - see
 #line_source_integration_imp). */
inline int sm_line_chord_gl_points = 2;

/** Smallest source extent the line path integrates, as a RATIO of the extent to the source-detector
 distance: a source shell thinner than this in any dimension is integrated as if it were exactly
 this thick (the SCALAR part of the dimension is floored, its derivative lane is kept - see
 #integrate_volumetric_calculators), so the value and derivative stay finite and continuous down to
 an extent of exactly zero.

 The ratio, not an absolute size, is what matters.  Computed from a detector-side origin, a line's
 chord through the source comes out of a ray/quadric intersection whose discriminant is a
 difference of terms of order (distance)^2, i.e. a relative error of eps*(distance/extent)^2
 (measured: the sphere sweep in test_ShieldingDimLimit swung 38% at extent/distance = 1e-8).
 #line_shell_intervals_imp therefore re-origins every line at its closest approach to the assembly
 origin first: the single subtraction that involves the detector-side origin then costs
 eps*(distance/extent) - about 1e-6 at this ratio - and nothing downstream sees a distance-sized
 term.  Below the floor the chord/volume quotient the integrand needs is a genuine 0/0 (both vanish
 together); its limit is what the floored source reproduces, to O(ratio). */
inline double sm_line_path_extent_ratio_floor = 1.0e-10;


/** Picks the argument with the larger / smaller SCALAR part (a kink, not a smoothing - the
 derivative follows whichever branch is active). */
template<typename T>
inline const T &select_max( const T &x, const T &y ) { return (scalar_of(x) >= scalar_of(y)) ? x : y; }
template<typename T>
inline const T &select_min( const T &x, const T &y ) { return (scalar_of(x) <= scalar_of(y)) ? x : y; }


/** Interval [a,b] of the line point(s) = o - s*d (s the distance BACK from `o` toward the source;
 d the unit photon direction) inside the slab |coord| <= half, on one axis.  False when the line
 misses the slab.  A ray parallel to the slab's planes is inside (unbounded interval, with large
 sentinels) or outside for every s.

 `D` is the direction's type: double for the cached detector-side lines, T for a ray whose direction
 depends on the fit parameters (#shell_path_to_point_imp). */
template<typename T, typename D>
inline bool line_slab_interval_imp( const T &o_c, const D &d_c, const T &half, T &a, T &b )
{
  static const double big = 1.0e300;
  if( scalar_of(d_c) == 0.0 )
  {
    if( std::fabs(scalar_of(o_c)) >= std::fabs(scalar_of(half)) )
      return false;
    a = T(-big);
    b = T(big);
    return true;
  }
  // coord(s) = o_c - s d_c in [-half, half]
  const T s1 = (o_c - half) / d_c;
  const T s2 = (o_c + half) / d_c;
  a = select_min( s1, s2 );
  b = select_max( s1, s2 );
  return scalar_of(a) < scalar_of(b);
}//line_slab_interval_imp(...)


/** Interval of the line inside the quadric |o_perp - s d_perp|^2 <= R^2 over the given axes
 (two axes: an infinite cylinder along the third; three: a sphere).  False when missed.

 Solved through the closest approach s* = (o.d)/(d.d) and the impact parameter |o_perp - s* d_perp|
 formed as a VECTOR, rather than through the textbook discriminant B^2 - A C: that discriminant is a
 difference of two terms of order |o|^2 whose result is of order R^2, so it loses
 eps*(|o|/R)^2 of relative precision, while the vector form's only cancellation is the linear
 o - s* d.  See #sm_line_path_extent_ratio_floor for why this matters. */
template<int NAxes, typename T, typename D>
inline bool line_quadric_interval_imp( const T o[3], const D d[3], const T &radius, T &a, T &b )
{
  using namespace std;
  using namespace ceres;
  static const double big = 1.0e300;

  // The direction may be scalar (cached lines) or T (a node ray); A follows it.
  using AType = std::conditional_t<std::is_same_v<D,double>,double,T>;
  AType A( 0.0 );
  T B( 0.0 );
  for( int i = 0; i < NAxes; ++i )
  {
    A += d[i]*d[i];
    B += o[i]*d[i];
  }

  if( scalar_of(A) < 1.0e-24 )
  {
    // Parallel to the cylinder axis: inside for every s, or never.
    T C( 0.0 );
    for( int i = 0; i < NAxes; ++i )
      C += o[i]*o[i];
    if( scalar_of(C) >= scalar_of(radius*radius) )
      return false;
    a = T(-big);
    b = T(big);
    return true;
  }

  // Closest approach, and the squared impact parameter from the vector o_perp - s* d_perp.
  const T s_star = B / A;
  T b2( 0.0 );
  for( int i = 0; i < NAxes; ++i )
  {
    const T q = o[i] - s_star*d[i];
    b2 += q*q;
  }

  const T disc = radius*radius - b2;
  if( scalar_of(disc) <= 0.0 )
    return false;
  const T root = sqrt( disc / A );
  a = s_star - root;
  b = s_star + root;
  return true;
}//line_quadric_interval_imp(...)


/** Interval of the line inside one shell of the given geometry (dims per #GeometryType:
 cylinders {radius, half_z}, box {hx, hy, hz}, sphere {radius}).  False when missed.
 The assembly frame is used throughout: cylinders along z, box axis-aligned, sphere at the origin. */
template<typename T, typename D>
inline bool line_shell_interval_imp( const GeometryType geometry, const std::array<T,3> &dims,
                                     const T o[3], const D d[3], T &a, T &b )
{
  switch( geometry )
  {
    case GeometryType::Spherical:
      return line_quadric_interval_imp<3>( o, d, dims[0], a, b )
             && (scalar_of(a) < scalar_of(b));

    case GeometryType::CylinderEndOn:
    case GeometryType::CylinderSideOn:
    {
      T az, bz, as, bs;
      if( !line_slab_interval_imp( o[2], d[2], dims[1], az, bz ) )
        return false;
      if( !line_quadric_interval_imp<2>( o, d, dims[0], as, bs ) )
        return false;
      a = select_max( az, as );
      b = select_min( bz, bs );
      return scalar_of(a) < scalar_of(b);
    }

    case GeometryType::Rectangular:
    {
      T ax, bx, ay, by, az, bz;
      if( !line_slab_interval_imp( o[0], d[0], dims[0], ax, bx )
          || !line_slab_interval_imp( o[1], d[1], dims[1], ay, by )
          || !line_slab_interval_imp( o[2], d[2], dims[2], az, bz ) )
        return false;
      a = select_max( ax, select_max( ay, az ) );
      b = select_min( bx, select_min( by, bz ) );
      return scalar_of(a) < scalar_of(b);
    }

    case GeometryType::NumGeometryType:
      break;
  }//switch( geometry )

  assert( 0 );
  return false;
}//line_shell_interval_imp(...)


/** Intervals [a_l, b_l] of the line point(s) = o - s d through every shell of a nested stack,
 outermost first: `crossed[l]` is set for each shell the line enters, a_l is clamped at zero (nothing
 on the far side of `o`, which is a detector-side point), and once a shell is missed every shell
 inside it is too.  The single place both the per-line integration and the per-point shell walk
 (#shell_path_to_point_imp) get their intervals from.

 RE-ORIGIN.  `o` is a whole source-detector distance from the source, and computed from it the
 intervals lose precision as the source shrinks - quadratically for the quadrics, linearly for the
 slabs (see #line_quadric_interval_imp).  So every interval is computed from the foot of the
 perpendicular from the assembly origin, o' = o - s* d with s* = (o.d)/(d.d), whose magnitude is of
 order the source extent for any line that hits the source, and shifted back by s*.  The only
 operation that still sees a distance-sized term is the one subtraction forming o', which costs
 eps*(distance/extent) - the linear sensitivity #sm_line_path_extent_ratio_floor is set by. */
template<typename T, typename D>
inline void line_shell_intervals_imp( const GeometryType geometry,
                                      const std::vector<typename DistributedSrcCalcT<T>::ShellInfo> &shells,
                                      const T o[3], const D d[3],
                                      std::vector<T> &a, std::vector<T> &b, std::vector<char> &crossed )
{
  const size_t num_shells = shells.size();
  a.resize( num_shells );
  b.resize( num_shells );
  crossed.assign( num_shells, 0 );

  using AType = std::conditional_t<std::is_same_v<D,double>,double,T>;
  AType dd( 0.0 );
  T od( 0.0 );
  for( int i = 0; i < 3; ++i )
  {
    dd += d[i]*d[i];
    od += o[i]*d[i];
  }
  const T s_star = od / dd;
  const T o_near[3] = { o[0] - s_star*d[0], o[1] - s_star*d[1], o[2] - s_star*d[2] };

  for( size_t l = num_shells; l-- > 0; )
  {
    T al, bl;
    if( !line_shell_interval_imp( geometry, shells[l].dims, o_near, d, al, bl ) )
      break;
    al += s_star;
    bl += s_star;
    // A shell entirely behind `o` cannot happen for a detector-side origin (pastDetector guard);
    //  treat it as missed anyway.
    if( scalar_of(bl) <= 0.0 )
      break;
    if( scalar_of(al) < 0.0 )
      al = T(0.0);
    a[l] = al;
    b[l] = bl;
    crossed[l] = 1;
  }
}//line_shell_intervals_imp(...)


/** The per-shell path a photon emitted at a point travels through on its way to the detector-side
 point `o` along the straight line between them - the CENTRE-RAY convention of the element path's
 shell walkers (eval_cylinder / eval_rect / eval_spherical), reproduced from the line intervals so
 that the cascade correction (#build_cascade_field) and the effective-shielding report get the same
 quantities on the line path that `record_path` / `record_generic` accumulate on the element path.
 Filled by #shell_path_to_point_imp; consumed by #record_shell_path_imp.

 With len_l = max(0, min(b_l, s_end) - a_l) the length of shell l's interval on the near side of the
 emission point (s_end from `o`), the element conventions map as

     element walker                                   here
     inner Generic shell counted only when crossed    crossed[l]  (len_l > 0)
     outer Generic shell unconditional                an outer shell always has len_l > 0
     outer Material shell = its near segment          own_len[l] = len_l - len_{l-1} = a_{l-1} - a_l
     inner Material shell = full nested chord         own_len[l] = len_l - len_{l-1}
     source shell = its own near piece                own_len[m] = (s_end - a_m) - len_{m-1}
     air = outermost shell's exit -> o                air = a[N-1]

 A point inside an inner core (a cascade-field node can be one) is handled by the same rule: its
 own shell's near piece is clipped at the point, and everything outside it is walked as usual.
 `ShellWalkMatchesElementCentreRay` (test_VolumetricLadder.cpp) pins this table against the element
 walkers directly - that test, not shared code, is what keeps the two in step. */
template<typename T>
struct ShellPathT
{
  std::vector<T> own_len;      //per shell: the path through that shell's OWN material
  std::vector<char> crossed;   //per shell: whether the path enters the shell's dims at all
  T air = T(0.0);              //from the outermost shell's exit to `o` (0 when nothing is crossed)
};//struct ShellPathT


/** Fills #ShellPathT for the point at distance `s_end` back along the line o - s d (d the unit
 photon direction, i.e. pointing from the source toward `o`; `d` may be T-valued). */
template<typename T>
inline void shell_path_to_point_imp( const GeometryType geometry,
                                     const std::vector<typename DistributedSrcCalcT<T>::ShellInfo> &shells,
                                     const T o[3], const T d[3], const T &s_end,
                                     ShellPathT<T> &out )
{
  const size_t num_shells = shells.size();
  std::vector<T> a, b;
  std::vector<char> crossed;
  line_shell_intervals_imp( geometry, shells, o, d, a, b, crossed );

  out.own_len.assign( num_shells, T(0.0) );
  out.crossed.assign( num_shells, 0 );
  out.air = T(0.0);

  T len_inner( 0.0 );
  for( size_t l = 0; l < num_shells; ++l )
  {
    if( !crossed[l] )
      continue;
    const T len = select_min( b[l], s_end ) - a[l];
    if( scalar_of(len) <= 0.0 )
      continue;   //the shell lies entirely beyond the emission point
    out.own_len[l] = len - len_inner;
    out.crossed[l] = 1;
    len_inner = len;
  }

  if( num_shells && crossed[num_shells-1] )
    out.air = a[num_shells-1];
}//shell_path_to_point_imp(...)


/** The same walk for a calculator's own geometry: from the emission point `point` (assembly frame)
 to the detector position, using the calculator's shells. */
template<typename T>
inline void shell_path_from_point_imp( const DistributedSrcCalcT<T> &calc, const T point[3],
                                       ShellPathT<T> &out )
{
  using namespace std;
  using namespace ceres;
  const T o[3] = { calc.m_detector.position[0], calc.m_detector.position[1], calc.m_detector.position[2] };
  T d[3] = { o[0] - point[0], o[1] - point[1], o[2] - point[2] };
  const T s_end = sqrt( d[0]*d[0] + d[1]*d[1] + d[2]*d[2] );
  for( int i = 0; i < 3; ++i )
    d[i] = d[i] / s_end;
  shell_path_to_point_imp( calc.m_geometry, calc.m_shells, o, d, s_end, out );
}//shell_path_from_point_imp(...)


/** Loads a #ShellPathT into the calculator's per-ray scratch exactly as the element walkers do on
 their centre ray: reset, then `record_path` / `record_generic` per crossed shell in index order,
 then the air distance.  After this `cascade_correction_factor` can be called for that point. */
template<typename T>
inline void record_shell_path_imp( const DistributedSrcCalcT<T> &calc, const ShellPathT<T> &path )
{
  calc.reset_ray_accumulators();
  for( size_t l = 0; l < calc.m_shells.size(); ++l )
  {
    if( !path.crossed[l] )
      continue;
    const typename DistributedSrcCalcT<T>::ShellInfo &shell = calc.m_shells[l];
    switch( shell.type )
    {
      case ShellType::Material:
        calc.record_path( shell, path.own_len[l] );
        break;
      case ShellType::Generic:
        calc.record_generic( shell );
        break;
    }
  }
  if( calc.m_cascade )
    calc.m_ray_air_dist = path.air;
}//record_shell_path_imp(...)


/** The cascade-summing correction C(r) of one calculator, tabulated on a coarse grid over the
 source and interpolated along every chord - how the line path applies `cascade_correction_factor`,
 which the element path evaluates once per element from a centre-ray walk.

 WHY A FIELD.  The correction needs, per emission point, the per-partner optical depths of the
 centre ray to the detector (the scratch `record_path` fills) and then a run of the cascade engine
 for this calculator's window; per line node that is far too expensive, but C is a smooth ratio of
 order 0.8-1 over the source, so a handful of nodes and trilinear interpolation carry it.  The grid
 is in NORMALISED source coordinates, so it rides on the (T-valued) dimensions and every node's C
 keeps its derivative lanes: the walk, the flat-disk factor and the engine all run in T.

 NODES are equally spaced INCLUDING the end points (interpolation, not quadrature: the face nearest
 the detector is where C changes fastest, and Gauss-Legendre nodes would leave it to clamped
 extrapolation).  On-axis detectors get the symmetric reductions the element integrator uses (no
 theta axis for a cylinder, |x|,|y| for a box), so 16-36 nodes; off-axis 64.  A hollow shell is
 gridded over the OUTER solid: nodes inside the core are walked like any other point (the walk is
 well defined from there) and only anchor the interpolation, since no chord node lies in the core.
 Trilinear interpolation smears the kink C has at the core's silhouette - second order in a
 correction that is itself a few percent; `LineVsElementCascadeField` measures it on a hollow case.
 */
template<typename T>
struct CascadeFieldT
{
  GeometryType geometry = GeometryType::NumGeometryType;
  bool on_axis = false;
  bool periodic1 = false;     //axis 1 is the cylinder azimuth: nodes at k/n, wrapping
  int n[3] = { 1, 1, 1 };
  std::vector<T> C;           //[i0][i1][i2]

  bool empty() const { return C.empty(); }

  size_t index( const int i0, const int i1, const int i2 ) const
  {
    return static_cast<size_t>( (i0*n[1] + i1)*n[2] + i2 );
  }

  /** Trilinear C at normalised coordinates u (each in [0,1]; the periodic axis wraps). */
  T eval( const T u[3] ) const
  {
    int i0[3], i1[3];
    T t[3];
    for( int ax = 0; ax < 3; ++ax )
    {
      const int na = n[ax];
      if( na <= 1 )
      {
        i0[ax] = i1[ax] = 0;
        t[ax] = T(0.0);
        continue;
      }
      if( periodic1 && (ax == 1) )
      {
        T v = u[ax] * T(static_cast<double>(na));
        const double vs = scalar_of( v );
        const double wrapped = vs - std::floor( vs/na )*na;
        const int lo = std::min( na - 1, std::max( 0, static_cast<int>( std::floor( wrapped ) ) ) );
        i0[ax] = lo;
        i1[ax] = (lo + 1) % na;
        t[ax] = v - T( vs - (wrapped - lo) );   //fractional part, derivative lane kept
        continue;
      }
      const T v = u[ax] * T(static_cast<double>(na - 1));
      const double vs = scalar_of( v );
      int lo = static_cast<int>( std::floor( vs ) );
      lo = std::max( 0, std::min( na - 2, lo ) );
      i0[ax] = lo;
      i1[ax] = lo + 1;
      t[ax] = v - T(static_cast<double>(lo));
      // Clamp (flat) outside the grid; the sources never put a chord node there by more than round-off.
      if( scalar_of(t[ax]) < 0.0 )
        t[ax] = T(0.0);
      else if( scalar_of(t[ax]) > 1.0 )
        t[ax] = T(1.0);
    }

    const T one( 1.0 );
    T val( 0.0 );
    for( int c = 0; c < 8; ++c )
    {
      const bool b0 = c & 1, b1 = c & 2, b2 = c & 4;
      const T w = (b0 ? t[0] : one - t[0]) * (b1 ? t[1] : one - t[1]) * (b2 ? t[2] : one - t[2]);
      val += w * C[index( b0 ? i1[0] : i0[0], b1 ? i1[1] : i0[1], b2 ? i1[2] : i0[2] )];
    }
    return val;
  }//eval(...)
};//struct CascadeFieldT


/** Normalised-coordinate <-> assembly-frame maps of #CascadeFieldT for a calculator's source shell
 (the OUTER solid of shell `m`; the sphere's radial axis spans the shell itself). */
template<typename T>
struct CascadeFieldFrameT
{
  GeometryType geometry;
  bool on_axis;
  std::array<T,3> outer;     //outer dims of the source shell
  T inner_radius;            //sphere only: inner radius (0 for a full sphere)
  T det_dir[3];              //unit direction to the detector (sphere polar axis)
  T det_perp[3];             //a unit vector perpendicular to it

  CascadeFieldFrameT( const DistributedSrcCalcT<T> &calc )
  {
    using namespace std;
    using namespace ceres;
    geometry = calc.m_geometry;
    outer = calc.m_shells[calc.m_materialIndex].dims;
    inner_radius = ((geometry == GeometryType::Spherical) && (calc.m_materialIndex > 0))
                     ? calc.m_shells[calc.m_materialIndex-1].dims[0] : T(0.0);
    on_axis = (scalar_of(calc.m_detector.position[0]) == 0.0) && (scalar_of(calc.m_detector.position[1]) == 0.0);
    if( geometry == GeometryType::CylinderSideOn )
      on_axis = (scalar_of(calc.m_detector.position[1]) == 0.0) && (scalar_of(calc.m_detector.position[2]) == 0.0);

    T dist( 0.0 );
    for( int i = 0; i < 3; ++i )
      dist += calc.m_detector.position[i]*calc.m_detector.position[i];
    dist = sqrt( dist );
    for( int i = 0; i < 3; ++i )
      det_dir[i] = calc.m_detector.position[i] / dist;
    // Any perpendicular: cross with the axis least aligned with det_dir.
    const int k = (std::fabs(scalar_of(det_dir[0])) < 0.6) ? 0 : ((std::fabs(scalar_of(det_dir[1])) < 0.6) ? 1 : 2);
    T e[3] = { T(0.0), T(0.0), T(0.0) };
    e[k] = T(1.0);
    T p[3] = { det_dir[1]*e[2] - det_dir[2]*e[1], det_dir[2]*e[0] - det_dir[0]*e[2], det_dir[0]*e[1] - det_dir[1]*e[0] };
    const T pn = sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] );
    for( int i = 0; i < 3; ++i )
      det_perp[i] = p[i] / pn;
  }

  /** Grid shape and periodicity for this frame. */
  void shape( int n[3], bool &periodic1 ) const
  {
    periodic1 = false;
    switch( geometry )
    {
      case GeometryType::Spherical:
        n[0] = 4; n[1] = 4; n[2] = 1;
        break;
      case GeometryType::CylinderEndOn:
      case GeometryType::CylinderSideOn:
        n[0] = 4; n[1] = on_axis ? 1 : 4; n[2] = 4;
        periodic1 = !on_axis;
        break;
      case GeometryType::Rectangular:
        n[0] = on_axis ? 3 : 4; n[1] = on_axis ? 3 : 4; n[2] = 4;
        break;
      case GeometryType::NumGeometryType:
        assert( 0 );
        n[0] = n[1] = n[2] = 1;
        break;
    }
  }

  /** Assembly-frame point of normalised coordinates u. */
  void point( const double u[3], T p[3] ) const
  {
    using namespace std;
    using namespace ceres;
    const double two_pi = 2.0*PhysicalUnits::pi;
    switch( geometry )
    {
      case GeometryType::Spherical:
      {
        const T r = inner_radius + T(u[0])*(outer[0] - inner_radius);
        const double ct = 2.0*u[1] - 1.0;
        const double st = std::sqrt( std::max( 0.0, 1.0 - ct*ct ) );
        for( int i = 0; i < 3; ++i )
          p[i] = r*(T(ct)*det_dir[i] + T(st)*det_perp[i]);
        break;
      }
      case GeometryType::CylinderEndOn:
      case GeometryType::CylinderSideOn:
      {
        const T r = T(u[0])*outer[0];
        const double th = on_axis ? 0.0 : two_pi*u[1];
        p[0] = r*T(std::cos(th));
        p[1] = r*T(std::sin(th));
        p[2] = T(2.0*u[2] - 1.0)*outer[1];
        break;
      }
      case GeometryType::Rectangular:
        p[0] = on_axis ? T(u[0])*outer[0] : T(2.0*u[0] - 1.0)*outer[0];
        p[1] = on_axis ? T(u[1])*outer[1] : T(2.0*u[1] - 1.0)*outer[1];
        p[2] = T(2.0*u[2] - 1.0)*outer[2];
        break;
      case GeometryType::NumGeometryType:
        assert( 0 );
        break;
    }
  }//point(...)

  /** Normalised coordinates of an assembly-frame point (T). */
  void coords( const T p[3], T u[3] ) const
  {
    using namespace std;
    using namespace ceres;
    const double two_pi = 2.0*PhysicalUnits::pi;
    u[0] = u[1] = u[2] = T(0.0);
    switch( geometry )
    {
      case GeometryType::Spherical:
      {
        const T r = sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] );
        const T dr = outer[0] - inner_radius;
        u[0] = (scalar_of(dr) > 0.0) ? (r - inner_radius)/dr : T(0.0);
        const T ct = (scalar_of(r) > 0.0) ? (p[0]*det_dir[0] + p[1]*det_dir[1] + p[2]*det_dir[2])/r : T(0.0);
        u[1] = T(0.5)*(ct + T(1.0));
        break;
      }
      case GeometryType::CylinderEndOn:
      case GeometryType::CylinderSideOn:
      {
        const T r = sqrt( p[0]*p[0] + p[1]*p[1] );
        u[0] = r / outer[0];
        if( !on_axis )
        {
          T th = atan2( p[1], p[0] );
          if( scalar_of(th) < 0.0 )
            th += T(two_pi);
          u[1] = th / T(two_pi);
        }
        u[2] = T(0.5)*(p[2]/outer[1] + T(1.0));
        break;
      }
      case GeometryType::Rectangular:
        if( on_axis )
        {
          u[0] = ((scalar_of(p[0]) < 0.0) ? -p[0] : p[0]) / outer[0];
          u[1] = ((scalar_of(p[1]) < 0.0) ? -p[1] : p[1]) / outer[1];
        }else
        {
          u[0] = T(0.5)*(p[0]/outer[0] + T(1.0));
          u[1] = T(0.5)*(p[1]/outer[1] + T(1.0));
        }
        u[2] = T(0.5)*(p[2]/outer[2] + T(1.0));
        break;
      case GeometryType::NumGeometryType:
        assert( 0 );
        break;
    }
  }//coords(...)
};//struct CascadeFieldFrameT


/** Builds the cascade fields of a line group: one #CascadeFieldT per calculator (empty for a
 calculator without `m_cascade`).  The per-node shell walk is shared by the group (identical shells
 and partner coefficients - `cascade_mu` is per material and partner energy, not per calculator);
 each calculator then runs its own window through the cascade engine at every node, in parallel
 over calculators when `multithread` (each has its own scratch; the engine's cascade memo is
 mutex-guarded). */
template<typename T>
std::vector<CascadeFieldT<T>> build_cascade_fields( const std::vector<DistributedSrcCalcT<T>*> &group,
                                                    const bool multithread )
{
  using namespace std;

  vector<CascadeFieldT<T>> fields( group.size() );
  vector<size_t> with_cascade;
  for( size_t c = 0; c < group.size(); ++c )
    if( group[c]->m_cascade )
      with_cascade.push_back( c );
  if( with_cascade.empty() )
    return fields;

  const DistributedSrcCalcT<T> &lead = *group.front();
  const CascadeFieldFrameT<T> frame( lead );
  int n[3];
  bool periodic1;
  frame.shape( n, periodic1 );
  const size_t num_nodes = static_cast<size_t>( n[0]*n[1]*n[2] );

  // Node positions, walks and flat-disk factors: once per group.
  vector<ShellPathT<T>> walks( num_nodes );
  vector<T> det_factors( num_nodes );
  for( int i0 = 0; i0 < n[0]; ++i0 )
  {
    for( int i1 = 0; i1 < n[1]; ++i1 )
    {
      for( int i2 = 0; i2 < n[2]; ++i2 )
      {
        const double u[3] = { (n[0] > 1) ? static_cast<double>(i0)/(n[0] - 1) : 0.0,
                              periodic1 ? static_cast<double>(i1)/n[1]
                                        : ((n[1] > 1) ? static_cast<double>(i1)/(n[1] - 1) : 0.0),
                              (n[2] > 1) ? static_cast<double>(i2)/(n[2] - 1) : 0.0 };
        T p[3];
        frame.point( u, p );
        const size_t idx = static_cast<size_t>( (i0*n[1] + i1)*n[2] + i2 );
        shell_path_from_point_imp( lead, p, walks[idx] );
        det_factors[idx] = detector_response_factor( lead.m_detector, p );
      }
    }
  }

  const auto build_one = [&]( const size_t c )
  {
    const DistributedSrcCalcT<T> &calc = *group[c];
    CascadeFieldT<T> &f = fields[c];
    f.geometry = frame.geometry;
    f.on_axis = frame.on_axis;
    f.periodic1 = periodic1;
    for( int ax = 0; ax < 3; ++ax )
      f.n[ax] = n[ax];
    f.C.resize( num_nodes );
    for( size_t idx = 0; idx < num_nodes; ++idx )
    {
      record_shell_path_imp( calc, walks[idx] );
      f.C[idx] = calc.cascade_correction_factor( det_factors[idx], 1.0 );
    }
  };

  if( multithread && (with_cascade.size() > 1) )
  {
    std::mutex error_mutex;
    std::exception_ptr first_error;
    SpecUtilsAsync::ThreadPool pool;
    for( const size_t c : with_cascade )
    {
      pool.post( [c,&build_one,&error_mutex,&first_error](){
        try
        {
          build_one( c );
        }catch( std::exception & )
        {
          std::lock_guard<std::mutex> lock( error_mutex );
          if( !first_error )
            first_error = std::current_exception();
        }
      } );
    }
    pool.join();
    if( first_error )
      std::rethrow_exception( first_error );
  }else
  {
    for( const size_t c : with_cascade )
      build_one( c );
  }

  return fields;
}//build_cascade_fields(...)


/** The response prefactor P(E; d, cos_theta, phi) = exp(ln_eta + ln_N + ln_k) tabulated on a
 crystal-frame grid for ONE energy, so the line integrand can look it up (in T, so the position
 derivative flows) instead of paying a PCHIP evaluation per chord node.  `d` is measured from the
 crystal-face ORIGIN, as ceelo::DetectorResponse::common_eval measures it - NOT from the reference
 point the assembly's detector position stands for.  Axial responses ignore phi (one node);
 Quadrant (box-crystal) responses carry phi over [0,90] degrees and are reflected into it. */
struct PrefactorGrid
{
  double energy_keV = 0.0;
  std::vector<double> ln_d;       //ascending, cm
  std::vector<double> cos_t;      //ascending, in [0,1]
  std::vector<double> phi_deg;    //ascending; {0} for an Axial response
  std::vector<double> ln_p;       //[id][ic][ip], d-major
  std::vector<ceelo::ResponseFlag> flags;   //per node, as ln_p
  ceelo::ResponseFlag worst_flag = ceelo::ResponseFlag::Ok;
  double max_frac_sigma = 0.0;

  /** The worst flag over the nodes the source's chords can actually reach: the grid is padded well
   beyond them in BOTH distance and incidence angle, and a flag raised only in that padding is not
   one the source ever sees.  (A collimated detector is shadowed at 90 degrees whatever it is
   looking at, so a distance-only restriction would flag every source.)  The interval is widened to
   the enclosing grid cell, so a source lying between two nodes still picks up the flag of the
   nodes it sits between.  Phi is not restricted - for a Quadrant response the same source can span
   the whole folded range. */
  ceelo::ResponseFlag worst_flag_between( const double d_lo_cm, const double d_hi_cm,
                                          const double cos_lo, const double cos_hi ) const
  {
    // The nodes bracketing each end of an ascending-node interval.
    const auto bracket = []( const std::vector<double> &nodes, const double lo, const double hi,
                             size_t &i_lo, size_t &i_hi ){
      i_lo = 0;
      i_hi = nodes.empty() ? 0 : (nodes.size() - 1);
      for( size_t i = 0; i < nodes.size(); ++i )
      {
        if( nodes[i] <= lo )
          i_lo = i;
        if( nodes[i] >= hi )
        {
          i_hi = i;
          break;
        }
      }
    };

    size_t id_lo, id_hi, ic_lo, ic_hi;
    bracket( ln_d, std::log( d_lo_cm ), std::log( d_hi_cm ), id_lo, id_hi );
    bracket( cos_t, cos_lo, cos_hi, ic_lo, ic_hi );

    ceelo::ResponseFlag worst = ceelo::ResponseFlag::Ok;
    for( size_t id = id_lo; id <= id_hi; ++id )
      for( size_t ic = ic_lo; ic <= ic_hi; ++ic )
        for( size_t ip = 0; ip < phi_deg.size(); ++ip )
          if( static_cast<int>(flags[index(id,ic,ip)]) > static_cast<int>(worst) )
            worst = flags[index(id,ic,ip)];
    return worst;
  }//worst_flag_between(...)

  size_t index( const size_t id, const size_t ic, const size_t ip ) const
  {
    return (id*cos_t.size() + ic)*phi_deg.size() + ip;
  }

  /** Locates `x` (scalar) in the ascending `nodes`, returning the lower node index and the
   T-valued fractional position, clamped to the grid (zero slope past either end). */
  template<typename T>
  static void locate( const std::vector<double> &nodes, const T &x, size_t &i0, T &t )
  {
    const double xs = scalar_of( x );
    if( nodes.size() < 2 )
    {
      i0 = 0;
      t = T(0.0);
      return;
    }
    if( xs <= nodes.front() )
    {
      i0 = 0;
      t = T(0.0);
      return;
    }
    if( xs >= nodes.back() )
    {
      i0 = nodes.size() - 2;
      t = T(1.0);
      return;
    }
    const size_t hi = std::upper_bound( begin(nodes), end(nodes), xs ) - begin(nodes);
    i0 = hi - 1;
    t = (x - nodes[i0]) / (nodes[hi] - nodes[i0]);
  }//locate(...)

  /** P at (ln d, cos_theta, phi_deg), trilinear in ln P. */
  template<typename T>
  T eval( const T &ln_d_q, const T &cos_t_q, const T &phi_q ) const
  {
    using namespace std;
    using namespace ceres;

    size_t id, ic, ip = 0;
    T td, tc, tp( 0.0 );
    locate( ln_d, ln_d_q, id, td );
    locate( cos_t, cos_t_q, ic, tc );
    const bool has_phi = (phi_deg.size() > 1);
    if( has_phi )
      locate( phi_deg, phi_q, ip, tp );

    const size_t id1 = std::min( id + 1, ln_d.size() - 1 );
    const size_t ic1 = std::min( ic + 1, cos_t.size() - 1 );
    const size_t ip1 = has_phi ? std::min( ip + 1, phi_deg.size() - 1 ) : ip;

    // Bilinear in (d, cos_theta) at each of the (up to two) phi planes, then linear in phi.
    const auto plane = [&]( const size_t p ) -> T {
      const T v00 = T(ln_p[index(id,  ic,  p)]);
      const T v10 = T(ln_p[index(id1, ic,  p)]);
      const T v01 = T(ln_p[index(id,  ic1, p)]);
      const T v11 = T(ln_p[index(id1, ic1, p)]);
      return (T(1.0) - td)*((T(1.0) - tc)*v00 + tc*v01) + td*((T(1.0) - tc)*v10 + tc*v11);
    };

    T ln_val = plane( ip );
    if( has_phi && (ip1 != ip) )
      ln_val = (T(1.0) - tp)*ln_val + tp*plane( ip1 );

    return exp( ln_val );
  }//eval(...)
};//struct PrefactorGrid


/** The aperture quadratures a COLLIMATED response's shadow gate needs, on a coarse position grid.

 `DetectorResponse::common_eval` gates a collimated response on the transmitted fraction of an
 aperture quadrature at the query position - below 0.3 it flags Shadowed and inflates sigma, below
 0.05 NeedsMc - and skips the gate silently when handed an empty quadrature.  The gate moves only
 the flag and the sigma, never the value (the collimator is part of the geometry every line is
 traced through, so the kernel already carries the shadow), and its thresholds do not need position
 resolution, so a handful of quadratures at a reduced ray count serve the whole prefactor grid:
 #build_prefactor_grid hands each node the nearest one.  Position-only, so built once per line set. */
struct CollimatorGateGrid
{
  std::vector<double> ln_d, cos_t, phi_deg;
  std::vector<ceelo::ApertureQuadrature> q;   //[id][ic][ip]

  size_t index( const size_t id, const size_t ic, const size_t ip ) const
  {
    return (id*cos_t.size() + ic)*phi_deg.size() + ip;
  }

  static size_t nearest_index( const std::vector<double> &nodes, const double x )
  {
    size_t best = 0;
    for( size_t i = 1; i < nodes.size(); ++i )
      if( std::fabs(nodes[i] - x) < std::fabs(nodes[best] - x) )
        best = i;
    return best;
  }

  const ceelo::ApertureQuadrature &nearest( const double ln_d_q, const double cos_q, const double phi_q ) const
  {
    return q[index( nearest_index( ln_d, ln_d_q ), nearest_index( cos_t, cos_q ),
                    nearest_index( phi_deg, phi_q ) )];
  }
};//struct CollimatorGateGrid


inline std::shared_ptr<const CollimatorGateGrid> build_collimator_gate_grid( const ceelo::DetectorResponse &resp,
                                                                            const double d_lo_cm,
                                                                            const double d_hi_cm,
                                                                            const size_t num_d = 6,
                                                                            const size_t num_cos = 5,
                                                                            const size_t num_phi_quadrant = 3,
                                                                            const int n_rays = 256 )
{
  auto grid = std::make_shared<CollimatorGateGrid>();
  const double lo = std::log( d_lo_cm ), hi = std::log( d_hi_cm );
  for( size_t i = 0; i < num_d; ++i )
    grid->ln_d.push_back( lo + (hi - lo)*static_cast<double>(i)/static_cast<double>(num_d - 1) );
  for( size_t i = 0; i < num_cos; ++i )
    grid->cos_t.push_back( static_cast<double>(i)/static_cast<double>(num_cos - 1) );
  if( resp.descriptor.symmetry == ceelo::ResponseSymmetry::Quadrant )
  {
    for( size_t i = 0; i < num_phi_quadrant; ++i )
      grid->phi_deg.push_back( 90.0*static_cast<double>(i)/static_cast<double>(num_phi_quadrant - 1) );
  }else
  {
    grid->phi_deg.push_back( 0.0 );
  }

  grid->q.resize( grid->ln_d.size() * grid->cos_t.size() * grid->phi_deg.size() );
  for( size_t id = 0; id < grid->ln_d.size(); ++id )
    for( size_t ic = 0; ic < grid->cos_t.size(); ++ic )
      for( size_t ip = 0; ip < grid->phi_deg.size(); ++ip )
      {
        const Eigen::Vector3d pos = ceelo::source_position( std::exp( grid->ln_d[id] ), grid->cos_t[ic],
                                                            grid->phi_deg[ip] * PhysicalUnits::pi / 180.0 );
        grid->q[grid->index(id,ic,ip)] = ceelo::make_aperture_quadrature( resp.geometry(), pos, n_rays );
      }
  return grid;
}//build_collimator_gate_grid(...)


/** Builds #PrefactorGrid for `energy_keV` over distances [d_lo_cm, d_hi_cm] (log-spaced) and
 cos_theta in [0,1]; phi over [0,90] deg for Quadrant responses.  `gates` (collimated responses)
 supplies the shadow gate's quadrature per node, see #CollimatorGateGrid. */
inline std::shared_ptr<const PrefactorGrid> build_prefactor_grid( const ceelo::DetectorResponse &resp,
                                                                  const double energy_keV,
                                                                  const double d_lo_cm,
                                                                  const double d_hi_cm,
                                                                  const CollimatorGateGrid *gates = nullptr,
                                                                  const size_t num_cos = 33,
                                                                  const size_t num_phi_quadrant = 9,
                                                                  const int nodes_per_octave = 12 )
{
  assert( gates || !resp.descriptor.collimator );

  assert( (d_lo_cm > 0.0) && (d_hi_cm > d_lo_cm) );
  auto grid = std::make_shared<PrefactorGrid>();
  grid->energy_keV = energy_keV;

  // Distance nodes on a FIXED logarithmic lattice, ln d = k*delta (12 per octave), covering
  //  [d_lo, d_hi]: a grid rebuilt over a wider range (VolumetricLineCache::ensure_prefactor_range,
  //  as the fitted source grows or shrinks) reproduces every node it shares with the old one bit
  //  for bit, so the integrand is continuous across the extension.
  assert( (num_cos >= 2) && (num_phi_quadrant >= 2) && (nodes_per_octave >= 1) );
  const double delta = std::log( 2.0 ) / static_cast<double>( nodes_per_octave );
  const long k_lo = static_cast<long>( std::floor( std::log( d_lo_cm )/delta ) );
  const long k_hi = std::max( static_cast<long>( std::ceil( std::log( d_hi_cm )/delta ) ), k_lo + 1 );
  for( long k = k_lo; k <= k_hi; ++k )
    grid->ln_d.push_back( static_cast<double>(k) * delta );
  for( size_t i = 0; i < num_cos; ++i )
    grid->cos_t.push_back( static_cast<double>(i)/static_cast<double>(num_cos - 1) );
  if( resp.descriptor.symmetry == ceelo::ResponseSymmetry::Quadrant )
  {
    for( size_t i = 0; i < num_phi_quadrant; ++i )
      grid->phi_deg.push_back( 90.0*static_cast<double>(i)/static_cast<double>(num_phi_quadrant - 1) );
  }else
  {
    grid->phi_deg.push_back( 0.0 );
  }

  // The quadrature argument of fep_prefactor only feeds the collimator shadow gate (flag and sigma,
  //  never the value); an empty one is correct for an uncollimated response.
  const ceelo::ApertureQuadrature no_quadrature;

  grid->ln_p.resize( grid->ln_d.size() * grid->cos_t.size() * grid->phi_deg.size() );
  grid->flags.assign( grid->ln_p.size(), ceelo::ResponseFlag::Ok );
  for( size_t id = 0; id < grid->ln_d.size(); ++id )
  {
    const double d = std::exp( grid->ln_d[id] );
    for( size_t ic = 0; ic < grid->cos_t.size(); ++ic )
    {
      const double ct = grid->cos_t[ic];
      for( size_t ip = 0; ip < grid->phi_deg.size(); ++ip )
      {
        const double phi = grid->phi_deg[ip] * PhysicalUnits::pi / 180.0;
        const Eigen::Vector3d pos = ceelo::source_position( d, ct, phi );
        const ceelo::ApertureQuadrature &gate_q = gates ? gates->nearest( grid->ln_d[id], ct, grid->phi_deg[ip] )
                                                        : no_quadrature;
        const ceelo::EffResult pre = resp.fep_prefactor( energy_keV, pos, gate_q );
        if( !(pre.value > 0.0) )
          throw std::runtime_error( "build_prefactor_grid: non-positive prefactor" );
        grid->ln_p[grid->index(id,ic,ip)] = std::log( pre.value );
        grid->flags[grid->index(id,ic,ip)] = pre.flag;
        if( static_cast<int>(pre.flag) > static_cast<int>(grid->worst_flag) )
          grid->worst_flag = pre.flag;
        grid->max_frac_sigma = std::max( grid->max_frac_sigma, pre.sigma / pre.value );
      }
    }
  }

  return grid;
}//build_prefactor_grid(...)


/** Fraction of the line set aimed at the source's SURFACES rather than its volume - the defensive
 component of the mixture proposal (see DIRECTION PROPOSAL in the file comment).  Part of the
 cache's key. */
inline double sm_default_volumetric_line_surface_frac = 0.3;

/** TEST HOOK - when set, #VolumetricLineCache::traced returns the most recently traced set whatever
 dimensions are asked for, and the integration carries the crystal kernel WITHOUT its direction
 gradient, so a finite difference of the estimator measures the chain rule of everything analytic
 (the chain-rule lane of LinePathGradientVsFiniteDifference).  Leave false in production. */
inline bool sm_line_trace_hold = false;

/** TEST HOOKS - resolution of the prefactor grid (#build_prefactor_grid): cos_theta nodes on [0,1],
 phi nodes over a quadrant, and ln(d) nodes per octave.  These are the production values; the
 tests sweep them to measure the interpolation error of P on its own. */
inline size_t sm_prefactor_grid_num_cos = 33;
inline size_t sm_prefactor_grid_num_phi_quadrant = 9;
inline int sm_prefactor_grid_nodes_per_octave = 12;

/** TEST HOOK - evaluate the response prefactor DIRECTLY (ceelo::DetectorResponse::fep_prefactor) at
 every chord node instead of interpolating it from #PrefactorGrid.  T = double only - the direct
 evaluation carries no derivative lane - so it throws for a Jet.  Isolates the grid's interpolation
 error from everything else in the line integrand. */
inline bool sm_prefactor_direct_eval = false;

/** Number of contiguous index blocks the line set is summed in (#line_source_integration_imp):
 fixes the reduction order independently of the thread count, and the spread of the block estimates
 is the quadrature's error estimate (`m_est_rel_error`).  Not a tuning knob - changing it changes
 the rounding of every line-path result. */
inline size_t sm_line_error_blocks = 32;

/** The HEMISPHERE component of the mixture proposal (see #VolumetricLineCache::frac_hemi): enabled
 at build for a WIDE source - one whose padded extent, seen from the detector, spans a far/near
 distance ratio above `sm_volumetric_line_hemi_ratio` (near floored at the crystal's transverse
 extent) - with `sm_volumetric_line_hemi_frac` of the lines, taken from the volume share.
 MEASURED (LineProposalMixtureSweep, replica rms of 8 Sobol' replicas at 65536 lines): on the 20 m
 in-situ disk at 1 m (ratio 16) a share of 0 / 0.15 / 0.3 / 0.5 / 0.65 gives 4.0e-3 / 6.0e-4 /
 5.3e-4 / 3.5e-4 / 2.9e-4, on the 100 m disk (ratio 76) 1.5e-2 / 9.4e-4 / 1.2e-3 / 6.9e-4 / 6.7e-4;
 on compact sources (contact, far field, a 2 m disk at 1 m: ratios 2.3-3.4) the same shares cost
 1.0-2.5x, hence the trigger, set above the contact box's 3.4.  An inverse-square area density on
 the facing face was tried for the same purpose and did less (2.4e-3 / 3.6e-3 on the two disks)
 while hurting boxes; removed. */
inline double sm_volumetric_line_hemi_ratio = 4.0;
inline double sm_volumetric_line_hemi_frac = 0.5;

// The sequence/replica selector of a line set, #LineSampleParams, is declared in
//  GammaInteractionCalc.h (ShieldingSourceChi2Fcn holds one per fit).


/** Right-handed orthonormal tangent basis (e1, e2) perpendicular to the unit vector `a` - a
 deterministic function of `a` alone, so the pass that traces a perturbed direction and the pass
 that projects a derivative onto it agree without sharing state. */
inline void line_tangent_frame( const Eigen::Vector3d &a, Eigen::Vector3d &e1, Eigen::Vector3d &e2 )
{
  e1 = (std::fabs(a.z()) < 0.9) ? a.cross( Eigen::Vector3d(0.0, 0.0, 1.0) ).normalized()
                                : a.cross( Eigen::Vector3d(1.0, 0.0, 0.0) ).normalized();
  e2 = a.cross( e1 );
}//line_tangent_frame(...)


/** Traces one line through the DETECTOR (housing, dead layer, crystal) and fills the CeeLo ray the
 per-energy kernel is evaluated on: the material segments in traversal order, the active-crystal
 chord, and the photon direction.  `x` is the hull point and `w_out` the OUTWARD unit direction
 (toward the source; the photon travels along -w_out).  `scratch` is reused across calls.

 Returns false - leaving the ray with no segments and a zero chord, which
 `DetectorResponse::fep_line_probabilities` scores as zero - when the direction points behind the
 face plane or the line misses the active crystal.  Both are decisions a moving line set has to be
 able to make per evaluation, so neither is an assert.

 This is `ceelo::append_etendue_line` without the set: same trace, same conventions, but writing
 into a caller-owned ray so the lines can be re-traced as the source dimensions move.  Kept here
 rather than added to CeeLo so the vendored library stays untouched; it uses only its public API
 (Geometry::trace_ray and the two extent queries).  `Geometry::trace_ray` is const and the geometry
 has no mutable state, so concurrent calls with distinct `out`/`scratch` are safe. */
inline bool trace_detector_line( const ceelo::Geometry &geom, const Eigen::Vector3d &x,
                                 const Eigen::Vector3d &w_out, ceelo::KernelRay &out,
                                 std::vector<ceelo::PathSegment> &scratch )
{
  out.segs.clear();
  out.active_len = 0.0f;
  out.omega_w = static_cast<float>( 1.0/(4.0*PhysicalUnits::pi) );
  out.cos_incidence = static_cast<float>( std::fabs( w_out.z() ) );
  out.dir = (-w_out).cast<float>();
  if( !(w_out.z() < 0.0) )
    return false;   //behind the face plane

  // The trace must start outside the outermost shell, on the source side.
  const std::pair<double,double> zext = geom.outer_z_extent();
  const double back = 2.0*(geom.outer_bounding_radius()
                           + std::max( std::fabs(zext.first), std::fabs(zext.second) ))
                      + x.norm() + 1.0;
  const Eigen::Vector3d dir_in = -w_out;
  geom.trace_ray( x + w_out*back, dir_in, scratch );
  if( scratch.empty() )
    return false;
  std::sort( begin(scratch), end(scratch),
             []( const ceelo::PathSegment &a, const ceelo::PathSegment &b ){
               return a.t_start < b.t_start;
             } );

  double active = 0.0;
  for( const ceelo::PathSegment &seg : scratch )
  {
    const double len = seg.length();
    if( len <= 1.0e-12 )
      continue;
    if( seg.is_scoring )
      active += len;
    if( seg.material )
      out.segs.push_back( { seg.material, static_cast<float>(len), seg.is_scoring } );
  }

  out.active_len = static_cast<float>( active );
  if( !(active > 0.0) || out.segs.empty() )
  {
    out.segs.clear();
    out.active_len = 0.0f;
    return false;   //can never contribute
  }
  return true;
}//trace_detector_line(...)


/** Sets #sm_line_trace_hold for a scope and restores it on exit. */
struct ScopedLineTraceHold
{
  const bool previous;
  ScopedLineTraceHold() : previous( sm_line_trace_hold ) { sm_line_trace_hold = true; }
  ~ScopedLineTraceHold() { sm_line_trace_hold = previous; }
  ScopedLineTraceHold( const ScopedLineTraceHold & ) = delete;
  ScopedLineTraceHold &operator=( const ScopedLineTraceHold & ) = delete;
};//struct ScopedLineTraceHold


/** Whether a T carries a non-zero derivative lane (false for double). */
template<typename T>
inline bool has_derivative_lane( const T &x )
{
  if constexpr( std::is_same_v<T,double> )
    return false;
  else
    return (x.v.squaredNorm() > 0.0);
}


/** The crystal -> assembly rotation for a detector whose axis (detector -> assembly) is `axis`,
 with the crystal's +x rotated by `azimuth` about the axis from the reference direction (assembly
 +x, or +y when the axis is along x).  CeeLo's detector axis is -z, so M e_z = -axis. */
inline void detector_frame_rotation( const double axis[3], const double azimuth, double M[3][3] )
{
  const double a[3] = { axis[0], axis[1], axis[2] };
  // Reference transverse direction: assembly +x unless the axis is along x.
  double u0[3] = { 1.0, 0.0, 0.0 };
  if( std::fabs(a[0]) > 0.9 )
  {
    u0[0] = 0.0;
    u0[1] = 1.0;
  }
  // Project out the axis component and normalise.
  const double ua = u0[0]*a[0] + u0[1]*a[1] + u0[2]*a[2];
  for( int i = 0; i < 3; ++i )
    u0[i] -= ua*a[i];
  const double un = std::sqrt( u0[0]*u0[0] + u0[1]*u0[1] + u0[2]*u0[2] );
  for( int i = 0; i < 3; ++i )
    u0[i] /= un;
  // v0 = (-a) x u0 so that [u0 v0 -a] is right-handed.
  const double na[3] = { -a[0], -a[1], -a[2] };
  const double v0[3] = { na[1]*u0[2] - na[2]*u0[1], na[2]*u0[0] - na[0]*u0[2], na[0]*u0[1] - na[1]*u0[0] };
  const double c = std::cos( azimuth ), s = std::sin( azimuth );
  for( int i = 0; i < 3; ++i )
  {
    M[i][0] = c*u0[i] + s*v0[i];
    M[i][1] = -s*u0[i] + c*v0[i];
    M[i][2] = na[i];
  }
}//detector_frame_rotation(...)


/** Uniform point in the (outer) source solid of the given geometry and dims (assembly frame). */
inline std::array<double,3> uniform_point_in_solid( const GeometryType geometry,
                                                    const std::array<double,3> &dims,
                                                    const double u1, const double u2, const double u3 )
{
  const double two_pi = 2.0*PhysicalUnits::pi;
  switch( geometry )
  {
    case GeometryType::Spherical:
    {
      const double r = dims[0] * std::cbrt( u1 );
      const double ct = 1.0 - 2.0*u2;
      const double st = std::sqrt( std::max( 0.0, 1.0 - ct*ct ) );
      const double ph = two_pi*u3;
      return { r*st*std::cos(ph), r*st*std::sin(ph), r*ct };
    }
    case GeometryType::CylinderEndOn:
    case GeometryType::CylinderSideOn:
    {
      const double r = dims[0] * std::sqrt( u1 );
      const double ph = two_pi*u2;
      return { r*std::cos(ph), r*std::sin(ph), (2.0*u3 - 1.0)*dims[1] };
    }
    case GeometryType::Rectangular:
      return { (2.0*u1 - 1.0)*dims[0], (2.0*u2 - 1.0)*dims[1], (2.0*u3 - 1.0)*dims[2] };
    case GeometryType::NumGeometryType:
      break;
  }
  assert( 0 );
  return { 0.0, 0.0, 0.0 };
}//uniform_point_in_solid(...)


/** Volume of the solid (assembly frame dims). */
template<typename T>
inline T solid_volume( const GeometryType geometry, const std::array<T,3> &dims )
{
  const T pi( PhysicalUnits::pi );
  switch( geometry )
  {
    case GeometryType::Spherical:      return T(4.0/3.0)*pi*dims[0]*dims[0]*dims[0];
    case GeometryType::CylinderEndOn:
    case GeometryType::CylinderSideOn: return T(2.0)*pi*dims[0]*dims[0]*dims[1];
    case GeometryType::Rectangular:    return T(8.0)*dims[0]*dims[1]*dims[2];
    case GeometryType::NumGeometryType: break;
  }
  assert( 0 );
  return T(0.0);
}//solid_volume(...)


/** Largest distance from the assembly origin to a point of the solid. */
inline double solid_bounding_radius( const GeometryType geometry, const std::array<double,3> &dims )
{
  switch( geometry )
  {
    case GeometryType::Spherical:      return dims[0];
    case GeometryType::CylinderEndOn:
    case GeometryType::CylinderSideOn: return std::hypot( dims[0], dims[1] );
    case GeometryType::Rectangular:    return std::sqrt( dims[0]*dims[0] + dims[1]*dims[1] + dims[2]*dims[2] );
    case GeometryType::NumGeometryType: break;
  }
  assert( 0 );
  return 0.0;
}//solid_bounding_radius(...)


/** The component-wise scale that maps the UNIT solid (every dim 1) onto the solid with `dims`:
 (R,R,R) for a sphere, (R,R,H) for a cylinder, (W,H,D) for a box.  #uniform_point_in_solid and
 #unit_surface_point both factor this way, so a point frozen in unit coordinates follows the
 fitted dimensions by one multiplication per axis. */
template<typename T>
inline std::array<T,3> solid_scale( const GeometryType geometry, const std::array<T,3> &dims )
{
  switch( geometry )
  {
    case GeometryType::Spherical:      return { dims[0], dims[0], dims[0] };
    case GeometryType::CylinderEndOn:
    case GeometryType::CylinderSideOn: return { dims[0], dims[0], dims[1] };
    case GeometryType::Rectangular:    return dims;
    case GeometryType::NumGeometryType: break;
  }
  assert( 0 );
  return dims;
}//solid_scale(...)


/** Faces of a solid's surface: sphere 1; cylinders 3 (cap +z, cap -z, side); box 6 (+x, -x, +y, -y,
 +z, -z).  The surface-aimed lines are sampled per face with FIXED probabilities and their density
 needs the area of the face a line crosses, so both index faces the same way. */
inline int surface_face_count( const GeometryType geometry )
{
  switch( geometry )
  {
    case GeometryType::Spherical:      return 1;
    case GeometryType::CylinderEndOn:
    case GeometryType::CylinderSideOn: return 3;
    case GeometryType::Rectangular:    return 6;
    case GeometryType::NumGeometryType: break;
  }
  assert( 0 );
  return 0;
}//surface_face_count(...)


/** Area of face `f` of the solid with `dims`. */
template<typename T>
inline T surface_face_area( const GeometryType geometry, const int f, const std::array<T,3> &dims )
{
  const T pi( PhysicalUnits::pi );
  switch( geometry )
  {
    case GeometryType::Spherical:
      return T(4.0)*pi*dims[0]*dims[0];
    case GeometryType::CylinderEndOn:
    case GeometryType::CylinderSideOn:
      return (f < 2) ? (pi*dims[0]*dims[0]) : (T(4.0)*pi*dims[0]*dims[1]);
    case GeometryType::Rectangular:
    {
      const int ax = f/2;
      return T(4.0)*dims[(ax + 1)%3]*dims[(ax + 2)%3];
    }
    case GeometryType::NumGeometryType:
      break;
  }
  assert( 0 );
  return T(0.0);
}//surface_face_area(...)


/** Point on face `f` of the UNIT solid, uniform in the face's area, from two unit-square
 coordinates; scaled by #solid_scale it lies on the same face of any solid of that geometry. */
inline std::array<double,3> unit_surface_point( const GeometryType geometry, const int f,
                                                const double u1, const double u2 )
{
  const double two_pi = 2.0*PhysicalUnits::pi;
  switch( geometry )
  {
    case GeometryType::Spherical:
    {
      const double ct = 1.0 - 2.0*u1;
      const double st = std::sqrt( std::max( 0.0, 1.0 - ct*ct ) );
      const double ph = two_pi*u2;
      return { st*std::cos(ph), st*std::sin(ph), ct };
    }
    case GeometryType::CylinderEndOn:
    case GeometryType::CylinderSideOn:
    {
      if( f < 2 )
      {
        const double r = std::sqrt( u1 );
        const double ph = two_pi*u2;
        return { r*std::cos(ph), r*std::sin(ph), (f == 0) ? 1.0 : -1.0 };
      }
      const double ph = two_pi*u1;
      return { std::cos(ph), std::sin(ph), 2.0*u2 - 1.0 };
    }
    case GeometryType::Rectangular:
    {
      std::array<double,3> p;
      const int ax = f/2;
      p[ax] = (f % 2 == 0) ? 1.0 : -1.0;
      p[(ax + 1)%3] = 2.0*u1 - 1.0;
      p[(ax + 2)%3] = 2.0*u2 - 1.0;
      return p;
    }
    case GeometryType::NumGeometryType:
      break;
  }
  assert( 0 );
  return { 0.0, 0.0, 0.0 };
}//unit_surface_point(...)


/** The face a point `p` ON the surface of the solid with `dims` lies on (scalar decision, by which
 normalised coordinate is at its limit), and the outward unit normal there (T: a cylinder side's
 normal turns with the point).  Edges are measure zero and go to whichever face wins the compare. */
template<typename T>
inline int surface_face_at( const GeometryType geometry, const std::array<T,3> &dims,
                            const T p[3], T n[3] )
{
  using namespace std;
  using namespace ceres;

  switch( geometry )
  {
    case GeometryType::Spherical:
    {
      const T r = sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] );
      for( int i = 0; i < 3; ++i )
        n[i] = p[i] / r;
      return 0;
    }
    case GeometryType::CylinderEndOn:
    case GeometryType::CylinderSideOn:
    {
      const T rho = sqrt( p[0]*p[0] + p[1]*p[1] );
      const double cap_ratio = std::fabs( scalar_of(p[2]) ) / scalar_of( dims[1] );
      const double side_ratio = scalar_of( rho ) / scalar_of( dims[0] );
      if( cap_ratio >= side_ratio )
      {
        const bool top = (scalar_of(p[2]) >= 0.0);
        n[0] = T(0.0);
        n[1] = T(0.0);
        n[2] = T( top ? 1.0 : -1.0 );
        return top ? 0 : 1;
      }
      n[0] = p[0] / rho;
      n[1] = p[1] / rho;
      n[2] = T(0.0);
      return 2;
    }
    case GeometryType::Rectangular:
    {
      int ax = 0;
      double best = -1.0;
      for( int i = 0; i < 3; ++i )
      {
        const double ratio = std::fabs( scalar_of(p[i]) ) / scalar_of( dims[i] );
        if( ratio > best )
        {
          best = ratio;
          ax = i;
        }
      }
      const bool pos = (scalar_of(p[ax]) >= 0.0);
      for( int i = 0; i < 3; ++i )
        n[i] = T(0.0);
      n[ax] = T( pos ? 1.0 : -1.0 );
      return 2*ax + (pos ? 0 : 1);
    }
    case GeometryType::NumGeometryType:
      break;
  }
  assert( 0 );
  return 0;
}//surface_face_at(...)


/** Interval [a,b] of the line o - s d inside ONE solid, computed from the line's closest approach
 to the origin (the same re-origining #line_shell_intervals_imp does, for the same precision
 reason) and NOT clamped at s = 0: `a < 0` says the origin is inside.  False when missed or when
 the solid lies entirely behind `o`. */
template<typename T, typename D>
inline bool line_solid_interval_imp( const GeometryType geometry, const std::array<T,3> &dims,
                                     const T o[3], const D d[3], T &a, T &b )
{
  using AType = std::conditional_t<std::is_same_v<D,double>,double,T>;
  AType dd( 0.0 );
  T od( 0.0 );
  for( int i = 0; i < 3; ++i )
  {
    dd += d[i]*d[i];
    od += o[i]*d[i];
  }
  const T s_star = od / dd;
  const T o_near[3] = { o[0] - s_star*d[0], o[1] - s_star*d[1], o[2] - s_star*d[2] };
  if( !line_shell_interval_imp( geometry, dims, o_near, d, a, b ) )
    return false;
  a += s_star;
  b += s_star;
  return (scalar_of(b) > 0.0);
}//line_solid_interval_imp(...)


/** The detector-side line set for one source shell of one fit.
 Built by #build_volumetric_line_cache; owned (shared) by #ShieldingSourceChi2Fcn and referenced by
 every #DistributedSrcCalcT of that source.  The candidate lines and the proposal are immutable
 after construction; the per-dimension traces and the per-energy memos are mutex-guarded.

 A candidate line is a hull point plus an aim point frozen in NORMALISED coordinates (the unit
 solid, or a unit face); at every evaluation the aim point is scaled by the CURRENT source
 dimensions (#line_direction_imp), so the line set deforms continuously with the source and the
 direction, weight and chords are smooth functions of the fitted dimensions.  The crystal kernel
 of each line depends on its direction and is re-traced per distinct set of scalar dimensions
 (#traced); see #TracedLines. */
struct VolumetricLineCache
{
  /** Key: what the set was built for. */
  std::shared_ptr<const ceelo::DetectorResponse> response;
  GeometryType geometry = GeometryType::NumGeometryType;
  size_t material_index = 0;
  std::array<double,3> det_position = { 0.0, 0.0, 0.0 };
  std::array<double,3> det_axis = { 0.0, 0.0, -1.0 };
  double det_azimuth = 0.0;
  int num_lines = 0;
  double pad = 1.5;            //padding factor of the volume component's proposal solid
  double surface_frac = 0.3;   //fraction of lines aimed at the surface(s)
  LineSampleParams sample;     //sequence and replica the lines were drawn from

#if( PERFORM_DEVELOPER_CHECKS )
  /** Diagnostics: chord nodes evaluated, and how many of them PrefactorGrid had to clamp (cos_theta
   below 0, i.e. behind the crystal face plane, or ln(d) outside the grid's range) - both should be
   zero for every supported source; the tests read them. */
  mutable std::atomic<uint64_t> diag_num_nodes{ 0 }, diag_nodes_cos_clamped{ 0 }, diag_nodes_d_clamped{ 0 };
#endif

  /** Crystal (CeeLo) -> assembly rotation, and the reference point in the crystal frame (cm). */
  double M[3][3] = { {1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,1.0} };
  std::array<double,3> ref_c = { 0.0, 0.0, 0.0 };

  /** One candidate line: a hull point and a frozen normalised aim point. */
  struct Candidate
  {
    std::array<double,3> point_c;    //hull point, crystal frame (cm)
    std::array<double,3> x_rel;      //hull point relative to the detector position, assembly frame (PhysicalUnits)
    std::array<double,3> n_a;        //outward hull normal, assembly frame
    double area_weight = 0.0;        //A_face/p_face of the hull face (cm^2)
    std::array<double,3> q_unit;     //aim point in the unit solid (volume) or on the unit surface
    int component = 0;               //0: padded outer volume; 1: outer surface
  };
  std::vector<Candidate> cand;

  /** Mixture weights (sum to 1) and the fixed per-face sampling probabilities of the surface
   component (summing to 1 over #surface_face_count faces), chosen at build from the hint
   dimensions' face areas.  The DENSITY of a line uses the same tables whatever the dimensions, so
   the estimator stays unbiased; only its variance depends on how well they still match.

   There is deliberately NO component on a hollow source's INNER surface, and this is the one place
   the symmetry with the outer surface breaks down.  On the outer surface the 1/|n.w| in the density
   is cancelled by the source chord, which vanishes on the same limb - that cancellation is the
   whole design.  A line grazing the CORE still crosses plenty of source, so nothing cancels there:
   the weight collapses across the core's silhouette, and its derivative is enormous.  Measured on a
   water shell around a steel core, fitting the shell's outer radius (HollowSourceGradientProbe):
   with an inner component the analytic gradient came out -2.10e-4 against a true +6.5e-5 - the
   WRONG SIGN - while the value stayed correct to 0.03%.  Without it the gradient is +5.7e-5, i.e.
   the right sign and 13% low.  Do not add one back. */
  double frac_volume = 1.0, frac_outer = 0.0;
  std::array<double,6> face_prob_outer = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  /** HEMISPHERE component (`Candidate::component == 2`): the direction uniform over the hemisphere
   above the hull point's outward normal, density 1/(2 pi) for every kept line.  It depends on
   nothing the fit moves, so it is trivially smooth in the dimensions, and it puts a FLOOR under the
   mixture density: every line's weight is bounded by 2 pi / frac_hemi times its hull share.  The
   aim-point components have no such floor - a line aimed at a point a distance s from the hull
   point has density ~ s^2, so lines aimed at the near skin of a CONTACT source (the ones that
   matter most) carry weights ~ 1/s^2 - which is the heavy tail behind an effective sample size of
   0.1 N there.  Lines it sends past the source cost nothing but their share of N. */
  double frac_hemi = 0.0;
  /** The hemisphere share the build was ASKED for (the knob), whether or not the trigger fired -
   what #matches compares, since the fired value depends on the hint dims. */
  double hemi_frac_request = 0.0;

  /** The lines traced through the detector at ONE set of scalar source dimensions: the crystal
   segments (for the FEP kernel per energy), optionally the same at two perturbed directions (for
   the kernel's direction gradient, see #kernel_set), and the distance/incidence range the padded
   chords span (for the response flags).  Immutable except the lazily built parts, which are
   guarded by `mutex`. */
  struct TracedLines
  {
    /** Scalar outer dims the lines were aimed with. */
    std::array<double,3> key = { 0.0, 0.0, 0.0 };

    /** Per candidate: leaves the hull toward the source side at these dims (else weight zero). */
    std::vector<char> kept;

    /** One KernelRay per candidate, empty when dropped or when the line misses the active crystal
     (`fep_line_probabilities` then gives 0, and indices stay aligned with `cand`). */
    ceelo::ApertureQuadrature q;

    /** The same lines re-traced at directions w_c + delta*e1 and w_c + delta*e2 (crystal frame,
     e1/e2 from #line_tangent_frame), for a forward difference of the kernel in direction.  Built
     on the first Jet request that needs it. */
    mutable ceelo::ApertureQuadrature q_fd[2];
    mutable bool have_fd = false;

    /** Step of that forward difference (rad).  Not smaller: KernelRay stores segment lengths as
     float, so the kernel carries ~1e-7 of noise and the gradient ~1e-4 at this step. */
    static constexpr double fd_delta = 1.0e-3;

    /** Distance (cm, from the crystal-face origin) and incidence-cosine ranges the padded source
     chords span at these dims (see #PrefactorGrid::worst_flag_between). */
    double chord_d_lo_cm = 0.0, chord_d_hi_cm = 0.0;
    double chord_cos_lo = 0.0, chord_cos_hi = 1.0;

    /** Per-line FEP kernel at one energy, and its forward differences (null until built). */
    struct KernelSet
    {
      std::shared_ptr<const std::vector<double>> k, k1, k2;
    };

    mutable std::mutex mutex;
    mutable std::map<double,KernelSet> kernels;

    /** The kernel at `energy_keV`, with the direction differences when `with_fd` (which requires
     `have_fd`).  Memoized. */
    KernelSet kernel_set( const ceelo::DetectorResponse &resp, const double energy_keV,
                          const bool with_fd ) const
    {
      std::lock_guard<std::mutex> lock( mutex );
      KernelSet &ks = kernels[energy_keV];
      if( !ks.k )
      {
        auto k = std::make_shared<std::vector<double>>();
        resp.fep_line_probabilities( energy_keV, q, *k );
        ks.k = k;
      }
      if( with_fd && !ks.k1 )
      {
        assert( have_fd );
        auto k1 = std::make_shared<std::vector<double>>();
        auto k2 = std::make_shared<std::vector<double>>();
        resp.fep_line_probabilities( energy_keV, q_fd[0], *k1 );
        resp.fep_line_probabilities( energy_keV, q_fd[1], *k2 );
        ks.k1 = k1;
        ks.k2 = k2;
      }
      return ks;
    }//kernel_set(...)
  };//struct TracedLines

  /** Traces kept (most recent first): the optimizer evaluates at a point and its trial step and
   may return to the point, so two entries cover Levenberg-Marquardt's pattern without paying a
   third set's memory (~10 MB per 65k lines, x3 with the direction differences). */
  static constexpr size_t sm_max_traces = 2;

  mutable std::mutex memo_mutex;
  mutable std::deque<std::shared_ptr<const TracedLines>> traces;
  mutable std::shared_ptr<const TracedLines> latest;

  /** Distance range (cm) the prefactor grids currently cover; widened by #ensure_prefactor_range.
   Nodes sit on a fixed logarithmic lattice (#build_prefactor_grid), so widening the range leaves
   every existing node's value unchanged and the integrand continuous across the extension. */
  mutable double grid_d_lo_cm = 0.0, grid_d_hi_cm = 0.0;

  /** Collimated responses only: the shadow gate's quadratures (see #CollimatorGateGrid). */
  mutable std::shared_ptr<const CollimatorGateGrid> gate_grid;
  mutable std::map<double,std::shared_ptr<const PrefactorGrid>> prefactor_by_energy;

  /** The prefactor grid at `energy_keV` (memoized; call #ensure_prefactor_range first). */
  std::shared_ptr<const PrefactorGrid> prefactor( const double energy_keV ) const
  {
    double lo, hi;
    std::shared_ptr<const CollimatorGateGrid> gates;
    {
      std::lock_guard<std::mutex> lock( memo_mutex );
      const auto pos = prefactor_by_energy.find( energy_keV );
      if( pos != end(prefactor_by_energy) )
        return pos->second;
      lo = grid_d_lo_cm;
      hi = grid_d_hi_cm;
      gates = gate_grid;
    }
    std::shared_ptr<const PrefactorGrid> g = build_prefactor_grid( *response, energy_keV, lo, hi, gates.get(),
                                                                   sm_prefactor_grid_num_cos,
                                                                   sm_prefactor_grid_num_phi_quadrant,
                                                                   sm_prefactor_grid_nodes_per_octave );
    std::lock_guard<std::mutex> lock( memo_mutex );
    // Another thread may have widened (and cleared) the range while this grid was building; a grid
    //  built for the narrower range must not be memoised under the wider one.
    if( (lo != grid_d_lo_cm) || (hi != grid_d_hi_cm) )
      return g;
    return prefactor_by_energy.emplace( energy_keV, g ).first->second;
  }//prefactor(...)

  /** Makes the prefactor grids cover [d_lo, d_hi] (cm) with a margin, dropping the memoized grids
   (and rebuilding the collimator gate) when the range has to grow. */
  void ensure_prefactor_range( const double d_lo_cm, const double d_hi_cm ) const
  {
    std::lock_guard<std::mutex> lock( memo_mutex );
    if( (grid_d_hi_cm > 0.0) && (d_lo_cm >= grid_d_lo_cm) && (d_hi_cm <= grid_d_hi_cm) )
      return;
    const double lo = std::max( 0.5*d_lo_cm, 1.0e-3 );
    const double hi = std::max( 2.0*d_hi_cm, 2.0*lo );
    grid_d_lo_cm = (grid_d_hi_cm > 0.0) ? std::min( grid_d_lo_cm, lo ) : lo;
    grid_d_hi_cm = std::max( grid_d_hi_cm, hi );
    prefactor_by_energy.clear();
    if( response->descriptor.collimator )
      gate_grid = build_collimator_gate_grid( *response, grid_d_lo_cm, grid_d_hi_cm );
  }//ensure_prefactor_range(...)

  /** The worst response flag the source's chords can meet at `energy_keV` (near-field floor,
   energy clamping, collimator shadowing), at the most recently evaluated dimensions - what the
   fit's warnings report for a volumetric source, which the point query at the source centre alone
   would miss (e.g. a collimated detector looking at a source that extends into the shadow). */
  ceelo::ResponseFlag worst_flag( const double energy_keV ) const
  {
    std::shared_ptr<const TracedLines> t;
    {
      std::lock_guard<std::mutex> lock( memo_mutex );
      t = latest;
    }
    if( !t )
      return ceelo::ResponseFlag::Ok;
    ensure_prefactor_range( t->chord_d_lo_cm, t->chord_d_hi_cm );
    return prefactor( energy_keV )->worst_flag_between( t->chord_d_lo_cm, t->chord_d_hi_cm,
                                                        t->chord_cos_lo, t->chord_cos_hi );
  }//worst_flag(...)

  /** Whether this cache can serve the given configuration: the detector placement, response, line
   count and proposal knobs must match exactly.  The source DIMENSIONS are not part of the key -
   the aim points follow them (see the struct comment), so one set serves a whole fit. */
  bool matches( const ceelo::DetectorResponse *resp, const GeometryType geom, const size_t mat_index,
                const std::array<double,3> &det_pos, const std::array<double,3> &axis,
                const double azimuth, const int n, const double pad_factor,
                const double surface_fraction, const LineSampleParams &sample_params,
                const double hemi_fraction ) const
  {
    return (response.get() == resp) && (geometry == geom) && (material_index == mat_index)
           && (det_position == det_pos) && (det_axis == axis) && (det_azimuth == azimuth)
           && (num_lines == n) && (pad == pad_factor) && (surface_frac == surface_fraction)
           && (sample == sample_params) && (hemi_frac_request == hemi_fraction);
  }

  /** The lines traced at the scalar source dims `dims_outer`, with the direction differences when
   `want_fd`.  Memoized (#sm_max_traces); under #sm_line_trace_hold the most recent trace is
   returned regardless. */
  std::shared_ptr<const TracedLines> traced( const std::array<double,3> &dims_outer,
                                             const bool want_fd, const bool multithread ) const;

  /** Fills `q_fd` of a trace (idempotent). */
  void trace_direction_differences( const TracedLines &t, const bool multithread ) const;
};//struct VolumetricLineCache


/** Direction of candidate line `j` of `cache` at the given (T-valued, floored) source dims:
 the OUTWARD unit direction `w` from the hull point toward the scaled aim point, the hull cosine
 `cos_n`, the crystal-frame direction `w_c`, and the distance `s_endcap` back along the line from
 the hull point to the endcap-front plane (PhysicalUnits).  False when the line leaves the hull
 elsewhere or points behind the face plane - weight zero, not renormalised (a scalar decision; the
 weight goes to zero continuously as cos_n does).  The scalar trace and the T integration both use
 this, so the traced direction is exactly the scalar part of the integrated one. */
template<typename T>
inline bool line_direction_imp( const VolumetricLineCache &cache, const size_t j,
                                const std::array<T,3> &dims_outer, const T det_pos[3],
                                T w[3], T &cos_n, T w_c[3], T &s_endcap )
{
  using namespace std;
  using namespace ceres;
  const double cm = PhysicalUnits::cm;
  const VolumetricLineCache::Candidate &c = cache.cand[j];

  if( c.component == 2 )
  {
    // Hemisphere component: the frozen direction is (cos = u1, phi = 2 pi u2) about the hull
    //  normal - independent of the dims, hence plain doubles.
    Eigen::Vector3d e1, e2;
    line_tangent_frame( Eigen::Vector3d( c.n_a[0], c.n_a[1], c.n_a[2] ), e1, e2 );
    const double ct = c.q_unit[0];
    const double st = std::sqrt( std::max( 0.0, 1.0 - ct*ct ) );
    const double ph = 2.0*PhysicalUnits::pi*c.q_unit[1];
    for( int i = 0; i < 3; ++i )
      w[i] = T( ct*c.n_a[i] + st*(std::cos(ph)*e1[i] + std::sin(ph)*e2[i]) );
  }else
  {
    // Aim point: the frozen unit-coordinate point scaled by the current dims.
    std::array<T,3> scale;
    switch( c.component )
    {
      case 0:
      {
        const std::array<T,3> padded = { T(cache.pad)*dims_outer[0], T(cache.pad)*dims_outer[1],
                                         T(cache.pad)*dims_outer[2] };
        scale = solid_scale( cache.geometry, padded );
        break;
      }
      default: scale = solid_scale( cache.geometry, dims_outer ); break;
    }

    T wn( 0.0 );
    for( int i = 0; i < 3; ++i )
    {
      w[i] = T(c.q_unit[i])*scale[i] - (det_pos[i] + T(c.x_rel[i]));
      wn += w[i]*w[i];
    }
    if( !(scalar_of(wn) > 0.0) )
      return false;
    wn = sqrt( wn );
    for( int i = 0; i < 3; ++i )
      w[i] /= wn;
  }

  cos_n = w[0]*T(c.n_a[0]) + w[1]*T(c.n_a[1]) + w[2]*T(c.n_a[2]);
  if( !(scalar_of(cos_n) > 0.0) )
    return false;

  for( int i = 0; i < 3; ++i )
    w_c[i] = T(cache.M[0][i])*w[0] + T(cache.M[1][i])*w[1] + T(cache.M[2][i])*w[2];   //M^T w
  if( !(scalar_of(w_c[2]) < 0.0) )
    return false;

  // Distance back along the line from the hull point to the endcap-front plane z = ref_c.z.
  s_endcap = select_max( T( c.point_c[2] - cache.ref_c[2] ) / (-w_c[2]), T(0.0) ) * T(cm);
  return true;
}//line_direction_imp(...)


/** The mixture proposal's density of a line, per unit direction solid angle at its hull point -
 what the line's weight is the reciprocal of.  `o` is the hull point and `d` the photon direction.

 The volume component is the direction density of a point uniform in the padded outer solid,
 (s1^3 - s0^3)/(3 V_p) over the forward chord [s0, s1] through it.  A surface component is the
 density of a point uniform (per face, with the fixed face probabilities) on that surface, summed
 over the forward crossings P of the line with it: p_face |P - o|^2 / (A_face |n.w|).  That
 1/|n.w| is the point: it is the Jacobian between the transverse line measure and the surface
 measure, i.e. exactly the factor by which a line grazing the surface has its chord's dimension
 derivative blow up, so with it in the density the weighted contribution of every line - value and
 derivative - stays bounded (see DIRECTION PROPOSAL in the file comment).

 Every intersection here is computed from the dims the AIM POINTS were scaled by, never from the
 calculator's shell intervals: the density has to describe the distribution the lines were actually
 drawn from, and a shell can be missed (a core collapsed to nothing) while its component is still
 being sampled.  Getting that wrong biases the estimate by the whole weight of the orphaned
 component - measured at +12% for a hollow source whose core reached exactly zero. */
template<typename T>
inline T line_proposal_density_imp( const VolumetricLineCache &cache,
                                    const std::array<T,3> &dims_outer,
                                    const T o[3], const T d[3] )
{
  using namespace std;
  using namespace ceres;
  const GeometryType geometry = cache.geometry;
  const T w[3] = { -d[0], -d[1], -d[2] };

  T p( 0.0 );

  // Volume component.
  if( cache.frac_volume > 0.0 )
  {
    const std::array<T,3> padded = { T(cache.pad)*dims_outer[0], T(cache.pad)*dims_outer[1],
                                     T(cache.pad)*dims_outer[2] };
    T s0, s1;
    if( line_solid_interval_imp( geometry, padded, o, d, s0, s1 ) )
    {
      s0 = select_max( s0, T(0.0) );
      const T s3 = s1*s1*s1 - s0*s0*s0;
      if( scalar_of(s3) > 0.0 )
        p += T(cache.frac_volume) * s3 / (T(3.0)*solid_volume( geometry, padded ));
    }
  }

  // Surface components: every forward crossing of the line with that surface.
  const auto surface_term = [&]( const std::array<T,3> &dims,
                                 const std::array<double,6> &face_prob ) -> T {
    T sum( 0.0 );
    T s_a, s_b;
    if( !line_solid_interval_imp( geometry, dims, o, d, s_a, s_b ) )
      return sum;
    for( const T *s : { &s_a, &s_b } )
    {
      if( !(scalar_of(*s) > 0.0) )
        continue;   //behind the hull point
      const T P[3] = { o[0] - (*s)*d[0], o[1] - (*s)*d[1], o[2] - (*s)*d[2] };
      T n[3];
      const int f = surface_face_at( geometry, dims, P, n );
      const T cos_pn = n[0]*w[0] + n[1]*w[1] + n[2]*w[2];
      const T abs_cos = (scalar_of(cos_pn) < 0.0) ? -cos_pn : cos_pn;
      if( !(scalar_of(abs_cos) > 0.0) )
        continue;
      const T dist2 = (*s)*(*s);   //|P - o| = s: d is a unit vector
      sum += T(face_prob[f]) * dist2 / (surface_face_area( geometry, f, dims ) * abs_cos);
    }
    return sum;
  };

  if( cache.frac_outer > 0.0 )
    p += T(cache.frac_outer) * surface_term( dims_outer, cache.face_prob_outer );

  // Hemisphere component: uniform above the hull normal; the caller only asks about kept lines,
  //  which all lie in that hemisphere.
  if( cache.frac_hemi > 0.0 )
    p += T( cache.frac_hemi / (2.0*PhysicalUnits::pi) );

  return p;
}//line_proposal_density_imp(...)


/** The unit coordinates every candidate line is built from - eight per line:
   [0] hull face, [1],[2] point on the hull face, [3] mixture component, [4] source face,
   [5],[6],[7] aim point in the unit solid / on the unit face.
 Generated by the sequence #LineSampleParams names.  Halton reproduces the original construction
 bit for bit (bases 2,3,5 | 17,19,7,11,13 at index_offset + i); Sobol' assigns its best-behaved
 low dimensions to the aim point and applies a random digital shift per dimension from `seed`, so
 each seed is an independent randomisation of the same (t,s)-net. */
struct LineSampleStream
{
  static constexpr int num_dims = 8;

  LineSampleStream( const LineSampleParams &params, const size_t num_lines )
    : m_params( params )
  {
    if( params.kind != LineSampleParams::Kind::Sobol )
      return;
    // Sobol' coordinates are precomputed (the engine is sequential); 8 x N doubles, transient.
    boost::random::sobol engine( num_dims );
    engine.seed( params.index_offset );   //the point index to start at
    std::mt19937_64 shift_rng( params.seed * 0x9E3779B97F4A7C15ull + 0x2545F4914F6CDD1Dull );
    uint64_t shift[num_dims];
    for( int d = 0; d < num_dims; ++d )
      shift[d] = shift_rng();
    // Sobol' dimension d (0 = best) -> coordinate slot: aim point first, then hull point, then the
    //  three discrete choices.
    static const int slot_of_dim[num_dims] = { 5, 6, 7, 1, 2, 3, 0, 4 };
    m_sobol.resize( num_lines * num_dims );
    const double norm = 1.0 / 18446744073709551616.0;   //2^-64
    for( size_t i = 0; i < num_lines; ++i )
      for( int d = 0; d < num_dims; ++d )
        m_sobol[i*num_dims + slot_of_dim[d]] = static_cast<double>( engine() ^ shift[d] ) * norm;
  }

  void point( const size_t i, double u[num_dims] ) const
  {
    if( m_params.kind == LineSampleParams::Kind::Sobol )
    {
      for( int d = 0; d < num_dims; ++d )
        u[d] = m_sobol[i*num_dims + d];
      return;
    }
    const uint64_t idx = m_params.index_offset + static_cast<uint64_t>(i);
    u[0] = ceelo::halton( idx, 2 );
    u[1] = ceelo::halton( idx, 3 );
    u[2] = ceelo::halton( idx, 5 );
    u[3] = ceelo::halton( idx, 17 );
    u[4] = ceelo::halton( idx, 19 );
    u[5] = ceelo::halton( idx, 7 );
    u[6] = ceelo::halton( idx, 11 );
    u[7] = ceelo::halton( idx, 13 );
  }

private:
  LineSampleParams m_params;
  std::vector<double> m_sobol;
};//struct LineSampleStream


/** Hull points from a #LineSampleStream - the host-side twin of `ceelo::sample_hull_points`, so
 the coordinates can come from any sequence (CeeLo's takes only a Halton index offset).  Same faces
 (front disc/rect + side wall(s); no back face), same projected-area allocation with a 5% floor,
 same per-face point construction, written only against the public Geometry accessors;
 `HullSamplerMatchesCeeLo` pins it bit for bit against CeeLo's for the Halton kind. */
inline void host_sample_hull_points( const ceelo::Geometry &geom, const int n, const Eigen::Vector3d &toward,
                                     const bool have_direction, const LineSampleStream &stream,
                                     std::vector<ceelo::HullPoint> &out )
{
  const double pi = PhysicalUnits::pi;
  struct Face { int kind; double area; Eigen::Vector3d normal; int axis; double sign; double prob; };
  enum { FrontDisc, FrontRect, CylSide, BoxSide };
  std::vector<Face> faces;
  const double L = geom.detector_length();
  if( geom.shape() == ceelo::DetectorShape::Cylinder )
  {
    const double R = geom.detector_radius();
    faces.push_back( { FrontDisc, pi*R*R, { 0.0, 0.0, -1.0 }, 0, 1.0, 0.0 } );
    faces.push_back( { CylSide, 2.0*pi*R*L, { 0.0, 0.0, 0.0 }, 0, 1.0, 0.0 } );
  }else
  {
    const double hx = geom.detector_half_x(), hy = geom.detector_half_y();
    faces.push_back( { FrontRect, 4.0*hx*hy, { 0.0, 0.0, -1.0 }, 0, 1.0, 0.0 } );
    for( const double sgn : { +1.0, -1.0 } )
    {
      faces.push_back( { BoxSide, 2.0*hy*L, { sgn, 0.0, 0.0 }, 0, sgn, 0.0 } );
      faces.push_back( { BoxSide, 2.0*hx*L, { 0.0, sgn, 0.0 }, 1, sgn, 0.0 } );
    }
  }
  double norm = 0.0;
  for( Face &f : faces )
  {
    double proj = f.area;
    if( have_direction )
    {
      if( f.kind == CylSide )
      {
        const double st = std::sqrt( std::max( 0.0, 1.0 - toward.z()*toward.z() ) );
        proj = 2.0*geom.detector_radius()*L*st;
      }else
      {
        proj = f.area * std::max( 0.0, toward.dot( f.normal ) );
      }
    }
    f.prob = std::max( proj, 0.05*f.area );
    norm += f.prob;
  }
  for( Face &f : faces )
    f.prob /= norm;

  out.clear();
  out.reserve( static_cast<size_t>(n) );
  for( int i = 0; i < n; ++i )
  {
    double u[LineSampleStream::num_dims];
    stream.point( static_cast<size_t>(i), u );
    const double u_face = u[0], u1 = u[1], u2 = u[2];
    size_t fi = 0;
    double acc = 0.0;
    for( ; fi + 1 < faces.size(); ++fi )
    {
      acc += faces[fi].prob;
      if( u_face < acc )
        break;
    }
    const Face &f = faces[fi];
    ceelo::HullPoint hp;
    switch( f.kind )
    {
      case FrontDisc:
      {
        const double R = geom.detector_radius();
        const double r = R*std::sqrt( u1 );
        const double ph = 2.0*pi*u2;
        hp.point = Eigen::Vector3d( r*std::cos(ph), r*std::sin(ph), 0.0 );
        hp.normal = Eigen::Vector3d( 0.0, 0.0, -1.0 );
        break;
      }
      case FrontRect:
        hp.point = Eigen::Vector3d( (2.0*u1 - 1.0)*geom.detector_half_x(), (2.0*u2 - 1.0)*geom.detector_half_y(), 0.0 );
        hp.normal = Eigen::Vector3d( 0.0, 0.0, -1.0 );
        break;
      case CylSide:
      {
        const double R = geom.detector_radius();
        const double ph = 2.0*pi*u1;
        const double c = std::cos( ph ), s = std::sin( ph );
        hp.point = Eigen::Vector3d( R*c, R*s, u2*L );
        hp.normal = Eigen::Vector3d( c, s, 0.0 );
        break;
      }
      case BoxSide:
      {
        const double hx = geom.detector_half_x(), hy = geom.detector_half_y();
        if( f.axis == 0 )
        {
          hp.point = Eigen::Vector3d( f.sign*hx, (2.0*u1 - 1.0)*hy, u2*L );
          hp.normal = Eigen::Vector3d( f.sign, 0.0, 0.0 );
        }else
        {
          hp.point = Eigen::Vector3d( (2.0*u1 - 1.0)*hx, f.sign*hy, u2*L );
          hp.normal = Eigen::Vector3d( 0.0, f.sign, 0.0 );
        }
        break;
      }
    }
    hp.area_weight = f.area / f.prob;
    out.push_back( hp );
  }
}//host_sample_hull_points(...)


/** Nearest and farthest distance from the point `c` (assembly frame) to the solid with `dims` -
 the exact point-to-solid distance for the near side, the bounding-sphere over-estimate for the far
 side.  Decides whether a source is "wide" for the hemisphere component. */
inline void padded_solid_distance_range( const GeometryType geometry, const std::array<double,3> &dims,
                                         const std::array<double,3> &c, double &near_dist, double &far_dist )
{
  double cc = 0.0;
  for( int i = 0; i < 3; ++i )
    cc += c[i]*c[i];
  cc = std::sqrt( cc );
  far_dist = cc + solid_bounding_radius( geometry, dims );

  switch( geometry )
  {
    case GeometryType::Spherical:
      near_dist = std::max( 0.0, cc - dims[0] );
      break;
    case GeometryType::CylinderEndOn:
    case GeometryType::CylinderSideOn:
    {
      const double rho = std::hypot( c[0], c[1] );
      const double dr = std::max( 0.0, rho - dims[0] );
      const double dz = std::max( 0.0, std::fabs(c[2]) - dims[1] );
      near_dist = std::hypot( dr, dz );
      break;
    }
    case GeometryType::Rectangular:
    {
      double d2 = 0.0;
      for( int i = 0; i < 3; ++i )
      {
        const double d = std::max( 0.0, std::fabs(c[i]) - dims[i] );
        d2 += d*d;
      }
      near_dist = std::sqrt( d2 );
      break;
    }
    case GeometryType::NumGeometryType:
      near_dist = 0.0;
      break;
  }
}//padded_solid_distance_range(...)


/** Builds the line set for one source shell.  `source_outer_dims` are the
 SCALAR cumulative outer dims of the source shell AT BUILD TIME - a hint only: they set the
 hull-face allocation, the mixture's face probabilities and the sanity trace, while the lines
 themselves follow whatever dims each evaluation carries.  `det_*` is the scalar detector geometry
 (position = the response's reference point in the assembly frame). */
inline std::shared_ptr<const VolumetricLineCache> build_volumetric_line_cache(
                                          std::shared_ptr<const ceelo::DetectorResponse> response,
                                          const GeometryType geometry,
                                          const size_t material_index,
                                          const std::array<double,3> &source_outer_dims,
                                          const std::array<double,3> &det_position,
                                          const std::array<double,3> &det_axis,
                                          const double det_azimuth,
                                          const int num_lines,
                                          const double pad = 1.5,
                                          const double surface_frac = sm_default_volumetric_line_surface_frac,
                                          const LineSampleParams &sample = LineSampleParams(),
                                          const double hemi_frac = sm_volumetric_line_hemi_frac )
{
  using namespace std;
  const double cm = PhysicalUnits::cm;

  if( !response || (num_lines <= 0) )
    throw runtime_error( "build_volumetric_line_cache: invalid inputs" );
  // Every line that crosses the source crosses its surface, and the hemisphere covers everything,
  //  so the volume share may be zero as long as some component remains.
  if( !(pad >= 1.0) || !(surface_frac >= 0.0) || !(hemi_frac >= 0.0) || !(surface_frac + hemi_frac <= 1.0)
      || !(surface_frac + hemi_frac + (1.0 - surface_frac - hemi_frac) > 0.0) )
    throw runtime_error( "build_volumetric_line_cache: invalid proposal knobs" );

  auto cache = make_shared<VolumetricLineCache>();
  cache->response = response;
  cache->geometry = geometry;
  cache->material_index = material_index;
  cache->det_position = det_position;
  cache->det_axis = det_axis;
  cache->det_azimuth = det_azimuth;
  cache->num_lines = num_lines;
  cache->pad = pad;
  cache->surface_frac = surface_frac;
  cache->sample = sample;

  detector_frame_rotation( det_axis.data(), det_azimuth, cache->M );
  // The assembly's detector placement is measured to the DETECTOR FACE (InterSpec's one distance
  //  convention); its image in the crystal frame is the face position, whatever the descriptor's
  //  own reference_point says - see CeeLoUtils::sourcePositionFromFace.
  const Eigen::Vector3d r0 = CeeLoUtils::detectorFacePosition( response->descriptor );
  cache->ref_c = { r0.x(), r0.y(), r0.z() };

  // Hint dims, floored at the same extent ratio the dispatcher floors the source itself at.
  const double det_dist = std::sqrt( det_position[0]*det_position[0] + det_position[1]*det_position[1]
                                     + det_position[2]*det_position[2] );
  const double ext_floor = sm_line_path_extent_ratio_floor * det_dist;
  const int ndims = (geometry == GeometryType::Spherical) ? 1 : ((geometry == GeometryType::Rectangular) ? 3 : 2);
  // All three components are floored, exactly as line_source_integration_imp floors its `dims_o`,
  //  so the sanity trace below is keyed identically to the first evaluation's and gets reused.
  std::array<double,3> hint_outer = source_outer_dims;
  for( int i = 0; i < 3; ++i )
    hint_outer[i] = std::max( std::fabs(hint_outer[i]), ext_floor );
  std::array<double,3> pdims = hint_outer;
  for( int i = 0; i < ndims; ++i )
    pdims[i] *= pad;

  // Mixture: the surface share split between the outer and inner surfaces by area, each surface
  //  sampled per face BY AREA.
  //
  //  Not by the area a face presents to the detector, which is what the hull-point allocation uses
  //  and what an eye on the VALUE would suggest.  MEASURED (LineProposalSurfaceFractionSweep, a
  //  shielded end-on source at 60 keV): a projected-area allocation with a 5% floor sent the
  //  2^18-line GRADIENT error to +3.0 / +4.9 / +7.3 / +10.5% at surface fractions of 0.1 / 0.2 /
  //  0.3 / 0.5, growing with the fraction, while plain area holds it at +0.9 / +1.0 / +0.6 /
  //  -1.3%.  The VALUE is unbiased either way (0.05-0.15%), so it is a derivative-only trap and a
  //  value-motivated allocation walks straight into it.
  //
  //  The reason to expect that, and the reason plain area is the safe default: a dimension's
  //  boundary term lives on whichever face MOVES when that dimension changes - a cylinder's side
  //  wall for its radius, its caps for its half-length - which is unrelated to which face the
  //  detector sees, and area weighting covers every face in proportion to its size.  (The precise
  //  mechanism by which the projected allocation biases the derivative was not pinned down; the
  //  measurement above is the evidence, not the explanation.)
  const int nfaces = surface_face_count( geometry );
  double area_outer = 0.0;
  for( int f = 0; f < nfaces; ++f )
  {
    cache->face_prob_outer[f] = surface_face_area( geometry, f, hint_outer );
    area_outer += cache->face_prob_outer[f];
  }
  for( int f = 0; f < nfaces; ++f )
    cache->face_prob_outer[f] /= area_outer;
  // The hemisphere share fires only for a WIDE source (see sm_volumetric_line_hemi_ratio): the
  //  nearest and farthest points of the padded solid from the detector point, the near distance
  //  floored at the crystal's transverse extent (closer than that the flux no longer falls as
  //  1/d^2 and "wide" stops meaning anything).
  double near_dist = 0.0, far_dist = 0.0;
  padded_solid_distance_range( geometry, pdims, det_position, near_dist, far_dist );
  const double crystal_half = response->transverse_half_extent()*cm;
  const bool wide = (far_dist > sm_volumetric_line_hemi_ratio*std::max( near_dist, crystal_half ));
  const double hemi_used = wide ? hemi_frac : 0.0;
  cache->frac_volume = 1.0 - surface_frac - hemi_used;
  cache->frac_outer = surface_frac;
  cache->frac_hemi = hemi_used;
  cache->hemi_frac_request = hemi_frac;

  const ceelo::Geometry &geom = response->geometry();

  // Assembly -> crystal helpers (double).
  const auto to_crystal_dir = [&]( const double v[3], Eigen::Vector3d &out ) {
    for( int i = 0; i < 3; ++i )
      out[i] = cache->M[0][i]*v[0] + cache->M[1][i]*v[1] + cache->M[2][i]*v[2];   //M^T v
  };
  const auto to_assembly_dir = [&]( const Eigen::Vector3d &v, double out[3] ) {
    for( int i = 0; i < 3; ++i )
      out[i] = cache->M[i][0]*v.x() + cache->M[i][1]*v.y() + cache->M[i][2]*v.z();
  };

  // Source centre (assembly origin) in the crystal frame, for the hull-face allocation.
  const double neg_pos[3] = { -det_position[0], -det_position[1], -det_position[2] };
  Eigen::Vector3d centre_c;
  to_crystal_dir( neg_pos, centre_c );
  centre_c = centre_c/cm + r0;
  const Eigen::Vector3d hull_centre( 0.0, 0.0, 0.5*geom.detector_length() );
  Eigen::Vector3d toward = centre_c - hull_centre;
  const double toward_norm = toward.norm();
  const bool have_dir = (toward_norm > solid_bounding_radius( geometry, pdims )/cm);
  if( toward_norm > 0.0 )
    toward /= toward_norm;

  const LineSampleStream stream( sample, static_cast<size_t>(num_lines) );
  std::vector<ceelo::HullPoint> hull;
  host_sample_hull_points( geom, num_lines, toward, have_dir, stream, hull );

  cache->cand.resize( num_lines );
  for( int i = 0; i < num_lines; ++i )
  {
    const ceelo::HullPoint &hp = hull[i];
    double u[LineSampleStream::num_dims];
    stream.point( static_cast<size_t>(i), u );
    VolumetricLineCache::Candidate &c = cache->cand[i];

    c.point_c = { hp.point.x(), hp.point.y(), hp.point.z() };
    const Eigen::Vector3d xc_rel = (hp.point - r0)*cm;
    to_assembly_dir( xc_rel, c.x_rel.data() );
    to_assembly_dir( hp.normal, c.n_a.data() );
    c.area_weight = hp.area_weight;

    // Component and aim point (see LineSampleStream for which coordinate is which).
    const double u_comp = u[3];
    if( u_comp >= 1.0 - cache->frac_hemi )
    {
      c.component = 2;
      c.q_unit = { u[5], u[6], 0.0 };   //(cos, phi/2pi) about the hull normal
    }
    else if( u_comp < cache->frac_outer )
    {
      c.component = 1;
      const std::array<double,6> &fp = cache->face_prob_outer;
      const double u_face = u[4];
      int f = 0;
      double acc = 0.0;
      for( ; f + 1 < nfaces; ++f )
      {
        acc += fp[f];
        if( u_face < acc )
          break;
      }
      c.q_unit = unit_surface_point( geometry, f, u[5], u[6] );
    }else
    {
      c.component = 0;
      c.q_unit = uniform_point_in_solid( geometry, { 1.0, 1.0, 1.0 }, u[5], u[6], u[7] );
    }
  }//for( loop over lines )

  // Sanity trace at the hint dims: the set must reach the crystal from where the fit starts.
  const std::shared_ptr<const VolumetricLineCache::TracedLines> t = cache->traced( hint_outer, false, true );
  bool any_active = false;
  for( const ceelo::KernelRay &r : t->q.rays )
    any_active = any_active || (r.active_len > 0.0f);
  if( !any_active )
    throw runtime_error( "build_volumetric_line_cache: no line reaches the active crystal" );
  cache->ensure_prefactor_range( t->chord_d_lo_cm, t->chord_d_hi_cm );

  return cache;
}//build_volumetric_line_cache(...)


inline std::shared_ptr<const VolumetricLineCache::TracedLines> VolumetricLineCache::traced(
                                             const std::array<double,3> &dims_outer,
                                             const bool want_fd, const bool multithread ) const
{
  using namespace std;
  const double cm = PhysicalUnits::cm;

  const std::array<double,3> key = dims_outer;
  std::shared_ptr<const TracedLines> hit;
  {
    std::lock_guard<std::mutex> lock( memo_mutex );
    if( sm_line_trace_hold && latest )
      return latest;
    for( const std::shared_ptr<const TracedLines> &t : traces )
    {
      if( t->key == key )
      {
        latest = t;
        hit = t;
        break;
      }
    }
  }
  if( hit )
  {
    // Outside the lock: tracing the difference pair is ~150 ms and would serialise every other
    //  caller of this cache behind it.
    if( want_fd )
      trace_direction_differences( *hit, multithread );
    return hit;
  }

  auto t = std::make_shared<TracedLines>();
  t->key = key;
  const size_t n = cand.size();
  t->kept.assign( n, 0 );
  t->q.n_rays_total = static_cast<int>( n );
  t->q.rays.resize( n );

  const ceelo::Geometry &geom = response->geometry();
  const std::array<double,3> padded = { pad*dims_outer[0], pad*dims_outer[1], pad*dims_outer[2] };
  const double det_pos[3] = { det_position[0], det_position[1], det_position[2] };

  const size_t num_chunks = (multithread && (n > 4096))
                              ? static_cast<size_t>( std::max( 1, SpecUtilsAsync::num_logical_cpu_cores() ) )
                              : size_t(1);
  std::vector<double> d_min( num_chunks, 1.0e300 ), d_max( num_chunks, 0.0 );
  std::vector<double> c_min( num_chunks, 1.0 ), c_max( num_chunks, 0.0 );

  const auto do_chunk = [&]( const size_t chunk )
  {
    const size_t lo = (n*chunk)/num_chunks, hi = (n*(chunk + 1))/num_chunks;
    std::vector<ceelo::PathSegment> scratch;
    for( size_t j = lo; j < hi; ++j )
    {
      double w[3], cos_n, w_c[3], s_endcap;
      if( !line_direction_imp<double>( *this, j, dims_outer, det_pos, w, cos_n, w_c, s_endcap ) )
        continue;
      t->kept[j] = 1;
      const Candidate &c = cand[j];
      const Eigen::Vector3d x( c.point_c[0], c.point_c[1], c.point_c[2] );
      const Eigen::Vector3d wc( w_c[0], w_c[1], w_c[2] );
      trace_detector_line( geom, x, wc, t->q.rays[j], scratch );

      // Distance range of the padded chord from the crystal origin (for the prefactor grids).
      const double x_a[3] = { det_pos[0] + c.x_rel[0], det_pos[1] + c.x_rel[1], det_pos[2] + c.x_rel[2] };
      const double d[3] = { -w[0], -w[1], -w[2] };
      double s0, s1;
      if( !line_solid_interval_imp( geometry, padded, x_a, d, s0, s1 ) )
        continue;
      s0 = std::max( s0, 0.0 );
      for( const double s : { s0, s1 } )
      {
        const Eigen::Vector3d p_c = x + (s/cm)*wc;
        const double dist = p_c.norm();
        d_min[chunk] = std::min( d_min[chunk], dist );
        d_max[chunk] = std::max( d_max[chunk], dist );
        if( dist > 0.0 )
        {
          const double ct = std::max( 0.0, std::min( 1.0, -p_c.z()/dist ) );
          c_min[chunk] = std::min( c_min[chunk], ct );
          c_max[chunk] = std::max( c_max[chunk], ct );
        }
      }
    }//for( lines in chunk )
  };//do_chunk

  if( num_chunks > 1 )
  {
    SpecUtilsAsync::ThreadPool pool;
    for( size_t chunk = 0; chunk < num_chunks; ++chunk )
      pool.post( [&do_chunk,chunk](){ do_chunk( chunk ); } );
    pool.join();
  }else
  {
    do_chunk( 0 );
  }

  t->chord_d_lo_cm = *std::min_element( begin(d_min), end(d_min) );
  t->chord_d_hi_cm = *std::max_element( begin(d_max), end(d_max) );
  t->chord_cos_lo = *std::min_element( begin(c_min), end(c_min) );
  t->chord_cos_hi = std::max( *std::max_element( begin(c_max), end(c_max) ), t->chord_cos_lo );
  if( !(t->chord_d_hi_cm > 0.0) )
  {
    t->chord_d_lo_cm = 1.0e-3;
    t->chord_d_hi_cm = 2.0e-3;
  }

  if( want_fd )
    trace_direction_differences( *t, multithread );

  std::lock_guard<std::mutex> lock( memo_mutex );
  traces.push_front( t );
  while( traces.size() > sm_max_traces )
    traces.pop_back();
  latest = t;
  return t;
}//VolumetricLineCache::traced(...)


inline void VolumetricLineCache::trace_direction_differences( const TracedLines &t, const bool multithread ) const
{
  std::lock_guard<std::mutex> lock( t.mutex );
  if( t.have_fd )
    return;

  const ceelo::Geometry &geom = response->geometry();
  const size_t n = cand.size();
  for( int s = 0; s < 2; ++s )
  {
    t.q_fd[s].n_rays_total = static_cast<int>( n );
    t.q_fd[s].rays.resize( n );
  }

  const std::array<double,3> dims_o = t.key;
  const double det_pos[3] = { det_position[0], det_position[1], det_position[2] };

  const size_t num_chunks = (multithread && (n > 4096))
                              ? static_cast<size_t>( std::max( 1, SpecUtilsAsync::num_logical_cpu_cores() ) )
                              : size_t(1);
  const auto do_chunk = [&]( const size_t chunk )
  {
    const size_t lo = (n*chunk)/num_chunks, hi = (n*(chunk + 1))/num_chunks;
    std::vector<ceelo::PathSegment> scratch;
    for( size_t j = lo; j < hi; ++j )
    {
      if( !t.kept[j] )
        continue;
      const Candidate &c = cand[j];
      const Eigen::Vector3d x( c.point_c[0], c.point_c[1], c.point_c[2] );
      // The direction in double (the ray stores it as float, too coarse against fd_delta).
      double w[3], cos_n, w_c[3], s_endcap;
      line_direction_imp<double>( *this, j, dims_o, det_pos, w, cos_n, w_c, s_endcap );
      const Eigen::Vector3d wc( w_c[0], w_c[1], w_c[2] );
      Eigen::Vector3d e1, e2;
      line_tangent_frame( wc, e1, e2 );
      const Eigen::Vector3d w1 = (wc + TracedLines::fd_delta*e1).normalized();
      const Eigen::Vector3d w2 = (wc + TracedLines::fd_delta*e2).normalized();
      trace_detector_line( geom, x, w1, t.q_fd[0].rays[j], scratch );
      trace_detector_line( geom, x, w2, t.q_fd[1].rays[j], scratch );
    }
  };

  if( num_chunks > 1 )
  {
    SpecUtilsAsync::ThreadPool pool;
    for( size_t chunk = 0; chunk < num_chunks; ++chunk )
      pool.post( [&do_chunk,chunk](){ do_chunk( chunk ); } );
    pool.join();
  }else
  {
    do_chunk( 0 );
  }

  t.have_fd = true;
}//VolumetricLineCache::trace_direction_differences(...)


/** Gauss-Legendre nodes/weights on [0,1]. */
inline void unit_gauss_legendre( const int n, const double *&x, const double *&w )
{
  static const double x2[2] = { 0.5 - 0.5/std::sqrt(3.0), 0.5 + 0.5/std::sqrt(3.0) };
  static const double w2[2] = { 0.5, 0.5 };
  static const double x3[3] = { 0.5 - 0.5*std::sqrt(0.6), 0.5, 0.5 + 0.5*std::sqrt(0.6) };
  static const double w3[3] = { 5.0/18.0, 8.0/18.0, 5.0/18.0 };
  static const double x4[4] = { 0.5 - 0.5*0.861136311594053, 0.5 - 0.5*0.339981043584856,
                                0.5 + 0.5*0.339981043584856, 0.5 + 0.5*0.861136311594053 };
  static const double w4[4] = { 0.5*0.347854845137454, 0.5*0.652145154862546,
                                0.5*0.652145154862546, 0.5*0.347854845137454 };
  switch( n )
  {
    case 2: x = x2; w = w2; return;
    case 3: x = x3; w = w3; return;
    case 4: x = x4; w = w4; return;
  }
  assert( 0 );   //only 2, 3 and 4 points are tabulated
  x = x4;
  w = w4;
}//unit_gauss_legendre(...)


/** Integrates a GROUP of calculators that share geometry (same source shell, dims, detector, line
 cache, normalization and in-situ settings) and differ only in energy-dependent coefficients:
 the per-line chord bookkeeping is done once, the energies innermost.  Fills `integral`,
 `m_num_evals` and `m_est_rel_error` on each.  The lines are summed in #sm_line_error_blocks fixed
 contiguous blocks (run on a pool when `multithread`) and reduced in block order, so the result is
 bit-identical whatever the thread count; the block spread is the error estimate.

 `eff_out` (T = double only): when given, receives one #EffShieldComponents per calculator,
 accumulated in the same pass - every chord node's share of the integral times the areal density,
 AN-weighted areal density, hydrogen areal density and the mu-weighted counterparts of ITS OWN path
 to the detector (source near piece, core when the far piece looks through it, outer shells; no
 air) - the per-line analogue of what `integrate_effective_shielding` accumulates on the element
 path from each element's centre ray.  `c[0]` equals `integral`. */
template<typename T>
void line_source_integration_imp( const std::vector<DistributedSrcCalcT<T>*> &group,
                                  const bool multithread,
                                  std::vector<EffShieldComponents> *eff_out = nullptr )
{
  using namespace std;
  using namespace ceres;

  if( group.empty() )
    return;

  const bool accumulate_eff = (eff_out != nullptr);
  if constexpr( !std::is_same_v<T,double> )
  {
    if( accumulate_eff )
      throw logic_error( "line_source_integration_imp: effective-shielding components are double-only" );
  }

  DistributedSrcCalcT<T> &lead = *group.front();
  const std::shared_ptr<const VolumetricLineCache> cache = lead.m_lineCache;
  if( !cache || !lead.m_effResponse )
    throw logic_error( "line_source_integration_imp: calculator has no line cache/response" );

  for( DistributedSrcCalcT<T> *calc : group )
  {
    calc->finalize_shell_coefficients();
    assert( calc->m_lineCache == cache );
    assert( calc->m_geometry == lead.m_geometry );
    assert( calc->m_materialIndex == lead.m_materialIndex );
    assert( calc->m_shells.size() == lead.m_shells.size() );
    assert( calc->m_isInSituExponential == lead.m_isInSituExponential );
    assert( calc->m_normalizeByVolume == lead.m_normalizeByVolume );
#if( PERFORM_DEVELOPER_CHECKS )
    // The whole group is integrated on the LEAD's chords - the per-line intervals and the cascade
    //  field's per-node walks are computed once from `lead.m_shells` - so the members must really
    //  share their geometry, and differ only in the energy-dependent coefficients.
    for( size_t l = 0; l < lead.m_shells.size(); ++l )
      for( int d = 0; d < 3; ++d )
        assert( scalar_of(calc->m_shells[l].dims[d]) == scalar_of(lead.m_shells[l].dims[d]) );
#endif
  }

  const size_t num_calc = group.size();
  const size_t num_lines = cache->cand.size();
  const size_t m = lead.m_materialIndex;
  const size_t num_shells = lead.m_shells.size();
  const GeometryType geometry = lead.m_geometry;
  const std::vector<typename DistributedSrcCalcT<T>::ShellInfo> &shells = lead.m_shells;
  const double cm = PhysicalUnits::cm;

  const bool in_situ = lead.m_isInSituExponential;
  const double relax = lead.m_inSituRelaxationLength;
  const bool radial_profile = in_situ && ((geometry == GeometryType::CylinderSideOn)
                                          || (geometry == GeometryType::Spherical));
  const bool normalize = lead.m_normalizeByVolume;
  if( in_situ )
    assert( normalize && (m == 0) && (relax > 0.0) );

  // Emission normalisation (T): 1/V for TotalActivity, 1/norm for in-situ, 1 otherwise.
  T rho_scale( 1.0 );
  if( normalize )
  {
    const std::array<T,3> &Do = shells[m].dims;
    const T pi( PhysicalUnits::pi );
    if( in_situ )
    {
      // Per unit EMITTING AREA: the activity is Bq/m^2 of the surface the depth is measured from
      //  (the contract at GammaInteractionCalc::TraceActivityType), so the weight is
      //  A_surf / (depth-integral over the volume).  Written as one expression per geometry with
      //  the area's factors cancelled into the norm, so the vanishing-dimension limits the
      //  volume-normalised forms have are kept (test_ShieldingDimLimit pins them).
      const T L( relax );
      switch( geometry )
      {
        case GeometryType::Spherical:            // 4 pi R^2 / (4 pi R^3 h_sph)
          rho_scale = T(1.0) / (Do[0]*sphere_exp_norm_factor( Do[0]/L ));
          break;
        case GeometryType::CylinderEndOn:        // pi R^2 / (pi R^2 2 L_o g)
          rho_scale = T(1.0) / (T(2.0)*Do[1]*one_minus_exp_neg_over_x( T(2.0)*Do[1]/L ));
          break;
        case GeometryType::CylinderSideOn:       // 2 pi R (2 L_o) / (4 pi L_o R^2 h_side)
          rho_scale = T(1.0) / (Do[0]*cyl_side_exp_norm_factor( Do[0]/L ));
          break;
        case GeometryType::Rectangular:          // (2W)(2H) / (8 W H D_o g)
          rho_scale = T(1.0) / (T(2.0)*Do[2]*one_minus_exp_neg_over_x( T(2.0)*Do[2]/L ));
          break;
        case GeometryType::NumGeometryType:
          assert( 0 );
          break;
      }
    }else if( m == 0 )
    {
      switch( geometry )
      {
        case GeometryType::Spherical:      rho_scale = T(1.0) / (T(4.0/3.0)*pi*Do[0]*Do[0]*Do[0]); break;
        case GeometryType::CylinderEndOn:
        case GeometryType::CylinderSideOn: rho_scale = T(1.0) / (T(2.0)*pi*Do[0]*Do[0]*Do[1]); break;
        case GeometryType::Rectangular:    rho_scale = T(1.0) / (T(8.0)*Do[0]*Do[1]*Do[2]); break;
        case GeometryType::NumGeometryType: assert( 0 ); break;
      }
    }else
    {
      // Annular volume via thickness deltas (no near-cancellation in the derivative lane) - the
      //  same expansions eval_* use.
      const std::array<T,3> &Di = shells[m-1].dims;
      switch( geometry )
      {
        case GeometryType::Spherical:
        {
          const T dR = Do[0] - Di[0];
          rho_scale = T(1.0) / (T(4.0/3.0)*pi*dR*(Do[0]*Do[0] + Do[0]*Di[0] + Di[0]*Di[0]));
          break;
        }
        case GeometryType::CylinderEndOn:
        case GeometryType::CylinderSideOn:
        {
          const T dR = Do[0] - Di[0], dL = Do[1] - Di[1];
          const T vol_factor = dL*(Di[0]*Di[0]) + dR*(T(2.0)*Di[0]*Di[1]) + T(2.0)*Di[0]*dR*dL
                             + (dR*dR)*Di[1] + (dR*dR)*dL;
          rho_scale = T(1.0) / (T(2.0)*pi*vol_factor);
          break;
        }
        case GeometryType::Rectangular:
        {
          const T dW = Do[0] - Di[0], dH = Do[1] - Di[1], dD = Do[2] - Di[2];
          const T vol_factor = dW*Di[1]*Di[2] + dH*Di[0]*Di[2] + dD*Di[0]*Di[1]
                             + dW*dH*Di[2] + dW*dD*Di[1] + dH*dD*Di[0] + dW*dH*dD;
          rho_scale = T(1.0) / (T(8.0)*vol_factor);
          break;
        }
        case GeometryType::NumGeometryType:
          assert( 0 );
          break;
      }
    }
  }//if( normalize )

  // Proposal dims (T): the source shell's outer dims, floored as the dispatcher floors them (the
  //  floor moves the scalar part only, the derivative lane is kept), and the shell inside it for a
  //  hollow source.  The lines are aimed with these, so they follow the fitted dimensions.
  const T det_pos[3] = { lead.m_detector.position[0], lead.m_detector.position[1],
                         lead.m_detector.position[2] };
  double det_dist2 = 0.0;
  for( int i = 0; i < 3; ++i )
    det_dist2 += scalar_of(det_pos[i])*scalar_of(det_pos[i]);
  const double ext_floor = sm_line_path_extent_ratio_floor * std::sqrt( det_dist2 );
  std::array<T,3> dims_o;
  bool dims_have_lane = false;
  for( int i = 0; i < 3; ++i )
  {
    T v = shells[m].dims[i];
    if( scalar_of(v) < 0.0 )
      v = -v;
    if( scalar_of(v) < ext_floor )
      v += T( ext_floor - scalar_of(v) );
    dims_o[i] = v;
    dims_have_lane = dims_have_lane || has_derivative_lane( v );
  }
  // The crystal kernel's direction gradient is only needed when a dimension is being differentiated
  //  (never for T = double, and never under the test hold, which freezes the kernel).
  const bool want_fd = dims_have_lane && !sm_line_trace_hold;
  const std::array<double,3> key_o = { scalar_of(dims_o[0]), scalar_of(dims_o[1]), scalar_of(dims_o[2]) };
  const std::shared_ptr<const VolumetricLineCache::TracedLines> trace
                                              = cache->traced( key_o, want_fd, multithread );
  cache->ensure_prefactor_range( trace->chord_d_lo_cm, trace->chord_d_hi_cm );

  // Per-calculator energy-dependent inputs.
  std::vector<VolumetricLineCache::TracedLines::KernelSet> kernels( num_calc );
  std::vector<std::shared_ptr<const PrefactorGrid>> grids( num_calc );
  for( size_t c = 0; c < num_calc; ++c )
  {
    kernels[c] = trace->kernel_set( *cache->response, group[c]->m_energy, want_fd );
    grids[c] = cache->prefactor( group[c]->m_energy );
  }

  // Cascade-summing corrections, tabulated once per group and interpolated along the chords.
  const std::vector<CascadeFieldT<T>> cascade_fields = build_cascade_fields( group, multithread );
  bool any_cascade = false;
  for( const CascadeFieldT<T> &f : cascade_fields )
    any_cascade = any_cascade || !f.empty();
  const CascadeFieldFrameT<T> cascade_frame( lead );

  const int n_gl = radial_profile ? 4 : std::max( 2, std::min( 4, sm_line_chord_gl_points ) );
  const double *gl_x = nullptr, *gl_w = nullptr;
  unit_gauss_legendre( n_gl, gl_x, gl_w );

  // Depth (T) of a point below the emitting face, per geometry (in-situ only).
  const auto depth_at = [&]( const T p[3] ) -> T {
    const std::array<T,3> &Do = shells[m].dims;
    switch( geometry )
    {
      case GeometryType::CylinderEndOn:  return Do[1] - p[2];
      case GeometryType::Rectangular:    return Do[2] - p[2];
      case GeometryType::CylinderSideOn: return Do[0] - sqrt( p[0]*p[0] + p[1]*p[1] );
      case GeometryType::Spherical:      return Do[0] - sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] );
      case GeometryType::NumGeometryType: break;
    }
    return T(0.0);
  };

  // The lines are split into a FIXED number of contiguous index blocks, whatever the thread count:
  //  each block is summed sequentially in index order and the blocks are reduced in block order,
  //  so the result is bit-reproducible on any machine.  The blocks double as the error estimate
  //  (below): a contiguous block of a low-discrepancy sequence is itself a usable quadrature set,
  //  so the spread of the block estimates measures the quadrature's own scatter.
  const size_t num_chunks = static_cast<size_t>( std::min<size_t>( sm_line_error_blocks, num_lines ) );
  std::vector<std::vector<T>> partial( num_chunks, std::vector<T>( num_calc, T(0.0) ) );
  std::vector<std::vector<EffShieldComponents>> partial_eff( accumulate_eff ? num_chunks : 0,
                                                             std::vector<EffShieldComponents>( num_calc ) );
#if( PERFORM_DEVELOPER_CHECKS )
  std::vector<std::array<uint64_t,3>> partial_diag( num_chunks, std::array<uint64_t,3>{ 0, 0, 0 } );
#endif

  const auto do_chunk = [&]( const size_t chunk )
  {
    const size_t lo = (num_lines*chunk)/num_chunks;
    const size_t hi = (num_lines*(chunk + 1))/num_chunks;
    std::vector<T> &acc = partial[chunk];
    std::vector<EffShieldComponents> * const acc_eff = accumulate_eff ? &partial_eff[chunk] : nullptr;
#if( PERFORM_DEVELOPER_CHECKS )
    std::array<uint64_t,3> &diag = partial_diag[chunk];
#endif

    std::vector<T> a( num_shells ), b( num_shells );
    std::vector<char> crossed( num_shells );
    // Per-line, energy-independent chord bookkeeping.
    struct Piece { T s0, len; };
    std::vector<std::pair<size_t,T>> outer_len;      //(shell, near-segment length) for material shells outside the source
    std::vector<size_t> outer_generic;               //generic shells outside the source
    std::vector<std::pair<size_t,T>> core_len;       //(shell, total length) for material shells inside the source
    std::vector<size_t> core_generic;                //generic shells inside the source, crossed

    for( size_t j = lo; j < hi; ++j )
    {
      if( !trace->kept[j] )
        continue;
      const VolumetricLineCache::Candidate &cand = cache->cand[j];
      T w[3], cos_n, w_c[3], s_endcap;
      if( !line_direction_imp( *cache, j, dims_o, det_pos, w, cos_n, w_c, s_endcap ) )
        continue;   //cannot happen: the trace made the same scalar decision
      const T o[3] = { det_pos[0] + T(cand.x_rel[0]), det_pos[1] + T(cand.x_rel[1]),
                       det_pos[2] + T(cand.x_rel[2]) };
      const T d[3] = { -w[0], -w[1], -w[2] };

      line_shell_intervals_imp( geometry, shells, o, d, a, b, crossed );
      if( !crossed[m] )
        continue;

      // This line's weight (T): hull area share times its cosine over the mixture density, per
      //  line, over 4 pi (see the file comment's UNITS paragraph).
      const T p_line = line_proposal_density_imp( *cache, dims_o, o, d );
      if( !(scalar_of(p_line) > 0.0) )
        continue;
      const T w_line = T(cand.area_weight) * cos_n
                       / (p_line * T(4.0*PhysicalUnits::pi*static_cast<double>(num_lines)))
                       * T(cm*cm);

      // The direction's derivative lanes projected on the tangent basis the kernel's forward
      //  differences were traced along (crystal frame; zero scalar part, since the trace was made
      //  at exactly these scalar dims).
      T proj_fd[2] = { T(0.0), T(0.0) };
      if( want_fd )
      {
        Eigen::Vector3d e1, e2;
        line_tangent_frame( Eigen::Vector3d( scalar_of(w_c[0]), scalar_of(w_c[1]), scalar_of(w_c[2]) ),
                            e1, e2 );
        for( int i = 0; i < 3; ++i )
        {
          const T dw = w_c[i] - T( scalar_of(w_c[i]) );
          proj_fd[0] += dw * T(e1[i]);
          proj_fd[1] += dw * T(e2[i]);
        }
      }

      // Source pieces.
      Piece pieces[2];
      int num_pieces = 0;
      const bool has_core = (m > 0) && crossed[m-1];
      if( has_core )
      {
        pieces[0] = { a[m], a[m-1] - a[m] };
        pieces[1] = { b[m-1], b[m] - b[m-1] };
        num_pieces = 2;
      }else
      {
        pieces[0] = { a[m], b[m] - a[m] };
        num_pieces = 1;
      }

      // Outer shells (near segments only) and the air path.
      outer_len.clear();
      outer_generic.clear();
      for( size_t l = m + 1; l < num_shells; ++l )
      {
        if( shells[l].type == ShellType::Generic )
          outer_generic.push_back( l );
        else
          outer_len.emplace_back( l, a[l-1] - a[l] );
      }
      const T air_len = select_max( a[num_shells-1] - s_endcap, T(0.0) );

      // Core (only the far piece looks through it): every inner shell's full chord.
      core_len.clear();
      core_generic.clear();
      if( has_core )
      {
        for( size_t l = 0; l < m; ++l )
        {
          if( !crossed[l] )
            continue;
          if( shells[l].type == ShellType::Generic )
          {
            core_generic.push_back( l );
            continue;
          }
          const bool inner_crossed = (l > 0) && crossed[l-1];
          const T len = inner_crossed ? ((a[l-1] - a[l]) + (b[l] - b[l-1])) : (b[l] - a[l]);
          core_len.emplace_back( l, len );
        }
      }

      // Effective-shielding components: the density-weighted path outside the source (and through
      //  the core, for the far piece) is the same for every calculator; the mu-weighted one is per
      //  calculator, below.
      double ad_out = 0.0, an_ad_out = 0.0, ad_h_out = 0.0, ad_core = 0.0, an_ad_core = 0.0, ad_h_core = 0.0;
      if constexpr( std::is_same_v<T,double> )
      {
        if( acc_eff )
        {
          for( const std::pair<size_t,T> &ol : outer_len )
          {
            const double ad = shells[ol.first].density * ol.second;
            ad_out += ad;
            an_ad_out += shells[ol.first].effective_an * ad;
            ad_h_out += shells[ol.first].hydrogen_mass_frac * ad;
          }
          for( const size_t l : outer_generic )
          {
            ad_out += shells[l].areal_density;
            an_ad_out += shells[l].effective_an * shells[l].areal_density;
          }
          for( const std::pair<size_t,T> &cl : core_len )
          {
            const double ad = shells[cl.first].density * cl.second;
            ad_core += ad;
            an_ad_core += shells[cl.first].effective_an * ad;
            ad_h_core += shells[cl.first].hydrogen_mass_frac * ad;
          }
          for( const size_t l : core_generic )
          {
            ad_core += shells[l].areal_density;
            an_ad_core += shells[l].effective_an * shells[l].areal_density;
          }
        }
      }

      for( size_t c = 0; c < num_calc; ++c )
      {
        const double k = (*kernels[c].k)[j];
        T k_T( k );
        if( want_fd )
        {
          // The kernel's direction gradient (forward difference) chained onto the direction's
          //  lanes: the pathwise derivative of a line that moves with the source.  This term is
          //  NOT zero in the continuum limit - it is of order the relative variation of k across
          //  the set times dI/d(dim), i.e. tens of percent of the dimension gradient for a contact
          //  source at high energy - see the file comment.
          const double g1 = ((*kernels[c].k1)[j] - k) / VolumetricLineCache::TracedLines::fd_delta;
          const double g2 = ((*kernels[c].k2)[j] - k) / VolumetricLineCache::TracedLines::fd_delta;
          k_T += proj_fd[0]*T(g1) + proj_fd[1]*T(g2);
        }
        if( !(k > 0.0) && !has_derivative_lane( k_T ) )
          continue;
        const DistributedSrcCalcT<T> &calc = *group[c];
        const std::vector<typename DistributedSrcCalcT<T>::ShellInfo> &cs = calc.m_shells;
        const T mu_src = cs[m].fep_trans_len_coef;

        // Per-calculator (TOTAL-mu weighted) effective-shielding sums, see record_path.
        double mud_out = 0.0, an_mud_out = 0.0, mud_core = 0.0, an_mud_core = 0.0;
        double rho_src = 0.0, an_src = 0.0, h_src = 0.0, mu_tot_src = 0.0;
        if constexpr( std::is_same_v<T,double> )
        {
          if( acc_eff )
          {
            rho_src = cs[m].density;
            an_src = cs[m].effective_an;
            h_src = cs[m].hydrogen_mass_frac;
            mu_tot_src = cs[m].trans_len_coef;
            for( const std::pair<size_t,T> &ol : outer_len )
            {
              const double mud = cs[ol.first].trans_len_coef * ol.second;
              mud_out += mud;
              an_mud_out += cs[ol.first].effective_an * mud;
            }
            for( const size_t l : outer_generic )
            {
              mud_out += cs[l].trans_len_coef;
              an_mud_out += cs[l].effective_an * cs[l].trans_len_coef;
            }
            for( const std::pair<size_t,T> &cl : core_len )
            {
              const double mud = cs[cl.first].trans_len_coef * cl.second;
              mud_core += mud;
              an_mud_core += cs[cl.first].effective_an * mud;
            }
            for( const size_t l : core_generic )
            {
              mud_core += cs[l].trans_len_coef;
              an_mud_core += cs[l].effective_an * cs[l].trans_len_coef;
            }
          }
        }

        // Transmission through everything outside the source, common to both pieces.
        T tau_out( 0.0 );
        for( const std::pair<size_t,T> &ol : outer_len )
          tau_out += cs[ol.first].fep_trans_len_coef * ol.second;
        for( const size_t l : outer_generic )
          tau_out += cs[l].fep_trans_len_coef;
        if( calc.m_attenuateForAir )
          tau_out += T(calc.m_airTransLenCoef) * air_len;

        T line_sum( 0.0 );
        for( int pc = 0; pc < num_pieces; ++pc )
        {
          const Piece &piece = pieces[pc];
          const T &L = piece.len;
          if( scalar_of(L) <= 0.0 )
            continue;

          T tau_beyond = tau_out;
          if( pc == 1 )
          {
            // The far piece looks through the core and the near source piece.
            for( const std::pair<size_t,T> &cl : core_len )
              tau_beyond += cs[cl.first].fep_trans_len_coef * cl.second;
            for( const size_t l : core_generic )
              tau_beyond += cs[l].fep_trans_len_coef;
            tau_beyond += mu_src * pieces[0].len;
          }

          // Effective shielding: everything beyond this piece's own near path.
          double ad_b = 0.0, an_ad_b = 0.0, ad_h_b = 0.0, mud_b = 0.0, an_mud_b = 0.0;
          if constexpr( std::is_same_v<T,double> )
          {
            if( acc_eff )
            {
              ad_b = ad_out;
              an_ad_b = an_ad_out;
              ad_h_b = ad_h_out;
              mud_b = mud_out;
              an_mud_b = an_mud_out;
              if( pc == 1 )
              {
                const double near_src = pieces[0].len;
                ad_b += ad_core + rho_src*near_src;
                an_ad_b += an_ad_core + an_src*rho_src*near_src;
                ad_h_b += ad_h_core + h_src*rho_src*near_src;
                mud_b += mud_core + mu_tot_src*near_src;
                an_mud_b += an_mud_core + an_src*mu_tot_src*near_src;
              }
            }
          }

          // Sub-pieces: an in-situ profile whose depth is NOT linear along the line (the radial
          //  ones) is carried by the quadrature rather than the exponent, so the chord is cut into
          //  enough pieces that the profile varies by at most ~one relaxation length across each.
          int nsub = 1;
          if( radial_profile )
            nsub = std::max( 1, std::min( 8, static_cast<int>( std::ceil( scalar_of(L)/relax ) ) ) );
          // The effective-shielding pass (post-fit, once) takes the prefactor as constant across a
          //  sub-piece when it weights the path length, so it cuts the chord finely enough for that
          //  to hold at contact (P varies on the cm scale there).
          if( accumulate_eff )
            nsub = std::max( nsub, std::min( 32, static_cast<int>( std::ceil( scalar_of(L)/(0.25*cm) ) ) ) );
          const T sub_len = L / T(static_cast<double>(nsub));

          for( int sub = 0; sub < nsub; ++sub )
          {
          // Sub-piece [s_lo, s_lo + sub_len]; the whole piece when nsub == 1.
          const T s_lo = piece.s0 + T(static_cast<double>(sub))*sub_len;
          const T &Lp = sub_len;

          // Emission at the near end, and (for a depth that IS linear along the line) its rate.
          const T p0[3] = { o[0] - s_lo*d[0], o[1] - s_lo*d[1], o[2] - s_lo*d[2] };
          T rho0 = rho_scale;
          T rate0( 0.0 );
          T depth0( 0.0 );
          if( in_situ )
          {
            depth0 = depth_at( p0 );
            rho0 = rho_scale * exp( -depth0/T(relax) );
            // End-on and rectangular: depth = const + s*d_z, exactly linear, so the profile folds
            //  into the analytic exponent.  Radial profiles (side-on cylinder, sphere) are not
            //  linear in s and would need a signed rate that can be negative and large - which
            //  turns the analytic factor into exp(+big) - so they keep rate0 = 0 and are carried
            //  by the sub-piece quadrature below instead.
            if( !radial_profile )
              rate0 = d[2] / T(relax);
          }//if( in_situ )

          const T mu_eff = mu_src + rate0;
          const T x = mu_eff * Lp;
          const T g = one_minus_exp_neg_over_x( x );
          const T y1 = exp( -x );

          // Everything nearer the detector than this sub-piece attenuates it as well.
          const T tau_sub = tau_beyond + mu_src*(s_lo - piece.s0);
          const T pref = exp( -tau_sub ) * rho0 * Lp * g;

          // Average of the smooth remainder (prefactor P, radial-profile residual) over y in [y1,1].
          T mean( 0.0 );
          for( int n = 0; n < n_gl; ++n )
          {
            const T y = y1 + (T(1.0) - y1)*T(gl_x[n]);
            const T ds = (std::fabs(scalar_of(x)) > 1.0e-6) ? (-log( y )/mu_eff) : (Lp*T(1.0 - gl_x[n]));
            const T s = s_lo + ds;
            const T p[3] = { o[0] - s*d[0], o[1] - s*d[1], o[2] - s*d[2] };

            // Crystal-frame coordinates of p: x_c + (s/cm) w_c, i.e. r0 + M^T (p - det.position)/cm.
            T pc[3];
            for( int i = 0; i < 3; ++i )
            {
              const T rel0 = p[0] - lead.m_detector.position[0];
              const T rel1 = p[1] - lead.m_detector.position[1];
              const T rel2 = p[2] - lead.m_detector.position[2];
              pc[i] = (cache->M[0][i]*rel0 + cache->M[1][i]*rel1 + cache->M[2][i]*rel2)/cm + T(cache->ref_c[i]);
            }
            // Once the source surface reaches the detector face, line_shell_intervals_imp clamps the
            //  near chord end to 0 (L401-402) and emission nodes ride onto the crystal, sending
            //  dist->0.  At the crystal origin the incidence is genuinely undefined: -pc[2]/dist is
            //  0/0 (a NaN VALUE that poisons every residual) and its lane, log(dist) and atan2's lane
            //  all diverge.  Pin the distance to the grid's own inner node d_lo ONLY when it falls
            //  below it (a conditional floor, not an additive term): healthy nodes keep their exact
            //  distance and full derivative lane, a degenerate node gets the constant d_lo -
            //  reproducing the grid's flat-region boundary value and zeroing only that node's lane,
            //  never evaluating sqrt on the degenerate value.  Mirrors the cache-build guard's
            //  `if(dist>0)` (~L2456) and eff_response_factor (GammaInteractionCalc_imp.hpp:1517-1523).
            const double d_lo = std::exp( grids[c]->ln_d.front() );  // grid min distance, cm (>0)
            const T dist2 = pc[0]*pc[0] + pc[1]*pc[1] + pc[2]*pc[2];
            const T dist = (scalar_of(dist2) > d_lo*d_lo) ? sqrt( dist2 ) : T(d_lo);
            // dist >= d_lo > 0 and dist >= |pc[2]|, so cos_t is finite and in [-1,1] with its full
            //  derivative lane; locate clamps any boundary/last-ULP case to the grid's cos axis.
            const T cos_t = -pc[2]/dist;
            T phi( 0.0 );
            if( grids[c]->phi_deg.size() > 1 )
            {
              // Quadrant symmetry: fold into [0,90] degrees.  atan2 is singular (0/0 value AND lane)
              //  only on the crystal axis where the transverse radius is exactly zero; guard on that
              //  alone so every off-axis node keeps the base atan2 value and its full lane.  (A node
              //  with a tiny-but-nonzero transverse radius has a large-but-finite lane, which Ceres
              //  handles and the base code already produced.)
              const T ax = (scalar_of(pc[0]) < 0.0) ? -pc[0] : pc[0];
              const T ay = (scalar_of(pc[1]) < 0.0) ? -pc[1] : pc[1];
              if( (scalar_of(ax)*scalar_of(ax) + scalar_of(ay)*scalar_of(ay)) > 0.0 )
                phi = atan2( ay, ax ) * T(180.0/PhysicalUnits::pi);
            }
            const T ln_dist = log( dist );  // dist >= d_lo > 0 => finite in value and every lane
            T val = grids[c]->eval( ln_dist, cos_t, phi );
#if( PERFORM_DEVELOPER_CHECKS )
            // Every quantity feeding the grid, and its result, must be finite in value AND every lane
            //  (this is where cos_t = -pc[2]/dist = 0/0 used to poison the residual and Jacobian).
            assert( all_finite_lanes(dist) && (scalar_of(dist) > 0.0) );
            assert( all_finite_lanes(cos_t) );
            assert( all_finite_lanes(ln_dist) );
            assert( all_finite_lanes(phi) );
            assert( all_finite_lanes(val) );
            diag[0] += 1;
            if( scalar_of(cos_t) < 0.0 )
              diag[1] += 1;
            if( (scalar_of(dist) < std::exp(grids[c]->ln_d.front())) || (scalar_of(dist) > std::exp(grids[c]->ln_d.back())) )
              diag[2] += 1;
#endif
            if constexpr( std::is_same_v<T,double> )
            {
              if( sm_prefactor_direct_eval )
              {
                // TEST HOOK: the response's prefactor at the node itself, no grid (see the flag).
                assert( !cache->response->descriptor.collimator );
                static const ceelo::ApertureQuadrature no_quadrature;
                const Eigen::Vector3d pos_c( pc[0], pc[1], pc[2] );
                val = cache->response->fep_prefactor( calc.m_energy, pos_c, no_quadrature ).value;
              }
            }else
            {
              if( sm_prefactor_direct_eval )
                throw logic_error( "sm_prefactor_direct_eval is double-only" );
            }
            if( radial_profile )
              val *= exp( -(depth_at( p ) - depth0)/T(relax) );
            if( any_cascade && !cascade_fields[c].empty() )
            {
              T u[3];
              cascade_frame.coords( p, u );
              val *= cascade_fields[c].eval( u );
            }
            mean += T(gl_w[n]) * val;
          }//for( GL nodes )

          line_sum += pref * mean;

          if constexpr( std::is_same_v<T,double> )
          {
            if( acc_eff )
            {
              // This sub-piece's share of the integral, weighting the path from its emission points
              //  to the detector: everything beyond the piece, the near source piece up to s_lo (the
              //  same split as tau_sub), and the contribution-weighted mean of the remaining depth
              //  ds.  That last one is the first moment of exp(-mu_eff ds) over [0, Lp] divided by
              //  its zeroth, M1/M0 = (1 - (1+x)e^-x) / (mu_eff (1 - e^-x)), taken ANALYTICALLY: in
              //  the y = exp(-mu_eff ds) substitution ds is -log(y)/mu_eff, and the two-point rule
              //  that is exact for the smooth remainder under-integrates log(y) by 10% once the
              //  chord is optically thick (measured: <AD> 10% low on a far steel source at 60 keV).
              //  P is taken constant across the sub-piece for this weighting - a second-order
              //  covariance the value itself does not have to make.
              const double ws = w_line * k * pref * mean;
              const double xs = x, mu = mu_eff;
              double mean_ds;
              if( std::fabs(xs) < 1.0e-4 )
                mean_ds = Lp * (0.5 - xs/12.0);
              else
                mean_ds = (1.0 - (1.0 + xs)*std::exp(-xs)) / (mu*(1.0 - std::exp(-xs)));
              const double src_path = (s_lo - piece.s0) + mean_ds;
              EffShieldComponents &ec = (*acc_eff)[c];
              ec.c[0] += ws;
              ec.c[1] += ws * (ad_b + rho_src*src_path);
              ec.c[2] += ws * (an_ad_b + an_src*rho_src*src_path);
              ec.c[3] += ws * (ad_h_b + h_src*rho_src*src_path);
              ec.c[4] += ws * (mud_b + mu_tot_src*src_path);
              ec.c[5] += ws * (an_mud_b + an_src*mu_tot_src*src_path);
            }
          }
          }//for( sub-pieces )
        }//for( pieces )

#if( PERFORM_DEVELOPER_CHECKS )
        assert( all_finite_lanes(line_sum) );
#endif
        acc[c] += w_line * k_T * line_sum;
      }//for( calculators )
    }//for( lines in chunk )
  };//do_chunk

  if( multithread && (num_chunks > 1) && (num_lines > 4096) )
  {
    SpecUtilsAsync::ThreadPool pool;
    for( size_t chunk = 0; chunk < num_chunks; ++chunk )
      pool.post( [&do_chunk,chunk](){ do_chunk( chunk ); } );
    pool.join();
  }else
  {
    for( size_t chunk = 0; chunk < num_chunks; ++chunk )
      do_chunk( chunk );
  }

  for( size_t c = 0; c < num_calc; ++c )
  {
    T total( 0.0 );
    for( size_t chunk = 0; chunk < num_chunks; ++chunk )
      total += partial[chunk][c];
#if( PERFORM_DEVELOPER_CHECKS )
    assert( all_finite_lanes(total) );
#endif
    group[c]->integral = total;
    group[c]->m_num_evals = num_lines;

    // Error estimate from the block partial sums at TWO scales.  Each contiguous block of the
    //  sequence, scaled up by the block count, is an unbiased estimate of the whole, so the scatter
    //  of the B block estimates measures the error of a set N/B lines long, and the scatter of the
    //  S super-block estimates (S groups of B/S consecutive blocks) that of a set N/S long.  A
    //  low-discrepancy set's error falls as n^-alpha with alpha between 1/2 (plain Monte Carlo) and
    //  ~1, and the two scales measure alpha: s_B/s_S = (B/S)^alpha.  The super-block scatter is then
    //  extrapolated to the full set, err(N) = s_S * S^-alpha.  With alpha = 1/2 this is exactly the
    //  plain-MC standard error s_B/sqrt(B); assuming 1/2 outright over-read by 5-10x on a wide disk
    //  (measured: replica rms 0.1-0.26% against 1.3-1.5% estimated).
    //  Contiguous blocks, not an even/odd split: the hull face is chosen by the base-2 Halton digit,
    //  so even and odd indices are (mostly) different hull faces and their difference measured the
    //  face split rather than the error.  `LineErrorEstimateCalibration` pins this against the
    //  spread of independent replicas (`LineSampleParams::index_offset`).
    const double tot = scalar_of( total );
    const size_t B = num_chunks;
    const size_t S = std::min<size_t>( 4, B );
    double ss_b = 0.0, ss_s = 0.0;
    std::vector<double> sup( S, 0.0 );
    for( size_t b = 0; b < B; ++b )
    {
      const double pb = scalar_of( partial[b][c] );
      const double block_est = static_cast<double>(B) * pb;
      ss_b += (block_est - tot)*(block_est - tot);
      sup[(b*S)/B] += pb;
    }
    for( size_t s = 0; s < S; ++s )
    {
      const double sup_est = static_cast<double>(S) * sup[s];
      ss_s += (sup_est - tot)*(sup_est - tot);
    }
    double std_err = 0.0;
    if( (B > 1) && (S > 1) && (std::fabs(tot) > 0.0) )
    {
      const double s_b = std::sqrt( ss_b / static_cast<double>(B - 1) );
      const double s_s = std::sqrt( ss_s / static_cast<double>(S - 1) );
      double alpha = 0.5;
      if( (s_b > 0.0) && (s_s > 0.0) && (B > S) )
        alpha = std::max( 0.5, std::min( 1.0, std::log( s_b/s_s ) / std::log( static_cast<double>(B)/static_cast<double>(S) ) ) );
      std_err = s_s * std::pow( static_cast<double>(S), -alpha );
    }
    group[c]->m_est_rel_error = (std::fabs(tot) > 0.0) ? (std_err / std::fabs(tot)) : 0.0;
  }

#if( PERFORM_DEVELOPER_CHECKS )
  for( size_t chunk = 0; chunk < num_chunks; ++chunk )
  {
    cache->diag_num_nodes += partial_diag[chunk][0];
    cache->diag_nodes_cos_clamped += partial_diag[chunk][1];
    cache->diag_nodes_d_clamped += partial_diag[chunk][2];
  }
#endif

  if( accumulate_eff )
  {
    eff_out->assign( num_calc, EffShieldComponents() );
    for( size_t c = 0; c < num_calc; ++c )
      for( size_t chunk = 0; chunk < num_chunks; ++chunk )
        (*eff_out)[c] += partial_eff[chunk][c];
  }
}//line_source_integration_imp(...)


/** Smallest scalar extent of a calculator's source shell (its thickness for a hollow shell). */
template<typename T>
double smallest_source_extent( const DistributedSrcCalcT<T> &calc )
{
  const size_t m = calc.m_materialIndex;
  const int ndims = (calc.m_geometry == GeometryType::Spherical) ? 1
                    : ((calc.m_geometry == GeometryType::Rectangular) ? 3 : 2);
  double smallest = 1.0e300;
  for( int i = 0; i < ndims; ++i )
  {
    double ext = std::fabs( scalar_of(calc.m_shells[m].dims[i]) );
    if( m > 0 )
      ext -= std::fabs( scalar_of(calc.m_shells[m-1].dims[i]) );
    smallest = std::min( smallest, ext );
  }
  return smallest;
}//smallest_source_extent(...)


/** Raises every source-shell extent of `calc` that is below `ext_floor` to it: the SCALAR part of
 the dimension moves, its derivative lane does not, and the same shift is applied to every shell
 outside the source so the nesting order and the outer thicknesses are unchanged.  Returns whether
 anything was floored.  See #sm_line_path_extent_ratio_floor. */
template<typename T>
bool floor_source_extent( DistributedSrcCalcT<T> &calc, const double ext_floor )
{
  const size_t m = calc.m_materialIndex;
  const int ndims = (calc.m_geometry == GeometryType::Spherical) ? 1
                    : ((calc.m_geometry == GeometryType::Rectangular) ? 3 : 2);
  bool floored = false;
  for( int i = 0; i < ndims; ++i )
  {
    const double outer = std::fabs( scalar_of(calc.m_shells[m].dims[i]) );
    const double inner = (m > 0) ? std::fabs( scalar_of(calc.m_shells[m-1].dims[i]) ) : 0.0;
    const double ext = outer - inner;
    if( ext >= ext_floor )
      continue;
    const T delta( ext_floor - ext );
    for( size_t l = m; l < calc.m_shells.size(); ++l )
      calc.m_shells[l].dims[i] += delta;
    floored = true;
  }
  return floored;
}//floor_source_extent(...)


/** Whether the line path can serve this calculator at all. */
template<typename T>
bool line_path_applicable( const DistributedSrcCalcT<T> &calc )
{
  if( sm_volumetric_integrator_override == VolumetricIntegrator::Element )
    return false;
  if( !calc.m_effResponse || !calc.m_lineCache )
    return false;
  return true;
}//line_path_applicable(...)


/** How #integrate_volumetric_calculators and #integrate_effective_shielding_all split a set of
 calculators: the element-path ones, the line groups (calculators sharing a line cache, in-situ
 settings and normalisation - the per-line chord bookkeeping is shared within a group), and the
 floored copies (#floor_source_extent) standing in for calculators whose source extent is below
 #sm_line_path_extent_ratio_floor.  Group members and `element_only` point at the calculators
 themselves, or at a floored copy; `floored` says which original each copy stands for. */
template<typename T>
struct VolumetricPartitionT
{
  std::vector<DistributedSrcCalcT<T>*> element_only;
  std::vector<std::vector<DistributedSrcCalcT<T>*>> line_groups;
  std::vector<std::pair<DistributedSrcCalcT<T>*,std::unique_ptr<DistributedSrcCalcT<T>>>> floored;
};//struct VolumetricPartitionT


template<typename T>
VolumetricPartitionT<T> partition_volumetric_calculators(
                        const std::vector<std::unique_ptr<DistributedSrcCalcT<T>>> &calculators )
{
  using namespace std;

  struct GroupKey
  {
    const VolumetricLineCache *cache;
    bool in_situ;
    double relax;
    bool normalize;
    bool operator<( const GroupKey &rhs ) const
    {
      return std::tie( cache, in_situ, relax, normalize )
             < std::tie( rhs.cache, rhs.in_situ, rhs.relax, rhs.normalize );
    }
  };

  VolumetricPartitionT<T> part;
  map<GroupKey,size_t> group_index;

  for( const unique_ptr<DistributedSrcCalcT<T>> &calc : calculators )
  {
    if( !line_path_applicable( *calc ) )
    {
      if( sm_volumetric_integrator_override == VolumetricIntegrator::Line )
        throw runtime_error( "integrate_volumetric_calculators: line path forced but not applicable" );
      part.element_only.push_back( calc.get() );
      continue;
    }

    DistributedSrcCalcT<T> *target = calc.get();

    // Vanishing extent: integrate a copy floored at the ratio the chords are still accurate at.
    double dist = 0.0;
    for( int i = 0; i < 3; ++i )
    {
      const double c = scalar_of( calc->m_detector.position[i] );
      dist += c*c;
    }
    const double ext_floor = sm_line_path_extent_ratio_floor * std::sqrt( dist );
    if( smallest_source_extent( *calc ) < ext_floor )
    {
      unique_ptr<DistributedSrcCalcT<T>> copy = make_unique<DistributedSrcCalcT<T>>( *calc );
      floor_source_extent( *copy, ext_floor );
      target = copy.get();
      part.floored.emplace_back( calc.get(), std::move(copy) );
    }

    const GroupKey key{ target->m_lineCache.get(), target->m_isInSituExponential,
                        target->m_inSituRelaxationLength, target->m_normalizeByVolume };
    const auto pos = group_index.find( key );
    if( pos == end(group_index) )
    {
      group_index[key] = part.line_groups.size();
      part.line_groups.push_back( { target } );
    }else
    {
      part.line_groups[pos->second].push_back( target );
    }
  }//for( calculators )

  return part;
}//partition_volumetric_calculators(...)


/** Runs the element path on each calculator (one task each when multithreaded), rethrowing the
 first failure after the pool drains. */
template<typename T>
void integrate_element_calculators( const std::vector<DistributedSrcCalcT<T>*> &element_only,
                                    const bool multithread )
{
  if( element_only.empty() )
    return;

  if( multithread && (element_only.size() > 1) )
  {
    std::mutex error_mutex;
    std::exception_ptr first_error;
    SpecUtilsAsync::ThreadPool pool;
    for( DistributedSrcCalcT<T> *calc : element_only )
    {
      pool.post( [calc,&error_mutex,&first_error](){
        try
        {
          self_shielding_integration_imp( *calc );
        }catch( std::exception & )
        {
          std::lock_guard<std::mutex> lock( error_mutex );
          if( !first_error )
            first_error = std::current_exception();
        }
      } );
    }
    pool.join();
    if( first_error )
      std::rethrow_exception( first_error );
  }else
  {
    for( DistributedSrcCalcT<T> *calc : element_only )
      self_shielding_integration_imp( *calc );
  }
}//integrate_element_calculators(...)


/** Integrates every calculator: line groups where the line path applies, the element path
 elsewhere.  A calculator whose source extent is below #sm_line_path_extent_ratio_floor is
 integrated through a floored copy (#floor_source_extent) and receives that copy's result.
 Replaces the per-calculator `self_shielding_integration_imp` loops in the fit and display paths. */
template<typename T>
void integrate_volumetric_calculators( const std::vector<std::unique_ptr<DistributedSrcCalcT<T>>> &calculators,
                                       const bool multithread )
{
  using namespace std;

  const VolumetricPartitionT<T> part = partition_volumetric_calculators( calculators );

  integrate_element_calculators( part.element_only, multithread );

  // Line groups (each parallel over line chunks internally).
  for( const vector<DistributedSrcCalcT<T>*> &g : part.line_groups )
    line_source_integration_imp( g, multithread );

  for( const pair<DistributedSrcCalcT<T>*,unique_ptr<DistributedSrcCalcT<T>>> &fl : part.floored )
  {
    fl.first->integral = fl.second->integral;
    fl.first->m_num_evals = fl.second->m_num_evals;
    fl.first->m_est_rel_error = fl.second->m_est_rel_error;
  }
}//integrate_volumetric_calculators(...)


/** The effective-shielding components (#EffShieldComponents) of every calculator, in the
 calculators' order - the post-fit diagnostic pass behind `computeEffectiveShielding`.  Line groups
 accumulate them in one line pass (#line_source_integration_imp with `eff_out`); everything the line
 path does not serve goes through the element path's #integrate_effective_shielding.  Same
 partition, floor and override rules as #integrate_volumetric_calculators.

 The two paths define the report differently, on purpose: the element path weights each element's
 CENTRE-RAY path by the element's whole contribution, the line path weights every line's OWN path by
 that line's contribution - the attenuation-weighted path of the photons actually detected.  They
 coincide far from the detector (measured 0.2%) and differ at contact for an optically thick source,
 where the line path's <AD> is LOWER (9-42% measured on steel at 1 cm; the short paths dominate
 exp(-tau)).  EffectiveShieldingLineVsElement in test_VolumetricLinePath.cpp records both. */
inline std::vector<EffShieldComponents> integrate_effective_shielding_all(
                          const std::vector<std::unique_ptr<DistributedSrcCalcT<double>>> &calculators,
                          const bool multithread )
{
  using namespace std;

  vector<EffShieldComponents> answer( calculators.size() );
  map<const DistributedSrcCalcT<double>*,size_t> index;
  for( size_t i = 0; i < calculators.size(); ++i )
    index[calculators[i].get()] = i;

  const VolumetricPartitionT<double> part = partition_volumetric_calculators( calculators );
  for( const pair<DistributedSrcCalcT<double>*,unique_ptr<DistributedSrcCalcT<double>>> &fl : part.floored )
    index[fl.second.get()] = index.at( fl.first );

  for( DistributedSrcCalcT<double> *calc : part.element_only )
    answer[index.at(calc)] = integrate_effective_shielding( *calc );

  for( const vector<DistributedSrcCalcT<double>*> &g : part.line_groups )
  {
    vector<EffShieldComponents> comps;
    line_source_integration_imp( g, multithread, &comps );
    for( size_t c = 0; c < g.size(); ++c )
      answer[index.at(g[c])] = comps[c];
  }

  return answer;
}//integrate_effective_shielding_all(...)

}//namespace GammaInteractionCalc

#endif //VolumetricLineIntegration_imp_hpp
