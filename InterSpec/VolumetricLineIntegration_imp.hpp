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
 2-4 point Gauss-Legendre average in y = exp(-mu (s1 - s)).

 What this buys over the element path (measured there: 5.8e-4 s per element at 512 rays, boxes
 needing 11k-67k elements PER ENERGY, nothing shared across energies or evaluations): the cost is
 (lines) x (energies) exponentials with the detector work hoisted out of the fit entirely, and
 with the lines fixed the chords are exact in T, so d(integral)/d(dims) no longer carries the
 frozen-aperture staircase error the element path documents at eval_cylinder.

 DIRECTION PROPOSAL.  Directions are aimed at points drawn UNIFORMLY IN THE (padded) SOURCE SOLID
 rather than into a cone around its bounding sphere: for a point q uniform in a volume V_p the
 induced direction density is p(w) = (s1^3 - s0^3)/(3 V_p) over the chord [s0,s1] of the line
 through V_p, and the line's weight carries 1/p.  Every line then crosses the padded solid by
 construction, whatever its aspect ratio - where a cone around the bounding SPHERE would waste most
 of its lines on empty space for a needle or an edge-on sheet - and as a fitted dimension shrinks
 the lines follow it.  The proposal is scalar (frozen for autodiff); that is exact, because the true
 integral does not depend on it (importance sampling with a fixed proposal is unbiased for every
 parameter value inside its support, and the padding keeps the support).  The set is rebuilt only
 when the scalar source dimensions leave the window #VolumetricLineCache::matches accepts, never
 per Jet pass - see that function for why holding one proposal across a neighbourhood matters more
 than aiming it perfectly.

 LIMIT.  A source extent that is a vanishing fraction of the standoff makes chord/volume a 0/0 here.
 The intervals are computed from each line's closest approach to the assembly origin
 (#line_shell_intervals_imp), which keeps the chords accurate down to an extent/distance ratio of
 ~1e-10; below that the dispatcher (#integrate_volumetric_calculators) integrates a copy of the
 source whose vanishing extent is floored at that ratio (#sm_line_path_extent_ratio_floor), so the
 value and derivative stay finite and continuous all the way to exactly zero and the line path owns
 the whole domain (test_ShieldingDimLimit pins this).

 UNITS.  CeeLo works in cm; everything here is in PhysicalUnits.  The etendue weight of a line is
 `omega_w` (cm^2, already divided by 4 pi) times cm^2 in PhysicalUnits; multiplied by an emission
 density (1/volume) and a chord (length) it is dimensionless.
 */

#include <map>
#include <array>
#include <cmath>
#include <mutex>
#include <tuple>
#include <memory>
#include <vector>
#include <cassert>
#include <algorithm>
#include <type_traits>

#include "io/DetectorEtendue.h"
#include "io/LowDiscrepancy.h"
#include "io/DetectorResponse.h"

namespace GammaInteractionCalc
{

/** Which quadrature a volumetric calculator is integrated with. */
enum class VolumetricIntegrator : int
{
  /** Production choice: the line path wherever it applies (a response is attached, no cascade
   correction, no effective-AN/AD accumulation, no collimator), else the element path. */
  Auto,
  /** Force the per-element aperture path (the reference implementation). */
  Element,
  /** Force the line path (throws for calculators it cannot serve). */
  Line
};//enum class VolumetricIntegrator

/** TEST HOOK - overrides the integrator choice for every calculator in the process.  Leave at
 Auto in production; the tests use it to A/B the two quadratures on identical calculators. */
inline VolumetricIntegrator sm_volumetric_integrator_override = VolumetricIntegrator::Auto;

/** Gauss-Legendre points on a source chord, in the y = exp(-mu_eff (s1 - s)) substitution
 (2 = exact for a linear remainder; radial in-situ profiles use 4 - see #line_source_integration_imp). */
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
  ceelo::ResponseFlag worst_flag = ceelo::ResponseFlag::Ok;
  double max_frac_sigma = 0.0;

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


/** Builds #PrefactorGrid for `energy_keV` over distances [d_lo_cm, d_hi_cm] (log-spaced) and
 cos_theta in [0,1]; phi over [0,90] deg for Quadrant responses. */
inline std::shared_ptr<const PrefactorGrid> build_prefactor_grid( const ceelo::DetectorResponse &resp,
                                                                  const double energy_keV,
                                                                  const double d_lo_cm,
                                                                  const double d_hi_cm,
                                                                  const size_t num_d = 48,
                                                                  const size_t num_cos = 33,
                                                                  const size_t num_phi_quadrant = 9 )
{
  assert( (d_lo_cm > 0.0) && (d_hi_cm > d_lo_cm) );
  auto grid = std::make_shared<PrefactorGrid>();
  grid->energy_keV = energy_keV;

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

  // The quadrature argument of fep_prefactor only feeds the collimator shadow gate; the
  //  dispatcher keeps collimated responses on the element path, so an empty one is correct here.
  const ceelo::ApertureQuadrature no_quadrature;

  grid->ln_p.resize( grid->ln_d.size() * grid->cos_t.size() * grid->phi_deg.size() );
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
        const ceelo::EffResult pre = resp.fep_prefactor( energy_keV, pos, no_quadrature );
        if( !(pre.value > 0.0) )
          throw std::runtime_error( "build_prefactor_grid: non-positive prefactor" );
        grid->ln_p[grid->index(id,ic,ip)] = std::log( pre.value );
        if( static_cast<int>(pre.flag) > static_cast<int>(grid->worst_flag) )
          grid->worst_flag = pre.flag;
        grid->max_frac_sigma = std::max( grid->max_frac_sigma, pre.sigma / pre.value );
      }
    }
  }

  return grid;
}//build_prefactor_grid(...)


/** The detector-side line set for one source shell of one fit, plus the per-energy memos.
 Built by #build_volumetric_line_cache; owned (shared) by #ShieldingSourceChi2Fcn and referenced by
 every #DistributedSrcCalcT of that source.  Immutable after construction except the memos, which
 are mutex-guarded. */
struct VolumetricLineCache
{
  /** Key: what the set was built for. */
  std::shared_ptr<const ceelo::DetectorResponse> response;
  GeometryType geometry = GeometryType::NumGeometryType;
  size_t material_index = 0;
  std::array<double,3> source_outer_dims = { 0.0, 0.0, 0.0 };   //scalar, PhysicalUnits
  std::array<double,3> det_position = { 0.0, 0.0, 0.0 };
  std::array<double,3> det_axis = { 0.0, 0.0, -1.0 };
  double det_azimuth = 0.0;
  int num_lines = 0;

  /** Crystal (CeeLo) -> assembly rotation, and the reference point in the crystal frame (cm). */
  double M[3][3] = { {1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,1.0} };
  std::array<double,3> ref_c = { 0.0, 0.0, 0.0 };

  /** The padded proposal solid the directions were aimed at (assembly frame, PhysicalUnits). */
  std::array<double,3> proposal_dims = { 0.0, 0.0, 0.0 };

  /** The lines, in the crystal frame (cm) - `lines.q` carries the CeeLo kernel segments. */
  ceelo::EtendueLineSet lines;

  /** Per kept line, assembly frame, PhysicalUnits: origin relative to the detector position,
   unit photon direction (toward the detector), etendue weight (area), and the distance back
   along the line from the hull point to the endcap-front plane (the air path ends there). */
  std::vector<std::array<double,3>> origin_rel;
  std::vector<std::array<double,3>> dir;
  std::vector<double> weight;
  std::vector<double> s_endcap;

  /** Distances (cm, from the crystal-face origin) the prefactor grids cover. */
  double prefactor_d_lo_cm = 0.0, prefactor_d_hi_cm = 0.0;

  mutable std::mutex memo_mutex;
  mutable std::map<double,std::shared_ptr<const std::vector<double>>> kernel_by_energy;
  mutable std::map<double,std::shared_ptr<const PrefactorGrid>> prefactor_by_energy;

  /** Per-line FEP interaction probability at `energy_keV` (memoized). */
  std::shared_ptr<const std::vector<double>> kernel( const double energy_keV ) const
  {
    {
      std::lock_guard<std::mutex> lock( memo_mutex );
      const auto pos = kernel_by_energy.find( energy_keV );
      if( pos != end(kernel_by_energy) )
        return pos->second;
    }
    auto k = std::make_shared<std::vector<double>>();
    response->fep_line_probabilities( energy_keV, lines.q, *k );
    std::lock_guard<std::mutex> lock( memo_mutex );
    return kernel_by_energy.emplace( energy_keV, k ).first->second;
  }//kernel(...)

  /** The prefactor grid at `energy_keV` (memoized). */
  std::shared_ptr<const PrefactorGrid> prefactor( const double energy_keV ) const
  {
    {
      std::lock_guard<std::mutex> lock( memo_mutex );
      const auto pos = prefactor_by_energy.find( energy_keV );
      if( pos != end(prefactor_by_energy) )
        return pos->second;
    }
    std::shared_ptr<const PrefactorGrid> g = build_prefactor_grid( *response, energy_keV,
                                                                   prefactor_d_lo_cm, prefactor_d_hi_cm );
    std::lock_guard<std::mutex> lock( memo_mutex );
    return prefactor_by_energy.emplace( energy_keV, g ).first->second;
  }//prefactor(...)

  /** Whether this cache can still serve the given scalar configuration.

   The detector placement, response and line count must match exactly, but the SOURCE dimensions
   only have to be close: the line set is an importance-sampling proposal, and a proposal aimed at
   a slightly different solid is still unbiased as long as it covers the real one (the padding, see
   #build_volumetric_line_cache).  Reusing it matters for the fit, not for the cost: a rebuilt set
   is a different quadrature, so the estimate moves by its own discretisation error (~0.2% at 65536
   lines) whenever the dimensions move at all - and an objective that jumps by 0.2% between
   optimizer steps is one Levenberg-Marquardt cannot take clean steps on.  Holding one proposal
   across a whole neighbourhood makes the objective smooth there, which is what the optimizer
   actually needs; the tolerance is well inside the 1.5x padding, so coverage is never at risk.

   MEASURED, and this is the known weakness of the scheme: the window makes rebuilds RARE, it does
   not make them free.  `LineCacheRebuildContinuity` sweeps a source radius across a boundary and
   finds a step of ~1.6e-3 in the integral there - the proposal is re-aimed by the whole width of
   the window at once, so the sets either side are effectively independent quadratures and the
   estimate moves by their discretisation.  A fit that walks a source dimension across a boundary
   therefore meets a ~0.16% discontinuity, which is exactly what this window was meant to avoid and
   only defers.  The fix is NOT a wider window - a fit crosses it eventually - but to build the
   proposal ONCE per fit, padded to the dimension parameters' upper BOUNDS instead of their current
   value, so it spans the whole search domain and never needs rebuilding.  Tracked in TODO.md. */
  bool matches( const ceelo::DetectorResponse *resp, const GeometryType geom, const size_t mat_index,
                const std::array<double,3> &outer_dims, const std::array<double,3> &det_pos,
                const std::array<double,3> &axis, const double azimuth, const int n ) const
  {
    if( (response.get() != resp) || (geometry != geom) || (material_index != mat_index)
        || (det_position != det_pos) || (det_axis != axis) || (det_azimuth != azimuth)
        || (num_lines != n) )
      return false;

    // Dimensions at or below the extent floor all describe the same proposal (the floor is what the
    //  sampler aimed at), so a shell collapsing to zero does not rebuild on every evaluation.
    const double ext_floor = sm_line_path_extent_ratio_floor
                             * std::sqrt( det_pos[0]*det_pos[0] + det_pos[1]*det_pos[1] + det_pos[2]*det_pos[2] );
    for( size_t i = 0; i < 3; ++i )
    {
      const double have = std::max( source_outer_dims[i], ext_floor );
      const double want = std::max( outer_dims[i], ext_floor );
      if( have == want )
        continue;
      const double ratio = want/have;
      if( (ratio < 0.8) || (ratio > 1.2) )
        return false;
    }
    return true;
  }
};//struct VolumetricLineCache


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
inline double solid_volume( const GeometryType geometry, const std::array<double,3> &dims )
{
  const double pi = PhysicalUnits::pi;
  switch( geometry )
  {
    case GeometryType::Spherical:      return (4.0/3.0)*pi*dims[0]*dims[0]*dims[0];
    case GeometryType::CylinderEndOn:
    case GeometryType::CylinderSideOn: return 2.0*pi*dims[0]*dims[0]*dims[1];
    case GeometryType::Rectangular:    return 8.0*dims[0]*dims[1]*dims[2];
    case GeometryType::NumGeometryType: break;
  }
  assert( 0 );
  return 0.0;
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


/** Builds the line set for one source shell.  `source_outer_dims` are the SCALAR cumulative outer
 dims of the source shell (the solid the directions are aimed at, padded by `pad`), `det_*` the
 scalar detector geometry (position = the response's reference point in the assembly frame). */
inline std::shared_ptr<const VolumetricLineCache> build_volumetric_line_cache(
                                          std::shared_ptr<const ceelo::DetectorResponse> response,
                                          const GeometryType geometry,
                                          const size_t material_index,
                                          const std::array<double,3> &source_outer_dims,
                                          const std::array<double,3> &det_position,
                                          const std::array<double,3> &det_axis,
                                          const double det_azimuth,
                                          const int num_lines,
                                          const double pad = 1.5 )
{
  using namespace std;
  const double cm = PhysicalUnits::cm;

  if( !response || (num_lines <= 0) )
    throw runtime_error( "build_volumetric_line_cache: invalid inputs" );

  auto cache = make_shared<VolumetricLineCache>();
  cache->response = response;
  cache->geometry = geometry;
  cache->material_index = material_index;
  cache->source_outer_dims = source_outer_dims;
  cache->det_position = det_position;
  cache->det_axis = det_axis;
  cache->det_azimuth = det_azimuth;
  cache->num_lines = num_lines;

  detector_frame_rotation( det_axis.data(), det_azimuth, cache->M );
  const Eigen::Vector3d r0 = response->reference_point_position();
  cache->ref_c = { r0.x(), r0.y(), r0.z() };

  // Proposal solid: the source padded in every dimension, floored at the same extent ratio the
  //  dispatcher floors the source itself at, so the proposal always covers the (floored) source.
  const double det_dist = std::sqrt( det_position[0]*det_position[0] + det_position[1]*det_position[1]
                                     + det_position[2]*det_position[2] );
  const double ext_floor = sm_line_path_extent_ratio_floor * det_dist;
  std::array<double,3> pdims = source_outer_dims;
  const int ndims = (geometry == GeometryType::Spherical) ? 1 : ((geometry == GeometryType::Rectangular) ? 3 : 2);
  for( int i = 0; i < ndims; ++i )
    pdims[i] = pad * std::max( std::fabs(pdims[i]), ext_floor );
  cache->proposal_dims = pdims;
  const double vol_p = solid_volume( geometry, pdims );

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

  std::vector<ceelo::HullPoint> hull;
  ceelo::sample_hull_points( geom, num_lines, toward, have_dir, 0, hull );

  cache->lines.q.n_rays_total = num_lines;
  cache->lines.target_centre = centre_c;
  cache->lines.target_radius = solid_bounding_radius( geometry, pdims )/cm;
  cache->lines.q.rays.reserve( num_lines );
  cache->lines.origin.reserve( num_lines );
  cache->lines.dir.reserve( num_lines );
  cache->origin_rel.reserve( num_lines );
  cache->dir.reserve( num_lines );
  cache->weight.reserve( num_lines );
  cache->s_endcap.reserve( num_lines );

  double d_min_cm = 1.0e300, d_max_cm = 0.0;   //source-point distance range from the crystal origin
  const double inv_n = 1.0/static_cast<double>(num_lines);

  for( int i = 0; i < num_lines; ++i )
  {
    const ceelo::HullPoint &hp = hull[i];
    const uint64_t idx = static_cast<uint64_t>(i);

    // Hull point and its normal in the assembly frame (PhysicalUnits).
    const Eigen::Vector3d xc_rel = (hp.point - r0)*cm;
    double x_rel[3], n_a[3];
    to_assembly_dir( xc_rel, x_rel );
    to_assembly_dir( hp.normal, n_a );
    const double x_a[3] = { det_position[0] + x_rel[0], det_position[1] + x_rel[1], det_position[2] + x_rel[2] };

    // Target point uniform in the padded source solid; outward direction toward it.
    const std::array<double,3> q = uniform_point_in_solid( geometry, pdims,
                                     ceelo::halton(idx,7), ceelo::halton(idx,11), ceelo::halton(idx,13) );
    double w[3] = { q[0] - x_a[0], q[1] - x_a[1], q[2] - x_a[2] };
    const double wn = std::sqrt( w[0]*w[0] + w[1]*w[1] + w[2]*w[2] );
    if( !(wn > 0.0) )
      continue;
    for( int k = 0; k < 3; ++k )
      w[k] /= wn;

    Eigen::Vector3d w_c;
    to_crystal_dir( w, w_c );
    const double cos_n = w[0]*n_a[0] + w[1]*n_a[1] + w[2]*n_a[2];
    if( !(cos_n > 0.0) || !(w_c.z() < 0.0) )
      continue;   //leaves the hull elsewhere / behind the face plane: weight zero, not renormalised

    // Chord of this line through the PADDED solid: the direction density p(w) = (s1^3 - s0^3)/(3 V_p).
    const double o_d[3] = { x_a[0], x_a[1], x_a[2] };
    const double d_a[3] = { -w[0], -w[1], -w[2] };   //photon direction
    double s0p, s1p;
    if( !line_shell_interval_imp<double>( geometry, pdims, o_d, d_a, s0p, s1p ) )
      continue;   //grazes the proposal boundary at round-off level
    s0p = std::max( s0p, 0.0 );
    if( !(s1p > s0p) )
      continue;
    const double s3 = s1p*s1p*s1p - s0p*s0p*s0p;
    if( !(s3 > 0.0) )
      continue;
    const double inv_p = 3.0*vol_p/s3;   //sr

    const double etendue_cm2_sr = hp.area_weight * cos_n * inv_p * inv_n;
    cache->lines.total_etendue += etendue_cm2_sr;

    if( !ceelo::append_etendue_line( cache->lines, geom, hp.point, w_c, etendue_cm2_sr ) )
      continue;

    cache->origin_rel.push_back( { x_rel[0], x_rel[1], x_rel[2] } );
    cache->dir.push_back( { d_a[0], d_a[1], d_a[2] } );
    cache->weight.push_back( cache->lines.q.rays.back().omega_w * cm * cm );
    // Distance back along the line from the hull point to the endcap-front plane z = r0.z.
    cache->s_endcap.push_back( std::max( 0.0, (hp.point.z() - r0.z())/(-w_c.z()) ) * cm );

    // Distance range of the padded chord from the crystal origin (for the prefactor grids).
    for( const double s : { s0p, s1p } )
    {
      const Eigen::Vector3d p_c = hp.point + (s/cm)*w_c;
      const double d = p_c.norm();
      d_min_cm = std::min( d_min_cm, d );
      d_max_cm = std::max( d_max_cm, d );
    }
  }//for( loop over lines )

  if( cache->lines.q.rays.empty() )
    throw runtime_error( "build_volumetric_line_cache: no line reaches the active crystal" );

  // Prefactor grids span the padded chords with a further margin, and never reach 0 distance.
  cache->prefactor_d_lo_cm = std::max( 0.5*d_min_cm, 1.0e-3 );
  cache->prefactor_d_hi_cm = std::max( 2.0*d_max_cm, cache->prefactor_d_lo_cm*2.0 );

  return cache;
}//build_volumetric_line_cache(...)


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
    default: x = x4; w = w4; return;
  }
}//unit_gauss_legendre(...)


/** Integrates a GROUP of calculators that share geometry (same source shell, dims, detector, line
 cache, normalization and in-situ settings) and differ only in energy-dependent coefficients:
 the per-line chord bookkeeping is done once, the energies innermost.  Fills `integral`,
 `m_num_evals` and `m_est_rel_error` on each.  Chunks of lines run on `pool` when given, and are
 reduced in a fixed order so the result does not depend on the thread count.

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
    assert( !calc->m_cascade );
  }

  const size_t num_calc = group.size();
  const size_t num_lines = cache->dir.size();
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
      const T L( relax );
      switch( geometry )
      {
        case GeometryType::Spherical:
          rho_scale = T(1.0) / (T(4.0)*pi*Do[0]*Do[0]*Do[0]*sphere_exp_norm_factor( Do[0]/L ));
          break;
        case GeometryType::CylinderEndOn:
          rho_scale = T(1.0) / (pi*Do[0]*Do[0]*T(2.0)*Do[1]*one_minus_exp_neg_over_x( T(2.0)*Do[1]/L ));
          break;
        case GeometryType::CylinderSideOn:
          rho_scale = T(1.0) / (T(4.0)*pi*Do[1]*Do[0]*Do[0]*cyl_side_exp_norm_factor( Do[0]/L ));
          break;
        case GeometryType::Rectangular:
          rho_scale = T(1.0) / (T(8.0)*Do[0]*Do[1]*Do[2]*one_minus_exp_neg_over_x( T(2.0)*Do[2]/L ));
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

  // Per-calculator energy-dependent inputs.
  std::vector<std::shared_ptr<const std::vector<double>>> kernels( num_calc );
  std::vector<std::shared_ptr<const PrefactorGrid>> grids( num_calc );
  for( size_t c = 0; c < num_calc; ++c )
  {
    kernels[c] = cache->kernel( group[c]->m_energy );
    grids[c] = cache->prefactor( group[c]->m_energy );
  }

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

  const size_t num_chunks = (multithread && (num_lines > 4096))
                              ? static_cast<size_t>( std::max( 1, SpecUtilsAsync::num_logical_cpu_cores() ) )
                              : size_t(1);
  // Partial sums: [chunk][calc], plus even/odd line halves for the error estimate.
  std::vector<std::vector<T>> partial( num_chunks, std::vector<T>( num_calc, T(0.0) ) );
  std::vector<std::vector<double>> partial_even( num_chunks, std::vector<double>( num_calc, 0.0 ) );
  std::vector<std::vector<double>> partial_odd( num_chunks, std::vector<double>( num_calc, 0.0 ) );
  std::vector<std::vector<EffShieldComponents>> partial_eff( accumulate_eff ? num_chunks : 0,
                                                             std::vector<EffShieldComponents>( num_calc ) );

  const auto do_chunk = [&]( const size_t chunk )
  {
    const size_t lo = (num_lines*chunk)/num_chunks;
    const size_t hi = (num_lines*(chunk + 1))/num_chunks;
    std::vector<T> &acc = partial[chunk];
    std::vector<double> &acc_even = partial_even[chunk];
    std::vector<double> &acc_odd = partial_odd[chunk];
    std::vector<EffShieldComponents> * const acc_eff = accumulate_eff ? &partial_eff[chunk] : nullptr;

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
      const std::array<double,3> &orel = cache->origin_rel[j];
      const std::array<double,3> &d = cache->dir[j];
      const T o[3] = { lead.m_detector.position[0] + orel[0],
                       lead.m_detector.position[1] + orel[1],
                       lead.m_detector.position[2] + orel[2] };

      line_shell_intervals_imp( geometry, shells, o, d.data(), a, b, crossed );
      if( !crossed[m] )
        continue;

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
      const T air_len = select_max( a[num_shells-1] - T(cache->s_endcap[j]), T(0.0) );

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

      const double w_line = cache->weight[j];
      const bool even = ((j & 1) == 0);

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
        const double k = (*kernels[c])[j];
        if( !(k > 0.0) )
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
              rate0 = T( d[2]/relax );
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
            const T dist = sqrt( pc[0]*pc[0] + pc[1]*pc[1] + pc[2]*pc[2] );
            const T cos_t = -pc[2]/dist;
            T phi( 0.0 );
            if( grids[c]->phi_deg.size() > 1 )
            {
              // Quadrant symmetry: fold into [0,90] degrees.
              const T ax = (scalar_of(pc[0]) < 0.0) ? -pc[0] : pc[0];
              const T ay = (scalar_of(pc[1]) < 0.0) ? -pc[1] : pc[1];
              phi = atan2( ay, ax ) * T(180.0/PhysicalUnits::pi);
            }
            T val = grids[c]->eval( log( dist ), cos_t, phi );
            if( radial_profile )
              val *= exp( -(depth_at( p ) - depth0)/T(relax) );
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

        const T contrib = T(w_line*k) * line_sum;
        acc[c] += contrib;
        if( even )
          acc_even[c] += scalar_of( contrib );
        else
          acc_odd[c] += scalar_of( contrib );
      }//for( calculators )
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

  for( size_t c = 0; c < num_calc; ++c )
  {
    T total( 0.0 );
    double even = 0.0, odd = 0.0;
    for( size_t chunk = 0; chunk < num_chunks; ++chunk )
    {
      total += partial[chunk][c];
      even += partial_even[chunk][c];
      odd += partial_odd[chunk][c];
    }
    group[c]->integral = total;
    group[c]->m_num_evals = num_lines;
    const double tot = std::fabs( scalar_of(total) );
    group[c]->m_est_rel_error = (tot > 0.0) ? (0.5*std::fabs(even - odd)/tot) : 0.0;
  }

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
  if( !calc.m_effResponse || !calc.m_lineCache || calc.m_cascade )
    return false;
  if( calc.m_effResponse->descriptor.collimator )
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
