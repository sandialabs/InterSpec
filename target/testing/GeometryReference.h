#ifndef GeometryReference_h
#define GeometryReference_h
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

/** Independent reference implementations of the ray-tracing primitives in GammaInteractionCalc,
 plus formulation-independent invariant checks on their results.

 Unlike LegacyGeometryRef.h - which is a *frozen, known-wrong* copy of what shipped - this file is
 maintained and authoritative.  Each reference here deliberately uses a *different formulation* than
 the production code, so agreement between the two is evidence about the math rather than a
 restatement of it:

   - reference_cyl_exit      : slab-interval intersection (the t-interval inside the infinite
                               cylinder, intersected with the |z| <= half_length slab), vs
                               production's "solve the quadratic, then fall back to an end cap".
   - reference_rect_exit     : Kay/Kajiya slab traversal (t_enter = max lo_i, t_exit = min hi_i), vs
                               production's "pick the face on the ray's side, take the min t".
   - reference_sphere_exit   : the *stable* quadratic root, `g = (S-r)*(S+r)` and `t = g/(h+k)`, vs
                               production's radius-scaled discriminant (and the legacy expansion).

 A NOTE ON PRECISION: these are typed `long double` because that is free extra precision on x86.
 It buys nothing on Apple silicon - `long double` is binary64 there (`__LDBL_MANT_DIG__ == 53`), so
 the references are exactly `double` precision on the primary dev machine.  Their accuracy therefore
 has to come from *well-conditioned algebra*, not from the type; see the comments on
 reference_sphere_exit for the two cancellations that are avoided by hand.  Any check that needs
 genuinely-better-than-double reference values must be gated on `LDBL_MANT_DIG > 53`.

 The invariant checkers at the bottom (on_sphere / in_box / on_box / collinear / dist_matches) are
 the strongest oracles in this file, since they share no algebra at all with the code under test.
 */

#include <cmath>
#include <cfloat>
#include <limits>
#include <algorithm>

namespace GeomRef
{

/** An independent reference implementation of #cylinder_line_intersection, done in long double.

 The cylinder is centered at the origin and oriented along z.  Parameterizing the ray as
 P(t) = source + t*(detector - source), t increases monotonically toward the detector, so the
 t-interval inside the volume is the intersection of the infinite-cylinder interval with the
 |z| <= half_length slab interval; "toward detector" is that interval's high end, and "away from
 detector" its low end.

 Only valid for a source strictly inside the cylinder - which is the case the volumetric-source
 integrands exercise, and the one where the answer is unambiguous.

 @returns Whether the ray intersects the cylinder at all.

 (Moved here from test_GammaInteractionCalc.cpp so the parity suite can share it.)
 */
inline bool reference_cyl_exit( const double radius, const double half_length,
                                const double source[3], const double detector[3],
                                const bool toward_detector,
                                double exit_point[3], double &dist )
{
  typedef long double ld;

  const ld inf = std::numeric_limits<ld>::infinity();
  const ld rad = radius, half_z = half_length;
  const ld sx = source[0], sy = source[1], sz = source[2];
  const ld dx = ld(detector[0]) - sx, dy = ld(detector[1]) - sy, dz = ld(detector[2]) - sz;

  // The t-interval inside the infinite cylinder
  ld t_cyl_lo = -inf, t_cyl_hi = inf;
  const ld a = dx*dx + dy*dy;
  if( a > 0.0L )
  {
    const ld b = 2.0L*(sx*dx + sy*dy);
    const ld c = sx*sx + sy*sy - rad*rad;
    const ld disc = b*b - 4.0L*a*c;

    if( disc < 0.0L )
      return false;

    const ld sqrt_disc = sqrtl( disc );
    t_cyl_lo = (-b - sqrt_disc) / (2.0L*a);
    t_cyl_hi = (-b + sqrt_disc) / (2.0L*a);
  }else if( (sx*sx + sy*sy) > rad*rad )
  {
    return false;  //parallel to z, and outside the radius
  }

  // The t-interval inside the |z| <= half_length slab
  ld t_slab_lo = -inf, t_slab_hi = inf;
  if( dz != 0.0L )
  {
    const ld t_a = (-half_z - sz) / dz;
    const ld t_b = ( half_z - sz) / dz;
    t_slab_lo = (std::min)( t_a, t_b );
    t_slab_hi = (std::max)( t_a, t_b );
  }else if( fabsl(sz) > half_z )
  {
    return false;  //perpendicular to z, and past an end cap
  }

  const ld t_lo = (std::max)( t_cyl_lo, t_slab_lo );
  const ld t_hi = (std::min)( t_cyl_hi, t_slab_hi );

  if( t_lo > t_hi )
    return false;

  const ld t = toward_detector ? t_hi : t_lo;

  exit_point[0] = static_cast<double>( sx + t*dx );
  exit_point[1] = static_cast<double>( sy + t*dy );
  exit_point[2] = static_cast<double>( sz + t*dz );
  dist = static_cast<double>( fabsl(t) * sqrtl(dx*dx + dy*dy + dz*dz) );

  return true;
}//reference_cyl_exit(...)


/** An independent reference implementation of #exit_point_of_sphere_z / #sphere_exit_distance_scaled.

 Sphere of radius `sphere_rad` centered at the origin, detector at (0,0,observation_dist), source
 strictly inside (or on) the sphere.  Returns the (non-negative) distance from source to the chosen
 intersection, and writes that intersection to `exit_point`.

 With `u` the unit direction from the source toward the detector, `k = P.u` and `r = |P|`, the ray
 parameter satisfies `t^2 + 2*k*t - g = 0` with `g = S^2 - r^2`, giving `t = -k +/- sqrt(k^2 + g)`.
 Two hand-avoided cancellations make this a genuinely better-conditioned oracle than either
 production formulation:

   1. `g = (S - r)*(S + r)`, NOT `S*S - r*r`.  The difference-of-squares form loses roughly
      log2( S/(S-r) ) bits as the source approaches the surface - which is exactly the thin-shell
      regime the refactor is about.  The factored form's relative error is bounded by the inputs'.
   2. Each root is taken from whichever of the two algebraically-equal forms does not cancel, keyed
      off the sign of `k` (with `h = sqrt(k*k + g) >= |k|`, so `h+k` and `h-k` are each small exactly
      when `k` has the matching sign and `g` is small):

        t_plus  = -k + h  ==  g / (h + k)      -> use `h - k` when k <= 0, else `g/(h+k)`
        t_minus = -k - h  == -g / (h - k)      -> use `-(h + k)` when k >= 0, else `-g/(h-k)`

      Getting this branch wrong is not academic: for a near-surface source with a long transverse
      chord, `k` is a large negative number and `h + k` cancels to ~1e-6 of its operands, which
      costs ten significant digits.

 `positiveSolution` selects the forward root (toward the detector), else the backward root behind the
 source; the returned distance is non-negative either way, matching the production convention.
 */
inline double reference_sphere_exit( const double source[3], const double sphere_rad,
                                     const double observation_dist, const bool positiveSolution,
                                     double exit_point[3] )
{
  typedef long double ld;

  const ld S = sphere_rad;
  const ld px = source[0], py = source[1], pz = source[2];

  // Unit direction from the source toward the detector at (0,0,observation_dist)
  ld ux = -px, uy = -py, uz = ld(observation_dist) - pz;
  const ld un = sqrtl( ux*ux + uy*uy + uz*uz );
  ux /= un;  uy /= un;  uz /= un;

  const ld r = sqrtl( px*px + py*py + pz*pz );
  const ld k = px*ux + py*uy + pz*uz;
  const ld g = (S - r)*(S + r);            // == S^2 - r^2, without the cancellation
  const ld h = sqrtl( k*k + g );           // >= |k| for an interior source (g >= 0)

  // Pick the non-cancelling form of the wanted root (see the header comment above).
  ld t;
  if( positiveSolution )
    t = (k <= 0.0L) ? (h - k) : ( ((h + k) != 0.0L) ? (g / (h + k)) : (h - k) );
  else
    t = (k >= 0.0L) ? -(h + k) : ( ((h - k) != 0.0L) ? (-g / (h - k)) : -(h + k) );

  exit_point[0] = static_cast<double>( px + t*ux );
  exit_point[1] = static_cast<double>( py + t*uy );
  exit_point[2] = static_cast<double>( pz + t*uz );

  return static_cast<double>( fabsl(t) );
}//reference_sphere_exit(...)


/** An independent reference implementation of #rectangle_exit_location, via a Kay/Kajiya slab
 traversal rather than production's "pick the face on the ray's side, take the min t".

 Axis-aligned box of half-extents `half[3]` centered at the origin; `source` inside (or on) it.
 Returns the distance to where the ray toward `detector` leaves the box.

 Handles a zero half-extent exactly, and in the two physically distinct ways: an oblique ray through
 a zero-thickness slab has a zero chord (that slab's interval collapses to a point), while a ray
 *parallel* to the degenerate slab's planes is unconstrained by them (interval is all of R), so
 another axis sets the chord.  This is precisely the behavior #rectangle_exit_location_imp has and
 the frozen legacy copy does not.
 */
inline double reference_rect_exit( const double half[3], const double source[3],
                                   const double detector[3], double exit_point[3] )
{
  typedef long double ld;

  const ld inf = std::numeric_limits<ld>::infinity();

  ld d[3] = { ld(detector[0]) - source[0], ld(detector[1]) - source[1], ld(detector[2]) - source[2] };
  const ld dn = sqrtl( d[0]*d[0] + d[1]*d[1] + d[2]*d[2] );
  d[0] /= dn;  d[1] /= dn;  d[2] /= dn;

  ld t_exit = inf;
  for( int i = 0; i < 3; ++i )
  {
    if( d[i] == 0.0L )
      continue;   //ray parallel to this axis' planes: never leaves through them

    const ld t_a = (-ld(half[i]) - ld(source[i])) / d[i];
    const ld t_b = ( ld(half[i]) - ld(source[i])) / d[i];
    t_exit = (std::min)( t_exit, (std::max)( t_a, t_b ) );
  }//for( each axis )

  if( t_exit == inf )
    t_exit = 0.0L;   //degenerate: ray parallel to all three (impossible for a real direction)

  exit_point[0] = static_cast<double>( ld(source[0]) + t_exit*d[0] );
  exit_point[1] = static_cast<double>( ld(source[1]) + t_exit*d[1] );
  exit_point[2] = static_cast<double>( ld(source[2]) + t_exit*d[2] );

  return static_cast<double>( t_exit );
}//reference_rect_exit(...)


/** An independent reference implementation of #rectangle_intersections, same slab traversal as
 #reference_rect_exit but for a source and detector both *outside* the box.

 Fills the entry and exit points; returns whether the segment's supporting ray actually crosses the
 box in front of the source (`t_enter <= t_exit` and `t_exit >= 0`).
 */
inline bool reference_rect_intersections( const double half[3], const double source[3],
                                          const double detector[3],
                                          double enter_point[3], double exit_point[3] )
{
  typedef long double ld;

  const ld inf = std::numeric_limits<ld>::infinity();

  ld d[3] = { ld(detector[0]) - source[0], ld(detector[1]) - source[1], ld(detector[2]) - source[2] };
  const ld dn = sqrtl( d[0]*d[0] + d[1]*d[1] + d[2]*d[2] );
  d[0] /= dn;  d[1] /= dn;  d[2] /= dn;

  ld t_enter = -inf, t_exit = inf;
  for( int i = 0; i < 3; ++i )
  {
    if( d[i] == 0.0L )
    {
      if( fabsl( ld(source[i]) ) > ld(half[i]) )
        return false;   //parallel to this axis' planes, and outside the slab
      continue;
    }//if( parallel to this axis' planes )

    const ld t_a = (-ld(half[i]) - ld(source[i])) / d[i];
    const ld t_b = ( ld(half[i]) - ld(source[i])) / d[i];
    t_enter = (std::max)( t_enter, (std::min)( t_a, t_b ) );
    t_exit  = (std::min)( t_exit,  (std::max)( t_a, t_b ) );
  }//for( each axis )

  if( (t_enter > t_exit) || (t_exit < 0.0L) )
    return false;

  for( int i = 0; i < 3; ++i )
  {
    enter_point[i] = static_cast<double>( ld(source[i]) + t_enter*d[i] );
    exit_point[i]  = static_cast<double>( ld(source[i]) + t_exit *d[i] );
  }

  return true;
}//reference_rect_intersections(...)


// ---------------------------------------------------------------------------------------------
//  Formulation-independent invariants.  These assert *geometric properties* of a result rather
//  than re-deriving the result, so they share no algebra with the code under test.
// ---------------------------------------------------------------------------------------------

/** Is `pt` on the surface of an origin-centered sphere of radius `S`, to a relative tolerance? */
inline bool on_sphere( const double pt[3], const double S, const double rel_tol )
{
  const double r = std::sqrt( pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2] );
  return ( std::fabs(r - S) <= rel_tol*std::max( S, r ) );
}

/** Is `pt` inside (or on) an origin-centered box of half-extents `half`, to an absolute tolerance? */
inline bool in_box( const double pt[3], const double half[3], const double abs_tol )
{
  for( int i = 0; i < 3; ++i )
  {
    if( std::fabs(pt[i]) > (half[i] + abs_tol) )
      return false;
  }
  return true;
}

/** Is `pt` on the *surface* of the box - i.e. inside it, and touching at least one face? */
inline bool on_box( const double pt[3], const double half[3], const double abs_tol )
{
  if( !in_box( pt, half, abs_tol ) )
    return false;

  for( int i = 0; i < 3; ++i )
  {
    if( std::fabs( std::fabs(pt[i]) - half[i] ) <= abs_tol )
      return true;
  }
  return false;
}

/** Does `mid` lie on the ray from `from` toward `toward`?  Checked as a normalized cross-product
 magnitude, so it is scale-free and independent of how any of the three points was computed. */
inline bool collinear( const double from[3], const double mid[3], const double toward[3],
                       const double rel_tol )
{
  const double a[3] = { mid[0]-from[0], mid[1]-from[1], mid[2]-from[2] };
  const double b[3] = { toward[0]-from[0], toward[1]-from[1], toward[2]-from[2] };

  const double an = std::sqrt( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] );
  const double bn = std::sqrt( b[0]*b[0] + b[1]*b[1] + b[2]*b[2] );

  if( (an == 0.0) || (bn == 0.0) )
    return true;   //degenerate; nothing to say

  const double cx = a[1]*b[2] - a[2]*b[1];
  const double cy = a[2]*b[0] - a[0]*b[2];
  const double cz = a[0]*b[1] - a[1]*b[0];

  // Also require `mid` to be on the correct side (not behind `from`)
  const double dot = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  if( dot < -rel_tol*an*bn )
    return false;

  return ( std::sqrt( cx*cx + cy*cy + cz*cz ) <= rel_tol*an*bn );
}

/** Does the returned scalar `t` actually equal |exit - source|?  (The production functions all
 promise this, but compute the two through different expressions.) */
inline bool dist_matches( const double t, const double source[3], const double exit_pt[3],
                          const double rel_tol )
{
  const double dx = exit_pt[0]-source[0], dy = exit_pt[1]-source[1], dz = exit_pt[2]-source[2];
  const double d = std::sqrt( dx*dx + dy*dy + dz*dz );
  return ( std::fabs(t - d) <= rel_tol*std::max( 1.0, std::max(std::fabs(t), d) ) );
}

}//namespace GeomRef

#endif //GeometryReference_h
