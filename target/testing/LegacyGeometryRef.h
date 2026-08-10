#ifndef LegacyGeometryRef_h
#define LegacyGeometryRef_h
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

/** FROZEN pre-refactor copies of the double-valued ray-tracing primitives that used to live in
 src/GammaInteractionCalc.cpp, kept ONLY as a cross-check oracle for the migration of those
 functions onto the shared templated implementations in InterSpec/GammaInteractionCalc_imp.hpp.

     *** DO NOT "FIX" ANYTHING IN THIS FILE. ***

 Two of these bodies are known to be WRONG (both in rectangle_exit_location, each documented on the
 function below and each pinned by a test case in test_GeometryLegacyParity.cpp that asserts the bug
 is STILL HERE).  A third (exit_point_of_sphere_z) merely differs in formulation, at comparable
 accuracy.  Their entire value is that they are what shipped: the parity suite sweeps them against
 the new implementations to prove the migration changed nothing except those documented divergences.
 If a function here disagrees with the production one somewhere new, that is a finding about the
 production code - not a licence to edit this file.

 Frozen from commit 4c85390ae9ea4322dc52a7edd2d223ebc0044e14, src/GammaInteractionCalc.cpp lines
 574-644 (exit_point_of_sphere_z), 1173-1395 (cylinder_line_intersection), 1657-1768
 (rectangle_exit_location) and 1772-1945 (rectangle_intersections).

 Transcription rules applied (and the ONLY edits made to the frozen bodies):
   - `inline`, and moved into namespace LegacyGeom.
   - `noexcept` dropped; an oracle is never on a hot path.
   - Every `assert(...)` replaced in-place by a `// contract (legacy): ...` comment.  The parity
     driver deliberately feeds inputs outside the legacy contract (zero half-extents, sources
     exactly on a face, near-surface sphere points), which a Debug build would otherwise abort on
     before producing any comparison.
   - The `#if( DEBUG_RAYTRACE_CALCS )` tracing blocks and the stdout mutex dropped.
   - `CylExitDir` qualified as `GammaInteractionCalc::CylExitDir`.
   - Two pieces of dead code in `rectangle_exit_location` that existed only to serve the removed
     asserts: the `#ifndef NDEBUG const double src_copy[3] = {...}; #endif` snapshot (otherwise
     unused once the asserts are comments), and the trailing `assert(0); return 0.0;` after the
     if/else-if/else (unreachable - every branch returns).  Neither carries arithmetic.
 Every other line, including the TODOs, is verbatim - they are part of the frozen record.

 `distance` is not copied (it is a trivially identical three-subtraction Euclidean norm), and
 `point_to_line_dist` is not copied (it had zero callers and was deleted, not migrated).

 Retire this file, GeometryReference.h, and test_GeometryLegacyParity.cpp once the migration has
 shipped a release cycle.
 */

#include <cmath>
#include <cfloat>
#include <stdexcept>

#include "InterSpec/GammaInteractionCalc.h"

namespace LegacyGeom
{

/** Frozen copy of the pre-refactor double `GammaInteractionCalc::exit_point_of_sphere_z`.

 DIVERGENCE (#3 of 3), and NOT a bug: the replacement (#exit_point_of_sphere_z_imp ->
 #sphere_exit_distance_scaled) factors the radius out of the discriminant, which is what keeps the
 ceres::Jet derivative finite as the radius -> 0, but it is not more accurate than this hand-expanded
 closed form.  Measured in `SphereExitNearSurfaceAccuracy`: near the shell surface the two trade
 wins, and both land within ~1e-13 of the sphere radius.  That case also explains why (both
 formulations cancel when |P| -> S, just in different places) and what would actually fix it.

 Note the `cerr` diagnostics that preceded each throw are dropped here (the production replacement
 folds those values into the exception message instead).
 */
inline double exit_point_of_sphere_z( const double source_point[3],
                                      double exit_point[3],
                                      double sphere_rad,
                                      double observation_dist,
                                      bool postiveSolution = true )
{
  /*
   *Makes x0,y0,z0 be the point of intersection of a sphere or radius
   *'sphere_rad' centered at the origin, with the line pointing from <x0,y0,z0>
   *towards <0,0,observation_dist>
   *Note: CAN be used with a 2-dimensional integral which holds phi constant
   */

  using namespace std;

  const double a = source_point[0];
  const double b = source_point[1];
  const double c = source_point[2];
  const double &S = sphere_rad;
  const double &R = observation_dist;

  const double r = sqrt( a*a + b*b + c*c );
  if( observation_dist < sphere_rad )
  {
    throw runtime_error( "exit_point_of_sphere_z(...): obs_dist < sphere_rad" );
  }//if( observation_dist < sphere_rad )

  if( r > sphere_rad )
  {
    if( ((r-sphere_rad)/sphere_rad) < 0.0001 )
    {
      exit_point[0] = source_point[0];
      exit_point[1] = source_point[1];
      exit_point[2] = source_point[2];
      return 0.0;
    }//if( this is just a rounding error )

    throw runtime_error( "exit_point_of_sphere_z(...): r > sphere_rad" );
  }//if( r > sphere_rad )

  //TODO: factor the math below to save CPU time
  if( postiveSolution )
  {
    exit_point[0] =-(a*sqrt((R*R-2*c*R+c*c+b*b+a*a)*S*S+(-b*b-a*a)*R*R)-a*R*R+a*c*R)/(R*R-2*c*R+c*c+b*b+a*a);
    exit_point[1] =(-b*sqrt(R*R*S*S-2*c*R*S*S+c*c*S*S+b*b*S*S+a*a*S*S-b*b*R*R-a*a*R*R)+b*R*R-b*c*R)/(R*R-2*c*R+c*c+b*b+a*a);
    exit_point[2] =(R*(sqrt(R*R*S*S-2*c*R*S*S+c*c*S*S+b*b*S*S+a*a*S*S-b*b*R*R-a*a*R*R)+b*b+a*a)-c*sqrt(R*R*S*S-2*c*R*S*S+c*c*S*S+b*b*S*S+a*a*S*S-b*b*R*R-a*a*R*R))/(R*R-2*c*R+c*c+b*b+a*a);

    //  const double factor_1 = R*R-2*c*R+c*c+b*b+a*a;
    //  const double factor_2 = R*R*S*S-2*c*R*S*S+c*c*S*S+b*b*S*S+a*a*S*S-b*b*R*R-a*a*R*R;
    //  const double x_pos_new =-(a*sqrt((factor_1)*S*S+(-b*b-a*a)*R*R)-a*R*R+a*c*R)/factor_1;
    //  const double y_pos_new =(-b*sqrt(factor_2)+b*R*R-b*c*R)/factor_1;
    //  const double z_pos_new =(R*(sqrt(factor_2)+b*b+a*a)-c*sqrt(factor_2))/(R*R-2*c*R+c*c+b*b+a*a);
    //  assert( x_pos_new == exit_point[0] );
    //  assert( y_pos_new == exit_point[1] );
    //  assert( z_pos_new == exit_point[2] );
  }else
  {
    exit_point[0] = (a*sqrt((R*R-2*c*R+c*c+b*b+a*a)*S*S+(-b*b-a*a)*R*R)+a*R*R-a*c*R)/(R*R-2*c*R+c*c+b*b+a*a);
    exit_point[1] = (b*sqrt(R*R*S*S-2*c*R*S*S+c*c*S*S+b*b*S*S+a*a*S*S-b*b*R*R-a*a*R*R)+b*R*R-b*c*R)/(R*R-2*c*R+c*c+b*b+a*a);
    exit_point[2] = (c*sqrt(R*R*S*S-2*c*R*S*S+c*c*S*S+b*b*S*S+a*a*S*S-b*b*R*R-a*a*R*R)+R*(-sqrt(R*R*S*S-2*c*R*S*S+c*c*S*S+b*b*S*S+a*a*S*S-b*b*R*R-a*a*R*R)+b*b+a*a))/(R*R-2*c*R+c*c+b*b+a*a);
  }//if( postiveSolution ) / else

  const double dx = a - exit_point[0];
  const double dy = b - exit_point[1];
  const double dz = c - exit_point[2];

  return sqrt( dx*dx + dy*dy + dz*dz );
}//exit_point_of_sphere_z


/** Frozen copy of the pre-refactor double `GammaInteractionCalc::cylinder_line_intersection`.

 NO known divergence: after commits c1f9983e and 4c85390a this was line-for-line identical to
 #cylinder_line_intersection_imp, so the parity suite checks it bit-exactly.
 */
inline double cylinder_line_intersection( const double radius, const double half_length,
                              const double source[3],
                              const double detector[3],
                              const GammaInteractionCalc::CylExitDir direction,
                              double exit_point[3] )
{
  using namespace std;
  using GammaInteractionCalc::CylExitDir;

  // TODO: this function should be broken into two separate functions.  One to handle finding the exit
  //  point when you know the source is inside the volume.  And one to find both intersection points
  //  (if any) of external points.  This would both increase the efficiency of the function, and also
  //  make the use cleaner/easier.

  // TODO: need to clearly define what happens for points exactly on boundary, or just over or
  //       whatever (I think things are set up so equality is counted as being inside volume,
  //       but there may be edge cases (pun realized) that dont obey this.

  // contract (legacy): radius >= 0.0
  // contract (legacy): half_length >= 0.0
  // contract (legacy): direction is TowardDetector or AwayFromDetector

  // A convenience function for handling case where the line never enters our volume.
  auto handle_line_outside_volume = [&exit_point,&source]() -> double {
    exit_point[0] = source[0];
    exit_point[1] = source[1];
    exit_point[2] = source[2];

    return 0.0;
  };//handle_line_outside_volume lamda


  if( (radius <= 0.0) || (half_length <= 0.0) )
    return handle_line_outside_volume();


  // Get unit direction vector from source to final position
  double unit[3] = { detector[0] - source[0], detector[1] - source[1], detector[2] - source[2] };

  {// begin scope on norm
    const double norm = sqrt( unit[0]*unit[0] + unit[1]*unit[1] + unit[2]*unit[2] );
    unit[0] /= norm;
    unit[1] /= norm;
    unit[2] /= norm;
  }// end scope on norm

  // Check if parallel to z-axis
  //  We should probably compare to DBL_MIN, but realistically anything less than DBL_EPSILON
  //  is close enough to zero for our purposes (the DRF will fail far before this assumption fails)
  if( (fabs(unit[0]) < DBL_EPSILON) && (fabs(unit[1]) < DBL_EPSILON) )
  {
    // contract (legacy): fabs(unit[2]) > DBL_EPSILON

    const double r = sqrt(source[0]*source[0] + source[1]*source[1]);

    if( r > radius )
      return handle_line_outside_volume();

    double exit_z = (unit[2] > 0.0) ? half_length : -half_length;
    switch( direction )
    {
      case CylExitDir::TowardDetector:
        break;

      case CylExitDir::AwayFromDetector:
        exit_z *= -1.0;
        break;
    }//switch( direction )


    const double distance = fabs( exit_z - source[2] );

    exit_point[0] = source[0];
    exit_point[1] = source[1];
    exit_point[2] = exit_z;

    return distance;
  }//if( (fabs(unit[0]) < DBL_EPSILON) && (fabs(unit[1]) < DBL_EPSILON) )


  // Make sure both source and detector z-coordinates are not on the same end of the cylinder,
  //  and larger than the half length.
  if( (signbit(source[2]) == signbit(detector[2]))
     && (fabs(source[2]) > half_length)
     && (fabs(detector[2]) > half_length) )
  {
    return handle_line_outside_volume();
  }


  // Check that the source point isnt on the same side of the circle as the detector, but outside
  //  of the circles radius.  This will prevent the case where the infinite line would intersect
  //  the cylinder, but not between source and detector (which we should return zero for)
  // TODO: get rid of these sqrts, and also probably other ones throughout this function
  const double src_rad = sqrt( source[0]*source[0] + source[1]*source[1] );
  const double det_rad = sqrt( detector[0]*detector[0] + detector[1]*detector[1] );

  // A detector on the cylinder axis is never "on the same side" as the source, so this check does
  //  not apply to it (and the direction it is in is undefined).
  if( (src_rad >= radius) && (det_rad > 0.0) )
  {
    //There is probably a better way to check if source and detector are within 90 degrees of each
    //  other
    const double src_unit_x = source[0] / src_rad;
    const double src_unit_y = source[1] / src_rad;

    const double det_unit_x = detector[0] / det_rad;
    const double det_unit_y = detector[1] / det_rad;

    const double unit_dx = src_unit_x - det_unit_x;
    const double unit_dy = src_unit_y - det_unit_y;

    const double unit_dist_2 = unit_dx*unit_dx + unit_dy*unit_dy;

    // TODO: for radius=1, half_length=1, source={1,-1,0}, detector{1,1,0}, because of numerical
    //       rounding (unit_dist_2 = (2 - 4E-16) - i.e. 2*DBL_EPSILON) we will return in this next
    //       statement, which we probably shouldnt since it is an exact intersection.  So we should
    //       probably implement a check that accounts for these rounding errors
    if( unit_dist_2 <= 2.0 )
      return handle_line_outside_volume();
  }//if( src_rad > radius )


  // Find where the ray crosses the infinite cylinder (x^2 + y^2 = radius^2), parameterizing the ray
  //  as P(t) = source + t*(detector - source).  So t==0 is the source, t==1 is the detector, and t
  //  increases monotonically toward the detector - which makes the "toward detector" crossing
  //  unambiguously the larger t, and "away from detector" the smaller one.
  //
  // Note: do *not* order the two crossings by their distance to the detector.  For end-on geometry
  //  the detector lies on the cylinder axis, so both crossings are exactly equidistant from it, and
  //  the pick becomes a rounding-noise coin flip that sends the ray out the wrong end cap.
  const double dx = detector[0] - source[0];
  const double dy = detector[1] - source[1];
  const double dz = detector[2] - source[2];

  // Substituting P(t) into x^2 + y^2 = radius^2 gives the quadratic a*t^2 + b*t + c = 0.
  //  a > 0 here, since the parallel-to-z case was already handled above.
  const double a = dx*dx + dy*dy;
  const double b = 2.0*(source[0]*dx + source[1]*dy);
  const double c = source[0]*source[0] + source[1]*source[1] - radius*radius;

  const double discriminant = b*b - 4.0*a*c;

  if( discriminant < 0.0 )
    return handle_line_outside_volume();

  // Use the numerically stable form of the quadratic formula; the textbook form loses precision in
  //  one of the two roots when |b| ~ sqrt(discriminant), i.e. for rays that just glance the cylinder.
  const double sqrt_disc = sqrt( discriminant );
  const double q = -0.5*(b + ((b < 0.0) ? -sqrt_disc : sqrt_disc));
  const double t_root_a = q / a;
  const double t_root_b = (q != 0.0) ? (c / q) : t_root_a; //q is only zero when b and discriminant are
  const double t_near = (t_root_a < t_root_b) ? t_root_a : t_root_b;
  const double t_far  = (t_root_a < t_root_b) ? t_root_b : t_root_a;

  const bool toward_det = (direction == CylExitDir::TowardDetector);
  const double t_exit  = toward_det ? t_far  : t_near;
  const double t_other = toward_det ? t_near : t_far;

  double x_exit = source[0] + t_exit*dx;
  double y_exit = source[1] + t_exit*dy;
  double z_exit = source[2] + t_exit*dz;

  // If z_exit is past the half-length, then for the ray to intersect our *finite* cylinder at all,
  //  the other crossing of the infinite cylinder must be either inside the volume, or on its far
  //  side; if it is past the same end, the ray misses us entirely.
  if( fabs(z_exit) > half_length )
  {
    const double other_z_exit = source[2] + t_other*dz;

    if( (fabs(other_z_exit) > half_length) && (signbit(z_exit) == signbit(other_z_exit)) )
      return handle_line_outside_volume();

    // We leave through one of the end caps - solve for where on that disk.
    // contract (legacy): dz != 0.0

    const double z_cap = ((z_exit < 0.0) ? -half_length : half_length);
    const double t_cap = (z_cap - source[2]) / dz;

    x_exit = source[0] + t_cap*dx;
    y_exit = source[1] + t_cap*dy;
    z_exit = z_cap;
  }//if( fabs(z_exit) > half_length )

  // If we are here, we are guaranteed the line does go through our volume

  const double exit_dx = (source[0] - x_exit);
  const double exit_dy = (source[1] - y_exit);
  const double exit_dz = (source[2] - z_exit);
  const double dist_scaler = sqrt( exit_dx*exit_dx + exit_dy*exit_dy + exit_dz*exit_dz );

  exit_point[0] = x_exit;
  exit_point[1] = y_exit;
  exit_point[2] = z_exit;

  return dist_scaler;
}//double cylinder_line_intersection(...)


/** Frozen copy of the pre-refactor double `GammaInteractionCalc::rectangle_exit_location`.

 KNOWN DIVERGENCE (#1 of 3) - the "zero-transverse-dimension bug": the `t = (intersect - source)*inv_slope`
 form below evaluates `(0 - 0)*DBL_MAX == 0` when a half-extent is zero, so the degenerate plane wins the
 nearest-exit comparison and a zero chord is returned.  Pinned by `RectExitZeroExtentNewIsCorrect`.

 KNOWN DIVERGENCE (#2 of 3) - same root cause, but reachable with strictly positive dimensions: when the
 ray is parallel to an axis' planes (`norm[i] == 0`) AND the source lies exactly on that face
 (`source[i] == +/-half_i`), the product is again `0*DBL_MAX == 0` and a zero chord is returned.  This is
 reachable from `eval_rect`'s outer-shell loop, where the source is by construction the previous shell's
 exit point - i.e. exactly on a face.  Pinned by `RectExitSourceOnParallelFace`.

 #rectangle_exit_location_imp fixes both by short-circuiting the whole `t` to DBL_MAX when `norm[i]` is
 zero, and never dividing by a half-extent.
 */
inline double rectangle_exit_location( const double half_width, const double half_height,
                               const double half_depth,
                               const double source[3],
                               const double detector[3],
                               double exit_point[3] )
{
  using namespace std;

  // contract (legacy): half_width > 0.0
  // contract (legacy): half_height > 0.0
  // contract (legacy): half_depth > 0.0

  // contract (legacy): fabs(source[0]) <= (half_width + 1.0E-12)
  // contract (legacy): fabs(source[1]) <= (half_height + 1.0E-12)
  // contract (legacy): fabs(source[2]) <= (half_depth + 1.0E-12)

  // The detector is normally outside the box, but the point-source attenuation path can
  //  pass a detector that sits inside a (degenerate) over-thick shield; in that case the
  //  plane-intersection logic below still returns the far-surface exit (t > source-detector
  //  distance), which the caller caps at the true detector distance.  So we do not assert
  //  the detector is outside here.

  // contract (legacy): source and detector are not the same point

  double norm[3] = { detector[0] - source[0], detector[1] - source[1], detector[2] - source[2] };

  // Currently (20211126) the detector will be [0.0,0.0,m_observationDist], so we know which face
  //  of the detector the ray will exit, so we could first check for this case with this commented
  //  out code, and then return an answer efficiently; need to benchmark before bothering to
  //  uncomment this special case.
  //if( (fabs(detector[0]) <= half_width) && (fabs(detector[1]) <= half_height) )
  //{
  //  const double z_frac_inside = fabs( (half_depth - source[2]) / (detector[2] - source[2]) );
  //  exit_point[0] = source[0] + z_frac_inside * norm[0];
  //  exit_point[1] = source[1] + z_frac_inside * norm[1];
  //  exit_point[2] = ((detector[2] > 0.0) ? half_depth : -half_depth); //std::copysign(half_depth,detector[2]);
  //
  //  return distance( source, exit_point );
  //}//if( we know ray is exiting the face on positive/negative z )


  const double total_dist = sqrt( norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2] );
  norm[0] /= total_dist;
  norm[1] /= total_dist;
  norm[2] /= total_dist;

  // We'll find the intersection of all the possible planes, and then choose the one that
  //  is the closest to source
  // Recall equation of a line in three-space is:
  //   x = x_0 + t*a --> source[0] + (t * norm[0])
  //   y = z_0 + t*b --> source[1] + (t * norm[1])
  //   z = z_0 + t*c --> source[2] + (t * norm[2])

  const double inv_slope_x = (norm[0] == 0.0) ? DBL_MAX : (1.0 / norm[0]);
  const double x_intersect = (inv_slope_x >= 0.0) ? half_width : -half_width;
  const double t_intersect_x = (x_intersect - source[0])*inv_slope_x;

  const double inv_slope_y = (norm[1] == 0.0) ? DBL_MAX : (1.0 / norm[1]);
  const double y_intersect = (inv_slope_y >= 0.0) ? half_height : -half_height;
  const double t_intersect_y = (y_intersect - source[1])*inv_slope_y;

  const double inv_slope_z = (norm[2] == 0.0) ? DBL_MAX : (1.0 / norm[2]);
  const double z_intersect = (inv_slope_z >= 0.0) ? half_depth : -half_depth;
  const double t_intersect_z = (z_intersect - source[2])*inv_slope_z;

  if( (t_intersect_x <= t_intersect_y) && (t_intersect_x <= t_intersect_z) )
  {
    // We are exiting through the plane perpendicular to x-axis
    exit_point[0] = ((norm[0] >= 0.0) ? half_width : -half_width);
    exit_point[1] = source[1] + (t_intersect_x * norm[1]);
    exit_point[2] = source[2] + (t_intersect_x * norm[2]);

    // contract (legacy): the x-plane crossing lands on the +/-half_width face
    // contract (legacy): t_intersect_x equals the source-to-exit distance

    return t_intersect_x;
  }else if( t_intersect_y <= t_intersect_z )
  {
    // We are exiting through the plane perpendicular to y-axis
    exit_point[0] = source[0] + (t_intersect_y * norm[0]);
    exit_point[1] = ((norm[1] >= 0.0) ? half_height : -half_height);
    exit_point[2] = source[2] + (t_intersect_y * norm[2]);

    // contract (legacy): the y-plane crossing lands on the +/-half_height face
    // contract (legacy): t_intersect_y equals the source-to-exit distance

    return t_intersect_y;
  }else
  {
    // We are exiting through the plane perpendicular to z-axis

    exit_point[0] = source[0] + (t_intersect_z * norm[0]);
    exit_point[1] = source[1] + (t_intersect_z * norm[1]);
    exit_point[2] = ((norm[2] >= 0.0) ? half_depth : -half_depth);

    // contract (legacy): the z-plane crossing lands on the +/-half_depth face
    // contract (legacy): t_intersect_z equals the source-to-exit distance

    return t_intersect_z;
  }// if / else figure out where we are exiting.
}//rectangle_exit_location(...)


/** Frozen copy of the pre-refactor double `GammaInteractionCalc::rectangle_intersections`.

 NO known divergence: #rectangle_intersections_imp is an algorithmically identical transcription,
 including the same `inv_slope * DBL_MAX` formulation - it did NOT receive the zero-extent fix that
 #rectangle_exit_location_imp did, so both share that latent bug and the parity suite checks them
 bit-exactly over the (strictly positive) dimensions the one production caller uses.
 */
inline bool rectangle_intersections( const double half_width, const double half_height,
                             const double half_depth,
                             const double source[3],
                             const double detector[3],
                             double enter_point[3],
                             double exit_point[3] )
{
  using namespace std;

  // Only checking inputs sanity on debug builds since this is a hot-path function, and is only
  //  called from one spot, so it should be good to only check inputs on development builds.

  // Make sure we arent passing garbage dimensions in ever
  // contract (legacy): (half_width > 0.0) && (half_height > 0.0) && (half_depth > 0.0)

  // Dev test that both the source and detector points are outside of the box; only checking on
  //  debug builds because this should really be the case
  // contract (legacy): source is outside the box
  // contract (legacy): detector is outside the box

  // Make sure detector and source arent in same position.
  // contract (legacy): source and detector are not the same point

  // See notes in #rectangle_exit_location about

  double norm[3] = { detector[0] - source[0], detector[1] - source[1], detector[2] - source[2] };
  const double total_dist = sqrt( norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2] );
  norm[0] /= total_dist;
  norm[1] /= total_dist;
  norm[2] /= total_dist;

  const double inv_slope_x = (norm[0] == 0.0) ? DBL_MAX : (1.0 / norm[0]);
  const double x_intersect = (inv_slope_x >= 0.0) ? half_width : -half_width;
  const double t_intersect_x_det = (x_intersect - source[0])*inv_slope_x;
  const double t_intersect_x_src = (-x_intersect - source[0])*inv_slope_x;


  const double inv_slope_y = (norm[1] == 0.0) ? DBL_MAX : (1.0 / norm[1]);
  const double y_intersect = (inv_slope_y >= 0.0) ? half_height : -half_height;
  const double t_intersect_y_det = (y_intersect - source[1])*inv_slope_y;
  const double t_intersect_y_src = (-y_intersect - source[1])*inv_slope_y;

  const double inv_slope_z = (norm[2] == 0.0) ? DBL_MAX : (1.0 / norm[2]);
  const double z_intersect = (inv_slope_z >= 0.0) ? half_depth : -half_depth;
  const double t_intersect_z_det = (z_intersect - source[2])*inv_slope_z;
  const double t_intersect_z_src = (-z_intersect - source[2])*inv_slope_z;

  const bool intersects_x_src = (t_intersect_x_src >= 0.0);
  const bool intersects_y_src = (t_intersect_y_src >= 0.0);
  const bool intersects_z_src = (t_intersect_z_src >= 0.0);

  const bool x_before_y_src = (!intersects_y_src || (t_intersect_x_src <= t_intersect_y_src));
  const bool x_before_z_src = (!intersects_z_src || (t_intersect_x_src <= t_intersect_z_src));
  const bool y_before_z_src = (!intersects_z_src || (t_intersect_y_src <= t_intersect_z_src));

  if( intersects_x_src && x_before_y_src && x_before_z_src )
  {
    // We are exiting through the plane perpendicular to x-axis
    const double src_intersect_x = ((norm[0] >= 0.0) ? -half_width : half_width);
    const double src_intersect_y = source[1] + (t_intersect_x_src * norm[1]);
    const double src_intersect_z = source[2] + (t_intersect_x_src * norm[2]);

    if( (fabs(src_intersect_y) > half_height)
       || (fabs(src_intersect_z) > half_depth) )
    {
      // contract (legacy): min(t_src) <= min(t_det) + 1.0E-9
      return false;
    }

    enter_point[0] = src_intersect_x;
    enter_point[1] = src_intersect_y;
    enter_point[2] = src_intersect_z;

    // contract (legacy): the x-plane entry lands on the +/-half_width face
  }else if( intersects_y_src && y_before_z_src )
  {
    // We are exiting through the plane perpendicular to y-axis
    const double src_intersect_x = source[0] + (t_intersect_y_src * norm[0]);
    const double src_intersect_y = ((norm[1] >= 0.0) ? -half_height : half_height);
    const double src_intersect_z = source[2] + (t_intersect_y_src * norm[2]);

    if( (fabs(src_intersect_x) > half_width)
       || (fabs(src_intersect_z) > half_depth) )
    {
      // contract (legacy): min(t_src) <= min(t_det) + 1.0E-9
      return false;
    }

    enter_point[0] = src_intersect_x;
    enter_point[1] = src_intersect_y;
    enter_point[2] = src_intersect_z;

    // contract (legacy): the y-plane entry lands on the +/-half_height face
  }else if( intersects_z_src )
  {
    // We are exiting through the plane perpendicular to z-axis
    const double src_intersect_x = source[0] + (t_intersect_z_src * norm[0]);
    const double src_intersect_y = source[1] + (t_intersect_z_src * norm[1]);
    const double src_intersect_z = ((norm[2] >= 0.0) ? -half_depth : half_depth);

    if( (fabs(src_intersect_x) > half_width)
       || (fabs(src_intersect_y) > half_height) )
    {
      // contract (legacy): min(t_src) <= min(t_det) + 1.0E-9
      return false;
    }

    enter_point[0] = src_intersect_x;
    enter_point[1] = src_intersect_y;
    enter_point[2] = src_intersect_z;

    // contract (legacy): the z-plane entry lands on the +/-half_depth face
  }else
  {
    return false;
  }// if / else figure out where we are exiting.

  const bool intersects_x_det = (t_intersect_x_det >= 0.0);
  const bool intersects_y_det = (t_intersect_y_det >= 0.0);
  const bool intersects_z_det = (t_intersect_z_det >= 0.0);

  const bool x_before_y_det = (!intersects_y_det || (t_intersect_x_det <= t_intersect_y_det));
  const bool x_before_z_det = (!intersects_z_det || (t_intersect_x_det <= t_intersect_z_det));
  const bool y_before_z_det = (!intersects_z_det || (t_intersect_y_det <= t_intersect_z_det));



  if( intersects_x_det && x_before_y_det && x_before_z_det )
  {
    // We are exiting through the plane perpendicular to x-axis
    exit_point[0] = ((norm[0] >= 0.0) ? half_width : -half_width);
    exit_point[1] = source[1] + (t_intersect_x_det * norm[1]);
    exit_point[2] = source[2] + (t_intersect_x_det * norm[2]);

    // contract (legacy): the x-plane exit lands on the +/-half_width face
  }else if( intersects_y_det && y_before_z_det )
  {
    // We are exiting through the plane perpendicular to y-axis
    exit_point[0] = source[0] + (t_intersect_y_det * norm[0]);
    exit_point[1] = ((norm[1] >= 0.0) ? half_height : -half_height);
    exit_point[2] = source[2] + (t_intersect_y_det * norm[2]);

    // contract (legacy): the y-plane exit lands on the +/-half_height face
  }else if( intersects_z_det )
  {
    // We are exiting through the plane perpendicular to z-axis
    exit_point[0] = source[0] + (t_intersect_z_det * norm[0]);
    exit_point[1] = source[1] + (t_intersect_z_det * norm[1]);
    exit_point[2] = ((norm[2] >= 0.0) ? half_depth : -half_depth);

    // contract (legacy): the z-plane exit lands on the +/-half_depth face
  }else
  {
    // contract (legacy): unreachable
  }// if / else figure out where we are exiting.


  // contract (legacy): |enter_point| <= half + 1.0E-9 on each axis
  // contract (legacy): |exit_point|  <= half + 1.0E-9 on each axis

  return true;
}//rectangle_intersections( ... )

}//namespace LegacyGeom

#endif //LegacyGeometryRef_h
