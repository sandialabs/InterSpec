#ifndef GammaInteractionCalc_imp_hpp
#define GammaInteractionCalc_imp_hpp
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

/** Templated (double or ceres::Jet<>) implementations of the volumetric-source
 ray-tracing and integration used by the activity/shielding fit.

 These are transcriptions of the double-valued functions in
 GammaInteractionCalc.cpp (exit_point_of_sphere_z, cylinder_line_intersection,
 rectangle_exit_location/intersections, and DistributedSrcCalc's eval_*), made
 generic on the scalar type so that, with T = ceres::Jet<>, the derivatives of
 the volumetric integrals with respect to fit parameters (shell dimensions,
 generic-shield AN/AD, source activity) come out of the integration
 automatically.

 The Cuba/Cuhre integrator only works on doubles, so integration here is a
 self-contained h-adaptive tensor-product Gauss-Legendre scheme over the unit
 cube (#adaptive_unit_cube_integrate), with subdivision driven by the scalar
 part of T.

 Conventions match the double path: the source/shielding assembly is centered
 at the origin; shells are stored as cumulative outer dimensions; generic
 (AN/AD) shells have zero physical extent and their `trans_len_coef` is the
 total (dimensionless) attenuation; integration coordinates `xx[]` are on
 [0,1]^ndim with the Jacobian included in the integrand value.

 To use with ceres::Jet, include ceres.h (or ceres/jet.h) before this header.
 */

#include <map>
#include <array>
#include <cmath>
#include <mutex>
#include <tuple>
#include <cfloat>
#include <limits>
#include <memory>
#include <vector>
#include <cassert>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <functional>
#include <type_traits>

#include "SpecUtils/SpecUtilsAsync.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/MassAttenuationTool_imp.hpp"

namespace ceres
{
  /* dummy namespace for when using this file for only doubles, and ceres.h hasnt been included */
}

namespace GammaInteractionCalc
{

/** The scalar part of a double or ceres::Jet<> value. */
template<typename T>
inline double scalar_of( const T &val )
{
  if constexpr ( std::is_same_v<T,double> )
    return val;
  else
    return val.a;
}


template<typename T>
inline T distance_imp( const T a[3], const T b[3] )
{
  using namespace std;
  using namespace ceres;

  const T dx = a[0] - b[0];
  const T dy = a[1] - b[1];
  const T dz = a[2] - b[2];

  return sqrt( dx*dx + dy*dy + dz*dz );
}//distance_imp(...)


/** Templated equivalent of DetectorPeakResponse::fractionalSolidAngle(diam,dist). */
template<typename T>
inline T fractional_solid_angle_imp( const double detector_diameter, const T &dist )
{
  using namespace std;
  using namespace ceres;

  const double r = 0.5 * detector_diameter;
  return 0.5*(1.0 - (dist/sqrt(dist*dist + r*r)));
}//fractional_solid_angle_imp(...)


/** Whether the detector response depends only on the source-to-detector
 distance (as currently), or also on the angle of incidence.

 Symmetry shortcuts that rotate coordinates - notably treating an off-axis
 sphere as an on-axis sphere at the true line-of-sight distance - are only
 valid while this returns true; when angular detector response is added, this
 must return false (or become a DetectorPeakResponse property), and those
 shortcuts will fall back to general 3D integration.
 */
constexpr bool detector_response_is_isotropic()
{
  return true;
}


/** The detector, as seen from the assembly-centered coordinate system the
 ray-tracing works in.

 `position` is templated on T so that a future option where the user-entered
 distance is measured to the assembly *face* (making the detector position
 depend on fit dimensions) works without re-templating.
 */
template<typename T>
struct DetectorGeomT
{
  /** Center of the detector face, in assembly-centered coordinates. */
  T position[3];

  /** Unit vector of the detector face normal, pointing from the detector
   toward the assembly.  Not used by the current (isotropic) response - stored
   so a future angular detector response has its data available.
   */
  double axis[3];

  /** Radius of detector face. */
  double radius;

  /** Distance from detector face to actual detection element. */
  double setback;
};//struct DetectorGeomT


/** The single place that maps the user inputs (distance, and eventually
 off-axis offsets and a distance-to-face option) to a detector location in
 assembly-centered coordinates.

 The detector is along +z, except for CylinderSideOn where it is along +x
 (matching DistributedSrcCalc::eval_cylinder).
 */
template<typename T>
DetectorGeomT<T> detector_geom_from_config( const GeometryType geometry,
                                            const T &distance,
                                            const double detector_radius,
                                            const double detector_setback )
{
  DetectorGeomT<T> det;
  det.radius = detector_radius;
  det.setback = detector_setback;

  const bool is_side_on = (geometry == GeometryType::CylinderSideOn);

  det.position[0] = is_side_on ? distance : T(0.0);
  det.position[1] = T(0.0);
  det.position[2] = is_side_on ? T(0.0) : distance;

  det.axis[0] = is_side_on ? -1.0 : 0.0;
  det.axis[1] = 0.0;
  det.axis[2] = is_side_on ? 0.0 : -1.0;

  return det;
}//detector_geom_from_config(...)


/** The per-ray detector response: the probability a gamma emitted at
 `eval_point` toward the detector is incident on the detector face.

 This is the single change-point for adding an angular dependence to the
 detector efficiency: the angle of incidence is the angle between the
 (eval_point - det.position) ray and det.axis.  Until then, only the
 line-of-sight distance matters (see #detector_response_is_isotropic).
 */
template<typename T>
T detector_response_factor( const DetectorGeomT<T> &det, const T eval_point[3] )
{
  const T dist = distance_imp( eval_point, det.position );
  return fractional_solid_angle_imp( 2.0*det.radius, dist + det.setback );
}//detector_response_factor(...)


/** Templated transcription of #exit_point_of_sphere_z.

 Makes `exit_point` the intersection of the sphere of radius `sphere_rad`
 (centered at origin) with the line from `source_point` toward
 {0,0,observation_dist}; returns the distance from source to exit point.
 Safe to pass the same array for source and exit point.
 */
template<typename T>
T exit_point_of_sphere_z_imp( const T source_point[3],
                              T exit_point[3],
                              const T &sphere_rad,
                              const T &observation_dist,
                              const bool postiveSolution = true )
{
  using namespace std;
  using namespace ceres;

  const T a = source_point[0];
  const T b = source_point[1];
  const T c = source_point[2];
  const T &S = sphere_rad;
  const T &R = observation_dist;

  const T r = sqrt( a*a + b*b + c*c );
  if( observation_dist < sphere_rad )
    throw std::runtime_error( "exit_point_of_sphere_z_imp(...): obs_dist < sphere_rad" );

  if( r > sphere_rad )
  {
    if( ((r-sphere_rad)/sphere_rad) < 0.0001 )
    {
      exit_point[0] = source_point[0];
      exit_point[1] = source_point[1];
      exit_point[2] = source_point[2];
      return T(0.0);
    }//if( this is just a rounding error )

    throw std::runtime_error( "exit_point_of_sphere_z_imp(...): r > sphere_rad" );
  }//if( r > sphere_rad )

  const T denom = R*R - 2.0*c*R + c*c + b*b + a*a;
  const T sqrt_arg = R*R*S*S - 2.0*c*R*S*S + c*c*S*S + b*b*S*S + a*a*S*S - b*b*R*R - a*a*R*R;
  const T sqrt_term = sqrt( sqrt_arg );

  if( postiveSolution )
  {
    exit_point[0] = -(a*sqrt_term - a*R*R + a*c*R) / denom;
    exit_point[1] = (-b*sqrt_term + b*R*R - b*c*R) / denom;
    exit_point[2] = (R*(sqrt_term + b*b + a*a) - c*sqrt_term) / denom;
  }else
  {
    exit_point[0] = (a*sqrt_term + a*R*R - a*c*R) / denom;
    exit_point[1] = (b*sqrt_term + b*R*R - b*c*R) / denom;
    exit_point[2] = (c*sqrt_term + R*(-sqrt_term + b*b + a*a)) / denom;
  }//if( postiveSolution ) / else

  const T dx = a - exit_point[0];
  const T dy = b - exit_point[1];
  const T dz = c - exit_point[2];

  return sqrt( dx*dx + dy*dy + dz*dz );
}//exit_point_of_sphere_z_imp(...)


/** Templated transcription of #cylinder_line_intersection.

 Cylinder is along the z-axis, centered at origin.  Returns the distance from
 `source` to the exit point in the requested direction, or 0.0 if the line
 from `source` to `detector` does not intersect the cylinder volume.
 */
template<typename T>
T cylinder_line_intersection_imp( const T &radius, const T &half_length,
                                  const T source[3],
                                  const T detector[3],
                                  const CylExitDir direction,
                                  T exit_point[3] )
{
  using namespace std;
  using namespace ceres;

  assert( scalar_of(radius) >= 0.0 );
  assert( scalar_of(half_length) >= 0.0 );

  // A convenience function for handling case where the line never enters our volume.
  auto handle_line_outside_volume = [&exit_point,&source]() -> T {
    exit_point[0] = source[0];
    exit_point[1] = source[1];
    exit_point[2] = source[2];
    return T(0.0);
  };//handle_line_outside_volume lamda

  if( (radius <= 0.0) || (half_length <= 0.0) )
    return handle_line_outside_volume();

  // Get unit direction vector from source to final position
  T unit[3] = { detector[0] - source[0], detector[1] - source[1], detector[2] - source[2] };

  {// begin scope on norm
    const T norm = sqrt( unit[0]*unit[0] + unit[1]*unit[1] + unit[2]*unit[2] );
    unit[0] /= norm;
    unit[1] /= norm;
    unit[2] /= norm;
  }// end scope on norm

  // Check if parallel to z-axis
  if( (abs(unit[0]) < DBL_EPSILON) && (abs(unit[1]) < DBL_EPSILON) )
  {
    const T r = sqrt(source[0]*source[0] + source[1]*source[1]);

    if( r > radius )
      return handle_line_outside_volume();

    T exit_z = (unit[2] > 0.0) ? half_length : -half_length;
    switch( direction )
    {
      case CylExitDir::TowardDetector:
        break;

      case CylExitDir::AwayFromDetector:
        exit_z *= -1.0;
        break;
    }//switch( direction )

    const T dist = abs( exit_z - source[2] );

    exit_point[0] = source[0];
    exit_point[1] = source[1];
    exit_point[2] = exit_z;

    return dist;
  }//if( (abs(unit[0]) < DBL_EPSILON) && (abs(unit[1]) < DBL_EPSILON) )

  // Make sure both source and detector z-coordinates are not on the same end of the cylinder,
  //  and larger than the half length.
  if( (std::signbit(scalar_of(source[2])) == std::signbit(scalar_of(detector[2])))
     && (abs(source[2]) > half_length)
     && (abs(detector[2]) > half_length) )
  {
    return handle_line_outside_volume();
  }

  // Check that the source point isnt on the same side of the circle as the detector, but outside
  //  of the circles radius.  This will prevent the case where the infinite line would intersect
  //  the cylinder, but not between source and detector (which we should return zero for)
  const T src_rad = sqrt( source[0]*source[0] + source[1]*source[1] );
  if( src_rad >= radius )
  {
    const T src_unit_x = source[0] / src_rad;
    const T src_unit_y = source[1] / src_rad;

    const T det_rad = sqrt( detector[0]*detector[0] + detector[1]*detector[1] );
    const T det_unit_x = detector[0] / det_rad;
    const T det_unit_y = detector[1] / det_rad;

    const T unit_dx = src_unit_x - det_unit_x;
    const T unit_dy = src_unit_y - det_unit_y;

    const T unit_dist_2 = unit_dx*unit_dx + unit_dy*unit_dy;

    if( unit_dist_2 <= 2.0 )
      return handle_line_outside_volume();
  }//if( src_rad > radius )

  T x_exit, y_exit, z_exit;  //intersection in the direction of detector

  if( detector[0] != source[0] )
  {
    T other_x_exit;

    {// begin scope to solve for x,y of intersection
      // Line parametrically: P(t) = source + t*(detector-source); substitute into
      //  the circle equation x^2 + y^2 = r^2 and solve the resulting quadratic.
      const T dx = detector[0] - source[0];
      const T dy = detector[1] - source[1];

      const T a = dx*dx + dy*dy;
      const T b = 2.0*(source[0]*dx + source[1]*dy);
      const T c = source[0]*source[0] + source[1]*source[1] - radius*radius;

      const T discriminant = b*b - 4.0*a*c;

      if( discriminant < 0 )
        return handle_line_outside_volume();

      const T sqrt_discriminant = sqrt(discriminant);
      const T t1 = (-b - sqrt_discriminant) / (2.0*a);
      const T t2 = (-b + sqrt_discriminant) / (2.0*a);

      const T point1[2] = {source[0] + t1*dx, source[1] + t1*dy};
      const T point2[2] = {source[0] + t2*dx, source[1] + t2*dy};

      // Pick the solution towards/away-from the detector.
      //  The line parameter t increases from source (t=0) toward the detector (t=1), so the
      //  intersection with the larger t (always t2, since a > 0) is the one toward the detector.
      //  (Comparing xy-distances to the detector is an exact floating-point tie when the
      //  detector is on the cylinder axis, so must not be used here.)
      switch( direction )
      {
        case CylExitDir::TowardDetector:
          x_exit       = point2[0];
          y_exit       = point2[1];
          other_x_exit = point1[0];
          break;

        case CylExitDir::AwayFromDetector:
          x_exit       = point1[0];
          y_exit       = point1[1];
          other_x_exit = point2[0];
          break;
      }//switch( direction )
    }// end scope to solve for x,y of intersection

    // Now find the z corresponding to {x_exit, y_exit}.
    const T m_zx = unit[2] / unit[0];
    const T c_zx = source[2] - m_zx*source[0];
    z_exit = m_zx*x_exit + c_zx;

    // If the z_exit is past the half_length, lets check that the other intersection to the
    //  infinite cylinder happens either in the volume, or on the other side of our volume
    if( abs(z_exit) > half_length )
    {
      const T other_z_exit = m_zx*other_x_exit + c_zx;
      if( (abs(other_z_exit) > half_length)
          && (std::signbit(scalar_of(z_exit)) == std::signbit(scalar_of(other_z_exit))) )
        return handle_line_outside_volume();
    }//if( abs(z_exit) > half_length )
  }else
  {
    // The x coordinate of the source and detector location are the same.
    if( detector[0] > radius )
      return handle_line_outside_volume();

    x_exit = detector[0];

    y_exit = sqrt(radius*radius - x_exit*x_exit);

    switch( direction )
    {
      case CylExitDir::TowardDetector:   break;
      case CylExitDir::AwayFromDetector: y_exit *= -1.0; break;
    }//switch( direction )

    if( unit[1] < 0.0 )
      y_exit *= -1.0;

    // If detector and source x and y are the same, the line is parallel to the z-axis,
    //  which was already dealt with above.
    const T m_zy = unit[2] / unit[1];
    const T c_zy = source[2] - m_zy*source[1];
    z_exit = m_zy*y_exit + c_zy;

    if( abs(z_exit) > half_length )
    {
      const T other_y_exit = -y_exit;
      const T other_z_exit = m_zy*other_y_exit + c_zy;

      if( (abs(other_z_exit) > half_length)
          && (std::signbit(scalar_of(z_exit)) == std::signbit(scalar_of(other_z_exit))) )
        return handle_line_outside_volume();
    }//if( abs(z_exit) > half_length )
  }//if( detector[0] != source[0] ) / else

  // If we are here, we are guaranteed the line does go through our volume

  // If we are exiting the cylinder on the ends, we need to figure out where on those disks we exit
  if( abs(z_exit) > half_length )
  {
    const T z = ((z_exit < 0.0) ? -half_length : half_length);

    // Just solve for x and y as a function of z, and fill those in
    const T m_xz = unit[0] / unit[2];
    const T c_xz = source[0] - m_xz*source[2];

    const T m_yz = unit[1] / unit[2];
    const T c_yz = source[1] - m_yz*source[2];

    x_exit = m_xz*z + c_xz;
    y_exit = m_yz*z + c_yz;
    z_exit = z;
  }//if( abs(z_exit) > half_length )

  const T dx = (source[0] - x_exit);
  const T dy = (source[1] - y_exit);
  const T dz = (source[2] - z_exit);
  const T dist_scaler = sqrt( dx*dx + dy*dy + dz*dz );

  exit_point[0] = x_exit;
  exit_point[1] = y_exit;
  exit_point[2] = z_exit;

  return dist_scaler;
}//cylinder_line_intersection_imp(...)


/** Templated transcription of #rectangle_exit_location.

 `source` must be inside (or on the surface of) the box; returns the distance
 from source to where the ray toward `detector` exits the box.  `exit_point`
 may be the same array as `source`.
 */
template<typename T>
T rectangle_exit_location_imp( const T &half_width, const T &half_height,
                               const T &half_depth,
                               const T source[3],
                               const T detector[3],
                               T exit_point[3] )
{
  using namespace std;
  using namespace ceres;

  T norm[3] = { detector[0] - source[0], detector[1] - source[1], detector[2] - source[2] };

  const T total_dist = sqrt( norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2] );
  norm[0] /= total_dist;
  norm[1] /= total_dist;
  norm[2] /= total_dist;

  // Find the intersection with each of the three candidate exit planes, and take the closest.
  // Line in three-space: {x,y,z} = source + t*norm
  const T inv_slope_x = (norm[0] == 0.0) ? T(DBL_MAX) : (1.0 / norm[0]);
  const T x_intersect = (inv_slope_x >= 0.0) ? half_width : -half_width;
  const T t_intersect_x = (x_intersect - source[0])*inv_slope_x;

  const T inv_slope_y = (norm[1] == 0.0) ? T(DBL_MAX) : (1.0 / norm[1]);
  const T y_intersect = (inv_slope_y >= 0.0) ? half_height : -half_height;
  const T t_intersect_y = (y_intersect - source[1])*inv_slope_y;

  const T inv_slope_z = (norm[2] == 0.0) ? T(DBL_MAX) : (1.0 / norm[2]);
  const T z_intersect = (inv_slope_z >= 0.0) ? half_depth : -half_depth;
  const T t_intersect_z = (z_intersect - source[2])*inv_slope_z;

  if( (t_intersect_x <= t_intersect_y) && (t_intersect_x <= t_intersect_z) )
  {
    // We are exiting through the plane perpendicular to x-axis
    exit_point[0] = ((norm[0] >= 0.0) ? half_width : -half_width);
    exit_point[1] = source[1] + (t_intersect_x * norm[1]);
    exit_point[2] = source[2] + (t_intersect_x * norm[2]);

    return t_intersect_x;
  }else if( t_intersect_y <= t_intersect_z )
  {
    // We are exiting through the plane perpendicular to y-axis
    exit_point[0] = source[0] + (t_intersect_y * norm[0]);
    exit_point[1] = ((norm[1] >= 0.0) ? half_height : -half_height);
    exit_point[2] = source[2] + (t_intersect_y * norm[2]);

    return t_intersect_y;
  }else
  {
    // We are exiting through the plane perpendicular to z-axis
    exit_point[0] = source[0] + (t_intersect_z * norm[0]);
    exit_point[1] = source[1] + (t_intersect_z * norm[1]);
    exit_point[2] = ((norm[2] >= 0.0) ? half_depth : -half_depth);

    return t_intersect_z;
  }// if / else figure out where we are exiting.
}//rectangle_exit_location_imp(...)


/** Templated transcription of #rectangle_intersections.

 Both `source` and `detector` must be outside the box; fills in where the ray
 enters and exits the box, returning false if it misses entirely.
 */
template<typename T>
bool rectangle_intersections_imp( const T &half_width, const T &half_height,
                                  const T &half_depth,
                                  const T source[3],
                                  const T detector[3],
                                  T enter_point[3],
                                  T exit_point[3] )
{
  using namespace std;
  using namespace ceres;

  T norm[3] = { detector[0] - source[0], detector[1] - source[1], detector[2] - source[2] };
  const T total_dist = sqrt( norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2] );
  norm[0] /= total_dist;
  norm[1] /= total_dist;
  norm[2] /= total_dist;

  const T inv_slope_x = (norm[0] == 0.0) ? T(DBL_MAX) : (1.0 / norm[0]);
  const T x_intersect = (inv_slope_x >= 0.0) ? half_width : -half_width;
  const T t_intersect_x_det = (x_intersect - source[0])*inv_slope_x;
  const T t_intersect_x_src = (-x_intersect - source[0])*inv_slope_x;

  const T inv_slope_y = (norm[1] == 0.0) ? T(DBL_MAX) : (1.0 / norm[1]);
  const T y_intersect = (inv_slope_y >= 0.0) ? half_height : -half_height;
  const T t_intersect_y_det = (y_intersect - source[1])*inv_slope_y;
  const T t_intersect_y_src = (-y_intersect - source[1])*inv_slope_y;

  const T inv_slope_z = (norm[2] == 0.0) ? T(DBL_MAX) : (1.0 / norm[2]);
  const T z_intersect = (inv_slope_z >= 0.0) ? half_depth : -half_depth;
  const T t_intersect_z_det = (z_intersect - source[2])*inv_slope_z;
  const T t_intersect_z_src = (-z_intersect - source[2])*inv_slope_z;

  const bool intersects_x_src = (t_intersect_x_src >= 0.0);
  const bool intersects_y_src = (t_intersect_y_src >= 0.0);
  const bool intersects_z_src = (t_intersect_z_src >= 0.0);

  const bool x_before_y_src = (!intersects_y_src || (t_intersect_x_src <= t_intersect_y_src));
  const bool x_before_z_src = (!intersects_y_src || (t_intersect_x_src <= t_intersect_z_src));
  const bool y_before_z_src = (!intersects_z_src || (t_intersect_y_src <= t_intersect_z_src));

  if( intersects_x_src && x_before_y_src && x_before_z_src )
  {
    // We are entering through the plane perpendicular to x-axis
    const T src_intersect_x = ((norm[0] >= 0.0) ? -half_width : half_width);
    const T src_intersect_y = source[1] + (t_intersect_x_src * norm[1]);
    const T src_intersect_z = source[2] + (t_intersect_x_src * norm[2]);

    if( (abs(src_intersect_y) > half_height)
       || (abs(src_intersect_z) > half_depth) )
    {
      return false;
    }

    enter_point[0] = src_intersect_x;
    enter_point[1] = src_intersect_y;
    enter_point[2] = src_intersect_z;
  }else if( intersects_y_src && y_before_z_src )
  {
    // We are entering through the plane perpendicular to y-axis
    const T src_intersect_x = source[0] + (t_intersect_y_src * norm[0]);
    const T src_intersect_y = ((norm[1] >= 0.0) ? -half_height : half_height);
    const T src_intersect_z = source[2] + (t_intersect_y_src * norm[2]);

    if( (abs(src_intersect_x) > half_width)
       || (abs(src_intersect_z) > half_depth) )
    {
      return false;
    }

    enter_point[0] = src_intersect_x;
    enter_point[1] = src_intersect_y;
    enter_point[2] = src_intersect_z;
  }else if( intersects_z_src )
  {
    // We are entering through the plane perpendicular to z-axis
    const T src_intersect_x = source[0] + (t_intersect_z_src * norm[0]);
    const T src_intersect_y = source[1] + (t_intersect_z_src * norm[1]);
    const T src_intersect_z = ((norm[2] >= 0.0) ? -half_depth : half_depth);

    if( (abs(src_intersect_x) > half_width)
       || (abs(src_intersect_y) > half_height) )
    {
      return false;
    }

    enter_point[0] = src_intersect_x;
    enter_point[1] = src_intersect_y;
    enter_point[2] = src_intersect_z;
  }else
  {
    return false;
  }// if / else figure out where we are entering.

  const bool intersects_x_det = (t_intersect_x_det >= 0.0);
  const bool intersects_y_det = (t_intersect_y_det >= 0.0);
  const bool intersects_z_det = (t_intersect_z_det >= 0.0);

  const bool x_before_y_det = (!intersects_y_det || (t_intersect_x_det <= t_intersect_y_det));
  const bool x_before_z_det = (!intersects_z_det || (t_intersect_x_det <= t_intersect_z_det));
  const bool y_before_z_det = (!intersects_z_det || (t_intersect_y_det <= t_intersect_z_det));

  if( intersects_x_det && x_before_y_det && x_before_z_det )
  {
    exit_point[0] = ((norm[0] >= 0.0) ? half_width : -half_width);
    exit_point[1] = source[1] + (t_intersect_x_det * norm[1]);
    exit_point[2] = source[2] + (t_intersect_x_det * norm[2]);
  }else if( intersects_y_det && y_before_z_det )
  {
    exit_point[0] = source[0] + (t_intersect_y_det * norm[0]);
    exit_point[1] = ((norm[1] >= 0.0) ? half_height : -half_height);
    exit_point[2] = source[2] + (t_intersect_y_det * norm[2]);
  }else if( intersects_z_det )
  {
    exit_point[0] = source[0] + (t_intersect_z_det * norm[0]);
    exit_point[1] = source[1] + (t_intersect_z_det * norm[1]);
    exit_point[2] = ((norm[2] >= 0.0) ? half_depth : -half_depth);
  }else
  {
    assert( 0 );
  }// if / else figure out where we are exiting.

  return true;
}//rectangle_intersections_imp(...)


/** Templated version of #DistributedSrcCalc - the description of one
 volumetric source (one source nuclide at one gamma energy) and its
 surrounding shells, with everything that depends on fit parameters
 (dimensions, transmission coefficients, activity) carried as T.
 */
template<typename T>
struct DistributedSrcCalcT
{
  using ShellType = DistributedSrcCalc::ShellType;

  struct ShellInfo
  {
    /** Cumulative outer dimensions (not thicknesses); meaning per geometry
     matches DistributedSrcCalc::m_dimensionsTransLenAndType.
     */
    std::array<T,3> dims;

    /** Transmission length coefficient (per length) for Material shells; the
     total (dimensionless) attenuation for zero-extent Generic shells.
     */
    T trans_len_coef;

    ShellType type = ShellType::Material;

    // Per-shell metadata for the effective-AN/AD/H accumulation (double-only
    //  post-fit pass); unused during fitting.
    double density = 0.0;
    double effective_an = 0.0;
    double hydrogen_mass_frac = 0.0;
  };//struct ShellInfo

  GeometryType m_geometry = GeometryType::NumGeometryType;

  /** Index into m_shells of the shell the source is distributed in. */
  size_t m_materialIndex = 0;

  DetectorGeomT<T> m_detector;

  bool m_attenuateForAir = false;
  double m_airTransLenCoef = 0.0;

  bool m_isInSituExponential = false;
  double m_inSituRelaxationLength = -1.0;

  std::vector<ShellInfo> m_shells;

  T m_srcVolumetricActivity;

  /** Gamma energy; informational only. */
  double m_energy = 0.0;
  const SandiaDecay::Nuclide *m_nuclide = nullptr;

  // Outputs of self_shielding_integration_imp
  T integral;
  size_t m_num_evals = 0;
  double m_est_rel_error = 0.0;


  T eval_spherical( const double xx[], const int ndim ) const;
  T eval_single_cyl_end_on( const double xx[], const int ndim ) const;
  T eval_cylinder( const double xx[], const int ndim ) const;
  T eval_rect( const double xx[], const int ndim ) const;
};//struct DistributedSrcCalcT


template<typename T>
T DistributedSrcCalcT<T>::eval_spherical( const double xx[], const int ndim ) const
{
  using namespace std;
  using namespace ceres;

  assert( m_geometry == GeometryType::Spherical );
  assert( m_materialIndex < m_shells.size() );
  assert( m_shells[m_materialIndex].type == ShellType::Material );
  assert( (ndim == 2) || (ndim == 3) );

  const double pi = PhysicalUnits::pi;

  // The sphere shells are rotation-invariant, and (while the detector response
  //  is isotropic - see #detector_response_is_isotropic) only the line-of-sight
  //  distance to the detector matters, so we work in a frame rotated such that
  //  the detector lies on the +z axis at its true distance.
  static_assert( detector_response_is_isotropic(),
                 "Off-axis spheres need general 3D treatment once the detector response is angular" );
  const T obs_dist = sqrt( m_detector.position[0]*m_detector.position[0]
                           + m_detector.position[1]*m_detector.position[1]
                           + m_detector.position[2]*m_detector.position[2] );

  const T source_inner_rad = (m_materialIndex > 0) ? m_shells[m_materialIndex-1].dims[0] : T(0.0);
  const T source_outer_rad = m_shells[m_materialIndex].dims[0];

  //integration coordinates go from 0 to one, so we have to scale the variables
  //r:     goes from source inner radius to outer radius
  //theta: goes from 0 to pi
  //phi:   goes from 0 to 2pi (or is fixed for the 2D azimuthal-symmetry case)
  const double max_theta = pi;
  const double max_phi = 2.0 * pi;

  const double x_r = xx[0];
  const double x_theta = xx[1];
  const double x_phi = (ndim < 3) ? 0.5 : xx[2];

  const T r = source_inner_rad + x_r * (source_outer_rad - source_inner_rad);
  const double theta = x_theta * max_theta;
  const double phi = x_phi * max_phi;

  const T j = (source_outer_rad - source_inner_rad) * max_theta * max_phi;
  const T dV = j * r * r * std::sin(theta);

  T source_point[3];
  source_point[0] = r * std::sin( theta ) * std::cos( phi );
  source_point[1] = r * std::sin( theta ) * std::sin( phi );
  source_point[2] = r * std::cos( theta );

  // Save the source location since 'source_point' will get modified.
  const T x = source_point[0], y = source_point[1], z = source_point[2];

  T trans(0.0);

  {//begin code-block to compute distance through source
    T exit_point[3];
    const T srcRad = m_shells[m_materialIndex].dims[0];
    const T srcTransCoef = m_shells[m_materialIndex].trans_len_coef;
    const T dist_in_src = exit_point_of_sphere_z_imp( source_point, exit_point, srcRad, obs_dist );

    T min_rad(0.0);
    bool needShellCompute = false;
    T inner_shell_point[3] = { T(0.0), T(0.0), T(0.0) };

    if( m_materialIndex > 0 )
    {
      exit_point_of_sphere_z_imp( source_point, inner_shell_point, srcRad, obs_dist, false );

      //Take advantage of symmetry: move the position halfway between the positive and negative
      //  exit points, then check if this point is both in a sub shell, and closer to the
      //  detector than the source point.
      inner_shell_point[0] = 0.5*(inner_shell_point[0] + exit_point[0]);
      inner_shell_point[1] = 0.5*(inner_shell_point[1] + exit_point[1]);
      inner_shell_point[2] = 0.5*(inner_shell_point[2] + exit_point[2]);
      min_rad = sqrt( inner_shell_point[0]*inner_shell_point[0]
                      + inner_shell_point[1]*inner_shell_point[1]
                      + inner_shell_point[2]*inner_shell_point[2] );

      const T subrad = m_shells[m_materialIndex-1].dims[0];
      needShellCompute = ( (min_rad < subrad) && (inner_shell_point[2] > source_point[2]) );
    }//if( m_materialIndex > 0 )

    if( needShellCompute )
    {
      T original_point[3] = { source_point[0], source_point[1], source_point[2] };
      source_point[0] = inner_shell_point[0];
      source_point[1] = inner_shell_point[1];
      source_point[2] = inner_shell_point[2];

      //Calc how far from the gammas original position it was, to the first inner shell it hit
      const T innerRad = m_shells[m_materialIndex-1].dims[0];
      exit_point_of_sphere_z_imp( source_point, exit_point, innerRad, obs_dist, false );
      T srcDist2(0.0);
      for( int i = 0; i < 3; ++i )
      {
        const T diff = exit_point[i] - original_point[i];
        srcDist2 += diff*diff;
      }//for( int i = 0; i < 3; ++i )
      trans += (srcTransCoef * sqrt(srcDist2));

      //find inner most sphere the ray passes through
      size_t start_index = 0;
      while( (start_index < m_shells.size()) && (m_shells[start_index].dims[0] < min_rad) )
        ++start_index;

      assert( m_materialIndex != 0 );
      assert( start_index < m_materialIndex );
      assert( start_index < m_shells.size() );

      if( start_index == m_shells.size() )
        throw std::runtime_error( "Logic error 1 in DistributedSrcCalcT::eval_spherical(...)" );
      if( start_index >= m_materialIndex )
        throw std::runtime_error( "Logic error 2 in DistributedSrcCalcT::eval_spherical(...)" );
      if( m_materialIndex == 0 )
        throw std::runtime_error( "Logic error 3 in DistributedSrcCalcT::eval_spherical(...)" );

      //calc distance it travels through the inner spheres
      for( size_t index = start_index; index <= m_materialIndex; ++index )
      {
        const ShellInfo &shell = m_shells[index];

        switch( shell.type )
        {
          case ShellType::Material:
          {
            const T shellRad = shell.dims[0];
            T dist = exit_point_of_sphere_z_imp( source_point, source_point, shellRad, obs_dist );
            if( index != m_materialIndex )
              dist = 2.0*dist;

            trans += (shell.trans_len_coef * dist);
            break;
          }//case ShellType::Material:

          case ShellType::Generic:
          {
            trans += shell.trans_len_coef;
            break;
          }//case ShellType::Generic:
        }//switch( type )
      }//for( loop over shells the ray dips into )
    }else
    {
      source_point[0] = exit_point[0];
      source_point[1] = exit_point[1];
      source_point[2] = exit_point[2];
      trans += (srcTransCoef * dist_in_src);
    }//if( line actually goes into child sphere ) / else
  }//end codeblock to compute distance through source

  for( size_t i = m_materialIndex+1; i < m_shells.size(); ++i )
  {
    const ShellInfo &shell = m_shells[i];

    switch( shell.type )
    {
      case ShellType::Material:
      {
        const T dist_in_sphere = exit_point_of_sphere_z_imp( source_point, source_point,
                                                             shell.dims[0], obs_dist );
        trans += (shell.trans_len_coef * dist_in_sphere);
        break;
      }//case ShellType::Material:

      case ShellType::Generic:
      {
        trans += shell.trans_len_coef;
        break;
      }//case ShellType::Generic:
    }//switch( type )
  }//for( loop over shells outside the source )

  trans = exp( -trans );

  if( m_attenuateForAir )
  {
    // source_point is now the exit point on the last of the shielding, and the (rotated)
    //  detector is at {0,0,obs_dist}
    const T dx = -source_point[0];
    const T dy = -source_point[1];
    const T dz = source_point[2] - obs_dist;

    const T air_dist = sqrt( dx*dx + dy*dy + dz*dz );

    trans *= exp( -m_airTransLenCoef * air_dist );
  }//if( m_attenuateForAir )

  if( m_isInSituExponential )
  {
    assert( m_inSituRelaxationLength > 0.0 );
    trans *= exp( -(source_outer_rad - r) / m_inSituRelaxationLength );
  }

  {// begin apply detector response
    DetectorGeomT<T> rotated_det = m_detector;
    rotated_det.position[0] = T(0.0);
    rotated_det.position[1] = T(0.0);
    rotated_det.position[2] = obs_dist;

    const T eval_point[3] = { x, y, z };
    trans *= detector_response_factor( rotated_det, eval_point );
  }// end apply detector response

  return trans * dV;
}//eval_spherical(...)


template<typename T>
T DistributedSrcCalcT<T>::eval_single_cyl_end_on( const double xx[], const int ndim ) const
{
  using namespace std;
  using namespace ceres;

  assert( m_materialIndex == 0 );
  assert( m_geometry == GeometryType::CylinderEndOn );
  assert( m_shells.size() == 1 );
  assert( m_shells[0].type == ShellType::Material );
  assert( (ndim == 2) || (ndim == 3) );

  // This fast path additionally relies on the detector being on the cylinder axis.
  assert( scalar_of(m_detector.position[0]) == 0.0 );
  assert( scalar_of(m_detector.position[1]) == 0.0 );

  const T source_outer_rad = m_shells[m_materialIndex].dims[0];
  const T source_half_z = m_shells[m_materialIndex].dims[1];
  const T total_height = 2.0 * source_half_z;
  const T trans_len_coef = m_shells[m_materialIndex].trans_len_coef;
  const T &obs_dist = m_detector.position[2];

  const T source_inner_rad(0.0);

  //integration coordinates go from 0 to one for each dimension:
  //  r:     goes from cylindrical inner radius, to outer radius
  //  theta: goes from 0 to 2pi
  //  z:     goes from the negative half-height to positive half-height
  const double max_theta = 2.0 * PhysicalUnits::pi;

  const double x_r = xx[0];
  const double x_z = xx[ ((ndim==3) ? 2 : 1) ];

  const T r = source_inner_rad + x_r * (source_outer_rad - source_inner_rad);
  const T z = total_height * (x_z - 0.5);

  const T j = (source_outer_rad - source_inner_rad) * max_theta * total_height;
  const T dV = j * r;

  const T eval_z_dist_to_det = obs_dist - z;

  T exit_radius(0.0);
  T trans(0.0);

  {//begin code-block to compute distance through source
    // Take advantage of theta symmetry here
    const T z_dist_in_src = source_half_z - z;
    const T r_dist_in_src = r * z_dist_in_src / eval_z_dist_to_det;
    exit_radius = r - r_dist_in_src;

    const T dist_in_src = sqrt(z_dist_in_src*z_dist_in_src + r_dist_in_src*r_dist_in_src);

    trans += (trans_len_coef * dist_in_src);
  }//end codeblock to compute distance through source

  trans = exp( -trans );

  if( m_attenuateForAir )
  {
    const T dz = obs_dist - source_half_z;

    const T air_dist = sqrt( exit_radius*exit_radius + dz*dz );

    trans *= exp( -m_airTransLenCoef * air_dist );
  }//if( m_attenuateForAir )

  if( m_isInSituExponential )
  {
    assert( m_inSituRelaxationLength > 0.0 );
    trans *= exp( -(source_half_z - z) / m_inSituRelaxationLength );
  }

  // Finally toss in the geometric factor (e.g., 1/r2 from where we are evaluating to detector).
  const T eval_dist_to_det = sqrt( r*r + eval_z_dist_to_det*eval_z_dist_to_det );
  trans *= fractional_solid_angle_imp( 2.0*m_detector.radius, eval_dist_to_det + m_detector.setback );

  return trans * dV;
}//eval_single_cyl_end_on(...)


template<typename T>
T DistributedSrcCalcT<T>::eval_cylinder( const double xx[], const int ndim ) const
{
  using namespace std;
  using namespace ceres;

  assert( m_materialIndex < m_shells.size() );
  assert( (m_geometry == GeometryType::CylinderSideOn)
          || (m_geometry == GeometryType::CylinderEndOn) );
  assert( m_shells[m_materialIndex].type == ShellType::Material );
  assert( ((m_geometry == GeometryType::CylinderEndOn) && ((ndim == 2) || (ndim == 3)))
         || ((m_geometry == GeometryType::CylinderSideOn) && (ndim == 3)) );

  const T source_outer_rad = m_shells[m_materialIndex].dims[0];
  const T source_half_z = m_shells[m_materialIndex].dims[1];
  const T total_height = 2.0 * source_half_z;
  const T trans_len_coef = m_shells[m_materialIndex].trans_len_coef;

  //integration coordinates go from 0 to one for each dimension:
  //  r:     goes from cylindrical inner radius, to outer radius
  //  theta: goes from 0 to 2pi
  //  z:     goes from the negative half-height to positive half-height
  const double max_theta = 2.0 * PhysicalUnits::pi;

  const T r = xx[0] * source_outer_rad;
  const double theta = ((ndim==3) ? xx[1] : 0.0) * max_theta;
  const T z = total_height * (xx[((ndim==3) ? 2 : 1)] - 0.5);

  // If point to evaluate is within inner-cylinder, set the source term to zero and return.
  if( m_materialIndex > 0 )
  {
    const ShellInfo &inner = m_shells[m_materialIndex-1];

    if( (r < inner.dims[0]) && (abs(z) < inner.dims[1]) )
      return T(0.0);
  }//if( m_materialIndex > 0 )

  const T j = source_outer_rad * max_theta * total_height;
  const T dV = j * r;

  const T eval_point[3] = { r * std::cos(theta), r * std::sin(theta), z };
  const T detector_pos[3] = { m_detector.position[0], m_detector.position[1], m_detector.position[2] };

  T exit_point[3];
  const T dist_in_cyl = cylinder_line_intersection_imp( source_outer_rad, source_half_z, eval_point,
                                              detector_pos, CylExitDir::TowardDetector, exit_point );

  T trans(0.0);

  // Do transport through inner cylinders, and also subtract off that distance through source
  //  cylinder
  T inner_distance(0.0);
  for( size_t i = 0; i < m_materialIndex; ++i )
  {
    const ShellInfo &shell = m_shells[i];

    T local_exit_point[3];
    const T dist_to_exit = cylinder_line_intersection_imp( shell.dims[0], shell.dims[1], eval_point,
                                      detector_pos, CylExitDir::TowardDetector, local_exit_point );

    if( dist_to_exit <= 0.0 ) //Line doesnt intersect cylinder
      continue;

    T local_enter_point[3];
    const T dist_to_enter = cylinder_line_intersection_imp( shell.dims[0], shell.dims[1], eval_point,
                                   detector_pos, CylExitDir::AwayFromDetector, local_enter_point );

    const T local_distance = dist_to_exit - dist_to_enter;

    switch( shell.type )
    {
      case ShellType::Material:
        trans += (shell.trans_len_coef * (local_distance - inner_distance));
        break;

      case ShellType::Generic:
        trans += shell.trans_len_coef;
        break;
    }//switch( type )

    inner_distance = local_distance;
  }//for( size_t i = 0; i < m_materialIndex; ++i )

  trans += (trans_len_coef * (dist_in_cyl - inner_distance));

  // Do transport through outer cylinders
  for( size_t i = m_materialIndex + 1; i < m_shells.size(); ++i )
  {
    const ShellInfo &shell = m_shells[i];

    switch( shell.type )
    {
      case ShellType::Generic:
      {
        trans += shell.trans_len_coef;
        break;
      }//case ShellType::Generic:

      case ShellType::Material:
      {
        T outer_exit_point[3];
        const T dist_in_shield = cylinder_line_intersection_imp( shell.dims[0], shell.dims[1],
                                                                 exit_point, detector_pos,
                                                                 CylExitDir::TowardDetector, outer_exit_point );

        trans += (shell.trans_len_coef * dist_in_shield);

        exit_point[0] = outer_exit_point[0];
        exit_point[1] = outer_exit_point[1];
        exit_point[2] = outer_exit_point[2];

        break;
      }//case ShellType::Material:
    }//switch( type )
  }//for( loop over outer shielding )

  trans = exp( -trans );

  if( m_attenuateForAir )
  {
    const T air_dist = distance_imp( exit_point, detector_pos );

    trans *= exp( -m_airTransLenCoef * air_dist );
  }//if( m_attenuateForAir )

  if( m_isInSituExponential )
  {
    assert( m_inSituRelaxationLength > 0.0 );
    if( m_geometry == GeometryType::CylinderSideOn )
      trans *= exp( -(source_outer_rad - r) / m_inSituRelaxationLength );
    else
      trans *= exp( -(source_half_z - z) / m_inSituRelaxationLength );
  }//if( m_isInSituExponential )

  // Finally toss in the geometric factor (e.g., 1/r2 from where we are evaluating to detector).
  trans *= detector_response_factor( m_detector, eval_point );

  return trans * dV;
}//eval_cylinder(...)


template<typename T>
T DistributedSrcCalcT<T>::eval_rect( const double xx[], const int ndim ) const
{
  using namespace std;
  using namespace ceres;

  assert( m_geometry == GeometryType::Rectangular );
  assert( m_materialIndex < m_shells.size() );
  assert( m_shells[m_materialIndex].type == ShellType::Material );
  assert( ndim == 3 );

  const T half_width  = m_shells[m_materialIndex].dims[0];
  const T half_height = m_shells[m_materialIndex].dims[1];
  const T half_depth  = m_shells[m_materialIndex].dims[2];

  const T trans_len_coef = m_shells[m_materialIndex].trans_len_coef;

  // Translate the [0,1.0] integration coordinates to the world coordinates.
  const T eval_x = (xx[0] - 0.5) * 2.0 * half_width;
  const T eval_y = (xx[1] - 0.5) * 2.0 * half_height;
  const T eval_z = (xx[2] - 0.5) * 2.0 * half_depth;

  // Check to see if [eval_x,eval_y,eval_z] is inside an inner volume, and if so the source
  //  term is zero
  if( m_materialIndex > 0 )
  {
    const ShellInfo &inner = m_shells[m_materialIndex-1];
    if( (abs(eval_x) < inner.dims[0]) && (abs(eval_y) < inner.dims[1]) && (abs(eval_z) < inner.dims[2]) )
      return T(0.0);
  }//if( m_materialIndex > 0 )

  const T dV = 8.0 * half_width * half_height * half_depth;

  const T eval_loc[3] = { eval_x, eval_y, eval_z };
  const T detector_loc[3] = { m_detector.position[0], m_detector.position[1], m_detector.position[2] };

  T exit_point[3];
  const T dist_in_src = rectangle_exit_location_imp( half_width, half_height, half_depth,
                                                     eval_loc, detector_loc, exit_point );

  T trans(0.0);

  // Do transport through inner shieldings, and also subtract off that distance through source
  T inner_rect_dist(0.0);
  for( size_t i = 0; i < m_materialIndex; ++i )
  {
    const ShellInfo &shell = m_shells[i];

    T enter_loc[3], exit_loc[3];
    const bool intersects = rectangle_intersections_imp( shell.dims[0], shell.dims[1], shell.dims[2],
                                                         eval_loc, detector_loc, enter_loc, exit_loc );

    if( intersects )
    {
      const T dist = distance_imp( enter_loc, exit_loc );

      switch( shell.type )
      {
        case ShellType::Generic:
        {
          trans += shell.trans_len_coef;
          break;
        }//case ShellType::Generic:

        case ShellType::Material:
        {
          trans += (shell.trans_len_coef * (dist - inner_rect_dist));
          break;
        }//case ShellType::Material:
      }//switch( type )

      inner_rect_dist = dist;
    }//if( intersects )
  }//for( size_t i = 0; i < m_materialIndex; ++i )

  trans += (trans_len_coef * (dist_in_src - inner_rect_dist));

  // Account for additional external shieldings
  for( size_t i = m_materialIndex + 1; i < m_shells.size(); ++i )
  {
    const ShellInfo &shell = m_shells[i];

    switch( shell.type )
    {
      case ShellType::Generic:
      {
        trans += shell.trans_len_coef;
        break;
      }

      case ShellType::Material:
      {
        const T dist_in_shield = rectangle_exit_location_imp( shell.dims[0], shell.dims[1],
                                              shell.dims[2], exit_point, detector_loc, exit_point );

        trans += (shell.trans_len_coef * dist_in_shield);
        break;
      }//case ShellType::Material:
    }//switch( type )
  }//for( loop over outer shieldings )

  trans = exp( -trans );

  if( m_attenuateForAir )
  {
    const T air_dist = distance_imp( exit_point, detector_loc );
    trans *= exp( -m_airTransLenCoef * air_dist );
  }//if( m_attenuateForAir )

  if( m_isInSituExponential )
  {
    assert( m_inSituRelaxationLength > 0.0 );
    trans *= exp( -(half_depth - eval_z) / m_inSituRelaxationLength );
  }

  trans *= detector_response_factor( m_detector, eval_loc );

  return trans * dV;
}//eval_rect(...)


namespace QuadDetail
{
  // Gauss-Legendre nodes and weights mapped onto [0,1] (weights sum to 1).
  //  Same abscissas as boost::math::quadrature::gauss<double,N>, hard-coded to
  //  keep this header self-contained and trivially usable with ceres::Jet.
  static const double gl3_x[3] = { 0.5*(1.0 - 0.774596669241483377), 0.5, 0.5*(1.0 + 0.774596669241483377) };
  static const double gl3_w[3] = { 0.5*0.555555555555555556, 0.5*0.888888888888888889, 0.5*0.555555555555555556 };

  static const double gl5_x[5] = { 0.5*(1.0 - 0.906179845938663993), 0.5*(1.0 - 0.538469310105683091), 0.5,
                                   0.5*(1.0 + 0.538469310105683091), 0.5*(1.0 + 0.906179845938663993) };
  static const double gl5_w[5] = { 0.5*0.236926885056189088, 0.5*0.478628670499366468, 0.5*0.568888888888888889,
                                   0.5*0.478628670499366468, 0.5*0.236926885056189088 };

  /** Tensor-product Gauss-Legendre estimate of the integral of `f` over the
   axis-aligned cell [lower, lower+width] in `ndim` (2 or 3) dimensions.
   `f` takes a const double[3] of coordinates.
   */
  template<typename T, typename Fcn>
  T gl_cell_estimate( const Fcn &f, const int ndim,
                      const double lower[3], const double width[3],
                      const int npts, const double *xs, const double *ws,
                      size_t &num_evals )
  {
    T sum(0.0);
    double xx[3] = { 0.5, 0.5, 0.5 };

    for( int i = 0; i < npts; ++i )
    {
      xx[0] = lower[0] + width[0]*xs[i];
      for( int j = 0; j < npts; ++j )
      {
        xx[1] = lower[1] + width[1]*xs[j];

        if( ndim == 2 )
        {
          sum += (ws[i]*ws[j]) * f( xx );
        }else
        {
          for( int k = 0; k < npts; ++k )
          {
            xx[2] = lower[2] + width[2]*xs[k];
            sum += (ws[i]*ws[j]*ws[k]) * f( xx );
          }
        }//if( ndim == 2 ) / else
      }//for( j )
    }//for( i )

    const double vol = width[0] * width[1] * ((ndim == 3) ? width[2] : 1.0);
    num_evals += static_cast<size_t>( (ndim == 3) ? npts*npts*npts : npts*npts );

    return sum * vol;
  }//gl_cell_estimate(...)


  /** One cell of the adaptive subdivision: its bounds, its GL5 estimate, and
   the embedded GL3/GL5 error estimate that prioritizes refinement.
   */
  template<typename T>
  struct QuadCell
  {
    double err;
    int depth;
    double lower[3];
    double width[3];
    T fine;

    bool operator<( const QuadCell &rhs ) const { return err < rhs.err; }
  };//struct QuadCell

  template<typename T, typename Fcn>
  QuadCell<T> make_quad_cell( const Fcn &f, const int ndim,
                              const double lower[3], const double width[3],
                              const int depth, size_t &num_evals )
  {
    QuadCell<T> cell;
    cell.depth = depth;
    for( int i = 0; i < 3; ++i )
    {
      cell.lower[i] = lower[i];
      cell.width[i] = width[i];
    }

    const T coarse = gl_cell_estimate<T>( f, ndim, lower, width, 3, gl3_x, gl3_w, num_evals );
    cell.fine = gl_cell_estimate<T>( f, ndim, lower, width, 5, gl5_x, gl5_w, num_evals );
    cell.err = std::fabs( scalar_of(cell.fine) - scalar_of(coarse) );

    return cell;
  }//make_quad_cell(...)
}//namespace QuadDetail


/** Globally-adaptive tensor-product Gauss-Legendre integration of `f` over the
 unit cube [0,1]^ndim (ndim = 2 or 3).

 Each cell is estimated with an embedded GL3/GL5 pair; the cell with the
 largest error estimate is repeatedly bisected (in every dimension) until the
 summed error estimate drops below `epsrel` times the current integral
 estimate, or the evaluation budget runs out.  Subdivision decisions use only
 the scalar part of T, so T may be a ceres::Jet<> and the derivative
 components are integrated alongside the value.

 `f` is called as f(const double xx[3]) and must return T (the Jacobian onto
 [0,1]^ndim included, exactly like the Cuhre integrands).

 The error control is relative to the integral magnitude, so wildly
 near-cancelling integrands should not use this; the volumetric-source
 integrands are smooth-ish and positive.
 */
template<typename T, typename Fcn>
T adaptive_unit_cube_integrate( const Fcn &f, const int ndim,
                                const double epsrel,
                                const size_t max_evals,
                                double *est_error = nullptr,
                                size_t *num_evals_out = nullptr )
{
  assert( (ndim == 2) || (ndim == 3) );
  if( (ndim != 2) && (ndim != 3) )
    throw std::runtime_error( "adaptive_unit_cube_integrate: ndim must be 2 or 3" );

  // Cells this small are dominated by floating-point noise; also bounds the
  //  effort spent resolving any single discontinuity.
  const int max_depth = 24;

  // Always refine at least this deep: the embedded GL3/GL5 error estimate can
  //  be deceived on a barely-sampled domain (e.g., a symmetric inner void that
  //  both rules miss the same way).
  const int min_depth = 2;

  size_t num_evals = 0;

  const double root_lower[3] = { 0.0, 0.0, 0.0 };
  const double root_width[3] = { 1.0, 1.0, 1.0 };

  std::vector<QuadDetail::QuadCell<T>> heap;
  heap.push_back( QuadDetail::make_quad_cell<T>( f, ndim, root_lower, root_width, 0, num_evals ) );

  T total = heap.front().fine;
  double active_err = heap.front().err;
  double frozen_err = 0.0;  //error of cells at max_depth, that wont be refined further
  size_t num_below_min_depth = 1;

  while( !heap.empty()
         && ( (num_below_min_depth > 0)
              || ((active_err + frozen_err) > epsrel*std::max(std::fabs(scalar_of(total)), DBL_MIN)) )
         && (num_evals < max_evals) )
  {
    // Pop the cell with the largest error estimate; if any cells havent reached
    //  the minimum depth yet, take those first
    if( num_below_min_depth > 0 )
    {
      // Find a below-min-depth cell (cheap linear scan; only happens for the first few dozen cells)
      size_t index = 0;
      while( (index < heap.size()) && (heap[index].depth >= min_depth) )
        ++index;
      assert( index < heap.size() );
      std::swap( heap[index], heap.back() );
      std::make_heap( heap.begin(), heap.end() - 1 );
      --num_below_min_depth;
    }else
    {
      std::pop_heap( heap.begin(), heap.end() );
    }

    const QuadDetail::QuadCell<T> cell = heap.back();
    heap.pop_back();

    if( cell.depth >= max_depth )
    {
      // Give up on refining this cell; keep its contribution and account its error separately
      active_err -= cell.err;
      frozen_err += cell.err;
      continue;
    }//if( cell.depth >= max_depth )

    total -= cell.fine;
    active_err -= cell.err;

    // Bisect every dimension
    const double half_width[3] = { 0.5*cell.width[0], 0.5*cell.width[1],
                                   (ndim == 3) ? 0.5*cell.width[2] : cell.width[2] };
    const int num_children = (1 << ndim);

    for( int child = 0; child < num_children; ++child )
    {
      const double child_lower[3] = {
        cell.lower[0] + ((child & 0x1) ? half_width[0] : 0.0),
        cell.lower[1] + ((child & 0x2) ? half_width[1] : 0.0),
        cell.lower[2] + (((child & 0x4) && (ndim == 3)) ? half_width[2] : 0.0)
      };

      heap.push_back( QuadDetail::make_quad_cell<T>( f, ndim, child_lower, half_width,
                                                     cell.depth + 1, num_evals ) );
      total += heap.back().fine;
      active_err += heap.back().err;
      if( (cell.depth + 1) < min_depth )
        ++num_below_min_depth;
      std::push_heap( heap.begin(), heap.end() );
    }//for( loop over children )
  }//while( error too large, and budget remains )

  if( est_error )
    *est_error = active_err + frozen_err;
  if( num_evals_out )
    *num_evals_out = num_evals;

  return total;
}//adaptive_unit_cube_integrate(...)


/** Templated equivalent of #ShieldingSourceChi2Fcn::selfShieldingIntegration:
 integrates `calculator` over its source volume, filling
 `calculator.integral` (plus m_num_evals / m_est_rel_error).

 Unlike the Cuhre version, exceptions from the integrand (e.g., inconsistent
 geometry) propagate to the caller, which decides the failure policy.
 */
template<typename T>
void self_shielding_integration_imp( DistributedSrcCalcT<T> &calculator,
                                     const double epsrel = 1.0E-4,
                                     const size_t max_evals = 5000000 )
{
  int ndim = -1;  //the number of dimensions of the integral.

  switch( calculator.m_geometry )
  {
    case GeometryType::Spherical:
    case GeometryType::CylinderEndOn:
      ndim = 2;
      break;

    case GeometryType::CylinderSideOn:
    case GeometryType::Rectangular:
      ndim = 3;
      break;

    case GeometryType::NumGeometryType:
      assert( 0 );
      throw std::runtime_error( "self_shielding_integration_imp: invalid geometry" );
  }//switch( calculator.m_geometry )

  calculator.integral = T(0.0);
  calculator.m_num_evals = 0;
  calculator.m_est_rel_error = 0.0;

  double est_error = 0.0;
  size_t num_evals = 0;

  switch( calculator.m_geometry )
  {
    case GeometryType::Spherical:
    {
      const auto f = [&calculator,ndim]( const double xx[3] ) -> T {
        return calculator.eval_spherical( xx, ndim );
      };
      calculator.integral = adaptive_unit_cube_integrate<T>( f, ndim, epsrel, max_evals, &est_error, &num_evals );
      break;
    }//case GeometryType::Spherical:

    case GeometryType::CylinderEndOn:
    {
      if( (calculator.m_shells.size() == 1)
          && (scalar_of(calculator.m_detector.position[0]) == 0.0)
          && (scalar_of(calculator.m_detector.position[1]) == 0.0) )
      {
        // For a single, on-axis, end-on cylinder we can use the ever-so-slightly faster evaluator
        const auto f = [&calculator,ndim]( const double xx[3] ) -> T {
          return calculator.eval_single_cyl_end_on( xx, ndim );
        };
        calculator.integral = adaptive_unit_cube_integrate<T>( f, ndim, epsrel, max_evals, &est_error, &num_evals );
      }else
      {
        const auto f = [&calculator,ndim]( const double xx[3] ) -> T {
          return calculator.eval_cylinder( xx, ndim );
        };
        calculator.integral = adaptive_unit_cube_integrate<T>( f, ndim, epsrel, max_evals, &est_error, &num_evals );
      }
      break;
    }//case GeometryType::CylinderEndOn:

    case GeometryType::CylinderSideOn:
    {
      const auto f = [&calculator,ndim]( const double xx[3] ) -> T {
        return calculator.eval_cylinder( xx, ndim );
      };
      calculator.integral = adaptive_unit_cube_integrate<T>( f, ndim, epsrel, max_evals, &est_error, &num_evals );
      break;
    }//case GeometryType::CylinderSideOn:

    case GeometryType::Rectangular:
    {
      const auto f = [&calculator,ndim]( const double xx[3] ) -> T {
        return calculator.eval_rect( xx, ndim );
      };
      calculator.integral = adaptive_unit_cube_integrate<T>( f, ndim, epsrel, max_evals, &est_error, &num_evals );
      break;
    }//case GeometryType::Rectangular:

    case GeometryType::NumGeometryType:
      assert( 0 );
      break;
  }//switch( calculator.m_geometry )

  calculator.m_num_evals = num_evals;
  const double integral_scalar = std::fabs( scalar_of(calculator.integral) );
  calculator.m_est_rel_error = (integral_scalar > 0.0) ? (est_error / integral_scalar) : est_error;
}//self_shielding_integration_imp(...)


/** Templated equivalent of #transmition_coefficient_generic: the total
 (dimensionless) attenuation through a generic shielding of the given atomic
 number and areal density, preserving derivative information of both.
 */
template<typename T>
T transmission_coefficient_generic_imp( const T &atomic_number, const T &areal_density,
                                        const float energy )
{
  return areal_density * MassAttenuation::mass_atten_coef_frac_an( atomic_number, energy );
}//transmission_coefficient_generic_imp(...)


// --- ShieldingSourceChi2Fcn templated parameter accessors ---

template<typename T>
T ShieldingSourceChi2Fcn::sphericalThickness_imp( const size_t materialNum, const std::vector<T> &params ) const
{
  if( isGenericMaterial( materialNum ) )
    throw std::runtime_error( "sphericalThickness_imp: cant be called for generic materials" );

  return params.at( 2*m_nuclides.size() + 3*materialNum );
}

template<typename T>
T ShieldingSourceChi2Fcn::cylindricalRadiusThickness_imp( const size_t materialNum, const std::vector<T> &params ) const
{
  assert( !isGenericMaterial( materialNum ) );
  return params.at( 2*m_nuclides.size() + 3*materialNum );
}

template<typename T>
T ShieldingSourceChi2Fcn::cylindricalLengthThickness_imp( const size_t materialNum, const std::vector<T> &params ) const
{
  assert( !isGenericMaterial( materialNum ) );
  return params.at( 2*m_nuclides.size() + 3*materialNum + 1 );
}

template<typename T>
T ShieldingSourceChi2Fcn::rectangularWidthThickness_imp( const size_t materialNum, const std::vector<T> &params ) const
{
  assert( !isGenericMaterial( materialNum ) );
  return params.at( 2*m_nuclides.size() + 3*materialNum );
}

template<typename T>
T ShieldingSourceChi2Fcn::rectangularHeightThickness_imp( const size_t materialNum, const std::vector<T> &params ) const
{
  assert( !isGenericMaterial( materialNum ) );
  return params.at( 2*m_nuclides.size() + 3*materialNum + 1 );
}

template<typename T>
T ShieldingSourceChi2Fcn::rectangularDepthThickness_imp( const size_t materialNum, const std::vector<T> &params ) const
{
  assert( !isGenericMaterial( materialNum ) );
  return params.at( 2*m_nuclides.size() + 3*materialNum + 2 );
}

template<typename T>
T ShieldingSourceChi2Fcn::arealDensity_imp( const size_t materialNum, const std::vector<T> &params ) const
{
  if( !isGenericMaterial( materialNum ) )
    throw std::runtime_error( "arealDensity_imp: only valid for generic materials" );

  return params.at( 2*m_nuclides.size() + 3*materialNum + 1 );
}

template<typename T>
T ShieldingSourceChi2Fcn::atomicNumber_imp( const size_t materialNum, const std::vector<T> &params ) const
{
  if( !isGenericMaterial( materialNum ) )
    throw std::runtime_error( "atomicNumber_imp: only valid for generic materials" );

  return params.at( 2*m_nuclides.size() + 3*materialNum );
}


template<typename T>
T ShieldingSourceChi2Fcn::age_imp( const SandiaDecay::Nuclide *nuclide, const std::vector<T> &params ) const
{
  if( params.size() != numExpectedFitParameters() )
    throw std::runtime_error( "age_imp(...) invalid params size" );

  const size_t ind = nuclideIndex( nuclide );

  //If parameter is zero or positive, we will use this age.
  const T &par = params[2*ind + 1];
  if( scalar_of(par) >= -0.00001 )
    return (scalar_of(par) > 0.0) ? par : T(0.0);

  //Else the parameter value indicates the "defining" nuclide whose age should be
  // specified (a negative, constant, parameter - so no derivative info lost here).
  const double findex = -1.0*scalar_of(par);
  const int nearFIndex = static_cast<int>( std::round(findex) );

  if( (std::fabs(findex - nearFIndex) > 0.01) || (nearFIndex < 1) )
    throw std::runtime_error( "age_imp: invalid age-defining-nuclide parameter" );

  const size_t defining_nuc_age_index = 2*(nearFIndex - 1) + 1;
  if( defining_nuc_age_index >= params.size() )
    throw std::runtime_error( "age_imp: age-defining-nuclide index to large" );

  const T &defining_age = params[defining_nuc_age_index];
  if( scalar_of(defining_age) < -0.00001 )
    throw std::runtime_error( "age_imp: defining age is also less than zero (shouldnt happen)" );

  return (scalar_of(defining_age) > 0.0) ? defining_age : T(0.0);
}//age_imp(...)


template<typename T>
T ShieldingSourceChi2Fcn::volumeOfMaterial_imp( const size_t matn, const std::vector<T> &params ) const
{
  if( matn >= m_initial_shieldings.size() )
    throw std::runtime_error( "volumeOfMaterial_imp: invalid material index" );

  if( !m_initial_shieldings[matn].m_material )
    throw std::runtime_error( "volumeOfMaterial_imp: cant be called for generic material" );

  T inner_dim_1(0.0), inner_dim_2(0.0), inner_dim_3(0.0);
  T outer_dim_1(0.0), outer_dim_2(0.0), outer_dim_3(0.0);

  for( size_t index = 0; index <= matn; ++index )
  {
    if( !m_initial_shieldings[index].m_material )  //if a generic shielding
      continue;

    inner_dim_1 = outer_dim_1;
    inner_dim_2 = outer_dim_2;
    inner_dim_3 = outer_dim_3;

    switch( m_geometry )
    {
      case GeometryType::Spherical:
        outer_dim_1 += sphericalThickness_imp( index, params );
        break;

      case GeometryType::CylinderEndOn:
      case GeometryType::CylinderSideOn:
        outer_dim_1 += cylindricalRadiusThickness_imp( index, params );
        outer_dim_2 += cylindricalLengthThickness_imp( index, params );
        break;

      case GeometryType::Rectangular:
        outer_dim_1 += rectangularWidthThickness_imp( index, params );
        outer_dim_2 += rectangularHeightThickness_imp( index, params );
        outer_dim_3 += rectangularDepthThickness_imp( index, params );
        break;

      case GeometryType::NumGeometryType:
        assert( 0 );
        break;
    }//switch( m_geometry )
  }//for( int index = 0; index < matn; ++index )

  switch( m_geometry )
  {
    case GeometryType::Spherical:
      return (4.0/3.0)*PhysicalUnits::pi*( outer_dim_1*outer_dim_1*outer_dim_1
                                           - inner_dim_1*inner_dim_1*inner_dim_1 );

    case GeometryType::CylinderEndOn:
    case GeometryType::CylinderSideOn:
    {
      const T inner_vol = PhysicalUnits::pi * inner_dim_1 * inner_dim_1 * 2.0 * inner_dim_2;
      const T outer_vol = PhysicalUnits::pi * outer_dim_1 * outer_dim_1 * 2.0 * outer_dim_2;

      return outer_vol - inner_vol;
    }

    case GeometryType::Rectangular:
    {
      const T inner_vol = 8.0 * inner_dim_1 * inner_dim_2 * inner_dim_3;
      const T outer_vol = 8.0 * outer_dim_1 * outer_dim_2 * outer_dim_3;

      return outer_vol - inner_vol;
    }

    case GeometryType::NumGeometryType:
      assert( 0 );
      break;
  }//switch( m_geometry )

  assert( 0 );
  return T(0.0);
}//volumeOfMaterial_imp(...)


template<typename T>
T ShieldingSourceChi2Fcn::massFractionOfElement_imp( const size_t material_index,
                                                     const SandiaDecay::Nuclide *nuc,
                                                     const std::vector<T> &pars ) const
{
  // Value-only transcription of the double massFractionOfElement(...) - see that
  //  function for the full commentary on the parameter layout.
  assert( nuc );
  if( !nuc )
    throw std::runtime_error( "massFractionOfElement_imp: nullptr nuc" );

  if( material_index >= m_initial_shieldings.size() )
    throw std::logic_error( "massFractionOfElement_imp: invalid material index" );

  const ShieldingSourceFitCalc::ShieldingInfo &shielding = m_initial_shieldings[material_index];

  if( !shielding.m_material )
    throw std::runtime_error( "massFractionOfElement_imp: invalid input" );

  //If we arent fitting this mass fraction, return its initial (constant) fraction
  const short int atomic_num = nuc->atomicNumber;
  for( const auto &el_nucs : shielding.m_nuclideFractions_ )
  {
    if( !el_nucs.first || (el_nucs.first->atomicNumber != atomic_num) )
      continue;

    for( const std::tuple<const SandiaDecay::Nuclide *,double,bool> &nuc_info : el_nucs.second )
    {
      if( std::get<0>(nuc_info) == nuc )
      {
        if( !std::get<2>(nuc_info) )
          return T( std::get<1>(nuc_info) );
        break;
      }
    }//for( loop over nuclide information for this Element )
  }//for( loop over to check if we are actually fitting this mass fraction )

  //Count mass-fraction parameters from materials before this one
  size_t num_pre_el_coefs = 0;
  for( size_t pre_mat_index = 0; pre_mat_index < material_index; ++pre_mat_index )
  {
    const ShieldingSourceFitCalc::ShieldingInfo &pre_shield = m_initial_shieldings[pre_mat_index];
    if( !pre_shield.m_material ) //skip if generic material
      continue;

    for( const auto &el_nucs : pre_shield.m_nuclideFractions_ )
    {
      size_t num_fit_this_el = 0;
      for( const auto &nuc_info : el_nucs.second )
        num_fit_this_el += std::get<2>(nuc_info);

      if( num_fit_this_el > 1 )
        num_pre_el_coefs += (num_fit_this_el - 1);
    }//for( loop over nuclide fractions )
  }//for( loop over materials before current one )

  const std::tuple<const SandiaDecay::Nuclide *,double,bool> *initial_nuc_info = nullptr;
  const std::vector<std::tuple<const SandiaDecay::Nuclide *,double,bool>> *initial_el_info = nullptr;

  size_t nuc_index_in_el = 0, num_variable_nucs_in_el = 0;

  for( const auto &el_nucs : shielding.m_nuclideFractions_ )
  {
    const SandiaDecay::Element * const this_el = el_nucs.first;
    if( !this_el )
      throw std::logic_error( "massFractionOfElement_imp: nullptr element" );

    if( this_el->atomicNumber != atomic_num )
    {
      size_t num_fit_fracs = 0;
      for( const auto &nuc_info : el_nucs.second )
        num_fit_fracs += std::get<2>(nuc_info);

      if( num_fit_fracs > 0 )
        num_pre_el_coefs += (num_fit_fracs - 1);

      continue;
    }//if( not the element of interest )

    initial_el_info = &el_nucs.second;

    for( const auto &nuc_info : el_nucs.second )
    {
      num_variable_nucs_in_el += std::get<2>(nuc_info);

      if( !initial_nuc_info )
      {
        if( std::get<0>(nuc_info) == nuc )
          initial_nuc_info = &nuc_info;
        else
          nuc_index_in_el += std::get<2>(nuc_info);
      }//if( !initial_nuc_info )
    }//for( const auto &nuc_info : el_nucs.second )

    break;
  }//for( const auto &el_nucs : self_atten_nucs )

  if( !initial_el_info || !initial_nuc_info )
    throw std::runtime_error( "massFractionOfElement_imp: " + nuc->symbol
                              + " was not a self-attenuating nuclide in "
                              + shielding.m_material->name );

  if( num_variable_nucs_in_el <= 1 )
    throw std::logic_error( "massFractionOfElement_imp: invalid number of variable nucs for an element" );

  // If we arent fitting mass-fraction, we can return here.
  if( !std::get<2>(*initial_nuc_info) )
    return T( std::get<1>(*initial_nuc_info) );

  double totalfrac = 0.0;  //The total mass fraction of this element being fit
  for( const std::tuple<const SandiaDecay::Nuclide *,double,bool> &i : *initial_el_info )
  {
    if( std::get<2>(i) )
      totalfrac += std::get<1>(i);
  }//for( loop over nuclides of this element )

  size_t matmassfracstart = 2*m_nuclides.size() + 3*m_initial_shieldings.size() + num_pre_el_coefs;

  const size_t numfraccoefs = num_variable_nucs_in_el - 1;
  const size_t thismatmassfrac = matmassfracstart + nuc_index_in_el;

  T frac(1.0);
  for( size_t index = matmassfracstart; index < thismatmassfrac; ++index )
    frac *= (1.0 - pars.at(index));

  if( nuc_index_in_el != numfraccoefs )
    frac *= pars.at(thismatmassfrac);

  const T massFrac = totalfrac * frac;

  if( (std::isnan)(scalar_of(massFrac)) || (std::isinf)(scalar_of(massFrac)) )
    throw std::runtime_error( "massFractionOfElement_imp: invalid parameters" );

  return massFrac;
}//massFractionOfElement_imp(...)


template<typename T>
T ShieldingSourceChi2Fcn::activityOfSelfAttenSource_imp( const SandiaDecay::Nuclide *nuclide,
                                                         const std::vector<T> &params ) const
{
  assert( nuclide );
  if( !nuclide )
    throw std::runtime_error( "activityOfSelfAttenSource_imp: called with null nuclide" );

  assert( isSelfAttenSource(nuclide) );

  bool foundSrc = false;
  T activity(0.0);

  const size_t num_materials = m_initial_shieldings.size();
  for( size_t material_index = 0; material_index < num_materials; ++material_index )
  {
    const ShieldingSourceFitCalc::ShieldingInfo &shield = m_initial_shieldings[material_index];

    const std::shared_ptr<const Material> &mat = shield.m_material;
    if( !mat )  //if a generic shielding
      continue;

    // Determine if nuclide is a self-attenuating source of this shielding
    const std::tuple<const SandiaDecay::Nuclide *,double,bool> *self_att_pos = nullptr;
    for( const auto &el_nucs : shield.m_nuclideFractions_ )
    {
      if( !el_nucs.first || (el_nucs.first->atomicNumber != nuclide->atomicNumber) )
        continue;

      for( const auto &nuc_info : el_nucs.second )
      {
        if( std::get<0>(nuc_info) == nuclide )
        {
          self_att_pos = &nuc_info;
          break;
        }
      }//for( const auto &nuc_info : el_nucs.second )
    }//for( const auto &el_nucs : shield.m_nuclideFractions_ )

    if( !self_att_pos )
      continue;

    foundSrc = true;

    T massFrac(0.0);
    if( std::get<2>(*self_att_pos) )
      massFrac = massFractionOfElement_imp( material_index, nuclide, params );
    else
      massFrac = T( std::get<1>(*self_att_pos) );

    const T vol = volumeOfMaterial_imp( material_index, params );
    const T mass_grams = massFrac * static_cast<double>(mat->density) * vol / PhysicalUnits::gram;
    activity += mass_grams * nuclide->activityPerGram();
  }//for( loop over materials )

  if( !foundSrc )
    throw std::runtime_error( "activityOfSelfAttenSource_imp: " + nuclide->symbol
                              + " is not a self-attenuating source" );

  return activity;
}//activityOfSelfAttenSource_imp(...)


template<typename T>
T ShieldingSourceChi2Fcn::activity_imp( const SandiaDecay::Nuclide *nuclide, const std::vector<T> &params ) const
{
  if( params.size() != numExpectedFitParameters() )
    throw std::runtime_error( "activity_imp(...) invalid params size" );

  if( isSelfAttenSource(nuclide) )
    return activityOfSelfAttenSource_imp( nuclide, params );

  const size_t ind = nuclideIndex( nuclide );
  return params[2*ind] * sm_activityUnits;
}//activity_imp(...)


template<typename T>
void ShieldingSourceChi2Fcn::cluster_peak_activities_imp( std::map<double,T> &energy_count_map,
                  const std::vector<std::pair<double,double>> &energie_widths,
                  SandiaDecay::NuclideMixture &mixture,
                  const T &act, const T &age,
                  const double photopeakClusterSigma,
                  const double energyToCluster,
                  const bool accountForDecayDuringMeas,
                  const double measDuration )
{
  if( energie_widths.empty() )
    return;

  if( energy_count_map.empty() )
  {
    for( const std::pair<double,double> &dp : energie_widths )
      energy_count_map[dp.first] = T(0.0);
  }//if( energy_count_map.empty() )

  if( mixture.numInitialNuclides() != 1 )
    throw std::runtime_error( "cluster_peak_activities_imp: mixture must have exactly one parent nuclide" );

  const SandiaDecay::Nuclide * const nuclide = mixture.initialNuclide(0);
  assert( nuclide );

  const double age_val = scalar_of( age );

  // The gamma yields, and the parent-activity normalization, as a function of age - the
  //  decay calculations are not differentiable, so when the age carries derivative
  //  information (it is being fit), the derivative of each yield with respect to age is
  //  computed by numeric differencing and chained back through the age Jet
  //  (same approach as RelActCalcAuto's peaks_for_energy_range_imp).
  auto gammas_at_age = [&mixture, accountForDecayDuringMeas, measDuration]( const double an_age )
                                                  -> std::vector<SandiaDecay::EnergyRatePair> {
    if( accountForDecayDuringMeas && (measDuration > 0.0)
        && !IsNan(measDuration) && !IsInf(measDuration) )
      return decay_during_meas_corrected_gammas( mixture, an_age, measDuration );
    return mixture.photons( an_age, SandiaDecay::NuclideMixture::OrderByEnergy );
  };//gammas_at_age

  const std::vector<SandiaDecay::EnergyRatePair> gammas = gammas_at_age( age_val );

  //'gammas' have the activity of the original parent decreased by aging by 'age', but we
  //  want the parent to have 'sm_activityUnits' activity at 'age' - so a correction factor.
  const double age_sf = sm_activityUnits / mixture.activity( age_val, nuclide );

  // Decide if we need the age derivative, and if so compute forward/backward differences
  bool age_has_deriv = false;
  if constexpr ( !std::is_same_v<T,double> )
  {
    for( int i = 0; i < T::DIMENSION; ++i )
      age_has_deriv = (age_has_deriv || (age.v[i] != 0.0));
  }

  std::vector<SandiaDecay::EnergyRatePair> gammas_fwd, gammas_bwd;
  double forward_dt = -1.0, backward_dt = -1.0;
  double age_sf_fwd = 0.0, age_sf_bwd = 0.0;

  if( age_has_deriv )
  {
    // Step size heuristics, adapted from RelActCalcAuto: at large ages the yields barely
    //  change, so larger steps avoid a zero derivative from float rounding.
    const double half_life = nuclide->halfLife;
    if( age_val > 0.1*half_life )
    {
      if( age_val > 5.0*half_life )
        forward_dt = backward_dt = 0.1*age_val;
      else if( age_val > 2.5*half_life )
        forward_dt = backward_dt = 0.01*age_val;
      else
        forward_dt = backward_dt = 0.001*age_val;
    }else
    {
      forward_dt = std::min( 0.01*age_val, 0.001*half_life );
      if( forward_dt <= 0.0 )
        forward_dt = 0.001*half_life;  //age is zero
      if( age_val > 0.0 )
        backward_dt = std::min( 0.1*age_val, 0.001*half_life );
    }//if( age_val > 0.1*half_life ) / else

    gammas_fwd = gammas_at_age( age_val + forward_dt );
    age_sf_fwd = sm_activityUnits / mixture.activity( age_val + forward_dt, nuclide );
    assert( gammas_fwd.size() == gammas.size() );

    if( backward_dt > 0.0 )
    {
      gammas_bwd = gammas_at_age( age_val - backward_dt );
      age_sf_bwd = sm_activityUnits / mixture.activity( age_val - backward_dt, nuclide );
      assert( gammas_bwd.size() == gammas.size() );
    }
  }//if( age_has_deriv )

  for( size_t gamma_index = 0; gamma_index < gammas.size(); ++gamma_index )
  {
    const SandiaDecay::EnergyRatePair &aep = gammas[gamma_index];

    // Find which fit-peak (if any) this gamma line clusters into - identical branch
    //  structure to the double cluster_peak_activities(...)
    const std::pair<double,double> epair( aep.energy, 0.0 );
    std::vector<std::pair<double,double>>::const_iterator epos
        = std::lower_bound( energie_widths.begin(), energie_widths.end(), epair,
            []( const std::pair<double,double> &lhs, const std::pair<double,double> &rhs )
              { return lhs.first < rhs.first; } );

    if( (photopeakClusterSigma <= 0.0)
        && ((epos == energie_widths.end()) || (epos->first != aep.energy)) )
      continue;

    if( epos == energie_widths.end() )
    {
      //aep.energy is larger than the largest fit peak energy
      if( (aep.energy - energie_widths.back().first)
          > (energie_widths.back().second*photopeakClusterSigma) )
        continue;
      epos--;
    }else if( (epos == energie_widths.begin()) && (epos->first != aep.energy) )
    {
      //aep.energy is less than the smallest fit peak energy
      if( (epos->first - aep.energy) > (epos->second*photopeakClusterSigma) )
        continue;
    }else if( epos->first != aep.energy )
    {
      //see if the nearest peaks are close enough; if so, assign to closest peak
      const double next_width = epos->second;
      const double next_energy = epos->first;
      const double prev_width = (epos-1)->second;
      const double prev_energy = (epos-1)->first;

      const double prev_sigma = (aep.energy - prev_energy) / prev_width;
      const double next_sigma = (next_energy - aep.energy) / next_width;

      if( (prev_sigma < next_sigma) && (prev_sigma < photopeakClusterSigma) )
      {
        epos--;
      }else if( next_sigma < photopeakClusterSigma )
      {
        //epos is already the correct value
      }else
      {
        continue;
      }
    }//if / else

    if( energyToCluster > 0.0 )
    {
      const double dist = std::fabs( (energyToCluster - epos->first) / epos->second );
      if( dist > photopeakClusterSigma )
        continue;
    }//if( energyToCluster > 0.0 )

    const double energy = epos->first;

    if( energy_count_map.count( energy ) == 0 )
      throw std::runtime_error( "cluster_peak_activities_imp: logic error - energy not in map" );

    // The age-dependent part of the contribution: yield(age) * age_sf(age)
    const double f_val = aep.numPerSecond * age_sf;

    T f_jet( f_val );

    if( age_has_deriv )
    {
      if constexpr ( !std::is_same_v<T,double> )
      {
        assert( gamma_index < gammas_fwd.size() );
        assert( std::fabs(gammas_fwd[gamma_index].energy - aep.energy) < 1.0E-6 );
        const double f_plus = gammas_fwd[gamma_index].numPerSecond * age_sf_fwd;

        double df_dage;
        if( backward_dt > 0.0 )
        {
          assert( gamma_index < gammas_bwd.size() );
          assert( std::fabs(gammas_bwd[gamma_index].energy - aep.energy) < 1.0E-6 );
          const double f_minus = gammas_bwd[gamma_index].numPerSecond * age_sf_bwd;

          // Second-order central difference, valid for unequal step sizes
          const double hf = forward_dt, hb = backward_dt;
          df_dage = (hb*hb*f_plus - hf*hf*f_minus + (hf*hf - hb*hb)*f_val)
                    / (hf*hb*(hf + hb));
        }else
        {
          df_dage = (f_plus - f_val) / forward_dt;
        }

        f_jet.v = df_dage * age.v;  //chain rule back through however age depends on the fit parameters
      }//if constexpr ( !std::is_same_v<T,double> )
    }//if( age_has_deriv )

    const T contribution = act * f_jet / sm_activityUnits;
    energy_count_map[energy] += contribution;
  }//for( loop over gammas )
}//cluster_peak_activities_imp(...)


template<typename T>
std::vector<T> ShieldingSourceChi2Fcn::expected_peak_counts_imp( const std::vector<T> &params,
                                                          NucMixtureCache &mixturecache ) const
{
  using namespace std;
  using namespace ceres;

  typedef std::map<double,T> EnergyCountMapT;

  assert( !m_options.attenuate_for_air || !m_detector || !m_detector->isFixedGeometry() );

  if( params.size() != numExpectedFitParameters() )
    throw std::runtime_error( "expected_peak_counts_imp: invalid params size" );

  EnergyCountMapT energy_count_map;
  const std::vector<std::pair<double,double>> energie_widths = observedPeakEnergyWidths( m_peaks );

  // 1) Un-attenuated contributions of the point sources
  if( m_options.multiple_nucs_contribute_to_peaks )
  {
    for( const SandiaDecay::Nuclide *nuclide : m_nuclides )
    {
      if( isVolumetricSource(nuclide) )
        continue;

      const T act = activity_imp( nuclide, params );
      const T thisage = age_imp( nuclide, params );

      if( mixturecache.find(nuclide) == mixturecache.end() )
        mixturecache[nuclide].addNuclideByActivity( nuclide, sm_activityUnits );

      cluster_peak_activities_imp( energy_count_map, energie_widths, mixturecache[nuclide],
                                   act, thisage, m_options.photopeak_cluster_sigma, -1.0,
                                   m_options.account_for_decay_during_meas, m_realTime );
    }//for( const SandiaDecay::Nuclide *nuclide : m_nuclides )
  }else
  {
    for( const PeakDef &peak : m_peaks )
    {
      const SandiaDecay::Nuclide *nuclide = peak.parentNuclide();
      const SandiaDecay::RadParticle *particle = peak.decayParticle();

      if( !nuclide
          || (!particle && (peak.sourceGammaType() != PeakDef::AnnihilationGamma))
          || isVolumetricSource(nuclide) )
        continue;

      const T act = activity_imp( nuclide, params );
      const T thisage = age_imp( nuclide, params );

      if( mixturecache.find(nuclide) == mixturecache.end() )
        mixturecache[nuclide].addNuclideByActivity( nuclide, sm_activityUnits );

      const float energy = peak.gammaParticleEnergy();
      cluster_peak_activities_imp( energy_count_map, energie_widths, mixturecache[nuclide],
                                   act, thisage, m_options.photopeak_cluster_sigma, energy,
                                   m_options.account_for_decay_during_meas, m_realTime );
    }//for( const PeakDef &peak : m_peaks )
  }//if( m_options.multiple_nucs_contribute_to_peaks ) / else

  // 2) Propagate the point-source gammas through each shielding
  T shield_outer_rad(0.0);
  const size_t nMaterials = m_initial_shieldings.size();

  for( size_t materialN = 0; materialN < nMaterials; ++materialN )
  {
    std::function<T(float)> att_coef_fcn;
    const ShieldingSourceFitCalc::ShieldingInfo &shielding = m_initial_shieldings[materialN];
    const std::shared_ptr<const Material> &material = shielding.m_material;

    if( !material )
    {
      // Generic material here.
      T atomic_number = atomicNumber_imp( materialN, params );
      T areal_density = arealDensity_imp( materialN, params );

      // Clamp way-out-of-bounds values (the optimizer should keep things in bounds, but the
      //  Minuit path historically hasnt always) - only the scalar part is clamped, so
      //  derivative information survives.
      if( scalar_of(atomic_number) < MassAttenuation::sm_min_xs_atomic_number )
        atomic_number = T( static_cast<double>(MassAttenuation::sm_min_xs_atomic_number) );
      if( scalar_of(atomic_number) > MassAttenuation::sm_max_xs_atomic_number )
        atomic_number = T( static_cast<double>(MassAttenuation::sm_max_xs_atomic_number) );

      const double ad_in_gcm2 = scalar_of(areal_density) * PhysicalUnits::cm2 / PhysicalUnits::g;
      if( ad_in_gcm2 < 0.0 )
        areal_density = T(0.0);
      if( ad_in_gcm2 > sm_max_areal_density_g_cm2 )
        areal_density = T( sm_max_areal_density_g_cm2 * PhysicalUnits::g / PhysicalUnits::cm2 );

      att_coef_fcn = [atomic_number, areal_density]( float energy ) -> T {
        return transmission_coefficient_generic_imp( atomic_number, areal_density, energy );
      };
    }else
    {
      T thickness(0.0);
      switch( m_geometry )
      {
        case GeometryType::Spherical:
          thickness = sphericalThickness_imp( materialN, params );
          break;

        case GeometryType::CylinderEndOn:
          thickness = cylindricalLengthThickness_imp( materialN, params );
          break;

        case GeometryType::CylinderSideOn:
          thickness = cylindricalRadiusThickness_imp( materialN, params );
          break;

        case GeometryType::Rectangular:
          thickness = rectangularDepthThickness_imp( materialN, params );
          break;

        case GeometryType::NumGeometryType:
          assert( 0 );
          break;
      }//switch( m_geometry )

      shield_outer_rad += thickness;

      att_coef_fcn = [mat = material.get(), thickness]( float energy ) -> T {
        return thickness * transmition_length_coefficient( mat, energy );
      };
    }//if( generic material ) / else

    for( typename EnergyCountMapT::value_type &energy_count : energy_count_map )
      energy_count.second *= exp( -1.0 * att_coef_fcn( static_cast<float>(energy_count.first) ) );
  }//for( size_t materialN = 0; materialN < nMaterials; ++materialN )

  // 3) Attenuation in the air between the shielding surface and the detector
  if( m_options.attenuate_for_air && (!m_detector || !m_detector->isFixedGeometry()) )
  {
    T air_dist = m_distance - shield_outer_rad;
    if( scalar_of(air_dist) < 0.0 )
      air_dist = T(0.0);

    for( typename EnergyCountMapT::value_type &energy_count : energy_count_map )
    {
      const double coef = transmission_length_coefficient_air( static_cast<float>(energy_count.first) );
      energy_count.second *= exp( -coef * air_dist );
    }
  }//if( m_options.attenuate_for_air )

  // 4) Detector response (at the fixed user-entered distance)
  if( m_detector && m_detector->isValid() )
  {
    const bool fixed_geom = m_detector->isFixedGeometry();
    for( typename EnergyCountMapT::value_type &energy_count : energy_count_map )
    {
      const double eff = fixed_geom
                  ? m_detector->intrinsicEfficiency( static_cast<float>(energy_count.first) )
                  : m_detector->efficiency( static_cast<float>(energy_count.first), m_distance );
      energy_count.second *= eff;
    }
  }else
  {
    const double detDiam = 1.0 * PhysicalUnits::cm;
    const double r = 0.5 * detDiam;
    const double D = m_distance;
    const double fracAngle = 0.5*(1.0 - (D/std::sqrt(D*D + r*r)));
    for( typename EnergyCountMapT::value_type &energy_count : energy_count_map )
      energy_count.second *= fracAngle;
  }//if( m_detector && m_detector->isValid() ) / else

  // 5) Contributions from self-attenuating and trace sources
  std::vector<std::unique_ptr<DistributedSrcCalcT<T>>> calculators;

  for( size_t material_index = 0; material_index < nMaterials; ++material_index )
  {
    const ShieldingSourceFitCalc::ShieldingInfo &shield = m_initial_shieldings[material_index];
    const std::shared_ptr<const Material> material = shield.m_material;

    if( !material )
      continue;

    std::vector<const SandiaDecay::Nuclide *> combined_srcs;
    const std::vector<ShieldingSourceFitCalc::TraceSourceInfo> &trace_srcs = shield.m_traceSources;
    for( const ShieldingSourceFitCalc::TraceSourceInfo &p : trace_srcs )
      combined_srcs.push_back( p.m_nuclide );

    for( const auto &el_nucs : shield.m_nuclideFractions_ )
    {
      for( const std::tuple<const SandiaDecay::Nuclide *,double,bool> &nuc : el_nucs.second )
      {
        if( std::get<0>(nuc) )
          combined_srcs.push_back( std::get<0>(nuc) );
      }
    }//for( const auto &el_nucs : shield.m_nuclideFractions_ )

    if( combined_srcs.empty() )
      continue;

    if( m_detector && m_detector->isFixedGeometry() )
      throw std::logic_error( "Self-attenuating and trace sources are not allowed for"
                              " fixed geometry detector response functions." );

    DistributedSrcCalcT<T> baseCalculator;
    baseCalculator.m_geometry = m_geometry;
    baseCalculator.m_materialIndex = material_index;
    baseCalculator.m_attenuateForAir = m_options.attenuate_for_air;
    baseCalculator.m_isInSituExponential = false;
    baseCalculator.m_inSituRelaxationLength = -1.0;

    const double det_radius = (m_detector ? 0.5*m_detector->detectorDiameter() : 0.5*PhysicalUnits::cm);
    const double det_setback = (m_detector ? m_detector->detectorSetback() : 0.0);
    baseCalculator.m_detector = detector_geom_from_config( m_geometry, T(m_distance),
                                                           det_radius, det_setback );

    for( size_t src_index = 0; src_index < combined_srcs.size(); ++src_index )
    {
      const SandiaDecay::Nuclide *src = combined_srcs[src_index];
      const bool is_trace = (src_index < trace_srcs.size());

      T actPerVol(0.0);

      if( is_trace )
      {
        const T act = activity_imp( src, params );

        switch( trace_srcs[src_index].m_type )
        {
          case TraceActivityType::TotalActivity:
            actPerVol = act / volumeOfMaterial_imp( material_index, params );
            break;

          case TraceActivityType::ActivityPerCm3:
            actPerVol = act / PhysicalUnits::cm3;
            break;

          case TraceActivityType::ActivityPerGram:
            actPerVol = act * static_cast<double>(material->density) / PhysicalUnits::g;
            break;

          case TraceActivityType::ExponentialDistribution:
          {
            const double relaxation_len = trace_srcs[src_index].m_relaxationDistance;
            assert( relaxation_len > 0.0 );

            //actually activity of the entire soil-column, per surface area
            actPerVol = act / PhysicalUnits::m2;

            // Normalize the activity by integrating over the entire volume
            const double L = relaxation_len;
            switch( m_geometry )
            {
              case GeometryType::Spherical:
              {
                const T R = sphericalThickness_imp( material_index, params );
                const T norm = 4*PhysicalUnits::pi * L * (L*L*(2.0 - 2.0*exp(-R/L)) - 2.0*L*R + R*R);
                actPerVol /= norm;
                break;
              }//case GeometryType::Spherical:

              case GeometryType::Rectangular:
              case GeometryType::CylinderEndOn:
              {
                T R(0.0);
                if( m_geometry == GeometryType::Rectangular )
                  R = 2.0 * rectangularDepthThickness_imp( material_index, params );
                else
                  R = 2.0 * cylindricalLengthThickness_imp( material_index, params );

                const T norm = L * (1.0 - exp(-R / L));
                actPerVol /= norm;
                break;
              }//case Rectangular or CylinderEndOn

              case GeometryType::CylinderSideOn:
              {
                const T R = cylindricalRadiusThickness_imp( material_index, params );
                const T norm = 2.0 * L * PhysicalUnits::pi * (L*(exp(-R/L) - 1.0) + R);
                actPerVol /= norm;
                break;
              }//case GeometryType::CylinderSideOn:

              case GeometryType::NumGeometryType:
                assert( 0 );
                break;
            }//switch( m_geometry )

            break;
          }//case TraceActivityType::ExponentialDistribution:

          case TraceActivityType::NumTraceActivityType:
            assert( 0 );
            throw std::runtime_error( "Invalid trace activity type" );
        }//switch( trace_srcs[src_index].m_type )
      }else
      {
        const double actPerMass = src->activityPerGram() / PhysicalUnits::gram;
        const T massFract = massFractionOfElement_imp( material_index, src, params );
        actPerVol = actPerMass * massFract * static_cast<double>(material->density);
      }//if( is_trace ) / else

      const T thisage = age_imp( src, params );

      if( mixturecache.find(src) == mixturecache.end() )
        mixturecache[src].addNuclideByActivity( src, sm_activityUnits );

      EnergyCountMapT local_energy_count_map;

      if( m_options.multiple_nucs_contribute_to_peaks )
      {
        cluster_peak_activities_imp( local_energy_count_map, energie_widths, mixturecache[src],
                                     actPerVol, thisage, m_options.photopeak_cluster_sigma, -1.0,
                                     m_options.account_for_decay_during_meas, m_realTime );
      }else
      {
        for( const PeakDef &peak : m_peaks )
        {
          if( (peak.parentNuclide() == src)
              && (peak.decayParticle() || (peak.sourceGammaType() == PeakDef::AnnihilationGamma)) )
          {
            cluster_peak_activities_imp( local_energy_count_map, energie_widths, mixturecache[src],
                                         actPerVol, thisage, m_options.photopeak_cluster_sigma,
                                         peak.gammaParticleEnergy(),
                                         m_options.account_for_decay_during_meas, m_realTime );
          }
        }//for( const PeakDef &peak : m_peaks )
      }//if( m_options.multiple_nucs_contribute_to_peaks ) / else

      for( const typename EnergyCountMapT::value_type &energy_count : local_energy_count_map )
      {
        if( scalar_of(energy_count.second) == 0.0 )
          continue;

        auto calculator = std::make_unique<DistributedSrcCalcT<T>>( baseCalculator );

        calculator->m_nuclide = src;
        calculator->m_energy = energy_count.first;
        calculator->m_srcVolumetricActivity = energy_count.second;

        if( m_options.attenuate_for_air )
          calculator->m_airTransLenCoef = transmission_length_coefficient_air( static_cast<float>(energy_count.first) );
        else
          calculator->m_airTransLenCoef = 0.0;

        if( is_trace && (trace_srcs[src_index].m_type == TraceActivityType::ExponentialDistribution) )
        {
          calculator->m_isInSituExponential = true;
          calculator->m_inSituRelaxationLength = trace_srcs[src_index].m_relaxationDistance;
          assert( calculator->m_inSituRelaxationLength > 0.0 );
        }//if( exponentially distributed trace source )

        std::array<T,3> outer_dims = { T(0.0), T(0.0), T(0.0) };

        for( size_t subMat = 0; subMat < nMaterials; ++subMat )
        {
          typename DistributedSrcCalcT<T>::ShellInfo shell;

          if( isGenericMaterial( subMat ) )
          {
            // A generic material at the very center has zero volume, so skip it
            if( subMat == 0 )
              continue;

            const T an = atomicNumber_imp( subMat, params );
            const T ad = arealDensity_imp( subMat, params );

            shell.dims = outer_dims;
            shell.trans_len_coef = transmission_coefficient_generic_imp( an, ad, static_cast<float>(calculator->m_energy) );
            shell.type = DistributedSrcCalc::ShellType::Generic;
            calculator->m_shells.push_back( shell );
            continue;
          }//if( isGenericMaterial( subMat ) )

          switch( m_geometry )
          {
            case GeometryType::Spherical:
              outer_dims[0] += sphericalThickness_imp( subMat, params );
              break;

            case GeometryType::CylinderEndOn:
            case GeometryType::CylinderSideOn:
              outer_dims[0] += cylindricalRadiusThickness_imp( subMat, params );
              outer_dims[1] += cylindricalLengthThickness_imp( subMat, params );
              break;

            case GeometryType::Rectangular:
              outer_dims[0] += rectangularWidthThickness_imp( subMat, params );
              outer_dims[1] += rectangularHeightThickness_imp( subMat, params );
              outer_dims[2] += rectangularDepthThickness_imp( subMat, params );
              break;

            case GeometryType::NumGeometryType:
              assert( 0 );
              break;
          }//switch( m_geometry )

          bool pastDetector = false;
          switch( m_geometry )
          {
            case GeometryType::Spherical:
            case GeometryType::CylinderSideOn:
              pastDetector = (scalar_of(outer_dims[0]) > m_distance);
              break;

            case GeometryType::CylinderEndOn:
              pastDetector = (scalar_of(outer_dims[1]) > m_distance);
              break;

            case GeometryType::Rectangular:
              pastDetector = (scalar_of(outer_dims[2]) > m_distance);
              break;

            case GeometryType::NumGeometryType:
              assert( 0 );
              break;
          }//switch( m_geometry )

          if( pastDetector && !(m_detector && m_detector->isFixedGeometry()) )
            throw std::runtime_error( "expected_peak_counts_imp: radius > distance" );

          const std::shared_ptr<const Material> &sub_material = m_initial_shieldings[subMat].m_material;

          shell.dims = outer_dims;
          shell.trans_len_coef = T( transmition_length_coefficient( sub_material.get(),
                                                  static_cast<float>(calculator->m_energy) ) );
          shell.type = DistributedSrcCalc::ShellType::Material;
          calculator->m_shells.push_back( shell );
        }//for( size_t subMat = 0; subMat < nMaterials; ++subMat )

        if( calculator->m_shells.empty() )
          throw std::logic_error( "No source/shielding sphere for calculator" );

        calculators.push_back( std::move(calculator) );
      }//for( loop over local_energy_count_map )
    }//for( loop over combined_srcs )
  }//for( loop over materials )

  if( !calculators.empty() )
  {
    if( m_options.multithread_self_atten && (calculators.size() > 1) )
    {
      std::mutex error_mutex;
      std::exception_ptr first_error;

      SpecUtilsAsync::ThreadPool pool;
      for( const std::unique_ptr<DistributedSrcCalcT<T>> &calculator : calculators )
      {
        pool.post( [&calculator, &error_mutex, &first_error](){
          try
          {
            self_shielding_integration_imp( *calculator );
          }catch( std::exception & )
          {
            std::lock_guard<std::mutex> lock( error_mutex );
            if( !first_error )
              first_error = std::current_exception();
          }
        } );
      }//for( loop over calculators )
      pool.join();

      if( first_error )
        std::rethrow_exception( first_error );
    }else
    {
      for( const std::unique_ptr<DistributedSrcCalcT<T>> &calculator : calculators )
        self_shielding_integration_imp( *calculator );
    }//if( multithread ) / else

    for( const std::unique_ptr<DistributedSrcCalcT<T>> &calculator : calculators )
    {
      T contrib = calculator->integral * calculator->m_srcVolumetricActivity;

      if( m_detector && m_detector->isValid() )
        contrib *= m_detector->intrinsicEfficiency( static_cast<float>(calculator->m_energy) );

      energy_count_map[calculator->m_energy] += contrib;
    }//for( loop over calculators )
  }//if( !calculators.empty() )

  // 6) Account for live time
  for( typename EnergyCountMapT::value_type &energy_count : energy_count_map )
    energy_count.second *= m_liveTime;

  // 7) Map clustered energies to per-peak expected counts - same peak inclusion and
  //    matching as expected_observed_chis(...)
  std::vector<T> answer;
  answer.reserve( m_peaks.size() );

  for( const PeakDef &peak : m_peaks )
  {
    if( !peak.decayParticle() && (peak.sourceGammaType() != PeakDef::AnnihilationGamma) )
      continue;

    const double energy = peak.gammaParticleEnergy();

    if( energy_count_map.count(energy) == 0 )
      throw std::runtime_error( "expected_peak_counts_imp: logic error - peak energy not in map" );

    T expected_counts(0.0);
    const double sigma = peak.gausPeak() ? peak.sigma() : 0.25*peak.roiWidth();
    for( const typename EnergyCountMapT::value_type &energy_count : energy_count_map )
    {
      if( std::fabs(energy_count.first - energy) < (0.1*sigma) )
        expected_counts += energy_count.second;
    }

    answer.push_back( expected_counts );
  }//for( const PeakDef &peak : m_peaks )

  return answer;
}//expected_peak_counts_imp(...)

}//namespace GammaInteractionCalc

#endif //GammaInteractionCalc_imp_hpp
