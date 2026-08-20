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
#include <tuple>
#include <cstdio>
#include <limits>
#include <random>
#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <algorithm>
#include <stdexcept>

#include "ceres/jet.h"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/Integrate.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/MassAttenuationTool_imp.hpp"

//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE GammaInteractionCalc_suite
//#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


using namespace std;
using namespace boost::unit_test;
using namespace GammaInteractionCalc;


std::string g_test_file_dir;

// We need to set the static data directory, so the code knows where
//  like sandia.decay.xml is located.
void set_data_dir()
{
  // We only need to initialize things once
  static bool s_have_set = false;
  if( s_have_set )
    return;
  
  s_have_set = true;
  
  int argc = boost::unit_test::framework::master_test_suite().argc;
  char **argv = boost::unit_test::framework::master_test_suite().argv;
  
  string datadir;
  
  for( int i = 1; i < argc; ++i )
  {
    cout << "Arg " << i << ": '" << argv[i] << "'" << endl;
    const string arg = argv[i];
    if( SpecUtils::istarts_with( arg, "--datadir=" ) )
      datadir = arg.substr( 10 );
    
    if( SpecUtils::istarts_with( arg, "--testfiledir=" ) )
      g_test_file_dir = arg.substr( 14 );
  }//for( int arg = 1; arg < argc; ++ arg )
  
  SpecUtils::ireplace_all( datadir, "%20", " " );
  SpecUtils::ireplace_all( g_test_file_dir, "%20", " " );
  
  // Search around a little for the data directory, if it wasnt specified
  if( datadir.empty() )
  {
    for( const auto &d : { "data", "../data", "../../data", "../../../data", "/Users/wcjohns/rad_ana/InterSpec/data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) )
      {
        datadir = d;
        break;
      }
    }//for( loop over candidate dirs )
  }//if( datadir.empty() )
  
  const string required_data_file = "findCharacteristics/202204_example_problem_1.n42";
  if( g_test_file_dir.empty() )
  {
    for( const auto &d : { "test_data", "../test_data", "../../test_data", "../target/testing/test_data/", "/Users/wcjohns/rad_ana/InterSpec/target/testing/test_data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, required_data_file) ) )
      {
        g_test_file_dir = d;
        break;
      }
    }//for( loop over candidate dirs )
  }
  
  const string sandia_deacay_file = SpecUtils::append_path(datadir, "sandia.decay.xml");
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( sandia_deacay_file ), "sandia.decay.xml not at '" << sandia_deacay_file << "'" );
  
  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );
  
  // Make sure we can actually init the decay database
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );
  BOOST_REQUIRE_MESSAGE( db->nuclide("U238"), "SandiaDecayDataBase empty?" );
  
  const string required_data_path = SpecUtils::append_path(g_test_file_dir, required_data_file);
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( required_data_path ), "'" << required_data_file << "' not at '" << required_data_path << "' - CWD='" << SpecUtils::get_working_path() << "'"  );
}//void set_data_dir()


namespace
{
  // A simple structure to represent a 2D point.
  struct Point {
    long double x, y;
  };
  
  /**
   * @brief Finds the intersection points of a line and a circle centered at the origin.
   * The points are ordered by their distance from source.
   * @param source The first point defining the line.
   * @param detector The second point defining the line.
   * @param r The radius of the circle.
   * @return A pair of intersection points.
   * @throws std::runtime_error if there are no intersections.
   */
  std::pair<Point, Point> findIntersections(Point source, Point detector, double r) {
    // Line parametrically: P(t) = source + t*(detector-source)
    // x = source.x + t*(detector.x - source.x)
    // y = source.y + t*(detector.y - source.y)
    
    long double dx = detector.x - source.x;
    long double dy = detector.y - source.y;
    
    // Substitute into circle equation x² + y² = r²:
    // (source.x + t*dx)² + (source.y + t*dy)² = r²
    //
    // Expanding and rearranging into quadratic form: at² + bt + c = 0
    long double a = dx*dx + dy*dy;
    long double b = 2.0*(source.x*dx + source.y*dy);
    long double c = source.x*source.x + source.y*source.y - r*r;
    
    // Calculate discriminant
    long double discriminant = b*b - 4.0*a*c;
    
    if (discriminant < 0) {
      throw std::runtime_error("No intersection points found between line and circle");
    }
    
    // Calculate the parameter values for intersection points
    long double sqrt_discriminant = sqrt(discriminant);
    long double t1 = (-b - sqrt_discriminant) / (2.0*a);
    long double t2 = (-b + sqrt_discriminant) / (2.0*a);
    
    // Calculate intersection points
    Point point1 = {source.x + t1*dx, source.y + t1*dy};
    Point point2 = {source.x + t2*dx, source.y + t2*dy};
    
    // Calculate distances from source to order the points
    long double dist1_sq = (point1.x - source.x)*(point1.x - source.x) + (point1.y - source.y)*(point1.y - source.y);
    long double dist2_sq = (point2.x - source.x)*(point2.x - source.x) + (point2.y - source.y)*(point2.y - source.y);
    
    // Return points ordered by distance from source (closest first)
    if (dist1_sq <= dist2_sq) {
      return std::make_pair(point1, point2);
    } else {
      return std::make_pair(point2, point1);
    }
  }

  /** An independent reference implementation of #cylinder_line_intersection, done in long double.

   The cylinder is centered at the origin and oriented along z.  Parameterizing the ray as
   P(t) = source + t*(detector - source), t increases monotonically toward the detector, so the
   t-interval inside the volume is the intersection of the infinite-cylinder interval with the
   |z| <= half_length slab interval; "toward detector" is that interval's high end, and "away from
   detector" its low end.

   Only valid for a source strictly inside the cylinder - which is the case the volumetric-source
   integrands exercise, and the one where the answer is unambiguous.

   @returns Whether the ray intersects the cylinder at all.
   */
  bool reference_cyl_exit( const double radius, const double half_length,
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
}


BOOST_AUTO_TEST_CASE( CylinderLineIntersection )
{
  // This function doesnt exhaustively test #cylinder_line_intersection, but I think it hits all the
  //  lines of its code, so is decent, although there is probably some permutation or possible test
  //  case left out.
  
  const double oneOverSqrt2 = 1.0 / sqrt(2.0);
  
  double dist;
  double radius = 1.0;
  double half_length = 0.5;
  double source_xyz[3], detector_xyz[3], exit_point[3], start_point[3];
  
  
  // Start tests of a point in the cylinder, and going to the detector
  source_xyz[0] = 0.0; source_xyz[1] = 0.0; source_xyz[2] = 0.0;
  detector_xyz[0] = 0.0; detector_xyz[1] = 2.0; detector_xyz[2] = 0.0;
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 0.0), 1.0E-14 );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 1.0), 1.0E-14 );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - 0.0), 1.0E-14 );
  BOOST_CHECK_SMALL( fabs(dist - 1.0), 1.0E-14 );
  
  
  source_xyz[0] = 0.5; source_xyz[1] = 0.5; source_xyz[2] = 0.0;
  detector_xyz[0] = 2.0; detector_xyz[1] = 2.0; detector_xyz[2] = 0.0;
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - oneOverSqrt2), 1.0E-14 );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - oneOverSqrt2), 1.0E-14 );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - 0.0), 1.0E-14 );
  BOOST_CHECK_SMALL( fabs(dist - 0.292893), 1.0E-5 );
  
  
  source_xyz[0] = 0.0; source_xyz[1] = 0.0; source_xyz[2] = 0.0;
  detector_xyz[0] = -2.0; detector_xyz[1] = 0.0; detector_xyz[2] = -2.0;
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - -0.5), 1.0E-14 );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 0), 1.0E-14 );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - -0.5), 1.0E-14 );
  BOOST_CHECK_SMALL( fabs(dist - oneOverSqrt2), 1.0E-14 );
  
  
  source_xyz[0] = 0.5; source_xyz[1] = 0.5; source_xyz[2] = 0.0;
  detector_xyz[0] = 0.5; detector_xyz[1] = 0.5; detector_xyz[2] = 1.0;
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 0.5), 1.0E-14 );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 0.5), 1.0E-14 );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - 0.5), 1.0E-14 );
  BOOST_CHECK_SMALL( fabs(dist - 0.5), 1.0E-14 );
  
  
  source_xyz[0] = 0.5; source_xyz[1] = 0.5; source_xyz[2] = 0.0;
  detector_xyz[0] = 0.0; detector_xyz[1] = 0.0; detector_xyz[2] = 1.0;
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 0.25), 1.0E-14 );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 0.25), 1.0E-14 );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - 0.5), 1.0E-14 );
  BOOST_CHECK_SMALL( fabs(dist - sqrt(2*0.25*0.25 + 0.5*0.5)), 1.0E-14 );
  
  
  radius = 1.0;
  half_length = 1.0;
  source_xyz[0] = 0.0; source_xyz[1] = 0.0; source_xyz[2] = 0.0;
  detector_xyz[0] = 0; detector_xyz[1] = 5; detector_xyz[2] = 5;
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 0.0), 1.0E-14 );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 1.0), 1.0E-14 );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - 1.0), 1.0E-14 );
  BOOST_CHECK_SMALL( fabs(dist - sqrt(2.0)), 1.0E-14 );
  
  
  radius = 1.0;
  half_length = 1.0;
  source_xyz[0] = 0.25; source_xyz[1] = 0.25; source_xyz[2] = 0.25;
  detector_xyz[0] = 0.25; detector_xyz[1] = 0.25; detector_xyz[2] = 5;
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 0.25), 1.0E-14 );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 0.25), 1.0E-14 );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - 1.0), 1.0E-14 );
  BOOST_CHECK_SMALL( fabs(dist - 0.75), 1.0E-14 );
  

  radius = 225000;
  half_length = 1000;
  source_xyz[0] = 112500.0; source_xyz[1] = 0.0; source_xyz[2] = 0.0;
  detector_xyz[0] = 0.0; detector_xyz[1] = 0.0; detector_xyz[2] = 2000;
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 0.5*source_xyz[0]), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - source_xyz[1]), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - half_length), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(dist - sqrt(half_length*half_length + 0.25*source_xyz[0]*source_xyz[0])), 1.0E-12 );
  
  
  // Start tests of a point outside the cylinder, and going through to the detector on other side
  radius = 1;
  half_length = 1;
  source_xyz[0] = 0.0; source_xyz[1] = 0.0; source_xyz[2] = -2.0;
  detector_xyz[0] = 0.0; detector_xyz[1] = 0.0; detector_xyz[2] = 2;
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - source_xyz[0]), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - source_xyz[1]), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - half_length), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(dist - 3.0), 1.0E-12 );
  
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::AwayFromDetector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - source_xyz[0]), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - source_xyz[1]), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[2] + half_length), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(dist - 1.0), 1.0E-12 );
  
  dist = cylinder_line_intersection( radius, half_length, exit_point, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - source_xyz[0]), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - source_xyz[1]), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - half_length), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(dist - 2.0), 1.0E-12 );
  
  
  radius = 1;
  half_length = 1;
  source_xyz[0] = 10.0; source_xyz[1] = 10.0; source_xyz[2] = 0.0;
  detector_xyz[0] = 0.0; detector_xyz[1] = 0.0; detector_xyz[2] = 10.0;
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - source_xyz[0]), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - source_xyz[1]), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - source_xyz[2]), 1.0E-12 );
  BOOST_CHECK_EQUAL( dist, 0.0 );
  
  radius = 1;
  half_length = 1;
  source_xyz[0] = -10.0; source_xyz[1] = 0.0; source_xyz[2] = -10.0;
  detector_xyz[0] = 10.0; detector_xyz[1] = 0.0; detector_xyz[2] = 10;
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 1.0), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 0.0), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - 1.0), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(dist - sqrt(11.0*11.0 + 11.0*11.0)), 1.0E-12 );
  
  
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::AwayFromDetector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - -1.0), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 0.0), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - -1.0), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(dist - sqrt(2*81.0)), 1.0E-12 );
  
  // {-1,0,-1} to {1,0,1}
  dist = cylinder_line_intersection( radius, half_length, exit_point, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 1.0), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 0.0), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - 1.0), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(dist - sqrt(8.0)), 1.0E-12 );
  
  
  radius = 1;
  half_length = 1;
  source_xyz[0] = -5; source_xyz[1] = 5; source_xyz[2] = 0.0;
  detector_xyz[0] = 5; detector_xyz[1] = 4; detector_xyz[2] = 10.0;
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - source_xyz[0]), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - source_xyz[1]), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - source_xyz[2]), 1.0E-12 );
  BOOST_CHECK_EQUAL( dist, 0.0 );
  
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::AwayFromDetector, exit_point );
  BOOST_CHECK_EQUAL( dist, 0.0 );
  
  radius = 1;
  half_length = 1;
  source_xyz[0] = 5; source_xyz[1] = 5; source_xyz[2] = -10.0;
  detector_xyz[0] = 5; detector_xyz[1] = 5; detector_xyz[2] = 10.0;
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - source_xyz[0]), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - source_xyz[1]), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - source_xyz[2]), 1.0E-12 );
  BOOST_CHECK_EQUAL( dist, 0.0 );
  
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::AwayFromDetector, exit_point );
  BOOST_CHECK_EQUAL( dist, 0.0 );
  
  
  radius = 0.0;
  half_length = 0.0;
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - source_xyz[0]), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - source_xyz[1]), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - source_xyz[2]), 1.0E-12 );
  BOOST_CHECK_EQUAL( dist, 0.0 );
  
  
  radius = 1;
  half_length = 1;
  source_xyz[0] = 2; source_xyz[1] = 2; source_xyz[2] = -10.0;
  detector_xyz[0] = 2; detector_xyz[1] = 3; detector_xyz[2] = 10.0;
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - source_xyz[0]), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - source_xyz[1]), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - source_xyz[2]), 1.0E-12 );
  BOOST_CHECK_EQUAL( dist, 0.0 );
  
  
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::AwayFromDetector, exit_point );
  BOOST_CHECK_EQUAL( dist, 0.0 );
  
  
  radius = 1;
  half_length = 1;
  source_xyz[0] = 0.5; source_xyz[1] = 0; source_xyz[2] = 2;
  detector_xyz[0] = 0.5; detector_xyz[1] = 2; detector_xyz[2] = 10.0;
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - source_xyz[0]), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - source_xyz[1]), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - source_xyz[2]), 1.0E-12 );
  BOOST_CHECK_EQUAL( dist, 0.0 );
  
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::AwayFromDetector, exit_point );
  BOOST_CHECK_EQUAL( dist, 0.0 );
  
  
  radius = 1;
  half_length = 1;
  source_xyz[0] = 0.0; source_xyz[1] = 2; source_xyz[2] = -2.0;
  detector_xyz[0] = 0.0; detector_xyz[1] = 1; detector_xyz[2] = 2.0;
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - source_xyz[0]), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - source_xyz[1]), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - source_xyz[2]), 1.0E-12 );
  BOOST_CHECK_EQUAL( dist, 0.0 );
  
  
  // Now do some cases to simulate what will happen when integrating over a volume and there is a
  //  sub-volume
  
  radius = 1;
  half_length = 1;
  source_xyz[0] = 1.5; source_xyz[1] = 0; source_xyz[2] = 0;
  detector_xyz[0] = -1.5; detector_xyz[1] = 0; detector_xyz[2] = 0;
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::AwayFromDetector, start_point );
  BOOST_CHECK_SMALL( fabs(start_point[0] - 1.0), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(start_point[1] - 0), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(start_point[2] - 0), 1.0E-12 );
  BOOST_CHECK_EQUAL( dist, 0.5 );
  
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] + 1.0), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 0), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - 0), 1.0E-12 );
  BOOST_CHECK_EQUAL( dist, 2.5 );
  
  dist = cylinder_line_intersection( radius, half_length, start_point, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( exit_point[0] + 1.0, 1.0E-12 );
  BOOST_CHECK_SMALL( exit_point[1] - 0, 1.0E-12 );
  BOOST_CHECK_SMALL( exit_point[2] - 0, 1.0E-12 );
  BOOST_CHECK_SMALL( dist - 2.0, 1.0E-12 );
  
  
  
  
  radius = 1;
  half_length = 1;
  source_xyz[0] = -1.5; source_xyz[1] = 0; source_xyz[2] = 0;
  detector_xyz[0] = 1.5; detector_xyz[1] = 0; detector_xyz[2] = 0;
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::AwayFromDetector, start_point );
  BOOST_CHECK_SMALL( fabs(start_point[0] + 1.0), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(start_point[1] - 0), 1.0E-12 );
  BOOST_CHECK_SMALL( fabs(start_point[2] - 0), 1.0E-12 );
  BOOST_CHECK_EQUAL( dist, 0.5 );
  
  dist = cylinder_line_intersection( radius, half_length, start_point, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( exit_point[0] - 1.0, 1.0E-12 );
  BOOST_CHECK_SMALL( exit_point[1] - 0, 1.0E-12 );
  BOOST_CHECK_SMALL( exit_point[2] - 0, 1.0E-12 );
  BOOST_CHECK_SMALL( dist - 2.0, 1.0E-12 );
  
  
  // Test on boundaries
  
  // Here we'll make the line pass *just* inside the circle to get around how exact intersections
  //  arent totally worked out yet.
  radius = 1.0;
  half_length = 1;
  source_xyz[0] = 1.0; source_xyz[1] = -1.0; source_xyz[2] = 0;
  // Note: if we set `detector_xyz[0] = (1.0 - 1.0E-6)` - or smaller delta, then tests will fail on Windows, but not macOS - have not yet checked if `findIntersections(...)` or `cylinder_line_intersection(...)` is giving the more accurate answer, but it seems like one (or likely both), are causing some reall numerical roundoffs and such...
  detector_xyz[0] = (1.0 - 5.0E-5); detector_xyz[1] = 1; detector_xyz[2] = 0;
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::TowardDetector, exit_point );
  pair<Point, Point> intersections = findIntersections({source_xyz[0],source_xyz[1]}, {detector_xyz[0],detector_xyz[1]}, radius);
  
  long double tollerance = 1.0E-6;
  BOOST_CHECK_SMALL( exit_point[0] - intersections.second.x, tollerance ); //0.999998998585
  BOOST_CHECK_SMALL( exit_point[1] - intersections.second.y, tollerance ); //0.001415
  BOOST_CHECK_SMALL( exit_point[2] - 0, 1.0E-9 );
  long double test_dist = sqrt(pow( (source_xyz[0] - intersections.second.x), 2.0 ) + pow( (source_xyz[1] - intersections.second.y), 2.0 ));
  BOOST_CHECK_SMALL( dist - test_dist, tollerance );
  
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::AwayFromDetector, exit_point );
  BOOST_CHECK_SMALL( exit_point[0] - intersections.first.x, tollerance ); //0.999999001413
  BOOST_CHECK_SMALL( exit_point[1] - intersections.first.y, tollerance ); //-0.001413
  BOOST_CHECK_SMALL( exit_point[2] - 0, 1.0E-9 );
  test_dist = sqrt(pow( (source_xyz[0] - intersections.first.x), 2.0 ) + pow( (source_xyz[1] - intersections.first.y), 2.0 ));
  BOOST_CHECK_SMALL( dist - test_dist, tollerance );
  
  dist = cylinder_line_intersection( radius, half_length, exit_point, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( exit_point[0] - intersections.second.x, tollerance );
  BOOST_CHECK_SMALL( exit_point[1] - intersections.second.y, tollerance );
  BOOST_CHECK_SMALL( exit_point[2] - 0, 1.0E-6 );
  test_dist = sqrt(pow( (intersections.first.x - intersections.second.x), 2.0 ) + pow( (intersections.first.y - intersections.second.y), 2.0 ));
  BOOST_CHECK_SMALL( dist - test_dist, tollerance );

  
  radius = 1;
  half_length = 1;
  source_xyz[0] = 1.0; source_xyz[1] = 0; source_xyz[2] = -5;
  detector_xyz[0] = 1.0; detector_xyz[1] = 0; detector_xyz[2] = 5;
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::AwayFromDetector, exit_point );
  BOOST_CHECK_SMALL( exit_point[0] - 1.0, 1.0E-12 );
  BOOST_CHECK_SMALL( exit_point[1] - 0, 1.0E-12 );
  BOOST_CHECK_SMALL( exit_point[2] - -1.0, 1.0E-12 );
  BOOST_CHECK_EQUAL( dist, 4.0 );
  
  
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( exit_point[0] - 1.0, 1.0E-12 );
  BOOST_CHECK_SMALL( exit_point[1] - 0, 1.0E-12 );
  BOOST_CHECK_SMALL( exit_point[2] - 1.0, 1.0E-12 );
  BOOST_CHECK_EQUAL( dist, 6.0 );
  
   
  // TODO: need more tests for exactly on boundary, or whatever
  
  // Check case where source is between cylinder and detector, so line segment doesnt actually
  //  intersect the cylinder, although it would if line was infinitely long.
  radius = 6.3499999999999996;
  half_length = 106.67999999999999;
  source_xyz[0] = 12.691863768989704; source_xyz[1] = 0.050459131221286244; source_xyz[2] = 0.0;
  detector_xyz[0] = 269.23999999999995; detector_xyz[1] = 0; detector_xyz[2] = 5;
  double to_enter_dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::AwayFromDetector, exit_point );
  double to_exit_dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_EQUAL( to_enter_dist, 0.0 );
  BOOST_CHECK_EQUAL( to_exit_dist, 0.0 );
  
  
  
  //radius = 6.3499999999999996;
  //half_length = 106.67999999999999;
  //// source is at radius 6.35002389242
  //source[0] = -2.430014265260013; source[1] = 5.8666714672787954; source[2] = 106.68019720704258;
  //detector[0] = 269.23999999999995; detector[1] = 0; detector[2] = 0;
  //dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::AwayFromDetector, exit_point );
  //dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );

  
  radius = 6.3499999999999996;
  half_length = 106.67999999999999;
  source_xyz[0] = 4.4901280605345768; source_xyz[1] = 4.4901280605345759; source_xyz[2] = 0;
  detector_xyz[0] = 269.23999999999995; detector_xyz[1] = 0; detector_xyz[2] = 0;
  
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_EQUAL( dist, 0.0 );
  
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::AwayFromDetector, exit_point );
  BOOST_CHECK_EQUAL( dist, 0.0 );
  
  
  radius = 6.3499999999999996;
  half_length = 106.67999999999999;
  source_xyz[0] = -2.430014265260013; source_xyz[1] = 5.8666714672787954; source_xyz[2] = 106.68019720704258;
  detector_xyz[0] = 269.23999999999995; detector_xyz[1] = 0 ;detector_xyz[2] = 0;
  
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::TowardDetector, exit_point );
  BOOST_CHECK_SMALL( radius - sqrt(exit_point[0]*exit_point[0] + exit_point[1]*exit_point[1]), 1.0E-9*(std::max)(1.0,radius) );
  
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::AwayFromDetector, exit_point );
  BOOST_CHECK_SMALL( exit_point[2] - half_length, 1.0E-9*(std::max)(1.0,half_length) ); //exit on end
}

BOOST_AUTO_TEST_CASE( CylinderEndOnExitCap )
{
  // For end-on geometry the detector sits on the cylinder axis, so the detector's xy-projection is
  //  the circle's *center*, and both crossings of the infinite cylinder are exactly equidistant from
  //  it.  Ordering the crossings by that distance (which #cylinder_line_intersection used to do) is
  //  therefore a rounding-noise coin flip, and when it flips the ray is sent out the wrong end cap.
  //  Every case below has an analytic answer, so this pins the ordering down.

  const double radius = 2.0;
  const double half_length = 0.15;
  const double det_dist = 100.0;
  const double detector[3] = { 0.0, 0.0, det_dist };

  const double radial_fracs[] = { 0.0, 0.01, 0.3, 0.7, 0.999 };
  const double z_fracs[] = { -0.99, -0.5, 0.0, 0.5, 0.99 };
  const size_t num_theta = 72;

  for( const double r_frac : radial_fracs )
  {
    const double r = r_frac * radius;

    for( const double z_frac : z_fracs )
    {
      const double z = z_frac * half_length;

      for( size_t theta_index = 0; theta_index < num_theta; ++theta_index )
      {
        const double theta = (2.0 * PhysicalUnits::pi * theta_index) / num_theta;
        const double source[3] = { r*cos(theta), r*sin(theta), z };

        // Heading toward an on-axis detector the ray's radius only shrinks, so it always leaves
        //  through the +half_length cap.
        const double t_toward = (half_length - z) / (det_dist - z);
        const double expected_toward[3] = {
          source[0]*(1.0 - t_toward), source[1]*(1.0 - t_toward), half_length
        };
        const double expected_toward_dist = sqrt( (r*t_toward)*(r*t_toward)
                                                 + (half_length - z)*(half_length - z) );

        double exit_toward[3];
        const double dist_toward = cylinder_line_intersection( radius, half_length, source, detector,
                                                        CylExitDir::TowardDetector, exit_toward );

        BOOST_CHECK_SMALL( exit_toward[0] - expected_toward[0], 1.0E-12*radius );
        BOOST_CHECK_SMALL( exit_toward[1] - expected_toward[1], 1.0E-12*radius );
        BOOST_CHECK_SMALL( exit_toward[2] - expected_toward[2], 1.0E-12*half_length );
        BOOST_CHECK_SMALL( dist_toward - expected_toward_dist, 1.0E-12*(std::max)(1.0,expected_toward_dist) );

        // Going away from the detector the ray leaves either through the -half_length cap, or - if
        //  it reaches the radius first - through the side wall.  The crossings of the infinite
        //  cylinder are at t = 1 -+ radius/r (for an on-axis detector), so the backward one is
        //  1 - radius/r; whichever of that and the -half_length cap comes first (largest t, since
        //  both are negative) is where we leave.
        const double t_cap_away = (-half_length - z) / (det_dist - z);
        const double t_side_away = (r > 0.0) ? (1.0 - radius/r)
                                             : -std::numeric_limits<double>::infinity();
        const double t_away = (std::max)( t_cap_away, t_side_away );

        const double expected_away[3] = {
          source[0]*(1.0 - t_away), source[1]*(1.0 - t_away), z + t_away*(det_dist - z)
        };
        const double expected_away_dist = fabs(t_away) * sqrt( r*r + (det_dist - z)*(det_dist - z) );

        double exit_away[3];
        const double dist_away = cylinder_line_intersection( radius, half_length, source, detector,
                                                        CylExitDir::AwayFromDetector, exit_away );

        BOOST_CHECK_SMALL( exit_away[0] - expected_away[0], 1.0E-12*radius );
        BOOST_CHECK_SMALL( exit_away[1] - expected_away[1], 1.0E-12*radius );
        BOOST_CHECK_SMALL( exit_away[2] - expected_away[2], 1.0E-12*half_length );
        BOOST_CHECK_SMALL( dist_away - expected_away_dist, 1.0E-12*(std::max)(1.0,expected_away_dist) );

        // Both exit points must be on the surface, and on opposite sides of the source, so the two
        //  distances must add up to the chord between them.
        const double chord = GammaInteractionCalc::distance( exit_toward, exit_away );
        BOOST_CHECK_SMALL( (dist_toward + dist_away) - chord, 1.0E-11*(std::max)(1.0,chord) );
      }//for( loop over theta )
    }//for( loop over z )
  }//for( loop over radius )
}//BOOST_AUTO_TEST_CASE( CylinderEndOnExitCap )


BOOST_AUTO_TEST_CASE( CylinderIntersectionRotationalInvariance )
{
  // With the detector on the cylinder axis the whole problem is symmetric about z, so rotating the
  //  source about the axis must not change the exit radius, exit z, or path length at all.  This is
  //  the cleanest statement of the coin-flip bug: it used to hold for only ~2/3 of the angles.

  const double radius = 6.35;
  const double half_length = 106.68;
  const double det_dist = 269.24;
  const double detector[3] = { 0.0, 0.0, det_dist };

  const double sources_r_z[][2] = {
    { 0.1*radius, 0.0 }, { 0.5*radius, 0.9*half_length }, { 0.99*radius, -0.5*half_length },
    { 0.25*radius, 0.999*half_length }, { 0.75*radius, -0.999*half_length }
  };

  for( const auto &r_and_z : sources_r_z )
  {
    const double r = r_and_z[0], z = r_and_z[1];

    for( const CylExitDir dir : { CylExitDir::TowardDetector, CylExitDir::AwayFromDetector } )
    {
      double ref_rad = -1.0, ref_z = 0.0, ref_dist = 0.0;

      for( size_t theta_index = 0; theta_index < 128; ++theta_index )
      {
        const double theta = (2.0 * PhysicalUnits::pi * theta_index) / 128;
        const double source[3] = { r*cos(theta), r*sin(theta), z };

        double exit_point[3];
        const double dist = cylinder_line_intersection( radius, half_length, source, detector,
                                                       dir, exit_point );
        const double exit_rad = sqrt( exit_point[0]*exit_point[0] + exit_point[1]*exit_point[1] );

        if( theta_index == 0 )
        {
          ref_rad = exit_rad;
          ref_z = exit_point[2];
          ref_dist = dist;
        }else
        {
          BOOST_CHECK_SMALL( exit_rad - ref_rad, 1.0E-12*(std::max)(1.0,ref_rad) );
          BOOST_CHECK_SMALL( exit_point[2] - ref_z, 1.0E-12*(std::max)(1.0,fabs(ref_z)) );
          BOOST_CHECK_SMALL( dist - ref_dist, 1.0E-12*(std::max)(1.0,ref_dist) );
        }

        // The azimuth of the exit point must also be preserved (the ray stays in the plane
        //  containing the axis and the source)
        BOOST_CHECK_SMALL( exit_point[0]*sin(theta) - exit_point[1]*cos(theta),
                          1.0E-11*(std::max)(1.0,ref_rad) );
      }//for( loop over theta )
    }//for( loop over direction )
  }//for( loop over source r,z )
}//BOOST_AUTO_TEST_CASE( CylinderIntersectionRotationalInvariance )


BOOST_AUTO_TEST_CASE( CylinderIntersectionReferenceSweep )
{
  // Compare against an independent long-double t-interval tracer, for sources strictly inside the
  //  cylinder (so both crossings always exist, and neither of the "line misses the volume"
  //  early-outs applies).  Covers end-on, side-on, and off-axis detector placements.

  std::mt19937 rng( 20260726 );
  std::uniform_real_distribution<double> uniform( 0.0, 1.0 );

  const double cyl_dims[][2] = { {2.0, 0.15}, {6.35, 106.68}, {225000.0, 1000.0}, {0.5, 0.5} };

  for( const auto &dims : cyl_dims )
  {
    const double radius = dims[0], half_length = dims[1];
    const double det_dist = 20.0 * (std::max)( radius, half_length );

    for( size_t iter = 0; iter < 2000; ++iter )
    {
      // Sources strictly inside; the 0.999 keeps us off the exact boundary, where the "on the
      //  surface counts as outside" convention is still not nailed down (see the TODO in
      //  cylinder_line_intersection).
      const double r = 0.999 * radius * uniform(rng);
      const double theta = 2.0 * PhysicalUnits::pi * uniform(rng);
      const double source[3] = { r*cos(theta), r*sin(theta), 0.999*half_length*(2.0*uniform(rng) - 1.0) };

      double detector[3];
      switch( iter % 3 )
      {
        case 0:  //end-on: detector on the axis
          detector[0] = 0.0; detector[1] = 0.0; detector[2] = det_dist;
          break;

        case 1:  //side-on: detector on the x-axis
          detector[0] = det_dist; detector[1] = 0.0; detector[2] = 0.0;
          break;

        default: //arbitrary direction
        {
          const double det_theta = 2.0 * PhysicalUnits::pi * uniform(rng);
          const double det_phi = PhysicalUnits::pi * uniform(rng);
          detector[0] = det_dist * sin(det_phi) * cos(det_theta);
          detector[1] = det_dist * sin(det_phi) * sin(det_theta);
          detector[2] = det_dist * cos(det_phi);
          break;
        }
      }//switch( iter % 3 )

      for( const bool toward : { true, false } )
      {
        const CylExitDir dir = toward ? CylExitDir::TowardDetector : CylExitDir::AwayFromDetector;

        double ref_exit[3], ref_dist;
        BOOST_REQUIRE( reference_cyl_exit( radius, half_length, source, detector, toward,
                                          ref_exit, ref_dist ) );

        double exit_point[3];
        const double dist = cylinder_line_intersection( radius, half_length, source, detector,
                                                       dir, exit_point );

        const double scale = (std::max)( radius, half_length );
        BOOST_CHECK_SMALL( exit_point[0] - ref_exit[0], 1.0E-9*scale );
        BOOST_CHECK_SMALL( exit_point[1] - ref_exit[1], 1.0E-9*scale );
        BOOST_CHECK_SMALL( exit_point[2] - ref_exit[2], 1.0E-9*scale );
        BOOST_CHECK_SMALL( dist - ref_dist, 1.0E-9*(std::max)(scale,ref_dist) );
      }//for( loop over direction )
    }//for( loop over random configurations )
  }//for( loop over cylinder dimensions )
}//BOOST_AUTO_TEST_CASE( CylinderIntersectionReferenceSweep )


BOOST_AUTO_TEST_CASE( CylinderThinOuterShellPathLength )
{
  // The symptom that made this matter: for a multi-layer end-on cylinder, `eval_cylinder` walks the
  //  ray out of the source cylinder and then charges each outer shell whatever
  //  #cylinder_line_intersection returns for it.  When the toward/away pick flipped, a shell only a
  //  micron thicker than the source got charged the *full stack length*.

  const double radius = 2.0 * PhysicalUnits::cm;
  const double half_length = 0.15 * PhysicalUnits::cm;
  const double shell_thickness = 1.0E-4 * PhysicalUnits::cm;  //1 micron
  const double det_dist = 100.0 * PhysicalUnits::cm;
  const double detector[3] = { 0.0, 0.0, det_dist };

  for( size_t theta_index = 0; theta_index < 64; ++theta_index )
  {
    const double theta = (2.0 * PhysicalUnits::pi * theta_index) / 64;

    for( const double r_frac : { 0.0, 0.05, 0.5, 0.95 } )
    {
      for( const double z_frac : { -0.9, 0.0, 0.9 } )
      {
        const double r = r_frac * radius;
        const double source[3] = { r*cos(theta), r*sin(theta), z_frac*half_length };

        double src_exit[3];
        const double dist_in_src = cylinder_line_intersection( radius, half_length, source, detector,
                                                             CylExitDir::TowardDetector, src_exit );
        BOOST_CHECK_GT( dist_in_src, 0.0 );
        BOOST_CHECK_SMALL( src_exit[2] - half_length, 1.0E-12*half_length );

        double shell_exit[3];
        const double dist_in_shell = cylinder_line_intersection( radius + shell_thickness,
                                                    half_length + shell_thickness, src_exit,
                                                    detector, CylExitDir::TowardDetector, shell_exit );

        // The ray is near-normal to the cap (radius/det_dist ~ 0.02), so the path through the shell
        //  must be the shell thickness to well within a percent - not the ~2*half_length it became
        //  when the ray was sent out the wrong end cap.
        BOOST_CHECK_SMALL( dist_in_shell - shell_thickness, 0.01*shell_thickness );
        BOOST_CHECK_SMALL( shell_exit[2] - (half_length + shell_thickness),
                          1.0E-9*(half_length + shell_thickness) );
      }//for( loop over z )
    }//for( loop over radius )
  }//for( loop over theta )
}//BOOST_AUTO_TEST_CASE( CylinderThinOuterShellPathLength )


BOOST_AUTO_TEST_CASE( RectangularIntersections )
{
  bool intersected;
  double half_width, half_height, half_depth, dist_in_shape;
  double source[3], detector[3], exit_point[3], enter_point[3];
  
  // First we'll test the simple case where we know the ray exits the volume on the plane at
  //  +-half_depth on the z-axis
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 0.0; detector[1] = 0.0; detector[2] = 10.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  BOOST_CHECK_SMALL( fabs(dist_in_shape - 1.0), 1.0E-9*(std::max)(1.0,dist_in_shape) );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 0.0), 1.0E-9*(std::max)(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 0.0), 1.0E-9*(std::max)(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - half_depth), 1.0E-9*(std::max)(half_depth,fabs(exit_point[2])) );
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 1.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 0.0; detector[1] = 0.0; detector[2] = 2.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  BOOST_CHECK_SMALL( fabs(dist_in_shape - sqrt(1.0*1.0 + 0.5*0.5)), 1.0E-9*(std::max)(1.0,dist_in_shape) );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 0.5), 1.0E-9*(std::max)(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 0.0), 1.0E-9*(std::max)(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - half_depth), 1.0E-9*(std::max)(half_depth,fabs(exit_point[2])) );
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.5; source[1] = -0.5; source[2] = -1.0;
  detector[0] = 0.0; detector[1] = 0.0; detector[2] = 3.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  BOOST_CHECK_SMALL( fabs(dist_in_shape - sqrt(0.25*0.25 + 0.25*0.25 + 2.0*2.0)), 1.0E-9*(std::max)(1.0,dist_in_shape) );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 0.25), 1.0E-9*(std::max)(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - -0.25), 1.0E-9*(std::max)(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - half_depth), 1.0E-9*(std::max)(half_depth,fabs(exit_point[2])) );
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.5; source[1] = -0.5; source[2] = -1.0;
  detector[0] = 0.5; detector[1] = -0.5; detector[2] = 3.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  BOOST_CHECK_SMALL( fabs(dist_in_shape - 2), 1.0E-9*(std::max)(1.0,dist_in_shape) );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 0.5), 1.0E-9*(std::max)(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - -0.5), 1.0E-9*(std::max)(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - half_depth), 1.0E-9*(std::max)(half_depth,fabs(exit_point[2])) );
  
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 0.5; detector[1] = -0.5; detector[2] = -2.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  BOOST_CHECK_SMALL( fabs(dist_in_shape - sqrt(0.25*0.25 + 0.25*0.25 + 1.0*1.0)), 1.0E-9*(std::max)(1.0,dist_in_shape) );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 0.25), 1.0E-9*(std::max)(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - -0.25), 1.0E-9*(std::max)(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - -half_depth), 1.0E-9*(std::max)(half_depth,fabs(exit_point[2])) );
  
  
  // Now test the other cases of the ray exiting the rectangle on an arbitrary face.
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 2.0; detector[1] = 0.0; detector[2] = 0.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  BOOST_CHECK_SMALL( fabs(dist_in_shape - 1.0), 1.0E-9*(std::max)(1.0,dist_in_shape) );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 1.0), 1.0E-9*(std::max)(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 0.0), 1.0E-9*(std::max)(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - 0.0), 1.0E-9*(std::max)(half_depth,fabs(exit_point[2])) );
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = -2.0; detector[1] = 0.0; detector[2] = 0.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  BOOST_CHECK_SMALL( fabs(dist_in_shape - 1.0), 1.0E-9*(std::max)(1.0,dist_in_shape) );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - -1.0), 1.0E-9*(std::max)(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 0.0), 1.0E-9*(std::max)(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - 0.0), 1.0E-9*(std::max)(half_depth,fabs(exit_point[2])) );
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 0.0; detector[1] = 2.0; detector[2] = 0.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  BOOST_CHECK_SMALL( fabs(dist_in_shape - 1.0), 1.0E-9*(std::max)(1.0,dist_in_shape) );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 0.0), 1.0E-9*(std::max)(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 1.0), 1.0E-9*(std::max)(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - 0.0), 1.0E-9*(std::max)(half_depth,fabs(exit_point[2])) );
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 0.0; detector[1] = -2.0; detector[2] = 0.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  BOOST_CHECK_SMALL( fabs(dist_in_shape - 1.0), 1.0E-9*(std::max)(1.0,dist_in_shape) );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 0.0), 1.0E-9*(std::max)(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - -1.0), 1.0E-9*(std::max)(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - 0.0), 1.0E-9*(std::max)(half_depth,fabs(exit_point[2])) );
  
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = -1.0; source[1] = 0.5; source[2] = -0.5;
  detector[0] = 3.0; detector[1] = 0.0; detector[2] = 0.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  BOOST_CHECK_SMALL( fabs(dist_in_shape - sqrt(0.25*0.25 + 0.25*0.25 + 2.0*2.0)), 1.0E-9*(std::max)(1.0,dist_in_shape) );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - half_width), 1.0E-9*(std::max)(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 0.25), 1.0E-9*(std::max)(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - -0.25), 1.0E-9*(std::max)(1.0,fabs(exit_point[2])) );
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 1.0; source[1] = -1.0; source[2] = -1.0;
  detector[0] = 0.0; detector[1] = 3.0; detector[2] = 0.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  BOOST_CHECK_SMALL( fabs(dist_in_shape - sqrt(0.5*0.5 + 0.5*0.5 + 2.0*2.0)), 1.0E-9*(std::max)(1.0,dist_in_shape) );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 0.5), 1.0E-9*(std::max)(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 1.0), 1.0E-9*(std::max)(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - -0.5), 1.0E-9*(std::max)(1.0,fabs(exit_point[2])) );
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = -0.5; source[1] = -0.5; source[2] = 0.5;
  detector[0] = 2.5; detector[1] = 1.0; detector[2] = 1.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 1.0), 1.0E-9*(std::max)(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 0.25), 1.0E-9*(std::max)(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - 0.75), 1.0E-9*(std::max)(1.0,fabs(exit_point[2])) );
  
  
  // Test rectangle_intersections function
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = -2.0; source[1] = 2.0; source[2] = 0.0;
  detector[0] = 2.0; detector[1] = 2.0; detector[2] = 0.0;
  intersected = rectangle_intersections( half_width, half_height, half_depth,
                                        source, detector, enter_point, exit_point );
  BOOST_CHECK( !intersected );
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = -2.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 2.0; detector[1] = 0.0; detector[2] = 0.0;
  intersected = rectangle_intersections( half_width, half_height, half_depth,
                                        source, detector, enter_point, exit_point );
  BOOST_CHECK( intersected );
  BOOST_CHECK_EQUAL( enter_point[0], -1.0 );
  BOOST_CHECK_EQUAL( enter_point[1], 0.0 );
  BOOST_CHECK_EQUAL( enter_point[2], 0.0 );
  BOOST_CHECK_EQUAL( exit_point[0], 1.0 );
  BOOST_CHECK_EQUAL( exit_point[1], 0.0 );
  BOOST_CHECK_EQUAL( exit_point[2], 0.0 );
  
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.0; source[1] = 0.0; source[2] = -2.0;
  detector[0] = 0.0; detector[1] = 0.0; detector[2] = 2.0;
  intersected = rectangle_intersections( half_width, half_height, half_depth,
                                        source, detector, enter_point, exit_point );
  BOOST_CHECK( intersected );
  BOOST_CHECK_EQUAL( enter_point[0], 0.0 );
  BOOST_CHECK_EQUAL( enter_point[1], 0.0 );
  BOOST_CHECK_EQUAL( enter_point[2], -1.0 );
  BOOST_CHECK_EQUAL( exit_point[0], 0.0 );
  BOOST_CHECK_EQUAL( exit_point[1], 0.0 );
  BOOST_CHECK_EQUAL( exit_point[2], 1.0 );
  
  half_width = 10.0; half_height = 10.0; half_depth = 10.0;
  source[0] = -10; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 0.0; detector[1] = 0.0; detector[2] = 250;
  intersected = rectangle_intersections( half_width, half_height, half_depth,
                                        source, detector, enter_point, exit_point );
  BOOST_CHECK( intersected );
} 


BOOST_AUTO_TEST_CASE( DistributedSrcCalcTests )
{
  set_data_dir();
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );
  const string materialfile = SpecUtils::append_path( InterSpec::staticDataDirectory(), "MaterialDataBase.txt" );
  
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  BOOST_REQUIRE( MaterialDB::initialized() );

  const std::shared_ptr<const MaterialDB> materialdb = MaterialDB::instance();
  BOOST_REQUIRE( materialdb );

  DistributedSrcCalc ObjectToIntegrate;

  const double energy = 185.0*PhysicalUnits::keV;
  ObjectToIntegrate.m_geometry = GeometryType::Spherical;
  ObjectToIntegrate.m_materialIndex = 0;
  ObjectToIntegrate.m_attenuateForAir = false;
  ObjectToIntegrate.m_airTransLenCoef = 0.0;
  ObjectToIntegrate.m_isInSituExponential = false;
  ObjectToIntegrate.m_inSituRelaxationLength = 0.0;
  ObjectToIntegrate.m_detectorRadius  = 2.0 * PhysicalUnits::cm;
  ObjectToIntegrate.m_observationDist = 400.0 * PhysicalUnits::cm;

  double sphereRad = 0.0, transLenCoef = 0.0;

  std::shared_ptr<const Material> material = materialdb->material( "void" );
  transLenCoef = GammaInteractionCalc::transmition_length_coefficient( material.get(), energy );
  sphereRad += 99.5* PhysicalUnits::cm;
#if( defined(__GNUC__) && __GNUC__ < 5 )
  ObjectToIntegrate.m_dimensionsTransLenAndType.push_back( tuple<array<double,3>,double,DistributedSrcCalc::ShellType>{{sphereRad,0.0,0.0},transLenCoef,DistributedSrcCalc::ShellType::Material} );
#else
  ObjectToIntegrate.m_dimensionsTransLenAndType.push_back( {{sphereRad,0.0,0.0},transLenCoef,DistributedSrcCalc::ShellType::Material} );
#endif

  material = materialdb->material( "U" );
  transLenCoef = GammaInteractionCalc::transmition_length_coefficient( material.get(), energy );
  sphereRad += 0.5 * PhysicalUnits::cm;
#if( defined(__GNUC__) && __GNUC__ < 5 )
  ObjectToIntegrate.m_dimensionsTransLenAndType.push_back( tuple<array<double,3>,double,DistributedSrcCalc::ShellType>{{sphereRad,0.0,0.0},transLenCoef,DistributedSrcCalc::ShellType::Material} );
#else
  ObjectToIntegrate.m_dimensionsTransLenAndType.push_back( {{sphereRad,0.0,0.0},transLenCoef,DistributedSrcCalc::ShellType::Material} );
#endif
  ObjectToIntegrate.m_materialIndex = ObjectToIntegrate.m_dimensionsTransLenAndType.size() - 1;

  material = materialdb->material( "Fe" );
  transLenCoef = GammaInteractionCalc::transmition_length_coefficient( material.get(), energy );
  sphereRad += 0.5 * PhysicalUnits::cm;
  ObjectToIntegrate.m_dimensionsTransLenAndType.push_back( {{sphereRad,0.0,0.0},transLenCoef,DistributedSrcCalc::ShellType::Material} );

  
  int ndim = 2;  //the number of dimensions of the integral.
  void *userdata = (void *)&ObjectToIntegrate;
  const double epsrel = 1e-5;  //the requested relative accuracy
  const double epsabs = -1.0;//1e-12; //the requested absolute accuracy
  const int mineval = 0; //the minimum number of integrand evaluations required.
  const int maxeval = 5000000; //the (approximate) maximum number of integrand evaluations allowed.

  int nregions, neval, fail;
  double integral, error, prob;

  ndim = 2;
  Integrate::CuhreIntegrate( ndim, DistributedSrcCalc_integrand_spherical, userdata, epsrel, epsabs,
                            Integrate::LastImportanceFcnt, mineval, maxeval, nregions, neval,
                            fail, integral, error, prob );

  //printf("ndim=%d CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n",
  //    ndim, nregions, neval, fail);
  //printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n", integral, error, prob);
  //printf("\n\n" );
  
  // Check that the integration succeeded
  BOOST_CHECK_EQUAL( fail, 0 );
  BOOST_CHECK_GT( neval, 0 );
  BOOST_CHECK_GT( nregions, 0 );
  
  // Check that the result is close to the expected value from the comment
  BOOST_CHECK_CLOSE( integral, 2.8626, 0.1 );  // within 0.01%
  
  ndim = 3;
  Integrate::CuhreIntegrate( ndim, DistributedSrcCalc_integrand_spherical, userdata, epsrel, epsabs,
                             Integrate::LastImportanceFcnt, mineval, maxeval, nregions, neval,
                            fail, integral, error, prob );
  
  //printf("ndim=%d CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n",
  //     ndim, nregions, neval, fail);
  //printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n", integral, error, prob );
  //cout << endl << endl;
  
  // Check that the 3D integration also succeeded
  BOOST_CHECK_EQUAL( fail, 0 );
  BOOST_CHECK_GT( neval, 0 );
  BOOST_CHECK_GT( nregions, 0 );
  
  // Check that the result is close to the expected value from the comment
  BOOST_CHECK_CLOSE( integral, 2.8626, 0.01 );  // within 0.1%
}//BOOST_AUTO_TEST_CASE( DistributedSrcCalcTests )


/** Tests #MassAttenuation::mass_atten_coef_frac_an - the C1-continuous
 fractional-atomic-number attenuation coefficient used (via
 massAttenuationCoefficientFracAN and RelActCalc::get_atten_coef_for_an) by
 generic-shielding transmission, and by auto-differentiated fits of atomic
 number.
 */
BOOST_AUTO_TEST_CASE( FracAtomicNumberAttenuationCoef )
{
  set_data_dir();

  const double test_energies[] = { 59.5, 185.7, 661.7, 2614.5 };  //keV

  // 1) Must reproduce the element values exactly at integer atomic number
  for( const double energy : test_energies )
  {
    for( const int an : { 1, 13, 26, 53, 74, 92, 98 } )
    {
      const double frac_mu = MassAttenuation::mass_atten_coef_frac_an<double>( an, energy );
      const double elem_mu = MassAttenuation::massAttenuationCoefficientElement( an, energy );
      BOOST_CHECK_CLOSE( frac_mu, elem_mu, 1.0E-4 );
    }
  }//for( loop over energies )

  // 2) Value and derivative continuity across integer atomic numbers.
  //    (the previous linear interpolation had a discontinuous derivative here)
  const double delta = 1.0E-4;
  for( const double energy : test_energies )
  {
    for( const double an : { 26.0, 74.0 } )
    {
      const double mu_below = MassAttenuation::mass_atten_coef_frac_an<double>( an - delta, energy );
      const double mu_at = MassAttenuation::mass_atten_coef_frac_an<double>( an, energy );
      const double mu_above = MassAttenuation::mass_atten_coef_frac_an<double>( an + delta, energy );

      BOOST_CHECK_CLOSE( mu_below, mu_at, 0.1 );
      BOOST_CHECK_CLOSE( mu_above, mu_at, 0.1 );

      const double slope_below = (mu_at - mu_below) / delta;
      const double slope_above = (mu_above - mu_at) / delta;
      BOOST_CHECK_MESSAGE( fabs(slope_above - slope_below) <= 0.01*std::max(fabs(slope_above),fabs(slope_below)),
                           "Derivative discontinuity at AN=" << an << ", E=" << energy
                           << " keV: slope below=" << slope_below << ", above=" << slope_above );
    }
  }//for( loop over energies )

  // 3) ceres::Jet derivative must match a numeric derivative of the double path,
  //    including exactly at an integer atomic number (a knot of the interpolation)
  using Jet1 = ceres::Jet<double,1>;
  const double diff_step = 1.0E-5;
  for( const double energy : { 59.5, 661.7 } )
  {
    for( const double an : { 25.5, 26.0, 73.2 } )
    {
      const Jet1 an_jet( an, 0 );
      const Jet1 mu_jet = MassAttenuation::mass_atten_coef_frac_an( an_jet, static_cast<float>(energy) );

      const double mu_val = MassAttenuation::mass_atten_coef_frac_an<double>( an, energy );
      const double mu_plus = MassAttenuation::mass_atten_coef_frac_an<double>( an + diff_step, energy );
      const double mu_minus = MassAttenuation::mass_atten_coef_frac_an<double>( an - diff_step, energy );
      const double numeric_deriv = (mu_plus - mu_minus) / (2.0*diff_step);

      BOOST_CHECK_CLOSE( mu_jet.a, mu_val, 1.0E-9 );
      BOOST_CHECK_MESSAGE( fabs(mu_jet.v[0] - numeric_deriv) <= 1.0E-4*fabs(numeric_deriv),
                           "Jet derivative (" << mu_jet.v[0] << ") doesnt match numeric ("
                           << numeric_deriv << ") at AN=" << an << ", E=" << energy << " keV" );
    }
  }//for( loop over energies )
}//BOOST_AUTO_TEST_CASE( FracAtomicNumberAttenuationCoef )
namespace
{
  /** Sets up a `CylinderEndOn` distributed source: a U cylinder of the given dimensions, optionally
   wrapped in an Fe shell of thickness `shell_thickness` (pass <= 0 for no shell).

   The source material is always the inner-most shell, e.g. `m_materialIndex` is 0.
   */
  DistributedSrcCalc make_end_on_cylinder( const double energy, const double radius,
                                          const double half_length, const double shell_thickness,
                                          const double obs_dist )
  {
    const std::shared_ptr<const MaterialDB> materialdb = MaterialDB::instance();
    BOOST_REQUIRE( materialdb );

    DistributedSrcCalc calc;
    calc.m_geometry = GeometryType::CylinderEndOn;
    calc.m_materialIndex = 0;
    calc.m_attenuateForAir = false;
    calc.m_airTransLenCoef = 0.0;
    calc.m_isInSituExponential = false;
    calc.m_inSituRelaxationLength = 0.0;
    calc.m_detectorRadius = 2.0 * PhysicalUnits::cm;
    calc.m_detectorSetback = 0.0;
    calc.m_observationDist = obs_dist;
    calc.m_energy = energy;
    calc.m_nuclide = nullptr;
    calc.m_srcVolumetricActivity = 1.0;
    calc.integral = 0.0;

    const std::shared_ptr<const Material> uranium = materialdb->material( "U" );
    BOOST_REQUIRE( uranium );
    const double u_coef = GammaInteractionCalc::transmition_length_coefficient( uranium.get(), energy );
    calc.m_dimensionsTransLenAndType.push_back( { {radius, half_length, 0.0}, u_coef,
                                                  DistributedSrcCalc::ShellType::Material } );

    if( shell_thickness > 0.0 )
    {
      const std::shared_ptr<const Material> iron = materialdb->material( "Fe" );
      BOOST_REQUIRE( iron );
      const double fe_coef = GammaInteractionCalc::transmition_length_coefficient( iron.get(), energy );
      calc.m_dimensionsTransLenAndType.push_back( { {radius + shell_thickness,
                                                     half_length + shell_thickness, 0.0}, fe_coef,
                                                    DistributedSrcCalc::ShellType::Material } );
    }//if( shell_thickness > 0.0 )

    return calc;
  }//make_end_on_cylinder(...)


  /** Cuhre-integrates a #DistributedSrcCalc using the general cylindrical integrand. */
  double integrate_cylindrical( const DistributedSrcCalc &calc, const int ndim )
  {
    void * const userdata = (void *)&calc;
    const double epsrel = 1e-6, epsabs = -1.0;
    const int mineval = 0, maxeval = 5000000;

    int nregions, neval, fail;
    double integral, error, prob;

    Integrate::CuhreIntegrate( ndim, DistributedSrcCalc_integrand_cylindrical, userdata, epsrel,
                              epsabs, Integrate::LastImportanceFcnt, mineval, maxeval, nregions,
                              neval, fail, integral, error, prob );

    BOOST_REQUIRE_EQUAL( fail, 0 );

    return integral;
  }//integrate_cylindrical(...)
}//namespace


BOOST_AUTO_TEST_CASE( CylinderEndOnSingleVsGeneralIntegrand )
{
  // For a single shielding, end-on cylinders are integrated with the closed-form
  //  `eval_single_cyl_end_on`, which never calls #cylinder_line_intersection.  The general
  //  `eval_cylinder` must give the same answer point-for-point; it is the discrepancy between the
  //  two that the (previously commented out) assert in `eval_single_cyl_end_on` checks for.
  set_data_dir();

  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  BOOST_REQUIRE( MaterialDB::initialized() );

  const double radius = 2.0 * PhysicalUnits::cm;
  const double half_length = 0.15 * PhysicalUnits::cm;
  const double obs_dist = 100.0 * PhysicalUnits::cm;

  for( const double energy : { 60.0*PhysicalUnits::keV, 186.0*PhysicalUnits::keV, 2614.0*PhysicalUnits::keV } )
  {
    for( const bool attenuate_air : { false, true } )
    {
      DistributedSrcCalc calc = make_end_on_cylinder( energy, radius, half_length, -1.0, obs_dist );
      calc.m_attenuateForAir = attenuate_air;
      calc.m_airTransLenCoef = attenuate_air ? 1.0E-4 : 0.0;

      const int two_dim = 2, three_dim = 3;

      for( size_t i = 0; i <= 8; ++i )
      {
        for( size_t j = 0; j <= 8; ++j )
        {
          double xx_2d[2] = { i/8.0, j/8.0 };  //{r, z}

          double single_ff[1] = { 0.0 }, general_ff[1] = { 0.0 };
          calc.eval_single_cyl_end_on( xx_2d, &two_dim, single_ff, nullptr );
          calc.eval_cylinder( xx_2d, &two_dim, general_ff, nullptr );

          const double scale = (std::max)( fabs(single_ff[0]), fabs(general_ff[0]) );
          BOOST_CHECK_SMALL( single_ff[0] - general_ff[0], 1.0E-9*(std::max)(1.0E-12,scale) );

          // End-on is symmetric about the cylinder axis, so the 3-dimensional integrand must give
          //  the same value at every theta as the theta-collapsed 2-dimensional one.
          for( size_t k = 0; k <= 8; ++k )
          {
            double xx_3d[3] = { i/8.0, k/8.0, j/8.0 };  //{r, theta, z}

            double three_dim_ff[1] = { 0.0 };
            calc.eval_cylinder( xx_3d, &three_dim, three_dim_ff, nullptr );

            BOOST_CHECK_SMALL( single_ff[0] - three_dim_ff[0], 1.0E-9*(std::max)(1.0E-12,scale) );
          }//for( loop over theta )
        }//for( loop over z )
      }//for( loop over r )
    }//for( loop over whether to attenuate for air )
  }//for( loop over energy )
}//BOOST_AUTO_TEST_CASE( CylinderEndOnSingleVsGeneralIntegrand )


BOOST_AUTO_TEST_CASE( CylinderEndOnMultiShellTransport )
{
  // A multi-layer end-on cylinder routes to `eval_cylinder`, which is where the toward/away pick
  //  mattered.  A thin outer shell must attenuate by (essentially) exp(-mu*thickness); previously
  //  rays were charged the full stack length instead, biasing these integrals low by 3-36%.
  set_data_dir();

  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  BOOST_REQUIRE( MaterialDB::initialized() );

  const std::shared_ptr<const MaterialDB> materialdb = MaterialDB::instance();
  BOOST_REQUIRE( materialdb );
  const std::shared_ptr<const Material> iron = materialdb->material( "Fe" );
  BOOST_REQUIRE( iron );

  const double radius = 2.0 * PhysicalUnits::cm;
  const double half_length = 0.15 * PhysicalUnits::cm;
  const double obs_dist = 100.0 * PhysicalUnits::cm;

  // The rays are near-normal to the end cap (radius/obs_dist ~ 0.02), so the path through the shell
  //  is its thickness to within ~0.02%; allow a bit of slop on top of the integration accuracy.
  const double max_obliquity = sqrt( 1.0 + (radius/(obs_dist - half_length))*(radius/(obs_dist - half_length)) );

  for( const double energy : { 60.0*PhysicalUnits::keV, 2614.0*PhysicalUnits::keV } )
  {
    const double fe_coef = GammaInteractionCalc::transmition_length_coefficient( iron.get(), energy );

    const DistributedSrcCalc bare = make_end_on_cylinder( energy, radius, half_length, -1.0, obs_dist );
    const double bare_integral = integrate_cylindrical( bare, 2 );
    BOOST_REQUIRE_GT( bare_integral, 0.0 );

    double prev_ratio = 1.0;

    for( const double thickness : { 1.0E-4, 1.0E-2, 0.1, 0.5 } )  //cm; first is 1 micron
    {
      const double shell_thickness = thickness * PhysicalUnits::cm;
      const DistributedSrcCalc shielded = make_end_on_cylinder( energy, radius, half_length,
                                                              shell_thickness, obs_dist );
      const double shielded_integral = integrate_cylindrical( shielded, 2 );
      const double ratio = shielded_integral / bare_integral;

      const double expected_hi = exp( -fe_coef * shell_thickness );
      const double expected_lo = exp( -fe_coef * shell_thickness * max_obliquity );

      BOOST_CHECK_MESSAGE( (ratio <= (1.0 + 1.0E-6)) && (ratio >= (expected_lo - 1.0E-3)),
        "Fe shell of " << thickness << " cm at " << energy/PhysicalUnits::keV << " keV attenuated by "
        << ratio << ", but expected within [" << expected_lo << ", " << expected_hi << "]" );

      // The shell is the only difference between the two integrals, so their ratio is just the
      //  (near-normal-incidence) slab attenuation through it.
      BOOST_CHECK_CLOSE( ratio, expected_hi, 0.1 );

      BOOST_CHECK_LE( ratio, prev_ratio + 1.0E-9 );  //more shielding can only attenuate more
      prev_ratio = ratio;
    }//for( loop over shell thickness )
  }//for( loop over energy )
}//BOOST_AUTO_TEST_CASE( CylinderEndOnMultiShellTransport )
