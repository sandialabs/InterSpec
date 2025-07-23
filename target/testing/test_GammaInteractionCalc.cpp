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

#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>
#include <cstdio>
#include <stdexcept>
#include <utility>

//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE GammaInteractionCalc_suite
//#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/Integrate.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/GammaInteractionCalc.h"

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
    for( const auto &d : { "test_data", "../test_data", "../../test_data", "/Users/wcjohns/rad_ana/InterSpec/target/testing/test_data" } )
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
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( required_data_path ), "'" << required_data_file << "' not at '" << required_data_path << "'" );
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
  detector_xyz[0] = (1.0 - 2.0E-6); detector_xyz[1] = 1; detector_xyz[2] = 0;
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
  BOOST_CHECK_SMALL( radius - sqrt(exit_point[0]*exit_point[0] + exit_point[1]*exit_point[1]), 1.0E-9*std::max(1.0,radius) );
  
  dist = cylinder_line_intersection( radius, half_length, source_xyz, detector_xyz, CylExitDir::AwayFromDetector, exit_point );
  BOOST_CHECK_SMALL( exit_point[2] - half_length, 1.0E-9*std::max(1.0,half_length) ); //exit on end
}

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
  BOOST_CHECK_SMALL( fabs(dist_in_shape - 1.0), 1.0E-9*std::max(1.0,dist_in_shape) );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 0.0), 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 0.0), 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - half_depth), 1.0E-9*std::max(half_depth,fabs(exit_point[2])) );
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 1.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 0.0; detector[1] = 0.0; detector[2] = 2.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  BOOST_CHECK_SMALL( fabs(dist_in_shape - sqrt(1.0*1.0 + 0.5*0.5)), 1.0E-9*std::max(1.0,dist_in_shape) );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 0.5), 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 0.0), 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - half_depth), 1.0E-9*std::max(half_depth,fabs(exit_point[2])) );
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.5; source[1] = -0.5; source[2] = -1.0;
  detector[0] = 0.0; detector[1] = 0.0; detector[2] = 3.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  BOOST_CHECK_SMALL( fabs(dist_in_shape - sqrt(0.25*0.25 + 0.25*0.25 + 2.0*2.0)), 1.0E-9*std::max(1.0,dist_in_shape) );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 0.25), 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - -0.25), 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - half_depth), 1.0E-9*std::max(half_depth,fabs(exit_point[2])) );
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.5; source[1] = -0.5; source[2] = -1.0;
  detector[0] = 0.5; detector[1] = -0.5; detector[2] = 3.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  BOOST_CHECK_SMALL( fabs(dist_in_shape - 2), 1.0E-9*std::max(1.0,dist_in_shape) );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 0.5), 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - -0.5), 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - half_depth), 1.0E-9*std::max(half_depth,fabs(exit_point[2])) );
  
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 0.5; detector[1] = -0.5; detector[2] = -2.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  BOOST_CHECK_SMALL( fabs(dist_in_shape - sqrt(0.25*0.25 + 0.25*0.25 + 1.0*1.0)), 1.0E-9*std::max(1.0,dist_in_shape) );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 0.25), 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - -0.25), 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - -half_depth), 1.0E-9*std::max(half_depth,fabs(exit_point[2])) );
  
  
  // Now test the other cases of the ray exiting the rectangle on an arbitrary face.
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 2.0; detector[1] = 0.0; detector[2] = 0.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  BOOST_CHECK_SMALL( fabs(dist_in_shape - 1.0), 1.0E-9*std::max(1.0,dist_in_shape) );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 1.0), 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 0.0), 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - 0.0), 1.0E-9*std::max(half_depth,fabs(exit_point[2])) );
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = -2.0; detector[1] = 0.0; detector[2] = 0.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  BOOST_CHECK_SMALL( fabs(dist_in_shape - 1.0), 1.0E-9*std::max(1.0,dist_in_shape) );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - -1.0), 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 0.0), 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - 0.0), 1.0E-9*std::max(half_depth,fabs(exit_point[2])) );
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 0.0; detector[1] = 2.0; detector[2] = 0.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  BOOST_CHECK_SMALL( fabs(dist_in_shape - 1.0), 1.0E-9*std::max(1.0,dist_in_shape) );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 0.0), 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 1.0), 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - 0.0), 1.0E-9*std::max(half_depth,fabs(exit_point[2])) );
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 0.0; detector[1] = -2.0; detector[2] = 0.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  BOOST_CHECK_SMALL( fabs(dist_in_shape - 1.0), 1.0E-9*std::max(1.0,dist_in_shape) );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 0.0), 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - -1.0), 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - 0.0), 1.0E-9*std::max(half_depth,fabs(exit_point[2])) );
  
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = -1.0; source[1] = 0.5; source[2] = -0.5;
  detector[0] = 3.0; detector[1] = 0.0; detector[2] = 0.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  BOOST_CHECK_SMALL( fabs(dist_in_shape - sqrt(0.25*0.25 + 0.25*0.25 + 2.0*2.0)), 1.0E-9*std::max(1.0,dist_in_shape) );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - half_width), 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 0.25), 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - -0.25), 1.0E-9*std::max(1.0,fabs(exit_point[2])) );
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 1.0; source[1] = -1.0; source[2] = -1.0;
  detector[0] = 0.0; detector[1] = 3.0; detector[2] = 0.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  BOOST_CHECK_SMALL( fabs(dist_in_shape - sqrt(0.5*0.5 + 0.5*0.5 + 2.0*2.0)), 1.0E-9*std::max(1.0,dist_in_shape) );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 0.5), 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 1.0), 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - -0.5), 1.0E-9*std::max(1.0,fabs(exit_point[2])) );
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = -0.5; source[1] = -0.5; source[2] = 0.5;
  detector[0] = 2.5; detector[1] = 1.0; detector[2] = 1.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  BOOST_CHECK_SMALL( fabs(exit_point[0] - 1.0), 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  BOOST_CHECK_SMALL( fabs(exit_point[1] - 0.25), 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  BOOST_CHECK_SMALL( fabs(exit_point[2] - 0.75), 1.0E-9*std::max(1.0,fabs(exit_point[2])) );
  
  
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
  
  MaterialDB materialdb;
  BOOST_REQUIRE_NO_THROW( materialdb.parseGadrasMaterialFile( materialfile, db, false ) );
  
  //const Material *material = materialdb.material( "Air" );
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

  const Material *material = materialdb.material( "void" );
  transLenCoef = GammaInteractionCalc::transmition_length_coefficient( material, energy );
  sphereRad += 99.5* PhysicalUnits::cm;
#if( defined(__GNUC__) && __GNUC__ < 5 )
  ObjectToIntegrate.m_dimensionsTransLenAndType.push_back( tuple<array<double,3>,double,DistributedSrcCalc::ShellType>{{sphereRad,0.0,0.0},transLenCoef,DistributedSrcCalc::ShellType::Material} );
#else
  ObjectToIntegrate.m_dimensionsTransLenAndType.push_back( {{sphereRad,0.0,0.0},transLenCoef,DistributedSrcCalc::ShellType::Material} );
#endif

  material = materialdb.material( "U" );
  transLenCoef = GammaInteractionCalc::transmition_length_coefficient( material, energy );
  sphereRad += 0.5 * PhysicalUnits::cm;
#if( defined(__GNUC__) && __GNUC__ < 5 )
  ObjectToIntegrate.m_dimensionsTransLenAndType.push_back( tuple<array<double,3>,double,DistributedSrcCalc::ShellType>{{sphereRad,0.0,0.0},transLenCoef,DistributedSrcCalc::ShellType::Material} );
#else
  ObjectToIntegrate.m_dimensionsTransLenAndType.push_back( {{sphereRad,0.0,0.0},transLenCoef,DistributedSrcCalc::ShellType::Material} );
#endif
  ObjectToIntegrate.m_materialIndex = ObjectToIntegrate.m_dimensionsTransLenAndType.size() - 1;

  material = materialdb.material( "Fe" );
  transLenCoef = GammaInteractionCalc::transmition_length_coefficient( material, energy );
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
