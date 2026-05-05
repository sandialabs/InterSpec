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
#include <vector>
#include <memory>
#include <iostream>


#define BOOST_TEST_MODULE GadrasShieldScatter_suite
#include <boost/test/included/unit_test.hpp>


#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/GadrasShieldScatter.h"


using namespace std;
using namespace boost::unit_test;


namespace
{
  string find_data_dir( const string &probe )
  {
    int argc = boost::unit_test::framework::master_test_suite().argc;
    char **argv = boost::unit_test::framework::master_test_suite().argv;

    string datadir;
    for( int i = 1; i < argc; ++i )
    {
      const string arg = argv[i];
      if( SpecUtils::istarts_with( arg, "--datadir=" ) )
        datadir = arg.substr( 10 );
    }
    SpecUtils::ireplace_all( datadir, "%20", " " );

    if( datadir.empty() )
    {
      const string trials[] = {
        "data", "../data", "../../data", "../../../data"
      };
      for( const string &d : trials )
      {
        if( SpecUtils::is_file( SpecUtils::append_path( d, probe ) ) )
        {
          datadir = d;
          break;
        }
      }
    }
    return datadir;
  }


}//anonymous namespace


static unique_ptr<GadrasShieldScatter> g_new;
static bool g_loaded = false;


static void ensure_loaded()
{
  if( g_loaded )
    return;
  g_loaded = true;

  const string data_dir = find_data_dir( "sandia.shieldscatter.db" );
  const string new_file = SpecUtils::append_path( data_dir, "sandia.shieldscatter.db" );

  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( new_file ), "Missing " << new_file );

  // GadrasShieldScatter needs MassAttenuationTool, which needs the static data dir.
  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( data_dir ) );

  BOOST_REQUIRE_NO_THROW( g_new.reset( new GadrasShieldScatter( new_file ) ) );
}


BOOST_AUTO_TEST_CASE( LoadFiles )
{
  ensure_loaded();
  BOOST_REQUIRE( g_new );
  BOOST_CHECK( g_new->groupCount() > 0 );
}


// Sanity: shapeFactor extremes (point=0, slab=1) should both produce
// non-negative finite results, and they should differ from each other for
// nontrivial shielding (otherwise the geometry knob does nothing).
BOOST_AUTO_TEST_CASE( ShapeFactorEndpoints )
{
  ensure_loaded();

  BOOST_REQUIRE( g_new );

  vector<double> point_scatter, slab_scatter;
  g_new->computeShieldScatter( 661.7, 26.0, 5.0, 0.0, 0.0, point_scatter );
  g_new->computeShieldScatter( 661.7, 26.0, 5.0, 0.0, 1.0, slab_scatter );

  BOOST_REQUIRE_EQUAL( point_scatter.size(), slab_scatter.size() );

  double point_sum = 0.0, slab_sum = 0.0, abs_diff = 0.0;
  for( size_t i = 0; i < point_scatter.size(); ++i )
  {
    BOOST_CHECK( std::isfinite( point_scatter[i] ) );
    BOOST_CHECK( std::isfinite( slab_scatter[i] ) );
    BOOST_CHECK( point_scatter[i] >= 0.0 );
    BOOST_CHECK( slab_scatter[i] >= 0.0 );
    point_sum += point_scatter[i];
    slab_sum += slab_scatter[i];
    abs_diff += std::abs( point_scatter[i] - slab_scatter[i] );
  }

  BOOST_CHECK_MESSAGE( point_sum > 0.0 || slab_sum > 0.0,
    "Both point and slab scatter are zero, table likely empty" );
  if( point_sum > 1e-6 && slab_sum > 1e-6 )
  {
    const double rel = abs_diff / ( point_sum + slab_sum );
    BOOST_CHECK_MESSAGE( rel > 1e-3,
      "Point and slab geometry give identical scatter -- shape blending broken" );
  }
}


// Sanity: zero areal density should give zero scatter.
BOOST_AUTO_TEST_CASE( ZeroShielding )
{
  ensure_loaded();

  BOOST_REQUIRE( g_new );

  vector<double> scatter;
  g_new->computeShieldScatter( 661.7, 26.0, 0.0, 0.0, 0.0, scatter );

  BOOST_REQUIRE_EQUAL( scatter.size(), static_cast<size_t>( g_new->groupCount() ) );
  for( double v : scatter )
    BOOST_CHECK_EQUAL( v, 0.0 );

  // Same via getContinuum: scatter answer should be all zero, transmission == intensity.
  vector<float> binning( 1024 );
  for( size_t i = 0; i < 1024; ++i )
    binning[i] = ( 800.0f * static_cast<float>( i ) ) / 1023.0f;

  vector<float> ans;
  const float trans = g_new->getContinuum( ans, 661.7f, 1.0f, 26.0f, 0.0f, 0.0f, binning );

  BOOST_CHECK_CLOSE( trans, 1.0f, 0.001f );
  for( float v : ans )
    BOOST_CHECK_EQUAL( v, 0.0f );
}
