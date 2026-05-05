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
#include "InterSpec/GadrasSpecFunc.h"
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


  // Sum the total scatter intensity in `bins`.
  float sum_bins( const vector<float> &bins )
  {
    float s = 0.0f;
    for( float v : bins )
      s += v;
    return s;
  }


  bool any_bad( const vector<float> &bins )
  {
    for( float v : bins )
    {
      if( !std::isfinite( v ) || v < 0.0f )
        return true;
    }
    return false;
  }
}//anonymous namespace


// One global pair, loaded once.
static unique_ptr<GadrasScatterTable> g_old;
static unique_ptr<GadrasShieldScatter> g_new;
static bool g_loaded = false;


static void ensure_loaded()
{
  if( g_loaded )
    return;
  g_loaded = true;

  const string data_dir_old = find_data_dir( "GadrasContinuum.lib" );
  const string data_dir_new = find_data_dir( "sandia.shieldscatter.db" );
  const string old_file = SpecUtils::append_path( data_dir_old, "GadrasContinuum.lib" );
  const string new_file = SpecUtils::append_path( data_dir_new, "sandia.shieldscatter.db" );

  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( old_file ), "Missing " << old_file );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( new_file ), "Missing " << new_file );

  // The new wrapper needs MassAttenuationTool, which needs the static data dir.
  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( data_dir_old ) );

  BOOST_REQUIRE_NO_THROW( g_old.reset( new GadrasScatterTable( old_file ) ) );
  BOOST_REQUIRE_NO_THROW( g_new.reset( new GadrasShieldScatter( new_file ) ) );
}


BOOST_AUTO_TEST_CASE( LoadFiles )
{
  ensure_loaded();
  BOOST_REQUIRE( g_old );
  BOOST_REQUIRE( g_new );
  BOOST_CHECK( g_new->groupCount() > 0 );
}


// Compare uncollided transmission and total scatter integral for a sweep over
// (material, areal density, source energy). Both implementations share the
// same TransmissionH() helper internally, so transmission should match very
// closely; the scatter algorithms differ substantially so we accept order-of-
// magnitude agreement.
BOOST_AUTO_TEST_CASE( CompareAirFePb )
{
  ensure_loaded();

  BOOST_REQUIRE( g_old );
  BOOST_REQUIRE( g_new );
  
  struct Mat
  {
    const char *name;
    float Z;
    vector<float> ads;  // g/cm^2
  };

  const Mat materials[] = {
    { "Air", 7.36f, { 0.1f, 0.5f, 2.0f } },
    { "Fe",  26.0f, { 1.0f, 5.0f, 20.0f, 50.0f } },
    { "Pb",  82.0f, { 1.0f, 5.0f, 20.0f, 50.0f } }
  };

  const float energies[] = { 60.0f, 122.0f, 356.0f, 661.7f, 1173.0f, 1332.0f, 2614.0f };

  size_t cells_compared = 0;
  for( const Mat &m : materials )
  {
    for( float ad : m.ads )
    {
      for( float E : energies )
      {
        // Build a 1024-bin linear binning to maxE (the same shape DoseCalc uses).
        const float maxE = E + 100.0f;
        const size_t nbins = 1024;
        vector<float> binning( nbins );
        for( size_t i = 0; i < nbins; ++i )
          binning[i] = ( maxE * static_cast<float>( i ) ) / static_cast<float>( nbins - 1 );

        vector<float> ans_old, ans_new;
        const float trans_old = g_old->getContinuum( ans_old, E, 1.0f, m.Z, ad, 0.0f, binning );
        const float trans_new = g_new->getContinuum( ans_new, E, 1.0f, m.Z, ad, 0.0f, binning );

        BOOST_CHECK_MESSAGE( !any_bad( ans_old ),
          "Old scatter has NaN/neg for " << m.name << " AD=" << ad << " E=" << E );
        BOOST_CHECK_MESSAGE( !any_bad( ans_new ),
          "New scatter has NaN/neg for " << m.name << " AD=" << ad << " E=" << E );

        // Transmission: should match within 2% relative; both share the same
        // TransmissionH helper so this would only fail if one implementation
        // changed how it scales the return value.
        if( trans_old > 1e-12f && trans_new > 1e-12f )
        {
          const float rel = std::abs( trans_old - trans_new )
                            / std::max( trans_old, trans_new );
          BOOST_CHECK_MESSAGE( rel < 0.02f,
            "Transmission disagree for " << m.name << " AD=" << ad << " E=" << E
              << ": old=" << trans_old << " new=" << trans_new );
        }

        // Scatter integral: algorithm differs significantly between the two
        // (different table, different energy mesh, different geometry model).
        // We just require they're in the same order of magnitude when both
        // produce nonzero output, and that neither is wildly larger than the
        // unattenuated source intensity (1.0).
        const float sum_old = sum_bins( ans_old );
        const float sum_new = sum_bins( ans_new );

        BOOST_CHECK_MESSAGE( sum_old < 10.0f,
          "Old scatter integral too large: " << m.name << " AD=" << ad
            << " E=" << E << " sum=" << sum_old );
        BOOST_CHECK_MESSAGE( sum_new < 10.0f,
          "New scatter integral too large: " << m.name << " AD=" << ad
            << " E=" << E << " sum=" << sum_new );

        if( sum_old > 1e-6f && sum_new > 1e-6f )
        {
          const float ratio = sum_new / sum_old;
          BOOST_CHECK_MESSAGE( ratio > 0.05f && ratio < 20.0f,
            "Scatter integrals differ by >20x for " << m.name << " AD=" << ad
              << " E=" << E << ": old=" << sum_old << " new=" << sum_new );
        }

        ++cells_compared;
      }//for energies
    }//for ads
  }//for materials

  BOOST_CHECK( cells_compared > 0 );
  BOOST_TEST_MESSAGE( "Compared " << cells_compared << " (material, AD, energy) cells." );
}


// Sanity: shapeFactor extremes (point=0, slab=1) should both produce
// non-negative finite results, and they should differ from each other for
// nontrivial shielding (otherwise the geometry knob does nothing).
BOOST_AUTO_TEST_CASE( ShapeFactorEndpoints )
{
  ensure_loaded();

  BOOST_REQUIRE( g_old );
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

  BOOST_REQUIRE( g_old );
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
