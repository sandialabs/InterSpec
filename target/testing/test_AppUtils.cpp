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

#include <string>
#include <iostream>

#include <Wt/Utils>

//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE QRSpectrum_suite
//#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


#include "InterSpec/AppUtils.h"



using namespace std;
using namespace boost::unit_test;


BOOST_AUTO_TEST_CASE( SplitUri )
{
  string host, path, query_str, fragment;

  string url = "interspec://decay/chain?nuclide=U238&activity=3uCi&initialage=20y&timespan=22y&actunits=ci";
  AppUtils::split_uri( url, host, path, query_str, fragment );
  BOOST_CHECK( host == "decay" );
  BOOST_CHECK( path == "chain" );
  BOOST_CHECK( query_str == "nuclide=U238&activity=3uCi&initialage=20y&timespan=22y&actunits=ci" );
  BOOST_CHECK( fragment == "" );

  url = "interspec://decay/chart?nuc=Ba133&act=3uCi&nuc=Cs137&act=2uci&actunits=ci&timespan=20y#somefragment";
  AppUtils::split_uri( url, host, path, query_str, fragment );
  BOOST_CHECK( host == "decay" );
  BOOST_CHECK( path == "chart" );
  BOOST_CHECK( query_str == "nuc=Ba133&act=3uCi&nuc=Cs137&act=2uci&actunits=ci&timespan=20y" );
  BOOST_CHECK( fragment == "somefragment" );


  url = "RADDATA://G0/000/NCF60UV/PGV2ER2*TFY1GQIN%ZQLOQMTVN4S0ZV6QDIQDJ$D...";
  AppUtils::split_uri( url, host, path, query_str, fragment );
  BOOST_CHECK( host == "G0" );
  BOOST_CHECK( path == "000/000/NCF60UV/PGV2ER2*TFY1GQIN%ZQLOQMTVN4S0ZV6QDIQDJ$D..." );
  BOOST_CHECK( query_str == "" );
  BOOST_CHECK( fragment == "" );

  url = "RADDATA://G0/000/NCF60UV/#";
  AppUtils::split_uri( url, host, path, query_str, fragment );
  BOOST_CHECK( host == "G0" );
  BOOST_CHECK( path == "000/NCF60UV/" );
  BOOST_CHECK( query_str == "" );
  BOOST_CHECK( fragment == "" );

  url = "RADDATA://G0/000/NCF60UV/?#";
  AppUtils::split_uri( url, host, path, query_str, fragment );
  BOOST_CHECK( host == "G0" );
  BOOST_CHECK( path == "000/NCF60UV/" );
  BOOST_CHECK( query_str == "" );
  BOOST_CHECK( fragment == "" );

  url = "RADDATA://G0/000/NCF6#0UV/?";
  AppUtils::split_uri( url, host, path, query_str, fragment );
  BOOST_CHECK( host == "G0" );
  BOOST_CHECK( path == "000/NCF6#0UV/" );
  BOOST_CHECK( query_str == "" );
  BOOST_CHECK( fragment == "" );


  url = "RADDATA://G0/000/NCF6#0UV/?query";
  AppUtils::split_uri( url, host, path, query_str, fragment );
  BOOST_CHECK( host == "G0" );
  BOOST_CHECK( path == "000/NCF6" );
  BOOST_CHECK( query_str == "" );
  BOOST_CHECK( fragment == "0UV/?query" );

  url = "RADDATA://G0/000/NCF6#frag";
  AppUtils::split_uri( url, host, path, query_str, fragment );
  BOOST_CHECK( host == "G0" );
  BOOST_CHECK( path == "000/NCF6#0UV/" );
  BOOST_CHECK( query_str == "" );
  BOOST_CHECK( fragment == "frag" );

  url = "interspec://decay/chart#somefragment?notaquery=3";
  AppUtils::split_uri( url, host, path, query_str, fragment );
  BOOST_CHECK( host == "decay" );
  BOOST_CHECK( path == "chart" );
  BOOST_CHECK( query_str == "" );
  BOOST_CHECK( fragment == "somefragment?notaquery=3" );

  url = "interspec://decay/chart#?notaquery=3";
  AppUtils::split_uri( url, host, path, query_str, fragment );
  BOOST_CHECK( host == "decay" );
  BOOST_CHECK( path == "chart" );
  BOOST_CHECK( query_str == "" );
  BOOST_CHECK( fragment == "?notaquery=3" );


  url = "interspec://decay/chart?#frag3";
  AppUtils::split_uri( url, host, path, query_str, fragment );
  BOOST_CHECK( host == "decay" );
  BOOST_CHECK( path == "chart" );
  BOOST_CHECK( query_str == "" );
  BOOST_CHECK( fragment == "frag3" );

  url = "interspec://decay/chart?#";
  AppUtils::split_uri( url, host, path, query_str, fragment );
  BOOST_CHECK( host == "decay" );
  BOOST_CHECK( path == "chart" );
  BOOST_CHECK( query_str == "" );
  BOOST_CHECK( fragment == "" );


  url = "interspec://decay?#";
  AppUtils::split_uri( url, host, path, query_str, fragment );
  BOOST_CHECK( host == "decay" );
  BOOST_CHECK( path == "" );
  BOOST_CHECK( query_str == "" );
  BOOST_CHECK( fragment == "" );

  url = "interspec://decay?#frag";
  AppUtils::split_uri( url, host, path, query_str, fragment );
  BOOST_CHECK( host == "decay" );
  BOOST_CHECK( path == "" );
  BOOST_CHECK( query_str == "" );
  BOOST_CHECK( fragment == "frag" );

  url = "interspec://decay?p=1#frag";
  AppUtils::split_uri( url, host, path, query_str, fragment );
  BOOST_CHECK( host == "decay" );
  BOOST_CHECK( path == "" );
  BOOST_CHECK( query_str == "p=1" );
  BOOST_CHECK( fragment == "frag" );


  url = "interspec://decay/?p=1#frag";
  AppUtils::split_uri( url, host, path, query_str, fragment );
  BOOST_CHECK( host == "decay" );
  BOOST_CHECK( path == "" );
  BOOST_CHECK( query_str == "p=1" );
  BOOST_CHECK( fragment == "frag" );


/*
  url = "";
  AppUtils::split_uri( url, host, path, query_str );
  BOOST_CHECK( host == "" );
  BOOST_CHECK( path == "" );
  BOOST_CHECK( query_str == "" );


  url = "";
  AppUtils::split_uri( url, host, path, query_str );
  BOOST_CHECK( host == "" );
  BOOST_CHECK( path == "" );
  BOOST_CHECK( query_str == "" );
*/
}//BOOST_AUTO_TEST_CASE( SimpleSpecEncode )

