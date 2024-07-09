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

#ifdef _WIN32
// For some reason, we need to include the following includes, before unit_test.hpp,
//  or we get a bunch of errors relating to winsock.k being included multiple times,
//  or something
#include "winsock2.h"
#include "Windows.h"
#endif

//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ShieldingSourceFitCalc_suite
//#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"


#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectionLimitCalc.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using namespace boost::unit_test;


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


BOOST_AUTO_TEST_CASE( BasicCurrieMda )
{
  using namespace DetectionLimitCalc;
  
  set_data_dir();
  
  auto spectrum = make_shared<SpecUtils::Measurement>();
  float liveTime = 1.0, realTime = 1.0;
  size_t num_channels = 1024;
  auto counts = make_shared<vector<float>>(num_channels, 10.0f);
  spectrum->set_gamma_counts( counts, liveTime, realTime );
  
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( num_channels, {0.0f, 3.0f}, {} );
  spectrum->set_energy_calibration( cal );
  
  CurrieMdaInput input;
  input.spectrum = spectrum;
  input.gamma_energy = 512.0f;     //channel: 170.667
  input.roi_lower_energy = 500.0f; //channel: 166.667
  input.roi_upper_energy = 524.0f; //channel: 174.667
  input.num_lower_side_channels = 4;
  input.num_upper_side_channels = 4;
  input.detection_probability = 0.95;
  input.additional_uncertainty = 0.0;
  
  CurrieMdaResult result;
  
  BOOST_REQUIRE_NO_THROW( result = currie_mda_calc( input ) );
  
  BOOST_CHECK_EQUAL( result.first_lower_continuum_channel, 163 );
  BOOST_CHECK_EQUAL( result.last_lower_continuum_channel, 166 );
  BOOST_CHECK_EQUAL( result.lower_continuum_counts_sum, 40.0f );
  BOOST_CHECK_EQUAL( result.first_upper_continuum_channel, 175 );
  BOOST_CHECK_EQUAL( result.last_upper_continuum_channel, 178 );
  BOOST_CHECK_EQUAL( result.upper_continuum_counts_sum, 40.0f );
  BOOST_CHECK_EQUAL( result.first_peak_region_channel, 167 );
  BOOST_CHECK_EQUAL( result.last_peak_region_channel, 174 );
  BOOST_CHECK_EQUAL( result.peak_region_counts_sum, 80.0f );
  BOOST_CHECK_CLOSE_FRACTION( result.continuum_eqn[0], 3.3333333333333335, 1.0E-5 );
  BOOST_CHECK_LT( fabs(result.continuum_eqn[1]), 1.0E-5 );
  BOOST_CHECK_EQUAL( result.estimated_peak_continuum_counts, 80.0f );
  BOOST_CHECK_CLOSE_FRACTION( result.estimated_peak_continuum_uncert, sqrt(80.0f), 1.0E-5 ); //4 side channels each side
  
  const double k = 1.64485;
  const double region_sigma = sqrt( 80 + 80 );
  BOOST_CHECK_CLOSE_FRACTION( result.decision_threshold, k*region_sigma, 1.0E-5 ); //continuum uncert^2 + expected amount uncert^2
  BOOST_CHECK_CLOSE_FRACTION( result.detection_limit, (2.0*k*region_sigma + k*k), 1.0E-5 );
  BOOST_CHECK_EQUAL( result.source_counts, 0.0f );
  BOOST_CHECK_CLOSE_FRACTION( result.lower_limit, -k*region_sigma, 1.0E-5 );
  BOOST_CHECK_CLOSE_FRACTION( result.upper_limit, k*region_sigma, 1.0E-5 );
}//BOOST_AUTO_TEST_CASE( BasicCurrieMda )


BOOST_AUTO_TEST_CASE( Table16OfAQ48 )
{
  //Example given on page 68 of:
  //  INTERNATIONAL ATOMIC ENERGY AGENCY, Determination and Interpretation of
  //  Characteristic Limits for Radioactivity Measurements, IAEA Analytical Quality
  //  in Nuclear Applications Series No. 48, IAEA, Vienna (2017)
  using namespace DetectionLimitCalc;
  
  set_data_dir();
  
  auto spectrum = make_shared<SpecUtils::Measurement>();
  float liveTime = 1.0, realTime = 1.0;
  size_t num_channels = 1024;
  auto counts = make_shared<vector<float>>(num_channels, 0.0f);
  (*counts)[992]  = 59;
  (*counts)[993]  = 56;
  (*counts)[994]  = 74;
  (*counts)[995]  = 72;
  (*counts)[996]  = 64;
  (*counts)[997]  = 46; //First channel of ROI
  (*counts)[998]  = 64;
  (*counts)[999]  = 67;
  (*counts)[1000] = 82;
  (*counts)[1001] = 95;
  (*counts)[1002] = 157;
  (*counts)[1003] = 398;
  (*counts)[1004] = 807;
  (*counts)[1005] = 1480;
  (*counts)[1006] = 1814;
  (*counts)[1007] = 1936;
  (*counts)[1008] = 1575;
  (*counts)[1009] = 940;
  (*counts)[1010] = 457;
  (*counts)[1011] = 207;
  (*counts)[1012] = 82;
  (*counts)[1013] = 49;
  (*counts)[1014] = 50;
  (*counts)[1015] = 45;
  (*counts)[1016] = 43; //Last channel of ROI
  (*counts)[1017] = 45;
  (*counts)[1018] = 51;
  (*counts)[1019] = 35;
  (*counts)[1020] = 45;
  (*counts)[1021] = 50;
  
  spectrum->set_gamma_counts( counts, liveTime, realTime );
  
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_polynomial( num_channels, {0.0f, 1.0f}, {} );
  spectrum->set_energy_calibration( cal );
  
  CurrieMdaInput input;
  input.spectrum = spectrum;
  input.gamma_energy = 1007.0f;
  input.roi_lower_energy = 997.0f;
  input.roi_upper_energy = 1017.0f;
  input.num_lower_side_channels = 5;
  input.num_upper_side_channels = 5;
  input.detection_probability = 0.95;
  input.additional_uncertainty = 0.0;
  
  CurrieMdaResult result;
  
  BOOST_REQUIRE_NO_THROW( result = currie_mda_calc( input ) );
  
  BOOST_CHECK_EQUAL( result.peak_region_counts_sum, 10394 );
  BOOST_CHECK_EQUAL( result.lower_continuum_counts_sum, 325 );  //n_1
  BOOST_CHECK_EQUAL( result.upper_continuum_counts_sum, 226 );  //n_2
  BOOST_CHECK_EQUAL( result.lower_continuum_counts_sum + result.upper_continuum_counts_sum, 551 );  //n_0
  BOOST_CHECK_EQUAL( 1 + result.last_lower_continuum_channel - result.first_lower_continuum_channel, 5 );
  BOOST_CHECK_EQUAL( 1 + result.last_upper_continuum_channel - result.first_upper_continuum_channel, 5 );
  BOOST_CHECK_EQUAL( 1 + result.last_peak_region_channel - result.first_peak_region_channel, 20 ); //n_g - num channels in peak reagion
  BOOST_CHECK_EQUAL( result.estimated_peak_continuum_counts, 1102.0f );
  const double w = 0.002458; //weight: counts to bq
  //const double w_uncert_rel = 0.0478; //We dont take this into account, yet
  const double k = 1.64485;
  const double a_star = w * result.decision_threshold;
  BOOST_CHECK_CLOSE_FRACTION( a_star, 0.232506, 1.0E-3 ); //We get: 0.232467
  
  const double a_pound = w * result.detection_limit;
  BOOST_CHECK_CLOSE_FRACTION( a_pound, 0.471705, 1.0E-3 ); //We get: 0.471583
}//BOOST_AUTO_TEST_CASE( Table16OfAQ48 )
