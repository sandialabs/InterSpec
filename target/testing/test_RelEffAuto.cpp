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

#include "rapidxml/rapidxml.hpp"

//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE RelEffManual_suite
//#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/RelActCalc.h"

#include "InterSpec/RelActCalcAuto_imp.hpp"


using namespace std;
using namespace boost::unit_test;

static_assert( USE_REL_ACT_TOOL, "Compile-time support for Rel Act tool is required." );

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



BOOST_AUTO_TEST_CASE( FitContinuum )
{
  //set_data_dir();
  
  
  vector<PeakDef> peaks;
  peaks.emplace_back( 103.5, 2.5, 150.0 );

  constexpr size_t nbin = 7;
  constexpr float energies[nbin+1] = {100.0f, 101.0f, 102.0f, 103.0f, 104.0f, 105.0f, 106.0f, 107.0f};
  constexpr float data[nbin] = {900.0f, 1090.0f, 990.0f, 1090.0f, 910.0f, 1090.0f, 950.0f};
  const PeakContinuum::OffsetType offset_type = PeakContinuum::OffsetType::Linear;
  const int num_polynomial_terms = PeakContinuum::num_parameters( offset_type );
  const bool step_continuum = PeakContinuum::is_step_continuum( offset_type );
  constexpr double ref_energy = energies[0];
  vector<double> continuum_coeffs(num_polynomial_terms, 0.0);
  double peak_counts[nbin];
    
  RelActCalcAuto::fit_continuum( energies, data, nbin, num_polynomial_terms, step_continuum,
                                  ref_energy, peaks, continuum_coeffs.data(), peak_counts );
    
    
  vector<double> dummy_amps, continuum_coeffs_old, dummy_amp_uncert, continuum_uncerts;
  fit_amp_and_offset( energies, data, nbin, num_polynomial_terms,
                       step_continuum, ref_energy, {}, {}, peaks,
                      PeakDef::SkewType::NoSkew, nullptr, dummy_amps,
                       continuum_coeffs_old, dummy_amp_uncert, continuum_uncerts );
    
  BOOST_REQUIRE( continuum_coeffs_old.size() == static_cast<size_t>(num_polynomial_terms) );
  for( size_t i = 0; i < num_polynomial_terms; ++i )
  {
    const double new_coef = continuum_coeffs[i];
    const double old_coef = continuum_coeffs_old[i];
    const double diff = fabs(new_coef - old_coef);
    BOOST_CHECK( (diff < 0.00001*std::max(fabs(new_coef), fabs(old_coef))) || (diff < 1.0E-12) );
  }
  
  double old_way_peak_counts[nbin] = { 0.0 };
  
  for( const PeakDef &peak : peaks )
    peak.gauss_integral( energies, old_way_peak_counts, nbin );

  std::shared_ptr<PeakContinuum> continuum = peaks[0].continuum();
  BOOST_REQUIRE( !!continuum );
  continuum->setType( offset_type );
  continuum->setParameters( ref_energy, continuum_coeffs, {} );

  for( size_t bin = 0; bin < nbin; ++bin )
    old_way_peak_counts[bin] += continuum->offset_integral(energies[bin], energies[bin+1], nullptr );

  for( size_t bin = 0; bin < nbin; ++bin )
  {
    const double new_val = peak_counts[bin];
    const double old_coef = old_way_peak_counts[bin];
    const double diff = fabs(new_val - old_coef);
    BOOST_CHECK( (diff < 0.00001*std::max(fabs(new_val), fabs(old_coef))) || (diff < 1.0E-12) );
  }

}//BOOST_AUTO_TEST_CASE( FitRelActManualToKnown )
