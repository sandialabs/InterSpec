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

#include "InterSpec/RelActCalcAuto_imp.hpp"

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
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/DecayDataBaseServer.h"

#include "InterSpec/PeakFit_imp.hpp"


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
  const int num_polynomial_terms = static_cast<int>( PeakContinuum::num_parameters( offset_type ) );
  const bool step_continuum = PeakContinuum::is_step_continuum( offset_type );
  constexpr double ref_energy = energies[0];
  vector<double> continuum_coeffs(num_polynomial_terms, 0.0);
  double peak_counts[nbin];
    
  PeakFit::fit_continuum( energies, data, nullptr, nbin, num_polynomial_terms, step_continuum,
                                  ref_energy, peaks, false, continuum_coeffs.data(), peak_counts );


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


BOOST_AUTO_TEST_CASE( test_pu242_by_correlation )
{
  // We will first roughly test PuCorrMethod::ByPu239Only to data given in paper.
  //
  // Fig 3 in Swinhoe 2010 gives Pu239 content vs Pu242 content; I manually
  //  extracted the following values from the fit line in the PDF.
  const vector<pair<double,double>> swinhoe_approx_fig_3_data = {
    {0.55496, 0.06998},
    {0.55894, 0.06821},
    {0.56438, 0.06619},
    {0.57073, 0.06399},
    {0.57829, 0.06154},
    {0.58453, 0.05952},
    {0.59068, 0.05756},
    {0.59471, 0.05634},
    {0.59844, 0.05524},
    {0.60146, 0.05444},
    {0.60579, 0.05322},
    {0.61526, 0.05077},
    {0.62101, 0.04906},
    {0.62605, 0.04802},
    {0.63088, 0.04691},
    {0.63542, 0.04594},
    {0.6402, 0.0449}
  };//swinhoe_approx_fig_3_data
  
  
  for( const auto x_y : swinhoe_approx_fig_3_data )
  {
    const double x = x_y.first;
    const double y = x_y.second;
    const double gamma_spec_pu239 = x / (1.0 - y);
    const double gamma_spec_pu_other = (1.0 - x - y)/(1.0 - y);
    
    // gamma_spec_pu239 plus gamma_spec_pu_other should sum to 1.0
    BOOST_CHECK( fabs(1.0 - (gamma_spec_pu239 + gamma_spec_pu_other)) < 0.001 );
    
    RelActCalc::Pu242ByCorrelationInput input;
    input.pu_age = 0.0;
    input.pu238_rel_mass = gamma_spec_pu_other;
    input.pu239_rel_mass = gamma_spec_pu239;
    // Pu240, and Pu241/Am241 are irrelevant, all that
    
    RelActCalc::Pu242ByCorrelationOutput output = RelActCalc::correct_pu_mass_fractions_for_pu242( input, RelActCalc::PuCorrMethod::ByPu239Only );
    
    //cout << "For Swinhoe [" << x << ", " << y << "]: Pu239: " << output.pu239_mass_frac
    //     << ", Pu242: " << output.pu242_mass_frac << " +- " << 100.0*output.pu242_uncert << "%\n";
    
    BOOST_CHECK( fabs(output.pu239_mass_frac - x) < 0.005 );
    BOOST_CHECK( fabs(output.pu242_mass_frac - y) < 0.0005 );
  }//for( const auto x_y : swinhoe_approx_fig_3_data )
  
  
  // For PuCorrMethod::Bignan95_BWR and PuCorrMethod::Bignan95_PWR, we dont have nearly as good
  //  of comparison data
  RelActCalc::Pu242ByCorrelationInput input;
  input.pu_age = 0.0;
  input.pu238_rel_mass = 0.0120424;
  input.pu239_rel_mass = 0.6649628;
  input.pu240_rel_mass = 0.2327493;
  input.pu241_rel_mass = 0.0501864;
  //input.pu241_rel_mass = 0.0361259;
  RelActCalc::Pu242ByCorrelationOutput output = RelActCalc::correct_pu_mass_fractions_for_pu242( input, RelActCalc::PuCorrMethod::Bignan95_BWR );
  //cout << "For Bignan95_BWR: Pu239: " << output.pu239_mass_frac
  //     << ", Pu242: " << output.pu242_mass_frac << " +- " << 100.0*output.pu242_uncert << "%\n";
  //For Bignan95_BWR: Pu239: 0.668767, Pu242: 0.0345679 +- 14%
  BOOST_CHECK( fabs(output.pu239_mass_frac - 0.668767) < 0.0001 );
  BOOST_CHECK( fabs(output.pu242_mass_frac - 0.0345679) < 0.0001 );
  BOOST_CHECK( fabs(output.pu242_uncert - 0.14) < 0.0001 );
  
  output = RelActCalc::correct_pu_mass_fractions_for_pu242( input, RelActCalc::PuCorrMethod::Bignan95_PWR );
  //cout << "For Bignan95_PWR: Pu239: " << output.pu239_mass_frac
  //     << ", Pu242: " << output.pu242_mass_frac << " +- " << 100.0*output.pu242_uncert << "%\n";
  //For Bignan95_PWR: Pu239: 0.664565, Pu242: 0.0406335 +- 6%
  BOOST_CHECK( fabs(output.pu239_mass_frac - 0.664565) < 0.0001 );
  BOOST_CHECK( fabs(output.pu242_mass_frac - 0.0406335) < 0.0001 );
  BOOST_CHECK( fabs(output.pu242_uncert - 0.06) < 0.0001 );
}//BOOST_AUTO_TEST_CASE( test_pu242_by_correlation )
