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
#include "InterSpec/RelActManualGui.h"
#include "InterSpec/RelActCalcManual.h"
#include "InterSpec/RelEffShieldWidget.h"
#include "InterSpec/DecayDataBaseServer.h"


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


BOOST_AUTO_TEST_CASE( FitRelActManualToKnown )
{
  set_data_dir();
  
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );
  
  MaterialDB matdb;
  
  const string materialfile = SpecUtils::append_path( InterSpec::staticDataDirectory(), "MaterialDataBase.txt" );
  BOOST_REQUIRE_NO_THROW( matdb.parseGadrasMaterialFile( materialfile, db, false ) );
  
  const Material * const iron = matdb.material("Fe (iron)");
  BOOST_REQUIRE( iron );
  
  const string test_n42_file = SpecUtils::append_path(g_test_file_dir, "manual_rel_eff/spec184_235U_12.9543.n42");
  BOOST_REQUIRE( SpecUtils::is_file(test_n42_file) );
  
  
  SpecMeas infile;
  const bool loaded = infile.load_file( test_n42_file, SpecUtils::ParserType::Auto );
  BOOST_REQUIRE( loaded );
  
  std::shared_ptr<const DetectorPeakResponse> det = infile.detector();
  BOOST_REQUIRE( det );
  
  shared_ptr<deque<shared_ptr<const PeakDef>>> orig_peaks = infile.peaks( {1} );
  BOOST_REQUIRE( orig_peaks && orig_peaks->size() );

  BOOST_REQUIRE( infile.detector_names().size() == 1 );
  const shared_ptr<const SpecUtils::Measurement> meas = infile.measurement( 1, infile.detector_names()[0] );
  BOOST_REQUIRE( meas );
  
  const rapidxml::xml_document<char> * const guiStateXml = infile.relActManualGuiState();
  BOOST_REQUIRE( guiStateXml );
  
  
  RelActManualGui::RelActCalcRawInput calc_raw_input;
  auto guiState = make_shared<RelActManualGui::GuiState>();
  BOOST_REQUIRE_NO_THROW( guiState->deSerialize( guiStateXml->first_node() ) );
  BOOST_REQUIRE( guiState->m_relEffEqnFormIndex == RelActCalc::RelEffEqnForm::FramPhysicalModel ); // Just to make sure
  
  calc_raw_input.state = guiState;
  calc_raw_input.fore_spec = meas;
  //calc_raw_input.back_spec = ...;
  calc_raw_input.peaks = *orig_peaks;
  calc_raw_input.detector = det;
  
  RelActCalcManual::RelEffInput calc_input;
  
  BOOST_REQUIRE_NO_THROW( RelActManualGui::prepare_calc_input( calc_raw_input, &matdb, calc_input ) );
  BOOST_CHECK( calc_input.eqn_form == guiState->m_relEffEqnFormIndex );
  
  // We could check that we see U235, U238, U234, and U232 in `calc_input.peaks`
  
  const RelActCalcManual::RelEffSolution solution = RelActCalcManual::solve_relative_efficiency( calc_input );
  
  BOOST_CHECK( solution.m_status == RelActCalcManual::ManualSolutionStatus::Success );

  const double u235_mass_frac = solution.mass_fraction( "U235" );
  
  BOOST_CHECK( (u235_mass_frac > 0.11) && (u235_mass_frac < 0.15) ); //Truth value is 12.9543, we are getting 0.11204
  // We could check uncertainties
  
  // Lets try using a different equation form
  guiState->m_physModelUseHoerl = false;
  guiState->m_selfAttenShield.reset();
  guiState->m_externalShields.clear();
  guiState->m_relEffEqnFormIndex = RelActCalc::RelEffEqnForm::LnX;
  guiState->m_relEffEqnOrderIndex = 4;
  guiState->m_addUncertIndex = RelActManualGui::AddUncert::FivePercent;
  
  RelActCalcManual::RelEffInput calc_input_lnx;
  BOOST_REQUIRE_NO_THROW( RelActManualGui::prepare_calc_input( calc_raw_input, &matdb, calc_input_lnx ) );
  BOOST_REQUIRE( calc_input_lnx.eqn_form == RelActCalc::RelEffEqnForm::LnX );
  
  // First lets try using least linear squares for relative eff. equation
  calc_input_lnx.use_ceres_to_fit_eqn = false;
  const RelActCalcManual::RelEffSolution lnx_lls_solution = RelActCalcManual::solve_relative_efficiency( calc_input_lnx );
  BOOST_CHECK( lnx_lls_solution.m_status == RelActCalcManual::ManualSolutionStatus::Success );
  
  // First lets try using Ceres to get relative eff. equation
  calc_input_lnx.use_ceres_to_fit_eqn = true;
  const RelActCalcManual::RelEffSolution lnx_ceres_solution = RelActCalcManual::solve_relative_efficiency( calc_input_lnx );
  BOOST_CHECK( lnx_ceres_solution.m_status == RelActCalcManual::ManualSolutionStatus::Success );
  
  const double lls_enrich = lnx_lls_solution.mass_fraction( "U235" );
  const double ceres_enrich = lnx_ceres_solution.mass_fraction( "U235" );
  BOOST_CHECK( (lls_enrich > 0.1) && (u235_mass_frac < 0.15) ); //Actual value 12.9543; this method seems to give 0.1085
  
  BOOST_CHECK( fabs(lls_enrich - ceres_enrich) < 0.0001 );
}//BOOST_AUTO_TEST_CASE( FitRelActManualToKnown )
