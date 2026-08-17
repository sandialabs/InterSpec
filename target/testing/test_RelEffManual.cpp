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

#include <set>
#include <cmath>
#include <string>
#include <iostream>
#include <functional>

#include "rapidxml/rapidxml.hpp"

//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE RelEffManual_suite
//#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/PhysicalUnits.h"
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
  
  BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  BOOST_REQUIRE( MaterialDB::initialized() );

  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  BOOST_REQUIRE( matdb );

  const std::shared_ptr<const Material> iron = matdb->material( "Fe (iron)" );
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
  // prepare_calc_input expects peaks already filtered to nuclide/reaction-assigned, manual-rel-eff peaks
  //  (normally done by get_raw_info_for_calc_input()); the file's peak list can include others, so filter here.
  for( const shared_ptr<const PeakDef> &p : *orig_peaks )
    if( p && (p->parentNuclide() || p->reaction()) && p->useForManualRelEff() )
      calc_raw_input.peaks.push_back( p );
  calc_raw_input.detector = det;
  
  RelActCalcManual::RelEffInput calc_input;
  
  BOOST_REQUIRE_NO_THROW( RelActManualGui::prepare_calc_input( calc_raw_input, calc_input ) );
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
  BOOST_REQUIRE_NO_THROW( RelActManualGui::prepare_calc_input( calc_raw_input, calc_input_lnx ) );
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


// Regression test for fixed mass-fraction constraints in the manual solver.  A *fixed* (lower==upper)
// mass-fraction constraint removes its nuclide from the free fit; the solver previously left such a nuclide
// with the act-ratio -1.0 norm sentinel while its own asserts (functor ctor / mass_fraction()) expected 1.0 -
// contradictory invariants that aborted assert-enabled builds.  Fixed mass-fraction nuclides now carry norm
// 1.0 (matching range constraints and release behavior).  spec184 has fitted peaks for U232/U234/U235/U238,
// so pinning U234 (other U isotopes free) exercises the path.
BOOST_AUTO_TEST_CASE( FitRelActManualFixedMassFractionConstraint )
{
  set_data_dir();

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );

  const string test_n42_file = SpecUtils::append_path( g_test_file_dir, "manual_rel_eff/spec184_235U_12.9543.n42" );
  BOOST_REQUIRE( SpecUtils::is_file(test_n42_file) );

  SpecMeas infile;
  BOOST_REQUIRE( infile.load_file( test_n42_file, SpecUtils::ParserType::Auto ) );

  const std::shared_ptr<const DetectorPeakResponse> det = infile.detector();
  BOOST_REQUIRE( det );

  const shared_ptr<deque<shared_ptr<const PeakDef>>> orig_peaks = infile.peaks( {1} );
  BOOST_REQUIRE( orig_peaks && orig_peaks->size() );

  BOOST_REQUIRE( infile.detector_names().size() == 1 );
  const shared_ptr<const SpecUtils::Measurement> meas = infile.measurement( 1, infile.detector_names()[0] );
  BOOST_REQUIRE( meas );

  const rapidxml::xml_document<char> * const guiStateXml = infile.relActManualGuiState();
  BOOST_REQUIRE( guiStateXml );

  RelActManualGui::RelActCalcRawInput calc_raw_input;
  auto guiState = make_shared<RelActManualGui::GuiState>();
  BOOST_REQUIRE_NO_THROW( guiState->deSerialize( guiStateXml->first_node() ) );
  calc_raw_input.state = guiState;
  calc_raw_input.fore_spec = meas;
  // prepare_calc_input expects peaks already filtered to nuclide/reaction-assigned, manual-rel-eff peaks
  //  (normally done by get_raw_info_for_calc_input()); the file's peak list can include others, so filter here.
  for( const shared_ptr<const PeakDef> &p : *orig_peaks )
    if( p && (p->parentNuclide() || p->reaction()) && p->useForManualRelEff() )
      calc_raw_input.peaks.push_back( p );
  calc_raw_input.detector = det;

  RelActCalcManual::RelEffInput calc_input;
  BOOST_REQUIRE_NO_THROW( RelActManualGui::prepare_calc_input( calc_raw_input, calc_input ) );

  // The constraint's m_specific_activities must list every nuclide of the constrained element present in the
  //  fit (including the constrained one), so gather the uranium isotopes actually in the problem's peaks.
  std::set<std::string> u_isotopes;
  for( const RelActCalcManual::GenericPeakInfo &peak : calc_input.peaks )
  {
    for( const RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas )
    {
      const SandiaDecay::Nuclide * const n = db->nuclide( line.m_isotope );
      if( n && (n->atomicNumber == 92) )
        u_isotopes.insert( line.m_isotope );
    }
  }
  BOOST_REQUIRE_MESSAGE( u_isotopes.count("U234"), "spec184 problem has no U234 peak to constrain." );
  BOOST_REQUIRE_MESSAGE( u_isotopes.size() >= 2, "Need >=2 U isotopes so at least one stays free." );

  // First solve UNCONSTRAINED, to learn U234's natural mass fraction (which we then pin) and the free-isotope
  //  result a consistent fixed-constraint solve should reproduce.
  RelActCalcManual::RelEffSolution unc_sol;
  BOOST_REQUIRE_NO_THROW( unc_sol = RelActCalcManual::solve_relative_efficiency( calc_input ) );
  BOOST_REQUIRE( unc_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );
  const double u234_natural = unc_sol.mass_fraction( "U234" );
  const double u235_unconstrained = unc_sol.mass_fraction( "U235" );
  BOOST_REQUIRE_MESSAGE( (u234_natural > 1.0e-4) && (u234_natural < 0.5),
                         "Unexpected unconstrained U234 mass fraction " << u234_natural );

  // Pin U234 at the value it already wants and re-solve; fixing a nuclide at its natural value must not move
  //  the free isotopes.  (Without the fix, this solve aborts in debug builds via contradictory norm asserts.)
  RelActCalcManual::RelEffInput constrained_input = calc_input;
  RelActCalcManual::MassFractionConstraint u234_fixed;
  u234_fixed.m_nuclide = "U234";
  u234_fixed.m_mass_fraction_lower = u234_natural;
  u234_fixed.m_mass_fraction_upper = u234_natural;
  for( const std::string &iso : u_isotopes )
    u234_fixed.m_specific_activities[iso] = db->nuclide(iso)->activityPerGram();
  constrained_input.mass_fraction_constraints.push_back( u234_fixed );

  RelActCalcManual::RelEffSolution con_sol;
  BOOST_REQUIRE_NO_THROW( con_sol = RelActCalcManual::solve_relative_efficiency( constrained_input ) );
  BOOST_CHECK_MESSAGE( con_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success,
                       "Manual solve with a fixed mass-fraction constraint did not succeed." );
  if( con_sol.m_status != RelActCalcManual::ManualSolutionStatus::Success )
    return;

  // U234 must decode back to its pinned value, and U235 must match the unconstrained result.
  BOOST_CHECK_MESSAGE( fabs(con_sol.mass_fraction("U234") - u234_natural) < 1.0e-3,
                       "Fixed U234 mass fraction " << con_sol.mass_fraction("U234")
                       << " is not pinned at its natural value " << u234_natural );
  BOOST_CHECK_MESSAGE( fabs(con_sol.mass_fraction("U235") - u235_unconstrained) < 0.02,
                       "U235 moved from " << u235_unconstrained << " (unconstrained) to "
                       << con_sol.mass_fraction("U235") << " when U234 was pinned at its natural value." );
}//BOOST_AUTO_TEST_CASE( FitRelActManualFixedMassFractionConstraint )


namespace
{
  /** Loads the spec184 test file (once) and prepares a `RelActCalcManual::RelEffInput` from its
   embedded "Isotopics from peaks" GUI state; `modify_state` lets a test tweak the deserialized
   GuiState (eqn form, shields, add. uncert., ...) before the input is prepared.
   */
  RelActCalcManual::RelEffInput spec184_calc_input(
                      const std::function<void(RelActManualGui::GuiState &)> &modify_state )
  {
    set_data_dir();

    if( !MaterialDB::initialized() )
    {
      BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
      BOOST_REQUIRE( MaterialDB::initialized() );
    }

    static std::shared_ptr<SpecMeas> s_infile;
    if( !s_infile )
    {
      const string test_n42_file = SpecUtils::append_path( g_test_file_dir, "manual_rel_eff/spec184_235U_12.9543.n42" );
      BOOST_REQUIRE( SpecUtils::is_file(test_n42_file) );

      std::shared_ptr<SpecMeas> infile = make_shared<SpecMeas>();
      BOOST_REQUIRE( infile->load_file( test_n42_file, SpecUtils::ParserType::Auto ) );
      s_infile = infile;
    }//if( !s_infile )

    const std::shared_ptr<const DetectorPeakResponse> det = s_infile->detector();
    BOOST_REQUIRE( det );

    const shared_ptr<deque<shared_ptr<const PeakDef>>> orig_peaks = s_infile->peaks( {1} );
    BOOST_REQUIRE( orig_peaks && orig_peaks->size() );

    BOOST_REQUIRE( s_infile->detector_names().size() == 1 );
    const shared_ptr<const SpecUtils::Measurement> meas = s_infile->measurement( 1, s_infile->detector_names()[0] );
    BOOST_REQUIRE( meas );

    const rapidxml::xml_document<char> * const guiStateXml = s_infile->relActManualGuiState();
    BOOST_REQUIRE( guiStateXml );

    auto guiState = make_shared<RelActManualGui::GuiState>();
    BOOST_REQUIRE_NO_THROW( guiState->deSerialize( guiStateXml->first_node() ) );

    modify_state( *guiState );

    RelActManualGui::RelActCalcRawInput calc_raw_input;
    calc_raw_input.state = guiState;
    calc_raw_input.fore_spec = meas;
    // prepare_calc_input expects peaks already filtered to nuclide/reaction-assigned, manual-rel-eff
    //  peaks (normally done by get_raw_info_for_calc_input()).
    for( const shared_ptr<const PeakDef> &p : *orig_peaks )
      if( p && (p->parentNuclide() || p->reaction()) && p->useForManualRelEff() )
        calc_raw_input.peaks.push_back( p );
    calc_raw_input.detector = det;

    RelActCalcManual::RelEffInput calc_input;
    BOOST_REQUIRE_NO_THROW( RelActManualGui::prepare_calc_input( calc_raw_input, calc_input ) );

    return calc_input;
  }//spec184_calc_input(...)


  /** Returns a GuiState modifier that switches to an empirical rel-eff form (no shields/Hoerl). */
  std::function<void(RelActManualGui::GuiState &)> empirical_state_mod(
              const RelActCalc::RelEffEqnForm form, const int order,
              const RelActManualGui::AddUncert add_uncert = RelActManualGui::AddUncert::FivePercent )
  {
    return [form, order, add_uncert]( RelActManualGui::GuiState &state )
    {
      state.m_physModelUseHoerl = false;
      state.m_selfAttenShield.reset();
      state.m_externalShields.clear();
      state.m_relEffEqnFormIndex = form;
      state.m_relEffEqnOrderIndex = order;
      state.m_addUncertIndex = add_uncert;
    };
  }//empirical_state_mod(...)
}//namespace


// Regression for the rel-eff uncertainty band with `use_ceres_to_fit_eqn = true` and an empirical
// equation form: RelEffSolution::rel_eff_eqn_uncert() used to have a fall-through switch where only
// the LnX case returned (LnY fell through into LnXLnY; FramEmpirical accumulated into a shadowed
// local that was then discarded), ending in a logic-error throw - so the GUI band silently
// disappeared for those forms.  All empirical forms now delegate to
// RelActCalc::eval_eqn_uncertainty().
BOOST_AUTO_TEST_CASE( CeresEmpiricalEqnUncertPerForm )
{
  using RelActCalc::RelEffEqnForm;

  for( const RelEffEqnForm form : { RelEffEqnForm::LnX, RelEffEqnForm::LnY,
                                    RelEffEqnForm::LnXLnY, RelEffEqnForm::FramEmpirical } )
  {
    RelActCalcManual::RelEffInput calc_input = spec184_calc_input( empirical_state_mod(form, 3) );
    calc_input.use_ceres_to_fit_eqn = true;

    RelActCalcManual::RelEffSolution solution;
    BOOST_REQUIRE_NO_THROW( solution = RelActCalcManual::solve_relative_efficiency( calc_input ) );
    BOOST_REQUIRE_MESSAGE( solution.m_status == RelActCalcManual::ManualSolutionStatus::Success,
                           "Ceres-mode solve failed for eqn form " << RelActCalc::to_str(form) );

    for( const double energy : { 150.0, 186.0, 400.0, 1001.0 } )
    {
      double value = -1.0, uncert = -1.0;
      BOOST_REQUIRE_NO_THROW( value = solution.rel_eff_eqn_value( energy ) );
      BOOST_CHECK_NO_THROW( uncert = solution.rel_eff_eqn_uncert( energy ) );
      BOOST_CHECK_MESSAGE( !std::isnan(uncert) && !std::isinf(uncert) && (uncert > 0.0),
                           "Invalid uncertainty band value " << uncert << " for eqn form "
                           << RelActCalc::to_str(form) << " at " << energy << " keV." );
      BOOST_CHECK_MESSAGE( uncert < 100.0*value,
                           "Implausibly large band " << uncert << " (curve value " << value
                           << ") for eqn form " << RelActCalc::to_str(form) << " at " << energy << " keV." );
    }//for( loop over energies )
  }//for( loop over empirical eqn forms )
}//BOOST_AUTO_TEST_CASE( CeresEmpiricalEqnUncertPerForm )


// Validates the physical-model uncertainty band against an independent finite-difference oracle.
// The band used to contract a hand-derived physical-space gradient (AD in internal PhysicalUnits,
// real atomic numbers) with the Ceres parameter-space covariance (AD as g/cm2 numbers, AN divided
// by ns_an_ceres_mult) - no chain rule applied - understating e.g. the AD contribution by
// g_per_cm2^2 ~ 3.9e9.  The band is now an auto-differentiated gradient in Ceres-parameter space;
// this test independently rebuilds the parameter->curve mapping and checks sqrt(grad^T*C*grad).
BOOST_AUTO_TEST_CASE( PhysicalModelBandMatchesFiniteDifference )
{
  RelActCalcManual::RelEffInput calc_input = spec184_calc_input( []( RelActManualGui::GuiState & ){} );
  BOOST_REQUIRE( calc_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel );

  RelActCalcManual::RelEffSolution solution;
  BOOST_REQUIRE_NO_THROW( solution = RelActCalcManual::solve_relative_efficiency( calc_input ) );
  BOOST_REQUIRE( solution.m_status == RelActCalcManual::ManualSolutionStatus::Success );
  BOOST_REQUIRE( !solution.m_rel_eff_eqn_covariance.empty() );

  const RelActCalcManual::RelEffInput &in = solution.m_input;
  const vector<double> &coefs = solution.m_rel_eff_eqn_coefficients;
  const vector<vector<double>> &cov = solution.m_rel_eff_eqn_covariance;
  BOOST_REQUIRE( cov.size() == coefs.size() );

  // Independent re-implementation of the Ceres-parameter -> curve mapping (see
  //  `make_phys_eqn_input(...)` in RelActCalcManual.cpp): an AN slot exists only for shields with
  //  no material and fit_atomic_number, carrying AN/ns_an_ceres_mult; AD slots hold the areal
  //  density as a g/cm2 number; the two trailing Hoerl slots are offset by
  //  ns_decay_hoerl_b_offset / ns_decay_hoerl_c_offset.
  const auto eval_curve = [&in]( const vector<double> &pars, const double energy ) -> double
  {
    size_t idx = 0;

    std::optional<RelActCalc::PhysModelShield<double>> self_atten;
    const std::shared_ptr<const RelActCalc::PhysicalModelShieldInput> &self_in = in.phys_model_self_atten;
    if( self_in && (self_in->material
                    || ((self_in->atomic_number >= 1.0) && (self_in->atomic_number <= 98.0))
                    || self_in->fit_atomic_number) )
    {
      RelActCalc::PhysModelShield<double> shield;
      shield.material = self_in->material;
      if( !shield.material )
        shield.atomic_number = self_in->fit_atomic_number
                                 ? pars.at(idx++) * RelActCalc::ns_an_ceres_mult
                                 : self_in->atomic_number;
      shield.areal_density = pars.at(idx++) * PhysicalUnits::g_per_cm2;
      self_atten = std::move(shield);
    }//if( there is a self-atten shield )

    vector<RelActCalc::PhysModelShield<double>> external_attens;
    for( const std::shared_ptr<const RelActCalc::PhysicalModelShieldInput> &ext_in : in.phys_model_external_attens )
    {
      if( !ext_in->material && ((ext_in->atomic_number < 1.0) || (ext_in->atomic_number > 98.0)) )
        continue;

      RelActCalc::PhysModelShield<double> shield;
      shield.material = ext_in->material;
      if( !shield.material )
        shield.atomic_number = ext_in->fit_atomic_number
                                 ? pars.at(idx++) * RelActCalc::ns_an_ceres_mult
                                 : ext_in->atomic_number;
      shield.areal_density = pars.at(idx++) * PhysicalUnits::g_per_cm2;
      external_attens.push_back( std::move(shield) );
    }//for( loop over external attenuators )

    std::optional<double> hoerl_b, hoerl_c;
    if( in.phys_model_use_hoerl )
    {
      hoerl_b = (pars.at(idx++) - RelActCalc::ns_decay_hoerl_b_offset) * RelActCalc::ns_decay_hoerl_b_multiple;
      hoerl_c = (pars.at(idx++) - RelActCalc::ns_decay_hoerl_c_offset) * RelActCalc::ns_decay_hoerl_c_multiple;
    }

    BOOST_REQUIRE( idx == pars.size() );

    return RelActCalc::eval_physical_model_eqn( energy, self_atten, external_attens,
                                                in.phys_model_detector.get(), hoerl_b, hoerl_c );
  };//eval_curve

  for( const double energy : { 122.0, 186.0, 400.0, 766.0, 1001.0 } )
  {
    vector<double> gradient( coefs.size(), 0.0 );
    for( size_t i = 0; i < coefs.size(); ++i )
    {
      const double step = 1.0e-6 * (std::max)( 1.0, fabs(coefs[i]) );
      vector<double> pars_plus = coefs, pars_minus = coefs;
      pars_plus[i] += step;
      pars_minus[i] -= step;
      gradient[i] = (eval_curve(pars_plus, energy) - eval_curve(pars_minus, energy)) / (2.0*step);
    }//for( loop over parameters )

    double uncert_sq = 0.0;
    for( size_t i = 0; i < coefs.size(); ++i )
      for( size_t j = 0; j < coefs.size(); ++j )
        uncert_sq += gradient[i] * cov[i][j] * gradient[j];

    const double expected = sqrt( (std::max)(uncert_sq, 0.0) );

    double band = -1.0;
    BOOST_REQUIRE_NO_THROW( band = solution.rel_eff_eqn_uncert( energy ) );

    BOOST_CHECK_MESSAGE( fabs(band - expected) <= (0.01*(std::max)(band, expected) + 1.0e-12),
                         "Band " << band << " disagrees with finite-difference expectation "
                         << expected << " at " << energy << " keV." );
  }//for( loop over energies )
}//BOOST_AUTO_TEST_CASE( PhysicalModelBandMatchesFiniteDifference )


// Validates the full-Jacobian relative-activity covariance (`m_rel_act_covariance`) invariants
// on an unconstrained problem: symmetry, per-isotope uncert == sqrt(diagonal), and agreement
// with the analytic unconstrained expectation norm_i^2 * m_nonlin_covariance[i][i] (the
// parameter-space covariance stays pristine; both matrices carry the same chi2/dof inflation).
BOOST_AUTO_TEST_CASE( RelActCovarianceInvariants )
{
  RelActCalcManual::RelEffInput calc_input = spec184_calc_input(
                              empirical_state_mod(RelActCalc::RelEffEqnForm::LnX, 4) );

  RelActCalcManual::RelEffSolution sol;
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcManual::solve_relative_efficiency( calc_input ) );
  BOOST_REQUIRE( sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );
  BOOST_REQUIRE( !sol.m_rel_act_covariance.empty() );
  BOOST_REQUIRE( sol.m_rel_act_covariance.size() == sol.m_rel_activities.size() );
  BOOST_REQUIRE( !sol.m_nonlin_covariance.empty() );
  BOOST_REQUIRE( sol.m_activity_norms.size() == sol.m_rel_activities.size() );

  for( size_t i = 0; i < sol.m_rel_act_covariance.size(); ++i )
  {
    for( size_t j = 0; j < sol.m_rel_act_covariance.size(); ++j )
    {
      BOOST_CHECK_MESSAGE( fabs(sol.m_rel_act_covariance[i][j] - sol.m_rel_act_covariance[j][i])
                             <= (1.0e-9*fabs(sol.m_rel_act_covariance[i][j]) + 1.0e-25),
                           "Rel-act covariance not symmetric at (" << i << "," << j << ")" );
    }

    const double diag = sol.m_rel_act_covariance[i][i];
    BOOST_CHECK( diag >= 0.0 );
    BOOST_CHECK_MESSAGE( fabs(sol.m_rel_activities[i].m_rel_activity_uncert - sqrt((std::max)(diag, 0.0)))
                           <= 1.0e-9*sqrt((std::max)(diag, 1.0e-30)),
                         "Activity uncert != sqrt(cov diagonal) for " << sol.m_rel_activities[i].m_isotope );

    // No constraints in this problem, so RelAct_i = norm_i * par_i exactly.
    const double expected = sol.m_activity_norms[i] * sol.m_activity_norms[i]
                            * sol.m_nonlin_covariance[i][i];
    BOOST_CHECK_MESSAGE( fabs(diag - expected) <= 1.0e-6*(std::max)(diag, expected),
                         "Jet-sweep covariance diagonal " << diag << " disagrees with analytic "
                         << expected << " for " << sol.m_rel_activities[i].m_isotope );
  }//for( loop over isotopes )
}//BOOST_AUTO_TEST_CASE( RelActCovarianceInvariants )


// Act-ratio constraint values and uncertainties (findings D(a)).  Ties U234 to U238 at the
// unconstrained best-fit ratio; the constrained pair's ratio uncertainty must be exactly zero,
// free-pair values must (nearly) reproduce the unconstrained solve, and the controlled nuclide
// inherits its controller's relative uncertainty exactly.
BOOST_AUTO_TEST_CASE( ActRatioConstraintUncerts )
{
  RelActCalcManual::RelEffInput calc_input = spec184_calc_input(
                              empirical_state_mod(RelActCalc::RelEffEqnForm::LnX, 4) );

  RelActCalcManual::RelEffSolution unc_sol;
  BOOST_REQUIRE_NO_THROW( unc_sol = RelActCalcManual::solve_relative_efficiency( calc_input ) );
  BOOST_REQUIRE( unc_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );

  const double unc_ratio_234_238 = unc_sol.activity_ratio( "U234", "U238" );
  const double unc_ratio_235_238 = unc_sol.activity_ratio( "U235", "U238" );
  double unc_uncert_235_238 = -1.0;
  BOOST_REQUIRE_NO_THROW( unc_uncert_235_238 = unc_sol.activity_ratio_uncert( "U235", "U238" ) );
  BOOST_REQUIRE( unc_uncert_235_238 > 0.0 );

  RelActCalcManual::RelEffInput con_input = calc_input;
  RelActCalcManual::ManualActRatioConstraint constraint;
  constraint.m_constrained_nuclide = "U234";
  constraint.m_controlling_nuclide = "U238";
  constraint.m_constrained_to_controlled_activity_ratio = unc_ratio_234_238;
  con_input.act_ratio_constraints.push_back( constraint );

  RelActCalcManual::RelEffSolution con_sol;
  BOOST_REQUIRE_NO_THROW( con_sol = RelActCalcManual::solve_relative_efficiency( con_input ) );
  BOOST_REQUIRE( con_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );

  // Constrained pair: exactly-fixed ratio, zero uncertainty.
  BOOST_CHECK_CLOSE_FRACTION( con_sol.activity_ratio("U234","U238"), unc_ratio_234_238, 1.0e-9 );
  double con_uncert_234_238 = -1.0;
  BOOST_REQUIRE_NO_THROW( con_uncert_234_238 = con_sol.activity_ratio_uncert( "U234", "U238" ) );
  BOOST_CHECK_SMALL( con_uncert_234_238, 1.0e-12 );

  // Constraining U234 at its natural best-fit value shouldnt move the free pair's value (much);
  //  its uncertainty legitimately shrinks (the constraint adds information), so only require the
  //  same rough scale.
  BOOST_CHECK_CLOSE_FRACTION( con_sol.activity_ratio("U235","U238"), unc_ratio_235_238, 0.02 );
  double con_uncert_235_238 = -1.0;
  BOOST_REQUIRE_NO_THROW( con_uncert_235_238 = con_sol.activity_ratio_uncert( "U235", "U238" ) );
  BOOST_CHECK( con_uncert_235_238 > 0.0 );
  BOOST_CHECK_MESSAGE( (con_uncert_235_238 > 0.4*unc_uncert_235_238)
                        && (con_uncert_235_238 < 2.0*unc_uncert_235_238),
                       "U235/U238 ratio uncert changed implausibly: " << unc_uncert_235_238
                       << " (unconstrained) vs " << con_uncert_235_238 << " (constrained)" );

  // The controlled nuclide inherits the controller's relative uncertainty exactly.
  const double u234_uncert = con_sol.relative_activity_uncertainty( "U234" );
  const double u238_uncert = con_sol.relative_activity_uncertainty( "U238" );
  BOOST_REQUIRE( (u234_uncert > 0.0) && (u238_uncert > 0.0) );
  BOOST_CHECK_CLOSE_FRACTION( u234_uncert / con_sol.relative_activity("U234"),
                              u238_uncert / con_sol.relative_activity("U238"), 1.0e-6 );
}//BOOST_AUTO_TEST_CASE( ActRatioConstraintUncerts )


// Range mass-fraction constraint uncertainties (findings D(b)).  U234 constrained to a window
// around its natural fraction: the solved fraction lands in the window; its activity uncertainty
// is positive (the old diagonal-only covariance transform is gone); the previously-unimplemented
// activity_ratio_uncert case (one nuclide mass-fraction constrained) now works; and the +-1 sigma
// mass-fraction variations bracket the nominal.
BOOST_AUTO_TEST_CASE( MassFracConstraintUncerts )
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );

  RelActCalcManual::RelEffInput calc_input = spec184_calc_input(
                              empirical_state_mod(RelActCalc::RelEffEqnForm::LnX, 4) );

  RelActCalcManual::RelEffSolution unc_sol;
  BOOST_REQUIRE_NO_THROW( unc_sol = RelActCalcManual::solve_relative_efficiency( calc_input ) );
  BOOST_REQUIRE( unc_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );
  const double u234_natural = unc_sol.mass_fraction( "U234" ); //spec184 is U-only, so == element fraction

  std::set<std::string> u_isotopes;
  for( const RelActCalcManual::GenericPeakInfo &peak : calc_input.peaks )
  {
    for( const RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas )
    {
      const SandiaDecay::Nuclide * const nuc = db->nuclide( line.m_isotope );
      if( nuc && (nuc->atomicNumber == 92) )
        u_isotopes.insert( line.m_isotope );
    }
  }
  BOOST_REQUIRE( u_isotopes.count("U234") && (u_isotopes.size() >= 2) );

  RelActCalcManual::RelEffInput con_input = calc_input;
  RelActCalcManual::MassFractionConstraint constraint;
  constraint.m_nuclide = "U234";
  constraint.m_mass_fraction_lower = 0.5 * u234_natural;
  constraint.m_mass_fraction_upper = (std::min)( 2.0 * u234_natural, 0.99 );
  for( const std::string &iso : u_isotopes )
    constraint.m_specific_activities[iso] = db->nuclide(iso)->activityPerGram();
  con_input.mass_fraction_constraints.push_back( constraint );

  RelActCalcManual::RelEffSolution con_sol;
  BOOST_REQUIRE_NO_THROW( con_sol = RelActCalcManual::solve_relative_efficiency( con_input ) );
  BOOST_REQUIRE( con_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );

  const double con_frac = con_sol.mass_fraction( "U234" );
  BOOST_CHECK( (con_frac >= (constraint.m_mass_fraction_lower - 1.0e-6))
               && (con_frac <= (constraint.m_mass_fraction_upper + 1.0e-6)) );

  double u234_uncert = -1.0;
  BOOST_REQUIRE_NO_THROW( u234_uncert = con_sol.relative_activity_uncertainty( "U234" ) );
  BOOST_CHECK_MESSAGE( u234_uncert > 0.0, "Range-constrained U234 has non-positive uncert" );

  // These two used to hit the "not implemented" -1.0 sentinel branches.
  double ratio_uncert = -1.0;
  BOOST_REQUIRE_NO_THROW( ratio_uncert = con_sol.activity_ratio_uncert( "U234", "U235" ) );
  BOOST_CHECK( ratio_uncert > 0.0 );
  BOOST_REQUIRE_NO_THROW( ratio_uncert = con_sol.activity_ratio_uncert( "U235", "U234" ) );
  BOOST_CHECK( ratio_uncert > 0.0 );

  const double nominal = con_sol.mass_fraction( "U234" );
  const double plus = con_sol.mass_fraction( "U234", 1.0 );
  const double minus = con_sol.mass_fraction( "U234", -1.0 );
  BOOST_CHECK( fabs(plus - minus) > 0.0 );
  BOOST_CHECK( ((std::min)(plus, minus) <= (nominal + 1.0e-9))
               && (nominal <= ((std::max)(plus, minus) + 1.0e-9)) );
}//BOOST_AUTO_TEST_CASE( MassFracConstraintUncerts )


// A FIXED mass-fraction nuclide in a MIXED-element problem used to get a hard-coded activity
// uncertainty of 0.0 (old TODO); the full-Jacobian covariance now propagates the element-total
// variation.  Synthetic problem: U235+U238 plus Cs137, flat true rel-eff (counts = act*yield),
// U235's fraction of the U element fixed at its true value.
BOOST_AUTO_TEST_CASE( FixedMassFracMixedElementUncert )
{
  set_data_dir();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );

  const double u235_act = 1.0e5, u238_act = 8.0e4, cs137_act = 1.2e5;

  const auto make_peak = []( const double energy, const double counts, const double yield,
                             const string &iso ) -> RelActCalcManual::GenericPeakInfo
  {
    RelActCalcManual::GenericPeakInfo peak;
    peak.m_energy = peak.m_mean = energy;
    peak.m_fwhm = 1.0;
    peak.m_counts = counts;
    peak.m_counts_uncert = sqrt( counts );
    peak.m_base_rel_eff_uncert = 0.01;
    peak.m_source_gammas.emplace_back( yield, iso );
    return peak;
  };

  RelActCalcManual::RelEffInput input;
  input.peaks.push_back( make_peak( 143.8, u235_act*0.1096, 0.1096, "U235" ) );
  input.peaks.push_back( make_peak( 185.7, u235_act*0.572, 0.572, "U235" ) );
  input.peaks.push_back( make_peak( 661.7, cs137_act*0.851, 0.851, "Cs137" ) );
  input.peaks.push_back( make_peak( 766.4, u238_act*0.00294, 0.00294, "U238" ) );
  input.peaks.push_back( make_peak( 1001.0, u238_act*0.00842, 0.00842, "U238" ) );
  input.eqn_form = RelActCalc::RelEffEqnForm::LnX;
  input.eqn_order = 1;

  const double u235_mass = u235_act / db->nuclide("U235")->activityPerGram();
  const double u238_mass = u238_act / db->nuclide("U238")->activityPerGram();
  const double u235_frac = u235_mass / (u235_mass + u238_mass);

  RelActCalcManual::MassFractionConstraint constraint;
  constraint.m_nuclide = "U235";
  constraint.m_mass_fraction_lower = u235_frac;
  constraint.m_mass_fraction_upper = u235_frac;
  constraint.m_specific_activities["U235"] = db->nuclide("U235")->activityPerGram();
  constraint.m_specific_activities["U238"] = db->nuclide("U238")->activityPerGram();
  input.mass_fraction_constraints.push_back( constraint );

  RelActCalcManual::RelEffSolution sol;
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcManual::solve_relative_efficiency( input ) );
  BOOST_REQUIRE( sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );

  // The fitted activities should reproduce truth (flat rel-eff, self-consistent counts).
  BOOST_CHECK_CLOSE_FRACTION( sol.activity_ratio("U235","Cs137"), u235_act/cs137_act, 0.02 );

  // Regression: fixed-mass-fraction nuclide in a mixed element must carry a positive activity
  //  uncertainty (its activity co-varies with the U-element total), not the old hard 0.0.
  double u235_uncert = -1.0;
  BOOST_REQUIRE_NO_THROW( u235_uncert = sol.relative_activity_uncertainty( "U235" ) );
  BOOST_CHECK_MESSAGE( u235_uncert > 0.0,
                       "Fixed U235 in a mixed element still has zero/negative uncertainty" );

  // Its fraction of ALL fitted mass still varies (the Cs mass does not co-vary proportionally).
  const double nominal = sol.mass_fraction( "U235" );
  double plus = nominal;
  BOOST_REQUIRE_NO_THROW( plus = sol.mass_fraction( "U235", 1.0 ) );
  BOOST_CHECK_MESSAGE( fabs(plus - nominal) > 1.0e-12,
                       "Fixed U235 all-mass fraction does not vary with +1 sigma" );
}//BOOST_AUTO_TEST_CASE( FixedMassFracMixedElementUncert )


// N1 cross-validation: for every empirical eqn form, `use_ceres_to_fit_eqn = true` and the
// default LLS method must give the same answers - point values tightly, uncertainties loosely
// (same reporting gauge, but different estimators: profiled-LLS vs held-coefficient covariance).
// The Ceres mode holds one coefficient to pin the scale gauge, and reports in the LLS convention
// (average measured rel. eff. == 1), so raw activities and the DOF are directly comparable.
BOOST_AUTO_TEST_CASE( CeresVsLlsEmpiricalAgreement )
{
  using RelActCalc::RelEffEqnForm;

  const pair<RelEffEqnForm,int> form_orders[] = {
    { RelEffEqnForm::LnX, 4 }, { RelEffEqnForm::LnY, 3 },
    { RelEffEqnForm::LnXLnY, 3 }, { RelEffEqnForm::FramEmpirical, 3 }
  };

  for( const pair<RelEffEqnForm,int> &form_order : form_orders )
  {
    const RelEffEqnForm form = form_order.first;
    const char * const form_name = RelActCalc::to_str( form );

    RelActCalcManual::RelEffInput calc_input = spec184_calc_input(
                              empirical_state_mod(form, form_order.second) );

    calc_input.use_ceres_to_fit_eqn = false;
    RelActCalcManual::RelEffSolution lls_sol;
    BOOST_REQUIRE_NO_THROW( lls_sol = RelActCalcManual::solve_relative_efficiency( calc_input ) );
    BOOST_REQUIRE_MESSAGE( lls_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success,
                           "LLS solve failed for " << form_name );

    calc_input.use_ceres_to_fit_eqn = true;
    RelActCalcManual::RelEffSolution ceres_sol;
    BOOST_REQUIRE_NO_THROW( ceres_sol = RelActCalcManual::solve_relative_efficiency( calc_input ) );
    BOOST_REQUIRE_MESSAGE( ceres_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success,
                           "Ceres solve failed for " << form_name );

    // Same effective-parameter bookkeeping in both modes -> identical DOF.
    BOOST_CHECK_MESSAGE( lls_sol.m_dof == ceres_sol.m_dof,
                         "DOF mismatch for " << form_name << ": " << lls_sol.m_dof
                         << " (LLS) vs " << ceres_sol.m_dof << " (Ceres)" );

    // Both modes report the SAME counts-space objective (m_chi2_fit_weights), but only the Ceres
    //  mode minimizes it over the curve coefficients directly, so it can never do worse:
    //   - LnX is fit un-transformed, and the LLS inner fit's weighted residual
    //     (C/S - f)/(sigma/S) is algebraically the counts residual (C - f*S)/sigma - i.e., an
    //     EXACT variable projection - so the two modes tie.
    //   - the exponential forms are fit in LOG space by the inner LLS, which is a different
    //     (approximate) criterion, so Ceres reaches a measurably lower chi2.
    //  Measured on this problem 20260816: LnX tie to ~1e-5; LnY 99.90 -> 98.59,
    //  LnXLnY 97.81 -> 95.60, FramEmpirical 99.75 -> 98.10.
    BOOST_CHECK_MESSAGE( ceres_sol.m_chi2_fit_weights
                           <= (lls_sol.m_chi2_fit_weights + 1.0e-3*lls_sol.m_chi2_fit_weights),
                         "Ceres mode should minimize the counts-space chi2 at least as well as the"
                         " LLS mode for " << form_name << ": " << lls_sol.m_chi2_fit_weights
                         << " (LLS) vs " << ceres_sol.m_chi2_fit_weights << " (Ceres)" );

    if( form == RelEffEqnForm::LnX ) //exact variable projection - must agree tightly
      BOOST_CHECK_CLOSE_FRACTION( lls_sol.m_chi2_fit_weights, ceres_sol.m_chi2_fit_weights, 1.0e-4 );

    // Scale-invariant point values.  LnX agrees to ~1e-4; for the exponential forms the two
    //  modes genuinely converge to slightly different minima (~1-2% in most quantities): the
    //  LLS mode's inner fit minimizes a log-space weighted rel-eff objective (an approximate
    //  variable projection), while the Ceres mode minimizes the true counts-space objective
    //  jointly.  The tolerances below document that intrinsic method difference.
    const double mass_frac_tol = (form == RelEffEqnForm::LnX) ? 1.0e-4 : 3.0e-3;
    BOOST_CHECK_MESSAGE( fabs(lls_sol.mass_fraction("U235") - ceres_sol.mass_fraction("U235")) < mass_frac_tol,
                         "U235 mass fraction mismatch for " << form_name << ": "
                         << lls_sol.mass_fraction("U235") << " vs " << ceres_sol.mass_fraction("U235") );
    BOOST_CHECK_CLOSE_FRACTION( lls_sol.activity_ratio("U235","U238"),
                                ceres_sol.activity_ratio("U235","U238"), 0.025 );

    // Gauge-dependent point values - comparable because both are in the LLS gauge.
    const double rel_act_tol = (form == RelEffEqnForm::LnX) ? 0.01 : 0.05;
    for( const RelActCalcManual::IsotopeRelativeActivity &lls_act : lls_sol.m_rel_activities )
    {
      const double ceres_act = ceres_sol.relative_activity( lls_act.m_isotope );
      BOOST_CHECK_MESSAGE( fabs(lls_act.m_rel_activity - ceres_act)
                             < rel_act_tol*(std::max)( fabs(lls_act.m_rel_activity), fabs(ceres_act) ),
                           "Rel act mismatch (" << form_name << ") for " << lls_act.m_isotope
                           << ": " << lls_act.m_rel_activity << " (LLS) vs " << ceres_act << " (Ceres)" );
    }

    const double curve_tol = (form == RelEffEqnForm::LnX) ? 0.01 : 0.04;
    for( const double energy : { 150.0, 400.0, 1001.0 } )
      BOOST_CHECK_CLOSE_FRACTION( lls_sol.rel_eff_eqn_value(energy),
                                  ceres_sol.rel_eff_eqn_value(energy), curve_tol );

    // Activity uncertainties: loose factor comparison (different estimators, same gauge).
    for( const RelActCalcManual::IsotopeRelativeActivity &lls_act : lls_sol.m_rel_activities )
    {
      const double lls_uncert = lls_act.m_rel_activity_uncert;
      const double ceres_uncert = ceres_sol.relative_activity_uncertainty( lls_act.m_isotope );
      BOOST_REQUIRE( lls_uncert > 0.0 );
      BOOST_CHECK_MESSAGE( (ceres_uncert > 0.5*lls_uncert) && (ceres_uncert < 2.0*lls_uncert),
                           "Rel act uncert mismatch (" << form_name << ") for " << lls_act.m_isotope
                           << ": " << lls_uncert << " (LLS) vs " << ceres_uncert << " (Ceres)" );
    }

    // Both modes report the band from the coefficient covariance CONDITIONAL on the fitted
    //  activities, so the band must not depend much on which method fit the coefficients.
    //  (Before 20260816 the Ceres mode used the marginal sub-block, which is dominated by the
    //  unidentifiable curve<->activity scale trade: bands were 2.4x to 528x the conditional
    //  ones on this problem.)
    for( const double energy : { 186.0, 400.0 } )
    {
      double lls_band = -1.0, ceres_band = -1.0;
      BOOST_REQUIRE_NO_THROW( lls_band = lls_sol.rel_eff_eqn_uncert(energy) );
      BOOST_REQUIRE_NO_THROW( ceres_band = ceres_sol.rel_eff_eqn_uncert(energy) );
      BOOST_REQUIRE( lls_band > 0.0 );
      BOOST_CHECK_MESSAGE( fabs(ceres_band - lls_band) < 0.25*(std::max)(lls_band, ceres_band),
                           "Band mismatch (" << form_name << ") at " << energy << " keV: "
                           << lls_band << " (LLS) vs " << ceres_band << " (Ceres)" );
    }
  }//for( loop over empirical eqn forms )
}//BOOST_AUTO_TEST_CASE( CeresVsLlsEmpiricalAgreement )


// Regression for the unweighted (m_base_rel_eff_uncert == -1) Ceres-empirical mode: the
// unweighted residual C_p/(k*S_p) - f/k scales as 1/k, so before the gauge pin the cost fell
// like 1/k^2 and the minimizer was driven toward infinite activities.  With one coefficient
// held, the solve must stay bounded and near the LLS-mode answer.
BOOST_AUTO_TEST_CASE( CeresEmpiricalUnweightedBounded )
{
  RelActCalcManual::RelEffInput calc_input = spec184_calc_input(
        empirical_state_mod( RelActCalc::RelEffEqnForm::LnY, 3, RelActManualGui::AddUncert::Unweighted ) );

  calc_input.use_ceres_to_fit_eqn = false;
  RelActCalcManual::RelEffSolution lls_sol;
  BOOST_REQUIRE_NO_THROW( lls_sol = RelActCalcManual::solve_relative_efficiency( calc_input ) );
  BOOST_REQUIRE( lls_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );

  calc_input.use_ceres_to_fit_eqn = true;
  RelActCalcManual::RelEffSolution ceres_sol;
  BOOST_REQUIRE_NO_THROW( ceres_sol = RelActCalcManual::solve_relative_efficiency( calc_input ) );
  BOOST_REQUIRE( ceres_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );

  for( const RelActCalcManual::IsotopeRelativeActivity &lls_act : lls_sol.m_rel_activities )
  {
    const double ceres_act = ceres_sol.relative_activity( lls_act.m_isotope );
    BOOST_CHECK_MESSAGE( !std::isnan(ceres_act) && !std::isinf(ceres_act) && (ceres_act > 0.0),
                         "Unweighted Ceres activity invalid for " << lls_act.m_isotope );
    BOOST_CHECK_MESSAGE( (ceres_act > 0.01*lls_act.m_rel_activity)
                          && (ceres_act < 100.0*lls_act.m_rel_activity),
                         "Unweighted Ceres activity for " << lls_act.m_isotope << " ("
                         << ceres_act << ") ran away from the LLS value (" << lls_act.m_rel_activity << ")" );
  }//for( loop over isotopes )
}//BOOST_AUTO_TEST_CASE( CeresEmpiricalUnweightedBounded )


// DOF bookkeeping (findings A5): empirical formula; Hoerl on/off differ by exactly 2 (the old
// bookkeeping charged one phantom parameter even with Hoerl off); a fixed mass-fraction
// constraint frees one parameter.
BOOST_AUTO_TEST_CASE( DofBookkeeping )
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );

  // Empirical LLS: dof = num_peaks - num_isotopes - eqn_order.
  RelActCalcManual::RelEffInput lnx_input = spec184_calc_input(
                              empirical_state_mod(RelActCalc::RelEffEqnForm::LnX, 4) );
  RelActCalcManual::RelEffSolution lnx_sol;
  BOOST_REQUIRE_NO_THROW( lnx_sol = RelActCalcManual::solve_relative_efficiency( lnx_input ) );
  BOOST_REQUIRE( lnx_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );
  const int num_peaks = static_cast<int>( lnx_input.peaks.size() );
  const int num_isos = static_cast<int>( lnx_sol.m_rel_activities.size() );
  BOOST_CHECK_EQUAL( lnx_sol.m_dof, num_peaks - num_isos - 4 );

  // Physical model: Hoerl on vs off must differ by exactly the two Hoerl parameters.
  RelActCalcManual::RelEffInput hoerl_on_input = spec184_calc_input(
                              []( RelActManualGui::GuiState &state ){ state.m_physModelUseHoerl = true; } );
  RelActCalcManual::RelEffInput hoerl_off_input = spec184_calc_input(
                              []( RelActManualGui::GuiState &state ){ state.m_physModelUseHoerl = false; } );
  RelActCalcManual::RelEffSolution hoerl_on_sol, hoerl_off_sol;
  BOOST_REQUIRE_NO_THROW( hoerl_on_sol = RelActCalcManual::solve_relative_efficiency( hoerl_on_input ) );
  BOOST_REQUIRE_NO_THROW( hoerl_off_sol = RelActCalcManual::solve_relative_efficiency( hoerl_off_input ) );
  BOOST_REQUIRE( hoerl_on_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );
  BOOST_REQUIRE( hoerl_off_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );
  BOOST_CHECK_EQUAL( hoerl_off_sol.m_dof, hoerl_on_sol.m_dof + 2 );

  // A fixed mass-fraction constraint holds one activity slot constant -> one more DOF.
  std::set<std::string> u_isotopes;
  for( const RelActCalcManual::GenericPeakInfo &peak : lnx_input.peaks )
    for( const RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas )
      if( db->nuclide(line.m_isotope) && (db->nuclide(line.m_isotope)->atomicNumber == 92) )
        u_isotopes.insert( line.m_isotope );

  RelActCalcManual::RelEffInput con_input = lnx_input;
  RelActCalcManual::MassFractionConstraint u234_fixed;
  u234_fixed.m_nuclide = "U234";
  u234_fixed.m_mass_fraction_lower = u234_fixed.m_mass_fraction_upper = lnx_sol.mass_fraction( "U234" );
  for( const std::string &iso : u_isotopes )
    u234_fixed.m_specific_activities[iso] = db->nuclide(iso)->activityPerGram();
  con_input.mass_fraction_constraints.push_back( u234_fixed );

  RelActCalcManual::RelEffSolution con_sol;
  BOOST_REQUIRE_NO_THROW( con_sol = RelActCalcManual::solve_relative_efficiency( con_input ) );
  BOOST_REQUIRE( con_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );
  BOOST_CHECK_EQUAL( con_sol.m_dof, lnx_sol.m_dof + 1 );
}//BOOST_AUTO_TEST_CASE( DofBookkeeping )


// Regression for the act-ratio/mass-fraction recursion hazard (findings A8): a nuclide whose
// act-ratio chain terminates on a mass-fraction constrained nuclide of its own element makes
// relative_activity() recurse without bound (decode -> relative_activity -> decode -> ...),
// crashing with a stack overflow rather than an exception.  The validator must reject it - and
// must NOT reject the safe variant where the chain terminates on another element.
BOOST_AUTO_TEST_CASE( A8ClosureConfigThrows )
{
  set_data_dir();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );

  const auto make_peak = []( const double energy, const double counts, const double yield,
                             const string &iso ) -> RelActCalcManual::GenericPeakInfo
  {
    RelActCalcManual::GenericPeakInfo peak;
    peak.m_energy = peak.m_mean = energy;
    peak.m_fwhm = 1.0;
    peak.m_counts = counts;
    peak.m_counts_uncert = sqrt( counts );
    peak.m_base_rel_eff_uncert = 0.01;
    peak.m_source_gammas.emplace_back( yield, iso );
    return peak;
  };

  RelActCalcManual::RelEffInput input;
  input.peaks.push_back( make_peak( 129.3, 1000.0, 6.31e-5, "Pu239" ) );
  input.peaks.push_back( make_peak( 160.3, 500.0, 4.02e-6, "Pu240" ) );
  input.peaks.push_back( make_peak( 185.7, 2000.0, 0.572, "U235" ) );
  input.peaks.push_back( make_peak( 413.7, 3000.0, 1.47e-5, "Pu239" ) );
  input.eqn_form = RelActCalc::RelEffEqnForm::LnX;
  input.eqn_order = 1;

  RelActCalcManual::MassFractionConstraint pu239_range;
  pu239_range.m_nuclide = "Pu239";
  pu239_range.m_mass_fraction_lower = 0.5;
  pu239_range.m_mass_fraction_upper = 0.7;
  pu239_range.m_specific_activities["Pu239"] = db->nuclide("Pu239")->activityPerGram();
  pu239_range.m_specific_activities["Pu240"] = db->nuclide("Pu240")->activityPerGram();
  input.mass_fraction_constraints.push_back( pu239_range );

  // Dangerous: Pu240's chain terminates on Pu239, which is mass-fraction constrained in a block
  //  whose decode needs Pu240's activity -> unbounded recursion.  Must throw (not crash).
  {
    RelActCalcManual::RelEffInput bad_input = input;
    RelActCalcManual::ManualActRatioConstraint constraint;
    constraint.m_constrained_nuclide = "Pu240";
    constraint.m_controlling_nuclide = "Pu239";
    constraint.m_constrained_to_controlled_activity_ratio = 0.3;
    bad_input.act_ratio_constraints.push_back( constraint );

    BOOST_CHECK_THROW( bad_input.check_nuclide_constraints(), std::exception );
  }

  // Safe: Pu240 controlled by U235 (chain terminates on an unconstrained, other-element nuclide).
  {
    RelActCalcManual::RelEffInput good_input = input;
    RelActCalcManual::ManualActRatioConstraint constraint;
    constraint.m_constrained_nuclide = "Pu240";
    constraint.m_controlling_nuclide = "U235";
    constraint.m_constrained_to_controlled_activity_ratio = 0.3;
    good_input.act_ratio_constraints.push_back( constraint );

    BOOST_CHECK_NO_THROW( good_input.check_nuclide_constraints() );
  }
}//BOOST_AUTO_TEST_CASE( A8ClosureConfigThrows )


// `add_nuclides_to_peaks` used to zip an energy-SORTED yield map positionally against the
// caller-ORDERED peak list, silently mis-assigning yields for unsorted input (findings A9).
BOOST_AUTO_TEST_CASE( UnsortedPeaksToAddNuclides )
{
  set_data_dir();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );

  const auto make_peak = []( const double energy ) -> RelActCalcManual::GenericPeakInfo
  {
    RelActCalcManual::GenericPeakInfo peak;
    peak.m_energy = peak.m_mean = energy;
    peak.m_fwhm = 2.0;
    peak.m_counts = 1000.0;
    peak.m_counts_uncert = sqrt( 1000.0 );
    peak.m_base_rel_eff_uncert = 0.01;
    return peak;
  };

  RelActCalcManual::SandiaDecayNuc u235, u238;
  u235.nuclide = db->nuclide( "U235" );
  u238.nuclide = db->nuclide( "U238" );
  u235.age = u238.age = 20.0 * PhysicalUnits::year;
  BOOST_REQUIRE( u235.nuclide && u238.nuclide );

  const vector<double> energies{ 143.8, 185.7, 766.4, 1001.0 };
  vector<RelActCalcManual::GenericPeakInfo> sorted_peaks, shuffled_peaks;
  for( const double energy : energies )
    sorted_peaks.push_back( make_peak(energy) );
  for( const size_t index : { 3, 0, 2, 1 } )
    shuffled_peaks.push_back( sorted_peaks[index] );

  vector<RelActCalcManual::GenericPeakInfo> sorted_result, shuffled_result;
  BOOST_REQUIRE_NO_THROW( sorted_result = RelActCalcManual::add_nuclides_to_peaks( sorted_peaks, {u235, u238}, 0.0, 1.5 ) );
  BOOST_REQUIRE_NO_THROW( shuffled_result = RelActCalcManual::add_nuclides_to_peaks( shuffled_peaks, {u235, u238}, 0.0, 1.5 ) );

  BOOST_REQUIRE( sorted_result.size() == shuffled_result.size() );

  for( const RelActCalcManual::GenericPeakInfo &sorted_peak : sorted_result )
  {
    bool found_match = false;
    for( const RelActCalcManual::GenericPeakInfo &shuffled_peak : shuffled_result )
    {
      if( shuffled_peak.m_energy != sorted_peak.m_energy )
        continue;

      found_match = true;
      BOOST_REQUIRE_MESSAGE( shuffled_peak.m_source_gammas.size() == sorted_peak.m_source_gammas.size(),
                             "Source gamma count differs at " << sorted_peak.m_energy << " keV" );
      for( size_t i = 0; i < sorted_peak.m_source_gammas.size(); ++i )
      {
        BOOST_CHECK_EQUAL( shuffled_peak.m_source_gammas[i].m_isotope, sorted_peak.m_source_gammas[i].m_isotope );
        BOOST_CHECK_MESSAGE( fabs(shuffled_peak.m_source_gammas[i].m_yield - sorted_peak.m_source_gammas[i].m_yield)
                               <= 1.0e-12*sorted_peak.m_source_gammas[i].m_yield,
                             "Yield mis-assigned at " << sorted_peak.m_energy << " keV for "
                             << sorted_peak.m_source_gammas[i].m_isotope );
      }
    }//for( find the same-energy peak in the shuffled result )

    BOOST_CHECK_MESSAGE( found_match, "No shuffled-result peak at " << sorted_peak.m_energy << " keV" );
  }//for( loop over sorted-result peaks )
}//BOOST_AUTO_TEST_CASE( UnsortedPeaksToAddNuclides )


// The public fit_rel_eff_eqn_lls overload forwards a caller-supplied isotope list into a
// lower_bound lookup: an unsorted list or a missing isotope used to silently attribute counts to
// the wrong nuclide (findings A10).  The overload now sorts (isotopes, rel_acts) together and
// the lookup checks equality.
BOOST_AUTO_TEST_CASE( LlsRejectsUnknownIsotope )
{
  set_data_dir();

  const auto make_peak = []( const double energy, const double counts, const double yield,
                             const string &iso ) -> RelActCalcManual::GenericPeakInfo
  {
    RelActCalcManual::GenericPeakInfo peak;
    peak.m_energy = peak.m_mean = energy;
    peak.m_fwhm = 1.0;
    peak.m_counts = counts;
    peak.m_counts_uncert = sqrt( counts );
    peak.m_base_rel_eff_uncert = 0.01;
    peak.m_source_gammas.emplace_back( yield, iso );
    return peak;
  };

  vector<RelActCalcManual::GenericPeakInfo> peaks;
  peaks.push_back( make_peak( 143.8, 1096.0, 0.1096, "U235" ) );
  peaks.push_back( make_peak( 185.7, 5720.0, 0.572, "U235" ) );
  peaks.push_back( make_peak( 766.4, 294.0, 0.00294, "U238" ) );
  peaks.push_back( make_peak( 1001.0, 842.0, 0.00842, "U238" ) );

  // A peak isotope missing from the list must throw, not silently take a neighbors activity.
  {
    vector<double> fit_pars;
    const vector<string> isotopes{ "U238" };
    const vector<double> rel_acts{ 1.0e5 };
    BOOST_CHECK_THROW( RelActCalcManual::fit_rel_eff_eqn_lls( RelActCalc::RelEffEqnForm::LnX, 1,
                              isotopes, rel_acts, peaks, fit_pars, nullptr ),
                       std::exception );
  }

  // An unsorted caller list must give identical coefficients to the sorted one.
  {
    vector<double> sorted_fit_pars, unsorted_fit_pars;
    const vector<string> sorted_isotopes{ "U235", "U238" }, unsorted_isotopes{ "U238", "U235" };
    const vector<double> sorted_rel_acts{ 1.0e5, 8.0e4 }, unsorted_rel_acts{ 8.0e4, 1.0e5 };

    BOOST_REQUIRE_NO_THROW( RelActCalcManual::fit_rel_eff_eqn_lls( RelActCalc::RelEffEqnForm::LnX, 1,
                              sorted_isotopes, sorted_rel_acts, peaks, sorted_fit_pars, nullptr ) );
    BOOST_REQUIRE_NO_THROW( RelActCalcManual::fit_rel_eff_eqn_lls( RelActCalc::RelEffEqnForm::LnX, 1,
                              unsorted_isotopes, unsorted_rel_acts, peaks, unsorted_fit_pars, nullptr ) );

    BOOST_REQUIRE( sorted_fit_pars.size() == unsorted_fit_pars.size() );
    for( size_t i = 0; i < sorted_fit_pars.size(); ++i )
      BOOST_CHECK_CLOSE_FRACTION( sorted_fit_pars[i], unsorted_fit_pars[i], 1.0e-12 );
  }
}//BOOST_AUTO_TEST_CASE( LlsRejectsUnknownIsotope )


// The `point_estimate_only` flag contract: same point estimates and chi2 as a full solve, but no
// covariance or derived uncertainties (used by the AN scan and profile-likelihood sub-solves).
BOOST_AUTO_TEST_CASE( PointEstimateOnlyFlag )
{
  RelActCalcManual::RelEffInput calc_input = spec184_calc_input(
                              empirical_state_mod(RelActCalc::RelEffEqnForm::LnX, 4) );

  RelActCalcManual::RelEffSolution full_sol;
  BOOST_REQUIRE_NO_THROW( full_sol = RelActCalcManual::solve_relative_efficiency( calc_input ) );
  BOOST_REQUIRE( full_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );

  calc_input.point_estimate_only = true;
  RelActCalcManual::RelEffSolution pe_sol;
  BOOST_REQUIRE_NO_THROW( pe_sol = RelActCalcManual::solve_relative_efficiency( calc_input ) );
  BOOST_REQUIRE( pe_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );

  BOOST_CHECK( pe_sol.m_nonlin_covariance.empty() );
  BOOST_CHECK( pe_sol.m_rel_act_covariance.empty() );
  BOOST_CHECK( pe_sol.m_rel_act_jacobian.empty() );
  for( const RelActCalcManual::IsotopeRelativeActivity &rel_act : pe_sol.m_rel_activities )
    BOOST_CHECK( rel_act.m_rel_activity_uncert < 0.0 );

  BOOST_CHECK_CLOSE_FRACTION( pe_sol.m_chi2, full_sol.m_chi2, 1.0e-9 );
  BOOST_CHECK_CLOSE_FRACTION( pe_sol.m_chi2_fit_weights, full_sol.m_chi2_fit_weights, 1.0e-9 );
  BOOST_CHECK_EQUAL( pe_sol.m_dof, full_sol.m_dof );

  for( const RelActCalcManual::IsotopeRelativeActivity &rel_act : full_sol.m_rel_activities )
    BOOST_CHECK_CLOSE_FRACTION( pe_sol.relative_activity(rel_act.m_isotope),
                                rel_act.m_rel_activity, 1.0e-9 );
}//BOOST_AUTO_TEST_CASE( PointEstimateOnlyFlag )


// Profile-likelihood on an easy, near-Gaussian case: the 1-sigma interval should roughly match
// the covariance-based mass_fraction(+-1) half-widths (spec184 is U-only, so element-relative
// and all-mass fractions coincide), and the 2-sigma interval must contain the 1-sigma one.
BOOST_AUTO_TEST_CASE( ProfileMatchesCovarianceEasyCase )
{
  RelActCalcManual::RelEffInput calc_input = spec184_calc_input(
                              empirical_state_mod(RelActCalc::RelEffEqnForm::LnX, 4) );

  RelActCalcManual::RelEffSolution sol;
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcManual::solve_relative_efficiency( calc_input ) );
  BOOST_REQUIRE( sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );

  RelActCalcManual::ProfileMassFractionOptions options;
  options.nuclide = "U235";

  RelActCalcManual::ProfileMassFractionResult profile;
  BOOST_REQUIRE_NO_THROW( profile = RelActCalcManual::profile_mass_fraction( calc_input, sol, options ) );

  BOOST_REQUIRE_EQUAL( profile.intervals.size(), size_t(2) );
  const RelActCalcManual::ProfileMassFractionInterval &one_sigma = profile.intervals[0];
  const RelActCalcManual::ProfileMassFractionInterval &two_sigma = profile.intervals[1];

  const double nominal = profile.nominal_mass_fraction;
  BOOST_CHECK_CLOSE_FRACTION( nominal, sol.mass_fraction("U235"), 1.0e-6 );

  BOOST_CHECK( !one_sigma.lower_at_bound && !one_sigma.upper_at_bound );
  BOOST_CHECK( one_sigma.lower_frac < nominal );
  BOOST_CHECK( one_sigma.upper_frac > nominal );

  // 2-sigma interval strictly contains the 1-sigma one.
  BOOST_CHECK( two_sigma.lower_frac < one_sigma.lower_frac );
  BOOST_CHECK( two_sigma.upper_frac > one_sigma.upper_frac );

  // Rough agreement with the covariance-based half-widths (both carry the same chi2/dof
  //  inflation): per side within a factor of two - the profile is genuinely asymmetric, which
  //  the covariance cannot express - and the total interval width within 40%.
  const double cov_plus = sol.mass_fraction( "U235", 1.0 ) - sol.mass_fraction( "U235" );
  const double cov_minus = sol.mass_fraction( "U235" ) - sol.mass_fraction( "U235", -1.0 );
  const double profile_plus = one_sigma.upper_frac - nominal;
  const double profile_minus = nominal - one_sigma.lower_frac;

  BOOST_CHECK_MESSAGE( (profile_plus > 0.5*cov_plus) && (profile_plus < 2.0*cov_plus),
                       "1-sigma upper half-width: profile " << profile_plus
                       << " vs covariance " << cov_plus );
  BOOST_CHECK_MESSAGE( (profile_minus > 0.5*cov_minus) && (profile_minus < 2.0*cov_minus),
                       "1-sigma lower half-width: profile " << profile_minus
                       << " vs covariance " << cov_minus );

  const double profile_width = profile_plus + profile_minus;
  const double cov_width = cov_plus + cov_minus;
  BOOST_CHECK_MESSAGE( fabs(profile_width - cov_width) < 0.4*(std::max)(profile_width, cov_width),
                       "1-sigma interval width: profile " << profile_width
                       << " vs covariance " << cov_width );
}//BOOST_AUTO_TEST_CASE( ProfileMatchesCovarianceEasyCase )


// Profile-likelihood end-to-end on the physical model (Hoerl on, per the file's embedded state):
// finite, bracketing, not at-bound intervals for U235 on the spec184 problem.
BOOST_AUTO_TEST_CASE( ProfilePhysicalModelTruth )
{
  RelActCalcManual::RelEffInput calc_input = spec184_calc_input( []( RelActManualGui::GuiState & ){} );
  BOOST_REQUIRE( calc_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel );

  RelActCalcManual::RelEffSolution sol;
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcManual::solve_relative_efficiency( calc_input ) );
  BOOST_REQUIRE( sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );

  RelActCalcManual::ProfileMassFractionOptions options;
  options.nuclide = "U235";

  RelActCalcManual::ProfileMassFractionResult profile;
  BOOST_REQUIRE_NO_THROW( profile = RelActCalcManual::profile_mass_fraction( calc_input, sol, options ) );

  BOOST_REQUIRE_EQUAL( profile.intervals.size(), size_t(2) );

  for( const RelActCalcManual::ProfileMassFractionInterval &interval : profile.intervals )
  {
    BOOST_CHECK( !std::isnan(interval.lower_frac) && !std::isnan(interval.upper_frac) );
    BOOST_CHECK_MESSAGE( interval.lower_frac < profile.nominal_mass_fraction,
                         "Lower " << interval.lower_frac << " not below nominal "
                         << profile.nominal_mass_fraction );
    BOOST_CHECK_MESSAGE( interval.upper_frac > profile.nominal_mass_fraction,
                         "Upper " << interval.upper_frac << " not above nominal "
                         << profile.nominal_mass_fraction );
    BOOST_CHECK( !interval.lower_at_bound && !interval.upper_at_bound );
  }//for( const auto &interval : profile.intervals )

  // Truth is 12.9543 wt% U235; the 2-sigma profile interval should be in the right neighborhood
  //  (this is a real spectrum with model error, so keep the check loose).
  BOOST_CHECK( profile.intervals[1].lower_frac < 0.15 );
  BOOST_CHECK( profile.intervals[1].upper_frac > 0.09 );
}//BOOST_AUTO_TEST_CASE( ProfilePhysicalModelTruth )


// A pre-existing range constraint on the profiled nuclide restricts the scan domain: end-points
// clip to the window with the at-bound flags set.
BOOST_AUTO_TEST_CASE( ProfileRespectsConstraintWindow )
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );

  RelActCalcManual::RelEffInput calc_input = spec184_calc_input(
                              empirical_state_mod(RelActCalc::RelEffEqnForm::LnX, 4) );

  RelActCalcManual::RelEffSolution unc_sol;
  BOOST_REQUIRE_NO_THROW( unc_sol = RelActCalcManual::solve_relative_efficiency( calc_input ) );
  BOOST_REQUIRE( unc_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );
  const double u235_nominal = unc_sol.mass_fraction( "U235" );

  std::set<std::string> u_isotopes;
  for( const RelActCalcManual::GenericPeakInfo &peak : calc_input.peaks )
    for( const RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas )
      if( db->nuclide(line.m_isotope) && (db->nuclide(line.m_isotope)->atomicNumber == 92) )
        u_isotopes.insert( line.m_isotope );

  // Constrain U235 to a window much narrower than its uncertainty, centered on the nominal.
  const double cov_sigma = 0.5*fabs( unc_sol.mass_fraction("U235",1.0) - unc_sol.mass_fraction("U235",-1.0) );
  BOOST_REQUIRE( cov_sigma > 0.0 );

  RelActCalcManual::RelEffInput con_input = calc_input;
  RelActCalcManual::MassFractionConstraint constraint;
  constraint.m_nuclide = "U235";
  constraint.m_mass_fraction_lower = u235_nominal - 0.2*cov_sigma;
  constraint.m_mass_fraction_upper = u235_nominal + 0.2*cov_sigma;
  for( const std::string &iso : u_isotopes )
    constraint.m_specific_activities[iso] = db->nuclide(iso)->activityPerGram();
  con_input.mass_fraction_constraints.push_back( constraint );

  RelActCalcManual::RelEffSolution con_sol;
  BOOST_REQUIRE_NO_THROW( con_sol = RelActCalcManual::solve_relative_efficiency( con_input ) );
  BOOST_REQUIRE( con_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );

  RelActCalcManual::ProfileMassFractionOptions options;
  options.nuclide = "U235";

  RelActCalcManual::ProfileMassFractionResult profile;
  BOOST_REQUIRE_NO_THROW( profile = RelActCalcManual::profile_mass_fraction( con_input, con_sol, options ) );
  BOOST_REQUIRE( !profile.intervals.empty() );

  for( const RelActCalcManual::ProfileMassFractionInterval &interval : profile.intervals )
  {
    BOOST_CHECK( interval.lower_frac >= (constraint.m_mass_fraction_lower - 1.0e-9) );
    BOOST_CHECK( interval.upper_frac <= (constraint.m_mass_fraction_upper + 1.0e-9) );
    BOOST_CHECK_MESSAGE( interval.lower_at_bound && interval.upper_at_bound,
                         "Expected both interval ends clipped to the constraint window (CL "
                         << interval.confidence_level << ")" );
  }//for( const auto &interval : profile.intervals )
}//BOOST_AUTO_TEST_CASE( ProfileRespectsConstraintWindow )


// The functor used to replace any initial activity estimate below 1.0 with 1.0, destroying the
// starting scale for genuinely small relative activities (e.g., CPS-scaled data); now only
// non-positive/non-finite estimates are replaced.  Synthetic flat-rel-eff problem where all true
// activities are ~1e-3.
BOOST_AUTO_TEST_CASE( NormFloorSmallActivities )
{
  set_data_dir();

  const double u235_act = 2.0e-3, u238_act = 1.0e-3;

  const auto make_peak = []( const double energy, const double counts, const double yield,
                             const string &iso ) -> RelActCalcManual::GenericPeakInfo
  {
    RelActCalcManual::GenericPeakInfo peak;
    peak.m_energy = peak.m_mean = energy;
    peak.m_fwhm = 1.0;
    peak.m_counts = counts;
    peak.m_counts_uncert = 0.02 * counts; //2% - as if CPS-scaled data
    peak.m_base_rel_eff_uncert = 0.01;
    peak.m_source_gammas.emplace_back( yield, iso );
    return peak;
  };

  RelActCalcManual::RelEffInput input;
  input.peaks.push_back( make_peak( 143.8, u235_act*0.1096, 0.1096, "U235" ) );
  input.peaks.push_back( make_peak( 185.7, u235_act*0.572, 0.572, "U235" ) );
  input.peaks.push_back( make_peak( 766.4, u238_act*0.00294, 0.00294, "U238" ) );
  input.peaks.push_back( make_peak( 1001.0, u238_act*0.00842, 0.00842, "U238" ) );
  input.eqn_form = RelActCalc::RelEffEqnForm::LnX;
  input.eqn_order = 1;

  RelActCalcManual::RelEffSolution sol;
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcManual::solve_relative_efficiency( input ) );
  BOOST_REQUIRE( sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );

  BOOST_CHECK_CLOSE_FRACTION( sol.activity_ratio("U235","U238"), u235_act/u238_act, 0.05 );

  // The old floor pushed a warning about replacing the small estimates; that must be gone.
  for( const std::string &warning : sol.m_warnings )
    BOOST_CHECK_MESSAGE( warning.find("will use 1.0 instead") == std::string::npos,
                         "Small (but valid) initial activity was floored: " << warning );
}//BOOST_AUTO_TEST_CASE( NormFloorSmallActivities )


// A single badly mis-fit peak (wrong continuum, bad skew, unmodeled interference) must not
// inflate every reported uncertainty.  The uncertainty scale uses the MEDIAN of the squared
// pulls, so it ignores a few outliers; the plain chi2/dof would be dragged up by the one bad
// peak.  The offending peak should also be named in the warnings so it can be re-fit.
BOOST_AUTO_TEST_CASE( RobustCovScaleIgnoresOutlierPeak )
{
  set_data_dir();

  const double u235_act = 1.0e5, u238_act = 8.0e4;

  const auto make_peak = []( const double energy, const double yield, const string &iso,
                             const double act ) -> RelActCalcManual::GenericPeakInfo
  {
    RelActCalcManual::GenericPeakInfo peak;
    peak.m_energy = peak.m_mean = energy;
    peak.m_fwhm = 1.0;
    peak.m_counts = act * yield;      //flat (==1) true rel. eff., so the model fits exactly
    peak.m_counts_uncert = sqrt( peak.m_counts );
    peak.m_base_rel_eff_uncert = 0.0; //statistical weights only, to keep the pulls clean
    peak.m_source_gammas.emplace_back( yield, iso );
    return peak;
  };

  RelActCalcManual::RelEffInput input;
  input.eqn_form = RelActCalc::RelEffEqnForm::LnX;
  input.eqn_order = 1;
  input.peaks.push_back( make_peak( 143.8, 0.1096, "U235", u235_act ) );
  input.peaks.push_back( make_peak( 163.4, 0.0508, "U235", u235_act ) );
  input.peaks.push_back( make_peak( 185.7, 0.5720, "U235", u235_act ) );
  input.peaks.push_back( make_peak( 202.1, 0.0108, "U235", u235_act ) );
  input.peaks.push_back( make_peak( 205.3, 0.0501, "U235", u235_act ) );
  input.peaks.push_back( make_peak( 766.4, 0.00294, "U238", u238_act ) );
  input.peaks.push_back( make_peak( 880.5, 0.00069, "U238", u238_act ) );
  input.peaks.push_back( make_peak( 1001.0, 0.00842, "U238", u238_act ) );
  input.peaks.push_back( make_peak( 1193.0, 0.00016, "U238", u238_act ) );
  input.peaks.push_back( make_peak( 1737.7, 0.00021, "U238", u238_act ) );

  RelActCalcManual::RelEffSolution clean_sol;
  BOOST_REQUIRE_NO_THROW( clean_sol = RelActCalcManual::solve_relative_efficiency( input ) );
  BOOST_REQUIRE( clean_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );
  BOOST_CHECK_MESSAGE( clean_sol.m_cov_scale < 1.05,
                       "Consistent data should not inflate uncertainties; scale was "
                       << clean_sol.m_cov_scale );

  // Now push ONE peak far off, as a mis-fit peak would be.  Deliberately NOT the strongest peak:
  //  that one carries most of its nuclide's weight, so the fitted activity simply follows it and
  //  the discrepancy is redistributed onto its siblings.  A mid-weight peak stays an outlier.
  RelActCalcManual::RelEffInput outlier_input = input;
  const double bad_energy = outlier_input.peaks[1].m_energy; //163.4 keV
  outlier_input.peaks[1].m_counts += 15.0 * outlier_input.peaks[1].m_counts_uncert;

  RelActCalcManual::RelEffSolution outlier_sol;
  BOOST_REQUIRE_NO_THROW( outlier_sol = RelActCalcManual::solve_relative_efficiency( outlier_input ) );
  BOOST_REQUIRE( outlier_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );

  // chi2/dof is badly degraded by the one peak...
  BOOST_REQUIRE( outlier_sol.m_dof > 0 );
  const double mean_based_scale = outlier_sol.m_chi2_fit_weights / outlier_sol.m_dof;
  BOOST_CHECK_MESSAGE( mean_based_scale > 3.0,
                       "Test setup: the outlier should dominate chi2/dof, got " << mean_based_scale );

  // ...but the robust scale, and hence everyone's uncertainty, must be left alone.
  BOOST_CHECK_MESSAGE( outlier_sol.m_cov_scale < 1.5,
                       "One mis-fit peak inflated the uncertainty scale to "
                       << outlier_sol.m_cov_scale << " (chi2/dof would give " << mean_based_scale << ")" );

  // And the peak should be called out, so the user can go re-fit it.
  bool warned_about_peak = false;
  for( const std::string &warning : outlier_sol.m_warnings )
    warned_about_peak |= (warning.find( SpecUtils::printCompact(bad_energy,5) ) != std::string::npos);
  BOOST_CHECK_MESSAGE( warned_about_peak,
                       "No warning identified the strongly deviating peak at " << bad_energy << " keV" );
}//BOOST_AUTO_TEST_CASE( RobustCovScaleIgnoresOutlierPeak )



// Profiling on top of an auto-estimated additional uncertainty must profile ONE likelihood: every
// trial has to keep the peak weights the nominal fit ended up with.  If a trial re-ran the
// estimate it would re-weight the peaks so its own scatter looked normal, flattening the profile
// into a "hit the bounds" non-answer (and costing several solves per trial).
BOOST_AUTO_TEST_CASE( ProfileWithAutoAddUncert )
{
  {
    RelActCalcManual::RelEffInput input = spec184_calc_input( []( RelActManualGui::GuiState & ){} );
    input.auto_estimate_add_uncert = true;
    for( RelActCalcManual::GenericPeakInfo &peak : input.peaks )
      peak.m_base_rel_eff_uncert = 0.0;

    const char * const form_name = "stat-multiple";

    RelActCalcManual::RelEffSolution sol;
    BOOST_REQUIRE_NO_THROW( sol = RelActCalcManual::solve_relative_efficiency( input ) );
    BOOST_REQUIRE( sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );
    BOOST_REQUIRE_MESSAGE( sol.m_auto_stat_uncert_multiple > 1.0,
                           "No automatic uncertainty estimate was made (" << form_name << ")" );

    RelActCalcManual::ProfileMassFractionOptions options;
    options.nuclide = "U235";

    RelActCalcManual::ProfileMassFractionResult profile;
    BOOST_REQUIRE_NO_THROW( profile = RelActCalcManual::profile_mass_fraction( input, sol, options ) );
    BOOST_REQUIRE( !profile.intervals.empty() );

    const RelActCalcManual::ProfileMassFractionInterval &one_sigma = profile.intervals.front();
    BOOST_CHECK_MESSAGE( !one_sigma.lower_at_bound && !one_sigma.upper_at_bound,
                         "Profile ran to the domain bounds (" << form_name << ") - the trial solves"
                         " are probably not on the same likelihood as the nominal one" );
    BOOST_CHECK( one_sigma.lower_frac < profile.nominal_mass_fraction );
    BOOST_CHECK( one_sigma.upper_frac > profile.nominal_mass_fraction );

    // A sane interval, not one spanning the whole physical range.
    BOOST_CHECK_MESSAGE( (one_sigma.upper_frac - one_sigma.lower_frac) < 0.5,
                         "Profile interval (" << form_name << ") spans "
                         << (one_sigma.upper_frac - one_sigma.lower_frac)
                         << " of the [0,1] mass-fraction range" );
  }
}//BOOST_AUTO_TEST_CASE( ProfileWithAutoAddUncert )


// The "Auto" additional uncertainty scales every peak's statistical uncertainty by ONE common
// factor.  The point of that form is that a uniform scaling factors straight out of the weighted
// least squares objective: the fit cannot move, only the uncertainties widen.  Validated against
// the reference materials 20260816 - it reproduces the statistics-only enrichment exactly, which
// is the most accurate setting we have, while covering truth (spec184 12.60 +- 1.71 vs truth
// 12.9543).  Contrast a fractional-of-peak-area term, which re-weights the peaks
// against each other and biases enrichment low.
BOOST_AUTO_TEST_CASE( AutoStatUncertMultipleDoesNotMoveFit )
{
  RelActCalcManual::RelEffInput stat_only = spec184_calc_input( []( RelActManualGui::GuiState & ){} );
  for( RelActCalcManual::GenericPeakInfo &peak : stat_only.peaks )
    peak.m_base_rel_eff_uncert = 0.0;

  RelActCalcManual::RelEffSolution stat_sol;
  BOOST_REQUIRE_NO_THROW( stat_sol = RelActCalcManual::solve_relative_efficiency( stat_only ) );
  BOOST_REQUIRE( stat_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );

  RelActCalcManual::RelEffInput auto_input = stat_only;
  auto_input.auto_estimate_add_uncert = true;

  RelActCalcManual::RelEffSolution auto_sol;
  BOOST_REQUIRE_NO_THROW( auto_sol = RelActCalcManual::solve_relative_efficiency( auto_input ) );
  BOOST_REQUIRE( auto_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );

  BOOST_REQUIRE_MESSAGE( auto_sol.m_auto_stat_uncert_multiple > 1.0,
                         "Expected the peaks to need widening; got multiple "
                         << auto_sol.m_auto_stat_uncert_multiple );

  // The fit must not have moved.
  BOOST_CHECK_CLOSE_FRACTION( auto_sol.mass_fraction("U235"), stat_sol.mass_fraction("U235"), 1.0e-4 );
  for( const RelActCalcManual::IsotopeRelativeActivity &act : stat_sol.m_rel_activities )
    BOOST_CHECK_CLOSE_FRACTION( auto_sol.relative_activity(act.m_isotope), act.m_rel_activity, 1.0e-4 );

  // Scaling every peak by k is the same thing the post-hoc robust covariance inflation was doing,
  //  so the estimated k should match sqrt(the statistics-only fit's variance inflation)...
  BOOST_CHECK_CLOSE_FRACTION( auto_sol.m_auto_stat_uncert_multiple,
                              sqrt(stat_sol.m_cov_scale), 0.05 );

  // ...which leaves nothing further to inflate, and lands chi2/dof (under the fit weights) near 1.
  BOOST_CHECK_MESSAGE( auto_sol.m_cov_scale < 1.5,
                       "Covariance was inflated a second time, scale = " << auto_sol.m_cov_scale );
  BOOST_REQUIRE( auto_sol.m_dof > 0 );
  const double fit_chi2_per_dof = auto_sol.m_chi2_fit_weights / auto_sol.m_dof;
  BOOST_CHECK_MESSAGE( (fit_chi2_per_dof > 0.2) && (fit_chi2_per_dof < 5.0),
                       "chi2/dof under the widened uncertainties should be near 1, got "
                       << fit_chi2_per_dof );

  // Reported uncertainties should agree with the statistics-only route (which inflates post-hoc).
  const double stat_uncert = 0.5*fabs( stat_sol.mass_fraction("U235",1.0) - stat_sol.mass_fraction("U235",-1.0) );
  const double auto_uncert = 0.5*fabs( auto_sol.mass_fraction("U235",1.0) - auto_sol.mass_fraction("U235",-1.0) );
  BOOST_CHECK_CLOSE_FRACTION( auto_uncert, stat_uncert, 0.05 );

  // The statistics-only chi2 must survive as a diagnostic - widening the uncertainties should not
  //  hide that the model does not describe the data within counting statistics.
  BOOST_CHECK_CLOSE_FRACTION( auto_sol.m_chi2, stat_sol.m_chi2, 0.05 );
}//BOOST_AUTO_TEST_CASE( AutoStatUncertMultipleDoesNotMoveFit )



// Data that already agrees with its own uncertainties needs no widening, so asking for the
// automatic estimate must be a no-op there.
BOOST_AUTO_TEST_CASE( AutoUncertOnConsistentData )
{
  set_data_dir();

  const double u235_act = 1.0e5, u238_act = 8.0e4;
  const auto make_peak = []( const double energy, const double yield, const string &iso,
                             const double act ) -> RelActCalcManual::GenericPeakInfo
  {
    RelActCalcManual::GenericPeakInfo peak;
    peak.m_energy = peak.m_mean = energy;
    peak.m_fwhm = 1.0;
    peak.m_counts = act * yield;   //flat (==1) true rel. eff., so the model fits exactly
    peak.m_counts_uncert = sqrt( peak.m_counts );
    peak.m_base_rel_eff_uncert = 0.0;
    peak.m_source_gammas.emplace_back( yield, iso );
    return peak;
  };

  RelActCalcManual::RelEffInput clean;
  clean.eqn_form = RelActCalc::RelEffEqnForm::LnX;
  clean.eqn_order = 1;
  clean.auto_estimate_add_uncert = true;
  clean.peaks.push_back( make_peak( 143.8, 0.1096, "U235", u235_act ) );
  clean.peaks.push_back( make_peak( 163.4, 0.0508, "U235", u235_act ) );
  clean.peaks.push_back( make_peak( 185.7, 0.5720, "U235", u235_act ) );
  clean.peaks.push_back( make_peak( 205.3, 0.0501, "U235", u235_act ) );
  clean.peaks.push_back( make_peak( 766.4, 0.00294, "U238", u238_act ) );
  clean.peaks.push_back( make_peak( 1001.0, 0.00842, "U238", u238_act ) );
  clean.peaks.push_back( make_peak( 1193.0, 0.00016, "U238", u238_act ) );
  clean.peaks.push_back( make_peak( 1737.7, 0.00021, "U238", u238_act ) );

  RelActCalcManual::RelEffSolution sol;
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcManual::solve_relative_efficiency( clean ) );
  BOOST_REQUIRE( sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );
  BOOST_CHECK_MESSAGE( sol.m_auto_stat_uncert_multiple < 1.05,
                       "Self-consistent data should need no widening, got multiple "
                       << sol.m_auto_stat_uncert_multiple );
}//BOOST_AUTO_TEST_CASE( AutoUncertOnConsistentData )


// "Stat. Only (no widening)" reports the raw statistical uncertainty - no widening for scatter -
// which on a real spectrum is far smaller than what the default reports.  The fit itself must be
// untouched either way; this only changes what is reported.
BOOST_AUTO_TEST_CASE( NoWideningOptionGivesRawStatisticalUncert )
{
  RelActCalcManual::RelEffInput widened = spec184_calc_input( []( RelActManualGui::GuiState & ){} );
  for( RelActCalcManual::GenericPeakInfo &peak : widened.peaks )
    peak.m_base_rel_eff_uncert = 0.0;

  RelActCalcManual::RelEffInput raw = widened;
  raw.widen_uncerts_for_scatter = false;

  RelActCalcManual::RelEffSolution widened_sol, raw_sol;
  BOOST_REQUIRE_NO_THROW( widened_sol = RelActCalcManual::solve_relative_efficiency( widened ) );
  BOOST_REQUIRE_NO_THROW( raw_sol = RelActCalcManual::solve_relative_efficiency( raw ) );
  BOOST_REQUIRE( widened_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );
  BOOST_REQUIRE( raw_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );

  BOOST_CHECK_CLOSE_FRACTION( raw_sol.mass_fraction("U235"), widened_sol.mass_fraction("U235"), 1.0e-9 );

  BOOST_CHECK_MESSAGE( raw_sol.m_cov_scale == 1.0,
                       "Uncertainties were widened despite being asked not to; scale "
                       << raw_sol.m_cov_scale );
  BOOST_REQUIRE_MESSAGE( widened_sol.m_cov_scale > 1.0,
                         "Test spectrum was expected to scatter beyond its uncertainties" );

  const double raw_uncert = 0.5*fabs( raw_sol.mass_fraction("U235",1.0) - raw_sol.mass_fraction("U235",-1.0) );
  const double widened_uncert = 0.5*fabs( widened_sol.mass_fraction("U235",1.0) - widened_sol.mass_fraction("U235",-1.0) );
  BOOST_REQUIRE( raw_uncert > 0.0 );
  BOOST_CHECK_CLOSE_FRACTION( widened_uncert / raw_uncert, sqrt(widened_sol.m_cov_scale), 0.05 );
}//BOOST_AUTO_TEST_CASE( NoWideningOptionGivesRawStatisticalUncert )


// Asking for a nuclide that is not in the solution must throw, not abort a debug build via an
// assert - callers (GUI chart title, LLM tool) legitimately ask for "U235" without knowing what
// was fit, inside a try/catch.
BOOST_AUTO_TEST_CASE( MassFractionUnknownNuclideThrows )
{
  RelActCalcManual::RelEffInput input = spec184_calc_input(
                              empirical_state_mod(RelActCalc::RelEffEqnForm::LnX, 4) );

  RelActCalcManual::RelEffSolution sol;
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcManual::solve_relative_efficiency( input ) );
  BOOST_REQUIRE( sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );

  BOOST_CHECK_THROW( sol.mass_fraction( "Pu239" ), std::exception );
  BOOST_CHECK_THROW( sol.mass_fraction( "Pu239", 1.0 ), std::exception );
  BOOST_CHECK_THROW( sol.mass_fraction( "Cs137", -1.0 ), std::exception );
}//BOOST_AUTO_TEST_CASE( MassFractionUnknownNuclideThrows )






