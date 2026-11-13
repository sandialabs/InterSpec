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

/* Element (fluorescence) x-ray sources in the "Isotopics from peaks" (RelActCalcManual) calculation:
 - with a Physical model whose self-attenuator contains the element, the x-rays get the same fluorescence
   source-depth correction RelActCalcAuto uses, so exact synthetic data built with it are fit exactly;
 - when the correction applies (`add_fluorescence_depth_models`);
 - the GUI's input preparation accepts element x-ray peaks;
 - an element x-ray peak's "use for isotopics" flag is opt-in, and survives an XML round trip.
 Needs the InterSpec "data" directory and test files (`-- --datadir=... --testfiledir=...`).
 */

#include "InterSpec_config.h"

#include <cmath>
#include <deque>
#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <filesystem>

#include "ceres/jet.h"

#include "rapidxml/rapidxml.hpp"

#define BOOST_TEST_MODULE RelActCalcManual_XrayFluorescence_suite
#include <boost/test/included/unit_test.hpp>

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
#include "InterSpec/RelActCalc_imp.hpp"
#include "InterSpec/RelEffShieldWidget.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using namespace boost::unit_test;

namespace
{
  string g_data_dir;
  string g_test_file_dir;

  void set_data_dir()
  {
    static bool s_have_set = false;
    if( s_have_set )
      return;
    s_have_set = true;

    const int argc = framework::master_test_suite().argc;
    char ** const argv = framework::master_test_suite().argv;

    string datadir;
    for( int i = 1; i < argc; ++i )
    {
      const string arg = argv[i];
      if( SpecUtils::istarts_with( arg, "--datadir=" ) )
        datadir = arg.substr( 10 );
      if( SpecUtils::istarts_with( arg, "--testfiledir=" ) )
        g_test_file_dir = arg.substr( 14 );
    }

    if( datadir.empty() )
    {
      for( const char *d : { "data", "../data", "../../data", "../../../data" } )
      {
        if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) )
        {
          datadir = d;
          break;
        }
      }
    }

    if( g_test_file_dir.empty() )
    {
      for( const char *d : { "test_data", "../test_data", "../../test_data" } )
      {
        if( SpecUtils::is_directory( SpecUtils::append_path(d, "manual_rel_eff") ) )
        {
          g_test_file_dir = d;
          break;
        }
      }
    }

    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( SpecUtils::append_path( datadir, "sandia.decay.xml" ) ),
                           "sandia.decay.xml not in '" << datadir << "'; pass -- --datadir=..." );
    BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );
    BOOST_REQUIRE( DecayDataBaseServer::database() );
    g_data_dir = datadir;

    if( !MaterialDB::initialized() )
      BOOST_REQUIRE_NO_THROW( MaterialDB::initialize() );
  }//set_data_dir()


  shared_ptr<const DetectorPeakResponse> generic_hpge()
  {
    auto drf = make_shared<DetectorPeakResponse>();
    drf->fromGadrasDirectory( SpecUtils::append_path( InterSpec::staticDataDirectory(),
                                                      "GenericGadrasDetectors/HPGe 40%" ) );
    BOOST_REQUIRE( drf->isValid() );
    return drf;
  }


  /** A U-metal (or `material_name`) Physical-model problem with U235, U238 and U fluorescence x-ray peaks
   (both K-alpha and K-beta, so the x-rays span more than 5 keV); peak counts are left at 1.0.
   */
  struct UProblem
  {
    RelActCalcManual::RelEffInput input;
    vector<RelActCalcManual::SandiaDecayNuc> sources;
  };

  UProblem make_u_problem( const string &material_name, const bool include_k_beta )
  {
    set_data_dir();
    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    const SandiaDecay::Nuclide * const u235 = db->nuclide( "U235" );
    const SandiaDecay::Nuclide * const u238 = db->nuclide( "U238" );
    const SandiaDecay::Element * const uranium = db->element( "U" );
    BOOST_REQUIRE( u235 && u238 && uranium );

    UProblem problem;
    for( const SandiaDecay::Nuclide *nuc : { u235, u238 } )
    {
      RelActCalcManual::SandiaDecayNuc src;
      src.nuclide = nuc;
      src.age = 20.0 * PhysicalUnits::year;
      problem.sources.push_back( src );
    }
    RelActCalcManual::SandiaDecayNuc el_src;
    el_src.element = uranium;
    problem.sources.push_back( el_src );

    vector<double> energies{ 143.767, 163.356, 185.715, 205.316, 766.361, 1001.026 };
    for( const SandiaDecay::EnergyIntensityPair &xray : uranium->xrays )
    {
      const bool k_alpha = (fabs(xray.energy - 94.654) < 0.01) || (fabs(xray.energy - 98.434) < 0.01);
      const bool k_beta = (fabs(xray.energy - 111.298) < 0.01);
      if( k_alpha || (include_k_beta && k_beta) )
        energies.push_back( xray.energy );
    }
    BOOST_REQUIRE_EQUAL( energies.size(), (include_k_beta ? 9 : 8) );

    vector<RelActCalcManual::GenericPeakInfo> peaks;
    for( const double energy : energies )
    {
      RelActCalcManual::GenericPeakInfo peak;
      peak.m_energy = peak.m_mean = energy;
      peak.m_fwhm = 0.7;
      peak.m_counts = 1.0;
      peak.m_counts_uncert = 1.0;
      peak.m_base_rel_eff_uncert = 0.0;
      peaks.push_back( peak );
    }

    vector<RelActCalcManual::PeakCsvInput::NucAndAge> isotopes;
    for( const RelActCalcManual::SandiaDecayNuc &src : problem.sources )
      isotopes.emplace_back( src.nuclide ? src.nuclide->symbol : src.element->symbol, src.age, false );

    const RelActCalcManual::PeakCsvInput::NucMatchResults matched
          = RelActCalcManual::PeakCsvInput::fill_in_nuclide_info( peaks,
                       RelActCalcManual::PeakCsvInput::NucDataSrc::SandiaDecay, {}, isotopes, 1.0f, {}, -1.0f );
    BOOST_REQUIRE( matched.unused_isotopes.empty() );
    BOOST_REQUIRE_EQUAL( matched.peaks_matched.size(), energies.size() );

    RelActCalcManual::RelEffInput &input = problem.input;
    input.peaks = matched.peaks_matched;
    input.eqn_form = RelActCalc::RelEffEqnForm::FramPhysicalModel;
    input.eqn_order = 0;
    input.use_ceres_to_fit_eqn = true;
    input.phys_model_use_hoerl = false;
    input.phys_model_detector = generic_hpge();

    auto self_atten = make_shared<RelActCalc::PhysicalModelShieldInput>();
    self_atten->material = MaterialDB::instance()->material( material_name );
    BOOST_REQUIRE( self_atten->material );
    self_atten->areal_density = 2.0 * PhysicalUnits::g_per_cm2;
    self_atten->fit_areal_density = true;
    self_atten->lower_fit_areal_density = 0.0;
    self_atten->upper_fit_areal_density = 50.0 * PhysicalUnits::g_per_cm2;
    input.phys_model_self_atten = self_atten;

    return problem;
  }//make_u_problem(...)


  /** Sets each peak's counts to the Physical model prediction - including the fluorescence correction of
   `input.phys_model_fluorescence` - for the given activities and self-attenuator areal density.
   */
  void set_exact_counts( RelActCalcManual::RelEffInput &input, const map<string,double> &activities,
                         const double areal_density )
  {
    RelActCalc::PhysModelShield<double> self_atten;
    self_atten.material = input.phys_model_self_atten->material;
    self_atten.areal_density = areal_density;

    for( RelActCalcManual::GenericPeakInfo &peak : input.peaks )
    {
      const double eff = RelActCalc::eval_physical_model_eqn( peak.m_energy, self_atten, {},
                                                 input.phys_model_detector.get(), std::nullopt, std::nullopt );
      double counts = 0.0;
      for( const RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas )
      {
        double corr = 1.0;
        for( const RelActCalcManual::FluorescenceDepthSource &fluor : input.phys_model_fluorescence )
        {
          if( fluor.m_element != line.m_isotope )
            continue;

          vector<double> weights;
          for( size_t i = 0; i < fluor.m_model.exciting_base_weight.size(); ++i )
            weights.push_back( fluor.m_model.exciting_base_weight[i] * activities.at( fluor.m_exciting_sources[i] ) );
          corr = RelActCalc::fluorescence_self_atten_correction( fluor.m_model.mu_at( peak.m_energy ),
                                                   fluor.m_model.exciting_mu, weights, areal_density );
        }//for( loop over fluorescing sources )

        counts += activities.at( line.m_isotope ) * line.m_yield * corr * eff;
      }//for( loop over lines )

      peak.m_counts = counts;
      peak.m_counts_uncert = sqrt( counts );
    }//for( loop over peaks )
  }//set_exact_counts(...)
}//namespace


BOOST_AUTO_TEST_CASE( fluorescence_depth_model_applicability )
{
  // U in a U-metal self-attenuator, with K-alpha and K-beta peaks: the model applies, excited by the nuclides.
  {
    UProblem problem = make_u_problem( "U", true );
    RelActCalcManual::add_fluorescence_depth_models( problem.input, problem.sources, 0.0 );
    BOOST_REQUIRE_EQUAL( problem.input.phys_model_fluorescence.size(), 1 );
    const RelActCalcManual::FluorescenceDepthSource &fluor = problem.input.phys_model_fluorescence[0];
    BOOST_CHECK_EQUAL( fluor.m_element, "U" );
    BOOST_REQUIRE( !fluor.m_exciting_sources.empty() );
    for( const string &src : fluor.m_exciting_sources )
      BOOST_CHECK_MESSAGE( (src == "U235") || (src == "U238"), "Unexpected exciting source " << src );
  }

  // Only the K-alpha lines (< 5 keV apart): the correction is ~constant over them, so it is not used.
  {
    UProblem problem = make_u_problem( "U", false );
    RelActCalcManual::add_fluorescence_depth_models( problem.input, problem.sources, 0.0 );
    BOOST_CHECK( problem.input.phys_model_fluorescence.empty() );
  }

  // A self-attenuator without uranium.
  {
    UProblem problem = make_u_problem( "Fe", true );
    RelActCalcManual::add_fluorescence_depth_models( problem.input, problem.sources, 0.0 );
    BOOST_CHECK( problem.input.phys_model_fluorescence.empty() );
  }

  // Not the Physical model.
  {
    UProblem problem = make_u_problem( "U", true );
    problem.input.eqn_form = RelActCalc::RelEffEqnForm::LnX;
    RelActCalcManual::add_fluorescence_depth_models( problem.input, problem.sources, 0.0 );
    BOOST_CHECK( problem.input.phys_model_fluorescence.empty() );
  }
}//fluorescence_depth_model_applicability


BOOST_AUTO_TEST_CASE( exact_data_with_fluorescence_is_fit_exactly )
{
  UProblem problem = make_u_problem( "U", true );
  RelActCalcManual::RelEffInput &input = problem.input;
  RelActCalcManual::add_fluorescence_depth_models( input, problem.sources, 0.0 );
  BOOST_REQUIRE_EQUAL( input.phys_model_fluorescence.size(), 1 );

  // Check the model's set-up against attenuation coefficients computed here, since the exact data below are
  //  generated with the model itself.
  {
    const RelActCalcManual::FluorescenceDepthSource &fluor = input.phys_model_fluorescence[0];
    BOOST_REQUIRE( input.phys_model_self_atten && input.phys_model_self_atten->material );
    const Material * const material = input.phys_model_self_atten->material.get();
    BOOST_CHECK( fluor.m_material == material );

    // Mass attenuation coefficient (per areal density, the units the correction takes), at the energy.
    const auto mu_of = [material]( const double energy ) -> double {
      return GammaInteractionCalc::transmition_length_coefficient( material, energy ) / material->density;
    };
    //  (`mu_at` gives the value at the element's nearest x-ray line, so allow for rounding of the energy.)
    BOOST_CHECK_CLOSE( fluor.m_model.mu_at( 98.434 ), mu_of( 98.434 ), 0.01 );
    BOOST_CHECK_CLOSE( fluor.m_model.mu_at( 111.298 ), mu_of( 111.298 ), 0.01 );

    // Only photons above the U K edge (115.6 keV) excite K x-rays, and U235's 185.7 keV line dominates its
    //  excitation (lines are merged into ~3 % wide bins of attenuation coefficient).
    const double max_mu = mu_of( 115.7 );
    double u235_max_weight = 0.0, u235_max_weight_mu = 0.0;
    for( size_t i = 0; i < fluor.m_model.exciting_mu.size(); ++i )
    {
      BOOST_CHECK( (fluor.m_model.exciting_mu[i] > 0.0) && (fluor.m_model.exciting_mu[i] <= max_mu) );
      BOOST_CHECK( fluor.m_model.exciting_base_weight[i] > 0.0 );
      if( (fluor.m_exciting_sources[i] == "U235") && (fluor.m_model.exciting_base_weight[i] > u235_max_weight) )
      {
        u235_max_weight = fluor.m_model.exciting_base_weight[i];
        u235_max_weight_mu = fluor.m_model.exciting_mu[i];
      }
    }//for( loop over exciting lines )
    BOOST_CHECK_CLOSE( u235_max_weight_mu, mu_of( 185.715 ), 3.0 );
  }

  const double true_ad = 10.0 * PhysicalUnits::g_per_cm2;
  const map<string,double> true_acts{ {"U235", 2.0E8}, {"U238", 1.0E9}, {"U", 1.0E8} };
  set_exact_counts( input, true_acts, true_ad );

  for( const RelActCalcManual::GenericPeakInfo &peak : input.peaks )
  {
    BOOST_TEST_MESSAGE( "Peak " << peak.m_energy << " keV: " << peak.m_counts << " counts, "
                        << peak.m_source_gammas.size() << " source line(s)" );
    BOOST_REQUIRE( peak.m_counts > 100.0 );
  }

  // The size of the effect being tested: K-beta1 relative to K-alpha1, versus a uniform source.
  {
    const RelActCalcManual::FluorescenceDepthSource &fluor = input.phys_model_fluorescence[0];
    vector<double> weights;
    for( size_t i = 0; i < fluor.m_model.exciting_base_weight.size(); ++i )
      weights.push_back( fluor.m_model.exciting_base_weight[i] * true_acts.at( fluor.m_exciting_sources[i] ) );
    const double ka1 = RelActCalc::fluorescence_self_atten_correction( fluor.m_model.mu_at( 98.434 ),
                                                     fluor.m_model.exciting_mu, weights, true_ad );
    const double kb1 = RelActCalc::fluorescence_self_atten_correction( fluor.m_model.mu_at( 111.298 ),
                                                     fluor.m_model.exciting_mu, weights, true_ad );
    BOOST_TEST_MESSAGE( "Fluorescence correction: K-alpha1 " << ka1 << ", K-beta1 " << kb1 );
    BOOST_CHECK( (kb1 / ka1) > 1.005 );
    BOOST_CHECK( (kb1 / ka1) < 1.15 );
  }

  RelActCalcManual::RelEffSolution sol = RelActCalcManual::solve_relative_efficiency( input );
  BOOST_REQUIRE_MESSAGE( sol.m_status == RelActCalcManual::ManualSolutionStatus::Success,
                         "Solve failed: " << sol.m_error_message );
  BOOST_REQUIRE( sol.m_phys_model_self_atten_shield );

  const double fit_ad = sol.m_phys_model_self_atten_shield->m_areal_density;
  BOOST_TEST_MESSAGE( "With the fluorescence model: chi2=" << sol.m_chi2 << ", AD="
                      << fit_ad / PhysicalUnits::g_per_cm2 << " g/cm2" );
  BOOST_CHECK_SMALL( sol.m_chi2, 1.0E-3 );
  BOOST_CHECK_CLOSE( fit_ad, true_ad, 0.1 );
  BOOST_CHECK_CLOSE( sol.relative_activity("U235") / sol.relative_activity("U238"),
                     true_acts.at("U235") / true_acts.at("U238"), 0.01 );

  // Post-fit accessors carry the same correction as the fit.
  BOOST_CHECK_EQUAL( sol.fluorescence_correction( "U235", 185.715 ), 1.0 );
  for( const RelActCalcManual::GenericPeakInfo &peak : sol.peaks_with_effective_yields() )
  {
    double src_counts = 0.0;
    for( const RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas )
      src_counts += sol.relative_activity( line.m_isotope ) * line.m_yield;
    BOOST_CHECK_CLOSE( src_counts * sol.rel_eff_eqn_value( peak.m_energy ), peak.m_counts, 0.05 );
  }

  // The text and HTML outputs handle the element source.
  {
    stringstream summary;
    BOOST_CHECK_NO_THROW( sol.print_summary( summary ) );
    BOOST_CHECK( !summary.str().empty() );

    // Without a WApplication, the HTML report reads its template from "InterSpec_resources/" relative to the
    //  working directory - the repository's base directory, which holds the "data" directory.
    namespace fs = std::filesystem;
    const fs::path orig_cwd = fs::current_path();
    const fs::path base_dir = fs::canonical( g_data_dir ).parent_path();
    if( fs::is_directory( base_dir / "InterSpec_resources" ) )
    {
      stringstream html;
      fs::current_path( base_dir );
      BOOST_CHECK_NO_THROW( sol.print_html_report( html, "exact U", nullptr, {}, nullptr, 0.0, {} ) );
      fs::current_path( orig_cwd );
      BOOST_CHECK( html.str().find( "U235" ) != string::npos );
    }else
    {
      BOOST_TEST_MESSAGE( "Skipping HTML report: no InterSpec_resources in " << base_dir );
    }
  }

  // The model is tied to the material it was made for.
  {
    RelActCalcManual::RelEffInput other_material_input = input;
    auto self_atten = make_shared<RelActCalc::PhysicalModelShieldInput>( *input.phys_model_self_atten );
    self_atten->material = make_shared<Material>( *self_atten->material );
    other_material_input.phys_model_self_atten = self_atten;
    const RelActCalcManual::RelEffSolution other_sol = RelActCalcManual::solve_relative_efficiency( other_material_input );
    BOOST_CHECK( other_sol.m_status != RelActCalcManual::ManualSolutionStatus::Success );
  }

  // The same data, fit without the correction, cannot be described by the uniform-source model.
  input.phys_model_fluorescence.clear();
  RelActCalcManual::RelEffSolution uniform_sol = RelActCalcManual::solve_relative_efficiency( input );
  BOOST_REQUIRE( uniform_sol.m_status == RelActCalcManual::ManualSolutionStatus::Success );
  BOOST_TEST_MESSAGE( "Without the fluorescence model: chi2=" << uniform_sol.m_chi2 << ", AD="
             << uniform_sol.m_phys_model_self_atten_shield->m_areal_density / PhysicalUnits::g_per_cm2 << " g/cm2" );
  BOOST_CHECK( uniform_sol.m_chi2 > 10.0 );
}//exact_data_with_fluorescence_is_fit_exactly


BOOST_AUTO_TEST_CASE( gui_input_accepts_element_xray_peaks )
{
  set_data_dir();
  BOOST_REQUIRE( !g_test_file_dir.empty() );

  const string n42 = SpecUtils::append_path( g_test_file_dir, "manual_rel_eff/spec184_235U_12.9543.n42" );
  auto infile = make_shared<SpecMeas>();
  BOOST_REQUIRE( infile->load_file( n42, SpecUtils::ParserType::Auto ) );

  const shared_ptr<deque<shared_ptr<const PeakDef>>> file_peaks = infile->peaks( {1} );
  BOOST_REQUIRE( file_peaks && !file_peaks->empty() );

  auto state = make_shared<RelActManualGui::GuiState>();
  BOOST_REQUIRE( infile->relActManualGuiState() );
  BOOST_REQUIRE_NO_THROW( state->deSerialize( infile->relActManualGuiState()->first_node() ) );
  BOOST_REQUIRE( state->m_relEffEqnFormIndex == RelActCalc::RelEffEqnForm::FramPhysicalModel );

  const SandiaDecay::Element * const uranium = DecayDataBaseServer::database()->element( "U" );

  // The file has U, Th, and Pb x-ray peaks, none of which are marked for use.
  size_t num_xray_peaks = 0;
  for( const shared_ptr<const PeakDef> &p : *file_peaks )
  {
    if( p->xrayElement() )
    {
      num_xray_peaks += 1;
      BOOST_CHECK_MESSAGE( !p->useForManualRelEff(), "Use-flag set on " << p->xrayElement()->symbol
                                                     << " x-ray peak at " << p->mean() << " keV" );
    }
  }//for( loop over peaks )
  BOOST_CHECK( num_xray_peaks >= 3 );

  // Opt the U K-alpha peaks in, and add a U K-beta1 peak so the x-rays span > 5 keV.
  RelActManualGui::RelActCalcRawInput raw;
  raw.state = state;
  raw.fore_spec = infile->measurement( 1, infile->detector_names().at(0) );
  raw.detector = infile->detector();
  double ka1_amp = 0.0;
  for( const shared_ptr<const PeakDef> &p : *file_peaks )
  {
    shared_ptr<PeakDef> peak = make_shared<PeakDef>( *p );
    if( peak->xrayElement() == uranium )
    {
      peak->useForManualRelEff( true );
      if( fabs(peak->xrayEnergy() - 98.434) < 0.01 )
        ka1_amp = peak->amplitude();
    }

    if( peak->useForManualRelEff() && (peak->parentNuclide() || peak->reaction() || peak->xrayElement()) )
      raw.peaks.push_back( peak );
  }//for( loop over peaks )
  BOOST_REQUIRE( ka1_amp > 0.0 );

  auto kb1 = make_shared<PeakDef>( 111.298, 0.3, 0.42 * ka1_amp );
  kb1->setAmplitudeUncert( sqrt( kb1->amplitude() ) );
  kb1->setXray( uranium, 111.298f );
  BOOST_CHECK( !kb1->useForManualRelEff() );  // x-rays are opt-in
  kb1->useForManualRelEff( true );
  BOOST_CHECK( kb1->useForManualRelEff() );
  raw.peaks.push_back( kb1 );

  RelActCalcManual::RelEffInput calc_input;
  BOOST_REQUIRE_NO_THROW( RelActManualGui::prepare_calc_input( raw, calc_input ) );

  size_t num_u_xray_lines = 0;
  for( const RelActCalcManual::GenericPeakInfo &peak : calc_input.peaks )
    for( const RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas )
      num_u_xray_lines += (line.m_isotope == "U");
  BOOST_CHECK_EQUAL( num_u_xray_lines, 3 );

  // UO2 self-attenuator contains U, so the depth correction is set up.
  BOOST_REQUIRE_EQUAL( calc_input.phys_model_fluorescence.size(), 1 );
  BOOST_CHECK_EQUAL( calc_input.phys_model_fluorescence[0].m_element, "U" );

  RelActCalcManual::RelEffSolution sol = RelActCalcManual::solve_relative_efficiency( calc_input );
  BOOST_REQUIRE_MESSAGE( sol.m_status == RelActCalcManual::ManualSolutionStatus::Success,
                         "Solve failed: " << sol.m_error_message );
  BOOST_TEST_MESSAGE( "spec184 with U x-rays: chi2/dof=" << sol.m_chi2 << "/" << sol.m_dof
                      << ", U235 enrichment " << sol.mass_fraction( "U235" ) );
  BOOST_CHECK( sol.relative_activity( "U" ) > 0.0 );

  // K-beta1 escapes the material more easily than K-alpha1, so is raised relative to a uniform source.
  const double kb1_ka1 = sol.fluorescence_correction( "U", 111.298 ) / sol.fluorescence_correction( "U", 98.434 );
  BOOST_TEST_MESSAGE( "spec184 K-beta1/K-alpha1 fluorescence correction: " << kb1_ka1 );
  BOOST_CHECK( (kb1_ka1 > 1.0) && (kb1_ka1 < 1.1) );
}//gui_input_accepts_element_xray_peaks


BOOST_AUTO_TEST_CASE( xray_peak_use_flag_is_opt_in_and_serialized )
{
  set_data_dir();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Element * const uranium = db->element( "U" );
  const SandiaDecay::Element * const lead = db->element( "Pb" );

  PeakDef peak( 98.434, 0.3, 1000.0 );
  BOOST_REQUIRE( peak.useForManualRelEffUserPreference() );  // the PeakDef default
  peak.setXray( uranium, 98.434f );
  BOOST_CHECK( !peak.useForManualRelEff() );

  peak.useForManualRelEff( true );
  peak.setXray( uranium, 94.654f );  // same element: keeps the users choice
  BOOST_CHECK( peak.useForManualRelEff() );
  peak.setXray( lead, 74.969f );     // different element: back to off
  BOOST_CHECK( !peak.useForManualRelEff() );
  peak.setXray( uranium, 98.434f );
  peak.useForManualRelEff( true );

  const auto round_trip = []( const PeakDef &input ) -> PeakDef {
    rapidxml::xml_document<char> doc;
    rapidxml::xml_node<char> * const base = doc.allocate_node( rapidxml::node_element, "Base" );
    doc.append_node( base );
    map<shared_ptr<PeakContinuum>,int> continuums;
    rapidxml::xml_node<char> * const peak_node = input.toXml( base, base, continuums );
    BOOST_REQUIRE( peak_node && (continuums.size() == 1) );
    map<int,shared_ptr<PeakContinuum>> id_to_continuum;
    for( const auto &cont_id : continuums )
      id_to_continuum[cont_id.second] = cont_id.first;
    PeakDef output;
    output.fromXml( peak_node, id_to_continuum );
    return output;
  };//round_trip

  BOOST_CHECK( round_trip( peak ).useForManualRelEff() );
  peak.useForManualRelEff( false );
  BOOST_CHECK( !round_trip( peak ).useForManualRelEff() );

  // Nuclide peaks are unaffected.
  PeakDef gamma( 185.7, 0.4, 1000.0 );
  gamma.setNuclearTransition( db->nuclide( "U235" ), nullptr, -1, PeakDef::NormalGamma );
  gamma.useForManualRelEff( true );
  BOOST_CHECK( round_trip( gamma ).useForManualRelEffUserPreference() );
}//xray_peak_use_flag_is_opt_in_and_serialized
