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
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include "InterSpec_config.h"

#include <cmath>
#include <memory>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>

// <boost/test/included/unit_test.hpp> includes <windows.h>; pulling winsock2.h in first keeps any
//  header that later needs it consistent, and protects against a future include pulling
//  <boost/asio> (MSVC C1189 "WinSock.h has already been included" otherwise).  NOTE this TU is
//  currently "light" (nothing in its include chain pulls boost/asio, so NOMINMAX is NOT defined):
//  any std::min/std::max added here must use the parenthesized (std::max)( a, b ) form - see
//  CLAUDE.md.  Same include approach as test_DetectorPeakResponse.cpp.
#ifdef _WIN32
  #define WIN32_LEAN_AND_MEAN
  #include <winsock2.h>
  #include <windows.h>
#endif

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/BatchRelActAuto.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

#define BOOST_TEST_MODULE RelActCalcAuto_MultiCurve_suite
#include <boost/test/included/unit_test.hpp>

using namespace std;
using namespace boost::unit_test;

static_assert( USE_REL_ACT_TOOL, "Compile-time support for the Rel. Act. tool is required for this test." );


namespace
{
  string g_test_file_dir, g_data_dir;

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

    SpecUtils::ireplace_all( datadir, "%20", " " );
    SpecUtils::ireplace_all( g_test_file_dir, "%20", " " );

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
      for( const char *d : { "test_data", "../test_data", "../../test_data",
                             "../../target/testing/test_data" } )
      {
        if( SpecUtils::is_directory( SpecUtils::append_path(d, "multicurve") ) )
        {
          g_test_file_dir = d;
          break;
        }
      }
    }

    g_data_dir = datadir;
    const string sandia_decay_file = SpecUtils::append_path( datadir, "sandia.decay.xml" );
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( sandia_decay_file ),
                           "sandia.decay.xml not at '" << sandia_decay_file << "'" );

    BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );

    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    BOOST_REQUIRE_MESSAGE( db, "Error initializing SandiaDecayDataBase" );
  }//set_data_dir()


  shared_ptr<const SpecUtils::Measurement> load_meas( const string &leaf_name )
  {
    const string path = SpecUtils::append_path(
                          SpecUtils::append_path( g_test_file_dir, "multicurve" ), leaf_name );
    auto meas = make_shared<SpecMeas>();
    BOOST_REQUIRE_MESSAGE( meas->load_file( path, SpecUtils::ParserType::Auto ),
                           "Failed to load '" << path << "'" );
    BOOST_REQUIRE( !meas->measurements().empty() );
    return meas->measurements().front();
  }//load_meas(...)


  shared_ptr<DetectorPeakResponse> load_detective_x()
  {
    const string tsv_path = SpecUtils::append_path( g_data_dir, "common_drfs.tsv" );
    ifstream tsv( tsv_path.c_str(), ios::binary );
    BOOST_REQUIRE_MESSAGE( tsv.good(), "Could not open '" << tsv_path << "'" );

    vector<string> credits;
    vector<shared_ptr<DetectorPeakResponse>> drfs;
    DetectorPeakResponse::parseMultipleRelEffDrfCsv( tsv, credits, drfs );

    for( const shared_ptr<DetectorPeakResponse> &drf : drfs )
    {
      if( drf && SpecUtils::icontains( drf->name(), "Detective-X_LANL_100cm" ) )
        return drf;
    }
    BOOST_REQUIRE_MESSAGE( false, "No Detective-X_LANL_100cm DRF in common_drfs.tsv" );
    return nullptr;
  }//load_detective_x()


  RelActCalcAuto::Options load_u_inside_u_options()
  {
    const string preset = SpecUtils::append_path( g_data_dir, "rel_act/HPGe U inside U.xml" );
    shared_ptr<RelActCalcAuto::RelActAutoGuiState> state;
    // This load path also exercises the M7 headless-validity fix (the preset used to pair
    //  FwhmForm NotApplicable with a non-Fixed estimation method, which solve() rejects).
    BOOST_REQUIRE_NO_THROW( state = BatchRelActAuto::load_state_from_xml_file( preset ) );
    BOOST_REQUIRE( !!state );

    RelActCalcAuto::Options options = state->options;
    options.skew_type = PeakDef::SkewType::NoSkew;  //the simulated spectra are pure Gaussian
    return options;
  }//load_u_inside_u_options()
}//namespace


// First multi-curve regression test (2026-07 review, M8e): the simulated "easy" two-disk case
//  (1 mm DU front disk, 3 mm 93.3% HEU back disk, both behind 1 g/cm2 Fe; see
//  test_data/multicurve/README.md) fit with the shipped "HPGe U inside U" preset from DEFAULT seeds
//  must reach the truth basin - this is the M1 seed-attribution fix working end-to-end (before it,
//  this exact fit crawled a false valley into a solver-budget hard failure).
//  NOTE: a full two-curve solve - takes a few tens of seconds in Release, several minutes in Debug.
BOOST_AUTO_TEST_CASE( two_disk_default_seeding_reaches_truth )
{
  set_data_dir();

  const RelActCalcAuto::Options options = load_u_inside_u_options();
  BOOST_REQUIRE_EQUAL( options.rel_eff_curves.size(), size_t(2) );

  const shared_ptr<const SpecUtils::Measurement> foreground = load_meas( "easy_two_disk_fore.n42" );
  const shared_ptr<const SpecUtils::Measurement> background = load_meas( "easy_background_long.n42" );
  const shared_ptr<DetectorPeakResponse> drf = load_detective_x();

  RelActCalcAuto::RelActAutoSolution sol;
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, foreground, background, drf, {}, nullptr ) );

  BOOST_REQUIRE_MESSAGE( sol.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success,
                         "Solve failed: " << sol.m_error_message );

  const double chi2_dof = (sol.m_dof_data > 0) ? sol.m_chi2_data/static_cast<double>(sol.m_dof_data) : -1.0;
  BOOST_CHECK_MESSAGE( (chi2_dof > 0.0) && (chi2_dof < 2.0),
                       "chi2/dof = " << chi2_dof << " (truth basin is ~0.76)" );

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const u235 = db->nuclide( "U235" );

  // The curve with `shielded_by_other_phys_model_curve_shieldings` is the inner (back, HEU) disk.
  for( size_t curve_index = 0; curve_index < options.rel_eff_curves.size(); ++curve_index )
  {
    const bool is_inner = !options.rel_eff_curves[curve_index]
                              .shielded_by_other_phys_model_curve_shieldings.empty();
    const double truth_enrich = is_inner ? 0.933 : 0.002;

    pair<double,optional<double>> enrich;
    BOOST_REQUIRE_NO_THROW( enrich = sol.mass_enrichment_fraction( u235, curve_index ) );
    BOOST_REQUIRE( enrich.second.has_value() );

    // Within the larger of 3 sigma or an absolute floor (Ceres is threaded, so exact values vary
    //  slightly run to run; the truth-basin values are 93.8 +- 3.4 and 0.13 +- 0.04 wt%).
    //  (std::max) parenthesized: this is a "light" boost-test TU (no boost/asio in its include
    //  chain), so <windows.h>'s max macro is live on MSVC - see CLAUDE.md.
    const double tolerance = (std::max)( 3.0*(*enrich.second), (is_inner ? 0.03 : 0.0015) );
    BOOST_CHECK_MESSAGE( std::fabs(enrich.first - truth_enrich) < tolerance,
                         "curve " << curve_index << (is_inner ? " (inner)" : " (outer)")
                         << " enrichment " << enrich.first << " +- " << *enrich.second
                         << " vs truth " << truth_enrich );
  }//for( each curve )

  // The separation diagnostics must be populated for a successful two-curve fit.
  BOOST_CHECK( sol.m_curve_separation_status
               != RelActCalcAuto::RelActAutoSolution::CurveSeparationStatus::NotApplicable );
  BOOST_CHECK( !sol.m_cross_curve_correlations.empty() );
  BOOST_CHECK( sol.m_cross_curve_max_corr.has_value() );
  BOOST_CHECK_EQUAL( sol.m_evidence_purity.size(), size_t(2) );
  BOOST_CHECK( !sol.m_enrichment_diff_z.empty() );

  // The two enrichments differ hugely, so the detection statistic must see it (review: z ~ 28).
  double max_z = 0.0;
  for( const RelActCalcAuto::RelActAutoSolution::EnrichmentDiffZ &diff : sol.m_enrichment_diff_z )
    max_z = (std::max)( max_z, diff.z );
  BOOST_CHECK_MESSAGE( max_z > 5.0, "max enrichment-difference z = " << max_z );

  // The automatic merged single-curve comparison ran, and clearly prefers two curves.
  BOOST_REQUIRE( sol.m_merged_single_curve_comparison.has_value() );
  if( sol.m_merged_single_curve_comparison->valid )
  {
    BOOST_CHECK_MESSAGE( sol.m_merged_single_curve_comparison->delta_chi2
                           > 10.0*sol.m_merged_single_curve_comparison->extra_dof_of_multi,
                         "delta chi2 = " << sol.m_merged_single_curve_comparison->delta_chi2
                         << " for " << sol.m_merged_single_curve_comparison->extra_dof_of_multi
                         << " extra DOF" );
  }

}//BOOST_AUTO_TEST_CASE( two_disk_default_seeding_reaches_truth )


// Mixed physical + non-physical multi-curve config with an areal-density prior enabled: developer
//  builds used to spuriously abort on the eval residual-count assert, which omitted
//  m_phys_model_param_priors (2026-07 review, M6/A29).  The point of this test is simply that a
//  debug build gets through eval without asserting; the fit quality is irrelevant (ROIs are trimmed
//  to keep it quick).
BOOST_AUTO_TEST_CASE( mixed_physical_and_empirical_curves_with_ad_prior )
{
  set_data_dir();

  RelActCalcAuto::Options options = load_u_inside_u_options();
  BOOST_REQUIRE_EQUAL( options.rel_eff_curves.size(), size_t(2) );

  // Keep only a few ROIs so the solve is quick - this test only needs eval to run.
  if( options.rois.size() > 4 )
    options.rois.resize( 4 );

  // Find the inner (shielded) curve; it stays physical and gets the AD-bias prior; the OTHER curve
  //  becomes the non-physical LAST curve the assert needed to count priors for.
  size_t inner_index = 0;
  for( size_t i = 0; i < options.rel_eff_curves.size(); ++i )
  {
    if( !options.rel_eff_curves[i].shielded_by_other_phys_model_curve_shieldings.empty() )
      inner_index = i;
  }
  const size_t outer_index = (inner_index == 0) ? size_t(1) : size_t(0);

  // Order the curves so the non-physical one is LAST (the assert only fired on the last curve).
  if( outer_index < inner_index )
    std::swap( options.rel_eff_curves[0], options.rel_eff_curves[1] );
  RelActCalcAuto::RelEffCurveInput &phys_curve = options.rel_eff_curves[0];
  RelActCalcAuto::RelEffCurveInput &empirical_curve = options.rel_eff_curves[1];

  phys_curve.shielded_by_other_phys_model_curve_shieldings.clear();
  if( phys_curve.phys_model_self_atten )
  {
    auto self_atten = make_shared<RelActCalc::PhysicalModelShieldInput>( *phys_curve.phys_model_self_atten );
    self_atten->ad_bias.use = true;  //enables the m_phys_model_param_priors residual row(s)
    phys_curve.phys_model_self_atten = self_atten;
  }

  empirical_curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::LnX;
  empirical_curve.rel_eff_eqn_order = 2;
  empirical_curve.phys_model_self_atten.reset();
  empirical_curve.phys_model_external_atten.clear();
  empirical_curve.shielded_by_other_phys_model_curve_shieldings.clear();
  empirical_curve.mass_fraction_constraints.clear();
  empirical_curve.act_ratio_constraints.clear();

  options.same_corr_fcn_for_all_rel_eff_curves = false;
  options.same_external_shielding_for_all_rel_eff_curves = false;

  const shared_ptr<const SpecUtils::Measurement> foreground = load_meas( "easy_two_disk_fore.n42" );
  const shared_ptr<const SpecUtils::Measurement> background = load_meas( "easy_background_long.n42" );
  const shared_ptr<DetectorPeakResponse> drf = load_detective_x();

  RelActCalcAuto::RelActAutoSolution sol;
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, foreground, background, drf, {}, nullptr ) );

  // Reaching here in a developer (assert-enabled) build IS the regression check; the solve status
  //  just must not indicate the problem could not even be set up.
  BOOST_CHECK( sol.m_status != RelActCalcAuto::RelActAutoSolution::Status::NotInitiated );
  BOOST_CHECK( sol.m_status != RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem );
}//BOOST_AUTO_TEST_CASE( mixed_physical_and_empirical_curves_with_ad_prior )
