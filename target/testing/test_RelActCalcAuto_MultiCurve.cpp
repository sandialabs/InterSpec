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
#include "InterSpec/PeakFitUtils.h"
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
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, foreground, background, drf, {},
                                    PeakFitUtils::coarse_det_type( foreground, nullptr ), nullptr ) );

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

  // The per-source model counts (weight denominators of the attributed shares) are filled alongside,
  //  and the cross-curve count-attribution fractions sum to 1 per shared source.
  BOOST_CHECK_EQUAL( sol.m_source_model_counts.size(), size_t(2) );
  const vector<RelActCalcAuto::RelActAutoSolution::SourceCountAttribution> attribs
                                                                = sol.source_count_attributions();
  BOOST_CHECK( !attribs.empty() );
  for( const RelActCalcAuto::RelActAutoSolution::SourceCountAttribution &attrib : attribs )
  {
    BOOST_CHECK( attrib.total_counts > 0.0 );
    double frac_sum = 0.0;
    for( const pair<size_t,double> &curve_frac : attrib.curve_fractions )
      frac_sum += curve_frac.second;
    BOOST_CHECK_MESSAGE( std::fabs(frac_sum - 1.0) < 1.0E-9,
                         RelActCalcAuto::to_name(attrib.source) << " fractions sum to " << frac_sum );
  }

  // Retired jargon must not reach any user-facing separation text.
  const string verdict = sol.curve_separation_verdict( false );
  BOOST_CHECK_MESSAGE( (verdict.find("purity") == string::npos)
                        && (verdict.find("blended") == string::npos),
                       "verdict: " << verdict );
  for( const string &warning : sol.m_warnings )
    BOOST_CHECK_MESSAGE( (warning.find("purity") == string::npos)
                          && (warning.find("blended") == string::npos),
                         "warning: " << warning );

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
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, foreground, background, drf, {},
                                    PeakFitUtils::coarse_det_type( foreground, nullptr ), nullptr ) );

  // Reaching here in a developer (assert-enabled) build IS the regression check; the solve status
  //  just must not indicate the problem could not even be set up.
  BOOST_CHECK( sol.m_status != RelActCalcAuto::RelActAutoSolution::Status::NotInitiated );
  BOOST_CHECK( sol.m_status != RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem );
}//BOOST_AUTO_TEST_CASE( mixed_physical_and_empirical_curves_with_ad_prior )


// ---------------------------------------------------------------------------------------------
// Synthetic-solution tests of the curve-separation classification (no fits needed): populate the
//  metric members directly, run finalize_curve_separation_status(), and assert on the status /
//  detection basis / wording.  The baseline numbers replicate a real enriched-U-inside-lower-
//  enriched-U measurement (2026-08): both enrichment-difference z's ~2 (marginal), merged
//  single-curve comparison delta-chi2 = 15.17 for 2 extra DOF at chi2/dof = 571/383, inner-curve
//  U235 fully absorbed by the outer shell (attributed share 0.011) but inner U238 partially
//  resolved (share 0.46).
// ---------------------------------------------------------------------------------------------
namespace
{
  using RelActCalcAuto::RelActAutoSolution;

  RelActAutoSolution make_u_inside_u_like_solution()
  {
    set_data_dir();
    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    const SandiaDecay::Nuclide * const u235 = db->nuclide( "U235" );
    const SandiaDecay::Nuclide * const u238 = db->nuclide( "U238" );
    BOOST_REQUIRE( u235 && u238 );

    RelActAutoSolution sol;
    sol.m_options.rel_eff_curves.resize( 2 );
    sol.m_options.rel_eff_curves[0].name = "Inner";
    sol.m_options.rel_eff_curves[1].name = "Outer";

    sol.m_chi2_data = 571.0;
    sol.m_dof_data = 383;
    sol.m_r2 = 0.999;
    sol.m_num_rank_deficient_dirs = 0;
    sol.m_jacobian_condition_number = 148.0;

    // Inner U235 is fully absorbed by the outer shell (expected physics); inner U238 is anchored.
    sol.m_evidence_purity.assign( 2, {} );
    sol.m_evidence_purity[0][u235] = 0.011;
    sol.m_evidence_purity[0][u238] = 0.46;
    sol.m_evidence_purity[1][u235] = 0.94;
    sol.m_evidence_purity[1][u238] = 0.54;

    sol.m_source_model_counts.assign( 2, {} );
    sol.m_source_model_counts[0][u235] = 500.0;
    sol.m_source_model_counts[0][u238] = 30000.0;
    sol.m_source_model_counts[1][u235] = 45000.0;
    sol.m_source_model_counts[1][u238] = 35000.0;

    RelActAutoSolution::CrossCurveCorrelation corr;
    corr.curve_a = 0;
    corr.curve_b = 1;
    corr.param_a = corr.param_b = "Act(U235)";
    corr.correlation = 0.5;
    sol.m_cross_curve_correlations.push_back( corr );
    sol.m_cross_curve_max_corr = corr;

    RelActAutoSolution::EnrichmentDiffZ u235_diff;
    u235_diff.curve_a = 0;
    u235_diff.curve_b = 1;
    u235_diff.nuclide = RelActCalcAuto::SrcVariant( u235 );
    u235_diff.enrichment_a = 0.775;
    u235_diff.sigma_a = 0.18;
    u235_diff.enrichment_b = 0.00896;
    u235_diff.sigma_b = 0.33;
    u235_diff.z = std::fabs( u235_diff.enrichment_a - u235_diff.enrichment_b )
                  / std::sqrt( u235_diff.sigma_a*u235_diff.sigma_a + u235_diff.sigma_b*u235_diff.sigma_b );
    u235_diff.reliable = true;

    RelActAutoSolution::EnrichmentDiffZ u238_diff = u235_diff;
    u238_diff.nuclide = RelActCalcAuto::SrcVariant( u238 );
    u238_diff.enrichment_a = 0.217;
    u238_diff.enrichment_b = 0.991;
    u238_diff.z = std::fabs( u238_diff.enrichment_a - u238_diff.enrichment_b )
                  / std::sqrt( u238_diff.sigma_a*u238_diff.sigma_a + u238_diff.sigma_b*u238_diff.sigma_b );

    sol.m_enrichment_diff_z.push_back( u235_diff );
    sol.m_enrichment_diff_z.push_back( u238_diff );

    RelActAutoSolution::MergedCurveComparison merged;
    merged.valid = true;
    merged.multi_chi2_data = 571.0;
    merged.multi_dof_data = 383;
    merged.merged_chi2_data = 571.0 + 15.17;
    merged.merged_dof_data = 385;
    merged.delta_chi2 = 15.17;
    merged.extra_dof_of_multi = 2;
    sol.m_merged_single_curve_comparison = merged;

    return sol;
  }//make_u_inside_u_like_solution()

  bool contains( const string &haystack, const string &needle )
  {
    return (haystack.find(needle) != string::npos);
  }
}//namespace


// The motivating real-world case: two marginal z's (~2) plus a decisive single-curve rejection
//  (delta-chi2 = 15.17 >> the 3-sigma-scaled bar of ~11.9) must count as detected-distinct via the
//  corroborated tier, and the expected-physics blocked inner U235 must NOT drive the headline to
//  "Poorly separated" - it becomes a per-nuclide caveat.
BOOST_AUTO_TEST_CASE( marginal_z_corroborated_by_merged_rejection_is_detected )
{
  RelActAutoSolution sol = make_u_inside_u_like_solution();
  sol.finalize_curve_separation_status();

  BOOST_CHECK( sol.curves_distinct_basis()
               == RelActAutoSolution::CurveDistinctBasis::ZCorroboratedByMerged );
  BOOST_CHECK( sol.curves_detected_distinct() );
  BOOST_CHECK( sol.m_curve_separation_status
               == RelActAutoSolution::CurveSeparationStatus::WellSeparated );
  BOOST_CHECK_EQUAL( string(sol.curve_separation_display()), string("Separated") );

  BOOST_REQUIRE( sol.m_merged_single_curve_comparison.has_value() );
  BOOST_CHECK( !sol.m_merged_single_curve_comparison->single_curve_adequate );

  const string verdict = sol.curve_separation_verdict( false );
  BOOST_CHECK_MESSAGE( contains( verdict, "genuinely different" ), "verdict: " << verdict );
  // The blocked inner U235 gets the per-nuclide caveat, with the count attribution...
  BOOST_CHECK_MESSAGE( contains( verdict, "peaks of its own" ), "verdict: " << verdict );
  BOOST_CHECK_MESSAGE( contains( verdict, "'Inner'" ), "verdict: " << verdict );
  // ...and the retired jargon stays retired.
  BOOST_CHECK_MESSAGE( !contains( verdict, "purity" ), "verdict: " << verdict );
  BOOST_CHECK_MESSAGE( !contains( verdict, "blended" ), "verdict: " << verdict );
  BOOST_CHECK_MESSAGE( !contains( verdict, "NOT well separated" ), "verdict: " << verdict );

  for( const string &warning : sol.m_warnings )
  {
    BOOST_CHECK_MESSAGE( !contains( warning, "purity" ) && !contains( warning, "blended" ),
                         "warning: " << warning );
  }
}//BOOST_AUTO_TEST_CASE( marginal_z_corroborated_by_merged_rejection_is_detected )


// S1 guard (2026-07 validation): an inflated z = 4.85 on genuinely identical stacked disks, with
//  the merged fit doing about as well (delta-chi2 = 5.5 for 2 DOF), must stay vetoed - removing the
//  old have-reliable-z early-out must not resurrect it through the new tiers.
BOOST_AUTO_TEST_CASE( s1_inflated_z_still_vetoed_by_adequate_merged_fit )
{
  RelActAutoSolution sol = make_u_inside_u_like_solution();
  sol.m_chi2_data = 100.0;
  sol.m_dof_data = 100;
  sol.m_enrichment_diff_z.resize( 1 );
  sol.m_enrichment_diff_z[0].z = 4.85;
  sol.m_merged_single_curve_comparison->delta_chi2 = 5.5;
  // Everything blended, as in the S1 collapse minimum.
  for( auto &curve_shares : sol.m_evidence_purity )
    for( auto &src_share : curve_shares )
      src_share.second = 0.016;

  sol.finalize_curve_separation_status();

  BOOST_CHECK( sol.curves_distinct_basis() == RelActAutoSolution::CurveDistinctBasis::None );
  BOOST_CHECK( !sol.curves_detected_distinct() );
  BOOST_CHECK( sol.z_detection_vetoed_by_merged() );
  BOOST_CHECK( sol.m_curve_separation_status
               == RelActAutoSolution::CurveSeparationStatus::Degenerate );

  // The z table row must carry the veto annotation (headline-vs-evidence consistency, S8).
  const string annotation = sol.z_row_annotation( sol.m_enrichment_diff_z[0] );
  BOOST_CHECK_MESSAGE( contains( annotation, "not treated as a detection" ),
                       "annotation: " << annotation );
}//BOOST_AUTO_TEST_CASE( s1_inflated_z_still_vetoed_by_adequate_merged_fit )


// A marginal z with the merged fit adequate (delta-chi2 within the 3-sigma-scaled bar) is NOT a
//  detection; same just below the bar.
BOOST_AUTO_TEST_CASE( marginal_z_without_merged_rejection_is_not_detected )
{
  {
    RelActAutoSolution sol = make_u_inside_u_like_solution();
    sol.m_chi2_data = 100.0;
    sol.m_dof_data = 100;
    sol.m_merged_single_curve_comparison->delta_chi2 = 7.0;  //bar is 2 + 3*sqrt(4) = 8
    sol.finalize_curve_separation_status();
    BOOST_CHECK( sol.curves_distinct_basis() == RelActAutoSolution::CurveDistinctBasis::None );
    BOOST_CHECK( sol.m_merged_single_curve_comparison->single_curve_adequate );
  }

  {
    // Boundary: just below the bar still fails.
    RelActAutoSolution sol = make_u_inside_u_like_solution();
    sol.m_chi2_data = 100.0;
    sol.m_dof_data = 100;
    sol.m_merged_single_curve_comparison->delta_chi2 = 7.99;
    sol.finalize_curve_separation_status();
    BOOST_CHECK( sol.curves_distinct_basis() == RelActAutoSolution::CurveDistinctBasis::None );
  }
}//BOOST_AUTO_TEST_CASE( marginal_z_without_merged_rejection_is_not_detected )


// The pre-existing fallback: no reliable z at all, delta-chi2 above the 5-sigma-scaled bar
//  (2 + 5*sqrt(4) = 12 at chi2/dof <= 1) detects as MergedOnly.
BOOST_AUTO_TEST_CASE( merged_only_detection_without_z_table )
{
  RelActAutoSolution sol = make_u_inside_u_like_solution();
  sol.m_chi2_data = 100.0;
  sol.m_dof_data = 100;
  sol.m_enrichment_diff_z.clear();
  sol.m_merged_single_curve_comparison->delta_chi2 = 13.0;
  sol.finalize_curve_separation_status();

  BOOST_CHECK( sol.curves_distinct_basis() == RelActAutoSolution::CurveDistinctBasis::MergedOnly );
  BOOST_CHECK( sol.curves_detected_distinct() );
}//BOOST_AUTO_TEST_CASE( merged_only_detection_without_z_table )


// A decisive merged rejection with only a tiny z (< 1.5) detects the CURVES as distinct, but the
//  verdict must say plainly the compositions themselves are not distinguished.
BOOST_AUTO_TEST_CASE( decisive_merged_with_tiny_z_detects_curves_not_composition )
{
  RelActAutoSolution sol = make_u_inside_u_like_solution();
  sol.m_chi2_data = 100.0;
  sol.m_dof_data = 100;
  for( RelActAutoSolution::EnrichmentDiffZ &diff : sol.m_enrichment_diff_z )
    diff.z = 1.0;
  sol.m_merged_single_curve_comparison->delta_chi2 = 20.0;
  sol.finalize_curve_separation_status();

  BOOST_CHECK( sol.curves_distinct_basis() == RelActAutoSolution::CurveDistinctBasis::MergedOnly );
  const string verdict = sol.curve_separation_verdict( false );
  BOOST_CHECK_MESSAGE( contains( verdict, "NOT distinguished" ), "verdict: " << verdict );
}//BOOST_AUTO_TEST_CASE( decisive_merged_with_tiny_z_detects_curves_not_composition )


// A whole curve with no resolved evidence (every attributed share < 0.3) still downgrades a
//  detected pair to PoorlySeparated, with the unanchored-curve trigger text.
BOOST_AUTO_TEST_CASE( unanchored_curve_downgrades_detected_pair )
{
  RelActAutoSolution sol = make_u_inside_u_like_solution();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const u238 = db->nuclide( "U238" );
  sol.m_evidence_purity[0][u238] = 0.05;  //now NO source on 'Inner' reaches 0.3
  for( RelActAutoSolution::EnrichmentDiffZ &diff : sol.m_enrichment_diff_z )
    diff.z = 5.0;  //clear Tier-1 detection (delta-chi2 = 15.17 is above the veto bar of 8)
  sol.finalize_curve_separation_status();

  BOOST_CHECK( sol.curves_distinct_basis() == RelActAutoSolution::CurveDistinctBasis::ZScore );
  BOOST_CHECK( sol.m_curve_separation_status
               == RelActAutoSolution::CurveSeparationStatus::PoorlySeparated );
  BOOST_CHECK_EQUAL( string(sol.curve_separation_display()),
                     string("Distinct curves - see per-nuclide notes") );

  const string trigger = sol.curve_separation_trigger_text();
  BOOST_CHECK_MESSAGE( contains( trigger, "no nuclide on 'Inner' has individually resolved peaks" ),
                       "trigger: " << trigger );

  BOOST_REQUIRE( !sol.m_warnings.empty() );
  bool warning_names_trigger = false;
  for( const string &warning : sol.m_warnings )
    warning_names_trigger = (warning_names_trigger || contains( warning, "individually resolved" ));
  BOOST_CHECK( warning_names_trigger );
}//BOOST_AUTO_TEST_CASE( unanchored_curve_downgrades_detected_pair )


// S4 guard: a fit that does not describe the data (R^2 below the floor) can never claim more than
//  PoorlySeparated, regardless of detection evidence.
BOOST_AUTO_TEST_CASE( poor_fit_quality_overrides_detection )
{
  RelActAutoSolution sol = make_u_inside_u_like_solution();
  sol.m_r2 = 0.5;
  sol.finalize_curve_separation_status();

  BOOST_CHECK( sol.poor_fit_quality() );
  BOOST_CHECK( sol.m_curve_separation_status
               == RelActAutoSolution::CurveSeparationStatus::PoorlySeparated );
  const string trigger = sol.curve_separation_trigger_text();
  BOOST_CHECK_MESSAGE( contains( trigger, "does not describe the data" ), "trigger: " << trigger );

  // A failed fit forfeits the "distinct" label/verdict even though a detection basis exists - the
  //  evidence comes from the same fit the verdict declares meaningless.  The verdict is the failed-
  //  fit trigger alone, with no "genuinely different" claim and no endorsement of uncertainties.
  BOOST_CHECK_EQUAL( string(sol.curve_separation_display()), string("Poorly separated") );
  const string verdict = sol.curve_separation_verdict( false );
  BOOST_CHECK_MESSAGE( contains( verdict, "does not describe the data" ), "verdict: " << verdict );
  BOOST_CHECK_MESSAGE( !contains( verdict, "genuinely" ) && !contains( verdict, "uncertainties" ),
                       "verdict: " << verdict );
  bool have_cannot_assess_warning = false;
  for( const string &warning : sol.m_warnings )
    have_cannot_assess_warning = (have_cannot_assess_warning
                                  || contains( warning, "cannot be assessed" ));
  BOOST_CHECK( have_cannot_assess_warning );
}//BOOST_AUTO_TEST_CASE( poor_fit_quality_overrides_detection )


// The per-source count-attribution accessor: fractions per source sum to 1, only shared sources
//  are reported, and the fractions match the model counts.
BOOST_AUTO_TEST_CASE( source_count_attributions_fractions )
{
  RelActAutoSolution sol = make_u_inside_u_like_solution();
  const vector<RelActAutoSolution::SourceCountAttribution> attribs = sol.source_count_attributions();
  BOOST_REQUIRE_EQUAL( attribs.size(), size_t(2) );  //U235 and U238, each on both curves

  for( const RelActAutoSolution::SourceCountAttribution &attrib : attribs )
  {
    BOOST_REQUIRE_EQUAL( attrib.curve_fractions.size(), size_t(2) );
    double frac_sum = 0.0;
    for( const pair<size_t,double> &curve_frac : attrib.curve_fractions )
      frac_sum += curve_frac.second;
    BOOST_CHECK_MESSAGE( std::fabs(frac_sum - 1.0) < 1.0E-9, "fractions sum to " << frac_sum );
  }

  // U235: 500 of 45500 total on 'Inner' (index 0).
  const RelActAutoSolution::SourceCountAttribution &u235_attrib = attribs[0];
  BOOST_CHECK_EQUAL( RelActCalcAuto::to_name(u235_attrib.source), string("U235") );
  BOOST_CHECK_MESSAGE( std::fabs(u235_attrib.total_counts - 45500.0) < 1.0E-6,
                       "U235 total = " << u235_attrib.total_counts );
  for( const pair<size_t,double> &curve_frac : u235_attrib.curve_fractions )
  {
    if( curve_frac.first == 0 )
      BOOST_CHECK_MESSAGE( std::fabs(curve_frac.second - 500.0/45500.0) < 1.0E-9,
                           "Inner U235 fraction = " << curve_frac.second );
  }
}//BOOST_AUTO_TEST_CASE( source_count_attributions_fractions )
