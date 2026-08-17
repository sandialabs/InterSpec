/* InterSpec: an application to analyze spectral gamma radiation data.

 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.

 This library is free software; you can redistribute it and/or modify it under
 the terms of the GNU Lesser General Public License as published by the Free
 Software Foundation; either version 2.1 of the License, or (at your option)
 any later version.
 */

#include "InterSpec_config.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#define BOOST_TEST_MODULE RelActCalcAuto_Pu610775_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/BatchRelActAuto.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/SpecMeas.h"

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

using namespace std;
using namespace boost::unit_test;

static_assert( USE_REL_ACT_TOOL, "Compile-time support for Rel Act tool is required." );

namespace
{
string g_data_dir;
string g_fixture_dir;
using json = nlohmann::json;

struct PuCase
{
  const char *id;
  const char *spectrum;
  const char *background;
  array<double,5> truth;
  bool correlation_in_range;
};

const array<PuCase,3> sm_cases{{
  { "Pu93", "48xmDT-CBNMPu93-10cm-Fe16-Cd00-0.spe",
    "mDT-BKG-2017-10-11_summed.n42",
    {{0.000091612433,0.935920939673,0.063100423177,0.000490931012,0.000396093704}},
    false },
  { "Pu70", "48xmDT-CBNMPu70-10cm-Fe16-Cd00-0.spe",
    "mDT-BKG-2017-09-20_summed.n42",
    {{0.006924214632,0.767646102488,0.191082614349,0.012580638132,0.021766430398}},
    true },
  { "Pu61", "48xmDT-CBNMPu61-10cm-Fe16-Cd00-0.spe",
    "mDT-BKG-2017-10-11_summed.n42",
    {{0.009854691411,0.662259574515,0.268432937204,0.015008318693,0.044444478178}},
    true }
}};

void initialize_paths()
{
  static bool initialized = false;
  if( initialized )
    return;

  string test_file_dir;
  const int argc = framework::master_test_suite().argc;
  char ** const argv = framework::master_test_suite().argv;
  for( int i = 1; i < argc; ++i )
  {
    const string arg = argv[i];
    if( SpecUtils::istarts_with(arg,"--datadir=") )
      g_data_dir = arg.substr(10);
    if( SpecUtils::istarts_with(arg,"--testfiledir=") )
      test_file_dir = arg.substr(14);
  }
  SpecUtils::ireplace_all(g_data_dir,"%20"," ");
  SpecUtils::ireplace_all(test_file_dir,"%20"," ");

  if( g_data_dir.empty() )
    for( const char *candidate : {"data","../data","../../data","../../../data"} )
      if( SpecUtils::is_file(SpecUtils::append_path(candidate,"sandia.decay.xml")) )
      {
        g_data_dir = candidate;
        break;
      }
  if( test_file_dir.empty() )
    for( const char *candidate : {"test_data","../test_data","../../test_data",
                                   "../../../target/testing/test_data"} )
      if( SpecUtils::is_directory(SpecUtils::append_path(candidate,
                                      "rel_act_auto/pu_610_775")) )
      {
        test_file_dir = candidate;
        break;
      }

  BOOST_REQUIRE( SpecUtils::is_file(SpecUtils::append_path(g_data_dir,"sandia.decay.xml")) );
  g_fixture_dir = SpecUtils::append_path(test_file_dir,"rel_act_auto/pu_610_775");
  BOOST_REQUIRE( SpecUtils::is_directory(g_fixture_dir) );
  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory(g_data_dir) );
  BOOST_REQUIRE( DecayDataBaseServer::database() );
  initialized = true;
}

const json &fixture_manifest()
{
  static json manifest;
  static bool loaded = false;
  if( !loaded )
  {
    const string path = SpecUtils::append_path(g_fixture_dir,"manifest.json");
    ifstream input(path);
    if( !input )
      throw runtime_error("Could not open Pu truth manifest " + path);
    manifest = json::parse(input);
    loaded = true;
  }
  return manifest;
}

double effective_age_years( const PuCase &test_case )
{
  for( const json &manifest_case : fixture_manifest().at("cases") )
  {
    if( manifest_case.at("spectrum_file").get<string>() != test_case.spectrum )
      continue;
    const double age = manifest_case.at("fit_acceptance").at("known_age_diagnostic")
                                    .at("effective_age_years_365_25d").get<double>();
    if( !std::isfinite(age) || (age <= 0.0) )
      throw runtime_error(string("Invalid manifest effective age for ") + test_case.id);
    return age;
  }
  throw runtime_error(string("Pu truth manifest has no case for ") + test_case.spectrum);
}

shared_ptr<const SpecUtils::Measurement> load_measurement( const string &filename )
{
  auto file = make_shared<SpecMeas>();
  const string path = SpecUtils::append_path(g_fixture_dir,filename);
  BOOST_REQUIRE_MESSAGE( file->load_file(path,SpecUtils::ParserType::Auto),
                         "Could not load " << path );
  for( const auto &measurement : file->measurements() )
    if( measurement && measurement->num_gamma_channels() )
      return measurement;
  BOOST_REQUIRE_MESSAGE( false, "No gamma measurement in " << path );
  return nullptr;
}

RelActCalcAuto::Options preset_options()
{
  const string path = SpecUtils::append_path(g_data_dir,
                                              "rel_act/HPGe Pu (610-775 keV).xml");
  shared_ptr<RelActCalcAuto::RelActAutoGuiState> state;
  BOOST_REQUIRE_NO_THROW( state = BatchRelActAuto::load_state_from_xml_file(path) );
  BOOST_REQUIRE( state );
  return state->options;
}

RelActCalcAuto::RelActAutoSolution solve_case( const PuCase &test_case,
                                               const bool known_age,
                                               const bool auto_profile )
{
  RelActCalcAuto::Options options = preset_options();
  options.auto_profile_weak_mass_fractions = auto_profile;
  // These grids are the acceptance case for the expensive path: multi-start basin search plus
  // bounded profile-likelihood intervals.  They are gated out of ctest for exactly that reason.
  options.robust_solve = auto_profile;
  if( known_age )
  {
    for( RelActCalcAuto::NucInputInfo &input : options.rel_eff_curves.at(0).nuclides )
    {
      input.age = effective_age_years(test_case) * PhysicalUnits::year;
      input.fit_age = false;
      input.fit_age_min.reset();
      input.fit_age_max.reset();
    }
  }

  const auto foreground = load_measurement(test_case.spectrum);
  const auto background = load_measurement(test_case.background);
  RelActCalcAuto::RelActAutoSolution solution;
  BOOST_REQUIRE_NO_THROW( solution = RelActCalcAuto::solve(
      options,foreground,background,nullptr,{},PeakFitUtils::CoarseResolutionType::High,nullptr) );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status(solution.m_status),
                         test_case.id << " failed: " << solution.m_error_message );
  return solution;
}

array<const SandiaDecay::Nuclide *,5> pu_nuclides()
{
  const auto *db = DecayDataBaseServer::database();
  return {{db->nuclide("Pu238"),db->nuclide("Pu239"),db->nuclide("Pu240"),
           db->nuclide("Pu241"),db->nuclide("Pu242")}};
}

array<double,5> fitted_fractions( const RelActCalcAuto::RelActAutoSolution &solution )
{
  const auto nuclides = pu_nuclides();
  array<double,5> result{};
  for( size_t i = 0; i < result.size(); ++i )
  {
    BOOST_REQUIRE( nuclides[i] );
    result[i] = solution.mass_enrichment_result(nuclides[i],0).fraction;
    BOOST_REQUIRE( std::isfinite(result[i]) );
    BOOST_CHECK_GE( result[i], 0.0 );
    BOOST_CHECK_LE( result[i], 1.0 );
  }
  return result;
}

struct FittedAgeDiagnostic
{
  double age = -1.0;
  double one_sigma = -1.0;
};

/** Check the high-energy preset's fitted age as a numerical diagnostic only.  The certificate
 acquisition-time propagation supplies isotope-composition truth, not a certified chemical-
 separation date, so this deliberately checks finite/bounded reporting rather than age accuracy. */
FittedAgeDiagnostic check_free_age_diagnostic(
                         const RelActCalcAuto::RelActAutoSolution &solution,
                         const string &case_id )
{
  BOOST_REQUIRE_EQUAL( solution.m_options.rel_eff_curves.size(),1u );
  BOOST_REQUIRE_EQUAL( solution.m_rel_activities.size(),1u );

  FittedAgeDiagnostic answer;
  size_t fitted_age_rows = 0;
  for( const RelActCalcAuto::NuclideRelAct &output : solution.m_rel_activities[0] )
  {
    if( !output.age_was_fit )
      continue;
    ++fitted_age_rows;

    const auto input = std::find_if(
        solution.m_options.rel_eff_curves[0].nuclides.begin(),
        solution.m_options.rel_eff_curves[0].nuclides.end(),
        [&output]( const RelActCalcAuto::NucInputInfo &candidate ) {
          return candidate.source == output.source;
        } );
    BOOST_REQUIRE( input != solution.m_options.rel_eff_curves[0].nuclides.end() );

    // With no explicit upper bound, solve setup caps fitted ages at 120 years.  Use the explicit
    // input bounds when present and the documented setup default otherwise.
    const double lower = input->fit_age_min.value_or(0.0);
    const double upper = input->fit_age_max.value_or(120.0*PhysicalUnits::year);
    const double bound_tolerance = 64.0*std::numeric_limits<double>::epsilon()
        *(std::max)({1.0,std::fabs(lower),std::fabs(upper)});
    BOOST_CHECK_MESSAGE( std::isfinite(output.age)
                           && (output.age >= lower-bound_tolerance)
                           && (output.age <= upper+bound_tolerance),
                         case_id << " fitted " << output.name() << " age "
                                 << output.age/PhysicalUnits::year
                                 << " y is outside [" << lower/PhysicalUnits::year << ", "
                                 << upper/PhysicalUnits::year << "] y" );
    BOOST_CHECK_MESSAGE( std::isfinite(output.age_uncertainty)
                           && (output.age_uncertainty >= 0.0),
                         case_id << " fitted " << output.name()
                                 << " age uncertainty is " << output.age_uncertainty );

    if( answer.age < 0.0 )
    {
      answer.age = output.age;
      answer.one_sigma = output.age_uncertainty;
    }else
    {
      // NucsOfElSameAge is enabled by the preset, so all four displayed Pu rows must decode the
      // same controlling coordinate and uncertainty.
      BOOST_CHECK_SMALL( output.age-answer.age,
                         1.0e-12*(std::max)(1.0,std::fabs(answer.age)) );
      BOOST_CHECK_SMALL( output.age_uncertainty-answer.one_sigma,
                         1.0e-12*(std::max)(1.0,std::fabs(answer.one_sigma)) );
    }
  }

  BOOST_REQUIRE_EQUAL( fitted_age_rows,4u );
  BOOST_REQUIRE( std::isfinite(answer.age) && (answer.age >= 0.0) );
  BOOST_REQUIRE( std::isfinite(answer.one_sigma) && (answer.one_sigma >= 0.0) );
  return answer;
}

bool html_report_resources_available()
{
  const char * const needed[] = {
    "static_text/auto_rel_act_report.tmplt.html",
    "d3.v3.min.js", "SpectrumChartD3.js", "SpectrumChartD3.css",
    "RelEffPlot.js", "RelEffPlot.css"
  };
  for( const char * const filename : needed )
    if( !SpecUtils::is_file(SpecUtils::append_path("InterSpec_resources",filename)) )
      return false;
  return true;
}

string html_table_with_caption( const string &html, const string &caption )
{
  const size_t caption_pos = html.find(caption);
  if( caption_pos == string::npos )
    return {};
  const size_t table_begin = html.rfind("<table",caption_pos);
  const size_t table_end = html.find("</table>",caption_pos);
  if( (table_begin == string::npos) || (table_end == string::npos) )
    return {};
  return html.substr(table_begin,table_end + string("</table>").size() - table_begin);
}

string html_nuclide_row( const string &table, const string &symbol )
{
  const size_t row_begin = table.find("<tr><td>" + symbol);
  if( row_begin == string::npos )
    return {};
  const size_t row_end = table.find("</tr>",row_begin);
  if( row_end == string::npos )
    return {};
  return table.substr(row_begin,row_end + string("</tr>").size() - row_begin);
}

void check_interval_values( const string &row,
                            const RelActCalcAuto::RelActAutoSolution::MassFractionProfileInterval &interval,
                            const string &where )
{
  const string lower = SpecUtils::printCompact(100.0*interval.lower,4);
  const string upper = SpecUtils::printCompact(100.0*interval.upper,4);
  BOOST_CHECK_MESSAGE( row.find(lower) != string::npos,
                       where << " omitted lower endpoint " << lower << "%: " << row );
  BOOST_CHECK_MESSAGE( row.find(upper) != string::npos,
                       where << " omitted upper endpoint " << upper << "%: " << row );
}

/** Exercises the two legacy, in-code report writers with a live Pu cost functor, but no additional
 optimization.  A synthetic profile makes the expected 68/95 values independent of optimizer
 behavior; a second copy has its covariance enlarged to simulate the hundreds-percent local
 Gaussian that first-party surfaces must classify instead of rendering as a symmetric error. */
void check_pu_report_uncertainty_contract(
                  const RelActCalcAuto::RelActAutoSolution &solution,
                  const SandiaDecay::Nuclide * const target )
{
  using Solution = RelActCalcAuto::RelActAutoSolution;
  BOOST_REQUIRE( target );

  const Solution::MassFractionResult baseline = solution.mass_enrichment_result(target,0);
  BOOST_REQUIRE( baseline.covariance_one_sigma.has_value() );
  BOOST_REQUIRE( std::isfinite(*baseline.covariance_one_sigma) );
  BOOST_REQUIRE_GT( *baseline.covariance_one_sigma,0.0 );

  Solution::MassFractionProfileInterval interval68;
  interval68.confidence_level = 0.6827;
  interval68.delta_chi2 = solution.m_cov_scale;
  interval68.lower = (std::max)(0.0,baseline.fraction - 0.007321);
  interval68.upper = (std::min)(1.0,baseline.fraction + 0.004567);
  interval68.lower_kind = Solution::MassFractionProfileEndpointKind::LikelihoodCrossing;
  interval68.upper_kind = Solution::MassFractionProfileEndpointKind::LikelihoodCrossing;

  Solution::MassFractionProfileInterval interval95;
  interval95.confidence_level = 0.95;
  interval95.delta_chi2 = solution.m_cov_scale * 3.841458820694124;
  interval95.lower = (std::max)(0.0,baseline.fraction - 0.01543);
  interval95.upper = (std::min)(1.0,baseline.fraction + 0.01234);
  interval95.lower_kind = Solution::MassFractionProfileEndpointKind::PhysicalLimit;
  interval95.upper_kind = Solution::MassFractionProfileEndpointKind::LikelihoodCrossing;

  Solution with_profile = solution;
  if( with_profile.m_mass_fraction_profiles.empty() )
    with_profile.m_mass_fraction_profiles.resize(1);
  Solution::MassFractionProfileResult profile;
  profile.status = Solution::MassFractionProfileStatus::BoundaryLimited;
  profile.reason = Solution::MassFractionProfileReason::Forced;
  profile.intervals = {interval68,interval95};
  profile.num_fits = 7;
  with_profile.m_mass_fraction_profiles[0][target->symbol] = profile;

  // The structured result itself is the schema oracle for the injected intervals.
  const Solution::MassFractionResult structured = with_profile.mass_enrichment_result(target,0);
  BOOST_REQUIRE( structured.profile.has_value() );
  BOOST_REQUIRE_EQUAL( structured.profile->intervals.size(),2u );
  BOOST_CHECK_EQUAL( structured.profile->intervals[0].delta_chi2,interval68.delta_chi2 );
  BOOST_CHECK( structured.profile->intervals[1].lower_kind
               == Solution::MassFractionProfileEndpointKind::PhysicalLimit );
  BOOST_CHECK( structured.profile->intervals[1].upper_kind
               == Solution::MassFractionProfileEndpointKind::LikelihoodCrossing );

  // Make the retained Gaussian about 275 percentage points wide, then remove the profile.  This
  // recreates the HPGe-Pu failure mode without depending on a particular optimizer basin.
  Solution unreliable = solution;
  BOOST_REQUIRE( !unreliable.m_covariance.empty() );
  const double wanted_sigma = 2.75;
  const double covariance_scale
      = std::pow(wanted_sigma / *baseline.covariance_one_sigma,2.0);
  BOOST_REQUIRE( std::isfinite(covariance_scale) );
  for( vector<double> &row : unreliable.m_covariance )
    for( double &entry : row )
      entry *= covariance_scale;
  if( !unreliable.m_mass_fraction_profiles.empty() )
    unreliable.m_mass_fraction_profiles[0].erase(target->symbol);

  const Solution::MassFractionResult weak = unreliable.mass_enrichment_result(target,0);
  BOOST_REQUIRE( weak.covariance_one_sigma.has_value() );
  BOOST_REQUIRE_GT( *weak.covariance_one_sigma,1.0 );
  BOOST_REQUIRE( weak.covariance_quality != Solution::MassFractionCovarianceQuality::Usable );
  BOOST_CHECK( !weak.profile.has_value() );

  stringstream summary;
  BOOST_REQUIRE_NO_THROW( unreliable.print_summary(summary) );
  const string summary_text = summary.str();
  const size_t pu_section = summary_text.find("Plutonium Enrichment");
  BOOST_REQUIRE( pu_section != string::npos );
  const string row_prefix = "    " + target->symbol + ":";
  const size_t summary_row_begin = summary_text.find(row_prefix,pu_section);
  BOOST_REQUIRE( summary_row_begin != string::npos );
  const size_t summary_row_end = summary_text.find('\n',summary_row_begin);
  const string summary_row = summary_text.substr(summary_row_begin,
      summary_row_end == string::npos ? string::npos : summary_row_end-summary_row_begin);
  BOOST_CHECK_MESSAGE( summary_row.find(" ± ") == string::npos,
                       "Structured unreliable covariance leaked as raw symmetric uncertainty: "
                       << summary_row );
  BOOST_CHECK_MESSAGE( SpecUtils::icontains(summary_row,"unreliable")
                         || SpecUtils::icontains(summary_row,"unavailable"),
                       "The suppressed local Gaussian was not classified in text output: "
                       << summary_row );

  if( !html_report_resources_available() )
  {
    BOOST_TEST_MESSAGE( "Skipping in-code HTML uncertainty assertions: report resources unavailable." );
    return;
  }

  stringstream weak_html_stream;
  BOOST_REQUIRE_NO_THROW( unreliable.print_html_report(weak_html_stream) );
  const string weak_table = html_table_with_caption(weak_html_stream.str(),
                                                     "Plutonium mass fractions");
  BOOST_REQUIRE( !weak_table.empty() );
  BOOST_CHECK_MESSAGE( !SpecUtils::icontains(weak_table,"1-σ Uncert"),
                       "The bounded 68% column is still labeled as a raw Gaussian 1-sigma." );
  BOOST_CHECK_MESSAGE( SpecUtils::icontains(weak_table,"68%"),
                       "Corrected Pu table has no 68% interval header." );
  BOOST_CHECK_MESSAGE( SpecUtils::icontains(weak_table,"95%"),
                       "Corrected Pu table has no detailed 95% profile header." );

  const string weak_row = html_nuclide_row(weak_table,target->symbol);
  BOOST_REQUIRE( !weak_row.empty() );
  const string raw_sigma_cell = "<td>"
      + SpecUtils::printCompact(100.0*(*weak.covariance_one_sigma),4) + "</td>";
  BOOST_CHECK_MESSAGE( weak_row.find(raw_sigma_cell) == string::npos,
                       "Structured unreliable covariance leaked into HTML as "
                       << raw_sigma_cell << ": " << weak_row );
  BOOST_CHECK_MESSAGE( SpecUtils::icontains(weak_row,"unreliable")
                         || SpecUtils::icontains(weak_row,"unavailable")
                         || (weak_row.find("--") != string::npos),
                       "The suppressed local Gaussian was not classified in HTML: " << weak_row );

  stringstream profiled_html_stream;
  BOOST_REQUIRE_NO_THROW( with_profile.print_html_report(profiled_html_stream) );
  const string profiled_html = profiled_html_stream.str();

  const string corrected_table = html_table_with_caption(profiled_html,
                                                          "Plutonium mass fractions");
  BOOST_REQUIRE( !corrected_table.empty() );
  BOOST_CHECK( SpecUtils::icontains(corrected_table,"68%") );
  BOOST_CHECK( SpecUtils::icontains(corrected_table,"95%") );
  const string corrected_row = html_nuclide_row(corrected_table,target->symbol);
  BOOST_REQUIRE( !corrected_row.empty() );
  check_interval_values(corrected_row,interval68,"corrected Pu 68% column");
  check_interval_values(corrected_row,interval95,"corrected Pu 95% column");

  const string detailed_table = html_table_with_caption(profiled_html,
                                           "Relative activities and mass fractions");
  BOOST_REQUIRE( !detailed_table.empty() );
  BOOST_CHECK_MESSAGE( SpecUtils::icontains(detailed_table,"68%"),
                       "Detailed mass-fraction table has no 68% interval header." );
  BOOST_CHECK_MESSAGE( SpecUtils::icontains(detailed_table,"95%"),
                       "Detailed mass-fraction table has no 95% profile header." );
  const string detailed_row = html_nuclide_row(detailed_table,target->symbol);
  BOOST_REQUIRE( !detailed_row.empty() );
  check_interval_values(detailed_row,interval68,"detailed 68% column");
  check_interval_values(detailed_row,interval95,"detailed 95% column");
}
}//namespace


BOOST_AUTO_TEST_CASE( known_effective_age_accuracy_grid )
{
  initialize_paths();
  static const array<double,5> max_error{{0.005,0.03,0.03,0.015,0.015}};
  // When ByPu239Only is outside its validation domain, every corrected component inherits the
  // unsupported Pu-242 extrapolation through renormalization.  The directly fitted pre-correlation
  // composition remains the meaningful conditional accuracy diagnostic (Pu-242 is not fitted).
  static const array<double,4> pu93_direct_truth{{
      0.000091648734,0.936291798960,0.063125426761,0.000491125544 }};
  for( const PuCase &test_case : sm_cases )
  {
    const auto solution = solve_case(test_case,true,false);
    const auto fitted = fitted_fractions(solution);
    double sum = 0.0, total_variation = 0.0;
    for( size_t i = 0; i < fitted.size(); ++i )
    {
      sum += fitted[i];
      total_variation += std::fabs(fitted[i] - test_case.truth[i]);
      if( (i != 4) || test_case.correlation_in_range )
        BOOST_CHECK_MESSAGE( std::fabs(fitted[i] - test_case.truth[i]) <= max_error[i],
          test_case.id << " isotope index " << i << " fit=" << fitted[i]
                       << " truth=" << test_case.truth[i] );
    }
    BOOST_CHECK_SMALL( sum - 1.0, 1.0e-8 );
    BOOST_CHECK_LE( 0.5*total_variation, 0.05 );

    if( !test_case.correlation_in_range )
    {
      const auto nuclides = pu_nuclides();
      double uncorrected_sum = 0.0;
      for( size_t i = 0; i < nuclides.size(); ++i )
      {
        const auto result = solution.mass_enrichment_result(nuclides[i],0);
        BOOST_REQUIRE( result.pu242_correlation_extrapolated );
        BOOST_REQUIRE( result.uncorrected_fraction.has_value() );
        BOOST_REQUIRE( std::isfinite(*result.uncorrected_fraction) );
        BOOST_CHECK_GE( *result.uncorrected_fraction,0.0 );
        BOOST_CHECK_LE( *result.uncorrected_fraction,1.0 );
        uncorrected_sum += *result.uncorrected_fraction;
        if( i < pu93_direct_truth.size() )
          BOOST_CHECK_MESSAGE(
              std::fabs(*result.uncorrected_fraction-pu93_direct_truth[i]) <= max_error[i],
              test_case.id << " pre-correlation isotope index " << i
                           << " fit=" << *result.uncorrected_fraction
                           << " truth=" << pu93_direct_truth[i] );
        else
          BOOST_CHECK_SMALL( *result.uncorrected_fraction,1.0e-14 );
      }
      BOOST_CHECK_SMALL( uncorrected_sum-1.0,1.0e-10 );
    }
  }
}


// CONV-01/03: a literal one-source Pu problem has no composition split, Pu-242 correlation, or
// semantic activity-rescue variants to hide behind.  It must remain a valid deterministic solve
// and report the sole modeled Pu isotope as a physical unit fraction.
BOOST_AUTO_TEST_CASE( single_pu_source_is_deterministic )
{
  initialize_paths();
  const PuCase &test_case = sm_cases[1];
  const SandiaDecay::Nuclide * const pu240
      = DecayDataBaseServer::database()->nuclide("Pu240");
  BOOST_REQUIRE( pu240 );

  RelActCalcAuto::Options options = preset_options();
  options.auto_profile_weak_mass_fractions = false;
  options.auto_simplify_model = false;
  options.energy_cal_type = RelActCalcAuto::EnergyCalFitType::NoFit;
  options.fwhm_form = RelActCalcAuto::FwhmForm::Polynomial_2;
  options.fwhm_estimation_method
      = RelActCalcAuto::FwhmEstimationMethod::FixedToAllPeaksInSpectrum;
  options.skew_type = PeakDef::SkewType::NoSkew;
  options.additional_br_uncert = 0.0;
  for( RelActCalcAuto::RoiRange &roi : options.rois )
    roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;

  BOOST_REQUIRE_EQUAL( options.rel_eff_curves.size(),1u );
  RelActCalcAuto::RelEffCurveInput &curve = options.rel_eff_curves.front();
  curve.pu242_correlation_method = RelActCalc::PuCorrMethod::NotApplicable;
  curve.mass_fraction_constraints.clear();
  curve.act_ratio_constraints.clear();
  curve.nuclides.erase( std::remove_if(curve.nuclides.begin(),curve.nuclides.end(),
      [pu240]( const RelActCalcAuto::NucInputInfo &input ) {
        return RelActCalcAuto::nuclide(input.source) != pu240;
      }),curve.nuclides.end() );
  BOOST_REQUIRE_EQUAL( curve.nuclides.size(),1u );
  curve.nuclides.front().age = effective_age_years(test_case)*PhysicalUnits::year;
  curve.nuclides.front().fit_age = false;
  curve.nuclides.front().fit_age_min.reset();
  curve.nuclides.front().fit_age_max.reset();
  curve.nuclides.front().force_profile_mass_fraction = false;

  const auto foreground = load_measurement(test_case.spectrum);
  const auto background = load_measurement(test_case.background);
  RelActCalcAuto::RelActAutoSolution solutions[2];
  for( size_t repeat = 0; repeat < 2; ++repeat )
  {
    BOOST_REQUIRE_NO_THROW( solutions[repeat] = RelActCalcAuto::solve(
        options,foreground,background,nullptr,{},
        PeakFitUtils::CoarseResolutionType::High,nullptr) );
    BOOST_REQUIRE_MESSAGE(
        RelActCalcAuto::RelActAutoSolution::is_usable_status(solutions[repeat].m_status),
        "Single-Pu repeat " << repeat << " failed: " << solutions[repeat].m_error_message );
    BOOST_REQUIRE_EQUAL( solutions[repeat].m_rel_activities.size(),1u );
    BOOST_REQUIRE_EQUAL( solutions[repeat].m_rel_activities.front().size(),1u );
    BOOST_CHECK( RelActCalcAuto::nuclide(
        solutions[repeat].m_rel_activities.front().front().source) == pu240 );
    BOOST_CHECK_GT( solutions[repeat].m_rel_activities.front().front().rel_activity,0.0 );
    const auto fraction = solutions[repeat].mass_enrichment_result(pu240,0);
    BOOST_CHECK_SMALL( fraction.fraction-1.0,1.0e-12 );
    BOOST_CHECK( !fraction.pu242_correlation_extrapolated );
  }

  BOOST_REQUIRE_NE( solutions[0].m_frozen_gamma_membership_hash,UINT64_C(0) );
  BOOST_REQUIRE_NE( solutions[0].m_frozen_layout_hash,UINT64_C(0) );
  BOOST_CHECK_EQUAL( solutions[0].m_frozen_gamma_membership_hash,
                     solutions[1].m_frozen_gamma_membership_hash );
  BOOST_CHECK_EQUAL( solutions[0].m_frozen_layout_hash,solutions[1].m_frozen_layout_hash );
  BOOST_CHECK_EQUAL( solutions[0].m_frozen_model_policy_hash,
                     solutions[1].m_frozen_model_policy_hash );
  BOOST_CHECK_SMALL( solutions[0].m_chi2-solutions[1].m_chi2,
                     1.0e-9*(std::max)(1.0,std::fabs(solutions[0].m_chi2)) );
  BOOST_CHECK_SMALL( solutions[0].m_chi2_data-solutions[1].m_chi2_data,
                     1.0e-9*(std::max)(1.0,std::fabs(solutions[0].m_chi2_data)) );
}


BOOST_AUTO_TEST_CASE( free_age_grid_reports_bounded_profile_uncertainty )
{
  initialize_paths();
  const auto nuclides = pu_nuclides();
  const char * const filter_value = std::getenv("INTERSPEC_REL_ACT_PU_CASE");
  const string case_filter = filter_value ? filter_value : "";
  const PuCase *filtered_case = nullptr;
  bool checked_report_contract = false;
  bool checked_generated_pu242_profile = false;
  size_t truth_coverage_checks = 0;
  size_t expected_truth_coverage_checks = 0;
  size_t cases_checked = 0;
  for( const PuCase &test_case : sm_cases )
  {
    if( !case_filter.empty() && (case_filter != test_case.id) )
      continue;
    filtered_case = &test_case;
    ++cases_checked;
    if( test_case.correlation_in_range )
      expected_truth_coverage_checks += test_case.truth.size();

    const auto solution = solve_case(test_case,false,true);
    if( !case_filter.empty() )
    {
      BOOST_TEST_MESSAGE( test_case.id << " selected solution: status="
          << static_cast<int>(solution.m_status) << ", full objective="
          << std::setprecision(17) << solution.m_chi2
          << ", data chi2=" << solution.m_chi2_data << ", cov scale=" << solution.m_cov_scale
          << ", kappa=" << solution.m_jacobian_condition_number
          << ", rank-deficient directions=" << solution.m_num_rank_deficient_dirs );
      for( const string &warning : solution.m_warnings )
        BOOST_TEST_MESSAGE( test_case.id << " solution warning: " << warning );
    }
    check_free_age_diagnostic(solution,test_case.id);
    const auto fitted = fitted_fractions(solution);
    double sum = 0.0;
    for( const double fraction : fitted )
      sum += fraction;
    BOOST_CHECK_SMALL( sum - 1.0, 1.0e-8 );

    // Include correlation-generated Pu-242: it is deliberately absent from NucInputInfo, but a
    // weak structured covariance must still make it an automatic profile target.
    double uncorrected_sum = 0.0;
    bool have_runtime_extrapolation_flag = false;
    bool runtime_extrapolated = false;
    for( size_t i = 0; i < nuclides.size(); ++i )
    {
      const auto result = solution.mass_enrichment_result(nuclides[i],0);
      // Pu-242 correlation renormalizes every reported Pu component.  If the correlation itself is
      // outside its validation domain, all corrected intervals are conditional on an uncalibrated
      // model extrapolation; applying a nominal 95% coverage gate to only Pu-238..241 would be as
      // misleading as applying it to generated Pu-242.  In-range Pu70/Pu61 retain all five gates.
      // Whether the certificate lies inside the correlations validated domain determines if a
      // corrected-quantity coverage assertion is scientifically calibrated.  The runtime flag is
      // intentionally separate: a free-age, non-identifiable central basin may itself leave that
      // domain even for an in-range reference material, and must then be reported conditionally.
      const bool require_truth_coverage = test_case.correlation_in_range;
      if( have_runtime_extrapolation_flag )
        BOOST_CHECK_EQUAL( result.pu242_correlation_extrapolated,runtime_extrapolated );
      else
      {
        runtime_extrapolated = result.pu242_correlation_extrapolated;
        have_runtime_extrapolation_flag = true;
      }
      if( result.pu242_correlation_extrapolated )
      {
        BOOST_CHECK( SpecUtils::icontains(result.diagnostic,"conditional") );
        BOOST_CHECK( SpecUtils::icontains(result.diagnostic,"systematic") );
      }
      BOOST_REQUIRE( result.uncorrected_fraction.has_value() );
      BOOST_REQUIRE( std::isfinite(*result.uncorrected_fraction) );
      BOOST_CHECK_GE( *result.uncorrected_fraction,0.0 );
      BOOST_CHECK_LE( *result.uncorrected_fraction,1.0 );
      uncorrected_sum += *result.uncorrected_fraction;
      if( result.covariance_quality
          == RelActCalcAuto::RelActAutoSolution::MassFractionCovarianceQuality::Usable )
      {
        // Automatic profiling is intentionally limited to weak quantities.  A healthy quantity's
        // structured local covariance still defines its conventional two-sided Gaussian 95% band;
        // clip that band to the physical fraction range and enforce the same truth gate.
        BOOST_REQUIRE_MESSAGE( result.covariance_one_sigma.has_value(),
                               test_case.id << " " << nuclides[i]->symbol
                                            << " usable covariance has no local sigma" );
        const double sigma = *result.covariance_one_sigma;
        BOOST_REQUIRE( std::isfinite(sigma) );
        BOOST_REQUIRE_GE( sigma,0.0 );
        constexpr double gaussian_95_quantile = 1.959963984540054;
        const double lower = (std::max)(0.0,result.fraction-gaussian_95_quantile*sigma);
        const double upper = (std::min)(1.0,result.fraction+gaussian_95_quantile*sigma);
        if( require_truth_coverage )
        {
          BOOST_CHECK_LE( lower,test_case.truth[i] );
          BOOST_CHECK_GE( upper,test_case.truth[i] );
          ++truth_coverage_checks;
        }
        continue;
      }

      BOOST_REQUIRE_MESSAGE( result.profile.has_value(),
                             test_case.id << " " << nuclides[i]->symbol
                                          << " weak covariance was not profiled" );
      const auto status = result.profile->status;
      const bool profile_classified
          = status == RelActCalcAuto::RelActAutoSolution::MassFractionProfileStatus::Complete
          || status == RelActCalcAuto::RelActAutoSolution::MassFractionProfileStatus::BoundaryLimited
          || status == RelActCalcAuto::RelActAutoSolution::MassFractionProfileStatus::NonIdentifiable;
      if( !case_filter.empty() && !profile_classified )
      {
        BOOST_TEST_MESSAGE( test_case.id << " " << nuclides[i]->symbol
            << " failed profile: status=" << static_cast<int>(status)
            << ", reason=" << static_cast<int>(result.profile->reason)
            << ", fits=" << result.profile->num_fits << ", fraction="
            << std::setprecision(17) << result.fraction
            << ", covariance quality=" << static_cast<int>(result.covariance_quality)
            << ", local sigma="
            << (result.covariance_one_sigma ? std::to_string(*result.covariance_one_sigma)
                                            : string("unavailable"))
            << ", full objective=" << solution.m_chi2
            << ", data chi2=" << solution.m_chi2_data
            << ", diagnostic=" << result.diagnostic
            << ", profile message=" << result.profile->message );
        for( const auto &sample : result.profile->sampled_delta_chi2 )
          BOOST_TEST_MESSAGE( test_case.id << " " << nuclides[i]->symbol
              << " failed-profile sample: fraction=" << std::setprecision(17) << sample.first
              << ", delta chi2=" << sample.second
              << ", implied objective=" << solution.m_chi2 + sample.second );
      }
      BOOST_CHECK_MESSAGE( profile_classified,
          test_case.id << " " << nuclides[i]->symbol << " profile status="
                       << static_cast<int>(status) << ": " << result.profile->message );
      const auto interval = std::find_if(result.profile->intervals.begin(),
                                          result.profile->intervals.end(),
        []( const auto &candidate ){
          return std::fabs(candidate.confidence_level - 0.95) < 0.01;
        });
      BOOST_REQUIRE_MESSAGE( interval != result.profile->intervals.end(),
          test_case.id << " " << nuclides[i]->symbol
                       << " profile supplied no 95% interval: " << result.profile->message );
      BOOST_CHECK_GE( interval->lower, 0.0 );
      BOOST_CHECK_LE( interval->upper, 1.0 );
      if( require_truth_coverage )
      {
        BOOST_CHECK_LE( interval->lower, test_case.truth[i] );
        BOOST_CHECK_GE( interval->upper, test_case.truth[i] );
        ++truth_coverage_checks;
      }
      BOOST_CHECK_LE( result.profile->num_fits, 32u );
      if( i == 4 )
      {
        BOOST_CHECK( result.profile->reason
            == RelActCalcAuto::RelActAutoSolution::MassFractionProfileReason::AutomaticWeak );
        checked_generated_pu242_profile = true;
      }
    }
    BOOST_CHECK_SMALL( uncorrected_sum-1.0,1.0e-10 );
    if( runtime_extrapolated )
    {
      const bool warned = std::any_of( solution.m_warnings.begin(),solution.m_warnings.end(),
        []( const string &warning ) {
          return SpecUtils::icontains(warning,"outside range")
              && SpecUtils::icontains(warning,"Pu242 correlation");
        } );
      BOOST_CHECK_MESSAGE( warned,test_case.id
          << " extrapolated the Pu-242 correlation without a solution-level warning" );
    }

    const auto pu242 = solution.mass_enrichment_result(nuclides[4],0);
    if( !test_case.correlation_in_range )
      BOOST_CHECK( pu242.pu242_correlation_extrapolated );
    BOOST_CHECK( std::isfinite(pu242.fraction) );

    // Reuse this already-computed fit for the hard-coded text/HTML formatter regression.  Pu-239
    // is present in every case and has a nonzero local covariance; run it once for the grid.
    if( !checked_report_contract )
    {
      check_pu_report_uncertainty_contract(solution,nuclides[1]);
      checked_report_contract = true;
    }
  }
  BOOST_REQUIRE_MESSAGE( cases_checked > 0,
                         "INTERSPEC_REL_ACT_PU_CASE='" << case_filter
                         << "' did not match Pu93, Pu70, or Pu61" );
  BOOST_CHECK( checked_report_contract );
  BOOST_CHECK_EQUAL( truth_coverage_checks,expected_truth_coverage_checks );
  BOOST_CHECK_MESSAGE( checked_generated_pu242_profile,
                       "No weak correlation-generated Pu-242 quantity was automatically profiled" );

  // Repeat one representative free-age case without profile scans.  This keeps the determinism
  // oracle focused on the selected physical basin/frozen objective and avoids repeating any of the
  // three expensive profile grids above.  The normal grid uses in-range Pu70; a diagnostic filter
  // repeats the selected material so the entire invocation remains case-specific.
  const PuCase &repeat_case = case_filter.empty() ? sm_cases[1] : *filtered_case;
  const auto repeat_a = solve_case(repeat_case,false,false);
  const auto repeat_b = solve_case(repeat_case,false,false);
  const FittedAgeDiagnostic age_a = check_free_age_diagnostic(repeat_a,repeat_case.id);
  const FittedAgeDiagnostic age_b = check_free_age_diagnostic(repeat_b,repeat_case.id);
  BOOST_REQUIRE_NE( repeat_a.m_frozen_gamma_membership_hash,UINT64_C(0) );
  BOOST_REQUIRE_NE( repeat_a.m_frozen_layout_hash,UINT64_C(0) );
  BOOST_CHECK_EQUAL( repeat_a.m_frozen_gamma_membership_hash,
                     repeat_b.m_frozen_gamma_membership_hash );
  BOOST_CHECK_EQUAL( repeat_a.m_frozen_layout_hash,repeat_b.m_frozen_layout_hash );
  const double objective_scale = (std::max)(1.0,std::fabs(repeat_a.m_chi2));
  BOOST_CHECK_SMALL( repeat_a.m_chi2-repeat_b.m_chi2,1.0e-9*objective_scale );
  BOOST_CHECK_SMALL( age_a.age-age_b.age,
                     1.0e-9*(std::max)(1.0,std::fabs(age_a.age)) );
  BOOST_CHECK_SMALL( age_a.one_sigma-age_b.one_sigma,
                     1.0e-8*(std::max)(1.0,std::fabs(age_a.one_sigma)) );
  const auto fractions_a = fitted_fractions(repeat_a);
  const auto fractions_b = fitted_fractions(repeat_b);
  for( size_t isotope = 0; isotope < fractions_a.size(); ++isotope )
    BOOST_CHECK_SMALL( fractions_a[isotope]-fractions_b[isotope],1.0e-8 );
}


BOOST_AUTO_TEST_CASE( pu242_augmented_lagrangian_targets_reported_fraction_only )
{
  initialize_paths();
  const PuCase &test_case = sm_cases[1]; // Pu70: inside the Pu-242 correlation validation range.
  const auto nuclides = pu_nuclides();
  const SandiaDecay::Nuclide * const pu242 = nuclides[4];
  BOOST_REQUIRE( pu242 );

  RelActCalcAuto::Options options = preset_options();
  options.auto_profile_weak_mass_fractions = false;
  for( RelActCalcAuto::RoiRange &roi : options.rois )
    roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
  for( RelActCalcAuto::NucInputInfo &input : options.rel_eff_curves.at(0).nuclides )
  {
    input.age = effective_age_years(test_case) * PhysicalUnits::year;
    input.fit_age = false;
    input.fit_age_min.reset();
    input.fit_age_max.reset();
    input.force_profile_mass_fraction = false;
  }

  const auto foreground = load_measurement(test_case.spectrum);
  const auto background = load_measurement(test_case.background);
  RelActCalcAuto::RelActAutoSolution baseline;
  BOOST_REQUIRE_NO_THROW( baseline = RelActCalcAuto::solve(
      options,foreground,background,nullptr,{},
      PeakFitUtils::CoarseResolutionType::High,nullptr) );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status(baseline.m_status),
                         baseline.m_error_message );

  const auto nominal = baseline.mass_enrichment_result(pu242,0);
  BOOST_REQUIRE( std::isfinite(nominal.fraction) );
  BOOST_REQUIRE( nominal.uncorrected_fraction.has_value() );
  BOOST_CHECK_SMALL( *nominal.uncorrected_fraction,1.0e-14 );
  const double requested = (std::min)(0.95,nominal.fraction + 0.002);
  BOOST_REQUIRE_GT( requested,nominal.fraction );

  vector<shared_ptr<const PeakDef>> frozen_peaks;
  frozen_peaks.reserve( baseline.m_fit_peaks_in_spectrums_cal.size() );
  for( const PeakDef &peak : baseline.m_fit_peaks_in_spectrums_cal )
    frozen_peaks.push_back( make_shared<const PeakDef>(peak) );

  // This is the production initial penalty, chosen to make most profile points one conditional fit.
  constexpr double penalty = 1.0e12;
  RelActCalcAuto::Options::ProfileOnlyMassFractionConstraint constraint;
  constraint.rel_eff_curve_index = 0;
  constraint.nuclide = pu242; // Correlation-generated: deliberately absent from fitted inputs.
  constraint.reported_fraction = requested;
  constraint.lagrange_multiplier = 0.0;
  constraint.penalty = penalty;
  options.profile_only_mass_fraction_constraint = constraint;

  RelActCalcAuto::RelActAutoSolution conditional;
  BOOST_REQUIRE_NO_THROW( conditional = RelActCalcAuto::solve(
      options,foreground,background,nullptr,frozen_peaks,
      PeakFitUtils::CoarseResolutionType::High,nullptr) );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status(conditional.m_status),
                         conditional.m_error_message );

  const auto fitted = conditional.mass_enrichment_result(pu242,0);
  BOOST_REQUIRE( std::isfinite(fitted.fraction) );
  BOOST_REQUIRE( std::isfinite(conditional.m_profile_constraint_violation) );
  BOOST_CHECK_SMALL( fitted.fraction-requested,3.0e-6 );
  BOOST_CHECK_SMALL( conditional.m_profile_constraint_violation
                       - (fitted.fraction-requested),1.0e-10 );

  // Ceres minimizes 0.5*sum(r^2); the production AL row is therefore
  // sqrt(mu)*(h+lambda/mu).  Its diagnostic contribution is excluded from physical m_chi2.
  const double expected_penalty_chi2
      = penalty*conditional.m_profile_constraint_violation
                *conditional.m_profile_constraint_violation;
  BOOST_CHECK_SMALL( conditional.m_profile_constraint_penalty_chi2
                       - expected_penalty_chi2,
                     1.0e-8*(1.0+expected_penalty_chi2) );
  BOOST_CHECK_SMALL( conditional.m_optimizer_chi2
                       - conditional.m_chi2
                       - conditional.m_profile_constraint_penalty_chi2,
                     1.0e-8*(1.0+conditional.m_optimizer_chi2) );
}
