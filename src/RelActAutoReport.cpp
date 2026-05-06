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
#include <chrono>
#include <ctime>
#include <random>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "Wt/WDateTime.h"
#include "Wt/WApplication.h"
#include "Wt/WLocalDateTime.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/DateTime.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/D3SpectrumExport.h"
#include "SpecUtils/EnergyCalibration.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/AppUtils.h"
#include "InterSpec/PeakDef.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/RelEffChart.h"
#include "InterSpec/BatchInfoLog.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/RelActAutoReport.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/RelActCalc_imp.hpp"

using namespace std;

namespace
{
#if( defined(_WIN32) )
  const char ns_path_sep = '\\';
#else
  const char ns_path_sep = '/';
#endif

  /** Resolve `InterSpec_resources/` whether or not we are running inside a Wt app. */
  string interspec_resources_dir()
  {
    Wt::WApplication * const app = Wt::WApplication::instance();
    if( app )
      return SpecUtils::append_path( app->docRoot(), "InterSpec_resources" );

    const string static_data_dir = InterSpec::staticDataDirectory().empty()
                                   ? string("./data") : InterSpec::staticDataDirectory();
    const string app_root = SpecUtils::append_path( static_data_dir, ".." );
    return SpecUtils::append_path( app_root, "InterSpec_resources" );
  }


  /** Format the m_status enum into a short human-readable string for HTML/text use. */
  string status_fail_reason( const RelActCalcAuto::RelActAutoSolution::Status status )
  {
    switch( status )
    {
      case RelActCalcAuto::RelActAutoSolution::Status::Success:               return "Success";
      case RelActCalcAuto::RelActAutoSolution::Status::NotInitiated:          return "Not initiated";
      case RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem:  return "Failed to setup problem";
      case RelActCalcAuto::RelActAutoSolution::Status::FailToSolveProblem:    return "Failed to solve problem";
      case RelActCalcAuto::RelActAutoSolution::Status::UserCanceled:          return "User canceled";
    }
    return "Unknown";
  }


  /** Energy-calibration parameter index -> physical-units conversion (matches legacy report). */
  double phys_energy_cal_value( const RelActCalcAuto::RelActAutoSolution &sol, const size_t i )
  {
    double value = (sol.m_energy_cal_adjustments[i]/RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset - 1.0);

    const auto cal = sol.m_spectrum ? sol.m_spectrum->energy_calibration() : nullptr;
    const size_t num_chan = (cal ? cal->num_channels() : 1);

    if( i == 0 )
      value *= RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple;
    else if( i == 1 )
      value *= RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple / num_chan;
    else if( i == 2 )
      value *= RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple / (num_chan * num_chan);

    return value;
  }


  /** Build the same `vector<RelEffChart::ReCurveInfo>` the legacy report builds, ready for
   `RelEffChart::jsonForData`.  Logic mirrors `print_html_report`'s rel-eff section.
   */
  vector<RelEffChart::ReCurveInfo> make_rel_eff_info_sets(
      const RelActCalcAuto::RelActAutoSolution &sol )
  {
    vector<RelEffChart::ReCurveInfo> rel_eff_info_sets;

    for( size_t rel_eff_index = 0; rel_eff_index < sol.m_rel_eff_forms.size(); ++rel_eff_index )
    {
      if( rel_eff_index >= sol.m_options.rel_eff_curves.size()
          || rel_eff_index >= sol.m_rel_eff_coefficients.size()
          || rel_eff_index >= sol.m_rel_activities.size() )
      {
        break;
      }

      const RelActCalcAuto::RelEffCurveInput &rel_eff = sol.m_options.rel_eff_curves[rel_eff_index];
      RelEffChart::ReCurveInfo info;
      info.live_time = sol.m_spectrum ? sol.m_spectrum->live_time() : 1.0;

      if( rel_eff_index < sol.m_obs_eff_for_each_curve.size() )
      {
        for( const RelActCalcAuto::RelActAutoSolution::ObsEff &obs_eff
             : sol.m_obs_eff_for_each_curve[rel_eff_index] )
        {
          if( obs_eff.observed_efficiency > 0.0
              && obs_eff.num_sigma_significance > 2.5
              && obs_eff.fraction_roi_counts > 0.05
              && obs_eff.within_roi )
          {
            info.obs_eff_data.push_back( obs_eff );
          }
        }
      }

      info.rel_acts = sol.m_rel_activities[rel_eff_index];
      info.re_curve_name = Wt::WString::fromUTF8( rel_eff.name );

      try
      {
        info.re_curve_eqn_txt = Wt::WString::fromUTF8( "y = " + sol.rel_eff_txt( false, rel_eff_index ) );
        info.js_rel_eff_eqn = sol.rel_eff_eqn_js_function( rel_eff_index );
      }catch( const std::exception &e )
      {
        cerr << "RelActAutoReport: error preparing rel-eff curve " << rel_eff_index << ": " << e.what() << endl;
      }

      rel_eff_info_sets.push_back( info );
    }

    return rel_eff_info_sets;
  }
}//namespace


namespace RelActAutoReport
{

// ---- Inja callback helpers ------------------------------------------------

void html_sanitize( string &val )
{
  // Same set of escapes as the legacy print_html_report's local lambda.
  SpecUtils::ireplace_all( val, "&", "&amp;" );
  SpecUtils::ireplace_all( val, "<", "&lt;" );
  SpecUtils::ireplace_all( val, ">", "&gt;" );
  SpecUtils::ireplace_all( val, "'", "&#39;" );
  SpecUtils::ireplace_all( val, "\"", "&quot;" );
}

string format_scientific( double value, int precision )
{
  ostringstream ss;
  ss << std::scientific << std::setprecision( precision ) << value;
  return ss.str();
}

string format_percentage( double value, int precision )
{
  ostringstream ss;
  ss << std::fixed << std::setprecision( precision ) << (value * 100.0) << "%";
  return ss.str();
}

string format_duration( int64_t microseconds )
{
  const double seconds = microseconds * 1.0e-6;
  return PhysicalUnits::printToBestTimeUnits( seconds );
}

string pct_callback( vector<const nlohmann::json *> &args )
{
  try
  {
    if( args.empty() || args[0]->is_null() )
      return "--";
    if( args[0]->is_string() )
      return args[0]->get<string>();
    if( !args[0]->is_number() )
      throw runtime_error( "pct: argument is not a number" );

    const double value = args[0]->get<double>();
    int precision = 4;
    if( args.size() > 1 && args[1]->is_number() )
      precision = std::max( 1, args[1]->get<int>() );

    return SpecUtils::printCompact( 100.0 * value, static_cast<size_t>(precision) );
  }catch( const std::exception &e )
  {
    cerr << "Error in 'pct': " << e.what() << endl;
    return string("ErrorPct{") + e.what() + "}";
  }
}

string safe_html_callback( vector<const nlohmann::json *> &args )
{
  if( args.empty() || args[0]->is_null() )
    return "";
  string val;
  if( args[0]->is_string() )
    val = args[0]->get<string>();
  else
    val = args[0]->dump();
  html_sanitize( val );
  return val;
}


// ---- Template directory + Inja environment helpers -----------------------

string default_template_dir()
{
  const string static_data_dir = InterSpec::staticDataDirectory().empty()
                                 ? string("./data") : InterSpec::staticDataDirectory();
  const string app_root = SpecUtils::append_path( static_data_dir, ".." );
  const string docroot  = SpecUtils::append_path( app_root, "InterSpec_resources" );
  const string static_txt = SpecUtils::append_path( docroot, "static_text" );
  return SpecUtils::append_path( static_txt, "IsotopicsByNuclidesReportTmplts" ) + ns_path_sep;
}

string template_include_dir( const string &user_path )
{
  if( user_path.empty() || SpecUtils::iequals_ascii( user_path, "none" ) )
    return "";
  if( SpecUtils::iequals_ascii( user_path, "default" ) )
    return default_template_dir();

  string tmplt_dir = user_path;
  if( !tmplt_dir.empty() && tmplt_dir.back() != ns_path_sep )
    tmplt_dir += ns_path_sep;
  return tmplt_dir;
}


string writable_template_dir()
{
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP || BUILD_AS_UNIT_TEST_SUITE )
  // InterSpec::writableDataDirectory() throws if not initialised (e.g. in unit
  // tests), so guard with a try.  Also return "" if the directory does not yet
  // exist, since enumeration would just return nothing anyway.
  string base;
  try
  {
    base = InterSpec::writableDataDirectory();
  }catch( std::exception & )
  {
    return "";
  }
  if( base.empty() )
    return "";
  const string tmplt_dir = SpecUtils::append_path( base, "IsotopicsByNuclidesReportTmplts" );
  if( !SpecUtils::is_directory(tmplt_dir) )
    return "";
  return tmplt_dir + ns_path_sep;
#else
  // Build configurations without a writable data directory (e.g. peak_fit_improve)
  //  cannot supply user-provided templates.
  return "";
#endif
}


vector<pair<string,string>> load_spectrum_chart_js_and_css()
{
  // Reuse BatchInfoLog's identical helper.
  return BatchInfoLog::load_spectrum_chart_js_and_css();
}

vector<pair<string,string>> load_rel_eff_plot_js_and_css()
{
  vector<pair<string,string>> answer;

#if( SpecUtils_ENABLE_D3_CHART )
  const string docroot = interspec_resources_dir();
  try
  {
    string plot_js  = AppUtils::file_contents( SpecUtils::append_path( docroot, "RelEffPlot.js" ) );
    string plot_css = AppUtils::file_contents( SpecUtils::append_path( docroot, "RelEffPlot.css" ) );
    answer.emplace_back( "RelEffPlot_JS", std::move(plot_js) );
    answer.emplace_back( "RelEffPlot_CSS", std::move(plot_css) );
  }catch( const std::exception &e )
  {
    cerr << "RelActAutoReport: failed to load RelEffPlot assets: " << e.what() << endl;
  }
#endif

  return answer;
}


inja::Environment get_default_inja_env( const string &tmplt_dir )
{
  if( !tmplt_dir.empty() && !SpecUtils::is_directory(tmplt_dir) )
    throw runtime_error( string("Template include directory, '") + tmplt_dir
                         + "', does not exist." );

  inja::Environment env{ tmplt_dir };
  env.set_trim_blocks( true );

#if( BUILD_FOR_WEB_DEPLOYMENT )
  env.set_search_included_templates_in_files( false );
#else
  env.set_search_included_templates_in_files( !tmplt_dir.empty() );
#endif

  // BatchInfoLog callbacks (delegated to its tested implementations).
  env.add_callback( "printFixed",   2, &BatchInfoLog::printFixed );
  env.add_callback( "printCompact", 2, &BatchInfoLog::printCompact );

  // RelActAutoReport-specific callbacks.
  env.add_callback( "pct",       1, &pct_callback );
  env.add_callback( "pct",       2, &pct_callback );
  env.add_callback( "safe_html", 1, &safe_html_callback );

  // Backward-compatibility callbacks (kept simple - prefer printFixed/printCompact instead).
  env.add_callback( "scientific", 2, []( vector<const nlohmann::json *> &args ) -> string {
    if( args.size() < 2 || !args[0]->is_number() || !args[1]->is_number() )
      return "";
    return format_scientific( args[0]->get<double>(), args[1]->get<int>() );
  } );

  env.add_callback( "format", 2, []( vector<const nlohmann::json *> &args ) -> string {
    if( args.size() < 2 || !args[0]->is_string() || !args[1]->is_number() )
      return "";
    char buffer[256] = { '\0' };
    snprintf( buffer, sizeof(buffer), args[0]->get<string>().c_str(), args[1]->get<double>() );
    return string(buffer);
  } );

  env.add_callback( "default", 2, []( vector<const nlohmann::json *> &args ) -> string {
    if( args.size() < 2 ) return "";
    if( args[0]->is_null() || args[0]->empty() )
      return args[1]->is_string() ? args[1]->get<string>() : args[1]->dump();
    return args[0]->is_string() ? args[0]->get<string>() : args[0]->dump();
  } );

  // Pre-register the default templates as named includes.
  // Mirrors BatchInfoLog::get_default_inja_env's pattern.
  try
  {
    inja::Environment sub_env;
    sub_env.add_callback( "printFixed",   2, &BatchInfoLog::printFixed );
    sub_env.add_callback( "printCompact", 2, &BatchInfoLog::printCompact );
    sub_env.add_callback( "pct",          1, &pct_callback );
    sub_env.add_callback( "pct",          2, &pct_callback );
    sub_env.add_callback( "safe_html",    1, &safe_html_callback );

    const string default_dir = default_template_dir();
    const string html_path = SpecUtils::append_path( default_dir, "std_rel_eff_summary.tmplt.html" );
    const string txt_path  = SpecUtils::append_path( default_dir, "std_rel_eff_summary.tmplt.txt" );
    const string multi_html_path = SpecUtils::append_path( default_dir, "std_multi_file_summary.tmplt.html" );

    if( SpecUtils::is_file(html_path) )
    {
      inja::Template html_tmplt = sub_env.parse_template( html_path );
      env.include_template( "default-rel-act-auto-html-results", html_tmplt );
    }
    if( SpecUtils::is_file(txt_path) )
    {
      inja::Template txt_tmplt = sub_env.parse_template( txt_path );
      env.include_template( "default-rel-act-auto-txt-results", txt_tmplt );
    }
    if( SpecUtils::is_file(multi_html_path) )
    {
      inja::Template multi_tmplt = sub_env.parse_template( multi_html_path );
      env.include_template( "default-rel-act-auto-html-multi-file-summary", multi_tmplt );
    }
  }catch( const std::exception &e )
  {
    throw runtime_error( "Error loading default RelActAuto report template: " + string(e.what()) );
  }

  return env;
}


// ---- solution_to_json ----------------------------------------------------

nlohmann::json solution_to_json( const RelActCalcAuto::RelActAutoSolution &sol )
{
  using nlohmann::json;
  json data;

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const bool success = (sol.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success);
  const float live_time = sol.m_spectrum ? sol.m_spectrum->live_time() : 1.0f;
  const bool have_multiple_rel_eff = (sol.m_rel_eff_forms.size() > 1);

  // ---- Status / metadata ----
  data["status"] = {
    { "success",       success },
    { "status_code",   static_cast<int>(sol.m_status) },
    { "fail_reason",   status_fail_reason( sol.m_status ) },
    { "error_message", sol.m_error_message }
  };

  data["spectrum_title"]        = sol.m_options.spectrum_title;
  data["have_multiple_rel_eff"] = have_multiple_rel_eff;
  data["live_time_s"]           = live_time;

  // ---- Goodness-of-fit ----
  data["chi2"]            = sol.m_chi2;
  data["dof"]             = static_cast<int64_t>(sol.m_dof);
  data["chi2_per_dof"]    = sol.m_dof > 0 ? sol.m_chi2 / static_cast<double>(sol.m_dof) : 0.0;
  {
    char buf[64] = { '\0' };
    snprintf( buf, sizeof(buf), "%.6G", sol.m_chi2 );
    data["chi2_str"] = string(buf);
    snprintf( buf, sizeof(buf), "%.6G", data["chi2_per_dof"].get<double>() );
    data["chi2_per_dof_str"] = string(buf);
  }

  // ---- Warnings ----
  data["warnings"] = json::array();
  data["warnings_html"] = json::array();
  for( const string &w : sol.m_warnings )
  {
    data["warnings"].push_back( w );
    string sanitized = w;
    html_sanitize( sanitized );
    data["warnings_html"].push_back( sanitized );
  }

  // ---- Rel-eff curves ----
  data["rel_eff_curves"] = json::array();
  for( size_t i = 0; i < sol.m_rel_eff_forms.size(); ++i )
  {
    json curve;
    curve["index"] = static_cast<int64_t>(i);
    curve["name"] = (i < sol.m_options.rel_eff_curves.size())
                    ? sol.m_options.rel_eff_curves[i].name : string();
    curve["rel_eff_eqn_type"] = (i < sol.m_options.rel_eff_curves.size())
                                ? RelActCalc::to_str( sol.m_options.rel_eff_curves[i].rel_eff_eqn_type )
                                : "";
    if( i < sol.m_options.rel_eff_curves.size()
        && sol.m_options.rel_eff_curves[i].rel_eff_eqn_type != RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      curve["rel_eff_eqn_order"] = static_cast<int64_t>(sol.m_options.rel_eff_curves[i].rel_eff_eqn_order);
    }

    if( i < sol.m_rel_eff_coefficients.size() )
      curve["coefficients"] = sol.m_rel_eff_coefficients[i];

    try
    {
      curve["equation_text"] = sol.rel_eff_txt( false, i );
      curve["equation_html"] = sol.rel_eff_txt( true, i );
      curve["js_rel_eff_eqn"] = sol.rel_eff_eqn_js_function( i );
    }catch( const std::exception &e )
    {
      curve["equation_error"] = e.what();
    }

    data["rel_eff_curves"].push_back( curve );
  }

  // ---- Rel-eff chart JSON (interactive D3 chart payload) ----
  if( success )
  {
    try
    {
      const vector<RelEffChart::ReCurveInfo> infos = make_rel_eff_info_sets( sol );
      const string rel_eff_chart_json = RelEffChart::jsonForData( infos );
      data["rel_eff_chart_json"] = rel_eff_chart_json;
    }catch( const std::exception &e )
    {
      cerr << "RelActAutoReport: rel_eff_chart_json failed: " << e.what() << endl;
      data["rel_eff_chart_json"] = "[]";
    }
  }else
  {
    data["rel_eff_chart_json"] = "[]";
  }

  // ---- Per-curve relative activities + Pu corrections ----
  // Pu corrections are also nested into each `curve_acts` object (key "pu"), so that templates
  // that just iterate `relative_activities` can render the corresponding Pu table without
  // doing array indexing across two parallel arrays (which Inja v3 does not support cleanly).
  data["relative_activities"] = json::array();
  data["pu_corrections"] = json::array();

  for( size_t curve_idx = 0; curve_idx < sol.m_rel_activities.size(); ++curve_idx )
  {
    const auto &rel_acts = sol.m_rel_activities[curve_idx];
    const RelActCalcAuto::RelEffCurveInput * rel_eff = nullptr;
    if( curve_idx < sol.m_options.rel_eff_curves.size() )
      rel_eff = &sol.m_options.rel_eff_curves[curve_idx];

    json curve_acts;
    curve_acts["curve_index"] = static_cast<int64_t>(curve_idx);

    double sum_rel_mass = 0.0;
    for( const RelActCalcAuto::NuclideRelAct &act : rel_acts )
    {
      const SandiaDecay::Nuclide *nuc = RelActCalcAuto::nuclide(act.source);
      if( nuc )
        sum_rel_mass += act.rel_activity / nuc->activityPerGram();
    }

    curve_acts["nuclides"] = json::array();
    for( const RelActCalcAuto::NuclideRelAct &act : rel_acts )
    {
      json nuc_info;
      nuc_info["name"] = act.name();
      nuc_info["rel_activity"] = act.rel_activity;
      nuc_info["rel_activity_uncertainty"] = act.rel_activity_uncertainty;
      nuc_info["age"] = act.age;
      nuc_info["age_uncertainty"] = act.age_uncertainty;
      nuc_info["age_was_fit"] = act.age_was_fit;

      // Pre-formatted "rel_activity ± uncertainty" matching the legacy raw stream output.
      {
        ostringstream oss;
        oss << act.rel_activity << " ± " << act.rel_activity_uncertainty;
        nuc_info["rel_activity_with_uncert_str"] = oss.str();
      }

      const SandiaDecay::Nuclide *nuc = RelActCalcAuto::nuclide(act.source);
      bool is_pu242_by_corr = false;
      if( nuc && rel_eff
          && nuc->atomicNumber == 94 && nuc->massNumber == 242
          && rel_eff->pu242_correlation_method != RelActCalc::PuCorrMethod::NotApplicable )
      {
        is_pu242_by_corr = true;
      }
      nuc_info["is_pu242_by_corr"] = is_pu242_by_corr;

      if( nuc )
      {
        const double rel_mass = act.rel_activity / nuc->activityPerGram();
        nuc_info["rel_mass"] = rel_mass;
        nuc_info["total_mass_fraction"] = sum_rel_mass > 0 ? (rel_mass / sum_rel_mass) : 0.0;

        try
        {
          const pair<double,std::optional<double>> enr = sol.mass_enrichment_fraction(nuc, curve_idx);
          nuc_info["enrichment"] = enr.first;
          if( enr.second.has_value() )
          {
            nuc_info["has_enrichment_uncert"]   = true;
            nuc_info["enrichment_uncert"]       = *enr.second;
            nuc_info["enrichment_minus_2sigma"] = enr.first - 2.0 * (*enr.second);
            nuc_info["enrichment_plus_2sigma"]  = enr.first + 2.0 * (*enr.second);
          }else
          {
            nuc_info["has_enrichment_uncert"] = false;
          }
        }catch( const std::exception &e )
        {
          nuc_info["enrichment_error"] = e.what();
        }

        try { nuc_info["detector_counts"] = sol.nuclide_counts(act.source, curve_idx); }
        catch( const std::exception &e ) { nuc_info["counts_error"] = e.what(); }

        // Age display (matches legacy formatting).
        const PhysicalUnits::UnitNameValuePair time_units = PhysicalUnits::bestTimeUnit( act.age );
        const double age_units = act.age / time_units.second;
        const double age_uncert_units = act.age_uncertainty / time_units.second;
        ostringstream age_oss;
        if( act.age_was_fit )
          age_oss << PhysicalUnits::printValueWithUncertainty( age_units, age_uncert_units, 5 ) << " " << time_units.first;
        else if( act.age == 0.0 )
          age_oss << "0";
        else
          age_oss << SpecUtils::printCompact( age_units, 5 ) << " " << time_units.first;
        nuc_info["age_str"] = age_oss.str();
      }

      curve_acts["nuclides"].push_back( nuc_info );
    }
    data["relative_activities"].push_back( curve_acts );

    // ---- Pu corrections (per curve) ----
    if( success
        && curve_idx < sol.m_corrected_pu.size()
        && curve_idx < sol.m_uncorrected_pu.size()
        && sol.m_corrected_pu[curve_idx]
        && sol.m_uncorrected_pu[curve_idx] )
    {
      const auto &pu_corr = *sol.m_corrected_pu[curve_idx];
      const auto &pu_unc  = *sol.m_uncorrected_pu[curve_idx];

      json pu_data;
      pu_data["curve_index"] = static_cast<int64_t>(curve_idx);

      // Set of which Pu nuclides actually appear in this curve's activities.
      set<string> pu_nuclide_names;
      set<double> pu_nuclide_ages;
      vector<tuple<const SandiaDecay::Nuclide *,double>> pu_rel_acts;
      double pu239_act = 0.0;
      for( const RelActCalcAuto::NuclideRelAct &act : rel_acts )
      {
        const SandiaDecay::Nuclide *nuc = RelActCalcAuto::nuclide(act.source);
        if( nuc && nuc->atomicNumber == 94 )
        {
          if( nuc->massNumber == 239 )
            pu239_act = act.rel_activity;
          pu_nuclide_names.insert( act.name() );
          pu_nuclide_ages.insert( act.age );
          pu_rel_acts.emplace_back( nuc, act.rel_activity );
        }
      }

      // Build ordered Pu rows (Pu238 -> Pu242).
      pu_data["rows"] = json::array();
      static const struct{ const char *iso; int A; } iso_table[] = {
        { "Pu238", 238 }, { "Pu239", 239 }, { "Pu240", 240 }, { "Pu241", 241 }, { "Pu242", 242 }
      };
      const double uncorr_for[5] = {
        pu_unc.pu238_rel_mass, pu_unc.pu239_rel_mass,
        pu_unc.pu240_rel_mass, pu_unc.pu241_rel_mass,
        0.0  //Pu242 by-correlation: blank in uncorrected column
      };

      for( size_t k = 0; k < 5; ++k )
      {
        if( !db ) break;
        const SandiaDecay::Nuclide *nuc = db->nuclide(iso_table[k].iso);
        if( !nuc )
          continue;
        if( !pu_nuclide_names.count(iso_table[k].iso) && iso_table[k].A != 242 )
          continue;

        double mass_frac = -1.0, mass_frac_uncert = -1.0;
        try
        {
          const pair<double,std::optional<double>> v = sol.mass_enrichment_fraction(nuc, curve_idx);
          mass_frac = v.first;
          if( v.second.has_value() )
            mass_frac_uncert = *v.second;
        }catch( const std::exception & )
        {
          // Fallback to corrected struct if covariance failed (matches legacy fallback).
          if( iso_table[k].A == 238 )      mass_frac = pu_corr.pu239_mass_frac;
          else if( iso_table[k].A == 239 ) mass_frac = pu_corr.pu239_mass_frac;
          else if( iso_table[k].A == 240 ) mass_frac = pu_corr.pu240_mass_frac;
          else if( iso_table[k].A == 241 ) mass_frac = pu_corr.pu241_mass_frac;
          else if( iso_table[k].A == 242 )
          {
            mass_frac = pu_corr.pu242_mass_frac;
            mass_frac_uncert = mass_frac * pu_corr.pu242_uncert;
          }
        }

        json prow;
        prow["nuclide"] = iso_table[k].iso;
        prow["mass_number"] = iso_table[k].A;
        prow["is_by_corr"] = (iso_table[k].A == 242);
        prow["uncorrected_mass_frac"] = (iso_table[k].A == 242) ? -1.0 : uncorr_for[k];
        prow["mass_frac"] = mass_frac;
        prow["mass_frac_uncert"] = mass_frac_uncert;
        pu_data["rows"].push_back( prow );
      }

      pu_data["mass_fractions"] = {
        { "pu238", pu_corr.pu238_mass_frac },
        { "pu239", pu_corr.pu239_mass_frac },
        { "pu240", pu_corr.pu240_mass_frac },
        { "pu241", pu_corr.pu241_mass_frac },
        { "pu242", pu_corr.pu242_mass_frac },
        { "pu242_uncert", pu_corr.pu242_uncert },
        { "is_within_range", pu_corr.is_within_range }
      };

      pu_data["uncorrected_mass_fractions"] = {
        { "pu238", pu_unc.pu238_rel_mass },
        { "pu239", pu_unc.pu239_rel_mass },
        { "pu240", pu_unc.pu240_rel_mass },
        { "pu241", pu_unc.pu241_rel_mass },
        { "pu_other", pu_unc.other_pu_mass },
        { "age",         PhysicalUnits::printToBestTimeUnits(pu_unc.pu_age, 6) },
        { "age_seconds", pu_unc.pu_age }
      };

      // T=0 back-decayed Pu mass fractions, when ages are uniform and there's >1 Pu nuclide.
      if( db && pu_nuclide_ages.size() == 1 && pu_nuclide_names.size() > 1 )
      {
        const SandiaDecay::Nuclide *pu239 = db->nuclide("Pu239");
        const SandiaDecay::Nuclide *pu242 = db->nuclide("Pu242");
        if( pu239 && pu242 && pu_corr.pu239_mass_frac > 0.0 )
        {
          const double pu242_rel_act = pu239_act
                                       * (pu_corr.pu242_mass_frac * pu242->activityPerGram())
                                       / (pu_corr.pu239_mass_frac * pu239->activityPerGram());
          pu_rel_acts.emplace_back( pu242, pu242_rel_act );

          const double back_decay_time = *begin(pu_nuclide_ages);
          try
          {
            vector<tuple<const SandiaDecay::Nuclide *,double,double>> at_t0
                = RelActCalc::back_decay_relative_activities( back_decay_time, pu_rel_acts );
            std::sort( at_t0.begin(), at_t0.end(),
                       []( const tuple<const SandiaDecay::Nuclide *,double,double> &l,
                           const tuple<const SandiaDecay::Nuclide *,double,double> &r ){
                         return get<0>(l)->symbol < get<0>(r)->symbol;
                       } );

            pu_data["back_decayed"] = json::array();
            for( const tuple<const SandiaDecay::Nuclide *,double,double> &nam : at_t0 )
            {
              json brow;
              brow["nuclide"] = get<0>(nam)->symbol;
              brow["mass_frac"] = get<2>(nam);
              pu_data["back_decayed"].push_back( brow );
            }
          }catch( const std::exception &e )
          {
            cerr << "RelActAutoReport: back_decay_relative_activities failed: " << e.what() << endl;
          }
        }
      }

      data["pu_corrections"].push_back( pu_data );
      // Also nest this Pu data inside the matching curve_acts entry for easier templating.
      data["relative_activities"].back()["pu"] = pu_data;
    }
  }

  // ---- Activity / mass ratios (per curve) ----
  data["ratios"] = json::array();
  for( size_t curve_idx = 0; curve_idx < sol.m_rel_activities.size(); ++curve_idx )
  {
    const auto &rel_activities = sol.m_rel_activities[curve_idx];
    json curve_ratios;
    curve_ratios["curve_index"] = static_cast<int64_t>(curve_idx);
    curve_ratios["pairs"] = json::array();

    for( size_t i = 1; i < rel_activities.size(); ++i )
    {
      for( size_t j = 0; j < i; ++j )
      {
        const RelActCalcAuto::NuclideRelAct &nuc_i = rel_activities[i];
        const RelActCalcAuto::NuclideRelAct &nuc_j = rel_activities[j];
        if( RelActCalcAuto::is_null(nuc_i.source) || RelActCalcAuto::is_null(nuc_j.source) )
          continue;

        const SandiaDecay::Nuclide *nuc_i_nuc = RelActCalcAuto::nuclide(nuc_i.source);
        const SandiaDecay::Nuclide *nuc_j_nuc = RelActCalcAuto::nuclide(nuc_j.source);

        json pair_data;
        pair_data["numerator"]   = nuc_i.name();
        pair_data["denominator"] = nuc_j.name();
        const double act_ij = (nuc_j.rel_activity != 0.0) ? (nuc_i.rel_activity / nuc_j.rel_activity) : 0.0;
        const double act_ji = (nuc_i.rel_activity != 0.0) ? (nuc_j.rel_activity / nuc_i.rel_activity) : 0.0;
        pair_data["activity_ratio"] = act_ij;
        pair_data["activity_ratio_inv"] = act_ji;
        {
          char buf[64] = { '\0' };
          snprintf( buf, sizeof(buf), "%.6G", act_ij ); pair_data["activity_ratio_str"] = string(buf);
          snprintf( buf, sizeof(buf), "%.6G", act_ji ); pair_data["activity_ratio_inv_str"] = string(buf);
        }

        if( nuc_i_nuc && nuc_j_nuc )
        {
          const double mass_i = nuc_i.rel_activity / nuc_i_nuc->activityPerGram();
          const double mass_j = nuc_j.rel_activity / nuc_j_nuc->activityPerGram();
          const double mr_ij = (mass_j != 0.0) ? (mass_i / mass_j) : 0.0;
          const double mr_ji = (mass_i != 0.0) ? (mass_j / mass_i) : 0.0;
          pair_data["mass_ratio"] = mr_ij;
          pair_data["mass_ratio_inv"] = mr_ji;
          char buf[64] = { '\0' };
          snprintf( buf, sizeof(buf), "%.6G", mr_ij ); pair_data["mass_ratio_str"] = string(buf);
          snprintf( buf, sizeof(buf), "%.6G", mr_ji ); pair_data["mass_ratio_inv_str"] = string(buf);
        }else
        {
          pair_data["mass_ratio_str"]     = "--";
          pair_data["mass_ratio_inv_str"] = "--";
        }

        // Uncertainty (both directions).
        try
        {
          const double u_ij = sol.activity_ratio_uncertainty( nuc_i.source, curve_idx, nuc_j.source, curve_idx );
          const double u_ji = sol.activity_ratio_uncertainty( nuc_j.source, curve_idx, nuc_i.source, curve_idx );
          pair_data["activity_ratio_uncertainty"]     = u_ij;
          pair_data["activity_ratio_uncertainty_inv"] = u_ji;
          if( act_ij != 0.0 && act_ji != 0.0 )
          {
            char buf[64] = { '\0' };
            snprintf( buf, sizeof(buf), "%.6G%%", 100.0 * u_ij / act_ij );
            pair_data["uncert_pct_str"] = string(buf);
            snprintf( buf, sizeof(buf), "%.6G%%", 100.0 * u_ji / act_ji );
            pair_data["uncert_inv_pct_str"] = string(buf);
          }else
          {
            pair_data["uncert_pct_str"]     = "--";
            pair_data["uncert_inv_pct_str"] = "--";
          }
        }catch( const std::exception & )
        {
          pair_data["uncert_pct_str"]     = "--";
          pair_data["uncert_inv_pct_str"] = "--";
        }

        curve_ratios["pairs"].push_back( pair_data );
      }
    }
    data["ratios"].push_back( curve_ratios );

    // Also nest into the matching curve_acts entry for easier templating.
    if( curve_idx < data["relative_activities"].size() )
      data["relative_activities"][curve_idx]["ratios"] = curve_ratios;
  }

  // ---- Per-peak table ----
  data["peaks"] = json::array();
  if( !sol.m_fit_peaks.empty() )
  {
    vector<pair<PeakDef,size_t>> peak_rel_eff;
    for( size_t r = 0; r < sol.m_fit_peaks_for_each_curve.size(); ++r )
    {
      for( const PeakDef &p : sol.m_fit_peaks_for_each_curve[r] )
        peak_rel_eff.emplace_back( p, r );
    }
    for( const PeakDef &p : sol.m_fit_peaks )
    {
      if( !p.parentNuclide() && !p.xrayElement() && !p.reaction() )
        peak_rel_eff.emplace_back( p, size_t(0) );
    }

    std::sort( peak_rel_eff.begin(), peak_rel_eff.end(),
               []( const pair<PeakDef,size_t> &l, const pair<PeakDef,size_t> &r ){
                 return l.first.mean() < r.first.mean();
               } );

    for( const pair<PeakDef,size_t> &pr : peak_rel_eff )
    {
      const PeakDef &peak = pr.first;
      const size_t rel_eff_index = pr.second;

      json p;
      p["energy"] = peak.mean();
      p["amplitude"] = peak.amplitude();
      p["amplitude_uncertainty"] = peak.amplitudeUncert();
      p["amplitude_uncertainty_percent"] = (peak.amplitude() != 0.0)
        ? 100.0 * peak.amplitudeUncert() / peak.amplitude() : 0.0;
      p["rel_eff_index"] = static_cast<int64_t>(rel_eff_index);

      const SandiaDecay::Nuclide *nuc = peak.parentNuclide();
      const SandiaDecay::Element *el  = peak.xrayElement();
      const ReactionGamma::Reaction *rx = peak.reaction();
      const bool is_floating = (!nuc && !el && !rx);
      p["is_floating"] = is_floating;

      if( nuc )      { p["nuclide"] = nuc->symbol; p["source_type"] = "nuclide"; }
      else if( el )  { p["nuclide"] = el->symbol;  p["source_type"] = "element"; }
      else if( rx )  { p["nuclide"] = rx->name();  p["source_type"] = "reaction"; }
      else           { p["nuclide"] = "";          p["source_type"] = "floating"; }

      // Yield, CPS/yield, fitted rel-eff (only meaningful for nuclide peaks).
      double yield = 0.0;
      double cps_over_yield = 0.0;
      double fit_rel_eff = -1.0, fit_rel_eff_uncert_pct = -1.0;
      if( nuc && !is_floating
          && rel_eff_index < sol.m_rel_activities.size() )
      {
        const RelActCalcAuto::NuclideRelAct *nuc_info = nullptr;
        for( const auto &act : sol.m_rel_activities[rel_eff_index] )
        {
          if( RelActCalcAuto::nuclide(act.source) == nuc )
            nuc_info = &act;
        }
        if( nuc_info )
        {
          for( const pair<double,double> &eb : nuc_info->gamma_energy_br )
          {
            if( fabs(eb.first - peak.gammaParticleEnergy()) < 1.0e-4 )
              yield += eb.second;
          }
          if( yield > 0.0 && live_time > 0.0 )
            cps_over_yield = peak.amplitude() / (yield * live_time);

          try
          {
            const pair<double,double> reu = sol.relative_efficiency_with_uncert( peak.mean(), rel_eff_index );
            fit_rel_eff = reu.first;
            if( reu.first != 0.0 )
              fit_rel_eff_uncert_pct = 100.0 * reu.second / reu.first;
          }catch( const std::exception & )
          {
            try { fit_rel_eff = sol.relative_efficiency( peak.mean(), rel_eff_index ); }
            catch( const std::exception & ) { /* leave as -1 */ }
          }
        }
      }
      p["yield"] = yield;
      p["cps_over_yield"] = cps_over_yield;
      p["fit_rel_eff"] = fit_rel_eff;
      p["fit_rel_eff_uncert_pct"] = fit_rel_eff_uncert_pct;

      if( peak.continuum() )
      {
        p["continuum"] = {
          { "type",         PeakContinuum::offset_type_label_tr( peak.continuum()->type() ) },
          { "lower_energy", peak.continuum()->lowerEnergy() },
          { "upper_energy", peak.continuum()->upperEnergy() }
        };
      }

      data["peaks"].push_back( p );
    }
  }//if( include_peak_details )

  // ---- Energy calibration ----
  bool ene_cal_fit = false;
  for( const bool &fit : sol.m_fit_energy_cal )
    ene_cal_fit = (ene_cal_fit || fit);

  json ene_cal;
  ene_cal["was_fit"] = ene_cal_fit;
  ene_cal["energy_cal_type"] = RelActCalcAuto::to_str( sol.m_options.energy_cal_type );
  ene_cal["adjustments"] = json::array();
  static const char * const ene_cal_label_long[3]  = { "offset adjustment", "gain adjustment", "cubic adjust" };
  static const char * const ene_cal_label_short[3] = { "offset", "gain", "cubic" };
  static const char * const ene_cal_units[3]       = { "keV", "keV/chnl", "keV/chnl2" };
  static const char * const ene_cal_units_long[3]  = { "keV", "keV/channel", "keV/channel\xc2\xb2" };
  for( size_t i = 0; i < sol.m_fit_energy_cal.size() && i < 3; ++i )
  {
    if( !sol.m_fit_energy_cal[i] )
      continue;
    const double phys_value = phys_energy_cal_value( sol, i );
    json adj;
    adj["parameter_index"] = static_cast<int64_t>(i);
    adj["raw_value"] = sol.m_energy_cal_adjustments[i];
    adj["was_fit"] = true;
    adj["physical_value"] = phys_value;
    adj["type"] = ene_cal_label_short[i];
    adj["type_long"] = ene_cal_label_long[i];
    adj["units"] = ene_cal_units[i];
    adj["units_long"] = ene_cal_units_long[i];
    {
      char buf[64] = { '\0' };
      snprintf( buf, sizeof(buf), "%.6G", phys_value );
      adj["physical_value_str"] = string(buf);
    }
    ene_cal["adjustments"].push_back( adj );
  }
  ene_cal["deviation_pair_offsets"] = json::array();
  for( const pair<double,double> &dev : sol.m_deviation_pair_offsets )
  {
    char buf[64] = { '\0' };
    snprintf( buf, sizeof(buf), "{%.1f,%.2f}", dev.first, dev.second );
    json d;
    d["anchor"] = dev.first;
    d["offset"] = dev.second;
    d["display"] = string(buf);
    ene_cal["deviation_pair_offsets"].push_back( d );
  }
  data["energy_calibration"] = ene_cal;

  // ---- Options ----
  data["options"] = {
    { "energy_cal_type",          RelActCalcAuto::to_str( sol.m_options.energy_cal_type ) },
    { "fwhm_form",                RelActCalcAuto::to_str( sol.m_options.fwhm_form ) },
    { "fwhm_estimation_method",   RelActCalcAuto::to_str( sol.m_options.fwhm_estimation_method ) },
    { "skew_type",                static_cast<int>(sol.m_options.skew_type) },
    { "skew_type_str",            PeakDef::to_string(sol.m_options.skew_type) },
    { "additional_br_uncert",     sol.m_options.additional_br_uncert },
    { "spectrum_title",           sol.m_options.spectrum_title }
  };

  // ---- Per-curve options ----
  data["curve_options"] = json::array();
  for( const RelActCalcAuto::RelEffCurveInput &c : sol.m_options.rel_eff_curves )
  {
    json co;
    co["name"] = c.name;
    co["nucs_of_el_same_age"] = c.nucs_of_el_same_age;
    co["rel_eff_eqn_type"] = RelActCalc::to_str( c.rel_eff_eqn_type );
    if( c.rel_eff_eqn_type != RelActCalc::RelEffEqnForm::FramPhysicalModel )
      co["rel_eff_eqn_order"] = static_cast<int64_t>(c.rel_eff_eqn_order);
    data["curve_options"].push_back( co );
  }

  // ---- Foreground / background metadata ----
  // Filename + sample numbers come from `m_options` (callers populate them via
  // `RelActCalcAuto::set_input_spec_info`).  `display_filename` is a pre-formatted
  // "Filename, Samples: {1,3,5}" string -- empty samples means just the filename, per the
  // trivial-omit rule (single-sample file, or fg+bg cover the whole file).
  const auto build_display_filename = []( const string &fname, const std::set<int> &samples ) -> string {
    if( fname.empty() )
      return "";
    if( samples.empty() )
      return fname;
    string s = fname + ", Samples: {";
    bool first = true;
    for( int n : samples ){ s += (first ? "" : ","); s += std::to_string(n); first = false; }
    s += "}";
    return s;
  };

  const auto samples_to_json_array = []( const std::set<int> &samples ) -> nlohmann::json {
    nlohmann::json arr = nlohmann::json::array();
    for( int n : samples )
      arr.push_back( n );
    return arr;
  };

  if( sol.m_foreground )
  {
    data["foreground"] = {
      { "filename",          sol.m_options.foreground_filename },
      { "title",             sol.m_foreground->title() },
      { "sample_numbers",    samples_to_json_array(sol.m_options.foreground_sample_numbers) },
      { "display_filename",  build_display_filename( sol.m_options.foreground_filename,
                                                     sol.m_options.foreground_sample_numbers ) },
      { "live_time",         sol.m_foreground->live_time() },
      { "real_time",         sol.m_foreground->real_time() },
      { "start_time",        SpecUtils::to_iso_string( sol.m_foreground->start_time() ) }
    };
  }
  if( sol.m_background )
  {
    data["background"] = {
      { "filename",          sol.m_options.background_filename },
      { "title",             sol.m_background->title() },
      { "sample_numbers",    samples_to_json_array(sol.m_options.background_sample_numbers) },
      { "display_filename",  build_display_filename( sol.m_options.background_filename,
                                                     sol.m_options.background_sample_numbers ) },
      { "live_time",         sol.m_background->live_time() },
      { "real_time",         sol.m_background->real_time() }
    };
  }

  // ---- Detector metadata ----
  if( sol.m_drf )
  {
    data["detector"] = {
      { "name",              sol.m_drf->name() },
      { "description",       sol.m_drf->description() },
      { "diameter_mm",       sol.m_drf->detectorDiameter() / PhysicalUnits::mm },
      { "has_intrinsic_eff", sol.m_drf->isValid() },
      { "has_fwhm",          sol.m_drf->hasResolutionInfo() }
    };
  }

  // ---- Timing ----
  {
    const int64_t cov_calls = static_cast<int64_t>(sol.m_num_function_eval_total)
                              - static_cast<int64_t>(sol.m_num_function_eval_solution);
    data["timing"] = {
      { "function_eval_solution",   static_cast<int64_t>(sol.m_num_function_eval_solution) },
      { "function_eval_total",      static_cast<int64_t>(sol.m_num_function_eval_total) },
      { "solve_calls",              static_cast<int64_t>(sol.m_num_function_eval_solution) },
      { "cov_calls",                cov_calls },
      { "microseconds_eval",        sol.m_num_microseconds_eval },
      { "microseconds_in_eval",     sol.m_num_microseconds_in_eval },
      { "duration_str",             format_duration( sol.m_num_microseconds_eval ) }
    };
  }

  // ---- Timestamps ----
  string utc_time, local_time;
  if( Wt::WApplication::instance() )
  {
    utc_time   = Wt::WDateTime::currentDateTime().toString("yyyyMMdd hh:mm:ss").toUTF8();
    local_time = Wt::WLocalDateTime::currentDateTime().toString("yyyyMMdd hh:mm:ss").toUTF8();
  }else
  {
    const auto utc_ts = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::now() );
    utc_time = SpecUtils::to_common_string( utc_ts, true );

    std::time_t current_time = std::chrono::system_clock::to_time_t(utc_ts);
    struct tm current_local_time;
#if( defined(WIN32) )
    localtime_s( &current_local_time, &current_time );
#else
    localtime_r( &current_time, &current_local_time );
#endif
    char buf[64] = { '\0' };
    std::strftime( buf, sizeof(buf), "%e-%b-%Y %r", &current_local_time );
    local_time = buf;
    SpecUtils::trim( local_time );
  }
  data["timestamps"] = { { "utc", utc_time }, { "local", local_time } };
  data["compile_timestamp"] = string(__TIMESTAMP__);

  // ---- Asset blobs (D3 + RelEffPlot + SpectrumChartD3 JS/CSS) ----
  // These are large strings (~hundreds of KB each); strip with
  //   jq 'del(.assets, .spectrum_chart.set_js)'
  // when dumping JSON for documentation.
  data["assets"] = json::object();
  for( const auto &kv : load_rel_eff_plot_js_and_css() )
    data["assets"][kv.first] = kv.second;
  for( const auto &kv : load_spectrum_chart_js_and_css() )
    data["assets"][kv.first] = kv.second;

  // ---- Embedded D3 spectrum chart payload ----
  if( sol.m_spectrum )
  {
    try
    {
      // The bundled `set_js` hard-codes this id throughout (resize observer, write_*_for_chart
      // calls, getElementById, etc.).  Templates must place a <div> with this id where they
      // want the chart to render -- the id is exposed as `spectrum_chart.div_id` so they don't
      // have to hard-code the literal string themselves.
      // We append a random suffix so multiple solutions can render side-by-side on the same
      // HTML page (e.g. the multi-file batch summary in §11 of the README) without their
      // spectrum charts colliding.
      string spec_div_id;
      {
        thread_local std::mt19937_64 rng( std::random_device{}() );
        std::uniform_int_distribution<uint64_t> dist;
        std::ostringstream oss;
        oss << "specchart_" << std::hex << std::setw(16) << std::setfill('0') << dist(rng);
        spec_div_id = oss.str();
      }

      stringstream set_js_str;
      D3SpectrumExport::write_js_for_chart( set_js_str, spec_div_id, "", "Energy (keV)", "Counts" );

      // Suffix the observer var so multiple per-call set_js blocks can coexist on the same page
      // without colliding on this top-level `let` (the chart-specific helpers from
      // D3SpectrumExport are already suffixed via spec_div_id).
      set_js_str <<
      "  let spec_observer_" << spec_div_id << " = new ResizeObserver(entries => {\n"
      "    for (let entry of entries) {\n"
      "      if (entry.target && (entry.target.id === \"" << spec_div_id << "\")) {\n"
      "        spec_chart_" << spec_div_id << ".handleResize(false);\n"
      "      }\n"
      "    }\n"
      "  });\n"
      "  spec_observer_" << spec_div_id << ".observe( document.getElementById(\"" << spec_div_id << "\") );\n";

      D3SpectrumExport::D3SpectrumChartOptions chart_options;
      chart_options.m_useLogYAxis = true;
      chart_options.m_legendEnabled = false;
      chart_options.m_compactXAxis = true;
      chart_options.m_allowDragRoiExtent = false;
      D3SpectrumExport::write_set_options_for_chart( set_js_str, spec_div_id, chart_options );
      set_js_str << "  spec_chart_" << spec_div_id << ".setShowLegend(false);\n";

      // Initial display energy: ROI bounds + 10% padding, clamped to spectrum range.
      float lower_energy = sol.m_spectrum->gamma_energy_max();
      float upper_energy = sol.m_spectrum->gamma_energy_min();
      for( const RelActCalcAuto::RoiRange &rr : sol.m_options.rois )
      {
        lower_energy = std::min( lower_energy, static_cast<float>(rr.lower_energy) );
        upper_energy = std::max( upper_energy, static_cast<float>(rr.upper_energy) );
      }
      const float range = upper_energy - lower_energy;
      lower_energy -= 0.1f * range;
      upper_energy += 0.1f * range;
      lower_energy = std::max( sol.m_spectrum->gamma_energy_min(), lower_energy );
      upper_energy = std::min( sol.m_spectrum->gamma_energy_max(), upper_energy );
      if( upper_energy <= lower_energy )
      {
        lower_energy = sol.m_spectrum->gamma_energy_min();
        upper_energy = sol.m_spectrum->gamma_energy_max();
      }
      set_js_str << "  spec_chart_" << spec_div_id
                 << ".setXAxisRange(" << lower_energy << ", " << upper_energy << ", false);\n";

      D3SpectrumExport::D3SpectrumOptions spec_options;
      spec_options.spectrum_type = SpecUtils::SpectrumType::Foreground;
      vector<shared_ptr<const PeakDef>> peaks;
      for( const PeakDef &p : sol.m_fit_peaks )
        peaks.push_back( make_shared<PeakDef>(p) );
      spec_options.peaks_json = PeakDef::peak_json( peaks, sol.m_spectrum, Wt::WColor(0,51,255), 255 );

      const SpecUtils::Measurement * const meas_ptr = sol.m_spectrum.get();
      D3SpectrumExport::write_and_set_data_for_chart( set_js_str, spec_div_id,
          { std::make_pair(meas_ptr, spec_options) } );

      data["spectrum_chart"] = {
        { "have_spectrum", true },
        { "div_id",        spec_div_id },
        { "set_js",        set_js_str.str() },
        { "lower_energy",  lower_energy },
        { "upper_energy",  upper_energy }
      };
    }catch( const std::exception &e )
    {
      cerr << "RelActAutoReport: spectrum_chart payload failed: " << e.what() << endl;
      data["spectrum_chart"] = { { "have_spectrum", false } };
    }
  }else
  {
    data["spectrum_chart"] = { { "have_spectrum", false } };
  }

  // ---- ROIs ----
  data["rois"] = json::array();
  for( const RelActCalcAuto::RoiRange &r : sol.m_options.rois )
  {
    data["rois"].push_back( {
      { "lower_energy",     r.lower_energy },
      { "upper_energy",     r.upper_energy },
      { "continuum_type",   PeakContinuum::offset_type_label_tr( r.continuum_type ) },
      { "range_limits_type", RelActCalcAuto::RoiRange::to_str( r.range_limits_type ) }
    } );
  }

  return data;
}//solution_to_json


// ---- render_template ------------------------------------------------------

namespace
{
  /** Turn the (possibly user-supplied) custom_template into a "html"/"txt"/"json" key,
   or - if it is a file - return the resolved path as the second element of the pair.
   First element is the format key used to pick the named default template.
   An empty `tmplt` defaults to the HTML template. */
  pair<string,string> resolve_template( const string &tmplt_in, const string &include_dir )
  {
    string t = tmplt_in;
    SpecUtils::trim( t );
    if( t.empty() ) return { "html", "" };

    if( SpecUtils::iequals_ascii(t, "html") ) return { "html", "" };
    if( SpecUtils::iequals_ascii(t, "txt") || SpecUtils::iequals_ascii(t, "text") ) return { "txt", "" };
    if( SpecUtils::iequals_ascii(t, "json") ) return { "json", "" };
    if( SpecUtils::iequals_ascii(t, "html-summary") || SpecUtils::iequals_ascii(t, "summary") ) return { "html-summary", "" };

    // Treat as a filename - try, in order: user-supplied include dir, the
    // user's writable IsotopicsByNuclidesReportTmplts/ dir, an absolute path,
    // then the bundled default dir.
    const string tmplt_dir = template_include_dir( include_dir );
    if( !tmplt_dir.empty() )
    {
      const string full = SpecUtils::append_path( tmplt_dir, t );
      if( SpecUtils::is_file(full) )
        return { "", full };
    }
    const string user_dir = writable_template_dir();
    if( !user_dir.empty() )
    {
      const string full = SpecUtils::append_path( user_dir, t );
      if( SpecUtils::is_file(full) )
        return { "", full };
    }
    if( SpecUtils::is_file(t) )
      return { "", t };

    const string default_dir = default_template_dir();
    const string full = SpecUtils::append_path( default_dir, t );
    if( SpecUtils::is_file(full) )
      return { "", full };

    throw runtime_error( "RelActAutoReport: could not find template '" + t + "'" );
  }
}//namespace


string render_template( inja::Environment &env,
                        const nlohmann::json &data,
                        const string &tmplt_in,
                        const string &include_dir )
{
  const pair<string,string> tmplt = resolve_template( tmplt_in, include_dir );

  // JSON dump: skip Inja entirely.
  if( tmplt.first == "json" )
    return data.dump( 2 );

  if( tmplt.first == "html" )
    return env.render( "{% include \"default-rel-act-auto-html-results\" %}", data );
  if( tmplt.first == "txt" )
    return env.render( "{% include \"default-rel-act-auto-txt-results\" %}", data );
  if( tmplt.first == "html-summary" )
    return env.render( "{% include \"default-rel-act-auto-html-multi-file-summary\" %}", data );

  // Filename path.  We can't just call `env.parse_template(full_path)`: Inja's
  // parse_template does `env.input_path + filename` unconditionally, so passing
  // a full path while the env was built with a non-empty input_path produces
  // garbage like "/include_dir//tmp/foo.tmplt.html".
  //
  // Two cases:
  //   1) The user explicitly named a `--report-template-include-dir`. We
  //      require the resolved template path to live under it (otherwise their
  //      `{% include %}` directives will resolve against the wrong directory),
  //      and pass the relative-to-include-dir name to `parse_template`.
  //   2) Otherwise, we read the file content ourselves and feed it to
  //      `env.parse(string)`.  This sidesteps Inja's path-prepending entirely;
  //      `{% include %}` directives still resolve against the env's input_path
  //      (which is whatever the caller built the env with -- typically the
  //      bundled default dir).
  // Canonicalize so the prefix-check below works on symlinked paths
  // (notably macOS, where `/tmp` -> `/private/tmp`).
  string full_path = tmplt.second;
  {
    string canon = full_path;
    if( SpecUtils::make_canonical_path(canon) )
      full_path = canon;
  }

  
  // Read the template file directly and feed to env.parse().  We avoid
  // env.parse_template(name) because Inja resolves it against env.input_path,
  // which the caller's env may not have set up to match the template's directory.
  std::vector<char> content;
  try
  {
    SpecUtils::load_file_data( full_path.c_str(), content );
  }catch( std::exception &e )
  {
    // Okay - we'll try using `env` to read the file - we're desperate here
    const bool explicit_include_dir
      = !include_dir.empty()
      && !SpecUtils::iequals_ascii(include_dir, "none")
      && !SpecUtils::iequals_ascii(include_dir, "default");

    if( explicit_include_dir )
    {
      string include_canon = include_dir;
      if( !SpecUtils::make_canonical_path(include_canon) )
        throw runtime_error( "RelActAutoReport: include directory '" + include_dir
                           + "' does not exist or is not accessible." );
      if( !include_canon.empty() && include_canon.back() != ns_path_sep )
        include_canon += ns_path_sep;

      if( !SpecUtils::starts_with(full_path, include_canon.c_str()) )
        throw runtime_error( "RelActAutoReport: template '" + full_path
                            + "' is not under the specified include directory '"
                            + include_canon + "'.  Either pass a template path under that"
                            " directory, or omit --report-template-include-dir to let it be"
                            " inferred from the template's parent directory." );
      const string relative_name = full_path.substr( include_canon.size() );
      inja::Template parsed = env.parse_template( relative_name );
      return env.render( parsed, data );
    }//if( explicit_include_dir )

    throw runtime_error( "RelActAutoReport: could not read template '" + full_path
                        + "': " + e.what() );
  }//try / catch 

  if( !content.empty() && content.back() == '\0' )
    content.pop_back();
  const std::string content_str( content.begin(), content.end() );

  inja::Template parsed = env.parse( content_str );
  return env.render( parsed, data );
}


void render_template( ostream &out,
                      inja::Environment &env,
                      const nlohmann::json &data,
                      const string &tmplt,
                      const string &include_dir )
{
  out << render_template( env, data, tmplt, include_dir );
}

}//namespace RelActAutoReport
