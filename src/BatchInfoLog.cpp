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
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <exception>


#include "external_libs/SpecUtils/3rdparty/inja/inja.hpp"

#include "SpecUtils/DateTime.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/D3SpectrumExport.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/BatchPeak.h"
#include "InterSpec/EnergyCal.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/BatchInfoLog.h"
#include "InterSpec/BatchSampleSelect.h"
#include "InterSpec/BatchActivity.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/ShowRiidInstrumentsAna.h"
#include "InterSpec/ShieldSourcePullTrend.h"
#include "InterSpec/ShieldingSourceFitPlot.h"

using namespace std;

namespace
{
#if( defined(_WIN32) )
  const char ns_path_sep = '\\';
#else
  const char ns_path_sep = '/';
#endif
}//namespace


namespace BatchInfoLog
{
  std::string default_template_dir()
  {
    const string static_data_dir = InterSpec::staticDataDirectory().empty() ? string("./data")
                                                                  : InterSpec::staticDataDirectory();
    //Also see: `WServer::instance()->appRoot()`, which isnt valid right now
    const string app_root = SpecUtils::append_path( static_data_dir, ".." );
    const string docroot  = SpecUtils::append_path( app_root, "InterSpec_resources" );
    const string static_txt = SpecUtils::append_path( docroot, "static_text" );
    
    // inja assumes trailing path separator, for its template path
    
    return SpecUtils::append_path( static_txt, "ShieldSourceFitLog" ) + ns_path_sep;
  }//std::string default_template_dir()
  
  
  std::string template_include_dir( const BatchPeak::BatchPeakFitOptions &options )
  {
    if( SpecUtils::iequals_ascii( options.template_include_dir, "none" ) )
      return "";
    
    if( SpecUtils::iequals_ascii( options.template_include_dir, "default" ) )
      return default_template_dir();
    
    string tmplt_dir = options.template_include_dir;
    if( !tmplt_dir.empty() && (tmplt_dir.back() != ns_path_sep) )
      tmplt_dir += ns_path_sep;
    
    return tmplt_dir;
  }//std::string template_include_dir( const BatchPeak::BatchPeakFitOptions &options )
  
  
  inja::Environment get_default_inja_env( const string &tmplt_dir )
  {
    if( !tmplt_dir.empty() && !SpecUtils::is_directory(tmplt_dir) )
      throw runtime_error( string("Template include directory, '") + tmplt_dir
                          + "', doesnt look to be a valid directory - not performing analysis." );
    
    
    inja::Environment env{ tmplt_dir };
    env.set_trim_blocks( true ); // remove the first newline after a block
    //env.set_lstrip_blocks( true ); //strip the spaces and tabs from the start of a line to a block
    
    
  #if( BUILD_FOR_WEB_DEPLOYMENT )
      env.set_search_included_templates_in_files( false );
  #else
    //To think about more: is there any security issues with allowing.
    //  It doesnt *look* like inja prevents using things like "../../../SomeSensitveFile.txt" as templates.
    env.set_search_included_templates_in_files( !tmplt_dir.empty() );
  #endif
    
    // Add some callbacks incase people want more control over the precision of their printouts
    env.add_callback( "printFixed", 2, &BatchInfoLog::printFixed );
    env.add_callback( "printCompact", 2, &BatchInfoLog::printCompact );
    
    try
    {
      // If we're using a custom include path, opening templates from the default template location
      //  is problematic, so we'll just create a new `inja::Environment` to open the default
      //  templates, and then add them to `env` - not really tested yet.
      inja::Environment sub_env;
      sub_env.add_callback( "printFixed", 2, &BatchInfoLog::printFixed );
      sub_env.add_callback( "printCompact", 2, &BatchInfoLog::printCompact );
      
      const string default_tmplt_dir = BatchInfoLog::default_template_dir();
      
      {
        const string def_txt_tmplt = SpecUtils::append_path( default_tmplt_dir, "std_fit_log.tmplt.txt" );
        inja::Template txt_tmplt = sub_env.parse_template( def_txt_tmplt );
        env.include_template( "default-act-fit-txt-results", txt_tmplt );
      }
      
      {
        const string def_html_tmplt = SpecUtils::append_path( default_tmplt_dir, "act_fit.tmplt.html" );
        inja::Template html_tmplt = sub_env.parse_template( def_html_tmplt );
        env.include_template( "default-act-fit-html-results", html_tmplt );
      }
      
      {
        const string def_csv_summary_tmplt = SpecUtils::append_path( default_tmplt_dir, "std_summary.tmplt.csv" );
        inja::Template csv_sum_tmplt = sub_env.parse_template( def_csv_summary_tmplt );
        env.include_template( "default-act-fit-csv-summary", csv_sum_tmplt );
      }
      
      {
        const string def_html_summary_tmplt = SpecUtils::append_path( default_tmplt_dir, "std_summary.tmplt.html" );
        inja::Template html_sum_tmplt = sub_env.parse_template( def_html_summary_tmplt );
        env.include_template( "default-act-fit-html-summary", html_sum_tmplt );
      }
      
      {
        const string tmplt_path = SpecUtils::append_path( default_tmplt_dir, "std_peak_fit_log.tmplt.txt" );
        inja::Template tmplt = sub_env.parse_template( tmplt_path );
        env.include_template( "default-peak-fit-txt-results", tmplt );
      }
      
      {
        const string tmplt_path = SpecUtils::append_path( default_tmplt_dir, "peak_fit.tmplt.html" );
        inja::Template tmplt = sub_env.parse_template( tmplt_path );
        env.include_template( "default-peak-fit-html-results", tmplt );
      }
      
      {
        const string tmplt_path = SpecUtils::append_path( default_tmplt_dir, "std_peak_fit_summary.tmplt.csv" );
        inja::Template tmplt = sub_env.parse_template( tmplt_path );
        env.include_template( "default-peak-fit-csv-summary", tmplt );
      }
      
      {
        const string tmplt_path = SpecUtils::append_path( default_tmplt_dir, "std_peak_fit_summary.tmplt.html" );
        inja::Template tmplt = sub_env.parse_template( tmplt_path );
        env.include_template( "default-peak-fit-html-summary", tmplt );
      }
    }catch( std::exception &e )
    {
      throw runtime_error( "Error loading default analysis report template: " + string(e.what()) );
    }
    
    return env;
  }//inja::Environment get_default_inja_env( const string &tmplt_dir )
  
  
  vector<pair<string,string>> load_spectrum_chart_js_and_css()
  {
    vector<pair<string,string>> answer;
    
#if( SpecUtils_ENABLE_D3_CHART )
    //InterSpec::setStaticDataDirectory( SpecUtils::append_path(datadir,"data") );
    assert( !InterSpec::staticDataDirectory().empty() );
    const string static_data_dir = InterSpec::staticDataDirectory().empty() ? string("./data")
    : InterSpec::staticDataDirectory();
    //Also see: `WServer::instance()->appRoot()`, which isnt valid right now
    const string app_root = SpecUtils::append_path( static_data_dir, ".." );
    const string docroot  = SpecUtils::append_path( app_root, "InterSpec_resources" );

    // d3.v3.min.js and SpectrumChartD3.{js,css} are deployed from
    // external_libs/SpecUtils/d3_resources/ into the build's InterSpec_resources/ at build
    // time (not committed to source).  In the unit-test build the staticDataDirectory points
    // at the source tree, so the deployed files may not be there yet — fall back to the
    // SpecUtils source location so tests can still embed the assets.
    const auto pick_path = [&docroot,&app_root]( const string &filename ) -> string {
      const string deployed = SpecUtils::append_path( docroot, filename );
#if( BUILD_AS_UNIT_TEST_SUITE )
      if( !SpecUtils::is_file(deployed) )
      {
        const string src = SpecUtils::append_path(
          SpecUtils::append_path( app_root, "external_libs/SpecUtils/d3_resources" ),
          filename );
        if( SpecUtils::is_file(src) )
          return src;
      }
#endif
      return deployed;
    };

    const string sc_js_fn = SpecUtils::is_file( SpecUtils::append_path( docroot, "SpectrumChartD3.min.js") )
                            ? "SpectrumChartD3.min.js" : "SpectrumChartD3.js";
    const string sc_css_fn = SpecUtils::is_file( SpecUtils::append_path( docroot, "SpectrumChartD3.min.css") )
                            ? "SpectrumChartD3.min.css" : "SpectrumChartD3.css";

    string d3_js  = AppUtils::file_contents( pick_path( "d3.v3.min.js" ) );

    string sc_js  = AppUtils::file_contents( pick_path( sc_js_fn ) );
    string sc_css = AppUtils::file_contents( pick_path( sc_css_fn ) );

    answer.emplace_back( "D3_JS", std::move(d3_js) );
    answer.emplace_back( "SpectrumChart_JS", std::move(sc_js) );
    answer.emplace_back( "SpectrumChart_CSS", std::move(sc_css) );
#endif // SpecUtils_ENABLE_D3_CHART
    
    return answer;
  }//load_spectrum_chart_js_and_css()


  vector<pair<string,string>> load_shielding_fit_plot_js_and_css()
  {
    vector<pair<string,string>> answer;

#if( SpecUtils_ENABLE_D3_CHART )
    //InterSpec::setStaticDataDirectory( SpecUtils::append_path(datadir,"data") );
    assert( !InterSpec::staticDataDirectory().empty() );
    const string static_data_dir = InterSpec::staticDataDirectory().empty() ? string("./data")
    : InterSpec::staticDataDirectory();
    //Also see: `WServer::instance()->appRoot()`, which isnt valid right now
    const string app_root = SpecUtils::append_path( static_data_dir, ".." );
    const string docroot  = SpecUtils::append_path( app_root, "InterSpec_resources" );

    // d3.v3.min.js is deployed from external_libs/SpecUtils/d3_resources/ at build time;
    // see load_spectrum_chart_js_and_css() for the BUILD_AS_UNIT_TEST_SUITE rationale.
    string d3_js_path = SpecUtils::append_path( docroot, "d3.v3.min.js" );
#if( BUILD_AS_UNIT_TEST_SUITE )
    if( !SpecUtils::is_file(d3_js_path) )
    {
      const string src = SpecUtils::append_path(
        SpecUtils::append_path( app_root, "external_libs/SpecUtils/d3_resources" ),
        "d3.v3.min.js" );
      if( SpecUtils::is_file(src) )
        d3_js_path = src;
    }
#endif
    string d3_js = AppUtils::file_contents( d3_js_path );

    string plot_js  = AppUtils::file_contents( SpecUtils::append_path( docroot, "ShieldingSourceFitPlot.js" ) );
    string plot_css = AppUtils::file_contents( SpecUtils::append_path( docroot, "ShieldingSourceFitPlot.css" ) );

    answer.emplace_back( "D3_JS", std::move(d3_js) );
    answer.emplace_back( "ShieldingSourceFitPlot_JS", std::move(plot_js) );
    answer.emplace_back( "ShieldingSourceFitPlot_CSS", std::move(plot_css) );
#endif // SpecUtils_ENABLE_D3_CHART

    return answer;
  }//load_shielding_fit_plot_js_and_css()

  string render_template( string tmplt, 
                         inja::Environment &env,
                         const TemplateRenderType type,
                         const BatchPeak::BatchPeakFitOptions &options,
                         const nlohmann::json &data )
  {
    string rpt;
    
    switch( type )
    {
      case TemplateRenderType::ActShieldIndividual:
      {
        if( SpecUtils::iequals_ascii(tmplt, "txt" ) )
          rpt = env.render("{% include \"default-act-fit-txt-results\" %}", data);
        else if( SpecUtils::iequals_ascii(tmplt, "html" ) )
          rpt = env.render("{% include \"default-act-fit-html-results\" %}", data);
        break;
      }//case TemplateRenderType::ActShieldIndividual:
        
      case TemplateRenderType::ActShieldSummary:
      {
        if( SpecUtils::iequals_ascii(tmplt, "csv" ) )
          rpt = env.render( "{% include \"default-act-fit-csv-summary\" %}", data );
        else if( SpecUtils::iequals_ascii(tmplt, "html" ) )
          rpt = env.render( "{% include \"default-act-fit-html-summary\" %}", data );
        break;
      }//case TemplateRenderType::ActShieldSummary:
        
      case TemplateRenderType::PeakFitIndividual:
      {
        if( SpecUtils::iequals_ascii(tmplt, "txt" ) )
          rpt = env.render("{% include \"default-peak-fit-txt-results\" %}", data);
        else if( SpecUtils::iequals_ascii(tmplt, "html" ) )
          rpt = env.render("{% include \"default-peak-fit-html-results\" %}", data);
        break;
      }//case TemplateRenderType::PeakFitIndividual:
      
      case TemplateRenderType::PeakFitSummary:
      {
        if( SpecUtils::iequals_ascii(tmplt, "csv" ) )
          rpt = env.render( "{% include \"default-peak-fit-csv-summary\" %}", data );
        else if( SpecUtils::iequals_ascii(tmplt, "html" ) )
          rpt = env.render( "{% include \"default-peak-fit-html-summary\" %}", data );
        break;
      }//case TemplateRenderType::PeakFitSummary:
    }//switch( type )
    
    
    if( rpt.empty() )
    {
      const string tmplt_dir = BatchInfoLog::template_include_dir( options );
      const bool is_in_inc = SpecUtils::is_file( SpecUtils::append_path(tmplt_dir, tmplt) );
      
      inja::Template injatmplt;
      if( is_in_inc )
      {
        injatmplt = env.parse_template( tmplt );
      }else
      {
        bool is_file = SpecUtils::is_file( tmplt );
        if( !is_file )
        {
          const string default_tmplt_dir = BatchInfoLog::default_template_dir();
          const string tmplt_in_def_path = SpecUtils::append_path(default_tmplt_dir, tmplt);
          
          is_file = SpecUtils::is_file( tmplt_in_def_path );
          
#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
          // TODO: consider using `AppUtils::locate_file(...)` to find the file
#endif
          
          if( !is_file )
            throw runtime_error( "Could not find template file '" + tmplt + "'."
                                " Please specify full path to file, or use the"
                                " 'report-template-include-dir' option to specify"
                                " directory where reports are located." );
          tmplt = tmplt_in_def_path;
        }
        
        inja::Environment sub_env;
        sub_env.add_callback( "printFixed", 2, &BatchInfoLog::printFixed );
        sub_env.add_callback( "printCompact", 2, &BatchInfoLog::printCompact );
        injatmplt = sub_env.parse_template( tmplt );
      }//
      
      rpt = env.render( injatmplt, data );
    }//if( default report format ) / else
    
    return rpt;
  };//render_template(...)
  
  
  std::string suggested_output_report_filename( const std::string &filename, 
                                               const std::string tmplt,
                                               const TemplateRenderType type,
                                               const BatchPeak::BatchPeakFitOptions &options )
  {
    string outname = SpecUtils::filename( filename );
    const string file_ext = SpecUtils::file_extension(outname);
    if( !file_ext.empty() )
      outname = outname.substr(0, outname.size() - file_ext.size());
    
    string tmplt_name = SpecUtils::filename( tmplt );
    string tmplt_ext = SpecUtils::file_extension(tmplt_name);
    
    if( SpecUtils::iequals_ascii(tmplt, "txt")
       || SpecUtils::iequals_ascii(tmplt, "text")
       || SpecUtils::iequals_ascii(tmplt, "csv")
       || SpecUtils::iequals_ascii(tmplt, "html") )
    {
      switch( type )
      {
        case TemplateRenderType::ActShieldIndividual:
          tmplt_name = "act_fit";
          break;
          
        case TemplateRenderType::ActShieldSummary:
          tmplt_name = "summary";
          break;
          
        case TemplateRenderType::PeakFitIndividual:
          tmplt_name = "peak_fit";
          break;
          
        case TemplateRenderType::PeakFitSummary:
          tmplt_name = "summary";
          break;
      }//switch( type )
      
      tmplt_ext = "." + tmplt;
      SpecUtils::to_lower_ascii( tmplt_ext );
    }//if( template is "txt", "text", "csv", or "html" )
    
    
    size_t pos = SpecUtils::ifind_substr_ascii(tmplt_name, "tmplt");
    if( pos == string::npos )
      pos = SpecUtils::ifind_substr_ascii(tmplt_name, "template");
    if( pos != string::npos )
      tmplt_name = tmplt_name.substr(0, pos);
    if( SpecUtils::iends_with(tmplt_name, "_") 
       || SpecUtils::iends_with(tmplt_name, ".")
       || SpecUtils::iends_with(tmplt_name, "-") )
    {
      tmplt_name = tmplt_name.substr(0, tmplt_name.size() - 1);
    }
    
    if( tmplt_ext.empty()
       || SpecUtils::iequals_ascii(tmplt_ext, "tmplt" )
       || SpecUtils::iequals_ascii(tmplt_ext, "template" ) )
      tmplt_ext = SpecUtils::file_extension(tmplt_name);
    
    if( tmplt_ext.empty() )
      tmplt_ext = ".txt";
    
    outname += (outname.empty() ? "" : "_") + tmplt_name + tmplt_ext;
    return SpecUtils::append_path(options.output_dir, outname );
  }//std::string suggested_output_report_filename(...)
  
  
  std::string printFixed( std::vector<const nlohmann::json *> &args )
  {
    try
    {
      if( args.empty() )
        return "";
      if( args[0]->is_null() )
        return "--";  //This happens if you try to put a inf or NaN as the value - the JSON will just have a null object, since JSON doesnt support these values
      if( args[0]->is_string() )
        return args[0]->get<string>();
      if( !args[0]->is_number() )
        throw runtime_error( "not a number, like expected." );

      const double val = args.at(0)->get<double>();
      const int numDecimal = std::max( 0, args.at(1)->get<int>() );
      
      char buffer[64] = { '\0' };
      snprintf( buffer, sizeof(buffer), "%.*f", numDecimal, val );
      
      return std::string(buffer);
    }catch( inja::InjaError &e )
    {
      const string msg = "Error converting 'printFixed' argument to number.\n"
      "line " + std::to_string(e.location.line) + ", column " + std::to_string(e.location.column)
      + "): " + e.message + ".";
      
      cerr << msg << endl;
      throw;
    }catch( std::exception &e )
    {
      cerr << "Error in 'printFixed': " << e.what() << endl;
      return "ErrorPrintingValue{" + string(e.what()) + "}";
    }
  };
  
  std::string printCompact( std::vector<const nlohmann::json *> &args )
  {
    try
    {
      if( args.empty() )
        return "";
      if( args[0]->is_null() )
        return "--";  //This happens if you try to put a inf or NaN as the value - the JSON will just have a null object, since JSON doesnt support these values
      if( args[0]->is_string() )
        return args[0]->get<string>();
      if( !args[0]->is_number() )
        throw runtime_error( "not a number, like expected." );
      
      const double val = args.at(0)->get<double>();
      const int numSigFig = args.at(1)->get<int>();
      if( numSigFig <= 1 )
        throw runtime_error( "printCompact: you must print at least one significant figures" );
      
      return SpecUtils::printCompact( val, static_cast<size_t>(numSigFig) );
    }catch( inja::InjaError &e )
    {
      const string msg = "Error converting 'printCompact' argument to number.\n"
      "line " + std::to_string(e.location.line) + ", column " + std::to_string(e.location.column)
      + "): " + e.message + ".";
      
      cerr << msg << endl;
      throw;
    }catch( std::exception &e )
    {
      cerr << "Error in 'printCompact': " << e.what() << endl;
      return "ErrorPrintingValue{" + string(e.what()) + "}";
    }
    return "";
  };
  
  
  // Adds the basic direct info on a source (nuclide name, activity, age, etc), but does not
  //  Add which peaks it contributes to, or any information on gammas
void add_basic_src_details( const GammaInteractionCalc::SourceDetails &src,
                          const std::shared_ptr<const DetectorPeakResponse> &drf,
                          const bool useBq,
                          const std::vector<GammaInteractionCalc::ShieldingDetails> *shield_details,
                          nlohmann::basic_json<> &src_json )
{
  const DetectorPeakResponse::EffGeometryType eff_type = drf ? drf->geometryType()
                                              : DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic;
  
  const string act_postfix = DetectorPeakResponse::det_eff_geom_type_postfix(eff_type);
  
  
  src_json["Nuclide"] = src.nuclide->symbol;
  src_json["Activity"] = PhysicalUnits::printToBestActivityUnits(src.activity,4,!useBq) + act_postfix;
  src_json["Activity_bq"] = src.activity / PhysicalUnits::bq;
  src_json["Activity_kBq"] = src.activity / PhysicalUnits::kBq;
  src_json["Activity_MBq"] = src.activity / PhysicalUnits::MBq;
  src_json["Activity_GBq"] = src.activity / PhysicalUnits::GBq;
  src_json["Activity_ci"] = src.activity / PhysicalUnits::ci;
  src_json["Activity_mCi"] = src.activity / PhysicalUnits::mCi;
  src_json["Activity_uCi"] = src.activity / PhysicalUnits::microCi;
  src_json["Activity_pCi"] = src.activity / PhysicalUnits::pCi;
  src_json["ActivityPostFix"] = act_postfix;
  
  bool activityIsFit = src.activityIsFit;
  if( src.isSelfAttenSource )
  {
    activityIsFit |= src.isSelfAttenVariableMassFrac;
   
    assert( shield_details );
    if( shield_details )
    {
      // TODO: check if mass fraction is being fit, or if a shielding dimension is being fit
      assert( src.selfAttenShieldIndex < shield_details->size() );
      const GammaInteractionCalc::ShieldingDetails &shield = shield_details->at(src.selfAttenShieldIndex);
      assert( shield.m_name == src.selfAttenShieldName );
      for( unsigned int i = 0; i < shield.m_num_dimensions; ++i )
        activityIsFit |= shield.m_fit_dimension[i];
    }//if( shield_details )
  }//if( src.isSelfAttenSource )
  
  src_json["ActivityIsFit"] = activityIsFit;
  
  
  if( activityIsFit )
  {
    src_json["ActivityUncert"] = PhysicalUnits::printToBestActivityUnits(src.activityUncertainty,4,!useBq) + act_postfix;
    
    src_json["ActivityUncert_bq"]  = src.activityUncertainty / PhysicalUnits::bq;
    src_json["ActivityUncert_kBq"] = src.activityUncertainty / PhysicalUnits::kBq;
    src_json["ActivityUncert_MBq"] = src.activityUncertainty / PhysicalUnits::MBq;
    src_json["ActivityUncert_GBq"] = src.activityUncertainty / PhysicalUnits::GBq;
    src_json["ActivityUncert_ci"]  = src.activityUncertainty / PhysicalUnits::ci;
    src_json["ActivityUncert_mCi"] = src.activityUncertainty / PhysicalUnits::mCi;
    src_json["ActivityUncert_uCi"] = src.activityUncertainty / PhysicalUnits::microCi;
    src_json["ActivityUncert_pCi"] = src.activityUncertainty / PhysicalUnits::pCi;
    
    const double act_uncert_percent = 100.0 * src.activityUncertainty / src.activity;
    src_json["ActivityUncertPercent"] = SpecUtils::printCompact( act_uncert_percent, 4 );
  }else
  {
    assert( src.activityUncertainty <= 0.0 );
  }
  
  src_json["NuclideMass"] = PhysicalUnits::printToBestMassUnits(src.nuclideMass, 4);
  src_json["Age"] = PhysicalUnits::printToBestTimeUnits(src.age, 4);
  src_json["AgeSeconds"] = src.age / PhysicalUnits::second;
  src_json["AgeDays"] = src.age / PhysicalUnits::day;
  src_json["AgeYears"] = src.age / PhysicalUnits::year;
  src_json["AgeIsFittable"] = src.ageIsFittable;
  src_json["AgeIsFit"] = src.ageIsFit;
  
  if( src.ageIsFit )
  {
    src_json["AgeUncert"] = PhysicalUnits::printToBestTimeUnits(src.ageUncertainty, 4);
    src_json["AgeUncertSeconds"] = src.ageUncertainty / PhysicalUnits::second;
    src_json["AgeUncertDays"] = src.ageUncertainty / PhysicalUnits::day;
    src_json["AgeUncertYears"] = src.ageUncertainty / PhysicalUnits::year;
  }else
  {
    assert( src.ageUncertainty <= 0.0 );
  }
  
  if( src.ageDefiningNuc )
    src_json["AgeDefiningNuclide"] = src.ageDefiningNuc->symbol;
  
  src_json["IsTraceSource"] = src.isTraceSource;
  
  if( src.isTraceSource )
  {
    string trace_src_postfix = "";
    
    switch( src.traceActivityType )
    {
      case GammaInteractionCalc::TraceActivityType::TotalActivity:
        trace_src_postfix = "";
        break;
        
      case GammaInteractionCalc::TraceActivityType::ActivityPerCm3:
        trace_src_postfix = "/cm^3";
        break;
      
      case GammaInteractionCalc::TraceActivityType::ExponentialDistribution:
        trace_src_postfix = "/m2 exp";
        break;
        
      case GammaInteractionCalc::TraceActivityType::ActivityPerGram:
        trace_src_postfix = "/g";
        break;
        
      case GammaInteractionCalc::TraceActivityType::NumTraceActivityType:
        break;
    }//switch( type )
    
    src_json["TraceActivityType"] = GammaInteractionCalc::to_str(src.traceActivityType);
    src_json["TraceDisplayActivity"] = PhysicalUnits::printToBestActivityUnits(src.traceSrcDisplayAct,4,!useBq) + trace_src_postfix;
    src_json["TraceDisplayActivity_bq"] = src.traceSrcDisplayAct / PhysicalUnits::bq;
    src_json["TraceDisplayActivity_kBq"] = src.traceSrcDisplayAct / PhysicalUnits::kBq;
    src_json["TraceDisplayActivity_MBq"] = src.traceSrcDisplayAct / PhysicalUnits::MBq;
    src_json["TraceDisplayActivity_GBq"] = src.traceSrcDisplayAct / PhysicalUnits::GBq;
    src_json["TraceDisplayActivity_ci"]  = src.traceSrcDisplayAct / PhysicalUnits::ci;
    src_json["TraceDisplayActivity_mCi"] = src.traceSrcDisplayAct / PhysicalUnits::mCi;
    src_json["TraceDisplayActivity_uCi"] = src.traceSrcDisplayAct / PhysicalUnits::microCi;
    src_json["TraceDisplayActivity_pCi"] = src.traceSrcDisplayAct / PhysicalUnits::pCi;
    
    src_json["TraceActivityPostFix"] = trace_src_postfix;
    
    if( src.activityIsFit )
    {
      src_json["TraceDisplayActivityUncert"] = PhysicalUnits::printToBestActivityUnits(src.traceSrcDisplayActUncertainty,4,!useBq) + trace_src_postfix;
      src_json["TraceDisplayActivityUncert_bq"] = src.traceSrcDisplayActUncertainty / PhysicalUnits::bq;
      src_json["TraceDisplayActivityUncert_kBq"] = src.traceSrcDisplayActUncertainty / PhysicalUnits::kBq;
      src_json["TraceDisplayActivityUncert_MBq"] = src.traceSrcDisplayActUncertainty / PhysicalUnits::MBq;
      src_json["TraceDisplayActivityUncert_GBq"] = src.traceSrcDisplayActUncertainty / PhysicalUnits::GBq;
      src_json["TraceDisplayActivityUncert_ci"]  = src.traceSrcDisplayActUncertainty / PhysicalUnits::ci;
      src_json["TraceDisplayActivityUncert_mCi"] = src.traceSrcDisplayActUncertainty / PhysicalUnits::mCi;
      src_json["TraceDisplayActivityUncert_uCi"] = src.traceSrcDisplayActUncertainty / PhysicalUnits::microCi;
      src_json["TraceDisplayActivityUncert_pCi"] = src.traceSrcDisplayActUncertainty / PhysicalUnits::pCi;
      
    }else
    {
      assert( src.traceSrcDisplayActUncertainty <= 0.0 );
    }
    
    if( src.traceActivityType == GammaInteractionCalc::TraceActivityType::ExponentialDistribution )
    {
      src_json["TraceRelaxationLength"] = PhysicalUnits::printToBestLengthUnits(src.traceRelaxationLength, 4);
      src_json["TraceRelaxationLength_mm"] = src.traceRelaxationLength / PhysicalUnits::mm;
      src_json["TraceRelaxationLength_cm"] = src.traceRelaxationLength / PhysicalUnits::cm;
      src_json["TraceRelaxationLength_m"] = src.traceRelaxationLength / PhysicalUnits::m;
      src_json["TraceRelaxationLength_inch"] = src.traceRelaxationLength / (2.54*PhysicalUnits::cm);
    }
  }//if( src.isTraceSource )
  
  src_json["IsSelfAttenSource"] = src.isSelfAttenSource;
  if( src.isSelfAttenSource )
  {
    src_json["SelfAttenShieldIndex"] = static_cast<int>( src.selfAttenShieldIndex );
    src_json["SelfAttenShieldName"] = src.selfAttenShieldName;
    src_json["SelfAttenIsVariableMassFrac"] = src.isSelfAttenVariableMassFrac;
    src_json["SelfAttenMassFrac"] = src.selfAttenMassFrac;
    if( src.isSelfAttenVariableMassFrac )
      src_json["SelfAttenMassFracUncert"] = src.selfAttenMassFracUncertainty;
  }
}//add_basic_src_details( SourceDetails &src, src_json )
  
  
  
  void add_act_shield_fit_options_to_json( const ShieldingSourceFitCalc::ShieldingSourceFitOptions &options,
                               const double distance,
                               const GammaInteractionCalc::GeometryType geometry,
                               const std::shared_ptr<const DetectorPeakResponse> &drf,
                               nlohmann::json &data )
  {
    const bool fixedGeom = (drf && drf->isFixedGeometry());
    
    data["FixedGeometryDetector"] = fixedGeom;
    
    auto &fit_setup = data["ActShieldFitSetup"];
    if( !fixedGeom )
    {
      fit_setup["Distance"] = PhysicalUnits::printToBestLengthUnits( distance, 3 );
      fit_setup["Distance_mm"] = distance / PhysicalUnits::mm;
      fit_setup["Distance_cm"] = distance / PhysicalUnits::cm;
      fit_setup["Distance_m"] = distance / PhysicalUnits::meter;
      fit_setup["Distance_km"] = distance / (1000.0*PhysicalUnits::meter);
      fit_setup["Distance_inch"] = distance / (2.54*PhysicalUnits::cm);
      fit_setup["Distance_feet"] = distance / (12.0*2.54*PhysicalUnits::cm);
      fit_setup["Geometry"] = GammaInteractionCalc::to_str( geometry );
    }else
    {
      string gem_desc;
      switch( drf->geometryType() )
      {
        case DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic:
        case DetectorPeakResponse::EffGeometryType::FarFieldAbsolute:
          assert( (drf->geometryType() != DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic)
                 && (drf->geometryType() != DetectorPeakResponse::EffGeometryType::FarFieldAbsolute) );
          break;
          
        case DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct:
          gem_desc = "total activity";
          break;
          
        case DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2:
          gem_desc = "activity per square centimeter";
          break;
          
        case DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2:
          gem_desc = "activity per square meter";
          break;
          
        case DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram:
          gem_desc = "activity per gram";
          break;
      }//switch( drf->geometryType() )
      
      fit_setup["FixedGeometryType"] = gem_desc;
    }//if( !fixedGeom ) / else
    
    auto &fit_options = fit_setup["FitOptions"];
    fit_options["InterferenceCorrection"]   = options.multiple_nucs_contribute_to_peaks;
    fit_options["AttenuateForAir"]          = (options.attenuate_for_air && (!drf || !drf->isFixedGeometry()));
    fit_options["DecayDuringMeasurement"]   = options.account_for_decay_during_meas;
    fit_options["MultithreadSelfAttenCalc"] = options.multithread_self_atten;
    fit_options["PhotopeakClusterSigma"]    = options.photopeak_cluster_sigma;
    fit_options["BackgroundPeakSubtract"]   = options.background_peak_subtract;
    fit_options["ElementNuclidesSameAge"]   = options.same_age_isotopes;
    fit_options["AccountForDrfUncert"]      = options.account_for_drf_uncert;
    fit_options["CorrectForCascadeSumming"] = options.correct_for_cascade_summing;
  }//void add_act_shield_fit_options_to_json(...)
  
  
  /** Adds basic information about a peak (energy, fwhm, counts, etc), but not any information
     about gammas that contribute to it, etc
   */
  void add_basic_peak_info( const GammaInteractionCalc::PeakDetail &peak, nlohmann::basic_json<> &peak_json )
  {
    char buffer[64] = { '\0' };
    snprintf( buffer, sizeof(buffer), "%.2f", peak.energy );
    peak_json["Energy"] = string(buffer);
    peak_json["Energy_keV"] = peak.energy;
    
    snprintf( buffer, sizeof(buffer), "%.3f", peak.decayParticleEnergy );
    peak_json["DecayParticleEnergy"] = string(buffer);
    peak_json["DecayParticleEnergy_keV"] = peak.decayParticleEnergy;
    peak_json["AssignedNuclide"] = peak.assignedNuclide;
    peak_json["FWHM"] = peak.fwhm;
    
    
    peak_json["Counts"] = peak.counts;
    peak_json["CountsStr"] = SpecUtils::printCompact(peak.counts, 4);
    peak_json["CountsUncert"] = peak.countsUncert;
    peak_json["CountsUncertStr"] = SpecUtils::printCompact(peak.countsUncert, 4);
    peak_json["Cps"] = peak.cps;
    peak_json["CpsStr"] = SpecUtils::printCompact(peak.cps,4);
    peak_json["CpsUncert"] = peak.cpsUncert;
    peak_json["ShieldAttenuations"] = peak.m_attenuations;
    
    peak_json["AttenuationByShieldingFactor"] = peak.m_totalShieldAttenFactor;
    peak_json["AttenuationByAirFactor"] = peak.m_airAttenFactor;
    peak_json["AttenuationTotalFactor"] = peak.m_totalAttenFactor;
    
    if( peak.backgroundCounts > 0.0 )
    {
      peak_json["BackgroundCounts"] = peak.backgroundCounts;
      peak_json["BackgroundCountsStr"] = SpecUtils::printCompact(peak.backgroundCounts, 4);
      
      if( peak.backgroundCountsUncert > 0.0 )
      {
        peak_json["BackgroundCountsUncert"] = peak.backgroundCountsUncert;
        peak_json["BackgroundCountsUncertStr"] = SpecUtils::printCompact(peak.backgroundCountsUncert, 4);
      }
    }//if( peak.backgroundCounts > 0.0 )
    
    peak_json["SignalCounts"] = peak.observedCounts;
    peak_json["SignalCountsStr"] = SpecUtils::printCompact(peak.observedCounts, 4);
    peak_json["SignalCountsUncert"] = peak.observedUncert;
    peak_json["SignalCountsUncertStr"] = SpecUtils::printCompact(peak.observedUncert, 4);
    
    peak_json["PredictedCounts"] = peak.expectedCounts;
    peak_json["PredictedNumSigmaOff"] = peak.numSigmaOff;
    
    peak_json["ObservedOverPredicted"] = peak.observedOverExpected;
    peak_json["ObservedOverPredictedUncert"] = peak.observedOverExpectedUncert;
    
    peak_json["DetectorSolidAngleFraction"] = peak.detSolidAngle;
    peak_json["DetectorIntrinsicEff"] = peak.detIntrinsicEff;
    peak_json["DetectorEff"] = peak.detEff;

    // Fractional 1-sigma DRF efficiency uncertainty at this energy; present
    //  only when the `account_for_drf_uncert` option was on and the DRF has
    //  uncertainty info (SignalCountsUncert then includes this component).
    if( peak.drfEffFracUncert > 0.0 )
    {
      peak_json["DrfEffFracUncert"] = peak.drfEffFracUncert;
      peak_json["DrfEffFracUncertPercentStr"]
                     = SpecUtils::printCompact( 100.0*peak.drfEffFracUncert, 3 );
    }

    // Validity flag of the detector-efficiency query at this energy/geometry;
    //  present only when the query fell outside the responses validated
    //  regime (near-field, refuse-grade off-axis, shadowed, clamped).
    if( peak.drfEffFlag != DetectorPeakResponse::EffFlag::Ok )
      peak_json["DrfEffFlag"] = DetectorPeakResponse::effFlagName( peak.drfEffFlag );

    // True-coincidence (cascade) summing correction applied to the predicted
    //  counts; present only when the `correct_for_cascade_summing` option was
    //  on and this peaks nuclide has coincident emissions.
    if( peak.cascadeCorrApplied )
    {
      auto &cascade = peak_json["CascadeCorrection"];
      cascade["NetMult"] = peak.cascadeNetMult;
      cascade["NetMultStr"] = SpecUtils::printCompact( peak.cascadeNetMult, 5 );
      cascade["SummingOut"] = peak.cascadeSummingOut;
      cascade["SummingOutStr"] = SpecUtils::printCompact( peak.cascadeSummingOut, 5 );
      cascade["SummingIn"] = peak.cascadeSummingIn;
      cascade["SummingInStr"] = SpecUtils::printCompact( peak.cascadeSummingIn, 5 );
    }//if( peak.cascadeCorrApplied )
  }//add_basic_peak_info( const PeakDetail &peak, nlohmann::basic_json<> &peak_json )
   
  
  void add_gamma_info_for_peak( const GammaInteractionCalc::PeakDetailSrc &ps,
                    const GammaInteractionCalc::SourceDetails * const src,
                    const std::shared_ptr<const DetectorPeakResponse> &drf,
                    const bool useBq,
                    const std::vector<GammaInteractionCalc::ShieldingDetails> * const shield_details,
                    nlohmann::basic_json<> &gamma_json )
  {
    assert( ps.nuclide );
    
    char buffer[64] = { '\0' };
    snprintf( buffer, sizeof(buffer), "%.3f", ps.energy );
    
    gamma_json["Nuclide"] = ps.nuclide ? ps.nuclide->symbol : string("null");
    gamma_json["Energy"] = string(buffer);
    gamma_json["Energy_keV"] = ps.energy;
    
    gamma_json["BranchingRatio"] = ps.br;
    gamma_json["BranchingRatioStr"] = SpecUtils::printCompact( ps.br, 5 );
    
    gamma_json["PredictedCounts"] = ps.modelContribToPeak;
    gamma_json["PredictedCountsStr"] = SpecUtils::printCompact(ps.modelContribToPeak, 4);
    
    gamma_json["SourcePhotonsCps"] = ps.cpsAtSource;
    gamma_json["SourcePhotonsCpsStr"] = SpecUtils::printCompact( ps.cpsAtSource, 5 );
    
    gamma_json["SourcePhotons"] = ps.countsAtSource;
    gamma_json["SourcePhotonsStr"] = SpecUtils::printCompact( ps.countsAtSource, 5 );
    
    gamma_json["HasDecayCorrection"] = (ps.decayCorrection > 0.0);
    if( ps.decayCorrection > 0.0 )
    {
      gamma_json["DecayCorrection"] = ps.decayCorrection;
      gamma_json["DecayCorrectionStr"] = SpecUtils::printCompact(ps.decayCorrection, 4);
    }
    
    // The other information in `PeakDetailSrc` should be repeat of
    //  `GammaInteractionCalc::SourceDetails`.  We _could_ access all this information
    //  from the templating code, since we are in a loop over `SourceDetails`, but
    //  to make things easier on people, we'll just re-include it here.
    if( src )
      add_basic_src_details( *src, drf, useBq, shield_details, gamma_json );
  };//add_gamma_info_for_peak(...)
  
  
  void shield_src_fit_results_to_json( const ShieldingSourceFitCalc::ModelFitResults &results,
                                      const std::shared_ptr<const DetectorPeakResponse> &drf,
                                      const bool useBq,
                                      nlohmann::json &data )
  {
    const DetectorPeakResponse::EffGeometryType eff_type = drf ? drf->geometryType()
                                                : DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic;
    
    const string act_postfix = DetectorPeakResponse::det_eff_geom_type_postfix(eff_type);
    
    add_act_shield_fit_options_to_json( results.options, results.distance, results.geometry, drf, data );
    
    // Now put in results
    data["FitChi2"] = results.chi2;
    data["EstimatedDistanceToMinimum"] = results.edm;
    data["NumberFcnCalls"] = results.num_fcn_calls;
    data["NumDof"] = results.numDOF;
    
    auto &fit_pars = data["RawFitParameter"];
    fit_pars["Values"] = results.paramValues;
    fit_pars["Errors"] = results.paramErrors;
    
    if( !results.errormsgs.empty() )
      data["ErrorMessages"] = results.errormsgs;

    // Non-fatal fit warnings (poor average deviation, x-ray peaks, supplemental-info failures, ...).
    //  Templates render these under {% if HasWarnings %} / {{ Warnings }}; the GUI and batch both
    //  reach this function, so populating them here covers every report.
    data["HasWarnings"] = !results.warnings.empty();
    data["Warnings"] = results.warnings;

    // The automatic pull-trend interpretation (too much/little shielding, wrong effective atomic
    //  number, ...) as a dedicated field, so a template can place it under the chi chart (like the
    //  GUI) and/or in the warnings area.  The report always reflects a completed fit, so it's valid.
    //  Always present (with HasConclusion) so template access is unconditionally safe.
    {
      auto &tj = data["PullTrend"];
      tj["HasConclusion"] = false;
      tj["Message"] = "";
      if( results.pull_trend )
      {
        const ShieldSourcePullTrend::TrendResult &trend = *results.pull_trend;
        const string trend_msg = ShieldSourcePullTrend::conclusion_report_text( trend );
        tj["Message"] = trend_msg;
        tj["HasConclusion"] = !trend_msg.empty();
        tj["SlopeT"] = trend.slopeT;
        tj["CurvatureT"] = trend.curvatureT;
        tj["ReducedChi2"] = trend.redChi2;
        tj["NumPeaksUsed"] = static_cast<int>( trend.numPointsUsed );
        tj["AtomicNumberDiscriminationPower"] = trend.anDiscriminationT;
        tj["AtomicNumberDiscriminable"] = trend.anDiscriminable;
      }
    }

    int num_sources = 0;
    bool hasAnyTraceSrc = false, hasAnyVolumetricSrc = false, hasFitAnyAge = false;
    if( results.source_calc_details )
    {
      for( const GammaInteractionCalc::SourceDetails &src : *results.source_calc_details )
      {
        data["Sources"].push_back( {} );
        auto &src_json = data["Sources"].back();
       
        ++num_sources;
        hasAnyTraceSrc |= src.isTraceSource;
        hasAnyVolumetricSrc |= src.isSelfAttenSource;
        hasFitAnyAge |= src.ageIsFit;
        
        add_basic_src_details( src, drf, useBq, results.shield_calc_details.get(), src_json );
        
        // We wont put peaks into this struct, but instead when we make the JSON, we'll
        //  insert peaks from `PeakDetail` as `PeakDetailSrc::nuclide` match this nuclide.
        if( results.peak_calc_details )
        {
          for( const GammaInteractionCalc::PeakDetail &peak : *results.peak_calc_details )
          {
            bool src_contribs_to_peak = false;
            for( const GammaInteractionCalc::PeakDetailSrc &ps :  peak.m_sources )
            {
              src_contribs_to_peak = (ps.nuclide == src.nuclide);
              if( src_contribs_to_peak )
                break;
            }
            if( !src_contribs_to_peak )
              continue;
            
            src_json["PeaksThisNucContributesTo"].push_back( {} );
            nlohmann::basic_json<> &peak_json = src_json["PeaksThisNucContributesTo"].back();
            
            add_basic_peak_info( peak, peak_json );
            
            for( const GammaInteractionCalc::PeakDetailSrc &ps :  peak.m_sources )
            {
              if( ps.nuclide != src.nuclide )
                continue;
              
              peak_json["ThisNucsGammasForPeak"].push_back( {} );
              nlohmann::basic_json<> &gamma_json = peak_json["ThisNucsGammasForPeak"].back();
              
              add_gamma_info_for_peak( ps, &src, drf, useBq, results.shield_calc_details.get(), gamma_json );
            }
          }//for( loop over results.peak_calc_details )
        }//if( results.peak_calc_details )
      }//for( loop over results.source_calc_details )
    }//if( results.source_calc_details )
    
    data["NumSources"] = num_sources;
    
    data["HasTraceSource"] = hasAnyTraceSrc;
    data["HasSelfAttenSource"] = hasAnyVolumetricSrc;
    data["HasVolumetricSource"] = (hasAnyTraceSrc || hasAnyVolumetricSrc);
    data["AnySourceAgeFit"] = hasFitAnyAge;
    
    auto &shieldings_json = data["Shieldings"];
    shieldings_json["Geometry"] = GammaInteractionCalc::to_str(results.geometry);
    
    bool fitAnyShielding = false;
    if( !results.shield_calc_details )
    {
      shieldings_json["NumberShieldings"] = 0;
    }else
    {
      const vector<GammaInteractionCalc::ShieldingDetails> &shield_details
                                                            = *results.shield_calc_details;
      shieldings_json["NumberShieldings"] = static_cast<int>( shield_details.size() );
      switch( results.geometry )
      {
        case GammaInteractionCalc::GeometryType::Spherical:
          shieldings_json["DimensionMeanings"] = vector<string>{"Radius"};
          shieldings_json["NumDimensions"] = 1;
          break;
          
        case GammaInteractionCalc::GeometryType::CylinderEndOn:
        case GammaInteractionCalc::GeometryType::CylinderSideOn:
          shieldings_json["DimensionMeanings"] = vector<string>{"Radius", "Length"};
          shieldings_json["NumDimensions"] = 2;
          break;
          
        case GammaInteractionCalc::GeometryType::Rectangular:
          shieldings_json["DimensionMeanings"] = vector<string>{"Width", "Height", "Depth"};
          shieldings_json["NumDimensions"] = 3;
          break;
          
        case GammaInteractionCalc::GeometryType::NumGeometryType:
          assert( 0 );
          break;
      }//switch( m_geometry )
      
      for( size_t shield_num = 0; shield_num < shield_details.size(); ++shield_num )
      {
        const GammaInteractionCalc::ShieldingDetails &shield = shield_details[shield_num];
        
        shieldings_json["Shields"].push_back({});
        auto &shield_json = shieldings_json["Shields"].back();
        
        shield_json["Name"] = shield.m_name;
        shield_json["ShieldingNumber"] = static_cast<int>( shield_num );
        shield_json["IsGeneric"] = shield.m_is_generic;
        if( !shield.m_is_generic )
        {
          shield_json["Formula"] = shield.m_chemical_formula;
          const double density = shield.m_density * PhysicalUnits::cm3 / PhysicalUnits::g;
          shield_json["Density_gPerCm3"] = density;
          shield_json["DensityStr"] = SpecUtils::printCompact(density, 5) + " g/cm3";
        }
        
        if( shield.m_an > 0 )
          shield_json["AN"] = SpecUtils::printCompact(shield.m_an, 4);
        if( shield.m_ad > 0 )
          shield_json["AD"] = SpecUtils::printCompact(shield.m_ad*PhysicalUnits::cm2/PhysicalUnits::g, 4);
        
        if( shield.m_is_generic )
        {
          const bool fitAn = shield.m_fit_dimension[0];
          const bool fitAD = shield.m_fit_dimension[1];
          shield_json["FitAN"] = fitAn;
          shield_json["FitAD"] = fitAD;
          fitAnyShielding |= fitAn;
          fitAnyShielding |= fitAD;
        }
        
        const vector<bool> fit_dim( shield.m_fit_dimension, shield.m_fit_dimension + shield.m_num_dimensions );
        shield_json["DimensionIsFit"] = fit_dim;
        
        shield_json["NumDimensions"] = static_cast<int>( shield.m_num_dimensions );
        shield_json["Geometry"] = GammaInteractionCalc::to_str(shield.m_geometry);
        
        shield_json["Thickness"] = PhysicalUnits::printToBestLengthUnits(shield.m_thickness,3);
        shield_json["Thickness_mm"] = shield.m_thickness / PhysicalUnits::mm;
        shield_json["Thickness_cm"] = shield.m_thickness / PhysicalUnits::cm;
        shield_json["Thickness_m"] = shield.m_thickness / PhysicalUnits::meter;
        shield_json["Thickness_inch"] = shield.m_thickness / (2.54*PhysicalUnits::cm);
        shield_json["Thickness_feet"] = shield.m_thickness / (12.0*2.54*PhysicalUnits::cm);
        shield_json["VolumeCm3"] = shield.m_volume / PhysicalUnits::cm3;
        shield_json["VolumeUncertCm3"] = shield.m_volume_uncert / PhysicalUnits::cm3;
        
        shield_json["InnerRadius"] = PhysicalUnits::printToBestLengthUnits(shield.m_inner_rad, 3);
        shield_json["OuterRadius"] = PhysicalUnits::printToBestLengthUnits(shield.m_inner_rad + shield.m_thickness, 3);
        
        const vector<double> inner_dims{ shield.m_inner_dimensions, shield.m_inner_dimensions + shield.m_num_dimensions };
        const vector<double> outer_dims{ shield.m_outer_dimensions, shield.m_outer_dimensions + shield.m_num_dimensions };
        vector<double> thicknesses( shield.m_num_dimensions );
        const vector<double> dim_uncerts( shield.m_dimension_uncert, shield.m_dimension_uncert + shield.m_num_dimensions );
        
        vector<double> inner_dims_mm = inner_dims;
        vector<double> outer_dims_mm = outer_dims;
        vector<double> thicknesses_mm = thicknesses;
        vector<double> dim_uncerts_mm = dim_uncerts;
        
        vector<double> inner_dims_cm = inner_dims, outer_dims_cm = outer_dims;
        vector<double> thicknesses_cm = thicknesses, dim_uncerts_cm = dim_uncerts;
        
        vector<double> inner_dims_m = inner_dims, outer_dims_m = outer_dims;
        vector<double> thicknesses_m = thicknesses, dim_uncerts_m = dim_uncerts;
        
        vector<double> inner_dims_inch = inner_dims, outer_dims_inch = outer_dims;
        vector<double> thicknesses_inch = thicknesses, dim_uncerts_inch = dim_uncerts;
        
        
        vector<string> inner_dims_strs( shield.m_num_dimensions );
        vector<string> outer_dims_strs( shield.m_num_dimensions );
        vector<string> thicknesses_strs( shield.m_num_dimensions );
        vector<string> dim_uncerts_strs( shield.m_num_dimensions );
        
        for( unsigned int dim = 0; dim < shield.m_num_dimensions; ++dim )
        {
          fitAnyShielding |= fit_dim[dim];
          thicknesses[dim] = outer_dims[dim] - inner_dims[dim];
          thicknesses_strs[dim] = PhysicalUnits::printToBestLengthUnits( thicknesses[dim], 5 );
          inner_dims_strs[dim] = PhysicalUnits::printToBestLengthUnits( inner_dims[dim], 5 );
          outer_dims_strs[dim] = PhysicalUnits::printToBestLengthUnits( outer_dims[dim], 5 );
          if( fit_dim[dim] )
            dim_uncerts_strs[dim] = PhysicalUnits::printToBestLengthUnits( dim_uncerts[dim], 5 );
          
          inner_dims_mm[dim] /= PhysicalUnits::mm;
          outer_dims_mm[dim] /= PhysicalUnits::mm;
          thicknesses_mm[dim] /= PhysicalUnits::mm;
          dim_uncerts_mm[dim] /= PhysicalUnits::mm;
          
          inner_dims_cm[dim] /= PhysicalUnits::cm;
          outer_dims_cm[dim] /= PhysicalUnits::cm;
          thicknesses_cm[dim] /= PhysicalUnits::cm;
          dim_uncerts_cm[dim] /= PhysicalUnits::cm;
          
          inner_dims_m[dim] /= PhysicalUnits::m;
          outer_dims_m[dim] /= PhysicalUnits::m;
          thicknesses_m[dim] /= PhysicalUnits::m;
          dim_uncerts_m[dim] /= PhysicalUnits::m;
          
          inner_dims_inch[dim] /= (2.54*PhysicalUnits::cm);
          outer_dims_inch[dim] /= (2.54*PhysicalUnits::cm);
          thicknesses_inch[dim] /= (2.54*PhysicalUnits::cm);
          dim_uncerts_inch[dim] /= (2.54*PhysicalUnits::cm);
        }//for( unsigned int dim = 0; dim < shield.m_num_dimensions; ++dim )
        
        shield_json["ThicknessesUncerts"] = dim_uncerts_strs;
        shield_json["InnerDims"]          = inner_dims_strs;
        shield_json["OuterDims"]          = outer_dims_strs;
        shield_json["Thicknesses"]        = thicknesses_strs;
        
        shield_json["InnerDims_mm"] = inner_dims_mm;
        shield_json["OuterDims_mm"] = outer_dims_mm;
        shield_json["Thicknesses_mm"] = thicknesses_mm;
        shield_json["ThicknessesUncerts_mm"] = dim_uncerts_mm;
        
        shield_json["InnerDims_cm"] = inner_dims_cm;
        shield_json["OuterDims_cm"] = outer_dims_cm;
        shield_json["Thicknesses_cm"] = thicknesses_cm;
        shield_json["ThicknessesUncerts_cm"] = dim_uncerts_cm;
        
        shield_json["InnerDims_m"] = inner_dims_m;
        shield_json["OuterDims_m"] = outer_dims_m;
        shield_json["Thicknesses_m"] = thicknesses_m;
        shield_json["ThicknessesUncerts_m"] = dim_uncerts_m;
        
        shield_json["InnerDims_inch"] = inner_dims_inch;
        shield_json["OuterDims_inch"] = outer_dims_inch;
        shield_json["Thicknesses_inch"] = thicknesses_inch;
        shield_json["ThicknessesUncerts_inch"] = dim_uncerts_inch;
        
        bool fittingAnyMassFrac = false;
        for( const GammaInteractionCalc::ShieldingDetails::SelfAttenComponent &comp : shield.m_mass_fractions )
        {
          fittingAnyMassFrac |= comp.m_is_fit;
          
          shield_json["SelfAttenSources"].push_back( {} );
          auto &self_atten_json = shield_json["SelfAttenSources"].back();
        
          assert( comp.m_nuclide );
          self_atten_json["Nuclide"] = comp.m_nuclide ? comp.m_nuclide->symbol : string("null");
          self_atten_json["IsFittingMassFraction"] = comp.m_is_fit;
          self_atten_json["MassFraction"] = comp.m_mass_frac;
          self_atten_json["MassFractionStr"] = SpecUtils::printCompact(comp.m_mass_frac, 5);
          if( comp.m_is_fit )
          {
            self_atten_json["MassFractionUncert"] = comp.m_mass_frac_uncert;
            self_atten_json["MassFractionUncertStr"] = SpecUtils::printCompact(comp.m_mass_frac_uncert, 5);
          }
          
          assert( results.source_calc_details );
          if( results.source_calc_details )
          {
            const vector<GammaInteractionCalc::SourceDetails> &srcs = *results.source_calc_details;
            auto pos = std::find_if( begin(srcs), end(srcs),
              [&comp]( const GammaInteractionCalc::SourceDetails &src ){
                return (src.nuclide == comp.m_nuclide);
            } );
            assert( pos != end(srcs) );
            if( pos != end(srcs) )
              add_basic_src_details( *pos, drf, useBq, results.shield_calc_details.get(), self_atten_json );
          }//if( results.source_calc_details )
        }//for( SelfAttenComponent & comp: shield.m_mass_fractions )
        
        shield_json["FitAnyMassFraction"] = fittingAnyMassFrac;
        
        
        for( const GammaInteractionCalc::ShieldingDetails::TraceSrcDetail &trace : shield.m_trace_sources )
        {
          shield_json["TraceSources"].push_back( {} );
          auto &trace_src_json = shield_json["TraceSources"].back();
          
          assert( trace.m_nuclide );
          trace_src_json["Nuclide"] = trace.m_nuclide ? trace.m_nuclide->symbol : string("null");
          trace_src_json["TraceSourceType"] = GammaInteractionCalc::to_str( trace.m_trace_type );
          if( trace.m_is_exp_dist )
          {
            trace_src_json["RelaxationLength"] = PhysicalUnits::printToBestLengthUnits( trace.m_relaxation_length, 4 );
            trace_src_json["RelaxationLength_mm"] = trace.m_relaxation_length / PhysicalUnits::mm;
            trace_src_json["RelaxationLength_cm"] = trace.m_relaxation_length / PhysicalUnits::cm;
            trace_src_json["RelaxationLength_m"] = trace.m_relaxation_length / PhysicalUnits::m;
            trace_src_json["RelaxationLength_inch"] = trace.m_relaxation_length / (2.54*PhysicalUnits::cm);
          }
          
          assert( results.source_calc_details );
          if( results.source_calc_details )
          {
            const vector<GammaInteractionCalc::SourceDetails> &srcs = *results.source_calc_details;
            auto pos = std::find_if( begin(srcs), end(srcs),
              [&trace]( const GammaInteractionCalc::SourceDetails &src ){
                return (src.nuclide == trace.m_nuclide);
            } );
            assert( pos != end(srcs) );
            if( pos != end(srcs) )
              add_basic_src_details( *pos, drf, useBq, results.shield_calc_details.get(), trace_src_json );
          }//if( results.source_calc_details )
        }//for( const TraceSrcDetail &src : shield.m_trace_sources )
        
      }//for( const GammaInteractionCalc::ShieldingDetails &shield : *results.shield_calc_details )
    }//if( results.shield_calc_details )
    
    data["AnyShieldingFit"] = fitAnyShielding;
    
    
    // Add Peak Details to JSON
    if( results.peak_calc_details )
    {
      auto &peaks_json = data["PeaksUsedForActivityFitting"];
      
      for( const GammaInteractionCalc::PeakDetail &peak : *results.peak_calc_details )
      {
        peaks_json["Peaks"].push_back( {} );
        nlohmann::basic_json<> &peak_json = peaks_json["Peaks"].back();
        add_basic_peak_info( peak, peak_json );
        
        for( const GammaInteractionCalc::PeakDetailSrc &pksrc : peak.m_sources )
        {
          assert( pksrc.nuclide && results.source_calc_details );
          if( !pksrc.nuclide || !results.source_calc_details )
            continue;
          
          const vector<GammaInteractionCalc::SourceDetails> &srcs = *results.source_calc_details;
          auto pos = std::find_if( begin(srcs), end(srcs),
            [&peak]( const GammaInteractionCalc::SourceDetails &src ){
              return src.nuclide && (src.nuclide->symbol == peak.assignedNuclide);
          } );
          assert( pos != end(srcs) );
          if( pos == end(srcs) )
            continue;
          
          
          peak_json["Sources"].push_back( {} );
          nlohmann::basic_json<> &src_json = peak_json["Sources"].back();
          
          add_basic_src_details( *pos, drf, useBq, results.shield_calc_details.get(), src_json );
          
          src_json["HasDecayCorrection"] = (pksrc.decayCorrection > 0.0);
          if( pksrc.decayCorrection > 0.0 )
          {
            peak_json["DecayCorrection"] = pksrc.decayCorrection;
            peak_json["DecayCorrectionStr"] = SpecUtils::printCompact(pksrc.decayCorrection, 4);
          }
          
          char buffer[64] = { '\0' };
          snprintf( buffer, sizeof(buffer), "%.2f", pksrc.energy );
          
          src_json["Energy"] = string(buffer);
          src_json["Energy_keV"] = pksrc.energy;
          
          src_json["SourcePhotonsCps"] = pksrc.cpsAtSource;
          src_json["SourcePhotonsCpsStr"] = SpecUtils::printCompact(pksrc.cpsAtSource, 4);
          
          src_json["BranchingRatio"] = pksrc.br;
          src_json["BranchingRatioStr"] = SpecUtils::printCompact( pksrc.br, 5 );
          
          snprintf( buffer, sizeof(buffer), "%.2f", pksrc.countsAtSource );
          src_json["SourcePhotons"] = pksrc.countsAtSource;
          src_json["SourcePhotonsStr"] = string(buffer);
          
          src_json["PredictedCounts"] = pksrc.modelContribToPeak;
          src_json["PredictedCountsStr"] = SpecUtils::printCompact(pksrc.modelContribToPeak, 4);
        }//for( const GammaInteractionCalc::PeakDetailSrc &pksrc : peak.m_sources )
        
        // We could loop over the volumetric sources, and add info, but...
        /*
        struct VolumeSrc
        {
          bool trace;
          double integral;
          double volume;
          double averageEfficiencyPerSourceGamma;
          double srcVolumetricActivity;
          bool inSituExponential;
          double inSituRelaxationLength;
          double detIntrinsicEff;
          std::string sourceName;
        };//struct VolumeSrc
        
        std::vector<VolumeSrc> m_volumetric_srcs;
         */
      }//for( loop over results.peak_calc_details )
    }//if( results.peak_calc_details )
                                      
    auto add_detector_to_json = []( nlohmann::json &data, const shared_ptr<const DetectorPeakResponse> &drf ){
      if( !drf || !drf->isValid() )
        return;
                                        
      auto &drf_obj = data["Detector"];
      drf_obj["Name"] = drf->name();
      drf_obj["Description"] = drf->description();
      drf_obj["Diameter_mm"] = drf->detectorDiameter() / PhysicalUnits::mm;
      drf_obj["Diameter_cm"] = drf->detectorDiameter() / PhysicalUnits::cm;
      drf_obj["Diameter_m"] = drf->detectorDiameter() / PhysicalUnits::m;
      drf_obj["Diameter_inch"] = drf->detectorDiameter() / (2.54*PhysicalUnits::cm);
      drf_obj["Radius_mm"] = 0.5*drf->detectorDiameter() / PhysicalUnits::mm;
      drf_obj["Radius_cm"] = 0.5*drf->detectorDiameter() / PhysicalUnits::cm;
      drf_obj["Radius_m"] = 0.5*drf->detectorDiameter() / PhysicalUnits::m;
      drf_obj["Radius_inch"] = 0.5*drf->detectorDiameter() / (2.54*PhysicalUnits::cm);
      drf_obj["Diameter"] = PhysicalUnits::printToBestLengthUnits( drf->detectorDiameter(), 3 );
      drf_obj["Radius"] = PhysicalUnits::printToBestLengthUnits( 0.5*drf->detectorDiameter(), 3 );
      if( drf->detectorSetback() > 0.0 )
      {
        drf_obj["Setback_cm"] = drf->detectorSetback() / PhysicalUnits::cm;
        drf_obj["Setback"] = PhysicalUnits::printToBestLengthUnits( drf->detectorSetback(), 3 );
      }
      drf_obj["FixedGeometry"] = drf->isFixedGeometry();
    };
             
      if( drf )
        add_detector_to_json( data, drf );
    
    if( results.peak_comparisons )
    {
      // This info is already in the Peaks info, but we'll put in anyways
      auto &energy_obj = data["PeakToModelComparison"];

      for( const GammaInteractionCalc::PeakResultPlotInfo &p : *results.peak_comparisons )
      {
        energy_obj["UsedPeaks"].push_back( {} );
        nlohmann::basic_json<> &p_json = energy_obj["UsedPeaks"].back();

        char buffer[64] = { '\0' };
        snprintf( buffer, sizeof(buffer), "%.2f", p.energy );
        p_json["Energy"] = string(buffer);
        p_json["NumSigmaOff"] = SpecUtils::printCompact( p.numSigmaOff, 6 ); //`(observed_counts - expected_counts) / observed_uncertainty`
        p_json["ObservedOverExpected"] = SpecUtils::printCompact( p.observedOverExpected, 6 );
        p_json["ObservedOverExpectedUncert"] = SpecUtils::printCompact( p.observedOverExpectedUncert, 6 );
      }
    }//if( results.peak_comparisons )

    // The per-peak detection limit checks, and the activities implied by peaks that were fit but
    //  not used by the model.  Done here so the GUIs fit log and the batch reports both get it.
    add_supplemental_peak_info_to_json( data, results.supplemental_peak_info, drf, useBq );

    // Add ShieldingSourceFitPlot chart data
    const string plot_json_str = ShieldingSourceFitPlot::jsonForData( results );
    try
    {
      data["ShieldingSourceFitPlotData"] = nlohmann::json::parse( plot_json_str );
    }catch( std::exception &e )
    {
      cerr << "Error parsing ShieldingSourceFitPlot JSON: " << e.what() << endl;
      data["ShieldingSourceFitPlotData"] = nlohmann::json::object();
    }

  }//void shield_src_fit_results_to_json()
  
  
  void add_hist_to_json( nlohmann::json &spec_obj,
                        const bool is_background,
                        const double display_scale_factor,
                       const shared_ptr<const SpecUtils::Measurement> &spec_ptr,
                       const shared_ptr<const SpecMeas> &spec_file,
                       const std::set<int> &sample_numbers,
                       const string &filename,
                        const deque<std::shared_ptr<const PeakDef>> * const peak_fit )
  {
    if( !spec_ptr )
      return;
    
    const SpecUtils::Measurement &spec = *spec_ptr;
    
    const double lt = spec.live_time();
    const double rt = spec.real_time();
    const string lt_str = PhysicalUnits::printToBestTimeUnits(lt,3);
    const string rt_str = PhysicalUnits::printToBestTimeUnits(rt,3);
    const vector<int> sample_nums( begin(sample_numbers), end(sample_numbers) );
    
    D3SpectrumExport::D3SpectrumOptions spec_json_options;
    if( peak_fit )
    {
      const deque<std::shared_ptr<const PeakDef>> &fore_peaks = *peak_fit;
      const vector<shared_ptr<const PeakDef> > inpeaks( begin(fore_peaks), end(fore_peaks) );
      spec_json_options.peaks_json = PeakDef::peak_json( inpeaks, spec_ptr, Wt::WColor(), 255 );
    }//if( fit_results.m_peak_fit_results )

    //spec_json_options.line_color = "rgb(0,0,0)"; //black
    spec_json_options.peak_color = "rgba(0,51,255,0.6)";
    spec_json_options.title = "";
    spec_json_options.display_scale_factor = ((display_scale_factor > 0.0) && !IsNan(display_scale_factor) && !IsInf(display_scale_factor))
                                               ? display_scale_factor
                                               : 1.0;
    spec_json_options.spectrum_type = is_background ? SpecUtils::SpectrumType::Background : SpecUtils::SpectrumType::Foreground;
    
    //We will only have foreground or background on this spectrum, so even if we say
    //  background ID is 2, instead of a negative number, things should be fine
    const size_t specID = is_background ? 2 : 1;
    const int backID = is_background ? -1 : 2;
    
    stringstream spec_strm;
    D3SpectrumExport::write_spectrum_data_js( spec_strm, spec, spec_json_options, specID, backID );
    const string spectrum_json_str = spec_strm.str();
    
    nlohmann::json spectrum_json;
    try
    {
      if( !spectrum_json_str.empty() )
        spectrum_json = nlohmann::json::parse( spectrum_json_str );
    }catch( std::exception &e )
    {
      cerr << "Failed to parse spectrum JSON: " << e.what()
      << "\n\nJSON: " << spectrum_json_str << endl << endl;
      assert( 0 );
    }
    
    spec_obj["LiveTime"] = lt_str;
    spec_obj["RealTime"] = rt_str;
    spec_obj["LiveTime_s"] = lt;
    spec_obj["RealTime_s"] = rt;
    spec_obj["DeadTime"] = PhysicalUnits::printToBestTimeUnits(rt - lt);
    spec_obj["DeadTime_s"] = (rt - lt)/PhysicalUnits::second;
    spec_obj["DeadTime_percent"] = 100.0*(rt - lt) / rt; //If rt is zero, will create a NaN, and then the JSON will have a `null` for the value
    spec_obj["StartTime"] = SpecUtils::to_extended_iso_string( spec.start_time() );
    spec_obj["StartTime_iso"] = SpecUtils::to_iso_string( spec.start_time() );
    spec_obj["StartTime_vax"] = SpecUtils::to_vax_string( spec.start_time() );
    spec_obj["StartTimeIsValid"] = !SpecUtils::is_special( spec.start_time() );
    spec_obj["LowerSpectrumEnergy"] = spec.gamma_channel_lower(0);
    spec_obj["UpperSpectrumEnergy"] = spec.gamma_channel_upper(spec.num_gamma_channels() - 1);
    spec_obj["NumberChannels"] = (int)spec.num_gamma_channels();
    spec_obj["GammaSum"] = spec.gamma_count_sum();
    spec_obj["GammaCps"] = (spec.live_time() > 0.0) ? (spec.gamma_count_sum() / spec.live_time()) : -1.0;
    spec_obj["SampleNumbers"] = sample_nums;
    //detector names?
    spec_obj["Filename"] = filename;
    spec_obj["spectrum"] = spectrum_json;
    spec_obj["HasGps"] = spec.has_gps_info();
    if( spec.has_gps_info() )
    {
      spec_obj["Longitude"] = spec.longitude();
      spec_obj["Latitude"] = spec.latitude();
    }
    spec_obj["HasNeutrons"] = spec.contained_neutron();
    if( spec.contained_neutron() )
    {
      spec_obj["NeutronCounts"] = spec.neutron_counts_sum();
      spec_obj["NeutronLiveTime"] = spec.neutron_live_time();
      spec_obj["NeutronCps"] = spec.neutron_counts_sum() / spec.neutron_live_time();
    }
    
    if( !spec.title().empty() )
      spec_obj["SpectrumTitle"] = spec.title();
    
    //The measured ambient radiation dose equivalent rate value, in microsieverts per hour (μSv/h).
    if( spec.dose_rate() >= 0.0 )
      spec_obj["DoseRate_uSvPerHour"] = spec.dose_rate();
    
    //The measured radiation exposure rate value, in milliroentgen per hour (mR/h).
    if( spec.exposure_rate() >= 0.0 )
      spec_obj["ExposureRate_mRPerHour"] = spec.exposure_rate();
    
    if( !spec.detector_type().empty() )
      spec_obj["DetectorTypeDesc"] = spec.detector_type();
    
    if( !spec.detector_name().empty() )
      spec_obj["DetectorName"] = spec.detector_name();
    
    //SourceType source_type() const;
    
    if( !spec.remarks().empty() )
      spec_obj["SpectrumRemarks"] = spec.remarks();
    
    if( !spec.parse_warnings().empty() )
      spec_obj["SpectrumParseWarnings"] = spec.parse_warnings();
    
    spec_obj["IsDerivedData"] = (spec.derived_data_properties() != 0);
    
    spec_obj["DerivedIoiSum"] = static_cast<bool>(spec.derived_data_properties()
                                                  & static_cast<uint32_t>(SpecUtils::Measurement::DerivedDataProperties::ItemOfInterestSum));
    spec_obj["DerivedUsedForRidAnalysis"] = static_cast<bool>(spec.derived_data_properties()
                                                              & static_cast<uint32_t>(SpecUtils::Measurement::DerivedDataProperties::UsedForAnalysis));
    spec_obj["DerivedProcessedFurther"] = static_cast<bool>(spec.derived_data_properties()
                                                            & static_cast<uint32_t>(SpecUtils::Measurement::DerivedDataProperties::ProcessedFurther));
    spec_obj["DerivedBackgroundSub"] = static_cast<bool>(spec.derived_data_properties()
                                                         & static_cast<uint32_t>(SpecUtils::Measurement::DerivedDataProperties::BackgroundSubtracted));
    spec_obj["DerivedIsBackground"] = static_cast<bool>(spec.derived_data_properties()
                                                        & static_cast<uint32_t>(SpecUtils::Measurement::DerivedDataProperties::IsBackground));
    
    add_energy_cal_json( spec_obj["EnergyCal"], spec.energy_calibration() );
    
    if( spec_file )
    {
      spec_obj["InstrumentModel"] = spec_file->instrument_model();
      spec_obj["SerialNumber"] = spec_file->instrument_id();
      spec_obj["Manufacturer"] = spec_file->manufacturer();
      spec_obj["InstrumentType"] = spec_file->instrument_type();
      spec_obj["DetectorType"] = SpecUtils::detectorTypeToString( spec_file->detector_type() );
      spec_obj["NumberRecordsInFile"] = static_cast<int>( spec_file->num_measurements() );
      spec_obj["RemarksInFile"] = spec_file->remarks();
      spec_obj["ParseWarningsForFile"] = spec_file->parse_warnings();
      
      spec_obj["HasInstrumentRid"] = !!spec_file->detectors_analysis();
      if( spec_file->detectors_analysis() )
        spec_obj["InstrumentRidSummary"] = riidAnaSummary( spec_file );
    }//if( spec_file )
  }//add_hist_to_json(...)
  
  
  void add_energy_cal_json( nlohmann::basic_json<> &cal_obj,
                           const std::shared_ptr<const SpecUtils::EnergyCalibration> &cal )
  {
    cal_obj["NumChannels"] = cal ? static_cast<int>(cal->num_channels()) : 0;
    const SpecUtils::EnergyCalType cal_type = cal ? cal->type() : SpecUtils::EnergyCalType::InvalidEquationType;
    switch( cal_type )
    {
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
        cal_obj["Type"] = "Polynomial";
        break;
        
      case SpecUtils::EnergyCalType::FullRangeFraction:
        cal_obj["Type"] = "FullRangeFraction";
        break;
        
      case SpecUtils::EnergyCalType::LowerChannelEdge:
        cal_obj["Type"] = "LowerChannelEdge";
        break;
        
      case SpecUtils::EnergyCalType::InvalidEquationType:
        cal_obj["Type"] = "Invalid";
        break;
    }//switch( cal->type() )
    
    if( cal_type == SpecUtils::EnergyCalType::InvalidEquationType )
      return;
    
    cal_obj["LowerEnergy"] = cal->lower_energy();
    cal_obj["UpperEnergy"] = cal->upper_energy();
    
    if( cal->type() != SpecUtils::EnergyCalType::LowerChannelEdge )
    {
      cal_obj["Coefficients"] = vector<double>{ begin(cal->coefficients()), end(cal->coefficients()) };
      
      if( !cal->deviation_pairs().empty() )
      {
        vector<vector<double>> pairs;
        for( const pair<float,float> &p : cal->deviation_pairs() )
          pairs.push_back( { static_cast<double>(p.first), static_cast<double>(p.second) } );
        cal_obj["DeviationPairs"] = pairs;
      }//if( dev pairs )
    }//if( polynomial or FRF )
  }//void add_energy_cal_json(...)
  
  
  void add_activity_fit_options_to_json( nlohmann::json &data,
                                      const BatchActivity::BatchActivityFitOptions &options )
  {
    add_peak_fit_options_to_json( data, options );
    
    auto &options_obj = data["PeakFitOptions"];
    
    const bool overiding_dist = options.distance_override.has_value();
    options_obj["IsSpecifyingDistance"] = overiding_dist;
    if( overiding_dist )
    {
      const double dist = options.distance_override.value();
      options_obj["SpecifiedDistance_m"] = dist / PhysicalUnits::m;
      options_obj["SpecifiedDistance_cm"] = dist / PhysicalUnits::cm;
      options_obj["SpecifiedDistance"] = PhysicalUnits::printToBestLengthUnits( dist, 6 );
    }
    
    options_obj["UseBq"] = options.use_bq;
    options_obj["IsSpecifyingDetector"] = !!options.drf_override;
    if( options.drf_override )
      options_obj["SpecifiedDetectorName"] = options.drf_override->name();
    
    options_obj["HardBackgroundSubtracted"] = !!options.hard_background_sub;
  }//add_activity_fit_options_to_json(...)
  
  
  void add_exe_info_to_json( nlohmann::json &data )
  {
    data["InterSpecCompileDate"] = string(__DATE__);
    data["InterSpecCompileDateIso"] = to_string( AppUtils::compile_date_as_int() );
    const auto now = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::now() );
    data["AnalysisTime"] = SpecUtils::to_iso_string( now );
    data["CurrentWorkingDirectory"] = SpecUtils::get_working_path();
  #if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
    try
    {
      const string exe_path = AppUtils::current_exe_path();
      data["InterSpecExecutablePath"] = exe_path;
    }catch( std::exception & )
    {
      assert( 0 );
    }
  #endif //if( not web or mobile )
  }//void add_exe_info_to_json( nlohmann::json &data )
  

  void add_peak_fit_options_to_json( nlohmann::json &data, const BatchPeak::BatchPeakFitOptions &options )
  {
    auto &options_obj = data["PeakFitOptions"];
    options_obj["FitAllPeaks"] = options.fit_all_peaks;
    options_obj["RefitEnergyCal"] = options.refit_energy_cal;
    options_obj["UseExemplarEnergyCal"] = options.use_exemplar_energy_cal;
    options_obj["WriteN42WithResults"] = options.write_n42_with_results;
    options_obj["ShowNonFitPeaks"] = options.show_nonfit_peaks;
    options_obj["OutputDir"] = options.output_dir;
    options_obj["CreateCsvOutput"] = options.create_csv_output;
    options_obj["CreateJsonOutput"] = options.create_json_output;
    options_obj["OverwriteOutputFiles"] = options.overwrite_output_files;
    options_obj["MultiSampleHandling"] = BatchSampleSelect::to_str( options.multi_sample_handling );
    
    if( !options.background_subtract_file.empty() || options.cached_background_subtract_spec )
    {
      string back_filename = options.background_subtract_file;
      if( back_filename.empty() && options.cached_background_subtract_spec )
        back_filename = options.cached_background_subtract_spec->filename();

      options_obj["BackgroundSubFile"] = back_filename;
      if( !options.background_subtract_samples.empty() )
      {
        options_obj["BackgroundSubSamples"] = vector<int>{ begin(options.background_subtract_samples),
          end(options.background_subtract_samples)
        };
      }
      options_obj["UsedExistingBackgroundPeak"] = options.use_existing_background_peaks;
      options_obj["UseExemplarEnergyCalForBackground"] = options.use_exemplar_energy_cal_for_background;
    }//if( !options.background_subtract_file.empty() || options.cached_background_subtract_spec )
    
    options_obj["ReportTemplateIncludeDir"] = options.template_include_dir;
    options_obj["ReportTemplates"] = options.report_templates;
    options_obj["ReportSummaryTemplates"] = options.summary_report_templates;
    options_obj["PeakStatThreshold"] = options.peak_stat_threshold;
    options_obj["PeakShapeThreshold"] = options.peak_hypothesis_threshold;
    options_obj["NotFitPeakMdaMethod"] = BatchPeak::to_str( options.not_fit_peak_mda );
    options_obj["MdaConfidenceLevel"] = options.mda_confidence_level;
    options_obj["MdaConfidenceLevelPercent"] = 100.0*options.mda_confidence_level;
    options_obj["MdaConfidenceLevelStr"] = DetectionLimitCalc::confidence_level_str( options.mda_confidence_level );
    options_obj["MdaNumSideChannels"] = static_cast<int>( options.mda_num_side_channels );
#if( ALLOW_SPECIFY_MDA_ROI_WIDTH )
    options_obj["MdaRoiNumFwhm"] = options.mda_roi_num_fwhm;
    options_obj["MdaSignalFractionInRoi"] = DetectionLimitCalc::gaussian_fraction_in_roi( options.mda_roi_num_fwhm );
#endif
  }//add_peak_fit_options_to_json(...)
  
  
  void add_peak_source_info_to_json( const PeakDef &peak, nlohmann::basic_json<> &peak_json )
  {
    peak_json["HasSourceAssigned"] = peak.hasSourceGammaAssigned();
    peak_json["SourceEnergy"] = peak.hasSourceGammaAssigned() ? peak.gammaParticleEnergy() : 0.0f;

    if( peak.parentNuclide() )
    {
      peak_json["SourceType"] = "Nuclide";
      peak_json["SourceName"] = peak.parentNuclide()->symbol;

      const SandiaDecay::Transition *trans = peak.nuclearTransition();
      if( trans && trans->parent )
        peak_json["SourceGammaParent"] = trans->parent->symbol;
      if( trans && trans->child )
        peak_json["SourceGammaChild"] = trans->child->symbol;
    }else if( peak.xrayElement() )
    {
      peak_json["SourceType"] = "X-Ray";
      peak_json["SourceName"] = peak.xrayElement()->name + " x-ray";
    }else if( peak.reaction() )
    {
      peak_json["SourceType"] = "Reaction";
      peak_json["SourceName"] = peak.reaction()->name();
    }else
    {
      peak_json["SourceType"] = "";
      peak_json["SourceName"] = "";
    }
  }//void add_peak_source_info_to_json(...)

  void add_peak_identity_to_json( const PeakDef &peak, nlohmann::basic_json<> &peak_json )
  {
    peak_json["PeakMean"] = peak.mean();
    peak_json["PeakMeanUncert"] = peak.meanUncert();
    peak_json["PeakFwhm"] = peak.gausPeak() ? peak.fwhm() : (peak.upperX() - peak.lowerX());
    peak_json["PeakAmplitude"] = peak.amplitude();
    peak_json["PeakAmplitudeUncert"] = peak.amplitudeUncert();
    peak_json["PeakMeanStr"] = SpecUtils::printCompact( peak.mean(), 5 );
    peak_json["PeakAmplitudeStr"] = SpecUtils::printCompact( peak.amplitude(), 4 );
    peak_json["PeakAmplitudeUncertStr"] = SpecUtils::printCompact( peak.amplitudeUncert(), 4 );

    add_peak_source_info_to_json( peak, peak_json );
  }//void add_peak_identity_to_json(...)



  void add_currie_check_to_json( nlohmann::basic_json<> &json,
                                 const DetectionLimitCalc::PeakCurrieCheck &check )
  {
    const DetectionLimitCalc::CurrieMdaResult &res = check.result;

    json["ConfidenceLevel"] = res.input.detection_probability;
    json["ConfidenceLevelPercent"] = 100.0*res.input.detection_probability;
    json["ConfidenceLevelStr"] = DetectionLimitCalc::confidence_level_str( res.input.detection_probability );
    json["NumSideChannels"] = static_cast<int>( res.input.num_lower_side_channels );
    json["SignalFractionInRoi"] = check.signal_fraction_in_roi;
    json["ShortDescription"] = check.short_description;
    json["ResultSummary"] = check.result_summary;
    json["RegionIsEmpty"] = check.region_is_empty;
    json["ResultType"] = DetectionLimitCalc::to_str( check.result_type );

    json["CurrieComputed"] = check.computed;
    if( !check.computed )
    {
      json["CurrieError"] = check.error_message;
      return;
    }

    json["DecisionThreshold_counts"] = res.decision_threshold;
    json["DetectionLimit_counts"] = res.detection_limit;
    json["SourceCounts"] = res.source_counts;
    json["UpperLimit_counts"] = res.upper_limit;
    json["LowerLimit_counts"] = res.lower_limit;
    json["PeakRegionCounts"] = res.peak_region_counts_sum;
    json["ContinuumCounts"] = res.estimated_peak_continuum_counts;
    json["ContinuumCountsUncert"] = res.estimated_peak_continuum_uncert;

    // ISO 11929-1:2019 clause 10.  `SourceCounts` is the primary result and is kept - it is what the
    //  detection decision is made on - while these are what ISO says to quote as the measured value.
    //  They are reported as a pair: the best estimate's own uncertainty is always the smaller of the
    //  two, so quoting it beside the primary result's would overstate the spread.
    json["BestEstimateCounts"] = res.best_estimate;
    json["BestEstimateCountsUncert"] = res.best_estimate_uncert;

    json["DecisionThresholdStr"] = SpecUtils::printCompact( res.decision_threshold, 4 );
    json["DetectionLimitStr"] = SpecUtils::printCompact( res.detection_limit, 4 );
    json["UpperLimitStr"] = SpecUtils::printCompact( res.upper_limit, 4 );
    json["SourceCountsStr"] = SpecUtils::printCompact( res.source_counts, 4 );
    json["BestEstimateCountsStr"] = SpecUtils::printCompact( res.best_estimate, 4 );
    json["BestEstimateCountsUncertStr"] = SpecUtils::printCompact( res.best_estimate_uncert, 4 );

    json["RoiLowerEnergy"] = res.input.roi_lower_energy;
    json["RoiUpperEnergy"] = res.input.roi_upper_energy;
    json["RoiWidth_keV"] = (res.input.roi_upper_energy - res.input.roi_lower_energy);
    json["RoiLowerChannel"] = static_cast<int>( res.first_peak_region_channel );
    json["RoiUpperChannel"] = static_cast<int>( res.last_peak_region_channel );
  }//void add_currie_check_to_json(...)


  void add_supplemental_peak_info_to_json( nlohmann::basic_json<> &data,
                    const std::vector<ShieldingSourceFitCalc::SupplementalPeakInfo> &supp_info,
                    const std::shared_ptr<const DetectorPeakResponse> &drf,
                    const bool useBq )
  {
    if( supp_info.empty() )
      return;

    const DetectorPeakResponse::EffGeometryType eff_type = drf ? drf->geometryType()
                                          : DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic;
    const bool use_curie = !useBq;
    const string act_postfix = DetectorPeakResponse::det_eff_geom_type_postfix( eff_type );

    nlohmann::basic_json<> &obj = data["SupplementalPeakInfo"];

    size_t num_currie_checks = 0, num_nominal_activities = 0, num_not_used = 0;

    for( const ShieldingSourceFitCalc::SupplementalPeakInfo &info : supp_info )
    {
      if( !info.peak )
        continue;

      num_currie_checks += info.currie.computed;

      // The peaks that were fit but not used by the model; the synthetic peaks stand in for peaks
      //  that were never observed, so they are not "peaks not used", they are not peaks at all.
      const bool is_not_used = (!info.used_for_fit && !info.synthetic);
      num_not_used += is_not_used;

      nlohmann::basic_json<> peak_json;
      add_peak_identity_to_json( *info.peak, peak_json );

      peak_json["UsedForFit"] = info.used_for_fit;
      peak_json["IsSyntheticPeak"] = info.synthetic;
      peak_json["HasCurrieCheck"] = info.currie.computed;
      add_currie_check_to_json( peak_json["CurrieCheck"], info.currie );

      peak_json["ModelEvaluated"] = info.model_evaluated;
      if( !info.model_evaluated )
      {
        peak_json["NotEvaluatedReason"] = info.not_evaluated_reason;
      }else
      {
        peak_json["ExpectedCounts"] = info.expected_counts;
        peak_json["ExpectedCountsStr"] = SpecUtils::printCompact( info.expected_counts, 4 );
        peak_json["NuclideExpectedCounts"] = info.nuclide_expected_counts;
        peak_json["OtherSourcesExpectedCounts"] = info.other_srcs_expected_counts;
        peak_json["SharedWithOtherSources"] = info.shared_with_other_sources;
        peak_json["CountsPerBq"] = info.counts_per_bq;
        peak_json["DetectorEff"] = info.det_efficiency;
        peak_json["ShieldingTransmission"] = info.shield_transmission;
        peak_json["AirTransmission"] = info.air_transmission;
        peak_json["IsVolumetricSource"] = info.is_volumetric_source;
      }//if( !info.model_evaluated ) / else

      peak_json["ShortDescription"] = info.short_description;
      peak_json["Description"] = info.description;
      peak_json["ResultSummary"] = info.result_summary;
      peak_json["Caveats"] = info.caveats;
      peak_json["HasCaveats"] = !info.caveats.empty();

      peak_json["HasNominalActivity"] = info.has_implied_activity;
      if( info.has_implied_activity )
      {
        num_nominal_activities += 1;

        /** Adds an activity, in a handful of units, under `<prefix>`, `<prefix>_bq`, etc.

         Follows the same conventions as `add_mda_to_json`: nothing is added for a non-positive
         activity, since `printToBestActivityUnits` renders those as an unreadable string of
         femtocuries.  Templates must guard with `existsIn(...)`.
         */
        const auto add_activity = [&peak_json,use_curie,&act_postfix]( const string &prefix,
                                                                       const double activity ){
          if( (activity <= 0.0) || IsNan(activity) || IsInf(activity) )
            return;

          peak_json[prefix] = PhysicalUnits::printToBestActivityUnits(activity,4,use_curie) + act_postfix;
          peak_json[prefix + "_bq"] = activity / PhysicalUnits::bq;
          peak_json[prefix + "_kBq"] = activity / PhysicalUnits::kBq;
          peak_json[prefix + "_MBq"] = activity / PhysicalUnits::MBq;
          peak_json[prefix + "_ci"] = activity / PhysicalUnits::ci;
          peak_json[prefix + "_mCi"] = activity / PhysicalUnits::mCi;
          peak_json[prefix + "_uCi"] = activity / PhysicalUnits::microCi;
        };//add_activity lambda

        add_activity( "ImpliedActivity", info.implied_activity );
        add_activity( "ImpliedActivityUncert", info.implied_activity_uncert );
        add_activity( "FitActivity", info.fit_activity );
        add_activity( "FitActivityUncert", info.fit_activity_uncert );

        peak_json["ActivityPostFix"] = act_postfix;
        peak_json["RatioToFitActivity"] = info.ratio_to_fit;
        peak_json["RatioToFitActivityUncert"] = info.ratio_to_fit_uncert;
        peak_json["RatioToFitActivityStr"] = SpecUtils::printCompact( info.ratio_to_fit, 4 );
        peak_json["NumSigmaOff"] = info.num_sigma_off;
        peak_json["NumSigmaOffStr"] = SpecUtils::printCompact( info.num_sigma_off, 3 );
      }//if( info.has_implied_activity )

      obj["AllPeaks"]["Peaks"].push_back( peak_json );

      if( is_not_used )
        obj["PeaksNotUsedForFit"]["Peaks"].push_back( peak_json );
    }//for( const SupplementalPeakInfo &info : supp_info )

    obj["HasCurrieChecks"] = (num_currie_checks > 0);
    obj["HasNominalActivities"] = (num_nominal_activities > 0);
    obj["AnyPeakNotUsedForFit"] = (num_not_used > 0);

    // Always present, so a template can loop without an `existsIn` guard on the array itself.
    if( !obj.contains("AllPeaks") )
      obj["AllPeaks"]["Peaks"] = nlohmann::json::array();
    if( !obj.contains("PeaksNotUsedForFit") )
      obj["PeaksNotUsedForFit"]["Peaks"] = nlohmann::json::array();
  }//void add_supplemental_peak_info_to_json(...)


  /** Names the continuum treatment used by a deconvolution limit, for the report.

   The treatment is per region of interest.  A batch limit uses one region per peak, but rather than
   report an arbitrary one if that ever stops being true, disagreement is reported as "mixed".
   `FixedByFullRange` is deprecated and unreachable from the UI, but a saved state can still carry
   it, so it is named rather than silently folded in with the others.
   */
  string decon_continuum_norm_str( const vector<DetectionLimitCalc::DeconRoiInfo> &roi_info )
  {
    if( roi_info.empty() )
      return "";

    const auto norm_name = []( const DetectionLimitCalc::DeconContinuumNorm norm ) -> string {
      switch( norm )
      {
        case DetectionLimitCalc::DeconContinuumNorm::Floating:         return "Floating";
        case DetectionLimitCalc::DeconContinuumNorm::FixedByEdges:     return "FixedByEdges";
        case DetectionLimitCalc::DeconContinuumNorm::FixedByFullRange: return "FixedByFullRange";
      }
      return "";
    };

    const string first = norm_name( roi_info.front().cont_norm_method );
    for( size_t i = 1; i < roi_info.size(); ++i )
    {
      if( norm_name( roi_info[i].cont_norm_method ) != first )
        return "mixed";
    }

    return first;
  }//string decon_continuum_norm_str(...)


  void add_mda_to_json( nlohmann::basic_json<> &mda_json, const BatchPeak::NotFitPeakMda &mda )
  {
    const DetectionLimitCalc::CurrieMdaResult &res = mda.currie.result;

    add_currie_check_to_json( mda_json, mda.currie );

    mda_json["Description"] = mda.description;
    mda_json["ActivitySummary"] = mda.activity_summary;
    mda_json["Caveats"] = mda.caveats;
    mda_json["HasCaveats"] = !mda.caveats.empty();
    mda_json["OverlapsFitPeak"] = mda.overlaps_fit_peak;
    mda_json["OverlapsOtherNotFitPeak"] = mda.overlaps_other_unfit_peak;
    mda_json["PeakWasFit"] = mda.peak_was_fit;

    if( mda.peak_was_fit && mda.fit_peak )
    {
      mda_json["FitPeakMean"] = mda.fit_peak->mean();
      mda_json["FitPeakFwhm"] = mda.fit_peak->fwhm();
      mda_json["FitPeakArea"] = mda.fit_peak->peakArea();
      mda_json["FitPeakAreaUncert"] = mda.fit_peak->peakAreaUncert();
      mda_json["FitPeakMeanStr"] = SpecUtils::printCompact( mda.fit_peak->mean(), 5 );
      mda_json["FitPeakAreaStr"] = SpecUtils::printCompact( mda.fit_peak->peakArea(), 4 );
      mda_json["FitPeakAreaUncertStr"] = SpecUtils::printCompact( mda.fit_peak->peakAreaUncert(), 4 );
    }//if( mda.peak_was_fit && mda.fit_peak )

    mda_json["HasActivity"] = mda.has_activity;
    if( !mda.has_activity )
    {
      if( !mda.no_activity_reason.empty() )
        mda_json["NoActivityReason"] = mda.no_activity_reason;
    }else
    {
      const bool use_curie = !mda.use_bq;
      const string &postfix = mda.activity_postfix;

      /** Adds an activity, in a handful of units, under `<prefix>`, `<prefix>_bq`, etc.

       Nothing is added for a non-positive activity: `PhysicalUnits::printToBestActivityUnits`
       renders those as an unreadable string of femtocuries, and an activity that works out
       negative (a deficit of counts, or a lower limit below zero) is not a number to put in a
       report.  Templates must guard with `existsIn(...)`.
       */
      const auto add_activity = [&mda_json,use_curie,&postfix]( const string &prefix, const double activity ){
        if( (activity <= 0.0) || IsNan(activity) || IsInf(activity) )
          return;

        mda_json[prefix] = PhysicalUnits::printToBestActivityUnits(activity,4,use_curie) + postfix;
        mda_json[prefix + "_bq"] = activity / PhysicalUnits::bq;
        mda_json[prefix + "_kBq"] = activity / PhysicalUnits::kBq;
        mda_json[prefix + "_MBq"] = activity / PhysicalUnits::MBq;
        mda_json[prefix + "_ci"] = activity / PhysicalUnits::ci;
        mda_json[prefix + "_mCi"] = activity / PhysicalUnits::mCi;
        mda_json[prefix + "_uCi"] = activity / PhysicalUnits::microCi;
      };//add_activity lambda

      mda_json["Nuclide"] = mda.nuclide ? mda.nuclide->symbol : string();
      mda_json["ShieldingTransmission"] = mda.shield_transmission;
      mda_json["AirTransmission"] = mda.air_transmission;
      mda_json["DetectorEff"] = mda.det_efficiency;
      mda_json["LiveTime_s"] = mda.live_time;
      mda_json["GammasPerBq"] = mda.gammas_per_bq;
      mda_json["ActivityPostFix"] = postfix;

      // Each of these is the activity form of the like-named counts field above; the minimum
      //  detectable activity ("MDA") is `DetectionLimitActivity`, i.e. Ld.
      add_activity( "DetectionLimitActivity", res.detection_limit / mda.gammas_per_bq );
      add_activity( "DecisionThresholdActivity", res.decision_threshold / mda.gammas_per_bq );

      // `add_activity` drops any of these that work out non-positive - e.g. the upper limit for
      //  a deficit of counts, or a lower limit that falls below zero.
      add_activity( "UpperLimitActivity", res.upper_limit / mda.gammas_per_bq );

      if( mda.currie.result_type == DetectionLimitCalc::PeakCurrieCheck::ResultType::Detected )
      {
        add_activity( "ObservedActivity", res.source_counts / mda.gammas_per_bq );
        add_activity( "ObservedActivityLower", res.lower_limit / mda.gammas_per_bq );
        add_activity( "ObservedActivityUpper", res.upper_limit / mda.gammas_per_bq );
      }//if( detected )

      // The activity form of the ISO best estimate and its uncertainty; see the counts fields above.
      //  Unlike `ObservedActivity` these are given whether or not a detection was declared, per ISO
      //  11929-1:2019 clause 10 NOTE 2 - subject, like every activity here, to `add_activity`
      //  dropping a value that is not positive, which for these means an empty region.
      add_activity( "BestEstimateActivity", res.best_estimate / mda.gammas_per_bq );
      add_activity( "BestEstimateActivityUncert", res.best_estimate_uncert / mda.gammas_per_bq );
    }//if( !mda.has_activity ) / else

    mda_json["DeconAttempted"] = mda.decon_attempted;
    if( mda.decon_attempted )
    {
      mda_json["DeconComputed"] = mda.decon_computed;
      mda_json["DeconQuantityIsCounts"] = mda.decon_quantity_is_counts;

      if( !mda.decon_error.empty() )
        mda_json["DeconError"] = mda.decon_error;

      if( mda.decon_computed && mda.decon_result )
      {
        const DetectionLimitCalc::DeconActivityOrDistanceLimitResult &decon = *mda.decon_result;

        mda_json["DeconFoundUpperLimit"] = decon.foundUpperCl;
        mda_json["MethodsDisagree"] = mda.methods_disagree;
        if( mda.decon_over_currie_ratio > 0.0 )
          mda_json["DeconOverCurrieRatio"] = mda.decon_over_currie_ratio;
        mda_json["DeconBestChi2"] = decon.overallBestChi2;
        mda_json["DeconLimitChi2"] = decon.upperLimitChi2;
        mda_json["DeconBestQuantity"] = decon.overallBestQuantity;

        // Which of the several distinct quantities this number is.  Without these a report cannot
        //  tell an upper bound on the loaded spectrum from a predicted sensitivity for a future
        //  one, and the two must never be worded the same way.
        mda_json["DeconContinuumNorm"] = decon_continuum_norm_str( decon.baseInput.roi_info );
        mda_json["DeconMeasurementModel"]
              = (decon.baseInput.measurement_model
                     == DetectionLimitCalc::DeconMeasurementModel::BackgroundReference)
                  ? "BackgroundReference" : "CurrentSpectrum";
        mda_json["DeconLimitType"]
              = (decon.limitType == DetectionLimitCalc::DeconLimitType::CentralInterval)
                  ? "CentralInterval" : "OneSidedUpperLimit";
        mda_json["DeconIsPredictedSensitivity"] = decon.is_predicted_sensitivity;
        // The REAL time, matching what the user entered and what every result string and the GUI
        //  quote.  `sampleExposure` is the live time the likelihood used, and quoting that here
        //  made batch and GUI disagree about the same calculation by the dead-time fraction.
        if( decon.is_predicted_sensitivity && (decon.sampleRealTime > 0.0) )
          mda_json["DeconSampleExposure_s"] = decon.sampleRealTime;

        if( decon.overallBestResults )
        {
          // Numerical-health diagnostics; iterations are normally small and restarts normally zero.
          mda_json["DeconOptimizerIterations"] = decon.overallBestResults->num_continuum_iterations;
          mda_json["DeconOptimizerRestarts"] = decon.overallBestResults->num_continuum_restarts;
          mda_json["DeconStatisticName"] = decon.overallBestResults->statistic_name;
        }

        // Things that qualify the limit without preventing it - most often that overlapping
        //  regions of interest were combined, or that the search range had to be extended.  These
        //  change how the number should be read, so they belong in the report.
        if( !decon.warnings.empty() )
        {
          nlohmann::json &warnings_json = mda_json["DeconWarnings"];
          for( const string &warning : decon.warnings )
            warnings_json.push_back( warning );
        }

        if( decon.foundUpperCl )
        {
          if( mda.decon_quantity_is_counts )
          {
            mda_json["DeconUpperLimit_counts"] = decon.upperLimit;
            mda_json["DeconUpperLimitStr"] = SpecUtils::printCompact( decon.upperLimit, 4 );
          }else
          {
            const bool use_curie = !mda.use_bq;
            mda_json["DeconUpperLimit"] = PhysicalUnits::printToBestActivityUnits(decon.upperLimit,4,use_curie)
                                          + mda.activity_postfix;
            mda_json["DeconUpperLimit_bq"] = decon.upperLimit / PhysicalUnits::bq;
            mda_json["DeconUpperLimit_ci"] = decon.upperLimit / PhysicalUnits::ci;
            mda_json["DeconUpperLimit_uCi"] = decon.upperLimit / PhysicalUnits::microCi;

            // These strings contain HTML, so are only suitable for HTML reports.
            mda_json["DeconLimitText"] = decon.limitText;
            mda_json["DeconBestChi2Text"] = decon.bestCh2Text;
          }//if( counts ) / else( activity )
        }//if( decon.foundUpperCl )
      }//if( mda.decon_computed && mda.decon_result )
    }//if( mda.decon_attempted )
  }//void add_mda_to_json(...)


  void add_peaks_to_json( nlohmann::json &json,
                          deque<shared_ptr<const PeakDef>> peaks,
                          const shared_ptr<const SpecUtils::Measurement> &spectrum,
                          const vector<BatchPeak::NotFitPeakMda> * const mdas )
  {
    // We will write down the peaks and continua in separate arrays, and give an index to link
    //  them.  By default the peaks and continua are sorted by energy, but we'll also through
    //  in some arrays of indexes so we can sort them in other orders when templating
    
    std::sort( begin(peaks), end(peaks), &PeakDef::lessThanByMeanShrdPtr );
    
    vector<shared_ptr<const PeakContinuum>> continua;
    for( const auto &p : peaks )
    {
      if( std::find(begin(continua), end(continua), p->continuum()) == end(continua) )
        continua.push_back( p->continuum() );
    }
    
    auto sorted_indices = [&peaks, spectrum]( Wt::SortOrder order, PeakModel::Columns column )
     -> vector<int> {
      deque<shared_ptr<const PeakDef>> peaks_copy = peaks;
      
      auto sortfcn = [column, order, spectrum]( const shared_ptr<const PeakDef> &lhs,
                                                const shared_ptr<const PeakDef> &rhs ) -> bool {
        return PeakModel::compare( lhs, rhs, column, order, spectrum );
      };
      stable_sort( begin(peaks_copy), end(peaks_copy), sortfcn );
      
      vector<int> indices;
      for( const shared_ptr<const PeakDef> &p : peaks )
      {
        const auto pos = std::find( begin(peaks_copy), end(peaks_copy), p );
        assert( pos != end(peaks_copy) );
        indices.push_back( static_cast<int>(pos - begin(peaks_copy)) );
      }
      return indices;
     };//sorted_indices lamda
    
    using So = Wt::SortOrder;
    using Col = PeakModel::Columns;
    json["PeakSortIndex_Energy_Ascend"] = sorted_indices( So::Ascending, Col::kMean );
    json["PeakSortIndex_Energy_Descend"] = sorted_indices( Wt::SortOrder::Descending, PeakModel::Columns::kMean );
    
    json["PeakSortIndex_Isotope_Ascend"] = sorted_indices( So::Ascending, Col::kIsotope );
    json["PeakSortIndex_Isotope_Descend"] = sorted_indices( Wt::SortOrder::Descending, PeakModel::Columns::kIsotope );
    
    json["PeakSortIndex_Mean_Ascend"] = sorted_indices( Wt::SortOrder::Ascending, PeakModel::Columns::kMean );
    json["PeakSortIndex_Mean_Descend"] = sorted_indices( Wt::SortOrder::Descending, PeakModel::Columns::kMean );
    
    json["PeakSortIndex_Amp_Ascend"] = sorted_indices( Wt::SortOrder::Ascending, PeakModel::Columns::kAmplitude );
    json["PeakSortIndex_Amp_Descend"] = sorted_indices( Wt::SortOrder::Descending, PeakModel::Columns::kAmplitude );
    
    json["PeakSortIndex_Fwhm_Ascend"] = sorted_indices( Wt::SortOrder::Ascending, PeakModel::Columns::kFwhm );
    json["PeakSortIndex_Fwhm_Descend"] = sorted_indices( Wt::SortOrder::Descending, PeakModel::Columns::kFwhm );
    
    json["PeakSortIndex_SrcEnergy_Ascend"] = sorted_indices( Wt::SortOrder::Ascending, PeakModel::Columns::kPhotoPeakEnergy );
    json["PeakSortIndex_SrcEnergy_Descend"] = sorted_indices( Wt::SortOrder::Descending, PeakModel::Columns::kPhotoPeakEnergy );
    
    json["PeakSortIndex_RoiCounts_Ascend"] = sorted_indices( Wt::SortOrder::Ascending, PeakModel::Columns::kRoiCounts );
    json["PeakSortIndex_RoiCounts_Descend"] = sorted_indices( Wt::SortOrder::Descending, PeakModel::Columns::kRoiCounts );
    
    json["PeakSortIndex_DistSrcEnergyToMean_Ascend"] = sorted_indices( Wt::SortOrder::Ascending, PeakModel::Columns::kDifference );
    json["PeakSortIndex_DistSrcEnergyToMean_Descend"] = sorted_indices( Wt::SortOrder::Descending, PeakModel::Columns::kDifference );
    
    json["PeakSortIndex_UseForActivity_Ascend"] = sorted_indices( Wt::SortOrder::Ascending, PeakModel::Columns::kUseForShieldingSourceFit );
    json["PeakSortIndex_UseForActivity_Descend"] = sorted_indices( Wt::SortOrder::Descending, PeakModel::Columns::kUseForShieldingSourceFit );
    
    json["PeakSortIndex_UseForEnergyCal_Ascend"] = sorted_indices( Wt::SortOrder::Ascending, PeakModel::Columns::kUseForCalibration );
    json["PeakSortIndex_UseForEnergyCal_Descend"] = sorted_indices( Wt::SortOrder::Descending, PeakModel::Columns::kUseForCalibration );
    
    
    for( int cont_index = 0; cont_index < static_cast<int>(continua.size()); ++cont_index )
    {
      const shared_ptr<const PeakContinuum> &cont = continua[cont_index];
      assert( cont );
      vector<int> peaks_with_cont;
      for( int peak_index = 0; peak_index < peaks.size(); ++peak_index )
      {
        if( peaks[peak_index]->continuum() == cont )
          peaks_with_cont.push_back( static_cast<int>(peak_index) );
      }
     
      json["Continua"].push_back( {} );
      auto &cont_json = json["Continua"].back();
      cont_json["PeakIndexes"] = peaks_with_cont;
      cont_json["ContinuumIndex"] = static_cast<int>( cont_index );
       
      const PeakContinuum::OffsetType type = cont->type();
      cont_json["ContinuumType"] = PeakContinuum::offset_type_str( type );
      const size_t npar = PeakContinuum::num_parameters( type );
      cont_json["NumberParameters"] = static_cast<int>( npar );
      cont_json["IsStepContinuum"] = PeakContinuum::is_step_continuum( type );
      cont_json["IsPolynomial"] = cont->isPolynomial();
      
      cont_json["IsEnergyRangeDefined"] = cont->energyRangeDefined();
      
      cont_json["LowerEnergy"] = cont->lowerEnergy();
      cont_json["UpperEnergy"] = cont->upperEnergy();
      
      assert( spectrum );
      if( spectrum && spectrum->energy_calibration() && spectrum->energy_calibration()->valid() )
      {
        const auto cal = spectrum->energy_calibration();
        const double lower_channel = cal->channel_for_energy( cont->lowerEnergy() );
        const double upper_channel = cal->channel_for_energy( cont->upperEnergy() );
        const int lower_channel_int = static_cast<int>( std::round(lower_channel) );
        const int upper_channel_int = static_cast<int>( std::round(upper_channel - 0.5) );
        const int num_channel = std::max( 0, 1 + upper_channel_int - lower_channel_int );
        
        cont_json["HasChannelRange"] = true;
        cont_json["LowerChannel"] = lower_channel;
        cont_json["UpperChannel"] = upper_channel;
        cont_json["LowerChannelInt"] = lower_channel_int;
        cont_json["UpperChannelInt"] = upper_channel_int;
        cont_json["NumberChannels"] = upper_channel - lower_channel;
        cont_json["NumberChannelsInt"] = num_channel;
        
        vector<int> channel_numbers( num_channel );
        vector<double> cont_area( num_channel, 0.0 ), channel_energies( num_channel + 1 );
        for( int i = lower_channel_int; i <= upper_channel_int; ++i )
        {
          const int index = i - lower_channel_int;
          const double channel_lower = cal->energy_for_channel( i );
          const double channel_upper = cal->energy_for_channel( i + 1 );
          
          channel_numbers[index] = i;
          channel_energies[index] = channel_lower;
          channel_energies[index + 1] = channel_upper;
          
          try
          {
            {
              vector<shared_ptr<const PeakDef>> roi_peaks;
              for( const int pi : peaks_with_cont )
                roi_peaks.push_back( peaks[pi] );
              cont_area[index] = cont->offset_integral( channel_lower, channel_upper, spectrum, roi_peaks );
            }
          }catch( std::exception & )
          {
            assert( 0 );
          }
        }//for( int i = lower_channel_int; i <= upper_channel_int; ++i )
        
        cont_json["ChannelNumbers"] = channel_numbers;
        cont_json["ChannelEnergies"] = channel_energies;
        cont_json["ChannelContinuumArea"] = cont_area;
      }else
      {
        cont_json["HasChannelRange"] = false;
        cont_json["LowerChannel"] = 0.0;
        cont_json["UpperChannel"] = 0.0;
        cont_json["NumberChannel"] = 0.0;
        cont_json["NumberChannelInt"] = 0;
        cont_json["LowerChannelInt"] = 0;
        cont_json["UpperChannelInt"] = 0;
      }
      
      try
      {
        {
          vector<shared_ptr<const PeakDef>> roi_peaks;
          for( const int pi : peaks_with_cont )
            roi_peaks.push_back( peaks[pi] );
          cont_json["ContinuumArea"] = cont->offset_integral( cont->lowerEnergy(), cont->upperEnergy(), spectrum, roi_peaks );
        }
      }catch( std::exception &e )
      {
        cont_json["ContinuumArea"] = 0.0;
        cont_json["Warnings"].push_back( "Error computing total integral area: " + string(e.what()) );
      }
      
      cont_json["ParameterReferenceEnergy"] = cont->referenceEnergy();
      cont_json["Parameters"] = cont->parameters();
      cont_json["ParameterUncertainties"] = cont->uncertainties();
      cont_json["ParameterIsForFitting"] = cont->fitForParameter();
    }//for( const shared_ptr<const PeakContinuum> &cont : continua )
    
    
    for( int peak_index = 0; peak_index < static_cast<int>(peaks.size()); ++peak_index )
    {
      const shared_ptr<const PeakDef> &p = peaks[peak_index];
      const auto cont_pos = std::find( begin(continua), end(continua), p->continuum() );
      assert( cont_pos != end(continua) );
      const auto continuum_index = cont_pos - begin(continua);
      assert( (continuum_index >= 0) && (continuum_index < continua.size()) );
      
      json["Peaks"].push_back( {} );
      auto &peak_json = json["Peaks"].back();
      
      peak_json["ContinuumIndex"] = static_cast<int>( continuum_index );
      
      peak_json["SkewType"] = PeakDef::to_string( p->skewType() );
      peak_json["NumSkewParameters"] = PeakDef::num_skew_parameters( p->skewType() );
      
      peak_json["PeakMean"] = p->mean();
      peak_json["PeakMeanUncert"] = p->meanUncert();
      
      peak_json["PeakAmplitude"] = p->amplitude();
      peak_json["PeakAmplitudeUncert"] = p->amplitudeUncert();
      
      peak_json["DataDefined"] = !p->gausPeak();
      peak_json["GaussianDefined"] = p->gausPeak();
      
      peak_json["Chi2Dof"] = p->chi2dof();
      peak_json["HasChi2Dof"] = p->chi2Defined();
      
      // Put peak info
      switch( p->type() )
      {
        case PeakDef::DefintionType::GaussianDefined:
        {
          peak_json["PeakSigma"] = p->sigma();
          peak_json["PeakSigmaUncert"] = p->sigmaUncert();
          peak_json["PeakFwhm"] = p->fwhm();
          peak_json["PeakFwhmUncert"] = 2.35482*p->sigmaUncert();
          
          break;
        }//case DefintionType::GaussianDefined:
          
        case PeakDef::DefintionType::DataDefined:
        {
          // We'll add in placeholders for sigma and FWHM, so templates that expect these values wont be messed up
          peak_json["PeakSigma"] = 0.0;
          peak_json["PeakSigmaUncert"] = 0.0;
          peak_json["PeakFwhm"] = 0.0;
          peak_json["PeakFwhmUncert"] = 0.0;
          break;
        }//case DefintionType::DataDefined:
      }//switch( p->type() )
      
      
      peak_json["LowerEnergy"] = p->lowerX();
      peak_json["UpperEnergy"] = p->upperX();
      peak_json["RoiWidth"] = p->roiWidth();

      
      if( spectrum && spectrum->energy_calibration() && spectrum->energy_calibration()->valid() )
      {
        peak_json["HasChannelRange"] = true;
        const auto cal = spectrum->energy_calibration();
        const double lower_channel = cal->channel_for_energy( p->lowerX() );
        const double upper_channel = cal->channel_for_energy( p->upperX() );
        const int lower_channel_int = static_cast<int>( std::round(lower_channel) );
        const int upper_channel_int = static_cast<int>( std::round(upper_channel - 0.5) );
        const int num_channel = std::max( 0, 1 + upper_channel_int - lower_channel_int );
        
        peak_json["LowerChannel"] = lower_channel;
        peak_json["UpperChannel"] = upper_channel;
        peak_json["LowerChannelInt"] = lower_channel_int;
        peak_json["UpperChannelInt"] = upper_channel_int;
        peak_json["NumberChannels"] = upper_channel - lower_channel;
        peak_json["NumberChannelsInt"] = num_channel;

        // Peak mean & widths expressed in channel numbers (channel-unit equivalents of the
        //  keV PeakMean / PeakSigma / PeakFwhm fields above).  Widths are converted by mapping
        //  the mean-centered energy interval of that width to a channel count, which stays
        //  correct for non-linear energy calibrations.
        try
        {
          const double mean = p->mean();
          peak_json["PeakMeanChannel"] = cal->channel_for_energy( mean );

          // Maps a keV width, centered on the peak mean, to a width in channels.
          const auto width_to_channels = [cal,mean]( const double width_kev ) -> double {
            return cal->channel_for_energy( mean + 0.5*width_kev )
                   - cal->channel_for_energy( mean - 0.5*width_kev );
          };

          peak_json["PeakMeanUncertChannel"] = width_to_channels( p->meanUncert() );

          if( p->gausPeak() )
          {
            peak_json["PeakSigmaChannel"]       = width_to_channels( p->sigma() );
            peak_json["PeakSigmaUncertChannel"] = width_to_channels( p->sigmaUncert() );
            peak_json["PeakFwhmChannel"]        = width_to_channels( p->fwhm() );
            peak_json["PeakFwhmUncertChannel"]  = width_to_channels( 2.35482*p->sigmaUncert() );
          }else
          {
            peak_json["PeakSigmaChannel"]       = 0.0;
            peak_json["PeakSigmaUncertChannel"] = 0.0;
            peak_json["PeakFwhmChannel"]        = 0.0;
            peak_json["PeakFwhmUncertChannel"]  = 0.0;
          }
        }catch( std::exception & )
        {
          // channel_for_energy can throw (lower-channel-edge cal, energy out of range);
          //  not expected in practice, so just emit sentinels.
          peak_json["PeakMeanChannel"]        = -1.0;
          peak_json["PeakMeanUncertChannel"]  = -1.0;
          peak_json["PeakSigmaChannel"]       = -1.0;
          peak_json["PeakSigmaUncertChannel"] = -1.0;
          peak_json["PeakFwhmChannel"]        = -1.0;
          peak_json["PeakFwhmUncertChannel"]  = -1.0;
        }

        // TODO: put in arrays of gaussian integrals, data counts, and continuum integral
        //double gauss_integral( const double x0, const double x1 ) const;
        //void gauss_integral( const float *energies, double *channels, const size_t nchannel ) const;
      }else
      {
        peak_json["HasChannelRange"] = false;
      }
      
    
      peak_json["AreaBetweenContinuumAndData"] = p->areaFromData( spectrum );
      peak_json["UseForEnergyCal"] = p->useForEnergyCalibration();
      peak_json["UseForActivityFit"] = p->useForShieldingSourceFit();
      peak_json["UseForIsotopicsFromPeaks"] = p->useForManualRelEff();
      peak_json["UseForDetEffFit"] = p->useForDrfIntrinsicEffFit();
      peak_json["UseForDetFwhmFit"] = p->useForDrfFwhmFit();
      peak_json["HasPeakUserLabel"] = !p->userLabel().empty();
      peak_json["PeakUserLabel"] = p->userLabel();
      add_peak_source_info_to_json( *p, peak_json );
    
      
      vector<string> coef_names;
      for( PeakDef::CoefficientType c = PeakDef::CoefficientType(0);
          c < PeakDef::CoefficientType::NumCoefficientTypes;
          c = PeakDef::CoefficientType( c + 1 ) )
      {
        coef_names.push_back( PeakDef::to_string(c) );
      }
      
      peak_json["CoefficientValues"] = vector<double>( p->coefficients(), p->coefficients() + PeakDef::CoefficientType::NumCoefficientTypes );
      peak_json["CoefficientUncerts"] = vector<double>( p->uncertainties(), p->uncertainties() + PeakDef::CoefficientType::NumCoefficientTypes );
      peak_json["CoefficientFit"] = vector<bool>( p->fitFors(), p->fitFors() + PeakDef::CoefficientType::NumCoefficientTypes );
      peak_json["CoefficientNames"] = coef_names;
      
      const Wt::WColor &color = p->lineColor();
      peak_json["PeakColor"] = color.cssText();

      // Detection limit ("MDA") info; only present for peaks that could not be fit, and only
      //  when limits were computed.  Never added for peaks that were fit, so that templates
      //  written before this existed see no change.
      if( mdas )
      {
        const BatchPeak::NotFitPeakMda *mda = nullptr;
        for( size_t i = 0; !mda && (i < mdas->size()); ++i )
        {
          if( (*mdas)[i].exemplar_peak == p )
            mda = &((*mdas)[i]);
        }

        peak_json["HasMda"] = !!mda;
        if( mda )
          add_mda_to_json( peak_json["Mda"], *mda );
      }//if( mdas )
    }//for( const shared_ptr<const PeakDef> &p : peaks )
  }//void add_peaks_to_json(...)


  /** Fills out `json` with the peaks that could not be fit, and their detection limits.

   Always sets `json["HasMdas"]`, so templates can use it without an `existsIn` guard.
   */
  void add_not_fit_peaks_to_json( nlohmann::json &json,
                                 const BatchPeak::BatchPeakFitResult &fit_results )
  {
    deque<shared_ptr<const PeakDef>> peaks( begin(fit_results.unfit_exemplar_peaks),
                                            end(fit_results.unfit_exemplar_peaks) );

    const vector<BatchPeak::NotFitPeakMda> &mdas = fit_results.not_fit_peak_mdas;

    add_peaks_to_json( json, peaks, fit_results.spectrum, mdas.empty() ? nullptr : &mdas );

    bool any_mda = false;
    for( const shared_ptr<const PeakDef> &p : peaks )
    {
      for( const BatchPeak::NotFitPeakMda &mda : mdas )
        any_mda = (any_mda || (mda.exemplar_peak == p));
    }

    json["HasMdas"] = any_mda;
  }//void add_not_fit_peaks_to_json(...)


  void add_peak_fit_results_to_json( nlohmann::basic_json<> &data,
                                    const BatchPeak::BatchPeakFitResult &fit_results )
  {
    data["Success"] = fit_results.success;
    data["Warnings"] = fit_results.warnings;
    data["HasWarnings"] = !fit_results.warnings.empty();

    data["HasSpectrum"] = !!fit_results.spectrum;
    if( fit_results.spectrum )
    {
      fit_results.spectrum->set_title( SpecUtils::filename(fit_results.file_path) );

      auto &spec_obj = data["foreground"];

      add_hist_to_json( spec_obj, true, 1.0,
                       fit_results.spectrum,
                       fit_results.measurement,
                       fit_results.sample_numbers,
                       SpecUtils::filename(fit_results.file_path),
                       &(fit_results.fit_peaks) );
    }//if( fit_results.spectrum )


    const bool successful_ene_refit = (fit_results.original_energy_cal && fit_results.refit_energy_cal);
    data["EnergyCalIsRefit"] = successful_ene_refit;

    if( successful_ene_refit )
    {
      auto &ene_cal_json = data["EnergyCalRefitResult"];

      if( fit_results.original_energy_cal )
        add_energy_cal_json( ene_cal_json["Original"], fit_results.original_energy_cal );

      if( fit_results.refit_energy_cal )
        add_energy_cal_json( ene_cal_json["Refit"], fit_results.refit_energy_cal );
    }//if( successful_ene_refit )


    data["FitAnyPeak"] = !fit_results.fit_peaks.empty();
    if( !fit_results.fit_peaks.empty() )
      add_peaks_to_json( data["FitPeaks"], fit_results.fit_peaks, fit_results.spectrum, nullptr );

    data["FitAllPeaks"] = fit_results.unfit_exemplar_peaks.empty();
    data["AnyNotFitPeakMda"] = !fit_results.not_fit_peak_mdas.empty();
    if( !fit_results.unfit_exemplar_peaks.empty() )
      add_not_fit_peaks_to_json( data["NotFitPeaks"], fit_results );

    data["ExemplarHasPeaks"] = !fit_results.exemplar_peaks.empty();
    if( !fit_results.exemplar_peaks.empty() )
    {
      // Every exemplar peak gets a detection limit, whether or not it was fit; for a peak that was
      //  fit this is a quality check.  Entries carry "PeakWasFit" to tell the two cases apart.
      const vector<BatchPeak::NotFitPeakMda> * const mdas
              = fit_results.exemplar_peak_mdas.empty() ? nullptr : &fit_results.exemplar_peak_mdas;

      add_peaks_to_json( data["ExemplarPeaks"], fit_results.exemplar_peaks,
                         fit_results.exemplar_spectrum, mdas );

      data["ExemplarPeaks"]["HasMdas"] = !!mdas;
    }//if( !fit_results.exemplar_peaks.empty() )

    data["AnyExemplarPeakMda"] = !fit_results.exemplar_peak_mdas.empty();

    
    // For peak searches, background subtraction are always a hard channel-by-channel subtraction,
    //  and `fit_results.spectrum` is after the subtraction
    //if( fit_results.background )
    //  add_hist_to_json( data["foreground"], true, display_scale_factor, fit_results.background, ... );
    //options.background_subtract_file;
    
     //std::shared_ptr<const SpecMeas> exemplar;
     //std::set<int> exemplar_sample_nums;
     
     //std::deque<std::shared_ptr<const PeakDef>> exemplar_peaks;
     //std::shared_ptr<const SpecUtils::Measurement> exemplar_spectrum;
     //std::vector<std::shared_ptr<const PeakDef>> unfit_exemplar_peaks;  //Exemplar peaks not found in the spectrum
  }//void add_peak_fit_results_to_json(...)


  void add_not_fit_peaks_to_act_shield_json( nlohmann::basic_json<> &data,
                                            const BatchPeak::BatchPeakFitResult &peak_fit_results )
  {
    data["AnyNotFitPeakMda"] = !peak_fit_results.not_fit_peak_mdas.empty();

    if( !peak_fit_results.unfit_exemplar_peaks.empty() )
      add_not_fit_peaks_to_json( data["NotFitPeaks"], peak_fit_results );
  }//void add_not_fit_peaks_to_act_shield_json(...)



  void write_json( const BatchPeak::BatchPeakFitOptions &options,
                  vector<string> &warnings,
                  const string &filename,
                  nlohmann::json json_copy )
  {
    string leaf_name = SpecUtils::filename(filename);
    if( leaf_name.empty() )
    {
      leaf_name = "summary.json";
    }else
    {
      const string file_ext = SpecUtils::file_extension(leaf_name);
      if( !file_ext.empty() )
        leaf_name = leaf_name.substr(0, leaf_name.size() - file_ext.size());
      leaf_name += "_results.json";
    }
    
    string out_json = SpecUtils::append_path(options.output_dir, leaf_name);
    
    if( SpecUtils::is_file(out_json) && !options.overwrite_output_files )
    {
      warnings.push_back( "Not writing '" + out_json + "', as it would overwrite a file."
                         " See the '--overwrite-output-files' option to force writing." );
    }else
    {
  #ifdef _WIN32
      const std::wstring wout_json = SpecUtils::convert_from_utf8_to_utf16(out_json);
      std::ofstream output_json( wout_json.c_str(), ios::binary | ios::out );
  #else
      std::ofstream output_json( out_json.c_str(), ios::binary | ios::out );
  #endif
      
      if( !output_json )
      {
        warnings.push_back( "Failed to open '" + out_json + "', for writing.");
      }else
      {
  #if( SpecUtils_ENABLE_D3_CHART )
        if( json_copy.count("D3_JS") )
          json_copy["D3_JS"] = "/* Removed for brevity - this string will have a value of the contents of the file InterSpec_resource/d3.v3.min.js during analysis in InterSpec_batch. */";
        if( json_copy.count("SpectrumChart_JS") )
          json_copy["SpectrumChart_JS"] = "/* Removed for brevity - this string will have a value of the contents of the file InterSpec_resource/SpectrumChartD3.js during analysis in InterSpec_batch.  */";
        if( json_copy.count("SpectrumChart_CSS") )
          json_copy["SpectrumChart_CSS"] = "/* Removed for brevity - this string will have a value of the contents of the file InterSpec_resource/SpectrumChartD3.css during analysis in InterSpec_batch. */";
        if( json_copy.count("ShieldingSourceFitPlot_JS") )
          json_copy["ShieldingSourceFitPlot_JS"] = "/* Removed for brevity - this string will have a value of the contents of the file InterSpec_resource/ShieldingSourceFitPlot.js during analysis in InterSpec_batch. */";
        if( json_copy.count("ShieldingSourceFitPlot_CSS") )
          json_copy["ShieldingSourceFitPlot_CSS"] = "/* Removed for brevity - this string will have a value of the contents of the file InterSpec_resource/ShieldingSourceFitPlot.css during analysis in InterSpec_batch. */";
  #endif // SpecUtils_ENABLE_D3_CHART
        
        output_json << std::setw(4) << json_copy << std::endl;
        cout << "Have written '" << out_json << "'" << endl;
      }
    }//if( SpecUtils::is_file( outcsv ) ) / else
  }//void write_json(...)
}//namespace BatchInfoLog
