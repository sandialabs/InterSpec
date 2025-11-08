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

#include "rapidxml/rapidxml.hpp"

#include "Minuit2/MnUserParameters.h"

#include "external_libs/SpecUtils/3rdparty/inja/inja.hpp"

#include "SpecUtils/DateTime.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "InterSpec/PeakDef.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/BatchPeak.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/BatchInfoLog.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/BatchActivity.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"


using namespace std;

namespace BatchActivity
{

const char *BatchActivityFitResult::to_str( const BatchActivityFitResult::ResultCode code )
{
  switch( code )
  {
    case ResultCode::CouldntInitializeStaticResources:     return "CouldntInitializeStaticResources";
    case ResultCode::NoExemplar:                           return "NoExemplar";
    case ResultCode::CouldntOpenExemplar:                  return "CouldntOpenExemplar";
    case ResultCode::ErrorPickingSpectrumFromExemplar:     return "ErrorPickingSpectrumFromExemplar";
    case ResultCode::CouldntOpenInputFile:                 return "CouldntOpenInputFile";
    case ResultCode::CouldntOpenBackgroundFile:            return "CouldntOpenBackgroundFile";
    case ResultCode::NoInputSrcShieldModel:                return "NoInputSrcShieldModel";
    case ResultCode::ForegroundSampleNumberUnderSpecified: return "ForegroundSampleNumberUnderSpecified";
    case ResultCode::BackgroundSampleNumberUnderSpecified: return "BackgroundSampleNumberUnderSpecified";
    case ResultCode::InvalidLiveTimeForHardBackSub:        return "InvalidLiveTimeForHardBackSub";
    case ResultCode::SpecifiedDistanceWithFixedGeomDet:    return "SpecifiedDistanceWithFixedGeomDet";
    case ResultCode::ErrorWithHardBackgroundSubtract:      return "ErrorWithHardBackgroundSubtract";
    case ResultCode::ErrorApplyingExemplarEneCalToFore:    return "ErrorApplyingExemplarEneCalToFore";
    case ResultCode::ForegroundPeakFitFailed:              return "ForegroundPeakFitFailed";
    case ResultCode::BackgroundPeakFitFailed:              return "BackgroundPeakFitFailed";
    case ResultCode::NoExistingBackgroundPeaks:            return "NoExistingBackgroundPeaks";
    case ResultCode::NoFitForegroundPeaks:                 return "NoFitForegroundPeaks";
    case ResultCode::NoDetEffFnct:                         return "NoDetEffFnct";
    case ResultCode::InvalidDistance:                      return "InvalidDistance";
    case ResultCode::InvalidGeometry:                      return "InvalidGeometry";
    case ResultCode::InvalidFitOptions:                    return "InvalidFitOptions";
    case ResultCode::ExemplarUsedBackSubButNoBackground:   return "ExemplarUsedBackSubButNoBackground";
    case ResultCode::NoShieldingsNode:                     return "NoShieldingsNode";
    case ResultCode::ErrorParsingShielding:                return "ErrorParsingShielding";
    case ResultCode::MissingNuclidesNode:                  return "MissingNuclidesNode";
    case ResultCode::InvalidNuclideNode:                   return "InvalidNuclideNode";
    case ResultCode::NoSourceNuclides:                     return "NoSourceNuclides";
    case ResultCode::FitNotSuccessful:                     return "FitNotSuccessful";
    case ResultCode::DidNotFitAllSources:                  return "DidNotFitAllSources";
    case ResultCode::FitThrewException:                    return "FitThrewException";
    case ResultCode::UnknownStatus:                        return "UnknownStatus";
    case ResultCode::Success:                              return "Success";
  }//switch( code )
  
  return "InvalidResultCode";
}//const char *BatchActivityFitResult::to_str( const ResultCode code )
  
shared_ptr<DetectorPeakResponse> init_drf_from_name( std::string drf_file, std::string drf_name )
{
  if( drf_file.empty() && drf_name.empty() )
    return nullptr;
    
  //  Check if we are over-riding the DRF to use
  if( !drf_file.empty() )
  {
    // If user specified a directory, assume its a GADRAS DRF
    if( SpecUtils::is_directory(drf_file) )
    {
      try
      {
        auto drf = make_shared<DetectorPeakResponse>();
        drf->fromGadrasDirectory( drf_file );
        return drf;
      }catch( std::exception &e )
      {
        throw runtime_error( "Directory '" + drf_file + "' wasnt a GADRAS DRF: " + string(e.what()) );
      }
    }//if( SpecUtils::is_directory(drf_file) )
    
    
    // Try a URI file
    if( SpecUtils::is_file(drf_file) && (SpecUtils::file_size(drf_file) < 64*1024) )
    {
      try
      {
        std::vector<char> data;
        SpecUtils::load_file_data( drf_file.c_str(), data );
        string url_query( begin(data), end(data) );
        SpecUtils::trim( url_query );
        auto drf = DetectorPeakResponse::parseFromAppUrl( url_query );
        if( drf )
          return drf;
      }catch( std::exception &e )
      {
        // Was not a file with URI as its contents
      }
    }else
    {
      try
      {
        string url_query = drf_file;
        SpecUtils::trim( url_query );
        auto drf = DetectorPeakResponse::parseFromAppUrl( url_query );
        if( drf )
          return drf;
      }catch( std::exception &e )
      {
        // Was not a URI
      }
    }//if( was a small file ) / else
    
    
    
    // Try a CSV/TSV file that may have multiple DRF in it
    fstream input( drf_file.c_str(), ios::binary | ios::in );
    if( !input )
      throw runtime_error( "Couldnt open DRF file '" + drf_file + "'." );
      
      
    vector<shared_ptr<DetectorPeakResponse>> drfs;
    try
    {
      vector<string> credits;
        
      // This should also cover files that could be parsed by
      //  `DetectorPeakResponse::parseSingleCsvLineRelEffDrf(string &)`
        
      DetectorPeakResponse::parseMultipleRelEffDrfCsv( input, credits, drfs );
    }catch(...)
    {
    }//try / catch for parseMultipleRelEffDrfCsv
      
      
    if( drfs.size() == 1 )
    {
      if( drf_name.empty() )
      {
        return drfs[0];
      }else
      {
        if( (drfs[0]->name() != drf_name) || (drfs[0]->description() != drf_name) )
          throw runtime_error( "The DRF in '" + drf_file + "' is named '"
                              + drfs[0]->name() + "', but you specified a name of '"
                                + drf_name + "'" );
        return drfs[0];
      }//if( drf_name.empty() ) / else
    }else if( drfs.size() > 1 )
    {
      string names;
      for( size_t i = 0; i < drfs.size(); ++i )
      {
        names += string(names.empty() ? "" : ", ") + "'" + drfs[i]->name() + "'";
        if( !drf_name.empty() && (drfs[i]->name() == drf_name) )
          return drfs[i];  // TODO: check if multiple DRFs with same name
      }
        
      if( drf_name.empty() )
        throw runtime_error( "No DRF name specified; possible DRFs in '" + drf_file
                              + "', are: " + names );
    }//if( drfs.size() == 1 ) / else
    
    
    // Try a Rel Eff CSV file
    try
    {
      input.seekg( 0, ios::beg );
      input.clear();
      
      auto det = DrfSelect::parseRelEffCsvFile( drf_file );
      if( det )
        return det;
    }catch( std::exception &e )
    {
      
    }//try catch
    
    
    // Try a ECC file
    try
    {
      input.seekg( 0, ios::beg );
      input.clear();
        
      auto answer = DetectorPeakResponse::parseEccFile( input );
      if( std::get<0>(answer) )
        return std::get<0>(answer);
    }catch( std::exception &e )
    {
      // Not a .ECC file
    }
  }else if( !drf_name.empty() )  //if( !drf_file.empty() )
  {
    try
    {
      string manufacturer = "";
      string model = drf_name;
      auto drf = DrfSelect::initARelEffDetector( SpecUtils::DetectorType::Unknown,
                                                manufacturer, model );
      if( drf )
        return drf;
    }catch( std::exception &e )
    {
      
    }//try /catch to load a standard Detector by name
    
    
    // TODO: Look for previous DRFs in the DB, that the user has used
      
  }//if( !drf_file.empty() ) / else if( !drf_name.empty() )


  // Try to open as a spectrum file, and see if there is a detector saved with it.
  //  (untested as of 20250527)
  {//begin try to open as spectrum file
    SpecMeas spec;
    const bool sucess = spec.load_file( drf_file, SpecUtils::ParserType::Auto, drf_file );
    if( sucess && spec.detector() )
      return spec.detector();
  }//end try to open as spectrum file

  
  throw runtime_error( "Could not load specified detector efficiency function."
                        + (drf_file.empty() ? string("") : " Filename='" + drf_file + "'.")
                        + (drf_name.empty() ? string("") : " Name='" + drf_name + "'.") );
  
  return nullptr;
}//shared_ptr<DetectorPeakResponse> init_drf_from_name( std::string drf_file, std::string drf-name )

  
void fit_activities_in_files( const std::string &exemplar_filename,
                          std::shared_ptr<const SpecMeas> cached_exemplar_n42,
                          const std::set<int> &exemplar_sample_nums,
                          const std::vector<std::string> &files,
                          std::vector<std::shared_ptr<SpecMeas>> optional_cached_files,
                          const BatchActivityFitOptions &options,
                          BatchActivityFitSummary * const summary_results )
{
  vector<string> warnings;
  
  if( files.empty() )
    throw runtime_error( "No input files specified." );

  if( !optional_cached_files.empty() && (optional_cached_files.size() != files.size()) )
    throw runtime_error( "If you specify cached files, you must specify the same number of files." );
  
  if( !options.output_dir.empty() && !SpecUtils::is_directory(options.output_dir) )
    throw runtime_error( "Output directory ('" + options.output_dir + "'), is not a directory." );
    
  if( options.write_n42_with_results && options.output_dir.empty() )
    throw runtime_error( "If you specify to write N42 files with the fit peaks, you must specify an output directory." );
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  if( !db )
  {
    const char *message = "Unable to load the nuclide decay database."
    " Run the executable from a working directory where the 'data' folder is "
    " located - the program expects to load data/sandia.decay.xml.";
    
    throw runtime_error( message );
  }//if( !db )
  
  
  const vector<pair<string,string>> spec_chart_js_and_css = BatchInfoLog::load_spectrum_chart_js_and_css();
  
  // Load report templates, and setup inja::environment
  inja::Environment env = BatchInfoLog::get_default_inja_env( options );
  
  bool set_setup_info_to_summary_json = false;
  nlohmann::json summary_json;
  
  BatchInfoLog::add_exe_info_to_json( summary_json );
  BatchInfoLog::add_activity_fit_options_to_json( summary_json, options );
  
  for( const pair<string,string> &key_val : spec_chart_js_and_css )
    summary_json[key_val.first] = key_val.second;
  
  summary_json["ExemplarFile"] = exemplar_filename;
  if( !exemplar_sample_nums.empty() )
    summary_json["ExemplarSampleNumbers"] = vector<int>{begin(exemplar_sample_nums), end(exemplar_sample_nums)};
  summary_json["InputFiles"] = files;
  
  if( summary_results )
  {
    summary_results->options = options;
    summary_results->exemplar_filename = exemplar_filename;
    summary_results->exemplar = cached_exemplar_n42;
    summary_results->exemplar_sample_nums = exemplar_sample_nums;
  }//if( summary_results )


  for( size_t file_index = 0; file_index < files.size(); ++file_index )
  {
    const string filename = files[file_index];
    string leaf_name = SpecUtils::filename(filename);
    const shared_ptr<SpecMeas> cached_file = optional_cached_files.empty() ? nullptr : optional_cached_files[file_index];

    const BatchActivityFitResult fit_results
                 = fit_activities_in_file( exemplar_filename, exemplar_sample_nums,
                                     cached_exemplar_n42, filename, cached_file, options );
    
    if( (fit_results.m_result_code == BatchActivityFitResult::ResultCode::CouldntOpenExemplar)
       || (fit_results.m_result_code == BatchActivityFitResult::ResultCode::CouldntOpenBackgroundFile) )
      throw runtime_error( fit_results.m_error_msg );
    
    if( !cached_exemplar_n42 )
      cached_exemplar_n42 = fit_results.m_exemplar_file;
    
    for( const string &warn : fit_results.m_warnings )
      warnings.push_back( "File '" + leaf_name + "': " + warn );
    
    if( !set_setup_info_to_summary_json && fit_results.m_fit_results )
    {
      std::shared_ptr<const DetectorPeakResponse> drf = options.drf_override;
      if( !drf && cached_exemplar_n42 )
        drf = cached_exemplar_n42->detector();
      
      const double distance = fit_results.m_fit_results->distance;
      const GammaInteractionCalc::GeometryType geometry = fit_results.m_fit_results->geometry;
      const ShieldingSourceFitCalc::ShieldingSourceFitOptions &fit_opts = fit_results.m_fit_results->options;
      BatchInfoLog::add_act_shield_fit_options_to_json( fit_opts, distance, geometry, drf, summary_json );
      set_setup_info_to_summary_json = true;
    }//if( !set_setup_info_to_summary_json )
    
    
    // The goal is to create a template that prints out the exact same info as as the current GUI
    //  log file, and then upgrade from this to hit HTML with the charts and peak fits and stuff.
    nlohmann::json data;
    
    data["ExemplarFile"] = exemplar_filename;
    if( !exemplar_sample_nums.empty() )
      data["ExemplarSampleNumbers"] = vector<int>{begin(exemplar_sample_nums), end(exemplar_sample_nums)};
    data["Filepath"] = filename;
    data["Filename"] = SpecUtils::filename( filename );
    data["ParentDir"] = SpecUtils::parent_path( filename );
    data["HasWarnings"] = !fit_results.m_warnings.empty();
    data["Warnings"] = fit_results.m_warnings;
    
    data["HasErrorMessage"] = !fit_results.m_error_msg.empty();
    if( !fit_results.m_error_msg.empty() )
      data["ErrorMessage"] = fit_results.m_error_msg;
    data["ResultCodeInt"] = static_cast<int>( fit_results.m_result_code );
    data["ResultCode"] = BatchActivityFitResult::to_str( fit_results.m_result_code );
    //m_background_file, m_background_sample_numbers
    
    
    if( fit_results.m_foreground )
    {
      auto &spec_obj = data["foreground"];
      
      BatchInfoLog::add_hist_to_json( spec_obj, false, fit_results.m_foreground,
                       fit_results.m_foreground_file,
                       fit_results.m_foreground_sample_numbers, fit_results.m_filename,
                       (fit_results.m_peak_fit_results ? &(fit_results.m_peak_fit_results->fit_peaks) : nullptr) );
    }//if( fit_results.m_foreground )
    
    if( fit_results.m_background && !fit_results.m_options.hard_background_sub )
    {
      string filename = fit_results.m_background_file ? fit_results.m_background_file->filename()
                        : string();
      
      // We'll only add background file info, if different than foreground file
      shared_ptr<const SpecMeas> back_file = (fit_results.m_background_file != fit_results.m_foreground_file)
                                                      ? fit_results.m_background_file : nullptr;
      
      auto &spec_obj = data["background"];
      BatchInfoLog::add_hist_to_json( spec_obj, true, fit_results.m_background,
                       back_file,
                       fit_results.m_background_sample_numbers, filename,
                       (fit_results.m_background_peak_fit_results ? &(fit_results.m_background_peak_fit_results->fit_peaks) : nullptr) );
      
      data["background"]["Normalization"] = fit_results.m_foreground->live_time() / fit_results.m_background->live_time();
    }//if( fit_results.m_background )
  
    
    BatchInfoLog::add_activity_fit_options_to_json( data, fit_results.m_options );
    
    for( const pair<string,string> &key_val : spec_chart_js_and_css )
      data[key_val.first] = key_val.second;
    
    const bool success = ((fit_results.m_result_code == BatchActivity::BatchActivityFitResult::ResultCode::Success)
                          && fit_results.m_fit_results);
    if( success )
    {
      cout << "Success analyzing '" << filename << "'!" << endl;
      assert( fit_results.m_fit_results );
      assert( fit_results.m_peak_fit_results && fit_results.m_peak_fit_results->measurement );
    }else
    {
      cout << "Failure analyzing '" << filename << "': " << fit_results.m_error_msg << endl;
      warnings.push_back( "Failed in analyzing '" + filename + "': " + fit_results.m_error_msg );
    }
    
    data["Success"] = success;
    
    assert( (fit_results.m_result_code != BatchActivity::BatchActivityFitResult::ResultCode::Success)
             || fit_results.m_fit_results );
    
    const bool useBq = UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    
    
    shared_ptr<const DetectorPeakResponse> drf = fit_results.m_options.drf_override;
    if( !drf && fit_results.m_exemplar_file )
      drf = fit_results.m_exemplar_file->detector();
    const bool use_bq = fit_results.m_options.use_bq;
      
    // Some info about the compiled application
    BatchInfoLog::add_exe_info_to_json( data );
    
    //Now add info about the analysis setup
    data["HasFitResults"] = !!fit_results.m_fit_results;
    if( fit_results.m_fit_results )
      BatchInfoLog::shield_src_fit_results_to_json( *fit_results.m_fit_results, drf, use_bq, data );
      
    if( summary_results )
    {
      summary_results->file_results.push_back( fit_results );
      summary_results->file_json.push_back( data.dump() );
      summary_results->file_peak_csvs.push_back( "" );              //We will fill this in below
      summary_results->file_reports.push_back( vector<string>() );  //We will fill this in below
      summary_results->exemplar = fit_results.m_exemplar_file;      //JIC we havent grabbed this yet
    }//if( summary_results )


    for( string tmplt : options.report_templates )
    {
      try
      {
        const string rpt = BatchInfoLog::render_template( tmplt, env,
                            BatchInfoLog::TemplateRenderType::ActShieldIndividual, options, data );
        
        if( summary_results )
          summary_results->file_reports.back().push_back( rpt );
        
        if( options.to_stdout && !SpecUtils::iequals_ascii(tmplt, "html" ) )
          cout << "\n\n" << rpt << endl << endl;
        
        if( !options.output_dir.empty() )
        {
          const string out_file
                    = BatchInfoLog::suggested_output_report_filename( filename, tmplt, 
                                  BatchInfoLog::TemplateRenderType::ActShieldIndividual, options );
          
          if( SpecUtils::is_file(out_file) && !options.overwrite_output_files )
          {
            warnings.push_back( "Not writing '" + out_file + "', as it would overwrite a file."
                               " See the '--overwrite-output-files' option to force writing." );
          }else
          {
#ifdef _WIN32
            const std::wstring woutcsv = SpecUtils::convert_from_utf8_to_utf16(out_file);
            std::ofstream output( woutcsv.c_str(), ios::binary | ios::out );
#else
            std::ofstream output( out_file.c_str(), ios::binary | ios::out);
#endif
            if( !output )
              warnings.push_back( "Failed to open report output '" + out_file + "', for writing.");
            else
              output.write( rpt.c_str(), rpt.size() );
          }//if( is file ) / else write file
        }//if( !options.output_dir.empty() )
      }catch( inja::InjaError &e )
      {
        const string msg = "Error templating results (" + e.type + ": line "
        + std::to_string(e.location.line) + ", column " + std::to_string(e.location.column)
        + " of '" + tmplt + "'): " + e.message + ". While processing '" + filename + "'.";
        
        cerr << msg << endl;
        warnings.push_back( msg );
      }catch( std::exception &e )
      {
        cerr << "Error templating results: " << e.what() << endl;
        warnings.push_back( "Error templating results: " + string(e.what()) );
      }
    }//for( const string &tmplt : options.report_templates )
    
    
    if( options.write_n42_with_results && fit_results.m_peak_fit_results
       && fit_results.m_peak_fit_results->measurement )
    {
      // TODO: need to have `fit_activities_in_file` add peaks and shielding model to a file,
      //       perhaps in `fit_results.m_peak_fit_results->measurement`, or maybe better yet,
      //       create whole new std::shared_ptr<SpecMeas> in BatchActivityFitResult
      const string message = "Written N42 file does not currently have Act/Shielding model"
      " written to it, only fit peaks, sorry - will.";
      cerr << message << endl;
      if( std::find(begin(warnings), end(warnings), message) == end(warnings) )
        warnings.push_back( message );
      
      const BatchPeak::BatchPeakFitResult &peak_fit_results = *fit_results.m_peak_fit_results;
      assert( peak_fit_results.measurement );
      
      string outn42 = SpecUtils::append_path(options.output_dir, SpecUtils::filename(filename) );
      if( !SpecUtils::iequals_ascii(SpecUtils::file_extension(filename), ".n42") )
        outn42 += ".n42";
      
      if( SpecUtils::is_file(outn42) && !options.overwrite_output_files )
      {
        warnings.push_back( "Not writing '" + outn42 + "', as it would overwrite a file."
                           " See the '--overwrite-output-files' option to force writing." );
      }else
      {
        if( !peak_fit_results.measurement->save2012N42File( outn42 ) )
          warnings.push_back( "Failed to write '" + outn42 + "'.");
        else
          cout << "Have written '" << outn42 << "' with peaks" << endl;
      }
    }//if( options.write_n42_with_results )
    
    if( !options.output_dir.empty() && options.create_json_output )
      BatchInfoLog::write_json( options, warnings, filename, data );

    if( fit_results.m_peak_fit_results )
    {
      const BatchPeak::BatchPeakFitResult &peak_fit_res = *fit_results.m_peak_fit_results;
      deque<shared_ptr<const PeakDef>> fit_peaks = peak_fit_res.fit_peaks;
      if( options.show_nonfit_peaks )
      {
        for( const auto p : peak_fit_res.unfit_exemplar_peaks )
        {
          auto np = make_shared<PeakDef>(*p);
          np->setAmplitude( 0.0 );
          np->setAmplitudeUncert( 0.0 );
          np->setSigmaUncert( 0.0 );
          auto cont = make_shared<PeakContinuum>( *np->continuum() );
          const size_t num_cont_pars = PeakContinuum::num_parameters(cont->type());
          for( size_t i = 0; i < num_cont_pars; ++i )
          {
            cont->setPolynomialCoef( i, 0.0 );
            cont->setPolynomialUncert( i, -1.0 );
          }
          np->setContinuum( cont );
          np->set_coefficient( -1.0, PeakDef::Chi2DOF );
          fit_peaks.push_back(np);
        }
        std::sort( begin(fit_peaks), end(fit_peaks), &PeakDef::lessThanByMeanShrdPtr );
      }//if( make_nonfit_peaks_zero )
      
      if( !options.output_dir.empty() && options.create_csv_output )
      {
        const string file_ext = SpecUtils::file_extension(leaf_name);
        if( !file_ext.empty() )
          leaf_name = leaf_name.substr(0, leaf_name.size() - file_ext.size());
        
        string outcsv = SpecUtils::append_path(options.output_dir, leaf_name) + "_peaks.CSV";
        
        if( SpecUtils::is_file(outcsv) && !options.overwrite_output_files )
        {
          warnings.push_back( "Not writing '" + outcsv + "', as it would overwrite a file."
                             " See the '--overwrite-output-files' option to force writing." );
        }else
        {
#ifdef _WIN32
          const std::wstring woutcsv = SpecUtils::convert_from_utf8_to_utf16(outcsv);
          std::ofstream output_csv( woutcsv.c_str(), ios::binary | ios::out );
#else
          std::ofstream output_csv( outcsv.c_str(), ios::binary | ios::out );
#endif
          
          if( !output_csv )
          {
            warnings.push_back( "Failed to open '" + outcsv + "', for writing.");
          }else
          {
            PeakModel::write_peak_csv( output_csv, leaf_name, PeakModel::PeakCsvType::Full,
                                      fit_peaks, peak_fit_res.spectrum );
            cout << "Have written '" << outcsv << "'" << endl;
          }
        }//if( SpecUtils::is_file( outcsv ) ) / else
      }//if( !options.output_dir.empty() && options.create_csv_output )
      
      if( options.to_stdout )
      {
        const string leaf_name = SpecUtils::filename(filename);
        cout << "peaks for '" << leaf_name << "':" << endl;
        PeakModel::write_peak_csv( cout, leaf_name, PeakModel::PeakCsvType::Full,
                                  fit_peaks, peak_fit_res.spectrum );
        cout << endl;
      }

      if( summary_results )
      {
        const string leaf_name = SpecUtils::filename(filename);
        stringstream ss;
        PeakModel::write_peak_csv( ss, leaf_name, PeakModel::PeakCsvType::Full,
                                      fit_peaks, peak_fit_res.spectrum );
        summary_results->file_peak_csvs.back() = ss.str();
      }
    }//if( fit_results.m_peak_fit_results )
    
    for( const pair<string,string> &key_val : spec_chart_js_and_css )
      data.erase(key_val.first);
    
    summary_json["Files"].push_back( data );
  }//for( const string filename : files )
  
  // Add any encountered errors to output summary JSON
  for( const string &warn : warnings )
    summary_json["Warnings"].push_back( warn );
  
  // Now write summary report(s)
  for( const string &summary_tmplt : options.summary_report_templates )
  {
    try
    {
      const string rpt = BatchInfoLog::render_template( summary_tmplt, env,
                       BatchInfoLog::TemplateRenderType::ActShieldSummary, options, summary_json );

     if( summary_results )
      summary_results->summary_reports.push_back( rpt );

      if( options.to_stdout && !SpecUtils::iequals_ascii(summary_tmplt, "html" ) )
        cout << "\n\n" << rpt << endl << endl;
      
      if( !options.output_dir.empty() )
      {
        const string out_file
                    = BatchInfoLog::suggested_output_report_filename( "", summary_tmplt,
                                    BatchInfoLog::TemplateRenderType::ActShieldSummary, options );
        
        if( SpecUtils::is_file(out_file) && !options.overwrite_output_files )
        {
          warnings.push_back( "Not writing '" + out_file + "', as it would overwrite a file."
                             " See the '--overwrite-output-files' option to force writing." );
        }else
        {
#ifdef _WIN32
          const std::wstring woutcsv = SpecUtils::convert_from_utf8_to_utf16(out_file);
          std::ofstream output( woutcsv.c_str(), ios::binary | ios::out );
#else
          std::ofstream output( out_file.c_str(), ios::binary | ios::out );
#endif
          if( !output )
            throw runtime_error( "Failed to open summary report output, '" + out_file + "'" );
          
          output.write( rpt.c_str(), rpt.size() );
        }//if( file exists and dont overwrite ) / else
      }
    }catch( inja::InjaError &e )
    {
      const string msg = "Error templating summary output (" + e.type + ": line "
      + std::to_string(e.location.line) + ", column " + std::to_string(e.location.column)
      + " of '" + summary_tmplt + "'): " + e.message + ".";
      
      cerr << msg << endl;
      warnings.push_back( msg );
    }catch( std::exception &e )
    {
      warnings.push_back( "Error making summary output: " + string(e.what()) );
    }
  }//if( !options.summary_report_template.empty() )
  
  
  if( !options.output_dir.empty() && options.create_json_output )
    BatchInfoLog::write_json( options, warnings, "", summary_json );
  
  if( !warnings.empty() )
    cerr << endl << endl;
  for( const auto warn : warnings )
    cerr << warn << endl;

  if( summary_results )
  {
    summary_results->summary_json = summary_json.dump();
    summary_results->warnings = warnings;
  }//if( summary_results )
}//fit_activities_in_files(...)
  

BatchActivityFitResult fit_activities_in_file( const std::string &exemplar_filename,
                          std::set<int> exemplar_sample_nums,
                          std::shared_ptr<const SpecMeas> cached_exemplar_n42,
                          const std::string &filename,
                          std::shared_ptr<SpecMeas> specfile,
                          const BatchActivityFitOptions &options )
{
  //  TODO: allow specifying, not just in the exemplar N42.  Also note InterSpec defines a URL encoding for model as well
  BatchActivityFitResult result;
  result.m_options = options;
  result.m_result_code = BatchActivityFitResult::ResultCode::UnknownStatus;
  result.m_filename = filename;
  result.m_exemplar_sample_numbers = exemplar_sample_nums;
  
  
  if( !cached_exemplar_n42 )
  {
    auto tmp = make_shared<SpecMeas>();
    if( !tmp->load_file( exemplar_filename, SpecUtils::ParserType::Auto ) )
    {
      result.m_error_msg = "Could not load exemplar '" + exemplar_filename + "'.";
      result.m_result_code = BatchActivityFitResult::ResultCode::CouldntOpenExemplar;
      return result;
    }//if( !meas->load_file( filename, ParserType::Auto ) )
    
    cached_exemplar_n42 = tmp;
  }//if( !cached_exemplar_n42 )
  
  result.m_exemplar_file = cached_exemplar_n42;
  
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
  {
    result.m_error_msg = "Could not initialize nuclide DecayDataBase.";
    result.m_result_code = BatchActivityFitResult::ResultCode::CouldntInitializeStaticResources;
    return result;
  }
  
  MaterialDB matdb;
  const string materialfile = SpecUtils::append_path( InterSpec::staticDataDirectory(), "MaterialDataBase.txt" );
  try
  {
    matdb.parseGadrasMaterialFile( materialfile, db, false );
  }catch( std::exception &e )
  {
    result.m_error_msg = "Could not initialize shielding material database.";
    result.m_result_code = BatchActivityFitResult::ResultCode::CouldntInitializeStaticResources;
    return result;
  }
  
  const Material * const iron = matdb.material("Fe (iron)");
  if( !iron )
  {
    result.m_error_msg = "Couldnt access materials from shielding database.";
    result.m_result_code = BatchActivityFitResult::ResultCode::CouldntInitializeStaticResources;
    return result;
  }
  
  if( !cached_exemplar_n42 )
  {
    result.m_error_msg = "No exemplar spectrum file provided";
    result.m_result_code = BatchActivityFitResult::ResultCode::NoExemplar;
    return result;
  }//if( !cached_exemplar_n42 )
  
  if( !specfile )
  {
    specfile = make_shared<SpecMeas>();
    if( !specfile->load_file( filename, SpecUtils::ParserType::Auto ) )
    {
      result.m_error_msg = "Could not load foreground '" + filename + "'.";
      result.m_result_code = BatchActivityFitResult::ResultCode::CouldntOpenInputFile;
      return result;
    }//if( !meas->load_file( filename, ParserType::Auto ) )
  }//if( !specfile )

  result.m_foreground_file = specfile;
  set<int> foreground_sample_numbers = result.m_foreground_sample_numbers;
  
  shared_ptr<SpecMeas> backfile;
  if( !options.background_subtract_file.empty() || options.cached_background_subtract_spec )
  {
    if( options.cached_background_subtract_spec )
    {
      backfile = options.cached_background_subtract_spec;
      result.m_background_file = backfile;
    }else if( options.background_subtract_file == filename )
    {
      backfile = specfile;
      result.m_background_file = specfile;
    }else if( options.background_subtract_file == exemplar_filename )
    {
      //backfile = cached_exemplar_n42;
      result.m_background_file = cached_exemplar_n42;
    }else
    {
      auto tmp = make_shared<SpecMeas>();
      if( !tmp->load_file(options.background_subtract_file, SpecUtils::ParserType::Auto) )
      {
        result.m_error_msg = "Could not load background '" + options.background_subtract_file + "'.";
        result.m_result_code = BatchActivityFitResult::ResultCode::CouldntOpenBackgroundFile;
        return result;
      }
      
      backfile = tmp;
      result.m_background_file = backfile;
    }//if( options.background_subtract_file == filename ) / else
  }//if( !options.background_subtract_file.empty() )
  
  
  // Find the sample number to use for either the foreground, or background measurement
  auto find_sample = []( shared_ptr<const SpecMeas> meas, const SpecUtils::SourceType wanted ) -> int {
    // This logic is probably repeated in a number of places throughout InterSpec now...
    assert( (wanted == SpecUtils::SourceType::Foreground)
           || (wanted == SpecUtils::SourceType::Background) );
    
    if( (wanted != SpecUtils::SourceType::Foreground)
           && (wanted != SpecUtils::SourceType::Background) )
    {
      throw std::logic_error( "Invalid src type" );
    }
    
    if( !meas )
      throw runtime_error( "No SpecMeas passed in." );
    
    const vector<string> &detectors = meas->detector_names();
    const set<int> &sample_nums = meas->sample_numbers();
    
    set<int> foreground_samples, background_samples;
    
    for( const int sample : sample_nums )
    {
      bool classified_sample_num = false;
      for( const string &det : detectors )
      {
        auto m = meas->measurement( sample, det );
        if( !m )
          continue;
        
        switch( m->source_type() )
        {
          case SpecUtils::SourceType::IntrinsicActivity:
          case SpecUtils::SourceType::Calibration:
            break;
            
          case SpecUtils::SourceType::Background:
            classified_sample_num = true;
            background_samples.insert(sample);
            break;
            
          case SpecUtils::SourceType::Foreground:
          case SpecUtils::SourceType::Unknown:
            classified_sample_num = true;
            foreground_samples.insert(sample);
            break;
        }//switch( m->source_type() )
        
        if( classified_sample_num )
          break;
      }//for( const string &det : detectors )
    }//for( const int sample : sample_nums )
    
    const set<int> &samples = (wanted == SpecUtils::SourceType::Foreground) ? foreground_samples 
                                                                            : background_samples;
    
    // If a file only has a single sample in it, and its marked background, but we are requesting
    //  foreground, use the single sample.  This happens, for example, when fitting peaks in a
    //  background file, where the sole spectrum has been marked as background.
    if( samples.empty()
       && (wanted == SpecUtils::SourceType::Foreground)
       && (background_samples.size() == 1) )
    {
      return *begin(background_samples);
    }
    
    if( samples.size() != 1 )
      throw runtime_error( "Sample number to use could not be uniquely identified." );
    
    return *begin(samples);
  };//int find_sample(...)
  
  
  shared_ptr<const SpecUtils::Measurement> foreground;
  try
  {
    const int sample_num = find_sample( specfile, SpecUtils::SourceType::Foreground );
    
    try
    {
      vector<shared_ptr<const SpecUtils::Measurement>> meass = specfile->sample_measurements(sample_num);
      if( meass.size() == 1 )
        foreground = meass[0];
      else
        foreground = specfile->sum_measurements( {sample_num}, specfile->detector_names(), nullptr );
      if( !foreground )
        throw runtime_error( "Missing measurement." );
    }catch( std::exception &e )
    {
      result.m_error_msg = "Could getting foreground measurement: " + string(e.what());
      result.m_result_code = BatchActivityFitResult::ResultCode::ForegroundSampleNumberUnderSpecified;
      return result;
    }//Try / catch get
    
    result.m_foreground = foreground;
    result.m_foreground_sample_numbers = set<int>{ sample_num };
    assert( result.m_foreground_file );
  }catch( std::exception &e )
  {
    result.m_error_msg = "Could not determine foreground: " + string(e.what());
    result.m_result_code = BatchActivityFitResult::ResultCode::ForegroundSampleNumberUnderSpecified;
    return result;
  }// try / catch (find foreground)

  
  set<int> background_sample_nums;
  shared_ptr<const SpecUtils::Measurement> background;
  try
  {
    if( (options.background_subtract_file == exemplar_filename) || (backfile == cached_exemplar_n42) )
    {
      if( options.background_subtract_samples.empty() )
      {
        const int sample_num = find_sample( cached_exemplar_n42, SpecUtils::SourceType::Background );
        background_sample_nums = set<int>{sample_num};
      }else
      {
        background_sample_nums = options.background_subtract_samples;
      }
      
      try
      {
        background = cached_exemplar_n42->sum_measurements( background_sample_nums, cached_exemplar_n42->detector_names(), nullptr );
      }catch( std::exception &e )
      {
        result.m_error_msg = "Error summing background measurements from exemplar: " + string(e.what());
        result.m_result_code = BatchActivityFitResult::ResultCode::BackgroundSampleNumberUnderSpecified;
        return result;
      }//Try / catch get
    }else if( backfile )
    {
      if( options.background_subtract_samples.empty() )
      {
        const SpecUtils::SourceType t = (backfile == specfile)
        ? SpecUtils::SourceType::Background
        : SpecUtils::SourceType::Foreground;
        
        const int sample_num = find_sample( backfile, t );
        
        background_sample_nums = set<int>{sample_num};
        
        try
        {
          vector<shared_ptr<const SpecUtils::Measurement>> meass = backfile->sample_measurements(sample_num);
          if( meass.size() == 1 )
            background = meass[0];
          else
            background = backfile->sum_measurements( background_sample_nums, backfile->detector_names(), nullptr );
          if( !background )
            throw runtime_error( "Missing measurement." );
        }catch( std::exception &e )
        {
          result.m_error_msg = "Could not get background measurement: " + string(e.what());
          result.m_result_code = BatchActivityFitResult::ResultCode::BackgroundSampleNumberUnderSpecified;
          return result;
        }//Try / catch get
      }else
      {
        background_sample_nums = options.background_subtract_samples;
        try
        {
          background = backfile->sum_measurements( background_sample_nums, backfile->detector_names(), nullptr );
          if( !background )
            throw runtime_error( "Missing measurement." );
        }catch( std::exception &e )
        {
          result.m_error_msg = "Could not sum background samples: " + string(e.what());
          result.m_result_code = BatchActivityFitResult::ResultCode::BackgroundSampleNumberUnderSpecified;
          return result;
        }//Try / catch get
      }//if( options.background_subtract_samples.empty() ) / else
    }//if( options.background_subtract_file == exemplar_filename ) / else
    
    result.m_background = background;
    result.m_background_sample_numbers = background_sample_nums;
    assert( (!!result.m_background_file) == (!!result.m_background) );
    assert( options.background_subtract_file.empty() == (!background) );
    assert( options.background_subtract_file.empty() == (!result.m_background_file) );
  }catch( std::exception &e )
  {
    result.m_error_msg = "Could not determine background: " + string(e.what());
    result.m_result_code = BatchActivityFitResult::ResultCode::BackgroundSampleNumberUnderSpecified;
    return result;
  }// try / catch (find foreground)
    
  
  // Apply exemplar energy cal - if wanted
  if( options.use_exemplar_energy_cal || options.use_exemplar_energy_cal_for_background )
  {
    
    set<int> exemplar_sample_nums;
    shared_ptr<const SpecUtils::Measurement> exemplar_spectrum;
    shared_ptr<const deque<std::shared_ptr<const PeakDef>>> exemplar_peaks;
    
    try
    {
      BatchPeak::get_exemplar_spectrum_and_peaks( exemplar_spectrum, exemplar_peaks,
                                                 exemplar_sample_nums, cached_exemplar_n42, true );
      if( !exemplar_spectrum )
        throw runtime_error( "Couldnt determine exemplar spectrum" );
      
      auto cal = exemplar_spectrum->energy_calibration();
      if( !cal || !cal->valid() )
        throw runtime_error( "Exemplar energy calibration invalid" );
    }catch( std::exception &e )
    {
      result.m_error_msg = "Error determining exemplar: " + string(e.what());
      result.m_result_code = BatchActivityFitResult::ResultCode::ErrorPickingSpectrumFromExemplar;
      return result;
    }
    
    shared_ptr<const SpecUtils::EnergyCalibration> exemplar_cal
                                                  = exemplar_spectrum->energy_calibration();
    // Make a copy of energy cal, jic, to keep things independednt (I dont think this is really necassary)
    if( options.use_exemplar_energy_cal || options.use_exemplar_energy_cal_for_background )
      exemplar_cal = make_shared<SpecUtils::EnergyCalibration>( *exemplar_cal );
    
    if( options.use_exemplar_energy_cal )
    {
      try
      {
        assert( foreground );
        auto fore = make_shared<SpecUtils::Measurement>( *foreground );
        
        BatchPeak::propagate_energy_cal( exemplar_cal, fore, specfile,
                                        result.m_foreground_sample_numbers );
        
        foreground = fore;
        result.m_foreground = fore;
      }catch( std::exception &e )
      {
        result.m_error_msg = "Error changing foreground energy calibration: " + string(e.what());
        result.m_result_code = BatchActivityFitResult::ResultCode::ErrorApplyingExemplarEneCalToFore;
        return result;
      }
    }//if( options.use_exemplar_energy_cal )
    
    if( background && options.use_exemplar_energy_cal_for_background )
    {
      try
      {
        assert( background );
        assert( options.background_subtract_file != exemplar_filename );
        auto back = make_shared<SpecUtils::Measurement>( *background );
        
        BatchPeak::propagate_energy_cal( exemplar_cal, back, backfile,
                                        result.m_background_sample_numbers );
        
        background = back;
        result.m_background = back;
      }catch( std::exception &e )
      {
        result.m_error_msg = "Error changing background energy calibration: " + string(e.what());
        result.m_result_code = BatchActivityFitResult::ResultCode::ErrorApplyingExemplarEneCalToFore;
        return result;
      }
    }
  }//if( options.use_exemplar_energy_cal || options.use_exemplar_energy_cal_for_background )
  
  
  if( options.hard_background_sub )
  {
    assert( background && result.m_foreground );
    try
    {
      if( !background )
        throw std::runtime_error( "no background spectrum" );
      if( !result.m_foreground )
        throw std::runtime_error( "no foreground spectrum" );
      
      // We _can_ subtract spectra with different number of channels, but its probably
      //  a mistake on the users part, so we'll throw an error
      if( background->num_gamma_channels() != result.m_foreground->num_gamma_channels() )
      {
        throw std::runtime_error( "Mis-match of channels between background ("
                                 + std::to_string(background->num_gamma_channels())
                                 + ") and foreground ("
                                 + std::to_string(result.m_foreground->num_gamma_channels())
                                 + ")" );
      }//if( background and foreground number of channels dont match )
  
      auto back_cal = background->energy_calibration();
      auto fore_cal = result.m_foreground->energy_calibration();
      
      if( !back_cal || !back_cal->valid() || !back_cal->channel_energies() )
        throw runtime_error( "Invalid background energy calibration" );
      
      if( !fore_cal || !fore_cal->valid() || !fore_cal->channel_energies() )
        throw runtime_error( "Invalid foreground energy calibration" );
      
  
      shared_ptr<const vector<float>> fore_counts = result.m_foreground->gamma_counts();
      shared_ptr<const vector<float>> back_counts = background->gamma_counts();
      
      // Make sure back_counts has the same energy calibration and fore_counts, so we can subtract
      //  on a bin-by-bin basis
      if( (back_cal != fore_cal) && (*back_cal) != (*fore_cal) )
      {
        auto new_backchan = make_shared<vector<float>>( fore_counts->size(), 0.0f );
        SpecUtils::rebin_by_lower_edge( *back_cal->channel_energies(), *back_counts,
                                       *fore_cal->channel_energies(), *new_backchan );
        back_counts = new_backchan;
      }
      
      // Create what will be the background subtracted foreground
      auto back_sub_counts = make_shared<vector<float>>( *fore_counts );
      
      //back_counts and fore_counts should always be the same size, but we'll be a bit verbose anyway
      assert( back_counts->size() == fore_counts->size() );
      const size_t nchann = std::min( back_counts->size(), fore_counts->size() );
      
      // Do the actual background subtraction
      const bool no_neg = true;
      const bool do_round = false;
      
      const float sf = result.m_foreground->live_time() / background->live_time();
      if( (sf <= 0.0) || IsNan(sf) || IsInf(sf) )
      {
        result.m_error_msg = "Could not determine live-times to perform hard background subtraction.";
        result.m_result_code = BatchActivityFitResult::ResultCode::InvalidLiveTimeForHardBackSub;
        return result;
      }
      
      for( size_t i = 0; i < nchann; ++i )
      {
        float &val = (*back_sub_counts)[i];
        val -= sf*(*back_counts)[i];
        
        if( no_neg )
          val = std::max( 0.0f, val );
        
        if( do_round )
          val = std::round( val );
      }//for( size_t i = 0; i < nchann; ++i )
      
      // Create a new Measurement object, based on the old foreground
      auto newspec = make_shared<SpecUtils::Measurement>( *result.m_foreground );
      newspec->set_gamma_counts( back_sub_counts, foreground->live_time(), foreground->real_time() );
      vector<string> remarks = foreground->remarks();
      remarks.push_back( "This spectrum has been background subtracted in InterSpec" );
      newspec->set_remarks( remarks );
      newspec->set_sample_number( 1 );
      
      // Create a new spectrum file object, and set new background subtracted Measurement as its only
      //  record
      auto newmeas = make_shared<SpecMeas>();
      
      // Copy just the SpecUtils::SpecFile stuff over to 'newmeas' so we dont copy things like
      //  displayed sample numbers, and uneeded peaks and stuff
      static_cast<SpecUtils::SpecFile &>(*newmeas) = static_cast<const SpecUtils::SpecFile &>(*result.m_foreground_file);
      newmeas->remove_measurements( newmeas->measurements() );
      newmeas->set_uuid( "" ); // Need to make sure UUID will get updated.
      newmeas->set_filename( "bkgsub_" + newmeas->filename() ); // Update filename
      newmeas->add_measurement( newspec, true ); // Actually add the measurement
      
      specfile = newmeas;
      foreground = newspec;
      result.m_foreground = newspec;
      foreground_sample_numbers.clear();
      foreground_sample_numbers.insert( newspec->sample_number() );
    }catch( std::exception &e )
    {
      result.m_error_msg = "Error performing hard background subtract: " + string(e.what());
      result.m_result_code = BatchActivityFitResult::ResultCode::ErrorWithHardBackgroundSubtract;
      return result;
    }//try / catch
  }//if( options.hard_background_sub )
  
  
  // Now need to figure out the foreground and possibly background sample numbers in the exemplar
  BatchPeak::BatchPeakFitOptions peak_fit_options = options;
  peak_fit_options.background_subtract_file = "";
  peak_fit_options.cached_background_subtract_spec = nullptr;
  peak_fit_options.use_exemplar_energy_cal = false;  //We've already applied this.
  peak_fit_options.use_exemplar_energy_cal_for_background = false;

  const BatchPeak::BatchPeakFitResult foreground_peak_fit_result
              = BatchPeak::fit_peaks_in_file( exemplar_filename, exemplar_sample_nums,
                                            cached_exemplar_n42, filename, specfile,
                                             foreground_sample_numbers, peak_fit_options );
  
  result.m_peak_fit_results = make_shared<BatchPeak::BatchPeakFitResult>( foreground_peak_fit_result );

  // Energy cal may have been modified - pick up these changes
  result.m_foreground = foreground_peak_fit_result.spectrum;
  result.m_foreground_file = foreground_peak_fit_result.measurement;
  

  if( !foreground_peak_fit_result.success )
  {
    result.m_error_msg = "Fitting of foreground peaks failed.";
    result.m_result_code = BatchActivityFitResult::ResultCode::ForegroundPeakFitFailed;
    
    return result;
  }//if( !foreground_peaks.success )
  
  
  if( (!options.background_subtract_file.empty() || background) && !options.hard_background_sub )
  {
    if( options.use_existing_background_peaks )
    {
      const set<int> &sample_nums = result.m_background_sample_numbers;
      const shared_ptr<const SpecMeas> &background = result.m_background_file;
      
      assert( background );
      assert( !options.hard_background_sub );
      
      shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks
                                              = background->peaks( sample_nums );
      
      if( !peaks || !peaks->size() )
      {
        result.m_error_msg = "It was requested to use existing background peaks, but there"
                            " were none for the background spectrum (please fit in InterSpec,"
                            " then export N42-2012 file).";
        result.m_result_code = BatchActivityFitResult::ResultCode::NoExistingBackgroundPeaks;
        
        return result;
      }//if( !peaks || !peaks->size() )
      
      auto background_peaks = make_shared<BatchPeak::BatchPeakFitResult>();
      background_peaks->file_path = exemplar_filename;
      background_peaks->options = options;
      background_peaks->exemplar = background;
      background_peaks->exemplar_sample_nums = result.m_background_sample_numbers;
      background_peaks->exemplar_peaks = *peaks;
      //background_peaks->unfit_exemplar_peaks;  //Exemplar peaks not found in the spectrum
      const vector<string> &det_names = background->detector_names();
      background_peaks->measurement = std::make_shared<SpecMeas>();
      static_cast<SpecUtils::SpecFile &>((*background_peaks->measurement)) = static_cast<const SpecUtils::SpecFile &>(*background);
      background_peaks->spectrum = background->sum_measurements( sample_nums, det_names, nullptr );
      assert( background_peaks->spectrum );
      if( !background_peaks->spectrum )
      {
        result.m_error_msg = "Failed to get background spectrum from file.";
        result.m_result_code = BatchActivityFitResult::ResultCode::NoExistingBackgroundPeaks;
        
        return result;
      }
      background_peaks->exemplar_spectrum = make_shared<SpecUtils::Measurement>( *background_peaks->spectrum );
      background_peaks->sample_numbers = sample_nums;
      background_peaks->fit_peaks = *peaks;
      
      background_peaks->background = nullptr;
      background_peaks->success = true;
      background_peaks->warnings = {};
      
      result.m_background_peak_fit_results = background_peaks;
    }else
    {
      const BatchPeak::BatchPeakFitResult background_peaks
      = BatchPeak::fit_peaks_in_file( exemplar_filename,
                                     exemplar_sample_nums,
                                     cached_exemplar_n42,
                                     options.background_subtract_file,
                                     backfile,
                                     result.m_background_sample_numbers,
                                     peak_fit_options );
      
      result.m_background_peak_fit_results = make_shared<BatchPeak::BatchPeakFitResult>( background_peaks );
      
      // Energy cal may have been modified - pick up these changes (as of 20250530, background is not )
      result.m_background = background_peaks.spectrum;
      result.m_background_file = background_peaks.measurement;
      
      if( !background_peaks.success )
      {
        result.m_error_msg = "Fitting of background peaks failed.";
        result.m_result_code = BatchActivityFitResult::ResultCode::BackgroundPeakFitFailed;
        
        return result;
      }//if( !foreground_peaks.success )
    }//if( options.use_existing_background_peaks )
  }//if( backfile )
  
  
  const rapidxml::xml_document<char> *model_xml = cached_exemplar_n42->shieldingSourceModel();
  const rapidxml::xml_node<char> *base_node = model_xml ? model_xml->first_node() : nullptr;
  if( !base_node )
  {
    result.m_error_msg = "Exemplar file did not have a Shielding/Source fit model in it";
    result.m_result_code = BatchActivityFitResult::ResultCode::NoInputSrcShieldModel;
    return result;
  }
  
  
  deque<shared_ptr<const PeakDef>> foreground_peaks = foreground_peak_fit_result.fit_peaks;
  if( foreground_peaks.empty() )
  {
    result.m_error_msg = "No foreground peaks fit.";
    result.m_result_code = BatchActivityFitResult::ResultCode::NoFitForegroundPeaks;
    return result;
  }
  
  shared_ptr<const DetectorPeakResponse> detector = options.drf_override;
  if( !detector )
    detector = cached_exemplar_n42->detector();
  
  if( !detector )
  {
    result.m_error_msg = "No detector efficiency function specified.";
    result.m_result_code = BatchActivityFitResult::ResultCode::NoDetEffFnct;
    return result;
  }
  
  const bool specified_dist = options.distance_override.has_value();
  if( specified_dist && detector && detector->isFixedGeometry() )
  {
    result.m_error_msg = "You can not specify a distance and use a fixed geometry detector.";
    result.m_result_code = BatchActivityFitResult::ResultCode::SpecifiedDistanceWithFixedGeomDet;
    return result;
  }
  
  
  vector<ShieldingSourceFitCalc::ShieldingInfo> shield_definitions;
  vector<ShieldingSourceFitCalc::SourceFitDef> src_definitions;
  ShieldingSourceFitCalc::ShieldingSourceFitOptions fit_options;
  GammaInteractionCalc::GeometryType geometry = GammaInteractionCalc::GeometryType::NumGeometryType;
  double distance = 1.0 * PhysicalUnits::meter;
  {
    GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceConfig config;
    try
    {
      config.deSerialize( base_node, &matdb );
    }catch( std::exception &e )
    {
      const string msg = e.what();
      result.m_error_msg = "Exemplar file '" + exemplar_filename + "' " + msg;
      if( msg.find("ShieldingSourceFitOptions") != string::npos )
        result.m_result_code = BatchActivityFitResult::ResultCode::InvalidFitOptions;
      else if( msg.find("Geometry") != string::npos )
        result.m_result_code = BatchActivityFitResult::ResultCode::InvalidGeometry;
      else if( msg.find("Distance") != string::npos )
        result.m_result_code = BatchActivityFitResult::ResultCode::InvalidDistance;
      else if( msg.find("Shieldings") != string::npos )
        result.m_result_code = (msg.find("missing") != string::npos)
                                ? BatchActivityFitResult::ResultCode::NoShieldingsNode
                                : BatchActivityFitResult::ResultCode::ErrorParsingShielding;
      else if( msg.find("Shielding") != string::npos )
        result.m_result_code = BatchActivityFitResult::ResultCode::ErrorParsingShielding;
      else if( msg.find("Nuclide") != string::npos || msg.find("sources XML") != string::npos )
        result.m_result_code = (msg.find("missing") != string::npos)
                                ? BatchActivityFitResult::ResultCode::MissingNuclidesNode
                                : BatchActivityFitResult::ResultCode::InvalidNuclideNode;
      else
        result.m_result_code = BatchActivityFitResult::ResultCode::FitThrewException;
      return result;
    }
    
    geometry = config.geometry;
    if( geometry == GammaInteractionCalc::GeometryType::NumGeometryType )
    {
      result.m_error_msg = "Exemplar file '" + exemplar_filename + "' has invalid geometry.";
      result.m_result_code = BatchActivityFitResult::ResultCode::InvalidGeometry;
      return result;
    }
    
    distance = config.distance;
    fit_options = config.options;
    shield_definitions = config.shieldings;
    src_definitions = config.sourceDefinitions();
  }
  
  if( specified_dist )
    distance = options.distance_override.value();

  shared_ptr<const deque<std::shared_ptr<const PeakDef>>> background_peaks;
  if( background 
     && result.m_background_peak_fit_results 
     && !options.hard_background_sub
     && result.m_background_peak_fit_results->success )
  {
    fit_options.background_peak_subtract = true; //Original exemplar may not have used background subtraction
    background_peaks = make_shared<deque<shared_ptr<const PeakDef>>>(result.m_background_peak_fit_results->fit_peaks);
  }
  
  
  // Check that if options said to use background-subtraction, that we actually can.
  //  I'm a little mixed here - perhaps if no background is specified, we should just ignore this,
  //  or just create a warning???
  //  For now - lets just warn if background subtraction isnt effective
  if( result.m_background_peak_fit_results 
     && !options.hard_background_sub
     && (!background_peaks || background_peaks->empty()) )
  {
    result.m_warnings.push_back( "No background peaks fit, although background subtraction requested" );
  }//if( result.m_background_peak_fit_results )
  
  
  if( fit_options.background_peak_subtract && (!background_peaks || background_peaks->empty()) )
  {
    fit_options.background_peak_subtract = false;
    string msg = "Exemplar file '" + exemplar_filename + "' Options indicated"
    " background-subtract, but either no background was given, or no background peaks fit.";
    result.m_warnings.push_back( msg );
    //result.m_error_msg = msg;
    //result.m_result_code = BatchActivityFitResult::ResultCode::ExemplarUsedBackSubButNoBackground;
    //return result;
  }
  if( src_definitions.empty() )
  {
    result.m_error_msg = "Exemplar file '" + exemplar_filename + "' no sources defined.";
    result.m_result_code = BatchActivityFitResult::ResultCode::NoSourceNuclides;
    return result;
  }//if( src_definitions.empty() )
  
  // We may not have fit peaks for all `src_definitions`.  Lets remove these sources,
  //  and add warnings about it
  {// Begin remove sources without peaks
    vector<ShieldingSourceFitCalc::SourceFitDef> srcs_with_peaks;
    for( const ShieldingSourceFitCalc::SourceFitDef &src : src_definitions )
    {
      const SandiaDecay::Nuclide * const nuclide = src.nuclide;
      assert( nuclide );
      if( !nuclide )
        continue;
      
      bool have_nuc_in_peak = false;
      for( const auto &p : foreground_peaks )
        have_nuc_in_peak |= (p->parentNuclide() == nuclide);
      
      if( have_nuc_in_peak )
        srcs_with_peaks.push_back( src );
      else
        result.m_warnings.push_back( "No peak assigned to nuclide " + nuclide->symbol
                                    + ", not using this nuclide." );
    }//for( const ShieldingSourceFitCalc::SourceFitDef &src : src_definitions )
    
    srcs_with_peaks.swap( src_definitions );
  }// End remove sources without peaks
  
  try
  {
    // We have all the parts, lets do the computation:
    GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceInput chi_input;
    chi_input.config.distance = distance;
    chi_input.config.geometry = geometry;
    chi_input.config.shieldings = shield_definitions;
  chi_input.config.setSourceDefinitions( src_definitions );
    chi_input.config.options = fit_options;
    chi_input.detector = detector;
    chi_input.foreground = foreground;
    chi_input.background = background;
    chi_input.foreground_peaks = foreground_peaks;
    chi_input.background_peaks = background_peaks;
    pair<shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> fcn_pars =
    GammaInteractionCalc::ShieldingSourceChi2Fcn::create( chi_input );
    
    auto inputPrams = make_shared<ROOT::Minuit2::MnUserParameters>();
    *inputPrams = fcn_pars.second;
    
    auto progress = make_shared<ShieldingSourceFitCalc::ModelFitProgress>();
    auto fit_results = make_shared<ShieldingSourceFitCalc::ModelFitResults>();
    result.m_fit_results = fit_results;
    
    auto progress_fcn = [progress](){
      // We dont really care too much probably right now
    };
    
    bool finished_fit_called = false;
    auto finished_fcn = [fit_results, &finished_fit_called](){
      finished_fit_called = true;
    };
    
    ShieldingSourceFitCalc::fit_model( "", fcn_pars.first, inputPrams, progress, progress_fcn, fit_results, finished_fcn );
    
    assert( finished_fit_called );
    
    result.m_result_code = BatchActivityFitResult::ResultCode::Success;
    
    if( fit_results->successful != ShieldingSourceFitCalc::ModelFitResults::FitStatus::Final )
    {
      result.m_error_msg = "Fit not successful.";
      result.m_result_code = BatchActivityFitResult::ResultCode::FitNotSuccessful;
      return result;
    }
    
    
    for( const ShieldingSourceFitCalc::SourceFitDef insrc : src_definitions )
    {
      assert( insrc.nuclide );
      bool found_in_output = false;
      for( const ShieldingSourceFitCalc::SourceFitDef &fitsrc : fit_results->fit_src_info )
      {
        assert( fitsrc.nuclide );
        found_in_output = (fitsrc.nuclide == insrc.nuclide);
        if( found_in_output )
          break;
      }//for( const ShieldingSourceFitCalc::SourceFitDef &fitsrc : fit_results->fit_src_info )
      
      if( !found_in_output )
      {
        result.m_warnings.push_back( "Fit results does not include " + insrc.nuclide->symbol );
        result.m_result_code = BatchActivityFitResult::ResultCode::DidNotFitAllSources;
      }//if( !found_in_output )
    }//for( const ShieldingSourceFitCalc::SourceFitDef insrc : src_definitions )
    
    // TODO: create an output file that has peaks, drf, and the fit model.  This will mean creating the XML that represents the fit model then turning it into a string
  }catch( std::exception &e )
  {
    result.m_error_msg = e.what();
    result.m_result_code = BatchActivityFitResult::ResultCode::FitThrewException;
  }//try / catch to do the fit
  
  cout << "Finished fitting for '" << filename << "'" << endl;
  
  return result;
}//fit_activities_in_file(...)
  
}//namespace BatchPeak

