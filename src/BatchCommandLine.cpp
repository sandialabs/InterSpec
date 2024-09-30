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

#include <regex>
#include <iostream>

#include <boost/program_options.hpp>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/AppUtils.h"
#include "InterSpec/BatchPeak.h"
#include "InterSpec/BatchActivity.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/BatchCommandLine.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;


namespace
{
  set<int> sequenceStrToSampleNums( string fulltxt )
  {
    //Lets replace numbers seperated by spaces, to be seperated by commas.
    //  We cant do a simple replace of spaces to commas, and using a regex would
    //  require using a lookahead or behind, and I dont think boost supports that
    //  always.  Note that below while loop is a little ineficient, but whatever
    std::smatch mtch;
    std::regex expr( ".*(\\d\\s+\\d).*" );
    while( std::regex_match( fulltxt, mtch, expr ) )
      fulltxt[mtch[1].first - fulltxt.begin() + 1] = ',';
    
    //Using a sregex_token_iterator doesnt work for single digit numbers, ex
    //  '1 2 3 45 983 193'  will go to '1,2 3,45,983,193'
    //  std::regex expr( "\\d\\s+\\d" );
    //  for( std::sregex_token_iterator iter(fulltxt.begin(),fulltxt.end(),expr,0);
    //      iter != std::sregex_token_iterator(); ++iter )
    //    fulltxt[iter->first - fulltxt.begin()+1] = ',';
      
    vector<string> sampleranges;
    SpecUtils::split( sampleranges, fulltxt, "," );

    set<int> sampleNumToLoad;
    
    for( string textStr : sampleranges )
    {
      SpecUtils::trim( textStr );
      if( textStr.empty() )
        continue;

      std::smatch matches;
      std::regex range_expression( "(\\d+)\\s*(\\-|to|through)\\s*(\\d+)",
                                  std::regex::ECMAScript | std::regex::icase );
      if( std::regex_match( textStr, matches, range_expression ) )
      {
        string firstStr = string( matches[1].first, matches[1].second );
        string lastStr = string( matches[3].first, matches[3].second );

        int first = std::stoi( firstStr );
        int last = std::stoi( lastStr );
        if( last < first )
          std::swap( last, first );

        for( int i = first; i <= last; ++i )
          sampleNumToLoad.insert( i );
      }else
      {
        const int sample = std::stoi( textStr );

        sampleNumToLoad.insert( sample );
      }//if( is sample range ) / else
    }//for( string textStr : sampleranges )
    
    return sampleNumToLoad;
  }//set<int> sequenceStrToSampleNums( string fulltxt )
  
}//namespace

namespace BatchCommandLine
{

int run_batch_command( int argc, char **argv )
{
  namespace po = boost::program_options;
  
  unsigned term_width = AppUtils::terminal_width();
  unsigned min_description_length = ((term_width < 80u) ? term_width/2u : 40u);
  
  po::options_description cl_desc("Allowed batch options", term_width, min_description_length);
  cl_desc.add_options()
  ("batch-peak-fit", "Batch-fit peaks.")
  ("batch-act-fit", "Batch shielding/source fit.")
  ;
  
  po::variables_map cl_vm;
  try
  {
    po::parsed_options parsed_opts
    = po::command_line_parser(argc,argv)
      .allow_unregistered()
      .options(cl_desc)
      .run();
    
    po::store( parsed_opts, cl_vm );
    po::notify( cl_vm );
  }catch( std::exception &e )
  {
    std::cerr << "Command line argument error: " << e.what() << std::endl << std::endl;
    std::cout << cl_desc << std::endl;
    return 1;
  }//try catch
  
  bool successful = true;
  
  const bool batch_peak_fit = cl_vm.count("batch-peak-fit");
  const bool batch_act_fit = cl_vm.count("batch-act-fit");
  
  if( batch_peak_fit || batch_act_fit )
  {
    try
    {
      if( batch_peak_fit && batch_act_fit )
        throw std::runtime_error( "You may not specify both 'batch-peak-fit' and 'batch-act-fit'." );
      
      bool output_stdout, refit_energy_cal, write_n42_with_results;
      bool use_exemplar_energy_cal, use_exemplar_energy_cal_for_background;
      bool show_nonfit_peaks, overwrite_output_files, create_csv_output, create_json_output;
      bool use_existing_background_peaks;
      double peak_stat_threshold = 0.0, peak_hypothesis_threshold = 0.0;
      vector<std::string> input_files, report_templates, summary_report_templates;
      string ini_file_path, exemplar_path, output_path, exemplar_samples, background_sub_file, background_samples, template_include_dir;
      
      po::options_description peak_cl_desc("Allowed batch peak-fit, and activity-fit options", term_width, min_description_length);
      peak_cl_desc.add_options()
      ("help,h",  "Produce help message")
      ("ini-file,i", po::value<std::string>(&ini_file_path)->default_value(""),
       "Path to INI file that can specify command line option defaults.\n"
       "(so you dont have to re-type things all the time)\n"
       "If not specified, will look for \"InterSpec_batch.ini\" in the current directory.")
      ("exemplar", po::value<std::string>(&exemplar_path),
       "File containing exemplar peaks, that will try to be fitted in the input spectra."
       " Can be a N42-2012 file save from InterSpec that contains peaks, or a peaks CSV file"
       " exported from the \"Peak Manager\" tab."
       )
      ("exemplar-sample-nums", po::value<std::string>(&exemplar_samples),
       "Only applicable if the exemplar file is an N42 file, and there are peaks for multiple"
       " sample numbers.  This string is interpreted similar to on the \"Spectrum Files\" tab,"
       " where a value such as \"1-3,8\" would use the peaks that were fit when samples 1,2,3, and 8"
       " where shown (e.g. summed for display), and peaks fit; usually you will just specify a"
       " single sample number though."
       )
      ("input-file", po::value<vector<std::string>>(&input_files)->multitoken(),
       "One or more spectrum files to fit peaks for.  If a directory, all files in it, recursively,"
       " will attempt to be used."
       )
      ("refit-energy-cal", po::value<bool>(&refit_energy_cal)->implicit_value(true)->default_value(false),
       "After initial peak fit, uses the those peaks (and their assigned nuclides) to adjust energy"
       " gain, then refits the peaks with the updated energy calibration."
       )
      ("use-exemplar-energy-cal", po::value<bool>(&use_exemplar_energy_cal)->implicit_value(true)->default_value(false),
       "Use the exemplar N42 energy calibration with the input foreground files.\n"
       " If refit-energy-cal is specified true, then the exemplar energy cal will be used for"
       " initial fit, and then refined and peaks refit.\n"
       "Only applicable if N42 file is used for exemplar."
       )
      ("use-exemplar-energy-cal-for-background",
       po::value<bool>(&use_exemplar_energy_cal_for_background)->implicit_value(true)->default_value(false),
       "Use the exemplar N42 energy calibration for the background file.\n"
       "Only applicable if N42 file is used for exemplar."
       )
      ("peak-stat-threshold",
       po::value<double>(&peak_stat_threshold)->default_value(2.0),
       "The improvement to the Chi2 of a peak fit required, over just fitting the continuum, to the ROI.\n"
       "A negative or zero value indicates no requirement (and default, since we are asserting peak"
       " is likely in the spectrum for batch analysis), and for general peak searching, reasonable"
       " values are between ~1 (a weak peak) and ~5 (a significant peak)."
       )
      ("peak-shape-threshold",
       po::value<double>(&peak_hypothesis_threshold)->default_value(1.0),
       "Requirement for how compatible the ROI must be to Gaussian peaks + continuum.\n"
       "It is the ratio of the null hypothesis chi2 (continuum only, no Gaussian),"
       "to the test hypothesis (continuum + Gaussian) chi2.\n"
       "A reasonable values for this seems to be in the 1 to 5 range.\n"
       "A zero or negative value will mean no requirement, and also no"
       "'peak-stat-threshold' requirement."
       )
      ("write-n42-with-results", po::value<bool>(&write_n42_with_results)->implicit_value(true)->default_value(false),
       "Adds the fit peaks to the input spectrum file , and then saves as a N42."
       " You must specify 'output_path' if you specify this option; also, will refuse to overwrite"
       " existing files."
       )
      ("overwrite-output-files", po::value<bool>(&overwrite_output_files)->implicit_value(true)->default_value(false),
       "Allows overwriting output N42, CSV, or report files.  By default will not overwrite files."
       )
      ("peak-csv-output", po::value<bool>(&create_csv_output)->implicit_value(true)->default_value(true),
       "Output peak fit CSV."
       )
      ("result-json-output", po::value<bool>(&create_json_output)->implicit_value(true)->default_value(false),
       "Writes the JSON used to create the report templates, out to file."
       )
      ("print", po::value<bool>(&output_stdout)->implicit_value(true)->default_value(false),
       "Print results to stdout."
       )
      ("out-dir", po::value<string>(&output_path)->default_value(""),
       "The directory to write peak CSV files (and optionally) N42 files to; if empty, CSV files will not be written.")
      ("back-sub-file", po::value<string>(&background_sub_file)->default_value(""),
       "File to use as the background.\n"
       "For activity/shielding fits, you can use the 'hard-background-subtract' option to specify if background"
       " subtraction should be on a peak-by-peak basis, or if a live-time-normalized, hard-background-subtraction (e.g., channel-by-channel).\n"
       "For peak search, a live-time-normalized, hard-background-subtraction (currently must be single record), will always be performed.")
      ("background-sample-nums", po::value<std::string>(&background_samples),
       "The sample numbers from the background file to use; if left empty will try to determine, and fail if not unique.\n\t"
       "Only applicable if the background subtraction file is specified."
       )
      ("use-existing-background-peaks", po::value<bool>(&use_existing_background_peaks)->default_value(false)->implicit_value(true),
       "Uses the fit peaks of the specified background, rather than re-fitting peaks.\n\t"
       "Peaks from InterSpec must exist in the spectrum file, and this option can not be"
       " used with hard background subtract."
       )
      ("include-nonfit-peaks", po::value<bool>(&show_nonfit_peaks)->default_value(false),
       "Include peaks that are not fit into the output CSV peak results; these peaks will have a zero amplitude set."
       )
      ("file-report-template", po::value<vector<string>>(&report_templates)->multitoken(),
       "One or more spectrum report template files - for each input file, a report will be produced with each template specified."
       " You can also specify \"default\" (=txt+html), \"txt\", \"html\", or \"none\"."
       )
      ("summary-report-template", po::value<vector<string>>(&summary_report_templates)->multitoken(),
       "One or more template files used to provide a summary of all input files."
       " You can also specify \"default\" (=csv+html), \"csv\", \"html\", or \"none\"."
       )
      ("report-template-include-dir", po::value<string>(&template_include_dir)->default_value("default"),
       "A directory containing report templates; specifying this allows templates to include other templates."
       " You can also specify \"default\", or \"none\"; specifying anything removes standard report template directory."
       )
      ;
      
      
      bool use_bq = false, hard_background_sub = false;
      try
      {
        use_bq = InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", nullptr );
      }catch( std::exception &)
      {
        assert( 0 );
      }
      
      string drf_file, drf_name, distance_override;
      po::options_description activity_cl_desc("Activity-fit options", term_width, min_description_length);
      activity_cl_desc.add_options()
      ("drf-file",  po::value<string>(&drf_file)->default_value(""),
                    "Path to a file containing DRF to use (overrides DRF contained in exemplar N42 file)")
      ("drf-name",  po::value<string>(&drf_name)->default_value(""),
                    "The name of a DRF to use; either within the file specified by 'drf-file'"
                    " (which may have multiple DRFs)"
                    ", or built-in, or previously used in InterSpec "
                    "(overrides DRF contained in exemplar N42 file).")
      ("distance", po::value<string>(&distance_override)->default_value(""),
                  "Optional distance to override default distance in exemplar.\n\t"
                  "Can not be specified with fixed geometry detector responses.\n\t"
                  "Ex: 100cm, 3ft, 10meters, 3.1inches, '1.1 cm + 1in', etc.")
      ("use-bq", po::value<bool>(&use_bq)->default_value(use_bq)->implicit_value(true),
       "Whether to use Curies or Becquerels for displaying activity.")
      ("hard-background-subtract", po::value<bool>(&hard_background_sub)->default_value(false)->implicit_value(true),
       "When true, the live-time normalized background is subtracted, channel-by-channel, "
       "from the foreground, before fitting for peaks.")
      ;
      
      po::variables_map cl_vm;
      try
      {
        po::parsed_options parsed_peak_opts
        = po::command_line_parser(argc,argv)
          .allow_unregistered()
          .options(peak_cl_desc)
          .run();
        
        po::store( parsed_peak_opts, cl_vm );
        
        if( batch_act_fit )
        {
          po::parsed_options parsed_act_opts
          = po::command_line_parser(argc,argv)
            .allow_unregistered()
            .options(activity_cl_desc)
            .run();
          po::store( parsed_act_opts, cl_vm );
        }//if( batch_act_fit )
        
        if( ini_file_path.empty() && SpecUtils::is_file("InterSpec_batch.ini") )
          ini_file_path = "InterSpec_batch.ini";
        
        if( !ini_file_path.empty() )
        {
          try 
          {
            ifstream input( ini_file_path.c_str() );
            if( !input )
              throw runtime_error( "Could not open config file '" + ini_file_path + "'" );
           
            // It *looks* like the INI file will NOT overwrite the values from the command line
            po::store( po::parse_config_file( input, peak_cl_desc, true), cl_vm );
            input.seekg( 0 );
            po::store( po::parse_config_file( input, activity_cl_desc, true), cl_vm );
            input.seekg( 0 );
            if( batch_act_fit )
            {
              input.seekg( 0 );
              po::store( po::parse_config_file(input, activity_cl_desc, true), cl_vm );
            }
            
            std::cout << "Using settings from '" << ini_file_path << "'" << endl;
          }catch( const std::exception & e )
          {
            std::cerr << "Error reading INI file: " << e.what() << " - skipping." << std::endl;
          }
        }//if( ini_file_path.empty() )
        
        po::notify( cl_vm );
      }catch( std::exception &e )
      {
        std::cerr << "Command line argument error: " << e.what() << std::endl << std::endl;
        std::cout << peak_cl_desc << std::endl;
        if( batch_act_fit )
          std::cout << activity_cl_desc << std::endl;
        return 1;
      }//try catch
      
      if( cl_vm.count("help") )
      {
        std::cout << "Available command-line options for batch peak or activity fitting are:\n";
        std::cout << peak_cl_desc << std::endl;
        if( batch_act_fit )
          std::cout << activity_cl_desc << std::endl;
        return 0;
      }//if( cl_vm.count("help") )
      
      set<int> exemplar_sample_nums;
      if( !exemplar_samples.empty() )
        exemplar_sample_nums = sequenceStrToSampleNums( exemplar_samples );
      
      // Expand directories to be files
      vector<string> expanded_input_files;
      for( const string &filename : input_files )
      {
        if( SpecUtils::is_directory(filename) )
        {
          const vector<string> subfiles = SpecUtils::recursive_ls(filename);
          expanded_input_files.insert( end(expanded_input_files), begin(subfiles), end(subfiles) );
        }else
        {
          expanded_input_files.push_back( filename );
        }
      }//for( string filename : input_files )
      
      
      set<int> background_sample_nums;
      if( !background_samples.empty() )
        background_sample_nums = sequenceStrToSampleNums( background_samples );
      
      if( !background_sample_nums.empty() && background_sub_file.empty() )
      {
        throw runtime_error( "You can not specify background sample numbers without specifying a background file" );
      }
      
      bool templates_contain_none = false, contain_not_none = false;
      for( const string &tmplt : report_templates )
      {
        templates_contain_none |= SpecUtils::iequals_ascii(tmplt, "none");
        contain_not_none |= !SpecUtils::iequals_ascii(tmplt, "none");
      }
      
      if( report_templates.empty() && (output_stdout || !output_path.empty()) )
        report_templates.push_back( "default" );
      
      if( summary_report_templates.empty() && (output_stdout || !output_path.empty()) )
        summary_report_templates.push_back( "default" );
      
      auto set_def = []( vector<string> &args, const vector<string> def_args ){
        const auto eq = []( const string &v ){ return SpecUtils::iequals_ascii(v, "default" ); };
        const size_t num_defaults = std::count_if( begin(args), end(args), eq);
        args.erase( std::remove_if( begin(args), end(args), eq), end(args) );
        if( num_defaults )
          args.insert( end(args), begin(def_args), end(def_args) );
        
        const auto want_none = [](const string &v){ return SpecUtils::iequals_ascii(v, "none"); };
        if( std::count_if( begin(args), end(args), want_none) )
          args.clear();
      };
      
      set_def( report_templates, {"txt", "html"} );
      set_def( summary_report_templates, {"csv", "html"} );
      
      if( contain_not_none && templates_contain_none )
        throw runtime_error( "You can not specify to use no report templates, and also specify to use other templates." );
      
      if( output_path.empty() && !summary_report_templates.empty() )
        throw runtime_error( "You must provide an output directory if specifying to use any templates for results." );
      
      if( hard_background_sub && background_sub_file.empty() )
        throw runtime_error( "You must specify a background spectrum when requesting a hard background subtract." );
      
      if( hard_background_sub && use_existing_background_peaks )
        throw runtime_error( "You can no use existing background peaks, and perform a hard background subtract." );
      
      if( use_existing_background_peaks && background_sub_file.empty() )
        throw runtime_error( "You must specify a background spectrum when requesting to use existing background peaks." );
      
      if( use_exemplar_energy_cal_for_background && background_sub_file.empty() )
        throw runtime_error( "You must specify a background spectrum when requesting to use exemplar energy calibration with background." );
      
      const string norm_background_path = SpecUtils::lexically_normalize_path(background_sub_file);
      const string norm_exemplar_path = SpecUtils::lexically_normalize_path(exemplar_path);
      if( (norm_background_path == norm_exemplar_path) && use_exemplar_energy_cal_for_background )
      {
        throw runtime_error( "You can not apply exemplar energy calibration to the background file,"
                            "when background file is the same as the exemplar (this could be a no-op,"
                            " but assuming an error in configuration)." );
      }
      
      boost::optional<double> distance_override_val;
      if( !distance_override.empty() )
      {
        double distance = -1.0;
        try
        {
          // If no units specified, default to cm
          if( distance_override.find_first_not_of( " \t0123456789.eE+-\n" ) == string::npos )
            distance_override += " cm";
          
          distance = PhysicalUnits::stringToDistance( distance_override );
        }catch( std::exception & )
        {
          throw runtime_error( "The distance string ('" + distance_override
                              + "') is not a valid distance" );
        }
        
        if( distance <= 0.0 )
          throw runtime_error( "Distance can not be negative or zero." );
        
        distance_override_val = distance;
      }//if( !distance_override.empty() )
      
      BatchActivity::BatchActivityFitOptions options;  //derived from BatchPeak::BatchPeakFitOptions
      options.to_stdout = output_stdout;
      options.refit_energy_cal = refit_energy_cal;
      options.use_exemplar_energy_cal = use_exemplar_energy_cal;
      options.write_n42_with_results = write_n42_with_results;
      options.show_nonfit_peaks = show_nonfit_peaks;
      options.overwrite_output_files = overwrite_output_files;
      options.create_csv_output = create_csv_output;
      options.create_json_output = create_json_output;
      options.output_dir = output_path;
      options.background_subtract_file = background_sub_file;
      options.background_subtract_samples = background_sample_nums;
      options.use_existing_background_peaks = use_existing_background_peaks;
      options.use_exemplar_energy_cal_for_background = use_exemplar_energy_cal_for_background;
      options.peak_stat_threshold = peak_stat_threshold;
      options.peak_hypothesis_threshold = peak_hypothesis_threshold;
      options.use_bq = use_bq;
      options.report_templates = report_templates;
      options.summary_report_templates = summary_report_templates;
      options.template_include_dir = template_include_dir;
      options.distance_override = distance_override_val;
      options.hard_background_sub = hard_background_sub;
      
      if( batch_peak_fit )
      {
        BatchPeak::fit_peaks_in_files( exemplar_path, exemplar_sample_nums, 
                                      expanded_input_files, options );
      }//if( batch_peak_fit )
      
      if( batch_act_fit )
      {
        options.drf_override = BatchActivity::init_drf_from_name( drf_file, drf_name );
        BatchActivity::fit_activities_in_files( exemplar_path, exemplar_sample_nums,
                                               expanded_input_files, options );
      }//if( batch_act_fit )
    }catch( std::exception &e )
    {
      successful = false;
      cerr << "Error performing batch analysis: " << e.what() << endl;
    }
  }else
  {
    successful = false;
    cerr << "Please specify either 'batch-peak-fit', or 'batch-act-fit'" << endl;
  }//if( batch_peak_fit || batch_act_fit )
  

  cout << "Done with batch processing." << endl;
  
  return successful ? EXIT_SUCCESS : EXIT_FAILURE;
}//int run_batch_command( int argc, const char **argv )
  
  
}//namespace BatchCommandLine
