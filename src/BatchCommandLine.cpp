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
      
      bool output_stdout, refit_energy_cal, use_exemplar_energy_cal, write_n42_with_peaks, show_nonfit_peaks;
      vector<std::string> input_files;
      string exemplar_path, output_path, exemplar_samples, background_sub_file, background_samples;
      
      po::options_description peak_cl_desc("Allowed batch peak-fit, and activity-fit options", term_width, min_description_length);
      peak_cl_desc.add_options()
      ("help,h",  "Produce help message")
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
      ("refit-energy-cal", po::value<bool>(&refit_energy_cal)->default_value(false),
       "After initial peak fit, uses the those peaks (and their assigned nuclides) to adjust energy"
       " gain, then refits the peaks with the updated energy calibration."
       )
      ("use-exemplar-energy-cal", po::value<bool>(&use_exemplar_energy_cal)->default_value(false),
       "Use the exemplar N42 energy calibration with the input files."
       " If refit-energy-cal is specified true, then the exemplar energy cal will be used for"
       " initial fit, and then refined and peaks refit.\n"
       "Only applicable if N42 file is used for exemplar."
       )
      ("write-n42-with-peaks", po::value<bool>(&write_n42_with_peaks)->default_value(false),
       "Adds the fit peaks to the input spectrum file , and then saves as a N42."
       " You must specify 'output_path' if you specify this option; also, will refuse to overwrite"
       " existing files."
       )
      ("print", po::value<bool>(&output_stdout)->default_value(true),
       "Print peak-fits to stdout."
       )
      ("out-dir", po::value<string>(&output_path)->default_value(""),
       "The directory to write peak CSV files (and optionally) N42 files to; if empty, CSV files will not be written.")
      ("back-sub-file", po::value<string>(&background_sub_file)->default_value(""),
       "File to use as the background, to perform a live-time-normalized, hard-background-subtraction with (currently must be single record).")
      ("include-nonfit-peaks", po::value<bool>(&show_nonfit_peaks)->default_value(false),
       "Include peaks that are not fit into the output CSV peak results."
       )
      ;
      
      
      string drf_file, drf_name;
      po::options_description activity_cl_desc("Activity-fit options", term_width, min_description_length);
      activity_cl_desc.add_options()
      ("drf-file",  po::value<string>(&drf_file)->default_value(""),
                    "Path to a file containing DRF to use (overrides DRF contained in exemplar N42 file)")
      ("drf-name",  po::value<string>(&drf_name)->default_value(""),
                    "The name of a DRF to use; either within the file specified by 'drf-file'"
                    " (which may have multiple DRFs)"
                    ", or built-in, or previously used in InterSpec "
                    "(overrides DRF contained in exemplar N42 file).")
      ("background-sample-nums", po::value<std::string>(&background_samples),
       "The sample numbers from the background file to use; if left empty will try to determine, and fail if not unique.\n\t"
       "Only applicable if the background subtraction file is specified."
       )
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
      
      BatchActivity::BatchActivityFitOptions options;  //derived from BatchPeak::BatchPeakFitOptions
      options.to_stdout = output_stdout;
      options.refit_energy_cal = refit_energy_cal;
      options.use_exemplar_energy_cal = use_exemplar_energy_cal;
      options.write_n42_with_peaks = write_n42_with_peaks;
      options.show_nonfit_peaks = show_nonfit_peaks;
      options.output_dir = output_path;
      options.background_subtract_file = background_sub_file;
      options.background_subtract_samples = background_sample_nums;
      
      if( batch_peak_fit )
      {
        BatchPeak::fit_peaks_in_files( exemplar_path, exemplar_sample_nums, 
                                      expanded_input_files, options );
      }//if( batch_peak_fit )
      
      if( batch_act_fit )
      {
        options.drf_override = BatchActivity::init_drf_from_name( drf_file, drf_name ); 
        
        //cerr << "batch-act-fit no implemented yet" << endl;
        //throw runtime_error( "Not implemented yet" );
        
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
