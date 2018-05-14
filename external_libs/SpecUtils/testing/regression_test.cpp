/* SpecUtils: a library to parse, save, and manipulate gamma spectrum data files.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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


#include "SpecUtils_config.h"

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>

#include <boost/program_options.hpp>

#include "SpecUtils/UtilityFunctions.h"
#include "SpecUtils/SpectrumDataStructs.h"

#if( !PERFORM_DEVELOPER_CHECKS )
#error PERFORM_DEVELOPER_CHECKS must be enabled to compile regression_test
#endif

using namespace std;
using boost::filesystem::path;
namespace po = boost::program_options;

//g_truth_n42_dirname: name of the sub-directory that truth N42 files will be
//  stored in.
const string g_truth_n42_dirname = "truth_n42s";

//g_automated_mode: specifies whether being ran in an automated mode, or with
//  user interaction.
bool g_automated_mode = false;

//g_sub_test_dir: specifies directory containing file types to be tested.  All
//  files if its empty.
string g_sub_test_dir = "";

//g_parse_time_filename: the name of the file written to the test base directory
//  that stores the parse times for files.  Not added to SVN.
const string g_parse_time_filename = "parsetimes.txt";

//open_spec_file(): uses the OS X 'open' command to open the spectrum file
//  with InterSpec running on localport:8080.
void open_spec_file( const path &p );

//open_spec_file_in_textmate(): : uses the OS X 'open' command to open the
//  spectrum file in textmate
void open_spec_file_in_textmate( const path &p );

//open_directory(): uses OS X 'open' command to open a finder window for the
//  specified directory; if a file is passed in, its parent directory is opened.
void open_directory( const path &p );

//handle_no_truth_files(): interactively creates truth files, based on prompting
//  user what actions should be taken (so they can decide it a truth file should
//  be created or not).
void handle_no_truth_files( const string basedir );

//FixableErrors: types of errors that are non-fatal for parsing that might be
//  expected to change over time, and hence you might want the truth level
//  information to be updated.
enum FixableErrors
{
  UuidError,
  NumFixableErrors
};//enum FixableErrors

//get_compare_error_type(): searches the what() given by the exception thrown
//  by MeasurmentInfo::equalEnough() to determine error type.  Returns
//  NumFixableErrors if it wasnt recognized or is not fixable.
FixableErrors get_compare_error_type( const string &msg );

//check_files_with_truth_n42(): check files that have a cooresponding truth n42
//  file, to be sure that the original and truth n42 files parse to similar
//  results using the 'equalEnough(...)' test.
void check_files_with_truth_n42( const string basedir );

//check_serialization_to_n42(): for test files (that have truth N42 files) this
//  function tests that the oringal file can be read in, written out to a 2011
//  N42, and then read back in and ensured its equalEnough() to the original.
void check_serialization_to_n42( const string basedir );

//add_truth_n42(): adds a truth n42 file for the MeasurementInfo and path
//  passed in.  Will fail if a truth n42 file already exists, unless force is
//  specified to be true.  Checks the created n42 file to be sure it can be read
//  back in and pass the 'equalEnough(...)' test, otherwise wont add truth n42
//  file.  Will add resulting added file (and possibly directory) to SVN.
//  Returns true if the truth N42 file was created.
bool add_truth_n42( const MeasurementInfo &m, const path &p, const bool force );

//check_parse_time(): compares the parse times of files with truth n42 files
//  against previous parse times.  Parses the file 10 times and takes the
//  minimum CPU time run as the parse time.
void check_parse_time( const string basedir );

//url_encode(): encodes a string appropriately to be passed as a URL argument.
string url_encode( const string &value );

//print_summary(): prints a reasonably brief summary to the provided stream.
void print_summary( const MeasurementInfo &info, std::ostream &out );

//print_one_line_summary(): prints a single line summary to the provided stream.
void print_one_line_summary( const Measurement &info, std::ostream &out );

//candidate_test_files(): return all candidate files, regardless if they have a
//  matching truth_n42.
vector<path> candidate_test_files( const string basedir );

//candidates_with_truth_n42_files(): returns only candidate files that have
//  truth information as well.
vector<path> candidates_with_truth_n42_files( const string basedir );

//truth_n42_files(): return all truth_n42 files.
vector<path> truth_n42_files( const string basedir );

//candidates_without_truth_n42_files(): return files that do not have truthfiles
vector<path> candidates_without_truth_n42_files( const string basedir );


int main( int argc, char **argv )
{
  //TODO: add test theat ensures serializing to N42 and back again pass
  //      equalEnough(...) test.
  
  //test_base_directory: the directory where the test file structure is based.
  string test_base_directory = "/Users/wcjohns/rad_ana/InterSpec/testing/file_format_test_spectra";

  
  po::options_description desc("Allowed options");
  
  desc.add_options()
  ( "help,h", "produce help message")
  ( "batch,b",
    po::value<bool>(&g_automated_mode)->default_value(false),
    "Run in non-interactive automated test mode.")
  ( "basedir,d",
    po::value<string>(&test_base_directory)->default_value(test_base_directory),
    "Directory where the test files are located.")
  ( "subdir,s",
    po::value<string>(&g_sub_test_dir)->default_value(""),
    "Sub-directory in 'basedir' of files to test.")
  ( "action,a",
   po::value< vector<string> >(),
   "Action to perform. Either 'n42test', 'regression' (or equivalently 'test'),"
   " 'addfiles', or 'timing'.  If blank defaults to 'test' if in automated mode"
   ", or 'n42test', 'addfiles', 'test', and 'timing' otherwise." )
  ;
  
  po::variables_map vm;
  
  try
  {
    po::store( po::parse_command_line( argc, argv, desc ), vm );
    po::notify( vm );
  }catch( std::exception &e )
  {
    cerr << "Invalid command line argument\n\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  
  if( vm.count( "help" ) )
  {
    cout << desc << "\n";
    return EXIT_SUCCESS;
  }

  
  if( !UtilityFunctions::is_directory( test_base_directory ) )
  {
    cerr << "Base directory '" << test_base_directory << "' is not a valid "
         << "directory\n";
    return EXIT_FAILURE;
  }//if( base directory is invalid )
  
  cout << "File base directory: '" << test_base_directory << "'\n";
  
  if( g_automated_mode )
    cout << "Testing in automated mode\n";
  else
    cout << "Testing in interactive mode\n";
  
  if( !g_sub_test_dir.empty() )
  {
    test_base_directory
       = UtilityFunctions::append_path( test_base_directory, g_sub_test_dir );
    if( !UtilityFunctions::is_directory(test_base_directory) )
    {
      cerr << "Test sub directory '" << g_sub_test_dir << "' is not a valid"
           << " directory\n";
      return EXIT_FAILURE;
    }//if( sub directory is invalid )
    
    cout << "Only testing files in the '" << g_sub_test_dir
         << "' subdirectory\n";
  }//if( !g_sub_test_dir.empty() )
  
  vector<string> actions;
  
  if( vm.count("action") )
    actions = vm["action"].as< vector<string> >();
  
  if( actions.empty() )
  {
    if( g_automated_mode )
      actions.push_back( "test" );
    else
      actions = {"addfiles", "test", "timing", "n42test"};
  }//if( vm.count("action") ) / else
  
  for( string action : actions )
  {
    if( action == "n42test" )
      check_serialization_to_n42( test_base_directory );
    else if( action == "regression" || action == "test" )
      check_files_with_truth_n42( test_base_directory );
    else if( action == "addfiles" )
      handle_no_truth_files( test_base_directory );
    else if( action == "timing" )
      check_parse_time( test_base_directory );
    else
    {
      cerr << "Invalid action type '" << action << "', valid options are "
           << "'n42test', 'regression', 'test', 'addfiles', 'timing', or blank\n";
      return EXIT_FAILURE;
    }
  }//for( string action : actions )
  
  return EXIT_SUCCESS;
}//int main( int argc, char **argv )



void open_spec_file( const path &p )
{
  const string command = "open http://localhost:8080/?specfilename="
                   + url_encode( p.string<string>() );
  system( command.c_str() );
  
//The following wine comand crashes Peak Easy, but would be nice.
//  command = "wine '/Users/wcjohns/Desktop/triage_spectra/PeakEasy 4.74/PeakEasy 4.74.exe' '" + path.string<string>() + "' &"
//  system( command.c_str() );
}//void open_spec_file( const path &p )



void open_spec_file_in_textmate( const path &p )
{
  //const string command = "open -a TextMate '" + p.string<string>() + "'";
  const string command = "/usr/local/bin/code '" + p.string<string>() + "'";

  system( command.c_str() );
}//void open_spec_file_in_textmate( const path &p )



void open_directory( const path &p )
{
  string command = "open '";
  
  if( !boost::filesystem::is_regular_file(p) )
    command += p.string<string>() + "'";
  else
    command += p.parent_path().string<string>() + "'";
  
  system( command.c_str() );
}//void open_directory( const path &p )



void check_parse_time( const string basedir )
{
  const int ntimes_parse = 10;
  map<path,double> cpu_parse_times, wall_parse_times;
  const vector<path> with_truth = candidates_with_truth_n42_files( basedir );
  
  const boost::posix_time::ptime start_time
                             = boost::posix_time::microsec_clock::local_time();
  
  for( const path &fpath : with_truth )
  {
    const string filename = fpath.string<string>();
    const string extention = fpath.extension().string<string>();
    
    for( int i = 0; i < ntimes_parse; ++i )
    {
      MeasurementInfo info;
      
      const double orig_wall_time = UtilityFunctions::get_wall_time();
      const double orig_cpu_time = UtilityFunctions::get_cpu_time();
    
      const bool parsed = info.load_file( filename, kAutoParser, extention );
    
      const double final_cpu_time = UtilityFunctions::get_cpu_time();
      const double final_wall_time = UtilityFunctions::get_wall_time();
    
      if( parsed && orig_cpu_time > 0.0 && final_cpu_time > 0.0 )
      {
        double cpu_dt = final_cpu_time - orig_cpu_time;
        double wall_dt = final_wall_time - orig_wall_time;
        const map<path,double>::const_iterator pos = cpu_parse_times.find( fpath );
        if( pos != cpu_parse_times.end() && (pos->second < cpu_dt) )
        {
          cpu_dt = pos->second;
          wall_dt = wall_parse_times[pos->first];
        }
        
        cpu_parse_times[fpath] = cpu_dt;
        wall_parse_times[fpath] = wall_dt;
      }//if( parsed && orig_cpu_time > 0.0 && final_cpu_time > 0.0 )
    }//for( int i = 0; i < ntimes_parse; ++i )
  }//for( const path &fpath : with_truth )
  
  string prevtimestr;
  map<path,double> previous_cpu_parse_times, previous_wall_parse_times;
  
  const string timingname
              = UtilityFunctions::append_path( basedir, g_parse_time_filename );
  
  {//begin read parsetimes
    ifstream parsetimes( timingname.c_str(), ios::in | ios::binary );
    if( parsetimes.is_open() && parsetimes.good() )
    {
      string filename;
      UtilityFunctions::safe_get_line( parsetimes, prevtimestr );
      
      while( UtilityFunctions::safe_get_line( parsetimes, filename ) )
      {
        if( filename.size() == 0 )
          continue;
        
        string times;
        UtilityFunctions::safe_get_line( parsetimes, times );
        
        double cputime, walltime;
        
        const int n = sscanf( times.c_str(), "%lf %lf", &cputime, &walltime );
        if( n != 2 )
        {
          cerr << "Error reading times for file '" << filename << "'\n"
               << "Stopping parsing timing file.\n";
          break;
        }//if( n != 2 )
        
        previous_cpu_parse_times[filename] = cputime;
        previous_wall_parse_times[filename] = walltime;
      }//while( UtilityFunctions::safeGetline( parsetimes, line ) )
    }//if( parsetimes.is_open() && parsetimes.good() )
  }//end read parsetimes
  
  
  //print comparison
  double prev_cpu_total = 0.0, prev_wall_total = 0.0;
  double current_cpu_total = 0.0, current_wall_total = 0.0;
  
  bool previous_had_all = true;
  cout << "Previous parse time: " << prevtimestr << endl;
  for( map<path,double>::iterator i = cpu_parse_times.begin();
      i != cpu_parse_times.end(); ++i )
  {
    const double cputime = i->second;
    const double walltime = wall_parse_times[i->first];
    
    current_cpu_total += cputime;
    current_wall_total +=walltime;
    
    string name = i->first.filename().string<string>();
    if( name.size() > 30 )
      name = name.substr( 0, 27 ) + "...";
    cout << setw(31) << std::left << name << ": {cpu: "
         << setprecision(6) << std::fixed << cputime << ", wall: "
         << setprecision(6) << walltime << "}\n";
    
    if( previous_cpu_parse_times.count(i->first) )
    {
      prev_cpu_total += previous_cpu_parse_times[i->first];
      prev_wall_total += previous_wall_parse_times[i->first];
      
      cout << "                      previous : {cpu: "
           << setprecision(6) << previous_cpu_parse_times[i->first]
           << ", wall: "
           << setprecision(6) << previous_wall_parse_times[i->first] << "}\n";
    }else
    {
      previous_had_all = false;
      cout << "          no previous          \n";
    }
    cout << "\n";
  }//for( loop over parse times )
  
  cout << "Current total  : {cpu: " << setprecision(6) << current_cpu_total
       << ", wall: " << setprecision(6) << prev_wall_total << "}\n";
  
  if( previous_had_all )
    cout << "Previous total : {cpu: " << setprecision(6) << prev_cpu_total
         << ", wall: " << setprecision(6) << prev_wall_total << "}\n\n";
  else
    cout << "Did not have previous timings for all the files\n\n";
  
  //Decide if we should save the current results.
  char action = (g_automated_mode ? 'n' : '\0');
  while( action != 'y' && action != 'n' )
  {
    cout << "Would you like to save these latest timings? y/n\n";
    cin >> action;
  }//while( user hasnt entered valid answer )
  
  if( action == 'y' )
  {
    ofstream file( timingname.c_str(), ios::out | ios::binary );
    if( file.is_open() )
    {
      file << UtilityFunctions::to_extended_iso_string( start_time ) << endl;
      for( map<path,double>::iterator i = cpu_parse_times.begin();
           i != cpu_parse_times.end(); ++i )
      {
        file << i->first.string<string>() << endl << i->second
            << " " << wall_parse_times[i->first] << endl;
      }
      
      cout << "Saved timings to '" << timingname << "'" << endl;
    }else
    {
      cerr << "Failed to open '" << timingname << "' for writing times\n";
    }
  }//if( action == 'y' )
  
}//void compare_parse_time( const string basedir )



FixableErrors get_compare_error_type( const string &msg )
{
  if( UtilityFunctions::icontains( msg, "UUID of LHS" ) )
    return UuidError;
  
  return NumFixableErrors;
}//FixableErrors get_compare_error_type( const string &msg )



void check_files_with_truth_n42( const string basedir )
{
  int initial = 0, initial_parsed = 0, failed_original_parsed = 0;
  int failed_truth_parsed = 0, initial_with_truth = 0, passed_tests = 0;
  int failed_tests = 0, updated_truths = 0, truths_failed_to_update = 0;
  
  const vector<path> with_truth = candidates_with_truth_n42_files( basedir );
  
  map<path,double> parse_times;
  
  for( const path &fpath : with_truth )
  {
    ++initial;
    
    const string filename = fpath.filename().string<string>();
    const string originalpath = fpath.string<string>();
    const string originalext = fpath.extension().string<string>();
    
    MeasurementInfo original;
    
    const bool originalstatus
                = original.load_file( originalpath, kAutoParser, originalext );
    
    if( !originalstatus )
    {
      ++failed_original_parsed;
      cerr << "Failed to parse original file " << fpath << "\n\tskipping; type 'c' and enter to continue.\n\n";
      
      char awk = (g_automated_mode ? 'c' : 'n');
      while( awk != 'c' )
        cin >> awk;
      
      continue;
    }
    
    ++initial_parsed;
    
    const path tpath = fpath.parent_path() / g_truth_n42_dirname
                       / (filename + ".n42");
    
    if( !boost::filesystem::is_regular_file(tpath) )
    {
      cerr << "Fatal error: " << fpath << " doesnt have truth file at "
           << tpath << "\n\n";
      exit( EXIT_FAILURE );
    }
    
    MeasurementInfo truth;
    const bool truthstat = truth.load_file( tpath.c_str(), K2011ICD1Parser, "" );
    
    if( !truthstat )
    {
      ++failed_truth_parsed;
      cerr << "Failed to parse truth file " << tpath << "\n\tskipping.\n\n";
      continue;
    }
    
    ++initial_with_truth;
    
    truth.set_filename( original.filename() );
    
    try
    {
      MeasurementInfo::equalEnough( original, truth );
      ++passed_tests;
    }catch( std::exception &e )
    {
      ++failed_tests;
      
      const string description = e.what();
      
      cerr << "\n" << fpath << "\nfailed comparison with previous parsing: "
           << description
           << "\n\t(Current parse is LHS, previous parse is RHS)\n"
           << "\n\tWhat would like to do?"
           << "\n";
      
      const FixableErrors errortype = get_compare_error_type( description );
      
      char action = (g_automated_mode ? 'i' : '\0');
      while( (action != 'i') && (action != 'u') )
      {
        cout << "\ti: ignore\n"
             << "\to: open original file\n"
             << "\tt: open truth n42\n"
             << "\td: open contiaing directory\n"
             << "\tp: print summary of current parsing\n"
             << "\tq: print summary of truth\n"
             << "\tu: update truth\n";
        if( errortype != NumFixableErrors )
          cout << "\ts: set new error value to old parsing and try again\n";
        
        
        cin >> action;
        
        switch( action )
        {
          case 'i':                                         break;
          case 'o': open_spec_file( fpath );                break;
          case 't': open_spec_file( tpath );                break;
          case 'd': open_directory( fpath );                break;
          case 'p': print_summary( original, cout );        break;
          case 'q': print_summary( truth, cout );           break;
          
          case 'u':
            if( add_truth_n42( original, fpath, true ) )
            {
              ++updated_truths;
              cout << "\nUpdated truth info file.\n\n";
            }else
            {
              ++truths_failed_to_update;
              cout << "\nFailed to updated truth info file.\n\n";
            }
          break;
            
          case 's':
          {
            try
            {
              switch( errortype )
              {
                case UuidError:
                  truth.set_uuid( original.uuid() );
                break;
                  
                case NumFixableErrors:
                  assert( 0 );
                  break;
              }//switch( errortype )
              
              MeasurementInfo::equalEnough( original, truth );
              
              cout << "\nFixing the issue allowed the comparison test to pass."
                   << "\nWould you like to update the truth level information?"
                   << " (y/n)\n";
              cin >> action;
              if( action == 'y' )
              {
                if( add_truth_n42( original, fpath, true ) )
                {
                  ++updated_truths;
                  cout << "\nUpdated truth info file.\n\n";
                }else
                {
                  ++truths_failed_to_update;
                  cout << "\nFailed to updated truth info file.\n\n";
                }
                action = 'u';
              }//if( action == 'y' )
            }catch( std::exception &e )
            {
              cout << "\nAfter fixing error, there was another error: \n\t"
                   << e.what() << "\nNot updating truth information.\n\n";
            }
            
            break;
          }//case 's':
            
          default:                                          break;
        }//switch( action )
      }//while( action != 'i' && action != 'u' )
    }//try / catch
    
  }//for( const path &path : with_truth )
  
  cout << "Of the " << initial << " initial test files " << initial_parsed
       << " were parsable (" << failed_original_parsed << " failed).\n"
       << failed_truth_parsed << " of the truth N42 files failed to parse.\n"
       << "Of the " << initial_with_truth << " parsible original files with "
       << "valid truth N42 files: \n"
       << "\t" << passed_tests << " passed comparison\n"
       << "\t" << failed_tests << " failed comaprison, with "
       << updated_truths << " truth N42 files updated.\n";
  if( truths_failed_to_update )
    cerr << truths_failed_to_update << " truth n42 files failed to update!\n";
}//void check_files_with_truth_n42( const string basedir )



void check_serialization_to_n42( const string basedir )
{
  size_t ninitial = 0, nOrigFileFailParse = 0,
         nFailToSerialize = 0, nSerializedFileFailParse = 0,
         npassed = 0, nfailed = 0;
  vector<path> failedcompare;
  const path tempdir = UtilityFunctions::temp_dir();
  const vector<path> with_truth = candidates_with_truth_n42_files( basedir );
  for( const path &fpath : with_truth )
  {
    ++ninitial;
    
    const string filename = fpath.filename().string<string>();
    const string originalpath = fpath.string<string>();
    const string originalext = fpath.extension().string<string>();
    
    MeasurementInfo info;
    
    bool status = info.load_file( originalpath, kAutoParser, originalext );

    if( !status )
    {
      ++nOrigFileFailParse;
      cerr << "N42 Serialization Test: Failed to parse input file " << fpath
           << "\n\n";
      continue;
    }//if( !status )
    
    
    const string tempname = UtilityFunctions::temp_file_name( filename, tempdir.c_str() );
    
    {//Begin codeblock to serialize to temporary file
      ofstream output( tempname.c_str() );
      if( !output )
      {
        ++nFailToSerialize;
        cerr << "N42 Serialization Test: Couldnt open temporary file "
             << tempname << "\n\n";
        UtilityFunctions::remove_file( tempname.c_str() );
        continue;
      }
    
      status = info.write_2011_N42( output );

      if( !status )
      {
        ++nFailToSerialize;
        cerr << "N42 Serialization Test: Couldnt serialize " << fpath
             << " to temp file " << tempname << "\n\n";
        UtilityFunctions::remove_file( tempname.c_str() );
        continue;
      }//if( !status )
    }//End codeblock to serialize to temporary file
    
    MeasurementInfo reread;
    status = reread.load_file( tempname.c_str(), K2011ICD1Parser, "" );
    
    if( !status )
    {
      ++nSerializedFileFailParse;
      cerr << "N42 Serialization Test: Couldnt parse serialized N42 file for "
           << fpath << "\n\n";
      UtilityFunctions::remove_file( tempname.c_str() );
      continue;
    }//if( !status )
    
    reread.set_filename( info.filename() );
    
    try
    {
      MeasurementInfo::equalEnough( info, reread );
      ++npassed;
    }catch( std::exception &e )
    {
      const string error_msg = e.what();
      ++nfailed;
      failedcompare.push_back( fpath );
      cerr << "N42 Serialization Test: comparison test for " << fpath
           << " failed with error:\n\t" << error_msg << "\n"
           << "\t(LHS is original parse, RHS is re-ad back in)\n\n";
      
      if( UtilityFunctions::contains(error_msg, "MeasurementInfo: Number of remarks in LHS") )
      {
        for( const string r : info.remarks() )
          cout << "\t\tLHS remark: '"  << r << "'" << endl;
        for( const string r : reread.remarks() )
          cout << "\t\tRHS remark: '"  << r << "'" << endl;
      }
      
    }//try / catch
    
    UtilityFunctions::remove_file( tempname.c_str() );
  }//for( const path &fpath : with_truth )
  
  
  cout << "N42 Serialization Test Results:\n"
       << "\tNumber of input files: " << ninitial << "\n"
       << "\tNumber of input files that failed to parse: " << nOrigFileFailParse
       << "\n\tNumber of files that failed to serialize to N42: "
       << nFailToSerialize
       << "\n\tNumber of serialed files that couldnt be parsed: "
       << nSerializedFileFailParse
       << "\n\tNumber of files that failed comparison: " << nfailed
       << "\n\tNumber of files that passed comparison: " << npassed
       << "\n\n";
  
  if( failedcompare.size() )
  {
    cout << "Files failing comparison:\n";
    for( const path &p : failedcompare )
      cout << "\t" << p << endl;
    cout << endl << endl;
  }
  
  if( !g_automated_mode && (nFailToSerialize || nSerializedFileFailParse || nfailed) )
  {
    char dummy;
    
    
    
    cout << "There was an error, enter 'x' to exit the app, or any other key"
         << " to continue.\n";
    cin >> dummy;
    
    if( dummy == 'x' )
      exit( EXIT_FAILURE );
  }//if( !g_automated_mode )
  
}//void check_serialization_to_n42( const string basedir )



bool add_truth_n42( const MeasurementInfo &info, const path &p,
                    const bool force )
{
  const path parentdir = p.parent_path();
  const path truthdir = parentdir / g_truth_n42_dirname;
  const path truth_n42 = truthdir / (p.filename().string<string>() + ".n42");
  path old_n42;
  const bool prevexist = boost::filesystem::is_regular_file( truth_n42 );
  if( !force && prevexist )
  {
    cerr << "File " << truth_n42 << " already exits, not re-creating\n";
    return false;
  }if( prevexist )
  {
    old_n42 = truth_n42.string<string>() + ".prev";
    boost::filesystem::rename( truth_n42, old_n42 );
  }
  
  try
  {
    if( !boost::filesystem::is_directory( truthdir ) )
    {
      try
      {
        boost::filesystem::create_directory( truthdir );
        const string command = "svn add '" + truthdir.string<string>() + "'";
        const int val = system( command.c_str() );
        if( val != 0 )
          cerr << "\n\nThere may have been an issue adding " << truthdir
               << " to the SVN repository.  Return code " << val << "\n";
      }catch( std::exception &e )
      {
        throw runtime_error( "Couldnt create directory "
                        + truthdir.string<string>() + ", so skipping file" );
      }
    }//if( !boost::filesystem::is_directory( truthdir ) )
  
    {//begin write file
      ofstream output( truth_n42.c_str() );
      if( !output )
        throw runtime_error( "Couldnt create file " + truth_n42.string<string>()
                             + ", so skipping file" );
  
      const bool status = info.write_2011_N42( output );
    
      if( !status )
        throw runtime_error( "Failed to write to file "
                        + truth_n42.string<string>() + ", so skipping file" );
    }//end write files
  
    
    MeasurementInfo reloadedinfo;
    const bool reloadstatus
                 = reloadedinfo.load_file( truth_n42.c_str(), K2011ICD1Parser, "" );
    if( !reloadstatus )
      throw runtime_error( "Failed to read in written n42 file" );
  
    reloadedinfo.set_filename( info.filename() );
    
    try
    {
      MeasurementInfo::equalEnough( info, reloadedinfo );
    }catch( std::exception &e )
    {
      char option = '\0';
      while( option != 'n' && option != 'y' )
      {
        cerr << "Writing " << truth_n42 << " to a file and then reading back in"
             << " resulted in\n\t" << e.what() << endl
             << "\t(LHS is original parse, RHS is re-ad back in)\n\n"
             << "What would you like to do:\n"
             << "\tn: skip this file\n"
             << "\ty: use this file anyway\n";
        cin >> option;
      }//while( option != 'n' && option != 'y' )
      
      if( option == 'n' )
        throw runtime_error( "Failed to make the MeasurementInfo--->N42"
                             "--->MeasurementInfo round trip" );
    }//try / catch for equalEnoughs
    
    if( !old_n42.empty() )
      try{ boost::filesystem::remove( old_n42 ); }catch(...){}
    
    const string command = "svn add '" + truth_n42.string<string>() + "'";
    const int rval = system( command.c_str() );
    if( rval != 0 )
      cerr << "\n\nThere may have been an issue adding " << truth_n42
           << " to the SVN repository.  Return code " << rval << "\n";
    
    cout << "Added truth n42 file: " << truth_n42 << "\n\n\n";
    
  }catch( std::exception &e )
  {
    cerr << e.what() << "\n\tskipping writing file\n";
    if( !force && boost::filesystem::is_regular_file( truth_n42 ) )
    {
      try{ boost::filesystem::remove( truth_n42 ); }catch(...){}
      if( !old_n42.empty() )
        try{ boost::filesystem::rename( old_n42, truth_n42 ); }catch(...){}
    }//if( a new file was written )
    
    return false;
  }//try / catch
  
  return true;
}//void add_truth_n42( MeasurementInfo &info, path &p )



void handle_no_truth_files( const string basedir )
{
  size_t nfailed_parse = 0, nadded = 0, nfail_add = 0, nignored = 0;
  const vector<path> no_truth = candidates_without_truth_n42_files( basedir );
  
  cout << "\nFound " << no_truth.size() << " files without truth N42 files\n\n";
  
  for( const path &path : no_truth )
  {
    const string filenamestr = path.string<string>();
    const string extention = path.extension().string<string>();
    
    MeasurementInfo info;
    const bool status = info.load_file( filenamestr, kAutoParser, extention );
    
    if( !status )
    {
      ++nfailed_parse;
      cerr << "\nFailed to parse file " << path << ", type 'c' and hit enter to continue" << endl;
      
      char awk = (g_automated_mode ? 'c' : 'n');
      while( awk != 'c' )
        cin >> awk;

      
      continue;
    }
    
    char action = (g_automated_mode ? 'c' : '\0');
    
    while( (action != 'i') && (action != 'c') )
    {
      cout << "File " << path << " does not have a truth N42 file, would you "
           << "like to:\n"
           << "\to: open\n"
           << "\tt: open file in textmate\n"
           << "\td: open containing directory\n"
           << "\tp: print summary\n"
           << "\tc: create truth N42 file\n"
           << "\ti: ignore file ?\n";
      cin >> action;
      
      switch( action )
      {
        case 'i':                                     break;
        case 'c':                                     break;
        case 'o': open_spec_file( path );             break;
        case 't': open_spec_file_in_textmate( path ); break;
        case 'd': open_directory( path );             break;
        case 'p': print_summary( info, cout );        break;
        default:                                      break;
      }//switch( action )
    }//while( (action != 'i') && (action != 'c') )
    
    if( action == 'c' )
    {
      const bool added = add_truth_n42( info, path, false );
      nadded += (added ? 1 : 0);
      nfail_add += (added ? 0 : 1);
    }

    if( action == 'i' )
      ++nignored;
  }//for( const path &path : no_truth )

  
  cout << "\n\nResults of trying to add truth N42 files:\n"
       << "\tAdded " << nadded << " truth N42 files.\n"
       << "\tFailed to add " << nfail_add << " truth N42 files due to N42 not parsing exactly like original.\n"
       << "\tIgnored " << nignored << " files.\n"
       << "\tFailed to parse " << nfailed_parse << " potential input files.\n";
}//void handle_no_truth_files( const string basedir )



void print_one_line_summary( const Measurement &meas, std::ostream &out )
{
  out << "Sample " << meas.sample_number() << " detector '"
      << meas.detector_name() << "', LT=" << meas.live_time()
      << ", RT=" << meas.real_time() << ", GammaSum=" << meas.gamma_count_sum();
  if( meas.contained_neutron() )
    out << ", NeutronSum=" << meas.neutron_counts_sum();
  else
    out << ", No neutron detector";
  
  switch( meas.source_type() )
  {
    case Measurement::Background:        out << ", Background";        break;
    case Measurement::Calibration:       out << ", Calibration";       break;
    case Measurement::Foreground:        out << ", Foreground";        break;
    case Measurement::IntrinsicActivity: out << ", IntrinsicActivity"; break;
    case Measurement::UnknownSourceType: out << ", UnknownSourceType"; break;
  }//switch( meas.source_type() )
  
  out << ", " << meas.start_time();
  
  if( meas.has_gps_info() )
    out << ", GPS(" << meas.latitude() << "," << meas.longitude() << "," << meas.position_time() << ")";
  
}//void print_one_line_summary( const Measurement &info, std::ostream &out )



void print_summary( const MeasurementInfo &info, std::ostream &out )
{
  vector< std::shared_ptr<const Measurement> > meass = info.measurements();

  const size_t ndet = info.detector_names().size();
  
  out << info.filename() << " successfully parsed to yeild " << meass.size()
      << " Measurements.\n"
      << "\tThere are " << ndet << " detectors: ";
  for( size_t i = 0; i < ndet; ++i )
    out << (i ? ", " : "") << info.detector_names()[i];
  out << "\n\tWith total live time " << info.gamma_live_time() << ", real time "
      << info.gamma_real_time() << ", and " << info.neutron_counts_sum()
      << " neutrons\n";
  
  if( info.manufacturer().size() )
    out << "\tmanufacturer: " << info.manufacturer() << "\n";
  if( info.instrument_model().size() )
    out << "\tinstrument_model: " << info.instrument_model() << "\n";
  if( info.instrument_id().size() )
    out << "\tinstrument_id (serial #): " << info.instrument_id() << "\n";
  if( info.uuid().size() )
    out << "\tuuid: " << info.uuid() << "\n";
  if( info.lane_number() > -1 )
    out << "\tlane_number: " << info.lane_number() << "\n";
  out << "\tAnd is " << (info.passthrough() ? "": "not ")
      << "passthrough/searchmode data.\n";
  
  //print out analysis info
  std::shared_ptr<const DetectorAnalysis> ana = info.detectors_analysis();
  if( !ana )
    out << "\tDoes not contain analysis results\n";
  else
    out << "\tContains analysis results with " << ana->results_.size() << "nuclides\n";


  for( size_t i = 0; i < meass.size(); ++i )
  {
    out << setw(4) << i << ": ";
    print_one_line_summary( *meass[i], out );
    out << "\n";
  }//for( size_t i = 0; i < meass.size(); ++i )
  
  out << "\n";
}//void print_summary( const MeasurementInfo &info, std::ostream out );



string url_encode( const string &value )
{
  //adapted from http://stackoverflow.com/questions/154536/encode-decode-urls-in-c
  ostringstream escaped;
  escaped.fill('0');
  escaped << hex;
  
  for( string::value_type c : value )
  {
    if (isalnum(c) || c == '-' || c == '_' || c == '.' || c == '~')
      escaped << c;
    else
      escaped << '%' << setw(2) << int((unsigned char) c);
  }//for( string::value_type c : value )
  
  return escaped.str();
}//string url_encode( const string &value )



vector<path> candidate_test_files( const string basedir )
{
  vector<path> filenames;
  
  const vector<string> allfiles = UtilityFunctions::recursive_ls( basedir );
  
  for( const string &filepath : allfiles )
  {
    const string filename = UtilityFunctions::filename(filepath);
    const string parentdir = path(filepath).parent_path().filename().string<string>();
    
    if( filename != "source.txt"
        && filename != g_parse_time_filename
        && parentdir != g_truth_n42_dirname
        && filename.size() && filename[0] != '.' )
      filenames.push_back( filepath );
  }//for( const path &filename : allfiles )
  
  return filenames;
}//vector<string> candidate_test_files( const string basedir )



vector<path> candidates_with_truth_n42_files( const string basedir )
{
  vector<path> results;
  const vector<path> truthfiles = truth_n42_files( basedir );
  const vector<path> candidates = candidate_test_files( basedir );
  
  for( const path &cand : candidates )
  {
    const string filename = cand.filename().string<string>() + ".n42";
    const path testpath = cand.parent_path() / g_truth_n42_dirname / filename;
    const vector<path>::const_iterator pos
                     = find( truthfiles.begin(), truthfiles.end(), testpath );
    if( pos != truthfiles.end() )
      results.push_back( cand );
  }//for( const path &cand : candidates )
  
  return results;
}//vector<path> candidates_with_truth_n42_files( const string basedir )



vector<path> truth_n42_files( const string basedir )
{
  vector<path> filenames;
  
  const vector<string> allfiles = UtilityFunctions::recursive_ls( basedir );
  
  for( const string &filestr : allfiles )
  {
    const string filename = path(filestr).filename().string<string>();
    const string parentdir = path(filestr).parent_path().filename().string<string>();
    
    if( filename != "source.txt"
        && filename != g_parse_time_filename
        && parentdir == g_truth_n42_dirname
        && filename.size() && filename[0] != '.'  )
      filenames.push_back( filestr );
  }//for( const path &filename : allfiles )
  
  return filenames;
}//vector<path> truth_n42_files( const string basedir )



vector<path> candidates_without_truth_n42_files( const string basedir )
{
  vector<path> results;
  const vector<path> truthfiles = truth_n42_files( basedir );
  const vector<path> candidates = candidate_test_files( basedir );
  
  for( const path &cand : candidates )
  {
    const string filename = cand.filename().string<string>() + ".n42";
    const path testpath = cand.parent_path() / g_truth_n42_dirname / filename;
    const vector<path>::const_iterator pos
                      = find( truthfiles.begin(), truthfiles.end(), testpath );
    if( pos == truthfiles.end() )
      results.push_back( cand );
  }//for( const path &cand : candidates )
  
  return results;
}//vector<path> candidates_without_truth_n42_files( const string basedir )

