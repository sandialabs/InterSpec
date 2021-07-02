//  Created by wcjohns on 20110322.
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

#include <string>
#include <iostream>

#include <Wt/WString>
#include <Wt/WApplication>
#include <Wt/WEnvironment>

#include "InterSpec/InterSpec.h"
#include "SpecUtils/Filesystem.h"
#include "InterSpec/InterSpecApp.h"
#include "SpecUtils/SerialToDetectorModel.h"
#include "InterSpec/DataBaseVersionUpgrade.h"

#if( ANDROID )
#include "android/AndroidUtils.hpp"
#endif


#if( BUILD_AS_OFFLINE_ANALYSIS_TEST_SUITE )
#include "InterSpec/SpectrumViewerTester.h"
#endif

#if( BUILD_AS_COMMAND_LINE_CODE_DEVELOPMENT )
#include "testing/developcode.h"
#endif

//Forward declaration
Wt::WApplication *createApplication( const Wt::WEnvironment &env );

#ifdef _WIN32
void getUtf8Args( int &argc, char ** &argv );
#endif

void processCustomArgs( int argc, char **argv );


int main( int argc, char **argv )
{
#ifdef _WIN32
  getUtf8Args( argc, argv );
#endif

  
#if( BUILD_AS_COMMAND_LINE_CODE_DEVELOPMENT )
  return developcode::run_development_code();
#endif
    
#if( BUILD_AS_OFFLINE_ANALYSIS_TEST_SUITE )
  SpectrumViewerTester::doOfflineTesting();
  return 1;
#endif
  
#if( ANDROID )
  AndroidUtils::androidbuf stdbuf( AndroidUtils::androidbuf::FromCout );
  AndroidUtils::androidbuf errbuf( AndroidUtils::androidbuf::FromCerr );
  AndroidUtils::set_anrdoid_cwd( argc, argv );
#endif //#if( ANDROID )
  
  std::cout << std::showbase << std::hex << "Running with Wt version "
            << WT_VERSION << std::dec << ", from executable compiled on "
            << __DATE__ << std::endl;
  
#if( PERFORM_DEVELOPER_CHECKS )
  std::cout << "Developer tests are being performed" << std::endl;
#endif

  std::cout << std::endl;
  
#if( WT_VERSION >= 0x3030300 )
  //Make it so WString defaults to assuming std::string or char * are UTF8
  //  encoded, rather than the system encoding.
  Wt::WString::setDefaultEncoding( Wt::UTF8 );
#endif
  
#if( BUILD_AS_LOCAL_SERVER )
  // For development we'll put in some default command line argument that assume the CWD is
  //  either the base InterSpec code directory, or the CMake build directory.
  std::string def_args[] = { argv[0], "--docroot=.", "--http-address=127.0.0.1",
    "--http-port=8080", "--config=./data/config/wt_config_localweb.xml", "--accesslog=-",
    "--no-compression"
  };
  const size_t default_argc = sizeof(def_args) / sizeof(def_args[0]);
  char *default_argv[default_argc] = { nullptr };
  
  if( argc == 1 )
  {
    std::cout << "Using default set of command line arguments\n" << std::endl;
    
    for( size_t i = 0; i < default_argc; ++i )
      default_argv[i] = &(def_args[i].front());
    
    argc = static_cast<int>( default_argc );
    argv = default_argv;
  }//if( no command line arguments were given ).
#endif //BUILD_AS_LOCAL_SERVER
  
  processCustomArgs( argc, argv );
  
  DataBaseVersionUpgrade::checkAndUpgradeVersion();
  
  return Wt::WRun( argc, argv, &createApplication );
}//int main( int argc, const char * argv[] )



Wt::WApplication *createApplication( const Wt::WEnvironment &env )
{
  return new InterSpecApp( env );
}// Wt::WApplication *createApplication(const Wt::WEnvironment& env)

#ifdef _WIN32

#define WIN32_LEAN_AND_MEAN 1
#include <windows.h>
#include <stdio.h>
#include <shellapi.h>

#include "SpecUtils/StringAlgo.h"

/** Get command line arguments encoded as UTF-8.
    This function just leaks the memory
 
 Note that environment variables are not in UTF-8, we could account for this
 similar to:
 wchar_t *wenvstrings = GetEnvironmentStringsW();
 ...
 */
void getUtf8Args( int &argc, char ** &argv )
{
  LPWSTR *argvw = CommandLineToArgvW( GetCommandLineW(), &argc );
  if( !argvw )
  {
    std::cout << "CommandLineToArgvW failed - good luck" << std::endl;
    return ;
  }
  
  argv = (char **)malloc(sizeof(char *)*argc);
    
  for( int i = 0; i < argc; ++i)
  {
    printf("Argument: %d: %ws\n", i, argvw[i]);
  
    const std::string asutf8 = SpecUtils::convert_from_utf16_to_utf8( argvw[i] );
    argv[i] = (char *)malloc( sizeof(char)*(asutf8.size()+1) );
    strcpy( argv[i], asutf8.c_str() );
  }//for( int i = 0; i < argc; ++i)

  // Free memory allocated for CommandLineToArgvW arguments.
  LocalFree(argvw);
}//void processCustomArgs()
#endif


void processCustomArgs( int argc, char **argv )
{
  bool set_serial_num_file = false;
  for( int i = 1; i < (argc-1); ++i )
  {
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || (BUILD_AS_LOCAL_SERVER && (defined(WIN32) || defined(__APPLE__)) ) )
    if( argv[i] == std::string("--userdatadir") )
    {
      try
      {
        InterSpec::setWritableDataDirectory( argv[i+1] );
        const std::vector<std::string> serial_db = SpecUtils::ls_files_in_directory( argv[i+1], "serial_to_model.csv" );
        if( !serial_db.empty() )
        {
          SerialToDetectorModel::set_detector_model_input_csv( serial_db[0] );
          set_serial_num_file = true;
        }
      }catch( std::exception & )
      {
        std::cerr << "Invalid userdatadir ('" << argv[i+1] << "') specified"
                  << std::endl;
      }
    }//if( argv[i] == std::string("--userdatadir") )
#endif  //if( not a webapp )
  }//for( int i = 1; i < (argc-1); ++i )
  
  
#if( BUILD_AS_LOCAL_SERVER )
  if( !set_serial_num_file )
  {
    const char * const searchdirs[] = { "data", "data_ouo" };
    for( const char *dir : searchdirs )
    {
      const std::vector<std::string> serial_db = SpecUtils::ls_files_in_directory( dir, "serial_to_model.csv" );
      if( !serial_db.empty() )
      {
        SerialToDetectorModel::set_detector_model_input_csv( serial_db[0] );
        break;
      }
    }
  }//if( !set_serial_num_file )
#endif
}//void processCustomArgs( int argc, char **argv )
