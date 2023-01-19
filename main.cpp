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

#if( USE_SQLITE3_DB )
#include "InterSpec/DataBaseUtils.h"
#endif


#if( BUILD_AS_COMMAND_LINE_CODE_DEVELOPMENT )
#include "testing/developcode.h"
#endif


#include "InterSpec/QRSpecDev.h"

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN 1
#include <windows.h>
#include <stdio.h>
#include <direct.h>
#include <shellapi.h>

#include "SpecUtils/StringAlgo.h"
#endif


//Forward declaration
Wt::WApplication *createApplication( const Wt::WEnvironment &env );

#ifdef _WIN32
void getUtf8Args( int &argc, char ** &argv );
#endif

void processCustomArgs( int argc, char **argv );


int main( int argc, char **argv )
{
  return QRSpecDev::dev_code();
  
#ifdef _WIN32
  getUtf8Args( argc, argv );
#endif
  
#if( BUILD_AS_COMMAND_LINE_CODE_DEVELOPMENT )
  return developcode::run_development_code();
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

#if(_WIN32)
    // I cant get MSVC to set CWD to anywhere besides InterSpec/out/build/x64-Debug/,
    //  so we'll look for our resources up to three levels up
    const std::string targetfile = "InterSpec_resources/InterSpec.css";
    std::string pardir = "";
    while (pardir.size() < 9)
    {
      const std::string testfile = SpecUtils::append_path(pardir, targetfile);
      if (SpecUtils::is_file(testfile))
      {
        if(!pardir.empty())
        {
          std::cout << "Will change cwd to '" << pardir << "'." << std::endl;
          if (_chdir(pardir.c_str()))
            std::cerr << "Failed to change CWD to '" << pardir << "'" << std::endl;
        }
        break;
      }
      pardir = SpecUtils::append_path(pardir, "..");
    }//while (pardir.size() < 6)
#endif
    
    // If there is a "user_data" directory in the CWD, we'll set this as the writeable data
    //  directory to simulate desktop app behavior of saving DRFs and similar
    const std::string cwd = SpecUtils::get_working_path();
    std::cout << "cwd='" << cwd << "'" << std::endl;

    const std::string dev_user_data = SpecUtils::append_path( cwd, "user_data" );
    if( SpecUtils::is_directory( dev_user_data ) )
      InterSpec::setWritableDataDirectory( dev_user_data );
    else
      std::cerr << "No '" << dev_user_data << "' - not setting writeable data directory.\n";
  }//if( no command line arguments were given ).
#endif //BUILD_AS_LOCAL_SERVER
  
  processCustomArgs( argc, argv );
  
  DataBaseVersionUpgrade::checkAndUpgradeVersion();
  
  // TODO: switch to using InterSpecServer::startServer(), etc
  return Wt::WRun( argc, argv, &createApplication );
}//int main( int argc, const char * argv[] )



Wt::WApplication *createApplication( const Wt::WEnvironment &env )
{
  return new InterSpecApp( env );
}// Wt::WApplication *createApplication(const Wt::WEnvironment& env)


#ifdef _WIN32
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
}//void getUtf8Args()
#endif


void processCustomArgs( int argc, char **argv )
{
  bool set_serial_num_file = false;
  for( int i = 1; i < (argc-1); ++i )
  {
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER )
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

#if( USE_SQLITE3_DB )
    if( argv[i] == std::string("--userdb") )
    {
      try
      {
        // TODO: make sure the filename is valid and we can write to it, and also let using = sign
        DataBaseUtils::setPreferenceDatabaseFile( argv[i+1] );
      }catch( std::exception & )
      {
        std::cerr << "Invalid userdb ('" << argv[i+1] << "') specified" << std::endl;
      }
    }//if( argv[i] == std::string("--userdatadir") )
#endif
   
    if( argv[i] == std::string("--static-data-dir") )
    {
      try
      {
        // Let user use the = sign
        InterSpec::setStaticDataDirectory( argv[i+1] );
      }catch( std::exception &e )
      {
        std::cerr << "Invalid static data directory ('" << argv[i+1] << "') specified: " << e.what() << std::endl;
      }
    }//if( argv[i] == std::string("--userdatadir") )

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
