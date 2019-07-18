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
#include "InterSpec/InterSpecApp.h"
#include "SpecUtils/UtilityFunctions.h"
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

#if( BUILD_AS_ELECTRON_APP )
#include "target/electron/ElectronUtils.h"
#endif

//Forward declaration
Wt::WApplication *createApplication( const Wt::WEnvironment &env );

#ifdef _WIN32
void getUtf8Args( int &argc, char ** &argv );
#endif

void processCustomArgs( int argc, char **argv );


//#include "InterSpec/GammaInteractionCalc.h"

int main( int argc, char **argv )
{
//  GammaInteractionCalc::example_integration();
//  return 1;
  
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
            << WT_VERSION << ", from executable compiled on "
            << std::dec << COMPILE_DATE_AS_INT << std::endl;
  
#if( PERFORM_DEVELOPER_CHECKS )
  std::cout << "Developer tests are being performed" << std::endl;
#endif

  std::cout << std::endl;
  
#if(WT_VERSION>=0x3030300)
  //Make it so WString defaults to assuming std::string or char * are UTF8
  //  encoded, rather than the system encoding.
  Wt::WString::setDefaultEncoding( Wt::UTF8 );
#endif
  
  
  processCustomArgs( argc, argv );

  
#if( BUILD_AS_ELECTRON_APP )
  return ElectronUtils::run_app(argc,argv);
#endif
  
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

/** Get command line arguments encoded as UTF-8.
    This function just leaks the memory
 */
void getUtf8Args( int &argc, char ** &argv );
{
  //ToDo:
  //  - on windows:
  //     - Convert command line arguments to UTF-8
  //     - Make all UtilityFunctions filesystem calls convert inputs from UTF8 to UTF16 and call wide versions of windows functions
  //     - Make all UtilityFunctions filesystem related functiosn return UTF8 encoded results
  //     - Make sure all instacces of DataBaseUtils::get/setPreferenceDatabaseFile() are aware of changes
  //     - Check that database upgrade functions use wide functions
  //     - Check target/electron/ElectronUtils.cpp for consistency
  //     - Check file query widget to make sure it is okay 
  //     - Check all instances of InterSpec::setWritableDataDirectory()/InterSpec::writableDataDirectory() are aware of changes
  //     - Check anywhere that calls cwd() or get_working_path() or temp_dir()
  //     - TO convert between utf16 and utf8 see from http://www.nubaria.com/en/blog/?p=289, or can do it cross-platform on c++11
  
  //
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
  
    const std::string asutf8 = UtilityFunctions::convert_from_utf16_to_utf8( argvw[i] );
    argv[i] = (char *)malloc( sizeof(char)*(asutf8.size()+1) );
    strcpy( argv[i], asutf8.c_str() );
  }//for( int i = 0; i < argc; ++i)

  // Free memory allocated for CommandLineToArgvW arguments.
  LocalFree(argvw);
}//void processCustomArgs()
#else
void processCustomArgs( int argc, char **argv )
{
  for( int i = 1; i < (argc-1); ++i )
  {
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || (BUILD_AS_LOCAL_SERVER && (defined(WIN32) || defined(__APPLE__)) ) )
    if( argv[i] == std::string("--userdatadir") )
    {
      try
      {
        InterSpec::setWritableDataDirectory( argv[i+1] );
        const std::vector<std::string> serial_db = UtilityFunctions::ls_files_in_directory( argv[i+1], "serial_to_model.csv" );
        if( !serial_db.empty() )
          SerialToDetectorModel::set_detector_model_input_csv( serial_db[0] );
      }catch( std::exception &e )
      {
        std::cerr << "Invalid userdatadir ('" << argv[i+1] << "') specified"
                  << std::endl;
      }
    }//if( argv[i] == std::string("--userdatadir") )
#endif  //if( not a webapp )
  }//for( int i = 1; i < (argc-1); ++i )
}//void processCustomArgs( int argc, char **argv )
#endif
