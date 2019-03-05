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

//#include "InterSpec/GammaInteractionCalc.h"

int main( int argc, char **argv )
{
//  GammaInteractionCalc::example_integration();
//  return 1;
  
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
  
  //TODO 20181009: make the following setting an option on startup for electron/iOS/Android/macOS
//"Have not yet implemented calling SerialToDetectorModel::set_detector_model_input_csv(...) properly"
  SerialToDetectorModel::set_detector_model_input_csv( "data/OUO_detective_serial_to_model.csv" );
  
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


