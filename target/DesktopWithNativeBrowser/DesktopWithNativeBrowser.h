#ifndef DesktopWithNativeBrowser_h
#define DesktopWithNativeBrowser_h
/* InterSpec: an application to analyze spectral gamma radiation data.
 
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

#include "InterSpec_config.h"

#include <map>
#include <string>
#include <stdlib.h>

#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>

#include <Wt/Dbo/Dbo>
#include <Wt/WString>
#include <Wt/WDateTime>
#include <Wt/WResource>
#include <Wt/WApplication>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
            
class InterSpecApp;


namespace DesktopWithNativeBrowser
{
  void experimental_startServer( int argc, char *argv[],
                                 Wt::WApplication::ApplicationCreator createApplication );
  void experimental_killServer();
	
  //runInDefaultBrowser(...): creates a server on a random local port.  Does not
  //  return untill the server is termenated
  //  argv[0]   (required): is the name of the application, e.g. "/Users/wcjohns/radana/InterSpec/build/bin/InterSpec.exe"
  //  argv[1]   (optional): if a file that exists on filesystem, this file is
  //                        opened.  No effect on iOS and Android systems
  //  argv[>=1]: 
  //             --basedir="/path/to/resources/and/data" will cause app to 
  //               change the current working directory to one specified, where
  //               its assumed Wt's "resources" folder is, as well as 
  //               InterSpec_resources, data, sandia.decay.xml, sandia.xray.xml,
  //               etc.
  //             --forceserve: if specified, will create a new server no matter 
  //               what
  //             --nobrowsertab: if specified no tab in the native browser will
  //               be opened, unless argv[1] specifies a file, then it will
  //               be openeind in a new tab no matter what
  int runInDefaultBrowser( int argc, char *argv[],
                           Wt::WApplication::ApplicationCreator createApplication );
        
  //portBeingServedOn(): will only be valid if this instance of the app is 
  //  serving the webpages.  Will be -1 if not serving.
  int portBeingServedOn();
  
  //urlBeingServedOn(): will only be valid if this instance of the app is 
  //  serving the webpages.  Will be empty if not serving.
  //  Example value returned: "http://127.0.0.1:7234"
  std::string urlBeingServedOn();
  
  
  //runningInstances(): returns pointers to the apps this proccess created, make
  //  sure to get a WApplication::UpdateLock before using the returned pointers
//  std::vector<const InterSpecApp *> runningInstances();
  
  int DesktopRunInDefaultBrowser( int argc, char *argv[],
                                 Wt::WApplication::ApplicationCreator creator );
	
  enum InstanceState
  {
    InstanceStateNotKnown,
    OtherInstanceRunning,
    NoOtherInstanceRunning,
  };//enum InstanceState
  
  bool changeToBaseDir( int argc, char *argv[] );
  bool shouldForceServe( int argc, char *argv[] );
  bool shouldOpenBlankTab( int argc, char *argv[] );
  std::string getWtConfigXml( int argc, char *argv[] );
  boost::filesystem::path fileToOpen( int argc, char *argv[] );
  boost::filesystem::path instanceInfoFile( const char *exename );
  
  //addFileToOpenToDatabase(...): returns -1 on error
  int addFileToOpenToDatabase( boost::filesystem::path file );
  
#if( !defined(IOS) && !defined(ANDROID) )
  void openUrl( std::string url, boost::filesystem::path file = "" );
#endif
  
  void handleHttpCallback( boost::system::error_code err,
                           const Wt::Http::Message &response,
                           Wt::Http::Client *client,
                           int argc, char *argv[],
                           std::string url,
                           Wt::WApplication::ApplicationCreator appCreator );                          
  
  bool checkForExistingInstanceAndCallIt( int argc, char *argv[],
                                 Wt::WApplication::ApplicationCreator appCreator  );
   
  class ServerAliveResource : public Wt::WResource
  {
  public:
    ServerAliveResource( Wt::WObject *parent = 0 )
     : WResource( parent ){}
    virtual ~ServerAliveResource(){ beingDeleted(); }
    virtual void handleRequest( const Wt::Http::Request &request,
                                Wt::Http::Response &response )
    {
      //Could use WTestEnvironment here to make sure an instance can actually be
      //  made
      response.addHeader( "Cache-Control", "no-cache, no-store, must-revalidate" );
      response.addHeader( "Pragma", "no-cache" );
      response.addHeader( "Expires", "0" );
      response.out() << "Alive\r\n";
      response.setStatus( 200 );
    }
  };//class DetectorUploadResource 
}//namespace DesktopWithNativeBrowser

#endif //DesktopWithNativeBrowser_h
