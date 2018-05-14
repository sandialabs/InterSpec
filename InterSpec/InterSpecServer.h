//
//  InterSpecServer.h
//  InterSpec
//
//  Created by Johnson, William C on 6/25/17.
//  Copyright © 2017 Sandia National Laboratories. All rights reserved.
//

#ifndef InterSpecServer_hpp
#define InterSpecServer_hpp

#include <string>
#include <stdio.h>

#include "InterSpec_config.h"

#include <Wt/WApplication>

namespace InterSpecServer
{
  void startServer( int argc, char *argv[],
                                Wt::WApplication::ApplicationCreator createApplication );
  void killServer();
  
  
  
  //portBeingServedOn(): will only be valid if this instance of the app is
  //  serving the webpages.  Will be -1 if not serving.
  int portBeingServedOn();
  
  //urlBeingServedOn(): will only be valid if this instance of the app is
  //  serving the webpages.  Will be empty if not serving.
  //  Example value returned: "http://127.0.0.1:7234"
  std::string urlBeingServedOn();
  
  
  bool changeToBaseDir( int argc, char *argv[] );
  
  std::string getWtConfigXml( int argc, char *argv[] );
}


#endif /* InterSpecServer_hpp */
