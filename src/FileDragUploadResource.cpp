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

// Disable some warnings in xutility
#pragma warning(disable:4996)

#include <string>
#include <fstream>
#include <iostream>
#include <iterator>

#include <Wt/Utils>
#include <Wt/WResource>
#include <Wt/Json/Value>
#include <Wt/Json/Parser>
#include <Wt/Json/Object>
#include <Wt/WJavaScript>
#include <Wt/WApplication>
#include <Wt/WPaintDevice>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/WPaintedWidget>
#include <Wt/WContainerWidget>
#include <Wt/WCssDecorationStyle>
#include <Wt/Chart/WAbstractChart>


#include "SpecUtils/Filesystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/FileDragUploadResource.h"

#ifdef _WIN32
#include "SpecUtils/StringAlgo.h"
#endif


using namespace Wt;
using namespace std;


FileDragUploadResource::FileDragUploadResource( WObject *parent  )
    : WResource( parent ),
      m_fileDrop( this )
{
  setUploadProgress( true );

  //If we enable this next call, then we dont get immediate action to
  //  the event loop when we upload a file, as we would expect.
  //  Not sure why.
  //  But this is fine as `handleRequest(...)` explicitly takes a `WApplication::UpdateLock`
  //setTakesUpdateLock( true );
}


FileDragUploadResource::~FileDragUploadResource()
{
  beingDeleted();
  
  clearSpooledFiles();
}


Wt::Signal<std::string,std::string > &FileDragUploadResource::fileDrop() //<display_name,spool_name>
{
  return m_fileDrop;
}


void FileDragUploadResource::handleRequest( const Http::Request& request,
                                            Http::Response& response )
{
  response.setMimeType("text/html; charset=utf-8");
  response.addHeader("Expires", "Sun, 14 Jun 2020 00:00:00 GMT");
  response.addHeader("Cache-Control", "max-age=315360000");

  std::ostream &output = response.out();

  if( request.tooLarge() )
  {
    response.setStatus( 413 );
    const int64_t max_size = wApp->maximumRequestSize();
    output << "Sorry, the maximum upload size is " << max_size/1024 << " kb";
    return;
  }//if( request.tooLarge() )

  auto app = WApplication::instance();
  if( !app )
  {
    cerr << "Uploaded file to a non Wt-Sesssion - bailing" << endl;
    response.setStatus( 403 ); //Forbidden
    return;
  }//if( !app )

  WApplication::UpdateLock lock( app );
  if( !lock )
  {
    response.setStatus( 500 );  //Internal Server error
    output << "App session may be over.";
    return;
  }//if( !lock )
  
  
#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP )
  // See FileUploadFcn in InterSpec.js
  const string isFilePath = request.headerValue("Is-File-Path");
  
  if( (isFilePath == "1") || (isFilePath == "true") || (isFilePath == "yes") )
  {
    assert( request.uploadedFiles().empty() );
    
    try
    {
      // We'll make sure this is the primary electron instance making this request, jic
      if( !InterSpecApp::isPrimaryWindowInstance() )
        throw runtime_error( "Opening via full path only allowed for primary instance" );
      
      // We'll also double check the request is from this computer (which should already be ensured
      //  by checking the primary instance, but jic)
      const string clientAddress = request.clientAddress();  //'127.0.0.1'
      //const string &clientAddress = wApp->environment().clientAddress();
      if( clientAddress.find("127.0.0.1") == std::string::npos )
        throw runtime_error( "Opening via full path only allowed from localhost" );

      std::istreambuf_iterator<char> eos;
      string body(std::istreambuf_iterator<char>(request.in()), eos);
      
      Json::Object result;
      Json::parse( body, result );
      if( !result.contains("fullpath") )
        throw std::runtime_error( "Body JSON did not contain a 'fullpath' entry." );
      
      const WString wfullpath = result.get("fullpath");
      const std::string fullpath = wfullpath.toUTF8();
      
      {// begin test if can read file
#ifdef _WIN32
        const std::wstring wname = SpecUtils::convert_from_utf8_to_utf16( fullpath );
        std::ifstream file( wname.c_str() );
#else
        std::ifstream file( fullpath.c_str() );
#endif
        if( !file.good() )
          throw std::runtime_error( "Could not read '" + fullpath + "'" );
      }// end test if can read file
      
      cout << "Will open spectrum file using path='" << fullpath << "'" << endl;

      m_spooledFiles.push_back( {SpecUtils::filename(fullpath), fullpath, false} );
      m_fileDrop.emit( SpecUtils::filename(fullpath), fullpath );
    }catch( std::exception &e )
    {
      cerr << "Failed to parse fullpath POST request: " << e.what() << " - returning status 406.\n";
      
      response.setStatus( 406 );
    }//try / catch to figure out how to interpret
    
    app->triggerUpdate();
    return;
  }//if( (fullpath == "1") || (fullpath == "true") || (fullpath == "yes") )
#endif  //BUILD_AS_ELECTRON_APP
  
  // TODO: Instead of having multiple Resources we upload files to, we could have a single
  //       upload Resource, and instead use "foreground", "background", "secondary", "multiple", "batch"
  //       as the parameter name, which would then let us upload multiple files at the same time.
  //       The native filesystem path stuff could also be easily modified.
  assert( request.uploadedFiles().size() <= 1 );
  
  for( const Wt::Http::UploadedFileMap::value_type &parname_file : request.uploadedFiles() )
  {
    assert( parname_file.first == "file" );
    
    const Wt::Http::UploadedFile &file = parname_file.second;
    assert( !file.clientFileName().empty() );
    
    if( !file.clientFileName().empty() )
    {
      const string &clientFileName = file.clientFileName();
      const string &spoolFileName = file.spoolFileName();
      file.stealSpoolFile(); //Take ownership of this file
      m_spooledFiles.emplace_back( clientFileName, spoolFileName, true );
      m_fileDrop.emit( clientFileName, spoolFileName );
    }
  }//for( const Wt::Http::UploadedFileMap::value_type &parname_file : request.uploadedFiles() )

  app->triggerUpdate();
}//void FileDragUploadResource::handleRequest(...)


vector<tuple<string,string,bool>> FileDragUploadResource::takeSpooledFiles()
{
  vector<tuple<string,string,bool>> files;
  files.swap( m_spooledFiles );
  return files;
}

const std::vector<std::tuple<std::string,std::string,bool>> &FileDragUploadResource::spooledFiles() const
{
  return m_spooledFiles;
}


void FileDragUploadResource::clearSpooledFiles()
{
  for( size_t i = 0; i < m_spooledFiles.size(); ++i )
  {
    const bool isTempFile = std::get<2>( m_spooledFiles[i] );

    if( isTempFile )
    {
      const string &spooledName = std::get<1>( m_spooledFiles[i] );
      const bool success = SpecUtils::remove_file( spooledName );
      if( !success )
        cerr << "FileDragUploadResource: Warning, could not delete file '" << spooledName << "'" << endl;
    }
  }//for( size_t i = 0; i < m_spooledFiles.size(); ++i )

  m_spooledFiles.clear();
}//void clearSpooledFiles()
