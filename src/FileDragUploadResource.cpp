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

/*
 * See also: http://www.webtoolkit.eu/wt/blog/2010/03/02/javascript_that_is_c__
 */
#define INLINE_JAVASCRIPT(...) #__VA_ARGS__


FileDragUploadResource::FileDragUploadResource( WObject *parent  )
    : WResource( parent ),
      m_fileDrop( this )
{
}


FileDragUploadResource::~FileDragUploadResource()
{
  beingDeleted();
  for( size_t i = 0; i < m_spooledFiles.size(); ++i )
  {
    const bool success = SpecUtils::remove_file( m_spooledFiles[i] );
    if( !success )
      cerr << "Warning, could not delete file '" << m_spooledFiles[i] << "'" << endl;
  }
}


Wt::Signal<std::string,std::string > &FileDragUploadResource::fileDrop() //<display_name,spool_name>
{
  return m_fileDrop;
}


void FileDragUploadResource::handleRequest( const Http::Request& request,
                                            Http::Response& response )
{
// TODO: It's assuming handleRequest(...) is only called once all the data is uploaded; this needs
//       to be verified.
  
// TODO: Currently the client side JS FileUploadFcn(...) function just puts file contents inside the
//  the POST body, so it doesnt look like files, for all the handling we do here.
//  We could, and maybe should, use a <form /> to upload the files to take advantage of the
//  WResource spooling of files, and uploading multiple files and all that, like could be done in
//  this next bit of commented out text.
//
//    std::vector<Http::UploadedFile> files;
//    std::pair<Http::UploadedFileMap::const_iterator, Http::UploadedFileMap::const_iterator> range
//               = request.uploadedFiles().equal_range(key);
//    for(Http::UploadedFileMap::const_iterator i = range.first; i != range.second; ++i)
//      files.push_back(i->second);
//    for (unsigned i = 0; i < files.size(); ++i)
//      if (!files[i].clientFileName().empty())
//      {
//        m_spooledFiles.push_back( files[i].spoolFileName() );
//        m_fileDrop.emit( files[i].clientFileName(), files[i].spoolFileName() );
//      }
  
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
  
  
#if( BUILD_AS_ELECTRON_APP )
  // See FileUploadFcn in InterSpec.js
  const string isFilePath = request.headerValue("Is-File-Path");
  
  if( (isFilePath == "1") || (isFilePath == "true") || (isFilePath == "yes") )
  {
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
      
      WApplication::UpdateLock lock( app );

      m_fileDrop.emit( SpecUtils::filename(fullpath), fullpath );
      
      app->triggerUpdate();
    }catch( std::exception &e )
    {
      cerr << "Failed to parse fullpath POST request: " << e.what() << " - returning status 406.\n";
      
      response.setStatus( 406 );
      return;
    }//try / catch to figure out how to interpret
    
    return;
  }//if( (fullpath == "1") || (fullpath == "true") || (fullpath == "yes") )
#endif  //BUILD_AS_ELECTRON_APP
  
  const int datalen = request.contentLength();

  if( datalen )
  {
    const string temp_name = SpecUtils::temp_file_name( wApp->sessionId(), InterSpecApp::tempDirectory() );
#ifdef _WIN32
    const wstring wtemp_name = SpecUtils::convert_from_utf8_to_utf16(temp_name);
    ofstream spool_file( wtemp_name.c_str(), ios::binary|ios::out );
#else
    ofstream spool_file( temp_name.c_str(), ios::binary|ios::out );
#endif
    
    if( spool_file.is_open() )
    {
      spool_file << request.in().rdbuf();
      spool_file.close();
      const string userNameEncoded = request.headerValue( "X-File-Name" );
      const string userName = Wt::Utils::urlDecode(userNameEncoded);
      
      //cerr << "\n\n\nuserName = '" << userName << "'\n\n" << endl;
      
      auto app = WApplication::instance();
      WApplication::UpdateLock lock( app );
      
      if( lock )
      {
        m_spooledFiles.push_back( temp_name );
        m_fileDrop.emit( userName, temp_name );
        app->triggerUpdate();
      }else
      {
        response.setStatus( 500 );  //Internal Server error
        output << "App session may be over.";
      }
    }else
    {
      response.setStatus( 500 );  //Internal Server error
      output << "Coulnt open temporary file on server";
    }//if( spool_file.is_open() ) / else
  }//if( request.contentLength() )
}//void FileDragUploadResource::handleRequest(...)
