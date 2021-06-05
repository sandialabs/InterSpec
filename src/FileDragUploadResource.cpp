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


#include <Wt/WResource>
#include <Wt/WJavaScript>
#include <Wt/WApplication>
#include <Wt/WPaintDevice>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/WPaintedWidget>
#include <Wt/WContainerWidget>
#include <Wt/WCssDecorationStyle>
#include <Wt/Chart/WAbstractChart>

#include "InterSpec/InterSpecApp.h"
#include "SpecUtils/Filesystem.h"
#if( !ANDROID && !IOS )
#include "InterSpec/FileDragUploadResource.h"
#endif

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
      m_fileDrop( NULL )
{
  m_fileDrop = new Wt::Signal<std::string, std::string>( this );
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

  if( m_fileDrop )
    delete m_fileDrop;
}


Wt::Signal<std::string,std::string > *FileDragUploadResource::fileDrop() //<display_name,spool_name>
{
  return m_fileDrop;
}


void FileDragUploadResource::handleRequest( const Http::Request& request,
                                            Http::Response& response )
{
//XXX - It's assuming handleRequest(...) is only called once all the data
//      is uploaded.
//XXX - Need to deal with case if the data was spooled to file ? - maybe not?
//    std::vector<Http::UploadedFile> files;
//    std::pair<Http::UploadedFileMap::const_iterator, Http::UploadedFileMap::const_iterator> range
//               = request.uploadedFiles().equal_range(key);
//    for(Http::UploadedFileMap::const_iterator i = range.first; i != range.second; ++i)
//      files.push_back(i->second);
//    for (unsigned i = 0; i < files.size(); ++i)
//      if (!files[i].clientFileName().empty())
//      {
//        m_spooledFiles.push_back( files[i].spoolFileName() );
//        m_dragCanvas->fileDrop().emit( files[i].clientFileName(), files[i].spoolFileName() );
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
      const string userName = request.headerValue( "X-File-Name" );
      
      auto app = WApplication::instance();
      WApplication::UpdateLock lock( app );
      
      if( lock )
      {
        m_spooledFiles.push_back( temp_name );
        if( m_fileDrop )
          m_fileDrop->emit( userName, temp_name );
      
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
