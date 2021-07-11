#ifndef FileDragUploadResource_h
#define FileDragUploadResource_h
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

#include <Wt/WEvent>
#include <Wt/WResource>
#include <Wt/Http/Request>
#include <Wt/Http/Response>


class FileDragUploadResource : public Wt::WResource
{
  //TODO - create own .h, .cpp, and .js  files for FileDragUploadResource
  //
  //FileDragUploadResource is designed to provide a URL to post a file to from
  //  the client side inorder to open the file.  It is not prticularly well
  //  tested as of yet, but seems to work; will not work on older IE..
  //Uploaded file are kept on disk until this object is destructed - which
  //  could be a source of optimization

public:
  FileDragUploadResource( Wt::WObject *parent  );
  virtual ~FileDragUploadResource();
  Wt::Signal<std::string,std::string > &fileDrop(); //<display_name,spool_name>

protected:
  //XXX - It's assuming handleRequest(...) is only called once all the data
  //      is uploaded.  See implementation for other possible issues.
  virtual void handleRequest( const Wt::Http::Request& request,
                              Wt::Http::Response& response );

private:
  std::vector<std::string> m_spooledFiles;
  Wt::Signal< std::string, std::string > m_fileDrop;  //<display_name,spool_name>
};//class FileDragUploadResource ...


#endif  //FileDragUploadResource_h

