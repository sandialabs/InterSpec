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

#include <tuple>
#include <string>

#include <Wt/WEvent>
#include <Wt/WResource>
#include <Wt/Http/Request>
#include <Wt/Http/Response>


class FileDragUploadResource : public Wt::WResource
{
  //FileDragUploadResource is intended to provide a URL to post a file to from
  //  the client side in order to open the file.
  //Uploaded file are kept on disk until this object is destructed - which
  //  could be a source of optimization

public:
  FileDragUploadResource( Wt::WObject *parent  );
  virtual ~FileDragUploadResource();
  Wt::Signal<std::string,std::string> &fileDrop(); //<display_name,spool_name>

  /** Take ownership of all files currently owned by this object.
   *  The first string is the display name, the second is the spooled file name,
   *  and the boolean is true than the caller must take care of deleting the file  
   * (e.g., it is a temporary file uploaded over http that needs deleting).  
   *  If false, the file is must not be deleted (e.g., on a native app and we 
   * revieved the full filesystem path, instead of doing an upload over http).
  */
  std::vector<std::tuple<std::string,std::string,bool>> takeSpooledFiles();

  /** Get currently spooled files.
   *  The first string is the display name, the second is the spooled file name,
   *  and the boolean is true if the file is owned by this object (e.g., it is a temporary file
   *  that will be deleted).  If false, the file is must not be deleted (e.g., on a native
   *  app and we revieved the full filesystem path, instead of doing an upload
   *  over http).
   * Does not transfer ownership of the files.
  */
  const std::vector<std::tuple<std::string,std::string,bool>> &spooledFiles() const;

  /** Clear all files currently owned by this object; e.g., deletes from disk. */
  void clearSpooledFiles();

protected:
  virtual void handleRequest( const Wt::Http::Request& request,
                              Wt::Http::Response& response );

private:
  /** Files currently owned by this object.  
   * The first string is the display name, the second is the spooled file name. 
  */
  std::vector<std::tuple<std::string,std::string,bool>> m_spooledFiles;

  /** Signal emitted when a file is dropped.
   *  The first string is the display name, the second is the spooled file name. 
  */
  Wt::Signal< std::string, std::string > m_fileDrop;  //<display_name,spool_name>
};//class FileDragUploadResource ...


#endif  //FileDragUploadResource_h

