#ifndef DbToFilesystemLink_h
#define DbToFilesystemLink_h
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

#include <Wt/Dbo/Dbo>
#include <Wt/WString>
#include <Wt/WDateTime>

//This namespace allows you to open a file from the filesystem by specifting
//  an index in the URL, which cooresponds to an entry in the sqlite3 database
//  (at filesystem location FileNumToFilePathDBName), which points to the
//  path on the filesystem of the file to be opened
//  e.x. https://hekili.ca.sandia.gov/srb/spectroscopy/InterSpec/?specfile=3
//Note that this class (or app) currently doesnt properly treat non-ascii
//  characters in filenames correctly.


namespace DbToFilesystemLink
{
  class FileIdToLocation;
  
  //setFileNumToFilePathDBNameBasePath(...): returns whether successful or not
  //  Default path is CWD.
  bool setFileNumToFilePathDBNameBasePath( const std::string &path );
  
  //addFileToOpenToDatabase(...): returns the index of the passed in filerequest
  int addFileToOpenToDatabase( const FileIdToLocation &filerequest );
  int addFileToOpenToDatabase( const std::string &foreground );
  
  //getFileToOpenToInfo(...): will throw an exception if index is invalid, or
  //  it runs into other trouble.
  //FileIdToLocation::m_fufilled will be marked as true in the database, but
  //  not for the returned object (to allow you to check this)
  FileIdToLocation getFileToOpenToInfo( const int index );
  
  //FileIdToLocation:
  class FileIdToLocation
  {
  public:
    //m_utcRequestTime:
    Wt::WDateTime m_utcRequestTime;
    
    //m_fufilled: in the future may be used to prevent a request from being
    //  fufilled multiple times
    bool m_fufilled;
    
    //m_userName: if specified, will need to match current username when this
    //  request is tried to be fufilled to open the file
    Wt::WString m_userName;
    
    //m_foregroundFilePath: must be a valid path in order for this request to
    //  be fulfilled, or it must be a url beginning with "interspec://".
    //
    // If this path begins with "interspec://", then InterSpecApp::handleAppUrl(...) will
    //  be called with this URL, and attempting to load the path or the background/secondary
    //  path will not be done.
    Wt::WString m_foregroundFilePath;
    
    //m_backgroundFilePath: if left blank or is invalid, wont be loaded
    Wt::WString m_backgroundFilePath;
    
    //m_secondForegroundFilePath: if left blank or is invalid, wont be loaded
    Wt::WString m_secondForegroundFilePath;
    
    
    FileIdToLocation()
    : m_utcRequestTime( Wt::WDateTime::currentDateTime() ),
      m_fufilled( false ),
      m_userName(),
      m_foregroundFilePath(),
      m_backgroundFilePath(),
      m_secondForegroundFilePath()
    {}
    
    template<class Action>
    void persist(Action& a)
    {
      Wt::Dbo::field( a, m_utcRequestTime, "RequestTime" );
      Wt::Dbo::field( a, m_fufilled,       "Fufilled" );
      Wt::Dbo::field( a, m_userName,       "UserName" );
      
      Wt::Dbo::field( a, m_foregroundFilePath, "ForegroundFilePath" );
      Wt::Dbo::field( a, m_backgroundFilePath, "BackgroundFilePath" );
      Wt::Dbo::field( a, m_secondForegroundFilePath, "SecondForegroundFilePath" );
    }//void persist(Action& a)
  };//class FileToTable
}//namespace DbToFilesystemLink

#endif  //DbToFilesystemLink_h
