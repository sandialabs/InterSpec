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

#include <mutex>
#include <string>
#include <iostream>

#include "InterSpec/DbToFilesystemLink.h"
#include "SpecUtils/UtilityFunctions.h"

#include <Wt/Dbo/Dbo>
#include <Wt/WString>
#include <Wt/WDateTime>
#include <Wt/WApplication>
#include <Wt/Dbo/WtSqlTraits>
#include <Wt/Dbo/backend/Sqlite3>


using namespace std;
using namespace Wt;

namespace DbToFilesystemLink
{
  std::mutex FileNumToFilePathDBNameMutex;
  std::string FileNumToFilePathDBName = "FileNumToFilePath.sqlite";
  
  
  bool setFileNumToFilePathDBNameBasePath( const std::string &path )
  {
    if( UtilityFunctions::is_directory( path ) )
    {
      string newpath = UtilityFunctions::append_path( path, "FileNumToFilePath.sqlite" );
      
      {
        std::lock_guard<std::mutex> lock( FileNumToFilePathDBNameMutex );
        FileNumToFilePathDBName = newpath;
      }
      //iOS 10.3 changed permissions model, so as a hack lets write this file so we can set permissions.
      if( !UtilityFunctions::is_file(newpath) )
        addFileToOpenToDatabase( "DummyFile" );
    }else
    {
      cerr << "setFileNumToFilePathDBNameBasePath: invalid input path: " << path << endl;
      return false;
    }
    
    return true;
  }//setFileNumToFilePathDBNameBasePath(...)
  
  FileIdToLocation getFileToOpenToInfo( const int file_id )
  {
    string dbfile;
    
    {
      std::lock_guard<std::mutex> lock( FileNumToFilePathDBNameMutex );
      dbfile = FileNumToFilePathDBName;
    }
    
    
    if( !UtilityFunctions::is_file(dbfile) )
      throw runtime_error( "Database file '" + dbfile + "' doesnt exist" );
    
    Wt::Dbo::Session session;
    Wt::Dbo::backend::Sqlite3 fileIdDB( dbfile );
    fileIdDB.setProperty( "show-queries", "false" );
    session.setConnection( fileIdDB );
    session.mapClass<FileIdToLocation>( "FileIdToLocation" );
    try{ session.createTables(); }catch(...){}
    
    FileIdToLocation info;
    Wt::Dbo::ptr<FileIdToLocation> file_entry;
    
    {//begin block to read database
      Dbo::Transaction transaction( session );
      file_entry = session.find<FileIdToLocation>().where( "id = ?" )
      .bind( file_id ).limit(1).resultValue();
      if( file_entry )
      {
        info = *file_entry;
        file_entry.modify()->m_fufilled = true;
      }
      transaction.commit();
    }//end block to read database

    if( !file_entry )
      throw runtime_error( "Unable to find file requested in URL in the database" );
    
    return info;
  }//FileIdToLocation getFileToOpenToInfo( const int index )
  
  int addFileToOpenToDatabase( const std::string &foreground )
  {
    FileIdToLocation info;
    info.m_foregroundFilePath = foreground;
    return addFileToOpenToDatabase( info );
  }//int addFileToOpenToDatabase( const std::string &foreground )
  
  
  int addFileToOpenToDatabase( const FileIdToLocation &info )
  {
    try
    {
      Wt::Dbo::Session session;
      std::string dbfilename;
      
      {
        std::lock_guard<std::mutex> lock( FileNumToFilePathDBNameMutex );
        dbfilename = FileNumToFilePathDBName;
      }
      
      Wt::Dbo::backend::Sqlite3 fileIdDB( dbfilename );
      fileIdDB.setProperty( "show-queries", "false" );
      session.setConnection( fileIdDB );
      session.mapClass<FileIdToLocation>( "FileIdToLocation" );
      try{ session.createTables(); }catch(...){}
      
      Wt::Dbo::Transaction transaction( session );
      Wt::Dbo::ptr<FileIdToLocation> dbfile( new FileIdToLocation( info ) );
      dbfile = session.add( dbfile );
      transaction.commit();
      
      return static_cast<int>( dbfile.id() );
    }catch( std::exception &e )
    {
      cerr << "Error creating FileIdToLocation database: " << e.what() << endl;
    }
    
    return -1;
  }//int addFileToOpenToDatabase( const std::string file );

  

  
  
  
}//namespace DbToFilesystemLink
