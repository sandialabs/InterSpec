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
#include <iterator>

#include <Wt/Dbo/Dbo>
#include <Wt/Dbo/SqlConnection>

#if( USE_MYSQL_DB )
#include <Wt/Dbo/backend/MySQL>
#endif

#if( USE_SQLITE3_DB )
#include <Wt/Dbo/backend/Sqlite3>
#endif

#include "InterSpec/InterSpecUser.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/DataBaseVersionUpgrade.h"


namespace DataBaseVersionUpgrade
{
  void checkAndUpgradeVersion()
  {
    //Synchronize database here
    //Set up the session; open the database.
    std::unique_ptr<Wt::Dbo::SqlConnection> database;
    
    //Synchronizes database schema is up to date...if not, does table modifications
    //to keep database in sync with what Interspec requires.
    //Incremental updates are done to upgrade database to latest version
    //Make sure query works both in sqlite and mysql
    //firefox/sqlite UI editor: https://code.google.com/p/sqlite-manager/
    //sqllite/mysql simulator: http://sqlfiddle.com/#!7/c6d5b
    //Mysql check: http://www.piliapp.com/mysql-syntax-check/
    
    
    int version = 0;
    
    {
      try
      {
        try
        {
          std::shared_ptr<Wt::Dbo::Session> sqlSession = getSession( database );
          mapDbClasses( sqlSession.get() );
          sqlSession->createTables();
          
          //If the previous statment didnt throw, then this is a brand new
          //  database, so lets set the version to current
          setDBVersion( DB_SCHEMA_VERSION, sqlSession );
        }catch( std::exception &e )
        {
          std::cout << "Trying initial table mapping, got exception: " << e.what() << std::endl;
        }
        
        std::shared_ptr<Wt::Dbo::Session> sqlSession = getSession( database );
        Wt::Dbo::Transaction transaction( *sqlSession );
        version = sqlSession->query<int>("select schema_version from InterSpecRegistry");
        transaction.commit();
      }catch( std::exception & )
      {
        //Does not exist, start from scratch version 0
        version = 0;
      } //catch
    } //get version
    
    std::cerr << "INTERSPEC DB VERSION: " << version << std::endl;
    
    if( version == DB_SCHEMA_VERSION )
    {
      std::cerr<<"No need to update database, everything in sync"<<std::endl;
      return;
    } //no need to update database

    
    if( version<4 && version<DB_SCHEMA_VERSION )
    {
      //20190717: removed support for upgrading versions less than 4 to remove
      //          requirement of needing sqlite3.h header available.
#if( USE_SQLITE3_DB )
      //ToDo: rename existing database file, and create a new database
      std::cerr << "Upgrading sqlite3 databases below version 4 no longer supported - things will likely fail" << std::endl;
#elif( USE_MYSQL_DB )
      std::cerr << "Upgrading MySQL databases below version 4 no longer supported - things will likely fail" << std::endl;
#endif
    }
    
    /* I *believe* there to be NO InterSpec instalations with a version less
       than 5, so we could safely delete everything above here, but I guess I'll
       leave in as template code for future upgrades.
     */
    if( version<5 && version<DB_SCHEMA_VERSION )
    {
      try
      {
        std::shared_ptr<Wt::Dbo::Session> sqlSession = getSession( database );

        executeSQL("ALTER TABLE UserState ADD COLUMN SnapshotTagParent_id bigint",sqlSession);
            
        version = 5;
        setDBVersion( version, sqlSession );
        std::cout << "DB_SCHEMA_VERSION has been upgraded to " << version << std::endl;
      }catch( std::exception &e )
      {
        const std::string errormsg = e.what();
        std::cerr << "DB_SCHEMA_VERSION 5 : Failed to add SnapshotTagParent to table UserState."
                  << errormsg << std::endl;
        
        if( errormsg.find("duplicate") != std::string::npos )
        {
          try
          {
            std::shared_ptr<Wt::Dbo::Session> sqlSession = getSession( database );
            version = 5;
            setDBVersion( version, sqlSession );
            std::cout << "DB_SCHEMA_VERSION has been set to " << version << std::endl;
          }catch( std::exception &ee )
          {
            std::cerr << "DB_SCHEMA_VERSION 5 : Failed to force version to 5: "
                      << ee.what() << std::endl;
          }
        }//if( the error message indicated column was already in schema )
      }//try / catch( add SnapshotTagParent_id column )
    }//if (version<5 && version<DB_SCHEMA_VERSION)
    
    
    /* Version 6 created 20181018.  Adds table "InterSpecGlobalSetting" to hold
       options that should be globally applied on application startup, like
       the mapping from serial number to detector version, or update to
       sandia.decay.xml.
     */
    if( version<6 && version<DB_SCHEMA_VERSION )
    {
      //For a while we had a table "InterSpecGlobalSetting" but this got removed 20191015, so the upgrade code also got removed.
      std::shared_ptr<Wt::Dbo::Session> sqlSession = getSession( database );
      version = 7;  //This should have been 6, so 
      setDBVersion( version, sqlSession );
    }//if( version<6 && version<DB_SCHEMA_VERSION )
    
    //Skipping version 7, due to accidentally marging version 6 as 7 into DB
    //  (although version 6 was never distributed so it would have been fine,
    //   but just to be thorough)
    
    if( version<8 && version<DB_SCHEMA_VERSION )
    {
      std::shared_ptr<Wt::Dbo::Session> sqlSession = getSession( database );
        
      executeSQL("DROP TABLE IF EXISTS ColorThemeInfo;",sqlSession);
        
      //The following has only been checked for SQLite3
      const char *sql_statement = R"sql(create table "ColorThemeInfo" (
      "id" integer primary key autoincrement,
      "version" integer not null,
      "InterSpecUser_id" bigint,
      "ThemeName" text not null,
      "ThemeDescription" text not null,
      "CreationTime" text,
      "ModifiedTime" text,
      "JsonData" text not null,
      constraint "fk_ColorThemeInfo_InterSpecUser" foreign key ("InterSpecUser_id") references "InterSpecUser" ("id") on delete cascade deferrable initially deferred
      ))sql";
        
      executeSQL( sql_statement, sqlSession );
      
      version = 8;
      setDBVersion( version, sqlSession );
    }//if( version<6 && version<DB_SCHEMA_VERSION )
    
    
    if( version<9 && version<DB_SCHEMA_VERSION )
    {
      std::shared_ptr<Wt::Dbo::Session> sqlSession = getSession( database );
      
      //The following has only been checked for SQLite3.  Positions 'ColorThemeJson'
      //  after the previous last column, witch I think Dbo needs.
      const char *sql_statement = "ALTER TABLE UserState ADD COLUMN ColorThemeJson text;";
      executeSQL( sql_statement, sqlSession );
      
      version = 9;
      setDBVersion( version, sqlSession );
    }//if( version<6 && version<DB_SCHEMA_VERSION )
    
    
    if( version<10 && version<DB_SCHEMA_VERSION )
    {
      /*
      try
      {
        std::shared_ptr<Wt::Dbo::Session> sqlSession = getSession( database );
        sqlSession->mapClass<UseDrfPref>( "UseDrfPref" );
        sqlSession->createTables();
        std::cout << "Created the UseDrfPref Table" << std::endl;
      }catch( std::exception &e )
      {
        //Get error Class 13InterSpecUser was not mapped
        std::cerr << "DB_SCHEMA_VERSION 10: Failed to create UseDrfPref table: "
        << e.what() << std::endl
        << "I guess we'll go on, but things will probably blow up." << std::endl;
      }
      */
      
      std::shared_ptr<Wt::Dbo::Session> sqlSession = getSession( database );
      
      const char *sql_statement = nullptr;
      
      sql_statement = R"Delim(create table "UseDrfPref" (
      "id" integer primary key autoincrement,
      "version" integer not null,
      "InterSpecUser_id" bigint,
      "MatchField" integer not null,
      "Flags" integer not null,
      "DrfIndex" bigint not null,
      "Criteria" text not null,
      constraint "fk_UseDrfPref_InterSpecUser" foreign key ("InterSpecUser_id") references "InterSpecUser" ("id") on delete cascade deferrable initially deferred
      ))Delim";
      executeSQL( sql_statement, sqlSession );
      

      sql_statement = "ALTER TABLE DetectorPeakResponse ADD COLUMN m_flags bigint default 0 not null;";
      executeSQL( sql_statement, sqlSession );
      sql_statement = "ALTER TABLE DetectorPeakResponse ADD COLUMN m_expOfLogPowerSeriesUncerts text default \"\" not null;";
      executeSQL( sql_statement, sqlSession );
      sql_statement = "ALTER TABLE DetectorPeakResponse ADD COLUMN m_resolutionUncerts text default \"\" not null;";
      executeSQL( sql_statement, sqlSession );
      sql_statement = "ALTER TABLE DetectorPeakResponse ADD COLUMN m_lowerEnergy double precision default 0 not null;";
      executeSQL( sql_statement, sqlSession );
      sql_statement = "ALTER TABLE DetectorPeakResponse ADD COLUMN m_upperEnergy double precision default 0 not null;";
      executeSQL( sql_statement, sqlSession );
      sql_statement = "ALTER TABLE DetectorPeakResponse ADD COLUMN m_createdUtc bigint default 0 not null;";
      executeSQL( sql_statement, sqlSession );
      sql_statement = "ALTER TABLE DetectorPeakResponse ADD COLUMN m_lastUsedUtc bigint default 0 not null;";
      executeSQL( sql_statement, sqlSession );
      
      version = 10;
      setDBVersion( version, sqlSession );
    }//if( version<6 && version<DB_SCHEMA_VERSION )
    
    
    if( version<11 && version<DB_SCHEMA_VERSION )
    {
      std::shared_ptr<Wt::Dbo::Session> sqlSession = getSession( database );
      
      //The following has only been checked for SQLite3.  Positions 'ColorThemeJson'
      //  after the previous last column, witch I think Dbo needs.
      const char *sql_statement = "ALTER TABLE UserState ADD COLUMN GammaXsToolUri text;";
      executeSQL( sql_statement, sqlSession );
      sql_statement = "ALTER TABLE UserState ADD COLUMN DoseCalcToolUri text;";
      executeSQL( sql_statement, sqlSession );
      sql_statement = "ALTER TABLE UserState ADD COLUMN OneOverR2ToolUri text;";
      executeSQL( sql_statement, sqlSession );
      sql_statement = "ALTER TABLE UserState ADD COLUMN UnitsConverterToolUri text;";
      executeSQL( sql_statement, sqlSession );
      sql_statement = "ALTER TABLE UserState ADD COLUMN NucDecayInfoUri text;";
      executeSQL( sql_statement, sqlSession );
      sql_statement = "ALTER TABLE UserState ADD COLUMN EnergyRangeSumUri text;";
      executeSQL( sql_statement, sqlSession );
      sql_statement = "ALTER TABLE UserState ADD COLUMN FluxToolUri text;";
      executeSQL( sql_statement, sqlSession );
      
      version = 11;
      setDBVersion( version, sqlSession );
    }//if( version<11 && version<DB_SCHEMA_VERSION )
    
    /// ******************************************************************
    /// DB_SCHEMA_VERSION is at 11.  Add Version 12 here.  Update InterSpecUser.h!
    /// ******************************************************************
  }//void checkAndUpgradeVersion()
  
  
  std::shared_ptr<Wt::Dbo::Session> getSession( std::unique_ptr<Wt::Dbo::SqlConnection> &database )
  {
    if( !database )
      database.reset( DataBaseUtils::getDatabaseConnection() );
    database->setProperty( "show-queries", "true" );
    
    std::shared_ptr<Wt::Dbo::Session> sqlSession;
    sqlSession.reset( new Wt::Dbo::Session() );
    sqlSession->setConnection( *database );
    
    return sqlSession;
  }//shared_ptr<Session> getSession( scoped_ptr<SqlConnection> &database )

  
  
  //Convenience function to execute SQL commands
  void executeSQL( const std::string &sql, std::shared_ptr<Wt::Dbo::Session> sqlSession)
  {
    Wt::Dbo::Transaction transaction( *sqlSession );
    sqlSession->execute(sql);
    transaction.commit();
  } //executeSQL(std::string sql, std::shared_ptr<Wt::Dbo::Session> m_sqlSession)
  
  //Sets the DB registry with the schema version
  void setDBVersion(int version, std::shared_ptr<Wt::Dbo::Session> sqlSession)
  {
    //clear
    std::cerr<<"Clear registry table"<<std::endl;
    try
    {
      executeSQL("DELETE FROM InterSpecRegistry",sqlSession);
    }catch( std::exception & )
    {
    }
    
    {//add version to registry
      Wt::Dbo::Transaction transaction( *sqlSession );
      std::cerr<<"Initializing registry to set DB version: "<<version<<std::endl;
      executeSQL("CREATE TABLE IF NOT EXISTS InterSpecRegistry(schema_version INTEGER);", sqlSession);
      sqlSession->execute("INSERT INTO InterSpecRegistry(schema_version) VALUES (?)").bind(version);
      transaction.commit();
    }//add version to registry
  } //setDBVersion(int version, std::shared_ptr<Wt::Dbo::Session> m_sqlSession)

  
}//namespace DataBaseVersionUpgrade


using namespace Wt;
using namespace std;
