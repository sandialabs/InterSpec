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
#include <cassert>
#include <iostream>
#include <iterator>

#include <Wt/Dbo/Dbo.h>
#include <Wt/Dbo/SqlConnection.h>

#if( USE_MYSQL_DB )
#include <Wt/Dbo/backend/MySQL.h>
#endif

#if( USE_SQLITE3_DB )
#include <Wt/Dbo/backend/Sqlite3.h>
#endif

#include "InterSpec/InterSpecUser.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/DataBaseVersionUpgrade.h"


namespace DataBaseVersionUpgrade
{
  /** Returns whether `table` currently has a column named `column`.

   Asks the database catalog (`pragma_table_info` for SQLite3, `information_schema.columns`
   for MySQL).  Returns false if `table` itself doesnt exist, and logs-and-returns false if
   the query fails - so callers must tolerate the subsequent `ALTER TABLE` failing, making a
   false negative a logged no-op rather than an aborted start-up.
   */
  bool tableHasColumn( const std::string &table, const std::string &column,
                       std::shared_ptr<Wt::Dbo::Session> sqlSession );

  /** Adds the `DetectorPeakResponse.m_drfExtra` column, if it is not already present.

   Used by both the version 13 and version 14 upgrade steps, so that neither aborts start-up
   on a database that already has the column.
   */
  void addDrfExtraColumn( std::shared_ptr<Wt::Dbo::Session> sqlSession );


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

    if( version > DB_SCHEMA_VERSION )
    {
      // All the upgrade steps below are guarded on `version < N`, so a database written by a newer
      //  build silently falls through them all, and then Dbo throws on the first query against a
      //  table whose columns dont match.  We cant downgrade, so just make the reason obvious.
      std::cerr << "Database is version " << version << ", but this build of InterSpec expects "
                << DB_SCHEMA_VERSION << " - it was written by a newer version of InterSpec, and"
                   " can not be downgraded.  Expect database errors." << std::endl;
      return;
    }//if( version > DB_SCHEMA_VERSION )

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
    
    if( version<12 && version<DB_SCHEMA_VERSION )
    {
      std::shared_ptr<Wt::Dbo::Session> sqlSession = getSession( database );
      
      //The following has only been checked for SQLite3.  Positions 'ColorThemeJson'
      //  after the previous last column, witch I think Dbo needs.
      const char *sql_statement = "ALTER TABLE UserState ADD COLUMN DetectionSensitivityToolUri text;";
      executeSQL( sql_statement, sqlSession );
      sql_statement = "ALTER TABLE UserState ADD COLUMN SimpleMdaUri text;";
      executeSQL( sql_statement, sqlSession );
      
      version = 12;
      setDBVersion( version, sqlSession );
    }//if( version<11 && version<DB_SCHEMA_VERSION )

    if( version<13 && version<DB_SCHEMA_VERSION )
    {
      std::shared_ptr<Wt::Dbo::Session> sqlSession = getSession( database );

      // Holds a `<DrfExtra>` XML fragment with optional DRF extensions
      //  (efficiency uncertainty, total efficiency, and future additions).
      //  The following has only been checked for SQLite3.
      addDrfExtraColumn( sqlSession );

      version = 13;
      setDBVersion( version, sqlSession );
    }//if( version<13 && version<DB_SCHEMA_VERSION )

    if( version<14 && version<DB_SCHEMA_VERSION )
    {
      // Version 13 added `m_drfExtra`, but some databases ended up stamped as version 13 without
      //  ever receiving the column - the registry is written in a few places based only on whether
      //  `Session::createTables()` threw, so a database can be marked current without its schema
      //  actually matching.  Re-check, and add the column if it is missing.  A database that
      //  upgraded normally through the version 13 block above already has it, so this is a no-op.
      std::shared_ptr<Wt::Dbo::Session> sqlSession = getSession( database );

      addDrfExtraColumn( sqlSession );

      version = 14;
      setDBVersion( version, sqlSession );
    }//if( version<14 && version<DB_SCHEMA_VERSION )

    if( version<15 && version<DB_SCHEMA_VERSION )
    {
      std::shared_ptr<Wt::Dbo::Session> sqlSession = getSession( database );

      // Holds the Batch Decay tool's URI-encoded state so it can be saved/restored with app state.
      //  The following has only been checked for SQLite3; ADD COLUMN appends after the last column,
      //  which Dbo needs.
      const char *sql_statement = "ALTER TABLE UserState ADD COLUMN DecayBatchUri text;";
      executeSQL( sql_statement, sqlSession );

      version = 15;
      setDBVersion( version, sqlSession );
    }//if( version<15 && version<DB_SCHEMA_VERSION )

    /// ******************************************************************
    /// DB_SCHEMA_VERSION is at 15.  Add Version 16 here.  Update InterSpecUser.h!
    /// ******************************************************************
  }//void checkAndUpgradeVersion()
  
  
  std::shared_ptr<Wt::Dbo::Session> getSession( std::unique_ptr<Wt::Dbo::SqlConnection> & /*database*/ )
  {
    // In Wt4, Session::setConnection() takes ownership of the connection via unique_ptr.
    // We create a fresh connection for each session, letting the session own it.
    std::unique_ptr<Wt::Dbo::SqlConnection> conn( DataBaseUtils::getDatabaseConnection() );
    conn->setProperty( "show-queries", "true" );

    std::shared_ptr<Wt::Dbo::Session> sqlSession = std::make_shared<Wt::Dbo::Session>();
    sqlSession->setConnection( std::move(conn) );

    return sqlSession;
  }//shared_ptr<Session> getSession( std::unique_ptr<SqlConnection> &database )

  
  
  //Convenience function to execute SQL commands
  void executeSQL( const std::string &sql, std::shared_ptr<Wt::Dbo::Session> sqlSession)
  {
    Wt::Dbo::Transaction transaction( *sqlSession );
    sqlSession->execute(sql);
    transaction.commit();
  } //executeSQL(std::string sql, std::shared_ptr<Wt::Dbo::Session> m_sqlSession)


  bool tableHasColumn( const std::string &table, const std::string &column,
                       std::shared_ptr<Wt::Dbo::Session> sqlSession )
  {
    // Note: dont be tempted to probe this with something like
    //   `select count("column") from "table"`
    //  - SQLite treats a double-quoted identifier that doesnt resolve as a string literal, so
    //  that query happily succeeds for a column that isnt there.  Ask the catalog instead.
#if( USE_MYSQL_DB && !USE_SQLITE3_DB )
    const std::string sql = "select count(*) from information_schema.columns"
                            " where table_schema=database() and table_name='" + table + "'"
                            " and column_name='" + column + "'";
#elif( USE_SQLITE3_DB )
    const std::string sql = "select count(*) from pragma_table_info('" + table + "')"
                            " where name='" + column + "'";
#else
#error "Either a SQLITE3 or MySQL database must be selected"
#endif

    try
    {
      Wt::Dbo::Transaction transaction( *sqlSession );
      const int num_matching = sqlSession->query<int>( sql );
      transaction.commit();

      return (num_matching > 0);
    }catch( std::exception &e )
    {
      std::cerr << "Failed to check if '" << table << "' has column '" << column << "': "
                << e.what() << std::endl;
      return false;
    }
  }//bool tableHasColumn( const string &table, const string &column, shared_ptr<Session> )


  void addDrfExtraColumn( std::shared_ptr<Wt::Dbo::Session> sqlSession )
  {
    if( tableHasColumn( "DetectorPeakResponse", "m_drfExtra", sqlSession ) )
      return;

    try
    {
      // Note the '' (a real SQL string literal) rather than "" - SQLite only accepts the latter
      //  as an empty string by way of its double-quoted-identifier fallback, which a build with
      //  SQLITE_DQS=0 turns off.
      const char *sql_statement = "ALTER TABLE DetectorPeakResponse ADD COLUMN m_drfExtra text default '' not null;";
      executeSQL( sql_statement, sqlSession );
    }catch( std::exception &e )
    {
      // Dont let this abort start-up; without the column, saving DRFs to the database will fail,
      //  but the rest of the application is still usable.
      std::cerr << "Failed to add 'm_drfExtra' column to DetectorPeakResponse: " << e.what() << std::endl;

      // Only a developer-side concern if the table is actually present - a database that doesnt
      //  have the table at all will get it, with the column, from `Session::createTables()`.
      assert( !tableHasColumn( "DetectorPeakResponse", "m_name", sqlSession ) );
    }
  }//void addDrfExtraColumn( std::shared_ptr<Wt::Dbo::Session> sqlSession )


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
