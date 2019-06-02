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

#include <fstream>
#include <iostream>
#include <iterator>
#include <string>

#include <Wt/Dbo/Dbo>
#include <Wt/Dbo/SqlConnection>

#if( USE_MYSQL_DB )
#include <Wt/Dbo/backend/MySQL>
#endif

#if( USE_SQLITE3_DB )
//sqlite3.h needed for the upgrade to verison 4 of the database
//#ifndef _WIN32
//Theres a good chance sqlite3.h might not be found on Windows.
//  wcjohns just hacked it by copying it from the Wt src directory
//  into the build folder
#include <sqlite3.h>
//#endif 

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
          std::cerr << "Trying initial table mapping, got exception: " << e.what() << std::endl;
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

    //Need to upgrade database
    if( version<1 && version<DB_SCHEMA_VERSION )
    {
/*
 //If we dont have have a version number, then perhaps we should just create all 
 //  the tables
      {//begin getting rid of all the tables
        std::cerr << "Database is from to old of an InterSpec version, "
                     "removing all user data" << std::endl;
        std::shared_ptr<Wt::Dbo::Session> sqlSession = getSession( database );
        mapDbClasses( sqlSession.get() );
        sqlSession->dropTables();
      }
      
      try
      {
        std::shared_ptr<Wt::Dbo::Session> sqlSession = getSession( database );
        mapDbClasses( sqlSession.get() );
        sqlSession->createTables();
        setDBVersion( DB_SCHEMA_VERSION, sqlSession );
        return;
      }catch( std::exception &e )
      {
        std::cerr << "Error crating database tables: " << e.what() << std::endl;
      }
*/
      
      std::shared_ptr<Wt::Dbo::Session> sqlSession = getSession( database );
      std::cerr<<"Upgrading to Version 1"<<std::endl;
        
      std::cerr<<"Creating InterSpecUser"<<std::endl;
#if( USE_SQLITE3_DB )
      executeSQL("CREATE TABLE IF NOT EXISTS InterSpecUser (id integer primary key autoincrement,userName text not null,isMobile boolean not null,accessCount integer not null,spectraFilesOpened integer not null,firstAccessUTC text,previousAccessUTC text,currentAccessStartUTC text,totalTimeInApp integer);",sqlSession);
#elif( USE_MYSQL_DB )
      executeSQL("CREATE TABLE IF NOT EXISTS InterSpecUser (id integer primary key AUTO_INCREMENT,userName text not null,isMobile boolean not null,accessCount integer not null,spectraFilesOpened integer not null,firstAccessUTC text,previousAccessUTC text,currentAccessStartUTC text,totalTimeInApp integer);",sqlSession);
#endif
        
      std::cerr<<"Creating ShieldingSourceModel"<<std::endl;
#if( USE_SQLITE3_DB )
      executeSQL("CREATE TABLE IF NOT EXISTS ShieldingSourceModel (id integer primary key autoincrement,version integer not null,InterSpecUser_id bigint,Name varchar(255) not null,Description varchar(511) not null,WriteProtected boolean not null,CreationTime text,SerializeTime text,XmlData text not null,constraint fk_ShieldingSourceModel_InterSpecUser foreign key (InterSpecUser_id) references InterSpecUser (id) on delete cascade);",sqlSession);
#elif( USE_MYSQL_DB )
      executeSQL("CREATE TABLE IF NOT EXISTS ShieldingSourceModel (id integer primary key AUTO_INCREMENT,version integer not null,InterSpecUser_id bigint,Name varchar(255) not null,Description varchar(511) not null,WriteProtected boolean not null,CreationTime text,SerializeTime text,XmlData text not null,constraint fk_ShieldingSourceModel_InterSpecUser foreign key (InterSpecUser_id) references InterSpecUser (id) on delete cascade);",sqlSession);
#endif
        
      std::cerr<<"Creating UserFileInDB"<<std::endl;
#if( USE_SQLITE3_DB )
      executeSQL("CREATE TABLE IF NOT EXISTS UserFileInDb (id integer primary key autoincrement,version integer not null,InterSpecUser_id bigint,SnapshotParent_id bigint,UploadTime text,UUID varchar(40) not null,Filename text not null,Description text not null,WriteProtected boolean not null,UserHasModified boolean not null,SessionID varchar(32) not null,NumSamples integer not null,IsPassthrough boolean not null,TotalLiveTime double precision not null,TotalRealTime double precision not null,TotalGammaCounts double precision not null,TotalNeutronCounts double precision not null,NumDetectors integer not null,HasNeutronDetector boolean not null,MeasurementsStartTime text,IsPartOfSaveState boolean not null,SerializeTime text,constraint fk_UserFileInDb_InterSpecUser foreign key (InterSpecUser_id) references InterSpecUser (id) on delete cascade,constraint fk_UserFileInDb_SnapshotParent foreign key (SnapshotParent_id) references UserFileInDb (id) on delete cascade);",sqlSession);
#elif( USE_MYSQL_DB )
      executeSQL("CREATE TABLE IF NOT EXISTS UserFileInDb (id integer primary key AUTO_INCREMENT,version integer not null,InterSpecUser_id bigint,SnapshotParent_id bigint,UploadTime text,UUID varchar(40) not null,Filename text not null,Description text not null,WriteProtected boolean not null,UserHasModified boolean not null,SessionID varchar(32) not null,NumSamples integer not null,IsPassthrough boolean not null,TotalLiveTime double precision not null,TotalRealTime double precision not null,TotalGammaCounts double precision not null,TotalNeutronCounts double precision not null,NumDetectors integer not null,HasNeutronDetector boolean not null,MeasurementsStartTime text,IsPartOfSaveState boolean not null,SerializeTime text,constraint fk_UserFileInDb_InterSpecUser foreign key (InterSpecUser_id) references InterSpecUser (id) on delete cascade,constraint fk_UserFileInDb_SnapshotParent foreign key (SnapshotParent_id) references UserFileInDb (id) on delete cascade);",sqlSession);
#endif
        
      std::cerr<<"Creating UserFileInDbData"<<std::endl;
#if( USE_SQLITE3_DB || USE_MYSQL_DB)
      executeSQL("CREATE TABLE IF NOT EXISTS UserFileInDbData (version integer not null,UserFileInDb_id bigint,gzipCompressed boolean not null,FileFormat integer not null,FileData blob not null,primary key (UserFileInDb_id),constraint fk_UserFileInDbData_UserFileInDb foreign key (UserFileInDb_id) references UserFileInDb (id) on delete cascade);",sqlSession);
#endif
        
      std::cerr<<"Creating UserOption"<<std::endl;
#if( USE_SQLITE3_DB )
      executeSQL("CREATE TABLE IF NOT EXISTS UserOption (id integer primary key autoincrement,name varchar(30) not null,value varchar(35) not null,type integer not null,InterSpecUser_id bigint,constraint fk_UserOption_InterSpecUser foreign key (InterSpecUser_id) references InterSpecUser (id) on delete cascade);",sqlSession);
#elif( USE_MYSQL_DB )
        executeSQL("CREATE TABLE IF NOT EXISTS UserOption (id integer primary key AUTO_INCREMENT,name varchar(30) not null,value varchar(35) not null,type integer not null,InterSpecUser_id bigint,constraint fk_UserOption_InterSpecUser foreign key (InterSpecUser_id) references InterSpecUser (id) on delete cascade);",sqlSession);
#endif
        
      std::cerr<<"Creating UserState"<<std::endl;
#if( USE_SQLITE3_DB )
      executeSQL("CREATE TABLE IF NOT EXISTS UserState (id integer primary key autoincrement,version integer not null,InterSpecUser_id bigint,StateType integer not null,Name varchar(127) not null,Description varchar(255) not null,WriteProtected boolean not null,CreationTime text,SerializeTime text,ForegroundId integer not null,BackgroundId integer not null,SecondForegroundId integer not null,ShieldSourceModelId integer not null,OtherSpectraCsvIds varchar(255) not null,ForegroundSampleNumsCsvIds varchar(127) not null,SecondForegroundSampleNumsCsvIds varchar(127) not null,BackgroundSampleNumsCsvIds varchar(127) not null,ShowingDetectorNumbersCsv varchar(127) not null,EnergyAxisMinimum double precision not null,EnergyAxisMaximum double precision not null,CountsAxisMinimum double precision not null,CountsAxisMaximum double precision not null,DisplayBinFactor integer not null,ShownDisplayFeatures integer not null,BackgroundSubMode integer not null,CurrentTab integer not null,GammaLinesXml text not null,IsotopeSearchEnergiesXml text not null,ShowingMarkers integer not null,DisabledNotifications integer not null,ShowingPeakLabels integer not null,ShowingWindows integer not null,UserOptions text not null,constraint fk_UserState_InterSpecUser foreign key (InterSpecUser_id) references InterSpecUser (id) on delete cascade);",sqlSession);
#elif( USE_MYSQL_DB )
      executeSQL("CREATE TABLE IF NOT EXISTS UserState (id integer primary key AUTO_INCREMENT,version integer not null,InterSpecUser_id bigint,StateType integer not null,Name varchar(127) not null,Description varchar(255) not null,WriteProtected boolean not null,CreationTime text,SerializeTime text,ForegroundId integer not null,BackgroundId integer not null,SecondForegroundId integer not null,ShieldSourceModelId integer not null,OtherSpectraCsvIds varchar(255) not null,ForegroundSampleNumsCsvIds varchar(127) not null,SecondForegroundSampleNumsCsvIds varchar(127) not null,BackgroundSampleNumsCsvIds varchar(127) not null,ShowingDetectorNumbersCsv varchar(127) not null,EnergyAxisMinimum double precision not null,EnergyAxisMaximum double precision not null,CountsAxisMinimum double precision not null,CountsAxisMaximum double precision not null,DisplayBinFactor integer not null,ShownDisplayFeatures integer not null,BackgroundSubMode integer not null,CurrentTab integer not null,GammaLinesXml text not null,IsotopeSearchEnergiesXml text not null,ShowingMarkers integer not null,DisabledNotifications integer not null,ShowingPeakLabels integer not null,ShowingWindows integer not null,UserOptions text not null,constraint fk_UserState_InterSpecUser foreign key (InterSpecUser_id) references InterSpecUser (id) on delete cascade);",sqlSession);
#endif
        
      std::cerr<<"Creating shielding_source_model"<<std::endl;
#if( USE_SQLITE3_DB || USE_MYSQL_DB)
      executeSQL("CREATE TABLE IF NOT EXISTS shielding_source_model (UserFileInDb_id bigint not null,ShieldingSourceModel_id bigint not null,primary key (UserFileInDb_id, ShieldingSourceModel_id),constraint fk_shielding_source_model_key1 foreign key (UserFileInDb_id) references UserFileInDb (id) on delete cascade,constraint fk_shielding_source_model_key2 foreign key (ShieldingSourceModel_id) references ShieldingSourceModel (id) on delete cascade);",sqlSession);
#endif
        
      version=1;
      setDBVersion( version, sqlSession );
    }// version 1
      
      
    //Each version after 1 should be updating columns, so it is incremental and
    //  NOT destroying existing tables
    if( version<2 && version<DB_SCHEMA_VERSION )
    {
      try
      {
        std::shared_ptr<Wt::Dbo::Session> sqlSession = getSession( database );
        
        try
        {
          executeSQL("DROP TABLE IF EXISTS DetectorPeakResponse;",sqlSession);
        }catch( std::exception & )
        {
        }
        
        std::cerr << "Added DetectorPeakResponse to database" << std::endl;
        sqlSession->mapClass<DetectorPeakResponse>( "DetectorPeakResponse" );
        sqlSession->createTables();
        version = 2;
        setDBVersion( version, sqlSession );
      }catch( std::exception &e )
      {
        std::cerr << "DB_SCHEMA_VERSION 2 : Failed to create DetectorPeakResponse table: "
                  << e.what() << std::endl;
      }
    }//if (version<2 && version<DB_SCHEMA_VERSION)
      
      
    if( version<3 && version<DB_SCHEMA_VERSION )
    {
      try
      {
        std::shared_ptr<Wt::Dbo::Session> sqlSession = getSession( database );
        executeSQL( "ALTER TABLE DetectorPeakResponse ADD Hash bigint not null DEFAULT '0'", sqlSession );
        executeSQL( "ALTER TABLE DetectorPeakResponse ADD ParentHash bigint not null DEFAULT '0'", sqlSession );
        version = 3;
        setDBVersion(version,sqlSession);
      }catch( std::exception &e )
      {
        std::cerr << "DB_SCHEMA_VERSION 3 : Failed to update DetectorPeakResponse table with Hash and PArentHash: "
                  << e.what() << std::endl;
      }//try / catch
    } //  if( version<3 && version<DB_VERSION )
      
      
    if( version<4 && version<DB_SCHEMA_VERSION )
    {
      //Change the InterSpecUser from having a BOOLEAN column named 'isMobile'
      //  to an INTEGER column named deviceType.
      //  This is to allow seperately tracking users phones and tablets.
      try
      {
        const char *sql = 0;

#if( USE_MYSQL_DB )
        {
          std::shared_ptr<Wt::Dbo::Session> sqlSession = getSession( database );
          sql = "ALTER TABLE InterSpecUser CHANGE isMobile deviceType INT";
          executeSQL( sql, sqlSession );
        }
#elif( USE_SQLITE3_DB )
		//There is a liniking issue for debug version of InterSpec when using shared libraries on Windows
#if( !defined(_DEBUG) || !(defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64)) )
        //Running with scissors:
        //http://stackoverflow.com/questions/805363/how-do-i-rename-a-column-in-a-sqlite-database-table
        
        sql = "PRAGMA writable_schema = 1;"
              "UPDATE SQLITE_MASTER SET SQL = 'CREATE TABLE InterSpecUser"
                "( id integer primary key autoincrement,"
                "userName text not null,"
                "deviceType int not null,"
                "accessCount integer not null,"
                "spectraFilesOpened integer not null,"
                "firstAccessUTC text,"
                "previousAccessUTC text,"
                "currentAccessStartUTC text,"
                "totalTimeInApp integer)"
              "' WHERE NAME = 'InterSpecUser';"
              "PRAGMA writable_schema = 0;";

        //For some reason if the above SQL is executaed using the Wt Dbo
        //  backend, it doesnt actually work - so here we'll drop down
        //  to using the SQLite3 API
        const std::string &dbLocation = DataBaseUtils::preferenceDatabaseFile();
        sqlite3 *db = 0;
        char *errmsg = 0;
        int rc = sqlite3_open( dbLocation.c_str(), &db );
        if( rc )
        {
          char buff[512];
          snprintf( buff, sizeof(buff), "Upgrade DB to V4: can't open DB: %s\n",
                    sqlite3_errmsg(db) );
          sqlite3_close(db);
          throw std::runtime_error( buff );
        }//if( rc )
          
        rc = sqlite3_exec( db, sql, NULL, 0, &errmsg );
        if( rc != SQLITE_OK )
        {
          char buff[512];
          snprintf( buff, sizeof(buff), "Upgrade DB to V4 error: %s", errmsg );
          sqlite3_free( errmsg );
          sqlite3_close( db );
          throw std::runtime_error( buff );
        }//if( rc != SQLITE_OK )
          
        sqlite3_close( db );
/*
 //The below is the "proper" way to upgrade the SQLite3 schema, however, it
 //  then gives a constraint failed error when adding user options into the
 //  the database in the InterSpec constructor - I assume due to the 
 //  autoincrementing of the primary key causing issues.
 //
 //SQLite3 cant rename columns or change types of existing tables, so
 //  we have to create a new table, and copy over the old one to the
 //  new one.

 //Update 20181106: Foreign key constraints would need to be temporarily turned off
 //  for the below "proper" method to work.  I.e., "PRAGMA foreign_keys = ON;"
 //  see https://www.sqlite.org/foreignkeys.html#fk_enable
 
          sql = "ALTER TABLE InterSpecUser RENAME TO preV4_InterSpecUser;";
          executeSQL( sql, sqlSession );
          
          sql = "CREATE TABLE InterSpecUser "
                "( id integer primary key autoincrement,"
                    "userName text not null,"
                    "deviceType int not null,"
                    "accessCount integer not null,"
                    "spectraFilesOpened integer not null,"
                    "firstAccessUTC text,"
                    "previousAccessUTC text,"
                    "currentAccessStartUTC text,"
                    "totalTimeInApp integer);";
          executeSQL( sql, sqlSession );
            
          sql = "INSERT INTO InterSpecUser("
                  "id,userName,deviceType,accessCount,spectraFilesOpened,"
                  "firstAccessUTC,previousAccessUTC,currentAccessStartUTC,"
                  "totalTimeInApp)"
                "SELECT id,userName,isMobile,accessCount,spectraFilesOpened,"
                  "firstAccessUTC,previousAccessUTC,currentAccessStartUTC,"
                  "totalTimeInApp "
                "FROM preV4_InterSpecUser;";
          executeSQL( sql, sqlSession );
          //Lets leave the old table in, JIC...
//          sql = "DROP TABLE preV4_InterSpecUser;";
//          executeSQL( sql, sqlSession );
*/
#endif  //if WIN32

#else
#error "checkAndUpgradeVersion(): only considers MySQL or SQLite3 database backends"
#endif
        
        std::shared_ptr<Wt::Dbo::Session> sqlSession = getSession( database );
        //Now make all the previous %_mobile user accounts into %_phone accounts,
        //  this assumes
        sql = "UPDATE InterSpecUser set userName=replace(userName,'mobile','phone')"
              " WHERE userName LIKE '%_mobile';";
        executeSQL( sql, sqlSession );
          
        version = 4;
        setDBVersion( version, sqlSession );
          
        std::cout << "Succesfully conveted InterSpecUser to have "
                       "INTEGER deviceType column, instead of "
                       "BOOLEAN isMobile columns" << std::endl
                    << "Database is now at version " << version << std::endl;
      }catch( std::exception &e )
      {
        std::cerr << "DB_SCHEMA_VERSION 4 : Failed to convert InterSpecUser: " << e.what() << std::endl;
      }//try / catch
    }//  if( version<4 && version<DB_VERSION )
    
    
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
      try
      {
        std::shared_ptr<Wt::Dbo::Session> sqlSession = getSession( database );
        executeSQL("DROP TABLE IF EXISTS InterSpecGlobalSetting;",sqlSession);
      }catch( std::exception & )
      {
        //We *should* probably get here, as the table shouldnt exist.
      }
      
      std::shared_ptr<Wt::Dbo::Session> sqlSession = getSession( database );
      sqlSession->mapClass<InterSpecGlobalSetting>( "InterSpecGlobalSetting" );
      sqlSession->createTables();
      
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
      
      /*
       create table "UseDrfPref" (
       "id" integer primary key autoincrement,
       "version" integer not null,
       "InterSpecUser_id" bigint,
       "MatchField" integer not null,
       "Flags" integer not null,
       "DrfIndex" bigint not null,
       "Criteria" text not null,
       constraint "fk_UseDrfPref_InterSpecUser" foreign key ("InterSpecUser_id") references "InterSpecUser" ("id") on delete cascade deferrable initially deferred
       )
       */
      
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
      
      
      /*
       create table "DetectorPeakResponse" (
       "id" integer primary key autoincrement,
       "version" integer not null,
       "m_name" varchar(255) not null,
       "m_description" varchar(255) not null,
       "m_detectorDiameter" real not null,
       "m_efficiencyEnergyUnits" real not null,
       "m_resolutionForm" integer not null,
       "m_resolutionCoeffs" text not null,
       "m_efficiencySource" integer not null,
       "m_efficiencyForm" integer not null,
       "m_energyEfficiencies" text not null,
       "m_efficiencyFormula" text not null,
       "m_expOfLogPowerSeriesCoeffs" text not null,
       "InterSpecUser_id" integer not null,
       "Hash" bigint not null,
       "ParentHash" bigint not null,
       "m_flags" bigint not null,
       "m_expOfLogPowerSeriesUncerts" text not null,
       "m_resolutionUncerts" text not null,
       "m_lowerEnergy" double precision not null,
       "m_upperEnergy" double precision not null,
       "m_createdUtc" bigint not null,
       "m_lastUsedUtc" bigint not null
       )
       */
      
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
    
    
    /// ******************************************************************
    /// DB_SCHEMA_VERSION is at 10.  Add Version 11 here.  Update InterSpecUser.h!
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
