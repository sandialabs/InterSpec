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
#include <atomic>
#include <iostream>

#if( !USE_GLOBAL_DATABASE_CONNECTION_POOL )
#include <Wt/WTimer>
#endif 

#include <Wt/WApplication>
#include <Wt/Dbo/Session>
#include <Wt/Dbo/Exception>
#include <Wt/Dbo/Transaction>
#include <Wt/Dbo/SqlConnection>
#include <Wt/Dbo/FixedSqlConnectionPool>

#if( USE_MYSQL_DB )
#include <Wt/Dbo/backend/MySQL>
#endif

#if( USE_SQLITE3_DB )
#include <Wt/Dbo/backend/Sqlite3>
#endif

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"
#include "rapidxml/rapidxml_utils.hpp"

#include "InterSpec/InterSpecUser.h"
#include "InterSpec/DataBaseUtils.h"

using namespace std;
using namespace Wt;


namespace
{
#if( USE_SQLITE3_DB )
  std::mutex PreferenceDbMutex;
  int PreferenceDatabaseFileSetCount = 0;
  std::string PreferenceDatabaseFile = "InterSpecUserData.db";
#endif

#if( USE_GLOBAL_DATABASE_CONNECTION_POOL )
  //We need a mechanism to ensure that the database connection stays alive until
  //  the database is finished being written; this is what DbPoolManager does.
  class DbPoolManager
  {
    public:
    DbPoolManager()
    {
      m_numconnection.store( 0, std::memory_order_seq_cst );
    }
    
    ~DbPoolManager()
    {
      std::lock_guard<std::mutex> lock( m_mutex );
        
      const int nconn = m_numconnection.exchange( -1, std::memory_order_seq_cst );
      if( nconn > 0 )
      {
        vector<Dbo::SqlConnection *> connections;
        
        //By requesting as many connections as we created we are waitging for
        //  all Dbo::Transactions to return their Dbo::Session object, which
        //  happened and the end of the transaction, thus assuring all database
        //  operations in progress are complete.
        for( int i = 0; i < nconn; ++i )
          connections.push_back( m_pool->getConnection() );
      
        for( size_t i = 0; i < connections.size(); ++i )
          m_pool->returnConnection( connections[i] );
        
        m_pool.reset();
      }//if( nconn > 0 )
    }//~DbPoolManager()
    
    
    Wt::Dbo::SqlConnectionPool *databaseConnectionPool()
    {
      if( m_numconnection.load(std::memory_order_consume) <= 0 )
      {
        std::lock_guard<std::mutex> lock( m_mutex );
        
        const int nprevconn = m_numconnection.load(std::memory_order_consume);
        
        if( nprevconn == 0 )
        {
#if( defined(IOS) || defined(ANDROID) )
          const int nconn = 1;
#else
          const int nconn = 5;
#endif
          Dbo::SqlConnection *connection = DataBaseUtils::getDatabaseConnection();
          connection->setProperty( "show-queries", "false" );
          
          Dbo::FixedSqlConnectionPool *pool
                    = new Wt::Dbo::FixedSqlConnectionPool( connection, nconn );
          m_numconnection.store( nconn, std::memory_order_seq_cst );
          m_pool.reset( pool );
          
          Dbo::Session sesh;
          sesh.setConnectionPool( *pool );
          mapDbClasses( &sesh );
          try
          {
            //NOTE: this will almost always throw exception because it only
            //creates Tables if none of them already exist.  We are managing
            //and keeping the database always updated in main.cpp
            sesh.createTables();
            
            {//add version to registry
              Wt::Dbo::Transaction transaction( sesh );
              sesh.execute( "DELETE FROM InterSpecRegistry" );
              sesh.execute("INSERT INTO InterSpecRegistry(schema_version) VALUES (?)").bind(DB_SCHEMA_VERSION);
              transaction.commit();
            }
          }catch( std::exception & )
          {
            //    cerr << "Failed to createTables (if even one exists): " << e.what() << endl;
          }//try / catch to create tables
        }else if( nprevconn < 0 )
        {
          throw runtime_error( "Attempt to get DB Pool during program exit" );
        }//if( nprevconn == 0 ) / else if( nprevconn < 0 )
      }//if( !m_numconnection )
      
      assert( !!m_pool );
      
      return m_pool.get();
    }//Wt::Dbo::SqlConnectionPool *databaseConnectionPool()
    
  protected:
    std::mutex m_mutex;
    std::atomic<int> m_numconnection;
    std::unique_ptr<Wt::Dbo::FixedSqlConnectionPool> m_pool;
  };//class DbPoolManager
  
  DbPoolManager DbConnectionPoolManager;
#endif
}//namespace



namespace DataBaseUtils
{

DbTransaction::DbTransaction( DbSession &session )
  : m_lock( session.m_mutex )
{
  m_transaction = new Dbo::Transaction( *(session.m_session) );
}
  
DbTransaction::~DbTransaction()
{
  delete m_transaction;
  m_transaction = 0;
}

bool DbTransaction::commit()
{
  return m_transaction->commit();
}
  
bool DbTransaction::isActive() const
{
  return m_transaction->isActive();
}
 
void DbTransaction::rollback()
{
  m_transaction->rollback();
}
  
  
DbSession::DbSession( DbSession &rhs )
  : m_session( new Dbo::Session() )
{
#if( USE_GLOBAL_DATABASE_CONNECTION_POOL )
  m_session->setConnectionPool( *databaseConnectionPool() );
#else
  m_keepAliveTimer = 0;
  m_connection = rhs.m_connection;
  m_session->setConnectionPool( *m_connection );
#endif
  
  mapDbClasses( m_session.get() );
}//DbSession( DbSession &rhs )
  
  
DbSession::DbSession()
  : m_session( new Dbo::Session() )
{
  try
  {
#if( !USE_GLOBAL_DATABASE_CONNECTION_POOL )
    m_keepAliveTimer = 0;
    Wt::Dbo::SqlConnection *connection = getDatabaseConnection();
    m_connection.reset( new Dbo::FixedSqlConnectionPool( connection, 1 ) );
    m_session->setConnectionPool( *m_connection );
    mapDbClasses( m_session.get() );
    
    try
    {
      m_session->createTables();  //will throw if any tables already exist
      
      {//add version to registry
        Wt::Dbo::Transaction transaction( *m_session );
        m_session->execute( "DELETE FROM InterSpecRegistry" );
        m_session->execute("INSERT INTO InterSpecRegistry(schema_version) VALUES (?)").bind(DB_SCHEMA_VERSION);
        transaction.commit();
      }
    }catch( std::exception & )
    {
//    cerr << "Failed to createTables (if even one exists): " << e.what() << endl;
    }//try / catch to create tables
    
    if( wApp )
    {
      //Only make a keep alive timer if in a Wt event loop
      m_keepAliveTimer = new WTimer();
      m_keepAliveTimer->setInterval( 1000*60*60 );
      m_keepAliveTimer->timeout().connect( boost::bind( &DbSession::keepDbConnectionAlive, this ) );
      m_keepAliveTimer->start();
    }//if( wApp )
#else
    m_session->setConnectionPool( *databaseConnectionPool() );
    mapDbClasses( m_session.get() );
#endif
  }
  catch( Wt::Dbo::Exception &e )
  {
    throw std::runtime_error( "Failed to open SqlConnection to " + DataBaseUtils::preferenceDatabaseFile()
                              + ": " + string(e.what()) + ": code " + e.code());
  }catch( std::exception &e )
  {
    throw std::runtime_error( "Failed to open SqlConnection to " + DataBaseUtils::preferenceDatabaseFile() + ": " + string(e.what()) );
  }
  
}

  
DbSession::~DbSession()
{
#if( !USE_GLOBAL_DATABASE_CONNECTION_POOL )
  if( m_keepAliveTimer )
  {
    m_keepAliveTimer->stop();  //redundant?
    delete m_keepAliveTimer;
  }//if( m_keepAliveTimer )
  
  m_session.reset();
  m_connection.reset();
#endif
}

#if( !USE_GLOBAL_DATABASE_CONNECTION_POOL )
void DbSession::keepDbConnectionAlive()
{
  try
  {
    Wt::Dbo::Transaction transaction( *m_session );
      int version = m_session->query<int>("select schema_version from InterSpecRegistry");
      m_session->execute("DELETE FROM InterSpecRegistry WHERE schema_version>0 limit 100");
      m_session->execute("INSERT INTO InterSpecRegistry(schema_version) VALUES (?)").bind(version);
    transaction.commit();
  }catch(...)
  {
    cerr << "Caught exception keeping connection alive" << endl;
  }
}//void keepDbConnectionAlive()
#endif
  
  
Wt::Dbo::Session *DbSession::session()
{
  return m_session.get();
}
  
  
#if( USE_GLOBAL_DATABASE_CONNECTION_POOL )
Wt::Dbo::SqlConnectionPool *databaseConnectionPool()
{
  return DbConnectionPoolManager.databaseConnectionPool();
}//Wt::Dbo::SqlConnectionPool *databaseConnectionPool();
#endif
  
  
Wt::Dbo::SqlConnection *getDatabaseConnection()
{
#if( USE_MYSQL_DB && !USE_SQLITE3_DB )
  try
  {
    const char* databaseXML = "interSpec"; //default prod db
    if (strcmp(MYSQL_DATABASE_TO_USE,"dev")==0)
      databaseXML = "interSpecDevel"; //devl db
    else if (strcmp(MYSQL_DATABASE_TO_USE,"qc")==0)
      databaseXML = "interSpecQA"; //qc db
    
    unsigned int port;
    string db, user, passwd, host, socket;
    getMySqlDatabaseConnectionInfo( db, user, passwd, host, port, socket, databaseXML );
    return new Dbo::backend::MySQL( db, user, passwd, host, port, socket );
  }catch( std::exception &e )
  {
    cerr << "\n\nError initializing MySQL database connection for InterSpec: "
         << e.what() << endl << endl << endl;
    throw runtime_error( "Couldnt open database connection for user "
                        "preferences and storing spectrum files." );
  }//try / catch
#elif( USE_SQLITE3_DB )
  const string &dbLocation = DataBaseUtils::preferenceDatabaseFile();
  //Note that on Windows sqlite3 takes UTF-8 encoded path (not whatever codepage
  //  is currently defined)
  return new Wt::Dbo::backend::Sqlite3( dbLocation );
#else
#error "Either a SQLITE3 or MySQL database must be selected to use to store user preferences and spectra"
#endif
}//Wt::Dbo::SqlConnection *getDatabaseConnection()

  
#if( USE_MYSQL_DB )
  void getMySqlDatabaseConnectionInfo( string &db, string &dbuser,
                                                      string &dbpasswd, string &dbhost,
                                                      unsigned int &port, string &dbsocket, const char* dbSource )
  {
    using rapidxml::internal::compare;
    typedef rapidxml::xml_node<char>      XmlNode;
    typedef rapidxml::xml_attribute<char> XmlAttribute;
    
    //will throw runtime_error upon failure
    std::vector<char> data;
    UtilityFunctions::load_file_data( DATABASE_PASSWORD_FILE, data );
    
    rapidxml::xml_document<char> doc;
    const int flags = rapidxml::parse_normalize_whitespace
    | rapidxml::parse_trim_whitespace;
    
    doc.parse<flags>( &data.front() );
    XmlNode *doc_node = doc.first_node( "databases", 9 );
    
    if( !doc_node || !doc_node->name() )
      throw runtime_error( string("no databases node in ")
                          + DATABASE_PASSWORD_FILE );
    
    for( const XmlNode *node = doc_node->first_node( "database", 8 );
        node;
        node = node->next_sibling( "database", 8 ) )
    {
      const XmlAttribute *name = node->first_attribute( "name", 4 );
      
      //Load the appropriate database depending whether we are prod/qc/dev
      
      //NOTE: THIS CAUSES OUTRAGE TO CRASH FATALLY
      if( !name
         || !compare( name->value(), name->value_size(), dbSource, strlen(dbSource), 1) )
        continue;
      
      db = name->value();
      
      //We found the <database name="interSpec" dbm="MySQL"> node
      
      //TODO: figure out how dbsocket should actually be handled
      //      dbsocket = "/var/run/mysqld/mysqld.sock";
      dbsocket = "";
      
      //TODO: figure out how dbhost should actually be treated
      //      const XmlNode *ip_node = node->first_node( "ip", 2 );
      //      if( !ip_node || !ip_node->value() );
      //        throw runtime_error( "Invalid or missing ip node" );
      //      dbhost = ip_node->value();
#if( APPLE )
      dbhost = "LOCALHOST";
#else
      dbhost = "localhost";
#endif
      
      const XmlNode *port_node = node->first_node( "port", 4 );
      if( !port_node || !port_node->value()
         || !(stringstream(port_node->value())>>port) )
        throw runtime_error( "Invalid or missing port node" );
      
      const XmlNode *users = node->first_node( "users", 5 );
      if( !users )
        throw runtime_error( "Invalid or missing users node" );
      
      for( const XmlNode *user = users->first_node( "user", 4 );
          user;
          user = user->next_sibling( "user", 4 ) )
      {
        const XmlAttribute *role = user->first_attribute( "role", 4 );
        if( !role || !role->value()
           || !compare( role->value(), role->value_size(), "default", 7, 1 ) )
          continue;
        
        const XmlAttribute *id = user->first_attribute( "id", 2 );
        if( !id || !id->value() )
          throw runtime_error( "Missing id attribute in user node" );
        
        const XmlAttribute *password = user->first_attribute( "password", 8 );
        if( !password || !password->value() )
          throw runtime_error( "Missing password attribute in user node" );
        
        dbuser = id->value();
        dbpasswd = password->value();
        
        return;
      }//for( loop over users )
      
      throw runtime_error( "Unable to find user node for with 'default' role "
                          "in " + string(DATABASE_PASSWORD_FILE) );
    }//for( loop over database node )
    
    throw runtime_error( "Unable to find database node for interSpec database "
                        "in " + string(DATABASE_PASSWORD_FILE) );
  }//getMySqlDatabaseConnectionInfo(...)
#endif //#if( USE_MYSQL_DB )

  
#if( USE_SQLITE3_DB )
  const std::string &preferenceDatabaseFile()
  {
    std::lock_guard<std::mutex> lock( PreferenceDbMutex );
    ++PreferenceDatabaseFileSetCount;
    return PreferenceDatabaseFile;
  }//std::string getWritableResourceDirectory()
  
  
  void setPreferenceDatabaseFile( const std::string &dir )
  {
    std::lock_guard<std::mutex> lock( PreferenceDbMutex );
    
    if( dir == PreferenceDatabaseFile )
      return;
    
    if( PreferenceDatabaseFileSetCount )
      throw runtime_error( "You may only call setPreferenceDatabaseFile once" );
    
    ++PreferenceDatabaseFileSetCount;
    PreferenceDatabaseFile = dir;
  }//void InterSpecApp::setWritableResourceDirectory( std::string &dir )
#endif  //#if( USE_SQLITE3_DB )

  
}//namespace DataBaseUtils
