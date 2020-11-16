#ifndef DataBaseUtils_h
#define DataBaseUtils_h
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
#include <memory>

namespace  Wt
{
  class WTimer;
  
  namespace Dbo
  {
    class Session;
    class Transaction;
    class SqlConnection;
    class SqlConnectionPool;
  }//namespace Dbo
}//namespace  Wt


#if( USE_MYSQL_DB )
//Before 20141014, we were having trouble keeping a single MySQL connection
//  pool open on Hekili, so we switched to creating a new connection everytime
//  a session is created.  This adds overhead whenever a new InterSpec
//  is created, but avoid the issue of stale connections.
//  Ideally we would like to fix the timeout issue with MySQL, so this
//  workaround doesnt have to be used, but its kinda a pain to diagnose and
//  debug this issue.
//  If you are using the SQLite3 database backend, then this problem doesnt
//  exist, so we will default to using a global connection pool to avoid this
//  extra overhead.
#define USE_GLOBAL_DATABASE_CONNECTION_POOL 0
#else
#define USE_GLOBAL_DATABASE_CONNECTION_POOL 1
#endif

namespace DataBaseUtils
{
  //Create a thread safe version of Wt::Dbo Session and Transaction.
  //  Note that not inheriting from the Dbo objects forces us to use this thread
  //  safe mechanism.
  class DbSession;

#if( USE_GLOBAL_DATABASE_CONNECTION_POOL )
  //databaseConnectionPool(): returns a pointer to the global database
  //  connection pool.
  //  The connection pool is lazily created the first time you call this
  //  function, and may throw an exception if the database cant be connected
  //  to.
  Wt::Dbo::SqlConnectionPool *databaseConnectionPool();
#endif
  
  //getDatabaseConnection(): creates a fresh database connection for you.
  //  Note that most of the time you will probably want to call
  //  create a new DbSession() object and not this function.
  //  You own the returned object, so make sure you delete it.
  //  Will throw exception if issue connecting to the database.
  Wt::Dbo::SqlConnection *getDatabaseConnection();
  
  
  class DbTransaction
  {
    //You can only have DbTransactions for a session in a single thread at a
    //  time
  public:
    DbTransaction( DbSession &session );
    ~DbTransaction() noexcept(false);
    
    bool commit();
    bool isActive() const;
    void rollback(); 
    
  protected:
    std::lock_guard<std::recursive_mutex> m_lock;
    Wt::Dbo::Transaction *m_transaction;
  };//class DbTransaction
  
  
  class DbSession
  {
    //DbSession is meant to be used instead of Wt::Dbo::Session everywhere.
    //  This class abstracts away the opening and connection to the database
    //  and the creation of the Wt::Dbo::Session.
  public:
    
    //DbSession(): constructs the DbSession object.
    //  If USE_GLOBAL_DATABASE_CONNECTION_POOL is true, then a the global
    //  database connection pool will be used, otherwise a new connection pool
    //  will be created (an share-owned) by this object.
    //Will throw exception on error.
    DbSession();
    
    //DbSession( DbSession & ): constructs a DbSession that shares the
    //  connection pool of 'rhs', but NOT the Wt::Dbo::Session.
    //Note, in the case USE_GLOBAL_DATABASE_CONNECTION_POOL is true, then this
    //  constructor functions identically to the default constructor.
    //This constructor is to avoid the overhead of creating a new connection
    //  to the database whenever a DetectorDisplay, DbStateBrowser,
    //  DbFileBrowser, or a UseInfoWindow instance is created, by allowing the
    //  use of the connection of the DbSession object the InterSpec owns;
    //  for the case of USE_GLOBAL_DATABASE_CONNECTION_POOL==0.
    DbSession( DbSession &rhs );
    
    ~DbSession();
    
    //session(): the returned Wt::Dbo::Session will already have all of our
    //  classes mapped and tables created; the returned object is owned by *this
    //  and should only be used when there is an active DbTransaction.
    Wt::Dbo::Session *session();
    
  protected:
    std::unique_ptr<Wt::Dbo::Session> m_session;
    std::recursive_mutex m_mutex;
    friend class DbTransaction;
    
#if( !USE_GLOBAL_DATABASE_CONNECTION_POOL )
    void keepDbConnectionAlive();
    Wt::WTimer *m_keepAliveTimer;  //should probably use a asio deadline_timer,
    std::shared_ptr<Wt::Dbo::SqlConnectionPool> m_connection;
#endif
  };//class DbSession
  
#if( USE_MYSQL_DB )
  //getMySqlDatabaseConnectionInfo(): A helper function to retrive
  void getMySqlDatabaseConnectionInfo( std::string &db,
                                       std::string &dbuser,
                                       std::string &dbpasswd,
                                       std::string &dbhost,
                                       unsigned int &port,
                                       std::string &dbsocket,
                                       const char* dbSource );
#endif  //#if( USE_MYSQL_DB )

  
#if( USE_SQLITE3_DB )
  //preferenceDatabaseFile(): location where the SQLITE3 database resides on the
  //  filesystem.  Note that a reference to a global variable is returned, since
  //  in principle it should only be changed by setPreferenceDatabaseFile(...)
  //  and only then before this function is ever called (this is dicey, but I
  //  think valid).
  const std::string &preferenceDatabaseFile();
  
  //setPreferenceDatabaseFile(...): should be set before any instances of
  //  InterSpecApp are created or preferenceDatabaseFile() or
  //  databaseConnectionPool() is called; if this is violated an exception will
  //  be thrown.
  void setPreferenceDatabaseFile( const std::string &dir );
#endif //#if( USE_SQLITE3_DB )
  
}//namespace DataBaseUtils

#endif //DataBaseUtils_h
