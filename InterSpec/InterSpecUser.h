#ifndef InterSpecUser_h
#define InterSpecUser_h
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

#include <map>
#include <chrono>
#include <string>
#include <vector>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <stdexcept>

#include <boost/any.hpp>
#include <boost/function.hpp>
#include <boost/date_time.hpp>

#include <Wt/WSignal>
#include <Wt/Dbo/Dbo>
//#include <Wt/Dbo/weak_ptr>  //usable with Wt 3.3.1, but not currently using
#include <Wt/Dbo/SqlTraits>
#include <Wt/Dbo/collection>
#include <Wt/Dbo/WtSqlTraits>

#include "InterSpec/DataBaseUtils.h"

class SpecMeas;
class InterSpec;
struct UserState;
struct UseDrfPref;
class UserFileInDb;
class ColorThemeInfo;
class UserFileInDbData;
struct ShieldingSourceModel;


class DataStr : public std::vector<unsigned char>{ public: };
typedef DataStr FileData_t;

class UserOption;
class InterSpecUser;

namespace rapidxml
{
  template<class Ch> class xml_node;
}//namespace rapidxml


//The database this InterSpec is using; if higher than database registry, will
//  automatically update tables at next execution
//  See DataBaseVersionUpgrade.cpp/.h
#define DB_SCHEMA_VERSION 12


namespace Wt
{
  class WSpinBox;
  class WCheckBox;
  class WRadioButton;
  class WApplication;
  class WDoubleSpinBox;
  
  namespace Dbo
  {
    //Remove optimistic concurrency for UserOption and InterSpecUser classes
    //  to avoid "Stale object" exceptions when the same user has multiple
    //  sessions of the app open.  A related issue not dealt with is that
    //  one session will overwrite whats saved in another one, and stuff wont be
    //  propogated between them
    template<>
    struct dbo_traits<UserOption> : public dbo_default_traits
    {
      static const char *versionField() { return 0; }
    };
  
    template<>
    struct dbo_traits<InterSpecUser> : public dbo_default_traits
    {
      static const char *versionField() { return 0; }
    };
    
    //Specialize for OnoeToOne mapping between UserFileInDb and UserFileInDbData
    //  Note we are not removing optimistic concurrency from these classes,
    //  so you should be ready to catch the Wt::Dbo::StaleObjectException
    //  when saving to the database
    template<>
    struct dbo_traits<UserFileInDbData> : public dbo_default_traits
    {
      typedef ptr<UserFileInDb> IdType;
      static IdType invalidId() { return ptr<UserFileInDb>(); }
      static const char *surrogateIdField() { return 0; }
    };
    
    //We have to specialize sql_value_traits<FileData_t, void> so that for MySQL
    //  not just a BLOB will be used, but a MEDIUMBLOB (a pain!)
    template<>
    struct sql_value_traits<FileData_t, void>
    {
      static const bool specialized = true;
      static const char *type(SqlConnection *conn, int size);
      static void bind(const FileData_t& v,
                       SqlStatement *statement, int column, int size);
      static bool read(FileData_t& v, SqlStatement *statement,
                       int column, int size);
    };
    
    //specialize UserFileInDb to be indexed off of "UUID" and "Filename" as well
  }//namespace Dbo
}//namespace Wt


/** Converts from std::chrono to boost ptime.
 
 A (temporary?) function to help smooth converting the code from using boost::posix_time to
 std::chrono.
 */
boost::posix_time::ptime to_ptime( const std::chrono::system_clock::time_point &rhs );


/** Converts from boost ptime to std::chrono.
 
 A (temporary?) function to help smooth converting the code from using boost::posix_time to
 std::chrono.
 */
std::chrono::system_clock::time_point to_time_point( const boost::posix_time::ptime &rhs );


//mapDbClasses(...) Maps the database classes to Dbo::Session
void mapDbClasses( Wt::Dbo::Session *session );


class UserOption
{
public:
  enum DataType{ String, Decimal, Integer, Boolean };

  //Right now MySQL is only used for web deployments, and there are no string
  //  option values that should be very long, since none of them should
  //  coorespond to file paths.
  //  Non-web deployments have a few options (like File Query Tool path,
  //  detector response folder location), that include file paths, so we should
  //  allow longer preference values.  Note, SQLite does not impose any
  //  length restrictions on string fields, even if they are declared with
  //  something like VARCHAR(255).
#if( USE_MYSQL_DB )
  static const size_t sm_max_name_str_len = 30;
  static const size_t sm_max_value_str_len = 35;
#elif( USE_SQLITE3_DB )
  static const size_t sm_max_name_str_len = 30;
  static const size_t sm_max_value_str_len = 4096;
#else
  static_assert( 0, "Compile-time database type not recognized" );
#endif
  
  Wt::Dbo::ptr<InterSpecUser> m_user;
  
  DataType m_type;
  std::string m_name;
  std::string m_value;
  
  //value(): throws runtime_error if m_name or m_value is empty;
  boost::any value() const;
  
  template<class Action>
  void persist( Action &a )
  {
    Wt::Dbo::field( a, m_name, "name", sm_max_name_str_len );
    Wt::Dbo::field( a, m_value, "value", sm_max_value_str_len );
    Wt::Dbo::field( a, m_type, "type" );
    Wt::Dbo::belongsTo( a, m_user, "InterSpecUser", Wt::Dbo::OnDeleteCascade );
  }//void persist( Action &a )
};//class UserOption


class InterSpecUser
{
  //Right now access to the InterSpecUser class should only be done from the
  //  the main gui thread or when you have a active DataBaseUtils::DbTransaction
  //  in order to avoid potential multithreading issue.
public:
  enum DeviceType
  {
    Desktop = 0x0,
    PhoneDevice = 0x1,
    TabletDevice = 0x2
  };

public:
  //InterSpecUser(): default constructor necessary for Wt::Dbo to be able to
  //  construct InterSpecUser objects to later call persist(...) on to read
  //  from the database.  The conseuence of this is if the userName() of a
  //  InterSpecUser object is empty, it is an invalid object.
  InterSpecUser();
  
  //InterSpecUser(): constructor for a new user.  No preferences will be
  //  associated for this user.  You should call initFromDefaultValues(...)
  //  in order to do this.
  InterSpecUser( const std::string &username, DeviceType type );
  
  const std::string &userName() const;
  
  /** The number of new InterSpec sessions that has been created for this user. */
  int accessCount() const;
  
  /** Approximate total amount of time user has actively used the app before the current session.
   
   Counts time between user generated events in the app, but only if less than 5 minutes between interactions.
   See #InterSpecApp::notify for implementation of when time is added.
   
   This value is only incremented #InterSpecApp::prepareForEndOfSession, so it does not include time from the current session.
   */
  std::chrono::steady_clock::time_point::duration totalTimeInApp() const;
  
  /** Approximate number of spectrum files opened.
   
   The count is only incremented if #SpecMeasManager::checkIfPreviouslyOpened is called after opening a file, which is only
   called some of the time from #SpecMeasManager::displayFile, specifically I think only when files are opened from the
   filesystem.
   */
  int numSpectraFilesOpened() const;
  
  /** The UTC time when this user was first created in the database. */
  std::chrono::system_clock::time_point firstAccessUTC() const;
  
  //userFromViewer: a simple helper function to return the user from a
  //  spectrum viewer pointer; necessary to break a "Member access into
  //  incomplete type 'InterSpec'" issue.
  static Wt::Dbo::ptr<InterSpecUser> &userFromViewer( InterSpec *viewer );
  
  //sqlFromViewer: same story as userFromViewer()
  static std::shared_ptr<DataBaseUtils::DbSession> sqlFromViewer( InterSpec *viewer );
  
  
  //preferenceValueAny(...): Retrieves preference value.  If value has not
  //  been previously set, will return default value as well as adding this
  //  value to to the users preferences.  If no preference by the passed in name
  //  is found, nor one with a default value, a runtime_error will be thrown.
  //  There must be a Dbo::Session associated with usr or else and exception
  //  will be thrown.
  //  The InterSpec is necessary in order to add the default value to the
  //  database in the case the user doesnt already have a preference by the 
  //  given name, but that preference does exist in the default values.
  static boost::any preferenceValueAny( const std::string &name, InterSpec *viewer );
  
  /** Convienince function to call for #preferenceValueAny */
  template<typename T>
  static T preferenceValue( const std::string &name,
                            InterSpec *viewer );
  
  //preferenceValue(...): similar to above, but if user doesnt already have a
  //  value for the desired preference, an exception will be thrown, even if
  //  a value would exist in the default values file. Can be called at any time.
  boost::any preferenceValueAny( const std::string &name ) const;
  
  template<typename T>
  T preferenceValue( const std::string &name ) const;
  
  //setPreferenceValue(): Sets preference value for named preference to both
  //  the InterSpecUser in memory and the database.
  //  If the preference isnt already in memory, and its not in
  //  data/default_preferences.xml then will throw an exception
  //  Note the reason for passing InterSpecUser and viewer pointers is so
  //  we dont have to #include InterSpec.h in this file.
  template<typename T>
  static void setPreferenceValue( Wt::Dbo::ptr<InterSpecUser> user,
                                  const std::string &name, const T &value,
                                  InterSpec *viewer );
  
  /** Specialization of #setPreferenceValue for boolean preference; does same as other variant
   of the function, but also calls functions in m_onBoolChangeSignals.
   */
  static void setPreferenceValue( Wt::Dbo::ptr<InterSpecUser> user,
                                 const std::string &name,
                                 const bool &value,
                                 InterSpec *viewer );
  
  /** Just a convenience function to call above variant, for binding to signals. */
  static void setBoolPreferenceValue( Wt::Dbo::ptr<InterSpecUser> user,
                                 const std::string &name,
                                 const bool &value,
                                 InterSpec *viewer );
  
protected:
  template<typename T>
  static void setPreferenceValueWorker( Wt::Dbo::ptr<InterSpecUser> user,
                                 const std::string &name, const T &value,
                                 InterSpec *viewer );

public:
  
  /** Add a function to callback when the preference changes, either through loading  a new state,
   or toggling the preference checkbox or whatever.
   
   This variant of #addCallbackWhenChanged is useful to book a slot (member function) of a
   Wt::WObject up to, so the callback lifetime will be limited to the lifetime of the WObject.
   
   This function must be called from the application thread such that InterSpec::instance() will be
   non-void.
   */
  template<class T, class V>
  static void addCallbackWhenChanged( Wt::Dbo::ptr<InterSpecUser> &user,
                                      const std::string &name, T *target,
                                      void (V::*method)(bool) );
  
  /** Adds a function to callback when the preference changes, either through loading
   a new state, or toggling the preference checkbox or whatever.
   
   This function must be called from the application thread such that InterSpec::instance() will be
   non-void
   
   Note, the function passed in must be sure that it will be safe to call until the end of the users
   session.  If you need limit lifetime of callbacks, either make sure your caller connects to
   a slot of a Wt::WObject, or is otherwise somehow safe for the entirety of the InterSpec class
   lifetime.
   */
  template<class T>
  static void addCallbackWhenChanged( Wt::Dbo::ptr<InterSpecUser> &user,
                                     const std::string &name,
                                     const T &fcn );
  
  /** Makes it so if the user changes the value via this GUI element, the value of the preference
   stored in memory and in the database will be correspondingly updated.  Also makes it so if
   a different GUI element updates this preference value, or #setPreferenceValue is called for
   this preference, this widgets state will be correspondingly updated as well.
   
   However, the signals you have hooked up to this widget wont be called when the value changes by
   another widget or #setPreferenceValue, to handle this case, use the #addCallbackWhenChanged
   function to add a callback for whenever the value of the preference changes.
   */
  static void associateWidget( Wt::Dbo::ptr<InterSpecUser> user,
                               const std::string &name,
                               Wt::WCheckBox *cb,
                               InterSpec *viewer );
  
  /*
  static void associateWidget( Wt::Dbo::ptr<InterSpecUser> user,
                               const std::string &name,
                               Wt::WSpinBox *spinner,
                               InterSpec *viewer );
  static void associateWidget( Wt::Dbo::ptr<InterSpecUser> user,
                               const std::string &name,
                               Wt::WDoubleSpinBox *spinner,
                               InterSpec *viewer );
   */
  
  
  
  //startingNewSession(): increments access counts, and updates previous and
  //  current timestamps.  Should be called at the begining of a new InterSpec.
  //  Should only be called from within an active transaction.
  void startingNewSession();
 
  //addUsageTimeDuration(): add a access time duration.
  //  Should only be called from within an active transaction.
  void addUsageTimeDuration( const std::chrono::steady_clock::time_point::duration &duration );
 
  //incrementSpectraFileOpened(): increments m_spectraFilesOpened, and is
  //  intended to be called whenever a new spectrum file is opened.
  //  Should only be called from within an active transaction.
  //  Currently only called from SpecMeasManager::checkIfPreviouslyOpened(...)
  //  (e.g. when a user selects to display a new file)
  void incrementSpectraFileOpened();
  
  //userFiles(): useful to a users previously saved files
  const Wt::Dbo::collection< Wt::Dbo::ptr<UserFileInDb> > &userFiles() const;
  
  //userOptionsToXml(): serializes the current user options to a node called
  //  <preferences> under parent passed in, with each preferences being placed
  //  in a <pref> node under it.  Uses same schema as 'default_preferences.xml'.
  //Returns the <preferences> node.
  ::rapidxml::xml_node<char> *userOptionsToXml(
                                            ::rapidxml::xml_node<char> *parent,
                                            InterSpec *viewer ) const;
  
  //restoreUserPrefsFromXml(): sets the currently active user preferences from
  //  passed in XML node.  'preferences' must be a node named <preferences>
  //  and should be in same format that userOptionsToXml() creates.
  //  If there is a current user option that is not in the passd in XML, its
  //  value will not be altered.  Any new options in the XML not currently in
  //  memory will be set.
  //Throws exception upon error.
  static void restoreUserPrefsFromXml( Wt::Dbo::ptr<InterSpecUser> user,
                                  const ::rapidxml::xml_node<char> *preferences,
                                InterSpec *viewer );
  
  template<class Action>
  void persist( Action &a )
  {
    Wt::Dbo::field( a, m_userName, "userName" );
    Wt::Dbo::field( a, m_deviceType, "deviceType" );
    Wt::Dbo::field( a, m_accessCount, "accessCount" );
    Wt::Dbo::field( a, m_spectraFilesOpened, "spectraFilesOpened" );
    Wt::Dbo::field( a, m_firstAccessUTC, "firstAccessUTC" );
    Wt::Dbo::field( a, m_previousAccessUTC, "previousAccessUTC" );
    Wt::Dbo::field( a, m_currentAccessStartUTC, "currentAccessStartUTC" );
    Wt::Dbo::field( a, m_totalTimeInApp, "totalTimeInApp" );  //shows up as milliseconds in the sqlite3 database it looks like.
    
    Wt::Dbo::hasMany( a, m_userFiles, Wt::Dbo::ManyToOne, "InterSpecUser" );
    Wt::Dbo::hasMany( a, m_dbPreferences, Wt::Dbo::ManyToOne, "InterSpecUser" );
    Wt::Dbo::hasMany( a, m_shieldSrcModels, Wt::Dbo::ManyToOne,"InterSpecUser");
    Wt::Dbo::hasMany( a, m_colorThemes, Wt::Dbo::ManyToOne,"InterSpecUser");
    Wt::Dbo::hasMany( a, m_drfPref, Wt::Dbo::ManyToOne,"InterSpecUser");
  }//void persist( Action &a )
  
  std::chrono::system_clock::time_point currentAccessStartUTC() const;
  
  const Wt::Dbo::collection< Wt::Dbo::ptr<ShieldingSourceModel> > &shieldSrcModels() const;
  const Wt::Dbo::collection< Wt::Dbo::ptr<UserState> > &userStates() const;
  const Wt::Dbo::collection< Wt::Dbo::ptr<ColorThemeInfo> >   &colorThemes() const;
  const Wt::Dbo::collection< Wt::Dbo::ptr<UseDrfPref> > &drfPrefs() const;
  
protected:
  void incrementAccessCount();
  void setCurrentAccessTime( const std::chrono::system_clock::time_point &utcTime );
  void setPreviousAccessTime( const std::chrono::system_clock::time_point &utcTime );
 
  //initFromDefaultValues(...): Set prefernces will be set to default values.
  //  There must be an active transaction associated with the session passed in.
  //  Will throw if default values XML file (m_defaultPreferenceFile) is not
  //  found, or is invalid or ill-formatted, or if the user already has any
  //  preferences.
  static void initFromDefaultValues( Wt::Dbo::ptr<InterSpecUser> user,
                          std::shared_ptr<DataBaseUtils::DbSession> session );
  
  static void initFromDbValues( Wt::Dbo::ptr<InterSpecUser> user,
                          std::shared_ptr<DataBaseUtils::DbSession> session );
  
  
  //getDefaultUserPreference(...): will throw exception upon error, otherwise
  //  results will always be valid.
  //Will search for user option specialized for DeviceType (represented by the
  //  int 'type') before returning the general option
  static UserOption *getDefaultUserPreference( const std::string &name,
                                               const int type );
 
  typedef std::map<std::string,boost::any> PreferenceMap;
  
  std::string m_userName;
  int m_deviceType;
  int m_accessCount;
  int m_spectraFilesOpened;
  boost::posix_time::ptime m_firstAccessUTC;
  boost::posix_time::ptime m_previousAccessUTC;
  boost::posix_time::ptime m_currentAccessStartUTC;
  boost::posix_time::time_duration m_totalTimeInApp;
  
  //These are mutable so Dbo::ptr<t>.modify() dont need to be called
  mutable PreferenceMap m_preferences;
  
  /** Holds callbacks set from #addCallbackWhenChanged, for boolean preferences.
   
   I believe using a Wt::Signals::signal allows makes it so we can connect widget slots (function
   calls), and then the connection will automatically get deleted when the widget is deleted, making
   things safe.
   
   Note that these function may be "safe" for widgets getting deleted (e.g., signals automatically
   disconnected), if the signal is connected to an object deriving from Wt::WObject.
   
   Currently only boolean preferences require callbacks in InterSpec; in the future, if other
   preference types (string, int, doubles) need callbacks, we will have to re-factor things, or
   add in analogous member variables.
   
   Mutable so Dbo::ptr<t>.modify() doesnt need to be called
   */
  mutable std::map<std::string,std::shared_ptr<Wt::Signals::signal<void(bool)>>> m_onBoolChangeSignals;
  
  
  Wt::Dbo::collection< Wt::Dbo::ptr<UserFileInDb> >         m_userFiles;
  Wt::Dbo::collection< Wt::Dbo::ptr<UserOption> >           m_dbPreferences;
  Wt::Dbo::collection< Wt::Dbo::ptr<ShieldingSourceModel> > m_shieldSrcModels;
  Wt::Dbo::collection< Wt::Dbo::ptr<ColorThemeInfo> >       m_colorThemes;
  Wt::Dbo::collection< Wt::Dbo::ptr<UserState> >            m_userStates;

  Wt::Dbo::collection< Wt::Dbo::ptr<UseDrfPref> >           m_drfPref;
  
  
//  Wt::Dbo::collection< Wt::Dbo::ptr<DbUserState::Spectrum> > m_spectrums;
  
  /** The file name of the default preferences XML file.  File name is relative
      to InterSpec::dataDirectory(), and has a default value of
      "default_preferences.xml".
   */
  static const std::string sm_defaultPreferenceFile;
  
  friend class InterSpec;
};//struct UserOption


//We will seperate UserFileInDb and UserFileInDbData classes to maybye speed
//  up the queries so the spectrum data wont be loaded during searches
class UserFileInDb
{
public:
  //sm_maxFileSizeBytes: maximum size of spectrum file we'll save to the
  //  database.  Dictated by MySQL medium blob size (a choice in
  //  sql_value_traits<FileData_t>).
  const static size_t sm_maxFileSizeBytes = 16777215;
  const static int sm_maxUuidLength = 40;
  const static int sm_maxSessionIdLength = 32;

public:
  UserFileInDb();
  
  Wt::Dbo::ptr<InterSpecUser> user;
  
  std::string uuid;
  std::string filename;
  std::string description;

  bool userHasModified;
  Wt::WDateTime uploadTime;
  Wt::WDateTime serializeTime;
  std::string sessionID;
  
  int numSamples;
  bool isPassthrough;
  double totalLiveTime;
  double totalRealTime;
  double totalGammaCounts;
  double totalNeutronCounts;
  int numDetectors;
  bool hasNeutronDetector;
  Wt::WDateTime measurementsStartTime;
  
  Wt::Dbo::ptr<UserFileInDb> snapshotParent;
  Wt::Dbo::collection< Wt::Dbo::ptr<UserFileInDb> > snapshots;
  
  Wt::Dbo::collection< Wt::Dbo::ptr<UserFileInDbData> > filedata;
//Could use below instead of collection with newer Wt
//  Wt::Dbo::weak_ptr<UserFileInDbData> filedata;
  
  Wt::Dbo::collection<Wt::Dbo::ptr<ShieldingSourceModel> > modelsUsedWith;
  
  //isWriteProtected(): returns if this UserFileInDb has been made write
  //  protected
  bool isWriteProtected() const;
  
  //makeWriteProtected(...): sets the writeprotected flag so object wont
  //  be re-saved to the database in an altered state.  Should only be called
  //  from within an active transaction
  static void makeWriteProtected( Wt::Dbo::ptr<UserFileInDb> ptr );
  
  //removeWriteProtection(...): sets the writeprotected flag so object can
  //  be re-saved to the database.  Should only be called from within an active
  //  transaction
  static void removeWriteProtection( Wt::Dbo::ptr<UserFileInDb> ptr );
  
  //makeDeepWriteProtectedCopyInDatabase(...): creates a new UserFileInDb object
  //  in the database, including a new UserFileInDbData entry, as well as a any
  //  ShieldingSourceModels which arent write proteced.
  static Wt::Dbo::ptr<UserFileInDb> makeDeepWriteProtectedCopyInDatabase(
                                            Wt::Dbo::ptr<UserFileInDb> orig,
                                            DataBaseUtils::DbSession &session,
                                            bool isSaveState );
  
  
  template<class Action>
  void persist( Action &a )
  {
    Wt::Dbo::belongsTo( a, user, "InterSpecUser", Wt::Dbo::OnDeleteCascade );
    
    Wt::Dbo::belongsTo( a, snapshotParent, "SnapshotParent", Wt::Dbo::OnDeleteCascade );
    Wt::Dbo::hasMany( a, snapshots, Wt::Dbo::ManyToOne, "SnapshotParent" );
    
    Wt::Dbo::hasMany( a, filedata, Wt::Dbo::ManyToOne, "UserFileInDb" );
//    Wt::Dbo::hasOne( a, filedata, "UserFileInDb" );
    Wt::Dbo::hasMany( a, modelsUsedWith,
                      Wt::Dbo::ManyToMany, "shielding_source_model" );
    
    Wt::Dbo::field( a, uploadTime,      "UploadTime" );
    Wt::Dbo::field( a, uuid,            "UUID", sm_maxUuidLength );
    Wt::Dbo::field( a, filename,        "Filename" );
    Wt::Dbo::field( a, description,     "Description" );
    Wt::Dbo::field( a, writeprotected,  "WriteProtected" );
    Wt::Dbo::field( a, userHasModified, "UserHasModified" );
    Wt::Dbo::field( a, sessionID,       "SessionID", sm_maxSessionIdLength );
    
    Wt::Dbo::field( a, numSamples,            "NumSamples" );
    Wt::Dbo::field( a, isPassthrough,         "IsPassthrough" );
    Wt::Dbo::field( a, totalLiveTime,         "TotalLiveTime" );
    Wt::Dbo::field( a, totalRealTime,         "TotalRealTime" );
    Wt::Dbo::field( a, totalGammaCounts,      "TotalGammaCounts" );
    Wt::Dbo::field( a, totalNeutronCounts,    "TotalNeutronCounts" );
    Wt::Dbo::field( a, numDetectors,          "NumDetectors" );
    Wt::Dbo::field( a, hasNeutronDetector,    "HasNeutronDetector" );
    Wt::Dbo::field( a, measurementsStartTime, "MeasurementsStartTime" );
    Wt::Dbo::field( a, isPartOfSaveState,     "IsPartOfSaveState" );
    
    
    if( a.getsValue() )
    {
      serializeTime = Wt::WDateTime::currentDateTime();
      Wt::Dbo::field( a, serializeTime, "SerializeTime" );
    }//if( saving to DB )
    
    if( a.setsValue() || a.isSchema() )
    {
      Wt::Dbo::field( a, serializeTime, "SerializeTime" );
    }//if( reading from DB )
  }//void persist( Action &a )
  
private:
  //writeprotected: not yet fully implemented, but intended to allow a way to
  //  ensure database entry doesnt get changed in the future.  This is mostly
  //  relevant for saving the applications state, you dont want the user to
  //  change the individual components, thus corrupting the state
  bool writeprotected;
  
  bool isPartOfSaveState;
};//class UserFileInDb


#if( HAS_ZLIB_SUPPORT && !IOS && !ANDROID )
//Some example compression ratios for boost binary files:
//  compressed: 367780 bytes  vs 659168   bytes; savings: 44% (Ba133 example file)
//  compressed: 4051514 bytes vs 70129153 bytes; savings: 94% (passthrough example file)
//  compressed: 4817 bytes    vs 7775     bytes; savings: 38% (1024 channel file)
//  compressed: 130361 bytes  vs 197899   bytes; savings: 34% ("2011_11_07_15_31_130 - airplane at altitude.SPE")
//  compressed: 65806 bytes   vs 99409    bytes; savings: 34% ("Alphas on Boron.Chn")
//  compressed: 94017 bytes   vs 198993   bytes; savings: 53% ("detector problem at 548.spc")
//  compressed: 66315 bytes   vs 99429    bytes; savings: 33% ("fertilizer_TexasA&M.chn")
//I have not yet looked at how much, or if, time the compressing saves when
//  placing into or reading from the database

//ALLOW_SAVE_TO_DB_COMPRESSION: use gzip compression to save to database.
//  Could probably allow for iOS and Android, but I havent tested this...
#define ALLOW_SAVE_TO_DB_COMPRESSION 0

#endif //HAS_ZLIB_SUPPORT

 
class FileToLargeForDbException : public std::exception
{
  std::string m_msg;
public:
  FileToLargeForDbException( const size_t saveSize, const size_t limit )
    : std::exception()
  {
    char msg[100];
    snprintf( msg, sizeof(msg),
              "Spectrum file save size is %i kb; the limit is %i kb.",
              static_cast<int>(saveSize/1024), static_cast<int>(limit/1024) );
    m_msg = msg;
  }

  
  virtual const char* what() const noexcept
  {
    return m_msg.c_str();
  }
};//class FileToLargeForDbException

class UserFileInDbData
{
public:
  enum SerializedFileFormat
  {
    k2012N42
  };//enum SerializedFileFormat
  
  const static SerializedFileFormat sm_defaultSerializationFormat;
  
public:
  UserFileInDbData();
  Wt::Dbo::ptr<UserFileInDb> fileInfo;
  
  //gzipCompressed: only set true when writing if ALLOW_SAVE_TO_DB_COMPRESSION
  //  has been set.
  bool gzipCompressed;
  
  //So we can choose to save files to the datbase in something besides native
  //  binary format
  SerializedFileFormat fileFormat;
  
  //fileData: the actual data of the serialized SpecMeas object, may be
  //  compressed.
  FileData_t fileData;
  
  //setFileData(...): serializes spectrumFile to fileData as a binary native
  //  file format.
  //  Will throw FileToLargeForDbException if serialization is larger than
  //  UserFileInDb::sm_maxFileSizeBytes.
  //  Will throw runtime_error if any other issues.
  void setFileData( std::shared_ptr<const SpecMeas> spectrumFile,
                    const SerializedFileFormat format );

  //setFileData(...): same as other setFileData(...) function, but instead
  //  sets the data from a file on the filesystem.  Will throw if the file
  //  reading fails for any reason.
  void setFileData( const std::string &path,
                    const SerializedFileFormat format );
  
  //decodeSpectrum(): de-serializes data currently in fileData.
  //  Will throw if de-serialization fails, otherwise will always return
  //  a valid SpecMeas object.
  std::shared_ptr<SpecMeas> decodeSpectrum() const;
  
  template<class Action>
  void persist( Action &a )
  {
    Wt::Dbo::id( a, fileInfo, "UserFileInDb", Wt::Dbo::OnDeleteCascade );
    Wt::Dbo::field( a, gzipCompressed, "gzipCompressed" );
    Wt::Dbo::field( a, fileFormat, "FileFormat" );
    Wt::Dbo::field( a, fileData, "FileData" );
  }//void persist( Action &a )
};//class UserFileInDbData


class ColorThemeInfo
{
  /** The current user color theme is determined by a (int) user preference,
      "ColorThemeIndex", which holds the index of of the of the entry in this
      table to use.  The 'InterSpecUser_id' must also match current user.
   */
public:
  Wt::Dbo::ptr<InterSpecUser> user;
  
  Wt::WString theme_name;
  Wt::WString theme_description;
  
  Wt::WDateTime creation_time;
  Wt::WDateTime modified_time;
  
  /** The JSON representation of the ColorTheme class. */
  std::string json_data;
  
  template<class Action>
  void persist( Action &a )
  {
    Wt::Dbo::belongsTo( a, user, "InterSpecUser", Wt::Dbo::OnDeleteCascade );
    Wt::Dbo::field( a, theme_name, "ThemeName" );
    Wt::Dbo::field( a, theme_description, "ThemeDescription" );
    
    //For some reason on Anddroid I've gotten an exception before - I guess the
    //  timestamps arent that important so we'll just live with it.
    try{
      Wt::Dbo::field( a, creation_time, "CreationTime" );
    }catch(...){
      std::cerr << "Caught exception getting ColorThemeInfo creation time from database" << std::endl;
    }
    try{
      Wt::Dbo::field( a, modified_time, "ModifiedTime" );
    }catch(...){
      std::cerr << "Caught exception getting ColorThemeInfo modified time from database" << std::endl;
    }
    
    Wt::Dbo::field( a, json_data, "JsonData" );
  }//void persist( Action &a )

  
};//class ColorThemeInfo


struct ShieldingSourceModel
{
  ShieldingSourceModel()
    : writeprotected(false)
  {}
  
  Wt::Dbo::ptr<InterSpecUser> user;
  Wt::WString name;
  Wt::WString description;
    
  Wt::WDateTime creationTime;
  Wt::WDateTime serializeTime;

  Wt::Dbo::collection<Wt::Dbo::ptr<UserFileInDb> > filesUsedWith;
  
  std::string xmlData;
  
  bool isWriteProtected() const;

  //makeWriteProtected(...): an active transaction must exist when calling this
  //  function
  static void makeWriteProtected( Wt::Dbo::ptr<ShieldingSourceModel> ptr );
  
  //removeWriteProtection(): an active transaction must exist when calling this
  //  function
  static void removeWriteProtection( Wt::Dbo::ptr<ShieldingSourceModel> ptr );

  
  void shallowEquals( const ShieldingSourceModel &rhs );
  
  template<class Action>
  void persist( Action &a )
  {
    Wt::Dbo::belongsTo( a, user, "InterSpecUser", Wt::Dbo::OnDeleteCascade );
    Wt::Dbo::hasMany( a, filesUsedWith,
                      Wt::Dbo::ManyToMany, "shielding_source_model" );
  
    Wt::Dbo::field( a, name, "Name", 255 );
    Wt::Dbo::field( a, description, "Description", 511 );
    Wt::Dbo::field( a, writeprotected, "WriteProtected" );
    
    Wt::Dbo::field( a, creationTime, "CreationTime" );
    Wt::Dbo::field( a, serializeTime, "SerializeTime" );
    Wt::Dbo::field( a, xmlData, "XmlData" );
  }//void persist( Action &a )
  
private:
  //writeprotected: not yet implemented, but intended to allow a way to ensure
  //  database entry doesnt get changed in the future.  This is mostly relevant
  //  for saving the applications state, you dont want the user to change the
  //  individual components, thus corrupting the state
  bool writeprotected;
};//struct ShieldingSourceModel


//UserState - struct to roughly handle saving the users state to the database,
//  and allow restoring the state.
//  Design desision: all variables, besides writeprotected are left to be public
//  variables since nothing significant would be gained by making them private,
//  and doing so would add a lot of boilerplate code.
struct UserState
{
  enum UserStateType
  {
    kEndOfSessionTemp,  //May be overwritten upon next end of session save
    kUserSaved, //the default state when this snapshot is saved
    kForTest,
    // should add a kPeriodicAutoSaveState
    kUndefinedStateType,
    kEndOfSessionHEAD  //Special case where the end of session should not be deleted on cleanup (ie if it is assocated with a snapshot, rather than a temporary kEndOfSession which gets cleaned up)
  };//enum UserStateType
  
  enum CurrentTab
  {
    kPeakInfo,
    kGammaLines,
    kCalibration,
    kIsotopeSearch,
    kFileTab,
#if( USE_TERMINAL_WIDGET )
    kTerminalTab,
#endif
#if( USE_REL_ACT_TOOL )
    kRelActManualTab,
#endif
    kNoTabs
  };//enum CurrentTab

  enum FeatureMarkers
  {
    kEscapePeaks = 0x1,
    kCompPeak    = 0x2,
    kComptonEdge = 0x4,
    kSumPeak     = 0x8
  };//enum FeatureMarkers
  
  enum ShowingWindows
  {
    kEnergyCalibration = 0x1,
    kDrfSelectSelect = 0x2
    //etc..
  };//enum ShowingWindows
  
  enum SpectrumSubtractMode
  {
    kNoSpectrumSubtract,
    kBackgorundSubtract,
    kContinuumSubtract
  };//enum SpectrumSubtractMode
  
  enum ShownDisplayFeatures
  {
    kDockedWindows          = 0x1,
    kLogSpectrumCounts      = 0x2,
    kVerticalGridLines      = 0x4,
    kHorizontalGridLines    = 0x8,
    kSpectrumLegend         = 0x10,
    kTimeSeriesLegend       = 0x20,
    kShowingShieldSourceFit = 0x40,
    kShowingEnergySearch    = 0x80
#if( USE_TERMINAL_WIDGET )
  , kShowingTerminalWidget  = 0x0100
#endif
#if( USE_REL_ACT_TOOL )
  , kShowingRelActManual    = 0x0200
  , kShowingRelActAuto      = 0x0400
#endif
  , kShowingMultimedia      = 0x000800
  , kShowingGammaXsTool     = 0x001000
  , kShowingDoseCalcTool    = 0x002000
  , kShowing1OverR2Tool     = 0x004000
  , kShowingUnitConvertTool = 0x008000
  , kShowingNucDecayInfo    = 0x010000
  , kShowingEnergyRangeSum  = 0x020000
  , kShowingFluxTool        = 0x040000
#if( USE_DETECTION_LIMIT_TOOL )
  , kShowingDetectionSens   = 0x080000
  , kShowingSimpleMda       = 0x100000
#endif //USE_DETECTION_LIMIT_TOOL
  };//enum ShownDisplayFeatures
  
  //UserState(): default constructor, initializes values to reasonable defaults
  UserState();
  
  Wt::Dbo::ptr<InterSpecUser> user;
    
  Wt::Dbo::ptr<UserState> snapshotTagParent;
  Wt::Dbo::collection< Wt::Dbo::ptr<UserState> > snapshotTags;
    
  //isWriteProtected(): returns if this UserState has been marked read only
  bool isWriteProtected() const;
  
  //makeWriteProtected(...): sets the writeprotected flag so object wont
  //  be re-saved to the database in an altered state.  Should only be called
  //  from within an active transaction.
  //  Also makes the foreground/background/second spectrum write protected, as
  //  well as the source/shielding model.
  static void makeWriteProtected( Wt::Dbo::ptr<UserState> ptr,
                                  Wt::Dbo::Session *session );

  //removeWriteProtection(...): removes write protection flag.
  static void removeWriteProtection( Wt::Dbo::ptr<UserState> ptr,
                                     Wt::Dbo::Session *session );
  
  UserStateType stateType;
  
  Wt::WDateTime creationTime;
  Wt::WDateTime serializeTime;
  
  Wt::WString name;
  Wt::WString description;
  
  //The following indices should probably be converted over to being proper
  //  pointers/collections (or maybye Dbo::weak_ptr), so as to not run into
  //  issues with Wt::Dbo optimizations (e.g. indices not normally assigned
  //  until last recursive Transaction is committed).
  int foregroundId;
  int backgroundId;
  int secondForegroundId;
  int shieldSourceModelId;
  std::string otherSpectraCsvIds;
  std::string foregroundSampleNumsCsvIds;
  std::string secondForegroundSampleNumsCsvIds;
  std::string backgroundSampleNumsCsvIds;
  std::string showingDetectorNumbersCsv;
  
  double energyAxisMinimum, energyAxisMaximum;
  double countsAxisMinimum, countsAxisMaximum;
  
  // TODO: should add time chart limits here - if showing
  
  /** deprecated */
  int displayBinFactor;
  
  int shownDisplayFeatures;  //bitwise or of ShownDisplayFeatures
  
  SpectrumSubtractMode backgroundSubMode;
  CurrentTab currentTab;
  std::string gammaLinesXml;
  std::string isotopeSearchEnergiesXml;
  
  std::string gammaXsToolUri;
  std::string doseCalcToolUri;
  std::string oneOverR2ToolUri;
  std::string unitsConverterToolUri;
  std::string nucDecayInfoUri;
  std::string energyRangeSumUri;
  std::string fluxToolUri;
  std::string detectionSensitivityToolUri;
  std::string simpleMdaUri;
  
  int showingMarkers;        //bitwise or of FeatureMarkers (not implemented yet)
  int disabledNotifications; //bitwise or of (0x1<<WarningWidget::WarningMsgLevel)
  int showingPeakLabels;     //bitwise or of PeakLabels
  int showingWindows;        //(not implemented yet)
  
  std::string userOptionsJson;

  // Note: we are including the whole color theme definition here, and not a link to the theme in
  //   the database.  
  std::string colorThemeJson;
  
  template<class Action>
  void persist( Action &a )
  {
    Wt::Dbo::belongsTo( a, user, "InterSpecUser", Wt::Dbo::OnDeleteCascade );
      
    Wt::Dbo::belongsTo( a, snapshotTagParent, "SnapshotTagParent", Wt::Dbo::OnDeleteCascade );
    Wt::Dbo::hasMany( a, snapshotTags, Wt::Dbo::ManyToOne, "SnapshotTagParent" );
      
      
    Wt::Dbo::field( a, stateType, "StateType" );
    Wt::Dbo::field( a, name, "Name", 127 );
    Wt::Dbo::field( a, description, "Description", 255 );
    Wt::Dbo::field( a, writeprotected, "WriteProtected" );
    Wt::Dbo::field( a, creationTime, "CreationTime" );
    Wt::Dbo::field( a, serializeTime, "SerializeTime" );
    Wt::Dbo::field( a, foregroundId, "ForegroundId" );
    Wt::Dbo::field( a, backgroundId, "BackgroundId" );
    Wt::Dbo::field( a, secondForegroundId, "SecondForegroundId" );
    Wt::Dbo::field( a, shieldSourceModelId, "ShieldSourceModelId" );
    Wt::Dbo::field( a, otherSpectraCsvIds, "OtherSpectraCsvIds", 255 );
    
    Wt::Dbo::field( a, foregroundSampleNumsCsvIds,
                                      "ForegroundSampleNumsCsvIds", 127 );
    Wt::Dbo::field( a, secondForegroundSampleNumsCsvIds,
                                      "SecondForegroundSampleNumsCsvIds", 127 );
    Wt::Dbo::field( a, backgroundSampleNumsCsvIds,
                                      "BackgroundSampleNumsCsvIds", 127 );
    Wt::Dbo::field( a, showingDetectorNumbersCsv,
                                      "ShowingDetectorNumbersCsv", 127 );
    
    Wt::Dbo::field( a, energyAxisMinimum, "EnergyAxisMinimum" );
    Wt::Dbo::field( a, energyAxisMaximum, "EnergyAxisMaximum" );
    Wt::Dbo::field( a, countsAxisMinimum, "CountsAxisMinimum" );
    Wt::Dbo::field( a, countsAxisMaximum, "CountsAxisMaximum" );
    Wt::Dbo::field( a, displayBinFactor,  "DisplayBinFactor" );
    
    Wt::Dbo::field( a, shownDisplayFeatures, "ShownDisplayFeatures" );
    Wt::Dbo::field( a, backgroundSubMode, "BackgroundSubMode" );
    Wt::Dbo::field( a, currentTab, "CurrentTab" );
    Wt::Dbo::field( a, gammaLinesXml, "GammaLinesXml" );
    Wt::Dbo::field( a, isotopeSearchEnergiesXml, "IsotopeSearchEnergiesXml" );
    Wt::Dbo::field( a, showingMarkers, "ShowingMarkers" );
    Wt::Dbo::field( a, disabledNotifications, "DisabledNotifications" );
    Wt::Dbo::field( a, showingPeakLabels, "ShowingPeakLabels" );
    Wt::Dbo::field( a, showingWindows, "ShowingWindows" );
    Wt::Dbo::field( a, userOptionsJson, "UserOptions" );
    Wt::Dbo::field( a, colorThemeJson, "ColorThemeJson" );
    
    Wt::Dbo::field( a, gammaXsToolUri, "GammaXsToolUri" );
    Wt::Dbo::field( a, doseCalcToolUri, "DoseCalcToolUri" );
    Wt::Dbo::field( a, oneOverR2ToolUri, "OneOverR2ToolUri" );
    Wt::Dbo::field( a, unitsConverterToolUri, "UnitsConverterToolUri" );
    Wt::Dbo::field( a, nucDecayInfoUri, "NucDecayInfoUri" );
    Wt::Dbo::field( a, energyRangeSumUri, "EnergyRangeSumUri" );
    Wt::Dbo::field( a, fluxToolUri, "FluxToolUri" );
    Wt::Dbo::field( a, detectionSensitivityToolUri, "DetectionSensitivityToolUri" );
    Wt::Dbo::field( a, simpleMdaUri, "SimpleMdaUri" );
  }//void persist( Action &a )

private:
  //writeprotected: not yet implemented, but intended to allow a way to ensure
  //  database entry doesnt get changed in the future.  This is mostly relevant
  //  for saving the applications state, you dont want the user to change the
  //  individual components, thus corrupting the state
  bool writeprotected;
  
  static void setWriteProtection( Wt::Dbo::ptr<UserState> ptr,
                                  Wt::Dbo::Session *session,
                                  bool protect );
};//struct UserState



struct UseDrfPref
{
  enum UseDrfType
  {
    UseDetectorSerialNumber = 0,
    UseDetectorModelName = 1
  };

  Wt::Dbo::ptr<InterSpecUser> m_user;
  UseDrfType m_field;
  std::string m_criteria;
  int m_flags;
  long long m_drfIndex;
  
  template<class Action>
  void persist( Action &a )
  {
    Wt::Dbo::belongsTo( a, m_user, "InterSpecUser", Wt::Dbo::OnDeleteCascade );
    Wt::Dbo::field( a, m_field, "MatchField" );
    Wt::Dbo::field( a, m_flags, "Flags" );
    Wt::Dbo::field( a, m_drfIndex, "DrfIndex" );
    Wt::Dbo::field( a, m_criteria, "Criteria" );
  }
};//struct UseDrfPref


//Implementation of inline and templated functions
#include <sstream>
#include <stdexcept>
#include <Wt/Dbo/Session>
#include <Wt/Dbo/Transaction>


template<typename T>
T InterSpecUser::preferenceValue( const std::string &name,
                                  InterSpec *viewer )
{
  boost::any value = InterSpecUser::preferenceValueAny( name, viewer );
  return boost::any_cast<T>( value );
}//preferenceValue(...)


template<typename T>
T InterSpecUser::preferenceValue( const std::string &name ) const
{
  return boost::any_cast<T>( preferenceValueAny(name) );
}//preferenceValue(...)


template<typename T>
void InterSpecUser::setPreferenceValueWorker( Wt::Dbo::ptr<InterSpecUser> user,
                                       const std::string &name, const T &value,
                                       InterSpec *viewer )
{
  //std::cout << "setPreferenceValue: " << name << std::endl;
  
  if( name.size() > UserOption::sm_max_name_str_len )
    throw std::runtime_error( "Invalid name for preference: " + name );
  
  std::shared_ptr<DataBaseUtils::DbSession> sql = sqlFromViewer( viewer );
  
  const T oldVal = preferenceValue<T>( name, viewer );
  if( value == oldVal )
    return;
  
  if( sql->session() != user.session() )
    throw std::runtime_error( "InterSpecUser::setPreferenceValue() must be called"
                             " with same db session as user" );
  
  std::vector< Wt::Dbo::ptr<UserOption> > options;
  {//begin interacting with the database
    DataBaseUtils::DbTransaction transaction( *sql );
    Wt::Dbo::collection< Wt::Dbo::ptr<UserOption> > optioncol
                                    = user->m_dbPreferences.find()
                                         .where( "name=?" ).bind( name );
    std::copy( optioncol.begin(), optioncol.end(), std::back_inserter(options) );
    transaction.commit();
  }//end interacting with the database
  
  Wt::Dbo::ptr<UserOption> option;
  const size_t noptions = options.size();
  if( noptions == 0 )
  {
    UserOption *newoption = new UserOption();
    newoption->m_name = name;
    newoption->m_user = user;
    
    if( typeid(value) == typeid(std::string) )
      newoption->m_type = UserOption::String;
    else if( typeid(value) == typeid(double) || typeid(value) == typeid(float) )
      newoption->m_type = UserOption::Decimal;
    else if( typeid(value) == typeid(int)
            || typeid(value) == typeid(unsigned int)
            || typeid(value) == typeid(long long) )
      newoption->m_type = UserOption::Integer;
    else if( typeid(value) == typeid(bool) )
      newoption->m_type = UserOption::Boolean;
    else
    {
      delete newoption;
      throw std::runtime_error( "setPreferenceValue(...): invalid type: "
                               + std::string(typeid(value).name()) );
    }
    
    DataBaseUtils::DbTransaction transaction( *sql );
    option = sql->session()->add( newoption );
    transaction.commit();
  }else if( noptions == 1 )
  {
    option = options.front();
  }else
  {
    //Hmmm, not sure how this happened, but I have made it here for "ColorThemeIndex".
    //  We will take the first option and
#if( PERFORM_DEVELOPER_CHECKS )
    char buffer[1024];
    snprintf( buffer, sizeof(buffer), "Invalid number of preferences (%i) for %s for user %s; will"
                                      " remove all of them after the first from the database.",
              static_cast<int>(noptions), name.c_str(), user->userName().c_str() );
    log_developer_error( __func__, buffer );
#endif
    
    option = options.front();
    for( size_t i = 1; i < noptions; ++i )
    {
      try
      {
        DataBaseUtils::DbTransaction transaction( *sql );
        options[i].remove();
      }catch(...)
      {
        std::cerr << "Caught exception removing duplicate preference from database" << std::endl;
      }
    }//for( size_t i = 1; i < noptions; ++i )
    
    //throw runtime_error( "Invalid number of preferences for " + name + " for user " + user->userName() );
  }//if( noptions == 0 ) / else / else
  
  std::stringstream valuestrm;
  valuestrm << value;
  std::string strval = valuestrm.str();
  if( strval.size() > UserOption::sm_max_value_str_len )
    strval = strval.substr( 0, UserOption::sm_max_value_str_len );
  
  {//begin interacting with the database
    DataBaseUtils::DbTransaction transaction( *sql );
    option.modify()->m_value = strval;
    transaction.commit();
//    std::cerr << "Commited " << option->m_name << " value of "
//              << option->m_value << " to user " << option->m_user.id()
//              << " with ID=" << option.id() << std::endl;
  }//end interacting with the database
  
  user->m_preferences[name] = option->value();

}//setPreferenceValueWorker(...)




template<typename T>
void InterSpecUser::setPreferenceValue( Wt::Dbo::ptr<InterSpecUser> user,
                                       const std::string &name, const T &value,
                                       InterSpec *viewer )
{
  setPreferenceValueWorker( user, name, value, viewer );
}


template<class T>
void InterSpecUser::addCallbackWhenChanged( Wt::Dbo::ptr<InterSpecUser> &user,
                                           const std::string &name, const T &fcn )
{
  assert( user );
  if( !user ) //shouldnt ever happen
    return;
  
  //Make sure a valid bool preference
  user->preferenceValue<bool>( name );
  
  std::shared_ptr<Wt::Signals::signal<void(bool)>> &signal = user->m_onBoolChangeSignals[name];
  if( !signal )
    signal = std::make_shared<Wt::Signals::signal<void(bool)>>();
  signal->connect( fcn );
}//addCallbackWhenChanged(...)


template<class T, class V>
void InterSpecUser::addCallbackWhenChanged( Wt::Dbo::ptr<InterSpecUser> &user,
                                            const std::string &name,
                                            T *target, void (V::*method)(bool) )
{
  assert( user );
  if( !user ) //shouldnt ever happen
    return;
  
  //Make sure a valid bool preference
  user->preferenceValue<bool>( name );
  
  // Retrieve (or create) the signal, and connect things up
  std::shared_ptr<Wt::Signals::signal<void(bool)>> &signal = user->m_onBoolChangeSignals[name];
  if( !signal )
    signal = std::make_shared<Wt::Signals::signal<void(bool)>>();
  
  signal->connect( boost::bind(method, target, boost::placeholders::_1) );
}//addCallbackWhenChanged(...)

#endif //InterSpecUser_h


