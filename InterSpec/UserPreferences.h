#ifndef UserPrefernces_h
#define UserPrefernces_h
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
#include <string>
#include <type_traits>

#include <Wt/WObject>
#include <Wt/WSignal>
#include <Wt/Dbo/ptr>

class InterSpec;
class UserOption;
class InterSpecUser;
namespace rapidxml{ template<class Ch> class xml_node; }


/**  This class tracks user preferences that get stored in the database.
 
 It is used to access or change user preferences, or hold callbacks to perform
 an action when a value is changed.
 
 This functionality can not be part of the `InterSpecUser` class, as the object
 pointed to by `InterSpec::m_user` is not stable, and will change when it
 gets refreshed from the database, so adding any extra non-database backed info
 there will get lost randomly.
 
 This class also caches database values to speed things up a bit - its not much
 (~60 ms on first load, on a M1 macbook pro, a few ms on subsequent loads),
 but death by a thousand paper cuts is always an issue.
 */
class UserPreferences : public Wt::WObject
{
public:
  /** You must pass in a valid pointer.
   */
  UserPreferences( InterSpec *parent );
  
  
  /** Returns preference value, for the named preference.
   
   If preference has not been previously set, will return default value as
   well as adding this value to to the users preferences.
   
   If no preference by the passed in name is found, nor one with a default value, a
   `runtime_error` will be thrown.
   
   First checks if the preference is in `UserPreferences::m_options`, and then if not,
   will check the database (incase another session has added it since this session began),
   and then finally gets the default value from the XML file, and both adds it to the database,
   and returns that answer.
  */
  static boost::any preferenceValueAny( const std::string &name, InterSpec *viewer );
  
  /** Convenience function to call for #preferenceValueAny */
  template<typename T>
  static T preferenceValue( const std::string &name,
                            InterSpec *viewer );
  
  /** Sets preference value for named preference to both the in-memory value (stored by `UserPreferences`)
   and the database.
   
   If a preference with the specified name isnt already in memory, and the name is not in
   `data/default_preferences.xml` then will throw an exception.
   
   This function merely calls the various `setPreferenceValueInternal(...)` functions,
   but is set up this way to make sure implicit conversions of types is not performed somewhere,
   leading to wonky results.
   */
  template<typename T>
  static typename std::enable_if<
          std::is_same<T, bool>::value ||
          std::is_same<T, int>::value ||
          std::is_same<T, double>::value ||
          std::is_same<T, std::string>::value
      >::type
  setPreferenceValue( const std::string &name, const T& value, InterSpec *viewer );
  
  
  /** Just a convenience function to call above variant, for binding to signals. */
  static void setBoolPreferenceValue( const std::string &name, const bool &value, InterSpec *viewer );
  
    
  /** Add a function to callback when the preference changes, either through loading  a new state,
   or toggling the preference checkbox or whatever.
   
   This variant of #addCallbackWhenChanged is useful to book a slot (member function) of a
   Wt::WObject up to, so the callback lifetime will be limited to the lifetime of the WObject.
   
   This function must be called from the application thread such that InterSpec::instance() will be
   non-void.
   */
  template<class T, class V>
  void addCallbackWhenChanged( const std::string &name, T *target,
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
  void addCallbackWhenChanged( const std::string &name,
                               const T &fcn );
  
  /** Makes it so if the user changes the value via this GUI element, the value of the preference
   stored in memory and in the database will be correspondingly updated.  Also makes it so if
   a different GUI element updates this preference value, or #setPreferenceValue is called for
   this preference, this widgets state will be correspondingly updated as well.
   
   However, the signals you have hooked up to this widget wont be called when the value changes by
   another widget or #setPreferenceValue, to handle this case, use the #addCallbackWhenChanged
   function to add a callback for whenever the value of the preference changes.
   */
  static void associateWidget( const std::string &name,
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

  //userOptionsToXml(): serializes the current user options to a node called
  //  <preferences> under parent passed in, with each preferences being placed
  //  in a <pref> node under it.  Uses same schema as 'default_preferences.xml'.
  //Returns the <preferences> node.
  rapidxml::xml_node<char> *userOptionsToXml( rapidxml::xml_node<char> *parent,
                                             InterSpec *viewer ) const;
  
  //restoreUserPrefsFromXml(): sets the currently active user preferences from
  //  passed in XML node.  'preferences' must be a node named <preferences>
  //  and should be in same format that userOptionsToXml() creates.
  //  If there is a current user option that is not in the passd in XML, its
  //  value will not be altered.  Any new options in the XML not currently in
  //  memory will be set.
  //Throws exception upon error.
  static void restoreUserPrefsFromXml( const rapidxml::xml_node<char> *preferences,
                                InterSpec *viewer );

protected:
  
  static void setPreferenceValueInternal( const std::string &name,
                                 const std::string &value,
                                 InterSpec *viewer );
  
  /** Similar to `setPreferenceValueInternal(string,string,InterSpec)`, but for boolean
   preferences, and also calls functions in `UserPreferences::m_onBoolChangeSignals`.
   */
  static void setPreferenceValueInternal( const std::string &name,
                                 const bool &value,
                                 InterSpec *viewer );
  
  /** Similar to `setPreferenceValueInternal(string,string,InterSpec)`, but for integer preferences. */
  static void setPreferenceValueInternal( const std::string &name,
                                 const int &value,
                                 InterSpec *viewer );
  
  /** Similar to `setPreferenceValueInternal(string,string,InterSpec)`, but for double preferences. */
  static void setPreferenceValueInternal( const std::string &name,
                                 const double &value,
                                 InterSpec *viewer );
  
  static void setPreferenceValueWorker( const std::string &name,
                                       std::string value_as_string,
                                       InterSpec *viewer );
  
protected:
  InterSpec *m_interspec;
  std::map<std::string,Wt::Dbo::ptr<UserOption>> m_options;
  
  /** Holds callbacks set from #addCallbackWhenChanged, for boolean preferences.
   
   I believe using a Wt::Signals::signal allows makes it so we can connect widget slots (function
   calls), and then the connection will automatically get deleted when the widget is deleted, making
   things safe.
   
   Note that these function may be "safe" for widgets getting deleted (e.g., signals automatically
   disconnected), if the signal is connected to an object deriving from Wt::WObject.
   
   Currently only boolean preferences require callbacks in InterSpec; in the future, if other
   preference types (string, int, doubles) need callbacks, we will have to re-factor things, or
   add in analogous member variables.
   */
  std::map<std::string,std::shared_ptr<Wt::Signals::signal<void(bool)>>> m_onBoolChangeSignals;
};//class UserPreferences



template<typename T>
T UserPreferences::preferenceValue( const std::string &name,
                                  InterSpec *viewer )
{
  boost::any value = UserPreferences::preferenceValueAny( name, viewer );
  return boost::any_cast<T>( value );
}//preferenceValue(...)


template<typename T>
typename std::enable_if<
        std::is_same<T, bool>::value ||
        std::is_same<T, int>::value ||
        std::is_same<T, double>::value ||
        std::is_same<T, std::string>::value
    >::type
UserPreferences::setPreferenceValue(  const std::string &name, const T& value, InterSpec *viewer  )
{
  setPreferenceValueInternal( name, value, viewer );
}

template<class T>
void UserPreferences::addCallbackWhenChanged( const std::string &name, const T &fcn )
{
  //Make sure a valid bool preference
  preferenceValue<bool>( name, m_interspec );
  
  std::shared_ptr<Wt::Signals::signal<void(bool)>> &signal = m_onBoolChangeSignals[name];
  if( !signal )
    signal = std::make_shared<Wt::Signals::signal<void(bool)>>();
  signal->connect( fcn );
}//addCallbackWhenChanged(...)


template<class T, class V>
void UserPreferences::addCallbackWhenChanged( const std::string &name,
                                            T *target, void (V::*method)(bool) )
{
  //Make sure a valid bool preference
  preferenceValue<bool>( name, m_interspec );
  
  // Retrieve (or create) the signal, and connect things up
  std::shared_ptr<Wt::Signals::signal<void(bool)>> &signal = m_onBoolChangeSignals[name];
  if( !signal )
    signal = std::make_shared<Wt::Signals::signal<void(bool)>>();
  
  signal->connect( boost::bind(method, target, boost::placeholders::_1) );
}//addCallbackWhenChanged(...)

#endif //UserPrefernces_h
