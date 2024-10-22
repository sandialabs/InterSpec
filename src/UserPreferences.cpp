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
#include <vector>
#include <functional>

#include <Wt/WObject>
#include <Wt/WSignal>
#include <Wt/Dbo/ptr>
#include <Wt/WCheckBox>

//#include <boost/scope_exit.hpp> //temporary for debug

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/UserPreferences.h"

using namespace Wt;
using namespace std;

namespace
{
  /** Throws exception if the value string doesnt represent the specified data type. */
  void check_string_gives_type( const string &val_str, const UserOption::DataType type, const string &name )
  {
    switch( type )
    {
      case UserOption::String:
        break;
        
      case UserOption::Decimal:
      {
        double val;
        if( !SpecUtils::parse_double(val_str.c_str(), val_str.size(), val) )
          throw std::logic_error( "The decimal preference '" + name
                                 + "' is invalid type; given value: '" + val_str + "'." );
        break;
      }
        
      case UserOption::Integer:
      {
        int val;
        if( !SpecUtils::parse_int(val_str.c_str(), val_str.size(), val) )
          throw std::logic_error( "The integer preference '" + name
                                 + "' is invalid type; given value: '" + val_str + "'." );
        break;
      }
        
      case UserOption::Boolean:
      {
        if( (val_str != "0") && (val_str != "1") )
          throw std::logic_error( "The boolean preference '" + name
                                 + "' must be '0' or '1'; given value: '" + val_str + "'." );
        break;
      }
    }//switch( ptr->m_type )
  };//check_string_gives_type lambda
  
  
  Wt::Dbo::ptr<UserOption> parseUserOption( const rapidxml::xml_node<char> *pref )
  {
    using rapidxml::internal::compare;
    typedef rapidxml::xml_attribute<char> XmlAttribute;
    
    UserOption *option = new UserOption();
    Wt::Dbo::ptr<UserOption> option_owner( option );
    try
    {
      const XmlAttribute *name_att = pref->first_attribute( "name", 4 );
      const XmlAttribute *type_att = pref->first_attribute( "type", 4 );
      if( !name_att || !name_att->value() || !type_att || !type_att->value()  )
        throw runtime_error( "Ill formatted default preferences file" );
      
      const char *typestr = type_att->value();
      std::size_t typestrlen = type_att->value_size();
      
      option->m_name = name_att->value();
      
      // Get rid of "_tablet" or "_phone"
      auto pos = option->m_name.find( "_tablet" );
      if( pos == string::npos )
        pos = option->m_name.find( "_phone" );
      if( pos != string::npos )
        option->m_name = option->m_name.substr(0, pos);
      
      const char *valuestr = pref->value();
      
      // We'll let "String" type values be empty, but no other types.
      if( !valuestr && !compare(typestr, typestrlen, "String", 6, true) )
        throw runtime_error( "Missing value for \"" + option->m_name + "\" in default preferences file" );
      
      if( valuestr )
        option->m_value = std::string( valuestr, valuestr + pref->value_size() );
      
      if( option->m_name.length() > UserOption::sm_max_name_str_len )
        throw runtime_error( "Pref \"" + option->m_name + "\" name to long" );
      if( option->m_value.length() > UserOption::sm_max_value_str_len )
        throw runtime_error( "Pref \"" + option->m_name + "\" value to long" );
      
      if( compare( typestr, typestrlen, "String", 6, true) )
        option->m_type = UserOption::String;
      else if( compare( typestr, typestrlen, "Decimal", 7, true) )
        option->m_type = UserOption::Decimal;
      else if( compare( typestr, typestrlen, "Integer", 7, true) )
        option->m_type = UserOption::Integer;
      else if( compare( typestr, typestrlen, "Boolean", 7, true) )
      {
        option->m_type = UserOption::Boolean;
        if( option->m_value == "true" )
          option->m_value = "1";
        else if( option->m_value == "false" )
          option->m_value = "0";
      }else
        throw runtime_error( "Invalid \"type\" (\"" + string(typestr) + "\") "
                            " for pref " + string(typestr) );
      
      check_string_gives_type( option->m_value, option->m_type, option->m_name );
    }catch( std::exception &e )
    {
      throw runtime_error( e.what() );
    }
    
    return option_owner;
  }//UserOption *parseUserOption( rapidxml::xml_node<char> *node )
  
  
  /** Returns a map from preference name, to default values, for all preferences.

   The map key will have "_phone" and "_tablet" appended to the preference names, for specializations
   for those devices, but the name in the `UserOption` will not have those postfixes.  The `Dbo::ptr`
   will NOT have been added to the database or anything.
   */
  map<string,Dbo::ptr<UserOption>> defaultUserPreferences()
  {
    map<string,Dbo::ptr<UserOption>> answer;
    
    const string filename = SpecUtils::append_path( InterSpec::staticDataDirectory(), InterSpecUser::sm_defaultPreferenceFile );
      
    std::vector<char> data;
    SpecUtils::load_file_data( filename.c_str(), data );
      
    rapidxml::xml_document<char> doc;
    const int flags = (rapidxml::parse_normalize_whitespace | rapidxml::parse_trim_whitespace);
      
    doc.parse<flags>( &data.front() );
    rapidxml::xml_node<char> *node = doc.first_node();
    if( !node || !node->name()
         || !rapidxml::internal::compare( node->name(), node->name_size(), "preferences", 11, true) )
      throw runtime_error( "defaultUserPreferences: invalid first node" );
      
    for( const rapidxml::xml_node<char> *pref = node->first_node( "pref", 4 );
        pref;
        pref = pref->next_sibling( "pref", 4 ) )
    {
      const rapidxml::xml_attribute<char> *name_att = pref->first_attribute( "name", 4 );
      assert( name_att && name_att->value() && name_att->value_size() );
      if( !name_att || !name_att->value() || !name_att->value_size() )
        continue;
        
      const char *prefname = name_att->value();
      const size_t prefnamelen = name_att->value_size();
      std::string prefnamestr(prefname, prefname + prefnamelen);
        
      Wt::Dbo::ptr<UserOption> option = parseUserOption( pref );
      
      // option->m_name will have had "_phone" and "_tablet" removed, even though `prefnamestr`
      //  will have these in them
      assert( option.get() );
      assert( !SpecUtils::icontains( option->m_name, "_phone") );
      assert( !SpecUtils::icontains( option->m_name, "_tablet") );
        
      answer[prefnamestr] = option;
    }//for( loop over preferences )
  
    return answer;
  }//map<string,Wt::Dbo::ptr<UserOption>> defaultUserPreferences()


  /** The returned option will NOT be in the database. 
   Returned value will be valid, and will throw exception if invalid preference name.
   */
  Wt::Dbo::ptr<UserOption> getDefaultUserPreference( const std::string &name, const int type )
  {
    const map<string,Wt::Dbo::ptr<UserOption>> prefs = defaultUserPreferences();
    
    const bool isphone = ((type & InterSpecUser::PhoneDevice) != 0);
    const bool istablet = ((type & InterSpecUser::TabletDevice) != 0);
    
    if( isphone )
    {
      const auto pos = prefs.find( name + "_phone" );
      assert( (pos == end(prefs)) || pos->second );
      if( pos != end(prefs) )
        return pos->second;
    }//
    
    if( istablet )
    {
      const auto pos = prefs.find( name + "_tablet" );
      assert( (pos == end(prefs)) || pos->second );
      if( pos != end(prefs) )
        return pos->second;
    }//if( istablet )
    
    
    const auto pos = prefs.find( name );
    assert( (pos == end(prefs)) || pos->second );
    if( pos != end(prefs) )
      return pos->second;
    
    //Note: the string "couldn't find preference by name" is currently used in
    //      restoreUserPrefsFromXml(...) to check if a preference with this name is no longer used.
    //      So dont change next string without also changing there - or really, a more robust
    //      indication should be used.
    throw runtime_error( "getDefaultUserPreference(...):"
                        " couldn't find preference by name " + name );
    
    return Wt::Dbo::ptr<UserOption>{};
  }//Wt::Dbo::ptr<UserOption> getDefaultUserPreference( const std::string &name, const int type )
}//namespace


UserPreferences::UserPreferences( InterSpec *parent )
 : Wt::WObject( parent ),
  m_interspec( parent ),
  m_options{}
{
  assert( m_interspec );
  if( !m_interspec )
    throw runtime_error( "UserPreferences got nullptr" );
  
  const Wt::Dbo::ptr<InterSpecUser> &user = m_interspec->user();
  std::shared_ptr<DataBaseUtils::DbSession> sql = m_interspec->sql();
  assert( sql && user );
  if( !sql || !user )
    throw runtime_error( "UserPreferences invalid sql" );
  
  //Takes about 200 microseconds to get all the preferences from the database, if the user
  //  already has them all
  //Takes about 1 milli-second to get default preference values and put them into the database
  //  if the havent already been put in.
  //const auto start_time = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::now() );
  //BOOST_SCOPE_EXIT(start_time){
  //  const auto end_time = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::now() );
  //  cout << "Took " << (end_time - start_time).count() << " microseconds to get initial preferences from "
  //  << "data-base" << "." << endl;
  //} BOOST_SCOPE_EXIT_END
  
  // Since we only create this class in one place, where there is already
  //  an active transaction, we dont need to do out own, but the effects
  //  of the transaction are recursive, so I dont think it adds any real
  //  extra overhead to just be safe
  DataBaseUtils::DbTransaction transaction( *sql );
  
  const Dbo::collection<Dbo::ptr<UserOption>> &dbprefs = user->preferences();
  for( Dbo::collection<Dbo::ptr<UserOption>>::const_iterator iter = dbprefs.begin();
      iter != dbprefs.end();
      ++iter )
  {
    Dbo::ptr<UserOption> option = *iter;
    assert( option );
    
    const string &name = option->m_name;
    
    if( m_options.count(name) )
    {
      cerr << "Found a duplicate preference for '" << name << "', will delete the former one" << endl;
      assert( 0 );
      m_options[name].remove();
    }
    
    try
    {
      check_string_gives_type( option->m_value, option->m_type, option->m_name );
    }catch( std::exception &e )
    {
      cerr << "User preference in database gave error: " << e.what() << endl;
      assert( 0 );
      
      try
      {
        option = getDefaultUserPreference( option->m_name,  user->deviceType() );
      }catch( std::exception & )
      {
        cerr << "Failed to get a default preference value for '" << option->m_name
        << "' - skipping." << endl;
        assert( 0 );
        continue;
      }
    }//try / catch - to make sure database value string represents purported data type
    
    m_options[name] = option;
  }//for( loop over preferences in the database );
  
  if( m_options.empty() )
  {
    auto app = dynamic_cast<InterSpecApp *>( wApp );
    const bool isphone = app && app->isPhone();
    const bool istablet = app && app->isTablet();
    
    const map<string,Dbo::ptr<UserOption>> defopts = defaultUserPreferences();
    
    for( const auto &name_opt : defopts )
    {
      const string &name = name_opt.first;
      Dbo::ptr<UserOption> opt = name_opt.second;
      assert( opt );
      if( !opt )
        continue;
      
      if( !istablet && (name.find("_tablet") != string::npos) )
        continue;
      
      if( !isphone && (name.find("_phone") != string::npos) )
        continue;
      
      if( m_options.count(opt->m_name) )
        continue;
      
      opt.modify()->m_user = user;
      sql->session()->add( opt );
      m_options[opt->m_name] = opt;
    }//for( const auto &name_opt : defopts )
  }//if( m_options.empty() )
  
  transaction.commit();
}//UserPreferences constructor


void UserPreferences::setBoolPreferenceValue( const std::string &name,
                                       const bool &value,
                                       InterSpec *viewer )
{
  setPreferenceValueInternal( name, value, viewer );
}


void UserPreferences::setPreferenceValueInternal( const std::string &name,
                                   const bool &value,
                                   InterSpec *viewer )
{
  const string value_as_str = value ? "1" : "0";
  const bool updated = UserPreferences::setPreferenceValueWorker( name, value_as_str, viewer );
  
  // If preference value was not updated, we wont call the callbacks - if we do, we may get stuck
  //  in an infinite recursion loop.
  if( !updated )
    return;
  
  UserPreferences * const self = viewer->preferences();
  assert( self );
  
  const auto callback_pos = self->m_onBoolChangeSignals.find(name);
  if( (callback_pos != end(self->m_onBoolChangeSignals)) && callback_pos->second )
  {
    Wt::Signals::signal<void(bool)> &signal = *callback_pos->second;
    signal(value);
  }
}//setPreferenceValueInternal(...)


void UserPreferences::setPreferenceValueInternal( const std::string &name,
                               const int &value,
                               InterSpec *viewer )
{
  const string value_as_str = std::to_string(value);
  UserPreferences::setPreferenceValueWorker( name, value_as_str, viewer );
}//void UserPreferences::setPreferenceValueInternal(int)


void UserPreferences::setPreferenceValueInternal( const std::string &name,
                               const double &value,
                               InterSpec *viewer )
{
  const std::string value_as_str = SpecUtils::printCompact( value, 12 );
  UserPreferences::setPreferenceValueWorker( name, value_as_str, viewer );
}//void UserPreferences::setPreferenceValueInternal(double)


void UserPreferences::setPreferenceValueInternal( const std::string &name,
                               const std::string &value,
                               InterSpec *viewer )
{
  UserPreferences::setPreferenceValueWorker( name, value, viewer );
}//void UserPreferences::setPreferenceValueInternal(string)




bool UserPreferences::setPreferenceValueWorker( const std::string &name,
                                               std::string value_as_string,
                                               InterSpec *viewer )
{
  if( name.size() > UserOption::sm_max_name_str_len )
    throw std::runtime_error( "Invalid name for preference: " + name );
  
  assert( viewer );
  if( !viewer )
    throw std::runtime_error( "setPreferenceValueWorker: called with invalid viewer ptr" );
  
  
  UserPreferences * const self = viewer->preferences();
  assert( self );
  auto old_pos = self->m_options.find(name);
  if( (old_pos != end(self->m_options)) && (old_pos->second->m_value == value_as_string) )
    return false;
  
  std::shared_ptr<DataBaseUtils::DbSession> sql = viewer->sql();
  const Wt::Dbo::ptr<InterSpecUser> &user = viewer->user();
  
  assert( sql && user );
  if( !sql || !user )
    throw runtime_error( "Invalid sql or user ptr" );
  
  if( sql->session() != user.session() )
    throw std::runtime_error( "UserPreferences::setPreferenceValue() must be called"
                             " with same db session as user" );
  
  if( value_as_string.size() > UserOption::sm_max_value_str_len )
    value_as_string = value_as_string.substr( 0, UserOption::sm_max_value_str_len );
  
  if( old_pos != end(self->m_options) )
  {
    try
    {
      check_string_gives_type( value_as_string, old_pos->second->m_type, name );
      
      DataBaseUtils::DbTransaction transaction( *sql );
      
      old_pos->second.reread(); // In case another session has changed the value, so we dont get a stale object exception
      old_pos->second.modify()->m_value = value_as_string;
      old_pos->second->value();
      
      transaction.commit();
    }catch( std::exception &e )
    {
      assert( 0 );
      std::cerr << "Caught exception setting preference value to database" << std::endl;
    }
    
    return true;
  }//if( old_pos != end(self->m_options) )
  
  
  // If we are here, we dont have a value in memory so we'll check if its been added
  //  to the database since we initialized this session, and if not, we'll start with
  //  default value
  DataBaseUtils::DbTransaction transaction( *sql );
  
  Dbo::collection< Wt::Dbo::ptr<UserOption> > options_in_db
                                    = user->preferences().find()
                                         .where( "name = ?" ).bind( name );
  const size_t noptions = options_in_db.size();
  
  assert( noptions < 2 );
    
  Wt::Dbo::ptr<UserOption> new_opt;
  if( noptions >= 1 )
  {
    vector< Dbo::ptr<UserOption> > options;
    std::copy( options_in_db.begin(), options_in_db.end(), std::back_inserter(options) );
    assert( options.size() == 1 );
    
    new_opt = options.back();
    new_opt.modify()->m_value = value_as_string;
    
    if( noptions > 1 )
    {
      //  We shouldnt actually ever make it here I think.
#if( PERFORM_DEVELOPER_CHECKS )
      char buffer[1024];
      snprintf( buffer, sizeof(buffer), "Invalid number of preferences (%i) for %s for user %s; will"
               " remove all of them after the first from the database.",
               static_cast<int>(noptions), name.c_str(), user->userName().c_str() );
      log_developer_error( __func__, buffer );
#endif
      
      // Remove all but the most recent option (assuming results are returned sorted by ID)
      for( size_t i = 0; (i+1) < options.size(); ++i )
        options[i].remove();
    }//if( noptions > 1 )
  }else
  {
    new_opt = getDefaultUserPreference( name, user->deviceType() );
    new_opt.modify()->m_user = user;
    new_opt.modify()->m_value = value_as_string;
    sql->session()->add( new_opt );
  }//if( optioncol.size() >= 1 ) / else
    
  transaction.commit();
  
  return true;
}//setPreferenceValueWorker(...)


boost::any UserPreferences::preferenceValueAny( const std::string &name, InterSpec *viewer )
{
  // These first few things are just for debug
  //bool from_default = false;
  //const auto start_time = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::now() );
  //BOOST_SCOPE_EXIT(start_time, name, from_default){
  //  const auto end_time = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::now() );
  //  cout << "Took " << (end_time - start_time).count() << " microseconds to get '" << name
  //  << "' preference from " << (from_default ? "default-xml" : "data-base") << "." << endl;
  //} BOOST_SCOPE_EXIT_END
  
  if( !viewer )
  {
    Dbo::ptr<UserOption> option = getDefaultUserPreference( name, InterSpecUser::DeviceType::Desktop );
    boost::any value = option->value();
    return value;
  }
  
  UserPreferences * const self = viewer->preferences();
  assert( self );
  
  auto known_pos = self->m_options.find(name);
  assert( (known_pos == end(self->m_options)) || known_pos->second.get() );
  if( known_pos != end(self->m_options) )
    return known_pos->second->value();
  
  const Dbo::ptr<InterSpecUser> &user = viewer->user();
  std::shared_ptr<DataBaseUtils::DbSession> sql = viewer->sql();
  
  if( !user || !sql )
    throw std::runtime_error( "preferenceValueAny(...): invalid usr or sql ptr" );
  
  //from_default = true;
  
  Wt::Dbo::ptr<UserOption> def_value = getDefaultUserPreference( name, user->deviceType() );
  
  DataBaseUtils::DbTransaction transaction( *sql );
  def_value.modify()->m_user = user;
  sql->session()->add( def_value );
  transaction.commit();
  
  self->m_options[name] = def_value;
   
  return def_value->value();
}//boost::any preferenceValue( const std::string &name, InterSpec *viewer );



void UserPreferences::associateWidget( const std::string &name,
                                     Wt::WCheckBox *cb,
                                     InterSpec *viewer )
{
  const bool value = preferenceValue<bool>( name, viewer );
  cb->setChecked( value );
  
  viewer->preferences()->addCallbackWhenChanged( name, cb, &WCheckBox::setChecked );

  cb->checked().connect( boost::bind( &UserPreferences::setBoolPreferenceValue, name, true, viewer ) );
  cb->unChecked().connect( boost::bind( &UserPreferences::setBoolPreferenceValue, name, false, viewer ) );
  
  /*
   //We need to emit the checked() and unChecked() signals so that any side-effect
   //  can happen from the change.  For actually changing the state of the widget
   //  we can safely do this (incase 'cb' gets deleted at some point) using
   //  WApplication::bind, but to call the emit functions I couldnt think of a
   //  safe way to do this, other than searching the widget tree and making sure
   //  that we can find the widget (hence, know it hasnt been deleted) before
   //  de-referening the 'cb' passed into this function; this is a bit of a hack
   //  but it works for now.
   const string cbid = cb->id();
   
  std::function<void(boost::any)> fcn = [=]( boost::any valueAny ){
    const bool value = boost::any_cast<bool>(valueAny);
    const bool setCbChecked = reverseValue ? !value : value;
    
    // The below doesnt seem to find widgets in AuxWindows (and maybe pop-ups)
    auto w = wApp->domRoot()->findById(cbid);
    if( !w && wApp->domRoot2() )
      w = wApp->domRoot2()->findById(cbid);
    if( !w && wApp->root() )
      w = wApp->root()->findById(cbid);
    
    if( w )
    {
      cb->changed().emit();
      if( value )
        cb->checked().emit();
      else
        cb->unChecked().emit();
    }else
    {
      cerr << "Couldnt find widget with cbid='" << cbid << "', so wont call any side-effect functions for pref '"
           << name << "'" << endl;
    }
  };//fcn
  
   UserPreferences::associateFunction( user, name, fcn, viewer );
   */
}//void associateWidget( )


/*
void UserPreferences::associateWidget( Wt::Dbo::ptr<InterSpecUser> user,
                                          const std::string &name,
                                          Wt::WDoubleSpinBox *sb,
                                         InterSpec *viewer )
{
  const double value = preferenceValue<double>( name, viewer );
  
  sb->setValue( value );
  sb->valueChanged().connect(
                      boost::bind( &UserPreferences::setPreferenceValue<double>,
                                   user, name, boost::placeholders::_1, viewer ) );
  
  const string sbid = sb->id();
  
  std::function<void(boost::any)> fcn = [=]( boost::any valueAny ){
    const double value = boost::any_cast<double>(valueAny);
  
    auto w = wApp->domRoot()->findById(sbid);
    if( !w && wApp->domRoot2() )
      w = wApp->domRoot2()->findById(sbid);
    
    if( w )
    {
      sb->setValue( value );
      sb->valueChanged().emit( value );
    }else
    {
      cerr << "Couldnt find WDoubleSpinBox with id='" << sbid << "', so wont call any side-effect functions" << endl;
    }
  };//fcn
  
 UserPreferences::associateFunction( user, name, fcn, viewer );
}//void associateWidget(...)


void UserPreferences::associateWidget( Wt::Dbo::ptr<InterSpecUser> user,
                                     const std::string &name,
                                     Wt::WSpinBox *sb,
                                    InterSpec *viewer )
{
  const int value = preferenceValue<int>( name, viewer );
  
  sb->setValue( value );
  
  sb->valueChanged().connect(
                             boost::bind( &UserPreferences::setPreferenceValue<int>,
                                         user, name, boost::placeholders::_1, viewer ) );
  
  const string sbid = sb->id();
  
  std::function<void(boost::any)> fcn = [=]( boost::any valueAny ){
    const int value = boost::any_cast<int>(valueAny);
    
    auto w = wApp->domRoot()->findById(sbid);
    if( !w && wApp->domRoot2() )
      w = wApp->domRoot2()->findById(sbid);
    
    if( w )
    {
      sb->setValue( value );
      sb->valueChanged().emit( value );
    }else
    {
      cerr << "Couldnt find WSpinBox with id='" << sbid << "', so wont call any side-effect functions" << endl;
    }
  };//fcn
  
 UserPreferences::associateFunction( user, name, fcn, viewer );
}//void UserPreferences::associateWidget(...)
*/


void UserPreferences::restoreUserPrefsFromXml( const rapidxml::xml_node<char> *prefs_node,
                                InterSpec *viewer )
{
  using namespace rapidxml;
  using rapidxml::internal::compare;
  
  if( !viewer || !prefs_node
      || !compare( prefs_node->name(), prefs_node->name_size(), "preferences", 11, true) )
    throw runtime_error( "restoreUserPrefsFromXml: invalid input" );
  
 
  for( const xml_node<char> *pref = prefs_node->first_node( "pref", 4 ); pref;
       pref = pref->next_sibling( "pref", 4 ) )
  {
    Wt::Dbo::ptr<UserOption> option;
    try
    {
      option = parseUserOption(pref);
      assert( option );
      
      switch( option->m_type )
      {
        case UserOption::String:
        {
          const string value = boost::any_cast<string>( option->value() );
          setPreferenceValueInternal( option->m_name, value, viewer );
          break;
        }//case String
          
        case UserOption::Decimal:
        {
          const double value = boost::any_cast<double>( option->value() );
          setPreferenceValueInternal( option->m_name, value, viewer );
          break;
        }//case Decimal
          
        case UserOption::Integer:
        {
          const int value = boost::any_cast<int>( option->value() );
          setPreferenceValueInternal( option->m_name, value, viewer );
          break;
        }//case Integer
          
        case UserOption::Boolean:
        {
          const bool value = boost::any_cast<bool>( option->value() );
          setPreferenceValueInternal( option->m_name, value, viewer );
          break;
        }//case Boolean
      }//switch( datatype )
    }catch( std::exception &e )
    {
      const string errmsg = e.what();
      if( SpecUtils::icontains( errmsg, "couldn't find preference by name" ) )
      {
        cerr << "Warning: couldnt find a preference named '"
             << (option ? option->m_name : string("N/A")) << "' that is in the"
             << " file being loaded, but apears to no longer be used." << endl;
      }else
      {
        throw runtime_error( "Failed to deserialize user prefernces from XML: " + errmsg );
      }
    }//try / catch
  }//for( loop over prefs )
}//void restoreUserPrefsFromXml(...)


rapidxml::xml_node<char> *UserPreferences::userOptionsToXml(
                                    rapidxml::xml_node<char> *parent_node,
                                              InterSpec *viewer ) const
{
  using namespace ::rapidxml;
  
  if( !parent_node )
    throw runtime_error( "userOptionsToXml: invalid input" );
  
  xml_document<char> *doc = parent_node->document();
  
  xml_node<char> *prefs_node = doc->allocate_node( node_element, "preferences" );
  parent_node->append_node( prefs_node );
  
  vector< Dbo::ptr<UserOption> > options;
  
  {//begin codeblock to retrieve preferences from database
    std::shared_ptr<DataBaseUtils::DbSession> sql = viewer->sql();
    const Wt::Dbo::ptr<InterSpecUser> &user = viewer->user();
    
    DataBaseUtils::DbTransaction transaction( *sql );

    std::copy( user->preferences().begin(), user->preferences().end(),
               std::back_inserter(options) );
    transaction.commit();
  }//end codeblock to retrieve prefernces from database
  
  for( vector< Dbo::ptr<UserOption> >::const_iterator iter = options.begin();
      iter != options.end(); ++iter )
  {
    Dbo::ptr<UserOption> option = *iter;
    
    const string &name  = option->m_name;
    const string &value = option->m_value;
    const char *namestr = doc->allocate_string( name.c_str(), name.size()+1 );
    
    const char *valstr  = 0;
    switch( option->m_type )
    {
      case UserOption::String: case UserOption::Decimal: case UserOption::Integer:
        valstr = doc->allocate_string( value.c_str(), value.size()+1 );
      break;
        
      case UserOption::Boolean:
        valstr = (boost::any_cast<bool>(option->value()) ? "true" : "false");
      break;
    }//switch( m_type )
    
    const char *typestr = 0;
    switch( option->m_type )
    {
      case UserOption::String:  typestr = "String";  break;
      case UserOption::Decimal: typestr = "Decimal"; break;
      case UserOption::Integer: typestr = "Integer"; break;
      case UserOption::Boolean: typestr = "Boolean"; break;
    }//switch( m_type )
    
    xml_node<char> *node = doc->allocate_node( node_element, "pref", valstr );
    xml_attribute<char> *name_att = doc->allocate_attribute( "name", namestr );
    xml_attribute<char> *type_att = doc->allocate_attribute( "type", typestr );
    node->append_attribute( name_att );
    node->append_attribute( type_att );
    prefs_node->append_node( node );
  }//for( loop over DB entries )
  
  return prefs_node;
}//xml_node<char> *userOptionsToXml( xml_node<char> * ) const



  
  
