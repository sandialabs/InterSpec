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

#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include <boost/iostreams/stream.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/back_inserter.hpp>


#include <Wt/Dbo/Dbo>
#include <Wt/WSpinBox>
#include <Wt/WCheckBox>
#include <Wt/WApplication>
#include <Wt/WRadioButton>
#include <Wt/WDoubleSpinBox>

#if( HAS_WTDBOSQLITE3 )
#include "Wt/Dbo/backend/Sqlite3"
#endif

#if( HAS_WTDBOMYSQL )
#include "Wt/Dbo/backend/MySQL"
#endif

#if( HAS_WTDBOPOSTGRES )
#include "Wt/Dbo/backend/Postgres"
#endif

#if( HAS_WTDBOFIREBIRD )
#include "Wt/Dbo/backend/Firebird"
#endif

#include "3rdparty/date/include/date/date.h"

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/PopupDiv.h"
#include "InterSpec/InterSpec.h"
#include "SpecUtils/Filesystem.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/DetectorPeakResponse.h"

#if( ALLOW_SAVE_TO_DB_COMPRESSION )
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#endif

using namespace Wt;
using namespace std;




namespace Wt
{
  namespace Dbo
  {
    const char *sql_value_traits<FileData_t>
    ::type(SqlConnection *conn, int size)
    {
      return conn->blobType();
    }
    
    void sql_value_traits< FileData_t >
    ::bind(const FileData_t& v, SqlStatement *statement,
           int column, int size)
    {
      statement->bind(column, v);
    }
    
    bool sql_value_traits<FileData_t >
    ::read( FileData_t &v, SqlStatement *statement, int column, int size )
    {
      return statement->getResult(column, &v, size);
    }
  }//namespace Dbo
}//namespace Wt

//Uhhhhg! this whole having to use vector<unsigned char> to write to the
//  database, but char based things is a pain!  So as a horrible work around
//  we will specialize boost::iostream::back_insert_device
namespace boost
{
  namespace iostreams
  {
    template<>
    class back_insert_device < FileData_t >
    {
    public:
      typedef char      char_type;
      typedef sink_tag  category;
      back_insert_device( FileData_t &cnt) : container( &cnt ) { }
      std::streamsize write( const char_type* s, std::streamsize n )
      {
        const FileData_t::value_type *start = (const FileData_t::value_type *)s;
        const FileData_t::value_type *end = start + n;
        container->insert( container->end(), start, end);
        return n;
      }
    protected:
       FileData_t* container;
    };
  }//namespace iostreams
}//namespace boost

const UserFileInDbData::SerializedFileFormat
      UserFileInDbData::sm_defaultSerializationFormat
= UserFileInDbData::k2012N42; //kNativeText;//kNativeBinary; //kNativeText;

const std::string InterSpecUser::sm_defaultPreferenceFile = "default_preferences.xml";

namespace
{
  UserOption *parseUserOption( const rapidxml::xml_node<char> *pref )
  {
    using rapidxml::internal::compare;
    typedef rapidxml::xml_attribute<char> XmlAttribute;
    
    UserOption *option = new UserOption;
    try
    {
      const XmlAttribute *name_att = pref->first_attribute( "name", 4 );
      const XmlAttribute *type_att = pref->first_attribute( "type", 4 );
      if( !name_att || !name_att->value() || !type_att || !type_att->value()  )
        throw runtime_error( "Ill formatted default preferences file" );
      
      const char *typestr = type_att->value();
      std::size_t typestrlen = type_att->value_size();
      
      option->m_name = name_att->value();
      
      const char *valuestr = pref->value();
      
      // We'll let "String" type values be empty, but no other types.
      if( !valuestr )
      {
        if( !compare(typestr, typestrlen, "String", 6, true) )
          throw runtime_error( "Missing value for \"" + option->m_name + "\" in default preferences file" );
        
        valuestr = "";
      }//if( value is empty in XML )
      
      option->m_value = valuestr;
      
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
    }catch( std::exception &e )
    {
      delete option;
      throw runtime_error( e.what() );
    }
    return option;
  }//UserOption *parseUserOption( rapidxml::xml_node<char> *node )
  
  
  /** Compares two boost::any objects to check if their underlying type is
   the same, and if so, if their values are equal.
   
   Currently only supports underlying types of std::string, bool, int, double
   to coorespond to UserOption::DataType.
   
   Throws exception if types are not supported.
   */
  bool boost_any_equal( const boost::any& lhs, const boost::any& rhs )
  {
   
    if( lhs.type() == typeid(std::string) )
    {
      if( rhs.type() != typeid(std::string) )
        return false;
      
      return boost::any_cast<std::string>(lhs) == boost::any_cast<std::string>(rhs);
    }//if( string )
    
    
    if( lhs.type() == typeid(double) || lhs.type() == typeid(float) )
    {
      if( rhs.type() != typeid(double) && rhs.type() != typeid(float) )
        return false;
      
      double lhsval, rhsval;
      if( lhs.type() == typeid(double) )
        lhsval = boost::any_cast<double>(lhs);
      else
        lhsval = boost::any_cast<float>(lhs);
      
      if( rhs.type() == typeid(double) )
        rhsval = boost::any_cast<double>(rhs);
      else
        rhsval = boost::any_cast<float>(rhs);
      
      return lhsval == rhsval;
    }//if( double or float )
    
    if( lhs.type() == typeid(int) || lhs.type() == typeid(unsigned int) || lhs.type() == typeid(long long) )
    {
      if( rhs.type() != typeid(int) && rhs.type() != typeid(unsigned int) && rhs.type() != typeid(long long) )
        return false;
      
      int64_t lhsval, rhsval;
      if( lhs.type() == typeid(int) )
        lhsval = boost::any_cast<int>(lhs);
      else if( lhs.type() == typeid(unsigned int) )
        lhsval = boost::any_cast<unsigned int>(lhs);
      else
        lhsval = boost::any_cast<long long>(lhs);
      
      if( rhs.type() == typeid(int) )
        rhsval = boost::any_cast<int>(rhs);
      else if( rhs.type() == typeid(unsigned int) )
        rhsval = boost::any_cast<unsigned int>(rhs);
      else
        rhsval = boost::any_cast<long long>(rhs);
      
      return lhsval == rhsval;
    }
    
    if( lhs.type() == typeid(bool) )
    {
      if( rhs.type() != typeid(bool) )
        return false;
      
      return boost::any_cast<bool>(lhs) == boost::any_cast<bool>(rhs);
    }//if( string )
    
    cerr << "boost_any_equal: unimplemented type encountered" << endl;
    
    throw std::runtime_error("boost_any_equal: unimplemented type encountered");
  }//boost_any_equal

}//namespace


boost::posix_time::ptime to_ptime( const std::chrono::system_clock::time_point &rhs )
{
  auto dp = ::date::floor<::date::days>(rhs);
  auto ymd = ::date::year_month_day{dp};
  auto time = ::date::make_time(rhs - dp);
  
  boost::posix_time::time_duration td( time.hours().count(),
                                      time.minutes().count(),
                                      time.seconds().count(),
                                      date::floor<chrono::microseconds>(time.subseconds()).count()
                                      );
  
  const unsigned month_num = static_cast<unsigned>(ymd.month());
  boost::gregorian::greg_month month = boost::date_time::months_of_year( month_num );
  /*
   switch( static_cast<unsigned>(ymd.month()) )
  {
    case 1: month = boost::gregorian::Jan; break;
    case 2: month = boost::gregorian::Feb; break;
    case 3: month = boost::gregorian::Mar; break;
    case 4: month = boost::gregorian::Apr; break;
    case 5: month = boost::gregorian::May; break;
    case 6: month = boost::gregorian::Jun; break;
    case 7: month = boost::gregorian::Jul; break;
    case 8: month = boost::gregorian::Aug; break;
    case 9: month = boost::gregorian::Sep; break;
    case 10: month = boost::gregorian::Oct; break;
    case 11: month = boost::gregorian::Nov; break;
    case 12: month = boost::gregorian::Dec; break;
    default: assert( 0 ); throw runtime_error( "wth: invalid month?" ); break;
  }
   */
  
  const boost::gregorian::date d( static_cast<int>(ymd.year()), month, (unsigned)ymd.day() );
  
  return boost::posix_time::ptime( d, td );
}//to_ptime(...)


std::chrono::system_clock::time_point to_time_point( const boost::posix_time::ptime &rhs )
{
  if( rhs.is_special() )
    return {};
  
  unsigned short year = static_cast<unsigned short>( rhs.date().year() );
  const unsigned month = rhs.date().month().as_number();
  const unsigned day = rhs.date().day().as_number();
  
  const int64_t nmicro = rhs.time_of_day().total_microseconds();
  date::year_month_day ymd{ date::year(year), date::month(month), date::day(day) };
  
  date::sys_days days = ymd;
  std::chrono::system_clock::time_point tp = days;
  tp += chrono::microseconds(nmicro);
  
  return tp;
}


void mapDbClasses( Wt::Dbo::Session *session )
{
  session->mapClass<InterSpecUser>( "InterSpecUser" );
  session->mapClass<UserOption>( "UserOption" );
  session->mapClass<UserFileInDb>( "UserFileInDb" );
  session->mapClass<UserFileInDbData>( "UserFileInDbData" );
  session->mapClass<ShieldingSourceModel>( "ShieldingSourceModel" );
  session->mapClass<UserState>( "UserState" );
  session->mapClass<DetectorPeakResponse>( "DetectorPeakResponse");
  session->mapClass<ColorThemeInfo>( "ColorThemeInfo");
  session->mapClass<UseDrfPref>( "UseDrfPref" );
}//void mapDbClasses( Wt::Dbo::Session *session )


UserState::UserState()
  : stateType( kUndefinedStateType ),
    creationTime( Wt::WDateTime::currentDateTime() ),
    serializeTime( Wt::WDateTime::currentDateTime() ),
    foregroundId( -1 ), backgroundId( -1 ), secondForegroundId( -1 ),
    shieldSourceModelId( -1 ),
    energyAxisMinimum( 0.0 ), energyAxisMaximum( 3000.0 ),
    countsAxisMinimum( -1.0 ), countsAxisMaximum( -1.0 ), displayBinFactor( 1 ),
    shownDisplayFeatures( 0 ), backgroundSubMode( kNoSpectrumSubtract ),
    currentTab( kNoTabs ), showingMarkers( 0 ), disabledNotifications( 0 ),
    showingPeakLabels( 0 ), showingWindows( 0 ),
    writeprotected( false )
{
}

void InterSpecUser::setBoolPreferenceValue( Wt::Dbo::ptr<InterSpecUser> user,
                                       const std::string &name,
                                       const bool &value,
                                       InterSpec *viewer )
{
  setPreferenceValue( user, name, value, viewer );
}


void InterSpecUser::setPreferenceValue( Wt::Dbo::ptr<InterSpecUser> user,
                                   const std::string &name,
                                   const bool &value,
                                   InterSpec *viewer )
{
  setPreferenceValueWorker( user, name, value, viewer );
  
  const auto callback_pos = user->m_onBoolChangeSignals.find(name);
  if( (callback_pos != end(user->m_onBoolChangeSignals)) && callback_pos->second )
    (*callback_pos->second)(value);
}//setPreferenceValue(...)


boost::any InterSpecUser::preferenceValueAny( const std::string &name, InterSpec *viewer )
{
  using namespace Wt;
  
  if( !viewer )
  {
    UserOption *option = getDefaultUserPreference( name, DeviceType::Desktop );
    boost::any value = option->value();
    delete option;
    return value;
  }

  //This next line is the only reason InterSpec.h needs to be included
  //  above
  Dbo::ptr<InterSpecUser> &user = userFromViewer(viewer);
  std::shared_ptr<DataBaseUtils::DbSession> sql = sqlFromViewer(viewer);
  
  if( !user || !sql )
    throw std::runtime_error( "preferenceValueAny(...): invalid usr or sql ptr" );
  
  PreferenceMap::const_iterator pos;
  const PreferenceMap &prefs = user->m_preferences;
  
  if( user->m_deviceType & PhoneDevice )
  {
    pos = prefs.find( name + "_phone" );
    if( pos != prefs.end() )
      return pos->second;
  }//if( user->m_deviceType & PhoneDevice )
  
  if( user->m_deviceType & TabletDevice )
  {
    pos = prefs.find( name + "_tablet" );
    if( pos != prefs.end() )
      return pos->second;
  }//if( user->m_deviceType & TabletDevice )
  
  pos = prefs.find( name );
  if( pos != prefs.end() )
    return pos->second;
  
  UserOption *option_raw = getDefaultUserPreference( name, user->m_deviceType );
  Wt::Dbo::ptr<UserOption> option( option_raw );
  option.modify()->m_user = user;
  
  DataBaseUtils::DbTransaction transaction( *sql );
  boost::any value;
  try
  {
    sql->session()->add( option );
    value = option->value();
    user->m_preferences[name] = value;
  }catch( std::exception &e )
  {
    //
    transaction.rollback();
    throw e;
  }//try / catch
  
  transaction.commit();
  
  return value;
}//boost::any preferenceValue( const std::string &name, InterSpec *viewer );


boost::any InterSpecUser::preferenceValueAny( const std::string &name ) const
{
  PreferenceMap::const_iterator pos;
  
  if( m_deviceType & PhoneDevice )
  {
    pos = m_preferences.find( name  + "_phone" );
    if( pos != m_preferences.end() )
      return pos->second;
  }//if( user->m_deviceType & PhoneDevice )
  
  if( m_deviceType & TabletDevice )
  {
    pos = m_preferences.find( name  + "_tablet" );
    if( pos != m_preferences.end() )
      return pos->second;
  }//if( user->m_deviceType & TabletDevice )
  
  pos = m_preferences.find( name );
  if( pos == m_preferences.end() )
    throw std::runtime_error( "No preference value for " + name + ", try"
                             " calling other InterSpecUser::preferenceValue(...) function" );
  return pos->second;
}//boost::any preferenceValueAny( const std::string &name ) const;


void InterSpecUser::associateWidget( Wt::Dbo::ptr<InterSpecUser> user,
                                     const std::string &name,
                                     Wt::WCheckBox *cb,
                                     InterSpec *viewer )
{
  const bool value = preferenceValue<bool>( name, viewer );
  cb->setChecked( value );
  
  InterSpecUser::addCallbackWhenChanged( user, name, cb, &WCheckBox::setChecked );

  cb->checked().connect( boost::bind( &InterSpecUser::setBoolPreferenceValue, user, name, true, viewer ) );
  cb->unChecked().connect( boost::bind( &InterSpecUser::setBoolPreferenceValue, user, name, false, viewer ) );
  
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
  
  InterSpecUser::associateFunction( user, name, fcn, viewer );
   */
}//void associateWidget( )


/*
void InterSpecUser::associateWidget( Wt::Dbo::ptr<InterSpecUser> user,
                                          const std::string &name,
                                          Wt::WDoubleSpinBox *sb,
                                         InterSpec *viewer )
{
  const double value = preferenceValue<double>( name, viewer );
  
  sb->setValue( value );
  sb->valueChanged().connect(
                      boost::bind( &InterSpecUser::setPreferenceValue<double>,
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
  
  InterSpecUser::associateFunction( user, name, fcn, viewer );
}//void associateWidget(...)


void InterSpecUser::associateWidget( Wt::Dbo::ptr<InterSpecUser> user,
                                     const std::string &name,
                                     Wt::WSpinBox *sb,
                                    InterSpec *viewer )
{
  const int value = preferenceValue<int>( name, viewer );
  
  sb->setValue( value );
  
  sb->valueChanged().connect(
                             boost::bind( &InterSpecUser::setPreferenceValue<int>,
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
  
  InterSpecUser::associateFunction( user, name, fcn, viewer );
}//void InterSpecUser::associateWidget(...)
*/


void UserFileInDbData::setFileData( const std::string &path,
                                    const SerializedFileFormat format )
{
  fileData.clear();
#ifdef _WIN32
  const std::wstring wpath = SpecUtils::convert_from_utf8_to_utf16(path);
  std::ifstream file( wpath.c_str() );
#else
  std::ifstream file( path.c_str() );
#endif
  
  if( !file )
    throw runtime_error( "UserFileInDbData::setFileData():"
                        " couldnt open the cached spectrum file." );
  
  const istream::pos_type orig_pos = file.tellg();
  file.seekg( 0, ios::end );
  const istream::pos_type eof_pos = file.tellg();
  file.seekg( orig_pos, ios::beg );
  const size_t filelen = 0 + eof_pos - orig_pos;

#if( ALLOW_SAVE_TO_DB_COMPRESSION )
  {
    namespace io = boost::iostreams;
    
    gzipCompressed = true;
    const size_t pre_mem_size_size
       = static_cast<size_t>( 2.2 * double(UserFileInDb::sm_maxFileSizeBytes) );
    if( filelen > pre_mem_size_size )
      throw FileToLargeForDbException( filelen, pre_mem_size_size );

    io::filtering_ostream compressor;
    compressor.push( io::gzip_compressor() );
    compressor.push( io::back_inserter( fileData ) );
    io::copy( file, compressor );  //default uses 4 kb buffer
    
    if( format == UserFileInDbData::k2012N42 )
      compressor << static_cast<unsigned char>(0);
    
    if( fileData.size() > UserFileInDb::sm_maxFileSizeBytes )
    {
      fileData.clear();
      throw FileToLargeForDbException( fileData.size(),
                                       UserFileInDb::sm_maxFileSizeBytes );
    }//if( file too large )

//    cerr << "UserFileInDbData::setFileData(string): compressed "
//         << filelen << " bytes to " << fileData.size() << " bytes for a savings "
//         << 100.0*((double(filelen)-double(fileData.size()))/double(filelen))
//         << "%" << endl;
  }
#else
  if( filelen > UserFileInDb::sm_maxFileSizeBytes )
    throw runtime_error( "UserFileInDbData::setFileData():"
                        " Spectrum file top large to serialize." );
  gzipCompressed = false;
  fileData.resize( filelen );
  if( !file.read( (char *)&fileData[0], filelen ) )
    throw runtime_error( "UserFileInDbData::setFileData():"
                        " couldnt fully read cached spectrum file." );
  
  if( format == UserFileInDbData::k2012N42 )
    fileData.push_back( static_cast<unsigned char>(0) );
#endif
}//void setFileData( const std::string &path )


UserFileInDb::UserFileInDb()
{
  writeprotected = false;
  isPartOfSaveState = false;
}


bool UserFileInDb::isWriteProtected() const
{
  return writeprotected;
}


void UserFileInDb::makeWriteProtected( Wt::Dbo::ptr<UserFileInDb> ptr )
{
  if( !ptr || ptr->writeprotected )
    return;
  ptr.modify()->writeprotected = true;
}//void UserFileInDb::makeWriteProtected( Wt::Dbo::ptr<UserFileInDb> ptr )


void UserFileInDb::removeWriteProtection( Wt::Dbo::ptr<UserFileInDb> ptr )
{
  if( !ptr || !ptr->writeprotected )
    return;
  ptr.modify()->writeprotected = false;
}//void removeWriteProtection( Wt::Dbo::ptr<UserFileInDb> ptr )


Dbo::ptr<UserFileInDb> UserFileInDb::makeDeepWriteProtectedCopyInDatabase(
                                                  Dbo::ptr<UserFileInDb> orig,
                                                  DataBaseUtils::DbSession &sqldb,
                                                  bool isSaveState )
{
  if( !orig )
    return orig;
  
  auto session = sqldb.session();
  if( !session )
    throw runtime_error( "UserFileInDb::makeDeepWriteProtectedCopyInDatabase():"
                          " invalid input.");
  
  UserFileInDb *newfile = new UserFileInDb();
  newfile->user = orig->user;
  newfile->uuid = orig->uuid;
  newfile->filename = orig->filename;
  
  newfile->userHasModified = orig->userHasModified;
  newfile->uploadTime = orig->uploadTime;
  newfile->serializeTime = WDateTime::currentDateTime();
  newfile->sessionID = orig->sessionID;
  
  newfile->numSamples = orig->numSamples;
  newfile->isPassthrough = orig->isPassthrough;
  newfile->totalLiveTime = orig->totalLiveTime;
  newfile->totalRealTime = orig->totalRealTime;
  newfile->totalGammaCounts = orig->totalGammaCounts;
  newfile->totalNeutronCounts = orig->totalNeutronCounts;
  newfile->numDetectors = orig->numDetectors;
  newfile->hasNeutronDetector = orig->hasNeutronDetector;
  newfile->measurementsStartTime = orig->measurementsStartTime;
  
  newfile->writeprotected = true;
  newfile->isPartOfSaveState = isSaveState;

  newfile->snapshotParent = orig->snapshotParent;
  
  DataBaseUtils::DbTransaction transaction( sqldb );
  
  try
  {
    Dbo::ptr<UserFileInDb> answer = session->add( newfile );
  
    
    for( Dbo::collection< Dbo::ptr<UserFileInDbData> >::const_iterator iter = orig->filedata.begin();
        iter != orig->filedata.end(); ++iter )
    {
      UserFileInDbData *newdata = new UserFileInDbData( **iter );
      newdata->fileInfo = answer;
      session->add( newdata );
    }
      
    for( Dbo::collection< Dbo::ptr<ShieldingSourceModel> >::const_iterator iter
                                                = orig->modelsUsedWith.begin();
        iter != orig->modelsUsedWith.end();
        ++iter )
    {
      if( (*iter)->isWriteProtected() )
      {
        answer.modify()->modelsUsedWith.insert( *iter );
      }else
      {
        ShieldingSourceModel *newmodel = new ShieldingSourceModel();
        newmodel->shallowEquals( **iter );
        Dbo::ptr<ShieldingSourceModel> modelptr = session->add( newmodel );
        modelptr.modify()->filesUsedWith.insert( answer );
        ShieldingSourceModel::makeWriteProtected( modelptr );
      }//if( (*iter)->isWriteProtected() ) / else
    }//for( loop over ShieldingSourceModel )
    
    transaction.commit();
    return answer;
  }catch( std::exception &e )
  {
    transaction.rollback();
    throw runtime_error( e.what() );
  }//try catch
  
  return Dbo::ptr<UserFileInDb>();
}//makeDeepWriteProtectedCopyInDatabase(...)


bool UserState::isWriteProtected() const
{
  return writeprotected;
}//bool isWriteProtected() const


void ShieldingSourceModel::shallowEquals( const ShieldingSourceModel &rhs )
{
  user = rhs.user;
  name = rhs.name;
  description = rhs.description;
  creationTime = rhs.creationTime;
  serializeTime = rhs.serializeTime;
  //filesUsedWith = rhs.filesUsedWith;
  xmlData = rhs.xmlData;
  writeprotected = rhs.writeprotected;
}


void UserState::setWriteProtection( Wt::Dbo::ptr<UserState> ptr,
                                    Wt::Dbo::Session *session,
                                    bool protect )
{
  if( !ptr || (ptr->writeprotected==protect) )
    return;
  
  ptr.modify()->writeprotected = protect;
  
  Dbo::ptr<ShieldingSourceModel> fitmodel;
  Dbo::ptr<UserFileInDb> dbforeground, dbsecond, dbbackground;
  
  if( ptr->foregroundId >= 0 )
    dbforeground = session->find<UserFileInDb>().where( "id = ?" )
    .bind( ptr->foregroundId );
  if( ptr->backgroundId >= 0 )
    dbbackground = session->find<UserFileInDb>().where( "id = ?" )
    .bind( ptr->backgroundId );
  if( ptr->secondForegroundId >= 0 )
    dbsecond = session->find<UserFileInDb>().where( "id = ?" )
    .bind( ptr->secondForegroundId );
  if( ptr->shieldSourceModelId >= 0 )
    fitmodel = session->find<ShieldingSourceModel>()
    .where( "id = ?" ).bind( ptr->shieldSourceModelId );
  
  if( dbforeground && protect )
    UserFileInDb::makeWriteProtected( dbforeground );
  else if( dbforeground )
    UserFileInDb::removeWriteProtection( dbforeground );
  
  if( dbsecond && protect )
    UserFileInDb::makeWriteProtected( dbsecond );
  else if( dbsecond )
    UserFileInDb::removeWriteProtection( dbsecond );
  
  if( dbbackground && protect )
    UserFileInDb::makeWriteProtected( dbbackground );
  else if( dbbackground )
    UserFileInDb::removeWriteProtection( dbbackground );
  
  if( fitmodel && protect )
    ShieldingSourceModel::makeWriteProtected( fitmodel );
  else if( fitmodel )
    ShieldingSourceModel::removeWriteProtection( fitmodel );
}//void setWriteProtection(...)


void UserState::makeWriteProtected( Wt::Dbo::ptr<UserState> ptr,
                                    Wt::Dbo::Session *session )
{
  setWriteProtection( ptr, session, true );
}//void makeWriteProtected( Wt::Dbo::ptr<UserFileInDb> ptr )


void UserState::removeWriteProtection( Wt::Dbo::ptr<UserState> ptr,
                                       Wt::Dbo::Session *session )
{
  setWriteProtection( ptr, session, false );
}//void removeWriteProtection( Wt::Dbo::ptr<UserFileInDb> ptr )


bool ShieldingSourceModel::isWriteProtected() const
{
  return writeprotected;
}

void ShieldingSourceModel::makeWriteProtected(
                                        Wt::Dbo::ptr<ShieldingSourceModel> ptr )
{
  if( !ptr || ptr->writeprotected )
    return;
  ptr.modify()->writeprotected = true;  
}


void ShieldingSourceModel::removeWriteProtection( Wt::Dbo::ptr<ShieldingSourceModel> ptr )
{
  if( !ptr || !ptr->writeprotected )
    return;
  ptr.modify()->writeprotected = false;
}

void UserFileInDbData::setFileData( std::shared_ptr<const SpecMeas> spectrumFile,
                          const UserFileInDbData::SerializedFileFormat format )
{
  namespace io = boost::iostreams;

  if( !spectrumFile )
    return;
  
  fileData.clear();
  const size_t memsize = spectrumFile->memmorysize();
    
  //XXX - below guess on how much memorry to reserve is based on almost
  //      nothing
#if( ALLOW_SAVE_TO_DB_COMPRESSION )
  const size_t reserved_size = memsize;
  const size_t pre_mem_size_size
       = static_cast<size_t>( 2.2 * double(UserFileInDb::sm_maxFileSizeBytes) );
#else
  const size_t reserved_size = 8*1024+static_cast<size_t>(1.2*double(memsize));
  const size_t pre_mem_size_size = UserFileInDb::sm_maxFileSizeBytes;
#endif
    
  if( memsize > pre_mem_size_size )
    throw FileToLargeForDbException( memsize, pre_mem_size_size );

  fileData.reserve( reserved_size );
    
  try
  {
#if( ALLOW_SAVE_TO_DB_COMPRESSION )
    gzipCompressed = true;
    io::stream_buffer< io::back_insert_device< FileData_t > > buff( fileData );
    io::filtering_stream<io::output> outStream;
    outStream.push( io::gzip_compressor() );
    outStream.push( buff );
#else
    gzipCompressed = false;
    io::stream_buffer< io::back_insert_device< FileData_t > > buff( fileData );
    std::ostream outStream( &buff );
#endif
      
    switch( format )
    {
      case UserFileInDbData::k2012N42:
        spectrumFile->write_2012_N42( outStream );
        outStream << static_cast<unsigned char>(0);
        break;
    }//switch( format )
      
    fileFormat = format;
  }catch( std::exception &e )
  {
    cerr << "Failed to serialize spectrum to the database: "
         << e.what() << endl;
    throw runtime_error( "Failed to serialize spectrum to the database" );
  }//try / catch to serialize spectrumFile

  const size_t actual = fileData.size();
  if( actual > UserFileInDb::sm_maxFileSizeBytes )
  {
    fileData.clear();
    throw FileToLargeForDbException( actual, UserFileInDb::sm_maxFileSizeBytes );
  }//if( file is too big to save to database )
}//void UserFileInDbData::setFileData( std::shared_ptr<SpecMeas> spectrumFile )


std::shared_ptr<SpecMeas> UserFileInDbData::decodeSpectrum() const
{
  namespace io = boost::iostreams;
  std::shared_ptr<SpecMeas> spectrumFile( new SpecMeas() );

  try
  {
    if( !fileData.size() )
      throw runtime_error( "no data" );
    
    const char *start = (const char *)&fileData[0];
    const char *end = start + fileData.size();
    
    std::unique_ptr<std::istream> instrm;
    
    if( !gzipCompressed )
    {
      if( fileFormat != UserFileInDbData::k2012N42 )
        instrm.reset( new io::filtering_istream( boost::make_iterator_range(start,end) ) );
    }else
    {
#if( ALLOW_SAVE_TO_DB_COMPRESSION )
      io::filtering_stream<io::input> *instrmptr = new io::filtering_stream<io::input>();
      instrm.reset( instrmptr );
      instrmptr->push( io::gzip_decompressor() );
      instrmptr->push( boost::make_iterator_range( start, end ) );
#else
      throw runtime_error( "InterSpec built without zlip support,"
                          " cant de-serialize gzip compressed spectrum" );
#endif
    }//if( !gzipCompressed ) / else
    
    switch( fileFormat )
    {
      case UserFileInDbData::k2012N42:
      {
        bool loaded = false;
        if( gzipCompressed )
        {
          loaded = spectrumFile->SpecMeas::load_from_N42( *instrm );
        }else
        {
          --end;
          while( end != start && *end )
            --end;
          
#if( PERFORM_DEVELOPER_CHECKS )
          if( end == start )
            log_developer_error( __func__, "data from database wasnt null terminated" );
#endif
          
          if( end == start )
            throw runtime_error( "data from database wasnt null terminated" );
          
          loaded = spectrumFile->SpecMeas::load_N42_from_data( (char *)start, (char *)end );
        }
        
        if( !loaded )
          throw runtime_error( "Failed to load file from N42 format serialized to the database." );
      }//case UserFileInDbData::k2012N42:
        
    }//switch( format )

  }catch( std::exception &e )
  {
    stringstream strm;
    strm << "UserFileInDbData::decodeSpectrum(): caught "  << e.what()
         << " trying to deserialize DB entry";
    cerr << strm.str() << endl;
    throw runtime_error( strm.str() );
  }//try / catch
  
  return spectrumFile;
}//std::shared_ptr<SpecMeas> decodeSpectrum()

UserFileInDbData::UserFileInDbData()
 : gzipCompressed( false ),
   fileFormat( UserFileInDbData::k2012N42 )
{
  
}


boost::any UserOption::value() const
{
  if( m_name.empty() || ((m_type != String) && m_value.empty()) )
    throw runtime_error( "UserOption::value(): un-initialized UserOption" );
  
  switch( m_type )
  {
    case String:
      return m_value;
    case Decimal:
      return std::stod( m_value );
    case Integer:
      return std::stoi( m_value );
    case Boolean:
      return (m_value=="true" || m_value=="1");
  }//switch( m_type )
  
  throw runtime_error( "UserOption::value(): invalid m_type" );
  return boost::any();
}//boost::any value() const


InterSpecUser::InterSpecUser()
{
}


InterSpecUser::InterSpecUser( const std::string &username,
                              InterSpecUser::DeviceType type )
 : m_userName( username ),
   m_deviceType( type ),
   m_accessCount( 0 ),
   m_spectraFilesOpened( 0 ),
   m_firstAccessUTC( boost::posix_time::second_clock::universal_time() )
{
}//InterSpecUser::InterSpecUser(...)


void InterSpecUser::initFromDefaultValues( Wt::Dbo::ptr<InterSpecUser> user,
                          std::shared_ptr<DataBaseUtils::DbSession> session )
{
  using rapidxml::internal::compare;
  typedef rapidxml::xml_node<char>      XmlNode;
  typedef rapidxml::xml_attribute<char> XmlAttribute;

  if( !session )
    throw runtime_error( "InterSpecUser::initFromDefaultValues(...):"
                         " no valid session associated with user ptr" );
  if( !user->m_preferences.empty() )
    throw runtime_error( "InterSpecUser::initFromDefaultValues(...):"
                         " you cant call this function if preferences have"
                         " already been initialized" );
  try
  {
    user.modify();
  }catch(...)
  {
    throw runtime_error( "InterSpecUser::initFromDefaultValues(...):"
                         " there is no active transaction." );
  }
  
  const string filename = SpecUtils::append_path( InterSpec::staticDataDirectory(), sm_defaultPreferenceFile );
  
  std::vector<char> data;
  SpecUtils::load_file_data( filename.c_str(), data );
    
  rapidxml::xml_document<char> doc;
  const int flags = rapidxml::parse_normalize_whitespace
                    | rapidxml::parse_trim_whitespace;
    
  doc.parse<flags>( &data.front() );
  XmlNode *node = doc.first_node();
  if( !node || !node->name()
      || !compare( node->name(), node->name_size(), "preferences", 11, true) )
    throw runtime_error( "InterSpecUser: invlaid first node" );

  DataBaseUtils::DbTransaction transaction( *session );
  
  //we actually need to go through here and eliminate options where a "phone"
  //  or "tablet" option is available, and also rename ish...
  try
  {
    for( const XmlNode *pref = node->first_node( "pref", 4 );
         pref;
         pref = pref->next_sibling( "pref", 4 ) )
    {
      UserOption *option = parseUserOption( pref );
      try
      {
        user.modify()->m_preferences[option->m_name] = option->value();
        option->m_user = user;
        session->session()->add( option );
      }catch( std::exception &e )
      {
        const string value = pref->value() ? pref->value() : "";
        const string msg = "Value \"" + value + "\" is not convertible to the"
                           " intended type in " + filename
                           + " for pref "
                           + option->m_name + "\n" + string(e.what());
        delete option;
        throw runtime_error( msg );
      }//try / catch
    }//for( loop over preferences )

    transaction.commit();
  }catch( std::exception &e )
  {
    cerr << "\n\nInterSpecUser::initFromDefaultValues(...)\t" << e.what() << endl;
    user.modify()->m_preferences.clear();
    transaction.rollback();
    throw e;
  }//try / catch
  
}//void initFromDefaultValues()


void InterSpecUser::initFromDbValues( Wt::Dbo::ptr<InterSpecUser> user,
                          std::shared_ptr<DataBaseUtils::DbSession> session )
{
  
  if( !session )
    throw runtime_error( "InterSpecUser::initFromDbValues(...):"
                        " no valid session associated with user ptr" );
  if( !user->m_preferences.empty() )
    throw runtime_error( "InterSpecUser::initFromDbValues(...):"
                        " you cant call this function if prefernces have"
                        " already been initialized" );

  DataBaseUtils::DbTransaction transaction( *session );  
  const Wt::Dbo::collection< Wt::Dbo::ptr<UserOption> > &prefs = user->m_dbPreferences;
  
  vector< Dbo::ptr<UserOption> > options;
  std::copy( prefs.begin(), prefs.end(), std::back_inserter(options) );
  
  InterSpecUser *usr = user.modify();
  
  for( vector< Dbo::ptr<UserOption> >::const_iterator iter = options.begin();
      iter != options.end(); ++iter )
  {
    Dbo::ptr<UserOption> option = *iter;
    switch( option->m_type )
    {
      case UserOption::String:
        usr->m_preferences[option->m_name] = boost::any( option->m_value );
      break;
        
      case UserOption::Decimal:
      {
        double val;
        if( !(stringstream(option->m_value) >> val) )
          throw runtime_error( "InterSpecUser::initFromDbValues(): invalid"
                               " double value '" + option->m_value
                               + "' for user '" + user->m_userName + "'" );
        usr->m_preferences[option->m_name] = val;
        break;
      }//case UserOption::Decimal:
      
      case UserOption::Integer:
      {
        int val;
        if( !(stringstream(option->m_value) >> val) )
          throw runtime_error( "InterSpecUser::initFromDbValues(): invalid"
                              " int value '" + option->m_value
                              + "' for user '" + user->m_userName + "'" );
        usr->m_preferences[option->m_name] = val;
        break;
      }//case UserOption::Integer:
        
      case UserOption::Boolean:
      {
        bool val;
        if( !(stringstream(option->m_value) >> val) )
          throw runtime_error( "InterSpecUser::initFromDbValues(): invalid"
                              " int value '" + option->m_value
                              + "' for user '" + user->m_userName + "'" );
        usr->m_preferences[option->m_name] = val;
        break;
      }//case UserOption::Boolean:
    }//switch( option->m_type )
  }//for( loop over DB entries )
  
  //The call to usr->modify() was just to get access to to the UserOption ptr,
  // no changes have been made, so theres no reason to write to the database.
  transaction.rollback();
}//void initFromDbValues( Wt::Dbo::ptr<InterSpecUser> user )



UserOption *InterSpecUser::getDefaultUserPreference( const std::string &name,
                                                     const int type )
{
  using rapidxml::internal::compare;
  typedef rapidxml::xml_attribute<char> XmlAttribute;
  
  const bool isphone = ((type & InterSpecUser::PhoneDevice) != 0);
  const bool istablet = ((type & InterSpecUser::TabletDevice) != 0);
  const char *nameptr = name.c_str();
  const size_t namelen = name.length();
  
  const string filename = SpecUtils::append_path( InterSpec::staticDataDirectory(), sm_defaultPreferenceFile );
  
  std::vector<char> data;
  SpecUtils::load_file_data( filename.c_str(), data );
  
  rapidxml::xml_document<char> doc;
  const int flags = rapidxml::parse_normalize_whitespace
                    | rapidxml::parse_trim_whitespace;
  
  doc.parse<flags>( &data.front() );
  rapidxml::xml_node<char> *node = doc.first_node();
  if( !node || !node->name()
     || !compare( node->name(), node->name_size(), "preferences", 11, true) )
    throw runtime_error( "InterSpecUser: invalid first node" );
  
  const rapidxml::xml_node<char> *defnode = 0, *tabnode = 0, *phonenode = 0;
  
  for( const rapidxml::xml_node<char> *pref = node->first_node( "pref", 4 );
      pref;
      pref = pref->next_sibling( "pref", 4 ) )
  {
    const XmlAttribute *name_att = pref->first_attribute( "name", 4 );
    const char *prefname = name_att ? name_att->value() : "";
    const size_t prefnamelen = name_att ? name_att->value_size() : 0;
    const size_t substrlen = std::min(prefnamelen,name.size());
    
    //Check if the substring up to the '_' in "PreferenceName_phone" matches
    if( !compare(prefname, substrlen, nameptr, namelen, true) )
      continue;
  
    if( prefnamelen == namelen )
    {
      //If the preference name length matches the wanted name length, we have
      //  default value
      defnode = pref;
      
      //break out of the loop if this is all we need
      if( type == InterSpecUser::Desktop )
        break;
    }else if( isphone && compare(prefname+substrlen, prefnamelen-substrlen, "_phone", 6, true) )
    {
      //The preference name ends in "_phone", and the device is a phone, we
      //  dont need to keep looping.
      phonenode = pref;
      break;
    }else if( istablet && compare(prefname+substrlen, prefnamelen-substrlen, "_tablet", 7, true) )
    {
      //The preference name ends in "_tablet", and the device is a tablet, we
      //  can break as long as the device isnt also a phone.
      tabnode = pref;
      if( !isphone )
        break;
    }
  }//for( loop over preferneces )
    
  if( phonenode )
    return parseUserOption( phonenode );
  if( tabnode )
    return parseUserOption( tabnode );
  if( defnode )
    return parseUserOption( defnode );
  
  //Note: the string "couldn't find preference by name" is currently used in
  //      restoreUserPrefsFromXml(...) to check if a preference with this name is no longer used.
  //      So dont change next string without also changing there - or really, a more robust
  //      indication should be used.
  throw runtime_error( "InterSpecUser::getDefaultUserPreference(...):"
                       " couldn't find preference by name " + name );
}//UserOption *getDefaultUserPreference( const std::string &name )


const std::string &InterSpecUser::userName() const
{
  return m_userName;
}


int InterSpecUser::accessCount() const
{
  return m_accessCount;
}


std::chrono::steady_clock::time_point::duration InterSpecUser::totalTimeInApp() const
{
  return chrono::microseconds( m_totalTimeInApp.total_microseconds() );
}

std::chrono::system_clock::time_point InterSpecUser::currentAccessStartUTC() const
{
  return to_time_point( m_currentAccessStartUTC );
}


int InterSpecUser::numSpectraFilesOpened() const
{
  return m_spectraFilesOpened;
}


std::chrono::system_clock::time_point InterSpecUser::firstAccessUTC() const
{
  return to_time_point( m_firstAccessUTC );
}


Wt::Dbo::ptr<InterSpecUser> &InterSpecUser::userFromViewer( InterSpec *viewer )
{
  static Dbo::ptr<InterSpecUser> dummy;
  if( !viewer )
    return dummy;
  return viewer->m_user;
}

std::shared_ptr<DataBaseUtils::DbSession> InterSpecUser::sqlFromViewer( InterSpec *viewer )
{
  if( !viewer )
    return std::shared_ptr<DataBaseUtils::DbSession>();
  return viewer->sql();
}

void InterSpecUser::startingNewSession()
{
  incrementAccessCount();
  setPreviousAccessTime( to_time_point( m_currentAccessStartUTC ) );
  setCurrentAccessTime( chrono::system_clock::now() );
}//void startingNewSession()
  

void InterSpecUser::addUsageTimeDuration( const std::chrono::steady_clock::time_point::duration &duration )
{
  m_totalTimeInApp += boost::posix_time::microseconds( chrono::duration_cast<chrono::microseconds>(duration).count() );
}//void addUsageTimeDuration(...)


void InterSpecUser::incrementSpectraFileOpened()
{
  ++m_spectraFilesOpened;
}//void incrementSpectraFileOpened()


void InterSpecUser::restoreUserPrefsFromXml(
                                Wt::Dbo::ptr<InterSpecUser> user,
                                const ::rapidxml::xml_node<char> *prefs_node,
                                InterSpec *viewer )
{
  using namespace rapidxml;
  using rapidxml::internal::compare;
  
  if( !user || !prefs_node
      || !compare( prefs_node->name(), prefs_node->name_size(), "preferences", 11, true) )
    throw runtime_error( "restoreUserPrefsFromXml: invalid input" );
  
 
  for( const xml_node<char> *pref = prefs_node->first_node( "pref", 4 ); pref;
       pref = pref->next_sibling( "pref", 4 ) )
  {
    std::unique_ptr<UserOption> option;
    try
    {
      option.reset( parseUserOption(pref) );
      
      switch( option->m_type )
      {
        case UserOption::String:
        {
          const string value = boost::any_cast<string>( option->value() );
          setPreferenceValue( user, option->m_name, value, viewer );
          break;
        }//case String
          
        case UserOption::Decimal:
        {
          const double value = boost::any_cast<double>( option->value() );
          setPreferenceValue( user, option->m_name, value, viewer );
          break;
        }//case Decimal
          
        case UserOption::Integer:
        {
          const int value = boost::any_cast<int>( option->value() );
          setPreferenceValue( user, option->m_name, value, viewer );
          break;
        }//case Integer
          
        case UserOption::Boolean:
        {
          const bool value = boost::any_cast<bool>( option->value() );
          setPreferenceValue( user, option->m_name, value, viewer );
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


::rapidxml::xml_node<char> *InterSpecUser::userOptionsToXml(
                                    ::rapidxml::xml_node<char> *parent_node,
                                              InterSpec *viewer ) const
{
  using namespace ::rapidxml;
  
  if( !parent_node )
    throw runtime_error( "userOptionsToXml: invalid input" );
  
  xml_document<char> *doc = parent_node->document();
  
  xml_node<char> *prefs_node = doc->allocate_node( node_element, "preferences" );
  parent_node->append_node( prefs_node );
  
  vector< Dbo::ptr<UserOption> > options;
  
  {//begin codeblock to retrieve prefernces from database
    std::shared_ptr<DataBaseUtils::DbSession> sql = viewer->sql();
    DataBaseUtils::DbTransaction transaction( *sql );

    std::copy( m_dbPreferences.begin(), m_dbPreferences.end(),
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




const Wt::Dbo::collection< Wt::Dbo::ptr<UserFileInDb> > &InterSpecUser::userFiles() const
{
  return m_userFiles;
}

const Wt::Dbo::collection< Wt::Dbo::ptr<ShieldingSourceModel> > &InterSpecUser::shieldSrcModels() const
{
  return m_shieldSrcModels;
}

const Wt::Dbo::collection< Wt::Dbo::ptr<UserState> > &InterSpecUser::userStates() const
{
  return m_userStates;
}

const Wt::Dbo::collection< Wt::Dbo::ptr<ColorThemeInfo> > &InterSpecUser::colorThemes() const
{
  return m_colorThemes;
}

const Wt::Dbo::collection< Wt::Dbo::ptr<UseDrfPref> > &InterSpecUser::drfPrefs() const
{
  return m_drfPref;
}

void InterSpecUser::incrementAccessCount()
{
  ++m_accessCount;
}//void incrementAccessCount();


void InterSpecUser::setCurrentAccessTime( const std::chrono::system_clock::time_point &utcTime )
{
  m_currentAccessStartUTC = to_ptime( utcTime );
}//void setCurrentAccessTime( const boost::posix_time::ptime &utcTime )


void InterSpecUser::setPreviousAccessTime( const std::chrono::system_clock::time_point &utcTime )
{
  m_previousAccessUTC = to_ptime( utcTime );
}//void setPreviousAccessTime( const boost::posix_time::ptime &utcTime );
