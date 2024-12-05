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
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <assert.h>
#include <iostream>

#include <boost/iostreams/stream.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/back_inserter.hpp>


#include <Wt/Dbo/Dbo>
#include <Wt/WString>
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


FileToLargeForDbException::FileToLargeForDbException( const size_t saveSize, const size_t limit )
  : std::exception(),
    m_msg{},
    m_saveSize( saveSize ), 
    m_limit( limit )
{
  m_msg = message().toUTF8();
}

Wt::WString FileToLargeForDbException::message() const
{
  return WString::tr( "err-save-to-large-for-db" )
    .arg( static_cast<int>(m_saveSize/1024) )
    .arg( static_cast<int>(m_limit/1024) );
}


UserState::UserState()
  : stateType( kUndefinedStateType ),
    creationTime( Wt::WDateTime::currentDateTime() ),
    serializeTime( Wt::WDateTime::currentDateTime() ),
    foreground(), background(), secondForeground(),
    energyAxisMinimum( 0.0 ), energyAxisMaximum( 3000.0 ),
    countsAxisMinimum( -1.0 ), countsAxisMaximum( -1.0 ), displayBinFactor( 1 ),
    shownDisplayFeatures( 0 ), backgroundSubMode( kNoSpectrumSubtract ),
    currentTab( kNoTabs ), showingMarkers( 0 ), disabledNotifications( 0 ),
    showingPeakLabels( 0 ), showingWindows( 0 )
{
}


void UserState::removeFromDatabase( Wt::Dbo::ptr<UserState> state,
                                          DataBaseUtils::DbSession &session )
{
  if( !state || (state.id() < 0) || !state->user )
    throw runtime_error( "UserFileInDbData::removeFromDatabase: invalid input." );
  
  DataBaseUtils::DbTransaction transaction( session );
  
  session.session()->execute( "PRAGMA foreign_keys = ON;" );
  
  // Delete all the tags of this state (I dont *think we need to do this....)
  /*
  for( auto iter = state->snapshotTags.begin(); iter != state->snapshotTags.end(); ++iter )
  {
    auto &kid = *iter;
    assert( kid != state );
    if( kid == state )
      continue;
    
    try
    {
      UserState::removeFromDatabase( kid, session );
    }catch( std::exception & )
    {
      transaction.rollback();
      throw;
    }
  }//for( const auto &kid : children )
   */
  
  try
  {
    state.remove();
    transaction.commit();
  }catch( std::exception &e )
  {
    cerr << "removeFromDatabase(...) caught error: " << e.what() << endl;
    transaction.rollback();
    throw;
  }//try / catch
}//removeFromDatabase( Dbo::ptr<UserState> state );


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
  isPartOfSaveState = false;
}


Dbo::ptr<UserFileInDb> UserFileInDb::makeDeepCopyOfFileInDatabase(
                                                  Dbo::ptr<UserFileInDb> orig,
                                                  DataBaseUtils::DbSession &sqldb,
                                                  bool isSaveState )
{
  if( !orig )
    return orig;
  
  auto session = sqldb.session();
  if( !session )
    throw runtime_error( "UserFileInDb::makeDeepCopyOfFileInDatabase():"
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
  
  newfile->isPartOfSaveState = isSaveState;

  newfile->snapshotParent = orig->snapshotParent;
  
  DataBaseUtils::DbTransaction transaction( sqldb );
  
  try
  {
    Dbo::ptr<UserFileInDb> answer = session->add( newfile );
  
    // The passed in `UserFileInDb` may not actually be in the database, so we have to check
    //  for this.
    if( orig.id() >= 0 )
    {
      Wt::Dbo::ptr<UserFileInDbData> data = orig->filedata.lock();
      if( !data )
        throw runtime_error( "No filedata db entry" );
      
      UserFileInDbData *newdata = new UserFileInDbData( *data );
      newdata->fileInfo = answer;
      session->add( newdata );
      
      for( Dbo::collection< Dbo::ptr<ShieldingSourceModel> >::const_iterator iter
          = orig->modelsUsedWith.begin();
          iter != orig->modelsUsedWith.end();
          ++iter )
      {
        ShieldingSourceModel *newmodel = new ShieldingSourceModel();
        newmodel->shallowEquals( **iter );
        Dbo::ptr<ShieldingSourceModel> modelptr = session->add( newmodel );
        modelptr.modify()->filesUsedWith.insert( answer );
      }//for( loop over ShieldingSourceModel )
    }//if( orig.id() >= 0 )
    
    transaction.commit();
    return answer;
  }catch( std::exception &e )
  {
    transaction.rollback();
    throw runtime_error( e.what() );
  }//try catch
  
  return Dbo::ptr<UserFileInDb>();
}//makeDeepCopyOfFileInDatabase(...)


void ShieldingSourceModel::shallowEquals( const ShieldingSourceModel &rhs )
{
  user = rhs.user;
  name = rhs.name;
  description = rhs.description;
  creationTime = rhs.creationTime;
  serializeTime = rhs.serializeTime;
  //filesUsedWith = rhs.filesUsedWith;
  xmlData = rhs.xmlData;
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
    {
      double val = 0.0;
      const bool parsed = SpecUtils::parse_double( m_value.c_str(), m_value.size(), val );
      assert( parsed );
      //if( !parsed )
      //  throw std::logic_error( "Invalid Decimal str: '" + m_value + "'" );
      
      return val;
      //return std::stod( m_value );
    }//case Decimal:
    
    case Integer:
    {
      int val = 0;
      const bool parsed = SpecUtils::parse_int( m_value.c_str(), m_value.size(), val );
      assert( parsed );
      //if( !parsed )
      //  throw std::logic_error( "Invalid Integer str: '" + m_value + "'" );
      
      return val;
      //return std::stoi( m_value );
    }//case Integer:
    
    case Boolean:
    {
      assert( (m_value == "0") || (m_value == "1") );
      return (m_value=="true" || m_value=="1");
    }
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


const std::string &InterSpecUser::userName() const
{
  return m_userName;
}

int InterSpecUser::deviceType() const
{
  return m_deviceType;
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


const Wt::Dbo::collection< Wt::Dbo::ptr<UserFileInDb> > &InterSpecUser::userFiles() const
{
  return m_userFiles;
}

const Wt::Dbo::collection< Wt::Dbo::ptr<UserOption> > &InterSpecUser::preferences() const
{
  return m_dbPreferences;
}

const Wt::Dbo::collection< Wt::Dbo::ptr<ShieldingSourceModel> > &InterSpecUser::shieldSrcModels() const
{
  return m_shieldSrcModels;
}

//const Wt::Dbo::collection< Wt::Dbo::ptr<UserState> > &InterSpecUser::userStates() const
//{
//  return m_userStates;
//}

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
