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

#include <deque>
#include <vector>
#include <string>
#include <memory>
#include <sstream>
#include <fstream>
#include <iostream>

// #include <fcntl.h>

#include <stdio.h>
#include <stdlib.h>


// Block out some warnings occurring in boost files.
#pragma warning(disable:4244) // warning C4244: 'initializing' : conversion from 'std::streamoff' to 'size_t', possible loss of data
#pragma warning(disable:4308) // warning C4308: negative integral constant converted to unsigned type

#include <boost/any.hpp>
#include <boost/system/error_code.hpp>

#include <Wt/WText>
#include <Wt/Utils>
#include <Wt/WLabel>
#include <Wt/WAnchor>
#include <Wt/WString>
#include <Wt/Dbo/Dbo>
#include <Wt/WServer>
#include <Wt/WIconPair>
#include <Wt/WTreeView>
#include <Wt/WGroupBox>
#include <Wt/WResource>
#include <Wt/WDateTime>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WModelIndex>
#include <Wt/WFileUpload>
#include <Wt/WButtonGroup>
#include <Wt/WRadioButton>
#include <Wt/WApplication>
#include <Wt/WProgressBar>
#include <Wt/Http/Response>
#include <Wt/WItemDelegate>
#include <Wt/WBorderLayout>
#include <Wt/Dbo/WtSqlTraits>
#include <Wt/WContainerWidget>
#include <Wt/WAbstractItemModel>

#include "InterSpec/PeakDef.h"
#include "SpecUtils/SpecFile.h"
#include "InterSpec/PopupDiv.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/WarningWidget.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/RowStretchTreeView.h"

#if( SpecUtils_ENABLE_D3_CHART )
#include "SpecUtils/D3SpectrumExport.h"
#endif

using namespace Wt;
using namespace std;

namespace
{
  bool is_close_to_int( double t )
  {
    return (fabs(t-floor(t)) < 0.01 );
  }
  
  void giveMessageToApp( const string msg,
                        const WarningWidget::WarningMsgLevel level )
  {
    passMessage( msg, level );
    if( wApp )  //wApp should aways be valid here
      wApp->triggerUpdate();
  }
  
  //postMessageToApp(...): thread safe message to call passMessage(...), incase
  //  the call arises after wApp no longer exists.
  void postMessageToApp( const string &msg,
                        const WarningWidget::WarningMsgLevel level,
                        const string &appId )
  {
    WServer *server = WServer::instance();  //can this ever be NULL?
    if( server )
      server->post( appId, boost::bind( &giveMessageToApp, msg, level ) );
  }
  
  const char * const error_saving_spectrum_msg
                           = "There was an error saving a spectrum to the"
                             " database from memory. It shouldnt effect"
                             " this session, but possibly future ones using"
                             " the spectrum.";
  
  const char * const error_saving_spectrum_size_msg
                              = "Spectrum too large to save to database.";
  
  const char * const error_saving_spectrum_stale_msg
                              = "You may be saving the same spectrum to the"
                                " database in mutliple sessions, becareful you"
                                " dont over-write work in one session from"
                                " another.";
}//namespace




SpectraHeader::SpectraHeader()
{
  live_time = real_time = 0.0;
  contained_neutron_ = false;
  sample_number = -1;
  gamma_counts_ = neutron_counts_ = 0.0;
  spectra_type = SpecUtils::SourceType::Unknown;
  is_derived_data = false;
}//SpectraHeader default constructor


SpectraHeader::SpectraHeader( const std::vector<std::shared_ptr<const SpecUtils::Measurement>> &sample_measurements )
{
  init( sample_measurements );
}


void SpectraHeader::init( const std::vector<std::shared_ptr<const SpecUtils::Measurement>> &measurements )
{
  live_time = real_time = 0.0;
  contained_neutron_ = false;
  sample_number = -1;
  gamma_counts_ = neutron_counts_ = 0.0;
  spectra_type = SpecUtils::SourceType::Unknown;
  is_derived_data = false;

  for( const std::shared_ptr<const SpecUtils::Measurement> &m : measurements )
  {
    live_time += m->live_time();
    real_time += m->live_time();
    contained_neutron_ |= m->contained_neutron();
    gamma_counts_ += m->gamma_count_sum();
    neutron_counts_ += m->neutron_counts_sum();
    detector_names.push_back( m->detector_name() );
    detector_numbers_.push_back( m->detector_number() );
    spectra_type = m->source_type();
    if( m->derived_data_properties() )
      is_derived_data = true;

    for( size_t i = 0; i < m->remarks().size(); ++i )
    {
      if( !remarks.empty() )
        remarks += ", ";
      string thisRemark = m->remarks()[i];

      size_t pos = thisRemark.find( m->detector_name() );
      if( pos != string::npos )
        thisRemark.erase( pos, m->detector_name().length() );

      pos = thisRemark.find( "Speed" );
      if( pos != string::npos )
      {
        size_t end = thisRemark.find_first_not_of( "= \t", pos+5 );
        thisRemark.erase( pos, end-pos );
      }//if( pos != string::npos )

      pos = thisRemark.find( "Survey" );
      if( pos != string::npos )
      {
        size_t end = thisRemark.find_first_not_of( " \t=", pos+6 );
        thisRemark.erase( pos, end-pos );
      }//if( pos != string::npos )
    }//for( const string &s : m->remarks )

    //In principle we should check that if we've alread filled in below,
    // we arent filling in a new different value

    char buffer[32];
    snprintf( buffer, sizeof(buffer), "%.1f m/s", m->speed() );
    speed_ = buffer;
    sample_number = m->sample_number();
    start_time = WDateTime::fromPosixTime( to_ptime( m->start_time() ) );
  }//for( std::shared_ptr<const SpecUtils::Measurement> &m : measurements )

  live_time /= measurements.size();
  real_time /= measurements.size();
}//SpectraHeader constructor


SpectraHeader::~SpectraHeader()
{
}


SpectraFileHeader::SpectraFileHeader( bool keepInMemmory,
                                      InterSpec *viewer )
{
  assert( viewer );
  m_viewer = viewer;
  m_sql = viewer->sql();
  m_fileSystemLocation = "";
  m_totalLiveTime = m_totalRealTime = 0.0;
  m_totalGammaCounts = m_totalNeutronCounts = 0.0;
  m_numDetectors = -1;
  m_hasNeutronDetector = false;
  m_keepCache = keepInMemmory;
  m_user = viewer->user();
  m_modifiedSinceDecode = false;
  m_candidateForSavingToDb = true;
  m_app = wApp;
  if( m_app )
    m_appId = m_app->sessionId();
}//SpectraFileHeader constructor


SpectraFileHeader::~SpectraFileHeader() noexcept(true)
{
  try
  {
    string fileSystemLocation;
    std::shared_ptr<SpecMeas> memObj;
    
    {
      RecursiveLock lock( m_mutex );
      m_aboutToBeDeletedConnection.disconnect();
      memObj = m_weakMeasurmentPtr.lock();
      fileSystemLocation = m_fileSystemLocation;
    }
  
    if( fileSystemLocation.size() )
    {
      const bool status = SpecUtils::remove_file( fileSystemLocation );;
      if( !status )
        throw runtime_error( m_fileSystemLocation + " didn't exist to delete" );
    }// if we should delete the file
  }catch( ... )
  {
    //cerr << "\n~SpectraFileHeader() caught: " << e.what() << endl;
  }
}//SpectraFileHeader destructor


#if( USE_DB_TO_STORE_SPECTRA )
std::shared_ptr<SpecMeas> SpectraFileHeader::resetFromDatabase(
                                                Dbo::ptr<UserFileInDb> info )
{
  Dbo::ptr<UserFileInDbData> dbdata;
  
  {//begin interaction with Database
    DataBaseUtils::DbTransaction transaction( *m_sql );
    
    if( !info || info->filedata.size() < 1 )
      throw runtime_error( "SpectraFileHeader::resetFromDatabase(...):"
                           "No UserFileInDbData for this UserFileInDb" );
  
    dbdata = m_sql->session()->find<UserFileInDbData>()
                             .where( "UserFileInDb_id = ?" )
                             .bind( info->filedata.front().id() );
    transaction.commit();
  }//end interaction with Database
 
  RecursiveLock lock( m_mutex );
  std::shared_ptr<SpecMeas> memobj = dbdata->decodeSpectrum();
  setMeasurmentInfo( memobj );
  
  m_displayName = SpecUtils::filename( info->filename );
  m_uploadTime = info->uploadTime;
  m_modifiedSinceDecode = info->userHasModified;
  
  m_fileDbEntry = info;

  return memobj;
}//std::shared_ptr<SpecMeas> resetFromDatabase( Dbo::ptr<UserFileInDb> dbfile )


void SpectraFileHeader::setNotACandiateForSavingToDb()
{
  m_candidateForSavingToDb = false;
}

bool SpectraFileHeader::candidateForSavingToDb() const
{
  return m_candidateForSavingToDb;
}

bool SpectraFileHeader::shouldSaveToDb() const
{
  return (m_candidateForSavingToDb && m_user );
}

void SpectraFileHeader::setDbEntry( Wt::Dbo::ptr<UserFileInDb> entry )
{
  Dbo::ptr<UserFileInDb> initialentry;
  
  {
    RecursiveLock lock( m_mutex );
    if( m_fileDbEntry == entry )
      return;
    initialentry = m_fileDbEntry;
    m_fileDbEntry.reset();
  }
  
  if( initialentry )
  {
    DataBaseUtils::DbTransaction transaction( *m_sql );
    initialentry.remove();
    transaction.commit();
  }//if( m_fileDbEntry )
  
  if( entry && entry.id() < 0 )
  {
    DataBaseUtils::DbTransaction transaction( *m_sql );
    m_sql->session()->add( entry );
    transaction.commit();
  }//if( entry && entry.isTransient() )
  
  {
    RecursiveLock lock( m_mutex );
    m_fileDbEntry = entry;
  }
}//void setDbEntry( Wt::Dbo::ptr<UserFileInDb> entry )


void SpectraFileHeader::clearDbEntry()
{
  RecursiveLock lock( m_mutex );
  m_fileDbEntry.reset();
}


Wt::Dbo::ptr<UserFileInDb> SpectraFileHeader::dbEntry()
{
  return m_fileDbEntry;
}//Wt::Dbo::ptr<UserFileInDb> dbEntry()


void SpectraFileHeader::setBasicFileInDbInfo( UserFileInDb *info ) const
{
  RecursiveLock lock( m_mutex );
  
  info->user = m_user;
  info->uuid = m_uuid;
  if( static_cast<int>(info->uuid.length()) > UserFileInDb::sm_maxUuidLength )
    info->uuid = info->uuid.substr( 0, UserFileInDb::sm_maxUuidLength );
  info->sessionID       = m_appId;
  if( static_cast<int>(info->uuid.length()) > UserFileInDb::sm_maxSessionIdLength )
    info->sessionID = info->sessionID.substr( 0, UserFileInDb::sm_maxSessionIdLength );
  
  info->filename      = m_displayName; //meas->filename();
  info->uploadTime    = m_uploadTime;
  info->serializeTime = WDateTime::currentDateTime();
  info->numSamples    = m_numSamples;
  info->isPassthrough = m_isPassthrough;
  info->totalLiveTime = m_totalLiveTime;
  info->totalRealTime = m_totalRealTime;
  info->totalGammaCounts      = m_totalGammaCounts;
  info->totalNeutronCounts    = m_totalNeutronCounts;
  info->numDetectors          = m_numDetectors;
  info->hasNeutronDetector    = m_hasNeutronDetector;
  info->measurementsStartTime = m_spectrumTime;
  info->userHasModified       = m_modifiedSinceDecode;
//  info->spectrumType
//  info->backgroundSpectrumId = -1;
}//void setBasicFileInDbInfo( UserFileInDb &entry )



void SpectraFileHeader::saveToDatabaseFromTempFileWorker() const
{
  string msg;
  WarningWidget::WarningMsgLevel level = WarningWidget::WarningMsgInfo;

  
  try
  {
    saveToDatabaseFromTempFile();
  }catch( FileToLargeForDbException &e )
  {
    msg = error_saving_spectrum_size_msg;
    msg += string(" ") + e.what();
    level = WarningWidget::WarningMsgInfo;
  }catch( Wt::Dbo::StaleObjectException & )
  {
    try
    {
      {
        DataBaseUtils::DbTransaction transaction( *m_sql );
        m_fileDbEntry.reread();
        transaction.commit();
      }
      
      saveToDatabaseFromTempFile();
      
      msg = error_saving_spectrum_stale_msg;
      level = WarningWidget::WarningMsgInfo;
    }catch( std::exception &e )
    {
      msg = error_saving_spectrum_msg;
      level = WarningWidget::WarningMsgHigh;
      cerr << "\n\nError saving spectrum to db from tmp file (2nd attempt): "
           << e.what() << endl;
    }//try / catch
  }catch( std::exception &e )
  {
    cerr << "\n\nThere was an error saving to the database from the "
         << "temporary file: " << e.what() << endl << endl;
    msg = error_saving_spectrum_msg;
    level = WarningWidget::WarningMsgHigh;
  }//try / catch
  
  
  if( !m_appId.empty() && msg.size() )
    postMessageToApp( msg, level, m_appId );
}//saveToDatabaseFromTempFileWorker(...)


void SpectraFileHeader::saveToDatabaseWorker(
                            std::shared_ptr<SpecMeas> measurment,
                            std::shared_ptr<SpectraFileHeader> header )
{
  string msg;
  WarningWidget::WarningMsgLevel level = WarningWidget::WarningMsgInfo;
  
  try
  {
    SpectraFileHeader::saveToDatabase( measurment, header );
  }catch( FileToLargeForDbException &e )
  {
    msg = error_saving_spectrum_size_msg;
    msg += string(" ") + e.what();
    level = WarningWidget::WarningMsgInfo;
  }catch( Wt::Dbo::StaleObjectException & )
  {
    try
    {
      {
        DataBaseUtils::DbTransaction transaction( *(header->m_sql) );
        
        header->m_fileDbEntry.reread();
        transaction.commit();
      }
    
      SpectraFileHeader::saveToDatabase( measurment, header );
      
      msg = error_saving_spectrum_stale_msg;
      level = WarningWidget::WarningMsgInfo;
    }catch( std::exception &e )
    {
      msg = error_saving_spectrum_msg;
      level = WarningWidget::WarningMsgHigh;
      cerr << "\n\nError saving spectrum to database (2nd attempt): "
           << e.what() << endl;
    }//try / catch
  }catch( std::exception &e )
  {
    msg = error_saving_spectrum_msg;
    level = WarningWidget::WarningMsgHigh;
    cerr << "\n\nError saving spectrum to database: " << e.what() << endl;
  }//try / catch
  
  if( !header->m_appId.empty() && msg.size() )
    postMessageToApp( msg, level, header->m_appId );
}//void saveToDatabaseWorker( std::shared_ptr<SpecMeas> measurment )


void SpectraFileHeader::saveToDatabaseFromTempFile() const
{
  Dbo::ptr<UserFileInDb> fileDbEntry;
  
  {//begin locked section
    RecursiveLock lock( m_mutex );
    if( !shouldSaveToDb() )
      return;

    if( m_fileSystemLocation.empty() )
      throw runtime_error( "SpectraFileHeader::saveToDatabaseFromTempFile():"
                           " no cached file");
    
    fileDbEntry = m_fileDbEntry;
  }//end locked section
  
  Dbo::ptr<UserFileInDbData> data;

  //There is a chance m_fileDbEntry has gone stale
  vector< Dbo::ptr<UserFileInDbData> > files;
  if( fileDbEntry )
  {
    DataBaseUtils::DbTransaction transaction( *m_sql );
    
    //grab a fresh copy of fileDbEntry to avoid any thread safety issues
    fileDbEntry = m_sql->session()->find<UserFileInDb>()
                                  .where( "id = ?" )
                                  .bind( fileDbEntry.id() );
    if( fileDbEntry )
      std::copy( fileDbEntry->filedata.begin(), fileDbEntry->filedata.end(),
                  std::back_inserter(files) );
    
    transaction.commit();
  }//if( m_fileDbEntry )

  if( files.size() )
  {
    data = files[0];
    DataBaseUtils::DbTransaction transaction( *m_sql );
    fileDbEntry.modify()->serializeTime = WDateTime::currentDateTime();
    transaction.commit();
  }else
  {
    DataBaseUtils::DbTransaction transaction( *m_sql );
    UserFileInDb *info = new UserFileInDb();
    setBasicFileInDbInfo( info );
    fileDbEntry = m_sql->session()->add( info );
    
    UserFileInDbData *dataptr = new UserFileInDbData();
    dataptr->fileInfo = fileDbEntry;
    data = m_sql->session()->add( dataptr );
    transaction.commit();
  }//if( m_fileDbEntry )
  
  
  try
  {
    DataBaseUtils::DbTransaction transaction( *m_sql );
    data = m_sql->session()->find<UserFileInDbData>()
                            .where( "UserFileInDb_id = ?" )
                            .bind( data.id() );
    if( data )
      data.modify()->setFileData( m_fileSystemLocation,
                               UserFileInDbData::sm_defaultSerializationFormat );
    transaction.commit();
  }catch( FileToLargeForDbException &e )
  {
    try
    {
      DataBaseUtils::DbTransaction transaction( *m_sql );
      if( data )
        data.remove();
      if( fileDbEntry )
        fileDbEntry.remove();
      fileDbEntry.reset();
      transaction.commit();
    }catch(...){}
    
    {
      RecursiveLock lock( m_mutex );
      m_fileDbEntry = fileDbEntry;
    }
    
    cerr << "SpectraFileHeader::saveToDatabaseFromTempFile: File too large to save to DB" << endl;
    throw e;
  }catch( std::exception &e )
  {
    try
    {
      DataBaseUtils::DbTransaction transaction( *m_sql );
      if( data )
        data.remove();
      if( fileDbEntry )
        fileDbEntry.remove();
      fileDbEntry.reset();
      transaction.commit();
    }catch(...){}
    
    {
      RecursiveLock lock( m_mutex );
      m_fileDbEntry = fileDbEntry;
    }
    
    cerr << "SpectraFileHeader::saveToDatabaseFromTempFile: caught: " << e.what() << endl;
    throw runtime_error( e.what() );
  }//try / catch

  
  {
    RecursiveLock lock( m_mutex );
    m_fileDbEntry = fileDbEntry;
  }
  
  //XXX - a Wt::Dbo::StaleObjectException may be thrown if modifying the same
  //      UserFileInDbData object in another session - but this is unlikely
  //      to happen, so we wont worry about it
}//void saveToDatabaseFromTempFile() const


void SpectraFileHeader::saveToDatabase( std::shared_ptr<SpecMeas> meas,
                                        std::shared_ptr<SpectraFileHeader> header )
{
  header->saveToDatabase( meas );
}

void SpectraFileHeader::saveToDatabase( std::shared_ptr<const SpecMeas> input ) const
{
  if( !input )
    throw runtime_error( "\n\n\nSpectraFileHeader::saveToDatabase(): !input" );
  
  Wt::log("info") << "Saving '" << input->filename() << " to DB.";
  
  std::shared_ptr<const SpecMeas> meas;
 
  {
    RecursiveLock lock( m_mutex );
    meas = m_weakMeasurmentPtr.lock();
  
    if( !m_weakMeasurmentPtr.expired() && meas != m_weakMeasurmentPtr.lock() )
      throw runtime_error( "SpectraFileHeader::saveToDatabase(): meas != m_weakMeasurmentPtr" );
  
    if( !m_user || !m_user.session() )
      throw runtime_error( "SpectraFileHeader::saveToDatabase(): !m_user || !m_user.session()" );
  }
  
  //In principle we should take a recursive mutex lock on meas for this entire function
  //  or do a deep copy...
  
#if( PERFORM_DEVELOPER_CHECKS && SpecUtils_ENABLE_EQUALITY_CHECKS )
  {//begin code block to do check
    string filepath = SpecUtils::temp_file_name("DevTestSaveToDb", InterSpecApp::tempDirectory() );
    const string filename = filepath + ".n42";
    
    bool written = false;
    
    {
#ifdef _WIN32
      const std::wstring wfilename = SpecUtils::convert_from_utf8_to_utf16(filename);
      ofstream output( wfilename.c_str(), ios::binary | ios::out );
#else
      ofstream output( filename.c_str(), ios::binary | ios::out );
#endif
      if( !output.is_open() )
      {
        log_developer_error( __func__, "Failed to open temp output file" );
      }else
      {
        written = meas->write_2012_N42( output );
        if( !written )
          log_developer_error( __func__, "Failed to write write_2012_N42" );
      }
    }
    
    if( written )
    {
      SpecMeas newmeas;
      const bool read = newmeas.load_N42_file( filename );
      
      if( !read )
      {
        log_developer_error( __func__, "Failed to re-read in 2012 N42 file" );
      }else
      {
        newmeas.set_filename( meas->filename() );
        try
        {
          SpecMeas::equalEnough( *meas, newmeas );
          
          cout << "Saving SpecMeas to file and reading it back in worked" << endl;
        }catch( std::exception &e )
        {
          log_developer_error( __func__,
                               ("Failed check comparing re-serialized SpecMeas object: " + string(e.what())).c_str() );
        }
      }
    }//if( written )
    
    SpecUtils::remove_file( filename );
  }//end code block to do check
#endif //#if( PERFORM_DEVELOPER_CHECKS )
  
  Wt::Dbo::ptr<UserFileInDb> fileDbEntry;
  
  bool modifiedSinceDecode;
  
  {
    RecursiveLock lock( m_mutex );
    modifiedSinceDecode = m_modifiedSinceDecode = meas->modified_since_decode();
    fileDbEntry = m_fileDbEntry;
  }
  
  
  if( !fileDbEntry )
  {
    UserFileInDb *info = new UserFileInDb();
    Dbo::ptr<UserFileInDb> info_dbo_ptr( info );
    setBasicFileInDbInfo( info );  //takes lock of m_mutex during execution
    UserFileInDbData *data = new UserFileInDbData();
    Dbo::ptr<UserFileInDbData> data_dbo_ptr( data );
    
    try
    {
      data->setFileData( meas, UserFileInDbData::sm_defaultSerializationFormat );
    }catch( FileToLargeForDbException &e )
    {
      throw e;
    }catch( std::exception &e )
    {
      throw runtime_error( e.what() );
    }//try / catch
    
    {
      DataBaseUtils::DbTransaction transaction( *m_sql );
      fileDbEntry = m_sql->session()->add( info_dbo_ptr );
      data->fileInfo = fileDbEntry;
      Dbo::ptr<UserFileInDbData> dataPtr = m_sql->session()->add( data_dbo_ptr );
      transaction.commit();
    }
  }else
  {
    //TODO: couldnt the below query just be accomplished by `fileDbEntry.reread()`?
    try
    {
      DataBaseUtils::DbTransaction transaction( *m_sql );
      fileDbEntry = m_sql->session()->find<UserFileInDb>().where("id = ?").bind(fileDbEntry.id()).resultValue();
      transaction.commit();
      
      if( !fileDbEntry )
        runtime_error("no entry");
    }catch( Wt::Dbo::Exception &e )
    {
      throw runtime_error( "SpectraFileHeader::saveToDatabase(), database error: " + string(e.what()) );
    }catch( std::exception & )
    {
      throw runtime_error( "SpectraFileHeader::saveToDatabase(): error re-reading UserFileInDb in DB" );
    }//try / catch
    
    typedef Dbo::collection<Dbo::ptr<UserFileInDbData> > Files;
    
    Dbo::ptr<UserFileInDbData> data;
    
    {//begin interact with database
      DataBaseUtils::DbTransaction transaction( *m_sql );
      for( Files::const_iterator i = fileDbEntry->filedata.begin();
          i != fileDbEntry->filedata.end(); ++i )
      {
        data = *i;
        break;
      }
      
      transaction.commit();
    }//end interact with database
    
    if( data )
    {
      bool modified;
      
      {
        RecursiveLock lock( m_mutex );
        modified = meas->modified();
      }

      if( !modified )
      {
        cerr << "Spectrum '" << input->filename() << "' has not been modified, not re-writing to database" << endl;
        return;
      }//if( !meas->modified() )
      
      Wt::log("info") << "Updating spectrum in database for '" << input->filename() << "'.";
      
      // Make a copy of the `UserFileInDbData`, so we can do the work of serializing things, while
      //  not holding a DB transaction
      UserFileInDbData newdata = *data;
      newdata.setFileData( meas, UserFileInDbData::sm_defaultSerializationFormat );
      
      try
      {
        DataBaseUtils::DbTransaction transaction( *m_sql );
        (*data.modify()) = newdata;
        fileDbEntry.modify()->userHasModified = modifiedSinceDecode;
        transaction.commit();
      }catch( std::exception &e )
      {
        try
        {
          DataBaseUtils::DbTransaction transaction( *m_sql );
          fileDbEntry.remove();
          transaction.commit();
        }catch(...){}
        
        {
          RecursiveLock lock( m_mutex );
          m_fileDbEntry.reset();
        }
        
        throw runtime_error( e.what() );
      }//try / catch
    }else
    {
      UserFileInDbData *dataptr = new UserFileInDbData();
      dataptr->fileInfo = fileDbEntry;
      try
      {
        dataptr->setFileData( meas, UserFileInDbData::sm_defaultSerializationFormat );
      }catch( std::exception &e )
      {
        delete dataptr;
        
        try
        {
          DataBaseUtils::DbTransaction transaction( *m_sql );
          fileDbEntry.remove();
          transaction.commit();
        }catch(...){}
        
        {
          RecursiveLock lock( m_mutex );
          m_fileDbEntry.reset();
        }

        throw runtime_error( e.what() );
      }//try / catch
      
      DataBaseUtils::DbTransaction transaction( *m_sql );
      data = m_sql->session()->add( dataptr );
      Wt::log("info") << "Adding '" << meas->filename() << "' spectrum to the database as id "
                      << data.id() << " to parent " << m_fileDbEntry.id();
      transaction.commit();
    }//if( !data )
    
    {
      DataBaseUtils::DbTransaction transaction( *m_sql );
      fileDbEntry.modify()->serializeTime = WDateTime::currentDateTime();
    
      //XXX - a Wt::Dbo::StaleObjectException may be thrown if modifying the same
      //      UserFileInDbData object in another session - but this is unlikely
      //      to happen, so we wont worry about it
      transaction.commit();
    }
  }//if( m_fileDbEntry )
  
  {
    DataBaseUtils::DbTransaction transaction( *m_sql );
    m_sql->session()->flush();  //needed so fileDbEntry.id() will be come not -1
    transaction.commit();
  }
  
  {
    RecursiveLock lock( m_mutex );
    m_fileDbEntry = fileDbEntry;
  }
  
  
  {
    std::shared_ptr<SpecMeas> ncmeas;
  
    {
      RecursiveLock lock( m_mutex );
      ncmeas = m_weakMeasurmentPtr.lock();
      if( !!ncmeas )
        ncmeas->reset_modified();
    }
  }
}//void SpectraFileHeader::saveToDatabase(...)


std::shared_ptr<SpecMeas> SpectraFileHeader::readFromDataBase() const
{
  long long int fileDbEntryID;
  
  {
    RecursiveLock lock( m_mutex );
    if( !m_fileDbEntry )
      throw runtime_error( "SpectraFileHeader::readFromDataBase(): "
                           "no database entry" );
    fileDbEntryID = m_fileDbEntry.id();
  }

  Wt::Dbo::ptr<UserFileInDb> fileDbEntry;
  
  {//Start interacting with DB
    DataBaseUtils::DbTransaction transaction( *m_sql );
    fileDbEntry = m_sql->session()->find<UserFileInDb>().where("id = ?").bind( fileDbEntryID );
    transaction.commit();
    
    if( !fileDbEntry )
      throw runtime_error( "SpectraFileHeader::readFromDataBase(): "
                           "no database entry" );
  }//End interacting with DB
  
  
  Dbo::ptr<UserFileInDbData> data;
  
  {//begin interaction with db
    DataBaseUtils::DbTransaction transaction( *m_sql );
  
    typedef Dbo::collection<Dbo::ptr<UserFileInDbData> > Files;
    if( fileDbEntry && fileDbEntry->filedata.size() )
      data = m_fileDbEntry->filedata.front();

    transaction.commit();
  }//end interaction with db
  
  std::shared_ptr<SpecMeas> meas = data->decodeSpectrum();

  if( meas )
    cerr << "Read spectrumFile from database" << endl;
  else
    cerr << "Failed to read spectrumFile from database" << endl;
  
  if( meas )
  {
    RecursiveLock lock( m_mutex );
    meas->reset_modified();
  }
  
  return meas;
}//std::shared_ptr<SpecMeas> readFromDataBase() const
#endif //#if( USE_DB_TO_STORE_SPECTRA )

void SpectraFileHeader::releaseCacheReference()
{
  m_cachedMeasurement.reset();
}//void releaseCacheReference()


void SpectraFileHeader::setKeepCachedInMemmorry( bool cache )
{
  m_keepCache = cache;
  if( !cache )
    releaseCacheReference();

  if( cache )
  {
    std::shared_ptr<SpecMeas> memObj = m_weakMeasurmentPtr.lock();
    if( memObj && !m_cachedMeasurement )
      m_cachedMeasurement = m_weakMeasurmentPtr.lock();
    else if( m_cachedMeasurement != memObj )
      m_weakMeasurmentPtr = m_cachedMeasurement;
  }//if( cache )
}//void setKeepCachedInMemmorry( bool cache )



std::shared_ptr<SpecMeas> SpectraFileHeader::initFile(
                                           const std::string &filename,
                                           const SpecUtils::ParserType parseType,
                                           std::string orig_file_ending ) throw()
{
  try
  {
    std::shared_ptr<SpecMeas> info = std::make_shared<SpecMeas>();

    const bool success = info->load_file( filename, parseType, orig_file_ending );

    if( success )
    {
      //Make sure the file we're opening doesn't have erroneous DB state indexes in it, which could
      //  cause us to overwrite a state.
      // I *think* all opening of spectrum files from filesystem go through here of `setFile(...)`.
      //Currently we are not writing this info to the SpecFile, so this call is not actually needed.
      info->clearAllDbStateId(); //JIV
      
      RecursiveLock lock( m_mutex );
      m_weakMeasurmentPtr = info;
      return info;
    }//if( success )

  }catch( const std::exception &e )
  {
    cerr << "SpectraFileHeader::initFile(...) caught exception:\n\t" << e.what() << endl;
  }catch(...)
  {
    cerr << "SpectraFileHeader::initFile(...) caught unknown exception!" << endl;
  }

  return std::shared_ptr<SpecMeas>();
}//std::shared_ptr<SpecMeas> initFile( const std::string &filename, ParserType parseType ) throw()


void SpectraFileHeader::errorSavingCallback( std::string fileLocation,
                                             std::shared_ptr<SpecMeas> meas ) const
{
  RecursiveLock lock( m_mutex );
  m_fileSystemLocation = "";
  cerr << "SpectraFileHeader::errorSavingCallback(...)\n\tThere was an error saving file to location: "
       << fileLocation << " - will attempt to cache it instead, but this is "
       << "no garuntee" << endl;

  //XXX - can we add a notification to the user here?

  //Lets cache this SpecMeas object ince we couldnt save it
  //  XXX - should check to see if meas is consistent with m_cachedMeasurement
  //        and m_weakMeasurmentPtr
  if( !meas )
    meas = m_weakMeasurmentPtr.lock();

  if( meas != m_cachedMeasurement )
    m_cachedMeasurement = meas;
  if( meas != m_weakMeasurmentPtr.lock() )
    m_weakMeasurmentPtr = meas;
}//void SpectraFileHeader::errorSavingCallback()




void SpectraFileHeader::saveToFileSystemImmediately( SpecMeas *meas ) const
{
  RecursiveLock lock( m_mutex );
  m_modifiedSinceDecode = meas->modified_since_decode();
  
  try
  {
    if( !meas )
      throw runtime_error( "You must pass in a SpecMeas object to save." );

    std::shared_ptr<SpecMeas> from_mem = measurementIfInMemory();
    if( from_mem && from_mem.get() != meas )
      throw runtime_error( "SpecMeas passed in is not the same as in memmorry" );

    if( !m_fileSystemLocation.empty() )
    {
      SpecUtils::remove_file( m_fileSystemLocation );
      m_fileSystemLocation = "";
    }//if( !m_fileSystemLocation.empty() )
    
    const string tempfile = SpecUtils::temp_file_name( m_displayName, InterSpecApp::tempDirectory() );
    m_fileSystemLocation = tempfile;
    
    {
#ifdef _WIN32
      const std::wstring wfilename = SpecUtils::convert_from_utf8_to_utf16(m_fileSystemLocation);
      ofstream output( wfilename.c_str(), ios::binary | ios::out );
#else
      ofstream output( m_fileSystemLocation.c_str(), ios::binary | ios::out );
#endif
      if( !output.is_open() )
        throw runtime_error( "Couldnt open file for writing: " + m_fileSystemLocation );
      if( !meas->write_2012_N42( output ) )
        throw runtime_error( "Failed to write 2012 N42 do to: " + m_fileSystemLocation );
    }

//#if( USE_DB_TO_STORE_SPECTRA )
//    if( shouldSaveToDb() )
//      saveToDatabaseFromTempFileWorker();
//#endif
  }catch( std::exception &e )
  {
    cerr << "SpectraFileHeader::saveToFileSystemImmediately(...)\n\tCaught: " << e.what() << endl;

    //XXX - should notify user of the error here?
  }//try / catch
}//void saveToFileSystemImmediately( SpecMeas *meas ) const



void SpectraFileHeader::saveToFileSystem( std::shared_ptr<SpecMeas> measurment ) const
{
  bool success = false;
  std::shared_ptr<SpecMeas> info;

  try
  {
    std::shared_ptr<SpecMeas> from_mem;
    
    {
      RecursiveLock lock( m_mutex );
      from_mem = measurementIfInMemory();
    }
    
    
    if( (from_mem && measurment) && (from_mem!=measurment) )
      throw runtime_error( "Measurement passed in is not same as currently in memory" );

    if( !measurment && !from_mem )
      from_mem = parseFile();

    info = from_mem ? from_mem : measurment;

    if( !info )
      throw runtime_error( "There is no SpecUtils::SpecFile anywhere to save" );

    {
      RecursiveLock lock( m_mutex );
      m_modifiedSinceDecode = info->modified_since_decode();
    }
    
    cerr << "In SpectraFileHeader::saveToFileSystem(...) for '"
         << info->filename() << "' m_modifiedSinceDecode="
         << m_modifiedSinceDecode << endl;

    try
    {
      RecursiveLock lock( m_mutex );
      if( !m_fileSystemLocation.empty() )
        SpecUtils::remove_file( m_fileSystemLocation );
      m_fileSystemLocation = "";
    }catch(...)
    {
    }

    const string tempfile = SpecUtils::temp_file_name( m_displayName, InterSpecApp::tempDirectory() );

    {
      RecursiveLock lock( m_mutex );
      m_fileSystemLocation = tempfile;
      
      m_userStateDbIndexes = info->dbUserStateIndexes();
    }
    
    success = true;
    
//    success = info->save_native_file( tempfile.generic_string() );
    boost::function<void()> error_callback = boost::bind( &SpectraFileHeader::errorSavingCallback, this, tempfile, info );
    SpecMeas::save2012N42FileInClientThread( info, tempfile, error_callback );
    
    
//#if( USE_DB_TO_STORE_SPECTRA )
//    if( m_app && shouldSaveToDb() )
//    {
//      boost::function<void(void)> worker
//                  = m_app->bind( boost::bind( &SpectraFileHeader::saveToDatabaseWorker, this, measurment ) );
//      WServer::instance()->post( m_app->sessionId(), worker );
//    }else if( m_candidateForSavingToDb )
//    {
//      cerr << "\n\nSpectraFileHeader::saveToFileSystem(...)\n\tSerious error: m_app is invalid"
//           << endl << endl;
//    }//if( m_app ) / else
//#endif
  }catch( std::exception &e )
  {
    cerr << "SpectraFileHeader::saveToFileSystem(...)\n\tCaught: " << e.what() << endl;

    //XXX - should notify user of the error here?

    success = false;
  } //try / catch

  if( !success )
  {
    cerr << "SpectraFileHeader::saveToFileSystem(...)\n\tWarning, couldnt save temp file "
         << " - am caching instead" << endl;
    if( info )
    {
      RecursiveLock lock( m_mutex );
      m_cachedMeasurement = info;
      m_weakMeasurmentPtr = info;
    }
  }//if( !success )
}//void saveToFileSystem()


struct SpectraHeaderMaker
{
  std::shared_ptr<const SpecMeas> m_info;
  vector<SpectraHeader>::iterator m_headerPos;  //expects m_headerPos through m_headerPos+m_samplenums->size()-1 to be allocated and assinable
  vector<int>::const_iterator m_sampleNumBegin;
  vector<int>::const_iterator m_sampleNumEnd;

  SpectraHeaderMaker( std::shared_ptr<const SpecMeas> info,
                      vector<int>::const_iterator sample_num_begin,
                      vector<int>::const_iterator sample_num_end,
                      vector<SpectraHeader>::iterator headerpos )
    : m_info( info ), m_headerPos( headerpos ),
      m_sampleNumBegin( sample_num_begin ), m_sampleNumEnd( sample_num_end )
  {
  }

  void operator()()
  {
    for( vector<int>::const_iterator iter = m_sampleNumBegin;
         iter != m_sampleNumEnd;
         ++iter )
    {
      const int sample = *iter;
      vector< std::shared_ptr<const SpecUtils::Measurement> > measurements
                                      = m_info->sample_measurements( sample );
      m_headerPos->init( measurements );
      ++m_headerPos;
    }//for( const int sample_number : sample_numbers )
  }//void operator()()
};//struct SpectraHeaderMaker


void SpectraFileHeader::setMeasurmentInfo( std::shared_ptr<SpecMeas> info )
{
  RecursiveLock lock( m_mutex );

  m_displayName = info->filename();

  m_uploadTime = WDateTime::currentDateTime();
  const set<int> sample_numbers = info->sample_numbers();
  const vector<int> &detector_numbers = info->detector_numbers();
  m_numSamples = static_cast<int>( sample_numbers.size() );
  m_isPassthrough = info->passthrough();

  m_totalLiveTime = info->gamma_live_time();
  m_totalRealTime = info->gamma_real_time();
  m_totalGammaCounts = info->gamma_count_sum();
  m_totalNeutronCounts = info->neutron_counts_sum();
  m_numDetectors = static_cast<int>( info->detector_numbers().size() );
  m_hasNeutronDetector = !info->neutron_detector_names().empty();
  
  m_uuid = info->uuid();
  if( static_cast<int>(m_uuid.length()) > UserFileInDb::sm_maxUuidLength )
    m_uuid = m_uuid.substr( 0, UserFileInDb::sm_maxUuidLength );
  
  m_modifiedSinceDecode = info->modified_since_decode();
  
  // Incase we are updating an already filled out entity, clear out some variables we will fill in
  m_samples.clear();
  m_hasNeutronDetector = false;
  m_spectrumTime = Wt::WDateTime();
  
  for( const int sample : sample_numbers )
  {
    for( const int detector : detector_numbers )
    {
      const std::shared_ptr<const SpecUtils::Measurement> meas = info->measurement( sample, detector );
      //This next line is not necessary, as long as we properly initialized 'info'
      m_hasNeutronDetector |= (meas && meas->contained_neutron());
      if( meas && (sample == (*sample_numbers.begin())) )
        m_spectrumTime = WDateTime::fromPosixTime( to_ptime( meas->start_time() ) );
    }//for( const int detector : detector_numbers )
  }//for( const int sample : sample_numbers )

  const size_t nsamplenums = sample_numbers.size();


  if( nsamplenums < 250 ) //250 is arbitrary
  {
    for( const int sample_number : sample_numbers )
    {
      vector< std::shared_ptr<const SpecUtils::Measurement> > measurements
                                      = info->sample_measurements( sample_number );
      m_samples.push_back( SpectraHeader( measurements ) );
    }//for( const int sample_number : sample_numbers )
  }else
  {
    //We'll do this multithreaded
//    const vector<int> sample_num_vec( sample_numbers.begin(), sample_numbers.end() );
    vector<int> sample_num_vec( nsamplenums );
    std::copy( sample_numbers.begin(), sample_numbers.end(),
               sample_num_vec.begin() );

    const int nthread = SpecUtilsAsync::num_logical_cpu_cores();
    const size_t meas_per_thread = nsamplenums / nthread;

//    vector<SpectraHeaderMaker> workers;
    SpecUtilsAsync::ThreadPool pool;
    m_samples.resize( nsamplenums );

    for( size_t pos = 0; pos < nsamplenums; pos += meas_per_thread )
    {
      vector<int>::const_iterator start, end;
      start = sample_num_vec.begin() + pos;

      if( (pos + meas_per_thread) <= nsamplenums )
        end = start + meas_per_thread;
      else
        end = sample_num_vec.begin() + nsamplenums;

      pool.post( SpectraHeaderMaker( info, start, end, m_samples.begin()+pos ) );
//      workers.push_back(
//                SpectraHeaderMaker( info, start, end, m_samples.begin()+pos ) );
    }//for( size_t pos = 0; pos < nstart; pos += meas_per_thread )

//    SpecUtils::do_asyncronous_work( workers, false );
    pool.join();
  }//if( sample_numbers.size() < 250 ) / else

  if( m_keepCache )
    m_cachedMeasurement = info;

  m_weakMeasurmentPtr = info;


  if( m_fileSystemLocation.size() )
  {
    SpecUtils::remove_file( m_fileSystemLocation );
    m_fileSystemLocation = "";
  }//if( m_fileSystemLocation.size() )


  if( info )
  {
    if( m_aboutToBeDeletedConnection.connected() )
      m_aboutToBeDeletedConnection.disconnect();
    m_aboutToBeDeletedConnection
             = info->aboutToBeDeleted().connect( boost::bind(
                                &SpectraFileHeader::saveToFileSystemImmediately,
                                                              this, info.get() ) );
  }
}//void setMeasurmentInfo( std::shared_ptr<SpecMeas> );


void SpectraFileHeader::setFile( const std::string &displayFileName, std::shared_ptr<SpecMeas> info )
{
  if( !info )
    throw runtime_error( "Not a valid spectrum file." );
  
  RecursiveLock lock( m_mutex );
    
  m_weakMeasurmentPtr = info;
  
  info->set_filename( displayFileName );
  
  info->reset_modified();
  info->reset_modified_since_decode();
  
  //Make sure the file we're opening doesn't have erroneous DB state indexes in it, which could
  //  cause us to overwrite a state.
  // I *think* all opening of spectrum files from filesystem go through here of `initFile(...)`.
  //Currently we are not writing this info to the SpecFile, so this call is not actually needed.
  info->clearAllDbStateId(); //JIC
  
  m_modifiedSinceDecode = false;
  
  setMeasurmentInfo( info );
}//void setFile( const std::string &displayFileName, std::shared_ptr<SpecMeas> meas );


std::shared_ptr<SpecMeas> SpectraFileHeader::setFile(
                                   const std::string &displayFileName,
                                   const std::string &filename,
                                   SpecUtils::ParserType parseType )
{
  Wt::log("debug") << "SpectraFileHeader::setFile('" << displayFileName << "', '" << filename << "')";
  try
  {
    if( !SpecUtils::is_file(filename) )
      throw runtime_error( "" );
  }catch(...)
  {
    Wt::log("error") << "Could not access the file '"
         << displayFileName << "', located at '" << filename << "'";
    
    throw runtime_error( "Could not access file '" + displayFileName + "'" );
  }//try / catch
  
  string orig_file_ending = SpecUtils::file_extension(displayFileName);
  shared_ptr<SpecMeas> info = initFile( filename, parseType, orig_file_ending );
  if( !info )
  {
    throw std::runtime_error( "Could not open '" + displayFileName
                             + "' with any of the available decoders, sorry." );
  }
  
  setFile( displayFileName, info );

  return info;
}//setFile(...)


std::shared_ptr<SpecMeas> SpectraFileHeader::parseFile() const
{
  string filesystemlocation;
 
  {//begin mutex protected code
    RecursiveLock lock( m_mutex );
    if( m_cachedMeasurement )
      return m_cachedMeasurement;
  
    std::shared_ptr<SpecMeas> memObj = m_weakMeasurmentPtr.lock();

    if( memObj )
      return memObj;
  
    filesystemlocation = m_fileSystemLocation;
  }//end mutex protected code
  
  cerr << "In parseFile() and not using weak ptr" << endl;

  bool success = false;
  auto info = std::make_shared<SpecMeas>();

  success = info->load_N42_file( filesystemlocation );

  if( !success )
  {
    stringstream msg;
    msg << "SpectraFileHeader::parseFile():\n\tfor some reason I couldnt parse "
        << m_fileSystemLocation << " using method using the expected native "
        << "binary method - this is odd and probably an internal logic error, "
        << "please report this error";
            cout<<"Exit 3 LOCK parseFile"<<endl;
    throw std::runtime_error( msg.str() );
  }//if( !success )  
  
  {//begin mutex protected code
    RecursiveLock lock( m_mutex );
    if( m_keepCache )
      m_cachedMeasurement = info;
    m_weakMeasurmentPtr = info;
    
    // TODO: Should consider using an on error callback for following connections
    if( info )
    {
      // Restore mappings of sample numbers in this file, to UserState indexes in the database
      info->clearAllDbStateId(); //JIC, but not needed
      for( const auto &sn_index : m_userStateDbIndexes )
        info->setDbStateId( sn_index.second, sn_index.first );
      m_userStateDbIndexes.clear();
      
      if( m_aboutToBeDeletedConnection.connected() )
        m_aboutToBeDeletedConnection.disconnect();
      m_aboutToBeDeletedConnection = info->aboutToBeDeleted().connect(
                  boost::bind( &SpectraFileHeader::saveToFileSystemImmediately,
                               this, info.get() ) );
    }//if( info )

  }//end mutex protected code
  
  
  return info;
}//std::shared_ptr<SpecMeas> parseFile() const


std::shared_ptr<SpecMeas> SpectraFileHeader::measurementIfInMemory() const
{
  return m_weakMeasurmentPtr.lock();
}//std::shared_ptr<SpecMeas> measurementIfInMemory()


int SpectraFileHeader::numSamples() const
{
  return m_numSamples;
}//int numSamples() const


Wt::WString SpectraFileHeader::displayName() const
{
  return m_displayName;
}


const Wt::WDateTime &SpectraFileHeader::uploadTime() const
{
  return m_uploadTime;
}


float SpectraFileHeader::totalLiveTime() const
{
  return m_totalLiveTime;
}


float SpectraFileHeader::totalRealTime() const
{
  return m_totalRealTime;
}


float SpectraFileHeader::totalGammaCounts() const
{
  return m_totalGammaCounts;
}


float SpectraFileHeader::totalNeutronCounts() const
{
  return m_totalNeutronCounts;
}


int SpectraFileHeader::numDetectors() const
{
  return m_numDetectors;
}


bool SpectraFileHeader::passthrough() const
{
  return m_isPassthrough;
}


bool SpectraFileHeader::hasNeutronDetector() const
{
  return m_hasNeutronDetector;
}


const Wt::WDateTime &SpectraFileHeader::spectrumTime() const
{
  return m_spectrumTime;
}


SpectraFileModel::SpectraFileModel( Wt::WObject *parent )
  : WAbstractItemModel( parent )
{
  //m_spectra
}


SpectraFileModel::~SpectraFileModel()
{
}


int SpectraFileModel::columnCount( const Wt::WModelIndex &parent ) const
{
  const Level indexLevel = level(parent);
  switch( indexLevel )
  {
    case FileHeaderLevel:
      return static_cast<int>( NumDisplayFields );
    case SampleLevel:
      return 0;
    case InvalidLevel:
      return static_cast<int>( NumDisplayFields );
  }//switch( level(index) )

  return 0;
//  if( parent.internalPointer() )
//    return 0;
//  return static_cast<int>( NumDisplayFields );
}//int columnCount(...)


int SpectraFileModel::rowCount( const Wt::WModelIndex &parent ) const
{  
  const Level indexLevel = level(parent);
  switch( indexLevel )
  {
    case FileHeaderLevel:
    {
      const int row = parent.row();
      if( row >= static_cast<int>( m_spectra.size() ) )
        return 0;

      std::shared_ptr<const SpectraFileHeader> header = m_spectra[row];
      if( !header )
        return 0;

      const int nspectra = header->numSamples();
      if( nspectra < 2 )
        return 0;

      return nspectra;
    }//case SpectraFileModel::FileHeaderLevel:

    case SampleLevel:
      return 0;
    case InvalidLevel:
      return static_cast<int>( m_spectra.size() );
  }//switch( level(index) )

//  if( parent.internalPointer() )
//    return 0;
//  if( !parent.isValid() )
//    return static_cast<int>( m_spectra.size() );
//  const int row = parent.row();
//  if( row >= static_cast<int>( m_spectra.size() ) )
//    return 0;
//  std::shared_ptr<SpectraFileHeader> header = m_spectra[row];
//  if( !header )
//    return 0;
//  const int nspectra = header->numSamples();
//  if( nspectra < 2 )
//    return 0;
//  return nspectra;
  return 0;
}//int rowCount(...)


boost::any SpectraFileModel::data( const Wt::WModelIndex &index,
                                   int role ) const
{
  if( role != Wt::DisplayRole )
  {
//    DisplayRole = 0,      //!< Role for textual representation
//    DecorationRole = 1,   //!< Role for the url of an icon
//    EditRole = 2,         //!< Role for the edited value
//    StyleClassRole = 3,   //!< Role for the style class
//    CheckStateRole = 4,
//    ToolTipRole = 5,      //!< Role for a tooltip
//    InternalPathRole = 6, //!< Role for an internal path activated when clicked
//    UrlRole = 7,          //!< Role for a url activated when clicked
//    LevelRole = 8,        //!< Level in aggregation, for header data.
//    MarkerPenColorRole = 16,  //!< Marker pen color (for Chart::WCartesianChart)
//    MarkerBrushColorRole = 17,//!< Marker brush color (for Chart::WCartesianChart)
//    cerr << "\trole != Wt::DisplayRole: " << role << endl;
    return boost::any();
  }//else cerr << "Asked for display role" << endl;

  const Level indexLevel = level(index);

  if( indexLevel == InvalidLevel )
  {
    cerr << "\t!index.isValid()" << endl;
    return boost::any();
  }//if( indexLevel == InvalidLevel )

  const int row = index.row();
  stringstream strm;
  strm.setf( ios::fixed );

  if( indexLevel == SampleLevel )
  {
    //the parent index will have a null pointer, but the index.internalPointer()
    //  will point to its parent
    SpectraFileHeader *fileHeader = (SpectraFileHeader *)index.internalPointer();
    if( !fileHeader )
      throw std::runtime_error( "SpectraFileModel::data(...): Serious logic"
                                " error, valid WModelIndex should always have a"
                                " valid internal pointer when it has a parent" );

    //now make sure this is still a valid pointer to a parent
    bool found = false;
    for( size_t i = 0; !found && i < m_spectra.size(); ++i )
      found = (found || (m_spectra[i].get()==fileHeader) );

    if( !found )
    {
      cerr << "SpectraFileModel::data(...):\n\tfailed to find spectra, returning empty answer"
           << endl;
      return boost::any();
    }//if( !found )

    if( row >= static_cast<int>( fileHeader->m_samples.size() ) )
    {
      cerr << "\trow >= static_cast<int>( fileHeader->m_samples.size() ) row="
           << row << ", fileHeader->m_samples.size()="
           << fileHeader->m_samples.size() << endl;
      return boost::any();
    }

    const SpectraHeader &spectra_header = fileHeader->m_samples[row];

    const DisplayFields field = DisplayFields( index.column() );
    
    switch( field )
    {
      case kDisplayName:
        strm << spectra_header.remarks;
      break;
        
      case kNumDetectors:
        strm << spectra_header.detector_numbers_.size();
      break;
      case kUploadTime:
      break;
        
      case kNumMeasurements:
        strm << (row+1);
      break;
        
      case kLiveTime:
        if( is_close_to_int(spectra_header.live_time) )
          strm << setprecision(1) << spectra_header.live_time;
        else
          strm << setprecision(2) << spectra_header.live_time;
      break;
        
      case kRealTime:
        if( is_close_to_int(spectra_header.real_time) )
          strm << setprecision(1) << spectra_header.real_time;
        else
          strm << setprecision(2) << spectra_header.real_time;
      break;
        
      case kGammaCounts:
        if( is_close_to_int(spectra_header.gamma_counts_) )
          strm << setprecision(0) << spectra_header.gamma_counts_;
        else
          strm << setprecision(2) << spectra_header.gamma_counts_;
      break;
        
      case kNeutronCounts:
        if( spectra_header.contained_neutron_ )
        {
          if( is_close_to_int(spectra_header.neutron_counts_) )
            strm << setprecision(0) << spectra_header.neutron_counts_;
          else
            strm << setprecision(2) << spectra_header.neutron_counts_;
        }//if( spectra_header.contained_neutron_ )
      break;
        
      case kSpectrumTime:
        strm << spectra_header.start_time.toString( DATE_TIME_FORMAT_STR );
      break;
        
      case NumDisplayFields:
      break;
    };//switch( field )
    
    string str = strm.str();
    
    //get rid of malicious stuff here
    
    if( !str.empty() )
      return boost::any( WString(str) );
    return boost::any();
  }//if( index.parent().isValid() )

//  indexLevel == FileHeaderLevel

  if( row >= static_cast<int>( m_spectra.size() ) )
  {
    cerr << "\trow >= static_cast<int>( m_spectra.size() )" << endl;
    return boost::any();
  }

  std::shared_ptr<SpectraFileHeader> header = m_spectra[row];
  
  const DisplayFields field = DisplayFields( index.column() );
  
  switch( field )
  {
    case kDisplayName:
      strm << header->displayName().toUTF8();  //return boost::any( header->displayName() );
    break;
      
    case kUploadTime:
      strm << header->uploadTime().toString( DATE_TIME_FORMAT_STR );
    break;
      
    case kNumMeasurements:
      strm << header->numSamples();
    break;
      
    case kNumDetectors:
      strm << header->numDetectors();
    break;
      
    case kLiveTime:
      if( is_close_to_int(header->totalLiveTime()) )
        strm << setprecision(1) << header->totalLiveTime();
      else
        strm << setprecision(2) << header->totalLiveTime();
    break;
      
    case kRealTime:
      if( is_close_to_int(header->totalRealTime()) )
        strm << setprecision(1) << header->totalRealTime();
      else
        strm << setprecision(2) << header->totalRealTime();
    break;
      
    case kGammaCounts:
      if( is_close_to_int( header->totalGammaCounts() ) )
        strm << setprecision(0) << header->totalGammaCounts();
      else
        strm << setprecision(2) << header->totalGammaCounts();
    break;
      
    case kNeutronCounts:
      if( header->hasNeutronDetector() )
      {
        if( is_close_to_int( header->totalNeutronCounts() ) )
          strm << setprecision(0) << header->totalNeutronCounts();
        else
          strm << setprecision(2) << header->totalNeutronCounts();
      }//if( header->hasNeutronDetector() )
    break;
      
    case kSpectrumTime:
      strm << header->spectrumTime().toString( DATE_TIME_FORMAT_STR );
    break;
      
    case NumDisplayFields:
    break;
  };//switch( field )

  string str = strm.str();
  
  //get rid of malicious stuff here
  
  if( !str.empty() )
    return boost::any( WString(str) );
  return boost::any();
}//boost::any SpectraFileModel::data(...)


std::shared_ptr<const SpectraFileHeader> SpectraFileModel::fileHeader( int row ) const
{
  if( (row >= 0) && (row < static_cast<int>(m_spectra.size())) )
    return m_spectra[row];
  return nullptr;
}


std::shared_ptr<SpectraFileHeader> SpectraFileModel::fileHeader( int row )
{
  if( (row >= 0) && (row < static_cast<int>(m_spectra.size())) )
    return m_spectra[row];
  return nullptr;
}

std::shared_ptr<SpectraFileHeader> SpectraFileModel::fileHeader(
                        std::shared_ptr<const SpecUtils::SpecFile> measinfo )
{
  WModelIndex in = index( measinfo );
  if( !in.isValid() )
    return std::shared_ptr<SpectraFileHeader>();
  return fileHeader( in.row() );
}//fileHeader(...)


#if( USE_DB_TO_STORE_SPECTRA )
Wt::Dbo::ptr<UserFileInDb> SpectraFileModel::dbEntry(
                        std::shared_ptr<const SpecUtils::SpecFile> measinfo )
{
  std::shared_ptr<SpectraFileHeader> head = fileHeader( measinfo );
  if( !head )
    return Dbo::ptr<UserFileInDb>();
  return head->dbEntry();
}//dbEntry(...)
#endif


Wt::WModelIndex SpectraFileModel::index( std::shared_ptr<SpectraFileHeader> header ) const
{
  for( size_t i = 0; i < m_spectra.size(); ++i )
  {
    const std::shared_ptr<SpectraFileHeader> &thisHeader = m_spectra[i];
    if( thisHeader == header )
      return index( static_cast<int>(i), 0, WModelIndex() );
  }//for spectra

  return WModelIndex();
}//Wt::WModelIndex index( std::shared_ptr<SpectraFileHeader> header ) const;


Wt::WModelIndex SpectraFileModel::index( std::shared_ptr<const SpecUtils::SpecFile> measinfo ) const
{
  if( !measinfo )
    return WModelIndex();
  
  for( size_t i = 0; i < m_spectra.size(); ++i )
  {
    const std::shared_ptr<SpectraFileHeader> &thisHeader = m_spectra[i];
    if( !thisHeader )
      continue;
    std::shared_ptr<SpecMeas> meas;
    meas = thisHeader->measurementIfInMemory();

    if( meas == measinfo )
      return index( static_cast<int>(i), 0, WModelIndex() );
  }//for spectra

  return WModelIndex();
}//Wt::WModelIndex SpectraFileModel::index( std::shared_ptr<SpecMeas> measinfo ) const


int SpectraFileModel::addRow( std::shared_ptr<SpectraFileHeader> fileInfo )
{
  beginInsertRows( WModelIndex(), 
	               static_cast<int>(m_spectra.size()), 
				   static_cast<int>(m_spectra.size()) );
  m_spectra.push_back( fileInfo );
  endInsertRows();

  return static_cast<int>( m_spectra.size()-1 );
}//int addRow(...)


SpectraFileModel::Level SpectraFileModel::level( const WModelIndex &index ) const
{

  if( !index.isValid() )
    return InvalidLevel;

  if( !index.parent().isValid() )
    return FileHeaderLevel;

  return SampleLevel;
}//Level level( const WModelIndex &index ) const


bool SpectraFileModel::removeRows( int row, int count, const WModelIndex &parent )
{
  if( parent.isValid() )
  {
    cerr << "SpectraFileModel::removeRows(...):\n\tyou can only remove SpectraFileHeader rows"
         << endl;
    return false;
  }//if( parent.isValid() )

  if( row >= rowCount() )
    return false;

  std::vector< std::shared_ptr<SpectraFileHeader> >::iterator startRow;
  startRow = m_spectra.begin() + row;

  count = (int)min( static_cast<ptrdiff_t>(count), m_spectra.end() - startRow);

  if( count <= 0 )
  {
    cerr << "SpectraFileModel::removeRows(...)\n\t: error, trying to remove zero rows" << endl;
    return false;
  }//if( count <= 0 )

  const int lastRow = row + count - 1;

  /*
  //I wouldn't have guessed we would have to do this, but if we don't also
  //  separately notify SpecMeasManager::m_treeView that we are
  //  deleting the children, all of the WModelIndex's pointing to the
  //  children become orphaned.
  for( int i = row; i <= lastRow; ++i )
  {
    WModelIndex thisIndex = index( i, 0 );
    std::shared_ptr<SpectraFileHeader> thisHeader = m_spectra.at(i);
    const int nsamples = rowCount( thisIndex );

    beginRemoveRows( thisIndex, 0, nsamples );
    thisHeader->m_numSamples = 0;
    thisHeader->m_samples.clear();
    endRemoveRows();
  }//for( int i = row; i <= lastRow; ++i )
*/
  //now actually remove this file
  beginRemoveRows( parent, row, lastRow );
  m_spectra.erase( startRow, startRow+count );
  endRemoveRows();

  return true;
}//bool removeRows(...)


WModelIndex SpectraFileModel::parent( const WModelIndex &index ) const
{
  SpectraFileHeader *parentPtr = (SpectraFileHeader *)index.internalPointer();
  if( !parentPtr )
    return WModelIndex();

  for( size_t row = 0; row < m_spectra.size(); ++row )
    if( parentPtr == m_spectra[row].get() )
      return createIndex( static_cast<int>(row), 0, (void *)NULL );

  cerr << "SpectraFileModel::parent(...):\n\tfailed to find parent for " << parentPtr << endl;

  return WModelIndex();
}//WModelIndex parent( const WModelIndex &index ) const


WModelIndex SpectraFileModel::index( int row, int column,
                                     const Wt::WModelIndex &parent ) const
{
  if( row<0 || column<0 )
    return WModelIndex();

  if( !parent.isValid() )
    return createIndex( row, column, (void *)0 );

  if( column >= NumDisplayFields || (column<0) )
    return WModelIndex();

  int parentRow = parent.row();
  if( parentRow >= static_cast<int>(m_spectra.size()) )
    return WModelIndex();

  return createIndex( row, column, m_spectra[parentRow].get() );
}//WModelIndex index(...)


boost::any SpectraFileModel::headerData( int section, Orientation orientation, int role ) const
{
  if( orientation == Horizontal && role==LevelRole )
    return 0;
  
  if( role != Wt::DisplayRole  )
    return boost::any();

  if( orientation == Horizontal )
  {
    // Make sure it's within the proper bounds
    if( section < 0 || section >= columnCount() )
      return boost::any();

    const DisplayFields field = DisplayFields( section );
    switch( field )
    {
      case kDisplayName:     return boost::any( WString("File") );
      case kNumMeasurements: return boost::any( WString("N-Samples") );
      case kLiveTime:        return boost::any( WString("Live Time (s)") );
      case kRealTime:        return boost::any( WString("Real Time (s)") );
      case kGammaCounts:     return boost::any( WString("Gam. Count") );
      case kNeutronCounts:   return boost::any( WString("Neut. Count") );
      case kSpectrumTime:    return boost::any( WString("Time Taken") );
      case kNumDetectors:    return boost::any( WString("N-Dets.") );
      case kUploadTime:      return boost::any( WString("Loaded") );
      case NumDisplayFields: break;
    };//switch( field )

  }else if( orientation == Vertical )
  {
    if( section < 0 || section >= rowCount() )
      return boost::any();
    return boost::any();
  }//if( orientation == Horizontal ) / else Vertical

  return boost::any();
}//any SpectraFileModel::headerData(...)


struct SortSpectraFileModel
{
  SortSpectraFileModel( SpectraFileModel::DisplayFields column, SortOrder order )
    :m_order(order), m_column( column )
  {
  }

  bool operator()( const std::shared_ptr<SpectraFileHeader> &lhs,
                   const std::shared_ptr<SpectraFileHeader> &rhs ) const
  {
    if( m_order == AscendingOrder )
    {
      switch( m_column )
      {
        case SpectraFileModel::kDisplayName:     return (lhs->displayName().toUTF8() < rhs->displayName().toUTF8());
        case SpectraFileModel::kNumMeasurements: return (lhs->numSamples() < rhs->numSamples());
        case SpectraFileModel::kNumDetectors:    return (lhs->numDetectors() < rhs->numDetectors());
        case SpectraFileModel::kLiveTime:        return (lhs->totalLiveTime() < rhs->totalLiveTime());
        case SpectraFileModel::kRealTime:        return (lhs->totalRealTime() < rhs->totalRealTime());
        case SpectraFileModel::kGammaCounts:     return (lhs->totalGammaCounts() < rhs->totalGammaCounts());
        case SpectraFileModel::kNeutronCounts:   return (lhs->totalNeutronCounts() < rhs->totalNeutronCounts());
        case SpectraFileModel::kSpectrumTime:    return (lhs->spectrumTime() < rhs->spectrumTime());
        case SpectraFileModel::kUploadTime:      return (lhs->uploadTime() < rhs->uploadTime());
        case SpectraFileModel::NumDisplayFields: return false;
      }//switch( m_column )
    }//if( order == AscendingOrder )

    switch( m_column )
    {
      case SpectraFileModel::kDisplayName:     return (lhs->displayName().toUTF8() >= rhs->displayName().toUTF8());
      case SpectraFileModel::kNumMeasurements: return (lhs->numSamples() >= rhs->numSamples());
      case SpectraFileModel::kNumDetectors:    return (lhs->numDetectors() >= rhs->numDetectors());
      case SpectraFileModel::kLiveTime:        return (lhs->totalLiveTime() >= rhs->totalLiveTime());
      case SpectraFileModel::kRealTime:        return (lhs->totalRealTime() >= rhs->totalRealTime());
      case SpectraFileModel::kGammaCounts:     return (lhs->totalGammaCounts() >= rhs->totalGammaCounts());
      case SpectraFileModel::kNeutronCounts:   return (lhs->totalNeutronCounts() >= rhs->totalNeutronCounts());
      case SpectraFileModel::kSpectrumTime:    return (lhs->spectrumTime() >= rhs->spectrumTime());
      case SpectraFileModel::kUploadTime:      return (lhs->uploadTime() >= rhs->uploadTime());
      case SpectraFileModel::NumDisplayFields: return true;
    }//switch( m_column )

    //shouldnt ever get here
    return false;
  }//bool operator()

  SortOrder m_order;
  SpectraFileModel::DisplayFields m_column;
};//struct SortSpectraFileModel


void SpectraFileModel::sort( int column, SortOrder order )
{
  layoutAboutToBeChanged().emit();
  std::sort( m_spectra.begin(), m_spectra.end(),
             SortSpectraFileModel(SpectraFileModel::DisplayFields(column),
                                  order) );
  layoutChanged().emit();
}//void sort(...);


void *SpectraFileModel::toRawIndex( const Wt::WModelIndex &index ) const
{
  const Level indLevel = level( index );
  if( indLevel == InvalidLevel )
    return NULL;

  const int row = index.row();
  const DisplayFields column = DisplayFields( index.column() );


  if( indLevel == FileHeaderLevel )
  {
    if( static_cast<size_t>(row) >= m_spectra.size() )
      return NULL;

    const SpectraFileHeader &header = *(m_spectra[row].get());
    switch( column )
    {
      case kDisplayName:     return (void *)&(header.m_displayName);
      case kNumMeasurements:  return (void *)&(header.m_numSamples);
      case kNumDetectors:    return (void *)&(header.m_numDetectors);
      case kLiveTime:        return (void *)&(header.m_totalLiveTime);
      case kRealTime:        return (void *)&(header.m_totalRealTime);
      case kGammaCounts:     return (void *)&(header.m_totalGammaCounts);
      case kNeutronCounts:   return (void *)&(header.m_totalNeutronCounts);
      case kSpectrumTime:    return (void *)&(header.m_spectrumTime);
      case kUploadTime:      return (void *)&(header.m_uploadTime);
      case NumDisplayFields: return NULL;
    }//switch( column )
  }//if( !index.parent().isValid() )

  //
//  indLevel == SampleLevel

  try
  {
    const SpectraFileHeader &header = *(m_spectra.at(index.parent().row()).get());
    if( row > 1 )
      return NULL;
    if( static_cast<size_t>(row) >= header.m_samples.size() )
      return NULL;

    const SpectraHeader &measurement = header.m_samples[row];

    switch( column )
    {
      case kDisplayName:     case kUploadTime:
      case kNumMeasurements:  case kNumDetectors:
        return NULL;

      case kLiveTime:        return (void *)&(measurement.live_time);
      case kRealTime:        return (void *)&(measurement.real_time);
      case kGammaCounts:     return (void *)&(measurement.gamma_counts_);
      case kNeutronCounts:   return (void *)&(measurement.neutron_counts_);
      case kSpectrumTime:    return (void *)&(measurement.start_time);
      case NumDisplayFields: return NULL;
    }//switch( column )
  }catch(...){ return NULL; }

  cerr << "SpectraFileModel::toRawIndex(...):\n\tshouldnt have got here" << endl;

  return NULL;
}//void *toRawIndex( const Wt::WModelIndex &index ) const


WModelIndex SpectraFileModel::fromRawIndex( void *rawIndex ) const
{
  for( int row = 0; row < static_cast<int>(m_spectra.size()); ++row )
  {
    const SpectraFileHeader &header = *(m_spectra[row].get());
    if( rawIndex == &(header.m_displayName) )
      return index( row, kDisplayName );
    if( rawIndex == &(header.m_uploadTime) )
      return index( row, kUploadTime );
    if( rawIndex == &(header.m_numSamples) )
      return index( row, kNumMeasurements );
    if( rawIndex == &(header.m_numDetectors) )
      return index( row, kNumDetectors );
    if( rawIndex == &(header.m_totalLiveTime) )
      return index( row, kLiveTime );
    if( rawIndex == &(header.m_totalRealTime) )
      return index( row, kRealTime );
    if( rawIndex == &(header.m_totalGammaCounts) )
      return index( row, kGammaCounts );
    if( rawIndex == &(header.m_totalNeutronCounts) )
      return index( row, kNeutronCounts );
    if( rawIndex == &(header.m_spectrumTime) )
      return index( row, kSpectrumTime );
  }//for( loop over rows to check...)

  for( size_t i = 0; i < m_spectra.size(); ++i )
  {
    const SpectraFileHeader &spectra = *(m_spectra[i].get());
    const int nsamples = static_cast<int>( spectra.m_samples.size() );

    for( int row = 0; row < nsamples; ++row )
    {
      const SpectraHeader &header = spectra.m_samples[row];
      if( rawIndex == &(header.live_time) )
        return createIndex( 0, kLiveTime,      (void *)&spectra );
      if( rawIndex == &(header.real_time) )
        return createIndex( 0, kRealTime,      (void *)&spectra );
      if( rawIndex == &(header.gamma_counts_) )
        return createIndex( 0, kGammaCounts,   (void *)&spectra );
      if( rawIndex == &(header.neutron_counts_) )
        return createIndex( 0, kNeutronCounts, (void *)&spectra );
      if( rawIndex == &(header.start_time) )
        return createIndex( 0, kSpectrumTime,  (void *)&spectra );
    }//for( int row = 0; row < nmeas; ++row )
  }//for( loop over rows to check...)

  cerr << "Failed to find WModelIndex for void*=" << rawIndex << endl;

  return WModelIndex();
}//WModelIndex fromRawIndex( void *rawIndex ) const


Wt::WFlags<Wt::ItemFlag> SpectraFileModel::flags( const WModelIndex &/*index*/ ) const
{
  return ItemIsSelectable;
}

