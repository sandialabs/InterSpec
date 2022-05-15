#ifndef SpectraFileModel_h
#define SpectraFileModel_h
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

#include <set>
#include <mutex>
#include <deque>
#include <memory>
#include <vector>
#include <string>
#include <sstream>

#include <boost/any.hpp>

#include <Wt/WString>
#include <Wt/Dbo/Dbo>
#include <Wt/Dbo/ptr>
#include <Wt/WDateTime>
#include <Wt/WResource>
#include <Wt/WDateTime>
#include <Wt/WModelIndex>
#include <Wt/Dbo/Session>
#include <Wt/Dbo/WtSqlTraits>
#include <Wt/WContainerWidget>
#include <Wt/WAbstractItemModel>

#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/SpecMeasManager.h"


//The basic goal of the classes in this file is to a way to manage user
//  uploaded files, with a way to manage them on disk, display summaries of them
//  on screen, and load them from disk into SpecUtils::SpecFile objects if they have
//  been flushed out of memorry.


namespace Wt
{
  class JSlot;
  class WObject;
  class WTreeView;
  class WFileUpload;
  class WApplication;
} // namespace Wt

namespace SpecUtils{ class SpecFile; }
namespace SpecUtils{ class Measurement; }
namespace SpecUtils{ enum class SourceType : int; }
namespace SpecUtils{ enum class SpectrumType : int; }
namespace SpecUtils{ enum class SaveSpectrumAsType : int; }

class PopupDiv;
class PopupDivMenu;
class InterSpec;
class SpecMeasManager;
class PopupDivMenuItem;
class SpectraFileHeader;



//SpectraHeader should be renamed MeasurementHeader
class SpectraHeader
{
  //Holds information about a single spectra for display purposes

public:
  SpectraHeader();
  SpectraHeader( const std::vector<std::shared_ptr<const SpecUtils::Measurement>> &sample_measurements );
  virtual ~SpectraHeader();

  void init( const std::vector<std::shared_ptr<const SpecUtils::Measurement>> &sample_measurements );

  float live_time;
  float real_time;
  bool contained_neutron_;
  int  sample_number;
  double gamma_counts_;
  double neutron_counts_;
  SpecUtils::SourceType spectra_type;
  Wt::WString speed_;
  std::vector<Wt::WString> detector_names;
  std::vector<int> detector_numbers_;
  std::string  remarks;
  Wt::WDateTime start_time;
  bool is_derived_data;
};//struct SpectraHeader


//SpectraFileHeader should be renamed SpecFileHeader
class SpectraFileHeader
{
  //Holds information about a spectra file, as well as keeps a copy of the
  //  SpecMeas in a temporary file if SpecMeasManager::sm_maxTempCacheSize has
  //  been exceeded causing the SpecMeas object to be removed from memory.
  //  If the user prefernces want it, on the destruction of the
  //  SpectraFileHeader, the SpecMeas will be placed into the database (same
  //  one as m_user is in) for access in subsequent sessions.
  //  (note that the database isnt just always used instread of temporary
  //   files, for performance - although I have yet to actually properly
  
  //
  //XXX - should convert shared_ptr<SpecMeas> to
  //      shared_ptr<const SpecUtils::SpecFile> wherever possible, as currently
  //      all caching mechanisms implicitly rely on SpecUtils::SpecFile not being
  //      changed out of this class
public:
  SpectraFileHeader( Wt::Dbo::ptr<InterSpecUser> user,
                     bool keepInMemmory,
                     InterSpec *viewer );

  virtual ~SpectraFileHeader() noexcept(true);

  //setFile thows std::runtime_error(..) on failure.  File passed in does
  //  not need to persist on the file system past calling this function.
  std::shared_ptr<SpecMeas> setFile( const std::string &displayFileName,
                                     const std::string &fileSystemLocation,
                                     SpecUtils::ParserType parseType = SpecUtils::ParserType::Auto );

  //setMeasurmentInfo(...) takes care of setting all the information,
  //  and serializing the the SpecMeas to a temporary location on file
  void setMeasurmentInfo( std::shared_ptr<SpecMeas> measurment );

  //saveToFileSystem(...): If there is currently a file pointed to by
  //  m_fileSystemLocation, than it will be deleted.
  //  If 'measurment' and measurementIfInMemory() are non-null, they must
  //  be equal (passing in a pointer can help keep object in
  //  m_weakMeasurmentPtr)
  //  If errors are encountered, setKeepCachedInMemmorry(true) is called, but
  //  no exceptions thrown
  //  The saving is actually done in a thread owned by WServer.
  //XXX - Should consider adding an on error callback!
  void saveToFileSystem( std::shared_ptr<SpecMeas> measurment )  const;

  //saveToFileSystemImmediately(...): saves (in calling thread) the passed in
  //  SpecMeas object; if a copy of the SpecMeas object in memmory can be
  //  obtained vai measurementIfInMemory(), then the passed in SpecMeas
  //  object _must_ be the same, or will throw an exception.
  //Note that this function is needed (as apposed to saveToFileSystem(...))
  //  because it appears by the time the SpecMeas destructor is called, the
  //  std::shared_ptr's have all already lost their references to the
  //  SpecMeas object, meaning m_weakMeasurmentPtr is already been reset to
  //  not to point to anything.
  //  -we could still save in a background thread if we made a copy of the
  //   SpecMeas object being deleted, and then saved it in another thread,
  //   its just not totally obvious if that would be worth the while (although
  //   making the copy is fairly cheap - I *think*)
  void saveToFileSystemImmediately( SpecMeas *meas ) const;

  //errorSavingCallback(...): gets called when there is an error saving
  //  to native file.  Attempts to instead cache the measurment by
  //  setting m_cachedMeasurement equal to meas (so make sure you pass in
  //  correct meas)
  void errorSavingCallback( std::string fileLocation,
                            std::shared_ptr<SpecMeas> meas ) const;

  std::shared_ptr<SpecMeas> parseFile() const;
  
  int numSamples() const;
  Wt::WString displayName() const;
  const Wt::WDateTime &uploadTime() const;
  float totalLiveTime() const;
  float totalRealTime() const;
  float totalGammaCounts() const;
  float totalNeutronCounts() const;
  int numDetectors() const;
  bool hasNeutronDetector() const;
  const Wt::WDateTime &spectrumTime() const;
  bool passthrough() const;

  void releaseCacheReference();
  void setKeepCachedInMemmorry( bool cache = true );

  std::shared_ptr<SpecMeas> measurementIfInMemory() const;

#if( USE_DB_TO_STORE_SPECTRA )
  //Experimental DB section
  //A lot more work hase to be done to this code, including:
  //  --Making sure all ways a SpecMeas can be modified are accounted for
  //  [DONE]--Create method so that either the SpecMeas get saved only when the users
  //    session is terminating, or a mechanism so that if a file is too large
  //    for the database, then it will be stored as a temporary file on the
  //    filesystem
  //  [DONE]--Add in an option for the user to allow saving spectra to the database
  //  [DONE]--Add in mechanism to check if the user has previously loaded a spectrum,
  //    and if they had made any changes to it; if they did, notify them of this
  //  [DONE]--Add in compile time option to choose between a MySQL or a SQLite3
  //    database, and of course deal with passwords and such
  //  [DONE]--Investigate if using compression would be helpful, see
  //    http://www.boost.org/doc/libs/1_44_0/libs/iostreams/doc/classes/gzip.html#examples
  //  --Do a performance check to make sure database interactions arent
  //    happening in the main GUI thread, and that DB operations arent too
  //    expensive.
  //  --Verify that when a UserFileInDb gets deleted, so does its
  //    UserFileInDbData (and same with InterSpecUser and UserFileInDb)
  //  [DONE]--Use temporary files for caching the measurments, and the database
  //    only for long term storage
  
  //saveToDatabase(...): saves the SpecMeas to the database if it has been
  //  modified since last save, or never saved at all.
  //  Will throw runtime_error if spectrum is too large to save in database, or
  //  there are other issues with the database.
  //  A Wt::Dbo::StaleObjectException may be thrown if you are concurrently
  //  saving the same file from another session (unlikely)
  //  Does not consult m_candidateForSavingToDb.
  void saveToDatabase( std::shared_ptr<const SpecMeas> meas ) const;
  
  static void saveToDatabase( std::shared_ptr<SpecMeas> meas,
                              std::shared_ptr<SpectraFileHeader> header );
  
  //saveToDatabaseWorker(...): an exception safe version of the above.  If
  //  the operation fails a message is printed out to the user.
  static void saveToDatabaseWorker( std::shared_ptr<SpecMeas> measurment,
                                    std::shared_ptr<SpectraFileHeader> header );
  
  //saveToDatabaseFromTempFile(): saves the SpecMeas to the database,
  //  essentially just coping the file to the UserFileInDbData object.
  //  This funtions is useful when the spectrum is no longer in memory, but you
  //  want to save it to the database.
  //  A Wt::Dbo::StaleObjectException may be thrown if you are concurrently
  //  saving the same file from another session (unlikely)
  //  Does not consult m_candidateForSavingToDb.
  void saveToDatabaseFromTempFile() const;
  
  //saveToDatabaseFromTempFileWorker(): an exception safe version of the above.
  //  If the operation fails a message is printed out to the user.
  void saveToDatabaseFromTempFileWorker() const;

  //readFromDataBase(): reads serialized data from the database and
  //  de-serializes it.  Returned pointer should always be valid.
  //  Throws exception on error.
  std::shared_ptr<SpecMeas> readFromDataBase() const;
  
  //setBasicFileInDbInfo(...): sets the simple (integers, floats, timestamps)
  //  to UserFileInDb; does not serialize the SpecMeas data or deal with
  //  the UserFileInDbData object
  void setBasicFileInDbInfo( UserFileInDb *info ) const;
  
  //setDbEntry(...): deletes database entry represented by m_fileDbEntry (if it
  //  exists and doesnt equal entry) and then associate the spectrum represented
  //  by this object to 'entry' via setting m_fileDbEntry ot 'entry'.  Does not
  //  check that 'entry' is a correct UserFileInDb object, and does not
  //  serialize the SpecMeas object, nor mark it as modified.  May throw
  //  exception if 'entry' is invalid, or there are errors righting to database.
  void setDbEntry( Wt::Dbo::ptr<UserFileInDb> entry );
  
  //dbEntry(): returns m_fileDbEntry, which may be empty
  Wt::Dbo::ptr<UserFileInDb> dbEntry();
  
  //setNotACandiateForSavingToDb(...): call this so that the SpecMeas object
  //  represetnted by *this will not be saved to the database when the
  //  SpectraFileHeader destructor is called.
  void setNotACandiateForSavingToDb();
  
  //candidateForSavingToDb(): returns m_candidateForSavingToDb and-ed with if
  //  the user preference says if it should be saved
  bool candidateForSavingToDb() const;

  //shouldSaveToDb(): returns m_candidateForSavingToDb and-ed with if the user
  //  preference says if it should be saved
  bool shouldSaveToDb() const;
  
  std::shared_ptr<SpecMeas> resetFromDatabase(
                                          Wt::Dbo::ptr<UserFileInDb> dbfile );
#endif
  
  
protected:
  //initFile(): will return empty pointers in case of failure. Does not throw.
  std::shared_ptr<SpecMeas> initFile( const std::string &filename,
                                        const SpecUtils::ParserType parseTypeconst,
                                        std::string orig_file_ending ) throw();

public:
//  protected:
  //The below must remain public for the sake of
  //  SpectraFileModel::toRawIndex()/fromRawIndex(),
  //  or SpectraFileModel be made a friend class
  mutable std::string m_fileSystemLocation;
  mutable Wt::Dbo::ptr<UserFileInDb> m_fileDbEntry;
  
  InterSpec *m_viewer;
  
  std::string m_displayName;
  std::string m_uuid;
  Wt::WDateTime m_uploadTime;
  int m_numSamples;
  bool m_isPassthrough;

  float m_totalLiveTime;
  float m_totalRealTime;
  float m_totalGammaCounts;
  float m_totalNeutronCounts;
  int m_numDetectors;
  bool m_hasNeutronDetector;
  Wt::WDateTime m_spectrumTime;
  std::vector<SpectraHeader>     m_samples;
  
  //m_modifiedSinceDecode: only updated when writing to file, or database, and
  //  exists to catch the edge case where the SpecMeas object has been written
  //  to file and not in memorry, but it was requested to save the SpecMeas
  //  object to the database (to fill out UserFileInDb::userHasModified).
  mutable bool m_modifiedSinceDecode;
  
  //If caching is enable the SpecUtils::SpecFile object will be kept in memorry
  //  even if no where else references the object
  bool m_keepCache;
  mutable std::shared_ptr<SpecMeas> m_cachedMeasurement;

  //m_candidateForSavingToDb: we may not always want to save the spectrum to
  //  the database, even if the user preference is to save it.  This variable
  //  indicates this.  An example of this would be if we opened a file spectrum
  //  file in a previous session, and then re-opened it in this one, and the
  //  user decides to use there previous version of the spectra cause maybe
  //  they had done some peak id or something.
  //  Also, just because this variable is true doesnt mean we will save the file
  //  to the database (user prefernce may be to not)
  bool m_candidateForSavingToDb;
  
  //This weak pointer tracks if the SpecUtils::SpecFile object cooresping to
  //  this SpectraFileHeader exists anywhere in memmory, if so parseFile()
  //  will return this, rather than re-parsing the spectrum
  mutable std::weak_ptr<SpecMeas> m_weakMeasurmentPtr;

  mutable std::recursive_mutex m_mutex; //XXX - right now only used in a couple select places
  typedef std::lock_guard<std::recursive_mutex>  RecursiveLock;

  Wt::Dbo::ptr<InterSpecUser> m_user;
  std::shared_ptr<DataBaseUtils::DbSession> m_sql;  //same as m_viewer uses
  
  
  std::string m_appId;
  mutable Wt::WApplication *m_app;

#ifdef WT_USE_BOOST_SIGNALS2
  mutable boost::signals2::connection m_aboutToBeDeletedConnection;
#else
  mutable boost::signals::connection m_aboutToBeDeletedConnection;
#endif
};//class SpectraFileHeader


//SpectraFileModel should be renamed SpecFileModel
class SpectraFileModel : public Wt::WAbstractItemModel
{
  //This class creates a 2 level model where the first level (labeled
  //  FileHeaderLevel) cooresponds to spectrum files (stored as
  //  SpectraFileHeader's), the sublevel to this (labeled SampleLevel)
  //  cooresponds to one sample Measurement (stored as SpectraHeader's), which
  //  may coorespond to more than one detector.
  //  If there is only one sample in a file, the SampleLevel will not exist.
  //
  //WModelIndex's for FileHeaderLevel have an invalid parent, and row
  //  cooresponding to its place in m_spectra.
  //WModelIndex's for SampleLevel have a parent cooresponding to its respective
  //  FileHeaderLevel with the internal pointer pointing to that parent
  //
  //Could consider adding another level to the model to coorespond to a spectrum
  //  (eg when a sample had multiple detectors).
  //

public:
  enum DisplayFields
  {
    kDisplayName,
    kNumMeasurements,
    kLiveTime,
    kRealTime,
    kGammaCounts,
    kNeutronCounts,
    kSpectrumTime,
    kNumDetectors,
    kUploadTime,
    NumDisplayFields
  };//enum DisplayFields

  enum Level
  {
    FileHeaderLevel,
//    MeasurementLevel,
    SampleLevel,
    InvalidLevel
  };

  SpectraFileModel( Wt::WObject *parent = NULL );
  ~SpectraFileModel();
  
  virtual int columnCount( const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
  virtual int rowCount( const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
  virtual boost::any data( const Wt::WModelIndex &index,
                           int role = Wt::DisplayRole ) const;
  std::shared_ptr<SpectraFileHeader> fileHeader( int row );
  Wt::WModelIndex index( std::shared_ptr<SpectraFileHeader> header ) const;

  //fileHeader: a convience function to avoid calling index(...) first.
  //  Results may not be valid pointer
  std::shared_ptr<SpectraFileHeader> fileHeader( std::shared_ptr<const SpecUtils::SpecFile> measinfo );
  
#if( USE_DB_TO_STORE_SPECTRA )
  //dbEntry: a convience function for calling SpectraFileHeader::dbEntry(...)
  Wt::Dbo::ptr<UserFileInDb> dbEntry( std::shared_ptr<const SpecUtils::SpecFile> measinfo );
#endif
  
  
  //Following fnct untested wcjohns 20120807
  Wt::WModelIndex index( std::shared_ptr<const SpecUtils::SpecFile> measinfo ) const;

  std::shared_ptr<const SpectraFileHeader> fileHeader( int row ) const;
  int addRow( std::shared_ptr<SpectraFileHeader> fileInfo );
  virtual bool removeRows( int row, int count,
                           const Wt::WModelIndex &parent = Wt::WModelIndex() );

  Level level( const Wt::WModelIndex &index ) const;

  virtual Wt::WModelIndex parent( const Wt::WModelIndex &index ) const;
  virtual Wt::WModelIndex index( int row, int column,
                             const Wt::WModelIndex &parent = Wt::WModelIndex()
                                                             ) const;
  virtual boost::any headerData( int section,
                                 Wt::Orientation orientation = Wt::Horizontal,
                                 int role = Wt::DisplayRole ) const;
  virtual void sort( int column, Wt::SortOrder order = Wt::AscendingOrder );
  virtual void *toRawIndex( const Wt::WModelIndex &index ) const;
  virtual Wt::WModelIndex fromRawIndex( void *rawIndex ) const;

  virtual Wt::WFlags<Wt::ItemFlag> flags( const Wt::WModelIndex &index ) const;  //should return ItemIsSelectable for SpectraFiles with one Measurement, or for SpectraHeaders

protected:
  std::vector< std::shared_ptr<SpectraFileHeader> > m_spectra;
};//class SpectraFileModel


//DownloadSpectrumResource should be renamed SpectrumDownloadResource
class DownloadSpectrumResource : public Wt::WResource
{
  //A simple class to stream a spectrum upon demand for which ever items are
  //  highlighted in the tree view.
public:
  DownloadSpectrumResource( SpecUtils::SaveSpectrumAsType type,
                            SpecMeasManager *display,
                            Wt::WObject *parent = NULL );
  virtual ~DownloadSpectrumResource();

private:
  SpecUtils::SaveSpectrumAsType m_type;
  SpecMeasManager *m_display;
  Wt::WApplication *m_app;

  virtual void handleRequest( const Wt::Http::Request& request,
                               Wt::Http::Response& response);

public:
  
  static void write_file( std::ostream &output,
                          const SpecUtils::SaveSpectrumAsType type,
                          const std::shared_ptr<const SpecMeas> Measurement,
                          const std::set<int> &samplenums,
                          const std::vector<std::string> &detectornums,
                          const InterSpec *viewer );
  
  //handle_resource_request(): does the actual streaming of the SpecMeas.
  //  'samplenums' and 'detectornums' are only used for file formats where SpecMeas must be
  //  collapsed down into a single spectrum (i.e., Chn, IntegerSpcType, SpcBinaryFloat, SpcAscii,
  //  SpeIaea, Cnf, Tka, will have all the spectra summed into a single spectrum); for other formats
  //  the entire SpecMeas is written out.
  //
  //  Specifying empty detector and sample numbers will default to all samples/detectors.
  static void handle_resource_request(
                                SpecUtils::SaveSpectrumAsType type,
                                std::shared_ptr<const SpecMeas> Measurement,
                                const std::set<int> &samplenums,
                                const std::vector<std::string> &detectornums,
                                const InterSpec *viewer,
                                const Wt::Http::Request& request,
                                Wt::Http::Response& response );
};// class DownloadSpectrumResource


//SpecificSpectrumResource should be renamed SpectrumCustomDownloadResource
class SpecificSpectrumResource : public Wt::WResource
{
  //A simple class to stream a spectrum upon demand, intendend for use with the
  //  "Save As" dialog.
public:
  SpecificSpectrumResource( SpecUtils::SaveSpectrumAsType type,
                            Wt::WObject *parent = NULL );
  virtual ~SpecificSpectrumResource();

  void setSpectrum( std::shared_ptr<const SpecMeas> spec,
                    const std::set<int> &samplenums,
                    const std::vector<std::string> &detnames );

  Wt::Signal<> &downloadStarting();
  Wt::Signal<> &downloadComplete();
  
#if( ANDROID )
  const std::set<int> &samplenums();
  const std::vector<std::string> &detnames();
  SpecUtils::SaveSpectrumAsType type();
  std::shared_ptr<const SpecMeas> spectrum();
#endif
  
private:
  std::set<int> m_samplenums;
  std::vector<std::string> m_detnames;
  SpecUtils::SaveSpectrumAsType m_type;
  std::shared_ptr<const SpecMeas> m_spectrum;
  Wt::WApplication *m_app;

  virtual void handleRequest( const Wt::Http::Request& request,
                               Wt::Http::Response& response);

  Wt::Signal<> m_downloadComplete;
}; // class SpecificSpectrumResource


//DownloadCurrentSpectrumResource: to download
class DownloadCurrentSpectrumResource : public Wt::WResource
{
public:
  DownloadCurrentSpectrumResource( SpecUtils::SpectrumType spectrum,
                                   SpecUtils::SaveSpectrumAsType format,
                                   InterSpec *viewer,
                                   Wt::WObject *parent = NULL );
  virtual ~DownloadCurrentSpectrumResource();
  
  Wt::Signal<> &downloadStarting();
  Wt::Signal<> &downloadComplete();
private:
  SpecUtils::SpectrumType m_spectrum;
  SpecUtils::SaveSpectrumAsType m_format;
  InterSpec *m_viewer;
  Wt::WApplication *m_app;
  
  virtual void handleRequest( const Wt::Http::Request& request,
                             Wt::Http::Response& response);
  
  Wt::Signal<> m_downloadComplete;
};// class SpecificSpectrumResource


#endif


