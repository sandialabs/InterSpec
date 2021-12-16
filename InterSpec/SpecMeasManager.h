#ifndef SpecMeasManager_h
#define SpecMeasManager_h
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
#include <mutex>
#include <memory>
#include <vector>
#include <string>

#include <boost/any.hpp>

#include <Wt/WString>
#include <Wt/WResource>
#include <Wt/WDateTime>
#include <Wt/WModelIndex>
#include <Wt/Dbo/Session>
#include <Wt/WContainerWidget>
#include <Wt/Dbo/SqlConnection>
#include <Wt/WAbstractItemModel>

#include "InterSpec/InterSpec.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/SpectraFileModel.h"

namespace SpecUtils{ enum class ParserType : int; }
namespace SpecUtils{ enum class SpectrumType : int; }

// Some forward declarations
namespace Wt
{
  class JSlot;
  class WObject;
  class WCheckBox;
  class WTreeView;
  class WFileUpload;
  class WButtonGroup;
  class WApplication;
  class WSelectionBox;
  class WContainerWidget;
}//namespace Wt

class PopupDiv;
class SpecMeas;
class UserFileInDb;
class PopupDivMenu;
class SpectraHeader;
class InterSpec;
class SimpleDialog;
class SpecMeasManager;
class SpectraFileModel;
class PopupDivMenuItem;
class SpectraFileHeader;
class RowStretchTreeView;
#if( !ANDROID && !IOS )
class FileDragUploadResource;
#endif
class DownloadSpectrumResource;
class SpecificSpectrumResource;

namespace SpecUtils{ enum class SpectrumType : int; }

namespace DataBaseUtils
{
  class DbSession;
}


class SpecMeasManager : public Wt::WObject
{
public:
  SpecMeasManager( InterSpec *viewer );
  void  startSpectrumManager();
  virtual ~SpecMeasManager();
  SpectraFileModel *model();
  const SpectraFileModel *model() const;
  RowStretchTreeView *treeView();
  const RowStretchTreeView *treeView() const;
  const InterSpec *viewer() const;

  //setFile(): returns the row in m_fileModel of the new file
  //  thows std::runtime_error(..) on failure
  int setFile( const std::string &displayName,
               const std::string &fileLocation,
                std::shared_ptr<SpectraFileHeader> &header,
                std::shared_ptr<SpecMeas> &Measurement,
                SpecUtils::ParserType parseType = SpecUtils::ParserType::Auto );

  static void displayInvalidFileMsg( std::string filename, std::string errormsg );
  
  //Adds measurment to the file model, and user database, returning a header for
  //  it.  The measurment is not checked to see if its appropriate for
  //  displaying or working with; that is, conditions such as multiple
  //  calibrations are not checked for.
  std::shared_ptr<SpectraFileHeader> addFile( const std::string &displayName,
                                      std::shared_ptr<SpecMeas> measurement );
  
  //cleanupQuickSaveAsDialog deletes the dialog and sets the m_specificResources
  //  data to null
  void cleanupQuickSaveAsDialog( AuxWindow *dialog, Wt::WApplication *app );

  static void fileTooLarge( const ::int64_t size_tried );

  //The dataUploaded(..) with two arguments is to help keep from having to
  //  parse file twice to display when SpectraHeader is not cahcing the
  //  SpecMeas obj.
  //XXX - The dataUploaded() with one argument is called whe a user uploads from
  //      the 'File Manager' screen - when they will probably want to load the
  //      file immediately - if SpectraHeader is not chaching, this class should
  //      cache the SpecMeas untill the 'File Manager' is,
  //      closed/minimized - w/ this modification I think it would be acceptable
  //      to not cache files in SpectraHeader (except for its weak_ptr<>).

  //Returns the row in m_fileModel of the new file, -1 on error
  int dataUploaded2( Wt::WFileUpload *upload , SpecUtils::SpectrumType type);
  int dataUploaded( Wt::WFileUpload *upload );
  int dataUploaded( Wt::WFileUpload *upload,
                    std::shared_ptr<SpecMeas> &meas_ptr );
  
  //loadFromFileSystem(...) loads a file from disk, as if it were uploaded.
  //  Returns whether or not file was loaded
  bool loadFromFileSystem( const std::string &filename, SpecUtils::SpectrumType type,
                           SpecUtils::ParserType parseType = SpecUtils::ParserType::Auto );
  
  
  void selectionChanged();
  void removeSelected();
  
  //removeSpecMeas(): removes the entry cooresponding to the passed in SpecMeas
  //  from the SpectraFileModel.   If 'undisplay' is specified, and the SpecMeas
  //  was being displayed, then the SpecMeas will also be un-displayed from the
  //  GUI.  Calling this function will also prevent the SpectraFileHeader
  //  cooresponding to 'meas' from saving the file to disk from the destructor
  //  of SpectraFileHeader
  void removeSpecMeas( std::shared_ptr<const SpecMeas> meas,
                       const bool undisplay );
  
#if( USE_DB_TO_STORE_SPECTRA )
  void removeSpecMeas( Wt::Dbo::ptr<UserFileInDb> row );
#endif
  void removeAllFiles();
  void renameSaveAsFile();
  void newFileFromSelection();
  std::shared_ptr<SpecMeas> selectedToSpecMeas() const;
  void unDisplay( SpecUtils::SpectrumType type );
  
  //the loadSelected() with second argument helps keep from having to parse
  //  the file twice when SpectraHeader is not cacheing the SpecMeas obj.
  void loadSelected( const SpecUtils::SpectrumType type, const bool doPreviousEnergyRangeCheck );
  void loadSelected( const SpecUtils::SpectrumType type, std::shared_ptr<SpecMeas> ptr,
                     const bool doPreviousEnergyRangeCheck  );

  void startQuickUpload();
  void finishQuickUpload( Wt::WFileUpload *upload, const SpecUtils::SpectrumType type );

  
  //handleNonSpectrumFile(): if a file couldnt be parsed, this function will
  //  look at the file, and try to suggest to the user what type of file it...
  //  If it returns true, then it could determine the file type, and notified
  //  the user.  If false, the file type is not known.
  //  Meant to be called from within the event loop.
  bool handleNonSpectrumFile( const std::string &displayName,
                              const std::string &fileLocation );
  
  bool handleCALpFile( std::istream &input, SimpleDialog *dialog, bool autoApply );
  
  enum class VariantChecksToDo
  {
    None = 0,
    MultipleEnergyCal = 1,
    DerivedDataAndEnergy = 2
  };
  
  // displayFile(...) displays the file passed in as specified type, if it can.
  //  --if kForground is specified and the measurment contains a background
  //    spectrum then the background will be loaded as well (as seperate graph).
  //  --if SpecUtils::SpectrumType::SecondForeground is specified, then only non-bakcground and
  //    non-calibration spectra will be displayed.
  //  --if SpecUtils::SpectrumType::Background is specified and the number of measurments is greater
  //    than 1, the first spectrum in the file will be loaded while warning user
  //For all cases, if the file contains a calibration measurment, it will not be
  //loaded and the user will be notified of this.
  //Function may throw runtime_error in case of fileModelRow doesnt coorespond
  //  to measement_ptr, or other inconsistencies or failed logic
  //
  //Note that the spectrum should not be assumed to have been loaded after this
  //  function call completes, but may be completed later after some
  //  (asyncronous) interactions with the user.
  //
  //If fileModelRow is negative AND measement_ptr is empty, then current
  //  spectrum will be unloaded
  //if checkIfPreviouslyOpened is specified, then the database will be searched
  //  to see if the user has previously worked with the same spectrum file, and
  //  if so either prompt them if they would like to resume with that work (if
  //  any had taken place) or associate the SpectraFileHeader with that DB
  //  entry if no work had taken place.  If the found matching SpectraFileHeader
  //  has been loaded in the current session, and has not been modified, then
  //  that SpectraFileHeader will be displayed rather than the one passed in.
  //If doPreviousEnergyRangeCheck is true then will be checked if this spectrum
  //  is a candidate to keep the energy calibration the same as the current
  //  foreground, and if it is, the user will be prompted about this.
  //If checkIfAppropriateForViewing is true, then SpecMeas will be checked if
  //  it has multiple binnings specified for the same data, and if so, prompt
  //  the user which one they would like to keep.
  //TODO: modify SpectraFileModel to return row of a given MeasurmentInfo
  //      pointer, which would remove need to pass in fileModelRow into this fcn
  void displayFile( int fileModelRow,
                    std::shared_ptr<SpecMeas> measement_ptr,
                    const SpecUtils::SpectrumType type,
                    bool checkIfPreviouslyOpened,
                    const bool doPreviousEnergyRangeCheck,
                    const VariantChecksToDo viewingChecks
                    );

  std::shared_ptr<SpectraFileHeader> selectedFile() const; //returns first file selected (I think)
  std::vector<std::shared_ptr<SpectraFileHeader> > getSelectedFiles() const; //returns all selected
  std::set<int> selectedSampleNumbers() const;
  void setDisplayedToSelected();
  
  
  // TODO There may be a race condition in the displayQuickSaveAsDialog()
  // ... I think its fine now, but could use some more double checking

  //  displayQuickSaveAsDialog() will eventually become depreciated I think
  //  after testing USE_SAVEAS_FROM_MENU section of code from InterSpec
  void displayQuickSaveAsDialog();


  void displayIsBeingShown();
  void displayIsBeingHidden();
  
  //clearTempSpectrumInfoCache(): clears the m_tempSpectrumInfoCache after
  //  trying to save them to disk
  void clearTempSpectrumInfoCache();
  
  //addToTempSpectrumInfoCache(...): places the SpecMeas into
  //  m_tempSpectrumInfoCache (if not already in) so it will remain in memorry.
  //  The function then goes through m_tempSpectrumInfoCache and deletes older
  //  SpecMeas objects until memorry taken up is less than sm_maxTempCacheSize.
  void addToTempSpectrumInfoCache( std::shared_ptr<const SpecMeas> meas ) const;
  
  //removeFromSpectrumInfoCache(...): removes the passed in spectrum from
  //  the cache, allowing it to be destructed and memorry be reclaimed.
  //  If saveToDisk is specified then it will be first saved to disk for
  //  later access (to avoid doing this from within the SpecMeas destructor).
  void removeFromSpectrumInfoCache( std::shared_ptr<const SpecMeas> meas,
                                    bool saveToDisk ) const;

  //serializeToTempFile(...): intended to be called right before a
  //  SpecMeas object is expected to be destructed in order to save it to a
  //  temporary file or the database, so that if the user loads it later it can
  //  be read from disk and be in its current state.  Note that this would
  //  happen anyway from the SpecMeas destructor, but this allows it to be more
  //  efficient and to happen in a seperate thread.
  void serializeToTempFile( std::shared_ptr<const SpecMeas> meas ) const;
  void uploadSpectrum();
  
#if( USE_DB_TO_STORE_SPECTRA )
  //setDbEntry(...): analagous to setFile(...)
  int setDbEntry( Wt::Dbo::ptr<UserFileInDb> dbfile,
                  std::shared_ptr<SpectraFileHeader> &header,
                  std::shared_ptr<SpecMeas> &measurement,
                  bool enforceUser );
  
  //checkIfPreviouslyOpened(...): checks if the user has previously opened this
  //  file, and if appropriate prompts them if they would like to use there
  //  last session with the spectrum.
  //The mutex and bool ptr are to avoid a race condition when this function is
  //  called from outside the event loop
  //  (e.g. WServer::instance()->ioService().post(...) ), where teh WApplication
  //  could get terminated while the function is executing.  The bool indicates
  //  if SpecMeasManager has destructed, and if not the mutex will keep it from
  //  destructing.
  void checkIfPreviouslyOpened( const std::string sessionID,
                                std::shared_ptr<SpectraFileHeader> header,
                                SpecUtils::SpectrumType type,
                                std::shared_ptr< std::mutex > mutex,
                                std::shared_ptr<bool> destructed );
  
  void createPreviousSpectraDialog( const std::string sessionID,
                                    std::shared_ptr<SpectraFileHeader> header,
                                    const SpecUtils::SpectrumType type,
                                    const std::vector< Wt::Dbo::ptr<UserFileInDb> > modifiedFiles,
                                   const std::vector< Wt::Dbo::ptr<UserFileInDb> > unModifiedFiles );
  
  void userCanceledResumeFromPreviousOpened( AuxWindow *window,
                                 std::shared_ptr<SpectraFileHeader> header );
  
  //saveToDatabase(...): saves the SpecMeas to the database in another thread;
  //  must be called from a thread where WApplication::instance() is available
  //  to ensure thread safety.
  void saveToDatabase( std::shared_ptr<const SpecMeas> meas ) const;
  
  void storeSpectraInDb();
  void finishStoreAsSpectrumInDb( Wt::WLineEdit *name,
                                  Wt::WTextArea *description,
                                  std::shared_ptr<SpecMeas> meas,
                                  AuxWindow *window );
  void storeSpectraSnapshotInDb( const std::string name = "" );
  void finishSaveSnapshotInDb(
                      const std::vector< std::shared_ptr<SpecMeas> > specs,
                      const std::vector< Wt::Dbo::ptr<UserFileInDb> > dbs,
                      const std::vector< Wt::WLineEdit * > edits,
                      const std::vector< Wt::WCheckBox * > cbs,
                               AuxWindow *window );
  void startStoreSpectraAsInDb();
  void browseDatabaseSpectrumFiles( SpecUtils::SpectrumType type );
  void showPreviousDatabaseSpectrumFiles( SpecUtils::SpectrumType type, std::shared_ptr<SpectraFileHeader> header );
#endif
  
#if( !ANDROID && !IOS )
  //Some resources to enable drag-n-drop of spectrum files to various widgets.
  //  Call FileDragUploadResource::addDragNDropToWidget( WWebWidget * ) to
  //  enable a widgets drag-n-drop support.
  FileDragUploadResource *dragNDrop( SpecUtils::SpectrumType type );
  FileDragUploadResource *foregroundDragNDrop();
  FileDragUploadResource *secondForegroundDragNDrop();
  FileDragUploadResource *backgroundDragNDrop();
#endif
  
#if( SUPPORT_ZIPPED_SPECTRUM_FILES )
  //handleZippedFile:  presents the user with a dialog to extract and use one
  //  of the spectrum files in a zip archive.  Returns true if a valid zip file.
  //  'name' is the display name of the original file, while 'spoolName' is
  //  is the actual file on disk.  If spectrum_type is not SpecUtils::SpectrumType::Foreground,
  //  SpecUtils::SpectrumType::SecondForeground, or SpecUtils::SpectrumType::Background, then user will be given option of what
  //  to open it as.
  // ToDo: having SpecUtils::SpectrumType potentially be invalid feels very wrong - should fix signature of function up
  bool handleZippedFile( const std::string &name,
                         const std::string &spoolName,
                         const SpecUtils::SpectrumType spectrum_type );
  
  //group should have three buttons coorespnding SpecUtils::SpectrumType::Foreground, SpecUtils::SpectrumType::Background, and
  //  SpecUtils::SpectrumType::SecondForeground.
  //if index is valid, open the file that row cooreponds to; if not, will open
  //  the selected row (if none selected open nothin, but that shouldnt ever
  //  happen)
  void extractAndOpenFromZip( const std::string &spoolName,
                              Wt::WButtonGroup *group,
                              Wt::WTreeView *table,
                              AuxWindow *window,
                              Wt::WModelIndex index );
  
#endif
  
  //Handles a file dropped onto the application, or finishes opening files from
  //  filesystem URL.
  //  Does not delete the file after opening.
  //  Note: This function may return immediately, posting doing the actual work to another thread.
  //        If you want to complete the parsing/opening of the file before returning, call
  //        #handleFileDropWorker.
  void handleFileDrop( const std::string &name,
                       const std::string &spoolName,
                       SpecUtils::SpectrumType type );
  
  void handleFileDropWorker( const std::string &name,
                       const std::string &spoolName,
                       SpecUtils::SpectrumType type,
                       SimpleDialog *dialog,
                       Wt::WApplication *app );

protected:
  //Called from inside displayFile(...) to see if there are options for
  //  displaying the file the user should select before continuing.
  //
  //Returns true if the user needs to be prompted, false if loading file can
  //  be continued.
  //
  //Throws exception if any input is invalid.
  bool checkForAndPromptUserForDisplayOptions( std::shared_ptr<SpectraFileHeader> header,
                                              std::shared_ptr<SpecMeas> measement_ptr,
                                              const SpecUtils::SpectrumType type,
                                              const bool checkIfPreviouslyOpened,
                                              const bool doPreviousEnergyRangeCheck,
                                              const VariantChecksToDo viewingChecks );

  void selectEnergyBinning( const std::string binning,
                            std::shared_ptr<SpectraFileHeader> header,
                            std::shared_ptr<SpecMeas> meas,
                            const SpecUtils::SpectrumType type,
                            const bool checkIfPreviouslyOpened,
                            const bool doPreviousEnergyRangeCheck );
  
  enum class DerivedDataToKeep{ All, RawOnly, DerivedOnly };
  void selectDerivedDataChoice( const DerivedDataToKeep tokeep,
                           std::shared_ptr<SpectraFileHeader> header,
                           std::shared_ptr<SpecMeas> meas,
                           const SpecUtils::SpectrumType type,
                           const bool checkIfPreviouslyOpened,
                           const bool doPreviousEnergyRangeCheck );
  
  
private:
  AuxWindow* m_spectrumManagerWindow;
  Wt::WContainerWidget *createButtonBar();
  void deleteSpectrumManager();
  Wt::WContainerWidget *createTreeViewDiv();

  
protected:
    
  RowStretchTreeView    *m_treeView;
  SpectraFileModel *m_fileModel;
  Wt::WFileUpload  *m_fileUpload;
  InterSpec   *m_viewer;
  Wt::WContainerWidget *m_spectrumManagertreeDiv;
    
//  Wt::WPushButton  *m_setPrimaryButton;
//  Wt::WPushButton  *m_setSecondaryButton;
//  Wt::WPushButton  *m_setBackgroundButton;
//  Wt::WPushButton  *m_unDisplayMenuButton;
  
  
  Wt::WPushButton *m_setButton;
  Wt::WPopupMenuItem  *m_setAsForeground;
  Wt::WPopupMenuItem  *m_setAsBackground;
  Wt::WPopupMenuItem  *m_setAsSecForeground;
  Wt::WPushButton *m_combineButton;
  Wt::WPushButton *m_saveButton; // replaces m_saveFileAsButton
  Wt::WPushButton *m_deleteButton;
  Wt::WPushButton *m_removeForeButton;
  Wt::WPushButton *m_removeBackButton;
  Wt::WPushButton *m_removeFore2Button;
  
  
  
//  PopupDivMenu     *m_unDisplayPopup;
//  PopupDivMenuItem *m_unDisplayPrimaryButton;
//  PopupDivMenuItem *m_unDisplaySecondaryButton;
//  PopupDivMenuItem *m_unDisplayBackgroundButton;
  
  Wt::WContainerWidget *m_infoHandler;
  std::vector< std::string > m_labelTags;

//  Wt::WPushButton  *m_makeNewFileButton;
//  Wt::WPushButton  *m_removeFileButton;
//  Wt::WPushButton  *m_saveFileAsButton;
//  PopupDivMenu     *m_saveAsPopup;
  
  DownloadSpectrumResource *m_downloadResources[static_cast<int>(SpecUtils::SaveSpectrumAsType::NumTypes)];
  SpecificSpectrumResource *m_specificResources[static_cast<int>(SpecUtils::SaveSpectrumAsType::NumTypes)];
  
#if( !ANDROID && !IOS )
  FileDragUploadResource *m_foregroundDragNDrop;
  FileDragUploadResource *m_secondForegroundDragNDrop;
  FileDragUploadResource *m_backgroundDragNDrop;
#endif
  
  //m_sql same as m_viewer->sql();
  std::shared_ptr<DataBaseUtils::DbSession> m_sql;
  
  std::shared_ptr< std::mutex > m_destructMutex;
  std::shared_ptr< bool > m_destructed;
  

#if( !defined(MAX_SPECTRUM_MEMMORY_SIZE_MB) ||  MAX_SPECTRUM_MEMMORY_SIZE_MB < 0 )
  static const size_t sm_maxTempCacheSize = 0;
#else
  static const size_t sm_maxTempCacheSize = 1024 * 1024 * MAX_SPECTRUM_MEMMORY_SIZE_MB;
#endif

  mutable std::deque< std::shared_ptr<const SpecMeas> > m_tempSpectrumInfoCache;
};//class SpecMeasManager

#endif
