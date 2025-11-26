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
#include <boost/asio/deadline_timer.hpp>

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

// Some forward declarations
namespace SpecUtils{ enum class ParserType : int; }
namespace SpecUtils{ enum class SpectrumType : int; }

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
class BatchGuiDialog;
class SpecMeasManager;
class SpectraFileModel;
class PopupDivMenuItem;
class SpectraFileHeader;
class RowStretchTreeView;
class FileDragUploadResource;
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

  static void fileTooLarge( const ::int64_t size_tried );

  //The dataUploaded(..) with two arguments is to help keep from having to
  //  parse file twice to display when SpectraHeader is not caching the
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
  void newFileFromSelection();
  void sumSelectedSpectra();
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
                              const std::string &fileLocation,
                              SpecUtils::SpectrumType type );
  
  /** Handles parsing multiple DRF CSV/TSV files when dropped onto the app.
   
   Will parse file and prompt user to select a DRF, and if they also want to save the file for later use.
   
   @param input The input CSV/TSV file containing one or more DRFs.  See #DetectorPeakResponse::parseMultipleRelEffDrfCsv
          for a description of the file contents.
   @param displayName The display name (on the users system) of the DRF file.
   @param fileLocation The spool file location of the \p input stream; used for copying the file into the apps writable data
          directory if the user chooses to save it for later.
   
   @returns True if file contained DRFs and option of selecting one will be presented to the user.
   */
  bool handleMultipleDrfCsv( std::istream &input,
                             const std::string &displayName,
                             const std::string &fileLocation );
 
  /** Reads in GammaQuant CSV of detector efficiencies.
   
   The first cell must be "Detector ID", then rows through "Coefficient h".
   Each column is a different detector, with the first column being the labels ("Detector ID", "Calibration Geometry",
   "Comments", ..., "Coefficient h").
   */
  bool handleGammaQuantDrfCsv( std::istream &input,
                             const std::string &displayName,
                             const std::string &fileLocation );
  
  /** Reads a CALp file from input stream and then either applies the CALp to current spectra, or prompts the user how to apply it.
   Function may return asynchronously to the CALp being applied, as the application may prompt user for options/input.
   
   @param input The input CALp stream; will usually be an std::ifstream of a CALp file
   @param dialog The dialog to use for interacting with the user.  Must be valid.  Dialog contents will be cleared and replaced, only
          if this function returns true.  Dialog will be accepted (e.g. `WDialog::done(Accepted)` called) if \p autoApply is true, and
          the application is unambiguous.
   @param autoApply For unambiguous cases (e.g., CALp matches data detector) if the CALp should be applied, without bothering
          to prompt the user.  If auto applied, the dialog will be accepted.
   
   @returns True if dialog was modified, or CALp is accepted, or the CALp _may_ be applied, or an error message about the file (such
          as it is an invalid CALp file) will be given to the user.  Returns false if the CALp is not even considered for application and no
          error message will be presented to the user, and input \p dialog was not modified.
          (not a very clean mechanism, but it works for the current two use cases of this function)
   */
  bool handleCALpFile( std::istream &input, SimpleDialog *dialog, bool autoApply );
  
#if( USE_REL_ACT_TOOL )
  bool handleRelActAutoXmlFile( std::istream &input, SimpleDialog *dialog );
#endif
  
  /** Handles the Shielding/Source fit XML file */
  bool handleShieldingSourceFile( std::istream &input, SimpleDialog *dialog );
  
  /** Handles Source.lib files dropped onto the app. */
  bool handleSourceLibFile( std::istream &input, SimpleDialog *dialog );
  
  /** Handles the user dropping a .ECC file produced from ISOCS. */
  bool handleEccFile( std::istream &input, SimpleDialog *dialog );
  
  /** Some input files contain duplicate data - we will ask the user how they want to handle
   this, first handling "Derived Data", then "Multiple Energy Calibration Types", then
   "Multiple Virtual Detectors"
   */
  enum class VariantChecksToDo
  {
    /** No checks for things like multiple energy cals, derived data, or virtual detectors. */
    None = 0,
    
    /** RSI systems may have multiple Virtual Detectors (ex "VD1", "VD2", etc) defined, that may be sums of
     multiple detection elements, or not, or be overlapping, or whatever - so lets alert the user, and let them
     choose what to do.
     */
    MultiVirtualDets = 1,
    
    /** Some spectrum files will include both linear and compressed spectra, or multiple energy ranges of
     the same data.  The user may want to only load one of these, so they arent seeing duplicate data.
     */
    MultiEnergyCalsAndMultiVirtualDets = 2,
    
    /** For some files, the user may just want to see the derived data.
     Or perhaps just the raw data, without the essentially duplicate derived data.
     */
    DerivedDataAndMultiEnergyAndMultipleVirtualDets = 3
  };//enum class VariantChecksToDo
  
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
  
  
  void startSaveSelected();
  
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
  //  (e.g. WServer::instance()->ioService().boost::asio::io_service::post(...) ),
  //  where the WApplication could get terminated while the function is executing.
  //  The bool indicates if SpecMeasManager has destructed, and if not the mutex
  //  will keep it from destructing.
  void checkIfPreviouslyOpened( const std::string sessionID,
                                std::shared_ptr<SpectraFileHeader> header,
                                SpecUtils::SpectrumType type,
                                std::shared_ptr< std::mutex > mutex,
                                std::shared_ptr<bool> destructed );
  
  
  void showPreviousSpecFileUsesDialog( std::shared_ptr<SpectraFileHeader> header,
                                   const SpecUtils::SpectrumType type,
                                   const std::vector<Wt::Dbo::ptr<UserFileInDb>> &modifiedFiles,
                                   const std::vector<Wt::Dbo::ptr<UserFileInDb>> &unModifiedFiles,
                                   const std::vector<Wt::Dbo::ptr<UserState>> &userStatesWithFile );
  
  
  
  void userCanceledResumeFromPreviousOpened( std::shared_ptr<SpectraFileHeader> header );
  
  //saveToDatabase(...): saves the SpecMeas to the database in another thread;
  //  must be called from a thread where WApplication::instance() is available
  //  to ensure thread safety.
  void saveToDatabase( std::shared_ptr<const SpecMeas> meas ) const;

  
  void browsePrevSpectraAndStatesDb();
#endif
  
  //Some resources to enable drag-n-drop of spectrum files to various widgets.
  //  Call FileDragUploadResource::addDragNDropToWidget( WWebWidget * ) to
  //  enable a widgets drag-n-drop support.
  FileDragUploadResource *dragNDrop( SpecUtils::SpectrumType type );
  FileDragUploadResource *foregroundDragNDrop();
  FileDragUploadResource *secondForegroundDragNDrop();
  FileDragUploadResource *backgroundDragNDrop();
#if( USE_BATCH_GUI_TOOLS )
  FileDragUploadResource *batchDragNDrop();
#endif

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
  
  /** Called by FileDragUploadResource, as files are being uploaded; if file is larger than `sm_minNumBytesShowUploadProgressDialog`,
   then will create a status message dialog, `m_processingUploadDialog`.
   */
  void handleDataRecievedStatus( uint64_t num_bytes_recieved, uint64_t num_bytes_total, SpecUtils::SpectrumType type );
  
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

#if( USE_BATCH_GUI_TOOLS )
  void showBatchDialog();
  void handleBatchDialogFinished();
#endif

#if( USE_QR_CODES )
  void handleSpectrumUrl( std::string &&url );
  void displaySpectrumQrCode( const SpecUtils::SpectrumType type );
  void multiSpectrumDialogDone();
#endif
  
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
                                              VariantChecksToDo viewingChecks );

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
  
  void selectVirtualDetectorChoice( const std::set<std::string> tokeep,
                           std::shared_ptr<SpectraFileHeader> header,
                           std::shared_ptr<SpecMeas> meas,
                           const SpecUtils::SpectrumType type,
                           const bool checkIfPreviouslyOpened,
                           const bool doPreviousEnergyRangeCheck );
  
  void handleCancelPreviousStatesDialog( AuxWindow *dialog );
  void handleClosePreviousStatesDialogAfterSelect( AuxWindow *dialog );
  
  /** Deletes the dialog, but only if the passed in dialog is the same as `m_processingUploadDialog` */
  void checkCloseUploadDialog( SimpleDialog *dialog, Wt::WApplication *app );
  
private:
  Wt::WContainerWidget *createButtonBar();
  void deleteSpectrumManager();
  Wt::WContainerWidget *createTreeViewDiv();
  void closeNonSpecFileDialog();
  
protected:
  AuxWindow *m_spectrumManagerWindow;
  RowStretchTreeView    *m_treeView;
  SpectraFileModel *m_fileModel;
  Wt::WFileUpload  *m_fileUpload;
  InterSpec   *m_viewer;
  Wt::WContainerWidget *m_spectrumManagertreeDiv;
  
  Wt::WPushButton *m_setButton;
  Wt::WPopupMenuItem  *m_setAsForeground;
  Wt::WPopupMenuItem  *m_setAsBackground;
  Wt::WPopupMenuItem  *m_setAsSecForeground;
  Wt::WPushButton *m_combineToNewFileButton;
  Wt::WPushButton *m_subsetOfMeasToNewFileButton;
  Wt::WPushButton *m_sumSpectraButton;
  Wt::WPushButton *m_saveButton;
  Wt::WPushButton *m_deleteButton;
  Wt::WPushButton *m_removeForeButton;
  Wt::WPushButton *m_removeBackButton;
  Wt::WPushButton *m_removeFore2Button;
  
  Wt::WContainerWidget *m_infoHandler;
  std::vector< std::string > m_labelTags;
  
  FileDragUploadResource *m_foregroundDragNDrop;
  FileDragUploadResource *m_secondForegroundDragNDrop;
  FileDragUploadResource *m_backgroundDragNDrop;
#if( USE_BATCH_GUI_TOOLS )
  FileDragUploadResource *m_batchDragNDrop;
#endif
  
  SimpleDialog *m_multiUrlSpectrumDialog;
  
  //m_sql same as m_viewer->sql();
  std::shared_ptr<DataBaseUtils::DbSession> m_sql;
  
  std::shared_ptr<std::mutex> m_destructMutex;
  std::shared_ptr<bool> m_destructed;
  
  /** Dialog shown when user loads a file that has been previously either explicitly saved state, or automatically saved. */
  AuxWindow *m_previousStatesDialog;
  
  /** Dialog created when a file is uploading, from `handleDataRecievedStatus()`, letting the user know that something is going on.
   On the client-side there is an upload status that is shown, but for large files the upload status in JS can be significantly off from what
   the server is seeing, leading to a potentially long gap between the client-side status, and when the file actually finishes uploading and
   is being parsed (which another dialog will be shown for parsing large files) - so this dialog covers this time gap to keep the user informed.
   */
  SimpleDialog *m_processingUploadDialog;
  
  /** Timer used to make sure the dialog is destroyed, even if upload stalls-out; reset in every call to `handleDataRecievedStatus()`,
   with a timeout of 30 seconds.
   
   Did not use a `Wt::WTimer` so we dont effect the DOM when creating and resetting timer, as the WTimer gets inserted into the DOM, and
   maybe also added to the client-side JS (but I didnt check if creating/resetting a WTimer would cause a network request, but a quick glance at
   the Wt source code did look like this might be the case, but also maybe not if its server-side only connection)
   */
  std::unique_ptr<boost::asio::deadline_timer> m_processingUploadTimer;

  /** Dialog created when a non-spectrum file is dropped on the app. */
  SimpleDialog *m_nonSpecFileDialog;

#if( USE_BATCH_GUI_TOOLS )
  /** Dialog created when a batch of files is dropped on the app. */
  BatchGuiDialog *m_batchDialog;
#endif
  
#if( !defined(MAX_SPECTRUM_MEMMORY_SIZE_MB) ||  MAX_SPECTRUM_MEMMORY_SIZE_MB < 0 )
  static const size_t sm_maxTempCacheSize = 0;
#else
  static const size_t sm_maxTempCacheSize = 1024 * 1024 * MAX_SPECTRUM_MEMMORY_SIZE_MB;
#endif

  const static size_t sm_minNumBytesShowUploadProgressDialog = 100 * 1024;
  
  mutable std::deque< std::shared_ptr<const SpecMeas> > m_tempSpectrumInfoCache;
};//class SpecMeasManager

#endif
