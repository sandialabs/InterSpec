#ifndef DrfSelect_h
#define DrfSelect_h
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

#include <atomic>
#include <memory>
#include <string>
#include <vector>
#include <cstdint>
#include <fstream>
#include <functional>

#include <Wt/WColor.h>
#include <Wt/WSignal.h>
#include <Wt/WTextArea.h>
#include <Wt/WContainerWidget.h>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpecUser.h"  //for UseDrfPref::UseDrfType enum only
#include "InterSpec/DetectorPeakResponse.h"


namespace Wt
{
  class WText;
  class WTimer;
  class WCheckBox;
  class WLineEdit;
  class WPushButton;
  class WStackedWidget;
}//namespace Wt

class InterSpec;
class InterSpecUser;
class DrfLibraryList;
class DrfFormulaForm;
class DrfPreviewPanel;
class DrfLocationsPage;
class DrfFoundList;
struct DrfDropBatch;

namespace DrfImport
{
  struct FoundDetector;
}
class DrfModifyWindow;
class DrfImportWidget;
class DetectorPeakResponse;

namespace DrfLibrary
{
  struct Entry;
  struct Location;
  struct Snapshot;
  struct LocationStatus;
}

namespace DataBaseUtils
{
  class DbSession;
}

namespace SpecUtils
{
  class SpecFile;
  enum class DetectorType : int;
}

/** The detector in use - its icon and name - shown in the tools that need one; clicking it (or
 Enter / Space when it has focus) opens the Detector Efficiency dialog.  It follows
 `InterSpec::detectorChanged()` / `detectorModified()`.
 */
class DetectorDisplay : public Wt::WContainerWidget
{
public:
  /** `specViewer` must be valid, and outlive this widget. */
  explicit DetectorDisplay( InterSpec *specViewer );

  void setDetector( std::shared_ptr<DetectorPeakResponse> det );
  void editDetector();

  std::shared_ptr<DetectorPeakResponse> detector();
  std::shared_ptr<const DetectorPeakResponse> detector() const;

protected:
  /** The generic detector icon (images/detector_generic.svg, as a mask in the text colour), shown
   unless the DRF has its own (#m_drfIcon). */
  Wt::WContainerWidget *m_genericIcon;
  Wt::WText *m_drfIcon;
  Wt::WText *m_text;

  InterSpec *m_interspec;

  /** Weak: the detector is owned by the foreground spectrum (or by the tool that set it). */
  std::weak_ptr<DetectorPeakResponse> m_currentDetector;
};//class DetectorDisplay


/** The "Detector Efficiency" dialog: a Library of every detector InterSpec knows about (built in,
 found in standard and user folders, and saved in the database), a preview of the one picked
 compared against the one in use, and an Import tab for bringing in new ones.

 Picking a detector only previews it; "Use" makes it the detector in use, app-wide (with an undo
 step), and closing the dialog never applies a detector that was only previewed.

 The detector list comes from `DrfLibrary`, which caches what it reads, so re-opening the dialog only
 re-reads files and database rows that changed.
 */
class DrfSelect : public Wt::WContainerWidget
{
public:
  /** DrfSelect Constructor
   @param currentDet The detector in use.
   @param parentWindow Used to place the Cancel / Use buttons and help into its footer.
   */
  DrfSelect( std::shared_ptr<DetectorPeakResponse> currentDet,
             InterSpec *specViewer,
             AuxWindow *parentWindow = nullptr );

  virtual ~DrfSelect();

  /** Decodes the URL, and then prompts the user if they want to use the DRF. */
  static void handle_app_url_drf( const std::string &url_query );

  /** Emitted when the dialog should close. */
  Wt::Signal<> &done();

  /** The detector in use (the foreground's). */
  std::shared_ptr<DetectorPeakResponse> detector();

  /** The detector previewed (fully loaded), or null. */
  std::shared_ptr<const DetectorPeakResponse> previewedDetector() const;

  /** Makes `drf` the detector in use, app-wide (with an undo step), records it in the database, and
   - if those boxes are ticked - as the detector to pick automatically for the foreground's serial
   number or model.  The dialog stays open, showing it as "Just changed".
   */
  void useDetector( std::shared_ptr<DetectorPeakResponse> drf );

  /** Leaves no detector in use. */
  void useNoDetector();

  /** Slot for `InterSpec::detectorChanged()` / `detectorModified()`: the "In use now" strip and the
   list follow the detector in use, however it changed. */
  void setDetector( std::shared_ptr<DetectorPeakResponse> det );

  /** Searches the database for DRFs from the specified user with hashes
      matching the passed in DRF, and updates the last used time for all of
      them (should only find at most one).  If the DRF does not exist in the DB
      then it is added to is.

      The appearance (icon / 3D model) is not in the hash, so a matching row is given the DRF's
      appearance if it differs; one is only ever removed from a row when `appearance_is_authoritative`
      is true (the user just edited the DRF), since a DRF without one may simply be an older copy.
   */
  static void updateLastUsedTimeOrAddToDb( std::shared_ptr<DetectorPeakResponse> drf,
                                           long long db_user_id,
                                           std::shared_ptr<DataBaseUtils::DbSession> sql,
                                           const bool appearance_is_authoritative = false );
  
  /** Checks the database to see if the measurement serial number or measurement
      model cooresponds to a user preference in the database for  a DRF to use.
      If found, returns the DRF, if not, returns null.

      Takes the user's database id rather than a `Wt::Dbo::ptr`, so it can be called from a worker
      thread (a `Dbo::ptr` copy touches a non-atomic reference count).
   */
  static std::shared_ptr<DetectorPeakResponse> getUserPreferredDetector(
                                            std::shared_ptr<DataBaseUtils::DbSession> sql,
                                            long long db_user_id,
                                            const std::string &serial_number,
                                            SpecUtils::DetectorType detType,
                                            const std::string &detector_model );
  
  /** Records `drf` as the detector to pick automatically for spectra with `meas`'s serial number
   or model (per `prefType`), replacing any earlier choice.  Call on the session thread: it copies
   the `Wt::Dbo::ptr` user, whose reference count is not thread safe.
   */
  static void setUserPrefferedDetector( std::shared_ptr<DetectorPeakResponse> drf,
                                        std::shared_ptr<DataBaseUtils::DbSession> sql,
                                        Wt::Dbo::ptr<InterSpecUser> user,
                                        UseDrfPref::UseDrfType prefType,
                                        std::shared_ptr<const SpecUtils::SpecFile> meas );

  /** The directories #initAGadrasDetector searches, separated by semicolons or newlines: the user
   preference "GadrasDRFPath", or the fixed data-directory locations on web deployments.

   Reads a user preference, so call it on the session thread and pass the result to any worker.
   */
  static std::string gadrasDrfSearchPaths( InterSpec *interspec );

  //Will init detector in the data/detector_responses folder, and return result.
  //  throws exception if there is an error.
  //Type should be a DetectorType enum (just not including that header to save
  //   on deplandcies)
  //  `searchPaths` is from #gadrasDrfSearchPaths; taking it as a value lets this run off the
  //  session thread.
  static std::shared_ptr<DetectorPeakResponse> initAGadrasDetector( const SpecUtils::DetectorType type,
                                                                     const std::string &searchPaths );

  /** Looks in `searchPaths` (see #gadrasDrfSearchPaths) for a detector matching specified name.
   Throws exception in error; returned detector should always be valid.
   */
  static std::shared_ptr<DetectorPeakResponse> initAGadrasDetector( const std::string &name,
                                                                     const std::string &searchPaths );
  
  /** Inits adetector from a directory that has a Detector.dat and Efficiency.csv file in it
      Throws exception in error; returned detector should always be valid.
   */
  static std::shared_ptr<DetectorPeakResponse> initAGadrasDetectorFromDirectory(
                              const std::string &directory,
                              const bool allow_missing_efficiency_csv = false );

  static std::shared_ptr<DetectorPeakResponse> initARelEffDetector( const SpecUtils::DetectorType type,
                                                                   std::string manufacturer,
                                                                   std::string model );
  
  /** Per-caller additions to the "Use DRF?" dialog #createChooseDrfDialog builds.
   
   Exists so one caller (the app-URL import, which offers detector-modeling choices) can add its
   own control and act on the result, without every other caller of that dialog growing the same
   parameters.  All three are optional.
   */
  struct ChooseDrfHooks
  {
    /** Called with the dialog's content area, after the DRF description and before the
     "use as default" checkboxes, to add caller-specific controls.
     */
    std::function<void(Wt::WContainerWidget *)> addContent;
    
    /** Called with the chosen DRF when the user accepts, *before* it is written to the database or
     applied to the session - so a hook that modifies the DRF (and therefore its hash) is what gets
     stored and used.
     */
    std::function<void(std::shared_ptr<DetectorPeakResponse>)> beforeAccept;
    
    /** Called with the chosen DRF after it has been stored and applied - for a follow-up that
     needs the detector already in use, such as opening another tool on it.
     */
    std::function<void(std::shared_ptr<DetectorPeakResponse>)> afterAccept;
  };//struct ChooseDrfHooks

  
  
  /** Creates a dialog to allow the user to select to use (one of) the passed in DRFs.
   
   @param drfs The DRFs to allow the user to choose from; if empty (or only null or invalid drfs),
          no dialog will be shown.  If a single DRF, then user wont be given a choice.
   @param mainMsgHtml The HTML to show as the description for the DRFs; this appears above the
          drf selection combo-box (if more than one DRF)
   @param creditsHtml The HTML to show below the combo-box.
   @param saveDrfsCallBack If non-empty, a checkbox to allow saving DRFs will be shown, and if
          checked when the user accepts dialog, will call this callback.
   @param hooks Optional per-caller additions; see #ChooseDrfHooks.
   */
  static void createChooseDrfDialog( std::vector<std::shared_ptr<DetectorPeakResponse>> drfs,
                                    Wt::WString mainMsgHtml,
                                    std::string creditsHtml,
                                    std::function<void()> saveDrfsCallBack = nullptr,
                                    ChooseDrfHooks hooks = ChooseDrfHooks() );
  
  /** Checks if file at passed in path is a TSV/CSV file that contains
   coefficients for the exp( c0 + c1*logx + c2*logx^2 + ...) equation.
   If so, returns detector.  If not, returns nullptr.

   ToDo: Currently reads in most of the information exported in the CSV from MakeDrf tool, but maybe not all.
   */
  static std::shared_ptr<DetectorPeakResponse> parseInterSpecRelEffCsvFile( const std::string filePath );

  /** Searches for potential relative efficiency detector files in static and user data directories.
   Returns a vector of file paths to potential RelEff detector definition files.
   */
  static std::vector<std::string> potential_rel_eff_det_files();

  /** Recursively searches a directory for GADRAS DRF directories.
   Returns a vector of directory paths that contain both Efficiency.csv and Detector.dat files.
   @param sourcedir The directory to search recursively
   */
  static std::vector<std::string> recursive_list_gadras_drfs( const std::string &sourcedir );

  /** Searches for potential GADRAS detector directories.
   Returns a vector of directory paths that may contain GADRAS DRF subdirectories.
   Note: The returned paths may include parent directories that contain many sub-directories
         with different detector efficiencies. Use recursive_list_gadras_drfs() to find
         the actual DRF directories within these paths.
   @param interspec Pointer to InterSpec instance; may be nullptr (will skip user preferences if null)
   */
  static std::vector<std::string> potential_gadras_det_dirs( InterSpec *interspec );

  /** The standalone FWHM-from-foreground and Monte-Carlo tools finished: record the detector they
   made, and close. */
  void handleFitFwhmFinished( std::shared_ptr<DetectorPeakResponse> drf );
  void handleMcResponseFinished( std::shared_ptr<DetectorPeakResponse> drf );

  /** Opens (or re-shows) the "Customize" dialog (rename, geometry + Monte-Carlo characterization,
   FWHM, uncertainty), seeded with the previewed detector - or, on the Import tab, with the file
   being imported (for a geometry-only file, to characterize it).
   */
  void openModifyWindow();

  /** Receives the customized detector from the Customize dialog and uses it. */
  void handleModifyFinished( std::shared_ptr<DetectorPeakResponse> drf );

  /** The "Customize" dialog this widget owns, or nullptr when it is not open.
   \sa InterSpec::drfModifyWidget
   */
  DrfModifyWindow *modifyWindow();

  /** Closes the "Customize" dialog without recording an undo/redo step - for use *from* an
   undo/redo step, whose counterpart already covers the close. */
  void programmaticallyCloseModifyWindow();

protected:
  enum class Tab { Library, Import, Locations };

  void showTab( const Tab tab );

  /** Re-reads the Library (off the session thread; cheap when nothing changed). */
  void startScan();
  void receiveSnapshot( std::shared_ptr<const DrfLibrary::Snapshot> snapshot, const bool final );

  /** The user picked a Library entry: previews it (see #selectEntry), and stops the preview
   following the detector in use. */
  void pickEntry( const std::string &entryId );

  /** Previews a Library entry (empty to clear); loads its full detector off the session thread. */
  void selectEntry( const std::string &entryId );
  const DrfLibrary::Entry *pickedEntry() const;

  /** The Library entry for a detector, by its stored hash; empty if none. */
  std::string entryIdFor( const std::shared_ptr<const DetectorPeakResponse> &drf ) const;

  void updateInUseStrip();
  void updatePreview();
  void updatePrimaryButton();
  void updateChips();
  void updateFilter();
  void updateAutoPick();

  void handlePrimary();

  /** On a phone: shows (or hides) the Details screen - the preview, full width, with a back button. */
  void setDetailsShown( const bool show );
  void handleImportChanged();
  void showInUseInList();

  /** After using a detector that may be new to the Library (imported, or customized): switch to
   the Library, with filters cleared, and pick it once it is listed. */
  void showUsedInLibrary();
  void confirmDeletePicked();
  void deleteEntry( const long long dbId, const std::string &name );
  void showQrCode();
  void showNotice( const Wt::WString &message );

  // Library › Locations (desktop builds)
  void updateLocationsPage();
  void setLocationShown( const std::string &id, const bool shown );
  void rescanLocation( const std::string &id );
  /** Asks before removing a location the user added (its files are not touched). */
  void confirmRemoveLocation( const std::string &id );
  void removeLocation( const std::string &id );
  void locateLocation( const std::string &id );
  void startAddLocationFolder();
  void startAddLocationFile();
  /** Adds a folder or file to the Library (saved in the preferences), and re-reads it.  Returns
   whether it is now in the Library (false if it could not be added). */
  bool addUserLocation( const std::string &path );

  /** The detector the Import tab would use now (a file's, a detector found in a dropped folder, or
   the formula's), or null. */
  std::shared_ptr<DetectorPeakResponse> importCandidate() const;

  /** For a geometry-only file on the Import tab, the detector to characterize; else null. */
  std::shared_ptr<DetectorPeakResponse> importSeed() const;

  // A drop of several detectors' files (a folder)
  void handleBatchDropped( std::shared_ptr<const DrfDropBatch> batch );
  void receiveFound( std::shared_ptr<const std::vector<DrfImport::FoundDetector>> found,
                     const size_t numSkipped );
  void pickFound( const size_t index );
  void clearFound();
  /** Copies the dropped detectors into the user-data `drfs` folder (so the Library lists them), and
   uses the picked one once the Library has read them. */
  void keepFoundAndUse();
  /** After a scan, uses the Library entry kept by #keepFoundAndUse (see #m_useAfterScanRoot). */
  void useKeptEntry();

  InterSpec *m_interspec;
  std::shared_ptr<DataBaseUtils::DbSession> m_sql;

  /** The detector in use (follows `InterSpec::detectorChanged()`). */
  std::shared_ptr<DetectorPeakResponse> m_inUse;
  /** Whether "Use" was pressed in this dialog for #m_inUse (whose hash #m_justUsedHash is). */
  bool m_justChanged;
  uint64_t m_justUsedHash;

  Tab m_tab;

  /** A phone (or narrow window): Library, Import and Details are separate screens. */
  bool m_narrow;
  bool m_showDetails;
  Wt::WContainerWidget *m_detailsHeader;
  Wt::WPushButton *m_detailsBack;
  /** Phone: the bar at the bottom holding the main button and the auto-pick boxes (else null). */
  Wt::WContainerWidget *m_actionBar;
  std::shared_ptr<const std::vector<DrfLibrary::Entry>> m_entries;
  std::string m_pickedId;
  /** Whether the user picked a detector (otherwise the one in use is picked when the Library is read). */
  bool m_userPicked;
  /** Pick the detector in use once the next scan finishes (it may have only just been saved). */
  bool m_pickInUseAfterScan;
  bool m_scrolledToPicked;
  std::shared_ptr<const DetectorPeakResponse> m_pickedDrf;
  Wt::WString m_pickedError;
  std::string m_suggestedNote;

  uint64_t m_scanGeneration, m_loadGeneration;
  std::shared_ptr<std::atomic<bool>> m_scanCancel;
  bool m_scanning;

  /** Library filter state. */
  int m_originChip;  ///< -1 for "All", else a `DrfLibrary::Origin`
  bool m_withGeometry;

  // The "In use now" strip
  Wt::WContainerWidget *m_stripInUse, *m_stripNone;
  Wt::WText *m_stripIcon, *m_stripJust, *m_stripName, *m_stripMeta;

  // Tabs
  Wt::WPushButton *m_libTab, *m_importTab;
  Wt::WStackedWidget *m_stack;

  // Library
  DrfLibraryList *m_list;
  Wt::WLineEdit *m_search;
  Wt::WTimer *m_searchTimer, *m_checkingTimer;
  std::vector<std::pair<int,Wt::WPushButton *>> m_chips;
  Wt::WPushButton *m_geoChip;

  // Import
  Wt::WPushButton *m_fileModeBtn, *m_formulaModeBtn;
  Wt::WStackedWidget *m_importStack;
  DrfImportWidget *m_importWidget;
  Wt::WContainerWidget *m_formats;  ///< the "Recognized formats" tiles
  DrfFoundList *m_foundList;
  DrfFormulaForm *m_formula;

  /** The detectors found in the last drop of several files, and the one picked (npos for none). */
  std::shared_ptr<const DrfDropBatch> m_batch;
  std::shared_ptr<const std::vector<DrfImport::FoundDetector>> m_found;
  size_t m_foundPick;
  std::string m_foundFolder;   ///< the dropped folder's name, or empty for loose files
  uint64_t m_batchGeneration;

  /** Set by #keepFoundAndUse: the folder the detectors were copied to, the picked one's name, and
   what to use if the Library does not list it after all. */
  std::string m_useAfterScanRoot, m_useAfterScanName, m_useAfterScanPath;
  std::shared_ptr<DetectorPeakResponse> m_useAfterScanFallback;
  /** Use this Library entry as soon as its detector is loaded. */
  std::string m_useWhenLoadedId;
  /** On a phone, using a detector closes the dialog - once it is used, for a kept folder's (which is
   used only after the Library has read it). */
  bool m_closeAfterUse;

  DrfPreviewPanel *m_preview;

  // Locations
  Wt::WPushButton *m_manageLocations;
  DrfLocationsPage *m_locationsPage;
  /** The locations the last scan read, and what it found in each. */
  std::shared_ptr<const std::vector<DrfLibrary::Location>> m_scanLocations;
  std::shared_ptr<const std::vector<DrfLibrary::LocationStatus>> m_locationStatuses;

  // Footer
  Wt::WContainerWidget *m_autoPick;
  Wt::WCheckBox *m_autoSerial, *m_autoModel;
  Wt::WPushButton *m_cancelBtn, *m_primaryBtn;

  Wt::WContainerWidget *m_notice;
  Wt::WText *m_noticeText;

  /** The "Customize" dialog, while open (auto-nulls on close). */
  Wt::Core::observing_ptr<DrfModifyWindow> m_modifyWindow;

  Wt::Signal<> m_doneSignal;
};//class DrfSelect


class DrfSelectWindow : public AuxWindow
{
  friend class AuxWindow;

public:
  virtual ~DrfSelectWindow();

  DrfSelect *widget();

protected:
  DrfSelectWindow( InterSpec *viewer );

  DrfSelect *m_edit;
  InterSpec *m_interspec;
};//class DrfSelectWindow


#endif //DrfSelect_h
