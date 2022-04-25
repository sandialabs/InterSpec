#ifndef SpectrumViewer_h
#define SpectrumViewer_h
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
#include <deque>
#include <memory>
#include <vector>

#include <Wt/Dbo/Dbo>
#include <Wt/WContainerWidget>

//#include "SpecUtils/SpecFile.h"
//Without including InterSpecUser.h here, we get some wierd issues with the
//  DB optomistic versioning...
#include "InterSpec/InterSpecUser.h"

class PeakDef;
class PopupDiv;
class AuxWindow;
class GoogleMap;
class PeakModel;
class MaterialDB;
class D3TimeChart;
class DecayWindow;
struct ColorTheme;
class UserFileInDb;
class PopupDivMenu;
class EnergyCalTool;
class SpectrumChart;
class UseInfoWindow;
class WarningWidget;
class TerminalWidget;
class WarningMessage;
class PeakEditWindow;
class PeakInfoDisplay;
class SpecMeasManager;
class GammaCountDialog;
class PopupDivMenuItem;
class SpectraFileHeader;
class PopupWarningWidget;
class D3SpectrumDisplayDiv;
class FeatureMarkerWindow;
class DetectorPeakResponse;
class IsotopeSearchByEnergy;
class ShieldingSourceDisplay;
class EnergyCalPreserveWindow;
class SimpleNuclideAssistPopup;
class ReferencePhotopeakDisplay;
class LicenseAndDisclaimersWindow;
namespace D3SpectrumExport{ struct D3SpectrumChartOptions; }

#if( INCLUDE_ANALYSIS_TEST_SUITE )
class SpectrumViewerTester;
#endif

namespace SpecUtils{ class SpecFile; }
namespace SpecUtils{ class Measurement; }
namespace SpecUtils{ enum class SpectrumType : int; }
namespace SpecUtils{ enum class DetectorType : int; }
namespace SpecUtils{ enum class OccupancyStatus : int; }


namespace DataBaseUtils
{
  class DbSession;
}

namespace Wt
{
  class WText;
  class WMenu;
  class WMenuItem;
  class WTextArea;
  class WLineEdit;
  class WTabWidget;
  class WGridLayout;
  class WPopupMenu;
  class WPushButton;
  class WSuggestionPopup;
    
  namespace Dbo
  {
    template<class Ch> class ptr;
  }
}//namespace Wt

/** A convience function to call #InterSpec::logMessage
 
 \deprecated Please directly call InterSpec::instance()->logMessage(message,source,priority).
 */
void log_error_message( const std::string &message, const std::string &source, const int priority );


enum class FeatureMarkerType : int
{
  /// \TODO: move this to some other header
  EscapePeakMarker,
  ComptonPeakMarker,
  ComptonEdgeMarker,
  SumPeakMarker,
  NumFeatureMarkers
};//enum FeatureMarkerType


class InterSpec : public Wt::WContainerWidget
{
  /****************************************************************************\
  | This class is the overall architecture for the InterSpec application.      |
  | It contains most of the widgets for the visible app and handles most of    |
  | the actions taken by the user that are not within a specific tool.         |
  |
  \****************************************************************************/

public:
  InterSpec( Wt::WContainerWidget *parent = 0 );

  ~InterSpec() noexcept(true);
    
  /** Returns the InterSpec instance corresponding to the current WApplication instance.
   Will return nullptr if WApplication::instance() is null (e.g., current thread is not in a
   WApplication event loop).
   */
  static InterSpec *instance();
  
  /** Sets the directory InterSpec "data" files are located, including cross
      sections, materials, detector response functions, nuclear decay info
      and similar.
   
      Must be called before creating any InterSpec class instances (this doesnt
      have enforcement implemented, but behaviour is undefined otherwise).
   
      Can be an absolute path or relative; if relative. it is to the CWD.
      Default is "data".
   
      Alls sets appropriate files or directories for DecayDataBaseServer,
      ReactionGammaServer, and MassAttenuation.  Other tools call
      InterSpec::dataDirectory() to find the appropriate directory.
   */
  static void setStaticDataDirectory( const std::string &dir );
 
  /** Directory where files like cross sections, materials, detector response
      function, nuclear decay information, and similar are stored.
   */
  static std::string staticDataDirectory();
 
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER )
  /** Sets the directory were we can write write the user preference database
   file (if using sqlite3); also the location will search for extra detector
   response functions.
   
   Will throw exception if not empty string and its an invalid directory.
   */
  static void setWritableDataDirectory( const std::string &dir );
  
  /** Returns the location you can write files to, such as user preferences.
      Will throw exception if hasnt been set (or set with empty string)
   */
  static std::string writableDataDirectory();
#endif  //if( not a webapp )
  
  /** Function called by Wt when the client reports a change in widget size. */
  virtual void layoutSizeChanged( int width, int height );

  //findAndSetExcludedSamples(...): excludes background and calibration samples
  //  in m_dataMeasurement - placing excluded sample numbers into
  //  m_excludedSamples.
  void findAndSetExcludedSamples( std::set<int> definetly_keep_samples );

#if( SpecUtils_ENABLE_D3_CHART )
  //print_d3_json(): Output data from current user chart data to JSON
  //  format to be used for running D3.js HTML files.
  std::string print_d3_json() const;
  
  //print_d3_reference_gammas(): Output data from current reference gammas
  //  displayed on the chart to JSON format to be used for rendering
  //  reference gammas in D3.js HTML files.
  std::string print_d3_reference_gammas() const;
  
  //getD3SpectrumOptions(): Output current chart options to a struct
  //  to provide the same "state" of the current user session to
  //  D3.js HTML files.
  D3SpectrumExport::D3SpectrumChartOptions getD3SpectrumOptions() const;
#endif //#if( SpecUtils_ENABLE_D3_CHART )
  
  
  //userOpenFileFromFilesystem(): appropriate to call if the user opens a file
  //  by double-clicking it in the finder.  Will check if the file might be
  //  a background to the current foreground, and if so, prompt the user how
  //  they would like to open it; if the file probably isnt a background to
  //  the current foreground, then the files is opened as foreground.  Returns
  //  status of if the file could be parsed or not.
  bool userOpenFileFromFilesystem( const std::string filepath, std::string displayFileName = "" );
  
  //promptUserHowToOpenFile(): method called by userOpenFileFromFilesystem()
  //  and when the file could possibly be a background to the current foreground.
  void promptUserHowToOpenFile( std::shared_ptr<SpecMeas> meas,
                                 std::shared_ptr<SpectraFileHeader> header );
  
  //finishLoadUserFilesystemOpenedFile(): is used by promptUserHowToOpenFile()
  //  if the user doesnt cancel the the opening of the file.
  void finishLoadUserFilesystemOpenedFile( std::shared_ptr<SpecMeas> meas,
                                           std::shared_ptr<SpectraFileHeader> header,
                                           const SpecUtils::SpectrumType type );

  // Options for setting the spectrum.
  //  E.g., what we should potentially prompt the user for.  Like if were opening up a spectrum file
  //    that belongs to the same detector as previous, then should we ask if the old energy
  //    calibration should be preserved.  However a lot of the times we're calling #setSpectrum, we
  //    know we shouldnt bug the user, so we dont want these things checked.
  enum SetSpectrumOptions
  {
    CheckToPreservePreviousEnergyCal = 0x01,
    CheckForRiidResults = 0x02
    
    // TODO: it seems both these options are the same everywhere - maybe go back and just use a bool
  };//SetSpectrumOptions
  
  
  //setSpectrum(...) is intended to be the only place m_dataMeasurement,
  //  m_secondDataMeasurement, or m_backgroundMeasurement are set.
  //  If sample_numbers is empty, then it will default to all sample numbers
  //  except in the case of a passthrough SpecUtils::SpectrumType::Foreground, then background and
  //  calibration samples will not be used (in order to show background/calib
  //  samples in passthrough foreground files, you must explicity pass these
  //  sample numbers in sample_numbers)
  //If you are passing in a new SpecUtils::SpectrumType::Foreground spectrum that has the same number
  //  of bins as the previous spectrum, and the new SpecMeas doesnt have a valid
  //  detector, then the new SpecMeas object will be assigned a pointer to the
  //  old detector (they will now share them) - there may be issues with
  //  the SpecMeasManager clearing some measurments out of memorry so spectra
  //  no longer share detector objects.
  //
  //emits the m_displayedSpectrumChangedSignal
  void setSpectrum( std::shared_ptr<SpecMeas> meas,
                    std::set<int> sample_numbers,
                    const SpecUtils::SpectrumType spec_type,
                    const Wt::WFlags<SetSpectrumOptions> options = 0 );

  //reloadCurrentSpectrum(...): reloads the specified spectrum.  This function
  //  is useful when you change teh SpecMeas object (e.x. live or real time),
  //  and would like to propogate this to to the GUI componenets
  void reloadCurrentSpectrum( SpecUtils::SpectrumType spec_type );

  /** Loads either the user-preferred DRF (e.g., they checked a box to always use a DRF for a
   specific serial number or model), or a DRF based on DetectorType/manufacturer/model from
   DRFs available in either the static and user data directories.
   
   The modified() and modified_since_decode() statuses of meas will not be changed by loading of a
   DRF.
   */
  void loadDetectorResponseFunction( std::shared_ptr<SpecMeas> meas,
                                     const SpecUtils::DetectorType type,
                                     const std::string serial_number,
                                     const std::string manufacturer,
                                     const std::string model,
                                     const bool tryDefaultDrf,
                                     const std::string sessionId );

  
  /** Handles "deep" urls.
   
   Meant to handle URLs with the scheme "interspec://...", that would be passed to the application
   by the OS, like when a QR code is scanned.
   
   The prefix "interspec://" is optional, and may be omitted.
   
   An example URL is "interspec://drf/specify?v=1"
   
   Throws std::exception if url cant be used.
   
   As of 20220405: Implementation and scope of use still being fleshed out.
   */
  void handleAppUrl( std::string url );
  
  
  
  //For the 'add*Menu(...)' functions, the menuDiv passed in *must* be a
  //  WContainerWidget or a PopupDivMenu
  void addFileMenu( Wt::WWidget *menuDiv, const bool isAppTitlebar );
  void addDisplayMenu( Wt::WWidget *menuDiv );
  void addDetectorMenu( Wt::WWidget *menuDiv );
  void addToolsMenu( Wt::WWidget *menuDiv );
  void addPeakLabelSubMenu( PopupDivMenu *parentWidget );
  void addAboutMenu( Wt::WWidget *menuDiv );

  //addPeak(): Adds a new peak to the peak model, and returns the models index
  //  of the new peak. If associateShownNucXrayRctn is specified true _and_ the
  //  user is showwing some reference gamma lines, than the new peak will be
  //  assigned to be from the shown lines if and are reasonably close.
  //  If the returned WModelIndex is not valid, then the peak was not added.
  Wt::WModelIndex addPeak( PeakDef peak, const bool associateShownNucXrayRctn );
  
  
  Wt::WContainerWidget *menuDiv();

  const std::set<int> &displayedSamples( SpecUtils::SpectrumType spectrum_type ) const;
  
  /** Foreground samples that might be displayed.
   
   TODO: Eliminate calling this function, in coordination with upgrading
         passthroughTimeToSampleNumber() (see TODOs for that function)
   */
  std::set<int> validForegroundSamples() const;
  
  std::vector<std::string> detectorsToDisplay( const SpecUtils::SpectrumType type ) const;

  
  std::shared_ptr<SpecMeas> measurment( SpecUtils::SpectrumType spectrum_type );
  std::shared_ptr<const SpecMeas> measurment( SpecUtils::SpectrumType spectrum_type ) const;

  /** Returns the spectrum displayed on the client side; will return a nullptr
      if the spectrum isnt displayed.
  */
  std::shared_ptr<const SpecUtils::Measurement> displayedHistogram( SpecUtils::SpectrumType spectrum_type ) const;

  /** Returns the sample numbers marked as Foreground, Background, or Secondary (i.e., intrinsic or
   known), in the foreground spectrum file.
   This is used for displaying the highlighted regions on the time series chart.
   */
  std::set<int> sampleNumbersForTypeFromForegroundFile(
                                                const SpecUtils::SpectrumType spec_type ) const;
  
#if( IOS )
  void exportSpecFile();
#endif
  
#if( USE_DB_TO_STORE_SPECTRA )
  //measurmentFromDb(...): returns the measurment that has been, or will be
  //  serialized to the database.  If 'update' is false, then just the last
  //  serialization will be returned and in fact may be null.  If 'update'
  //  is true, then the measurment will be saved to the database first (unless
  //  it hasnt been modified since last saving) and then returned; if the user
  //  preference is to not save spectra to the database, then it will not be
  //  be saved.  In the case the user preference is to not save to the database,
  //  but the file is already in there, but has been modified in memmorry, then
  //  the file will not be re-saved to the database.
  //  Function wont throw, but may return a null pointer
  Wt::Dbo::ptr<UserFileInDb> measurmentFromDb( SpecUtils::SpectrumType type, bool update );
  
  //saveStateToDb( entry ): saves the application state to the pointer passed
  //  in.  Note that pointer passed in must be a valid pointer associated with
  //  this m_user.  If entry->stateType is specified for testing, then every
  //  thing will be copied and write protected so it cant be changed in the
  //  future. May throw exception.  std::runtime_error exceptions will contain
  //  messages that are reasonably okay to send to users.
  void saveStateToDb( Wt::Dbo::ptr<UserState> entry );
  
  //serializeStateToDb(...): a convience function for saveStateToDb(...).
  Wt::Dbo::ptr<UserState> serializeStateToDb( const Wt::WString &name,
                                              const Wt::WString &desc,
                                              const bool forTesting,
                                              Wt::Dbo::ptr<UserState> parent);
  
  //loadStateFromDb( entry ): sets the applications state to that in 'entry'
  void loadStateFromDb( Wt::Dbo::ptr<UserState> entry );

  //testLoadSaveState(): a temporary function to help develop the loading and
  //  saving of the applications state
//  void testLoadSaveState();
  
    void stateSave();
    void stateSaveAs();
    void stateSaveTag();
    void stateSaveAsAction( Wt::WLineEdit *nameedit,
                         Wt::WTextArea *descriptionarea,
                    AuxWindow *window,
                         const bool forTesting);
    void stateSaveTagAction(Wt::WLineEdit *nameedit,
                          AuxWindow *window);
  
  //startStoreStateInDb(...): If state hasnt been saved yet (m_currentStateID>0)
  //  then a dialog will be popped up allowing user to enter a title and
  //  description.
  //  If the state is already in the database, and either the state is not
  //  readonly, or 'allowOverWrite' is specified true, then the database will be
  //  updated.
  //  If 'asNewState' is specified, then even if the current state is already in
  //  the database, the process (e.g. user dialog) for saving a new state in the
  //  database will be started.
  //  Specifying 'forTesting' will cause the for testing flag to be set, as well
  //  as for the state to be marked as read only.
  void startStoreStateInDb( const bool forTesting,
                            const bool asNewState,
                            const bool allowOverWrite,
                            const bool endOfSessionNoDelete );
  
  void finishStoreStateInDb( Wt::WLineEdit *name, Wt::WTextArea *description,
                            AuxWindow *window, const bool forTesting , Wt::Dbo::ptr<UserState> parent);
  void browseDatabaseStates( bool testStates );
  
#if( INCLUDE_ANALYSIS_TEST_SUITE )
  void startStateTester();
  
  //startStoreTestStateInDb(): a convience function that makes a dialog to give
  //  the user an option to overwrite their current state or create a new one.
  void startStoreTestStateInDb();
  
  //storeTestStateToN42(): stores foreground, background, and shielding/source
  //  model into a since 2012 N42 file.  Does not throw, but will notify the
  //  user via the GUI upon error.
  void storeTestStateToN42( std::ostream &output,
                            const std::string &name, const std::string &descr );
  
  //loadTestStateFromN42(): Attempts to load a state previously saved to an
  //  XML file via storeTestStateToN42
  void loadTestStateFromN42( const std::string filename );
  
  //startN42TestStates(): creates dialog that lists files ending in '.n42' in
  //  the 'analysis_tests' folder so the user can select and then load one.
  void startN42TestStates();
#endif
  
  
#endif //#if( USE_DB_TO_STORE_SPECTRA )
  

  /** Converts the current chart to a PNG and triggers a download of the
     resulting file.
   
   Note: Currently not fully working for D3 based spectrum chart.
   
   @param spectrum If true, make a PNG or SVG for the spcetrum chart.  If false make
          PNG or SVG for the time-series chart
   @param asPng if true, save as a PNG.  If false, save as a SVG
   
   Note: currently only PNG is supported for time chart
   */
  void saveChartToImg( const bool spectrum, const bool asPng );
  
  //displayScaleFactor(...): the live time scale factor used by to display
  //  the histogram returned by displayedHistogram( spectrum_type ).
  double displayScaleFactor( SpecUtils::SpectrumType spectrum_type ) const;
  
  
  //setDisplayScaleFactor(...): sets the display scale factor. SpecUtils::SpectrumType
  //  can not be SpecUtils::SpectrumType::Foreground (an exception will be thrown).
  //This will make it so the effective live time of the spectrum will be
  //  'sf' times the original live time of the spectrum.
  void setDisplayScaleFactor( const double sf,
                              const SpecUtils::SpectrumType spectrum_type );
  
  //changeDisplayedSampleNums(...): called by both changeTimeRange() and
  //  sampleNumbersAddded() when the user wants to change displayed sample numbers
  void changeDisplayedSampleNums( const std::set<int> &samples,
                                 const SpecUtils::SpectrumType type );

  
  //liveTime(): returns the live time of the desired spectrum.  Will return 0.0
  //  if not defined.
  float liveTime( SpecUtils::SpectrumType spectrum_type ) const;
  
  PopupDivMenu *displayOptionsPopupDiv();

  /** Emits a message to the user.
   
    The priority corresponds to the WarningWidget::WarningMsgLevel enum, namely:
   - Info: 0
   - Low: 1
   - Medium: 2
   - High: 3
   - Save: 4
   
   \sa log_error_message
   */
  void logMessage( const Wt::WString& message, const Wt::WString& source, int priority );
  
  virtual Wt::Signal< Wt::WString, Wt::WString, int > &messageLogged();
  
  void toggleToolTip( const bool showToolTips );
  
  /** Width of the apps window, in pixels. Will return a value of zero if (not yet) available. */
  int renderedWidth() const;
  /** Height of the apps window, in pixels. Will return a value of zero if (not yet) available. */
  int renderedHeight() const;

  //ToDo: add a signal that fires when the app size changes, to allow adjusting AuxWindow's and stuff.

  
  //displayedSpectrumChanged(): is a singal emitted whenever a spectrum or
  //  sample numbers of a spectrum is changed.
  //  Note that this signal is _not_ emitted if only the display scale factor
  //  is changed.
  Wt::Signal<SpecUtils::SpectrumType,                //which spectrum changed
             std::shared_ptr<SpecMeas>, //A pointer to the new spectrum
             std::set<int>,             //The sample numbers displayed
             std::vector<std::string>   //The detectors displayed
             >& displayedSpectrumChanged();

  /** Signal emmitted when user changes a spectrums scale factor in the "Spectrum Files" tab, or using graphical y-axis scalers.
   Is not emmitted when new spectra or sample numbers are loaded/changed.
   */
  Wt::Signal<SpecUtils::SpectrumType,double> &spectrumScaleFactorChanged();
  
  
  //addHighlightedEnergyRange(): Adds highlighted range to the energy spectrum.
  //  Returns the ID of the highlight region, which you will need to remove
  //  the highlight region.
  size_t addHighlightedEnergyRange( const float lowerEnergy,
                                   const float upperEnergy,
                                   const Wt::WColor &color );
  //removeHighlightedEnergyRange(): removes the region cooresponding to the ID
  //  passed in from the energy chart.  Returns succesfulness of the removing
  //  the region.
  bool removeHighlightedEnergyRange( const size_t regionid );
  
  //setDisplayedEnergyRange(): sets the displayed energy range that should be
  //  shown.  Range is not checked for validity. E.g. you should not ask to zoom
  //  in to less than one bin.  You can also specify a lower or higher energy
  //  than data is available for.  Definitely dont input the same energy for both
  //  values, or bad things will happen.
  void setDisplayedEnergyRange( float lowerEnergy, float upperEnergy );
  
  /** Sets the Y-axis range.
   
   Returns true if successful, or false if the request couldnt be fully honored
   (e.g., a negative lower counts was specified, but its currently log-scale)
   */
  bool setYAxisRange( float lower_counts, float upper_counts );
  
  //displayedSpectrumRange(): Grab the ranges in y and y that are currently
  //  displayed to the user.
  void displayedSpectrumRange( double &xmin, double &xmax, double &ymin, double &ymax ) const;
  
  void handleShiftAltDrag( double lowEnergy, double upperEnergy );
  
  void setToolTabsVisible( bool show );
  bool toolTabsVisible() const;
  
  void showMakeDrfWindow();
  void showDrfSelectWindow();
  void showCompactFileManagerWindow();
  
  //Nuclide Search
  void closeNuclideSearchWindow();
  void showNuclideSearchWindow();
  
  
  //initMaterialDbAndSuggestions(): initializes m_materialDB and
  //  m_shieldingSuggestion, and posts to the WServer threadpool a call to
  //  fillMaterialDbAndPushSuggestionsToUsers() which does the actual parsing
  //  of the material database.  Its done in two stages so as to not slow
  //  the initial loading of the app; it would be nice to not have to parse the
  //  materia DB file in the event loop at all, but this makes it easy to avoid
  //  race conditions, or whatever.
  void initMaterialDbAndSuggestions();
  

  //fillMaterialDb(): populates materialDB, and then calls the provided
  //  'update' function by posting to the WServer thread pool for 'sessionid'
  //  so it will be executed in its event loop.
  static void fillMaterialDb( std::shared_ptr<MaterialDB> materialDB,
                              const std::string sessionid,
                              boost::function<void(void)> update );
  
  //pushMaterialSuggestionsToUsers(): should be called from application loop, to
  //  fill m_shieldingSuggestion (from m_materialDB) and then push to the user.
  void pushMaterialSuggestionsToUsers();
  
  void showGammaXsTool();
  void showDoseTool();
  void showShieldingSourceFitWindow();
  void closeShieldingSourceFitWindow();
  void saveShieldingSourceModelToForegroundSpecMeas();
  

  //Note: I am not a big fan of how I treat m_referencePhotopeakLines in this class,
  //      but it is how it is due to technical difficulties that I dont want
  //      to figure out how to fix - right now I'm willing to accept the
  //      bandwidth and rendering in-efficiencies

  //When showGammaLinesWindow() is called, m_referencePhotopeakLines is deleted,
  //  all displayed photopeak lines cleared; an AuxWindow is the created, as
  //  well as a new ReferencePhotopeakDisplay (which is assigned to
  //  m_referencePhotopeakLines).
  void showGammaLinesWindow();

  //When closeGammaLinesWindow(...) is called, the current m_referencePhotopeakLines
  //  is deleted after clearing all lines, the a new m_referencePhotopeakLines
  //  is created and added to the tabs widget at bottom tool tabs.
  void closeGammaLinesWindow();

  //handleToolTabChanged(...): sets focus to isotope line edit
  void handleToolTabChanged( int tabSwitchedTo );


  SpecMeasManager *fileManager();
  
  PeakModel *peakModel();

  /** The material database.
   
   Object will be alive as long as *this
   */
  MaterialDB *materialDataBase();
  
  /** The suggestion pop-up widget for shielding names; used globally for all shielding name inputs
   so that there is not duplicate copies of the widget in the DOM.
   
   Object will be alive as long as *this.
   */
  Wt::WSuggestionPopup *shieldingSuggester();
  
  
  //detectorChanged(): signal emited when the detector is changed to a
  //  completely new detector.  Note that the object pointed to stays the same
  //  but it has been totally redifined.
  Wt::Signal<std::shared_ptr<DetectorPeakResponse> > &detectorChanged();

  //detectorModified(): signal emited when a property of the current detector
  //  object is modified.
  Wt::Signal<std::shared_ptr<DetectorPeakResponse> > &detectorModified();

  void showEnergyCalWindow();
  void handEnergyCalWindowClose();

  EnergyCalTool *energyCalTool();
  
  
  void showWarningsWindow();
  void handleWarningsWindowClose( bool closeWindowTo );

  void showPeakInfoWindow();
  void handlePeakInfoClose();

#if( USE_GOOGLE_MAP )
  void createMapWindow( SpecUtils::SpectrumType spectrum_type );
  void displayOnlySamplesWithinView( GoogleMap *map,
                                     const SpecUtils::SpectrumType targetSamples,
                                     const SpecUtils::SpectrumType fromSamples );
#endif
  
#if( USE_SEARCH_MODE_3D_CHART )
  void create3DSearchModeChart();
#endif

  /** Show the RIID results included in the spectrum file. */
  void showRiidResults( const SpecUtils::SpectrumType type );
  
#if( USE_TERMINAL_WIDGET )
  void createTerminalWidget();
  void handleTerminalWindowClose();
#endif

  /** Will show the disclaimer, license, and statment window, setting
      m_licenseWindow pointer with its value.
      If is_awk is true, the confirm button will say "Awcknowledge" and the
      window wont be otherwise closeable; if false button will say "close" and
      the close icon in top bar will be avaialble.
   */
  void showLicenseAndDisclaimersWindow( const bool is_awk,
                                      std::function<void()> finished_callback );
  
  /** Brings up a dialog asking the user to confirm starting a new session, and if they select so, will start new session. */
  void startClearSession();
    
  /** Deletes disclaimer, license, and statment window and sets m_licenseWindow
      to nullptr;
   */
  void deleteLicenseAndDisclaimersWindow();
  
  //showWelcomeDialog(...): shows the welcome dialog.  If force==false, and the
  //  user has opted to not show the dialog in a previous display, then it
  //  wont be shown; if force==true, than it is shown no matter what and the
  //  "do not show again" check box wont be shown to effect users prefference.
  //  After creation m_useInfoWindow will point to the dialog.
  void showWelcomeDialog( bool force = false );

  //deleteWelcomeCountDialog(): deletes the currently showing dialog (if its
  //  showing), and sets m_useInfoWindow to null
  void deleteWelcomeCountDialog();
  
  //deleteEnergyCalPreserveWindow(): deletes the currently showing diaolg (if its
  //  showing), and sets m_preserveCalibWindow to null
  void deleteEnergyCalPreserveWindow();
  
  //showIEWarningDialog(): returns NULL if user previously specified to not show
  //  again, otherwise it returns the AuxWIndow it is displaying.  The dialog
  //  is by default shown visible, and will deleted when user is done with it.
  AuxWindow *showIEWarningDialog();
  
  // The user itself gets to be public--no need to protect access to it.
  //Note 20130116: m_user should be made protected, but really the whole
  //  preference thing should be re-done, see README
  Wt::Dbo::ptr<InterSpecUser> m_user;
  
  //sql returns the DbSession (which holds the Wt::Dbo::Session) associated
  //  with m_user.  The reason for using an indirection via
  //  DataBaseUtils::DbSession is to ensure thread safety when interacting with
  //  the database.  A shared pointer is used to ensure lifetime of the Dbo::
  //  Session object (important for writing state or spectrum file to the
  //  database after destruction of the InterSpec object)
  std::shared_ptr<DataBaseUtils::DbSession> sql();
  

  /** Hard re-displays the foreground, background, and 2nd datas.
   Displayed energy range will only be changed if the currently displayed range goes past what the data now displays.
   Will keep current display scale factors for foreground/background
   
   Useful after calibrations.
   */
  void refreshDisplayedCharts();
  
protected:

  //doFinishupSetSpectrumWork(...): intended to do things like calculate the
  //  continuum in the background, or guess and load the detector, from a
  //  background thread (in the WServer thread pool).  The reason this function
  //  is necessary is so that after doing these standard operations, the
  //  modified and modifiedSinceDecode statuses can be reset since we shouldnt
  //  consider these standard things modifications.
  //  Note that the WApplication::UpdateLock is not taken while the workers are
  //  doing their thing, so you should explicitly do this if necassarry, or call
  //  this function from within the main loop.
  void doFinishupSetSpectrumWork( std::shared_ptr<SpecMeas> meas,
                            std::vector<boost::function<void(void)> > workers );
  
  void createOneOverR2Calculator();
  void createUnitsConverterTool();
  void createFluxTool();
  void createDecayInfoWindow();
  void deleteDecayInfoWindow();
  void createFileParameterWindow();
  
  void updateGuiForPrimarySpecChange( std::set<int> display_sample_nums );
  
  static std::set<int> sampleRangeToSet( int start_sample,  int end_sample,
                                         std::shared_ptr<const SpecMeas> meas,
                                         const std::set<int> &excluded_samples );
  
  //sample_real_time_increment(): returns the real time for a given sample
  //  number; this is what is on the x-axis of the time series plot.
  //  Returns the maximum real time of any measurment in the sample, for any of
  //  the specified detecor numbers.
  static float sample_real_time_increment( const std::shared_ptr<const SpecMeas> &meas,
                                           const int sample,
                                           const std::vector<std::string> &detector_names );
  
  /** Function to call when the time chart is clicked or tapped on.
   @param sample_start The starting sample number (as defined by the SpecFile) that was drug.
   @param sample_end The ending sample number (as defined by the SpecFile) that was drug.
   @param modifiers The keyoard modifier pressed.  Ex.
          Shift key means: add specified samples to displayed samples
          Alt (mac option) key means operation being performed is for the background spectrum
          Conrol means remove the specified samples from being displayed
   */
  void timeChartDragged( const int sample_start, const int sample_end,
                             Wt::WFlags<Wt::KeyboardModifier> modifiers );
  
  /** A simple passthrough to #timeChartDragged to handle when time chart is clicked on or tapped.
   */
  void timeChartClicked( const int sample_number, Wt::WFlags<Wt::KeyboardModifier> modifiers );
    
  void displayForegroundData( const bool keep_current_energy_range );
  void displaySecondForegroundData();
  void displayBackgroundData();
  void displayTimeSeriesData();

  
  //detectorsToDisplayChanged(): a callback function for when the user selects a
  //  detector to be displayed or not.
  void detectorsToDisplayChanged();
  
  
  /** Returns a vector of pairs that indicate the cumulative chart starting time of each interval, and
    the sample number of that interval, for use in the time series chart.
    The cumulative times are strictly increasing in value, and may not correspond to the integrated
    real-times of all previous samples, because backgrounds are placed at negative times, and also
    may be compressed if they are really long.
    The last entry contains the upper edge of the last time segment, and a garbage sample number.
    Background samples will always come first, followed by foreground, then "derived" data.  Derived
    data are not sorted according to type (foreground/background/known/intrinsic)
  
   TODO: Make it so this function returns a tuple with:
         {cumulative time, sample number, gamma counts, neutron counts, time, back/fore/second, occupancy}
         and have the time history chart store this info, so this function only ever gets called
         when file is loaded or number detectors changed.
   TODO: Integrate the logic of validForegroundSamples() into this function, and maybe eliminate
         other places
   */
  std::vector<std::pair<float,int> > passthroughTimeToSampleNumber() const;

  //showNewWelcomeDialog(): see notes for showWelcomeDialog().  This function
  //  will eventually replace showWelcomeDialog().
  void showNewWelcomeDialog( bool force = false );

  // Cookie management~  /*should be placed into user options*.
  void setShowIEWarningDialogCookie( bool show );

  /** Adds menu items to "Tools" menu for the tools that are ordinarily shown
      as tabs.  This function is called when the "Hide Tool Tabs" menu option
      is triggered (either manually or by user preference when InterSpec is
      started; also, always called when started on a phone).
   
   \sa removeToolsTabToMenuItems
   \sa m_tabToolsMenuItems
   */
  void addToolsTabToMenuItems();
  
  /** Remove the tool tab items from the "Tools" menu.  This is called when you
      are showing the tool tabs.
   
   \sa addToolsTabToMenuItems
   \sa m_tabToolsMenuItems
   */
  void removeToolsTabToMenuItems();
  
  void showGammaCountDialog();
  void deleteGammaCountDialog();

#if( USE_SPECRUM_FILE_QUERY_WIDGET )
  void showFileQueryDialog();
  void deleteFileQueryDialog();
#endif
  
  void deletePeakEdit();
  void createPeakEdit( double energy );
  void handleRightClick( const double energy, const double counts,
                         const double pageX, const double pageY );
  void handleLeftClick( const double energy, const double counts,
                        const double pageX, const double pageY );
  void rightClickMenuClosed();
  
  void peakEditFromRightClick();
  void refitPeakFromRightClick();
  void deletePeakFromRightClick();
  void addPeakFromRightClick();
  void makePeakFromRightClickHaveOwnContinuum();
  void shareContinuumWithNeighboringPeak( const bool shareWithLeft );
  void handleChangeContinuumTypeFromRightClick( const int peak_continuum_offset_type );
  
  //updateRightClickNuclidesMenu(): meant to be called from within the
  //  application loop (so wApp is valid).  Sets the contents of the
  //  "Change Nuclide" sub-menu of the right click menu, to nuclides
  void updateRightClickNuclidesMenu( const std::shared_ptr<const PeakDef> p,
                      std::shared_ptr< std::vector<std::string> > nuclides );
  //setPeakNuclide(): changes the nuclide of the peak passed in to the nuclide
  //  named; if the nuclide has energy information associated with it, it will
  //  be taken into account (see PeakModel::setNuclideXrayReaction() for details
  //  of how).
  //The details of xrays, vs reactions, vs nuclides hasnt been worked out yet!
  void setPeakNuclide( const std::shared_ptr<const PeakDef> peak,
                       std::string nuclide );
  
  std::shared_ptr<const PeakDef> nearestPeak( const double energy ) const;
  
  void setIsotopeSearchEnergy( double energy );
  
//If we are using D3 to render the spectrum chart, we need to have feature markers available
//  to be able to display/hide them on the chart
public:
  //Tracking of which feature markers are being shown on the c++ side of things
  //  is currently only used for export to the D3 chart...
  
  void setFeatureMarkerOption( const FeatureMarkerType option, const bool show );
  bool showingFeatureMarker( const FeatureMarkerType option );
  void setComptonPeakAngle( const int angle );
  void toggleFeatureMarkerWindow();
  void deleteFeatureMarkerWindow();
  
public:

#if( USE_DB_TO_STORE_SPECTRA )
  //getCurrentStateID() open access to the current state ID
  int getCurrentStateID() { return m_currentStateID; }
#endif
  
  //Peak finding functions
  void searchForSinglePeak( const double x );
  
  /** Function to call when the automated search for peaks (throughout the
      entire spectrum) begins.  Currently all this function does is disable the
      peak search button.
   */
  void automatedPeakSearchStarted();
  
  /** Function to call when the automated search for peak finishes.
     Re-enables peak search button and forces a chart re-drawing if D3 based
     dispay.
   */
  void automatedPeakSearchCompleted();
  
  /** Returns if the current color theme says you should make peaks the same
     color as the reference photopeak lines assigned to the peak.
   */
  bool colorPeaksBasedOnReferenceLines() const;
  
  //searchForHintPeaks(): launches the job to search for peaks (single threaded)
  //  which will call setHintPeaks(...) when done.
  void searchForHintPeaks( const std::shared_ptr<SpecMeas> &data,
                           const std::set<int> &samples );
  
  //setHintPeaks(): sets the hint peaks (SpecMeas::m_autoSearchPeaks and
  //  SpecMeas::m_autoSearchInitialPeaks) if spectrum.lock() yeilds a valid ptr.
  //  If the user has changed peaks from existingPeaks, then results will be
  //  merged.
  //  This function should be called from the main event loop.
  void setHintPeaks( std::weak_ptr<SpecMeas> spectrum,
                     std::set<int> samplenums,
                     std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > existingPeaks,
                     std::shared_ptr<std::vector<std::shared_ptr<const PeakDef> > > resultpeaks );
  
  
  //findPeakFromUserRange(): Depreciated 20150204 by wcjohns in favor of calling
  // InterSpec::findPeakFromControlDrag().  Keeping around JIC for a little
  // while.
  //void findPeakFromUserRange( double x0, double x1 );
  
  
  void excludePeaksFromRange( double x0, double x1 );
  
  //guessIsotopesForPeaks(): makes a best guess for which isotopes are
  //  responsible for the identified photopeaks.  Does not modify peaks which
  //  have already been assigned an isotope.
  //If you wish to perform the work outside of the event loop (e./g. in a
  //  background thread), you should pass in the WApplication pointer, however
  //  the work takes an application lock, so gui will become unresponsive anyway
  //  The InterSpec pointer is actually only necessary for the
  //  DetectorPeakResponse.
  void guessIsotopesForPeaks( Wt::WApplication *app );
  
    
  bool isMobile() const;
  bool isPhone() const;
  bool isTablet() const;
  bool isDesktop() const;
  bool isDedicatedApp() const;
  
  
  //Some functions that effect the display options
  void setLogY( bool logy );
  void setBackgroundSub( bool sub );
  void setVerticalLines( bool show );
  void setHorizantalLines( bool show );
  
  /** A "hard" background subtraction alters the data, subtracting the background counts from foreground counts on a channel by
   channel basis, with the resulting spectrum now having incorrect variances.
   */
  void startHardBackgroundSub();
  void finishHardBackgroundSub( std::shared_ptr<bool> truncate_neg, std::shared_ptr<bool> round_counts );
  
  void setXAxisSlider( bool show );
  void setXAxisCompact( bool compact );
  void setShowYAxisScalers( bool show );
  
  ReferencePhotopeakDisplay *referenceLinesWidget();

#if( defined(WIN32) && BUILD_AS_ELECTRON_APP )
  //When users drag files from Outlook on windows into the app
  //  you can call the following functions
  //  - NOT TESTED for ELECTRON (holdover form old Qt target).
  void dragEventWithFileContentsStarted();
  void dragEventWithFileContentsFinished();
#endif
  

  /** Returns the current color theme.
   Returned ptr may not be in the database yet.
   Will always return a valid pointer, and will not throw.
   Will set m_colorTheme if not previously set.
  */
  std::shared_ptr<const ColorTheme> getColorTheme();

  /** Apply the color theme passed in, to the app; if theme is nullptr then the
      current color theme is retirieved fromm the DB and applied.
      Does not alter the users preference, but does set m_colorTheme.
  */
  void applyColorTheme( std::shared_ptr<const ColorTheme> theme );

  
  /** Signal emitted when the color theme is changed.
   Useful when sub-windows are created to keep them in-sync if the color theme changes
   */
  Wt::Signal< std::shared_ptr<const ColorTheme> > &colorThemeChanged();

  
//Applying the color theme from JS seems to work, but only tested on single
//  browser and operating system, so leaving disabled for the moment, until I
//  at least make sure it doesnt cause problems on older borwsers.
//  Also, currently dont use JS for OSX App build, but there isnt really a
//  reason why not, other than I just havent checked it wont cause issues.
#define APPLY_OS_COLOR_THEME_FROM_JS 1
  
#if( APPLY_OS_COLOR_THEME_FROM_JS && !BUILD_AS_OSX_APP && !IOS && !BUILD_AS_ELECTRON_APP )
  /** Sets up client-side JS to detect the operating system color-theme.  Will
   trigger the signal to be emitted right after initial load.
   */
  void initOsColorThemeChangeDetect();
#endif
  
#if( BUILD_AS_OSX_APP || APPLY_OS_COLOR_THEME_FROM_JS || IOS || BUILD_AS_ELECTRON_APP )
  /** Notification that the operating system changed themes (e.g., macOS went
   into dark mode, or light mode).
   
   If app is currently in either default "Default" or "Dark" then this function
   will switch between these two themes to match the OS.  If in a custom theme
   then no action will be take.
   
   Currently only supports themes 'dark', 'light', 'no-preference', and 'no-support'
   */
  void osThemeChange( std::string name );
#endif
  
#if( USE_DB_TO_STORE_SPECTRA )
  /** Returns the database index of the current UserState.
   If the app state is not based off of one saved to the database, will return -1
   (e.g., user has not loaded current app state form database, and has not
   selected "Store As...").
   */
  int currentAppStateDbId();
  
  /** Disasociates apps current state from any state saved to the database.
   Doesnt remove state from database, or change the spectrum or any thing, just
   disables the tagging or saving over current state in the database.
   */
  void resetCurrentAppStateDbId();
#endif
  
protected:
  
#if( USE_DB_TO_STORE_SPECTRA )
  /** Enables or disables menu items related to saving state. */
  void updateSaveWorkspaceMenu();
#endif
  
  
  /** Handles setting the color theme options to the nuclide reference line
  widget.  If you pass in a nullptr, then current theme will be retrieved from
  database, otherwise will use one passed in.
  */
  void setReferenceLineColors( std::shared_ptr<const ColorTheme> theme );
  
  /** Creates the ColorThemeWidget Window that allows users to alter themes/colors. */
  void showColorThemeWindow();
  
#if( !ANDROID && !IOS )
  void initDragNDrop();
#endif //#if( !ANDROID && !IOS )
  
  void initHotkeySignal();
  void hotKeyPressed( const unsigned int value );
  void arrowKeyPressed( const unsigned int value );
  
  //detectClientDeviceType(): goes through and sets the various bits of
  //  m_clientDeviceType, according to the user agent string, as well as
  //  other environment variables.  Not excedingly well implemented/tested
  //  as of 20140110, but still pretty reasonable.
  void detectClientDeviceType();
  
  
protected:
  PeakModel *m_peakModel;
  D3SpectrumDisplayDiv *m_spectrum;
  D3TimeChart *m_timeSeries;
  
  PopupDivMenu *m_detectorToShowMenu;
  Wt::WPushButton *m_mobileMenuButton;
  Wt::WContainerWidget *m_mobileBackButton;
  Wt::WContainerWidget *m_mobileForwardButton;
  Wt::WContainerWidget *m_notificationDiv; //has id="qtip-growl-container"
  
  void handleUserIncrementSampleNum( SpecUtils::SpectrumType type, bool increment);

  
 /* Start widgets this class keeps track of */

  Wt::Signal< Wt::WString, Wt::WString, int > m_messageLogged;
  
  WarningWidget          *m_warnings;
  AuxWindow              *m_warningsWindow;
  
  SpecMeasManager        *m_fileManager; // The file manager
  
  
#if( USE_CSS_FLEX_LAYOUT )
  Wt::WContainerWidget *m_chartResizer;
  Wt::WContainerWidget *m_toolsResizer;
#else
  Wt::WGridLayout        *m_layout;
  
  Wt::WContainerWidget   *m_charts;
  Wt::WContainerWidget   *m_chartResizer;
  Wt::WGridLayout        *m_toolsLayout;
#endif
  
  Wt::WContainerWidget   *m_menuDiv; // The top menu bar.

  //m_peakInfoWindow is deleted when tool tabs are shown because the layout
  //  gets all messed up for some reason when m_peakInfoDisplay is removed
  //  and placed back in it.
  PeakInfoDisplay        *m_peakInfoDisplay;
  AuxWindow              *m_peakInfoWindow;

  //m_peakEditWindow: used to ensure only one peak editor window is open at a
  //  time.  Will be null if no peak editor is open; valid if one is open.
  PeakEditWindow *m_peakEditWindow;
  
  
  //m_currentToolsTab: used to track which tab is currently showing when the
  //  tools tab is shown.  This variable is necessary so that handleToolTabChanged(...)
  //  can know what tab is being changed from, and not just what tab is being
  //  changed to.  In everyother function this variable will always be up to
  //  date when in tool tabs are shown.
  int m_currentToolsTab;
  
  //m_toolsTabs: will be null when not tool tabs are hidden, and non-null in
  //  when they are visible
  Wt::WTabWidget *m_toolsTabs;

  EnergyCalTool          *m_energyCalTool;
  AuxWindow              *m_energyCalWindow;
  GammaCountDialog       *m_gammaCountDialog;
  AuxWindow              *m_specFileQueryDialog;

  Wt::WSuggestionPopup   *m_shieldingSuggestion;
  ShieldingSourceDisplay *m_shieldingSourceFit;
  AuxWindow              *m_shieldingSourceFitWindow;
  std::shared_ptr<MaterialDB> m_materialDB;

  //m_nuclideSearchWindow: only valid when in tool tabs are hidden, and the user
  //  currently has the window nuclide search window open.
  AuxWindow              *m_nuclideSearchWindow;
  
  //m_nuclideSearchContainer: holds the nuclide search content in tab when tool
  //  tabs are visible.  Will be valid when in tool tabs visible, and null when
  //  not.
  WContainerWidget       *m_nuclideSearchContainer;
  
  //m_nuclideSearch: Nuclide Search widget.  Will always be a valid pointer,
  //  although not always in the DOM (specifically when tool tabs are hidden).
  IsotopeSearchByEnergy  *m_nuclideSearch;
  
  //DataBaseUtils::DbSession is an indirect way to holds the Wt::Dbo::Session
  //  object associated with m_user.  This indirection forces you to use
  //  DataBaseUtils::DbTransaction to interact with the database, which is all
  //  a thread safe way to use the Wt::Dbo::Session associated with m_user (so
  //  we can post jobs to the thread pool that require writing to the database).
  //  A shared pointer is used to ensure lifetime of the Dbo::Session object
  //  when writing to the database, so it wont get deleted when the
  //  InterSpec class is deleted.
  std::shared_ptr<DataBaseUtils::DbSession> m_sql;
  
  //This menu implementation uses somethng that visually looks like a WPopupMenuItem.
  PopupDivMenu         *m_fileMenuPopup;
  PopupDivMenu         *m_toolsMenuPopup;
  PopupDivMenu         *m_helpMenuPopup;
  PopupDivMenu         *m_displayOptionsPopupDiv;
  
#if( USE_DB_TO_STORE_SPECTRA )
  PopupDivMenuItem *m_saveState;
  PopupDivMenuItem *m_saveStateAs;
  PopupDivMenuItem *m_createTag;
  
  //m_currentStateID: if the user has saved a state, this will be the database
  //  ID of that state.  If there is no saved state, this variable will be -1.
  //  It will also be reset to -1 when a new foreground is loaded.
  int m_currentStateID;
#endif
  
  enum RightClickItems
  {
    kPeakEdit,
    kRefitPeak,
    kRefitROI,
    kChangeNuclide,
    kChangeContinuum,
    kDeletePeak,
    kAddPeak,
    kShareContinuumWithLeftPeak,
    kShareContinuumWithRightPeak,
    kMakeOwnContinuum,
    kNumRightClickItems
  };//enum RightClickItems
  
  PopupDivMenu         *m_rightClickMenu;
  double                m_rightClickEnergy;
  Wt::WMenuItem        *m_rightClickMenutItems[kNumRightClickItems];
  PopupDivMenu         *m_rightClickNuclideSuggestMenu;
  PopupDivMenu         *m_rightClickChangeContinuumMenu;
  
  Wt::WMenuItem        *m_showPeakManager;
  
#if( !IOS && !ANDROID )
#define USE_SAVEAS_FROM_MENU 1
  //m_downloadMenu:
  PopupDivMenu *m_downloadMenu;
  //m_downloadMenus: menu items that allow user to download current showing
  //  spectrums (SpecUtils::SpectrumType::Foreground, SpecUtils::SpectrumType::SecondForeground, SpecUtils::SpectrumType::Background)
  PopupDivMenu *m_downloadMenus[3];
#else
#define USE_SAVEAS_FROM_MENU 0
#endif
  
  //If I ever get the preference tracking stuff working better, I could probably
  //  eliminate the following variables
  PopupDivMenuItem *m_logYItems[2];
  PopupDivMenuItem *m_toolTabsVisibleItems[2];
  PopupDivMenuItem *m_backgroundSubItems[2];
  PopupDivMenuItem *m_hardBackgroundSub;
  PopupDivMenuItem *m_verticalLinesItems[2];
  PopupDivMenuItem *m_horizantalLinesItems[2];
  PopupDivMenuItem *m_showXAxisSliderItems[2];
  PopupDivMenuItem *m_showYAxisScalerItems[2];
  PopupDivMenuItem *m_compactXAxisItems[2];
  
  enum class ToolTabMenuItems
  {
    EnergyCal,
    RefPhotopeaks,
    PeakManager,
    NuclideSearch,
    AutoPeakSearch,
    Seperator,
    NumItems
  };//enum ToolTabMenuItems
  
  /** When the tool tabs are hidden, these menu items will show the respective
   tool that would normally be in a tab.  These menu items are deleted when tool
   tabs are shown.
   
   \sa addToolsTabToMenuItems
   \sa removeToolsTabToMenuItems
   */
  Wt::WMenuItem *m_tabToolsMenuItems[static_cast<int>(ToolTabMenuItems::NumItems)];
  
  /** Tracks which features are being used
   TODO: remove this member variable and, and use m_featureMarkers to track things
   */
  bool m_featureMarkersShown[static_cast<int>(FeatureMarkerType::NumFeatureMarkers)];
  
  /** A window that controls if S.E., D.E., Compton Peak, Compton Edge, or Sum
   Peaks are shown.  Is null when window is not showing.
   */
  FeatureMarkerWindow *m_featureMarkers;
  
  PopupDivMenuItem *m_featureMarkerMenuItem;

  
#if( USE_GOOGLE_MAP )
  PopupDivMenuItem *m_mapMenuItem;
#endif
  
#if( USE_SEARCH_MODE_3D_CHART )
  PopupDivMenuItem *m_searchMode3DChart;
#endif

  PopupDivMenuItem *m_showRiidResults;
  
#if( USE_TERMINAL_WIDGET )
  PopupDivMenuItem *m_terminalMenuItem;
  TerminalWidget   *m_terminal;
  AuxWindow        *m_terminalWindow;
#endif
  
  std::set<int> m_excludedSamples;//these are samples that should not be displayed for the primary spectrum
  std::set<int> m_displayedSamples;
  std::set<int> m_backgroundSampleNumbers;
  std::set<int> m_sectondForgroundSampleNumbers;
  
  enum ClientDeviceType
  {
    DesktopBrowserClient = 0x01,
    PhoneClient          = 0x02,
    TabletClient         = 0x04,
    MobileClient         = 0x08,
    HighBandwithClient   = 0x10,
    DedicatedAppClient   = 0x20,
    NumClientDeviceType  = 0x40
  };//enum ClientDeviceType
  
  //m_targetDevice: a variable to cache the type of client we have, so we dont
  //  have to do string compares (user agent, ip address, etc) each time we want
  //  to know something.  Bits are set according to ClientDeviceType enum.
  unsigned int m_clientDeviceType;

  //m_referencePhotopeakLines: is a pointer to the widget where you can type in
  //  nuclides, reactions or elements (ex "U235", "W", "Ge(n,n)") to see the
  //  reference photpeaks on the energy spectrum chart.
  ReferencePhotopeakDisplay *m_referencePhotopeakLines;
  AuxWindow                 *m_referencePhotopeakLinesWindow;

  /** m_licenseWindow: pointer to window showing disclaimers, licenses, and
      other statements.  Will be nullptr if not showing.
   */
  LicenseAndDisclaimersWindow *m_licenseWindow;
  
  //m_useInfoWindow: a pointer to the window initially shown when app starts up
  //  that allows you to select a previous work state or spectrum, or to watch
  //  the how-to movies.  Will be null if a window is currently not showing,
  //  and non-null if one is showing.  Is used to ensure only one window is open
  //  at a time.
  UseInfoWindow *m_useInfoWindow;
  
  /** Used to make sure only one decay info window is shwoing, and also apply
   color changes when color theme is updated
   */
  DecayWindow *m_decayInfoWindow;
  
  //m_preserveCalibWindow: a pointer to the window that prompts the user if they
  //  would like to use a calibration from a previously used spectrum if the one
  //  they just uploaded is from the same detector as the previous one.
  EnergyCalPreserveWindow *m_preserveCalibWindow;
  
  
  //Current width and height are set in layoutSizeChanged(...).
  int m_renderedWidth;
  int m_renderedHeight;

  //The below pointers may not be valid, so make sure to check them
//  std::shared_ptr<SpecMeas> m_measurments[SpecUtils::SpectrumType::Background+1]; - should convert the below to use the following
  std::shared_ptr<SpecMeas> m_dataMeasurement;
  std::shared_ptr<SpecMeas> m_secondDataMeasurement;
  std::shared_ptr<SpecMeas> m_backgroundMeasurement;  //If a "passthrough" spectrum is uploaded
                                           //  for background, then by default all
                                           //  time slices are displayed

  //TODO: m_displayedSpectrumChangedSignal is kinda heavy duty, I could probably
  //      due with just emmitting SpecUtils::SpectrumType and not sampleNumbers or
  //      measurments
  Wt::Signal<SpecUtils::SpectrumType,
             std::shared_ptr<SpecMeas>/*measurment*/,
             std::set<int>, /*sample_numbers*/
             std::vector<std::string> /*detectors*/
             > m_displayedSpectrumChangedSignal;

  Wt::Signal<SpecUtils::SpectrumType,double> m_spectrumScaleFactorChanged;
  

  //Signals for when the current detector is changed or modified
  Wt::Signal<std::shared_ptr<DetectorPeakResponse> > m_detectorChanged;
  Wt::Signal<std::shared_ptr<DetectorPeakResponse> > m_detectorModified;

  //Connections to the current foreground SpecMeas object that need to be
  //  connected and disconnected when changing spectrums.
#ifdef WT_USE_BOOST_SIGNALS2
  boost::signals2::connection m_detectorChangedConnection;
  boost::signals2::connection m_detectorModifiedConnection;
  boost::signals2::connection m_displayedSpectrumChanged;
#else
  boost::signals::connection m_detectorChangedConnection;
  boost::signals::connection m_detectorModifiedConnection;
  boost::signals::connection m_displayedSpectrumChanged;
#endif

  //m_unNamedJSlots: consists of ptrs to JSlots we want to keep around, but
  //  dont care enough to make them member variables so we can be sure to delete
  //  them eventually.
  std::vector< std::shared_ptr<Wt::JSlot> > m_unNamedJSlots;

  std::unique_ptr< Wt::JSignal<unsigned int> > m_hotkeySignal;
  
#if( APPLY_OS_COLOR_THEME_FROM_JS && !BUILD_AS_OSX_APP && !IOS )
  /** Signal emitted from JS when the operating systems color theme changes, or
      right after application starts.
   
   The string can have values: 'dark', 'light', 'no-preference', and 'no-support'
   */
  std::unique_ptr<Wt::JSignal<std::string> > m_osColorThemeChange;
#endif
  
  /** Determine if peak color should be assigned based on refernce line (if
      being shown) it gets associated with.  Controlled as part of the
      #ColorTheme.
   */
  bool m_colorPeaksBasedOnReferenceLines;
  
  std::string m_currentColorThemeCssFile;
  
  std::shared_ptr<const ColorTheme> m_colorTheme;
  
  Wt::Signal<std::shared_ptr<const ColorTheme>> m_colorThemeChanged;
  
  bool m_findingHintPeaks;
  std::deque<boost::function<void()> > m_hintQueue;
  
  static std::mutex sm_staticDataDirectoryMutex;
  static std::string sm_staticDataDirectory;
  
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER )
  static std::mutex sm_writableDataDirectoryMutex;
  static std::string sm_writableDataDirectory;
#endif  //if( not a webapp )

#if( INCLUDE_ANALYSIS_TEST_SUITE )
  friend class SpectrumViewerTester;
  
  D3SpectrumDisplayDiv *spectrum(){ return m_spectrum; }
#endif
};//class InterSpec

#endif


