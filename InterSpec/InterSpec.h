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
#include <tuple>
#include <memory>
#include <vector>

#include <Wt/Dbo/Dbo>
#include <Wt/WContainerWidget>

//Without including InterSpecUser.h here, we get some weird issues with the
//  DB optimistic versioning...
#include "InterSpec/InterSpecUser.h"

class PeakDef;
class PopupDiv;
class SpecMeas;
class AuxWindow;
class GoogleMap;
class PeakModel;
class MaterialDB;
class D3TimeChart;
class DecayWindow;
struct ColorTheme;
class UserFileInDb;
class PopupDivMenu;
class SimpleDialog;
class EnergyCalTool;
class GammaXsWindow;
class MakeDrfWindow;
class OneOverR2Calc;
class SpectrumChart;
class UseInfoWindow;
class WarningWidget;
class DoseCalcWindow;
class FluxToolWindow;
class PeakEditWindow;
class RefLineDynamic;

#if( USE_LLM_INTERFACE )
class LlmToolGui;
#endif
class WarningMessage;
class DrfSelectWindow;
class PeakInfoDisplay;
class SpecMeasManager;
class UndoRedoManager;
class UserPreferences;
class GammaCountDialog;
class PopupDivMenuItem;
class SpectraFileHeader;
class PopupWarningWidget;
class UnitsConverterTool;
struct ExternalRidResults;
class FeatureMarkerWindow;
class D3SpectrumDisplayDiv;
class DetectionLimitWindow;
class DetectorPeakResponse;
class ExportSpecFileWindow;
class MakeFwhmForDrfWindow;
class IsotopeSearchByEnergy;
class ShieldingSourceDisplay;
class EnergyCalPreserveWindow;
class ReferencePhotopeakDisplay;
class DetectionLimitSimpleWindow;
class SimpleActivityCalcWindow;
class LicenseAndDisclaimersWindow;
namespace HelpSystem{ class HelpWindow; }
namespace D3SpectrumExport{ struct D3SpectrumChartOptions; }

#if( USE_TERMINAL_WIDGET )
class TerminalWidget;
#endif

#if( USE_REMOTE_RID )
class RemoteRid;
#endif

#if( INCLUDE_ANALYSIS_TEST_SUITE )
class SpectrumViewerTester;
#endif

#if( USE_REL_ACT_TOOL )
class RelActAutoGui;
class RelActManualGui;
#endif

#if( USE_LEAFLET_MAP )
class LeafletRadMapWindow;
#endif

namespace SpecUtils{ class SpecFile; }
namespace SpecUtils{ class Measurement; }
namespace SpecUtils{ enum class SpectrumType : int; }
namespace SpecUtils{ enum class DetectorType : int; }
namespace SpecUtils{ enum class OccupancyStatus : int; }
namespace PeakSearchGuiUtils{ enum class RefitPeakType : int; }


namespace DataBaseUtils
{
  class DbSession;
}

namespace Wt
{
  class WText;
  class WMenu;
  class WDialog;
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
  static InterSpec_API void setStaticDataDirectory( const std::string &dir );
 
  /** Directory where files like cross sections, materials, detector response
      function, nuclear decay information, and similar are stored.
   */
  static InterSpec_API std::string staticDataDirectory();
 
  /** Returns if the staticDataDirectory has been explicitly set. */
  static InterSpec_API bool haveSetStaticDataDirectory();
  
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP || BUILD_AS_UNIT_TEST_SUITE )
  /** Sets the directory were we can write write the user preference database
   file (if using sqlite3); also the location will search for extra detector
   response functions.
   
   Will throw exception if not empty string and its an invalid directory.
   */
  static InterSpec_API void setWritableDataDirectory( const std::string &dir );
  
  /** Returns the location you can write files to, such as user preferences.
      
   For desktop and mobile apps, this will be the standard application data directory specified by the OS:
   - macOS usually: /Users/username/Library/Containers/gov.sandia.macOS.InterSpec/Data/Library/Application Support/sandia.InterSpec
   - Win32 usually: C:\Users\username\AppData\Roaming\InterSpec
   For development builds, this may be the "${cwd}/user_data"
   
   Will throw exception if hasnt been set (or set with empty string)
   */
  static InterSpec_API std::string writableDataDirectory();
#endif  //if( not a webapp )
  
  /** Function called by Wt when the client reports a change in widget size. */
  virtual void layoutSizeChanged( int width, int height );

  //findAndSetExcludedSamples(...): excludes background and calibration samples
  //  in m_dataMeasurement - placing excluded sample numbers into
  //  m_excludedSamples.
  void findAndSetExcludedSamples( std::set<int> definetly_keep_samples );

#if( SpecUtils_ENABLE_D3_CHART )
  
  //getD3SpectrumOptions(): Output current chart options to a struct
  //  to provide the same "state" of the current user session to
  //  D3.js HTML files.
  D3SpectrumExport::D3SpectrumChartOptions getD3SpectrumOptions() const;
#endif //#if( SpecUtils_ENABLE_D3_CHART )
  
  
  /** Appropriate to call if the user opens a file by double-clicking it in the finder.
  
   Parses the file, and then calls #userOpenFile, so will behave as that function describes.
  
   */
  bool userOpenFileFromFilesystem( const std::string filepath, std::string displayFileName = "" );
  
  /** Function to call when the user opens a spectrum file using the operating system (e.g.,
   double click on file, open a QR code or URL, etc).
   
   Checks if a background to the current foreground, and if so, prompt the user how
   they would like to open it; if the file probably isnt a background to the current foreground,
   then the files is opened as foreground. 
   */
  void userOpenFile( std::shared_ptr<SpecMeas> meas, std::string displayFileName = "" );
  
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
    CheckForRiidResults = 0x02,
    SkipParseWarnings = 0x04
#if( USE_REMOTE_RID )
    , SkipExternalRid = 0x08
#endif
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
   If a spectrum URL, the "RADDATA://G0/" prefix is required.
   
   An example URL is "interspec://drf/specify?v=1"
   
   The passed in string is assumed to be url-encoded.
   
   Throws std::exception if url cant be used.
   */
  void handleAppUrl( const std::string &url_encoded_url );
  
  /** Creates a window users can copy/paste a URL into; just passes the URL to #handleAppUrl */
  SimpleDialog *makeEnterAppUrlWindow();
  void handleAppUrlClosed();
  
  /** Loads the XML file for current locale, to use for localizing strings.
   
   Simply passes through to the `InterSpecApp` function of the same name
  \sa InterSpecApp::useMessageResourceBundle
   */
  void useMessageResourceBundle( const std::string &name );
  
  //For the 'add*Menu(...)' functions, the menuDiv passed in *must* be a
  //  WContainerWidget or a PopupDivMenu
  void addFileMenu( Wt::WWidget *menuDiv, const bool isAppTitlebar );
  void addEditMenu( Wt::WWidget *menuDiv );
  void addViewMenu( Wt::WWidget *menuDiv );
  void addDetectorMenu( Wt::WWidget *menuDiv );
  void addToolsMenu( Wt::WWidget *menuDiv );
  void addPeakLabelSubMenu( PopupDivMenu *parentWidget );
  void addAboutMenu( Wt::WWidget *menuDiv );

  /** Adds a new peak to the peak model, and returns the models index of the new peak.
   
     Does NOT add a undo/redo action.
   
   @param peak the peak to add
   @param associateShownNucXrayRctn is specified true _and_ the user is showing some
          reference gamma lines, OR `ref_line_name` is not empty, and the peak doesnt already have
          a source assigned to it, then the new peak will be assigned to be from the shown lines if
          and are reasonably close.
   @param ref_line_name optional string to specify the source name to assign to the peak.  if
          specified, this string will be tried as a source, before the reference lines.
   
   @returns the WModelIndex of the added peak if the SpectrumType was a foreground and the peak was added
            (may have failed to be added because it is outside the spectrums energy range, or the peak was not initialized)
            If not for the foreground, or the peak was not added, then the returned index will be invalid.
   */
  Wt::WModelIndex addPeak( PeakDef peak, const bool associateShownNucXrayRctn,
                          const SpecUtils::SpectrumType spec_type,
                          const std::string &ref_line_name = "" );
  
  /** Sets the peaks for the given spectrum.  If foreground, you should consider instead to use the PeakModel.
   
   @param spectrum Which spectrum to set the peaks for.
   @param peaks The peaks to set.  Must not be nullptr, or the currently set deque of peaks.
   
   Will trigger update of displayed spectrum.
   
   Throws exception if spectrum type specified is not displayed, or peaks is nullptr.
   */
  void setPeaks( const SpecUtils::SpectrumType spectrum, std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> peaks );
  
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
  
  
#if( USE_DB_TO_STORE_SPECTRA )
  /** Note: see comments in the top of InterSpecUser.h for an explanation of the use-model for saving to the database. */
  
  /**  Returns the measurement that has been, or will be serialized to the database.
   
   @param update If false, then just the last serialization will be returned and in fact may be null.  If true
          then the measurement will be saved to the database first (unless it hasnt been modified since
          last saving) and then returned
  
   Function may throw exception (e.g. FileToLargeForDbException, Wt::Dbo::Exception, std::exception) if it runs into an
   error (file too large to be in DB, or DB issue), or may return a null pointer (the measurement type isnt loaded.
   */
  Wt::Dbo::ptr<UserFileInDb> measurementFromDb( const SpecUtils::SpectrumType type,
                                               const bool update );
  
  //saveStateToDb( entry ): saves the application state to the pointer passed
  //  in.  Note that pointer passed in must be a valid pointer associated with
  //  this m_user.  If `entry->stateType` is specified for as end of session,
  //  or `entry->snapshotTagParent` is non-null, then things will be copied
  //  to unique entries in the database.
  //  May throw exception.  std::runtime_error exceptions will contain
  //  messages that are reasonably okay to send to users.
  void saveStateToDb( Wt::Dbo::ptr<UserState> entry );
  
  //loadStateFromDb( entry ): sets the applications state to that in 'entry'
  void loadStateFromDb( Wt::Dbo::ptr<UserState> entry );

  
  /** Called from the "Store As..." menu item.
   If state is already saved in the DB, updates it.
   This function should only be called when we are already connected to a state in the DB,
   but if not, will then call `stateSaveAs()` (but again, this shouldn't happen).
   */
  void stateSave();
  /** Called from the "Store As..." menu item - creates a dialog to store current state under new name/desc. */
  void stateSaveAs();
  /** Called from the "Store As..." dialog to save things to DB. */
  void stateSaveAsFinish( Wt::WLineEdit *name, Wt::WTextArea *desc, AuxWindow *window );
  /** Called from the "Tag..." menu item, to create a tag of the current state. */
  void stateSaveTag();
  /** Called from the "Tag..." dialog to actually create the tag in the DB. */
  void stateSaveTagFinish( Wt::WLineEdit *name, AuxWindow *window );
  
  /** Saves the Act/Shield and Rel Eff tool states to the in-memory `SpecMeas` objects, then updates
   the database with either the current app state, or the current SpecMeas object, depending if we are
   connected to a app-state or not.
   If connected to an app state, will create, or replace the states `kUserStateAutoSavedWork` state in DB.
   
   @param doAsync If true, saving to the database will be performed asynchronously. If false, will do all work
   during the call.
   */
  void saveStateAtForegroundChange( const bool doAsync );
  
  /** Removes all previous `kEndOfSessionTemp` sessions for the user from the database, and then
   if the "AutoSaveSpectraToDb" preference is true, will create the new `kEndOfSessionTemp` state
   that will be loaded next time the app is started.
   */
  void saveStateForEndOfSession();
  
  /** Called by #InterSpec::stateSaveTagFinish and #InterSpec::stateSaveAsFinish, to actually create
   the db user state, and then call `saveStateToDb(...)` and update menu items.
   */
  void finishStoreStateInDb( const Wt::WString &name,
                            const Wt::WString &description,
                            Wt::Dbo::ptr<UserState> parent );
  
#if( INCLUDE_ANALYSIS_TEST_SUITE )
  void startStateTester();
  
  //Creates a dialog that allows user to name and describe the state; when
  // user then clicks save, will pass off to `storeTestStateToN42(...)`.
  void startStoreTestState();
  
  //storeTestStateToN42(): stores foreground, background, and shielding/source
  //  model into a since 2012 N42 file.  Does not throw, but will notify the
  //  user via the GUI upon error.
  void storeTestStateToN42( const Wt::WString name, const Wt::WString descr );
  
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
                              const SpecUtils::SpectrumType spectrum_type,
                             const bool addUndoRedoStep );
  
  /** This function is called when the user slides the slider on the spectrum, through the
   D3SpectrumDisplayDiv::yAxisScaled() signal.
   */
  void handleDisplayScaleFactorChangeFromSpectrum( const double sf, const double prevSF,
                                                  const SpecUtils::SpectrumType spec_type );
  
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
   */
  void logMessage( const Wt::WString& message, int priority );
  
  /** In case you need to display a custom Toast item, you can access the WarningWidget
   using this function.
   */
  WarningWidget *warningWidget();
  
  virtual Wt::Signal< Wt::WString, int > &messageLogged();
  
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
  
  /** Signal emited when the hint-peaks (the automatic search peaks) for the spectrum is set.
   This signal will always be called after `displayedSpectrumChanged()` for foreground
   background spectra, if either a fresh search for peaks is done, or previously found peaks are
   re-used.
   */
  Wt::Signal<SpecUtils::SpectrumType> &hintPeaksSet();
  
  /** Signal emitted when new external RID results are recieved.  See #RemoteRid. */
  Wt::Signal<std::shared_ptr<const ExternalRidResults>> &externalRidResultsRecieved();
  
  //addHighlightedEnergyRange(): Adds highlighted range to the energy spectrum.
  //  Returns the ID of the highlight region, which you will need to remove
  //  the highlight region.
  size_t addHighlightedEnergyRange( const float lowerEnergy,
                                   const float upperEnergy,
                                   const Wt::WColor &color );
  //removeHighlightedEnergyRange(): removes the region corresponding to the ID
  //  passed in from the energy chart.  Returns successfulness of the removing
  //  the region.
  bool removeHighlightedEnergyRange( const size_t regionid );
  
  //setDisplayedEnergyRange(): sets the displayed energy range that should be
  //  shown.  Range is not checked for validity. E.g. you should not ask to zoom
  //  in to less than one bin.  You can also specify a lower or higher energy
  //  than data is available for.  Definitely dont input the same energy for both
  //  values, or bad things will happen.
  void setDisplayedEnergyRange( float lowerEnergy, float upperEnergy );
  
  /** Sets the Y-axis range.
   
   Returns set range, and empty string if successful, or otherwise a message explaining why the request couldnt be fully honored
   (e.g., a negative lower counts was specified, but its currently log-scale)
   */
  std::tuple<double,double,Wt::WString> setYAxisRange( double lower_counts, double upper_counts );

  /** When displaying the spectrum in log-y, sets the lower y-axis value to show, if there are channels with zero counts. */
  bool setLogYAxisMin( const double lower_value );
  
  //displayedSpectrumRange(): Grab the ranges in y and y that are currently
  //  displayed to the user.
  void displayedSpectrumRange( double &xmin, double &xmax, double &ymin, double &ymax ) const;
  
  void handleShiftAltDrag( double lowEnergy, double upperEnergy );
  
  void setToolTabsVisible( bool show );
  bool toolTabsVisible() const;
  
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
  /** Function called when user changes the tool tab height.
   We will remember this, so we can set it back
   */
  void toolTabContentHeightChanged( int height );
#endif
  
  /** Makes a MakeDrf Window and returns it, or if one was already present, returns it. */
  MakeDrfWindow *showMakeDrfWindow();
  /** Returns the pointer to current MakeDrf Window.  Will by nullptr if not currently showing */
  MakeDrfWindow *makeDrfWindow();
  void handleCloseMakeDrfWindow( MakeDrfWindow *window );
  
  DrfSelectWindow *showDrfSelectWindow();
  void closeDrfSelectWindow();
  
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
  
  GammaXsWindow *showGammaXsTool();
  void deleteGammaXsTool();
  
  DoseCalcWindow *showDoseTool();
  void deleteDoseCalcTool();
  
  /** If "Activity/Shielding Fit" window is not showing, creates the tool/window, and returns the tool.
   If it is already showing, no changes are made, and a pointer to the tool is returned.
   When the window shown to the user is closed, the #closeShieldingSourceFit function will automatically be called, and tool deleted.
   */
  ShieldingSourceDisplay *shieldingSourceFit();
  void closeShieldingSourceFit();
  void saveShieldingSourceModelToForegroundSpecMeas();
  
  /** Shows the help window; if a preselect (ex., "energy-calibration", or "getting-started"), then that topic will be set to be selected.
   */
  void showHelpWindow( const std::string &preselect = "" );

  /** Closes the help window. */
  void closeHelpWindow();
  
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
  
  /** Returns `m_materialDB`.  Note `MaterialDB` should be thread-safe, so you can use this for long running jobs where
   this instance of InterSpec may disappear.
   */
  std::shared_ptr<MaterialDB> materialDataBaseShared();
  
  /** The suggestion pop-up widget for shielding names; used globally for all shielding name inputs
   so that there is not duplicate copies of the widget in the DOM.
   
   Object will be alive as long as *this.
   */
  Wt::WSuggestionPopup *shieldingSuggester();
  
  /** The RefLineDynamic class. */
  RefLineDynamic *refLineDynamic();
  
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
  
  UndoRedoManager *undoRedoManager();
  
  PeakEditWindow *peakEdit();
  
  void showWarningsWindow();
  void handleWarningsWindowClose();

  void showPeakInfoWindow();
  void handlePeakInfoClose();

  GammaCountDialog *showGammaCountDialog();
  void deleteGammaCountDialog();
  
#if( USE_GOOGLE_MAP )
  void createMapWindow( SpecUtils::SpectrumType spectrum_type );
  void displayOnlySamplesWithinView( GoogleMap *map,
                                     const SpecUtils::SpectrumType targetSamples,
                                     const SpecUtils::SpectrumType fromSamples );
#elif( USE_LEAFLET_MAP )
  /** Create a Leaflet map window - first showing a warning to users about fetching data from arbgis, if they havent said
   to not show that warning.
   @param spectrum_type Which spectrum file to display GPS coordinates for.
   @param noWarning If true, the warning dialog will not be shown - should be true only for executing undo/redo.
   */
  void createMapWindow( const SpecUtils::SpectrumType spectrum_type, const bool noWarning );
  
  /** Handles the closing of the Leaflet map warning window (sets `m_leafletWarning` to nullptr). */
  void handleLeafletWarningWindowClose();
  
  
  void handleLeafletMapOpen( LeafletRadMapWindow *window );
  
  /** Handles the closing of the Leaflet map window (when the AuxWIndow gets hidden), and inserts a undo/redo step. */
  void handleLeafletMapClosed();
  
  /** Closes the Leaflet map warning window, without creating an undo/redo step. */
  void programmaticallyCloseLeafletWarning();
  
  /** Closes the Leaflet map window, without inserting undo/redo step. */
  void programmaticallyCloseLeafletMap();
#endif
  
#if( USE_SEARCH_MODE_3D_CHART )
  void create3DSearchModeChart();
  void programmaticallyClose3DSearchModeChart();
  void handle3DSearchModeChartClose( AuxWindow *window );
#endif

  /** Show the RIID results included in the spectrum file; and sets `m_riidDisplay` to window pointer */
  void showRiidResults( const SpecUtils::SpectrumType type );
  
  /** Handles the user closing the RIID Display; sets `m_riidDisplay` to nullptr, and tries adds undo/redo step. . */
  void handleRiidResultsClose();
  
  /** Closes the RIID Display, without setting undo/redo step. */
  void programmaticallyCloseRiidResults();
  
  /** Show images included in the spectrum file. */
  void showMultimedia( const SpecUtils::SpectrumType type );
  
  void handleMultimediaClose( SimpleDialog *dialog );
  
  void programmaticallyCloseMultimediaWindow();
  
#if( USE_TERMINAL_WIDGET )
  void createTerminalWidget();
  void handleTerminalWindowClose();
#endif
  

#if( USE_REMOTE_RID )
  RemoteRid *remoteRid();
  void createRemoteRidWindow();
  void deleteRemoteRidWindow();
  
  void setAutoRemoteRidResultDialog( SimpleDialog *dialog );
  void handleAutoRemoteRidResultDialogClose();
  void programaticallyCloseAutoRemoteRidResultDialog();
#endif


#if( USE_REL_ACT_TOOL )
  RelActAutoGui *showRelActAutoWindow();
  void handleRelActAutoClose();
  
  RelActManualGui *createRelActManualWidget();
  void handleRelActManualClose();
  
  void saveRelActManualStateToForegroundSpecMeas();
  void saveRelActAutoStateToForegroundSpecMeas();
#endif
  
#if( USE_TERMINAL_WIDGET || USE_REL_ACT_TOOL )
  void handleToolTabClosed( const int tabnum );
#endif
  
  
  OneOverR2Calc *createOneOverR2Calculator();
  void deleteOneOverR2Calc();
  
  UnitsConverterTool *createUnitsConverterTool();
  void deleteUnitsConverterTool();
  
  FluxToolWindow *createFluxTool();
  void deleteFluxTool();
  
  DecayWindow *createDecayInfoWindow();
  void deleteDecayInfoWindow();

  /** If no `MakeFwhmForDrfWindow` is currently showing, will create one and return a pointer to it.
   If a window is currently showing, will return a pointer to it.
   @param use_auto_fit_peaks_too Whether to use auto-search for peaks and include those in the table of
          peaks that can be used.  Only has an effect if the window is created; if a window already exists, doesnt
          have an effect.
   */
  MakeFwhmForDrfWindow *fwhmFromForegroundWindow( const bool use_auto_fit_peaks_too );
  
  /** If a `MakeFwhmForDrfWindow` is showing, deletes it, and sets `m_addFwhmTool` to nullptr.
   */
  void deleteFwhmFromForegroundWindow();

#if( USE_LLM_INTERFACE )
  /** Create and show the LLM tool widget in the tools tab. */
  void createLlmTool();
  
  /** Handle cleanup when LLM tool is closed. */
  void handleLlmToolClose();
#endif
  
  /** Will show the disclaimer, license, and statment window, setting
      m_licenseWindow pointer with its value.
   */
  void showLicenseAndDisclaimersWindow();
  
  /** Deletes disclaimer, license, and statment window and sets m_licenseWindow
   to nullptr;
   */
  void deleteLicenseAndDisclaimersWindow();
  
  /** Shows the welcome dialog.
   
   @param force If true, the window is shown immediately, and regardless of the "ShowSplashScreen" user preference, and also
          an undo/redo step will be attempted to be added.  If false then window will only be shown if "ShowSplashScreen"
          user preference is true, and wont be shown immediately, but rather posted to be shown next event loop, and with a
          "do not show again" checkbox.
   */
  void showWelcomeDialog( const bool force );
  void showWelcomeDialogWorker( const bool force );

  /** Deletes the currently showing dialog (if its showing), and sets m_useInfoWindow to null.
   
   @param addUndoRedoStep If true, will attempt to add an undo step to re-open the Welcome dialog
   */
  void deleteWelcomeDialog( const bool addUndoRedoStep );
  
  //deleteEnergyCalPreserveWindow(): deletes the currently showing dialog (if its
  //  showing), and sets m_preserveCalibWindow to null
  void deleteEnergyCalPreserveWindow();
  
  /** Shows the spectrum file export dialog.
   
   Sets the `m_exportSpecFileWindow` member variable to the created window.
   */
  ExportSpecFileWindow *createExportSpectrumFileDialog();
  
  /** Deletes the spectrum file export dialog.
   
   Sets the `m_exportSpecFileWindow` member variable to nullptr.
   */
  void handleExportSpectrumFileDialogClose();

#if( USE_DETECTION_LIMIT_TOOL )
  /** If `query_str` is not empty, the handle app URI function will be called. */
  void showDetectionLimitTool( const std::string &query_str );
  DetectionLimitWindow *createDetectionLimitTool();
  void handleDetectionLimitWindowClose();
  void programmaticallyCloseDetectionLimit();
  
  DetectionLimitSimpleWindow *showSimpleMdaWindow();
  void handleSimpleMdaWindowClose();
  void programmaticallyCloseSimpleMda();
#endif //USE_DETECTION_LIMIT_TOOL
  
  SimpleActivityCalcWindow *showSimpleActivityCalcWindow();
  void handleSimpleActivityCalcWindowClose();
  void programmaticallyCloseSimpleActivityCalc();
  void startSimpleActivityCalcFromRightClick();
  
  /** Brings up a dialog asking the user to confirm starting a new session, and if they select so, will start new session. */
  void startClearSession();
  
  /** Pointer to class to access user preferences. */
  UserPreferences *preferences();
  
  /** The user information in the database. */
  const Wt::Dbo::ptr<InterSpecUser> &user();
  
  /** Calls the `reread()` function on `m_user`, which refreshes the information pointed to by
   `m_user` to match what is currently in the database.  Please note, this may (and maybe always)
   cause the `m_user` pointer to point to a different location in memory.
   */
  void reReadUserInfoFromDb();
  
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
  
  
  /** Enables the "View" -> "Show Images" menu item, based on if the foreground
   spectrum file contains any images.
   */
  void checkEnableViewImageMenuItem();
  
  void createFileParameterWindow( const SpecUtils::SpectrumType type );
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
  
#if( USE_DETECTION_LIMIT_TOOL )
  void fitNewPeakNotInRoiFromRightClick();
  void startAddPeakFromRightClick();
  void searchOnEnergyFromRightClick();
  void startSimpleMdaFromRightClick();
#endif //USE_DETECTION_LIMIT_TOOL
  
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
  

#if( USE_SPECRUM_FILE_QUERY_WIDGET )
  void showFileQueryDialog();
  void deleteFileQueryDialog();
#endif
  
  void deletePeakEdit();
  void createPeakEdit( double energy );
  
  
  void handleRightClick( const double energy, const double counts,
                         const double pageX, const double pageY,
                        const std::string &ref_line_name  );
  void handleLeftClick( const double energy, const double counts,
                        const double pageX, const double pageY,
                       const std::string &ref_line_name );
  void rightClickMenuClosed();
  
  void peakEditFromRightClick();
  
  void refitPeakFromRightClick( const PeakSearchGuiUtils::RefitPeakType type );
  
  void setMeanToRefPhotopeak();
  void deletePeakFromRightClick();
  void addPeakFromRightClick();
  void makePeakFromRightClickHaveOwnContinuum();
  void shareContinuumWithNeighboringPeak( const bool shareWithLeft );
  void handleChangeContinuumTypeFromRightClick( const int peak_continuum_offset_type );
  void handleChangeSkewTypeFromRightClick( const int peak_continuum_offset_type );
  
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
  
  void setIsotopeSearchEnergy( double energy );
  
//If we are using D3 to render the spectrum chart, we need to have feature markers available
//  to be able to display/hide them on the chart
public:

  std::shared_ptr<const PeakDef> nearestPeak( const double energy ) const;

  //Tracking of which feature markers are being shown on the c++ side of things
  //  is currently only used for export to the D3 chart...
  
  void setFeatureMarkerOption( const FeatureMarkerType option, const bool show );
  bool showingFeatureMarker( const FeatureMarkerType option );
  void setComptonPeakAngle( const int angle );
  void toggleFeatureMarkerWindow();
  void displayFeatureMarkerWindow( const bool show );
  
public:

  
  /** Function to add undo/redo step for when the user changes x-axis range of spectrum. */
  void handleSpectrumChartXRangeChange( const double xmin, const double xmax,
                                       const double oldXmin, const double oldXmax,
                                       const bool user_interaction);
  
  //Peak finding functions
  void searchForSinglePeak( const double x, const std::string &ref_line_name, Wt::WFlags<Wt::KeyboardModifier> mods );
  
  
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
                           const std::set<int> &samples,
                          const bool isHPGe );
  
  //setHintPeaks(): sets the hint peaks (SpecMeas::m_autoSearchPeaks and
  //  SpecMeas::m_autoSearchInitialPeaks) if spectrum.lock() yeilds a valid ptr.
  //  If the user has changed peaks from existingPeaks, then results will be
  //  merged.
  //  This function should be called from the main event loop.
  void setHintPeaks( std::weak_ptr<SpecMeas> spectrum,
                     std::set<int> samplenums,
                     std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > existingPeaks,
                     std::shared_ptr<std::vector<std::shared_ptr<const PeakDef> > > resultpeaks );
  
  
  void excludePeaksFromRange( double x0, double x1 );
  
  /** Returns if detected this is a mobile device, based on user-agent string, or compile-time
   options.
   
   The function will return false on tablets when the "TabletUseDesktopMenus" preference is true.
   If you want to know the true detected value, call InterSpecApp::isMobile(), which does not
   account for the user preference.
   */
  bool isMobile() const;
  
  /** Returns if detected this is a mobile device, based on user-agent string. */
  bool isPhone() const;
  
  /** Returns if detected this is a tablet, based on user-agent string.
   
   The function will return false on tablets when the "TabletUseDesktopMenus" preference is true.
   If you want to know the true detected value, call InterSpecApp::isMobile(), which does not
   account for the user preference.
   */
  bool isTablet() const;
  
  /** Returns true id this is not detected to be a mobile device. */
  bool isDesktop() const;
  
  /** Returns !BUILD_FOR_WEB_DEPLOYMENT */
  bool isDedicatedApp() const;
  
  /** Returns true if built for Android, or "Android" found in user-agent string. */
  bool isAndroid() const;
  
  
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
  
  void setXAxisSlider( const bool show, const bool addUndoRedo );
  void setXAxisCompact( bool compact );
  void setShowYAxisScalers( bool show );
  
  ReferencePhotopeakDisplay *referenceLinesWidget();
  
  IsotopeSearchByEnergy *nuclideSearch();
  
  PeakInfoDisplay *peakInfoDisplay();

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
  
#if( APPLY_OS_COLOR_THEME_FROM_JS && !BUILD_AS_OSX_APP && !BUILD_AS_ELECTRON_APP )
  /** Sets up client-side JS to detect the operating system color-theme.  Will
   trigger the signal to be emitted right after initial load.
   */
  void initOsColorThemeChangeDetect();
#endif
  
#if( BUILD_AS_OSX_APP || APPLY_OS_COLOR_THEME_FROM_JS || IOS || BUILD_AS_ELECTRON_APP || BUILD_AS_WX_WIDGETS_APP  )
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
  long long int currentAppStateDbId();
  
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
  
  void initDragNDrop();
  
  void hotKeyPressed( const unsigned int value );
  void arrowKeyPressed( const unsigned int value );
  
  //detectClientDeviceType(): goes through and sets the various bits of
  //  m_clientDeviceType, according to the user agent string, as well as
  //  other environment variables.  Not excedingly well implemented/tested
  //  as of 20140110, but still pretty reasonable.
  void detectClientDeviceType();
  
  /** Changes the local, to the language specified. 
   
   @param languageCode The code for the language, ex. "nl" for Dutch, "fr" for French, "en" for English, or "en_GB" for Great Britain.
   */
  void changeLocale( std::string languageCode );
  
protected:
  Wt::Dbo::ptr<InterSpecUser> m_user;
  UserPreferences *m_preferences;
  
  PeakModel *m_peakModel;
  D3SpectrumDisplayDiv *m_spectrum;
  D3TimeChart *m_timeSeries;
  
  PopupDivMenu *m_detectorToShowMenu;
  Wt::WPushButton *m_mobileMenuButton;
  Wt::WContainerWidget *m_mobileBackButton;
  Wt::WContainerWidget *m_mobileForwardButton;
  Wt::WContainerWidget *m_notificationDiv; //has id="qtip-growl-container"
  
  void handleUserIncrementSampleNum( SpecUtils::SpectrumType type, bool increment);

  Wt::Signal< Wt::WString, int > m_messageLogged;
  
  WarningWidget          *m_warnings;
  AuxWindow              *m_warningsWindow;
  
  SpecMeasManager        *m_fileManager; // The file manager
  
  Wt::WGridLayout        *m_layout;
  
  Wt::WContainerWidget   *m_charts;
  Wt::WContainerWidget   *m_chartResizer;
  Wt::WGridLayout        *m_toolsLayout;
  
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
  //  changed to.  In every other function this variable will always be up to
  //  date when in tool tabs are shown.
  int m_currentToolsTab;
  
  //m_toolsTabs: will be null when not tool tabs are hidden, and non-null in
  //  when they are visible
  Wt::WTabWidget *m_toolsTabs;
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
  int m_toolsTabsContentHeight;
#endif

  EnergyCalTool          *m_energyCalTool;
  AuxWindow              *m_energyCalWindow;
  GammaCountDialog       *m_gammaCountDialog;
  AuxWindow              *m_specFileQueryDialog;

  Wt::WSuggestionPopup   *m_shieldingSuggestion;
  ShieldingSourceDisplay *m_shieldingSourceFit;
  AuxWindow              *m_shieldingSourceFitWindow;
  std::shared_ptr<MaterialDB> m_materialDB;
  
#if( USE_REL_ACT_TOOL )
  RelActAutoGui          *m_relActAutoGui;
  AuxWindow              *m_relActAutoWindow;
  PopupDivMenuItem       *m_relActAutoMenuItem;
  
  RelActManualGui        *m_relActManualGui;
  AuxWindow              *m_relActManualWindow;
  PopupDivMenuItem       *m_relActManualMenuItem;
#endif

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
  
  //This menu implementation uses something that visually looks like a WPopupMenuItem.
  PopupDivMenu         *m_fileMenuPopup;
  PopupDivMenu         *m_editMenuPopup;
  PopupDivMenu         *m_toolsMenuPopup;
  PopupDivMenu         *m_helpMenuPopup;
  PopupDivMenu         *m_displayOptionsPopupDiv;
  
#if( USE_DB_TO_STORE_SPECTRA )
  PopupDivMenuItem *m_saveState;
  PopupDivMenuItem *m_saveStateAs;
  PopupDivMenuItem *m_createTag;
#endif
  
  PopupDivMenu *m_languagesSubMenu;
  
  enum RightClickItems
  {
    kPeakEdit,
    kRefitPeakStandard,
    kRefitRoiStandard,
    kRefitRoiAgressive,
    kRefitPeakWithDrfFwhm,
    kSetMeanToNucOrRefLinePhotopeak,
    kChangeNuclide,
    kChangeContinuum,
    kChangeSkew,
    kDeletePeak,
    kAddPeakToRoi,
    kShareContinuumWithLeftPeak,
    kShareContinuumWithRightPeak,
    kMakeOwnContinuum,

#if( USE_DETECTION_LIMIT_TOOL )
    kFitNewPeakNotInRoi,
    kAddPeakNotInRoi,
    kSearchEnergy,
    kSimpleMda,
    kSimpleActivityCalc,
#endif
    
    kNumRightClickItems
  };//enum RightClickItems
  
  PopupDivMenu         *m_rightClickMenu;
  double                m_rightClickEnergy;
  /** The ref-line info from the client-side when there is a right-click.  This may be a displayed reference line, or it could be a kinematic reference line, and is in the
   form like "Th232;S.E. of 2614.5 keV".
   */
  std::string           m_rightClickRefLineHint;
  Wt::WMenuItem        *m_rightClickMenutItems[kNumRightClickItems];
  PopupDivMenu         *m_rightClickNuclideSuggestMenu;
  PopupDivMenu         *m_rightClickChangeContinuumMenu;
  PopupDivMenu         *m_rightClickChangeSkewMenu;
  
  Wt::WMenuItem        *m_showPeakManager;
    
  PopupDivMenuItem *m_exportSpecFileMenu;
  ExportSpecFileWindow *m_exportSpecFileWindow;
  
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
   TODO: remove this member variable and, and use m_featureMarkersWindow to track things
   */
  bool m_featureMarkersShown[static_cast<int>(FeatureMarkerType::NumFeatureMarkers)];
  
  /** A window that controls if S.E., D.E., Compton Peak, Compton Edge, or Sum
   Peaks are shown.  Is null when window is not showing.
   */
  FeatureMarkerWindow *m_featureMarkersWindow;
  
  PopupDivMenuItem *m_featureMarkerMenuItem;
  PopupDivMenuItem *m_dynamicRefLineEnableMenuItem;
  PopupDivMenuItem *m_dynamicRefLineDisableMenuItem;

  SimpleDialog *m_multimedia;

#if( USE_REMOTE_RID )
  /** When the user has selected spectra to be sent off to external RID analysis, and results to
   be displayed in a dialog, instead of a Toast message, this pointer will keep track of the dialog.
   
   \sa getAutoRemoteRidResultDialog
   \sa handleAutoRemoteRidResultDialogClose
   */
  SimpleDialog *m_autoRemoteRidResultDialog;
#endif
  
  GammaXsWindow *m_gammaXsToolWindow;
  DoseCalcWindow *m_doseCalcWindow;
  OneOverR2Calc *m_1overR2Calc;
  UnitsConverterTool *m_unitsConverter;
  FluxToolWindow *m_fluxTool;
  MakeDrfWindow *m_makeDrfTool;
  
  
#if( USE_GOOGLE_MAP || USE_LEAFLET_MAP )
  PopupDivMenuItem *m_mapMenuItem;
  
#if( USE_LEAFLET_MAP )
  SimpleDialog *m_leafletWarning;
  LeafletRadMapWindow *m_leafletWindow;
#endif
#endif
  
  SimpleDialog *m_enterUri;
  
#if( USE_SEARCH_MODE_3D_CHART )
  PopupDivMenuItem *m_searchMode3DChart;
#endif

  PopupDivMenuItem *m_showRiidResults;
  PopupDivMenuItem *m_showMultimedia;
  
#if( USE_TERMINAL_WIDGET )
  PopupDivMenuItem *m_terminalMenuItem;
  TerminalWidget   *m_terminal;
  AuxWindow        *m_terminalWindow;
#endif
  
#if( USE_REMOTE_RID )
  PopupDivMenuItem *m_remoteRidMenuItem;
  RemoteRid        *m_remoteRid;
  AuxWindow        *m_remoteRidWindow;
#endif

#if( USE_DETECTION_LIMIT_TOOL )
  DetectionLimitSimpleWindow *m_simpleMdaWindow;
  DetectionLimitWindow *m_detectionLimitWindow;
#endif
  SimpleActivityCalcWindow *m_simpleActivityCalcWindow;
  
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
    DedicatedAppClient   = 0x10,
    AndroidClient        = 0x20,
    NumClientDeviceType  = 0x40
  };//enum ClientDeviceType
  
  //m_targetDevice: a variable to cache the type of client we have, so we dont
  //  have to do string compares (user agent, ip address, etc) each time we want
  //  to know something.  Bits are set according to ClientDeviceType enum.
  unsigned int m_clientDeviceType;

  //m_referencePhotopeakLines: is a pointer to the widget where you can type in
  //  nuclides, reactions or elements (ex "U235", "W", "Ge(n,n)") to see the
  //  reference photopeaks on the energy spectrum chart.
  ReferencePhotopeakDisplay *m_referencePhotopeakLines;
  AuxWindow                 *m_referencePhotopeakLinesWindow;
  RefLineDynamic            *m_refLineDynamic;

  HelpSystem::HelpWindow *m_helpWindow;
  
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
  
  MakeFwhmForDrfWindow *m_addFwhmTool;
  
  //m_preserveCalibWindow: a pointer to the window that prompts the user if they
  //  would like to use a calibration from a previously used spectrum if the one
  //  they just uploaded is from the same detector as the previous one.
  EnergyCalPreserveWindow *m_preserveCalibWindow;

#if( USE_LLM_INTERFACE )
  /** Menu item for opening the LLM tool. */
  PopupDivMenuItem *m_llmToolMenuItem;
  /** LLM tool widget for user interaction. */
  LlmToolGui          *m_llmTool;
#endif
  
#if( USE_SEARCH_MODE_3D_CHART )
  /** Pointer to window showing the Search Mode 3D data view. */
  AuxWindow *m_3dViewWindow;
#endif
  
  /** Pointer to window created by the #showRiidResults function. */
  SimpleDialog *m_riidDisplay;
  
  DrfSelectWindow *m_drfSelectWindow;
  
  UndoRedoManager *m_undo;
  
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
  
#if( APPLY_OS_COLOR_THEME_FROM_JS && !BUILD_AS_OSX_APP )
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
  Wt::Signal<SpecUtils::SpectrumType> m_hintPeaksSet;
  
  Wt::Signal<std::shared_ptr<const ExternalRidResults>> m_externalRidResultsRecieved;
  
  /** Some informational messages should only be shown once, like when you click on the
   energy tab, so we'll keep track of if we have shown a message.
   */
  std::set<std::string> m_infoNotificationsMade;
  
#if( INCLUDE_ANALYSIS_TEST_SUITE )
  friend class SpectrumViewerTester;
  
  D3SpectrumDisplayDiv *spectrum(){ return m_spectrum; }
#endif
};//class InterSpec

#endif


