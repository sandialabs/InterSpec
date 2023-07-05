#ifndef EnergyCalTool_h
#define EnergyCalTool_h
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
#include <string>
#include <vector>
#include <memory>

#include <Wt/WContainerWidget>

#include "InterSpec/AuxWindow.h"

class SpecMeas;
class AuxWindow;
class InterSpec;
class PeakModel;
class RowStretchTreeView;
class EnergyCalAddActionsWindow;
class EnergyCalGraphicalConfirm;

namespace Wt
{
  class WText;
  class WCheckBox;
  class WGridLayout;
  class WPushButton;
  class WStackedWidget;
}//namespace Wt

namespace SpecUtils
{
  struct EnergyCalibration;
  enum class SpectrumType : int;
}


/** Previous to 20211213 the "Fit Coeffs" button was way over to the left, in the "More Actions" column.  I also crammed the download
 CALp link over there as well.
 
 However, it makes more sense to have this button and link near the actual coefficients, but the problem is that this area is owned by the
 EnergyCalImp::CalDisplay widget - so to have the buttons/links there, it requires some additional boilerplate.  For them moment, we'll
 use this compile-time define to dictate where things are, until we decide on things.
 */
#define IMP_COEF_FIT_BTN_NEAR_COEFS 1
#define IMP_CALp_BTN_NEAR_COEFS 0


namespace EnergyCalImp
{
  class CalDisplay;
  class CALpDownloadResource;
}


/** This class handles energy calibration by the user.
 
 Energy calibration can be a suprising complicated or nuanced operation, but the user almost never
 cares about the details, and they just want to make a peak lineup with a reference line.  So the
 goal of this class is facilitate the most complex energy calibrations, while defaulting to letting
 the user not worry about the details.
 
 There can be a wide variation of energy calibration scenerious in spectrum files:
 - Single record with defined energy calibration
 - Multiple samples, all from the same detector, all with the same energy calibration
 - Multiple samples, all from the same detector, but calibration changes between samples
 - Multiple samples, each with different detectors, but all detectors and samples use same energy
   calbration
 - Multiple samples, each with different detectors, but each detector uses a different energy
   calibration, but that remains same for all samples.
 - Multiple samples, each with different detectors, and some or all detectors calibration may change
   throughout with sample.
 - Some sample/detectors may be lower channel, while others polynomial or FRF
 - Different deviation pairs for each detector, but constant throughout the file
 - Different deviation pairs for each detector, but may change throughout the file
 
 And the user may want to apply calibration differences in a few ways:
 - Apply change in energy calibration to every detector and every sample in file
 - Apply change to only certain detectors, but for every sample in file
 - Apply change to only the sample numbers currently contributing to spectrum, but all detectors
 - Apply change to only the sample numbers and detectors currently contributing to spectrum
 - Apply change to all detectors, but only visible samples
 Where calibration differences means propogating the change in linear term for detector Aa1, to
 Aa2, and not just simply setting Aa1 and Aa2's linear coefficients equal to each other.
 
 Users may want to:
 - re-bin all spectrum to the same energy calibration before/after doing recalibration
 - linearize spectra
 - change from polynomial to FRF, etc
 - Edit deviation pairs
 - Convert lower channel energies to polynomial/FRF coefficients, or other way
 - Revert energy calibration to original calbiration given in file
 - Revert to a previous energy calibration
 
 The ways users may want to perform a calbiration include:
 - Fit one or more peaks that have an assigned nuclide, and use that to fit coefficients
 - Drag spectrum using mouse to where it should be; they may do this for, offset, gain, or deviation
   pairs
 - Use multiple files to fit coefficients (useful for lower resolution detectors that measured one
   nuclide at a time).
 
 
 Current assumptions made by this tool:
 - ...
 
 TODO Tues:
 - Add a "Only Changed Detector option" to "Apply Changes To" colum
 - [implemented but not enabled] By default dont show deviation pairs display unless a detector has
   them defined; maybe a "+ deviation pairs" button or something incase user wants to show them
 - Implement undu...
 */

enum class MoreActionsIndex : int
{
  Linearize,
  Truncate,
  CombineChannels,
  ConvertToFrf,
  ConvertToPoly,
  MultipleFilesCal,
  //Maybe add sum spectra, although this isnt a energy cal related thing
  //Revert, and revert back to original calibrations.
  NumMoreActionsIndex
};//enum MoreActionsIndex

/** A struct that indicates what SpecUtils::Measurement's to apply a coefficent change to.
  \TODO: if (or hopefully when) the InterSpec class allows selecting detectors seperately for
         foreground/back/sec., we will need to consider upgrading how we indicate things \
         because there is an edge-case where detectors wanted will differ by sample number
 
 \sa MeasToApplyCoefChangeToWeak (for undo/redo)
 */
struct MeasToApplyCoefChangeTo
{
  std::shared_ptr<SpecMeas> meas;
  std::set<int> sample_numbers;
  std::set<std::string> detectors;
};//struct MeasToApplyCoefChangeTo



class EnergyCalTool : public Wt::WContainerWidget
{
public:
  EnergyCalTool( InterSpec *viewer, PeakModel *peakModel, Wt::WContainerWidget *parent = nullptr );
  virtual ~EnergyCalTool();
  
  void setWideLayout();
  void setTallLayout();
  
  void refreshGuiFromFiles();
  
  void handleGraphicalRecalRequest( double xstart, double xfinish );
  
  void deleteGraphicalRecalConfirmWindow();
  
  void updateFitButtonStatus();
  
#if( !IMP_CALp_BTN_NEAR_COEFS )
  void updateCALpButtonsStatus();
#endif
  
  void displayedSpecChangedCallback( const SpecUtils::SpectrumType,
                                     const std::shared_ptr<SpecMeas> meas,
                                     const std::set<int> samples,
                                     const std::vector<std::string> detectors );
  
  void userChangedCoefficient( const size_t coefnum, EnergyCalImp::CalDisplay *display );
  
  /** Applies the change in coefficients from 'orig' to 'updated' to all the measurements the
   user has selected in the GUI (e.g., Back/For/Sec, visible detectors, entire file, etc).
   
   Only propagates for the coefficients, not changes to deviation pairs.
   Changes to deviation pairs will cause spectra that dont originally have the 'orig' calibration to
   have their coefficients altered to account for the change in deviation pairs, which is probably
   never the wanted behavior.
   
   Throws std::exception on error.
   */
  void applyCalChange( std::shared_ptr<const SpecUtils::EnergyCalibration> orig,
                       std::shared_ptr<const SpecUtils::EnergyCalibration> updated,
                       const std::vector<MeasToApplyCoefChangeTo> &changemeas,
                       const bool isOffsetOnly );
  
  /** Sets the energy calibration of the specified SpecMea/SampleNumber/Detector to the given energy calibration.
   
   Updates peaks only if `adjust_peaks` is true, and one of the calibrations being set corresponds to a calibration being set.
   
   Note, does not call #EnergyCalTool::doRefreshFromFiles or #InterSpec::refreshDisplayedCharts (which #applyCalChange does),
   so you should call these functions after you call into this function all the times you currently will.
   
   Throws std::exception on error.
   */
  void setEnergyCal( std::shared_ptr<const SpecUtils::EnergyCalibration> new_cal,
                     const MeasToApplyCoefChangeTo &changemeas,
                     const bool adjust_peaks );
  
  /** Adds a deviation pair from the graphical calibration.  Adds the deviation pair measurements
   the user has selected in the GUI (e.g., Back/For/Sec, visible detectors, entire file, etc).
   */
  void addDeviationPair( const std::pair<float,float> &new_pair );
  
  
  /** Callback function that gets called when the user changes the deviation pairs from the gui.
   
   @param display The calibration display that got altered.
   @param devPairFieldType Which part of deviation pair got altered; ss value of enum type
          #DeviationPairDisplay::UserFieldChanged.
   
   */
  void userChangedDeviationPair( EnergyCalImp::CalDisplay *display, const int devPairFieldType );
  
  void displayedSpectrumChanged( const SpecUtils::SpectrumType type,
                                 const std::shared_ptr<SpecMeas> &meas,
                                 const std::set<int> &samples,
                                 const std::vector<std::string> &detectors );
  
  
  void setShowNoCalInfo( const bool nocal );
  
  //setWasGraphicalRecal( int type ): 'type' is of type
  //  GraphicalRecalConfirm::RecalTypes,
  void setWasGraphicalRecal( int type, float energy );
  
  /** Gives a text summary of what spectra changes will be applied to.
   Will be blank if unambigous.
   */
  std::string applyToSummaryTxt() const;
  
  /** Returns which SpecUtils::Measurement need to be updated, based on what files are loaded and
   what options the user has choosen.
   */
  std::vector<MeasToApplyCoefChangeTo> measurementsToApplyCoeffChangeTo();
  
  
  
  /** Applies an energy calibration to the currently displayed spectrum.
   This is the function called when the user drops a CALp file onto InterSpec.
   
   Effects all currently showing detectors, but if there is multiple detectors, and only a single calibration passed in, will take displayed
   spectrum as reference and display differences to the rest of the detectors; handling may change in the future.
   
   Throws exception on error.
   */
  void applyCALpEnergyCal( std::map<std::string,std::shared_ptr<const SpecUtils::EnergyCalibration>> det_to_cal,
                            const SpecUtils::SpectrumType specfile,
                            const bool all_detectors, const bool all_samples );
  
  /** The spectrum type (foreground, background, secondary) of the coefficients currently showing.
   
   i.e., in the "Detector" column, which tab ("For.", "Back", "Sec.") is currently selected.  If there is only a single detector, and only a
   single file, this whole column will likely not be shown, and Foreground will be returned.
   */
  SpecUtils::SpectrumType typeOfCurrentlyShowingCoefficients() const;
  
  /** The detector name currently selected in the currently selected in the "Detector" column.
   If there is only a single detector, then the option for the user to select different detectors wont be shown, and this function will return
   the relevant detectors name.
   
   Currently commented out since it is unused and untested.
   */
  //std::string detectorNameOfCurrentlyShowingCoefficients() const;
  
  /** Returns the pointer to the CALp resource (derived from Wt::WResource) that can be used to export CALp files.
   */
  EnergyCalImp::CALpDownloadResource *calpResources();

  /** Makes a dialog the user can then use to upload a CALp file */
  void handleRequestToUploadCALp();
  
  
  void moreActionBtnClicked( const MoreActionsIndex index );
  
  void cancelMoreActionWindow();
  
protected:
  enum class LayoutType{ Tall, Wide };
  void initWidgets( LayoutType layout );
  
  
  void doRefreshFromFiles();
  
  void specTypeToDisplayForChanged();
  
  void fitCoefficients();
  bool canDoEnergyFit();
  
  /** Have #refreshGuiFromFiles only ever get called from #render to avoid current situation
   of calling refreshGuiFromFiles multiple times for a single render (the looping over files can
   get expensive probably)
  */
  
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags) override;
  
  
  enum ApplyToCbIndex
  {
    ApplyToForeground = 0,
    ApplyToBackground,
    ApplyToSecondary,
    ApplyToDisplayedDetectors,
    ApplyToAllDetectors,
    ApplyToDisplayedSamples,
    ApplyToAllSamples,
    NumApplyToCbIndex
  };//enum ApplyToCbIndex
  
  void applyToCbChanged( const ApplyToCbIndex index );
  
  
  /** Returns the gamma detector names that are available for display, given the displayed samples.
   
   Note: this does not return what is actually displayed, use m_interspec->detectorsToDisplay(type)
         for that.
   
   \TODO: Move this function into a static function that takes a SpecFile as input too, and use
         within the InterSpec class.  Probably also extend to neutron detector names.
   */
  std::set<std::string> gammaDetectorsForDisplayedSamples( const SpecUtils::SpectrumType type );
  
#if( IMP_COEF_FIT_BTN_NEAR_COEFS )
  std::vector<EnergyCalImp::CalDisplay *> calDisplays();
#endif
  
  
  enum EnergyCalToolRenderFlags
  {
    FullGuiUpdate = 0x1
  };//EnergyCalToolRenderFlags
  
  Wt::WFlags<EnergyCalToolRenderFlags> m_renderFlags;
  
  InterSpec *m_interspec;
  PeakModel *m_peakModel;
  
  EnergyCalImp::CALpDownloadResource *m_calpResource;
  
  WContainerWidget *m_tallLayoutContent;
  
  RowStretchTreeView *m_peakTable;  /** \TODO: can maybe remove this member variable */
  
  /** Menu that lets you choose which spectrum (For., Back, Sec.) to show the detectors to select */
  Wt::WMenu *m_specTypeMenu;
  
  /** Stack that holds the menus that let you select which detector to show calibration info for. */
  Wt::WStackedWidget *m_specTypeMenuStack;
  
  /** The menus that let you select which detector to show calibration info for.
   These menues are held in the #m_specTypeMenuStack stack.
   
   These entries will be nullptr if there is not a SpecFile for its respective {for, back, sec}
   */
  Wt::WMenu *m_detectorMenu[3];
  
  /** The stack that holds all the calibration info displays for all the spectrum files and
   detectors (i.e., this stack is shared by all m_detectorMenu[]'s).
   */
  Wt::WStackedWidget *m_calInfoDisplayStack;
  
  Wt::WText *m_noCalTxt;
  Wt::WContainerWidget *m_moreActionsColumn;
  Wt::WContainerWidget *m_applyToColumn;
  Wt::WContainerWidget *m_detColumn;
  Wt::WGridLayout *m_detColLayout;
  Wt::WContainerWidget *m_calColumn;
  Wt::WContainerWidget *m_peakTableColumn;
  
  Wt::WGridLayout *m_layout;
  
  Wt::WCheckBox *m_applyToCbs[ApplyToCbIndex::NumApplyToCbIndex];
  
  Wt::WAnchor *m_moreActions[static_cast<int>(MoreActionsIndex::NumMoreActionsIndex)];
  
#if( !IMP_COEF_FIT_BTN_NEAR_COEFS )
  Wt::WPushButton *m_fitCalBtn;
#endif
  
#if( !IMP_CALp_BTN_NEAR_COEFS )
#if( BUILD_AS_OSX_APP || IOS )
  Wt::WAnchor *m_downloadCALp;
#else
  Wt::WPushButton *m_downloadCALp;
#endif
  
  Wt::WPushButton *m_uploadCALp;
#endif
  
  time_t m_lastGraphicalRecal;  // \TODO: switch this to std::chrono::timepoint
  int m_lastGraphicalRecalType;
  float m_lastGraphicalRecalEnergy;
  EnergyCalGraphicalConfirm *m_graphicalRecal;
  
  std::shared_ptr<SpecMeas> m_currentSpecMeas[3];
  std::set<int> m_currentSampleNumbers[3];
  
  EnergyCalAddActionsWindow *m_addActionWindow;
  
  /*
   We could hook up to SpecMeas::aboutToBeDeleted() to know whene file is about to go-away, but this
   could be when the file is being cached on disk or, or unloaded.  So maybe a new signal is needed
   in SpecMeasManager that more properly emits, and tells you whats going on (also, probably need
   to check that memory is actually reclaimed).
  struct CalibrationInformation
  {
    void reset();
    CalibrationInformation();
    
   //should probably use weak_ptrs below
    std::map<std::shared_ptr<const SpecUtils::EnergyCalibration>,std::vector<std::shared_ptr<const SpecUtils::Measurement>>> m_cal_to_meas;
    
    std::shared_ptr<const SpecUtils::EnergyCalibration> m_display_cal; //needed when shifting peaks back
    
    
    SpecUtils::EnergyCalType type;
    std::vector<float> coefficients;
    std::vector< std::pair<float,float> > deviationpairs;
    
    std::set<int> sample_numbers;
    std::set<int> detectors_numbers;
    
    //could also hold original SpecMeas::SampleNumsToPeakMap ...
  };//struct CalibrationInformation
  */
};//class EnergyCalTool

#endif //EnergyCalTool_h

