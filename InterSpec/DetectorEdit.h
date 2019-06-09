#ifndef DetectorEdit_h
#define DetectorEdit_h
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

#include <vector>
#include <string>
#include <fstream>

#include <Wt/WColor>
#include <Wt/WSignal>
#include <Wt/WTextArea>
#include <Wt/WContainerWidget>
#include <Wt/Chart/WCartesianChart>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/DetectorPeakResponse.h"


namespace Wt
{
  class WText;
  class WMenu;
  class WLabel;
  class WCheckBox;
  class WLineEdit;
  class WComboBox;
  class WTabWidget;
  class WPushButton;
  class WFileUpload;
  class WPushButton;
  class WTextArea;
  class WStandardItemModel;
}//namespace Wt

class InterSpec;
class InterSpecUser;
class RelEffDetSelect;
class GadrasDetSelect;
class SpectraFileModel;
class DetectorPeakResponse;

namespace DataBaseUtils
{
  class DbSession;
}

class DetectorDisplay : public Wt::WContainerWidget
{
/*
  This class provides a way to consistently display information about the
  currently defined detector, and also allow the user to change this by clicking
  on this widget.

  TODO as of 20121115:
    -Need to style this widget okay, and indicate the user can click on it to
     edit the detector
*/

public:
  //DetectorDisplay constructor: the InterSpec and SpectraFileModel objects
  //  passed in are assumed to be valid and longer lived than this
  //  DetectorDisplay
  DetectorDisplay( InterSpec *specViewer,
                   SpectraFileModel *fileModel,
                   Wt::WContainerWidget *parent = 0 );

  virtual ~DetectorDisplay();

  void setDetector( std::shared_ptr<DetectorPeakResponse> det );
  void editDetector();

  std::shared_ptr<DetectorPeakResponse> detector();

protected:
  static const char * const sm_noDetectorTxt;
  static const char * const sm_noDetectorTxtMbl;

  Wt::WText *m_text;

  InterSpec *m_interspec;
  SpectraFileModel *m_fileModel;
  std::weak_ptr<DetectorPeakResponse> m_currentDetector;
};//class DetectorDisplay


class DetectorEdit : public Wt::WContainerWidget
{
/*
  This class is meant to allow the user to change or edit the detector
  associated with a file in a consistent way across all windows that might use
  a detector.
*/
  
public:
  /** DetectorEdit Constructor
   @param currentDet Current DRF; will be used to help select what DRF to show.
   @param parentWindow Used only to place the close, cancel and help buttons
          into the footer
   */
  DetectorEdit( std::shared_ptr<DetectorPeakResponse> currentDet,
                InterSpec *specViewer,
                SpectraFileModel *fileModel,
                AuxWindow *parentWindow = 0
 );
  virtual ~DetectorEdit();

  
  //verify manual definition, set Apply to enabled if ok
  void verifyManualDefinition();
  
  //User defined detector in functional form
  void setDefineDetector();
  
  //Action when the user clicks on the ButtonGroup to select
  //defined/upload/functional form definition for the detector
  void selectButton( Wt::WStackedWidget * stack,
                     Wt::WMenu *group,
                     bool activateCallBack = true);
    
  //done(): signal emmitted when the user is done using this widget.  By the
  //  time this signal is emmtitted, all detectorChanged or detectorModified
  //  signals that need to be emitted, have been.
  Wt::Signal<> &done();

  //emitChangedSignal(): causes the detectorChanged signal to be emitted only if
  //  it is necessary (eg different detector object than last time an emit was
  //  done)
  void emitChangedSignal();

  //emitModifiedSignal(): causes the detectorModified signal to be emitted only
  //  if it is necessary (eg if the current detector does not equal previously
  //  emitted detector)
  void emitModifiedSignal();


  /** Searches the database for DRFs from the specified user with hashes
      matching the passed in DRF, and updates the last used time for all of
      them (should only find at most one).  If the DRF does not exist in the DB
      then it is added to is.
   */
  static void updateLastUsedTimeOrAddToDb( std::shared_ptr<DetectorPeakResponse> drf,
                                           long long db_user_id,
                                           std::shared_ptr<DataBaseUtils::DbSession> sql );
  
  /** Checks the database to see if the measurement serial number or measurement
      model cooresponds to a user preference in the database for  a DRF to use.
      If found, returns the DRF, if not, returns null.
   */
  static std::shared_ptr<DetectorPeakResponse> getUserPrefferedDetector(
                                std::shared_ptr<DataBaseUtils::DbSession> sql,
                                Wt::Dbo::ptr<InterSpecUser> user,
                                std::shared_ptr<const MeasurementInfo> meas );
  
  /**
   */
  static void setUserPrefferedDetector( std::shared_ptr<DetectorPeakResponse> drf,
                                        std::shared_ptr<DataBaseUtils::DbSession> sql,
                                        Wt::Dbo::ptr<InterSpecUser> user,
                                        UseDrfPref::UseDrfType prefType,
                                        std::shared_ptr<const MeasurementInfo> meas );
  
  void acceptAndFinish();
  void cancelAndFinish();

  //Looks through directories specified by "GadrasDRFPath" for directories that
  //  contain Detector.dat and Efficiency.csv; doesnt test that they are valid.
  //  First each pair is the path to the detector, second in the pair is a
  //  reasonable display name.
  std::vector< std::pair<std::string,std::string> > avaliableGadrasDetectors() const;

  //Will init detector in the data/detector_responses folder, and return result.
  //  throws exception if there is an error.
  //Type should be a DetectorType enum (just not including that header to save
  //   on deplandcies)
  static std::shared_ptr<DetectorPeakResponse> initAGadrasDetector( const int type, InterSpec *interspec );
  
  /** Looks to path in user prefernce "GadrasDRFPath" for a detector matching specified name.
   Throws exception in error; returned detector should always be valid.
   */
  static std::shared_ptr<DetectorPeakResponse> initAGadrasDetector( const std::string &name, InterSpec *interspec );
  
  /** Inits adetector from a directory that has a Detector.dat and Efficiency.csv file in it
      Throws exception in error; returned detector should always be valid.
   */
  static std::shared_ptr<DetectorPeakResponse> initAGadrasDetectorFromDirectory( const std::string &directory, InterSpec *interspec );

  //Will init detector in the data/detector_responses folder, and return result.
  //  throws exception if there is an error.
  //Type should be a DetectorType enum (just not including that header to save
  //   on deplandcies)
  static std::shared_ptr<DetectorPeakResponse> initARelEffDetector( const int type, InterSpec *interspec );
  
  //Callbacks for when detector is changed or modified
  void gadrasDetectorSelectCallback();
  void relEffDetectorSelectCallback();
  
  
  /** Reason #fileUpladedCallback is being called */
  enum class UploadCallbackReason
  {
    ImportTabChosen,
    DetectorDiameterChanged,
    EfficiencyCsvUploaded,
    DetectorDotDatUploaded
  };//enum class UploadCallbackReason
  
  /** Function called when "Import" tab is chosen, detector diameter is changed,
      or if a Efficiency.csv, or Detector.dat file is uplaoded.
      @param context Has value 0 
   */
  void fileUploadedCallback( const UploadCallbackReason context );

  //updates energy efficient chart
  void updateChart();
  
  //Listens to change in detector from outside this dialog, like when
  //a spectrum is created.  But be careful of feedback loop.
  void setDetector( std::shared_ptr<DetectorPeakResponse> det );
 
  //initializes the dialog with appropriate detector selected
  void init();
  
  //detector(): returns m_detector.
  std::shared_ptr<DetectorPeakResponse> detector();
  
  //deletes a row from DB table
  void deleteDBTableSelected();
  
  //a row is selected, so update det and charts
  void dbTableSelectionChanged();

protected:
  void setAcceptButtonEnabled( const bool enable );
  
  /** Checks if file at passed in path is a TSV/CSV file that contains
      coeffeicents for the exp( c0 + c1*logx + c2*logx^2 + ...) equation.
      If so, returns detector.  If not, returns nullptr.
   */
  static std::shared_ptr<DetectorPeakResponse> checkIfFileIsRelEff( const std::string tsvfilepath );
  
protected:
  WContainerWidget *m_footer;
  InterSpec *m_interspec;
  SpectraFileModel *m_fileModel;
  Wt::WColor m_chartEnergyLineColor;
  Wt::WColor m_chartFwhmLineColor;
  Wt::Chart::WCartesianChart *m_chart;
  Wt::WStandardItemModel* m_efficiencyModel;
  
  //m_detector: the detector that is currently defined for the user to edit.
  //  When ever there is a substantial edit, this detector wil be emitted via
  //  the appropriate detectorChanged or detectorModified signal.
  std::shared_ptr<DetectorPeakResponse> m_detector;

  //m_previousDetectorDef: this is equal to the previously emmitted detector
  //  so that the detectorModified signals wont be emitted when its not
  //  necessary (prob faster to check this, than push updates to clients -
  //  unverified though)
  DetectorPeakResponse m_previousDetectorDef;

  //m_previousEmmittedDetector: pointer to the previously emitted detector so
  //  un-needed detectorChanged signals wont be emmitted
  std::shared_ptr<DetectorPeakResponse> m_previousEmmittedDetector;


  //m_originalDetector: A pointer to the original detector passed into the
  //  constructor.
  std::shared_ptr<DetectorPeakResponse> m_originalDetector;

  //m_originalDetectorCopy: a _copy_ of the detector passed in.  If the user
  //  hits "cancel", then the detectorChanged signal is emitted with the
  //  original detector pointer, but first set equal to m_originalDetectorCopy
  //  incase there has been any changes made
  std::shared_ptr<DetectorPeakResponse> m_originalDetectorCopy;
  
  Wt::Signal<> m_doneSignal;


  Wt::WTabWidget *m_tabs;

  Wt::WLineEdit *m_detectorDiameter;
  Wt::WContainerWidget *m_detectrDiameterDiv;
  Wt::WFileUpload *m_efficiencyCsvUpload;
  Wt::WFileUpload *m_detectorDotDatUpload;

  Wt::WPushButton *m_acceptButton;
  Wt::WPushButton *m_cancelButton;
  
  Wt::WLineEdit    *m_detectorManualFunctionName;
  Wt::WTextArea    *m_detectorManualFunctionText;
  Wt::WLineEdit    *m_detectorManualDescription;
  Wt::WButtonGroup *m_eqnEnergyGroup;
  Wt::WButtonGroup *m_absOrIntrinsicGroup;
  Wt::WLineEdit    *m_detectorManualDiameterText;
  Wt::WLineEdit    *m_detectorManualDistText;
  Wt::WLabel       *m_detectorManualDistLabel;
  Wt::WPushButton  *m_manualSetButton;
  
  GadrasDetSelect *m_gadrasDetSelect;
  RelEffDetSelect *m_relEffSelect;
  
  Wt::WMenu *m_drfTypeMenu;
  Wt::WStackedWidget *m_drfTypeStack;
  
  Wt::WPushButton* m_deleteButton;
  Wt::WTreeView* m_DBtable;
  
  std::shared_ptr<DataBaseUtils::DbSession> m_sql;
  Wt::Dbo::QueryModel< Wt::Dbo::ptr<DetectorPeakResponse> >    *m_model;

  Wt::WCheckBox *m_defaultForSerialNumber;
  Wt::WCheckBox *m_defaultForDetectorModel;
};//class DetectorEdit


class DetectorEditWindow : public AuxWindow
{
/*
  Provides a window that contains a DetectorEdit
*/
public:
  DetectorEditWindow( std::shared_ptr<DetectorPeakResponse> currentDet,
                InterSpec *specViewer,
                SpectraFileModel *fileModel );
  virtual ~DetectorEditWindow();

protected:
  static void acceptAndDelete( DetectorEditWindow *window );

  DetectorEdit *m_edit;
};//class DetectorEditWindow


#endif //DetectorEdit_h
