#ifndef ShieldingSourceDisplay_h
#define ShieldingSourceDisplay_h
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

#include <map>
#include <tuple>
#include <mutex>
#include <utility>

#if( INCLUDE_ANALYSIS_TEST_SUITE )
#include <boost/optional.hpp>
#endif

#include <Wt/WRectF>
#include <Wt/WColor>
#include <Wt/WPainter>
#include <Wt/WModelIndex>
#include <Wt/WGridLayout>
#include <Wt/WContainerWidget>
#include <Wt/WAbstractItemModel>
#include <Wt/Chart/WCartesianChart>

#include "InterSpec/DetectorPeakResponse.h" //DetectorPeakResponse::EffGeometryType
#include "InterSpec/ShieldingSourceFitCalc.h"

//Forward declarations
class PeakDef;
struct Material;
class PeakModel;
class AuxWindow;
class MaterialDB;
class ColorTheme;
class PopupDivMenu;
class SwitchCheckbox;
class DetectorDisplay;
class PopupDivMenuItem;
struct ShieldingSourceModel;
#if( INCLUDE_ANALYSIS_TEST_SUITE )
class SpectrumViewerTester;
#endif

namespace SandiaDecay
{
  struct Nuclide;
  struct Element;
}//namespace SandiaDecay

namespace ROOT
{
  namespace Minuit2
  {
    class MnUserParameters;
  }//namespace Minuit2
}//namespace ROOT

//A Forward declaration
namespace rapidxml
{
  template<class Ch> class xml_node;
  template<class Ch> class xml_document;
}//namespace rapidxml

namespace Wt
{
  class WText;
  class WLabel;
  class WAnchor;
  class WSvgImage;
  class WCheckBox;
  class WLineEdit;
  class WTreeView;
  class WComboBox;
  class WGridLayout;
  class WFileUpload;
  class WSelectionBox;
  class WSuggestionPopup;
  class WStandardItemModel;
//  namespace Chart
//  {
//    class WCartesianChart;
//  }//namespace Chart
}//namespace Wt


namespace GammaInteractionCalc
{
  enum class GeometryType : int;
  enum class ModelSourceType : int;
  class ShieldingSourceChi2Fcn;
  struct SourceFitDef;
}//namespace GammaInteractionCalc

class InterSpec;
class SourceFitModel;
class ShieldingSelect;
class NativeFloatSpinBox;
class ShieldingSourceDisplay;


/*
class ShieldingSourceConstraint : public Wt::WContainerWidget
{
public:
  enum Type
  {
    Undefined,
    ActivityToActivity,
    MassToMass,
    ShieldingThickness
  };//enum typedef


protected:
};//class ShieldingSourceConstraint : public Wt::WContainerWidget
*/







class SourceFitModel: public Wt::WAbstractItemModel
{
public:
  enum Columns
  {
    kIsotope, kFitActivity, kFitAge, kActivity, kAge, kIsotopeMass,
    kActivityUncertainty, kAgeUncertainty,
#if( INCLUDE_ANALYSIS_TEST_SUITE )
    kTruthActivity, kTruthActivityTolerance, kTruthAge, kTruthAgeTolerance,
#endif
    kNumColumns
  };//enum Columns


  //SourceFitModel constructor: the peakModel pointer must be valid.
  //  'sameAgeIsotopes' specifies whether isotopes of the same element should be
  //  forced to have the same age.
  SourceFitModel( PeakModel *peakModel,
                  const bool sameAgeIsotopes,
                  Wt::WObject *parent = 0 );
  virtual ~SourceFitModel();

  //numNuclides(): returns same as rowCount()
  int numNuclides() const;

  //row(...): get the row cooresponfding to the nuclide passed in.  Returns -1
  //  if the nuclide is not in this model, or is invalid
  int row( const SandiaDecay::Nuclide *nuclide ) const;

  //Accessor functions all taking row number according to currently sorted order
  //  All of the below will throw std::runtime_error if invald index.
  const SandiaDecay::Nuclide *nuclide( int nuc ) const;
  double activity( int nuc ) const;
  double activityUncert( int nuc ) const;
  
  /** Returns m_nuclides[nuc].fitActivity for point and trace sources, and always false for self-atten sources; note that
   this is different than data( {row,kFitActivity} ) which only returns a non-empty value for point sources.
   */
  bool fitActivity( int nuc ) const;
  
  //age(): returns IsoFitStruct::age, which is marked age for this nuclide,
  //  which if there is a defining age nuclide set for this nuclide, you should
  //  get the age for that nuclide. See ageDefiningNuclide(...) for this case.
  double age( int nuc ) const;
  double ageUncert( int nuc ) const;

#if( INCLUDE_ANALYSIS_TEST_SUITE )
  boost::optional<double> truthActivity( int nuc ) const;
  boost::optional<double> truthActivityTolerance( int nuc ) const;
  boost::optional<double> truthAge( int nuc ) const;
  boost::optional<double> truthAgeTolerance( int nuc ) const;
#endif
  
  //fitAge(): returns IsoFitStruct::fitAge, which is if it is marked to fit age
  //  for this nuclide, not if you should fit for the age of this nuclide since
  //  it may have another defining age nuclide that controlls the age, see
  //  ageDefiningNuclide(...) for this case
  bool fitAge( int nuc ) const;
  
  
  ShieldingSourceFitCalc::ModelSourceType sourceType( int nuc ) const;
  
  /** Returns true if a shielding determined source, or a trace source. */
  bool isVolumetricSource( int nuc ) const;
  
  int nuclideIndex( const SandiaDecay::Nuclide *nuc ) const;
 
  //setSharedAgeNuclide(): in order to make it so all isotopes of an element
  //  can be made to have the same age, we'll have it so one of the isotopes
  //  controls the age (and if it can be fit) for all of them in a nuclide.
  //The dependantNucs' age will be controlled by definingNuc.  Setting the definingNuc
  //  to NULL (or to same value as dependantNuc) will disable having another nuclide
  //  control this age.
  //If dependantNuc->atomicNumber != definingNuc->atomicNumber, an exception will
  //  be thrown
  void setSharedAgeNuclide( const SandiaDecay::Nuclide *dependantNuc,
                             const SandiaDecay::Nuclide *definingNuc );
  
  //ageDefiningNuclide(...): returns the nuclide that controlls the age for the
  //  passed in nuclide.  If the passed in nuclide does not have another nuclide
  //  that controls its age, it returns the same nuclide that was passed in.
  const SandiaDecay::Nuclide *ageDefiningNuclide( const SandiaDecay::Nuclide *dependantNuc ) const;
  
  
  /** Sets the source type for a nuclide.
   Sources default to ShieldingSourceFitCalc::ModelSourceType::Point, which lets the user edit the activity and whether or not to fit the activity in the table.
   
   For ShieldingSourceFitCalc::ModelSourceType::Intrinsic or ShieldingSourceFitCalc::ModelSourceType::Trace the activity field in the table will become non-editable, and the option to
   fit this quantity will become non-visible, as the ShieldingSelect controls this option via fitting for thickness and/or trace activity.
   */
  void setSourceType( const SandiaDecay::Nuclide *nuc, ShieldingSourceFitCalc::ModelSourceType type );
  
  void makeActivityNonEditable( const SandiaDecay::Nuclide *nuc );

  void peakModelRowsInsertedCallback( Wt::WModelIndex index,
                                      int firstRow, int lastRow );
  void peakModelRowsRemovedCallback( Wt::WModelIndex index,
                                     int firstRow, int lastRow );

  //peakModelDataChangedCallback(...): will throw std::runtime_error if
  //  firstRow != lastRow
  void peakModelDataChangedCallback( Wt::WModelIndex topLeft,
                                     Wt::WModelIndex bottomRight );
  void peakModelResetCallback();

  //insertPeak(...): Will not insert peak if isotope is already accounted for
  void insertPeak( const std::shared_ptr<const PeakDef> peak );

  //repopulateIsotopes(): go through and and sync this models isotopes with
  // its parent m_peakModel
  void repopulateIsotopes();

  //setUseSameAgeForIsotopes(...): sets whether isotopes of an element should
  //  all share the same age
  void setUseSameAgeForIsotopes( bool useSame );
  
  //index(...): a shorthand method of getting an index for a specfici nuclide.
  //  Returns invalid index if nuclide passed in is not in this model
  Wt::WModelIndex index( const SandiaDecay::Nuclide *n, Columns column ) const;

  //index(...): a shorthand method of getting an index for a specfici nuclide.
  //  Returns invalid index if nuclide passed in is not in this model
  Wt::WModelIndex index( const std::string &nuclide, Columns column ) const;

  //Functions defined Wt::WAbstractItemModel that this class customizes
  virtual int columnCount( const Wt::WModelIndex &p = Wt::WModelIndex() ) const;
  virtual int rowCount( const Wt::WModelIndex &p = Wt::WModelIndex() ) const;
  virtual Wt::WFlags<Wt::ItemFlag> flags( const Wt::WModelIndex &index ) const;
  virtual Wt::WModelIndex parent( const Wt::WModelIndex &index ) const;
  virtual boost::any data( const Wt::WModelIndex &index,
                           int role = Wt::DisplayRole ) const;
  boost::any headerData( int section,
                         Wt::Orientation orientation, int role ) const;
  virtual Wt::WModelIndex index( int row, int column,
                      const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;

  //setData(...) note that if activity or age is set, then the associated
  //  uncertainty for that quantity will be reset as well.
  virtual bool setData( const Wt::WModelIndex &index,
                        const boost::any &value, int role = Wt::EditRole );
  virtual void sort( int column, Wt::SortOrder order = Wt::AscendingOrder );


  //compare(...): function to compare IsoFitStruct according to relevant Column;
  //  functions similar to operator<
  static bool compare( const ShieldingSourceFitCalc::IsoFitStruct &lhs,
                      const ShieldingSourceFitCalc::IsoFitStruct &rhs,
                       Columns sortColumn, Wt::SortOrder order );

  void displayUnitsChanged( bool displayBq );
  
  const std::vector<ShieldingSourceFitCalc::IsoFitStruct> &underlyingData() const;
  
  void setUnderlyingData( const std::vector<ShieldingSourceFitCalc::IsoFitStruct> &data );
  
  /** Set the DetectorPeakResponse::EffGeometryType, so correct label for activity can be shown.
   */
  void setDetectorType( const DetectorPeakResponse::EffGeometryType det_type );
  DetectorPeakResponse::EffGeometryType detType() const;
  
protected:
  Wt::SortOrder m_sortOrder;
  Columns m_sortColumn;
  bool m_displayCuries;
  PeakModel *m_peakModel;
  std::vector<ShieldingSourceFitCalc::IsoFitStruct> m_nuclides;
  bool m_sameAgeForIsotopes;
  DetectorPeakResponse::EffGeometryType m_det_type;
  
  //m_previousResults: when a isotope gets removed from this model, we'll cache
  //  its current value, since it will often times get added again and be
  //  intended to be the same value
  std::map<const SandiaDecay::Nuclide *, ShieldingSourceFitCalc::IsoFitStruct> m_previousResults;
  
  friend class ShieldingSourceDisplay;
};//class SourceFitModel


class ShieldingSourceDisplay : public Wt::WContainerWidget
{
public:
  ShieldingSourceDisplay( PeakModel *peakModel,
                          InterSpec *specViewer,
                          Wt::WSuggestionPopup *materialSuggest,
                          MaterialDB *materialDB,
                          Wt::WContainerWidget *parent = 0 );
  
  /** Creates a AuxWindow with a ShieldingSourceDisplay in it.
   
   @returns the created ShieldingSourceDisplay and AuxWindow. If for some reason there was an issue
   making the widgets, the returned pair will be nullptr's (and also an error message will be
   displayed to the user).
   */
  static std::pair<ShieldingSourceDisplay *,AuxWindow *> createWindow( InterSpec *viewer  );
  
  virtual ~ShieldingSourceDisplay() noexcept(true);

#if( INCLUDE_ANALYSIS_TEST_SUITE )
  /** Creates a window that lets you enter truth-values and tolerances for all the quantities
   currently marked to be fit for.
   */
  void showInputTruthValuesWindow();
  
  /** Sets activitities to like 1mCi, and thicknesses to 1 cm, so this was test fits wont be
   starting at an already correct value.
   */
  void setFitQuantitiesToDefaultValues();
  
  /** Returns <num values specified, num fit values, is valid> */
  std::tuple<int,int,bool> numTruthValuesForFitValues();
  
  /** Renders the Chi2Chart to a SVG image */
  void renderChi2Chart( Wt::WSvgImage &image );

  /** Tests the current values, for all quantities being fit for, against truth-level values.
   
   @returns <IfTestSuccesful,NumCorrect,NumTested,TextInfoLines>
   */
  std::tuple<bool,int,int,std::vector<std::string>> testCurrentFitAgainstTruth();
#endif
  
  //add generic shielding
  void addGenericShielding();
  
  /** Creates a new shielding, and hooks up all the necessary signals.
   
   @param before If nullptr, the shielding will be added after the last ShieldingSelect; if non-null, the new ShieldingSelect will be added
                before the specified ShieldingSelect.
   @param updateChiChartAndAddUndoRedo If true, the Chi2 chart will be updated after adding, AND a undo/redo point will be
                added.
  */
  ShieldingSelect *addShielding( ShieldingSelect *before, const bool updateChiChartAndAddUndoRedo );
  
  //doAddShielding(), doAddShieldingBefore(), doAddShieldingAfter:
  //  convience functions that calls addShielding(...) appropriately.
  void doAddShielding();
  void doAddShieldingBefore( ShieldingSelect *select );
  void doAddShieldingAfter( ShieldingSelect *select );
  size_t numberShieldings() const;
  
  void removeShielding( ShieldingSelect *select );
  void materialModifiedCallback( ShieldingSelect *select );
  void materialChangedCallback( ShieldingSelect *select );
  void updateActivityOfShieldingIsotope( ShieldingSelect *select,
                                         const SandiaDecay::Nuclide *nuc );
  void isotopeIsBecomingVolumetricSourceCallback( ShieldingSelect *select,
                                           const SandiaDecay::Nuclide *nuc,
                                           const ShieldingSourceFitCalc::ModelSourceType type );
  void isotopeRemovedAsVolumetricSourceCallback( ShieldingSelect *select,
                                             const SandiaDecay::Nuclide *nuc,
                                             const ShieldingSourceFitCalc::ModelSourceType type );
  
  /** Function called whenever a ShieldingSelect is changed by the user.
   
   This function is not called when values are changed programmatically (e.g., as a side-effect of another setting being changed).
   
   Storing the XML for the ShieldingSelect, as apposed to the whole ShieldingSourceDisplay, takes up only about 3% as much
   memmory, but at the cost of not currently getting things exactly right for trace and self-atten sources.
   
   @param select The ShieldingSelect that was changed.
   @param prev_state The XML representation of the ShieldingSelect before it was changed.
   @param current_state The XML representation of the updated ShieldingSelect state.
   */
  void handleShieldingUndoRedoPoint( const ShieldingSelect * const select,
                                    const std::shared_ptr<const std::string> &prev_state,
                                    const std::shared_ptr<const std::string> &current_state );
  
  /** Performs the actual fit of shielding, activities, and ages;
   called when user clicks "Perform Model Fit" button.
   
   @param fitInBackground If true, the fit is performed in the background, and function will return
          immediately (i.e., before fit is performed).  If false, function will not return until
          fit is finished (i.e., blocking).
   @returns A pointer to the fitting results; will be nullptr if fit was not started.  Note that
            if fitting is being performed in the background, the returned object will be updated as
            fitting is being performed (use m_mutex to ensure safe access).
   */
  std::shared_ptr<ShieldingSourceFitCalc::ModelFitResults> doModelFit( const bool fitInBackground );
  
  /** Cancels the current model fit happening */
  void cancelModelFit();
  
  /** Cancels the current model fit happening, and will make it so the `gui_updater` argument of
   `ShieldingSourceFitCalc::fit_model` will not be called, but the GUI will be put back into
   a state that the user can use it.
   
   This method is intended for the undo/redo mechanism, in case a user hits undo while a fit is happening.
   */
  void cancelModelFitWithNoUpdate();
  
  /** Function that must be called from within the main Wt event loop thread,
   that updates the GUI with the fit progress.
   */
  void updateGuiWithModelFitProgress( std::shared_ptr<ShieldingSourceFitCalc::ModelFitProgress> progress );
  
  /** Function that must be called from within the main Wt event loop thread,
      that updates the GUI with the fit results.
   */
  void updateGuiWithModelFitResults( std::shared_ptr<ShieldingSourceFitCalc::ModelFitResults> results );
  
  //updateCalcLogWithFitResults(): adds in fit for values and their
  //  uncertainties in the calculation log.
  void updateCalcLogWithFitResults( std::shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn> chi2,
                                    std::shared_ptr<ShieldingSourceFitCalc::ModelFitResults> results );
  
  //initialSizeHint(...) gives the initial hint about the size so which widgets
  //  should be shown can be decided on.  Furthermore, calling this function
  //  resets m_nResizeSinceHint, so that the next 4 calls to
  //  layoutSizeChanged(...) will be ignored, since this is how many times this
  //  function typically gets called when creating an AuxWindow with a
  //  ShieldingSourceDisplay as its main widget.
  void initialSizeHint( int width, int height );
  
  void handleUserDistanceChange();
  
  GammaInteractionCalc::GeometryType geometry() const;
  void handleGeometryTypeChange();
  
  //checkAndWarnZeroMassFraction(): makes sure the user hasnt asked to use
  //  a shielding as a source, but specified a zero mass fraction.  Throws
  //  std::exception on error, with a message approprtiate for output to user.
  void checkAndWarnZeroMassFraction();
  
  //checkDistanceAndThicknessConsistent(): called when the m_distanceEdit is
  //  changed; makes sure radius of outer shielding is less than this distnace,
  //  and updates chi2 chart
  void checkDistanceAndThicknessConsistent();

  /** Checks to see if fitting for more than one atomic number of generic shielding.
   TODO: Could probably tighten things up a bit to avoid degeneracies
   */
  void checkForMultipleGenericMaterials();
  
  void handleShieldingChange();
  
  void handleDetectorChanged( std::shared_ptr<DetectorPeakResponse> new_det );
  
  void updateChi2Chart();
  
  void showCalcLog();
  void closeCalcLogWindow();
  
  /** Returns the inner ShieldingSelect of the one passed in; e.g., returns the ShieldingSelect that is contained by the one passed in.
   
   Returns nullptr if the inner-most shielding.
   
   Will throw exception if the ShieldingSelect is not a valid pointer owned by this ShieldingSourceDisplay.
   */
  const ShieldingSelect *innerShielding( const ShieldingSelect * const select ) const;
  
  
  //testSerialization(): simply tries to round-trip this ShieldingSourceDisplay
  //  to and then back from XML.  Does not actually check it was correctly done
  //  but does print the XML to cerr, and will trigger exceptions being thrown
  //  if XML isnt as expected; intended to be used for development.
  void testSerialization();
  
  //serialize(): creates (and returns) a node "ShieldingSourceFit" under
  //  'parent' which contains the XML for this object.
  ::rapidxml::xml_node<char> *serialize( ::rapidxml::xml_node<char> *parent );
  
  //deSerialize(): takes in a "ShieldingSourceFit" node and sets the state of
  //  this object to match the XML.
  //Throws when it runs into an unexpected situation, or invalid parent_node
  void deSerialize( const ::rapidxml::xml_node<char> *parent_node );
  
  //some helper functions for deserialization:
  void deSerializePeaksToUse( const ::rapidxml::xml_node<char> *peaks );
  void deSerializeSourcesToFitFor( const ::rapidxml::xml_node<char> *sources );
  void deSerializeShieldings( const ::rapidxml::xml_node<char> *shiledings );
  
  void startModelUpload();
  void modelUploadError( const ::int64_t size_tried );
  void finishModelUpload( Wt::WFileUpload *upload );
  void closeModelUploadWindow();
  
#if( USE_DB_TO_STORE_SPECTRA )
  void startSaveModelToDatabase( bool promptNewTitleOnExistingEntry );
  void closeSaveModelToDatabaseWindow();
  
  //finishSaveModelToDatabase(...): Wont throw exception, and instead returns
  //  the success status.
  //  Will set m_modelInDb to the result; if saving fails m_modelInDb will be
  //  set to NULL.
  bool finishSaveModelToDatabase( const Wt::WString &name,
                                  const Wt::WString &desc );
  void finishGuiSaveModelToDatabase( Wt::WLineEdit *name_edit,
                                  Wt::WLineEdit *desc_edit );
  void saveCloneModelToDatabase();
  
  //saveModelIfAlreadyInDatabase(...): wont throw, but may silently fail
  void saveModelIfAlreadyInDatabase();
#endif //#if( USE_DB_TO_STORE_SPECTRA )
  
  std::string defaultModelName() const;
  std::string defaultModelDescription() const;
  
#if( USE_DB_TO_STORE_SPECTRA )
  void startBrowseDatabaseModels();
  void closeBrowseDatabaseModelsWindow();
  void finishLoadModelFromDatabase( Wt::WSelectionBox *selected,
                                    Wt::WSelectionBox *other_select );
  
  //loadModelFromDb(...):  deserializes the model passed in.
  //  Will not throw, but may fail to load the passed in entry, in which
  //  case it will attempt to not make any changes (although not strictly
  //  garunteed) to current model, and inform user of the error and return
  //  the state of loading the passed in entry.
  bool loadModelFromDb( Wt::Dbo::ptr<ShieldingSourceModel> entry );
  
  void removeModelFromDb( Wt::WSelectionBox *selec1,
                          Wt::WSelectionBox *selec2 );
  
  //modelInDb(): will resave m_modelInDb (if valid) and return it.  Otherwise
  //  if m_modelInDb is not valid, and it appears the current model is
  //  non-empty, the current model will be saved (with default name/description)
  //  and returned.
  //  Shouldnt throw.
  Wt::Dbo::ptr<ShieldingSourceModel> modelInDb();
#endif //#if( USE_DB_TO_STORE_SPECTRA )
  
  //reset(): removes all shielding and sources; removes all peaks from fit
  void reset();

  //newForegroundSet(): call this function to reset the 'm_modifiedThisForeground'
  //  variable, so we can kinda track if the current state of this model should
  //  be associated with the current foregorund.
  void newForegroundSet();
  
  //userChangedDuringCurrentForeground(): lets you know if the user has changed
  //  this foreground at all by adding or removing shielding, fitting for
  //  shielding, uploading a model, or if they modified a shielding since the
  //  last time newForegroundSet() has been called.
  bool userChangedDuringCurrentForeground() const;
  
  //addSourceIsotopesToShieldings(): called _after_ an isotope is added to
  //  m_sourceModel
  void addSourceIsotopesToShieldings( Wt::WModelIndex,
                                      int firstRow, int lastRow );

  //removeSourceIsotopesFromShieldings(): called _before_ the isotopes is
  //  removed from m_sourceModel
  void removeSourceIsotopesFromShieldings( Wt::WModelIndex,
                                           int firstRow, int lastRow );

  void multiNucsPerPeakChanged();
  void attenuateForAirChanged();
  void backgroundPeakSubChanged();
  void sameIsotopesAgeChanged();
  void decayCorrectChanged();
  void showGraphicTypeChanged();
  

  std::pair<std::shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>,
            ROOT::Minuit2::MnUserParameters> shieldingFitnessFcn();
  

  //toggle checkbox/chart
  void toggleUseAll(Wt::WCheckBox* button);
  void updateAllPeaksCheckBox(  Wt::WCheckBox *but);
    
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags ) override;

  
  SourceFitModel *sourceFitModel();
  
  ShieldingSourceFitCalc::ShieldingSourceFitOptions fitOptions() const;
protected:
  /** Disables fit button and other elements, and hides some stuff, etc. */
  void setWidgetStateForFitStarting();
  
  /** Undoes the changes from #setWidgetStateForFitStarting */
  void setWidgetStateForFitBeingDone();
  
  void updateChi2ChartActual();
  virtual void layoutSizeChanged( int width, int height ) override;
  
protected:
  class Chi2Graphic;  //forward declaration

  bool m_chi2ChartNeedsUpdating;
  
  int m_width, m_height, m_nResizeSinceHint;

  //m_modifiedThisForeground: tracks if the user has modified this
  //  ShieldingSourceDisplay for the current foreground.
  //As of 20140706 it hasnt been tested that this state is actually properly
  //  tracked.
  bool m_modifiedThisForeground;

  PeakModel *m_peakModel;       //belongs to m_specViewer, not *this
  InterSpec *m_specViewer;
  SourceFitModel *m_sourceModel;

  Wt::WTreeView *m_peakView;
  Wt::WTreeView *m_sourceView;

  DetectorDisplay *m_detectorDisplay;
  
  Wt::WLabel *m_distanceLabel;
  /** Is set after validating a user entered distance string.  If user entered
   *  string is invalid, then this value is used to go back to.
   */
  std::string m_prevDistStr;
  Wt::WLineEdit *m_distanceEdit;

  Wt::WPushButton *m_addMaterialShielding;
  Wt::WPushButton *m_addGenericShielding;
  
  Wt::WGridLayout *m_layout;
  PopupDivMenu *m_addItemMenu;
  
#if( USE_DB_TO_STORE_SPECTRA )
  PopupDivMenuItem *m_saveAsNewModelInDb;
#endif
  
  Wt::Dbo::ptr<ShieldingSourceModel> m_modelInDb;
  
  Wt::WSuggestionPopup *m_materialSuggest;

  //m_shieldingSelects: contains objects of class ShieldingSelect
  Wt::WContainerWidget *m_shieldingSelects;

  Wt::WLabel *m_geometryLabel;
  GammaInteractionCalc::GeometryType m_prevGeometry;
  Wt::WComboBox *m_geometrySelect;

  Wt::WText *m_fixedGeometryTxt;
  
  Wt::WText *m_showChi2Text;
  Wt::WStandardItemModel *m_chi2Model;
  Chi2Graphic *m_chi2Graphic;
  

  Wt::WCheckBox  *m_multiIsoPerPeak;
  Wt::WCheckBox  *m_attenForAir;
  Wt::WCheckBox  *m_backgroundPeakSub;
  Wt::WCheckBox  *m_sameIsotopesAge;
  Wt::WCheckBox  *m_decayCorrect;
  SwitchCheckbox *m_showChiOnChart;
  Wt::WContainerWidget *m_optionsDiv;
  double m_photopeak_cluster_sigma;
  bool m_multithread_computation;
  
  PopupDivMenuItem *m_showLog;
  
  AuxWindow *m_logDiv;
  std::vector<std::string> m_calcLog;
  
  AuxWindow *m_modelUploadWindow;
#if( USE_DB_TO_STORE_SPECTRA )
  AuxWindow *m_modelDbBrowseWindow;
  AuxWindow *m_modelDbSaveWindow;
#endif
  
  //m_materialDB: not owned by this object, but passed in at construction.
  MaterialDB *m_materialDB;


  Wt::WPushButton *m_fitModelButton;
  Wt::WText *m_fitProgressTxt;
  Wt::WPushButton *m_cancelfitModelButton;
  
  std::mutex m_currentFitFcnMutex;  //protects the shared_ptr only, not the object it points to
  std::shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn> m_currentFitFcn;
  
  //A class to draw the chi2 distribution of the fit to activity/shielding.
  //  We have to overide the Paint(...) method to draw some text on chart
  //  indicating the chi2
  class Chi2Graphic : public Wt::Chart::WCartesianChart
  {
  public:
    Chi2Graphic( Wt::WContainerWidget *parent = 0 );
    virtual ~Chi2Graphic();
    virtual void paint( Wt::WPainter &painter,
                        const Wt::WRectF &rectangle = Wt::WRectF() ) const;
    virtual void paintEvent( Wt::WPaintDevice *paintDevice );
    void setNumFitForParams( unsigned int npar );
    
    void setShowChiOnChart( const bool show_chi );
    void setTextPenColor( const Wt::WColor &color );
    void setColorsFromTheme( std::shared_ptr<const ColorTheme> theme );
  protected:
    void calcAndSetAxisRanges();
    void calcAndSetAxisPadding( double yHeightPx );
    
    int m_nFitForPar;
    bool m_showChi;
    Wt::WColor m_textPenColor;
  };//class WCartesianChart

  static const int sm_xmlSerializationMajorVersion;
  static const int sm_xmlSerializationMinorVersion;
  
#if( INCLUDE_ANALYSIS_TEST_SUITE )
  // So the tester can call updateChi2ChartActual
  friend class SpectrumViewerTester;
#endif
};//class ShieldingSourceDisplay



#endif //ifndef ShieldingSourceDisplay_h
