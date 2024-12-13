#ifndef ShieldingSelect_h
#define ShieldingSelect_h
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

#include "InterSpec/ShieldingSourceFitCalc.h"

//Forward declarations
class PeakDef;
struct Material;
class PeakModel;
class AuxWindow;
class MaterialDB;
class PopupDivMenu;
class SwitchCheckbox;
class DetectorDisplay;
class PopupDivMenuItem;
class DetectorPeakResponse;
struct ShieldingSourceModel;
#if( INCLUDE_ANALYSIS_TEST_SUITE )
class SpectrumViewerTester;
#endif

namespace SandiaDecay
{
  struct Nuclide;
  struct Element;
}//namespace SandiaDecay


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
  class WCheckBox;
  class WLineEdit;
  class WGridLayout;
  class WSuggestionPopup;
}//namespace Wt


class InterSpec;
class SourceFitModel;
class TraceSrcDisplay;
class NativeFloatSpinBox;
class ShieldingSourceDisplay;

namespace GammaInteractionCalc
{
  enum class ModelSourceType : int;
}

namespace GammaInteractionCalc
{
  enum class GeometryType : int;
  enum class TraceActivityType : int;
}//namespace GammaInteractionCalc


class SourceCheckbox : public Wt::WContainerWidget
{
public:
  SourceCheckbox( const SandiaDecay::Nuclide *nuclide,
                  double massFrac, Wt::WContainerWidget *parent = 0 );
  virtual ~SourceCheckbox();

  double massFraction() const;
  
  /** Sets the mass fraction and uncertainty.
   If there is no uncertainty, pass in 0.0 or less.
   Uncertainty is (currently) only displayed in the tool-tip of this widget.
   */
  void setMassFraction( double frac, double uncert );
  const SandiaDecay::Nuclide *isotope() const;
  bool useAsSource() const;
  void setUseAsSource( bool use );

  bool fitMassFraction() const;
  void setFitMassFraction( const bool fit );
  
  Wt::EventSignal<> &checked();
  Wt::EventSignal<> &unChecked();
  Wt::Signal<float> &massFractionChanged();

  Wt::EventSignal<> &fitMassFractionChecked();
  Wt::EventSignal<> &fitMassFractionUnChecked();
  
protected:
  void handleUseCbChange();
  void handleFitMassFractionChanged();
  
  Wt::WCheckBox *m_useAsSourceCb;
  Wt::WLabel *m_label;
  NativeFloatSpinBox *m_massFraction;
  Wt::WCheckBox *m_fitFraction;
  const SandiaDecay::Nuclide *m_nuclide;
};//class SourceCheckbox


class ShieldingSelect : public Wt::WContainerWidget
{
public:
  /** ShieldingSelect constructor for use when you arent fitting material thickness or AN/AD - e.g., everywhere except in
   ShieldingSourceDisplay.
   */
  ShieldingSelect( MaterialDB *materialDB,
                   Wt::WSuggestionPopup *materialSuggest,
                   Wt::WContainerWidget *parent = 0 );
  
  //ShieldingSelect constructor: if forFitting==true, then the checkboxes
  //  that tell if a quantity should be fit for, are constructed.  if
  //  forFitting==false, then any calls; also if true then cursor focus will
  //  to the material edit
  //
  //  If materialSuggest does not have an object name, one will be assigned to
  //  it, to allow safe removing of the suggest from the edit in the destructor.
  ShieldingSelect( MaterialDB *materialDB,
                   SourceFitModel *sourceModel,
                   Wt::WSuggestionPopup *materialSuggest,
                   const ShieldingSourceDisplay *shieldSource,
                   Wt::WContainerWidget *parent = 0 );
  
  virtual ~ShieldingSelect();

  /** Sets the current geometry.
   If a generic material will have no visual effect.
   */
  void setGeometry( GammaInteractionCalc::GeometryType type );
  
  /** @returns the current geometry. */
  GammaInteractionCalc::GeometryType geometry() const;
  
  /** Sets if trace or self-attenuating sources are allowed.
   
   @param fixed_geom If true, then any trace or self-attenuating sources, as well as menu
          item to add a trace source, will be removed, and geometry will be set to spherical.
          If false, will make sure menu item to add trace-source is enabled.
   
   If `m_forFitting` is false, and this function is called, will throw exception.
   */
  void setFixedGeometry( const bool fixed_geom );
  
  //isGenericMaterial(): tells you if the material is defined by areal density
  //  and atomic number, or if a pre-defined material
  bool isGenericMaterial() const;

  /** @returns sphere thickness in SandiaDecay/PhysicalUnits units.  Note that this is the thickness, and not the radius (which is the
               sum of all thicknesses of ShieldingSelects before this one.
  
   Throws std::runtime_error if a GenericMaterial, non-spherical, or invalid text entered.
   */
  double thickness() const;

  
  double cylindricalRadiusThickness() const;
  double cylindricalLengthThickness() const;
  double rectangularWidthThickness() const;
  double rectangularHeightThickness() const;
  double rectangularDepthThickness() const;
  
  //atomicNumber(): returns values of 1.0 through 100.0
  //  throws std::runtime_error if not a GenericMaterial or unable to convert
  //  text into a number.  If text is blank returns 26.0.
  double atomicNumber() const;

  //arealDensity(): retuns dimension of [Mass/Length^2] in
  //  SandiaDecay/PhysicalUnits units
  //  throws std::runtime_error if not a GenericMaterial or unable to convert
  //  text into a number.  If text is blank, returns 0.
  double arealDensity() const;

  //fitAtomicNumber():
  //  throws std::runtime_error if not a GenericMaterial or unable to convert
  //  text into a number, or if m_forFitting==false
  bool fitAtomicNumber() const;
  
  //fitArealDensity():
  //  throws std::runtime_error if not a GenericMaterial or unable to convert
  //  text into a number, or if m_forFitting==false
  bool fitArealDensity() const;

  
  //fitThickness():
  //  throws std::runtime_error if a GenericMaterial, m_forFitting==false, or not spherical
  bool fitThickness() const;

  bool fitCylindricalRadiusThickness() const;
  bool fitCylindricalLengthThickness() const;
  
  bool fitRectangularWidthThickness() const;
  bool fitRectangularHeightThickness() const;
  bool fitRectangularDepthThickness() const;
  
  
  /** Sometime we know we dont have the power to fit for shielding dimensions, so in these cases we wont the user have the option.
   
   These setFit**Enabled(...) functions will unset the fit option and hide the checkbox.
   
   Will throw exception if you call the function for a geometry different than the currently set one, or if not for fitting.
   */
  void setFitCylindricalRadiusEnabled( const bool allow );
  void setFitCylindricalLengthEnabled( const bool allow );
  
  void setFitRectangularHeightEnabled( const bool allow );
  void setFitRectangularWidthEnabled( const bool allow );
  
  //material(): returns the currently selected Material.  Note that if the
  //  GenericMaterial is selected, will return NULL.
  //  If the material description has been changed, then a new Material will
  //  be created in memorry.  Otherwise, the previously existing Material in
  //  memmorry will be returned - unmodified.  Call handleMaterialChange() and
  //  handleIsotopicChange(...) to deal with modifying materials.
  // TODO: should the returned pointer be made const, or maybe unique?
  std::shared_ptr<const Material> material();
  
  /** Returns the current m_currentMaterial.
   
   Does not update material if description has been changed.
   
   Note that although the pointer returned is const, it points to the same object as m_currentMaterial,
   which is non-const, so may be modified.
   */
  std::shared_ptr<const Material> currentMaterial() const;

  //material( string ) returns the material for the text passed in.  Does not
  //  change the material represented by the widget, but will add the material
  //  passed in to the possible suggestions if it is a chemical formula.
  //  Returns a pointer to an object managed by m_materialDB.
  //  Returns null on error.
  const Material *material( const std::string &text );
  
  //remove() is the signal emitted when the user clicks the close button
  Wt::Signal<ShieldingSelect *> &remove();

  //Signals emitted when the user requests to add another shielding. Can only be
  //  fired if this ShieldingSelect is for fitting.  Argument will be *this.
  Wt::Signal<ShieldingSelect *> &addShieldingBefore();
  Wt::Signal<ShieldingSelect *> &addShieldingAfter();
  
  //materialModified() is signal emitted when the material is modified.
  Wt::Signal<ShieldingSelect *> &materialModified();

  //materialChanged() is signal emitted when the material is changed.
  Wt::Signal<ShieldingSelect *> &materialChanged();

  //  The removingIsotopeAsSource() signal will be emitted for all
  //  checked source isotopes if a new material has been selected

  //activitiesFromThicknessNeedUpdating(): signal emitted when the activity
  //  of the isotope (with symbol given by second argument) needs updating
  //  to to a change in thickness or something
  Wt::Signal<ShieldingSelect *,const SandiaDecay::Nuclide *>
                                           &activityFromVolumeNeedUpdating();

  //addingIsotopeAsSource(): Signal emitted when isotope is checked
  Wt::Signal<const SandiaDecay::Nuclide *,ShieldingSourceFitCalc::ModelSourceType> &addingIsotopeAsSource();

  //removingIsotopeAsSource(): Signal emitted when isotope is unchecked
  Wt::Signal<const SandiaDecay::Nuclide *,ShieldingSourceFitCalc::ModelSourceType> &removingIsotopeAsSource();

  /** Signal to indicate when the user has chnaged any of the values displayed.
   Used for undo/redo support.
   */
  Wt::Signal<ShieldingSelect *,
              std::shared_ptr<const std::string> /* Previous state xml */,
              std::shared_ptr<const std::string> /* New state */>
             &userChangedStateSignal();
  
  //Will un-check the checkbox for the nuclide to be a self-attenuating source, or a trace source.
  //  Does nothing if there is no isotope cb, or trace source with the nuclide.
  //Note: does not emit the removingIsotopeAsSource() signal
  void uncheckSourceIsotopeCheckBox( const SandiaDecay::Nuclide *nuc );

  //sourceRemovedFromModel(...): removes the checkbox option for
  //  given isotope, and updates possible trace source options.
  void sourceRemovedFromModel( const SandiaDecay::Nuclide *nuc );

  //setMassFractionDisplaysToMaterial(): updates displayed mass fractions to that of
  //  the Material passed in.
  void setMassFractionDisplaysToMaterial( std::shared_ptr<const Material> mat );

  /** Adjusts the "other" non-src nuclide row entry, to reflect one mins checked source nuclides fraction */
  void updateSelfAttenOtherNucFraction();
  
  /** Returns the isotopes currently checked for use as self-attenuating sources. 
   
   Note: this returns nuclides that are selected to fit mass fraction, as well as not.
   */
  std::vector<const SandiaDecay::Nuclide *> selfAttenNuclides() const;
  
  /** Returns the elements who the "other" component is selected to be fit. */
  std::vector<const SandiaDecay::Element *> elementsFittingNonSourceComponent() const;
  
  /** Gets the (self-attenuating) source nuclides (not trace nuclides), and their respective mass fractions.
   
   Mass-fractions are for the element, not the entire material.
   Note: mass-fractions are retrieved from the GUI state - they are not the values in `m_currentMaterial`.
   Sources not being used as self-atten sources are not included.
   */
  typedef std::tuple<const SandiaDecay::Nuclide *,double,bool> NucMasFrac;
  std::map<const SandiaDecay::Element *,std::vector<NucMasFrac>> sourceNuclideMassFractions() const;

  //setClosableAndAddable(...): by default widget will have close icon and emit
  //  the remove() signal, but setClosable() lets you change this
  void setClosableAndAddable( bool closeable , Wt::WGridLayout* layout);

  //fitForAnyMassFractions(): returns if the user has asked to use this material
  //  as at least two sources, and checked that they want to fit for mass fracs
  bool fitForAnyMassFractions() const;
  
  //setMassFraction(...): set the mass fraction for a given nuclide.  If nuclide
  //  is not currently being used as a source, or fraction is not between
  //  0 and 1.0, an exception will be thrown.
  //  This function does not check if all mass fractions add up to 1.0, or
  //  anything else like that, so be careful
  void setMassFraction( const SandiaDecay::Nuclide * const nuc,
                       const double fraction, const double uncert );
  
  /** How to hold information about each nuclide, of an element, for setting mass fraction information. */
  struct MassFracInfo
  {
    /**  If nuclide is nullptr, than that indicates, it is the "other" (non-source) portion of the element. */
    const SandiaDecay::Nuclide *m_nuclide = nullptr;
    double m_fraction = 0.0;
    /** Uncertainty on mass fraction; if not applicable set as zero or negative value. */
    double m_frac_uncert = 0.0;
    bool m_use_as_source = false;
    bool m_fit_mass_frac = false;
  };//struct MassFracInfo
  
  
  /** Sets all the source mass fractions at once; will normalize total mass fractions for each element to 1.0.
   If mass fractions are less than zero, or if no current material, will through exception.
   */
  void setMassFractions( std::map<const SandiaDecay::Element *,std::vector<MassFracInfo>> fractions );
  
  
  //setMaterialNameAndThickness(...): sets the current material name and
  //  thickness to the specified strings.  If necessary will make it so this is
  //  not a GenericMaterial.  Throws exception if the material cant be found
  //  in the database, or the thickness is an invalid thickness (gui wont be
  //  changed in that case).  Also calls handle material change.
  void setMaterialNameAndThickness( const std::string &name,
                                    const std::string &thickness );
  
  /** Toggles the widget to be a generic widget (if it wasnt already) and sets
     the atomic number and areal density.
     AN sould be between 1 and 100 (inclusive), and AD in PhysicalUnits g/cm2,
     and a positive value.  If errors, exception will be thrown.
   */
  void setAtomicNumberAndArealDensity( const double an, const double ad );
  
  /** Similar to other variant of this function, but sets the display string
   to the value passed in so exact value is preserved.
   
   Strings must either be blank, or valid numbers.
   Areal density is in units of g/cm2 (i.g., between 0 and 500).
   */
  void setAtomicNumberAndArealDensity( const std::string &an, const std::string &ad );
  
  void setToNoShielding();
  
  /** Returns true if material field is empty, or [thickness|AD] is empty or zero.
   
   For cylindrical and rectangular geometries, only checks dimension is empty or zero, for the one between center and detector.
   */
  bool isNoShielding();
  
  //Simple accessors
  Wt::WLineEdit *materialEdit();
  Wt::WLineEdit *thicknessEdit();
  
  const Wt::WLineEdit *materialEdit() const;
  const Wt::WLineEdit *thicknessEdit() const;
  
  /** Returns all thickness edits, for the current geometry. */
  std::vector<Wt::WLineEdit *> distanceEdits();

  /** Sets the spherical thickness value.
   
   If a negative value is passed in, then thickness input text is set to blank.
   
   Note the thickness is rounded to have at most three places after the decimal point.
   
   Note: #handleMaterialChange function will not be called; if you need to trigger the appropriate materialChanged() or
   materialModified() signals.
   
   Will throw exception if non-spherical geometry.
   */
  void setSphericalThickness( const double thickness );
  
  void setCylindricalRadiusThickness( const double radius );
  void setCylindricalLengthThickness( const double length );
  void setRectangularWidthThickness( const double width );
  void setRectangularHeightThickness( const double height );
  void setRectangularDepthThickness( const double depth );
  
  void setSphericalThicknessEditEnabled( const bool enabled );
  
  NativeFloatSpinBox *arealDensityEdit();
  NativeFloatSpinBox *atomicNumberEdit();
  
  /** Returns all the nuclides this shielding is a trace source for. */
  std::vector<const SandiaDecay::Nuclide *> traceSourceNuclides() const;
  
  /** Returns if this shielding is a trace source for the specified nuclide. */
  bool isTraceSourceForNuclide( const SandiaDecay::Nuclide *nuc ) const;
  
  /** Returns the current total activity for the specified trace source nuclide.
   
   Throws exception if not a trace source for the nuclide.
   */
  double traceSourceTotalActivity( const SandiaDecay::Nuclide *nuc ) const;
  
  /** Sets the specified nuclides activity, and causes correct signals to be emitted so SourceFitMode
   
   You always pass in total activity of trace source, which will then be converted to the appropriate display activity.
   
   Throws exception if not a trace source for the nuclide.
   */
  void setTraceSourceTotalActivity( const SandiaDecay::Nuclide *nuc, const double activity );
  
  
  /** Returns the current display activity for the specified trace source nuclide; e.g., might be total activity, or activity per cm^3, or
   per gram.
   
   Throws exception if not a trace source for the nuclide.
   
   \sa traceSourceType
   */
  double traceSourceDisplayActivity( const SandiaDecay::Nuclide *nuc ) const;
  
  
  /** Returns if the activity for this source should be fit for; this returns the same value as the model should return for this quantity. */
  bool fitTraceSourceActivity( const SandiaDecay::Nuclide *nuc ) const;
  
  
  /** Recalculates total activity, based on currently displayed activity and current geometry size.
   
   @returns the updated activity value.
   
   Throws exception if not a trace source for the nuclide.
   */
  double updateTotalTraceSourceActivityForGeometryChange( const SandiaDecay::Nuclide *nuc );
  
  
  /** Returns the current total activity for the specified trace source nuclide.
   
   Throws exception if not a trace source for the nuclide.
   */
  GammaInteractionCalc::TraceActivityType traceSourceType( const SandiaDecay::Nuclide *nuc ) const;
  
  /** Returns the relaxation length for the specified trace source.
   
   Throws exception if not a trace source for the nuclide, or TraceActivityType != TraceActivityType::ExponentialDistribution.
   */
  double relaxationLength( const SandiaDecay::Nuclide *nuc ) const;
  
  
  //serialize(...): saves state as a <Shielding />  node.
  void serialize( rapidxml::xml_node<char> *parent_node ) const;
  
  /** Returns the widget to the state specified by the XML node
   
   However,  `m_forFitting` must be the same in the XML as object this function
   is being called on (or else an exception is thrown).
    
   @param shielding_node The XML node saved by `serialize(...)`
   @param is_fixed_geom_det If the current DetectorPeakResponse is for fixed geometry or not.
          If false, then intrinsic and trace sources will not be allowed, which prevents some callback that will hit asserts.
   
   Throws exception on invalid XML or `m_forFitting` mismatch
   */
  void deSerialize( const rapidxml::xml_node<char> *shielding_node,
                   const bool is_fixed_geom_det );
  
  /** Encodes current tool state to app-url format.  Returned string is just the query portion of URL;
   so will look something like "V=1&G=S&D1=1.2cm", and it will not be url-encoded.
   
   TODO: Currently does not encode self-attenuating source, trace-source information, or "truth" values.
   */
  std::string encodeStateToUrl() const;

  void handleAppUrl( std::string query_str );
  
  ShieldingSourceFitCalc::ShieldingInfo toShieldingInfo() const;
  void fromShieldingInfo( const ShieldingSourceFitCalc::ShieldingInfo &info );
  
#if( INCLUDE_ANALYSIS_TEST_SUITE )
  boost::optional<double> truthThickness; //Spherical thickness, radial thickness (cylindrical), or rectagular width
  boost::optional<double> truthThicknessTolerance;
  boost::optional<double> truthAD;
  boost::optional<double> truthADTolerance;
  boost::optional<double> truthAN;
  boost::optional<double> truthANTolerance;
  
  boost::optional<double> truthThicknessD2;  //Cylinder length, or rectangular height
  boost::optional<double> truthThicknessD2Tolerance;
  boost::optional<double> truthThicknessD3;  //Rectangular depth
  boost::optional<double> truthThicknessD3Tolerance;
  
  std::map<const SandiaDecay::Element *,std::map<const SandiaDecay::Nuclide *,std::pair<double,double>>> truthFitMassFractions;
#endif

  
protected:
  void init();

  //emitRemoveSignal():  not only emits the remove() signal, but first unchecks
  //  all checkboxes in m_asSourceCBs, and emits the removingIsotopeAsSource()
  //  signal
  void emitRemoveSignal();

  //emitAddBeforeSignal(): emits the m_addShieldingBefore signal.
  void emitAddBeforeSignal();
  
  //emitAddAfterSignal(): emits the m_emitAddAfterSignal signal
  void emitAddAfterSignal();
  
  /** Checks #m_geometry is compatible with the desired geometry, and that it is not a generic material, and throws an exception if
   not compatible.
   GeometryType::CylinderEndOn and GeometryType::CylinderSideOn are treated as the same.
   GeometryType::NumGeometryType checks that this is a generic material
   */
  void checkIsCorrectCurrentGeometry( const GammaInteractionCalc::GeometryType wanted,
                                      const char *fcn_name ) const;
  
  void addTraceSource();
  void removeTraceSourceWidget( TraceSrcDisplay *src );
  void setTraceSourceMenuItemStatus();
  void handleTraceSourceNuclideChange( TraceSrcDisplay *src, const SandiaDecay::Nuclide *oldNuc );
  void handleTraceSourceWidgetAboutToBeRemoved( TraceSrcDisplay *nuc );
  void handleTraceSourceActivityChange( const SandiaDecay::Nuclide *nuc, const double activity );
  
  /** Returns the volume (in PhysicalUnits units) of this shielding.  This value does not include the volumes of sub-shieldings. */
  double shieldingVolume() const;
  
  /** Returns the mass (in PhysicalUnits units) of this shielding.  This value does not include the mass of sub-shieldings. */
  double shieldingMass() const;
  
  /** Returns the surface area for the surface that an in-situ exponential surface-contamination would be applicable to. */
  double inSituSurfaceArea() const;
  
  
  //This simply toggles the generic, and calls handleMaterialChange()
  void handleToggleGeneric();
  
  void displayInputsForCurrentGeometry();
  
  //handleMaterialChange(): handles when the user changes or modifies the
  //  current material.  If the material is changed, then possible source
  //  isotopes whic are no longer valid will be reomoved.
  //  Emits either the materialChanged() or materialModified() signals depending
  //  on what has been done.
  void handleMaterialChange();
  
  /** Handles when `m_materialEdit` emits the `changed()` or `enterPressed()` signals.
   
   Does a simple guess to correct for the material name not being valid if the user clicked on one
   of the popup suggestion (this first `m_materialEdit->changed()` signal will have only the
   partial text the user entered to get the suggestions).
   */
  void handleUserChangedMaterialName();
  
  /** Updates #m_currentMaterial if `m_materialEdit->text() != m_currentMaterialDescrip` */
  void updateMaterialFromUserInputTxt();
  
  /** Handles updating #m_prevState, and emitting userChangedStateSignal(). */
  void handleUserChangeForUndoRedo();
  void handleUserChangeForUndoRedoWorker( const bool emit_change );

  /** Removes uncertainty number from a distance edit text; useful when user changes text of a WLineEdit that is also fit for, so may
   have an uncertainty value in it
   */
  void removeUncertFromDistanceEdit( Wt::WLineEdit *edit );
  
  /** Changes the mass fraction of 'nuclide'.
   
   @param fraction The fraction of the element, the specified nuclide, accounts for.
          Value should be between 0.0, and 1.0.
   @param nuclide The nuclide to set the mass-fraction for.
  
   
  This function tries to adjust the other isotope fractions in a way that is kinda intuitive to
  the user, while enforcing consistency.
   
   The overall mass-fraction for the respective element stays the same, this function just
   changes the isotopic composition of the element.  (i.e. the mass-fraction of the element
   will stay the same).
   */
  void handleIsotopicChange( const float fraction, 
                            const SandiaDecay::Nuclide * const nuclide );

  /** Called when fit mass fraction checkbox is changed for a nuclide.
   If it is the "other nucs" checkbox, then nuclide will be nullptr, and element will say which element this is for.
   */
  void handleFitMassFractionChanged( const bool fit, const SandiaDecay::Nuclide * const nuclide,
                                    const SandiaDecay::Element *el );
  
  //modelNuclideAdded(...): adds a checkbox for Nuclide and connects
  //  appropriate signals, and updates trace sources
  void modelNuclideAdded( const SandiaDecay::Nuclide *iso );
  
  /** Returns the mass-fraction of the nuclide, of the entire material.
 
   If the material specifies the isotopic fractions, that value (the mass fraction relative to the element) will be returned.
   If the material specifies just the element of the nuclide, and the nuclide is a naturally present nuclide, then
   the function returns the elements mass fraction.
   If the element is specified, but the nuclide is not a naturally present nuclide, returns the element mass fraction time
   1 over the total number of nuclides in element (this is a fairly arbitrary behavior).
   
   Throws exception if material does not contain the nuclide, or the element of the nuclide.
   Also throws exception if nuclide or material are nullptr.
   */
  static double nuclidesMassFractionInElementOfMaterial( const SandiaDecay::Nuclide * const iso,
                                       const std::shared_ptr<const Material> &mat );

  static double elementsMassFractionInMaterial( const SandiaDecay::Element * const el,
                                               const std::shared_ptr<const Material> &mat );
  
  //isotopeCheckedCallback(...): emits the addingIsotopeAsSource() signal
  void isotopeCheckedCallback( const SandiaDecay::Nuclide *iso );

  //isotopeUnCheckedCallback(...): emits the removingIsotopeAsSource() signal
  void isotopeUnCheckedCallback( const SandiaDecay::Nuclide *iso );

  /** Checks if only a single isotope is fitting mass-fraction, then will also make sure "other" non-src fraction will be fit.
   Will also hide/show "other" non-src fraction, depending if any nuclides are being used.
   */
  void checkAndUpdateMassFractionCanFit();

  /** Returns the trace source widget for the specfied nuclide, or nullptr if nuclide is not a trace source. */
  const TraceSrcDisplay *traceSourceWidgetForNuclide( const SandiaDecay::Nuclide *nuc ) const;
  
  /** Non-const version of #traceSourceWidgetForNuclide */
  TraceSrcDisplay *traceSourceWidgetForNuclide( const SandiaDecay::Nuclide *nuc );
  
  /** Serializes the widgets state to a #rapidxml::xml_document object.
   
   Intended to update #m_prevState during the undo/redo step stuff.
   */
  std::shared_ptr<const std::string> getStateAsXml() const;
  
protected:
  /** Pointer to the ShieldingSourceDisplay this ShieldingSelect belongs to - if it belongs to this tool, otherwise will be nullptr. */
  const ShieldingSourceDisplay *m_shieldSrcDisp;
  
  Wt::WImage* m_toggleImage;
  const bool m_forFitting;
  MaterialDB *m_materialDB;
  SourceFitModel *m_sourceModel;
  
  GammaInteractionCalc::GeometryType m_geometry;

  Wt::WSuggestionPopup *m_materialSuggest;  //Not owned by this object
  
  Wt::WLineEdit *m_materialEdit;
  bool m_isGenericMaterial;
  Wt::WText *m_materialSummary;
  Wt::WPushButton *m_closeIcon;
  Wt::WPushButton *m_addIcon;
  PopupDivMenuItem *m_addTraceSourceItem;
  
  /** If true, then `m_addTraceSourceItem` should be hidden, spherical geometry, and there should be no trace sources,
   and no self-attenuating sources.
   */
  bool m_fixedGeometry;
  
  /** Stack to hold the dimension edits for spherical, cylindrical, and rectangular geometries, as well as m_genericMaterialDiv.
   First entry is m_genericMaterialDiv, then spherical, then cylindrical, then rectangular
   */
  Wt::WStackedWidget *m_dimensionsStack;
  
  Wt::WContainerWidget *m_genericDiv;
  NativeFloatSpinBox *m_arealDensityEdit;
  Wt::WCheckBox *m_fitArealDensityCB;
  NativeFloatSpinBox *m_atomicNumberEdit;
  Wt::WCheckBox *m_fitAtomicNumberCB;
  
  
  Wt::WContainerWidget *m_sphericalDiv;
  Wt::WLineEdit *m_thicknessEdit;
  Wt::WCheckBox *m_fitThicknessCB;
  
  
  Wt::WContainerWidget *m_cylindricalDiv;
  Wt::WLineEdit *m_cylRadiusEdit;
  Wt::WCheckBox *m_fitCylRadiusCB;
  Wt::WLineEdit *m_cylLengthEdit;
  Wt::WCheckBox *m_fitCylLengthCB;
  
  
  Wt::WContainerWidget *m_rectangularDiv;
  Wt::WLineEdit *m_rectWidthEdit;
  Wt::WCheckBox *m_fitRectWidthCB;
  Wt::WLineEdit *m_rectHeightEdit;
  Wt::WCheckBox *m_fitRectHeightCB;
  Wt::WLineEdit *m_rectDepthEdit;
  Wt::WCheckBox *m_fitRectDepthCB;
  
  
  //Wt::WCheckBox *m_fitMassFrac;  //is nullptr if !m_forFitting

  Wt::WContainerWidget *m_asSourceCBs;
  
  typedef std::map<const SandiaDecay::Element *,Wt::WContainerWidget *> ElementToNuclideMap;
  ElementToNuclideMap m_sourceIsotopes;

  Wt::WContainerWidget *m_traceSources;
  
  
  std::string m_currentMaterialDescrip;
  
  // TODO: should m_currentMaterial be made a const pointer?
  std::shared_ptr<const Material> m_currentMaterial;

  Wt::Signal<ShieldingSelect *> m_removeSignal;
  Wt::Signal<ShieldingSelect *> m_addShieldingBefore;
  Wt::Signal<ShieldingSelect *> m_addShieldingAfter;
  Wt::Signal<ShieldingSelect *> m_materialModifiedSignal;
  Wt::Signal<ShieldingSelect *> m_materialChangedSignal;
  Wt::Signal<const SandiaDecay::Nuclide *,ShieldingSourceFitCalc::ModelSourceType> m_addingIsotopeAsSource;
  Wt::Signal<const SandiaDecay::Nuclide *,ShieldingSourceFitCalc::ModelSourceType> m_removingIsotopeAsSource;
  Wt::Signal<ShieldingSelect *,const SandiaDecay::Nuclide *>
                                            m_activityFromVolumeNeedUpdating;
  
  //To add undo/redo support, I think we need to add:
  Wt::Signal<ShieldingSelect *,
             std::shared_ptr<const std::string>,
             std::shared_ptr<const std::string>> m_userChangedStateSignal;
  std::shared_ptr<const std::string> m_prevState;

  friend class TraceSrcDisplay;
  friend class ShieldingSourceDisplay;
};//struct ShieldingSelect

#endif //ShieldingSelect_h
