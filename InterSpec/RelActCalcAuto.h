#ifndef RelActCalcAuto_h
#define RelActCalcAuto_h
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
#include <atomic>
#include <string>
#include <memory>
#include <vector>
#include <ostream>
#include <variant>
#include <optional>

#include "InterSpec/PeakDef.h" //for PeakContinuum::OffsetType and PeakDef::SkewType
#include "InterSpec/RelActCalc.h"

// Forward declarations
class DetectorPeakResponse;
struct Material;
class MaterialDB;

namespace rapidxml
{
  template<class Ch> class xml_node;
}

namespace SpecUtils
{
  class Measurement;
  struct EnergyCalibration;
}//namespace SpecUtils

namespace SandiaDecay
{
  struct Nuclide;
}

namespace RelActCalc
{
  enum class PuCorrMethod : int;
  enum class RelEffEqnForm : int;
  struct PhysicalModelShieldInput;
  template<typename T> struct Pu242ByCorrelationOutput;
}//namespace RelActCalc

namespace RelActCalcAutoImp
{
  struct RelActAutoCostFcn;
}

/*
 Further things to consider:
 - For the relative efficiency chart, currently all markers for peaks info fall exactly on the
    RelEff line. What we should do is, after finding the solution, cluster the peaks together
    (similar to Act/Shield fit), and then make a template for each clustered peak (sum all peaks into
    a channel count array, with separate array for each clustered peak), with some threshold (like
    99% detection threshold) that if the clustered peaks area is less than, it just gets subtracted
    from the data, unless its the only peak for that ROI.
    Then fit the templates to data, and use the fit values to make the markers on the RelEff plot.
    I think the plotting JS can handle showing info for multiple nuclides in a single plots
    mouse-over.  We will need to create new variable that is analogous to
    #RelActAutoSolution::m_fit_peaks, but maps multiple source nuclides to a single peak.
 - Sum peaks - not totally sure how to handle - can there just be one parameter that controls
    strength of peak random summing, and then another for cascade?  For cascades we dont have
    gamma-x-ray coincidences (or is it the case we can just assume constant fraction between any
    gamma and any nuclide?), so we might need two parameters, one for gamams, one for gamma-neutrons
    and I guess maybe xray-xray?
 */

namespace RelActCalcAuto
{

int run_test();

/** A typdef for passing either Nuclide, Element, or Reaction to functions. 
 
 Note that it also includes using a monostate, because we would like to semi-enforse that if a pointer is
 used, then it is valid; if we didnt include the monostate, then default construction of the variant would
 use a nullptr for the first pointer type.
*/
typedef std::variant<std::monostate,const SandiaDecay::Nuclide *, const SandiaDecay::Element *, const ReactionGamma::Reaction *> SrcVariant;

/** Returns the #SandiaDecay::Nuclide from the #SrcVariant, returning nullptr if not a #SandiaDecay::Nuclide. */
const SandiaDecay::Nuclide *nuclide( const SrcVariant &src );

/** Returns the #SandiaDecay::Element from the #SrcVariant, returning nullptr if not a #SandiaDecay::Element. */
const SandiaDecay::Element *element( const SrcVariant &src );

/** Returns the #ReactionGamma::Reaction from the #SrcVariant, returning nullptr if not a #ReactionGamma::Reaction. */
const ReactionGamma::Reaction *reaction( const SrcVariant &src );

/** Returns true if the pointer in the #SrcVariant is a nullptr. */
bool is_null( const SrcVariant &src );

/** Returns a string name for the #SrcVariant. 
   Either SandiaDecay::Nuclide::symbol, SandiaDecay::Element::symbol, or ReactionGamma::Reaction::name(). 
   
   Throws exception if monostate applicable pointer is nullptr.
*/
std::string to_name( const SrcVariant &src );

/** Returns a #SrcVariant from a string name.
 * 
   If a Nuclide, Element, or Reaction can not be found for name, return monostate.
*/
SrcVariant source_from_string( const std::string &name );


/** Struct to specify an energy range to consider for doing relative-efficiency/activity calc.
 */
struct RoiRange
{
  double lower_energy = -1.0;
  double upper_energy = -1.0;
  
  /** The continuum type to use.
   
   TODO: Allow setting to #PeakContinuum::OffsetType::External to auto-choose this.
   */
  PeakContinuum::OffsetType continuum_type = PeakContinuum::OffsetType::Quadratic;

  /** Specifies how the energy range limits should be interpreted during fitting.

   This controls whether the ROI boundaries are strictly enforced, can expand to accommodate
   peaks at the edges, or can be broken into smaller ranges based on peak locations.
   */
  enum class RangeLimitsType : int
  {
    /** The lower and upper energies are strictly enforced - the ROI will not be
     adjusted or broken up based on peak locations.
     */
    Fixed = 0,

    /** The lower and upper energy values can be expanded to accommodate the FWHM of peaks
     at or near the ROI edges. Not well tested yet.
     */
    CanExpandForFwhm = 1,

    /** The ROI can be broken up into multiple smaller ROIs based on expected significant peaks.
     The range may be subdivided to optimize fitting.
     */
    CanBeBrokenUp = 2
  };//enum class RangeLimitsType

  /** Specifies how the ROI limits should be interpreted. */
  RangeLimitsType range_limits_type = RangeLimitsType::CanBeBrokenUp;

  /** Returns string representation of the RangeLimitsType.
   String returned is a static string, so do not delete it.
   */
  static const char *to_str( const RangeLimitsType type );

  /** Converts from the string representation of RangeLimitsType to enumerated value.

   Throws exception if invalid string (i.e., any string not returned by to_str(RangeLimitsType) ).
   */
  static RangeLimitsType range_limits_type_from_str( const char *str );

  static const int sm_xmlSerializationVersion = 0;
  void toXml( ::rapidxml::xml_node<char> *parent ) const;
  void fromXml( const ::rapidxml::xml_node<char> *parent );

#if( PERFORM_DEVELOPER_CHECKS )
  static void equalEnough( const RoiRange &lhs, const RoiRange &rhs );
#endif
};//struct RoiRange


/** Struct to specify a nuclide to use for doing relative-efficiency/activity calc.
 
 TODO: min/max rel act to be a constraint
 */
struct NucInputInfo
{
  SrcVariant source;
  
  /** Age in units of PhysicalUnits (i.e., 1.0 == second), for the nuclide.
   
   If the age is being fit, this will be used as the initial age to start with.
   
   Must not be negative if this is a nuclide, but must be zero or negative for element or reaction.
   */
  double age = -1.0;
  
  /** If nuclide age should be fit.  Must be false for x-ray or reaction. */
  bool fit_age = false;
  
  /** Minimum age to fit; if specified, must be zero or greater. 
   
   In units of PhysicalUnits (i.e., 1.0 == second).
   
   May not be specified if `fit_age` is false, or a x-ray or reaction.
  */
  std::optional<double> fit_age_min;
  
  /** Maximum age to fit; if specified, must be greater than #fit_age_min (if specified). 
  
   In units of PhysicalUnits (i.e., 1.0 == second).

   May not be specified if `fit_age` is false, or a x-ray or reaction.
  */
  std::optional<double> fit_age_max;

  /** Minimum relative activity to fit; if specified, must be zero or greater. 
   
   May not be specified if this nuclide is subject to a ActRatioConstraint (as the controlled nuclide).

   Setting #min_rel_act equal to #max_rel_act will cause the fit to be constrained to that value (
   and if #starting_rel_act is set, it must be exactly equal to that value).
  */
  std::optional<double> min_rel_act;
  
  /** Maximum relative activity to fit; if specified, must be greater than #min_rel_act (if specified). 
   
   May not be specified if this nuclide is subject to a ActRatioConstraint (as the controlled nuclide).

   Setting #min_rel_act equal to #max_rel_act will cause the fit to be constrained to that value (
   and if #starting_rel_act is set, it must be exactly equal to that value).
  */
  std::optional<double> max_rel_act;

  /** Starting relative activity; if specified, must be zero or greater. 
   
   May not be specified if this nuclide is subject to a ActRatioConstraint (as the controlled nuclide).
  */
  std::optional<double> starting_rel_act;

  /** Energy corresponding to SandiaDecay::EnergyRatePair::energy, or equivalent for x-ray or reactions. */
  std::vector<double> gammas_to_exclude;
  
  std::string peak_color_css;
  
  const std::string name() const;

  bool operator==( const NucInputInfo &rhs ) const;

  static const int sm_xmlSerializationVersion = 0;
  void toXml( ::rapidxml::xml_node<char> *parent ) const;
  void fromXml( const ::rapidxml::xml_node<char> *parent );

#if( PERFORM_DEVELOPER_CHECKS )
  static void equalEnough( const NucInputInfo &lhs, const NucInputInfo &rhs );
#endif
};//struct NucInputInfo


/** A peak at a specific energy, that has a free-floating amplitude, and may also have a floating
 FWHM.
 */
struct FloatingPeak
{
  /** Energy (in keV) of the peak.
   
   Note that if energy calibration is being fit for, this energy will not have that correction
   applied when fitting things; this is because these peaks are nominally expected to be added
   by a user who gets the energy based on looking at the spectrum.
   */
  double energy = -1.0;
  
  /** If true, the FWHM of the peak will be allowed to freely vary from 0.25 to 4.0 times (both
   numbers arbitrarily chosen) the FWHM that would be predicted by the functional FWHM
   for the energy.
   */
  bool release_fwhm = false;
  
  /** Whether the energy correction should be applied while performing the fit for the answer.
 
   That is, if #FloatingPeak::energy represents a known true gamma energy, set this value to true;
   if this peak represents an observed peak in the spectrum (maybe from unknown origins), set
   this value to false.
   
   Only has an effect when the energy calibration is being fit for as well.
   */
  bool apply_energy_cal_correction = true;
  
  // TODO: should maybe have a max-width in FloatingPeak
  
  static const int sm_xmlSerializationVersion = 0;
  void toXml( ::rapidxml::xml_node<char> *parent ) const;
  void fromXml( const ::rapidxml::xml_node<char> *parent );

#if( PERFORM_DEVELOPER_CHECKS )
  static void equalEnough( const FloatingPeak &lhs, const FloatingPeak &rhs );
#endif
};//struct FloatingPeak


/** The FWHM functional form to use.
 
 See #DetectorPeakResponse::ResolutionFnctForm for details.
 
 
 See #DetectorPeakResponse::Polynomial for details
 
 FWHM = sqrt( Sum_i{A_i*pow(x/1000,i)} );
 */
enum class FwhmForm : int
{
  /** See #DetectorPeakResponse::peakResolutionFWHM implementation for details.
   Will use three parameters to fit FWHM as a function of energy.
   */
  Gadras,
  
  /** The equation `FWHM = sqrt( A_0 + A_1*Energy + A_2/Energy )`.
   This is what FRAM uses.
   */
  SqrtEnergyPlusInverse,
  
  /** The equation `FWHM = A_0 + A_1*sqrt(Energy)`
   This is what some commercial applications use - but should probably only be used for cases of interfacing to the.
   */
  ConstantPlusSqrtEnergy,
  
  /** The linear polynomial: e.g. FWHM = sqrt(A_0 + A_1*1*0.001*Energy)
   The "2" refers to how many parameters are being fit
   */
  Polynomial_2,
  
  /** The quadratic polynomial: e.g. FWHM = sqrt(A_0 + A_1*0.001*Energy + A_2*0.000001*Energy*Energy */
  Polynomial_3,
  Polynomial_4,
  Polynomial_5,
  Polynomial_6,
  
  /** Berstein polynomial with 2 coefficients (order 1) */
  Berstein_2,
  
  /** Berstein polynomial with 3 coefficients (order 2) */
  Berstein_3,
  
  /** Berstein polynomial with 4 coefficients (order 3) */
  Berstein_4,
  
  /** Berstein polynomial with 5 coefficients (order 4) */
  Berstein_5,
  
  /** Berstein polynomial with 6 coefficients (order 5) */
  Berstein_6,

  /** Do not fit the FWHM equation - use the FWHM from the detector efficiency function.
   
   #Options::fwhm_estimation_method must be set to FwhmEstimationMethod::FixedToDetectorEfficiency,
   and you must provide a DetectorPeakResponse, with FWHM information, to `RelActCalcAuto::solve(...)`.
  */
  NotApplicable
};//enum FwhmForm


/** Returns string representation of the #FwhmForm.
 String returned is a static string, so do not delete it.
 */
const char *to_str( const FwhmForm form );

/** Converts from the string representation of #FwhmForm to enumerated value.
 
 Throws exception if invalid string (i.e., any string not returned by #to_str(FwhmForm) ).
 */
FwhmForm fwhm_form_from_str( const char *str );


/** How the FWHM of peaks should be determined.
 
 We could add in ability to specify parameters and the ranges they can vary over, but for now
 we'll keep it simple - it isnt clear complicating things will actually be useful to users.
 */
enum class FwhmEstimationMethod : int
{
  /** Will use the FWHM of the DetectorPeakResponse, if available, or if not will fit the
   peaks in the spectrum, and use the derived FWHM to start with, refining it in the non-linear
   fit for the relative activities.
   */
  StartFromDetEffOrPeaksInSpectrum,

  /** Use the FWHM equation fit from all peaks in the spectrum, 
    but further refine this during the non-linear fit for the relative activities. 
  */
  StartingFromAllPeaksInSpectrum,
  
  /** Use the FWHM equation fit from all peaks in the spectrum, 
    and do not further refine this during the non-linear fit for the relative activities. 
   */
  FixedToAllPeaksInSpectrum,

  /** Use the detector efficiency function to estimate the FWHM,
    but further refine this during the non-linear fit for the relative activities. 

    Will throw exception if detector efficiency function does not have a FWHM
   */
  StartingFromDetectorEfficiency,
  
  /** Use the detector efficiency function to estimate the FWHM,
    and do not further refine this during the non-linear fit for the relative activities. 

    Options::fwhm_form
   */
  FixedToDetectorEfficiency
};//enum FwhmEstimationMethod


/** Returns string representation of the #FwhmEstimationMethod.
 String returned is a static string, so do not delete it.
 */
const char *to_str( const FwhmEstimationMethod form );


/** Converts from the string representation of #FwhmEstimationMethod to enumerated value.
 
 Throws exception if invalid string (i.e., any string not returned by #to_str(FwhmEstimationMethod) ).
 */
FwhmEstimationMethod fwhm_estimation_method_from_str( const char *str );


/** Type of energy calibration fitting to perform.

 Controls whether and how energy calibration is adjusted during the fit.
 */
enum class EnergyCalFitType : int
{
  /** Do not fit any energy calibration corrections. */
  NoFit,

  /** Fit linear energy calibration corrections (offset and gain adjustments). */
  LinearFit,

  /** Fit non-linear energy calibration corrections using deviation pairs.

   Requires at least 3 ROIs; if less than 3 ROIs, will use linear fit instead.
   Each ROI gets an anchor point at its largest peak (or midpoint 
   if no peaks).  The first and last ROIs have their deviation offsets fixed to zero.  Middle
   ROIs have their deviation offsets fitted.  The deviation pairs are interpolated using a
   cubic spline to provide smooth energy corrections.
   */
  NonLinearFit
};//enum class EnergyCalFitType


/** Returns string representation of the #EnergyCalFitType.
 String returned is a static string, so do not delete it.
 */
const char *to_str( const EnergyCalFitType type );


/** Converts from the string representation of #EnergyCalFitType to enumerated value.

 Throws exception if invalid string (i.e., any string not returned by #to_str(EnergyCalFitType) ).
 */
EnergyCalFitType energy_cal_fit_type_from_str( const char *str );


size_t num_parameters( const FwhmForm eqn_form );


struct RelEffCurveInput
{
  RelEffCurveInput();

  /** The name of the relative efficiency curve. 
   
   This value is not used within calculations, but is used for the reporting and 
   display of results.
  */
  std::string name;

  /** The nuclides that apply to this relative efficiency curve.
   * You may specify the same nuclide for multiple different RelEffCurveInput inputs,
   * but becareful of degeneracies in the problem.
   */
  std::vector<NucInputInfo> nuclides;
  
  /** If true, all nuclides of a given element will be constrained to be the same age.
   
   If true, all #NucInputInfo::age of nuclides of the same element must be the same, or an exception
   will be thrown (even if age is being fit, this is .
   */
  bool nucs_of_el_same_age;
  
  /** The functional form of the relative efficiency equation to use.
   */
  RelActCalc::RelEffEqnForm rel_eff_eqn_type;
  
  /** The number of energy dependent terms to use for the relative efficiency equation.  I.e., the
   number of terms fit for will be one more than this value.
   
   Not used for `RelEffEqnForm::FramPhysicalModel`; see `phys_model_self_atten_an`
   and `phys_model_ext_atten_an`
   */
  size_t rel_eff_eqn_order;

  /** Only used for `RelEffEqnForm::FramPhysicalModel` - shielding definitions for self-attenuating sources.
   
   May be nullptr if physical model.
   */
  std::shared_ptr<const RelActCalc::PhysicalModelShieldInput> phys_model_self_atten;
  
  /** Only used for `RelEffEqnForm::FramPhysicalModel` - shielding definitions for external attenuators. */
  std::vector<std::shared_ptr<const RelActCalc::PhysicalModelShieldInput>> phys_model_external_atten;

  /** The indices of other Physical Model Rel. Eff. curves, whose self-attenuation and external attuators act as external shieldings for
   this Rel. Eff. curve.  If any indices are provided, a number of requirements must be met.
   - The indices provided must be for Physical Model defined rel. eff. curves, and the indix of this curve may not be provided.
   - The `Options::same_external_shielding_for_all_rel_eff_curves` option must be false.
   - The curve(s) shielding this curve may not also be shielded by this curve.
   */
  std::set<size_t> shielded_by_other_phys_model_curve_shieldings;

  /** If true, fit the modified Hoerl equation form for the physical model.
   * If false, do not fit the modified Hoerl equation form (its value will be constant value of 1.0).
   *
   * Ignored if not using `RelActCalc::RelEffEqnForm::FramPhysicalModel`.
  */
  bool phys_model_use_hoerl = true;

  /** The method to use for Pu-242 mass-enrichment estimation (all methods use correlation to other
   Pu isotopics).
   
   When specified #RelActAutoSolution::m_corrected_pu will be filled-out (unless an error occurred
   while doing the correction - which shouldnt really ever happen).
   
   Defaults to #PuCorrMethod::NotApplicable; if any other value is specified, and sufficient
   plutonium isotopes for that method are not in the problem, then finding the solution will
   fail.
   */
  RelActCalc::PuCorrMethod pu242_correlation_method;

  /** A constraint on the activity of a nuclide, relative to another nuclide. 
   
   TODO: move controlled/constrained nuclide to be `SrcVariant` (i.e., `std::variant<Nuc,El,Rctn>`)
   TODO: currently this constraint only applies to a single Rel. Eff. curve; see TODO note in
        MassFractionConstraint, and do similar for this.
  */
  struct ActRatioConstraint
  {
    /** The source whose relative activity is varied in the fit. */
    SrcVariant controlling_source;

    /** The source that is fixed to the `controlling_source`. */
    SrcVariant constrained_source;

    /** The ratio of the activity of the `constrained_source` to the `controlling_source`. 
     
     e.g., `act(constrained_source) = constrained_to_controlled_activity_ratio * act(controlling_source)`
    */
    double constrained_to_controlled_activity_ratio = 0.0;
    
    static ActRatioConstraint from_mass_ratio( const SandiaDecay::Nuclide *constrained, 
                                              const SandiaDecay::Nuclide *controlling, 
                                            const double mass_ratio_constrained_to_controlling );

    // TODO: add ability to specify an uncertainty.
    //double constrained_to_controlled_activity_ratio_uncertainty = 0.0;

  static const int sm_xmlSerializationVersion = 0;
  void toXml( ::rapidxml::xml_node<char> *parent ) const;
  void fromXml( const ::rapidxml::xml_node<char> *constraint_node );

#if( PERFORM_DEVELOPER_CHECKS )
  static void equalEnough( const ActRatioConstraint &lhs, const ActRatioConstraint &rhs );
#endif
  };//struct ActRatioConstraint

  /** Constraints on nuclide activities. 
    
    All input nuclides must be valid, and be found in `nuclides`.

    The constrained nuclides must not be cyclical; e.g., if Nuc1 is constrained to Nuc2, 
    then there can not be a constraint that Nuc2 is constrained to Nuc1.
  */
  std::vector<ActRatioConstraint> act_ratio_constraints;


  /** A constraint on the mass fraction of an nuclide within an element. 
   

   TODO: currently this constraint only applies to a single Rel. Eff. curve; if you have multiple
          curves with same constraint/nuclide, then each constraint will be applied seperately.
          We should add in a `set<size_t>` that specifies the indexes of all the Rel. Eff. curves
          that this constraint applies to (if only single curve, then allow it to be empty), and 
          then move the constraint up a level to the `Options` struct.
  */
  struct MassFractionConstraint
  {
    const SandiaDecay::Nuclide *nuclide = nullptr;
    double lower_mass_fraction = 0.0;
    double upper_mass_fraction = 0.0;

    static const int sm_xmlSerializationVersion = 0;
    void toXml( ::rapidxml::xml_node<char> *parent ) const;
    void fromXml( const ::rapidxml::xml_node<char> *parent );

#if( PERFORM_DEVELOPER_CHECKS )
  static void equalEnough( const MassFractionConstraint &lhs, const MassFractionConstraint &rhs );
#endif
  };//struct MassFractionConstraint

  /** Constraints on nuclide mass fractions. 
    
    All input nuclides must be valid, and be found in `nuclides`.

    The constrained nuclides must not be a controlled nuclide by ActRatioConstraint,
    or have the min or max relative activity specified.

    If the constrained nuclide controlls the activity of another nuclide, then that other nuclide
    must be of a different element.
    (this could be changed in the future, but for the initial implementation, we'll enforce this)

    The lower mass fractions, for a particular element, must sum to less than 1.0.

    There must be at least one nuclide for the element that does not have a mass fraction constraint.
    (this _could_ be relaxed in the future, but for the initial implementation, this is required)
  */
  std::vector<MassFractionConstraint> mass_fraction_constraints;


  /** Checks that the nuclide constraints are valid.

   Checks for cyclical constraints, and that all constrained nuclides are found in #nuclides.
   
   Throws an exception if they are not valid.
   */
  void check_nuclide_constraints() const; 

  
  static const int sm_xmlSerializationVersion = 0;

  /** Puts this object to XML
   * 
   * @param parent The parent node to add this object to.
   * @param into_parent_node If true, this objects fields will be added to the parent node, otherwise a new node will be created.
   * @returns The node the fields were written into (so the added node, or the parent node passed in).
   */
  rapidxml::xml_node<char> *toXml( ::rapidxml::xml_node<char> *parent ) const;
  void fromXml( const ::rapidxml::xml_node<char> *parent, MaterialDB *materialDB );

#if( PERFORM_DEVELOPER_CHECKS )
  static void equalEnough( const RelEffCurveInput &lhs, const RelEffCurveInput &rhs );
#endif
};//struct RelEffCurveInput



struct Options
{
  Options();

  /** Type of energy calibration fitting to perform.

   Controls whether and how energy calibration is adjusted during the fit.
   */
  EnergyCalFitType energy_cal_type;

  /** The functional form of the FWHM equation to use.  The coefficients of this equation will
   be fit for across all energy ranges.
   */
  FwhmForm fwhm_form;
  
  /** How the FWHM of peaks should be determined. */
  FwhmEstimationMethod fwhm_estimation_method;

  /** Optional title of the spectrum; used as title in HTML and text summaries. */
  std::string spectrum_title;
  
  
  /** Peak skew to apply to the entire spectrum.
   
   Under development: currently, if total energy range being fit is less than 100 keV, then all peaks will share the same skew.
   Otherwise a linear energy dependence will be assumed, where the fitting parameters will be for the spectrums lower
   energy, and the spectrums upper energy; not all skew parameters are allowed to vary with energy; e.g., the Crystal Ball
   power law is not allowed to have an energy dependence (see `PeakDef::is_energy_dependent(...)`).
   */
  PeakDef::SkewType skew_type;

  /** Whether to use Lorentzian (Voigt) peak shapes for x-ray peaks.
   *
   * When true, x-ray peaks will use VoigtPlusBortel skew type with the
   * Lorentzian HWHM set to the natural x-ray linewidth (plus thermal/recoil
   * Doppler broadening for decay x-rays).
   *
   * Only compatible with skew_type == NoSkew or skew_type == GaussPlusBortel.
   * If set to true with an incompatible skew_type, then the problem will fail to setup (an exception
   * in `check_same_hoerl_and_external_shielding_specifications()` will be thrown)
   */
  bool lorentzian_xrays;

  /** An additional uncertainty applied to each roughly independent peak.
   * 
   * Contributing gammas are clustered into roughly independent peaks based on their energy and the 
   * initially estimated FWHM from the spectrum.  Then for each clustered energy range, the gammas 
   * within that range are allowed to vary in amplitude, with the deviation from nominally predicted
   * amplitude punished according to this uncertanty.
   * 
   * Its not perfect, but its something.
   */
  double additional_br_uncert;
  

  /** The relative efficiency curves to use.
   */
  std::vector<RelActCalcAuto::RelEffCurveInput> rel_eff_curves;


  std::vector<RelActCalcAuto::RoiRange> rois;


  std::vector<RelActCalcAuto::FloatingPeak> floating_peaks;


  /** If true, use the same Hoerl equation form for all relative efficiency curves.
   *
   * If false, each relative efficiency curve will have its own Hoerl equation fit (if applicable).
   * 
   * if true, and you do not have multiple physical models with `RelEffCurveInput::phys_model_use_hoerl`
   * set to true, then an exception will be thrown.
   */
  bool same_hoerl_for_all_rel_eff_curves;

  /** If true, use the same external shielding for all relative efficiency curves.
   *
   * If true, you must have multiple physical models defined, each with the same 
   * external shielding defined for each relative efficiency curve - or an exception 
   * will be thrown.
   * 
   * If the initial starting solution for the AutoRelEff curve is to be determined using
   * the ManualRelEff method, then the first RelEffCurveInput estimate will be what is 
   * used.
   */
  bool same_external_shielding_for_all_rel_eff_curves;

  /** If using the same Hoerl function, or external shielding for all relative efficiency curves,
   * this will check that the specifications are consistent.
   *
   * Throws an exception if they are not consistent.
   */
  void check_same_hoerl_and_external_shielding_specifications() const;

  /** Version history:
   - 20250117: incremented to 1 to handle FramPhysicalModel; if not this model, will still write version 0.
   - 20250130: incremented to 2 to handle multiple rel eff curves; can read backward compatible, but not write.
   */
  static const int sm_xmlSerializationVersion = 2;
  rapidxml::xml_node<char> *toXml( ::rapidxml::xml_node<char> *parent ) const;
  
  /** Sets the member variables from an XML element created by `toXml(...)`.
   @param parent An XML element with name "Options".
   @param materialDB The material database to use to retrieve a material from for the #PhysicalModelShieldInput;
          for other equation types, or for AN/AD defined shields, this isnt used/required.
          `same_hoerl_for_all_rel_eff_curves` and `same_external_shielding_for_all_rel_eff_curves`
          were also added, and are optional.
   */
  void fromXml( const ::rapidxml::xml_node<char> *parent, MaterialDB *materialDB );

#if( PERFORM_DEVELOPER_CHECKS )
  static void equalEnough( const Options &lhs, const Options &rhs );
#endif
};//struct Options

/** Struct to a sources fit relative activity.*/
struct NuclideRelAct
{
  SrcVariant source;
  
  std::string name() const;

  double age = -1.0;
  double age_uncertainty = -1.0;
  bool age_was_fit = false;
  
  double rel_activity = -1.0;
  double rel_activity_uncertainty = -1.0;
  
  /** The energy and its respective number of gammas per decay for this nuclide, at its age. 
   Note: this has not had the branching ratio uncertainty correction applied (if they were fit for).
  */
  std::vector<std::pair<double,double>> gamma_energy_br;
};//struct NuclideRelAc


struct FloatingPeakResult
{
  double energy;
  double amplitude;
  double amplitude_uncert;
  double fwhm;
  double fwhm_uncert;
};//struct FloatingPeakResult


/** This struct hold the `RelActAutoGui` state, but is defined here because we may want to use the XML
 from the GUI for batch processing, or whatever.
*/
struct RelActAutoGuiState
{
  RelActAutoGuiState();

  std::string note;

  /** Description field for this configuration. Not currently exposed in GUI. */
  std::string description;

  RelActCalcAuto::Options options;
  
    
  bool background_subtract;
  
  bool show_ref_lines;
  double lower_display_energy;
  double upper_display_energy;
  
  /** Returns XML node added; i.e., will have name "RelActCalcAuto" */
  ::rapidxml::xml_node<char> *serialize( ::rapidxml::xml_node<char> *parent ) const;
  void deSerialize( const rapidxml::xml_node<char> *base_node, MaterialDB *materialDb );

#if( PERFORM_DEVELOPER_CHECKS )
  static void equalEnough( const RelActAutoGuiState &lhs, const RelActAutoGuiState &rhs );
#endif
};//struct RelActAutoGuiState
  
  
struct RelActAutoSolution
{
  struct ObsEff; //forward decleration
  
  RelActAutoSolution();
  
  std::ostream &print_summary( std::ostream &strm ) const;
  
  void print_html_report( std::ostream &strm ) const;
  
  /** Prints the txt version of relative eff eqn. */
  std::string rel_eff_txt( const bool html_format, const size_t rel_eff_index ) const;
  
  /** Returns the fractional amount of an element (by mass), the nuclide composes, as well as (hopefully) the uncertainty.
   
   @param nuclide The nuclide of interest.
   @param rel_eff_index The relative efficiency index for the nuclide
   @returns The fraction of mass (so between 0.0 and 1.0), the input nuclide is responsible for
            in the solution, as well as the uncertainty, if the covariance matrix is defined.
            Ex., if Eu152 is passed in, and that was the only europium isotope
            then 1.0 will be returned.  If it was a natural uranium problem, and U235 was passed
            in, then the returned value would likely be something like 0.0072
   
   Note: if input is a Pu nuclide, and a Pu242-by-correlation is being applied, the returned value will take this into account; however,
   the uncertainty in Pu242 amount will not be accounted for (currently).
   
   Throws exception if \c nuclide was not in the specified relative efficiency curve, or other errors encountered.
   */
  std::pair<double,std::optional<double>> mass_enrichment_fraction( const SandiaDecay::Nuclide *nuclide,
                                                                   const size_t rel_eff_index ) const;

  /** Returns the mass ratio of two nuclides.
   
   Throws exception if either input \c nuclide is nullptr, or was not in the problem.
   
   TODO: add uncertainty, via returning pair<double,double>
   */
  double mass_ratio( const SandiaDecay::Nuclide *numerator, const SandiaDecay::Nuclide *denominator, const size_t rel_eff_index ) const;
  

  const NuclideRelAct &nucinfo( const SrcVariant src, const size_t rel_eff_index ) const;
  

  /** Returns the activity ratio of two nuclides, and if possible, the uncertainty of the ratio.

   The uncertainty will not be provided if covariance matrix is empty, or if it is Inf or NaN.

   Throws exception if either input \c nuclide is nullptr, or was not in the problem.
   */
  std::pair<double,std::optional<double>> activity_ratio( const SrcVariant &numerator, const size_t num_rel_eff_index,
                                              const SrcVariant &denominator, const size_t denom_rel_eff_index ) const;


  /** Walks the chain of activity ratio constraints to find the controlling source.

   @param [in,out] src The source to start from.
   @param [in] rel_eff_index The relative efficiency curve index to use.
   @param [out] multiple The multiple to multiply the activity ratio by.
   @returns True if a controlling source was found, false otherwise.
  
   throws exception if \c rel_eff_index or \c src is invalid.
  */
  bool walk_to_controlling_nuclide( SrcVariant &src, const size_t rel_eff_index, double &multiple ) const;

  /** Returns the uncertainty on the activity ratio of two nuclides.
  
  Throws exception if:
  - covariance was not able to be fit
  - either source passed in was not in the problem
  */
  double activity_ratio_uncertainty( SrcVariant numerator, size_t numerator_rel_eff_index, 
                                      SrcVariant denominator, size_t denominator_rel_eff_index ) const;
  
  /** Returns the relative activity of a nuclide.
  
  Please note that the interpretation of this value is a little tortured, as it depends on the relative efficiency curve.

  Throws exception if \c nuclide is nullptr, or was not in the problem.
  
  TODO: add uncertainty, via returning pair<double,double>
  */
  double rel_activity( const SrcVariant &src, const size_t rel_eff_index ) const;

  std::pair<double,double> rel_activity_with_uncert( const SrcVariant &src, const size_t rel_eff_index ) const;


  /** Returns the counts in all peaks for a source in the relative efficency curve.

  Note: it actually returns the sum of peak amplitudes, for all gammas within `m_final_roi_ranges`.

  Throws exception if \c src is nullptr, or was not in the problem, or any counts were inf or NaN.
  */
  double nuclide_counts( const SrcVariant &src, const size_t rel_eff_index ) const;

  /**  Gives the relative efficiency for a given energy. 
  */
  double relative_efficiency( const double energy, const size_t rel_eff_index ) const;
  
  /** Gives the relative efficiency for a given energy, as well as the uncertainty.
   
   Throws exception if covariance matrix is invalid.
   */
  std::pair<double,double> relative_efficiency_with_uncert( const double energy, const size_t rel_eff_index ) const;
  

  /** Get the index of specified nuclide within #m_rel_activities and #m_nonlin_covariance. */
  size_t nuclide_index( const SrcVariant &src, const size_t rel_eff_index ) const;

  /** Returns the index of the source activity within `m_final_parameters`.
   
   Age is +1.
  */
  size_t fit_parameters_index_for_source( const SrcVariant &src, const size_t rel_eff_index ) const;
  
  /** Returns result of `RelActCalc::rel_eff_eqn_js_function(...)` or
   `RelActCalc::physical_model_rel_eff_eqn_js_function(...)`
   */
  std::string rel_eff_eqn_js_function( const size_t rel_eff_index ) const;
  
  /** Returns JS function that represents the error bars of RelEffEqn, or "null" string if they arent avaiable. */
  std::string rel_eff_eqn_js_uncert_fcn( const size_t rel_eff_index ) const;
  
  /** Returns the updated energy calibration.
   
   Note: You could instead call `m_spectrum->energy_calibration()` to avoid computing from scratch.
   
   Will throw exception if runs into any issues.
   */
  std::shared_ptr<SpecUtils::EnergyCalibration> get_adjusted_energy_cal() const;
  
  static std::vector<std::vector<RelActCalcAuto::RelActAutoSolution::ObsEff>>
  fit_free_peak_amplitudes( const RelActCalcAuto::Options &options,
                            const RelActCalcAutoImp::RelActAutoCostFcn *cost_functor,
                            const std::vector<double> &parameters,
                           const RelActCalcAuto::RelActAutoSolution &solution );
  
  enum class Status : int
  {
    Success,
    NotInitiated,
    FailedToSetupProblem,
    FailToSolveProblem,
    UserCanceled
  };//
  
  RelActAutoSolution::Status m_status;
  std::string m_error_message;
  std::vector<std::string> m_warnings;
  
  Options m_options;
  
  std::shared_ptr<const RelActCalcAutoImp::RelActAutoCostFcn> m_cost_functor;
  
  /** The original foreground passed into `solve_ceres(...)`. */
  std::shared_ptr<const SpecUtils::Measurement> m_foreground;
  
  /** The original background passed into `solve_ceres(...)`. */
  std::shared_ptr<const SpecUtils::Measurement> m_background;
  
  /** The final fit parameters. */
  std::vector<double> m_final_parameters;
  
  /** The uncertainties on the final fit parameters; the sqrt of covariance matrix diagonal.
   Will be empty if covariance matrix failed to compute; otherwise same ordering as
   `m_final_parameters`.
   */
  std::vector<double> m_final_uncertainties;

  /** The scale factors used between what was fit for, and and the physical values (e.g., activity, age, etc.).
   
   If a paramater isnt used (e.g., age of nuclide when age isnt fit), then the scale factor may be 0.0; however,
   the scale factor may not be zero for these cases - (e.g., for energy calibration the scale factor that would 
   be used is returned, event if it is not fit for).
   */
  std::vector<double> m_parameter_scale_factors;
  
  /** The covariance matrix of the fit.
 
   Will be empty if computation of the covariance failed.
   the rows/columns have same size and ordering as `m_final_parameters`.

   Indexed as `m_covariance[row][col]`, by convention.
   */
  std::vector<std::vector<double>> m_covariance;

  /** The covariance matrix, converted to physical units - i.e., multiplied by parameter scale factors, or converted to physical units.

   Will be empty if computation of the covariance failed.

   Indexed as `m_phys_units_cov[row][col]`, by convention.
   */
  std::vector<std::vector<double>> m_phys_units_cov;

  /** The short names of the parameters. */
  std::vector<std::string> m_parameter_names;

  std::vector<bool> m_parameter_were_fit;
  
  /** This is a spectrum that will be background subtracted and energy calibrated, if those options
   where wanted.
   */
  std::shared_ptr<SpecUtils::Measurement> m_spectrum;
  
  std::vector<float> m_channel_counts;
  std::vector<float> m_channel_counts_uncerts;
  
  /** These are the peaks passed into the #solve function, or if no peaks were passed in, the
   auto-search peaks from the (potentially background subtracted) foreground.
   */
  std::vector<std::shared_ptr<const PeakDef>> m_spectrum_peaks;
  
  std::vector<RelActCalc::RelEffEqnForm> m_rel_eff_forms;
  
  /** The coefficients for the relative efficiency curve(s)
   * 
   * If Hoerl or external shieldings were shared between Physical Model curves, then the
   * coefficents fit for these values will be copied from the first physical model curve
   * to all the others.
   */
  std::vector<std::vector<double>> m_rel_eff_coefficients;

  /** The covariance matrix of the relative efficiency curve(s)
   * 
   * If Hoerl or external shieldings were shared between Physical Model curves, then the
   * covariance matrix entries these values will NOT have been copied from the first 
   * physical model curve to all the others - since I dont think this makes sense.
   */
  std::vector<std::vector<std::vector<double>>> m_rel_eff_covariance;
  
  /** The relative activities of each of the input nuclides.
   
   If a Pu242 correlation method was specified in #Options::pu242_correlation_method, then Pu242
   will NOT be in this variable, see #m_corrected_pu.

   Note: the gamma yields are nominal for the activity/age, and have not had the branching ratio
   uncertainty correction applied (if they were fit for).
   */
  std::vector<std::vector<NuclideRelAct>> m_rel_activities;
  
  std::vector<std::vector<double>> m_rel_act_covariance;
  
  FwhmForm m_fwhm_form;
  std::vector<double> m_fwhm_coefficients;
  std::vector<std::vector<double>> m_fwhm_covariance;
  
  std::vector<FloatingPeakResult> m_floating_peaks;
  
  /** The fit peaks.
   
   If energy calibration fitting was selected, these peaks will have been adjusted back to "true" energy.
   
   You would use these peaks for `m_spectrum`.
   
   Will include free-floating peaks (they will have nullptr sources).
   */
  std::vector<PeakDef> m_fit_peaks;
  
  /** Same as `m_fit_peaks`, but with the peaks seperated by relative efficiency curve (they will share peak continua between curves).
   Will not include free-floating peaks.
   */
  std::vector<std::vector<PeakDef>> m_fit_peaks_for_each_curve;
  
  /** The fit peaks, in the energy calibration of the spectrum (i.e., what you would display,
   if you are not updating the displayed spectrum for the adjusted energy cal).
   
   You would use these peaks for `m_foreground`.

   Note: if background subtraction is used, and you want to display the peaks on 
     non-background-subtracted data, please use `m_peaks_without_back_sub`.
   You would use these peaks for `m_foreground - m_background` (i.e., if you do background subtraction, then the continuums would not be correct,
   unless you background subtract from `m_foreground`).
   */
  std::vector<PeakDef> m_fit_peaks_in_spectrums_cal;
  
  /** Same as `m_fit_peaks_in_spectrums_cal`, but seperated by rel eff curve, and not including free-floating peaks. */
  std::vector<std::vector<PeakDef>> m_fit_peaks_in_spectrums_cal_for_each_curve;
  
  /** The fit peaks, in the spectrums cal, but with the continuums adjusted to display the peaks in 
   non-background-subtracted data.  If background subtraction wasnt used to do the fit, then these
   peaks will be the same as `m_fit_peaks_in_spectrums_cal`.
  */
  std::vector<PeakDef> m_peaks_without_back_sub;

  /** The measured rel. eff. for a energy; these are determined after the fit finishes, peaks are clustered (peaks within 1.5 sigma of
   each other are combined), then holding all other things (FWHM, energy cal, etc) constant, the amplitudes are allowed to freely float;
   then the relative efficieciency is determined from this.
   */
  struct ObsEff
  {
    /** The effective mean of all the peaks clutered. */
    double energy;
    /** The relative efficiency, as determined by the full solution to the problem. */
    double orig_solution_eff;
    /** The relative efficiency from the amplitude-unconstrained fit. */
    double observed_efficiency;
    double observed_efficiency_uncert;
    /** The scale factor to multiply the peak amplitude fit in the solution, to get what was freely fit for. */
    double observed_scale_factor;
    /** The unconstrained fit peak area, over its uncertainty.  Should be about 4, before we bother showing on chart */
    double num_sigma_significance;
    double cluster_lower_energy, cluster_upper_energy;
    double roi_lower_energy, roi_upper_energy;
    /** The unconstrained fit peak amplitude.  This amplitude includes contributions from all peaks, and all relative efficiency curves. */
    double fit_clustered_peak_amplitude;
    double fit_clustered_peak_amplitude_uncert;
    /** The starting amiplit*/
    double initial_clustered_peak_amplitude;
    /** The effective sigma, after clustering all the input peaks together that are within 1.5 sigma of each other. */
    double effective_sigma;
    double fraction_roi_counts;
    /* The mean +- 1-sigma is fully within the ROI */
    bool within_roi;

    /** The peaks, who have had their amplitudes scaled to the unconstrined fit value, who where clustered together.
     Ordered by largest peak first.
     */
    std::vector<PeakDef> fit_peaks;
  };//struct ObsEff

  /** The fit efficiencies for the peaks fit to data, allowing the amplitudes to freely float.
   Each `ObsEff` may be from one to a number of peaks that have been clustered together because they were effectively
   indistiguishable (means within 1.5 sigma is the criteria used).

   Peaks are clustered together (e.g., peaks less than 1.5 sigma from each other are combined), then the amplitudes and continuums
   are re-fit to the data, with all other parameters fixed at thier final fit values.  These peaks give a measure of how far off from the
   relative efficiency curve the data is from the actual solution.

   Will be empty if re-fitting process fails.

   These peaks are seperated by rel eff curve, do not include free-floating peaks, an are in "true" energy calibration of the spectrum.
   */
  std::vector<std::vector<ObsEff>> m_obs_eff_for_each_curve;


  /** When a ROI is #RoiRange::force_full_range is false, independent energy ranges will
   be assessed based on peak localities and expected counts; this variable holds the ROI
   ranges that were assessed and used to compute final answer.
   If all input RoiRanges had #RoiRange::force_full_range as true, and computation was
   successful then this variable will be equal to #m_options.rois.
   If computation is not successful, this variable may, or may not, be empty.

   Note: these ROIs are in true energy, not the energy of the spectrum.
   */
  std::vector<RoiRange> m_final_roi_ranges;

  /** These ROIs are in the energy of the spectrum. */
  std::vector<RoiRange> m_final_roi_ranges_in_spectrum_cal;
  
  
  /** This DRF will be the input DRF you passed in, if it was valid and had energy resolution info.
     
   If you passed in a valid DRF, but no energy resolution info, this DRF will have
   the resolution info added.
   
   If you didnt pass in a valid detector, then this will be a generic detector, with rough
   energy resolution derived from the spectrum, or at least the best we could do, before
   the fit.
   
   It may be useful to re-use the DRF on subsequent calls to  avoid the overhead of searching
   for all the peaks and such (although this isnt a garuntee this wont need to be done).
   */
  std::shared_ptr<const DetectorPeakResponse> m_drf;
  
  
  /** If #Options::pu242_correlation_method is not #PuCorrMethod::NotApplicable, then
   the result of the correction for Pu242 will be stored here.
   
   Note: this is a shared ptr just to avoid creating a copy constructor that would be needed
         if we make it a unique ptr...
   */
  std::vector<std::shared_ptr<const RelActCalc::Pu242ByCorrelationOutput<double>>> m_corrected_pu;
  
  /** Provides the input to Pu242 corrections, so you can compare against corrected values.
   
   Note that `mass_enrichment_fraction(...)` returns the corrected Pu fractions, so this is the only way to
   access the uncorrected mass fractions (unless you manually compute from relative activities)
   */
  std::vector<std::shared_ptr<const RelActCalc::Pu242ByCorrelationInput<double>>> m_uncorrected_pu;
  
  /** We will allow corrections to the first following number of energy calibration coefficients.
 
   Implemented for values of 2 or 3; lower channel energy calibration will only ever use up to 2.
   */
  constexpr static size_t sm_num_energy_cal_pars = 2;
  
  /** It seems parameters centered around zero maybe don't work quite as well, we'll center the energy calibration
   parameters at 1, and allow to vary between 0 and 2.
   */
  constexpr static double sm_energy_par_offset = 1.0;
  
  /** We will allow the energy offset to vary by the minimum of +-15 keV or `sm_energy_offset_range_fwhm` times lower energy FWHM . This is arbitrarily chosen. */
  constexpr static double sm_energy_offset_range_keV = 15.0;  // Allow +-15 keV offset adjust
  /** If we are off by more than 1 FWHM on the offset, we are likely to be falling into some false minima - this has not been investigated. */
  constexpr static double sm_energy_offset_range_fwhm = 1.0;
  
  /** We will allow the energy gain to vary by +-20 keV, or 4 FWHM - whichever is less, at the right side of the spectrum.
   This is arbitrarily chosen.
   */
  constexpr static double sm_energy_gain_range_keV   = 20.0; // Allow up to 20 keV energy range gain adjust
  constexpr static double sm_energy_gain_range_fwhm  = 4.0; //Allow up to 4 FWHM range gain adjust
  
  /** If `sm_num_energy_cal_pars == 3` we will allow up to 5 keV adjustment to quadratic energy range, at the right-hand
   side of the spectrum.
   */
  constexpr static double sm_energy_quad_range_keV   = 10.0;  // Allow 10 keV energy range quad correction
  
  /** With numeric differentiation, we run into an issue where the initial small delta use (e.x. 1.0E-6 of parameter value)
   isnt quite enough to actually have a change (and might even get lost in float precision of energy cal pars) - so we will "fix" this by
   multiplying the change of the value away from `sm_energy_par_offset`, by a multiple like 100, so this way the parameter will
   actually fit.
   We then adjust the limits of the parameter (`sm_energy_offset_range_keV` and `sm_energy_gain_range_keV`) to take
   this into account, so the gain parameter will go between 0.8 and 1.2, for +- 20 keV effect on the spectrum
   */
  constexpr static double sm_energy_cal_multiple    = 100.0;
  
  /** The fit energy calibration adjustments.
   
   To get the effect, in keV to the spectrum (e.g., on right-hand side of spectrum), you would do:
     `(m_energy_cal_adjustments[i]/sm_energy_par_offset - 1)*sm_energy_cal_multiple`
   */
  std::array<double,sm_num_energy_cal_pars> m_energy_cal_adjustments;
  
  /** Whether each energy cal parameter was fit.
   Which parameters are fit are subject to number and energy range of ROIs.
   */
  std::array<bool,sm_num_energy_cal_pars> m_fit_energy_cal;

  size_t m_num_deviations_fit = 0;
  
  /** The fitted deviation pair offsets for non-linear energy calibration.

   Each entry is (anchor_energy_keV, offset_keV).  The offset at lower and upper ROIs is fixed to 0.
   Entries are sorted by anchor energy.

   Only valid if `Options::energy_cal_type` is `EnergyCalFitType::NonLinearFit`.
   */
  std::vector<std::pair<double, double>> m_deviation_pair_offsets;

  /** The index of the first parameter that will be used to adjust the peak amplitude. 
   
   Will be valid only if `m_options.additional_br_uncert > 0.0`.
  */
  size_t m_add_br_uncert_start_index = std::numeric_limits<size_t>::max();
  
  /** The energy ranges cooresponding to the additional peak amplitude uncertainty paramaters.
   
   Will only be valid/filled-out if `m_options.additional_br_uncert > 0.0`.
   
   \sa m_add_br_uncert_start_index
   */
  std::vector<std::pair<double,double>> m_peak_ranges_with_uncert;
  
  /** */
  double m_chi2;
  
  /** The number of degrees of freedom in the fit for equation parameters.
   
   Note: this is the number data channels used, minus the number of (estimated effective) fit paramaters.
   */
  size_t m_dof;
  
  /** A struct to hold information about the Physical Model result of the fit. */
  struct PhysicalModelFitInfo
  {
    /** Info on an individual shielding. */
    struct ShieldInfo
    {
      std::shared_ptr<const Material> material;
      std::optional<double> atomic_number;        /// Will not be present if material
      std::optional<double> atomic_number_uncert; /// Only present if not a material, and AN was fit; sqrt of diagonal of covariance matrix
      bool atomic_number_was_fit;                 /// If AN was fit
      double areal_density;                       /// In PhysicalUnits units
      std::optional<double> areal_density_uncert; /// Only present if AD was fit; sqrt of diagonal of covariance matrix
      bool areal_density_was_fit;                 /// If AD was fit
    };//struct ShieldInfo
    
    /** Self attenuator info will only be present if was input into fit. */
    std::optional<ShieldInfo> self_atten;
    std::vector<ShieldInfo> ext_shields;

    /** If `RelEffCurveInput::shielded_by_other_phys_model_curve_shieldings` was specified, the "outer" curves shieldings
     (both self-atten and external atten) will copied to this variable.
     */
    std::vector<ShieldInfo> shields_from_other_curves;

    // Modified Hoerl corrections only present if fitting the Hoerl function was selected
    //  Uncertainties are just sqrt of covariance diagnal
    std::optional<double> hoerl_b, hoerl_b_uncert;
    std::optional<double> hoerl_c, hoerl_c_uncert;
  };//struct PhysicalModelFitInfo
  
  /** The curve information for Physical Model curves.
   */
  std::vector<std::optional<PhysicalModelFitInfo>> m_phys_model_results;
  
  
  /** The number of evaluation calls it took L-M to reach a solution.
   Only useful for debugging and curiosity.
   */
  int m_num_function_eval_solution;
  
  /** The number of evaluation calls it took to reach a solution, and compute final covariance. */
  int m_num_function_eval_total;
  
  /** This is the total time spent solving the problem. */
  int m_num_microseconds_eval;
  
  /** This is the total time spent in the `eval(...)` function. */
  int m_num_microseconds_in_eval;
};//struct RelEffSolution


/**
 
 @param rel_eff_order The number of energy dependent terms to have in the relative efficiency
        equation (e.g., one more parameter than this will be fit for).
 */
RelActAutoSolution solve( const Options options,
                         std::shared_ptr<const SpecUtils::Measurement> foreground,
                         std::shared_ptr<const SpecUtils::Measurement> background,
                         std::shared_ptr<const DetectorPeakResponse> drf,
                         std::vector<std::shared_ptr<const PeakDef>> all_peaks,
                         std::shared_ptr<std::atomic_bool> cancel_calc = nullptr
                         );



}//namespace RelActCalcAuto

#endif //RelActCalcAuto_h
