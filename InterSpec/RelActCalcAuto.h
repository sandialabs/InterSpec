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

#include <atomic>
#include <string>
#include <memory>
#include <vector>
#include <ostream>
#include <optional>

#include "InterSpec/PeakDef.h" //for PeakContinuum::OffsetType and PeakDef::SkewType


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
}//namespace SpecUtils

namespace SandiaDecay
{
  struct Nuclide;
}

namespace RelActCalc
{
  enum class PuCorrMethod : int;
  enum class RelEffEqnForm : int;
  struct Pu242ByCorrelationOutput;
  struct PhysicalModelShieldInput;
}//namespace RelActCalc



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
 - Should add option to give each branching-ratio an independent uncertainty (like 1% or something),
   as well the same for FWHM
 */

namespace RelActCalcAuto
{

int run_test();

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
  
  /** Dont allow cutting range down based on peaks not being anywhere reasonably near edge. */
  bool force_full_range = false;
  
  /** If we have a peak right-on the edge, allow the ROI to extend out to a few FWHM.
   If #force_full_range is true, this value must be false.
   */
  bool allow_expand_for_peak_width = false;
  
  static const int sm_xmlSerializationVersion = 0;
  void toXml( ::rapidxml::xml_node<char> *parent ) const;
  void fromXml( const ::rapidxml::xml_node<char> *parent );
};//struct RoiRange


/** Struct to specify a nuclide to use for doing relative-efficiency/activity calc.
 */
struct NucInputInfo
{
  const SandiaDecay::Nuclide *nuclide = nullptr;
  
  /** Age in units of PhysicalUnits (i.e., 1.0 == second).
   
   If the age is being fit, this will be used as the initial age to start with.
   
   Must not be negative.
   */
  double age = -1.0;
  
  
  bool fit_age = false;
  
  /** Energy corresponding to SandiaDecay::EnergyRatePair::energy */
  std::vector<double> gammas_to_exclude;
  
  std::string peak_color_css;
  
  static const int sm_xmlSerializationVersion = 0;
  void toXml( ::rapidxml::xml_node<char> *parent ) const;
  void fromXml( const ::rapidxml::xml_node<char> *parent );
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
  Polynomial_6
};//enum ResolutionFnctForm

/** Returns string representation of the #FwhmForm.
 String returned is a static string, so do not delete it.
 */
const char *to_str( const FwhmForm form );

/** Converts from the string representation of #FwhmForm to enumerated value.
 
 Throws exception if invalid string (i.e., any string not returned by #to_str(FwhmForm) ).
 */
FwhmForm fwhm_form_from_str( const char *str );

/** Evaluates the FWHM equation for the input energy, returning the FWHM. */
float eval_fwhm( const float energy, const FwhmForm form, const std::vector<float> &coeffs );

size_t num_parameters( const FwhmForm eqn_form );


struct Options
{
  Options();
 
  // TODO: 
  // - Fix RelEff to flat 1.0 - dont fit
  //  - Do not fit FWHM - use DRF, or if that doesnt have it, use FWHM eqn fit from all peaks in spectrum
  //    This should allow getting a rough idea of which gamma lines _could_ be statistically significant, and hence be used.  And also, for generic nuclides, it will allow fitting small numbers of peaks; so like we could still use the "auto" rel act stuff to fit peaks, even for Cs137
  
  /** Whether to allow making small adjustments to the gain and/or offset of the energy calibration.
   
   Which coefficients are fit will be determined based on energy ranges used.
   */
  bool fit_energy_cal;
  
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
  
  /** The functional form of the FWHM equation to use.  The coefficients of this equation will
   be fit for across all energy ranges.
   */
  FwhmForm fwhm_form;
  
  /** Optional title of the spectrum; used as title in HTML and text summaries. */
  std::string spectrum_title;
  
  /** The method to use for Pu-242 mass-enrichment estimation (all methods use correlation to other
   Pu isotopics).
   
   When specified #RelActAutoSolution::m_corrected_pu will be filled-out (unless an error occurred
   while doing the correction - which shouldnt really ever happen); this will then be used by the
   #RelActAutoSolution::mass_enrichment_fraction function, and activity ratio function.
   
   Defaults to #PuCorrMethod::NotApplicable; if any other value is specified, and sufficient
   plutonium isotopes for that method are not in the problem, then finding the solution will
   fail.
   */
  RelActCalc::PuCorrMethod pu242_correlation_method;
  
  
  /** Peak skew to apply to the entire spectrum.
   
   Under development: currently, if total energy range being fit is less than 100 keV, then all peaks will share the same skew.
   Otherwise a linear energy dependance will be assumed, where the fitting parameters will be for the spectrums lower
   energy, and the spectrums upper energy; not all skew parameters are allowed to vary with energy; e.g., the Crystal Ball
   power law is not allowed to have an energy dependence (see `PeakDef::is_energy_dependent(...)`).
   */
  PeakDef::SkewType skew_type;
  
  /** Only used for `RelEffEqnForm::FramPhysicalModel` - shielding definitions for self-attenuating sources.
   
   May be nullptr if physical model.
   */
  std::shared_ptr<const RelActCalc::PhysicalModelShieldInput> phys_model_self_atten;
  
  /** Only used for `RelEffEqnForm::FramPhysicalModel` - shielding definitions for external attenuators. */
  std::vector<std::shared_ptr<const RelActCalc::PhysicalModelShieldInput>> phys_model_external_atten;

  /** If true, fit the modified Hoerl equation form for the physical model.
   * If false, do not fit the modified Hoerl equation form (its value will be constant value of 1.0).
   *
   * Ignored if not using `RelActCalc::RelEffEqnForm::FramPhysicalModel`.
  */
  bool phys_model_use_hoerl = true;
  
  /** Version history:
   - 20250117: incremented to 1 to handle FramPhysicalModel; if not this model, will still write version 0.
   */
  static const int sm_xmlSerializationVersion = 1;
  rapidxml::xml_node<char> *toXml( ::rapidxml::xml_node<char> *parent ) const;
  
  /** Sets the member variables from an XML element created by `toXml(...)`.
   @param parent An XML element with name "Options".
   @param materialDB The material database to use to retrieve a material from for the #PhysicalModelShieldInput;
          for other equation types, or for AN/AD defined shields, this isnt used/required.
   */
  void fromXml( const ::rapidxml::xml_node<char> *parent, MaterialDB *materialDB );
};//struct Options


struct NuclideRelAct
{
  const SandiaDecay::Nuclide *nuclide;
  
  double age;
  double age_uncertainty;
  bool age_was_fit;
  
  double rel_activity;
  double rel_activity_uncertainty;
  
  /** The energy and its respective number of gammas per decay for this nuclide, at its age. */
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
  
  RelActCalcAuto::Options options;
  std::vector<RelActCalcAuto::RoiRange> rois;
  std::vector<RelActCalcAuto::NucInputInfo> nuclides;
  std::vector<RelActCalcAuto::FloatingPeak> floating_peaks;
    
  bool background_subtract;
  RelActCalc::PuCorrMethod pu_correlation_method;
  
  bool show_ref_lines;
  double lower_display_energy;
  double upper_display_energy;
  
  /** Returns XML node added; i.e., will have name "RelActCalcAuto" */
  ::rapidxml::xml_node<char> *serialize( ::rapidxml::xml_node<char> *parent ) const;
  void deSerialize( const rapidxml::xml_node<char> *base_node, MaterialDB *materialDb );
};//struct RelActAutoGuiState
  
  
struct RelActAutoSolution
{
  RelActAutoSolution();
  
  std::ostream &print_summary( std::ostream &strm ) const;
  
  void print_html_report( std::ostream &strm ) const;
  
  /** Prints out the JSON data the JS the RelEff chart accepts.
   
   Note: currently this code largely duplicates #RelEffChart::setData, so need to refactor.
   */
  void rel_eff_json_data( std::ostream &json, std::ostream &css ) const;
  
  /** Prints the txt version of relative eff eqn. */
  std::string rel_eff_txt( const bool html_format ) const;
  
  /** Returns the fractional amount of an element (by mass), the nuclide composes.
   
   Throws exception if \c nuclide was not in the problem.
   
   @param nuclide The nuclide of interest.
   @returns The fraction of mass (so between 0.0 and 1.0), the input nuclide is responsible for
            in the solution.  Ex., if Eu152 is passed in, and that was the only europium isotope
            then 1.0 will be returned.  If it was a natural uranium problem, and U235 was passed
            in, then the returned value would likely be something like 0.0072
   
   TODO: add uncertainty, via returning pair<double,double>
   */
  double mass_enrichment_fraction( const SandiaDecay::Nuclide *nuclide ) const;
  
  /** Returns the mass ratio of two nuclides.
   
   Throws exception if either input \c nuclide is nullptr, or was not in the problem.
   
   TODO: add uncertainty, via returning pair<double,double>
   */
  double mass_ratio( const SandiaDecay::Nuclide *numerator, const SandiaDecay::Nuclide *denominator ) const;
  
  /** Returns the activity ratio of two nuclides.
   
   Throws exception if either input \c nuclide is nullptr, or was not in the problem.
   
   TODO: add uncertainty, via returning pair<double,double>
   */
  double activity_ratio( const SandiaDecay::Nuclide *numerator, const SandiaDecay::Nuclide *denominator ) const;
  
  /** Get the index of specified nuclide within #m_rel_activities and #m_nonlin_covariance. */
  size_t nuclide_index( const SandiaDecay::Nuclide *nuclide ) const;
  
  /** Returns result of `RelActCalc::rel_eff_eqn_js_function(...)` or
   `RelActCalc::physical_model_rel_eff_eqn_js_function(...)`
   */
  std::string rel_eff_eqn_js_function() const;
  
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
  
  std::shared_ptr<const SpecUtils::Measurement> m_foreground;
  std::shared_ptr<const SpecUtils::Measurement> m_background;
  
  /** The final fit parameters. */
  std::vector<double> m_final_parameters;
  
  /** The uncertainties on the final fit parameters; the sqrt of covariance matrix diagonal.
   Will be empty if covariance matrix failed to compute; otherwise same ordering as
   `m_final_parameters`.
   */
  std::vector<double> m_final_uncertainties;
  
  /** The covariance matrix of the fit.
 
   Will be empty if computation of the covariance failed.
   the rows/columns have same size and ordering as `m_final_parameters`.
   */
  std::vector<std::vector<double>> m_covariance;
  
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
  
  RelActCalc::RelEffEqnForm m_rel_eff_form;
  
  std::vector<double> m_rel_eff_coefficients;
  std::vector<std::vector<double>> m_rel_eff_covariance;
  
  /** The relative activities of each of the input nuclides.
   
   If a Pu242 correlation method was specified in #Options::pu242_correlation_method, then Pu242
   will NOT be in this variable, see #m_corrected_pu.
   */
  std::vector<NuclideRelAct> m_rel_activities;
  
  std::vector<std::vector<double>> m_rel_act_covariance;
  
  FwhmForm m_fwhm_form;
  std::vector<double> m_fwhm_coefficients;
  std::vector<std::vector<double>> m_fwhm_covariance;
  
  std::vector<FloatingPeakResult> m_floating_peaks;
  
  std::vector<PeakDef> m_fit_peaks;
  
  std::vector<RoiRange> m_input_roi_ranges;
  
  /** When a ROI is #RoiRange::force_full_range is false, independent energy ranges will
   be assessed based on peak localities and expected counts; this variable holds the ROI
   ranges that were assessed and used to compute final answer.
   If all input RoiRanges had #RoiRange::force_full_range as true, and computation was
   successful then this variable will be equal to #m_input_roi_ranges.
   If computation is not successful, this variable may, or may not, be empty.
   */
  std::vector<RoiRange> m_final_roi_ranges;
  
  
  /** This DRF will be the input DRF you passed in, if it was valid and had energy resolution info.
   Otherwise, this will be a "FLAT" DRF with a rough energy resolution function fit from the
   foreground spectrum; in this case, it may be useful to re-use the DRF on subsequent calls to
   avoid the overhead of searching for all the peaks and such.
   */
  std::shared_ptr<const DetectorPeakResponse> m_drf;
  
  
  /** If #Options::pu242_correlation_method is not #PuCorrMethod::NotApplicable, then
   the result of the correction for Pu242 will be stored here.
   
   Note: this is a shared ptr just to avoid creating a copy constructor that would be needed
         if we make it a unique ptr...
   */
  std::shared_ptr<const RelActCalc::Pu242ByCorrelationOutput> m_corrected_pu;
  
  double m_energy_cal_adjustments[2];
  bool m_fit_energy_cal[2];
  
  
  /** */
  double m_chi2;
  
  /** The number of degrees of freedom in the fit for equation parameters.
   
   Note: this is currently just a
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
    
    // Modified Hoerl corrections only present if fitting the Hoerl function was selected
    //  Uncertainties are just sqrt of covariance diagnal
    std::optional<double> hoerl_b, hoerl_b_uncert;
    std::optional<double> hoerl_c, hoerl_c_uncert;
  };//struct PhysicalModelFitInfo
  
  
  std::optional<PhysicalModelFitInfo> m_phys_model_result;
  
  
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
                         const std::vector<RoiRange> energy_ranges,
                         const std::vector<NucInputInfo> nuclides,
                         const std::vector<FloatingPeak> extra_peaks,
                         std::shared_ptr<const SpecUtils::Measurement> foreground,
                         std::shared_ptr<const SpecUtils::Measurement> background,
                         std::shared_ptr<const DetectorPeakResponse> drf,
                         std::vector<std::shared_ptr<const PeakDef>> all_peaks,
                         std::shared_ptr<std::atomic_bool> cancel_calc = nullptr
                         );


}//namespace RelActCalcAuto

#endif //RelActCalcAuto_h
