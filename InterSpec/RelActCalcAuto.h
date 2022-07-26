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

#include "InterSpec/PeakDef.h" //for PeakContinuum::OffsetType


// Forward declarations
class DetectorPeakResponse;

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
  enum class RelEffEqnForm : int;
}

/*
 We need to take in:
 
 foreground_spectrum;
 background_spectrum; //optional
 
 struct RoiRanges{ double lower_energy, upper_energy; continuum_type; bool force_full_range; }
 vector<RoiRanges>: roi_ranges;
 
 struct NucInputInfo{ const SandiaDecay::Nuclide *nuc; double age; bool fit_age; vector<double> gammas_to_release; }
 vector<NucInputInfo> nuclides;
 
 struct FreeFloatPeak{ double energy; bool release_fwhm; }
 vector<FreeFloatPeak> free_floating_peaks;
 
 options{ bool fit_energy_cal; RelEffEqnType; size_t RelEffEqnOrder; FwhmFcnForm/Order; }
 
 
 For deciding peak ROIs:
 - First use whole specified ROI, and fit for a solution
 - Eliminate statistically insignificant gamma lines; then use remaining gammas to decide peak ROI widths (and if should be further broken up?)
 
 To get rid of the degeneracy of rel-eff and rel-act:
 - Add an extra term to the Chi2 that forces the rel-act equation to be 1.0 at the lowest energy; maybe multiple difference of current rel-eff equation from 1.0 by the area of the largest peak in any of the ROIs
 
 Parameters to fit:
 - Energy offset; Energy gain adjust //will be fixed to begin with
 - Fwhm Parameters  // Will be `num_parameters(FwhmForm)` number of these
 - RelEff coefs     // There will be one more than the "order" of relative eff eqn
 - {Activity,Age}   // One pair for each nuclide.  Age will be mostly fixed, but can also be tied to another nuclide of the same elements age
 - Free float peaks {amplitude, width-multiple} // The width multiple is of what the Fwhm Parameters predict (e.g., 1.0 says use parameter value).
 
 
 Further things to consider:
 - Sum peaks - not totally sure how to handle - can there just be one parameter that controls strength of peak random summing, and then another for cascade?
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
  
  bool release_fwhm = false;
  
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
  
  /** The equation FWHM = sqrt( A_0 + A_1*Energy + A_2/Energy ).
   This is what FRAM uses.
   */
  SqrtEnergyPlusInverse,
  
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
   */
  size_t rel_eff_eqn_order;
  
  /** The functional form of the FWHM equation to use.  The coefficients of this equation will
   be fit for across all energy ranges.
   */
  FwhmForm fwhm_form;
  
  /** Optional title of the spectrum; used as title in HTML and text summaries. */
  std::string spectrum_title;
  
  static const int sm_xmlSerializationVersion = 0;
  rapidxml::xml_node<char> *toXml( ::rapidxml::xml_node<char> *parent ) const;
  void fromXml( const ::rapidxml::xml_node<char> *parent );
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


struct RelActAutoSolution
{
  RelActAutoSolution();
  
  std::ostream &print_summary( std::ostream &strm ) const;
  
  void print_html_report( std::ostream &strm ) const;
  
  /** Prints out the JSON data the JS the RelEff chart accepts.
   
   Note: currently this code largely duplicates #RelEffChart::setData, so need to refactor.
   */
  void rel_eff_json_data( std::ostream &json, std::ostream &css ) const;
  
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
  
  /** Get the index of specified nuclide within #m_rel_activities and #m_rel_act_covariance. */
  size_t nuclide_index( const SandiaDecay::Nuclide *nuclide ) const;
  
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
  
  std::vector<double> m_final_parameters;
  
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
  std::vector<RelActCalcAuto::RoiRange> m_final_roi_ranges;
  
  
  /** This DRF will be the input DRF you passed in, if it was valid and had energy resolution info.
   Otherwise, this will be a "FLAT" DRF with a rough energy resolution function fit from the
   foreground spectrum; in this case, it may be useful to re-use the DRF on subsequent calls to
   avoid the overhead of searching for all the peaks and such.
   */
  std::shared_ptr<const DetectorPeakResponse> m_drf;
  
  
  double m_energy_cal_adjustments[2];
  bool m_fit_energy_cal[2];
  
  
  /** */
  double m_chi2;
  
  /** The number of degrees of freedom in the fit for equation parameters.
   
   Note: this is currently just a
   */
  size_t m_dof;
  
  /** The number of evaluation calls it took L-M to reach a solution.
   Only useful for debugging and curiosity.
   */
  int m_num_function_eval_solution;
  
  /** The number of evaluation calls it took to reach a solution, and compute final covariance. */
  int m_num_function_eval_total;
  
  int m_num_microseconds_eval;
};//struct RelEffSolution

/**
 
 @param rel_eff_order The number of energy dependent terms to have in the relative efficiency
        equation (e.g., one more parameter than this will be fit for).
 */
RelActAutoSolution solve( Options options,
                         std::vector<RoiRange> energy_ranges,
                         std::vector<NucInputInfo> nuclides,
                         std::vector<FloatingPeak> extra_peaks,
                         std::shared_ptr<const SpecUtils::Measurement> foreground,
                         std::shared_ptr<const SpecUtils::Measurement> background,
                         std::shared_ptr<const DetectorPeakResponse> drf,
                         std::vector<std::shared_ptr<const PeakDef>> all_peaks,
                         std::shared_ptr<std::atomic_bool> cancel_calc = nullptr
                         );


}//namespace RelActCalcAuto

#endif //RelActCalcAuto_h
