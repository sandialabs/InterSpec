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

#include <string>
#include <memory>
#include <vector>

#include "InterSpec/PeakDef.h" //for PeakContinuum::OffsetType


// Forward declarations
class DetectorPeakResponse;

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
 - Fwhm Parameters
 - RelEff coefs
 - Activities
 - Ages //will mostly be fixed
 - Free float peaks {amplitude, width} // if width is negative, then use Fwhm eqn
 
 */

namespace RelActCalcAuto
{

/** Struct to specify an energy range to consider for doing relative-efficiency/activity calc.
 */
struct RoiRange
{
  double lower_energy = -1.0;
  double upper_energy = -1.0;
  PeakContinuum::OffsetType continuum_type = PeakContinuum::OffsetType::Quadratic;
  
  /** Dont allow cutting range down based on peaks not being anywhere reasonably near edge. */
  bool force_full_range = false;
  
  /** If we have a peak right-on the edge, allow the ROI to extend out to a few FWHM.
   If #force_full_range is true, this value must be false.
   */
  bool allow_expand_for_peak_width = false;
};//struct RoiRange


/** Struct to specify a nuclide to use for doing relative-efficiency/activity calc.
 */
struct NucInputInfo
{
  const SandiaDecay::Nuclide *nuclide = nullptr;
  
  /** Must not be negative. */
  double age = -1.0;
  
  
  bool fit_age = false;
  
  /** Energy corresponding to SandiaDecay::EnergyRatePair::energy */
  std::vector<double> gammas_to_exclude;
};//struct NucInputInfo


/** A peak at a specific energy, that has a free-floating amplitude, and may also have a floating
 FWHM.
 */
struct FloatingPeak
{
  double energy = -1.0;
  bool release_fwhm = false;
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
  
  /** The linear polynomial: e.g. FWHM = sqrt(A_0 + A_1*1/1000) */
  Polynomial_2,
  
  /** The quadratic polynomial: e.g. FWHM = sqrt(A_0 + A_1*0.001*Energy + A_2*0.000001*Energy*Energy */
  Polynomial_3,
  Polynomial_4,
  Polynomial_5,
  Polynomial_6
};//enum ResolutionFnctForm



struct Options
{
  Options();
  
  bool fit_energy_cal;
  
  RelActCalc::RelEffEqnForm rel_eff_eqn_type;
  
  /** The number of energy dependent terms to use for the relative efficiency equation.  I.e., the
   number of terms fit for will be one more than this value.
   */
  size_t rel_eff_eqn_order;
  
  FwhmForm fwhm_form;
};//struct Options


struct NuclideRelAct
{
  const SandiaDecay::Nuclide *nuclide;
  
  double age;
  double age_uncertainty;
  bool age_was_fit;
  
  double activity;
  double activity_uncertainty;
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
  
  enum class Status : int
  {
    Success,
    NotInitiated,
    FailedToSetupProblem,
    FailToSolveProblem
  };//
  
  RelActAutoSolution::Status m_status;
  std::string m_error_message;
  
  std::shared_ptr<const SpecUtils::Measurement> m_foreground;
  std::shared_ptr<const SpecUtils::Measurement> m_background;
  
  /** This is a spectrum that will be background subtracted and energy calibrated, if those options
   where wanted.
   */
  std::shared_ptr<SpecUtils::Measurement> m_spectrum;
  
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
  
  /** This DRF will be the input DRF you passed in, if it was valid and had energy resolution info.
   Otherwise, this will be a "FLAT" DRF with a rough energy resolution function fit from the
   foreground spectrum; in this case, it may be useful to re-use the DRF on subsequent calls to
   avoid the overhead of searching for all the peaks and such.
   */
  std::shared_ptr<const DetectorPeakResponse> m_drf;
  
  double m_energy_cal_adjustments[2];
  bool m_fit_energy_cal_adjustments;
  
  
  /** */
  double m_chi2;
  
  /** The number of degrees of freedom in the fit for equation parameters.
   
   Note: this is currently just a
   */
  size_t m_dof;
  
  /** The number steps it took L-M to reach a solution. Only useful for debugging and curiosity. */
  int m_num_function_eval;
  
  
  int m_num_microseconds_eval;
};//struct RelEffSolution


RelActAutoSolution solve( Options options,
                         std::vector<RoiRange> energy_ranges,
                         std::vector<NucInputInfo> nuclides,
                         std::vector<FloatingPeak> extra_peaks,
                         FwhmForm fwhm_form,
                         RelActCalc::RelEffEqnForm rel_eff_form,
                         size_t rel_eff_order,
                         std::shared_ptr<const SpecUtils::Measurement> foreground,
                         std::shared_ptr<const SpecUtils::Measurement> background,
                         const std::shared_ptr<const DetectorPeakResponse> drf
                         );


}//namespace RelActCalcAuto

#endif //RelActCalcAuto_h
