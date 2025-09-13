#ifndef RelActCalc_h
#define RelActCalc_h
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

#include <tuple>
#include <vector>
#include <cstddef>
#include <optional>
#include <functional>

// Forward declarations
struct Material;
class MaterialDB;
class DetectorPeakResponse;
class PeakDef;

namespace SpecUtils
{
  class Measurement;
}

namespace rapidxml
{
  template<class Ch> class xml_node;
}

/** The \c RelActCalc namespace captures things that are common between the "Auto" and "Manual"
 relative eff/act calculations.
 */
namespace RelActCalc
{

/** The available forms of relative efficiency equations.
  
  Where x in energy (in keV), and y is area/br.
  */
enum class RelEffEqnForm : int
{
  /** y = a + b*ln(x) + c*(ln(x))^2 + d*(ln(x))^3 + ...
   */
  LnX,
  
  /** y = exp( a + b*x + c/x + d/x^2 + e/x^3 + ... )
   */
  LnY,
  
  
  /** y = exp( a  + b*(lnx) + c*(lnx)^2 + d*(lnx)^3 + ... )
   */
  LnXLnY,
  
  /** y = exp( a + b/x^2 + c*(lnx) + d*(lnx)^2 + e*(lnx)^3 )
   
   See LAUR-11-03005.
   Note: FRAM uses MeV, instead of keV, but it was verified using either
   units gives the same relative activities, so we this library will
   stick to using keV for consistency.
   */
  FramEmpirical,
  
  /** The (slightly modified) FRAM Physical Model.
   
   See pages 98 and 99 of "FRAM Version 6.1, User Manual "
   
      
   Vo, Duc T., and Sampson, Thomas E. FRAM, Version 6.1 User Manual. United States: N. p., 2020. Web. doi:10.2172/1599022.
   https://www.osti.gov/biblio/1599022
   */
  FramPhysicalModel
};//enum class RelEffEqnForm

/** Returns string representation of RelEffEqnForm, i.e., "LnX", "LnY", "LnXLnY", or "FRAM Empirical" */
const char *to_str( const RelEffEqnForm form );

RelEffEqnForm rel_eff_eqn_form_from_str( const char *str );

/** Returns a human-readable representation of the relative efficiency equation.
 @param eqn_form The form of the equation used.
 @param coefs The coefficients for the equation.
 */
std::string rel_eff_eqn_text( const RelEffEqnForm eqn_form,
                             const std::vector<double> &coefs );

/** Returns a javascript expression of the relative efficiency equation.
 
 Example return value: "function(x){ return Math.exp( 1.2 + x + ... ); }"
 */
std::string rel_eff_eqn_js_function( const RelEffEqnForm eqn_form, const std::vector<double> &coefs );

/** Evaluate the relative efficiency equation at a given energy, for a specific form of the equation, and coefficients.
 
 Can not be called for `RelEffEqnForm::FramPhysicalModel`, will throw exception - see `eval_physical_model_eqn()`.
 */
double eval_eqn( const double energy, const RelEffEqnForm eqn_form,
                const double * const coefs, const size_t num_coefs );

/** A convenience call signature for the above #eval_eqn */
double eval_eqn( const double energy, const RelEffEqnForm eqn_form, const std::vector<double> &coefs );

/** Evaluate the uncertainty in the relative efficiency equation.
 
 @param coefs The coefficients of the nominal equation.
 @param covariance The covariance matrix of the ``log'' transformed relative
        efficiency equation coefficients; i.e., as computed using the least 
        linear squares method - not as computed by ceres.

 Can not be called for `RelEffEqnForm::FramPhysicalModel`, will throw exception 
 - see `RelEffSolution::rel_eff_eqn_uncert(...)`.
 
 TODO: this function can return NaN values - when the square uncertainty is negative...
 */
double eval_eqn_uncertainty( const double energy, const RelEffEqnForm eqn_form,
                             const std::vector<double> &coefs,
                             const std::vector<std::vector<double>> &covariance );


/** Function to back-calculate Relative Activities and masses of a set of nuclides, a specified time interval previous.

 @param back_decay_time A positive time interval to back-decay nuclides to.
 @param nuclide_rel_acts The nuclides, and their relative activities at the time of measurement.

 Returns vector of `{nuclide, RelActivityAtTimeZero, RelMassAtTimeZero}`, where the RelativeMasses are as a fraction of all nuclides passed in.
 */
// Implemented in RelActCalc_imp.hpp
template<typename T>
inline std::vector<std::tuple<const SandiaDecay::Nuclide *,T,T>> back_decay_relative_activities(
                                    const T back_decay_time,
                                    std::vector<std::tuple<const SandiaDecay::Nuclide *,T>> &nuclide_rel_acts );


/** A struct to specify the relative masses of the Pu and Am241 nuclides that
 are observable by gamma spec.
 
 Masses do not need to be normalized to 1, but their ratios to each other do
 need to be correct.
 */
template<typename T>
struct Pu242ByCorrelationInput
{
  T pu238_rel_mass = T(0.0);
  T pu239_rel_mass = T(0.0);
  T pu240_rel_mass = T(0.0);
  T pu241_rel_mass = T(0.0);
  /** A value to capture all other non-Pu242 plutonium isotopes (I dont expect to ever be used) */
  T other_pu_mass  = T(0.0);

  /** The age of the Pu, so this way the mass fractions can be back-corrected to T=0, to ever so slightly increase the accuracy of the correction. */
  T pu_age = T(0.0);
};//struct Pu242ByCorrelationInput


/** The nuclide mass fractions (so between 0 and 1 - with all nuclides summing to 1.0)
 as determined after accounting for the predicted amount (by correlation) of Pu242.
 */
template<typename T>
struct Pu242ByCorrelationOutput
{
  T pu238_mass_frac;
  T pu239_mass_frac;
  T pu240_mass_frac;
  T pu241_mass_frac;
  T pu242_mass_frac;
  
  /** If enrichment is within literatures specified range.
   This is a very coarse indicator if its a valid method.
   */
  bool is_within_range;
  
  /** The very coarsely estimated fractional uncertainty on the Pu242 component
   from the correlation estimate; roughly taken from the literature, or some
   wild guess, but I guess better than nothing.  Does not cover additional
   uncerts from using the method outside the literatures specified conditions.
   */
  T pu242_uncert;
};//struct Pu242ByCorrelationOutput


/** Different methods of determining Pu-242 from correlations with the other plutonium nuclides. */
enum class PuCorrMethod : int
{
  /** Pu242 estimator based on correlation with Pu239 only
   
   The fitted power function can be used to predict the 242Pu
   content of the samples with an r.m.s. deviation of 1.2% over the
   range 55–64% 239Pu.
   Tests on other miscellaneous plutonium isotopic compositions
   indicate that the correlation works well (1% relative) for 239Pu
   content less than 70% and moderately well (4% relative) up to
   80% 239Pu
   
   From paper referenced below:
   In FRAM, this function is the equivalent of setting 'Pu242-correlation' value to  9.66E-03 and
   the 'Pu239_exponent' to -3.83 and the other coefficients to zero.
   
   See:
   M.T. Swinhoe, T. Iwamoto, T. Tamura
   'Determination of 242Pu by correlation with 239Pu only',
   Nuclear Instruments and Methods in Physics Research A 615 (2010) 136–137
   https://www.sciencedirect.com/science/article/pii/S0168900210000045#bib2
   */
  ByPu239Only,
  
  /** Bignan 95 Pressurized Water Reactor
   
   Bignan 95, referenced below, uses Pu238/Pu239 and Pu240/Pu239 mass ratios, to estimate
   the Pu242/Pu239 mass ratio.  The Pu238, Pu239, and Pu240 can be determined from gamma
   spectroscopy (Pu242 typically cant be determined by gamma spec), as given by G. Bignan et al 98.
   
   From the paper cited below, the average error on Pu242 evaluation of about 3% for PWR, and
   7% for BWR.
   
   See:
   G. Bignan, W. Ruther, H. Ottmar, A. Schubert, and C. Zimmerman,
   'Plutonium isotopic determination by gamma spectrometry: Recommendations for the 242Pu content
   evaluation using a new algorithm',
   Bulletin ESARDA Nr 28 (1998), 1.
   https://publications.jrc.ec.europa.eu/repository/handle/JRC85594
   */
  Bignan95_PWR,

  /** Same as #PuCorrMethod::Bignan95_PWR, but for Boiling Water Reactors.
   */
  Bignan95_BWR,
  
  // TODO: implement Gunnink80, or values from FRAM (at least for comparison)
  
  /** Not applicable. */
  NotApplicable
};//enum class PuCorrMethod

const std::string &to_str( const PuCorrMethod method );
const std::string &to_description( const PuCorrMethod method );
/** Given the isotopics determined from gamma-spec, returns mass fractions of those
 isotopes, as well as Pu242, as determined from the specified correlation estimate
 method.

 This function is implemented in RelActCalc_imp.hpp, so you will need to include that if you would like to use it.
 */
template<typename T>
Pu242ByCorrelationOutput<T> correct_pu_mass_fractions_for_pu242( Pu242ByCorrelationInput<T> input, PuCorrMethod method );


    

/** Convert a mass ratio to an activity ratio.
 
 E.g., if you want to get the activity ratio of 0.72% enriched U235 to U238, you would:
```
double act_ratio = mass_ratio_to_act_ratio( u235, u238, 0.0072/(1-0.0072) );
```
 */
double mass_ratio_to_act_ratio( const SandiaDecay::Nuclide * const numerator_nuclide, 
                                const SandiaDecay::Nuclide * const denominator_nuclide, 
                                const double mass_ratio );
    
/** A structure to allow inputting self-attenuating, and external attenuating shielding. */
struct PhysicalModelShieldInput
{
  /** The atomic number of the shielding.
    Must be in range [1,98], unless you are defining a material, and then it must be 0.0.
    If you wish to instead specify a material, this number must be 0.0f.
       
    If this is an input, and you are fitting atomic number, then this is the starting value for the fit.
    */
  double atomic_number = 0.0;
      
  /** The shielding material to to use - if non-null, then atomic number must be 0. */
  std::shared_ptr<const Material> material = nullptr;
    
    
  /** The areal density, in units of PhysicalUnits, e.g., `areal_density = 18*PhysicalUnits::g/PhysicalUnits::cm2` *.
    If you are fitting the AD, this will be the starting value for the fit.
  */
  double areal_density = 0.0;
    
  /** If the atomic number of the shielding should be fit.
    Fitting atomic number max of one self-attenuating shielding, and one external shielding is allowed.
       
    Normally false.  Must be false if material is specified.
  */
  bool fit_atomic_number = false;
      
  /** The lower atomic number you want to allow. */
  double lower_fit_atomic_number = 1.0;
  /** The upper atomic number you want to allow. */
  double upper_fit_atomic_number = 98.0;
      
  /** If the shielding thickness should be fit.  Note this is areal density, and not thickness.
    if false, set the AD value by setting BOTH lower and upper AD limit values to be greater than zero, and the same.
    If true, and lower AD and upper AD limits are not equal, then the specified limits will be applied.
      
    Normally true.
  */
  bool fit_areal_density = true;
      
  /** The lower AD for the shielding.  Must be in range [0,500] g/cm2, if fitting AD. */
  double lower_fit_areal_density = 0.0;
      
  /** The upper AD for the shielding.  Must be greater than to lower AD, if fitting, and must be in range [0,500] g/cm2. */
  double upper_fit_areal_density = 0.0;
    
  /** Checks specified constraints are obeyed - throwing an exception if not. */
  void check_valid() const;
  
  static const int sm_xmlSerializationVersion = 0;
  rapidxml::xml_node<char> *toXml( ::rapidxml::xml_node<char> *parent ) const;
  void fromXml( const ::rapidxml::xml_node<char> *parent, MaterialDB *materialDB );
  
  static const double sm_upper_allowed_areal_density_in_g_per_cm2; //Set to 500
  
#if( PERFORM_DEVELOPER_CHECKS )
  static void equalEnough( const PhysicalModelShieldInput &lhs, const PhysicalModelShieldInput &rhs );
#endif
};//struct PhysicalModelShieldInput


/** Since AN ranges from 1 to ~100, we'll scale AN by this amount in the Ceres solver. 
 * Center around 50
 * Aim to have the parameter between0.1 and 10.0.
 * 
*/
const double ns_an_ceres_mult = 50;

/**
 Offsets and multipliers for Hoerl function paramaters, when going from Ceres to value to be evaluated, e.x.: 
  b = (x[b_index] - ns_decay_hoerl_b_offset) * ns_decay_hoerl_b_multiple
  c = (x[c_index] - ns_decay_hoerl_c_offset) * ns_decay_hoerl_c_multiple
 */
const double ns_decay_hoerl_b_offset = 1.0;   
const double ns_decay_hoerl_b_multiple = 1.0;
const double ns_decay_hoerl_c_offset = 0.0;
const double ns_decay_hoerl_c_multiple = 1.0;


template<typename T = double>
struct PhysModelShield
{
  std::shared_ptr<const Material> material;
  T atomic_number{0.0};
  T areal_density{0.0};
};//struct PhysModelShield
  
  
double eval_physical_model_eqn( const double energy,
                               const std::optional<PhysModelShield<double>> &self_atten,
                               const std::vector<PhysModelShield<double>> &external_attens,
                               const DetectorPeakResponse * const drf,
                               std::optional<double> hoerl_b,
                               std::optional<double> hoerl_c );
  
/** Please note, that the
 */
std::function<double(double)> physical_model_eff_function( const std::optional<PhysModelShield<double>> &self_atten,
                                                          const std::vector<PhysModelShield<double>> &external_attens,
                                                          const std::shared_ptr<const DetectorPeakResponse> &drf,
                                                          std::optional<double> hoerl_b,
                                                          std::optional<double> hoerl_c );
  
  
std::string physical_model_rel_eff_eqn_text( const std::optional<PhysModelShield<double>> &self_atten,
                                            const std::vector<PhysModelShield<double>> &external_attens,
                                            const std::shared_ptr<const DetectorPeakResponse> &drf,
                                            std::optional<double> hoerl_b,
                                            std::optional<double> hoerl_c,
                                            const bool html_format );

std::string physical_model_rel_eff_eqn_js_function( const std::optional<PhysModelShield<double>> &self_atten,
                                                   const std::vector<PhysModelShield<double>> &external_attens,
                                                   const DetectorPeakResponse * const drf,
                                                   std::optional<double> hoerl_b,
                                                   std::optional<double> hoerl_c );

/** Refit the continuums for ROIs (peaks grouped by shared continuum) in polynomial-based continuums.
 
 @param solution_peaks The fit peaks from solution.m_fit_peaks_in_spectrums_cal
 @param foreground The foreground spectrum measurement
 @returns A vector of all peaks (both modified and unmodified) sorted by energy
 
 For polynomial-based continuums (excluding NoOffset and External), this function groups peaks
 by their shared continuum and refits the continuum coefficients. If there are more than 4 ROIs,
 the refitting is done using ThreadPool for parallel processing.
 */
std::vector<PeakDef> refit_roi_continuums( const std::vector<PeakDef> &solution_peaks,
                                          std::shared_ptr<const SpecUtils::Measurement> foreground );

}//namespace RelActCalc


#endif //RelActCalc_h
