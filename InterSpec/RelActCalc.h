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

/** Caps the number of threads an individual relative-activity/efficiency solve uses
 internally — both the ROI/atomic-number thread-pools and Ceres' own thread count,
 which are otherwise sized to the hardware concurrency.

 When many solves run concurrently — e.g. the peak-fit optimization GA on a many-core
 machine — that `(outer parallelism) x (hardware_concurrency)` thread count can exhaust
 the OS thread/process limit (`pthread_create` failing with EAGAIN, i.e. "Resource
 temporarily unavailable").

 Pass 0 (the default) for "auto" (use the hardware concurrency); pass a small value
 (e.g. 1) when you provide your own outer parallelism.  Process-wide and thread-safe.
 */
void set_max_solve_threads( const unsigned num_threads );

/** The resolved per-solve thread count: the value given to #set_max_solve_threads
 clamped to [1, hardware_concurrency], or the hardware concurrency when set to 0.
 Always >= 1.
 */
unsigned max_solve_threads();

/** Number of empirical-correction basis coefficients (replaces the Hoerl b,c) - kept at 2 so the
 physical-model parameter-block layout is identical for every correction type. */
constexpr int ns_num_basis_correction_terms = 2;
static_assert( ns_num_basis_correction_terms == 2,
              "The correction always uses exactly 2 parameter slots ('Hoerl b,c'); changing this breaks the"
              " invariant parameter-block layout shared by all correction types." );

/** The single "max window-swing" knob R for the empirical-correction bounds: the largest factor any one
 correction basis term may reshape the relative-efficiency curve across the fit energy window.  Every
 correction form derives symmetric, window-aware bounds about its identity from R, so the start (identity)
 is always interior - never pinned on a bound:
   - Hoerl b   (E_MeV^b):     |b|     <= ln(R)/ln(Ehi/Elo)
   - Hoerl c   (gamma=ln c):  |ln c|  <= ln(R)/(1/Elo_MeV - 1/Ehi_MeV)
   - Chebyshev a1,a2:         |a|     <= ln(R)/2          (basis already normalized to u in [-1,1])
 (with narrow-window clamps b<=2, c in [1/3,3]).  R=10 is a generous guard rail (~7x above the corrections
 actually seen in fits); a 10..30 sweep showed larger R does not improve accuracy and only inflates the
 rel-eff uncertainty band.  This bound is the SOLE guard on the correction now that the coefficient prior
 has been removed. */
constexpr double ns_corr_max_window_swing = 10.0;

/** The empirical correction term applied on top of (DRF x attenuation) in the Physical Model rel-eff.

  - None: no correction.
  - Hoerl: the modified Hoerl function `E_MeV^b * c^(1/E_MeV)` (2 params b,c).
  - Chebyshev: a pivot-anchored exponential of a low-order Chebyshev polynomial in log-energy,
    `correction(E) = exp( sum_{k=1,2} a_k * (T_k(u(E)) - T_k(u(E0))) )`; `u` maps the fit range to [-1,1] and
    the pivot anchoring forces `correction(E0)=1` (adds SHAPE, not scale).  Far better conditioned than the
    Hoerl (whose b,c basis functions are ~96% collinear over a typical range).
 All non-None forms use exactly 2 parameters, so the parameter-block layout is identical. */
enum class PhysModelCorrFcn : int
{
  None,
  Hoerl,
  Chebyshev
};//enum class PhysModelCorrFcn

/** Returns "None"/"Hoerl"/"Chebyshev". */
const char *to_str( const PhysModelCorrFcn corr );
/** Inverse of `to_str`; throws std::runtime_error on an unrecognized string. */
PhysModelCorrFcn phys_model_corr_fcn_from_str( const char *str );

/** A weak Gaussian prior ("biasing"/regularization) applied to a Physical-Model parameter to stabilize the
 near-degenerate shielding/correction covariance.

  - use == false: prior is OFF (no residual added).
  - use == true, value == nullopt: prior ON using the call-site default weight.
  - use == true, value set: prior ON using the specified weight.
 It is an invalid configuration to specify a non-zero `value` while `use` is false (callers throw). */
struct PriorWeightOption
{
  bool use = false;
  std::optional<double> value;

  /** The effective weight: 0.0 if off, else `value` or `default_weight`. */
  double effective_weight( const double default_weight ) const
  {
    return use ? value.value_or(default_weight) : 0.0;
  }

  bool operator==( const PriorWeightOption &rhs ) const
  {
    return (use == rhs.use) && (value == rhs.value);
  }
  bool operator!=( const PriorWeightOption &rhs ) const { return !(*this == rhs); }
};//struct PriorWeightOption

/** Default weights for the Physical-Model areal-density regularization priors (used when a shield's
 `ad_bias` is enabled without an explicit value).  Acts like a Gaussian with sigma = 1/weight in g/cm^2,
 pulling each fit areal density toward 0. */
constexpr double ns_default_ext_atten_ad_prior_weight  = 0.1;   // sigma ~ 10 g/cm^2
constexpr double ns_default_self_atten_ad_prior_weight = 0.05;  // sigma ~ 20 g/cm^2 (self-atten is physical)

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

 Tiny-negative squared uncertainties (numerical noise from the covariance
 matrix) are clamped to zero rather than producing NaN.
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
      
  /** The lower AD for the shielding.  Must be in range [0,500] g/cm2, if fitting AD.

   In units of `PhysicalUnits` - i.e., you need to divide by `PhysicalUnits::g_per_cm2` to printout to g/cm2 human values.
   */
  double lower_fit_areal_density = 0.0;
      
  /** The upper AD for the shielding.  Must be greater than to lower AD, if fitting, and must be in range [0,500] g/cm2.

   In units of `PhysicalUnits` - i.e., you need to divide by `PhysicalUnits::g_per_cm2` to printout to g/cm2 human values.
   */
  double upper_fit_areal_density = 0.0;

  /** Optional weak regularization prior pulling this shield's FIT areal density toward 0 g/cm^2 (off by
   default).  Only has an effect when `fit_areal_density` is true.  Used to break the
   correction-vs-attenuation degeneracy when the data does not strongly constrain the shielding. */
  RelActCalc::PriorWeightOption ad_bias;

  /** Checks specified constraints are obeyed - throwing an exception if not. */
  void check_valid() const;

  static const int sm_xmlSerializationVersion = 0;
  rapidxml::xml_node<char> *toXml( ::rapidxml::xml_node<char> *parent ) const;
  void fromXml( const ::rapidxml::xml_node<char> *parent );
  
  static const double sm_upper_allowed_areal_density_in_g_per_cm2; //Set to 500
  
  bool operator==( const PhysicalModelShieldInput &rhs ) const;
  bool operator!=( const PhysicalModelShieldInput &rhs ) const;

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

/** Areal density (g/cm^2) is rescaled and offset in Ceres parameter space so its magnitude
 is comparable to the other fit parameters, and its value is kept away from zero (Ceres is
 poorly behaved when a parameter sits at exactly zero / at a bound).

   physical_AD[g/cm^2] = (x[ad_index] - ns_ad_ceres_offset) * ns_ad_ceres_mult
   x[ad_index]         =  physical_AD / ns_ad_ceres_mult + ns_ad_ceres_offset

 The common case is AD < ~50 g/cm^2, which with these values lives in ceres [0.5, 2.5].
 */
const double ns_ad_ceres_offset = 0.5;
const double ns_ad_ceres_mult = 25.0;

/** Deviation-pair offsets (keV) rescaled and offset similarly.  Parameter is symmetric and
 may be negative; the offset just keeps the (physical=0) start point away from zero in
 Ceres space.

   physical_offset[keV] = (x[par] - ns_dev_ceres_offset) * ns_dev_ceres_mult
   x[par]               =  physical_offset / ns_dev_ceres_mult + ns_dev_ceres_offset
 */
const double ns_dev_ceres_offset = 0.5;
const double ns_dev_ceres_mult = 3.0;


template<typename T = double>
struct PhysModelShield
{
  std::shared_ptr<const Material> material;
  T atomic_number{0.0};
  T areal_density{0.0};
};//struct PhysModelShield
  
  
/** Evaluates the Physical-Model relative efficiency at `energy`.

 `hoerl_b`/`hoerl_c` carry the 2 correction parameters (Hoerl b,c, or the basis coefficients a1,a2).
 `corr_fcn` selects the correction form.  `corr_lower_energy`/`corr_upper_energy`/`corr_pivot_energy` define
 the energy reference frame for a basis (Chebyshev) correction; they are ignored for None/Hoerl. */
double eval_physical_model_eqn( const double energy,
                               const std::optional<PhysModelShield<double>> &self_atten,
                               const std::vector<PhysModelShield<double>> &external_attens,
                               const DetectorPeakResponse * const drf,
                               std::optional<double> hoerl_b,
                               std::optional<double> hoerl_c,
                               const double corr_lower_energy = 0.0,
                               const double corr_upper_energy = 0.0,
                               const double corr_pivot_energy = 0.0,
                               const PhysModelCorrFcn corr_fcn = PhysModelCorrFcn::Hoerl );

/** Please note, that the
 */
std::function<double(double)> physical_model_eff_function( const std::optional<PhysModelShield<double>> &self_atten,
                                                          const std::vector<PhysModelShield<double>> &external_attens,
                                                          const std::shared_ptr<const DetectorPeakResponse> &drf,
                                                          std::optional<double> hoerl_b,
                                                          std::optional<double> hoerl_c,
                                                          const double corr_lower_energy = 0.0,
                                                          const double corr_upper_energy = 0.0,
                                                          const double corr_pivot_energy = 0.0,
                                                          const PhysModelCorrFcn corr_fcn = PhysModelCorrFcn::Hoerl );


std::string physical_model_rel_eff_eqn_text( const std::optional<PhysModelShield<double>> &self_atten,
                                            const std::vector<PhysModelShield<double>> &external_attens,
                                            const std::shared_ptr<const DetectorPeakResponse> &drf,
                                            std::optional<double> hoerl_b,
                                            std::optional<double> hoerl_c,
                                            const bool html_format,
                                            const double corr_lower_energy = 0.0,
                                            const double corr_upper_energy = 0.0,
                                            const double corr_pivot_energy = 0.0,
                                            const PhysModelCorrFcn corr_fcn = PhysModelCorrFcn::Hoerl );

std::string physical_model_rel_eff_eqn_js_function( const std::optional<PhysModelShield<double>> &self_atten,
                                                   const std::vector<PhysModelShield<double>> &external_attens,
                                                   const DetectorPeakResponse * const drf,
                                                   std::optional<double> hoerl_b,
                                                   std::optional<double> hoerl_c,
                                                   const double corr_lower_energy = 0.0,
                                                   const double corr_upper_energy = 0.0,
                                                   const double corr_pivot_energy = 0.0,
                                                   const PhysModelCorrFcn corr_fcn = PhysModelCorrFcn::Hoerl );

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
