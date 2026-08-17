#ifndef RelActCalcManual_h
#define RelActCalcManual_h
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
#include <vector>
#include <memory>
#include <cstdint>
#include <istream>
#include "ReactionGamma.h"

// Forward declarations
class PeakDef;

class DetectorPeakResponse;

namespace RelActCalc
{
  enum class RelEffEqnForm : int;
  struct PhysicalModelShieldInput;
}

namespace SandiaDecay
{
  struct Nuclide;
  struct Element;
}

namespace SpecUtils
{
  class Measurement;
}

namespace ceres
{
  class Problem;
}

/** The structs/functions in this namespace facilitate performing a relative activity calculation
 on user-fit peaks (hence the term "Manual").
 
 In order to facilitate using nuclear data from both SandiaDecay, and other "generic" sources, there
 are alternative functions to accept input either from a {SandiaDecay::Nuclide,double(age)} pair, or
 from a {std::string(nuclide),double(yield)} source; this is in anticipation of allowing fitting for
 age when the SandiaDecay input method is used (which isnt implemented yet).
 */
namespace RelActCalcManual
{

/** A SandiaDecay defined nuclide source at a fixed age, for use in adding GenericLineInfo to peaks
 from nuclides.
 */
struct SandiaDecayNuc
{
  const SandiaDecay::Nuclide *nuclide = nullptr;
  double age = -1;
  bool correct_for_decay_during_meas = false;
  const SandiaDecay::Element *element = nullptr;
  const ReactionGamma::Reaction *reaction = nullptr;
};//struct SandiaDecayNuc


/** Struct to specify the yield of a nuclide, that will be associated with a specific peak. */
struct GenericLineInfo
{
  /** The number of gammas produced per decay, for this gamma. */
  double m_yield;
  
  /** The name of the parent isotope of the decay chain this gamma is for.
   
   This value is used to group gamma lines together, to figure out the relative amplitudes of each
   nuclides activity.
   
   The actual value of this string is unimportant, and doesnt need to correspond to a nuclide;
   it is only used for grouping informations together.
   
   For example, for the U238 decay chain, a lot of the gammas you will use come from the decay of
   Pa234m in the U238 decay chain - you would label all of these lines as "U238".
   */
  std::string m_isotope;
  
  
  /** Convenience constructor. */
  GenericLineInfo( const double yield, const std::string &isotope );
  
  GenericLineInfo();
};//struct GenericLineInfo


/** Struct to represent a peak. */
struct GenericPeakInfo
{
  /** Energy of peak, in keV
   
   All source gammas in \c m_source_gammas are assumed at this energy.
   
   This value will be set to the energy of the assigned gamma, and is the energy used for clustering
   and calculations.
   */
  double m_energy;
  
  /** This is the fit peak mean.
   
   Only used for plotting.
   */
  double m_mean;
  
  /** The FWHM of the peak; not use for relative activity or efficiency calculations, but useful
   for assigning source gamma terms.
   */
  double m_fwhm;
  
  /** The peak amplitude.
   
   You can also use counts per second for this value, as long as all peaks
   are consistent, and you provide the correct 1-sigma uncertainty of the cps.
   */
  double m_counts;
  
  /** The uncertainty on the peak amplitude. */
  double m_counts_uncert;
  
  /** Additional uncertainty that is independent of the relative uncertainty value,
   and is added in quadrature to the statistical uncertainty (i.e., #m_counts_uncert),
   to get the total uncertainty used for each relative efficiency point.
   
   Must be in range [0,1]
   
   This keeps a very high-statistics peak from forcing the fit relative efficiency curve off of the
   general trend.  A value of 0.5 produces reasonable results for HPGe Uranium spectra; the exact
   value used isnt super important, just that there is a non-zero value.
   
   As a special case, if a value of -1 is specified (and must be specified for every peak),
   an unweighted fit will be performed (i.e., no stat uncert taken into account - all peaks
   contribute equally, despite their size).
   
   Note: replacing this with a `ceres::LossFunction` was tried and rejected - Ceres applies the
   loss to a whole residual block (this problem uses one block for every peak), so it acts on the
   total chi2 and cannot move the minimum.  Per-peak Huber influence was also implemented and
   measured, and moved the answer the wrong way.  See the comment at the `lossfcn` declaration in
   RelActCalcManual.cpp for the measurements and the spectra they were made on.
   */
  double m_base_rel_eff_uncert;
  
  /** The list of gammas that contributes to this peak.
   
   Must be at least one entry, and duplicate entries are not allowed.
   */
  std::vector<GenericLineInfo> m_source_gammas;
  
  
  /** Simple constructor that initializes all the member variables to zero. */
  GenericPeakInfo();
};//struct GenericPeakInfo

/** A constraint on the activity ratio of two nuclides.
 
 This is used to pin the activity ratio of one nuclide to another.
 */
struct ManualActRatioConstraint
{
  std::string m_constrained_nuclide;
  std::string m_controlling_nuclide;
  double m_constrained_to_controlled_activity_ratio;

/*
  static const int sm_xmlSerializationVersion = 0;
  void toXml( ::rapidxml::xml_node<char> *parent ) const;
  void fromXml( const ::rapidxml::xml_node<char> *constraint_node );

#if( PERFORM_DEVELOPER_CHECKS )
  static void equalEnough( const ManualActRatioConstraint &lhs, const ManualActRatioConstraint &rhs );
#endif
*/
};//struct ManualActRatioConstraint


// Mass-fraction constraints solve through a per-element "sigma-block" (an exact reparameterization
//  shared with RelActCalcAuto - see RelActCalc::MassFracBlockSpec in RelActCalc_imp.hpp), and are
//  exercised by the RelActCalcAuto mass-fraction tests (every constrained Auto solve pre-solves
//  through this Manual path).  The flag is kept only to make the feature easy to locate/disable.
#define USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT 1

#if( USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT )
/** A constraint on the mass fraction of a nuclide within its element.

 Fixed (lower == upper) constraints pin the fraction exactly; range constraints fit it within
 [lower, upper].  An element may have every one of its nuclides constrained, in which case the
 windows must be able to sum to exactly 1, and the elements total (relative mass) is fit instead
 of any individual activity.
 */
struct MassFractionConstraint
{
  std::string m_nuclide;
  double m_mass_fraction_lower;
  double m_mass_fraction_upper;

  /** You must specify the specific activities of each nuclide for this same element,
   including for this nuclide itself.
   */
  std::map<std::string, double> m_specific_activities;
};//struct MassFractionConstraint
#endif //USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT

/** Adds the `GenericLineInfo` info (e.g. nuclides and their BR) to input `peaks` by clustering gamma lines of
 provided nuclides.
 
 @param peaks Input peaks, with all info except `GenericPeakInfo::m_source_gammas` filled out.
        Note: there must not be any duplicate energies, or else an exception will be thrown, as assignment would be ambigious..
 @param nuclides The input nuclides to cluster and assign to peaks
 @param real_time The real time of the measurement - only used if correcting for the nuclides decay during measurement
        \sa `SandiaDecayNuc::correct_for_decay_during_meas`
 @param cluster_sigma The number of peak sigma to cluster gamma energies as being responsible for the peaks.
 */
std::vector<GenericPeakInfo> add_nuclides_to_peaks( const std::vector<GenericPeakInfo> &peaks,
                                                   const std::vector<SandiaDecayNuc> &nuclides,
                                                   const double real_time,
                                                   const double cluster_sigma = 1.5
                                                   );


/** Uses weighted least squares (i.e., matrix based solution) to determine relative efficiency
 equation parameters and uncertainties.
 
 The resulting relative efficiency curve will be "pinned", within uncertainties, so the
 lowest-energy will have a relative efficiency of 1.0.
 
 Note that for RelEffEqnForm::LnY and RelEffEqnForm::LnXLnY, the log of each side is taken before
 solving for the parameters, meaning the uncertainty estimates become skewed (using a L-M based
 fit was implemented to compare against, but is not currently included, as it is basically just
 duplicate code).
 
 Throws exception on input or solution error.
 
 @par fcn_form The form of equation to fit coefficients for.
 @par order The order (e.g., number of energy dependent terms) the relative efficiency equation
      should have.  The resulting \c fit_pars will have one more entry than this value.
 @par isotopes Names of the nuclides that will contribute to the peaks; must exactly match nuclides
      specified in the \c peak_infos peaks; string value has no meaning, just used to group
      activities across multiple peaks.
 @par peak_infos The input peaks information (e.g., peaks from PeakEasy/InterSpec/Gadras with
      nuclide info associated with them).
 @par fit_pars [out] The fit coefficients will be placed into this vector; will have \c order+1
      entries; these coefficients are what you will pass into #eval_eqn to get the relative
      efficiency of an energy.
 @par covariance [out] If non-nullptr, the covariance matrix will be placed into this vector of
      vectors; result will be of size order by order.
      Note that the full covariance matrix is required to estimate the function uncertainty as
      assuming coefficients are uncorrelated leads to wildly crazy uncertainties.
      Also, for RelEffEqnForm::LnY, RelEffEqnForm::LnXLnY, and RelEffEqnForm::FramEmpirical the
      covariance matrix is for the log of the equation.
 */
void fit_rel_eff_eqn_lls( const RelActCalc::RelEffEqnForm fcn_form,
                           const size_t order,
                           const std::vector<std::string> &isotopes,
                           const std::vector<double> &rel_acts,
                           const std::vector<GenericPeakInfo> &peak_infos,
                           std::vector<double> &fit_pars,
                           std::vector<std::vector<double>> *covariance );

/** A convenience method for calling the above using SandiaDecay defined nuclides.
 
 Throws exception on input or solution error.
 
 @par fcn_form The form of equation to fit coefficients for.
 @par order The order (e.g., number of energy dependent terms) the rell eff equation should have.
 @par nuclides The nuclides that will contribute to the peaks.
 @par base_rel_eff_uncert This is an additional uncertainty added in quadrature to the statistical
      uncertainty of the relative efficiency equation.  It is treated as uncorrelated between each
      peak.  The purpose is to keep high-stat peaks from dominating the rel eff curve away from the
      general trend.  This variable must be in the range [0,1], or the special value of -1.
      A reasonable value seems to be 0.5, with the exact value not mattering much, just as long as
      there is something.  Also, the special value of -1.0 will cause the fit to be totally
      un-weighted (e.g., each peak contribute same to fit, no matter their statistical uncertainty).
 @par peak_infos The input peaks information
 @par fit_pars [out] The fit coefficients will be placed into this vector; see notes above.
 @par covariance [out] If non-nullptr, the covariance matrix will be placed into this vector of
      vectors; see notes above
 */
/*
void fit_rel_eff_eqn_lls( const RelActCalc::RelEffEqnForm fcn_form,
                           const size_t order,
                           const std::vector<SandiaDecayNucRelAct<double>> &nuclides,
                           const double base_rel_eff_uncert,
                           const std::vector<std::shared_ptr<const PeakDef>> &peak_infos,
                           std::vector<double> &fit_pars,
                           std::vector<std::vector<double>> *covariance );
*/

/** Function to do actual LLS work.  Split out from the above for testing purposes.  */
void fit_rel_eff_eqn_lls( const RelActCalc::RelEffEqnForm fcn_form,
                         const size_t order,
                         const std::vector<double> &energies,
                         const std::vector<double> &data_values,
                         const std::vector<double> &data_uncertainties,
                         std::vector<double> &fit_pars,
                         std::vector<std::vector<double>> *covariance );


/** For a given relative efficiency curve, will fit the relative activities of the isotopes to it.
 
 \param fit_rel_acts The relative activities of the isotopes. Corresponds to \p isotopes on an index-by-index basis.
 
 */
void fit_act_to_rel_eff( const std::function<double(double)> &eff_fcn,
                        const std::vector<std::string> &isotopes,
                        const std::vector<GenericPeakInfo> &peak_infos,
                        std::vector<double> &fit_rel_acts,
                        std::vector<double> &fit_rel_act_uncerts );

/** convenience method for calling above `fit_act_to_rel_eff(...)`. 
 * @param eqn_form The form of the relative efficiency equation to use.
 *                  May not be `RelActCalc::RelEffEqnForm::FramPhysicalModel`.
 */
void fit_act_to_rel_eff( const RelActCalc::RelEffEqnForm eqn_form,
                        const std::vector<double> &eqn_pars,
                        const std::vector<std::string> &isotopes,
                        const std::vector<GenericPeakInfo> &peak_infos,
                        std::vector<double> &fit_rel_acts,
                        std::vector<double> &fit_rel_act_uncerts );

struct RelEffInput;  //Forward declaration
/** Fits the activities, for a given PhysicalModel Relative Efficiency equation.
 * 
 * The input `RelEffInput` must have valid DRF and self attenuation set.  Parameter values
 * are determined similar to the starting values used in the full Relative Activity fit.
 */
void fit_act_to_phys_rel_eff( const RelEffInput &input,
                        const std::vector<std::string> &isotopes,
                        const std::vector<GenericPeakInfo> &peak_infos,
                        std::vector<double> &fit_rel_acts,
                        std::vector<double> &fit_rel_act_uncerts );


/** A simple struct to hold the determined relative activity of a isotope. */
struct IsotopeRelativeActivity
{
  /** The isotope this information is for; corresponds to a value provided by #GenericLineInfo. */
  std::string m_isotope;
  
  /** The relative activity for the isotope.  That is, this value, times relative efficiency (and
   possible times live-time if the peaks you passed were in CPS instead of raw counts), times
   the gammas yield, will give you the predicted peak amplitude.
   */
  double m_rel_activity;
  
  /** The uncertainty of the activity; this value is just the sqrt of the diagonal of the relative
   activity covariance matrix.  I.e., it doesnt account for uncertainty of the relative efficiency;
   correlations between isotope activities could also be important, and not taken into account here.
   */
  double m_rel_activity_uncert;
};//struct IsotopeRelativeActivity


/** The input for the manual relative efficiency calculation. */
struct RelEffInput
{
  /** The required input peak information.
 
   Each peak must have at least one source isotope, have counts above zero, have counts
   uncertainty greater than zero (or equal to -1.0), and each isotope must have yields
   at each energy greater than 0.
   
   If you are decay-correcting yields during measurement, the `GenericPeakInfo::m_source_gammas::m_yield`
   values should have already had this applied.
   */
  std::vector<GenericPeakInfo> peaks;
  
  /** If you are decay-correcting peaks, then the non-decay-corrected version of the peaks
   are put here.   Only peaks with corrections should be put here, and should only include info
   for only gamma lines that have been corrected.
   
   Only includes peaks that have at least one correction, and the only entries in
   `GenericPeakInfo::m_source_gammas`, are the ones with corrections.
   
   The peaks are _not_ used during computation, and placed here only so
   you can list effects of decay-correction.  Not super clean, but maybe the best way
   to carry all the information together.
   */
  std::vector<GenericPeakInfo> peaks_before_decay_correction;
  
  /** Any warnings or notes you encountered while preparing input.
   
   These are not used for computation, but will be copied into `RelEffSolution::m_warnings`.
   */
  std::vector<std::string> prep_warnings;
  
  /** The equation form to use; i.e., EqnForm::LnX, EqnForm::LnY, or EqnForm::LnXLn */
  RelActCalc::RelEffEqnForm eqn_form;

  /** The order of equation to be fit for.  i.e., the number of energy-dependent terms
  to be fit for (e.g., the total number of terms will be this plus one). 
  
  If eqn_form is `RelActCalc::RelEffEqnForm::FramPhysicalModel`, then this value must
  be 0, as it is not used.
  */
  size_t eqn_order;

  /** If true, use Ceres to fit the relative efficiency equation.
   * If false, use LLS to fit the relative efficiency equation.
   *
   * For `RelActCalc::RelEffEqnForm::FramPhysicalModel`, this must be true.
   */
  bool use_ceres_to_fit_eqn = false;

  /** Skip the covariance / Jacobian / derived-uncertainty computations: \p m_nonlin_covariance,
   \p m_rel_act_covariance / \p m_rel_act_jacobian and \p m_nonlin_jacobian come back empty, and
   the relative-activity uncertainties are the -1 "not available" sentinel.
   Point estimates, chi2, DOF, and shield fit values are still produced.
   (In the LLS fit mode \p m_rel_eff_eqn_covariance is still filled - it falls out of the same
   linear fit that produces the curve coefficients, so there is nothing to skip.)
   Used by the atomic-number-scan multi-start and by profile-likelihood sub-solves, where only
   the chi2 of the re-fit is needed.
   */
  bool point_estimate_only = false;

  /** Skip the SCAN_AN_FOR_BEST_FIT multi-start over self-attenuator atomic number (AN is still
   fit locally when `fit_atomic_number` is set).  Used by profile-likelihood sub-solves, which
   warm-start near the nominal solution.
   */
  bool skip_an_scan = false;

  /** If true, widen every peak's statistical uncertainty by a single common factor estimated
   from how far the peaks actually sit from the fitted model - see #estimate_stat_uncert_multiple.

   Because every peak is scaled by the SAME factor, the 1/k^2 factors straight out of the
   weighted-least-squares objective: the fitted values cannot move, only the uncertainties widen.
   That is the point of this form.  Widening peaks by differing amounts - which is what a
   fractional-of-peak-area term such as #GenericPeakInfo::m_base_rel_eff_uncert does - re-weights
   the peaks against each other and shifts the answer (see RelActManualGui::AddUncert for the
   measured size of that shift).

   The estimate is outlier-insensitive (median based), so a few badly mis-fit peaks do not drive
   it.  The factor used is reported in RelEffSolution::m_auto_stat_uncert_multiple.
   Ignored for unweighted fits.
   */
  bool auto_estimate_add_uncert = false;

  /** If true (the default), the reported uncertainties are widened when the peaks scatter about
   the fitted model by more than their uncertainties allow - see RelEffSolution::m_cov_scale.

   Set false to get the purely statistical uncertainty, i.e. what the covariance says if the model
   is taken to be exactly right.  That is the correct quantity to quote only when the fit really
   does describe the data within its uncertainties; on real uranium spectra it is typically far
   too small (measured on CBNM446: +-0.04 wt% reported while sitting 0.58 wt% from the certified
   value), because the peaks disagree with any achievable relative-efficiency curve by much more
   than counting statistics.
   */
  bool widen_uncerts_for_scatter = true;


  /** If true, fit the modified Hoerl equation form for the physical model.
   * If false, do not fit the modified Hoerl equation form (its value will be constant value of 1.0).
   * 
   * Ignored if not using `RelActCalc::RelEffEqnForm::FramPhysicalModel`.
  */
  bool phys_model_use_hoerl = true;

  /** The detector peak response to use as part of the relative efficiency equation for 
   `RelActCalc::RelEffEqnForm::FramPhysicalModel`, and if so, the detector must be non-null 
   and valid.
   Not used for other equation forms.
  */
  std::shared_ptr<const DetectorPeakResponse> phys_model_detector;
  
  /** The self attenuation if the equation form is `RelActCalc::RelEffEqnForm::FramPhysicalModel`.
    If specified for any other equation form, will throw an exception.
  */
  std::shared_ptr<const RelActCalc::PhysicalModelShieldInput> phys_model_self_atten;
  
  /** The external attenuations for equation form `RelActCalc::RelEffEqnForm::FramPhysicalModel`. 
    If specified for any other equation form, will throw an exception.
  */
  std::vector<std::shared_ptr<const RelActCalc::PhysicalModelShieldInput>> phys_model_external_attens;

  /** The activity ratio constraints. */
  std::vector<ManualActRatioConstraint> act_ratio_constraints;

#if( USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT )
  /** The mass fraction constraints. 
   
   Note: this initial implementation is less than ideal, as the initial relative activities
   of the non-mass-fraction constrained isotopes are determined by just ignoring the 
   mass-fraction-constrained isotopes exist at all.  This can lead to say the 185 keV peak
   being attributed to the tiny (BR=1.3E-4) peak of U238 at 184.8 keV, instead of the
   main peak of U235.
  */
  std::vector<MassFractionConstraint> mass_fraction_constraints;
#endif //USE_REL_ACT_MANUAL_MASS_FRACTION_CONSTRAINT
  /** Checks that the nuclide constraints are valid.

   Checks for cyclical constraints, and that all constrained nuclides are found in #nuclides.
   
   Throws an exception if they are not valid.
   */
  void check_nuclide_constraints() const; 
};//struct RelEffInput

/** The status of fitting for a solution. */
enum class ManualSolutionStatus : int
{
  NotInitialized,
  ErrorInitializing,
  ErrorFindingSolution,
  ErrorGettingSolution,
  Success
};//enum class ManSolutionStatus


/** Options for a profile-likelihood scan of one nuclide's mass fraction - see
 #profile_mass_fraction.
 */
struct ProfileMassFractionOptions
{
  /** The nuclide whose (element-relative) mass fraction to profile; must resolve to a
   SandiaDecay nuclide present in the problem's peaks.
   */
  std::string nuclide;

  /** The confidence levels to find interval end-points for, as CENTRAL two-sided intervals;
   each maps to a delta-chi2 threshold of quantile(chi_squared(1), CL).
   Defaults to 1-sigma and 2-sigma.
   */
  std::vector<double> confidence_levels{ 0.682689492, 0.954499736 };

  /** Maximum bisection solves per interval side, per confidence level. */
  size_t max_solves_per_side = 12;

  /** Absolute tolerance on the mass fraction at each interval crossing. */
  double frac_tolerance = 1.0e-4;
};//struct ProfileMassFractionOptions


/** One profile-likelihood interval, at a single confidence level. */
struct ProfileMassFractionInterval
{
  double confidence_level = 0.0;

  /** The chi2 rise defining the interval end-points: quantile(chi_squared(1), CL), times the
   nominal solution's max(1, m_cov_scale) - so the profile carries the same chi2/dof error
   inflation the covariance-based uncertainties do, and the two agree in the Gaussian limit.
   */
  double delta_chi2 = 0.0;

  /** Interval end-points, as ELEMENT-relative mass fractions (the constraint-native quantity;
   equals RelEffSolution::mass_fraction(nuclide) for single-element problems).
   */
  double lower_frac = 0.0;
  double upper_frac = 0.0;

  /** Set when the corresponding end-point ran into the physical [0,1] boundary (or a
   pre-existing mass-fraction constraint window edge) before crossing the delta-chi2 threshold;
   the end-point then holds the boundary value.
   */
  bool lower_at_bound = false;
  bool upper_at_bound = false;
};//struct ProfileMassFractionInterval


/** Result of a profile-likelihood scan of one nuclide's mass fraction - see
 #profile_mass_fraction.
 */
struct ProfileMassFractionResult
{
  std::string nuclide;

  /** The nominal ELEMENT-relative mass fraction (fraction of the summed mass of the nuclides
   in the constraint roster / element), from the nominal solution.
   */
  double nominal_mass_fraction = 0.0;

  /** The fit-weight chi2 (see RelEffSolution::m_chi2_fit_weights) the delta-chi2 is referenced
   to; normally the nominal solution's, but re-anchored to the scan minimum if the profile
   found a better fit (a warning is added in that case).
   */
  double nominal_chi2 = 0.0;

  /** One interval per requested confidence level, in the requested order. */
  std::vector<ProfileMassFractionInterval> intervals;

  /** The (element-relative fraction, fit-weight chi2) points evaluated during the scan, sorted
   by fraction - useful for plotting or diagnosing the profile shape.
   */
  std::vector<std::pair<double,double>> scan_points;

  std::vector<std::string> warnings;
};//struct ProfileMassFractionResult

/** A struct to hold the information about the solution to fitting the relative activities and
 efficiency curves.
 */
struct RelEffSolution
{
  /** The input used to compute this solution.
   
   Note that if you fit the atomic number of a physical model shielding, this values in this
   struct may not exactly match your input, because we are currently scanning the AN
   range to get the an.  See portions of the source code marked by
   `SCAN_AN_FOR_BEST_FIT`.
   */
  RelEffInput m_input;

  std::vector<double> m_rel_eff_eqn_coefficients;

  /** The covariance matrix of the relative efficiency equation coefficients.

   For the empirical equation forms this is always the covariance CONDITIONAL on the fitted
   relative activities (from the linear-least-squares fit of the curve to the measured rel. eff.
   points), no matter which method fit the coefficients.  This is deliberate: for these forms
   only the product (curve x activities) is determined by the data, so the marginal coefficient
   covariance is dominated by that unidentifiable trade - and the measured rel. eff. points the
   band is drawn against are themselves computed with the fitted activities, so they move with
   the curve under it.

   For `RelActCalc::RelEffEqnForm::FramPhysicalModel` this is the rel-eff sub-block of the full
   Ceres parameter covariance (marginal); there the DRF sets the absolute scale, so the shield
   and Hoerl parameters are identifiable in their own right.

   Includes the \p m_cov_scale inflation.
  */
  std::vector<std::vector<double>> m_rel_eff_eqn_covariance;
  
  /** The relative activities for each of the input nuclides. */
  std::vector<IsotopeRelativeActivity> m_rel_activities;

  /** The parameters fit by Ceres.

   Note: kept in the raw solver gauge - for Ceres-fit empirical forms the reported
   \p m_rel_activities and \p m_rel_eff_eqn_coefficients have additionally been re-expressed in
   the LLS-mode normalization convention (average measured rel. eff. == 1), so they will not
   reproduce exactly from these raw parameters.
   */
  std::vector<double> m_fit_parameters;

  /** Covariance matrix of nonlinear parameters fit by Ceres, in raw (Ceres) parameter space.

   If not empty, the first `m_rel_activities.size()` indices are the activity multiple
   parameters, in the same index ordering as \p m_rel_activities.
   But also see `m_activity_norms`, as you need to multiple the relative activities by these
   scale factors before using with this covariance matrix; for act-ratio controlled or
   mass-fraction constrained nuclides the parameter slot is not a simple activity multiple at
   all (constant sentinel, or a sigma-block carrier/distribution slot in [0.5,1.5]) - use
   \p m_rel_act_covariance for anything expressed in relative activities.
   (Prior to 20250816 the rows/columns of range-constrained mass-fraction nuclides were
   pre-scaled by a diagonal d(RelAct)/d(par) derivative; they no longer are.)

   If the equation form is `RelActCalc::RelEffEqnForm::PhysicalModel`, then the following
   indeices are for the shielding parameters:
   - {self-atten AN}
   - {self-atten AD}
   - {external-atten 0 AN} (if >=1 external attens specified)
   - {external-atten 0 AD} (if >=1 external attens specified)
   - {external-atten 1 AN} (if >=2 external attens specified)
   - {external-atten 1 AD} (if >=2 external attens specified)
   - ...
   - {Modified Hoerl b}
   - {Modified Hoerl c}
   */
  std::vector<std::vector<double>> m_nonlin_covariance;

  /** Jacobian of the final reported relative activities with respect to the fit parameters,
   `m_rel_act_jacobian[i][j] = d(RelAct_i)/d(par_j)`, computed by automatic differentiation
   through the cost functor at the solution.  Rows index \p m_rel_activities; columns index
   \p m_fit_parameters.  These are derivatives of the FINAL relative activities - activity
   norms, act-ratio constraint chains, and mass-fraction block decoding are all included -
   so no \p m_activity_norms scaling applies.  Empty if covariance computation failed.
   */
  std::vector<std::vector<double>> m_rel_act_jacobian;

  /** Covariance of the final reported relative activities:
   `m_rel_act_jacobian * m_nonlin_covariance * m_rel_act_jacobian^T`.
   Same index ordering as \p m_rel_activities.  This is the matrix every derived uncertainty
   (per-isotope activity sigmas, activity ratios, mass-fraction variations) is computed from.
   Empty if covariance computation failed.
   */
  std::vector<std::vector<double>> m_rel_act_covariance;


  /** The Jacobian matrix of the nonlinear parameters fit by Ceres.
  
  i.e, `m_nonlin_jacobian[k][i] = d residual[k] / d parameters[i]`
  
  To access the Jacobian for the k'th residual and i'th parameter, 
  use `m_nonlin_jacobian[k][i]`.  The k-index cooresponds to the index 
  of the peak in `m_input.peaks`.  The i-index cooresponds to the index 
  of the parameter in `m_fit_parameters`.  You might want to think of it as
  `m_nonlin_jacobian[peak][parameter]`.

  Note: if the compile time option `USE_RESIDUAL_TO_BREAK_DEGENERACY` is true, 
  then there will be one more residual than the number of peaks.
  */
  std::vector<std::vector<double>> m_nonlin_jacobian;

  /** The Chi2 summed over all peaks between their actual and fit relative efficiencies.
   
   Note that this always uses the peaks statistical uncertainties, and does not include
   #GenericPeakInfo::m_base_rel_eff_uncert - see \p m_chi2_fit_weights for the chi2 under the
   weights the fit actually minimized.
   */
  double m_chi2 = 0.0;

  /** The chi2 under the weights the fit actually minimized - i.e., including any
   #GenericPeakInfo::m_base_rel_eff_uncert - equal to 2x the final Ceres cost.
   Compare \p m_chi2, which is statistical-only and intended for display.
   This is the chi2 used for the \p m_cov_scale covariance inflation, and the one a
   profile-likelihood scan should difference.
   Will be -1.0 if the solve did not complete.
   */
  double m_chi2_fit_weights = -1.0;

  /** Variance inflation applied (once) to the stored covariances - and hence to every derived
   uncertainty - when the data scatter about the model exceeds the assumed uncertainties.
   1.0 when m_dof <= 0, for unweighted fits, or on covariance failure; never less than 1.

   This uses an OUTLIER-INSENSITIVE estimate of the scatter,
   `median(pull^2) * (num_peaks/m_dof) / 0.4549`, rather than the usual chi2/dof: the excess in
   these fits is typically a few individually mis-fit peaks (continuum/skew/interference) that
   land many sigma out without much affecting the composition, and the mean of the squared
   pulls would let them inflate every reported uncertainty.  For a broad excess the two agree.

   Divide the covariances by this to recover the raw Ceres values.
   \sa RelActCalcAuto::RelActAutoSolution::m_cov_scale, which uses the chi2/dof form.
   */
  double m_cov_scale = 1.0;

  /** When RelEffInput::auto_estimate_add_uncert was used, the common multiple `k` that every
   peak's statistical uncertainty was scaled by (>= 1).  -1.0 when it was not used.

   Note this leaves the fitted values untouched - it only widens uncertainties - so the reported
   chi2/dof lands near 1 by construction and \p m_cov_scale has nothing left to inflate.
   */
  double m_auto_stat_uncert_multiple = -1.0;

  /** The number of degrees of freedom in the fit: the number of peaks minus the effective
   number of fitted parameters.

   The effective parameter count comes from the solver bookkeeping: all Ceres parameters, minus
   parameters held constant (act-ratio controlled and fixed mass-fraction activity slots,
   non-fit shield parameters, and - for Ceres-fit empirical forms - the gauge-pinned curve
   coefficient), plus, for the LLS fit mode, the (eqn_order + 1) profiled curve coefficients
   less the one combination fixed by the average-measured-rel-eff normalization convention.

   So with no constraints this is num_peaks - num_isotopes - eqn_order for the empirical forms,
   and num_peaks - num_isotopes - num_fit_shield_pars - (use_hoerl ? 2 : 0) for the Physical
   Model (which has no normalization degeneracy, since the DRF provides an absolute scale).

   Can legitimately be 0 - guard before dividing by this.
   */
  int m_dof = 0;
  
  /** The number of evaluation calls it took L-M to reach a solution.
   Only useful for debugging and curiosity.
   */
  int m_num_function_eval_solution = 0;
  
  /** The number of evaluation calls it took to reach a solution, and compute final covariance. */
  int m_num_function_eval_total = 0;
  
  /** How long it took to compute the answer (only for curiosity).
   (64-bit, since an int of microseconds overflows past ~35 minutes)
   */
  std::int64_t m_num_microseconds_eval = 0;
  
  /** As an internal detail of fitting the relative efficiencies, we normalize the activities to a flat relative efficiency curve of 1.0, and then
   fit for the multiple of the normalization that yields the best solution.  This member variable keeps track of these normalizations; we
   are keeping them around to help in computing the correlation compensated ratios (although we could instead modify
   m_nonlin_covariance - but we'll just keep an extra variable around to make things a little easier to debug).

   The entries in this vector correspond to \p m_rel_activities on an index-by-index basis.

   To state it another way, these are the relative activities if the relative efficiency curve is 1.0.

   If a nuclide is controlled by another nuclide, its value will be -1.
   If a nuclide is mass-fraction constrained (either fixed to a specific value, or constrained to a
   range), its value in this entry will be 1.0.

   NOTE: for Ceres-fit empirical forms the reported \p m_rel_activities have additionally been
   re-expressed in the LLS-mode gauge, so `m_rel_activities[i] / m_activity_norms[i]` does NOT
   recover the fitted parameter there (it is off by that gauge multiple), and neither does
   combining these norms with \p m_nonlin_covariance.  Use \p m_rel_act_covariance for anything
   expressed in relative activities.
   */
  std::vector<double> m_activity_norms;
  
  /** The status of creating the solution. */
  ManualSolutionStatus m_status = ManualSolutionStatus::NotInitialized;
  
  /** Error message if computation wasn't successful. */
  std::string m_error_message;
  
  /** Warnings about the setup or solution of the problem; by no means comprehensive of potential
   issues!

   Note: contents of `RelEffInput::prep_warnings` are copied into this variable.
   */
  std::vector<std::string> m_warnings;

  /** Optional profile-likelihood mass-fraction results.

   NOT filled by #solve_relative_efficiency - the caller (GUI background job, LLM tool, tests)
   runs #profile_mass_fraction on demand and stores the results here so the reporting surfaces
   (chart title, mass-fraction table, HTML report, JSON) can pick them up.
   */
  std::vector<ProfileMassFractionResult> m_profile_mass_fractions;



  /** A struct to hold the self attenuation shield fit results. 
   * This is fine for simple accesses, but not if you need to take into account the correlations, which 
   * you really should do a lot of the time.
  */
  struct PhysModelShieldFit
  {
    std::shared_ptr<const Material> m_material;
    double m_areal_density = 0.0;
    double m_areal_density_uncert = -1.0;
    double m_atomic_number = 0.0;
    double m_atomic_number_uncert = -1.0;
  };//struct PhysModelShieldFit
  std::unique_ptr<PhysModelShieldFit> m_phys_model_self_atten_shield;
  std::vector<std::unique_ptr<PhysModelShieldFit>> m_phys_model_external_atten_shields;

  /////////////////////////////////////////////////////////////////////////////////////////////////
  //            Below here are member functions for simplified access to information             //
  /////////////////////////////////////////////////////////////////////////////////////////////////
  
  /** Returns the index in #m_rel_activities for the specified nuclide.
   
   Throws std::exception if an invalid nuclide.
   */
  size_t nuclide_index( const std::string &nuc ) const;

  /** Walks the constraints to find the controlling nuclide for the specified nuclide.
   * 
   * @param iso_index The index of the nuclide to walk to the controlling nuclide of.
   * @param multiple The multiple to multiply the activity by to get to the controlling nuclide.
   * 
   * @return True if a controlling nuclide was found, false otherwise.
   * 
   * \sa RelEffInput::act_ratio_constraints
   */
  bool walk_to_controlling_nuclide( size_t &iso_index, double &multiple ) const;
  
  /** The relative activity of a nuclide.
   
   Note that this assumes the peak counts are in CPS (or equivalently the data is 1-second long),
   so you will need to divide by live-time to compare to the "auto" relative activity.
   
   Throws std::exception if an invalid nuclide.
   */
  double relative_activity( const std::string &nuclide ) const;

  /** Returns the uncertainty of a nuclide.

      Will throw exception if covariances were not computed, or invlaid nuclide.
   */
  double relative_activity_uncertainty( const std::string &nuclide ) const;

  /** The fit relative efficiency curve value; the curve is shifted so its centered around 1
   over all your input points.
   */
  double relative_efficiency( const double energy ) const;
  

  /** Returns the activity ratio between the two isotopes at index \p iso1 and index \p iso2, where the
   indexes are int #RelEffSolution::m_rel_activities.
   */
  double activity_ratio( const size_t iso1, const size_t iso2 ) const;
  
  /** A convenience method for the above #ratio function,
   
   If either isotope is invalid, will throw std::exception.
   */
  double activity_ratio( const std::string &iso1, const std::string &iso2 ) const;
  
  /** Returns the activity ratio uncertainty between the two isotopes at index \p iso1 and index
   \p iso2, taking into account correlations of the relative activities.
   
   Will throw exception if covariances were not computed, or invalid indexes.
   
   Note: it appears that taking into account correlations usually makes the uncertainty _smaller_ than not taking them into account.
   */
  double activity_ratio_uncert( const size_t iso1, const size_t iso2 ) const;
  
  /** A convenience method for the above #ratio_uncert function,
   
   If either isotope is invalid, or covariances not computed, will throw std::exception.
   */
  double activity_ratio_uncert( const std::string &iso1, const std::string &iso2 ) const;
  
  /** Returns the mass fraction of the specified nuclide.
   
   Will throw exception if invalid nuclide name (e.g., a reaction), or negative mass fraction.
   */
  double mass_fraction( const std::string &iso ) const;
  
  /** Returns the mass fraction of the specified nuclide, with it varied the specified number of sigma away.
   
   @param num_sigma The number of sigma away from nominal, to vary the nuclide in questions activity.
   
   Should be properly taking into account the covariance matrix for relative activities (but not rel-eff curve).
   
   Will throw exception if invalid nuclide name (e.g., a reaction), or negative mass fraction.
   */
  double mass_fraction( const std::string &iso, const double num_sigma ) const;
  
  /** Returns a short parameter description.
   
   Same indexing as `m_fit_parameters`; throws exception if invalid index.
   */
  std::string parameter_name( const size_t par_num ) const;
  
  /** Prints out a summary of the results to the provided stream; for development/debug. */
  std::ostream &print_summary( std::ostream &strm ) const;
  
  /** Creates a self-contained HTML report of the results.
   
   @param strm The stream to place the HTML file into.
   @param spectrum_title The title to display on the HTML pace
   @param spectrum The optional spectrum to display on the HTML page (may be nullptr)
   @param spectrum_display_peaks The peaks to display on the spectrum; may be empty.
   @param background The optional background spectrum to display on the HTML page (may be nullptr)
   @param background_normalization The background normalization; if zero or negative, and a background
          is provided, then spectrum live-times will be used for normalization.
   */
  void print_html_report( std::ostream &strm,
                         std::string spectrum_title,
                         std::shared_ptr<const SpecUtils::Measurement> spectrum,
                         std::vector<std::shared_ptr<const PeakDef>> spectrum_display_peaks,
                         std::shared_ptr<const SpecUtils::Measurement> background,
                         double background_normalization,
                         std::vector<std::shared_ptr<const PeakDef>> background_peaks
                         ) const;
  
  /** Makes a HTML table of the activity and mass fractions of all the nuclides.
   
   The table has CSS style class "nuctable resulttable" and columns: "Nuclide", "Rel. Act",
   and "Mass Fraction".
   */
  void get_mass_fraction_table( std::ostream &strm ) const;
  
  /** Makes a HTML table of the activity and mass ratios for all the nuclides.
   
   The table has CSS style class "massratiotable resulttable" and columns: "Nuclide", "Mass Ratio",
   and "Activity Ratio".
   */
  void get_mass_ratio_table( std::ostream &strm ) const;

  /** Writes a short HTML summary of the fitted Physical Model shieldings - the material or
   atomic number, and the areal density, each with its 1-sigma uncertainty when available.
   Writes nothing for non-physical-model solutions, or when no shieldings were used.
   Used by both the GUI results tab and print_html_report().
   */
  void get_phys_model_shield_text( std::ostream &strm ) const;
  
  /** Returns the value of the relative efficiency equation at the specified energy.
   */
  double rel_eff_eqn_value( const double energy ) const;

  /** Returns the uncertainty of the relative efficiency equation at the specified energy.
    
     throws std::exception if covariances are not available, or there is another error.
   */
  double rel_eff_eqn_uncert( const double energy ) const;

  /** Returns a JavaScript function that returns the uncertainty of the relative efficiency equation at the specified energy.
   * 
   * The function interpolates between values computed by the C++.
   * 
   * If couldnt evaluate the function, returns "null".
   */
  std::string rel_eff_eqn_js_uncert_fcn() const;


  /** Returns the equation text.
   
   @param html_format Currently only applicable to Physical Model.
   
   This is a convenience wrapper over `RelActCalc::rel_eff_eqn_text(...)` and `RelActCalc::physical_model_rel_eff_eqn_text(...)`
   */
  std::string rel_eff_eqn_txt( const bool html_format ) const;
  
  /** A convenience wrapper over `RelActCalc::rel_eff_eqn_js_function(...)` and
   `RelActCalc::physical_model_rel_eff_eqn_js_function(...)`
   */
  std::string rel_eff_eqn_js_function() const;
};//struct RelEffSolution


/** Solve for the relative efficiency equation and relative activities for all isotopes
 
 @param input The input peaks and options for the calculation.
 @returns The solution to the problem.  Make sure to check RelEffSolution::m_status,
          RelEffSolution::m_warnings, and RelEffSolution::m_error_message.
 */
RelEffSolution solve_relative_efficiency( const RelEffInput &input );


/** Estimates the common multiple `k` by which every peak's statistical uncertainty should be
 scaled for the peaks' scatter about the fitted model to look statistically reasonable.

 Uses an outlier-insensitive statistic (median of the squared pulls, matched to a chi-squared(1)
 median and corrected by num_peaks/dof), so a few badly mis-fit peaks cannot inflate it: `k = sqrt(max(1, median(pull^2)*num_peaks/(dof*0.4549)))`.

 Scaling every peak by one factor does not move the weighted-least-squares solution, so this
 widens the uncertainties without disturbing the fit.

 Returns 1.0 when the peaks are already consistent with their uncertainties, and -1.0 if no
 estimate could be made (unweighted fit, too few peaks, or an unsuccessful solution).
 */
double estimate_stat_uncert_multiple( const RelEffSolution &solution, const double max_multiple = 25.0 );


/** Profile-likelihood interval for one nuclide's mass fraction (of its element).

 For a grid of trial fractions f around the nominal, the problem is re-solved with the target
 nuclide's mass fraction FIXED at f (a lower == upper #MassFractionConstraint), and the
 fit-weight chi2 (see RelEffSolution::m_chi2_fit_weights - the objective the solver actually
 minimizes, so all nuisance parameters are re-fit at each trial value) is recorded; interval
 end-points are the chi2(f) = chi2_min + delta_chi2 crossings, found by bisection.  Compared to
 the covariance-based RelEffSolution::mass_fraction(nuclide, num_sigma), this gives honest
 asymmetric intervals, stays valid near the physical [0,1] bounds and constraint-window edges,
 and does not rely on inverting a possibly ill-conditioned covariance.

 Sub-solves use RelEffInput::point_estimate_only and RelEffInput::skip_an_scan (a physical-model
 self-atten AN is re-fit locally, seeded from the nominal solution), so the whole scan is
 typically well under a second for empirical forms, and a few seconds for the physical model.

 If the target nuclide already has a range mass-fraction constraint, the scan is restricted to
 its window (end-points then flagged `*_at_bound` when clipped).

 Throws std::exception on invalid input (unknown nuclide, element with only one isotope in the
 problem, non-Success nominal solution).

 @param input The same input the nominal solution was solved with.
 @param nominal_solution The nominal solution (must be Success).
 @param options The nuclide to profile, confidence levels, and scan controls.
 */
ProfileMassFractionResult profile_mass_fraction( const RelEffInput &input,
                                                 const RelEffSolution &nominal_solution,
                                                 const ProfileMassFractionOptions &options );


/** Functions in this namespace are for importing peak data from CSV files, and then matching
 up nuclide info, if it wasnt in the CSV files.
 Accepts CSV files from InterSpec, PeakEasy, and GADRAS-DRF.
 */
namespace PeakCsvInput
{

/** Nominally InterSpec uses SandiaDecay for gamma yields everywhere, but for benchmarking
 their are some slightly different Uranium yields available.
 */
enum class NucDataSrc : int
{
  Icrp107_U,
  Lanl_U,
  IcrpLanlGadras_U,
  SandiaDecay,
  Undefined
};//enum class NucDataSrc

const char *to_str( const NucDataSrc src );


struct NuclideInfo
{
  std::string parent, source_nuclide;
  float energy, yield;
  bool optional;
  
  NuclideInfo( const char *p, const char *nuc, bool opt, float kev, float br );
};//struct NuclideInfo



/** Struct to hold results of matching peaks from CSV to nuclides.
 
 For logging purposes peaks or requested nuclides that werent used are also stored.
 */
struct NucMatchResults
{
  /** Peaks (originally from CSV file - e.g., InterSpec, PeakEasy or GADRAS) that were successfully
   matched to at least one source data gamma line.
   */
  std::vector<RelActCalcManual::GenericPeakInfo> peaks_matched;
  
  /** Peaks that were not matched to any source gammas */
  std::vector<RelActCalcManual::GenericPeakInfo> peaks_not_matched;
  
  /** Peaks specifically excluded from the analysis via the 'exclude-peak' command line argument. */
  std::vector<RelActCalcManual::GenericPeakInfo> peaks_excluded;
  
  /** The isotopes (with names normalized to the form like U238, Co60, Pa234m, etc) used for the
   matching.
   */
  std::vector<std::string> used_isotopes;
  
  /** Isotopes that were requested to be used, but no peaks matched to them. */
  std::vector<std::string> unused_isotopes;
  
  /** The ages of the isotopes used; the entries in this vector will correspond to the
   entries in #NucMatchResults::used_isotopes, on an index-by-index basis.
   If age was not applicable (i.e., a Uranium isotope with a non SandiaDecay NucDataSrc),
   then the age will be negative.
   */
  std::vector<double> used_isotope_ages;
  
  /** Where the nuclear source data (e.g., gamma energy and branch ratios) came from. */
  NucDataSrc data_source;
  
  /** The energy tolerance used to match gamma lines to fit peaks; the number of the peaks Gaussian
   sigmas the gamma line must be within the peak mean.
   */
  float match_sigma_tolerance;
  
  /** The energy ranges used; if empty, all energies used. */
  std::vector<std::pair<float,float>> energy_ranges;
  
  /** Source gammas that are matched to at least one peak, and within the allowed energy range, and
   not specifically excluded.
   */
  std::vector<NuclideInfo> source_gammas_used;
  
  /** Source gammas that are in a valid energy range, not specifically excluded, but didnt match to
   any of the peaks fit in data.
   */
  std::vector<NuclideInfo> source_gammas_not_used;
  
  /** The non-decay-corrected values for peaks that received decay-during-measurement corrections.
   Will only have entries in `RelActCalcManual::GenericPeakInfo::m_source_gammas`
   that received corrections.
   */
  std::vector<RelActCalcManual::GenericPeakInfo> not_decay_corrected_peaks;
};//struct NucMatchResults


/** Reads in a peak CSV from InterSpec, PeakEasy, or GADRAS. */
std::vector<RelActCalcManual::GenericPeakInfo> peak_csv_to_peaks( std::istream &csv );

/** Simple struct to hold a nuclides name and its age. */
struct NucAndAge
{
  /** The nuclide name, in a format that SandiaDecay will understand. */
  std::string nuclide;
  
  /** Age of nuclide; a negative value will cause the default age for the nuclide to be used. */
  double age = -1.0;
  
  /** If the nuclides decay during the measurement should be accounted for. */
  bool decay_during_measurement;
  
  NucAndAge( const std::string &nuc, const double the_age, const bool decay_correct_during_meas )
  : nuclide(nuc),
    age( the_age ),
    decay_during_measurement( decay_correct_during_meas )
  {}
};//NucAndAge

/** Matches peaks up to source nuclides, and filters peaks based on energy ranges and, not matching
 an input nuclide gamma, and explicitly not-wanted peaks.
 
 @param peaks Input peaks from the CSV file.
 @param nuc_data_src The nuclear data source to use.  If any non-uranium nuclides are specified,
        this must be #NucDataSrc::SandiaDecay.  When #NucDataSrc::SandiaDecay is used, then default
        ages are assumed.
 @param energy_ranges The energy ranges of peaks to keep.  If empty, will not filter on this.
 @param isotopes The names of isotopes to potentially match up.
 @param energy_tolerance_sigma The matching tolerance, in the peaks sigma; all source gammas within
        this energy range of the peak mean will be attributed to the peak.
 @param excluded_peak_energies Peaks to explicitly exclude from the analysis.
 @param measurement_duration The duration of the measurement; only used if correcting activities for decay during measurement.
        \sa `NucAndAge::decay_during_measurement`
 */
NucMatchResults fill_in_nuclide_info( const std::vector<RelActCalcManual::GenericPeakInfo> peaks,
                                     const NucDataSrc nuc_data_src,
                                     const std::vector<std::pair<float,float>> energy_ranges,
                                     std::vector<NucAndAge> isotopes,
                                     const float energy_tolerance_sigma,
                                     const std::vector<float> excluded_peak_energies,
                                     const float measurement_duration );
}//namespace PeakCsvInput

}//namespace RelActCalcManual

#endif //RelActCalcManual_h
