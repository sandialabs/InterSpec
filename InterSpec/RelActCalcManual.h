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

#include <vector>
#include <memory>
#include <istream>

// Forward declarations
class PeakDef;

namespace RelActCalc
{
  enum class RelEffEqnForm : int;
}

namespace SandiaDecay
{
  struct Nuclide;
}

namespace SpecUtils
{
  class Measurement;
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
};//struct SandiaDecayNucRelAct


/** Information about a SandiaDecay defined nuclide for input into relative efficiency calculation.
 */
struct SandiaDecayNucRelAct
{
  const SandiaDecay::Nuclide *nuclide = nullptr;
  double age = -1.0;
  double rel_activity = -1.0;
};//struct SandiaDecayNucRelAct


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
   */
  double m_energy;
  
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
   
   TODO: Try replacing this value by using a ceres::LossFunction.
   */
  double m_base_rel_eff_uncert;
  
  /** The list of gammas that contributes to this peak.
   
   Must be at least one entry, and duplicate entries are not allowed.
   */
  std::vector<GenericLineInfo> m_source_gammas;
  
  
  /** Simple constructor that initializes all the member variables to zero. */
  GenericPeakInfo();
};//struct GenericPeakInfo



std::vector<GenericPeakInfo> add_nuclides_to_peaks( const std::vector<GenericPeakInfo> &peaks,
                                                   const std::vector<SandiaDecayNuc> &nuclides,
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
void fit_rel_eff_eqn_lls( const RelActCalc::RelEffEqnForm fcn_form,
                           const size_t order,
                           const std::vector<SandiaDecayNucRelAct> &nuclides,
                           const double base_rel_eff_uncert,
                           const std::vector<std::shared_ptr<const PeakDef>> &peak_infos,
                           std::vector<double> &fit_pars,
                           std::vector<std::vector<double>> *covariance );


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
void fit_act_to_rel_eff( const RelActCalc::RelEffEqnForm eqn_form,
                        const std::vector<double> &eqn_pars,
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
   
   TODO: check correlations between isotopes to see if we need to take this into account to have a reasonable uncertainty
   */
  double m_rel_activity_uncert;
};//struct IsotopeRelativeActivity



/** The status of fitting for a solution. */
enum class ManualSolutionStatus : int
{
  NotInitialized,
  ErrorInitializing,
  ErrorFindingSolution,
  ErrorGettingSolution,
  Success
};//enum class ManSolutionStatus

/** A struct to hold the information about the solution to fitting the relative activities and
 efficiency curves.
 */
struct RelEffSolution
{
  RelActCalc::RelEffEqnForm m_rel_eff_eqn_form;
  size_t m_rel_eff_eqn_order = 0;
  std::vector<double> m_rel_eff_eqn_coefficients;
  std::vector<std::vector<double>> m_rel_eff_eqn_covariance;
  
  /** The relative activities for each of the input nuclides. */
  std::vector<IsotopeRelativeActivity> m_rel_activities;

  /** Covariance matrix of relative activities; same index ordering as \p m_rel_activities. */
  std::vector<std::vector<double>> m_rel_act_covariance;

  /** The Chi2 summed over all peaks between their actual and fit relative efficiencies.
   
   Note that this always uses the peaks statistical uncertainties, and does not include
   #GenericPeakInfo::m_base_rel_eff_uncert
   */
  double m_chi2 = 0.0;
  
  /** The number of degrees of freedom in the fit for equation parameters.
   
   Note: right now just have as just the number of peaks used minus number of equation parameters fit for, minus one less than the
   number of isotopes - probably off by one, or more - need to think on this and come back to it
   */
  size_t m_dof = 0;
  
  /** The number of evaluation calls it took L-M to reach a solution.
   Only useful for debugging and curiosity.
   */
  int m_num_function_eval_solution = 0;
  
  /** The number of evaluation calls it took to reach a solution, and compute final covariance. */
  int m_num_function_eval_total = 0;
  
  /** How long it took to compute the answer (only for curiosity) */
  int m_num_microseconds_eval = 0;
  
  /** As an internal detail of fitting the relative efficiencies, we normalize the activities to a flat relative efficiency curve of 1.0, and then
   fit for the multiple of the normalization that yields the best solution.  This member variable keeps track of these normalizations; we
   are keeping them around to help in computing the correlation compensated ratios (although we could instead modify
   m_rel_act_covariance - but we'll just keep an extra variable around to make things a little easier to debug).
   
   The entries in this vector correspond to \p m_rel_activities on an index-by-index basis.
   
   To state it another way, these are the relative activities if the relative efficiency curve is 1.0.
   */
  std::vector<double> m_activity_norms;
  
  /** The status of creating the solution. */
  ManualSolutionStatus m_status = ManualSolutionStatus::NotInitialized;
  
  /** Error message if computation wasn't successful. */
  std::string m_error_message;
  
  /** Warnings about the setup or solution of the problem; by no means comprehensive of potential
   issues!
   */
  std::vector<std::string> m_warnings;
  
  /** The input peaks that were used to create the solution. */
  std::vector<GenericPeakInfo> m_input_peak;
  
  
  /////////////////////////////////////////////////////////////////////////////////////////////////
  //            Below here are member functions for simplified access to information             //
  /////////////////////////////////////////////////////////////////////////////////////////////////
  
  /** Returns the index in #m_rel_activities for the specified nuclide.
   
   Throws std::exception if an invalid nuclide.
   */
  size_t nuclide_index( const std::string &nuc ) const;
  
  /** The relative activity of a nuclide.
   
   Note that this assumes the peak counts are in CPS (or equivalently the data is 1-second long),
   so you will need to divide by live-time to compare to the "auto" relative activity.
   
   Throws std::exception if an invalid nuclide.
   */
  double relative_activity( const std::string &nuclide ) const;
  
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
  
  /** Returns the mass fraction of the specified nuclide. */
  double mass_fraction( const std::string &iso ) const;
  
  /** Prints out a summary of the results to the provided stream; for development/debug. */
  std::ostream &print_summary( std::ostream &strm ) const;
  
  /** Creates a self-contained HTML report of the results.
   
   @param strm The stream to place the HTML file into.
   @param spectrum_title The title to display on the HTML pace
   @param spectrum The optional spectrum to display on the HTML page (may be nullptr)
   @param spectrum_display_peaks The peaks to display on the spectrum; may be empty.
   */
  void print_html_report( std::ostream &strm,
                         std::string spectrum_title,
                         std::shared_ptr<const SpecUtils::Measurement> spectrum,
                         std::vector<std::shared_ptr<const PeakDef>> spectrum_display_peaks
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
};//struct RelEffSolution

/** Solve for the relative efficiency equation and relative activities for all isotopes
 
 @param peak_infos The input peak information.
 Each peak must have at least one source isotope, have counts above zero, have counts
 uncertainty greater than zero (or equal to -1.0), and each isotope must have yields
 at each energy greater than 0.
 @param eqn_form The equation form to use; i.e., EqnForm::LnX, EqnForm::LnY, or EqnForm::LnXLn
 @param eqn_order The order of equation to be fit for.  i.e., the number of energy-dependent terms
 to be fit for (e.g., the total number of terms will be this plus one).

 @returns The solution to the problem.  Make sure to check RelEffSolution::m_status,
          RelEffSolution::m_warnings, and RelEffSolution::m_error_message.
 */
RelEffSolution solve_relative_efficiency( const std::vector<GenericPeakInfo> &peak_infos,
                                         const RelActCalc::RelEffEqnForm eqn_form,
                                         const size_t eqn_order );


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
  
  NucAndAge( const std::string &nuc, const double the_age )
  : nuclide(nuc),
    age( the_age )
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
 @param energy_tolerance_fwhm The matching tolerance, in the peaks FWHM; all source gammas within
        this energy range of the peak mean will be attributed to the peak.  A value of 0.57 will
        include all gammas within 1.5 sigma of the mean.
 @param excluded_peak_energies Peaks to explicitly exclude from the analysis.
 */
NucMatchResults fill_in_nuclide_info( const std::vector<RelActCalcManual::GenericPeakInfo> peaks,
                                     const NucDataSrc nuc_data_src,
                                     const std::vector<std::pair<float,float>> energy_ranges,
                                     std::vector<NucAndAge> isotopes,
                                     const float energy_tolerance_sigma,
                                     const std::vector<float> excluded_peak_energies );
}//namespace PeakCsvInput

}//namespace RelActCalcManual

#endif //RelActCalcManual_h
