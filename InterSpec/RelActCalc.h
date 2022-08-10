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

#include <vector>
#include <cstddef>

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
};//enum class RelEffEqnForm

/** Returns string representation of RelEffEqnForm, i.e., "LnX", "LnY", "LnXLnY", or "FRAM Empirical" */
const char *to_str( const RelEffEqnForm form );

RelEffEqnForm rel_eff_eqn_form_from_str( const char *str );

/** Returns a human-readable representation of the relative efficiency equation. */
std::string rel_eff_eqn_text( const RelEffEqnForm eqn_form, const std::vector<double> &coefs );

/** Returns a javascript expression of the relative efficiency equation.
 
 Example return value: "function(x){ return Math.exp( 1.2 + x + ... ); }"
 */
std::string rel_eff_eqn_js_function( const RelEffEqnForm eqn_form, const std::vector<double> &coefs );

/** Evaluate the relative efficiency equation at a given energy, for a specific form of the equation, and coefficients. */
double eval_eqn( const double energy, const RelEffEqnForm eqn_form,
                const double * const coefs, const size_t num_coefs );

/** A convenience call signature for the above #eval_eqn */
double eval_eqn( const double energy, const RelEffEqnForm eqn_form, const std::vector<double> &coefs );

/** Evaluate the uncertainty in the relative efficiency equation - assuming all uncertainties are uncorrelated - this will tend to
 over-estimate the errors.
 
 TODO: this function can return NaN values - when the square uncertainty is negative...
 */
double eval_eqn_uncertainty( const double energy, const RelEffEqnForm eqn_form,
                             const std::vector<std::vector<double>> &covariance );



/** A struct to specify the relative masses of the Pu and Am241 nuclides that
 are observable by gamma spec.
 
 Masses do not need to be normalized to 1, but their ratios to each other do
 need to be correct.
 */
struct Pu242ByCorrelationInput
{
  float pu238_rel_mass = 0.0f;
  float pu239_rel_mass = 0.0f;
  float pu240_rel_mass = 0.0f;
  float pu241_rel_mass = 0.0f;
  float am241_rel_mass = 0.0f;
  /** A value to capture all other non-Pu242 plutonium isotopes (I dont expect to ever be used) */
  float other_pu_mass  = 0.0f;
};//struct Pu242ByCorrelationInput


/** The nuclide mass fractions (so between 0 and 1 - with all nuclides summing to 1.0)
 as determined after accounting for the predicted amount (by correlation) of Pu242.
 */
struct Pu242ByCorrelationOutput
{
  float pu238_mass_frac;
  float pu239_mass_frac;
  float pu240_mass_frac;
  float pu241_mass_frac;
  float am241_mass_frac;
  float pu242_mass_frac;
  
  /** If enrichment is within literatures specified range.
   This is a very coarse indicator if its a valid method.
   */
  bool is_within_range;
  
  /** The very coarsely estimated fractional uncertainty on the Pu242 component
   from the correlation estimate; roughly taken from the literature, or some
   wild guess, but I guess better than nothing.  Does not cover additional
   uncerts from using the method outside the literatures specified conditions.
   */
  float pu242_uncert;
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
   In FRAM, this function is the equivalent of setting ‘‘Pu242-correlation’’ value to  9.66E-03 and
   the ‘‘Pu239_exponent’’ to -3.83 and the other coefficients to zero.
   
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
 
 Currently doesnt take into account decay of Am241 rel to Pu241, beyond just adding
 Am241 mass to Pu241.
 */
Pu242ByCorrelationOutput correct_pu_mass_fractions_for_pu242( Pu242ByCorrelationInput input, PuCorrMethod method );

void test_pu242_by_correlation();
}//namespace RelActCalc


#endif //RelActCalc_h
