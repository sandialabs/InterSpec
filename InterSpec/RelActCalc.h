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

}//namespace RelActCalc


#endif //RelActCalc_h
