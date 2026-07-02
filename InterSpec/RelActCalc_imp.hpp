#ifndef RelActCalc_imp_hpp
#define RelActCalc_imp_hpp
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

#include <tuple>
#include <vector>
#include <cassert>
#include <utility>
#include <optional>
#include <algorithm>
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/GammaInteractionCalc.h"

#ifdef _MSC_VER
#undef isinf
#undef isnan
#undef isfinite
#undef isnormal
#undef ERROR
#endif

namespace ceres
{
  /* dummy namespace for when using this file for only doubles, and ceres.h hasnt been included */
}

/** This .hpp file is is to allow evaluating Relative Efficiency functions as either `double` or `ceres::Jet<>`.
 
 Its a little messy, but you need to include this .hpp file after including ceres.h, and your std headers.
 */
namespace RelActCalc
{
  
/** Fraction of an element's post-fixed-constraint budget (`1 - Σ fixed`) kept in reserve for the
 element's unconstrained nuclides: the range-constrained nuclides' total mass fraction is
 hard-bounded (a Ceres box constraint on the carrier parameter) at
 `1 - fixed_sum - ns_mass_frac_min_remainder_frac*(1 - fixed_sum)`, so the `1/(1 - sum_constrained)`
 factor in the relative-activity decode can never blow up, with no smooth-cap warping of the
 feasible region (the old soft-cap knee compressed everything above 95% of budget - a single [0,1]
 constraint could never decode above ~0.98).  When every nuclide of the element is constrained
 there is no reserve (the fractions must sum to exactly 1).  See #MassFracBlockSpec. */
constexpr double ns_mass_frac_min_remainder_frac = 1.0e-6;

/** Relative half-width of the quadratic-hinge smoothing zones used when distributing the total
 constrained mass fraction among the individual windows (#decode_mass_frac_block): each window's
 conditional feasible interval gets its endpoints smoothed over
 `ns_mass_frac_stick_radius_frac * (upper - lower)` of that window, so the decode is C1 in every
 parameter (safe through `ceres::Jet`) while the windows stay exactly respected (see the one-sided
 bound properties on #qmax_hinge).  Always RELATIVE to the window involved - tiny windows (e.g.
 U232 in [0, 0.9E-9]) are legitimate, and must never be swamped by an absolute epsilon. */
constexpr double ns_mass_frac_stick_radius_frac = 1.0e-3;


/** C1 smooth approximation of `max(x, a)` with compact support: exact for `|x - a| >= r`, a
 quadratic blend inside.  One-sided-bound properties (relied on by #decode_mass_frac_block):
   - `qmax_hinge(x,a,r) >= a` and `>= x` everywhere (never below the true max);
   - the excess over `max(x,a)` is at most `r/4` (largest at `x == a`);
   - value and first derivative are continuous at the seams, so it is safe to evaluate and
     differentiate through `ceres::Jet` (the branch compares only the scalar part, and both value
     and slope match across the branch points).
 */
template<typename T>
T qmax_hinge( const T &x, const double a, const double r )
{
  assert( r > 0.0 );

  const T d = x - a;
  if( d >= r )
    return x;
  if( d <= -r )
    return T(a);
  return a + (d + r)*(d + r) / (4.0*r);
}


/** C1 smooth approximation of `min(x, b)` - the mirror of #qmax_hinge, via
 `min(x,b) = x + b - max(x,b)`; `qmin_hinge(x,b,r) <= min(x,b)` everywhere, deficit at most `r/4`. */
template<typename T>
T qmin_hinge( const T &x, const double b, const double r )
{
  return x + b - qmax_hinge( x, b, r );
}


/** Per-element specification of the exact "sigma-block" decode for mass-fraction constraints.

 An element's mass-fraction-constrained nuclides split into FIXED constraints (lower == upper;
 their parameter stays constant and they contribute only `fixed_sum`), and RANGE constraints
 (lower < upper; `num_range = lower.size()` of them, in a deterministic order with the "carrier"
 first).  Each range-constrained nuclide owns exactly one activity slot in the fit:

   - The carrier slot holds `t in [0,1]`, mapped to the range-constrained nuclides' TOTAL mass
     fraction `sigma = sig_lo + t*(sig_hi - sig_lo)`.  `sig_hi <= 1 - fixed_sum - delta` is a hard
     Ceres box bound, so the element's constrained sum `fixed_sum + sigma` is structurally below 1
     with EXACT gradients - no throw, no soft-cap warp; an infeasible demand surfaces as the
     parameter pinning at its bound (which the existing at-bound warning machinery reports).
     When `all_constrained` (the element has no unconstrained nuclide), `sigma` is instead the
     constant `1 - fixed_sum`, and the carrier slot holds the element's total (relative mass)
     scale.
   - The remaining `num_range - 1` slots hold `g_k in [0,1]`, distributing `sigma` among the
     individual windows through sequential conditional intervals (#decode_mass_frac_block); the
     carrier nuclide takes the exact remainder, in-window by construction.
 */
struct MassFracBlockSpec
{
  /** Σ of the element's fixed (lower == upper) constrained mass fractions. */
  double fixed_sum = 0.0;

  /** True when every nuclide of the element carries a mass-fraction constraint. */
  bool all_constrained = false;

  /** Hard margin below 1 reserved for the element's unconstrained nuclides
   (`ns_mass_frac_min_remainder_frac * (1 - fixed_sum)`); 0 when `all_constrained`. */
  double delta = 0.0;

  /** Box for the range-constrained total: `sig_lo = Σ lower`,
   `sig_hi = min( Σ upper, 1 - fixed_sum - delta )`; both equal `1 - fixed_sum` when
   `all_constrained`. */
  double sig_lo = 0.0, sig_hi = 0.0;

  /** Windows of the range-constrained nuclides, carrier first. */
  std::vector<double> lower, upper;

  /** Suffix window sums over the yet-unassigned set after assigning window k (the carrier plus
   windows k+1..): `after_lower[k] = lower[0] + Σ_{j>k} lower[j]`, likewise `after_upper[k]`;
   entries 1..num_range-1 are what the decode uses. */
  std::vector<double> after_lower, after_upper;

  /** Hinge half-width for window k: `ns_mass_frac_stick_radius_frac * (upper[k] - lower[k])`. */
  std::vector<double> radius;

  /** Hinge half-width of the pinched-conditional-interval (width) guard - relative to the
   smallest window in the block, capped at 1E-12. */
  double width_radius = 1.0e-12;
};//struct MassFracBlockSpec


/** Builds the #MassFracBlockSpec for one element.

 @param range_windows   (lower,upper) of the element's range-constrained (lower < upper) nuclides,
                        in the deterministic block order - the first entry is the carrier.
 @param fixed_sum       Σ of the element's fixed (lower == upper) constrained fractions.
 @param all_constrained Whether every nuclide of the element is mass-fraction constrained.

 Feasibility preconditions are enforced by `check_nuclide_constraints()`:
 mixed: `fixed_sum + Σ lower < 1`; all-constrained: `fixed_sum + Σ lower <= 1 <= fixed_sum + Σ upper`.
 */
inline MassFracBlockSpec make_mass_frac_block_spec( const std::vector<std::pair<double,double>> &range_windows,
                                                    const double fixed_sum,
                                                    const bool all_constrained )
{
  MassFracBlockSpec spec;
  spec.fixed_sum = fixed_sum;
  spec.all_constrained = all_constrained;

  const size_t num_range = range_windows.size();
  spec.lower.resize( num_range );
  spec.upper.resize( num_range );

  double lower_sum = 0.0, upper_sum = 0.0, min_width = 1.0;
  for( size_t k = 0; k < num_range; ++k )
  {
    spec.lower[k] = range_windows[k].first;
    spec.upper[k] = range_windows[k].second;
    assert( spec.lower[k] < spec.upper[k] ); //fixed (lower == upper) constraints belong in `fixed_sum`
    lower_sum += spec.lower[k];
    upper_sum += spec.upper[k];
    min_width = std::min( min_width, spec.upper[k] - spec.lower[k] );
  }

  // For a mixed element check_nuclide_constraints() guarantees `fixed_sum < 1`; an all-constrained
  //  all-fixed element has `fixed_sum == 1` (to ~1E-6), so clamp the tiny/negative leftover to 0.
  const double budget = std::max( 1.0 - fixed_sum, 0.0 );
  assert( all_constrained || (budget > 0.0) );

  if( all_constrained )
  {
    // No unconstrained nuclides: the range fractions must sum to exactly the leftover budget.
    spec.delta = 0.0;
    spec.sig_lo = spec.sig_hi = budget;
  }else
  {
    spec.delta = ns_mass_frac_min_remainder_frac * budget;
    spec.sig_lo = lower_sum;
    spec.sig_hi = std::min( upper_sum, budget - spec.delta );
    // check_nuclide_constraints() guarantees `lower_sum < budget`; keep the box non-inverted in
    //  edge cases anyway (a point box gets its parameter pinned const by the setup code).
    spec.sig_hi = std::max( spec.sig_hi, spec.sig_lo );
  }

  if( num_range > 0 )
  {
    spec.after_lower.resize( num_range, 0.0 );
    spec.after_upper.resize( num_range, 0.0 );
    double l_tail = spec.lower[0], u_tail = spec.upper[0];
    for( size_t k = num_range - 1; k > 0; --k )
    {
      spec.after_lower[k] = l_tail;
      spec.after_upper[k] = u_tail;
      l_tail += spec.lower[k];
      u_tail += spec.upper[k];
    }
    spec.after_lower[0] = l_tail; //== Σ all lowers; unused by the decode, but handy for asserts
    spec.after_upper[0] = u_tail; //== Σ all uppers

    spec.radius.resize( num_range, 0.0 );
    for( size_t k = 0; k < num_range; ++k )
      spec.radius[k] = ns_mass_frac_stick_radius_frac * (spec.upper[k] - spec.lower[k]);

    spec.width_radius = std::min( 1.0e-12, ns_mass_frac_stick_radius_frac * min_width );
  }//if( num_range > 0 )

  return spec;
}//make_mass_frac_block_spec(...)


/** Decodes the element's range-constrained mass fractions from the block parameters.

 @param spec      The block specification (#make_mass_frac_block_spec).
 @param sigma     Total mass fraction of the range-constrained nuclides, in
                  `[spec.sig_lo, spec.sig_hi]` (from the carrier parameter, or the constant
                  `1 - fixed_sum` when `spec.all_constrained`).
 @param gs        The `num_range - 1` distribution values, each in [0,1] (may be nullptr when
                  `num_range <= 1`).
 @param fractions Receives the `num_range` decoded fractions, in block order (carrier first).

 Guarantees: `Σ fractions == sigma` exactly; `fractions[k] >= lower[k]` exactly;
 `fractions[k] <= upper[k] + width_radius/4` (exact to ~2.5E-13); C1 in (sigma, gs), so safe
 through `ceres::Jet`.
 */
template<typename T>
void decode_mass_frac_block( const MassFracBlockSpec &spec, const T &sigma, const T * const gs,
                             T * const fractions )
{
  const size_t num_range = spec.lower.size();
  if( num_range == 0 )
    return;

  T rem = sigma;
  for( size_t k = 1; k < num_range; ++k )
  {
    // The exact conditional feasible interval for f_k, given the remainder and that the
    //  yet-unassigned windows must still be satisfiable, is
    //  [max(lower_k, rem - after_upper_k), min(upper_k, rem - after_lower_k)] - non-empty for any
    //  sigma inside its box; the hinges round its corners so the decode is C1 everywhere.
    const T lo = qmax_hinge( rem - spec.after_upper[k], spec.lower[k], spec.radius[k] );
    const T hi = qmin_hinge( rem - spec.after_lower[k], spec.upper[k], spec.radius[k] );
    const T width = qmax_hinge( hi - lo, 0.0, spec.width_radius );

    fractions[k] = lo + gs[k-1]*width;
    rem -= fractions[k];
  }//for( size_t k = 1; k < num_range; ++k )

  fractions[0] = rem; //the carrier absorbs the exact remainder (in-window by construction)
}//decode_mass_frac_block(...)


/** Inverts #decode_mass_frac_block for starting values: given target fractions (e.g., from the
 manual-solution estimate), computes a `sigma` and `gs` that decode to (approximately) those
 fractions, each kept strictly inside its box by `margin` (Ceres' projected-gradient steps behave
 poorly for a parameter that starts exactly on a bound).

 Runs the same smoothed forward decode incrementally, so the start is exactly self-consistent;
 the clamps absorb any infeasibility of the targets.
 */
inline void invert_mass_frac_block( const MassFracBlockSpec &spec, const double * const target_fractions,
                                    double &sigma, double * const gs, const double margin = 1.0e-3 )
{
  const size_t num_range = spec.lower.size();
  if( num_range == 0 )
  {
    sigma = spec.sig_lo;
    return;
  }

  double target_sigma = 0.0;
  for( size_t k = 0; k < num_range; ++k )
    target_sigma += target_fractions[k];

  if( spec.all_constrained || (spec.sig_hi <= spec.sig_lo) )
  {
    sigma = spec.sig_hi; //the all-constrained constant, or a pinched box
  }else
  {
    const double t_raw = (target_sigma - spec.sig_lo) / (spec.sig_hi - spec.sig_lo);
    const double t = std::min( 1.0 - margin, std::max( margin, t_raw ) );
    sigma = spec.sig_lo + t*(spec.sig_hi - spec.sig_lo);
  }

  double rem = sigma;
  for( size_t k = 1; k < num_range; ++k )
  {
    const double lo = qmax_hinge( rem - spec.after_upper[k], spec.lower[k], spec.radius[k] );
    const double hi = qmin_hinge( rem - spec.after_lower[k], spec.upper[k], spec.radius[k] );
    const double width = qmax_hinge( hi - lo, 0.0, spec.width_radius );

    // Solve for g only when the conditional interval has meaningful width RELATIVE to this
    //  window (never vs an absolute epsilon - tiny windows are legitimate); else the midpoint.
    double g = 0.5;
    if( width > 1.0e-9*(spec.upper[k] - spec.lower[k]) )
      g = std::min( 1.0 - margin, std::max( margin, (target_fractions[k] - lo)/width ) );

    gs[k-1] = g;
    rem -= (lo + g*width);
  }//for( size_t k = 1; k < num_range; ++k )
}//invert_mass_frac_block(...)


template<typename T>
T eval_eqn_imp( const double energy, const RelActCalc::RelEffEqnForm eqn_form,
                  const T * const coeffs, const size_t num_coefs )
{
  using namespace std;
  using namespace ceres;
  
  if( energy <= 0.0 )
    throw std::runtime_error( "eval_eqn: energy must be greater than zero." );
  
  if( num_coefs < 1 )
    throw std::runtime_error( "eval_eqn: need at least one coefficients." );
  
  
  T answer( 0.0 );
  
  switch( eqn_form )
  {
    case RelActCalc::RelEffEqnForm::LnX:
    {
      // y = a + b*ln(x) + c*(ln(x))^2 + d*(ln(x))^3 + ...
      const double log_energy = std::log(energy);
      
      for( size_t order = 0; order < num_coefs; ++order )
      {
        switch( order )
        {
          case 0:  answer += coeffs[order];                                    break;
          case 1:  answer += coeffs[order] * log_energy;                       break;
          default: answer += coeffs[order] * pow( log_energy, (double)order ); break;
        }//switch( order )
      }//for( loop over coeffs )
      
      break;
    }//case RelEffEqnForm::LnX:
      
    case RelActCalc::RelEffEqnForm::LnY:
    {
      // y = exp( a + b*x + c/x + d/x^2 + e/x^3 + ... )
      // Note that it would be a little more convenient to have y = exp(a*x + b + c/x + ...), but
      //  I want to keep the leading term independent of energy.
      for( size_t order = 0; order < num_coefs; ++order )
      {
        switch( order )
        {
          case 0:  answer += coeffs[order];                              break;
          case 1:  answer += coeffs[order] * energy;                     break;
          case 2:  answer += coeffs[order] / energy;                     break;
          default: answer += coeffs[order] / pow( energy, order - 1.0 ); break;
        }//switch( order )
      }//for( loop over coeffs )
      
      answer = exp( answer );
      
      break;
    }//case RelEffEqnForm::LnY:
      
    case RelActCalc::RelEffEqnForm::LnXLnY:
    {
      // y = exp (a  + b*(lnx) + c*(lnx)^2 + d*(lnx)^3 + ... )
      const double log_energy = std::log(energy);
      
      for( size_t order = 0; order < num_coefs; ++order )
      {
        switch( order )
        {
          case 0:  answer += coeffs[order]; break;
          case 1:  answer += coeffs[order] * log_energy; break;
          default: answer += coeffs[order] * pow( log_energy, (double)order ); break;
        }//switch( order )
      }//for( loop over coeffs )
      
      answer = exp( answer );
      
      break;
    }//case RelEffEqnForm::LnXLnY:
      
    case RelActCalc::RelEffEqnForm::FramEmpirical:
    {
      // y = exp( a + b/x^2 + c*(lnx) + d*(lnx)^2 + e*(lnx)^3 )
      const double log_energy = log(energy);
      
      for( size_t order = 0; order < num_coefs; ++order )
      {
        switch( order )
        {
          case 0:  answer += coeffs[order];  break;
          case 1:  answer += coeffs[order] / (energy*energy); break;
          default: answer += coeffs[order] * pow( log_energy, order - 1.0 ); break;
        }//switch( order )
      }//for( loop over coeffs )
      
      answer = exp( answer );
      
      break;
    }//case RelEffEqnForm::FramEmpirical:
      
    case RelActCalc::RelEffEqnForm::FramPhysicalModel:
      throw std::runtime_error( "RelActCalc::eval_eqn() can not be called for FramPhysicalModel." );
  }//switch( eqn_form )
  
  return answer;
}//eval_eqn(...)
    

/** This function is the equivalent of `mass_attenuation_coef(...)`, but preserves ceres::Jet
    derivative information by linearly interpolating between floor(AN) and ceil(AN).
*/
template<typename T>
T get_atten_coef_for_an( const T &an, const float energy )
{
  // `T` is either a `ceres::Jet<>` or a `double` here, so we'll use compile time if
  //  to get the scalar part
  //  (the better way to do this is to use `auto an_scalar = ceres::internal::AsScalar(an);`, but
  //   that requires including "ceres/internal/jet_traits.h", so we'll skip it)
  double an_scalar;
  if constexpr ( !std::is_same_v<T, double> )
    an_scalar = an.a;
  else
    an_scalar = an;

  // Allow a small margin past [1,98]: callers no longer hard-clamp the (Jet) atomic number to the
  //  bound (that would zero its derivative), so the optimizer may evaluate a hair outside while the
  //  table lookup below stays valid via std::clamp.
  assert( (an_scalar >= 1.0 - 1.0e-3) && (an_scalar <= 98.0 + 1.0e-3) );

  // Clamp to [1, 97] so we can always interpolate to lower_an+1 without exceeding element 98
  const int lower_an = std::clamp( static_cast<int>(std::floor(an_scalar)), 1, 97 );
  const int upper_an = lower_an + 1;

  const double lower_mu = MassAttenuation::massAttenuationCoefficientElement(lower_an, energy);
  const double upper_mu = MassAttenuation::massAttenuationCoefficientElement(upper_an, energy);
  const T anfrac = an - static_cast<double>(lower_an);  //This preserves the Jet derivative
  const T mu = (1.0 - anfrac)*lower_mu + anfrac*upper_mu;

  return mu;
}//T get_atten_coef_for_an( const T &an )


/** Evaluates the pivot-anchored Chebyshev empirical correction in log-energy:

   correction(E) = exp( sum_{k=1..N} coeffs[k-1] * ( T_k(u(E)) - T_k(u(E0)) ) )

 with `u(E) = 2*(log E - log Elo)/(log Ehi - log Elo) - 1` and `T_k` the Chebyshev polynomial of degree k.
 The basis polynomials depend only on energy (not the fit coefficients), so they are evaluated in `double`;
 only the linear combination with `coeffs` (which may be a `ceres::Jet`) carries derivatives.  The pivot
 anchor (subtracting `T_k(u(E0))`) forces `correction(E0) = 1`, so the correction contributes only SHAPE, not
 overall scale - removing the scale degeneracy with the relative activities. */
template<typename T>
T physical_model_basis_correction( const double energy,
                                   const std::vector<T> &coeffs,
                                   const double lower_energy,
                                   const double upper_energy,
                                   const double pivot_energy )
{
  if( coeffs.empty() )
    return T( 1.0 );

  // If the energy reference frame was not supplied (e.g. a caller that does not yet plumb it through),
  //  skip the correction (identity) rather than producing garbage.
  if( !((lower_energy > 0.0) && (upper_energy > lower_energy) && (pivot_energy > 0.0)) )
    return T( 1.0 );

  const double log_lo = std::log( lower_energy );
  const double log_span = std::log( upper_energy ) - log_lo;
  assert( log_span > 0.0 );

  // Map energy -> [-1,1] in log-energy (efficiency/attenuation are ~power laws -> ~linear in log E).
  const double u  = 2.0*(std::log( energy )       - log_lo)/log_span - 1.0;
  const double u0 = 2.0*(std::log( pivot_energy ) - log_lo)/log_span - 1.0;

  // Chebyshev recurrence (in double; basis is independent of the fit coefficients):
  //   T_0=1, T_1=x, T_{k+1} = 2x*T_k - T_{k-1}
  double Tkm1_u = 1.0, Tk_u = u;    // T_0, T_1 evaluated at u
  double Tkm1_0 = 1.0, Tk_0 = u0;   // T_0, T_1 evaluated at u0

  T exponent( 0.0 );
  for( size_t k = 1; k <= coeffs.size(); ++k )
  {
    exponent += coeffs[k-1] * (Tk_u - Tk_0);  // anchored basis term T_k(u) - T_k(u0)

    const double Tnext_u = 2.0*u *Tk_u - Tkm1_u;
    const double Tnext_0 = 2.0*u0*Tk_0 - Tkm1_0;
    Tkm1_u = Tk_u; Tk_u = Tnext_u;
    Tkm1_0 = Tk_0; Tk_0 = Tnext_0;
  }//for( k = 1 .. N )

  using std::exp;
  return exp( exponent );
}//physical_model_basis_correction(...)


template<typename T>
T eval_physical_model_eqn_imp( const double energy,
                                std::optional<RelActCalc::PhysModelShield<T>> self_atten,
                                const std::vector<RelActCalc::PhysModelShield<T>> &external_attens,
                                const DetectorPeakResponse * const drf,
                                std::optional<T> b,
                                std::optional<T> c,
                                const double corr_lower_energy = 0.0,
                                const double corr_upper_energy = 0.0,
                                const double corr_pivot_energy = 0.0,
                                const RelActCalc::PhysModelCorrFcn corr_fcn = RelActCalc::PhysModelCorrFcn::Hoerl )
{
  // Note 20250114: this function currently interpolates between floor(AN) and one above it.
  //                However, I dont think this is correct, we should be interpolating a bit better
  //                so at least the derivative comes out better (think of an integer AN, we should
  //                take into account the AN one below, as well as the one above)
  
  using namespace std;
  using namespace ceres; //So we can use the math functions defined in the `ceres` namespace for Jets
  
  const float energyf = static_cast<float>( energy );
  
  assert( b.has_value() == c.has_value() );
  
  T answer( 1.0 );
  
  assert( b.has_value() == c.has_value() );
  if( b.has_value() != c.has_value() )
    throw std::logic_error( "hoerl_b.has_value() != hoerl_c.has_value()" );
  
  if( (corr_fcn != RelActCalc::PhysModelCorrFcn::None) && b.has_value() && c.has_value() )
  {
    if( corr_fcn == RelActCalc::PhysModelCorrFcn::Chebyshev )
    {
      // Pivot-anchored Chebyshev correction; `b`,`c` carry the basis coefficients a_1, a_2.
      const std::vector<T> corr_coeffs{ *b, *c };
      answer = physical_model_basis_correction( energy, corr_coeffs,
                                  corr_lower_energy, corr_upper_energy, corr_pivot_energy );
    }else
    {
      // Modified Hoerl function E_MeV^b * c^(1/E_MeV); evaluated in exp/log form for AD/Jet stability.
      const T b_val = *b;
      const T c_val = *c;
      // The b/c bounds are sized to keep the exponent's swing bounded over the correction's fit
      //  window [corr_lower_energy, corr_upper_energy] (see RelActCalcAuto.cpp ~4104).  The correction
      //  is meaningless outside that window, and below it the 1/E_MeV term grows without bound, blowing
      //  the exponent up until exp() overflows (-> the throw below).  Clamp the evaluation energy into
      //  the window so the exponent stays in the range the bounds were sized for - preventing the
      //  overflow at its source.  `energy` is a plain double here, so this does not touch any
      //  ceres::Jet derivative.
      double eval_energy = energy;
      if( (corr_upper_energy > corr_lower_energy) && (corr_lower_energy > 0.0) )
        eval_energy = std::max( corr_lower_energy, std::min( corr_upper_energy, energy ) );
      const double energy_mev = 0.001*eval_energy;
      answer = exp( b_val*log(energy_mev) + log(c_val)/energy_mev );
    }//if( basis ) / else( Hoerl )

    assert( !isnan(answer) && !isinf(answer) );
    if( isnan(answer) || isinf(answer) )
      throw std::logic_error( "physical-model correction gives eqn NaN or Inf" );
  }//if( correction active )
  
  const double det_part = drf ? drf->intrinsicEfficiency( energyf ) : 1.0;
  answer *= det_part;
  
  assert( !isnan(answer) && !isinf(answer) );
  
  if( self_atten.has_value() )
  {
    T mu( 0.0 );
    if( self_atten->material )
      mu = T( GammaInteractionCalc::transmition_length_coefficient( self_atten->material.get(), energyf ) / self_atten->material->density );
    else
      mu = get_atten_coef_for_an( self_atten->atomic_number, energyf );
    
    T areal_density = self_atten->areal_density;

    // `setup_physical_model_shield_par` loosens the Ceres-level lower bound on AD by 1e-5 g/cm^2
    //  to escape an active-bound LM trap, so the optimizer may evaluate at AD a hair below zero.
    //  We deliberately do NOT clamp `areal_density` to >= 0 here: the self-attenuation factor
    //  (1 - e^{-x})/x is smooth and accurate for the tiny negative x reachable inside that margin,
    //  so letting the true value/gradient flow keeps the AD Jacobian column non-zero (clamping a
    //  `ceres::Jet` with `fmax(jet,0)` would zero its derivative and recreate the trap - see A2).
    //  The converged AD is clamped to >= 0 where the solution is extracted (`get_shield_info` in
    //  RelActCalcAuto.cpp), so the user never sees a negative AD; here we only guard grossly
    //  negative values. `areal_density` is in PhysicalUnits, so the -1e-3 g/cm^2 throw threshold
    //  (well below the -1e-5 loosening) compares in physical units.
    assert( !isinf(areal_density) );
    if( (areal_density <= -1.0e-3 * PhysicalUnits::g_per_cm2)
        || isnan(areal_density) || isinf(areal_density) )
      throw std::runtime_error( "eval_physical_model_eqn: areal density must be >= 0 - got value " );

    assert( mu >= 0.0 );
    if( mu < 0.0 )
      mu = fmax(mu, 0.0);

    assert( areal_density >= -1.0e-3 * PhysicalUnits::g_per_cm2 );

    //The simple computation to handel self-attenuation would be:
    //  answer *= (1.0 - exp(-mu * areal_density)) / (mu * areal_density);
    //But if -mu*areal_density is zero, we get NaN, and even as it goes towards zero, we will greatly lose precision

    const T epsilon = T(1.0e-8);
    const T x = mu * areal_density;

    if( x < epsilon )
    {
      // Use the Taylor series expansion for small x.
      // f(x) = 1 - x/2 + x^2/6
      answer *= (T(1.0) - 0.5*x + ((x*x) / 6.0)); // + -(x*x*x)/24.0 + ...
    }else
    {
      // For larger x, use the standard formula, but with expm1
      // to avoid catastrophic cancellation near zero.
      // 1 - exp(-x) = -expm1(-x)
      answer *= -expm1(-x) / x;
    }//if( (mu*areal_density) < epsilon ) / else

    assert( !isnan(answer) && !isinf(answer) );
  }//if( self_atten->has_value() )
  
  for( const RelActCalc::PhysModelShield<T> &ext_atten : external_attens )
  {
    // TODO: `GammaInteractionCalc::transmition_length_coefficient` can be a real bottleneck of computation - at soem point we should memoise its results
    T mu( 0.0 );
    if( ext_atten.material )
      mu = T( GammaInteractionCalc::transmition_length_coefficient( ext_atten.material.get(), energyf ) / ext_atten.material->density );
    else
      mu = get_atten_coef_for_an( ext_atten.atomic_number, energyf );
    
    T areal_density = ext_atten.areal_density;

    // See note in self-atten branch above: we do NOT clamp `areal_density` to >= 0 (that would zero
    //  the Jet gradient across the deliberately-loosened bound margin - A2); e^{-mu*AD} is smooth for
    //  the tiny negative AD reachable inside the margin, and the converged value is clamped >= 0 at
    //  extraction. The -1e-3 g/cm^2 throw (in PhysicalUnits, below the -1e-5 loosening) only catches
    //  grossly negative values.
    assert( !isinf(areal_density) );
    if( (areal_density <= -1.0e-3 * PhysicalUnits::g_per_cm2)
        || isnan(areal_density) || isinf(areal_density) )
      throw std::runtime_error( "eval_physical_model_eqn: areal density must be >= 0." );

    assert( mu >= 0.0 );
    if( mu < 0.0 )
      mu = fmax(mu, 0.0);

    assert( areal_density >= -1.0e-3 * PhysicalUnits::g_per_cm2 );

    // Apply unconditionally: `mu >= 0` (clamped above) and `AD >= -1e-5` is a valid exponent.
    //  Guarding on `areal_density >= 0` would skip the term for AD in the loosened margin and
    //  re-zero the external-AD gradient there.
    answer *= exp( -mu * areal_density );
    
    assert( !isnan(answer) && !isinf(answer) );
  }//for( size_t i = 0; i < external_attens.size(); ++i )
  
  assert( !isnan(answer) && !isinf(answer) );
  
  return answer;
}//eval_physical_model_eqn_imp(...)
  

  template<typename T>
Pu242ByCorrelationOutput<T> correct_pu_mass_fractions_for_pu242( Pu242ByCorrelationInput<T> input, PuCorrMethod method )
{
  using namespace std;
  using namespace ceres;
  
  // First, lets normalize the the input relative mass, jic it isnt already
  const bool corr_for_age = (input.pu_age > 0.0);
  const SandiaDecay::SandiaDecayDataBase * const db = corr_for_age ? DecayDataBaseServer::database() : nullptr;
  const SandiaDecay::Nuclide * const pu238 = corr_for_age ? db->nuclide( "Pu238" ) : nullptr;
  const SandiaDecay::Nuclide * const pu239 = corr_for_age ? db->nuclide( "Pu239" ) : nullptr;
  const SandiaDecay::Nuclide * const pu240 = corr_for_age ? db->nuclide( "Pu240" ) : nullptr;
  const SandiaDecay::Nuclide * const pu241 = corr_for_age ? db->nuclide( "Pu241" ) : nullptr;
  const SandiaDecay::Nuclide * const pu242 = corr_for_age ? db->nuclide( "Pu242" ) : nullptr;

  if( corr_for_age )
  {
    //Need to back-decay to T=0, to apply correction.  We will then re-decay to the final age for the answer.
    assert( db );
    assert( pu238 && pu239 && pu240 && pu241 );

    vector<tuple<const SandiaDecay::Nuclide *,T>> input_nuclide_rel_acts;
    if( input.pu238_rel_mass > 0.0 )
      input_nuclide_rel_acts.emplace_back( pu238, input.pu238_rel_mass * pu238->activityPerGram() );
    if( input.pu239_rel_mass > 0.0 )
      input_nuclide_rel_acts.emplace_back( pu239, input.pu239_rel_mass * pu239->activityPerGram() );
    if( input.pu240_rel_mass > 0.0 )
      input_nuclide_rel_acts.emplace_back( pu240, input.pu240_rel_mass * pu240->activityPerGram() );
    if( input.pu240_rel_mass > 0.0 )
      input_nuclide_rel_acts.emplace_back( pu241, input.pu241_rel_mass * pu241->activityPerGram() );

    const vector<tuple<const SandiaDecay::Nuclide *,T,T>> time_zero_vals
                                  = back_decay_relative_activities( input.pu_age, input_nuclide_rel_acts );
    for( const tuple<const SandiaDecay::Nuclide *,T,T> &nuc_act_mass : time_zero_vals )
    {
      const SandiaDecay::Nuclide * const nuc = get<0>(nuc_act_mass);
      const T rel_act = get<1>(nuc_act_mass);
      if( nuc == pu238 )
        input.pu238_rel_mass = rel_act / nuc->activityPerGram();
      else if( nuc == pu239 )
        input.pu239_rel_mass = rel_act / nuc->activityPerGram();
      else if( nuc == pu240 )
        input.pu240_rel_mass = rel_act / nuc->activityPerGram();
      else if( nuc == pu241 )
        input.pu241_rel_mass = rel_act / nuc->activityPerGram();
      else { assert( 0 ); }
    }//for( loop over back-decayed nuclides )

    // We will just let `input.other_pu_mass` stay the same - we dont expect to use it anyway
  }//if( input.pu_age > 0.0 )

  T sum_input_mass = T(0.0);
  sum_input_mass += input.pu238_rel_mass;
  sum_input_mass += input.pu239_rel_mass;
  sum_input_mass += input.pu240_rel_mass;
  sum_input_mass += input.pu241_rel_mass;
  sum_input_mass += input.other_pu_mass;

  input.pu238_rel_mass /= sum_input_mass;
  input.pu239_rel_mass /= sum_input_mass;
  input.pu240_rel_mass /= sum_input_mass;
  input.pu241_rel_mass /= sum_input_mass;
  input.other_pu_mass  /= sum_input_mass;

  T pu242_mass_frac = T(0.0), fractional_uncert = T(0.0);
  switch( method )
  {
    case PuCorrMethod::Bignan95_PWR:
    case PuCorrMethod::Bignan95_BWR:
    {
      // Note: Bignan 98 provides a nice "database" in Table 1 that could be used to improve
      //  the value of c_0 used, based on ratios of Pu238, Pu240, and Pu242 to Pu239, but
      //  realistically this is way to fine-meshed for the calculations we could hope to do
      //  in InterSpec
      
      const T c_0 = ((method == PuCorrMethod::Bignan95_PWR) ? T(1.313) : T(1.117));
      const T pu238_to_pu239 = input.pu238_rel_mass / input.pu239_rel_mass;
      const T pu240_to_pu239 = input.pu240_rel_mass / input.pu239_rel_mass;
      const T pu242_to_pu239 = c_0 * pow( pu238_to_pu239, T(0.33) ) * pow( pu240_to_pu239, T(1.7) );
      
      pu242_mass_frac = pu242_to_pu239 * input.pu239_rel_mass;
      break;
    }//case PuCorrMethod::Bignan95_BWR:
    
    case PuCorrMethod::ByPu239Only:
    {
      const T A = T(9.66E-3);
      const T C = T(-3.83);
      
      pu242_mass_frac = A * pow( input.pu239_rel_mass, C );
      break;
    }//case PuCorrMethod::ByPu239Only:
      
    case PuCorrMethod::NotApplicable:
      pu242_mass_frac = T(0.0);
      break;
  }//switch( method )
  
  
  Pu242ByCorrelationOutput<T> answer;
  
  // We need to correct for the Pu242 mass fraction.
  //  See equation 8-14 (page 249) in:
  //    "Plutonium Isotopic Composition by Gamma-Ray Spectroscopy"
  //    T. E. Sampson
  // https://www.lanl.gov/org/ddste/aldgs/sst-training/_assets/docs/PANDA/Plutonium%20Isotopic%20Composition%20by%20Gamma-Ray%20Spectroscopy%20Ch.%208%20p.%20221-272.pdf
  
  answer.pu238_mass_frac = input.pu238_rel_mass * (T(1.0) - pu242_mass_frac);
  answer.pu239_mass_frac = input.pu239_rel_mass * (T(1.0) - pu242_mass_frac);
  answer.pu240_mass_frac = input.pu240_rel_mass * (T(1.0) - pu242_mass_frac);
  answer.pu241_mass_frac = input.pu241_rel_mass * (T(1.0) - pu242_mass_frac);
  answer.pu242_mass_frac = pu242_mass_frac;

  if( input.pu_age > 0.0 )
  {
    answer.pu238_mass_frac *= exp( -input.pu_age * pu238->decayConstant() );
    answer.pu239_mass_frac *= exp( -input.pu_age * pu239->decayConstant() );
    answer.pu240_mass_frac *= exp( -input.pu_age * pu240->decayConstant() );
    answer.pu241_mass_frac *= exp( -input.pu_age * pu241->decayConstant() );
    answer.pu242_mass_frac *= exp( -input.pu_age * pu242->decayConstant() );

    const T norm_amount = answer.pu238_mass_frac
                        + answer.pu239_mass_frac
                        + answer.pu240_mass_frac
                        + answer.pu241_mass_frac
                        + answer.pu242_mass_frac;

    answer.pu238_mass_frac /= norm_amount;
    answer.pu239_mass_frac /= norm_amount;
    answer.pu240_mass_frac /= norm_amount;
    answer.pu241_mass_frac /= norm_amount;
    answer.pu242_mass_frac /= norm_amount;
  }//if( input.pu_age > 0.0 )

  
  switch( method )
  {
    case PuCorrMethod::Bignan95_PWR:
    case PuCorrMethod::Bignan95_BWR:
    {
      const T pu238_pu239 = answer.pu238_mass_frac / answer.pu239_mass_frac;
      const T pu240_pu239 = answer.pu240_mass_frac / answer.pu239_mass_frac;
      const T pu242_pu239 = answer.pu242_mass_frac / answer.pu239_mass_frac;
      
      answer.is_within_range = ((pu238_pu239 >= T(0.007851)) && (pu238_pu239 <= T(0.02952)))
                                && ((pu240_pu239 >= T(0.2688)) && (pu240_pu239 <= T(0.4586)))
                                && ((pu242_pu239 >= T(0.03323)) && (pu242_pu239 <= T(0.1152)));
      
      if( method == PuCorrMethod::Bignan95_PWR )
        answer.pu242_uncert = T(0.03);
      else
        answer.pu242_uncert = T(0.07);
      
      break;
    }//case PuCorrMethod::Bignan95_BWR:
      
    case PuCorrMethod::ByPu239Only:
    {
      answer.is_within_range = (answer.pu239_mass_frac >= T(0.55)) && (answer.pu239_mass_frac <= T(0.80));
      
      if( (answer.pu239_mass_frac >= T(0.55)) && (answer.pu239_mass_frac <= T(0.64)) )
        answer.pu242_uncert = T(0.012);
      else if( answer.pu239_mass_frac < T(0.55) )
        answer.pu242_uncert = T(0.05); //totally made up
      else if( answer.pu239_mass_frac < T(0.70) )
        answer.pu242_uncert = T(0.01);
      else if( answer.pu239_mass_frac < T(0.80) )
        answer.pu242_uncert = T(0.04);
      else
        answer.pu242_uncert = T(0.05); //totally made up
      
      break;
    }//case PuCorrMethod::ByPu239Only
      
    case PuCorrMethod::NotApplicable:
      answer.is_within_range = true;
      answer.pu242_uncert = T(0.0);
    break;
  }//switch( method )

  // Except for the totally made up uncertainties (which are out of validated ranges), the
  //  actual errors are likely much larger than reported in the paper, as their data probably
  //  came from similar sources, or at least dont include nearly all the ways Pu is made, so
  //  we'll throw an arbitrary factor of 2 onto the uncertainty.
  const T engineering_uncert_multiple = T(2.0);
  answer.pu242_uncert *= engineering_uncert_multiple;

#ifndef NDEBUG
  const T total_pu = answer.pu238_mass_frac
    + answer.pu239_mass_frac
    + answer.pu240_mass_frac
    + answer.pu241_mass_frac
    + answer.pu242_mass_frac;

  assert( (total_pu > 0.999) && (total_pu < 1.001) );
#endif //NDEBUG

  return answer;
}//correct_pu_mass_fractions_for_pu242( ... )

template<typename T>
std::vector<std::tuple<const SandiaDecay::Nuclide *,T,T>> back_decay_relative_activities(
    const T back_decay_time,
    std::vector<std::tuple<const SandiaDecay::Nuclide *,T>> &nuclide_rel_acts )
{
  using namespace std;
  using namespace ceres;
  
  assert( back_decay_time >= T(0.0) );
  std::vector<std::tuple<const SandiaDecay::Nuclide *,T,T>> answer;

  T rel_mass_sum = T(0.0);
  for( const auto &nuc_activity : nuclide_rel_acts )
  {
    const SandiaDecay::Nuclide * const nuc = std::get<0>(nuc_activity);
    const T final_act = std::get<1>(nuc_activity);
    assert( nuc );
    if( !nuc || (final_act <= T(0.0)) )
      continue;

    const T decrease_factor = exp( -back_decay_time * static_cast<T>(nuc->decayConstant()) );
    const T initial_activity = final_act / decrease_factor;

    answer.emplace_back( nuc, initial_activity, T(0.0) );

    const T initial_rel_mass = initial_activity / static_cast<T>(nuc->activityPerGram());
    rel_mass_sum += initial_rel_mass;
  }

  for( size_t i = 0; i < answer.size(); ++i )
  {
    const SandiaDecay::Nuclide * const nuc = std::get<0>(answer[i]);
    const T initial_activity = std::get<1>(answer[i]);
    const T initial_rel_mass = initial_activity / static_cast<T>(nuc->activityPerGram());

    std::get<2>(answer[i]) = initial_rel_mass / rel_mass_sum;
  }

    return answer;
}

}//namespace RelActCalc

#endif //RelActCalc_imp_hpp
