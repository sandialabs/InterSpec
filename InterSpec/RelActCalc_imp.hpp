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
#include <optional>
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
    

/** This function is the equivalent of `mass_attenuation_coef(...)`, but
   I blindly *think/hope/wish* it preserve the ceres::Jet derivative information.
   TODO: check that this is actually proper way to propagate derivate info for fitting AN
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
  
  assert( (an_scalar >= 1.0) && (an_scalar <= 98.0) );
  
  const int lower_an = std::max( 1, static_cast<int>( std::floor(an_scalar) ) );
  const int sign = (lower_an < 98) ? 1 : -1;
  const int upper_an = lower_an + sign;

  const double lower_mu = MassAttenuation::massAttenuationCoefficientElement(lower_an, energy);
  const double upper_mu = MassAttenuation::massAttenuationCoefficientElement(upper_an, energy);
  const T anfrac = an - static_cast<double>(lower_an);  //Looks like this preserves the derivative
  const T mu = (1.0 - anfrac)*lower_mu + anfrac*upper_mu;
  
  return mu * static_cast<double>(sign);
}//T get_atten_coef_for_an( const T &an )

    
template<typename T>
T eval_physical_model_eqn_imp( const double energy,
                                std::optional<RelActCalc::PhysModelShield<T>> self_atten,
                                const std::vector<RelActCalc::PhysModelShield<T>> &external_attens,
                                const DetectorPeakResponse * const drf,
                                std::optional<T> b,
                                std::optional<T> c )
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
  
  if( b.has_value() && c.has_value() )
  {
    const T b_val = *b;
    const T c_val = *c;
    const T b_part = pow(0.001*energy, b_val);
    const T c_part = pow(c_val, 1000.0/energy);
    answer = b_part * c_part;
    
    assert( !isnan(answer) && !isinf(answer) );
    if( isnan(answer) || isinf(answer) )
      throw std::logic_error( "hoerl_b or hoerl_c gives eqn NaN or Inf" );
  }
  
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
    
    assert( (areal_density >= -1.0E-3) && !isinf(areal_density) );
    if( (areal_density <= -1.0E-3) || isnan(areal_density) || isinf(areal_density) )
      throw std::runtime_error( "eval_physical_model_eqn: areal density must be >= 0." );
    
    if( areal_density < 0.0 )
      areal_density = fmax(areal_density, 0.0);
    
    if( (mu > 0.0) && (areal_density > 0.0) )
      answer *= (1.0 - exp(-mu * areal_density)) / (mu * self_atten->areal_density);
    
    assert( !isnan(answer) && !isinf(answer) );
  }//if( self_atten->has_value() )
  
  for( const RelActCalc::PhysModelShield<T> &ext_atten : external_attens )
  {
    // TODO: `GammaInteractionCalc::transmition_length_coefficient` can be a real bottleneck of computation - at soem point we should memoise its results
    T mu;
    if( ext_atten.material )
      mu = T( GammaInteractionCalc::transmition_length_coefficient( ext_atten.material.get(), energyf ) / ext_atten.material->density );
    else
      mu = get_atten_coef_for_an( ext_atten.atomic_number, energyf );
    
    T areal_density = ext_atten.areal_density;
    
    assert( (areal_density >= -1.0E-3) && !isinf(areal_density) );
    if( (areal_density <= -1.0E-3) || isnan(areal_density) || isinf(areal_density) )
      throw std::runtime_error( "eval_physical_model_eqn: areal density must be >= 0." );
    
    if( areal_density < 0.0 )
      areal_density = fmax(areal_density, 0.0);
    
    if( (mu > 0.0) && (ext_atten.areal_density > 0.0) )
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
