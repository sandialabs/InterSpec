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

#include <optional>
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/GammaInteractionCalc.h"

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
  
}//namespace RelActCalc

#endif //RelActCalc_imp_hpp
