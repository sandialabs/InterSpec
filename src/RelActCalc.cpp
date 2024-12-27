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

#include <cmath>
#include <string>
#include <vector>
#include <assert.h>
#include <iostream> //for cout, only for debug
#include <cstring>
#include <stdexcept>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/MaterialDB.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"

using namespace std;

namespace 
{
  /** Make "Pu (plutonium)" into "Pu", for printing PhysicalModelShieldInput materials. */
  string cleanup_mat_name( string text )
  {
    size_t first = text.find('(');
    if( first != std::string::npos )
    {
      size_t second = text.find(')', first);
      if( second != std::string::npos )
        text.erase(first, second - first + 1);
    }
    
    SpecUtils::trim( text );
    return text;
  };//cleanup_mat_name(...)
}//namespace 

namespace RelActCalc
{
const double PhysicalModelShieldInput::sm_upper_allowed_areal_density_in_g_per_cm2 = 500.0;

const char *to_str( const RelEffEqnForm form )
{
  switch( form )
  {
    case RelEffEqnForm::LnX:    return "LnX";
    case RelEffEqnForm::LnY:    return "LnY";
    case RelEffEqnForm::LnXLnY: return "LnXLnY";
    case RelEffEqnForm::FramEmpirical: return "FRAM Empirical";
    case RelEffEqnForm::FramPhysicalModel: return "FRAM Physical";
      
  }
  
  assert( 0 );
  throw runtime_error( "to_str(RelEffEqnForm): invalid input" );
  return "";
}//to_str(...)


RelEffEqnForm rel_eff_eqn_form_from_str( const char *str )
{
  // Not sure of a good way to get a warning if we change RelEffEqnForm enum arbitrarily
  const RelEffEqnForm eqn_forms[] = {
    RelEffEqnForm::LnX,
    RelEffEqnForm::LnY,
    RelEffEqnForm::LnXLnY,
    RelEffEqnForm::FramEmpirical,
    RelEffEqnForm::FramPhysicalModel
  };
  
  for( const RelEffEqnForm eqn : eqn_forms )
  {
    const char *this_eqn_str = to_str( eqn );
    if( SpecUtils::iequals_ascii( this_eqn_str, str ) )
      return eqn;
  }
  
  throw runtime_error( "String '" + std::string(str) + "' not a valid RelEffEqnForm" );
}//rel_eff_eqn_rorm_from_str(...)


string rel_eff_eqn_text( const RelEffEqnForm eqn_form, const std::vector<double> &coefs )
{
  const size_t nsigfig = 5;
  string rel_eff_eqn_str;
  switch( eqn_form )
  {
    case RelActCalc::RelEffEqnForm::LnX:
    {
      for( size_t i = 0; i < coefs.size(); ++i )
      {
        const auto val = coefs[i];
        if( i == 0 )
        {
          rel_eff_eqn_str += SpecUtils::printCompact(val,nsigfig);
        }else
        {
          rel_eff_eqn_str += (val < 0.0 ? " - " : " + " )
                             +  SpecUtils::printCompact(fabs(val),nsigfig)
                             + "*ln(x)^" + to_string(i);
        }
      }//for( size_t i = 0; i < coefs.size(); ++i )
      
      break;
    }//case RelActCalc::RelEffEqnForm::LnX:
      
    case RelActCalc::RelEffEqnForm::LnY:
    {
      //y = exp( a + b*x + c/x + d/x^2 + e/x^3 + ... )
      rel_eff_eqn_str += "exp( ";
      
      for( size_t i = 0; i < coefs.size(); ++i )
      {
        const auto val = coefs[i];
        if( i == 0 )
        {
          rel_eff_eqn_str += SpecUtils::printCompact(val,nsigfig);
        }else
        {
          rel_eff_eqn_str += (val < 0.0 ? " - " : " + " )
                             + SpecUtils::printCompact(fabs(val),nsigfig);
          if( i == 1 )
            rel_eff_eqn_str += "*x";
          else if( i == 2 )
            rel_eff_eqn_str += "/x";
          else
            rel_eff_eqn_str += "/x^" +  to_string( i - 1 );
        }
      }//for( size_t i = 0; i < coefs.size(); ++i )
      rel_eff_eqn_str += " )";
      
      break;
    }//case RelActCalc::RelEffEqnForm::LnY:
      
    case RelActCalc::RelEffEqnForm::LnXLnY:
    {
      rel_eff_eqn_str += "exp( ";
      
      for( size_t i = 0; i < coefs.size(); ++i )
      {
        const auto val = coefs[i];
        if( i == 0 )
        {
          rel_eff_eqn_str += to_string(val);
        }else
        {
          rel_eff_eqn_str += (val < 0.0 ? " - " : " + " )
                             + SpecUtils::printCompact(fabs(val),nsigfig)
                             + "*ln(x)^" + to_string(i);
        }
      }//for( size_t i = 0; i < coefs.size(); ++i )
      
      rel_eff_eqn_str += " )";
      break;
    }//case RelActCalc::RelEffEqnForm::LnXLnY:
      
    case RelActCalc::RelEffEqnForm::FramEmpirical:
    {
      rel_eff_eqn_str += "exp( ";
      
      for( size_t i = 0; i < coefs.size(); ++i )
      {
        const auto val = coefs[i];
        if( i == 0 )
        {
          rel_eff_eqn_str += SpecUtils::printCompact(val,nsigfig);
        }else if( i == 1 )
        {
          rel_eff_eqn_str += (val < 0.0 ? " - " : " + " )
                             + SpecUtils::printCompact(fabs(val),nsigfig) + "/(x*x)";
        }else
        {
          rel_eff_eqn_str += (val < 0.0 ? " - " : " + " )
                             + SpecUtils::printCompact(fabs(val),nsigfig)
                             + "*ln(x)^" + to_string(i-1);
        }
      }//for( size_t i = 0; i < coefs.size(); ++i )
      
      rel_eff_eqn_str += " )";
      
      break;
    }//case RelActCalc::RelEffEqnForm::FramEmpirical:
      
    case RelActCalc::RelEffEqnForm::FramPhysicalModel:
      throw runtime_error( "rel_eff_eqn_text: can not be called for FramPhysicalModel" );
  }//switch( answer.m_eqn_form )
  
  return rel_eff_eqn_str;
}//rel_eff_eqn_text(...)


std::string rel_eff_eqn_js_function( const RelEffEqnForm eqn_form, const std::vector<double> &coefs )
{
  string rel_eff_fcn;
  rel_eff_fcn = "function(x){ return ";
  switch( eqn_form )
  {
    case RelActCalc::RelEffEqnForm::LnX:
    {
      //y = a + b*ln(x) + c*(ln(x))^2 + d*(ln(x))^3 + ...
      for( size_t i = 0; i < coefs.size(); ++i )
      {
        if( i == 0 )
          rel_eff_fcn += to_string(coefs[i]);
        else if( i == 1 )
          rel_eff_fcn += " + " + to_string(coefs[i]) + "*Math.log(x)";
        else
          rel_eff_fcn += " + " + to_string(coefs[i]) + "*Math.pow( Math.log(x), " + to_string(i) + " )";
      }//for( size_t i = 0; i < coefs.size(); ++i )
      
      break;
    }//case RelActCalc::RelEffEqnForm::LnX:
      
    case RelActCalc::RelEffEqnForm::LnY:
    {
      //y = exp( a + b*x + c/x + d/x^2 + e/x^3 + ... )
      rel_eff_fcn += "Math.exp( ";
      
      for( size_t i = 0; i < coefs.size(); ++i )
      {
        rel_eff_fcn += (i ? " + " : "") + to_string(coefs[i]);
        if( i == 1 )
          rel_eff_fcn += "*x";
        else if( i > 1 )
          rel_eff_fcn += "/Math.pow(x," + to_string(i-1) + ")";
      }//for( size_t i = 0; i < coefs.size(); ++i )
      
      rel_eff_fcn += " )";
      break;
    }//case RelActCalc::RelEffEqnForm::LnY:
      
    case RelActCalc::RelEffEqnForm::LnXLnY:
    {
      //y = exp(a  + b*(lnx) + c*(lnx)^2 + d*(lnx)^3 + ... )
      
      rel_eff_fcn += "Math.exp( ";
      
      for( size_t i = 0; i < coefs.size(); ++i )
      rel_eff_fcn += (i ? " + " : "") + to_string(coefs[i]) + "*Math.pow( Math.log(x), " + to_string(i) + " )";
      
      rel_eff_fcn += " )";
      break;
    }//case RelActCalc::RelEffEqnForm::LnXLnY:
      
    case RelActCalc::RelEffEqnForm::FramEmpirical:
    {
      rel_eff_fcn += "Math.exp( ";
      
      for( size_t i = 0; i < coefs.size(); ++i )
      {
        rel_eff_fcn += (i ? " + " : "") + to_string(coefs[i]);
        if( i == 1 )
          rel_eff_fcn += "/(x*x)";
        else if( i > 1 )
          rel_eff_fcn += "*Math.pow( Math.log(x)," + to_string(i-1) + ")";
      }
      
      rel_eff_fcn += " )";
      break;
      
    case RelActCalc::RelEffEqnForm::FramPhysicalModel:
      throw runtime_error( "rel_eff_eqn_js_function: can not be called for FramPhysicalModel" );
    }//case RelActCalc::RelEffEqnForm::FramEmpirical:
  }//switch( answer.m_eqn_form )
  
  rel_eff_fcn += "; }";
  
  return rel_eff_fcn;
}//rel_eff_eqn_js_function(...);


double eval_eqn( const double energy, const RelEffEqnForm eqn_form, const vector<double> &coeffs )
{
  return eval_eqn( energy, eqn_form, &(coeffs[0]), coeffs.size() );
}


double eval_eqn( const double energy, const RelEffEqnForm eqn_form,
                const double * const coeffs, const size_t num_coefs )
{
  if( energy <= 0.0 )
    throw runtime_error( "eval_eqn: energy must be greater than zero." );
  
  if( num_coefs < 1 )
    throw runtime_error( "eval_eqn: need at least one coefficients." );
  
  
  double answer = 0.0;
  
  switch( eqn_form )
  {
    case RelEffEqnForm::LnX:
    {
      // y = a + b*ln(x) + c*(ln(x))^2 + d*(ln(x))^3 + ...
      const double log_energy = std::log(energy);
      
      for( size_t order = 0; order < num_coefs; ++order )
      {
        switch( order )
        {
          case 0:  answer += coeffs[order];                                         break;
          case 1:  answer += coeffs[order] * log_energy;                            break;
          default: answer += coeffs[order] * std::pow( log_energy, (double)order ); break;
        }//switch( order )
      }//for( loop over coeffs )
      
      break;
    }//case RelEffEqnForm::LnX:
      
    case RelEffEqnForm::LnY:
    {
      // y = exp( a + b*x + c/x + d/x^2 + e/x^3 + ... )
      // Note that it would be a little more convenient to have y = exp(a*x + b + c/x + ...), but
      //  I want to keep the leading term independent of energy.
      for( size_t order = 0; order < num_coefs; ++order )
      {
        switch( order )
        {
          case 0:  answer += coeffs[order];                                   break;
          case 1:  answer += coeffs[order] * energy;                          break;
          case 2:  answer += coeffs[order] / energy;                          break;
          default: answer += coeffs[order] / std::pow( energy, order - 1.0 ); break;
        }//switch( order )
      }//for( loop over coeffs )
      
      answer = std::exp( answer );
      
      break;
    }//case RelEffEqnForm::LnY:
      
    case RelEffEqnForm::LnXLnY:
    {
      // y = exp (a  + b*(lnx) + c*(lnx)^2 + d*(lnx)^3 + ... )
      const double log_energy = std::log(energy);
      
      for( size_t order = 0; order < num_coefs; ++order )
      {
        switch( order )
        {
          case 0:  answer += coeffs[order]; break;
          case 1:  answer += coeffs[order] * log_energy; break;
          default: answer += coeffs[order] * std::pow( log_energy, (double)order ); break;
        }//switch( order )
      }//for( loop over coeffs )
      
      answer = std::exp( answer );
      
      break;
    }//case RelEffEqnForm::LnXLnY:
      
    case RelEffEqnForm::FramEmpirical:
    {
      // y = exp( a + b/x^2 + c*(lnx) + d*(lnx)^2 + e*(lnx)^3 )
      const double log_energy = std::log(energy);
      
      for( size_t order = 0; order < num_coefs; ++order )
      {
        switch( order )
        {
          case 0:  answer += coeffs[order];  break;
          case 1:  answer += coeffs[order] / (energy*energy); break;
          default: answer += coeffs[order] * std::pow( log_energy, order - 1.0 ); break;
        }//switch( order )
      }//for( loop over coeffs )
      
      answer = std::exp( answer );
      
      break;
    }//case RelEffEqnForm::FramEmpirical:
      
    case RelEffEqnForm::FramPhysicalModel:
      throw runtime_error( "RelActCalc::eval_eqn() can not be called for FramPhysicalModel." );
  }//switch( eqn_form )
  
  return answer;
}//eval_eqn(...)


double eval_eqn_uncertainty( const double energy, const RelEffEqnForm eqn_form,
                            const std::vector<std::vector<double>> &covariance )
{
  if( covariance.empty() )
    throw runtime_error( "eval_eqn_uncertainty: empty coefficients passed in." );
  
  switch( eqn_form )
  {
    case RelEffEqnForm::LnX:
    {
      // y = a + b*ln(x) + c*(ln(x))^2 + d*(ln(x))^3 + ...
      const double log_energy = std::log(energy);
      
      double uncert_sq = 0.0;
      
      for( size_t i = 0; i < covariance.size(); ++i )
      {
        assert( covariance[i].size() == covariance.size() );
        if( covariance[i].size() != covariance.size() )  //JIC for release builds
          throw runtime_error( "eval_eqn_uncertainty: covariance not a square matrix." );
        
        for( size_t j = 0; j <= covariance.size(); ++j )
        uncert_sq += std::pow(log_energy,1.0*i) * covariance[i][j] * std::pow(log_energy,1.0*j);
      }//for( size_t i = 0; i < coefs.size(); ++i )
      
      assert( uncert_sq >= 0.0 );
      
      return sqrt(uncert_sq);
    }//case RelEffEqnForm::LnX:
      
      
    case RelEffEqnForm::LnY:
    {
      // y = exp( a + b*x + c/x + d/x^2 + e/x^3 + ... )
      const auto eval_term_log = [energy]( const size_t order ) -> double {
        switch( order )
        {
          case 0:  return 1.0;
          case 1:  return energy;
          default:
            break;
        }//switch( order )
        
        return std::pow( energy, 1.0 - order );
      };//eval_term_log
      
      double log_uncert_sq = 0.0;
      
      for( size_t i = 0; i < covariance.size(); ++i )
      {
        assert( covariance[i].size() == covariance.size() );
        if( covariance[i].size() != covariance.size() )  //JIC for release builds
          throw runtime_error( "eval_eqn_uncertainty: covariance not a square matrix." );
        
        for( size_t j = 0; j <= covariance.size(); ++j )
        log_uncert_sq += eval_term_log(i) * covariance[i][j] * eval_term_log(j);
      }//for( size_t i = 0; i < coefs.size(); ++i )
      
      assert( log_uncert_sq >= 0.0 );
      
      return exp( sqrt(log_uncert_sq) );
    }//case RelEffEqnForm::LnY:
      
    case RelEffEqnForm::LnXLnY:
    {
      // y = exp(a  + b*(lnx) + c*(lnx)^2 + d*(lnx)^3 + ... )
      const double log_energy = std::log(energy);
      
      double log_uncert_sq = 0.0;
      
      for( size_t i = 0; i < covariance.size(); ++i )
      {
        assert( covariance[i].size() == covariance.size() );
        if( covariance[i].size() != covariance.size() )  //JIC for release builds
          throw runtime_error( "eval_eqn_uncertainty: covariance not a square matrix." );
        
        for( size_t j = 0; j <= covariance.size(); ++j )
        log_uncert_sq += std::pow(log_energy,1.0*i) * covariance[i][j] * std::pow(log_energy,1.0*j);
      }//for( size_t i = 0; i < coefs.size(); ++i )
      
      assert( log_uncert_sq >= 0.0 );
      
      return exp( sqrt(log_uncert_sq) );
    }//case RelEffEqnForm::LnXLnY:
      
    case RelEffEqnForm::FramEmpirical:
    {
      // y = exp( a + b/x^2 + c*(lnx) + d*(lnx)^2 + e*(lnx)^3 )
      const double log_energy = std::log(energy);
      
#ifdef _MSC_VER
#pragma message( "eval_eqn_uncertainty with RelEffEqnForm::FramEmpirical not tested!" )
#else
#warning "eval_eqn_uncertainty with RelEffEqnForm::FramEmpirical not tested!"
#endif
      
      double log_uncert_sq = 0.0;
      
      for( size_t i = 0; i < covariance.size(); ++i )
      {
        assert( covariance[i].size() == covariance.size() );
        if( covariance[i].size() != covariance.size() )  //JIC for release builds
          throw runtime_error( "eval_eqn_uncertainty: covariance not a square matrix." );
        
        double i_component = 0.0;
        switch( i )
        {
          case 0:  i_component = 1.0;  break;
          case 1:  i_component = 1.0 / (energy*energy); break;
          default: i_component = std::pow(log_energy, i - 1.0); break;
        }//switch( i )
        
        for( size_t j = 0; j <= covariance.size(); ++j )
        {
          double j_component = 0.0;
          switch( j )
          {
            case 0:  j_component = 1.0;  break;
            case 1:  j_component = 1.0 / (energy*energy); break;
            default: j_component = std::pow(log_energy, j - 1.0); break;
          }//switch( i )
          
          log_uncert_sq += i_component * covariance[i][j] * j_component;
        }
      }//for( size_t i = 0; i < coefs.size(); ++i )
      
      assert( log_uncert_sq >= 0.0 );
      
      return exp( sqrt(log_uncert_sq) );
    }//case RelEffEqnForm::FramEmpirical:
      
    case RelEffEqnForm::FramPhysicalModel:
      assert( 0 );
      throw runtime_error( "RelActCalc::eval_eqn_uncertainty() can not be called for FramPhysicalModel." );
  }//switch( eqn_form )
  
  
  assert( 0 );   //shouldnt ever get here
  return -999.9; //avoid compiler warning
}//double eval_eqn_uncertainty(...)


const std::string &to_str( const PuCorrMethod method )
{
  const static std::string s_Bignan95_PWR{ "Bignan95_PWR" };
  const static std::string s_Bignan95_BWR{ "Bignan95_BWR" };
  const static std::string s_ByPu239Only{ "ByPu239Only" };
  const static std::string s_NotApplicable{ "NotApplicable" };
  
  switch( method )
  {
    case PuCorrMethod::Bignan95_PWR:  return s_Bignan95_PWR;
    case PuCorrMethod::Bignan95_BWR:  return s_Bignan95_BWR;
    case PuCorrMethod::ByPu239Only:   return s_ByPu239Only;
    case PuCorrMethod::NotApplicable: return s_NotApplicable;
  }//switch( method )
  
  assert( 0 );
  throw std::logic_error( "invalid PuCorrMethod" );
  
  return s_NotApplicable;
}//const std::string &to_str( const PuCorrMethod method )


const std::string &to_description( const PuCorrMethod method )
{
  const static std::string s_desc_Bignan95_PWR{ "By 238,239,240 ratio PWR" };
  const static std::string s_desc_Bignan95_BWR{ "By 238,239,240 ratio BWR" };
  const static std::string s_desc_ByPu239Only{ "By Pu239 Frac." };
  const static std::string s_desc_NotApplicable{ "None" };
  
  switch( method )
  {
    case PuCorrMethod::Bignan95_PWR:  return s_desc_Bignan95_PWR;
    case PuCorrMethod::Bignan95_BWR:  return s_desc_Bignan95_BWR;
    case PuCorrMethod::ByPu239Only:   return s_desc_ByPu239Only;
    case PuCorrMethod::NotApplicable: return s_desc_NotApplicable;
  }//switch( method )
  
  assert( 0 );
  throw std::logic_error( "invalid PuCorrMethod" );
  
  return s_desc_NotApplicable;
}


Pu242ByCorrelationOutput correct_pu_mass_fractions_for_pu242( Pu242ByCorrelationInput input, PuCorrMethod method )
{
  // First, lets normalize the the input relative mass, jic it isnt already
  float sum_input_mass = 0.0f;
  sum_input_mass += input.pu238_rel_mass;
  sum_input_mass += input.pu239_rel_mass;
  sum_input_mass += input.pu240_rel_mass;
  sum_input_mass += input.pu241_rel_mass;
  sum_input_mass += input.other_pu_mass;
  
  // TODO: in principle we want to account for the decay of Am241 and Pu241 better, but for now we'll just be really gross about it and equate the two nuclides
  sum_input_mass += input.am241_rel_mass;
  
  input.pu238_rel_mass /= sum_input_mass;
  input.pu239_rel_mass /= sum_input_mass;
  input.pu240_rel_mass /= sum_input_mass;
  input.pu241_rel_mass /= sum_input_mass;
  input.am241_rel_mass /= sum_input_mass;
  input.other_pu_mass  /= sum_input_mass;
  
  double pu242_mass_frac = 0.0, fractional_uncert = 0.0;
  switch( method )
  {
    case PuCorrMethod::Bignan95_PWR:
    case PuCorrMethod::Bignan95_BWR:
    {
      // Note: Bignan 98 provides a nice "database" in Table 1 that could be used to improve
      //  the value of c_0 used, based on ratios of Pu238, Pu240, and Pu242 to Pu239, but
      //  realistically this is way to fine-meshed for the calculations we could hope to do
      //  in InterSpec
      
      const double c_0 = ((method == PuCorrMethod::Bignan95_PWR) ? 1.313 : 1.117);
      const double pu238_to_pu239 = input.pu238_rel_mass / input.pu239_rel_mass;
      const double pu240_to_pu239 = input.pu240_rel_mass / input.pu239_rel_mass;
      const double pu242_to_pu239 = c_0 * std::pow( pu238_to_pu239, 0.33 ) * std::pow( pu240_to_pu239, 1.7 );
      
      pu242_mass_frac = pu242_to_pu239 * input.pu239_rel_mass;
      break;
    }//case PuCorrMethod::Bignan95_BWR:
    
    case PuCorrMethod::ByPu239Only:
    {
      const double A = 9.66E-3;
      const double C = -3.83;
      
      pu242_mass_frac = A * std::pow( input.pu239_rel_mass, C );
      break;
    }//case PuCorrMethod::ByPu239Only:
      
    case PuCorrMethod::NotApplicable:
      pu242_mass_frac = 0.0;
      break;
  }//switch( method )
  
  
  Pu242ByCorrelationOutput answer;
  
  // We need to correct for the Pu242 mass fraction.
  //  See equation 8-14 (page 249) in:
  //    "Plutonium Isotopic Composition by Gamma-Ray Spectroscopy"
  //    T. E. Sampson
  // https://www.lanl.gov/org/ddste/aldgs/sst-training/_assets/docs/PANDA/Plutonium%20Isotopic%20Composition%20by%20Gamma-Ray%20Spectroscopy%20Ch.%208%20p.%20221-272.pdf
  
  answer.pu238_mass_frac = input.pu238_rel_mass * (1.0 - pu242_mass_frac);
  answer.pu239_mass_frac = input.pu239_rel_mass * (1.0 - pu242_mass_frac);
  answer.pu240_mass_frac = input.pu240_rel_mass * (1.0 - pu242_mass_frac);
  answer.pu241_mass_frac = input.pu241_rel_mass * (1.0 - pu242_mass_frac);
  answer.am241_mass_frac = input.am241_rel_mass * (1.0 - pu242_mass_frac);
  answer.pu242_mass_frac = pu242_mass_frac;
  
  
  switch( method )
  {
    case PuCorrMethod::Bignan95_PWR:
    case PuCorrMethod::Bignan95_BWR:
    {
      const double pu238_pu239 = answer.pu238_mass_frac / answer.pu239_mass_frac;
      const double pu240_pu239 = answer.pu240_mass_frac / answer.pu239_mass_frac;
      const double pu242_pu239 = answer.pu242_mass_frac / answer.pu239_mass_frac;
      
      answer.is_within_range = ((pu238_pu239 >= 0.007851) && (pu238_pu239 <= 0.02952))
                                && ((pu240_pu239 >= 0.2688) && (pu240_pu239 <= 0.4586))
                                && ((pu242_pu239 >= 0.03323) && (pu242_pu239 <= 0.1152));
      
      if( method == PuCorrMethod::Bignan95_PWR )
        answer.pu242_uncert = 0.03;
      else
        answer.pu242_uncert = 0.07;
      
      break;
    }//case PuCorrMethod::Bignan95_BWR:
      
    case PuCorrMethod::ByPu239Only:
    {
      answer.is_within_range = (answer.pu239_mass_frac >= 0.55) && (answer.pu239_mass_frac <= 0.80);
      
      if( (answer.pu239_mass_frac >= 0.55) && (answer.pu239_mass_frac <= 0.64) )
        answer.pu242_uncert = 0.012;
      else if( answer.pu239_mass_frac < 0.55 )
        answer.pu242_uncert = 0.05; //totally made up
      else if( answer.pu239_mass_frac < 0.70 )
        answer.pu242_uncert = 0.01;
      else if( answer.pu239_mass_frac < 0.80 )
        answer.pu242_uncert = 0.04;
      else
        answer.pu242_uncert = 0.05; //totally made up
      
      break;
    }//case PuCorrMethod::ByPu239Only
      
    case PuCorrMethod::NotApplicable:
      answer.is_within_range = true;
      answer.pu242_uncert = 0.0;
    break;
  }//switch( method )
  
  
  // Except for the totally made up uncertainties (which are out of validated ranges), the
  //  actual errors are likely much larger than reported in the paper, as their data probably
  //  came from similar sources, or at least dont include nearly all the ways Pu is made, so
  //  we'll throw an arbitrary factor of 2 onto the uncertainty.
  const double engineering_uncert_multiple = 2.0;
  answer.pu242_uncert *= engineering_uncert_multiple;
  
  return answer;
}//correct_pu_mass_fractions_for_pu242( ... )



void test_pu242_by_correlation()
{
  // We will first roughly test PuCorrMethod::ByPu239Only to data given in paper.
  //
  // Fig 3 in Swinhoe 2010 gives Pu239 content vs Pu242 content; I manually
  //  extracted the following values from the fit line in the PDF.
  const vector<pair<double,double>> swinhoe_approx_fig_3_data = {
    {0.55496, 0.06998},
    {0.55894, 0.06821},
    {0.56438, 0.06619},
    {0.57073, 0.06399},
    {0.57829, 0.06154},
    {0.58453, 0.05952},
    {0.59068, 0.05756},
    {0.59471, 0.05634},
    {0.59844, 0.05524},
    {0.60146, 0.05444},
    {0.60579, 0.05322},
    {0.61526, 0.05077},
    {0.62101, 0.04906},
    {0.62605, 0.04802},
    {0.63088, 0.04691},
    {0.63542, 0.04594},
    {0.6402, 0.0449}
  };//swinhoe_approx_fig_3_data
  
  
  for( const auto x_y : swinhoe_approx_fig_3_data )
  {
    const double x = x_y.first;
    const double y = x_y.second;
    const double gamma_spec_pu239 = x / (1.0 - y);
    const double gamma_spec_pu_other = (1.0 - x - y)/(1.0 - y);
    
    // gamma_spec_pu239 plus gamma_spec_pu_other should sum to 1.0
    assert( fabs(1.0 - (gamma_spec_pu239 + gamma_spec_pu_other)) < 0.001 );
    
    Pu242ByCorrelationInput input;
    input.pu238_rel_mass = gamma_spec_pu_other;
    input.pu239_rel_mass = gamma_spec_pu239;
    // Pu240, and Pu241/Am241 are irrelevant, all that
    
    Pu242ByCorrelationOutput output = correct_pu_mass_fractions_for_pu242( input, PuCorrMethod::ByPu239Only );
    
    //cout << "For Swinhoe [" << x << ", " << y << "]: Pu239: " << output.pu239_mass_frac
    //     << ", Pu242: " << output.pu242_mass_frac << " +- " << 100.0*output.pu242_uncert << "%\n";
    
    assert( fabs(output.pu239_mass_frac - x) < 0.005 );
    assert( fabs(output.pu242_mass_frac - y) < 0.0005 );
  }//for( const auto x_y : swinhoe_approx_fig_3_data )
  
  
  // For PuCorrMethod::Bignan95_BWR and PuCorrMethod::Bignan95_PWR, we dont have nearly as good
  //  of comparison data
  Pu242ByCorrelationInput input;
  input.pu238_rel_mass = 0.0120424;
  input.pu239_rel_mass = 0.6649628;
  input.pu240_rel_mass = 0.2327493;
  input.pu241_rel_mass = 0.0501864;
  input.pu241_rel_mass = 0.0361259;
  Pu242ByCorrelationOutput output = correct_pu_mass_fractions_for_pu242( input, PuCorrMethod::Bignan95_BWR );
  cout << "For Bignan95_BWR: Pu239: " << output.pu239_mass_frac
       << ", Pu242: " << output.pu242_mass_frac << " +- " << 100.0*output.pu242_uncert << "%\n";
  
  output = correct_pu_mass_fractions_for_pu242( input, PuCorrMethod::Bignan95_PWR );
  cout << "For Bignan95_PWR: Pu239: " << output.pu239_mass_frac
       << ", Pu242: " << output.pu242_mass_frac << " +- " << 100.0*output.pu242_uncert << "%\n";
}//void test_pu242_by_correlation()


void PhysicalModelShieldInput::check_valid() const
{
  const PhysicalModelShieldInput &input = *this;
    
  double ad_gcm2 = input.areal_density / PhysicalUnits::g_per_cm2;
  const double ad_ll_gcm2 = input.lower_fit_areal_density / PhysicalUnits::g_per_cm2;
  double ad_ul_gcm2 = input.upper_fit_areal_density / PhysicalUnits::g_per_cm2;
  if( ad_ul_gcm2 == 0.0 )
    ad_ul_gcm2 = sm_upper_allowed_areal_density_in_g_per_cm2;
  
  if( input.fit_areal_density )
  {
    if( ad_gcm2 == 0.0 )
      ad_gcm2 = 0.5*(ad_ll_gcm2 + ad_ul_gcm2);
    
    if( ad_ul_gcm2 <= ad_ll_gcm2 )
      throw runtime_error( "PhysicalModelShieldInput: upper AD must be greater than lower AD." );
    
    if( (ad_gcm2 < ad_ll_gcm2) || (ad_gcm2 > ad_ul_gcm2) )
      throw runtime_error( "PhysicalModelShieldInput: starting AD must be within limits." );
  }//if( input.fit_areal_density )
  
  if( (ad_gcm2 < 0.0) || (ad_gcm2 > sm_upper_allowed_areal_density_in_g_per_cm2) )
    throw runtime_error( "PhysicalModelShieldInput: AD must be in range [0,500] g/cm2." );
  
  if( input.material )
  {
    if( input.atomic_number != 0.0 )
      throw runtime_error( "PhysicalModelShieldInput: when material is specified"
                          " the atomic number must be set to zero");
    if( input.fit_atomic_number )
      throw runtime_error( "PhysicalModelShieldInput: when material is specified, must set fit"
                          " atomic number to false." );
  }else
  {
    if( input.fit_atomic_number )
    {
      if( (input.atomic_number != 0.0)
         && ((input.atomic_number < 1.0) || (input.atomic_number > 98)) )
        throw runtime_error( "PhysicalModelShieldInput: atomic number must be in [1,98], or 0 when"
                            " fitting it." );
      
      if( (input.lower_fit_atomic_number != 0.0)
         && ((input.lower_fit_atomic_number < 1.0) || (input.lower_fit_atomic_number > 98)) )
        throw runtime_error( "PhysicalModelShieldInput: lower atomic number must be in [1,98]." );
      
      if( (input.upper_fit_atomic_number != 0.0)
         && ((input.upper_fit_atomic_number < 1.0) || (input.upper_fit_atomic_number > 98)) )
        throw runtime_error( "PhysicalModelShieldInput: upper atomic number must be in [1,98]." );
      
      if( (input.upper_fit_atomic_number <= input.lower_fit_atomic_number)
         && (input.upper_fit_atomic_number != 0.0) )
        throw runtime_error( "PhysicalModelShieldInput: upper atomic number must be less than lower." );
    }else
    {
      if( (input.atomic_number < 1.0) || (input.atomic_number > 98) )
        throw runtime_error( "PhysicalModelShieldInput: invalid atomic number; must be in [1,98]" );
    }
  }//if( input.material ) / else.
    
};//PhysicalModelShieldInput::check_valid()
  


double eval_physical_model_eqn( const double energy,
                                 const shared_ptr<const Material> &self_atten,
                                 const vector<shared_ptr<const Material>> &external_attens,
                                 const DetectorPeakResponse &drf,
                                 const double * const paramaters,
                                 const size_t num_pars )
{
  const size_t num_expect_pars = 2 + 2*external_attens.size() + 2;
  assert( num_expect_pars == num_pars );
  if( num_expect_pars != num_pars )
    throw runtime_error( "eval_physical_model_eqn: invalid number of parameters." );
 
  const auto sanity_check_shield = [paramaters,num_pars]( const shared_ptr<const Material> &material,
                                                const size_t par_start ){
    double atomic_number = paramaters[par_start];
    double areal_density = paramaters[par_start + 1];
    
    assert( (areal_density >= -1.0E-6) && !IsInf(areal_density) );
    if( (areal_density <= -1.0E-6) || IsNan(areal_density) || IsInf(areal_density) )
      throw runtime_error( "eval_physical_model_eqn: areal density must be >= 0." );
    areal_density = std::max( 0.0, areal_density );

    if( material )
    {
      assert( atomic_number == 0.0f );
      if( atomic_number != 0.0f )
        throw runtime_error( "eval_physical_model_eqn: atomic number must be zero if material defined." );
    }else
    {
      if( IsNan(atomic_number) || IsInf(atomic_number) )
        throw runtime_error( "eval_physical_model_eqn: atomic number is inf or NaN." );

      if( (atomic_number < 0.9) || (atomic_number > 98.1) )
        throw runtime_error( "eval_physical_model_eqn: atomic number must in in range [1,98]" );

      atomic_number = std::max( 1.0, std::min( 98.0, atomic_number ) );
    }
  };//sanity_check_shield lamda
  
  sanity_check_shield( self_atten, 0 );
  for( size_t i = 0; i < external_attens.size(); ++i )
    sanity_check_shield( external_attens[i], 2 + i*2 );
  
  using GammaInteractionCalc::mass_attenuation_coef;
  using GammaInteractionCalc::transmition_length_coefficient;
    
    
  const float energyf = static_cast<float>( energy );
    
  double self_atten_an = paramaters[0];
  double self_atten_ad = paramaters[1] * PhysicalUnits::g / PhysicalUnits::cm2;
  
  assert( self_atten_ad >= -1.0E-3 );
  if( self_atten_ad < -1.0E-3 )
    throw runtime_error( "eval_physical_model_eqn: self atten areal-density may not be less than zero." );
  self_atten_ad = std::max( 0.0, self_atten_ad );

  assert( !IsNan(self_atten_an) && !IsInf(self_atten_an) && (self_atten_an > 0.9) && (self_atten_an < 98.1) );

  if( IsNan(self_atten_an) || IsInf(self_atten_an) || (self_atten_an < 0.9) && (self_atten_an > 98.1) )
    throw runtime_error( "eval_physical_model_eqn: self atten atomic number must be in range [1,98]." );
  self_atten_an = std::max( 1.0, std::min( 98.0, self_atten_an ) );

  const double self_atten_mu = self_atten
            ? transmition_length_coefficient( self_atten.get(), energyf ) / self_atten->density
            :  mass_attenuation_coef( self_atten_an, energyf );
    
  assert( self_atten_mu > 0.0 );
  
  const double sa_part = (self_atten_ad == 0.0)
                          ? 1.0
                          : (1 - exp(-self_atten_mu * self_atten_ad)) / (self_atten_mu * self_atten_ad);
    
  double ext_atten_part = 1.0;
  for( size_t i = 0; i < external_attens.size(); ++i )
  {
    const double atten_an = paramaters[2 + 2*i];
    const double atten_ad = paramaters[2 + 2*i + 1] * PhysicalUnits::g / PhysicalUnits::cm2;
      
    const shared_ptr<const Material> &mat = external_attens[i];
    const double self_atten_mu = mat
              ? transmition_length_coefficient( mat.get(), energyf ) / mat->density
              :  mass_attenuation_coef( atten_an, energyf );
      
    ext_atten_part *= exp( -self_atten_mu * atten_ad );
  }//for( size_t i = 0; i < external_attens.size(); ++i )
    
    
  const double det_part = drf.intrinsicEfficiency( energyf );
    
//#warning "Not including correction factor in Physical Model Eff"
//static int ntimeshere = 0;
//if( ntimeshere++ < 5 )
//  cerr << "eval_physical_model_eqn: corr_factor not included in Physical Model Eff" << endl;
//  double corr_factor = 1.0;
  const double b = paramaters[2 + 2*external_attens.size() + 0];
  const double c = paramaters[2 + 2*external_attens.size() + 1];
  const double corr_factor = std::pow(0.001*energy, b) * std::pow( c, 1000.0/energy );
    
  return sa_part * ext_atten_part * det_part * corr_factor;
}//eval_physical_model_eqn(...)
  
    
double eval_physical_model_eqn_uncertainty( const double energy,
                                 const std::shared_ptr<const Material> &self_atten,
                                 const std::vector<std::shared_ptr<const Material>> &external_attens,
                                 const DetectorPeakResponse &drf,
                                 const std::vector<std::vector<double>> &covariance )
{
#warning "eval_physical_model_eqn_uncertainty not implemented."
  cerr << "eval_physical_model_eqn_uncertainty not implemented." << endl;
  return 0.0;
}
  
std::function<double(double)> physical_model_eff_function( const shared_ptr<const Material> &self_atten,
                                  const vector<shared_ptr<const Material>> &external_attens,
                                  const DetectorPeakResponse &drf,
                                  const double * const paramaters,
                                  const size_t num_pars )
{
  const size_t num_expect_pars = 2 + 2*external_attens.size() + 2;
  assert( num_expect_pars == num_pars );
  if( num_expect_pars != num_pars )
    throw runtime_error( "physical_model_eff_function: invalid number of parameters." );
 
  const auto sanity_check_shield = [paramaters]( const shared_ptr<const Material> &material,
                                                const size_t par_start ){
    double atomic_number = paramaters[par_start];
    double areal_density = paramaters[par_start + 1];
    
    assert( (areal_density >= -1.0E-6) && !IsInf(areal_density) );
    if( (areal_density < -1.0E-6) || IsNan(areal_density) || IsInf(areal_density) )
      throw runtime_error( "physical_model_eff_function: areal density must be >= 0." );
    areal_density = std::max( 0.0, areal_density );
    
    assert( !IsNan(atomic_number) && !IsInf(atomic_number) && (atomic_number > 0.9) && (atomic_number < 98.1) );
    if( IsNan(atomic_number) || IsInf(atomic_number) || (atomic_number < 0.9) || (atomic_number > 98.1) )
      throw runtime_error( "physical_model_eff_function: atomic number must be in range [1,98]." );
    atomic_number = std::max( 1.0, std::min( 98.0, atomic_number ) );

    if( material )
    {
      assert( atomic_number == 0.0f );
      if( atomic_number != 0.0f )
        throw runtime_error( "physical_model_eff_function: atomic number must be zero if material defined." );
    }else
    {
      if( IsNan(atomic_number) || IsInf(atomic_number) )
        throw runtime_error( "physical_model_eff_function: atomic number is inf or NaN." );
      if( !((atomic_number >= 1.0) && (atomic_number <= 98.0)) ) //should catch Inf an NaN
        throw runtime_error( "physical_model_eff_function: atomic number must in in range [1,98]" );
    }
  };//sanity_check_shield lamda
  
  sanity_check_shield( self_atten, 0 );
  for( size_t i = 0; i < external_attens.size(); ++i )
    sanity_check_shield( external_attens[i], 2 + i*2 );
  
  if( paramaters[0] < 0.0 )
    throw runtime_error( "physical_model_eff_function: self atten areal-density may not be zero." );
  
  const function<float( float )> drffcn = drf.intrinsicEfficiencyFcn();
  const vector<double> pars( paramaters, paramaters + num_pars );
  
  
  return [drffcn, pars, self_atten, external_attens]( double energy ) -> double {
    using GammaInteractionCalc::mass_attenuation_coef;
    using GammaInteractionCalc::transmition_length_coefficient;
    
    
    const float energyf = static_cast<float>( energy );
    
    const double self_atten_an = pars[0];
    const double self_atten_ad = pars[1] * PhysicalUnits::g / PhysicalUnits::cm2;
    
    const double self_atten_mu = self_atten
            ? transmition_length_coefficient( self_atten.get(), energyf ) / self_atten->density
            :  mass_attenuation_coef( self_atten_an, energyf );
    
    const double sa_part = (1 - exp(-self_atten_mu * self_atten_ad)) / (self_atten_mu * self_atten_ad);
    
    double ext_atten_part = 1.0;
    for( size_t i = 0; i < external_attens.size(); ++i )
    {
      const double atten_an = pars[2 + 2*i];
      const double atten_ad = pars[2 + 2*i + 1] * PhysicalUnits::g / PhysicalUnits::cm2;
      
      const shared_ptr<const Material> &mat = external_attens[i];
      const double self_atten_mu = mat
              ? transmition_length_coefficient( mat.get(), energyf ) / mat->density
              :  mass_attenuation_coef( atten_an, energyf );
      
      ext_atten_part *= exp( -self_atten_mu * atten_ad );
    }//for( size_t i = 0; i < external_attens.size(); ++i )
    
    
    const double det_part = drffcn( energyf );
    
    const double b = pars[2 + 2*external_attens.size() + 0];
    const double c = pars[2 + 2*external_attens.size() + 1];
    const double corr_factor = std::pow(0.001*energy, b) * std::pow( c, 0.001*energy );
    
    return sa_part * ext_atten_part * det_part * corr_factor;
  };
}//physical_model_eff_function(...)
  


string physical_model_rel_eff_eqn_text( const std::shared_ptr<const Material> &self_atten,
                                                const vector<shared_ptr<const Material>> &external_attens,
                                                const DetectorPeakResponse &drf,
                                                const double * const paramaters,
                                                const size_t num_pars,
                                                const bool html_format )
{
  string eqn;
  if( !paramaters || (num_pars == 2) )
    return eqn;
  
  const size_t sa_an_index = 0;
  const size_t sa_ad_index = 1;
  if( self_atten || (paramaters[sa_ad_index] >= 1) )
  {
    const double sa_an = paramaters[sa_an_index];
    const double sa_ad = paramaters[sa_ad_index];
    const string sa_ad_str = SpecUtils::printCompact(sa_ad, 3);
    const string mu_name = self_atten ? cleanup_mat_name(self_atten->name) : SpecUtils::printCompact(sa_an, 3);
    if( html_format )
    {
      eqn +=
      "<span style=\"display: inline-block; vertical-align: middle;\">\n"
        "<span style=\"display: block; text-align: center;\">\n"
          "(1 - exp(-" + sa_ad_str + "*μ<sub>" + mu_name + "</sub>))\n"
        "</span>\n"
        "<span style=\"display: block; border-top: 1px solid black; text-align: center;\">\n"
          "(" + sa_ad_str + "*μ<sub>" + mu_name + "</sub>)\n"
        "</span>\n"
      "</span>\n";
    }else
    {
      eqn += "(1 - exp(-" + sa_ad_str + "*μ_" + mu_name + "))/(" + sa_ad_str + "*μ_" + mu_name + ")";
    }//if( html_format ) / else
  }//if( using self attenuation )
  
  if( !external_attens.empty() )
  {
    eqn += " * [";
    
    for( size_t i = 0; i < external_attens.size(); ++i )
    {
      const size_t an_index = 2 + 2*i + 0;
      const size_t ad_index = 2 + 2*i + 1;
      if( ad_index >= num_pars )
        continue;

      eqn += (i ? " * " : "");
      const shared_ptr<const Material> &mat = external_attens[i];
      const string mu_name = mat ? cleanup_mat_name(mat->name) : SpecUtils::printCompact(paramaters[an_index], 3);
      const string ad_str = SpecUtils::printCompact(paramaters[ad_index], 3);
      if( html_format )
      { 
        eqn += "exp(-" + ad_str + "*μ<sub>" + mu_name + "</sub>)";
      }else
      {
        eqn += "exp(-" + ad_str + "*μ_" + mu_name + ")";
      }//if( html_format ) / else
    }//for( size_t i = 0; i < external_attens.size(); ++i )
    
    eqn += "]";
  }//if( !external_attens.empty() )
       
  const size_t b_index = 2 + 2*external_attens.size() + 0;
  const size_t c_index = 2 + 2*external_attens.size() + 1;
  if( c_index < num_pars )
  {
    if( html_format )
    {
      eqn += " * [Detector Efficiency] * [E<sup>" + SpecUtils::printCompact(paramaters[b_index], 3) + "</sup> * "
        + SpecUtils::printCompact(paramaters[c_index], 3) + "<sup>1/E</sup>]";
    }else
    {
      eqn += " * [Detector Efficiency] * [E^" + SpecUtils::printCompact(paramaters[b_index], 3) + " * "
         + SpecUtils::printCompact(paramaters[c_index], 3) + "^(1/E)]";
    }//if( html_format ) / else
  }//if( c_index < num_pars )
   
  return eqn;
}//string physical_model_rel_eff_eqn_text(...)


std::string physical_model_rel_eff_eqn_js_function( const std::shared_ptr<const Material> &self_atten,
                                                       const std::vector<std::shared_ptr<const Material>> &external_attens,
                                                       const DetectorPeakResponse &drf,
                                                       const double * const paramaters,
                                                       const size_t num_pars )
{
  cerr << "physical_model_rel_eff_eqn_js_function: TODO: implement better than just interpolating" << endl;
  string fcn = "function(x){\n"
  "  const points = [";
  for( int x = 20; x < 3000; )
  {
    double y = eval_physical_model_eqn( x, self_atten, external_attens, drf, paramaters, num_pars );
    if( (y < 0.0) || IsNan(y) || IsInf(y) )
      y = 0.0;
    
    fcn += (x == 20) ? "" : ",";
    fcn += "[" + std::to_string(x) + "," + SpecUtils::printCompact(y, 4) + "]";
    
    if( x < 130 )
      x += 1;
    else if( x < 300 )
      x += 5;
    else
      x += 15;
  }//for( int x = 20; x < 3000; )
  fcn += "];\n"
  "  if( x <= points[0][0] )\n"
  "    return points[0][1];\n"
  "  if( x >= points[points.length - 1][0] )\n"
  "    return points[points.length - 1][1];\n"
  "  for (let i = 0; i < points.length - 1; i++) {\n"
  "    const [x1, y1] = points[i];\n"
  "    const [x2, y2] = points[i + 1];\n"
  "    if( x >= x1 && x <= x2) {\n"
  "      const t = (x - x1) / (x2 - x1);\n"
  "      return y1 + t * (y2 - y1);\n"
  "    }\n"
  "  }\n"
  "console.assert(0,'Shouldnt get here in interpolating');\n"
  "return points[points.length - 1][1];"
  "}";
  
  return fcn;
}

  
}//namespace RelActCalc
