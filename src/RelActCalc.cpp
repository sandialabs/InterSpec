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
#include "InterSpec/XmlUtils.hpp"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"

#include "InterSpec/RelActCalc_imp.hpp"

using namespace std;

namespace 
{
  /** Make "Pu (plutonium)" into "Pu", for printing PhysicalModelShieldInput materials. */
  string cleanup_mat_name( const std::shared_ptr<const Material> &material )
  {
    assert( material );
    if( !material )
      return "null";
    
    string text = material->name;
    if( !material->description.empty() && (material->description.size() < material->name.size()) )
      text = material->description;
    
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
  return RelActCalc::eval_eqn_imp<double>( energy, eqn_form, coeffs, num_coefs );
}//eval_eqn(...)


double eval_eqn_uncertainty( const double energy, const RelEffEqnForm eqn_form,
                            const std::vector<double> &coefs,
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
        
        for( size_t j = 0; j < covariance.size(); ++j )
          uncert_sq += std::pow(log_energy,1.0*i) * covariance[i][j] * std::pow(log_energy,1.0*j);
      }//for( size_t i = 0; i < coefs.size(); ++i )
      
      assert( uncert_sq >= 0.0 );
      
      return sqrt(uncert_sq);
    }//case RelEffEqnForm::LnX:
      
      
    case RelEffEqnForm::LnY:
    {
      // y = exp( a + b*x + c/x + d/x^2 + e/x^3 + ... )
      const double eval_val = eval_eqn( energy, eqn_form, coefs );
      
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
        
        for( size_t j = 0; j < covariance.size(); ++j )
          log_uncert_sq += eval_term_log(i) * covariance[i][j] * eval_term_log(j);
      }//for( size_t i = 0; i < coefs.size(); ++i )
      
      assert( log_uncert_sq >= 0.0 );
      
      return eval_val * sqrt(log_uncert_sq);
    }//case RelEffEqnForm::LnY:
      
    case RelEffEqnForm::LnXLnY:
    {
      // y = exp(a  + b*(lnx) + c*(lnx)^2 + d*(lnx)^3 + ... )
      const double log_energy = std::log(energy);
      const double eval_val = eval_eqn( energy, eqn_form, coefs );
      
      double log_uncert_sq = 0.0;
      
      for( size_t i = 0; i < covariance.size(); ++i )
      {
        assert( covariance[i].size() == covariance.size() );
        if( covariance[i].size() != covariance.size() )  //JIC for release builds
          throw runtime_error( "eval_eqn_uncertainty: covariance not a square matrix." );
        
        for( size_t j = 0; j < covariance.size(); ++j )
        log_uncert_sq += std::pow(log_energy,1.0*i) * covariance[i][j] * std::pow(log_energy,1.0*j);
      }//for( size_t i = 0; i < coefs.size(); ++i )
      
      assert( log_uncert_sq >= 0.0 );
      
      return eval_val * sqrt(log_uncert_sq);
    }//case RelEffEqnForm::LnXLnY:
      
    case RelEffEqnForm::FramEmpirical:
    {
      // y = exp( a + b/x^2 + c*(lnx) + d*(lnx)^2 + e*(lnx)^3 )
      const double log_energy = std::log(energy);
      const double eval_val = eval_eqn( energy, eqn_form, coefs );
      
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
        
        for( size_t j = 0; j < covariance.size(); ++j )
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
      
      return eval_val * sqrt(log_uncert_sq);
    }//case RelEffEqnForm::FramEmpirical:
      
    case RelEffEqnForm::FramPhysicalModel:
      assert( 0 );
      throw runtime_error( "RelActCalc::eval_eqn_uncertainty() can not be called for FramPhysicalModel." );
  }//switch( eqn_form )
  
  
  assert( 0 );   //shouldnt ever get here
  return -999.9; //avoid compiler warning
}//double eval_eqn_uncertainty(...)


std::vector<std::tuple<const SandiaDecay::Nuclide *,double,double>> back_decay_relative_activities(
                                            const double back_decay_time,
                                            std::vector<std::tuple<const SandiaDecay::Nuclide *,double>> &nuclide_rel_acts )
{
  assert( back_decay_time >= 0.0 );
  vector<tuple<const SandiaDecay::Nuclide *,double,double>> answer;

  double rel_mass_sum = 0.0;
  for( const tuple<const SandiaDecay::Nuclide *,double> &nuc_activity : nuclide_rel_acts )
  {
    const SandiaDecay::Nuclide * const nuc = get<0>(nuc_activity);
    const double final_act = get<1>(nuc_activity);
    assert( nuc );
    if( !nuc || (final_act <= 0.0) )
      continue;

    const double decrease_factor = std::exp( -back_decay_time * nuc->decayConstant() );
    const double initial_activity = final_act / decrease_factor;

    answer.emplace_back( nuc, initial_activity, 0.0 );

    const double initial_rel_mass = initial_activity / nuc->activityPerGram();
    rel_mass_sum += initial_rel_mass;
  }//for( const tuple<const SandiaDecay::Nuclide *,double> &nuc_activity : nuclide_rel_acts )


  for( size_t i = 0; i < answer.size(); ++i )
  {
    const SandiaDecay::Nuclide * const nuc = get<0>(answer[i]);
    const double initial_activity = get<1>(answer[i]);
    const double initial_rel_mass = initial_activity / nuc->activityPerGram();

    get<2>(answer[i]) = initial_rel_mass / rel_mass_sum;
  }

  return answer;
}//back_decay_relative_activities(...)


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

    vector<tuple<const SandiaDecay::Nuclide *,double>> input_nuclide_rel_acts;
    if( input.pu238_rel_mass > 0.0f )
      input_nuclide_rel_acts.emplace_back( pu238, input.pu238_rel_mass * pu238->activityPerGram() );
    if( input.pu239_rel_mass > 0.0f )
      input_nuclide_rel_acts.emplace_back( pu239, input.pu239_rel_mass * pu239->activityPerGram() );
    if( input.pu240_rel_mass > 0.0f )
      input_nuclide_rel_acts.emplace_back( pu240, input.pu240_rel_mass * pu240->activityPerGram() );
    if( input.pu240_rel_mass > 0.0f )
      input_nuclide_rel_acts.emplace_back( pu241, input.pu241_rel_mass * pu241->activityPerGram() );

    const vector<tuple<const SandiaDecay::Nuclide *,double,double>> time_zero_vals
                                  = back_decay_relative_activities( input.pu_age, input_nuclide_rel_acts );
    for( const tuple<const SandiaDecay::Nuclide *,double,double> &nuc_act_mass : time_zero_vals )
    {
      const SandiaDecay::Nuclide * const nuc = get<0>(nuc_act_mass);
      const double rel_act = get<1>(nuc_act_mass);
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

  float sum_input_mass = 0.0f;
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
  answer.pu242_mass_frac = pu242_mass_frac;

  if( input.pu_age > 0.0 )
  {
    const double pre_decay_sum_mass_frac = input.pu238_rel_mass + input.pu239_rel_mass + input.pu240_rel_mass
                                                + input.pu241_rel_mass + answer.pu242_mass_frac;

    answer.pu238_mass_frac *= std::exp( -input.pu_age * pu238->decayConstant() );
    answer.pu239_mass_frac *= std::exp( -input.pu_age * pu239->decayConstant() );
    answer.pu240_mass_frac *= std::exp( -input.pu_age * pu240->decayConstant() );
    answer.pu241_mass_frac *= std::exp( -input.pu_age * pu241->decayConstant() );
    answer.pu242_mass_frac *= std::exp( -input.pu_age * pu242->decayConstant() );


    const double post_decay_sum_mass_frac = input.pu238_rel_mass + input.pu239_rel_mass + input.pu240_rel_mass
                                              + input.pu241_rel_mass + answer.pu242_mass_frac;

    const double overall_decay = pre_decay_sum_mass_frac / post_decay_sum_mass_frac;
    answer.pu238_mass_frac *= overall_decay;
    answer.pu239_mass_frac *= overall_decay;
    answer.pu240_mass_frac *= overall_decay;
    answer.pu241_mass_frac *= overall_decay;
    answer.pu242_mass_frac *= overall_decay;
  }//if( input.pu_age > 0.0 )

  
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






double mass_ratio_to_act_ratio( const SandiaDecay::Nuclide * const numerator_nuclide, 
                                const SandiaDecay::Nuclide * const denominator_nuclide, 
                                const double mass_ratio )
{
  assert( numerator_nuclide );
  assert( denominator_nuclide );
  assert( mass_ratio > 0.0 );
/*
Nuc1: 3 bq/g
Nuc2: 2 bq/g

constraining Nuc1, controlling Nuc2

want mass ratio 1:
  mass_ratio_to_act_ratio( Nuc1, Nuc2, 1 ) --> (3/2)*1 = 3/2 (3 bq Nuc1, 2 bq Nuc2)

want mass ratio 0.5:
  mass_ratio_to_act_ratio( Nuc1, Nuc2, 0.5 ) --> (3/2)*0.5 = 3/4 (3 bq Nuc1, 4 bq Nuc2 --> 1g Nuc1, 2g Nuc2)
*/

  const double num_act_per_g = numerator_nuclide->activityPerGram();
  const double denom_act_per_g = denominator_nuclide->activityPerGram();
  const double act_ratio = mass_ratio * num_act_per_g / denom_act_per_g;

  return act_ratio;
}

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
  
  
#if( PERFORM_DEVELOPER_CHECKS )
void PhysicalModelShieldInput::equalEnough( const PhysicalModelShieldInput &lhs, const PhysicalModelShieldInput &rhs )
{
  if( fabs(lhs.atomic_number - rhs.atomic_number) > 1.0E-6 )
    throw runtime_error( "PhysicalModelShieldInput::equalEnough: atomic number mismatch" );
  
  if( !lhs.material != !rhs.material )
    throw runtime_error( "PhysicalModelShieldInput::equalEnough: material being present mismatch" );
  
  if( lhs.material && (lhs.material->name != rhs.material->name) )
    throw runtime_error( "PhysicalModelShieldInput::equalEnough: material name mismatch" );
  
  if( fabs(lhs.areal_density - rhs.areal_density) > 1.0E-3 )
    throw runtime_error( "PhysicalModelShieldInput::equalEnough: areal density mismatch" );
  
  if( lhs.fit_atomic_number != rhs.fit_atomic_number )
    throw runtime_error( "PhysicalModelShieldInput::equalEnough: fit atomic number mismatch" );
  
  if( fabs(lhs.lower_fit_atomic_number - rhs.lower_fit_atomic_number) > 1.0E-4 )
    throw runtime_error( "PhysicalModelShieldInput::equalEnough: lower fit atomic number mismatch" );
  
  if( fabs(lhs.upper_fit_atomic_number - rhs.upper_fit_atomic_number) > 1.0E-4 )
    throw runtime_error( "PhysicalModelShieldInput::equalEnough: upper fit atomic number mismatch" );
  
  if( lhs.fit_areal_density != rhs.fit_areal_density )
    throw runtime_error( "PhysicalModelShieldInput::equalEnough: fit areal density mismatch" );
  
  if( fabs(lhs.lower_fit_areal_density - rhs.lower_fit_areal_density) > 1.0E-6 )
    throw runtime_error( "PhysicalModelShieldInput::equalEnough: lower fit areal density mismatch" );
  
  if( fabs(lhs.upper_fit_areal_density - rhs.upper_fit_areal_density) > 1.0E-3 )
    throw runtime_error( "PhysicalModelShieldInput::equalEnough: upper fit areal density mismatch" );
}//void PhysicalModelShieldInput::equalEnough(...)
#endif
  
rapidxml::xml_node<char> *PhysicalModelShieldInput::toXml( ::rapidxml::xml_node<char> *parent ) const
{
  using namespace rapidxml;
  
  assert( parent );
  if( !parent || !parent->document() )
    throw runtime_error( "PhysicalModelShieldInput::toXml: invalid parent." );
  
  xml_document<char> *doc = parent->document();
  xml_node<char> *base_node = doc->allocate_node( node_element, "PhysicalModelShield" );
  parent->append_node( base_node );
  
  XmlUtils::append_version_attrib( base_node, PhysicalModelShieldInput::sm_xmlSerializationVersion );
  
  XmlUtils::append_float_node( base_node, "AtomicNumber", atomic_number );
  if( material )
    XmlUtils::append_string_node( base_node, "Material", material->name );
  XmlUtils::append_float_node( base_node, "ArealDensity", areal_density / PhysicalUnits::g_per_cm2 );
  XmlUtils::append_bool_node( base_node, "FitAtomicNumber", fit_atomic_number );
  XmlUtils::append_float_node( base_node, "LowerFitAtomicNumber", lower_fit_atomic_number );
  XmlUtils::append_float_node( base_node, "UpperFitAtomicNumber", upper_fit_atomic_number );
  XmlUtils::append_bool_node( base_node, "FitArealDensity", fit_areal_density );
  XmlUtils::append_float_node( base_node, "LowerFitArealDensity", lower_fit_areal_density / PhysicalUnits::g_per_cm2 );
  XmlUtils::append_float_node( base_node, "UpperFitArealDensity", upper_fit_areal_density / PhysicalUnits::g_per_cm2 );

  return base_node;
}//rapidxml::xml_node<char> *toXml( ::rapidxml::xml_node<char> *parent ) const
  
  
void PhysicalModelShieldInput::fromXml( const ::rapidxml::xml_node<char> *parent, MaterialDB *materialDB )
{
  try
  {
    if( !parent )
      throw runtime_error( "invalid input" );
    
    if( !rapidxml::internal::compare( parent->name(), parent->name_size(), "PhysicalModelShield", 19, false ) )
      throw std::logic_error( "invalid input node name" );
    
    // A reminder double check these logics when changing RoiRange::sm_xmlSerializationVersion
    static_assert( PhysicalModelShieldInput::sm_xmlSerializationVersion == 0,
                  "needs to be updated for new serialization version." );
    
    XmlUtils::check_xml_version( parent, PhysicalModelShieldInput::sm_xmlSerializationVersion );
    
    const rapidxml::xml_node<char> *material_node = XML_FIRST_NODE( parent, "Material" );
    
    
    areal_density = 0.0;
    lower_fit_areal_density = 0.0;
    upper_fit_areal_density = 0.0;
    fit_areal_density = XmlUtils::get_bool_node_value( parent, "FitArealDensity" );
    
    if( fit_areal_density )
    {
      // We require AD limits if we are fitting AD
      lower_fit_areal_density = XmlUtils::get_float_node_value( parent, "LowerFitArealDensity" );
      upper_fit_areal_density = XmlUtils::get_float_node_value( parent, "UpperFitArealDensity" );
      
      // Wont require AD if we are fitting it
      try
      {
        areal_density = XmlUtils::get_float_node_value( parent, "ArealDensity" );
      }catch( std::exception & )
      {
      }
    }else
    {
      // We wont really require AD limits if we are not fitting AD
      try
      {
        lower_fit_areal_density = XmlUtils::get_float_node_value( parent, "LowerFitArealDensity" );
      }catch( std::exception & )
      {
      }
      
      try
      {
        upper_fit_areal_density = XmlUtils::get_float_node_value( parent, "UpperFitArealDensity" );
      }catch( std::exception & )
      {
      }
      
      // Require AD if we are fitting it
      areal_density = XmlUtils::get_float_node_value( parent, "ArealDensity" );
    }//if( fit_areal_density ) / else
    
    areal_density *= PhysicalUnits::g_per_cm2;
    lower_fit_areal_density *= PhysicalUnits::g_per_cm2;
    upper_fit_areal_density *= PhysicalUnits::g_per_cm2;
    
    
    fit_atomic_number = false;
    lower_fit_atomic_number = 1.0;
    upper_fit_atomic_number = 98.0;
    
    if( material_node )
    {
      if( !materialDB )
        throw runtime_error( "PhysicalModelShieldInput::fromXml: need MaterialDB" );
      
      const string name = SpecUtils::xml_value_str(material_node);
      const Material * mat = materialDB->material( name );
      if( !mat )
      {
        try
        {
          const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
          mat = materialDB->parseChemicalFormula( name, db );
        }catch( std::exception & )
        {
        }
      }//if( !mat )
      
      if( !mat )
        throw runtime_error( "Invalid material name '" + name + "'" );
      
      material = make_shared<Material>( *mat );
      
      try
      {
        fit_atomic_number = XmlUtils::get_bool_node_value( parent, "FitAtomicNumber" );
      }catch( std::exception & )
      {
      }
      
      try
      {
        atomic_number = XmlUtils::get_float_node_value( parent, "AtomicNumber" );
        lower_fit_atomic_number = XmlUtils::get_float_node_value( parent, "LowerFitAtomicNumber" );
        upper_fit_atomic_number = XmlUtils::get_float_node_value( parent, "UpperFitAtomicNumber" );
      }catch( std::exception & )
      {
      }
    }else
    {
      fit_atomic_number = XmlUtils::get_bool_node_value( parent, "FitAtomicNumber" );
      atomic_number = XmlUtils::get_float_node_value( parent, "AtomicNumber" );
      lower_fit_atomic_number = XmlUtils::get_float_node_value( parent, "LowerFitAtomicNumber" );
      upper_fit_atomic_number = XmlUtils::get_float_node_value( parent, "UpperFitAtomicNumber" );
      
      double lower_an = 1.0, upper_an = 98.0;
      if( fit_atomic_number && (lower_fit_atomic_number == upper_fit_atomic_number) )
      {
        lower_fit_atomic_number = 1.0;
        upper_fit_atomic_number = 98.0;
      }
      
      if( fit_atomic_number )
      {
        if( (lower_fit_atomic_number < 1.0) || (lower_fit_atomic_number > 98.0) )
          throw runtime_error( "invalid lower atomic number fit range." );
        
        if( (upper_fit_atomic_number < 1.0) || (upper_fit_atomic_number > 98.0) )
          throw runtime_error( "invalid lower atomic number fit range." );
        
        if( lower_fit_atomic_number >= upper_fit_atomic_number )
          throw runtime_error( "lower atomic number is larger than upper atomic number." );
      }//if( fit_atomic_number )
      
      if( (atomic_number < 1.0) || (atomic_number > 98.0) )
      {
        assert( atomic_number == 0.0 );
        if( fit_atomic_number )
          atomic_number = 0.5*(lower_fit_atomic_number + upper_fit_atomic_number);
        else
          atomic_number = 26.0;
      }//
    }// if( material_node ) / else
  }catch( std::exception &e )
  {
    throw runtime_error( "PhysicalModelShieldInput::fromXml(): " + string(e.what()) );
  }
}//void fromXml( const ::rapidxml::xml_node<char> *parent )


double eval_physical_model_eqn( const double energy,
                               const std::optional<PhysModelShield<double>> &self_atten,
                               const std::vector<PhysModelShield<double>> &external_attens,
                               const DetectorPeakResponse * const drf,
                               std::optional<double> hoerl_b,
                               std::optional<double> hoerl_c )
{
  return eval_physical_model_eqn_imp<double>( energy, self_atten, external_attens, drf, hoerl_b, hoerl_c );
}//eval_physical_model_eqn(...)
  
    
double eval_physical_model_eqn_uncertainty( const double energy,
                               const std::optional<PhysModelShield<double>> &self_atten,
                               const std::vector<PhysModelShield<double>> &external_attens,
                               const DetectorPeakResponse * const drf,
                               std::optional<double> hoerl_b,
                               std::optional<double> hoerl_c,
                               const std::vector<std::vector<double>> &covariance )
{
#pragma message( "eval_physical_model_eqn_uncertainty not implemented." )
static int ntimeshere = 0;
if( ntimeshere++ < 5 )  
  cerr << "eval_physical_model_eqn_uncertainty not implemented. " << endl;
return 0.0;
}
  
std::function<double(double)> physical_model_eff_function( const std::optional<PhysModelShield<double>> &self_atten,
                                                          const std::vector<PhysModelShield<double>> &external_attens,
                                                          const std::shared_ptr<const DetectorPeakResponse> &drf,
                                                          std::optional<double> hoerl_b,
                                                          std::optional<double> hoerl_c )
{
  const auto sanity_check_shield = []( const PhysModelShield<double> &shield ){
    double atomic_number = shield.atomic_number;
    double areal_density = shield.areal_density;
    
    assert( (areal_density >= -1.0E-6) && !IsInf(areal_density) );
    if( (areal_density <= -1.0E-6) || IsNan(areal_density) || IsInf(areal_density) )
      throw runtime_error( "physical_model_eff_function: areal density must be >= 0." );
    areal_density = std::max( 0.0, areal_density );

    if( shield.material )
    {
      assert( atomic_number == 0.0f );
      if( atomic_number != 0.0f )
        throw runtime_error( "physical_model_eff_function: atomic number must be zero if material defined." );
    }else if( atomic_number == 0.0 )
    {
      assert( areal_density == 0.0 );
      if( areal_density != 0.0 )
        throw runtime_error( "physical_model_eff_function: areal density must be zero if atomic number is zero." );
    }else
    {
      assert( !IsNan(atomic_number) && !IsInf(atomic_number) && (atomic_number > 0.9) && (atomic_number < 98.1) );
      
      if( IsNan(atomic_number) || IsInf(atomic_number) )
        throw runtime_error( "physical_model_eff_function: atomic number is inf or NaN." );

      if( (atomic_number < 0.9) || (atomic_number > 98.1) )
        throw runtime_error( "physical_model_eff_function: atomic number must in in range [1,98]" );
    }
  };//sanity_check_shield lamda
  
  if( self_atten.has_value() )
    sanity_check_shield( *self_atten );
  for( size_t i = 0; i < external_attens.size(); ++i )
    sanity_check_shield( external_attens[i] );
  
  assert( hoerl_b.has_value() == hoerl_c.has_value() );
  if( hoerl_b.has_value() != hoerl_c.has_value() )
    throw std::logic_error( "hoerl_b.has_value() != hoerl_c.has_value()" );
  
  //const function<float(float)> drffcn = drf ? drf->intrinsicEfficiencyFcn() : []( float ){ return 1.0f; };
  
  return [drf, self_atten, external_attens, hoerl_b, hoerl_c]( double energy ) -> double {
    return eval_physical_model_eqn_imp<double>( energy, self_atten, external_attens, drf.get(), hoerl_b, hoerl_c );
  };
}//physical_model_eff_function(...)
  


string physical_model_rel_eff_eqn_text( const std::optional<PhysModelShield<double>> &self_atten,
                                       const std::vector<PhysModelShield<double>> &external_attens,
                                       const std::shared_ptr<const DetectorPeakResponse> &drf,
                                       std::optional<double> hoerl_b,
                                       std::optional<double> hoerl_c,
                                                const bool html_format )
{
  string eqn;
  
  if( self_atten.has_value() && (self_atten->material || ((self_atten->atomic_number > 0.9) && (self_atten->atomic_number < 98.1))) )
  {
    const double sa_an = self_atten->atomic_number;
    const double sa_ad = self_atten->areal_density / PhysicalUnits::g_per_cm2;
    const string sa_ad_str = SpecUtils::printCompact(sa_ad, 3);
    const string mu_name = self_atten->material ? cleanup_mat_name(self_atten->material) : SpecUtils::printCompact(sa_an, 3);
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
      const PhysModelShield<double> &atten = external_attens[i];
      if( !atten.material && ((atten.atomic_number < 0.9) || (atten.atomic_number > 98.1)) )
        continue;
        
      eqn += (i ? " * " : "");
      const shared_ptr<const Material> &mat = atten.material;
      const string mu_name = mat ? cleanup_mat_name(mat) : SpecUtils::printCompact(atten.atomic_number, 3);
      const string ad_str = SpecUtils::printCompact(atten.areal_density/PhysicalUnits::g_per_cm2, 3);
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

  eqn += (eqn.empty() ? "" : " * ");
  eqn += "[Det. Eff.]";

  
  if( hoerl_b.has_value() && hoerl_c.has_value() && ((hoerl_b.value() != 0.0) || (hoerl_c.value() != 1.0)) )
  {
    if( html_format )
    {
      eqn += " * [E<sup>" + SpecUtils::printCompact(hoerl_b.value(), 3) + "</sup> * "
        + SpecUtils::printCompact(hoerl_c.value(), 3) + "<sup>1/E</sup>]";
    }else
    {
      eqn += " * [E^" + SpecUtils::printCompact(hoerl_b.value(), 3) + " * "
         + SpecUtils::printCompact(hoerl_c.value(), 3) + "^(1/E)]";
    }//if( html_format ) / else
  }//if( c_index < num_pars )
   
  return eqn;
}//string physical_model_rel_eff_eqn_text(...)


std::string physical_model_rel_eff_eqn_js_function( const std::optional<PhysModelShield<double>> &self_atten,
                                                   const std::vector<PhysModelShield<double>> &external_attens,
                                                   const DetectorPeakResponse * const drf,
                                                   std::optional<double> hoerl_b,
                                                   std::optional<double> hoerl_c )
{
  cerr << "physical_model_rel_eff_eqn_js_function: TODO: implement better than just interpolating" << endl;
  string fcn = "function(x){\n"
  "  const points = [";
  for( int x = 20; x < 3000; )
  {
    double y = eval_physical_model_eqn( x, self_atten, external_attens, drf, hoerl_b, hoerl_c );
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
