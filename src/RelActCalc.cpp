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
#include <cstring>
#include <stdexcept>

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/RelActCalc.h"
#include "InterSpec/PhysicalUnits.h"

using namespace std;


namespace RelActCalc
{

const char *to_str( const RelEffEqnForm form )
{
  switch( form )
  {
    case RelEffEqnForm::LnX:    return "LnX";
    case RelEffEqnForm::LnY:    return "LnY";
    case RelEffEqnForm::LnXLnY: return "LnXLnY";
    case RelEffEqnForm::FramEmpirical: return "FRAM Empirical";
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
    RelEffEqnForm::FramEmpirical
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
          rel_eff_eqn_str += PhysicalUnits::printCompact(val,nsigfig);
        }else
        {
          rel_eff_eqn_str += (val < 0.0 ? " - " : " + " )
                             +  PhysicalUnits::printCompact(fabs(val),nsigfig)
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
          rel_eff_eqn_str += PhysicalUnits::printCompact(val,nsigfig);
        }else
        {
          rel_eff_eqn_str += (val < 0.0 ? " - " : " + " )
                             + PhysicalUnits::printCompact(fabs(val),nsigfig);
          if( i == 1 )
            rel_eff_eqn_str += "*x";
          else if( i == 2 )
            rel_eff_eqn_str += "/x";
          else
            rel_eff_eqn_str += "/x^" +  to_string( i - 1 );
        }
      }//for( size_t i = 0; i < coefs.size(); ++i )
      rel_eff_eqn_str += ")";
      
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
                             + PhysicalUnits::printCompact(fabs(val),nsigfig)
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
          rel_eff_eqn_str += PhysicalUnits::printCompact(val,nsigfig);
        }else if( i == 1 )
        {
          rel_eff_eqn_str += (val < 0.0 ? " - " : " + " )
                             + PhysicalUnits::printCompact(fabs(val),nsigfig) + "/(x*x)";
        }else
        {
          rel_eff_eqn_str += (val < 0.0 ? " - " : " + " )
                             + PhysicalUnits::printCompact(fabs(val),nsigfig)
                             + "*ln(x)^" + to_string(i-1);
        }
      }//for( size_t i = 0; i < coefs.size(); ++i )
      
      rel_eff_eqn_str += " )";
      
      break;
    }//case RelActCalc::RelEffEqnForm::FramEmpirical:
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
  }//switch( eqn_form )
  
  
  assert( 0 );   //shouldnt ever get here
  return -999.9; //avoid compiler warning
}//double eval_eqn_uncertainty(...)
}//namespace RelActCalc
