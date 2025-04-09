//
//  TerminalModel.cpp
//  InterSpec
//
//  Created by Christian Kenneth Morte on 7/15/16.
//

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

// Need to make sure we have M_PI, etc
#define _USE_MATH_DEFINES

#include "InterSpec_config.h"

#include <regex>
#include <tuple>
#include <vector>
#include <memory>
#include <stdio.h>
#include <sstream>
#include <iostream>
#include <functional>

#include <boost/math/distributions/chi_squared.hpp>

#include "Wt/WModelIndex"

#include "mpParser.h"

#include "InterSpec/PeakEdit.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/InterSpec.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/TerminalModel.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"



// The regex in GCC 4.8.x does not have working regex
#if( defined(__GLIBCXX__) && (__cplusplus < 201402L) )
static_assert( defined(_GLIBCXX_REGEX_DFS_QUANTIFIERS_LIMIT) \
              || defined(_GLIBCXX_REGEX_STATE_LIMIT) \
              || (defined(_GLIBCXX_RELEASE) && _GLIBCXX_RELEASE > 4), "GCC 4.8 is not supported due to buggy regex implementation" );
#endif


/*
 Helper Namespace:
 Define any useful functions here for helping in creating general functions for your parser. For example, I used this namespace to create the math functions
 such as modulus, division, and calculating p-values. In addition, I use this namespace to define various regular expressions that are used throughout usage
 of the Terminal.
 */
namespace {
    /* Regex Expressions */
    //The variableAssignmentRegex and validVariableRegex are invalid for gcc 4.8.4
    const char * const ns_variableAssignmentRegexArg = "^(?:\\s)*([A-Za-z]\\w*)(?:\\s)*?=(?:\\s)*([^=]+)$";
    const char * const ns_validVariableRegexArg = "^([a-zA-Z]\\w*|_[a-zA-Z]\\w*)+$";
    const char * const ns_keywordRegexArg = "^(cos|sin|tan|acos|asin|atan|sinh|cosh|tanh|asinh|acosh|atanh|log2|log10|log|ln|exp|sqrt|sign|rint|abs|min|max|sum|avg|default|empty)$";

    /* Internal Parser helper methods */
    std::string convertDoubleTypeToString(const double value)
    {
        std::ostringstream os;
        os << value;
        return os.str();
    }
    
    std::string errorMessage(const mup::ParserError& e) {
        std::ostringstream message;
        message << "Error code " << e.GetCode() << ": " << e.GetMsg();
        return message.str();
    }
    
    // Math functions
    double modulus( const int arg1, const int arg2 ) {
        if (arg2 < 0)   return (int) arg1 % - (int) arg2;
        double result = (int) arg1 % (int) arg2;
        if ( result < 0 ) return result + arg2;
        return result;
    }
    double pValueFromChiSquare( const double chiSquareScore, const double df ) {
        try {
            boost::math::chi_squared distribution( df );
            return 1 - boost::math::cdf( distribution, chiSquareScore );
        } catch ( const std::exception & ) {
            throw mup::ParserError( "Domain error calculating the p-value." );
        }
    }
    double toAngle( const double radians ) { return radians * (180/M_PI); }
    double toRadians( const double angle ) { return angle * (M_PI/180); }
}

/*  Define non-built-in parser functions within a default namespace under class declarations.
 This namespace carries all the functions that are to be defined within the parser of the TerminalModel class.
 To define a function:
 1. Create a class that INHERITS FROM mup::ICallback.
 a) Initialize the ICallback member with the following parameters ( mup::cmFUNC, <function_name>, <number of arguments> ) and any other members.
 where <function_name> is a string for the callback the user inputs, and <number of arguments> is the number of arguments the function needs.
 b) If you want to have the function use/compute a value from the internal spectrum, then you MUST also include and initialize a TerminalModel member within the
 class.
 
 2. Overload the method:   virtual void Eval( mup::ptr_val_type& ret, const mup::ptr_val_type* argv, int a_iArgc )
 where "ret" is a reference to the return value, "argv" is a pointer array to the arguments, and "a_iArgc" is the number of arguments.
 a) To successfully have a parser function return the result, assign the "ret" argument to the final result of the function evaluation.
 b) If you want to have the function use/compute a value from the internal spectrum, then you MUST also define a function within the TerminalModel class that
 executes whatever result/operation you want to have, and call that within the Eval method. See LiveTime and RealTime classes for an example.
 
 3. (OPTIONAL) Overload the GetDesc() and Clone() functions.
 a)  const mup::char_type* GetDesc() const   --> Returns callback description of function.
 b)  mup::IToken* Clone() const { return new <FUNCTION_CLASS>(*this); }
 where <FUNCTION_CLASS> is the class being implemented for the function.
 *** If you would like to include the function callback description in the function drop-down list, you MUST define the GetDesc() method ***
 
 4. (OPTIONAL) Create a tag function inside the function class. This function should return a string that represents any search tags that would be used for
    searching this function inside the Command Helper tool.
 a) For consistency, you should define this function as "std::string tags() const".
 b) Each search tag is separated by a string (so if you want to have certain tags together as one word, you must combine using other characters).
 
 5. (OPTIONAL) Create a tooltip function inside the funciton class. This function should rturn a string that represents the text inside the tooltip popup
    when the user highlights a function in the command helper menu.
 a) For consistency, you should define this function as "std::string toolTip() const".
 
 6. Inside the constructor for TerminalModel, dynamically allocate an instance of your function class. (create a pointer to the class and assign
 with the 'new' operator). Pass this pointer into the TerminalModel's addFunction( functionPointer, tags, toolTip  ) method. If you had defined the tags
 and/or tooltips function within your class, then pass a call to that method into the 'tags' argument and/or 'toolTip of "addFunction", respectively.

 For more help in defining parser functions/operators for the TerminalModel class, please follow the link:
 **************************************************************************************
 http://beltoforion.de/article.php?a=muparserx&hl=en&p=extending&s=idPageTop#idPageTop
 **************************************************************************************
 */

namespace {
    // Operators: both binary and unary operators. These act as separators between operands to execute evaluations.
    class Modulus : public mup::IOprtBin
    {
    public:
        Modulus(const char* sign) :IOprtBin(sign, (int)mup::prMUL_DIV, mup::oaNONE) {};
        virtual void Eval( mup::ptr_val_type& ret, const mup::ptr_val_type* argv, int a_iArgc ) {
            const mup::IValue *arg1 = argv[0].Get();
            const mup::IValue *arg2 = argv[1].Get();
            if ( !arg1->IsNonComplexScalar() )      throw mup::ParserError( "Modulus operator must precede a value of type scalar." );
            else if ( !arg2->IsNonComplexScalar() ) throw mup::ParserError( "Modulus operator must succede a value of type scalar." );
            *ret = modulus( arg1->GetInteger(), arg2->GetInteger() );
        }
        const mup::char_type* GetDesc() const { return "Modulus operator."; }
        mup::IToken* Clone() const { return new Modulus(*this); }
    };
    class Division : public mup::IOprtBin
    {
    public:
        Division(const char* sign) :IOprtBin(sign, (int)mup::prMUL_DIV, mup::oaNONE) {};
        virtual void Eval( mup::ptr_val_type& ret, const mup::ptr_val_type* argv, int a_iArgc ) {
            const mup::IValue *arg1 = argv[0].Get();
            const mup::IValue *arg2 = argv[1].Get();
            if ( !arg1->IsNonComplexScalar() )      throw mup::ParserError( "Integer division operator must precede a value of type scalar." );
            else if ( !arg2->IsNonComplexScalar() ) throw mup::ParserError( "Integer division operator must succede a value of type scalar." );
            
            *ret = div( arg1->GetInteger(),  arg2->GetInteger() ).quot;
        }
        const mup::char_type* GetDesc() const { return "Integer division operator."; }
        mup::IToken* Clone() const { return new Division(*this); }
    };
    
    // Trigonometric Functions
    class ToAngle : public mup::ICallback
    {
    public:
        ToAngle() :ICallback(mup::cmFUNC, "toAngle", 1) {};
        
        virtual void Eval( mup::ptr_val_type& ret, const mup::ptr_val_type* argv, int a_iArgc ) {
            if ( !argv[0]->IsNonComplexScalar() )
                throw mup::ParserError( "Argument for function 'toAngle' must be of type 'scalar value'." );
            *ret = toAngle( argv[0]->GetFloat() );
        }
        const mup::char_type* GetDesc() const { return "toAngle( radians )"; }
        std::string tags() const { return "trig trigonometry trigonometric sin cos sin cosine angle radians to_angle to angle convert conversion pi 180 degrees"; }
        std::string toolTip() const { return "Returns the <b>angle</b> value for a given <i>radian</i>."; }
        mup::IToken* Clone() const { return new ToAngle(*this); }
    };
    class ToRadians : public mup::ICallback
    {
    public:
        ToRadians() :ICallback(mup::cmFUNC, "toRadians", 1) {};
        
        virtual void Eval( mup::ptr_val_type& ret, const mup::ptr_val_type* argv, int a_iArgc ) {
            if ( !argv[0]->IsNonComplexScalar() )
                throw mup::ParserError( "Argument for function 'toRadians' must be of type 'scalar value'." );
            *ret = toRadians( argv[0]->GetFloat() );
        }
        const mup::char_type* GetDesc() const { return "toRadians( angle )"; }
        std::string tags() const { return "trig trigonometry trigonometric sin cos sin cosine angle radians to_angle to angle convert conversion pi 180 degrees"; }
        std::string toolTip() const { return "Returns the <b>radian</b> value for a given <i>angle</i>."; }
        mup::IToken* Clone() const { return new ToRadians(*this); }
    };
    
    
    
    // Nuclide functions
    class NuclideIntensity : public mup::ICallback
    {
    public:
        NuclideIntensity(TerminalModel *model, bool withArg=false)
        :ICallback(mup::cmFUNC, withArg ? "nuclideIntensityForParticle" : "nuclideIntensity", withArg ? 2 : 1), tm(model), withParticleArgument(withArg) {};
        
        virtual void Eval( mup::ptr_val_type& ret, const mup::ptr_val_type* argv, int a_iArgc ) {
            if ( withParticleArgument ) {
                if ( !argv[0]->IsString() ) throw mup::ParserError( "Argument for particle in function 'nuclideIntensityForParticle'" + std::string("' must be of type 'string'.") );
                else if ( !argv[1]->IsNonComplexScalar() )  throw mup::ParserError( "Argument for energy in function 'nuclideIntensityForParticle'" + std::string("' must be of type 'scalar value'.") );
                *ret = tm->nuclideIntensityForParticle( argv[0]->GetString(), argv[1]->GetFloat() );
                
            } else {
                if ( !argv[0]->IsNonComplexScalar() ) throw mup::ParserError( "Argument for 'nuclideIntensity'" + std::string("' must be of type 'scalar value'.") );
                *ret = tm->nuclideIntensity(  argv[0]->GetFloat() );
            }
        }
        const mup::char_type* GetDesc() const { return withParticleArgument ? "nuclideIntensityForParticle( particle, energy )" : "nuclideIntensity( energy )"; }
        std::string tags() const { return "br B.R. energy nuclides age particles highest"; }
        std::string toolTip() const { return withParticleArgument ?
            "Gets nuclide’s <b>intensity value</b> at specified <i>energy</i> and <i>particle</i>. Returns <b><font color='red'>error message</font></b> if both of these parameters could not be matched." :
            "Gets the nuclide’s <b>intensity value</b> at a specified <i>energy</i>. Returns <b><font color='red'>error message</font></b> if no energy value could be matched with user’s specification or no nuclide is shown."; }
        mup::IToken* Clone() const { return new NuclideIntensity(*this); }
    private: TerminalModel *tm; bool withParticleArgument;
    };
    
    
    // Live/Real time functions
    class LiveTime : public mup::ICallback
    {
    public:
        LiveTime(TerminalModel *model, bool withArg=false) :ICallback(mup::cmFUNC, withArg ? "liveTimeOf" : "liveTime", withArg ? 1 : 0), tm(model), withArgument(withArg) {};
        
        virtual void Eval( mup::ptr_val_type& ret, const mup::ptr_val_type* argv, int a_iArgc ) {
            if ( withArgument && !argv[0]->IsString() ) throw mup::ParserError( "Argument '" + argv[0]->GetString() + "' is not of type string/id." );
            if ( withArgument ) *ret = tm->liveTime( std::string( argv[0]->GetString() ) );
            else *ret = tm->liveTimeWithoutArgument();
        }
        const mup::char_type* GetDesc() const { return withArgument ? "liveTimeOf( spectrum )" : "liveTime()"; }
        std::string tags() const { return "live_time live time spectra foreground fg background bg secondary sfg"; }
        std::string toolTip() const { return withArgument ?
            "Returns the <b>live time</b> (in seconds) of <i>spectrum</i>, or <b>0</b> if not known. Returns <b><font color='red'>error message</font></b> if <i>spectrum</i> was not detected." :
            "Automatically detects single spectrum and returns its <b>live time value</b>. Returns <b><font color='red'>error message</font></b> if multiple or no spectra detected."; }
        mup::IToken* Clone() const { return new LiveTime(*this); }
    private: TerminalModel *tm; bool withArgument;
    };
    
    class RealTime : public mup::ICallback
    {
    public:
        RealTime(TerminalModel *model, bool withArg=false) :ICallback(mup::cmFUNC, withArg ? "realTimeOf" : "realTime", withArg ? 1 : 0), tm(model), withArgument(withArg) {};
        
        virtual void Eval( mup::ptr_val_type& ret, const mup::ptr_val_type* argv, int a_iArgc ) {
            if ( withArgument && !argv[0]->IsString() ) throw mup::ParserError( "Argument '" + argv[0]->GetString() + "' is not of type string/id." );
            if ( withArgument ) *ret = tm->realTime( std::string( argv[0]->GetString() ) );
            else *ret = tm->realTimeWithoutArgument();
        }
        const mup::char_type* GetDesc() const { return withArgument ? "realTimeOf( spectrum )" : "realTime()"; }
        std::string tags() const { return "real_time real time spectra foreground fg background bg secondary sfg"; }
        std::string toolTip() const { return withArgument ?
            "Returns the <b>real time</b> (in seconds) of <i>spectrum</i>, or <b>0</b> if not known. Returns <b><font color='red'>error message</font></b> if <i>spectrum</i> was not detected." :
            "Automatically detects single spectrum and returns its <b>real time value</b>. Returns <b><font color='red'>error message</font></b> if multiple or no spectra detected."; }
        mup::IToken* Clone() const { return new RealTime(*this);}
    private: TerminalModel *tm; bool withArgument;
    };
    
    
    // Peak Functions
    class PeakFunction : public mup::ICallback
    {
    public:
        PeakFunction(TerminalModel *model, const char* funName)
          : ICallback(mup::cmFUNC, funName, std::string(funName) == "peakGauss" ? 3 : 1),
            tm(model),
            functionName(funName),
            functionName_energy_arg( std::string(funName) + "( energy )")
      {};
        
        virtual void Eval( mup::ptr_val_type& ret, const mup::ptr_val_type* argv, int a_iArgc ) {
            if ( functionName == "peakGauss" ) {
                if ( !argv[0]->IsNonComplexScalar() || !argv[1]->IsNonComplexScalar() || !argv[2]->IsNonComplexScalar() )
                    throw mup::ParserError( "Arguments for function '" + functionName + "' must be all of type 'scalar value'." );
                *ret = tm->peakGaussianIntegral( argv[0]->GetFloat(), argv[1]->GetFloat(), argv[2]->GetFloat() );
                return;
            }
            if ( !argv[0]->IsNonComplexScalar() ) throw mup::ParserError( "Argument for function '" + functionName + "' must be of type 'scalar value'." );
            const double arg = argv[0]->GetFloat();
            if      ( functionName == "peakArea" )    *ret = tm->peakArea( arg );
            else if ( functionName == "peakMean" )    *ret = tm->peakMean( arg );
            else if ( functionName == "peakSigma" )   *ret = tm->peakSigma( arg );
            else if ( functionName == "peakFwhm" )    *ret = tm->peakFwhm( arg );
            else if ( functionName == "peakAmp" )     *ret = tm->peakAmp( arg );
            else if ( functionName == "peakChi2Dof" ) *ret = tm->peakChi2dof( arg );
        }
        const mup::char_type* GetDesc() const {
            if( functionName == "peakGauss" )
              return "peakGauss( energy, x0, x1 )";
            else
              return functionName_energy_arg.c_str();
        }
        std::string tags() const {
            if      ( functionName == "peakArea" )    return "areas peaks peak area";
            else if ( functionName == "peakMean" )    return "mean peaks average middle mid-point top";
            else if ( functionName == "peakSigma" )   return "sigma peaks gaussian gauss";
            else if ( functionName == "peakFwhm" )    return "fwhm peaks fullwidth halfmaximum full width at half maximum";
            else if ( functionName == "peakAmp" )     return "amplitude peaks amps period periodic";
            else if ( functionName == "peakChi2Dof" ) return "chi squared peaks chisquared chi-squared degree degrees of freedom dof chi^2";
            else if ( functionName == "peakGauss" )   return "gaussian integral peaks euler-poisson";
            else                                      return "";
        }
        std::string toolTip() const {
            if      ( functionName == "peakArea" )    return "Returns a peak’s <b>area</b> at a specified <i>energy</i>. Returns <b>0</b> if no peak was found.";
            else if ( functionName == "peakMean" )    return "Returns a peak’s <b>mean</b> at a specified <i>energy</i>. Returns <b>0</b> if no peak was found.";
            else if ( functionName == "peakSigma" )   return "Returns a peak’s <b>sigma</b> at a specified <i>energy</i>. Returns <b>0</b> if no peak found.</td><td><i>*Peak must be gaussian or else error is returned.";
            else if ( functionName == "peakFwhm" )    return "Returns a peak’s <b>full width at half-maximum (FWHM)</b> at a specified <i>energy</i>. Returns <b>0</b> if no peak found.";
            else if ( functionName == "peakAmp" )     return "Returns a peak’s <b>amplitude</b> at a specified <i>energy</i>. Returns <b>0</b> if no peak found.";
            else if ( functionName == "peakChi2Dof" ) return "Returns a peak’s <b>chi-squared degrees of freedom</b> at a specified <i>energy</i>. Returns <b>0</b> if no peak found.";
            else if ( functionName == "peakGauss" )   return "Returns a peak’s <b>gaussian integral</b> at a specified <i>energy</i> and <i>bounds</i>. Returns <b>0</b> if no peak found or x-boundaries are equivalent.";
            else                                      return "";
        }   
        mup::IToken* Clone() const { return new PeakFunction(*this); }
    
    private:
      TerminalModel *tm;
      std::string functionName;
      std::string functionName_energy_arg;
    };
    
    
    /* 
     Gamma Functions - Obtain specific gamma information from the spectrum.
     
     GammaAtFunctions: These are functions that can be called without having to specify a specific spectrum. In other words, as long as ONLY ONE 
     spectrum is present, then this function can be called (providing the rest of the arguments passed are correct) without having to specify a specific
     spectrum. However, if there are multiple or no spectra detetected, the function returns an error message on the Terminal window.
     
     GammaForFunctions: These are functions that REQUIRE the first argument to be an ID/string of the spectrum to gather specific gamma information.
     If the spectrum specified in the first argument is not detected, then an error message is returned on the Terminal window.
    */
    class GammaAtFunction : public mup::ICallback
    {
    public:
        GammaAtFunction(TerminalModel *model, const std::string &funName)
        :ICallback( mup::cmFUNC, funName.c_str(),
                    (funName == "gammaSum" || funName == "gammaIntegral" || funName == "drfEfficiency") ? 2
                      : (funName == "numGammas" || funName == "gammaMin" || funName == "gammaMax") ? 0 : 1),
        tm(model),
        functionName(funName),
        functionName_no_arg( std::string(funName) + "()" ),
        functionName_energy_arg( std::string(funName) + "( energy )" )
      {};
      
        virtual void Eval( mup::ptr_val_type& ret, const mup::ptr_val_type* argv, int a_iArgc ) {
            const bool functionHasTwoArgs = functionName == "gammaSum" || functionName == "gammaIntegral" || functionName == "drfEfficiency";
            const bool functionHasNoArgs = functionName == "numGammas" || functionName == "gammaMin" || functionName == "gammaMax";
            
          if( functionName == "drfEfficiency" )
          {
            if ( !argv[0]->IsNonComplexScalar() )
              throw mup::ParserError( "First argument for function '" + functionName + "' must be of type 'scalar value'." );
            
            if ( !argv[1]->IsString() )
              throw mup::ParserError( "Second argument for function '" + functionName + "' must be a string giving distance (ex \"1m\")." );
            
            *ret = tm->drfEfficiency( argv[0]->GetFloat(), argv[1]->GetString() );
            return;
          } else if( functionHasTwoArgs )
          {
                if ( !argv[0]->IsNonComplexScalar() || !argv[1]->IsNonComplexScalar() )
                    throw mup::ParserError( "Arguments for function '" + functionName + "' must be all of type 'scalar value'." );
                *ret = (functionName == "gammaSum") ?  tm->gammaSumAt( argv[0]->GetFloat(), argv[1]->GetFloat() ) : tm->gammaIntegralAt( argv[0]->GetFloat(), argv[1]->GetFloat() );
                return;
            }
            else if ( functionHasNoArgs ) {
                if      ( functionName == "numGammas" )  *ret = tm->numGammaChannels( );
                else if ( functionName == "gammaMin" )   *ret = tm->gammaMin( );
                else if ( functionName == "gammaMax" )   *ret = tm->gammaMax( );
                return;
            }
          
          if ( functionName == "drfGeometricEff" )
          {
            if ( !argv[0]->IsString() )
              throw mup::ParserError( "Argument for function '" + functionName + "' must be of a string giving distance (ex \"1.2m\")." );
            *ret = tm->drfGeometricEff( argv[0]->GetString() );
          }
          
          if( !argv[0]->IsNonComplexScalar() )
            throw mup::ParserError( "Argument for function '" + functionName + "' must be of type 'scalar value'." );
          
          const double arg = argv[0]->GetFloat();
          if      ( functionName == "gammaChannel" )   *ret = tm->gammaChannelAt( arg );
            else if ( functionName == "gammaContent" )   *ret = tm->gammaChannelCountAt( arg );
            else if ( functionName == "gammaLower" )     *ret = tm->gammaChannelLowerEnergyAt( arg );
            else if ( functionName == "gammaCenter" )    *ret = tm->gammaChannelCentralEnergyAt( arg );
            else if ( functionName == "gammaUpper" )     *ret = tm->gammaChannelHigherEnergyAt( arg );
            else if ( functionName == "gammaWidth" )     *ret = tm->gammaChannelWidthAt( arg );
            else if ( functionName == "gammaEnergyForChannel" )     *ret = tm->gammaEnergyForChannelAt( arg );
            else if ( functionName == "drfFWHM" )     *ret = tm->drfFWHM( arg );
            else if ( functionName == "drfIntrinsicEff" )     *ret = tm->drfIntrinsicEff( arg );
          
          
          
          
        }
        const mup::char_type* GetDesc() const {
            const bool functionHasNoArgs = functionName == "numGammas" || functionName == "gammaMin" || functionName == "gammaMax";
            
            if ( functionName == "gammaSum" )              return "gammaSum( start_bin(int), end_bin(int) )";
            else if ( functionName == "gammaIntegral" )    return "gammaIntegral( energy_low, energy_high )";
            else if ( functionName == "drfGeometricEff" )  return "drfGeometricEff( distance )";
            else if ( functionName == "drfEfficiency" )    return "drfEfficiency( energy, distance )";
            else if ( functionHasNoArgs )                   return functionName_no_arg.c_str();
            else                                            return functionName_energy_arg.c_str();
        }
        std::string tags() const {
            if      ( functionName == "gammaChannel" )    return "gammas channels gamma_channels";
            else if ( functionName == "gammaContent" )    return "gammas contents gamma_contents channel channels";
            else if ( functionName == "gammaLower" )      return "gammas lower_energy gamma_energy energies energy lower";
            else if ( functionName == "gammaCenter" )     return "gammas center_energy gamma channels central energy energies";
            else if ( functionName == "gammaUpper" )      return "gammas upper_energy gamma_channels channels upper higher bound energy energies";
            else if ( functionName == "gammaWidth" )      return "gammas widths gamma_channels channels width";
            else if ( functionName == "gammaEnergyForChannel" ) return "gammas gamma energy channels";
            else if ( functionName == "numGammas" )       return "number of gammas num_gammas gamma channels numbers";
            else if ( functionName == "gammaMin" )        return "gammas gamma channels minimum min energy energies";
            else if ( functionName == "gammaMax" )        return "gammas gamma channels maximum max energy energies";
            else if ( functionName == "gammaSum" )        return "gammas gamma channels sum add added total two arguments";
            else if ( functionName == "gammaIntegral" )   return "gammas gamma channels integral integrate total lower bound upper higher x sum";
            else if ( functionName == "drfFWHM" )         return "drf detector response FWHM";
            else if ( functionName == "drfIntrinsicEff" ) return "drf detector response intrinsic efficiency";
            else if ( functionName == "drfGeometricEff" ) return "drf detector response geometric efficiency";
            else if ( functionName == "drfEfficiency" )   return "drf detector response efficiency";
            else                                          return "drf detector ";
        }
      
      std::string toolTip() const
      {
        if( functionName == "gammaChannel" )
          return "Returns <b>gamma channel</b> containing energy. If <i>energy</i> is below zero, then 0 is returned. If the <i>energy</i> is above the last channel, the last channel is returned. Automatically detects displayed spectrum, returns <b><font color='red'>error message</font></b> if multiple or no spectra detected.";
        else if ( functionName == "gammaContent" )
          return "Returns <b>gamma channel contents</b> for a specified spectrum in a specified channel. Returns 0 if channel <i>energies</i> is not defined or <i>channel</i> is invalid (too large). Automatically detects displayed spectrum, returns <b><font color='red'>error message</font></b> if multiple or no spectra detected.";
        else if ( functionName == "gammaLower" )
          return "Returns <b>lower energy</b> of specified <i>gamma channel</i> for a specified spectrum. Automatically detects displayed spectrum, returns <b><font color='red'>error message</font></b> if multiple or no spectra detected.";
        else if ( functionName == "gammaCenter" )
          return "Returns <b>central energy</b> of specified <i>gamma channel</i> for a specified spectrum. For last channel, returns width of second-to-last channel. Automatically detects displayed spectrum, returns <b><font color='red'>error message</font></b> if multiple or no spectra detected.";
        else if ( functionName == "gammaUpper" )
          return "Returns <b>energy</b> for a spectra just past energy range the specified <i>channel</i> contains. Returns error if channel is invalid. Automatically detects displayed spectrum, returns <b><font color='red'>error message</font></b> if multiple or no spectra detected.";
        else if ( functionName == "gammaWidth" )
          return "Returns <b>energy width</b> of a channel. If at last channel, then <b>width of second-to-last channel</b> is returned. Automatically detects displayed spectrum, returns <b><font color='red'>error message</font></b> if multiple or no spectra detected.";
        else if( functionName == "gammaEnergyForChannel" )
          return "Returns the energy corresponding to the provided fractional channel; e.x., if you pass in integer, will return lower energy of the channel. Automatically detects displayed spectrum, returns <b><font color='red'>error message</font></b> if multiple or no spectra detected.";
        else if ( functionName == "numGammas" )
          return "Returns <b>minimum number of channels</b> of channel energies or gamma counts for spectrum. Returns 0 if neither is defined. Automatically detects displayed spectrum, returns <b><font color='red'>error message</font></b> if multiple or no spectra detected.";
        else if ( functionName == "gammaMin" )
          return "Returns <b>minimum gamma energy</b>. Automatically detects displayed spectrum, returns <b><font color='red'>error message</font></b> if multiple or no spectra detected.";
        else if ( functionName == "gammaMax" )
          return "Returns <b>maximum gamma energy</b>. Automatically detects displayed spectrum, returns <b><font color='red'>error message</font></b> if multiple or no spectra detected.";
        else if ( functionName == "gammaSum" )
          return "Get the <b>sum of gamma channel contents</b> for all channels in between (inclusive) <i>start_bin</i> and <i>end_bin</i> for a spectrum. Returns 0 if start_bin too large or gamma counts invalid. If end_bin too large, then it will be clamped to number of channels. Automatically detects displayed spectrum, returns <b><font color='red'>error message</font></b> if multiple or no spectra detected.";
        else if ( functionName == "gammaIntegral" )
          return "Get <b>integral of gamma counts</b> between <i>energy_low</i> and <i>energy_high</i> for a spectrum. Returns <b>0</b> if channel energies or gamma counts invalid. Automatically detects displayed spectrum, returns <b><font color='red'>error message</font></b> if multiple or no spectra detected.";
        else if( functionName == "drfFWHM" )
          return "Returns the full-width-at-half-maximum according to the current detector response function, for the specified energy.";
        else if( functionName == "drfIntrinsicEff" )
          return "Returns the intrinsic efficiency (i.e., probability of gamma striking face of detector contributing towards a peak) for the specified energy.";
        else if( functionName == "drfGeometricEff" )
          return "Returns the fraction of gammas emitted from a point source at the specified distance, that will impinge upon the detector face.";
        else if( functionName == "drfEfficiency" )
          return "Returns the fraction of gammas, at the specified energy and distance, that will contribute to a photo-peak.";
        
        return "";
      }
      
        mup::IToken* Clone() const { return new GammaAtFunction(*this); }
    private:
      TerminalModel *tm;
      std::string functionName;
      std::string functionName_no_arg;
      std::string functionName_energy_arg;
    };
    
    class GammaForFunction : public mup::ICallback
    {
    public:
        GammaForFunction(TerminalModel *model, const char* funName)
        :ICallback(mup::cmFUNC, funName,
                   std::string(funName) == "gammaSumFor" || std::string(funName) == "gammaIntegralFor" ? 3 : std::string(funName) == "numGammasFor" || std::string(funName) == "gammaMinFor" || std::string(funName) == "gammaMaxFor" ? 1 : 2),
        tm(model),
      functionName_spectrum_arg( std::string(funName) + "( spectrum )" ),
      functionName_spec_energy_arg( std::string(funName) + "( spectrum, energy )" ),
      functionName(funName)
      {};
        
        virtual void Eval( mup::ptr_val_type& ret, const mup::ptr_val_type* argv, int a_iArgc ) {
            const bool functionHasThreeArgs = functionName == "gammaSumFor" || functionName == "gammaIntegralFor";
            const bool functionHasOneArg = functionName == "numGammasFor" || functionName == "gammaMinFor" || functionName == "gammaMaxFor";
            
            if ( !argv[0]->IsString() ) throw mup::ParserError( "Argument 'spectrum' for function '" + functionName + "' must be of type 'string'." );
            const std::string& arg0 = argv[0]->GetString();
            
            if ( functionHasThreeArgs ) {
                if ( !argv[1]->IsNonComplexScalar() || !argv[2]->IsNonComplexScalar() )
                    throw mup::ParserError( "Arguments 2 and 3 for function '" + functionName + "' must be all of type 'scalar value'." );
                *ret = (functionName == "gammaSumFor") ?  tm->gammaSumFor( arg0, argv[1]->GetFloat(), argv[2]->GetFloat() ) : tm->gammaIntegralFor( arg0, argv[1]->GetFloat(), argv[2]->GetFloat() );
                return;
            }
            else if ( functionHasOneArg ) {
                if      ( functionName == "numGammasFor" )  *ret = tm->numGammaChannelsFor( arg0 );
                else if ( functionName == "gammaMinFor" )   *ret = tm->gammaMinFor( arg0 );
                else if ( functionName == "gammaMaxFor" )   *ret = tm->gammaMaxFor( arg0 );
                return;
            }
            if ( !argv[1]->IsNonComplexScalar() )
                throw mup::ParserError( "Arguments 2 and 3 for function '" + functionName + "' must be of type 'scalar value'." );
            const double arg1 = argv[1]->GetFloat();
            if      ( functionName == "gammaChannelFor" )   *ret = tm->gammaChannelFor( arg0, arg1 );
            else if ( functionName == "gammaContentFor" )   *ret = tm->gammaChannelCountFor( arg0, arg1 );
            else if ( functionName == "gammaLowerFor" )     *ret = tm->gammaChannelLowerEnergyFor( arg0, arg1 );
            else if ( functionName == "gammaCenterFor" )    *ret = tm->gammaChannelCentralEnergyFor( arg0, arg1 );
            else if ( functionName == "gammaUpperFor" )     *ret = tm->gammaChannelHigherEnergyFor( arg0, arg1 );
            else if ( functionName == "gammaWidthFor" )     *ret = tm->gammaChannelWidthFor( arg0, arg1 );
        }
        const mup::char_type* GetDesc() const {
            const bool functionHasOneArg = functionName == "numGammasFor" || functionName == "gammaMinFor" || functionName == "gammaMaxFor";
            
            if ( functionName == "gammaSumFor" )              return "gammaSumFor( spectrum, start_bin, end_bin )";
            else if ( functionName == "gammaIntegralFor" )    return "gammaIntegralFor( spectrum, energy_low, energy_high )";
            else if ( functionHasOneArg )                   return functionName_spectrum_arg.c_str();
            else                                            return functionName_spec_energy_arg.c_str();
        }
        std::string tags() const {
            if      ( functionName == "gammaChannelFor" )   return "gammas channels gamma_channels for spectrum spectra";
            else if ( functionName == "gammaContentFor" )   return "gammas contents gamma_contents channel channels for spectrum spectra";
            else if ( functionName == "gammaLowerFor" )     return "gammas lower_energy gamma_energy energies energy lower for spectrum spectra";
            else if ( functionName == "gammaCenterFor" )    return "gammas center_energy gamma channels central energy energies for spectrum spectra";
            else if ( functionName == "gammaUpperFor" )     return "gammas upper_energy gamma_channels channels upper higher bound energy energies for spectrum spectra";
            else if ( functionName == "gammaWidthFor" )     return "gammas widths gamma_channels channels width for spectrum spectra";
            else if ( functionName == "numGammasFor" )      return "number of gammas num_gammas gamma channels numbers for spectrum spectra";
            else if ( functionName == "gammaMinFor" )       return "gammas gamma channels minimum min energy energies for spectrum spectra";
            else if ( functionName == "gammaMaxFor" )       return "gammas gamma channels maximum max energy energies for spectrum spectra";
            else if ( functionName == "gammaSumFor" )       return "gammas gamma channels sum add added total two arguments for spectrum spectra";
            else if ( functionName == "gammaIntegralFor" )  return "gammas gamma channels integral integrate total lower bound upper higher x sum for spectrum spectra";
            else                                             return "";
        }
        std::string toolTip() const {
            if      ( functionName == "gammaChannelFor" )   return "Returns <b>gamma channel</b> containing energy. If <i>energy</i> is below zero, then 0 is returned. If the <i>energy</i> is above the last channel, the last channel is returned.";
            else if ( functionName == "gammaContentFor" )   return "Returns <b>gamma channel contents</b> for a specified spectrum in a specified channel. Returns 0 if channel <i>energies</i> is not defined or <i>channel</i> is invalid (too large).";
            else if ( functionName == "gammaLowerFor" )     return "Returns <b>lower energy</b> of specified <i>gamma channel</i> for a specified spectrum.";
            else if ( functionName == "gammaCenterFor" )    return "Returns <b>central energy</b> of specified <i>gamma channel</i> for a specified spectrum. For last channel, returns width of second-to-last channel.";
            else if ( functionName == "gammaUpperFor" )     return "Returns <b>energy</b> for a spectra just past energy range the specified <i>channel</i> contains. Returns error if channel is invalid.";
            else if ( functionName == "gammaWidthFor" )     return "Returns <b>energy width</b> of a channel. If at last channel, then <b>width of second-to-last channel</b> is returned.";
            else if ( functionName == "numGammasFor" )      return "Returns <b>minimum number of channels</b> of channel energies or gamma counts for spectrum. Returns 0 if neither is defined.";
            else if ( functionName == "gammaMinFor" )       return "Returns <b>minimum gamma energy</b>.";
            else if ( functionName == "gammaMaxFor" )       return "Returns <b>maximum gamma energy</b>.";
            else if ( functionName == "gammaSumFor" )       return "Get the <b>sum of gamma channel contents</b> for all channels in between (inclusive) <i>start_bin</i> and <i>end_bin</i> for a spectrum. Returns 0 if start_bin too large or gamma counts invalid. If end_bin too large, then it will be clamped to number of channels.";
            else if ( functionName == "gammaIntegralFor" )  return "Get <b>integral of gamma counts</b> between <i>energy_low</i> and <i>energy_high</i> for a spectrum. Returns <b>0</b> if channel energies or gamma counts invalid.";
            else                                         return "";
        }
        mup::IToken* Clone() const { return new GammaForFunction(*this); }
    private:
      TerminalModel *tm;
      std::string functionName;
      std::string functionName_spectrum_arg;
      std::string functionName_spec_energy_arg;
    };
    
    
    // Statistic Functions
    class PValueFromChiSquare : public mup::ICallback
    {
    public:
        PValueFromChiSquare() :ICallback(mup::cmFUNC, "pValueFromChiSquare", 2) {};
        
        virtual void Eval( mup::ptr_val_type& ret, const mup::ptr_val_type* argv, int a_iArgc ) {
            if ( !argv[0]->IsNonComplexScalar() || !argv[1]->IsNonComplexScalar() )
                throw mup::ParserError( "Arguments for function 'pValueFromChiSquare' must be all of type 'scalar value'." );
            *ret = pValueFromChiSquare( argv[0]->GetFloat(), argv[1]->GetFloat() );
        }
        const mup::char_type* GetDesc() const { return "pValueFromChiSquare( chiScore, df )"; }
        std::string tags() const { return "p-value pvals chi2 chi^2 chisquared chi-squared pvalues p value values stats statistics"; }
        std::string toolTip() const { return "Returns the p-value from a specific chi-score and degrees of freedom."; }
        mup::IToken* Clone() const { return new PValueFromChiSquare(*this); }
    };
}



// Main Methods
// Constructor
TerminalModel::TerminalModel( InterSpec* viewer )
{
    m_parser.reset( new mup::ParserX() );
    m_variables = VariableMap();
    
    m_viewer = viewer;
    
  m_foregroundHistogram = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  m_backgroundHistogram = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Background );
  m_secondaryHistogram  = m_viewer->displayedHistogram( SpecUtils::SpectrumType::SecondForeground );
    
    /* Note: The order you add functions/commands WILL MATTER for how it is displayed in the search menu!
             Commands are listed by the order the addCommand/addFunction methods are called. */
    
    m_commandMap = CommandMap();
    addDropDownListHeader( "Peak Commands" );
    addCommand( "searchPeak", SearchPeakCommand );
    addCommand( "deletePeak", DeletePeakCommand );
    addCommand( "refitPeak" , RefitPeakCommand  );
    
    addDropDownListHeader( "Spectrum Viewer Commands" );
    addCommand( "setRange"  , SetEnergyRangeCommand );
    addCommand( "setYRange"  , SetYaxisRangeCommand );
    addCommand( "setLogYAxisMin"  , SetLogYAxisMinCommand );
    addDropDownListSeparator();
    
    addDropDownListHeader( "Nuclide Commands" );
    addCommand( "setNuclide", SetNuclideCommand );
    addDropDownListSeparator();
    
    addDropDownListHeader( "Variable Commands" );
    addCommand( "showDefinedVariables"    , VarmapCommand );
    addCommand( "clearVariable"  , ClearVarCommand   );
    addDropDownListSeparator();
    
    addDropDownListHeader( "Terminal Commands" );
    addCommand( "darken"    , DarkenCommand     );
    addCommand( "lighten"   , LightenCommand    );
    addDropDownListSeparator();
    
    addCommand( "saveFile"  , SaveCommand       );      // currently unfinished
    
    // Define non-built-in operators
    m_parser->RemoveOprt( "//" );   // remove the original "//" operator (string concatenation), I don't think it's that useful
    m_parser->DefineOprt( new Modulus( "%" ) );
    m_parser->DefineOprt( new Modulus( "mod" ) );
    m_parser->DefineOprt( new Division( "//" ) );
    m_parser->DefineOprt( new Division( "div" ) );
 
    // Define non-built-in functions
    // Define Live/Real time functions
    addDropDownListHeader( "Live Time / Real Time Functions" );
    LiveTime* liveTimeFun = new LiveTime(this);              addFunction( liveTimeFun, liveTimeFun->tags(), liveTimeFun->toolTip() );       // live time for fg, bg, secondary fg
    LiveTime* liveTimeFunWArg = new LiveTime(this, true);  addFunction( liveTimeFunWArg, liveTimeFunWArg->tags(), liveTimeFunWArg->toolTip() ); // live time function that auto detects live time value
    RealTime* realTimeFun = new RealTime(this);              addFunction( realTimeFun, realTimeFun->tags(), realTimeFun->toolTip() );       // real time for fg, bg, secondary fg
    RealTime* realTimeFunWArg = new RealTime(this, true);  addFunction( realTimeFunWArg, realTimeFunWArg->tags(), realTimeFunWArg->toolTip() ); // real time funciton that auto detects real time value
    
    addDropDownListHeader( "Nuclide Functions" );
    NuclideIntensity* nuclideIntensityFun = new NuclideIntensity(this);  addFunction( nuclideIntensityFun, nuclideIntensityFun->tags(), nuclideIntensityFun->toolTip() );
    NuclideIntensity* nuclideIntensityFunWArg = new NuclideIntensity(this, true);   addFunction( nuclideIntensityFunWArg, nuclideIntensityFunWArg->tags(), nuclideIntensityFunWArg->toolTip() );
    
    // Define Gamma functions
    addDropDownListHeader( "Gamma Functions - with Spectrum argument" );
    GammaForFunction* numGammasForFunc = new GammaForFunction(this, "numGammasFor");            addFunction( numGammasForFunc, numGammasForFunc->tags(), numGammasForFunc->toolTip() );
    GammaForFunction* gammaChannelForFunc = new GammaForFunction(this, "gammaChannelFor");      addFunction( gammaChannelForFunc, gammaChannelForFunc->tags(), gammaChannelForFunc->toolTip() );
    GammaForFunction* gammaContentForFunc = new GammaForFunction(this, "gammaContentFor");      addFunction( gammaContentForFunc, gammaContentForFunc->tags(), gammaContentForFunc->toolTip() );
    GammaForFunction* gammaLowerForFunc = new GammaForFunction(this, "gammaLowerFor");          addFunction( gammaLowerForFunc, gammaLowerForFunc->tags(), gammaLowerForFunc->toolTip() );
    GammaForFunction* gammaCenterForFunc = new GammaForFunction(this, "gammaCenterFor");        addFunction( gammaCenterForFunc, gammaCenterForFunc->tags(), gammaCenterForFunc->toolTip() );
    GammaForFunction* gammaUpperForFunc = new GammaForFunction(this, "gammaUpperFor");          addFunction( gammaUpperForFunc, gammaUpperForFunc->tags(), gammaUpperForFunc->toolTip() );
    GammaForFunction* gammaWidthForFunc = new GammaForFunction(this, "gammaWidthFor");          addFunction( gammaWidthForFunc, gammaWidthForFunc->tags(), gammaWidthForFunc->toolTip() );
    GammaForFunction* gammaIntegralForFunc = new GammaForFunction(this, "gammaIntegralFor");    addFunction( gammaIntegralForFunc, gammaIntegralForFunc->tags(), gammaIntegralForFunc->toolTip() );
    GammaForFunction* gammaSumForFunc = new GammaForFunction(this, "gammaSumFor");              addFunction( gammaSumForFunc, gammaSumForFunc->tags(), gammaSumForFunc->toolTip() );
    GammaForFunction* gammaMinForFunc = new GammaForFunction(this, "gammaMinFor");              addFunction( gammaMinForFunc, gammaMinForFunc->tags(), gammaMinForFunc->toolTip() );
    GammaForFunction* gammaMaxForFunc = new GammaForFunction(this, "gammaMaxFor");              addFunction( gammaMaxForFunc, gammaMaxForFunc->tags(), gammaMaxForFunc->toolTip() );
    addDropDownListSeparator();
    
    addDropDownListHeader( "Gamma Functions - without Spectrum argument" );
    GammaAtFunction* numGammasFunc = new GammaAtFunction(this, "numGammas");              addFunction( numGammasFunc, numGammasFunc->tags(), numGammasFunc->toolTip() );
    GammaAtFunction* gammaChannelFunc = new GammaAtFunction(this, "gammaChannel");        addFunction( gammaChannelFunc, gammaChannelFunc->tags(), gammaChannelFunc->toolTip() );
    GammaAtFunction* gammaContentFunc = new GammaAtFunction(this, "gammaContent");        addFunction( gammaContentFunc, gammaContentFunc->tags(), gammaContentFunc->toolTip() );
    GammaAtFunction* gammaLowerFunc = new GammaAtFunction(this, "gammaLower");            addFunction( gammaLowerFunc, gammaLowerFunc->tags(), gammaLowerFunc->toolTip() );
    GammaAtFunction* gammaCenterFunc = new GammaAtFunction(this, "gammaCenter");          addFunction( gammaCenterFunc, gammaCenterFunc->tags(), gammaCenterFunc->toolTip() );
    GammaAtFunction* gammaUpperFunc = new GammaAtFunction(this, "gammaUpper");            addFunction( gammaUpperFunc, gammaUpperFunc->tags(), gammaUpperFunc->toolTip() );
  
  GammaAtFunction* gammaWidthFunc = new GammaAtFunction(this, "gammaWidth");
  addFunction( gammaWidthFunc, gammaWidthFunc->tags(), gammaWidthFunc->toolTip() );
  
  GammaAtFunction* gammaEnergyForChannelFunc = new GammaAtFunction(this, "gammaEnergyForChannel");
  addFunction( gammaEnergyForChannelFunc, gammaEnergyForChannelFunc->tags(), gammaEnergyForChannelFunc->toolTip() );
  
  
    GammaAtFunction* gammaIntegralFunc = new GammaAtFunction(this, "gammaIntegral");      addFunction( gammaIntegralFunc, gammaIntegralFunc->tags(), gammaIntegralFunc->toolTip() );
    GammaAtFunction* gammaSumFunc = new GammaAtFunction(this, "gammaSum");                addFunction( gammaSumFunc, gammaSumFunc->tags(), gammaSumFunc->toolTip() );
    GammaAtFunction* gammaMinFunc = new GammaAtFunction(this, "gammaMin");                addFunction( gammaMinFunc, gammaMinFunc->tags(), gammaMinFunc->toolTip() );
    GammaAtFunction* gammaMaxFunc = new GammaAtFunction(this, "gammaMax");                addFunction( gammaMaxFunc, gammaMaxFunc->tags(), gammaMaxFunc->toolTip() );
    
    // Define Peak functions
    addDropDownListHeader( "Peak Functions" );
    PeakFunction* peakAreaFun = new PeakFunction(this, "peakArea");         addFunction( peakAreaFun, peakAreaFun->tags(), peakAreaFun->toolTip() );
    PeakFunction* peakMeanFun = new PeakFunction(this, "peakMean");         addFunction( peakMeanFun, peakMeanFun->tags(), peakMeanFun->toolTip() );
    PeakFunction* peakSigmaFun = new PeakFunction(this, "peakSigma");       addFunction( peakSigmaFun, peakSigmaFun->tags(), peakSigmaFun->toolTip() );
    PeakFunction* peakFwhmFun = new PeakFunction(this, "peakFwhm");         addFunction( peakFwhmFun, peakFwhmFun->tags(), peakFwhmFun->toolTip() );
    PeakFunction* peakAmpFun = new PeakFunction(this, "peakAmp");           addFunction( peakAmpFun, peakAmpFun->tags(), peakAmpFun->toolTip() );
    PeakFunction* peakChi2DofFun = new PeakFunction(this, "peakChi2Dof");   addFunction( peakChi2DofFun, peakChi2DofFun->tags(), peakChi2DofFun->toolTip() );
    PeakFunction* peakGaussFun = new PeakFunction(this, "peakGauss");       addFunction( peakGaussFun, peakGaussFun->tags(), peakGaussFun->toolTip() );
    
  
  //Define FWHM functions
  addDropDownListHeader( "Detector Response Functions" );
  GammaAtFunction *drfFWHMFunc = new GammaAtFunction(this, "drfFWHM");
  addFunction( drfFWHMFunc, drfFWHMFunc->tags(), drfFWHMFunc->toolTip() );
  
  GammaAtFunction *drfIntrinsicEffFunc = new GammaAtFunction(this, "drfIntrinsicEff");
  addFunction( drfIntrinsicEffFunc, drfIntrinsicEffFunc->tags(), drfIntrinsicEffFunc->toolTip() );
  
  GammaAtFunction *drfGeometricEffFunc = new GammaAtFunction(this, "drfGeometricEff");
  addFunction( drfGeometricEffFunc, drfGeometricEffFunc->tags(), drfGeometricEffFunc->toolTip() );
  
  GammaAtFunction *drfEfficiencyFunc = new GammaAtFunction(this, "drfEfficiency");
  addFunction( drfEfficiencyFunc, drfEfficiencyFunc->tags(), drfEfficiencyFunc->toolTip() );
  
  
    // Define non-built-in statistical functions
    addDropDownListHeader( "Statistical Functions" );
    PValueFromChiSquare* pValFun = new PValueFromChiSquare;      addFunction( pValFun, pValFun->tags(), pValFun->toolTip() );
    
    // Define non-built-in statistical functions
    addDropDownListHeader( "Trigonometric Functions" );
    ToAngle* toAngleFun = new ToAngle;      addFunction( toAngleFun, toAngleFun->tags(), toAngleFun->toolTip() );
    ToRadians* toRadiansFun = new ToRadians;      addFunction( toRadiansFun, toRadiansFun->tags(), toRadiansFun->toolTip() );
    
    // Define non-built-in constants
    m_parser->DefineConst( "_sqrt2", mup::Value(M_SQRT2) );
    
    // Define non-built-in string constants
    m_parser->DefineConst("foreground", mup::Value("fg"));             m_parser->DefineConst("fg", mup::Value("fg"));     // These are string constants for the arguments for the live time
    m_parser->DefineConst("secondaryforeground", mup::Value("sfg"));   m_parser->DefineConst("sfg", mup::Value("sfg"));   // and real time functions; avoids use of having to use double-quotes
    m_parser->DefineConst("background", mup::Value("bg"));             m_parser->DefineConst("bg", mup::Value("bg"));     // in the arguments for calling the functions
    
    m_parser->DefineConst("xray", mup::Value("xray"));
    m_parser->DefineConst("gamma", mup::Value("gamma"));
    m_parser->DefineConst("alpha", mup::Value("alpha"));
    m_parser->DefineConst("beta", mup::Value("beta"));
}

// Destructor
TerminalModel::~TerminalModel()
{
    for (auto i : m_variables)       // delete all the pointers to doubles inside the variable map
        delete i.second;
}

// Set the Spectrum Viewer to extract data from it
void TerminalModel::setViewer( InterSpec* spectViewer )
{
    m_viewer = spectViewer;
}

// Determines input type of a user-inputted string
TerminalModel::InputType TerminalModel::inputType(const std::string& input)
{
  std::regex variableAssignmentRegex( ns_variableAssignmentRegexArg );
  
    if (std::regex_match(input, variableAssignmentRegex))
        return VariableAssignment;
    
    else if (std::regex_match(input, std::regex(commandRegex()) )) {
        std::string command = input.substr(0, input.find("("));
        
        // If the command being called is a keyword for an operation, then it is an OPERATION
      std::regex keywordRegex( ns_keywordRegexArg );
        if (std::regex_match(command, keywordRegex)) return Operation;
        else return Command;
    }
    // Default input
    return Operation;
}

// Evalutes inputted string into a result string
std::string TerminalModel::evaluate(std::string input)
{
    SpecUtils::trim(input);
    
    if (!input.empty()) {
        const TerminalModel::InputType type = inputType(input);             // Get input-type of user's input
        
        if (type == Command)
            return doCommand(input);
            
        else if (type == VariableAssignment)
            return assignVariable(input);
            
        else if (type == Operation) {
            m_parser->SetExpr( input );                                         // Sets the expression for the parser into the current input
            
            try {
                const double result = m_parser->Eval().GetFloat();              // Get the value of the evaluated input
                return convertDoubleTypeToString(result);                       // Convert result into string format  
            } catch ( const mup::ParserError& e ) {
                // These are errors specific to the parser
              return "Error code " + std::to_string(e.GetCode()) + ": " + e.GetMsg();
            } catch (std::runtime_error re) {
              return "Runtime error: " + std::string(re.what());
            } catch (std::exception e) {                                        // catch any other possible exceptions (non-parser specific)
              return "Error: " + std::string(e.what());
            }
        }
    }
    return "";                                                              // print out empty line if input is empty
}

void TerminalModel::addDropDownListItem( const std::string& item, const std::string& tags, const std::string& toolTip ) {
    m_dropDownList.push_back( CommandHelperTuple( item, tags, toolTip ) );
}
void TerminalModel::addDropDownListHeader( const std::string& header ) { addDropDownListItem( "<header:" + header ); }
void TerminalModel::addDropDownListSeparator()                         { addDropDownListItem( "<separator>" ); }
TerminalModel::CommandHelperList TerminalModel::commandsFunctionsList()               { return m_dropDownList; }


// Non-built-in Function Methods (for parser)
// Adds function definition to the parser
void TerminalModel::addFunction( mup::ICallback* function, const std::string& tags, const std::string& toolTip )
{
    addDropDownListItem( function->GetDesc(), tags, toolTip );
    m_parser->DefineFun( function );
}

// Function helpers
// Updates all the spectra
void TerminalModel::updateHistograms()
{
    m_foregroundHistogram = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
    m_backgroundHistogram = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Background );
    m_secondaryHistogram  = m_viewer->displayedHistogram( SpecUtils::SpectrumType::SecondForeground );
}

// These methods help extract the live/real time of the foreground, secondary foreground, background
double TerminalModel::foregroundLiveTime() {
    m_foregroundHistogram = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
    if (m_foregroundHistogram == nullptr) throw mup::ParserError( "Foreground could not be detected" );
    return m_foregroundHistogram->live_time();
}
double TerminalModel::foregroundRealTime() {
    m_foregroundHistogram = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
    if (m_foregroundHistogram == nullptr) throw mup::ParserError( "Foreground could not be detected" );
    return m_foregroundHistogram->real_time();
}
double TerminalModel::secondaryForegroundLiveTime() {
    m_foregroundHistogram = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
    if (m_secondaryHistogram == nullptr) throw mup::ParserError( "Secondary foreground could not be detected" );
    return m_secondaryHistogram->live_time();
}
double TerminalModel::secondaryForegroundRealTime() {
    m_secondaryHistogram = m_viewer->displayedHistogram( SpecUtils::SpectrumType::SecondForeground );
    if (m_secondaryHistogram == nullptr) throw mup::ParserError( "Secondary foreground could not be detected" );
    return m_secondaryHistogram->real_time();
}
double TerminalModel::backgroundLiveTime() {
    m_backgroundHistogram = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Background );
    if (m_backgroundHistogram == nullptr) throw mup::ParserError( "Background could not be detected" );
    return m_backgroundHistogram->live_time();
}
double TerminalModel::backgroundRealTime() {
    m_backgroundHistogram = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Background );
    if (m_backgroundHistogram == nullptr) throw mup::ParserError( "Background could not be detected" );
    return m_backgroundHistogram->real_time();
}

double TerminalModel::nuclideIntensity( const double energy )
{
  ReferencePhotopeakDisplay *nuclideReference = m_viewer->referenceLinesWidget();
  if( !nuclideReference )
    throw mup::ParserError( "ReferencePhotopeakDisplay is currently NULL." );
  
  const ReferenceLineInfo& nuclideInfo = nuclideReference->currentlyShowingNuclide();
        
  if( nuclideInfo.m_validity != ReferenceLineInfo::InputValidity::Valid )
    throw mup::ParserError( "No nuclide is currently being shown." );
  
  
  double nearest_intensity = -1.0, nearest_distance = 1.0E6;
  
  for( const ReferenceLineInfo::RefLine &line : nuclideInfo.m_ref_lines )
  {
    const double dist = fabs(line.m_energy - energy);
    if( dist < nearest_distance )
    {
      nearest_intensity = line.m_normalized_intensity; // TODO: do we want line.m_decay_intensity here?
      nearest_distance = dist;
    }
  }//
  
  if( nearest_distance > 0.1 )
    throw mup::ParserError( "Intensity value not found for nuclide at specified intensity." );
  
  return nearest_intensity;
}//TerminalModel::nuclideIntensity(...)


double TerminalModel::nuclideIntensityForParticle(  std::string particle, const double energy )
{
  ReferencePhotopeakDisplay *nuclideReference = m_viewer->referenceLinesWidget();
  if( !nuclideReference )
    throw mup::ParserError( "ReferencePhotopeakDisplay is currently NULL." );
  
  const ReferenceLineInfo& nuclideInfo = nuclideReference->currentlyShowingNuclide();
  
  if( nuclideInfo.m_validity != ReferenceLineInfo::InputValidity::Valid )
    throw mup::ParserError( "No nuclide is currently being shown." );
  
  SpecUtils::to_lower_ascii( particle );
  if( (particle == "x-ray") || (particle == "x ray") )
    particle = "xray";
  
  ReferenceLineInfo::RefLine::Particle particle_type;
  
  if( particle == "xray" )
    particle_type = ReferenceLineInfo::RefLine::Particle::Xray;
  else if( particle == "alpha" )
    particle_type = ReferenceLineInfo::RefLine::Particle::Alpha;
  else if( particle == "beta" )
    particle_type = ReferenceLineInfo::RefLine::Particle::Beta;
  else if( particle == "gamma")
    particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
  else
    throw mup::ParserError( "Particle type must be one of {'xray', 'alpha', 'beta', 'gamma'}" );
  
  double nearest_intensity = -1.0, nearest_distance = 1.0E6;
  
  for( const ReferenceLineInfo::RefLine &line : nuclideInfo.m_ref_lines )
  {
    if( line.m_particle_type != particle_type )
      continue;
    
    const double dist = fabs(line.m_energy - energy);
    if( dist < nearest_distance )
    {
      nearest_intensity = line.m_normalized_intensity; // TODO: do we want line.m_decay_intensity here?
      nearest_distance = dist;
    }
  }//
  
  if( nearest_distance > 0.1 )
    throw mup::ParserError( "Intensity value not found for nuclide at specified intensity." );
  
  return nearest_intensity;
}//nuclideIntensityForParticle(...)


double TerminalModel::nuclideEnergy( const double intensity )
{
  ReferencePhotopeakDisplay *nuclideReference = m_viewer->referenceLinesWidget();
  if( !nuclideReference )
    throw mup::ParserError( "ReferencePhotopeakDisplay is currently NULL." );
  
  const ReferenceLineInfo& nuclideInfo = nuclideReference->currentlyShowingNuclide();
  
  if( nuclideInfo.m_validity != ReferenceLineInfo::InputValidity::Valid )
    throw mup::ParserError( "No nuclide is currently being shown." );
  
  double nearest_energy = -1.0, nearest_distance = 1.0E6;
  for( const ReferenceLineInfo::RefLine &line : nuclideInfo.m_ref_lines )
  {
    const double dist = fabs(line.m_normalized_intensity - intensity); // TODO: do we want line.m_decay_intensity here?
    if( dist < nearest_distance )
    {
      nearest_energy = line.m_energy;
      nearest_distance = dist;
    }
  }//
  
  if( nearest_distance > 0.01 )
    throw mup::ParserError( "Intensity value within 0.01 of " + std::to_string(intensity) + " not found." );
  
  return nearest_energy;
}//nuclideEnergy(...)


// Return the live time of either the foreground, background, or secondary foreground (determined in argument by user)
// Throws error if invalid name is provided.
double TerminalModel::liveTime( const std::string& argument )
{
    const std::string& arg (argument);
    
    m_foregroundHistogram = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
    m_backgroundHistogram = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Background );
    m_secondaryHistogram = m_viewer->displayedHistogram( SpecUtils::SpectrumType::SecondForeground );
    
    if ( std::regex_match( arg, std::regex( "^(\\s*(foreground|fg)\\s*)$", std::regex::icase ) ) )    // foreground
        return foregroundLiveTime();
    else if ( std::regex_match( arg, std::regex( "^(\\s*(secondary|sfg|secondaryforeground|secondforeground)\\s*)$", std::regex::icase ) ) )  // second foreground
        return secondaryForegroundLiveTime();
    else if ( std::regex_match( arg, std::regex( "^(\\s*(background|bg|background|back)\\s*)$", std::regex::icase ) ) )  // background
        return backgroundLiveTime();
    
    // Throw an ecINVALID_NAME error if argument user provided is not valid
    throw mup::ParserError( "Invalid argument for function 'live_time'" );
}

// Return the the live time of only single live time value.
// Throws error if multiple live time values are found.
double TerminalModel::liveTimeWithoutArgument()
{
    updateHistograms();
    
    if (!m_foregroundHistogram && !m_secondaryHistogram && !m_backgroundHistogram)
        throw mup::ParserError( "No spectra detected. Please add a spectrum to use." );
    
    if (m_foregroundHistogram && !m_secondaryHistogram && !m_backgroundHistogram)
        return m_foregroundHistogram->live_time();
    else if (!m_foregroundHistogram && m_secondaryHistogram && !m_backgroundHistogram)
        return m_secondaryHistogram->live_time();
    else if (!m_foregroundHistogram && !m_secondaryHistogram && m_backgroundHistogram)
        return m_backgroundHistogram->live_time();
    
    else if ( m_foregroundHistogram && m_secondaryHistogram && m_backgroundHistogram &&
                m_foregroundHistogram->live_time() == m_secondaryHistogram->live_time() && m_secondaryHistogram->live_time() == m_backgroundHistogram->live_time() )
        return m_foregroundHistogram->live_time();
    else if ( m_foregroundHistogram && !m_secondaryHistogram && m_backgroundHistogram && m_foregroundHistogram->live_time() == m_backgroundHistogram->live_time() )
        return m_foregroundHistogram->live_time();
    else if ( m_foregroundHistogram && m_secondaryHistogram && !m_backgroundHistogram && m_foregroundHistogram->live_time() == m_secondaryHistogram->live_time() )
        return m_foregroundHistogram->live_time();
    else if ( !m_foregroundHistogram && m_secondaryHistogram && m_backgroundHistogram && m_backgroundHistogram->live_time() == m_secondaryHistogram->live_time() )
        return m_backgroundHistogram->live_time();
    
    else
        throw mup::ParserError( "Multiple spectra detected. Please specify a spectrum to use." );
}

// Return the real time of either the foreground, background, or secondary foreground (determined in argument by user)
// Throws error if multiple live time values are found
double TerminalModel::realTime( const std::string& argument )
{
    const std::string& arg (argument);
    
    if ( std::regex_match( arg, std::regex( "^(\\s*(foreground|fg)\\s*)$", std::regex::icase ) ) )    // foreground
        return foregroundRealTime();
    else if ( std::regex_match( arg, std::regex( "^(\\s*(secondary|sfg|secondaryforeground|secondforeground)\\s*)$", std::regex::icase ) ) )  // second foreground
        return secondaryForegroundRealTime();
    else if ( std::regex_match( arg, std::regex( "^(\\s*(background|bg|background|back)\\s*)$", std::regex::icase ) ) )  // background
        return backgroundRealTime();
    
    // Throw an ecINVALID_NAME error if argument user provided is not valid
    throw mup::ParserError ( "Invalid argument for function 'real_time'" );
}

// Return the the real time of only single real time value.
// Throws error if multiple real time values are found.
double TerminalModel::realTimeWithoutArgument()
{
    updateHistograms();
    
    if (!m_foregroundHistogram && !m_secondaryHistogram && !m_backgroundHistogram)
        throw mup::ParserError( "No spectra detected. Please add a spectrum to use." );
    
    if (m_foregroundHistogram && !m_secondaryHistogram && !m_backgroundHistogram)
        return m_foregroundHistogram->real_time();
    else if (!m_foregroundHistogram && m_secondaryHistogram && !m_backgroundHistogram)
        return m_secondaryHistogram->real_time();
    else if (!m_foregroundHistogram && !m_secondaryHistogram && m_backgroundHistogram)
        return m_backgroundHistogram->real_time();
    
    else if ( m_foregroundHistogram && m_secondaryHistogram && m_backgroundHistogram &&
             m_foregroundHistogram->real_time() == m_secondaryHistogram->real_time() && m_secondaryHistogram->real_time() == m_backgroundHistogram->real_time() )
        return m_foregroundHistogram->real_time();
    else if ( m_foregroundHistogram && !m_secondaryHistogram && m_backgroundHistogram && m_foregroundHistogram->real_time() == m_backgroundHistogram->real_time() )
        return m_foregroundHistogram->real_time();
    else if ( m_foregroundHistogram && m_secondaryHistogram && !m_backgroundHistogram && m_foregroundHistogram->real_time() == m_secondaryHistogram->real_time() )
        return m_foregroundHistogram->live_time();
    else if ( !m_foregroundHistogram && m_secondaryHistogram && m_backgroundHistogram && m_backgroundHistogram->real_time() == m_secondaryHistogram->real_time() )
        return m_backgroundHistogram->real_time();
    
    else
        throw mup::ParserError( "Multiple spectra detected. Please specify a spectrum to use." );
}

float TerminalModel::gammaChannel( std::shared_ptr<const SpecUtils::Measurement> histogram, const double energy )
{
  auto cal = histogram ? histogram->energy_calibration() : nullptr;
  
  if( !cal || !cal->valid() )
    throw mup::ParserError( "No spectrum detected to gather corresponding gamma channel with energy." );
  
  return cal->channel_for_energy(energy);  //returns double
  //return histogram->find_gamma_channel( energy ); //returns int
}

// Gets the gamma channel count of a specific spectrum
float TerminalModel::gammaChannelContent( std::shared_ptr<const SpecUtils::Measurement> histogram, const double energy )
{
    if ( !histogram )
        throw mup::ParserError( "No spectrum detected to gather gamma channel count." );
    return histogram->gamma_channel_content( histogram->find_gamma_channel( energy ) );
}

// Gets the lower energy of a gamma channel of a specific spectrum
float TerminalModel::gammaChannelLowerEnergy( std::shared_ptr<const SpecUtils::Measurement> histogram, const double channel )
{
    if ( !histogram )
        throw mup::ParserError( "No spectrum detected to gather gamma channel lower energy." );
    return histogram->gamma_channel_lower( channel );
}

// Gets the central energy of a gamma channel of a specific spectrum
float TerminalModel::gammaChannelCentralEnergy( std::shared_ptr<const SpecUtils::Measurement> histogram, const double channel )
{
    if ( !histogram )
        throw mup::ParserError( "No spectrum detected to gather gamma channel central energy." );
    return histogram->gamma_channel_center( channel );
}

// Gets the upper energy of a gamma channel of a specific spectrum
float TerminalModel::gammaChannelHigherEnergy( std::shared_ptr<const SpecUtils::Measurement> histogram, const double channel )
{
    if ( !histogram )
        throw mup::ParserError( "No spectrum detected to gather gamma channel lower energy." );
    return histogram->gamma_channel_upper( channel );
}

float TerminalModel::gammaChannelWidth( std::shared_ptr<const SpecUtils::Measurement> histogram, const double channel )
{
    if ( !histogram )
        throw mup::ParserError( "No spectrum detected to gather gamma channel width." );
    return histogram->gamma_channel_width( channel );
}


float TerminalModel::gammaEnergyForChannel( std::shared_ptr<const SpecUtils::Measurement> histogram, const double channel )
{
  auto cal = histogram ? histogram->energy_calibration() : nullptr;
  if ( !cal || !cal->valid() )
    throw mup::ParserError( "No spectrum or invalid energy calibration to convert channel to energy." );
  
  return cal->energy_for_channel( channel );
}


// Gets the gamma integral of a specific spectrum
float TerminalModel::gammaIntegral( std::shared_ptr<const SpecUtils::Measurement> histogram, double energyLow, double energyHigh )
{
    if ( !histogram )
        throw mup::ParserError( "No spectrum detected to gather gamma channel lower energy." );
    
    if (energyLow > energyHigh)
        std::swap( energyLow, energyHigh );
    return histogram->gamma_integral( energyLow, energyHigh );
}

// Gets the gamma channel sum of a specific spectrum
float TerminalModel::gammaSum( std::shared_ptr<const SpecUtils::Measurement> histogram, const double startBin, const double endBin )
{
    if ( !histogram )
        throw mup::ParserError( "No spectrum detected to get gamma channel sum." );
  
  //TODO: gamma_channels_sum(...) take integer startBin and endBin - need to rectify and such
  
    return histogram->gamma_channels_sum( startBin, endBin );
}

double TerminalModel::numGammaChannels()
{
    updateHistograms();
    
    if (!m_foregroundHistogram && !m_secondaryHistogram && !m_backgroundHistogram)
        throw mup::ParserError( "No spectrum detected. Please add a spectrum to use." );
    
    if (m_foregroundHistogram && !m_secondaryHistogram && !m_backgroundHistogram)
        return m_foregroundHistogram->num_gamma_channels();
    else if (!m_foregroundHistogram && m_secondaryHistogram && !m_backgroundHistogram)
        return m_secondaryHistogram->num_gamma_channels();
    else if (!m_foregroundHistogram && !m_secondaryHistogram && m_backgroundHistogram)
        return m_backgroundHistogram->num_gamma_channels();
    
    else if ( m_foregroundHistogram && m_secondaryHistogram && m_backgroundHistogram &&
             m_foregroundHistogram->num_gamma_channels() == m_secondaryHistogram->num_gamma_channels() &&
             m_secondaryHistogram->num_gamma_channels() == m_backgroundHistogram->num_gamma_channels() )
        return m_foregroundHistogram->num_gamma_channels();
    else if ( m_foregroundHistogram && !m_secondaryHistogram && m_backgroundHistogram &&
             m_foregroundHistogram->num_gamma_channels() == m_backgroundHistogram->num_gamma_channels() )
        return m_foregroundHistogram->num_gamma_channels();
    else if ( m_foregroundHistogram && m_secondaryHistogram && !m_backgroundHistogram &&
             m_foregroundHistogram->num_gamma_channels() == m_secondaryHistogram->num_gamma_channels() )
        return m_foregroundHistogram->num_gamma_channels();
    else if ( !m_foregroundHistogram && m_secondaryHistogram && m_backgroundHistogram &&
             m_backgroundHistogram->num_gamma_channels() == m_secondaryHistogram->num_gamma_channels() )
        return m_backgroundHistogram->num_gamma_channels();
    
    else
        throw mup::ParserError( "Multiple spectra detected. Please specify a spectrum to use." );
}

// Automatically gets the gamma channel count for a specific spectrum, throws error if spectrum undetected
double TerminalModel::numGammaChannelsFor( const std::string& histogram )
{
    const std::string& hist (histogram);
    updateHistograms();
    
    if ( std::regex_match( hist, std::regex( "^(\\s*(foreground|fg)\\s*)$", std::regex::icase ) ) ) { // foreground
        if ( !m_foregroundHistogram ) throw mup::ParserError( "Foreground not detected." );
        else return m_foregroundHistogram->num_gamma_channels();
        
    } else if ( std::regex_match( hist, std::regex( "^(\\s*(secondary|sfg|secondaryforeground|secondforeground)\\s*)$", std::regex::icase ) ) ) { // second foreground
        if ( !m_secondaryHistogram ) throw mup::ParserError( "Secondary foreground not detected." );
        else return m_secondaryHistogram->num_gamma_channels();
        
    } else if ( std::regex_match( hist, std::regex( "^(\\s*(background|bg|background|back)\\s*)$", std::regex::icase ) ) ) { // background
        if ( !m_backgroundHistogram ) throw mup::ParserError( "Background not detected." );
        else return m_backgroundHistogram->num_gamma_channels();
    }
    
    throw mup::ParserError ( "Invalid argument for function 'numGammaChannelsFor( spectrum )'" );
}

double TerminalModel::gammaChannelAt( const double energy ) { return gammaFunctionOneArg( energy, &TerminalModel::gammaChannel ); }

double TerminalModel::gammaChannelFor( const std::string& histogram, const double energy )
{ return gammaFunctionOneArgFor(histogram, energy, &TerminalModel::gammaChannel ); }

// Automatically gets the gamma channel count of a spectrum, throws error if multiple or no spectra detected
double TerminalModel::gammaChannelCountAt ( const double energy ) { return gammaFunctionOneArg( energy, &TerminalModel::gammaChannelContent ); }
double TerminalModel::gammaChannelCountFor( const std::string& histogram, const double energy )
{ return gammaFunctionOneArgFor(histogram, energy, &TerminalModel::gammaChannelContent ); }

// Automatically gets the lower energy of a gamma channel of a spectrum, throws error if multiple or no spectra detected
double TerminalModel::gammaChannelLowerEnergyAt( const double energy ) { return gammaFunctionOneArg( energy, &TerminalModel::gammaChannelLowerEnergy ); }

// Automatically gets the lower energy for a gamma channel for a specific spectrum, throws error if spectrum undetected
double TerminalModel::gammaChannelLowerEnergyFor( const std::string& histogram, const double energy )
{ return gammaFunctionOneArgFor(histogram, energy, &TerminalModel::gammaChannelLowerEnergy ); }

double TerminalModel::gammaChannelCentralEnergyAt( const double energy ) { return gammaFunctionOneArg( energy, &TerminalModel::gammaChannelCentralEnergy ); }

// Automatically gets the lower energy for a gamma channel for a specific spectrum, throws error if spectrum undetected
double TerminalModel::gammaChannelCentralEnergyFor( const std::string& histogram, const double energy )
{ return gammaFunctionOneArgFor(histogram, energy, &TerminalModel::gammaChannelCentralEnergy ); }

// Automatically gets the higher energy of a gamma channel of a spectrum, throws error if multiple or no spectra detected
double TerminalModel::gammaChannelHigherEnergyAt( const double energy ) { return gammaFunctionOneArg( energy, &TerminalModel::gammaChannelHigherEnergy ); }

// Automatically gets the higher energy for a gamma channel for a specific spectrum, throws error if spectrum undetected
double TerminalModel::gammaChannelHigherEnergyFor( const std::string& histogram, const double energy )
{ return gammaFunctionOneArgFor(histogram, energy, &TerminalModel::gammaChannelHigherEnergy ); }

double TerminalModel::gammaChannelWidthAt( const double energy ) { return gammaFunctionOneArg( energy, &TerminalModel::gammaChannelWidth ); }

double TerminalModel::gammaEnergyForChannelAt( const double energy ) { return gammaFunctionOneArg( energy, &TerminalModel::gammaEnergyForChannel ); }


// Automatically gets the higher energy for a gamma channel for a specific spectrum, throws error if spectrum undetected
double TerminalModel::gammaChannelWidthFor( const std::string& histogram, const double energy )
{ return gammaFunctionOneArgFor(histogram, energy, &TerminalModel::gammaChannelWidth ); }

// Automatically gets the integral of a gamma channel of a spectrum, throws error if multiple or no spectra detected
double TerminalModel::gammaIntegralAt( const double energyLow, const double energyHigh )
{ return gammaFunctionTwoArg( energyLow, energyHigh, &TerminalModel::gammaIntegral ); }

// Automatically gets the integral for a gamma channel for a specific spectrum, throws error if spectrum undetected
double TerminalModel::gammaIntegralFor( const std::string& histogram, const double energyLow, const double energyHigh )
{ return gammaFunctionTwoArgFor( histogram, energyLow, energyHigh, &TerminalModel::gammaIntegral ); }

// Automatically gets the sum of a gamma channels of a spectrum, throws error if multiple or no spectra detected
double TerminalModel::gammaSumAt( const double startBin, const double endBin )
{ return gammaFunctionTwoArg( startBin, endBin, &TerminalModel::gammaSum ); }

// Automatically gets the sum of gamma channels for a specific spectrum, throws error if spectrum undetected
double TerminalModel::gammaSumFor( const std::string& histogram, const double startBin, const double endBin )
{ return gammaFunctionTwoArgFor( histogram, startBin, endBin, &TerminalModel::gammaIntegral ); }

double TerminalModel::gammaMin()
{
    updateHistograms();
    
    if (!m_foregroundHistogram && !m_secondaryHistogram && !m_backgroundHistogram)
        throw mup::ParserError( "No spectra detected. Please add a spectrum to use." );
    
    if (m_foregroundHistogram && !m_secondaryHistogram && !m_backgroundHistogram)
        return m_foregroundHistogram->gamma_energy_min();
    else if (!m_foregroundHistogram && m_secondaryHistogram && !m_backgroundHistogram)
        return m_secondaryHistogram->gamma_energy_min();
    else if (!m_foregroundHistogram && !m_secondaryHistogram && m_backgroundHistogram)
        return m_backgroundHistogram->gamma_energy_min();
    
    else if ( m_foregroundHistogram && m_secondaryHistogram && m_backgroundHistogram &&
             m_foregroundHistogram->gamma_energy_min() == m_secondaryHistogram->gamma_energy_min() &&
             m_secondaryHistogram->gamma_energy_min() == m_backgroundHistogram->gamma_energy_min() )
        return m_foregroundHistogram->num_gamma_channels();
    else if ( m_foregroundHistogram && !m_secondaryHistogram && m_backgroundHistogram &&
             m_foregroundHistogram->gamma_energy_min() == m_backgroundHistogram->gamma_energy_min() )
        return m_foregroundHistogram->gamma_energy_min();
    else if ( m_foregroundHistogram && m_secondaryHistogram && !m_backgroundHistogram &&
             m_foregroundHistogram->gamma_energy_min() == m_secondaryHistogram->gamma_energy_min() )
        return m_foregroundHistogram->gamma_energy_min();
    else if ( !m_foregroundHistogram && m_secondaryHistogram && m_backgroundHistogram &&
             m_backgroundHistogram->gamma_energy_min() == m_secondaryHistogram->gamma_energy_min() )
        return m_backgroundHistogram->gamma_energy_min();
    
    else
        throw mup::ParserError( "Multiple spectra detected. Please specify a spectrum to use." );
}

// Automatically gets the gamma energy minimum for a specific spectrum, throws error if spectrum undetected
double TerminalModel::gammaMinFor( const std::string& histogram )
{
    const std::string& hist (histogram);
    updateHistograms();
    
    if ( std::regex_match( hist, std::regex( "^(\\s*(foreground|fg)\\s*)$", std::regex::icase ) ) ) { // foreground
        if ( !m_foregroundHistogram ) throw mup::ParserError( "Foreground not detected." );
        else return m_foregroundHistogram->gamma_energy_min();
        
    } else if ( std::regex_match( hist, std::regex( "^(\\s*(secondary|sfg|secondaryforeground|secondforeground)\\s*)$", std::regex::icase ) ) ) { // second foreground
        if ( !m_secondaryHistogram ) throw mup::ParserError( "Secondary foreground not detected." );
        else return m_secondaryHistogram->gamma_energy_min();
        
    } else if ( std::regex_match( hist, std::regex( "^(\\s*(background|bg|background|back)\\s*)$", std::regex::icase ) ) ) { // background
        if ( !m_backgroundHistogram ) throw mup::ParserError( "Background not detected." );
        else return m_backgroundHistogram->gamma_energy_min();
    }
    
    throw mup::ParserError ( "Invalid argument for function 'gammaMinFor( spectrum )'" );
}

double TerminalModel::gammaMax()
{
    updateHistograms();
    
    if (!m_foregroundHistogram && !m_secondaryHistogram && !m_backgroundHistogram)
        throw mup::ParserError( "No spectra detected. Please add a spectrum to use." );
    
    if (m_foregroundHistogram && !m_secondaryHistogram && !m_backgroundHistogram)
        return m_foregroundHistogram->gamma_energy_max();
    else if (!m_foregroundHistogram && m_secondaryHistogram && !m_backgroundHistogram)
        return m_secondaryHistogram->gamma_energy_max();
    else if (!m_foregroundHistogram && !m_secondaryHistogram && m_backgroundHistogram)
        return m_backgroundHistogram->gamma_energy_max();
    
    else if ( m_foregroundHistogram && m_secondaryHistogram && m_backgroundHistogram &&
             m_foregroundHistogram->gamma_energy_max() == m_secondaryHistogram->gamma_energy_max() &&
             m_secondaryHistogram->gamma_energy_max() == m_backgroundHistogram->gamma_energy_max() )
        return m_foregroundHistogram->gamma_energy_max();
    else if ( m_foregroundHistogram && !m_secondaryHistogram && m_backgroundHistogram &&
             m_foregroundHistogram->gamma_energy_max() == m_backgroundHistogram->gamma_energy_max() )
        return m_foregroundHistogram->gamma_energy_max();
    else if ( m_foregroundHistogram && m_secondaryHistogram && !m_backgroundHistogram &&
             m_foregroundHistogram->gamma_energy_max() == m_secondaryHistogram->gamma_energy_max() )
        return m_foregroundHistogram->gamma_energy_max();
    else if ( !m_foregroundHistogram && m_secondaryHistogram && m_backgroundHistogram &&
             m_backgroundHistogram->gamma_energy_max() == m_secondaryHistogram->gamma_energy_max() )
        return m_backgroundHistogram->gamma_energy_max();
    
    else
        throw mup::ParserError( "Multiple spectra detected. Please specify a spectrum to use." );
}

// Automatically gets the gamma energy minimum for a specific spectrum, throws error if spectrum undetected
double TerminalModel::gammaMaxFor( const std::string& histogram )
{
    const std::string& hist (histogram);
    updateHistograms();
    
    if ( std::regex_match( hist, std::regex( "^(\\s*(foreground|fg)\\s*)$", std::regex::icase ) ) ) { // foreground
        if ( !m_foregroundHistogram ) throw mup::ParserError( "Foreground not detected." );
        else return m_foregroundHistogram->gamma_energy_max();
        
    } else if ( std::regex_match( hist, std::regex( "^(\\s*(secondary|sfg|secondaryforeground|secondforeground)\\s*)$", std::regex::icase ) ) ) { // second foreground
        if ( !m_secondaryHistogram ) throw mup::ParserError( "Secondary foreground not detected." );
        else return m_secondaryHistogram->gamma_energy_max();
        
    } else if ( std::regex_match( hist, std::regex( "^(\\s*(background|bg|background|back)\\s*)$", std::regex::icase ) ) ) { // background
        if ( !m_backgroundHistogram ) throw mup::ParserError( "Background not detected." );
        else return m_backgroundHistogram->gamma_energy_max();
    }
    
    throw mup::ParserError ( "Invalid argument for function 'gammaMaxFor( spectrum )'" );
}



double TerminalModel::drfFWHM( const double energy )
{
  std::shared_ptr<SpecMeas> foreground = m_viewer->measurment(SpecUtils::SpectrumType::Foreground);
  if( !foreground )
    throw mup::ParserError( "Foreground not foreground loaded." );
  
  std::shared_ptr<const DetectorPeakResponse> det = foreground->detector();
  if( !det )
    throw mup::ParserError( "No detector response loaded." );
  
  if( !det->isValid() )
    throw mup::ParserError( "Detector response loaded is not valid." );
  
  if( !det->hasResolutionInfo() )
    throw mup::ParserError( "Detector response does not contain peak-resolution information." );
  
  try
  {
    return det->peakResolutionFWHM( static_cast<float>(energy) );
  }catch( std::exception &e )
  {
    throw mup::ParserError( "Error getting peak FWHM: " + std::string(e.what()) );
  }
  
  return 0.0;
}//drfFWHM(energy)


double TerminalModel::drfIntrinsicEff( const double energy )
{
  std::shared_ptr<SpecMeas> foreground = m_viewer->measurment(SpecUtils::SpectrumType::Foreground);
  if( !foreground )
    throw mup::ParserError( "Foreground not foreground loaded." );
  
  std::shared_ptr<const DetectorPeakResponse> det = foreground->detector();
  if( !det )
    throw mup::ParserError( "No detector response loaded." );
  
  if( !det->isValid() )
    throw mup::ParserError( "Detector response loaded is not valid." );
  
  try
  {
    return det->intrinsicEfficiency( static_cast<float>(energy) );
  }catch( std::exception &e )
  {
    throw mup::ParserError( "Error getting intrinsic efficiency: " + std::string(e.what()) );
  }
  
  return 0.0;
}//drfIntrinsicEff(energy)


double TerminalModel::drfGeometricEff( const std::string &distance_str )
{
  std::shared_ptr<SpecMeas> foreground = m_viewer->measurment(SpecUtils::SpectrumType::Foreground);
  if( !foreground )
    throw mup::ParserError( "Foreground not foreground loaded." );
  
  std::shared_ptr<const DetectorPeakResponse> det = foreground->detector();
  if( !det || !det->isValid() )
    throw mup::ParserError( "No valid detector response loaded." );
  
  double distance = 0;
  try
  {
    distance = PhysicalUnits::stringToDistance(distance_str);
  }catch(...)
  {
    throw mup::ParserError( "Could not convert '" + distance_str + "' to a distance" );
  }
  
  if( distance < 0.0 )
    throw mup::ParserError( "Distance must not be negative." );
  
  try
  {
    return DetectorPeakResponse::fractionalSolidAngle( det->detectorDiameter(), distance );
  }catch( std::exception &e )
  {
    throw mup::ParserError( "Error getting geometric efficiency: " + std::string(e.what()) );
  }
  
  return 0.0;
}//drfGeometricEff(energy)


double TerminalModel::drfEfficiency( const double energy, const std::string &distance_str )
{
  std::shared_ptr<SpecMeas> foreground = m_viewer->measurment(SpecUtils::SpectrumType::Foreground);
  if( !foreground )
    throw mup::ParserError( "Foreground not foreground loaded." );
  
  std::shared_ptr<const DetectorPeakResponse> det = foreground->detector();
  if( !det || !det->isValid() )
    throw mup::ParserError( "No valid detector response loaded." );
  
  const bool fixed_geom = det->isFixedGeometry();
  
  double distance = 0;
  try
  {
    distance = PhysicalUnits::stringToDistance(distance_str);
  }catch(...)
  {
    if( !fixed_geom )
      throw mup::ParserError( "Could not convert '" + distance_str + "' to a distance" );
  }
  
  if( fixed_geom && (distance > 0.0) )
    throw mup::ParserError( "DRF is for fixed geometry, so you must specify a distance of 0 or negative" );
  
  if( distance < 0.0 )
    throw mup::ParserError( "Distance must not be negative." );
  
  try
  {
    return det->efficiency( energy, distance );
  }catch( std::exception &e )
  {
    throw mup::ParserError( "Error getting geometric efficiency: " + std::string(e.what()) );
  }
  
  return 0.0;
}//double drfEfficiency( const double energy, const std::string &distance )



float TerminalModel::gammaFunctionOneArg( const double energy, float (TerminalModel::*func)(std::shared_ptr<const SpecUtils::Measurement> histogram, const double arg) )
{
    updateHistograms();
    
    if (!m_foregroundHistogram && !m_secondaryHistogram && !m_backgroundHistogram)
        throw mup::ParserError( "No spectra detected. Please add a spectrum to use." );
    
    const double energy_value ( energy );
    
    if (m_foregroundHistogram && !m_secondaryHistogram && !m_backgroundHistogram)
        return (this->*func)( m_foregroundHistogram, energy_value );
    else if (!m_foregroundHistogram && m_secondaryHistogram && !m_backgroundHistogram)
        return (this->*func)( m_secondaryHistogram, energy_value );
    else if (!m_foregroundHistogram && !m_secondaryHistogram && m_backgroundHistogram)
        return (this->*func)( m_backgroundHistogram, energy_value );
    
    else if ( m_foregroundHistogram && m_secondaryHistogram && m_backgroundHistogram &&
             (this->*func)( m_foregroundHistogram, energy_value ) == (this->*func)( m_secondaryHistogram, energy_value ) &&
             (this->*func)( m_secondaryHistogram, energy_value ) == (this->*func)( m_backgroundHistogram, energy_value ) )
        return (this->*func)( m_foregroundHistogram, energy_value );
    else if ( m_foregroundHistogram && !m_secondaryHistogram && m_backgroundHistogram &&
             (this->*func)( m_foregroundHistogram, energy_value ) == (this->*func)( m_backgroundHistogram, energy_value ) )
        return (this->*func)( m_foregroundHistogram, energy_value );
    else if ( m_foregroundHistogram && m_secondaryHistogram && !m_backgroundHistogram &&
             (this->*func)( m_foregroundHistogram, energy_value ) == (this->*func)( m_secondaryHistogram, energy_value ) )
        return (this->*func)( m_foregroundHistogram, energy_value );
    else if ( !m_foregroundHistogram && m_secondaryHistogram && m_backgroundHistogram &&
             (this->*func)( m_backgroundHistogram, energy_value ) == (this->*func)( m_secondaryHistogram, energy_value ) )
        return (this->*func)( m_backgroundHistogram, energy_value );
    
    else
        throw mup::ParserError( "Multiple spectra detected. Please specify a spectrum to use." );
}

float TerminalModel::gammaFunctionTwoArg( const double arg1, const double arg2, float (TerminalModel::*func)(std::shared_ptr<const SpecUtils::Measurement> histogram,
                                                                                                              const double arg1, const double arg2 ) )
{
    updateHistograms();
    
    if (!m_foregroundHistogram && !m_secondaryHistogram && !m_backgroundHistogram)
        throw mup::ParserError( "No spectra detected. Please add a spectrum to use." );
    
    const double argument1 ( arg1 );
    const double argument2 ( arg2 );
    
    if (m_foregroundHistogram && !m_secondaryHistogram && !m_backgroundHistogram)
        return (this->*func)( m_foregroundHistogram, argument1, argument2 );
    else if (!m_foregroundHistogram && m_secondaryHistogram && !m_backgroundHistogram)
        return (this->*func)( m_secondaryHistogram, argument1, argument2 );
    else if (!m_foregroundHistogram && !m_secondaryHistogram && m_backgroundHistogram)
        return (this->*func)( m_backgroundHistogram, argument1, argument2 );
    
    else if ( m_foregroundHistogram && m_secondaryHistogram && m_backgroundHistogram &&
             (this->*func)( m_foregroundHistogram, argument1, argument2 ) == (this->*func)( m_secondaryHistogram, argument1, argument2 ) &&
             (this->*func)( m_secondaryHistogram, argument1, argument2 ) == (this->*func)( m_backgroundHistogram, argument1, argument2 ) )
        return (this->*func)( m_foregroundHistogram, argument1, argument2 );
    else if ( m_foregroundHistogram && !m_secondaryHistogram && m_backgroundHistogram &&
             (this->*func)( m_foregroundHistogram, argument1, argument2 ) == (this->*func)( m_backgroundHistogram, argument1, argument2 ) )
        return (this->*func)( m_foregroundHistogram, argument1, argument2 );
    else if ( m_foregroundHistogram && m_secondaryHistogram && !m_backgroundHistogram &&
             (this->*func)( m_foregroundHistogram, argument1, argument2 ) == (this->*func)( m_secondaryHistogram, argument1, argument2 ) )
        return (this->*func)( m_foregroundHistogram, argument1, argument2 );
    else if ( !m_foregroundHistogram && m_secondaryHistogram && m_backgroundHistogram &&
             (this->*func)( m_backgroundHistogram, argument1, argument2 ) == (this->*func)( m_secondaryHistogram, argument1, argument2 ) )
        return (this->*func)( m_backgroundHistogram, argument1, argument2 );
    
    else
        throw mup::ParserError( "Multiple spectra detected. Please specify a spectrum to use." );
}

float TerminalModel::gammaFunctionOneArgFor( const std::string& histogram, const double energy,
                                            float (TerminalModel::*func)(std::shared_ptr<const SpecUtils::Measurement> histogram, const double arg) )
{
    const std::string& hist (histogram);
    const double energy_value ( energy );
    
    updateHistograms();
    
    if ( std::regex_match( hist, std::regex( "^(\\s*(foreground|fg)\\s*)$", std::regex::icase ) ) ) { // foreground
        if ( !m_foregroundHistogram ) throw mup::ParserError( "Foreground not detected." );
        else return (this->*func)( m_foregroundHistogram, energy_value );
        
    } else if ( std::regex_match( hist, std::regex( "^(\\s*(secondary|sfg|secondaryforeground|secondforeground)\\s*)$", std::regex::icase ) ) ) { // second foreground
        if ( !m_secondaryHistogram ) throw mup::ParserError( "Secondary foreground not detected." );
        else return (this->*func)( m_secondaryHistogram, energy_value );
        
    } else if ( std::regex_match( hist, std::regex( "^(\\s*(background|bg|background|back)\\s*)$", std::regex::icase ) ) ) { // background
        if ( !m_backgroundHistogram ) throw mup::ParserError( "Background not detected." );
        else return (this->*func)( m_backgroundHistogram, energy_value );
        
    }
    throw mup::ParserError ( "Invalid argument for function." );
}

float TerminalModel::gammaFunctionTwoArgFor( const std::string& histogram, const double arg1, const double arg2,
                             float (TerminalModel::*func)(std::shared_ptr<const SpecUtils::Measurement> histogram, const double arg1, const double arg2) )
{
    const std::string& hist (histogram);
    const double arg_1 ( arg1 );
    const double arg_2 ( arg2 );
    
    updateHistograms();
    
    if ( std::regex_match( hist, std::regex( "^(\\s*(foreground|fg)\\s*)$", std::regex::icase ) ) ) { // foreground
        if ( !m_foregroundHistogram ) throw mup::ParserError( "Foreground not detected." );
        else return (this->*func)( m_foregroundHistogram, arg_1, arg_2 );
        
    } else if ( std::regex_match( hist, std::regex( "^(\\s*(secondary|sfg|secondaryforeground|secondforeground)\\s*)$", std::regex::icase ) ) ) { // second foreground
        if ( !m_secondaryHistogram ) throw mup::ParserError( "Secondary foreground not detected." );
        else return (this->*func)( m_secondaryHistogram, arg_1, arg_2 );
        
    } else if ( std::regex_match( hist, std::regex( "^(\\s*(background|bg|background|back)\\s*)$", std::regex::icase ) ) ) { // background
        if ( !m_backgroundHistogram ) throw mup::ParserError( "Background not detected." );
        else return (this->*func)( m_backgroundHistogram, arg_1, arg_2 );
        
    }
    throw mup::ParserError ( "Invalid argument for function." );
}


// Command Methods
// Return a regex object of the current regex-string for checking commands
std::string TerminalModel::commandRegex() { return commandRegexStr;  }

void TerminalModel::addCommand(const std::string& command, CommandType type)
{
  commandRegexStr.insert(2, command + "|");              // insert new command into the command-checking regex
  m_commandMap.insert ( CommandPair (command, type) );     // add new command to the map of known commands
    
  switch( m_commandMap.at(command) )
  {                                  // add to helper list of commands
    case VarmapCommand:
      addDropDownListItem( "showDefinedVariables( optional_variable )",
                           "defined vars map varmap",
                           "Output <b>information about defined variables</b>."
                           " Can specify a <i>single</i> variable, if given as an argument,"
                           " or output <b>all variable values</b> if no argument is specified." );
    break;
      
    case SearchPeakCommand:
      addDropDownListItem( "searchPeak( energy )",
                           "addpeak add peaks searchfor",
                           "Searches for and <b>adds peak</b> onto the spectrum at the"
                           " specified <i>energy</i> value." );
    break;
    
    case DeletePeakCommand:
      addDropDownListItem( "deletePeak( energy )",
                           "removepeak remove peaks deletepeaks",
                           "Searches for and <b>deletes peak</b> at the specified <i>energy</i>"
                           " value in the spectrum." );
    break;
    
    case RefitPeakCommand:
      addDropDownListItem( "refitPeak( energy )",
                           "refitpeaks refitting",
                           "Searches for and <b>refits peak</b> at the specified <i>energy</i>"
                           " value in the spectrum." );
    break;
    
    case DarkenCommand:
      addDropDownListItem( "darken()",
                           "command line cli black interface",
                           "<b>Darkens</b> the text areas for the Terminal and"
                           " <b><font color='white'>whitens</font></b> the text font." );
    break;
    
    case LightenCommand:
      addDropDownListItem( "lighten()",
                           "command line cli white interface",
                           "<b><font color='white'>Whitens</font></b> the text areas for the"
                           " Terminal and <b>darkens</b> the font. Currently set to default." );
    break;
    
    case SetNuclideCommand:
      addDropDownListItem( "setNuclide( nuclide, optional_age )",
                           "nuclides age reference photopeaks",
                           "Sets the <i>nuclide</i> and <i>age</i> inside the"
                           " ‘Reference Photopeaks Tab’. If the <i>age</i> is not specified,"
                           " then the default age is provided. If the nuclide could not be found,"
                           " then nothing happens." );
    break;
    
    case ClearVarCommand:
      addDropDownListItem( "clearVariable( variable )",
                           "variables delete map vars variablemap varmap",
                           "Deletes a stored <i>variable</i>. If the <i>variable</i>"
                           " does not exist, then an <b><font color='red'>error message</font></b>"
                           " is returned." );
    break;
    
    case SetEnergyRangeCommand:
      addDropDownListItem( "setRange( lower_energy, upper_energy )",
                           "x zooming energy range energyrange axis zoom",
                           "Sets the displayed x-axis range of the current spectrum."
                           " Outputs if action was a success or not." );
    break;
        
        
    case SetYaxisRangeCommand:
      addDropDownListItem( "setYRange( lower_y_counts, upper_y_counts )",
                           "y axis zooming counts range yrange axis zoom",
                           "Sets the displayed y-axis range of the current spectrum."
                           " Outputs if action was a success or not." );
    break;
      
    case SetLogYAxisMinCommand:
      addDropDownListItem( "setLogYAxisMin( lower_y_counts )",
                          "y axis zooming counts range yrange axis zoom",
                           "Sets the minimum displayed y-axis value when there is zero counts." );
    break;
        
    default: break;
  }//switch( m_commandMap.at(command) )
  
}//void TerminalModel::addCommand(const std::string& command, CommandType type)

// Executes inputted command from user
std::string TerminalModel::doCommand(const std::string& input)
{
  std::smatch match;
  std::regex_match(input, match, std::regex(commandRegex()) );
    
    const std::string& command = match[1];
    const std::string& arguments = match[2];
    
    switch ( m_commandMap.at(command) ) {
        case VarmapCommand:          return variableMapStr( arguments );
        case SetEnergyRangeCommand:  return setEnergyRange( arguments );
        case SetYaxisRangeCommand:   return setYRange( arguments );
        case SetLogYAxisMinCommand:  return setLogYAxisMin( arguments );
        case SearchPeakCommand:      return searchforPeak ( arguments );
        case DeletePeakCommand:      return deletePeak    ( arguments );
        case RefitPeakCommand:       return refitPeak     ( arguments );
        case DarkenCommand:          return darken        ( arguments );
        case LightenCommand:         return lighten       ( arguments );
        case SetNuclideCommand:      return setNuclide    ( arguments );
        case SaveCommand:            return saveFile      ( arguments );
        case ClearVarCommand:        return clearVar      ( arguments );
        default: break;
    }
    return "Error: unsupported command";
}

std::string TerminalModel::saveFile( const std::string& arguments )
{
    // Takes in destination path and output type of file
    const std::regex argumentRegex = std::regex ( "^\\s*([\\s*\\S+\\s*]*)\\s*\\,\\s*([\\s*\\S+\\s*]*)\\s*$" );
    std::ostringstream os;
    
    if (!std::regex_match(arguments, argumentRegex)) {      // User provided invalid arguments, display proper error message
        os << "The arguments (" << arguments << ") are not valid arguments for the function saveFile( path, output_type )";
        return os.str();
    }
  std::smatch match;
  std::regex_match(arguments, match, argumentRegex);
    
    const std::string& path = match[1], outputType = match[2];
    
    return std::string("[:savefile:") + path + std::string(":") + outputType + std::string(":]");
}

// Command: darken - darkens the text input and historytxtdiv fields and whitens the text colors if arg set to 'true'
// Syntax: darken(bool)
// Return value: a string determining if the command was successful or not
std::string TerminalModel::darken( const std::string& arguments ) { return "[:darken:]"; }

// Command: lighten - lightens the text input and historytxtdiv fields and darkens the text colors if arg set to 'true'
// Syntax: lighten(bool)
// Return value: a string determining if the command was successful or not
std::string TerminalModel::lighten( const std::string& arguments ) { return "[:lighten:]"; }

// Command: varmap - prints out information about defined m_variables. Can specify a single variable, or all m_variables.
// Syntax: varmap() -->      returns string of all declared m_variables
//         varmap( var ) --> returns value of a single variable (in string format, so you can't use its result in an expression)
// Return value: a string-format of the variable map inside model
std::string TerminalModel::variableMapStr( const std::string& args )
{
    std::ostringstream os;
    std::string arguments (args);
    SpecUtils::trim(arguments);
    
    if (arguments.empty() || arguments.find_first_not_of(' ') == std::string::npos) {              // Print out all m_variables if the argument is empty
        os << "VariableMap = {";
        if (!m_variables.empty()) {
            auto i = m_variables.begin();
            os << i->first << ":" << (i->second)->GetFloat();                                                 // Format: VariableMap = { var1:val1, var2:val2, ... }
            
            for (++i; i != m_variables.end(); ++i)
                os << ", " << i->first << ":" << (i->second)->GetFloat();
        }
        os << "}";
        
    } else {                                                          // Provide value for a specified variable
      
      std::regex validVariableRegex( ns_validVariableRegexArg );
      
        if ( isVariable(arguments) )                                  os << "Current value for variable (" << arguments << ") = " << m_variables.at(arguments);
        else if ( std::regex_match(arguments, validVariableRegex) )   os << "\"" << arguments << "\" is not a declared variable in the map";    // User specified undeclared variable
        else                                                          os << "\"" << arguments << "\" is not a valid variable";   // User specified an invalid variable
    }

    return os.str();
}

std::string TerminalModel::clearVar( const std::string& args )
{
    std::ostringstream os;
    std::string arguments (args);
    SpecUtils::trim(arguments);
    
    if ( arguments.empty() )
        return "Error: Please enter a variable to be removed.";
      
    if ( isVariable(arguments) ) {
        os << "Removing variable (" << arguments << ") that had value (" << m_variables.at(arguments)->GetFloat() << ")";
        m_parser->RemoveVar( arguments );
        delete m_variables[ arguments ];
        m_variables.erase( arguments );
    }else if ( std::regex_match(arguments, std::regex(ns_validVariableRegexArg)) ) {
      os << "\"" << arguments << "\" is not a declared variable in the map";    // User specified undeclared variable
    } else {
      os << "\"" << arguments << "\" is not a valid variable";   // User specified an invalid variable
    }
    
    return os.str();
}

// Command: setrange - zooms in or out of a graph inside a specific range of keV (or any other used unit)
// Syntax: setrange( lower, upper )
// Return value: string stating whether zoom was successful or unsuccessful
std::string TerminalModel::setEnergyRange( const std::string& arguments )
{
    std::ostringstream os;                      // Regex matching arguments: must be two valid expressions separated by commas inside parentheses
    const std::regex argumentRegex = std::regex ( "^\\s*([\\s*\\S*\\s*]*)\\s*\\,\\s*([\\s*\\S*\\s*]*)\\s*$" );
    
    
    if (!std::regex_match(arguments, argumentRegex)) {      // User provided invalid arguments, display proper error message
        os << "The arguments (" << arguments << ") are not valid arguments for the function setRange( lowerBound, upperBound )";
        return os.str();
    }
    try {
      std::smatch match;
      std::regex_match(arguments, match, argumentRegex);
        
        float lower, upper;
        m_parser->SetExpr(match[1]);  lower = m_parser->Eval().GetFloat();
        m_parser->SetExpr(match[2]);  upper = m_parser->Eval().GetFloat();
        
        if (lower > upper) std::swap(lower, upper);
        if (  std::abs(upper-lower) <= 1 )
            throw mup::ParserError( "Invalid arguments for function (setRange). Lower and upper bound must not have difference less than or equal to one." );
        
        m_viewer->setDisplayedEnergyRange(lower, upper);
        
        
        os << "Now setting energy range with lower bound (" << lower << ") and upper bound (" << upper << ").";
        
    } catch ( const mup::ParserError& e ) {
        os << "Error code " << e.GetCode() << ": " << e.GetMsg();
        
    } catch ( const std::exception& e ) {
        os << "Error: " << e.what();
    }
    return os.str();
}


std::string TerminalModel::setYRange( const std::string& arguments )
{
  // Regex matching arguments: must be two valid expressions separated by commas inside parentheses
  std::ostringstream os;
  const std::regex argumentRegex = std::regex ( "^\\s*([\\s*\\S*\\s*]*)\\s*\\,\\s*([\\s*\\S*\\s*]*)\\s*$" );
  
  if (!std::regex_match(arguments, argumentRegex)) {
    // User provided invalid arguments, display proper error message
    os << "The arguments (" << arguments << ") are not valid arguments for the function setYRange( lower_counts, upper_counts )";
    return os.str();
  }
  
  try
  {
    std::smatch match;
    std::regex_match(arguments, match, argumentRegex);
    
    float lower, upper;
    m_parser->SetExpr(match[1]);
    lower = m_parser->Eval().GetFloat();
    
    m_parser->SetExpr(match[2]);
    upper = m_parser->Eval().GetFloat();
    
    if (lower > upper)
      std::swap(lower, upper);
    
    //const float diff = upper - lower;
    //if( diff <= 0.0 )
    //  throw mup::ParserError( "Invalid arguments for function (setYRange). Lower and upper bound"
    //                           " must not have difference less than or equal to 0.01." );
    
    const std::tuple<double,double,Wt::WString> result = m_viewer->setYAxisRange(lower, upper);
    
    if( std::get<2>(result).empty() )
      os << "Setting y-range to [" << lower << ", " << upper << "].";
    else
      os << "Setting y-range to [" << lower << ", " << upper << "] not fully be fulfilled; set to ["
         << std::get<0>(result) << ", " << std::get<1>(result) << "]: "
         << std::get<2>(result).toUTF8();
  } catch ( const mup::ParserError& e ) {
    os << "Error code " << e.GetCode() << ": " << e.GetMsg();
    
  } catch ( const std::exception& e ) {
    os << "Error: " << e.what();
  }
  return os.str();
  
}//string TerminalModel::setYRange(string)


std::string TerminalModel::setLogYAxisMin( const std::string& arguments )
{
  try 
  {
    m_parser->SetExpr( arguments );
    const double counts = m_parser->Eval().GetFloat();
    
    if( counts <= 0 )
      throw std::runtime_error( "value must be larger than zero" );
    
    if( !m_viewer->setLogYAxisMin(counts) )
      return "Failed to set lower value.";
  }catch( const mup::ParserError &e )
  {
    return errorMessage(e);
  }catch( const std::exception &e )
  {
    return "Error: " + std::string(e.what());
  }
  
  return "";
}//string TerminalModel::setYRange(string)


// Command: searchpeak - does the work of double-clicking a peak inside the spectrum
// Syntax: searchpeak( x )
// Return value: none, basically highlights a peak on the spectrum and outputs whether or not it was successful or not  
std::string TerminalModel::searchforPeak( const std::string& arguments )
{
    std::ostringstream os;
    m_parser->SetExpr( arguments );
    
    try {
        const double energy = m_parser->Eval().GetFloat();
        
        PeakModel *peakModel = m_viewer->peakModel();
        
        if ( !peakModel )
            throw mup::ParserError( "Peak model is not set." );
        if ( !(peakModel->peaks()) )
            throw mup::ParserError( "No peaks detected for peak model." );
        if ( energy == 0 )
            throw mup::ParserError( "Unable to add peak at energy (0)." );      // throw domain error if user is adding peak at 0 keV
        
        
        const int previousPeakVecSize = static_cast<int>( peakModel->peaks()->size() );
    
        m_viewer->searchForSinglePeak( energy );
        const int currentPeakVecSize = static_cast<int>( peakModel->peaks()->size() );
        PeakModel::PeakShrdPtr peak = peakModel->nearestPeak( energy );
        
        if ( !peak || currentPeakVecSize <= previousPeakVecSize || energy == 0 )
            os << "Could not find a peak near " << energy << " keV";
        else
            os << "Successfully added peak near " << energy << " keV; actual mean value was at " << peak->mean() << " keV";
        
    } catch (const mup::ParserError& e) {
        return errorMessage(e);
        
    } catch ( const std::exception& e ) {
        os << "Error: " << e.what();
    }
    
    return os.str();
}

// Command: deletepeak - does the work of double-clicking a peak inside the spectrum
// Syntax: deletepeak( x )
// Return value: none, basically deletes the specified peak on the spectrum and outputs whether or not it was successful or not
std::string TerminalModel::deletePeak( const std::string& arguments )
{
    std::ostringstream os;
    m_parser->SetExpr( arguments );
    
    try {
        const double energy = m_parser->Eval().GetFloat();
        
        PeakModel *peakModel = m_viewer->peakModel();
        if (!peakModel)
            return "Error: Peak model is not set.";
        
        PeakModel::PeakShrdPtr peak = peakModel->nearestPeak( energy );
        
        if ( !peak )
            os << "Error: Could not delete peak at " << energy << " keV; no peaks available to delete";
        else if ( !energyIsWithinPeak( peak, energy ) )
            os << "Could not delete peak at " << energy << " keV; please try a different energy level";
        else {
            peakModel->removePeak( peak );
            os << "Successfully deleted peak near " << energy << " keV; " << "actual mean value was at " << peak->mean() << " keV";
        }
        
    } catch ( const mup::ParserError& e ) {
        return errorMessage(e);
        
    } catch ( const std::exception& e ) {
        os << "Error: " << e.what();
    }
    
    return os.str();
}

// Command: refitpeak - does the work of refitting the peak through mouse-clicking
// Syntax: refitpeak( x )
// Return value: none, basically refits the specified peak on the spectrum and outputs whether or not it was successful or not
std::string TerminalModel::refitPeak ( const std::string& argument )
{
    std::ostringstream os;
    m_parser->SetExpr( argument );
    
    try {
        const double energy = m_parser->Eval().GetFloat();
        
        PeakModel *peakModel = m_viewer->peakModel();
        PeakModel::PeakShrdPtr peak = peakModel->nearestPeak( energy );
        
        if ( !peak )
            os << "Error: Could not refit peak at " << energy << " keV; no peaks available to delete";
        else if ( !energyIsWithinPeak( peak, energy ) )
            os << "Could not refit peak at " << energy << " keV; no peaks detected near " << energy << " keV";
        else {
            os << "Peak refitting not yet implemented.";
//            os << "Successfully refitted peak near " << energy << " keV; " << "actual mean value was at " << peak->mean() << " keV";
        }
        
    } catch ( const mup::ParserError& e ) {
        return errorMessage(e);
        
    } catch ( const std::exception& e ) {
        os << "Error: " << e.what();
    }
    
    return os.str();
}

std::string TerminalModel::setNuclide( const std::string& arguments )
{
    ReferencePhotopeakDisplay *nuclideReference = m_viewer->referenceLinesWidget();
    const std::regex argumentRegex = std::regex ( "^(.+)\\,(.+)$" );
    const std::regex argumentRegexWithoutAge = std::regex ( "^([^,]+)$" );
    
    if ( !std::regex_match(arguments, argumentRegex) && !std::regex_match(arguments, argumentRegexWithoutAge) )
        return "Error: failed to set nuclude and age because of invalid arguments.";
    
  std::smatch match;
  std::regex_match(arguments, match, argumentRegex);
    
    const bool argumentContainsAge = std::regex_match(arguments, argumentRegex);
    
    const std::string& nuclide = argumentContainsAge ? match[1] : arguments,
                        age = argumentContainsAge ? match[2] : std::string();
    
    try { nuclideReference->setNuclideAndAge( nuclide, age ); }   // not sure if this method throws an error, so just being on the safe side
    catch ( const std::exception& e ) { return std::string("Error: ") + e.what(); }
    
    return std::string("Setting nuclide (") + SpecUtils::trim_copy(nuclide) + ") and age ("
           + (argumentContainsAge ? SpecUtils::trim_copy(age) : std::string("default")) + ")";
}

// Command: peak_area - does the work of double-clicking a peak inside the spectrum
// Syntax: peak_area( x )
// Return value: area of the specified peak, 0 if no peak was found
double TerminalModel::peakArea( const double argument )
{
    const double energy = argument;
    
    PeakModel *peakModel = m_viewer->peakModel();
    PeakModel::PeakShrdPtr peak = peakModel->nearestPeak( energy );
    
    return peak && energyIsWithinPeak( peak, energy ) ? peak->peakArea() : 0;
}

// Command: peak_mean - gets the mean value of a peak
// Syntax: peak_mean( x )
// Return value: mean of the specified peak, 0 if no peak was found
double TerminalModel::peakMean( const double argument )
{
    const double energy = argument;
    
    PeakModel *peakModel = m_viewer->peakModel();
    PeakModel::PeakShrdPtr peak = peakModel->nearestPeak( energy );
    
    return peak && energyIsWithinPeak( peak, energy ) ? peak->mean() : 0;
}

// Command: peak_fwhm - gets the fwhm value of a peak
// Syntax: peak_fwhm( x )
// Return value: fwhm of the specified peak, 0 if no peak was found
double TerminalModel::peakFwhm( const double argument )
{
    const double energy = argument;
    
    PeakModel *peakModel = m_viewer->peakModel();
    PeakModel::PeakShrdPtr peak = peakModel->nearestPeak( energy );
    
    return peak && energyIsWithinPeak( peak, energy ) ? peak->fwhm() : 0;
}

// Command: peak_amp - gets the amplitude value of a peak
// Syntax: peak_amp( x )
// Return value: amplitude of the specified peak, 0 if no peak was found
double TerminalModel::peakAmp( const double argument )
{
    const double energy = argument;
    
    PeakModel *peakModel = m_viewer->peakModel();
    PeakModel::PeakShrdPtr peak = peakModel->nearestPeak( energy );
    
    return peak && energyIsWithinPeak( peak, energy ) ? peak->amplitude() : 0;
}

// Command: peak_sigma - gets the sigma value of a peak
// Syntax: peak_sigma( x )
// Return value: sigma value of the specified peak, 0 if no peak was found, throws error if peak is not Gaussian
double TerminalModel::peakSigma( const double argument )
{
    const double energy = argument;
    
    PeakModel *peakModel = m_viewer->peakModel();
    PeakModel::PeakShrdPtr peak = peakModel->nearestPeak( energy );
    
    if ( peak && energyIsWithinPeak( peak, energy ) ) {
        if ( !peak->gausPeak() )
            throw mup::ParserError( "Peak type must be Gaussian" );
        return peak->sigma();
    }
    return 0;
}

// Command: peak_chi2 - gets the chi2dof value of a peak
// Syntax: peak_chi2( x )
// Return value: chi2dof of the specified peak, 0 if no peak was found
double TerminalModel::peakChi2dof ( const double argument )
{
    const double energy = argument;
    
    PeakModel *peakModel = m_viewer->peakModel();
    PeakModel::PeakShrdPtr peak = peakModel->nearestPeak( energy );
    
    return peak && energyIsWithinPeak( peak, energy ) ? peak->chi2dof() : 0;
}

// Command: peak_gauss - gets the gaussian integral of a peak
// Syntax: peak_gauss( energy, x0, x1 )
// Return value: gaussian integral of the specified peak, 0 if no peak was found
double TerminalModel::peakGaussianIntegral ( const double argument, const double x0, const double x1 )
{
    const double energy = argument;
    
    PeakModel *peakModel = m_viewer->peakModel();
    PeakModel::PeakShrdPtr peak = peakModel->nearestPeak( energy );
    
    return peak && energyIsWithinPeak( peak, energy ) ? peak->gauss_integral( x0, x1 ) : 0;
}






// Variable Methods
// Assigns variable into the map, proper error handling done
// How it's done: capture variable on left side of '=' sign
//  Evaluate expression on right side of '=' sign
//      Set the variable value to that evaluated expression
std::string TerminalModel::assignVariable(const std::string& input)
{
  std::smatch match;
    if( !std::regex_match(input, match, std::regex(ns_variableAssignmentRegexArg)) ){
      return "No variable to assign";
    }
    
    std::ostringstream os;
    
    const std::string& variable = match[1];
    if ( std::regex_match(variable, std::regex(ns_keywordRegexArg) ) ) {                                              // invalid variable error (built-in error)
        const mup::ParserError e = mup::ParserError( "Invalid name. Variable must not be a keyword.");
        return errorMessage(e);
    }
    
    double value = 0;
    const std::string& expression = match[2];
    
    m_parser->SetExpr(expression);
    try { value = double( m_parser->Eval().GetFloat() ); }                               // catch error on right side of variable assignment
    catch ( const mup::ParserError& e ) {
        os << "Error code " << e.GetCode() << ": " << e.GetMsg();
        return os.str();
        
    } catch ( const std::exception& e ) {
        return std::string("Error: ") + std::string(e.what());
    }
    
    if ( m_parser->IsVarDefined(variable) ) {    // if variable is already in Variable Map
        
        // if re-assigning variable with same value, don't assign again
        if (m_variables.find(variable) != m_variables.end() && *(m_variables.at(variable)) == mup::Value(value) ) {
            os << "Variable " << variable << " already initialized with value(" << value << ")";
            return os.str();
        }
        
        // replace old value in variable map with new value
        const VariableMap::iterator& i = m_variables.find(variable);
        os << "Re-assigned variable " << variable << " from old value(" << *(i->second) << ") to new value(" << value << ")";
        
//        mup::Value* to_delete = i->second;
        *m_variables.at( variable ) = mup::Value(value);
//        delete to_delete;
//        m_parser->DefineVar(variable, mup::Variable( mup::Value( i->second ) ) );
        
        return os.str();
    }
    
    // add new variable inside the variable map
    m_variables.insert ( std::pair<std::string,mup::Value*> (variable, new mup::Value(value) ) );
    m_parser->DefineVar(variable, mup::Variable( m_variables.at( variable ) ) );
    os << "Assigned variable " << variable << " to value(" << value << ")";
    
    return os.str();
}


// Can use this to allow implicit variable declarations inside calculator
// Eg.  x is not defined yet by the user
//      User inputs: 'x + 2'    ----> result is 2, x implicitly defined and equals 0
double* TerminalModel::addImplicitVariable(const char* variable, void* pUserData)
{
    static double valBuff[100];     // array of m_variables held inside class
    static int iVal = -1;           // amount of m_variables allowed
    
    std::cout << " Adding variable: " << variable << " \\ (slots left: " << 99-(++iVal) << ")" << std::endl;
    valBuff[iVal] = 0;
//    m_variables.insert ( std::pair<std::string,mup::Value*> ( variable, new mup::Value(0) ) );
    
    if (iVal >= 99)
        throw mup::ParserError("Variable buffer overflow.");
    
    return &valBuff[iVal];
}



// - Helper Methods
// Determines if a given energy is within a certain range in or near a peak.
bool TerminalModel::energyIsWithinPeak( PeakModel::PeakShrdPtr peak, const double energy )
{
    const double lowerBound =  peak->mean() - ( 3 * (peak->mean() - peak->lowerX()) );
    const double upperBound = peak->mean() + ( 3 * (peak->upperX() - peak->mean()) );
    return lowerBound <= energy && energy <= upperBound;
}

// Determines if argument string is was declared and defined in the calculator
bool TerminalModel::isVariable(const std::string& inputVariable)
{
    if (m_variables.count(inputVariable) <= 0)
        return false;
    
    return true;
}
