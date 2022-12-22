//
//  TerminalModel.h
//  InterSpec
//
//  Created by Christian Kenneth Morte on 7/15/16.
//
//

#ifndef TerminalModel_h
#define TerminalModel_h
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

#include <map>
#include <array>
#include <tuple>
#include <string>
#include <vector>
#include <memory>
#include <algorithm>


#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"


//Forward declarations
class InterSpec;
namespace SpecUtils{ class Measurement; }

namespace mup
{
  class Value;
  class ParserX;
  class ICallback;
}


// Useful typedefs
typedef std::map<std::string, mup::Value*> VariableMap;
typedef std::map<std::string, int>  CommandMap;

class TerminalModel {
public:     // Main methods
    TerminalModel( InterSpec* viewer );
    ~TerminalModel();
    
    std::string evaluate( std::string input );
    void setViewer( InterSpec* spectViewer );
    
    typedef std::tuple<std::string,std::string,std::string> CommandHelperTuple;
    typedef std::vector<CommandHelperTuple> CommandHelperList;
    CommandHelperList commandsFunctionsList();
    
protected:  // Input Types
    enum InputType { Command, VariableAssignment, Operation };
    InputType inputType(const std::string& input);
    
    
public:  // Function methods (for parser, can combine with expressions)
    // To add a new function in the Terminal tool, please see the instructions underneath the TerminalModel class declaration.
    void addFunction( mup::ICallback* function, const std::string& tags="", const std::string& toolTip="" );
    
    double nuclideIntensity(  const double energy );
    double nuclideIntensityForParticle( std::string particle, const double energy );
    double nuclideEnergy( const double intensity );
    
    double liveTime( const std::string& argument );
    double realTime( const std::string& argument );
    double liveTimeWithoutArgument();
    double realTimeWithoutArgument();
    
    double peakArea     ( const double argument );
    double peakMean     ( const double argument );
    double peakSigma    ( const double argument );
    double peakFwhm     ( const double argument );
    double peakAmp      ( const double argument );
    double peakChi2dof  ( const double argument );
    
    double peakGaussianIntegral ( const double argument, const double x0, const double x1 );
    
    double numGammaChannels();
    double numGammaChannelsFor( const std::string& histogram );
    double gammaChannelAt( const double energy );
    double gammaChannelFor( const std::string& histogram, const double energy );
    double gammaChannelCountAt ( const double energy );
    double gammaChannelCountFor( const std::string& histogram, const double energy );
    double gammaChannelLowerEnergyFor( const std::string& histogram, const double channel );
    double gammaChannelLowerEnergyAt ( const double channel );
    double gammaChannelCentralEnergyFor( const std::string& histogram, const double channel );
    double gammaChannelCentralEnergyAt ( const double channel );
    double gammaChannelHigherEnergyFor( const std::string& histogram, const double channel );
    double gammaChannelHigherEnergyAt( const double channel );
    double gammaChannelWidthFor( const std::string& histogram, const double channel );
    double gammaChannelWidthAt( const double channel );
    double gammaEnergyForChannelAt( const double channel );
    double gammaIntegralAt( const double energyLow, const double energyHigh );
    double gammaIntegralFor( const std::string& histogram, const double energyLow, const double energyHigh );
    double gammaSumAt( const double startBin, const double endBin );
    double gammaSumFor( const std::string& histogram, const double startBin, const double endBin );
    double gammaMin();
    double gammaMinFor( const std::string& histogram );
    double gammaMax();
    double gammaMaxFor( const std::string& histogram );
  
  double drfFWHM( const double energy );
  double drfIntrinsicEff( const double energy );
  double drfGeometricEff( const std::string &distance );
  double drfEfficiency( const double energy, const std::string &distance );
  
protected:  /* Command methods (complete actions on Spectrum, cannot be used with parser)
             To add a new command in the Terminal tool:
                 1. Add the corresponding CommandType enum inside the TerminalModel Header file.
                 2. Call the 'addCommand( command, type )' method inside the main constructor of TerminalModel.
                        - command:  arg for what the command name should be (this is what the user calls)
                        - type:     arg for the CommandType enum of the command
                 3. Declare/define the main function your command will execute when it is called by the user.
                        - Ensure that the return value is a string. Commands have no return value; they execute commands inside InterSpec.
                 4. Assign the proper CommandType enum of your new command in the switch statement of TerminalModel's 'doCommand( input )'
                    with the function you just declared in step 3.
                 5. Assure that the command works properly inside the Terminal tool. If it doesn't, repeat the steps below to see if any mistakes
                    were made. Contact Christian Morte or Will for help or suggestions.
             
             Extra notes for commands: Methods for commands should NOT throw errors, instead return a proper an error message.
     */
    
    enum CommandType {
      VarmapCommand,
      ClearVarCommand,
      SetEnergyRangeCommand,
      SetYaxisRangeCommand,
      SearchPeakCommand,
      DeletePeakCommand,
      RefitPeakCommand,
      DarkenCommand,
      LightenCommand,
      SetNuclideCommand,
      SaveCommand
    };
    typedef std::pair<std::string, CommandType> CommandPair;
    std::string commandRegexStr = std::string("^(default*)\\(([\\s*\\S*\\s*]*)\\)$");  // initial command regex; this changes everytime new command is issued
    
    static double* addImplicitVariable(const char* variable, void* pUserData);
    std::string assignVariable(const std::string& input);
    
    // Commands (these have no return values; instead output strings of info after the command)
    std::string saveFile      ( const std::string& arguments );
    std::string darken        ( const std::string& arguments );
    std::string lighten       ( const std::string& arguments );
    std::string setEnergyRange( const std::string& arguments );
    std::string setYRange     ( const std::string& arguments );
    std::string searchforPeak ( const std::string& arguments );
    std::string deletePeak    ( const std::string& arguments );
    std::string refitPeak     ( const std::string& arguments );
    std::string setNuclide    ( const std::string& arguments );
    std::string variableMapStr( const std::string& arguments );
    std::string clearVar      ( const std::string& arguments );
    
    void addCommand(const std::string& command, CommandType type);
    std::string commandRegex();
    std::string doCommand(const std::string& input);
    
    
protected:  // Internal helper methods
    // Drop-down list helpers
    void addDropDownListHeader( const std::string& header );
    void addDropDownListSeparator();
    void addDropDownListItem( const std::string& command, const std::string& tags="", const std::string& tooltip="" );
    
    // Variable helpers
    bool isVariable(const std::string& inputVariable);
    
    // Non-built-in Parser Function Helpers
    void updateHistograms();
    
    double foregroundLiveTime();
    double secondaryForegroundLiveTime();
    double backgroundLiveTime();
    
    double foregroundRealTime();
    double secondaryForegroundRealTime();
    double backgroundRealTime();
    
    float gammaFunctionOneArg( const double energy, float (TerminalModel::*func)(std::shared_ptr<const SpecUtils::Measurement> histogram, const double arg) );
    float gammaFunctionTwoArg( const double arg1, const double arg2, float (TerminalModel::*func)(std::shared_ptr<const SpecUtils::Measurement> histogram, const double arg1, const double arg2 ) );
    
    float gammaFunctionOneArgFor( const std::string& histogram, const double arg,
                                        float (TerminalModel::*func)(std::shared_ptr<const SpecUtils::Measurement> histogram, const double arg) );
    float gammaFunctionTwoArgFor( const std::string& histogram, const double arg1, const double arg2,
                                        float (TerminalModel::*func)(std::shared_ptr<const SpecUtils::Measurement> histogram, const double arg1, const double arg2) );
    
    float gammaChannel( std::shared_ptr<const SpecUtils::Measurement> histogram, const double energy );
    float gammaChannelContent( std::shared_ptr<const SpecUtils::Measurement> histogram, const double energy );
    float gammaChannelLowerEnergy( std::shared_ptr<const SpecUtils::Measurement> histogram, const double channel );
    float gammaChannelCentralEnergy( std::shared_ptr<const SpecUtils::Measurement> histogram, const double channel );
    float gammaChannelHigherEnergy( std::shared_ptr<const SpecUtils::Measurement> histogram, const double channel );
    float gammaChannelWidth( std::shared_ptr<const SpecUtils::Measurement> histogram, const double channel );
    float gammaEnergyForChannel( std::shared_ptr<const SpecUtils::Measurement> histogram, const double channel );
    float gammaIntegral( std::shared_ptr<const SpecUtils::Measurement> histogram, const double energyLow, const double energyHigh );  // continous sum (energy)
    float gammaSum( std::shared_ptr<const SpecUtils::Measurement> histogram, const double startBin, const double endBin );            // discrete sum (bins)
    
    bool energyIsWithinPeak( PeakModel::PeakShrdPtr peak, const double energy );
    
private:
    // Terminal Model Members
    std::unique_ptr<mup::ParserX>   m_parser;
    CommandMap     m_commandMap;
    VariableMap    m_variables;
    
    InterSpec *m_viewer;
    std::shared_ptr<const SpecUtils::Measurement> m_foregroundHistogram;
    std::shared_ptr<const SpecUtils::Measurement> m_backgroundHistogram;
    std::shared_ptr<const SpecUtils::Measurement> m_secondaryHistogram;
    
    /* To add a new command/function into the Drop-Down helper list:
     FOR COMMANDS ONLY:
            1.  Inside the void TerminalModel::addCommand(const std::string& command, CommandType type) method, create a new case inside the 
                switch statement with your new command's CommandType enum value.
            2.  Call the TerminalModel's addDropDownListItem( item, tags, toolTip ) method with your desired function specifications passed in the 'item' argument. Make
                sure that the case block ends with a 'break' statement.
     
    FOR FUNCTIONS ONLY:
            1.  In order to provide a function callback description inside the drop down list, you MUST overload a function class' GetDesc() method.
                Return the desired callback description for the function inside this method.
            2.  For including search tags for your function, you must define a tags function inside the function class. This function, for consistency, should be defined as "std::string tags() const" and return a string value for all the search tags for that specific function. Tags are separated by spaces.
            3. For including a tool-tip for a function, you must define a tool-tip function inside the function class. This function, for consistency, should be defined as "std::string toolTip() const" and return a string value for the text inside the command's tool-tip.
            4.  Calling the void TerminalModel::addFunction( mup::ICallback* function, const std::string& tags ) will automatically add the callback description into the drop down list.
            
            For more information on creating functions and their tags, refer to the instructions comment above the namespace where all class functions and commands are located in the TerminalModel.cpp.
    */
    CommandHelperList m_dropDownList;
};

#endif /* TerminalModel_h */
