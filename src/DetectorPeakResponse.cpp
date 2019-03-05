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

#include <mutex>
#include <cmath>
#include <memory>
#include <cctype>
#include <string>
#include <vector>
#include <istream>
#include <fstream>
#include <sstream>
#include <utility>
#include <sstream>
#include <numeric>
#include <stdexcept>

#include <boost/functional/hash.hpp>

#include "rapidxml/rapidxml.hpp"

//Roots Minuit2 includes
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
//#include "Minuit2/Minuit2Minimizer.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/SimplexMinimizer.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinimize.h"

#include "mpParser.h"

#include "InterSpec/PeakModel.h"
#include "InterSpec/PhysicalUnits.h"
#include "SpecUtils/UtilityFunctions.h"
#include "SpecUtils/SpectrumDataStructs.h"
#include "InterSpec/DetectorPeakResponse.h"


using namespace std;
const int DetectorPeakResponse::sm_xmlSerializationVersion = 0;

namespace
{
  //calcA(...) is for use from DetectorPeakResponse::akimaInterpolate(...).
  //  This function does not to any error checking that the input is valid.
  float calcA( const size_t i,
              const std::vector<DetectorPeakResponse::EnergyEfficiencyPair> &xy )
  {
    const size_t n = xy.size();
    const float x_i  = xy[i].energy;
    const float x_m1 = xy[i-1].energy;
    const float x_m2 = xy[i-2].energy;
    const float x_p1 = xy[i+1].energy;
    const float x_p2 = xy[i+2].energy;
    const float y_i  = xy[i].efficiency;
    const float y_m1 = xy[i-1].efficiency;
    const float y_m2 = xy[i-2].efficiency;
    const float y_p1 = xy[i+1].efficiency;
    const float y_p2 = xy[i+2].efficiency;
    const float p_i = (y_p1 - y_i) / (x_p1 - x_i);
    
    float A_i;
    if( i == 0 )
    {
      A_i = p_i;
    }else if( i == 1 )
    {
      const float p_0 = (xy[1].efficiency - xy[0].efficiency) / (xy[1].energy - xy[0].energy);
      A_i = (p_0+p_i)/2.0f;
    }else if( i == (n-1) )
    {
      A_i = p_i;
    }else if( i == (n-2) )
    {
      const float p = (xy[n-2].efficiency - xy[n-3].efficiency) / (xy[n-2].energy - xy[n-3].energy);
      A_i = (p_i + p)/2.0f;
    }else
    {
      const float p_1 = (y_p2 - y_p1) / (x_p2 - x_p1);
      const float p_m1 = (y_i - y_m1) / (x_i - x_m1);
      const float p_m2 = (y_m1 - y_m2) / (x_m1 - x_m2);
      
      const float w1 = fabs( p_1 - p_i );
      const float w2 = fabs( p_m1 - p_m2 );
      if( (w1+w2) == 0.0f )
        A_i = (p_m1 + p_i) / 2.0f;
      else
        A_i = (w1*p_m1 + w2*p_i) / (w1 + w2);
    }
    return A_i;
  }//float calcA(...)
}//namespace


  
FormulaWrapper::FormulaWrapper( const std::string &fcnstr )
  : m_fcnstr( fcnstr ), m_var_name( "x" )
{
  UtilityFunctions::to_lower( m_fcnstr );
  UtilityFunctions::erase_any_character( m_fcnstr, "\t\n\r" );
  
  m_var_name = find_variable_name( m_fcnstr );
  
  m_parser.reset( new mup::ParserX() );
  m_value.reset( new mup::Value(100.0) );
  
  if( m_fcnstr.empty() )
    m_fcnstr = "1";
  
  try
  {
    m_parser->DefineVar( m_var_name, mup::Variable(m_value.get()) );
    m_parser->SetExpr( m_fcnstr );
  
    //Trigger an initial evaluation - jic.
    m_parser->Eval();
  }catch( mup::ParserError &e ) //if no Formula is set or in case of any other error related to the formula.
  {
    //stringstream msg;
    //msg << "Error evaluating expression \"" << e.GetExpr() << "\": " << e.GetMsg();
    
    std::string msg = e.GetMsg();
    msg = "<span style=\"color:black;font-weight:bold;\">" + msg + ": </span>";
    
    const int errorpos = e.GetPos();
    std::string preeqn, posteqn, errorstr;
    if( errorpos+1 < int(m_fcnstr.size()) && errorpos >= 0 )
      posteqn = m_fcnstr.substr(errorpos+1);
    if( posteqn.size() > 12 )
      posteqn = posteqn.substr(0,12) + "...";
    if( errorpos < int(m_fcnstr.size()) && errorpos >= 0 )
      errorstr = m_fcnstr.substr(errorpos,1);
    
    if( errorpos > 14 && errorpos < int(m_fcnstr.size()) )
      preeqn = "..." + m_fcnstr.substr(errorpos-11,11);
    else if( errorpos < int(m_fcnstr.size()) )
      preeqn = m_fcnstr.substr(0,errorpos);
    
    msg += "<span style=\"font-family:monospace;color:black;\">"
    + preeqn
    + "<span style=\"color:red;text-decoration:underline\">"
    + errorstr + "</span><span style=\"color:#2F4F4F;\">"
    + posteqn+ "</span></span>";
    
    throw std::runtime_error( msg );
    
    
    //throw std::runtime_error( "Error evaluating expression \"" + e.GetExpr() + "\": " + e.GetMsg() );
    
    //const string_type& GetExpr() const;
    //string_type GetMsg() const;
    //int GetPos() const;
    //const string_type& GetToken() const;
    //EErrorCodes GetCode() const;
    //ErrorContext& GetContext();
    
  }catch( std::exception &e )
  {
    throw std::runtime_error( "Error parsing expression: " + string(e.what()) );
  }//try / catch
  
}//FormulaWrapper
  
FormulaWrapper::~FormulaWrapper()
{
}

std::string FormulaWrapper::find_variable_name( std::string eqn )
{
  //Make case-insensitive
  UtilityFunctions::to_lower( eqn );
  
  //Get rid of whitespaces
  //eqn.erase( std::remove_if(eqn.begin(), eqn.end(), boost::algorithm::is_any_of(" \t\n\r") ), eqn.end());
  
  auto is_whitespace_fcn = [](const char &c) -> bool {
    return (c==' ' || c=='\t' || c=='\n' || c=='\r');
  };
  
  eqn.erase( std::remove_if(eqn.begin(), eqn.end(), is_whitespace_fcn), eqn.end());
  
  //Leave unicode alone, dont erase non-print chars
  //eqn.erase( std::remove_if(eqn.begin(), eqn.end(), [](char c){ !isprint(static_cast<unsigned char>(c)); } ), eqn.end());
  
  //grab all inner-most parenthesis
  std::map<std::string,int> argcounts;
  size_t open_pos = eqn.find('(');
  while( open_pos != string::npos )
  {
    const size_t close_pos = eqn.find( ')', open_pos + 1 );
    if( close_pos == string::npos )
      break;
    
    bool found_open = false;
    for( size_t i = open_pos+1; !found_open && i < close_pos; ++i )
      found_open = (eqn[i] == '(');
    
    if( !found_open )
    {
      const string substr = eqn.substr(open_pos+1, close_pos-open_pos-1);
      auto iter = argcounts.find(substr);
      if( iter == end(argcounts) )
        argcounts[substr] = 1;
      else
        iter->second += 1;
    }//
    
    open_pos = eqn.find( '(', open_pos + 1 );
  }//while( open_pos != string::npos )
  
  int max_occurances = 0;
  for( const auto x : argcounts )
    max_occurances = std::max( max_occurances, x.second );
  
  for( const auto x : argcounts )
    if( x.second == max_occurances )
      return x.first;
  
  //Default
  return "x";
}//void find_variable_name( std::string eqn )
  
  
float FormulaWrapper::efficiency( const float x )
{
  try
  {
    std::lock_guard<std::mutex> lock( m_mutex );
    *m_value = x;
    return static_cast<float>( m_parser->Eval().GetFloat() );
  }catch( mup::ParserError &e )
  {
    throw std::runtime_error( "Error evaluating expression \"" + e.GetExpr() + "\": " + e.GetMsg() );
  }//try / catch
}//float efficiency( const float x )
  
double FormulaWrapper::operator()( const float x )
{
  return efficiency(x);
}


DetectorPeakResponse::DetectorPeakResponse()
  : m_name( "DetectorPeakResponse" ),
    m_description( "" ),
    m_detectorDiameter( 0.0f ),
    m_efficiencyEnergyUnits( static_cast<float>(PhysicalUnits::keV) ),
    m_resolutionForm( DetectorPeakResponse::kNumResolutionFnctForm ),
    m_efficiencySource( DetectorPeakResponse::kUnknownEfficiencySource ),
    m_efficiencyForm( DetectorPeakResponse::kNumEfficiencyFnctForms ),
    m_hash( 0 ),
    m_parentHash( 0 )
{
} //DetectorPeakResponse()

DetectorPeakResponse::DetectorPeakResponse( const std::string &name,
                                            const std::string &descrip )
  : m_user( -1 ),
    m_name( name ),
    m_description( descrip ),
    m_detectorDiameter( 0.0 ),
    m_efficiencyEnergyUnits( static_cast<float>(PhysicalUnits::keV) ),
    m_resolutionForm( DetectorPeakResponse::kNumResolutionFnctForm ),
    m_efficiencySource( DetectorPeakResponse::kUnknownEfficiencySource ),
    m_efficiencyForm( DetectorPeakResponse::kNumEfficiencyFnctForms ),
    m_hash( 0 ),
    m_parentHash( 0 )
{
} //DetectorPeakResponse( const std::string &name, const std::string &descrip )


DetectorPeakResponse::~DetectorPeakResponse()
{
} //~DetectorPeakResponse()


uint64_t DetectorPeakResponse::hashValue() const
{
  return m_hash;
}

uint64_t DetectorPeakResponse::parentHashValue() const
{
  return m_parentHash;
}


void DetectorPeakResponse::setParentHashValue( const uint64_t val )
{
  m_parentHash = val;
}


void DetectorPeakResponse::computeHash()
{
  std::size_t seed = 0;
  
  boost::hash_combine( seed, m_name );
  boost::hash_combine( seed, m_description );
  
  boost::hash_combine( seed, m_detectorDiameter );
  boost::hash_combine( seed, m_efficiencyEnergyUnits );
  boost::hash_combine( seed, m_resolutionForm );
  for( const float val : m_resolutionCoeffs )
    boost::hash_combine( seed, val );
  
  boost::hash_combine( seed, m_efficiencySource );
  boost::hash_combine( seed, m_efficiencyForm );
  
  for( const EnergyEfficiencyPair &a : m_energyEfficiencies )
  {
    boost::hash_combine( seed, a.energy );
    boost::hash_combine( seed, a.efficiency );
  }
  
  boost::hash_combine( seed, m_efficiencyFormula );

  for( const float val : m_expOfLogPowerSeriesCoeffs )
    boost::hash_combine( seed, val );
  
  m_hash = seed;
}//void computeHash()


//isValid(): Tells wether or not this DetectorPeakResponse has been properly
//  defined and can be used
bool DetectorPeakResponse::isValid() const
{
  switch( m_efficiencyForm )
  {
    case kEnergyEfficiencyPairs:
      return !m_energyEfficiencies.empty();
    case kFunctialEfficienyForm:
      return !m_efficiencyFormula.empty();
    case kExpOfLogPowerSeries:
      return !m_expOfLogPowerSeriesCoeffs.empty();
    case kNumEfficiencyFnctForms:
      break;
  }//switch( m_efficiencyForm )
  
  return false;
}//bool isValid() const


//hasResolutionInfo(): Tells whether or not this DetectorPeakResponse has
//  (valid) detector resolution information
bool DetectorPeakResponse::hasResolutionInfo() const
{
  switch( m_resolutionForm )
  {
    case kGadrasResolutionFcn:
      return (m_resolutionCoeffs.size() == 3);
    case kSqrtPolynomial:
      return (m_resolutionCoeffs.size() > 1);
    case kNumResolutionFnctForm:
      return false;
  }//switch( m_resolutionForm )
  
  return false;
}//bool hasResolutionInfo() const


//reset(): clears all data, and causes isValid() to return false
void DetectorPeakResponse::reset()
{
  m_name.clear();
  m_description.clear();
  m_detectorDiameter = 0.0;
  m_efficiencyEnergyUnits = static_cast<float>( PhysicalUnits::keV );
  m_resolutionForm = kNumResolutionFnctForm;
  m_efficiencyForm = kNumEfficiencyFnctForms;
  m_efficiencySource = kUnknownEfficiencySource;
  m_resolutionCoeffs.clear();
  m_energyEfficiencies.clear();
  m_hash = m_parentHash = 0;
  
/*
  m_name = "Flat";
  m_description = "Flat response function; no resolution information";
  m_detectorDiameter = 1.0*PhysicalUnits::cm;
  EnergyEfficiencyPair p;
  p.energy = 0.0;
  p.efficiency = 1.0;
  m_energyEfficiencies.push_back( p );
  p.energy = 10000.0 * PhysicalUnits::keV;
  p.efficiency = 1.0;
  m_energyEfficiencies.push_back( p );
*/
}//void reset()


bool DetectorPeakResponse::operator==( const DetectorPeakResponse &rhs ) const
{
  if( isValid() != rhs.isValid() )
    return false;
  return (m_name==rhs.m_name
          && m_description==rhs.m_description
          && fabs(m_detectorDiameter-rhs.m_detectorDiameter)<0.000001
          && fabs(m_efficiencyEnergyUnits-rhs.m_efficiencyEnergyUnits)<0.000001f
          && m_resolutionForm==rhs.m_resolutionForm
          && m_efficiencyForm==rhs.m_efficiencyForm
          && m_resolutionCoeffs==rhs.m_resolutionCoeffs
          && m_energyEfficiencies==rhs.m_energyEfficiencies
          && m_hash==rhs.m_hash
          && m_parentHash==rhs.m_parentHash );
}//operator==


DetectorPeakResponse::EfficiencyDefinitionSource
                                DetectorPeakResponse::efficiencySource() const
{
  return m_efficiencySource;
}//EfficiencyDefinitionSource efficiencySource() const

float DetectorPeakResponse::efficiencyEnergyUnits() const
{
  return m_efficiencyEnergyUnits;
}

void DetectorPeakResponse::setEfficiencySource(
                    const DetectorPeakResponse::EfficiencyDefinitionSource src )
{
  m_efficiencySource = src;
  computeHash();
}//void setEfficiencySource( const EfficiencyDefinitionSource src )


void DetectorPeakResponse::fromEnergyEfficiencyCsv( std::istream &input,
                                    const float detectorDiameter,
                                    const float energyUnits )
{
  m_energyEfficiencies.clear();
  m_detectorDiameter = detectorDiameter;
  m_efficiencySource = kUserUploadedEfficiencyCsv;
  m_efficiencyEnergyUnits = energyUnits;

  if( m_detectorDiameter <= 0.0 || IsInf(m_detectorDiameter) || IsNan(m_detectorDiameter) )
    throw runtime_error( "Detector diameter must be greater than 0.0" );
  
  if( m_efficiencyEnergyUnits <= 0.0 || IsInf(m_efficiencyEnergyUnits) || IsNan(m_efficiencyEnergyUnits) )
    throw runtime_error( "Energy units must be greater than 0.0" );

  
  bool gotline = false;
  int linen = 0;
  string line;
  while( UtilityFunctions::safe_get_line( input, line ) )
  {
    ++linen;
    UtilityFunctions::trim( line );

    if( line.empty() || line[0] == '#' )
      continue;
    
    if( !isdigit(line[0]) )
    {
      if( !gotline )
        continue;
      throw runtime_error( "Invalid efficiency file.  After the CSV file header"
                           ", all lines must begin with a number (the energy"
                           " and then the % efficiency)" );
    }//if( !isdigit(line[0]) )

    vector<float> fields;
    const bool ok = UtilityFunctions::split_to_floats( line.c_str(),
                                                       line.size(), fields );
    
    if( (fields.size() < 2) || !ok )
    {
      if( !gotline )
        continue;
      
      throw runtime_error( "Invalid efficiency file.  After the CSV file header"
                           ", all lines must have at least two numbers: "
                           "energy and % efficiency. Further comma seperated "
                           "numbers are alllowed, but no other text" );
    }//if( fields.size() < 2 )
    
    EnergyEfficiencyPair point;
    point.energy = fields[0] * SandiaDecay::keV;
    point.efficiency = fields[1] / 100.0; //files have it listed as a percentage
    
    gotline = true;
    m_energyEfficiencies.push_back( point );
    
    if( linen > 5000 || m_energyEfficiencies.size() > 3000 )
      throw runtime_error( "Efficiency CSV files can have a max of 5000 total lines or 3000 energy efficiency pairs" );
    
  }//while( UtilityFunctions::safe_get_line(input) )
  
  
  
  if( m_energyEfficiencies.size() < 2 )
    throw runtime_error( "File was not a valid efficiency file, need at least"
                        " two efficiency points." );
  
  //Lets check if the numbers were strickly increasing or strickly decreasing
  const bool increasing = (m_energyEfficiencies[1].energy > m_energyEfficiencies[0].energy);
  for( size_t i = 2; i < m_energyEfficiencies.size(); ++i )
  {
    //We'll let two be energied be the exact same in a row
    const float laste = m_energyEfficiencies[i-1].energy;
    const float thise = m_energyEfficiencies[i].energy;
    
    if( IsNan(laste) || IsNan(thise)
        || IsInf(laste) || IsInf(thise) )
      throw runtime_error( "NaN or Inf energy detected from CSV file" );
    
    if( laste<0.0f || thise<0.0f )
      throw runtime_error( "Energy can not be less than zero in CSV efficiecny file" );
    
    if( m_energyEfficiencies[i-1].efficiency < 0.0f
       || m_energyEfficiencies[i].efficiency < 0.0f
       || m_energyEfficiencies[i-1].efficiency > 1.0f
       || m_energyEfficiencies[i].efficiency > 1.0f )
      throw runtime_error( "Intrinsic efficiency can not be less than zero or"
                           " greater than 100 in CSV efficiency file" );
    
    if( fabs(laste - thise) < 1.0E-5*std::max(fabs(thise),fabs(laste)) )
      continue;
    
    const bool bigger = (thise > laste);
    if( increasing != bigger )
      throw runtime_error( "Energies in efficiency CSV must be strickly "
                           "increasing or strickly decreasing" );
  }//for( size_t i = 2; i < m_energyEfficiencies.size(); ++i )
  

  if( !increasing )
    std::reverse( m_energyEfficiencies.begin(), m_energyEfficiencies.end() );
  
  m_efficiencyForm = kEnergyEfficiencyPairs;
  
  computeHash();
}//void fromEnergyEfficiencyCsv(...)



void DetectorPeakResponse::setIntrinsicEfficiencyFormula( const string &fcnstr,
                                                          const float diameter,
                                                          const float energyUnits )
{
  std::shared_ptr<FormulaWrapper > expression
                                = std::make_shared<FormulaWrapper>( fcnstr );
  
  m_efficiencyForm = kFunctialEfficienyForm;
  m_efficiencyFormula = fcnstr;
  m_efficiencySource = kUserEfficiencyEquationSpecified;
  m_efficiencyFcn = boost::bind( &FormulaWrapper::efficiency, expression, _1  );
  m_detectorDiameter = diameter;
  m_efficiencyEnergyUnits = energyUnits;
  
  computeHash();
}//void setIntrinsicEfficiencyFormula( const std::string &fcn )



void DetectorPeakResponse::fromGadrasDefinition( std::istream &csvFile,
                                                 std::istream &detDatFile )
{
  float detWidth = 0.0, heightToWidth = 0.0;

  m_efficiencyEnergyUnits = static_cast<float>( PhysicalUnits::keV );
  m_resolutionCoeffs.clear();
  m_resolutionForm = kGadrasResolutionFcn;
  m_resolutionCoeffs.resize( 3, 0.0 );

  string line;
  while( UtilityFunctions::safe_get_line( detDatFile, line ) )
  {
    vector<string> parts;
    UtilityFunctions::trim( line );
    if( line.empty() || !isdigit(line[0]) )
      continue;

    UtilityFunctions::split( parts, line, " \t" );
    
    try
    {
      const int parnum = std::stoi( parts.at(0) );
      const float value = static_cast<float>( std::stod( parts.at(1) ) );
//      const int value2 = std::stoi( parts.at(2) );
//      const string descrip = parts.at(3);
      switch( parnum )
      {
        case 6: m_resolutionCoeffs[0] = value; break;
        case 7: m_resolutionCoeffs[1] = value; break;
        case 8: m_resolutionCoeffs[2] = value; break;
        case 11: detWidth = value;      break;
        case 12: heightToWidth = value; break;
      }//switch( parnum )
    }catch(...)
    {
      cerr << "\nError reading line \"" << line << "\"" << endl;
      continue;
    }
  }//while( UtilityFunctions::safe_get_line( detDatFile, line ) )
  
  if( detWidth<=0.0 || heightToWidth<=0.0 )
    throw runtime_error( "Couldnt find detector dimensions in the Detector.dat file" );
  
  const float surfaceArea = detWidth * detWidth * heightToWidth;
  //const float diam = (4.0f/3.14159265359f) * sqrt(surfaceArea) * static_cast<float>(PhysicalUnits::cm);
  const float diam = 2.0f*sqrt(surfaceArea/3.14159265359f) * static_cast<float>(PhysicalUnits::cm);
  
  fromEnergyEfficiencyCsv( csvFile, diam, static_cast<float>(PhysicalUnits::keV) );
  
  m_efficiencySource = kGadrasEfficiencyDefintion;
  
  computeHash();
}//void fromGadrasDefinition(...)


void DetectorPeakResponse::fromGadrasDirectory( const std::string &dir )
{
  string eff_file, dat_file;
  
  //We will go through a little bit of trouble to not be case-sensitive
  //  For efficiency since some of the Detector.dat files that come with
  //  GADRAS are all uppercase.
  const vector<string> files = UtilityFunctions::ls_files_in_directory( dir );
  for( size_t i = 0; (eff_file.empty() || dat_file.empty()) && i < files.size(); ++i )
  {
    const auto &file = files[i];
    auto filename = UtilityFunctions::filename(file);
    if( UtilityFunctions::iequals( filename, "Efficiency.csv") )
      eff_file = filename;
    else if( UtilityFunctions::iequals( filename, "Detector.dat") )
      dat_file = filename;
  }
  
  if( dat_file.empty() )
    throw runtime_error( "Directory '" + dir + "' did not contain a Detector.dat file." );
  
  if( eff_file.empty() )
    throw runtime_error( "Directory '" + dir + "' did not contain a Efficiency.csv file." );
  
  ifstream csv( eff_file.c_str(), ios_base::binary|ios_base::in );
  ifstream datFile( dat_file.c_str(), ios_base::binary|ios_base::in );
  
  if( !csv.is_open() )
    throw runtime_error( "Could not open file '" + eff_file + "'." );
  
  if( !datFile.is_open() )
    throw runtime_error( "Could not open file '" + dat_file + "'." );
  
  fromGadrasDefinition( csv, datFile );
}//void fromGadrasDirectory( const std::string &dir )


void DetectorPeakResponse::fromExpOfLogPowerSeriesAbsEff(
                                   const std::vector<float> &coefs,
                                   const float charactDist,
                                   const float det_diam,
                                   const float energyUnits )
{
  if( coefs.empty() )
    throw runtime_error( "fromExpOfLogPowerSeriesAbsEff(...): invlaid input" );
  
  m_energyEfficiencies.clear();
  m_expOfLogPowerSeriesCoeffs.clear();
  m_efficiencyFcn = std::function<float(float)>();
  m_efficiencyFormula.clear();
  
  m_efficiencyForm = kExpOfLogPowerSeries;
  m_efficiencyEnergyUnits = energyUnits;
  m_detectorDiameter = det_diam;
  m_efficiencySource = kRelativeEfficiencyDefintion;
  m_expOfLogPowerSeriesCoeffs = coefs;
  
  //now we need to account for characterizationDist
  if( charactDist > 0.0f )
  {
    //x^n * x^m = x^(n+m)
    const float gfactor = fractionalSolidAngle( det_diam, charactDist );
    m_expOfLogPowerSeriesCoeffs[0] += log( 1.0f/gfactor );
  }//if( characterizationDist > 0.0f )
  
  computeHash();
}//void fromExpOfLogPowerSeriesAbsEff


void DetectorPeakResponse::toXml( ::rapidxml::xml_node<char> *parent, 
                                  ::rapidxml::xml_document<char> *doc ) const
{
  using namespace rapidxml;
  
  char buffer[128];
  
  xml_node<char> *base_node = doc->allocate_node( node_element, "DetectorPeakResponse" );
  parent->append_node( base_node );
  
  snprintf( buffer, sizeof(buffer), "%i", sm_xmlSerializationVersion );
  const char *value = doc->allocate_string( buffer );
  xml_attribute<char> *attr = doc->allocate_attribute( "version", value );
  base_node->append_attribute( attr );
  
  const char *val = doc->allocate_string( m_name.c_str() );
  xml_node<char> *node = doc->allocate_node( node_element, "Name", val );
  base_node->append_node( node );
  
  val = doc->allocate_string( m_description.c_str() );
  node = doc->allocate_node( node_element, "Description", val );
  base_node->append_node( node );
  
  snprintf( buffer, sizeof(buffer), "%1.8E", m_detectorDiameter );
  val = doc->allocate_string( buffer );
  node = doc->allocate_node( node_element, "DetectorDiameter", val );
  base_node->append_node( node );
  
  switch( m_efficiencySource )
  {
    case kGadrasEfficiencyDefintion:       val = "GadrasEfficiencyDefintion";       break;
    case kRelativeEfficiencyDefintion:     val = "RelativeEfficiencyDefintion";     break;
    case kUserUploadedEfficiencyCsv:       val = "UserUploadedEfficiencyCsv";       break;
    case kUserEfficiencyEquationSpecified: val = "UserEfficiencyEquationSpecified"; break;
    case kUnknownEfficiencySource:         val = "UnknownEfficiencySource";         break;
  }//switch( m_efficiencySource )
  
  node = doc->allocate_node( node_element, "EfficiencySource", val );
  base_node->append_node( node );
  
  snprintf( buffer, sizeof(buffer), "%g", m_efficiencyEnergyUnits );
  val = doc->allocate_string( buffer );
  node = doc->allocate_node( node_element, "EfficiencyEnergyUnits", val );
  base_node->append_node( node );
  
  switch( m_resolutionForm )
  {
    case kGadrasResolutionFcn:   val = "GadrasResolutionFcn"; break;
    case kSqrtPolynomial:        val = "SqrtPolynomial";      break;
    case kNumResolutionFnctForm: val = "Undefined";           break;
  }//switch( m_resolutionForm )

  node = doc->allocate_node( node_element, "ResolutionForm", val );
  base_node->append_node( node );

  
  switch( m_efficiencyForm )
  {
    case kEnergyEfficiencyPairs:  val = "EnergyEfficiencyPairs"; break;
    case kFunctialEfficienyForm:  val = "FunctialEfficienyForm"; break;
    case kExpOfLogPowerSeries:    val = "ExpOfLogPowerSeries";   break;
    case kNumEfficiencyFnctForms: val = "Undefined";             break;
  }//switch( m_efficiencyForm )
  
  node = doc->allocate_node( node_element, "EfficiencyForm", val );
  base_node->append_node( node );
  
  
  if( m_resolutionCoeffs.size() )
  {
    stringstream valstrm;
    for( size_t i = 0; i < m_resolutionCoeffs.size(); ++i )
      valstrm << (i?" ":"") << m_resolutionCoeffs[i];
    val = doc->allocate_string( valstrm.str().c_str() );
    node = doc->allocate_node( node_element, "ResolutionCoefficients", val );
    base_node->append_node( node );
  }//if( m_resolutionCoeffs.size() )
  
  if( m_energyEfficiencies.size() )
  {
    stringstream valstrm;
    for( size_t i = 0; i < m_energyEfficiencies.size(); ++i )
      valstrm << (i?" ":"") << m_energyEfficiencies[i].energy 
              << " " << m_energyEfficiencies[i].efficiency;
    val = doc->allocate_string( valstrm.str().c_str() );
    node = doc->allocate_node( node_element, "EnergyEfficiencies", val );
    base_node->append_node( node );    
  }//if( m_energyEfficiencies.size() )
  
  if( m_efficiencyFormula.size() )
  {
    val = doc->allocate_string( m_efficiencyFormula.c_str() );
    node = doc->allocate_node( node_element, "EfficiencyFormula", val );
    base_node->append_node( node );    
  }//if( m_efficiencyFormula.size() )
  
  if( m_expOfLogPowerSeriesCoeffs.size() )
  {
    stringstream valstrm;
    for( size_t i = 0; i < m_expOfLogPowerSeriesCoeffs.size(); ++i )
      valstrm << (i?" ":"") << m_expOfLogPowerSeriesCoeffs[i];
    val = doc->allocate_string( valstrm.str().c_str() );
    node = doc->allocate_node( node_element, "ExpOfLogPowerSeriesCoeffs", val );
    base_node->append_node( node );
  }//if( m_expOfLogPowerSeriesCoeffs.size() )
  
  stringstream hashstrm, parenthashstrm;
  hashstrm << m_hash;
  parenthashstrm << m_parentHash;
  
  val = doc->allocate_string( hashstrm.str().c_str() );
  node = doc->allocate_node( node_element, "Hash", val );
  base_node->append_node( node );
  
  val = doc->allocate_string( parenthashstrm.str().c_str() );
  node = doc->allocate_node( node_element, "ParentHash", val );
  base_node->append_node( node );
}//toXml(...)


void DetectorPeakResponse::fromXml( const ::rapidxml::xml_node<char> *parent )
{
  using namespace rapidxml;
  using ::rapidxml::internal::compare;
  
  if( !parent )
    throw runtime_error( "DetectorPeakResponse::fromXml(...): invalid input" );
  
  if( !compare( parent->name(), parent->name_size(), "DetectorPeakResponse", 20, false ) )
    throw std::logic_error( "DetectorPeakResponse::fromXml(...): invalid input node name" );
  
  const xml_attribute<char> *att = parent->first_attribute( "version", 7 );
  
  int version;
  if( !att || !att->value() || (sscanf(att->value(), "%i", &version)!=1) )
    throw runtime_error( "DetectorPeakResponse invalid version" );
  
  if( version != sm_xmlSerializationVersion )
    throw runtime_error( "Invalid DetectorPeakResponse version" );
  

  const xml_node<char> *node = parent->first_node( "Name", 4 );
  if( !node || !node->value() )
    throw runtime_error( "DetectorPeakResponse missing Name node" );
  m_name = node->value();
  
  
  node = parent->first_node( "Description", 11 );
  if( !node || !node->value() )
    throw runtime_error( "DetectorPeakResponse missing Description node" );
  m_description = node->value();
  
  node = parent->first_node( "DetectorDiameter", 16 );
  if( !node || !node->value() )
    throw runtime_error( "DetectorPeakResponse missing DetectorDiameter node" );
  m_detectorDiameter = atof( node->value() );
  
  node = parent->first_node( "EfficiencySource", 16 );
  if( !node || !node->value() )
    throw runtime_error( "DetectorPeakResponse missing EfficiencySource node" );

  if( compare(node->value(),node->value_size(),"GadrasEfficiencyDefintion",25,false) )
    m_efficiencySource = kGadrasEfficiencyDefintion;
  else if( compare(node->value(),node->value_size(),   "SimpleMassEfficiencyDefintion",29,false)
           || compare(node->value(),node->value_size(),"RelativeEfficiencyDefintion",27,false) )
    m_efficiencySource = kRelativeEfficiencyDefintion;
  else if( compare(node->value(),node->value_size(),"UserUploadedEfficiencyCsv",25,false) )
    m_efficiencySource = kUserUploadedEfficiencyCsv;
  else if( compare(node->value(),node->value_size(),"UserEfficiencyEquationSpecified",31,false) )
    m_efficiencySource = kUserEfficiencyEquationSpecified;
  else if( compare(node->value(),node->value_size(),"UnknownEfficiencySource",23,false) )
    m_efficiencySource = kUnknownEfficiencySource;
  else 
    throw runtime_error( "DetectorPeakResponse: invalid EfficiencySource value" );

  
  node = parent->first_node( "EfficiencyEnergyUnits", 21 );
  if( !node || !node->value() )
    throw runtime_error( "DetectorPeakResponse missing EfficiencyEnergyUnits node" );
  m_efficiencyEnergyUnits = atof( node->value() );
  
  
  node = parent->first_node( "ResolutionForm", 14 );
  if( !node || !node->value() )
    throw runtime_error( "DetectorPeakResponse missing ResolutionForm node" );
  
  if( compare(node->value(),node->value_size(),"GadrasResolutionFcn",19,false) )
    m_resolutionForm = kGadrasResolutionFcn;
  else if( compare(node->value(),node->value_size(),"SqrtPolynomial",14,false) )
    m_resolutionForm = kSqrtPolynomial;
  else if( compare(node->value(),node->value_size(),"Undefined",9,false) )
    m_resolutionForm = kNumResolutionFnctForm;
  else
    throw runtime_error( "DetectorPeakResponse: invalid ResolutionForm value" );
  
  
  node = parent->first_node( "EfficiencyForm", 14 );
  if( !node || !node->value() )
    throw runtime_error( "DetectorPeakResponse missing EfficiencyForm node" );
  
  if( compare(node->value(),node->value_size(),"EnergyEfficiencyPairs",21,false) )
    m_efficiencyForm = kEnergyEfficiencyPairs;
  else if( compare(node->value(),node->value_size(),"FunctialEfficienyForm",21,false) )
    m_efficiencyForm = kFunctialEfficienyForm;
  else if( compare(node->value(),node->value_size(),"ExpOfLogPowerSeries",19,false) )
    m_efficiencyForm = kExpOfLogPowerSeries;
  else if( compare(node->value(),node->value_size(),"Undefined",9,false) )
    m_efficiencyForm = kNumEfficiencyFnctForms;
  else
    throw runtime_error( "DetectorPeakResponse: invalid EfficiencyForm value" );
  
  m_resolutionCoeffs.clear();
  node = parent->first_node( "ResolutionCoefficients", 22 );
  if( node && node->value() )
    UtilityFunctions::split_to_floats( node->value(), node->value_size(), m_resolutionCoeffs );
  
  m_energyEfficiencies.clear();
  node = parent->first_node( "EnergyEfficiencies", 18 );
  if( node && node->value() )
  {
    vector<float> values;
    UtilityFunctions::split_to_floats( node->value(), node->value_size(), values );
    if( (values.size()%2) != 0 )
      throw runtime_error( "DetectorPeakResponse: invalid number of energy efficiency pairs" );
    
    for( size_t i = 0; i < values.size(); i += 2 )
    {
      EnergyEfficiencyPair p;
      p.energy = values[i];
      p.efficiency = values[i+1];
      m_energyEfficiencies.push_back( p );
    }
  }//if( node && node->value() )
  
  m_efficiencyFcn = std::function<float(float)>();
  m_efficiencyFormula.clear();
  node = parent->first_node( "EfficiencyFormula", 17 );
  if( node && node->value() )
  {
    m_efficiencyFormula = node->value();
    UtilityFunctions::trim( m_efficiencyFormula );
    
    if( m_efficiencyFormula.size() )
    {
      //if( m_efficiencySource != kUserEfficiencyEquationSpecified )
        //throw runtime_error( "An detector efficiency formula was specified but the the EfficiencySource had a different value" );
      try
      {
        auto expression = std::make_shared<FormulaWrapper>( m_efficiencyFormula );
        m_efficiencyFcn = boost::bind( &FormulaWrapper::efficiency, expression, _1  );
      }catch( std::exception &e )
      {
        throw runtime_error( "Invalid detector efficiency formula in XML: " + string(e.what()) );
      }
    }//if( m_efficiencyFormula.size() )
  }//if( node && node->value() )
  
  m_expOfLogPowerSeriesCoeffs.clear();
  node = parent->first_node( "ExpOfLogPowerSeriesCoeffs", 25 );
  if( node && node->value() )
    UtilityFunctions::split_to_floats( node->value(), node->value_size(), m_expOfLogPowerSeriesCoeffs );
  
  
  node = parent->first_node( "Hash", 4 );
  if( !node || !node->value() )
    throw runtime_error( "DetectorPeakResponse missing Hash node" );
  if( !(stringstream(node->value()) >> m_hash) )
    throw runtime_error( "DetectorPeakResponse invalid Hash" );
  
  node = parent->first_node( "ParentHash", 10 );
  if( !node || !node->value() )
    throw runtime_error( "DetectorPeakResponse missing ParentHash node" );
  if( !(stringstream(node->value()) >> m_parentHash) )
    throw runtime_error( "DetectorPeakResponse invalid ParentHash" );
}//void fromXml( ::rapidxml::xml_node<char> *parent )



#if( PERFORM_DEVELOPER_CHECKS )
void DetectorPeakResponse::equalEnough( const DetectorPeakResponse &lhs,
                                        const DetectorPeakResponse &rhs )
{
  char buffer[512];
  
  if( lhs.m_name != rhs.m_name )
  {
    snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: name of LHS ('%s')"
             " doesnt match RHS ('%s')",
             lhs.m_name.c_str(), rhs.m_name.c_str() );
    throw runtime_error(buffer);
  }

  if( lhs.m_description != rhs.m_description )
  {
    snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: name of LHS ('%s')"
             " doesnt match RHS ('%s')",
             lhs.m_description.c_str(), rhs.m_description.c_str() );
    throw runtime_error(buffer);
  }
  
  if( fabs(lhs.m_detectorDiameter - rhs.m_detectorDiameter) > 0.001 )
  {
    snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: detector diameter"
              " of LHS (%1.8E) doesnt match RHS (%1.8E)",
             lhs.m_detectorDiameter, rhs.m_detectorDiameter );
    throw runtime_error(buffer);
  }

  if( fabs(lhs.m_efficiencyEnergyUnits - rhs.m_efficiencyEnergyUnits) > 0.001 )
  {
    snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: efficiency units"
             " of LHS (%1.8E) doesnt match RHS (%1.8E)",
             lhs.m_detectorDiameter, rhs.m_detectorDiameter );
    throw runtime_error(buffer);
  }

  if( lhs.m_resolutionForm != rhs.m_resolutionForm )
  {
    snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: resolution"
              " functional form of LHS (%i) doesnt match RHS (%i)",
             int(lhs.m_resolutionForm), int(rhs.m_resolutionForm) );
    throw runtime_error(buffer);
  }
  
  
  if( lhs.m_resolutionCoeffs.size() != rhs.m_resolutionCoeffs.size() )
  {
    snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: size of"
             " resolution coefficients of LHS (%i) doesnt match RHS (%i)",
             int(lhs.m_resolutionCoeffs.size()),
             int(rhs.m_resolutionCoeffs.size()) );
    throw runtime_error(buffer);
  }
  
  for( size_t i = 0; i < lhs.m_resolutionCoeffs.size(); ++i )
  {
    const float a = lhs.m_resolutionCoeffs[i];
    const float b = rhs.m_resolutionCoeffs[i];
    if( fabs(a-b) > (1.0E-5 * std::max(fabs(a),fabs(b))) )
    {
      snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: resolution"
      " coefficient %i of LHS (%1.8E) doesnt match RHS (%1.8E)", int(i), a, b );
      throw runtime_error( buffer );
    }
  }
  
  if( lhs.m_efficiencySource != rhs.m_efficiencySource )
  {
    snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: efficiency source"
             " of LHS (%i) doesnt match RHS (%i)",
             int(lhs.m_efficiencySource),
             int(rhs.m_efficiencySource) );
    throw runtime_error(buffer);
  }
  
  if( lhs.m_efficiencyForm != rhs.m_efficiencyForm )
  {
    snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: efficiency"
              " functional form of LHS (%i) doesnt match RHS (%i)",
             int(lhs.m_efficiencyForm),
             int(rhs.m_efficiencyForm) );
    throw runtime_error(buffer);
  }
  
  
  if( lhs.m_energyEfficiencies.size() != rhs.m_energyEfficiencies.size() )
  {
    snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: size of"
             " efficiency coefficients of LHS (%i) doesnt match RHS (%i)",
             int(lhs.m_energyEfficiencies.size()),
             int(rhs.m_energyEfficiencies.size()) );
    throw runtime_error(buffer);
  }
  
  for( size_t i = 0; i < lhs.m_energyEfficiencies.size(); ++i )
  {
    float a = lhs.m_energyEfficiencies[i].energy;
    float b = rhs.m_energyEfficiencies[i].energy;
    if( fabs(a-b) > (1.0E-5 * std::max(fabs(a),fabs(b))) )
    {
      snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: efficiency"
               " energy %i of LHS (%1.8E) doesnt match RHS (%1.8E)", int(i), a, b );
      throw runtime_error( buffer );
    }
    
    a = lhs.m_energyEfficiencies[i].efficiency;
    b = rhs.m_energyEfficiencies[i].efficiency;
    if( fabs(a-b) > 1.0E-5 )
    {
      snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: efficiency"
               " %i of LHS (%1.8E) doesnt match RHS (%1.8E)", int(i), a, b );
      throw runtime_error( buffer );
    }
  }
  
  
  if( lhs.m_efficiencyFormula != rhs.m_efficiencyFormula )
  {
    snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: efficiency formula"
              " of LHS ('%s') doesnt match RHS ('%s')",
             lhs.m_efficiencyFormula.c_str(), rhs.m_efficiencyFormula.c_str() );
    throw runtime_error(buffer);
  }
  
  if( (!lhs.m_efficiencyFcn) != (!rhs.m_efficiencyFcn) )
  {
    snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: availability of"
             " efficiency formula of LHS (%s) doesnt match RHS (%s)",
             (!lhs.m_efficiencyFcn?"missing":"available"),
             (!rhs.m_efficiencyFcn?"missing":"available") );
    throw runtime_error(buffer);
  }
  
  
  if( lhs.m_expOfLogPowerSeriesCoeffs.size() != rhs.m_expOfLogPowerSeriesCoeffs.size() )
  {
    snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: size of"
             " exponential of log power series coefficients of LHS (%i) doesnt"
             " match RHS (%i)",
             int(lhs.m_expOfLogPowerSeriesCoeffs.size()),
             int(rhs.m_expOfLogPowerSeriesCoeffs.size()) );
    throw runtime_error(buffer);
  }
  
  for( size_t i = 0; i < lhs.m_expOfLogPowerSeriesCoeffs.size(); ++i )
  {
    const float a = lhs.m_expOfLogPowerSeriesCoeffs[i];
    const float b = rhs.m_expOfLogPowerSeriesCoeffs[i];
    if( fabs(a-b) > (1.0E-5 * std::max(fabs(a),fabs(b))) )
    {
      snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: exponential of"
                " log power series coefficient %i of LHS (%1.8E)"
                " doesnt match RHS (%1.8E)", int(i), a, b );
      throw runtime_error( buffer );
    }
  }
  
  if( lhs.m_hash != rhs.m_hash )
  {
    //should use PRIu64 instead of casting, but whatever
    snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: hash of LHS (%llud)"
              " doesnt match RHS (%llud)",
              (long long unsigned int)lhs.m_hash,
              (long long unsigned int)rhs.m_hash );
    throw runtime_error(buffer);
  }

  if( lhs.m_parentHash != rhs.m_parentHash )
  {
    snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: hash of LHS (%llud)"
              " doesnt match RHS (%llud)",
              (long long unsigned int)lhs.m_parentHash,
              (long long unsigned int)rhs.m_parentHash );
    throw runtime_error(buffer);
  }
}//void equalEnough(...)
#endif //PERFORM_DEVELOPER_CHECKS



float DetectorPeakResponse::fractionalSolidAngle( const float detDiam,
                                                   const float D,
                                                   const float source_rad )
{
  /*
    See also page 119 in Knoll for approximations where source and detector
    diameters are somewhat comparible to distance of source from detector
    An unverified implementation of that is:
  */
  const float alpha = pow(source_rad/D,2.0f);
  const float beta = pow(0.5f*detDiam/D,2.0f);
  const float F1 = ((5.0f/16.0f)*beta/pow(1.0f+beta,3.5f)) - ((35.0f/64.0f)*beta*beta/pow(1.0f+beta,4.5f));
  const float F2 = ((35.0f/128.0f)*beta/pow(1.0f+beta,4.5f)) - ((315.0/256.0)*beta*beta/pow(1.0f+beta,11.0f/12.0f)) + ((1155.0f/1024.0f)*beta*beta*beta/pow(1.0f+beta,6.5f));

  const float omega = 0.5f * ( 1.0f - (1.0f/sqrt(1.0f+beta)) - ((3.0f/8.0f)*alpha*beta/pow(1.0f+beta,2.5f)) + alpha*alpha*F1 - alpha*alpha*alpha*F2 );

  return omega;
}//float DetectorPeakResponse::fractionalSolidAngle(...)


float DetectorPeakResponse::fractionalSolidAngle( const float detDiam,
                                                  const float D )
{
  const float r = 0.5f * detDiam;
  return 0.5f*(1.0f - (D/sqrt(D*D+r*r)));
}//float fractionalSolidAngle(...)


float DetectorPeakResponse::efficiency( const float energy,
                                        const float dist ) const
{
  const float fracSolidAngle = fractionalSolidAngle( m_detectorDiameter, dist );
  return fracSolidAngle * intrinsicEfficiency( energy );
}//float efficiency( const float energy ) const


const vector<DetectorPeakResponse::EnergyEfficiencyPair> &DetectorPeakResponse::getEnergyEfficiencyPair() const
{
  return m_energyEfficiencies;
}//std::vector<EnergyEfficiencyPair> DetectorPeakResponse::getEnergyEfficiencyPair()


float DetectorPeakResponse::intrinsicEfficiencyFromPairs( float energy ) const
{
  energy /= m_efficiencyEnergyUnits;
  
  if( m_energyEfficiencies.size() < 2 )
    throw runtime_error("DetectorPeakResponse objects must be initialized "
                        "before calling intrinsicEfficiencyFromPairs(...)");
  return akimaInterpolate( energy, m_energyEfficiencies );
}//double DetectorPeakResponse::absoluteEfficiencyFromPairs( const float energy ) const




float DetectorPeakResponse::akimaInterpolate( const float z,
            const std::vector<DetectorPeakResponse::EnergyEfficiencyPair> &xy )
{
  //adapted from http://users-phys.au.dk/fedorov/nucltheo/Numeric/now/interp.pdf
  vector<DetectorPeakResponse::EnergyEfficiencyPair>::const_iterator pos;
  pos = lower_bound( xy.begin(), xy.end(), z );
  const size_t n = xy.size();
  
  if( pos == xy.begin() )
    return xy.front().efficiency;
  if( pos == xy.end() )
    return xy.back().efficiency;
  
  const size_t i = (pos-1) - xy.begin();
  
  const float x_i = xy[i].energy;
  const float x_p1 = xy[i+1].energy;
  const float y_i = xy[i].efficiency;
  const float y_p1 = xy[i+1].efficiency;
  
  if( n < 6 || i < 2 || (n-i) < 3  ) //We'll just use linear interpolation here
  {
    const float d = (z - x_i) / (x_p1 - x_i);
    return y_i + d*(y_p1 - y_i);
  }//if( n < 6 )
  
  const float dx = z - x_i;
  const float h_i = (x_p1 - x_i);
  const float A_i = calcA( i, xy );
  const float A_p1 = calcA( i+1, xy );
  const float p_i = (y_p1 - y_i) / h_i;
  const float c_i = (3.0f*p_i - 2.0f*A_i - A_p1) / h_i;
  const float d_i = (A_p1 + A_i - 2.0f*p_i) / h_i / h_i;
  
  return y_i + dx*(A_i + dx*(c_i + dx*d_i));
}//akimaInterpolate(...)


float DetectorPeakResponse::intrinsicEfficiencyFromFcn( float energy ) const
{
  energy /= m_efficiencyEnergyUnits;
  
  if( !m_efficiencyFcn )
    throw runtime_error( "DetectorPeakResponse objects must be initialized "
                         "before calling intrinsicEfficiencyFromFcn(...)" );
  return m_efficiencyFcn( energy );
}//double intrinsicEfficiencyFromFcn( const float energy ) const


float DetectorPeakResponse::intrinsicEfficiencyFromExpLnEqn( float energy ) const
{
  if( m_expOfLogPowerSeriesCoeffs.empty() )
    throw runtime_error( "DetectorPeakResponse objects must be initialized "
                         "before calling intrinsicEfficiencyFromExpLnEqn(...)" );
 
  energy /= m_efficiencyEnergyUnits;
  return expOfLogPowerSeriesEfficiency( energy, m_expOfLogPowerSeriesCoeffs );
}//float intrinsicEfficiencyFromExpLnEqn( float energy ) const


float DetectorPeakResponse::expOfLogPowerSeriesEfficiency( const float energy,
                                           const std::vector<float> &coefs )
{
  float exparg = 0.0f;
  const float x = log( energy );
  for( size_t i = 0; i < coefs.size(); ++i )
    exparg += coefs[i] * pow(x,static_cast<float>(i));
  return exp( exparg );
}//float expOfLogPowerSeriesEfficiency(...)


float DetectorPeakResponse::intrinsicEfficiency( const float energy ) const
{
  
  switch( m_efficiencyForm )
  {
    case kEnergyEfficiencyPairs:
      return intrinsicEfficiencyFromPairs( energy );
      
    case kFunctialEfficienyForm:
      return intrinsicEfficiencyFromFcn( energy );
      
    case kExpOfLogPowerSeries:
      return intrinsicEfficiencyFromExpLnEqn( energy );
      
    case kNumEfficiencyFnctForms:
      break;
  }//switch( m_efficiencyForm )

  
  throw runtime_error( "DetectorPeakResponse::intrinsicEfficiency:"
                       " undefined efficiency" );
}//float intrinsicEfficiency( const float energy ) const;



float DetectorPeakResponse::peakResolutionFWHM( float energy,
                                                ResolutionFnctForm fcnFrm,
                                                const std::vector<float> &pars )
{
  /*
   //The below is to compare for the Anthony NM detector, what FWHM we should
   //  get, verses what we actually get reading the Detector.dat from GADRAS
   std::vector<float> coeffs;
   coeffs.push_back( 0.1401463  );
   coeffs.push_back( 0.187422  );
   coeffs.push_back( -2.07732e-09 );
   
   float channel = find_bin_from_polynomial( energy, coeffs, 16384 );
   
   double expected_fwhm_channels = 6.220468 + 0.0005996688*channel + channel*channel*-1.25688e-10;
   double fwhm_kev = (0.187422 + channel*-2.07732e-09)* expected_fwhm_channels;
   cerr << "FWHM at energy=" << energy << " to be Expected=" << fwhm_kev
   << ", Actual = " << 2.35482 * peakResolutionSigma( energy )
   << "  - channel=" << channel
   << endl;
   */
  /*
   //From Greg T.:
   IF (E.GE.661.OR.A(6).EQ.0) THEN
   GetFWHM=6.61*A(7)*(E/661.)**A(8)
   ELSE IF (A(6).LT.0) THEN
   P=A(8)**(1./LOG(1.-A(6)))
   GetFWHM=6.61*A(7)*(E/661.)**P
   ELSE
   IF (A(6).GT.6.61*A(7)) THEN
   GetFWHM=A(6)
   ELSE
   A7=SQRT((6.61*A(7))**2-A(6)**2)/6.61
   GetFWHM=SQRT(A(6)**2+(6.61*A7*(E/661.)**A(8))**2)
   END IF
   END IF
   */
  
  
  switch( fcnFrm )
  {
    case kGadrasResolutionFcn:
    {
      if( pars.size() != 3 )
        throw std::runtime_error( "DetectorPeakResponse::peakResolutionSigma():"
                                 " pars not defined" );
      
      energy /= PhysicalUnits::keV;
      const float &a = pars[0];
      const float &b = pars[1];
      const float &c = pars[2];
      
      if( energy >= 661.0f || fabs(a)<float(1.0E-6) )
        return 6.61f * b * pow(energy/661.0f, c);
      
      if( a < 0.0 )
      {
        const float p = pow( c, float(1.0f/log(1.0f-a)) );
        return 6.61f * b * pow(energy/661.0f, p);
      }//if( a < 0.0 )
      
      if( a > 6.61f*b )
        return a;
      
      const float A7 = sqrt( pow(float(6.61f*b), float(2.0f))-a*a )/6.61f;
      return sqrt(a*a + pow(float(6.61f * A7 * pow(energy/661.0f, c)), 2.0f));
    }//case kGadrasResolutionFcn:
      
    case kSqrtPolynomial:
    {
      //FWHM = sqrt(A1+A2*energy+A3*energy*energy + A4/energy)
      if( pars.size() < 2 )
        throw runtime_error( "DetectorPeakResponse::peakResolutionSigma():"
                             " pars not defined" );
      const float A1 = pars[0];
      const float A2 = pars[1];
//      const float A3 = pars.size()>2 ? pars[2] : 0.0f;
//      const float A4 = pars.size()>3 ? pars[3] : 0.0f;
      
      return A1 + A2 *std::sqrt(energy);
//      return std::sqrt( A1 + A2*energy + A3*energy*energy + A4/energy);
    }//case kSqrtPolynomial:
      
    case kNumResolutionFnctForm:
      throw std::runtime_error( "DetectorPeakResponse::peakResolutionSigma():"
                                " Resolution not defined" );
    break;
  }//switch( m_resolutionForm )
  
  //Lets keep MSVS happy
  assert(0);
  return 0.0f;
}//static float peakResolutionFwhmGadras(...)


float DetectorPeakResponse::peakResolutionSigma( const float energy,
                                                  ResolutionFnctForm fcnFrm,
                                              const std::vector<float> &pars )
{
  return peakResolutionFWHM( energy, fcnFrm, pars ) / 2.35482f;
}//static double peakResolutionSigma(...)


float DetectorPeakResponse::peakResolutionFWHM( const float energy ) const
{
  return peakResolutionFWHM( energy, m_resolutionForm, m_resolutionCoeffs );
}//double peakResolutionFWHM( float energy ) const


float DetectorPeakResponse::peakResolutionSigma( const float energy ) const
{
  return peakResolutionFWHM(energy) / 2.35482f;
}//double peakResolutionSigma( float energy ) const


float DetectorPeakResponse::detectorDiameter() const
{
  return m_detectorDiameter;
}

const std::string &DetectorPeakResponse::efficiencyFormula() const
{
  return m_efficiencyFormula;
}

const std::string &DetectorPeakResponse::name() const
{
  return m_name;
}

const std::string &DetectorPeakResponse::description() const
{
  return m_description;
}

void DetectorPeakResponse::setName( const std::string &name )
{
  m_name = name;
  computeHash();
}

void DetectorPeakResponse::setDescription( const std::string &descrip )
{
  m_description = descrip;
  computeHash();
}

namespace
{
  template<class T> struct index_compare_descend
  {
    index_compare_descend(const T arr) : arr(arr) {} //pass the actual values you want sorted into here
    bool operator()(const size_t a, const size_t b) const
    {
      return arr[a] > arr[b];
    }
    const T arr;
  };//struct index_compare_descend
}//namespace


DetectorPeakResponse::PeakInput_t DetectorPeakResponse::removeOutlyingWidthPeaks(
                           DetectorPeakResponse::PeakInput_t peaks,
                           DetectorPeakResponse::ResolutionFnctForm fnctnlForm,
                                        const vector<float> &coefficients )
{
  const double npeaks = static_cast<double>( peaks->size() );
  if( npeaks < 5 )
    return peaks;
  const size_t ndel_max = static_cast<size_t>( floor(0.2*npeaks) );
  
  vector<double> weights;
  double mean_weight = 0.0;
  for( const PeakModel::PeakShrdPtr peak : *peaks )
  {
    const double predicted_sigma = DetectorPeakResponse::peakResolutionSigma( peak->mean(), fnctnlForm, coefficients );
    const double chi2 = DetectorPeakResponse::peak_width_chi2( predicted_sigma, *peak );
    mean_weight += chi2/npeaks;
    weights.push_back( chi2 );
  }//for( const EnergySigma &es : m_energies_and_sigmas )
  
  vector<size_t> indices( npeaks );
  for( size_t i = 0; i < npeaks; ++i )
    indices[i] = i;
  std::sort( indices.begin(), indices.end(), index_compare_descend<vector<double>&>(weights) );
  
  size_t lastInd = 0;
  while( lastInd <= ndel_max && weights[indices[lastInd]] > 2.5*mean_weight )
    ++lastInd;
  indices.erase( indices.begin() + lastInd, indices.end() );
  
  std::shared_ptr< deque< std::shared_ptr<const PeakDef> > > reduced_peaks( new deque< PeakModel::PeakShrdPtr >() );
  for( size_t i = 0; i < npeaks; ++i )
    if( find( indices.begin(), indices.end(), i) == indices.end() )
      reduced_peaks->push_back( peaks->operator[](i) );
  
  return reduced_peaks;
}//removeOutlyingWidthPeaks(...)


void DetectorPeakResponse::fitResolution( DetectorPeakResponse::PeakInput_t peaks,
                    const std::shared_ptr<const Measurement> meas,
                    const DetectorPeakResponse::ResolutionFnctForm fnctnlForm )
{
  if( !peaks )
    throw runtime_error( "DetectorPeakResponse::fitResolution(...): invalid input" );
  
  vector<float> coefficients = m_resolutionCoeffs;
  performResolutionFit( peaks, meas, fnctnlForm, coefficients );
  peaks = removeOutlyingWidthPeaks( peaks, fnctnlForm, coefficients );
  performResolutionFit( peaks, meas, fnctnlForm, coefficients );
  
  m_resolutionCoeffs = coefficients;
  m_resolutionForm = fnctnlForm;
  
  computeHash();
}//void fitResolution(...)


double DetectorPeakResponse::peak_width_chi2( double predicted_sigma,
                                              const PeakDef &peak )
{
  double measured_sigma = peak.sigma();
  double measured_sigma_uncert = peak.sigmaUncert();
  if( measured_sigma_uncert <= 0.01 )
    measured_sigma_uncert = 1.0;
  
  const double chi = (measured_sigma - predicted_sigma)/measured_sigma_uncert;
  return chi*chi;
}//peak_width_chi2(...)



namespace
{
class DetectorResolutionFitness
: public ROOT::Minuit2::FCNBase
{
protected:
  std::deque<PeakModel::PeakShrdPtr> m_peaks;
  DetectorPeakResponse::ResolutionFnctForm m_form;
  
public:
  DetectorResolutionFitness( const std::deque<PeakModel::PeakShrdPtr> &peaks,
                             DetectorPeakResponse::ResolutionFnctForm form )
  : m_peaks( peaks ),
    m_form( form )
  {
  }
  
  virtual ~DetectorResolutionFitness(){}
  
  
  virtual double DoEval( const std::vector<double> &x ) const
  {
    for( size_t i = 0; i < x.size(); ++i )
      if( IsInf(x[i]) || IsNan(x[i]) )
        return 999999.9;
   
    vector<float> fx( x.size() );
    for( size_t i = 0; i < x.size(); ++i )
      fx[i] = static_cast<float>( x[i] );
    
    double chi2 = 0.0;
    for( const PeakModel::PeakShrdPtr &peak : m_peaks )
    {
      if( !peak )  //probably isnt needed
        continue;      
      const float predicted = DetectorPeakResponse::peakResolutionSigma(
                                                    peak->mean(), m_form, fx );
      if( predicted <= 0.0 || IsNan(predicted) || IsInf(predicted) )
        return 999999.9;
      
      chi2 += DetectorPeakResponse::peak_width_chi2( predicted, *peak );
    }//for( const EnergySigma &es : m_energies_and_sigmas )
    
    return chi2 / x.size();
  }//DoEval();
  
  virtual double operator()( const std::vector<double> &x ) const
  {
    return DoEval( x );
  }
  
  DetectorResolutionFitness&	operator=( const DetectorResolutionFitness &rhs )
  {
    if( &rhs != this )
    {
      m_peaks = rhs.m_peaks;
      m_form = rhs.m_form;
    }
    return *this;
  }
  
  virtual double Up() const{ return 1.0; }
};//class DetectorResolutionFitness
}//namespace

void DetectorPeakResponse::performResolutionFit(
                      DetectorPeakResponse::PeakInput_t peaks,
                      const std::shared_ptr<const Measurement> meas,
                      const DetectorPeakResponse::ResolutionFnctForm fnctnlForm,
                      std::vector<float> &answer )
{
  if( !peaks || peaks->empty() )
    throw runtime_error( "DetectorPeakResponse::performResolutionFit(...):"
                        " no input peaks" );
  
  answer.clear();

  
  bool highres;
  if( !!meas )
  {
    highres = (meas->num_gamma_channels() > 3000);
  }else
  {
    const size_t index = peaks->size() / 2;
    const double ratio = (*peaks)[index]->fwhm() / (*peaks)[index]->mean();
    highres = (ratio < 0.03); //whatever, this is JIC anyway
  }//if( !!meas ) / else
  
  /*
High res.:
   Detective-EX100:         {1.54000, 0.26412, 0.33000}
   Detective-EX:            {1.34000, 0.27759, 0.31000}
   Detctive:                {1.77168, 0.23563, 0.57223}
   Falcon 5000:             {1.00000, 0.20028, 0.56000}
   MicroDetective:          {1.30000, 0.23731, 0.54000}
   Ortec IDM Portal         {1.11000, 0.23598, 0.39000}
   Transpec                 {1.33000, 0.25862, 0.28000}
Low res.:
   FieldSec:                {-2.00000, 6.80000, 0.70000}
   Gr130                    {-5.00000, 7.00000, 0.60000}
   GR135:                   {-7.00000, 6.84000, 0.55000}
   GR135 Plus:              {-7.00000, 6.84000, 0.55000}
   Identifinder-LaBr3:      { 3.00000, 3.30000, 0.60000}
   Identifinder-N:          {-7.00000, 7.50000, 0.55000}
   Identifinder-NGH         {-6.00000, 7.80000, 0.55000}
   Identifinder-NGH         {-5.00000, 7.30000, 0.59000}
   InSpector 1000 (LaBr3):  { 0.00000, 3.60000, 0.51000}
   InSpector 1000 (NaI)     {-3.00000, 7.60000, 0.58000}
   Interceptor (CZT)        {-0.52000, 2.13000, 0.68000}
   LaBr3 PNNL               { 2.00000, 2.60000, 0.52000}
   NaI 3x3                  {-4.00000, 6.30000, 0.60000}
   Pager-S                  { 0.00000, 8.50000, 0.50000}
   Rad-Pack                 { 0.00000, 8.50000, 0.50000}
   RadSeeker (LaBr3 or NaI) { 4.99023, 2.99883, 0.56174}
   Raider                   { 2.40000, 2.50000, 0.20000}
   Ranger                   {-6.00000, 7.80000, 0.50000}
   Sam-935                  {-6.00000, 7.30000, 0.67000}
   Sam-EagleLaBr3           { 7.00000, 2.60000, 0.52000}
   Sam-Eagle Nai            {-3.00000, 7.00000, 0.65000}
   Spir-ID                  { 7.34303, 3.05497, 0.65990}
   Spir-ID NaI              { 7.44000, 7.61000, 0.63000}
   */
  
  double a_initial, b_initial, c_initial, d_initial;
  double lowerA, upperA, lowerB, upperB, lowerC, upperC;
  
  switch( fnctnlForm )
  {
    case kGadrasResolutionFcn:
    {
      if( highres )
      {
        lowerA = 0.75*1.0;     //Falcon 5000 is 1.0
        upperA = 2.0*1.77;     //Detective-DX 1.77
        lowerB = 0.75*0.20028; //Falcon 5000 is 0.20028
        upperB = 2.0*0.27759;  //Detective-EX is 0.27759
        lowerC = 0.75*0.31;    //Detective-EX
        upperC = 1.5*0.57223;  //Detective-DX
      }else
      {
        lowerA = 1.5*-7.0;     //GR135
        upperA = 1.5*7.44000;  //Spir-ID NaI
        lowerB = 0.75*2.13000; //Interceptor (CZT)
        upperB = 1.5*8.50000;  //Rad-Pack
        lowerC = 0.75*0.20000; //Raider
        upperC = 1.25*0.70000; //FieldSec:
      }
      
      if( answer.size() == 3 )
      { //caller has provided some default values, lets use them
        a_initial = static_cast<double>( answer[0] );
        b_initial = static_cast<double>( answer[1] );
        c_initial = static_cast<double>( answer[2] );
      }else
      {
        if( highres )
        {
          a_initial = 0.5*(1.0 + 1.77);
          b_initial = 0.5*(0.20028 + 0.27759);
          c_initial = 0.5*(0.31 + 0.57223);
        }else
        {
          a_initial = 0.0;//0.5*(-7.0 + 7.44000);
          b_initial = 0.5*(2.13000 + 8.50000);
          c_initial = 0.5*(0.20000 + 0.70000);
        }
      }//if( answer.size() == 3 ) / else
      
      break;
    }//case kGadrasResolutionFcn:
      
    case kSqrtPolynomial:
    {
      if( highres )
      {
        a_initial = 2.0;  //-2.68,  3.96,   2.54,   2.45,    3.3,    2.2,   4.63, 3.8
        b_initial = 0.05;  //0.0579, 0.0488, 0.0697, 0.0459, 0.0669, 0.0605, 0.1, 0.137
        lowerA = -5.0;
        upperA = 7.0;
        lowerB = 0.02;
        upperB = 0.2;
      }else
      {
        //Based on very few detectors checked
        a_initial = -7.0; //-7.76, -7.85, -7.77, 4.78
        b_initial = 2.0;  //1.87, 1.46,  2.1, 2.62
        lowerA = -10.0;
        upperA = 5.0;
        lowerB = 0.5;
        upperB = 4.0;
      }//if( highres ) / else

      c_initial = 0.0;
      d_initial = 0.0;
      
      break;
    }//case kSqrtPolynomial:
      
    case kNumResolutionFnctForm:
      //We should never make it here if the detector FWHM functional form has
      //  not been explicitly set, so throw and error if we got here
      throw runtime_error( "DetectorPeakResponse::performResolutionFit(...):"
                          " invalid ResolutionFnctForm" );
      break;
  }//switch( fnctnlForm )
  
  
  //For high res we should only have to to this one fit
  //But for lowres detectors we should consider the cases of
  //  A==0, A<0, and A>0 when fitting for kGadrasResolutionFcn
  
  DetectorResolutionFitness fitness( *peaks, fnctnlForm );
  
  ROOT::Minuit2::MnUserParameters inputPrams;
  switch( fnctnlForm )
  {
    case kGadrasResolutionFcn:
      if( peaks->size() == 1 )
      {
        inputPrams.Add( "A", a_initial );
        inputPrams.Add( "B", b_initial, 0.1*(upperB-lowerB), lowerB, upperB );
        inputPrams.Add( "C", c_initial );
      }else if( peaks->size() == 2 )
      {
        inputPrams.Add( "A", a_initial );
        inputPrams.Add( "B", b_initial, 0.1*(upperB-lowerB), lowerB, upperB );
        inputPrams.Add( "C", c_initial, 0.1*(upperC-lowerC), lowerC, upperC );
      }else
      {
        inputPrams.Add( "A", a_initial, 0.1*(upperA-lowerA), lowerA, upperA );
        inputPrams.Add( "B", b_initial, 0.1*(upperB-lowerB), lowerB, upperB );
        inputPrams.Add( "C", c_initial, 0.1*(upperC-lowerC), lowerC, upperC );
      }//if( all_peaks.size() == 1 )
      break;
      
    case kSqrtPolynomial:
    {
      if( peaks->size() == 1 || peaks->size() == 2 )
      {
        inputPrams.Add( "A", a_initial );
        inputPrams.Add( "B", b_initial, 0.1*(upperB-lowerB) /*, lowerB, upperB*/ );
      }else
      {
        inputPrams.Add( "A", a_initial, 0.1*(upperA-lowerA), lowerA, upperA );
        inputPrams.Add( "B", b_initial, 0.1*(upperB-lowerB), lowerB, upperB );
      }
      break;
    }//case kSqrtPolynomial:
      
    case kNumResolutionFnctForm:
      break;
  }//switch
  
  
  ROOT::Minuit2::MnUserParameterState inputParamState( inputPrams );
  ROOT::Minuit2::MnStrategy strategy( 2 ); //0 low, 1 medium, >=2 high
  ROOT::Minuit2::MnMinimize fitter( fitness, inputParamState, strategy );

  double tolerance = 1.0;
  unsigned int maxFcnCall = 50000;
  ROOT::Minuit2::FunctionMinimum minimum = fitter( maxFcnCall, tolerance );
  const ROOT::Minuit2::MnUserParameters fitParams = fitter.Parameters();
  
  cerr << "FWHM final chi2=" << fitness.DoEval( fitParams.Params() ) << endl;
  
  if( !minimum.IsValid() )
  {
    stringstream msg;
    msg  << BOOST_CURRENT_FUNCTION
         << ": status is not valid"
         << "\n\tHasMadePosDefCovar: " << minimum.HasMadePosDefCovar()
         << "\n\tHasAccurateCovar: " << minimum.HasAccurateCovar()
         << "\n\tHasReachedCallLimit: " << minimum.HasReachedCallLimit()
         << "\n\tHasValidCovariance: " << minimum.HasValidCovariance()
         << "\n\tHasValidParameters: " << minimum.HasValidParameters()
         << "\n\tIsAboveMaxEdm: " << minimum.IsAboveMaxEdm()
         << endl;
    if( minimum.IsAboveMaxEdm() )
      msg << "\t\tEDM=" << minimum.Edm() << endl;
    cerr << endl << msg.str() << endl;
    throw std::runtime_error( msg.str() );
  }//if( !minimum.IsValid() )
  
  for( size_t i = 0; i < fitParams.Params().size(); ++i )
    cout << "\tFit Par" << i << "=" << fitParams.Params()[i] << endl;
  
  for( const double p : fitParams.Params() )
    answer.push_back( static_cast<float>(p) );
}//std::vector<float> performResolutionFit(...)


