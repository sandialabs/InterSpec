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
#include <ctime>
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
#include <iostream>
#include <iomanip>

#include <boost/tokenizer.hpp>
#include <boost/functional/hash.hpp>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"

#include "mpParser.h"

#include "Wt/Utils"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "InterSpec/AppUtils.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"


using namespace std;
using SpecUtils::Measurement;

const int DetectorPeakResponse::sm_xmlSerializationVersion = 2;

namespace
{
  /** Mutex to protect the "generic" detectors, returned by like `DetectorPeakResponse::getGenericHPGeDetector()` */
  std::mutex s_generic_det_mutex;
  
  
  
  //calcA(...) is for use from DetectorPeakResponse::akimaInterpolate(...).
  //  This function does not to any error checking that the input is valid.
  float calcA( const size_t i,
              const std::vector<DetectorPeakResponse::EnergyEfficiencyPair> &xy )
  {
    assert( i >= 2 );
    assert( (i + 2) < xy.size() );
    
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
  
  
  //removeOutlyingWidthPeaks(...): removes peaks whos width does not agree
  //  well with the functional form passed in.  Returns surviving peaks.
  DetectorPeakResponse::PeakInput_t removeOutlyingWidthPeaks( DetectorPeakResponse::PeakInput_t peaks,
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
      double predicted_sigma = DetectorPeakResponse::peakResolutionSigma( peak->mean(), fnctnlForm, coefficients );
      if( IsNan(predicted_sigma) || IsInf(predicted_sigma) )
        predicted_sigma = 0.0;
      
      const double chi2 = MakeDrfFit::peak_width_chi2( predicted_sigma, *peak );
      
      mean_weight += chi2/npeaks;
      weights.push_back( chi2 );
    }//for( const EnergySigma &es : m_energies_and_sigmas )
    
    vector<size_t> indices( npeaks );
    for( size_t i = 0; i < npeaks; ++i )
      indices[i] = i;
    
    std::sort( begin(indices), end(indices), [&weights](const size_t &a, const size_t &b){
      return weights[b] < weights[a];
    } );
    
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

  //Allow the file to be comma or tab delimited, but if an individual field
  //  contains a comma or tab, then the field must be quoted by a double
  //  quote.  Note that if you just copy cells from Microsoft Excel, that
  //  contain a comma, and then past into a text editor, fields with a comma
  //  are not quoted.
  void split_escaped_csv( vector<string> &fields, const string &line )
  {
    fields.clear();
  
    typedef boost::tokenizer<boost::escaped_list_separator<char> > Tokeniser;
    boost::escaped_list_separator<char> separator("\\",",\t", "\"");
    Tokeniser t( line, separator );
    for( Tokeniser::iterator it = t.begin(); it != t.end(); ++it )
      fields.push_back(*it);
  }//void split_escaped_csv(...)
  
  string to_url_flt_array( const std::vector<float> &vals, const size_t num_sig_fig )
  {
    string answer;
    for( size_t i = 0; i < vals.size(); ++i )
      answer += (i ? "*" : "") + SpecUtils::printCompact( vals[i], num_sig_fig );
    return answer;
  }//string to_url_flt_array

  vector<float> from_url_flt_array( const string &val )
  {
    vector<float> answer;
    vector<string> parts;
    SpecUtils::split( parts, val, "*" );
    for( const string &v : parts )
    {
      float vf;
      if( (stringstream(v) >> vf) )
      {
        answer.push_back( vf );
      }else
      {
        assert( 0 );
        cerr << "Failed in from_url_flt_array: '" << val << "'" << endl;
        return {};
      }
    }//for( const string &v : parts )
    
    return answer;
  }//from_url_flt_array(...)

  std::string url_encode( const std::string& url, const std::string &not_allowed, const bool qr_ascii_only )
  {
    auto to_hex = [](int n) -> char {
      return "0123456789ABCDEF"[(n & 0xF)];
    };
    
    const std::string unsafe_chars = " $&+,:;=?@'\"<>#%{}|\\^~[]`/";
    const string qr_ascii_allowed = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ $%*+-./:";
    
    std::stringstream result;
    
    for( size_t i = 0; i < url.size(); ++i )
    {
      unsigned char c = (unsigned char)url[i];
      
      bool allowable = ((c <= 31) || (c >= 127) || (unsafe_chars.find(c) != std::string::npos));
      
      if( qr_ascii_only && allowable )
        allowable = (qr_ascii_allowed.find(url[i]) != string::npos);
      
      if( allowable )
        result << '%' << to_hex(c >> 4) << to_hex(c);
      else
        result << (char)c;
    }//for( size_t i = 0; i < url.size(); ++i )
    
    return result.str();
  }//url_encode(...)
}//namespace


/** Struct that a string formated equation for evaluation of efficiency.
 
Three expression evaluators have been tried; all have benefits and drawbacks.
  CLHEP/Evaluator: only added 21 kb to final binary size, but slowest option
  exprtk.hpp: added 5.3 Mb to final binary size, but 30x faster than CLHEP
  muparserx: adds 430 kb to binary size, and about 10x faster than CLHEP

 Since muparserx is what TerminalWidget uses, it was decided on 20180408 to
 only use muparserx.

100000 evaluations of:
  exp(-343.6330974237 + 269.1023287277*log(x) + -83.8077567526*log(x)^2
  + 12.9980559362*log(x)^3 + -1.0068649823*log(x)^4 + 0.0311640084*log(x)^5)
  muparser x took:  cpu=0.143294s, wall=0.14341s  (1.4 us/eval)
  evaluator x took: cpu=1.10172s, wall=1.10209s   (11 us/eval)
*/
struct FormulaWrapper
{
  /** Constructor that takes an equation as a string, and creates a callable
   object to evaluate that equation.
   
   \param fcnstr The function to use.  Ex: "exp(-1.2 + 3*lox(x) + ...)"
   \param isMeV Only used to determine which energy value to use to test if the
   function is valid or not;  If isMeV is true, uses 0.1, else 100.
   */
  FormulaWrapper( const std::string &fcnstr, const bool isMev );
  ~FormulaWrapper();
  
  float efficiency( const float x );
  double operator()( const float x );
  
  /** Finds the variable the user most likely intended to be the energy variable
   for detector response functions, by looking for arguments inside
   paranthesis.
   Returned answer is always lower case.  Spaces, tabs, and newlines are all
   removed from input string before searching.
   Returns "x" if it doesnt find any other candidates.
   
   e.x., assert(find_variable_name("exp(2.0*log(x) - 1.1*log(x)^2)") == "x");
   assert(find_variable_name("exp(2.0*log(E) - 1.1*log(E)^2)") == "e");
   assert(find_variable_name("3.2*exp(energy)") == "energy");
   assert(find_variable_name("") == "x");
   assert(find_variable_name("1-energy") == "x"); //no (...), so defaults to "x"
   */
  static std::string find_variable_name( std::string eqn );
  
  std::mutex m_mutex;
  std::string m_fcnstr;
  std::string m_var_name;
  
  std::unique_ptr<mup::ParserX> m_parser;
  std::unique_ptr<mup::Value> m_value;
};//struct FormulaWrapper

  
FormulaWrapper::FormulaWrapper( const std::string &fcnstr, const bool isMev )
  : m_fcnstr( fcnstr ), m_var_name( "x" )
{
  SpecUtils::to_lower_ascii( m_fcnstr );
  SpecUtils::erase_any_character( m_fcnstr, "\t\n\r" );
  
  
  m_parser.reset( new mup::ParserX() );
  m_value.reset( new mup::Value(isMev ? 0.1 : 100.0) );
  
  if( m_fcnstr.empty() )
    m_fcnstr = "1";
  
  
  
  try
  {
    m_var_name = "x";
    mup::ParserX trialparser;
    trialparser.DefineVar( m_var_name, mup::Variable(m_value.get()) );
    trialparser.SetExpr( m_fcnstr );
    trialparser.Eval();  //Trigger evaluation; if we get an exception, we'll try somethign besides "x".
  }catch( mup::ParserError &e )
  {
    const mup::EErrorCodes code = e.GetCode();
    
    if( code == mup::EErrorCodes::ecUNASSIGNABLE_TOKEN )
      m_var_name = e.GetToken();
    else
      m_var_name = find_variable_name( m_fcnstr );
  }//try catch to figure out if 'x' works for varaiable name, or we need something else
  
  
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
  SpecUtils::to_lower_ascii( eqn );
  
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
      
      bool maybeOtherStuff = false;
      //Really weak test to make sure parenthesis only contains the variable name
      //  ToDo: improve this test!  Wouldnt catc should be y in 'pow(y,2)', etc.
      if( substr.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_") != string::npos )
        maybeOtherStuff = true;
      
      if( open_pos == 0 )  //Check equation isnt starting with this parenthesis
        maybeOtherStuff = true;
      else if( !isalpha(eqn[open_pos-1]) )  //check this parenthesis is somethign lile 'log(y)', or 'sqrt(y)', and not '1/(x)^2'.  ... actualy sure this test is necassary...
        maybeOtherStuff = true;
      
      if( !maybeOtherStuff )
      {
        auto iter = argcounts.find(substr);
        if( iter == end(argcounts) )
          argcounts[substr] = 1;
        else
          iter->second += 1;
      }//if( !maybeOtherStuff )
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

const std::string &DetectorPeakResponse::det_eff_geom_type_postfix( const DetectorPeakResponse::EffGeometryType type )
{
  static const string s_empty{}, s_cm2{"/cm2"}, s_m2{"/m2"}, s_gram{"/g"};
  switch( type )
  {
    case DetectorPeakResponse::EffGeometryType::FarField:
    case DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct:
      return s_empty;
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2:
      return s_cm2;
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2:
      return s_m2;
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram:
      return s_gram;
  }//switch( m_det_type )
  assert( 0 );
  return s_empty;
}//string det_eff_geom_type_postfix( DetectorPeakResponse::EffGeometryType )


DetectorPeakResponse::DetectorPeakResponse()
  : m_name( "DetectorPeakResponse" ),
    m_description( "" ),
    m_detectorDiameter( 0.0f ),
    m_efficiencyEnergyUnits( static_cast<float>(PhysicalUnits::keV) ),
    m_resolutionForm( DetectorPeakResponse::kNumResolutionFnctForm ),
    m_efficiencySource( DrfSource::UnknownDrfSource ),
    m_efficiencyForm( DetectorPeakResponse::kNumEfficiencyFnctForms ),
    m_hash( 0 ),
    m_parentHash( 0 ),
    m_flags( 0 ),
    m_lowerEnergy( 0.0 ),
    m_upperEnergy( 0.0 ),
    m_createdUtc( 0 ),
    m_lastUsedUtc( 0 ),
    m_geomType(EffGeometryType::FarField)
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
    m_efficiencySource( DrfSource::UnknownDrfSource ),
    m_efficiencyForm( DetectorPeakResponse::kNumEfficiencyFnctForms ),
    m_hash( 0 ),
    m_parentHash( 0 ),
    m_flags( 0x0 ),
    m_lowerEnergy( 0.0 ),
    m_upperEnergy( 0.0 ),
    m_createdUtc( 0 ),
    m_lastUsedUtc( 0 ),
    m_geomType(EffGeometryType::FarField)
{
}//DetectorPeakResponse( const std::string &name, const std::string &descrip )


DetectorPeakResponse::~DetectorPeakResponse()
{
}//~DetectorPeakResponse()


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
  
  for( const float val : m_expOfLogPowerSeriesUncerts )
    boost::hash_combine( seed, val );
  
  //For legacy DRFs (pre 1.0.4 / 20190527) that dont have DRF source info,
  //  lower/upper energies, and m_flags, dont change the hash.
  if( m_lowerEnergy != 0.0 || m_upperEnergy != 0.0 )
  {
    boost::hash_combine( seed, m_lowerEnergy );
    boost::hash_combine( seed, m_upperEnergy );
  }
  
  //Dont hash based on m_createdUtc or m_lastUsed
  
  if( m_geomType != EffGeometryType::FarField )
    boost::hash_combine( seed, m_geomType );
  
  m_hash = seed;
}//void computeHash()


//isValid(): Tells wether or not this DetectorPeakResponse has been properly
//  defined and can be used
bool DetectorPeakResponse::isValid() const
{
  switch( m_efficiencyForm )
  {
    case kEnergyEfficiencyPairs:
      return (m_energyEfficiencies.size() >= 2);
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
    case kSqrtEnergyPlusInverse:
      return (m_resolutionCoeffs.size() == 3);
    case kConstantPlusSqrtEnergy:
      return (m_resolutionCoeffs.size() == 2);
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
  m_resolutionCoeffs.clear();
  m_energyEfficiencies.clear();
  m_hash = m_parentHash = 0;
  m_flags = 0;
  m_lowerEnergy = m_upperEnergy = 0.0;
  m_efficiencySource = DrfSource::UnknownDrfSource;
  m_createdUtc = m_lastUsedUtc = 0;
  m_geomType = EffGeometryType::FarField;
  
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
          && m_parentHash==rhs.m_parentHash
          && m_lowerEnergy==rhs.m_lowerEnergy
          && m_upperEnergy==rhs.m_upperEnergy
          //m_createdUtc
          //m_lastUsedUtc
          && m_geomType==rhs.m_geomType
          );
}//operator==


DetectorPeakResponse::DrfSource DetectorPeakResponse::drfSource() const
{
  return m_efficiencySource;
}


float DetectorPeakResponse::efficiencyEnergyUnits() const
{
  return m_efficiencyEnergyUnits;
}

void DetectorPeakResponse::setDrfSource( const DetectorPeakResponse::DrfSource source )
{
  if( m_efficiencySource == source )
    return;
  
  m_efficiencySource = source;
  computeHash();
}//void setDrfSource( const DrfSource source );


void DetectorPeakResponse::setEnergyRange( const double lower, const double upper )
{
  m_lowerEnergy = lower;
  m_upperEnergy = upper;
  computeHash();
}//void setEnergyRange( const double lower, const double upper )


double DetectorPeakResponse::lowerEnergy() const
{
  return m_lowerEnergy;
}


double DetectorPeakResponse::upperEnergy() const
{
  return m_upperEnergy;
}


bool DetectorPeakResponse::isFixedGeometry() const
{
  return (m_geomType != EffGeometryType::FarField);
}//bool isFixedGeometry() const;


DetectorPeakResponse::EffGeometryType DetectorPeakResponse::geometryType() const
{
  return m_geomType;
}


void DetectorPeakResponse::updateLastUsedTimeToNow()
{
  m_lastUsedUtc = std::time(nullptr);
}//void updateLastUsedTimeToNow()


void DetectorPeakResponse::fromEnergyEfficiencyCsv( std::istream &input,
                                    const float detectorDiameter,
                                    const float energyUnits,
                                    //const bool is_fixed_geometry
                                    const EffGeometryType geometry_type )
{
  if( // !is_fixed_geometry
     (geometry_type == EffGeometryType::FarField)
     && ((detectorDiameter <= 0.0) || IsInf(detectorDiameter) || IsNan(detectorDiameter)) )
    throw runtime_error( "Detector diameter must be greater than 0.0" );
  
  if( energyUnits <= 0.0 || IsInf(energyUnits) || IsNan(energyUnits) )
    throw runtime_error( "Energy units must be greater than 0.0" );

  
  bool gotline = false;
  int linen = 0;
  string line;
  std::vector<EnergyEfficiencyPair> energyEfficiencies;
  
  while( SpecUtils::safe_get_line( input, line ) )
  {
    ++linen;
    SpecUtils::trim( line );

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
    const bool ok = SpecUtils::split_to_floats( line.c_str(),
                                                       line.size(), fields );
    
    if( (fields.size() < 2) || !ok )
    {
      if( !gotline )
        continue;
      
      throw runtime_error( "Invalid efficiency file.  After the CSV file header"
                           ", all lines must have at least two numbers: "
                           "energy and % efficiency. Further comma separated "
                           "numbers are allowed, but no other text" );
    }//if( fields.size() < 2 )
    
    EnergyEfficiencyPair point;
    point.energy = fields[0] * PhysicalUnits::keV;
    point.efficiency = fields[1] / 100.0; //files have it listed as a percentage
    
    gotline = true;
    energyEfficiencies.push_back( point );
    
    if( linen > 5000 || energyEfficiencies.size() > 3000 )
      throw runtime_error( "Efficiency CSV files can have a max of 5000 total lines or 3000 energy efficiency pairs" );
    
  }//while( SpecUtils::safe_get_line(input) )

  
  if( energyEfficiencies.size() < 2 )
    throw runtime_error( "File was not a valid efficiency file, need at least"
                        " two efficiency points." );
  
  //Lets check if the numbers were strickly increasing or strickly decreasing
  const bool increasing = (energyEfficiencies[1].energy > energyEfficiencies[0].energy);
  for( size_t i = 2; i < energyEfficiencies.size(); ++i )
  {
    //We'll let two be energied be the exact same in a row
    const float laste = energyEfficiencies[i-1].energy;
    const float thise = energyEfficiencies[i].energy;
    
    if( IsNan(laste) || IsNan(thise)
        || IsInf(laste) || IsInf(thise) )
      throw runtime_error( "NaN or Inf energy detected from CSV file" );
    
    if( laste<0.0f || thise<0.0f )
      throw runtime_error( "Energy can not be less than zero in CSV efficiecny file" );
    
    if( energyEfficiencies[i-1].efficiency < 0.0f
       || energyEfficiencies[i].efficiency < 0.0f
       || energyEfficiencies[i-1].efficiency > 1.0f
       || energyEfficiencies[i].efficiency > 1.0f )
      throw runtime_error( "Intrinsic efficiency can not be less than zero or"
                           " greater than 100 in CSV efficiency file" );
    
    if( fabs(laste - thise) < 1.0E-5*std::max(fabs(thise),fabs(laste)) )
      continue;
    
    const bool bigger = (thise > laste);
    if( increasing != bigger )
      throw runtime_error( "Energies in efficiency CSV must be strickly "
                           "increasing or strickly decreasing" );
  }//for( size_t i = 2; i < energyEfficiencies.size(); ++i )
  

  if( !increasing )
    std::reverse( energyEfficiencies.begin(), energyEfficiencies.end() );
  
  m_detectorDiameter = detectorDiameter;
  m_energyEfficiencies = energyEfficiencies;
  m_efficiencyEnergyUnits = energyUnits;
  
  m_efficiencySource = DrfSource::UserImportedIntrisicEfficiencyDrf;
  m_efficiencyForm = kEnergyEfficiencyPairs;
  
  m_flags = 0;
  
  m_lowerEnergy = m_energyEfficiencies.front().energy;
  m_upperEnergy = m_energyEfficiencies.back().energy;
  
  m_lastUsedUtc = m_createdUtc = std::time(nullptr);
  m_geomType = geometry_type;
  
  computeHash();
}//void fromEnergyEfficiencyCsv(...)



void DetectorPeakResponse::setIntrinsicEfficiencyFormula( const string &fcnstr,
                                                          const float diameter,
                                                          const float energyUnits,
                                                          const float lowerEnergy,
                                                          const float upperEnergy,
                                                          //const bool fixedGeometry
                                                          const EffGeometryType geometry_type )
{
  const bool isMeV = (energyUnits > 10.f);
  std::shared_ptr<FormulaWrapper > expression
                                = std::make_shared<FormulaWrapper>( fcnstr, isMeV );
  
  m_efficiencyForm = kFunctialEfficienyForm;
  m_efficiencyFormula = fcnstr;
  m_efficiencyFcn = boost::bind( &FormulaWrapper::efficiency, expression,
                                boost::placeholders::_1  );
  m_detectorDiameter = diameter;
  m_efficiencyEnergyUnits = energyUnits;
  
  m_flags = 0;
  
  m_lowerEnergy = lowerEnergy;
  m_upperEnergy = upperEnergy;
  
  m_efficiencySource = DrfSource::UserSpecifiedFormulaDrf;
  
  m_lastUsedUtc = m_createdUtc = std::time(nullptr);
  m_geomType = geometry_type;
  
  computeHash();
}//void setIntrinsicEfficiencyFormula( const std::string &fcn )


function<float(float)> DetectorPeakResponse::makeEfficiencyFunctionFromFormula(
                                                                     const std::string &formula,
                                                                     const bool isMeV )
{
  shared_ptr<FormulaWrapper> expression = make_shared<FormulaWrapper>( formula, isMeV );
  return [expression](float a) -> float { return expression->efficiency(a); };
}


void DetectorPeakResponse::fromGadrasDefinition( std::istream &csvFile,
                                                 std::istream &detDatFile )
{
  float detWidth = 0.0, heightToWidth = 0.0;
  float lowerEnergy = 0.0, upperEnergy = 0.0;

  m_efficiencyEnergyUnits = static_cast<float>( PhysicalUnits::keV );
  m_resolutionCoeffs.clear();
  m_resolutionForm = kGadrasResolutionFcn;
  m_resolutionCoeffs.resize( 3, 0.0 );

  
  // First decide if old style Detector.dat, or the newer XML Detector.dat
  const istream::pos_type orig_pos = detDatFile.tellg();
  
  string line;
  if( !SpecUtils::safe_get_line( detDatFile, line ) )
    throw std::runtime_error( "Couldnt read first line of Detector.dat" );
  
  const bool is_xml = (line.find("xml") != string::npos);
  
  if( is_xml )
  {
    try
    {
      // A helper function to grab float values from Detector.dat xml.
      auto get_float_value = []( const rapidxml::xml_node<char> * const parent, const string &name ) -> float {
        const auto target = parent->first_node( name.c_str(), name.size() );
        const auto value = SpecUtils::xml_first_node( target, "value" );
        const string val_str = SpecUtils::xml_value_str( value );
        
        float val;
        //will fail if val_str empty (which we want), but will parse strings like " 19.2 as", which
        //  maybe its debatable if we want this behavior, but I dont think this matters much.
        if( !(stringstream(val_str) >> val) )
          throw runtime_error( "Missing <" + name + "> node." );
        return val;
      };//get_float_value lambda
      
      
      detDatFile.seekg(orig_pos); //probably not really needed
      rapidxml::file<char> input_file( detDatFile );
      
      rapidxml::xml_document<char> doc;
      doc.parse<rapidxml::parse_trim_whitespace>( input_file.data() );
      
      const auto gamma_detector = doc.first_node( "gamma_detector" );
      if( !gamma_detector )
        throw runtime_error( "Missing <gamma_detector> node." );
      
      lowerEnergy = get_float_value( gamma_detector, "weight_range_lower" );  //template error / min energy, keV
      upperEnergy = get_float_value( gamma_detector, "weight_range_upper" );  //chi-square / max energy, keV
      
      const auto dimensions = XML_FIRST_NODE( gamma_detector, "dimensions" );
      if( !dimensions )
        throw runtime_error( "Missing <dimensions> node." );
      
      detWidth = get_float_value( dimensions, "width" );
      heightToWidth = get_float_value( dimensions, "height_to_width_ratio" );
      
      const auto peak_shape = XML_FIRST_NODE( gamma_detector, "peak_shape" );
      if( !peak_shape )
        throw runtime_error( "Missing <peak_shape> node." );
      
      m_resolutionCoeffs[0] = get_float_value( peak_shape, "fwhm_offset" );
      m_resolutionCoeffs[1] = get_float_value( peak_shape, "fwhm_at_661keV" );
      m_resolutionCoeffs[2] = get_float_value( peak_shape, "fwhm_power" );
    }catch( std::exception &e )
    {
      throw std::runtime_error( "Failed to read XML Detector.dat: " + std::string(e.what()) );
    }//try / catch
  }else
  {
    
    do  //We've already got the first line
    {
      vector<string> parts;
      SpecUtils::trim( line );
      if( line.empty() || !isdigit(line[0]) )
        continue;
      
      SpecUtils::split( parts, line, " \t" );
      
      try
      {
        const int parnum = std::stoi( parts.at(0) );
        const float value = static_cast<float>( std::stod( parts.at(1) ) );
        //      const int value2 = std::stoi( parts.at(2) );
        //      const string descrip = parts.at(3);
        switch( parnum )
        {
          case 6:  m_resolutionCoeffs[0] = value; break;
          case 7:  m_resolutionCoeffs[1] = value; break;
          case 8:  m_resolutionCoeffs[2] = value; break;
          case 11: detWidth = value;              break;
          case 12: heightToWidth = value;         break;
          //case 35: lowerEnergy = value;           break;  //LLD(keV)
          //case ??: upperEnergy = value;           break;  //chi-square / max energy, keV
        }//switch( parnum )
      }catch(...)
      {
        cerr << "\nError reading line \"" << line << "\"" << endl;
        continue;
      }
    }while( SpecUtils::safe_get_line( detDatFile, line ) );
  }//if( is_xml ) / else
  
  
  if( detWidth<=0.0 || heightToWidth<=0.0 )
    throw runtime_error( "Couldnt find detector dimensions in the Detector.dat file" );
  
  const float surfaceArea = detWidth * detWidth * heightToWidth;
  //const float diam = (4.0f/3.14159265359f) * sqrt(surfaceArea) * static_cast<float>(PhysicalUnits::cm);
  const float diam = 2.0f*sqrt(surfaceArea/3.14159265359f) * static_cast<float>(PhysicalUnits::cm);
  
  //const bool fixed_geom = false;
  const EffGeometryType geometry_type = EffGeometryType::FarField;
  fromEnergyEfficiencyCsv( csvFile, diam, static_cast<float>(PhysicalUnits::keV), geometry_type );
  
  m_flags = 0;
  
  m_lowerEnergy = lowerEnergy;
  m_upperEnergy = upperEnergy;

  m_efficiencySource = DrfSource::UserAddedGadrasDrf;
  
  m_lastUsedUtc = m_createdUtc = std::time(nullptr);
  m_geomType = geometry_type;
  
  computeHash();
}//void fromGadrasDefinition(...)


void DetectorPeakResponse::fromGadrasDirectory( const std::string &dir )
{
  string eff_file, dat_file;
  
  //We will go through a little bit of trouble to not be case-sensitive
  //  For efficiency since some of the Detector.dat files that come with
  //  GADRAS are all uppercase.
  const vector<string> files = SpecUtils::ls_files_in_directory( dir );
  for( size_t i = 0; (eff_file.empty() || dat_file.empty()) && i < files.size(); ++i )
  {
    const auto &file = files[i];
    auto filename = SpecUtils::filename(file);
    if( SpecUtils::iequals_ascii( filename, "Efficiency.csv") )
      eff_file = SpecUtils::append_path( dir, filename );
    else if( SpecUtils::iequals_ascii( filename, "Detector.dat") )
      dat_file = SpecUtils::append_path( dir, filename );
  }
  
  if( dat_file.empty() )
    throw runtime_error( "Directory '" + dir + "' did not contain a Detector.dat file." );
  
  if( eff_file.empty() )
    throw runtime_error( "Directory '" + dir + "' did not contain a Efficiency.csv file." );
  
#ifdef _WIN32
  const std::wstring weff_file = SpecUtils::convert_from_utf8_to_utf16(eff_file);
  const std::wstring wdat_file = SpecUtils::convert_from_utf8_to_utf16(dat_file);
  ifstream csv( weff_file.c_str(), ios_base::binary|ios_base::in );
  ifstream datFile( wdat_file.c_str(), ios_base::binary|ios_base::in );
#else
  ifstream csv( eff_file.c_str(), ios_base::binary|ios_base::in );
  ifstream datFile( dat_file.c_str(), ios_base::binary|ios_base::in );
#endif
  
  if( !csv.is_open() )
    throw runtime_error( "Could not open file '" + eff_file + "'." );
  
  if( !datFile.is_open() )
    throw runtime_error( "Could not open file '" + dat_file + "'." );
  
  fromGadrasDefinition( csv, datFile );
}//void fromGadrasDirectory( const std::string &dir )


void DetectorPeakResponse::fromExpOfLogPowerSeriesAbsEff(
                                   const std::vector<float> &coefs,
                                   const std::vector<float> &uncerts,
                                   const float charactDist,
                                   const float det_diam,
                                   const float equationEnergyUnits,
                                   const float lowerEnergy,
                                   const float upperEnergy,
                                   //const bool fixedGeometry
                                   const EffGeometryType geometry_type
)
{
  if( coefs.empty() )
    throw runtime_error( "fromExpOfLogPowerSeriesAbsEff(...): invalid input" );
  
  if( !uncerts.empty() && (uncerts.size() != coefs.size()) )
    throw runtime_error( "DetectorPeakResponse::fromExpOfLogPowerSeriesAbsEff: uncertainties"
                        " must either be empty, or same size as coefficients." );
  
  m_energyEfficiencies.clear();
  m_expOfLogPowerSeriesCoeffs.clear();
  m_expOfLogPowerSeriesUncerts.clear();
  
  m_efficiencyFcn = std::function<float(float)>();
  m_efficiencyFormula.clear();
  
  m_efficiencyForm = kExpOfLogPowerSeries;
  m_efficiencyEnergyUnits = equationEnergyUnits;
  m_detectorDiameter = det_diam;
  m_expOfLogPowerSeriesCoeffs = coefs;
  m_expOfLogPowerSeriesUncerts = uncerts;
  
  //now we need to account for characterizationDist
  if( (charactDist > 0.0f) && (geometry_type == EffGeometryType::FarField) )
  {
    //x^n * x^m = x^(n+m)
    const float gfactor = fractionalSolidAngle( det_diam, charactDist );
    // TODO: Note that this does not account for attenuation in the air, which is something like 2.2% at 59 keV
    m_expOfLogPowerSeriesCoeffs[0] += log( 1.0f/gfactor );
  }//if( characterizationDist > 0.0f )
  
  m_flags = 0;
  
  m_lowerEnergy = lowerEnergy;
  m_upperEnergy = upperEnergy;
  
  m_efficiencySource = DrfSource::UserAddedRelativeEfficiencyDrf;
  
  m_lastUsedUtc = m_createdUtc = std::time(nullptr);
  m_geomType = geometry_type;
  
  computeHash();
}//void fromExpOfLogPowerSeriesAbsEff


std::shared_ptr<DetectorPeakResponse> DetectorPeakResponse::parseSingleCsvLineRelEffDrf( std::string &line )
{
  std::shared_ptr<DetectorPeakResponse> det;
  
  vector<string> fields;
  fields.reserve( 20 );
  
  split_escaped_csv( fields, line );
  
  // If the line is a "URL" encoded DRF, try parsing it.  Right now I'm being a little frugal
  //  about only letting there be 4 fields, and there are probably some more smaller issues to
  //  be taken care of.
  if( (fields.size() == 4) && SpecUtils::iequals_ascii(fields[2], "UrlEncoded") )
  {
    try
    {
      const string unencoded = Wt::Utils::urlDecode( fields[3] );
      return parseFromAppUrl( unencoded );
    } catch( std::exception &e )
    {
      cerr << "Failed to parse 'UrlEncoded' DRF: " << e.what() << endl;
    }
    
    return det;
  }//if( a UrlEncoded DRF )
  
  if( fields.size() < 16 )
    return det;
  
  try
  {
    string name = fields[0] + " (" + fields[1] + ")";
    SpecUtils::ireplace_all( name, "%20", " " );
    SpecUtils::trim( name );

    vector<float> coefs;
    for( int i = 3; i < 11; ++i )
      coefs.push_back( static_cast<float>( std::stod( fields[i] ) ) );
    
    //Get rid of the zero coefficients
    for( size_t i = coefs.size()-1; i > 0; --i )
    {
      if( fabs(coefs[i]) > 1.0E-14 ) //1.0E-14 chosen arbitrarily
      {
        coefs.erase( coefs.begin()+i+1, coefs.end() );
        break;
      }
    }//for( size_t i = coefs.size()-1; i > 0; --i )
    
    const float dist = static_cast<float>( std::stod(fields[14])*PhysicalUnits::cm );
    const float diam = 2.0f*static_cast<float>( std::stod(fields[15])*PhysicalUnits::cm );
    const float eunits = static_cast<float>( PhysicalUnits::MeV );
    
    string description = fields[2] + " - from Relative Eff. File";
    SpecUtils::ireplace_all( description, "%20", " " );
    SpecUtils::trim( description );
    det.reset( new DetectorPeakResponse( name, description ) );
    det->fromExpOfLogPowerSeriesAbsEff( coefs, {}, dist, diam, eunits, 0.0f, 0.0f, EffGeometryType::FarField );
    det->setDrfSource( DetectorPeakResponse::DrfSource::UserAddedRelativeEfficiencyDrf );
  }catch( std::exception &e )
  {
    det.reset();
    cerr << "DetectorPeakResponse::parseSingleCsvLineRelEffDrf() caught: '" << e.what()
         << "', on line='" << line << endl;
  }//try /catch
  
  return det;
}//shared_ptr<DetectorPeakResponse> parseSingleCsvLineRelEffDrf( std::string &line )


void DetectorPeakResponse::parseMultipleRelEffDrfCsv( std::istream &input,
                                                     std::vector<std::string> &credits,
                                      std::vector<std::shared_ptr<DetectorPeakResponse>> &drfs )
{
  credits.clear();
  drfs.clear();
  
  const std::streampos start_pos = input.tellg();
  
  vector<vector<float>> detcoefs;
  string line;
  while( SpecUtils::safe_get_line( input, line, 2048 ) )
  {
    SpecUtils::trim( line );
    
    if( SpecUtils::istarts_with( line, "#credit:") )
      credits.push_back( SpecUtils::trim_copy( line.substr(8) ) );
    
    if( line.empty() || line[0]=='#' )
      continue;
    
    std::shared_ptr<DetectorPeakResponse> det = parseSingleCsvLineRelEffDrf( line );
    
    if( det )
      drfs.push_back( det );
  }//while( SpecUtils::safe_get_line( input, line ) )
  
  if( drfs.empty() )
  {
    // If we didnt read and DRFs, lets be nice and reset file location to where we started
    credits.clear();
    input.seekg( start_pos, ios::beg );
  }
}//parseMultipleRelEffDrfCsv(...)



void DetectorPeakResponse::parseGammaQuantRelEffDrfCsv( std::istream &input,
                                      std::vector<std::shared_ptr<DetectorPeakResponse>> &drfs,
                                      std::vector<std::string> &credits,
                                      std::vector<std::string> &warnings )
{
  drfs.clear();
  warnings.clear();
  
  const std::streampos start_pos = input.tellg();
  
  try
  {
    string line;
    
    // Allow for some empty lines at the top
    while( SpecUtils::safe_get_line( input, line, 128*1024 ) )
    {
      SpecUtils::trim( line );
      
      // Remove any leading non-ascii characters - e.g., UTF-8 BOM
      while( line.size() && (static_cast<unsigned char>(line.front()) > 127) )
        line.erase( begin(line) );
      
      if( SpecUtils::istarts_with( line, "#credit:") )
        credits.push_back( SpecUtils::trim_copy( line.substr(8) ) );
      
      if( line.empty() || (line[0] == '#') )
        continue;
      
      const size_t pos = line.find_first_not_of(" \t,");
      if( pos == string::npos )
        continue;
      
      break;
    }//while( SpecUtils::safe_get_line( input, line, 128*1024 ) )
    
    if( line.empty() )
      throw runtime_error( "No content" );
    
    size_t row_num = 0;
    vector<string> detector_id;
    split_escaped_csv( detector_id, line );
    const size_t num_col = detector_id.size();
    if( (num_col < 2) || !SpecUtils::icontains(detector_id[0], "Detector ID") )
      throw runtime_error( "First column must be row headers, and start with \"Detector ID\"" );
    
    // Define a utility lambda to help get each line
    auto get_row = [num_col, &line, &row_num, &input, &warnings]( const bool is_empty, const string &expected_title )
      -> vector<string> {
      row_num += 1;
      if( !SpecUtils::safe_get_line( input, line, 128*1024 ) )
        throw runtime_error( "Not enough rows - only " + std::to_string(row_num) + " rows." );
      
      // TODO: check if last cell has newline in it, and if so, get the next line, and append to current.
        
      SpecUtils::trim( line );
      vector<string> cols;
      split_escaped_csv( cols, line );
        
      if( !is_empty && (cols.size() != num_col) )
        throw runtime_error( "Inconsistent number of columns on line " + std::to_string(row_num)
                            + " - expected " + std::to_string(num_col) + " and received "
                            + std::to_string(cols.size()) );
        
      SpecUtils::ireplace_all(cols[0], "\n", " ");
      SpecUtils::ireplace_all(cols[0], "\r", " ");
      SpecUtils::ireplace_all(cols[0], "\t", " ");
      while( SpecUtils::icontains(cols[0], "  ") )
        SpecUtils::ireplace_all(cols[0], "  ", " ");
      
      if( !expected_title.empty() && !SpecUtils::icontains(cols[0], expected_title) )
        warnings.push_back( "Warning: row " + std::to_string(row_num)
                           + " doesnt have expected header '" + expected_title + "'" );
        
      return cols;
    };//get_row lambda
    
    
    const vector<string> calibration_geometry = get_row(false, "Calibration Geometry");
    const vector<string> comments = get_row(false, "Comments");
    const vector<string> is_far_field = get_row(false, "Far-Field Point Source Eff Cal?");
    get_row(true, "");
    const vector<string> det_radius = get_row(false, "Detector Radius - a (cm)");
    const vector<string> distance = get_row(false, "Source to Detector Face - d (cm)");
    get_row(false, "Source Radius - r (cm)"); //We dont use this info, yet
    const vector<string> detector_setbacks = get_row(false, "Detector Setback (cm)");  //We dont use this info, yet
    get_row(false, "Detector Area - (cm2)");
    get_row(false, "Source to Crystal Face - d (cm)");
    get_row(false, "Relative Efficiency (%)");
    get_row(true, "");
    get_row(true, "EFFICIENCY CALIBRATION AND CURVE FIT INFORMATION");
    const vector<string> lower_energy = get_row(false, "Lower Energy Cutoff (keV)");
    const vector<string> upper_energy = get_row(false, "Upper Energy Cutoff (keV)");
    get_row(true, "");
    get_row(true, ""); //Actual header:  "y = exp (a  + b(lnx) + c(lnx)2 + d(lnx)3 + e(lnx)4 + f(lnx)5 + g(lnx)6 + h(lnx)7)"
    const vector<string> keV_or_MeV = get_row(false, "keV or MeV");
    const vector<string> coef_a = get_row(false, "Coefficient a");
    const vector<string> coef_b = get_row(false, "Coefficient b");
    const vector<string> coef_c = get_row(false, "Coefficient c");
    const vector<string> coef_d = get_row(false, "Coefficient d");
    const vector<string> coef_e = get_row(false, "Coefficient e");
    const vector<string> coef_f = get_row(false, "Coefficient f");
    vector<string> coef_g, coef_h;
    try
    {
      coef_g = get_row(false, "Coefficient g");
      coef_h = get_row(false, "Coefficient h");
    }catch( std::exception & )
    {
      //G and H arent used that much anyway
    }
    
    for( size_t col = 1; col < num_col; ++col )
    {
      const string &name = detector_id[col];
      if( name.empty() )
        continue;
      
      const string &comment = comments[col];
      const string &cal_geom = calibration_geometry[col];
      
      const string &far_field = is_far_field[col];
      if( far_field.empty() )
        continue;
      
      const bool is_fixed_geom = SpecUtils::iequals_ascii( far_field, "No" );
      const bool is_far_field_val = SpecUtils::iequals_ascii( far_field, "Yes" );
      
      if( !is_fixed_geom && !is_far_field_val )
      {
        warnings.push_back( "Far-Field row of column " + std::to_string(col) 
                           + " should be either Yes or No, was '" + far_field + "'" );
        continue;
      }
      
      double det_rad = 0, src_to_det_face = 0, det_setback = 0, low_energy = 0, up_energy = 0;
      
      const string &det_rad_str = det_radius[col];
      if( !SpecUtils::parse_double(det_rad_str.c_str(), det_rad_str.size(), det_rad) )
      {
        warnings.push_back( "Failed to parse detector radius, '" + det_rad_str
                           + "', to a number, skipping detector '" + name + "'" );
        continue;
      }
      det_rad *= PhysicalUnits::cm;
      if( det_rad < 0.0 )
      {
        warnings.push_back( "Detector radius negative ('" + det_rad_str
                           + "') - skipping detector '" + name + "'" );
        continue;
      }
      
      const string &dist = distance[col];
      if( !SpecUtils::parse_double(dist.c_str(), dist.size(), src_to_det_face) )
      {
        warnings.push_back( "Failed to parse distance to detector, '" + dist
                           + "', to a number, skipping detector '" + name + "'" );
        continue;
      }
      src_to_det_face *= PhysicalUnits::cm;
      
      if( src_to_det_face < 0.0 )
      {
        warnings.push_back( "Src to face distance, '" + dist
                           + "', is negative - skipping detector '" + name + "'" );
        continue;
      }
        
      const string &setback = detector_setbacks[col];
      if( !SpecUtils::parse_double(setback.c_str(), setback.size(), det_setback) )
      {
        warnings.push_back( "Failed to parse detector setback, '" + setback
                           + "', to a number, skipping detector '" + name + "'" );
        continue;
      }
        
      
      const string &low_energy_str = lower_energy[col];
      if( !SpecUtils::parse_double(low_energy_str.c_str(), low_energy_str.size(), low_energy) )
      {
        low_energy = 0;
        warnings.push_back( "Failed to parse lower energy, '" + low_energy_str
                           + "', to a number, for detector '" + name + "'" );
        //continue;
      }
      
      const string &up_energy_str = upper_energy[col];
      if( !SpecUtils::parse_double(up_energy_str.c_str(), up_energy_str.size(), up_energy) )
      {
        up_energy = 0;
        low_energy = 0;
        warnings.push_back( "Failed to parse upper energy, '" + up_energy_str
                           + "', to a number, for detector '" + name + "'" );
        //continue;
      }
        
      const string &keV_or_MeV_str = keV_or_MeV[col];
      const bool is_keV = SpecUtils::iequals_ascii(keV_or_MeV_str, "keV");
      const bool is_MeV = SpecUtils::iequals_ascii(keV_or_MeV_str, "MeV");
      if( !is_keV && !is_MeV )
      {
        warnings.push_back( "The keV/MeV field for detector '" + name + "', is invalid, '" 
                           + keV_or_MeV_str+ "' (must be either keV or MeV_ - skipping detector" );
        continue;
      }
      const double energy_units = is_keV ? PhysicalUnits::keV : PhysicalUnits::MeV;
      
      vector<string> coef_strs;
      coef_strs.push_back( coef_a[col] );
      coef_strs.push_back( coef_b[col] );
      coef_strs.push_back( coef_c[col] );
      coef_strs.push_back( coef_d[col] );
      coef_strs.push_back( coef_e[col] );
      coef_strs.push_back( coef_f[col] );
      if( (coef_g.size() > col) && (!coef_g[col].empty()) )
      {
        coef_strs.push_back( coef_g[col] );
        if( (coef_h.size() > col) && (!coef_h[col].empty()) )
          coef_strs.push_back( coef_h[col] );
      }//if( we have coef g or h )
      
      bool success_parse_all_coef = true;
      size_t last_non_zero_coef = 0;
      vector<float> coef_values;
      for( size_t coef_index = 0; coef_index < coef_strs.size(); ++coef_index )
      {
        const string &str = coef_strs[coef_index];
        float coef;
        if( !SpecUtils::parse_float(str.c_str(), str.size(), coef) )
        {
          success_parse_all_coef = false;
          break;
        }
        
        coef_values.push_back( coef );
        
        if( coef != 0.0f )
          last_non_zero_coef = coef_index;
      }//for( loop over coefficients )
        
      if( !success_parse_all_coef )
      {
        warnings.push_back( "Failed to parse all eqn coefficients for detector '" + name + "' - skipping detector" );
        continue;
      }
      
      if( last_non_zero_coef < 1 )
      {
        warnings.push_back( "No non-zero eqn coefficients for detector '" + name + "' - skipping detector" );
        continue;
      }
      
      coef_values.resize( last_non_zero_coef + 1 );
      
      EffGeometryType geometry_type = EffGeometryType::FarField;
      if( is_fixed_geom )
      {
        if( SpecUtils::icontains(name, "Bq-cm2") )
        {
          geometry_type = EffGeometryType::FixedGeomActPerCm2;
        }else if( SpecUtils::icontains(name, "OnContact") || SpecUtils::icontains(comment, "OnContact") )
        {
          geometry_type = EffGeometryType::FixedGeomTotalAct;
        }else if( SpecUtils::icontains(name, "Bq-g") )
        {
          geometry_type = EffGeometryType::FixedGeomActPerGram;
        }else
        {
          warnings.push_back( "Couldnt determine fixed-geometry type for detector '" + name + "' - skipping detector" );
          continue;
        }
      }//if( is_fixed_geom )
      
      try
      {
        auto det = make_shared<DetectorPeakResponse>();
        det->fromExpOfLogPowerSeriesAbsEff( coef_values, {}, src_to_det_face, 2*det_rad,
                                           energy_units, low_energy, up_energy, geometry_type );
        if( !det->isValid() )
          throw runtime_error( "Not valid" );
        
        det->setName( name );
        string desc = comment + ((!comment.empty() && !cal_geom.empty()) ? ", " : "") + cal_geom;
        det->setDescription( desc );
        det->setDrfSource( DetectorPeakResponse::DrfSource::UserImportedIntrisicEfficiencyDrf );
        
        drfs.push_back( det );
      }catch( std::exception &e )
      {
        warnings.push_back( "Detector '" + name + "' was invalid ('" + string(e.what())
                           + "') - skipping detector" );
      }
    }//for( size_t col = 1; col < num_col; ++col )
    
    if( drfs.empty() )
      throw runtime_error( "No valid detector efficiencies found." );
  }catch( std::exception &e )
  {
    cerr << "Failed to parse GammaQuant Det Effs: " << e.what() << endl;
    input.seekg( start_pos, ios::beg );
    input.clear( ios::failbit );
    throw;
  }
}//void DetectorPeakResponse::parseGammaQuantRelEffDrfCsv(...)


std::shared_ptr<DetectorPeakResponse> DetectorPeakResponse::parseFromAppUrl( const std::string &url_query )
{
  auto drf = make_shared<DetectorPeakResponse>();
  drf->fromAppUrl( url_query );
  assert( drf->isValid() );
  
  return drf;
}//shared_ptr<DetectorPeakResponse> parseFromAppUrl( const std::string &url_query )


std::string DetectorPeakResponse::toAppUrl() const
{
  // QR codes can hold:
  //  Numeric only: Max. 7,089 characters (0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  //  Alphanumeric: Max. 4,296 characters (0–9, A–Z [upper-case only], space, $, %, *, +, -, ., /, :)
  //  Binary/byte:  Max. 2,953 characters (8-bit bytes) (23624 bits)
  //
  // So we will try to fit the DRF inside of this limitation, by first not changing anything, and
  // strlen("interspec://drf/specify?") == 24, so we'll round up to 30.
  const size_t max_binary_num = 2953 - 30;
  //const size_t max_alpha_num = 4296 - 30;
  
  // Currently we are just trying to fit everything in the binary size
  
  //*********************************************************************
  // A note for the future:
  //  We are actually pretty close to being an QR Alphanumeric code.
  //  The '&' and '=' would have to be replaced by say like '/' and ':'
  //  and also the '?' that indicates the query string would have to
  //  become say '%3F', but we're mostly there.
  //  Also we would need to change calls to url_encode to be like:
  //    url_encode(...,"",true);
  //*********************************************************************
  
  
  if( !isValid() )
    throw runtime_error( "Invalid DRF." );
  
  map<string,string> parts;
  parts["VER"] = "1";
  
  if( !m_name.empty() )
    parts["NAME"] = url_encode( m_name, "", false );
  
  if( !m_description.empty() )
    parts["DESC"] = url_encode( m_description, "", false );
  
  if( (m_geomType == EffGeometryType::FarField) || (m_detectorDiameter > 0.0) )
    parts["DIAM"] = SpecUtils::printCompact( m_detectorDiameter, 5 );
  
  // We'll assume units are MeV, unless stated otherwise.
  const bool notKeV = (fabs(m_efficiencyEnergyUnits - PhysicalUnits::keV) > 0.001);
  if( notKeV )
    parts["EUNIT"] = SpecUtils::printCompact( m_efficiencyEnergyUnits, 4 );
  
  switch( m_efficiencyForm )
  {
    case EfficiencyFnctForm::kEnergyEfficiencyPairs:
    {
      assert( m_energyEfficiencies.size() );
      
      parts["EFFT"] = "P";
      vector<float> energies, effs;
      
      for( size_t i = 0; i < m_energyEfficiencies.size(); ++i )
      {
        const EnergyEfficiencyPair &v = m_energyEfficiencies[i];
        energies.push_back( v.energy );
        effs.push_back( v.efficiency );
      }
      
      parts["EFFX"] = to_url_flt_array( energies, 4 );
      parts["EFFY"] = to_url_flt_array( effs, 5 );
      break;
    }//case EfficiencyFnctForm::kEnergyEfficiencyPairs:
      
    case EfficiencyFnctForm::kFunctialEfficienyForm:
    {
      assert( m_efficiencyFormula.size() );
      
      parts["EFFT"] = "F";
      parts["EFFE"] = url_encode(m_efficiencyFormula, "", false);
      break;
    }
      
    case EfficiencyFnctForm::kExpOfLogPowerSeries:
    {
      assert( m_expOfLogPowerSeriesCoeffs.size() );
      
      parts["EFFT"] = "E";
      parts["EFFC"] = to_url_flt_array( m_expOfLogPowerSeriesCoeffs, 7 );
      if( !m_expOfLogPowerSeriesUncerts.empty() )
        parts["EFFU"] = to_url_flt_array( m_expOfLogPowerSeriesUncerts, 7 );
      break;
    }
      
    case EfficiencyFnctForm::kNumEfficiencyFnctForms:
      assert( 0 );
      throw std::logic_error("m_efficiencyForm");
      break;
  }//switch( m_efficiencyForm )
  
  
  switch( m_resolutionForm )
  {
    case ResolutionFnctForm::kGadrasResolutionFcn:
      parts["FWHMT"] = "GAD";
      break;
    case kSqrtEnergyPlusInverse:
      parts["FWHMT"] = "FRAM";
      break;
    case kConstantPlusSqrtEnergy:
      parts["FWHMT"] = "GENIE";
      break;
    case ResolutionFnctForm::kSqrtPolynomial:
      parts["FWHMT"] = "SQRTPOLY";
      break;
      
    case ResolutionFnctForm::kNumResolutionFnctForm:
      break;
  }//switch( m_resolutionForm )
  
  if( m_resolutionCoeffs.size() )
    parts["FWHMC"] = to_url_flt_array( m_resolutionCoeffs, 7 );
  
  if( m_resolutionUncerts.size() )
    parts["FWHMU"] = to_url_flt_array( m_resolutionUncerts, 7 );
  
  parts["ORIGIN"] = std::to_string( static_cast<int>(m_efficiencySource) );
  
  if( m_hash )
    parts["HASH"] = std::to_string( m_hash );
  
  if( m_parentHash )
    parts["HASHP"] = std::to_string( m_parentHash );
  
  if( (fabs(m_lowerEnergy - m_upperEnergy) > 1.0) && (m_upperEnergy > 0.0) )
  {
    parts["LOWE"] = SpecUtils::printCompact( m_lowerEnergy, 4 );
    parts["HIGHE"] = SpecUtils::printCompact( m_upperEnergy, 4 );
  }
  
  if( m_createdUtc )
    parts["CREATED"] = std::to_string( m_createdUtc );
  
  if( m_lastUsedUtc )
    parts["LASTUSED"] = std::to_string( m_lastUsedUtc );
  
  switch( m_geomType )
  {
    case EffGeometryType::FarField:
      // We wont note this, but if decide to in the future, use "FAR-FIELD"
      //parts["GEOM"] = "FAR-FIELD";
      break;
      
    case EffGeometryType::FixedGeomTotalAct:
      parts["GEOM"] = "FIXED-TOTAL";
      parts["FIXGEOM"] = "1"; //Left-over from when fixed geometry was just a boolean - can be removed probably
      break;
      
    case EffGeometryType::FixedGeomActPerCm2:
      parts["GEOM"] = "FIXED-PER-CM2";
      break;
      
    case EffGeometryType::FixedGeomActPerM2:
      parts["GEOM"] = "FIXED-PER-M2";
      break;
      
    case EffGeometryType::FixedGeomActPerGram:
      parts["GEOM"] = "FIXED-PER-GRAM";
      break;
  }//switch( m_geomType )
  
  auto current_url_len = [&parts]() -> size_t {
    size_t nchar = (2 * parts.size()) - (parts.size() ? 1 : 0);
    for( const auto &p : parts )
      nchar += p.first.size() + p.second.size();
    return nchar;
  };//current_url_len(...)

  
  auto combine_parts = [&parts]() -> string {
    string answer;
    for( const auto &p : parts )
      answer += (answer.size() ? "&" : "") + p.first + "=" + p.second;
    
    /*
     //If we wanted to check if it was QR ASCII
    for( const auto &p : parts )
    {
      //0–9, A–Z [upper-case only], space, $, %, *, +, -, ., /, :
      const string allowed = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ $%*+-./:";
      for( auto c : parts.first )
      {
        assert( allowed.find(c) != string::npos );
      }
      
      for( auto c : parts.second )
      {
        assert( allowed.find(c) != string::npos );
      }
    }//for( const auto &p : parts )
    */
    
    return answer;
  };//
  
  
  //  For the DRFs I have, as of 20220408, they all come in under 1000 bytes, so we wont really
  //  usually need to shorten things, but we'll add this code in anyway
  
  if( current_url_len() < max_binary_num )
    return combine_parts();
  
  // Remove part of URL, and return if now short engough
  auto remove_part = [&parts,&current_url_len,max_binary_num]( const string &part ) -> bool {
    if( parts.count(part) )
      parts.erase( part );
    
    return (current_url_len() < max_binary_num);
  };//remove_part
  
  if( remove_part("LASTUSED") )
    return combine_parts();
  
  if( remove_part("ORIGIN") )
    return combine_parts();
  
  if( remove_part("HASHP") )
    return combine_parts();
  
  if( remove_part("CREATED") )
    return combine_parts();
  
  remove_part("LOWE");
  if( remove_part("HIGHE") )
    return combine_parts();
  
  if( remove_part("HASH") )
    return combine_parts();
  
  // TODO: try to remove as little of the description as possible... but shorten the string
  //       before URL encoding so you dont have hanging encoding (e.g., "%2" would make things
  //       invalid)
  if( remove_part("DESC") )
    return combine_parts();

  if( remove_part("FWHMU") )
    return combine_parts();
  
  if( remove_part("EFFU") )
    return combine_parts();
  
  remove_part("FWHMT");
  if( remove_part("FWHMC") )
    return combine_parts();
  
  if( remove_part("NAME") )
    return combine_parts();
  
  throw runtime_error( "toAppUrl: Unable to shorten information enough" );
  
  return "";
}//std::string toAppUrl() const


void DetectorPeakResponse::fromAppUrl( std::string url_query )
{
  map<string,string> parts = AppUtils::query_str_key_values( url_query );
  
  if( !parts.count("VER") || (parts["VER"] != "1") )
    throw runtime_error( "fromAppUrl: missing or invalid 'VER'" );
  
  string name, desc, eqn;
  float detectorDiameter = 0.0, efficiencyEnergyUnits = 1.0;
  EfficiencyFnctForm eff_form;
  vector<EnergyEfficiencyPair> energyEfficiencies;
  std::function<float(float)> efficiencyFcn;
  vector<float> expOfLogPowerSeriesCoeffs, expOfLogPowerSeriesUncerts;
  ResolutionFnctForm resolutionForm = ResolutionFnctForm::kNumResolutionFnctForm;
  vector<float> resolutionCoeffs, resolutionUncerts;
  DrfSource drf_source = DrfSource::UnknownDrfSource;
  uint64_t hash = 0, parent_hash = 0;
  int64_t createdUtc = 0, lastUsedUtc = 0;
  double lowerEnergy = 0.0, upperEnergy = 0.0;
  //bool fixedGeometry = false;
  EffGeometryType geom_type = EffGeometryType::FarField;
  
  if( parts.count("NAME") )
    name = parts["NAME"];
  
  if( parts.count("DESC") )
    desc = parts["DESC"];

  SpecUtils::ireplace_all( name, "%20", " " );
  SpecUtils::trim( name );

  SpecUtils::ireplace_all( desc, "%20", " " );
  SpecUtils::trim( desc );

  if( parts.count("FIXGEOM") )
  {
    // This is vestigial code that can probably be deleted - it was only used between
    //  20230916 and 20231010 on development builds (maybe it also on a released Android version?)
    const string &v = parts["FIXGEOM"];
    if( (v != "0") && (v != "1")
       && !SpecUtils::iequals_ascii(v, "true")
       && !SpecUtils::iequals_ascii(v, "false") )
      throw runtime_error( "fromAppUrl: FIXGEOM not boolean." );
     
    if( (v == "1") || SpecUtils::iequals_ascii(v, "true") )
      geom_type = EffGeometryType::FixedGeomTotalAct;
  }//if( parts.count("FIXGEOM") )
  
  if( parts.count("GEOM") )
  {
    const string &v = parts["GEOM"];
    if( SpecUtils::iequals_ascii(v, "FAR-FIELD") )
      geom_type = EffGeometryType::FarField;
    else if( SpecUtils::iequals_ascii(v, "FIXED-TOTAL") )
      geom_type = EffGeometryType::FixedGeomTotalAct;
    else if( SpecUtils::iequals_ascii(v, "FIXED-PER-CM2") )
      geom_type = EffGeometryType::FixedGeomActPerCm2;
    else if( SpecUtils::iequals_ascii(v, "FIXED-PER-M2") )
      geom_type = EffGeometryType::FixedGeomActPerM2;
    else if( SpecUtils::iequals_ascii(v, "FIXED-PER-GRAM") )
      geom_type = EffGeometryType::FixedGeomActPerGram;
    else
      throw runtime_error( "fromAppUrl: invalid GEOM field: '" + v + "'" );
  }//if( parts.count("GEOM") )
  
  //if( !parts.count("DIAM") && !fixedGeometry )
  if( !parts.count("DIAM") && (geom_type == EffGeometryType::FarField) )
    throw runtime_error( "fromAppUrl: missing required DIAM component" );
  
  if( parts.count("DIAM") )
  {
    if( !(stringstream(parts["DIAM"]) >> detectorDiameter) || (detectorDiameter <= 0.0) )
      throw runtime_error( "fromAppUrl: invalid DIAM component" );
  }else
  {
    detectorDiameter = 0.0;
  }//if( parts.count("DIAM") )
  
  if( parts.count("EUNIT") )
  {
    if( !(stringstream(parts["EUNIT"]) >> efficiencyEnergyUnits) || (efficiencyEnergyUnits <= 0.0) )
      throw runtime_error( "fromAppUrl: invalid EUNIT component" );
  }
  
  if( !parts.count("EFFT") )
    throw runtime_error( "fromAppUrl: missing required EFFT component" );
  
  
  SpecUtils::to_upper_ascii( parts["EFFT"] );
  if( parts["EFFT"] == "P" )
  {
    eff_form = EfficiencyFnctForm::kEnergyEfficiencyPairs;
  
    if( !parts.count("EFFX") || !parts.count("EFFY") )
      throw runtime_error( "fromAppUrl: missing required EFFX or EFFY component for Eff Pair" );
      
    const vector<float> x = from_url_flt_array( parts["EFFX"] );
    const vector<float> y = from_url_flt_array( parts["EFFY"] );
    
    if( (x.size() < 2) || (x.size() != y.size()) )
      throw runtime_error( "fromAppUrl: missing required EFFX or EFFY component for Eff Pair" );
    
    for( size_t i = 0; i < x.size(); ++i )
    {
      EnergyEfficiencyPair p;
      p.energy = x[i];
      p.efficiency = y[i];
      energyEfficiencies.push_back( p );
    }
  }else if( parts["EFFT"] == "F" )
  {
    eff_form = EfficiencyFnctForm::kFunctialEfficienyForm;
    
    if( !parts.count("EFFE") )
      throw runtime_error( "fromAppUrl: missing required EFFE component for Eff Eqn" );
    
    eqn = parts["EFFE"];
    
    try
    {
      const bool isMeV = (efficiencyEnergyUnits > 10.0f);
      auto expression = std::make_shared<FormulaWrapper>( m_efficiencyFormula, isMeV );
      efficiencyFcn = boost::bind( &FormulaWrapper::efficiency, expression, boost::placeholders::_1  );
    }catch( std::exception &e )
    {
      throw runtime_error( "fromAppUrl: Invalid detector efficiency formula: " + string(e.what()) );
    }
  }else if( parts["EFFT"] == "E" )
  {
    eff_form = EfficiencyFnctForm::kExpOfLogPowerSeries;
    
    if( !parts.count("EFFC") )
      throw runtime_error( "fromAppUrl: missing EFFC for " );
      
    if( !parts.count("EFFC") )
      throw runtime_error( "fromAppUrl: missing required EFFC component for Eff Series" );
      
    expOfLogPowerSeriesCoeffs = from_url_flt_array( parts["EFFC"] );
    if( expOfLogPowerSeriesCoeffs.empty() )
      throw runtime_error( "fromAppUrl: invalid EFFC component for Eff Series" );
    
    if( parts.count("EFFU") )
    {
      expOfLogPowerSeriesUncerts = from_url_flt_array( parts["EFFU"] );
      if( !expOfLogPowerSeriesUncerts.empty() )
      {
        if( expOfLogPowerSeriesUncerts.size() != expOfLogPowerSeriesCoeffs.size() )
          throw runtime_error( "fromAppUrl: EFFC and EFFU are different lengths" );
      }
    }
  }else
  {
    throw runtime_error( "fromAppUrl: invalid EFFT value: " + parts["EFFT"]  );
  }
  
                              
  if( parts.count("FWHMT") )
  {
    if( parts["FWHMT"] == "GAD" )
    {
      resolutionForm = ResolutionFnctForm::kGadrasResolutionFcn;
    }else if( parts["FWHMT"] == "FRAM" )
    {
      resolutionForm = ResolutionFnctForm::kSqrtEnergyPlusInverse;
    }else if( parts["FWHMT"] == "GENIE" )
    {
      resolutionForm = ResolutionFnctForm::kConstantPlusSqrtEnergy;
    }else if( parts["FWHMT"] == "SQRTPOLY" )
    {
      resolutionForm = ResolutionFnctForm::kSqrtPolynomial;
    }else
    {
      throw runtime_error( "fromAppUrl: invalid FWHMT value '" + parts["FWHMT"] + "'" );
    }
    
    resolutionCoeffs = from_url_flt_array( parts["FWHMC"] );
            
    if( parts.count("FWHMU") )
      resolutionUncerts = from_url_flt_array( parts["FWHMU"] );
      
    if( resolutionUncerts.size() && (resolutionUncerts.size() != resolutionCoeffs.size()) )
       throw runtime_error( "fromAppUrl: FWHMC and FWHMU different lengths" );
      
    bool invalidNumCoef = false;
    switch( resolutionForm )
    {
      case ResolutionFnctForm::kGadrasResolutionFcn:
      case ResolutionFnctForm::kSqrtEnergyPlusInverse:
        invalidNumCoef = (resolutionCoeffs.size() != 3);
         
        break;
        
      case ResolutionFnctForm::kConstantPlusSqrtEnergy:
        invalidNumCoef = (resolutionCoeffs.size() != 2);
         break;
        
      case ResolutionFnctForm::kSqrtPolynomial:
        invalidNumCoef = (resolutionCoeffs.size() < 2);
        break;
       
      case ResolutionFnctForm::kNumResolutionFnctForm:
        break;
    }//switch( resolutionForm )
    
    if( invalidNumCoef )
      throw runtime_error( "fromAppUrl: invalid number resolution coefs" );
  }//if( parts.count("FWHMT") )
      
  if( parts.count("ORIGIN") )
  {
    int val;
    if( !(stringstream(parts["ORIGIN"]) >> val) )
      throw runtime_error( "fromAppUrl: ORIGIN must be an integer value." );
    
    drf_source = static_cast<DrfSource>( val );
    switch( drf_source )
    {
      case DrfSource::UnknownDrfSource:
      case DrfSource::DefaultGadrasDrf:
      case DrfSource::UserAddedGadrasDrf:
      case DrfSource::UserAddedRelativeEfficiencyDrf:
      case DrfSource::DefaultRelativeEfficiencyDrf:
      case DrfSource::UserImportedIntrisicEfficiencyDrf:
      case DrfSource::UserImportedGadrasDrf:
      case DrfSource::UserSpecifiedFormulaDrf:
      case DrfSource::UserCreatedDrf:
      case DrfSource::FromSpectrumFileDrf:
      case DrfSource::IsocsEcc:
        break;
      
      default:
        throw runtime_error( "fromAppUrl: invalid ORIGIN value." );
    }//switch( static_cast<DrfSource>( val ) )
  }//if( parts.count("ORIGIN") )
       
  if( parts.count("HASH") )
  {
     if( !(stringstream(parts["HASH"]) >> hash) )
       throw runtime_error( "fromAppUrl: HASH not an integer." );
  }//if( parts.count("HASH") )
       
  if( parts.count("HASHP") )
  {
    if( !(stringstream(parts["HASHP"]) >> parent_hash) )
      throw runtime_error( "fromAppUrl: HASHP not an integer." );
  }//if( parts.count("HASH") )
    
  if( parts.count("CREATED") )
  {
    if( !(stringstream(parts["CREATED"]) >> createdUtc) )
      throw runtime_error( "fromAppUrl: CREATED not an integer." );
  }//if( parts.count("CREATED") )
       
  if( parts.count("LASTUSED") )
  {
    if( !(stringstream(parts["LASTUSED"]) >> lastUsedUtc) )
      throw runtime_error( "fromAppUrl: LASTUSED not an integer." );
  }//if( parts.count("LASTUSED") )
  
  if( parts.count("LOWE") )
  {
    if( !(stringstream(parts["LOWE"]) >> lowerEnergy) )
      throw runtime_error( "fromAppUrl: LOWE not an float." );
  }
  
  if( parts.count("HIGHE") )
  {
    if( !(stringstream(parts["HIGHE"]) >> upperEnergy) )
      throw runtime_error( "fromAppUrl: HIGHE not an float." );
      
    // We wont actually check this, as its not required to be isValid()
    //if( upperEnergy < lowerEnergy )
    //  throw runtime_error( "fromAppUrl: HIGHE is lower than LOWE." );
  }
 
  // We should be good to go here
  m_name = name;
  m_description = desc;
       
  m_detectorDiameter = detectorDiameter;
  m_efficiencyEnergyUnits = efficiencyEnergyUnits;
       
  m_efficiencyForm = eff_form;
  m_efficiencyFormula = eqn;
  m_efficiencyFcn = efficiencyFcn;
  m_energyEfficiencies = energyEfficiencies;
  m_expOfLogPowerSeriesCoeffs = expOfLogPowerSeriesCoeffs;
  m_expOfLogPowerSeriesUncerts = expOfLogPowerSeriesUncerts;
       
  m_resolutionForm = resolutionForm;
  m_resolutionCoeffs = resolutionCoeffs;
  m_resolutionUncerts = resolutionUncerts;
       
  m_efficiencySource = drf_source;
       
  m_hash = hash;
  m_parentHash = parent_hash;
  m_createdUtc = createdUtc;
  m_lastUsedUtc = lastUsedUtc;
       
  m_lowerEnergy = lowerEnergy;
  m_upperEnergy = upperEnergy;
  
  m_geomType = geom_type;
  
  if( !isValid() )
    throw runtime_error( "fromAppUrl: DRF is invalid - even though it shouldnt be - logic error in this function." );
}//void fromAppUrl( std::string url_query )


tuple<shared_ptr<DetectorPeakResponse>,double,double>
  DetectorPeakResponse::parseEccFile( std::istream &input )
{
  /*
   SGI_template: SPHERE
   ISOCS_file_name: PointSource.gis
   Detector_name: UserDetectorName
   Collimator_name: no_collimator

   Convrgence_[%]: 1.0000
   Test_description:
   Comment:
   Date_Time: Tue_Sep_19_16:30:38_2023
   Source_area_cm2:  3.14159e-4
   Source_grams:  5.23599e-9
   keV_eff_%err_effw_%cnvrg(i)_%cnvrg(i-1)_pntsN:    45.00    1.07165e-4   15.0    5.61115e-13   -0.000001   -0.000003     1022
   keV_eff_%err_effw_%cnvrg(i)_%cnvrg(i-1)_pntsN:    60.00    1.64972e-4   10.0    8.63793e-13   -0.000000   -0.000003     1022
   keV_eff_%err_effw_%cnvrg(i)_%cnvrg(i-1)_pntsN:    80.00    1.97791e-4   10.0    1.03563e-12   -0.000001   -0.000002     1022
   ...
   */
  
  // We dont currently handle DRF uncertainties (other than m_expOfLogPowerSeriesUncerts, which
  //  we dont actually use anywhere anyway...), but in the future we hopefully will, so we'll
  //  parse them here to.
  string line;
  vector<pair<float,float>> energy_error;
  vector<EnergyEfficiencyPair> energy_efficiencies;
  
  double source_area = 0.0, source_mass = 0.0;
  string det_name, ISOCS_fname, coll_name, test_desc, comment, date_time, src_area, src_grams;
  
  while( SpecUtils::safe_get_line( input, line, 8192 ) )
  {
    SpecUtils::trim( line );
    
    string label;
    const size_t label_pos = line.find(':');
    if( label_pos != string::npos )
    {
      label = line.substr(0, label_pos);
      line = line.substr( label_pos + 1 );
      SpecUtils::trim( line );
    }
    
    if( SpecUtils::istarts_with(label, "Detector_name") )
    {
      det_name = line;
    }else if( SpecUtils::istarts_with(label, "ISOCS_file_name") )
    {
      ISOCS_fname = line;
    }else if( SpecUtils::istarts_with(label, "Collimator_name") )
    {
      if( line != "no_collimator" )
        coll_name = line;
    }else if( SpecUtils::istarts_with(label, "Test_description") )
    {
      test_desc = line;
    }else if( SpecUtils::istarts_with(label, "Comment") )
    {
      comment = line;
    }else if( SpecUtils::istarts_with(label, "Date_Time") )
    {
      date_time = line;
    }else if( SpecUtils::istarts_with(label, "Source_area_cm2") )
    {
      src_area = line + " cm2";  //It looks like "Source_area_cm2" is always used, so we can assume cm2 always.
    }else if( SpecUtils::istarts_with(label, "Source_grams") )
    {
      src_grams = line + " g";
    }else if( SpecUtils::istarts_with(label, "keV_eff_%err") )
    {
      // line == "45.00    1.07165e-4   15.0    5.61115e-13   -0.000001   -0.000003     1022"
      vector<float> values;
      if( !SpecUtils::split_to_floats( line, values ) )
        cerr << "Warning: didnt parse ECC line entirely to floats: '" << line << "'" << endl;
      if( values.size() != 7 )
        cerr << "Warning: didnt parse ECC line to 7 floats: '" << line << "' (got " << values.size()
             << " floats)" << endl;
      
      if( values.size() > 3 )
      {
        EnergyEfficiencyPair ene_eff;
        ene_eff.energy = values[0];
        ene_eff.efficiency = values[1];
        energy_efficiencies.push_back( ene_eff );
        
        if( ene_eff.energy <= 10.0 ) //ISOCS only goes down to 45 keV, but we'll be generous
          throw runtime_error( "parseEccFile: energy <10 keV (" + to_string(ene_eff.energy) + ")" );
        if( ene_eff.energy > 8000.0 ) //ISOCS only goes up to 7 MeV, but we'll be generous
          throw runtime_error( "parseEccFile: energy >8 MeV (" + to_string(ene_eff.energy) + ")" );
        
        if( (ene_eff.efficiency < 0.0) || IsNan(ene_eff.efficiency) || IsInf(ene_eff.efficiency) )
          throw runtime_error( "parseEccFile: efficiency <0 ("
                              + to_string(ene_eff.efficiency) + " at " + to_string(ene_eff.energy)
                              + " keV)" );
        
        energy_error.emplace_back( values[0], 0.001f * values[2] );
      }//if( values.size() > 3 )
    }//if( label is some value ) / else if( ... )
  }//while( SpecUtils::safe_get_line( input, line, 8192 ) )
  
  if( energy_efficiencies.size() < 4 )
    throw runtime_error( "parseEccFile: not enough energy efficiency pairs found." );
  
  // Make sure efficiency points are sorted - even though they already should be.
  std::sort( begin(energy_efficiencies), end(energy_efficiencies),
    []( const EnergyEfficiencyPair &lhs, const EnergyEfficiencyPair &rhs) -> bool {
    return lhs.energy < rhs.energy;
  } );
  
  shared_ptr<DetectorPeakResponse> answer = make_shared<DetectorPeakResponse>();
  answer->m_name = det_name;
  if( !ISOCS_fname.empty() )
    answer->m_name += (answer->m_name.empty() ? "" : " - ") + ISOCS_fname;
  if( !coll_name.empty() )
    answer->m_name += (answer->m_name.empty() ? "" : " - ") + ISOCS_fname;
  
  answer->m_description = comment;
  if( !test_desc.empty() )
    answer->m_description += (answer->m_description.empty() ? "" : ". ") + test_desc;
  //if( !date_time.empty() )
  //  answer->m_description += (answer->m_description.empty() ? "" : " - ") + date_time;
  if( !src_area.empty() )
  {
    try
    {
      //Ex. "300 cm2"
      // Instead we could use regex from PhysicalUnits to extract the number.
      const size_t pos = src_area.find("cm2");
      if( pos == string::npos )
        throw runtime_error( "no 'cm2' in src area string." );
      const string number = SpecUtils::trim_copy( src_area.substr(0,pos) );
      
      if( !SpecUtils::parse_double( number.c_str(), number.length(), source_area ) )
        throw runtime_error( "Failed to parse '" + number + "' as a number." );
      
      source_area *= PhysicalUnits::cm2;
    }catch( std::exception &e )
    {
      cerr << "Failed to interpret '" << src_area << "' as a surface area: " << e.what() << "." << endl;
    }
    
    answer->m_description += (answer->m_description.empty() ? "" : ". Src area: ") + src_area;
  }
  if( !src_grams.empty() )
  {
    try
    {
      source_mass = PhysicalUnits::stringToMass( src_grams );
    }catch( std::exception &e )
    {
      cerr << "Failed to interpret '" << src_grams << "' as a mass: " << e.what() << "." << endl;
    }
    answer->m_description += (answer->m_description.empty() ? "" : ". Src mass: ") + src_grams;
  }
  
  answer->m_detectorDiameter = 0.0f;
  answer->m_efficiencyEnergyUnits = static_cast<float>(PhysicalUnits::keV);
  answer->m_resolutionForm = ResolutionFnctForm::kNumResolutionFnctForm;
  answer->m_efficiencySource = DrfSource::IsocsEcc;
  answer->m_efficiencyForm = EfficiencyFnctForm::kEnergyEfficiencyPairs;
  answer->m_energyEfficiencies = energy_efficiencies;
  answer->m_flags = 0;
  answer->m_lowerEnergy = energy_efficiencies.front().energy;
  answer->m_upperEnergy = energy_efficiencies.back().energy;
  answer->m_createdUtc = std::time(nullptr);
  answer->m_lastUsedUtc = answer->m_createdUtc;
  answer->m_geomType = EffGeometryType::FixedGeomTotalAct;
  answer->m_parentHash = 0;
  answer->computeHash();
  
  return {answer, source_area, source_mass};
}//shared_ptr<DetectorPeakResponse> parseEccFile( std::istream &input )


shared_ptr<DetectorPeakResponse> DetectorPeakResponse::convertFixedGeometryToFarField(
                                                        const double diameter,
                                                        const double distance,
                                                        const bool correct_for_air_atten ) const
{
  if( !isValid() )
    throw runtime_error( "DetectorPeakResponse::convertFixedGeometryToFarField: Invalid input DRF" );
  
  if( m_geomType != EffGeometryType::FixedGeomTotalAct )
    throw runtime_error( "DetectorPeakResponse::convertFixedGeometryToFarField:"
                        " Input DRF not fixed-geometry" );
  
  if( (distance < 0.0) || IsInf(distance) || IsNan(distance) )
    throw runtime_error( "DetectorPeakResponse::convertFixedGeometryToFarField: Invalid distance" );
  
  if( (diameter <= 0.0) || IsInf(diameter) || IsNan(diameter) )
    throw runtime_error( "DetectorPeakResponse::convertFixedGeometryToFarField: Invalid diameter" );
  
  if( correct_for_air_atten &&
     ( (m_efficiencyForm != EfficiencyFnctForm::kEnergyEfficiencyPairs)
      && (m_efficiencyForm != EfficiencyFnctForm::kFunctialEfficienyForm)) )
    throw runtime_error( "DetectorPeakResponse::convertFixedGeometryToFarField: air attenuation"
                        " correction only allowed if DRF defined using energy-efficiency pairs"
                        " or a functional efficiency form" );
  
  
  shared_ptr<DetectorPeakResponse> answer = make_shared<DetectorPeakResponse>(*this);
  answer->m_geomType = EffGeometryType::FarField;
  answer->m_detectorDiameter = diameter;
  
  const double energy_units = answer->m_efficiencyEnergyUnits;
  const float distancef = static_cast<float>( distance );
  
  
  const double frac_angle = fractionalSolidAngle( diameter, distance );
  if( (frac_angle <= 0.0) || IsInf(frac_angle) || IsNan(frac_angle) )
    throw runtime_error( "DetectorPeakResponse::convertFixedGeometryToFarField: invalid"
                         " fractional solid angle" );
  
  switch( answer->m_efficiencyForm )
  {
    case kEnergyEfficiencyPairs:
    {
      for( EnergyEfficiencyPair &ene_eff : answer->m_energyEfficiencies )
      {
        ene_eff.efficiency /= frac_angle;
        
        if( correct_for_air_atten && (ene_eff.efficiency > 0.0) )
        {
          const float energy = static_cast<float>( energy_units * ene_eff.energy );
          const double mu = GammaInteractionCalc::transmission_coefficient_air( energy, distancef );
          const double transmission_frac = exp( -mu );
          if( (transmission_frac <= 0.0) || IsInf(transmission_frac) || IsNan(transmission_frac) )
            throw runtime_error( "DetectorPeakResponse::convertFixedGeometryToFarField: air"
                                " attenuation correction at "
                                + std::to_string(floor(100*energy + 0.5)/100.0)
                                + " keV was not possible." );
          
          ene_eff.efficiency /= transmission_frac;
        }//if( correct_for_air_atten )
      }//for( loop over answer->m_energyEfficiencies )
      
      break;
    }//case kEnergyEfficiencyPairs:
    
      
    case kFunctialEfficienyForm:
    {
      if( !answer->m_efficiencyFcn )
        throw runtime_error( "DetectorPeakResponse::convertFixedGeometryToFarField: "
                            "no function defined, which is unexpected." );
      function<float(float)> old_fnctn = answer->m_efficiencyFcn;
      
      answer->m_efficiencyFcn
      = [old_fnctn, frac_angle, correct_for_air_atten, distancef]( float energy ) -> float {
        const float fixed_geom_eff = old_fnctn( energy );
        double eff = fixed_geom_eff / frac_angle;
        if( correct_for_air_atten && (eff > 0.0) ){
          const double mu = GammaInteractionCalc::transmission_coefficient_air( energy, distancef );
          const double transmission_frac = exp( -mu );
          // We should probably check `transmission_frac` is not zero - but I guess we'll have to
          //  look for this NaN later...
          eff /= transmission_frac;
        }
        
        return static_cast<float>( eff );
      };//answer->m_efficiencyFcn lambda defintion
      
      break;
    }//case kFunctialEfficienyForm:
      
    case kExpOfLogPowerSeries:
    {
      assert( !correct_for_air_atten );
      if( m_expOfLogPowerSeriesCoeffs.empty() )
        throw runtime_error( "DetectorPeakResponse::convertFixedGeometryToFarField: "
                            "no coefficients defined, which is unexpected." );
      
      answer->m_expOfLogPowerSeriesCoeffs[0] += static_cast<float>( log(1.0 / frac_angle) );
      
      // TODO: if correct_for_air_atten is true, we could try to re-fit for the equation, or something...
      
      break;
    }//case kExpOfLogPowerSeries:
      
    case kNumEfficiencyFnctForms:
      assert( 0 );
      throw runtime_error( "DetectorPeakResponse::convertFixedGeometryToFarField: invalid function form" );
      break;
  }//switch( m_efficiencyForm )
  
  answer->computeHash();

  return answer;
}//shared_ptr<DetectorPeakResponse> convertFixedGeometryToFarField(...)


shared_ptr<DetectorPeakResponse> DetectorPeakResponse::convertFixedGeometryType( const double quantity,
                                                            const EffGeometryType to_type ) const
{
  if( !isValid() )
    throw runtime_error( "DetectorPeakResponse::convertFixedGeometryType: Invalid input DRF" );
  
  if( m_geomType == EffGeometryType::FarField )
    throw runtime_error( "DetectorPeakResponse::convertFixedGeometryType:"
                        " Input DRF not fixed-geometry" );
  
  if( (quantity <= 0.0) || IsInf(quantity) || IsNan(quantity) )
    throw runtime_error( "DetectorPeakResponse::convertFixedGeometryType: Invalid surface_area" );

  shared_ptr<DetectorPeakResponse> answer = make_shared<DetectorPeakResponse>(*this);
  if( answer->m_geomType == to_type )
    return answer;

  answer->m_geomType = to_type;
  
  double correction = quantity;
  switch( to_type )
  {
    case EffGeometryType::FarField:
      assert( 0 );
      break;
      
    case EffGeometryType::FixedGeomTotalAct:
      switch( answer->m_geomType )
      {
        case EffGeometryType::FarField:
        case EffGeometryType::FixedGeomTotalAct:
          assert( 0 );
          break;
          
        case EffGeometryType::FixedGeomActPerCm2:
          correction /= PhysicalUnits::cm2;
          break;
          
        case EffGeometryType::FixedGeomActPerM2:
          correction /= PhysicalUnits::m2;
          break;
          
        case EffGeometryType::FixedGeomActPerGram:
          correction /= PhysicalUnits::gram;
          break;
      }//switch( answer->m_geomType )
      
      correction = 1.0 / correction;
      break;
      
    case EffGeometryType::FixedGeomActPerCm2:
      correction /= PhysicalUnits::cm2;
      break;
      
    case EffGeometryType::FixedGeomActPerM2:
      correction /= PhysicalUnits::m2;
      break;
      
    case EffGeometryType::FixedGeomActPerGram:
      correction /= PhysicalUnits::gram;
      break;
  }//switch( to_type )
  
  if( (correction <= 0.0) || IsNan(correction) || IsInf(correction) )
    throw runtime_error( "DetectorPeakResponse::convertFixedGeometryType: "
                         "The correction became invalid (" + to_string(correction) + ")." );
  
  switch( answer->m_efficiencyForm )
  {
    case kEnergyEfficiencyPairs:
    {
      for( EnergyEfficiencyPair &ene_eff : answer->m_energyEfficiencies )
      {
        ene_eff.efficiency *= correction;
      }//for( loop over answer->m_energyEfficiencies )
      
      break;
    }//case kEnergyEfficiencyPairs:
    
      
    case kFunctialEfficienyForm:
    {
      if( !answer->m_efficiencyFcn )
        throw runtime_error( "DetectorPeakResponse::convertFixedGeometryType: "
                            "no function defined, which is unexpected." );
      function<float(float)> old_fnctn = answer->m_efficiencyFcn;
      
      answer->m_efficiencyFcn
      = [old_fnctn, correction]( float energy ) -> float {
        return static_cast<float>( correction * old_fnctn(energy) );
      };//answer->m_efficiencyFcn lambda defintion
      
      break;
    }//case kFunctialEfficienyForm:
      
    case kExpOfLogPowerSeries:
    {
      if( m_expOfLogPowerSeriesCoeffs.empty() )
        throw runtime_error( "DetectorPeakResponse::convertFixedGeometryTypeToFarField: "
                            "no coefficients defined, which is unexpected." );
      
      answer->m_expOfLogPowerSeriesCoeffs[0] += static_cast<float>( log(correction) );
      
      break;
    }//case kExpOfLogPowerSeries:
      
    case kNumEfficiencyFnctForms:
      assert( 0 );
      throw runtime_error( "DetectorPeakResponse::convertFixedGeometryTypeToFarField: invalid function form" );
      break;
  }//switch( m_efficiencyForm )
  
  answer->computeHash();

  return answer;
}//std::shared_ptr<DetectorPeakResponse> convertFixedGeometryType( const double quanitity, const EffGeometryType to_type ) const;


void DetectorPeakResponse::setFwhmCoefficients( const std::vector<float> &coefs,
                         const ResolutionFnctForm form )
{
  switch( form )
  {
    case ResolutionFnctForm::kSqrtPolynomial:
      if( coefs.empty() )
        throw runtime_error( "setFwhmCoefficients: Sqrt polynomial equation must have at least one coefficient." );
      break;
      
    case ResolutionFnctForm::kGadrasResolutionFcn:
      if( coefs.size() != 3 )
        throw runtime_error( "setFwhmCoefficients: GADRAS equation must have three coefficients." );
      break;
      
    case ResolutionFnctForm::kSqrtEnergyPlusInverse:
      if( coefs.size() != 3 )
        throw runtime_error( "setFwhmCoefficients: sqrt(A0+A1*E+A2/E) equation must have three coefficients." );
      break;
    
    case ResolutionFnctForm::kConstantPlusSqrtEnergy:
      if( coefs.size() != 2 )
        throw runtime_error( "setFwhmCoefficients: A0 + A1*sqrt(E) equation must have two coefficients." );
      break;
      
    case ResolutionFnctForm::kNumResolutionFnctForm:
      if( !coefs.empty() )
        throw runtime_error( "setFwhmCoefficients: NumResolutionFnctForm must not have any coefficients." );
      break;
  }//switch( form )
  
  m_resolutionForm = form;
  m_resolutionCoeffs = coefs;
  
  computeHash();
}//void setFwhmCoefficients(...)


void DetectorPeakResponse::toXml( ::rapidxml::xml_node<char> *parent,
                                  ::rapidxml::xml_document<char> *doc ) const
{
  using namespace rapidxml;
  
  char buffer[128];
  
  xml_node<char> *base_node = doc->allocate_node( node_element, "DetectorPeakResponse" );
  parent->append_node( base_node );
  
  // We will write XML version 0, if m_geomType is not far-field (this was only change between 0 and 1)
  static_assert( sm_xmlSerializationVersion == 2, "Update DetectorPeakResponse sm_xmlSerializationVersion");
  
  if( m_resolutionForm == ResolutionFnctForm::kConstantPlusSqrtEnergy )
  {
    snprintf( buffer, sizeof(buffer), "%i", sm_xmlSerializationVersion );
  }else
  {
    snprintf( buffer, sizeof(buffer), "%i", ((m_geomType != EffGeometryType::FarField) ? 1 : 0) );
  }
  
  const char *value = doc->allocate_string( buffer );
  xml_attribute<char> *attr = doc->allocate_attribute( "version", value );
  base_node->append_attribute( attr );
  
  const char *val = doc->allocate_string( m_name.c_str() );
  xml_node<char> *node = doc->allocate_node( node_element, "Name", val );
  base_node->append_node( node );
  
  val = doc->allocate_string( m_description.c_str() );
  node = doc->allocate_node( node_element, "Description", val );
  base_node->append_node( node );
  
  if( (m_geomType == EffGeometryType::FarField) || (m_detectorDiameter > 0.0) )
  {
    snprintf( buffer, sizeof(buffer), "%1.8E", m_detectorDiameter );
    val = doc->allocate_string( buffer );
    node = doc->allocate_node( node_element, "DetectorDiameter", val );
    base_node->append_node( node );
  }//if( (m_geomType == EffGeometryType::FarField) || (m_detectorDiameter > 0.0) )
  
  val = "UnknownDrfSource";
  switch( m_efficiencySource )
  {
    case UnknownDrfSource:                                                             break;
    case DefaultGadrasDrf:                  val = "DefaultGadrasDrf";                  break;
    case UserAddedGadrasDrf:                val = "UserAddedGadrasDrf";                break;
    case UserAddedRelativeEfficiencyDrf:    val = "UserAddedRelativeEfficiencyDrf";    break;
    case DefaultRelativeEfficiencyDrf:      val = "DefaultRelativeEfficiencyDrf";      break;
    case UserImportedIntrisicEfficiencyDrf: val = "UserImportedIntrisicEfficiencyDrf"; break;
    case UserImportedGadrasDrf:             val = "UserImportedGadrasDrf";             break;
    case UserSpecifiedFormulaDrf:           val = "UserSpecifiedFormulaDrf";           break;
    case UserCreatedDrf:                    val = "UserCreatedDrf";                    break;
    case FromSpectrumFileDrf:               val = "FromSpectrumFileDrf";               break;
    case DrfSource::IsocsEcc:               val = "ISOCS";                             break;
  }//switch( m_efficiencySource )
  
  node = doc->allocate_node( node_element, "EfficiencySource", val );
  base_node->append_node( node );
  
  snprintf( buffer, sizeof(buffer), "%g", m_efficiencyEnergyUnits );
  val = doc->allocate_string( buffer );
  node = doc->allocate_node( node_element, "EfficiencyEnergyUnits", val );
  base_node->append_node( node );
  
  switch( m_resolutionForm )
  {
    case kGadrasResolutionFcn:    val = "GadrasResolutionFcn";   break;
    case kSqrtEnergyPlusInverse:  val = "SqrtEnergyPlusInverse"; break;
    case kConstantPlusSqrtEnergy: val = "ConstantPlusSqrtEnergy"; break;
    case kSqrtPolynomial:         val = "SqrtPolynomial";        break;
    case kNumResolutionFnctForm:  val = "Undefined";             break;
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
  
  
  if( m_expOfLogPowerSeriesUncerts.size() )
  {
    stringstream valstrm;
    for( size_t i = 0; i < m_expOfLogPowerSeriesUncerts.size(); ++i )
      valstrm << (i?" ":"") << m_expOfLogPowerSeriesUncerts[i];
    val = doc->allocate_string( valstrm.str().c_str() );
    node = doc->allocate_node( node_element, "ExpOfLogPowerSeriesUncerts", val );
    base_node->append_node( node );
  }//if( m_expOfLogPowerSeriesUncerts.size() )
  
  
  stringstream hashstrm, parenthashstrm, flagsstrm;
  hashstrm << m_hash;
  parenthashstrm << m_parentHash;
  flagsstrm << m_flags;
  
  val = doc->allocate_string( hashstrm.str().c_str() );
  node = doc->allocate_node( node_element, "Hash", val );
  base_node->append_node( node );
  
  val = doc->allocate_string( parenthashstrm.str().c_str() );
  node = doc->allocate_node( node_element, "ParentHash", val );
  base_node->append_node( node );
  
  val = doc->allocate_string( flagsstrm.str().c_str() );
  node = doc->allocate_node( node_element, "Flags", val );
  base_node->append_node( node );
  
  if( m_lowerEnergy != 0.0 || m_upperEnergy != 0.0 )
  {
    snprintf( buffer, sizeof(buffer), "%1.8E", m_lowerEnergy );
    val = doc->allocate_string( buffer );
    node = doc->allocate_node( node_element, "LowerEnergy", val );
    base_node->append_node( node );
    
    snprintf( buffer, sizeof(buffer), "%1.8E", m_upperEnergy );
    val = doc->allocate_string( buffer );
    node = doc->allocate_node( node_element, "UpperEnergy", val );
    base_node->append_node( node );
  }
  
  
  if( m_geomType == EffGeometryType::FixedGeomTotalAct )
  {
    // Added 20230916, e.g., for InterSpec v1.0.12
    //  - but then made irrelevant before v1.0.12 released, on 20231011
    // This section of code can be removed
    node = doc->allocate_node( node_element, "FixedGeometry", "1" );
    base_node->append_node( node );
  }//if( m_geomType == EffGeometryType::FixedGeomTotalAct )
  
  
  // Added 20231011 to describe geometry type
  const char *geom_type_str = "FAR-FIELD";
  switch( m_geomType )
  {
    case EffGeometryType::FarField:            break;
    case EffGeometryType::FixedGeomTotalAct:   geom_type_str = "FIXED-TOTAL";    break;
    case EffGeometryType::FixedGeomActPerCm2:  geom_type_str = "FIXED-PER-CM2";  break;
    case EffGeometryType::FixedGeomActPerM2:   geom_type_str = "FIXED-PER-M2";   break;
    case EffGeometryType::FixedGeomActPerGram: geom_type_str = "FIXED-PER-GRAM"; break;
  }//switch( m_geomType )
  
  node = doc->allocate_node( node_element, "Geometry", geom_type_str );
  base_node->append_node( node );
  
  
  if( m_createdUtc )
  {
    stringstream strm;
    strm << m_createdUtc;
    val = doc->allocate_string( strm.str().c_str() );
    node = doc->allocate_node( node_element, "CreationTimeUtc", val );
    base_node->append_node( node );
  }
  
  if( m_lastUsedUtc )
  {
    stringstream strm;
    strm << m_lastUsedUtc;
    val = doc->allocate_string( strm.str().c_str() );
    node = doc->allocate_node( node_element, "LastUsedTimeUtc", val );
    base_node->append_node( node );
  }
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
  
  if( version > sm_xmlSerializationVersion )
    throw runtime_error( "Invalid DetectorPeakResponse version" );

  const xml_node<char> *node = parent->first_node( "Name", 4 );
  if( !node || !node->value() )
    throw runtime_error( "DetectorPeakResponse missing Name node" );
  m_name = node->value();
  
  
  node = parent->first_node( "Description", 11 );
  if( !node || !node->value() )
    throw runtime_error( "DetectorPeakResponse missing Description node" );
  m_description = node->value();
  
  
  m_geomType = EffGeometryType::FarField;
  
  node = parent->first_node( "FixedGeometry", 13 );
  if( node )
  {
    // Added 20230916, e.g., for InterSpec v1.0.12
    if( compare(node->value(), node->value_size(), "1", 1, false)
       || compare(node->value(), node->value_size(), "true", 4, false) )
    {
      m_geomType = EffGeometryType::FixedGeomTotalAct;
    }else if( !compare(node->value(), node->value_size(), "0", 1, false)
            && !compare(node->value(), node->value_size(), "false", 5, false) )
    {
      throw runtime_error( "DetectorPeakResponse invalid FixedGeometry ('"
                          + SpecUtils::xml_value_str(node) + "')" );
    }
  }//if( not "FixedGeometry" node )
  
  
  // Geometry node added 20231011 to describe geometry type
  node = parent->first_node( "Geometry", 8 );
  if( !node && (version >= 1) && (m_geomType != EffGeometryType::FixedGeomTotalAct) )
    throw runtime_error( "DetectorPeakResponse no Geometry node" );
  
  if( node )
  {
    if( compare(node->value(), node->value_size(), "FAR-FIELD", 9, false) )
      m_geomType = EffGeometryType::FarField;
    else if( compare(node->value(), node->value_size(), "FIXED-TOTAL", 11, false) )
      m_geomType = EffGeometryType::FixedGeomTotalAct;
    else if( compare(node->value(), node->value_size(), "FIXED-PER-CM2", 13, false) )
      m_geomType = EffGeometryType::FixedGeomActPerCm2;
    else if( compare(node->value(), node->value_size(), "FIXED-PER-M2", 12, false) )
      m_geomType = EffGeometryType::FixedGeomActPerM2;
    else if( compare(node->value(), node->value_size(), "FIXED-PER-GRAM", 14, false) )
      m_geomType = EffGeometryType::FixedGeomActPerGram;
    else
      throw runtime_error( "DetectorPeakResponse has Geometry value: "
                          + string(node->value(), node->value() + node->value_size()) );
  }//if( node )
  
  node = parent->first_node( "DetectorDiameter", 16 );
  if( (!node || !node->value()) && (m_geomType == EffGeometryType::FarField) )
    throw runtime_error( "DetectorPeakResponse missing DetectorDiameter node" );
  if( node && node->value() )
    m_detectorDiameter = atof( node->value() );
  else
    m_detectorDiameter = 0.0;
  
  node = parent->first_node( "EfficiencySource", 16 );
  if( !node || !node->value() )
    throw runtime_error( "DetectorPeakResponse missing EfficiencySource node" );
  
  if( compare(node->value(),node->value_size(),"GadrasEfficiencyDefintion",25,false)
     || compare(node->value(),node->value_size(),"DefaultGadrasDrf",16,false) )
    m_efficiencySource = DrfSource::DefaultGadrasDrf;
  else if( compare(node->value(),node->value_size(),   "SimpleMassEfficiencyDefintion",29,false)
           || compare(node->value(),node->value_size(),"RelativeEfficiencyDefintion",27,false)
           || compare(node->value(),node->value_size(),"UserAddedRelativeEfficiencyDrf",30,false) )
    m_efficiencySource = DrfSource::UserAddedRelativeEfficiencyDrf;
  else if( compare(node->value(),node->value_size(),"DefaultRelativeEfficiencyDrf",28,false) )
    m_efficiencySource = DrfSource::DefaultRelativeEfficiencyDrf;
  else if( compare(node->value(),node->value_size(),"UserUploadedEfficiencyCsv",25,false)
          || compare(node->value(),node->value_size(),"UserImportedIntrisicEfficiencyDrf",33,false) )
    m_efficiencySource = DrfSource::UserImportedIntrisicEfficiencyDrf;
  else if( compare(node->value(),node->value_size(),"UserEfficiencyEquationSpecified",31,false)
          || compare(node->value(),node->value_size(),"UserSpecifiedFormulaDrf",23,false) )
    m_efficiencySource = DrfSource::UserSpecifiedFormulaDrf;
  else if( compare(node->value(),node->value_size(),"UnknownEfficiencySource",23,false)
          || compare(node->value(),node->value_size(),"UnknownDrfSource",16,false) )
    m_efficiencySource = DrfSource::UnknownDrfSource;
  else if( compare(node->value(),node->value_size(),"UserAddedGadrasDrf",18,false) )
    m_efficiencySource = DrfSource::UserAddedGadrasDrf;
  else if( compare(node->value(),node->value_size(),"UserImportedGadrasDrf",21,false) )
    m_efficiencySource = DrfSource::UserImportedGadrasDrf;
  else if( compare(node->value(),node->value_size(),"UserCreatedDrf",14,false) )
    m_efficiencySource = UserCreatedDrf;
  else if( compare(node->value(),node->value_size(),"FromSpectrumFileDrf",19,false) )
    m_efficiencySource = FromSpectrumFileDrf;
  else if( compare(node->value(),node->value_size(),"ISOCS",5,false) )
    m_efficiencySource = IsocsEcc;
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
  else if( compare(node->value(),node->value_size(),"SqrtEnergyPlusInverse",21,false) )
    m_resolutionForm = kSqrtEnergyPlusInverse;
  else if( compare(node->value(),node->value_size(),"ConstantPlusSqrtEnergy",22,false) )
    m_resolutionForm = kConstantPlusSqrtEnergy;
  else if( compare(node->value(),node->value_size(),"SqrtPolynomial",14,false) )
    m_resolutionForm = kSqrtPolynomial;
  else if( compare(node->value(),node->value_size(),"Undefined",9,false) )
    m_resolutionForm = kNumResolutionFnctForm;
  else
    m_resolutionForm = kNumResolutionFnctForm;
    //throw runtime_error( "DetectorPeakResponse: invalid ResolutionForm value" );
  
  
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
  if( m_resolutionForm != kNumResolutionFnctForm )
  {
    node = parent->first_node( "ResolutionCoefficients", 22 );
    if( node && node->value() )
      SpecUtils::split_to_floats( node->value(), node->value_size(), m_resolutionCoeffs );
  }//if( m_resolutionForm != kNumResolutionFnctForm )
  
  m_energyEfficiencies.clear();
  node = parent->first_node( "EnergyEfficiencies", 18 );
  if( node && node->value() )
  {
    vector<float> values;
    SpecUtils::split_to_floats( node->value(), node->value_size(), values );
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
    SpecUtils::trim( m_efficiencyFormula );
    
    if( m_efficiencyFormula.size() )
    {
      //if( m_efficiencySource != kUserEfficiencyEquationSpecified )
        //throw runtime_error( "An detector efficiency formula was specified but the the EfficiencySource had a different value" );
      const bool isMeV = (m_efficiencyEnergyUnits > 10.0f);
      
      try
      {
        auto expression = std::make_shared<FormulaWrapper>( m_efficiencyFormula, isMeV );
        m_efficiencyFcn = boost::bind( &FormulaWrapper::efficiency, expression,
                                      boost::placeholders::_1  );
      }catch( std::exception &e )
      {
        throw runtime_error( "Invalid detector efficiency formula in XML: " + string(e.what()) );
      }
    }//if( m_efficiencyFormula.size() )
  }//if( node && node->value() )
  
  m_expOfLogPowerSeriesCoeffs.clear();
  m_expOfLogPowerSeriesUncerts.clear();
  
  node = XML_FIRST_NODE(parent, "ExpOfLogPowerSeriesCoeffs");
  if( node && node->value() )
    SpecUtils::split_to_floats( node->value(), node->value_size(), m_expOfLogPowerSeriesCoeffs );
  
  
  node = XML_FIRST_NODE(parent, "ExpOfLogPowerSeriesUncerts");
  if( node && node->value() )
    SpecUtils::split_to_floats( node->value(), node->value_size(), m_expOfLogPowerSeriesUncerts );
  
  if( !m_expOfLogPowerSeriesUncerts.empty()
     && (m_expOfLogPowerSeriesUncerts.size() != m_expOfLogPowerSeriesCoeffs.size()) )
    throw runtime_error( "DetectorPeakResponse number eff coeffs doesnt match number of uncerts" );
  
  
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
  
  
  
  node = parent->first_node( "Flags", 5 );
  if( node && node->value_size() )
  {
    if( !(stringstream(node->value()) >> m_flags) )
      throw runtime_error( "DetectorPeakResponse invalid Flags" );
  }else
  {
    m_flags = 0;
  }
  
  
  node = parent->first_node( "LowerEnergy", 11 );
  if( node && node->value_size() )
    m_lowerEnergy = atof( node->value() );
  else
    m_lowerEnergy = 0.0f;
  
  node = parent->first_node( "UpperEnergy", 11 );
  if( node && node->value_size() )
    m_upperEnergy = atof( node->value() );
  else
    m_upperEnergy = 0.0f;
  
  node = parent->first_node( "CreationTimeUtc", 15 );
  if( node && node->value_size() )
  {
    if( !(stringstream(node->value()) >> m_createdUtc) )
      throw runtime_error( "DetectorPeakResponse invalid CreationTimeUtc" );
  }else
  {
    m_createdUtc = 0;
  }
  
  node = parent->first_node( "LastUsedTimeUtc", 15 );
  if( node && node->value_size() )
  {
    if( !(stringstream(node->value()) >> m_lastUsedUtc) )
      throw runtime_error( "DetectorPeakResponse invalid LastUsedTimeUtc" );
  }else
  {
    m_lastUsedUtc = 0;
  }
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
             lhs.m_efficiencyEnergyUnits, rhs.m_efficiencyEnergyUnits );
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
  
  if( lhs.m_expOfLogPowerSeriesUncerts.size() != rhs.m_expOfLogPowerSeriesUncerts.size() )
  {
    snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: size of"
             " exponential of log power series uncertainties of LHS (%i) doesnt"
             " match RHS (%i)",
             int(lhs.m_expOfLogPowerSeriesUncerts.size()),
             int(rhs.m_expOfLogPowerSeriesUncerts.size()) );
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
  
  for( size_t i = 0; i < lhs.m_expOfLogPowerSeriesUncerts.size(); ++i )
  {
    const float a = lhs.m_expOfLogPowerSeriesUncerts[i];
    const float b = rhs.m_expOfLogPowerSeriesUncerts[i];
    const float diff = fabs(a-b);
    if( fabs(a-b) > (1.0E-5 * std::max(fabs(a),fabs(b))) && (diff > 1.0E-8) )
    {
      snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: exponential of"
               " log power series uncertainty %i of LHS (%1.8E)"
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
  
  if( lhs.m_flags != rhs.m_flags )
  {
    snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: flags of LHS (%llud)"
             " doesnt match RHS (%llud)",
             (long long unsigned int)lhs.m_flags,
             (long long unsigned int)rhs.m_flags );
    throw runtime_error(buffer);
  }
  
  if( fabs(lhs.m_lowerEnergy - rhs.m_lowerEnergy) > 0.01 )
  {
    snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: lower energy"
             " of LHS (%1.8E) doesnt match RHS (%1.8E)",
             lhs.m_lowerEnergy, rhs.m_lowerEnergy );
    throw runtime_error(buffer);
  }
  
  if( fabs(lhs.m_upperEnergy - rhs.m_upperEnergy) > 0.01 )
  {
    snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: upper energy"
             " of LHS (%1.8E) doesnt match RHS (%1.8E)",
             lhs.m_upperEnergy, rhs.m_upperEnergy );
    throw runtime_error(buffer);
  }
  
  if( lhs.m_createdUtc != rhs.m_createdUtc )
  {
    snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: CreationTimeUtc"
             " of LHS (%lld) doesnt match RHS (%llud)",
             static_cast<long long int>(lhs.m_createdUtc),
             static_cast<long long int>(rhs.m_createdUtc) );
    throw runtime_error(buffer);
  }
  
  if( lhs.m_lastUsedUtc != rhs.m_lastUsedUtc )
  {
    snprintf( buffer, sizeof(buffer), "DetectorPeakResponse: LastUsedTimeUtc"
             " of LHS (%lld) doesnt match RHS (%llud)",
             static_cast<long long int>(lhs.m_lastUsedUtc),
             static_cast<long long int>(rhs.m_lastUsedUtc) );
    throw runtime_error(buffer);
  }
}//void equalEnough(...)
#endif //PERFORM_DEVELOPER_CHECKS



double DetectorPeakResponse::fractionalSolidAngle( const double detDiam,
                                                   const double D,
                                                   const double source_rad )
{
  /*
    See also page 119 in Knoll for approximations where source and detector
    diameters are somewhat comparable to distance of source from detector
    An unverified implementation of that is:
  */
  const double alpha = pow(source_rad/D,2.0);
  const double beta = pow(0.5*detDiam/D,2.0);
  const double F1 = ((5.0/16.0)*beta/pow(1.0+beta,3.5)) - ((35.0/64.0)*beta*beta/pow(1.0+beta,4.5));
  const double F2 = ((35.0/128.0)*beta/pow(1.0+beta,4.5)) - ((315.0/256.0)*beta*beta/pow(1.0+beta,11.0/12.0)) + ((1155.0/1024.0)*beta*beta*beta/pow(1.0+beta,6.5));

  const double omega = 0.5 * ( 1.0 - (1.0/sqrt(1.0+beta)) - ((3.0/8.0)*alpha*beta/pow(1.0+beta,2.5))
                       + alpha*alpha*F1 - alpha*alpha*alpha*F2 );

  return omega;
}//float DetectorPeakResponse::fractionalSolidAngle(...)


double DetectorPeakResponse::fractionalSolidAngle( const double detDiam, const double D ) noexcept
{
  // Using a metric that if you increase the distance by 1%, but get an identical results, you have
  //  hit the limit of numerical accuracy (not perfect, but whatever), using a 5cm Detector:
  //  - original float-based implementation: 12.65 m
  //  - upgrading to doubles: 194698 m
  //  - dividing the "D" through into the sqrt argument, reduced numerical accuracy for both float
  //    and double, as did moving the 0.5 through the parenthesis
  
  const double r = 0.5 * detDiam;
  return 0.5*(1.0 - (D/sqrt(D*D+r*r)));
  
  // The below should be a numerically more stable computation... not tested yet.
  //const double t_0 = D / sqrt((D * D) + (r * r));
  //const double t_1 = fma(D, (-0.5 / hypot(D, r)), 0.5);
  //if (t_0 <= 0.9674535662344166)
  //  return log(exp(t_1));
  //if (t_0 <= 1.0)
  //  return 0.25 * (pow(r, 2.0) / pow(D, 2.0));
  //return log1p(expm1(t_1));
}//double fractionalSolidAngle(...)


double DetectorPeakResponse::efficiency( const float energy, const double dist ) const
{
  const double fracSolidAngle = fractionalSolidAngle( m_detectorDiameter, dist );
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
  
  if( (n < 6) || (i < 2) || ((n-i) <= 3)  ) //We'll just use linear interpolation here
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
  double exparg = 0.0;
  const double x = log( static_cast<double>(energy) );
  for( size_t i = 0; i < coefs.size(); ++i )
    exparg += coefs[i] * pow(x,static_cast<double>(i));
  const double answer = exp( exparg );
  return (answer >= 0.0) ? static_cast<float>(answer) : 0.0f;
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



std::function<float( float )> DetectorPeakResponse::intrinsicEfficiencyFcn() const
{
  const double energy_units = m_efficiencyEnergyUnits;
  const double lowerEnergy = m_lowerEnergy;
  const double upperEnergy = m_upperEnergy;

  switch( m_efficiencyForm )
  {
    case kEnergyEfficiencyPairs:
    {
      if( m_energyEfficiencies.size() < 2 )
        return nullptr;

      const vector<EnergyEfficiencyPair> energy_effs = m_energyEfficiencies;

      return [energy_units, energy_effs]( float energy ) -> float {
        return akimaInterpolate( energy / energy_units, energy_effs );
      };
    }//case kEnergyEfficiencyPairs:

    case kFunctialEfficienyForm:
    {
      const std::function<float( float )> eff_fcnt = m_efficiencyFcn;
      if( !eff_fcnt )
        return nullptr;

      return [energy_units, eff_fcnt, lowerEnergy, upperEnergy]( float energy ) -> float {
        if( (lowerEnergy > 10.0) && (energy < lowerEnergy) )
          return eff_fcnt( lowerEnergy / energy_units );
        if( (upperEnergy > 10.0) && (energy > upperEnergy) )
          return eff_fcnt( upperEnergy / energy_units );
        
        return eff_fcnt( energy / energy_units );
      };
    }//case kFunctialEfficienyForm:


    case kExpOfLogPowerSeries:
    {
      if( m_expOfLogPowerSeriesCoeffs.empty() )
        return nullptr;

      const vector<float> coeffs = m_expOfLogPowerSeriesCoeffs;
      return [energy_units, coeffs, lowerEnergy, upperEnergy]( float energy ) -> float {
        if( (lowerEnergy > 10.0) && (energy < lowerEnergy) )
          return expOfLogPowerSeriesEfficiency( lowerEnergy / energy_units, coeffs );
        if( (upperEnergy > 10.0) && (energy > upperEnergy) )
          return expOfLogPowerSeriesEfficiency( upperEnergy / energy_units, coeffs );
        return expOfLogPowerSeriesEfficiency( energy / energy_units, coeffs );
      };
    }//case kExpOfLogPowerSeries:


    case kNumEfficiencyFnctForms:
      return nullptr;
  }//switch( m_efficiencyForm )

  return nullptr;
}//std::function<float( float )> DetectorPeakResponse::intrinsicEfficiencyFcn() const



float DetectorPeakResponse::peakResolutionFWHM( float energy,
                                                ResolutionFnctForm fcnFrm,
                                                const std::vector<float> &pars )
{
  switch( fcnFrm )
  {
    case kGadrasResolutionFcn:
    {
      if( pars.size() != 3 )
        throw std::runtime_error( "DetectorPeakResponse::peakResolutionSigma():"
                                 " pars not defined" );
      assert( PhysicalUnits::keV == 1.0 );
      
      /*
       // 20231223: a updated straight-forward translation from the fortran is:
       const float &resolutionOffset = pars[0];
       const float &resolution661 = pars[1];
       const float &resolutionPower = pars[2];
       energy = std::max( 1.0f, energy );
       if( energy > 661.0f )
         return 6.61f * resolution661 * std::pow( energy/661.0f, resolutionPower );
       
       if( resolutionOffset >= 0.0 )
       {
         const float ZeroLimit = std::max( 0.0f, fabs(resolutionOffset)*(661.0f - energy)/661.0f );
         const float FWHM = 6.61f * resolution661 * std::pow( energy/661.0f, resolutionPower );
         return std::sqrt( ZeroLimit*ZeroLimit + FWHM*FWHM );
       }
       
       const float p = std::pow(resolutionPower, (1.0f/std::log(1.0f - resolutionOffset)) );
       return 6.61f * resolution661 * std::pow(std::max(30.0f,energy)/661.0f, p );
       */
      
      const double a = pars[0];
      const double b = pars[1];
      const double c = pars[2];
      
      if( energy >= 661.0f || fabs(a)<float(1.0E-6) )
        return static_cast<float>( 6.61 * b * pow(energy/661.0, c) );

      if( a < 0.0 )
      {
        const double p = pow( c, 1.0/log(1.0-a) );
        return static_cast<float>( 6.61 * b * pow(energy/661.0, p) );
      }//if( a < 0.0 )
      
      if( a > 6.61*b )
        return static_cast<float>( a );

      const double A7 = sqrt( pow( 6.61*b, 2.0) - a*a ) / 6.61;
      return static_cast<float>( sqrt(a*a + pow( 6.61 * A7 * pow(energy/661.0, c), 2.0)) );
    }//case kGadrasResolutionFcn:
    
    case kSqrtEnergyPlusInverse:
    {
      if( pars.size() != 3 )
        throw std::runtime_error( "DetectorPeakResponse::peakResolutionSigma():"
                                 " pars not defined" );
      energy /= PhysicalUnits::keV;
      
      return sqrt(pars[0] + pars[1]*energy + pars[2]/energy);
    }//case kSqrtEnergyPlusInverse:
      
    case kConstantPlusSqrtEnergy:
    {
      if( pars.size() != 2 )
        throw std::runtime_error( "DetectorPeakResponse::peakResolutionSigma():"
                                 " pars not defined" );
      energy /= PhysicalUnits::keV;
      
      return pars[0] + pars[1]*sqrt(energy);
    }//case kConstantPlusSqrtEnergy:
      
    case kSqrtPolynomial:
    {
      if( pars.size() < 1 )
        throw runtime_error( "DetectorPeakResponse::peakResolutionSigma():"
                             " pars not defined" );

      energy /= PhysicalUnits::MeV;
      //return  A1 + A2*std::pow( energy + A3*energy*energy, A4 );

      // Use Horner's method to evaluate the polynomial - more stable.
      double val = pars.back();
      for( int i = static_cast<int>(pars.size()) - 2; i >= 0; i -= 1 )
        val = val * energy + pars[i]; // Multiply by x and add the next coefficient

#ifndef NDEBUG
      double nonstable_val = pars[0];
      for( size_t i = 1; i < pars.size(); ++i )
        nonstable_val += pars[i] * pow(static_cast<double>(energy), static_cast<double>(i) );
      const double diff = fabs(nonstable_val - val);
      assert( (diff < 1.0E-4) || (diff < 1.0E-3*max(fabs(nonstable_val), fabs(val))) );
#endif

      return sqrt( val );
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


std::string DetectorPeakResponse::javaScriptFwhmFunction() const
{
  if( !hasResolutionInfo() )
  {
    throw std::runtime_error( "DetectorPeakResponse::javaScriptFwhmFunction(): "
                              "detector does not have resolution information" );
  }
  
  return javaScriptFwhmFunction( m_resolutionCoeffs, m_resolutionForm );
}//std::string javaScriptFwhmFunction() const


std::string DetectorPeakResponse::javaScriptFwhmFunction( const std::vector<float> &coeffs,
                                                          const ResolutionFnctForm form )
{
  std::string js_fwhm_fcn;
  
  switch( form )
  {
    case kGadrasResolutionFcn:
      if( coeffs.size() >= 3 )
      {
        js_fwhm_fcn = "function(e) { var a=" + std::to_string(coeffs[0]) +
        ", b=" + std::to_string(coeffs[1]) +
        ", c=" + std::to_string(coeffs[2]) +
        "; e = Math.max(1.0, e);" +
        " if(e > 661.0) return 6.61 * b * Math.pow(e/661.0, c);" +
        " if(a >= 0.0) { var ZeroLimit = Math.max(0.0, Math.abs(a)*(661.0-e)/661.0);" +
        " var FWHM = 6.61 * b * Math.pow(e/661.0, c);" +
        " return Math.sqrt(ZeroLimit*ZeroLimit + FWHM*FWHM); }" +
        " if(a > 6.61*b) return a;" +
        " var A7 = Math.sqrt(Math.pow(6.61*b, 2.0) - a*a) / 6.61;" +
        " return Math.sqrt(a*a + Math.pow(6.61 * A7 * Math.pow(e/661.0, c), 2.0)); }";
      }
      break;
      
    case kSqrtEnergyPlusInverse:
      if( coeffs.size() >= 3 )
      {
        js_fwhm_fcn = "function(e) { return Math.sqrt(" + std::to_string(coeffs[0]) +
        " + " + std::to_string(coeffs[1]) + "*e + " + std::to_string(coeffs[2]) + "/e); }";
      }
      break;
      
    case kConstantPlusSqrtEnergy:
      if( coeffs.size() >= 2 )
      {
        js_fwhm_fcn = "function(e) { return " + std::to_string(coeffs[0]) +
        " + " + std::to_string(coeffs[1]) + "*Math.sqrt(e); }";
      }
      break;
      
    case kSqrtPolynomial:
      if( !coeffs.empty() )
      {
        js_fwhm_fcn = "function(e) { e = e/1000.0; var val = " + std::to_string(coeffs.back()) + ";";
        for( int i = static_cast<int>(coeffs.size()) - 2; i >= 0; i-- )
        {
          js_fwhm_fcn += " val = val * e + " + std::to_string(coeffs[i]) + ";";
        }
        js_fwhm_fcn += " return Math.sqrt(val); }";
      }
      break;
      
    case kNumResolutionFnctForm:
      break;
  }
  
  if( js_fwhm_fcn.empty() )
  {
    throw std::runtime_error( "DetectorPeakResponse::javaScriptFwhmFunction(): "
                              "unable to generate JavaScript function for resolution form" );
  }
  
  return js_fwhm_fcn;
}//static std::string javaScriptFwhmFunction(...)


float DetectorPeakResponse::detectorDiameter() const
{
  return m_detectorDiameter;
}


void DetectorPeakResponse::printDetectorParameterizationToStdout() const
{
  using std::cout;
  using std::endl;
  
  cout << "\n========== DETECTOR PARAMETERIZATION ===========" << endl;
  cout << "// Copy and paste this into hardcoded detector creation code" << endl;
  cout << "Name: \"" << m_name << "\"" << endl;
  cout << "Description: \"" << m_description << "\"" << endl;
  
  // Detector diameter
  if( m_detectorDiameter > 0.0 )
  {
    const std::string diam_str = PhysicalUnits::printToBestLengthUnits( m_detectorDiameter );
    cout << "Detector Diameter: " << diam_str << " (" << m_detectorDiameter << " mm)" << endl;
  }
  
  // Geometry type
  cout << "Geometry Type: ";
  switch( m_geomType )
  {
    case EffGeometryType::FarField: cout << "FarField"; break;
    case EffGeometryType::FixedGeomTotalAct: cout << "FixedGeomTotalAct"; break;
    case EffGeometryType::FixedGeomActPerCm2: cout << "FixedGeomActPerCm2"; break;
    case EffGeometryType::FixedGeomActPerM2: cout << "FixedGeomActPerM2"; break;
    default: cout << "Unknown(" << static_cast<int>(m_geomType) << ")"; break;
  }
  cout << endl;
  
  // Efficiency information
  cout << "\n--- EFFICIENCY ---" << endl;
  cout << "Efficiency Function Type: ";
  switch( m_efficiencyForm )
  {
    case kEnergyEfficiencyPairs: cout << "EnergyEfficiencyPairs"; break;
    case kFunctialEfficienyForm: cout << "FunctionalEfficiencyForm"; break;
    case kExpOfLogPowerSeries: cout << "ExpOfLogPowerSeries"; break;
    case kNumEfficiencyFnctForms: cout << "NumEfficiencyFnctForms"; break;
  }
  cout << endl;
  
  cout << "Energy Units: " << m_efficiencyEnergyUnits;
  if( std::abs(m_efficiencyEnergyUnits - static_cast<float>(PhysicalUnits::keV)) < 0.001 )
    cout << " (keV)";
  else if( std::abs(m_efficiencyEnergyUnits - static_cast<float>(PhysicalUnits::MeV)) < 0.001 )
    cout << " (MeV)";
  cout << endl;
  
  if( m_efficiencyForm == kExpOfLogPowerSeries && !m_expOfLogPowerSeriesCoeffs.empty() )
  {
    cout << "ExpOfLogPowerSeries Coefficients: {";
    for( size_t i = 0; i < m_expOfLogPowerSeriesCoeffs.size(); ++i )
    {
      if( i > 0 ) cout << ", ";
      cout << std::scientific << std::setprecision(8) << m_expOfLogPowerSeriesCoeffs[i];
    }
    cout << "}" << endl;
    
    if( !m_expOfLogPowerSeriesUncerts.empty() )
    {
      cout << "ExpOfLogPowerSeries Uncertainties: {";
      for( size_t i = 0; i < m_expOfLogPowerSeriesUncerts.size(); ++i )
      {
        if( i > 0 ) cout << ", ";
        cout << std::scientific << std::setprecision(8) << m_expOfLogPowerSeriesUncerts[i];
      }
      cout << "}" << endl;
    }
  }
  
  if( m_efficiencyForm == kFunctialEfficienyForm && !m_efficiencyFormula.empty() )
  {
    cout << "Efficiency Formula: \"" << m_efficiencyFormula << "\"" << endl;
  }
  
  if( m_efficiencyForm == kEnergyEfficiencyPairs && !m_energyEfficiencies.empty() )
  {
    cout << "Energy-Efficiency Pairs (" << m_energyEfficiencies.size() << " points):" << endl;
    for( size_t i = 0; i < std::min(size_t(10), m_energyEfficiencies.size()); ++i )
    {
      cout << "  {" << m_energyEfficiencies[i].energy << ", " 
           << std::scientific << std::setprecision(6) << m_energyEfficiencies[i].efficiency << "}";
      if( i + 1 < std::min(size_t(10), m_energyEfficiencies.size()) ) cout << ",";
      cout << endl;
    }
    if( m_energyEfficiencies.size() > 10 )
      cout << "  ... (" << (m_energyEfficiencies.size() - 10) << " more points)" << endl;
    
    
    const int fcnOrder = 6;
    std::vector<float> fit_pars, fit_pars_incert;
    std::vector<MakeDrfFit::DetEffDataPoint> eff_points;
    for( const EnergyEfficiencyPair &eff : m_energyEfficiencies )
    {
      if( eff.energy < 45 || eff.energy > 3000 )
        continue;
      
      MakeDrfFit::DetEffDataPoint p;
      p.energy = eff.energy / 1000.0;
      p.efficiency = eff.efficiency;
      p.efficiency_uncert = 0.01*eff.efficiency;
      eff_points.push_back( p );
    }
    const double chi2 = MakeDrfFit::performEfficiencyFit( eff_points, fcnOrder, fit_pars, fit_pars_incert );
    
    cout << "ExpOfLogPowerSeries (chi2/dof=" << chi2 << "/" << (eff_points.size() - fcnOrder) << ") equivalent Coefficients: {";
    for( size_t i = 0; i < fit_pars.size(); ++i )
    {
      if( i > 0 ) cout << ", ";
      cout << std::scientific << std::setprecision(8) << fit_pars[i];
    }
    cout << "}" << endl;
    
    cout << "\n\nComparison to original:" << endl;
    double max_error = 0.0, min_error = 100000, average_error = 0.0, num_points = 0;
    for( const EnergyEfficiencyPair &eff : m_energyEfficiencies )
    {
      const double fit_eff = expOfLogPowerSeriesEfficiency( 0.001*eff.energy, fit_pars);
      const double percent_error = (100.0*fabs(eff.efficiency - fit_eff)/eff.efficiency);
      if( eff.energy >= 49 && eff.energy < 2650.0 )
      {
        num_points += 1;
        max_error = std::max( max_error, percent_error );
        min_error = std::min( min_error, percent_error );
        average_error += percent_error;
      }
      
      cout << "\tEnergy: " << std::fixed << eff.energy << " keV, OrigEff: " << std::scientific << eff.efficiency
      << ", FitEff: " << fit_eff
      << ", error: " << std::fixed << percent_error << "%" << endl;
    }
    average_error /= num_points;
    cout << "Between 50 and 2600 keV, average error is " << average_error << "%, with max error " << max_error << "%" << endl;
    cout << "----" << endl;
  }
  
  // Resolution (FWHM) information
  if( hasResolutionInfo() )
  {
    cout << "\n--- RESOLUTION (FWHM) ---" << endl;
    cout << "Resolution Function Type: ";
    switch( m_resolutionForm )
    {
      case kGadrasResolutionFcn: cout << "GadrasResolutionFcn"; break;
      case kSqrtPolynomial: cout << "SqrtPolynomial"; break;
      case kSqrtEnergyPlusInverse: cout << "SqrtEnergyPlusInverse"; break;
      default: cout << "Unknown(" << static_cast<int>(m_resolutionForm) << ")"; break;
    }
    cout << endl;
    
    if( !m_resolutionCoeffs.empty() )
    {
      cout << "Resolution Coefficients: {";
      for( size_t i = 0; i < m_resolutionCoeffs.size(); ++i )
      {
        if( i > 0 ) cout << ", ";
        cout << std::scientific << std::setprecision(8) << m_resolutionCoeffs[i];
      }
      cout << "}" << endl;
      
      if( !m_resolutionUncerts.empty() )
      {
        cout << "Resolution Uncertainties: {";
        for( size_t i = 0; i < m_resolutionUncerts.size(); ++i )
        {
          if( i > 0 ) cout << ", ";
          cout << std::scientific << std::setprecision(8) << m_resolutionUncerts[i];
        }
        cout << "}" << endl;
      }
    }
  }
  else
  {
    cout << "\n--- RESOLUTION (FWHM) ---" << endl;
    cout << "No resolution information available" << endl;
  }
  
  cout << "\n--- EXAMPLE CODE ---" << endl;
  cout << "// Create detector like this:" << endl;
  cout << "auto detector = std::make_shared<DetectorPeakResponse>(\"" << m_name << "\", \"" << m_description << "\");" << endl;
  
  if( m_efficiencyForm == kExpOfLogPowerSeries && !m_expOfLogPowerSeriesCoeffs.empty() )
  {
    cout << "std::vector<float> eff_coeffs = {";
    for( size_t i = 0; i < m_expOfLogPowerSeriesCoeffs.size(); ++i )
    {
      if( i > 0 ) cout << ", ";
      cout << m_expOfLogPowerSeriesCoeffs[i] << "f";
    }
    cout << "};" << endl;
    cout << "detector->fromExpOfLogPowerSeriesAbsEff(eff_coeffs, {}, 0.0f, " 
         << m_detectorDiameter << "f, " << m_efficiencyEnergyUnits << "f, 0.0f, 0.0f, "
         << "DetectorPeakResponse::EffGeometryType::FarField);" << endl;
  }
  
  if( hasResolutionInfo() && !m_resolutionCoeffs.empty() )
  {
    cout << "std::vector<float> res_coeffs = {";
    for( size_t i = 0; i < m_resolutionCoeffs.size(); ++i )
    {
      if( i > 0 ) cout << ", ";
      cout << m_resolutionCoeffs[i] << "f";
    }
    cout << "};" << endl;
    cout << "detector->setResolutionFunctionCoeffs(res_coeffs, DetectorPeakResponse::";
    switch( m_resolutionForm )
    {
      case kGadrasResolutionFcn: cout << "kGadrasResolutionFcn"; break;
      case kSqrtPolynomial: cout << "kSqrtPolynomial"; break;
      case kSqrtEnergyPlusInverse: cout << "kSqrtEnergyPlusInverse"; break;
      default: cout << "kGadrasResolutionFcn"; break;
    }
    cout << ");" << endl;
  }
  
  cout << "=================================================" << endl << endl;
}



std::shared_ptr<const DetectorPeakResponse> DetectorPeakResponse::getGenericHPGeDetector()
{
  static std::shared_ptr<const DetectorPeakResponse> s_detector;
  
  std::lock_guard<std::mutex> lock( s_generic_det_mutex );
  
  if( !s_detector )
  {
    auto detector = std::make_shared<DetectorPeakResponse>( "Generic HPGe", "" );
    
    // Eff coefficients from "ORTEC Detective-X_LANL_100cm (59%)"
    std::vector<float> eff_coeffs = {-1.58287716e+00, -7.39438057e-01, 3.77921872e-02, -4.59978022e-02, -4.18332592e-02};
    detector->fromExpOfLogPowerSeriesAbsEff( eff_coeffs, {}, 0.0f, 6.50*PhysicalUnits::cm, PhysicalUnits::MeV,
                                            55.0f, 3000.0f, EffGeometryType::FarField );
    
    // Use FWHM from the generic 40% HPGe GADRAS detector included with InterSpec
    std::vector<float> res_coeffs = {1.54999995e+00, 2.50000000e-01, 3.49999994e-01};
    
    detector->setFwhmCoefficients( res_coeffs, kGadrasResolutionFcn );
    detector->setDrfSource( DrfSource::UnknownDrfSource );
    
    s_detector = std::const_pointer_cast<const DetectorPeakResponse>( detector );
  }//if( !s_detector )
  
  return s_detector;
}


std::shared_ptr<const DetectorPeakResponse> DetectorPeakResponse::getGenericNaIDetector()
{
  static std::shared_ptr<const DetectorPeakResponse> s_detector;
  
  std::lock_guard<std::mutex> lock( s_generic_det_mutex );
  
  if( !s_detector )
  {
    auto detector = std::make_shared<DetectorPeakResponse>( "Generic NaI(Tl)", "12% Efficiency, 6.9% FWHM @661 keV" );
    
    // From "Canberra Inspector 1000 (12%)"
    vector<float> eff_coeffs = {-2.22847652e+00f, -1.89709997e+00f, -1.02779996e+00f, -4.63600010e-01f, -1.06799997e-01f};
    detector->fromExpOfLogPowerSeriesAbsEff(eff_coeffs, {}, 0.0f, 5.08*PhysicalUnits::cm, PhysicalUnits::MeV,
                                            45.0f, 3000.0f, DetectorPeakResponse::EffGeometryType::FarField );
    
    // FWHM from "RadSeeker-NaI"
    std::vector<float> res_coeffs = {-8.18999958e+00, 6.94000006e+00, 5.64000010e-01};
    
    detector->setFwhmCoefficients( res_coeffs, kGadrasResolutionFcn );
    detector->setDrfSource( DrfSource::UnknownDrfSource );
    
    s_detector = std::const_pointer_cast<const DetectorPeakResponse>( detector );
  }
  
  return s_detector;
}


std::shared_ptr<const DetectorPeakResponse> DetectorPeakResponse::getGenericLaBrDetector()
{
  static std::shared_ptr<const DetectorPeakResponse> s_detector;
  
  std::lock_guard<std::mutex> lock( s_generic_det_mutex );
  
  if( !s_detector )
  {
    auto detector = std::make_shared<DetectorPeakResponse>( "Generic LaBr3(Ce)", "" );
    
    // Fit to GADRAS efficiency of Sam-Eagle-LaBr
    // Between 50 and 2600 keV, average error is 1.40918552%, with max error 3.70071643%
    std::vector<float> eff_coeffs = {-1.66245103e+00, -1.07331991e+00, 9.46287289e-02, 3.29341553e-02, -7.69617930e-02, -1.86134428e-02};
    detector->fromExpOfLogPowerSeriesAbsEff( eff_coeffs, {}, 0.0f, 3.70*PhysicalUnits::cm, PhysicalUnits::MeV,
                                            50.0f, 2650.0f, EffGeometryType::FarField );
    
    // From GADRAS Sam-Eagle-LaBr
    std::vector<float> res_coeffs = {7.0, 2.6e+00, 5.2e-01};
    
    detector->setFwhmCoefficients( res_coeffs, kGadrasResolutionFcn /* TODO: verify resolution function type */ );
    detector->setDrfSource( DrfSource::UnknownDrfSource );
    
    s_detector = std::const_pointer_cast<const DetectorPeakResponse>( detector );
  }
  
  return s_detector;
}


std::shared_ptr<const DetectorPeakResponse> DetectorPeakResponse::getGenericCZTGeneralDetector()
{
  static std::shared_ptr<const DetectorPeakResponse> s_detector;
  
  std::lock_guard<std::mutex> lock( s_generic_det_mutex );
  
  if( !s_detector )
  {
    auto detector = std::make_shared<DetectorPeakResponse>( "Generic CZT", "" );
    
    // Fit to the GADRAS efficiencies between 45 and 3000 of the Kromek GR1 detector included with InterSpec
    //Between 50 and 2600 keV, average error is 1.36569567%, with max error 4.19282222%
    std::vector<float> eff_coeffs = {-3.36118150e+00, -1.55648577e+00, 2.62750953e-01, -2.36228734e-01, -2.88489103e-01, -5.32837957e-02};
    detector->fromExpOfLogPowerSeriesAbsEff( eff_coeffs, {}, 0.0f, 1.07*PhysicalUnits::cm, PhysicalUnits::MeV,
                                            50.0f, 2650.0f, EffGeometryType::FarField );
    
    // FWHM from "Kromek GR1"
    std::vector<float> res_coeffs = {8.94999981e+00f, 2.39000010e+00f, 3.44000012e-01f};
    
    detector->setFwhmCoefficients( res_coeffs, kGadrasResolutionFcn );
    detector->setDrfSource( DrfSource::UnknownDrfSource );
    
    s_detector = std::const_pointer_cast<const DetectorPeakResponse>( detector );
  }
  
  return s_detector;
}


std::shared_ptr<const DetectorPeakResponse> DetectorPeakResponse::getGenericCZTGoodDetector()
{
  static std::shared_ptr<const DetectorPeakResponse> s_detector;
  
  std::lock_guard<std::mutex> lock( s_generic_det_mutex );
  
  if( !s_detector )
  {
    auto detector = std::make_shared<DetectorPeakResponse>( "Generic CZT High-Quality", "" );
    
    // A fit to a GADRAS eff points for a M400
    // Between 50 and 2600 keV, average error is 1.50666823%, with max error 3.77653729%
    vector<float> eff_coeffs = {-3.35549235e+00, -1.39031291e+00, 2.49679491e-01, -2.98282117e-01, -2.85447180e-01, -4.70053926e-02};
    detector->fromExpOfLogPowerSeriesAbsEff( eff_coeffs, {}, 0.0f, 2.20*PhysicalUnits::cm, PhysicalUnits::MeV,
                                            50.0f, 2650.0f, EffGeometryType::FarField );
    // From M400 CZT
    vector<float> res_coeffs = {2.7, 0.699, 0.753};
    
    detector->setFwhmCoefficients( res_coeffs, kGadrasResolutionFcn /* TODO: verify resolution function type */ );
    detector->setDrfSource( DrfSource::DefaultGadrasDrf );
    
    s_detector = std::const_pointer_cast<const DetectorPeakResponse>( detector );
  }
  
  return s_detector;
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



void DetectorPeakResponse::fitResolution( DetectorPeakResponse::PeakInput_t peaks,
                    const std::shared_ptr<const Measurement> meas,
                    const DetectorPeakResponse::ResolutionFnctForm fnctnlForm )
{
  if( !peaks )
    throw runtime_error( "DetectorPeakResponse::fitResolution(...): invalid input" );
  
  int sqrtEqnOrder = static_cast<int>(peaks->size()) / 2;
  if( sqrtEqnOrder > 6 )
    sqrtEqnOrder = 6;
  else if( peaks->size() < 3 )
    sqrtEqnOrder = static_cast<int>( peaks->size() );
  
  vector<float> coefficients = m_resolutionCoeffs, uncerts;
  MakeDrfFit::performResolutionFit( peaks, fnctnlForm, sqrtEqnOrder, coefficients, uncerts );
  peaks = removeOutlyingWidthPeaks( peaks, fnctnlForm, coefficients );
  MakeDrfFit::performResolutionFit( peaks, fnctnlForm, sqrtEqnOrder, coefficients, uncerts );
  
  m_resolutionCoeffs = coefficients;
  m_resolutionForm = fnctnlForm;
  
  computeHash();
}//void fitResolution(...)


