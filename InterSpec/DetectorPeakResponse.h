#ifndef DetectorPeakResponse_h
#define DetectorPeakResponse_h
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

#include <deque>
#include <cmath>
#include <thread>
#include <memory>
#include <cctype>
#include <string>
#include <vector>
#include <cstring>
#include <istream>
#include <sstream>
#include <utility>
#include <sstream>
#include <stdexcept>
#include <functional>

#include <Wt/Dbo/Dbo>
#include <Wt/WDateTime>
#include <Wt/Dbo/SqlTraits>
#include <Wt/Dbo/WtSqlTraits>

#include "InterSpec/PeakModel.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/WarningWidget.h"
#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/DetectorPeakResponse.h"


class PeakDef;
class PeakModel;
class Measurement;
class InterSpecUser;

namespace rapidxml
{
  template<class Ch> class xml_node;
  template<class Ch> class xml_document;
}//namespace rapidxml

//Three expression evaluators have been tried; all have benefits and drawbacks.
//  CLHEP/Evaluator: only added 21 kb to final binary size, but slowest option
//  exprtk.hpp: added 5.3 Mb to final binary size, but 30x faster than CLHEP
//  muparserx: adds 430 kb to binary size, and about 10x faster than CLHEP
//
// Since muparserx is what TerminalWidget uses, it was decided on 20180408 to
// only use muparserx.

//100000 evaluations of:
//  exp(-343.6330974237 + 269.1023287277*log(x) + -83.8077567526*log(x)^2
//  + 12.9980559362*log(x)^3 + -1.0068649823*log(x)^4 + 0.0311640084*log(x)^5)
//  muparser x took:  cpu=0.143294s, wall=0.14341s  (1.4 us/eval)
//  evaluator x took: cpu=1.10172s, wall=1.10209s   (11 us/eval)

//Forward declarations
namespace mup
{
  class Value;
  class ParserX;
}

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


/*
 //ToDo: Add table to database to track which DRF to use for a given detector
 //      model, or serial number,
 struct DrfToUsePreference
 {
 Wt::Dbo::ptr<InterSpecUser> user;
 long int / <DetectorPeakResponse>  //Dbo::weak_ptr<DetectorPeakResponse>
 enum MatchType{ MatchByDetectorModel, MatchBySerialNumber };
 MatchType match_type;
 string detector_dentifier;
 
 //could add info such as detector name, or whatever
 }
 */

class DetectorPeakResponse
{
  /*
   //An example of defining a detectors efficiency via functional form:
   DetectorPeakResponse det( "DetectorName", "Detecotor Description" );
   
   const std::string fcn = "exp(-343.6330974237 + 269.1023287277*log(x)"
                           "+ -83.8077567526*log(x)^2  + 12.9980559362*log(x)^3"
                           "+ -1.0068649823*log(x)^4 + 0.0311640084*log(x)^5)";
   try
   {
     float detector_diameter = PhysicalUnits::stringToDistance( "2.2 cm" );
     det.setIntrinsicEfficiencyFormula( fcn, detector_diameter, 1.0f, 0.0f, 0.0f );
   }catch( std::exception &e )
   {
     std::cerr << "Caught: " << e.what() << std::endl;
     return;
   }//try / catch
   
   std::cout << det.intrinsicEfficiency(121.78f) << std::endl; //0.625191
   std::cout << det.intrinsicEfficiency(411.02f) << std::endl; //0.333307
   std::cout << det.intrinsicEfficiency(500.0f) << std::endl;    //0.285503
   std::cout << det.intrinsicEfficiency(700.0f) << std::endl;    //0.219004
   std::cout << det.intrinsicEfficiency(800.0f) << std::endl;    //0.197117
  */

  
  
public:
  enum ResolutionFnctForm
  {
    kGadrasResolutionFcn, //See peakResolutionFWHM() implementation
    kSqrtPolynomial,  //FWHM = sqrt( Sum_i{A_i*pow(x/1000,i)} );  //previously implemented as A1 + A2*pow( energy/1E3 + A3*energy*energy/1E6, A4 ); Implementation not finalized.
    kNumResolutionFnctForm
  };//enum ResolutionFnctForm

  
  enum EfficiencyFnctForm
  {
    kEnergyEfficiencyPairs,
    kFunctialEfficienyForm,
    kExpOfLogPowerSeries,
    kNumEfficiencyFnctForms
  };//enum EfficiencyFnctForm
  
  /** Enum used to indicate where the DRF came from.  This is used primarily to
      help decide what DRFs to show when user browses database.
   */
  enum DrfSource
  {
    UnknownDrfSource = 4,
    
    /** One of the GADRAS based DRFs that come with InterSpec. */
    DefaultGadrasDrf = 0,
    
    /** A GADRAS DRF from the filesystem (like C:\Gadras\Detectors, or
     InterSpecs user data directory.
     */
    UserAddedGadrasDrf = 5,
    
    /** Relative (or intrinsic) efficiency DRF from a CSV or TSV from a user
     specified directory, like InterSpecs user data directory.
     */
    UserAddedRelativeEfficiencyDrf = 1,
    
    /** User used the "Import" of Detector Select tool to upload Efficiency.csv
     file, and specified the detector diameter.
    */
    UserImportedIntrisicEfficiencyDrf = 2,
    
    /** User used the "Import" of DetectorSelect tool to upload an
        Efficiency.csv and Detector.dat file.
     */
    UserImportedGadrasDrf = 6,
    
    /** User specified a formula in the Detector Select Tool. */
    UserSpecifiedFormulaDrf = 3,
    
    /** User used the MakeDrf tool to create the DRF. */
    UserCreatedDrf = 7,
    
    /** The DRF was included in a spectrum (N42) file. */
    FromSpectrumFileDrf = 8
  };//enum DrfSource
  
public:
  //DetectorPeakResponse(): constructs with no defined resolution, efficiency,
  //  or description; detector name is "DetectorPeakResponse".
  DetectorPeakResponse();
  
  
  //DetectorPeakResponse(...): constructs with no defined resolution or
  //  efficiency, but with given name and description.
  DetectorPeakResponse( const std::string &name,
                        const std::string &descrip = "" );
  
  
  
  //Destructor
  virtual ~DetectorPeakResponse();
  
  //The hash value gives a fast/unique way to compare detectors.
  uint64_t hashValue() const;
  uint64_t parentHashValue() const;
  void setParentHashValue( const uint64_t val );

  //isValid(): Tells wether or not this DetectorPeakResponse has been properly
  //  defined and can be used.  Right now this is solely based on whether or
  //  not 'm_energyEfficiencies' has been populated.
  bool isValid() const;
  
  
  //hasResolutionInfo(): Tells whether or not this DetectorPeakResponse has
  //  (valid) detector resolution information
  bool hasResolutionInfo() const;
  

  //reset(): clears all data, and causes isValid() to return false
  void reset();


  //operator==: compares all member varialbles.  Currently doesnt account for
  //  small precision errors on floating point variables.
  bool operator==( const DetectorPeakResponse &rhs ) const;


  //fromEnergyEfficiencyCsv(...): first field of each line should be centroid,
  //  and second the absolute efficiency (so efficiency of gammas hitting
  //  detector face) in percents, of the detector.
  //  All lines with non-digit first non-whitespace characters, will be ignored,
  //  as will lines with only one number on them.
  //  Input can be space, comma, or tab seperated.
  //  Will throw std::runtime_error() on input error.
  //  Detector diameter is in units of PhysicalUnits
  //  Recomputes hash value.
  void fromEnergyEfficiencyCsv( std::istream &input,
                const float detectorDiameter, const float energyUnits );
  
  
  //setIntrinsicEfficiencyFormula(): sets m_efficiencyForm, m_efficiencyFormula,
  //  and m_efficiencyFcn to correspond to the fucntional form passed in.
  //  Formula should be in the a general functional form, where energy is
  //  represented by the letter 'x' and is in keV.  The formula is for the
  //  intrinsic efficiency of the detector, and should range from 0.0 to 1.0,
  //  in the valid energy range.
  //  Recomputes hash value.
  //e.x.
  //  fcn = "exp(-343.6330974237 + 269.1023287277*log(x)"
  //        "+ -83.8077567526*log(x)^2  + 12.9980559362*log(x)^3"
  //        "+ -1.0068649823*log(x)^4 + 0.0311640084*log(x)^5)";
  //Should give:
  //  intrinsicEfficiency(121.78)==0.625191;
  //  intrinsicEfficiency(411.02)==0.333307;
  //  intrinsicEfficiency(500)   ==0.285503;
  //  intrinsicEfficiency(700)   ==0.219004;
  //  intrinsicEfficiency(800)   ==0.197117
  //
  //The 'detector_diameter' is in units of PhysicalUnits (e.g. mm).
  //
  //The lower and upper energy parameters specify the energy range the function
  //  is valid over; if unknown specify 0.0f for both.
  //
  //Throws std::runtime_error if invalid expression, with a message that is
  //  valid HTML
  //Note, valid function names are:
  //   abs min, max, sqrt, pow, sin, cos, tan, asin, acos, atan, atan2,
  //    sinh, cosh, tanh, exp, log, ln (synonym for log), log10
  void setIntrinsicEfficiencyFormula( const std::string &fcn,
                                      const float detector_diameter,
                                      const float eqnEnergyUnits,
                                      const float lowerEnergy,
                                      const float upperEnergy );
  
  //fromGadrasDefinition(...): accepts the Efficiency.csv and Detector.dat
  //  files from GADRAS to define the detector.
  //Note that just the Detector.dat file contains enough info to define the
  //  detector, so I could try to only accept this one file, but this isnt
  //  too important right now, since Efficiency.csv is also provided.
  //  Recomputes hash value.
  void fromGadrasDefinition( std::istream &efficiencyCsvFile,
                             std::istream &detDatFile );

  /** Convience function that calls #fromGadrasDefinition with the
      Efficiency.csv and Detector.dat files in the specified directory.
      Throws exception on issue.
   */
  void fromGadrasDirectory( const std::string &dir );
  
  
  /** Sets the detectors as a kExpOfLogPowerSeries effiency detector, using the
      fit absolute efficiency coefficients.
   
   \param coefs The coefficients for the exp(A0 + A1*log(x) + A2*log(x)^2...)
          equation.  If coefficients are not for intrinsic efficiency (e.g.,
          characterizationDist is not 0.0), they coeffcients will be converted
          to intrinsic for internal use
   \param characterizationDist Distance used when fitting coefficents.  If
          coefficients are for intrinsic efficiency, then this value will be
          0.0f.
   \param equationEnergyUnits The energy units equation should be evalueated
          in. If equn is in MeV, this will be 1000.  If in keV, will be 1.0.
   \param lowerEnergy Lower energy, in keV, the equation is good for; if
          unknown pass a value of 0.0f for this parameter and upperEnergy.
   \param upperEnergy Upper energy, in keV, the equation is good for; if
          unknown pass a value of 0.0f for this parameter and lowerEnergy.
   
   Recomputes hash value.
  */
  void fromExpOfLogPowerSeriesAbsEff( const std::vector<float> &coefs,
                                      const float characterizationDist,
                                      const float detector_diameter,
                                      const float equationEnergyUnits,
                                      const float lowerEnergy,
                                      const float upperEnergy );

  /**
   if form==kGadrasResolutionFcn then coefs must have 3 entries
   if form==kSqrtPolynomial then coefs must not be empty,
      and coefficients must have been fit for energy in MeV
   if form==kNumResolutionFnctForm then coefs must be empty
   */
  void setFwhmCoefficients( const std::vector<float> &coefs,
                            const ResolutionFnctForm form );
  
  //efficiency(...): Currently just linearly interpolates between surrounding
  //  efficiency points for m_efficiencyForm=kEnergyEfficiencyPairs.  Will throw
  //  std::runtime_exception if the object has not been initialized.  Above or
  //  below maximum energies of the efficiency will return upper or lower
  //  efficiencies, respectively.
  //Returns efficiency per decay measured at `distance`
  //Energy should be in units of SandiaDecay (e.g. keV=1.0).
  float efficiency( const float energy, const float distance ) const;


  //intrinsicEfficiency(...): Currently just linearly interpolates between
  //  surrounding efficiency points.  Will throw std::runtime_exception if the
  //  object has not been initialized.  Above or below maximum energies of the
  //  efficiency will return upper or lower efficiencies, respectively.
  float intrinsicEfficiency( const float energy ) const;


  //fractionalSolidAngle(...) returns the fraction of gamma rays from a point
  //  source that would strike the detector face of a detector with diameter
  //  'detectorDiameter' at a distance from source of distance.
  static float fractionalSolidAngle( const float detector_diameter,
                                     const float observation_distance );


  //fractionalSolidAngle(...): similar to the above, but takes into account
  //  the source radius (assumed to be flat round plane); see pg 119 in Knoll
  //  for details on approxiation used.
  static float fractionalSolidAngle( const float detector_diameter,
                                     const float observation_distance,
                                     const float source_radius );


  //peakResolutionFWHM(...): returns the full width at half max of the detector.
  //  If the m_resolutionCoeffs are not defined, an exception is thrown.
  float peakResolutionFWHM( const float energy ) const;
  static float peakResolutionFWHM( float energy,
                                   ResolutionFnctForm fcnFrm,
                                   const std::vector<float> &pars );

  
  //peakResolutionSigma(...): returns the resolution sigma of the detector,
  //  that is FWHM/2.35482.
  //  If the m_resolutionCoeffs are not defined, an exception is thrown.
  float peakResolutionSigma( const float energy ) const;
  static float peakResolutionSigma( const float energy,
                                    ResolutionFnctForm fcnFrm,
                                    const std::vector<float> &pars );

  
  //Simple accessors
  float detectorDiameter() const;
  const std::string &efficiencyFormula() const;
  const std::string &name() const;
  const std::string &description() const;
  DrfSource drfSource() const;
  
  float efficiencyEnergyUnits() const;
  
  //Simple setters (all recompute hash value)
  void setName( const std::string &name );
  void setDescription( const std::string &descrip );

  //Search for: setFwhmCoefficients, fromGadrasDirectory, fromGadrasDefinition, setIntrinsicEfficiencyFormula, fromEnergyEfficiencyCsv, fromExpOfLogPowerSeriesAbsEff
  //  And maybe consider making them take a DrfSource argument
  void setDrfSource( const DrfSource source );
  
  /** Sets the energy range the DRF is valid over; currently not enforced or
      indicated anywhere in InterSpec, but may be in the future.
   */
  void setEnergyRange( const double lower, const double upper );
  
  /** Updated the #m_lastUsedUtc member variable to current time.  Does not
      save to database.
   */
  void updateLastUsedTimeToNow();
  
  //Some temporary accessors for debugging 2019050
  ResolutionFnctForm resolutionFcnType() const { return m_resolutionForm; }
  const std::vector<float> &resolutionFcnCoefficients() const { return m_resolutionCoeffs; }
  EfficiencyFnctForm efficiencyFcnType() const { return m_efficiencyForm; }
  //std::vector<EnergyEfficiencyPair> m_energyEfficiencies;
  //std::string m_efficiencyFormula;
  //std::function<float(float)> m_efficiencyFcn;
  const std::vector<float> &efficiencyExpOfLogsCoeffs() const { return m_expOfLogPowerSeriesCoeffs; }
  
  
  
  //Some methods to fit detector resolution from data.
  //    Only kinda tested as of 20130428
  typedef std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > PeakInput_t;
  
  
  //fitResolution(...): performs a fit to the passed in peaks, and assigns the
  //  found coeffficients to this DetectorPeakResponse.  Note that this funtion
  //  may iterately discard outlier peaks, who contribute more than 2.5 times
  //  the mean per-peak contribution, down to a minum of 80% of previous peaks,
  //  or a minum of 5 peaks; currently a max of one iteration will be performed,
  //  but this may be changed in the future.
  //  The provided Measurement should be same one used to fit the peaks for.
  //Re-computes hash as well.
  //Throws exception on error or failure.
  void fitResolution( PeakInput_t peaks,
                      const std::shared_ptr<const Measurement> meas,
                      const ResolutionFnctForm fnctnlForm );
  
  //expOfLogPowerSeriesEfficiency(...): evalutaes absolute efficiency for
  // coefficients of the form:
  //  x=log(energy);  //natural log of energy in keV
  //  abs_eff = exp(coefs[0]+coefs[1]*x+coefs[2]*x^2+coefs[3]*x^3+coefs[4]*x^4+...)
  //It is expected energy is in the same units (e.g. keV, MeV, etc) as the coefs
  static float expOfLogPowerSeriesEfficiency( const float energy,
                                              const std::vector<float> &coefs );
  
  
  void toXml( ::rapidxml::xml_node<char> *parent, 
              ::rapidxml::xml_document<char> *doc ) const;
  void fromXml( const ::rapidxml::xml_node<char> *parent );
  
#if( PERFORM_DEVELOPER_CHECKS )
  //equalEnough(...): tests whether the passed in Measurement objects are
  //  equal, for most intents and purposes.  Allows some small numerical
  //  rounding to occur.
  //Throws an std::exception with a brief explanaition when an issue is found.
  static void equalEnough( const DetectorPeakResponse &lhs,
                           const DetectorPeakResponse &rhs );
#endif

public:
  struct EnergyEfficiencyPair
  {
    float energy;
    float efficiency;
    bool operator<( const EnergyEfficiencyPair &rhs ) const
    { return energy < rhs.energy; }
    bool operator<( const float rhs ) const
    { return energy < rhs; }
    bool operator==( const EnergyEfficiencyPair &rhs ) const
    { return energy==rhs.energy && efficiency==rhs.efficiency; }
  };//struct EnergyEfficiencyPair

  //Returns the energy efficiency pairs
  const std::vector<EnergyEfficiencyPair> &getEnergyEfficiencyPair() const;

  static float akimaInterpolate( const float energy,
                                const std::vector<EnergyEfficiencyPair> &xy );
  
  //20190525: Why is m_user a raw index, and not a Wt::Dbo::ptr<InterSpecUser>?
  //          Maybe for database upgrade so Wt::Dbo doesnt need a reference to
  //          InterSpecUser?
  int m_user;
  
protected:
  void computeHash();
  
protected:
  //intrinsicEfficiencyFrom...(...) functions assume energy is input in keV
  float intrinsicEfficiencyFromPairs( float energy ) const;
  float intrinsicEfficiencyFromFcn( float energy ) const;
  float intrinsicEfficiencyFromExpLnEqn( float energy ) const;
  
  std::string m_name;
  std::string m_description;

  //m_detectorDiameter: assumes a round detector face, if other shape you will
  //  need to find the equivalent diameter for the face surface area of
  //  detector
  float m_detectorDiameter;

  //m_efficiencyEnergyUnits: units the absolute energy efficinecy formula,
  //  equation, or EnergyEfficiencyPairs are expecting.  Defaults to
  //  PhysicalUnits::keV
  float m_efficiencyEnergyUnits;

  ResolutionFnctForm m_resolutionForm;
  std::vector<float> m_resolutionCoeffs;
  /** Valid only if same size as m_resolutionCoeffs. */
  std::vector<float> m_resolutionUncerts;
  
  DrfSource m_efficiencySource;
  
  //
  EfficiencyFnctForm m_efficiencyForm;
  
  //Design decision: I dont like have member variables that are only used for
  //  certain m_efficiencyForm types, and would have preffered to just have a
  //  single std::function<double(double)> object that would abstract awway
  //  the differences, however, the ability to serialize the detector response
  //  function easily, made this difficult to accomplish, so I'm doing it the
  //  bone-headed way I'm not entirely satisfied with at the moment.
  
  //m_energyEfficiencies: the raw intrinsic energy to efficincy pairs, only
  //  filled out if (m_efficiencyForm==kEnergyEfficiencyPairs), which also
  //  implies the efficiency came from a CSV file.
  std::vector<EnergyEfficiencyPair> m_energyEfficiencies;
  
   
  //m_efficiencyFormula: the raw functional form for the intrinsic effeiciency,
  // e.x. "-343.6 + 269.1*ln(x) + -83.8*ln(x)^2  + 13.0*ln(x)^3".
  //  Only filled out if m_efficiencyForm==kFunctialEfficienyForm.
  // 'x' should is energy in keV, and function should vary between 0.0 and 1.0
  //  for valid energy range.
  std::string m_efficiencyFormula;
  
  //m_efficiencyFcn: the actual function that returns the intrinsic efficiency
  //  when m_efficiencyForm==kFunctialEfficienyForm
  std::function<float(float)> m_efficiencyFcn;
  
  //m_expOfLogPowerSeriesCoeffs: the coefficients for the intrinsic efficiency
  //  when m_efficiencyForm==kExpOfLogPowerSeries
  std::vector<float> m_expOfLogPowerSeriesCoeffs;
 
  /** Valid if same size as m_expOfLogPowerSeriesCoeffs. */
  std::vector<float> m_expOfLogPowerSeriesUncerts;
  
  //In order to keep track of lineage and uniqness of detectors, we will use
  //  hash values.  All non-serialization related non-const member functions
  //  should recompute the m_hash value when/if it makes any changes to the
  //  detector.
  uint64_t m_hash;
  
  //The m_parentHash is supposed to help track lineage of the detector.  If
  //  retrieve a detector from a database or spectrum file, and then alter and
  //  use it, you should set m_parentHash to the original value of m_hash you
  //  loaded the detector from.  This is not managed in this class; you must
  //  do this.
  //  (we may work out a better system once we actually implement modifying
  //   detectors in InterSpec)
  uint64_t m_parentHash;
  
  /** Not currently used, but in place for future upgrades, so DB schema wont
   have to be changed.
   */
  uint64_t m_flags;
  
  /** The lower energy (in keV) the DRF is good to; not (currently) enforced
   anywhere in InterSpec, but is good to know.
   */
  double m_lowerEnergy;
  
  /** The upper energy (in keV) the DRF is good to; not (currently) enforced
   anywhere in InterSpec, but is good to know.
   */
  double m_upperEnergy;
  
  /** Time when DRF was created. */
  int64_t m_createdUtc;
  
  /** Last time the DRF was used. */
  int64_t m_lastUsedUtc;
  
  static const int sm_xmlSerializationVersion;
  
public:
  
  template<class Action>
  void saveFloatVectorToDB( std::vector<float> &vFloat,
                            const std::string &dbname, Action &a )
  {
    try
    {
      std::stringstream ssv;
      for( size_t i = 0; i < vFloat.size(); ++i )
        ssv << (i?" ":"") << vFloat[i];
      std::string result = ssv.str();
      Wt::Dbo::field( a, result, dbname );
    }catch( std::exception & )
    {
      throw std::runtime_error( "Error saving field " + dbname
                                + " of DetectorPeakResponse to database." );
    } //catch
  }//saveFloatVectorToDB(...)
  
  template<class Action>
  void loadDBToFloatVector( std::vector<float> &vFloat,
                            const std::string &dbname, Action &a )
  {
    vFloat.clear();
    
    std::string result;
    Wt::Dbo::field( a, result, dbname );
    
    UtilityFunctions::split_to_floats( result.c_str(), result.size(), vFloat );
  } //loadDBToFloatVector(std::vector<float> vFloat, std::string dbname, Action a)
  
  template<class Action>
  void saveEnergyEfficiencyPairVectorToDB(
                                  std::vector<EnergyEfficiencyPair> &effPair,
                                  const std::string &dbname, Action &a )
  {
    try
    {
      std::stringstream ssv;
      for( size_t i = 0; i < effPair.size(); ++i )
        ssv << (i?" ":"") << effPair[i].energy << " " << effPair[i].efficiency;
      std::string result = ssv.str();
      Wt::Dbo::field( a, result, dbname );
    }catch( std::exception & )
    {
      throw std::runtime_error( "Error saving field " + dbname + " to database" );
    } //catch
  }//saveEnergyEfficiencyPairVectorToDB(...)
  
  
  template<class Action>
  void loadDBToEnergyEfficiencyPairVector(
                                    std::vector<EnergyEfficiencyPair> &effPair,
                                    const std::string &dbname, Action &a )
  {
    effPair.clear();
    
    std::string result;
    Wt::Dbo::field( a, result, dbname );
    
    std::vector<float> vals;
    UtilityFunctions::split_to_floats( result.c_str(), result.size(), vals );
    
    if( (vals.size()%2) != 0 )
      throw std::runtime_error( dbname + " field of DetectorPeakResponse "
                                "database entry had invalid number of fields" );
    
    for( size_t i = 0; i < vals.size(); i += 2 )
	{
	  EnergyEfficiencyPair val;
	  val.energy = vals[i];
	  val.efficiency = vals[i+1];
	  effPair.push_back( val );
      //effPair.push_back( EnergyEfficiencyPair{vals[i],vals[i+1]} );
	}
  }// void loadDBToEnergyEfficiencyPairVector(...)
  
  template<class Action>
  void persist( Action &a )
  {
    Wt::Dbo::field( a, m_name, "m_name", 255 );
    Wt::Dbo::field( a, m_description, "m_description", 255 );
    Wt::Dbo::field( a, m_detectorDiameter, "m_detectorDiameter" );
    Wt::Dbo::field( a, m_efficiencyEnergyUnits, "m_efficiencyEnergyUnits" );
    Wt::Dbo::field( a, m_resolutionForm, "m_resolutionForm" );
    
    if( a.getsValue() )
      saveFloatVectorToDB(m_resolutionCoeffs, "m_resolutionCoeffs", a);
    if( a.setsValue() || a.isSchema() )
      loadDBToFloatVector(m_resolutionCoeffs, "m_resolutionCoeffs", a);
    
    Wt::Dbo::field( a, m_efficiencySource, "m_efficiencySource" );
    Wt::Dbo::field( a, m_efficiencyForm, "m_efficiencyForm" );

    if( a.getsValue() )
      saveEnergyEfficiencyPairVectorToDB(m_energyEfficiencies, "m_energyEfficiencies", a);
    if( a.setsValue() || a.isSchema() )
      loadDBToEnergyEfficiencyPairVector(m_energyEfficiencies, "m_energyEfficiencies", a);
    
    Wt::Dbo::field( a, m_efficiencyFormula, "m_efficiencyFormula" );
    
    if( (a.setsValue() || a.isSchema()) && m_efficiencyFormula.size() )
    {
      const bool isMeV = (m_efficiencyEnergyUnits > 10.0f);
      try
      {
        std::shared_ptr<FormulaWrapper> expression = std::make_shared<FormulaWrapper>(m_efficiencyFormula,isMeV);
        m_efficiencyFcn = [expression](float a) -> float { return expression->efficiency(a); };
      }catch( std::exception & )
      {
        //In principle this shouldnt happen - in practice it might
        m_efficiencyFcn = std::function<float(float)>();
        if( m_efficiencyFormula.find( "invalid formula:" ) == std::string::npos )
          m_efficiencyFormula = "invalid formula: " + m_efficiencyFormula;
      }//try / catch
    }//if( reading from DB )
    
    if( a.getsValue() )
      saveFloatVectorToDB(m_expOfLogPowerSeriesCoeffs, "m_expOfLogPowerSeriesCoeffs", a);
    if( a.setsValue() || a.isSchema() )
      loadDBToFloatVector(m_expOfLogPowerSeriesCoeffs, "m_expOfLogPowerSeriesCoeffs", a);
    
    Wt::Dbo::field( a, m_user, "InterSpecUser_id" );
    
    //Wt::Dbo doesnt support unsigned integers, so we got a little workaround
    int64_t hash, parentHash, flags;
    if( a.getsValue() )
    {
      hash = reinterpret_cast<int64_t&>(m_hash);
      parentHash = reinterpret_cast<int64_t&>(m_parentHash);
      flags = reinterpret_cast<int64_t&>(m_flags);
    }
    
    
    Wt::Dbo::field( a, hash, "Hash" );
    Wt::Dbo::field( a, parentHash, "ParentHash" );
    Wt::Dbo::field( a, flags, "m_flags" );
    
    if( a.setsValue() || a.isSchema() )
    {
      m_hash = reinterpret_cast<uint64_t&>(hash);
      m_parentHash = reinterpret_cast<uint64_t&>(parentHash);
      m_flags = reinterpret_cast<uint64_t&>(flags);
    }

    if( a.getsValue() )
      saveFloatVectorToDB(m_expOfLogPowerSeriesUncerts, "m_expOfLogPowerSeriesUncerts", a);
    if( a.setsValue() || a.isSchema() )
      loadDBToFloatVector(m_expOfLogPowerSeriesUncerts, "m_expOfLogPowerSeriesUncerts", a);
    
    if( a.getsValue() )
      saveFloatVectorToDB(m_resolutionUncerts, "m_resolutionUncerts", a);
    if( a.setsValue() || a.isSchema() )
      loadDBToFloatVector(m_resolutionUncerts, "m_resolutionUncerts", a);
    
    Wt::Dbo::field( a, m_lowerEnergy, "m_lowerEnergy" );
    Wt::Dbo::field( a, m_upperEnergy, "m_upperEnergy" );
    Wt::Dbo::field( a, m_createdUtc, "m_createdUtc" );
    Wt::Dbo::field( a, m_lastUsedUtc, "m_lastUsedUtc" );
  } //void persist( Action &a )
};//class DetectorPeakResponse


#endif  //DetectorPeakResponse_h
