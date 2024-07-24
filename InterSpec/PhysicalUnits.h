//Note these units are NOT CLHEP compatible!

#ifndef PhysicalUnits_h
#define PhysicalUnits_h 1
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

#include <vector>
#include <string>
#include <utility>

/*
  Note: This file evovled from the SimpleSimAnaTools project,
        and named "SimpleSimUnits.hh" but I copied it to this project,
        changed the definitions of second and eV, and have made several
        additional functions for input/output -  wcjohns 20120812
*/



namespace PhysicalUnits
{
  // Constants
  static const double pi  = 3.14159265358979323846;
  static const double mole = 1.0;
  static const double Avogadro = 6.02214179e+23/mole;

  // Length / areas
  static const double mm  = 1.;
  static const double cm  = 10.*mm;
  static const double meter  = 1000.0*mm;
  static const double m  = meter;
  static const double nm  = 1.e-9 *m;
  static const double um  = 1.e-6 *m;
  static const double angstrom  = 1.e-10*m;
  static const double fermi     = 1.e-15*m;
  static const double m2 = m*m;
  static const double cm2 = cm*cm;
  static const double cm3 = cm*cm*cm;
  static const double mm2 = mm*mm;
  static const double barn = 1.e-28*m2;
  static const double millibarn = 1.e-3 *barn;

  // Angle
  static const double radian = 1.0;
  static const double degree = pi*radian, deg =pi*radian;

  // Time
  static const double second = 1.0; // in CLHEP 1.e+9;  //a nano-second is basic unit of time
  static const double minute = 60.0 * second;
  static const double hour   = 60.0 * minute;
  static const double day    = 24.0 * hour;
  static const double year   = 365.2425 * day;  //Not int CLHEP!!
  static const double picosecond  = 1.0e-12 * second;
  static const double nanosecond  = 1.0e-9 * second;
  static const double microsecond = 1.0e-6 * second;
  static const double millisecond = 1.0e-3 * second;
  static const double month       = year / 12.0;

  // Charge
  static const double e_SI   = 1.602176487e-19;// positron charge in coulomb

  // Energy
  static const double MeV = 1.0E3;
  static const double  eV = 1.e-6*MeV;
  static const double keV = 1.e-3*MeV;
  static const double joule  = eV/e_SI;// joule = 6.24150 e+12 * MeV

  // Mass
  static const double  kilogram = joule*second*second/m2, kg = joule*second*second/m2;
  static const double      gram = 1.e-3*kg, g = 1.e-3*kg;
  static const double milligram = 1.e-3*gram, mg =  1.e-3*gram;

  // Activity
  static const double becquerel = 1./second ;
  static const double curie = 3.7e+10 * becquerel;
  static const double bq = 1.0 * becquerel;
  static const double kBq = 1.0e+3 * becquerel;
  static const double MBq = 1.0e+6 * becquerel;
  static const double GBq = 1.0e+9 * becquerel;
  static const double TBq = 1.0e+12 * becquerel;
  static const double fCi = 1.0e-15 * curie;
  static const double pCi = 1.0e-12 * curie;
  static const double nCi = 1.0e-9 * curie;
  static const double microCi = 1.0e-6 * curie;
  static const double mCi = 1.0e-3 * curie;
  static const double ci = curie;
  
  // Equivalent radiation dose
  static const double sievert = joule / kilogram;
  static const double rem = 0.01*sievert;

  // Absorbed dose
  static const double gray = joule / kilogram;
  static const double rad = 0.01 * gray;
  

  // Physical Constants
  static const double c_light   = 2.99792458e+8 * m/second;
  static const double c_squared = c_light * c_light;
  static const double amu_c2 = 931.494028 * MeV;
  static const double amu = amu_c2/c_squared;


  // TODO: consider making user literals, for example to write something like 25.0_cm, define:
  // long double operator"" _cm ( long double dist ){ return dist * PhysicalUnits::cm; }

  /** The number of Gaussian sigmas in a Full-Width at Half-maximum, for the normal distribution. */
  static const float fwhm_nsigma = 2.35482f;

  //sm_distanceRegex: a javascript regex string to validate user input distances
  //  stringToDistance.  Allows positive and negative values.
  extern const char * const sm_distanceRegex;
  extern const char * const sm_distanceUnitOptionalRegex;
  extern const char * const sm_distanceUncertaintyRegex;
  extern const char * const sm_distanceUncertaintyUnitsOptionalRegex;
  
  extern const char * const sm_activityRegex;
  extern const char * const sm_activityUnitOptionalRegex;

  extern const char * const sm_absorbedDoseRegex;
  extern const char * const sm_equivalentDoseRegex;

  extern const char * const sm_timeDurationRegex;
  extern const char * const sm_timeDurationHalfLiveOptionalRegex;
  extern const char * const sm_timeDurationHalfLiveOptionalPosOrNegRegex;
  
  extern const char * const sm_decimalRegex;
  extern const char * const sm_positiveDecimalRegex;

  //printToBest____Units(...): makes an attempt to turn the input into a user
  //  readable string.  The optional defineitions of units are convience
  //  variables for when you're working in units that are not from this
  //  namespace.
  //  TODO: instead of specifying number of digits past the decimal, should convert to using number of significant figures.
  std::string printToBestLengthUnits( double length,
                                      int maxNpostDecimal = 2,
                                      double cm_definition = cm );
  std::wstring printToBestLengthUnits( double length, double uncert,
                                      int maxNpostDecimal = 2,
                                      const double cm_definition = cm );
  std::string printToBestActivityUnits( double activity,
                                        int maxNpostDecimal = 2,
                                        bool useCuries = true, //or else use becquerel
                                        double bq_definition = becquerel );
  std::string printToBestTimeUnits( double time,
                                    int maxNpostDecimal = 2,
                                    double second_definition = second );
  std::string printToBestMassUnits( const double mass,
                                    const int maxNpostDecimal = 2,
                                    const double gram_definition = gram );
  
  std::string printToBestEquivalentDoseRateUnits( double dosePerTime,
                                   const int maxNpostDecimal = 2,
                                   const bool useSievertPerHr = true,
                                   const double svPerHrDefinition = sievert/hour );

  std::string printToBestEquivalentDoseUnits( double dose,
                                            const int maxNpostDecimal = 2,
                                            const bool use_sievert = true,
                                            const double sievert_definition = sievert );

  std::string printToBestAbsorbedDoseUnits( double dose,
                                          const int maxNpostDecimal = 2,
                                          const bool use_gray = true,
                                          const double gray_definition = gray );

  std::string printToBestSpecificActivityUnits( double activityPerMass,
                                                const int numSignificantFigures = 3,
                                                const bool useCurie = true );
  
  //stringToTimeDuration(...): Takes a user string, and attempts to interpret
  //  if as a time duration. Takes for example '5.2 y', '52.3 s', '00:01:2.1',
  //  '3.2d 15h', etc.
  // Also supports the extended format of ISO 8601, i.e., PnYnMnDTnHnMnS,
  //   Examples are: P1Y2M3DT10H30M, -P1347M, PT1M
  //
  //  throws std::runtime_error on failure
  //Note: time periods may be negative.
  //Note: the ISO 8601 format has only barely been tested as of 20190605
  double stringToTimeDuration( std::string str, double second_def = second );
  
  //stringToTimeDurationPossibleHalfLife(): similar to stringToTimeDuration()
  //  but time units can also be specified in half lifes (hl|halflife|halflives
  //  |half\\-life|half\\-lives|half\\slives|half\\slife), which also must be
  //  specified and greater than zero.
  double stringToTimeDurationPossibleHalfLife( std::string str,
                                               const double halflife,
                                                double second_def = second );

  //stringToActivity(...): Takes a user string, and attempts to interpret
  //  if as an activity.  Ex. '52.3 bq', '88 MCi', etc
  //  throws std::runtime_error on failure
  double stringToActivity( std::string str, double bq_def = becquerel );


  //stringToDistance(...): Accepts a floating point number followed by: m,
  //  meter, cm, mm, ft, feet, ', in, inches, and ".
  //  If more than one distance string is specified, then results of the
  //  individual string are summed, e.g.: '33.2ft 12 in' = 1042.42 cm, or
  //  '41.9 m 13cm' = 4203 cm.
  //  Uncertainties may be specified inside of of paranthesis, and are ignored,
  //  e.x. '41.3 (+-19.9) cm 1.0 (+-0.1) cm' returns 42.3 cm
  //Throws std::runtime_exception on failure.
  double stringToDistance( std::string str, double cm_definition = cm );

  /** Converts a string into absorbed dose.
   
   Accepted physical units: "gray", "Gy", "rad", "erg" (assumes per gram), "erg/g", "erg per gram"
   
   Example inputs: "1.2 gray", "1.23gy", "100 urad", "100 mrad", "100 milli-rad", "100 millirad",
                   "13 erg/g", "1 gray 1 Gy".
   
   Throws std::runtime_exception on failure.
   */
  double stringToAbsorbedDose( const std::string &str, const double gray_definition = gray );

  /** Converts a string into equivalent dose.
   
   Accepted physical units: "sievert", "Sv", "rem", "roentgen", "r\xC3\xB6entgen"
   
   Example inputs: "1.2 sv", "1.23 urem", "100 kilo-rem", "1 gigasv", "1gsv", "1.2mrem", (milli-rem)
                   "1.2Mrem" (mega-rem)
   
   Throws std::runtime_exception on failure.
   */
  double stringToEquivalentDose( const std::string &str, const double sievert_def = sievert );

  /** Converts a string into a mass.
   
   Currently only accepts units of grams (g), pounds (lb, lbs), ounces (oz), and grains, with prefixes of pico (p), nano (n), micro (u),
   milli (m), and kilo (k).
   
   Example valid input: "1.2 gram", "1.2 grams", "3g", "5 kg", "5 kilogram", "5 kilo-gram", "2lb", "2lb 3oz", "2lb - 3oz".

   Throws std::runtime_exception on failure.
   */
  double stringToMass( const std::string &str, const double gram_def = gram );


  /** Converts a string into an energy in keV.
    
    Currently only accepts units of eV, keV, or MeV.

    Example valid input: "185 keV", "1.011MeV", "2614 keV - 511 keV".

    Throws std::runtime_exception on failure.
  */
  double stringToEnergy( const std::string &str, const double keV_def = keV );


  /** Prints the passed in value and uncertainty in a form like "1.23 \xC2\xB1 0.123", where \xC2\xB1 is the plus-or-minus character.
   Output will use scientific notation decimals, depending on fewest number of characters (i.e., the '%g' printf flag).
     
   Currently, If value uncertainty and uncertainty "overlap" with more than one digit, then we will only print uncertainty out to the same
   decimal order as the value; if they only overlap by one digit, or not at all, then we will print out one more decimal point than value will
   be printed to.
   
   \TODO: Decide if above policy is a good way to print things out, and if so, finalize
   \TODO: trailing zeros may-or-may-not be appended after decimal point, this should be decided on
          and fixed
   
   Some example values (value,uncert,nsigfig) -> answer;
   - (4.66712, 0.0457612, 4) -> '4.667 ± 0.046'
   - (4.66712e-16, 4.57612e-18, 4) -> '4.667e-16 ± 4.6e-18'
   - (4.66712, 0.0457612, 3) -> '4.67 ± 0.046'
   - (4.66712, 0.0457612, 2) -> '4.7 ± 0.05'
   - (4.66712, 0.0457612, 1) -> '5 ± 0'
   - (0.1, 0.001, 2) -> '0.1 ± 0.001'
   - (0.1, 1.001, 3) -> '0.1 ± 1'
   - (0.1, 1.001, 4) -> '0.1 ± 1.001'
   - (1.23457, 0.897654, 1) -> '1 ± 0.9'
   - (1.23457, 0.897654, 2) -> '1.2 ± 0.9'
   - (1.23457, 0.897654, 3) -> '1.23 ± 0.9'
   - (1.23457, 0.897654, 4) -> '1.235 ± 0.898'
   - (1.23457, 8.97654e-06, 1) -> '1 ± 0'
   */
  std::string printValueWithUncertainty( double value, double uncert, size_t nsigfig );

  /** A simple passthrough to `SpecUtils::printCompact`, for compatibility. */
  [[deprecated("Use SpecUtils::printCompact(double,size_t).")]]
  std::string printCompact( const double value, const size_t sig_figs );

  //Functions below here were imported 20121014 from some other code I have
  //  written, and are not yet consistent with the above functions (e.g.
  //  will have different available prefixes on the units and such)
  typedef std::pair<std::string,double>    UnitNameValuePair;
  typedef std::vector< UnitNameValuePair > UnitNameValuePairV;

  extern const UnitNameValuePairV sm_lengthUnitNameValues;
  extern const UnitNameValuePairV sm_activityUnitNameValues;
  extern const UnitNameValuePairV sm_lengthUnitHtmlNameValues;
  extern const UnitNameValuePairV sm_activityUnitHtmlNameValues;
  extern const UnitNameValuePairV sm_doseRateUnitHtmlNameValues;

  //Returns a reference to an element in sm_activityUnitNameValues.
  const UnitNameValuePair &bestActivityUnit( const double activity,
                                             bool useCuries = true );
  
  //Returns a reference to an element in sm_activityUnitHtmlNameValues
  const UnitNameValuePair &bestActivityUnitHtml( const double activity,
                                                 bool useCuries = true );
  
  const UnitNameValuePair &bestDoseUnitHtml( const double dose,
                                                bool useRem = true );
  
  const UnitNameValuePair &bestLengthUnit( const double length );
  const UnitNameValuePair &bestLengthUnitHtml( const double length );
  
  UnitNameValuePair bestTimeUnit( double time );
}//namespace PhysicalUnits
#endif
