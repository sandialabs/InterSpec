#ifndef SpectrumDataStructs_h
#define SpectrumDataStructs_h
/* SpecUtils: a library to parse, save, and manipulate gamma spectrum data files.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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

#include "SpecUtils_config.h"

#if( SPECTRUM_DATA_STRUCTS_NO_BOOST_LIB )
#define BOOST_DATE_TIME_NO_LIB
#endif

#include <set>
#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <stdexcept>
#include <functional>

#include <boost/date_time/posix_time/posix_time.hpp>

#include "SpecUtils/SpecUtilsAsync.h"


/*
Shortcommings that wcjohns should be addressed
  -Many of the ICD1 fields possible are not checked for
    -comments for multiple different tags, ...
 -Energy calibration should be made its own object and shared between 
  Measurments.
 -Neutron meausruemtns should have their own live and real time
 -Should add a DetectorInfo object that Measurement objects point to and share.
   -Should add things like dimention and RadDetectorKindCode to this object
 -Should eliminate the <InterSpec:DetectorType> tag in written N42 2012 files.
 -Should consider adding explicit dead_time field 
 -Should add elevation and uncertainties to GPS coordinates
 -Should add detector and item orientation/positions
 -Should change Measurement::detector_type_ to detector_description_
 -Should consider removing measurment_operator and inspection and make location 
  part of detector location object
 -Should add in InstrumentToItemDistance and InstrumentToItemBearing to
 -There is a degeneracy in MeasurementInfo between: detector_type_,
  instrument_type_, manufacturer_, and instrument_model_ - so this should 
  be sorted out.
 -the generated UUID should maybe be more stable with respect to just the 
  spectroscopic information.
 -Should rename DetectorType to SystemType or DetectionSystemType
 -Should add in tag that indicates original file type, which will survive 
  serialization to N42-2012 and back in
 -For analysis result should add information on what isotopes where in the 
  alarm templates
 -Should add a Dose field to Measurement; see CountDose for a starting point
 -Need to reduce the compilation memory requirments to allow compiling on 
  devices with only 1 GB of ram.
 -Should break implementation up into many files (ex SpectrumDataStructs_pcf.cpp
  SpectrumDataStructs_2012N42.cpp, etc.
*/


enum ParserType
{
  k2006Icd1Parser,
  K2012ICD1Parser,
  kSpcParser,  //ascii or binary
  kGR135Parser,
  kPcfParser,
  kChnParser,
  kIaeaParser,
  kTxtOrCsvParser,
  kCanberraCnfParser,
  kTracsMpsParser,
  kAramParser,
  kSPMDailyFile, //SpectroscopicPortalMonitor
  kAmptekMca,
  kMicroRaider,
  kOrtecListMode,
  kAutoParser
};//enum ParserType


enum SpectrumType
{
  kForeground,
  kSecondForeground,
  kBackground
};//enum SpectrumType


enum SaveSpectrumAsType
{
  /** See #MeasurementInfo::write_txt for details. */
  kTxtSpectrumFile,
  
  /** See #MeasurementInfo::write_csv for details. */
  kCsvSpectrumFile,
  
  /** See #MeasurementInfo::write_pcf for details. */
  kPcfSpectrumFile,
  
  /** See #MeasurementInfo::write_2006_N42_xml for details. */
  kXmlSpectrumFile,
  
  /** See #MeasurementInfo::write_2012_N42 for details. */
  k2012N42SpectrumFile,
  
  /** See #MeasurementInfo::write_integer_chn for details. */
  kChnSpectrumFile,
  
  /** See #MeasurementInfo::write_binary_spc for details. */
  kBinaryIntSpcSpectrumFile,
  
  /** See #MeasurementInfo::write_binary_spc for details. */
  kBinaryFloatSpcSpectrumFile,
  
  /** See #MeasurementInfo::write_ascii_spc for details. */
  kAsciiSpcSpectrumFile,
  
  /** See #MeasurementInfo::write_binary_exploranium_gr130v0 for details. */
  kExploraniumGr130v0SpectrumFile,
  
  /** See #MeasurementInfo::write_binary_exploranium_gr135v2 for details. */
  kExploraniumGr135v2SpectrumFile,
  
  /** See #MeasurementInfo::write_iaea_spe for details. */
  kIaeaSpeSpectrumFile,

#if( ENABLE_D3_CHART_EXPORTING )
  /** See #MeasurementInfo::write_d3_html for details. */
  kD3HtmlSpectrumFile,
#endif
  
  kNumSaveSpectrumAsType
};//enum SaveSpectrumAsType


const char *descriptionText( const SpectrumType type );

//suggestedNameEnding(): returns suggested lowercase file name ending for type
//  passed in.  Does not contain the leading '.' for extentions
const char *suggestedNameEnding( const SaveSpectrumAsType type );


//spectrumTypeFromDescription(..): the inverse of descriptionText(SpectrumType)
//  throw runtiem_exception if doesnt match
SpectrumType spectrumTypeFromDescription( const char *descrip );

const char *descriptionText( const SaveSpectrumAsType type );


//DetectorType is intended to identify the detector system used to aquire the
//  spectra in a spectrum file.  It typically refers to a whole system which
//  may be comprised of multpile subdetectors, like for portals, or gamma
//  detectors with neutron detectors as well.
//Currently only opening files in MeasurementInfo that are specialezed (spc,
//  dat, etc) support filling out detector_type_, and even then I've been a
//  bit lazy and not added all systems, just the most common ones I work with.
enum DetectorType
{
  kGR135Detector,
  kIdentiFinderDetector,       //First gen identiFINDER with smaller crystal
                               //than NGs; note sometimes called identiFINDER-N.
                               // I dont have any examples of this
  kIdentiFinderNGDetector,     //Used for both the NG and NGH since same crystal
                               //  size (NGH has neutron tube)
  kIdentiFinderLaBr3Detector,  // Probably not detected (ever?)
  
  //The kDetectiveDetector is a default for when the type of detective cant be
  //  determined, not an actual detector type.  This enum doesnt consider the
  //  difference between the EX and DX series; the DX are same gamma crystal,
  //  but do not have a neutron detector.
  kDetectiveDetector,          
  kDetectiveExDetector,
  kDetectiveEx100Detector,
  kOrtecIDMPortalDetector,     //only identified from N42 files
  kSAIC8Detector,              //only identified from N42 files
  kFalcon5000,
  kUnknownDetector,
  kMicroDetectiveDetector,
  kMicroRaiderDetector,
  kRadHunterNaI,
  kRadHunterLaBr3,
  kRsi701,
  kRsi705,
  kAvidRsi, //unspecified RSI/Avid system, usually model is stated as RS???
  kOrtecRadEagleNai,
  kOrtecRadEagleCeBr2Inch,
  kOrtecRadEagleCeBr3Inch,
  kOrtecRadEagleLaBr,
  kSam940LaBr3,  //The LaBr3 may not always be detector, and then it will be assigned kSame940
  kSam940,
  kSam945
  //RadSeekerLaBr1.5x1.5
  //RadSeekerNaI2x2 (although should check for this, see MeasurementInfo::set_n42_2006_instrument_info_node_info
};//enum DetectorType



//A Forward declaration
namespace rapidxml
{
  template<class Ch> class xml_node;
  template<class Ch> class xml_document;
  template<class Ch> class xml_attribute;
}//namespace rapidxml

class InterSpec;
class SpecMeas;
class Measurement;
class MeasurementInfo;
class DetectorAnalysis;
struct MeasurementCalibInfo; //defined in SpectrumDataStructs.cpp (used for parsing N42 2006/2012 files and rebinning)
struct SpectrumNodeDecodeWorker;
struct GrossCountNodeDecodeWorker;
#if( ENABLE_D3_CHART_EXPORTING )
namespace D3SpectrumExport{ struct D3SpectrumChartOptions; }
#endif


//Some typedefs and enums used for decode_2012_N42_rad_measurment_node(...)
enum DetectionType{ GammaDetection, NeutronDetection, GammaAndNeutronDetection, OtherDetection };
typedef std::map<std::string,std::pair<DetectionType,std::string> > IdToDetectorType;
typedef std::map<std::string,MeasurementCalibInfo>                  DetectorToCalibInfo;


typedef std::vector< std::pair<float,float> > DeviationPairVec;
typedef std::shared_ptr< std::vector<float> > ShrdFVecPtr;
typedef std::shared_ptr< const std::vector<float> > ShrdConstFVecPtr;

typedef std::shared_ptr<Measurement> MeasurementShrdPtr;
typedef std::shared_ptr<const Measurement> MeasurementConstShrdPtr;

//gamma_integral(): get the integral of gamma counts between lowenergy and
//  upperunergy; linear approximation is used for fractions of channels.
double gamma_integral( const std::shared_ptr<const Measurement> &hist,
                       const float lowenergy, const float upperunergy );


// XXX - all these functions should be placed in a namespace, and
//      possible in a seperate source/header file
 
//expand_counted_zeros(...): requires zeros to be identically 0.0f, in order
//  be expanded by the next value.  The value following a zero is rounded to
//  nearest integer (no integer check is performed).
void expand_counted_zeros( const std::vector<float> &data,
                           std::vector<float> &results );

//compress_to_counted_zeros(...): contents less than 1E-8 are assumed to be
//  zeros
void compress_to_counted_zeros( const std::vector<float> &data,
                                std::vector<float> &results );

//polynomial_coef_to_fullrangefraction(..): only works for up to the first
//  four coefficints, as the fifth one for FRF doesnt translate easily.
std::vector<float> polynomial_coef_to_fullrangefraction(
                          const std::vector<float> &coeffs, const size_t nbin );

//fullrangefraction_coef_to_polynomial(..): only considers the first four
//  coefficients, as the fifth coefficient of FRF coorespods to a
//  a term like: C_4 / (1.0f+60.0*x)
std::vector<float> fullrangefraction_coef_to_polynomial(
                          const std::vector<float> &coeffs, const size_t nbin );


//mid_channel_polynomial_to_fullrangeFraction(): converts coefficients from
//  polynomial equation, which uses the convention the energy given by the
//  equation is the middle of the channel, to standard full range fraction.
std::vector<float> mid_channel_polynomial_to_fullrangeFraction(
                          const std::vector<float> &coeffs, const size_t nbin );

//polynomial_cal_remove_first_channels(): changes the polynomial calibration
//  coefficients, 'orig_coefs' to remove the first 'n' channels of the binning.
//  Would probably work for adding channels, but untested.
//Truncates binning to 6th order polynomial.
std::vector<float> polynomial_cal_remove_first_channels( const int n,
                                        const std::vector<float> &orig_coefs );

//polynomial_binning(...): returns lower channel energies from input polynomial
//  calibration equation and deviation pairs.  Note that this function uses
//  the convention that the energy of the lower edge of the channel is given by:
//  E_i = C_0 + i*C_1 + i*i*C_2 + ...
//Note that Canberra and Raytheon ASP systems may use the convention that
//  the energy given by the polynomial if for the bin mid-point.
std::shared_ptr< const std::vector<float> >
                 polynomial_binning( const std::vector<float> &coeffs,
                                     const size_t nbin,
                                     const DeviationPairVec &deviation_pairs );

//fullrangefraction_binning(...): returns lower channel energies from input
//  full width fraction calibration equation and deviation pairs.
//Uses the definition for the i'th channel:
//  x = i / nbin;
//  E_i = C_0 + x*C_1 + x*x*C_2 + x*x*x*C_3 + C_4/(1+60*x);
std::shared_ptr< const std::vector<float> >
                 fullrangefraction_binning( const std::vector<float> &coeffs,
                                            const size_t nbin,
                                            const DeviationPairVec &dev_pairs );

//fullrangefraction_energy(...): returns the energy cooresponding to the
//  passed in bin number; note that the bin_number may be non-integer. An
//  integer value gives you the energy of lower edge of that channel.
float fullrangefraction_energy( float bin_number,
                                const std::vector<float> &coeffs,
                                const size_t nbin,
                                const DeviationPairVec &deviation_pairs );

//find_bin_fullrangefraction(...):  returns the bin (including fractional
//  portion) corresponding to 'energy'.  If deviation_pairs is empty then a
//  algabraic approach is used, otherwise a binary search is performed to find
//  the bin that comes within 'accuracy' of 'energy'
float find_bin_fullrangefraction( const double energy,
                                  const std::vector<float> &coeffs,
                                  const size_t nbin,
                                  const DeviationPairVec &deviation_pairs,
                                  const double accuracy = 0.001 );

float offset_due_to_deviation_pairs( float original_energy,
                                     const DeviationPairVec &dps );

std::shared_ptr< const std::vector<float> > apply_deviation_pair(
                       std::shared_ptr< const std::vector<float> > binning,
                       const std::vector< std::pair<float,float> > &dev_pairs );


//rebin_by_lower_edge(...): not incredably well testted, but appears to be
//  better than the previous Measurment::rebin_by_lower_edge(...), but it
//  so it currently used.  There is some code in the function that
//  does test this function if PERFORM_DEVELOPER_CHECKS is true.
//  Also, a timing benchmark should be done between the alternate
//  implementations to decide which to use (although
//  Measurment::rebin_by_lower_edge does have some potential known bugs...)
//  There are some tests started for this function in
//  testing/testRebinByLowerEnergy.cpp.
void rebin_by_lower_edge( const std::vector<float> &original_energies,
                          const std::vector<float> &original_counts,
                          const std::vector<float> &new_energies,
                          std::vector<float> &resulting_counts );

/*
//returns _fractional_ bin energy belowngs in.
//Depreciated since it doesnt take into account DeviationPairs
float find_bin_from_polynomial( const float energy,
                                const std::vector<float> &coeffs,
                                const size_t nbin );
*/


//returns energy of _fractional_ channel number; an integer channel number
//  cooresponds to the lower edge of the channel.
float bin_number_to_energy_polynomial( const float bin,
                                       const std::vector<float> &coeffs,
                                       const size_t nbin );

int sample_num_from_remark( const std::string &remark ); //returns -1 on error
float speed_from_remark( std::string remark );
std::string detector_name_from_remark( const std::string &remark );


//time_duration_in_seconds(...): Reads times in formats similar to "PT16M44S"
//  or "13H82M49.33S".  Returns partial answer upon failure (and thus 0.0 on
//  complete failure), that is "PT16M44AS" would return 16 minutes, 0 seconds.
//  Shouldnt ever throw.
float time_duration_in_seconds( const std::string &duration );
float time_duration_in_seconds( const char *duration_str, const size_t length );

//dose_units_usvPerH(...): returns the dose units indicated by the string, in
//  units such that ia micro-sievert per hour is equal to 1.0.
//Currently only handles the label "uSv" and "uRem/h", e.g. fnctn not really
//  impleneted.
//Returns 0.0 on error.
float dose_units_usvPerH( const char *str, const size_t str_length );

//int detector_num_from_name( std::string name );

//Returns true if the file is likely a spectrum file, based off of file
//  extenstion, file size, etc..  By no means definitive, but useful when
//  looping through a large amount of files in order to filter out files likely
//  to not be spectrum files (but will also filter out a small amount of actual
//  spectrum files in practice).
bool likely_not_spec_file( const std::string &file );

//setAnalysisInformation(...): adds to analysis the information
//  in AnalysisResults node.  Currently only adds the NuclideAnalysis to
//  that is in identiFINDER files to analysis
void setAnalysisInformation( const rapidxml::xml_node<char> *analysis_node,
                             std::shared_ptr<DetectorAnalysis> analysis );


//detectorTypeToString(): returns string which cooresponds to the convention
//  InterSpec is using to represent detector response functions on disk.
const std::string &detectorTypeToString( const DetectorType type );

//We have to convert from 2006 N42 to 2012 instrument typs
const std::string &convert_n42_instrument_type_from_2006_to_2012(
                                                    const std::string &input );
    
/*
//TODO: start using (something like) this EnergyCalibration struct to
//      represent calibration.
struct EnergyCalibration
{
  enum CalibrationModel
  {
    Polynomial,
    FullRangeFraction,
    LowerChannelEdge,
    UnknownCalibrationModel
  };//From ICD1 Spec Polynomial Pade Exponential PolyLogarithmic

  CalibrationModel equation_type_;
  
  size_t nbin_;
  std::vector<float> coefficients_;
  std::vector< std::pair<float,float> > deviation_pairs_;
  std::shared_ptr< const std::vector<float> > binning_;
};//struct EnergyCalibration
*/

class Measurement
{
public:

  enum QualityStatus
  {
    //The detector status reported in the file; not applicable to all formats,
    //  in which case should be marked as Missing, although some formats
    //  (notable N42 and MPS) default to Good.
    Good, Suspect, Bad, Missing
  };

  enum OccupancyStatus
  {
    //Reported occupancy status; not be applicable to all systems/formats, in
    //  which case is marked to UnknownOccupancyStatus.
    NotOccupied, Occupied, UnknownOccupancyStatus
  };

  enum EquationType
  {
    //The energy (or FWHM) calibration type that the calibration coefficients
    //  should be interpreted as; typically also the type in the file.
    Polynomial, FullRangeFraction, LowerChannelEdge, UnknownEquationType
  };
  
  enum SourceType
  {
    //Reported source type for a record; marked as UnknownSourceType unless
    //  file format explicitly specifies, or can reasonably be infered.
    IntrinsicActivity, Calibration, Background, Foreground, UnknownSourceType
  };
  

  Measurement();
  
  //operator=: operates as expected for most member variables.  Smart pointer
  //  to const objects (channel data and channel energies) are shallow copied.
  const Measurement &operator=( const Measurement &rhs );
  
  //memmorysize(): calculates the approximate amount of memorry this Measurment
  //  is taking up in memmory, including all of the objects which it owns (like
  //  pointers to float arrays and stuff).
  size_t memmorysize() const;

  //Simple accessor functions (cheap to call):
  
  //live_time(): returned in units of seconds.  Will be 0 if not known.
  inline float live_time() const;
  
  //real_time(): returned in units of seconds.  Will be 0 if not known.
  inline float real_time() const;
  
  //contained_neutron(): returns whether or not the measurment is thought to
  //  contain the possibility to detect neutrons (e.g. if a neutron detector was
  //  also present).  This may be true even if zero neutrons were detected.  For
  //  some detector types this value is infered through previous hardware
  //  knowledge.
  inline bool contained_neutron() const;
  
  //sample_number(): the sample number assigned to this Measurment.  If the
  //  'DontChangeOrReorderSamples' flag wasnt passed to the
  //  MeasurementInfo::cleanup_after_load() function, then this value may be
  //  assigned during file parsing.
  inline int sample_number() const;
  
  //title(): some formats such as .PCF or .DAT files will contain a title for
  //  the spectrum, or this may be set through set_title(...).
  inline const std::string &title() const;
  
  //occupied(): returns the occupancy status.  Detectors which do not contain
  //  this capability will return 'UnknownOccupancyStatus'
  inline OccupancyStatus occupied() const;
  
  //gamma_count_sum(): returns the sum of channel data counts for gamma data.
  inline double gamma_count_sum() const;
  
  //neutron_counts_sum(): returns the sum of neutron counts.
  inline double neutron_counts_sum() const;
  
  //speed(): returns the speed of the vehicle, object or detector, in m/s if
  //  known.  Otherwise return 0.0.
  inline float speed() const;
  
  //latitude(): returns the latitude of the measurment, in degrees, if known.
  //  Returns -999.9 otherwise.
  inline double latitude() const;
  
  //longitude(): returns the longitude, in degrees, of the measurment if known.
  //  Returns -999.9 otherwise.
  inline double longitude() const;
  
  //has_gps_info(): returns true only if both latitude and longitude are valid.
  inline bool has_gps_info() const;
  
  //valid_longitude(): checks if abs(longitude) is less than or equal to 180.
  static bool valid_longitude( const double longitude );
  
  //valid_latitude(): checks if abs(latitude) is less or equal to 90.
  static bool valid_latitude( const double latitude );
  
  //position_time(): returns the (local, or detector) time of the GPS fix, if
  //  known.  Returns boost::posix_time::not_a_date_time otherwise.
  inline const boost::posix_time::ptime &position_time() const;
  
  //detector_name(): returns the name of the detector within the device.
  //  May be empty string for single detector systems, or otherwise.
  //  ex: Aa1, Ba1, etc.
  inline const std::string &detector_name() const;
  
  //detector_number(): returns the detector number of the detector within the
  //  detection system.  Will have a 1 to 1 coorespondence with detector_name().
  inline int detector_number() const;
  
  //detector_type():  If the file specifies the detector type string, it _may_
  //  be retrieved here.  Note there is not much consistency between file
  //  formats in what you should expect to get from this function.
  //  e.x. "HPGe 50%", "NaI"
  inline const std::string &detector_type() const;
  
  //quality_status():  If not specified in file, will have value of 'Missing'.
  inline QualityStatus quality_status() const;
  
  //source_type():  Returns the source type if known.  For some formats (notably
  //  PCF snd spectroscopic daily files), anything not background will be marked
  //  as foreground (this behaviour may be cahnged in the future).  For other
  //  formats if not known, 'UnknownSourceType' is returned.
  inline SourceType source_type() const;
  
  //energy_calibration_model(): returns calibration model used for energy
  //  binning.  If a value of 'UnknownEquationType' is returned, then
  //  channel_energies() may or may not return a valid pointer; otherwise, if
  //  this Measurment is part of a MeasurmentInfo object constructed by parsing
  //  a file, then channel_energies() pointer _should_ be valid.
  inline EquationType energy_calibration_model() const;
  
  //remarks(): the list of remarks found while parsing this record from the file
  //  that pertain to this record specifically.  See also
  //  MeasurmentInformation::remarks().
  inline const std::vector<std::string> &remarks() const;
  
  //start_time(): start time of the measurement.  Returns
  //  boost::posix_time::not_a_date_time if could not be determined.
  inline const boost::posix_time::ptime &start_time() const;

  //start_time_copy(): start time of the measurement.  Returns
  //  boost::posix_time::not_a_date_time if could not be determined.
  inline const boost::posix_time::ptime start_time_copy() const;
  
  //calibration_coeffs(): returns the energy calibration coeificients.
  //  Returned vector should have at least two elements for Polynomial and
  //  FullRangeFraction energy calibration models.  Polynomial may have an
  //  arbitrary number of coefficients, while FullRangeFraction may have up to
  //  five.  For LowerChannelEdge calibration model the returned vector will
  //  most likely be empty (a memory optimization, may be changed in the future)
  //  so you should instead call channel_energies().
  inline const std::vector<float> &calibration_coeffs() const;
  
  //deviation_pairs(): returns the energy deviation pairs.  Sometimes also
  // refered to as nonlinear deviation pairs.
  //  TODO: insert description of how to actually use these.
  inline const DeviationPairVec &deviation_pairs() const;
  
  //channel_energies(): returns a vector containing the starting (lower) energy
  //  of the gamma channels, calculated using the energy calibration
  //  coefficients as well as the deviation pairs.  These channel energies are
  //  calculated during file parsing, or any subsequent re-calibrations;  the
  //  owining MeasurmentInfo object will make an attempt so that multiple
  //  Measurments that it owns, that have the same calibration, will also return
  //  pointers from channel_energies() that point to the same spot in memory
  //  (this is primarily a memory usage optimization).
  //  Typically the vector returned by channel_energies() will have the same
  //  number of channels as gamma_counts(), however for most LowerChannelEdge
  //  calibration model files, channel_energies() may have 1 more channels (to
  //  indicate end of last channel energy).
  //  Returned pointer may be null if energy calibration is unknown/unspecified.
  inline const std::shared_ptr< const std::vector<float> > &
                                                       channel_energies() const;
  
  //gamma_counts(): the channel counts of the gamma data.
  //  Returned pointer may be null if no gamma data present, or not thie
  //  Measurment is not properly initialized.
  inline const std::shared_ptr< const std::vector<float> > &
                                                           gamma_counts() const;
  
  //neutron_counts(): the channel counts of neutron data.  Currently none of
  //  the file formats give channelized neutron data, so this function may
  //  be removed in the future; use neutron_counts_sum() instead.  Currently
  //  returned vector will have a size 1 if the file contained neutron counts.
  inline const std::vector<float> &neutron_counts() const;

  //compare_by_sample_det_time: compares by sample_number_, and then
  //  detector_number_, then by start_time_, then source_type_
  static bool compare_by_sample_det_time( const MeasurementConstShrdPtr &lhs,
                                          const MeasurementConstShrdPtr &rhs );

  //set_title(): sets the title property.
  inline void set_title( const std::string &title );
  
  //set_gamma_counts(...): XXX - should deprecate!
  //  reset real and live times, updates total gamma counts
  inline void set_gamma_counts( ShrdConstFVecPtr counts,
                                const float livetime, const float realtime );
  
  //set_neutron_counts(): XXX - should deprecate!
  //   updates total nuetron counts.  Marks containing neutrons based on
  //   if input has any entries or not.
  inline void set_neutron_counts( const std::vector<float> &counts );
  
  //set_channel_energies(...): XXX - should deprecate!
  //  if channel_energies_ or gamma_counts_ must be same number channels as
  //  before
  inline void set_channel_energies( ShrdConstFVecPtr counts );

  //popuplate_channel_energies_from_coeffs(): uses calibration_coeffs_ and 
  //  deviation_pairs_ to populate channel_energies_.
  //This function should not be used when this Measurment is part of a
  //  MeasurementInfo object (since this could waste memorry), but is intentded 
  //  for the case when this Measurment is saved all by itself to XML by 
  //  write_2006_N42_xml() and then restored using
  //  set_2006_N42_spectrum_node_info().
  //Throws if gamma_counts_ or calibration_coeffs_ is empty, or if 
  //  channel_energies_ is already populated, or if UnknownEquationType.
  //channel_energies_ is garunteed to be valid after calling this function.
  void popuplate_channel_energies_from_coeffs();

  //To set real and live times, see MeasurementInfo::set_live_time(...)
  
  //Functions that mimmic ROOTs TH1 class
  //  XXX - note that I should _really_ convert these functions to be same
  //        style as the rest of the member functions, and to use a 1-based
  //        index instead of the stupid 1 based index they are using to mimic
  //        ROOT cerns functionality.
  //  Note that in following TH1, 'bin' is a 1-based system
  //  e.g. the first binn is bin 1, NOT bin 0
  //
  //All of these functions are depreciated!  (dont use them anywhere new)
  inline float GetBinCenter( int bin ) const;  //depreciated
  inline float GetBinContent( int bin ) const; //depreciated
  inline float GetBinLowEdge( int bin ) const; //depreciated
  inline float GetBinWidth( int bin ) const;   //depreciated
  inline int GetNbinsX() const;                //depreciated
  inline float Integral( int binx1=0, int binx2 = -1 ) const; //depreciated
  inline int FindFixBin( float x ) const;      //depreciated
  inline bool CheckBinRange( int bin ) const;  //depreciated

  //I want to get rid of all the CERN ROOT inspired functions (who uses 1 based
  //  indexing?), so I am slowly re-writing them in a (more) sane manner bellow.
  
  //num_gamma_channels(): returns the minimum number of channels of either
  //  this->channel_energies_ or this->gamma_counts_; if either is not defined
  //  0 is returned.
  inline size_t num_gamma_channels() const;
  
  //find_gamma_channel(): returns gamma channel containing 'energy'.  If
  //  'energy' is bellow the zeroth channel energy, 0 is returned; if 'energy'
  //  is above last channels energy, the index for the last channel is returned.
  //  If this->channel_energies_ is not defined, an exception is thrown.
  inline size_t find_gamma_channel( const float energy ) const;
  
  //gamma_channel_content(): returns the gamma channel contents for the
  //  specified channel.  Returns 0 if this->gamma_counts_ is not defined, or
  //  if specified channel is invalid (to large).
  inline float gamma_channel_content( const size_t channel ) const;
  
  //gamma_channel_lower(): returns the lower energy of the specified gamma
  //  channel.
  //Throws exception if invalid channel number, or  !this->channel_energies_.
  inline float gamma_channel_lower( const size_t channel ) const;
  
  //gamma_channel_center(): returns the central energy of the specified
  //  channel.  For the last channel, it is assumed its width is the same as
  //  the second to last channel.
  //Throws exception if invalid channel number, or  !this->channel_energies_.
  inline float gamma_channel_center( const size_t channel ) const;

  //gamma_channel_upper(): returns the energy just past the energy range the
  //  specified channel contains (e.g. the lower edge of the next bin).  For the
  //  last bin, the bin width is assumed to be the same as the previous bin.
  //Throws exception if invalid channel number, or  !this->channel_energies_.
  inline float gamma_channel_upper( const size_t channel ) const;

  //gamma_channel_width(): returns the energy width of the specified channel.
  //  For the last channel, it is assumed its width is the same as the second
  //  to last channel.
  //Throws exception if invalid channel number, or  !this->channel_energies_.
  inline float gamma_channel_width( const size_t channel ) const;
  
  //gamma_integral(): get the integral of gamma counts between lower_energy and
  //  upper_energy; linear approximation is used for fractions of channels.
  //If this->channel_energies_ or this->gamma_counts_ is invalid, 0.0 is
  //  returned.
  double gamma_integral( float lower_energy, float upper_energy ) const;
  
  //gamma_channels_sum(): get the sum of gamma gamma channel contents for all
  //  channels, both including and inbetween, 'startbin' and 'endbin'.
  //  If 'startbin' is invalid (to large), 0.0 is returned, regardless of value
  //  of 'endbin'.
  //  If 'endbin' is to large, then it will be clamped to number of channels.
  //  'startbin' and 'endbin' will be swapped if endbin < startbin.
  //  Returns 0 if this->gamma_counts_ is invalid (doesnt throw).
  double gamma_channels_sum( size_t startbin, size_t endbin ) const;
  
  //gamma_channel_energies(): returns channel_energies_, the lower edge
  //  energies of the gamma channels).
  //  The returned shared pointer may be null.
  //  The exact same as channel_energies(), just renamed to be consistent with
  //  above accessors.
  inline const std::shared_ptr< const std::vector<float> > &
                                                 gamma_channel_energies() const;
  
  //gamma_channel_contents(): returns gamma_counts_, the gamma channel data
  //  (counts in each energy bin).
  //  The returned shared pointer may be null.
  //  The exact same as gamma_counts(), just renamed to be consistent with
  //  above accessors.
  inline const std::shared_ptr< const std::vector<float> > &
                                                 gamma_channel_contents() const;
  
  inline float gamma_energy_min() const;
  inline float gamma_energy_max() const;
  
  
  //Functions to write this Measurement object out to external formats
  
  //write_2006_N42_xml(...): similar schema as Cambio uses, but with some added
  //  tags
  bool write_2006_N42_xml( std::ostream& ostr ) const;
  
  //write_csv(...): kinda similar to Cambio, except multiple spectra in 1 files
  bool write_csv( std::ostream& ostr ) const;
  
  //write_txt(...): kinda similar to Cambio, but contians most possible info
  bool write_txt( std::ostream& ostr ) const;

  
  //Shared pointers (such as gamma_counts_) are reset in the following function,
  //  meaning that shared_ptr's outside of this object will not be effected
  //  by reset() or set(...)
  void reset();

protected:
  
  inline void set_start_time( const boost::posix_time::ptime &timestamp );
  inline void set_remarks( const std::vector<std::string> &remar );
  inline void set_source_type( const SourceType type );

  
  //Functions to set the information in this Measurement object from external
  //  sources
  //set_n42_2006_spectrum_calibration_from_id(...) added for FLIR identiFINDER N42 files,
  //  but also other N42 files which use the 'CalibrationIDs' attribute in the
  //  <Spectrum> node
  void set_n42_2006_spectrum_calibration_from_id(
                                                 const rapidxml::xml_node<char> *doc_node,
                                                 const rapidxml::xml_node<char> *spec_node );
public:
  //add_calibration_to_2012_N42_xml(): writes calibration information to the
  //  xml document, with the "id" == caliId for the <EnergyCalibration> tag.
  //  Note, this funciton should probably be protected or private, but isnt
  //  currently due to an implementation detail.
  void add_calibration_to_2012_N42_xml(
                                       ::rapidxml::xml_node<char> *RadInstrumentData,
                                       std::mutex &xmldocmutex,
                                       const int caliId ) const;

protected:
  
public:
  //set_2006_N42_spectrum_node_info(...): sets information from a 2006 N42
  //  <Spectrum> node.
  //  Note channel_energies_ and deviation pairs will note be set.  Additionally
  //  calibration information may not be set if the node did not contain it.
  void set_2006_N42_spectrum_node_info( const rapidxml::xml_node<char> *node );

protected:
  
  void set_n42_2006_count_dose_data_info( const rapidxml::xml_node<char> *dose_data,
                                std::shared_ptr<DetectorAnalysis> analysis_info,
                                std::mutex *analysis_mutex );
  
  //set_gross_count_node_info(...): throws exception on error
  //  XXX - only implements nuetron gross gounts right now
  void set_n42_2006_gross_count_node_info(
                                 const rapidxml::xml_node<char> *gross_count_measu_node );

  
  
  //combine_gamma_channels(): combines every 'nchann' channel gamma channels
  //  together.  Will throw exception if (gamma_counts_->size() % nchann) != 0.
  //  If gamma_counts_ is undefined or empty, nothing will be done.
  void combine_gamma_channels( const size_t nchann );
  
  //truncate_gamma_channels(): removes channels bellow 'keep_first_channel'
  //  and above 'keep_last_channel'.  If 'keep_under_over_flow' is true, then
  //  adds the removed gamma counts to the first/last remaining bin.
  //Throws exception if keep_first_channel>=keep_last_channel, or if
  //  keep_last_channel>=num_gamma_channels()
  void truncate_gamma_channels( const size_t keep_first_channel,
                                const size_t keep_last_channel,
                                const bool keep_under_over_flow );
  
  
  //decode_n42_2006_binning(): parses coefficients and equation type from the input
  //  xml node.
  //  If equation type is not specified, but coefficients are found, the
  //  equation type may be tried to be guessed.
  //Throws exception if an invalid node, or un-recognized calibration type.
  static void decode_n42_2006_binning( const rapidxml::xml_node<char> *calibration_node,
                             std::vector<float> &coeffs,
                             Measurement::EquationType &eqnmodel );
  
  //set_n42_2006_detector_data_node_info(): silently returns if information isnt found
  static void set_n42_2006_detector_data_node_info(
                          const rapidxml::xml_node<char> *det_data_node,
                          std::vector<MeasurementShrdPtr> &measurs_to_update );

  
  //set_info_from_txt_or_csv(...): throws upon failure
  //  XXX - currently doesnt make use of all the information written out by
  //         write_txt(...)
  //  XXX - could be made more rebust to parse more inputs
  //  XXX - currently very CPU innefiecient
  //  XXX - should be documented data formats it accepts
  void set_info_from_txt_or_csv( std::istream &istr );
  
  //set_info_from_avid_mobile_txt(): throws upon failure.
  //  Called from set_info_from_txt_or_csv().
  void set_info_from_avid_mobile_txt( std::istream &istr );
  
  //Rebin the gamma_spectrum according to new_binning. new_binning should
  //  be the lower edges of bins, in keV - this does not shift the energy of the
  //  counts (eg peaks stay at same energy/width, just might have more or less
  //  bins)
  //Throws exception if (old or new) binning or channel energies arent defined
  //  or have less than 4 or so bins.
  void rebin_by_lower_edge( ShrdConstFVecPtr new_binning );
  
  //rebin_by_eqn(...) should probably just be renamed rebin(....)
  
  //rebin_by_eqn(...): Does not shift the energy of the counts (eg
  //  peaks stay at same energy).
  //  Does not store the eqn for the case of LowerChannelEdge binning
  //  Will throw std:runtime_error(...) if gamma_counts_ is not defined.
  //  or called with UnknownEquationType
  //Throws exception if 'type' is UnknownEquationType, or gamma_counts_ is empty
  //  or not defined, or potentialy (but not necassirly) if input equation is
  //  ill-specified.
  void rebin_by_eqn( const std::vector<float> &eqn,
                    const DeviationPairVec &dev_pairs,
                    EquationType type );
  
  //If the new binning has already been calculated when rebinning by polynomial
  //  equation, you can use this to save re-caclulating the lower edeg energy
  //  values, and also possibly save memmorry by sharing the vector among many
  //  Measurement objects.
  //No checks are made that new_binning is consistent with eqn!!!
  //Throws exception if 'type' is UnknownEquationType, or gamma_counts_ is empty
  //  or not defined, or potentialy (but not necassirly) if input equation is
  //  ill-specified.
  void rebin_by_eqn( const std::vector<float> &eqn,
                    const DeviationPairVec &dev_pairs,
                    ShrdConstFVecPtr new_binning,
                    EquationType type );
  
  
  //recalibrate_by_eqn(...) passing in a valid binning
  //  ShrdConstFVecPtr object will save recomputing this quantity, however no
  //  checks are performed to make sure 'eqn' actually cooresponds to 'binning'
  //  For type==LowerChannelEdge, eqn should coorespond to the energies of the
  //  lower edges of the bin.
  //Throws exception if 'type' is UnknownEquationType, or potentially (but not
  //  necassirily) if input is ill-specified, or if the passed in binning doesnt
  //  have the proper number of channels.
  void recalibrate_by_eqn( const std::vector<float> &eqn,
                          const DeviationPairVec &dev_pairs,
                          Measurement::EquationType type,
                          ShrdConstFVecPtr binning
#if( !SpecUtils_JAVA_SWIG )
                          = ShrdConstFVecPtr()  //todo: fix this more better
#endif
                          );
  
public:
  
#if( PERFORM_DEVELOPER_CHECKS )
  //equalEnough(...): tests whether the passed in Measurement objects are
  //  equal, for most intents and purposes.  Allows some small numerical
  //  rounding to occur.
  //Throws an std::exception with a brief explanaition when an issue is found.
  static void equalEnough( const Measurement &lhs, const Measurement &rhs );
#endif
  
  
protected:
  
  //live_time_: in units of seconds.  Typially 0.0f if not specified.
  float live_time_;
  
  //real_time_: in units of seconds.  Typially 0.0f if not specified.
  float real_time_;

  //contained_neutron_: used to specify if there was a neutron detector, but
  //  0 counts were actually detected.
  bool contained_neutron_;

  //sample_number_: first sample is typically 1 (not zero like in c++)
  int sample_number_;
  
  OccupancyStatus  occupied_;
  double gamma_count_sum_;
  double neutron_counts_sum_;
  float speed_;  //in m/s
  std::string detector_name_;
  int detector_number_;
  std::string detector_type_;  //e.x. "HPGe 50%". Roughly the equivalent N42 2012 "RadDetectorDescription" node
  QualityStatus quality_status_;

  //The following objects are for the energy calibration, note that there are
  //  similar quantities for resolution and such, that arent kept track of
  SourceType     source_type_;
  EquationType   energy_calibration_model_;

  std::vector<std::string>  remarks_;
  boost::posix_time::ptime  start_time_;

  std::vector<float>        calibration_coeffs_;  //should consider making a shared pointer (for the case of LowerChannelEdge binning)
  DeviationPairVec          deviation_pairs_;     //<energy,offset>

  std::shared_ptr< const std::vector<float> >   channel_energies_;    //gamma channel_energies_[energy_channel]
  std::shared_ptr< const std::vector<float> >   gamma_counts_;        //gamma_counts_[energy_channel]
  std::vector<float>        neutron_counts_;      //neutron_counts_[energy_channel]  //XXX - I have yet to see a detector file with more than one neutron channel - maybye this variable should be eliniated and only nuterton_counts_sum_ used

  //could add Alt, Speed, and Dir here, as well as explicit valid flag
  double latitude_;  //set to -999.9 if not specified
  double longitude_; //set to -999.9 if not specified
//  double elevation_;
  
  boost::posix_time::ptime position_time_;
  
  std::string title_;  //Actually used for a number of file formats

  friend class SpecMeas;
  friend class MeasurementInfo;
  friend struct SpectrumNodeDecodeWorker;
  friend struct GrossCountNodeDecodeWorker;
};//class Measurement

/*
class CountDose
{
public:
  enum DoseType
  {
    GammaDose,
    NeutronDose,
    TotalDose,
    OtherDoseType
  };//enum DoseType


  DoseType m_doseType;
  std::string m_remark;
  boost::posix_time::ptime m_startTime;
  float m_realTime;
  float m_doseRate;
};//class CountDose
*/

//Should MeasurementInfo be renamed to something better?
class MeasurementInfo
{
public:
  //ToDo:
  //  -wcjohns needs to document how these classes are structured, now that
  //   there development has stabilized
  //     e.g. note things like how all Measurments in measurements_
  //          may or may not have the same binning, or channel_energies()
  //          may return a null pointer, or passthrough vehicles are
  //          rebinned to a consistent FullWidthFraction binning
  //  -clean code up
  //  -see XXX notes below and in source
  //  -Maybe add member function to save spectrum to various file types into
  //   a buffer (string, stream, etc...), although is this best to do versus
  //  -Allow initializing from istreams, instead of just files
  //  -Make reporting of errors consistent (eg specify which exceptions can be
  //   thrown, or if none can be thrown, specify this)
  //
  //  -The ASP 16k channels, 133 samples, 8 detectors takes up 66 MB memmorry

public:
  MeasurementInfo();
  MeasurementInfo( const MeasurementInfo &rhs );  //calls operator=
  virtual ~MeasurementInfo();
 
  //operator=(...) copies all the 'rhs' information and creates a new set of
  //  Measurment objects so that if you apply changes to *this, it will not
  //  effect the lhs; this is however fairly effieient as the Measurment
  //  objects copy shallow copy of all ShrdConstFVecPtr instances.
  //  Since it is assumed the 'rhs' is in good shape, recalc_total_counts() and
  //  cleanup_after_load() are NOT called.
  const MeasurementInfo &operator=( const MeasurementInfo &rhs );

  //load_file(...): returns true when file is successfully loaded, false
  //  otherwise. Callling this function with parser_type==kAutoParser is
  //  the easiest way to load a spectrum file when you dont know the type of
  //  file.  The file_ending_hint is only used in the case of kAutoParser
  //  and uses the file ending to effect the order of parsers tried, example
  //  values for this mught be: "n24", "pcf", "chn", etc. The entire filename
  //  can be passed in since only the letters after the last period are used
  bool load_file( const std::string &filename,
                  ParserType parser_type,
                  std::string file_ending_hint
#if( !SpecUtils_JAVA_SWIG )
                 = ""  //todo: fix this more better
#endif
                 );

  //modified(): intended to indicate if object has been modified since last save
  inline bool modified() const;
  
  //reset_modified(): intended to be called after saving object
  inline void reset_modified();
  
  //modified_since_decode(): returns if object has been modified since decoding
  inline bool modified_since_decode() const;

  //reset_modified_since_decode(): intended to be called right after any initial
  //  adjustments following openeing of an object
  inline void reset_modified_since_decode();

  //simple accessors (no thread locks are aquired)
  inline float gamma_live_time() const;
  inline float gamma_real_time() const;
  inline double gamma_count_sum() const;
  inline double neutron_counts_sum() const;
  inline const std::string &filename() const;
  inline const std::vector<std::string> &detector_names() const;
  inline const std::vector<int> &detector_numbers() const;
  inline const std::vector<std::string> &neutron_detector_names() const;
  inline const std::string &uuid() const;
  inline const std::vector<std::string> &remarks() const;
  inline int lane_number() const;
  inline const std::string &measurement_location_name() const;
  inline const std::string &inspection() const;
  inline const std::string &measurment_operator() const;
  inline const std::set<int> &sample_numbers() const;
  inline size_t num_measurements() const;
  inline DetectorType detector_type() const;
  //instrument_type(): From ICD1 Specs InstrumentType can be:
  //   PortalMonitor, SpecPortal, RadionuclideIdentifier,
  //   PersonalRadiationDetector, SurveyMeter, Spectrometer, Other
  inline const std::string &instrument_type() const;
  inline const std::string &manufacturer() const;
  inline const std::string &instrument_model() const;
  inline const std::string &instrument_id() const;
  inline std::vector< std::shared_ptr<const Measurement> > measurements() const;
  inline std::shared_ptr<const Measurement> measurement( size_t num ) const;
  inline std::shared_ptr<const DetectorAnalysis> detectors_analysis() const;
  inline bool has_gps_info() const; //mean longitude/latitude are valid gps coords
  inline double mean_latitude() const;
  inline double mean_longitude() const;

  //passthrough(): returns true if it looked like this data was from a portal
  //  or search mode data.  Not 100% right always, but pretty close.
  bool passthrough() const;
  
  //simple setters (no thread locks are aquired)
  inline void set_filename( const std::string &n );
  inline void set_remarks( const std::vector<std::string> &n );
  inline void set_uuid( const std::string &n );
  inline void set_lane_number( const int num );
  inline void set_measurement_location_name( const std::string &n );
  inline void set_inspection( const std::string &n );
  inline void set_instrument_type( const std::string &n );
  inline void set_detector_type( const DetectorType type );
  inline void set_manufacturer( const std::string &n );
  inline void set_instrument_model( const std::string &n );
  inline void set_instrument_id( const std::string &n );


  //A little more complex setters:
  //set_live_time(...) and set_real_time(...) update both the measurment
  //  you pass in, as well as *this.  If measurment is invalid, or not in
  //  measurments_, than an exception is thrown.
  void set_live_time( const float lt, MeasurementConstShrdPtr measurment );
  void set_real_time( const float rt, MeasurementConstShrdPtr measurment );

  //set_start_time(...), set_remarks(...), set_spectra_type(...) allow
  //  setting the relevant variables of the 'measurment' passed in.  The reason
  //  you have to set these variables from MeasurmentInfo class, instead of
  //  directly from the Measurement class is because you should only be dealing
  //  with const pointers to these object for both the sake of the modified_
  //  flag, but also to ensure some amount of thread safety.
  void set_start_time( const boost::posix_time::ptime &timestamp,
                       const MeasurementConstShrdPtr measurment  );
  void set_remarks( const std::vector<std::string> &remarks,
                    const MeasurementConstShrdPtr measurment  );
  void set_source_type( const Measurement::SourceType type,
                         const MeasurementConstShrdPtr measurment );
  void set_position( double longitude, double latitude,
                     boost::posix_time::ptime position_time,
                     const MeasurementConstShrdPtr measurment );
  void set_title( const std::string &title,
                  const MeasurementConstShrdPtr measurment );
  
  //set_contained_neutrons(...): sets the specified measurment as either having
  //  contained neutron counts, or not.  If specified to be false, then counts
  //  is ignored.  If true, then the neutron sum counts is set to be as
  //  specified.
  void set_contained_neutrons( const bool contained, const float counts,
                               const MeasurementConstShrdPtr measurment );

  
  //add_measurment(...): adds the measurment to this MeasurmentInfo object and
  //  if 'doCleanup' is specified, then all sums will be recalculated, and
  //  binnings made consistent.  If you do not specify 'doCleanup' then
  //  things will be roughly updated, but the more thorough cleanup_after_load()
  //  is not called.
  //Will throw if 'meas' is already in this MeasurementInfo object, otherwise
  //  will ensure it has a unique sample number/ detector number combination
  //  (which means 'meas' may be modified).  Note that if a sample number
  //  needs to be assigned, it is chosen to be the very last available sample
  //  number available if that detector does not already have that one or else
  //  its assigned to be one larger sample number - this by no means captures
  //  all use cases, but good enough for now.
  void add_measurment( MeasurementShrdPtr meas, const bool doCleanup );
  
  //remove_measurment(...): removes the measurment from this MeasurmentInfo
  //  object and if 'doCleanup' is specified, then all sums will be
  //  recalculated.  If you do not specify 'doCleanup' then make sure to call
  //  cleanup_after_load() once you are done adding/removing Measurments if a
  //  rough fix up isnt good enough.
  //Will throw if 'meas' isnt currently in this MeasurementInfo.
  void remove_measurment( MeasurementConstShrdPtr meas, const bool doCleanup );
  
  //remove_measurments(): similar to remove_measurment(...), but more efficient
  //  for removing large numbers of measurments.  This function assumes
  //  the internal state of this MeasurmentInfo object is consistent
  //  (e.g. no measurments have been added or removed without 'cleaningup').
  void remove_measurments( const std::vector<MeasurementConstShrdPtr> &meas );
  
  //combine_gamma_channels(): combines 'ncombine' gamma channels for every
  //  Measurement that has exactly nchannels.  Returns number of Measurments
  //  modified.
  //  Throws exception if( (nchannels % ncombine) != 0 )
  size_t combine_gamma_channels( const size_t ncombine, const size_t nchannels );
  
  //combine_gamma_channels(): calls equivalent function on non-const version of
  //  'm'. See Measurement::combine_gamma_channels() for comments.
  //Throws exception if 'm' isnt owned by *this.
  void combine_gamma_channels( const size_t ncombine,
                               const std::shared_ptr<const Measurement> &m );
  
  //truncate_gamma_channels(): removes all channels bellow 'keep_first_channel'
  //  and above 'keep_last_channel', for every measurment that has 'nchannels'.
  //  If keep_under_over_flow is true, then removed channel counts will be added
  //  to the first/last channel of the remaing data.
  //  Returns number of modified Measurments.
  //Throws exception if keep_last_channel>=nchannels, or if
  //  keep_first_channel>=keep_last_channel.
  size_t truncate_gamma_channels( const size_t keep_first_channel,
                                  const size_t keep_last_channel,
                                  const size_t nchannels,
                                  const bool keep_under_over_flow );
  
  //truncate_gamma_channels(): removes all channels bellow 'keep_first_channel'
  //  and above 'keep_last_channel', for specified measurment.
  //  If keep_under_over_flow is true, then removed channel counts will be added
  //  to the first/last channel of the remaing data.
  //Throws exception if invalid Measurement, or if
  //  keep_last_channel>=m->num_gamma_channels(), or if
  //  keep_first_channel>=keep_last_channel.
  void truncate_gamma_channels( const size_t keep_first_channel,
                                const size_t keep_last_channel,
                                const bool keep_under_over_flow,
                                const std::shared_ptr<const Measurement> &m );
  
  
  //Not so simple accessors
  //occupancy_number_from_remarks(): tries to find the occupancy number from
  //  remarks_.
  //  returns -1 if one cannot be found.
  //  Currently only succeds if a remark starts with the string
  //    "Occupancy number = ", where this should be improved in the future
  //  In Principle there should be an EventNumber under the Event
  //  tag in ICD1 files, but I dont havea any exammples of this, so implementing
  //  EventNumber will have to wait
  int occupancy_number_from_remarks() const;

  //get all measurements cooresponding to 'sample_number', where
  //  sample_number may not be a zero based index (eg may start at one)
  std::vector< std::shared_ptr<const Measurement> > sample_measurements(
                                                const int sample_number ) const;

  std::shared_ptr<const Measurement> measurement( const int sample_number,
                                           const std::string &det_name ) const;
  std::shared_ptr<const Measurement> measurement( const int sample_number,
                                           const int detector_number ) const;

  //suggested_gamma_binning_index(...): returns the index of measurments_ to use
  //  as the binning, when you are summing over the specified sample numbers and
  //  detectors.
  //  This function chooses the Measurment with the largest number of gamma
  //  channels, if the Measurements have varying number of gamma channels.
  //'det_to_use' must be same size as, and coorespond 1:1 with detector_numbers_
  //Throws exception if 'det_to_use' is wrong size, no measurments available, or
  //  other errors.
  size_t suggested_gamma_binning_index( const std::set<int> &sample_numbers,
                                    const std::vector<bool> &det_to_use ) const;
  
  //sum_measurements(...): sums the specified sample_numbers and det_to_use
  //  (logically and-ed).
  //'det_to_use' must be same size as, and coorespond 1:1 with detector_numbers_
  //  or else an exception may be thrown.
  //Requires the selected samples and detectors to have at least one spectrum
  //  that can serve as gamma binning, or else nullptr will be returned.
  MeasurementShrdPtr sum_measurements( const std::set<int> &sample_numbers,
                                    const std::vector<bool> &det_to_use ) const;
  
  //sum_measurements(...): a convienience function for calling the the other
  //  form of this function.  'det_numbers' should contain the numbers of
  //  the detectors you would like included in the sum.  If any invalid detector
  //  numbers are included, an exception will be thrown.
  //Requires the selected samples and detectors to have at least one spectrum
  //  that can serve as gamma binning, or else nullptr will be returned.
  MeasurementShrdPtr sum_measurements( const std::set<int> &sample_numbers,
                                    const std::vector<int> &det_numbers ) const;
  
  //sum_measurements(...): sums measurments similar to the other variants by
  //  the same name, but uses the 'binTo' Measurment passed in as the bassis
  //  for the energy calibration binning.
  //'det_to_use' must be same size as, and coorespond 1:1 with detector_numbers_
  //  or else an exception may be thrown.
  //Requires the selected samples and detectors to have at least one spectrum
  //  that can serve as gamma binning, or else nullptr will be returned.
  MeasurementShrdPtr sum_measurements( const std::set<int> &sample_numbers,
                                       const std::vector<bool> &det_to_use,
                                       const MeasurementConstShrdPtr binTo ) const;
  
  MeasurementShrdPtr sum_measurements( const std::set<int> &sample_numbers,
                                      const std::vector<int> &det_numbers,
                                      const MeasurementConstShrdPtr binTo ) const;
  
  
  
  //memmorysize(): should be reasonbly accurate, but could definetly be off by a
  //  little bit.  Tries to take into account the shared float vectors may be
  //  shared between Measurment objects.
  size_t memmorysize() const; //in bytes

  //gamma_channel_counts(): loops over the Measurments and returns a set<size_t>
  //  containing all the channel_energies_->size() results
  std::set<size_t> gamma_channel_counts() const;

  //num_gamma_channels(): loops over the Measurments, and returns the size of
  //  the first channel_energies_ encountered.
  size_t num_gamma_channels() const;

  //keep_n_bin_spectra_only(..): return number of removed spectra
  //  doesnt remove neutron detectors
  size_t keep_n_bin_spectra_only( size_t nbin );

  
  //energy_cal_variants(): Some N42 files may have the same spectra, binned
  //  differently multiple times (ex. "Lin" and "Sqrt"), and when this is
  //  detected the detector they are assigned to is appended with a prefix
  //  such as "_intercal_LIN"
  std::set<std::string> energy_cal_variants() const;
  
  
  //keep_energy_cal_variant(): When #energy_cal_variants() returns multiple
  //  variants, you can use this function to remove all energy calibration
  //  variants, besides the one you specify, from the measurment.  If a spectrum
  //  is not part of a variant, it is kept.
  //Returns return number of removed spectra.
  //Throws exception if you pass in an invalid variant.
  size_t keep_energy_cal_variant( const std::string variant );
  
  
  //rremove_neutron_measurments() only removes neutron measurments that do not
  //  have a gamma binning defined
  size_t remove_neutron_measurments();

  //background_sample_number() returns numeric_limits<int>::min() if no
  // background is found; behavior undefined for more than one background sample
  int background_sample_number() const;
  
  //reset(): resets all variables to same state as just after construction
  void reset();

  //The following do not throw, but return false on falure, as well as call
  //  reset().  None of these functions are tested very well.
  virtual bool load_N42_file( const std::string &filename );
  bool load_pcf_file( const std::string &filename );
  bool load_spc_file( const std::string &filename );
  bool load_chn_file( const std::string &filename );
  bool load_iaea_file( const std::string &filename );
  bool load_binary_exploranium_file( const std::string &file_name );
  bool load_micro_raider_file( const std::string &filename );
  bool load_txt_or_csv_file( const std::string &filename );
  bool load_camberra_cnf_file( const std::string &filename );
  bool load_tracs_mps_file( const std::string &filename );
  bool load_aram_file( const std::string &filename );
  
  bool load_spectroscopic_daily_file( const std::string &filename );
  bool load_amptek_file( const std::string &filename );
  bool load_ortec_listmode_file( const std::string &filename );

  
  //load_from_N42: loads spectrum from a stream.  If failure, will return false
  //  and set the stream position back to original position.
  virtual bool load_from_N42( std::istream &istr );
 
  //load_N42_from_data(...): raw xml file data - must be 0 terminated
  virtual bool load_N42_from_data( char *data );
  
  //load_from_iaea_spc(...) and load_from_binary_spc(...) reset the
  //  input stream to original position and sets *this to reset state upon
  //  failure
  //  20120910: not well tested yet - only implemented for files with a single
  //  spectrum (I didnt have any multi-spectrum files to view structure)
  bool load_from_iaea_spc( std::istream &input );
  bool load_from_binary_spc( std::istream &input );
  
  //load_from_N42_document(...): loads information from the N42 document, either
  //  2006 or 2012 variants.
  //  May throw.  Returns success status.
  bool load_from_N42_document( const rapidxml::xml_node<char> *document_node );
  
  //load_from_micro_raider_from_data(...): loads information from the Micro
  //  Raider XML format.  This function takes as input the files contents.
  bool load_from_micro_raider_from_data( const char *data );
  
  //load_from_binary_exploranium(...): returns success status; calls reset()
  //  upon falure, and puts istr to original location
  bool load_from_binary_exploranium( std::istream &istr );
  
  //load_from_pcf(...): Set info from GADRAS PCF files.  Returns success
  //  status.  Calls reset() upon failure, and puts istr to original location.
  bool load_from_pcf( std::istream &istr );

  //load_from_txt_or_csv(...): returns state of successfulness, and typically
  //  wont throw.  Resets the object/stream upon failure.
  //  XXX - currently very CPU innefiecient
  //  XXX - should be documented data formats it accepts
  //  XXX - should add some protections to make sure file is not a binary file
  //        (currently this is done in loadTxtOrCsvFile(...) )
  bool load_from_txt_or_csv( std::istream &istr );

  //load_from_Gr135_txt(...): called from load_from_txt_or_csv(...) if it
  //  thinks it could be a GR135 text file
  bool load_from_Gr135_txt( std::istream &istr );
  
  //setInfoFromSpectroscopicDailyFiles(...): file format used for SLD
  //  spectroscopic portals, to reduce bandwidth; a txt based format
  bool load_from_spectroscopic_daily_file( std::istream &input );
  
  //setInfoFromAmetekMcaFiles(...): asncii file format used by Ametek MCA
  //  devices;
  bool load_from_amptek_mca( std::istream &input );
  
  //load_from_ortec_listmode(...): listmode data from ORTEC digiBASE (digibase-E
  //  and PRO list format not supported yet).
  bool load_from_ortec_listmode( std::istream &input );
  
  //bool load_from_iaea(...): an ASCII format standardized by the IAEA; not all
  //  portions of the standard have been implemented, since they are either
  //  not applicable, or I havent seen an example use of them. Also, only
  //  tested to work with multiple spectra in one file from PeakEasy (didnt
  //  seem like standard explain how it should be done, but I didnt look very
  //  closely).
  //see
  //  http://www.gbs-elektronik.de/old_page/mca/appendi1.pdf
  //  and
  //  http://www.gbs-elektronik.de/fileadmin/download/datasheets/spe-file-format_datasheet.pdf
  //  for specifications of IAEA file standards.
  //This function is computationally slower than it could be, to allow for
  //  a little more diverse set of intput.
  bool load_from_iaea( std::istream &istr );

  //bool load_from_chn(...): Load information from ORTECs binary CHN file.
  bool load_from_chn( std::istream &input );
  
  //bool load_from_camberra_cnf(...): loads info from Camberra CNF files.  Not
  //  particularly well tested, but seems to work okay
  bool load_from_camberra_cnf( std::istream &input );
  bool load_from_tracs_mps( std::istream &input );
  
  bool load_from_aram( std::istream &input );
  
  
  //cleanup_after_load():  Fixes up inconsistent calibrations, binnings and such,
  //  May throw exception on error.
  enum CleanupAfterLoadFlags
  {
    //RebinToCommonBinning: ensures all spectrums in measurements_ have the same
    //  binning.  Will cause the kHasCommonBinning flag to be set in
    //  properties_flags_ (this flag may get set anyway if all same spectrum in
    //  the file).
    RebinToCommonBinning = 0x1,
    
    //DontChangeOrReorderSamples: makes it so cleanup_after_load() wont change
    //  currently assigned sample or detector numbers, or change order of
    //  measurements_.
    DontChangeOrReorderSamples = 0x2,
    
    //Before 20141110 the rebin option was not available (always done), but for
    //  a few non-InterSpec applications, this isnt always desirable, so as a
    //  hack, we will define StandardCleanup, that you can change depending on
    //  the application; currently in the code cleanup_after_load() is always
    //  called without an argument, meaning it will default to StandardCleanup.
    //This is a bit of a hack, but time is limited...
#if( REBIN_FILES_TO_SINGLE_BINNING )
    StandardCleanup = RebinToCommonBinning
#else
    StandardCleanup = 0x0
#endif
  };//enum CleanupFlags
  
  virtual void cleanup_after_load( const unsigned int flags = StandardCleanup );
  
  //recalc_total_counts(): does NOT recalculate totals for Measurements, just
  //  from them!  Is always called by cleanup_after_load().
  void recalc_total_counts();
  
  //Rebin the gamma_spectrum according to new_binning. new_binning should
  //  be the lower edges of bins, in keV - note: alters individual bin contents.
  //  Necassarily information losing, so use sparingly and dont compound calls.
  void rebin_by_eqn( const std::vector<float> &eqn,
                     const DeviationPairVec &dev_pairs,
                     Measurement::EquationType type );

  //Recalibrate the spectrum so the existing channel counts coorespond
  //  to the energies of the new binning - note: does not alter bin contents.
  //  Also not that the recalibration is applied to all gamma detectors
  void recalibrate_by_lower_edge( ShrdConstFVecPtr binning );
  void recalibrate_by_eqn( const std::vector<float> &eqn,
                           const DeviationPairVec &dev_pairs,
                           Measurement::EquationType type );
  
  //If only certain detectors are specified, then those detectors will be
  //  recalibrated (channel contents not changed).  If rebin_other_detectors
  //  is true, then the other detectors will be rebinned (channel contents
  //  changed) so that they have the same calibration as the rest of the
  //  detectors.
  //Will throw exception if an empty set of detectors is passed in, or if none
  //  of the passed in names match any of the available names, since these are
  //  both likely a mistakes
  void recalibrate_by_eqn( const std::vector<float> &eqn,
                           const DeviationPairVec &dev_pairs,
                           Measurement::EquationType type,
                           const std::vector<std::string> &detectors,
                           const bool rebin_other_detectors );

  //calibration_is_valid(...): checks to make sure passed in calibration is
  //  valid.  Polynomial and FullRangeFraction types are checked to make sure
  //  that energy of the first two and last two bins are increasing left to
  //  right. LowerChannelEdge type is check that each bin is increasing over
  //  previous, and that it has at least as many bins as nbin.
  //  UnknownEquationType always returns false.
  static bool calibration_is_valid( const Measurement::EquationType type,
                                  const std::vector<float> &eqn,
                          const std::vector< std::pair<float,float> > &devpairs,
                                  size_t nbin );
  
  //Functions to export to various file formats

  //write_to_file(...): Writes the contents of this object to the specified
  //  file in the specified format.  For file formats such as CHN and SPC files
  //  that can only have one record, all Measurement's are summed togoether
  //  and written to the file; for other formats that can have multiple records,
  //  all Measurement will be written to the output file.
  //  Throws exception if an error is encountered.  Will not overwrite  an
  //  existing file.  Assumes no other programs/threads will attempt to access
  //  the same file while this function is being called.
  //  If no exception is thrown the specified file will exist and contain the
  //  relevant information/format.
  //  If an exception is thrown, there are no garuntees as to if the file will
  //  exist, or what its contents will be.
  void write_to_file( const std::string filename,
                      const SaveSpectrumAsType format ) const;
  
  //write_to_file(...): Similar to above write_to_file(...), except only the
  //  specified sample and detector numbers will be written.  When writing to
  //  the output file, if the ouput format supports multiple records, then
  //  each Measurement will be written to its own record.  If the format
  //  only supports a single record, then the specified sample and detector
  //  numbers will be summed together.
  //  Throws exception on error.  Will not overwrite existing file.
  void write_to_file( const std::string filename,
                      const std::set<int> sample_nums,
                      const std::set<int> det_nums,
                      const SaveSpectrumAsType format ) const;

  //write_to_file(...): Similar to above write_to_file(...), except only the
  //  specified sample and detector numbers will be written.  When writing to
  //  the output file, if the ouput format supports multiple records, then
  //  each Measurement will be written to its own record.  If the format
  //  only supports a single record, then the specified sample and detector
  //  numbers will be summed together.
  //  Throws exception on error.  Will not overwrite existing file.
  void write_to_file( const std::string filename,
                      const std::vector<int> sample_nums,
                      const std::vector<int> det_nums,
                      const SaveSpectrumAsType format ) const;
  
  //write(...): Wites the specified sample and detector numbers to the provided
  //  stream, in the format specified.  If the output format allows multiple
  //  records, each Measurement will be placed in its own record.  If the output
  //  format only allows a single records, the specified sample and detector
  //  numbers will be summed.
  //
  //  Throws exception on error.
  void write( std::ostream &strm,
              std::set<int> sample_nums,
              const std::set<int> det_nums,
              const SaveSpectrumAsType format ) const;
  
  //write_pcf(...): writes to GADRAS format, using the convention of what looks
  //  to be newer GADRAS files that include "DeviationPairsInFile" information.
  //  Tries to preserve speed, sample number, foreground/background information
  //  into the title if there is room
  bool write_pcf( std::ostream &ostr ) const;
  
  //write_2006_N42(...): writes 2006 N42 format, similar to that of Cambio (for
  //  compatibility), but tries to includes other information like GPS
  //  coordinates, speed, occupied, sample number, etc.
  //  XXX - deviation pairs are not tested, and probably not up to standard.
  virtual bool write_2006_N42( std::ostream &ostr ) const;
  
  
  /** The spectra in the current file are written out in a two column
   format (seperated by a comma); the first column is gamma channel
   lower edge energy, the second column is channel counts.  Each
   spectrum in the file are written out contiguously and seperated
   by a header that reads \"Energy, Data\".  Windows style line
   endings are used (\\n\\r).
   This format loses all non-spectral information, including live
   and real times, and is intended to be an easy way to import the
   spectral information into other programs like Excel.
   */
  bool write_csv( std::ostream &ostr ) const;
  
  
  /** Spectrum(s) will be written to an ascii text file.  At the
   beggining of the file the original file name, total live and
   real times, sum gamma counts, sum neutron counts, and any file
   level remarks will be written on seperate labeled lines.
   Then after two blank lines each spectrum in the current file
   will be written, seperated by two blank lines.
   Each spectrum will contain all remarks, measurment start time
   (if valid), live and real times, sample number, detector name,
   detector type, GPS coordinates/time (if valid), serial number
   (if present), energy
   calibration type and coefficient values, and neutron counts (if
   valid); the
   channel number, channel lower energy, and channel counts is then
   provided with each channel being placed on a seperate line and
   each field being seperated by a space.
   Any detector provided analysis in the original program, as well
   manufacturer, UUID, deviation pairs, lane information,
   location name, or spectrum title is lost.
   Cambio or other programs may not be able to read back in all
   information written to the txt file.
   The Windows line ending convention is used (\\n\\r).
   This is not a standard format commonly read by other programs,
   and is intended as a easily human readable summary of the
   spectrum file information
   */
  bool write_txt( std::ostream &ostr ) const;
  
  //write_integer_chn(): Sums over the passed in sample_nums and det_nums to
  //  write a single spectrum with integer count bins in CHN format to the
  //  stream.  If sample_nums and/or det_nums is empty, then all sample and/or
  //  detector numbers are assumed to be wanted.  Values in det_nums coorespond
  //  to values in detector_numbers_.
  // This format preserves the gamma spectrum, measurment start time, spectrum
  //  title (up to 63 characters)," detector description, and energy
  //  calibration.
  //  Energy deviation pairs and neutron counts, as well as any other meta
  //  information is not preserved.
  //  Also, wcjohns reverse engineered the file format to begin with, so this
  //  function may not be writing all the information it could to the file.
  bool write_integer_chn( std::ostream &ostr, std::set<int> sample_nums,
                          const std::set<int> &det_nums ) const;
  
  /** Enum to specify the type of binary SPC file to write. */
  enum SpcBinaryType{ IntegerSpcType, FloatSpcType };
  
  /** This format allows a single spectrum in the output, so the sample and 
   detector numbers specified will be summed to produce a single spectrum.
   If sample_nums and/or det_nums is empty, then all sample and/or detector
   numbers are assumed to be wanted.
   This format preserves the gamma spectrum, neutron counts, gps info, 
   measurment start time, detector serial number, and energy calibration (if 
   polynomnial or FWHM).
   Energy deviation pairs, analysis results, and other meta information will be
   lost.
   */
  bool write_binary_spc( std::ostream &ostr,
                         const SpcBinaryType type,
                         std::set<int> sample_nums,
                         const std::set<int> &det_nums ) const;
  
  
  /** This format allows a single spectrum in the output, so the sample and 
      detector numbers specified will be summed to produce a single spectrum.
      If sample_nums and/or det_nums is empty, then all sample and/or detector
      numbers are assumed to be wanted.
      This format preserves the gamma spectrum, neutron counts, gps info,
      measurment start time, detector serial number, and energy calibration (if
      polynomnial or FWHM).
      Energy deviation pairs, some analysis analysis results, and possibly some, 
      but not all, meta information will be lost.
   */
  bool write_ascii_spc( std::ostream &output,
                       std::set<int> sample_nums,
                       const std::set<int> &det_nums ) const;
  
  /** This format allows multiple spectra to be written to the file, however
      information on detector will be lost (records ordered by detector number,
      then by sample number), and spectra will be rebinned to 256 channels, and
      energy calibration information will be lost.  All analysis results,
      meta-information, and neutron counts will be lost.  Channel counts are 
      stored as 16 bit unsigned ints; counts over 32,767 will be wrapped 
      moduloed.
      Currently, if spectra are not a multiple of 256 channels, they will be
      linearized to 256 channels; note that the full range of the spectrum will
      be preserved, although it looks like GR130 spectra should only be ranged
      from 0 to 1.5, or 1 to 3.0 MeV.
   */
  bool write_binary_exploranium_gr130v0( std::ostream &output ) const;
  
  /** This format allows multiple spectra to be written to the file, however
   information on detector will be lost (records ordered by detector number,
   then by sample number), and spectra will be rebinned to 1024 channels, and 
   neutron counts preserved (marked zero if original didnt have any counts).
   Energy calibration will be converted to third order polynomial if polynomial 
   or full width fraction; spectra with energy calibration of by lower energy 
   will be converted to second order polynomial.
   All analysis results and meta-information will be lost.
   Channel counts over 65,535 will be trucated to this value.
   Currently, if spectra are not a multiple of 1024 channels, they will be 
   linearized to 1024 channels.
   */
  bool write_binary_exploranium_gr135v2( std::ostream &output ) const;
  
  
  virtual bool write_iaea_spe( std::ostream &output,
                               std::set<int> sample_nums,
                               const std::set<int> &det_nums ) const;

#if( ENABLE_D3_CHART_EXPORTING )
  bool write_d3_html( std::ostream &output,
                      const D3SpectrumExport::D3SpectrumChartOptions &options,
                      std::set<int> sample_nums,
                      const std::set<int> &det_nums ) const;
#endif
  
  //Incase InterSpec specific changes are made, please change this number
  //  Version 4: Made it so portal data that starts with a long background as
  //             its first sample will have the 'id' attribute of the
  //             <RadMeasurement> be "Background", and subsequent samples have
  //             ids similar to "Survey 1", "Survey 2", etc.  The id attrib
  //             of <spectrum> tags also start with these same strings.
  //             (Hack to be compatible with GADRAS)
  #define MeasurementInfo_2012N42_VERSION 4
  
  //The goal of create_2012_N42_xml(...) is that when read back into
  //  load_2012_N42_from_doc(...), the results should be identical (i.e. can
  //  make round trip).  May return null pointer if something goes drasitcally
  //  wrong.
  virtual std::shared_ptr< ::rapidxml::xml_document<char> > create_2012_N42_xml() const;
  
  //write_2012_N42(): a convience function that uses create_2012_N42_xml(...) to
  //  do the actual work
  virtual bool write_2012_N42( std::ostream& ostr ) const;


#if( PERFORM_DEVELOPER_CHECKS )
  //equalEnough(...): tests whether the passed in MeasurementInfo objects are
  //  equal, for most intents and purposes.  Allows some small numerical
  //  rounding to occur.
  //Throws an std::exception with a brief explanaition when an issue is found.
  static void equalEnough( const MeasurementInfo &lhs, const MeasurementInfo &rhs );
  
  double deep_gamma_count_sum() const;
  double deep_neutron_count_sum() const;
#endif
  
protected:
  
  //measurment(...): converts a const Measurement ptr to a non-const Measurement
  // ptr, as well as checking that the Measurement actually belong to this
  //  MeasurementInfo object. Returns empty pointer on error.
  //  Does not obtain a thread lock.
  std::shared_ptr<Measurement> measurment( std::shared_ptr<const Measurement> meas );
  
  //generate_psuedo_uuid(): uses things like gamma and nuetrons sums, times,
  //  number of samples, and calibration to generate a psuedo-UUID.
  //  Is called from cleanup_after_load() if a UUID doesnt already exist.
  //  Results kinda conform to the expected UUID v4 format (a random
  //  UUID) of YYMMDDHH-MMSS-4FFx-axxx-xxxxxxxxxxxx, where {Y,M,D,H,M,S,F}
  //  relate to the time of first measurment, and the x's are based off of
  //  hashing the various properties of the spectrum.
  //  Does not obtain a thread lock.
  std::string generate_psuedo_uuid() const;
  
  //find_detector_names(): looks through measurements_ to find all detector
  //  names.
  //  Does not obtain a thread lock.
  std::set<std::string> find_detector_names() const;
  
  //set_detector_type_from_other_info(): sets detector_type_ using
  //  manufacturer_, instrument_model_ instrument_id_, or other info, only if
  //  detector_type_ isnt set.  Called from cleanup_after_load()
  void set_detector_type_from_other_info();
  
  //set_n42_2006_instrument_info_node_info(...): called from load_2006_N42_from_doc(...)
  //  Does not obtain a thread lock.
  void set_n42_2006_instrument_info_node_info( const rapidxml::xml_node<char> *info_node );
  
  //set_n42_2006_deviation_pair_info(...): called from load_2006_N42_from_doc(...)
  //  Does not obtain a thread lock.
  void set_n42_2006_deviation_pair_info( const rapidxml::xml_node<char> *info_node,
                            std::vector<MeasurementShrdPtr> &measurs_to_update );
  
  //ensure_unique_sample_numbers(): unsures uniqu detector-name samplenumber
  //  combos
  //  -Not stable in terms of garunteeing measurments started at same time
  //   will have same sample number
  //  -ensures measurments_ sorted according to Measurement::compare_by_sample_det_time
  //  -if UUID is empty will generate a deterministic pseudo-UUID
  //Called from cleanup_after_load()
  //TODO: function should be parralized for measurments with many samples
  //      - currently measurments with large numbers of measurments (>500)
  //        dont ensure dense sample numbers as a computationally fast
  //        workaround.
  void ensure_unique_sample_numbers();
  
  //has_unique_sample_and_detector_numbers(): checks to make sure 
  bool has_unique_sample_and_detector_numbers() const;
  
  //setSampleNumbersByTimeStamp():
  //For files with < 500 samples, doesn't guarantee unique detctor-name
  //  sample-number combinations, but does try to preserve initial sample-number
  //  assignments.
  //For files with >= 500 samples, it does garuntee unique detector-name
  //  sample-number combinations, but it does not preserve initial sample-number
  //  assignments.
  //Called from cleanup_after_load()
  void set_sample_numbers_by_time_stamp();
  
  
  //load20XXN42File(...): loads the specialized type on N42 file.
  //  Throws std::exception on error loading.
  //  Both functions assume they are being called from loadN42FileData(...).
  void load_2006_N42_from_doc( const rapidxml::xml_node<char> *document_node );
  
  //ICD1 2012 info at http://www.nist.gov/pml/div682/grp04/upload/n42.xsd
  //  Example files at: http://www.nist.gov/pml/div682/grp04/n42-2012.cfm
  //  2012 ICD1 support still not very tested, and doesnt support non-ASP portals
  //  Dose rates and gross sums not parsed
  //  GPS and detector statuses are also not supported
  void load_2012_N42_from_doc( const rapidxml::xml_node<char> *document_node );

  
  static void add_analysis_results_to_2012_N42( const DetectorAnalysis &ana,
                                               ::rapidxml::xml_node<char> *RadInstrumentData,
                                               std::mutex &xmldocmutex );

  
  //add_spectra_to_measurment_node_in_2012_N42_xml(...): Adds the given
  //  spectra to the specified RadMeasurementNode.  All measurments should
  //  have the sample sample number, and the entries in calibid should
  //  coorespond one to one to the entries in measurments.
  //  If something drastically goes wrong, and an exception is thrown somewhere
  //  this function will not throw, it will print an error to stderror and not
  //  insert itself into the DOM; this is so this function is safe to call in
  //  its own thread with no error handling.  I expect this to never happen, so
  //  I'm not bothing with any better error handling.
  static void add_spectra_to_measurment_node_in_2012_N42_xml(
          ::rapidxml::xml_node<char> *RadMeasurementNode,
          const std::vector< std::shared_ptr<const Measurement> > measurments,
          const std::vector<size_t> calibids,
          std::mutex &xmldocmutex );
  
  
  //2012 N42 helper functions for loading (may throw exceptions)
  void set_2012_N42_instrument_info( const rapidxml::xml_node<char> *inst_info_node );
  static std::string concat_2012_N42_characteristic_node( const rapidxml::xml_node<char> *char_node );
  
  //decode_2012_N42_rad_measurment_node: a function to help decode 2012 N42
  //  RadMeasurment nodes in a mutlithreaded fashion.  To avoid having to put
  //  some class definitions in this header, I am using boost::any, which is
  //  horrble and lame, but in the end is typesafe.  Also, this helper function
  //  has to be a member function in order to access the member variables.
  //  I would preffer you didnt awknowledge the existence of this function.
  //id_to_dettypeany and calibrations_any can be empty, but you wont get the
  //  information for calibration
  static void decode_2012_N42_rad_measurment_node(
                                std::vector< MeasurementShrdPtr > &measurments,
                                const rapidxml::xml_node<char> *meas_node,
                                const IdToDetectorType *id_to_dettypeany_ptr,
                                DetectorToCalibInfo *calibrations_ptr,
                                std::mutex &meas_mutex,
                                std::mutex &calib_mutex );

  //decode_2012_N42_detector_state_and_quality: gets GPS and detector quality
  //  status as well as InterSpec specific RadMeasurementExtension infor (title)
  static void decode_2012_N42_detector_state_and_quality( MeasurementShrdPtr meas,
                                   const rapidxml::xml_node<char> *meas_node );
  
  //Gets N42 2012 <RadDetectorKindCode> element value
  std::string determine_rad_detector_kind_code() const;
  
  //setMeasurmentLocationInformation(...):  sets the measurment information
  //  for a particular <Measurement> section of N42 data.  The parced data
  //  sets both MeasurmentInfo member variables, as well as member variables
  //  (particularly gps) of the Measurnment's passed in (that should belong to
  //  the same <Measurement> section of the N42 file, since there may be
  //  multiple spectrums per measurment).
  void set_n42_2006_measurment_location_information(
                    const rapidxml::xml_node<char> *measured_item_info_node,
                    std::vector<MeasurementShrdPtr> measurements_applicable );

  //do_channel_data_xform(): utility function called by
  //  truncate_gamma_channels() and combine_gamma_channels().  For each
  //  Measurment with 'nchannels', xform is called for it.  Returns
  //  number of modified Measurments.  Will appropriately modify the
  //  kHasCommonBinning and kAllSpectraSameNumberChannels bits of
  //  properties_flags_, as well as set modified_ and modifiedSinceDecode_.
  size_t do_channel_data_xform( const size_t nchannels,
                std::function< void(std::shared_ptr<Measurement>) > xform );
  
  //Data members
  float gamma_live_time_;      //sum over all measurments
  float gamma_real_time_;      //sum over all measurments
  double gamma_count_sum_;      //sum over all measurments
  double neutron_counts_sum_;   //sum over all measurments
  std::string                 filename_;
  std::vector<std::string>    detector_names_;          //Names may have "_intercal_..." appended to them to account for multiple binnings of the same data.
  std::vector<int>            detector_numbers_;        //in same order as detector_names_
  std::vector<std::string>    neutron_detector_names_;  //These are the names of detectors that may hold nuetron information

  std::string uuid_;
  std::vector<std::string> remarks_;

  //These go under the <MeasuredItemInformation> node in N42 file
  //  There are more fields were not currently checking for (also remarks under
  //  this node are just put into remarks_), and not necassarily set for other
  //  file formats.
  //  Also, if multiple locations are specified in the file, the last one
  //  of the files what gets put here.
  //In the future all this info should be placed in a shared_ptr, and only
  //  pouplated if it actually exists in the file
  //  Should also consider moving to the Measurment class
  int lane_number_;
  std::string measurement_location_name_;
  std::string inspection_;
  std::string measurment_operator_;


  //Start dealing with sample numbers
  std::set<int> sample_numbers_;

  // map from sample_number to a vector with indices of measurments_ of all
  //  Measurement with that sample_number
  std::map<int, std::vector<size_t> > sample_to_measurments_;


  DetectorType detector_type_;  //This is deduced from the file

  //instrument_type_:
  //  From 2006 ICD1 specs (under node InstrumentType), can be:
  //    PortalMonitor, SpecPortal, RadionuclideIdentifier,
  //    PersonalRadiationDetector, SurveyMeter, Spectrometer, Other
  //  From 2012 ICD1 specs (under node RadInstrumentClassCode) can be:
  //    "Backpack", "Dosimeter",
  //    "Electronic Personal Emergency Radiation Detector", "Mobile System",
  //    "Network Area Monitor", "Neutron Handheld", "Personal Radiation Detector",
  //    "Radionuclide Identifier", "Portal Monitor",
  //    "Spectroscopic Portal Monitor",
  //    "Spectroscopic Personal Radiation Detector", "Gamma Handheld",
  //    "Transportable System", "Other"
  //I am currently using "RadionuclideIdentifier" for handheld detectors read in
  //  from not n42 files
  std::string instrument_type_;
  
  std::string manufacturer_;
  std::string instrument_model_;
  
  std::string instrument_id_;     //often times the serial number
  
  /** Equivalent of N42 2012 RadInstrumentVersion tag that has a name (first)
      and verson (second) field
   */
  std::vector<std::pair<std::string,std::string> > component_versions_;
  
  std::vector< std::shared_ptr<Measurement> > measurements_;

  double mean_latitude_, mean_longitude_;
  
  
  //This should actually be a map<> or something so there can be multiple
  // DetectorAnalysis objects for a file, each indexed by
  // radMeasurementGroupReferences, or radMeasurementReferences
  std::shared_ptr<const DetectorAnalysis> detectors_analysis_;

  
  //properties_flags_: intenteded to indicate boolean things about the
  //  measurment style, origin, properties, or other values.
  //  These flags are calculated and set in the cleanup_after_load() function.
  //  This value is also not included in the computation of the hash of this
  //  object in generate_psuedo_uuid().
  enum MeasurmentPorperties
  {
    //kPassthroughOrSearchMode: gretaer than 5 samples, with average real time
    //  less than 2.5 seconds.  May be improved in the future to ensure
    //  measurments are sequential.
    kPassthroughOrSearchMode = (1 << 0),
    
    //kHasCommonBinning: ensures that all spectrums in measurements_ share the
    //  same binning.
    kHasCommonBinning = (1 << 1),
    
    //kRebinnedToCommonBinning: spectra had to be rebinned in order to make it
    //  so all spectrums in measurements_ share the same binning.
    kRebinnedToCommonBinning = (1 << 2),
    
    //kAllSpectraSameNumberChannels: lets you know that all spectra have the
    //  same number of channels.
    kAllSpectraSameNumberChannels = (1 << 3),
    
    //kNotTimeSortedOrder: lets you know measurements_ is not sorted by start
    //  time.  Only possibly marked when cleanup_after_load() is called with the
    //  DontChangeOrReorderSamples flag.
    kNotTimeSortedOrder = (1 << 4),
    
    //kNotSampleDetectorTimeSorted: if set, then indicates measurements_
    //  is not sorted by sample number, then detector number, then time.
    kNotSampleDetectorTimeSorted = (1 << 5),
    
    //kNotUniqueSampleDetectorNumbers: marked when the the combination of
    //  sample and detector numbers does not uniquly identify a Measurment.
    //  May be marked when cleanup_after_load() is called with the
    //  DontChangeOrReorderSamples flag.  If not set, then each measurment for
    //  a sample number has a the same start time.
    kNotUniqueSampleDetectorNumbers = (1 << 6)
  };//enum MeasurmentPorperties
  
  uint32_t properties_flags_;
  
  //modified_: intended to be used to determine if changes have been
  //  made to oject since the last time the spectrum was saved to disk.
  //  Set to true by (hopefully) all MeasurementInfo functions which are
  //  modifying, except cleanup_after_load() and reset() which sets it to false.
  //  Is _not_ set to false by any of the "saving" functions, this is your
  //  responsibility.
  bool modified_;
  
  //modifiedSinceDecode_: indicates if file has been modified at all since
  // cleanup_after_load() was called.  Set to true by (hopefully) all
  //  MeasurementInfo functions which are modifying, except cleanup_after_load()
  //  and reset() which sets it to false.
  bool modifiedSinceDecode_;
 
protected:
  mutable std::recursive_mutex mutex_;
public:
  std::recursive_mutex &mutex() const { return mutex_; };
};//class MeasurementInfo


//DetectorAnalysisResult and DetectorAnalysis are a first hack at recording
//  analysis information from the detectors file in order to use it later.
//  I'm not happy with the current design, and this information is not
//  retrieved from all (most) file types...  But I needed this in, and was
//  under a serious time crunch.
//
//Should seperate DetectorAnalysisResult into Nuclide and Dose rate analysis
//  types (have not seen RadAlarm, SpectrumPeakAnalysisResults or
//  GrossCountAnalysisResults results).
//  To nuclide result could add:
//    NuclideShieldingArealDensityValue, NuclideShieldingAtomicNumber,
//    NuclideActivityUncertaintyValue, SourcePosition
//    NuclideSourceGeometryCode, NuclideIDConfidenceValue, etc
//    as well as more better recording out SourcePosition.
//
//  To dose result could add:
//    AverageDoseRateUncertaintyValue, MaximumDoseRateValue,
//    MinimumDoseRateValue, BackgroundDoseRateValue,
//    BackgroundDoseRateUncertaintyValue, TotalDoseValue,
//    as well as more better recording out SourcePosition.
class DetectorAnalysisResult
{
public:
  std::string remark_;
  std::string nuclide_;
  float activity_;            //in units of becquerel (eg: 1.0 == 1 bq)
  std::string nuclide_type_;
  std::string id_confidence_;
  float distance_;            //in units of mm (eg: 1.0 == 1 mm )
  
  float dose_rate_;           //in units of micro-sievert per hour
                              //   (eg: 1.0 = 1 u-sv)
  
  float real_time_;           //in units of seconds (eg: 1.0 = 1 s)
  
  //20171220: removed start_time_, since there is no cooresponding equivalent
  //  in N42 2012.  Instead should link to which detectors this result is
  //  applicable to.
  //boost::posix_time::ptime start_time_;  //start time of the cooresponding data (its a little unclear in some formats if this is the case...)
  
  std::string detector_;
  
  DetectorAnalysisResult();
  void reset();

#if( PERFORM_DEVELOPER_CHECKS )
  static void equalEnough( const DetectorAnalysisResult &lhs,
                           const DetectorAnalysisResult &rhs );
#endif
};//struct DetectorAnalysisResult


/** A class that aims to eventually be about equivalent of the N42 2012 
    <AnalysisResults> tag.
 */
class DetectorAnalysis
{
public:
  //Need to make an association of this with the sample/detector number, and
  //  allow MeasurementInfo to have multiple of these...

  /** Remarks included with the AnalysisResults */
  std::vector<std::string> remarks_;
  
  /** A unique name of the analysis algorithm. */
  std::string algorithm_name_;
  
  /** Information describing the version of a particular analysis algorithm component. */
  std::string algorithm_version_;
  
  /** Creator or implementer of the analysis algorithm. */
  std::string algorithm_creator_;
  
  /** Free-form text describing the analysis algorithm. */
  std::string algorithm_description_;
  
  /** Free-form text describing the overall conclusion of the analysis regarding 
      the source of concern.  Equivalent to <AnalysisResultDescription> or 
      <ThreatDescription> tag of N42 2012 or 2006 respectively
   */
  std::string algorithm_result_description_;
  
  
  //N42 2012 also has the following fields:
  //AnalysisAlgorithmSetting
  //AnalysisResultStatusCode
  //AnalysisConfidenceValue
  //RadAlarm
  //NuclideAnalysisResults
  //SpectrumPeakAnalysisResults
  //GrossCountAnalysisResults
  //DoseAnalysisResults
  //ExposureAnalysisResults
  //Fault
  //AnalysisStartDateTime
  //AnalysisComputationDuration
  
  std::vector<DetectorAnalysisResult> results_;
  
  DetectorAnalysis();
  void reset();
  
#if( PERFORM_DEVELOPER_CHECKS )
  static void equalEnough( const DetectorAnalysis &lhs,
                           const DetectorAnalysis &rhs );
#endif
};//struct DetectorAnalysisResults



//implementation of inlined functions
bool MeasurementInfo::modified() const
{
  return modified_;
}

void MeasurementInfo::reset_modified()
{
  std::unique_lock<std::recursive_mutex> lock( mutex_ );
  modified_ = false;
}

void MeasurementInfo::reset_modified_since_decode()
{
  std::unique_lock<std::recursive_mutex> lock( mutex_ );
  modifiedSinceDecode_ = false;
}

bool MeasurementInfo::modified_since_decode() const
{
  return modifiedSinceDecode_;
}

float MeasurementInfo::gamma_live_time() const
{
  return gamma_live_time_;
}

float MeasurementInfo::gamma_real_time() const
{
  return gamma_real_time_;
}

double MeasurementInfo::gamma_count_sum() const
{
  return gamma_count_sum_;
}

double MeasurementInfo::neutron_counts_sum() const
{
  return neutron_counts_sum_;
}

const std::string &MeasurementInfo::filename() const
{
  return filename_;
}

const std::vector<std::string> &MeasurementInfo::detector_names() const
{
  return detector_names_;
}

const std::vector<int> &MeasurementInfo::detector_numbers() const
{
  return detector_numbers_;
}

const std::vector<std::string> &MeasurementInfo::neutron_detector_names() const
{
  return neutron_detector_names_;
}

const std::string &MeasurementInfo::uuid() const
{
  return uuid_;
}

const std::vector<std::string> &MeasurementInfo::remarks() const
{
  return remarks_;
}

int MeasurementInfo::lane_number() const
{
  return lane_number_;
}

const std::string &MeasurementInfo::measurement_location_name() const
{
  return measurement_location_name_;
}

const std::string &MeasurementInfo::inspection() const
{
  return inspection_;
}

double Measurement::latitude() const
{
  return latitude_;
}

double Measurement::longitude() const
{
  return longitude_;
}

bool Measurement::has_gps_info() const
{
  return (valid_longitude(longitude_) && valid_latitude(latitude_));
}

const boost::posix_time::ptime &Measurement::position_time() const
{
  return position_time_;
}

const std::string &MeasurementInfo::measurment_operator() const
{
  return measurment_operator_;
}

const std::set<int> &MeasurementInfo::sample_numbers() const
{
//  std::unique_lock<std::recursive_mutex> lock( mutex_ );
//  set<int> answer;
//  fore( const MeasurementShrdPtr &meas : measurements_ )
//    answer.insert( meas->sample_number_ );
//  return answer;
  return sample_numbers_;
}//set<int> sample_numbers() const


size_t MeasurementInfo::num_measurements() const
{
  size_t n;

  {
    std::unique_lock<std::recursive_mutex> lock( mutex_ );
    n = measurements_.size();
  }

  return n;
}//size_t num_measurements() const


std::shared_ptr<const Measurement> MeasurementInfo::measurement(
                                                             size_t num ) const
{
  std::unique_lock<std::recursive_mutex> lock( mutex_ );
  const size_t n = measurements_.size();
 
  if( num >= n )
    throw std::runtime_error( "MeasurementInfo::measurement(size_t): invalid index" );
  
  return measurements_[num];
}


DetectorType MeasurementInfo::detector_type() const
{
  return detector_type_;
}

const std::string &MeasurementInfo::instrument_type() const
{
  return instrument_type_;
}

const std::string &MeasurementInfo::manufacturer() const
{
  return manufacturer_;
}

const std::string &MeasurementInfo::instrument_model() const
{
  return instrument_model_;
}

const std::string &MeasurementInfo::instrument_id() const
{
  return instrument_id_;
}

std::vector< std::shared_ptr<const Measurement> > MeasurementInfo::measurements() const
{
  std::unique_lock<std::recursive_mutex> lock( mutex_ );
  
  std::vector< std::shared_ptr<const Measurement> > answer;
  for( size_t i = 0; i < measurements_.size(); ++i )
    answer.push_back( measurements_[i] );
  return answer;
}//std::vector< std::shared_ptr<const Measurement> > measurements() const


std::shared_ptr<const DetectorAnalysis> MeasurementInfo::detectors_analysis() const
{
  return detectors_analysis_;
}

bool MeasurementInfo::has_gps_info() const
{
  return (Measurement::valid_longitude(mean_longitude_)
           && Measurement::valid_latitude(mean_latitude_));
}

double MeasurementInfo::mean_latitude() const
{
  return mean_latitude_;
}

double MeasurementInfo::mean_longitude() const
{
  return mean_longitude_;
}

void MeasurementInfo::set_filename( const std::string &n )
{
  filename_ = n;
  modified_ = modifiedSinceDecode_ = true;
}

void MeasurementInfo::set_remarks( const std::vector<std::string> &n )
{
  remarks_ = n;
  modified_ = modifiedSinceDecode_ = true;
}

void MeasurementInfo::set_uuid( const std::string &n )
{
  uuid_ = n;
  modified_ = modifiedSinceDecode_ = true;
}

void MeasurementInfo::set_lane_number( const int num )
{
  lane_number_ = num;
  modified_ = modifiedSinceDecode_ = true;
}

void MeasurementInfo::set_measurement_location_name( const std::string &n )
{
  measurement_location_name_ = n;
  modified_ = modifiedSinceDecode_ = true;
}

void MeasurementInfo::set_inspection( const std::string &n )
{
  inspection_ = n;
  modified_ = modifiedSinceDecode_ = true;
}

void MeasurementInfo::set_instrument_type( const std::string &n )
{
  instrument_type_ = n;
  modified_ = modifiedSinceDecode_ = true;
}

void MeasurementInfo::set_detector_type( const DetectorType type )
{
  detector_type_ = type;
  modified_ = modifiedSinceDecode_ = true;
}

void MeasurementInfo::set_manufacturer( const std::string &n )
{
  manufacturer_ = n;
  modified_ = modifiedSinceDecode_ = true;
}

void MeasurementInfo::set_instrument_model( const std::string &n )
{
  instrument_model_ = n;
  modified_ = modifiedSinceDecode_ = true;
}

void MeasurementInfo::set_instrument_id( const std::string &n )
{
  instrument_id_ = n;
  modified_ = modifiedSinceDecode_ = true;
}




//implementation of inline Measurment functions
float Measurement::live_time() const
{
  return live_time_;
}

float Measurement::real_time() const
{
  return real_time_;
}

bool Measurement::contained_neutron() const
{
  return contained_neutron_;
}

int Measurement::sample_number() const
{
  return sample_number_;
}

Measurement::OccupancyStatus Measurement::occupied() const
{
  return occupied_;
}

double Measurement::gamma_count_sum() const
{
  return gamma_count_sum_;
}

double Measurement::neutron_counts_sum() const
{
  return neutron_counts_sum_;
}

float Measurement::speed() const
{
  return speed_;
}

const std::string &Measurement::detector_name() const
{
  return detector_name_;
}

int Measurement::detector_number() const
{
  return detector_number_;
}

const std::string &Measurement::detector_type() const
{
  return detector_type_;
}

Measurement::QualityStatus Measurement::quality_status() const
{
  return quality_status_;
}

Measurement::SourceType Measurement::source_type() const
{
  return source_type_;
}


Measurement::EquationType Measurement::energy_calibration_model() const
{
  return energy_calibration_model_;
}


const std::vector<std::string> &Measurement::remarks() const
{
  return remarks_;
}

const boost::posix_time::ptime &Measurement::start_time() const
{
  return start_time_;
}

const boost::posix_time::ptime Measurement::start_time_copy() const
{
  return start_time_;
}


const std::vector<float> &Measurement::calibration_coeffs() const
{
  return calibration_coeffs_;
}

const DeviationPairVec &Measurement::deviation_pairs() const
{
  return deviation_pairs_;
}

const std::shared_ptr< const std::vector<float> > &Measurement::channel_energies() const
{
  return channel_energies_;
}

const std::shared_ptr< const std::vector<float> > &Measurement::gamma_counts() const
{
  return gamma_counts_;
}

void Measurement::set_start_time( const boost::posix_time::ptime &time )
{
  start_time_ = time;
}

void Measurement::set_remarks( const std::vector<std::string> &remar )
{
  remarks_ = remar;
}

void Measurement::set_source_type( const SourceType type )
{
  source_type_ = type;
}


void Measurement::set_gamma_counts( ShrdConstFVecPtr counts,
                                    const float livetime, const float realtime )
{
  if( !counts )
  {
    if( !gamma_counts_ )
      return;
    gamma_counts_.reset( new std::vector<float>( gamma_counts_->size(), 0.0f ) );
    return;
  }//if( !counts )

  const size_t size = counts->size();
  gamma_counts_ = counts;
  live_time_ = livetime;
  real_time_ = realtime;
  gamma_count_sum_ = 0.0;

  const std::vector<float> &rhs = *counts;
  for( size_t i = 0; i < size; ++i )
    gamma_count_sum_ += rhs[i];
}


void Measurement::set_neutron_counts( const std::vector<float> &counts )
{
  neutron_counts_ = counts;
  neutron_counts_sum_ = 0.0;
  contained_neutron_ = !counts.empty();
  const size_t size = counts.size();
  for( size_t i = 0; i < size; ++i )
    neutron_counts_sum_ += counts[i];
}

void Measurement::set_channel_energies( ShrdConstFVecPtr counts )  //if channel_energies_ or gamma_counts_ must be same number channels as before
{
  if( !counts )
    return;

  if( (channel_energies_ && !channel_energies_->empty() && channel_energies_->size() != counts->size())
      || ( gamma_counts_ && !gamma_counts_->empty() && gamma_counts_->size() != counts->size() ) )
  {
    throw std::runtime_error( "Measurement::set_channel_energies(...): number of bin mismatch" );
  }//if( nbins mismatch )

  channel_energies_ = counts;
}

const std::vector<float> &Measurement::neutron_counts() const
{
  return neutron_counts_;
}


float Measurement::GetBinCenter( int bin ) const
{
  return (float)( GetBinLowEdge(bin) + 0.5*GetBinWidth(bin) );
}

float Measurement::GetBinContent( int bin ) const
{
  if( !CheckBinRange(bin) )
    return 0.0;
  return (*gamma_counts_)[bin-1];
}

float Measurement::GetBinLowEdge( int bin ) const
{
  if( !CheckBinRange(bin) )
    return (float)-999.9;
  return  (*channel_energies_)[bin-1];
}

float Measurement::GetBinWidth( int bin ) const
{
  if( !CheckBinRange(bin) )
    return 0.0;
  if( channel_energies_->size()==1 )
    return 0.0;
  if( bin == static_cast<int>(channel_energies_->size()) )
    --bin;
  return ((*channel_energies_)[bin]) - ((*channel_energies_)[bin-1]);
}

int Measurement::GetNbinsX() const
{
  if( !channel_energies_ )
  {
    if( gamma_counts_ )
      return static_cast<int>( gamma_counts_->size() );
    return 0;
  }//

  return static_cast<int>( channel_energies_->size() );
}


float Measurement::Integral( int binx1, int binx2 ) const
{
  float answer = 0.0;
  if( !gamma_counts_ )
    return answer;

  --binx1;
  --binx2;

  const int nbinsx = static_cast<int>( gamma_counts_->size() );
  
  if( binx1 < 0 )
    binx1 = 0;
  if( binx2 < 0 || binx2 < binx1 || binx2 >= nbinsx )
    binx2 = nbinsx - 1;

  for( int i = binx1; i <= binx2; ++i )
    answer += gamma_counts_->operator[](i);
  return answer;
}



size_t Measurement::num_gamma_channels() const
{
  if( !gamma_counts_ || !channel_energies_ )
    return 0;
  
  return std::min( gamma_counts_->size(), channel_energies_->size() );
}


size_t Measurement::find_gamma_channel( const float x ) const
{
  if( !channel_energies_ || channel_energies_->empty() )
    throw std::runtime_error( "Measurement::find_gamma_channel(): "
                              "channel_energies_ not defined" );
  
  //Using upper_bound instead of lower_bound to properly handle the case
  //  where x == bin lower energy.
  const std::vector<float>::const_iterator pos_iter
                              = std::upper_bound( channel_energies_->begin(),
                                                  channel_energies_->end(), x );
  
  if( pos_iter == channel_energies_->begin() )
    return 0;
  
  if( pos_iter == channel_energies_->end() )
    return channel_energies_->size() - 1;
  
  return (pos_iter - channel_energies_->begin()) - 1;
}//size_t find_gamma_channel( const float energy ) const


float Measurement::gamma_channel_content( const size_t channel ) const
{
  if( !gamma_counts_ || channel >= gamma_counts_->size() )
    return 0.0f;
  
  return gamma_counts_->operator[]( channel );
}//float gamma_channel_content( const size_t channel ) const


float Measurement::gamma_channel_lower( const size_t channel ) const
{
  if( !channel_energies_ || channel >= channel_energies_->size() )
    throw std::runtime_error( "Measurement::gamma_channel_lower(): "
                              "channel_energies_ not defined" );
  
  return channel_energies_->operator[]( channel );
}//float gamma_channel_lower( const size_t channel ) const


float Measurement::gamma_channel_center( const size_t channel ) const
{
  return gamma_channel_lower( channel ) + 0.5f*gamma_channel_width( channel );
}//float gamma_channel_center( const size_t channel ) const


float Measurement::gamma_channel_upper( const size_t channel ) const
{
  if( !channel_energies_ || channel >= channel_energies_->size()
     || channel_energies_->size() < 2 )
    throw std::runtime_error( "Measurement::gamma_channel_upper(): "
                             "channel_energies_ not defined" );
  
  if( channel < (channel_energies_->size()-1) )
    return (*channel_energies_)[channel+1];
  
  return gamma_channel_lower(channel) + gamma_channel_width(channel);
}//float gamma_channel_upper( const size_t channel ) const


inline const std::shared_ptr< const std::vector<float> > &
                                     Measurement::gamma_channel_energies() const
{
  return channel_energies_;
}


const std::shared_ptr< const std::vector<float> > &
                                    Measurement::gamma_channel_contents() const
{
  return gamma_counts_;
}


float Measurement::gamma_channel_width( const size_t channel ) const
{
  if( !channel_energies_ || channel >= channel_energies_->size()
      || channel_energies_->size() < 2 )
    throw std::runtime_error( "Measurement::gamma_channel_width(): "
                              "channel_energies_ not defined" );
  
  if( channel == (channel_energies_->size()-1) )
    return (*channel_energies_)[channel] - (*channel_energies_)[channel-1];
  
  return (*channel_energies_)[channel+1] - (*channel_energies_)[channel];
}//float gamma_channel_width( const size_t channel ) const



int Measurement::FindFixBin( float x ) const
{
  //Note, this function returns an index one-above what you would use to access
  //  the channel_energies_ or gamma_counts_ arrays.  This is to be compatible
  //  with CERNs ROOT package functions.
  if( !channel_energies_ )
    return -1;
  
  //Using upper_bound instead of lower_bound to properly handle the case
  //  where x == bin lower energy.
  const std::vector<float>::const_iterator pos_iter
                             = std::upper_bound( channel_energies_->begin(),
                                                 channel_energies_->end(), x );

  if( pos_iter == channel_energies_->end() )
  {
    float width = GetBinWidth( static_cast<int>(channel_energies_->size() ) );
    if( x < (channel_energies_->back()+width) )
      return static_cast<int>( channel_energies_->size() );
    else return static_cast<int>( channel_energies_->size() + 1 );
  }//if( pos_iter ==  channel_energies_->end() )

  if( pos_iter == channel_energies_->begin() )
  {
    if( x < channel_energies_->front() )
      return 0;
    return 1;
  }//if( pos_iter ==  channel_energies_->begin() )
  
  return static_cast<int>( pos_iter - channel_energies_->begin() );
}//int FindFixBin( float x ) const


bool Measurement::CheckBinRange( int bin ) const
{
  if( !channel_energies_
      || bin < 1 || bin > static_cast<int>(channel_energies_->size()) )
    return false;
  return true;
}


const std::string &Measurement::title() const
{
  return title_;
}


void Measurement::set_title( const std::string &title )
{
  title_ = title;
}




float Measurement::gamma_energy_min() const
{
  if( !channel_energies_ || channel_energies_->empty() )
    return 0.0f;
  return (*channel_energies_)[0];
}

float Measurement::gamma_energy_max() const
{
  if( !channel_energies_ || channel_energies_->empty() )
    return 0.0f;
  
  const size_t nbin = channel_energies_->size();
  if( nbin < 2 )
    return (*channel_energies_)[0];
  
  return 2.0f*(*channel_energies_)[nbin-1] - (*channel_energies_)[nbin-2];
}

#endif
