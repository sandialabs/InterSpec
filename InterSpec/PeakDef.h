#ifndef PeakDef_h
#define PeakDef_h
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
#include <memory>
#include <vector>
#include <utility>
#include <assert.h>
#include <stdexcept>

#include <Wt/WColor>

#include "InterSpec/ReactionGamma.h"

//Forward declaration
class PeakDef;
namespace SpecUtils{ class Measurement; }
namespace SandiaDecay
{
  struct Nuclide;
  struct Element;
  struct Transition;
  struct RadParticle;
  struct EnergyIntensityPair;
}//namespace SandiaDecay


namespace rapidxml
{
  template<class Ch> class xml_node;
  template<class Ch> class xml_document;
}//namespace rapidxml


/**
 
 TODO: When PeakDefs are copied, their continuum still points to the same PeakContinuum object,
       meaning if the copies modify their continuum, the originals continuum also gets modified.
       This is currently necessary since multiple PeakDefs can share a PeakContinuum.  To fix this
       a ROI (e.g., continuum) should own the peaks, not the other way around
 */
struct PeakContinuum
{
  enum OffsetType : int
  {
    /** The gaussian will sit on top of y=0.
     No parameters.
     */
    NoOffset,
    
    /** The gaussian sits on top of a flat, constant y-offset.
     One parameter - the number of counts in the continuum per keV.  E.g., to get counts in a channel, multiply the continuum parameter
     value by the channel width, in keV.
     
     Note: This, and the other polynomial enums, have a value equal to the expected number of parameters to describe them
     */
    Constant = 1,
    
    /** The gaussian sits on top of a straight, but sloped line.
     Two parameters.
     The density of continuum counts (e.g., counts/keV) is given by dy = par[0] + (energy - ref_energy)*par[1], where 'ref_energy' is
     usually the mean of the first peak in the ROI.  A 'ref_energy' is used, instead of absolute energy, for numerical stability.
     To get the counts in a channel you must integrate the density, e.g.,
       y_chan = par[0]*(E_1 - E_0) + par[1]*(E_1^2 - E_0^2 + 2*ref_energy*(E_0 - E_1))
               where E_0 and E_1 is the lower and upper energies of the channel respectively.
     */
    Linear,
    
    /** Similar to #OffsetType::Linear, but a polynomial with three parameters. */
    Quadratic,
    
    /** Similar to #OffsetType::Linear, but a polynomial with four parameters. */
    Cubic,
    
    /** Two parameters; first gives offset as #OffsetType::Constant, with the second parameter giving change in dy (density of counts
     per keV) between the left and right side of the ROI - the change in count density varies across the ROI as the fraction of data in the
     ROI up to that point (if that makes any sense).
     */
    FlatStep,
    
    LinearStep,
    
    BiLinearStep,
    
    /** A continuum is determined algorithmically for the entire spectrum; not recommended to use. */
    External
  };//enum OffsetType
  
  /** Text appropriate for use as a label for the continuum type in the gui. */
  static const char *offset_type_label( const OffsetType type );
  
  /** Returns string to be used for XML or JSON identification of continuum type. */
  static const char *offset_type_str( const OffsetType type );
  
  /** Returns the number of parameters for a specified offset type. */
  static size_t num_parameters( const OffsetType type );
  
  /** Returns true if continuum type is FlatStep, LinearStep, or BiLinearStep and otherwise false */
  static bool is_step_continuum( const OffsetType type );
  
  /** Throws exception if string does not match a string returned by #offset_type_str.
   @param str String to be tested.  Must not be a null pointer, but string does not need to be null terminated.
   @param len The length of the string to be tested.
   
   Throws exception if an invalid string.
   */
  static OffsetType str_to_offset_type_str( const char * const str, const size_t len );
  
  PeakContinuum();
  
  //setType: makes sure things are in a consistent state.  E.g. proper number
  //  of polynomial coefficients, or that the external continuum is kept around
  //  when not needed.
  void setType( OffsetType type );
  OffsetType type() const { return m_type; }
  
  //setParameters: throws if NoOffset, External, or wrong number of parameters.
  //  If uncertainties is empty, will reset uncertainties to zero
  // - should consider renaming setPolynomialParameters(...)
  void setParameters( double referenceEnergy,
                     const std::vector<double> &parameters,
                     const std::vector<double> &uncertainties );

  //setParameters: a convenience function which actually calls other form of
  //  this member function.  This function assumes both passed in arrays are
  //  of a length equal to OffsetType; uncertainties may be a null pointer.
  void setParameters( double referenceEnergy,
                      const double *parameters,
                      const double *uncertainties );

  //setPolynomialCoefFitFor: sets whether the polynomial coefficient should be
  //  fit for.  Returns if the coefficient was set; returns false if you call it
  //  with an invalid index
  bool setPolynomialCoefFitFor( size_t polyCoefNum, bool fit );

  //setPolynomialCoef: sets the polynomial coefficient value. Returns if the
  //  coefficient was set; returns false if you call it with an invalid index
  bool setPolynomialCoef( size_t polyCoef, double val );
  
  //setPolynomialUncert: analagous to setPolynomialCoef(...), but for uncert.
  bool setPolynomialUncert( size_t polyCoef, double val );
  
  double referenceEnergy() const { return m_referenceEnergy; }
  const std::vector<double> &parameters() const { return m_values; }
  const std::vector<double> &uncertainties() const { return m_uncertainties; }
  std::vector<bool> fitForParameter() const { return m_fitForValue; }
  std::shared_ptr<const SpecUtils::Measurement> externalContinuum() const { return m_externalContinuum; }
  
  //setGlobalContinuum: throws if not a External OffsetType.
  void setExternalContinuum( const std::shared_ptr<const SpecUtils::Measurement> &data );
 
  //setRange: sets the energy range this continuum is applicable for
  void setRange( const double lowerenergy, const double upperenergy );
  

  /** Sets this to be a Linear OffsetType continuum, and the range to be from x0 to x1, with the resulting linear polynomial
     relative to m_referenceEnergy == x0.
   
   \param data The spectrum to use
   \param reference_energy The reference energy the equation will be based on.  Must be equal to, or between \p roi_lower
                           and \p roi_upper.
   \param roi_lower The lower energy of the ROI
   \param roi_upper The upper energy of the ROI.  Must be greater than \p roi_lower.
   \param num_lower_channels The number of channels below the ROI to use for calculating the equation; must be 1 or larger
   \param num_upper_channels The number of channels above the ROI to use for calculating the equation; must be 1 or larger
   
   Note that if \p roi_lower and/or \p roi_upper dont correspond to a bin edge, then the bordering bin will be skipped if the
   ROI extends more than 10% of the way into that bin.
   
   Will throw exception on input error, such as
   */
  void calc_linear_continuum_eqn( const std::shared_ptr<const SpecUtils::Measurement> &data,
                                 const double reference_energy,
                                 const double roi_lower, const double roi_upper,
                                 const size_t num_lower_channels,
                                 const size_t num_upper_channels );
  
  //offset_integral: returns the area of the continuum from x0 to x1.  If m_type
  //  is NoOffset, then will return 0.  If a polynomial, then the integral of
  //  the contonuum density will be returned.  If globally defined continuum,
  //  then the integral will be returned, using linear interpolation for
  //  ranges not exactly on the bin edges.
  //  If FlatStep, then data histogram is necassary
  double offset_integral( const double x0, const double x1,
                          const std::shared_ptr<const SpecUtils::Measurement> &data ) const;
  
  /** Returns true if a _valid_ polynomial, step, or external continuum type.
   Where valid polynomial and step continuums means any of the coefficients are non zero, not that they actually make sense.
   */
  bool parametersProbablySet() const;
  
  //energyRangeDefined: returns if an energy range has explicitely been set
  bool energyRangeDefined() const;
  
  /** Returns true if Constant, Linear, Quadratic, Cubic, or FlatStep: */
  bool isPolynomial() const;
  
  double lowerEnergy() const { return m_lowerEnergy; }
  double upperEnergy() const { return m_upperEnergy; }
  
  //eqn_from_offsets: uses bin heights at lowbin and highbin to calculate the
  //  linear continuum density coefficients, relative to referenceEnergy.
  //  nbinEachSide: the number of bin on each side of lowbin/highbin to use
  //                toward the continuum hight as well (0 means use just
  //                lowbin/highbin, 1 means use 3 bins centered around
  //                lowbin/highbin, etc)
  static void eqn_from_offsets( size_t lowchannel,
                                size_t highchannel,
                                const double referenceEnergy,
                                const std::shared_ptr<const SpecUtils::Measurement> &data,
                                const size_t num_lower_channels,
                                const size_t num_upper_channels,
                                double &m, double &b );
  
  //offset_eqn_integral: integrates the continuum density, specified by coefs,
  //  to return the continuum area.
  static double offset_eqn_integral( const double *coefs,
                                     OffsetType type,
                                     double x0, double x1,
                                     const double reference_energy );
  
  //translate_offset_polynomial: if you want the exact same polynomial
  //  line, but would like it with reference to a different energy,
  //  use this function.
  //  TODO: Does NOT support Quadratic, Cubic, or External continuum types 
  //  TODO: Currently untested
  static void translate_offset_polynomial( double *new_coefs,
                                           const double *old_coefs,
                                           OffsetType type,
                                           const double new_reference_energy,
                                           const double old_reference_energy );
  bool operator==( const PeakContinuum &rhs ) const;
 
  //toXml: XXX = does not support serializing of m_externalContinuum is pretty 
  //  bad currently!  Currently serializes a new continuum for each peak, even
  //  if they should be shared across peaks
  void toXml( rapidxml::xml_node<char> *parent, const int contId ) const;
  void fromXml( const rapidxml::xml_node<char> *node, int &contId );

  
#if( PERFORM_DEVELOPER_CHECKS )
  //equalEnough(...): tests whether the passed in PeakDef objects are
  //  equal, for most intents and purposes.  Allows some small numerical
  //  rounding to occur.
  //Throws an std::exception with a brief explanaition when an issue is found.
  static void equalEnough( const PeakContinuum &lhs, const PeakContinuum &rhs );
#endif

  
protected:
  OffsetType m_type;
  
  //m_lowerEnergy, m_upperEnergy: will be 0.0 if undefined.  Also if
  //  m_lowerEnergy==m_upperEnergy, then will be considered undefined as well.
  double m_lowerEnergy, m_upperEnergy;
  
  //m_referenceEnergy: the energy polynomials are evaluated relative to. Ex.
  //  continuum_density = m_values[0] + (energy - m_referenceEnergy)*m_values[1];
  double m_referenceEnergy;
  
  //m_values, m_uncertainties: the polynomial coefficients and their
  //  uncertainties.  The polynomial is actually a density per keV of the
  //  continuum, so needs to be integrated over the width of the bin.  Also,
  //  the energy is relative to m_referenceEnergy. E.g. at the m_referenceEnergy
  //  the continuum will have a density of m_values[0]
  std::vector<double> m_values, m_uncertainties;
  std::vector<bool> m_fitForValue;
  
  //m_externalContinuum: since the global continuum may be changed after fitting
  //  for a peak (thus possibly resulting in a poorly fit peak), we'll keep a
  //  refernce to the original.  This is slightly inefficient, but shouldnt be
  //  to bad.
  //  \TODO: need to verify when writing a SpecMeas object to a native file, if
  //         multiple peaks share a global continuum, it is only written once in
  //         the file.
  //  \TODO: use #Measurement::set_energy_calibration when translating a peaks mean for energy
  //         calibrations.
  std::shared_ptr<const SpecUtils::Measurement> m_externalContinuum;
  
  static const int sm_xmlSerializationVersion;
  
  friend std::ostream &operator<<( std::ostream &, const PeakContinuum & );
};//struct PeakContinuum



/** TODO: Define and implement this function
 
 @param fwhm The real or estimated FWHM of the peak
 @param energy The mean of the peak
 @param nchannel The number of channels in the spectrum.  If specified as zero, will not be used.
 @returns Whether this peak is likely to be from a HPGe detector or not.
 */
//bool is_high_resolution_peak( const float fwhm, const float energy, const size_t nchannel );
//
// OR
/** Returns true if the energy per channel is more inline with a high resolution system.
 */
//bool is_high_resolution( const std::shared_ptr<const SpecUtils::Measurement> &data );
//


//findROILimit(...): returns the channel number that should be the extent of the
//  region of interest, for the peak with the given mean and sigma.  If
//  high==true, then the bin higher in energy than the mean will be returned,
//  otherwise it will be the bin on the lower energy side of the mean.
//  The basic idea is to include a maximum of 11.75 sigma away from the mean
//  of the peak, but then start at ~1.5 sigma from mean, and try to detect if
//  a new feature is occuring, and if so, stop the region of interest there.
//  A feature is "detected" if the value of the bin contents exceeds 2.5 sigma
//  from the "expected" value, where the expected value starts off being the
//  smallest bin value so far (well, this bin averaged with the bins on either
//  side of it), and then each preceeding bin is added to this background
//  value.
//  This is kinda similar in principle to how PCGAP does it, but with further
//  modifications, and I'm sure some differences; see,
//  http://www.inl.gov/technicalpublications/Documents/3318133.pdf
size_t findROILimit( const PeakDef &peak, const std::shared_ptr<const SpecUtils::Measurement> &data, bool high);


//findROIEnergyLimits(...): a convience function that calls findROILimit(...)
//  inorder to set 'lowerEnengy' and 'upperEnergy'.  IF data is an invalid
//  pointer, then PeakDef::lowerX() and PeakDef::upperX() are used.
void findROIEnergyLimits( double &lowerEnengy, double &upperEnergy,
                         const PeakDef &peak, const std::shared_ptr<const SpecUtils::Measurement> &data );

//Returns if area from start2 to end2 is greater than or equal to the
//  area of (start1 to end1 minus nsigma).
//Answers the question: is region 2 similar or greater in area as region 1
bool isStatisticallyGreaterOrEqual( const size_t start1, const size_t end1,
                                    const size_t start2, const size_t end2,
                                    const std::shared_ptr<const SpecUtils::Measurement> &dataH,
                                    const double nsigma );

//If a polynomial continuum and the peakdoesnt specify the range, then it
//  will look at the data histogram for when data starts increasing, and set
//  this as the limit.
//Does nothing if the histogram is null.
void estimatePeakFitRange( const PeakDef &peak, const std::shared_ptr<const SpecUtils::Measurement> &dataH,
                          size_t &lower_channel, size_t &upper_channel );


class PeakDef
{
  // If I was to re-write this class from scratch I would define peak mean by channel, and FWHM
  //  by the fraction of the peak mean, and similarly for the other quantities so that it is
  //  invariant to the energy calibration - although this isnt without issues.
  // I would also make the PeakContinuum the primary ROI with a list of peaks belonging to it,
  //  instead of the other way around.
  
public:
  /** The skew type applied with Gaussian defined peaks.
   
   The skew types are listed roughly in order of what should be preferred, if it can describe the data.
   */
  enum SkewType
  {
    /** No skew - i.e., a purely Gaussian distribution. */
    NoSkew,
  
    /** The Bortel function, from the paper referenced below, is:
     
     Convolution of Gaussian with an left-hand exponential multiplied by a step function
     that goes to zero above the peak mean.  The Bortel paper cited below uses two
     exponentials, but we use only one for gamma spectroscopy.
     
     See: Analytical function for fitting peaks in alpha-particle spectra from Si detectors
     G. Bortels, P. Collaers
     International Journal of Radiation Applications and Instrumentation. Part A. Applied Radiation and Isotopes
     Volume 38, Issue 10, 1987, Pages 831-837
     https://doi.org/10.1016/0883-2889(87)90180-8
     
     Uses one skew parameter.
     
     Maybe call exGaussian?
     */
    Bortel,
    
    /** A model of an Doniach Sunjic asymmetric lineshape.
     
     Commonly use in X-ray Photoelectron Spectroscopy (XPS) peak fitting.
     
     See https://lmfit.github.io/lmfit-py/builtin_models.html#doniachmodel
     
     Uses one skew parameter.
     
     Not currently defined because the amplitude can be infinite for non-zero skew, meaning we need some-way to specify the reported range - likely is
     */
    // Doniach,
    
    /** An exponential tail stitched to a Gaussian core.
     
     See: A simple alternative to the Crystal Ball function.
     Souvik Das, arXiv:1603.08591
     https://arxiv.org/abs/1603.08591
     
     Uses one skew parameter.
     */
    GaussExp,
    
    /** A Gaussian core portion and a power-law low-end tail, below a threshold.
     
     See https://en.wikipedia.org/wiki/Crystal_Ball_function
     
     Uses two skew parameters.
     The first is `alpha`, which defines the threshold.
     The second is `n` which defines the power-law.
     */
    CrystalBall,
    
    /** The GausExp extended to a exponential tail on each side of the Gaussian.
     
     See same reference as for GaussExp.
     
     Uses two skew parameters.
     The first one is the skew on the left.
     The second one is the skew on the right.
     */
    ExpGaussExp,
    
    /** A "double-sided" version of the Crystal Ball distribution to account for high-energy skew.
     
     See: Search for resonances in diphoton events at $$ \sqrt{s}=13 $$ TeV with the ATLAS detector
     Aaboud, M., G. Aad, B. Abbott, J. Abdallah, O. Abdinov, B. Abeloos, R. Aben, et al. 2016
     http://nrs.harvard.edu/urn-3:HUL.InstRepos:29362185
     Chapter 6
     
     Note also that CERNs ROOT package implements this function, with likely with some consideration
     for power laws between 1 and 1.00001, that we dont currently do.  However, this implementation
     is independent of the ROOT implementation, and treats the distribution as having unit area
     when integrated over all `x`, which it doesnt appear the ROOT implementation does (but I didnt
     check)
     https://root.cern.ch/doc/master/RooCrystalBall_8cxx_source.html
     
     Uses four skew parameters.
     The first is `alpha_low`, which defines the threshold, on the left side.
     The second is `n_low` which defines the power-law on the left side.
     The third is `alpha_high`, which defines the threshold, on the right side.
     The second is `n_high` which defines the power-law on the right side.
     */
    DoubleSidedCrystalBall
  };//enum SkewType
  
  enum DefintionType
  {
    GaussianDefined,
    DataDefined
  };//enum PeakType
  
  enum CoefficientType
  {
    Mean,
    Sigma,
    GaussAmplitude,
    SkewPar0,
    SkewPar1,
    SkewPar2,
    SkewPar3,
    Chi2DOF,           //for peaks that share a ROI/Continuum, this values is for entire ROI/Continuum
    NumCoefficientTypes
  };//enum CoefficientType

  enum SourceGammaType
  {
    NormalGamma,
    AnnihilationGamma,
    SingleEscapeGamma,
    DoubleEscapeGamma,
    XrayGamma
  };//enum SourceGammaType
  
  //gammaTypeFromUserInput: guesses SourceGammaType from txt string.
  //  After calling txt will not contain the dingle/double escape peak text.
  //  strings recognized as something other than DecayGamma, are:
  //  'se ', 's.e.', 'de', 'de ', 'single escape', 'double escape'
  static void gammaTypeFromUserInput( std::string &txt, SourceGammaType &type );
  
  /** Function to extract energy from a peaks source string, such as "U238 185.7 keV".
   On success, returns energy extracted, as well as modifies the input string to remove the portion of the string that specified the energy.
   On failure, returns -1.0 and does not modify input string.
   Input is case insensitive.
   
   Examples input strings that will be successful:
    "fe xray 98.2 kev"         -> {98.2, "fe xray"},
    "5.34e+2 keV"              -> {534, ""},
    "hf178m 5.34e-3 Mev" -> {5.34, "hf178m"},
    "8.0E+02 kev hf178m" -> {800, "hf178m"},
    "hf178m2 574."            -> {574, "hf178m2"},
    "u232 98"                     -> {98, "u232"},
    "3.3mev be(a,n)"          -> {3300, "be(a,n)"},
    "co60 1173.23"            -> {1173.23, "co60"},
    "co60 1173.23 kev"     -> {1173.23, "co60"},
    "Pb 98.2"                     -> {98.2,"Pb"}
   Examples input strings that will NOT be successful:
    "98 u232", "u-232", "321 u-232", "1173.0 CO60", "Pb 98"
   */
  static double extract_energy_from_peak_source_string( std::string &str );
  
  
  static const char *to_string( const CoefficientType type );
  static const char *to_string( const SkewType type );
  static SkewType skew_from_string( const std::string &skew_str );
  
  /** Gives reasonable range for skew parameter values, as well as a reasonable starting value.
   
   Returns true if limits if the #CoefficientType is applicable to #SkewType, or else returns false and
   sets values to all zero
   */
  static bool skew_parameter_range( const SkewType skew_type, const CoefficientType coefficient,
                                   double &lower_value, double &upper_value,
                                   double &starting_value, double &step_size );
  
  /** Returns the number of parameters required to describe the skew type.
   Will return 0, 1, 2, 3, or 4
   */
  static size_t num_skew_parameters( const SkewType skew_type );
  
  /** Returns if a skew parameter has an energy dependence across the spectrum, or if the parameter should have
   the same value, no matter the energy of the peak.
   
   This is currently experimental - I dont actually know!!!
   */
  static bool is_energy_dependent( const SkewType skew_type, const CoefficientType coefficient );
  
public:
  PeakDef();

  //The following constructor constructs a gaussian peak with no skew or
  //  self defined background/continuum.
  PeakDef( double m, double s, double a );

  //The following constructor make a peak were m_PeakDef==false and
  //  instead the peak will be drawn as the difference between data and
  //  background.  If a background histogram is not proivided than a
  //  strait line spanning lowerx and upperx will be used.
  //This usage of PeakDef has not been well tested
  PeakDef( double lowerx, double upperx, double mean,
            std::shared_ptr<const SpecUtils::Measurement> data, std::shared_ptr<const SpecUtils::Measurement> background );

//  PeakDef( const PeakDef &rhs );
//  const PeakDef &operator=( const PeakDef &rhs );

  //continuum(): access the continuum
  std::shared_ptr<PeakContinuum> continuum();
  std::shared_ptr<const PeakContinuum> continuum() const;
  
  //getContinuum(): garunteed to be a valid pointer
  std::shared_ptr<PeakContinuum> getContinuum();
  
  //setContinuum(...): input must be a valid pointer or exception will be thrown
  void setContinuum( std::shared_ptr<PeakContinuum> continuum );
  
  //makeUniqueNewContinuum(): clones the current continuum, and sets it to be
  //  the one used; this is handy for instances when you want to be sure no
  //  other peaks share the continuum
  void makeUniqueNewContinuum();
  
  void reset();

  inline double mean() const;
  inline double sigma() const;  //throw exception if not a gausPeak
  inline double fwhm() const;
  inline double amplitude() const;
  inline bool gausPeak() const;
  
  //The below should in principle take care of gaussian area and the skew area
  double peakArea() const;
  double peakAreaUncert() const;
  
  //setPeakArea(): sets the area of the Gaussian + Skew (not just the gaussian
  //  component)
  void setPeakArea( const double a );
  void setPeakAreaUncert( const double a );
  
  /** Returns the area between the continuum and data, everywhere data is above
     continuum in the ROI.
   */
  double areaFromData( std::shared_ptr<const SpecUtils::Measurement> data ) const;
  
  
  inline double meanUncert() const;
  inline double sigmaUncert() const;
  inline double amplitudeUncert() const;


  inline void setMean( const double m );
  inline void setSigma( const double s );
  inline void setAmplitude( const double a );

  inline void setMeanUncert( const double m );
  inline void setSigmaUncert( const double s );
  inline void setAmplitudeUncert( const double a );


  inline bool useForEnergyCalibration() const;
  inline void useForEnergyCalibration( const bool use );

  inline bool useForShieldingSourceFit() const;
  inline void useForShieldingSourceFit( const bool use );
  
  inline bool useForManualRelEff() const;
  inline void useForManualRelEff( const bool use );
  
  /** Returns if should use for DRF intrinsic efficiency fit.  Note that this does not check that the nuclide and transition has actually
   been defined.
   */
  inline bool useForDrfIntrinsicEffFit() const;
  inline void setUseForDrfIntrinsicEffFit( const bool use );
  inline bool useForDrfFwhmFit() const;
  inline void setUseForDrfFwhmFit( const bool use );
  inline bool useForDrfDepthOfInteractionFit() const;
  inline void setUseForDrfDepthOfInteractionFit( const bool use );
  
  
  inline double chi2dof() const;
  inline bool chi2Defined() const;

  inline const std::string &userLabel() const;
  inline void setUserLabel( const std::string &utf8Label );
  
  /** Remove any Nuclide, Reaction, or Xray associated with the peak.
   After calling this function parentNuclide(), decayParticle(),
   nuclearTransition(), xrayElement() and reaction() will all return nullptr.
   */
  void clearSources();
  
  /** Check if nuclide, xray, or reaction has been set. */
  bool hasSourceGammaAssigned() const;
  
  //setNuclearTransition(...): sets the inforation about what nuclide is
  //  responsible for this peak.  The transition and radParticle index specifies
  //  which decay is responsible for the gamma, however, if isAnnihilationGamma
  //  is specified, then this is a 511 keV gamma from positron annihilation.
  //  Not that the particle may be an xray (caused by decay of nucleus)
  void setNuclearTransition(  const SandiaDecay::Nuclide *parentNuclide,
                              const SandiaDecay::Transition *transition,
                              const int radParticle,
                              const SourceGammaType sourceType );
  
  //parentNuclide(): returns the pointer to the nuclide of the parent
  //  responsible for creating this peak, although one of its prodigeny may
  //  have been the acutal nuclide that gave off the gamma (which you can get
  //  by decayParticle()->parent).
  //  If the parent nuclide, then the decayParticle() will give a non-null
  //  result _or_ isAnnihilationGamma() will be true.
  const SandiaDecay::Nuclide *parentNuclide() const;
  
  //decayParticle(): may return null even if parentNuclide() did not in the case
  //  of the gamma being an annihilation gamma (from a positron), in which
  //  case sourceGammaType() should return AnnihilationGamma.
  //  This being null and sourceGammaType() giving AnnihilationGamma type
  //  indicates there were mutliple positrons given off; if there is only one
  //  positron in the decay, than decayParticle() should be non-null, but
  //  you should not use SandiaDecay::RadParticle::energy as the photon energy,
  //  but 511 keV (which is what gammaParticleEnergy() will return) instead.
  //  if sourceGammaType() is Songle/DoubleEscapeGamma then you should subtract
  //  511/1022 keV from SandiaDecay::RadParticle::energy (which is what
  //  gammaParticleEnergy() will return)
  const SandiaDecay::RadParticle *decayParticle() const;
  
  //nuclearTransition(): returns the actual transition responsible for creating
  //  this peak.  May be null even if parentNuclide() is not, in the case
  //  of multiple electron capture decays from the parent, and the 511 keV gamma
  //  is responsible for this peak.
  const SandiaDecay::Transition *nuclearTransition() const;
  
  //decayParticleIndex(): this index in SandiaDecay::Transition::products of
  //  the particle responsible for this peak.
  int decayParticleIndex() const;
  
  //sourceGammaType():
  SourceGammaType sourceGammaType() const;
  
  void setXray( const SandiaDecay::Element *el, const float energy );
  const SandiaDecay::Element *xrayElement() const;
  float xrayEnergy() const;
  
  //setReaction:  If a valid rctn, then xray and isotopes will be un-set.
  //  An exception will be thrown if sourceType is PeakDef::AnnihilationGamma.
  void setReaction( const ReactionGamma::Reaction *rctn, const float energy,
                    const SourceGammaType sourceType );
  
  const ReactionGamma::Reaction *reaction() const;
  float reactionEnergy() const;
  
  //gammaParticleEnergy(): returns the energy of the gamma/xray responsible for
  //  this peak, despite if it is a normal nuclear transition, an annihilation,
  //  an xray, or a reaction gamma.  For single and double excape peaks, the 511
  //  or 1022 keV is subtracted off of the SandiaDecay::RadParticle::energy
  //  or reactionEnergy() that gave rise to this escape peak.
  //Throws an exception if there is no gamma associated with this peak.
  float gammaParticleEnergy() const;

  struct CandidateNuclide
  {
    const SandiaDecay::Nuclide *nuclide;
    const SandiaDecay::Transition *transition;
    int radparticleIndex; //Index of object in transition->products
    float weight;
    PeakDef::SourceGammaType sourceGammaType;
  };//struct CandidateNuclide

  const std::vector< CandidateNuclide > &candidateNuclides() const;
  void addCandidateNuclide( const CandidateNuclide &candidate );
  void addCandidateNuclide( const SandiaDecay::Nuclide * const nuclide,
                            const SandiaDecay::Transition * const transition,
                            const int radparticleIndex,
                            const SourceGammaType sourceGammaType,
                            const float weight );
  void setCandidateNuclides( const std::vector<CandidateNuclide> &candidates );


  /** Sets values of quantities/items that are not fit from the data.
   So, whether to use for energy calibration, peak color, wether a quantity should be fit for, etc.
   
   @param parent The peak whose values should be copied
   @param inheritNonFitForValues If true, then quantities that currently are selected to not be fit from data, but in principle could
          be, will also be copied.  So for example, if the peak-mean is selected to not be fit for, and this parameter is true, then the
          mean from \p parent will also be copied to *this.
   
   This function does not modify continuum extent, reference energy, or polynomial type, but does copy if polynomial coefficients should
   be fit for, only (or if different polynomial order continuums, copies up to the lesser order of the two continuums).
   
   TODO: consider if more continuum information should be copied, like polynomial values not being
         fit for, or continuum type.
   */
  void inheritUserSelectedOptions( const PeakDef &parent,
                                  const bool inheritNonFitForValues );

  double lowerX() const;
  double upperX() const;
  inline double roiWidth() const;
  
  /** Returns the the area of the Gaussian and Skew (if applicable) components of the peak,
   between x0 and x1.
   */
  double gauss_integral( const double x0, const double x1 ) const;
  
  /** Adds each channels contribution from the Gaussian and Skew (if applicable) into a channel
   array.
   
   The computation of the Gaussian integral calls the `erf` function twice (once for lower energy,
   and once for upper energy of each channel), which even when using an optimized function, takes
   up much of the CPU time of fitting peaks.  This call to get the peak area contributions to each
   channel effectively cuts the number of calls to the `erf` function in half.
   
   \param energies Array of lower channel energies; must have at least one more entry than
          `nchannel`
   \param channels Channel count array integrals of Gaussian and Skew will be _added_ to (e.g.,
          will not be zeroed); must have at least `nchannel` entries
   \param nchannel The number of channels to do the integration over.
   */
  void gauss_integral( const float *energies, double *channels, const size_t nchannel ) const;
  
  
  //offset_integral(): gives area of the continuum component between x0 and x1.
  double offset_integral( const double x0, const double x1,
                          const std::shared_ptr<const SpecUtils::Measurement> &data ) const;

  inline bool fitFor( CoefficientType type ) const;
  inline void setFitFor( CoefficientType type, bool fit );
  inline const bool *fitFors() const;
  
  inline double coefficient( CoefficientType type ) const;
  inline double uncertainty( CoefficientType type ) const;
  
  inline const double *coefficients() const;
  inline const double *uncertainties() const;
  
  inline void set_coefficient( double val, CoefficientType type );
  inline void set_uncertainty( double val, CoefficientType type );

  inline DefintionType type() const;
  
  inline SkewType skewType() const;
  inline void setSkewType( SkewType );
  
  /** Returns currently assigned color of the peak. Wt::WColor::isDefault()==true
      indicates no color set.
   */
  const Wt::WColor &lineColor() const;
  
  /** Sets the CSS style color of peak.
   */
  void setLineColor( const Wt::WColor &color );

  static bool lessThanByMean( const PeakDef &lhs, const PeakDef &rhs );
  static bool lessThanByMeanShrdPtr( const std::shared_ptr<const PeakDef> &lhs,
                                     const std::shared_ptr<const PeakDef> &rhs );
  
  bool operator==( const PeakDef &rhs ) const;

  
  static bool causilyConnected( const PeakDef &lower_peak,
                                const PeakDef &upper_peak,
                                const double ncausality,
                                const bool useRoiAsWell );
  
  //causilyDisconnected: checks if peaks plus/minus ncausality sigma overlaps.
  //  If peak is not gaussian peak, then uses ROI extent and ignores ncausality
  //  If 'useRoiAsWell' is specified, then this function will use the greater
  //  of ncausality*sigma, or continuum()->lowerX()/upperX().
  static bool causilyDisconnected( const PeakDef &lower_peak,
                                     const PeakDef &upper_peak,
                                     const double ncausality,
                                     const bool useRoiAsWell );
  
  /** Use heuristics to guess at a reasonable age a nuclide might be if
   * encountered in the real world.  
   *
   * \param nuc Nuclide to return the age for.
   * \param strAnswer If provided, will assign a suggested string for the age
   *        to this string.  May yield either an actual unit of time (ex.
   *        "2.31 years") or a multiple of half lives (ex. "3.5 HL").
   *        Will have a maximum of 2 decimal points.
   * \returns Suggested age= value in PhysicalUnits units.
   *
   * \sa defaultDecayTimeString
   */
  static double defaultDecayTime( const SandiaDecay::Nuclide *nuc,
                                  std::string *strAnswer = NULL );
  
  /** Lets not both changing/fitting ages for nuclides like
  * Cs137 that dont effectively change spectrum over the course of a few days.
  * Right now it catches simple cases where we shouldnt try to fit for age,
  * but it definetly misses a number of nuclides probably
  */
  static bool ageFitNotAllowed( const SandiaDecay::Nuclide *nuc );
  
  
  //XXX - Right now findNearestPhotopeak(...) decays the nuclide you pass in
  //      to 5 {prompt, if unexistent secular, if unexistent regular} half lifes
  //      then uses then matches the energy passed in to one of the lines - this
  //      could be improved by removing lines due to protigen from after the
  //      nuclide which defines the prompt quilibrium decay endpoint - but
  //      is this what we reeally want?
  //
  //findNearestPhotopeak(): doenst find the nearest photopeak, but instead
  //  the one with smallest
  //    FOM=(0.1*windowHalfWidth + abs(photopeak-energy) ) / photopeak_intensity
  //
  //  windowHalfWidth: specifies the maximum distnace in energy away from
  //                   'energy' the candidate can be.  If a value <= 0.0 is
  //                   passed in, then an infinite range is assumed, and the
  //                   the actual closest in energy of photopeak is found.
  //  xraysOnly: If true, reulting photopeak will only be able to have a x-ray
  //             as a source.  If false, can be xray, gamma, annih, S.E., D.E..
  //  sourceGammaType: in the case of gammas produced by positron
  //                   annihilation, this will be set to AnnihilationGamma.
  //                   If so, then transition and transition_index will be set
  //                   only if there was exactly one positron in the decay
  //                   series, otherwise these will be set to null and zero
  //                   respectively.
  //  The relative branching ratio must be at least 1.0e-10 in order for a match
  //  to be successful.
  static void findNearestPhotopeak( const SandiaDecay::Nuclide *nuclide,
                                   const double energy,
                                   const double windowHalfWidth,
                                   const bool xraysOnly,
                                   const SandiaDecay::Transition *&transition,
                                   size_t &transition_index,
                                   SourceGammaType &sourceGammaType );
  
  //findNearestXray(...): returns NULL on error
  static const SandiaDecay::EnergyIntensityPair *findNearestXray(
                              const SandiaDecay::Element *el, double energy );
  
  //toXml(...) serializes the peak and its continuums to the XML document, under
  //  the parent node.  The continuum is only serialized if it is not contained
  //  in the 'continuums' map, otherwise the peaks is assigned the the continuum
  //  ID indicated by 'continuums'.  Note that 'continuums' is not thread safe.
  //A new node named 'Peak' and potentially one named 'Continuum' will be
  //  created under parent.  The node corresponding to "Peak" will be returned.
  rapidxml::xml_node<char> *toXml( rapidxml::xml_node<char> *peak_parent,
              rapidxml::xml_node<char> *continuum_parent,
              std::map<std::shared_ptr<PeakContinuum>,int> &continuums ) const;
  
  //fromXml(...) reverses toXml, with the exception of if the continuum ID isnt
  //  in continuums, an exeption will be thrown.
  //An exception will be thrown if errors are ran into
  void fromXml( const rapidxml::xml_node<char> *node,
            const std::map<int,std::shared_ptr<PeakContinuum> > &continuums );
  

#if( SpecUtils_ENABLE_D3_CHART )
  static std::string gaus_peaks_to_json( const std::vector<std::shared_ptr<const PeakDef> > &peaks,
                                  const std::shared_ptr<const SpecUtils::Measurement> &foreground );
  static std::string peak_json( const std::vector<std::shared_ptr<const PeakDef> > &inpeaks,
                                const std::shared_ptr<const SpecUtils::Measurement> &foreground );
#endif
  

  friend std::ostream &operator<<( std::ostream &stream, const PeakDef &peak );

#if( PERFORM_DEVELOPER_CHECKS )
  //equalEnough(...): tests whether the passed in PeakDef objects are
  //  equal, for most intents and purposes.  Allows some small numerical
  //  rounding to occur.
  //Throws an std::exception with a brief explanaition when an issue is found.
  static void equalEnough( const PeakDef &lhs, const PeakDef &rhs );
#endif
  
public:
  std::string m_userLabel;  //Encoded as UTF8
  
  DefintionType m_type;
  SkewType m_skewType;
  double m_coefficients[NumCoefficientTypes];
  double m_uncertainties[NumCoefficientTypes];
  bool m_fitFor[NumCoefficientTypes];
  
  std::shared_ptr<PeakContinuum> m_continuum;
  
  /** For versioning the XML we will use a major.minor notation.
   Changes made to the XML that are backward compatible (e.g., will read into older version of InterSpec just fine - so this is like if you
   just add a field), increment the minor version.  Breaking changes (e.x., a change to the name of a XML tag) increment the major version
   and will cause older versions of InterSpec to not read in the peak.
   */
  static const int sm_xmlSerializationMajorVersion;
  static const int sm_xmlSerializationMinorVersion;
  
  //A maximum of one of the following will be valid: m_parentNuclide,
  //  m_xrayElement, or m_reaction
  //
  //m_parentNuclide may not be the actual nuclide the decay came from, but
  // instead the _ultimate_ parent the user assigned to this peak
  const SandiaDecay::Nuclide *m_parentNuclide;   //
  const SandiaDecay::Transition *m_transition;   // The general decay
  int m_radparticleIndex; //This specific peak, an object in m_transition->products
  SourceGammaType m_sourceGammaType;
  
  std::vector< CandidateNuclide > m_candidateNuclides;
  
  const SandiaDecay::Element *m_xrayElement;
  float m_xrayEnergy;
  const ReactionGamma::Reaction *m_reaction;
  float m_reactionEnergy;
  
  bool m_useForEnergyCal;
  bool m_useForShieldingSourceFit;
  bool m_useForManualRelEff;
  
  
  /** Fif this peak should be used in the detector response function model fit for intrinsic efficiency */
  bool m_useForDrfIntrinsicEffFit;
  bool m_useForDrfFwhmFit;
  bool m_useForDrfDepthOfInteractionFit;
  
  static const bool sm_defaultUseForDrfIntrinsicEffFit;
  static const bool sm_defaultUseForDrfFwhmFit;
  static const bool sm_defaultUseForDrfDepthOfInteractionFit;
  
  /** Line color of the peak.  Will also (currently) be used to set the fill of
      the peak, just with the alpha channel lowered.
      Alpha not currently allowed, partially because of Wt bug parsing strings
      with alpha specified.
   */
  Wt::WColor m_lineColor;
};//struct PeakDef


std::ostream &operator<<( std::ostream &stream, const PeakDef &peak );


double PeakDef::mean() const
{
//  assert( m_type==PeakDef::GaussianDefined );
  return m_coefficients[PeakDef::Mean];
}

double PeakDef::sigma() const
{
  if( m_type != PeakDef::GaussianDefined )
    throw std::runtime_error( "PeakDef::sigma(): not gaus peak" );
  return m_coefficients[PeakDef::Sigma];
}

double PeakDef::fwhm() const
{
  if( m_type != PeakDef::GaussianDefined )
    throw std::runtime_error( "PeakDef::fwhm(): not gaus peak" );
  return 2.35482*m_coefficients[PeakDef::Sigma];
}

double PeakDef::amplitude() const
{
  if( m_type != PeakDef::GaussianDefined )
    throw std::runtime_error( "PeakDef::amplitude(): not gaus peak" );
  return m_coefficients[PeakDef::GaussAmplitude];
}

double PeakDef::chi2dof() const
{
  return m_coefficients[PeakDef::Chi2DOF];
}

const std::string &PeakDef::userLabel() const
{
  return m_userLabel;
}

void PeakDef::setUserLabel( const std::string &utf8Label )
{
  m_userLabel = utf8Label;
}

bool PeakDef::chi2Defined() const
{
  return (m_coefficients[PeakDef::Chi2DOF] > 0.0);
}


bool PeakDef::gausPeak() const
{
  return (m_type==PeakDef::GaussianDefined);
}

double PeakDef::roiWidth() const
{
  return (upperX() - lowerX());
}

double PeakDef::coefficient( CoefficientType type ) const
{
  return m_coefficients[type];
}

bool PeakDef::fitFor( CoefficientType type ) const
{
  return m_fitFor[type];
}

void PeakDef::setFitFor( CoefficientType type, bool fit )
{
  m_fitFor[type] = fit;
}

const bool *PeakDef::fitFors() const
{
  return m_fitFor;
}

const double *PeakDef::coefficients() const
{
  return m_coefficients;
}

const double *PeakDef::uncertainties() const
{
  return m_uncertainties;
}

void PeakDef::set_coefficient( double val, CoefficientType type )
{
  m_coefficients[type] = val;
}

double PeakDef::uncertainty( CoefficientType type ) const
{
  return m_uncertainties[type];
}

void PeakDef::set_uncertainty( double val, CoefficientType type )
{
  m_uncertainties[type] = val;
}

PeakDef::DefintionType PeakDef::type() const
{
  return m_type;
}

PeakDef::SkewType PeakDef::skewType() const
{
  return m_skewType;
}

void PeakDef::setSkewType( PeakDef::SkewType t )
{
  m_skewType = t;
}

double PeakDef::meanUncert() const
{
  return m_uncertainties[PeakDef::Mean];
}

double PeakDef::sigmaUncert() const
{
  return m_uncertainties[PeakDef::Sigma];
}

double PeakDef::amplitudeUncert() const
{
  return m_uncertainties[PeakDef::GaussAmplitude];
}


void PeakDef::setMean( const double m )
{
  assert( !IsInf(m) && !IsNan(m) );
  m_coefficients[PeakDef::Mean] = m;
}

void PeakDef::setSigma( const double s )
{
  assert( !IsInf(s) && !IsNan(s) );
  
  if( m_type != PeakDef::GaussianDefined )
    throw std::runtime_error( "PeakDef::setSigma(): not gaus peak" );
  m_coefficients[PeakDef::Sigma] = s;
}

void PeakDef::setAmplitude( const double a )
{
  assert( !IsInf(a) && !IsNan(a) );
  
  if( m_type != PeakDef::GaussianDefined )
    throw std::runtime_error( "PeakDef::setAmplitude(): not gaus peak" );
  m_coefficients[PeakDef::GaussAmplitude] = a;
}

void PeakDef::setMeanUncert( const double m )
{
  m_uncertainties[PeakDef::Mean] = m;
}

void PeakDef::setSigmaUncert( const double s )
{
  m_uncertainties[PeakDef::Sigma] = s;
}

void PeakDef::setAmplitudeUncert( const double a )
{
  m_uncertainties[PeakDef::GaussAmplitude] = a;
}


bool PeakDef::useForEnergyCalibration() const
{
  return ( m_useForEnergyCal &&
           ( (m_radparticleIndex>=0 && m_transition && m_parentNuclide)
            || (m_parentNuclide && (m_sourceGammaType==PeakDef::AnnihilationGamma))
            || (m_xrayElement && m_xrayEnergy>0.0)
            || (m_reaction && m_reactionEnergy>0.0)
           ) );
}//void useForEnergyCalibration() const


void PeakDef::useForEnergyCalibration( const bool use )
{
  m_useForEnergyCal = use;
}//void useForEnergyCalibration( const bool use )


bool PeakDef::useForShieldingSourceFit() const
{
  if( !m_useForShieldingSourceFit )
    return m_useForShieldingSourceFit;
  
  switch( m_sourceGammaType )
  {
    case PeakDef::NormalGamma:
    case PeakDef::XrayGamma:
      return (m_radparticleIndex>=0 && m_transition && m_parentNuclide);
    case PeakDef::AnnihilationGamma:
      return (m_parentNuclide && (m_sourceGammaType==PeakDef::AnnihilationGamma));
    case PeakDef::SingleEscapeGamma:
    case PeakDef::DoubleEscapeGamma:
      break;
  }//switch( srcType )
  
  return false;
}//bool useForShieldingSourceFit() const


void PeakDef::useForShieldingSourceFit( const bool use )
{
  m_useForShieldingSourceFit = use;
}//void useForShieldingSourceFit( const bool use )


bool PeakDef::useForManualRelEff() const
{
  if( !m_useForManualRelEff || (!m_parentNuclide && !m_reaction) )
    return false;
  
  switch( m_sourceGammaType )
  {
    case PeakDef::NormalGamma:
      return m_useForManualRelEff;
      
    case PeakDef::XrayGamma:
    case PeakDef::AnnihilationGamma:
    case PeakDef::SingleEscapeGamma:
    case PeakDef::DoubleEscapeGamma:
      break;
  }//switch( srcType )
  
  return false;
}//bool useForManualRelEff() const


void PeakDef::useForManualRelEff( const bool use )
{
  m_useForManualRelEff = use;
}//void useForManualRelEff( const bool use )


bool PeakDef::useForDrfIntrinsicEffFit() const
{
  //if( !m_parentNuclide || !m_transition || (m_sourceGammaType != SourceGammaType::NormalGamma) )
  //  return false;
  
  return m_useForDrfIntrinsicEffFit;
}


void PeakDef::setUseForDrfIntrinsicEffFit( const bool use )
{
  m_useForDrfIntrinsicEffFit = use;
}


bool PeakDef::useForDrfFwhmFit() const
{
  return m_useForDrfFwhmFit;
}


void PeakDef::setUseForDrfFwhmFit( const bool use )
{
  m_useForDrfFwhmFit = use;
}


bool PeakDef::useForDrfDepthOfInteractionFit() const
{
  return m_useForDrfDepthOfInteractionFit;
}


void PeakDef::setUseForDrfDepthOfInteractionFit( const bool use )
{
  m_useForDrfDepthOfInteractionFit = use;
}

#endif
