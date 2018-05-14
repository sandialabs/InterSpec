#ifndef PeakDef_h
#define PeakDef_h
/* InterSpec: an application to analyze spectral gamma radiation data.
 
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

#include "InterSpec_config.h"

#include <map>
#include <memory>
#include <vector>
#include <utility>
#include <stdexcept>

#include "InterSpec/ReactionGamma.h"
#include "sandia_decay/SandiaDecay.h"
#include "InterSpec/DecayDataBaseServer.h"


//Forward declaration
class PeakDef;
class Measurement;
namespace SandiaDecay
{
  struct Nuclide;
  struct Element;
  struct Transition;
  struct RadParticle;
}//namespace SandiaDecay


namespace rapidxml
{
  template<class Ch> class xml_node;
  template<class Ch> class xml_document;
}//namespace rapidxml



struct PeakContinuum
{
  enum OffsetType
  {
    NoOffset,
    Constant = 1,  //purposely set to be the size of the expected num paramters
    Linear,
    Quardratic,
    Cubic,
    External
  };//enum OffsetType
  
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

  //setParameters: a convienience function which actually calls other form of
  //  this member function.  This funtion assumes both passed in arrays are
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
  
  double referenceEnergy() const { return m_refernceEnergy; }
  const std::vector<double> &parameters() const { return m_values; }
  const std::vector<double> &unertainties() const { return m_uncertainties; }
  std::vector<bool> fitForParameter() const { return m_fitForValue; }
  std::shared_ptr<const Measurement> externalContinuum() const { return m_externalContinuum; }
  
  //setGlobalContinuum: throws if not a External OffsetType.
  void setExternalContinuum( const std::shared_ptr<const Measurement> &data );
 
  //setRange: sets the energy range this continuum is applicable for
  void setRange( const double lowerenergy, const double upperenergy );
  
  //calc_linear_continuum_eqn: sets this to be a Linear OffsetType continuum,
  //  and the range to be from x0 to x1, with the resulting linear polynomial
  //  relative to m_refernceEnergy == x0
  //  nSideBinToAverage is the number of bins on each side of x0_bin and x1_bin
  //  that should be used to find data y-hieght for that respective side of
  //  continuum
  void calc_linear_continuum_eqn( const std::shared_ptr<const Measurement> &data,
                                  const double x0, const double x1,
                                  const int nSideBinToAverage );

  //offset_integral: returns the area of the continuum from x0 to x1.  If m_type
  //  is NoOffset, then will return 0.  If a polynomial, then the integral of
  //  the contonuum density will be returned.  If globally defined continuum,
  //  then the integral will be returned, using linear interpolation for
  //  ranges not exactly on the bin edges.
  double offset_integral( const double x0, const double x1 ) const;
  
  //defined: returns true if a _valid_ polynomial or external continuum type
  bool defined() const;
  
  //energyRangeDefined: returns if an energy range has explicitely been set
  bool energyRangeDefined() const;
  
  //isPolynomial:
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
                                const std::shared_ptr<const Measurement> &data,
                                const size_t nbinEachSide,
                                double &m, double &y0 );
  
  //offset_eqn_integral: integrates the continuum density, specified by coefs,
  //  to return the continuum area.
  static double offset_eqn_integral( const double *coefs,
                                     OffsetType type,
                                     double x0, double x1,
                                     const double reference_energy );
  
  //translate_offset_polynomial: if you want the exact same polynomial
  //  line, but would like it with reference to a different energy,
  //  use this funtion.
  //XXX - Currently only supports constant and linear polynomials.
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
  
  //m_refernceEnergy: the energy polynomials are evaluated relative to. Ex.
  //  continuum_density = m_values[0] + (energy - m_refernceEnergy)*m_values[1];
  double m_refernceEnergy;
  
  //m_values, m_uncertainties: the polynomial coefficients and their
  //  uncertainties.  The polynomial is actually a density per keV of the
  //  continuum, so needs to be integrated over the width of the bin.  Also,
  //  the energy is relative to m_refernceEnergy. E.g. at the m_refernceEnergy
  //  the continuum will have a density of m_values[0]
  std::vector<double> m_values, m_uncertainties;
  std::vector<bool> m_fitForValue;
  
  //m_externalContinuum: since the global continuum may be changed after fitting
  //  for a peak (thus possibly resulting in a poorly fit peak), we'll keep a
  //  refernce to the original.  This is slightly inefficient, but shouldnt be
  //  to bad.
  //  XXX - need to verify when writing a SpecMeas object to a native file, if
  //        multiple peaks share a global continuum, it is only written once in
  //        the file.
  std::shared_ptr<const Measurement> m_externalContinuum;
  
  static const int sm_xmlSerializationVersion;
  
  friend std::ostream &operator<<( std::ostream &, const PeakContinuum & );
};//struct PeakContinuum



//skewedGaussianIntegral(...):
//A gaussian convoluted with an exponential:
//  http://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
double skewedGaussianIntegral( double xlow, double xhigh,
                               double mean, double sigma,
                               double lambda );
//The chromotagraphy form:
double skewedGaussianIntegral( double x0,  //x-value to start integrating at
                               double x1,  //x-value to stop integrating at
                               double amplitude,  //Area of entire gaussian
                               double peak_center,//
                               double width, //
                               double skewness  //must be >= 0.03*width
                              );


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
size_t findROILimit( const PeakDef &peak, const std::shared_ptr<const Measurement> &data, bool high);


//findROIEnergyLimits(...): a convience function that calls findROILimit(...)
//  inorder to set 'lowerEnengy' and 'upperEnergy'.  IF data is an invalid
//  pointer, then PeakDef::lowerX() and PeakDef::upperX() are used.
void findROIEnergyLimits( double &lowerEnengy, double &upperEnergy,
                         const PeakDef &peak, const std::shared_ptr<const Measurement> &data );

//Returns if area from start2 to end2 is greater than or equal to the
//  area of (start1 to end1 minus nsigma).
//Answers the question: is region 2 similar or greater in area as region 1
bool isStatisticallyGreaterOrEqual( const size_t start1, const size_t end1,
                                    const size_t start2, const size_t end2,
                                    const std::shared_ptr<const Measurement> &dataH,
                                    const double nsigma );

//If a polynomial continuum and the peakdoesnt specify the range, then it
//  will look at the data histogram for when data starts increasing, and set
//  this as the limit.
//Does nothing if the histogram is null.
void estimatePeakFitRange( const PeakDef &peak, const std::shared_ptr<const Measurement> &dataH,
                          int &xlowbin, int &xhighbin );


class PeakDef
{
public:
  enum SkewType
  {
    NoSkew, LandauSkew
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
    LandauAmplitude,   //multiplies peak ampliture (so is between 0.0 and ~0.2)
    LandauMode,
    LandauSigma,
//    OffsetPolynomial0,
//    OffsetPolynomial1,
//    OffsetPolynomial2,
//    OffsetPolynomial3,
//    OffsetPolynomial4,
//    RangeStartEnergy,
//    RangeEndEnergy,
    Chi2DOF,  //for peaks that share a ROI/Continuum, this values is for entire ROI/Continuum
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
  
  static const char *to_string( const CoefficientType type );
  
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
            std::shared_ptr<const Measurement> data, std::shared_ptr<const Measurement> background );

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
  
  //The bellow should in principle take care of gaussian area and the skew area
  double peakArea() const;
  double peakAreaUncert() const;
  
  //setPeakArea(): sets the area of the Gaussian + Skew (not just the gaussian
  //  component)
  void setPeakArea( const double a );
  void setPeakAreaUncert( const double a );
  
  
  inline double meanUncert() const;
  inline double sigmaUncert() const;
  inline double amplitudeUncert() const;


  inline void setMean( const double m );
  inline void setSigma( const double s );
  inline void setAmplitude( const double a );

  inline void setMeanUncert( const double m );
  inline void setSigmaUncert( const double s );
  inline void setAmplitudeUncert( const double a );


  inline bool useForCalibration() const;
  inline void useForCalibration( const bool use );

  inline bool useForShieldingSourceFit() const;
  inline void useForShieldingSourceFit( const bool use );

  inline double chi2dof() const;
  inline bool chi2Defined() const;

  inline const std::string &userLabel() const;
  inline void setUserLabel( const std::string &utf8Label );
  
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
  //  or 1022 keV is subrtacted off of the SandiaDecay::RadParticle::energy
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


  //inheritUserSelectedOptions(): sets m_fitFor[] values, useForCalibration,
  //  useForShieldingSourceFit, and nuclide/xray/reaction values.
  //  If 'inheritNonFitForValues' is specified true, then mean, width, and or
  //  amplitude will also be inherited (only) IF they were specfified to not be
  //  fit for.
  //This function does not modify the continuum values (including ROI extent) or
  //  type of peak.
  void inheritUserSelectedOptions( const PeakDef &parent,
                                  const bool inheritNonFitForValues );

  double lowerX() const;
  double upperX() const;
  inline double roiWidth() const;
  
  //gauss_integral(): gives the area of the Gaussian and Skew (if applicable)
  //  components of the peak, between x0 and x1.
  //  Results have approximately 9 decimal digits of accuracy.
  double gauss_integral( const double x0, const double x1 ) const;
  
  //offset_integral(): gives area of the continuum component between x0 and x1.
  double offset_integral( const double x0, const double x1 ) const;

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
  
  //skew_integral(): gives the difference in aarea between the gaussian, and
  //  gaussian with skew applied, between x0 and x1.
  double skew_integral( const double x0, const double x1 ) const;

  static double landau_potential_lowerX( const double peak_mean, const double peak_sigma );
  static double landau_potential_upperX( const double peak_mean, const double peak_sigma );
  static double landau_integral( const double x0, const double x1,
                                 const double peak_mean,
                                 const double land_amp,
                                 const double land_mode,
                                 const double land_sigma );
  static double skew_integral( const double xlow,
                               const double xhigh,
                               const double peak_amplitude,
                               const double peak_mean,
                               const double p0,
                               const double p1,
                               const double p2 );

  static bool lessThanByMean( const PeakDef &lhs, const PeakDef &rhs );
  static bool lessThanByMeanShrdPtr( const std::shared_ptr<const PeakDef> &lhs,
                                     const std::shared_ptr<const PeakDef> &rhs );
  
  bool operator==( const PeakDef &rhs ) const;


  //gaus_integral(): Calculates the area of a Gaussian with specified mean,
  //  sigma, and amplitude, between x0 and x1.
  //  Results have approximately 9 decimal digits of accuracy.
  static double gaus_integral( const double peak_mean,
                               const double peak_sigma,
                               const double peak_amplitude,
                               const double x0, const double x1 );

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
  
  /** Use hyristics to guess at a reasonable age a nuclide might be if
   * encountered in the real world.  
   *
   * \param nuc Nuclide to return the age for.
   * \param strAnswer If provided, will assign a suggested string for the age
   *        to this string.  May yeild either a actual unit of time (ex. 
   *        "2.31 years") or a mutliple of half lives (ex. "3.5 HL").
   * Is set to have a maximum of 2 decimal points.
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
  
  //The Landau is evaluated as:
  //  landau_amplitude*TMath::Landau(mean-x,landau_mode,landau_sigma,true)
  //  reasonable
  //landau_mode is how far to the left of the peak mean the most probable value
  //  of the landau should be (a larger landau_mode shifts further to left)
  //landau_sigma is the width of the landau distribution
  //landau_amplitude is the total area of the landau distribution if integrated
  //  from negative to positive infinity
  //The maximum value of the landau distribution for this peak is at
  //  mean-landau_mode+(0.22278*landau_sigma)
  //  and by 2.622 landau_sigma to the right of this is down 0.01
  //  times the maximum amplitude.  It takes about 25 landau_sigma to the left
  //  of this point an amplitude of 0.01 times the maximum value
  //  double landau_amplitude, landau_mode, landau_sigma; //landau_amplitude==0.0 indicates dont draw landau
  DefintionType m_type;
  SkewType m_skewType;
  double m_coefficients[NumCoefficientTypes];
  double m_uncertainties[NumCoefficientTypes];
  bool m_fitFor[NumCoefficientTypes];
  
  std::shared_ptr<PeakContinuum> m_continuum;
  
  static const int sm_xmlSerializationVersion;
  
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
  
  bool m_useForCalibration;
  bool m_useForShieldingSourceFit;
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
  m_coefficients[PeakDef::Mean] = m;
}

void PeakDef::setSigma( const double s )
{
  if( m_type != PeakDef::GaussianDefined )
    throw std::runtime_error( "PeakDef::setSigma(): not gaus peak" );
  m_coefficients[PeakDef::Sigma] = s;
}

void PeakDef::setAmplitude( const double a )
{
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


bool PeakDef::useForCalibration() const
{
  return ( m_useForCalibration &&
           ( (m_radparticleIndex>=0 && m_transition && m_parentNuclide)
            || (m_parentNuclide && (m_sourceGammaType==PeakDef::AnnihilationGamma))
            || (m_xrayElement && m_xrayEnergy>0.0)
            || (m_reaction && m_reactionEnergy>0.0)
           ) );
}//void useForCalibration() const


void PeakDef::useForCalibration( const bool use )
{
  m_useForCalibration = use;
}//void useForCalibration( const bool use )


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

#endif
