#ifndef QLPeakDef_h
#define QLPeakDef_h

/* QLPeakDef is a version of InterSpecs PeakDef class that was
 branched 20171227 to create a version slightly more acceptable for QuickLook
 functionality.
 */

#include <memory>
#include <vector>
#include <utility>
#include <stdexcept>


//Forward declaration
class QLPeakDef;
class Measurement;

namespace rapidxml
{
  template<class Ch> class xml_node;
  template<class Ch> class xml_document;
}//namespace rapidxml


struct QLPeakContinuum
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
  
  QLPeakContinuum();
  
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
  bool operator==( const QLPeakContinuum &rhs ) const;
 
  //toXml: XXX = does not support serializing of m_externalContinuum is pretty 
  //  bad currently!  Currently serializes a new continuum for each peak, even
  //  if they should be shared across peaks
  void toXml( rapidxml::xml_node<char> *parent, const int contId ) const;
  void fromXml( const rapidxml::xml_node<char> *node, int &contId );


  
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
  //  XXX - need to verify when writing a QLSpecMeas object to a native file, if
  //        multiple peaks share a global continuum, it is only written once in
  //        the file.
  std::shared_ptr<const Measurement> m_externalContinuum;
  
  static const int sm_xmlSerializationVersion;
};//struct QLPeakContinuum



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



class QLPeakDef
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
  
  static const char *to_string( const CoefficientType type );
  
public:
  QLPeakDef();

  //The following constructor constructs a gaussian peak with no skew or
  //  self defined background/continuum.
  QLPeakDef( double m, double s, double a );

  //The following constructor make a peak were m_PeakDef==false and
  //  instead the peak will be drawn as the difference between data and
  //  background.  If a background histogram is not proivided than a
  //  strait line spanning lowerx and upperx will be used.
  //This usage of QLPeakDef has not been well tested
  QLPeakDef( double lowerx, double upperx, double mean,
            std::shared_ptr<const Measurement> data, std::shared_ptr<const Measurement> background );

//  QLPeakDef( const QLPeakDef &rhs );
//  const QLPeakDef &operator=( const QLPeakDef &rhs );

  //continuum(): access the continuum
  std::shared_ptr<QLPeakContinuum> continuum();
  std::shared_ptr<const QLPeakContinuum> continuum() const;
  
  //getContinuum(): garunteed to be a valid pointer
  std::shared_ptr<QLPeakContinuum> getContinuum();
  
  //setContinuum(...): input must be a valid pointer or exception will be thrown
  void setContinuum( std::shared_ptr<QLPeakContinuum> continuum );
  
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
  
  std::string parentNuclide() const;
  float gammaParticleEnergy() const;
  SourceGammaType sourceGammaType() const;
  
  std::string xrayElement() const;
  float xrayEnergy() const;
  
  std::string reaction() const;
  float reactionEnergy() const;

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

  static bool lessThanByMean( const QLPeakDef &lhs, const QLPeakDef &rhs );
  static bool lessThanByMeanShrdPtr( const std::shared_ptr<const QLPeakDef> &lhs,
                                     const std::shared_ptr<const QLPeakDef> &rhs );
  
  bool operator==( const QLPeakDef &rhs ) const;


  //gaus_integral(): Calculates the area of a Gaussian with specified mean,
  //  sigma, and amplitude, between x0 and x1.
  //  Results have approximately 9 decimal digits of accuracy.
  static double gaus_integral( const double peak_mean,
                               const double peak_sigma,
                               const double peak_amplitude,
                               const double x0, const double x1 );
 
  //toXml(...) serializes the peak and its continuums to the XML document, under
  //  the parent node.  The continuum is only serialized if it is not contained
  //  in the 'continuums' map, otherwise the peaks is assigned the the continuum
  //  ID indicated by 'continuums'.  Note that 'continuums' is not thread safe.
  //A new node named 'Peak' and potentially one named 'Continuum' will be
  //  created under parent.  The node corresponding to "Peak" will be returned.
  rapidxml::xml_node<char> *toXml( rapidxml::xml_node<char> *peak_parent,
              rapidxml::xml_node<char> *continuum_parent,
              std::map<std::shared_ptr<QLPeakContinuum>,int> &continuums ) const;
  
  //fromXml(...) reverses toXml, with the exception of if the continuum ID isnt
  //  in continuums, an exeption will be thrown.
  //An exception will be thrown if errors are ran into
  void fromXml( const rapidxml::xml_node<char> *node,
            const std::map<int,std::shared_ptr<QLPeakContinuum> > &continuums );
  
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
  
  std::shared_ptr<QLPeakContinuum> m_continuum;
  
  static const int sm_xmlSerializationVersion;
  
  //A maximum of one of the following will be valid: m_parentNuclide,
  //  m_xrayElement, or m_reaction
  //
  //m_parentNuclide may not be the actual nuclide the decay came from, but
  // instead the _ultimate_ parent the user assigned to this peak
  std::string m_parentNuclide;   //
  float m_gammaEnergy;
  SourceGammaType m_sourceGammaType;
  
  std::string m_xrayElement;
  float m_xrayEnergy;
  std::string m_reaction;
  float m_reactionEnergy;
  
  bool m_useForCalibration;
  bool m_useForShieldingSourceFit;
};//struct QLPeakDef



double QLPeakDef::mean() const
{
//  assert( m_type==QLPeakDef::GaussianDefined );
  return m_coefficients[QLPeakDef::Mean];
}

double QLPeakDef::sigma() const
{
  if( m_type != QLPeakDef::GaussianDefined )
    throw std::runtime_error( "QLPeakDef::sigma(): not gaus peak" );
  return m_coefficients[QLPeakDef::Sigma];
}

double QLPeakDef::fwhm() const
{
  if( m_type != QLPeakDef::GaussianDefined )
    throw std::runtime_error( "QLPeakDef::fwhm(): not gaus peak" );
  return 2.35482*m_coefficients[QLPeakDef::Sigma];
}

double QLPeakDef::amplitude() const
{
  if( m_type != QLPeakDef::GaussianDefined )
    throw std::runtime_error( "QLPeakDef::amplitude(): not gaus peak" );
  return m_coefficients[QLPeakDef::GaussAmplitude];
}

double QLPeakDef::chi2dof() const
{
  return m_coefficients[QLPeakDef::Chi2DOF];
}

const std::string &QLPeakDef::userLabel() const
{
  return m_userLabel;
}

void QLPeakDef::setUserLabel( const std::string &utf8Label )
{
  m_userLabel = utf8Label;
}

bool QLPeakDef::chi2Defined() const
{
  return (m_coefficients[QLPeakDef::Chi2DOF] > 0.0);
}


bool QLPeakDef::gausPeak() const
{
  return (m_type==QLPeakDef::GaussianDefined);
}

double QLPeakDef::roiWidth() const
{
  return (upperX() - lowerX());
}

double QLPeakDef::coefficient( CoefficientType type ) const
{
  return m_coefficients[type];
}

bool QLPeakDef::fitFor( CoefficientType type ) const
{
  return m_fitFor[type];
}

void QLPeakDef::setFitFor( CoefficientType type, bool fit )
{
  m_fitFor[type] = fit;
}

const bool *QLPeakDef::fitFors() const
{
  return m_fitFor;
}

const double *QLPeakDef::coefficients() const
{
  return m_coefficients;
}

const double *QLPeakDef::uncertainties() const
{
  return m_uncertainties;
}

void QLPeakDef::set_coefficient( double val, CoefficientType type )
{
  m_coefficients[type] = val;
}

double QLPeakDef::uncertainty( CoefficientType type ) const
{
  return m_uncertainties[type];
}

void QLPeakDef::set_uncertainty( double val, CoefficientType type )
{
  m_uncertainties[type] = val;
}

QLPeakDef::DefintionType QLPeakDef::type() const
{
  return m_type;
}

QLPeakDef::SkewType QLPeakDef::skewType() const
{
  return m_skewType;
}

void QLPeakDef::setSkewType( QLPeakDef::SkewType t )
{
  m_skewType = t;
}

double QLPeakDef::meanUncert() const
{
  return m_uncertainties[QLPeakDef::Mean];
}

double QLPeakDef::sigmaUncert() const
{
  return m_uncertainties[QLPeakDef::Sigma];
}

double QLPeakDef::amplitudeUncert() const
{
  return m_uncertainties[QLPeakDef::GaussAmplitude];
}


void QLPeakDef::setMean( const double m )
{
  m_coefficients[QLPeakDef::Mean] = m;
}

void QLPeakDef::setSigma( const double s )
{
  if( m_type != QLPeakDef::GaussianDefined )
    throw std::runtime_error( "QLPeakDef::setSigma(): not gaus peak" );
  m_coefficients[QLPeakDef::Sigma] = s;
}

void QLPeakDef::setAmplitude( const double a )
{
  if( m_type != QLPeakDef::GaussianDefined )
    throw std::runtime_error( "QLPeakDef::setAmplitude(): not gaus peak" );
  m_coefficients[QLPeakDef::GaussAmplitude] = a;
}

void QLPeakDef::setMeanUncert( const double m )
{
  m_uncertainties[QLPeakDef::Mean] = m;
}

void QLPeakDef::setSigmaUncert( const double s )
{
  m_uncertainties[QLPeakDef::Sigma] = s;
}

void QLPeakDef::setAmplitudeUncert( const double a )
{
  m_uncertainties[QLPeakDef::GaussAmplitude] = a;
}


bool QLPeakDef::useForCalibration() const
{
  return ( m_useForCalibration &&
           ( (m_parentNuclide.size() && m_gammaEnergy>0.0)
            || (m_parentNuclide.size() && (m_sourceGammaType==QLPeakDef::AnnihilationGamma))
            || (m_xrayElement.size() && m_xrayEnergy>0.0)
            || (m_reaction.size() && m_reactionEnergy>0.0)
           ) );
}//void useForCalibration() const


void QLPeakDef::useForCalibration( const bool use )
{
  m_useForCalibration = use;
}//void useForCalibration( const bool use )


bool QLPeakDef::useForShieldingSourceFit() const
{
  if( !m_useForShieldingSourceFit )
    return m_useForShieldingSourceFit;
  
  switch( m_sourceGammaType )
  {
    case QLPeakDef::NormalGamma:
    case QLPeakDef::XrayGamma:
      return (m_gammaEnergy>0.0 && m_parentNuclide.size());
    case QLPeakDef::AnnihilationGamma:
      return (m_parentNuclide.size() && (m_sourceGammaType==QLPeakDef::AnnihilationGamma));
    case QLPeakDef::SingleEscapeGamma:
    case QLPeakDef::DoubleEscapeGamma:
      break;
  }//switch( srcType )
  
  return false;
}//bool useForShieldingSourceFit() const


void QLPeakDef::useForShieldingSourceFit( const bool use )
{
  m_useForShieldingSourceFit = use;
}//void useForShieldingSourceFit( const bool use )

#endif
