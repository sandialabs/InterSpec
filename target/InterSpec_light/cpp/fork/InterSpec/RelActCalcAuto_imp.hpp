#ifndef RelActCalcAuto_imp_h
#define RelActCalcAuto_imp_h

#include <utility>
#include <vector>

#include <thread>


#include "Eigen/Dense"

#include "InterSpec/PeakDists.h"
#include "InterSpec/PeakDists_imp.hpp"

#include "SpecUtils/SpecUtilsAsync.h"

#include "InterSpec/PeakDists_imp.hpp"  //for `check_jet_for_NaN(...)`


namespace RelActCalcAutoImp
{
/** Maps an unconstrained trial FWHM to a strictly positive, resolvable width.

 The mapping is exactly the identity above twice the minimum width.  Below that join it is a
 monotone C1 exponential continuation with a positive asymptote, so invalid optimizer trials never
 make the peak model throw and normal persisted FWHM semantics remain unchanged. */
template<typename T>
T resolvable_fwhm_continuation( const T &raw_fwhm, const double channel_width,
                                const double minimum_channels = 1.25 )
{
  assert( std::isfinite(channel_width) && (channel_width > 0.0) );
  assert( std::isfinite(minimum_channels) && (minimum_channels > 0.0) );

  const double lower_asymptote = minimum_channels * channel_width;
  const double join = 2.0 * lower_asymptote;
  double scalar = 0.0;
  if constexpr ( std::is_same_v<T,double> )
    scalar = raw_fwhm;
  else
    scalar = raw_fwhm.a;

  if( scalar >= join )
    return raw_fwhm;

  const double scale = join - lower_asymptote;
  return T(lower_asymptote) + T(scale)*exp((raw_fwhm - T(join))/T(scale));
}
}//namespace RelActCalcAutoImp


namespace RelActCalcAuto
{
/** A stand-in for the `PeakDef` class to allow auto-differentiation, and also simplify things  */
template<typename T>
struct PeakDefImp
{
  T m_mean = T(0.0);
  T m_sigma = T(0.0);
  T m_amplitude = T(0.0);
  T m_skew_pars[6] = { T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0) };

  /** True energy of the source gamma or x-ray. */
  double m_src_energy = 0.0;

  PeakDef::SkewType m_skew_type = PeakDef::SkewType::NoSkew;
  PeakDef::SourceGammaType m_gamma_type = PeakDef::SourceGammaType::NormalGamma;

  size_t m_rel_eff_index = std::numeric_limits<size_t>::max();

  // Some functions to be compatible with PeakDef, in templated functions
  const T &mean() const { return m_mean; }
  const T &sigma() const { return m_sigma; }
  const T &amplitude() const { return m_amplitude; }
  PeakDef::SkewType skewType() const { return m_skew_type; }
  const T *skew_parameters() const { return m_skew_pars; }

  void setMean( const T &mean ) { m_mean = mean; }
  void setSigma( const T &sigma ) { m_sigma = sigma; }
  void setAmplitude( const T &amp ) { m_amplitude = amp; }

  inline void setSkewType( const PeakDef::SkewType &type )
  {
    m_skew_type = type;
  }

  inline void set_coefficient( T val, const PeakDef::CoefficientType &coef )
  {
    const int index = static_cast<int>( coef - PeakDef::CoefficientType::SkewPar0 );
    if( index >= (sizeof(m_skew_pars) / sizeof(m_skew_pars[0])) )
    {
      assert( coef == PeakDef::CoefficientType::Chi2DOF );
      return; //Like Chi2Dof
    }

    assert( (index >= 0) && (index < 6) );
    m_skew_pars[index] = val;
  }

  inline void set_uncertainty( T val, const PeakDef::CoefficientType &coef )
  {
    // Dummy so `PeakFitDiffCostFunction::parametersToPeaks(...)` can set uncertainties and compile.
  }
  
  void gauss_integral( const float *energies, T *channels, const size_t nchannel ) const
  {
    check_jet_array_for_NaN( channels, nchannel );

    switch( m_skew_type )
    {
      case PeakDef::NumSkewType:
        assert( 0 );
        //throw runtime_error( "RelActCalcAuto gauss_integral: NumSkewType is not a valid skew type" );
        // Fall through to NoSkew for non-debug builds

      case PeakDef::SkewType::NoSkew:
        PeakDists::gaussian_integral( m_mean, m_sigma, m_amplitude, energies, channels, nchannel );
        break;

      case PeakDef::SkewType::Bortel:
        PeakDists::bortel_integral( m_mean, m_sigma, m_amplitude, m_skew_pars[0], energies, channels, nchannel );
        break;

      case PeakDef::SkewType::CrystalBall:
        PeakDists::crystal_ball_integral( m_mean, m_sigma, m_amplitude, m_skew_pars[0], m_skew_pars[1], energies, channels, nchannel );
        break;

      case PeakDef::SkewType::DoubleSidedCrystalBall:
        PeakDists::double_sided_crystal_ball_integral( m_mean, m_sigma, m_amplitude,
                                           m_skew_pars[0], m_skew_pars[1],
                                           m_skew_pars[2], m_skew_pars[3],
                                           energies, channels, nchannel );
        break;

      case PeakDef::SkewType::GaussExp:
        PeakDists::gauss_exp_integral( m_mean, m_sigma, m_amplitude, m_skew_pars[0], energies, channels, nchannel );
        break;

      case PeakDef::SkewType::ExpGaussExp:
        PeakDists::exp_gauss_exp_integral( m_mean, m_sigma, m_amplitude, m_skew_pars[0], m_skew_pars[1], energies, channels, nchannel );
        break;

      case PeakDef::SkewType::VoigtPlusBortel:
        PeakDists::voigt_exp_integral( m_mean, m_sigma, m_amplitude, m_skew_pars[0], m_skew_pars[1], m_skew_pars[2], energies, channels, nchannel );
        break;

      case PeakDef::SkewType::GaussPlusBortel:
        PeakDists::gauss_plus_bortel_integral( m_mean, m_sigma, m_amplitude, m_skew_pars[0], m_skew_pars[1], energies, channels, nchannel );
        break;

      case PeakDef::SkewType::DoubleBortel:
        PeakDists::double_bortel_integral( m_mean, m_sigma, m_amplitude, m_skew_pars[0], m_skew_pars[1], m_skew_pars[2], energies, channels, nchannel );
        break;

      case PeakDef::SkewType::GadrasGeneric:
        PeakDists::gadras_integral( m_mean, m_sigma, m_amplitude, m_skew_pars,
                                    PeakDists::GadrasMaterial::Generic, energies, channels, nchannel );
        break;

      case PeakDef::SkewType::GadrasCZT:
        PeakDists::gadras_integral( m_mean, m_sigma, m_amplitude, m_skew_pars,
                                    PeakDists::GadrasMaterial::CZT_CdTe, energies, channels, nchannel );
        break;
    }//switch( skew_type )

    check_jet_array_for_NaN( channels, nchannel );
  }//void gauss_integral( const float *energies, double *channels, const size_t nchannel ) const

};//struct PeakDefImp

template<typename T>
struct PeakContinuumImp
{
  PeakContinuum::OffsetType m_type = PeakContinuum::OffsetType::NoOffset;

  T m_lower_energy = T(0.0);
  T m_upper_energy = T(0.0);
  T m_reference_energy = T(0.0);
  std::array<T,5> m_values = { T(0.0), T(0.0), T(0.0), T(0.0), T(0.0) };

  // Some functions for compatibility with `PeakContinuum` in PeakDef.h.
  void setRange( const T lowerenergy, const T upperenergy )
  {
    m_lower_energy = lowerenergy;
    m_upper_energy = upperenergy;
  }

  T lowerEnergy() const { return m_lower_energy; }
  T upperEnergy() const { return m_upper_energy; }

  void setType( PeakContinuum::OffsetType type )
  {
    m_type = type;
  }

  PeakContinuum::OffsetType type() const
  {
    return m_type;
  }

  T referenceEnergy() const { return m_reference_energy; }
  const std::array<T,5> &parameters() const { return m_values; }

  void setParameters( T referenceEnergy, const T *parameters, const T *uncertainties [[maybe_unused]] )
  {
    m_reference_energy = referenceEnergy;
    const size_t npar = PeakContinuum::num_parameters(m_type);
    for( size_t i = 0; i < npar; ++i )
      m_values[i] = parameters[i];
  }

  std::shared_ptr<const SpecUtils::Measurement> externalContinuum() const { return nullptr; }
};//struct PeakContinuumImp


template<typename T>
struct PeaksForEnergyRangeImp
{
  std::vector<PeakDefImp<T>> peaks;

  PeakContinuumImp<T> continuum;

  size_t first_channel;
  size_t last_channel;
  bool no_gammas_in_range;
  bool forced_full_range;

  /** Peak plus continuum counts for [first_channel, last_channel] */
  std::vector<T> peak_counts;
};//struct PeaksForEnergyRangeImp


}//namespace RelActCalcAuto

#endif //RelActCalcAuto_imp_h
