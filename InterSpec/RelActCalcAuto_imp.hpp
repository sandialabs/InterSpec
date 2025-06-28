#ifndef RelActCalcAuto_imp_h
#define RelActCalcAuto_imp_h

#include <vector>

#include <thread>

#include <boost/math/distributions/normal.hpp>

#include "Eigen/Dense"

#include "InterSpec/PeakDists.h"
#include "InterSpec/PeakDists_imp.hpp"

#include "SpecUtils/SpecUtilsAsync.h"

#include "InterSpec/PeakDists_imp.hpp"  //for `check_jet_for_NaN(...)`


namespace RelActCalcAuto
{
/** A stand-in for the `PeakDef` class to allow auto-differentiation, and also simplify things  */
template<typename T>
struct PeakDefImp
{
  T m_mean = T(0.0);
  T m_sigma = T(0.0);
  T m_amplitude = T(0.0);
  T m_skew_pars[4] = { T(0.0), T(0.0), T(0.0), T(0.0) };

  const SandiaDecay::Nuclide *m_parent_nuclide = nullptr;
  const SandiaDecay::Transition *m_transition = nullptr;
  size_t m_rad_particle_index = 0;

  const SandiaDecay::Element *m_xray_element = nullptr;
  const ReactionGamma::Reaction *m_reaction = nullptr;

  /** True energy of the source gamma or x-ray. */
  double m_src_energy = 0.0;

  PeakDef::SkewType m_skew_type = PeakDef::SkewType::NoSkew;
  PeakDef::SourceGammaType m_gamma_type = PeakDef::SourceGammaType::NormalGamma;

  size_t m_rel_eff_index = std::numeric_limits<size_t>::max();

  // Some functions to be compatible with PeakDef, in templated functions
  const T &mean() const { return m_mean; }
  const T &sigma() const { return m_sigma; }
  const T &amplitude() const { return m_amplitude; }

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

    assert( (index >= 0) && (index < 4) );
    m_skew_pars[index] = val;
  }

  void gauss_integral( const float *energies, T *channels, const size_t nchannel ) const
  {
    check_jet_array_for_NaN( channels, nchannel );

    switch( m_skew_type )
    {
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
    }//switch( skew_type )

    check_jet_array_for_NaN( channels, nchannel );
  }//void gauss_integral( const float *energies, double *channels, const size_t nchannel ) const


  /** Gives lower and upper energies that contain `1.0 - missing_frac` of the peak,

   For skewed distributions, particularly Crystal Bal, the energy range given by the desired
   peak coverage can get huge, so you can also specify the max number of FWHM to use, which will
   truncate the range in these cases of very large skew.

   @param missing_frac The fraction of the peaks area you are okay not including; half this
          amount will be missing from both lower and upper distribution tails.
          A typical value might be like 1.0E-4 to get 99.99% of the peak.
   @param max_num_fwhm If greater than zero, the energy range will be truncated to be this
          amount of FWHMs from the mean, if the coverage would have the energy range go really
          far out.
   */
  std::pair<double,double> peak_coverage_limits( const double missing_frac, const double max_num_fwhm )
  {
    using namespace std;

    double skew_pars[4] = { 0.0 };

    double mean, sigma;
    if constexpr ( !std::is_same_v<T, double> )
    {
      mean = m_mean.a;
      sigma = m_sigma.a;
    }else
    {
      mean = m_mean;
      sigma = m_sigma;
    }

    const size_t nskew_par = PeakDef::num_skew_parameters(m_skew_type);
    for( size_t i = 0; i < nskew_par; ++i )
    {
      if constexpr ( !std::is_same_v<T, double> )
        skew_pars[i] = m_skew_pars[i].a;
      else
        skew_pars[i] = m_skew_pars[i];
    }


    pair<double,double> vis_limits;

    switch( m_skew_type )
    {
      case PeakDef::SkewType::NoSkew:
      {
        const boost::math::normal_distribution gaus_dist( 1.0 );
        vis_limits.first = mean +  sigma*boost::math::quantile( gaus_dist, 0.5*missing_frac );
        vis_limits.second = mean + sigma*boost::math::quantile( gaus_dist, 1.0 - 0.5*missing_frac );
        break;
      }

      case PeakDef::SkewType::Bortel:
        vis_limits = PeakDists::bortel_coverage_limits( mean, sigma, skew_pars[0], missing_frac );
        break;

      case PeakDef::SkewType::GaussExp:
        vis_limits = PeakDists::gauss_exp_coverage_limits( mean, sigma, skew_pars[0], missing_frac );
        break;

      case PeakDef::SkewType::CrystalBall:
        try
      {
        vis_limits = PeakDists::crystal_ball_coverage_limits( mean, sigma, skew_pars[0], skew_pars[1], missing_frac );
      }catch( std::exception & )
      {
        // CB dist can have really long tail, causing the coverage limits to fail, because
        //  of unreasonable values - in this case we'll just go way out
        vis_limits.first = mean - 15.0*sigma;

        const boost::math::normal_distribution gaus_dist( 1.0 );
        vis_limits.second = mean + sigma*boost::math::quantile( gaus_dist, 1.0 - missing_frac );
      }
        break;

      case PeakDef::SkewType::ExpGaussExp:
        vis_limits = PeakDists::exp_gauss_exp_coverage_limits( mean, sigma, skew_pars[0],
                                                              skew_pars[1], missing_frac );
        break;

      case PeakDef::SkewType::DoubleSidedCrystalBall:
        try
      {
        vis_limits = PeakDists::double_sided_crystal_ball_coverage_limits( mean, sigma, skew_pars[0],
                                                                          skew_pars[1], skew_pars[2], skew_pars[3], missing_frac );
      }catch( std::exception & )
      {
        // CB dist can have really long tail, causing the coverage limits to fail, because
        //  of unreasonable values - in this case we'll just go way out
        vis_limits.first  = mean - 15.0*sigma;
        vis_limits.second = mean + 15.0*sigma;
      }//try / catch
        break;
    }//switch( m_skew_type )

    if( max_num_fwhm > 0.0 )
    {
      vis_limits.first = std::max( vis_limits.first, mean - 2.35482*max_num_fwhm*sigma );
      vis_limits.second = std::min( vis_limits.second, mean + 2.35482*max_num_fwhm*sigma );
    }//if( max_num_fwhm > 0.0 )

    return vis_limits;
  }//pair<double,double> peak_coverage_limits( const double missing_frac = 1.0E-6 )

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
