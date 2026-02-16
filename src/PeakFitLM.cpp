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


#ifdef _MSC_VER
#undef isinf
#undef isnan
#undef isfinite
#undef isnormal
#endif

#include "ceres/ceres.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "SpecUtils/SpecFile.h"
#include "InterSpec/PeakDists.h"
#include "InterSpec/PeakFitLM.h"

#if( PEAK_FIT_LM_PARALLEL_ROIS )
#include <future>
#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>
#include "SpecUtils/SpecUtilsAsync.h"
#endif
#include "InterSpec/PeakFitUtils.h"
#include "SpecUtils/EnergyCalibration.h"
#include "InterSpec/DetectorPeakResponse.h"

#include "InterSpec/PeakFit_imp.hpp"
#include "InterSpec/PeakDists_imp.hpp"
#include "InterSpec/RelActCalcAuto_imp.hpp"

// Undefine isnan and isinf macros
#undef isnan
#undef isinf


using namespace std;

namespace
{
/** To make the code of `PeakFitDiffCostFunction` work with either `PeakDef`
 or `RelActCalcAuto::PeakDefImp<T>`, we'll add in a few interface functions
 */
template<typename T>
struct PeakContinuumImp : public RelActCalcAuto::PeakContinuumImp<T>
{
  std::shared_ptr<const SpecUtils::Measurement> m_external_continuum;
  std::shared_ptr<const SpecUtils::Measurement> externalContinuum() const { return m_external_continuum; }
  void setExternalContinuum( const std::shared_ptr<const SpecUtils::Measurement> &data ){ m_external_continuum = data; }
};

template<typename T>
struct PeakDefImpWithCont : public RelActCalcAuto::PeakDefImp<T>
{
  std::shared_ptr<PeakContinuumImp<T>> m_continuum;
  std::shared_ptr<PeakContinuumImp<T>> getContinuum(){ return m_continuum; }
  void setContinuum( const std::shared_ptr<PeakContinuumImp<T>> &continuum ) { m_continuum = continuum;}
  void setMeanUncert( const T &mean_uncert ){ }
  void setAmplitudeUncert( const T &amp_uncert ){ }
  void setSigmaUncert( const T &sigma_uncert ){ }
  void inheritUserSelectedOptions( const PeakDef &parent, const bool inheritNonFitForValues [[maybe_unused]] )
  {
    //We dont really need this function, since we never use any of this info, but whatever.
    this->m_parent_nuclide = parent.parentNuclide();
    this->m_transition = parent.nuclearTransition();
    this->m_rad_particle_index = parent.decayParticleIndex();
    this->m_xray_element = parent.xrayElement();
    this->m_reaction = parent.reaction();
    this->m_src_energy = (this->m_parent_nuclide || this->m_xray_element || this->m_reaction)
                          ? parent.gammaParticleEnergy() : 0.0f;
    this->m_gamma_type = parent.sourceGammaType();
    //m_rel_eff_index = ...
  }

  static bool lessThanByMean( const PeakDefImpWithCont<T> &lhs, const PeakDefImpWithCont<T> &rhs )
  {
    return (lhs.m_mean < rhs.m_mean);
  }
};//struct PeakDefImpWithCont


/** Replaces contents of input peaks vector with completely new peaks, that crucually have had new continuums allocated. */
void local_unique_copy_continuum( vector<shared_ptr<const PeakDef>> &input_peaks )
{
  map<std::shared_ptr<const PeakContinuum>,vector<shared_ptr<PeakDef>>> contToPeaks;
  for( auto &p : input_peaks )
    contToPeaks[p->continuum()].push_back( make_shared<PeakDef>(*p) );

  for( auto &pp : contToPeaks )
  {
    pp.second[0]->makeUniqueNewContinuum();
    auto newcont = pp.second[0]->continuum();
    for( size_t i = 1; i < pp.second.size(); ++i )
      pp.second[i]->setContinuum( newcont );
  }

  input_peaks.clear();
  for( auto &pp : contToPeaks )
  {
    for( auto p : pp.second )
      input_peaks.push_back( p );
  }
  std::sort( begin(input_peaks), end(input_peaks), &PeakDef::lessThanByMeanShrdPtr );
}//unique_copy_continuum(...)
}//namespace


namespace PeakFitLM
{

/** Info about a single ROI (region of interest), used internally by PeakFitDiffCostFunction.
 A ROI is defined by peaks that share the same PeakContinuum pointer.
*/
struct RoiInfo
{
  std::vector<std::shared_ptr<const PeakDef>> peaks;  // sorted by mean
  double lower_energy;
  double upper_energy;
  size_t lower_channel;
  size_t upper_channel;
  double ref_energy;
  PeakContinuum::OffsetType offset_type;
  bool use_lls_for_cont;
  double max_initial_sigma;   // max sigma across peaks in this ROI
  size_t num_fit_sigmas;      // number of peaks with fitFor(Sigma)==true
};//struct RoiInfo


/** PeakFitDiffCostFunction fits peaks from one or more ROIs simultaneously using Ceres.
 Equivalent of the MultiPeakFitChi2Fcn class, but using Levenberg-Marquardt differentiation.

 When multiple ROIs are present, skew parameters may be energy-dependent (interpolated between
 anchor energies at the spectrum lower/upper bounds), if the ROIs span more than 100 keV and
 the skew type has energy-dependent parameters.

 Parameter layout (default / shared-skew mode):
   [skew_lower_pars (num_skew values) | skew_upper_pars* (M values, only energy-dep params)]
   | ROI_0_cont | ROI_0_sigma | ROI_0_mean | ROI_1_cont | ROI_1_sigma | ROI_1_mean | ...
   * only present when m_fit_skew_energy_dependence == true

 Parameter layout (IndependentSkewValues option):
   ROI_0_cont | ROI_0_sigma | ROI_0_mean | ROI_0_skew (num_skew values)
   | ROI_1_cont | ROI_1_sigma | ROI_1_mean | ROI_1_skew | ...
   There is no shared skew block; each ROI carries its own num_skew parameters at the end of
   its per-ROI block.
*/
struct PeakFitDiffCostFunction
{
  /** A struct with the info to feed into Ceres. */
  struct ProblemSetup
  {
    /** The starting parameters to use. */
    vector<double> m_parameters;
    /** The indexes of paramters that need to be held constant in the fit.  E.g. if a mean is held constant for a peak whose FWHM/amp is being fit. */
    vector<int> m_constant_parameters;
    /** Lower bounds parameters should be allowed to go to; will not have a value if it shouldnt be restricted. Will be same size as `m_parameters`. */
    vector<std::optional<double>> m_lower_bounds;
    /** Upper bounds parameters should be allowed to go to; will not have a value if it shouldnt be restricted. Will be same size as `m_parameters`. */
    vector<std::optional<double>> m_upper_bounds;
  };//struct ProblemSetup


  /** Build the vector of RoiInfo structs from starting peaks and data.
   Groups peaks by shared PeakContinuum pointer, fills channel bounds, sigma, etc.
   Static so it can be called from the constructor initializer list via a lambda.
  */
  static std::vector<RoiInfo> make_rois(
      const std::shared_ptr<const SpecUtils::Measurement> &data,
      const std::vector<std::shared_ptr<const PeakDef>> &starting_peaks,
      const Wt::WFlags<PeakFitLM::PeakFitLMOptions> & /*options*/ )
  {
    // Group peaks by continuum pointer
    std::map<std::shared_ptr<const PeakContinuum>, std::vector<std::shared_ptr<const PeakDef>>> cont_to_peaks;
    for( const auto &p : starting_peaks )
      cont_to_peaks[p->continuum()].push_back( p );

    std::vector<RoiInfo> rois;
    rois.reserve( cont_to_peaks.size() );

    for( auto &kv : cont_to_peaks )
    {
      RoiInfo roi;
      roi.peaks = kv.second;
      std::sort( roi.peaks.begin(), roi.peaks.end(), &PeakDef::lessThanByMeanShrdPtr );

      const std::shared_ptr<const PeakContinuum> &cont = kv.first;
      roi.lower_energy = cont->lowerEnergy();
      roi.upper_energy = cont->upperEnergy();
      roi.lower_channel = data->find_gamma_channel( static_cast<float>( roi.lower_energy ) );
      roi.upper_channel = data->find_gamma_channel( static_cast<float>( roi.upper_energy ) );
      roi.offset_type = cont->type();

      const double prev_ref = cont->referenceEnergy();
      roi.ref_energy = ((prev_ref >= roi.lower_energy) && (prev_ref <= roi.upper_energy))
                       ? prev_ref : roi.lower_energy;

      // Use LLS for continuum if all continuum params are being fit (or NoOffset/External)
      if( (roi.offset_type == PeakContinuum::OffsetType::NoOffset)
          || (roi.offset_type == PeakContinuum::OffsetType::External) )
      {
        roi.use_lls_for_cont = true;
      }
      else
      {
        roi.use_lls_for_cont = true;
        for( const bool fit_par : cont->fitForParameter() )
        {
          if( !fit_par )
          {
            roi.use_lls_for_cont = false;
            break;
          }
        }
      }

      // Max sigma across all peaks in this ROI
      roi.max_initial_sigma = 1.0;
      for( const auto &p : roi.peaks )
      {
        const double s = p->sigma();
        if( !isinf(s) && !isnan(s) && (s > 0.01) )
          roi.max_initial_sigma = std::max( roi.max_initial_sigma, s );
      }

      // Count how many peaks have Sigma being fit
      roi.num_fit_sigmas = 0;
      for( const auto &p : roi.peaks )
        roi.num_fit_sigmas += p->fitFor( PeakDef::Sigma ) ? 1 : 0;

      rois.push_back( std::move(roi) );
    }//for( auto &kv : cont_to_peaks )

    // Sort ROIs by lower energy
    std::sort( rois.begin(), rois.end(), []( const RoiInfo &a, const RoiInfo &b ){
      return a.lower_energy < b.lower_energy;
    } );

    return rois;
  }//make_rois(...)


  PeakFitDiffCostFunction( const std::shared_ptr<const SpecUtils::Measurement> data,
                           const std::vector<std::shared_ptr<const PeakDef>> &starting_peaks,
                           const double roi_lower_energy,    // kept for API compat, ignored (ROI info comes from peaks)
                           const double roi_upper_energy,    // kept for API compat, ignored
                           const double continuum_ref_energy, // kept for API compat, ignored
                           const PeakDef::SkewType skew_type,
                           const bool isHPGe,
                           const Wt::WFlags<PeakFitLM::PeakFitLMOptions> options )
  : m_data( data ),
    m_skew_type( skew_type ),
    m_isHPGe( isHPGe ),
    m_options( options ),
    m_ncalls( 0 ),
    m_rois( make_rois( data, starting_peaks, options ) ),
    m_total_num_peaks( ([this]() -> size_t {
      size_t n = 0;
      for( const RoiInfo &r : m_rois )
        n += r.peaks.size();
      return n;
    })() ),
    m_fit_skew_energy_dependence( ([this, &skew_type]() -> bool {
      // IndependentSkewValues uses per-ROI skew blocks, not a shared energy-dependent block
      if( m_options.testFlag( PeakFitLM::PeakFitLMOptions::IndependentSkewValues ) )
        return false;
      if( m_rois.size() <= 1 )
        return false;
      if( PeakDef::num_skew_parameters( skew_type ) == 0 )
        return false;
      // Check if any skew parameters are energy-dependent
      bool any_energy_dep = false;
      for( size_t i = 0; i < PeakDef::num_skew_parameters( skew_type ); ++i )
      {
        const auto ct = PeakDef::CoefficientType( static_cast<int>(PeakDef::SkewPar0) + static_cast<int>(i) );
        if( PeakDef::is_energy_dependent( skew_type, ct ) )
        {
          any_energy_dep = true;
          break;
        }
      }
      if( !any_energy_dep )
        return false;
      // Only energy-dependent if ROIs span more than 100 keV
      const double lower = m_rois.front().lower_energy;
      const double upper = m_rois.back().upper_energy;
      return (upper - lower) >= 100.0;
    })() ),
    m_skew_anchor_lower_energy( data->gamma_channel_lower( 0 ) ),
    m_skew_anchor_upper_energy( data->gamma_channel_upper( data->num_gamma_channels() - 1 ) ),
    m_external_continuum( ([this]() -> std::shared_ptr<const SpecUtils::Measurement> {
      for( const RoiInfo &roi : m_rois )
      {
        if( roi.offset_type == PeakContinuum::OffsetType::External )
        {
          for( const auto &p : roi.peaks )
          {
            if( p->continuum()->externalContinuum() )
              return p->continuum()->externalContinuum();
          }
        }
      }
      return nullptr;
    })() ),
    m_num_parameters( number_parameters() ),
    m_num_residuals( number_residuals() )
#if( PEAK_FIT_LM_PARALLEL_ROIS )
    , m_thread_pool( (m_rois.size() > 2)
        ? std::make_unique<boost::asio::thread_pool>(
            std::min( m_rois.size(),
                      static_cast<size_t>( std::max( 1, SpecUtilsAsync::num_physical_cpu_cores() ) ) ) )
        : nullptr )
#endif
  {
    if( !m_data || !m_data->gamma_counts() || m_data->gamma_counts()->empty() )
      throw runtime_error( "PeakFitDiffCostFunction: !m_data or !gamma_counts()" );

    if( m_rois.empty() )
      throw runtime_error( "PeakFitDiffCostFunction: no ROIs to fit" );

    const std::vector<float> &counts = *m_data->gamma_counts();
    const std::shared_ptr<const vector<float>> &energies = m_data->channel_energies();

    for( const RoiInfo &roi : m_rois )
    {
      if( roi.lower_energy >= roi.upper_energy )
        throw runtime_error( "PeakFitDiffCostFunction: ROI lower_energy >= upper_energy" );
      if( roi.upper_channel < roi.lower_channel )
        throw runtime_error( "PeakFitDiffCostFunction: ROI upper_channel < lower_channel" );
      if( roi.upper_channel >= counts.size() )
        throw runtime_error( "PeakFitDiffCostFunction: ROI upper_channel >= counts.size()" );
      if( !energies || (roi.upper_channel >= energies->size()) )
        throw runtime_error( "PeakFitDiffCostFunction: ROI upper_channel >= energies.size()" );
      if( roi.peaks.empty() )
        throw runtime_error( "PeakFitDiffCostFunction: ROI has no peaks" );
      if( roi.peaks.size() >= (roi.upper_channel - roi.lower_channel) )
        throw runtime_error( "PeakFitDiffCostFunction: too many peaks for ROI channel range" );

      if( (roi.offset_type == PeakContinuum::OffsetType::External)
          && (!m_external_continuum
              || (m_external_continuum->num_gamma_channels() < 7)
              || !m_external_continuum->energy_calibration()
              || !m_external_continuum->energy_calibration()->valid()) )
        throw runtime_error( "PeakFitDiffCostFunction: external continuum wanted but not valid" );
    }//for( const RoiInfo &roi : m_rois )

    // Check all ROIs have positive DOF
    for( size_t i = 0; i < m_rois.size(); ++i )
    {
      if( dof_for_roi( i ) <= 0.0 )
        throw runtime_error( "PeakFitDiffCostFunction: not enough data channels to fit parameters in ROI" );
    }
  }//PeakFitDiffCostFunction constructor


  // Returns count of sigma parameters for one ROI
  size_t roi_sigma_parameter_count( const RoiInfo &roi ) const
  {
    if( m_options.testFlag( PeakFitLM::PeakFitLMOptions::AllPeakFwhmIndependent ) )
      return roi.num_fit_sigmas;
    else
      return std::min( roi.num_fit_sigmas, size_t(2) );
  }

  // Returns total parameter count for one ROI (cont + sigma + mean params, plus per-ROI skew
  // when IndependentSkewValues is set).
  size_t roi_parameter_count( const RoiInfo &roi ) const
  {
    const size_t cont_pars = roi.use_lls_for_cont
                             ? size_t(0) : PeakContinuum::num_parameters( roi.offset_type );
    size_t n = cont_pars + roi_sigma_parameter_count( roi ) + roi.peaks.size();
    if( m_options.testFlag( PeakFitLM::PeakFitLMOptions::IndependentSkewValues ) )
      n += PeakDef::num_skew_parameters( m_skew_type );
    return n;
  }

  // Returns total residual count for one ROI (data channels + punishment residuals)
  size_t roi_residual_count( const RoiInfo &roi ) const
  {
    const size_t n_chan = roi.upper_channel - roi.lower_channel + 1;
    const bool punish_close = !m_options.testFlag( PeakFitLM::PeakFitLMOptions::DoNotPunishForBeingToClose );
    const bool punish_insig = m_options.testFlag( PeakFitLM::PeakFitLMOptions::PunishForPeakBeingStatInsig );
    size_t n = n_chan;
    if( (punish_close || punish_insig) && roi.peaks.size() > 1 )
      n += roi.peaks.size() - 1;
    if( punish_insig )
      n += 1;
    return n;
  }

  // Returns count of shared skew parameters (doubles for energy-dep params when multi-ROI).
  // Returns 0 when IndependentSkewValues is set, since skew params are folded into each ROI block.
  size_t skew_parameter_count() const
  {
    if( m_options.testFlag( PeakFitLM::PeakFitLMOptions::IndependentSkewValues ) )
      return 0;

    const size_t num_skew = PeakDef::num_skew_parameters( m_skew_type );
    if( !m_fit_skew_energy_dependence )
      return num_skew;

    // Layout: [sp0..spN-1 (lower/base values)] [upper values for energy-dep params only]
    size_t num_energy_dep = 0;
    for( size_t i = 0; i < num_skew; ++i )
    {
      const auto ct = PeakDef::CoefficientType( static_cast<int>(PeakDef::SkewPar0) + static_cast<int>(i) );
      if( PeakDef::is_energy_dependent( m_skew_type, ct ) )
        num_energy_dep += 1;
    }
    return num_skew + num_energy_dep;
  }

  // Returns degrees of freedom for a single ROI.
  // Shared skew parameters are not counted here (they span all ROIs).
  // When IndependentSkewValues is set, per-ROI skew parameters are counted here.
  double dof_for_roi( const size_t roi_index ) const
  {
    assert( roi_index < m_rois.size() );
    const RoiInfo &roi = m_rois[roi_index];
    const double num_channels = static_cast<double>( 1 + roi.upper_channel - roi.lower_channel );
    const size_t num_fit_cont = roi.use_lls_for_cont
                                ? size_t(0) : PeakContinuum::num_parameters( roi.offset_type );

    size_t num_fixed = 0;
    for( const auto &p : roi.peaks )
    {
      num_fixed += p->fitFor( PeakDef::GaussAmplitude ) ? 0 : 1;
      num_fixed += p->fitFor( PeakDef::Mean ) ? 0 : 1;
    }
    if( !roi.use_lls_for_cont )
    {
      for( const bool fit : roi.peaks[0]->continuum()->fitForParameter() )
        num_fixed += fit ? 0 : 1;
    }

    // Count fitted per-ROI skew parameters when IndependentSkewValues is active.
    // Mirrors setup_roi_parameters exactly: a parameter is counted as fitted if ANY
    // matching-type peak in the ROI has fitFor==true for it (OR semantics).
    size_t num_fit_skew = 0;
    if( m_options.testFlag( PeakFitLM::PeakFitLMOptions::IndependentSkewValues ) )
    {
      const size_t num_skew_pars = PeakDef::num_skew_parameters( m_skew_type );
      if( num_skew_pars > 0 )
      {
        vector<bool> fit_skew( num_skew_pars, false );
        for( const auto &p : roi.peaks )
        {
          if( p->skewType() != m_skew_type )
            continue;
          for( size_t i = 0; i < num_skew_pars; ++i )
          {
            const auto ct = PeakDef::CoefficientType( static_cast<int>(PeakDef::SkewPar0) + static_cast<int>(i) );
            if( p->fitFor( ct ) )
              fit_skew[i] = true;
          }
        }
        for( const bool fit : fit_skew )
          num_fit_skew += fit ? 1 : 0;
      }
    }

    return num_channels
           - 2.0*static_cast<double>( roi.peaks.size() )
           + static_cast<double>( num_fixed )
           - static_cast<double>( roi_sigma_parameter_count( roi ) )
           - static_cast<double>( num_fit_cont )
           - static_cast<double>( num_fit_skew );
  }

  // Returns all starting peaks flattened across all ROIs (for skew param initialization)
  std::vector<std::shared_ptr<const PeakDef>> all_starting_peaks() const
  {
    std::vector<std::shared_ptr<const PeakDef>> all;
    for( const RoiInfo &roi : m_rois )
      for( const auto &p : roi.peaks )
        all.push_back( p );
    return all;
  }

  // Computes total parameter count from scratch (called from initializer list)
  size_t number_parameters() const
  {
    size_t n = skew_parameter_count();
    for( const RoiInfo &roi : m_rois )
      n += roi_parameter_count( roi );
    return n;
  }

  // Computes total residual count from scratch (called from initializer list)
  size_t number_residuals() const
  {
    size_t n = 0;
    for( const RoiInfo &roi : m_rois )
      n += roi_residual_count( roi );
    return n;
  }


  //
  template<typename T>
  static std::shared_ptr<T> create_continuum( const std::shared_ptr<T> &other_cont [[maybe_unused]]) {
      return std::make_shared<T>();
  }

  /** Apply the shared skew parameters to all peaks in the given vector.
   For single-ROI (or no energy dependence), all peaks get the same skew values.
   For multi-ROI with energy-dependent params, the skew is interpolated between
   the lower and upper anchor energies based on each peak's mean energy.

   Skew parameter layout in params[0..skew_block_size-1]:
     params[0..num_skew-1]:        base (lower-anchor) values for all skew params
     params[num_skew..num_skew+M-1]: upper-anchor values for energy-dep params only (M = count of energy-dep params)

   uncertainties: if non-null, has the same layout as params (diagonal of covariance); used to
     set skew uncertainty on each peak for non-energy-dependent params.
   covariance: if non-null, the full num_total_pars x num_total_pars row-major covariance matrix
     for the entire parameter vector; needed for proper error propagation of interpolated
     energy-dependent skew params.  The indices in this matrix correspond to the global parameter
     array, so params[0] corresponds to covariance[0][0], etc.
   params_offset: the index of params[0] in the global parameter array (needed to index covariance).
  */
  template<typename PeakType, typename T>
  void apply_skew_to_peaks( vector<PeakType> &peaks, const T * const params,
                            const T * const uncertainties = nullptr,
                            const double * const covariance = nullptr,
                            const size_t num_total_pars = 0,
                            const size_t params_offset = 0 ) const
  {
    const size_t num_skew = PeakDef::num_skew_parameters( m_skew_type );
    if( num_skew == 0 )
      return;

    const double energy_span = m_skew_anchor_upper_energy - m_skew_anchor_lower_energy;

    for( PeakType &peak : peaks )
    {
      peak.setSkewType( m_skew_type );

      if( !m_fit_skew_energy_dependence )
      {
        // Single ROI or no energy dependence: all peaks share the same skew values
        for( size_t i = 0; i < num_skew; ++i )
        {
          const auto ct = PeakDef::CoefficientType( static_cast<int>(PeakDef::SkewPar0) + static_cast<int>(i) );
          peak.set_coefficient( params[i], ct );
          if( uncertainties )
            peak.set_uncertainty( uncertainties[i], ct );
        }
      }
      else
      {
        // Multi-ROI with energy dependence: interpolate per-peak based on mean
        T mean_val;
        if constexpr ( std::is_same_v<T, double> )
          mean_val = T( peak.mean() );
        else
          mean_val = peak.mean();

        const T mean_frac = (mean_val - T(m_skew_anchor_lower_energy)) / T(energy_span);

        size_t upper_offset = num_skew; // index of first upper-anchor param
        for( size_t i = 0; i < num_skew; ++i )
        {
          const auto ct = PeakDef::CoefficientType( static_cast<int>(PeakDef::SkewPar0) + static_cast<int>(i) );
          T val;
          if( PeakDef::is_energy_dependent( m_skew_type, ct ) )
          {
            // Interpolate between lower (params[i]) and upper (params[upper_offset])
            val = params[i] + mean_frac * (params[upper_offset] - params[i]);
            peak.set_coefficient( val, ct );

            // Uncertainty propagation only applies when T=double (post-solve), not during Jet-based solve.
            if constexpr ( std::is_same_v<T, double> )
            {
              // Proper error propagation: val = (1-f)*p_lower + f*p_upper
              // sigma^2 = (1-f)^2*Var[lower] + f^2*Var[upper] + 2*(1-f)*f*Cov[lower,upper]
              if( covariance && (num_total_pars > 0) )
              {
                const double f = mean_frac;
                const double one_minus_f = 1.0 - f;
                const size_t gi = params_offset + i;              // global index of lower-anchor param
                const size_t gu = params_offset + upper_offset;   // global index of upper-anchor param
                const double var_lower = covariance[gi * num_total_pars + gi];
                const double var_upper = covariance[gu * num_total_pars + gu];
                const double cov_lu    = covariance[gi * num_total_pars + gu];
                const double variance  = one_minus_f*one_minus_f*var_lower
                                         + f*f*var_upper
                                         + 2.0*one_minus_f*f*cov_lu;
                if( variance > 0.0 )
                  peak.set_uncertainty( sqrt(variance), ct );
              }
              else if( uncertainties )
              {
                // Fallback: propagate in quadrature ignoring correlation
                const double f = mean_frac;
                const double sigma_lower = uncertainties[i];
                const double sigma_upper = uncertainties[upper_offset];
                const double sigma = sqrt( (1.0-f)*(1.0-f)*sigma_lower*sigma_lower
                                           + f*f*sigma_upper*sigma_upper );
                if( sigma > 0.0 )
                  peak.set_uncertainty( sigma, ct );
              }
            }//if constexpr T==double

            upper_offset += 1;
          }
          else
          {
            val = params[i];
            peak.set_coefficient( val, ct );
            if( uncertainties )
              peak.set_uncertainty( uncertainties[i], ct );
          }
        }
      }
    }//for( PeakType &peak : peaks )
  }//apply_skew_to_peaks(...)


  template<typename PeakType,typename T>
  vector<PeakType> parametersToPeaks( const T * const params, const T * const uncertainties, T *residuals,
                                      const double * const covariance = nullptr,
                                      const size_t num_total_pars = 0 ) const
  {
    std::unique_ptr<vector<T>> local_residuals;
    if( !residuals )
    {
      local_residuals.reset( new vector<T>( number_residuals(), T(0.0) ) );
      residuals = local_residuals->data();
    }
    else
    {
      const size_t nresids = number_residuals();
      for( size_t i = 0; i < nresids; ++i )
        residuals[i] = T(0.0);
    }

    const size_t num_skew = PeakDef::num_skew_parameters( m_skew_type );
    const size_t skew_block_size = skew_parameter_count();

    // Pre-compute per-ROI parameter and residual offsets so we can launch ROIs in parallel.
    const size_t nrois = m_rois.size();
    vector<size_t> roi_param_offsets( nrois ), roi_residual_offsets( nrois );
    {
      size_t poff = skew_block_size, roff = 0;
      for( size_t i = 0; i < nrois; ++i )
      {
        roi_param_offsets[i]    = poff;
        roi_residual_offsets[i] = roff;
        poff += roi_parameter_count( m_rois[i] );
        roff += roi_residual_count( m_rois[i] );
      }
    }

    // Lambda that processes one ROI and returns its peaks; writes residuals into roi_residuals.
    // All inputs are read-only except roi_residuals (which is per-ROI private storage).
    // params[] is read-only (shared skew block + ROI params), safe to read from multiple threads.
    const auto process_one_roi = [&]( const size_t roi_idx,
                                      const size_t param_offset,
                                      T * const roi_residuals ) -> vector<PeakType>
    {
      const RoiInfo &roi = m_rois[roi_idx];
      const size_t num_roi_peaks = roi.peaks.size();
      const size_t num_sigmas_fit = roi_sigma_parameter_count( roi );
      const size_t num_fit_cont = roi.use_lls_for_cont
                                  ? size_t(0) : PeakContinuum::num_parameters( roi.offset_type );
      const size_t nchannel = roi.upper_channel - roi.lower_channel + 1;
      const double range = roi.upper_energy - roi.lower_energy;

      // Pointer into the parameter array for this ROI
      const T * const roi_params = params + param_offset;

      vector<PeakType> peaks, fixed_amp_peaks;

      // Compute min/max means across this ROI for sigma interpolation
      T min_mean = T( std::numeric_limits<double>::infinity() );
      T max_mean = T( -std::numeric_limits<double>::infinity() );
      for( size_t i = 0; i < num_roi_peaks; ++i )
      {
        const size_t mean_par_index = num_fit_cont + num_sigmas_fit + i;
        const T frac = roi_params[mean_par_index] - T(0.5);
        const T mean = T(roi.lower_energy) + frac * T(range);
        min_mean = min( min_mean, mean );
        max_mean = max( max_mean, mean );
      }

      size_t fit_sigma_num = 0;
      for( size_t i = 0; i < num_roi_peaks; ++i )
      {
        const std::shared_ptr<const PeakDef> &src_peak = roi.peaks[i];
        const size_t mean_par_index = num_fit_cont + num_sigmas_fit + i;
        assert( mean_par_index < roi_parameter_count( roi ) );

        const bool fit_amp = src_peak->fitFor( PeakDef::GaussAmplitude );

        const T frac = roi_params[mean_par_index] - T(0.5);
        const T mean = T(roi.lower_energy) + frac * T(range);
        const T amp = fit_amp ? T(1.0) : T(src_peak->amplitude());

        T sigma, sigma_uncert;

        if( !src_peak->fitFor( PeakDef::Sigma ) || (num_sigmas_fit == 0) )
        {
          sigma = T( src_peak->sigma() );
          sigma_uncert = T( src_peak->sigmaUncert() );
        }
        else if( m_options.testFlag( PeakFitLM::PeakFitLMOptions::AllPeakFwhmIndependent ) )
        {
          assert( fit_sigma_num < roi.num_fit_sigmas );
          const size_t sigma_index = num_fit_cont + fit_sigma_num;
          assert( sigma_index < roi_parameter_count( roi ) );

          sigma = roi_params[sigma_index] * T(roi.max_initial_sigma);
          sigma_uncert = uncertainties
                         ? (uncertainties[param_offset + sigma_index] * T(roi.max_initial_sigma))
                         : T(0.0);

          if( isnan(sigma) || isinf(sigma) )
            throw runtime_error( "Inf or NaN sigma (AllPeakFwhmIndependent)" );
        }
        else
        {
          // First sigma param = absolute sigma for lowest-energy peak;
          // second sigma param (if present) = relative multiplier for highest-energy peak.
          const size_t sigma_index = num_fit_cont;
          sigma = roi_params[sigma_index] * T(roi.max_initial_sigma);
          sigma_uncert = uncertainties
                         ? (uncertainties[param_offset + sigma_index] * T(roi.max_initial_sigma))
                         : T(0.0);

          if( (i > 0) && (num_sigmas_fit > 1) )
          {
            T frac_dist;
            const T mean_dists = max_mean - min_mean;
            if( mean_dists > T(1.0E-6) )
              frac_dist = (mean - min_mean) / mean_dists;
            else
              frac_dist = T(0.5);

            const T first_sigma = sigma;
            const T last_sigma  = sigma * roi_params[sigma_index + 1];
            sigma = first_sigma + frac_dist * (last_sigma - first_sigma);

            if( uncertainties )
            {
              const T first_uncert = sigma_uncert;
              const T last_uncert  = sigma_uncert * uncertainties[param_offset + sigma_index + 1];
              sigma_uncert = first_uncert + frac_dist * (last_uncert - first_uncert);
            }
          }//if( (i > 0) && (num_sigmas_fit > 1) )

          if( isnan(sigma) || isinf(sigma) )
            throw runtime_error( "Inf or NaN sigma (shared sigma)" );
        }//sigma selection

        // Protect against sigma smaller than one channel width (only when fitting sigma)
        if( (sigma < 0.2) && src_peak->fitFor( PeakDef::Sigma ) && (num_sigmas_fit != 0) )
        {
          const float left_val  = m_data->gamma_channel_lower( roi.lower_channel );
          const float upper_val = m_data->gamma_channel_upper( roi.upper_channel );
          const double avrg_chnl_width = (upper_val - left_val) / static_cast<double>( nchannel );
          const double min_reasonable = (0.42*0.99) * avrg_chnl_width;

          double sigma_d;
          if constexpr ( std::is_same_v<T, double> )
            sigma_d = sigma;
          else
            sigma_d = sigma.a;

          if( (sigma_d <= min_reasonable) || (sigma_d < 0.00025) )
            throw runtime_error( "peak sigma (" + std::to_string(sigma_d)
                                 + ") smaller than reasonable ("
                                 + std::to_string( min_reasonable ) + ")" );
        }

        if( src_peak->fitFor( PeakDef::Sigma ) )
          fit_sigma_num += 1;

        PeakType peak;
        peak.setMean( mean );
        peak.setSigma( sigma );
        peak.setAmplitude( amp );
        // Skew applied later via apply_skew_to_peaks
        peak.setSkewType( m_skew_type );

        if( uncertainties )
        {
          T mean_uncert = uncertainties[param_offset + mean_par_index];
          if( !src_peak->fitFor( PeakDef::Mean ) )
            mean_uncert = T( src_peak->meanUncert() );
          if( mean_uncert > 0.0 )
            peak.setMeanUncert( mean_uncert );
          if( sigma_uncert > 0.0 )
            peak.setSigmaUncert( sigma_uncert );
        }

        if( fit_amp )
          peaks.push_back( peak );
        else
        {
          peak.setAmplitudeUncert( T( src_peak->amplitudeUncert() ) );
          fixed_amp_peaks.push_back( peak );
        }
      }//for( size_t i = 0; i < num_roi_peaks; ++i )

      assert( fit_sigma_num == roi.num_fit_sigmas );
      std::sort( begin(peaks), end(peaks), &PeakType::lessThanByMean );

      // --- Compute predicted channel counts for this ROI ---
      const shared_ptr<const vector<float>> &energies_ptr = m_data->channel_energies();
      const vector<float> &counts_vec = *m_data->gamma_counts();
      const float * const energies       = energies_ptr->data() + roi.lower_channel;
      const float * const channel_counts = counts_vec.data()    + roi.lower_channel;

      const int  num_polynomial_terms = static_cast<int>( PeakContinuum::num_parameters( roi.offset_type ) );
      const bool is_step_continuum    = PeakContinuum::is_step_continuum( roi.offset_type );

      assert( !peaks.empty() || !fixed_amp_peaks.empty() );
      auto continuum = create_continuum( (!peaks.empty())
                                         ? peaks.front().getContinuum()
                                         : fixed_amp_peaks.front().getContinuum() );
      continuum->setRange( T(roi.lower_energy), T(roi.upper_energy) );
      continuum->setType( roi.offset_type );

      vector<T> peak_counts( nchannel, T(0.0) );

      // Build the skew parameter vector to pass to fit_amp_and_offset_imp.
      // In IndependentSkewValues mode the skew params are at the tail of this ROI's own block;
      // otherwise they live at params[0..num_skew-1] (the shared / energy-dependent block).
      const T *roi_skew_ptr;
      if( m_options.testFlag( PeakFitLM::PeakFitLMOptions::IndependentSkewValues ) )
        roi_skew_ptr = roi_params + num_fit_cont + num_sigmas_fit + num_roi_peaks;
      else
        roi_skew_ptr = params;
      const vector<T> skew_pars( roi_skew_ptr, roi_skew_ptr + num_skew );

      if( roi.offset_type == PeakContinuum::OffsetType::External )
      {
        assert( m_external_continuum );
        continuum->setExternalContinuum( m_external_continuum );

        if( !peaks.empty() )
        {
          // Subtract external continuum from data for amplitude fitting
          vector<float> data_copy( channel_counts, channel_counts + nchannel );
          vector<float> data_variances( channel_counts, channel_counts + nchannel );
          for( size_t i = 0; i < nchannel; ++i )
          {
            data_copy[i] -= m_external_continuum->gamma_integral( energies[i], energies[i+1] );
            data_copy[i] = std::max( 0.0f, data_copy[i] );
            if( data_variances[i] < PEAK_FIT_MIN_CHANNEL_UNCERT )
              data_variances[i] = PEAK_FIT_MIN_CHANNEL_UNCERT;
          }

          vector<T> means, sigmas;
          for( const auto &p : peaks )
          {
            means.push_back( p.mean() );
            sigmas.push_back( p.sigma() );
          }

          vector<T> amplitudes, cont_coeffs, amp_uncerts, cont_uncerts;
          PeakFit::fit_amp_and_offset_imp( energies, &data_copy[0], &data_variances[0], nchannel,
                                           num_polynomial_terms, is_step_continuum, T(roi.ref_energy),
                                           means, sigmas, fixed_amp_peaks, m_skew_type, skew_pars.data(),
                                           amplitudes, cont_coeffs, amp_uncerts, cont_uncerts,
                                           &peak_counts[0] );
          assert( peaks.size() == amplitudes.size() );
          for( size_t pi = 0; pi < peaks.size(); ++pi )
          {
            peaks[pi].setAmplitude( amplitudes[pi] );
            if( !amp_uncerts.empty() )
              peaks[pi].setAmplitudeUncert( amp_uncerts[pi] );
          }
        }
        else
        {
          for( PeakType &fp : fixed_amp_peaks )
            fp.gauss_integral( energies, &peak_counts[0], nchannel );
        }

        for( size_t i = 0; i < nchannel; ++i )
          peak_counts[i] += T( m_external_continuum->gamma_integral( energies[i], energies[i+1] ) );
      }
      else if( !roi.use_lls_for_cont )
      {
        // Non-LLS continuum: continuum parameters are in roi_params[0..num_fit_cont-1]
        continuum->setParameters( T(roi.ref_energy), roi_params, nullptr );

        for( PeakType &p : peaks )
          p.gauss_integral( energies, &peak_counts[0], nchannel );
        for( PeakType &fp : fixed_amp_peaks )
          fp.gauss_integral( energies, &peak_counts[0], nchannel );

        PeakDists::offset_integral( *continuum, energies, &peak_counts[0], nchannel, m_data );
      }
      else
      {
        // LLS continuum: fit continuum coefficients simultaneously with amplitudes
        vector<T> means, sigmas;
        for( const auto &p : peaks )
        {
          means.push_back( p.mean() );
          sigmas.push_back( p.sigma() );
        }

        vector<T> amplitudes, cont_coeffs, amp_uncerts, cont_uncerts;
        PeakFit::fit_amp_and_offset_imp( energies, channel_counts, nullptr, nchannel,
                                         num_polynomial_terms, is_step_continuum, T(roi.ref_energy),
                                         means, sigmas, fixed_amp_peaks, m_skew_type, skew_pars.data(),
                                         amplitudes, cont_coeffs, amp_uncerts, cont_uncerts,
                                         &peak_counts[0] );
        assert( peaks.size() == amplitudes.size() );
        for( size_t pi = 0; pi < peaks.size(); ++pi )
        {
          peaks[pi].setAmplitude( amplitudes[pi] );
          if( !amp_uncerts.empty() )
            peaks[pi].setAmplitudeUncert( amp_uncerts[pi] );
        }

        if( (roi.offset_type != PeakContinuum::OffsetType::NoOffset)
            && (roi.offset_type != PeakContinuum::OffsetType::External) )
        {
          continuum->setParameters( T(roi.ref_energy), cont_coeffs.data(), cont_uncerts.data() );
        }
      }//continuum handling

      // Merge fixed-amp peaks back into peaks
      if( !fixed_amp_peaks.empty() )
      {
        for( PeakType &fp : fixed_amp_peaks )
          peaks.push_back( fp );
        std::sort( begin(peaks), end(peaks), &PeakType::lessThanByMean );
      }

      // Set continuum and apply skew parameters to all peaks in this ROI.
      // In IndependentSkewValues mode, pass the per-ROI skew pointer; otherwise pass params
      // (which points to the shared/energy-dependent skew block at the start of the global array).
      for( PeakType &p : peaks )
        p.setContinuum( continuum );

      // Pass uncertainties at the same offset as the skew params, so uncertainties are set too.
      const T *roi_skew_uncert_ptr = uncertainties ? (uncertainties + (roi_skew_ptr - params)) : nullptr;
      const size_t skew_params_offset = static_cast<size_t>( roi_skew_ptr - params );
      apply_skew_to_peaks<PeakType, T>( peaks, roi_skew_ptr, roi_skew_uncert_ptr,
                                        covariance, num_total_pars, skew_params_offset );

      // --- Compute residuals for this ROI ---
      T chi2( 0.0 );
      for( size_t ch = 0; ch < nchannel; ++ch )
      {
        const double ndata = channel_counts[ch];
        if( ndata >= PEAK_FIT_MIN_CHANNEL_UNCERT )
          roi_residuals[ch] = (T(ndata) - peak_counts[ch]) / T(sqrt(ndata));
        else
          roi_residuals[ch] = peak_counts[ch]; // ad-hoc, follows PeakFitChi2Fcn
        chi2 += roi_residuals[ch] * roi_residuals[ch];
      }

      const double dof_val = dof_for_roi( roi_idx );
      const T chi2_dof = (dof_val > 0.0) ? (chi2 / T(dof_val)) : T(0.0);
      for( PeakType &p : peaks )
        p.set_coefficient( chi2_dof, PeakDef::Chi2DOF );

      // --- Punishment residuals (within this ROI only) ---
      // Punishment residuals immediately follow the channel residuals for this ROI.
      // Indexing: roi_residuals[nchannel + i - 1] for peaks[i] vs peaks[i-1].
      const double punishment_factor = 0.25 * static_cast<double>(nchannel) / static_cast<double>(num_roi_peaks);

      for( size_t i = 1; i < peaks.size(); ++i )
      {
        T avg_sigma;
        if constexpr ( !std::is_same_v<T, double> )
        {
          avg_sigma = T(0.5) * (peaks[i-1].sigma() + peaks[i].sigma());
        }
        else
        {
          const double s0 = peaks[i-1].gausPeak() ? peaks[i-1].sigma() : 0.25*peaks[i-1].roiWidth();
          const double s1 = peaks[i].gausPeak()   ? peaks[i].sigma()   : 0.25*peaks[i].roiWidth();
          avg_sigma = 0.5 * (s0 + s1);
        }

        if( !m_options.testFlag( PeakFitLM::PeakFitLMOptions::DoNotPunishForBeingToClose ) )
        {
          const T dist = abs( peaks[i-1].mean() - peaks[i].mean() );
          T reldist = dist / avg_sigma;
          if( (reldist < 0.01) || isinf(reldist) || isnan(reldist) )
          {
            if constexpr ( !std::is_same_v<T, double> )
            {
              reldist.a = 0.01;
              reldist.v = peaks[i-1].mean().v - peaks[i].mean().v;
            }
            else
            {
              reldist = 0.01;
            }
          }

          const size_t punish_idx = nchannel + i - 1;
          assert( punish_idx < roi_residual_count( roi ) );
          if( reldist < 1.25 )
            roi_residuals[punish_idx] = T(punishment_factor / reldist);
          else
            roi_residuals[punish_idx] = T(0.0);
        }//if( punish for peaks too close )

        if( m_options.testFlag( PeakFitLM::PeakFitLMOptions::PunishForPeakBeingStatInsig ) )
        {
          // As of 20250415, this section is untested.
          // Use last slot of punishment residuals for the insignificance penalty
          const size_t last_punish_idx = nchannel + peaks.size() - 1;
          assert( last_punish_idx < roi_residual_count( roi ) );
          roi_residuals[last_punish_idx] = T(0.0);  // initialized once; summed below

          for( size_t pi = 0; pi < peaks.size(); ++pi )
          {
            const PeakType &pk = peaks[pi];

            if( m_options.testFlag( PeakFitLM::PeakFitLMOptions::DoNotPunishForBeingToClose ) )
            {
              const size_t idx = nchannel + pi;
              assert( idx < roi_residual_count( roi ) );
              roi_residuals[idx] = T(0.0);
            }

            double pk_mean, pk_sigma, pk_amp;
            if constexpr ( std::is_same_v<T, double> )
            {
              pk_mean  = pk.mean();
              pk_sigma = pk.sigma();
              pk_amp   = pk.amplitude();
            }
            else
            {
              pk_mean  = pk.mean().a;
              pk_sigma = pk.sigma().a;
              pk_amp   = pk.amplitude().a;
            }

            const float lower_energy = static_cast<float>( pk_mean - 1.75*pk_sigma );
            const float upper_energy = static_cast<float>( pk_mean + 1.75*pk_sigma );
            const size_t lower_ch = std::max( roi.lower_channel, m_data->find_gamma_channel(lower_energy) );
            const size_t upper_ch = std::min( roi.upper_channel, m_data->find_gamma_channel(upper_energy) );

            double dataarea = 0.0;
            for( size_t bin = lower_ch; (bin < counts_vec.size()) && (bin <= upper_ch); ++bin )
              dataarea += counts_vec[bin];

            const double punishment_chi2 = 2.0 * static_cast<double>( 1 + upper_ch - lower_ch );
            const size_t this_punish_idx = nchannel + pi;

            if( pk_amp < std::max( 2.0*sqrt(dataarea), 1.0 ) )
            {
              if( pk_amp <= 1.0 )
              {
                if constexpr ( std::is_same_v<T, double> )
                {
                  roi_residuals[this_punish_idx] += 100.0 * punishment_chi2;
                }
                else
                {
                  T punish( 1.0 );
                  punish.v = pk.amplitude().v / pk.amplitude().a;
                  punish *= 100.0 * punishment_chi2;
                  roi_residuals[this_punish_idx] += punish;
                }
              }
              else
              {
                roi_residuals[this_punish_idx] += (T(-1.0) + 2.0*sqrt(dataarea)/pk.amplitude()) * T(punishment_chi2);
              }
            }//if( amp < 2*sqrt(dataarea) )
          }//for( size_t pi = 0; pi < peaks.size(); ++pi )
        }//if( PunishForPeakBeingStatInsig )
      }//for( size_t i = 1; i < peaks.size(); ++i )

      // --- Inherit user-selected options from matching starting peaks ---
      assert( roi.peaks.size() == peaks.size() );
      map<size_t,size_t> old_to_new;

      // Fixed-mean peaks keep their original index
      for( size_t idx = 0; idx < roi.peaks.size(); ++idx )
      {
        if( !roi.peaks[idx]->fitFor( PeakDef::Mean ) )
          old_to_new[idx] = idx;
      }

      for( size_t new_idx = 0; new_idx < peaks.size(); ++new_idx )
      {
        bool already = false;
        for( const auto &kv : old_to_new )
          already = (already || (kv.second == new_idx));
        if( already )
          continue;

        int closest = -1;
        const T new_mean = peaks[new_idx].mean();
        for( size_t orig_idx = 0; orig_idx < roi.peaks.size(); ++orig_idx )
        {
          if( old_to_new.count( orig_idx ) )
            continue;
          if( closest < 0 )
          {
            closest = static_cast<int>( orig_idx );
          }
          else
          {
            const double prev_m = roi.peaks[static_cast<size_t>(closest)]->mean();
            const double orig_m = roi.peaks[orig_idx]->mean();
            double nm;
            if constexpr ( std::is_same_v<T, double> )
              nm = new_mean;
            else
              nm = new_mean.a;
            if( std::abs(nm - orig_m) < std::abs(nm - prev_m) )
              closest = static_cast<int>( orig_idx );
          }
        }
        assert( closest >= 0 );
        old_to_new[static_cast<size_t>(closest)] = new_idx;
      }

      assert( old_to_new.size() == peaks.size() );
      for( const auto &mapping : old_to_new )
      {
        assert( mapping.first < peaks.size() );
        assert( mapping.second < peaks.size() );
        peaks[mapping.second].inheritUserSelectedOptions( *roi.peaks[mapping.first], true );
      }

      // VoigtPlusBortel special case: when global skew is Bortel or GaussPlusBortel
      // but an individual starting peak had VoigtPlusBortel type, remap the skew coefficients.
      // gamma_lor (SkewPar0) and R (SkewPar1, for Bortel global) are taken directly from the
      // stored original peak values (not fitted separately), since they are not in the global
      // parameter set.
      if( (m_skew_type == PeakDef::SkewType::Bortel)
          || (m_skew_type == PeakDef::SkewType::GaussPlusBortel) )
      {
        for( const auto &mapping : old_to_new )
        {
          const std::shared_ptr<const PeakDef> &orig = roi.peaks[mapping.first];
          PeakType &fitted = peaks[mapping.second];

          if( orig->skewType() != PeakDef::SkewType::VoigtPlusBortel )
            continue;

          if( m_skew_type == PeakDef::SkewType::Bortel )
          {
            // Bortel: SkewPar0=tau
            // VoigtPlusBortel: SkewPar0=gamma_lor, SkewPar1=R, SkewPar2=tau
            // roi_skew_ptr[0] is always tau for Bortel, whether shared or per-ROI.
            const T tau = roi_skew_ptr[0];
            fitted.setSkewType( PeakDef::SkewType::VoigtPlusBortel );
            fitted.set_coefficient( T(orig->coefficient(PeakDef::SkewPar0)), PeakDef::SkewPar0 ); // gamma_lor from orig
            fitted.set_coefficient( T(orig->coefficient(PeakDef::SkewPar1)), PeakDef::SkewPar1 ); // R from orig
            fitted.set_coefficient( tau, PeakDef::SkewPar2 );  // tau from skew block
            if( roi_skew_uncert_ptr )
              fitted.set_uncertainty( roi_skew_uncert_ptr[0], PeakDef::SkewPar2 );
          }
          else // GaussPlusBortel
          {
            // GaussPlusBortel: SkewPar0=R, SkewPar1=tau
            // VoigtPlusBortel: SkewPar0=gamma_lor, SkewPar1=R, SkewPar2=tau
            // roi_skew_ptr[0]/[1] are R and tau, whether shared or per-ROI.
            const T R_global   = roi_skew_ptr[0]; // GaussPlusBortel SkewPar0=R
            const T tau_global = roi_skew_ptr[1]; // GaussPlusBortel SkewPar1=tau
            fitted.setSkewType( PeakDef::SkewType::VoigtPlusBortel );
            fitted.set_coefficient( T(orig->coefficient(PeakDef::SkewPar0)), PeakDef::SkewPar0 ); // gamma_lor from orig
            fitted.set_coefficient( R_global,   PeakDef::SkewPar1 );
            fitted.set_coefficient( tau_global, PeakDef::SkewPar2 );
            if( roi_skew_uncert_ptr )
            {
              fitted.set_uncertainty( roi_skew_uncert_ptr[0], PeakDef::SkewPar1 );
              fitted.set_uncertainty( roi_skew_uncert_ptr[1], PeakDef::SkewPar2 );
            }
          }
        }//for( const auto &mapping : old_to_new )
      }//if( global skew is Bortel or GaussPlusBortel )

      return peaks;
    };//process_one_roi lambda

    // Collect all peaks from all ROIs into the final result
    vector<PeakType> all_peaks;
    all_peaks.reserve( m_total_num_peaks );

#if( PEAK_FIT_LM_PARALLEL_ROIS )
    // For >2 ROIs, evaluate each ROI on a thread from m_thread_pool.
    // The pool is sized to physical cores and lives for the lifetime of this cost function,
    // so thread creation cost is paid once at construction, not per Ceres evaluation.
    // SpecUtilsAsync::ThreadPool (GCD backend) is deliberately avoided: creating pools
    // nested inside another pool's workers causes severe delays on macOS (see PeakFit_imp.hpp).
    // Applies to both double and Jet instantiations: Jet arithmetic is ~8x more expensive
    // per channel than double, so it benefits at least as much from parallelism.
    if( nrois > 2 )
    {
      // Per-ROI residual buffers — each thread writes to its own buffer, eliminating races.
      vector<vector<T>> roi_residual_bufs( nrois );
      for( size_t i = 0; i < nrois; ++i )
        roi_residual_bufs[i].assign( roi_residual_count( m_rois[i] ), T(0.0) );

      // One promise/future pair per ROI carries the peak vector (and any exception).
      vector<std::promise<vector<PeakType>>> promises( nrois );
      vector<std::future<vector<PeakType>>>  futures;
      futures.reserve( nrois );
      for( size_t i = 0; i < nrois; ++i )
        futures.push_back( promises[i].get_future() );

      for( size_t i = 0; i < nrois; ++i )
      {
        boost::asio::post( *m_thread_pool,
          [&, i]()
          {
            try
            {
              promises[i].set_value(
                process_one_roi( i, roi_param_offsets[i], roi_residual_bufs[i].data() ) );
            }
            catch( ... )
            {
              promises[i].set_exception( std::current_exception() );
            }
          } );
      }

      // Collect results in ROI order (preserves peak ordering)
      for( size_t i = 0; i < nrois; ++i )
      {
        vector<PeakType> roi_peaks = futures[i].get(); // re-throws any stored exception
        for( PeakType &p : roi_peaks )
          all_peaks.push_back( std::move(p) );

        // Copy ROI residuals into the shared output array (sequential, no race)
        const size_t roff = roi_residual_offsets[i];
        const size_t rcount = roi_residual_count( m_rois[i] );
        for( size_t j = 0; j < rcount; ++j )
          residuals[roff + j] = roi_residual_bufs[i][j];
      }
    }else
#endif //PEAK_FIT_LM_PARALLEL_ROIS
    {
      // Serial path: used when PEAK_FIT_LM_PARALLEL_ROIS==0 or nrois<=2.
      // Each ROI writes directly into the correct slice of `residuals`.
      for( size_t roi_idx = 0; roi_idx < nrois; ++roi_idx )
      {
        T * const roi_res = residuals + roi_residual_offsets[roi_idx];
        vector<PeakType> roi_peaks = process_one_roi( roi_idx, roi_param_offsets[roi_idx], roi_res );
        for( PeakType &p : roi_peaks )
          all_peaks.push_back( std::move(p) );
      }
    }


#if( !defined(NDEBUG) && PERFORM_DEVELOPER_CHECKS && defined(CERES_PUBLIC_JET_H_) )
    const size_t nresids = number_residuals();
    for( size_t i = 0; i < nresids; ++i )
    {
      if constexpr ( std::is_same_v<T, double> )
      {
        assert( !isnan(residuals[i]) && !isinf(residuals[i]) );
      }
      else
      {
        check_jet_for_NaN( residuals[i] );
      }
    }
#endif

    return all_peaks;
  }//parametersToPeaks(...)


  template<typename T>
  bool operator()( T const *const *parameters, T *residuals ) const
  {
    m_ncalls += 1;

    try
    {
      T const * const pars = parameters[0];
      const T * const uncertainties = nullptr;
      if constexpr ( std::is_same_v<T, double> )
      {
        parametersToPeaks<PeakDef,double>( pars, uncertainties, residuals );
      }else
      {
        parametersToPeaks<PeakDefImpWithCont<T>,T>( pars, uncertainties, residuals );
      }
    }catch( std::exception &e )
    {
      cerr << "PeakFitDiffCostFunction: caught exception: '" << e.what() << "'" << endl;
      return false;
    }

    return true;
  }//operator()


  static PeakDef::SkewType skew_type_from_prev_peaks( const vector<shared_ptr<const PeakDef> > &inpeaks )
  {
    PeakDef::SkewType skew_type = PeakDef::SkewType::NoSkew;

    for( const auto &p : inpeaks )
    {
      // PeakDef::SkewType is defined roughly in order of how we should prefer them.
      //  However, if one of the peaks in the ROI required a higher-level of skew function
      //  to describe it, we will use that skew for the entire ROI.
      // Here we will use the skew value of the first peak we encounter, with the higher-level
      //  of skew function, as the starting value, and wether we should fit for the parameter
      //  at all.
      if( p->skewType() > skew_type )
        skew_type = p->skewType();
    }//for( const auto &p : inpeaks )

    return skew_type;
  }//void skew_type_from_prev_peaks(...)


  /** Set up the shared skew parameter block (params[0..skew_block_size-1]).
   Skew params are at the front of the global parameter vector.
   Layout:
     params[0..num_skew-1]:        base (lower-anchor) values for all skew params
     params[num_skew..num_skew+M-1]: upper-anchor values for energy-dep params only
       (only present when m_fit_skew_energy_dependence == true)
   When IndependentSkewValues is set, this function is a no-op (skew is per-ROI).
  */
  void setup_skew_parameters( double *pars,
                              vector<int> &constant_parameters,
                              vector<std::optional<double>> &lower_bounds,
                              vector<std::optional<double>> &upper_bounds,
                              const vector<shared_ptr<const PeakDef>> &inpeaks ) const
  {
    // With IndependentSkewValues, there is no shared skew block; each ROI sets up its own.
    if( m_options.testFlag( PeakFitLM::PeakFitLMOptions::IndependentSkewValues ) )
      return;

    const size_t num_skew_pars = PeakDef::num_skew_parameters( m_skew_type );

    if( num_skew_pars == 0 )
      return;

    vector<bool> fit_parameter( num_skew_pars, true );
    vector<double> starting_value( num_skew_pars, 0.0 );
    vector<double> lower_values( num_skew_pars, 0.0 );
    vector<double> upper_values( num_skew_pars, 0.0 );

    for( size_t i = 0; i < num_skew_pars; ++i )
    {
      const auto ct = PeakDef::CoefficientType( static_cast<int>(PeakDef::SkewPar0) + static_cast<int>(i) );
      double lower, upper, start, dx;
      const bool use = PeakDef::skew_parameter_range( m_skew_type, ct, lower, upper, start, dx );
      assert( use );
      if( !use )
        throw logic_error( "Inconsistent skew par val" );
      starting_value[i] = start;
      lower_values[i]   = lower;
      upper_values[i]   = upper;
    }

    // Use values from the first matching peak (matching global skew type)
    for( const auto &p : inpeaks )
    {
      if( p->skewType() != m_skew_type )
        continue;

      for( size_t i = 0; i < num_skew_pars; ++i )
      {
        const auto ct = PeakDef::CoefficientType( static_cast<int>(PeakDef::SkewPar0) + static_cast<int>(i) );
        fit_parameter[i] = p->fitFor( ct );
        double val = p->coefficient( ct );

        if( IsInf(val) || IsNan(val) || (val < lower_values[i]) || (val > upper_values[i]) )
          val = starting_value[i];

        // Sanity-clamp Crystal Ball power-law params that can drift high
        switch( m_skew_type )
        {
          case PeakDef::NumSkewType:
            assert( 0 );
            // Fall through to NoSkew for non-debug builds
          case PeakDef::NoSkew:   case PeakDef::Bortel:
          case PeakDef::DoubleBortel: case PeakDef::GaussPlusBortel:
          case PeakDef::GaussExp: case PeakDef::ExpGaussExp:
            break;

          case PeakDef::CrystalBall:
          case PeakDef::DoubleSidedCrystalBall:
          {
            switch( ct )
            {
              case PeakDef::Mean:           case PeakDef::Sigma:
              case PeakDef::GaussAmplitude: case PeakDef::NumCoefficientTypes:
              case PeakDef::Chi2DOF:
              case PeakDef::SkewPar0:
              case PeakDef::SkewPar2:
                if( (val > 3.0) && fit_parameter[i] )
                  val = starting_value[i];
                break;
              case PeakDef::SkewPar1:
              case PeakDef::SkewPar3:
                if( (val > 6.0) && fit_parameter[i] )
                  val = starting_value[i];
                break;
            }
            break;
          }

          case PeakDef::VoigtPlusBortel:
            break;
        }//switch( m_skew_type )

        starting_value[i] = val;
      }
    }//for( const auto &p : inpeaks )

    // For small/medium refinement, restrict skew range near starting value to prevent large shifts;
    // otherwise allow the full parameter range so fitting from default starting values can succeed.
    const bool restrict_skew_range = m_options.testFlag( PeakFitLM::PeakFitLMOptions::SmallRefinementOnly )
                                     || m_options.testFlag( PeakFitLM::PeakFitLMOptions::MediumRefinementOnly );

    // Set up the base (lower-anchor) skew parameters at params[0..num_skew_pars-1]
    for( size_t skew_index = 0; skew_index < num_skew_pars; ++skew_index )
    {
      pars[skew_index] = starting_value[skew_index];

      if( fit_parameter[skew_index] )
      {
        if( restrict_skew_range )
        {
          lower_bounds[skew_index] = std::max( lower_values[skew_index], 0.5*starting_value[skew_index] );
          upper_bounds[skew_index] = std::min( upper_values[skew_index], 1.5*fabs(starting_value[skew_index]) );
        }else
        {
          lower_bounds[skew_index] = lower_values[skew_index];
          upper_bounds[skew_index] = upper_values[skew_index];
        }
      }
      else
      {
        constant_parameters.push_back( static_cast<int>(skew_index) );
      }
    }

    // For energy-dependent multi-ROI fitting, set up upper-anchor params
    if( m_fit_skew_energy_dependence )
    {
      size_t upper_idx = num_skew_pars;
      for( size_t i = 0; i < num_skew_pars; ++i )
      {
        const auto ct = PeakDef::CoefficientType( static_cast<int>(PeakDef::SkewPar0) + static_cast<int>(i) );
        if( !PeakDef::is_energy_dependent( m_skew_type, ct ) )
          continue;

        // Initialize upper-anchor to same value as lower-anchor
        pars[upper_idx] = starting_value[i];
        if( fit_parameter[i] )
        {
          lower_bounds[upper_idx] = lower_bounds[i];
          upper_bounds[upper_idx] = upper_bounds[i];
        }
        else
        {
          constant_parameters.push_back( static_cast<int>(upper_idx) );
        }
        upper_idx += 1;
      }
    }//if( m_fit_skew_energy_dependence )
  }//setup_skew_parameters(...)


  /** Set up parameters for a single ROI at param_offset in the global parameter array. */
  void setup_roi_parameters( const RoiInfo &roi,
                             const size_t param_offset,
                             double *pars,
                             vector<int> &constant_parameters,
                             vector<std::optional<double>> &lower_bounds,
                             vector<std::optional<double>> &upper_bounds ) const
  {
    const size_t num_fit_cont = roi.use_lls_for_cont
                                ? size_t(0) : PeakContinuum::num_parameters( roi.offset_type );
    const size_t num_sigmas_fit = roi_sigma_parameter_count( roi );
    const double range = roi.upper_energy - roi.lower_energy;

    // Continuum parameters (if not LLS)
    if( !roi.use_lls_for_cont )
    {
      const shared_ptr<const PeakContinuum> initial_continuum = roi.peaks.front()->continuum();
      assert( initial_continuum );
      const vector<double> &cont_pars    = initial_continuum->parameters();
      const vector<bool>    par_fit_for  = initial_continuum->fitForParameter();
      assert( cont_pars.size() == num_fit_cont );

      for( size_t i = 0; i < num_fit_cont; ++i )
      {
        pars[param_offset + i] = cont_pars[i];
        if( !par_fit_for[i] )
          constant_parameters.push_back( static_cast<int>(param_offset + i) );
      }
    }

    // Compute mean range bin width for sigma constraints
    const size_t midbin = m_data->find_gamma_channel( static_cast<float>(0.5*(roi.lower_energy + roi.upper_energy)) );
    const float binwidth = m_data->gamma_channel_width( midbin );
    const double avrg_bin_width = [&](){
      const double left_energy  = m_data->gamma_channel_lower( roi.lower_channel );
      const double right_energy = m_data->gamma_channel_upper( roi.upper_channel );
      return (right_energy - left_energy) / static_cast<double>( 1 + roi.upper_channel - roi.lower_channel );
    }();

    double minsigma = binwidth, max_input_sigma = binwidth, min_input_sigma = binwidth;
    double maxsigma = 0.5*range;

    // Mean parameters and compute minsigma/maxsigma
    for( size_t i = 0; i < roi.peaks.size(); ++i )
    {
      const std::shared_ptr<const PeakDef> &peak = roi.peaks[i];
      const double mean  = peak->mean();
      const double sigma = peak->sigma();
      const size_t mean_par_index = param_offset + num_fit_cont + num_sigmas_fit + i;

      const double rel_mean = 0.5 + (mean - roi.lower_energy) / range;
      pars[mean_par_index] = rel_mean;

      if( !peak->fitFor( PeakDef::Mean ) )
      {
        constant_parameters.push_back( static_cast<int>(mean_par_index) );
      }
      else
      {
        lower_bounds[mean_par_index] = 0.5;
        upper_bounds[mean_par_index] = 1.5;

        if( m_options.testFlag( PeakFitLM::PeakFitLMOptions::MediumRefinementOnly ) )
        {
          lower_bounds[mean_par_index] = 0.5 + (mean - 0.5*sigma - roi.lower_energy) / range;
          upper_bounds[mean_par_index] = 0.5 + (mean + 0.5*sigma - roi.lower_energy) / range;
          assert( rel_mean >= *lower_bounds[mean_par_index] );
          assert( rel_mean <= *upper_bounds[mean_par_index] );
        }

        if( m_options.testFlag( PeakFitLM::PeakFitLMOptions::SmallRefinementOnly ) )
        {
          lower_bounds[mean_par_index] = 0.5 + (mean - 0.15*sigma - roi.lower_energy) / range;
          upper_bounds[mean_par_index] = 0.5 + (mean + 0.15*sigma - roi.lower_energy) / range;
          assert( rel_mean >= *lower_bounds[mean_par_index] );
          assert( rel_mean <= *upper_bounds[mean_par_index] );
        }

        if( rel_mean <= *lower_bounds[mean_par_index] )
          lower_bounds[mean_par_index] = rel_mean - std::max( 0.1, fabs(0.25*rel_mean) );
        if( rel_mean >= *upper_bounds[mean_par_index] )
          upper_bounds[mean_par_index] = 1.2*rel_mean;

        assert( rel_mean >= *lower_bounds[mean_par_index] );
        assert( rel_mean <= *upper_bounds[mean_par_index] );
      }

      min_input_sigma = std::min( min_input_sigma, sigma );
      max_input_sigma = std::max( max_input_sigma, sigma );

      if( !peak->fitFor( PeakDef::Sigma ) )
      {
        minsigma = std::min( minsigma, sigma );
        maxsigma = std::max( maxsigma, sigma );
      }
      else
      {
        float lowersigma, uppersigma;
        expected_peak_width_limits( mean, m_isHPGe, m_data, lowersigma, uppersigma );
        if( i == 0 )
          minsigma = lowersigma;
        if( i == (roi.peaks.size() - 1) )
          maxsigma = uppersigma;
      }
    }//for each peak

    // Sigma parameters
    if( num_sigmas_fit == 0 )
    {
      // No sigma parameters to set up
    }
    else if( m_options.testFlag( PeakFitLM::PeakFitLMOptions::AllPeakFwhmIndependent ) )
    {
      // One parameter per peak with fitted sigma
      size_t fit_sigma_index = 0;
      for( const std::shared_ptr<const PeakDef> &peak : roi.peaks )
      {
        if( !peak->fitFor( PeakDef::Sigma ) )
          continue;

        const double sigma = peak->sigma();
        const size_t sigma_par_idx = param_offset + num_fit_cont + fit_sigma_index;

        pars[sigma_par_idx] = sigma / roi.max_initial_sigma;
        lower_bounds[sigma_par_idx] = (0.5*std::min(minsigma,min_input_sigma)) / roi.max_initial_sigma;
        upper_bounds[sigma_par_idx] = (1.5*std::max(maxsigma,max_input_sigma)) / roi.max_initial_sigma;

        if( m_options.testFlag( PeakFitLM::PeakFitLMOptions::MediumRefinementOnly ) )
        {
          lower_bounds[sigma_par_idx] = 0.5 * sigma / roi.max_initial_sigma;
          upper_bounds[sigma_par_idx] = 1.5 * sigma / roi.max_initial_sigma;
        }
        if( m_options.testFlag( PeakFitLM::PeakFitLMOptions::SmallRefinementOnly ) )
        {
          lower_bounds[sigma_par_idx] = 0.85 * sigma / roi.max_initial_sigma;
          upper_bounds[sigma_par_idx] = 1.15 * sigma / roi.max_initial_sigma;
        }

        assert( pars[sigma_par_idx] >= *lower_bounds[sigma_par_idx] );
        assert( pars[sigma_par_idx] <= *upper_bounds[sigma_par_idx] );
        fit_sigma_index += 1;
      }
      assert( fit_sigma_index == num_sigmas_fit );
    }
    else
    {
      // One or two shared sigma parameters.
      // The parameter represents sigma / roi.max_initial_sigma.
      // roi.max_initial_sigma has a 1.0 keV floor, so the starting value may be
      // less than 1.0 when the actual peak sigmas are smaller than 1.0 keV.
      const size_t fit_sigma_idx = param_offset + num_fit_cont;
      pars[fit_sigma_idx] = max_input_sigma / roi.max_initial_sigma;

      lower_bounds[fit_sigma_idx] = (0.5*std::min(minsigma, 0.75*min_input_sigma)) / roi.max_initial_sigma;
      upper_bounds[fit_sigma_idx] = (1.5*std::max(maxsigma, 1.25*max_input_sigma)) / roi.max_initial_sigma;

      assert( pars[fit_sigma_idx] >= *lower_bounds[fit_sigma_idx] );
      assert( pars[fit_sigma_idx] <= *upper_bounds[fit_sigma_idx] );

      if( m_options.testFlag( PeakFitLM::PeakFitLMOptions::MediumRefinementOnly ) )
      {
        lower_bounds[fit_sigma_idx] = (0.5*min_input_sigma) / roi.max_initial_sigma;
        upper_bounds[fit_sigma_idx] = (1.5*max_input_sigma) / roi.max_initial_sigma;
        assert( pars[fit_sigma_idx] >= *lower_bounds[fit_sigma_idx] );
        assert( pars[fit_sigma_idx] <= *upper_bounds[fit_sigma_idx] );
      }
      if( m_options.testFlag( PeakFitLM::PeakFitLMOptions::SmallRefinementOnly ) )
      {
        lower_bounds[fit_sigma_idx] = (0.85*min_input_sigma) / roi.max_initial_sigma;
        upper_bounds[fit_sigma_idx] = (1.15*max_input_sigma) / roi.max_initial_sigma;
        assert( pars[fit_sigma_idx] >= *lower_bounds[fit_sigma_idx] );
        assert( pars[fit_sigma_idx] <= *upper_bounds[fit_sigma_idx] );
      }

      // Enforce minimum of 0.525 channels per sigma (narrowest commercially observed detector)
      bool any_peak_fixed_narrow = false, all_widths_fixed = true;
      for( const auto &peak : roi.peaks )
      {
        const bool fit_width = peak->fitFor( PeakDef::Sigma );
        any_peak_fixed_narrow |= (!fit_width && (peak->sigma() < 0.525*avrg_bin_width));
        all_widths_fixed = (all_widths_fixed && !fit_width);
      }

      if( !any_peak_fixed_narrow )
      {
        const double rel_bin_width = avrg_bin_width / roi.max_initial_sigma;
        lower_bounds[fit_sigma_idx] = std::max( 0.525*rel_bin_width, lower_bounds[fit_sigma_idx].value() );
        upper_bounds[fit_sigma_idx] = std::max( 0.525*rel_bin_width, upper_bounds[fit_sigma_idx].value() );
        pars[fit_sigma_idx]         = std::max( 0.525*rel_bin_width, pars[fit_sigma_idx] );
      }

      assert( pars[fit_sigma_idx] >= *lower_bounds[fit_sigma_idx] );
      assert( pars[fit_sigma_idx] <= *upper_bounds[fit_sigma_idx] );
      assert( all_widths_fixed || (lower_bounds[fit_sigma_idx].has_value()
              && (lower_bounds[fit_sigma_idx].value() >= 0.525*(avrg_bin_width/roi.max_initial_sigma))) );

      if( num_sigmas_fit > 1 )
      {
        // Second sigma param = multiplier for the high-energy end; range +-20%
        const size_t upper_sigma_idx = fit_sigma_idx + 1;
        pars[upper_sigma_idx] = 1.0;
        lower_bounds[upper_sigma_idx] = 0.8;
        upper_bounds[upper_sigma_idx] = 1.2;
        assert( pars[upper_sigma_idx] >= *lower_bounds[upper_sigma_idx] );
        assert( pars[upper_sigma_idx] <= *upper_bounds[upper_sigma_idx] );
      }
    }//sigma parameter setup

    // Per-ROI skew parameters (only when IndependentSkewValues option is set)
    if( m_options.testFlag( PeakFitLM::PeakFitLMOptions::IndependentSkewValues ) )
    {
      const size_t num_skew_pars = PeakDef::num_skew_parameters( m_skew_type );
      if( num_skew_pars > 0 )
      {
        // Determine default ranges and starting values for each skew parameter
        vector<bool>   fit_parameter( num_skew_pars, false ); // OR'd across matching peaks below
        vector<double> starting_value( num_skew_pars, 0.0 );
        vector<double> lower_values( num_skew_pars, 0.0 );
        vector<double> upper_values( num_skew_pars, 0.0 );

        for( size_t i = 0; i < num_skew_pars; ++i )
        {
          const auto ct = PeakDef::CoefficientType( static_cast<int>(PeakDef::SkewPar0) + static_cast<int>(i) );
          double lower, upper, start, dx;
          const bool use = PeakDef::skew_parameter_range( m_skew_type, ct, lower, upper, start, dx );
          assert( use );
          if( !use )
            throw logic_error( "Inconsistent skew par val (IndependentSkewValues)" );
          starting_value[i] = start;
          lower_values[i]   = lower;
          upper_values[i]   = upper;
        }

        // A skew parameter is fitted if ANY matching-type peak in the ROI has fitFor==true for it;
        // it is held constant only when ALL matching peaks agree it should be fixed.
        // fit_parameter[] starts false; we OR in each matching peak's fitFor flag.
        // Initial coefficient values come from the first matching peak.
        bool found_first = false;
        for( const auto &p : roi.peaks )
        {
          if( p->skewType() != m_skew_type )
            continue;

          for( size_t i = 0; i < num_skew_pars; ++i )
          {
            const auto ct = PeakDef::CoefficientType( static_cast<int>(PeakDef::SkewPar0) + static_cast<int>(i) );

            // OR across peaks: mark as fit if any peak requests it
            if( p->fitFor( ct ) )
              fit_parameter[i] = true;

            if( !found_first )
            {
              // Coefficient starting value: use first matching peak only
              double val = p->coefficient( ct );

              if( IsInf(val) || IsNan(val) || (val < lower_values[i]) || (val > upper_values[i]) )
                val = starting_value[i];

              // Sanity-clamp Crystal Ball power-law params that can drift high
              switch( m_skew_type )
              {
                case PeakDef::NumSkewType:
                  assert( 0 );
                  // Fall through to NoSkew for non-debug builds
                case PeakDef::NoSkew:   case PeakDef::Bortel:
                case PeakDef::DoubleBortel: case PeakDef::GaussPlusBortel:
                case PeakDef::GaussExp: case PeakDef::ExpGaussExp:
                  break;

                case PeakDef::CrystalBall:
                case PeakDef::DoubleSidedCrystalBall:
                {
                  switch( ct )
                  {
                    case PeakDef::Mean:           case PeakDef::Sigma:
                    case PeakDef::GaussAmplitude: case PeakDef::NumCoefficientTypes:
                    case PeakDef::Chi2DOF:
                    case PeakDef::SkewPar0:
                    case PeakDef::SkewPar2:
                      if( val > 3.0 )
                        val = starting_value[i];
                      break;
                    case PeakDef::SkewPar1:
                    case PeakDef::SkewPar3:
                      if( val > 6.0 )
                        val = starting_value[i];
                      break;
                  }
                  break;
                }

                case PeakDef::VoigtPlusBortel:
                  break;
              }//switch( m_skew_type )

              starting_value[i] = val;
            }//if( !found_first )
          }//for( size_t i = 0; i < num_skew_pars; ++i )

          found_first = true;
        }//for( const auto &p : roi.peaks )

        // Write the per-ROI skew params into the global parameter array
        const size_t skew_base_idx = param_offset + num_fit_cont + num_sigmas_fit + roi.peaks.size();
        const bool restrict_skew_range = m_options.testFlag( PeakFitLM::PeakFitLMOptions::SmallRefinementOnly )
                                         || m_options.testFlag( PeakFitLM::PeakFitLMOptions::MediumRefinementOnly );
        for( size_t skew_index = 0; skew_index < num_skew_pars; ++skew_index )
        {
          const size_t abs_idx = skew_base_idx + skew_index;
          pars[abs_idx] = starting_value[skew_index];

          if( fit_parameter[skew_index] )
          {
            if( restrict_skew_range )
            {
              lower_bounds[abs_idx] = std::max( lower_values[skew_index], 0.5*starting_value[skew_index] );
              upper_bounds[abs_idx] = std::min( upper_values[skew_index], 1.5*fabs(starting_value[skew_index]) );
            }else
            {
              lower_bounds[abs_idx] = lower_values[skew_index];
              upper_bounds[abs_idx] = upper_values[skew_index];
            }
          }
          else
          {
            constant_parameters.push_back( static_cast<int>(abs_idx) );
          }
        }
      }//if( num_skew_pars > 0 )
    }//if( IndependentSkewValues )
  }//setup_roi_parameters(...)


  ProblemSetup get_problem_setup() const
  {
    const size_t num_fit_pars = number_parameters();

    vector<int> constant_parameters;
    vector<double> parameters( num_fit_pars, 0.0 );
    vector<std::optional<double>> lower_bounds( num_fit_pars );
    vector<std::optional<double>> upper_bounds( num_fit_pars );

    // Shared skew parameters at the front
    setup_skew_parameters( &parameters[0], constant_parameters, lower_bounds, upper_bounds,
                           all_starting_peaks() );

    // Per-ROI parameters following the skew block
    size_t param_offset = skew_parameter_count();
    for( const RoiInfo &roi : m_rois )
    {
      setup_roi_parameters( roi, param_offset, &parameters[0], constant_parameters,
                            lower_bounds, upper_bounds );
      param_offset += roi_parameter_count( roi );
    }

    assert( param_offset == num_fit_pars );

    ProblemSetup prob_setup;
    prob_setup.m_parameters         = parameters;
    prob_setup.m_constant_parameters = constant_parameters;
    prob_setup.m_lower_bounds        = lower_bounds;
    prob_setup.m_upper_bounds        = upper_bounds;

    return prob_setup;
  }//get_problem_setup()

public:
  const std::shared_ptr<const SpecUtils::Measurement> m_data;
  const PeakDef::SkewType m_skew_type;
  const std::shared_ptr<const SpecUtils::Measurement> m_external_continuum;
  const bool m_isHPGe;
  const Wt::WFlags<PeakFitLMOptions> m_options;

  mutable std::atomic<unsigned int> m_ncalls;

  const std::vector<RoiInfo> m_rois;
  const size_t m_total_num_peaks;
  const size_t m_num_parameters;
  const size_t m_num_residuals;
  const bool m_fit_skew_energy_dependence;
  const double m_skew_anchor_lower_energy;
  const double m_skew_anchor_upper_energy;

#if( PEAK_FIT_LM_PARALLEL_ROIS )
  // Non-null only when m_rois.size() > 2; sized to min(nrois, physical_cores).
  const std::unique_ptr<boost::asio::thread_pool> m_thread_pool;
#endif
};//struct PeakFitDiffCostFunction



/** All peaks passed in must share a PeakContinuum.
 */
vector<shared_ptr<const PeakDef>> fit_peaks_in_roi_LM( const vector<shared_ptr<const PeakDef>> coFitPeaks,
                                                      const std::shared_ptr<const SpecUtils::Measurement> &dataH,
                                                      const bool isHPGe,
                                                      const Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options )
{
  /** For this first go, we will have Ceres fit for things.

   This is only tested to the "seems to work, for simple cases" level.

   In the future:
     - Need to deal with errors during solving, and during Covariance computation
     - Deal with setting number of threads reasonably.
     - Try out using a ceres::LossFunction
     - Optimize/pcik the Ceres-based paramaters to get reasonable fits.
   */

  if( coFitPeaks.empty() )
    throw runtime_error( "fit_peaks_in_roi_LM: empty input peaks." );

  // Lets make sure all coFitPeaks share a continuum.
  for( size_t i = 1; i < coFitPeaks.size(); ++i )
  {
    if( coFitPeaks[i]->continuum() != coFitPeaks[0]->continuum() )
      throw runtime_error( "fit_peak_for_user_click_LM: input peaks all must share a single continuum" );
  }//for( size_t i = 1; i < coFitPeaks.size(); ++i )

  try
  {
    // TODO: check that repeaded calls to this function wont cause adding/removing channels due to find_gamma_channel(...) rounding type things - e.g., do we need to subtract off a tiny bit from roiUpperEnergy

    const shared_ptr<const PeakContinuum> input_continuum = coFitPeaks[0]->continuum();

    const float roiLowerEnergy = static_cast<float>( input_continuum->lowerEnergy() );
    const float roiUpperEnergy = static_cast<float>( input_continuum->upperEnergy() );
    assert( roiLowerEnergy < roiUpperEnergy );
    if( roiLowerEnergy >= roiUpperEnergy )
      throw runtime_error( "Invalid energy range (" + std::to_string(roiLowerEnergy) + ", " + std::to_string(roiUpperEnergy) + ")" );


    const size_t lower_channel = dataH->find_gamma_channel(roiLowerEnergy);
    const size_t upper_channel = dataH->find_gamma_channel(roiUpperEnergy);
    const double roi_lower = dataH->gamma_channel_lower( lower_channel );
    const double roi_upper = dataH->gamma_channel_upper( upper_channel );
    
    assert( (lower_channel < 2) || (roi_lower <= roiLowerEnergy) );
    assert( (roi_upper >= roiUpperEnergy) || ((upper_channel + 3) >= dataH->num_gamma_channels()) );

    if( coFitPeaks.size() >= (upper_channel - lower_channel) ) //PeakFitDiffCostFunction will check for this in an assert
      throw runtime_error( "Invalid energy channels (" + std::to_string(lower_channel) + ", "
                          + std::to_string(upper_channel) + ") for " + std::to_string(coFitPeaks.size()) + " peaks." );
    assert( coFitPeaks.size() < (upper_channel - lower_channel) );

    const size_t nchannels = dataH->num_gamma_channels();
    const size_t midbin = dataH->find_gamma_channel( 0.5*(roiLowerEnergy + roiUpperEnergy) );
    const float binwidth = dataH->gamma_channel_width( midbin );
    const size_t num_fit_peaks = coFitPeaks.size() + 1;

    const PeakContinuum::OffsetType offset_type = input_continuum->type();
    const size_t num_continuum_pars = PeakContinuum::num_parameters( offset_type );

    const double prev_ref_energy = input_continuum->referenceEnergy();
    const bool prev_ref_energy_valid = ((prev_ref_energy >= roiLowerEnergy) && (prev_ref_energy <= roiUpperEnergy));
    const double reference_energy = prev_ref_energy_valid ? prev_ref_energy : roiLowerEnergy;

    const vector<bool> parFitFors = input_continuum->fitForParameter();
    assert( parFitFors.size() <= num_continuum_pars );

    const PeakDef::SkewType skew_type = PeakFitDiffCostFunction::skew_type_from_prev_peaks( coFitPeaks );
    const size_t num_skew = PeakDef::num_skew_parameters( skew_type );


    auto cost_functor = make_unique<PeakFitDiffCostFunction>( dataH, coFitPeaks, roiLowerEnergy, roiUpperEnergy,
                                                    reference_energy, skew_type, isHPGe, fit_options );

    //Choosing 8 paramaters to include in the `ceres::Jet<>` is 4 peaks in ROI, which covers most cases
    //  without introducing a ton of extra overhead.
    auto cost_function = new ceres::DynamicAutoDiffCostFunction<PeakFitDiffCostFunction,8>( cost_functor.get(), ceres::Ownership::DO_NOT_TAKE_OWNERSHIP );

    const size_t num_fit_pars = cost_functor->number_parameters();

    cost_function->AddParameterBlock( static_cast<int>(num_fit_pars) );
    cost_function->SetNumResiduals( static_cast<int>(cost_functor->number_residuals()) );

    const PeakFitDiffCostFunction::ProblemSetup prob_setup = cost_functor->get_problem_setup();

    assert( prob_setup.m_parameters.size() == num_fit_pars );
    assert( prob_setup.m_lower_bounds.size() == num_fit_pars );
    assert( prob_setup.m_upper_bounds.size() == num_fit_pars );

    vector<double> parameters = prob_setup.m_parameters;
    double * const pars = &parameters[0];


    std::unique_ptr<ceres::Problem> problem = make_unique<ceres::Problem>();

    // A brief look at a dataset of ~4k HPGe spectra with known truth-value peaks areas shows
    //  that using a loss function doesnt seem improve outcomes, when measured by sucessful
    //  peak fits, and by comparison of fit to truth peak areas.
    ceres::LossFunction *lossfcn = nullptr;

    problem->AddResidualBlock( cost_function, lossfcn, pars ); //Note: problem takes ownership of `cost_function`


    if( !prob_setup.m_constant_parameters.empty() )
    {
      ceres::Manifold *subset_manifold = new ceres::SubsetManifold( static_cast<int>(num_fit_pars), prob_setup.m_constant_parameters );
      problem->SetManifold( pars, subset_manifold );
    }

    for( size_t i = 0; i < num_fit_pars; ++i )
    {
      if( prob_setup.m_lower_bounds[i].has_value() )
      {
        assert( *prob_setup.m_lower_bounds[i] <= parameters[i] );
        problem->SetParameterLowerBound(pars, static_cast<int>(i), *prob_setup.m_lower_bounds[i] );
      }

      if( prob_setup.m_upper_bounds[i].has_value() )
      {
        assert( *prob_setup.m_upper_bounds[i] >= parameters[i] );
        problem->SetParameterUpperBound(pars, static_cast<int>(i), *prob_setup.m_upper_bounds[i] );
      }
    }//for( size_t i = 0; i < num_fit_pars; ++i )


    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_type = ceres::TRUST_REGION; //ceres::LINE_SEARCH
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT; //ceres::DOGLEG
    options.use_nonmonotonic_steps = true;
    options.max_consecutive_nonmonotonic_steps = 10;

    //options.max_num_consecutive_invalid_steps = 10; //default 5

    // Initial trust region was very coursely optimized using a fit over ~4k HPGe spectra, using the median peak area
    //  error from truth value (whose value was 1.65*sqrt(area)); but also the average error was also about optimal
    //  (value 3.12*sqrt(area)).
    options.initial_trust_region_radius = 35*(cost_functor->number_parameters() - prob_setup.m_constant_parameters.size());
    //options.initial_trust_region_radius = 1e4; //
    options.max_trust_region_radius = 1e16;

    // Minimizer terminates when the trust region radius becomes smaller than this value.
    options.min_trust_region_radius = 1e-32;
    // Lower bound for the relative decrease before a step is accepted.
    options.min_relative_decrease = 1e-3; //pseudo optimized based on success rate of fitting peaks - but unknown effect on accuracy of fits.

    // It looks like it usually takes well below 100 calls to solve most problems, but a tail up to ~500, and then rare
    //  (pathelogical?) cases that can take >50k
    //  However, for these pathological cases, if we just try again starting from where it left off, it will solve
    //  it super quick - I guess because there will be more "momentum" to find the minimum, or something, not totally
    //  sure.
    //  We could probably reduce this number of iterations to ~100 - but the effects or impact of this hasnt been
    //  evaluated.
    //  Also, have not checked dependence on number of peaks at all
    options.max_num_iterations = (coFitPeaks.size() > 2) ? 1000 : 500;
#ifndef NDEBUG
    options.max_solver_time_in_seconds = (coFitPeaks.size() > 2) ? 180.0 : 120;
#else
    options.max_solver_time_in_seconds = (coFitPeaks.size() > 2) ? 30.0 : 150.0;
#endif

#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
    options.minimizer_progress_to_stdout = true;
    options.logging_type = ceres::PER_MINIMIZER_ITERATION;
#else
    options.minimizer_progress_to_stdout = false;
    options.logging_type = ceres::SILENT;
#endif

    /** Termination when `(new_cost - old_cost) < function_tolerance * old_cost`
     Default 1e-9;  Setting this to 1e-6 from 1e-9 is a ~80x speedup, and slgiht success fitting more peaks.
     Values from 1E-5 to 1E-9 dont seem to have a big effect on accuracy, but 1E-7 had best results for a HPGe dataset

     Value      AvrgAccuracyNSigma   MedianAccuracyNSigma   CPU-Time
     1e-5      3.25068                           1.66096                             465.296 s
     1e-6      3.20428                           1.69388                             581.317 s
     1e-7      3.16456                           1.6934                               1864.5 s
     1e-8      3.29325                           1.6936                               16833.9 s
     */
    options.function_tolerance = 1e-7;
    options.parameter_tolerance = 1e-11; //Default value is 1e-8.  Using 1e-11, so its usually the function tolerance that terminates things.
    options.num_threads = 1; //Probably wont have much/any effect

    // Default value of `max_num_consecutive_invalid_steps` is 5, however, if we are re-fitting a peak who already has
    //  near perfect values, occasionally the fit fails with the devault values - I guess the trust-region is just so
    //  far off (but I dont really know - this is just a guess), but if we increase this value to 10, the fit seems
    //  to be sucessful, for a limited number of test cases.
    options.max_num_consecutive_invalid_steps = 10;

    // Just to note for the future, the following values are enabled by the default options, and are what we want.
    //options.preconditioner_type = ceres::JACOBI;
    //options.jacobi_scaling = true; // Let Ceres estimate scale based on Jacobian diagonal

    ceres::Solver::Summary summary;
    ceres::Solve(options, problem.get(), &summary);

#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
    //std::cout << summary.BriefReport() << "\n";
    std::cout << summary.FullReport() << "\n";
    cout << "Took " << cost_functor->m_ncalls.load() << " calls to solve." << endl;
#endif

    string failure_reason;

    switch( summary.termination_type )
    {
      case ceres::CONVERGENCE:
      case ceres::USER_SUCCESS:
        break;
        
      case ceres::NO_CONVERGENCE:
      {
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
        cerr << "The L-M ceres::Solver solving failed - NO_CONVERGENCE:\n" << summary.FullReport() << endl;
#endif
        // We will give it another go - for most cases it seems the solution will now be succesffuly found really
        //  quickly.  The guess is this is from momentum, or whatever, but not really sure.
        //cerr << "Initial L-M ceres::Solver failed with " << cost_functor->m_ncalls.load() << " calls and "
        //  << summary.num_successful_steps << " steps - will try again" << endl;

        options.max_num_iterations = 5000;
        // TODO: see if we should loosed any of the following up.
        //options.function_tolerance = 1e-6;
        //options.parameter_tolerance = 1e-9;
        //options.initial_trust_region_radius = 35*(cost_functor->number_parameters() - prob_setup.m_constant_parameters.size());
        cost_functor->m_ncalls = 0;

        summary = ceres::Solver::Summary();

        ceres::Solve(options, problem.get(), &summary);

        // If things worked out, we're good here and we can `break` from this case statement
        if( (summary.termination_type == ceres::CONVERGENCE)
            || (summary.termination_type == ceres::USER_SUCCESS) )
        {
          break;
        }

#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
        cerr << "Retry of L-M ceres::Solver Failed with " << cost_functor->m_ncalls.load() << " additional calls" << endl;
#endif

        // If we have failed here, we will re-try, but with numerical differentiation - havent explicitly found any
        //  cases where this is necassary, or helpful
        failure_reason = "The L-M ceres::Solver solving failed - NO_CONVERGENCE.";

        [[fallthrough]]; //Fallthrough intentional
      }//case ceres::NO_CONVERGENCE:

      case ceres::FAILURE:
      {
        if( failure_reason.empty() )
          failure_reason = "The L-M ceres::Solver solving failed - FAILURE.";

#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
        cerr << "The L-M ceres::Solver solving failed - FAILURE:\n" << summary.FullReport() << endl;
#endif
        // We will re-try with numerical differntiation - the one case that the author has observed this being
        //  necassary is with peak-refits, where all peak paramaeters are already about perfect, very rarely, it seems
        //  to throw off trust-region searching, or something - however, numerical differentiating seems to work for
        //  these (rare!) cases that this is observed.

        auto numeric_cost_fnct = new ceres::DynamicNumericDiffCostFunction<PeakFitDiffCostFunction>( cost_functor.get(),
                                                                            ceres::Ownership::DO_NOT_TAKE_OWNERSHIP );
        numeric_cost_fnct->AddParameterBlock( static_cast<int>(num_fit_pars) );
        numeric_cost_fnct->SetNumResiduals( static_cast<int>(cost_functor->number_residuals()) );

        std::unique_ptr<ceres::Problem> numerical_problem = make_unique<ceres::Problem>();
        numerical_problem->AddResidualBlock( numeric_cost_fnct, lossfcn, pars ); //Note: numerical_problem takes ownership of `numeric_cost_fnct`

        if( !prob_setup.m_constant_parameters.empty() )
        {
          ceres::Manifold *subset_manifold = new ceres::SubsetManifold( static_cast<int>(num_fit_pars), prob_setup.m_constant_parameters );
          numerical_problem->SetManifold( pars, subset_manifold );
        }

        for( size_t i = 0; i < num_fit_pars; ++i )
        {
          if( prob_setup.m_lower_bounds[i].has_value() )
            numerical_problem->SetParameterLowerBound(pars, static_cast<int>(i), *prob_setup.m_lower_bounds[i] );
          if( prob_setup.m_upper_bounds[i].has_value() )
            numerical_problem->SetParameterUpperBound(pars, static_cast<int>(i), *prob_setup.m_upper_bounds[i] );
        }//for( size_t i = 0; i < num_fit_pars; ++i )

        cost_functor->m_ncalls = 0;
        summary = ceres::Solver::Summary();
        ceres::Solve(options, numerical_problem.get(), &summary);

        switch( summary.termination_type )
        {
          case ceres::CONVERGENCE:
          case ceres::USER_SUCCESS:
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
            cout << "Was able to solved problem using numerical differentiation!" << endl;
#endif
            problem = std::move(numerical_problem);
            break;

          case ceres::NO_CONVERGENCE:
          case ceres::FAILURE:
          case ceres::USER_FAILURE:
          {
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
            cerr << "Failed to get successful fit with numerical differentiation either - giving up." << endl;
#endif
            throw runtime_error( failure_reason );
          }
        }

        break;
      }//case ceres::FAILURE:

      case ceres::USER_FAILURE:
      {
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
        cerr << "The L-M ceres::Solver solving failed - USER_FAILURE:\n" << summary.FullReport() << endl;
#endif
        throw runtime_error( "The L-M ceres::Solver solving failed - USER_FAILURE." );
      }//case ceres::USER_FAILURE:
    }//switch( summary.termination_type )

    cost_functor->m_ncalls = 0;

    ceres::Covariance::Options cov_options;
    cov_options.algorithm_type = ceres::CovarianceAlgorithmType::DENSE_SVD; //SPARSE_QR;

    // Some terms we are fitting may not matter much, so when computing the inverse of J'J, we will opt
    //  to drop all terms where the eigen value of that term, divided by the maximum eigenvalue is
    //  less than the condition number.
    //  TODO: Currently leaving condition number at default 1e-14 - but what value is reasonable should be investigated
    cov_options.null_space_rank = -1;
    cov_options.min_reciprocal_condition_number = 1e-14;

    vector<double> uncertainties( num_fit_pars, 0.0 );
    double *uncertainties_ptr = uncertainties.data();
    // Row-major covariance matrix (num_fit_pars x num_fit_pars); valid only when uncertainties_ptr != nullptr.
    vector<double> row_major_covariance;
    
    ceres::Covariance covariance(cov_options);
    vector<pair<const double*, const double*> > covariance_blocks;
    covariance_blocks.push_back( make_pair( pars, pars) );
    
    if( !covariance.Compute(covariance_blocks, problem.get()) )
    {
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
      cerr << "Failed to compute covariance!" << endl;
#endif
      uncertainties_ptr = nullptr;
    }else
    {
      row_major_covariance.resize( num_fit_pars * num_fit_pars );
      const vector<const double *> const_par_blocks( 1, pars );

      const bool success = covariance.GetCovarianceMatrix( const_par_blocks, row_major_covariance.data() );
      assert( success );
      if( success )
      {
        for( size_t i = 0; i < num_fit_pars; ++i )
        {
          // TODO: compare uncertainties with Minuit method - should maybe check out.
          if( row_major_covariance[i*num_fit_pars + i] > 0.0 )
            uncertainties[i] = sqrt( row_major_covariance[i*num_fit_pars + i] );
        }
      }else 
      {
        uncertainties_ptr = nullptr;
        row_major_covariance.clear();
      }//
    }//if( we failed to computer covariance ) / else

#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
    // Using numeric differentiation, should only be a single call to get covariance.
    cout << "Took " << cost_functor->m_ncalls.load() << " calls to get covariance." << endl;
#endif

    const double *cov_ptr = row_major_covariance.empty() ? nullptr : row_major_covariance.data();
    vector<double> residuals( cost_functor->number_residuals(), 0.0 );
    auto final_peaks = cost_functor->parametersToPeaks<PeakDef,double>( &parameters[0], uncertainties_ptr,
                                                                         residuals.data(), cov_ptr, num_fit_pars );

    vector<shared_ptr<const PeakDef>> results( final_peaks.size() );
    for( size_t i = 0; i < final_peaks.size(); ++i )
      results[i] = make_shared<PeakDef>( final_peaks[i] );

    return results;
  }catch( std::exception &e )
  {
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
    cout << "fit_peak_for_user_click_LM caught: " << e.what() << endl;
#endif
    //assert( 0 );
    throw;
  }//try / catch

  assert( 0 );
  throw runtime_error( "shouldnt have gotten here" );
  return {};
}//void fit_peaks_in_roi_LM(...)



void fit_peak_for_user_click_LM( PeakShrdVec &results,
                             const std::shared_ptr<const SpecUtils::Measurement> &dataH,
                             const vector<shared_ptr<const PeakDef>> &coFitPeaksInput,
                             const double mean0, const double sigma0,
                             const double area0,
                             const float roiLowerEnergy,
                             const float roiUpperEnergy,
                             const bool isHPGe )
{
  vector<shared_ptr<const PeakDef>> coFitPeaks = coFitPeaksInput;
  
  const PeakDef::SkewType skew_type = PeakFitDiffCostFunction::skew_type_from_prev_peaks( coFitPeaks );
  const size_t num_skew = PeakDef::num_skew_parameters( skew_type );

  std::shared_ptr<PeakDef> candidatepeak = std::make_shared<PeakDef>(mean0, sigma0, area0);

  //The below should probably go off the number of bins in the ROI
  const size_t num_fit_peaks = coFitPeaksInput.size() + 1;
  PeakContinuum::OffsetType offset_type = (num_fit_peaks < (isHPGe ? 3 : 2)) ? PeakContinuum::Linear : PeakContinuum::Quadratic;

  if( coFitPeaks.size() )
    offset_type = std::max( offset_type, coFitPeaks[0]->continuum()->type() );


  // Make sure all initial peaks share a continuum
  if( coFitPeaks.size() )
  {
    shared_ptr<PeakContinuum> initial_continuum = make_shared<PeakContinuum>( *coFitPeaks[0]->continuum() );
    initial_continuum->setRange( roiLowerEnergy, roiUpperEnergy );
    candidatepeak->setContinuum( initial_continuum );
    for( size_t i = 0; i < coFitPeaks.size(); ++i )
    {
      auto newpeak = make_shared<PeakDef>( *coFitPeaks[i] );
      newpeak->setContinuum( initial_continuum );
      coFitPeaks[i] = newpeak;
    }

    if( initial_continuum->type() != offset_type )
    {
      const vector<bool> parFitFors = initial_continuum->fitForParameter();
      assert( parFitFors.size() <= PeakContinuum::num_parameters( initial_continuum->type() ) );

      bool any_continuum_par_fixed = false;
      for( size_t fit_for_index = 0; fit_for_index < parFitFors.size(); ++fit_for_index )
        any_continuum_par_fixed |= !parFitFors[fit_for_index];

      if( !any_continuum_par_fixed )
      {
        const double prev_ref_energy = initial_continuum->referenceEnergy();
        const bool prev_ref_energy_valid = ((prev_ref_energy >= roiLowerEnergy) && (prev_ref_energy <= roiUpperEnergy));
        const double reference_energy = prev_ref_energy_valid ? prev_ref_energy : roiLowerEnergy;

        initial_continuum->calc_linear_continuum_eqn( dataH, reference_energy, roiLowerEnergy, roiUpperEnergy, 2, 2 );
      }
      initial_continuum->setType( offset_type );
    }
  }else
  {
    shared_ptr<PeakContinuum> initial_continuum = candidatepeak->continuum();
    initial_continuum->setRange( roiLowerEnergy, roiUpperEnergy );
    const double reference_energy = roiLowerEnergy;
    initial_continuum->calc_linear_continuum_eqn( dataH, reference_energy, roiLowerEnergy, roiUpperEnergy, 2, 2 );
    initial_continuum->setType( offset_type );
  }//if( coFitPeaks.size() )

  coFitPeaks.push_back( candidatepeak );
  std::sort( coFitPeaks.begin(), coFitPeaks.end(), &PeakDef::lessThanByMeanShrdPtr );


  const Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options( 0 );

  try
  {
    results = fit_peaks_in_roi_LM( coFitPeaks, dataH, isHPGe, fit_options );
  }catch( std::exception &e )
  {
    results.clear();
  }
}//void fit_peak_for_user_click_LM(...)


void fit_peaks_LM( vector<shared_ptr<const PeakDef>> &results,
                  const vector<shared_ptr<const PeakDef>> input_peaks,
                  shared_ptr<const SpecUtils::Measurement> data,
                  const double stat_threshold,
                  const double hypothesis_threshold,
                  const bool is_refit,
                  const bool isHPGe ) throw()
{
  try
  {
    // Check all input peaks share a ROI
    for( size_t i = 1; i < input_peaks.size(); ++i )
    {
      if( input_peaks[i]->continuum() != input_peaks[0]->continuum() )
        throw runtime_error( "fit_peaks_LM: all input peaks must share ROI." );
    }

    results.clear();

    //We have to separate out non-gaussian peaks since they cant enter the
    //  fitting methods
    vector<shared_ptr<const PeakDef>> near_peaks, datadefined_peaks;

    shared_ptr<const SpecUtils::Measurement> ext_continuum;

    near_peaks.reserve( input_peaks.size() );

    for( const shared_ptr<const PeakDef> &p : input_peaks )
    {
      if( p->gausPeak() )
        near_peaks.push_back( p );
      else
        datadefined_peaks.push_back( p );

      if( !!p->continuum()->externalContinuum() )
        ext_continuum = p->continuum()->externalContinuum();
    }//for( const PeakDef &p : all_near_peaks )


    if( near_peaks.empty() )
      return;

    // Completely replace contents of `near_peaks` with new peaks, and new `PeakContinuum`
    //  so we dont affect the input peaks.
    local_unique_copy_continuum( near_peaks );

    //Need to make sure near_peaks and fixedpeaks are all gaussian (if not
    //  seperate them out, and add them in later).  If fitpeaks is non-gaussian
    //  ignore it or throw an exception.

    double lowx( 0.0 ), highx( 0.0 );

    {
      const shared_ptr<const PeakDef> &lowgaus = near_peaks.front();
      const shared_ptr<const PeakDef> &highgaus = near_peaks.back();

      double dummy = 0.0;
      findROIEnergyLimits( lowx, dummy, *lowgaus, data, isHPGe );
      findROIEnergyLimits( dummy, highx, *highgaus, data, isHPGe );
    }


    Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options( 0 );
    if( is_refit )
      fit_options |= PeakFitLM::PeakFitLMOptions::MediumRefinementOnly;

    results = fit_peaks_in_roi_LM( near_peaks, data, isHPGe, fit_options );

    // I'm pretty sure things have been sorted, but lets check
    assert( std::is_sorted(begin(results), end(results), &PeakDef::lessThanByMeanShrdPtr ) );

    //
    //set_chi2_dof( data, results, 0, near_peaks.size() );


    for( size_t i = 1; i <= results.size(); ++i ) //Note weird convntion of index
    {
      const shared_ptr<const PeakDef> &peak = (results[i-1]);

      // Dont remove peaks whos amplitudes we arent fitting
      if( !peak->fitFor(PeakDef::GaussAmplitude) )
        continue;

      // Dont enforce a significance test if we are refitting the peak, and
      //  Sigma and Mean are fixed - the user probably knows what they
      //  are are doing.
      if( is_refit
         && !peak->fitFor(PeakDef::Sigma)
         && !peak->fitFor(PeakDef::Mean) )
      {
        continue;
      }

      const double num_sigma = (peak->amplitude() / peak->amplitudeUncert());
      const bool is_sig = (stat_threshold <= 0.0) || (num_sigma >= stat_threshold);

      const double dummy_stat_thresh = 0.0;
      const bool significant = chi2_significance_test( *peak, dummy_stat_thresh, hypothesis_threshold, {}, data );
      if( !is_sig || !significant )
      {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
        DebugLog(cerr) << "\tPeak at mean=" << peak->mean()
        << "is being discarded for not being significant"
        << "\n";
#endif
        results.erase( results.begin() + --i );
      }//if( !significant )
    }//for( size_t i = 1; i < fitpeaks.size(); ++i )

    bool removed_peak = false;
    for( size_t i = 1; i < results.size(); ++i ) //Note weird convention of index
    {
      const shared_ptr<const PeakDef> &this_peak = results[i-1];
      const shared_ptr<const PeakDef> &next_peak = results[i+1-1];

      // Dont remove peaks whos amplitudes we arent fitting
      if( !this_peak->fitFor(PeakDef::GaussAmplitude) )
        continue;

      const double min_sigma = min( this_peak->sigma(), next_peak->sigma() );
      const double mean_diff = next_peak->mean() - this_peak->mean();

      //In order to remove a gaussian, the peaks must both be within a sigma
      //  of eachother.  Note that this proccess doesnt care about the widths
      //  of the gaussians because we are assuming that the width of the gaussian
      //  should only be dependant on energy, so should only have one width of
      //  gaussian for a given energy
      if( (mean_diff/min_sigma) < 1.0 ) //XXX 1.0 chosen arbitrarily, and not checked
      {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
        DebugLog(cerr) << "Removing duplicate peak at x=" << this_peak->mean() << " sigma="
        << this_peak->sigma() << " in favor of mean=" << next_peak->mean()
        << " sigma=" << next_peak->sigma() << "\n";
#endif

        removed_peak = true;

        //Delete the peak with the worst chi2
        if( this_peak->chi2dof() > next_peak->chi2dof() )
          results.erase( begin(results) + i - 1 );
        else
          results.erase( begin(results) + i );

        i = i - 1; //incase we have multiple close peaks in a row
      }//if( (mean_diff/min_sigma) < 1.0 ) / else
    }//for( size_t i = 1; i < fitpeaks.size(); ++i )

    if( removed_peak )
      results = fit_peaks_in_roi_LM( results, data, isHPGe, fit_options );

    if( datadefined_peaks.size() )
    {
      results.insert( end(results), begin(datadefined_peaks), end(datadefined_peaks) );
      std::sort( begin(results), end(results), &PeakDef::lessThanByMeanShrdPtr );
    }

    return;
  }catch( std::exception &e )
  {
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO || !defined(NDEBUG) )
    cerr << "fit_peaks_LM: caught exception '" << e.what() << "'" << endl;
#endif
    results.clear();
    return;
  }//try / catch


  //We will only reach here if there was no exception, so since never expect
  //  this to actually happen, just assign the results to be same as the input
  assert( 0 );
  results = input_peaks;
}//vector<PeakDef> fit_peaks_LM(...);


vector<shared_ptr<const PeakDef>> fit_peaks_in_range_LM( const double x0, const double x1,
                                      const double ncausalitysigma,
                                      const double stat_threshold,
                                      const double hypothesis_threshold,
                                      const std::vector<std::shared_ptr<const PeakDef>> input_peaks,
                                      const std::shared_ptr<const SpecUtils::Measurement> data,
                                      const bool isRefit,
                                      const bool isHPGe )
{
  if( !data || (x1 < x0) )
    return input_peaks;

  vector<shared_ptr<const PeakDef>> all_peaks = input_peaks;
  std::sort( begin(all_peaks), end(all_peaks), &PeakDef::lessThanByMeanShrdPtr );

  const vector<shared_ptr<const PeakDef>> peaks_in_range = peaksInRange( x0, x1, ncausalitysigma, all_peaks );
  vector<shared_ptr<const PeakDef>> peaks_not_in_range = all_peaks;

  peaks_not_in_range.erase( std::remove_if( begin(peaks_not_in_range), end(peaks_not_in_range),
    [&peaks_in_range]( const shared_ptr<const PeakDef> &p ) -> bool {
      // Return false if the peak is in the energy range of interest
      return std::find(begin(peaks_in_range), end(peaks_in_range), p) != end(peaks_in_range);
    } ),
    end(peaks_not_in_range) );

  assert( (peaks_not_in_range.size() + peaks_in_range.size()) == input_peaks.size() );


  // The returned pointers all point to the original peaks.
  const vector<vector<shared_ptr<const PeakDef>>> seperated_peaks
                                          = causilyDisconnectedPeaks( ncausalitysigma, false, peaks_in_range );

  assert( ([seperated_peaks](){ size_t i = 0; for(auto &p:seperated_peaks) i += p.size(); return i; })() == peaks_in_range.size() );

  //Fit each of the ranges
  vector<vector<shared_ptr<const PeakDef>>> fit_peak_ranges( seperated_peaks.size() );
  SpecUtilsAsync::ThreadPool threadpool;
  for( size_t peakn = 0; peakn < seperated_peaks.size(); ++peakn )
  {
    threadpool.post( boost::bind( &fit_peaks_LM,
                                 boost::ref(fit_peak_ranges[peakn]),
                                 boost::cref(seperated_peaks[peakn]),
                                 data,
                                 stat_threshold,
                                 hypothesis_threshold,
                                 isRefit,
                                 isHPGe ) );
  }//for( size_t peakn = 0; peakn < seperated_peaks.size(); ++peakn )
  threadpool.join();

  vector<shared_ptr<const PeakDef>> results = peaks_not_in_range;

  //put the fit peaks back into 'input_peaks' so we can return all the peaks
  //  passed in, not just ones in X range of interest
  for( size_t peakn = 0; peakn < fit_peak_ranges.size(); ++peakn )
  {
    const vector<shared_ptr<const PeakDef>> &fit_peaks = fit_peak_ranges[peakn];
    results.insert( end(results), begin(fit_peaks), end(fit_peaks) );
  }//for( size_t peakn = 0; peakn < fit_peak_ranges.size(); ++peakn )

  std::sort( begin(results), end(results), &PeakDef::lessThanByMeanShrdPtr );

  //Now make sure peaks from two previously causally disconnected regions
  //  didn't migrate towards each other, causing the regions to become
  //  causally connected now
  bool migration = false;
  for( size_t peakn = 1; peakn < fit_peak_ranges.size(); ++peakn )
  {
    if( fit_peak_ranges[peakn-1].empty() || fit_peak_ranges[peakn].empty() )
      continue;

    const shared_ptr<const PeakDef> &last_peak = fit_peak_ranges[peakn-1].back();
    const shared_ptr<const PeakDef> &this_peak = fit_peak_ranges[peakn][0];
    migration |= PeakDef::causilyConnected( *last_peak, *this_peak, ncausalitysigma, false );
  }//for( size_t peakn = 0; peakn < seperated_peaks.size(); ++peakn )

  if( migration )
  {
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
    cerr << "fitPeaksInRange(...)\n\tWarning: Migration happened!" << endl;
#endif

    return fit_peaks_in_range_LM( x0, x1, ncausalitysigma,
                           stat_threshold, hypothesis_threshold,
                                 results, data, isRefit, isHPGe );
  }//if( migration )

  //  cout << "Fit took: " << timer.format() << endl;

  return results;
}//vector<shared_ptr<const PeakDef>> fit_peaks_in_range_LM(...)

std::vector<std::shared_ptr<const PeakDef>> refitPeaksThatShareROI_LM(
                                   const std::shared_ptr<const SpecUtils::Measurement> &data,
                                   const std::shared_ptr<const DetectorPeakResponse> &detector,
                                   const std::vector<std::shared_ptr<const PeakDef>> &inpeaks,
                                   const Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options )
{
  vector<shared_ptr<const PeakDef>> answer;

  try
  {
    if( inpeaks.empty() )
      return answer;

    std::shared_ptr<const PeakContinuum> origCont = inpeaks[0]->continuum();

    for( const shared_ptr<const PeakDef> &p : inpeaks )
      if( origCont != p->continuum() )
        throw runtime_error( "refitPeaksThatShareROI_LM: all input peaks must share a ROI" );


    //const bool isHPGe = PeakFitUtils::is_high_res( spec );
    const auto resType = PeakFitUtils::coarse_resolution_from_peaks(inpeaks);
    const bool isHPGe = (resType == PeakFitUtils::CoarseResolutionType::High);

    answer = fit_peaks_in_roi_LM( inpeaks, data, isHPGe, fit_options );


    //now we need to go through and make sure the peaks we're adding are both
    //  significant, and improve the chi2/dof from before the fit.
    const double lx = origCont->lowerEnergy();
    const double ux = origCont->upperEnergy();
    const int lower_channel = static_cast<int>( data->find_gamma_channel( lx ) );
    const int upper_channel = static_cast<int>( data->find_gamma_channel( ux ) );
    const int nbin = (upper_channel > lower_channel) ? (upper_channel - lower_channel) : 1;
    const double prechi2Dof = chi2_for_region( inpeaks, data, lower_channel, upper_channel ) / nbin;
    const double postchi2Dof = chi2_for_region( answer, data, lower_channel, upper_channel ) / nbin;


    if( prechi2Dof < postchi2Dof )
    {
      answer.clear();

      if( true )
      {
        const double ncausalitysigma = 0.0;
        const double stat_threshold  = 0.0;
        const double hypothesis_threshold = 0.0;

        const bool isRefit = true;
        const vector<shared_ptr<const PeakDef>> refit_peaks
                     = fit_peaks_in_range_LM( lx, ux, ncausalitysigma, stat_threshold, hypothesis_threshold,
                                             inpeaks, data, isRefit, isHPGe );


        if( refit_peaks.size() == inpeaks.size() )
        {
          const double refit_chi2Dof = chi2_for_region( refit_peaks, data, lower_channel, upper_channel ) / nbin;
          cout << "refitPeaksThatShareROI_LM: refit_chi2Dof=" << refit_chi2Dof << ", prechi2Dof=" << prechi2Dof << ", postchi2Dof=" << postchi2Dof << endl;
          if( refit_chi2Dof <= prechi2Dof )
          {
            //cout << "Using re-fit peaks!" << endl;
            answer = refit_peaks;
          }
        }//if( output_peak.size() == inpeaks.size() )
      }

      if( answer.empty() )
        return answer;
    }//if( prechi2Dof >= 0.99*postchi2Dof )


    for( const shared_ptr<const PeakDef> &peak : answer )
    {
      const double mean = peak->mean();
      const double sigma = peak->sigma();
      const double data_counts = SpecUtils::gamma_integral( data, mean-0.5*sigma, mean+0.5*sigma );
      const double area = peak->gauss_integral( mean-0.5*sigma, mean+0.5*sigma );

      if( area < 5.0 || area < sqrt(data_counts) )
      {
        answer.clear();
        break;
      }//if( area < 5.0 || area < sqrt(data) )

      if( answer.size() == 1 )
      {
        static int ntimes = 0;
        if( ntimes++ < 3 )
          cerr << "refitPeaksThatShareROI_LM: Should perform some additional chi2 checks for the "
          "case answer.first.size() == 1" << endl;
      }//if( answer.first.size() == 1 )
    }//for( const PeakDefShrdPtr peak : answer.first )

  }catch( std::exception & )
  {
    answer.clear();
    cerr << "refitPeaksThatShareROI_LM: failed to find a better fit" << endl;
  }

  return answer;
}//refitPeaksThatShareROI_LM(...)

}//namespace PeakFitLM
