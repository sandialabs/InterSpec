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
#include "InterSpec/PeakFitUtils.h"
#include "SpecUtils/EnergyCalibration.h"
#include "InterSpec/DetectorPeakResponse.h"

#include "InterSpec/PeakFit_imp.hpp"
#include "InterSpec/PeakDists_imp.hpp"
#include "InterSpec/RelActCalcAuto_imp.hpp"

using namespace std;

namespace PeakFitLM
{

/** PeakFitDiffCostFunction is a not exactly apples-to-apples experiment to see how well Ceres does fitting peaks, relative to Minuit2
    This class is somewhat equivalent of the MultiPeakFitChi2Fcn class.
    Or maybe more directly analogous to LinearProblemSubSolveChi2Fcn

 Current shortcommings of this class
 - Should add choice of a free-for-all for peak means, or to have it as a function of energy, like now.

 Takes 1 parameter blocks, of size number_parameters() - e.g., we are putting all parameters together:
 Emits `number_residuals()` residuals

 TODO: currently will uses the entire first/last bin or ROI energy range, no matter how-little the ROI extends into those channels - should probably do some rounding or something.
 TODO: PeakFitChi2Fcn also punishes for peak being statistically insignificant
 */
struct PeakFitDiffCostFunction
{
  PeakFitDiffCostFunction(const std::shared_ptr<const SpecUtils::Measurement> data,
                          const std::vector< std::shared_ptr<const PeakDef> > &starting_peaks,
                          const double roi_lower_energy,
                          const double roi_upper_energy,
                          const PeakContinuum::OffsetType offset_type,
                          const double continuum_ref_energy,
                          const PeakDef::SkewType skew_type )
  : m_data( data ),
    m_lower_channel( data->find_gamma_channel( roi_lower_energy ) ),
    m_upper_channel( data->find_gamma_channel( roi_upper_energy ) ),
    m_roi_lower_energy( roi_lower_energy ),
    m_roi_upper_energy( roi_upper_energy ),
    m_num_peaks( starting_peaks.size() ),
    // Make sure m_starting_peaks is sorted - use lambda trick so all member variables are const
    m_starting_peaks( ([starting_peaks]() -> vector<shared_ptr<const PeakDef>>{
      vector<shared_ptr<const PeakDef>> sorted_peaks = starting_peaks;
      std::sort( begin(sorted_peaks), end(sorted_peaks), &PeakDef::lessThanByMeanShrdPtr );
      return sorted_peaks;
    } )() ),
    // The number of parameters depend on how many sigmas are being fit; again, use lambda trick
    m_num_fit_sigmas( ([&starting_peaks]() -> size_t{
      size_t n_fit_sigma = 0;
      for( const auto &p : starting_peaks )
        n_fit_sigma += p->fitFor(PeakDef::Sigma);
      return n_fit_sigma;
    })() ),
    m_use_lls_for_cont( ([&starting_peaks,offset_type]() -> bool {
      auto continuum = starting_peaks.empty() ? nullptr : starting_peaks.front()->continuum();
      assert( continuum && (continuum->type() == offset_type) );
      const vector<bool> par_fit_for = continuum ? continuum->fitForParameter() : vector<bool>{};
      for( const bool do_fit : par_fit_for )
      {
        if( !do_fit )
          return false;
      }
      return true;
    })() ),
    m_offset_type( offset_type ),
    m_ref_energy( continuum_ref_energy ),
    m_skew_type( skew_type ),
    m_ncalls( 0 )
  {
    assert( m_data );
    assert( m_data->gamma_counts() && !m_data->gamma_counts()->empty() );
    assert( m_upper_channel >= m_lower_channel );
    assert( m_upper_channel < m_data->num_gamma_channels() );
    assert( m_roi_lower_energy < m_roi_upper_energy );
    
    assert( m_num_peaks > 0 );
    assert( m_num_peaks < (m_upper_channel - m_lower_channel) );
    
    if( !m_data || !m_data->gamma_counts() || m_data->gamma_counts()->empty() )
      throw runtime_error( "PeakFitDiffCostFunction: !m_data or !gamma_counts()" );
    
    if( m_roi_lower_energy >= m_roi_upper_energy )
      throw runtime_error( "PeakFitDiffCostFunction: m_roi_lower_energy >= m_roi_upper_energy" );
    
    if( m_upper_channel < m_lower_channel )
      throw runtime_error( "PeakFitDiffCostFunction: m_upper_channel < m_lower_channel" );
      
    const std::vector<float> &counts = *m_data->gamma_counts();
    if( m_upper_channel >= counts.size() )
      throw runtime_error( "PeakFitDiffCostFunction: m_upper_channel >= counts.size()" );

    const shared_ptr<const vector<float>> &energies = m_data->channel_energies();
    if( !energies || (m_upper_channel > energies->size()) )
      throw runtime_error( "PeakFitDiffCostFunction: m_upper_channel >= energies.size()" );

    if( !m_num_peaks )
      throw runtime_error( "PeakFitDiffCostFunction: no peaks to fit" );
    
    for( size_t i = 1; i < starting_peaks.size(); ++i )
    {
      if( starting_peaks[i]->continuum() != starting_peaks[0]->continuum() )
        throw runtime_error( "PeakFitDiffCostFunction: input peaks must all share a continuum" );
    }
    
    auto continuum = m_starting_peaks[0]->continuum();
    // TODO: get rid of m_roi_lower_energy, m_roi_upper_energy, m_offset_type, and m_ref_energy - since they are all duplicated
    assert( continuum->lowerEnergy() == m_roi_lower_energy );
    assert( continuum->upperEnergy() == m_roi_upper_energy );
    assert( continuum->type() == m_offset_type );
    assert( continuum->referenceEnergy() == m_ref_energy );
    
    
    
    bool valid_offset_type = false;
    switch( offset_type )
    {
      case PeakContinuum::NoOffset:     case PeakContinuum::Constant:
      case PeakContinuum::Linear:       case PeakContinuum::Quadratic:
      case PeakContinuum::Cubic:        case PeakContinuum::FlatStep:
      case PeakContinuum::LinearStep:   case PeakContinuum::BiLinearStep:
      case PeakContinuum::External:
        valid_offset_type = true;
        break;
    }//switch( offset_type )
    
    assert( valid_offset_type );
    
    if( !valid_offset_type )
      throw runtime_error( "PeakFitDiffCostFunction: !valid_offset_type" );
    
    // It wouldnt be that hard to support external continuum, but we'll leave that for later, at
    //  the moment.
    assert( m_offset_type != PeakContinuum::OffsetType::External );
    if( m_offset_type == PeakContinuum::OffsetType::External )
      throw runtime_error( "PeakFitDiffCostFunction: PeakContinuum::OffsetType::External not supported" );
    
    if( dof() <= 0.0 )
      throw runtime_error( "PeakFitDiffCostFunction: not enough data channels to fit all the parameters" );
  }//PeakFitDiffCostFunction constructor
  
  double dof() const
  {
    const double num_channels = 1.0 + m_upper_channel - m_lower_channel;
    const double num_skew = PeakDef::num_skew_parameters( m_skew_type );
    const double num_continuum_pars = PeakContinuum::num_parameters( m_offset_type );
    const double num_sigmas_fit = number_sigma_parameters();
    
    // fixed sigmas are handled by just not having them included in the parameters, so we wont
    //  count fixed sigmas in this next loop.
    size_t num_fixed_pars = 0;
    for( const auto &p : m_starting_peaks )
    {
      num_fixed_pars += !p->fitFor(PeakDef::GaussAmplitude);
      num_fixed_pars += !p->fitFor(PeakDef::Mean);
    }
    
    // Get the number of fixed continuum parameters.
    for( const bool do_fit_par : m_starting_peaks[0]->continuum()->fitForParameter() )
      num_fixed_pars += !do_fit_par;
    
    // TODO: MultiPeakFitChi2Fcn::dof() adds 1 to the degrees of freedom - not quite sure why that is, at the moment, or if we should do this here.
    return num_channels - 2.0*m_num_peaks + num_fixed_pars - num_sigmas_fit - num_continuum_pars;
  }//double dof() const
  
  size_t number_parameters() const
  {
    const size_t num_skew = PeakDef::num_skew_parameters( m_skew_type );
    const size_t num_continuum_pars = PeakContinuum::num_parameters( m_offset_type );

    return num_continuum_pars + num_skew + number_sigma_parameters() + m_num_peaks;
  }
  
  size_t number_sigma_parameters() const
  {
    return std::min( m_num_fit_sigmas, size_t(2) );
  }
  
  size_t number_residuals() const
  {
    return m_num_peaks + m_upper_channel - m_lower_channel;
  }//size_t number_residuals() const


  template<typename T>
  struct PeakDefImpWithCont : public RelActCalcAuto::PeakDefImp<T>
  {
    std::shared_ptr<RelActCalcAuto::PeakContinuumImp<T>> m_continuum;

    std::shared_ptr<RelActCalcAuto::PeakContinuumImp<T>> getContinuum(){ return m_continuum; }


    void setContinuum( std::shared_ptr<RelActCalcAuto::PeakContinuumImp<T>> continuum )
    {
      m_continuum = continuum;
    }

    void setMeanUncert( const T mean_uncert )
    {
      //
    }

    void setAmplitudeUncert( const T amp_uncert )
    {

    }

    void setSigmaUncert( const T sigma_uncert )
    {

    }


    static bool lessThanByMean( const PeakDefImpWithCont<T> &lhs, const PeakDefImpWithCont<T> &rhs )
    {
      return (lhs.m_mean < rhs.m_mean);
    }
  };//struct PeakDefImpWithCont


  //
  template<typename T>
  static std::shared_ptr<T> create_continuum( const std::shared_ptr<T> &other_cont [[maybe_unused]]) {
      return std::make_shared<T>();
  }

  template<typename PeakType,typename T>
  vector<PeakType> parametersToPeaks( const T * const params, const T * const uncertainties, T *residuals ) const
  {
    const size_t num_continuum_pars = PeakContinuum::num_parameters( m_offset_type );
    const size_t last_data_residual = m_upper_channel - m_lower_channel;

    auto starting_continuum = m_starting_peaks[0]->continuum();
    assert( starting_continuum->lowerEnergy() == m_roi_lower_energy );
    assert( starting_continuum->upperEnergy() == m_roi_upper_energy );
    assert( starting_continuum->type() == m_offset_type );
    assert( starting_continuum->referenceEnergy() == m_ref_energy );


    std::unique_ptr<vector<T>> local_residuals;
    if( !residuals )
    {
      local_residuals.reset( new vector<T>( number_residuals(), T(0.0) ) );
      residuals = local_residuals->data();
    }else
    {
      // Not sure if we need to do this or not
      const size_t nresids = number_residuals();
      for( size_t i = 0; i < nresids; ++i )
        residuals[i] = T(0.0);
    }

    const size_t num_sigmas_fit = number_sigma_parameters();
    const size_t num_skew = PeakDef::num_skew_parameters( m_skew_type );

    vector<PeakType> peaks, fixed_amp_peaks;

    for( size_t i = 0; i < m_num_peaks; ++i )
    {
      const size_t start_par = num_continuum_pars + num_skew + num_sigmas_fit + i;
      assert( start_par < number_parameters() );

      const bool fit_amp = m_starting_peaks[i]->fitFor(PeakDef::GaussAmplitude);

      const T mean = params[start_par + 0];
      const T amp = fit_amp ? T(1.0) : T(m_starting_peaks[i]->amplitude());

      T sigma, sigma_uncert;

      if( !m_starting_peaks[i]->fitFor(PeakDef::Sigma) || (num_sigmas_fit == 0) )
      {
        sigma = T(m_starting_peaks[i]->sigma());
        sigma_uncert = T(m_starting_peaks[i]->sigmaUncert());
      }else
      {
        sigma = params[num_continuum_pars + num_skew];
        sigma_uncert = uncertainties ? uncertainties[num_continuum_pars + num_skew] : T(0.0);

        if( num_sigmas_fit > 1 )
        {
          assert( (num_continuum_pars + num_skew + 1) < number_parameters() );

          const T frac = (mean - m_roi_lower_energy) / (m_roi_upper_energy - m_roi_lower_energy);
          sigma = params[num_continuum_pars + num_skew] + (frac * params[num_continuum_pars + num_skew + 1]);

          if( uncertainties )
          {
            //TODO: I'm not entirely conviced uncert is correct (note: assuming 100% correlation)
            //TODO: use the covariance between these two parameters to better assign sigma uncert
            sigma_uncert += frac * uncertainties[num_continuum_pars + num_skew + 1];
          }
        }//if( num_sigmas_fit > 1 )
      }//if( num_sigmas_fit == 0 ) / else

      PeakType peak;
      peak.setMean( mean );
      peak.setSigma( sigma );
      peak.setAmplitude( amp );

      peak.setSkewType( m_skew_type );
      for( size_t skew_index = 0; skew_index < num_skew; ++skew_index )
      {
        const PeakDef::CoefficientType coeff_num = static_cast<PeakDef::CoefficientType>( PeakDef::CoefficientType::SkewPar0 + static_cast<int>(skew_index) );
        const T skew_coef = params[num_continuum_pars + skew_index];
        peak.set_coefficient( skew_coef, coeff_num );
      }

      if( uncertainties )
      {
        T mean_uncert = uncertainties[start_par + 0];
        T amp_uncert = uncertainties[start_par + 1];

        if( !m_starting_peaks[i]->fitFor(PeakDef::Mean) )
          mean_uncert = T(m_starting_peaks[i]->meanUncert());

        if( !m_starting_peaks[i]->fitFor(PeakDef::GaussAmplitude) )
          amp_uncert = T(m_starting_peaks[i]->amplitudeUncert());

        if( mean_uncert > 0.0 )
          peak.setMeanUncert( mean_uncert );

        if( amp_uncert > 0.0 )
          peak.setAmplitudeUncert( amp_uncert );

        if( sigma_uncert > 0.0 )
          peak.setSigmaUncert( sigma_uncert );
      }//if( uncertainties )

      if( fit_amp )
        peaks.push_back( peak );
      else
        fixed_amp_peaks.push_back( peak );
    }//for( size_t i = 0; i < m_num_peaks; ++i )

    // Sort the peaks to make calculating the punishment for peaks being to close easier
    std::sort( begin(peaks), end(peaks), &PeakType::lessThanByMean );


    const size_t nchannel = (m_upper_channel - m_lower_channel) + 1;
    const shared_ptr<const vector<float>> &energies_ptr = m_data->channel_energies();
    assert( energies_ptr && (energies_ptr->size() > (m_lower_channel + nchannel)) );

    assert( !!m_data->gamma_counts() );
    const vector<float> &counts_vec = *m_data->gamma_counts();
    assert( counts_vec.size() > (m_lower_channel + nchannel) );

    vector<T> peak_counts( nchannel, T(0.0) );
    const float * const energies = energies_ptr->data() + m_lower_channel;
    const float * const channel_counts = counts_vec.data() + m_lower_channel;

    const int num_polynomial_terms = static_cast<int>( PeakContinuum::num_parameters( m_offset_type ) );
    const bool is_step_continuum = PeakContinuum::is_step_continuum( m_offset_type );

    auto continuum = create_continuum( peaks.front().getContinuum() );
    continuum->setRange( T(m_roi_lower_energy), T(m_roi_upper_energy) );
    continuum->setType( m_offset_type );

    if( !m_use_lls_for_cont )
    {
      continuum->setParameters( T(m_ref_energy), params, nullptr );

      for( size_t peak_index = 0; peak_index < peaks.size(); ++peak_index )
        peaks[peak_index].gauss_integral( energies, &(peak_counts[0]), nchannel );

      if constexpr ( !std::is_same_v<T, double> )
      {
        PeakDists::offset_integral( *continuum, energies, &(peak_counts[0]), nchannel, m_data );
      }else
      {
        continuum->offset_integral( energies, &(peak_counts[0]), nchannel, m_data );
      }

    }else
    {
      const bool multithread = true;
      T continuum_coeffs[5] = { T(0.0) };
      assert( num_polynomial_terms <= 5 );

      switch( m_offset_type )
      {
        case PeakContinuum::Constant:
        case PeakContinuum::Linear:       case PeakContinuum::Quadratic:
        case PeakContinuum::Cubic:        case PeakContinuum::FlatStep:
        case PeakContinuum::LinearStep:   case PeakContinuum::BiLinearStep:
        {
          vector<T> means, sigmas;
          for( const auto &peak : peaks )
          {
            means.push_back( peak.mean() );
            sigmas.push_back( peak.sigma() );
          }

          const vector<T> skew_pars( params + num_continuum_pars, params + num_continuum_pars + num_skew );

          vector<T> amplitudes, continuum_coeffs, amp_uncerts, cont_uncerts;

          PeakFit::fit_amp_and_offset_imp(energies, channel_counts, nchannel, num_polynomial_terms, is_step_continuum,
                                          T(m_ref_energy), means, sigmas, fixed_amp_peaks, m_skew_type, skew_pars.data(),
                                          amplitudes, continuum_coeffs, amp_uncerts, cont_uncerts,
                                          &(peak_counts[0]) );

          //PeakFit::fit_continuum( energies, channel_counts, nchannel, num_polynomial_terms,
          //                       is_step_continuum, T(m_ref_energy), peaks, multithread,
          //                       continuum_coeffs, &(peak_counts[0]) );

          assert( peaks.size() == amplitudes.size() );
          assert( amp_uncerts.empty() || (amp_uncerts.size() == peaks.size()));
          for( size_t peak_index = 0; peak_index < peaks.size(); ++peak_index )
          {
            peaks[peak_index].setAmplitude( amplitudes[peak_index] );
            if( !amp_uncerts.empty() )
              peaks[peak_index].setAmplitudeUncert( amp_uncerts[peak_index] );
          }

          continuum->setParameters( T(m_ref_energy), continuum_coeffs.data(), cont_uncerts.data() );
          break;
        }


        case PeakContinuum::NoOffset:
        case PeakContinuum::External:
          for( size_t peak_index = 0; peak_index < peaks.size(); ++peak_index )
            peaks[peak_index].gauss_integral( energies, &(peak_counts[0]) , nchannel );

          if( m_offset_type == PeakContinuum::OffsetType::External )
          {
            //continuum->setExternalContinuum( const std::shared_ptr<const SpecUtils::Measurement> &data );
            assert( 0 );
            throw runtime_error( "External continuum not implemented for PeakFitDiffCostFunction yet" );
          }
          break;
      }//switch( offset_type )
    }//if( m_use_lls_for_cont )


    // Put fixed amplitude peaks into `peaks`.
    if( !fixed_amp_peaks.empty() )
    {
      for( auto &peak : fixed_amp_peaks )
        peaks.push_back( peak );
      std::sort( begin(peaks), end(peaks), &PeakType::lessThanByMean );
    }


    for( auto &peak : peaks )
      peak.setContinuum( continuum );


    T chi2( 0.0 );
    for( size_t residual_num = 0; residual_num < nchannel; ++residual_num )
    {
      const double ndata = channel_counts[residual_num];

      assert( residual_num < number_residuals() );

      if( ndata > 0.000001f )  // The cutoff is arbitrary, but needs to be greater than zero
        residuals[residual_num] = (ndata - peak_counts[residual_num]) / sqrt(ndata);
      else
        residuals[residual_num] = peak_counts[residual_num]; // TODO: This is ad-hoc, and just following PeakFitChi2Fcn implementation for the moment

      chi2 += residuals[residual_num] * residuals[residual_num];
    }//for( size_t channel = m_lower_channel; channel <= m_upper_channel; ++channel )

    const T chi2_dof = chi2 / dof();
    for( size_t i = 0; i < m_num_peaks; ++i )
      peaks[i].set_coefficient( chi2_dof, PeakDef::Chi2DOF );

    // Punish if the peaks are too close.
    for( size_t i = 1; i < peaks.size(); ++i )
    {
      T sigma;
      if constexpr ( !std::is_same_v<T, double> )
      {
        const T sigma_i = peaks[i].sigma();
        const T sigma_j = peaks[i-1].sigma();
        sigma = 0.5*(sigma_j + sigma_i);
      }else
      {
        const double sigma_i = peaks[i].gausPeak() ? peaks[i].sigma() : 0.25*peaks[i].roiWidth();
        const double sigma_j = peaks[i-1].gausPeak() ? peaks[i-1].sigma() : 0.25*peaks[i-1].roiWidth();
        sigma = 0.5*(sigma_j + sigma_i);
      }


      const T dist = abs(peaks[i-1].mean() - peaks[i].mean());
      T reldist = dist / sigma;
      if( (reldist < 0.01) || isinf(reldist) || isnan(reldist) )
      {
        if constexpr ( !std::is_same_v<T, double> )
        {
          reldist.a = 0.01;
          reldist.v = peaks[i-1].mean().v - peaks[i].mean().v;
        }else
        {
          reldist = 0.01;
        }
      }

      const double punishment_factor = 2 * (m_upper_channel - m_lower_channel) / m_num_peaks;

      assert( (last_data_residual + i) < number_residuals() );

      if( reldist < 1.0 )
      {
        residuals[last_data_residual + i] = (punishment_factor / reldist);
      }else if( reldist < 1.5 )
      {
        // At reldist==1, this will be equal to the above case,
        //  and by reldist==1.5: exp( -pow(5*0.5,3.0) )== 1.64E-7
        residuals[last_data_residual + i] = punishment_factor * exp( -pow(5.0*(reldist - 1.0),3.0) );
      }else
      {
        residuals[last_data_residual + i] = T(0.0);
      }
    }//for( size_t i = 1; i < peaks.size(); ++i )

    // TODO: PeakFitChi2Fcn also punishes for peak being statistically insignificant


    if constexpr ( std::is_same_v<PeakType, PeakDef> )
    {
      // Now inherit all the assigned nuclides, colors, whatever from the starting peaks.
      //  However, ordering of peaks could swap, etc, so we need to match the original
      //  peaks up to the closest current peak.
      assert( m_starting_peaks.size() == peaks.size() );
      map<size_t,size_t> old_to_new_indexs;

      // If the mean of a peak is fixed, then we "know" it wont swap positions
      for( size_t index = 0; index < m_starting_peaks.size(); ++index )
      {
        if( !m_starting_peaks[index]->fitFor(PeakDef::Mean) )
          old_to_new_indexs[index] = index;
      }

      for( size_t new_index = 0; new_index < peaks.size(); ++new_index )
      {
        const double new_mean = peaks[new_index].mean();

        // Make sure we havent already assigned this new peak to an old peak (because mean was fixed)
        bool already_assigned = false;
        for( const auto &t : old_to_new_indexs )
          already_assigned = (already_assigned || (t.second == new_index));

        if( already_assigned )
          continue;

        int closest_orig = -1;
        for( size_t orig_index = 0; orig_index < m_starting_peaks.size(); ++orig_index )
        {
          // Check if we've already assigned this index
          if( old_to_new_indexs.count(orig_index) )
            continue;

          if( closest_orig < 0 )
          {
            closest_orig = static_cast<int>( orig_index );
          }else
          {
            const double prev_mean = m_starting_peaks[closest_orig]->mean();
            const double orig_mean = m_starting_peaks[orig_index]->mean();
            if( fabs(new_mean - orig_mean) < fabs(new_mean - prev_mean) )
              closest_orig = static_cast<int>( orig_index );
          }
        }//for( loop over original peaks )

        assert( closest_orig >= 0 );

        old_to_new_indexs[static_cast<size_t>(closest_orig)] = new_index;
      }//for( loop over new peaks )

      assert( old_to_new_indexs.size() == peaks.size() );

      for( const auto &old_to_new : old_to_new_indexs )
      {
        assert( old_to_new.first < peaks.size() );
        assert( old_to_new.second < peaks.size() );
        const bool inheritNonFitForValues = true; //I dont think this really matters
        peaks[old_to_new.second].inheritUserSelectedOptions( *m_starting_peaks[old_to_new.first],  inheritNonFitForValues );
      }//for( const auto &old_to_new : old_to_new_indexs )
    }//

    return peaks;
  }


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

  void setup_skew_parameters( double *pars,
                             vector<int> &constant_parameters,
                             vector<std::optional<double>> &lower_bounds,
                             vector<std::optional<double>> &upper_bounds,
                          const vector<shared_ptr<const PeakDef> > &inpeaks ) const
  {
    const size_t num_skew_pars = PeakDef::num_skew_parameters( m_skew_type );

    if( num_skew_pars == 0 )
      return;

    vector<bool> fit_parameter( num_skew_pars, true );
    vector<double> starting_value( num_skew_pars, 0.0 );
    vector<double> lower_values( num_skew_pars, 0.0 ), upper_values( num_skew_pars, 0.0 );

    for( size_t i = 0; i < num_skew_pars; ++i )
    {
      const auto ct = PeakDef::CoefficientType( PeakDef::CoefficientType::SkewPar0 + i);

      double lower, upper, start, dx;
      const bool use = PeakDef::skew_parameter_range( m_skew_type, ct, lower, upper, start, dx );
      assert( use );
      if( !use )
        throw logic_error( "Inconsistent skew par val" );

      starting_value[i] = start;
      lower_values[i] = lower;
      upper_values[i] = upper;
    }//for( size_t i = 0; i < num_skew_pars; ++i )


    for( const auto &p : inpeaks )
    {
      if( p->skewType() != m_skew_type )
        continue;

      for( size_t i = 0; i < num_skew_pars; ++i )
      {
        const auto ct = PeakDef::CoefficientType( PeakDef::CoefficientType::SkewPar0 + i);
        fit_parameter[i] = p->fitFor( ct );
        double val = p->coefficient( ct );

        if( IsInf(val) || IsNan(val) || (val < lower_values[i]) || (val > upper_values[i]) )
          val = starting_value[i];

        // When we are dragging a ROI edge, sometimes the power-law can get stuck up around 100,
        //  so lets try to avoid this
        switch( m_skew_type )
        {
          case PeakDef::NoSkew:   case PeakDef::Bortel:
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
            }//switch( ct )
            break;
          }//A Crystal Ball distribution
        }//switch( m_skew_type )

        starting_value[i] = val;
      }//for( size_t i = 0; i < num_skew_pars; ++i )
    }//for( const auto &p : inpeaks )

    assert( num_skew_pars == fit_parameter.size() );
    assert( num_skew_pars == starting_value.size() );

    const size_t num_continuum_pars = PeakContinuum::num_parameters( m_offset_type );

    for( size_t skew_index = 0; skew_index < num_skew_pars; ++skew_index )
    {
      const size_t par_index = num_continuum_pars + skew_index;
      pars[par_index] = starting_value[skew_index];

      if( fit_parameter[skew_index] )
      {
        lower_bounds[par_index] = std::max( lower_values[skew_index], 0.5*starting_value[skew_index] );
        upper_bounds[par_index] = std::min( upper_values[skew_index], 1.5*fabs(starting_value[skew_index]) );
      }else
      {
        constant_parameters.push_back( static_cast<int>(par_index) );
      }
    }//for( size_t skew_index = 0; skew_index < num_skew_pars; ++skew_index )
  }//setup_skew_parameters(...)

public:
  const std::shared_ptr<const SpecUtils::Measurement> m_data;
  const size_t m_lower_channel;
  const size_t m_upper_channel;
  const double m_roi_lower_energy;
  const double m_roi_upper_energy;
  
  const size_t m_num_peaks;
  const std::vector< std::shared_ptr<const PeakDef> > m_starting_peaks;
  const size_t m_num_fit_sigmas;
  const bool m_use_lls_for_cont;
  const PeakContinuum::OffsetType m_offset_type;
  const double m_ref_energy;
  const PeakDef::SkewType m_skew_type;

  mutable std::atomic<unsigned int> m_ncalls;
};//struct PeakFitDiffCostFunction



void fit_peak_for_user_click_LM( PeakShrdVec &results,
                             double &chi2Dof,
                             const std::shared_ptr<const SpecUtils::Measurement> &dataH,
                             PeakShrdVec coFitPeaks,
                             const double mean0, const double sigma0,
                             const double area0,
                             const float roiLowerEnergy,
                             const float roiUpperEnergy,
                             const bool isHPGe )
{
  /** For this first go, we will have Ceres fit for things.

   This is only tested to the "seems to work, for simple cases" level.
   
   We will:
     - Constrains the widths of the peaks to have a linearly related width
   
   In the future:
     - Add punishment for stat insignificant peaks
     - Need to deal with errors during solving, and during Covariance computation
     - Deal with setting number of threads reasonably.
     - Try out using a ceres::LossFunction
     - ...
   
   */
  
  typedef std::shared_ptr<const PeakDef> PeakDefShrdPtr;
  
  // Lets make sure all coFitPeaks share a continuum.
  for( size_t i = 1; i < coFitPeaks.size(); ++i )
  {
    if( coFitPeaks[i]->continuum() != coFitPeaks[0]->continuum() )
      throw runtime_error( "fit_peak_for_user_click_LM: input peaks all must share a single continuum" );
  }//for( size_t i = 1; i < coFitPeaks.size(); ++i )
  
  chi2Dof = DBL_MAX;
  results.clear();
  
  try
  {
    // TODO: check that repeaded calls to this function wont cause adding/removing channels due to find_gamma_channel(...) rounding type things - e.g., do we need to subtract off a tiny bit from roiUpperEnergy
    const size_t lower_channel = dataH->find_gamma_channel(roiLowerEnergy);
    const size_t upper_channel = dataH->find_gamma_channel(roiUpperEnergy);
    const double roi_lower = dataH->gamma_channel_lower( lower_channel );
    const double roi_upper = dataH->gamma_channel_upper( upper_channel );
    
    assert( roi_lower <= roiLowerEnergy );
    assert( roi_upper >= roiUpperEnergy );
    
    const size_t nchannels = dataH->num_gamma_channels();
    const size_t midbin = dataH->find_gamma_channel( mean0 );
    const float binwidth = dataH->gamma_channel_width( midbin );
    const size_t num_fit_peaks = coFitPeaks.size() + 1;
    
    double reference_energy = roiLowerEnergy;
    
    //The below should probably go off the number of bins in the ROI
    PeakContinuum::OffsetType offset_type;
    if( isHPGe )
      offset_type = (num_fit_peaks < 3) ? PeakContinuum::Linear : PeakContinuum::Quadratic;
    else
      offset_type = (num_fit_peaks < 2) ? PeakContinuum::Linear : PeakContinuum::Quadratic;
    
    //for( size_t i = 0; i < coFitPeaks.size(); ++i )
    if( coFitPeaks.size() )
    {
      //const auto continuum = coFitPeaks[i]->continuum();
      const auto continuum = coFitPeaks[0]->continuum();
      offset_type = std::max( offset_type, continuum->type() );
    }
    
    const size_t num_continuum_pars = PeakContinuum::num_parameters( offset_type );
    
    bool any_continuum_par_fixed = false;
    // TODO: we can get rid of continuum_parameter_fit... as this will be contained in `initial_continuum` as well.
    vector<bool> continuum_parameter_fit( num_continuum_pars, true );
    //for( size_t i = 0; i < coFitPeaks.size(); ++i )
    if( coFitPeaks.size() )
    {
      //const auto continuum = coFitPeaks[i]->continuum();
      const auto continuum = coFitPeaks[0]->continuum();
      const std::vector<bool> parFitFors = continuum->fitForParameter();
      assert( parFitFors.size() <= num_continuum_pars );
      
      for( size_t fit_for_index = 0; fit_for_index < parFitFors.size(); ++fit_for_index )
      {
        if( !parFitFors[fit_for_index] )
        {
          reference_energy = continuum->referenceEnergy();
          any_continuum_par_fixed = true;
          continuum_parameter_fit[fit_for_index] = false;
        }
      }//
    }//for( size_t i = 0; i < coFitPeaks.size(); ++i )

    const PeakDef::SkewType skew_type = PeakFitDiffCostFunction::skew_type_from_prev_peaks( coFitPeaks );
    const size_t num_skew = PeakDef::num_skew_parameters( skew_type );

    std::shared_ptr<PeakDef> candidatepeak = std::make_shared<PeakDef>(mean0, sigma0, area0);
    
    // Make sure all initial peaks share a continuum
    //  TODO: this common continuum shares a good amount of information we also pass to PeakFitDiffCostFunction ... should consolidate
    shared_ptr<PeakContinuum> initial_continuum;
    if( coFitPeaks.size() )
    {
      initial_continuum = make_shared<PeakContinuum>( *coFitPeaks[0]->continuum() );
      
      candidatepeak->setContinuum( initial_continuum );
      for( size_t i = 0; i < coFitPeaks.size(); ++i )
      {
        auto newpeak = make_shared<PeakDef>( *coFitPeaks[i] );
        newpeak->setContinuum( initial_continuum );
        coFitPeaks[i] = newpeak;
      }
    }else
    {
      initial_continuum = candidatepeak->continuum();
    }//if( coFitPeaks.size() )
    
    assert( initial_continuum );
    
    initial_continuum->setRange( roiLowerEnergy, roiUpperEnergy );
    if( !any_continuum_par_fixed )
      initial_continuum->calc_linear_continuum_eqn( dataH, reference_energy, roiLowerEnergy, roiUpperEnergy, 2, 2 );
    initial_continuum->setType( offset_type );


    coFitPeaks.push_back( candidatepeak );
    std::sort( coFitPeaks.begin(), coFitPeaks.end(), &PeakDef::lessThanByMeanShrdPtr );

    
    assert( roiLowerEnergy < roiUpperEnergy );
    if( roiLowerEnergy >= roiUpperEnergy )
      throw runtime_error( "Invalid energy range (" + std::to_string(roiLowerEnergy) + ", " + std::to_string(roiUpperEnergy) + ")" );


    assert( num_continuum_pars == initial_continuum->parameters().size() );


    const double range = (roi_upper - roi_lower);
      
    double minsigma = binwidth;
    double maxsigma = 0.5*range;
    
    auto cost_functor = new PeakFitDiffCostFunction( dataH, coFitPeaks, roiLowerEnergy, roiUpperEnergy, offset_type, reference_energy, skew_type );

    auto cost_function = new ceres::DynamicAutoDiffCostFunction<PeakFitDiffCostFunction,12>( cost_functor );

    const size_t num_fit_pars = cost_functor->number_parameters();
    const size_t num_sigmas_fit = cost_functor->number_sigma_parameters();

    cost_function->AddParameterBlock( static_cast<int>(num_fit_pars) );
    cost_function->SetNumResiduals( static_cast<int>(cost_functor->number_residuals()) );
    
    vector<double> parameters( num_fit_pars, 0.0 );
    double *pars = &parameters[0];


    vector<int> constant_parameters;
    vector<std::optional<double>> lower_bounds( num_fit_pars ), upper_bounds( num_fit_pars );

    if( cost_functor->m_use_lls_for_cont )
    {
      for( size_t i = 0; i < num_continuum_pars; ++i )
      {
        parameters[i] = -1.0;
        constant_parameters.push_back( static_cast<int>(i) );
      }
    }else
    {
      assert( initial_continuum );
      const vector<double> &cont_pars = initial_continuum->parameters();
      const vector<bool> par_fit_for = initial_continuum->fitForParameter();
      assert( cont_pars.size() == num_continuum_pars );
      assert( par_fit_for.size() == num_continuum_pars );

      for( size_t i = 0; i < num_continuum_pars; ++i )
      {
        parameters[i] = cont_pars[i];
        if( !continuum_parameter_fit[i] )
          constant_parameters.push_back( static_cast<int>(i) );
      }
    }//if( cost_function->m_use_lls_for_cont )


    cost_functor->setup_skew_parameters( pars, constant_parameters, lower_bounds, upper_bounds, coFitPeaks );


    for( size_t i = 0; i < coFitPeaks.size(); ++i )
    {
      const double mean = coFitPeaks[i]->mean();
      const double sigma = coFitPeaks[i]->sigma();
      const double amp = coFitPeaks[i]->amplitude();
      
      const size_t start_par = num_continuum_pars + num_skew + num_sigmas_fit + i;
      parameters[start_par] = mean;
      
      if( !coFitPeaks[i]->fitFor(PeakDef::Mean) )
      {
        constant_parameters.push_back( static_cast<int>(start_par + 0) );
      }else
      {
        if( coFitPeaks[i] == candidatepeak )
        {
          lower_bounds[start_par + 0] = mean - sigma;
          upper_bounds[start_par + 0] = mean + sigma;
        }else
        {
          lower_bounds[start_par + 0] = mean - 0.25*sigma;
          upper_bounds[start_par + 0] = mean + 0.25*sigma;
        }
      }//


      if( !coFitPeaks[i]->fitFor(PeakDef::Sigma) )
      {
        minsigma = std::min( minsigma, coFitPeaks[i]->sigma() );
        maxsigma = std::max( maxsigma, coFitPeaks[i]->sigma() );
      }else
      {
        float lowersigma, uppersigma;
        expected_peak_width_limits( mean, isHPGe, dataH, lowersigma, uppersigma );
        if( !i )
          minsigma = lowersigma;
        if( i == (coFitPeaks.size()-1) )
          maxsigma = 1.33*uppersigma;
      }
    }//for( size_t i = 0; i < coFitPeaks.size(); ++i )
    
    if( num_sigmas_fit > 0 )
    {
      lower_bounds[num_continuum_pars + num_skew] = 0.5*minsigma;
      upper_bounds[num_continuum_pars + num_skew] = 1.5*maxsigma;

      if( num_sigmas_fit > 1 )
      {
        lower_bounds[num_continuum_pars + num_skew + 1] = -0.15;
        upper_bounds[num_continuum_pars + num_skew + 1] = +0.15;
      }
    }//if( num_sigmas_fit > 0 )


    ceres::Problem problem;

    // TODO: investigate using a LossFunction; from a brief investigation it maybe really affects the amplitude uncertainty if not chosen carefully
    ceres::LossFunction *lossfcn = nullptr;
    //ceres::LossFunction *lossfcn = new ceres::CauchyLoss(0.5) or new HuberLoss(...), or ...

    problem.AddResidualBlock( cost_function, lossfcn, pars );


    if( !constant_parameters.empty() )
    {
      ceres::Manifold *subset_manifold = new ceres::SubsetManifold( static_cast<int>(num_fit_pars), constant_parameters );
      problem.SetManifold( pars, subset_manifold );
    }

    for( size_t i = 0; i < num_fit_pars; ++i )
    {
      if( lower_bounds[i].has_value() )
        problem.SetParameterLowerBound(pars, static_cast<int>(i), *lower_bounds[i] );
      if( upper_bounds[i].has_value() )
        problem.SetParameterUpperBound(pars, static_cast<int>(i), *upper_bounds[i] );
    }//for( size_t i = 0; i < num_fit_pars; ++i )


    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_type = ceres::TRUST_REGION; //ceres::LINE_SEARCH
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT; //ceres::DOGLEG
    options.use_nonmonotonic_steps = true;
    options.max_consecutive_nonmonotonic_steps = 5;
    options.initial_trust_region_radius = 1e4;
    options.max_trust_region_radius = 1e16;

    // Minimizer terminates when the trust region radius becomes smaller than this value.
    options.min_trust_region_radius = 1e-32;
    // Lower bound for the relative decrease before a step is accepted.
    options.min_relative_decrease = 1e-3;
    options.max_num_iterations = 50000;
#ifndef NDEBUG
    options.max_solver_time_in_seconds = 300.0;
#else
    options.max_solver_time_in_seconds = 120.0;
#endif
    options.function_tolerance = 1e-7; //default 1e-9;
    options.parameter_tolerance = 1e-11; //Default value is 1e-8.  Using 1e-11
    // TODO: there are a ton of ceres::Solver::Options that might be useful for us to set
    options.minimizer_progress_to_stdout = true;
    options.num_threads = 4;


    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    //std::cout << summary.BriefReport() << "\n";
    std::cout << summary.FullReport() << "\n";
    cout << "Took " << cost_functor->m_ncalls.load() << " calls to solve." << endl;
    cost_functor->m_ncalls = 0;
    
    switch( summary.termination_type )
    {
      case ceres::CONVERGENCE:
      case ceres::USER_SUCCESS:
        break;
        
      case ceres::NO_CONVERGENCE:
      case ceres::FAILURE:
      case ceres::USER_FAILURE:
        throw runtime_error( "The L-M ceres::Solver solving failed." );
        break;
    }//switch( summary.termination_type )
    
    
    //ceres::Solve(options, &problem, nullptr );
    
    
    ceres::Covariance::Options cov_options;
    cov_options.algorithm_type = ceres::CovarianceAlgorithmType::SPARSE_QR; //
    //cov_options.num_threads = 1;
    
    vector<double> uncertainties( num_fit_pars, 0.0 );
    double *uncertainties_ptr = uncertainties.data();
    
    
    ceres::Covariance covariance(cov_options);
    vector<pair<const double*, const double*> > covariance_blocks;
    covariance_blocks.push_back( make_pair( pars, pars) );
    
    if( !covariance.Compute(covariance_blocks, &problem) )
    {
      cerr << "Failed to compute covariance!" << endl;
      uncertainties_ptr = nullptr;
    }else
    {
      vector<double> row_major_covariance( num_fit_pars * num_fit_pars );
      const vector<const double *> const_par_blocks( 1, pars );

      const bool success = covariance.GetCovarianceMatrix( const_par_blocks, row_major_covariance.data() );
      assert( success );
      if( success )
      {
        //Obviously we dont need to recopy things, just to get the diagonal, but eventually we will save the full matrix....
        vector<vector<double>> covariance_matrix( num_fit_pars, vector<double>(num_fit_pars, 0.0) );

        for( size_t row = 0; row < num_fit_pars; ++row )
        {
          for( size_t col = 0; col < num_fit_pars; ++col )
            covariance_matrix[row][col] = row_major_covariance[row*num_fit_pars + col];
        }//for( size_t row = 0; row < num_fit_pars; ++row )

        for( size_t i = 0; i < num_fit_pars; ++i )
        {
          // TODO: compare uncertainties with Minuit method - should maybe check out.
          if( covariance_matrix[i][i] > 0.0 )
            uncertainties[i] = sqrt( covariance_matrix[i][i] );
        }
      }else
      {
        uncertainties_ptr = nullptr;
      }//
    }//if( we failed to computer covariance ) / else
    
    
    vector<double> residuals( cost_functor->number_residuals(), 0.0 );
    auto final_peaks = cost_functor->parametersToPeaks<PeakDef,double>( &parameters[0], uncertainties_ptr, residuals.data() );

    results.resize( final_peaks.size() );
    for( size_t i = 0; i < final_peaks.size(); ++i )
    {
      chi2Dof = final_peaks[i].chi2dof();
      results[i] = make_shared<PeakDef>( final_peaks[i] );
    }
    
    // One example (two peaks) with 46 iterations had 2529 to solve, and 17 calls to get covariance.
    cout << "Took " << cost_functor->m_ncalls.load() << " calls to get covariance." << endl;
  }catch( std::exception &e )
  {
    cout << "fit_peak_for_user_click_LM caught: " << e.what() << endl;
    //assert( 0 );
  }//try / catch
}//void fit_peak_for_user_click_LM(...)


}//namespace PeakFitLM
