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

#include <chrono>
#include <string>
#include <fstream>
#include <iostream>


//Roots Minuit2 includes
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnScan.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnSimplex.h"
#include "Minuit2/MinosError.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/CombinedMinimizer.h"
#include "Minuit2/SimplexMinimizer.h"

#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp> // For boost::asio::post

#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/SpecUtilsAsync.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/DetectorPeakResponse.h"

#include <Wt/WString>

#include "PeakFitImprove.h"
#include "FinalFit_GA.h"
#include "InitialFit_GA.h"

using namespace std;



namespace
{
std::function<double(const FinalPeakFitSettings &)> ns_ga_eval_fcn;

bool sm_has_been_called = false;

std::atomic<size_t> ns_generation_num( 0 );
std::atomic<size_t> ns_individuals_eval( 0 );


/** Fit the ROIs in parralele, using a SpecUtilsAsync::ThreadPool, in `final_peak_fit(...)`*/
const bool ns_final_peak_fit_parrallel = false;




class FitFinalPeakFitSettingsChi2
: public ROOT::Minuit2::FCNBase
{
  std::function<double( const FinalPeakFitSettings &)> m_eval_fcn;


public:

  FitFinalPeakFitSettingsChi2( std::function<double( const FinalPeakFitSettings &)> eval_fcn )
    : m_eval_fcn( eval_fcn )
  {
  }

  virtual double Up() const { return 1.0; }
  size_t nfitPars() const { return 1; }

  static FinalPeakFitSettings params_to_settings( const std::vector<double> &params )
  {
    assert( params.size() == 15 );
    if( params.size() != 15 )
      throw std::runtime_error( "params_to_settings: invalid number of parameters." );

    FinalPeakFitSettings settings;

    size_t index = 0;
    settings.require_combine_num_fwhm_near = params[index++];
    settings.not_allow_combine_num_fwhm_near = params[index++];
    //settings.combine_ROI_overlap_frac = params[index++];
    settings.cont_type_peak_nsigma_threshold = params[index++];
    settings.cont_type_left_right_nsigma = params[index++];
    settings.cont_poly_order_increase_chi2dof_required = params[index++];
    settings.cont_step_type_increase_chi2dof_required = params[index++];
    settings.skew_nsigma = params[index++];
    settings.left_residual_sum_min_to_try_skew = params[index++];
    settings.right_residual_sum_min_to_try_skew = params[index++];
    settings.skew_improve_chi2_dof_threshold = params[index++];
    settings.roi_extent_low_num_fwhm_base_highstat = params[index++];
    settings.roi_extent_high_num_fwhm_base_highstat = params[index++];
    settings.roi_extent_low_num_fwhm_base_lowstat = params[index++];
    settings.roi_extent_high_num_fwhm_base_lowstat = params[index++];
    settings.high_stat_threshold = params[index++];
    settings.roi_extent_low_num_fwhm_extra = params[index++];
    settings.roi_extent_high_num_fwhm_extra = params[index++];
    settings.roi_end_second_deriv_thresh = params[index++];
    settings.break_multi_roi_up_continuum_away_sigma = params[index++];
    settings.break_multi_roi_up_required_chi2dof_improve = params[index++];

    assert( index == params.size() );

    return settings;
  }

  virtual double operator()( const std::vector<double> &params ) const
  {
    const FinalPeakFitSettings settings = params_to_settings( params );

    try
    {
      return m_eval_fcn( settings );
    }catch( std::exception &e )
    {
      std::cerr << "FitFinalPeakFitSettingsChi2::operator() caught: " << e.what() << endl;
    }

    return std::numeric_limits<double>::max();
  }//operator()
};//class FitFinalPeakFitSettingsChi2


}//namespace


namespace FinalFit_GA
{

std::string FinalFitScore::print( const std::string &varname ) const
{
  string answer;
  answer += varname + ".area_score =                        " + std::to_string(area_score) + "\n";
  answer += varname + ".width_score =                       " + std::to_string(width_score) + "\n";
  answer += varname + ".position_score =                    " + std::to_string(position_score) + "\n";
  answer += varname + ".ignored_unexpected_peaks =          " + std::to_string(ignored_unexpected_peaks) + "\n";
  answer += varname + ".unexpected_peaks_sum_significance = " + std::to_string(unexpected_peaks_sum_significance) + "\n";
  answer += varname + ".num_peaks_used =                    " + std::to_string(num_peaks_used) + "\n";
  answer += varname + ".total_weight =                      " + std::to_string(total_weight) + "\n";
  return answer;
}


vector<PeakDef> final_peak_fit_for_roi( const vector<PeakDef> &pre_fit_peaks,
                                       const FinalPeakFitSettings &final_fit_settings,
                                       const bool isHPGe,
                                       const shared_ptr<const SpecUtils::Measurement> &data )
{
#ifndef NDEBUG
  // Make sure all peaks are from the same ROI - only need to do in debug builds
  for( size_t i = 1; i < pre_fit_peaks.size(); ++i )
  {
    assert( pre_fit_peaks[0].continuum() == pre_fit_peaks[i].continuum() );
  }
#endif

  assert( data );
  if( !data )
    return vector<PeakDef>();

  assert( !pre_fit_peaks.empty() );
  if( pre_fit_peaks.empty() )
    return vector<PeakDef>();

  const shared_ptr<const vector<float>> &channel_counts_ptr = data->gamma_channel_contents();
  const std::shared_ptr< const std::vector<float> > &energies_ptr = data->gamma_channel_energies();
  assert( channel_counts_ptr && !channel_counts_ptr->empty() && energies_ptr && !energies_ptr->empty() );
  if( !channel_counts_ptr || channel_counts_ptr->empty() || !energies_ptr || energies_ptr->empty() )
    return vector<PeakDef>{};
  const vector<float> &channel_counts = *channel_counts_ptr;
  const vector<float> &energies = *energies_ptr;

  const std::shared_ptr<const PeakContinuum> orig_continuum = pre_fit_peaks.front().continuum();

  if( 0.5*(pre_fit_peaks.front().fwhm() + pre_fit_peaks.back().fwhm()) <= 0.0 )
    return vector<PeakDef>();

  const bool is_high_stat = ([&final_fit_settings, &pre_fit_peaks](){
    for( const PeakDef &p : pre_fit_peaks )
    {
      if( (p.amplitude() > 0.0) && (sqrt(p.amplitude()) > final_fit_settings.high_stat_threshold) )
        return true;
    }
    return false;
  })();

  const bool lower_extend_base_nfwhm = is_high_stat ? final_fit_settings.roi_extent_low_num_fwhm_base_highstat
                                                    : final_fit_settings.roi_extent_low_num_fwhm_base_lowstat;
  const bool upper_extend_base_nfwhm = is_high_stat ? final_fit_settings.roi_extent_high_num_fwhm_base_highstat
                                                    : final_fit_settings.roi_extent_high_num_fwhm_base_lowstat;



  const auto refit_peaks = [data, isHPGe]( const vector<PeakDef> &these_input_peaks ) -> vector<PeakDef> {

    // Do an initial fit
    const double initial_stat_threshold = 0.0;
    const double initial_hypothesis_threshold = 0.0;
    const bool amplitudeOnly = false;

    vector<PeakDef> these_fit_peaks;
#if( USE_LM_PEAK_FIT )
      vector<shared_ptr<const PeakDef>> results_tmp, input_peaks_tmp;
      for( const PeakDef &p : these_input_peaks )
        input_peaks_tmp.push_back( make_shared<PeakDef>(p) );
      PeakFitLM::fit_peaks_LM( results_tmp, input_peaks_tmp, data,
                          initial_stat_threshold, initial_hypothesis_threshold,  amplitudeOnly, isHPGe );
    for( const shared_ptr<const PeakDef> &p : results_tmp )
      these_fit_peaks.push_back( *p );
#else
      fitPeaks( these_fit_peaks, initial_stat_threshold, initial_hypothesis_threshold,
                data, initial_peaks, amplitudeOnly, isHPGe );
#endif

    return these_fit_peaks;
  };

  const bool debug_peak = ((pre_fit_peaks.front().mean() > PeakFitImprove::debug_lower_energy && pre_fit_peaks.front().mean() < PeakFitImprove::debug_upper_energy)
                           || (pre_fit_peaks.back().mean() > PeakFitImprove::debug_lower_energy && pre_fit_peaks.back().mean() < PeakFitImprove::debug_upper_energy));

  if( debug_peak )
  {
    cout << "Final fit debugging peaks {";
    for( const auto &p : pre_fit_peaks )
      cout << p.mean() << ", ";
    cout << "}" << endl;
  }

  vector<PeakDef> initial_peaks = pre_fit_peaks;

  const auto find_ROI_limits = [&initial_peaks,&data,&final_fit_settings,&channel_counts,isHPGe,lower_extend_base_nfwhm,upper_extend_base_nfwhm, debug_peak]() -> pair<double,double> {

    const PeakDef &first_peak = initial_peaks.front();
    const PeakDef &last_peak = initial_peaks.back();

    const double first_mean = first_peak.mean();
    const double last_mean = last_peak.mean();
    const double first_fwhm = first_peak.fwhm();
    const double last_fwhm = last_peak.fwhm();
    const double avrg_fwhm = 0.5*(first_fwhm + last_fwhm);


    double nchannel_half_fwhm = -1.0;
    double effective_fwhm = avrg_fwhm;

    // We will force the FWHM to be at least 2.4 channels.
    {// Begin check `effective_fwhm` and compute `nchannel_half_fwhm`
      const double lower_roi_channel = data->energy_calibration()->channel_for_energy( first_peak.lowerX() );
      const double upper_roi_channel = data->energy_calibration()->channel_for_energy( first_peak.upperX() );
      const double num_roi_channels = upper_roi_channel - lower_roi_channel;
      const double num_roi_keV = first_peak.upperX() - first_peak.lowerX();
      double num_roi_fwhm = num_roi_keV / effective_fwhm;
      double num_channel_per_fwhm = num_roi_channels / num_roi_fwhm;

      const double min_channel_per_fwhm = 2.4; //2.4 picked fairly arbitrarily
      if( num_channel_per_fwhm < min_channel_per_fwhm )
      {
        effective_fwhm *= (min_channel_per_fwhm / num_channel_per_fwhm);
        num_roi_fwhm = num_roi_keV / effective_fwhm;
        num_channel_per_fwhm = num_roi_channels / num_roi_fwhm;
      }

      nchannel_half_fwhm = 0.5*num_channel_per_fwhm;
    }// End check `effective_fwhm` and compute `nchannel_half_fwhm`

    assert( nchannel_half_fwhm > 1.0 );

    const double base_roi_lower = first_mean - (lower_extend_base_nfwhm * effective_fwhm);
    const double base_roi_upper = last_mean + (upper_extend_base_nfwhm * effective_fwhm);

    const double min_roi_lower = base_roi_lower - effective_fwhm * final_fit_settings.roi_extent_low_num_fwhm_extra;
    int min_roi_lower_channel = static_cast<int>( data->find_gamma_channel( static_cast<float>(min_roi_lower) ) );

    const double max_roi_upper = base_roi_upper + effective_fwhm * final_fit_settings.roi_extent_high_num_fwhm_extra;
    int max_roi_upper_channel = static_cast<int>( data->find_gamma_channel( static_cast<float>(max_roi_upper) ) );

    assert( max_roi_upper_channel >= min_roi_lower_channel );

    const int min_possible_roi_channels = 20;
    const int initial_max_roi_nchannels = 1 + max_roi_upper_channel - min_roi_lower_channel;
    if( debug_peak )
      cout << "Initial ROI Width = " << initial_max_roi_nchannels << " channels; min_roi_lower=" << min_roi_lower
      << ", max_roi_upper=" << max_roi_upper << ", Type=" << PeakContinuum::offset_type_str(first_peak.continuum()->type()) << endl;
    if( initial_max_roi_nchannels < min_possible_roi_channels  )
    {
      const int num_ch_add = static_cast<int>( std::round(0.5*(min_possible_roi_channels - initial_max_roi_nchannels)) );
      if( debug_peak )
        cout << "num_ch_add=" << num_ch_add << endl;
      min_roi_lower_channel = std::max( 0, min_roi_lower_channel - num_ch_add );
      max_roi_upper_channel = std::min( static_cast<int>(channel_counts.size() - 1), max_roi_upper_channel + num_ch_add );
    }

    const int smooth_nchannel = std::max( 3, static_cast<int>(std::round(nchannel_half_fwhm)) );

    // Sometimes the input peak is artificially narrow, and the data is still dropping away at
    //  min_roi_lower_channel/max_roi_upper_channel, so we will override this, and allow for the ROI to expand beyond
    //  these values if the data is still dropping.
    const int min_num_extra_channel = 15; // Totally arbitrary
    const int base_roi_lower_ch = static_cast<int>( data->find_gamma_channel( static_cast<float>(base_roi_lower) ) );
    const int sanity_smooth_lower = std::max(base_roi_lower_ch - min_num_extra_channel, 0);
    const int smooth_lower_index = std::min( sanity_smooth_lower, std::max( min_roi_lower_channel - static_cast<int>(std::round(2.0*nchannel_half_fwhm)), 0 ) );

    const int base_roi_upper_ch = static_cast<int>( data->find_gamma_channel( static_cast<float>(base_roi_upper) ) );
    const int sanity_smooth_upper = std::max( base_roi_upper_ch + min_num_extra_channel, static_cast<int>(channel_counts.size()) - 1 );

    const int smooth_upper_index = std::max( sanity_smooth_upper, std::min( max_roi_upper_channel + static_cast<int>(std::round(2.0*nchannel_half_fwhm)),
                                            static_cast<int>(channel_counts.size()) - 1 ) );
    const int num_smooth_channel = 1 + smooth_upper_index - smooth_lower_index;

    vector<float> possible_roi_channels( num_smooth_channel );
    for( int i = 0; i < num_smooth_channel; ++i )
    {
      assert( static_cast<size_t>(i + smooth_lower_index) < channel_counts.size() );
      possible_roi_channels[i] = channel_counts[i + smooth_lower_index];
    }

    const int order = isHPGe ? 3 : 2;

    SavitzyGolayCoeffs first_deriv_sgcoeffs( smooth_nchannel, smooth_nchannel, order, 1 );
    vector<float> smoothed_first_derivs, first_derivs_variance;
    first_deriv_sgcoeffs.smooth_with_variance( possible_roi_channels, smoothed_first_derivs, first_derivs_variance );

    SavitzyGolayCoeffs second_deriv_sgcoeffs( smooth_nchannel, smooth_nchannel, order, 2 );
    vector<float> smoothed_second_derivs, second_derivs_variance;
    second_deriv_sgcoeffs.smooth_with_variance( possible_roi_channels, smoothed_second_derivs, second_derivs_variance );


    int roi_lower_channel = static_cast<int>( data->find_gamma_channel( static_cast<float>(base_roi_lower) ) );

    if( debug_peak )
      cout << "mean=" << first_mean << "/" << last_mean << ", nchannel_half_fwhm=" << nchannel_half_fwhm
      << " - ROI=[" << first_peak.lowerX() << ", " << first_peak.upperX() << "]" << endl;


    if( debug_peak )
      cout << "Base lower channel " << roi_lower_channel << ", " << base_roi_lower << " keV" << endl;

    // While derivative is positive, or not very negative, keep going away
    //  Note: this may go past `min_roi_lower_channel`
    for( int channel = roi_lower_channel; channel > 0; channel -= 1 )
    {
      const int index = channel - smooth_lower_index;
      if( index <= 0 )
        break;

      const double sigma = sqrt( first_derivs_variance[index] );
      const double significance = (first_derivs_variance[index] > 0.0) ? (smoothed_first_derivs[index] / sigma) : 0.0;

      const double next_sigma = sqrt( first_derivs_variance[index - 1] );
      const double next_significance = (first_derivs_variance[index - 1] > 0.0) ? (smoothed_first_derivs[index - 1] / next_sigma) : 0.0;

      if( debug_peak )
        cout << "  First derivative (low-side) at channel " << channel << " is deriv=" << smoothed_first_derivs[index]
        << " var=" << sigma << ", sig=" << significance
        << " next_deriv=" << smoothed_first_derivs[index - 1] << ", next_var=" << next_sigma << ", next_sig=" << next_significance
        << endl;


      if( (significance < 0.0) && (next_significance < 0.0) )
      {
        if( debug_peak )
          cout << "  stopping first derivative decent on low-side at channel " << channel << endl;

        // The derivative is basically the previous bin to the current bin, so we will decrement one channel so the ROI will include the lowest count channel
        if( (smoothed_first_derivs[index + 1] > 0.0) && (channel > min_roi_lower_channel) )
          roi_lower_channel = channel - 1;

        break;
      }

      roi_lower_channel = channel;
    }

    if( debug_peak )
      cout << "Derivative no longer positive on low side for channel " << roi_lower_channel << ", " << data->gamma_channel_upper(roi_lower_channel) << " keV" << endl;

    // Now we will sum the second derivative significance over ~0.5 FWHM to try and detect features
    for( int channel = roi_lower_channel; channel > min_roi_lower_channel; channel -= 1 )
    {
      const int index = channel - smooth_lower_index;
      const double sigma = sqrt( second_derivs_variance[index] );
      const double significance = (second_derivs_variance[index] > 0.0) ? (smoothed_second_derivs[index] / sigma) : 0.0;

      if( significance <= 0.0 )
      {
        roi_lower_channel = channel;
        continue;
      }

      // Now we'll sum over the next nchannel_half_fwhm channels, and if
      double sum_d2N = 0.0, sum_variance = 0.0;
      const int end_index = std::max( 0, index );
      const int start_index = std::max( 0, index - static_cast<int>(std::round(nchannel_half_fwhm)) );
      for( size_t j = start_index; j < end_index; ++j )
      {
        sum_d2N += smoothed_second_derivs[j];
        sum_variance += second_derivs_variance[j];
      }

      const double integrated_sigma = (sum_variance > 0.0) ? std::sqrt( sum_variance ) : 1.0;
      const double integrated_significance = sum_d2N / integrated_sigma;
      if( fabs(integrated_significance) > final_fit_settings.roi_end_second_deriv_thresh )
      {
        if( debug_peak )
          cout << "Detected feature on low-side at channel " << channel << ", " << data->gamma_channel_upper(channel) << " keV - significance=" << integrated_significance << endl;
        break;
      }
      roi_lower_channel = channel;
    }

    if( debug_peak )
      cout << "Second deriv feature det on low side stopped on channel " << roi_lower_channel << ", "
      << data->gamma_channel_lower(roi_lower_channel) << " keV"
      << " (min_roi_lower_channel=" << min_roi_lower_channel << ", base_roi_lower=" << base_roi_lower << ", min_roi_lower=" << min_roi_lower
      << ", first_fwhm=" << first_fwhm << ", effective_fwhm=" << effective_fwhm << ", nchannel_half_fwhm=" << nchannel_half_fwhm <<")"
      << endl;


    int roi_upper_channel = static_cast<int>( data->find_gamma_channel( static_cast<float>(base_roi_upper) ) );
    if( debug_peak )
      cout << "Base upper channel " << roi_upper_channel << ", " << base_roi_upper << " keV" << endl;

    // While derivative is negative, or not very positive, keep going away from peak mean
    //  Note: this may go past `max_roi_upper_channel`
    for( int channel = roi_upper_channel; channel < static_cast<int>(channel_counts.size() - 1); channel += 1 )
    {
      const int index = channel - smooth_lower_index;
      if( index >= num_smooth_channel )
        break;

      const double sigma = sqrt( first_derivs_variance[index] );
      const double significance = (first_derivs_variance[index] > 0.0) ? (smoothed_first_derivs[index] / sigma) : 0.0;

      const double next_sigma = sqrt( first_derivs_variance[index - 1] );
      const double next_significance = (first_derivs_variance[index - 1] > 0.0) ? (smoothed_first_derivs[index - 1] / next_sigma) : 0.0;

      if( (significance > 0.0) && (next_significance > 0.0) )
      {
        break;
      }

      roi_upper_channel = channel;
    }

    if( debug_peak )
      cout << "Derivative no longer negative on high side for channel " << roi_upper_channel << ", "
      << data->gamma_channel_upper(roi_upper_channel) << " keV"
      << ", smoothed_second_derivs=" << smoothed_second_derivs[roi_upper_channel-smooth_lower_index]
      << ", second_derivs_variance=" << second_derivs_variance[roi_upper_channel-smooth_lower_index]
      << endl;

    // Now we will sum the second derivative significance over ~0.5 FWHM to try and detect features
    for( int channel = roi_upper_channel; channel < max_roi_upper_channel; channel += 1 )
    {
      const int index = channel - smooth_lower_index;
      const double sigma = sqrt( second_derivs_variance[index] );
      const double significance = (second_derivs_variance[index] > 0.0) ? (smoothed_second_derivs[index] / sigma) : 0.0;

      if( debug_peak )
        cout << "  channel=" << channel << ", significanc=" << significance << ", " << smoothed_second_derivs[index] << " / " << sigma << endl;

      if( significance <= 0.0 )
      {
        roi_upper_channel = channel;
        continue;
      }

      // Now we'll sum over the next nchannel_half_fwhm channels, and if
      double sum_d2N = 0.0, sum_variance = 0.0;
      const int start_index = index + 1;
      const int end_index = std::min( start_index + static_cast<int>(std::round(nchannel_half_fwhm)) + 1, static_cast<int>(smoothed_second_derivs.size()) );
      for( size_t j = start_index; j < end_index; ++j )
      {
        sum_d2N += smoothed_second_derivs[j];
        sum_variance += second_derivs_variance[j];
      }

      const double integrated_sigma = (sum_variance > 0.0) ? std::sqrt( sum_variance ) : 1.0;
      const double integrated_significance = sum_d2N / integrated_sigma;
      if( debug_peak )
        cout << "  integrated_significance=" << integrated_significance << ", roi_end_second_deriv_thresh=" << final_fit_settings.roi_end_second_deriv_thresh << endl;

      if( fabs(integrated_significance) > final_fit_settings.roi_end_second_deriv_thresh )
      {
        if( debug_peak )
          cout << "Detected feature on high-side at channel " << channel << ", " << data->gamma_channel_lower(channel) << " keV - significance=" << integrated_significance << endl;
        break;
      }
      roi_upper_channel = channel;
    }

    if( debug_peak )
      cout << "After second deriv feature det on high side stopped on channel " << roi_upper_channel << ", " << data->gamma_channel_upper(roi_upper_channel) << " keV" << endl;

    //cout << "first_mean=" << first_mean << ", last_mean=" << last_mean << endl;
    /*
    for( int i = 0; i < num_smooth_channel; ++i )
    {
      const size_t channel = static_cast<size_t>(i + smooth_lower_index);
      const double sigma = sqrt( second_derivs_variance[i] );
      const double significance = (second_derivs_variance[i] > 0.0) ? (smoothed_second_derivs[i] / sigma) : 0.0;
      const double lenergy = data->gamma_channel_lower(channel);
      const double renergy = data->gamma_channel_upper(channel);
      const double energy = 0.5*(lenergy + renergy);

      const double nfwhm_away = (energy < first_mean) ? ((energy - first_mean) / first_fwhm) : ((energy - last_mean) / last_fwhm);
      if( debug_peak )
        cout << "Channel " << channel << " significance=" << significance
        << ", e=[" << energy << ", fwhm_away=" << nfwhm_away << "]"
        << endl;
    }
     */



    double roi_lower = data->gamma_channel_lower( static_cast<size_t>(roi_lower_channel) );
    double roi_upper = data->gamma_channel_upper( static_cast<size_t>(roi_upper_channel) ) - 0.001;

    // The original ROI definition (based off smoothed second derivative of data) is actually pretty close to
    //  what we want "by eye" - maybe a touch shorter than my preference, but reasonable.
    // If things look reasonable, and `roi_lower` or `roi_upper` is past the original ROI defintion, lets average
    //  average things
    const double orig_lower_energy = first_peak.lowerX();
    const double orig_upper_energy = first_peak.upperX();

    if( (orig_lower_energy < base_roi_lower) && (orig_lower_energy > min_roi_lower) && (roi_lower < orig_lower_energy) )
    {
      const double lower_avrg = 0.5*(orig_lower_energy + roi_lower);
      const size_t lower_avrg_channel = data->find_gamma_channel( lower_avrg );
      //cout << "Moving roi_lower=" << roi_lower << " ---> " << data->gamma_channel_lower( lower_avrg_channel ) << endl;

      roi_lower_channel = static_cast<int>(lower_avrg_channel);
      roi_lower = data->gamma_channel_lower( lower_avrg_channel );
    }

    if( (orig_upper_energy > base_roi_upper) && (orig_upper_energy < max_roi_upper) && (roi_upper > orig_upper_energy) )
    {
      const double upper_avrg = 0.5*(orig_upper_energy + roi_upper);
      const size_t upper_avrg_channel = data->find_gamma_channel( upper_avrg );
      //cout << "Moving roi_upper=" << roi_upper << " ---> " << data->gamma_channel_lower( upper_avrg_channel ) << endl;

      roi_upper_channel = static_cast<int>(upper_avrg_channel);
      roi_upper = data->gamma_channel_upper( upper_avrg_channel ) - 0.001;
    }

#if( PERFORM_DEVELOPER_CHECKS )
    const size_t test_roi_lower_ch = data->find_gamma_channel(static_cast<float>(roi_lower));
    const size_t test_roi_upper_ch = data->find_gamma_channel(static_cast<float>(roi_upper));
    assert( test_roi_lower_ch == static_cast<size_t>(roi_lower_channel) );
    assert( test_roi_upper_ch == static_cast<size_t>(roi_upper_channel) );
#endif
    return make_pair(roi_lower, roi_upper);
  };//find_ROI_limits lambda


  pair<double,double> roi_limits = find_ROI_limits();

  const double max_roi_width_change_without_refit = 0.2; //totally arbitrary
  const double orig_range = (pre_fit_peaks.front().upperX() - pre_fit_peaks.front().lowerX());
  if( debug_peak )
    cout << "orig_range=" << orig_range << ", updated_range=" << (roi_limits.second - roi_limits.first) << ", first_fwhm=" << initial_peaks.front().fwhm()
    << ", new ROI=[" << roi_limits.first << ", " << roi_limits.second << "]"
    << endl;

  if( fabs(orig_range - (roi_limits.second - roi_limits.first)) > max_roi_width_change_without_refit*orig_range )
  {
    shared_ptr<PeakContinuum> updated_continuum = make_shared<PeakContinuum>( *orig_continuum );
    updated_continuum->setRange( roi_limits.first, roi_limits.second );

    try
    {
      vector<PeakDef> updated_peaks = refit_peaks( initial_peaks );
      if( !updated_peaks.empty() )
        initial_peaks.swap( updated_peaks );
    }catch( std::exception &e )
    {
      if( debug_peak )
        cerr << "Initial peak refit failed: " << e.what() << endl;
    }

    roi_limits = find_ROI_limits();

    if( debug_peak )
      cout << "updated updated_range=" << (roi_limits.second - roi_limits.first) << " with first_fwhm=" << initial_peaks.front().fwhm()
      << " and extents [" << roi_limits.first << ", " << roi_limits.second << "] keV" << endl;
  }//if( fabs(orig_range - (roi_limits.second - roi_limits.first)) > 0.2*orig_range )

  const double roi_lower = roi_limits.first;
  const double roi_upper = roi_limits.second;

  if( debug_peak )
    cout << "  updated trial ROI=[" << roi_lower << ", " << roi_upper << "]" << endl;

  shared_ptr<PeakContinuum> initial_continuum = make_shared<PeakContinuum>( *orig_continuum );
  initial_continuum->setRange( roi_lower, roi_upper );
  for( PeakDef &p : initial_peaks )
    p.setContinuum( initial_continuum );

  const size_t original_roi_lower_channel = data->find_gamma_channel( orig_continuum->lowerEnergy() );
  const size_t original_roi_upper_channel = data->find_gamma_channel( orig_continuum->upperEnergy() );

  const size_t updated_roi_lower_channel = data->find_gamma_channel( roi_lower );
  const size_t updated_roi_upper_channel = data->find_gamma_channel( roi_upper );

  // For Chi2, we will use the intersection of original and updated ROIs...
  const size_t compare_lower_channel = std::max( original_roi_lower_channel, updated_roi_lower_channel );
  const size_t compare_upper_channel = std::min( original_roi_upper_channel, updated_roi_upper_channel );


 const auto chi2_for_region_wrapper = []( vector<PeakDef> peaks, const std::shared_ptr<const SpecUtils::Measurement> &data, const size_t xlowbin, const size_t xhighbin ) -> double {
  vector<shared_ptr<const PeakDef>> peaks_tmp;
  for( const PeakDef &p : peaks )
    peaks_tmp.push_back( make_shared<PeakDef>(p) );
  return chi2_for_region( peaks_tmp, data, static_cast<int>(xlowbin), static_cast<int>(xhighbin) );
 };

  const double intial_chi2 = chi2_for_region_wrapper( initial_peaks, data, compare_lower_channel, compare_upper_channel );

  vector<PeakDef> first_fit_peaks = refit_peaks( initial_peaks );

  if( first_fit_peaks.empty() )
    throw std::runtime_error( "final_peak_fit_for_roi: first_fit_peaks is empty" );

  bool peak_wider_than_roi = false;
  for( const auto &p : first_fit_peaks )
    peak_wider_than_roi |= (p.fwhm() > (p.upperX() - p.lowerX()));

  if( peak_wider_than_roi )
  {
    if( debug_peak )
      cout << "Swapping orig peaks back because first fit made FWHM->" << first_fit_peaks.front().fwhm()
      << " (from " << initial_peaks.front().fwhm()
      << ") while ROI is only [" << first_fit_peaks.front().lowerX() << ", " << first_fit_peaks.front().upperX() << endl;
    first_fit_peaks = initial_peaks;
  }

  const double initial_fit_chi2 = chi2_for_region_wrapper( first_fit_peaks, data, compare_lower_channel, compare_upper_channel );
  const double approx_dof = (compare_upper_channel - compare_lower_channel);

  if( debug_peak || PeakFitImprove::debug_printout )
    cout << "Initial ROI [" << roi_lower << ", " << roi_upper << "] fit chi2: " << initial_fit_chi2 << " compared to intial_chi2: " << intial_chi2
    << " - new mean=" << first_fit_peaks.front().mean() << " and fwhm=" << first_fit_peaks.front().fwhm() << endl;

  // We will actually take the new ROI width, no matter what, since we are optimizing things.

  bool explore_cont_type = false;
  for( const PeakDef &p : pre_fit_peaks )
  {
    switch( p.continuum()->type() )
    {
      case PeakContinuum::OffsetType::External:    // We wont be tuning external continuums, although we shouldnt actually get here
        assert( false );

      case PeakContinuum::OffsetType::FlatStep:    // We'll check for stepped continuums later
      case PeakContinuum::OffsetType::LinearStep:
      case PeakContinuum::OffsetType::BiLinearStep:
      case PeakContinuum::OffsetType::Cubic:        // We cant go any higher than cubic
        continue;

      // The current continuum types that we will now try other types to check for improvements.
      case PeakContinuum::OffsetType::NoOffset:
      case PeakContinuum::OffsetType::Constant:
      case PeakContinuum::OffsetType::Linear:
      case PeakContinuum::OffsetType::Quadratic:
      break;
    }//switch( p.continuum()->type() )

    const double peak_uncert = p.amplitudeUncert() > 0.0 ? p.amplitudeUncert() : sqrt(p.amplitude());
    const double nsigma = p.amplitude() / peak_uncert;
    explore_cont_type |= (nsigma > final_fit_settings.cont_type_peak_nsigma_threshold);
  }//for( const PeakDef &p : pre_fit_peaks )

  vector<PeakContinuum::OffsetType> cont_types_to_check;
  if( explore_cont_type )
  {
    for( PeakContinuum::OffsetType type = PeakContinuum::OffsetType(orig_continuum->type() + 1);
         type <= PeakContinuum::OffsetType::Cubic;
         type = PeakContinuum::OffsetType(type + 1) )
    {
      cont_types_to_check.push_back( type );
    }
  }//if( explore_cont_type )


  bool explore_stepped_cont = false;
  const bool starting_with_step_continuum = PeakContinuum::is_step_continuum( orig_continuum->type() );
  if( starting_with_step_continuum )
  {// Begin of check if we should explore stepped cont
    // We will sum data starting 1.5 FWHM below first peak, and 1.5 FWHM above last peak, for 0.5 FWHM each.
    const PeakDef &first_peak = pre_fit_peaks.front();
    const PeakDef &last_peak = pre_fit_peaks.back();
    const double first_peak_lower = first_peak.mean() - 1.5*first_peak.fwhm();
    const double roi_lower = first_peak.mean() - 2.0*first_peak.fwhm();
    const double roi_upper = last_peak.mean() + 2.0*last_peak.fwhm();
    const double last_peak_upper = last_peak.mean() + 1.5*last_peak.fwhm();

    const size_t roi_lower_channel = data->find_gamma_channel( roi_lower );
    const size_t roi_upper_channel = data->find_gamma_channel( roi_upper );
    const size_t first_peak_lower_channel = data->find_gamma_channel( first_peak_lower );
    const size_t last_peak_upper_channel = data->find_gamma_channel( last_peak_upper );
    assert( roi_lower_channel <= first_peak_lower_channel );
    assert( last_peak_upper_channel <= roi_upper_channel );
    const size_t num_channels_to_use = std::max( size_t(2), std::max( first_peak_lower_channel - roi_lower_channel, roi_upper_channel - last_peak_upper_channel ) );

    double lower_sum = 0.0, upper_sum = 0.0;
    for( size_t i = roi_lower_channel;  i < (roi_lower_channel + num_channels_to_use); ++i )
      lower_sum += data->gamma_channel_content(i);
    for( size_t i = last_peak_upper_channel; i < (last_peak_upper_channel + num_channels_to_use); ++i )
      upper_sum += data->gamma_channel_content(i);

    const double average_lower_sum = lower_sum / num_channels_to_use;
    const double average_upper_sum = upper_sum / num_channels_to_use;
    const double average_lower_uncert = sqrt( average_lower_sum );
    const double average_upper_uncert = sqrt( average_upper_sum );
    const double average_uncert = sqrt( average_lower_sum + average_upper_sum );
    if( average_uncert > 0.0 )
    {
      const double average_nsigma = (average_lower_sum - average_upper_sum) / average_uncert;
      explore_stepped_cont = (average_nsigma >= final_fit_settings.cont_type_left_right_nsigma);

      if( debug_peak )
        cout << "average_nsigma=" << average_nsigma << ", explore_stepped_cont=" << explore_stepped_cont << ", final_fit_settings.cont_type_left_right_nsigma=" << final_fit_settings.cont_type_left_right_nsigma << endl;
    }
  }// End of check if we should explore stepped cont

  if( explore_stepped_cont )
  {
    cont_types_to_check.push_back( PeakContinuum::OffsetType::FlatStep );
    cont_types_to_check.push_back( PeakContinuum::OffsetType::LinearStep );

    // I dont find bilinear steps to be very good, at least visually, so we will not explore them
    //cont_types_to_check.push_back( PeakContinuum::OffsetType::BiLinearStep );
  }else if( starting_with_step_continuum )
  {
    if( orig_continuum->type() != PeakContinuum::OffsetType::FlatStep )
      cont_types_to_check.push_back( PeakContinuum::OffsetType(orig_continuum->type() - 1) );

    if( orig_continuum->type() != PeakContinuum::OffsetType::BiLinearStep )
      cont_types_to_check.push_back( PeakContinuum::OffsetType(orig_continuum->type() + 1) );
  }

  vector<tuple<PeakContinuum::OffsetType,vector<PeakDef>,double>> trial_peak_fits;
  trial_peak_fits.push_back( make_tuple( orig_continuum->type(), first_fit_peaks, initial_fit_chi2/approx_dof )  );

  for( const PeakContinuum::OffsetType type : cont_types_to_check )
  {
    vector<PeakDef> these_input_peaks = first_fit_peaks;
    shared_ptr<PeakContinuum> cont = make_shared<PeakContinuum>( *orig_continuum );
    cont->calc_linear_continuum_eqn( data, cont->referenceEnergy(), cont->lowerEnergy(), cont->upperEnergy(), 3, 3 );
    cont->setType( type );
    for( PeakDef &p : these_input_peaks )
      p.setContinuum( cont );

    vector<PeakDef> these_fit_peaks;
    try
    {
      these_fit_peaks = refit_peaks( these_input_peaks );
    }catch( const std::exception &e )
    {
      // TODO/Note: It doesn look like we can get here, neither `PeakFitLM::fit_peaks_LM(...)` or `fitPeaks(...)` will through
      cout << "Error fitting peaks: " << e.what() << endl;
      continue;
    }//try / catch to fit peaks

    // We could have not fint any peaks; if `fit_peaks_LM(...)` fails, it will return empty results.
    if( these_fit_peaks.empty() )
      continue;

    const double these_fit_chi2 = chi2_for_region_wrapper( these_fit_peaks, data, compare_lower_channel, compare_upper_channel );
    trial_peak_fits.push_back( make_tuple( type, these_fit_peaks, these_fit_chi2/approx_dof ) );
  }//for( const PeakContinuum::OffsetType type : cont_types_to_check )

  // I'm pretty sure trial_peak_fits is sorted by continuum type, but just to be sure, we will sort it
  std::sort( begin(trial_peak_fits), end(trial_peak_fits), []( const auto &a, const auto &b ){
    return static_cast<int>(std::get<0>(a)) < static_cast<int>(std::get<0>(b));
  } );


  // We will seperately track the chosen index, and the best index, to prevent situations like where we
  //  started with a linear continuum, but the improvement to quadratic wasnt great than `cont_poly_order_increase_chi2dof_required`
  //  but the improvement of Cubic over linear would be `cont_poly_order_increase_chi2dof_required`
  size_t chosen_index = 0, best_index = 0;

  bool found_original_type = false;
  for( size_t i = 0; i < trial_peak_fits.size(); ++i )
  {
    if( std::get<0>(trial_peak_fits[i]) == orig_continuum->type() )
    {
      best_index = i;
      chosen_index = i;
      found_original_type = true;
      break;
    }
  }
  assert( found_original_type );

  if( debug_peak )
    cout << "orig_continuum->type()=" << PeakContinuum::offset_type_str(orig_continuum->type())
    << " and nest_index type=" << PeakContinuum::offset_type_str(std::get<0>(trial_peak_fits[chosen_index])) << endl;

  for( size_t i = 0; i < trial_peak_fits.size(); ++i )
  {
    if( i == chosen_index )
      continue;

    const bool prev_was_step = PeakContinuum::is_step_continuum( std::get<0>(trial_peak_fits[best_index]) );
    const bool this_is_step = PeakContinuum::is_step_continuum( std::get<0>(trial_peak_fits[i]) );

    const double prev_chi2 = std::get<2>(trial_peak_fits[best_index]);
    const double this_chi2 = std::get<2>(trial_peak_fits[i]);

    if( !prev_was_step && !this_is_step )
    {
      if( debug_peak )
        cout << "Will check continuum type " << PeakContinuum::offset_type_str(std::get<0>(trial_peak_fits[i]))
        << ", this_chi2=" << this_chi2 << ", prev_chi2=" << prev_chi2
        << endl;


      if( (this_chi2 + final_fit_settings.cont_poly_order_increase_chi2dof_required) < prev_chi2 )
      {
        if( debug_peak )
          cout << "Updating to non-step type from "
          << PeakContinuum::offset_type_str(std::get<0>(trial_peak_fits[chosen_index]))
          << " to "
          << PeakContinuum::offset_type_str(std::get<0>(trial_peak_fits[i]))
          << " because Chi2 went from " << prev_chi2 << " to " << this_chi2 << ", with step required of "
          << final_fit_settings.cont_poly_order_increase_chi2dof_required
          << endl;

        chosen_index = i;
      }
    }else if( prev_was_step && !this_is_step )
    {
      assert( 0 );  // We should never get here
    }else if( (!prev_was_step && this_is_step) || (prev_was_step && this_is_step) )
    {
      if( (this_chi2 + final_fit_settings.cont_step_type_increase_chi2dof_required) < prev_chi2 )
      {
        if( debug_peak )
          cout << "Udating to step type from "
          << PeakContinuum::offset_type_str(std::get<0>(trial_peak_fits[chosen_index]))
          << " to "
          << PeakContinuum::offset_type_str(std::get<0>(trial_peak_fits[i]))
          << " because Chi2 went from " << prev_chi2 << " to " << this_chi2 << ", with step required of "
          << final_fit_settings.cont_step_type_increase_chi2dof_required
          << endl;

        chosen_index = i;
      }
    }

    if( this_chi2 < prev_chi2 )
      best_index = i;

    if( debug_peak && (chosen_index != i) )
      cout << "Didnt update peak from " << PeakContinuum::offset_type_str(std::get<0>(trial_peak_fits[chosen_index]))
      << " to " << PeakContinuum::offset_type_str(std::get<0>(trial_peak_fits[i])) << endl;
  }//for( size_t i = 1; i < trial_peak_fits.size(); ++i )

  assert( std::get<0>(trial_peak_fits[0]) == orig_continuum->type() );

  if( (debug_peak || PeakFitImprove::debug_printout) && (trial_peak_fits.size() > 1) )
  {
    double sum_peak_area = 0.0;
    for( const auto &peak : std::get<1>(trial_peak_fits[chosen_index]) )
      sum_peak_area += peak.peakArea();

    cout << "Post continuum type fit: Chosen index: " << chosen_index << " gives Chi2Dof: " << std::get<2>(trial_peak_fits[chosen_index])
    << " (ContType: " << PeakContinuum::offset_type_str( std::get<0>(trial_peak_fits[chosen_index]) ) << ")"
    << ",\n\tcompared to initial Chi2Dof: " << initial_fit_chi2/approx_dof
    << " (ContType: " << PeakContinuum::offset_type_str( std::get<0>(trial_peak_fits[0]) ) << "),\n\t"
    "All Chi2s: [";
    for( size_t i = 0; i < trial_peak_fits.size(); ++i )
    {
      cout << "\n\t\t{" << PeakContinuum::offset_type_str( std::get<0>(trial_peak_fits[i]) )
      << ", " << std::get<2>(trial_peak_fits[i]) << "}";
    }//for( size_t i = 1; i < trial_peak_fits.size(); ++i )
    cout << "\n\t]\n\tpoly_order_inc_chi2dof_req=" << final_fit_settings.cont_poly_order_increase_chi2dof_required
    << ", and step_inc_chi2dof_req=" << final_fit_settings.cont_step_type_increase_chi2dof_required
    << "\n\tPeaks have area: " << sum_peak_area << " --- [";
    for( const auto &peak : std::get<1>(trial_peak_fits[chosen_index]) )
      cout << peak.peakArea() << ", ";
    cout << "]";
    cout << endl;
  }//if( PeakFitImprove::debug_printout && (trial_peak_fits.size() > 1) )

  const double after_cont_type_chi2dof = std::get<2>(trial_peak_fits[chosen_index]);
  const double after_cont_type_chi2 = after_cont_type_chi2dof * approx_dof;
  const vector<PeakDef> best_post_cont_type_peaks = std::get<1>(trial_peak_fits[chosen_index]);
  vector<PeakDef> answer = best_post_cont_type_peaks;

  double post_cont_max_nsigma = 0.0;
  size_t post_cont_most_sig_peak_index = 0;
  for( size_t i = 0; i < best_post_cont_type_peaks.size(); ++i )
  {
    const PeakDef &p = best_post_cont_type_peaks[i];
    const double peak_uncert = p.amplitudeUncert() > 0.0 ? p.amplitudeUncert() : sqrt(p.amplitude());
    const double nsigma = p.amplitude() / peak_uncert;
    if( nsigma > post_cont_max_nsigma )
    {
      post_cont_max_nsigma = nsigma;
      post_cont_most_sig_peak_index = i;
    }
  }//for( size_t i = 0; i < best_post_cont_type_peaks.size(); ++i )

  // Now we'll check if we should try skew; this isnt perfect, but we'll use just the most significant peak
  //  in the ROI, and check its residual sums.
  if( post_cont_max_nsigma > final_fit_settings.skew_nsigma )
  {
    if( debug_peak )
      cout << "Will check for skew" << endl;

    //Using 1.5 sigma below is arbitrary, but we'll use it for now
    const PeakDef &most_sig_peak = best_post_cont_type_peaks[post_cont_most_sig_peak_index];
    const double peak_lower_energy = most_sig_peak.mean() - 1.5*most_sig_peak.sigma();
    const double left_lower_energy = peak_lower_energy - 1.5*most_sig_peak.sigma();
    const double peak_upper_energy = most_sig_peak.mean() + 1.5*most_sig_peak.sigma();
    const double right_upper_energy = peak_upper_energy + 1.5*most_sig_peak.sigma();

    const size_t left_start_channel = data->find_gamma_channel( left_lower_energy );
    const size_t left_end_channel = data->find_gamma_channel( peak_lower_energy );

    const size_t right_start_channel = data->find_gamma_channel( peak_upper_energy );
    const size_t right_end_channel = data->find_gamma_channel( right_upper_energy );

    const bool too_narrow = (((left_start_channel + 2) >= left_end_channel) || ((right_start_channel + 2) >= right_end_channel));
    //assert( !too_narrow );

    double left_residual_sum = 0.0, right_residual_sum = 0.0;
    if( !too_narrow && (most_sig_peak.skewType() == PeakDef::SkewType::NoSkew) ) //We'll only try skew, if we dont currently have a user-preffered skew type
    {
      const auto sum_residuals = [&]( const size_t start_channel, const size_t end_channel) -> double {
        double sum = 0.0;
        for( size_t i = start_channel; i <= end_channel; ++i )
        {
          const double data_counts = data->gamma_channel_content(i);
          const float lenergy = data->gamma_channel_lower(i);
          const float uenergy = data->gamma_channel_upper(i);
          const double expected_continuum_counts = most_sig_peak.continuum()->offset_integral( lenergy, uenergy, data );
          const double expected_gauss_counts = most_sig_peak.gauss_integral( lenergy, uenergy );
          const double expected_counts = expected_continuum_counts + expected_gauss_counts;
          const double residual = (data_counts - expected_counts) / sqrt( std::max( 1.0, data_counts ) );
          sum += residual;
        }
        return sum / (end_channel - start_channel + 1);
      };//sum_residuals lamda

      left_residual_sum = sum_residuals( left_start_channel, left_end_channel );
      right_residual_sum = sum_residuals( right_start_channel, right_end_channel );

      const bool try_left_skew = left_residual_sum > final_fit_settings.left_residual_sum_min_to_try_skew;
      const bool try_right_skew = right_residual_sum > final_fit_settings.right_residual_sum_min_to_try_skew;

      vector<PeakDef::SkewType> skew_types_to_try;
      if( try_left_skew )
      {
        skew_types_to_try.push_back( PeakDef::SkewType::Bortel );
        skew_types_to_try.push_back( PeakDef::SkewType::GaussExp );
        skew_types_to_try.push_back( PeakDef::SkewType::CrystalBall );
      }//if( try left-sided skew )

      if( try_left_skew && try_right_skew )
      {
        skew_types_to_try.push_back( PeakDef::SkewType::ExpGaussExp );
        skew_types_to_try.push_back( PeakDef::SkewType::DoubleSidedCrystalBall );
      }//if( try both-sided skews )

      // Don't duplicate fitting the same skew type as the most significant peak
      skew_types_to_try.erase( std::remove_if(begin(skew_types_to_try), end(skew_types_to_try), [&most_sig_peak](const PeakDef::SkewType &a){
        return (a == most_sig_peak.skewType());
      }), end(skew_types_to_try) );


      vector<tuple<PeakDef::SkewType,vector<PeakDef>,double>> trial_skew_fits;
      assert( most_sig_peak.skewType() == PeakDef::SkewType::NoSkew );
      trial_skew_fits.emplace_back( most_sig_peak.skewType(), best_post_cont_type_peaks, after_cont_type_chi2dof );

      if( debug_peak )
      {
        cout << "  Checking: {";
        for( const PeakDef::SkewType skew_type : skew_types_to_try )
          cout << PeakDef::to_string(skew_type) << ", ";
        cout << "}" << endl;
      }

      for( const PeakDef::SkewType skew_type : skew_types_to_try )
      {
        vector<PeakDef> these_input_peaks = best_post_cont_type_peaks;
        shared_ptr<PeakContinuum> cont = make_shared<PeakContinuum>( *orig_continuum );
        for( PeakDef &p : these_input_peaks )
        {
          p.setSkewType( skew_type );
          p.setContinuum( cont ); //so we dont mess with the original continuum
        }

        vector<PeakDef> these_fit_peaks;

        try
        {
          these_fit_peaks = refit_peaks( these_input_peaks );

          if( these_fit_peaks.empty() )
            throw runtime_error( "failed fitting" );
        }catch( const std::exception &e )
        {
          if( debug_peak || PeakFitImprove::debug_printout )
            cout << "Error adding skew to peaks: " << e.what()  << " for skew " << PeakDef::to_string(skew_type) << endl;
          continue;
        }//try / catch to fit peaks

        const double these_fit_chi2 = chi2_for_region_wrapper( these_fit_peaks, data, compare_lower_channel, compare_upper_channel );
        trial_skew_fits.emplace_back( skew_type, these_fit_peaks, these_fit_chi2/approx_dof );
      }//for( const PeakDef::SkewType skew_type : skew_types_to_try )

      // Sort the skew fits by the skew type (I'm pretty sure they should already be sorted)
      std::sort( begin(trial_skew_fits), end(trial_skew_fits), []( const auto &a, const auto &b ){
        return std::get<0>(a) < std::get<0>(b);
      } );

      size_t best_skew_index = 0;
      bool found_orig_skew_index = false;
      for( size_t i = 0; i < trial_skew_fits.size(); ++i )
      {
        if( std::get<0>(trial_skew_fits[best_skew_index]) == most_sig_peak.skewType() )
        {
          best_skew_index = i;
          found_orig_skew_index = true;
          break;
        }
      }
      assert( found_orig_skew_index );
      assert( best_skew_index == 0 );

      for( size_t i = 0; i < trial_skew_fits.size(); ++i )
      {
        if( i == best_skew_index )
          continue;

        const double prev_chi2 = std::get<2>(trial_skew_fits[best_skew_index]);
        const double this_chi2 = std::get<2>(trial_skew_fits[i]);

        if( (this_chi2 - final_fit_settings.skew_improve_chi2_dof_threshold) < prev_chi2 )
          best_skew_index = i;
      }//for( size_t i = 1; i < trial_skew_fits.size(); ++i )

      if( (debug_peak || PeakFitImprove::debug_printout) && (trial_skew_fits.size() > 1) )
      {
        double sum_peak_area = 0.0;
        for( const auto &peak : std::get<1>(trial_skew_fits[best_skew_index]) )
          sum_peak_area += peak.peakArea();

        cout << "Post skew fit: Chosen index: " << best_skew_index << " gives Chi2Dof: " << std::get<2>(trial_skew_fits[best_skew_index])
        << " (SkewType: " << PeakDef::to_string( std::get<0>(trial_skew_fits[best_skew_index]) ) << ")"
        << "\n\tcompared to initial Chi2Dof: " << initial_fit_chi2/approx_dof
        << " (SkewType: " << PeakDef::to_string( std::get<0>(trial_skew_fits[0]) ) << ")"
        << "\n\tAll fits: [";

        for( size_t i = 0; i < trial_skew_fits.size(); ++i )
        {
          cout << "\n\t\t{" << PeakDef::to_string( std::get<0>(trial_skew_fits[i]) )
          << ", " << std::get<2>(trial_skew_fits[i]) << "}, ";
        }//for( size_t i = 1; i < trial_peak_fits.size(); ++i )
        cout << "\n\t]\n\tskew_improve_chi2_dof_threshold=" << final_fit_settings.skew_improve_chi2_dof_threshold
        << "\n\tPeaks have area: " << sum_peak_area << " and post_cont_max_nsigma=" << post_cont_max_nsigma;
        cout << endl;
      }
      answer = std::get<1>(trial_skew_fits[best_skew_index]);
    }else
    {
      if( debug_peak )
        cout << "  Skipping skew check: too_narrow=" << too_narrow << ", PrevSkewType=" << PeakDef::to_string(most_sig_peak.skewType())
        << ", left_start_channel=" << left_start_channel << ", left_end_channel=" << left_end_channel
        << ", right_start_channel=" << right_start_channel << ", right_end_channel=" << right_end_channel
        << endl;
    }//if( !too_narrow && (most_sig_peak.skew_type() == PeakDef::SkewType::NoSkew) )
  }else
  {
    if( debug_peak )
      cout << "Not checking for skew." << endl;
  }//if( post_cont_type_max_nsigma > final_fit_settings.skew_nsigma ) / else


  // Now if there are multiple peaks in the ROI, we'll check if we should break it up into multiple ROIs
  //  TODO: It would be nice to do this at an earlier stage (like right after the different types of continuums are fit, and to do this for each type of continuum), but each call to fit peaks assumes a single ROI - so logic would need adjusting for this in all parts that might fit the peaks afterwards
  if( (answer.size() > 1) && (final_fit_settings.break_multi_roi_up_continuum_away_sigma > 1.0) )
  {
    assert( std::count_if(begin(answer), end(answer), [&answer](const PeakDef &p) {
      return (p.continuum() != answer.front().continuum());
    }) == 0 );

    //Start with first peak mean, and add `roi_extent_low_num_fwhm_base_` FWHM, then check for each channel if we
    //  should break them until we get to within `roi_extent_high_num_fwhm_base_` FWHM of second peak, etc.
    //  If we should break up a ROI
    //final_fit_settings.break_multi_roi_up_required_chi2dof_improve
    vector<vector<PeakDef>> splitup_rois;
    const vector<PeakDef> &orig_roi = answer;
    assert( std::is_sorted(begin(orig_roi), end(orig_roi), []( const PeakDef &lhs, const PeakDef &rhs ){
      return lhs.mean() < rhs.mean();
    }) );

    // We will only break the ROI up at 1 point - where the deviation is max
    size_t max_dev_channel = 0, last_left_peak_index = 0;
    double max_deviation = -1.0;
    for( size_t orig_peak_index = 0; (orig_peak_index + 1) < orig_roi.size(); ++orig_peak_index )
    {
      const PeakDef &this_peak = orig_roi[orig_peak_index];
      const PeakDef &next_peak = orig_roi[orig_peak_index + 1];
      const double start_energy = this_peak.mean() + this_peak.fwhm() * lower_extend_base_nfwhm;
      const double end_energy = next_peak.mean() - next_peak.fwhm() * upper_extend_base_nfwhm;
      const size_t start_check_channel = data->find_gamma_channel( static_cast<float>(start_energy) );
      const size_t end_check_channel = data->find_gamma_channel( static_cast<float>(end_energy) );

      //We will use two data channels at a time to compute things... this is arbitrarily chosen and not checked
      for( size_t i = start_check_channel; i < end_check_channel; ++i )
      {
        const float data_cnts = data->gamma_channel_content(i)
        + data->gamma_channel_content(i + 1);
        const double uncert = (data_cnts > 0.1) ? sqrt(data_cnts) : 1.0;
        const float ch_lower_energy = data->gamma_channel_lower(i);
        const float ch_upper_energy = data->gamma_channel_upper(i + 1);

        double fit_amp = 0.0;
        for( const PeakDef &p : orig_roi )
          fit_amp += p.gauss_integral(ch_lower_energy, ch_upper_energy);
        fit_amp += orig_roi.front().continuum()->offset_integral(ch_lower_energy, ch_upper_energy, data);

        const double deviation = fabs(fit_amp - data_cnts) / uncert;
        if( deviation > max_deviation )
        {
          max_deviation = deviation;
          max_dev_channel = i;
          last_left_peak_index = orig_peak_index;
        }
      }//for( size_t i = start_check_channel; i < end_check_channel; ++i )
    }//for( size_t orig_peak_index = 0; orig_peak_index < orig_roi.size(); ++orig_peak_index )


    if( max_deviation > final_fit_settings.break_multi_roi_up_continuum_away_sigma )
    {
      const double break_energy = data->gamma_channel_center( max_dev_channel );
      const double orig_roi_lower = orig_roi.front().continuum()->lowerEnergy();
      const double orig_roi_upper = orig_roi.front().continuum()->upperEnergy();

      assert( break_energy > orig_roi_lower );
      assert( break_energy < orig_roi_upper );

      vector<PeakDef> left_peaks;
      shared_ptr<PeakContinuum> left_cont = make_shared<PeakContinuum>( *orig_roi.front().continuum() );
      left_cont->setRange( orig_roi_lower, break_energy );
      left_cont->setParameters( orig_roi_lower, left_cont->parameters(), left_cont->uncertainties() );

      for( size_t i = 0; i <= last_left_peak_index; ++i )
      {
        PeakDef p = orig_roi[i];
        p.setContinuum( left_cont );
        left_peaks.push_back( p );
      }


      vector<PeakDef> right_peaks;
      shared_ptr<PeakContinuum> right_cont = make_shared<PeakContinuum>( *orig_roi.front().continuum() );
      right_cont->setRange( break_energy, orig_roi_upper );
      right_cont->setParameters( break_energy, right_cont->parameters(), right_cont->uncertainties() );

      for( size_t i = last_left_peak_index + 1; i < orig_roi.size(); ++i )
      {
        PeakDef p = orig_roi[i];
        p.setContinuum( right_cont );
        right_peaks.push_back( p );
      }

      assert( !left_peaks.empty() );
      assert( !right_peaks.empty() );


      vector<PeakDef> left_fit_peaks, right_fit_peaks;
      try
      {
        left_fit_peaks = refit_peaks( left_peaks );
        right_fit_peaks = refit_peaks( right_peaks );

        if( left_fit_peaks.empty() || (left_fit_peaks.size() != left_peaks.size()) )
          throw runtime_error( "Diff number of left peaks - prev: " + to_string(left_peaks.size()) + ", fit: " + to_string(left_fit_peaks.size()) );

        if( right_fit_peaks.empty() || (right_fit_peaks.size() != right_peaks.size()) )
          throw runtime_error( "Diff number of right peaks - prev: " + to_string(right_peaks.size()) + ", fit: " + to_string(right_fit_peaks.size()) );

        const double lower_frac = (break_energy - orig_roi_lower) / (orig_roi_upper - orig_roi_lower);
        const double upper_frac = (orig_roi_upper - break_energy) / (orig_roi_upper - orig_roi_lower);
        assert( fabs((lower_frac + upper_frac) - 1.0) < 1.0E-6 );

        const double left_chi2dof = left_fit_peaks.front().chi2dof();
        const double right_chi2dof = right_fit_peaks.front().chi2dof();
        const double orig_dof = orig_roi.front().chi2dof();
        const double new_weighted_dof = lower_frac*left_chi2dof + upper_frac*right_chi2dof;
        if( (new_weighted_dof + final_fit_settings.break_multi_roi_up_required_chi2dof_improve) < orig_dof )
        {
          answer.clear();
          answer.insert( end(answer), begin(left_fit_peaks), end(left_fit_peaks) );
          answer.insert( end(answer), begin(right_fit_peaks), end(right_fit_peaks) );

          if( debug_peak || PeakFitImprove::debug_printout )
          {
            cout << "Split ROI UP!  Chi2Dof " << orig_dof << " --> " << left_chi2dof << " + " << right_chi2dof
            << " (weighted: " << new_weighted_dof << ") with ROIs:"
            << " [[" << orig_roi_lower << ", " << break_energy << "], [" << break_energy << ", " << orig_roi_upper << "]]"
            << " - on file title='" << data->title() << "'" << endl;
          }
        }
      }catch( const std::exception &e )
      {
        // TODO/Note: It doesn look like we can get here, neither `PeakFitLM::fit_peaks_LM(...)` or `fitPeaks(...)` will through
        if( debug_peak || PeakFitImprove::debug_printout )
          cout << "Error fitting broken-up ROI peaks: " << e.what() << endl;
      }//try / catch to fit peaks
    }//if( max_deviation > final_fit_settings.break_multi_roi_up_continuum_away_sigma )
  }//if( answer.size() > 1 )


  return answer;
}//vector<PeakDef> final_peak_fit_for_roi(...)


vector<PeakDef> final_peak_fit( const vector<PeakDef> &pre_fit_peaks,
                               const FinalPeakFitSettings &final_fit_settings,
                               const bool isHPGe,
                               const std::shared_ptr<const SpecUtils::Measurement> &data,
                               const bool multithread )
{
  // This function deals with `FinalPeakFitSettings::combine_nsigma_near` and
  //  `FinalPeakFitSettings::combine_ROI_overlap_frac`, and then passes the peaks
  //  off to `final_peak_fit_for_roi(...)`
  vector<PeakDef> answer;

  // Create a vector of ROIs, sorted by the ROI middle energy.
  vector<pair<shared_ptr<const PeakContinuum>,vector<PeakDef>>> rois;
  for( const PeakDef &p : pre_fit_peaks )
  {
    const auto comp_eq = [&p]( const pair<shared_ptr<const PeakContinuum>,vector<PeakDef>> &other ){
      return other.first == p.continuum();
    };
    const auto comp_lt = []( const pair<shared_ptr<const PeakContinuum>,vector<PeakDef>> &lhs,
                       const pair<shared_ptr<const PeakContinuum>,vector<PeakDef>> &rhs ){
      return (lhs.first->lowerEnergy() + lhs.first->upperEnergy()) < (rhs.first->lowerEnergy() + rhs.first->upperEnergy());
    };

    vector<pair<shared_ptr<const PeakContinuum>,vector<PeakDef>>>::iterator pos
                                           = std::find_if( begin(rois), end(rois), comp_eq );

    if( pos == end(rois) )
    {
      pair<shared_ptr<const PeakContinuum>,vector<PeakDef>> val{ p.continuum(), {p} };
      const auto lb_pos = std::lower_bound( begin(rois), end(rois), val, comp_lt );
      pos = rois.insert( lb_pos, std::move(val) );
    }else
    {
      const auto lb_pos = std::lower_bound( begin(pos->second), end(pos->second), p, &PeakDef::lessThanByMean );
      pos->second.insert( lb_pos, p );
      assert( std::is_sorted( begin(pos->second), end(pos->second), &PeakDef::lessThanByMean ) );
    }
  }//for( const PeakDef &p : pre_fit_peaks )


  bool combined_any_rois = false;
  for( size_t roi_index = 1; roi_index < rois.size(); ++roi_index )
  {
    const pair<shared_ptr<const PeakContinuum>,vector<PeakDef>> &prev_roi = rois[roi_index-1];
    const shared_ptr<const PeakContinuum> &prev_cont = prev_roi.first;

    const pair<shared_ptr<const PeakContinuum>,vector<PeakDef>> &this_roi = rois[roi_index];
    const shared_ptr<const PeakContinuum> &this_cont = this_roi.first;
    const vector<PeakDef> &this_peaks = this_roi.second;

    //bool combine_for_overlap = false;
    bool combine_rois = false;
    const double prev_lower = prev_cont->lowerEnergy();
    const double prev_upper = prev_cont->upperEnergy();

    const double this_lower = this_cont->lowerEnergy();
    const double this_upper = this_cont->upperEnergy();

    assert( (prev_lower + prev_upper) <= (this_lower + this_upper) );
    const double overlap_start = std::max(prev_lower, this_lower); //unnessary step since things are already sorted
    const double overlap_end = std::min(prev_upper, this_upper);
    if( overlap_start <= overlap_end )
    {
      const double overlap = overlap_end - overlap_start;
      assert( overlap >= 0.0 );
      const double lesser_roi_extent = std::min( (prev_upper - prev_lower), (this_upper - this_lower) );
      //combine_for_overlap = ((overlap / lesser_roi_extent) > final_fit_settings.combine_ROI_overlap_frac);
    }//if( overlap_start <= overlap_end )

    const PeakDef &last_prev_peak = prev_roi.second.back();
    const PeakDef &first_this_peak = this_roi.second.front();
    const double avrg_fwhm = 0.5*(last_prev_peak.fwhm() + first_this_peak.fwhm());
    const double energy_diff = fabs( last_prev_peak.mean() - first_this_peak.mean() );
    const double num_fwhm_near = energy_diff / avrg_fwhm;
    combine_rois = (num_fwhm_near < final_fit_settings.require_combine_num_fwhm_near);

    if( !combine_rois && (num_fwhm_near < final_fit_settings.not_allow_combine_num_fwhm_near) && (overlap_start <= overlap_end) )
    {
#warning "blah blah blah - do a better job trying to combine"
      const double overlap = overlap_end - overlap_start;
      assert( overlap >= 0.0 );
      const double lesser_roi_extent = std::min( (prev_upper - prev_lower), (this_upper - this_lower) );
      if( (overlap / lesser_roi_extent) > 0.15 )
      {
        //cout << "Combining [" << prev_lower << ", " << prev_upper << "] and [" << this_lower << ", " << this_upper << "]" << endl;
        combine_rois = true;
      }
    }//if( !combine_rois && (num_fwhm_near < final_fit_settings.not_allow_combine_num_fwhm_near) )


    if( combine_rois )
    {
      combined_any_rois = true;

      //Combine ROIs into a single ROI
      shared_ptr<PeakContinuum> new_cont = make_shared<PeakContinuum>( *prev_cont );
      const double new_lower = std::min(prev_lower, this_lower); //just to make sure nothings wonky
      const double new_upper = std::max(prev_upper, this_upper); //just to make sure nothings wonky
      new_cont->setRange( new_lower, std::max(prev_upper, new_upper) );
      const PeakContinuum::OffsetType new_offset_type = std::max( prev_cont->type(), this_cont->type() );
      if( new_offset_type != new_cont->type() )
      {
        const double reference_energy = 0.5*(new_lower + new_upper);
        new_cont->calc_linear_continuum_eqn( data, reference_energy, new_lower, new_upper, 3, 3 );
        new_cont->setType( new_offset_type );
      }//if( new_offset_type != new_cont->type() )

      // Update `rois` to use new continuum
      rois[roi_index-1].first = new_cont;

      // Add peaks from
      vector<PeakDef> &prev_peaks = rois[roi_index-1].second;
      prev_peaks.insert( end(prev_peaks), begin(this_peaks), end(this_peaks) );
      for( PeakDef &p : prev_peaks )
        p.setContinuum( new_cont );
      std::sort( begin(prev_peaks), end(prev_peaks), &PeakDef::lessThanByMean );

      // Erase current peaks, now that we have added them to the previous ROI
      rois.erase( begin(rois) + roi_index );

      assert( roi_index > 0 );

      // Lets check to make sure the peaks arent some-how duplicates - if they are, we'll erase the peaks that are
      //  about the same energy.
      for( int index = 0; (index + 1) < static_cast<int>(prev_peaks.size()); ++index )
      {
        const double avrg_sigma = 0.5*(prev_peaks[index].sigma() + prev_peaks[index+1].sigma());
        const double mean_diff = fabs(prev_peaks[index].mean() - prev_peaks[index+1].mean());

        const double sigma_near_threshold = 0.75;  // Arbitrarily chosen - could probably even be made larger to like 1.0
        if( mean_diff < sigma_near_threshold*avrg_sigma )
        {
          prev_peaks.erase( begin(prev_peaks) + index + 1 );
          index -= 1;
        }
      }//for( check if ROIs are duplicates )

      // Decrement the loop index, since now the next ROI is one element closer
      roi_index -= 1;
    }//if( combine_rois )
  }//for( size_t size_t roi_index = 1; roi_index < rois.size(); ++i )


  if( combined_any_rois && PeakFitImprove::debug_printout )
  {
    cout << "rois after combining:" << endl;
    for( const auto &p : rois )
    {
      cout << "    [" << p.first->lowerEnergy() << ", " << p.first->upperEnergy() << "]: ";
      for( const auto &peak : p.second )
        cout << "{" << peak.mean() << ", " << peak.fwhm() << "}, ";
      cout << endl;
    }
    cout << endl;
  }//if( debug_peak || PeakFitImprove::debug_printout )

  if( multithread )
  {
    std::mutex answer_mutex;

    SpecUtilsAsync::ThreadPool pool;

    for( const auto &cont_peaks : rois )
    {
      const vector<PeakDef> *in_peaks = &cont_peaks.second;

      pool.post( [&answer, &answer_mutex, in_peaks, &final_fit_settings, &data, isHPGe](){
        try
        {
          const vector<PeakDef> roi_answer = final_peak_fit_for_roi( *in_peaks, final_fit_settings, isHPGe, data );

          std::lock_guard<std::mutex> lock( answer_mutex );
          answer.insert( end(answer), begin(roi_answer), end(roi_answer) );
        }catch( std::exception &e )
        {
          if( PeakFitImprove::debug_printout )
          {
            cerr << "Fit of ROI=[" << (*in_peaks)[0].continuum()->lowerEnergy() << ", "
            << (*in_peaks)[0].continuum()->upperEnergy() << "] failed - will use input peaks as output.  Err: " << e.what()
            << endl;
          }

          std::lock_guard<std::mutex> lock( answer_mutex );
          answer.insert( end(answer), begin(*in_peaks), end(*in_peaks) );
        }
      } );
    }//for( const auto &cont_peaks : rois )
  }else
  {
    for( const auto &cont_peaks : rois )
    {
      const vector<PeakDef> &in_peaks = cont_peaks.second;

      try
      {
        const vector<PeakDef> roi_answer = final_peak_fit_for_roi( in_peaks, final_fit_settings, isHPGe, data );
        answer.insert( end(answer), begin(roi_answer), end(roi_answer) );
      }catch( std::exception &e )
      {
        if( PeakFitImprove::debug_printout )
        {
          cerr << "Fit of ROI=[" << in_peaks[0].continuum()->lowerEnergy() << ", "
          << in_peaks[0].continuum()->upperEnergy() << "] failed - will use input peaks as output.  Err: " << e.what()
          << endl;
        }
        answer.insert( end(answer), begin(in_peaks), end(in_peaks) );
      }
    }
  }//if( final_peak_fit_parrallel ) / else

  std::sort( begin(answer), end(answer), &PeakDef::lessThanByMean );

  return answer;
}//vector<PeakDef> final_peak_fit

FinalFitScore calculate_final_fit_score(
  const std::vector<PeakDef> &fit_peaks,
  const std::vector<ExpectedPhotopeakInfo> &expected_photopeaks,
  const double num_sigma_contribution )
{
  using namespace std;

  FinalFitScore score;

  for( const PeakDef &found_peak : fit_peaks )
  {
    const double found_energy = found_peak.mean();
    const double found_fwhm = found_peak.fwhm();
    const double peak_lower_contrib = found_energy - num_sigma_contribution*found_peak.sigma();
    const double peak_upper_contrib = found_energy + num_sigma_contribution*found_peak.sigma();

    //Note that ExpectedPhotopeakInfo cooresponds to a grouping of gamma lines we expect to
    // detect as a single peak.
    const ExpectedPhotopeakInfo * nearest_expected_peak = nullptr;

    for( const ExpectedPhotopeakInfo &expected_peak : expected_photopeaks )
    {
      const double expected_energy = expected_peak.effective_energy;

      if( (expected_energy >= peak_lower_contrib) && (expected_energy <= peak_upper_contrib) )
      {
        if( !nearest_expected_peak || (fabs(expected_energy - found_energy) < fabs(nearest_expected_peak->effective_energy - found_energy)) )
          nearest_expected_peak = &expected_peak;
      }
    }//for( const ExpectedPhotopeakInfo &expected_peak : expected_photopeaks )

    if( nearest_expected_peak )
    {
      score.num_peaks_used += 1;

      const double expected_energy = nearest_expected_peak->effective_energy;
      const double expected_fwhm = nearest_expected_peak->effective_fwhm;
      const double expected_sigma = expected_fwhm/2.35482;
      const double expected_area = nearest_expected_peak->peak_area;

      const double area_diff = fabs(found_peak.amplitude() - expected_area);
      const double area_score = area_diff / sqrt( (expected_area < 1.0) ? 1.0 : expected_area );
      const double width_score = fabs( expected_fwhm - found_fwhm ) / expected_sigma;
      const double position_score = fabs(found_energy - expected_energy) / expected_sigma;

      score.area_score += std::min( area_score, 20.0 );
      score.width_score += std::min( width_score, 1.0 );
      score.position_score += std::min( position_score, 1.5 );
    }else
    {
      // Found an extra peak we didn't expect
      if( found_peak.amplitude() < 1.0 ) //Not a real peak, so ignore it
        continue;

      const double sqrt_area = sqrt( found_peak.amplitude() );
      const double area_uncert = found_peak.amplitudeUncert() > 0.0 ? (std::max)(found_peak.amplitudeUncert(), sqrt_area) : sqrt_area;

      // If this is a significant peak that we didn't expect, count it
      // (Note: In eval_final_peak_fit, escape peaks are excluded, but here we don't have that info)
      score.ignored_unexpected_peaks += 1;
      score.unexpected_peaks_sum_significance += std::min( 7.5, found_peak.amplitude() / area_uncert ); //Cap at 7.5 sigma (arbitrary)
    }//if( nearest_expected_peak ) / else
  }//for( const PeakDef &found_peak : fit_peaks )

  // Calculate total weight for scoring
  if( score.num_peaks_used <= 1 )
    score.total_weight = score.area_score;
  else
    score.total_weight = score.area_score / score.num_peaks_used;

  return score;
}//calculate_final_fit_score(...)


FinalFitScore eval_final_peak_fit( const FinalPeakFitSettings &final_fit_settings,
                           const DataSrcInfo &src_info,
                           const vector<PeakDef> &intial_peaks,
                           const bool write_n42 )
{
  FinalFitScore score;

  const vector<shared_ptr<const SpecUtils::Measurement>> &src_spectra = src_info.src_info.src_spectra;
  assert( !src_spectra.empty() );
  const shared_ptr<const SpecUtils::Measurement> &data = src_spectra.front();
  assert( data );

#warning "final_peak_fit_for_roi: always assuming HPGe right now - for dev"
  const bool isHPGe = true;

  vector<PeakDef> fit_peaks = final_peak_fit( intial_peaks, final_fit_settings, isHPGe, data, ns_final_peak_fit_parrallel );

  // - Calculate score based on fit_peaks, we will judge on:
  //   - Area of each peaks
  //   - Peak width
  //   - Peak position

  // Single and double escape fraction for 20% generic HPGe
  const auto single_escape_sf = []( const double x ) -> double {
    return std::max( 0.0, (1.8768E-11 *x*x*x) - (9.1467E-08 *x*x) + (2.1565E-04 *x) - 0.16367 );
  };

  const auto double_escape_sf = []( const double x ) -> double {
    return std::max( 0.0, (1.8575E-11 *x*x*x) - (9.0329E-08 *x*x) + (2.1302E-04 *x) - 0.16176 );
  };

  vector<tuple<double,double,double,double>> possible_escape_peak_parents;
  for( const auto &p : src_info.expected_photopeaks )
  {
    if( p.effective_energy > 1060 )
    {
      const double se_area = p.peak_area * single_escape_sf( p.effective_energy );
      const double de_area = p.peak_area * double_escape_sf( p.effective_energy );
      possible_escape_peak_parents.emplace_back( p.effective_energy, p.peak_area, se_area, de_area );
    }//if( p.effective_energy > 1060 )
  }//for( const auto &p : src_info.expected_photopeaks )

  // Sort by S.E. amplitude
  std::sort( begin(possible_escape_peak_parents), end(possible_escape_peak_parents),
    []( const tuple<double,double,double,double> &lhs, const tuple<double,double,double,double> &rhs ){
      return get<2>(lhs) > get<2>(rhs);
  });

  // Pick only top 10 (arbitrary) energies to create single escape peaks
  if( possible_escape_peak_parents.size() > 10 )
    possible_escape_peak_parents.resize( 10 );


  // The number of sigma away from expected, that we will use to calculate the peak area.
  const double num_sigma_contribution = 1.5;

  for( PeakDef &found_peak : fit_peaks )
  {
    const double found_energy = found_peak.mean();
    const double found_fwhm = found_peak.fwhm();
    const double peak_lower_contrib = found_energy - num_sigma_contribution*found_peak.sigma();
    const double peak_upper_contrib = found_energy + num_sigma_contribution*found_peak.sigma();

    //Note that ExpectedPhotopeakInfo cooresponds to a grouping of gamma lines we expect to
    // detect as a single peak.
    const ExpectedPhotopeakInfo * nearest_expected_peak = nullptr;

    for( const ExpectedPhotopeakInfo &expected_peak : src_info.expected_photopeaks )
    {
      const double expected_energy = expected_peak.effective_energy;

      if( (expected_energy >= peak_lower_contrib) && (expected_energy <= peak_upper_contrib) )
      {
        if( !nearest_expected_peak || (fabs(expected_energy - found_energy) < fabs(nearest_expected_peak->effective_energy - found_energy)) )
          nearest_expected_peak = &expected_peak;
      }
    }//for( const ExpectedPhotopeakInfo &expected_peak : src_info.expected_photopeaks )

    if( nearest_expected_peak )
    {
      score.num_peaks_used += 1;

      const double expected_energy = nearest_expected_peak->effective_energy;
      const double expected_fwhm = nearest_expected_peak->effective_fwhm;
      const double expected_sigma = expected_fwhm/2.35482;
      const double expected_area = nearest_expected_peak->peak_area;

      const double area_diff = fabs(found_peak.amplitude() - expected_area);
      const double area_score = area_diff / sqrt( (expected_area < 1.0) ? 1.0 : expected_area );
      const double width_score = fabs( expected_fwhm - found_fwhm ) / expected_sigma;
      const double position_score = fabs(found_energy - expected_energy) / expected_sigma;

      score.area_score += std::min( area_score, 20.0 );
      score.width_score += std::min( width_score, 1.0 );
      score.position_score += std::min( position_score, 1.5 );

      if( write_n42 )
      {
        const double nsgima_off = (found_peak.amplitude() - expected_area) / found_peak.amplitudeUncert();
        //string label = std::format("TrueArea: {:.0f}, σ_off: {:.1f}", expected_area, nsgima_off );
        char label[128] = { '\0' };
        snprintf( label, sizeof(label), "TrueArea: %.0f, σ_off: %.1f", expected_area, nsgima_off );
        found_peak.setUserLabel( label );
        found_peak.setLineColor( Wt::GlobalColor::darkGreen );
      }//if( write_n42 )
    }else
    {
      // We found an extra peak we dont want - lets check if its a single or double escape peak
      //  and if so, ignore it; otherwise I guess we have to punish for this peak being here.
      bool is_escape_peak = false;
      for( const auto &ep : possible_escape_peak_parents )
      {
        if( (fabs( (get<0>(ep)-511.0) - found_energy) < 0.2*found_fwhm)
            ||  (fabs( (get<0>(ep)-1022.0) - found_energy) < 0.2*found_fwhm) )
        {
          is_escape_peak = true;
          break;
        }
      }

      // If it is an escape peak, ignore it.
      if( is_escape_peak )
      {
        found_peak.setUserLabel( "Escape Peak" );
        continue;
      }

      // If this is a significant peak, we messed up, so lets ignore it...
      const double area_uncert = found_peak.amplitudeUncert() > 0.0 ? found_peak.amplitudeUncert() : 1/sqrt(found_peak.amplitude());
      if( (found_peak.amplitude() > 8.0*area_uncert) || (fabs(found_peak.mean() - 511.0) < 5.0) )
      {
        if( PeakFitImprove::debug_printout )
          cout << "Found peak {mean: " << found_peak.mean()
          << ", fwhm: " << found_peak.fwhm()
          << ", area: " << found_peak.amplitude()
          << ", area_uncert: " << found_peak.amplitudeUncert() << "} that is unexpected, but too large, so ignoring." << endl;

        found_peak.setUserLabel( "Unexpected Large Peak" );
        found_peak.setLineColor( Wt::GlobalColor::magenta );
        continue;
      }

      if( found_peak.amplitude() < 1.0 ) //Not a real peak, so ignore it.
        continue;

      score.ignored_unexpected_peaks += 1;
      score.unexpected_peaks_sum_significance += std::min( 7.5, found_peak.amplitude() / area_uncert ); //Cap at 7.5 sigma (arbitrary)

      if( write_n42 )
      {
        found_peak.setUserLabel( "Unexpected Peak" );
        found_peak.setLineColor( Wt::GlobalColor::darkRed );
      }//if( write_n42 )
    }//if( !nearest_expected_peak )
  }//for( PeakDef &found_peak : fit_peaks )

  // TODO: we could/should add in weight for peak area_score, width_score, and position_score
  if( score.num_peaks_used <= 1 )
    score.total_weight = score.area_score;
  else
    score.total_weight = score.area_score / score.num_peaks_used;

  if( write_n42 )
  {
    string outdir = "output_n42";
    if( !SpecUtils::is_directory(outdir) && !SpecUtils::create_directory(outdir) )
      cerr << "Failed to create directory '" << outdir << "'" << endl;

    outdir = SpecUtils::append_path( outdir, src_info.detector_name );
    if( !SpecUtils::is_directory(outdir) && SpecUtils::create_directory(outdir) )
      cerr << "Failed to create directory '" << outdir << "'" << endl;

    outdir = SpecUtils::append_path( outdir, src_info.location_name );
    if( !SpecUtils::is_directory(outdir) && !SpecUtils::create_directory(outdir) )
      cerr << "Failed to create directory '" << outdir << "'" << endl;

    outdir = SpecUtils::append_path( outdir, src_info.live_time_name );
    if( !SpecUtils::is_directory(outdir) && !SpecUtils::create_directory(outdir) )
      cerr << "Failed to create directory '" << outdir << "'" << endl;

    const string out_n42 = SpecUtils::append_path( outdir, src_info.src_info.src_name ) + "_final_fit.n42";

    SpecMeas output;

    output.add_remark( score.print( "score") );

    output.set_instrument_model( src_info.detector_name );
    if( SpecUtils::icontains(src_info.detector_name, "Detective" ) )
      output.set_manufacturer( "ORTEC" );
    else if( SpecUtils::icontains(src_info.detector_name, "Falcon 5000" ) )
      output.set_manufacturer( "Canberra" );
    else if( SpecUtils::icontains(src_info.detector_name, "Fulcrum" ) )
      output.set_manufacturer( "PHDS" );

    output.set_measurement_location_name( src_info.location_name );

    shared_ptr<SpecUtils::Measurement> out_with_cand = make_shared<SpecUtils::Measurement>( *src_info.src_info.src_spectra.front() );
    out_with_cand->set_sample_number( 1 );

    const string title = src_info.detector_name + "/" + src_info.src_info.src_name;
    out_with_cand->set_title( title + " + final fit peaks" );
    auto now = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());
    out_with_cand->set_start_time( now );
    output.add_measurement( out_with_cand, false );

    deque<shared_ptr<const PeakDef>> peaks;
    for( const PeakDef &p : fit_peaks )
      peaks.push_back( make_shared<PeakDef>(p) );

    output.setPeaks( peaks, {1} );

    //ofstream outstrm( out_n42.c_str(), ios::out | ios::binary );
    //if( !outstrm )
    //  cerr << "Failed to open '" << out_n42 << "'" << endl;
    //output.write( outstrm, output.sample_numbers(), output.detector_names(), SpecUtils::SaveSpectrumAsType::N42_2012 );

    output.save2012N42File( out_n42, [=](){
      cerr << "Failed to write '" << out_n42 << "'" << endl;
    });

    //cout << "Wrote '" << out_n42 << endl;
  }//if( write_n42 )

  return score;
}//double eval_final_peak_fit(...)


FinalPeakFitSettings minuit_fit_final_pars( std::function<double( const FinalPeakFitSettings &)> ga_eval_fcn )
{
  /* Run on 20250511, using just Detective-X (14 hours on M4 - population of 100, and max of 1000 generations)
   best_final_fit_settings.combine_nsigma_near = 13.288225;
   //best_final_fit_settings.combine_ROI_overlap_frac = 0.802212;
   best_final_fit_settings.cont_type_peak_nsigma_threshold = 45.340986;
   best_final_fit_settings.cont_type_left_right_nsigma = 9.624457;
   best_final_fit_settings.cont_poly_order_increase_chi2dof_required = 1.157370;
   best_final_fit_settings.cont_step_type_increase_chi2dof_required = 0.549378;
   best_final_fit_settings.skew_nsigma = 5.420989;
   best_final_fit_settings.left_residual_sum_min_to_try_skew = 1.482605;
   best_final_fit_settings.right_residual_sum_min_to_try_skew = 4.213049;
   best_final_fit_settings.skew_improve_chi2_dof_threshold = 3.138754;
   best_final_fit_settings.roi_extent_low_num_fwhm_base = 3.519398;
   best_final_fit_settings.roi_extent_high_num_fwhm_base = 8.373682;
   best_final_fit_settings.roi_extent_mult_type = RoiExtentMultType::Linear;
   best_final_fit_settings.roi_extent_lower_side_stat_multiple = 0.266548;
   best_final_fit_settings.roi_extent_upper_side_stat_multiple = 0.724294;
   best_final_fit_settings.multi_roi_extent_lower_side_fwhm_mult = 0.613304;
   best_final_fit_settings.multi_roi_extent_upper_side_fwhm_mult = -0.680065;

   */

  FitFinalPeakFitSettingsChi2 chi2Fcn( ga_eval_fcn );

  ROOT::Minuit2::MnUserParameters inputPrams;
  inputPrams.Add( "require_combine_num_fwhm_near", 2.0, 0.5, 1.0, 8.5 );
  inputPrams.Add( "not_allow_combine_num_fwhm_near", 4.0, 0.5, 2.0, 15.0 );
  //inputPrams.Add( "combine_ROI_overlap_frac", 0.802212, 0.25, -1.0, 1.0 );
  inputPrams.Add( "cont_type_peak_nsigma_threshold", 45.340986, 5, 10.0, 100 );
  inputPrams.Add( "cont_type_left_right_nsigma", 9.624457, 3, 1.0, 40.0 );
  inputPrams.Add( "cont_poly_order_increase_chi2dof_required", 1.157370, 0.2, 0.0, 4.0 );
  inputPrams.Add( "cont_step_type_increase_chi2dof_required", 0.549378, 0.2, 0.0, 4.0 );
  inputPrams.Add( "skew_nsigma", 5.420989, 5, 0.0, 25.0 );
  inputPrams.Add( "left_residual_sum_min_to_try_skew", 1.482605, 0.5, 0.0, 10.0 );
  inputPrams.Add( "right_residual_sum_min_to_try_skew", 4.213049, 2.0, 0.0, 10.0 );
  inputPrams.Add( "skew_improve_chi2_dof_threshold", 3.138754, 0.5, 0.0, 5.0 );
  inputPrams.Add( "roi_extent_low_num_fwhm_base", 2.5, 0.5, 0.5, 9.0 );
  inputPrams.Add( "roi_extent_high_num_fwhm_base", 3.5, 0.5, 0.5, 9.0 );
  inputPrams.Add( "roi_extent_low_num_fwhm_extra", 1.5, 0.5, 0.0, 6.0 );
  inputPrams.Add( "roi_extent_high_num_fwhm_extra", 1.5, 0.5, 0.0, 6.0 );
  inputPrams.Add( "roi_end_second_deriv_thresh", 3, 0.0, 0, 20 );
  inputPrams.Add( "break_multi_roi_up_continuum_away_sigma", 4, 1.0, 0.0, 20.0 );
  inputPrams.Add( "break_multi_roi_up_required_chi2dof_improve", 0.5, 0.1, 0.1, 2.0 );

  cerr << "Returning ititial paramaters..." << endl;
  return FitFinalPeakFitSettingsChi2::params_to_settings( inputPrams.Params() );


  ROOT::Minuit2::MnUserParameterState inputParamState( inputPrams );
  ROOT::Minuit2::MnStrategy strategy( 2 ); //0 low, 1 medium, >=2 high

  const unsigned int maxFcnCall = 2000;
  const double tolerance = 0.001;
  ROOT::Minuit2::CombinedMinimizer fitter;
  ROOT::Minuit2::FunctionMinimum minimum
  = fitter.Minimize( chi2Fcn, inputParamState,
                    strategy, maxFcnCall, tolerance );

  //Not sure why Minuit2 doesnt like converging on the minumum verry well, but
  //  rather than showing the user an error message, we'll give it anither try
  if( minimum.IsAboveMaxEdm() )
  {
    ROOT::Minuit2::MnMigrad fitter( chi2Fcn, inputParamState, strategy );
    minimum = fitter( maxFcnCall, tolerance );
  }//if( minimum.IsAboveMaxEdm() )

  if( !minimum.IsValid() )
    throw runtime_error( Wt::WString::tr("dcw-err-failed-fit-AD").toUTF8() );

  const ROOT::Minuit2::MnUserParameters params = minimum.UserState().Parameters();
  const vector<double> pars = params.Params();
  cerr << "Fit " << pars[0] << " g/cm2 with EDM " << minimum.Edm() << endl;


  FinalPeakFitSettings final_settings = FitFinalPeakFitSettingsChi2::params_to_settings( pars );

  return final_settings;
}//FinalPeakFitSettings minuit_fit_final_pars( std::function<double( const FinalPeakFitSettings &)> ga_eval_fcn )


void do_final_peak_fit_ga_optimization( const FindCandidateSettings &candidate_settings,
                           const InitialPeakFindSettings &initial_fit_settings,
                           const vector<DataSrcInfo> &input_srcs )
{
  vector<vector<PeakDef>> initial_peak_fits( input_srcs.size(), vector<PeakDef>{} );

  {// Begin get initial peak fit
    cout << "Starting to fit final peaks." << endl;
    const double start_wall = SpecUtils::get_wall_time();
    const double start_cpu = SpecUtils::get_cpu_time();

    boost::asio::thread_pool pool(PeakFitImprove::sm_num_optimization_threads);
    std::atomic<long long> num_completed{0};
    const size_t num_total = input_srcs.size();

    for( size_t input_src_index = 0; input_src_index < num_total; ++input_src_index )
    {
      const DataSrcInfo &src = input_srcs[input_src_index];
      vector<PeakDef> *result = &(initial_peak_fits[input_src_index]);
      shared_ptr<const SpecUtils::Measurement> data = src.src_info.src_spectra[0];

      auto fit_initial_peaks_worker = [&candidate_settings, &initial_fit_settings, &num_completed, num_total, data, result](){
        size_t dummy1, dummy2;
        *result = InitialFit_GA::initial_peak_find_and_fit( initial_fit_settings, candidate_settings, data, false, dummy1, dummy2 );

        long long val = num_completed.fetch_add(1, std::memory_order_relaxed);
        if( ((val+1) % 100) == 0 )
          cout << "Completed " << (val + 1) << " spectra of " << num_total << " for initial peak fits." << endl;
      }; //fit_initial_peaks_worker

      boost::asio::post(pool, fit_initial_peaks_worker);
    }//for( loop over input_srcs )

    cout << "Submitted initial peak fits to the queue" << endl;
    pool.join();
    cout << "Done fitting initial peaks" << endl;

    const double end_wall = SpecUtils::get_wall_time();
    const double end_cpu = SpecUtils::get_cpu_time();
    cout << "Finished fitting initial peaks - took"
    << " {wall: " << (end_wall - start_wall)
    << " s, cpu: " << (end_cpu - start_cpu) << " s}."
    << endl;
  }// end get initial peak fit

  cout << "Will use candidate and initial fit settings:" << endl;
  cout << candidate_settings.print( "\tcandidate_settings" ) << endl;
  cout << initial_fit_settings.print( "\tinitial_fit_settings" ) << endl;
  cout << "\n" << endl;

  // Now go through and setup genetic algorithm, as well implement final peak fitting code
  const double start_wall = SpecUtils::get_wall_time();
  const double start_cpu = SpecUtils::get_cpu_time();

  const vector<vector<PeakDef>> * const intial_peaks_ptr = &initial_peak_fits;
  const vector<DataSrcInfo> * const data_srcs_ptr = &input_srcs;

  std::function<double( const FinalPeakFitSettings &)> ga_eval_fcn
        = [&intial_peaks_ptr, &data_srcs_ptr]( const FinalPeakFitSettings &settings ) -> double {

    assert( intial_peaks_ptr && data_srcs_ptr );
    const vector<vector<PeakDef>> &intial_peaks_ref = *intial_peaks_ptr;
    const vector<DataSrcInfo> &data_srcs_ref = *data_srcs_ptr;

    assert( intial_peaks_ref.size() == data_srcs_ref.size() );

    FinalFitScore sum_score;

    static std::mutex sm_score_mutex;
    static FinalFitScore sm_best_score{ 0.0, 0.0, 0.0, 0u, 0.0, DBL_MAX, 0u };

#warning "final_peak_fit_for_roi: always assuming HPGe right now - for dev"
  const bool isHPGe = true;

    if( PeakFitImprove::sm_num_threads_per_individual > 1 )
    {
      boost::asio::thread_pool pool(PeakFitImprove::sm_num_threads_per_individual);
      sm_best_score.total_weight = DBL_MAX;

      for( size_t src_index = 0; src_index < data_srcs_ref.size(); ++src_index )
      {
        const DataSrcInfo * const info = &(data_srcs_ref[src_index]);
        const vector<PeakDef> * const initial_peaks = &(intial_peaks_ref[src_index]);

        boost::asio::post(pool, [src_index,info,initial_peaks,settings,&sum_score,isHPGe](){

          const FinalFitScore score = FinalFit_GA::eval_final_peak_fit( settings, *info, *initial_peaks, false );

          std::lock_guard<std::mutex> lock( sm_score_mutex );
          sum_score.area_score += score.area_score;
          sum_score.width_score += score.width_score;
          sum_score.position_score += score.position_score;
          sum_score.ignored_unexpected_peaks += score.ignored_unexpected_peaks;
          sum_score.unexpected_peaks_sum_significance += score.unexpected_peaks_sum_significance;
          sum_score.total_weight += score.total_weight;//Better fits have lower score
          sum_score.num_peaks_used += score.num_peaks_used;
        } );
      }//for( size_t src_index = 0; src_index < data_srcs_ref.size(); ++src_index )

      pool.join();

      cout << "Finished individual"
      << settings.print( "  indiv" ) << endl
      << "With score " << sum_score.print("  score") << endl << endl;
    }else
    {
      for( size_t src_index = 0; src_index < data_srcs_ref.size(); ++src_index )
      {
        const DataSrcInfo &info = data_srcs_ref[src_index];
        const vector<PeakDef> &initial_peaks = intial_peaks_ref[src_index];

        const FinalFitScore score = FinalFit_GA::eval_final_peak_fit( settings, info, initial_peaks, false );

        sum_score.area_score += score.area_score;
        sum_score.width_score += score.width_score;
        sum_score.position_score += score.position_score;
        sum_score.ignored_unexpected_peaks += score.ignored_unexpected_peaks;
        sum_score.unexpected_peaks_sum_significance += score.unexpected_peaks_sum_significance;
        sum_score.total_weight += score.total_weight;//Better fits have lower score
        sum_score.num_peaks_used += score.num_peaks_used;
      }//for( size_t src_index = 0; src_index < data_srcs_ref.size(); ++src_index )
    }//if( do_multithread ) / else

    bool is_best_so_far = false;

    {//begin lock on sm_score_mutex
      std::lock_guard<std::mutex> lock( sm_score_mutex );
      is_best_so_far = (sum_score.total_weight < sm_best_score.total_weight);
      if( is_best_so_far )
      {
        cout << "Best score of:" << endl
        << sum_score.print( "\tbest_score_so_far" ) << endl
        << "For current best settings:\n"
        << settings.print("\tbest_settings_so_far") << endl
        << "For individual " << ns_individuals_eval.load()
        << " of gen " << ns_generation_num.load()
        << " weight over " << data_srcs_ref.size() << " spectra.\n"
        << "Previous best score was:\n"
        << sm_best_score.print("\tprev_best_settings") << endl
        << "\n(will write to N42 files)"
        << endl << endl;

        sm_best_score = sum_score;
      }
    }//end lock on sm_score_mutex

    // If best so far, we'll write out N42 files
    // Note however, that we arent taking a lock, so output files could be currupted garbage
    if( is_best_so_far )
    {
      //Write out best peaks so far to N42 files so we can inspect, as things are optimizing
      //we'll use a few threads to make it a little faster, but not too many, because we may be doing this from multiple threads already
      boost::asio::thread_pool pool(4);

      for( size_t src_index = 0; src_index < data_srcs_ref.size(); ++src_index )
      {
        const DataSrcInfo * const info = &(data_srcs_ref[src_index]);
        const vector<PeakDef> * const initial_peaks = &(intial_peaks_ref[src_index]);
        boost::asio::post(pool, [info,initial_peaks,&settings,isHPGe](){
          FinalFit_GA::eval_final_peak_fit( settings, *info, *initial_peaks, true );
        } );
      }//for( size_t src_index = 0; src_index < data_srcs_ref.size(); ++src_index )
      pool.join();
      std::lock_guard<std::mutex> lock( sm_score_mutex );
      cout << "Wrote output N42 files." << endl;
    }//if( sum_score.total_weight < sm_best_score.total_weight )

    cout << "Individual " << ns_individuals_eval.load()
        << " of gen " << ns_generation_num.load()
        << " weight over " << data_srcs_ref.size() << " spectra:\n"
        << sum_score.print("\tsum_score")
        << "--------" << endl
        << endl;

    ns_individuals_eval += 1;

    return sum_score.total_weight; //We are optimizing the cost so that lower is better.
  };// set InitialFit_GA::ns_ga_eval_fcn


  //cout << "Doing minuit_fit_final_pars..." << endl;
  //const FinalPeakFitSettings best_final_fit_settings = minuit_fit_final_pars( ga_eval_fcn );

  //cerr << "\n\n\nWarning - not doign genetic optimzation!!!\n\n" << endl;
  const FinalPeakFitSettings best_final_fit_settings = FinalFit_GA::do_ga_eval( ga_eval_fcn );

  //eval_candidate_settings_fcn( best_settings, best_initial_fit_settings, true );
  //cout << "Wrote N42s with best settings." << endl;

  const double end_wall = SpecUtils::get_wall_time();
  const double end_cpu = SpecUtils::get_cpu_time();
  auto now = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());

  for( size_t i = 0; i < initial_peak_fits.size(); ++i )
  {
    const DataSrcInfo src_info = input_srcs[i];
    const vector<PeakDef> &intial_peaks = initial_peak_fits[i];

    const bool write_n42 = true;
    eval_final_peak_fit( best_final_fit_settings, src_info, intial_peaks, write_n42 );
  }

  stringstream outmsg;
  outmsg << "Finish Time: " << SpecUtils::to_iso_string(now) << endl
  << endl
  << "Final settings used:\n"
  << candidate_settings.print("\tcandidate_settings") << endl
  << endl
  << initial_fit_settings.print("\tinitial_fit_settings") << endl
  << endl
#if( !WRITE_ALL_SPEC_TO_HTML ) //quick hack to get to compile
  << best_final_fit_settings.print("\tbest_final_fit_settings")
#endif
  //<< "Best settings had score " << best_sum_score << ", with values:" << endl
  << endl << endl
  << "Ran in {wall=" << (end_wall - start_wall)
  << ", cpu=" << (end_cpu - start_cpu) << "} seconds" << endl;

  {
    ofstream best_settings_file( "best_final_settings.txt" );
    best_settings_file << outmsg.str();
  }

  cout << outmsg.str()<< endl;


  //TODO:
  // - Should time using LinearProblemSubSolveChi2Fcn - may be faster/better?
  // - Also, try using Ceres instead of Minuit.  Look for `USE_LM_PEAK_FIT`
  // -Compare using PeakFitLM::fit_peak_for_user_click_LM(...) vs
  //  `refitPeaksThatShareROI(...)`
  // - See what works better for peaks - `LinearProblemSubSolveChi2Fcn` or `PeakFitChi2Fcn` - e.g. more accurate answers for some default fit scenarios.

}//void do_final_peak_fit_ga_optimization(...)



string FinalFitSolution::to_string( const string &separator ) const
{
  return
  string("require_combine_num_fwhm_near: ") + std::to_string(require_combine_num_fwhm_near)
  + separator + "not_allow_combine_num_fwhm_near: " + std::to_string(not_allow_combine_num_fwhm_near)
  //+ separator + "combine_ROI_overlap_frac: " + std::to_string(combine_ROI_overlap_frac)
  + separator + "cont_type_peak_nsigma_threshold: " + std::to_string(cont_type_peak_nsigma_threshold)
  + separator + "cont_type_left_right_nsigma: " + std::to_string(cont_type_left_right_nsigma)
  //+ separator + "stepped_roi_extent_lower_side_stat_multiple: " + std::to_string(stepped_roi_extent_lower_side_stat_multiple)
  //+ separator + "stepped_roi_extent_upper_side_stat_multiple: " + std::to_string(stepped_roi_extent_upper_side_stat_multiple)
  + separator + "cont_poly_order_increase_chi2dof_required: " + std::to_string(cont_poly_order_increase_chi2dof_required)
  + separator + "cont_step_type_increase_chi2dof_required: " + std::to_string(cont_step_type_increase_chi2dof_required)
  + separator + "skew_nsigma: " + std::to_string(skew_nsigma)
  + separator + "left_residual_sum_min_to_try_skew: " + std::to_string(left_residual_sum_min_to_try_skew)
  + separator + "right_residual_sum_min_to_try_skew: " + std::to_string(right_residual_sum_min_to_try_skew)
  + separator + "skew_improve_chi2_dof_threshold: " + std::to_string(skew_improve_chi2_dof_threshold)
  + separator + "roi_extent_low_num_fwhm_base_highstat: " + std::to_string(roi_extent_low_num_fwhm_base_highstat)
  + separator + "roi_extent_high_num_fwhm_base_highstat: " + std::to_string(roi_extent_high_num_fwhm_base_highstat)
  + separator + "roi_extent_low_num_fwhm_base_lowstat: " + std::to_string(roi_extent_low_num_fwhm_base_lowstat)
  + separator + "roi_extent_high_num_fwhm_base_lowstat: " + std::to_string(roi_extent_high_num_fwhm_base_lowstat)
  + separator + "high_stat_threshold: " + std::to_string(high_stat_threshold)
  + separator + "roi_extent_low_num_fwhm_extra: " + std::to_string(roi_extent_low_num_fwhm_extra)
  + separator + "roi_extent_high_num_fwhm_extra: " + std::to_string(roi_extent_high_num_fwhm_extra)
  + separator + "roi_end_second_deriv_thresh: " + std::to_string(roi_end_second_deriv_thresh)
  + separator + "break_multi_roi_up_continuum_away_sigma: " + std::to_string(break_multi_roi_up_continuum_away_sigma)
  + separator + "break_multi_roi_up_required_chi2dof_improve: " + std::to_string(break_multi_roi_up_required_chi2dof_improve)
  ;
}//to_string( separator )


void init_genes(FinalFitSolution& p,const std::function<double(void)> &rnd01)
{
  // rnd01() gives a random number in 0~1
  p.require_combine_num_fwhm_near = 1.0 + 7.5*rnd01();
  p.not_allow_combine_num_fwhm_near = 2.0 + 13*rnd01();
  //p.combine_ROI_overlap_frac = -1.0 + 2.0*rnd01();  // Range -1 to 1
  p.cont_type_peak_nsigma_threshold = 10.0 + 65.0*rnd01();  // Range 10-75
  p.cont_type_left_right_nsigma = 1.0 + 19.0*rnd01();  // Range 1-20
  //p.stepped_roi_extent_lower_side_stat_multiple = 0.0 + 2.0*rnd01();  // Range 0-20
  //p.stepped_roi_extent_upper_side_stat_multiple = 0.0 + 2.0*rnd01();  // Range 0-20
  p.cont_poly_order_increase_chi2dof_required = 0.0 + 1.5*rnd01();  // Range 0-1.5
  p.cont_step_type_increase_chi2dof_required = 0.0 + 1.5*rnd01();  // Range 0-1.5
  p.skew_nsigma = 4.0 + 8.0*rnd01();  // Range 0-12
  p.left_residual_sum_min_to_try_skew = 0.0 + 10.0*rnd01();  // Range 0-10
  p.right_residual_sum_min_to_try_skew = 0.0 + 5.0*rnd01();  // Range 0-5
  p.skew_improve_chi2_dof_threshold = 0.5 + 4.5*rnd01();  // Range 0.5-4.5
  p.roi_extent_low_num_fwhm_base_highstat  = 0.5 + 4.0*rnd01();  // Range 0.5-4
  p.roi_extent_high_num_fwhm_base_highstat = 0.5 + 4.0*rnd01();  // Range 0.5-4
  p.roi_extent_low_num_fwhm_base_lowstat  = 0.5 + 4.0*rnd01();  // Range 0.5-4
  p.roi_extent_high_num_fwhm_base_lowstat = 0.5 + 4.0*rnd01();  // Range 0.5-4
  p.high_stat_threshold = 5.0 + 95.0*rnd01();  // Range 5-100
  p.roi_extent_low_num_fwhm_extra  = 0.0 + 5.0*rnd01();  // Range 0 - 5
  p.roi_extent_high_num_fwhm_extra = 0.0 + 5.0*rnd01();  // Range 0 - 5
  p.roi_end_second_deriv_thresh    = 0.0 + 20.0*rnd01();  // Range 0 - 20 - total guess
  p.break_multi_roi_up_continuum_away_sigma = 0.0 + 15.0*rnd01();  // Range 0 - 15 - total guess - if less than 1, then requirement wont be applied
  p.break_multi_roi_up_required_chi2dof_improve = 0.1 + 2.9*rnd01(); //Range 0.1 - 3 - total guess
}


FinalPeakFitSettings genes_to_settings( const FinalFitSolution &solution )
{
  FinalPeakFitSettings settings;

  settings.require_combine_num_fwhm_near = solution.require_combine_num_fwhm_near;
  settings.not_allow_combine_num_fwhm_near = solution.not_allow_combine_num_fwhm_near;
  //settings.combine_ROI_overlap_frac = solution.combine_ROI_overlap_frac;
  settings.cont_type_peak_nsigma_threshold = solution.cont_type_peak_nsigma_threshold;
  settings.cont_type_left_right_nsigma = solution.cont_type_left_right_nsigma;
  //settings.stepped_roi_extent_lower_side_stat_multiple = solution.stepped_roi_extent_lower_side_stat_multiple;
  //settings.stepped_roi_extent_upper_side_stat_multiple = solution.stepped_roi_extent_upper_side_stat_multiple;
  settings.cont_poly_order_increase_chi2dof_required = solution.cont_poly_order_increase_chi2dof_required;
  settings.cont_step_type_increase_chi2dof_required = solution.cont_step_type_increase_chi2dof_required;
  settings.skew_nsigma = solution.skew_nsigma;
  settings.left_residual_sum_min_to_try_skew = solution.left_residual_sum_min_to_try_skew;
  settings.right_residual_sum_min_to_try_skew = solution.right_residual_sum_min_to_try_skew;
  settings.skew_improve_chi2_dof_threshold = solution.skew_improve_chi2_dof_threshold;
  settings.roi_extent_low_num_fwhm_base_highstat = solution.roi_extent_low_num_fwhm_base_highstat;
  settings.roi_extent_high_num_fwhm_base_highstat = solution.roi_extent_high_num_fwhm_base_highstat;
  settings.roi_extent_low_num_fwhm_base_lowstat = solution.roi_extent_low_num_fwhm_base_lowstat;
  settings.roi_extent_high_num_fwhm_base_lowstat = solution.roi_extent_high_num_fwhm_base_lowstat;
  settings.high_stat_threshold = solution.high_stat_threshold;
  settings.roi_extent_low_num_fwhm_extra = solution.roi_extent_low_num_fwhm_extra;
  settings.roi_extent_high_num_fwhm_extra = solution.roi_extent_high_num_fwhm_extra;
  settings.roi_end_second_deriv_thresh = solution.roi_end_second_deriv_thresh;
  settings.break_multi_roi_up_continuum_away_sigma = solution.break_multi_roi_up_continuum_away_sigma;
  settings.break_multi_roi_up_required_chi2dof_improve = solution.break_multi_roi_up_required_chi2dof_improve;

  return settings;
}


bool eval_solution( const FinalFitSolution &p, FinalFitCost &c )
{
  const FinalPeakFitSettings settings = genes_to_settings( p );

  try
  {
    c.objective1 = ns_ga_eval_fcn( settings );

    if( IsInf(c.objective1) || IsNan(c.objective1) )
    {
      cerr << "Value of c.objective1 is " << c.objective1 << " - will set to large positive value" << endl
      << settings.print("    bad_settings") << endl;

      c.objective1 = std::numeric_limits<double>::max();
      return false;
    }else
    {
      cout << "\n\nEvaluation gave " << c.objective1 << " for\n"
      << settings.print("    settings") << endl << endl;
    }
  }catch( std::exception &e )
  {
    cerr << "Failed to eval genes: " << e.what() << endl;
    return false;
  }

  return true;
}


FinalFitSolution mutate(
                        const FinalFitSolution& X_base,
                        const std::function<double(void)> &rnd01,
                        double shrink_scale)
{
  FinalFitSolution X_new;
  const double mu = 0.2*shrink_scale;  // mutation radius

  auto mutate_value = [&rnd01,mu](double val, double min_val, double max_val) {

    const double mutate_threshold = PeakFitImprove::sm_ga_mutate_threshold;
    if( rnd01() > mutate_threshold )
      return val;

    const double r = 2.0*(rnd01() - 0.5);
    val += r*mu*(max_val - min_val);
    return std::max(min_val, std::min(max_val, val));
  };

  X_new.require_combine_num_fwhm_near = mutate_value(X_base.require_combine_num_fwhm_near, 1.0, 8.5);
  X_new.not_allow_combine_num_fwhm_near = mutate_value(X_base.not_allow_combine_num_fwhm_near, 2.0, 15.0);
  //X_new.combine_ROI_overlap_frac = mutate_value(X_base.combine_ROI_overlap_frac, -1.0, 1.0);
  X_new.cont_type_peak_nsigma_threshold = mutate_value(X_base.cont_type_peak_nsigma_threshold, 10.0, 75.0);
  X_new.cont_type_left_right_nsigma = mutate_value(X_base.cont_type_left_right_nsigma, 1.0, 20.0);
  //X_new.stepped_roi_extent_lower_side_stat_multiple = mutate_value(X_base.stepped_roi_extent_lower_side_stat_multiple, 0.0, 2.0);
  //X_new.stepped_roi_extent_upper_side_stat_multiple = mutate_value(X_base.stepped_roi_extent_upper_side_stat_multiple, 0.0, 2.0);
  X_new.cont_poly_order_increase_chi2dof_required = mutate_value(X_base.cont_poly_order_increase_chi2dof_required, 0.0, 1.5);
  X_new.cont_step_type_increase_chi2dof_required = mutate_value(X_base.cont_step_type_increase_chi2dof_required, 0.0, 1.5);
  X_new.skew_nsigma = mutate_value(X_base.skew_nsigma, 4.0, 12.0);
  X_new.left_residual_sum_min_to_try_skew = mutate_value(X_base.left_residual_sum_min_to_try_skew, 0.0, 10.0);
  X_new.right_residual_sum_min_to_try_skew = mutate_value(X_base.right_residual_sum_min_to_try_skew, 0.0, 5.0);
  X_new.skew_improve_chi2_dof_threshold = mutate_value(X_base.skew_improve_chi2_dof_threshold, 0.5, 5.0);
  X_new.roi_extent_low_num_fwhm_base_highstat = mutate_value(X_base.roi_extent_low_num_fwhm_base_highstat, 0.5, 4.0);
  X_new.roi_extent_high_num_fwhm_base_highstat = mutate_value(X_base.roi_extent_high_num_fwhm_base_highstat, 0.5, 4.0);
  X_new.roi_extent_low_num_fwhm_base_lowstat = mutate_value(X_base.roi_extent_low_num_fwhm_base_lowstat, 0.5, 4.0);
  X_new.roi_extent_high_num_fwhm_base_lowstat = mutate_value(X_base.roi_extent_high_num_fwhm_base_lowstat, 0.5, 4.0);
  X_new.high_stat_threshold = mutate_value(X_base.high_stat_threshold, 5.0, 100.0);
  X_new.roi_extent_low_num_fwhm_extra = mutate_value(X_base.roi_extent_low_num_fwhm_extra, 0.0, 5.0);
  X_new.roi_extent_high_num_fwhm_extra = mutate_value(X_base.roi_extent_high_num_fwhm_extra, 0.0, 5.0);
  X_new.roi_end_second_deriv_thresh = mutate_value(X_base.roi_end_second_deriv_thresh, 0.0, 20.0);
  X_new.break_multi_roi_up_continuum_away_sigma = mutate_value(X_base.break_multi_roi_up_continuum_away_sigma, 0.0, 15.0);
  X_new.break_multi_roi_up_required_chi2dof_improve = mutate_value(X_base.break_multi_roi_up_required_chi2dof_improve, 0.1, 3.0);

  return X_new;
}


FinalFitSolution crossover( const FinalFitSolution& X1,
                            const FinalFitSolution& X2,
                            const std::function<double(void)> &rnd01 )
{
  FinalFitSolution X_new;

  auto cross_value = [&rnd01](double v1, double v2) {
    const double crossover_threshold = PeakFitImprove::sm_ga_crossover_threshold;
    if( rnd01() > crossover_threshold )
      return (rnd01() < 0.5) ? v1 : v2;

    const double r = rnd01();
    return r*v1 + (1.0-r)*v2;
  };

  X_new.require_combine_num_fwhm_near = cross_value(X1.require_combine_num_fwhm_near, X2.require_combine_num_fwhm_near);
  X_new.not_allow_combine_num_fwhm_near = cross_value(X1.not_allow_combine_num_fwhm_near, X2.not_allow_combine_num_fwhm_near);
  //X_new.combine_ROI_overlap_frac = cross_value(X1.combine_ROI_overlap_frac, X2.combine_ROI_overlap_frac);
  X_new.cont_type_peak_nsigma_threshold = cross_value(X1.cont_type_peak_nsigma_threshold, X2.cont_type_peak_nsigma_threshold);
  X_new.cont_type_left_right_nsigma = cross_value(X1.cont_type_left_right_nsigma, X2.cont_type_left_right_nsigma);
  //X_new.stepped_roi_extent_lower_side_stat_multiple = cross_value(X1.stepped_roi_extent_lower_side_stat_multiple, X2.stepped_roi_extent_lower_side_stat_multiple);
  //X_new.stepped_roi_extent_upper_side_stat_multiple = cross_value(X1.stepped_roi_extent_upper_side_stat_multiple, X2.stepped_roi_extent_upper_side_stat_multiple);
  X_new.cont_poly_order_increase_chi2dof_required = cross_value(X1.cont_poly_order_increase_chi2dof_required, X2.cont_poly_order_increase_chi2dof_required);
  X_new.cont_step_type_increase_chi2dof_required = cross_value(X1.cont_step_type_increase_chi2dof_required, X2.cont_step_type_increase_chi2dof_required);
  X_new.skew_nsigma = cross_value(X1.skew_nsigma, X2.skew_nsigma);
  X_new.left_residual_sum_min_to_try_skew = cross_value(X1.left_residual_sum_min_to_try_skew, X2.left_residual_sum_min_to_try_skew);
  X_new.right_residual_sum_min_to_try_skew = cross_value(X1.right_residual_sum_min_to_try_skew, X2.right_residual_sum_min_to_try_skew);
  X_new.skew_improve_chi2_dof_threshold = cross_value(X1.skew_improve_chi2_dof_threshold, X2.skew_improve_chi2_dof_threshold);
  X_new.roi_extent_low_num_fwhm_base_highstat = cross_value(X1.roi_extent_low_num_fwhm_base_highstat, X2.roi_extent_low_num_fwhm_base_highstat);
  X_new.roi_extent_high_num_fwhm_base_highstat = cross_value(X1.roi_extent_high_num_fwhm_base_highstat, X2.roi_extent_high_num_fwhm_base_highstat);
  X_new.roi_extent_low_num_fwhm_base_lowstat = cross_value(X1.roi_extent_low_num_fwhm_base_lowstat, X2.roi_extent_low_num_fwhm_base_lowstat);
  X_new.roi_extent_high_num_fwhm_base_lowstat = cross_value(X1.roi_extent_high_num_fwhm_base_lowstat, X2.roi_extent_high_num_fwhm_base_lowstat);
  X_new.high_stat_threshold = cross_value(X1.high_stat_threshold, X2.high_stat_threshold);
  X_new.roi_extent_low_num_fwhm_extra = cross_value(X1.roi_extent_low_num_fwhm_extra, X2.roi_extent_low_num_fwhm_extra);
  X_new.roi_extent_high_num_fwhm_extra = cross_value(X1.roi_extent_high_num_fwhm_extra, X2.roi_extent_high_num_fwhm_extra);
  X_new.roi_end_second_deriv_thresh = cross_value(X1.roi_end_second_deriv_thresh, X2.roi_end_second_deriv_thresh);
  X_new.break_multi_roi_up_continuum_away_sigma = cross_value(X1.break_multi_roi_up_continuum_away_sigma, X2.break_multi_roi_up_continuum_away_sigma);
  X_new.break_multi_roi_up_required_chi2dof_improve = cross_value(X1.break_multi_roi_up_required_chi2dof_improve, X2.break_multi_roi_up_required_chi2dof_improve);

  return X_new;
}


double calculate_SO_total_fitness(const GA_Type::thisChromosomeType &X)
{
  // finalize the cost
  double final_cost = 0.0;
  final_cost = X.middle_costs.objective1;

  if( IsInf(final_cost) || IsNan(final_cost) )
  {
    cerr << "calculate_SO_total_fitness: final-cost of is " << final_cost << endl;
  }

  return final_cost;
}


void SO_report_generation(
                          int generation_number,
                          const EA::GenerationType<FinalFitSolution,FinalFitCost> &last_generation,
                          const FinalFitSolution& best_genes)
{
  cout<<"\n\nGeneration ["<<generation_number<<"], ";
  cout<<"Best = "<<last_generation.best_total_cost<<", ";
  cout<<"Average = "<<last_generation.average_cost<<", ";

  cout<<", genes = "<<best_genes.to_string(",\n\t") << endl << endl;

  ns_generation_num += 1;
  ns_individuals_eval = 0;
}


FinalPeakFitSettings do_ga_eval( std::function<double(const FinalPeakFitSettings &)> ga_eval_fcn )
{
  assert( !sm_has_been_called );
  if( sm_has_been_called )
  {
    cerr << "You may only call FinalFit_GA::do_ga_eval(...) once per program execution." << endl;
    exit(1);
  }

  sm_has_been_called = true;

  EA::Chronometer timer;
  timer.tic();

  GA_Type ga_obj;
  ga_obj.problem_mode = EA::GA_MODE::SOGA;

  ga_obj.multi_threading = true;
  if( PeakFitImprove::debug_printout )
  {
    ga_obj.multi_threading = false;
  }else
  {
    ga_obj.N_threads = static_cast<int>(PeakFitImprove::sm_num_optimization_threads);
    cout << "Will use " << ga_obj.N_threads << " threads." << endl;
  }

  ga_obj.idle_delay_us = 1; // switch between threads quickly
  ga_obj.verbose = false;
  ga_obj.population = static_cast<unsigned int>(PeakFitImprove::sm_ga_population);
  ga_obj.generation_max = static_cast<int>(PeakFitImprove::sm_ga_generation_max);
  ga_obj.calculate_SO_total_fitness = calculate_SO_total_fitness;
  ga_obj.init_genes = init_genes;
  ga_obj.eval_solution = eval_solution;
  ga_obj.mutate = mutate;
  ga_obj.crossover = crossover;
  ga_obj.SO_report_generation = SO_report_generation;
  ga_obj.best_stall_max = static_cast<int>(PeakFitImprove::sm_ga_best_stall_max);
  ga_obj.elite_count = static_cast<int>(PeakFitImprove::sm_ga_elite_count);
  ga_obj.crossover_fraction = PeakFitImprove::sm_ga_crossover_fraction;
  ga_obj.mutation_rate = PeakFitImprove::sm_ga_mutation_rate;

  ns_ga_eval_fcn = ga_eval_fcn;

  ga_obj.solve();

  cout<<"The problem is optimized in "<<timer.toc()<<" seconds."<<endl;

  return genes_to_settings( ga_obj.last_generation.chromosomes[ga_obj.last_generation.best_chromosome_index].genes );
}
};//namespace FinalFit_GA
