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

#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/SpecUtilsAsync.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/DetectorPeakResponse.h"


#include "PeakFitImprove.h"
#include "FinalFit_GA.h"
#include "InitialFit_GA.h"

using namespace std;



namespace
{
std::function<double(const FinalPeakFitSettings &)> ns_ga_eval_fcn;

bool sm_has_been_called = false;
}


namespace FinalFit_GA
{
vector<PeakDef> final_peak_fit_for_roi( const vector<PeakDef> &pre_fit_peaks,
                               const FinalPeakFitSettings &final_fit_settings,
                               const DataSrcInfo &src_info )
{
#ifndef NDEBUG
  // Make sure all peaks are from the same ROI - only need to do in debug builds
  for( size_t i = 1; i < pre_fit_peaks.size(); ++i )
  {
    assert( pre_fit_peaks[0].continuum() == pre_fit_peaks[i].continuum() );
  }
#endif

  assert( !pre_fit_peaks.empty() );
  if( pre_fit_peaks.empty() )
    return vector<PeakDef>();

  const vector<shared_ptr<const SpecUtils::Measurement>> &src_spectra = src_info.src_info.src_spectra;
  assert( !src_spectra.empty() );
  const shared_ptr<const SpecUtils::Measurement> &data = src_spectra.front();
  assert( data );

  const std::shared_ptr<const PeakContinuum> orig_continuum = pre_fit_peaks.front().continuum();


  // Set ROI widthn based on
   //Single peak ROIS: Width will be modified based on stat uncert of initial fit.
  //double roi_extent_low_num_fwhm_base, roi_extent_high_num_fwhm_base;

  const PeakDef &first_peak = pre_fit_peaks.front();
  const PeakDef &last_peak = pre_fit_peaks.back();
  const double first_mean = first_peak.mean();
  const double last_mean = last_peak.mean();
  const double first_fwhm = first_peak.fwhm();
  const double last_fwhm = last_peak.fwhm();
  const double base_roi_lower = first_mean - final_fit_settings.roi_extent_low_num_fwhm_base*first_fwhm;
  const double base_roi_upper = last_mean + final_fit_settings.roi_extent_high_num_fwhm_base*last_fwhm;

  double roi_lower, roi_upper;

  if( pre_fit_peaks.size() == 1 )
  {
    // Single peak ROI
    const double central_continuum_area = orig_continuum->offset_integral( first_mean - 0.5*first_fwhm, first_mean + 0.5*last_fwhm, data );
    const double central_peak_area = first_peak.gauss_integral( first_mean - 0.5*first_fwhm, first_mean + 0.5*last_fwhm );
    double central_ratio = central_continuum_area / central_peak_area;
    central_ratio = std::max( 0.0, central_ratio );
    central_ratio = std::min( 1.0, central_ratio );

    switch( final_fit_settings.roi_extent_mult_type )
    {
      case FinalPeakFitSettings::RoiExtentMultType::Linear:
        roi_lower = base_roi_lower - final_fit_settings.roi_extent_lower_side_stat_multiple*central_ratio*first_peak.fwhm();
        roi_upper = base_roi_upper + final_fit_settings.roi_extent_upper_side_stat_multiple*central_ratio*first_peak.fwhm();
        break;

      case FinalPeakFitSettings::RoiExtentMultType::Sqrt:
        roi_lower = base_roi_lower - final_fit_settings.roi_extent_lower_side_stat_multiple*sqrt(central_ratio)*first_peak.fwhm();
        roi_upper = base_roi_upper + final_fit_settings.roi_extent_upper_side_stat_multiple*sqrt(central_ratio)*first_peak.fwhm();
        break;
    }//switch( final_fit_settings.roi_extent_mult_type )
  }else
  {
    // Multi-peak ROI
    roi_lower = base_roi_lower - final_fit_settings.multi_roi_extent_lower_side_fwhm_mult*first_peak.fwhm();
    roi_upper = base_roi_upper + final_fit_settings.multi_roi_extent_upper_side_fwhm_mult*last_peak.fwhm();
  }//if( pre_fit_peaks.size() == 1 ) / else

  vector<PeakDef> initial_peaks = pre_fit_peaks;
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



  // Do an initial fit
  const double initial_stat_threshold = 0.0;
  const double initial_hypothesis_threshold = 0.0;
  const bool amplitudeOnly = false;
  const bool isHPGe = true;
  #warning "final_peak_fit_for_roi: always assuming HPGe right now - for dev"

  vector<PeakDef> first_fit_peaks;
#if( USE_CERES_PEAK_FITTING )
  vector<shared_ptr<const PeakDef>> results_tmp, input_peaks_tmp;
  for( const PeakDef &p : initial_peaks )
    input_peaks_tmp.push_back( make_shared<PeakDef>(p) );
  PeakFitLM::fit_peaks_LM( results_tmp, input_peaks_tmp, data,
                          initial_stat_threshold, initial_hypothesis_threshold,  amplitudeOnly, isHPGe );
  for( const shared_ptr<const PeakDef> &p : results_tmp )
    first_fit_peaks.push_back( *p );
#else
  fitPeaks( first_fit_peaks, initial_stat_threshold, initial_hypothesis_threshold,
               data, initial_peaks, amplitudeOnly, isHPGe );
#endif

  assert( !first_fit_peaks.empty() );
  if( first_fit_peaks.empty() )
    throw std::runtime_error( "final_peak_fit_for_roi: first_fit_peaks is empty" );

  const double initial_fit_chi2 = chi2_for_region_wrapper( first_fit_peaks, data, compare_lower_channel, compare_upper_channel );
  const double approx_dof = (compare_upper_channel - compare_lower_channel);
  cout << "Initial fit chi2: " << initial_fit_chi2 << " compared to intial_chi2: " << intial_chi2 << endl;

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
    explore_cont_type |= ( nsigma > final_fit_settings.cont_type_peak_nsigma_threshold );
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
  if( !starting_with_step_continuum )
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
    const double average_uncert = sqrt( average_lower_sum*average_lower_sum + average_upper_sum*average_upper_sum );
    if( average_uncert > 0.0 )
    {
      const double average_nsigma = (average_lower_sum - average_upper_sum) / average_uncert;
      explore_stepped_cont = (average_nsigma >= final_fit_settings.cont_type_left_right_nsigma);
    }
  }// End of check if we should explore stepped cont

  if( explore_stepped_cont )
  {
    cont_types_to_check.push_back( PeakContinuum::OffsetType::FlatStep );
    cont_types_to_check.push_back( PeakContinuum::OffsetType::LinearStep );

    // I dont find bilinear steps to be very good, at least visually, so we will not explore them
    //cont_types_to_check.push_back( PeakContinuum::OffsetType::BiLinearStep );
  }else if( starting_with_step_continuum && (orig_continuum->type() != PeakContinuum::OffsetType::BiLinearStep) )
  {
    cont_types_to_check.push_back( PeakContinuum::OffsetType(orig_continuum->type() + 1) );
  }

  vector<tuple<PeakContinuum::OffsetType,vector<PeakDef>,double>> trial_peak_fits;
  trial_peak_fits.push_back( make_tuple( orig_continuum->type(), pre_fit_peaks, initial_fit_chi2/approx_dof )  );

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
#if( USE_CERES_PEAK_FITTING )
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
    }catch( const std::exception &e )
    {
      cout << "Error fitting peaks: " << e.what() << endl;
      continue;
    }//try / catch to fit peaks

    assert( !these_fit_peaks.empty() );
    if( these_fit_peaks.empty() )
      continue;

    const double these_fit_chi2 = chi2_for_region_wrapper( these_fit_peaks, data, compare_lower_channel, compare_upper_channel );
    trial_peak_fits.push_back( make_tuple( type, these_input_peaks, these_fit_chi2/approx_dof ) );
  }//for( const PeakContinuum::OffsetType type : cont_types_to_check )

  // I'm pretty sure trial_peak_fits is sorted by continuum type, but just to be sure, we will sort it
  std::sort( begin(trial_peak_fits), end(trial_peak_fits), []( const auto &a, const auto &b ){
    return static_cast<int>(std::get<0>(a)) < static_cast<int>(std::get<0>(b));
  } );

  size_t best_index = 0;
  for( size_t i = 1; i < trial_peak_fits.size(); ++i )
  {
    const bool prev_was_step = PeakContinuum::is_step_continuum( std::get<0>(trial_peak_fits[best_index]) );
    const bool this_is_step = PeakContinuum::is_step_continuum( std::get<0>(trial_peak_fits[i]) );

    const double prev_chi2 = std::get<2>(trial_peak_fits[best_index]);
    const double this_chi2 = std::get<2>(trial_peak_fits[i]);

    if( !prev_was_step && !this_is_step )
    {
      if( (this_chi2 - final_fit_settings.cont_poly_order_increase_chi2dof_required) > prev_chi2 )
        best_index = i;
    }else if( prev_was_step && !this_is_step )
    {
      assert( 0 );  // We should never get here
    }else if( (!prev_was_step && this_is_step) || (prev_was_step && this_is_step) )
    {
      if( (this_chi2 - final_fit_settings.cont_step_type_increase_chi2dof_required) > prev_chi2 )
        best_index = i;
    }
  }//for( size_t i = 1; i < trial_peak_fits.size(); ++i )

  assert( std::get<0>(trial_peak_fits[0]) == orig_continuum->type() );

  cout << "Post continuum type fit: Best index: " << best_index << " gives Chi2Dof: " << std::get<2>(trial_peak_fits[best_index])
  << " (ContType: " << PeakContinuum::offset_type_str( std::get<0>(trial_peak_fits[best_index]) ) << ")"
  << ", compared to initial Chi2Dof: " << initial_fit_chi2/approx_dof
  << " (ContType: " << PeakContinuum::offset_type_str( std::get<0>(trial_peak_fits[0]) ) << ")" << endl;

  const double after_cont_type_chi2dof = std::get<2>(trial_peak_fits[best_index]);
  const double after_cont_type_chi2 = after_cont_type_chi2dof * approx_dof;
  const vector<PeakDef> best_post_cont_type_peaks = std::get<1>(trial_peak_fits[best_index]);
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
    assert( !too_narrow );

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
#if( USE_CERES_PEAK_FITTING )
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
        }catch( const std::exception &e )
        {
          cout << "Error adding skew to peaks: " << e.what() << endl;
          assert( 0 );
          continue;
        }//try / catch to fit peaks

        const double these_fit_chi2 = chi2_for_region_wrapper( these_fit_peaks, data, compare_lower_channel, compare_upper_channel );
        trial_skew_fits.emplace_back( skew_type, these_input_peaks, these_fit_chi2/approx_dof );
      }//for( const PeakDef::SkewType skew_type : skew_types_to_try )

      // Sort the skew fits by the skew type (I'm pretty sure they should already be sorted)
      std::sort( begin(trial_skew_fits), end(trial_skew_fits), []( const auto &a, const auto &b ){
        return std::get<0>(a) < std::get<0>(b);
      } );

      size_t best_skew_index = 0;
      for( size_t i = 1; i < trial_skew_fits.size(); ++i )
      {
        const double prev_chi2 = std::get<2>(trial_skew_fits[best_skew_index]);
        const double this_chi2 = std::get<2>(trial_skew_fits[i]);

        if( (this_chi2 - final_fit_settings.skew_improve_chi2_dof_threshold) > prev_chi2 )
          best_skew_index = i;
      }//for( size_t i = 1; i < trial_skew_fits.size(); ++i )

      cout << "Post skew fit: Best index: " << best_skew_index << " gives Chi2Dof: " << std::get<2>(trial_skew_fits[best_skew_index])
      << " (SkewType: " << PeakDef::to_string( std::get<0>(trial_skew_fits[best_skew_index]) ) << ")"
      << ", compared to initial Chi2Dof: " << initial_fit_chi2/approx_dof
      << " (SkewType: " << PeakDef::to_string( std::get<0>(trial_skew_fits[0]) ) << ")" << endl;

      answer = std::get<1>(trial_skew_fits[best_skew_index]);
    }//if( !too_narrow && (most_sig_peak.skew_type() == PeakDef::SkewType::NoSkew) )
  }//if( post_cont_type_max_nsigma > final_fit_settings.skew_nsigma )

  return answer;
}//vector<PeakDef> final_peak_fit_for_roi(...)


vector<PeakDef> final_peak_fit( const vector<PeakDef> &pre_fit_peaks,
                               const FinalPeakFitSettings &final_fit_settings,
                               const DataSrcInfo &src_info )
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


  vector<pair<shared_ptr<const PeakContinuum>,vector<PeakDef>>> connected_rois;
  for( size_t roi_index = 1; roi_index < rois.size(); ++roi_index )
  {
    const pair<shared_ptr<const PeakContinuum>,vector<PeakDef>> &prev_roi = rois[roi_index-1];
    const shared_ptr<const PeakContinuum> &prev_cont = prev_roi.first;

    const pair<shared_ptr<const PeakContinuum>,vector<PeakDef>> &this_roi = rois[roi_index];
    const shared_ptr<const PeakContinuum> &this_cont = this_roi.first;
    const vector<PeakDef> &this_peaks = this_roi.second;

    bool combine_for_overlap = false, combine_for_nsigma_near = false;
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
      combine_for_overlap = ((overlap / lesser_roi_extent) > final_fit_settings.combine_ROI_overlap_frac);
    }//if( overlap_start <= overlap_end )

    const PeakDef &last_prev_peak = prev_roi.second.back();
    const PeakDef &first_this_peak = prev_roi.second.front();
    const double avrg_sigma = 0.5*(last_prev_peak.sigma() + first_this_peak.sigma());
    const double energy_diff = fabs( last_prev_peak.mean() - first_this_peak.mean() );
    combine_for_nsigma_near = ((energy_diff / avrg_sigma) < final_fit_settings.combine_nsigma_near);

    if( combine_for_overlap || combine_for_nsigma_near )
    {
      //Combine ROIs into a single ROI
      shared_ptr<PeakContinuum> new_cont = make_shared<PeakContinuum>( *prev_cont );
      const double new_lower = std::min(prev_lower, this_lower); //just to make sure nothings wonky
      const double new_upper = std::max(prev_upper, this_upper); //just to make sure nothings wonky
      new_cont->setRange( new_lower, std::max(prev_upper, new_upper) );
      const PeakContinuum::OffsetType new_offset_type = std::max( prev_cont->type(), this_cont->type() );
      if( new_offset_type != new_cont->type() )
      {
        assert( !src_info.src_info.src_spectra.empty() );
        const shared_ptr<const SpecUtils::Measurement> &data = src_info.src_info.src_spectra.front();
        assert( data );

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

      // Decrement the loop index, since now the next ROI is one element closer
      roi_index -= 1;
    }//if( combine_for_overlap || combine_for_nsigma_near )
  }//for( size_t size_t roi_index = 1; roi_index < rois.size(); ++i )


  for( const auto &cont_peaks : rois )
  {
    const vector<PeakDef> &in_peaks = cont_peaks.second;

    const vector<PeakDef> roi_answer = final_peak_fit_for_roi( in_peaks, final_fit_settings, src_info );
    answer.insert( end(answer), begin(roi_answer), end(roi_answer) );
  }

  std::sort( begin(answer), end(answer), &PeakDef::lessThanByMean );

  return answer;
}//vector<PeakDef> final_peak_fit


FinalFitScore eval_final_peak_fit( const FinalPeakFitSettings &final_fit_settings,
                           const DataSrcInfo &src_info,
                           const vector<PeakDef> &intial_peaks )
{
  FinalFitScore score;
  const vector<PeakDef> fit_peaks = final_peak_fit( intial_peaks, final_fit_settings, src_info );

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

  for( const PeakDef &found_peak : fit_peaks )
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

    if( !nearest_expected_peak )
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
        continue;

      // If this is a significant peak, we messed up, so lets ignore it...
      const double area_uncert = found_peak.amplitudeUncert() > 0.0 ? found_peak.amplitudeUncert() : 1/sqrt(found_peak.amplitude());
      if( found_peak.amplitude() > 8.0*area_uncert )
      {
        cout << "Found peak {mean: " << found_peak.mean()
        << ", fwhm: " << found_peak.fwhm()
        << ", area: " << found_peak.amplitude()
        << ", area_uncert: " << found_peak.amplitudeUncert() << "} that is unexpected, but too large, so ignoring." << endl;
        continue;
      }

      if( found_peak.amplitude() < 1.0 ) //Not a real peak, so ignore it.
        continue;

      score.ignored_unexpected_peaks += 1;
      score.unexpected_peaks_sum_significance += std::min( 7.5, found_peak.amplitude() / area_uncert ); //Cap at 7.5 sigma (arbitrary)
      assert( 0 ); // To force a check here...
    }//if( !nearest_expected_peak )

    score.num_peaks_used += 1;

    const double expected_energy = nearest_expected_peak->effective_energy;
    const double expected_fwhm = nearest_expected_peak->effective_fwhm;
    const double expected_sigma = expected_fwhm/2.35482;
    const double expected_area = nearest_expected_peak->peak_area;

    const double area_score = fabs(found_peak.amplitude() - expected_area) / sqrt(expected_area);
    const double width_score = fabs( expected_fwhm - found_fwhm ) / expected_sigma;
    const double position_score = fabs(found_energy - expected_energy) / expected_sigma;

    score.area_score += std::min( area_score, 20.0 );
    score.width_score += std::min( width_score, 1.0 );
    score.position_score += std::min( position_score, 1.5 );
  }//for( const PeakDef &found_peak : initial_peaks )

  score.total_weight = score.area_score / score.num_peaks_used;
  // TODO: we could/should add in weight for peak area_score, width_score, and position_score

  return score;
}//double eval_final_peak_fit(...)


void do_final_peak_fit_ga_optimization( const FindCandidateSettings &candidate_settings,
                           const InitialPeakFindSettings &initial_fit_settings,
                           const vector<DataSrcInfo> &input_srcs )
{
  vector<vector<PeakDef>> initial_peak_fits( input_srcs.size(), vector<PeakDef>{} );

  {// Begin get initial peak fit
    cout << "Starting to fit initial peaks." << endl;
    const double start_wall = SpecUtils::get_wall_time();
    const double start_cpu = SpecUtils::get_cpu_time();
    SpecUtilsAsync::ThreadPool pool;

    for( size_t input_src_index = 0; input_src_index < input_srcs.size(); ++input_src_index )
    {
      const DataSrcInfo &src = input_srcs[input_src_index];
      vector<PeakDef> *result = &(initial_peak_fits[input_src_index]);
      shared_ptr<const SpecUtils::Measurement> data = src.src_info.src_spectra[0];

      auto fit_initial_peaks_worker = [&candidate_settings, &initial_fit_settings, data, result](){
        size_t dummy1, dummy2;
        *result = InitialFit_GA::initial_peak_find_and_fit( initial_fit_settings, candidate_settings, data, dummy1, dummy2 );
      }; //fit_initial_peaks_worker

      pool.post( fit_initial_peaks_worker );
    }//for( loop over input_srcs )

    pool.join();

    const double end_wall = SpecUtils::get_wall_time();
    const double end_cpu = SpecUtils::get_cpu_time();
    cout << "Finished fitting initial peaks - took"
    << " {wall: " << (end_wall - start_wall)
    << " s, cpu: " << (end_cpu - start_cpu) << " s}."
    << endl;
  }// end get initial peak fit

  cout << "Will use candidate and initial fit settings:" << endl;
  candidate_settings.print( "\tcandidate_settings" );
  initial_fit_settings.print( "\tinitial_fit_settings" );
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

    double score_sum = 0.0;
    for( size_t src_index = 0; src_index < data_srcs_ref.size(); ++src_index )
    {
      const DataSrcInfo &info = data_srcs_ref[src_index];
      const vector<PeakDef> &initial_peaks = intial_peaks_ref[src_index];

      const double score = FinalFit_GA::eval_final_peak_fit( settings, info, initial_peaks ).total_weight; //Better fits have lower score
      score_sum += score;
    }

    return score_sum; //We are optimizing the cost so that lower is better.
  };// set InitialFit_GA::ns_ga_eval_fcn

#if( !WRITE_ALL_SPEC_TO_HTML ) //quick hack to get to compile
  const FinalPeakFitSettings best_final_fit_settings = FinalFit_GA::do_ga_eval( ga_eval_fcn );
#endif

  //eval_candidate_settings_fcn( best_settings, best_initial_fit_settings, true );
  //cout << "Wrote N42s with best settings." << endl;

  const double end_wall = SpecUtils::get_wall_time();
  const double end_cpu = SpecUtils::get_cpu_time();
  auto now = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());


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
  string("combine_nsigma_near: ") + std::to_string(combine_nsigma_near)
  + separator + "combine_ROI_overlap_frac: " + std::to_string(combine_ROI_overlap_frac)
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
  + separator + "roi_extent_low_num_fwhm_base: " + std::to_string(roi_extent_low_num_fwhm_base)
  + separator + "roi_extent_high_num_fwhm_base: " + std::to_string(roi_extent_high_num_fwhm_base)
  + separator + "roi_extent_mult_type: " + std::to_string(roi_extent_mult_type)
  + separator + "roi_extent_lower_side_stat_multiple: " + std::to_string(roi_extent_lower_side_stat_multiple)
  + separator + "roi_extent_upper_side_stat_multiple: " + std::to_string(roi_extent_upper_side_stat_multiple)
  + separator + "multi_roi_extent_lower_side_fwhm_mult: " + std::to_string(multi_roi_extent_lower_side_fwhm_mult)
  + separator + "multi_roi_extent_upper_side_fwhm_mult: " + std::to_string(multi_roi_extent_upper_side_fwhm_mult)
  ;
}//to_string( separator )


void init_genes(FinalFitSolution& p,const std::function<double(void)> &rnd01)
{
  // rnd01() gives a random number in 0~1
  p.combine_nsigma_near = 1.0 + 9.0*rnd01();  // Range 1-10
  p.combine_ROI_overlap_frac = -1.0 + 2.0*rnd01();  // Range -1 to 1
  p.cont_type_peak_nsigma_threshold = 5.0 + 55.0*rnd01();  // Range 5-60
  p.cont_type_left_right_nsigma = 1.0 + 19.0*rnd01();  // Range 1-20
  //p.stepped_roi_extent_lower_side_stat_multiple = 0.0 + 2.0*rnd01();  // Range 0-20
  //p.stepped_roi_extent_upper_side_stat_multiple = 0.0 + 2.0*rnd01();  // Range 0-20
  p.cont_poly_order_increase_chi2dof_required = 0.0 + 4.0*rnd01();  // Range 0-4
  p.cont_step_type_increase_chi2dof_required = 0.0 + 4.0*rnd01();  // Range 0-4
  p.skew_nsigma = 0.0 + 10.0*rnd01();  // Range 0-10
  p.left_residual_sum_min_to_try_skew = 0.0 + 5.0*rnd01();  // Range 0-5 (estimated)???
  p.right_residual_sum_min_to_try_skew = 0.0 + 5.0*rnd01();  // Range 0-5 (estimated)???
  p.skew_improve_chi2_dof_threshold = 0.0 + 4.0*rnd01();  // Range 0-4
  p.roi_extent_low_num_fwhm_base = 0.5 + 6.5*rnd01();  // Range 0.5-7
  p.roi_extent_high_num_fwhm_base = 0.5 + 6.5*rnd01();  // Range 0.5-7
  p.roi_extent_mult_type = static_cast<int>(rnd01() < 0.5);  // 0 or 1 for Linear/Sqrt
  p.roi_extent_lower_side_stat_multiple = 0.0 + 1.0*rnd01();  // Range 0-1
  p.roi_extent_upper_side_stat_multiple = 0.0 + 1.0*rnd01();  // Range 0-1
  p.multi_roi_extent_lower_side_fwhm_mult = -1.0 + 3.0*rnd01();  // Range -1 to 2
  p.multi_roi_extent_upper_side_fwhm_mult = -1.0 + 3.0*rnd01();  // Range -1 to 2
}


FinalPeakFitSettings genes_to_settings( const FinalFitSolution &solution )
{
  FinalPeakFitSettings settings;

  settings.combine_nsigma_near = solution.combine_nsigma_near;
  settings.combine_ROI_overlap_frac = solution.combine_ROI_overlap_frac;
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
  settings.roi_extent_low_num_fwhm_base = solution.roi_extent_low_num_fwhm_base;
  settings.roi_extent_high_num_fwhm_base = solution.roi_extent_high_num_fwhm_base;
  settings.roi_extent_mult_type = (solution.roi_extent_mult_type == 0) ?
  FinalPeakFitSettings::RoiExtentMultType::Linear :
  FinalPeakFitSettings::RoiExtentMultType::Sqrt;
  settings.roi_extent_lower_side_stat_multiple = solution.roi_extent_lower_side_stat_multiple;
  settings.roi_extent_upper_side_stat_multiple = solution.roi_extent_upper_side_stat_multiple;
  settings.multi_roi_extent_lower_side_fwhm_mult = solution.multi_roi_extent_lower_side_fwhm_mult;
  settings.multi_roi_extent_upper_side_fwhm_mult = solution.multi_roi_extent_upper_side_fwhm_mult;

  return settings;
}


bool eval_solution( const FinalFitSolution &p, FinalFitCost &c )
{
  const FinalPeakFitSettings settings = genes_to_settings( p );

  c.objective1 = ns_ga_eval_fcn( settings );

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
    const double r = 2.0*(rnd01() - 0.5);
    val += r*mu*(max_val - min_val);
    return std::max(min_val, std::min(max_val, val));
  };

  X_new.combine_nsigma_near = mutate_value(X_base.combine_nsigma_near, 1.0, 10.0);
  X_new.combine_ROI_overlap_frac = mutate_value(X_base.combine_ROI_overlap_frac, -1.0, 1.0);
  X_new.cont_type_peak_nsigma_threshold = mutate_value(X_base.cont_type_peak_nsigma_threshold, 5.0, 60.0);
  X_new.cont_type_left_right_nsigma = mutate_value(X_base.cont_type_left_right_nsigma, 1.0, 20.0);
  //X_new.stepped_roi_extent_lower_side_stat_multiple = mutate_value(X_base.stepped_roi_extent_lower_side_stat_multiple, 0.0, 2.0);
  //X_new.stepped_roi_extent_upper_side_stat_multiple = mutate_value(X_base.stepped_roi_extent_upper_side_stat_multiple, 0.0, 2.0);
  X_new.cont_poly_order_increase_chi2dof_required = mutate_value(X_base.cont_poly_order_increase_chi2dof_required, 0.0, 4.0);
  X_new.cont_step_type_increase_chi2dof_required = mutate_value(X_base.cont_step_type_increase_chi2dof_required, 0.0, 4.0);
  X_new.skew_nsigma = mutate_value(X_base.skew_nsigma, 0.0, 10.0);
  X_new.left_residual_sum_min_to_try_skew = mutate_value(X_base.left_residual_sum_min_to_try_skew, 0.0, 5.0);
  X_new.right_residual_sum_min_to_try_skew = mutate_value(X_base.right_residual_sum_min_to_try_skew, 0.0, 5.0);
  X_new.skew_improve_chi2_dof_threshold = mutate_value(X_base.skew_improve_chi2_dof_threshold, 0.0, 4.0);
  X_new.roi_extent_low_num_fwhm_base = mutate_value(X_base.roi_extent_low_num_fwhm_base, 0.5, 7.0);
  X_new.roi_extent_high_num_fwhm_base = mutate_value(X_base.roi_extent_high_num_fwhm_base, 0.5, 7.0);
  X_new.roi_extent_mult_type = (rnd01() < 0.1) ? (1 - X_base.roi_extent_mult_type) : X_base.roi_extent_mult_type;
  X_new.roi_extent_lower_side_stat_multiple = mutate_value(X_base.roi_extent_lower_side_stat_multiple, 0.0, 1.0);
  X_new.roi_extent_upper_side_stat_multiple = mutate_value(X_base.roi_extent_upper_side_stat_multiple, 0.0, 1.0);
  X_new.multi_roi_extent_lower_side_fwhm_mult = mutate_value(X_base.multi_roi_extent_lower_side_fwhm_mult, -1.0, 2.0);
  X_new.multi_roi_extent_upper_side_fwhm_mult = mutate_value(X_base.multi_roi_extent_upper_side_fwhm_mult, -1.0, 2.0);

  return X_new;
}


FinalFitSolution crossover(
                           const FinalFitSolution& X1,
                           const FinalFitSolution& X2,
                           const std::function<double(void)> &rnd01)
{
  FinalFitSolution X_new;

  auto cross_value = [&rnd01](double v1, double v2) {
    const double r = rnd01();
    return r*v1 + (1.0-r)*v2;
  };

  X_new.combine_nsigma_near = cross_value(X1.combine_nsigma_near, X2.combine_nsigma_near);
  X_new.combine_ROI_overlap_frac = cross_value(X1.combine_ROI_overlap_frac, X2.combine_ROI_overlap_frac);
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
  X_new.roi_extent_low_num_fwhm_base = cross_value(X1.roi_extent_low_num_fwhm_base, X2.roi_extent_low_num_fwhm_base);
  X_new.roi_extent_high_num_fwhm_base = cross_value(X1.roi_extent_high_num_fwhm_base, X2.roi_extent_high_num_fwhm_base);
  X_new.roi_extent_mult_type = (rnd01() < 0.5) ? X1.roi_extent_mult_type : X2.roi_extent_mult_type;
  X_new.roi_extent_lower_side_stat_multiple = cross_value(X1.roi_extent_lower_side_stat_multiple, X2.roi_extent_lower_side_stat_multiple);
  X_new.roi_extent_upper_side_stat_multiple = cross_value(X1.roi_extent_upper_side_stat_multiple, X2.roi_extent_upper_side_stat_multiple);
  X_new.multi_roi_extent_lower_side_fwhm_mult = cross_value(X1.multi_roi_extent_lower_side_fwhm_mult, X2.multi_roi_extent_lower_side_fwhm_mult);
  X_new.multi_roi_extent_upper_side_fwhm_mult = cross_value(X1.multi_roi_extent_upper_side_fwhm_mult, X2.multi_roi_extent_upper_side_fwhm_mult);

  return X_new;
}


double calculate_SO_total_fitness(const GA_Type::thisChromosomeType &X)
{
  // finalize the cost
  double final_cost = 0.0;
  final_cost = X.middle_costs.objective1;
  return final_cost;
}


void SO_report_generation(
                          int generation_number,
                          const EA::GenerationType<FinalFitSolution,FinalFitCost> &last_generation,
                          const FinalFitSolution& best_genes)
{
  cout<<"Generation ["<<generation_number<<"], ";
  cout<<"Best = "<<last_generation.best_total_cost<<", ";
  cout<<"Average = "<<last_generation.average_cost<<", ";

  cout<<", genes = "<<best_genes.to_string(", ")<<endl;
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
  ga_obj.idle_delay_us = 1; // switch between threads quickly
  ga_obj.verbose = false;
  ga_obj.population = 40;
  ga_obj.generation_max = 100;
  ga_obj.calculate_SO_total_fitness = calculate_SO_total_fitness;
  ga_obj.init_genes = init_genes;
  ga_obj.eval_solution = eval_solution;
  ga_obj.mutate = mutate;
  ga_obj.crossover = crossover;
  ga_obj.SO_report_generation = SO_report_generation;
  ga_obj.best_stall_max = 10;
  ga_obj.elite_count = 10;
  ga_obj.crossover_fraction = 0.7;
  ga_obj.mutation_rate = 0.4;

  ns_ga_eval_fcn = ga_eval_fcn;

  ga_obj.solve();

  cout<<"The problem is optimized in "<<timer.toc()<<" seconds."<<endl;

  return genes_to_settings( ga_obj.last_generation.chromosomes[ga_obj.last_generation.best_chromosome_index].genes );
}
};//namespace FinalFit_GA
