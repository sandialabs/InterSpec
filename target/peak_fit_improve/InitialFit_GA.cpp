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

#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp> // For boost::asio::post

#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/DetectorPeakResponse.h"

#include "PeakFitImprove.h"
#include "InitialFit_GA.h"
#include "CandidatePeak_GA.h"

using namespace std;

namespace
{
atomic<size_t> ns_num_evals_this_generation( 0 );

std::function<double( const InitialPeakFindSettings &)> ns_ga_eval_fcn;

std::ofstream sm_output_file;

bool sm_has_been_called = false;
bool sm_set_best_genes = false;
InitialFit_GA::InitialFitSolution sm_best_genes;
double sm_best_total_cost = 1.0E99;

}//namespace


namespace InitialFit_GA
{

vector<PeakDef> initial_peak_find_and_fit( const InitialPeakFindSettings &fit_settings,
                              const FindCandidateSettings &candidate_settings,
                              const shared_ptr<const SpecUtils::Measurement> &data,
                                          const bool multithread,
                              size_t &num_add_candidates_fit_for,  //Only for eval purposes
                              size_t &num_add_candidates_accepted //Only for eval purposes
                              )
{
  assert( RETURN_PeakDef_Candidates ); //dont want to do a static_assert - since I am still compiling all this code both ways...
  if( !RETURN_PeakDef_Candidates )
  {
    cerr << "Please change 'RETURN_PeakDef_Candidates' to true and recompile (needed for evaluating peak-fits)." << endl;
    exit(1);
  }

  num_add_candidates_fit_for = 0;
  num_add_candidates_accepted = 0;

  assert( data );
  assert( data->num_gamma_channels() > 256 );
  const size_t num_channels = data->num_gamma_channels();


  const shared_ptr<const SpecUtils::EnergyCalibration> cal = data->energy_calibration();
  assert( cal && cal->valid() );
  if( !cal || !cal->valid() || !cal->channel_energies() )
    return {};

  const vector<float> &channel_energies = *cal->channel_energies();
  assert( channel_energies.size() == (num_channels + 1) );
  if( channel_energies.size() <= num_channels )
    return {};

  assert( data->gamma_counts() );
  assert( data->gamma_counts()->size() == num_channels );
  if( !data->gamma_counts() || data->gamma_counts()->empty() )
    return {};

  const vector<float> &channel_counts = *data->gamma_counts();

  vector<tuple<float,float,float>> dummy;
  const vector<PeakDef> candidate_peaks = CandidatePeak_GA::find_candidate_peaks( data, 0, 0, dummy, candidate_settings );


#warning "initial_peak_find_and_fit: always assuming HPGe right now - for dev"
  //bool isHPGe = (data->num_gamma_channels() > 5000); // This will be wrong a lot...
  const bool isHPGe = true;
  // TODO: improve selection of HPGe here

  bool amplitudeOnly = false;
  vector<PeakDef> zeroth_fit_results, initial_fit_results;

  if( false )
  {//Begin optimization block - delte when done optimizing `fitPeaks(...)`
    const double start_wall = SpecUtils::get_wall_time();
    const double start_cpu = SpecUtils::get_cpu_time();

    for( size_t i = 0; i < 20; ++i )
    //while( true )
    {
      vector<PeakDef> zeroth_fit_results, initial_fit_results;
#if( USE_LM_PEAK_FIT )
      vector<shared_ptr<const PeakDef>> results_tmp, input_peaks_tmp;
      for( const auto &p : candidate_peaks )
        input_peaks_tmp.push_back( make_shared<PeakDef>(p) );

      PeakFitLM::fit_peaks_LM( results_tmp, input_peaks_tmp, data,
                              fit_settings.initial_stat_threshold, fit_settings.initial_hypothesis_threshold,
                              amplitudeOnly, isHPGe );
      for( const auto &p : results_tmp )
        zeroth_fit_results.push_back( *p );
#else
      fitPeaks( candidate_peaks,
               fit_settings.initial_stat_threshold, fit_settings.initial_hypothesis_threshold,
               data, zeroth_fit_results, amplitudeOnly, isHPGe );
#endif
    }

    const double end_wall = SpecUtils::get_wall_time();
    const double end_cpu = SpecUtils::get_cpu_time();

    cout << "Eval took " << (end_wall - start_wall) << " s, wall and " << (end_cpu - start_cpu) << " s cpu" << endl;
    exit(1);
  }//End optimization block

  /* For some reason calling `fitPeaksInRange(...)` works a lot better than directly calling
   `fitPeaks(...)`, even though it is what `fitPeaksInRange` calls/uses.
   Should probably look into.
   */
  //fitPeaks( candidate_peaks,
  //         fit_settings.initial_stat_threshold, fit_settings.initial_hypothesis_threshold,
  //         data, zeroth_fit_results, dummy_fixedpeaks, amplitudeOnly, isHPGe );
  //const double ncausality = 1.5;
  //zeroth_fit_results = fitPeaksInRange( 0.0, data->gamma_energy_max(), ncausality,
  //                fit_settings.initial_stat_threshold, fit_settings.initial_hypothesis_threshold,
  //                                     candidate_peaks, data, dummy_fixedpeaks, amplitudeOnly, isHPGe );

  auto fit_peaks_per_roi = [&fit_settings, &data, isHPGe, multithread]( const vector<PeakDef> &input_peaks ) -> vector<PeakDef> {

    const bool amplitudeOnly = false;

    map<const PeakContinuum *,vector<PeakDef>> roi_to_peaks_map;
    for( const PeakDef &p : input_peaks )
      roi_to_peaks_map[p.continuum().get()].push_back( p );

    vector<vector<PeakDef>> fit_rois( roi_to_peaks_map.size() );
#if( USE_LM_PEAK_FIT )
    vector<vector<shared_ptr<const PeakDef>>> fit_rois_tmp( roi_to_peaks_map.size() );
#endif

    if( multithread && (PeakFitImprove::sm_num_threads_per_individual > 1) )
    {
      boost::asio::thread_pool pool(PeakFitImprove::sm_num_threads_per_individual);

      size_t fit_rois_index = 0;
      for( auto &roi_peaks : roi_to_peaks_map )
      {
        const vector<PeakDef> &peaks = roi_peaks.second;

#if( USE_LM_PEAK_FIT )
        vector<shared_ptr<const PeakDef>> &results = fit_rois_tmp[fit_rois_index];
        vector<shared_ptr<const PeakDef>> input_peaks_tmp;
        for( const auto &p : peaks )
          input_peaks_tmp.push_back( make_shared<PeakDef>(p) );

        boost::asio::post(pool, boost::bind( &PeakFitLM::fit_peaks_LM,
                                     boost::ref(results),
                                     input_peaks_tmp,
                                     data,
                                     fit_settings.initial_stat_threshold,
                                     fit_settings.initial_hypothesis_threshold,
                                     amplitudeOnly,
                                     isHPGe ) );
#else
        vector<PeakDef> &results = fit_rois[fit_rois_index];

        threadpool.post( boost::bind( &fitPeaks,
                                     boost::cref(peaks),
                                     fit_settings.initial_stat_threshold,
                                     fit_settings.initial_hypothesis_threshold,
                                     data,
                                     boost::ref( results ),
                                     amplitudeOnly,
                                     isHPGe ) );
#endif
        fit_rois_index += 1;
      }//for( auto &roi_peaks : roi_to_peaks_map )

      pool.join();
    }else
    {
      size_t fit_rois_index = 0;
      for( auto &roi_peaks : roi_to_peaks_map )
      {
        const vector<PeakDef> &peaks = roi_peaks.second;

#if( USE_LM_PEAK_FIT )
        vector<shared_ptr<const PeakDef>> &results = fit_rois_tmp[fit_rois_index];
        vector<shared_ptr<const PeakDef>> input_peaks_tmp;
        for( const auto &p : peaks )
          input_peaks_tmp.push_back( make_shared<PeakDef>(p) );

        PeakFitLM::fit_peaks_LM( results, input_peaks_tmp, data,
                                fit_settings.initial_stat_threshold, fit_settings.initial_hypothesis_threshold,
                                amplitudeOnly, isHPGe );
#else
        vector<PeakDef> &results = fit_rois[fit_rois_index];
        fitPeaks( peaks, fit_settings.initial_stat_threshold, fit_settings.initial_hypothesis_threshold,
                 data, results, amplitudeOnly, isHPGe );
#endif
        fit_rois_index += 1;
      }//for( auto &roi_peaks : roi_to_peaks_map )
    }//if( multithread ) / else

#if( USE_LM_PEAK_FIT )
    for( size_t i = 0; i < fit_rois.size(); ++i )
    {
      vector<PeakDef> &to_fil = fit_rois[i];
      assert( to_fil.empty() );
      const vector<shared_ptr<const PeakDef>> &to_fill_from = fit_rois_tmp[i];
      for( const auto &p : to_fill_from )
        to_fil.push_back( *p );
    }
#endif

    vector<PeakDef> fit_results;
    fit_results.reserve( input_peaks.size() );

    for( const vector<PeakDef> &fit_peaks : fit_rois )
      fit_results.insert( end(fit_results), begin(fit_peaks), end(fit_peaks) );

    std::sort( begin(fit_results), end(fit_results), &PeakDef::lessThanByMean );

    return fit_results;
  };//fit_peaks_per_roi lambda

  zeroth_fit_results = fit_peaks_per_roi( candidate_peaks );

  // Now go through and make sure there is at least `0.5*initial_min_nsigma_roi` on left and right,
  //  and no more than `0.5*initial_max_nsigma_roi`, since the mean may have changed, etc
  bool modified_roi = false;
  for( size_t i = 0; i < zeroth_fit_results.size(); ++i )
  {
    PeakDef &peak = zeroth_fit_results[i];

    // We wont mess with multi peak ROIs right now
    if( (i > 0) && (peak.continuum() == zeroth_fit_results[i].continuum()) )
      continue;

    if( ((i + 1) < zeroth_fit_results.size()) && (peak.continuum() == zeroth_fit_results[i+1].continuum()) )
      continue;

    const double min_around = 0.5*fit_settings.initial_min_nsigma_roi * peak.sigma();
    const double max_around = 0.5*fit_settings.initial_max_nsigma_roi * peak.sigma();

    if( (peak.mean() - min_around) < peak.continuum()->lowerEnergy() )
    {
      modified_roi = true;
      peak.continuum()->setRange(peak.mean() - min_around, peak.continuum()->upperEnergy() );
    }

    if( (peak.mean() + min_around) > peak.continuum()->upperEnergy() )
    {
      modified_roi = true;
      peak.continuum()->setRange( peak.continuum()->lowerEnergy(), peak.mean() + min_around  );
    }


    if( (peak.mean() - max_around) > peak.continuum()->lowerEnergy() )
    {
      modified_roi = true;
      peak.continuum()->setRange(peak.mean() - max_around, peak.continuum()->upperEnergy() );
    }

    if( (peak.mean() + max_around) < peak.continuum()->upperEnergy() )
    {
      modified_roi = true;
      peak.continuum()->setRange( peak.continuum()->lowerEnergy(), peak.mean() + max_around  );
    }
  }//for( PeakDef &peak : zeroth_fit_results )

  //fitPeaks( zeroth_fit_results,
  //         fit_settings.initial_stat_threshold, fit_settings.initial_hypothesis_threshold,
  //         data, initial_fit_results, dummy_fixedpeaks, amplitudeOnly, isHPGe );

  //initial_fit_results = fitPeaksInRange( 0.0, data->gamma_energy_max(), ncausality,
  //                fit_settings.initial_stat_threshold, fit_settings.initial_hypothesis_threshold,
  //                zeroth_fit_results, data, dummy_fixedpeaks, amplitudeOnly, isHPGe );
  initial_fit_results = fit_peaks_per_roi( zeroth_fit_results );

  //vector<PeakDef> refit_results;
  //fitPeaks( initial_fit_results,
  //         fit_settings.initial_stat_threshold, fit_settings.initial_hypothesis_threshold,
  //         data, refit_results, dummy_fixedpeaks, amplitudeOnly, isHPGe );
  //refit_results.swap( initial_fit_results );

  //for( const auto &p : initial_fit_results )
  //{
  //  const bool significant = chi2_significance_test( p, fit_settings.initial_stat_threshold, fit_settings.initial_hypothesis_threshold,
  //                                                  {}, data );
  //  cout << "Peak at " << p.mean() << " significant=" << significant << " with nsigma=" << p.peakArea()/p.peakAreaUncert() << endl;
  //}


  // We will try to estimate the FWHM functional form to help search for peaks that weren't
  //  found by the candidate peak algorithm

  // Lets decide the FWHM equation form - we'll try to use what is requested, but if we dont
  //  have enough peaks, we'll downgrade things.
  int fwhm_poly_eqn_order = 0;
  vector<float> fwhm_coeffs, fwhm_coeff_uncerts;
  DetectorPeakResponse::ResolutionFnctForm fwhm_eqn_form = DetectorPeakResponse::ResolutionFnctForm::kNumResolutionFnctForm;
  switch( fit_settings.fwhm_fcn_form )
  {
    case InitialPeakFindSettings::FwhmFcnForm::Gadras:
      fwhm_eqn_form = DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn;
      break;

    case InitialPeakFindSettings::FwhmFcnForm::SqrtPolynomialTwoCoefs:
      if( initial_fit_results.size() > 2 )
      {
        fwhm_poly_eqn_order = 2;
        fwhm_eqn_form = DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;
      }
      break;

    case InitialPeakFindSettings::FwhmFcnForm::SqrtPolynomialThreeCoefs:
      fwhm_eqn_form = DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;
      if( initial_fit_results.size() > 3 )
      {
        fwhm_poly_eqn_order = 3;
        fwhm_eqn_form = DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;
      }else if( initial_fit_results.size() > 2 )
      {
        fwhm_poly_eqn_order = 2;
        fwhm_eqn_form = DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;
      }
      break;

    case InitialPeakFindSettings::FwhmFcnForm::SqrtEnergyPlusInverse:
      if( initial_fit_results.size() > 4 )
      {
        fwhm_eqn_form = DetectorPeakResponse::ResolutionFnctForm::kSqrtEnergyPlusInverse;
      }else if( initial_fit_results.size() > 2 )
      {
        fwhm_poly_eqn_order = 2;
        fwhm_eqn_form = DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;
      }
      break;

    case InitialPeakFindSettings::FwhmFcnForm::NumFwhmFcnForm:
      assert( 0 );
      break;
  }//switch( fit_settings.fwhm_fcn_form )


  try
  {
    if( initial_fit_results.empty() )
      throw runtime_error( "No peaks" );

    auto peaks = make_shared< deque<shared_ptr<const PeakDef>> >();
    for( const auto &p : initial_fit_results )
      peaks->push_back( make_shared<PeakDef>(p) );


    const PeakFitUtils::CoarseResolutionType coarse_res_type
                                  = PeakFitUtils::coarse_resolution_from_peaks( *peaks );
    if( coarse_res_type != PeakFitUtils::CoarseResolutionType::High )
    {
#ifndef NDEBUG
      cerr << "Found coarse_res_type=" << static_cast<int>(coarse_res_type) << " with " << initial_fit_results.size() << " peaks." << endl;
      cerr << endl;
#endif
    }

#warning "initial_peak_find_and_fit: always assuming HPGe right now - for dev"
    //const bool isHPGe = (coarse_res_type == PeakFitUtils::CoarseResolutionType::High)

    DetectorPeakResponse::ResolutionFnctForm form = DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;
    if( initial_fit_results.size() <= 2 )
      fwhm_eqn_form = DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn;

    MakeDrfFit::performResolutionFit( peaks, fwhm_eqn_form, fwhm_poly_eqn_order,
                                     fwhm_coeffs, fwhm_coeff_uncerts );

    double cs137_fwhm = DetectorPeakResponse::peakResolutionFWHM( 661, fwhm_eqn_form, fwhm_coeffs );
    //cout << "Initial FWHM=" << 100*cs137_fwhm/661 << "%" << endl;
    if( IsNan(cs137_fwhm) && (initial_fit_results.size() > 1) )
    {
      // This can happen if we have 2 or 3 peaks, below 661 keV, so that by the time we
      //  get up to 661, the equation is invalid
      peaks->erase( begin(*peaks) );

      if( initial_fit_results.size() <= 2 )
        fwhm_eqn_form = DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn;

      MakeDrfFit::performResolutionFit( peaks, fwhm_eqn_form, fwhm_poly_eqn_order,
                                       fwhm_coeffs, fwhm_coeff_uncerts );

      cs137_fwhm = DetectorPeakResponse::peakResolutionFWHM( 661, form, fwhm_coeffs );
    }//if( IsNan(cs137_fwhm) && (initial_fit_results.size() > 1) )

    if( IsNan(cs137_fwhm) || ((100*cs137_fwhm/661) < 0.05) )
      throw std::runtime_error( "Failed to fit valid FWHM function." );

    // Go through and remove outliers....
    if( peaks->size() > 5 )
    {
      auto new_peaks = make_shared< deque<shared_ptr<const PeakDef>> >();
      vector<pair<double,shared_ptr<const PeakDef>>> distances;
      for( const auto &p : *peaks )
      {
        double pred_fwhm = 0.0;
        const double mean_MeV = 0.001 * p->mean();  //kSqrtPolynomial
        for( size_t i = 0; i < fwhm_coeffs.size(); ++i )
          pred_fwhm += fwhm_coeffs[i] * std::pow( mean_MeV, 1.0*i );
        pred_fwhm = sqrt( pred_fwhm );
        //pred_fwhm = fwhm_coeffs[0] + fwhm_coeffs[1]*sqrt( p->mean() ); //

        const double frac_diff = fabs( p->fwhm() - pred_fwhm ) / p->fwhm();
        if( !IsNan(frac_diff) && !IsInf(frac_diff) )
          distances.emplace_back( frac_diff, p );
      }//for( const auto &p : initial_fit_peaks )

      std::sort( begin(distances), end(distances),
                []( const pair<double,shared_ptr<const PeakDef>> &lhs,
                   const pair<double,shared_ptr<const PeakDef>> &rhs) -> bool {
        return lhs.first > rhs.first;
      } );

      // Limit to un-selecting max of 20% of peaks (arbitrarily chosen), if the deviate
      // more than 25% from the fit (again, arbitrarily chosen).
      const size_t max_remove = static_cast<size_t>( std::ceil( 0.2*distances.size() ) );
      for( size_t index = 0; index < distances.size(); ++index )
      {
        if( (distances[index].first < 0.25) || (index > max_remove) )
          new_peaks->push_back( distances[index].second );
      }

      MakeDrfFit::performResolutionFit( new_peaks, fwhm_eqn_form, fwhm_poly_eqn_order,
                                       fwhm_coeffs, fwhm_coeff_uncerts );
    }//if( peaks->size() > 5 )

    cs137_fwhm = DetectorPeakResponse::peakResolutionFWHM( 661, form, fwhm_coeffs );

    if( IsNan(cs137_fwhm) || ((100*cs137_fwhm/661) < 0.05) )
      throw std::runtime_error( "Failed to fit valid FWHM function." );

    //cout << "Final FWHM=" << 100*cs137_fwhm/661 << "%" << endl;
  }catch( std::exception &e )
  {
    if( PeakFitImprove::debug_printout )
      cerr << ("Caught exception: " + string(e.what()) + "\n");

    fwhm_coeffs.clear();
    fwhm_coeff_uncerts.clear();
    fwhm_poly_eqn_order = 0;
    fwhm_eqn_form = DetectorPeakResponse::ResolutionFnctForm::kNumResolutionFnctForm;
  }//try / catch to fit a FWHM functional form


  // Now go through and find peaks we may have missed.
  //  We have the following three variables to help us with this:
  //    - fit_settings.search_roi_nsigma_deficit
  //    - fit_settings.search_stat_threshold
  //    - fit_settings.search_hypothesis_threshold
  //    - fit_settings.search_stat_significance
  vector<PeakDef> peaks_added_for_debug;
  try
  {
    auto ecal = data->energy_calibration();
    assert( ecal && ecal->valid() );

    const vector<float> &counts = *data->gamma_counts();


    // We dont want to identify the initial "turn-on" of data as a peak, so, fairly arbitrarily,
    //  we'll fast-forward until we see three slightly statistically significant decreases in
    //  channel counts.  This has not been optimized, and just sanity checked on a few spectra.
    size_t start_channel = 3;
    for( size_t num_down_chnls = 0; (num_down_chnls < 3) && start_channel < (num_channels - 3); ++start_channel )
    {
      assert( start_channel > 1 );
      const float &prev_contents = counts[start_channel - 1];
      const float &this_contents = counts[start_channel];
      const float nsigma_decrease_required = 0.75f; //Choosen fairly arbitrarily
      num_down_chnls += ((prev_contents > 0.0)
                          && (this_contents < (prev_contents - nsigma_decrease_required* sqrt(prev_contents))));
    }

    for( size_t center_channel = start_channel; center_channel < (num_channels - 3); ++center_channel )
    {
      // We will only check this channel if `center_channel` is higher than the three channels
      //  on either side
      bool highest_around = true;
      for( int i = -3; highest_around && (i <= 3); ++i )
        highest_around = ((i == 0) || (counts[center_channel] > counts[center_channel + i]));
      if( !highest_around )
        continue;

      float fwhm = std::numeric_limits<float>::quiet_NaN();
      const float energy = data->gamma_channel_center( center_channel );

      if( fwhm_eqn_form != DetectorPeakResponse::ResolutionFnctForm::kNumResolutionFnctForm )
        fwhm = DetectorPeakResponse::peakResolutionFWHM( energy, fwhm_eqn_form, fwhm_coeffs );

      if( (fwhm < 0.0001) || IsNan(fwhm) )
      {
        if( isHPGe )
          fwhm = PeakFitUtils::hpge_fwhm_fcn( energy );
        else
          fwhm = PeakFitUtils::nai_fwhm_fcn( energy );
      }//if( we dont have a valid FWHM )

      // `roi_lchannel` and `roi_uchannel` are the first and last channels of the ROI - the
      //  ROI includes these channels.
      //  Note: kinda arbitrarily using a ROI width of 4 FWHM - this has not been optimized beyond
      //        noting for a few ROIs, using only 2 FWHM lead to pretty poor results.
      //  TODO: after optimization of final fit parameters, should use those to decide ROI range, and re-optimize this stage
      int roi_lchannel = static_cast<int>( round( ecal->channel_for_energy(energy - 2*fwhm) ) );
      int roi_uchannel = static_cast<int>( round( ecal->channel_for_energy(energy + 2*fwhm) ) );

      // Make sure the ROI includes at least two channels both above and below the center channel.
      if( (roi_lchannel+1) >= static_cast<int>(center_channel) )
        roi_lchannel = static_cast<int>(center_channel) - 2;

      if( (roi_uchannel-1) <= static_cast<int>(center_channel) )
        roi_uchannel = static_cast<int>(center_channel) + 2;

      if( (roi_lchannel <= 5) || ((roi_uchannel + 5) >= static_cast<int>(num_channels)) )
        continue;

      const float roi_lower_energy = data->gamma_channel_lower( static_cast<size_t>(roi_lchannel) );
      const float roi_upper_energy = data->gamma_channel_upper( static_cast<size_t>(roi_uchannel) );

      // Check to make sure we arent overlapping with any existing ROIs.
      //  But allow some overlap if the peaks would be greater than two FWHM from each other
      double nearest_left_peak_upper = -1.0E6, nearest_right_peak_lower = 1.0E6;

      bool overlaps = false;
      for( size_t i = 0; !overlaps && (i < initial_fit_results.size()); ++i )
      {
        const PeakDef &p = initial_fit_results[i];
        const double peak_lower = std::max( p.lowerX(), p.mean() - p.fwhm() );
        const double peak_upper = std::min( p.upperX(), p.mean() + p.fwhm() );

        if( peak_upper < energy )
          nearest_left_peak_upper = std::max( nearest_left_peak_upper, peak_upper );
        if( peak_lower > energy )
          nearest_right_peak_lower = std::min( nearest_right_peak_lower, peak_lower );

        overlaps = ((roi_lower_energy <= p.upperX()) && (roi_upper_energy >= p.lowerX()))
                    && (fabs(energy - p.mean()) < 2*p.fwhm());
      }

      if( overlaps )
        continue;

      //if( energy > 967 && energy < 972 )
      //  cout << "For trial channel at energy=" << energy
      //  << " [" << roi_lower_energy << ", " << roi_upper_energy << "],"
      //  << " nearest_left_peak_upper=" << nearest_left_peak_upper
      //  << ", nearest_right_peak_lower=" << nearest_right_peak_lower
      //  << endl;

      // We'll fairly arbitrarily choose number of side channels to use to fit the continuum from
      int num_side_channels = (roi_uchannel - roi_lchannel) / 4;
      num_side_channels = std::max(num_side_channels, 4);
      num_side_channels = std::min(num_side_channels, 16);

      if( roi_lchannel <= num_side_channels )
        continue;

      if( (roi_uchannel + num_side_channels) >= static_cast<int>(num_channels) )
        continue;

      int nlower_channel = num_side_channels;
      const float side_lower_energy = data->gamma_channel_lower(roi_lchannel - nlower_channel);
      if( side_lower_energy < nearest_left_peak_upper )
      {
        const size_t prev_roi_uchannel = data->find_gamma_channel( nearest_left_peak_upper );
        if( static_cast<int>(prev_roi_uchannel) >= (roi_lchannel - 2) )
        {
          // Fairly arbitrary cutting down of ROI
          if( (roi_uchannel - roi_lchannel) > 3 )
            roi_lchannel += 1;
          nlower_channel = 2;
        }else
        {
          assert( roi_lchannel > prev_roi_uchannel );
          nlower_channel = roi_lchannel - static_cast<int>(prev_roi_uchannel);
          nlower_channel = std::max(nlower_channel, 2); //isnt needed, but jic
        }
      }//if( side_lower_energy < nearest_left_peak_upper )

      int nupper_channel = num_side_channels;
      const float side_upper_energy = data->gamma_channel_upper(roi_uchannel + nupper_channel);
      if( side_upper_energy > nearest_right_peak_lower )
      {
        const size_t prev_roi_lchannel = data->find_gamma_channel( nearest_right_peak_lower );
        if( prev_roi_lchannel < (roi_uchannel + 2) )
        {
          // Fairly arbitrary cutting down of ROI
          if( (roi_uchannel - roi_lchannel) > 3 )
            roi_uchannel -= 1;
          nupper_channel = 2;
        }else
        {
          assert( prev_roi_lchannel > roi_uchannel );
          nupper_channel = static_cast<int>(prev_roi_lchannel) - roi_uchannel;
          nupper_channel = std::max(nupper_channel, 2); //isnt needed, but jic
        }
      }//if( side_upper_energy > nearest_right_peak_lower )


      vector<double> cont_equation{0.0, 0.0};
      PeakContinuum::eqn_from_offsets( roi_lchannel, roi_uchannel, energy, data,
                                      nlower_channel, nupper_channel,
                                      cont_equation[1], cont_equation[0] );

      const double cont_area = PeakContinuum::offset_eqn_integral( &(cont_equation[0]),
                                                PeakContinuum::OffsetType::Linear,
                                                roi_lower_energy, roi_upper_energy, energy );

      const double data_area = data->gamma_integral( roi_lower_energy, roi_upper_energy );
      if( data_area <= 0.0 )
        continue;

      const double data_uncert = sqrt( data_area );

      const double area_diff = data_area - cont_area;

      if( area_diff < fit_settings.search_roi_nsigma_deficit*data_uncert )
        continue;

      num_add_candidates_fit_for += 1;

      PeakDef peak( energy, fwhm/PhysicalUnits::fwhm_nsigma, area_diff );
      peak.continuum()->setType( PeakContinuum::OffsetType::Linear );
      peak.continuum()->setParameters( energy, cont_equation, {} );
      peak.continuum()->setRange( roi_lower_energy, roi_upper_energy );

      vector<PeakDef> trial_peak_fit_results;
      //fitPeaks( {peak},
      //         fit_settings.search_stat_threshold, fit_settings.search_hypothesis_threshold,
      //         data, trial_peak_fit_results, dummy_fixedpeaks, amplitudeOnly, isHPGe );

      //amplitudeOnly = false;
      //trial_peak_fit_results = fitPeaksInRange( peak.mean() - 1.0, peak.mean() + 1.0, ncausality,
      //                fit_settings.search_stat_threshold, fit_settings.search_hypothesis_threshold,
      //                                         {peak}, data, dummy_fixedpeaks, amplitudeOnly, isHPGe );
      trial_peak_fit_results = fit_peaks_per_roi( {peak} );


      bool is_stat_sig = false;
      for( const auto &p : trial_peak_fit_results )
      {
        const double area = p.peakArea(), uncert = p.peakAreaUncert();
        is_stat_sig |= ((area > 0.0) && (uncert > 0.0) && ((area/uncert) > fit_settings.search_stat_significance));
      }

      if( is_stat_sig && trial_peak_fit_results.size() )
      {
        assert( trial_peak_fit_results.size() == 1 );
        num_add_candidates_accepted += 1;

        // Now fast forward to edge of new ROI
        assert( roi_uchannel > center_channel );
        center_channel = std::max( center_channel, static_cast<size_t>(roi_uchannel) );

        // We will add this peak to `initial_fit_results` right now, so that way the next potential
        //  peak can take its ROI into account.
        initial_fit_results.push_back( trial_peak_fit_results.front() );
        std::sort( begin(initial_fit_results), end(initial_fit_results), &PeakDef::lessThanByMean );

        peaks_added_for_debug.push_back( trial_peak_fit_results.front() );
      }
    }//for( loop over center_channel )

    if( peaks_added_for_debug.size() )
    {
    }
  }catch( std::exception &e )
  {
    // I dont think this should ever really happen - but we'll check.
    cerr << "initial_peak_find_and_fit: caught exception searching for more candidate peaks: " << e.what() << endl;
    assert( 0 );
  }


  // Now lets try to add peaks in existing ROIs
  //  for ~+-2.5 FWHM, look for a bin to go statistically significantly down, and if so, try adding a peak to the left/right of
  //  that, with the mean at the highest bin within the next ~1 FWHM.
  //
  // Lets look on the left and right side of each ROI to see if adding a peak there could be
  //  helpful.  Smaller features within ~`num_smooth_side_channels` of larger features will
  //  tend to get bulldozed, so we'll take this as our characteristic distance, and go out 1.5
  //  times (arbitrarily chosen) this length
  const size_t num_extra_chnl = (3*candidate_settings.num_smooth_side_channels) / 2;

  bool added_any_peaks = false;
  vector<shared_ptr<PeakDef>> peaks_with_add;
  for( const PeakDef &p : initial_fit_results )
    peaks_with_add.push_back( make_shared<PeakDef>(p) );

  vector< pair<shared_ptr<const PeakContinuum>,vector<shared_ptr<const PeakDef>>> > roi_to_peaks;

  {// begin be lazy and use a map to find ROIs to put in `roi_to_peaks`
    map<shared_ptr<const PeakContinuum>,vector<shared_ptr<const PeakDef>>> roi_to_peaks_map;
    for( shared_ptr<PeakDef> &p : peaks_with_add )
      roi_to_peaks_map[p->continuum()].push_back( p );
    for( const auto &c_p : roi_to_peaks_map )
      roi_to_peaks.push_back( c_p );

    // We dont really need to sort things, but we will since its nicer for debugging
    std::sort( begin(roi_to_peaks), end(roi_to_peaks),
        []( const pair<shared_ptr<const PeakContinuum>,vector<shared_ptr<const PeakDef>>> &lhs,
            const pair<shared_ptr<const PeakContinuum>,vector<shared_ptr<const PeakDef>>> &rhs ) -> bool {
      return lhs.first->referenceEnergy() < rhs.first->referenceEnergy();
    } );
  }// end be lazy and use a map to find ROIs to put in `roi_to_peaks`

  // When we add a peak, we will decrement the loop index, so it will try to add another peak
  //  to the same ROI - but to avoid some crazy cycle, we will track this.
  map<int,int> num_iters;
  for( int cont_peak_index = 0; cont_peak_index < static_cast<int>(roi_to_peaks.size()); ++cont_peak_index )
  {
    // I'm like 99.99% sure we can simply do `num_iters[cont_peak_index] += 1` and the counter will be initialized to zero - but I'll be extra sure
    auto num_iter_pos = num_iters.find(cont_peak_index);
    if( num_iter_pos == num_iters.end() )
      num_iter_pos = num_iters.insert( make_pair(cont_peak_index, 0) ).first;
    num_iter_pos->second += 1;
    if( num_iter_pos->second > 5 ) //arbitrarily allow adding up to 5 peaks
      continue;

    auto &c_p = roi_to_peaks[cont_peak_index];
    const shared_ptr<const PeakContinuum> orig_continuum = c_p.first;
    const vector<shared_ptr<const PeakDef>> &peaks_in_roi = c_p.second;

    assert( !peaks_in_roi.empty() );

    const double old_lower_energy = orig_continuum->lowerEnergy();
    const double old_upper_energy = orig_continuum->upperEnergy();

    size_t lowest_mean_channel = data->num_gamma_channels(), highest_mean_channel = 0;
    for( const auto &p : peaks_in_roi )
    {
      const double mean = p->mean();
      const double sigma = p->sigma();
      const size_t mean_channel = data->find_gamma_channel(mean);
      lowest_mean_channel = std::min( lowest_mean_channel, mean_channel );
      highest_mean_channel = std::max( highest_mean_channel, mean_channel );
    }

    const size_t old_roi_lower_chnl = data->find_gamma_channel( old_lower_energy );
    const size_t old_roi_upper_chnl = data->find_gamma_channel( old_upper_energy );


    // We'll fairly arbitrarily choose number of side channels to use to fit the continuum from
    size_t num_side_channels = (old_roi_upper_chnl - old_roi_lower_chnl) / 4;
    num_side_channels = std::max( num_side_channels, size_t(4) );
    num_side_channels = std::min( num_side_channels, size_t(16) );

    if( ((num_extra_chnl + num_side_channels + 1) >= old_roi_lower_chnl)
       || ((old_roi_upper_chnl + num_extra_chnl + num_side_channels + 2) >= num_channels) )
    {
      continue;
    }

    const size_t start_chnl = old_roi_lower_chnl - num_extra_chnl;
    const size_t end_chnl = old_roi_upper_chnl + num_extra_chnl;


    double avrg_gaus_sigma = 0.0, max_significance = 0.0;
    vector<double> means, sigmas;
    for( const shared_ptr<const PeakDef> &p : peaks_in_roi )
    {
      avrg_gaus_sigma += p->sigma();
      means.push_back( p->mean() );
      sigmas.push_back( p->sigma() );

      assert( (p->peakArea() < 0.0) || (p->peakAreaUncert() > 0.0) );
      if( p->peakAreaUncert() > 0.0 )
        max_significance = std::max( max_significance, p->peakArea() / p->peakAreaUncert() );
    }

    avrg_gaus_sigma /= peaks_in_roi.size();


    // Variable to track which channel is best to have new peak mean at, by using statistical
    //  significance of new peak.
    double max_stat_sig = 0.0;

    // Some variables just for debug printout
    double max_stat_chi2_improve = 0.0, max_stat_chi2_dof_improve = 0.0;

    // The set of peaks, where the added peak is most significant
    vector<shared_ptr<PeakDef>> max_sig_peaks;

    // Find other nearby peaks - we'll just kinda WAG it, since it doesnt matter a ton, I dont think
    vector<PeakDef> nearbyOtherRoiPeaks;
    for( const auto &other_c_p : roi_to_peaks )
    {
      const shared_ptr<const PeakContinuum> &other_cont = other_c_p.first;

      if( other_cont == orig_continuum )
        continue;

      const double other_lower_roi = other_cont->lowerEnergy() - 7*avrg_gaus_sigma;
      const double other_upper_roi = other_cont->lowerEnergy() + 7*avrg_gaus_sigma;

      if( (other_lower_roi < old_upper_energy) && (other_upper_roi > old_lower_energy) )
      {
        for( const auto p : other_c_p.second )
          nearbyOtherRoiPeaks.push_back( *p );
      }
    }//for( const auto &other_c_p : roi_to_peaks )


    for( size_t chnl = start_chnl; chnl < end_chnl; ++chnl )
    {
      const float added_mean = data->gamma_channel_center( chnl );

      // If we are closer than ~1.5 gaussian sigma widths to an existing peak in the ROI, I dont
      //  think we can fairly fit both peaks.  But we'll be a bit conservative and use 2 sigma.
      bool to_close = false;
      for( size_t i = 0; i < means.size(); ++i )
      {
        const double nsigma_indistinguishable = 2;
        to_close |= ((fabs(added_mean - means[i]) / sigmas[i]) < nsigma_indistinguishable);
      }

      if( to_close )
        continue;

      // We can probably use some simple logic to significantly reduce CPU time, but we'll
      //  leave this for later

      // We'll make sure the ROI goes at least 1.5 (chosen arbitrarily) FWHM to each side of added peak
      const double min_fwhm_roi_extend = 1.5;
      size_t trial_start_chnl = std::min( old_roi_lower_chnl, data->find_gamma_channel(added_mean - min_fwhm_roi_extend*avrg_gaus_sigma*PhysicalUnits::fwhm_nsigma) );
      size_t trial_end_chnl = std::max( old_roi_upper_chnl, data->find_gamma_channel(added_mean + min_fwhm_roi_extend*avrg_gaus_sigma*PhysicalUnits::fwhm_nsigma) );

      assert( trial_start_chnl < num_channels || trial_end_chnl < num_channels );
      if( trial_start_chnl >= num_channels || trial_end_chnl >= num_channels )
        continue;

      assert( trial_end_chnl > trial_start_chnl );
      if( (trial_start_chnl + 4) > trial_end_chnl )
        continue;


      // Lets check the other ROIs and see if we are overlapping with them
      //double edge_ish = (chnl < lowest_mean_channel) ? added_mean - 2.0*avrg_gaus_sigma : added_mean + 2.0*avrg_gaus_sigma;

      bool overlaps = false;
      for( const auto &other_c_p : roi_to_peaks )
      {
        const shared_ptr<const PeakContinuum> &other_cont = other_c_p.first;
        if( other_cont == orig_continuum )
          continue;

        //if( (edge_ish > other_cont->lowerEnergy())
        //    && (edge_ish < other_cont->upperEnergy()) )
        if( (added_mean > other_cont->lowerEnergy())
           && (added_mean < other_cont->upperEnergy()) )
        {
          overlaps = true;
          break;
        }
      }//for( const auto &other_c_p : roi_to_peaks )

      if( overlaps )
      {
        if( chnl > lowest_mean_channel )
          break;
        continue;
      }

      const size_t num_trial_channel = (1 + trial_end_chnl - trial_start_chnl);

      try
      {
        assert( orig_continuum->isPolynomial() );

        const int num_prev_poly = orig_continuum->isPolynomial() ? static_cast<int>(orig_continuum->parameters().size()) : 2;

        const float * const energies = &(channel_energies[trial_start_chnl]);
        const float * const data_cnts = &(channel_counts[trial_start_chnl]);
        int num_polynomial_terms = std::max( num_prev_poly, 2 );

        //Really high stat HPGe peaks that have a "step" in them are susceptible to getting a bunch
        //  of peaks added on in the low-energy range.  So for these we'll make the continuum
        //  step-functions - either flat, or if higher-stat, linear.
        //  These threshold are pure guesses - but its maybe not worth the effort to include in the optimization?
        const double significance_for_flat_step_continuum = 40;
        const double significance_for_linear_step_continuum = 60;

        bool step_continuum = (PeakContinuum::is_step_continuum( orig_continuum->type() )
                                    || (max_significance > significance_for_flat_step_continuum));
        if( significance_for_linear_step_continuum > 50 )
          num_polynomial_terms = std::max( num_polynomial_terms, 3 );

        const double ref_energy = orig_continuum->referenceEnergy();
        vector<double> trial_means = means, trial_sigmas = sigmas;


        const double *skew_parameters = nullptr;

        vector<double> amplitudes, continuum_coeffs, amplitudes_uncerts, continuum_coeffs_uncerts;

        const double without_new_peak_chi2 = fit_amp_and_offset( energies, data_cnts,
                              num_trial_channel, num_polynomial_terms, step_continuum, ref_energy,
                              trial_means, trial_sigmas, nearbyOtherRoiPeaks, PeakDef::SkewType::NoSkew,
                              skew_parameters, amplitudes, continuum_coeffs,
                              amplitudes_uncerts, continuum_coeffs_uncerts );

        const double initial_peak_area = std::accumulate( begin(amplitudes), end(amplitudes), 0.0 );

        trial_means.push_back( added_mean );
        trial_sigmas.push_back( avrg_gaus_sigma );
        const double with_new_peak_chi2 = fit_amp_and_offset( energies, data_cnts,
                              num_trial_channel, num_polynomial_terms, step_continuum, ref_energy,
                              trial_means, trial_sigmas, nearbyOtherRoiPeaks, PeakDef::SkewType::NoSkew,
                              skew_parameters, amplitudes, continuum_coeffs,
                              amplitudes_uncerts, continuum_coeffs_uncerts );

        assert( (means.size() + 1) == amplitudes.size() );
        if( (means.size() + 1) != amplitudes.size() )
          throw runtime_error( "Unexpectedly didnt fit peak...wtf?" );

        const double after_peak_area = std::accumulate( begin(amplitudes), end(amplitudes), 0.0 );

        // Lets not let total peak area decrease
        if( after_peak_area < initial_peak_area )
          continue;

        // If any peak in the ROI is insignificant, we'll reject this hypothesis
        bool a_peak_is_insignificant = false;
        for( size_t existing_index = 0; existing_index < amplitudes.size(); ++existing_index )
        {
          const double peak_area = amplitudes[existing_index];
          const double peak_area_uncert = amplitudes_uncerts[existing_index];

          assert( peak_area_uncert >= 0.0 );

          if( (peak_area < 0.0) || (peak_area_uncert < 0.0)
             || (peak_area < fit_settings.ROI_add_nsigma_required*peak_area_uncert) )
          {
            a_peak_is_insignificant = true;
            break;
          }
        }//for( loop over existing_index )

        if( a_peak_is_insignificant )
          continue;

        const double new_peak_area = amplitudes.back();
        const double new_peak_uncert = amplitudes_uncerts.back();
        const double nsigma = new_peak_area / new_peak_uncert;

        const size_t ndof = (num_trial_channel - amplitudes.size()) - num_polynomial_terms;
        const double chi2_improve = (without_new_peak_chi2 - with_new_peak_chi2);
        const double chi2_dof_improve = chi2_improve / ndof;
        // Note/TODO: The Chi2/Dof wont improve as much when adding a peak to a wide ROI - maybe there is a better measure to use?

        if( (nsigma > max_stat_sig)
           && (nsigma > fit_settings.ROI_add_nsigma_required)
           && (chi2_dof_improve > fit_settings.ROI_add_chi2dof_improve)
           && ((nsigma < significance_for_linear_step_continuum) || (chi2_dof_improve > 2*fit_settings.ROI_add_chi2dof_improve) ) //This is totally ad-hoc, and not tested
           )
        {
          max_stat_sig = nsigma;

          max_stat_chi2_improve = chi2_improve;
          max_stat_chi2_dof_improve = chi2_dof_improve;


          // Its wasteful to construct the peaks here/now, but it saves some code, so we'll do it anyway.
          max_sig_peaks.clear();
          for( size_t i = 0; i < amplitudes.size(); ++i )
          {
            auto new_peak = make_shared<PeakDef>(trial_means[i], trial_sigmas[i], amplitudes[i]);
            new_peak->setAmplitudeUncert( amplitudes_uncerts[i] ); // We dont know sigma and mean uncertainties

            max_sig_peaks.push_back( new_peak );
            if( i )
              new_peak->setContinuum( max_sig_peaks.front()->continuum() );
          }//for( size_t i = 0; i < amplitudes.size(); ++i )

          auto cont = max_sig_peaks.front()->continuum();
          cont->setRange( channel_energies[trial_start_chnl], channel_energies[trial_end_chnl+1] );

          if( step_continuum )
          {
            if( num_polynomial_terms == 2 )
              cont->setType( PeakContinuum::OffsetType::FlatStep );
            else if( num_polynomial_terms == 3 )
              cont->setType( PeakContinuum::OffsetType::LinearStep );
            else if( num_polynomial_terms == 4 )
              cont->setType( PeakContinuum::OffsetType::BiLinearStep );
            else
            {
              assert( 0 );
              throw logic_error( "Unknown step cont type" );
            }
          }else
          {
            if( num_polynomial_terms == 0 )
              cont->setType( PeakContinuum::OffsetType::NoOffset );
            else if( num_polynomial_terms == 1 )
              cont->setType( PeakContinuum::OffsetType::Constant );
            else if( num_polynomial_terms == 2 )
              cont->setType( PeakContinuum::OffsetType::Linear );
            else if( num_polynomial_terms == 3 )
              cont->setType( PeakContinuum::OffsetType::Quadratic );
            else if( num_polynomial_terms == 4 )
              cont->setType( PeakContinuum::OffsetType::Cubic );
            else
            {
              assert( 0 );
              throw logic_error( "Unknown cont type" );
            }
          }//if( step_continuum ) / else

          cont->setParameters( ref_energy, continuum_coeffs, continuum_coeffs_uncerts );

          // Anything else to set?
        }//if( nsigma > max_stat_sig )

        //if( (mean > 950 && mean < 975) || (mean > 1080 && mean < 1100) )
        if( false )
        {
          cout << "At " << added_mean << " keV, before adding peak, chi2=" << without_new_peak_chi2
          << ", and after adding chi2=" << with_new_peak_chi2
          << " - new-peak amp is " << amplitudes.back() << " +- " << amplitudes_uncerts.back()
          << " (" << ((nsigma > fit_settings.ROI_add_nsigma_required) ? "accept" : "reject") << ")" << endl
          << "\tThe chi2 improvements was " << chi2_improve
          << " and impr/dof=" << chi2_dof_improve << " (" << ((chi2_dof_improve > fit_settings.ROI_add_chi2dof_improve) ? "accept" : "reject") << ")" << endl
          << endl;
          cout << endl;
        }

        // Need to check that amplitude of new peak is positive, greater than a few sigma significant
        //  Also check that overall peak area went up, and anything else?
      }catch( std::exception &e )
      {
#ifndef NDEBUG
        cerr << "Caught exception trying to add a peak: " << e.what() << endl;
#endif
      }//try / catch
    }//for( size_t chnl = start_chnl; chnl < end_chnl; ++chnl )

    if( (max_stat_sig > 0.0) && !max_sig_peaks.empty() )
    {
      // We already checked that this candidate set of peaks is good to go, so we'll use them.
      added_any_peaks = true;

      if( PeakFitImprove::debug_printout )
      {
        shared_ptr<PeakDef> added_peak = max_sig_peaks.back();

        cout << "Adding peak at " << added_peak->mean() << " with area=" << added_peak->amplitude()
        << " +- " << added_peak->amplitudeUncert() << ", stat_sig=" << max_stat_sig
        << ", chi2_improve=" << max_stat_chi2_improve
        << ", chi2_dof_improve=" << max_stat_chi2_dof_improve
        << ". Other peaks in ROI: {";
        for( size_t i = 0; (i+1) < max_sig_peaks.size(); ++i )
          cout << (i ? ", " : "") << max_sig_peaks[i]->mean()
               << " (nsig=" << (max_sig_peaks[i]->amplitude() / max_sig_peaks[i]->amplitudeUncert())
               << ")";
        cout << "}" << endl;
      }//if( PeakFitImprove::debug_printout )

      for( const shared_ptr<const PeakDef> &peak : peaks_in_roi )
      {
        auto pos = std::find( begin(peaks_with_add), end(peaks_with_add), peak );
        assert( pos != end(peaks_with_add) );
        peaks_with_add.erase( pos );
      }

      roi_to_peaks[cont_peak_index].first = max_sig_peaks.front()->continuum();
      roi_to_peaks[cont_peak_index].second.clear();
      for( const shared_ptr<PeakDef> &p : max_sig_peaks )
      {
        peaks_with_add.push_back( p );
        roi_to_peaks[cont_peak_index].second.push_back( p );
      }

      // Lets go through and try to add another peak
      cont_peak_index -= 1;
    }//if( max_stat_sig > 0.0 && !max_sig_peaks.empty() )
  }//for( const auto &c_p : roi_to_peaks )


  if( added_any_peaks )
  {
    initial_fit_results.clear();
    for( const auto &p : peaks_with_add )
      initial_fit_results.push_back( *p );
    std::sort( begin(initial_fit_results), end(initial_fit_results), &PeakDef::lessThanByMean );
  }//if( added_any_peaks )

  return initial_fit_results;
}//vector<PeakDef> initial_peak_find_and_fit(...)


PeakFindAndFitWeights eval_initial_peak_find_and_fit( const InitialPeakFindSettings &fit_settings,
                                      const FindCandidateSettings &candidate_settings,
                                      const DataSrcInfo &src_info,
                                                     const bool multithread )
{
  const InjectSourceInfo &src = src_info.src_info;
  assert( src.src_spectra.size() >= 12 );
  if( src.src_spectra.size() < 12 )
    throw runtime_error( "Unexpected number of measurements." );

  const shared_ptr<const SpecUtils::Measurement> &data = src.src_spectra.front(); // TODO: we cold loop over all 16 of these histograms

  // Next two counters tell us how successful adding peaks were using the sliding linear continuum
  //  line, after the initial peak candidate find and fit.
  size_t num_add_candidates_fit_for = 0, num_add_candidates_accepted = 0;

  // `initial_peaks` not const for initial development
  vector<PeakDef> initial_peaks
                        = InitialFit_GA::initial_peak_find_and_fit( fit_settings, candidate_settings, data, multithread,
                                          num_add_candidates_fit_for, num_add_candidates_accepted);

  const vector<ExpectedPhotopeakInfo> &expected_photopeaks = src_info.expected_photopeaks;

  vector<PeakDef> found_not_expected;

  map<const ExpectedPhotopeakInfo *,vector<PeakDef>> found_maybe_wanted, found_def_wanted;


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


  // We'll just pick the closest expected peak in energy, to the detected peak.
  //   This isnt perfect, but close enough for our purposes, maybe
  for( PeakDef &found_peak : initial_peaks )
  {
    const double found_energy = found_peak.mean();
    const ExpectedPhotopeakInfo *nearest_expected_peak = nullptr;
    for( const ExpectedPhotopeakInfo &expected_peak : src_info.expected_photopeaks )
    {
      const double expected_energy = expected_peak.effective_energy;

      if( (found_energy >= expected_peak.roi_lower) && (found_energy <= expected_peak.roi_upper) )
      {
        if( !nearest_expected_peak
           || fabs(nearest_expected_peak->effective_energy - found_energy) > fabs(expected_energy - found_energy) )
        {
          nearest_expected_peak = &expected_peak;
        }
      }
    }//for( const ExpectedPhotopeakInfo &expected_peak : src_info.expected_photopeaks )


    if( nearest_expected_peak )
    {
      if( (nearest_expected_peak->nsigma_over_background >= JudgmentFactors::def_want_nsigma)
         && (nearest_expected_peak->peak_area > JudgmentFactors::min_def_wanted_counts) )
      {
        found_peak.setLineColor( Wt::GlobalColor::blue );
        found_peak.setUserLabel( "Found - def wanted: A=" + to_string( nearest_expected_peak->peak_area ) );
        found_def_wanted[nearest_expected_peak].push_back( found_peak );
      }else
      {
        found_peak.setLineColor( Wt::GlobalColor::cyan );
        found_peak.setUserLabel( "Found - opt Wanted: A=" + to_string( nearest_expected_peak->peak_area ) );
        found_maybe_wanted[nearest_expected_peak].push_back( found_peak );
      }
    }else
    {
      found_peak.setUserLabel( "Not Expected" );
      found_peak.setLineColor( Wt::GlobalColor::red );


      static const Wt::WColor purple( 235, 33, 188 );


      // If HPGe, check if this is a single or double escape peak - we dont have these in the
      //  "truth" lines.
      if( (found_energy > 508) && (found_energy < 514) )
      {
        found_peak.setUserLabel( "Annih." );
        found_peak.setLineColor( purple );
      }else if( data->num_gamma_channels() > 4098 )
      {
        const double range = found_peak.fwhm();
        const double se_parent = found_energy + 510.9989;
        const double de_parent = found_energy + 2*510.9989;

        for( const auto &ep : possible_escape_peak_parents )
        {
          if( fabs(get<0>(ep) - se_parent) < range )
          {
            found_peak.setLineColor( purple );
            found_peak.setUserLabel( "Possibly S.E. of " + to_string( get<0>(ep) ) );
            break;
          }else if( fabs(get<0>(ep) - de_parent) < range )
          {
            found_peak.setLineColor( purple );
            found_peak.setUserLabel( "Possibly D.E. of " + to_string( get<0>(ep) ) );
            break;
          }
        }//for( loop over possible escape peak parents )
      }//if( data->num_gamma_channels() > 4098 )

      found_not_expected.push_back( found_peak );
    }//if( nearest_expected_peak )
  }//for( const PeakDef &found_peak : initial_peaks )

  vector<const ExpectedPhotopeakInfo *> def_wanted_not_found, maybe_wanted_not_found;
  for( const ExpectedPhotopeakInfo &epi : expected_photopeaks )
  {
    if( (epi.nsigma_over_background >= JudgmentFactors::def_want_nsigma)
       && (epi.peak_area >= JudgmentFactors::min_def_wanted_counts)
       && !found_def_wanted.count(&epi) )
    {
      def_wanted_not_found.push_back( &epi );
    }//

    if( (epi.nsigma_over_background < JudgmentFactors::def_want_nsigma)
       && !found_maybe_wanted.count(&epi) )
    {
      maybe_wanted_not_found.push_back( &epi );
    }//
  }//for( const ExpectedPhotopeakInfo &epi : expected_photopeaks )

  if( PeakFitImprove::debug_printout )
  {
    cout << "Found " << found_def_wanted.size() << " of " << (def_wanted_not_found.size() + found_def_wanted.size())
    << " peaks definetly wanted.  Found " << found_maybe_wanted.size() << " of " << (maybe_wanted_not_found.size() + found_maybe_wanted.size())
    << " maybe wanted peak.  Found " << found_not_expected.size() << " unexpected peaks." << endl;
  }

  double score = found_def_wanted.size();

  for( const auto &pp : found_maybe_wanted )
  {
    const ExpectedPhotopeakInfo * epi = pp.first;
    if( epi->nsigma_over_background > JudgmentFactors::lower_want_nsigma )
    {
      if( epi->peak_area > JudgmentFactors::min_def_wanted_counts )
      {
        const double amount_short = JudgmentFactors::def_want_nsigma - epi->nsigma_over_background;
        const double fraction_short = amount_short / (JudgmentFactors::def_want_nsigma - JudgmentFactors::lower_want_nsigma);
        assert( (fraction_short >= 0.0) && (fraction_short <= 1.0) );
        score += JudgmentFactors::min_initial_fit_maybe_want_score + (1.0 - JudgmentFactors::min_initial_fit_maybe_want_score)*(1 - fraction_short);
      }else
      {
        const double fraction_short = 1.0 - (epi->peak_area / JudgmentFactors::min_def_wanted_counts);
        assert( (fraction_short >= 0.0) && (fraction_short <= 1.0) );
        score += JudgmentFactors::min_initial_fit_maybe_want_score + (1.0 - JudgmentFactors::min_initial_fit_maybe_want_score)*(1 - fraction_short);
      }
    }else
    {
      score += 0.0; //No punishment, but no reward.
    }
  }//for( const ExpectedPhotopeakInfo * epi : found_maybe_wanted )

  size_t num_found_not_expected = found_not_expected.size();

  // The 511 Annih. line isnt in our truth line, unless the source makes them, so lets
  //  remove this one from the score
  for( const auto &p : found_not_expected )
  {
    // we dont have "truth" lines for annihilation or escape peaks
    if( SpecUtils::icontains(p.userLabel(), "S.E.")
             || SpecUtils::icontains(p.userLabel(), "D.E.")
             || SpecUtils::icontains(p.userLabel(), "Annih.") )
    {
      if( num_found_not_expected > 0 )
        num_found_not_expected -= 1;
    }
  }

  score -= JudgmentFactors::initial_fit_extra_peak_punishment * num_found_not_expected;


  // Now punish for extra time fitting more peaks
  assert( num_add_candidates_fit_for >= num_add_candidates_accepted );
  if( num_add_candidates_fit_for > 0 )
    score -= JudgmentFactors::extra_add_fits_punishment * (num_add_candidates_fit_for - num_add_candidates_accepted) / num_add_candidates_fit_for;

  if( PeakFitImprove::debug_printout )
  {
    SpecMeas output;
    auto new_meas = make_shared<SpecUtils::Measurement>( *data );
    new_meas->set_sample_number( 1 );
    output.add_measurement( new_meas, false );

    deque<shared_ptr<const PeakDef>> peaks;
    for( const PeakDef &p : initial_peaks )
    {
      auto peak = make_shared<PeakDef>( p );
      peaks.push_back( peak );
    }

    //Now add missing definitely (and optionally) wanted peaks, and color them orange (yellow)
    vector<const ExpectedPhotopeakInfo *> missing_peaks = def_wanted_not_found;
    missing_peaks.insert( end(missing_peaks), begin(maybe_wanted_not_found), end(maybe_wanted_not_found) );

    for( size_t missing_index = 0; missing_index < missing_peaks.size(); ++missing_index )
    {
      const ExpectedPhotopeakInfo &missing = *missing_peaks[missing_index];
      const bool def_wanted = (missing_index < def_wanted_not_found.size());

      const double mean = missing.effective_energy;
      const double sigma = missing.gamma_lines.front().fwhm/2.35482f;
      auto peak = make_shared<PeakDef>( mean, sigma, missing.peak_area );
      peak->setFitFor( PeakDef::CoefficientType::Mean, false );
      peak->setFitFor( PeakDef::CoefficientType::Sigma, false );
      peak->setFitFor( PeakDef::CoefficientType::GaussAmplitude, false );
      peak->continuum()->setType( PeakContinuum::OffsetType::Linear );
      peak->continuum()->setRange( mean - 3*sigma, mean + 3*sigma );
      peak->continuum()->calc_linear_continuum_eqn( data, mean, mean - 4*sigma, mean + 4*sigma, 5, 5 );
      const Wt::WColor orange(252, 94, 3), yellow(168, 150, 50);
      peak->setLineColor( def_wanted ? orange : yellow );
      string label = def_wanted ? "Not det. - wanted" : "Not det. - poss.";
      label += ", N_Sigma=" + SpecUtils::printCompact(missing.nsigma_over_background, 3);
      peak->setUserLabel( label );

      peaks.push_back( peak );
    }//for( const ExpectedPhotopeakInfo &missing : not_detected )


    output.setPeaks( peaks, {1} );
    output.cleanup_after_load( SpecUtils::SpecFile::DontChangeOrReorderSamples );

    string outdir = "output_n42";
    if( !SpecUtils::is_directory(outdir) && SpecUtils::create_directory(outdir) )
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

    const string out_n42 = SpecUtils::append_path( outdir, src_info.src_info.src_name ) + "_initial_fit.n42";

    output.save2012N42File( out_n42, [=](){
      cerr << "Failed to write '" << out_n42 << "'" << endl;
    });

    cout << "Wrote: " << out_n42 << endl;
  }//if( write output N42 file to inspect


  const auto sum_area_weight = []( const map<const ExpectedPhotopeakInfo *,vector<PeakDef>> &data ) -> pair<double,double> {
    if( data.empty() )
      return {0.0, 0.0};
    vector<double> weights;
    weights.reserve( data.size() );

    double area_weight = 0;
    for( const auto &pp : data )
    {
      const ExpectedPhotopeakInfo * epi = pp.first;
      double fit_area = 0;
      for( const PeakDef &p : pp.second )
        fit_area += p.amplitude();
      const double sigma = (epi->peak_area < 10000.0) ? sqrt(epi->peak_area) : 0.01*epi->peak_area;
      const double weight = fabs(epi->peak_area - fit_area) / sigma;
      area_weight += weight;
      weights.push_back(weight);
    }

    std::sort( begin(weights), end(weights) );

    return {area_weight / data.size(), weights[weights.size()/2] };
  };//sum_area_weight lambda


  vector<double> not_wanted_area_weights;
  double not_wanted_area_weight = 0.0;
  for( const PeakDef &p : found_not_expected )
  {
    if( ((p.peakAreaUncert() > 0.001) && (p.amplitude() < 10.0*p.peakAreaUncert())) || (p.amplitude() > 25) )
    {
      not_wanted_area_weights.push_back( sqrt(p.peakArea()) );
      not_wanted_area_weight += sqrt( p.peakArea() );
    }
  }


  PeakFindAndFitWeights answer;
  answer.find_weight = score;
  const pair<double,double> def_want_weights = sum_area_weight( found_def_wanted );
  answer.def_wanted_area_weight = def_want_weights.first;
  answer.def_wanted_area_median_weight = def_want_weights.second;

  const pair<double,double> maybe_want_weights = sum_area_weight( found_maybe_wanted );
  answer.maybe_wanted_area_weight = maybe_want_weights.first;
  answer.maybe_wanted_area_median_weight = maybe_want_weights.second;

  if( not_wanted_area_weights.size() )
  {
    answer.not_wanted_area_weight = not_wanted_area_weight;
    answer.not_wanted_area_median_weight = not_wanted_area_weights[not_wanted_area_weights.size()/2];
  }else
  {
    answer.not_wanted_area_weight = 0.0;
    answer.not_wanted_area_median_weight = 0.0;
  }

  return answer;
}//double eval_initial_peak_find_and_fit(...)



string InitialFitSolution::to_string( const string &separator ) const
{
  return
  string("initial_stat_threshold: ")           + std::to_string(initial_stat_threshold)
  + separator + "initial_hypothesis_threshold: " + std::to_string(initial_hypothesis_threshold)
  + separator + "initial_min_nsigma_roi: "       + std::to_string(initial_min_nsigma_roi)
  + separator + "initial_max_nsigma_roi: "       + std::to_string(initial_max_nsigma_roi)
  + separator + "fwhm_fcn_form: "                + std::to_string(fwhm_fcn_form)
  + separator + "search_roi_nsigma_deficit: "    + std::to_string(search_roi_nsigma_deficit)
  + separator + "search_stat_threshold: "        + std::to_string(search_stat_threshold)
  + separator + "search_hypothesis_threshold: "  + std::to_string(search_hypothesis_threshold)
  + separator + "search_stat_significance: "     + std::to_string(search_stat_significance)
  + separator + "ROI_add_nsigma_required: "      + std::to_string(ROI_add_nsigma_required)
  + separator + "ROI_add_chi2dof_improve: "      + std::to_string(ROI_add_chi2dof_improve)
  ;
}//to_string( separator )



InitialPeakFindSettings genes_to_settings( const InitialFitSolution &solution )
{
  InitialPeakFindSettings settings;

  settings.initial_stat_threshold = solution.initial_stat_threshold;
  settings.initial_hypothesis_threshold = solution.initial_hypothesis_threshold;

  settings.initial_min_nsigma_roi = solution.initial_min_nsigma_roi;
  settings.initial_max_nsigma_roi = solution.initial_max_nsigma_roi;

  settings.fwhm_fcn_form = static_cast<InitialPeakFindSettings::FwhmFcnForm>( solution.fwhm_fcn_form );
  assert( settings.fwhm_fcn_form == InitialPeakFindSettings::FwhmFcnForm::Gadras
         || settings.fwhm_fcn_form == InitialPeakFindSettings::FwhmFcnForm::SqrtPolynomialTwoCoefs
         || settings.fwhm_fcn_form == InitialPeakFindSettings::FwhmFcnForm::SqrtPolynomialThreeCoefs
         || settings.fwhm_fcn_form == InitialPeakFindSettings::FwhmFcnForm::SqrtEnergyPlusInverse
         //|| settings.fwhm_fcn_form == InitialPeakFindSettings::FwhmFcnForm::NumFwhmFcnForm
         );

  settings.search_roi_nsigma_deficit = solution.search_roi_nsigma_deficit;
  settings.search_stat_threshold = solution.search_stat_threshold;
  settings.search_hypothesis_threshold = solution.search_hypothesis_threshold;
  settings.search_stat_significance = solution.search_stat_significance;
  settings.ROI_add_nsigma_required = solution.ROI_add_nsigma_required;
  settings.ROI_add_chi2dof_improve = solution.ROI_add_chi2dof_improve;

  return settings;
}//InitialPeakFindSettings genes_to_settings( const InitialFitSolution &solution )


void init_genes(InitialFitSolution& p,const std::function<double(void)> &rnd01)
{
  // rnd01() gives a random number in 0~1
  p.initial_stat_threshold       = 0.0 + 8*rnd01();
  p.initial_hypothesis_threshold = 0.25 + 2.75*rnd01(); //best values tend to be about 1.0
  p.initial_min_nsigma_roi       = 2.0 + 6*rnd01();
  p.initial_max_nsigma_roi       = 4.0 + 8*rnd01();
  p.fwhm_fcn_form                = static_cast<int>( InitialPeakFindSettings::FwhmFcnForm::SqrtPolynomialTwoCoefs ); //0+3*rnd01();
  p.search_roi_nsigma_deficit    = 3.0 + 7*rnd01();
  p.search_stat_threshold        = 1.0 + 8*rnd01();
  p.search_hypothesis_threshold  = -0.1+10.1*rnd01();
  p.search_stat_significance     = 1.0 + 5*rnd01();
  p.ROI_add_nsigma_required      = 1.0 + 7*rnd01();
  p.ROI_add_chi2dof_improve      = 0.0 + 8*rnd01();
}


bool eval_solution( const InitialFitSolution &p, InitialFitCost &c )
{
  const InitialPeakFindSettings settings = genes_to_settings( p );

  c.objective1 = ns_ga_eval_fcn( settings );

  ns_num_evals_this_generation += 1;
  if( (ns_num_evals_this_generation % 10) == 0 )
    cout << "Have evaluated " << ns_num_evals_this_generation.load() << " individuals this generation." << endl;


  if( IsNan(c.objective1) || IsInf(c.objective1) )
  {
    cerr << "Got an objective of " << c.objective1 << " for " << p.to_string(", ") << endl;
    return false;
  }

  return true; // solution is accepted
}

InitialFitSolution mutate( const InitialFitSolution& X_base,
                          const std::function<double(void)> &rnd01,
                          double shrink_scale )
{
  InitialFitSolution X_new;
  const double mu = 0.2*shrink_scale; // mutation radius (adjustable)
  const double mutate_threshold = PeakFitImprove::sm_ga_mutate_threshold;  //one paper suggest 20%
  bool in_range;

  size_t num_tries = 0;

  do{
    num_tries += 1;
    if( num_tries > 1000 )
    {
      cerr << "Has taken over " << num_tries << " tries to find some genes.\nX_new={\n"
      << X_new.to_string( "\n\t" )
      << "}\nX_base={\n"
      << X_base.to_string("\n\t" )
      << "}\nWill return new randomly inited gene.\n"
      << endl;

      std::random_device rng;
      std::uniform_real_distribution<double> unif_dist( 0.0, 1.0 );

      init_genes( X_new, [&](){return unif_dist(rng);} );
      return X_new;
    }

    in_range = true;
    X_new = X_base;
    if( rnd01() < mutate_threshold )
    {
      X_new.initial_stat_threshold += mu*(rnd01()-rnd01());
      in_range = in_range&&(X_new.initial_stat_threshold>=0.0 && X_new.initial_stat_threshold<8);
    }

    if( rnd01() < mutate_threshold )
    {
      X_new.initial_hypothesis_threshold+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.initial_hypothesis_threshold>=0.25 && X_new.initial_hypothesis_threshold<3.0);
    }

    if( rnd01() < mutate_threshold )
    {
      X_new.initial_min_nsigma_roi+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.initial_min_nsigma_roi>=2 && X_new.initial_min_nsigma_roi<8);
    }

    if( rnd01() < mutate_threshold )
    {
      X_new.initial_max_nsigma_roi+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.initial_max_nsigma_roi>=4 && X_new.initial_max_nsigma_roi<12);
    }

    X_new.fwhm_fcn_form = X_base.fwhm_fcn_form; //unnecessary, but to be explicit
    //X_new.fwhm_fcn_form+=mu*(rnd01()-rnd01());
    //in_range=in_range&&(X_new.fwhm_fcn_form>=0 && X_new.fwhm_fcn_form<3);

    if( rnd01() < mutate_threshold )
    {
      X_new.search_roi_nsigma_deficit+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.search_roi_nsigma_deficit>=3 && X_new.search_roi_nsigma_deficit < 10.0);
    }

    if( rnd01() < mutate_threshold )
    {
      X_new.search_stat_threshold+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.search_stat_threshold>=1 && X_new.search_stat_threshold<9);
    }

    if( rnd01() < mutate_threshold )
    {
      X_new.search_hypothesis_threshold+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.search_hypothesis_threshold>=-0.1 && X_new.search_hypothesis_threshold<10.0);
    }

    if( rnd01() < mutate_threshold )
    {
      X_new.search_stat_significance+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.search_stat_significance>=1 && X_new.search_stat_significance<6);
    }

    if( rnd01() < mutate_threshold )
    {
      X_new.ROI_add_nsigma_required+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.ROI_add_nsigma_required>=1 && X_new.ROI_add_nsigma_required<8);
    }

    if( rnd01() < mutate_threshold )
    {
      X_new.ROI_add_chi2dof_improve+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.ROI_add_chi2dof_improve>=0.0 && X_new.ROI_add_chi2dof_improve<8);
    }
  } while(!in_range);
  return X_new;
}

InitialFitSolution crossover( const InitialFitSolution& X1,
                             const InitialFitSolution& X2,
                             const std::function<double(void)> &rnd01 )
{
  const double crossover_threshold = PeakFitImprove::sm_ga_crossover_threshold;  //One paper suggest 80% crossover rate

  InitialFitSolution X_new;
  X_new = X1;

  double r;
  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.initial_stat_threshold=r*X1.initial_stat_threshold+(1.0-r)*X2.initial_stat_threshold;
  }else if( rnd01() < 0.5 )
  {
    X_new.initial_stat_threshold = X2.initial_stat_threshold;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.initial_hypothesis_threshold=r*X1.initial_hypothesis_threshold+(1.0-r)*X2.initial_hypothesis_threshold;
  }else if( rnd01() < 0.5 )
  {
    X_new.initial_hypothesis_threshold = X2.initial_hypothesis_threshold;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.initial_min_nsigma_roi=r*X1.initial_min_nsigma_roi+(1.0-r)*X2.initial_min_nsigma_roi;
  }else if( rnd01() < 0.5 )
  {
    X_new.initial_min_nsigma_roi = X2.initial_min_nsigma_roi;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.initial_max_nsigma_roi=r*X1.initial_max_nsigma_roi+(1.0-r)*X2.initial_max_nsigma_roi;
  }else if( rnd01() < 0.5 )
  {
    X_new.initial_max_nsigma_roi = X2.initial_max_nsigma_roi;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.fwhm_fcn_form=r*X1.fwhm_fcn_form+(1.0-r)*X2.fwhm_fcn_form;
  }else if( rnd01() < 0.5 )
  {
    X_new.fwhm_fcn_form = X2.fwhm_fcn_form;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.search_roi_nsigma_deficit=r*X1.search_roi_nsigma_deficit+(1.0-r)*X2.search_roi_nsigma_deficit;
  }else if( rnd01() < 0.5 )
  {
    X_new.search_roi_nsigma_deficit = X2.search_roi_nsigma_deficit;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.search_stat_threshold=r*X1.search_stat_threshold+(1.0-r)*X2.search_stat_threshold;
  }else if( rnd01() < 0.5 )
  {
    X_new.search_stat_threshold = X2.search_stat_threshold;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.search_hypothesis_threshold=r*X1.search_hypothesis_threshold+(1.0-r)*X2.search_hypothesis_threshold;
  }else if( rnd01() < 0.5 )
  {
    X_new.search_hypothesis_threshold = X2.search_hypothesis_threshold;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.search_stat_significance=r*X1.search_stat_significance+(1.0-r)*X2.search_stat_significance;
  }else if( rnd01() < 0.5 )
  {
    X_new.search_stat_significance = X2.search_stat_significance;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.ROI_add_nsigma_required=r*X1.ROI_add_nsigma_required+(1.0-r)*X2.ROI_add_nsigma_required;
  }else if( rnd01() < 0.5 )
  {
    X_new.ROI_add_nsigma_required = X2.ROI_add_nsigma_required;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.ROI_add_chi2dof_improve=r*X1.ROI_add_chi2dof_improve+(1.0-r)*X2.ROI_add_chi2dof_improve;
  }else if( rnd01() < 0.5 )
  {
    X_new.ROI_add_chi2dof_improve = X2.ROI_add_chi2dof_improve;
  }

  
  return X_new;
}

double calculate_SO_total_fitness(const GA_Type::thisChromosomeType &X)
{
  // finalize the cost
  double final_cost=0.0;
  final_cost+=X.middle_costs.objective1;
  return final_cost;
}




void SO_report_generation( int generation_number,
                          const EA::GenerationType<InitialFitSolution,InitialFitCost> &last_generation,
                          const InitialFitSolution& best_genes)
{
  bool best_yet = false;
  if( !sm_set_best_genes || (last_generation.best_total_cost < sm_best_total_cost) )
  {
    best_yet = true;
    sm_set_best_genes = true;
    sm_best_genes = best_genes;
    sm_best_total_cost = last_generation.best_total_cost;
  }

  sm_output_file
  << generation_number <<"\t"
  << last_generation.average_cost <<"\t"
  << last_generation.best_total_cost <<"\t"
  << "{" << best_genes.to_string(", ") << "}"
  << "\t" << (best_yet ? "BestYet" : "NoImprovement")
  << "\n\n";

  cout
  <<"Generation ["<<generation_number<<"], "
  <<"Best="<<last_generation.best_total_cost<<", "
  <<"Average="<<last_generation.average_cost<<", "
  << ", " << (best_yet ? "Best generation yet" : "no improvement")
  <<"\n  Best genes: {\n\t" <<best_genes.to_string("\n\t")  << "\n}\n"
  <<"Exe_time="<<last_generation.exe_time
  << endl;

  ns_num_evals_this_generation = 0;
}//SO_report_generation( ... )


InitialPeakFindSettings do_ga_eval( std::function<double( const InitialPeakFindSettings &)> ga_eval_fcn )
{
  assert( !sm_has_been_called );
  if( sm_has_been_called )
  {
    cerr << "You may only call InitialFit_GA::do_ga_eval(...) once per program execution." << endl;
    exit(1);
  }

  sm_has_been_called = true;

  assert( RETURN_PeakDef_Candidates ); //dont want to do a static_assert - since I am still compiling all this code both ways...
  if( !RETURN_PeakDef_Candidates )
  {
    cerr << "Please change 'RETURN_PeakDef_Candidates' to true and recompile." << endl;
    exit(1);
  }

  assert( !!ga_eval_fcn );
  if( !ga_eval_fcn )
    throw runtime_error( "Invalid eval function passed in." );

  ns_ga_eval_fcn = ga_eval_fcn;

  sm_output_file.open( "initial_fit_results.txt" );
  sm_output_file << "step" << "\t" << "cost_avg" << "\t" << "cost_best" << "\t" << "solution_best" << "\n";

  EA::Chronometer timer;
  timer.tic();

  GA_Type ga_obj;
  ga_obj.problem_mode=EA::GA_MODE::SOGA;
  ga_obj.multi_threading = true;
  ga_obj.idle_delay_us = 1; // switch between threads quickly
  ga_obj.dynamic_threading = (PeakFitImprove::sm_ga_population > PeakFitImprove::sm_num_optimization_threads);
  ga_obj.verbose = false;
  ga_obj.population = static_cast<unsigned int>(PeakFitImprove::sm_ga_population);
  ga_obj.N_threads = static_cast<int>(PeakFitImprove::sm_num_optimization_threads);
  ga_obj.generation_max = static_cast<int>(PeakFitImprove::sm_ga_generation_max);
  ga_obj.calculate_SO_total_fitness=calculate_SO_total_fitness;
  ga_obj.init_genes=init_genes;
  ga_obj.eval_solution=eval_solution;
  ga_obj.mutate=mutate;
  ga_obj.crossover=crossover;
  ga_obj.SO_report_generation=SO_report_generation;
  ga_obj.crossover_fraction = PeakFitImprove::sm_ga_crossover_fraction; //This *looks* to be the fraction of individuals that will have crossover applied to them
  ga_obj.mutation_rate=PeakFitImprove::sm_ga_mutation_rate; //This *looks* to be the fraction of individuals that will have mutattion applied to them
  ga_obj.best_stall_max = static_cast<int>(PeakFitImprove::sm_ga_best_stall_max);
  ga_obj.elite_count = static_cast<int>(PeakFitImprove::sm_ga_elite_count);
  EA::StopReason stop_reason = ga_obj.solve();

  cout<<"The problem is optimized in "<<timer.toc()<<" seconds."<<endl;

  sm_output_file.close();

  cout << "Stop reason was: " << ga_obj.stop_reason_to_string( stop_reason) << endl;
  cout << "The problem is optimized in " << timer.toc() << " seconds." << endl;

  sm_output_file.close();

  return genes_to_settings( sm_best_genes );
}
}//namespace InitialFit_GA
