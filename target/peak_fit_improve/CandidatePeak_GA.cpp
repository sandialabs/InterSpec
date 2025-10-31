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
#include <limits>
#include <string>
#include <fstream>
#include <iostream>

#include <boost/math/constants/constants.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"

#include "PeakFitImprove.h"
#include "CandidatePeak_GA.h"

using namespace std;


std::mutex EA::mtx_rand;


namespace
{
std::function<double(const FindCandidateSettings &)> ns_ga_eval_fcn;

std::ofstream sm_output_file;

bool sm_set_best_genes = false;
CandidatePeak_GA::CandidatePeakSolution sm_best_genes;
double sm_best_total_cost = 1.0E99;

bool sm_has_been_called = false;
}


namespace CandidatePeak_GA
{

std::vector<PeakDef> find_candidate_peaks( const std::shared_ptr<const SpecUtils::Measurement> data,
                          size_t start_channel,
                          size_t end_channel,
                          const FindCandidateSettings &settings )
{
  // Plan:
  //  - Should we split the spectrum up into N regions?  So we can refine have fewer smoothing bins lower in energy?
  //    - For HPGe it looks like we can avoid this (but out of time to investigate)
  std::vector<PeakDef> result_peaks;

  const size_t nchannel = data->num_gamma_channels();

  //We will let one bin fluctate negative to avoid fluctations near threshold_FOM
  //Untested for HPGe data.
  //  Currently (20141209) for low res data, this should be kept the same as
  //  in find_roi_for_2nd_deriv_candidate(...) {although all this is under
  //  development}.
  const size_t nFluxuate = settings.num_chan_fluctuate;

  if( nchannel < (2*nFluxuate + 2) || (nchannel < 4) )
    return {};

  if( start_channel >= nchannel )
    start_channel = 0;
  if( end_channel <= start_channel || end_channel >= (nchannel-1) )
    end_channel = nchannel - nFluxuate - 1;


  const double threshold_FOM = settings.threshold_FOM;
  const float pos_sum_threshold_sf = settings.pos_sum_threshold_sf;


  const int order = settings.smooth_polynomial_order; //highres ? 3 : 2;
  const size_t side_bins = settings.num_smooth_side_channels; //highres ? 4 : std::max( size_t(5), static_cast<size_t>( 0.022f*midenergy/midbinwidth + 0.5f ) );

  const vector<float> &spectrum = *data->gamma_counts();

  vector<float> second_deriv, smoothed_data, coarser_data, rougher_first_deriv, rougher_second_deriv, rougher_second_deriv_var;

  // Create a smoothed second order derivative; its negative-most values will correspond
  //  to peak centers
  smoothSpectrum( spectrum, static_cast<int>(side_bins), order, 2, second_deriv );

  // We will smooth the data a bit, so we can sum it to help determine stat significance
  smoothSpectrum( spectrum, static_cast<int>(side_bins), 1, 0, smoothed_data );

  smoothSpectrum( spectrum, 3, 3, 1, rougher_first_deriv );

  // We will use `rougher_second_deriv` to evaluate sigma, and ROI extents, because if
  // `side_bins`, than the predicted sigma and ROI extents will be a lot wider than
  // expected.
  {
    SavitzyGolayCoeffs sgcoeffs( 3, 3, 3, 2 );
    sgcoeffs.smooth_with_variance( spectrum, rougher_second_deriv, rougher_second_deriv_var );
  }


  //XXX: the below 1.5 is empiracally found, I'm not entirely sure where
  //     comes from...and infact might be higher
  const double amp_fake_factor = 1.0;


  const vector<float> &energies = *data->gamma_channel_energies();

  size_t minbin = 0, firstzero = 0, secondzero = 0;
  float secondsum = 0.0f, minval = 9999999999.9f;

  // Usually there is some turn-on for the spectrum, so lets skip past that
  //  by requiring 2nd derivative to go positive, then negative, then positive again
  if( start_channel == 0 )
  {
    while( (start_channel+nFluxuate) < end_channel )
    {
      bool all_neg = true;
      for( size_t ind = start_channel; ind < (start_channel+nFluxuate); ++ind )
        all_neg = (second_deriv[ind] < -1.0E-9);

      if( all_neg )
        break;
      ++start_channel;
    }

    while( (start_channel+nFluxuate) < end_channel )
    {
      bool all_pos = true;
      for( size_t ind = start_channel; ind < (start_channel+nFluxuate); ++ind )
        all_pos = (second_deriv[ind] > 1.0E-9);

      if( all_pos )
        break;
      ++start_channel;
    }
  }//if( start_channel == 0 )


  for( size_t channel = start_channel; channel <= end_channel; ++channel )
  {
    const double channel_energy = data->gamma_channel_center(channel);
    const bool debug_channel = (channel_energy > 1963 && channel_energy < 1996);

    const float secondDeriv = second_deriv[channel]; //Not dividing by binwidth^2 here,

    bool secondSumPositive = true;
    float positive2ndDerivSum = 0.0f, positive2ndDataSum = 0.0f;
    for( size_t i = 0; (i < nFluxuate) && ((channel+i) <= end_channel); ++i )
    {
      const bool above = (second_deriv[channel+i] > 0.0f);
      if( above )
      {
        positive2ndDerivSum += second_deriv[channel+i];
        positive2ndDataSum += smoothed_data[channel+i];
      }
      secondSumPositive &= above;
    }

    //Rather than using 'pos_sum_threshold_sf*secondsum', it might be better to
    //  use something involving sqrt(secondsum) since this makes a bit more
    //  sense.
    //Also, positive2ndDerivSum can probably also thresholded off of some sort of more
    //  absolute quantity
    secondSumPositive &= (positive2ndDerivSum > pos_sum_threshold_sf*secondsum);

    if( secondSumPositive && (minval < 99999999999.9)
       && (secondsum!=0.0) && (firstzero>0)
       && ((channel-firstzero)>2) )
    {
      /*
      float pos2ndDerivAreaSum = 0.0f;
      for( size_t ind = channel; ind <= end_channel; ++ind )
      {
        pos2ndDerivAreaSum += second_deriv[channel];

        bool all_neg = true;
        for( size_t fluc_ind = ind + 1; fluc_ind < (ind + nFluxuate + 1); ++fluc_ind )
          all_neg = (second_deriv[fluc_ind] <= 0.0);

        if( all_neg )
          break;
      }
       */

      if( pos_sum_threshold_sf < 0.0 )
      {
        // If `pos_sum_threshold_sf` is negative (which it usually is), then `channel` will be past
        //  the zero point, so lets try to go back to the most recent zero point.
        //  We want `secondzero` to be the first positive second-derivative channel,
        //  after the negative region
        size_t trial_second_zero = channel;
        while( (trial_second_zero > minbin)
              && (trial_second_zero > 1)
              && (second_deriv[trial_second_zero-1]  > 0.0) )
        {
          trial_second_zero -= 1;
        }

        secondzero = ((trial_second_zero > minbin) && second_deriv[trial_second_zero-1] < 0.0)
                        ? trial_second_zero
                        : channel;
      }else
      {
        secondzero = channel;
      }

      double mean = data->gamma_channel_center(minbin);

      // TODO: More carefully estimate lower and upper sigma-zeros, to have sub-channel accuracy
      //       Maybe fit a tine over the surrounding two (or at least <nFluxuate) bins on either side?
      //       Also, maybe use a less-smoothed 3rd order poly to estimate it
      const double sigma_smoothed = 0.5*(data->gamma_channel_center(secondzero)
                                         - data->gamma_channel_center(firstzero));
      double sigma = sigma_smoothed;


      /*
       // skewness tends to be reasonably near zero for legit peaks, but havent looked into
       //  any more detail about using `mean_channel`, `variance`, or `skewness` to help
       //  eliminate poor candidates
      //first, second, and third moments
      double mean_channel = 0, variance = 0, skewness = 0, i_sum = 0.0, second_deriv_sum = 0.0;
      for( size_t i = firstzero; i <= secondzero; ++i )
      {
        i_sum += i;
        second_deriv_sum += second_deriv[i];
        mean_channel += (i+0.5)*second_deriv[i];
      }
      mean_channel /= second_deriv_sum;

      for( size_t i = firstzero; i <= secondzero; ++i )
      {
        variance += second_deriv[i] * ((i+0.5) - mean_channel) * ((i+0.5) - mean_channel);
        skewness += second_deriv[i] * std::pow((i+0.5) - mean_channel, 3);
      }

      variance /= second_deriv_sum;
      skewness = skewness / (second_deriv_sum * std::pow(variance, 1.5));

      cerr << "mean=" << mean << ", minbin=" << minbin
       << ", firstzero=" << firstzero << ", secondzero=" << secondzero
       << ", mean_channel=" << mean_channel << ", sqrt(variance)=" << sqrt(variance)
       << " (s=" << sigma_smoothed << ")"
       << ", skewness=" << skewness
       << endl;
      */


      /*
      if( second_deriv[firstzero+1] > 0.0 )
      {
        for( int i = -5; i <= 5; ++i )
          cerr << "second_deriv[firstzero" << (i>=0 ? "+" : "") << i << "] = " << second_deriv[firstzero+i] << endl;
      }

      if( second_deriv[secondzero-1] > 0.0 )
      {
        for( int i = -5; i <= 5; ++i )
          cerr << "second_deriv[secondzero" << (i>=0 ? "+" : "") << i << "] = " << second_deriv[secondzero+i] << endl;
      }
       */

      assert( second_deriv[firstzero+1] <= 0.0 );
      //assert( second_deriv[secondzero-1] <= 0.0 );  For peaks with min-bin right on the edge, this test would fail


      // The more-coarse smoothed value at `minbin` may not be negative, due to a large-ish
      //  fluctuation, so we'll also look on either side of that channel (although I expect this
      //  is really rare).
      size_t rough_index = (rougher_second_deriv[minbin] < rougher_second_deriv[minbin-1]) ? minbin : minbin-1;
      rough_index = (rougher_second_deriv[rough_index] < rougher_second_deriv[minbin+1]) ? rough_index : minbin+1;

      double rougher_FOM = 0.0;

      // We'll take the ROI to include the first negative bins, after 2nd deriv gos positive
      size_t roi_begin_index = 0, roi_end_index = 0;

      bool rougher_confident_multi_roi = false;
      const double rougher_scnd_drv = rougher_second_deriv[rough_index];

      if( debug_channel )
        cout << "Candidate: mean=" << mean << ", rougher_scnd_drv=" << rougher_scnd_drv
        << ", minbin=" << minbin << ", firstzero=" << firstzero
        << ", secondzero=" << secondzero << endl;


      if( rougher_scnd_drv < 0.0 )
      {
        //size_t nFluxuate = 1; //
        size_t lower_pos_index = rough_index, upper_pos_index = rough_index;
        while( lower_pos_index > nFluxuate )
        {
          bool all_pos = true;
          for( size_t ind = lower_pos_index; all_pos && (ind > (lower_pos_index - nFluxuate)); --ind )
            all_pos = (rougher_second_deriv[ind] >= -std::numeric_limits<float>::epsilon());
          if( all_pos )
            break;

          lower_pos_index -= 1;
        }


        while( upper_pos_index < (end_channel - nFluxuate) )
        {
          bool all_pos = true;
          for( size_t ind = upper_pos_index; all_pos && (ind < (upper_pos_index + nFluxuate)); ++ind )
            all_pos = (rougher_second_deriv[ind] >= -std::numeric_limits<float>::epsilon());
          if( all_pos )
            break;

          upper_pos_index += 1;
        }

        assert( (rougher_second_deriv[lower_pos_index] >= 0)   || (lower_pos_index <= (nFluxuate+1)) );
        assert( (rougher_second_deriv[lower_pos_index+1] <= 0) || (lower_pos_index <= (nFluxuate+1)) );
        assert( (rougher_second_deriv[upper_pos_index] >= 0)   || ((upper_pos_index + nFluxuate + 1) >= end_channel) );
        assert( (rougher_second_deriv[upper_pos_index-1] <= 0) || (upper_pos_index <= (nFluxuate+1)) );

        //lower_pos_index and upper_pos_index are the first positive channels, so lets estimate,
        //  very crudely, where this crossing took place
        const float lower_crossing = lower_pos_index
          + (-rougher_second_deriv[lower_pos_index+1] / (rougher_second_deriv[lower_pos_index] - rougher_second_deriv[lower_pos_index+1]));
        const float upper_crossing = upper_pos_index
          + (rougher_second_deriv[upper_pos_index-1] / (rougher_second_deriv[upper_pos_index] - rougher_second_deriv[upper_pos_index-1]));

        const double lower_sigma_energy = data->energy_calibration()->energy_for_channel(lower_crossing);
        const double upper_sigma_energy = data->energy_calibration()->energy_for_channel(upper_crossing);

        const double rougher_sigma = 0.5*(upper_sigma_energy - lower_sigma_energy);


        double rougher_secondsum = 0.0;
        size_t min_rough_index = minbin;
        for( size_t index = lower_pos_index + 1; index < upper_pos_index; ++index )
        {
          if( rougher_second_deriv[index] < 0.0 )
            rougher_secondsum += rougher_second_deriv[index];
          if( rougher_second_deriv[index] < rougher_second_deriv[min_rough_index] )
            min_rough_index = index;
        }

        const double rougher_deriv_sigma = 0.5*( upper_crossing - lower_crossing );
        const double rougher_part = sqrt( 2.0 / ( boost::math::constants::pi<double>() *
                                         boost::math::constants::e<double>() ) )
                              / ( rougher_deriv_sigma * rougher_deriv_sigma );
        const double rougher_amplitude = -amp_fake_factor * rougher_secondsum / rougher_part;


        double rough_data_area = 0.0;
        for( size_t i = lower_pos_index+1; i <= upper_pos_index-1; ++i )
          rough_data_area += smoothed_data[i];
        const double rough_est_std_dev = sqrt( std::max(rough_data_area,1.0) );


        rougher_FOM = 0.68*rougher_amplitude / rough_est_std_dev;


        // We'll take the ROI to include the first negative bins, after 2nd deriv gos positive
        roi_begin_index = lower_pos_index;
        roi_end_index = upper_pos_index;

        assert( rougher_second_deriv[lower_pos_index] >= 0 || (lower_pos_index < (2*nFluxuate + 1) ) );
        assert( (rougher_second_deriv[upper_pos_index] >= 0) || (upper_pos_index > (end_channel - 2*nFluxuate - 1) ) );
        //double lower_data_counts = 0.0, upper_data_counts = 0.0, mid_data_counts = 0.0;
        //size_t tester = 0;
        for( ; (roi_begin_index > 0) && (rougher_second_deriv[roi_begin_index] > 0.0); --roi_begin_index )
        {
          //++tester;
          //lower_data_counts += spectrum[roi_begin_index];
          if( (roi_begin_index > 1) && (rougher_first_deriv[roi_begin_index] < 0) && (rougher_first_deriv[roi_begin_index-1] < 0) )
            break;
        }
        //assert( tester == (lower_pos_index - roi_begin_index) );
        //lower_data_counts /= (lower_pos_index - roi_begin_index);

        //tester = 0;
        for( ; ((roi_end_index+1) < end_channel) && (rougher_second_deriv[roi_end_index] > 0.0); ++roi_end_index )
        {
          //++tester;
          //upper_data_counts += spectrum[roi_end_index];
          if( ((roi_end_index + 1) < end_channel) && (rougher_first_deriv[roi_end_index] > 0.0) && (rougher_first_deriv[roi_end_index+1] > 0.0) )
            break;
        }
        //assert( tester == (roi_end_index - upper_pos_index) );
        //upper_data_counts /= (roi_end_index - upper_pos_index);

        //tester = 0;
        //for( size_t index = lower_pos_index + 1; index < upper_pos_index; ++index )
        //{
          //++tester;
          //mid_data_counts += spectrum[index];
        //}
        //assert( tester == (upper_pos_index - lower_pos_index - 1) );
        //mid_data_counts /= (upper_pos_index - lower_pos_index - 1);

        //cout << "midCPC: " << mid_data_counts << ", low: " << lower_data_counts << ", up: " << upper_data_counts << endl;

        //cout << "mean=" << mean << " -> sigma=" << sigma
        //<< ", rougher_sigma=" << rougher_sigma
        //<< ", rougher_amp=" << rougher_amplitude
        //<< ", rougher_FOM=" << rougher_FOM
        //<< ", rougher_ROI=[" << data->gamma_channel_lower(roi_begin_index) << ", " << data->gamma_channel_upper(roi_end_index) << "]"
        //<< ", energy(lower_pos_index)=" << data->gamma_channel_lower(lower_pos_index)
        //<< ", mean_second_deriv=" << rougher_second_deriv[data->find_gamma_channel(mean)]
        //<< endl;


        // Sum data in rougher +-sigma range
        // Sume data in ROI, excluding +-sigma range, and compare average counts per channel

        if( (rougher_sigma < sigma)
           && !IsNan(rougher_sigma)
           && !IsInf(rougher_sigma)
           && ((upper_crossing - lower_crossing) > 1.5f) )
        {
          sigma = rougher_sigma;
          if( rougher_second_deriv[min_rough_index] < 0.0 )
          {
            // TODO: we could use
            mean = data->gamma_channel_center(min_rough_index);
          }
        }

        // TODO: There is something to do using `rougher_second_deriv` to make sure there are at least a few more negative bins between `firstzero` and `secondzero` than posivite bins, to really get rid of larger-scale features
      }else //if( rougher_second_deriv[minbin] < 0 )
      {
        // If we are here, the `rougher_second_deriv` of the mean bin, is positive - indicating it is not the middle of the peak
        //
        // We tend to get here when:
        //  1) It is an insignificant peak
        //  2) Two HPGe peaks are close together, so that with the smoothing and everything, we cant resolve them
        //
        // We will still use the smoother second derivative to set the ROI bounds.
        // We will check the significance of the rougher second derivative, and if at least two channels next to
        //  each other pass a threshold of a positive value, then we its probably two peaks - in which case we will
        //  update the `mean` and `sigma` values using the rougher second derivative, for the largest of the peaks
        //  in the ROI.  We will NOT add more than one peak for this ROI - we will rely on later stages of peak fitting
        //  to pick this up.

        // We will use `second_deriv` to find the ROI extent - so thi
        //  TODO: Maybe use `rougher_first_deriv` to help with this.
        roi_begin_index = firstzero;
        roi_end_index = secondzero;

        //assert( second_deriv[firstzero] >= 0.0 ); // if sf is positive, I dont think this would be true
        //assert( second_deriv[secondzero] >= 0.0 ); // Not strictly true, especially if minbin is on edge of ROI

        // Make sure second_deriv[roi_begin_index] and second_deriv[roi_end_index] are positive,
        //  although they almost always should already be
        for( ; (roi_begin_index > 0) && (second_deriv[roi_begin_index] < 0.0); --roi_begin_index )
        {
        }

        for( ; ((roi_end_index+1) < end_channel) && (second_deriv[roi_end_index] < 0.0); ++roi_end_index )
        {
        }

        //cout << "1) roi_begin_index=" << data->gamma_channel_lower(roi_begin_index) << ", roi_end_index=" << data->gamma_channel_upper(roi_end_index) << endl;

        assert( second_deriv[roi_begin_index] >= 0.0 );
        assert( second_deriv[roi_end_index] >= 0.0 );

        // We are now in areas where second-derivative is positive - go until first negative bin,
        //  and this will define our ROI.
        for( ; (roi_begin_index > 0) && (second_deriv[roi_begin_index] > 0.0); --roi_begin_index )
        {
        }

        for( ; ((roi_end_index+1) < end_channel) && (second_deriv[roi_end_index] > 0.0); ++roi_end_index )
        {
        }

        if( debug_channel )
        {
          cout << "Looking at coarser (positive) second-deriv for mean=" << mean << " - minbin=" << minbin << endl;
          cout << "firstzero=" << data->gamma_channel_lower(firstzero) << ", secondzero=" << data->gamma_channel_upper(secondzero) << endl;

          cout << "Channel, Energy, 2nd, 2nd-rough, 2nd-rough-sqrt(var)" << endl;
          for( size_t i = roi_begin_index; i <= roi_end_index; ++i )
          {
            cout << i << ", " << data->gamma_channel_lower(i) << ", " << second_deriv[i] << ", " << rougher_second_deriv[i] << ", " << sqrt(rougher_second_deriv_var[i]);
            if( i == firstzero )
              cout << "  <-- firstzero";
            else if( i == secondzero )
              cout << "  <-- secondzero";
            cout << endl;
          }
          cout << endl;
        }

        // If the rougher second derivative being negative is really significant, we will update mean and sigma
        //  to match what it says (for the most significant peak in this ROI

        // Make lambda to give second-derivative significance for channels relative to `minbin`
        const auto second_deriv_sig = [minbin, &rougher_second_deriv, &rougher_second_deriv_var]( const int i ) -> double {
          if( (minbin < i) || ((minbin + i) >= rougher_second_deriv_var.size()) || (rougher_second_deriv_var[minbin + i] <= 0.0) )
            return 0.0;
          return rougher_second_deriv[minbin+i]/sqrt(rougher_second_deriv_var[minbin+i]);
        };

        int num_sig_pos = 0;
        const double second_deriv_sig_thresh = 1.0;  //This is totall arbitrary, and not vetted - just from a single ROI of a single file
        num_sig_pos += (second_deriv_sig(-1) > second_deriv_sig_thresh);
        num_sig_pos += (second_deriv_sig(0) > second_deriv_sig_thresh);
        num_sig_pos += (second_deriv_sig(1) > second_deriv_sig_thresh);

        const bool roi_center_is_dip = (num_sig_pos >= 2);

        if( roi_center_is_dip )
        {
          // We will verify we can see to peaks (means) inside the ROI, and if so, set `rougher_confident_multi_roi` to
          //  true, so we will know when we are making decisions.
          bool peak_lower = false, peak_upper = false;
          const double required_convex_limit = -1.5; // Totally arbitrarily chosen - not investigated - just something we are reasonably sure is the peak mea

          for( size_t i = roi_begin_index; !peak_lower && (i < minbin); ++i )
          {
            const double val = rougher_second_deriv[i];
            const double var = (rougher_second_deriv_var[i] > 0.0) ? sqrt(rougher_second_deriv_var[i]) : 1.0;
            peak_lower = ((val / var) < required_convex_limit);
          }

          for( size_t i = minbin + 1; !peak_upper && (i <= roi_end_index); ++i )
          {
            const double val = rougher_second_deriv[i];
            const double var = (rougher_second_deriv_var[i] > 0.0) ? sqrt(rougher_second_deriv_var[i]) : 1.0;
            peak_upper = ((val / var) < required_convex_limit);
          }

          rougher_confident_multi_roi = (peak_lower && peak_upper);
        }

        if( roi_center_is_dip )
        {
          //cout << "Will update for...mean=" << mean << endl;

          // Scan from `firstzero` to `secondzero` and find the most negative `rougher_second_deriv` - call this
          size_t min_2nd_deriv_index = firstzero;
          for( size_t i = firstzero + 1; i < secondzero; ++i )
          {
            if( rougher_second_deriv[i] < rougher_second_deriv[min_2nd_deriv_index] )
              min_2nd_deriv_index = i;
          }

          size_t rough_firstzero = firstzero, rough_secondzero = secondzero;
          for( size_t i = min_2nd_deriv_index + 1; i < secondzero; ++i )
          {
            if( std::signbit(rougher_second_deriv[i-1]) != std::signbit(rougher_second_deriv[i]) )
            {
              rough_secondzero = i;
              break;
            }
          }

          for( size_t i = min_2nd_deriv_index - 1; (i > 0) && (i > firstzero); ++i )
          {
            if( std::signbit(rougher_second_deriv[i]) != std::signbit(rougher_second_deriv[i+1]) )
            {
              rough_firstzero = i;
              break;
            }
          }

          const double rough_mean = data->gamma_channel_center(min_2nd_deriv_index);
          const double rough_sigma = 0.5*(data->gamma_channel_center(rough_secondzero)
                                             - data->gamma_channel_center(rough_firstzero));
          if( rough_sigma > 0.01 )
          {
            mean = rough_mean;
            sigma = rough_sigma;
          }
        }//if( num_sig_pos >= 2 )

        if( debug_channel )
          cout << "Coarser look updated mean=" << mean << ", sigma=" << sigma << endl;
      }//if( rougher_second_deriv can be used to better estimate sigma ) / else

      assert( roi_begin_index != roi_end_index );

      const float roi_start_energy = data->gamma_channel_lower( roi_begin_index );
      const float roi_end_energy = data->gamma_channel_upper( roi_end_index );

      const double deriv_sigma = 0.5*( secondzero - firstzero );
      const double part = sqrt( 2.0 / ( boost::math::constants::pi<double>() *
                                       boost::math::constants::e<double>() ) )
                            / ( deriv_sigma * deriv_sigma );
      const double amplitude = -amp_fake_factor * secondsum / part;


      // The following quantitiy looks to have some power, with it usually having a value ~0.1
      //  for legit peaks, and higher for non-legit.
      //  Not deeply looked into
      //const double min_to_area = second_deriv[minbin] / secondsum;
      //cout << "Energy=" << mean << " minval/sum=" << min_to_area << endl;


      double data_area = 0.0;
      for( size_t i = firstzero; i <= secondzero; ++i )
        data_area += smoothed_data[i];


      double est_std_dev = sqrt( std::max(data_area,1.0) );

      //In principle we would want to use the (true) continuums area to derive
      //  the est_std_dev from, but for practical purposes I think this can give
      //  us false positives fairly often

      const double figure_of_merit = 0.68*amplitude / est_std_dev;

      if( debug_channel )
      {
        cout << "Candidate mean=" << mean << ", est_std_dev=" << est_std_dev << ", amplitude=" << amplitude << ", figure_of_merit=" << figure_of_merit
        << ", threshold_FOM=" << threshold_FOM << ", more_scrutiny_FOM_threshold=" << settings.more_scrutiny_FOM_threshold
        << endl;
      }

      bool passed_higher_scrutiny = true;
      if( figure_of_merit < settings.more_scrutiny_FOM_threshold
         && (settings.more_scrutiny_FOM_threshold >= threshold_FOM) )
      {
        passed_higher_scrutiny = ( (rougher_scnd_drv < 0.0) && (rougher_FOM >= settings.more_scrutiny_coarser_FOM) );

        //rougher_amplitude

        //passed_higher_scrutiny
      }//if( figure_of_merit < settings.more_scrutiny_FOM_threshold )


      // For peaks in low-stat area, a decent way of separating noise peaks out, is to just use the
      //  ROI neighboring regions to estimate a straight line, and then check for at least one
      //  channel in the ROI, that goes above that line, by like ~4 sigma.
      // Right now, the amplitude below which we do this check is hard-coded, just because there
      //  are to many variables to loop over to optimize - after getting a better idea of optimal
      //  values of other variables, we'll add this one into the optimization, when fixing other.
      if( passed_higher_scrutiny && ((figure_of_merit < settings.more_scrutiny_FOM_threshold) || (amplitude < settings.amp_to_apply_line_test_below)) )
      {
        //Using the un-smoothed data - lets estimate the continuum using the area on
        //  either side of the ROI, then compute a chi2 for a straight line, and use this
        //  to help decide if a peak _could_ be useful
        const size_t cont_est_channels = 1 + roi_end_index - roi_begin_index;
        const size_t lower_edge_start = (cont_est_channels < roi_begin_index)
                                      ? (roi_begin_index - cont_est_channels) : 0;
        const size_t upper_edge_end = ((roi_end_index + cont_est_channels) < end_channel)
                                      ? (roi_end_index + cont_est_channels)
                                      : ((end_channel - cont_est_channels) > roi_end_index) ? (end_channel - cont_est_channels)
                                                                                            : std::min( roi_end_index + side_bins, end_channel - 1 );


        size_t lower_nchannel = 0, upper_nchannel = 0;
        double lower_cnts_per_chnl = 0.0, upper_cnts_per_chnl = 0.0;
        for( size_t i = lower_edge_start; i < roi_begin_index; ++i, ++lower_nchannel  )
          lower_cnts_per_chnl += spectrum[i];
        for( size_t i = roi_end_index + 1; i <= upper_edge_end; ++i, ++upper_nchannel )
          upper_cnts_per_chnl += spectrum[i];

        assert( lower_nchannel == (roi_begin_index - lower_edge_start) );
        assert( (upper_nchannel == (upper_edge_end - roi_end_index)) || (upper_edge_end < roi_end_index) ); //Second case can happen at the very end of the spectrum

        //if( mean > 640.0 && mean < 641 )
        //{
        //  cout << "lower_cnts_per_chnl=" << lower_cnts_per_chnl << ", upper_cnts_per_chnl=" << upper_cnts_per_chnl
        //  << ", lower_nchannel=" << lower_nchannel << ", upper_nchannel=" << upper_nchannel << endl;
        //}

        lower_cnts_per_chnl /= (lower_nchannel ? lower_nchannel : size_t(1));
        upper_cnts_per_chnl /= (upper_nchannel ? upper_nchannel : size_t(1));
        const double diff_per_chnl = (upper_cnts_per_chnl - lower_cnts_per_chnl) / (1 + roi_end_index - roi_begin_index);




        double chi2_dof = 0.0, max_chi2 = 0.0;
        // TODO: could/should probably shrink from checking entire ROI, +-1 sigma from mean
        for( size_t i = roi_begin_index; i <= roi_end_index; ++i )
        {
          const double pred_chnl_nts = lower_cnts_per_chnl + ((i-roi_begin_index) + 0.5) * diff_per_chnl;
          const double dc = spectrum[i] - pred_chnl_nts;

          //if( mean > 640.0 && mean < 641 )
          //{
          //  cout << "For energy " << mean << ", channel " << i << " has energy "
          //  << data->gamma_channel_lower(i) << " with counts "
          //  << spectrum[i] << " and predicted counts " << pred_chnl_nts << endl;
          //}

          const double uncert2 = std::max( (pred_chnl_nts > 0.0) ? pred_chnl_nts : static_cast<double>(spectrum[i]), 1.0 );
          const double val = dc*dc / uncert2;

          if( spectrum[i] > pred_chnl_nts )
            max_chi2 = std::max( max_chi2, val );

          chi2_dof += val;
        }
        chi2_dof /= (roi_end_index - roi_begin_index);

        if( PeakFitImprove::debug_printout )
        {
          //cout << "Straight-line Chi2 = " << chi2_dof << " and max_chi2=" << max_chi2 << " for energy=" << mean << endl;
        }

        //vector<double> cont_equation{0.0, 0.0};
        //PeakContinuum::eqn_from_offsets( roi_begin_index, roi_end_index, mean, data,
        //                                cont_est_channels, cont_est_channels,
        //                                cont_equation[1], cont_equation[0] );

        passed_higher_scrutiny = (max_chi2 > settings.more_scrutiny_min_dev_from_line);
      }//if( check Chi2 of region )


      if( (figure_of_merit < 5) && (rougher_scnd_drv >= 0.0) && !rougher_confident_multi_roi )
        passed_higher_scrutiny = false;

      if( debug_channel )
        cout << "  passed_higher_scrutiny=" << passed_higher_scrutiny << ", figure_of_merit=" << figure_of_merit << ", threshold_FOM=" << threshold_FOM << endl;

      //const double min_required_data = settings.min_counts_per_channel*2*deriv_sigma;

      if( (figure_of_merit > threshold_FOM) && passed_higher_scrutiny )
      {
        if( debug_channel || PeakFitImprove::debug_printout )
        {
           cout << "Accepted: energy=" << mean << ", FOM=" << figure_of_merit
                << ", amp=" << amplitude << ", FWHM=" << sigma*2.35482f
                << ", data_area=" << data_area << ", rougher_FOM=" << rougher_FOM
                << ", ROI=[" << roi_start_energy << ", " << roi_end_energy
                << endl;
        }

        PeakDef peak( mean, sigma, amplitude );
        peak.continuum()->setRange( roi_start_energy, roi_end_energy );

        vector<double> cont_equation{0.0, 0.0};
        PeakContinuum::eqn_from_offsets( roi_begin_index, roi_end_index, mean, data,
                                        settings.num_smooth_side_channels,
                                        settings.num_smooth_side_channels,
                                        cont_equation[1], cont_equation[0] );
        peak.continuum()->setType( PeakContinuum::OffsetType::Linear );
        peak.continuum()->setParameters( mean, cont_equation, {} );

        result_peaks.push_back( std::move(peak) );
      }else
      {
        if( debug_channel || PeakFitImprove::debug_printout )
       {
         cout << "Rejected: energy=" << mean << ", FOM=" << figure_of_merit
         << ", amp=" << amplitude << ", FWHM=" << sigma*2.35482f
         << ", data_area=" << data_area  << ", rougher_FOM=" << rougher_FOM
         << ", ROI=[" << roi_start_energy << ", " << roi_end_energy
         << endl;
       }
      }//if( region we were just in passed_threshold )

      secondsum = 0.0;
      minval = 9999999999.9f;
      minbin = secondzero = firstzero = 0;
    }else
    {
      bool belowzero = true, goingnegative = true, abovezero = true;
      for( size_t i = 0; i < nFluxuate; ++i )
      {
        if( (channel+i+1) < end_channel )
          goingnegative &= (second_deriv[channel+i+1] < 0.0f);
        if( channel >= i )
        {
          belowzero &= (second_deriv[channel-i] <= 0.0f);
          abovezero &= (second_deriv[channel-i] > 0.0f);
        }
      }//for( size_t i = 0; i < nFluxuate; ++i )

      if( channel && !firstzero && goingnegative )
      {
        firstzero = channel;
        minbin = channel;
        minval = secondDeriv;

        for( size_t i = 1; i < nFluxuate; ++i )
          if( channel >= i )
            secondsum += second_deriv[channel-i];
      }else if( secondSumPositive )
      {
        secondsum = 0.0;
        minval = 9999999999.9f;
        minbin = secondzero = firstzero = 0;
      }

      if( firstzero > 0 )
      {
        secondsum += secondDeriv;

        if( secondDeriv < minval )
        {
          minbin = channel;
          minval = secondDeriv;
        }
      }//if( firstzero > 0 )
    }//if( we are out of region of interest) / else( in region of )
  }//for( loop over bins )

  return result_peaks;
}//find_candidate_peaks(...)



//<score, num_peaks_found, num_possibly_accepted_peaks_not_found, num_extra_peaks>
tuple<double,size_t,size_t,size_t,size_t,size_t> eval_candidate_settings( const FindCandidateSettings settings, const vector<DataSrcInfo> &input_srcs, const bool write_n42 )
{
  double sum_score = 0.0;
  size_t num_possibly_accepted_peaks_not_found = 0, num_def_wanted_not_found = 0;
  size_t num_extra_peaks = 0, num_peaks_found = 0;
  size_t num_def_wanted_peaks_found = 0;

  for( const DataSrcInfo &info : input_srcs )
  {
    shared_ptr<const SpecUtils::Measurement> src_spectrum = info.src_info.src_spectra.front();
    const vector<PeakDef> peaks = CandidatePeak_GA::find_candidate_peaks( src_spectrum, 0, 0, settings );

    vector<tuple<float,float,float>> detected_peaks; //{mean, sigma, amplitude}
    for( const PeakDef &p : peaks )
      detected_peaks.emplace_back( p.mean(), p.sigma(), p.amplitude() );

    const vector<tuple<float,float,float>> orig_peak_candidates = detected_peaks;

    double score = 0.0;
    size_t num_detected_expected = 0, num_detected_not_expected = 0;
    size_t num_possibly_accepted_but_not_detected = 0, num_def_wanted_but_not_detected = 0;
    size_t num_def_wanted_detected = 0;

    // First, go through expected peaks, and match up to candidates
    vector<ExpectedPhotopeakInfo> possibly_expected_but_not_detected, expected_and_was_detected, def_expected_but_not_detected;
    possibly_expected_but_not_detected.reserve( info.expected_photopeaks.size() );
    def_expected_but_not_detected.reserve( info.expected_photopeaks.size() );
    expected_and_was_detected.reserve( info.expected_photopeaks.size() );

    vector<bool> detected_matched_to_expected( detected_peaks.size(), false );

    for( size_t expected_index = 0; expected_index < info.expected_photopeaks.size(); ++expected_index )
    {
      const ExpectedPhotopeakInfo &expected = info.expected_photopeaks[expected_index];

      vector<pair<tuple<float,float,float>,size_t>> detected_matching_expected; //{mean, sigma, amplitude}
      for( size_t det_index = 0; det_index < detected_peaks.size(); ++det_index )
      {
        const tuple<float,float,float> &det_peak = detected_peaks[det_index];
        const float mean = get<0>(det_peak);
        //const float sigma = get<1>(det_peak);
        //const float amp = get<1>(det_peak);

        if( (mean > expected.roi_lower) && (mean < expected.roi_upper) )
        {
          detected_matching_expected.push_back( make_pair(det_peak,det_index) );
          detected_matched_to_expected[det_index] = true;
        }
      }//for( size_t det_index = 0; det_index < detected_peaks.size(); ++i )

      if( detected_matching_expected.empty() )
      {
        num_possibly_accepted_but_not_detected += 1;

        if( (expected.nsigma_over_background > JudgmentFactors::def_want_nsigma)
           && (expected.peak_area > JudgmentFactors::min_def_wanted_counts) )
        {
          num_def_wanted_but_not_detected += 1;
          def_expected_but_not_detected.push_back( expected );
          if( PeakFitImprove::debug_printout )
            cerr << "Def. wanted m=" << expected.effective_energy << ", nsigma=" << expected.nsigma_over_background << ", but didnt find." << endl;
        }else
        {
          possibly_expected_but_not_detected.push_back( expected );
        }

        /*
         if( expected.nsigma_over_background < JudgmentFactors::lower_want_nsigma )
         score -= 0; //No punishment
         else if( expected.nsigma_over_background < JudgmentFactors::def_want_nsigma )
         score -= ( (JudgmentFactors::def_want_nsigma - expected.nsigma_over_background) / (expected.nsigma_over_background - JudgmentFactors::lower_want_nsigma) );
         else
         score -= 1;
         */
      }else
      {
        // we got a match
        num_detected_expected += 1;
        expected_and_was_detected.push_back( expected );

        if( (expected.nsigma_over_background > JudgmentFactors::def_want_nsigma)
           && (expected.peak_area > JudgmentFactors::min_def_wanted_counts) )
        {
          num_def_wanted_detected += 1;
        }

        if( expected.nsigma_over_background > JudgmentFactors::def_want_nsigma )
        {
          score += 1.0;
        }else if( expected.nsigma_over_background > JudgmentFactors::lower_want_nsigma )
        {
          const double amount_short = JudgmentFactors::def_want_nsigma - expected.nsigma_over_background;
          const double fraction_short = amount_short / (JudgmentFactors::def_want_nsigma - JudgmentFactors::lower_want_nsigma);
          assert( (fraction_short >= 0.0) && (fraction_short <= 1.0) );
          score += (1 - fraction_short);
        }else
        {
          score += 0.0; //No punishment, but no reward.
        }
      }//if( detected_matching_expected.empty() ) / else
    }//for( loop over expected_photopeaks )

    assert( expected_and_was_detected.size() == num_detected_expected );

    vector<tuple<float,float,float>> detected_expected, detected_not_expected;
    detected_expected.reserve( info.expected_photopeaks.size() );
    detected_not_expected.reserve( detected_peaks.size() );

    for( size_t i = 0; i < detected_matched_to_expected.size(); ++i )
    {
      if( detected_matched_to_expected[i] )
        detected_expected.push_back( detected_peaks[i] );
      else
        detected_not_expected.push_back( detected_peaks[i] );
    }//for( size_t i = 0; i < detected_matched_to_expected.size(); ++i )

    num_detected_not_expected = detected_not_expected.size();

    // The 511 Annih. line isnt in our truth line, unless the source makes them, so lets
    //  remove this one from the score
    for( const auto &p : detected_not_expected )
    {
      if( (get<0>(p) > 508) && (get<0>(p) < 514) )
      {
        num_detected_not_expected -= 1;
        break;
      }
    }


    score -= JudgmentFactors::found_extra_punishment * num_detected_not_expected;

    /*
    cout << "For " << info.src_info.file_base_path
    //<< info.src_info.src_name
    //<< "/" << info.detector_name
    //<< "/" << info.location_name
    //<< "/" << info.live_time_name
    << ":\n"
    << "\tnum_detected_expected=" << num_detected_expected
    << ", num_detected_not_expected=" << num_detected_not_expected
    << ", num_possibly_accepted_but_not_detected=" << num_possibly_accepted_but_not_detected
    << ", score=" << score << endl;
    */

    sum_score += score;
    num_peaks_found += num_detected_expected;
    num_possibly_accepted_peaks_not_found += num_possibly_accepted_but_not_detected;
    num_def_wanted_not_found += num_def_wanted_but_not_detected;
    num_def_wanted_peaks_found += num_def_wanted_detected;
    num_extra_peaks += num_detected_not_expected;

    if( write_n42 )
    {
      //const string src_dir = SpecUtils::parent_path( info.src_info.file_base_path );

      string outdir = "output_n42";
      if( !SpecUtils::is_directory(outdir) && !SpecUtils::create_directory(outdir) )
        cerr << "Failed to create directory '" << outdir << "'" << endl;

      outdir = SpecUtils::append_path( outdir, info.detector_name );
      if( !SpecUtils::is_directory(outdir) && SpecUtils::create_directory(outdir) )
        cerr << "Failed to create directory '" << outdir << "'" << endl;

      outdir = SpecUtils::append_path( outdir, info.location_name );
      if( !SpecUtils::is_directory(outdir) && !SpecUtils::create_directory(outdir) )
        cerr << "Failed to create directory '" << outdir << "'" << endl;

      outdir = SpecUtils::append_path( outdir, info.live_time_name );
      if( !SpecUtils::is_directory(outdir) && !SpecUtils::create_directory(outdir) )
        cerr << "Failed to create directory '" << outdir << "'" << endl;

      const string out_n42 = SpecUtils::append_path( outdir, info.src_info.src_name ) + "_peak_candidates.n42";

      SpecMeas output;

      output.add_remark( "Failed to find " + std::to_string(num_possibly_accepted_but_not_detected) + " peaks that are possibly accepted." );
      output.add_remark( "Found " + std::to_string(num_detected_expected) + " peaks that were expected." );
      output.add_remark( "Found " + std::to_string(num_detected_not_expected) + " peaks that were NOT expected." );
      output.add_remark( settings.print("settings") );
      output.set_instrument_model( info.detector_name );
      if( SpecUtils::icontains(info.detector_name, "Detective" ) )
        output.set_manufacturer( "ORTEC" );
      else if( SpecUtils::icontains(info.detector_name, "Falcon 5000" ) )
        output.set_manufacturer( "Canberra" );
      else if( SpecUtils::icontains(info.detector_name, "Fulcrum" ) )
        output.set_manufacturer( "PHDS" );

      output.set_measurement_location_name( info.location_name );

      shared_ptr<SpecUtils::Measurement> out_with_cand = make_shared<SpecUtils::Measurement>( *src_spectrum );
      out_with_cand->set_sample_number( 1 );

      const string title = info.detector_name + "/" + info.src_info.src_name;
      out_with_cand->set_title( title + " + candidate peaks" );
      auto now = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());
      out_with_cand->set_start_time( now );
      output.add_measurement( out_with_cand, false );

      deque<shared_ptr<const PeakDef>> peaks;

      vector<tuple<float,float,float>> all_candidates = detected_expected;
      all_candidates.insert( end(all_candidates), begin(detected_not_expected), end(detected_not_expected) );
      for( size_t peak_index = 0; peak_index < all_candidates.size(); ++peak_index )
      {
        const tuple<float,float,float> &p = all_candidates[peak_index]; //{mean, sigma, amplitude}
        const float mean = get<0>(p);
        const float sigma = get<1>(p);
        const float amp = get<2>(p);

        const bool is_expected = (peak_index < detected_expected.size());

        auto peak = make_shared<PeakDef>( mean, sigma, amp );
        peak->setFitFor( PeakDef::CoefficientType::Mean, false );
        peak->setFitFor( PeakDef::CoefficientType::Sigma, false );
        peak->setFitFor( PeakDef::CoefficientType::GaussAmplitude, false );
        peak->continuum()->setType( PeakContinuum::OffsetType::Linear );
        peak->continuum()->setRange( mean - 3*sigma, mean + 3*sigma );
        peak->continuum()->calc_linear_continuum_eqn( src_spectrum, mean, mean - 4*sigma, mean + 4*sigma, 5, 5 );
        peak->setLineColor( is_expected ? Wt::GlobalColor::darkGreen : Wt::GlobalColor::red );

        if( is_expected )
        {
          //Put truth-level info to user-label field
          string info_str;
          for( const ExpectedPhotopeakInfo &exp_info : info.expected_photopeaks )
          {
            if( mean > exp_info.roi_lower && mean < exp_info.roi_upper )
            {
              if( !info_str.empty() )
                info_str += ".\n";

              info_str += "E: A=" + SpecUtils::printCompact(exp_info.peak_area, 3)
              + ", W=" + SpecUtils::printCompact(exp_info.gamma_lines.front().fwhm, 3)
              + ", #s=" + SpecUtils::printCompact(exp_info.nsigma_over_background, 3);
            }
          }

          peak->setUserLabel( info_str );
        }else
        {
          peak->setUserLabel( "Not Expected" );
        }

        peaks.push_back( peak );
      }//for( const tuple<float,float,float> &p : all_candidates )


      vector<ExpectedPhotopeakInfo> missing_peaks = def_expected_but_not_detected;
      missing_peaks.insert( end(missing_peaks), begin(possibly_expected_but_not_detected), end(possibly_expected_but_not_detected) );

      for( size_t missing_index = 0; missing_index < missing_peaks.size(); ++missing_index )
      {
        const ExpectedPhotopeakInfo &missing = missing_peaks[missing_index];
        const bool def_wanted = (missing_index < def_expected_but_not_detected.size());

        const double mean = missing.effective_energy;
        const double sigma = missing.gamma_lines.front().fwhm/2.35482f;
        auto peak = make_shared<PeakDef>( mean, sigma, missing.peak_area );
        peak->setFitFor( PeakDef::CoefficientType::Mean, false );
        peak->setFitFor( PeakDef::CoefficientType::Sigma, false );
        peak->setFitFor( PeakDef::CoefficientType::GaussAmplitude, false );
        peak->continuum()->setType( PeakContinuum::OffsetType::Linear );
        peak->continuum()->setRange( mean - 3*sigma, mean + 3*sigma );
        peak->continuum()->calc_linear_continuum_eqn( src_spectrum, mean, mean - 4*sigma, mean + 4*sigma, 5, 5 );
        const Wt::WColor orange(252, 94, 3), yellow(168, 150, 50);
        peak->setLineColor( def_wanted ? orange : yellow );
        string label = def_wanted ? "Not det. - wanted" : "Not det. - poss.";
        label += ", N_Sigma=" + SpecUtils::printCompact(missing.nsigma_over_background, 3);
        peak->setUserLabel( label );

        peaks.push_back( peak );
      }//for( const ExpectedPhotopeakInfo &missing : not_detected )





      output.setPeaks( peaks, {1} );

      {// begin add second derivative to N42
        auto second_deriv = make_shared<vector<float>>();
        const int side_bins = settings.num_smooth_side_channels;
        const int poly_order = settings.smooth_polynomial_order;
        smoothSpectrum( src_spectrum, side_bins, poly_order, 2, *second_deriv );

        shared_ptr<SpecUtils::Measurement> second_deriv_meas = make_shared<SpecUtils::Measurement>();
        second_deriv_meas->set_gamma_counts( second_deriv, 1.0, 1.0 );
        second_deriv_meas->set_energy_calibration( src_spectrum->energy_calibration() );
        second_deriv_meas->set_sample_number( 2 );
        second_deriv_meas->set_title( title + " + smoothed second derivative." );
        output.add_measurement( second_deriv_meas, true );

        deque<shared_ptr<const PeakDef>> candidates;
        for( const tuple<float,float,float> &p : orig_peak_candidates )
        {
          auto peak = make_shared<PeakDef>( get<0>(p), get<1>(p), get<2>(p) );
          peak->continuum()->setType( PeakContinuum::OffsetType::Constant );
          peak->continuum()->setPolynomialCoef(0, 0.0);
          candidates.push_back( peak );
        }

        output.setPeaks( candidates, {2} );
      }// end add second derivative to N42


      {
        // I think this will be good for estimating sigma and area
        const int side_bins = 3;
        const int poly_order = 3;
        vector<float> rougher_second_deriv;
        smoothSpectrum( src_spectrum, side_bins, poly_order, 2, rougher_second_deriv );
        shared_ptr<SpecUtils::Measurement> rough_second_deriv_meas = make_shared<SpecUtils::Measurement>();
        auto rough_second_deriv_counts = make_shared<vector<float>>( rougher_second_deriv );
        rough_second_deriv_meas->set_gamma_counts( rough_second_deriv_counts, 1.0, 1.0 );
        rough_second_deriv_meas->set_energy_calibration( src_spectrum->energy_calibration() );
        rough_second_deriv_meas->set_sample_number( 3 );
        rough_second_deriv_meas->set_title( title + " + rougher smoothed second derivative." );
        output.add_measurement( rough_second_deriv_meas, true );
      }


      {
        const int side_bins = settings.num_smooth_side_channels;
        const int poly_order = 1;
        vector<float> smoothed_data;
        smoothSpectrum( src_spectrum, side_bins, poly_order, 0, smoothed_data );
        shared_ptr<SpecUtils::Measurement> meas = make_shared<SpecUtils::Measurement>();
        auto rough_second_deriv_counts = make_shared<vector<float>>( smoothed_data );
        meas->set_gamma_counts( rough_second_deriv_counts, out_with_cand->live_time(), out_with_cand->real_time() );
        meas->set_energy_calibration( src_spectrum->energy_calibration() );
        meas->set_sample_number( 4 );
        meas->set_title( title + " + smoothed data." );
        output.add_measurement( meas, true );
      }

      //ofstream outstrm( out_n42.c_str(), ios::out | ios::binary );
      //if( !outstrm )
      //  cerr << "Failed to open '" << out_n42 << "'" << endl;
      //output.write( outstrm, output.sample_numbers(), output.detector_names(), SpecUtils::SaveSpectrumAsType::N42_2012 );

      output.save2012N42File( out_n42, [=](){
        cerr << "Failed to write '" << out_n42 << "'" << endl;
      });

      //cout << "Wrote '" << out_n42 << endl;
    }//if( write_n42 )
  }//for( const DataSrcInfo &info : input_srcs )

  //cout << "Avrg score: " << sum_score/input_srcs.size() << endl;

  return {sum_score / input_srcs.size(), num_peaks_found, num_def_wanted_not_found, num_def_wanted_peaks_found, num_possibly_accepted_peaks_not_found, num_extra_peaks};
};//eval_candidate_settings lambda


string CandidatePeakSolution::to_string( const string &separator ) const
{
  return
  string("num_smooth_side_channels: ") + std::to_string(num_smooth_side_channels)
  + separator + "smooth_polynomial_order: " + std::to_string(smooth_polynomial_order)
  + separator + "threshold_FOM: " + std::to_string(threshold_FOM)
  + separator + "more_scrutiny_FOM_threshold_delta: " + std::to_string(more_scrutiny_FOM_threshold_delta)
  + separator + "pos_sum_threshold_sf: " + std::to_string(pos_sum_threshold_sf)
  + separator + "num_chan_fluctuate: " + std::to_string(num_chan_fluctuate)
  + separator + "more_scrutiny_coarser_FOM_delta: " + std::to_string(more_scrutiny_coarser_FOM_delta)
  + separator + "more_scrutiny_min_dev_from_line: " + std::to_string(more_scrutiny_min_dev_from_line)
  + separator + "amp_to_apply_line_test_below: " + std::to_string(amp_to_apply_line_test_below)
  ;
}


void init_genes(CandidatePeakSolution& p,const std::function<double(void)> &rnd01)
{
  // rnd01() gives a random number in 0~1
  p.num_smooth_side_channels          = 6+7*rnd01(); //8;                   //2+13*rnd01();
  p.smooth_polynomial_order           = 2;                   //2+1*rnd01();
  p.threshold_FOM                     = 0.75+1.25*rnd01(); //0.75 + 0.75*rnd01(); //0.8+2.7*rnd01();
  p.more_scrutiny_FOM_threshold_delta = -0.1+2.1*rnd01(); //-0.1 + 1.2*rnd01();  //-0.5+4*rnd01();
  p.pos_sum_threshold_sf              = -0.15+0.3*rnd01(); // -0.1 + 0.11*rnd01();   //-0.1+0.2*rnd01();
  p.num_chan_fluctuate                = 1;                   //1+3*rnd01();
  p.more_scrutiny_coarser_FOM_delta   = -0.1+4.1*rnd01(); //2 + 4*rnd01();       //-0.1+5.1*rnd01();
  p.more_scrutiny_min_dev_from_line   = 0 + 7*rnd01(); //0 + 10*rnd01();      //0.0+10*rnd01();
  p.amp_to_apply_line_test_below      = 0.0+100*rnd01(); //0 + 75*rnd01();      //0.0+120*rnd01();
}

FindCandidateSettings genes_to_settings( const CandidatePeakSolution &p )
{
  FindCandidateSettings settings;

  settings.num_smooth_side_channels = p.num_smooth_side_channels;
  settings.smooth_polynomial_order = p.smooth_polynomial_order;
  settings.threshold_FOM = p.threshold_FOM;
  settings.more_scrutiny_FOM_threshold = p.threshold_FOM + p.more_scrutiny_FOM_threshold_delta;
  settings.pos_sum_threshold_sf = p.pos_sum_threshold_sf;
  settings.num_chan_fluctuate = p.num_chan_fluctuate;
  settings.more_scrutiny_coarser_FOM = p.threshold_FOM + p.more_scrutiny_coarser_FOM_delta;
  settings.more_scrutiny_min_dev_from_line = p.more_scrutiny_min_dev_from_line;
  settings.amp_to_apply_line_test_below = p.amp_to_apply_line_test_below;

  return settings;
}

bool eval_solution( const CandidatePeakSolution &p, CandidatePeakCost &c )
{
  const FindCandidateSettings settings = genes_to_settings( p );

  assert(ns_ga_eval_fcn);

  try
  {
    c.objective1 = -1.0 * ns_ga_eval_fcn( settings );
  }catch( std::exception & )
  {
    return false; //reject solution
  }

  return true; // solution is accepted
}


CandidatePeakSolution mutate(
    const CandidatePeakSolution& X_base,
    const std::function<double(void)> &rnd01,
    double shrink_scale)
{
  CandidatePeakSolution X_new;
  const double mu = 0.2*shrink_scale; // mutation radius (adjustable)

  const double mutate_threshold = PeakFitImprove::sm_ga_mutate_threshold;

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

    if( rnd01() > mutate_threshold )
      X_new.num_smooth_side_channels += shrink_scale*(rnd01()-rnd01()); //not multiplying by `mu`, because we can get stuck in a single int
    //in_range=in_range&&(X_new.num_smooth_side_channels>=2 && X_new.num_smooth_side_channels<15);
    if( rnd01() > mutate_threshold )
      X_new.num_smooth_side_channels = std::max( X_new.num_smooth_side_channels, 6 );
    if( rnd01() > mutate_threshold )
    {
      X_new.num_smooth_side_channels = std::min( X_new.num_smooth_side_channels, 13 );
      in_range=in_range&&(X_new.num_smooth_side_channels>=6 && X_new.num_smooth_side_channels<=13);
    }

    //X_new.smooth_polynomial_order+=mu*(rnd01()-rnd01()); //This is an int we could get stuck in...
    //in_range=in_range&&(X_new.smooth_polynomial_order>=2 && X_new.smooth_polynomial_order<3);

    if( rnd01() > mutate_threshold )
    {
      X_new.threshold_FOM += mu * (rnd01() - rnd01());
      //in_range=in_range&&(X_new.threshold_FOM>=0.75 && X_new.threshold_FOM<3.5);
      in_range=in_range&&(X_new.threshold_FOM>=0.75 && X_new.threshold_FOM<=2);
    }

    if( rnd01() > mutate_threshold )
    {
      X_new.more_scrutiny_FOM_threshold_delta += mu*(rnd01() - rnd01());
      //in_range=in_range&&(X_new.more_scrutiny_FOM_threshold_delta>=-0.5 && X_new.more_scrutiny_FOM_threshold_delta<3.5);
      in_range=in_range&&(X_new.more_scrutiny_FOM_threshold_delta>=-0.1 && X_new.more_scrutiny_FOM_threshold_delta<=2);
    }

    if( rnd01() > mutate_threshold )
    {
      X_new.pos_sum_threshold_sf += 0.05 * mu * (rnd01() - rnd01()); //mu*(rnd01() - rnd01());
      //in_range=in_range&&(X_new.pos_sum_threshold_sf>=-0.1 && X_new.pos_sum_threshold_sf<0.1);
      in_range=in_range&&(X_new.pos_sum_threshold_sf>=-0.15 && X_new.pos_sum_threshold_sf<=0.15);
    }

    //X_new.num_chan_fluctuate+=mu*(rnd01()-rnd01());  //THis is an int we could get stuck in...
    //in_range=in_range&&(X_new.num_chan_fluctuate>=1 && X_new.num_chan_fluctuate<4);

    if( rnd01() > mutate_threshold )
    {
      X_new.more_scrutiny_coarser_FOM_delta += mu*(rnd01() - rnd01());
      //in_range=in_range&&(X_new.more_scrutiny_coarser_FOM_delta>=-0.1 && X_new.more_scrutiny_coarser_FOM_delta<5);
      in_range=in_range&&(X_new.more_scrutiny_coarser_FOM_delta>=-0.1 && X_new.more_scrutiny_coarser_FOM_delta<4);
    }

    if( rnd01() > mutate_threshold )
    {
      X_new.more_scrutiny_min_dev_from_line += 2 * mu*(rnd01()-rnd01());
      //in_range=in_range&&(X_new.more_scrutiny_min_dev_from_line>=0.0 && X_new.more_scrutiny_min_dev_from_line<10.0);
      in_range=in_range&&(X_new.more_scrutiny_min_dev_from_line>=0.0 && X_new.more_scrutiny_min_dev_from_line<7.0);
    }

    if( rnd01() > mutate_threshold )
      X_new.amp_to_apply_line_test_below += 5 * mu*(rnd01()-rnd01());
    if( rnd01() > mutate_threshold )
      X_new.amp_to_apply_line_test_below = std::max( X_new.amp_to_apply_line_test_below, 0 );
    if( rnd01() > mutate_threshold )
    {
      X_new.amp_to_apply_line_test_below = std::min( X_new.amp_to_apply_line_test_below, 100 );
      in_range=in_range&&(X_new.amp_to_apply_line_test_below>=0.0 && X_new.amp_to_apply_line_test_below<=100);
      //in_range=in_range&&(X_new.amp_to_apply_line_test_below>=0.0 && X_new.amp_to_apply_line_test_below<=75);
    }
  } while(!in_range);
  return X_new;
}


CandidatePeakSolution crossover(
    const CandidatePeakSolution& X1,
    const CandidatePeakSolution& X2,
    const std::function<double(void)> &rnd01)
{
  const double crossover_threshold = PeakFitImprove::sm_ga_crossover_threshold;


  CandidatePeakSolution X_new = X1;


  double r;
  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.num_smooth_side_channels=r*X1.num_smooth_side_channels+(1.0-r)*X2.num_smooth_side_channels;
    //X_new.num_smooth_side_channels = X1.num_smooth_side_channels;
  }else if( rnd01() < 0.5 )
  {
    X_new.num_smooth_side_channels = X2.num_smooth_side_channels;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    //X_new.smooth_polynomial_order=r*X1.smooth_polynomial_order+(1.0-r)*X2.smooth_polynomial_order;
    X_new.smooth_polynomial_order = X1.smooth_polynomial_order;
  }else if( rnd01() < 0.5 )
  {
    X_new.smooth_polynomial_order = X2.smooth_polynomial_order;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.threshold_FOM=r*X1.threshold_FOM+(1.0-r)*X2.threshold_FOM;
  }else if( rnd01() < 0.5 )
  {
    X_new.threshold_FOM = X2.threshold_FOM;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.more_scrutiny_FOM_threshold_delta=r*X1.more_scrutiny_FOM_threshold_delta+(1.0-r)*X2.more_scrutiny_FOM_threshold_delta;
  }else if( rnd01() < 0.5 )
  {
    X_new.more_scrutiny_FOM_threshold_delta = X2.more_scrutiny_FOM_threshold_delta;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.pos_sum_threshold_sf=r*X1.pos_sum_threshold_sf+(1.0-r)*X2.pos_sum_threshold_sf;
  }else if( rnd01() < 0.5 )
  {
    X_new.pos_sum_threshold_sf = X2.pos_sum_threshold_sf;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    //X_new.num_chan_fluctuate=r*X1.num_chan_fluctuate+(1.0-r)*X2.num_chan_fluctuate;
    X_new.num_chan_fluctuate = X1.num_chan_fluctuate;
  }else if( rnd01() < 0.5 )
  {
    X_new.num_chan_fluctuate = X2.num_chan_fluctuate;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.more_scrutiny_coarser_FOM_delta=r*X1.more_scrutiny_coarser_FOM_delta+(1.0-r)*X2.more_scrutiny_coarser_FOM_delta;
  }else if( rnd01() < 0.5 )
  {
    X_new.more_scrutiny_coarser_FOM_delta = X2.more_scrutiny_coarser_FOM_delta;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.more_scrutiny_min_dev_from_line=r*X1.more_scrutiny_min_dev_from_line+(1.0-r)*X2.more_scrutiny_min_dev_from_line;
  }else if( rnd01() < 0.5 )
  {
    X_new.more_scrutiny_min_dev_from_line = X2.more_scrutiny_min_dev_from_line;
  }

  if( rnd01() < crossover_threshold )
  {
    r=rnd01();
    X_new.amp_to_apply_line_test_below=r*X1.amp_to_apply_line_test_below+(1.0-r)*X2.amp_to_apply_line_test_below;
  }else if( rnd01() < 0.5 )
  {
    X_new.amp_to_apply_line_test_below = X2.amp_to_apply_line_test_below;
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


void SO_report_generation(
    int generation_number,
    const EA::GenerationType<CandidatePeakSolution,CandidatePeakCost> &last_generation,
    const CandidatePeakSolution& best_genes)
{
  bool best_yet = false;
  if( !sm_set_best_genes || (last_generation.best_total_cost < sm_best_total_cost) )
  {
    best_yet = true;
    sm_set_best_genes = true;
    sm_best_genes = best_genes;
    sm_best_total_cost = last_generation.best_total_cost;
  }

  cout <<"Generation ["<<generation_number<<"], "
  <<"Best="<<last_generation.best_total_cost << ", "
  <<"Average="<<last_generation.average_cost << ", "
  <<"Best genes: {\n\t" <<best_genes.to_string("\n\t")  << "\n}\n"
  <<"Exe_time="<<last_generation.exe_time
  << endl << endl;

  sm_output_file
  <<generation_number<<"\t"
  <<last_generation.average_cost<<"\t"
  <<last_generation.best_total_cost<<"\t"
  << "{" << best_genes.to_string(", ") << "}\n\n";
}


FindCandidateSettings do_ga_eval( std::function<double(const FindCandidateSettings &)> ga_eval_fcn )
{
  assert( !sm_has_been_called );
  if( sm_has_been_called )
  {
    cerr << "You should only call CandidatePeak_GA::do_ga_eval(...) once per program execution!!!" << endl;
    exit(1);
  }

  sm_has_been_called = true;

  assert( !!ga_eval_fcn );
  if( !ga_eval_fcn )
    throw runtime_error( "Invalid eval function passed in." );

  ns_ga_eval_fcn = ga_eval_fcn;


  sm_output_file.open("results.txt");
  sm_output_file<<"step"<<"\t"<<"cost_avg"<<"\t"<<"cost_best"<<"\t"<<"solution_best"<<"\n";

  EA::Chronometer timer;
  timer.tic();

  GA_Type ga_obj;
  ga_obj.problem_mode=EA::GA_MODE::SOGA; //Single objective genetic algorithm
  ga_obj.multi_threading=true;
  ga_obj.idle_delay_us=1; // switch between threads quickly
  ga_obj.dynamic_threading = true; //If false,  thread responsibilities are divided at the beginning
  ga_obj.verbose=false;
  ga_obj.population = static_cast<unsigned int>(PeakFitImprove::sm_ga_population);
  ga_obj.generation_max = static_cast<int>(PeakFitImprove::sm_ga_generation_max);
  ga_obj.calculate_SO_total_fitness=calculate_SO_total_fitness;
  ga_obj.init_genes=init_genes;
  ga_obj.eval_solution=eval_solution;
  ga_obj.mutate=mutate;
  ga_obj.crossover=crossover;
  ga_obj.SO_report_generation=SO_report_generation;
  ga_obj.crossover_fraction=PeakFitImprove::sm_ga_crossover_fraction;
  ga_obj.mutation_rate=PeakFitImprove::sm_ga_mutation_rate;
  ga_obj.best_stall_max = static_cast<int>(PeakFitImprove::sm_ga_best_stall_max);
  ga_obj.elite_count = static_cast<int>(PeakFitImprove::sm_ga_elite_count);
  ga_obj.N_threads = static_cast<int>(PeakFitImprove::sm_num_optimization_threads);
  EA::StopReason stop_reason = ga_obj.solve();

  cout << "Stop reason was: " << ga_obj.stop_reason_to_string( stop_reason) << endl;
  cout << "The problem is optimized in " << timer.toc() << " seconds." << endl;

  sm_output_file.close();

  return genes_to_settings( sm_best_genes );
}
}//namespace CandidatePeak_GA
