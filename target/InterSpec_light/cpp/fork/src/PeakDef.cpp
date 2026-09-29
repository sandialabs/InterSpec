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

#include <map>
#include <regex>
#include <cfloat>
#include <memory>
#include <numeric>
#include <iostream>


#include "nlohmann/json.hpp"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"


#include "InterSpec/PeakDef.h"
#include "InterSpec/LightMath.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/PeakDists.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PhysicalUnits.h"

#include "InterSpec/PeakDists_imp.hpp"

using namespace std;
using SpecUtils::Measurement;

/** Peak XML version changes:
 20260120: Added GaussPlusBortel and DoubleBortel skew models. Renamed VoigtWithExpTail to VoigtPlusBortel.
 20260113: Major version incremented from 0 to 1 to account for VoigtPlusBortel skew model (breaking change).
           For backward compatibility, peaks without VoigtPlusBortel/GaussPlusBortel/DoubleBortel still write version 0.2.
 20231031: Minor version incremented from 1 to 2 to account for new peak skew models, and removal of Landau skew model.
 */

// For backward compatibility, peaks without VoigtPlusBortel/GaussPlusBortel/DoubleBortel write the pre-Voigt version
const int pre_Voigt_serialization_major = 0;
const int pre_Voigt_serialization_minor = 2;

const bool PeakDef::sm_defaultUseForDrfIntrinsicEffFit = true;
const bool PeakDef::sm_defaultUseForDrfFwhmFit = true;
const bool PeakDef::sm_defaultUseForDrfDepthOfInteractionFit = false;




void findROIEnergyLimits( double &lowerEnengy, double &upperEnergy,
                         const PeakDef &peak, 
                         const std::shared_ptr<const Measurement> &data,
                         const bool isHPGe )
{
  std::shared_ptr<const PeakContinuum> continuum = peak.continuum();
  if( continuum->energyRangeDefined() )
  {
    lowerEnengy  = continuum->lowerEnergy();
    upperEnergy = continuum->upperEnergy();
    return;
  }//if( continuum->energyRangeDefined() )
  
  if( !data || (data->num_gamma_channels() < 2) )
  {
    lowerEnengy = peak.lowerX();
    upperEnergy = peak.upperX();
    return;
  }//if( !data )
  
  const size_t lowbin = findROILimit( peak, data, false, isHPGe );
  const size_t upbin  = findROILimit( peak, data, true, isHPGe );
  if( lowbin == 0 )
    lowerEnengy = data->gamma_channel_center( lowbin );
  else
    lowerEnengy = data->gamma_channel_lower( lowbin );
  
  if( (upbin+1) >= data->num_gamma_channels() )
    upperEnergy = data->gamma_channel_center( std::min(upbin,data->num_gamma_channels()-1) );
  else
    upperEnergy = data->gamma_channel_upper( upbin );
}//void findROIEnergyLimits(...)


#define PRINT_ROI_DEBUG_INFO 0

#if( PRINT_ROI_DEBUG_INFO )
namespace
{
  class DebugLog
  {
    //Make it so cout/cerr statments always end up non-interleaved when multiple
    // threads are calling cout/cerr
  public:
    explicit DebugLog( std::ostream &os ) : os(os) {}
    ~DebugLog() { os << ss.rdbuf() << std::flush; }
    template <typename T>
    DebugLog& operator<<(T const &t){ ss << t; return *this;}
  private:
    std::ostream &os;
    std::stringstream ss;
  };
}//namespace
#endif



size_t findROILimitHighRes( const PeakDef &peak, const std::shared_ptr<const Measurement> &dataH, bool high )
{
  if( !dataH || !dataH->energy_calibration() || !dataH->energy_calibration()->valid() )
    return 0;
  
  const float mean = static_cast<float>( peak.mean() );
  const float fwhm = static_cast<float>( peak.fwhm() );
  
  // The plan is to iterate over the V&V test files and adjust these next quantities to most closely
  //  match the human selected ROI widths.
  //  Also, need to investigate highres_shrink_roi(...) to work better, for example the Cs137 peak of a shielded specturm
  const double min_width_mult = 0.75;
  const float nominal_fwhm_mult = 2.0f;
  const float steep_fwhm_mult = 3.5f;
  const float steep_continuum_limit = 0.35f;
  const float feature_nsigma_limit = 2.25f;
  
  
  const int direction = high ? 1 : -1;
  float nominal_low_edge = mean - (nominal_fwhm_mult * fwhm);
  size_t nominal_low_channel = dataH->find_gamma_channel( nominal_low_edge );
  
  float nominal_up_edge = mean + (nominal_fwhm_mult * fwhm);
  size_t nominal_up_channel = dataH->find_gamma_channel( nominal_up_edge );
  
  double nominal_data_area = dataH->gamma_channels_sum( nominal_low_channel, nominal_up_channel );
  
  double coefficients[2];
  const size_t nSideChanel = 3;
  PeakContinuum::eqn_from_offsets( nominal_low_channel, nominal_up_channel,
                                   mean, dataH, nSideChanel, nSideChanel, coefficients[1], coefficients[0] );
  
  // cout << "coefficients[0]=" << coefficients[0] << ", coefficients[1]=" << coefficients[1]
  //  << "\n\tdataH->gamma_channel_lower(nominal_low_channel)=" << dataH->gamma_channel_lower(nominal_low_channel)
  //  << "\n\tdataH->gamma_channel_upper(nominal_up_channel)=" << dataH->gamma_channel_upper(nominal_up_channel)
  //  << endl << endl;
  
  double nominal_cont_area = PeakContinuum::offset_eqn_integral( coefficients,
                                                 PeakContinuum::OffsetType::Linear,
                                                 dataH->gamma_channel_lower(nominal_low_channel),
                                                 dataH->gamma_channel_upper(nominal_up_channel),
                                                 mean );
  
  const double nominal_peak_area = (nominal_data_area - nominal_cont_area);
  const double nominal_uncert_ratio = nominal_peak_area / sqrt(nominal_data_area);
  
  //cout << "\n\n\nmean=" << mean << ", nominal_data_area=" << nominal_data_area << endl;
  //cout << "nominal_peak_area=" << nominal_peak_area << ", nominal_cont_area=" << nominal_cont_area
  //     << ", nominal_uncert_ratio=" << nominal_uncert_ratio << endl;
  
  const double width = dataH->gamma_channel_upper(nominal_up_channel) - dataH->gamma_channel_lower(nominal_low_channel);
  const double offset_contrib = coefficients[0]*width;
  const double slope_contrib = 0.5*coefficients[1]*width*width;
  
  //cout << "For " << mean << " kev, offset contributes " << offset_contrib << ", while slope contributes "
  //     << slope_contrib << "--->" << slope_contrib/offset_contrib << endl;
  
  if( (offset_contrib <= 0.0) || ((direction*slope_contrib/offset_contrib) > steep_continuum_limit) ) // 0.35 is arbitrary
  {
    if( high )
    {
      nominal_up_edge = mean + (steep_fwhm_mult * fwhm);
      nominal_up_channel = dataH->find_gamma_channel( nominal_up_edge );
    }else
    {
      nominal_low_edge = mean - (steep_fwhm_mult * fwhm);
      nominal_low_channel = dataH->find_gamma_channel( nominal_low_edge );
    }
  }//if( coefficients[1] < -30 )
  
  
  const size_t feature_detect_start = dataH->find_gamma_channel( mean + direction*min_width_mult*fwhm );
  const size_t feature_detect_stop = high ? nominal_up_channel : nominal_low_channel;

  
  // Feature detect start FWHM mult
  const size_t nchannel = dataH->num_gamma_channels();
  const size_t nprev_avrg = 4;
  for( size_t channel = feature_detect_start; channel != feature_detect_stop; channel += direction )
  {
    if( (channel <= nprev_avrg) || ((channel + nprev_avrg + 1) >= nchannel) )
      continue;

    const size_t prev_sum_start = channel - direction*nprev_avrg;
    const size_t prev_sum_end = channel - direction;

    const float prev_sum = dataH->gamma_channels_sum( prev_sum_start, prev_sum_end );
    const float val = dataH->gamma_channel_content( channel );
    const float val_next = dataH->gamma_channel_content( channel + direction );

    const float prev_avrg = prev_sum / nprev_avrg;
    const float prev_avrg_uncert = std::max( 1.0f*nprev_avrg, std::sqrt(prev_sum) ) / nprev_avrg;

    const float max_allowable = std::ceil(prev_avrg + feature_nsigma_limit*prev_avrg_uncert) +  0.001f;

    //cout << "mean=" << mean << ", channel=" << channel << ", prev_avrg(" << prev_sum_start
    //     <<  "," << prev_sum_end<< ")=" << prev_avrg
    //     << ", val=" << val << ", nextval=" << val_next << ", max_allowable=" << max_allowable << ""
    //     << endl;

    // We'll require two bins to be outside of tolerance
    if( (val > max_allowable && val_next > max_allowable) )
    {
      // When the peaks FWHM is smaller than the channel width (e.g., a synthetic delta-like
      //  peak), `feature_detect_start` can be the peaks own channel, in which case backing off
      //  by one channel here would put the ROI limit on the wrong side of the peaks mean -
      //  giving an inverted (lower above upper) ROI once both sides do this - so clamp the limit
      //  to the peaks mean channel.
      const size_t limit = channel - direction;
      const size_t mean_channel = dataH->find_gamma_channel( mean );
      return high ? std::max(limit, mean_channel) : std::min(limit, mean_channel);
    }
  }//for( int bin = minBin; bin > lastbin; --bin )


  return high ? nominal_up_channel : nominal_low_channel;
}//findROILimitHighRes(...)
 

size_t findROILimit( const PeakDef &peak, 
                    const std::shared_ptr<const Measurement> &dataH,
                    const bool high,
                    const bool isHPGe )
{
  if( !peak.gausPeak() )
    return dataH->find_gamma_channel( (high ? peak.upperX() : peak.lowerX()) );
  

  //This implemntation is an adaptation of how PCGAP defines a region of
  //  interest.
  //  The basic idea is to include a maximum of 11.75 sigma away from the mean
  //  of the peak, but then start at ~1.5 sigma from mean, and try to detect if
  //  a new feature is occuring, and if so, stop the region of interest there.
  //  A feature is "detected" if the value of the bin contents exceeds 2.5 sigma
  //  from the "expected" value, where the expected value starts off being the
  //  smallest bin value so far (well, this bin averaged with the bins on either
  //  side of it), and then each preceeding bin is added to this background
  //  value.
  //
  //References for PCGAP are at:
  //  http://www.inl.gov/technicalpublications/Documents/3318133.pdf
  //  http://www.osti.gov/bridge/servlets/purl/800710-Zc4iYJ/native/800710.pdf
  
  typedef int indexing_t;
  //  typedef size_t indexing_t;
  
  const indexing_t nchannel = (!dataH ? indexing_t(0): static_cast<indexing_t>(dataH->num_gamma_channels()));
  
  if( nchannel<128 )
    throw runtime_error( "findROILimit(...): Invalid input" );
  
  const vector<float> &contents = *dataH->gamma_channel_contents();
  
  
  std::shared_ptr<const PeakContinuum> continuum = peak.continuum();
  double lowxrange  = continuum->lowerEnergy();
  double highxrange = continuum->upperEnergy();
  
  const bool definedRange = (lowxrange != highxrange);
  
  if( definedRange )
    return dataH->find_gamma_channel( (high ? (highxrange-0.00001) : lowxrange) );
  
  if( isHPGe )
    return findROILimitHighRes( peak, dataH, high );
  
  
  const double mean = peak.mean();
  const double sigma = peak.sigma();

  const int direction = high ? 1 : -1;
  lowxrange = mean + direction * 7.5*sigma;  //2.3 FWHM
  
  const indexing_t nSideChannel = 1;
  
  indexing_t startchannel = static_cast<indexing_t>( dataH->find_gamma_channel( mean + direction*1.5*sigma ) );
  if( !high && startchannel < nSideChannel )
    startchannel = nSideChannel;
  
  indexing_t minChannel = startchannel;
  double minVal = contents[startchannel];
  indexing_t nbackbin = 1 + 2*nSideChannel;
  float backval = dataH->gamma_channels_sum( startchannel - nSideChannel,
                                             startchannel + nSideChannel );
  backval = std::max( backval, static_cast<float>(nbackbin) );
  
  const indexing_t meanchannel = static_cast<indexing_t>( dataH->find_gamma_channel( mean ) );
  indexing_t lastchannel  = static_cast<indexing_t>( dataH->find_gamma_channel( lowxrange + direction*0.0001 ) );
  
  //Make sure were not looking to far and that loop will terminate
  if( high ) lastchannel = std::min( lastchannel, nchannel-2 );
  else       lastchannel = std::max( lastchannel, indexing_t(1) );
  if( (high && (startchannel > lastchannel)) || (!high && (startchannel < lastchannel)) )
    startchannel = lastchannel;
  if( high || lastchannel>=direction )
    lastchannel += direction;
  
  
#if( PRINT_ROI_DEBUG_INFO )
  const bool debug_this_peak = (fabs(peak.mean() - 12.15) < 1.0); // && direction>0;
  const vector<float> &energies = *dataH->gamma_channel_energies();
  
  if( debug_this_peak )
    DebugLog(cerr) << "\n\n\n\nTo start with, lowxrange=" << lowxrange << " lastbin="
    << lastchannel << ", x(lastbin)=" << energies[lastchannel]
    << " for peak at mean=" << mean << ", sigma=" << sigma << "\n";
#endif
  
  //Lets find the bin with the smallest contents, and
  for( indexing_t channel = startchannel + direction*nSideChannel;
       channel != lastchannel && channel>nSideChannel && channel < nchannel; channel += direction )
  {
    assert( channel < (dataH->num_gamma_channels() + 100) );
    const float val = contents[channel];
    
    if( val <= minVal && (!isHPGe || (contents[channel+direction] <= minVal)) )
    {
      minVal = val;
      minChannel = channel;
      nbackbin = 1 + 2*nSideChannel;
      backval = dataH->gamma_channels_sum( channel-nSideChannel,
                                           channel+nSideChannel );
      backval = std::max( backval, static_cast<float>(nbackbin) );
      channel += direction;
      if( channel == lastchannel )
        break;
#if( PRINT_ROI_DEBUG_INFO )
      if( debug_this_peak )
        DebugLog(cerr) << "FindMin: New background val at " << energies[channel]
        << ", backval/nbackbin=" << backval/nbackbin << "\n";
#endif
    }else
    {
      //If val is greater than backval, we _may_ be hitting a new feature, like
      //  a new peak, if we are, lets stop going any further away from the mean.
      //  Not detecting this may cause us to extend _past_ the new feature and
      //  find a global minimum we clearly dount want
      
      float background = backval / nbackbin;
      float background_sigma = sqrt(backval) / nbackbin;
      
      //For low resolution spectra, lets estimate the slope of the continuum, and
      //  use this to correct maximum_allowable number of counts; without doing this
      //  the lower ROI range for peaks on a falling continuum can be much to short.
      //  (not tested for HPGe)

      if( (nbackbin > 2) && !isHPGe && channel > nbackbin )
      {
        try
        {
          const float *x = &(*dataH->channel_energies())[0];
          const float *y = &(*dataH->gamma_counts())[0];
          x = x + channel - direction - ((direction < 0) ? 0 : nbackbin);
          y = y + channel - direction - ((direction < 0) ? 0 : nbackbin);
          
          vector<double> coeffs, uncerts;
          fit_to_polynomial( x, y, nbackbin, 1, coeffs, uncerts );
          
          const float thisx = ((direction < 0) ? x[direction] : x[nbackbin+1]);
          background = std::max( coeffs[0] + thisx*coeffs[1], 1.0*nbackbin );
        }catch(...)
        {
        }
      }//if( nbackbin > 2 )

      const float sigma = sqrt( background_sigma*background_sigma + background );
      //      double max_allowable = std::ceil(background + 2.575829*sigma) + 0.001;
      float max_allowable = std::ceil(background + 3.0f*sigma) + 0.001f;
      
      if( background < 20.0f )
      {
        max_allowable = static_cast<float>( LightMath::poisson_quantile( background, 0.99f ) );
      }//if( background < 20 )
      
      
      
      if( val > max_allowable && (!isHPGe || contents[channel+direction] > max_allowable) )
      {
        //XXX - the below 3 is purely empircal, and meant to help avoid
        //      contamination due to the new feature
        if( channel >= 3*direction )
          lastchannel = channel - 3*direction;
        
#if( PRINT_ROI_DEBUG_INFO )
        if( debug_this_peak )
          DebugLog(cerr) << "FindMin: Found last bin at " << dataH->gamma_channel_center(lastchannel)
          << ", val=" << val << ", max_allowable=" << max_allowable << "\n";
#endif
        break;
      }else
      {
        ++nbackbin;
        backval += val;
#if( PRINT_ROI_DEBUG_INFO )
        if( debug_this_peak )
          DebugLog(cerr) << "FindMin: bin " << channel << ", x(bin)=" << energies[channel]
          << ", val=" << val << ", max_allowable=" << max_allowable
          << ", background=" << background << ", sigma=" << sigma << "\n";
#endif
      }
    }//if( val <= minVal ) / else
  }//for( ; bin > lastbin; --bin )
  
  nbackbin = 1 + 2*nSideChannel;
  backval = dataH->gamma_channels_sum( minChannel - nSideChannel,
                                       minChannel + nSideChannel );
  backval = std::max( backval, static_cast<float>(nbackbin) );
  
#if( PRINT_ROI_DEBUG_INFO )
  if( debug_this_peak )
    DebugLog(cerr) << "1) Bin with smallest contends at "
    << dataH->gamma_channel_center(minChannel)
    << ", backval=" << (backval/3.0) << ", lastbin=" << lastchannel
    << " x(lastbin)=" << dataH->gamma_channel_center(lastchannel) << "\n";
#endif
  
  //Make sure were not looking to far and that loop will terminate
  if( high ) lastchannel = std::min( lastchannel, nchannel-2 );
  else       lastchannel = std::max( lastchannel, indexing_t(1) );
  if( (high && (minChannel > lastchannel)) || (!high && (minChannel < lastchannel)) )
    minChannel = lastchannel;
  
  if( lastchannel || direction>0 )
    lastchannel += direction;
  
#if( PRINT_ROI_DEBUG_INFO )
  if( debug_this_peak )
    DebugLog(cerr) << "2) Bin with smallest contends at " << energies[minChannel]
    << ", backval=" << (backval/3.0) << ", lastbin=" << lastchannel
    << " x(lastbin)=" << energies[lastchannel] << "\n";
#endif
  
  
  if( direction < 0 && ((float(lastchannel)/float(nchannel)) < 0.04) )
  {
    size_t lower_channel = 0, upper_channel = 0;
    PeakFitUtils::find_spectroscopic_extent( dataH, lower_channel, upper_channel );
    if( static_cast<int>(lower_channel) >= lastchannel )
    {
      lastchannel = lower_channel ? static_cast<indexing_t>(lower_channel - 1) : 0;
      minChannel = std::max( minChannel, lastchannel );
    }
  }//if( direction < 0 && dataH->GetBinCenter(lastbin) < 100.0 )
  
  
  // Bound `channel` to the valid channel range.  The termination test only checks for reaching
  // lastchannel/meanchannel, so on a degenerate low-resolution spectrum where the scan direction
  // points away from both, `channel` would run to INT_MAX and `contents[channel]` below would read
  // far out of bounds (heap OOB -> SIGBUS).  For valid inputs the loop stops at lastchannel/
  // meanchannel first, so this bound only affects the runaway case.
  for( indexing_t channel = minChannel + direction*nSideChannel;
      channel != lastchannel && channel != meanchannel
        && channel >= 0 && channel < nchannel; channel += direction )
  {
    const float val = contents[channel];
    const float nextval = (channel>1 && (nchannel-channel)>0)  //probably is fine, but we'll check JIC
                          ? contents[channel+direction]
                          : contents[channel];
    
    float back = backval / nbackbin;
    float back_uncert = std::sqrt( backval ) / nbackbin;
    const float sigma = std::sqrt( back + back_uncert*back_uncert );
    //    float max_allowable = std::ceil(back + 2.575829*sigma) +  0.001;
    //    float min_allowable = std::floor(back - 2.575829*sigma) - 0.001;
    float max_allowable = std::ceil(back + 2.8f*sigma) +  0.001f;
    float min_allowable = std::floor(back - 2.8f*sigma) - 0.001f;
    
    if( back < 20.0f )
    {
      max_allowable = static_cast<float>( LightMath::poisson_quantile( back, 0.99f ) );
      min_allowable = static_cast<float>( LightMath::poisson_quantile( back, 0.01f ) );
    }//if( val < 20.0 )
    
    //For high resolution spectra we'll require two bins to be outside of
    //  tollerance, since this will preserve the intent, but allow single bin
    //  spikes (which I swear are more common than poisson!) to not mess up the
    //  ROI.  It could probably be done for low resolution spectra, but I havent
    //  tested if (since single lowres peaks ROIs typically get calculated by
    //  find_roi_for_2nd_deriv_candidate(...) anyway
    if( (val>max_allowable && (!isHPGe ||nextval>max_allowable))
        || (val<min_allowable && (!isHPGe || nextval<min_allowable)) )
    {
#if( PRINT_ROI_DEBUG_INFO )
      if( debug_this_peak )
        DebugLog(cerr) << "Setting bin to " << (channel - direction) << " from expected "
        << lastchannel << " at energy=" << energies[channel-direction] << " kev"
        << ", this bring limit to " << (mean-energies[channel-direction])/sigma
        << " sigma from mean, val=" << val << ", min_allowable="
        << min_allowable << ", max_allowable=" << max_allowable << "\n";
#endif
      lastchannel = channel - (channel>0 ? direction : 0);
      break;
    }//if( val>max_allowable || val<min_allowable )
#if( PRINT_ROI_DEBUG_INFO )
    else
    {
      if( debug_this_peak )
        DebugLog(cerr) << "FindLimit: bin " << channel << ", x(bin)=" << energies[channel]
        << ", val=" << val
        << ", min_allowable=" << min_allowable
        << ", max_allowable=" << max_allowable
        << ", back=" << back << ", sigma=" << sigma << "\n";
    }
#endif
    
    nbackbin++;
    backval += val;
  }//for( int bin = minBin; bin > lastbin; --bin )
  
  
  //In principle, lastbin is the furthest from the mean we can end up, with
  //  the maximum being 11.75*sigma, or whereever a new feature was detected
  
#if( PRINT_ROI_DEBUG_INFO )
  if( debug_this_peak )
    DebugLog(cerr) << "minBin was " << minChannel << ", lastbin=" << lastchannel
    << " x(lastbin)=" << energies[lastchannel] << "\n";
#endif
  
  
  //Try to detect if there is a signficant skew on the peak by comparing
  //  4 to 7 sigma, to 7 to 11.75 sigma (or wherever is lastbin) to see if they
  //  are statistically compatible; if they are, just have ROI go to 7 sigma
  const int mean_channel = static_cast<int>( dataH->find_gamma_channel( mean ) );
  const int good_cont_channel = static_cast<int>( dataH->find_gamma_channel( mean + direction*7.05*sigma ) );
  if( (std::abs(int(lastchannel)-mean_channel) > std::abs(good_cont_channel-mean_channel)) )
  {
    const indexing_t nearest_channel = static_cast<indexing_t>( dataH->find_gamma_channel( mean + direction*3.5*sigma ) );
    const bool isNotDecreasing = isStatisticallyGreaterOrEqual( nearest_channel, good_cont_channel,
                                            good_cont_channel, lastchannel, dataH, 3.0 );

    
    
    if( high || isNotDecreasing )
    {
#if( PRINT_ROI_DEBUG_INFO )
      if( debug_this_peak )
        DebugLog(cerr) << "Setting lastchannel to " << lastchannel << "\n";
#endif
      lastchannel = good_cont_channel;
    }
    //    else
    //    {
    //      //now check from 11.75 to to 16.5 sigma, to see if we should include down
    //      //  to there.
    //      const int start = good_cont_bin;
    //      const size_t end_channel = dataH->find_gamma_channel( mean + direction*16.5*sigma );
    //      isNotDecreasing = isStatisticallyGreaterOrEqual( good_cont_bin, lastbin,
    //                                                      start, end_channel, dataH, 2.0 );
    //      if( !isNotDecreasing )
    //        lastbin = end;
    //    }
  }//if( abs(lastbin-mean_bin) > abs(good_cont_bin-mean_bin) )
  
  if( lastchannel < 0 )
    lastchannel = 0;
  else if( lastchannel >= static_cast<indexing_t>(dataH->num_gamma_channels()) )
    lastchannel = static_cast<indexing_t>(dataH->num_gamma_channels()) - 1;
  
  float val = dataH->gamma_channel_center(lastchannel);
  if( direction < 0 )
  {
    if( ((mean-val)/sigma) < 1.75 )
      lastchannel = static_cast<indexing_t>( dataH->find_gamma_channel( mean - 1.75*sigma ) );
  }else
  {
    if( ((val-mean)/sigma) < 1.75 )
      lastchannel = static_cast<indexing_t>( dataH->find_gamma_channel( mean + 1.75*sigma ) );
  }
  
#if( PRINT_ROI_DEBUG_INFO )
  if( debug_this_peak )
    DebugLog(cerr) << "Returning bin " << lastchannel
    << ", x=" << dataH->gamma_channel_center(lastchannel) << "\n";
#endif
  
  return lastchannel;
}//int findROILimit(...)




bool isStatisticallyGreaterOrEqual( const size_t start1, const size_t end1,
                                    const size_t start2, const size_t end2,
                                    const std::shared_ptr<const Measurement> &dataH,
                                    const double nsigma )
{
  size_t lowerbin = std::min( start1, end1 );
  size_t upperbin = std::max( start1, end1 );
  const double lower_area = dataH->gamma_channels_sum( lowerbin, upperbin );
  const int num_near_mean_bins = static_cast<int>( upperbin - lowerbin + 1 );
  const double avrg_near_mean_area = lower_area / num_near_mean_bins;
  const double avrg_near_mean_uncert = sqrt(lower_area) / num_near_mean_bins;
  
  lowerbin = std::min( start2, end2 );
  upperbin = std::max( start2, end2 );
  const double upper_area = dataH->gamma_channels_sum( lowerbin, upperbin );
  const int num_tail_bins = static_cast<int>( upperbin - lowerbin + 1 );
  const double tail_area = upper_area / num_tail_bins;
  const double avrg_tail_uncert = sqrt(upper_area) / num_tail_bins;
  
  const double uncert = sqrt(avrg_near_mean_uncert*avrg_near_mean_uncert
                             + avrg_tail_uncert*avrg_tail_uncert);
  
/*
#if( PRINT_ROI_DEBUG_INFO )
    DebugLog(cerr) << "Found avrg_near_mean_area=" << avrg_near_mean_area << "+-" << avrg_near_mean_uncert
    << ", tail_area=" << tail_area << "+-" << avrg_near_mean_uncert
    << ", uncert=" << uncert
    << " ---> (avrg_near_mean_area-2.0*uncert)=" << (avrg_near_mean_area-nsigma*uncert)
    << ", nsigma=" << ((avrg_near_mean_area-tail_area)/uncert)
    << "\n";
#endif
*/
  
  return (tail_area > (avrg_near_mean_area-nsigma*uncert));
}//isStatisticallyGreaterOrEqual( ... )


void estimatePeakFitRange( const PeakDef &peak, const std::shared_ptr<const Measurement> &dataH,
                           const bool isHPGe,
                           size_t &lower_channel, size_t &upper_channel )
{
  const size_t nchannel = dataH ? dataH->num_gamma_channels() : size_t(0);
  if( !nchannel )
    return;
  
  std::shared_ptr<const PeakContinuum> continuum = peak.continuum();
  double lowxrange  = continuum->lowerEnergy();
  double highxrange = continuum->upperEnergy();
  
  const bool definedRange = (lowxrange != highxrange);
  if( definedRange )
  {
    lower_channel  = max( dataH->find_gamma_channel(lowxrange), size_t(0) );
    upper_channel = min( dataH->find_gamma_channel(highxrange-0.00001), nchannel-1 );
    return;
  }//if( definedRange )
  
  const double mean = peak.mean();
  const double sigma = peak.gausPeak() ? peak.sigma() : 0.5*0.25*peak.roiWidth();
  
  
  if( continuum->type() == PeakContinuum::External )
  {
    lowxrange  = mean - 4.0*sigma;
    highxrange = mean + 4.0*sigma;
    
    lower_channel = dataH->find_gamma_channel(lowxrange);
    upper_channel = dataH->find_gamma_channel(highxrange);
    return;
  }//if( peak.m_offsetType == PeakDef::External )
  
  
  const bool polyContinuum = continuum->isPolynomial();
  if( polyContinuum )
  {
    lower_channel = findROILimit( peak, dataH, false, isHPGe );
    upper_channel = findROILimit( peak, dataH, true, isHPGe );
  }else
  {
    lower_channel = dataH->find_gamma_channel( mean - 4.0*sigma );
    upper_channel = dataH->find_gamma_channel( mean + 4.0*sigma );
  }
  
  if( lower_channel > upper_channel )
      std::swap( lower_channel, upper_channel );
      
  //Lets avoid some wierd going to too small of peak widths
  const size_t numfitbin = upper_channel - lower_channel;
  if( numfitbin <= 9 )  //9 chosen arbitrarily
  {
    lower_channel -= (10-numfitbin)/2;
    upper_channel += (10-numfitbin)/2;
  }//if( numfitbin <= 9 )
}//void setPeakXLimitsFromData( PeakDef &peak, const std::shared_ptr<const Measurement> &dataH )


ostream &operator<<( std::ostream &stream, const PeakContinuum &cont )
{
  switch( cont.type() )
  {
    case PeakContinuum::NoOffset:
      stream << "Underfined continuum";
    break;
      
    case PeakContinuum::External:
      stream << "Globally defined continuum";
    break;
      
    case PeakContinuum::FlatStep:
    case PeakContinuum::LinearStep:
    case PeakContinuum::BiLinearStep:
    {
      const char * const names[] = {"Flat", "Linear", "Bi-linear"};
      stream << names[cont.type() - PeakContinuum::FlatStep] << " step with coefficients {";
      for( size_t i = 0; i < cont.m_values.size(); ++i )
        stream << (i?", ":"") << cont.m_values[i];
      stream << "} relative to " << cont.m_referenceEnergy << " keV";
      break;
    }

    case PeakContinuum::FlatStepCDF:
    case PeakContinuum::LinearStepCDF:
    {
      const char * const names[] = {"Flat", "Linear"};
      stream << names[cont.type() - PeakContinuum::FlatStepCDF] << " step (CDF) with coefficients {";
      for( size_t i = 0; i < cont.m_values.size(); ++i )
        stream << (i?", ":"") << cont.m_values[i];
      stream << "} relative to " << cont.m_referenceEnergy << " keV";
      break;
    }

    case PeakContinuum::BiLinearStepCDF:
    {
      stream << "Bi-linear step (CDF) with coefficients {";
      for( size_t i = 0; i < cont.m_values.size(); ++i )
        stream << (i?", ":"") << cont.m_values[i];
      stream << "} relative to " << cont.m_referenceEnergy << " keV";
      break;
    }

    case PeakContinuum::Constant:   case PeakContinuum::Linear:
    case PeakContinuum::Quadratic: case PeakContinuum::Cubic:
      stream << "Polynomial continuum with values {";
      for( size_t i = 0; i < cont.m_values.size(); ++i )
        stream << (i?", ":"") << cont.m_values[i];
      stream << "} relative to " << cont.m_referenceEnergy << " keV";
    break;
  }//switch( m_type )

  stream << ", valid from " << cont.lowerEnergy()
         << " to " << cont.upperEnergy() << " keV";
  
  return stream;
}//operator<<( std::ostream &stream, const PeakContinuum &cont )


std::ostream &operator<<( std::ostream &stream, const PeakDef &peak )
{
  stream << "mean=" << peak.m_coefficients[PeakDef::Mean];
  if( peak.m_uncertainties[PeakDef::Mean] > 0.0 )
    stream << "+-" << peak.m_uncertainties[PeakDef::Mean];
  stream << ", sigma=" << peak.m_coefficients[PeakDef::Sigma];
  if( peak.m_uncertainties[PeakDef::Sigma] > 0.0 )
    stream << "+-" << peak.m_uncertainties[PeakDef::Sigma];
  stream << ", amplitude=" << peak.m_coefficients[PeakDef::GaussAmplitude];
  if( peak.m_uncertainties[PeakDef::GaussAmplitude] > 0.0 )
    stream << "+-" << peak.m_uncertainties[PeakDef::GaussAmplitude];

  if( peak.m_source.kind != PeakDef::Source::Kind::None )
    stream << ", source=" << peak.m_source.name << " " << peak.m_source.particle_energy << " keV";
  
  stream << ", " << *peak.m_continuum
         << ", chi2=" << peak.m_coefficients[PeakDef::Chi2DOF]
         << ", Skew0=" << peak.m_coefficients[PeakDef::SkewPar0]
         << ", Skew1=" << peak.m_coefficients[PeakDef::SkewPar1]
         << ", Skew2=" << peak.m_coefficients[PeakDef::SkewPar2]
         << ", Skew3=" << peak.m_coefficients[PeakDef::SkewPar3]
         << ", Skew4=" << peak.m_coefficients[PeakDef::SkewPar4]
         << ", Skew5=" << peak.m_coefficients[PeakDef::SkewPar5]
         << std::flush;
  return stream;
}//std::ostream &operator<<( std::ostream &stream, const PeakDef &peak )





PeakDef::PeakDef()
{
  reset();
}


void PeakDef::reset()
{
  m_userLabel                 = "";
  m_type                      = GaussianDefined;
  m_skewType                  = PeakDef::NoSkew;

  m_source                    = Source();
  m_useForEnergyCal           = true;
  m_useForShieldingSourceFit  = false;
  m_useForManualRelEff        = true;
  
  m_useForDrfIntrinsicEffFit       = PeakDef::sm_defaultUseForDrfIntrinsicEffFit;
  m_useForDrfFwhmFit               = PeakDef::sm_defaultUseForDrfFwhmFit;
  m_useForDrfDepthOfInteractionFit = PeakDef::sm_defaultUseForDrfDepthOfInteractionFit;
  
  m_lineColor                 = Wt::WColor();
  
  std::shared_ptr<PeakContinuum> newcont = std::make_shared<PeakContinuum>();
  m_continuum = newcont;
  
  for( CoefficientType t = CoefficientType(0);
       t < NumCoefficientTypes; t = CoefficientType(t+1) )
  {
    m_coefficients[t] = 0.0;
    m_uncertainties[t] = -1.0;
    
    switch( t )
    {
      case PeakDef::Mean:
      case PeakDef::Sigma:
      case PeakDef::GaussAmplitude:
      case PeakDef::SkewPar0:
      case PeakDef::SkewPar1:
      case PeakDef::SkewPar2:
      case PeakDef::SkewPar3:
      case PeakDef::SkewPar4:
      case PeakDef::SkewPar5:
        m_fitFor[t] = true;
      break;

      case PeakDef::Chi2DOF:
      case PeakDef::NumCoefficientTypes:
        m_fitFor[t] = false;
      break;
    }//switch( type )
  }//for( loop over coefficients )
}//void PeakDef::reset()


PeakDef::PeakDef( double m, double s, double a )
{
  reset();
  m_coefficients[PeakDef::Mean] = m;
  m_coefficients[PeakDef::Sigma] = s;
  m_coefficients[PeakDef::GaussAmplitude] = a;
}



PeakDef::PeakDef( double xlow, double xhigh, double mean,
                    std::shared_ptr<const Measurement> data, std::shared_ptr<const Measurement> background )
{
  reset();
  m_type = PeakDef::DataDefined;
  m_coefficients[PeakDef::Mean] = mean;
  m_continuum->setRange( xlow, xhigh );
  
  if( !data )
    return;

  m_continuum->setType( PeakContinuum::External );
  m_continuum->setExternalContinuum( background );
  m_continuum->setRange( xlow, xhigh );
  
  m_coefficients[PeakDef::GaussAmplitude] = gamma_integral( data, xlow, xhigh );
  if( background )
    m_coefficients[PeakDef::GaussAmplitude] -= gamma_integral( background, xlow, xhigh );
}//PeakDef( constructor )



const char *PeakDef::to_string( const CoefficientType type )
{
  switch( type )
  {
    case PeakDef::Mean:                return "Centroid";
    case PeakDef::Sigma:               return "Width";
    case PeakDef::GaussAmplitude:      return "Amplitude";
    case PeakDef::SkewPar0:            return "Skew0";
    case PeakDef::SkewPar1:            return "Skew1";
    case PeakDef::SkewPar2:            return "Skew2";
    case PeakDef::SkewPar3:            return "Skew3";
    case PeakDef::SkewPar4:            return "Skew4";
    case PeakDef::SkewPar5:            return "Skew5";
    case PeakDef::Chi2DOF:             return "Chi2";
    case PeakDef::NumCoefficientTypes: return "";
  }//switch( type )

  return "";
}//const char *PeakDef::to_string( const CoefficientType type )


const char *PeakDef::to_string( const SkewType type )
{
  switch( type )
  {
    case PeakDef::NumSkewType:
      assert( 0 );
      //throw runtime_error( "PeakDef::to_string(SkewType): NumSkewType is not a valid skew type" );
      // Fall through to NoSkew for non-debug builds

    case PeakDef::NoSkew:                 return "NoSkew";
    case PeakDef::Bortel:                 return "ExGauss";
    case PeakDef::DoubleBortel:           return "DoubleExGauss";
    case PeakDef::GaussPlusBortel:        return "GaussPlusExGauss";
    case PeakDef::GaussExp:               return "GaussExp";
    case PeakDef::CrystalBall:            return "CrystalBall";
    case PeakDef::ExpGaussExp:            return "ExpGaussExp";
    case PeakDef::DoubleSidedCrystalBall: return "DoubleSidedCrystalBall";
    case PeakDef::VoigtPlusBortel:        return "VoigtPlusExGauss";
    case PeakDef::GadrasGeneric:          return "GadrasGeneric";
    case PeakDef::GadrasCZT:              return "GadrasCZT";
  }//switch( skew_type )

  assert( 0 );
  throw runtime_error( "PeakDef::to_string(SkewType): invalid SkewType" );
  return "";
}//const char *to_string( const SkewType type )


const char *PeakDef::to_label( const SkewType type )
{
  switch( type )
  {
    case PeakDef::NumSkewType:
      assert( 0 );
      //throw runtime_error( "PeakDef::to_label(SkewType): NumSkewType is not a valid skew type" );
      // Fall through to NoSkew for non-debug builds

    case PeakDef::NoSkew:                 return "None";
    case PeakDef::Bortel:                 return "Exp*Gauss";
    case PeakDef::DoubleBortel:           return "Double Exp*Gauss";
    case PeakDef::GaussPlusBortel:        return "Gauss+Exp*Gauss";
    case PeakDef::GaussExp:               return "GaussExp";
    case PeakDef::CrystalBall:            return "Crystal Ball";
    case PeakDef::ExpGaussExp:            return "ExpGaussExp";
    case PeakDef::DoubleSidedCrystalBall: return "Double Crystal Ball";
    case PeakDef::VoigtPlusBortel:        return "Voigt+Exp*Gauss";
    case PeakDef::GadrasGeneric:          return "GADRAS (generic)";
    case PeakDef::GadrasCZT:              return "GADRAS (CZT/CdTe)";
  }//switch( skew_type )

  assert( 0 );
  throw runtime_error( "PeakDef::to_label(SkewType): invalid SkewType" );
  return "";
}//const char *to_label( const SkewType type );


PeakDef::SkewType PeakDef::skew_from_string( const string &skew_type_str )
{
  if( SpecUtils::iequals_ascii(skew_type_str,"NoSkew")
          || SpecUtils::iequals_ascii(skew_type_str,"None") )
    return PeakDef::SkewType::NoSkew;

  // Exp*Gauss (Bortel) - accept multiple aliases, including display label "Exp*Gauss"
  if( SpecUtils::iequals_ascii(skew_type_str,"ExGauss")
          || SpecUtils::iequals_ascii(skew_type_str,"Bortel")
          || SpecUtils::iequals_ascii(skew_type_str,"EMG")
          || SpecUtils::iequals_ascii(skew_type_str,"ExpGauss")
          || SpecUtils::iequals_ascii(skew_type_str,"Exp*Gauss") )
    return PeakDef::SkewType::Bortel;

  // Double Exp*Gauss (DoubleBortel) - accept aliases, including display label
  if( SpecUtils::iequals_ascii(skew_type_str,"DoubleExGauss")
          || SpecUtils::iequals_ascii(skew_type_str,"DoubleBortel")
          || SpecUtils::iequals_ascii(skew_type_str,"Double Exp*Gauss") )
    return PeakDef::SkewType::DoubleBortel;

  // Gauss+Exp*Gauss (GaussPlusBortel) - accept aliases, including display label
  if( SpecUtils::iequals_ascii(skew_type_str,"GaussPlusExGauss")
          || SpecUtils::iequals_ascii(skew_type_str,"GaussPlusBortel")
          || SpecUtils::iequals_ascii(skew_type_str,"GaussBortel")
          || SpecUtils::iequals_ascii(skew_type_str,"Gauss+Exp*Gauss") )
    return PeakDef::SkewType::GaussPlusBortel;

  if( SpecUtils::iequals_ascii(skew_type_str,"GaussExp") )
    return PeakDef::SkewType::GaussExp;

  if( SpecUtils::iequals_ascii(skew_type_str,"CrystalBall")
          || SpecUtils::iequals_ascii(skew_type_str,"CB")
          || SpecUtils::iequals_ascii(skew_type_str,"Crystal Ball") )
    return PeakDef::SkewType::CrystalBall;

  if( SpecUtils::iequals_ascii(skew_type_str,"ExpGaussExp") )
    return PeakDef::SkewType::ExpGaussExp;

  if( SpecUtils::iequals_ascii(skew_type_str,"DoubleSidedCrystalBall")
          || SpecUtils::iequals_ascii(skew_type_str,"DSCB")
          || SpecUtils::iequals_ascii(skew_type_str,"Double Crystal Ball") )
    return PeakDef::SkewType::DoubleSidedCrystalBall;

  // Voigt+Exp*Gauss (VoigtPlusBortel) - accept old names and display label
  if( SpecUtils::iequals_ascii(skew_type_str,"VoigtPlusExGauss")
          || SpecUtils::iequals_ascii(skew_type_str,"VoigtPlusBortel")
          || SpecUtils::iequals_ascii(skew_type_str,"VoigtWithExpTail")
          || SpecUtils::iequals_ascii(skew_type_str,"VoigtExp")
          || SpecUtils::iequals_ascii(skew_type_str,"VoigtBortel")
          || SpecUtils::iequals_ascii(skew_type_str,"Voigt+Exp*Gauss") )
    return PeakDef::SkewType::VoigtPlusBortel;

  // GADRAS peak shapes
  if( SpecUtils::iequals_ascii(skew_type_str,"GadrasGeneric")
          || SpecUtils::iequals_ascii(skew_type_str,"GADRAS (generic)")
          || SpecUtils::iequals_ascii(skew_type_str,"GadrasGen") )
    return PeakDef::SkewType::GadrasGeneric;

  if( SpecUtils::iequals_ascii(skew_type_str,"GadrasCZT")
          || SpecUtils::iequals_ascii(skew_type_str,"GADRAS (CZT/CdTe)")
          || SpecUtils::iequals_ascii(skew_type_str,"GadrasCdTe") )
    return PeakDef::SkewType::GadrasCZT;


  throw runtime_error( "Invalid peak skew type: " + string(skew_type_str) );

  return PeakDef::SkewType::NoSkew;
}//SkewType skew_from_string( const char *skew_str )


bool PeakDef::skew_parameter_range( const SkewType skew_type, const CoefficientType coef,
                                   double &lower_value, double &upper_value,
                                   double &starting_value, double &step_size )
{
  lower_value = upper_value = starting_value = step_size = 0.0;
  
  switch( skew_type )
  {
    case NumSkewType:
      assert( 0 );
      //throw runtime_error( "PeakDef::skew_parameter_range: NumSkewType is not a valid skew type" );
      // Fall through to NoSkew for non-debug builds

    case NoSkew:
      return false;
      
    case SkewType::Bortel:
    {
      if( coef != CoefficientType::SkewPar0 )
        return false;

      // The smaller the skew paramater, the less skew there is
      starting_value = 0.5;
      step_size = 1.0;
      lower_value = 0.0; //Below 0.005 would be numerically bad, but the Bortel function should protect against it.
      upper_value = 15;
      
      break;
    }//case SkewType::Bortel:
      
    case SkewType::CrystalBall:
    case SkewType::DoubleSidedCrystalBall:
    {
      switch( coef )
      {
        case CoefficientType::SkewPar2: //alpha (right)
          if( skew_type == SkewType::CrystalBall )
            return false;
          // fall-though intentional
        case CoefficientType::SkewPar0: //alpha (left)
          starting_value = 2; // Saying skew becomes significant after 2 sigma, is maybe reasonable
          step_size = 0.5;
          lower_value = 0.5;  // You should at least be gaussian for half a sigma
          // Pure Gaussian only as alpha->inf; the power-law tail is <1e-5 of the area by ~5 sigma,
          //  so 5.0 is the practical "no skew" ceiling (a larger alpha also shrinks the (n/alpha)^n
          //  term, so widening this does NOT worsen the CrystalBall overflow corner, which is at
          //  small alpha).  The auto-simplify de-pin fixes alpha here to drop the skew DOF when the
          //  data prefers no tail.
          upper_value = 5.0;
          break;
        
          
        case CoefficientType::SkewPar3: //n (right)
          if( skew_type == SkewType::CrystalBall )
            return false;
          // fall-though intentional
        case CoefficientType::SkewPar1: //n (left)
          // `n` is the power-law tail exponent.  These bounds are the single source of truth shared
          //  by every peak-fit path and by RelActCalcAuto (both read this function), so the
          //  CrystalBall range is unified across them.
          //   - lower 1.05 keeps the 1/(n-1) tail-normalization pole bounded (|1/(n-1)| <= 20);
          //     1.0 is a hard divide-by-zero.
          //   - upper 100 is a deliberately generous ceiling.  It used to also cap the (n/alpha)^n
          //     constant that overflowed for large n / small alpha, but the CrystalBall tail is now
          //     evaluated in a cancellation-free form (PeakDists::crystal_ball_tail_indefinite_t)
          //     that never forms that constant, so 100 is numerically safe.  (The likelihood in `n`
          //     flattens for large n - the tail approaches a fixed exponential - so a tail-heavy peak
          //     may pin here; that is handled by the bound-aware uncertainty path, not this bound.)
          starting_value = 2;
          step_size = 0.75;
          lower_value = 1.05;
          upper_value = 100;
          break;
          
        default:
          return false;
      }//switch( coef )
      
      break;
    }//case SkewType::CrystalBall, SkewType::DoubleSidedCrystalBall:
      
    case SkewType::GaussExp:
    case SkewType::ExpGaussExp:
    {
      switch( coef )
      {
        case CoefficientType::SkewPar1:
          if( skew_type == SkewType::GaussExp )
            return false;
          // fall-though intentional
        case CoefficientType::SkewPar0:
          starting_value = 1;  //A pretty good amount of skew
          step_size = 0.2;
          lower_value = 0.15;  //this is a huge amount of skew (~82% of the area is in the tail)
          // GaussExp/ExpGaussExp reach a pure Gaussian only as skew->inf, so this upper bound is a
          //  finite proxy for "no skew".  The exponential tail is ~0.06% of the peak area at 3.25,
          //  but ~0.003% at 4.0 (tail = (sigma/s)*exp(-s^2/2) over the total `gauss_exp_norm`), so
          //  4.0 is the practical "no skew" ceiling.  When the data wants even less tail, the
          //  RelActAuto auto-simplify de-pin fixes skew here and drops the DOF (see
          //  PeakDef::skew_no_skew_value) rather than the fit slamming into - and pinning at - the bound.
          upper_value = 4.0;
          break;
          
        default:
          return false;
      }//switch( coef )

      break;
    }//case SkewType::ExpGaussExp: case SkewType::GaussExp:

    case SkewType::VoigtPlusBortel:
    {
      switch( coef )
      {
        case CoefficientType::SkewPar0: // gamma_lor - Lorentzian HWHM in keV
          starting_value = 0.1; // 5 eV typical for flourescent x-rays, but decay x-rays can be like 80 eV
          step_size = 0.02;      // 2 eV steps
          lower_value = 0.0;    // 1 eV minimum
          upper_value = 10.0;      // Not sure if this is valid...
          break;

        case CoefficientType::SkewPar1: // R - mixing ratio (0=Voigt, 1=Bortel)
          starting_value = 0.1;   // 10% Bortel is typical
          step_size = 0.05;       // 5% steps
          lower_value = 0.0;      // Pure Voigt
          upper_value = 1.0;      // Pure Bortel
          break;

        case CoefficientType::SkewPar2: // tau - Bortel skew parameter
          starting_value = 1.0;   // Similar to Bortel skew
          step_size = 0.2;
          lower_value = 0.01;     // Gentle tail
          upper_value = 10;       // Steep tail
          break;

        default:
          return false;
      }//switch( coef )

      break;
    }//case SkewType::VoigtPlusBortel:

    case SkewType::GaussPlusBortel:
    {
      switch( coef )
      {
        case CoefficientType::SkewPar0: // R - mixing ratio (0=Gaussian, 1=Bortel)
          starting_value = 0.5;   // 50% mix
          step_size = 0.1;
          lower_value = 0.0;      // Pure Gaussian
          upper_value = 1.0;      // Pure Bortel
          break;

        case CoefficientType::SkewPar1: // tau - Bortel skew parameter
          starting_value = 0.5;
          step_size = 1.0;
          lower_value = 0.01;
          upper_value = 15.0;
          break;

        default:
          return false;
      }//switch( coef )

      break;
    }//case SkewType::GaussPlusBortel:

    case SkewType::DoubleBortel:
    {
      switch( coef )
      {
        case CoefficientType::SkewPar0: // tau1 - first exponential decay constant
          starting_value = 0.5;
          step_size = 1.0;
          lower_value = 0.01;
          upper_value = 15.0;
          break;

        case CoefficientType::SkewPar1: // tau2_delta - positive delta (tau2 = tau1 + tau2_delta)
          starting_value = 1.0;
          step_size = 0.5;
          lower_value = 0.0;      // When 0, reduces to single Bortel
          upper_value = 10.0;
          break;

        case CoefficientType::SkewPar2: // eta - weight of second exponential
          starting_value = 0.5;
          step_size = 0.1;
          lower_value = 0.0;      // Only first exponential (tau1)
          upper_value = 1.0;      // Only second exponential (tau2)
          break;

        default:
          return false;
      }//switch( coef )

      break;
    }//case SkewType::DoubleBortel:

    case SkewType::GadrasGeneric:
    case SkewType::GadrasCZT:
    {
      // SkewPar0=low_skew, SkewPar1=high_skew (fittable amplitudes, GADRAS magnitudes @ 661 keV).
      // SkewPar2/3 = low/high skew power (energy-dependence exponent, fixed detector characteristic).
      // SkewPar4/5 = low/high skew extent (tail slope shaping, fixed detector characteristic).
      switch( coef )
      {
        case CoefficientType::SkewPar0: // low_skew amplitude
        case CoefficientType::SkewPar1: // high_skew amplitude
          starting_value = 5.0;
          step_size = 0.5;
          lower_value = 0.0;
          upper_value = 100.0;
          break;

        case CoefficientType::SkewPar2: // low_skew_power
        case CoefficientType::SkewPar3: // high_skew_power
          starting_value = 0.0;
          step_size = 0.1;
          lower_value = 0.0;
          upper_value = 5.0;
          break;

        case CoefficientType::SkewPar4: // low_skew_extent
        case CoefficientType::SkewPar5: // high_skew_extent
          starting_value = 0.0;
          step_size = 0.5;
          lower_value = -10.0;
          upper_value = 10.0;
          break;

        default:
          return false;
      }//switch( coef )

      break;
    }//case SkewType::GadrasGeneric, SkewType::GadrasCZT:
  }//switch( skew_type )

  return true;
}//void skew_parameter_range(...)


bool PeakDef::skew_no_skew_value( const SkewType skew_type, const CoefficientType coef,
                                  double &no_skew_value )
{
  no_skew_value = 0.0;

  // Only meaningful for coefficients this skew type actually uses; piggy-back on the bounds
  //  function both to reject inapplicable coefficients and to read the (possibly widened) bounds
  //  so the asymptotic "no skew" value tracks `skew_parameter_range` automatically.
  double lower = 0.0, upper = 0.0, starting = 0.0, step = 0.0;
  if( !skew_parameter_range( skew_type, coef, lower, upper, starting, step ) )
    return false;

  switch( skew_type )
  {
    case NumSkewType:
    case NoSkew:
      return false;

    case SkewType::Bortel:
      // tau -> 0 is exactly a pure Gaussian (bortel_indefinite_integral returns pure erf for skew<=0).
      no_skew_value = lower;  // 0.0
      return true;

    case SkewType::GaussExp:
    case SkewType::ExpGaussExp:
      // skew -> inf approaches a pure Gaussian; the (negligible-tail) upper bound is the proxy.
      no_skew_value = upper;
      return true;

    case SkewType::CrystalBall:
    case SkewType::DoubleSidedCrystalBall:
      // alpha -> inf approaches a pure Gaussian; the power `n` is then irrelevant, so leave it neutral.
      if( (coef == CoefficientType::SkewPar0) || (coef == CoefficientType::SkewPar2) )
        no_skew_value = upper;     // alpha (left / right)
      else
        no_skew_value = starting;  // n (does not matter once alpha is at "no skew")
      return true;

    case SkewType::GaussPlusBortel:
      // R == 0 is a pure Gaussian; tau is then irrelevant.
      no_skew_value = (coef == CoefficientType::SkewPar0) ? lower : starting;
      return true;

    case SkewType::VoigtPlusBortel:
      // gamma_lor == 0 and R == 0 give a pure Gaussian; tau is then irrelevant.
      no_skew_value = ((coef == CoefficientType::SkewPar0) || (coef == CoefficientType::SkewPar1))
                        ? lower : starting;
      return true;

    case SkewType::DoubleBortel:
      // Always a sum of two exponential-Gaussian tails - there is no pure-Gaussian limit, so do not
      //  offer a skew removal for it.
      return false;

    case SkewType::GadrasGeneric:
    case SkewType::GadrasCZT:
      // GADRAS shapes carry an intrinsic, energy-dependent skew (the low/high tail amplitudes are
      //  fixed detector characteristics, not a generic Gaussian-plus-skew add-on), and do not use
      //  InterSpec's generic skew machinery - so there is no pure-Gaussian value to snap to.  Do
      //  not offer generic skew removal for them.
      return false;
  }//switch( skew_type )

  return false;
}//bool skew_no_skew_value(...)


size_t PeakDef::num_skew_parameters( const SkewType skew_type )
{
  switch ( skew_type )
  {
    case NumSkewType:
      assert( 0 );
      //throw std::logic_error( "NumSkewType is not a valid skew type" );
      // Fall through to NoSkew for non-debug builds

    case NoSkew:                 return 0;
    case Bortel:                 return 1;
    case DoubleBortel:           return 3; // tau1, tau2_delta, eta
    case GaussPlusBortel:        return 2; // R, tau
    case GaussExp:               return 1;
    case CrystalBall:            return 2;
    case ExpGaussExp:            return 2;
    case DoubleSidedCrystalBall: return 4;
    case VoigtPlusBortel:        return 3; // gamma_lor, R, tau
    case GadrasGeneric:          return 6; // low_skew, high_skew, low/high power, low/high extent
    case GadrasCZT:              return 6; // low_skew, high_skew, low/high power, low/high extent
  }//switch ( skew_type )

  assert( 0 );
  throw std::logic_error( "Invalid peak skew" );
  return 0;
}//size_t num_skew_parameters( const SkewType skew_type );


bool PeakDef::is_energy_dependent( const SkewType skew_type, const CoefficientType coef )
{
  switch( skew_type )
  {
    case NumSkewType:
      assert( 0 );
      //throw std::logic_error( "NumSkewType is not a valid skew type" );
      // Fall through to NoSkew for non-debug builds

    case NoSkew:
      return false;

    case SkewType::Bortel:
      assert( coef == CoefficientType::SkewPar0 );
      return (coef == CoefficientType::SkewPar0);

    case SkewType::DoubleBortel:
      // tau1 (SkewPar0) and tau2_delta (SkewPar1) are energy dependent
      // eta (SkewPar2) may or may not be energy dependent
      return ((coef == CoefficientType::SkewPar0) || (coef == CoefficientType::SkewPar1));

    case SkewType::GaussPlusBortel:
      // Both R (SkewPar0) and tau (SkewPar1) are energy dependent
      return ((coef == CoefficientType::SkewPar0) || (coef == CoefficientType::SkewPar1));

    case SkewType::GaussExp:
      assert( coef == CoefficientType::SkewPar0 );
      return (coef == CoefficientType::SkewPar0);

    case SkewType::CrystalBall:
      assert( (coef == CoefficientType::SkewPar0) || (coef == CoefficientType::SkewPar1) );
      return (coef == CoefficientType::SkewPar0);

    case SkewType::ExpGaussExp:
      assert( (coef == CoefficientType::SkewPar0) || (coef == CoefficientType::SkewPar1) );
      return ((coef == CoefficientType::SkewPar0) || (coef == CoefficientType::SkewPar1));

    case SkewType::DoubleSidedCrystalBall:
      return ((coef == CoefficientType::SkewPar0) || (coef == CoefficientType::SkewPar2));

    case SkewType::VoigtPlusBortel:
      // gamma_lor (SkewPar0) may be energy dependent if from lookup table
      // R (SkewPar1) and tau (SkewPar2) are energy dependent for RelActAuto
      return ((coef == CoefficientType::SkewPar0) || (coef == CoefficientType::SkewPar1)
              || (coef == CoefficientType::SkewPar2));

    case SkewType::GadrasGeneric:
    case SkewType::GadrasCZT:
      // The GADRAS shape applies its energy dependence intrinsically (via the skew "power" terms and
      // each peak's own energy), so we do NOT use InterSpec's generic energy-dependent-skew machinery -
      // doing so would double-count the energy dependence.
      return false;
  }//switch( skew_type )

  assert( 0 );
  return false;
}//bool is_energy_dependent( const SkewType skew_type, const CoefficientType coefficient )


bool PeakDef::skew_parameter_fit_by_default( const SkewType skew_type, const CoefficientType coef )
{
  switch( skew_type )
  {
    case SkewType::GadrasGeneric:
    case SkewType::GadrasCZT:
      // Only the two amplitude parameters are fit by default; the powers/extents are fixed
      // detector characteristics.
      return ((coef == CoefficientType::SkewPar0) || (coef == CoefficientType::SkewPar1));

    case NoSkew:
    case Bortel:
    case GaussExp:
    case CrystalBall:
    case ExpGaussExp:
    case DoubleSidedCrystalBall:
    case VoigtPlusBortel:
    case GaussPlusBortel:
    case DoubleBortel:
    case NumSkewType:
      break;
  }//switch( skew_type )

  // Default: any parameter within the skew type's parameter count is fit.
  const size_t par_num = static_cast<size_t>(coef) - static_cast<size_t>(CoefficientType::SkewPar0);
  return ((coef >= CoefficientType::SkewPar0) && (coef <= CoefficientType::SkewPar5)
          && (par_num < num_skew_parameters(skew_type)));
}//bool skew_parameter_fit_by_default(...)


double PeakDef::extract_energy_from_peak_source_string( std::string &str )
{  
  std::smatch energy_match;
  const std::regex energy_regexp("(?:^|\\s)((((\\d+(\\.\\d*)?)|(\\.\\d*))\\s*(?:[Ee][+\\-]?\\d+)?)\\s*(kev|mev|ev|$))",
                                 regex::ECMAScript | regex::icase );
  
  if( !std::regex_search(str, energy_match, energy_regexp) )
    return -1.0;
  
  
  const string &val = energy_match[2];
  const string &units = energy_match[7];
  const string &total_match = energy_match[0];
  
  //cout << "Match for '"  << str << "': ";
  //for( const auto i : energy_match )
  //  cout << "'" << i << "', ";
  //cout << ", val='" << val << ", units='" << units << "'" << endl;
  
  double energy = -1.0;
  if( !(stringstream(val) >> energy) )
  {
    assert( 0 ); //should ever get here if regex is well formed
    return -1.0;
  }
  
  if( SpecUtils::iequals_ascii(units, "kev") )
  {
    // Nothing to do here
  }else if( SpecUtils::iequals_ascii(units, "mev") )
  {
    energy *= 1000.0;
  }else if( SpecUtils::iequals_ascii(units, "ev") )
  {
    energy /= 1000.0;
  }else if( units.empty() )
  {
    if( (energy < 295.0) && (str.find('.') == string::npos) )
    {
      // If value is less than 295, and there is no decimal, then its possible we've picked up on
      //  isotope number (e.g., the 235 from U235), so for the moment, we'll reject this match,
      //  unless we unambiguously know there were numbers or reaction before what we think is the
      //  energy
      const auto matchpos = str.find(total_match);
      auto first_num_pos = str.find_first_of( "0123456789)" );
      //if( first_num_pos > matchpos )
      //  first_num_pos = str.find_first_of( "x-ray" );
      //if( first_num_pos > matchpos )
      //  first_num_pos = str.find_first_of( "xray" );
      //if( first_num_pos > matchpos )
      //  first_num_pos = str.find_first_of( "x ray" );
      
      if( first_num_pos < matchpos )
      {
        // We have an x-ray, or there were numbers, or closing parenthesis before our match
      }else
      {
        return -1.0;
      }
    }//if( (energy < 295.0) && (str.find('.') == string::npos) )
  }
  
  SpecUtils::ireplace_all( str, total_match.c_str(), " " );
  SpecUtils::trim( str );
  
  return energy;
};//extract_energy_from_peak_source_string


const char *PeakDef::to_str( const DefintionType type )
{
  switch( type )
  {
    case GaussianDefined:
      return "GaussianDefined";
    case DataDefined:
      return "DataDefined";
  }//switch( type )
  
  assert( 0 );
  return "InvalidDefintionType";
}//const char *to_str( const DefintionType type )


PeakDef::DefintionType PeakDef::peak_type_from_str( const char * const str )
{
  if( str && SpecUtils::icontains( str, "GaussianDefined" ) )
    return PeakDef::DefintionType::GaussianDefined;
  
  if( str && SpecUtils::icontains( str, "DataDefined" ) )
    return PeakDef::DefintionType::DataDefined;
  
  throw runtime_error( "Invalid PeakDef::DefintionType: '" + string(str ? str : "") );
}//DefintionType peak_type_from_str( const char * const str )


void PeakDef::gammaTypeFromUserInput( std::string &txt,
                                      PeakDef::SourceGammaType &type )
{
  
  type = PeakDef::NormalGamma;
  
  if( SpecUtils::icontains( txt, "s.e." ) )
  {
    type = PeakDef::SingleEscapeGamma;
    SpecUtils::ireplace_all( txt, "s.e.", "" );
  }else if( SpecUtils::icontains( txt, "single escape" ) || SpecUtils::icontains( txt, "single-escape" ) )
  {
    type = PeakDef::SingleEscapeGamma;
    SpecUtils::ireplace_all( txt, "single escape", "" );
    SpecUtils::ireplace_all( txt, "single-escape", "" );
  }else if( SpecUtils::iequals_ascii( txt, "s.e." )
     || SpecUtils::iequals_ascii( txt, "se" )
     || SpecUtils::iequals_ascii( txt, "escape" ) )
  {
    type = PeakDef::SingleEscapeGamma;
    txt = "";
  }else if( SpecUtils::icontains( txt, "se " ) && txt.size() > 5 )
  {
    type = PeakDef::SingleEscapeGamma;
    SpecUtils::ireplace_all( txt, "se ", "" );
  }else if( SpecUtils::icontains( txt, " se" ) && txt.size() > 5 )
  {
    type = PeakDef::SingleEscapeGamma;
    SpecUtils::ireplace_all( txt, " se", "" );
  }else if( SpecUtils::icontains( txt, "d.e." ) )
  {
    type = PeakDef::DoubleEscapeGamma;
    SpecUtils::ireplace_all( txt, "d.e.", "" );
  }else if( SpecUtils::icontains( txt, "double escape" ) || SpecUtils::icontains( txt, "double-escape" ) )
  {
    type = PeakDef::DoubleEscapeGamma;
    SpecUtils::ireplace_all( txt, "double escape", "" );
    SpecUtils::ireplace_all( txt, "double-escape", "" );
  }else if( SpecUtils::icontains( txt, "de " ) && txt.size() > 5 )
  {
    type = PeakDef::DoubleEscapeGamma;
    SpecUtils::ireplace_all( txt, "de ", "" );
  }else if( SpecUtils::icontains( txt, " de" ) && txt.size() > 5 )
  {
    type = PeakDef::DoubleEscapeGamma;
    SpecUtils::ireplace_all( txt, " de", "" );
  }else if( SpecUtils::iequals_ascii( txt, "d.e." )
           || SpecUtils::iequals_ascii( txt, "de" )  )
  {
    type = PeakDef::DoubleEscapeGamma;
    txt = "";
  }else if( SpecUtils::icontains( txt, "x-ray" )
      || SpecUtils::icontains( txt, "xray" )
     || SpecUtils::icontains( txt, "x ray" ) )
  {
    type = PeakDef::XrayGamma;
    SpecUtils::ireplace_all( txt, "xray", "" );
    SpecUtils::ireplace_all( txt, "x-ray", "" );
    SpecUtils::ireplace_all( txt, "x ray", "" );
  }
  
  SpecUtils::trim( txt );
}//PeakDef::SourceGammaType gammaType( std::string txt )


const Wt::WColor &PeakDef::lineColor() const
{
  return m_lineColor;
}

void PeakDef::setLineColor( const Wt::WColor &color )
{
  m_lineColor = color;
}


#if( SpecUtils_ENABLE_D3_CHART )
std::string PeakDef::gaus_peaks_to_json(const std::vector<std::shared_ptr<const PeakDef> > &peaks,
                                        const std::shared_ptr<const SpecUtils::Measurement> &foreground,
                                        const Wt::WColor &defaultPeakColor,
                                        int override_alpha )
{
  //Need to check all numbers to make sure not inf or nan
  
  stringstream answer;
  
  (void)override_alpha; //Light version does not override peak alpha
  
  if (peaks.empty())
    return answer.str();

  std::shared_ptr<const PeakContinuum> continuum = peaks[0]->continuum();
  if (!continuum)
    throw runtime_error("gaus_peaks_to_json: invalid continuum");
  
  
  if( IsInf(continuum->lowerEnergy()) || IsNan(continuum->lowerEnergy()) )
    throw runtime_error( "Continuum lower energy is invalid" );
  
  if( IsInf(continuum->upperEnergy()) || IsNan(continuum->upperEnergy()) )
    throw runtime_error( "Continuum upper energy is invalid" );
  
  const char *q = "\""; // for creating valid json format

  answer << "{" << q << "type" << q << ":" << q;
  switch( continuum->type() )
  {
    case PeakContinuum::NoOffset:     answer << "NoOffset";     break;
    case PeakContinuum::Constant:     answer << "Constant";     break;
    case PeakContinuum::Linear:       answer << "Linear";       break;
    case PeakContinuum::Quadratic:    answer << "Quadratic";    break;
    case PeakContinuum::Cubic:        answer << "Cubic";        break;
    case PeakContinuum::FlatStep:     answer << "FlatStep";     break;
    case PeakContinuum::LinearStep:   answer << "LinearStep";   break;
    case PeakContinuum::BiLinearStep: answer << "BiLinearStep"; break;
    case PeakContinuum::FlatStepCDF:  answer << "FlatStepCDF";  break;
    case PeakContinuum::LinearStepCDF:    answer << "LinearStepCDF";    break;
    case PeakContinuum::BiLinearStepCDF: answer << "BiLinearStepCDF"; break;
    case PeakContinuum::External:        answer << "External";        break;
  }//switch( continuum->type() )
  
  
  // We use the peaks defined range, and not the continuum, as the continuum may not
  //  have the range defined (but normally should).
  //  This next statement assumes peaks are sorted in increasing mean (but this only
  //  matters if ROI range is not defined) - we could check this, but maybe not
  //  worth the overhead for the edge case, that we might never encounter.
  answer << q << "," << q << "lowerEnergy" << q << ":" << peaks.front()->lowerX()
         << "," << q << "upperEnergy" << q << ":" << peaks.back()->upperX();
  
  
  if( foreground && foreground->channel_energies() && foreground->channel_energies()->size() > 2 )
  {
    // If we want integer counts (for simple spectra only), should use gamma_channels_sum, otherwise
    //   gamma_integral(...) will almost always return fractional.
    //const size_t lower_channel = foreground->find_gamma_channel( continuum->lowerEnergy() );
    //const size_t upper_channel = foreground->find_gamma_channel( continuum->upperEnergy() );
    //const double sum = foreground->gamma_channels_sum( lower_channel, upper_channel );
    const double sum = foreground->gamma_integral( continuum->lowerEnergy(), continuum->upperEnergy() );
    answer << "," << q << "roiCounts" << q << ":" << sum;
  }//if( foreground )
  
  
  switch( continuum->type() )
  {
    case PeakContinuum::NoOffset:
      break;
      
    case PeakContinuum::Constant:
    case PeakContinuum::Linear:
    case PeakContinuum::Quadratic:
    case PeakContinuum::Cubic:
    case PeakContinuum::FlatStep:
    case PeakContinuum::LinearStep:
    case PeakContinuum::BiLinearStep:
    case PeakContinuum::FlatStepCDF:
    case PeakContinuum::LinearStepCDF:
    case PeakContinuum::BiLinearStepCDF:
    {
      if( IsInf(continuum->referenceEnergy()) || IsNan(continuum->referenceEnergy()) )
        throw runtime_error( "Continuum reference energy is invalid" );

      answer << "," << q << "referenceEnergy" << q << ":" << continuum->referenceEnergy();
      const vector<double> &values = continuum->parameters();
      const vector<double> &uncerts = continuum->uncertainties();
      answer << "," << q << "coeffs" << q << ":[";
      for (size_t i = 0; i < values.size(); ++i)
      {
        if( IsInf(values[i]) || IsNan(values[i]) )
          throw runtime_error( "Continuum coef is invalid" );

        answer << (i ? "," : "") << values[i];
      }
      answer << "]," << q << "coeffUncerts" << q << ":[";
      for (size_t i = 0; i < uncerts.size(); ++i)
        answer << (i ? "," : "") << ((IsInf(uncerts[i]) || IsNan(uncerts[i])) ? -1.0 : uncerts[i]); //we'll let uncertainties slide since we dont use them
      answer << "]";

      answer << "," << q << "fitForCoeff" << q << ":[";
      for (size_t i = 0; i < continuum->fitForParameter().size(); ++i)
        answer << (i ? "," : "") << (continuum->fitForParameter()[i] ? "true" : "false");
      answer << "]";
      
      if( PeakContinuum::is_step_continuum( continuum->type() )
         && foreground && foreground->num_gamma_channels() )
      {
        const size_t nchannel = foreground->num_gamma_channels();
        
        //We'll put in the coefficients, but also the values to make things easy in JS
        size_t firstbin = foreground->find_gamma_channel( continuum->lowerEnergy() );
        size_t lastbin = foreground->find_gamma_channel( continuum->upperEnergy() );
        
        firstbin = (firstbin > 0) ? (firstbin - 1) : firstbin;
        firstbin = (firstbin > 0) ? (firstbin - 1) : firstbin;
        lastbin = (lastbin < (nchannel - 1)) ? (lastbin + 1) : lastbin;
        lastbin = (lastbin < (nchannel - 1)) ? (lastbin + 1) : lastbin;
        
        assert( firstbin <= lastbin );
        assert( lastbin <= (nchannel - 1) );
        
        // When the JSON of the spectrum chart is defined, D3SpectrumExport.cpp/write_spectrum_data_js(...)
        //  will send the energy calibration as coefficients sometimes, and lower channel energies other
        //  times.  If coefficients are sent, then the JS computes channel bounds using doubles.
        //  And somewhat surprisingly, rounding causes visual artifacts of the continuum when
        //  the 'continuumEnergies' and 'continuumCounts' arrays are used, which are only accurate
        //  to float levels.  So we will increase accuracy of the 'answer' stream here, which appears
        //  to be enough to avoid these artifacts, but also commented out is how we could compute
        //  to double precision to match what happens in the JS.
        const auto oldprecision = answer.precision();  //probably 6 always
        answer << std::setprecision(std::numeric_limits<float>::digits10 + 1);
        
        answer << "," << q << "continuumEnergies" << q << ":[";
        for (size_t i = firstbin; i <= lastbin; ++i)
          answer << ((i!=firstbin) ? "," : "") << foreground->gamma_channel_lower(i);
        answer << "]," << q << "continuumCounts" << q << ":[";
        
        /*
         //Implementation to compute energies to double precision - doesnt appear to be necessary.
        const SpecUtils::EnergyCalType caltype = foreground->energy_calibration_model();
        if(  (caltype == SpecUtils::EnergyCalType::Polynomial
             || caltype == SpecUtils::EnergyCalType::FullRangeFraction
             || caltype == SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial )
           && foreground->deviation_pairs().empty() )
        {
          const size_t nchannel = foreground->num_gamma_channels();
          const vector<float> &coefs = foreground->calibration_coeffs();
          const vector<pair<float,float>> &dev_pairs = foreground->deviation_pairs();
         
          answer << "," << q << "continuumEnergies" << q << ":[";
          for (size_t i = firstbin; i <= lastbin; ++i)
          {
            double energy;
            if( caltype == SpecUtils::EnergyCalType::FullRangeFraction )
              energy = SpecUtils::fullrangefraction_energy( i, coefs, nchannel, dev_pairs );
            else
              energy = SpecUtils::polynomial_energy( i, coefs, dev_pairs );
            answer << ((i!=firstbin) ? "," : "") << energy;
          }
          answer << "]," << q << "continuumCounts" << q << ":[";
        }else
        {
          answer << "," << q << "continuumEnergies" << q << ":[";
          for (size_t i = firstbin; i <= lastbin; ++i)
            answer << ((i!=firstbin) ? "," : "") << foreground->gamma_channel_lower(i);
          answer << "]," << q << "continuumCounts" << q << ":[";
        }
         */
        
        //Slow implementation pre 20251105
        //for (size_t i = firstbin; i <= lastbin; ++i)
        //{
        //  const float lower_x = foreground->gamma_channel_lower( i );
        //  const float upper_x = foreground->gamma_channel_upper( i );
        //  const float cont_counts = continuum->offset_integral( lower_x, upper_x, foreground );
        //  answer << ((i!=firstbin) ? "," : "") << cont_counts;
        //}
        
        {// Begin write continuum counts for each channel
          const shared_ptr<const vector<float>> &channel_energies_ptr = foreground->channel_energies();
          assert( channel_energies_ptr && (lastbin < channel_energies_ptr->size()) && (firstbin <= lastbin) );
          const size_t num_roi_channels = (lastbin > firstbin) ? (1 + lastbin - firstbin) : size_t(1);
          assert( num_roi_channels > 0 );

          const float * const energies = &(channel_energies_ptr->at(firstbin));
          vector<double> continuum_counts( nchannel, 0.0 );

          continuum->offset_integral( energies, &(continuum_counts[0]), num_roi_channels, foreground, peaks );

          for( size_t i = 0; i < num_roi_channels; ++i )
          {
            answer << (i ? "," : "") << continuum_counts[i];
#if( PERFORM_DEVELOPER_CHECKS )
            {
              // The batch and single-channel implementations must agree for every continuum type,
              //  including the peak-CDF step ones - they share `cdf_step_anchor_energies(...)`, so
              //  a divergence here means one of them drifted.
              const size_t channel = i + firstbin;
              const float lower_x = foreground->gamma_channel_lower( channel );
              const float upper_x = foreground->gamma_channel_upper( channel );
              const double cont_counts_test = continuum->offset_integral( lower_x, upper_x, foreground, peaks );
              assert( fabs(continuum_counts[i] - cont_counts_test) < std::max( 0.0001*std::max( fabs(continuum_counts[i]), fabs(cont_counts_test) ), 1.0E-6) );
            }
#endif
          }
        }// End write continuum counts for each channel
        
        answer << "]";
        
        answer << std::setprecision(9);
      }//if( continuum->type() == FlatStep/LinearStep/BiLinearStep )
      
      break;
    }//polynomial continuum
      
    case PeakContinuum::External:
    {
      if( continuum->externalContinuum()
          && continuum->externalContinuum()->num_gamma_channels() )
      {
        std::shared_ptr<const Measurement> hist = continuum->externalContinuum();
        size_t firstbin = hist->find_gamma_channel(continuum->lowerEnergy());
        size_t lastbin = hist->find_gamma_channel(continuum->upperEnergy());
        
        const size_t nchannel = hist->num_gamma_channels();
        firstbin = (firstbin > 0) ? (firstbin - 1) : firstbin;
        firstbin = (firstbin > 0) ? (firstbin - 1) : firstbin;
        lastbin = (lastbin < (nchannel - 1)) ? (lastbin + 1) : nchannel;
        lastbin = (lastbin < (nchannel - 1)) ? (lastbin + 1) : nchannel;
        
        //see comments above in FlatStep section
        const auto oldprecision = answer.precision();
        answer << std::setprecision(std::numeric_limits<float>::digits10 + 1);

        answer << "," << q << "continuumEnergies" << q << ":[";
        for (size_t i = firstbin; i <= lastbin; ++i)
          answer << ((i!=firstbin) ? "," : "") << hist->gamma_channel_lower(i);
        answer << "]," << q << "continuumCounts" << q << ":[";
        for (size_t i = firstbin; i <= lastbin; ++i)
          answer << ((i!=firstbin) ? "," : "") << hist->gamma_channel_content(i);
        answer << "]";
        
        answer << std::setprecision(9);
      }
    }//case PeakContinuum::External:
  }//switch( continuum->type() )


  answer << "," << q << "peaks" << q << ":[";
  for (size_t i = 0; i < peaks.size(); ++i)
  {
    const PeakDef &p = *peaks[i];
    if (continuum != p.continuum())
      throw runtime_error("gaus_peaks_to_json: peaks all must share same continuum");
    answer << (i ? "," : "") << "{";

    if( !p.userLabel().empty() )
    {
      string label = p.userLabel();
      
      answer << q << "userLabel" << q << ":" << nlohmann::json( label ).dump() << ",";
    }

    const Wt::WColor &peak_color = p.lineColor().isDefault() ? defaultPeakColor : p.lineColor();
    if( !peak_color.isDefault() )
      answer << q << "lineColor" << q << ":" << nlohmann::json( peak_color.cssText() ).dump() << ",";
    
    answer << q << "type" << q << ":" << q << PeakDef::to_str(p.type()) << q << ",";
    
    double dist_norm = 0.0;
    answer << q << "skewType" << q << ":";
    
    switch( p.skewType() )
    {
      case NumSkewType:
        assert( 0 );
        //throw runtime_error( "PeakDef JSON serialization: NumSkewType is not a valid skew type" );
        // Fall through to NoSkew for non-debug builds

      case PeakDef::NoSkew:
        answer << q << "NoSkew" << q << ",";
        break;
      
      case Bortel:
        answer << q << "ExGauss" << q << ",";
        break;
        
      case CrystalBall:
        answer << q << "CB" << q << ",";
        dist_norm = PeakDists::crystal_ball_norm( p.coefficient(CoefficientType::Sigma),
                                          p.coefficient(CoefficientType::SkewPar0),
                                      p.coefficient(CoefficientType::SkewPar1) );
        break;
        
      case DoubleSidedCrystalBall:
        answer << q << "DSCB" << q << ",";
        dist_norm = PeakDists::DSCB_norm( p.coefficient(CoefficientType::SkewPar0),
                                                   p.coefficient(CoefficientType::SkewPar1),
                                                   p.coefficient(CoefficientType::SkewPar2),
                                                   p.coefficient(CoefficientType::SkewPar3) );
        
        break;
        
      case GaussExp:
        answer << q << "GaussExp" << q << ",";
        dist_norm = PeakDists::gauss_exp_norm( p.coefficient(CoefficientType::Sigma),
                                   p.coefficient(CoefficientType::SkewPar0) );
        break;
        
      case ExpGaussExp:
        answer << q << "ExpGaussExp" << q << ",";
        dist_norm = PeakDists::exp_gauss_exp_norm( p.coefficient(CoefficientType::Sigma),
                                       p.coefficient(CoefficientType::SkewPar0),
                                       p.coefficient(CoefficientType::SkewPar1) );
        break;
        
      case VoigtPlusBortel:
        answer << q << "VoigtBortel" << q << ",";
        dist_norm = 1.0;  //Distribution should already be normed
        break;

      case GaussPlusBortel:
        answer << q << "GaussBortel" << q << ",";
        dist_norm = 1.0;  //Distribution should already be normed
        break;

      case DoubleBortel:
        answer << q << "DoubleBortel" << q << ",";
        dist_norm = 1.0;  //Distribution should already be normed
        break;

      case GadrasGeneric:
        answer << q << "GadrasGeneric" << q << ",";
        dist_norm = 1.0;  //Distribution should already be normed
        break;

      case GadrasCZT:
        answer << q << "GadrasCZT" << q << ",";
        dist_norm = 1.0;  //Distribution should already be normed
        break;
    }//switch( p.type() )

    if( (p.skewType() != PeakDef::NoSkew) && (!IsNan(dist_norm) && !IsInf(dist_norm)) )
      answer << q << "DistNorm" << q << ":" << dist_norm << ",";
    
    if( (p.type() == PeakDef::GaussianDefined) && (p.skewType() != PeakDef::NoSkew) )
    {
      double hidden_frac = 1.0E-6; //
      
      try
      {
        pair<double,double> vis_limits;
        
        switch( p.skewType() )
        {
          case NumSkewType:
            assert( 0 );
            //throw runtime_error( "PeakDef coverage limits: NumSkewType is not a valid skew type" );
            // Fall through to NoSkew for non-debug builds

          case NoSkew:
            vis_limits.first = p.mean() - 5.0*p.sigma();
            vis_limits.second = p.mean() + 5.0*p.sigma();
            break;
            
          case Bortel:
            vis_limits = PeakDists::bortel_coverage_limits( p.mean(), p.sigma(),
                                                           p.coefficient(CoefficientType::SkewPar0),
                                                           hidden_frac );
            break;
            
          case GaussExp:
            vis_limits = PeakDists::gauss_exp_coverage_limits( p.mean(), p.sigma(),
                                                              p.coefficient(CoefficientType::SkewPar0),
                                                              hidden_frac );
            break;
          case CrystalBall:
            try
            {
              vis_limits = PeakDists::crystal_ball_coverage_limits( p.mean(), p.sigma(),
                                                                 p.coefficient(CoefficientType::SkewPar0),
                                                                 p.coefficient(CoefficientType::SkewPar1),
                                                                 hidden_frac );
            }catch( std::exception & )
            {
              // CB dist can have really long tail, causing the coverage limits to fail, because
              //  of unreasonable values - in this case we'll use the entire ROI.
              vis_limits.first = p.lowerX();
              vis_limits.second = p.upperX();
            }
            break;
          case ExpGaussExp:
            vis_limits = PeakDists::exp_gauss_exp_coverage_limits( p.mean(), p.sigma(),
                                                                  p.coefficient(CoefficientType::SkewPar0),
                                                                  p.coefficient(CoefficientType::SkewPar1),
                                                                  hidden_frac );
            break;
          case DoubleSidedCrystalBall:
            // For largely skewed peaks, going out t 1E-6 is unreasonable, so we could limit
            //  this to 1E-3, to improve chances of success.
            //if( p.coefficient(CoefficientType::SkewPar1) < 5.0
            //   || p.coefficient(CoefficientType::SkewPar3) < 5.0)
            //  hidden_frac = 1.0E-3; //DSCB doesnt behave well for power-law peaks...
            
            try
            {
              const double mean = p.mean();
              const double sigma = p.sigma();
              const double left_skew = p.coefficient(CoefficientType::SkewPar0);
              const double left_n = p.coefficient(CoefficientType::SkewPar1);
              const double right_skew = p.coefficient(CoefficientType::SkewPar2);
              const double right_n = p.coefficient(CoefficientType::SkewPar3);
              
              vis_limits = PeakDists::double_sided_crystal_ball_coverage_limits( mean, sigma,
                                              left_skew, left_n, right_skew, right_n, hidden_frac );
            }catch( std::exception & )
            {
              // DSCB dist can have really long tail, causing the coverage limits to fail, because
              //  of unreasonable values - in this case we'll use the entire ROI.  This throw is the
              //  expected, cheap signal for a heavy tail (a power-law `n` near its 1.05 bound), not
              //  a numerical accident - see `PeakDists::sk_max_coverage_limit_nsigma`.
              vis_limits.first = p.lowerX();
              vis_limits.second = p.upperX();
            }//try / catch
            break; //case DoubleSidedCrystalBall:
            
            
          case VoigtPlusBortel:
            try
            {
              const double mean = p.mean();
              const double sigma = p.sigma();
              const double gamma_lor = p.coefficient(CoefficientType::SkewPar0);
              const double R = p.coefficient(CoefficientType::SkewPar1);
              const double tau = p.coefficient(CoefficientType::SkewPar2);

              vis_limits = PeakDists::voigt_exp_coverage_limits( mean, sigma, gamma_lor,
                                                                 R, tau, hidden_frac );
            }catch( std::exception & )
            {
              // Voigt plus Bortel dist can have really long tail, causing the coverage limits to fail,
              //  because of unreasonable values - in this case we'll use the entire ROI.
              vis_limits.first = p.lowerX();
              vis_limits.second = p.upperX();
            }//try / catch
            break; //case VoigtPlusBortel:

          case GaussPlusBortel:
            try
            {
              const double mean = p.mean();
              const double sigma = p.sigma();
              const double R = p.coefficient(CoefficientType::SkewPar0);
              const double tau = p.coefficient(CoefficientType::SkewPar1);

              vis_limits = PeakDists::gauss_plus_bortel_coverage_limits( mean, sigma, R, tau, hidden_frac );
            }catch( std::exception & )
            {
              vis_limits.first = p.lowerX();
              vis_limits.second = p.upperX();
            }//try / catch
            break; //case GaussPlusBortel:

          case DoubleBortel:
            try
            {
              const double mean = p.mean();
              const double sigma = p.sigma();
              const double tau1 = p.coefficient(CoefficientType::SkewPar0);
              const double tau2_delta = p.coefficient(CoefficientType::SkewPar1);
              const double eta = p.coefficient(CoefficientType::SkewPar2);

              vis_limits = PeakDists::double_bortel_coverage_limits( mean, sigma, tau1, tau2_delta, eta, hidden_frac );
            }catch( std::exception & )
            {
              vis_limits.first = p.lowerX();
              vis_limits.second = p.upperX();
            }//try / catch
            break; //case DoubleBortel:

          case GadrasGeneric:
          case GadrasCZT:
            try
            {
              const PeakDists::GadrasMaterial mat = (p.skewType() == GadrasCZT)
                                    ? PeakDists::GadrasMaterial::CZT_CdTe
                                    : PeakDists::GadrasMaterial::Generic;
              double skew[6];
              for( int i = 0; i < 6; ++i )
                skew[i] = p.coefficient( CoefficientType(CoefficientType::SkewPar0 + i) );

              vis_limits = PeakDists::gadras_coverage_limits( p.mean(), p.sigma(), skew, mat, hidden_frac );
            }catch( std::exception & )
            {
              vis_limits.first = p.lowerX();
              vis_limits.second = p.upperX();
            }//try / catch
            break; //case GadrasGeneric, GadrasCZT:
        }//switch( p.skewType() )
        
        answer << q << "visRange" << q << ":[" << vis_limits.first << "," << vis_limits.second << "],";
      }catch( std::exception &e )
      {
        cerr << "Failed to get limits for peak type " << PeakDef::to_string(p.skewType())
            << ": " << e.what() << endl
            << "\tFor peak with skew=[" << p.coefficient(CoefficientType::SkewPar0)
            << ", " << p.coefficient(CoefficientType::SkewPar1)
            << ", " << p.coefficient(CoefficientType::SkewPar2)
            << ", " << p.coefficient(CoefficientType::SkewPar3)
            << "], mean=" << p.mean() << ", and sigma=" << p.sigma()
            << " with prob=" << hidden_frac << endl;
      }//try / catch
    }//if( we should give a hint to the JS about range to draw skewed peaks )
    
    
    
    for (PeakDef::CoefficientType t = PeakDef::CoefficientType(0);
      t < PeakDef::NumCoefficientTypes; t = PeakDef::CoefficientType(t + 1))
    {
      // Dont include non-useful skew coefficients
      switch( t )
      {
        case SkewPar0:
        case SkewPar1:
        case SkewPar2:
        case SkewPar3:
        case SkewPar4:
        case SkewPar5:
        {
          double lower, upper, starting, step;
          if( !PeakDef::skew_parameter_range( p.m_skewType, t, lower, upper, starting, step ) )
            continue;

          break;
        }//case A Skew Parameter
          
        default:
          break;
      }//switch( t )
      
      double coef = p.coefficient(t), uncert = p.uncertainty(t);
      assert( !IsInf(coef) && !IsNan(coef) || (t == PeakDef::Chi2DOF) ); // PeakDef::Chi2DOF may be +inf, if it wasnt evaluated
      if( IsInf(coef) || IsNan(coef) )
      {
        //throw runtime_error( "Peak ceoff is inf or nan" );
        coef = 0.0;
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, "Peak coefficient is Inf or NaN" );
#endif//#if( PERFORM_DEVELOPER_CHECKS )
      }
      
      if( IsInf(uncert) || IsNan(uncert) )
        uncert = 0.0;
      
      answer << q << PeakDef::to_string(t) << q << ":[" << coef
        << "," << uncert << "," << (p.fitFor(t) ? "true" : "false")
        << "],";
    }//for(...)
    
    // The peak amplitude CPS, is only used for the little info box when you hover-over/tap a peak,
    //   so we'll just format the text here; may change in the future.
    if( (p.type() == PeakDef::GaussianDefined) && (p.amplitudeUncert() > 0.0f)
       && foreground && (foreground->live_time() > 0.0f) )
    {
      const float lt = foreground->live_time();
      const float cps = p.amplitude() / lt;
      const float cpsUncert = p.amplitudeUncert() / lt;
      const string uncertstr = PhysicalUnits::printValueWithUncertainty( cps, cpsUncert, 4 );
      answer << q << "cpsTxt" << q << ":" << q << uncertstr << q << ",";
    }//if( we can give peak CPS )
    
    answer << q << "useForEnergyCalibration" << q << ":" << (p.useForEnergyCalibration() ? "true" : "false")
           << "," << q << "forSourceFit" << q << ":" << (p.useForShieldingSourceFit() ? "true" : "false");

    if( p.useForDrfIntrinsicEffFit() != PeakDef::sm_defaultUseForDrfIntrinsicEffFit )
      answer << "," << q << "useForDrfIntrinsicEffFit" << q << ":" << (p.useForDrfIntrinsicEffFit() ? "true" : "false");
    
    if( p.useForDrfFwhmFit() != PeakDef::sm_defaultUseForDrfFwhmFit )
      answer << "," << q << "useForDrfFwhmFit" << q << ":" << (p.useForDrfFwhmFit() ? "true" : "false");
    
    if( p.useForDrfDepthOfInteractionFit() != PeakDef::sm_defaultUseForDrfDepthOfInteractionFit)
      answer << "," << q << "useForDrfDepthOfInteractionFit" << q << ":" << (p.useForDrfDepthOfInteractionFit() ? "true" : "false");
    
    
    /*
    const char *gammaTypeVal = 0;
    switch( p.sourceGammaType() )
    {
      case PeakDef::NormalGamma:       gammaTypeVal = "NormalGamma";       break;
      case PeakDef::AnnihilationGamma: gammaTypeVal = "AnnihilationGamma"; break;
      case PeakDef::SingleEscapeGamma: gammaTypeVal = "SingleEscapeGamma"; break;
      case PeakDef::DoubleEscapeGamma: gammaTypeVal = "DoubleEscapeGamma"; break;
      case PeakDef::XrayGamma:         gammaTypeVal = "XrayGamma";         break;
    }//switch( p.sourceGammaType() )

    if( p.parentNuclide() || p.xrayElement() || p.reaction() )
      answer << "," << q << "sourceType" << q << ":" << q << gammaTypeVal << q;
    */
    
    const auto srcTypeJSON = [q]( const PeakDef::SourceGammaType gamtype ) -> std::string{
      const char *type = "";
      switch( gamtype )
      {
        case NormalGamma:       return ""; break;
        case AnnihilationGamma: type = "annih."; break;
        case XrayGamma:         type = "x-ray";  break;
        case DoubleEscapeGamma: type = "D.E.";   break;
        case SingleEscapeGamma: type = "S.E.";   break;
      }//switch( p.sourceGammaType() )
      
      return string(",") + q + "type" + q + ":" + q + type + q;
    };//srcTypeJSON lamda
    
    
    const PeakDef::Source &src = p.source();
    const float src_energy = (IsInf(src.particle_energy) || IsNan(src.particle_energy)) ? 0.0f : src.particle_energy;
    switch( src.kind )
    {
      case PeakDef::Source::Kind::None:
        break;

      case PeakDef::Source::Kind::Nuclide:
        answer << "," << q << "nuclide" << q << ":{" << q << "name" << q << ":" << nlohmann::json( src.name ).dump();
        if( !src.decay_parent.empty() )
          answer << "," << q << "decayParent" << q << ":" << nlohmann::json( src.decay_parent ).dump();
        answer << "," << q << "energy" << q << ":" << src_energy << srcTypeJSON( src.gamma_type ) << "}";
        break;

      case PeakDef::Source::Kind::Xray:
      case PeakDef::Source::Kind::Reaction:
        answer << "," << q << ((src.kind == PeakDef::Source::Kind::Xray) ? "xray" : "reaction") << q
               << ":{" << q << "name" << q << ":" << nlohmann::json( src.name ).dump() << ","
               << q << "energy" << q << ":" << src_energy << srcTypeJSON( src.gamma_type ) << "}";
        break;
    }//switch( src.kind )

    answer << "}";
  }//for( size_t i = 0; i < peaks.size(); ++i )

  answer << "]}";

  return answer.str();
}//std::string PeakDef::gaus_peaks_to_json(...)


string PeakDef::peak_json(const vector<std::shared_ptr<const PeakDef> > &inpeaks,
                          const std::shared_ptr<const SpecUtils::Measurement> &foreground,
                          const Wt::WColor &defaultPeakColor,
                          int override_alpha )
{
  if (inpeaks.empty())
    return "[]";

  typedef std::map< std::shared_ptr<const PeakContinuum>, vector<std::shared_ptr<const PeakDef> > > ContinuumToPeakMap_t;

  ContinuumToPeakMap_t continuumToPeaks;
  for (size_t i = 0; i < inpeaks.size(); ++i)
    continuumToPeaks[inpeaks[i]->continuum()].push_back(inpeaks[i]);

  // Output ROIs in energy order (the map is ordered by pointer, which is not reproducible)
  vector<const ContinuumToPeakMap_t::value_type *> rois;
  for( const ContinuumToPeakMap_t::value_type &vt : continuumToPeaks )
    rois.push_back( &vt );
  std::sort( begin(rois), end(rois), []( const auto *lhs, const auto *rhs ){
    return lhs->second.front()->mean() < rhs->second.front()->mean();
  } );

  string json = "[";
  for( const ContinuumToPeakMap_t::value_type *vt : rois )
    json += ((json.size()>2) ? "," : "") + gaus_peaks_to_json(vt->second,foreground,defaultPeakColor,override_alpha);

  json += "]";
  
  if( json.find("nan") != string::npos )
  {
    cerr << "Found nan: " << json << endl;
    cerr << endl;
  }
  
  return json;
}//string peak_json( inpeaks )
#endif //#if( SpecUtils_ENABLE_D3_CHART )

const std::shared_ptr<PeakContinuum> &PeakDef::continuum()
{
  return m_continuum;
}

std::shared_ptr<const PeakContinuum> PeakDef::continuum() const
{
  return m_continuum;
}

std::shared_ptr<PeakContinuum> PeakDef::getContinuum()
{
  return m_continuum;
}

void PeakDef::setContinuum( std::shared_ptr<PeakContinuum> continuum )
{
  if( !continuum )
    throw runtime_error( "PeakDef::setContinuum(...): invalid input" );
  m_continuum = continuum;
}//void setContinuum(...)


void PeakDef::makeUniqueNewContinuum()
{
  m_continuum = std::make_shared<PeakContinuum>(*m_continuum);
}


//The below should in principle take care of gaussian area and the skew area
double PeakDef::peakArea() const
{
  return m_coefficients[PeakDef::GaussAmplitude];
}//double peakArea() const


double PeakDef::peakAreaUncert() const
{
  return m_uncertainties[PeakDef::GaussAmplitude];
}//double peakAreaUncert() const


void PeakDef::setPeakArea( const double a )
{
  m_coefficients[PeakDef::GaussAmplitude] = a;
}//void PeakDef::setPeakArea( const double a )


void PeakDef::setPeakAreaUncert( const double uncert )
{
  m_uncertainties[PeakDef::GaussAmplitude] = uncert;
}//void setPeakAreaUncert( const double a )


double PeakDef::areaFromData( std::shared_ptr<const Measurement> data ) const
{
  double sumval = 0.0;
  if( !data || !m_continuum || !m_continuum->energyRangeDefined() )
    return sumval;
  
  try
  {
    const float energyStart = (float)m_continuum->lowerEnergy();
    const float energyEnd = (float)m_continuum->upperEnergy();
    
    const size_t lower_channel = data->find_gamma_channel( energyStart );
    const size_t upper_channel = data->find_gamma_channel(energyEnd);
    
    for( size_t i = lower_channel; i <= upper_channel; ++i )
    {
      const float e0 = std::max( energyStart, data->gamma_channel_lower(i) );
      const float e1 = std::min( energyEnd, data->gamma_channel_upper(i) );
      const double data_area_i = data->gamma_integral(e0, e1);
      // `this` is passed as the ROI's only peer - see the limitation noted at the declaration.
      const PeakDef *self = this;
      const double cont_area_1 = m_continuum->offset_integral( e0, e1, data, &self, 1 );
      if( data_area_i > cont_area_1 )
        sumval += (data_area_i - cont_area_1);
    }
  }catch( std::exception &e )
  {
    cerr << "Caught exception in PeakDef::areaFromData(): " << e.what() << endl;
    return 0.0;
  }
  
  return sumval;
}//double areaFromData( std::shared_ptr<const Measurement> data ) const;

bool PeakDef::lessThanByMean( const PeakDef &lhs, const PeakDef &rhs )
{
  return (lhs.m_coefficients[PeakDef::Mean] < rhs.m_coefficients[PeakDef::Mean]);
}//lessThanByMean(...)

bool PeakDef::lessThanByMeanShrdPtr( const std::shared_ptr<const PeakDef> &lhs,
                                  const std::shared_ptr<const PeakDef> &rhs )
{
  if( !lhs || !rhs )
    return (lhs < rhs);
  return lessThanByMean( *lhs, *rhs );
}


bool PeakDef::causilyDisconnected( const PeakDef &lower_peak,
                                    const PeakDef &upper_peak,
                                    const double ncausality,
                                   const bool useRoiAsWell )
{
  if( lower_peak.continuum() == upper_peak.continuum() )
    return false;
  
  double lower_upper( 0.0 ), upper_lower( 0.0 );

  if( lower_peak.mean() < upper_peak.mean() )
  {
    lower_upper = lower_peak.gausPeak() ? lower_peak.mean() + ncausality*lower_peak.sigma() : lower_peak.upperX();
    upper_lower = upper_peak.gausPeak() ? upper_peak.mean() - ncausality*upper_peak.sigma() : upper_peak.lowerX();
    if( useRoiAsWell )
    {
      lower_upper = std::max( lower_upper, lower_peak.upperX() );
      upper_lower = std::min( upper_lower, upper_peak.lowerX() );
    }
  }else
  {
    lower_upper = upper_peak.gausPeak() ? upper_peak.mean() + ncausality*upper_peak.sigma() : upper_peak.upperX();
    upper_lower = lower_peak.gausPeak() ? lower_peak.mean() - ncausality*lower_peak.sigma() : lower_peak.lowerX();
    if( useRoiAsWell )
    {
      lower_upper = std::max( lower_upper, upper_peak.upperX() );
      upper_lower = std::min( upper_lower, lower_peak.lowerX() );
    }
  }//if( lower_peak.mean() < upper_peak.mean() ) / else

  return (upper_lower > lower_upper);
}//bool PeakDef::causilyDisconnected


bool PeakDef::causilyConnected( const PeakDef &lower_peak,
                                   const PeakDef &upper_peak,
                                   const double ncausality,
                                   const bool useRoiAsWell )
{
  return !causilyDisconnected( lower_peak, upper_peak, ncausality, useRoiAsWell );
}


bool PeakDef::operator==( const PeakDef &rhs ) const
{
  for( CoefficientType t = CoefficientType(0);
       t < NumCoefficientTypes; t = CoefficientType(t+1) )
  {
    if( m_coefficients[t] != rhs.m_coefficients[t] )
      return false;
  }

  return m_type==rhs.m_type
      && (*m_continuum == *rhs.m_continuum)
      && (m_source == rhs.m_source)
      && m_useForEnergyCal==rhs.m_useForEnergyCal
      && m_useForShieldingSourceFit==rhs.m_useForShieldingSourceFit
      && m_useForManualRelEff==rhs.m_useForManualRelEff
      && m_useForDrfIntrinsicEffFit == rhs.m_useForDrfIntrinsicEffFit
      && m_useForDrfFwhmFit == rhs.m_useForDrfFwhmFit
      && m_useForDrfDepthOfInteractionFit == rhs.m_useForDrfDepthOfInteractionFit
      && m_lineColor==rhs.m_lineColor
      ;
}//PeakDef::operator==


bool PeakDef::Source::operator==( const Source &rhs ) const
{
  return (kind == rhs.kind) && (name == rhs.name) && (decay_parent == rhs.decay_parent)
         && (decay_child == rhs.decay_child) && (particle_energy == rhs.particle_energy)
         && (gamma_type == rhs.gamma_type);
}


void PeakDef::clearSources()
{
  m_source = Source();
}


bool PeakDef::hasSourceGammaAssigned() const
{
  return (m_source.kind != Source::Kind::None);
}


std::string PeakDef::sourceName() const
{
  return m_source.name;
}


const PeakDef::Source &PeakDef::source() const
{
  return m_source;
}


void PeakDef::setSource( const Source &src )
{
  m_source = (src.kind == Source::Kind::None) ? Source() : src;
}


PeakDef::SourceGammaType PeakDef::sourceGammaType() const
{
  return m_source.gamma_type;
}


float PeakDef::gammaParticleEnergy() const
{
  if( m_source.kind == Source::Kind::None )
    throw runtime_error( "Peak doesnt have a gamma associated with it" );

  switch( m_source.gamma_type )
  {
    case NormalGamma:
    case XrayGamma:
      break;
    case AnnihilationGamma:
      return 510.99891f;
    case SingleEscapeGamma:
      return m_source.particle_energy - 510.9989f;
    case DoubleEscapeGamma:
      return m_source.particle_energy - 2.0f*510.9989f;
  }//switch( m_source.gamma_type )

  return m_source.particle_energy;
}//float gammaParticleEnergy() const


void PeakDef::inheritUserSelectedOptions( const PeakDef &parent,
                                          const bool inheritNonFitForValues )
{
  m_source = parent.m_source;
  
  m_lineColor = parent.lineColor();
  
  m_userLabel = parent.m_userLabel;
  m_useForEnergyCal = parent.m_useForEnergyCal;
  m_useForShieldingSourceFit = parent.m_useForShieldingSourceFit;
  m_useForManualRelEff = parent.m_useForManualRelEff;
  m_useForDrfIntrinsicEffFit = parent.m_useForDrfIntrinsicEffFit;
  m_useForDrfFwhmFit = parent.m_useForDrfFwhmFit;
  m_useForDrfDepthOfInteractionFit = parent.m_useForDrfDepthOfInteractionFit;
  
  
  for( CoefficientType t = CoefficientType(0);
      t < NumCoefficientTypes; t = CoefficientType(t+1) )
  {
    m_fitFor[t] = parent.m_fitFor[t];
    
    if( inheritNonFitForValues && !m_fitFor[t] )
    {
      switch( t )
      {
        case PeakDef::Mean:             case PeakDef::Sigma:
        case PeakDef::GaussAmplitude:
        case PeakDef::SkewPar0: case PeakDef::SkewPar1:
        case PeakDef::SkewPar2: case PeakDef::SkewPar3:
        case PeakDef::SkewPar4: case PeakDef::SkewPar5:
          m_coefficients[t]  = parent.m_coefficients[t];
          m_uncertainties[t] = parent.m_uncertainties[t];
        break;
        
        case PeakDef::Chi2DOF: case PeakDef::NumCoefficientTypes:
          break;
      }//switch( t )
    }//if( inheritNonFitForValues )
  }//for( loop over PeakDef::CoefficientType )
  
  
  const auto &rhs_cont = parent.m_continuum;
  if( m_continuum && rhs_cont )
  {
    //Currently not copying the following values of continuum:
    //double m_lowerEnergy, m_upperEnergy;
    //double m_referenceEnergy;
    //std::vector<double> m_values, m_uncertainties;
    
    const vector<bool> rhs_fit_fors = rhs_cont->fitForParameter();
    const vector<bool> orig_fit_fors = m_continuum->fitForParameter();
    for( size_t i = 0; (i < rhs_fit_fors.size()) && (i < orig_fit_fors.size()); ++i )
    {
      m_continuum->setPolynomialCoefFitFor( i, rhs_fit_fors[i] );
    }
    
    // If we wanted to copy values of coefficients, we could use:
    /*
    if( (m_continuum->type() == rhs_cont->type())
       && (m_continuum->type() != PeakContinuum::OffsetType::External ) )
    {
      switch( m_continuum->type() )
      {
        case PeakContinuum::NoOffset:
          break;
        case PeakContinuum::External:
          //if( inheritNonFitForValues )
          //  m_continuum->setExternalContinuum( rhs_cont->externalContinuum() );
          break;
          
        case PeakContinuum::Constant:
        case PeakContinuum::Linear:
        case PeakContinuum::Quadratic:
        case PeakContinuum::Cubic:
        case PeakContinuum::FlatStep:
        case PeakContinuum::LinearStep:
        case PeakContinuum::BiLinearStep:
        {
          const vector<bool> rhs_fit_fors = rhs_cont->fitForParameter();
          assert( rhs_fit_fors.size() == m_continuum->fitForParameter().size() );
          
          for( size_t i = 0; i < rhs_fit_fors.size(); ++i )
          {
            m_continuum->setPolynomialCoefFitFor( i, rhs_fit_fors[i] );
          
            if( inheritNonFitForValues && !rhs_fit_fors[i] )
            {
              const vector<double> &rhs_pars = rhs_cont->parameters();
              const vector<double> &rhs_uncerts = rhs_cont->uncertainties();
              assert( rhs_pars.size() == rhs_fit_fors.size() );
              assert( rhs_pars.size() == rhs_uncerts.size() );
              
              m_continuum->setPolynomialCoef( i, rhs_pars[i] );
              m_continuum->setPolynomialUncert( i, rhs_uncerts[i] );
            }
          }//for( size_t i = 0; i < rhs_fit_fors.size(); ++i )
          
          break;
        }
      }//switch( m_continuum->type() )
    }
     */
  }//if( m_continuum )
  
}//void inheritUserSelectedOptions(...)


bool PeakContinuum::parametersProbablySet() const
{
  switch( m_type )
  {
    case NoOffset:
      return true;
    
    case Constant:     case Linear:
    case Quadratic:    case Cubic:
    case FlatStep:     case LinearStep:
    case BiLinearStep:
    case FlatStepCDF:  case LinearStepCDF:
    case BiLinearStepCDF:
    {
      for( const auto &v : m_values )
      {
        if( v != 0.0 )
          return true;
      }
      break;
    }// polynomial or step continuum

    case External:
      return !!m_externalContinuum;
    break;
  }//switch( m_type )
  
  return false;
}//bool defined() const



double PeakDef::lowerX() const
{
  if( m_continuum->PeakContinuum::energyRangeDefined() )
    return m_continuum->lowerEnergy();

  const double mean = m_coefficients[PeakDef::Mean];
  const double sigma = m_coefficients[PeakDef::Sigma];
  
  switch( m_skewType )
  {
    case NumSkewType:
      assert( 0 );
      //throw runtime_error( "PeakDef::lowerX: NumSkewType is not a valid skew type" );
      // Fall through to NoSkew for non-debug builds

    case NoSkew:
      return mean - 4.0*sigma;
      
    case Bortel:
    {
      // Find `x` so that only 2.5% of the distribution is to the left, but dont go any more
      //  than 8 sigma to the left of the mean
      // TODO: turn this next loop into a definite equation
      const double skew = m_coefficients[PeakDef::SkewPar0];
      double lowx = mean - 8.0*sigma;
      while( lowx < (mean - 4.0*sigma) )
      {
        const double cumulative = PeakDists::bortel_indefinite_integral(lowx, mean, sigma, skew);
        if( cumulative > 0.025 )
          return lowx;
        lowx += 0.1*sigma;
      }
      return mean - 4.0*sigma;
    }//case Bortel
      
    case CrystalBall:
    {
      // Find `x` so that only 2.5% of the distribution is to the left, but dont go any more
      //  than 8 sigma to the left of the mean
      // TODO: turn this next loop into a definite equation
      const double alpha = m_coefficients[PeakDef::SkewPar0];
      const double n = m_coefficients[PeakDef::SkewPar1];
      double lowx = mean - 8.0*sigma;
      while( (lowx < (mean - 4.0*sigma)) && (((lowx - mean)/sigma) < -alpha) )
      {
        const double cumulative = PeakDists::crystal_ball_tail_indefinite_t( sigma, alpha, n, (lowx - mean)/sigma );
        if( cumulative > 0.025 )
          return lowx;
        lowx += 0.1*sigma;
      }
      return mean - 4.0*sigma;
    }//case CrystalBall
      
    case DoubleSidedCrystalBall:
    {
      // Find `x` so that only 2.5% of the distribution is to the left, but dont go any more
      //  than 8 sigma to the left of the mean
      // TODO: turn this next loop into a definite equation, or search for 0.05, or something
      const double alpha_left = m_coefficients[PeakDef::SkewPar0];
      const double n_left = m_coefficients[PeakDef::SkewPar1];
      const double alpha_right = m_coefficients[PeakDef::SkewPar2];
      const double n_right = m_coefficients[PeakDef::SkewPar3];
      const double norm = PeakDists::DSCB_norm( alpha_left, n_left, alpha_right, n_right );
      
      double lowx = mean - 8.0*sigma;
      while( (lowx < (mean - 4.0*sigma)) && (((lowx - mean)/sigma) < -alpha_left) )
      {
        const double t = (lowx - mean) / sigma;
        const double cumulative = norm * PeakDists::DSCB_left_tail_indefinite_non_norm_t( alpha_left, n_left, t );
        
        if( cumulative > 0.025 )
          return lowx;
        lowx += 0.1*sigma;
      }
      return mean - 4.0*sigma;
    }//case CrystalBall
      
    case GaussExp:
    {
      // Find `x` so that only 2.5% of the distribution is to the left, but dont go any more
      //  than 8 sigma to the left of the mean
      // TODO: turn this next loop into a definite equation, or search for 0.05, or something
      const double skew = m_coefficients[PeakDef::SkewPar0];
      
      double lowx = mean - 8.0*sigma;
      while( (lowx < (mean - 4.0*sigma)) && (((lowx - mean)/sigma) < -skew) )
      {
        const double cumulative = PeakDists::gauss_exp_indefinite(mean, sigma, skew, lowx);
        if( cumulative > 0.025 )
          return lowx;
        lowx += 0.1*sigma;
      }
      return mean - 4.0*sigma;
    }//case GaussExp:
      
      
    case ExpGaussExp:
    {
      // Find `x` so that only 2.5% of the distribution is to the left, but dont go any more
      //  than 8 sigma to the left of the mean
      // TODO: turn this next loop into a definite equation, or search for 0.05, or something
      const double skew_left = m_coefficients[PeakDef::SkewPar0];
      const double skew_right = m_coefficients[PeakDef::SkewPar1];
      
      double lowx = mean - 8.0*sigma;
      while( (lowx < (mean - 4.0*sigma)) && (((lowx - mean)/sigma) < -skew_left) )
      {
        const double cumulative = PeakDists::exp_gauss_exp_indefinite( mean, sigma, skew_left, skew_right, lowx);
        if( cumulative > 0.025 )
          return lowx;
        lowx += 0.1*sigma;
      }
      return mean - 4.0*sigma;
    }//case ExpGaussExp:

    case VoigtPlusBortel:
    {
      // Use coverage_limits function to find lower bound
      const double gamma_lor = m_coefficients[PeakDef::SkewPar0];
      const double R = m_coefficients[PeakDef::SkewPar1];
      const double tau = m_coefficients[PeakDef::SkewPar2];

      try
      {
        const auto limits = PeakDists::voigt_exp_coverage_limits( mean, sigma, gamma_lor,
                                                                  R, tau, 0.05 );
        return limits.first;
      }
      catch( std::exception & )
      {
        // Fallback if coverage_limits fails
        return mean - (8.0*sigma + R*tau*sigma*10.0);
      }
    }//case VoigtPlusBortel:

    case GaussPlusBortel:
    {
      const double R = m_coefficients[PeakDef::SkewPar0];
      const double tau = m_coefficients[PeakDef::SkewPar1];

      try
      {
        const auto limits = PeakDists::gauss_plus_bortel_coverage_limits( mean, sigma, R, tau, 0.05 );
        return limits.first;
      }
      catch( std::exception & )
      {
        return mean - (8.0*sigma + R*tau*sigma*10.0);
      }
    }//case GaussPlusBortel:

    case DoubleBortel:
    {
      const double tau1 = m_coefficients[PeakDef::SkewPar0];
      const double tau2_delta = m_coefficients[PeakDef::SkewPar1];
      const double eta = m_coefficients[PeakDef::SkewPar2];

      try
      {
        const auto limits = PeakDists::double_bortel_coverage_limits( mean, sigma, tau1, tau2_delta, eta, 0.05 );
        return limits.first;
      }
      catch( std::exception & )
      {
        const double tau2 = tau1 + tau2_delta;
        return mean - (8.0*sigma + std::max(tau1, tau2)*sigma*10.0);
      }
    }//case DoubleBortel:

    case GadrasGeneric:
    case GadrasCZT:
    {
      const PeakDists::GadrasMaterial mat = (m_skewType == GadrasCZT)
                            ? PeakDists::GadrasMaterial::CZT_CdTe
                            : PeakDists::GadrasMaterial::Generic;
      double skew[6];
      for( int i = 0; i < 6; ++i )
        skew[i] = m_coefficients[PeakDef::SkewPar0 + i];

      try
      {
        const auto limits = PeakDists::gadras_coverage_limits( mean, sigma, skew, mat, 0.05 );
        return limits.first;
      }catch( std::exception & )
      {
        return mean - 8.0*sigma;
      }
    }//case GadrasGeneric, GadrasCZT:
  }//switch( m_skewType )

  assert( 0 );
  return mean - 4.0*sigma;
}//double lowerX() const


double PeakDef::upperX() const
{
  if( m_continuum->PeakContinuum::energyRangeDefined() )
    return m_continuum->upperEnergy();
  
  const double mean = m_coefficients[PeakDef::Mean];
  const double sigma = m_coefficients[PeakDef::Sigma];
  
  switch( m_skewType )
  {
    case NumSkewType:
      assert( 0 );
      //throw runtime_error( "PeakDef::upperX: NumSkewType is not a valid skew type" );
      // Fall through to NoSkew for non-debug builds

    case NoSkew:
    case GaussExp:
    case Bortel:
    case CrystalBall:
      return mean + 4.0*sigma;
      break;
      
    case DoubleSidedCrystalBall:
    {
      // Find `x` so that only 2.5% of the distribution is to the left, but dont go any more
      //  than 8 sigma to the left of the mean
      // TODO: turn this next loop into a definite equation, or search for 0.05, or something
      const double alpha_left = m_coefficients[PeakDef::SkewPar0];
      const double n_left = m_coefficients[PeakDef::SkewPar1];
      const double alpha_right = m_coefficients[PeakDef::SkewPar2];
      const double n_right = m_coefficients[PeakDef::SkewPar3];
      
      const double oneOverSqrt2 = LightMath::one_div_root_two; //0.70710678118654752440
      const double sqrtPiOver2 = LightMath::root_half_pi;      //1.2533141373155002512078826424
      
      
      double start_cumulative = 0.0;
      start_cumulative += PeakDists::DSCB_left_tail_indefinite_non_norm_t( alpha_left,
                                                      n_left, -alpha_left );
      start_cumulative += PeakDists::DSCB_gauss_indefinite_non_norm_t( alpha_right );
      start_cumulative -= PeakDists::DSCB_gauss_indefinite_non_norm_t( -alpha_left );
      start_cumulative -= PeakDists::DSCB_right_tail_indefinite_non_norm_t( alpha_right, n_right, alpha_right );

      const double norm = PeakDists::DSCB_norm( alpha_left, n_left, alpha_right, n_right );
      start_cumulative *= norm;
      
      double upper_t = alpha_right + 1;
      while( upper_t < 8 )
      {
        const double cumulative = start_cumulative
                + norm*PeakDists::DSCB_right_tail_indefinite_non_norm_t( alpha_right, n_right, upper_t );
        if( cumulative >= 0.975 )
          return (upper_t*sigma + mean);
        upper_t += 0.1;
      }
      return mean + 8.0*sigma;
    }//case DoubleSidedCrystalBall:
      
    case ExpGaussExp:
    {
      // Find `x` so that only 2.5% of the distribution is to the right, but dont go any more
      //  than 8 sigma to the left of the mean
      // TODO: turn this next loop into a definite equation, or search for 0.05, or something
      const double skew_left = m_coefficients[PeakDef::SkewPar0];
      const double skew_right = m_coefficients[PeakDef::SkewPar1];
      
      double upperx = mean + skew_right*sigma;
      while( upperx < (mean + 8.0*sigma) )
      {
        const double cumulative = PeakDists::exp_gauss_exp_indefinite( mean, sigma, skew_left, skew_right, upperx);
        if( cumulative >= 0.975 )
          return upperx;
        upperx += 0.1*sigma;
      }
      return mean + 8.0*sigma;
      break;
    }//case DoubleSidedCrystalBall:, case ExpGaussExp:

    case VoigtPlusBortel:
    {
      // Use coverage_limits function to find upper bound
      const double gamma_lor = m_coefficients[PeakDef::SkewPar0];
      const double R = m_coefficients[PeakDef::SkewPar1];
      const double tau = m_coefficients[PeakDef::SkewPar2];

      try
      {
        const auto limits = PeakDists::voigt_exp_coverage_limits( mean, sigma, gamma_lor,
                                                                  R, tau, 0.05 );
        return limits.second;
      }
      catch( std::exception & )
      {
        // Fallback if coverage_limits fails - Voigt has heavy tails due to Lorentzian
        return mean + 8.0*sigma + 5.0*gamma_lor;
      }
    }//case VoigtPlusBortel:

    case GaussPlusBortel:
    {
      const double R = m_coefficients[PeakDef::SkewPar0];
      const double tau = m_coefficients[PeakDef::SkewPar1];

      try
      {
        const auto limits = PeakDists::gauss_plus_bortel_coverage_limits( mean, sigma, R, tau, 0.05 );
        return limits.second;
      }
      catch( std::exception & )
      {
        return mean + 8.0*sigma;
      }
    }//case GaussPlusBortel:

    case DoubleBortel:
    {
      const double tau1 = m_coefficients[PeakDef::SkewPar0];
      const double tau2_delta = m_coefficients[PeakDef::SkewPar1];
      const double eta = m_coefficients[PeakDef::SkewPar2];

      try
      {
        const auto limits = PeakDists::double_bortel_coverage_limits( mean, sigma, tau1, tau2_delta, eta, 0.05 );
        return limits.second;
      }
      catch( std::exception & )
      {
        return mean + 8.0*sigma;
      }
    }//case DoubleBortel:

    case GadrasGeneric:
    case GadrasCZT:
    {
      const PeakDists::GadrasMaterial mat = (m_skewType == GadrasCZT)
                            ? PeakDists::GadrasMaterial::CZT_CdTe
                            : PeakDists::GadrasMaterial::Generic;
      double skew[6];
      for( int i = 0; i < 6; ++i )
        skew[i] = m_coefficients[PeakDef::SkewPar0 + i];

      try
      {
        const auto limits = PeakDists::gadras_coverage_limits( mean, sigma, skew, mat, 0.05 );
        return limits.second;
      }catch( std::exception & )
      {
        return mean + 8.0*sigma;
      }
    }//case GadrasGeneric, GadrasCZT:
  }//switch( m_skewType )

  assert( 0 );
  return mean + 4.0*sigma;
}//double upperX() const


double PeakDef::gauss_integral( const double x0, const double x1 ) const
{
  const double &mean = m_coefficients[CoefficientType::Mean];
  const double &sigma = m_coefficients[CoefficientType::Sigma];
  const double &amp = m_coefficients[CoefficientType::GaussAmplitude];
  
  switch( m_skewType )
  {
    case NumSkewType:
      assert( 0 );
      //throw runtime_error( "PeakDef::gauss_integral: NumSkewType is not a valid skew type" );
      // Fall through to NoSkew for non-debug builds

    case SkewType::NoSkew:
      return amp*PeakDists::gaussian_integral( mean, sigma, x0, x1 );
          
    case SkewType::Bortel:
      return amp*PeakDists::bortel_integral( mean, sigma, m_coefficients[CoefficientType::SkewPar0], x0, x1 );
      
    //case SkewType::Doniach:
    //  return amp*doniach_integral( x0, x1, mean, sigma, m_coefficients[CoefficientType::SkewPar0] );
      
    case SkewType::CrystalBall:
      return amp*PeakDists::crystal_ball_integral( mean, sigma,
                            m_coefficients[CoefficientType::SkewPar0],
                            m_coefficients[CoefficientType::SkewPar1],
                            x0, x1 );
      
    case SkewType::DoubleSidedCrystalBall:
      return amp*PeakDists::double_sided_crystal_ball_integral( mean, sigma,
                                       m_coefficients[CoefficientType::SkewPar0],
                                       m_coefficients[CoefficientType::SkewPar1],
                                       m_coefficients[CoefficientType::SkewPar2],
                                       m_coefficients[CoefficientType::SkewPar3],
                                                         x0, x1 );
      break;
      
    case SkewType::GaussExp:
      return amp*PeakDists::gauss_exp_integral( mean, sigma,
                                         m_coefficients[CoefficientType::SkewPar0], x0, x1 );
      break;
      
    case SkewType::ExpGaussExp:
      return amp*PeakDists::exp_gauss_exp_integral( mean, sigma,
                                             m_coefficients[CoefficientType::SkewPar0],
                                             m_coefficients[CoefficientType::SkewPar1], x0, x1 );
      break;

    case SkewType::VoigtPlusBortel:
      return amp*PeakDists::voigt_exp_integral( mean, sigma,
                                                m_coefficients[CoefficientType::SkewPar0],
                                                m_coefficients[CoefficientType::SkewPar1],
                                                m_coefficients[CoefficientType::SkewPar2],
                                                x0, x1 );
      break;

    case SkewType::GaussPlusBortel:
      return amp*PeakDists::gauss_plus_bortel_integral( mean, sigma,
                                                        m_coefficients[CoefficientType::SkewPar0],
                                                        m_coefficients[CoefficientType::SkewPar1],
                                                        x0, x1 );
      break;

    case SkewType::DoubleBortel:
      return amp*PeakDists::double_bortel_integral( mean, sigma,
                                                    m_coefficients[CoefficientType::SkewPar0],
                                                    m_coefficients[CoefficientType::SkewPar1],
                                                    m_coefficients[CoefficientType::SkewPar2],
                                                    x0, x1 );
      break;

    case SkewType::GadrasGeneric:
    case SkewType::GadrasCZT:
    {
      const PeakDists::GadrasMaterial mat = (m_skewType == SkewType::GadrasCZT)
                            ? PeakDists::GadrasMaterial::CZT_CdTe
                            : PeakDists::GadrasMaterial::Generic;
      return amp*PeakDists::gadras_integral( mean, sigma,
                                             m_coefficients + CoefficientType::SkewPar0,
                                             mat, x0, x1 );
    }
  };//enum SkewType

  assert( 0 );
  throw runtime_error( "Invalid skew type" );
  return 0.0;
}//double gauss_integral( const double x0, const double x1 ) const;


void PeakDef::gauss_integral( const float *energies, double *channels, const size_t nchannel ) const 
{
  const double mean = m_coefficients[PeakDef::Mean];
  const double sigma = m_coefficients[PeakDef::Sigma];
  const double amp = m_coefficients[PeakDef::GaussAmplitude];
  
  PeakDists::photopeak_function_integral( mean, sigma, amp, m_skewType,
                                       m_coefficients + PeakDef::SkewPar0,
                                       nchannel, energies, channels );
}//void gauss_integral( const float * const energies, double *channels, const size_t nchannel )




////////////////////////////////////////////////////////////////////////////////

/** Text appropriate for use as a label for the continuum type in the gui. */
const char *PeakContinuum::offset_type_label_tr( const PeakContinuum::OffsetType type )
{
  switch( type )
  {
    case PeakContinuum::NoOffset:     return "pct-none";
    case PeakContinuum::Constant:     return "pct-constant";
    case PeakContinuum::Linear:       return "pct-linear";
    case PeakContinuum::Quadratic:    return "pct-quadratic";
    case PeakContinuum::Cubic:        return "pct-cubic";
    case PeakContinuum::FlatStep:      return "pct-flat-step";
    case PeakContinuum::LinearStep:    return "pct-linear-step";
    case PeakContinuum::BiLinearStep:  return "pct-bilinear-step";
    case PeakContinuum::FlatStepCDF:   return "pct-flat-step-cdf";
    case PeakContinuum::LinearStepCDF:    return "pct-linear-step-cdf";
    case PeakContinuum::BiLinearStepCDF: return "pct-bilinear-step-cdf";
    case PeakContinuum::External:        return "pct-global";
  }//switch( type )
  
  assert( 0 );
  return "pct-invalid";
}


const char *PeakContinuum::offset_type_str( const PeakContinuum::OffsetType type )
{
  switch( type )
  {
    case NoOffset:     return "NoOffset";
    case Constant:     return "Constant";
    case Linear:       return "Linear";
    case Quadratic:    return "Quardratic"; //Note mispelling of "Quardratic" is left for backwards compatibility (InterSpec v1.0.6 and before), but should eventually be fixed...
    case Cubic:        return "Cubic";
    case FlatStep:      return "FlatStep";
    case LinearStep:    return "LinearStep";
    case BiLinearStep:  return "BiLinearStep";
    case FlatStepCDF:      return "FlatStepCDF";
    case LinearStepCDF:    return "LinearStepCDF";
    case BiLinearStepCDF:  return "BiLinearStepCDF";
    case External:         return "External";
  }//switch( m_type )
  
  return "InvalidOffsetType";
}


size_t PeakContinuum::num_parameters( const PeakContinuum::OffsetType type )
{
  switch( type )
  {
    case OffsetType::NoOffset:
    case OffsetType::External:
      return 0;
      
    case OffsetType::Constant:
    case OffsetType::Linear:
    case OffsetType::Quadratic:
    case OffsetType::Cubic:
      return static_cast<size_t>(type);
      
    case OffsetType::FlatStep:
    case OffsetType::LinearStep:
    case OffsetType::BiLinearStep:
      return 2 + (type - FlatStep);

    case OffsetType::FlatStepCDF:
      return 2;
    case OffsetType::LinearStepCDF:
      return 3;
    case OffsetType::BiLinearStepCDF:
      return 4;
  }//switch( type )

  assert( 0 );
  throw std::runtime_error( "Somehow invalid continuum polynomial type." );

  return 0;
}//size_t num_parameters( const OffsetType type );


size_t PeakContinuum::num_linear_fit_pars( const PeakContinuum::OffsetType type )
{
  switch( type )
  {
    case OffsetType::NoOffset:
    case OffsetType::External:
      return 0;

    case OffsetType::Constant:
    case OffsetType::Linear:
    case OffsetType::Quadratic:
    case OffsetType::Cubic:
      return static_cast<size_t>( type );

    case OffsetType::FlatStep:
    case OffsetType::LinearStep:
    case OffsetType::BiLinearStep:
      return 2 + (type - FlatStep);

    case OffsetType::FlatStepCDF:
      return 1;  // Constant polynomial only; step coefficient is optimized by non-linear solver
    case OffsetType::LinearStepCDF:
      return 2;  // Linear polynomial; step coefficient is optimized by non-linear solver
    case OffsetType::BiLinearStepCDF:
      return 2;  // Linear polynomial; the two step coefficients are optimized by non-linear solver
  }//switch( type )

  assert( 0 );
  throw std::runtime_error( "Somehow invalid continuum type in num_linear_fit_pars." );

  return 0;
}//size_t num_linear_fit_pars( const OffsetType type );


size_t PeakContinuum::num_cdf_step_pars( const PeakContinuum::OffsetType type )
{
  switch( type )
  {
    case OffsetType::NoOffset:     case OffsetType::External:
    case OffsetType::Constant:     case OffsetType::Linear:
    case OffsetType::Quadratic:    case OffsetType::Cubic:
    case OffsetType::FlatStep:     case OffsetType::LinearStep:
    case OffsetType::BiLinearStep:
      return 0;

    case OffsetType::FlatStepCDF:
    case OffsetType::LinearStepCDF:
      return 1;

    case OffsetType::BiLinearStepCDF:
      return 2;
  }//switch( type )

  assert( 0 );
  throw std::runtime_error( "Somehow invalid continuum type in num_cdf_step_pars." );

  return 0;
}//size_t num_cdf_step_pars( const OffsetType type );


bool PeakContinuum::is_step_continuum( const OffsetType type )
{
  switch( type )
  {
    case PeakContinuum::NoOffset: case PeakContinuum::External:
    case PeakContinuum::Constant: case PeakContinuum::Linear:
    case PeakContinuum::Quadratic: case PeakContinuum::Cubic:
      return false;
      
    case PeakContinuum::FlatStep:
    case PeakContinuum::LinearStep:
    case PeakContinuum::BiLinearStep:
    case PeakContinuum::FlatStepCDF:
    case PeakContinuum::LinearStepCDF:
    case PeakContinuum::BiLinearStepCDF:
      return true;
  }//switch( cont->type() )

  assert( 0 );
  throw std::runtime_error( "Somehow invalid continuum polynomial type." );
  return false;
}//bool is_step_continuum( const OffsetType type );


bool PeakContinuum::is_peak_cdf_step_continuum( const OffsetType type )
{
  switch( type )
  {
    case PeakContinuum::NoOffset:     case PeakContinuum::External:
    case PeakContinuum::Constant:     case PeakContinuum::Linear:
    case PeakContinuum::Quadratic:    case PeakContinuum::Cubic:
    case PeakContinuum::FlatStep:     case PeakContinuum::LinearStep:
    case PeakContinuum::BiLinearStep:
      return false;

    case PeakContinuum::FlatStepCDF:
    case PeakContinuum::LinearStepCDF:
    case PeakContinuum::BiLinearStepCDF:
      return true;
  }//switch( type )

  assert( 0 );
  return false;
}//bool is_peak_cdf_step_continuum( const OffsetType type )


PeakContinuum::OffsetType PeakContinuum::str_to_offset_type_str( const char * const str, const size_t len )
{
  const auto compare = []( const char *a, const size_t alen, const char *b, const size_t blen, const bool ) -> bool {
    return SpecUtils::iequals_ascii( std::string(a, alen), std::string(b, blen) );
  };
  
  /*
   // Alternative implementation for this function that is untested
  for( OffsetType type = OffsetType(0); type <= External; type = OffsetType(type+1) )
  {
    const char * const teststr = offset_type_str(type);
    const size_t teststr_len = strlen(teststr);
    if( compare(str,len,teststr,teststr_len,false) )
      return type;
  }
  if( compare(str,len,"Quadratic",9,false) )
    return PeakContinuum::Quadratic;
  throw runtime_error( "Invalid continuum type" );
  */
  
  if( compare(str,len,"NoOffset",8,false) )
    return PeakContinuum::NoOffset;
  
  if( compare(str,len,"Constant",8,false) )
    return PeakContinuum::Constant;
  
  if( compare(str,len,"Linear",6,false) )
    return PeakContinuum::Linear;
  
  if( compare(str,len,"Quardratic",10,false) || compare(str,len,"Quadratic",9,false) )
    return PeakContinuum::Quadratic;
  
  if( compare(str,len,"Cubic",5,false) )
    return PeakContinuum::Cubic;
  
  if( compare(str,len,"FlatStep",8,false) )
    return PeakContinuum::FlatStep;
  
  if( compare(str,len,"LinearStep",10,false) )
    return PeakContinuum::LinearStep;
  
  if( compare(str,len,"BiLinearStepCDF",15,false) )
    return PeakContinuum::BiLinearStepCDF;

  if( compare(str,len,"BiLinearStep",12,false) )
    return PeakContinuum::BiLinearStep;

  if( compare(str,len,"FlatStepCDF",11,false) )
    return PeakContinuum::FlatStepCDF;

  if( compare(str,len,"LinearStepCDF",13,false) )
    return PeakContinuum::LinearStepCDF;

  if( compare(str,len,"External",8,false) )
    return PeakContinuum::External;
  
  throw runtime_error( "Invalid continuum type" );
}//str_to_offset_type_str(...)



PeakContinuum::PeakContinuum()
: m_type( PeakContinuum::NoOffset ),
  m_lowerEnergy( 0.0 ),
  m_upperEnergy( 0.0 ),
  m_referenceEnergy( 0.0 )
{
}//PeakContinuum constructor

bool PeakContinuum::operator==( const PeakContinuum &rhs ) const
{
  return m_type==rhs.m_type
         && m_lowerEnergy == rhs.m_lowerEnergy
         && m_lowerEnergy == rhs.m_lowerEnergy
         && m_upperEnergy == rhs.m_upperEnergy
         && m_referenceEnergy == rhs.m_referenceEnergy
         && m_values == rhs.m_values
         && m_uncertainties == rhs.m_uncertainties
         && m_fitForValue == rhs.m_fitForValue
         && m_externalContinuum == rhs.m_externalContinuum;
}

void PeakContinuum::setParameters( double referenceEnergy,
                                   const std::vector<double> &x,
                                   const std::vector<double> &uncertainties )
{
  // First check size of inputs are valid
  const size_t num_expected_pars = PeakContinuum::num_parameters( m_type );
  
  if( x.size() != num_expected_pars )
    throw runtime_error( "PeakContinuum::setParameters invalid parameter size" );
  
  if( !uncertainties.empty() && (uncertainties.size() != num_expected_pars) )
    throw runtime_error( "PeakContinuum::setParameters invalid uncert size" );
  
  m_values = x;
  m_referenceEnergy = referenceEnergy;
  m_fitForValue.resize( m_values.size(), true );
  m_uncertainties = uncertainties;
  m_uncertainties.resize( m_values.size(), 0.0 );
}//void setParameters(...)


void PeakContinuum::setParameters( double referenceEnergy,
                                   const double *parameters,
                                   const double *uncertainties )
{
  if( !parameters )
    throw runtime_error( "PeakContinuum::setParameters invalid parameters" );
  
  switch( m_type )
  {
    case NoOffset: case External:
      throw runtime_error( "PeakContinuum::setParameters(): called for external or no-offset continuum - not allowed" );
      
    case Constant:   case Linear:
    case Quadratic: case Cubic:
    case FlatStep:
    case LinearStep:
    case BiLinearStep:
    case FlatStepCDF:
    case LinearStepCDF:
    case BiLinearStepCDF:
    {
      const size_t npar = num_parameters(m_type);

      m_values.resize( npar );
      m_uncertainties.resize( npar );
      m_fitForValue.resize( npar, true );
      m_referenceEnergy = referenceEnergy;

      for( size_t i = 0; i < npar; ++i )
      {
        m_values[i] = parameters[i];
        m_uncertainties[i] = uncertainties ? uncertainties[i] : 0.0;
      }

      break;
    }//case - polynomial continuum
  };//switch( m_type )
}//setParameters


bool PeakContinuum::setPolynomialCoefFitFor( size_t polyCoefNum, bool fit )
{
  if( polyCoefNum >= m_fitForValue.size() )
    return false;
  
  m_fitForValue[polyCoefNum] = fit;
  return true;
}//bool setPolynomialCoefFitFor( size_t polyCoefNum, bool fit )


bool PeakContinuum::setPolynomialCoef( size_t polyCoef, double val )
{
  if( polyCoef >= m_values.size() )
    return false;
  
  m_values[polyCoef] = val;
  return true;
}


bool PeakContinuum::setPolynomialUncert( size_t polyCoef, double val )
{
  if( polyCoef >= m_uncertainties.size() )
    return false;
  
  m_uncertainties[polyCoef] = val;
  return true;
}


void PeakContinuum::setExternalContinuum( const std::shared_ptr<const Measurement> &data )
{
  if( m_type != External )
    throw runtime_error( "PeakContinuum::setExternalContinuum invalid m_type" );
 
  m_externalContinuum = data;
}//setExternalContinuum(...)


void PeakContinuum::setRange( const double lower, const double upper )
{
  m_lowerEnergy = lower;
  m_upperEnergy = upper;
  if( m_lowerEnergy > m_upperEnergy )
    std::swap( m_lowerEnergy, m_upperEnergy );
}//void setRange( const double lowerenergy, const double upperenergy )


bool PeakContinuum::energyRangeDefined() const
{
  return (m_lowerEnergy != m_upperEnergy);
}//bool energyRangeDefined() const


bool PeakContinuum::isPolynomial() const
{
  switch( m_type )
  {
    case NoOffset:   case External:
      return false;
      
    case Constant:      case Linear:
    case Quadratic:     case Cubic:
    case FlatStep:      case LinearStep:
    case BiLinearStep:
    case FlatStepCDF:   case LinearStepCDF:
    case BiLinearStepCDF:
      return true;

    break;
  }//switch( m_type )

  return false;
}//bool isPolynomial() const


namespace
{
/** What a single `PeakContinuum` parameter slot *means*.

 `PeakContinuum::setType(...)` uses this to decide which stored values may be carried across a
 type change: a value may only be copied into a slot with the identical meaning.  This matters
 most for the two kinds of "step", which are numerically incompatible:
   - `DataStep` (FlatStep/LinearStep) multiplies the in-ROI data fraction, so it is a count
     density in counts/keV - typically O(1) to O(100).
   - `CdfStep0` (FlatStepCDF/LinearStepCDF) multiplies the ROI peak area, so it is in keV^-1 -
     typically O(1e-3) to O(1e-6).
 The two differ by the total ROI peak area, i.e. a factor of 1e3 to 1e6.  Carrying one into the
 other collapses the fitted peak amplitude to ~0 and inflates the continuum.
 */
enum class ContParMeaning
{
  Poly0, Poly1, Poly2, Poly3,  //!< Polynomial terms, relative to the reference energy
  RightPoly0, RightPoly1,      //!< Right-hand line of a bi-linear step
  DataStep,                    //!< Step magnitude of FlatStep/LinearStep; counts/keV
  CdfStep0,                    //!< Step coefficient of the CDF step types; keV^-1
  CdfStep1                     //!< Energy-slope of BiLinearStepCDF's step coefficient; keV^-2
};//enum class ContParMeaning


/** Returns the meaning of each parameter slot of `type`.

 The returned size always equals `PeakContinuum::num_parameters(type)`; this table is the
 authoritative statement of what every continuum type's parameters are, and hence of which
 `setType(...)` transitions may preserve which values.
 */
const std::vector<ContParMeaning> &continuum_parameter_meanings( const PeakContinuum::OffsetType type )
{
  using M = ContParMeaning;

  static const std::vector<M> s_empty{};
  static const std::vector<M> s_constant{ M::Poly0 };
  static const std::vector<M> s_linear{ M::Poly0, M::Poly1 };
  static const std::vector<M> s_quadratic{ M::Poly0, M::Poly1, M::Poly2 };
  static const std::vector<M> s_cubic{ M::Poly0, M::Poly1, M::Poly2, M::Poly3 };
  static const std::vector<M> s_flat_step{ M::Poly0, M::DataStep };
  static const std::vector<M> s_linear_step{ M::Poly0, M::Poly1, M::DataStep };
  static const std::vector<M> s_bilinear_step{ M::Poly0, M::Poly1, M::RightPoly0, M::RightPoly1 };
  static const std::vector<M> s_flat_step_cdf{ M::Poly0, M::CdfStep0 };
  static const std::vector<M> s_linear_step_cdf{ M::Poly0, M::Poly1, M::CdfStep0 };
  static const std::vector<M> s_bilinear_step_cdf{ M::Poly0, M::Poly1, M::CdfStep0, M::CdfStep1 };

  switch( type )
  {
    case PeakContinuum::NoOffset:        return s_empty;
    case PeakContinuum::Constant:        return s_constant;
    case PeakContinuum::Linear:          return s_linear;
    case PeakContinuum::Quadratic:       return s_quadratic;
    case PeakContinuum::Cubic:           return s_cubic;
    case PeakContinuum::FlatStep:        return s_flat_step;
    case PeakContinuum::LinearStep:      return s_linear_step;
    case PeakContinuum::BiLinearStep:    return s_bilinear_step;
    case PeakContinuum::FlatStepCDF:     return s_flat_step_cdf;
    case PeakContinuum::LinearStepCDF:   return s_linear_step_cdf;
    case PeakContinuum::BiLinearStepCDF: return s_bilinear_step_cdf;
    case PeakContinuum::External:        return s_empty;
  }//switch( type )

  assert( 0 );
  throw std::runtime_error( "Somehow invalid continuum type in continuum_parameter_meanings." );
}//continuum_parameter_meanings( OffsetType )


/** Index of `meaning` within `meanings`, or `meanings.size()` if not present. */
size_t index_of_meaning( const std::vector<ContParMeaning> &meanings, const ContParMeaning meaning )
{
  for( size_t i = 0; i < meanings.size(); ++i )
  {
    if( meanings[i] == meaning )
      return i;
  }

  return meanings.size();
}//index_of_meaning(...)
}//namespace


void PeakContinuum::setType( PeakContinuum::OffsetType type )
{
  const PeakContinuum::OffsetType oldType = m_type;

  const std::vector<ContParMeaning> &old_meanings = continuum_parameter_meanings( oldType );
  const std::vector<ContParMeaning> &new_meanings = continuum_parameter_meanings( type );

#if( PERFORM_DEVELOPER_CHECKS )
  // The meaning table and num_parameters(...) are independent switches that must stay in step; if
  //  they drift, every transition through here silently mis-sizes the continuum.
  assert( old_meanings.size() == num_parameters(oldType) );
  assert( new_meanings.size() == num_parameters(type) );
#endif

  // Only trust as many old slots as the old type actually declared *and* we actually hold.
  const size_t num_old = std::min( old_meanings.size(), m_values.size() );

  std::vector<double> new_values( new_meanings.size(), 0.0 );
  std::vector<double> new_uncerts( new_meanings.size(), 0.0 );
  std::vector<bool> new_fit_for( new_meanings.size(), true );

  // Carry over every value whose meaning is unchanged; anything else starts at zero.
  for( size_t i = 0; i < new_meanings.size(); ++i )
  {
    const size_t j = index_of_meaning( old_meanings, new_meanings[i] );
    if( j >= num_old )
      continue;

    new_values[i] = m_values[j];
    if( j < m_uncertainties.size() )
      new_uncerts[i] = m_uncertainties[j];
    if( j < m_fitForValue.size() )
      new_fit_for[i] = m_fitForValue[j];
  }//for( size_t i = 0; i < new_meanings.size(); ++i )

  // Coming from a type with no right-hand line, seed it from the left one; starting a bi-linear
  //  step from right==left avoids an artificial discontinuity part way through the ROI.
  const size_t right0 = index_of_meaning( new_meanings, ContParMeaning::RightPoly0 );
  if( (right0 < new_meanings.size())
     && (index_of_meaning( old_meanings, ContParMeaning::RightPoly0 ) >= num_old) )
  {
    const size_t right1 = index_of_meaning( new_meanings, ContParMeaning::RightPoly1 );
    const size_t left0 = index_of_meaning( new_meanings, ContParMeaning::Poly0 );
    const size_t left1 = index_of_meaning( new_meanings, ContParMeaning::Poly1 );

    assert( (right1 < new_meanings.size()) && (left0 < new_meanings.size())
           && (left1 < new_meanings.size()) );

    new_values[right0] = new_values[left0];
    new_values[right1] = new_values[left1];
    new_uncerts[right0] = new_uncerts[left0];
    new_uncerts[right1] = new_uncerts[left1];
    // ...including whether they are being fit; a copy of a pinned line must not be free to move.
    new_fit_for[right0] = new_fit_for[left0];
    new_fit_for[right1] = new_fit_for[left1];
  }//if( the new type has a right-hand line the old type did not )

  m_type = type;
  m_values.swap( new_values );
  m_uncertainties.swap( new_uncerts );
  m_fitForValue.swap( new_fit_for );

  if( type == External )
    m_referenceEnergy = 0.0;
  else
    m_externalContinuum.reset();
}//void setType( PeakContinuum::OffsetType type )


void PeakContinuum::calc_linear_continuum_eqn( const std::shared_ptr<const SpecUtils::Measurement> &data,
                                              const double reference_energy,
                                              const double roi_start, const double roi_end,
                               const size_t num_lower_channels,
                               const size_t num_upper_channels )
{
  assert( data );
  assert( reference_energy >= roi_start );
  assert( reference_energy <= roi_end );
  assert( roi_end >= roi_start );
  assert( num_lower_channels > 0 );
  assert( num_upper_channels > 0 );
  
  if( !data || !data->energy_calibration() || !data->energy_calibration()->valid() )
    throw runtime_error( "PeakContinuum::calc_linear_continuum_eqn: invalid input data" );
  
  if( (reference_energy < roi_start) || (reference_energy > roi_end) )
    throw runtime_error( "PeakContinuum::calc_linear_continuum_eqn: reference energy must be within ROI" );
  
  if( roi_end < roi_start )
    throw runtime_error( "PeakContinuum::calc_linear_continuum_eqn: lower energy greater than upper energy" );
  
  if( !num_lower_channels || !num_upper_channels )
    throw runtime_error( "PeakContinuum::calc_linear_continuum_eqn: number of above/below channels must not be zero" );
  
  m_referenceEnergy = reference_energy;
  m_lowerEnergy = roi_start;
  m_upperEnergy = roi_end;
  m_type = PeakContinuum::Linear;
  m_values.resize( 2 );
  m_fitForValue.resize( 2, true );
  
  const shared_ptr<const SpecUtils::EnergyCalibration> cal = data->energy_calibration();
  assert( cal && cal->valid() );
  
  // We'll try to make the side channels be completely independent from the channels in the ROI, but
  //  we'll let there be up 10% of a channel overlap.
  double lower_cont_bound = cal->channel_for_energy( roi_start );
  double upper_cont_bound = cal->channel_for_energy( roi_end );
  
  double frac_part = lower_cont_bound - std::floor(lower_cont_bound);
  if( frac_part > 0.9 )
    lower_cont_bound = std::round( lower_cont_bound );
  else
    lower_cont_bound = std::floor( lower_cont_bound );
  
  frac_part = upper_cont_bound - std::floor(upper_cont_bound);
  if( frac_part < 0.1 )
    upper_cont_bound = std::round( upper_cont_bound );
  else
    upper_cont_bound = std::floor( upper_cont_bound );
  upper_cont_bound = std::max( upper_cont_bound, 1.0 );
  
  const size_t roi_first_channel = static_cast<size_t>( lower_cont_bound );
  // Upper cont bound give the channel number whos lower edge defines the upper bound on ROI, so we
  //  will to subtract 1 from it to find the last channel _in_ the ROI, however, we dont want to
  //  make the upper ROI channel be lower than the low one
  const size_t roi_last_channel = std::max( roi_first_channel, static_cast<size_t>( upper_cont_bound - 1.0 ) );
  
  double &m = m_values[1];
  double &b = m_values[0];
  
  PeakContinuum::eqn_from_offsets( roi_first_channel, roi_last_channel, reference_energy, data,
                                   num_lower_channels, num_upper_channels, m, b );
}//PeakContinuum::calc_linear_continuum_eqn


double PeakContinuum::offset_integral_non_cdf( const double x0, const double x1,
                                               const std::shared_ptr<const SpecUtils::Measurement> &data ) const
{
  
  // A lambda to integrate data for step function continuums
  auto integrate_for_step = [this,x0,x1,data]( double &roi_lower, double &roi_upper, double &cumulative_data, double &roi_data_sum ) {
    
    if( !data || !data->num_gamma_channels() )
      throw runtime_error( "PeakContinuum::offset_integral: invalid data spectrum passed in" );
    
    roi_lower = lowerEnergy();
    roi_upper = upperEnergy();
    
    // To be consistent with how fit_amp_and_offset(...) handles things, we will do our own
    //  summing here, rather than calling Measurement::gamma_integral.
    const vector<float> &counts = *data->gamma_counts();
    const size_t mid_channel = data->find_gamma_channel( 0.5*(x0 + x1) );
    const size_t lower_channel = data->find_gamma_channel( roi_lower );
    const size_t upper_channel = data->find_gamma_channel( roi_upper );
    
    // Adjusting roi_lower and roi_upper for rounding to channel edges to match fit_amp_and_offset(...)
    roi_lower = data->gamma_channel_lower(lower_channel);
    roi_upper = data->gamma_channel_upper(upper_channel);
    
    //assert( mid_channel >= lower_channel ); //PeakDef::gaus_peaks_to_json will call a few channels below and above ROI extents
    //assert( mid_channel <= upper_channel );
    assert( lower_channel <= upper_channel );
    assert( upper_channel < counts.size() );
    
#if( PEAK_CONTINUUM_DATA_STEP_SUBTRACT )
    // Find the minimum channel count in the ROI to subtract before forming the CDF;
    //  this concentrates the CDF shape around the peak rather than the flat continuum.
    float min_count = counts[lower_channel];
    for( size_t i = lower_channel + 1; i <= upper_channel; ++i )
      min_count = std::min( min_count, counts[i] );
    min_count = std::max( min_count, 0.0f );
#else
    const float min_count = 0.0f;
#endif

    cumulative_data = 0.0;
    roi_data_sum = 0.0;
    for( size_t i = lower_channel; i <= upper_channel; ++i )
    {
      roi_data_sum += (counts[i] - min_count);
      cumulative_data += (i < mid_channel) ? static_cast<double>( counts[i] - min_count ) : 0.0;
    }

    if( (mid_channel >= lower_channel) && (mid_channel <= upper_channel) )
      cumulative_data += 0.5 * (counts[mid_channel] - min_count);
  };//integrate_for_step lamda
  
  
  switch( m_type )
  {
    case NoOffset:
      return 0.0;
    
    case Constant: case Linear: case Quadratic: case Cubic:
      return offset_eqn_integral( &(m_values[0]), m_type, x0, x1, m_referenceEnergy );
    
    case FlatStep:
    case LinearStep:
    {
      assert( m_values.size() == (2 + (m_type - FlatStep)) );
      
      if( !data || !data->num_gamma_channels() )
        throw runtime_error( "PeakContinuum::offset_integral: invalid data spectrum passed in" );
      
      // If you change any of this, make sure you update fit_amp_and_offset(...), as well as
      //  the as well `PeakContinuum::offset_integral( const float *, double *, size_t, data )`.

      double roi_lower, roi_upper, cumulative_data, roi_data_sum;
      integrate_for_step( roi_lower, roi_upper, cumulative_data, roi_data_sum );
      
      const double x0_rel = x0 - m_referenceEnergy;
      const double x1_rel = x1 - m_referenceEnergy;
      
      const double frac_data = (roi_data_sum > 0.0) ? (cumulative_data / roi_data_sum) : 0.5;

      const double offset = m_values[0]*(x1_rel - x0_rel);
      const double linear = ((m_type == FlatStep) ? 0.0 :  0.5*m_values[1]*(x1_rel*x1_rel - x0_rel*x0_rel));
      const size_t step_index = ((m_type == FlatStep) ? 1 : 2);
      const double step_contribution = m_values[step_index] * frac_data * (x1_rel - x0_rel);
       
      const double answer = offset + linear + step_contribution;
      
      return std::max( answer, 0.0 );
    }//case FlatStep and LinearStep

      
    case BiLinearStep:
    {
      assert( m_values.size() == 4 );
      
      double roi_lower, roi_upper, cumulative_data, roi_data_sum;
      integrate_for_step( roi_lower, roi_upper, cumulative_data, roi_data_sum );

      const double x0_rel = x0 - m_referenceEnergy;
      const double x1_rel = x1 - m_referenceEnergy;

      const double frac_data = (roi_data_sum > 0.0) ? (cumulative_data / roi_data_sum) : 0.5;

      const double left_poly = m_values[0]*(x1_rel - x0_rel) + 0.5*m_values[1]*(x1_rel*x1_rel - x0_rel*x0_rel);
      const double right_poly = m_values[2]*(x1_rel - x0_rel) + 0.5*m_values[3]*(x1_rel*x1_rel - x0_rel*x0_rel);
      
      const double contrib = ((1.0 - frac_data) * left_poly) + (frac_data * right_poly);
      
      return std::max( 0.0, contrib );
    }//case BiLinearStep:
      
    case FlatStepCDF:
    case LinearStepCDF:
    case BiLinearStepCDF:
      assert( 0 && "offset_integral_non_cdf called for CDF step type" );
      throw runtime_error( "PeakContinuum::offset_integral_non_cdf: CDF step types require"
                          " offset_integral_cdf_step instead" );

    case External:
      if( !m_externalContinuum )
        return 0.0;
      return gamma_integral( m_externalContinuum, x0, x1 );
    break;
  };//enum OffsetType

  return 0.0;
}//double offset_integral_non_cdf( x0, x1, data )


double PeakContinuum::offset_integral( const double x0, const double x1,
                                       const std::shared_ptr<const SpecUtils::Measurement> &data,
                                       const PeakDef * const *roi_peaks,
                                       const size_t num_peaks ) const
{
  if( is_peak_cdf_step_continuum( m_type ) )
    return offset_integral_cdf_step( x0, x1, data, roi_peaks, num_peaks );
  return offset_integral_non_cdf( x0, x1, data );
}//double offset_integral( x0, x1, data, roi_peaks, num_peaks )


void PeakContinuum::offset_integral( const float *energies, double *channels, const size_t nchannel,
                                     const std::shared_ptr<const SpecUtils::Measurement> &data,
                                     const PeakDef * const *roi_peaks,
                                     const size_t num_peaks ) const
{
  if( is_peak_cdf_step_continuum( m_type ) )
    offset_integral_cdf_step( energies, channels, nchannel, data, roi_peaks, num_peaks );
  else
    offset_integral_non_cdf( energies, channels, nchannel, data );
}//void offset_integral( energies, channels, nchannel, data, roi_peaks, num_peaks )


double PeakContinuum::offset_integral( const double x0, const double x1,
                                       const std::shared_ptr<const SpecUtils::Measurement> &data,
                                       const std::vector<std::shared_ptr<const PeakDef>> &roi_peaks ) const
{
  // Delegate to the raw-pointer overload
  std::vector<const PeakDef *> ptrs( roi_peaks.size() );
  for( size_t i = 0; i < roi_peaks.size(); ++i )
    ptrs[i] = roi_peaks[i].get();
  return offset_integral( x0, x1, data, ptrs.data(), ptrs.size() );
}


void PeakContinuum::offset_integral( const float *energies, double *channels, const size_t nchannel,
                                     const std::shared_ptr<const SpecUtils::Measurement> &data,
                                     const std::vector<std::shared_ptr<const PeakDef>> &roi_peaks ) const
{
  std::vector<const PeakDef *> ptrs( roi_peaks.size() );
  for( size_t i = 0; i < roi_peaks.size(); ++i )
    ptrs[i] = roi_peaks[i].get();
  offset_integral( energies, channels, nchannel, data, ptrs.data(), ptrs.size() );
}


namespace
{
/** Peak j's ROI-anchored CDF for the channel spanning [x0,x1].

 The fitters build their step basis by cumulative-summing unit-area channel integrals across the
 ROI (`PeakDists::unit_pdf_to_cdf`), which for channel i gives
   sum_{k<i}(pdf_k) + 0.5*pdf_i  ==  0.5*(CDF(x0) + CDF(x1)) - CDF(roi_lower)
 i.e. the *average of the channel's edge CDFs*, not the CDF at its center.  Evaluating the analytic
 CDF the same way makes the two agree to floating-point rather than only to O(channel width^2) -
 using the midpoint instead leaves a trapezoid-vs-midpoint difference of about
 `(dx^2/8)*pdf'(center)`, which for a 20k-count peak on 0.35 keV channels is a few hundredths of a
 count per channel.

 Clamping to [f0,f1] anchors the result at the ROI's lower edge and saturates it at the upper edge,
 matching where the fitters' cumulative sum starts and stops; see
 `PeakContinuum::cdf_step_anchor_energies(...)`.
 */
double anchored_peak_cdf( const double x0, const double x1, const PeakDef &peak,
                          const double f0, const double f1 )
{
  const double * const skew_pars = peak.coefficients() + PeakDef::CoefficientType::SkewPar0;
  const double cdf0 = PeakDists::peak_cdf( x0, peak.mean(), peak.sigma(), peak.skewType(), skew_pars );
  const double cdf1 = PeakDists::peak_cdf( x1, peak.mean(), peak.sigma(), peak.skewType(), skew_pars );
  const double cdf = 0.5*(cdf0 + cdf1);

  return (std::min)( (std::max)( cdf, f0 ), f1 ) - f0;
}//anchored_peak_cdf(...)
}//namespace


bool PeakContinuum::cdf_step_anchor_energies( const std::shared_ptr<const SpecUtils::Measurement> &data,
                                              double &lower, double &upper ) const
{
  lower = 0.0;
  upper = 0.0;

  // With no ROI defined there is nothing to anchor to; caller then uses the un-anchored CDF.
  if( !energyRangeDefined() )
    return false;

  lower = m_lowerEnergy;
  upper = m_upperEnergy;

  // Prefer the channel edges the fitters actually start/stop their cumulative sums at, so the
  //  drawn continuum matches the fitted one.  `data` is optional - e.g. LeafletRadMap deliberately
  //  passes nullptr for multi-Measurement selections - so fall back to the raw ROI bounds, which
  //  differ by at most a fraction of a channel of CDF mass at the ROI edge.
  //  This mirrors what `offset_integral_non_cdf(...)` already does for the data-step types.
  //
  // `find_gamma_channel(...)` floors, which matches how every ROI in InterSpec is turned into a
  //  channel range - except `RelActCalcAuto`s `RoiRangeChannels::channel_range()`, which rounds to
  //  the nearest channel and then stores the caller's unrounded energies here.  For a RelAct ROI
  //  the anchors below can therefore be one channel away from the ones the fit used; see the note
  //  on `channel_range()` for the magnitude (<0.1% of the continuum) and for why the obvious fix
  //  is wrong.
  if( data && data->num_gamma_channels() )
  {
    const std::shared_ptr<const SpecUtils::EnergyCalibration> cal = data->energy_calibration();
    if( cal && cal->valid() )
    {
      const size_t lower_channel = data->find_gamma_channel( static_cast<float>(m_lowerEnergy) );
      const size_t upper_channel = data->find_gamma_channel( static_cast<float>(m_upperEnergy) );
      lower = data->gamma_channel_lower( lower_channel );
      upper = data->gamma_channel_upper( upper_channel );
    }
  }//if( data && data->num_gamma_channels() )

  return true;
}//bool cdf_step_anchor_energies( data, lower, upper ) const


double PeakContinuum::offset_integral_cdf_step( const double x0, const double x1,
                                                 const std::shared_ptr<const SpecUtils::Measurement> &data,
                                                 const PeakDef * const *roi_peaks, const size_t num_peaks ) const
{
  assert( is_peak_cdf_step_continuum( m_type ) );
  if( !is_peak_cdf_step_continuum( m_type ) )
    throw runtime_error( "PeakContinuum::offset_integral_cdf_step: only valid for CDF step types" );

  const size_t num_poly = num_linear_fit_pars( m_type );
  const size_t num_step = num_cdf_step_pars( m_type );
  assert( m_values.size() == (num_poly + num_step) );

  const double x0_rel = x0 - m_referenceEnergy;
  const double x1_rel = x1 - m_referenceEnergy;
  const double dx = x1_rel - x0_rel;
  // Average the *relative* edges rather than subtracting the reference from the absolute centre:
  //  the absolute energies are ~1e3 while E' is ~1e0, so the latter throws away several digits.
  const double center_rel = 0.5 * (x0_rel + x1_rel);

  // Polynomial part, integrated exactly across the channel
  double poly = m_values[0] * dx;
  if( num_poly > 1 )
    poly += 0.5 * m_values[1] * (x1_rel * x1_rel - x0_rel * x0_rel);

  double anchor_lower = 0.0, anchor_upper = 0.0;
  const bool anchored = cdf_step_anchor_energies( data, anchor_lower, anchor_upper );

  // Step part: (SUM_k s_k * E'^k) * SUM_j( amp_j * CDFbar_j(center) ) * dx
  double cdf_weighted_sum = 0.0;
  for( size_t j = 0; j < num_peaks; ++j )
  {
    const PeakDef *peak = roi_peaks[j];
    if( !peak )
      continue;

    const double f0 = anchored ? PeakDists::peak_cdf( anchor_lower, peak->mean(), peak->sigma(),
                                    peak->skewType(),
                                    peak->coefficients() + PeakDef::CoefficientType::SkewPar0 ) : 0.0;
    const double f1 = anchored ? PeakDists::peak_cdf( anchor_upper, peak->mean(), peak->sigma(),
                                    peak->skewType(),
                                    peak->coefficients() + PeakDef::CoefficientType::SkewPar0 ) : 1.0;

    cdf_weighted_sum += peak->amplitude() * anchored_peak_cdf( x0, x1, *peak, f0, f1 );
  }//for( size_t j = 0; j < num_peaks; ++j )

  double step_coeff = m_values[num_poly];
  if( num_step > 1 )
    step_coeff += m_values[num_poly + 1] * center_rel;

  return (std::max)( 0.0, poly + (step_coeff * cdf_weighted_sum * dx) );
}//double offset_integral_cdf_step( x0, x1, data, roi_peaks, num_peaks )


void PeakContinuum::offset_integral_cdf_step( const float *energies, double *channels, const size_t nchannel,
                                               const std::shared_ptr<const SpecUtils::Measurement> &data,
                                               const PeakDef * const *roi_peaks, const size_t num_peaks ) const
{
  assert( is_peak_cdf_step_continuum( m_type ) );
  if( !is_peak_cdf_step_continuum( m_type ) )
    throw runtime_error( "PeakContinuum::offset_integral_cdf_step: only valid for CDF step types" );

  assert( nchannel > 0 );
  if( !nchannel )
    return;

  const size_t num_poly = num_linear_fit_pars( m_type );
  const size_t num_step = num_cdf_step_pars( m_type );
  assert( m_values.size() == (num_poly + num_step) );

  double anchor_lower = 0.0, anchor_upper = 0.0;
  const bool anchored = cdf_step_anchor_energies( data, anchor_lower, anchor_upper );

  // Precompute per-peak info; the anchors in particular must be hoisted out of the channel loop -
  //  this batch overload exists because stepped continua are otherwise ~5x slower.
  struct PeakInfo
  {
    const PeakDef *peak;
    double amp;
    double f0, f1;                 //!< CDF at the ROI's lower/upper anchor energy
    std::vector<double> edge_cdf;  //!< CDF at each channel edge; `nchannel + 1` entries
  };

  std::vector<PeakInfo> peak_infos;
  peak_infos.reserve( num_peaks );
  for( size_t j = 0; j < num_peaks; ++j )
  {
    const PeakDef *peak = roi_peaks[j];
    if( !peak )
      continue;

    const double * const skew_pars = peak->coefficients() + PeakDef::CoefficientType::SkewPar0;

    PeakInfo info;
    info.peak = peak;
    info.amp = peak->amplitude();
    info.f0 = anchored ? PeakDists::peak_cdf( anchor_lower, peak->mean(), peak->sigma(),
                                              peak->skewType(), skew_pars ) : 0.0;
    info.f1 = anchored ? PeakDists::peak_cdf( anchor_upper, peak->mean(), peak->sigma(),
                                              peak->skewType(), skew_pars ) : 1.0;

    // Adjacent channels share an edge, so evaluate the (erf-heavy) CDF once per edge rather than
    //  twice per channel - this batch overload exists precisely to keep stepped continua fast.
    info.edge_cdf.resize( nchannel + 1 );
    for( size_t i = 0; i <= nchannel; ++i )
      info.edge_cdf[i] = PeakDists::peak_cdf( energies[i], peak->mean(), peak->sigma(),
                                              peak->skewType(), skew_pars );

    peak_infos.push_back( std::move(info) );
  }//for( size_t j = 0; j < num_peaks; ++j )

  for( size_t i = 0; i < nchannel; ++i )
  {
    const double x0_rel = energies[i] - m_referenceEnergy;
    const double x1_rel = energies[i + 1] - m_referenceEnergy;
    const double dx = x1_rel - x0_rel;
    // NB: `energies` is float, so `0.5*(energies[i] + energies[i+1])` would sum in float at ~1e3
    //  and lose ~6e-5 keV of E'; average the relative (double) edges instead.  Only the peak CDF,
    //  which genuinely needs absolute energies, uses the channel edges directly.
    const double center_rel = 0.5 * (x0_rel + x1_rel);

    double poly = m_values[0] * dx;
    if( num_poly > 1 )
      poly += 0.5 * m_values[1] * (x1_rel * x1_rel - x0_rel * x0_rel);

    double cdf_weighted_sum = 0.0;
    for( const PeakInfo &info : peak_infos )
    {
      const double cdf = 0.5*(info.edge_cdf[i] + info.edge_cdf[i+1]);
      cdf_weighted_sum += info.amp * ((std::min)( (std::max)( cdf, info.f0 ), info.f1 ) - info.f0);
    }

    double step_coeff = m_values[num_poly];
    if( num_step > 1 )
      step_coeff += m_values[num_poly + 1] * center_rel;

    channels[i] += (std::max)( 0.0, poly + (step_coeff * cdf_weighted_sum * dx) );
  }//for( size_t i = 0; i < nchannel; ++i )
}//void offset_integral_cdf_step( energies, channels, nchannel, data, roi_peaks, num_peaks )


void PeakContinuum::offset_integral_non_cdf( const float *energies, double *channels, const size_t nchannel,
                                             const std::shared_ptr<const SpecUtils::Measurement> &data ) const
{
  // This function should give the same answer as
  //  `PeakContinuum::offset_integral( double x0, const double x1, data)`, on a channel-by-channel
  //  basis, just be a little faster computationally, especially for stepped continua.
  PeakDists::offset_integral<PeakContinuum,double>( *this, energies, channels, nchannel, data );

  /*
  assert( nchannel > 0 );
  if( !nchannel )
    return;
  
  switch( m_type )
  {
    case NoOffset:
      return;
      
    case Constant: case Linear: case Quadratic: case Cubic:
    {
      for( size_t i = 0; i < nchannel; ++i )
      {
        const double x0 = energies[i] - m_referenceEnergy;
        const double x1 = energies[i+1] - m_referenceEnergy;
        
        double answer = 0.0;
        switch( m_type )
        {
          case NoOffset: case External: case FlatStep: case LinearStep: case BiLinearStep:
            assert( 0 );
            
          case Cubic:
            assert( m_values.size() == 4 );
            answer += 0.25*m_values[3]*(x1*x1*x1*x1 - x0*x0*x0*x0);
            //fall-through intentional
            
          case Quadratic:
            assert( m_values.size() >= 3 );
            answer += 0.333333333333333*m_values[2]*(x1*x1*x1 - x0*x0*x0);
            //fall-through intentional
            
          case Linear:
            assert( m_values.size() >= 2 );
            answer += 0.5*m_values[1]*(x1*x1 - x0*x0);
            //fall-through intentional
            
          case Constant:
            assert( m_values.size() >= 1 );
            answer += m_values[0]*(x1 - x0);
            break;
        };//switch( type )
        
        assert( std::max(answer,0.0) == offset_eqn_integral( &(m_values[0]), m_type, energies[i], energies[i+1], m_referenceEnergy ) );
        
        channels[i] += std::max( answer, 0.0 );
      }//for( size_t i = 0; i < nchannel; ++i )
      
      break;
    }//case Constant: case Linear: case Quadratic: case Cubic:
      
      
    case FlatStep:
    case LinearStep:
    case BiLinearStep:
    {
      if( !data || !data->num_gamma_channels() )
        throw runtime_error( "PeakContinuum::offset_integral: invalid data spectrum passed in" );
      
      // To be consistent with how fit_amp_and_offset(...) handles things, we will do our own
      //  summing here, rather than calling Measurement::gamma_integral.
      const size_t roi_lower_channel = data->find_gamma_channel( m_lowerEnergy );
      const size_t roi_upper_channel = data->find_gamma_channel( m_upperEnergy );
      
      //const double roi_lower = data->gamma_channel_lower(roi_lower_channel);
      //const double roi_upper = data->gamma_channel_upper(roi_upper_channel);
      
      const vector<float> &counts = *data->gamma_counts();
      
      assert( roi_lower_channel < counts.size() );
      assert( roi_upper_channel < counts.size() );
      
      
#if( PEAK_CONTINUUM_DATA_STEP_SUBTRACT )
      float min_count = counts[roi_lower_channel];
      for( size_t i = roi_lower_channel + 1; i <= roi_upper_channel; ++i )
        min_count = std::min( min_count, counts[i] );
      min_count = std::max( min_count, 0.0f );
#else
      const float min_count = 0.0f;
#endif
      // Compute roi_data_sum by per-element subtraction of min_count to match the
      //  single-channel offset_integral_non_cdf computation (avoid floating-point
      //  discrepancy from subtracting min_count*N at the end).
      double roi_data_sum = 0.0;
      for( size_t i = roi_lower_channel; i <= roi_upper_channel; ++i )
        roi_data_sum += (counts[i] - min_count);

      const size_t begin_channel = data->find_gamma_channel( energies[0] );
      assert( energies[0] == data->gamma_channel_lower(begin_channel) );
      const size_t end_channel = begin_channel + nchannel; //one past last channel we want

      // Lets check that `energies` points into the lower channel energies of `data`.
      //  We actually only care that the values of the array are the same, but we'll be
      //  a little tighter for development.
      const vector<float> &data_energies = *data->channel_energies();

      if( (energies[0] != data->gamma_channel_lower(begin_channel))
         || (energies[nchannel] != data->gamma_channel_lower(begin_channel+nchannel)) )
        throw std::logic_error( "PeakContinuum::offset_integral: for stepped continua" );

      double cumulative_data = 0.0;

      // In case we are starting part-way into the ROI
      if( begin_channel > roi_lower_channel )
      {
        for( size_t i = roi_lower_channel; i < begin_channel; ++i )
          cumulative_data += (counts[i] - min_count);
      }//if( begin_channel > roi_lower_channel )
      
      
      for( size_t i = begin_channel; i < end_channel; ++i )
      {
        const size_t input_index = i - begin_channel;
        assert( data_energies[i] == energies[input_index] );
        
        const double x0_rel = data_energies[i] - m_referenceEnergy;
        const double x1_rel = data_energies[i+1] - m_referenceEnergy;
        
        if( i >= roi_lower_channel && i <= roi_upper_channel )
          cumulative_data += 0.5 * (counts[i] - min_count);

        const double frac_data = (roi_data_sum > 0.0) ? (cumulative_data / roi_data_sum) : 0.5;

        switch( m_type )
        {
          case FlatStep:
          case LinearStep:
          {
            const double offset = m_values[0]*(x1_rel - x0_rel);
            const double linear = ((m_type == FlatStep) ? 0.0 :  0.5*m_values[1]*(x1_rel*x1_rel - x0_rel*x0_rel));
            const size_t step_index = ((m_type == FlatStep) ? 1 : 2);
            const double step_contribution = m_values[step_index] * frac_data * (x1_rel - x0_rel);

            const double answer = std::max( 0.0, offset + linear + step_contribution );

            assert( answer == offset_integral( data_energies[i], data_energies[i+1], data ) );

            channels[input_index] += answer;
            break;
          }//case FlatStep: case LinearStep:

          case BiLinearStep:
          {
            assert( m_values.size() == 4 );
            const double left_poly = m_values[0]*(x1_rel - x0_rel) + 0.5*m_values[1]*(x1_rel*x1_rel - x0_rel*x0_rel);
            const double right_poly = m_values[2]*(x1_rel - x0_rel) + 0.5*m_values[3]*(x1_rel*x1_rel - x0_rel*x0_rel);
            const double contrib = std::max( 0.0, ((1.0 - frac_data) * left_poly) + (frac_data * right_poly) );

            assert( contrib == offset_integral( data_energies[i], data_energies[i+1], data ) );

            channels[input_index] += contrib;
            break;
          }//case BiLinearStep:

          case NoOffset: case Constant: case Linear: case Quadratic: case Cubic: case External:
            assert( 0 );
            break;
        }//switch( m_type )

        if( i >= roi_lower_channel && i <= roi_upper_channel )
          cumulative_data += 0.5 * (counts[i] - min_count);
      }//for( size_t i = 0; i < channels; ++i )
      
      break;
    }//case FlatStep/LinearStep/BiLinearStep
      
    case External:
    {
      assert( m_externalContinuum );
      if( !m_externalContinuum )
        break;
      
      for( size_t i = 0; i < nchannel; ++i )
        channels[i] += gamma_integral( m_externalContinuum, energies[i], energies[i+1] );
      break;
    }//case External:
  }//switch( m_type )
  */
}//void PeakContinuum::offset_integral_non_cdf( energies, channels, nchannel, data )


void PeakContinuum::eqn_from_offsets( size_t lowchannel,
                             size_t highchannel,
                             const double reference_energy,
                             const std::shared_ptr<const SpecUtils::Measurement> &data,
                             const size_t num_lower_channels,
                             const size_t num_upper_channels,
                             double &m, double &b )
{
  // y = m*x + b, where y is density of continuum counts, per keV, at energy x
  
  
  if( !data || !data->energy_calibration() || !data->energy_calibration()->valid() )
    throw runtime_error( "PeakContinuum::calc_linear_continuum_eqn: invalid input data" );
  
  if( lowchannel > highchannel )
    throw runtime_error( "PeakContinuum::eqn_from_offsets: lower channel greater than upper channel" );
  
  if( !num_lower_channels || !num_upper_channels )
    throw runtime_error( "PeakContinuum::eqn_from_offsets: number of above/below channels must not be zero" );
  
  const shared_ptr<const SpecUtils::EnergyCalibration> cal = data->energy_calibration();
  assert( cal && cal->valid() );
  const size_t nchannel = cal->num_channels();
  const size_t last_channel = nchannel ? nchannel - 1 : size_t(0);
  
  if( (num_upper_channels >= nchannel) || (num_lower_channels >= nchannel) )
    throw runtime_error( "PeakContinuum::eqn_from_offsets: number of above/below channels is to large" );
  
  lowchannel = std::max( lowchannel, num_lower_channels );
  if( (highchannel + num_upper_channels) >= nchannel )
    highchannel = nchannel - num_upper_channels - 1;
    
  
  // Now get the first and last channels for regions above/below ROI - note that these will be
  //  inclusive; e.g., if num_lower_channels==1, then first and last channel will be equal.
  //  We will also make sure to not go below zero, but if we go above last channel in spectrum,
  //  we'll just be lazy and assume they have zero values anyway.
  
  const size_t lower_cont_first_channel = lowchannel - num_lower_channels;
  const size_t lower_cont_last_channel =  lowchannel - 1;
  
  const size_t upper_cont_first_channel = highchannel + 1;
  const size_t upper_cont_last_channel = upper_cont_first_channel + num_upper_channels - 1;
  
  const double lower_count_sum = data->gamma_channels_sum( lower_cont_first_channel, lower_cont_last_channel );
  const double upper_count_sum = data->gamma_channels_sum( upper_cont_first_channel, upper_cont_last_channel );
  
  const double lower_low_energy = cal->energy_for_channel( std::min(lower_cont_first_channel, last_channel) );
  const double lower_up_energy = cal->energy_for_channel( std::min(lower_cont_last_channel + 1, last_channel) );
  
  const double upper_low_energy = cal->energy_for_channel( std::min(upper_cont_first_channel, last_channel) );
  const double upper_up_energy = cal->energy_for_channel( std::min(upper_cont_last_channel + 1, last_channel) );
  
  const double lower_dx = lower_up_energy - lower_low_energy;
  const double upper_dx = upper_up_energy - upper_low_energy;
  
  const double y1 = lower_count_sum / (1.0 + lower_cont_last_channel - lower_cont_first_channel);
  const double y2 = upper_count_sum / (1.0 + upper_cont_last_channel - upper_cont_first_channel);
  
  const double lower_x1_rel   = lower_low_energy - reference_energy;
  const double lower_x2_rel   = lower_up_energy - reference_energy;
  const double upper_x1_rel   = upper_low_energy - reference_energy;
  const double upper_x2_rel   = upper_up_energy - reference_energy;
  const double lower_sqr_diff = lower_x2_rel*lower_x2_rel - lower_x1_rel*lower_x1_rel;
  const double upper_sqr_diff = upper_x2_rel*upper_x2_rel - upper_x1_rel*upper_x1_rel;
  
  // The equation "mx + b" gives us the density of continuum counts at energy = x, so we will
  //  integrate this equation over the region below the ROI, and above the ROI, and do a little
  //  algebra so that we solve for m and b, such that we get the exact amount of counts from the
  //  equation, as is in data
  //
  //lower_count_sum = lower_dx * b + 0.5*m*lower_sqr_diff
  //upper_count_sum = upper_dx * b + 0.5*m*upper_sqr_diff
  //  ==> b = (lower_count_sum - 0.5*m*lower_sqr_diff)/lower_dx
  //  ==> upper_count_sum = upper_dx * ((lower_count_sum - 0.5*m*lower_sqr_diff)/lower_dx) + 0.5*m*upper_sqr_diff
  //      upper_count_sum = upper_dx*lower_count_sum/lower_dx - upper_dx*0.5*m*lower_sqr_diff/lower_dx + 0.5*m*upper_sqr_diff
  //      upper_count_sum = upper_dx*lower_count_sum/lower_dx + 0.5*m*(upper_sqr_diff - upper_dx*lower_sqr_diff/lower_dx)
  //      upper_count_sum - upper_dx*lower_count_sum/lower_dx = 0.5*m*(upper_sqr_diff - upper_dx*lower_sqr_diff/lower_dx)
  //  ==> m = 2*(upper_count_sum - upper_dx*lower_count_sum/lower_dx) / (upper_sqr_diff - upper_dx*lower_sqr_diff/lower_dx)
  
  if( fabs( (upper_sqr_diff - (upper_dx*lower_sqr_diff/lower_dx)) ) < FLT_EPSILON )
  {
    // Extraordinarily unlikely to happen; just assign the continuum to be flat, and average density
    //  of below and above ROI
    m = 0.0;
    b = (0.5 * lower_count_sum / lower_dx) + (0.5 * upper_count_sum / upper_dx);
  }else
  {
    const double ratio_dx = upper_dx / lower_dx;
    const double mult_m = 2.0 / (upper_sqr_diff - ratio_dx*lower_sqr_diff);
    
    // Calc m first, since b is dependent on m.
    m = mult_m * (upper_count_sum - (upper_dx*lower_count_sum/lower_dx));
    b = (lower_count_sum - 0.5*m*lower_sqr_diff) / lower_dx;
    
    // TODO: also give uncertainty of parameters, dont think its to hard, e.g. (but needs to be
    //       double checked, I did the differentiation/summing correctly in my head):
    //  const double uncert_m = mult_m * sqrt( upper_count_sum + ratio_dx*ratio_dx*lower_count_sum );
    //  const double uncert_b = (1.0 / lower_dx) * sqrt( lower_count_sum + 0.5*0.5*lower_sqr_diff*lower_sqr_diff*uncert_m*uncert_m )
  }
  
#if( PERFORM_DEVELOPER_CHECKS )
  const double coefs[2] = { b, m };
  const double eqn_lower_counts = PeakContinuum::offset_eqn_integral( coefs,
                                                                     PeakContinuum::OffsetType::Linear,
                                                                     lower_low_energy, lower_up_energy,
                                                                     reference_energy );
  
  const double eqn_upper_counts = PeakContinuum::offset_eqn_integral( coefs,
                                                                     PeakContinuum::OffsetType::Linear,
                                                                     upper_low_energy, upper_up_energy,
                                                                     reference_energy );
  
  if( fabs(eqn_lower_counts - lower_count_sum) > 1.0E-4*std::max(eqn_lower_counts,lower_count_sum)
     && (lower_count_sum > 1.0E-3)  )
  {
    char buffer[512] = { '\0' };
    snprintf( buffer, sizeof(buffer), "Region below ROI data integral not equal to equation integral."
             "\n\tReference energy: %f"
             "\n\tlower_low_energy: %f"
             "\n\tlower_up_energy: %f"
             "\n\teqn_lower_counts: %f, data lower_count_sum: %f"
             "\n\tcoefficients: %f, %f\n",
             reference_energy, lower_low_energy, lower_up_energy,
             eqn_lower_counts, lower_count_sum, b, m );
    log_developer_error( __func__, buffer );
    assert( 0 );
  }
  
  if( (fabs(eqn_upper_counts - upper_count_sum) > 1.0E-4*std::max(eqn_upper_counts,upper_count_sum))
     && (upper_count_sum > 1.0E-3) )
  {
    char buffer[512] = { '\0' };
    snprintf( buffer, sizeof(buffer), "Region Above ROI data integral not equal to equation integral."
             "\n\tReference energy: %f"
             "\n\tupper_low_energy: %f"
             "\n\tupper_up_energy: %f"
             "\n\teqn_upper_counts: %f, data upper_count_sum: %f"
             "\n\tcoefficients: %f, %f\n",
             reference_energy, upper_low_energy, upper_up_energy,
             eqn_upper_counts, upper_count_sum, b, m );
    log_developer_error( __func__, buffer );
    assert( 0 );
  }
#endif
  
  if( IsNan(m) || IsInf(m) || IsNan(b) || IsInf(b) )
  {
    // Shouldnt ever get here.
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( __func__, "PeakContinuum::eqn_from_offsets(...): Invalid results" );
#else
    cerr << "PeakContinuum::eqn_from_offsets(...): Invalid results" << endl;
#endif
    m = b = 0.0;
  }//if( an invalid value of m or b )
  
  //cout << "For reference energy " << reference_energy << " fit ceofs {" << b << ", " << m << "}" << endl;
}//void PeakContinuum::eqn_from_offsets(...)


double PeakContinuum::offset_eqn_integral( const double *coefs,
                                          PeakContinuum::OffsetType type,
                                          double x0, double x1,
                                          const double peak_mean )
{
  double answer = 0.0;
  x0 -= peak_mean;
  x1 -= peak_mean;
  
  //Explicitly evaluating the polynomial speeds up peak fitting by about a factor of two - suprising!
  //const int maxorder = static_cast<int>( type );
  //for( int order = 0; order < maxorder; ++order )
  //{
  //  const double exp = order + 1.0;
  //  answer += (coefs[order]/exp) * (std::pow(x1,exp) - std::pow(x0,exp));
  //}//for( int order = 0; order < maxorder; ++order )
  
  switch( type )
  {
    case NoOffset: case External:
    case FlatStep: case LinearStep: case BiLinearStep:
    case FlatStepCDF: case LinearStepCDF: case BiLinearStepCDF:
      throw runtime_error( "PeakContinuum::offset_eqn_integral(...) may only be"
                          " called for polynomial backgrounds" );
      
    case Cubic:
      answer += 0.25*coefs[3]*(x1*x1*x1*x1 - x0*x0*x0*x0);
      //fallthrough intentional
      
    case Quadratic:
      answer += 0.333333333333333*coefs[2]*(x1*x1*x1 - x0*x0*x0);
      //fallthrough intentional
      
    case Linear:
      answer += 0.5*coefs[1]*(x1*x1 - x0*x0);
      //fallthrough intentional
      
    case Constant:
      answer += coefs[0]*(x1 - x0);
      break;
  };//enum OffsetType
  
  return std::max( answer, 0.0 );
}//offset_eqn_integral(...)




void PeakContinuum::translate_offset_polynomial( double *new_coefs,
                                                 const double *old_coefs,
                                                 PeakContinuum::OffsetType type,
                                                 const double new_center,
                                                 const double old_center )
{
  switch( type )
  {
    case NoOffset:
      return;
    
    case Constant:
      new_coefs[0] = old_coefs[0];
      break;
      
    case Linear:
      new_coefs[0] = old_coefs[0] + old_coefs[1] * (new_center - old_center);
      new_coefs[1] = old_coefs[1];
      break;
      
    case Quadratic:
    case Cubic:
      throw runtime_error( "translate_offset_polynomial does not yet support "
                           "quadratic or cubic polynomials" );
    
    case LinearStep:
    case LinearStepCDF:
      new_coefs[0] = old_coefs[0] + old_coefs[1] * (new_center - old_center);
      new_coefs[1] = old_coefs[1];
      new_coefs[2] = old_coefs[2];
      break;
      
    case BiLinearStep:
    case BiLinearStepCDF:
      new_coefs[0] = old_coefs[0] + old_coefs[1] * (new_center - old_center);
      new_coefs[1] = old_coefs[1];
      new_coefs[2] = old_coefs[2] + old_coefs[3] * (new_center - old_center);
      new_coefs[3] = old_coefs[3];
      break;
    
    case FlatStep:
    case FlatStepCDF:
      new_coefs[0] = old_coefs[0];
      new_coefs[1] = old_coefs[1];
      break;      

    case External:
      throw runtime_error( "translate_offset_polynomial does not support external continuum" );
  }//switch( type )
}//void translate_offset_polynomial(...)

