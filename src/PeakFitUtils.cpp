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

#include <vector>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;

namespace PeakFitUtils
{
  
float nai_fwhm_fcn( const float energy )
{
  static const vector<float> nai_fwhm_coefs{ -4.0f, 6.3f, 0.6f };   //"NaI 3x3"
  return DetectorPeakResponse::peakResolutionFWHM( energy,
                  DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn, nai_fwhm_coefs );
}//float nai_fwhm_fcn( const float energy )


float labr_fwhm_fcn( const float energy )
{
  static const vector<float> labr_fwhm_coefs{ 5.0f, 3.0f, 0.55f };  //"LaBr 10%"
  return DetectorPeakResponse::peakResolutionFWHM( energy,
                  DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn, labr_fwhm_coefs );
}//float labr_fwhm_fcn( const float energy )

  
float hpge_fwhm_fcn( const float energy )
{
  static const vector<float> hpge_fwhm_coefs{ 1.55f, 0.25f, 0.35f };//"HPGe 40%"
  return DetectorPeakResponse::peakResolutionFWHM( energy,
                  DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn, hpge_fwhm_coefs );
}//float hpge_fwhm_fcn( const float energy )

  
CoarseResolutionType coarse_resolution_from_peaks( const vector<shared_ptr<const PeakDef>> &peaks )
{
  size_t num_peaks = 0;
  double max_sig = 0.0;
  double low_w = 0, med_w = 0, high_w = 0, all_w = 0;
  
  for( const auto &p : peaks )
  {
    if( !p || !p->gausPeak() || (p->mean() > 3000) || (p->mean() < 50) )
      continue;
    
    num_peaks += 1;
    
    const double drf_low_fwhm = nai_fwhm_fcn( p->mean() );
    const double drf_med_fwhm = labr_fwhm_fcn( p->mean() );
    const double drf_high_fwhm = hpge_fwhm_fcn( p->mean() );
    
    const double chi_dof = p->chi2dof();
    const double stat_sig = p->peakArea() / p->peakAreaUncert();
    max_sig = std::max( max_sig, stat_sig );
    
    const double w = std::min( stat_sig, 10.0 );// / std::max( chi_dof, 0.5 );
    
    all_w += w;
    
    const double low_diff = fabs(p->fwhm() - drf_low_fwhm);
    const double med_diff = fabs(p->fwhm() - drf_med_fwhm);
    const double high_diff = fabs(p->fwhm() - drf_high_fwhm);
    if( (low_diff < med_diff) && (low_diff < high_diff) )
      low_w += w;
    else if( med_diff < high_diff )
      med_w += w;
    else
      high_w += w;
  }//for( const auto &p : peak_candidates )
  
  if( (num_peaks == 1) && (max_sig < 5) )
    return CoarseResolutionType::Unknown;
  
  if( all_w <= 0.0 )
    return CoarseResolutionType::Unknown;
  
  //if( low_w > high_w )
  //  return CoarseResolutionType::Low;
  
  
  if( (low_w > med_w) && (low_w > high_w) )
    return CoarseResolutionType::Low;
  
  if( (med_w > high_w) )
    return CoarseResolutionType::Medium;
  
  return CoarseResolutionType::High;
}//CoarseResolutionType coarse_resolution_from_peaks( const vector<shared_ptr<const PeakDef>> &peaks )

  
CoarseResolutionType coarse_resolution_from_peaks( const deque<std::shared_ptr<const PeakDef>> &inp )
{
  return coarse_resolution_from_peaks( vector<shared_ptr<const PeakDef>>{begin(inp), end(inp)} );
}
  
bool is_high_res( const std::shared_ptr<const SpecUtils::Measurement> &meas )
{
  // I dont think I've seen a HPGe spectrum straight from the MCA with less than 4096, but we'll be
  //  conservative and allow HPGe to go down to 2048 channels
  const size_t min_high_res_chan = 2048;
  
  // I have seen some poorly tuned/defined, one-off type MCAs have a ton of channels for a low
  //  resolution scintillator, but we'll be reasonable and say low-res has max of 4096 channels,
  //  which should account for (nearly all?) the CZT and LaBr detectors
  const size_t max_low_res_chan = 4096;
  
  // If channel counts leaves things ambiguous, we'll look at the width of channels to make a choice
  //   9 MeV spectrum with 9k channels: 1.098 kev/chan.
  //   Which, 1024 channels * 1.1 --> 1127 keV, which is less than 1400 keV almost all dets go to.
  const float max_high_res_chann_width = 1.1;
                                               
  
  const size_t nchannel = meas ? meas->num_gamma_channels() : 0;
  
  if( nchannel < min_high_res_chan )
    return false;
  
  if( nchannel > max_low_res_chan )
    return true;
  
  const auto cal = meas->energy_calibration();
  if( !cal || !cal->valid()
     || (cal->type() == SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial) )
    return (nchannel > max_low_res_chan);
  
  const float energy_range = cal->upper_energy() - cal->lower_energy();
  const float avrg_channel_width = energy_range / nchannel;
  
  // Note: some (usually low resolution) systems have a sqrt(energy) scale, which invalidates
  //       the average energy assumptions, and havent been looked into yet.
  
  return (avrg_channel_width <= max_high_res_chann_width);
}//bool is_high_res(...)

  
bool is_likely_high_res( InterSpec *viewer )
{
  assert( wApp );
  assert( viewer );
  
  PeakModel *peakModel = viewer->peakModel();
  assert( peakModel );
  if( !peakModel )
    return false;
  
  std::shared_ptr<const SpecMeas> meas = viewer->measurment(SpecUtils::SpectrumType::Foreground);
  assert( meas );
  
  shared_ptr<const SpecUtils::Measurement> foreground = viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  assert( foreground );
  
  if( !meas || !foreground || (foreground->num_gamma_channels() < 512) )
    return false;
  
  switch( meas->detector_type() )
  {
    case SpecUtils::DetectorType::Fulcrum:
    case SpecUtils::DetectorType::Fulcrum40h:
    case SpecUtils::DetectorType::DetectiveUnknown:
    case SpecUtils::DetectorType::DetectiveEx:
    case SpecUtils::DetectorType::DetectiveEx100:
    case SpecUtils::DetectorType::DetectiveEx200:
    case SpecUtils::DetectorType::DetectiveX:
    case SpecUtils::DetectorType::MicroDetective:
      return true;
      
      
    case SpecUtils::DetectorType::Exploranium:
    case SpecUtils::DetectorType::IdentiFinder:
    case SpecUtils::DetectorType::IdentiFinderNG:
    case SpecUtils::DetectorType::IdentiFinderLaBr3:
    case SpecUtils::DetectorType::IdentiFinderTungsten:
    case SpecUtils::DetectorType::IdentiFinderR425NaI:
    case SpecUtils::DetectorType::IdentiFinderR425LaBr:
    case SpecUtils::DetectorType::IdentiFinderR500NaI:
    case SpecUtils::DetectorType::IdentiFinderR500LaBr:
    case SpecUtils::DetectorType::IdentiFinderUnknown:
    case SpecUtils::DetectorType::SAIC8:
    case SpecUtils::DetectorType::MicroRaider:
    case SpecUtils::DetectorType::RadiaCode:
    case SpecUtils::DetectorType::RadHunterNaI:
    case SpecUtils::DetectorType::RadHunterLaBr3:
    case SpecUtils::DetectorType::Rsi701:
    case SpecUtils::DetectorType::Rsi705:
    case SpecUtils::DetectorType::AvidRsi:
    case SpecUtils::DetectorType::OrtecRadEagleNai:
    case SpecUtils::DetectorType::OrtecRadEagleCeBr2Inch:
    case SpecUtils::DetectorType::OrtecRadEagleCeBr3Inch:
    case SpecUtils::DetectorType::OrtecRadEagleLaBr:
    case SpecUtils::DetectorType::Sam940LaBr3:
    case SpecUtils::DetectorType::Sam940:
    case SpecUtils::DetectorType::Sam945:
    case SpecUtils::DetectorType::Srpm210:
    case SpecUtils::DetectorType::RIIDEyeNaI:
    case SpecUtils::DetectorType::RIIDEyeLaBr:
    case SpecUtils::DetectorType::RadSeekerNaI:
    case SpecUtils::DetectorType::RadSeekerLaBr:
    case SpecUtils::DetectorType::VerifinderNaI:
    case SpecUtils::DetectorType::VerifinderLaBr:
    case SpecUtils::DetectorType::KromekD3S:
    case SpecUtils::DetectorType::Sam950:
      return false;
      
    case SpecUtils::DetectorType::Falcon5000: //Any Canberra/Mirion system will be classified as a Falcon 5k
    case SpecUtils::DetectorType::Interceptor:
    case SpecUtils::DetectorType::Unknown:
      break;
  }//switch( meas->detector_type() )
  
  
  shared_ptr<const deque<shared_ptr<const PeakDef>>> fwhmPeaks = peakModel->peaks();
  if( !fwhmPeaks || fwhmPeaks->empty() )
  {
    const set<int> &foreSamples = viewer->displayedSamples(SpecUtils::SpectrumType::Foreground);
    fwhmPeaks = meas->automatedSearchPeaks(foreSamples);
  }
  
  if( fwhmPeaks && fwhmPeaks->size() )
  {
    const vector<shared_ptr<const PeakDef>> peakv( begin(*fwhmPeaks), end(*fwhmPeaks) );
    const auto type = PeakFitUtils::coarse_resolution_from_peaks(peakv);
    return (type == PeakFitUtils::CoarseResolutionType::High);
  }//if( fwhmPeaks && fwhmPeaks->size() )
  
  try
  {
    const double lower_energy = foreground->gamma_channel_lower( 0 );
    const double upper_energy = foreground->gamma_channel_upper( foreground->num_gamma_channels() - 1 );
    const double keV_per_channel = (upper_energy - lower_energy) / foreground->num_gamma_channels();
      
    // 9 MeV spectrum with 9k channels: 1.098 kev/chan.
    // Which, 1024 channels * 1.1 --> 1127 keV, which is less than 1400 keV almost all dets go to.
    const float max_high_res_chann_width = 1.1;
    return (keV_per_channel < max_high_res_chann_width);
  }catch( std::exception & )
  {
    assert( 0 );
  }
  
  return true;
}//bool is_likely_high_res( InterSpec *viewer )
  
}//namespace PeakFitUtils
