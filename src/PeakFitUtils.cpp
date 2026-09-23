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

#include <cmath>
#include <limits>
#include <memory>
#include <vector>
#include <algorithm>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/EnergyCalibration.h"

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakFitDetPrefs.h"
#include "InterSpec/PeakFitSpecImp.h"
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


float czt_fwhm_fcn( const float energy )
{
  // From Kromek GR1 generic CZT detector
  static const vector<float> czt_fwhm_coefs{ 8.95f, 2.39f, 0.344f };
  return DetectorPeakResponse::peakResolutionFWHM( energy,
                  DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn, czt_fwhm_coefs );
}//float czt_fwhm_fcn( const float energy )


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
  double low_w = 0, labr_w = 0, czt_w = 0, high_w = 0, all_w = 0;

  for( const auto &p : peaks )
  {
    if( !p || !p->gausPeak() || (p->mean() > 3000) || (p->mean() < 50) )
      continue;

    num_peaks += 1;

    const double drf_low_fwhm = nai_fwhm_fcn( p->mean() );
    const double drf_labr_fwhm = labr_fwhm_fcn( p->mean() );
    const double drf_czt_fwhm = czt_fwhm_fcn( p->mean() );
    const double drf_high_fwhm = hpge_fwhm_fcn( p->mean() );

    const double stat_sig = p->peakArea() / p->peakAreaUncert();
    max_sig = std::max( max_sig, stat_sig );

    const double w = std::min( stat_sig, 10.0 );
    all_w += w;

    // Vote for the detector type whose expected FWHM is closest to the measured peak
    const double low_diff = fabs( p->fwhm() - drf_low_fwhm );
    const double labr_diff = fabs( p->fwhm() - drf_labr_fwhm );
    const double czt_diff = fabs( p->fwhm() - drf_czt_fwhm );
    const double high_diff = fabs( p->fwhm() - drf_high_fwhm );

    const double min_diff = std::min( {low_diff, labr_diff, czt_diff, high_diff} );
    if( min_diff == low_diff )
      low_w += w;
    else if( min_diff == labr_diff )
      labr_w += w;
    else if( min_diff == czt_diff )
      czt_w += w;
    else
      high_w += w;
  }//for( const auto &p : peak_candidates )

  if( (num_peaks == 1) && (max_sig < 5) )
    return CoarseResolutionType::Unknown;

  if( all_w <= 0.0 )
    return CoarseResolutionType::Unknown;

  // Return the type with the highest weighted votes
  const double max_w = std::max( {low_w, labr_w, czt_w, high_w} );
  if( max_w == high_w )
    return CoarseResolutionType::High;
  if( max_w == czt_w )
    return CoarseResolutionType::CZT;
  if( max_w == labr_w )
    return CoarseResolutionType::LaBr;

  return CoarseResolutionType::Low;
}//CoarseResolutionType coarse_resolution_from_peaks( const vector<shared_ptr<const PeakDef>> &peaks )

  
CoarseResolutionType coarse_resolution_from_peaks( const deque<std::shared_ptr<const PeakDef>> &inp )
{
  return coarse_resolution_from_peaks( vector<shared_ptr<const PeakDef>>{begin(inp), end(inp)} );
}


CoarseResolutionType classify_det_type(
  const vector<PeakDef> &peaks,
  const shared_ptr<const SpecUtils::Measurement> & )
{
  // Two-stage FWHM-based classification.
  // Stage 1: High-resolution (HPGe) vs NotHigh
  // Stage 2: Low (NaI/CsI) vs MedRes (CZT/LaBr)

  // Helper: 2-way FWHM soft voting with outlier rejection.
  // Returns vote fraction for type A (ref_a).  1-fraction is type B.
  // confidence_out receives the winning vote fraction.
  struct TwoWayResult { double a_vote; double b_vote; int num_passing; };

  const auto two_way_classify = []( const vector<PeakDef> &pks,
    double ref_a_bias, double ref_b_bias,
    double ref_a_fwhm_at_e_fcn_nai_mult, double ref_a_fwhm_at_e_fcn_hpge_mult,
    double ref_b_fwhm_at_e_fcn_nai_mult, double ref_b_fwhm_at_e_fcn_hpge_mult,
    double dist_w_a, double dist_w_b,
    double amp_denom, double wt_clamp,
    double min_e, double max_e, double min_sig ) -> TwoWayResult
  {
    (void)ref_a_bias; (void)ref_b_bias; // unused, folded into the fwhm mults

    struct PV { double w; double a_frac; double b_frac; };
    vector<PV> pvs;
    pvs.reserve( pks.size() );

    double a_w = 0.0, b_w = 0.0;
    int np = 0;

    for( const PeakDef &p : pks )
    {
      if( !p.gausPeak() )
        continue;
      const double e = p.mean();
      if( (e < min_e) || (e > max_e) )
        continue;
      const double amp = p.amplitude();
      if( amp < min_sig )
        continue;
      np += 1;

      const double w = std::min( sqrt( amp ) / amp_denom, wt_clamp );

      const double ref_nai  = nai_fwhm_fcn( static_cast<float>( e ) );
      const double ref_hpge = hpge_fwhm_fcn( static_cast<float>( e ) );
      const double ref_a = ref_nai * ref_a_fwhm_at_e_fcn_nai_mult + ref_hpge * ref_a_fwhm_at_e_fcn_hpge_mult;
      const double ref_b = ref_nai * ref_b_fwhm_at_e_fcn_nai_mult + ref_hpge * ref_b_fwhm_at_e_fcn_hpge_mult;

      const double d_a = fabs( log( p.fwhm() / ref_a ) ) * dist_w_a;
      const double d_b = fabs( log( p.fwhm() / ref_b ) ) * dist_w_b;

      const double eps = 1.0e-6;
      const double inv_a = 1.0 / (d_a * d_a + eps);
      const double inv_b = 1.0 / (d_b * d_b + eps);
      const double inv_s = inv_a + inv_b;

      const double af = inv_a / inv_s;
      const double bf = inv_b / inv_s;

      a_w += w * af;
      b_w += w * bf;
      pvs.push_back( { w, af, bf } );
    }//for( peaks )

    // Outlier rejection (>4 peaks, max 25%)
    if( np > 4 )
    {
      const double mw = std::max( a_w, b_w );
      const int winner = (mw == a_w) ? 0 : 1;
      const size_t max_ex = static_cast<size_t>( np ) / 4;

      vector<pair<double, size_t>> outliers;
      for( size_t i = 0; i < pvs.size(); ++i )
      {
        const double wf = (winner == 0) ? pvs[i].a_frac : pvs[i].b_frac;
        const double bf = std::max( pvs[i].a_frac, pvs[i].b_frac );
        if( wf < bf )
          outliers.push_back( { wf, i } );
      }

      if( !outliers.empty() )
      {
        std::sort( outliers.begin(), outliers.end() );
        const size_t ne = std::min( outliers.size(), max_ex );

        a_w = 0.0;
        b_w = 0.0;
        vector<bool> excl( pvs.size(), false );
        for( size_t i = 0; i < ne; ++i )
          excl[outliers[i].second] = true;

        for( size_t i = 0; i < pvs.size(); ++i )
        {
          if( excl[i] )
            continue;
          a_w += pvs[i].w * pvs[i].a_frac;
          b_w += pvs[i].w * pvs[i].b_frac;
        }
      }
    }//if( outlier rejection )

    return { a_w, b_w, np };
  };//two_way_classify

  // --- Stage 1: HPGe vs NaI ---
  // Type A = HPGe (ref = hpge_bias * hpge_fwhm)
  // Type B = NaI (ref = nai_bias * nai_fwhm)
  // Parameters: defaults, will be updated after GA optimization
  const double s1_hpge_bias = 1.0;
  const double s1_nai_bias = 1.0;
  const double s1_hpge_dist_w = 1.0;
  const double s1_nai_dist_w = 1.0;
  const double s1_amp_denom = 10.0;
  const double s1_wt_clamp = 10.0;
  const double s1_min_e = 50.0;
  const double s1_max_e = 3000.0;
  const double s1_unknown_thr = 0.4;
  const int s1_min_peaks = 2;
  const double s1_min_sig = 2.5;

  const TwoWayResult s1 = two_way_classify( peaks,
    0.0, 0.0,
    0.0, s1_hpge_bias,           // type A (HPGe): 0*nai + hpge_bias*hpge
    s1_nai_bias, 0.0,            // type B (NaI): nai_bias*nai + 0*hpge
    s1_hpge_dist_w, s1_nai_dist_w,
    s1_amp_denom, s1_wt_clamp,
    s1_min_e, s1_max_e, s1_min_sig );

  if( s1.num_passing < s1_min_peaks )
    return CoarseResolutionType::Unknown;

  const double s1_total = s1.a_vote + s1.b_vote;
  if( s1_total <= 0.0 )
    return CoarseResolutionType::Unknown;

  const double s1_max = std::max( s1.a_vote, s1.b_vote );
  const double s1_conf = s1_max / s1_total;

  if( s1_conf < s1_unknown_thr )
    return CoarseResolutionType::Unknown;

  // If HPGe wins Stage 1, we're done
  if( s1_max == s1.a_vote )
    return CoarseResolutionType::High;

  // --- Stage 2: Low (NaI/CsI) vs MedRes (CZT/LaBr) ---
  // Type A = NaI (ref = nai_bias * nai_fwhm)
  // Type B = MedRes (ref = (1-frac)*nai*nai_bias + frac*hpge*hpge_bias)
  // Parameters: defaults, will be updated after GA optimization
  const double s2_nai_bias = 1.0;
  const double s2_medres_frac = 0.5;
  const double s2_hpge_bias = 1.0;
  const double s2_nai_dist_w = 1.0;
  const double s2_medres_dist_w = 1.0;
  const double s2_amp_denom = 10.0;
  const double s2_wt_clamp = 10.0;
  const double s2_min_e = 50.0;
  const double s2_max_e = 3000.0;
  const double s2_unknown_thr = 0.4;
  const int s2_min_peaks = 2;
  const double s2_min_sig = 2.5;

  const TwoWayResult s2 = two_way_classify( peaks,
    0.0, 0.0,
    s2_nai_bias, 0.0,  // type A (NaI): nai_bias*nai + 0*hpge
    (1.0 - s2_medres_frac) * s2_nai_bias, s2_medres_frac * s2_hpge_bias, // type B (MedRes)
    s2_nai_dist_w, s2_medres_dist_w,
    s2_amp_denom, s2_wt_clamp,
    s2_min_e, s2_max_e, s2_min_sig );

  if( s2.num_passing < s2_min_peaks )
    return CoarseResolutionType::Unknown;

  const double s2_total = s2.a_vote + s2.b_vote;
  if( s2_total <= 0.0 )
    return CoarseResolutionType::Unknown;

  const double s2_max = std::max( s2.a_vote, s2.b_vote );
  const double s2_conf = s2_max / s2_total;

  if( s2_conf < s2_unknown_thr )
    return CoarseResolutionType::Unknown;

  return (s2_max == s2.a_vote) ? CoarseResolutionType::Low : CoarseResolutionType::MedRes;
}//classify_det_type


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

  const std::shared_ptr<PeakModel> peakModel = viewer->peakModel();
  assert( peakModel );
  if( !peakModel )
    return false;

  std::shared_ptr<const SpecMeas> meas = viewer->measurment(SpecUtils::SpectrumType::Foreground);
  //assert( meas );

  shared_ptr<const SpecUtils::Measurement> foreground = viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  //assert( foreground );

  if( !meas || !foreground || (foreground->num_gamma_channels() < 512) )
    return false;

  // If the user has set PeakFitDetPrefs, use that as the definitive answer
  const shared_ptr<const PeakFitDetPrefs> fitPrefs = meas->peakFitDetPrefs();
  if( fitPrefs && (fitPrefs->m_det_type != CoarseResolutionType::Unknown) )
    return (fitPrefs->m_det_type == CoarseResolutionType::High);



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
    case SpecUtils::DetectorType::RadiaCodeCsI10:
    case SpecUtils::DetectorType::RadiaCodeCsI14:
    case SpecUtils::DetectorType::RadiaCodeGAGG10:
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
    case SpecUtils::DetectorType::H3D400:
    case SpecUtils::DetectorType::KromekD3S:
    case SpecUtils::DetectorType::Sam950:
    case SpecUtils::DetectorType::KromekD5:
    case SpecUtils::DetectorType::KromekGR1:
    case SpecUtils::DetectorType::Raysid:
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
  
CoarseResolutionType coarse_type_for_detector_type( const SpecUtils::DetectorType type )
{
  switch( type )
  {
    // HPGe detectors
    case SpecUtils::DetectorType::DetectiveUnknown:
    case SpecUtils::DetectorType::DetectiveEx:
    case SpecUtils::DetectorType::DetectiveEx100:
    case SpecUtils::DetectorType::DetectiveEx200:
    case SpecUtils::DetectorType::DetectiveX:
    case SpecUtils::DetectorType::MicroDetective:
    case SpecUtils::DetectorType::Falcon5000:
    case SpecUtils::DetectorType::Fulcrum:
    case SpecUtils::DetectorType::Fulcrum40h:
      return CoarseResolutionType::High;

    // NaI/CsI detectors (Low)
    case SpecUtils::DetectorType::IdentiFinder:
    case SpecUtils::DetectorType::IdentiFinderNG:
    case SpecUtils::DetectorType::IdentiFinderR425NaI:
    case SpecUtils::DetectorType::IdentiFinderR500NaI:
    case SpecUtils::DetectorType::IdentiFinderUnknown:
    case SpecUtils::DetectorType::Exploranium:
    case SpecUtils::DetectorType::SAIC8:
    case SpecUtils::DetectorType::RadHunterNaI:
    case SpecUtils::DetectorType::Rsi701:
    case SpecUtils::DetectorType::Rsi705:
    case SpecUtils::DetectorType::AvidRsi:
    case SpecUtils::DetectorType::OrtecRadEagleNai:
    case SpecUtils::DetectorType::Sam940:
    case SpecUtils::DetectorType::Sam945:
    case SpecUtils::DetectorType::Srpm210:
    case SpecUtils::DetectorType::RIIDEyeNaI:
    case SpecUtils::DetectorType::RadSeekerNaI:
    case SpecUtils::DetectorType::VerifinderNaI:
    case SpecUtils::DetectorType::KromekD3S:
    case SpecUtils::DetectorType::RadiaCodeCsI10:
    case SpecUtils::DetectorType::RadiaCodeCsI14:
      return CoarseResolutionType::Low;

    // LaBr3 detectors
    case SpecUtils::DetectorType::IdentiFinderLaBr3:
    case SpecUtils::DetectorType::IdentiFinderR425LaBr:
    case SpecUtils::DetectorType::IdentiFinderR500LaBr:
    case SpecUtils::DetectorType::RadHunterLaBr3:
    case SpecUtils::DetectorType::OrtecRadEagleLaBr:
    case SpecUtils::DetectorType::Sam940LaBr3:
    case SpecUtils::DetectorType::RIIDEyeLaBr:
    case SpecUtils::DetectorType::RadSeekerLaBr:
    case SpecUtils::DetectorType::VerifinderLaBr:
      return CoarseResolutionType::LaBr;

    // CZT detectors
    case SpecUtils::DetectorType::KromekGR1:
    case SpecUtils::DetectorType::MicroRaider:
    case SpecUtils::DetectorType::Interceptor:
    case SpecUtils::DetectorType::H3D400:
      return CoarseResolutionType::CZT;

    // CeBr3 (medium resolution)
    case SpecUtils::DetectorType::OrtecRadEagleCeBr2Inch:
    case SpecUtils::DetectorType::OrtecRadEagleCeBr3Inch:
      return CoarseResolutionType::MedRes;

    // Detectors whose type is ambiguous or unknown — fall through to string/FWHM check
    case SpecUtils::DetectorType::Unknown:
    case SpecUtils::DetectorType::Raysid:
    case SpecUtils::DetectorType::Sam950:
    case SpecUtils::DetectorType::KromekD5:
    case SpecUtils::DetectorType::IdentiFinderTungsten:
    case SpecUtils::DetectorType::RadiaCodeGAGG10:
      break;
  }//switch( type )

  return CoarseResolutionType::Unknown;
}//coarse_type_for_detector_type(...)


CoarseResolutionType coarse_det_type(
  const shared_ptr<const SpecUtils::Measurement> &meas,
  const shared_ptr<const SpecMeas> &spec )
{
  // Lambda to perform FWHM-based classification (Tier 3)
  const auto fwhm_classify = [&meas]() -> CoarseResolutionType {
    if( meas && (meas->num_gamma_channels() >= 16) )
    {
      const PeakFitSpec::SpecClassType spec_class = PeakFitSpec::initial_lowres_highres_classify( meas );
      switch( spec_class )
      {
        case PeakFitSpec::SpecClassType::High:
          return CoarseResolutionType::High;
        case PeakFitSpec::SpecClassType::LowOrMedRes:
          return CoarseResolutionType::LowOrMedRes;
        case PeakFitSpec::SpecClassType::Unknown:
          return is_high_res( meas ) ? CoarseResolutionType::High : CoarseResolutionType::LowOrMedRes;
      }
    }
    return CoarseResolutionType::Unknown;
  };//fwhm_classify lambda
  
  // If no SpecMeas provided, can only do FWHM-based classification
  if( !spec )
    return fwhm_classify();

  // Tier 1: Check SpecUtils::DetectorType from file parsing (most reliable)
  {
    const CoarseResolutionType from_type = coarse_type_for_detector_type( spec->detector_type() );
    if( from_type != CoarseResolutionType::Unknown )
      return from_type;
  }


  // Tier 2: Search instrument metadata strings for detector-type keywords
  {
    const string &inst_type = spec->instrument_type();
    const string &inst_model = spec->instrument_model();
    const string &inst_id = spec->instrument_id();

    // Collect all remarks into one string for searching
    string all_remarks;
    for( const string &r : spec->remarks() )
    {
      all_remarks += " ";
      all_remarks += r;
    }

    // Check each metadata field for detector type keywords (case-insensitive)
    const auto check_strings = [&]( const string &keyword ) -> bool {
      return SpecUtils::icontains( inst_type, keyword )
          || SpecUtils::icontains( inst_model, keyword )
          || SpecUtils::icontains( inst_id, keyword )
          || SpecUtils::icontains( all_remarks, keyword );
    };

    // HPGe indicators
    if( check_strings( "HPGe" ) || check_strings( "Germanium" )
       || check_strings( "coaxial" ) )
    {
      return CoarseResolutionType::High;
    }

    // LaBr indicators (check before NaI since some models have both)
    if( check_strings( "LaBr" ) || check_strings( "lanthanum bromide" ) )
      return CoarseResolutionType::LaBr;

    // CZT indicators
    if( check_strings( "CZT" ) || check_strings( "CdZnTe" )
       || check_strings( "cadmium zinc telluride" ) )
    {
      return CoarseResolutionType::CZT;
    }

    // NaI indicators
    if( check_strings( "NaI" ) || check_strings( "sodium iodide" )
       || check_strings( "Nal" ) )
    {
      return CoarseResolutionType::Low;
    }

    // CsI indicators
    if( check_strings( "CsI" ) || check_strings( "cesium iodide" ) )
      return CoarseResolutionType::Low;
  }


  // Tier 3: FWHM-based classification using GA-optimized settings
  return fwhm_classify();
}//coarse_det_type


CoarseResolutionType effective_det_type(
  const std::shared_ptr<const PeakFitDetPrefs> &prefs,
  const std::shared_ptr<const SpecUtils::Measurement> &meas,
  const std::shared_ptr<const SpecMeas> &spec )
{
  if( prefs && (prefs->m_det_type != CoarseResolutionType::Unknown) )
    return prefs->m_det_type;
  return coarse_det_type( meas, spec );
}


namespace
{
/** Lower spectroscopic-extent estimator ("C1"): a resolution-aware, statistically-
 thresholded refinement of the old fixed-window 2nd-derivative heuristic.

 The energy at which a gamma spectrum "turns on" (below it is electronic-noise / low-energy
 junk, above it is real spectroscopic data) is found in three stages, all expressed in
 resolution-aware units (multiples of the expected low-energy FWHM in channels) rather than
 magic channel constants:

  1. State machine on the Savitzky-Golay 2nd derivative (window ~ 0.5 FWHM), thresholded in
     units of the propagated Poisson sigma.  Finds the first significant negative curvature
     (the turn-on), advances to where it levels off, and distinguishes a genuine x-ray peak
     sitting on the turn-on (back up before it) from an electronic junk bump (a real dip
     followed by a higher persistent continuum - step past it).
  2. climb_steps: a low flat junk floor followed by a SHARP (>= climb_r x) persistent step
     up means everything below the step is junk (dominant on planar HPGe) - climb it.
  3. persistence_walk: advance up through the extended low-energy junk band (fluorescence /
     backscatter clutter) to the sustained continuum, stopping where the spectrum rises back
     up (the continuum) rather than dipping further.

 The extent is capped at an absolute energy (the junk band is bounded in absolute energy),
 and a walk that climbs into the spectrum's main body (a broad bremsstrahlung / scatter
 hump with no clean turn-on) is rejected back to the turn-on top.  Tuned and cross-validated
 against ~960 hand-picked lower-extent bands across six detector types (HPGe coax & planar,
 NaI, LaBr3, CZT, D3S): pooled band hit-rate ~77% vs ~29% for the previous heuristic.

 Returns the lower-extent channel index, or `nbin` on failure (no turn-on found in the
 lower third of the spectrum).  Search is limited to the first nbin/3 channels, preserving
 the historical failure contract.
*/
struct LowerExtentFinder
{
  // Tunables (cross-validated defaults; see function-level comment).
  const double m_f = 0.5;          // SG window = f * fwhm_low_ch
  const double m_t1 = 1.5;         // phase-1 negative-curvature threshold (sigma)
  const double m_t2 = 1.0;         // level-off threshold (sigma)
  const double m_t3 = 2.0;         // oscillation / peak threshold (sigma)
  const double m_tb = 1.0;         // peak left-extent back-up threshold (sigma)
  const double m_persist_q = 0.85; // persistence-walk dip factor
  const double m_climb_r = 3.0;    // junk-plateau sharp-step ratio
  const double m_e_cap_kev = 80.0; // absolute-energy ceiling on the extent
  const double m_body_frac = 0.60; // reject a walk landing above this fraction of the peak

  const vector<float> &m_counts;
  const size_t m_n;
  bool m_hr = false;               // is_high_res (gates the junk-plateau climb)
  size_t m_limit = 0;              // min( nbin/5, channel(e_cap) )
  double m_fw = 1.0;               // fwhm_low_ch
  int m_Wb = 2;                    // box-average half-context window
  int m_Hb = 6;                    // look-ahead horizon (channels)
  vector<float> m_smoothed, m_sigma, m_sm_box;

  LowerExtentFinder( const shared_ptr<const SpecUtils::Measurement> &meas )
  : m_counts( *meas->gamma_counts() ), m_n( meas->gamma_counts()->size() )
  {
    // Resolution scale available from just the Measurement.
    m_hr = PeakFitUtils::is_high_res( meas );
    const bool hr = m_hr;
    const PeakFitUtils::CoarseResolutionType det_type = hr
        ? PeakFitUtils::CoarseResolutionType::High
        : PeakFitUtils::CoarseResolutionType::LowOrMedRes;
    float sig_lo = 0.0f, sig_hi = 0.0f;
    expected_peak_width_limits( 60.0f, det_type, meas, sig_lo, sig_hi );
    const double fwhm_kev = 2.35482 * std::sqrt( std::max( sig_lo*sig_hi, 0.0f ) );

    const size_t ch60 = meas->find_gamma_channel( 60.0f );
    const double kev_per_ch = std::max( 1.0e-6,
        static_cast<double>( meas->gamma_channel_upper(ch60) - meas->gamma_channel_lower(ch60) ) );
    m_fw = std::max( fwhm_kev / kev_per_ch, 1.0 );

    // Clamp the SG half-window so the full window (2*W+1) never exceeds the spectrum
    // length - smooth_with_variance throws when nCoeffs > nSamples.  Only binds for
    // pathologically short spectra (nbin < 2*W+1); real spectra are far larger.
    const int W = std::min( clamp_win( m_f * m_fw, 2, 16 ),
                            static_cast<int>( (m_n - 1)/2 ) );
    SavitzyGolayCoeffs sg( W, W, 2, 2 );
    vector<float> var;
    sg.smooth_with_variance( m_counts, m_smoothed, var );
    m_sigma.resize( m_n );
    for( size_t i = 0; i < m_n; ++i )
      m_sigma[i] = (var[i] > 0.0f) ? std::sqrt(var[i]) : 1.0f;

    m_Wb = clamp_win( m_fw, 2, 16 );
    m_Hb = std::max( static_cast<int>( std::lround(10.0*m_fw) ), 3*m_Wb );
    m_sm_box = box_average_same( m_counts, m_Wb );

    const size_t cap_ch = meas->find_gamma_channel( static_cast<float>(m_e_cap_kev) );
    m_limit = std::min( m_n/5, cap_ch );
  }

  static int clamp_win( double x, int lo, int hi )
  {
    const int v = static_cast<int>( std::lround(x) );
    return std::max( lo, std::min( hi, v ) );
  }

  // numpy np.convolve(counts, ones(W)/W, mode="same"): centred moving average with zero
  // padding; for even W the window is W/2 channels left and W/2-1 right of the point.
  static vector<float> box_average_same( const vector<float> &a, int W )
  {
    const int N = static_cast<int>( a.size() );
    const int left = W/2, right = W-1-W/2;
    vector<float> out( N, 0.0f );
    for( int k = 0; k < N; ++k )
    {
      double s = 0.0;
      for( int d = -left; d <= right; ++d )
      {
        const int j = k + d;
        if( j >= 0 && j < N )
          s += a[j];
      }
      out[k] = static_cast<float>( s / W );
    }
    return out;
  }

  // Median of m_sm_box over [lo, hi) (numpy averages the two central values for even n).
  double sm_box_median( size_t lo, size_t hi ) const
  {
    if( hi <= lo )
      return 0.0;
    vector<float> tmp( m_sm_box.begin()+lo, m_sm_box.begin()+hi );
    std::sort( tmp.begin(), tmp.end() );
    const size_t m = tmp.size();
    return (m % 2) ? tmp[m/2] : 0.5*(tmp[m/2 - 1] + tmp[m/2]);
  }

  // One pass of the 2nd-derivative state machine from `start`.  Returns the lower channel
  // (or m_n on no-detection) and sets `peak_path` when a peak/spike branch was taken.
  size_t scan( size_t start, bool &peak_path, int depth = 0 ) const
  {
    peak_path = false;
    const size_t n = m_n;

    // Phase 1: first significant negative curvature at/after start.
    size_t ch = start;
    while( ch < m_limit && !(m_smoothed[ch] < -m_t1*m_sigma[ch]) )
      ++ch;
    if( ch >= m_limit )
      return n;
    const size_t first_neg = ch;

    // Phase 2: advance until curvature levels off, or swings significantly positive (peak).
    bool found_osc = false;
    while( ch < n )
    {
      if( std::fabs(m_smoothed[ch]) < m_t2*m_sigma[ch] )
        break;
      if( m_smoothed[ch] > m_t3*m_sigma[ch] )
      {
        found_osc = true;
        break;
      }
      ++ch;
    }

    // Verify the oscillation is a real peak (counts fall after) vs a turn-on kink.
    if( found_osc && (ch > first_neg) )
    {
      const size_t hw = std::min<size_t>( 3, (ch - first_neg)/2 );
      double sum_neg = 0.0, sum_osc = 0.0;
      for( size_t i = 0; i <= 2*hw; ++i )
      {
        const size_t ni = std::min( n-1, (first_neg >= hw) ? (first_neg - hw + i) : i );
        const size_t oi = std::min( n-1, (ch >= hw) ? (ch - hw + i) : i );
        sum_neg += m_counts[ni];
        sum_osc += m_counts[oi];
      }
      if( sum_osc >= sum_neg )
        found_osc = false;
    }

    if( found_osc )
    {
      // Distinguish a real x-ray peak on the turn-on (towers over what follows) from an
      // electronic junk bump (a real dip followed by a higher persistent continuum).
      const size_t osc_ch = ch;
      const size_t pk_lo = (first_neg >= size_t(m_Wb)) ? (first_neg - m_Wb) : 0;
      double peak_top = 0.0;
      for( size_t i = pk_lo; i <= osc_ch && i < n; ++i )
        peak_top = std::max( peak_top, static_cast<double>(m_sm_box[i]) );

      // The junk-bump-vs-peak recursion is a high-resolution phenomenon (electronic bumps
      // below the real continuum on HPGe).  On coarse detectors the turn-on is sharp with
      // real X-ray peaks right above it; cascading the recursion past them overshoots the
      // extent, so gate it to high-resolution spectra.
      const size_t v_hi = std::min( osc_ch + m_Hb, m_limit );
      if( m_hr && (depth < 3) && (v_hi > osc_ch) )
      {
        size_t v_idx = osc_ch;
        for( size_t i = osc_ch; i <= v_hi; ++i )
          if( m_sm_box[i] < m_sm_box[v_idx] )
            v_idx = i;
        const double valley = m_sm_box[v_idx];
        const size_t a_hi = std::min( v_idx + m_Hb, m_limit );
        const double after_med = (a_hi > v_idx) ? sm_box_median( v_idx, a_hi + 1 ) : 0.0;
        if( (valley < 0.5*peak_top) && (after_med >= peak_top) )
          return scan( v_idx, peak_path, depth + 1 );  // junk bump; real spectrum is above
      }

      // Real peak at turn-on: back up to before the peak.
      size_t left = first_neg;
      while( left > 0 && (m_smoothed[left] <= -m_tb*m_sigma[left]) )
        --left;
      const size_t neg_hw = (first_neg > left) ? (first_neg - left) : 1;
      const size_t buffered = (left > neg_hw) ? (left - neg_hw) : 0;

      size_t first_inf = start;
      while( first_inf < first_neg && (std::fabs(m_smoothed[first_inf]) <= m_t3*m_sigma[first_inf]) )
        ++first_inf;
      while( (first_inf + 3) < n
             && !( (m_counts[first_inf] > 0.0f) && (m_counts[first_inf+1] > 0.0f)
                   && (m_counts[first_inf+2] > 0.0f) && (m_counts[first_inf+3] > 0.0f) ) )
        ++first_inf;

      peak_path = true;
      return std::max( first_inf, buffered );
    }

    const size_t leveled = std::min( ch, n-1 );
    // Narrow turn-on spike with a valley just after (scintillator LLD artifact).
    const size_t p_lo = (first_neg >= 2) ? (first_neg - 2) : 0;
    const size_t p_hi = std::min( first_neg + 2, n-1 );
    if( (leveled > first_neg) && ((leveled - first_neg) < 10) )
    {
      double max_near = 0.0, min_after = std::numeric_limits<double>::max();
      for( size_t i = p_lo; i <= p_hi; ++i )
        max_near = std::max( max_near, static_cast<double>(m_counts[i]) );
      for( size_t i = first_neg + 1; i <= leveled && i < n; ++i )
        min_after = std::min( min_after, static_cast<double>(m_counts[i]) );
      if( (max_near > 5.0) && (min_after < 0.5*max_near) )
      {
        size_t pk = p_lo;
        for( size_t i = p_lo; i <= p_hi; ++i )
          if( m_counts[i] > m_counts[pk] )
            pk = i;
        size_t lower = pk;
        for( int i = 0; i < 3; ++i )
        {
          if( (lower == 0) || (m_counts[lower-1] == 0.0f) )
            break;
          --lower;
        }
        peak_path = true;
        return lower;
      }
    }
    return leveled;
  }

  // Junk-plateau climb: a low flat floor followed by a sharp, persistent step up.
  size_t climb_steps( size_t lower ) const
  {
    for( int iter = 0; iter < 3; ++iter )
    {
      const double level = std::max( sm_box_median( lower, std::min(lower + m_Wb + 1, m_n) ), 1.0e-9 );
      const size_t hi = std::min( lower + m_Hb, m_limit );
      bool have = false;
      size_t found = 0;
      for( size_t k = lower + m_Wb; k <= hi; ++k )
      {
        const double lk = sm_box_median( k, std::min(k + m_Wb + 1, m_n) );
        const size_t pre_lo = (k >= size_t(2*m_Wb)) ? (k - 2*m_Wb) : 0;
        const size_t pre_hi = std::max( (k >= size_t(m_Wb)) ? (k - m_Wb) : size_t(1), size_t(1) );
        const double pre = sm_box_median( pre_lo, pre_hi );
        if( (lk >= m_climb_r*level) && (lk >= m_climb_r*std::max(pre,1.0e-9)) )
        {
          const size_t seg_hi = std::min( k + m_Hb, m_limit );
          double seg_min = std::numeric_limits<double>::max();
          for( size_t i = k; i <= seg_hi; ++i )
            seg_min = std::min( seg_min, static_cast<double>(m_sm_box[i]) );
          if( (seg_hi >= k) && (seg_min >= 0.5*lk) )
          {
            have = true;
            found = k;
            break;
          }
        }
      }
      if( !have )
        break;
      const double lk = sm_box_median( found, std::min(found + m_Wb + 1, m_n) );
      while( (found + 1 <= m_limit) && (m_sm_box[found] < 0.9*lk) )
        ++found;
      lower = std::min( found + static_cast<size_t>(std::lround(m_fw)), m_limit );
    }
    return lower;
  }

  // Advance through the low-energy junk band to the sustained continuum.
  size_t persistence_walk( size_t k0 ) const
  {
    const double r_gate = 1.5;
    const double peak_t = 4.0;    // stop before crossing a significant peak (real data)
    const size_t limit = m_n/5;   // matches the standalone prototype's internal bound
    size_t k = k0;
    for( int iter = 0; iter < 200; ++iter )
    {
      const size_t hi = std::min( k + m_Hb, limit );
      if( k + 1 >= hi )
        break;
      // Segment (k, hi]; truncate at the first channel that rises above r_gate*level.
      size_t seg_hi = hi;
      const double thresh = r_gate * std::max( static_cast<double>(m_sm_box[k]), 1.0e-9 );
      for( size_t i = k+1; i <= hi; ++i )
      {
        if( m_sm_box[i] > thresh )
        {
          seg_hi = i - 1;     // pre-rise region only
          break;
        }
      }
      if( seg_hi < k+1 )
        break;
      size_t m_idx = k+1;
      for( size_t i = k+1; i <= seg_hi; ++i )
        if( m_sm_box[i] < m_sm_box[m_idx] )
          m_idx = i;
      // Peak-stop: never cross a statistically-significant peak (2nd deriv < -peak_t*sigma)
      // between the current channel and the candidate - a real low-energy source line is
      // spectroscopic data to keep, so the extent must stay below it.
      bool peak_ahead = false;
      for( size_t i = k+1; i <= m_idx; ++i )
      {
        if( m_smoothed[i] < -peak_t * m_sigma[i] )
        {
          peak_ahead = true;
          break;
        }
      }
      if( peak_ahead )
        break;
      if( m_sm_box[k] <= m_sm_box[m_idx] / m_persist_q )
        break;
      k = m_idx;
    }
    return std::min( k, limit );
  }

  size_t find()
  {
    bool peak_path = false;
    size_t lower = scan( 0, peak_path );
    if( lower >= m_n )
      return m_n;

    // The climb handles a low flat junk plateau followed by a sharp step up - a genuinely
    // HIGH-RESOLUTION phenomenon (planar HPGe).  On coarse detectors a steep but smooth
    // turn-on rise trips the sharp-step test spuriously, so gate the climb to high-res.
    if( (m_climb_r > 0.0) && m_hr )
      lower = climb_steps( lower );

    if( !peak_path && (m_persist_q > 0.0) )
    {
      const size_t near_gap = m_Hb;  // NB: 'near' is a <windows.h> macro on MSVC, so avoid it
      const size_t start_lower = lower;
      bool clamped = false;
      for( int iter = 0; iter < 3; ++iter )
      {
        const size_t raw = persistence_walk( lower );
        clamped = (raw > m_limit);          // walk had no clean stop below the cap
        const size_t walked = std::min( raw, m_limit );
        if( walked <= lower )
          break;
        bool nxt_peak = false;
        const size_t nxt = scan( walked, nxt_peak );
        // Reject a re-detected edge that leveled off PAST a statistically-significant peak
        // (a real source line the extent must stay below) - otherwise the loop climbs onto
        // that peak and marches to the cap.
        bool crossed_peak = false;
        if( (nxt < m_n) && (nxt > walked) )
        {
          for( size_t i = walked; (i <= nxt) && (i < m_n); ++i )
          {
            if( m_smoothed[i] < -4.0 * m_sigma[i] )
            {
              crossed_peak = true;
              break;
            }
          }
        }
        if( (nxt >= m_n) || nxt_peak || (nxt <= walked) || ((nxt - walked) > near_gap) || crossed_peak )
        {
          lower = walked;
          break;
        }
        lower = nxt;
      }
      // Reject a walk that climbed into the spectrum's main body (broad hump, no turn-on).
      if( lower > start_lower )
      {
        // A walk clamped at the absolute cap never found a clean stopping point below it
        // (the broad-hump / bremsstrahlung signature) - revert to the turn-on top.
        if( clamped )
        {
          lower = start_lower;
        }else
        {
          // Peak level is estimated over a stable wide range (nbin/3), NOT the extent
          // search limit - otherwise lowering the absolute cap would shrink this window
          // and make the body-fraction test misfire, reverting good long HPGe walks.
          const size_t body_hi = std::min( m_n/3, m_n-1 );
          double peak_level = 0.0;
          for( size_t i = 0; i <= body_hi; ++i )
            peak_level = std::max( peak_level, static_cast<double>(m_sm_box[i]) );
          if( (peak_level > 0.0) && (m_sm_box[lower] > m_body_frac*peak_level) )
            lower = start_lower;
        }
      }
    }

    if( lower > m_limit )
      return m_n;
    return lower;
  }
};//struct LowerExtentFinder
}//namespace


bool find_spectroscopic_extent( std::shared_ptr<const SpecUtils::Measurement> meas,
                               size_t &lower_channel,
                               size_t &upper_channel )
{
  // A valid energy calibration is required: the estimator sizes its windows from the
  // expected FWHM at 60 keV and caps the extent at an absolute energy, both of which need
  // channel<->energy mapping (find_gamma_channel throws without it).
  if( !meas || !meas->gamma_counts()
     || !meas->energy_calibration() || !meas->energy_calibration()->valid() )
    return false;

  const vector<float> &channel_counts = *meas->gamma_counts();
  const size_t nbin = channel_counts.size();
  if( nbin < 7 )
    return false;

  // Lower extent: resolution-aware "C1" estimator (see LowerExtentFinder above).
  const size_t detected_lower = LowerExtentFinder( meas ).find();

  if( detected_lower >= nbin || detected_lower > (nbin/5) )
  {
    lower_channel = upper_channel = 0;
    return false;
  }//
  lower_channel = detected_lower;

  size_t upperlastchannel = nbin - 1;
  while( upperlastchannel > 0 && channel_counts[upperlastchannel] < 5.0f )
    --upperlastchannel;
  
  //Start at the turn on energy, but make sure were in a region of the spectra
  //  with decent statistics (e.g. at least 20 counts per bin)
  size_t lastchannel = lower_channel;
  while( lastchannel < nbin && channel_counts[lastchannel] < 20.0 )
    ++lastchannel;
  
  //Now find where the spectrum drops below X counts per channel for a number
  //  of channels in a row
  int numbelow = 0;
  const int nbelowlimit = std::max( int(std::ceil(0.0015*nbin)), 3 );
  while( lastchannel < (nbin-1) && numbelow <= nbelowlimit )
    numbelow = (channel_counts[++lastchannel] > 0.1) ? 0 : numbelow + 1;
  
  upper_channel = (lastchannel == nbin) ? size_t(nbin-1) : lastchannel;
  upper_channel = std::max( upper_channel, upperlastchannel+1 );
  upper_channel = std::min( upper_channel, nbin-1 );
  
  return true;
}//bool find_spectroscopic_extent(...)

}//namespace PeakFitUtils
