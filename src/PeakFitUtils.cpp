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

#include "SpecUtils/StringAlgo.h"

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
  switch( spec->detector_type() )
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
  }//switch( spec->detector_type() )


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


}//namespace PeakFitUtils
