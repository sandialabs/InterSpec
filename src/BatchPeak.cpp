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

#include <set>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <exception>
#include <chrono>
#include <iomanip>
#include <sstream>

#include "external_libs/SpecUtils/3rdparty/inja/inja.hpp"

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/BatchPeak.h"
#include "InterSpec/EnergyCal.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/BatchInfoLog.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/BatchSampleSelect.h"
#include "InterSpec/DetectionLimitCalc.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"


using namespace std;

const char * const BatchPeak::BatchPeakFitOptions::sm_report_display_name_marker = ":--DisplayName--:";

namespace
{
  std::string get_report_template_name( const std::string &template_path )
  {
    const size_t pos = template_path.find( BatchPeak::BatchPeakFitOptions::sm_report_display_name_marker );
    if( pos == std::string::npos )
      return template_path;
    return template_path.substr( pos + strlen( BatchPeak::BatchPeakFitOptions::sm_report_display_name_marker ) );
  }

  std::string get_report_template_path( const std::string &template_path )
  {
    const size_t pos = template_path.find( BatchPeak::BatchPeakFitOptions::sm_report_display_name_marker );
    if( pos == std::string::npos )
      return template_path;
    return template_path.substr( 0, pos );
  }

  /** Returns the FWHM to use for an exemplar peak; for peaks that arent Gaussian (i.e., "data
   defined" peaks), the ROI is assumed to be the usual 2.5 FWHM wide.

   Returns zero if a width could not be determined.
   */
  double exemplar_peak_fwhm( const PeakDef &peak )
  {
    if( peak.gausPeak() && (peak.fwhm() > 0.0) )
      return peak.fwhm();

    const double roi_width = peak.upperX() - peak.lowerX();

    return (roi_width > 0.0) ? (roi_width / 2.5) : 0.0;
  }//double exemplar_peak_fwhm( const PeakDef & )


  /** Returns the energy range the Currie-style limit used, including the side channels used to
   estimate the continuum.
   */
  std::pair<double,double> currie_energy_range( const DetectionLimitCalc::CurrieMdaResult &result )
  {
    const std::shared_ptr<const SpecUtils::Measurement> &spec = result.input.spectrum;
    assert( spec && spec->energy_calibration() );

    const size_t first_channel = (result.input.num_lower_side_channels > 0)
                                    ? result.first_lower_continuum_channel
                                    : result.first_peak_region_channel;
    const size_t last_channel = (result.input.num_upper_side_channels > 0)
                                    ? result.last_upper_continuum_channel
                                    : result.last_peak_region_channel;

    return { spec->gamma_channel_lower(first_channel), spec->gamma_channel_upper(last_channel) };
  }//std::pair<double,double> currie_energy_range(...)


  /** Creates a fixed-geometry detector response with unity intrinsic efficiency, so that the
   "activity" limited by `DetectionLimitCalc::get_activity_or_distance_limits` is simply the
   number of peak counts.

   The FWHM curve is taken from `base_drf` if it has one, otherwise it is fit from the exemplar
   peaks, and failing that it is set to a constant.  In practice which curve is used doesnt matter,
   since each peaks FWHM is passed into the calculation explicitly; the curve just has to exist for
   the DRF to be accepted.

   Returns nullptr if no usable peak width could be determined.
   */
  std::shared_ptr<const DetectorPeakResponse> make_unity_fixed_geom_drf(
                        const std::vector<std::shared_ptr<const PeakDef>> &exemplar_peaks,
                        const std::shared_ptr<const SpecUtils::Measurement> &spectrum,
                        const std::shared_ptr<const DetectorPeakResponse> &base_drf )
  {
    auto drf = std::make_shared<DetectorPeakResponse>();

    try
    {
      drf->setIntrinsicEfficiencyFormula( "1.0", 2.54f*PhysicalUnits::cm, PhysicalUnits::keV,
                                     0.0f, 0.0f,
                                     DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct );
    }catch( std::exception & )
    {
      return nullptr;
    }

    if( base_drf && base_drf->hasResolutionInfo() )
    {
      try
      {
        drf->setFwhmCoefficients( base_drf->resolutionFcnCoefficients(), base_drf->resolutionFcnType() );
      }catch( std::exception & )
      {
      }
    }//if( base_drf && base_drf->hasResolutionInfo() )

    if( !drf->hasResolutionInfo() )
    {
      try
      {
        auto peaks = std::make_shared<std::deque<std::shared_ptr<const PeakDef>>>(
                                              begin(exemplar_peaks), end(exemplar_peaks) );
        drf->fitResolution( peaks, spectrum, DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial );
      }catch( std::exception & )
      {
      }
    }//if( !drf->hasResolutionInfo() )

    if( !drf->hasResolutionInfo() )
    {
      // Fall back to a constant width, taken as the median of the exemplar peak widths.
      std::vector<float> fwhms;
      for( const std::shared_ptr<const PeakDef> &p : exemplar_peaks )
      {
        const double fwhm = p ? exemplar_peak_fwhm(*p) : 0.0;
        if( fwhm > 0.0 )
          fwhms.push_back( static_cast<float>(fwhm) );
      }

      if( fwhms.empty() )
        return nullptr;

      std::sort( begin(fwhms), end(fwhms) );
      const float fwhm = fwhms[fwhms.size()/2];

      try
      {
        // For `kSqrtPolynomial`, FWHM = sqrt(a0 + a1*E_MeV + ...), so this gives a constant width.
        drf->setFwhmCoefficients( std::vector<float>{ fwhm*fwhm, 0.0f },
                                  DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial );
      }catch( std::exception & )
      {
        return nullptr;
      }
    }//if( !drf->hasResolutionInfo() )

    if( !drf->isValid() || !drf->hasResolutionInfo() )
      return nullptr;

    return drf;
  }//make_unity_fixed_geom_drf(...)
}//namespace

namespace BatchPeak
{

const char *to_str( const NotFitPeakMdaMethod method )
{
  switch( method )
  {
    case NotFitPeakMdaMethod::None:           return "none";
    case NotFitPeakMdaMethod::Currie:         return "currie";
    case NotFitPeakMdaMethod::CurrieAndDecon: return "currie-and-decon";
  }//switch( method )

  assert( 0 );
  throw runtime_error( "Invalid NotFitPeakMdaMethod value." );
}//const char *to_str( const NotFitPeakMdaMethod )


NotFitPeakMdaMethod not_fit_peak_mda_method_from_str( const string &str )
{
  string value = str;
  SpecUtils::trim( value );
  SpecUtils::to_lower_ascii( value );

  if( value == "none" )
    return NotFitPeakMdaMethod::None;

  if( value == "currie" )
    return NotFitPeakMdaMethod::Currie;

  if( (value == "currie-and-decon") || (value == "currie-and-deconvolution") )
    return NotFitPeakMdaMethod::CurrieAndDecon;

  throw runtime_error( "Invalid not-fit peak MDA method '" + str + "'; must be one of 'none',"
                       " 'currie', or 'currie-and-decon'." );
}//NotFitPeakMdaMethod not_fit_peak_mda_method_from_str( const string & )


void update_description( NotFitPeakMda &mda )
{
  mda.description = mda.result_summary;

  if( !mda.activity_summary.empty() )
    mda.description += "  " + mda.activity_summary;

  if( !mda.caveats.empty() )
    mda.description += "  " + mda.caveats;
}//void update_description( NotFitPeakMda & )


double gaussian_fraction_in_roi( const double num_fwhm )
{
  if( num_fwhm <= 0.0 )
    return 0.0;

  // Half-width, in units of sigma, is 0.5*num_fwhm*2.35482; the fraction within +-x sigma of the
  //  mean is erf( x / sqrt(2) ).
  const double half_width_sigma = 0.5 * num_fwhm * 2.35482;

  return erf( half_width_sigma / sqrt(2.0) );
}//double gaussian_fraction_in_roi( const double )


string confidence_level_str( const double confidence_level )
{
  char buffer[64] = { '\0' };

  if( confidence_level < 0.999 )
    snprintf( buffer, sizeof(buffer), "%.4g%%", 100.0*confidence_level );
  else
    snprintf( buffer, sizeof(buffer), "1-%.2G", (1.0 - confidence_level) );

  return buffer;
}//string confidence_level_str( const double )


void compute_decon_limit( NotFitPeakMda &mda,
                          const DetectionLimitCalc::DeconComputeInput &input,
                          const double gammas_per_bq,
                          const bool quantity_is_counts,
                          const double confidence_level,
                          const bool use_curie )
{
  mda.decon_attempted = true;
  mda.decon_quantity_is_counts = quantity_is_counts;
  mda.decon_computed = false;
  mda.decon_result = nullptr;

  if( !mda.currie_computed )
  {
    mda.decon_error = "the Currie-style limit, used to seed the search range, was not computed";
    return;
  }

  if( (gammas_per_bq <= 0.0) || IsNan(gammas_per_bq) || IsInf(gammas_per_bq) )
  {
    mda.decon_error = "invalid number of counts per unit activity";
    return;
  }

  const DetectionLimitCalc::CurrieMdaResult &res = mda.currie_result;

  // Exaggerate the range implied by the Currie-style limit, so we can be confident of bracketing
  //  the answer; the multiple is arbitrary, and matches what the GUI "Simple MDA" tool uses.
  const double diff_multiple = 50.0;
  const double smallest_range = 1.0 / gammas_per_bq;
  double min_quantity = 0.0, max_quantity = 0.0;

  if( res.source_counts > res.decision_threshold )
  {
    const double nominal = res.source_counts / gammas_per_bq;
    const double lower_diff = fabs( nominal - (res.lower_limit / gammas_per_bq) );
    const double upper_diff = fabs( (res.upper_limit / gammas_per_bq) - nominal );

    min_quantity = std::max( 0.0, nominal - diff_multiple*lower_diff );
    max_quantity = std::max( smallest_range, nominal + diff_multiple*upper_diff );
  }else if( res.upper_limit < 0.0f )
  {
    // Many fewer counts in the peak region than the sides predict; fall back to the Poisson
    //  uncertainty of the peak region to set the scale.
    const double poisson_uncert = sqrt( std::max( 0.0f, res.peak_region_counts_sum ) );
    min_quantity = 0.0;
    max_quantity = std::max( smallest_range, diff_multiple*poisson_uncert/gammas_per_bq );
  }else
  {
    min_quantity = 0.0;
    max_quantity = std::max( smallest_range, diff_multiple*res.upper_limit/gammas_per_bq );
  }//if( detected ) / else if( deficit ) / else

  if( (max_quantity <= min_quantity) || IsNan(max_quantity) || IsInf(max_quantity) )
  {
    mda.decon_error = "could not determine a valid range of values to search over";
    return;
  }

  try
  {
    const auto base_input = make_shared<const DetectionLimitCalc::DeconComputeInput>( input );
    const DetectionLimitCalc::DeconActivityOrDistanceLimitResult result
        = DetectionLimitCalc::get_activity_or_distance_limits( static_cast<float>(confidence_level),
                                                              base_input, false,
                                                              min_quantity, max_quantity, use_curie );

    mda.decon_result = make_shared<const DetectionLimitCalc::DeconActivityOrDistanceLimitResult>( result );
    mda.decon_computed = true;

    // The two methods look at the same data, so a large disagreement says something about the
    //  data, rather than about the methods - most often a curved continuum, or an interfering
    //  peak in the region.  Worth telling the user, rather than quietly reporting two numbers.
    const double currie_limit = res.upper_limit / gammas_per_bq;
    if( result.foundUpperCl && (currie_limit > 0.0) && (result.upperLimit > 0.0) )
    {
      mda.decon_over_currie_ratio = result.upperLimit / currie_limit;
      mda.methods_disagree = ((mda.decon_over_currie_ratio > 2.0)
                              || (mda.decon_over_currie_ratio < 0.5));

      if( mda.methods_disagree )
      {
        const double factor = (mda.decon_over_currie_ratio > 1.0)
                                ? mda.decon_over_currie_ratio : (1.0/mda.decon_over_currie_ratio);
        char buffer[256] = { '\0' };
        snprintf( buffer, sizeof(buffer),
                 "Note: the gross-counts and deconvolution limits differ by a factor of %.1f,"
                 " which usually means the continuum under the peak is not straight, or another"
                 " peak falls within the region.", factor );

        mda.caveats += string(mda.caveats.empty() ? "" : "  ") + buffer;
        update_description( mda );
      }//if( mda.methods_disagree )
    }//if( both limits are usable )
  }catch( std::exception &e )
  {
    mda.decon_error = e.what();
  }//try / catch
}//void compute_decon_limit(...)


void add_counts_decon_limits( vector<NotFitPeakMda> &mdas,
                              const shared_ptr<const SpecUtils::Measurement> &spectrum,
                              const shared_ptr<const DetectorPeakResponse> &drf,
                              const BatchPeakFitOptions &options )
{
  // Anything that already has a limit (e.g. the activity-space one an activity/shielding fit
  //  computes) is left alone.
  vector<shared_ptr<const PeakDef>> peaks;
  for( const NotFitPeakMda &mda : mdas )
  {
    if( !mda.decon_attempted && mda.exemplar_peak )
      peaks.push_back( mda.exemplar_peak );
  }

  if( peaks.empty() || !spectrum )
    return;

  const size_t num_side_channels = std::max( size_t(1), options.mda_num_side_channels );

  // The unity-efficiency detector response makes the "activity" the calculation limits be simply
  //  the number of peak counts.
  const shared_ptr<const DetectorPeakResponse> unity_drf
                                     = make_unity_fixed_geom_drf( peaks, spectrum, drf );

  for( NotFitPeakMda &mda : mdas )
  {
    if( mda.decon_attempted || !mda.exemplar_peak )
      continue;

    mda.decon_attempted = true;
    mda.decon_quantity_is_counts = true;

    const double fwhm = exemplar_peak_fwhm( *mda.exemplar_peak );

    if( !unity_drf )
    {
      mda.decon_error = "could not determine peak widths to use";
      continue;
    }

    if( fwhm <= 0.0 )
    {
      mda.decon_error = "exemplar peak has no width";
      continue;
    }

    if( !mda.currie_computed )
    {
      mda.decon_error = "the peak region could not be determined";
      continue;
    }

    DetectionLimitCalc::DeconRoiInfo roi;
    roi.roi_start = mda.currie_result.input.roi_lower_energy;
    roi.roi_end = mda.currie_result.input.roi_upper_energy;
    roi.continuum_type = PeakContinuum::OffsetType::Linear;
    roi.cont_norm_method = DetectionLimitCalc::DeconContinuumNorm::Floating;
    roi.num_lower_side_channels = num_side_channels;
    roi.num_upper_side_channels = num_side_channels;

    DetectionLimitCalc::DeconRoiInfo::PeakInfo peak_info;
    peak_info.energy = mda.currie_result.input.gamma_energy;
    peak_info.fwhm = static_cast<float>( fwhm );
    // With a unity-efficiency fixed-geometry DRF, a "counts per bq" of one makes the quantity
    //  being limited be the number of peak counts.
    peak_info.counts_per_bq_into_4pi = 1.0;
    roi.peak_infos.push_back( peak_info );

    DetectionLimitCalc::DeconComputeInput decon_input;
    decon_input.distance = 0.0;
    decon_input.activity = 0.0;
    decon_input.include_air_attenuation = false;
    decon_input.shielding_thickness = 0.0;
    decon_input.drf = unity_drf;
    decon_input.measurement = spectrum;
    decon_input.roi_info.push_back( roi );

    compute_decon_limit( mda, decon_input, 1.0, true, options.mda_confidence_level, true );
  }//for( NotFitPeakMda &mda : mdas )
}//void add_counts_decon_limits(...)


vector<NotFitPeakMda> compute_not_fit_peak_mdas(
                          const vector<shared_ptr<const PeakDef>> &unfit_exemplar_peaks,
                          const deque<shared_ptr<const PeakDef>> &fit_peaks,
                          const shared_ptr<const SpecUtils::Measurement> &spectrum,
                          const shared_ptr<const DetectorPeakResponse> &exemplar_drf,
                          const BatchPeakFitOptions &options )
{
  vector<NotFitPeakMda> answer;

  if( unfit_exemplar_peaks.empty()
     || (options.not_fit_peak_mda == NotFitPeakMdaMethod::None)
     || !spectrum
     || (spectrum->num_gamma_channels() < 4) )
  {
    return answer;
  }

  const string cl_str = confidence_level_str( options.mda_confidence_level );
  const bool want_decon = (options.not_fit_peak_mda == NotFitPeakMdaMethod::CurrieAndDecon);

  // A side-channel count of zero puts `currie_mda_calc` into its "the spectrum is asserted to be
  //  background" mode, which is a different calculation than we are documenting here.
  const size_t num_side_channels = std::max( size_t(1), options.mda_num_side_channels );

  for( const shared_ptr<const PeakDef> &peak : unfit_exemplar_peaks )
  {
    NotFitPeakMda mda;
    mda.exemplar_peak = peak;

    // Record the options used up front, so they are still reported for peaks whose limit cant
    //  be computed.
    mda.currie_result.input.detection_probability = options.mda_confidence_level;
    mda.currie_result.input.num_lower_side_channels = num_side_channels;
    mda.currie_result.input.num_upper_side_channels = num_side_channels;
    mda.signal_fraction_in_roi = gaussian_fraction_in_roi( options.mda_roi_num_fwhm );

    if( !peak )
    {
      assert( 0 );
      mda.currie_error = "invalid peak";
      mda.short_description = "Not computed";
      mda.result_summary = "Detection limit could not be computed: invalid peak.";
      update_description( mda );
      answer.push_back( mda );
      continue;
    }//if( !peak )

    const double mean = peak->mean();
    const double fwhm = exemplar_peak_fwhm( *peak );

    DetectionLimitCalc::CurrieMdaInput input;
    input.spectrum = spectrum;
    input.gamma_energy = static_cast<float>( mean );

    if( fwhm > 0.0 )
    {
      // The default region width of 2.5 FWHM (i.e., +-1.25 FWHM) is what ISO 11929:2010
      //  recommends.  For peaks that arent Gaussian, `exemplar_peak_fwhm(...)` gives back a width
      //  such that the default reproduces the peaks own ROI width.
      const double half_width = 0.5 * options.mda_roi_num_fwhm * fwhm;
      input.roi_lower_energy = static_cast<float>( mean - half_width );
      input.roi_upper_energy = static_cast<float>( mean + half_width );
    }else
    {
      input.roi_lower_energy = static_cast<float>( peak->lowerX() );
      input.roi_upper_energy = static_cast<float>( peak->upperX() );
    }

    input.num_lower_side_channels = num_side_channels;
    input.num_upper_side_channels = num_side_channels;
    input.detection_probability = options.mda_confidence_level;
    input.additional_uncertainty = 0.0f;

    try
    {
      mda.currie_result = DetectionLimitCalc::currie_mda_calc( input );
      mda.currie_computed = true;
    }catch( std::exception &e )
    {
      // Keep the input around, so the options actually used are still reported for this peak
      mda.currie_result.input = input;
      mda.currie_error = e.what();
      mda.result_type = NotFitPeakMda::MdaResultType::Error;
      mda.short_description = "Not computed";
      mda.result_summary = "Detection limit could not be computed: " + mda.currie_error + ".";
      update_description( mda );
      answer.push_back( mda );
      continue;
    }//try / catch

    const DetectionLimitCalc::CurrieMdaResult &res = mda.currie_result;

    // With (essentially) no counts anywhere in the region, the decision threshold and the upper
    //  limit both collapse to zero - the Gaussian statistics the calculation uses have broken
    //  down.  Left alone, that makes a single stray count read as a detection, and reports an
    //  upper limit of "less than 0 counts".  Ld does not degenerate (it tends to k^2, which is
    //  the right order for a zero-background Poisson limit), so fall back to quoting just that.
    mda.region_is_empty = (res.decision_threshold <= 0.0f);

    // Classify the result, using the same criteria as the GUI detection limit tools.
    if( mda.region_is_empty )
      mda.result_type = NotFitPeakMda::MdaResultType::NotDetected;
    else if( res.source_counts > res.decision_threshold )
      mda.result_type = NotFitPeakMda::MdaResultType::Detected;
    else if( res.upper_limit < 0.0f )
      mda.result_type = NotFitPeakMda::MdaResultType::Deficit;
    else
      mda.result_type = NotFitPeakMda::MdaResultType::NotDetected;

    // Flag limits whose continuum estimate may be contaminated by a neighboring peak.
    const pair<double,double> energy_range = currie_energy_range( res );
    for( const shared_ptr<const PeakDef> &p : fit_peaks )
    {
      if( p && (p->upperX() > energy_range.first) && (p->lowerX() < energy_range.second) )
        mda.overlaps_fit_peak = true;
    }

    for( const shared_ptr<const PeakDef> &p : unfit_exemplar_peaks )
    {
      if( p && (p != peak) && (p->mean() > energy_range.first) && (p->mean() < energy_range.second) )
        mda.overlaps_other_unfit_peak = true;
    }

    switch( mda.result_type )
    {
      case NotFitPeakMda::MdaResultType::NotDetected:
        if( mda.region_is_empty )
        {
          mda.short_description = "No counts in region";
          mda.result_summary = "Not detected; the peak region contains essentially no counts."
                               "  Minimum reliably detectable signal (Ld) is "
                               + SpecUtils::printCompact(res.detection_limit, 4) + " counts.";
          break;
        }//if( mda.region_is_empty )

        mda.short_description = "Less than Lc";
        mda.result_summary = "Not detected.  Observed signal is below the decision threshold"
                          " (Lc = " + SpecUtils::printCompact(res.decision_threshold, 4) + " counts);"
                          " signal is less than " + SpecUtils::printCompact(res.upper_limit, 4)
                          + " counts at the " + cl_str + " confidence level."
                          "  Minimum reliably detectable signal (Ld) is "
                          + SpecUtils::printCompact(res.detection_limit, 4) + " counts.";
        break;

      case NotFitPeakMda::MdaResultType::Detected:
        mda.short_description = "Greater than Lc";
        mda.result_summary = "Signal present but peak not fit: observed "
                          + SpecUtils::printCompact(res.source_counts, 4) + " counts is above the"
                          " decision threshold (Lc = "
                          + SpecUtils::printCompact(res.decision_threshold, 4) + " counts); signal"
                          " is between " + SpecUtils::printCompact(res.lower_limit, 4) + " and "
                          + SpecUtils::printCompact(res.upper_limit, 4) + " counts at the "
                          + cl_str + " confidence level.";
        break;

      case NotFitPeakMda::MdaResultType::Deficit:
        mda.short_description = "Fewer counts than expected";
        mda.result_summary = "Fewer counts than expected from the neighboring continuum; observed"
                          " signal is consistent with less than 0 counts at the " + cl_str
                          + " confidence level.";
        break;

      case NotFitPeakMda::MdaResultType::Error:
        assert( 0 );
        break;
    }//switch( mda.result_type )

    answer.push_back( mda );
  }//for( const shared_ptr<const PeakDef> &peak : unfit_exemplar_peaks )

  // The caveats apply to both limits, and go at the end of the description, after any activity
  //  information the activity/shielding fit may fill in later.
  for( NotFitPeakMda &mda : answer )
  {
    if( mda.overlaps_fit_peak )
      mda.caveats += string(mda.caveats.empty() ? "" : "  ")
                     + "Note: a fit peak overlaps the region used to estimate the continuum, so"
                       " this limit may be unreliable.";

    if( mda.overlaps_other_unfit_peak )
      mda.caveats += string(mda.caveats.empty() ? "" : "  ")
                     + "Note: another exemplar peak that was not fit lies within the evaluated"
                       " region.";

    // The calculation assumes the whole peak is inside the peak region, so a narrow region gives
    //  a limit on the counts within it, rather than on the peaks total area.  Negligible at the
    //  default region width, but worth saying once it is not.
    if( mda.currie_computed && (mda.signal_fraction_in_roi < 0.99) )
    {
      char buffer[256] = { '\0' };
      snprintf( buffer, sizeof(buffer),
               "Note: the peak region holds only %.0f%% of the peak, and the limit is not"
               " corrected for the signal outside it, so it is optimistic by about %.0f%%.",
               100.0*mda.signal_fraction_in_roi,
               100.0*(1.0/mda.signal_fraction_in_roi - 1.0) );
      mda.caveats += string(mda.caveats.empty() ? "" : "  ") + buffer;
    }//if( a meaningful part of the peak is outside the region )

    update_description( mda );
  }//for( NotFitPeakMda &mda : answer )

  // Done after the caveats above, so any note the comparison of the two methods adds goes last
  if( want_decon )
    add_counts_decon_limits( answer, spectrum, exemplar_drf, options );

  return answer;
}//compute_not_fit_peak_mdas(...)


void propagate_energy_cal( const shared_ptr<const SpecUtils::EnergyCalibration> &energy_cal,
                           shared_ptr<SpecUtils::Measurement> &to_spectrum,
                           shared_ptr<SpecMeas> &to_specfile,
                           const set<int> &used_sample_nums )
{
  assert( to_spectrum );
  assert( energy_cal );
  assert( to_specfile );
  
  shared_ptr<const SpecUtils::EnergyCalibration> original_cal = to_spectrum->energy_calibration();
  assert( original_cal );
  if( !energy_cal )
    throw runtime_error( "Missing energy in from spectrum." );
  
  const size_t num_spec_chan = to_spectrum->num_gamma_channels();
  shared_ptr<const SpecUtils::EnergyCalibration> updated_cal;
  if( energy_cal == original_cal )
  {
    // We already have this energy cal, nothing to do here
    updated_cal = energy_cal;
  }else
  {
    if( energy_cal->num_channels() == num_spec_chan )
    {
      to_spectrum->set_energy_calibration( energy_cal );
      updated_cal = energy_cal;
    }else
    {
      auto new_cal = make_shared<SpecUtils::EnergyCalibration>();
      
      switch( energy_cal->type() )
      {
        case SpecUtils::EnergyCalType::Polynomial:
        case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
          new_cal->set_polynomial( num_spec_chan, energy_cal->coefficients(), energy_cal->deviation_pairs() );
          break;
          
        case SpecUtils::EnergyCalType::FullRangeFraction:
          new_cal->set_full_range_fraction( num_spec_chan, energy_cal->coefficients(), energy_cal->deviation_pairs() );
          break;
          
        case SpecUtils::EnergyCalType::LowerChannelEdge:
          if( num_spec_chan > energy_cal->num_channels() )
          {
            throw std::runtime_error( " its lower energy channel calibration, and exemplar has"
                                     " fewer channels." );
            new_cal.reset();
          }else
          {
            vector<float> channel_energies = *new_cal->channel_energies();
            channel_energies.resize( num_spec_chan + 1 );
            new_cal->set_lower_channel_energy( num_spec_chan, channel_energies );
          }
          break;
          
        case SpecUtils::EnergyCalType::InvalidEquationType:
          assert( 0 );
          break;
      }//switch( energy_cal->type() )
      
      updated_cal = new_cal;
      if( new_cal )
        to_spectrum->set_energy_calibration( new_cal );
    }//if( num channels match exemplar ) / else
  }//if( energy_cal == original_cal )
  
  // Update `to_specfile` as well, as we may write it back out as a N42
  if( updated_cal && to_specfile )
  {
    const vector<string> &det_names = to_specfile->detector_names();
    
    if( used_sample_nums.size() > 1 )
    {
      // Translate peaks from old energy, to new energy
      auto peaks = to_specfile->peaks( used_sample_nums );
      if( peaks && peaks->size() && (original_cal != updated_cal) )
      {
        const deque<shared_ptr<const PeakDef>> new_peaks
               = EnergyCal::translatePeaksForCalibrationChange( *peaks, original_cal, updated_cal );
        to_specfile->setPeaks( new_peaks, used_sample_nums );
      }
    }//if( used_sample_nums.size() > 1 )
    
    for( const int sample : used_sample_nums )
    {
      shared_ptr<const SpecUtils::EnergyCalibration> prev_cal;
      for( const string &det : det_names )
      {
        auto m = to_specfile->measurement( sample, det );
        if( m && (m->num_gamma_channels() == updated_cal->num_channels()) )
        {
          prev_cal = m->energy_calibration();
          if( prev_cal != updated_cal )
            to_specfile->set_energy_calibration( updated_cal, m );
        }
      }//for( const string &det : det_names )
      
      if( prev_cal && (used_sample_nums.size() > 1) && (prev_cal != updated_cal) )
      {
        // Translate peaks from old energy, to new energy, only if we havent already done it
        auto peaks = to_specfile->peaks( {sample} );
        if( peaks && peaks->size() )
        {
          const deque<shared_ptr<const PeakDef>> new_peaks
          = EnergyCal::translatePeaksForCalibrationChange( *peaks, prev_cal, updated_cal );
          to_specfile->setPeaks( new_peaks, {sample} );
        }//if( peaks && peaks->size() )
      }//if( prev_cal && (used_sample_nums.size() > 1) )
    }//for( const int sample : used_sample_nums )
  }//if( updated_cal )
}//propagate_energy_cal(...)
  
void fit_energy_cal_from_fit_peaks( shared_ptr<SpecUtils::Measurement> &raw, vector<PeakDef> peaks )
{
  if( !raw || raw->num_gamma_channels() < 16 )
    throw runtime_error( "update_gain_from_peak: invalid input spectrum" );
  
  vector<EnergyCal::RecalPeakInfo> peakinfos;
  
  shared_ptr<const SpecUtils::EnergyCalibration> orig_cal = raw->energy_calibration();
  assert( orig_cal && orig_cal->valid() && (orig_cal->coefficients().size() > 1) );
  
  const SpecUtils::EnergyCalType energy_cal_type = (orig_cal && orig_cal->valid()) 
                                                    ? orig_cal->type()
                                                    : SpecUtils::EnergyCalType::InvalidEquationType;
  switch( energy_cal_type )
  {
    case SpecUtils::EnergyCalType::Polynomial:
    case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
    case SpecUtils::EnergyCalType::FullRangeFraction:
      if( orig_cal->coefficients().size() < 2 )
        throw std::logic_error( "Somehow the energy calibration has less than two coefficients." );
      break;
      
    case SpecUtils::EnergyCalType::LowerChannelEdge:
    case SpecUtils::EnergyCalType::InvalidEquationType:
      assert( 0 );
      throw std::logic_error( "The original calibration must be either polynomial or full-range-fraction." );
      break;
  }//switch( exemplar_cal->type() )
  
  
  for( const PeakDef &peak : peaks )
  {
    if( !peak.hasSourceGammaAssigned() )
    {
      cerr << "Warning: peak at " << peak.mean() << " keV doesnt have a source assigned, so will"
      << " not be used for energy calibration" << endl;
      continue;
    }
    
    
    if( !peak.useForEnergyCalibration() )  //This defaults to true, so if this is false, then user selected it
      continue;
    
    EnergyCal::RecalPeakInfo recalpeak;
    recalpeak.peakMean = peak.mean();
    recalpeak.peakMeanUncert = peak.meanUncert();
    recalpeak.peakMeanBinNumber = orig_cal->channel_for_energy( peak.mean() );
    recalpeak.photopeakEnergy = peak.gammaParticleEnergy();
    
    peakinfos.push_back( recalpeak );
  }//for( const double peak_energy : peak_energies )
  
  std::sort( begin(peakinfos), end(peakinfos), 
    []( const EnergyCal::RecalPeakInfo &lhs, const EnergyCal::RecalPeakInfo &rhs ) -> bool {
    return lhs.peakMean < rhs.peakMean;
  } );
  
  if( peakinfos.empty() )
    throw runtime_error( "No peaks selected for use in energy calibration." );
  
  const size_t nchannels = raw->num_gamma_channels();
  
  const vector<float> &orig_coefs = orig_cal->coefficients();
  const vector<pair<float,float>> &dev_pairs = orig_cal->deviation_pairs();
  const size_t num_coefs = orig_coefs.size();
  assert( num_coefs >= 2 );
  
  vector<float> fit_coefs_uncert;
  vector<float> fit_coefs = orig_coefs;
  vector<bool> fitfor( num_coefs, false );
  
  // If we only have a couple peaks, then we cant fit for like 4 coefficients; we'll
  //  just hardcode a rough heuristic for what coefficients to fit for, because I think
  //  bothering the user with this level of detail is probably a bit too much.
  //
  //  But also note that we have already fit peaks, so we must be reasonably close
  //  right now anyway, so we can be a bit more liberal in terms of fitting more
  //  coefficients than a user might normally do during an interactive session.
  //
  //  We dont want to update both gain and offset from like the two Co60 peaks, so
  //  we'll count the effective number of peaks, requiring them to be separated a bit.
  //  We will require the separation to be max( 100, min(0.1*total-energy-range, 200)) keV,
  //  with the 0.1, 100 and 200, all being entirely arbitrary, but hopefully reasonable.
  size_t num_effective_peaks = 1;
  const double energy_range = (orig_cal->upper_energy() - orig_cal->lower_energy());
  const double sep_dist = std::max( 100.0, std::min(0.1*energy_range, 200.0) );
  size_t last_peak = 0;
  for( size_t peak_index = 1; peak_index < peakinfos.size(); ++peak_index )
  {
    const double delta_energy = peakinfos[peak_index].peakMean - peakinfos[last_peak].peakMean;
    assert( delta_energy >= 0.0 );
    if( delta_energy >= sep_dist )
    {
      num_effective_peaks += 1;
      last_peak = peak_index;
    }
  }//for( loop over peaks to count how many are separated by at least `sep_dist` )
  
  if( num_effective_peaks == 1 )
  {
    fitfor[1] = true; // Only fit gain
  }else
  {
    for( size_t index = 0; (index < num_effective_peaks) && (index < fitfor.size()); ++index )
      fitfor[index] = true;
    // TODO: we could consider not fitting for offset if `peakinfos[0].peakMean` is less than 
    //       ~200 keV or something, but for the moment we wont, because we know we are close
    //       in energy calibration, and we know the peaks the user is interested in, so overfitting
    //       isnt as large of a concern as not lining up the ROI edges as much
  }
  
  shared_ptr<SpecUtils::EnergyCalibration> updated_cal = make_shared<SpecUtils::EnergyCalibration>();
  
  switch( energy_cal_type )
  {
    case SpecUtils::EnergyCalType::Polynomial:
    case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      EnergyCal::fit_energy_cal_poly( peakinfos, fitfor, nchannels, dev_pairs, fit_coefs, fit_coefs_uncert );
      updated_cal->set_polynomial( nchannels, fit_coefs, dev_pairs );
      break;
      
    case SpecUtils::EnergyCalType::FullRangeFraction:
      EnergyCal::fit_energy_cal_frf( peakinfos, fitfor, nchannels, dev_pairs, fit_coefs, fit_coefs_uncert );
      updated_cal->set_full_range_fraction( nchannels, fit_coefs, dev_pairs );
      break;
      
    case SpecUtils::EnergyCalType::LowerChannelEdge:
    case SpecUtils::EnergyCalType::InvalidEquationType:
      assert( 0 );
      break;
  }//switch( exemplar_cal->type() )
  
  raw->set_energy_calibration( updated_cal );
  
#if( PERFORM_DEVELOPER_CHECKS )
  cout << "Updated energy calibration using ROIs from exemplar.\n\tCoefficients:\n";
  assert( fit_coefs.size() == orig_cal->coefficients().size() );
  for( size_t i = 0; i < fit_coefs.size(); ++i )
    cout << "\t\t" << std::setprecision(6) << std::setw(12) << orig_cal->coefficients().at(i)
    << "\t-->\t" << std::setprecision(6) << std::setw(12) << fit_coefs[i] << endl;
  
  cout << "This moved peak energies:" << endl;
  for( const auto &peak : peaks )
  {
    const double energy = peak.mean();
    const double orig = orig_cal->channel_for_energy( energy );
    const double now = updated_cal->channel_for_energy( energy );
    cout << "\t" << std::setprecision(6) << std::setw(12) << energy << " keV from channel "
    << std::setprecision(6) << std::setw(12) << orig << " to "
    << std::setprecision(6) << std::setw(12) << now << endl;
  }
  
  cout << endl << endl;
#endif
}//void fit_energy_cal_from_fit_peaks(...)

  
std::shared_ptr<SpecMeas> write_concatenated_n42( const std::vector<ConcatRecord> &records,
                                                  const BatchPeakFitOptions &options,
                                                  std::vector<std::string> &warnings )
{
  if( records.empty() || options.output_dir.empty() )
    return nullptr;

  try
  {
    auto concatenated_spec = make_shared<SpecMeas>();
    int current_sample = 0;

    // Add remarks to indicate when this file was created
    auto now = chrono::system_clock::now();
    auto time_t = chrono::system_clock::to_time_t(now);
    stringstream time_str;
    time_str << put_time(localtime(&time_t), "%Y-%m-%d %H:%M:%S");
    concatenated_spec->add_remark( "Concatenated N42 file created on " + time_str.str() );

    // Only input files that produced more than one analysis need their sample numbers spelled out
    //  in the record title; adding it unconditionally would change the output of runs that arent
    //  using multi-sample handling at all.
    map<string,size_t> analyses_per_input;
    for( const ConcatRecord &record : records )
      analyses_per_input[record.source_file_path] += 1;

    for( const ConcatRecord &record : records )
    {
      if( !record.spectrum )
        continue;

      // Create a copy of the spectrum measurement
      auto new_meas = make_shared<SpecUtils::Measurement>( *record.spectrum );
      current_sample += 1;
      new_meas->set_sample_number( current_sample );

      // Add source file name to the measurement title.  When one input file was split into several
      //  analyses, include the sample numbers so the records stay distinguishable.
      string source_filename = SpecUtils::filename( record.source_file_path );
      if( analyses_per_input[record.source_file_path] > 1 )
      {
        string samples;
        for( const int sample : record.sample_numbers )
          samples += (samples.empty() ? "" : ",") + std::to_string( sample );
        if( !samples.empty() )
          source_filename += " (sample" + string(record.sample_numbers.size() > 1 ? "s " : " ") + samples + ")";
      }//if( this input file produced more than one analysis )

      // Note the summed spectrum often already carries the source file name as its title (see
      //  `SpecFile::sum_measurements`), so only append it when it actually adds something.
      const string original_title = new_meas->title();
      if( original_title.empty() || (original_title == SpecUtils::filename(record.source_file_path)) )
        new_meas->set_title( source_filename );
      else
        new_meas->set_title( source_filename + " - " + original_title );

      concatenated_spec->add_measurement( new_meas, false );

      if( !record.peaks.empty() )
      {
        const set<int> sample_nums = { current_sample };
        deque<shared_ptr<const PeakDef>> peaks_copy;
        for( const shared_ptr<const PeakDef> &peak : record.peaks )
          peaks_copy.push_back( peak );
        concatenated_spec->setPeaks( peaks_copy, sample_nums );
      }//if( !record.peaks.empty() )
    }//for( const ConcatRecord &record : records )

    concatenated_spec->cleanup_after_load( SpecUtils::SpecFile::ReorderSamplesByTime );

    const string concatenated_file = SpecUtils::append_path( options.output_dir, "concatenated.n42" );

    if( SpecUtils::is_file( concatenated_file ) && !options.overwrite_output_files )
    {
      warnings.push_back( "Not writing '" + concatenated_file + "', as it would overwrite a file."
                         " See the '--overwrite-output-files' option to force writing." );
    }else
    {
      if( !concatenated_spec->save2012N42File( concatenated_file ) )
        warnings.push_back( "Failed to write concatenated N42 file '" + concatenated_file + "'." );
      else
        cout << "Have written concatenated N42 file '" << concatenated_file << "'" << endl;
    }//if( would overwrite ) / else

    return concatenated_spec;
  }catch( const std::exception &e )
  {
    warnings.push_back( "Error creating concatenated N42 file: " + string(e.what()) );
  }//try / catch

  return nullptr;
}//std::shared_ptr<SpecMeas> write_concatenated_n42(...)


void fit_peaks_in_files( const std::string &exemplar_filename,
                        std::shared_ptr<const SpecMeas> cached_exemplar_n42,
                          const std::set<int> &exemplar_sample_nums,
                          const std::vector<std::string> &files,
                        std::vector<std::shared_ptr<SpecMeas>> cached_files,
                        const ::BatchPeak::BatchPeakFitOptions &options,
                        BatchPeakFitSummary * const results )
{
  vector<string> warnings;
  
  if( files.empty() )
    throw runtime_error( "No input files specified." );
  
  if( !options.output_dir.empty() && !SpecUtils::is_directory(options.output_dir) )
    throw runtime_error( "Output directory ('" + options.output_dir + "'), is not a directory." );
    
  if( options.write_n42_with_results && options.output_dir.empty() )
    throw runtime_error( "If you specify to write N42 files with the fit peaks, you must specify an output directory." );

  if( !cached_files.empty() && (cached_files.size() != files.size()) )
    throw runtime_error( "If you provide any cached intput files, the input vector of SpecMeas objects"
                        " must be the same size as input filenames" );

  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  if( !db )
  {
    const char *message = "Unable to load the nuclide decay database."
    " Run the executable from a working directory where the 'data' folder is "
    " located - the program expects to load data/sandia.decay.xml.";
    
    throw runtime_error( message );
  }//if( !db )
  
  const vector<pair<string,string>> spec_chart_js_and_css = BatchInfoLog::load_spectrum_chart_js_and_css();
  
  // Load report templates, and setup inja::environment
  const string tmplt_dir = BatchInfoLog::template_include_dir( options );
  inja::Environment env = BatchInfoLog::get_default_inja_env( tmplt_dir );
  
  
  nlohmann::json summary_json;
  
  BatchInfoLog::add_exe_info_to_json( summary_json );
  BatchInfoLog::add_peak_fit_options_to_json( summary_json, options );
  
  summary_json["ExemplarFile"] = exemplar_filename;
  if( !exemplar_sample_nums.empty() )
    summary_json["ExemplarSampleNumbers"] = vector<int>{begin(exemplar_sample_nums), end(exemplar_sample_nums)};
  summary_json["InputFiles"] = files;
  for( const pair<string,string> &key_val : spec_chart_js_and_css )
    summary_json[key_val.first] = key_val.second;

  // Records for the optional concatenated N42; accumulated as we go, rather than pulled out of
  //  `results` afterwards, so `--concatenate-to-n42` works even when the caller doesnt ask for the
  //  full per-file results to be retained (which the command line doesnt).
  vector<ConcatRecord> concat_records;

  // Each input file becomes one or more work items; more than one when the user has asked for a
  //  file holding several foreground records to be analyzed record-by-record.  The expansion is
  //  done a file at a time so that only one input file is held parsed in memory at once.
  size_t file_index = 0;  //index of the analysis, which is not the input file index when splitting

  for( size_t input_index = 0; input_index < files.size(); input_index += 1 )
  {
   const shared_ptr<SpecMeas> input_cached_file
                      = (input_index < cached_files.size()) ? cached_files[input_index] : nullptr;
   vector<BatchSampleSelect::InputWorkItem> work_items
      = BatchSampleSelect::expand_input_file( files[input_index], input_index, input_cached_file,
                                              options.multi_sample_handling );

   for( size_t item_index = 0; item_index < work_items.size(); item_index += 1, file_index += 1 )
   {
    BatchSampleSelect::InputWorkItem &item = work_items[item_index];
    const string &filename = item.filename;

    if( results )
    {
      results->file_results.resize( file_index + 1 );
      results->file_json.resize( file_index + 1 );
      results->file_peak_csvs.resize( file_index + 1 );
      results->file_reports.resize( file_index + 1 );
    }//if( results )

    // `fit_peaks_in_file` modifies the file handed to it, so work items that share a parsed file
    //  each need their own copy.  Only one copy is alive at a time.
    std::shared_ptr<SpecMeas> cached_file = item.source;
    if( item.needs_private_copy && item.source )
    {
      cached_file = make_shared<SpecMeas>();
      cached_file->uniqueCopyContents( *item.source );
    }

    const BatchPeak::BatchPeakFitResult fit_results
                 = fit_peaks_in_file( exemplar_filename, exemplar_sample_nums,
                                     cached_exemplar_n42, filename, cached_file,
                                     item.foreground_sample_numbers, options );

    // Release our reference to the parsed input file now that it has been analyzed, so that a run
    //  over many files doesnt hold every one of them in memory at once.
    item.source.reset();
    cached_file.reset();

    if( !cached_exemplar_n42 )
      cached_exemplar_n42 = fit_results.exemplar;

    // Only label the warnings when the input file was split or summed, so that output for the
    //  default handling stays exactly as it was.
    const bool label_warnings = (item.output_base_name != SpecUtils::filename(filename));
    for( const string &warn : fit_results.warnings )
      warnings.push_back( label_warnings ? ("File '" + item.label + "': " + warn) : warn );

    nlohmann::json data;
    BatchInfoLog::add_exe_info_to_json( data );
    BatchInfoLog::add_peak_fit_options_to_json( data, options );
    data["ExemplarFile"] = exemplar_filename;
    if( !exemplar_sample_nums.empty() )
      data["ExemplarSampleNumbers"] = vector<int>{begin(exemplar_sample_nums), end(exemplar_sample_nums)};
    data["Filepath"] = filename;
    // `Filename` identifies this analysis, and matches the output files written for it; for a file
    //  split into per-sample analyses it carries the same "_sampleN" infix the output files do.
    //  `SourceFilename` is always the unmodified input file leaf name.
    data["Filename"] = item.output_base_name;
    data["SourceFilename"] = SpecUtils::filename( filename );
    data["AnalysisLabel"] = item.label;
    data["IsSplitFromMultiSampleFile"] = (work_items.size() > 1);
    // Report the samples actually used, which are known even in `Auto` mode
    if( !fit_results.sample_numbers.empty() )
      data["ForegroundSampleNumbers"] = vector<int>{ begin(fit_results.sample_numbers),
                                                     end(fit_results.sample_numbers) };
    data["ParentDir"] = SpecUtils::parent_path( filename );

    BatchInfoLog::add_peak_fit_results_to_json( data, fit_results );
    
    summary_json["Files"].push_back( data );
    
    for( const pair<string,string> &key_val : spec_chart_js_and_css )
      data[key_val.first] = key_val.second;

    if( results )
    {
      results->file_json[file_index] = data.dump(2);
      results->file_results[file_index] = fit_results;
      results->file_reports[file_index].resize( options.report_templates.size() );
    }

    for( size_t tmplt_index = 0; tmplt_index < options.report_templates.size(); ++tmplt_index )
    {
      const std::string tmplt_path = get_report_template_path( options.report_templates[tmplt_index] );
      const std::string tmplt_name = get_report_template_name( options.report_templates[tmplt_index] );

      try
      {
        const string rpt = BatchInfoLog::render_template( tmplt_path, env,
                            BatchInfoLog::TemplateRenderType::PeakFitIndividual, options, data );

        if( results )
          results->file_reports[file_index][tmplt_index] = rpt;

        if( options.to_stdout && !SpecUtils::iequals_ascii(tmplt_path, "html" ) )
          cout << "\n\n" << rpt << endl << endl;
        
        if( !options.output_dir.empty() )
        {
          const string out_file
                    = BatchInfoLog::suggested_output_report_filename( item.output_base_name, tmplt_name,
                                  BatchInfoLog::TemplateRenderType::PeakFitIndividual, options );

          if( SpecUtils::is_file(out_file) && !options.overwrite_output_files )
          {
            warnings.push_back( "Not writing '" + out_file + "', as it would overwrite a file."
                               " See the '--overwrite-output-files' option to force writing." );
          }else
          {
#ifdef _WIN32
            const std::wstring woutcsv = SpecUtils::convert_from_utf8_to_utf16(out_file);
            std::ofstream output( woutcsv.c_str(), ios::binary | ios::out );
#else
            std::ofstream output( out_file.c_str(), ios::binary | ios::out);
#endif
            if( !output )
              warnings.push_back( "Failed to open report output '" + out_file + "', for writing.");
            else
              output.write( rpt.c_str(), rpt.size() );
          }//if( is file ) / else write file
        }//if( !options.output_dir.empty() )
      }catch( inja::InjaError &e )
      {
        const string msg = "Error templating results (" + e.type + ": line "
        + std::to_string(e.location.line) + ", column " + std::to_string(e.location.column)
        + " of '" + tmplt_path + "'): " + e.message + ". While processing '" + filename + "'.";
        
        cerr << msg << endl;
        warnings.push_back( msg );
      }catch( std::exception &e )
      {
        cerr << "Error templating results: " << e.what() << endl;
        warnings.push_back( "Error templating results: " + string(e.what()) );
      }
    }//for( const string &tmplt : options.report_templates )
    
    if( !options.output_dir.empty() && options.create_json_output )
      BatchInfoLog::write_json( options, warnings, item.output_base_name, data );
    
    if( !fit_results.success )
      continue;
    
    assert( fit_results.measurement );
    
    if( options.write_n42_with_results && fit_results.measurement )
    {
      string outn42 = SpecUtils::append_path(options.output_dir, item.output_base_name );
      if( !SpecUtils::iequals_ascii(SpecUtils::file_extension(item.output_base_name), ".n42") )
        outn42 += ".n42";
      
      if( SpecUtils::is_file( outn42 ) && !options.overwrite_output_files )
      {
        warnings.push_back( "Not writing '" + outn42 + "', as it would overwrite a file."
                           " See the '--overwrite-output-files' option to force writing." );
      }else
      {
        if( !fit_results.measurement->save2012N42File( outn42 ) )
          warnings.push_back( "Failed to write '" + outn42 + "'.");
        else
          cout << "Have written '" << outn42 << "' with peaks" << endl;
      }
    }//if( options.write_n42_with_peaks )
    
    deque<shared_ptr<const PeakDef>> fit_peaks = fit_results.fit_peaks;
    if( options.show_nonfit_peaks )
    {
      for( const auto p : fit_results.unfit_exemplar_peaks )
      {
        auto np = make_shared<PeakDef>(*p);
        np->setAmplitude( 0.0 );
        np->setAmplitudeUncert( 0.0 );
        np->setSigmaUncert( 0.0 );
        auto cont = make_shared<PeakContinuum>( *np->continuum() );
        const size_t num_cont_pars = PeakContinuum::num_parameters(cont->type());
        for( size_t i = 0; i < num_cont_pars; ++i )
        {
          cont->setPolynomialCoef( i, 0.0 );
          cont->setPolynomialUncert( i, -1.0 );
        }
        np->setContinuum( cont );
        np->set_coefficient( -1.0, PeakDef::Chi2DOF );
        fit_peaks.push_back(np);
      }
      std::sort( begin(fit_peaks), end(fit_peaks), &PeakDef::lessThanByMeanShrdPtr );
    }//if( make_nonfit_peaks_zero )
    
    
    if( !options.output_dir.empty() && options.create_csv_output )
    {
      const string &leaf_name = item.output_base_name;
      string outcsv = SpecUtils::append_path(options.output_dir, leaf_name) + ".CSV";
      
      if( SpecUtils::is_file(outcsv) && !options.overwrite_output_files )
      {
        warnings.push_back( "Not writing '" + outcsv + "', as it would overwrite a file."
                           " See the '--overwrite-output-files' option to force writing." );
      }else
      {
#ifdef _WIN32
        const std::wstring woutcsv = SpecUtils::convert_from_utf8_to_utf16(outcsv);
        std::ofstream output_csv( woutcsv.c_str() );
#else
        std::ofstream output_csv( outcsv.c_str() );
#endif
        
        if( !output_csv )
        {
          warnings.push_back( "Failed to open '" + outcsv + "', for writing.");
        }else
        {
          PeakModel::write_peak_csv( output_csv, leaf_name, 
                                    PeakModel::PeakCsvType::Full, fit_peaks, fit_results.spectrum );
          cout << "Have written '" << outcsv << "'" << endl;
        }
      }//if( SpecUtils::is_file( outcsv ) ) / else
    }//if( !options.output_dir.empty() )

    if( results )
    {
      stringstream out_csv;
      PeakModel::write_peak_csv( out_csv, item.output_base_name,
                                PeakModel::PeakCsvType::Full, fit_peaks, fit_results.spectrum );
      results->file_peak_csvs[file_index] = out_csv.str();
    }

    if( options.concatenate_to_n42 && !options.output_dir.empty() && fit_results.spectrum )
    {
      ConcatRecord record;
      record.source_file_path = filename;
      record.sample_numbers = fit_results.sample_numbers;
      record.spectrum = fit_results.spectrum;
      record.peaks = fit_peaks;
      concat_records.push_back( record );
    }//if( options.concatenate_to_n42 ... )

    if( options.to_stdout )
    {
      const string &leaf_name = item.output_base_name;
      cout << "peaks for '" << item.label << "':" << endl;
      PeakModel::write_peak_csv( cout, leaf_name, PeakModel::PeakCsvType::Full,
                                fit_peaks, fit_results.spectrum );
      cout << endl;
    }
   }//for( loop over work items of this input file )
  }//for( loop over input files )

  // `Files` holds one entry per analysis performed; with multi-sample handling this may be more
  //  entries than there were input files.
  summary_json["NumInputFiles"] = static_cast<int>( files.size() );
  summary_json["NumAnalyses"] = static_cast<int>( file_index );

  // Add any encountered errors to output summary JSON
  for( const string &warn : warnings )
    summary_json["Warnings"].push_back( warn );
  
  // Now write summary report(s)
  for( const string &summary_tmplt : options.summary_report_templates )
  {
    const std::string tmplt_path = get_report_template_path( summary_tmplt );
    const std::string tmplt_name = get_report_template_name( summary_tmplt );

    try
    {
      const string rpt = BatchInfoLog::render_template( tmplt_path, env,
                       BatchInfoLog::TemplateRenderType::PeakFitSummary, options, summary_json );

      if( results )
        results->summary_reports.push_back( rpt );

      if( options.to_stdout && !SpecUtils::iequals_ascii(summary_tmplt, "html" ) )
        cout << "\n\n" << rpt << endl << endl;
      
      if( !options.output_dir.empty() )
      {
        const string out_file
                    = BatchInfoLog::suggested_output_report_filename( "", tmplt_name,
                                    BatchInfoLog::TemplateRenderType::PeakFitSummary, options );
        
        if( SpecUtils::is_file(out_file) && !options.overwrite_output_files )
        {
          warnings.push_back( "Not writing '" + out_file + "', as it would overwrite a file."
                             " See the '--overwrite-output-files' option to force writing." );
        }else
        {
#ifdef _WIN32
          const std::wstring woutcsv = SpecUtils::convert_from_utf8_to_utf16(out_file);
          std::ofstream output( woutcsv.c_str(), ios::binary | ios::out );
#else
          std::ofstream output( out_file.c_str(), ios::binary | ios::out );
#endif
          if( !output )
            throw runtime_error( "Failed to open summary peak fit report output, '" + out_file + "'" );
          
          output.write( rpt.c_str(), rpt.size() );
        }//if( file exists and dont overwrite ) / else
      }
    }catch( inja::InjaError &e )
    {
      const string msg = "Error templating summary peak fit output (" + e.type + ": line "
      + std::to_string(e.location.line) + ", column " + std::to_string(e.location.column)
      + " of '" + tmplt_path + "'): " + e.message + ".";
      
      cerr << msg << endl;
      warnings.push_back( msg );
    }catch( std::exception &e )
    {
      warnings.push_back( "Error making summary peak fit output: " + string(e.what()) );
    }
  }//if( !options.summary_report_template.empty() )
  
  
  if( !options.output_dir.empty() && options.create_json_output )
    BatchInfoLog::write_json( options, warnings, "", summary_json );
  
  // Create concatenated N42 file if requested
  if( options.concatenate_to_n42 && !options.output_dir.empty() )
  {
    shared_ptr<SpecMeas> concatenated = write_concatenated_n42( concat_records, options, warnings );
    if( results )
      results->concatenated_results = concatenated;
  }//if( options.concatenate_to_n42 && !options.output_dir.empty() )
  
  if( !warnings.empty() )
    cerr << endl << endl;
  for( const auto warn : warnings )
    cerr << warn << endl;

  if( results )
  {
    results->options = options;
    results->exemplar_filename = exemplar_filename;
    results->exemplar = cached_exemplar_n42;
    results->exemplar_sample_nums = exemplar_sample_nums;

    results->summary_json = summary_json.dump(2);
    results->warnings = warnings;
  }
}//fit_peaks_in_files(...)
  
  
void get_exemplar_spectrum_and_peaks(
            std::shared_ptr<const SpecUtils::Measurement> &exemplar_spectrum,
            std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> &exemplar_peaks,
            std::set<int> &exemplar_sample_nums,
            const std::shared_ptr<const SpecMeas> &exemplar_n42,
            const bool require_peaks )
{
  const vector<string> det_names = exemplar_n42->detector_names();
  const set<set<int>> withPeakSampleNums = exemplar_n42->sampleNumsWithPeaks();
    
  if( exemplar_n42->measurements().empty() )
  {
    throw logic_error( "Exemplar spectrum file did not have any measurements." );
  }else if( !exemplar_sample_nums.empty() )
  {
    const set<set<int>> withPeakSampleNums = exemplar_n42->sampleNumsWithPeaks();
    if( !withPeakSampleNums.count(exemplar_sample_nums) )
      throw runtime_error( "The specified exemplar sample numbers did not have peaks associated with them." );
      
    exemplar_peaks = exemplar_n42->peaks( exemplar_sample_nums );
    exemplar_spectrum = exemplar_n42->sum_measurements( exemplar_sample_nums, det_names, nullptr );
      
    if( !exemplar_peaks || (require_peaks && exemplar_peaks->empty()) || !exemplar_spectrum )
      throw runtime_error( "The specified exemplar sample numbers did not have peaks, or spectra couldnt be summed." );
  }else if( exemplar_n42->measurements().size() == 1 )
  {
    exemplar_spectrum = exemplar_n42->measurements().front();
    if( !exemplar_spectrum )
      throw logic_error( "Unexpected invalid exemplar spectrum." );
      
    exemplar_sample_nums.insert( exemplar_spectrum->sample_number() );
    exemplar_peaks = exemplar_n42->peaks( exemplar_sample_nums );
    if( require_peaks && (!exemplar_peaks || exemplar_peaks->empty()) )
      throw logic_error( "Exemplar spectrum did not contain any peaks." );
  }else
  {
    set<set<int>> foregroundPeaks, backgroundPeaks, otherPeaks;
    for( const set<int> &samples : withPeakSampleNums )
    {
      auto peaks = exemplar_n42->peaks( samples );
      if( !peaks || peaks->empty() )
        continue;
        
      auto m = exemplar_n42->sum_measurements( samples, det_names, nullptr );
      if( !m )
        continue;
      
      switch( m->source_type() )
      {
        case SpecUtils::SourceType::IntrinsicActivity:
        case SpecUtils::SourceType::Calibration:
          otherPeaks.insert( samples );
          break;
          
        case SpecUtils::SourceType::Background:
          backgroundPeaks.insert( samples );
          break;
          
        case SpecUtils::SourceType::Foreground:
        case SpecUtils::SourceType::Unknown:
          foregroundPeaks.insert( samples );
          break;
      }//switch( m->source_type() )
    }//for( const set<int> &samples : withPeakSampleNums )
    
    if( foregroundPeaks.size() > 1 )
      throw runtime_error( "Ambiguous which peaks to use from exemplar file" );
    
    if( foregroundPeaks.empty() && backgroundPeaks.empty() )
      throw runtime_error( "No valid peaks exemplar file"
                          + string(otherPeaks.empty() ? "." : " (intrinsic and calibration spectra in files are ignored).") );
    
    if( foregroundPeaks.empty() && (backgroundPeaks.size() != 1) )
      throw runtime_error( "Ambiguous which peaks to use from exemplar file; multiple background spectra with peaks." );
    
    if( foregroundPeaks.size() == 1 )
    {
      exemplar_sample_nums = *begin(foregroundPeaks);
    }else if( backgroundPeaks.size() == 1 )
    {
      exemplar_sample_nums = *begin(backgroundPeaks);
    }else
    {
      throw logic_error( "Error getting peaks from exemplar." );
    }
    
    exemplar_peaks = exemplar_n42->peaks( exemplar_sample_nums );
    exemplar_spectrum = exemplar_n42->sum_measurements( exemplar_sample_nums, det_names, nullptr );
  }//if( sample nums specified ) / else ( single meas ) / else ( search for peaks )

  if( !exemplar_peaks )
    return;

  //We need to make sure that the exemplar peaks have unique continuums
  //  This is because the peak fitter will modify the continuums, and we need to make sure that
  //  the exemplar peaks have the same continuums as the fit peaks
  shared_ptr<deque<shared_ptr<const PeakDef>>> exemplar_peaks_copy = make_shared<deque<shared_ptr<const PeakDef>>>();
  
  std::map<std::shared_ptr<const PeakContinuum>, std::shared_ptr<PeakContinuum>> cont_map;
  for( size_t peak_index = 0; peak_index < exemplar_peaks->size(); peak_index += 1 )
  {
    const shared_ptr<const PeakDef> &peak = exemplar_peaks->at( peak_index );
    assert( peak );
    const shared_ptr<const PeakContinuum> cont = peak->continuum();
    auto pos = cont_map.find( cont );
    if( pos == end(cont_map) )
      pos = cont_map.insert( { cont, std::make_shared<PeakContinuum>( *cont ) } ).first;
    shared_ptr<PeakDef> new_peak = make_shared<PeakDef>( *peak );
    new_peak->setContinuum( pos->second );
    
    exemplar_peaks_copy->push_back( new_peak );
  }//for( size_t peak_index = 0; peak_index < exemplar_peaks->size(); peak_index += 1 )

  exemplar_peaks = exemplar_peaks_copy;
}//void get_exemplar_spectrum_and_peaks(...)
  

BatchPeak::BatchPeakFitResult fit_peaks_in_file( const std::string &exemplar_filename,
                          std::set<int> exemplar_sample_nums,
                          std::shared_ptr<const SpecMeas> cached_exemplar_n42,
                          const std::string &filename,
                          std::shared_ptr<SpecMeas> cached_spectrum,
                          std::set<int> foreground_sample_numbers,
                          const BatchPeakFitOptions &options )
{
  shared_ptr<const SpecMeas> exemplar_n42 = cached_exemplar_n42;

  if( !exemplar_n42 && !exemplar_filename.empty() )
  {
    // Read in the N42 file exported from InterSpec.
    //  This file should have good energy calibration applied, have the peaks fit that you care
    //  about, and have exactly one spectrum
    auto exemplar = make_shared<SpecMeas>();
    const bool exemplar_is_n42 = exemplar->load_N42_file( exemplar_filename );
        
    if( !exemplar_is_n42 && (options.use_exemplar_energy_cal || options.use_exemplar_energy_cal_for_background) )
      throw runtime_error( "Exemplar file wasnt an N42 file, but using its energy cal was specified - not allowed." );

    if( exemplar_is_n42 )
      

    if( exemplar_is_n42 )
      exemplar_n42 = exemplar;
  }//if( !cached_exemplar_n42 )
  
  std::shared_ptr<const SpecUtils::Measurement> exemplar_spectrum;
  std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> exemplar_peaks;

  if( options.fit_all_peaks && !options.use_exemplar_energy_cal )
  {
    // Nothing to do here
  }else if( options.fit_all_peaks )
  {
    // We want to use the exemplar energy cal
    if( !exemplar_n42 )
      throw runtime_error( "Exemplar file not provided, but using its energy cal was requested." );

    get_exemplar_spectrum_and_peaks( exemplar_spectrum, exemplar_peaks,
                                    exemplar_sample_nums, exemplar_n42, false );

  }else if( exemplar_n42 )
  {
    get_exemplar_spectrum_and_peaks( exemplar_spectrum, exemplar_peaks,
                                    exemplar_sample_nums, exemplar_n42, true );

    assert( exemplar_peaks );
    assert( exemplar_spectrum );
    if( !exemplar_peaks || !exemplar_spectrum )
      throw logic_error( "Logic error retrieving spectrum from exemplar." );
    
    if( !exemplar_peaks || !exemplar_spectrum )
      throw logic_error( "Logic error retrieving peaks from exemplar." );
    
    if( (options.use_exemplar_energy_cal || options.use_exemplar_energy_cal_for_background)
       && (exemplar_spectrum->energy_calibration_model() == SpecUtils::EnergyCalType::InvalidEquationType) )
      throw runtime_error( "Exemplar spectrum doesnt have a valid energy calibration." );
    
    if( exemplar_spectrum->num_gamma_channels() < 64  )
      throw runtime_error( "Exemplar spectrum doesnt have enough gamma channels." );

    assert( !exemplar_n42 || exemplar_spectrum );
    assert( !exemplar_n42 || (exemplar_peaks && !exemplar_peaks->empty()) );
  }//if( exemplar_is_n42 )

  BatchPeakFitResult results;
  results.file_path = filename;
  results.options = options;
  results.exemplar = exemplar_n42;
  results.exemplar_sample_nums = exemplar_sample_nums;
  //if( exemplar_peaks )
  //  results.exemplar_peaks = *exemplar_peaks;
  results.exemplar_spectrum = exemplar_spectrum;
  results.success = false;
  
  {
    shared_ptr<SpecMeas> specfile;
    
    if( cached_spectrum )
    {
      specfile = cached_spectrum;
    }else
    {
      specfile = make_shared<SpecMeas>();
      const bool loaded = specfile->load_file( filename, SpecUtils::ParserType::Auto, filename );
      if( !loaded || !specfile->num_measurements() )
      {
        if( SpecUtils::is_file(filename) )
          results.warnings.push_back( "Couldnt read in '" + filename + "' as a spectrum file -- skipping." );
        else
          results.warnings.push_back( "Could not access '" + filename + "' -- skipping." );

        return results;
      }//if( !loaded )
    }//if( cached_spectrum ) / else
    
    results.measurement = specfile;
    
    shared_ptr<SpecUtils::Measurement> spec;
    const vector<string> det_names = specfile->detector_names();
    
    set<int> used_sample_nums;
    if( !foreground_sample_numbers.empty() )
    {
      try
      {
        spec = specfile->sum_measurements( foreground_sample_numbers, det_names, nullptr );
        used_sample_nums = foreground_sample_numbers;
      }catch( std::exception &e )
      {
        results.warnings.push_back( "Invalid sample numbers specified to sum: " + string(e.what()) );
        results.success = false;
        return results;
      }
    }else
    {
      if( specfile->sample_numbers().size() == 1 )
      {
        used_sample_nums = specfile->sample_numbers();
        spec = specfile->sum_measurements( used_sample_nums, det_names, nullptr );
      }else
      {
        // Note: this ladder is deliberately more permissive than the one `BatchActivity` and
        //  `BatchRelActAuto` use - Foreground and Unknown are kept separate, and a lone
        //  calibration/intrinsic spectrum is accepted, since fitting peaks in one is a sensible
        //  thing to do even though fitting an activity to one is not.
        const BatchSampleSelect::SampleBuckets buckets
                                                = BatchSampleSelect::classify_samples( *specfile );
        const set<int> &foregroundSamples = buckets.foreground;
        const set<int> &backgroundSamples = buckets.background;
        const set<int> &unknownSamples = buckets.unknown;
        const set<int> &otherSamples = buckets.other;

        if( foregroundSamples.size() == 1 )
        {
          used_sample_nums = foregroundSamples;
        }else if( unknownSamples.size() == 1 )
        {
          used_sample_nums = unknownSamples;
        }else if( backgroundSamples.size() == 1 )
        {
          used_sample_nums = backgroundSamples;
        }else if( otherSamples.size() == 1 )
        {
          used_sample_nums = otherSamples;
        }else
        {
          results.warnings.push_back( "Spectrum file '" + filename + "' was ambiguous of which spectrum to use for peak fitting." );
          return results;
        }
        
        spec = specfile->sum_measurements( used_sample_nums, det_names, nullptr );
      }//if( specfile.sample_numbers().size() == 1 ) / else
    }//if( !foreground_sample_numbers.empty() ) / else
    
    assert( spec );
    if( !spec )
    {
      results.warnings.push_back( "Spectrum file '" + filename + "' failed to extract wanted spectrum." );
      return results;
    }
    
    
    if( options.use_exemplar_energy_cal )
    {
      try
      {
        // We will make a copy of the energy calibration in case we modify it further in the future; we dont
        //  want to modify the exemplar energy calibration
        shared_ptr<const SpecUtils::EnergyCalibration> exemplar_cal = exemplar_spectrum->energy_calibration();
        if( !exemplar_cal || !exemplar_cal->valid() )
          throw runtime_error( "Exemplar spectrum doesnt have a valid energy calibration." );
        
        shared_ptr<SpecUtils::EnergyCalibration> spec_cal = make_shared<SpecUtils::EnergyCalibration>( *exemplar_cal );

        propagate_energy_cal( spec_cal, spec, specfile, used_sample_nums );
      }catch( std::exception &e )
      {
        results.warnings.push_back( "Not using exemplar energy calibration for '" + filename + "': "
                                   + std::string(e.what()) );
      }
    }//if( options.use_exemplar_energy_cal )
    
    
    if( !spec || (spec->num_gamma_channels() < 64)
       || (spec->energy_calibration_model() == SpecUtils::EnergyCalType::InvalidEquationType) )
    {
      results.warnings.push_back( "Failed to get spectrum from file '" + filename + "' -- skipping." );
      return results;
    }
    
    set<int> back_sample_nums;
    shared_ptr<SpecMeas> background_n42;
    if( !options.background_subtract_file.empty() || options.cached_background_subtract_spec )
    {
      if( options.cached_background_subtract_spec  )
      {
        background_n42 = options.cached_background_subtract_spec;
      }else
      {
        background_n42 = make_shared<SpecMeas>();
        if( !background_n42->load_file( options.background_subtract_file, SpecUtils::ParserType::Auto ) )
          throw runtime_error( "Couldnt open background file '" + options.background_subtract_file + "'" );
      }
      
      if( options.background_subtract_samples.empty() )
      {
        back_sample_nums = background_n42->sample_numbers();
        if( back_sample_nums.size() != 1 )
          throw runtime_error( "There should only be a single sample in background subtract file." );
      }else
      {
        back_sample_nums = options.background_subtract_samples;
      }
      
      if( background_n42->num_measurements() == 1 )
      {
        assert( !back_sample_nums.empty() );
        if( !back_sample_nums.empty()
           && ((*begin(back_sample_nums)) != background_n42->measurements()[0]->sample_number() ) )
        {
          results.warnings.push_back( "Specified background sample number invalid." );
          return results;
        }
        
        results.background = make_shared<SpecUtils::Measurement>( *(background_n42->measurements()[0]) );
      }else
      {
        try
        {
          const vector<string> &dets = background_n42->detector_names();
          results.background = background_n42->sum_measurements( back_sample_nums, dets, nullptr );
        }catch( std::exception &e )
        {
          results.warnings.push_back( "Failed to sum spectrum from background '"
                                     + options.background_subtract_file + "' -- skipping '"
                                     + filename + "'." );
          return results;
        }//try / catch to sum background
        
        if( !results.background->energy_calibration()
           || !results.background->energy_calibration()->valid()
           || (results.background->num_gamma_channels() < 16)
           || (results.background->live_time() <= 1.0E-3) )
        {
          results.warnings.push_back( "Background '"
                                     + options.background_subtract_file
                                     + "' didnt have energy calibration, too few channels, or no live-time"
                                     " -- skipping '" + filename + "'." );
          return results;
        }//
        
        if( options.use_exemplar_energy_cal_for_background )
        {
          try
          {
             shared_ptr<const SpecUtils::EnergyCalibration> exemplar_cal = exemplar_spectrum->energy_calibration();
            if( !exemplar_cal || !exemplar_cal->valid() )
              throw runtime_error( "Exemplar spectrum doesnt have a valid energy calibration." );        
            shared_ptr<SpecUtils::EnergyCalibration> spec_cal = make_shared<SpecUtils::EnergyCalibration>( *exemplar_cal );

            propagate_energy_cal( spec_cal, results.background, background_n42, back_sample_nums );
          }catch( std::exception &e )
          {
            results.warnings.push_back( "Not using exemplar energy calibration for background of '" + filename + "': "
                                       + std::string(e.what()) );
          }
        }//if( options.use_exemplar_energy_cal_for_background )
      }//if( background->num_measurements() == 1 ) / else
      
      
      try
      {
        const bool no_neg = true;
        const bool do_round = false;
        
        double sf = spec->live_time() / results.background->live_time();
        if( IsInf(sf) || IsNan(sf) )
          sf = 1.0; //e.g., if we dont have live-time

        shared_ptr<const vector<float>> fore_counts = spec->gamma_counts();
        shared_ptr<const vector<float>> back_counts = results.background->gamma_counts();
        
        // Make sure back_counts has the same energy calibration and fore_counts, so we can subtract
        //  on a bin-by-bin basis
        if( results.background->energy_calibration() != spec->energy_calibration()
           && (*results.background->energy_calibration()) != (*spec->energy_calibration()) )
        {
          auto new_backchan = make_shared<vector<float>>( fore_counts->size(), 0.0f );
          SpecUtils::rebin_by_lower_edge( *results.background->channel_energies(), *back_counts,
                                         *spec->channel_energies(), *new_backchan );
          back_counts = new_backchan;
        }
        
        // Create what will be the background subtracted foreground
        auto back_sub_counts = make_shared<vector<float>>( *fore_counts );
        
        //back_counts and fore_counts should always be the same size, but we'll be a bit verbose anyway
        assert( back_counts->size() == fore_counts->size() );
        const size_t nchann = std::min( back_counts->size(), fore_counts->size() );
        
        // Do the actual background subtraction
        for( size_t i = 0; i < nchann; ++i )
        {
          float &val = (*back_sub_counts)[i];
          val -= sf*(*back_counts)[i];
          
          if( no_neg )
            val = std::max( 0.0f, val );
          
          if( do_round )
            val = std::round( val );
        }//for( size_t i = 0; i < nchann; ++i )
        
        // Create a new Measurement object, based on the old foreground
        auto newspec = make_shared<SpecUtils::Measurement>( *spec );
        newspec->set_gamma_counts( back_sub_counts, spec->live_time(), spec->real_time() );
        vector<string> remarks = spec->remarks();
        remarks.push_back( "This spectrum has been background subtracted in InterSpec" );
        newspec->set_remarks( remarks );
        newspec->set_sample_number( 1 );
        
        spec = newspec;
        
        results.measurement->removeAllPeaks();
        results.measurement->remove_measurements( results.measurement->measurements() );
        results.measurement->add_measurement( spec, true );
      }catch( std::exception &e )
      {
        results.warnings.push_back( "Error background subtracting: '" + string(e.what())
                                   + "' -- skipping '" + filename + "'." );
        return results;
      }//try / catch
    }//if( !options.background_subtract_file.empty() )
    
    
    results.sample_numbers = used_sample_nums;
    results.spectrum = spec;

    deque<shared_ptr<const PeakDef>> fit_peaks_ptrs;
    set<shared_ptr<const PeakDef>> unused_exemplar_peaks;

    if( options.fit_all_peaks )
    {
      std::shared_ptr<const DetectorPeakResponse> det;
      if( exemplar_n42 )
        det = exemplar_n42->detector();

      // TODO: 'peak-stat-threshold' and 'peak-shape-threshold' are not currently taken into account.
      if( (options.peak_stat_threshold != 2.0) || (options.peak_hypothesis_threshold != 1.0) )
        results.warnings.push_back( "peak-stat-threshold and peak-hypothesis-threshold are currently not used when fitting for all peaks." );

      const bool highres = PeakFitUtils::is_high_res( spec );
      vector<shared_ptr<const PeakDef>> fit_peaks = ExperimentalAutomatedPeakSearch::search_for_peaks( spec, det, nullptr, false, highres );
      fit_peaks_ptrs.insert( end(fit_peaks_ptrs), begin(fit_peaks), end(fit_peaks) );
    }else
    {
      vector<PeakDef> starting_peaks;

      if( exemplar_peaks )
      {
        for( const auto p : *exemplar_peaks )
          starting_peaks.push_back( *p );
      }else
      {
        // Open the input CSV file - and get those peaks.
        assert( !exemplar_n42 );

#ifdef _WIN32
        const std::wstring wfilename = SpecUtils::convert_from_utf8_to_utf16(exemplar_filename);
        std::ifstream input( wfilename.c_str() );
#else
        std::ifstream input( exemplar_filename.c_str() );
#endif

        if( !input.is_open() )
          throw runtime_error( "Exemplar peak file, '" + exemplar_filename + "', could not be opened." );

        try
        {
          starting_peaks = PeakModel::csv_to_candidate_fit_peaks( spec, input );
        }catch( std::exception &e )
        {
          throw runtime_error( "Invalid input exemplar CSV peak file: " + string(e.what()) );
        }

        if( starting_peaks.empty() )
        {
          results.warnings.push_back( "No candidate peaks for file '" + filename + "', perhaps peaks from CSV were not good candidate peaks -- skipping file." );

          return results;
        }

      }//if( !exemplar_peaks )

      if( starting_peaks.empty() && !options.fit_all_peaks )
      {
        results.warnings.push_back( "No candidate peaks for '" + filename + "' - maybe a programming logic error?" );

        return results;
      }

      // Sort the peaks by mean even though its probably already sorted
      std::sort( begin(starting_peaks), end(starting_peaks), &PeakDef::lessThanByMean );

      results.exemplar_peaks.clear();
      vector<shared_ptr<const PeakDef>> exemplar_peaks;
      for( const auto &p : starting_peaks )
      {
        auto ep = make_shared<PeakDef>(p);
        results.exemplar_peaks.push_back( ep );
        exemplar_peaks.push_back( ep );
        unused_exemplar_peaks.insert( ep );
        results.unfit_exemplar_peaks.push_back( ep ); //We will clear this later on if unsuccessful
      }

      //  We are re-using functions called by the GUI InterSpec, so there are some extra arguments
      //  that arent totally applicable to us right now.
      const double lower_energy = starting_peaks.front().mean() - 0.1;
      const double uppper_energy = starting_peaks.back().mean() + 0.1;

      // We are not refitting peaks, because the areas may be wildly different.
      const bool isRefit = false;

      const double ncausalitysigma = 0.0;
      const double stat_threshold = options.peak_stat_threshold;
      const double hypothesis_threshold = options.peak_hypothesis_threshold;

      vector<PeakDef> candidate_peaks;
      for( const auto &p : starting_peaks )
      {
        PeakDef peak = p;

        // If you had selected for certain peak quantities to not be fit for in exemplar spectrum
        //  in InterSpec, then by default those quantities also wouldnt be fit for here.
        //  You can alter this similar to below.
        //  For the moment we'll just make sure the peak amplitude will be fit for.

        //peak.setFitFor( PeakDef::Mean, true );
        //peak.setFitFor( PeakDef::Sigma, false );
        peak.setFitFor( PeakDef::GaussAmplitude, true );


        // Lets also make sure starting amplitude is something halfway reasonable,
        //  and continuum coefficients are reasonable starting values
        if( peak.gausPeak() )
        {
          std::shared_ptr<PeakContinuum> continuum = peak.continuum();
          assert( continuum );
          if( continuum && continuum->isPolynomial() )
          {
            const PeakContinuum::OffsetType origType = continuum->type();
            // calc_linear_continuum_eqn sets type to Linear and resizes m_values to 2 polynomial
            //  coefficients.  We compute the continuum integral while the type is still Linear
            //  (just polynomial, no step), then restore the original type.
            //  For step types (including CDF step), this gives a reasonable polynomial-only estimate
            //  for the starting amplitude — the step coefficient will be fit later.
            //  Note: setType(origType) will init the step coefficient to zero, so the only-one-peak
            //  vs all-ROI-peers distinction for CDF step types is moot here.
            continuum->calc_linear_continuum_eqn( spec, peak.mean(), peak.lowerX(), peak.upperX(), 2, 2 );

            const double mean = peak.mean(), fwhm = peak.fwhm();
            const double data_area = spec->gamma_integral( mean - fwhm, mean + fwhm );

            if( (data_area > 1) && (peak.amplitude() > data_area) )
            {
              const double cont_area = continuum->offset_integral( mean - fwhm, mean + fwhm,
                                                                    spec, nullptr, 0 );
              if( (cont_area > 0.0) && (cont_area < data_area) )
              {
                peak.setAmplitude( data_area - cont_area );
              }else
              {
                peak.setAmplitude( 0.25*data_area );
              }
            }//if( exemplar peak is clearly much larger than data )

            continuum->setType( origType );
          }//if( continuum )
        }//if( peak.gausPeak() )

        candidate_peaks.push_back( peak );
      }//for( const auto &p : exemplar_peaks )

      results.original_energy_cal = spec ? spec->energy_calibration() : nullptr;
      //if( !results.original_energy_cal )
      //  results.original_energy_cal = make_shared<SpecUtils::EnergyCalibration>( *spec->energy_calibration() ); //Make a copy, just to make sure it doesnt get messed up

      if( options.refit_energy_cal )
      {
        // We will refit the energy calibration - maybe a few times - to really hone in on things
        for( size_t i = 0; i < 4; ++i )
        {
          vector<PeakDef> energy_cal_peaks = candidate_peaks;
          for( auto &peak : energy_cal_peaks )
          {
            peak.setFitFor( PeakDef::Mean, true );
            peak.setFitFor( PeakDef::Sigma, true );
            peak.setFitFor( PeakDef::GaussAmplitude, true );
          }

          vector<PeakDef> peaks = fitPeaksInRange( lower_energy, uppper_energy, ncausalitysigma,
                                                  stat_threshold, hypothesis_threshold,
                                                  energy_cal_peaks, spec, {}, isRefit );

          try
          {
            fit_energy_cal_from_fit_peaks( spec, peaks );
          }catch( std::exception &e )
          {
            const string msg = "Energy calibration not performed: " + string(e.what());
            auto pos = std::find( begin(results.warnings), end(results.warnings), msg );
            if( pos == end(results.warnings) )
              results.warnings.push_back( msg );
          }
        }//for( size_t i = 0; i < 1; ++i )

        // Propagate the updated energy cal to the result file
        assert( spec && spec->energy_calibration() && spec->energy_calibration()->valid() );
        shared_ptr<const SpecUtils::EnergyCalibration> new_cal = spec ? spec->energy_calibration() : nullptr;
        if( new_cal && new_cal->valid() )
        {
          try
          {
            // Must pass the sample numbers being fit - `propagate_energy_cal` only updates the
            //  measurements of the samples given to it, so an empty set updates nothing, and the
            //  written out N42 would keep the original calibration, while the peaks in it were fit
            //  at the new calibration.  Sample numbers chosen same as the `setPeaks(...)` below.
            const set<int> cal_sample_nums = (options.background_subtract_file.empty()
                                              && !options.cached_background_subtract_spec)
                                             ? used_sample_nums : results.measurement->sample_numbers();
            propagate_energy_cal( new_cal, spec, results.measurement, cal_sample_nums );
          }catch( std::exception &e )
          {
            results.warnings.push_back( "Failed to propagate fit energy calibration in '" + filename + "'." );
          }
        }else
        {
          results.warnings.push_back( "Failed to fit an appropriate energy calibration in '" + filename + "'." );
        }

        results.refit_energy_cal = spec ? spec->energy_calibration() : nullptr;
      }//if( options.refit_energy_cal )

      vector<PeakDef> fit_peaks = fitPeaksInRange( lower_energy, uppper_energy, ncausalitysigma,
                                                  stat_threshold, hypothesis_threshold,
                                                  candidate_peaks, spec, {}, isRefit );

#if( !USE_LM_PEAK_FIT )
      // will re-fit the peaks again to make sure have the best solution.
      //  Note: the Ceres/LM-based fitting does not need this - it seems to always fit the best solution on first try
      for( size_t i = 0; i < 3; ++i )
      {
        const vector<PeakDef> prev_peaks = fit_peaks;
        fit_peaks = fitPeaksInRange( lower_energy, uppper_energy, ncausalitysigma,
                                    stat_threshold, hypothesis_threshold,
                                    fit_peaks, spec, {}, true );

        if( fit_peaks.size() != prev_peaks.size() )
          fit_peaks = prev_peaks;
      }
#endif // !USE_LM_PEAK_FIT

      //cout << "Fit for the following " << fit_peaks.size() << " peaks (the exemplar file had "
      //<< starting_peaks.size() <<  ") from the raw spectrum:"
      //<< endl;

      for( PeakDef &p: fit_peaks )
      {
        // Find nearest exemplar peak, and we'll use this to set nuclides, colors, and such
        const double fit_mean = p.mean();

        shared_ptr<const PeakDef> exemplar_parent;
        for( const auto &exemplar : exemplar_peaks )
        {
          const double exemplar_mean = exemplar->mean();

          const double energy_diff = fabs( fit_mean - exemplar_mean );

          // Calculate the overlap fraction of the ROI
          auto overlap_frac_fcn = []( const shared_ptr<const PeakDef> &peak_1, const PeakDef &p_2 ) -> double {
            const double overlapLow = std::max(peak_1->lowerX(), p_2.lowerX());
            const double overlapHigh = std::min(peak_1->upperX(), p_2.upperX() );
            const double intersection = std::max(0.0, overlapHigh - overlapLow);
            if( intersection <= 0.0 )
              return 0.0;
            const double length1 = peak_1->upperX() - peak_1->lowerX();
            const double length2 = p_2.upperX() - p_2.lowerX();
            const double unionLength = length1 + length2 - intersection;
            if( unionLength <= 0.0 )
              return unionLength;
            return intersection / unionLength;
          };

          const double overlap_frac = overlap_frac_fcn( exemplar, p );

          // We will require the fit peak to be within 1.5 FWHM (arbitrarily chosen distance)
          //  of the exemplar peak, and we will use the exemplar peak closest in energy to the
          //  fit peak
          const double fwhm_match_multiple = 1.5;
          const double min_overlap_frac = 0.75; //arbitrary
          if( ((energy_diff < fwhm_match_multiple*p.fwhm())
               || (energy_diff < fwhm_match_multiple*exemplar->fwhm())
               || (overlap_frac > min_overlap_frac) )
             && (!exemplar_parent || (energy_diff < fabs(exemplar_parent->mean() - fit_mean))) )
          {
            exemplar_parent = exemplar;
          }
        }//for( loop over exemplars to find original peak corresponding to the fit peak 'p' )

        if( exemplar_parent )
        {
          unused_exemplar_peaks.erase( exemplar_parent );
          p.inheritUserSelectedOptions( *exemplar_parent, false );
        }else
        {
          results.warnings.push_back( "In '" + filename + "', failed to find exemplar peak"
                                     " corresponding to peak fit with mean=" + std::to_string( p.mean() ) + " keV." );
          //cout << "\tmean=" << p.mean() << ", FWHM=" << p.fwhm() << ", area=" << p.amplitude() << endl;
        }//if( exemplar_parent ) / else
      }//for( PeakDef &p: fit_peaks )

      //cout << endl << endl;

      // Now we will associate the fit peaks to the spectrum and save an N42 file you can open up in
      //  InterSpec and inspect the fits.
      map<shared_ptr<const PeakContinuum>, shared_ptr<PeakContinuum>> peak_continuum_map;
      for( const auto &p: fit_peaks )
      {
        shared_ptr<PeakDef> new_peak = make_shared<PeakDef>(p);
        shared_ptr<const PeakContinuum> cont = new_peak->continuum();
        auto pos = peak_continuum_map.find( cont );
        if( pos == end(peak_continuum_map) )
          pos = peak_continuum_map.insert( { cont, std::make_shared<PeakContinuum>( *cont ) } ).first;
        new_peak->setContinuum( pos->second );
        fit_peaks_ptrs.push_back( new_peak );
      }

      if( options.background_subtract_file.empty() && !options.cached_background_subtract_spec )
      {
        specfile->setPeaks( fit_peaks_ptrs, used_sample_nums );
      }else
      {
        assert( specfile->num_measurements() == 1 );
        specfile->setPeaks( fit_peaks_ptrs, specfile->sample_numbers() );
      }
    }//if( options.fit_all_peaks ) / else

    results.success = true;
    results.fit_peaks = fit_peaks_ptrs;
    results.unfit_exemplar_peaks.clear();
    results.unfit_exemplar_peaks.insert( end(results.unfit_exemplar_peaks),
                                        begin(unused_exemplar_peaks), end(unused_exemplar_peaks) );
    std::sort( begin(results.unfit_exemplar_peaks), end(results.unfit_exemplar_peaks),
              &PeakDef::lessThanByMeanShrdPtr );

    if( !options.fit_all_peaks
       && (options.not_fit_peak_mda != NotFitPeakMdaMethod::None)
       && !results.unfit_exemplar_peaks.empty() )
    {
      const shared_ptr<const DetectorPeakResponse> exemplar_drf = exemplar_n42 ? exemplar_n42->detector()
                                                                               : nullptr;
      results.not_fit_peak_mdas = compute_not_fit_peak_mdas( results.unfit_exemplar_peaks,
                                                     results.fit_peaks, spec, exemplar_drf, options );
    }//if( we should compute detection limits for the peaks that were not fit )
  }
  
  return results;
}//fit_peaks_in_file(...)
  
}//namespace BatchPeak

