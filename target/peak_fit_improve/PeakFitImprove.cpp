//  Created by wcjohns on 20110322.
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

/** This file is for developing peak-fit improving code. */

#include "InterSpec_config.h"

#include <map>
#include <deque>
#include <mutex>
#include <atomic>
#include <chrono>
#include <string>
#include <fstream>
#include <numeric>
#include <iostream>

#include <boost/program_options.hpp>
#include <boost/math/constants/constants.hpp>

//#include "/external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"
#include <Wt/Json/Value>
#include <Wt/Json/Array>
#include <Wt/Json/Parser>
#include <Wt/Json/Object>
#include <Wt/Json/Serializer>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/DateTime.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "SpecUtils/D3SpectrumExport.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakFit.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakDists.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/InterSpecServer.h"
#include "InterSpec/ReferenceLineInfo.h"

#include "InitialFit_GA.h"
#include "PeakFitImprove.h"
#include "CandidatePeak_GA.h"
#include "FinalFit_GA.h"
#include "PeakFitImproveData.h"

using namespace std;






std::string FindCandidateSettings::print( const string &var_name ) const
{
  return var_name + ".num_smooth_side_channels = "   + std::to_string(num_smooth_side_channels)        + ";\n"
  + var_name + ".smooth_polynomial_order = "         + std::to_string(smooth_polynomial_order)         + ";\n"
  + var_name + ".threshold_FOM = "                   + std::to_string(threshold_FOM)                   + ";\n"
  + var_name + ".more_scrutiny_FOM_threshold = "     + std::to_string(more_scrutiny_FOM_threshold)     + ";\n"
  + var_name + ".pos_sum_threshold_sf = "            + std::to_string(pos_sum_threshold_sf)            + ";\n"
  + var_name + ".num_chan_fluctuate = "              + std::to_string(num_chan_fluctuate)              + ";\n"
  //+ var_name + ".min_counts_per_channel = "        + std::to_string(min_counts_per_channel)          + ";\n"
  + var_name + ".more_scrutiny_coarser_FOM = "       + std::to_string(more_scrutiny_coarser_FOM)       + ";\n"
  + var_name + ".more_scrutiny_min_dev_from_line = " + std::to_string(more_scrutiny_min_dev_from_line) + ";\n"
  + var_name + ".amp_to_apply_line_test_below = "    + std::to_string(amp_to_apply_line_test_below)    + ";\n"
  ;
}//std::string print( const string &var_name ) const


std::string FindCandidateSettings::to_json() const
{
    //nlohmann::json data;
  Wt::Json::Object data;
  data["FindCandidateSettingsVersion"] = 1;
  data["NumSmoothChannels"] = num_smooth_side_channels;
  data["SmoothPolyOrder"] = smooth_polynomial_order;
  data["MinFOM"] = threshold_FOM;
  data["MoreScrutinyFOM"] = more_scrutiny_FOM_threshold;
  data["SecondDerivSummingThresh"] = pos_sum_threshold_sf;
  data["NumChannelFluxuate"] = static_cast<int>(num_chan_fluctuate);
  //data["MinCountsPerChannel"] = min_counts_per_channel;
  data["MoreScrutinyCoarseFOM"] = more_scrutiny_coarser_FOM;
  data["MoreScrutinyMinDevFromLine"] = more_scrutiny_min_dev_from_line;
  data["AmpToApplyLineTestBelow"] = amp_to_apply_line_test_below;


  return Wt::Json::serialize( data );
}//std::string to_json() const



// Eval function to:
//  - Takes in a spectrum
//  - Uses optimum `FindCandidateSettings` that is already set, to find candidates
//  - Fits those candidates.
//  - Filters, combines, adds, and modifies initial fits, based on input settings
//    - Filter significance using `chi2_significance_test(...)` function:
//      - "initial_stat_threshold", "initial_hypothesis_threshold"
//      !these two parameters should be optimized separate from everything else!
//    - Combine based on how near means are, and how much initial ROIs are overlapping
//      - "combine_nsigma_near", "combine_ROI_overlap_frac"
//    - Determine coarse FWHM form
//      - "FWHM_interpolation_type" (functional form and order)
//    - Adds peaks based on now kinda known FWHM functional form.
//      Lets just do something stupid, and slide along non-ROI areas, and have an ROI ~5 sigma wide,
//      and draw a line using surrounding channels to fit a flat continuum, and if the deficit is
//      some reasonable value (2.5 sigma? - need to do a time-tradeoff comparison to get something
//      reasonable), then call `fit_amp_and_offset(...)` and `chi2_significance_test(...)`
//      "search_roi_nsigma_deficit", "search_stat_threshold", "search_hypothesis_threshold"
//    - Add peaks to ROIs, based on deficits in the ~1.5 to ~7.5 sigma ranges. Can probably do something
//      stupid like scan a ~2sigma wide area, and look for a deficit of counts.  Or could just go
//      slightly crazy and try calling `fit_amp_and_offset(...)` and look at improvements.
//      "ROI_add_nsigma_required", "ROI_add_chi2dof_improve"
//    - Modify continuum type, based on peak area, and difference between left-and-right side heights,
//      and slope on either side of ROI
//      "cont_type_peak_nsigma_threshold" (how many nsigma the peak is, to consider),
//      "cont_type_left_right_nsigma" (how many stat higher the left is than the right, where stat is based of side edge area),
//      "cont_type_sum_slopes" (to see if linear or bilinear continuum)
//      "cont_type_step_ncrease_chi2dof_required" (refit peak with different steps, and this is the chi2/dof needed to accept solution)
//    - Check if should add skew or not.
//      Right now maybe just just ratio of area between continuum and data for -1 to -4 sigma from mean,
//      compared to peak area, or maybe the ratio of sqrt of the areas.
//      "skew_nsigma"
//    - Determine ROI left and right widths for single peak ROIs
//      Have a base-width, and then add/subtract based off some multiple.
//      The multiple can be peak-area within +-1 FWHM, and continuum-area in +- 1 FWHM
//      Have an int to select subtract type, maybe, linear, sqrt, exp, and log
//      "roi_extent_low_num_fwhm_base", "roi_extent_high_num_fwhm_base"
//      "roi_extent_mult_type" (linear, sqrt),
//      "roi_extent_low_stat_multiple", "roi_extent_high_stat_multiple"
//    - Determine ROI left and right widths for multiple peak ROIs
//      Use single peak peak value for starting, then have an additional add/subtract on each side
//      "multi_roi_extent_low_fwhm_mult", "multi_roi_extent_high_fwhm_mult"
//    - A bool value to see how many times each ROI should be refit
//
//
//
// An example where we should have a stepped continuum
//  503 to 507: data=125.1
//  517 to 521: data=20.5
//  Peak area 19679.5, FWHM=2.1, Mean=511.9

std::string InitialPeakFindSettings::print( const string &var_name ) const
{
  string fwhm_fcn_form_string;
  switch( fwhm_fcn_form )
  {
    case FwhmFcnForm::Gadras:                   fwhm_fcn_form_string = "Gadras";                   break;
    case FwhmFcnForm::SqrtPolynomialTwoCoefs:   fwhm_fcn_form_string = "SqrtPolynomialTwoCoefs";   break;
    case FwhmFcnForm::SqrtPolynomialThreeCoefs: fwhm_fcn_form_string = "SqrtPolynomialThreeCoefs"; break;
    case FwhmFcnForm::SqrtEnergyPlusInverse:    fwhm_fcn_form_string = "SqrtEnergyPlusInverse";    break;
    case FwhmFcnForm::NumFwhmFcnForm:           fwhm_fcn_form_string = "NumFwhmFcnForm";           break;
  }//switch( fwhm_fcn_form )

  return var_name + ".initial_stat_threshold = "  + std::to_string(initial_stat_threshold)       + ";\n"
  + var_name + ".initial_hypothesis_threshold = " + std::to_string(initial_hypothesis_threshold) + ";\n"
  + var_name + ".initial_min_nsigma_roi = "       + std::to_string(initial_min_nsigma_roi)       + ";\n"
  + var_name + ".initial_max_nsigma_roi = "       + std::to_string(initial_max_nsigma_roi)       + ";\n"
  + var_name + ".fwhm_fcn_form_string = "         + fwhm_fcn_form_string                         + ";\n"
  + var_name + ".search_roi_nsigma_deficit = "    + std::to_string(search_roi_nsigma_deficit)    + ";\n"
  + var_name + ".search_stat_threshold = "        + std::to_string(search_stat_threshold)        + ";\n"
  + var_name + ".search_hypothesis_threshold = "  + std::to_string(search_hypothesis_threshold)  + ";\n"
  + var_name + ".search_stat_significance = "     + std::to_string(search_stat_significance)     + ";\n"
  + var_name + ".ROI_add_nsigma_required = "      + std::to_string(ROI_add_nsigma_required)      + ";\n"
  + var_name + ".ROI_add_chi2dof_improve = "      + std::to_string(ROI_add_chi2dof_improve)      + ";\n"
  ;
}//std::string print( const string &var_name ) const
  

std::string FinalPeakFitSettings::print( const string &var_name ) const
{
  return var_name + ".combine_nsigma_near = "                    + std::to_string(combine_nsigma_near)                         + ";\n"
  + var_name + ".combine_ROI_overlap_frac = "                    + std::to_string(combine_ROI_overlap_frac)                    + ";\n"
  + var_name + ".cont_type_peak_nsigma_threshold = "             + std::to_string(cont_type_peak_nsigma_threshold)             + ";\n"
  + var_name + ".cont_type_left_right_nsigma = "                 + std::to_string(cont_type_left_right_nsigma)                 + ";\n"
  //+ var_name + ".stepped_roi_extent_lower_side_stat_multiple = " + std::to_string(stepped_roi_extent_lower_side_stat_multiple) + ";\n"
  //+ var_name + ".stepped_roi_extent_upper_side_stat_multiple = " + std::to_string(stepped_roi_extent_upper_side_stat_multiple) + ";\n"
  + var_name + ".cont_poly_order_increase_chi2dof_required = "   + std::to_string(cont_poly_order_increase_chi2dof_required)   + ";\n"
  + var_name + ".cont_step_type_increase_chi2dof_required = "    + std::to_string(cont_step_type_increase_chi2dof_required)    + ";\n"
  + var_name + ".skew_nsigma = "                                 + std::to_string(skew_nsigma)                                 + ";\n"
  + var_name + ".left_residual_sum_min_to_try_skew = "           + std::to_string(left_residual_sum_min_to_try_skew)           + ";\n"
  + var_name + ".right_residual_sum_min_to_try_skew = "          + std::to_string(right_residual_sum_min_to_try_skew)          + ";\n"
  + var_name + ".skew_improve_chi2_dof_threshold = "             + std::to_string(skew_improve_chi2_dof_threshold)             + ";\n"
  + var_name + ".roi_extent_low_num_fwhm_base = "                + std::to_string(roi_extent_low_num_fwhm_base)                + ";\n"
  + var_name + ".roi_extent_high_num_fwhm_base = "               + std::to_string(roi_extent_high_num_fwhm_base)               + ";\n"
  + var_name + ".roi_extent_mult_type = "                        + std::string((roi_extent_mult_type == RoiExtentMultType::Linear)
                                                                               ? "RoiExtentMultType::Linear" : "RoiExtentMultType::Sqrt") + ";\n"
  + var_name + ".roi_extent_lower_side_stat_multiple = "         + std::to_string(roi_extent_lower_side_stat_multiple)         + ";\n"
  + var_name + ".roi_extent_upper_side_stat_multiple = "         + std::to_string(roi_extent_upper_side_stat_multiple)         + ";\n"
  + var_name + ".multi_roi_extent_lower_side_fwhm_mult = "       + std::to_string(multi_roi_extent_lower_side_fwhm_mult)       + ";\n"
  + var_name + ".multi_roi_extent_upper_side_fwhm_mult = "       + std::to_string(multi_roi_extent_upper_side_fwhm_mult)       + ";\n";
}//std::string print( const string &var_name ) const




PeakTruthInfo::PeakTruthInfo()
: energy( 0.0 ), cps( 0.0 ), area( 0.0 ), fwhm( 0.0 ), full_width( 0.0 ), label()
{
}

PeakTruthInfo::PeakTruthInfo( const std::string &line )
{
  vector<string> fields;
  SpecUtils::split( fields, line, "," );
  assert( fields.size() == 6 );
  if( fields.size() != 6 )
    throw runtime_error( "Unexpected number of fileds in line: '" + line + "'" );

  //"Energy (keV),CPS,Area,FWHM,FullWidth,Label"
  if( !SpecUtils::parse_float( fields[0].c_str(), fields[0].size(), energy ) )
    throw runtime_error( "Failed to parse energy of line: '" + line + "'" );

  if( !SpecUtils::parse_float( fields[1].c_str(), fields[1].size(), cps ) )
    throw runtime_error( "Failed to parse cps of line: '" + line + "'" );

  if( !SpecUtils::parse_float( fields[2].c_str(), fields[2].size(), area ) )
    throw runtime_error( "Failed to parse area of line: '" + line + "'" );

  if( !SpecUtils::parse_float( fields[3].c_str(), fields[3].size(), fwhm ) )
    throw runtime_error( "Failed to parse fwhm of line: '" + line + "'" );

  if( !SpecUtils::parse_float( fields[4].c_str(), fields[4].size(), full_width ) )
    throw runtime_error( "Failed to parse full_width of line: '" + line + "'" );

  label = fields[5];
  assert( (label == "S.E.") || (label == "D.E.")
         || (label == "EscapeXRay") || (label == "Peak")
         || (label == "EscapeXRay+Peak") || (label == "D.E.+EscapeXRay")
         || (label == "D.E.+Peak") || (label == "D.E.+S.E.")
         || (label == "Peak+S.E.")
         || (label == "EscapeXRay+S.E.") );
}//PeakTruthInfo(...)








int main( int argc, char **argv )
{
  const double start_wall = SpecUtils::get_wall_time();
  const double start_cpu = SpecUtils::get_cpu_time();
  
  string datadir;
  if( datadir.empty() )
  {
    string targetfile = "data/CharacteristicGammas.txt";
    if( !AppUtils::locate_file( targetfile, false, 5, false ) )
    {
      cerr << "Unable to find '" << targetfile << "' directory." << endl;
      return -24;
    }
    
    datadir = SpecUtils::parent_path( targetfile );
  }//if( docroot.empty() )
  
  assert( SpecUtils::is_directory(datadir) );
  
  InterSpec::setStaticDataDirectory( datadir );
  
  /*
   Known issues with "truth" CSV files
   - The truth lines are before random-summing, so the observed peaks will be smaller than the CSV says
   - The annihilation 511 keV will not be included (sources like Na22 will have 511 keV, but "truth" will claim be slightly smaller than actual)
   - Reaction lines not be accounted for, e.g. 488 keV Alpha-Li reaction
   - Single and double escape peaks not present in the "truth" lines.
   - I dont think sum lines are accounted for
   */

  const string base_dir = "/Users/wcjohns/rad_ana/peak_area_optimization/peak_fit_accuracy_inject/";

  vector<string> hpges {
    "Detective-X",
    "Detective-EX",
    "Detective-X_noskew",
    "Falcon 5000",
    "Fulcrum40h",
    "LANL_X",
    "HPGe_Planar_50%"
  };

  if( PeakFitImprove::debug_printout )
    hpges = vector<string>{ "Detective-X" };

#if( WRITE_ALL_SPEC_TO_HTML )
  hpges = vector<string>{ "Detective-X", "IdentiFINDER-R500-NaI" };
#endif


  vector<string> live_times{
    "30_seconds",
    "300_seconds",
    "1800_seconds"
  };

#if( WRITE_ALL_SPEC_TO_HTML )
  live_times = {"300_seconds"};
#endif


  const std::tuple<std::vector<DetectorInjectSet>,std::vector<DataSrcInfo>> &loaded_data
      = PeakFitImproveData::load_inject_data_with_truth_info( base_dir, hpges, live_times );

  const vector<DetectorInjectSet> &inject_sets = std::get<0>(loaded_data);
  const vector<DataSrcInfo> &input_srcs = std::get<1>(loaded_data);



#if( WRITE_ALL_SPEC_TO_HTML )
  PeakFitImproveData::write_html_summary( input_srcs );
  return 1;
#endif


  
  
  /**
   Default settings:
   */
  FindCandidateSettings best_settings;
  best_settings.num_smooth_side_channels = 4; // low res more
  best_settings.smooth_polynomial_order = 3;  // highres 3, lowres 2
  best_settings.threshold_FOM = 1.3;
  best_settings.more_scrutiny_FOM_threshold = 3.5;
  best_settings.more_scrutiny_coarser_FOM = 5.0;
  best_settings.pos_sum_threshold_sf = -0.01f;
  best_settings.num_chan_fluctuate = 2;
  //best_settings.min_counts_per_channel = 1.5;
  best_settings.more_scrutiny_min_dev_from_line = 4.0;
  best_settings.amp_to_apply_line_test_below = 40;
  
  
  double best_sum_score = -1.0E-6;
  size_t best_score_num_peaks_found = 0, best_score_num_possibly_accepted_peaks_not_found = 0, best_score_num_extra_peaks = 0, best_score_def_wanted_peaks_not_found = 0, best_score_def_wanted_peaks_found = 0;
  
  const double eval_start_wall = SpecUtils::get_wall_time();
  
  bool done_posting = false;
  size_t num_evaluations = 0, last_percent_printed = 0, num_posted = 0;
  std::mutex score_mutex;
  auto eval_candidate_settings_fcn = [&]( const FindCandidateSettings settings, const vector<DataSrcInfo> &input_srcs, const bool write_n42 ){

    const tuple<double,size_t,size_t,size_t,size_t,size_t> result = CandidatePeak_GA::eval_candidate_settings( settings, input_srcs, write_n42 );

    const double score = std::get<0>(result);
    const size_t num_peaks_found = std::get<1>(result);
    const size_t num_def_wanted_not_found = std::get<2>(result);
    const size_t def_wanted_peaks_found = std::get<3>(result);
    const size_t num_possibly_accepted_peaks_not_found = std::get<4>(result);
    const size_t num_extra_peaks = std::get<5>(result);
    
    std::lock_guard<std::mutex> lock( score_mutex );
    
    num_evaluations += 1;
    const size_t percent_done = (100*num_evaluations) / num_posted;
    
    if( done_posting && (percent_done > last_percent_printed) )
    {
      const double now_wall = SpecUtils::get_wall_time();
      
      last_percent_printed = percent_done;
      
      const double elapsed = now_wall - eval_start_wall;
      const double frac_eval = num_evaluations / static_cast<double>(num_posted);
      const double est_total_time = elapsed / frac_eval;
      cout << "\t" << percent_done << "% done. " << elapsed << " of "
      << est_total_time << " (estimated) seconds through." << endl;
    }//if( done_posting && (percent_done > last_percent_printed) )
    
    if( score > best_sum_score )
    {
      best_sum_score = score;
      best_settings = settings;
      best_score_num_peaks_found = num_peaks_found;
      best_score_def_wanted_peaks_not_found = num_def_wanted_not_found;
      best_score_def_wanted_peaks_found = def_wanted_peaks_found;
      best_score_num_possibly_accepted_peaks_not_found = num_possibly_accepted_peaks_not_found;
      best_score_num_extra_peaks = num_extra_peaks;
    }
  };
  
  
  //best_settings.num_smooth_side_channels = 13; // low res more
  //best_settings.smooth_polynomial_order = 2;  // highres 3, lowres 2
  //best_settings.threshold_FOM = 1.3;
  //best_settings.pos_sum_threshold_sf = 0.16f;
  //eval_candidate_settings_fcn( best_settings, input_srcs, false );

  
  enum class OptimizationAction : int
  {
    Candidate,
    InitialFit,
    FinalFit,
    CodeDev,
    AccuracyFromCsvsStudy
  };//enum class OptimizationAction : int
  
  const OptimizationAction action = OptimizationAction::CodeDev; //OptimizationAction::InitialFit;

  switch( action )
  {
    case OptimizationAction::Candidate:
    {
      // code to run candidate peak optimization
      const auto ga_eval = [&input_srcs](const FindCandidateSettings &settings) -> double {
        const tuple<double,size_t,size_t,size_t,size_t,size_t> score = CandidatePeak_GA::eval_candidate_settings( settings, input_srcs, false );
        return get<0>( score );
      };
      
      best_settings = CandidatePeak_GA::do_ga_eval( ga_eval );
      
      /*
      SpecUtilsAsync::ThreadPool pool;
      
      // With this many nested loops, its really easy for the number of iterations to explode to an
      // intractable amount
      //for( int num_side = 8; num_side < 11; ++num_side )
      for( int num_side = 7; num_side < 12; ++num_side )
      {
        for( int poly_order = 2; poly_order < 3; ++poly_order )
        {
          //for( double threshold_FOM = 1.3; threshold_FOM < 2.0; threshold_FOM += 0.2 )
          for( double threshold_FOM = 1.1; threshold_FOM < 1.8; threshold_FOM += 0.1 )
          {
            for( double more_scrutiny_FOM_threshold = threshold_FOM-0.001; more_scrutiny_FOM_threshold < (threshold_FOM + 1.25); more_scrutiny_FOM_threshold += 0.05 )
            //for( double more_scrutiny_FOM_threshold = 0; more_scrutiny_FOM_threshold == 0; more_scrutiny_FOM_threshold += 0.25 )
            {
              //for( float pos_sum_threshold_sf = -0.06f; pos_sum_threshold_sf < 0.01; pos_sum_threshold_sf += 0.02 )
              for( float pos_sum_threshold_sf = 0.0f; pos_sum_threshold_sf < 0.01; pos_sum_threshold_sf += 0.02 )
              {
                for( float more_scrutiny_coarser_FOM = threshold_FOM; more_scrutiny_coarser_FOM < (threshold_FOM + 1.0); more_scrutiny_coarser_FOM += 0.2 )
                //for( float more_scrutiny_coarser_FOM = threshold_FOM; more_scrutiny_coarser_FOM <= threshold_FOM; more_scrutiny_coarser_FOM += 1.0 )
                {
                  for( float more_scrutiny_min_dev_from_line = 2; more_scrutiny_min_dev_from_line < 4.0; more_scrutiny_min_dev_from_line += 0.5 )
                  {
                    for( float amp_to_apply_line_test_below = 50; amp_to_apply_line_test_below < 110; amp_to_apply_line_test_below += 10 )
                    {
                      {
                        std::lock_guard<std::mutex> lock( score_mutex );
                        ++num_posted;
                      }
                      
                      FindCandidateSettings settings;
                      settings.num_smooth_side_channels = num_side; // low res more
                      settings.smooth_polynomial_order = poly_order;  // highres 3, lowres 2
                      settings.threshold_FOM = threshold_FOM;
                      settings.more_scrutiny_FOM_threshold = more_scrutiny_FOM_threshold;
                      settings.pos_sum_threshold_sf = pos_sum_threshold_sf;
                      //settings.min_counts_per_channel = min_counts_per_ch;
                      settings.more_scrutiny_coarser_FOM = more_scrutiny_coarser_FOM;
                      settings.more_scrutiny_min_dev_from_line = more_scrutiny_min_dev_from_line;
                      settings.amp_to_apply_line_test_below = amp_to_apply_line_test_below;
                      
                      pool.post( [&input_srcs,settings,&score_mutex](){
                        try
                        {
                          eval_candidate_settings_fcn( settings, input_srcs, false );
                        }catch( std::exception &e )
                        {
                          //std::lock_guard<std::mutex> lock( score_mutex );
                          //cerr << "Caught exception: " << e.what() << ", for:" << endl
                          //<< settings.print("\tsettings") << endl;
                        }
                      } );
                    }//for( loop over amp_to_apply_line_test_below )
                  }//for( loop over more_scrutiny_min_dev_from_line )
                }//for( loop over more_scrutiny_coarser_FOM )
              }//for( loop over pos_sum_threshold_sf )
            }//for( loop over more_scrutiny_FOM_threshold )
          }//for( loop over threshold_FOM )
        }//for( loop over poly_order )
      }//for( loop over num_side )
      
      {
        std::lock_guard<std::mutex> lock( score_mutex );
        done_posting = true;
      }
      cout << "Posted " << num_posted << " evaluations." << endl;
      
      pool.join();
       */
      
      eval_candidate_settings_fcn( best_settings, input_srcs, true );
      cout << "Wrote N42s with best settings." << endl;
      
      const double end_wall = SpecUtils::get_wall_time();
      const double end_cpu = SpecUtils::get_cpu_time();
      auto now = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());
      
      stringstream outmsg;
      outmsg << "Finish Time: " << SpecUtils::to_iso_string(now) << endl
      << endl
      << best_settings.print("\tbest_settings") << endl
      << endl
      << "found_extra_punishment = " << JudgmentFactors::found_extra_punishment << endl
      << "Best settings had score " << best_sum_score << ", with values:" << endl
      << best_settings.print("\tsettings")
      << "And:\n"
      << "\tnum_peaks_found:" << best_score_num_peaks_found << endl
      << "\tdef_wanted_peaks_found:" << best_score_def_wanted_peaks_found << endl
      << "\tdef_wanted_peaks_not_found:" << best_score_def_wanted_peaks_not_found << endl
      << "\tnum_possibly_accepted_peaks_not_found:" << best_score_num_possibly_accepted_peaks_not_found << endl
      << "\tnum_extra_peaks:" << best_score_num_extra_peaks << endl
      << endl << endl
      << "Ran in {wall=" << (end_wall - start_wall)
      << ", cpu=" << (end_cpu - start_cpu) << "} seconds" << endl;

      {
        ofstream best_settings_file( "best_candidate_settings.txt" );
        best_settings_file << outmsg.str();
      }
      
      cout << outmsg.str()<< endl;
      break;
    }//case OptimizationAction::Candidate:
      
    case OptimizationAction::InitialFit:
    {
      // FindCandidateSettings:
      //   Apply best settings found; Generation 491, population best score: -18.1076, population average score: -14.9656
      best_settings.num_smooth_side_channels = 9;
      best_settings.smooth_polynomial_order = 2;
      best_settings.threshold_FOM = 1.127040;
      best_settings.more_scrutiny_FOM_threshold = best_settings.threshold_FOM + 0.498290;
      best_settings.pos_sum_threshold_sf = 0.081751;
      best_settings.num_chan_fluctuate = 1;
      best_settings.more_scrutiny_coarser_FOM = best_settings.threshold_FOM + 1.451548;
      best_settings.more_scrutiny_min_dev_from_line = 5.866464;
      best_settings.amp_to_apply_line_test_below = 6;
      
      std::function<double( const InitialPeakFindSettings &)> ga_eval_fcn 
            = [best_settings, &input_srcs]( const InitialPeakFindSettings &settings ) -> double {
        double score_sum = 0.0;
        for( const DataSrcInfo &info : input_srcs )
        {
          const double score = InitialFit_GA::eval_initial_peak_find_and_fit( settings, best_settings, info ).find_weight;
          score_sum += score;
        }
        
        return -score_sum;
      };// set InitialFit_GA::ns_ga_eval_fcn
      
      const InitialPeakFindSettings best_initial_fit_settings = InitialFit_GA::do_ga_eval( ga_eval_fcn );
      
      
      //eval_candidate_settings_fcn( best_settings, best_initial_fit_settings, input_srcs, true );
      //cout << "Wrote N42s with best settings." << endl;
      
      const double end_wall = SpecUtils::get_wall_time();
      const double end_cpu = SpecUtils::get_cpu_time();
      auto now = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());
      
      stringstream outmsg;
      outmsg << "Finish Time: " << SpecUtils::to_iso_string(now) << endl
      << endl
      << "Candidate settings used:\n"
      << best_settings.print("\tbest_settings") << endl
      << endl
      << best_initial_fit_settings.print("\tbest_initial_fit_settings")
      //<< "Best settings had score " << best_sum_score << ", with values:" << endl
      << endl << endl
      << "Ran in {wall=" << (end_wall - start_wall)
      << ", cpu=" << (end_cpu - start_cpu) << "} seconds" << endl;

      {
        ofstream best_settings_file( "best_candidate_settings.txt" );
        best_settings_file << outmsg.str();
      }
      
      cout << outmsg.str()<< endl;
      
      break;
    }//case OptimizationAction::InitialFit:
      
      
    case OptimizationAction::FinalFit:
    {
      FindCandidateSettings candidate_settings;
      candidate_settings.num_smooth_side_channels = 9;
      candidate_settings.smooth_polynomial_order = 2;
      candidate_settings.threshold_FOM = 1.127040;
      candidate_settings.more_scrutiny_FOM_threshold = best_settings.threshold_FOM + 0.498290;
      candidate_settings.pos_sum_threshold_sf = 0.081751;
      candidate_settings.num_chan_fluctuate = 1;
      candidate_settings.more_scrutiny_coarser_FOM = best_settings.threshold_FOM + 1.451548;
      candidate_settings.more_scrutiny_min_dev_from_line = 5.866464;
      candidate_settings.amp_to_apply_line_test_below = 6;
      best_settings = candidate_settings;
      
      
      InitialPeakFindSettings initial_fit_settings;
      //20240806: Generation [84], Best=-62960.8, Average=-58257.3, , Best generation yet
      initial_fit_settings.initial_stat_threshold = 2.890059;
      initial_fit_settings.initial_hypothesis_threshold = 0.746688;
      initial_fit_settings.initial_min_nsigma_roi = 3.801792;
      initial_fit_settings.initial_max_nsigma_roi = 7.315852;
      initial_fit_settings.fwhm_fcn_form = InitialPeakFindSettings::FwhmFcnForm::SqrtPolynomialTwoCoefs;
      initial_fit_settings.search_roi_nsigma_deficit = 5.642556;
      initial_fit_settings.search_stat_threshold = 2.951359;
      initial_fit_settings.search_hypothesis_threshold= 4.328742;
      initial_fit_settings.search_stat_significance = 1.814627;
      initial_fit_settings.ROI_add_nsigma_required = 3.654512;
      initial_fit_settings.ROI_add_chi2dof_improve = 0.374546;
      
      
      FinalFit_GA::do_final_peak_fit_ga_optimization( candidate_settings, initial_fit_settings, input_srcs );
      
      break;
    }//case OptimizationAction::FinalFit
      
      
    case OptimizationAction::AccuracyFromCsvsStudy:
    {
      const string result_base_path = "/Users/wcjohns/rad_ana/InterSpec/target/peak_fit_improve/build_xcode/PeakSearchCsvs";
      
      shared_ptr<DetectorPeakResponse> drf = DrfSelect::initARelEffDetector( SpecUtils::DetectorType::DetectiveX, "Ortec", "Detective-X" );
      
      size_t total_num_expected = 0;
      size_t peak_easy_num_found = 0, peak_easy_num_extra = 0;
      size_t gad4s_num_found = 0, gad4s_num_extra = 0;
      size_t gad5s_num_found = 0, gad5s_num_extra = 0;
      size_t g2k_num_found = 0, g2k_num_extra = 0;
      size_t interspec_num_found = 0, interspec_num_extra = 0;
      
      vector<float> peak_easy_nsigma_away, gadras_nsigma_away, g2k_nsigma_away, interspec_nsigma_away;
      
      vector<float> significances;
      for( const DataSrcInfo &info : input_srcs )
      {
        string result_path = SpecUtils::append_path( result_base_path, info.detector_name );
        result_path = SpecUtils::append_path( result_path, info.location_name );
        result_path = SpecUtils::append_path( result_path, info.live_time_name );

        const string peak_easy_result_dir = SpecUtils::append_path( result_path, "PeakEasy" );
        const string gadras_4sigma_result_dir = SpecUtils::append_path( result_path, "GADRAS_4sigma" );
        const string gadras_5sigma_result_dir = SpecUtils::append_path( result_path, "GADRAS_5sigma" );
        const string g2k_result_dir = SpecUtils::append_path( result_path, "G2k_reports/NBSSTD" );
        
        const string peak_easy_file = SpecUtils::append_path( peak_easy_result_dir, "PeakEasy_" + info.src_info.src_name + ".CSV" );
        const string gad4s_file = SpecUtils::append_path( gadras_4sigma_result_dir, info.src_info.src_name + "-Pks.CSV" );
        const string gad5s_file = SpecUtils::append_path( gadras_5sigma_result_dir, info.src_info.src_name + "-Pks.CSV" );
        
        
        string g2k_file;
        const vector<string> g2k_results = SpecUtils::recursive_ls( g2k_result_dir );
        for( const string &filename : g2k_results )
        {
          if( SpecUtils::iequals_ascii( SpecUtils::filename(filename) , info.src_info.src_name + ".RPT") )
          {
            assert( g2k_file.empty() );
            g2k_file = filename;
            //break;
          }
        }//for( const string &filename : g2k_results )
        assert( !g2k_file.empty() );
        
        assert( SpecUtils::is_file(peak_easy_file) );
        assert( SpecUtils::is_file(gad4s_file) );
        assert( SpecUtils::is_file(gad5s_file) );
        assert( SpecUtils::is_file(g2k_file) );
        
        ifstream peakeasy_strm( peak_easy_file.c_str(), (ios::binary | ios::in) );
        assert( peakeasy_strm.good() );
        if( !peakeasy_strm.good() )
          throw runtime_error( "Failed to open CSV: " + peak_easy_file );
        
        ifstream gad4s_strm( gad4s_file.c_str(), (ios::binary | ios::in) );
        assert( gad4s_strm.good() );
        if( !gad4s_strm.good() )
          throw runtime_error( "Failed to open CSV: " + gad4s_file );
        
        ifstream gad5s_strm( gad5s_file.c_str(), (ios::binary | ios::in) );
        assert( gad5s_strm.good() );
        if( !gad5s_strm.good() )
          throw runtime_error( "Failed to open CSV: " + gad5s_file );
        
        ifstream g2k_strm( g2k_file.c_str(), (ios::binary | ios::in) );
        assert( g2k_strm.good() );
        if( !g2k_strm.good() )
          throw runtime_error( "Failed to open RPT: " + g2k_file );
        
        shared_ptr<const SpecUtils::Measurement> data = info.src_info.src_spectra[0];
        
        //
        //vector<PeakDef> PeakModel::csv_to_candidate_fit_peaks( std::shared_ptr<const SpecUtils::Measurement> meas, std::istream &csv );
        
        const vector<PeakDef> peak_easy_peaks = PeakModel::gadras_peak_csv_to_peaks( data, peakeasy_strm );
        const vector<PeakDef> gad4s_peaks = PeakModel::gadras_peak_csv_to_peaks( data, gad4s_strm );
        const vector<PeakDef> gad5s_peaks = PeakModel::gadras_peak_csv_to_peaks( data, gad5s_strm );
        const vector<PeakFitImproveData::G2k::G2kPeak> g2k_peaks_initial = PeakFitImproveData::G2k::g2k_peaks_from_file( g2k_strm );
        vector<PeakDef> g2k_peaks;
        for( const PeakFitImproveData::G2k::G2kPeak &p : g2k_peaks_initial )
        {
          g2k_peaks.emplace_back( p.Energy, p.FWHM/PhysicalUnits::fwhm_nsigma, p.NetPeakArea );
          g2k_peaks.back().setAmplitudeUncert( p.NetAreaError ); //We get 2-sigma errors from G2k
        }
        // As of 20250418, `search_for_peaks(...)` still partially uses Minuit2 based peak-fitting,
        //  even when USE_CERES_PEAK_FITTING is true.
        vector<shared_ptr<const PeakDef> > interspec_peaks_initial
            = ExperimentalAutomatedPeakSearch::search_for_peaks( data, drf, nullptr, false, true );
        vector<PeakDef> interspec_peaks;
        for( const auto &p : interspec_peaks_initial )
          interspec_peaks.push_back( *p );
        
        
        // Get dist of stat sig of peaks
        for( const ExpectedPhotopeakInfo &expected : info.expected_photopeaks )
        {
          //expected.effective_energy
          significances.push_back( expected.nsigma_over_background );
        }
        
        const float min_nsigma = 5, min_counts = 15;
        size_t num_expected = 0;
        for( const ExpectedPhotopeakInfo &expected : info.expected_photopeaks )
        {
          if( expected.nsigma_over_background >= min_nsigma && expected.peak_area > min_counts )
            num_expected += 1;
        }
        
        total_num_expected += num_expected;
        
        auto num_found_or_not = [&info,min_nsigma,min_counts]( const vector<PeakDef> &peaks, vector<float> &nsigma_away ) -> pair<size_t,size_t> {
          vector<bool> fit_peaks_matched_to_expected( peaks.size(), false );
          vector<bool> expected_peaks_matched_to_fit( info.expected_photopeaks.size(), false );
          
          for( size_t i = 0; i < peaks.size(); ++i )
          {
            float expected_area = 0.0f;
            const PeakDef &peak = peaks[i];
            // Avoid Alpha-Li reaction and 511
            if( (fabs(peak.mean() - 478.0) < 3) || (fabs(peak.mean() - 511.0) < 2.0) )
            {
              continue;
            }
            
            if( info.expected_photopeaks.size() )
            {
              double nearest_de = DBL_MAX, largest_sig = 0;
              size_t nearest_expected = info.expected_photopeaks.size();
              for( size_t j = 0; j < info.expected_photopeaks.size(); ++j )
              {
                const ExpectedPhotopeakInfo &expected = info.expected_photopeaks[j];
                const double exp_energy = expected.effective_energy;
                const double delta_e = exp_energy - peak.mean();
                
                if( peak.mean() >= expected.roi_lower
                   && peak.mean() <= expected.roi_upper
                  && (delta_e < 0.5*expected.gamma_lines.front().fwhm)
                   //&& expected.nsigma_over_background > 4
                   //&& expected.nsigma_over_background < 10
                   )
                {
                  if( largest_sig < expected.nsigma_over_background )
                  {
                    largest_sig = expected.nsigma_over_background;
                    expected_area = expected.peak_area;
                  }
                }
                
                if( peak.mean() >= expected.roi_lower
                   && peak.mean() <= expected.roi_upper
                   && fabs(delta_e) < fabs(nearest_de) )
                {
                  fit_peaks_matched_to_expected[i] = true;
                  
                  if( expected.nsigma_over_background > min_nsigma
                     && expected.peak_area > min_counts )
                  {
                    nearest_expected = j;
                    nearest_de = delta_e;
                  }
                }
              }
              
              if( nearest_expected < info.expected_photopeaks.size() )
              {
                expected_peaks_matched_to_fit[nearest_expected] = true;
              }
            }//if( info.expected_photopeaks.size() )
            
            
            if( expected_area > 0 )
            {
              nsigma_away.push_back( (peak.amplitude() - expected_area) / peak.amplitudeUncert() );
            }
            
          }//for( const PeakDef &peak : peaks )
          
          
          size_t num_found = 0, num_not_expected = 0;
          
          for( size_t j = 0; j < info.expected_photopeaks.size(); ++j )
          {
            const ExpectedPhotopeakInfo &expected = info.expected_photopeaks[j];
            if( expected.nsigma_over_background > min_nsigma
               && expected.peak_area > min_counts )
            {
              if( expected_peaks_matched_to_fit[j] )
                num_found += 1;
            }
          }
          
          for( size_t i = 0; i < peaks.size(); ++i )
          {
            const PeakDef &peak = peaks[i];
            if( (fabs(peak.mean() - 478.0) > 3) && (fabs(peak.mean() - 511.0) > 2.0) )
            {
              if( !fit_peaks_matched_to_expected[i] )
                num_not_expected += 1;
            }
          }
          
          return make_pair( num_found, num_not_expected );
        };//num_found_or_not(...)
        
        
        const pair<size_t,size_t> peak_easy_found_not = num_found_or_not(peak_easy_peaks, peak_easy_nsigma_away);
        peak_easy_num_found += peak_easy_found_not.first;
        peak_easy_num_extra += peak_easy_found_not.second;
        
        const pair<size_t,size_t> gad_4s_found_not = num_found_or_not(gad4s_peaks, gadras_nsigma_away);
        gad4s_num_found += gad_4s_found_not.first;
        gad4s_num_extra += gad_4s_found_not.second;
        
        vector<float> dummy;
        const pair<size_t,size_t> gad_5s_found_not = num_found_or_not(gad5s_peaks, dummy);
        gad5s_num_found += gad_5s_found_not.first;
        gad5s_num_extra += gad_5s_found_not.second;
        
        const pair<size_t,size_t> g2k_found_not = num_found_or_not(g2k_peaks, g2k_nsigma_away);
        g2k_num_found += g2k_found_not.first;
        g2k_num_extra += g2k_found_not.second;
        
        const pair<size_t,size_t> interspec_found_not = num_found_or_not(interspec_peaks, interspec_nsigma_away);
        interspec_num_found += interspec_found_not.first;
        interspec_num_extra += interspec_found_not.second;
        
      }//for( const DataSrcInfo &info : input_srcs )
     
      ofstream output_significances( "truth_level_significances.csv" );
      for( const float val : significances )
        output_significances << val << endl;
      
      
      cout << "Expected total = " << total_num_expected << endl;
      cout << "Program, Num Expected Found, Not Found, Num Not Expected" << endl;
      cout << "GADRAS 4-sigma, " << gad4s_num_found << ", " << (total_num_expected - gad4s_num_found) << ", " << gad4s_num_extra << endl;
      cout << "GADRAS 5-sigma, " << gad5s_num_found << ", " << (total_num_expected - gad5s_num_found)<< ", " << gad5s_num_extra << endl;
      cout << "PeakEasy 5.21, " << peak_easy_num_found << ", " << (total_num_expected - peak_easy_num_found)<< ", " << peak_easy_num_extra << endl;
      cout << "InterSpec 1.0.12, " << interspec_num_found << ", " << (total_num_expected - interspec_num_found)<< ", " << interspec_num_extra << endl;
      cout << "Genie 2k, " << g2k_num_found << ", " << (total_num_expected - g2k_num_found) << ", " << g2k_num_extra << endl;
      
      
      {
        ofstream output_nsigma_off( "peakeasy_nsigma_away.csv" );
        for( const float val : peak_easy_nsigma_away )
          output_nsigma_off << val << endl;
      }
      
      {
        ofstream output_nsigma_off( "gadras_nsigma_away.csv" );
        for( const float val : gadras_nsigma_away )
          output_nsigma_off << val << endl;
      }
      
      {
        ofstream output_nsigma_off( "interspec_nsigma_away.csv" );
        for( const float val : interspec_nsigma_away )
          output_nsigma_off << val << endl;
      }
      
      {
        ofstream output_nsigma_off( "g2k_nsigma_away.csv" );
        for( const float val : g2k_nsigma_away )
          output_nsigma_off << val << endl;
      }
      
      if( false ) //output some N42 files
      {
        FindCandidateSettings candidate_settings;
        candidate_settings.num_smooth_side_channels = 9;
        candidate_settings.smooth_polynomial_order = 2;
        candidate_settings.threshold_FOM = 1.127040;
        candidate_settings.more_scrutiny_FOM_threshold = best_settings.threshold_FOM + 0.498290;
        candidate_settings.pos_sum_threshold_sf = 0.081751;
        candidate_settings.num_chan_fluctuate = 1;
        candidate_settings.more_scrutiny_coarser_FOM = best_settings.threshold_FOM + 1.451548;
        candidate_settings.more_scrutiny_min_dev_from_line = 5.866464;
        candidate_settings.amp_to_apply_line_test_below = 6;
        CandidatePeak_GA::eval_candidate_settings( candidate_settings, input_srcs, true );
      }//if( make N42 files )
      
      
      break;
    }//case OptimizationAction::AccuracyFromCsvsStudy:
      
    case OptimizationAction::CodeDev:
    {
      cerr << "Setting best_settings instead of actually finding them!" << endl;
      // Some non-optimal settings, to use for just development.
      best_settings.num_smooth_side_channels = 9;
      best_settings.smooth_polynomial_order = 2;
      best_settings.threshold_FOM = 0.914922;
      best_settings.more_scrutiny_FOM_threshold = best_settings.threshold_FOM + 0.643457;
      best_settings.pos_sum_threshold_sf = 0.048679;
      best_settings.num_chan_fluctuate = 1;
      best_settings.more_scrutiny_coarser_FOM = best_settings.threshold_FOM + 2.328440;
      best_settings.more_scrutiny_min_dev_from_line = 5.875653;
      best_settings.amp_to_apply_line_test_below = 6;
      
      if( false )
      {
        for( const DataSrcInfo &info : input_srcs )
        {
          const double start_wall = SpecUtils::get_wall_time();
          const double start_cpu = SpecUtils::get_cpu_time();
          
          for( size_t i = 0; i < 5000; ++i )
          {
            std::vector< std::tuple<float,float,float> > results;
            CandidatePeak_GA::find_candidate_peaks( info.src_info.src_spectra[0], 0, 0, results, best_settings );
          }
          
          const double end_wall = SpecUtils::get_wall_time();
          const double end_cpu = SpecUtils::get_cpu_time();
          
          cout << "Eval took " << (end_wall - start_wall) << " s, wall and " << (end_cpu - start_cpu) << " s cpu" << endl;
          exit( 0 );
        }
      }//if( benchmark find_candidate_peaks )
      
      if( false ) //benchmark peak-fitting
      {
        for( const DataSrcInfo &info : input_srcs )
        {
          const shared_ptr<const SpecUtils::Measurement> &data = info.src_info.src_spectra.front();
          
          vector<tuple<float,float,float>> dummy;
          vector<PeakDef> candidate_peaks = CandidatePeak_GA::find_candidate_peaks( data, 0, 0, dummy, best_settings );

          for( PeakDef &p : candidate_peaks )
          {
            p.continuum()->setType( PeakContinuum::OffsetType::FlatStep );
            //p.continuum()->setType( PeakContinuum::OffsetType::Linear );
          }

          const bool isHPGe = true, amplitudeOnly = false;
          vector<PeakDef> zeroth_fit_results, initial_fit_results;
          
          
          {//Begin optimization block - delta when done optimizing `fitPeaks(...)`
            const double start_wall = SpecUtils::get_wall_time();
            const double start_cpu = SpecUtils::get_cpu_time();
            
            //for( size_t i = 0; i < 1; ++i )
            for( size_t i = 0; i < 20; ++i )
              //while( true )
            {
              assert( 0 );
              
              // Check why one of these seems to be multithreaded, but not the other
              
              //We need to actually score using these two methods; using `fitPeaks(...)` takes
              //  {11.3/11.3}s, while `fitPeaksInRange(...)` takes {0.071/0.45}s.
              //  Although `fitPeaksInRange(...)` actually uses `fitPeaks(...)` to do the
              //  work - maybe fitPeaks doesnt separate ROIs???
              //  Also, need to compare to using `LinearProblemSubSolveChi2Fcn`, which is used by
              //  `refitPeaksThatShareROI(...)` and `refit_for_new_roi(...)` and `fit_peak_for_user_click(...)`
              //And also could try using Ceres to see if it works better than Minuit.
              vector<PeakDef> peaks;
#if( USE_CERES_PEAK_FITTING )
              vector<shared_ptr<const PeakDef>> results_tmp, input_peaks_tmp;
              for( const auto &p : candidate_peaks )
                input_peaks_tmp.push_back( make_shared<PeakDef>(p) );
              PeakFitLM::fit_peaks_LM( results_tmp, input_peaks_tmp, data, 0.0, 0.0, amplitudeOnly, isHPGe );
              for( const auto &p : results_tmp )
                peaks.push_back( *p );
#else
              fitPeaks( candidate_peaks, 0.0, 0.0, data, peaks, amplitudeOnly, isHPGe );
#endif

              //vector<PeakDef> peaksInRange = fitPeaksInRange( 0.0, data->gamma_energy_max(), 1.5, 0.0, 0.0, candidate_peaks, data, dummy_fixedpeaks, amplitudeOnly, isHPGe );
              
              //assert( peaksInRange.size() == peaks.size() );
              //for( size_t i = 0; i < peaksInRange.size(); ++i )
              //  cout << i << "\t" << peaks[i].mean() <<  ":\t"
              //  << candidate_peaks[i].amplitude() << "\t"
              //  << peaks[i].amplitude() << "\t" << peaksInRange[i].amplitude() << endl;
            }
            
            //Do this timing - and then disbale continua integral optimization and re-due
            const double end_wall = SpecUtils::get_wall_time();
            const double end_cpu = SpecUtils::get_cpu_time();
            
            cout << "Eval of step continua took " << (end_wall - start_wall) << " s, wall and " << (end_cpu - start_cpu) << " s cpu" << endl;
      exit(1);
          }//End optimization block
        }//for( const DataSrcInfo &info : input_srcs )
      }//if( benchmark peak-fitting )
      
      // Some real guess settings, to use just for development
      InitialPeakFindSettings fit_settings;
      fit_settings.initial_stat_threshold = 2.5; //  A reasonable search range of values is maybe between 0 and 8.
      fit_settings.initial_hypothesis_threshold = 0;  //A reasonable search range of values is maybe between -0.1 and 10.
      fit_settings.initial_min_nsigma_roi = 4;
      fit_settings.initial_max_nsigma_roi = 12;
      fit_settings.fwhm_fcn_form = InitialPeakFindSettings::FwhmFcnForm::SqrtPolynomialTwoCoefs;
      fit_settings.search_roi_nsigma_deficit = 2;//Reasonable search range of 1 to 10.
      fit_settings.search_stat_threshold = 3;//Reasonable search range of 0 and 8.
      fit_settings.search_hypothesis_threshold = 2;//Reasonable search range of 0 and 8. -0.1 and 10.
      fit_settings.search_stat_significance = 2.5;
      
      fit_settings.ROI_add_nsigma_required = 3;
      fit_settings.ROI_add_chi2dof_improve = 0.5;

      std::mutex score_mutex;
      double sum_find_weight = 0.0, sum_def_wanted_area_weight = 0.0, sum_maybe_wanted_area_weight = 0.0, sum_not_wanted_area_weight = 0.0;
      double sum_def_wanted_area_median_weight = 0.0, sum_maybe_wanted_area_median_weight = 0.0, sum_not_wanted_area_median_weight = 0.0;
      const double start_wall = SpecUtils::get_wall_time();
      const double start_cpu = SpecUtils::get_cpu_time();

      SpecUtilsAsync::ThreadPool pool;

      size_t num_posted = 0;
      for( const DataSrcInfo &info : input_srcs )
      {
        num_posted += 1;
        pool.post( [&](){
          if( PeakFitImprove::debug_printout )
          {
            std::lock_guard<std::mutex> lock( score_mutex );
            cout << "Evaluating " << info.location_name << "/" << info.live_time_name << "/" << info.detector_name << "/" << info.src_info.src_name << endl;
          }

          const InitialFit_GA::PeakFindAndFitWeights weight = InitialFit_GA::eval_initial_peak_find_and_fit( fit_settings, best_settings, info );

          std::lock_guard<std::mutex> lock( score_mutex );
          sum_find_weight += weight.find_weight;
          sum_def_wanted_area_weight += weight.def_wanted_area_weight;
          sum_maybe_wanted_area_weight += weight.maybe_wanted_area_weight;
          sum_not_wanted_area_weight += weight.not_wanted_area_weight;
          sum_def_wanted_area_median_weight += weight.def_wanted_area_median_weight;
          sum_maybe_wanted_area_median_weight += weight.maybe_wanted_area_median_weight;
          sum_not_wanted_area_median_weight += weight.not_wanted_area_median_weight;

          if( PeakFitImprove::debug_printout )
            cout << "For " << info.location_name << "/" << info.live_time_name << "/" << info.detector_name << "/" << info.src_info.src_name
            << ", got weight=" << sum_find_weight << endl;
        } );

        if( (num_posted % 100) == 0 )
        {
          pool.join();

          if( (num_posted % 1000) == 0 )
            cout << "Completed " << num_posted << " spectra of " << input_srcs.size() << endl;
        }//if( (num_posted % 50) == 0 )
      }//for( const DataSrcInfo &info : input_srcs )

      cout << "Have posted " << num_posted << " jobs - will wait on them to finish." << endl;
      pool.join();

      const double end_wall = SpecUtils::get_wall_time();
      const double end_cpu = SpecUtils::get_cpu_time();

      cout << "Average weights of " << input_srcs.size() << " files is:\n"
      << "\tFindWeight:          " << (sum_find_weight/input_srcs.size())              << endl
      << "\tDefWantAreaWeight:   " << (sum_def_wanted_area_weight/input_srcs.size())   << endl
      << "\tDefWantAreaMedian:   " << (sum_def_wanted_area_median_weight/input_srcs.size())   << endl
      << "\tMaybeWantAreaWeight: " << (sum_maybe_wanted_area_weight/input_srcs.size()) << endl
      << "\tMaybeWantAreaMedian: " << (sum_maybe_wanted_area_median_weight/input_srcs.size()) << endl
      << "\tDontWantAreaWeight:  " << (sum_not_wanted_area_weight/input_srcs.size())   << endl
      << "\tDontWantAreaMedian:  " << (sum_not_wanted_area_median_weight/input_srcs.size())   << endl
      << "Eval took " << (end_wall - start_wall) << " s, wall and " << (end_cpu - start_cpu) << " s cpu." << endl;

      break;
    }//case OptimizationAction::CodeDev:
  }//switch( action )

  
  return EXIT_SUCCESS;
}//int main( int argc, const char * argv[] )


