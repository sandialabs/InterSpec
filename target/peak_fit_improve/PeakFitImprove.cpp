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
#include <thread>

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

namespace PeakFitImprove
{
  size_t sm_num_optimization_threads = 8; // Will be set by command line or default calculation
  size_t sm_num_threads_per_individual = 1;
  size_t sm_ga_population = 1000;
  size_t sm_ga_generation_max = 100;
  size_t sm_ga_best_stall_max = 10;
  size_t sm_ga_elite_count = 10;
  double sm_ga_crossover_fraction = 0.7;
  double sm_ga_mutation_rate = 0.4;
  double sm_ga_mutate_threshold = 0.15;
  double sm_ga_crossover_threshold = -1.0; // Will be set based on action if not specified
}



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
  //+ var_name + ".search_stat_threshold = "        + std::to_string(search_stat_threshold)        + ";\n"
  + var_name + ".search_hypothesis_threshold = "  + std::to_string(search_hypothesis_threshold)  + ";\n"
  + var_name + ".search_stat_significance = "     + std::to_string(search_stat_significance)     + ";\n"
  + var_name + ".ROI_add_nsigma_required = "      + std::to_string(ROI_add_nsigma_required)      + ";\n"
  + var_name + ".ROI_add_chi2dof_improve = "      + std::to_string(ROI_add_chi2dof_improve)      + ";\n"
  ;
}//std::string print( const string &var_name ) const
  

std::string FinalPeakFitSettings::print( const string &var_name ) const
{
  return var_name + ".require_combine_num_fwhm_near = "          + std::to_string(require_combine_num_fwhm_near)               + ";\n"
  + var_name + ".not_allow_combine_num_fwhm_near = "             + std::to_string(not_allow_combine_num_fwhm_near)             + ";\n"
 // + var_name + ".combine_ROI_overlap_frac = "                    + std::to_string(combine_ROI_overlap_frac)                    + ";\n"
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
  + var_name + ".roi_extent_low_num_fwhm_base_highstat = "       + std::to_string(roi_extent_low_num_fwhm_base_highstat)       + ";\n"
  + var_name + ".roi_extent_high_num_fwhm_base_highstat = "      + std::to_string(roi_extent_high_num_fwhm_base_highstat)      + ";\n"
  + var_name + ".roi_extent_low_num_fwhm_base_lowstat = "        + std::to_string(roi_extent_low_num_fwhm_base_lowstat)        + ";\n"
  + var_name + ".roi_extent_high_num_fwhm_base_lowstat = "       + std::to_string(roi_extent_high_num_fwhm_base_lowstat)       + ";\n"
  + var_name + ".high_stat_threshold = "               + std::to_string(high_stat_threshold)                                   + ";\n"
  + var_name + ".roi_extent_low_num_fwhm_extra = "               + std::to_string(roi_extent_low_num_fwhm_extra)               + ";\n"
  + var_name + ".roi_extent_high_num_fwhm_extra = "              + std::to_string(roi_extent_high_num_fwhm_extra)              + ";\n"
  + var_name + ".roi_end_second_deriv_thresh = "                 + std::to_string(roi_end_second_deriv_thresh)                 + ";\n"
  + var_name + ".break_multi_roi_up_continuum_away_sigma = "     + std::to_string(break_multi_roi_up_continuum_away_sigma)     + ";\n"
  + var_name + ".break_multi_roi_up_required_chi2dof_improve = "     + std::to_string(break_multi_roi_up_required_chi2dof_improve)     + ";\n"
  ;
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


void create_n42_peak_fits_for_dir( const string &dir )
{
  const bool isHPGe = true;

  //Generation [48], Best=-19.7537, Average=-19.7535, Best genes:
  //Using all live time and HPGe detectors, and cities - 20250912
  const FindCandidateSettings hpge_candidate_settings = ([](){
    FindCandidateSettings settings;
    settings.num_smooth_side_channels = 3;
    settings.smooth_polynomial_order = 2;
    settings.threshold_FOM = 1.25;
    settings.more_scrutiny_FOM_threshold = 2.5;
    settings.pos_sum_threshold_sf = 0.119178;
    settings.num_chan_fluctuate = 1;
    settings.more_scrutiny_coarser_FOM = 3.001943;
    settings.more_scrutiny_min_dev_from_line = 6.816465;
    settings.amp_to_apply_line_test_below = 6.000000;
    return settings;
  })();

  const InitialPeakFindSettings hpge_initial_fit_settings = ([](){
    InitialPeakFindSettings settings;
    settings.initial_stat_threshold = 3.5;
    settings.initial_hypothesis_threshold = 0.5;
    settings.initial_min_nsigma_roi = 2.246770;
    settings.initial_max_nsigma_roi = 6.378162;
    settings.fwhm_fcn_form = InitialPeakFindSettings::FwhmFcnForm::SqrtPolynomialTwoCoefs;
    settings.search_roi_nsigma_deficit = 4.241748;
    //settings.search_stat_threshold = 8.051485;
    settings.search_hypothesis_threshold= 3.342207;
    settings.search_stat_significance = 2.025582;
    settings.ROI_add_nsigma_required = 3.526017;
    settings.ROI_add_chi2dof_improve = 0.45;
    return settings;
  })();

  // Not final settings
  const FinalPeakFitSettings hpge_final_fit_settings = ([](){
    FinalPeakFitSettings settings;

    settings.require_combine_num_fwhm_near = 2.0;
    settings.not_allow_combine_num_fwhm_near = 4.0;
    settings.cont_type_peak_nsigma_threshold = 40.450658;
    settings.cont_type_left_right_nsigma = 8.423588;
    settings.cont_poly_order_increase_chi2dof_required = 0.671733;
    settings.cont_step_type_increase_chi2dof_required = 0.312514;
    settings.skew_nsigma = 8.359958;
    settings.left_residual_sum_min_to_try_skew = 2.878255;
    settings.right_residual_sum_min_to_try_skew = 0.158640;
    settings.skew_improve_chi2_dof_threshold = 1.016935;
    settings.roi_extent_low_num_fwhm_base_highstat = 1.5;
    settings.roi_extent_high_num_fwhm_base_highstat = 1.5;
    settings.roi_extent_low_num_fwhm_base_lowstat = 1.5;
    settings.roi_extent_high_num_fwhm_base_lowstat = 1.5;
    settings.high_stat_threshold = 15.0;
    settings.roi_extent_low_num_fwhm_extra = 1.5;
    settings.roi_extent_high_num_fwhm_extra = 1.5;
    settings.roi_end_second_deriv_thresh = 3.5;
    settings.break_multi_roi_up_continuum_away_sigma = 5;
    settings.break_multi_roi_up_required_chi2dof_improve = 0.5;

    return settings;
  })();



  auto do_fit_peaks = [&hpge_candidate_settings, &hpge_initial_fit_settings, &hpge_final_fit_settings, isHPGe]
                        ( const shared_ptr<const SpecUtils::Measurement> &spectrum )
                      -> vector<PeakDef>
  {
    const bool multithread = true;
    size_t num_add_candidates_fit_for = 0, num_add_candidates_accepted = 0; //Only for eval purposes
    const vector<PeakDef> initial_peaks = InitialFit_GA::initial_peak_find_and_fit( hpge_initial_fit_settings,
                                                                                   hpge_candidate_settings,
                                                                                   spectrum,
                                                                                   multithread,
                                                                                   num_add_candidates_fit_for,
                                                                                   num_add_candidates_accepted);
    //return initial_peaks;

    cout << "Input into final_peak_fit:" << endl;
    for( const PeakDef &p : initial_peaks )
      cout << "    {" << p.mean() << ", " << p.fwhm() << ", " << p.lowerX() << "-" << p.upperX() << "}" << endl;
    cout << endl;

    const vector<PeakDef> fit_peaks = FinalFit_GA::final_peak_fit( initial_peaks, hpge_final_fit_settings, isHPGe,
                                                                  spectrum, multithread );

    return fit_peaks;
  };//do_fit_peaks


  ofstream output( dir + "/peak_fits.html" );

  // Lu177m_Unsh is probably the best test of settings...

  D3SpectrumExport::write_html_page_header( output, "Peak Fits", "InterSpec_resources" );

  output << "<body>" << endl;

  //Hack putting <style></style> block in HTML body, but wahtever, seems ot work
  output <<"<style>"
  << ".TopLinesTable{ margin-top: 25px; margin-left: auto; margin-right: auto; border-collapse: collapse; border: 1px solid black; }" << endl
  << "table, th, td{ border: 1px solid black; }" << endl
  << "fieldset{width: 90vw; margin-left: auto; margin-right: auto; margin-top: 20px;}" << endl
  << "</style>" << endl;

  vector<string> files = SpecUtils::recursive_ls( dir, ".csv" );

  for( size_t i = 0; i < files.size(); ++i )
  {
    const string filename = files[i];

    if( SpecUtils::icontains(filename, ".peaks.CSV") )
      continue;

    //if( !SpecUtils::icontains(filename, "spec_8h.csv") )
    //  continue;


    SpecMeas meas;
    const bool loaded = meas.load_file(filename, SpecUtils::ParserType::Auto );
    if( !loaded || (meas.num_measurements() != 1) )
    {
      cerr << "Failed to load '" << filename << "'" << endl;
      exit(1);
    }

    const shared_ptr<const SpecUtils::Measurement> spectrum_orig = meas.measurement_at_index( 0 );
    assert( spectrum_orig );
    //meas.set_live_time( 1.0, spectrum );
    //meas.set_real_time( 1.0, spectrum );

    const shared_ptr<const vector<float>> &counts_ptr = spectrum_orig->gamma_channel_contents();
    const vector<float> &counts_orig = *counts_ptr;
    float min_counts = 1.0E8;
    for( const float &c : counts_orig )
    {
      if( c > 0.5f )
        min_counts = std::min( min_counts, c );
    }
    shared_ptr<vector<float>> scaled_cnts = make_unique<vector<float>>( counts_orig );
    for( float &c : *scaled_cnts )
      c /= min_counts;

    shared_ptr<SpecUtils::Measurement> spectrum = make_shared<SpecUtils::Measurement>( *spectrum_orig );
    spectrum->set_gamma_counts( scaled_cnts, 1.0f, 1.0f );

    meas.remove_measurement( spectrum_orig, false );
    meas.add_measurement( spectrum, true );


    const vector<PeakDef> fit_peaks = do_fit_peaks( spectrum );
    vector<shared_ptr<const PeakDef>> fit_peaks_ptrs;
    std::deque<std::shared_ptr<const PeakDef>> peaks_dq;
    for( const PeakDef &p : fit_peaks )
    {
      auto np = make_shared<PeakDef>(p);
      fit_peaks_ptrs.push_back( np );
      peaks_dq.push_back( np );
    }


    meas.setPeaks( peaks_dq, {spectrum->sample_number()} );

    meas.save2012N42File( filename + ".n42" );


    {
      ofstream peak_out( (filename + ".peaks.CSV").c_str(), ios::out | ios::binary );
      PeakModel::write_peak_csv( peak_out, SpecUtils::filename(filename), PeakModel::PeakCsvType::Full, peaks_dq, spectrum );
    }


    string ref_line_json;
    std::string title = "";
    std::string dataTitle = "";
    bool useLogYAxis = true, showVerticalGridLines = false, showHorizontalGridLines = false;
    bool legendEnabled = true, compactXAxis = true;
    bool showPeakUserLabels = false, showPeakEnergyLabels = false, showPeakNuclideLabels = false, showPeakNuclideEnergyLabels = false;
    bool showEscapePeakMarker = false, showComptonPeakMarker = false, showComptonEdgeMarker = false, showSumPeakMarker = false;
    bool backgroundSubtract = false;
    float xMin = 0, xMax = 3000;
    std::map<std::string,std::string> refernce_lines_json;
    refernce_lines_json["TruthPeaks"] = ref_line_json;

    D3SpectrumExport::D3SpectrumChartOptions options( title, "Energy (keV)", "Counts/Channel",
                                                     dataTitle, useLogYAxis,
                                                     showVerticalGridLines, showHorizontalGridLines,
                                                     legendEnabled, compactXAxis,
                                                     showPeakUserLabels, showPeakEnergyLabels, showPeakNuclideLabels,
                                                     showPeakNuclideEnergyLabels, showEscapePeakMarker, showComptonPeakMarker,
                                                     showComptonEdgeMarker, showSumPeakMarker, backgroundSubtract,
                                                     xMin, xMax, refernce_lines_json );

    D3SpectrumExport::D3SpectrumOptions foreground_opts;
    foreground_opts.line_color = "black";
    foreground_opts.title = SpecUtils::filename(filename);
    foreground_opts.display_scale_factor = 1.0;
    foreground_opts.spectrum_type = SpecUtils::SpectrumType::Foreground;
    foreground_opts.peaks_json = PeakDef::peak_json( fit_peaks_ptrs, spectrum, Wt::WColor(), false );

    const string div_id = "chart_" + std::to_string(i);


    output << "<fieldset style=\"\">" << endl
    << "<legend>" << SpecUtils::filename(filename) << "</legend>" << endl;

    output << "<div id=\"" << div_id << "\" class=\"chart\" oncontextmenu=\"return false;\"></div>" << endl;  // Adding the main chart div


    output << "<script>" << endl;

    D3SpectrumExport::write_js_for_chart( output, div_id, options.m_dataTitle, options.m_xAxisTitle, options.m_yAxisTitle );

    std::vector< std::pair<const SpecUtils::Measurement *,D3SpectrumExport::D3SpectrumOptions> > measurements;
    measurements.emplace_back( spectrum.get(), foreground_opts );

    write_and_set_data_for_chart( output, div_id, measurements );


    output << R"delim(
    const resizeChart)delim" << i << R"delim( = function(){
      let height = window.innerHeight;
      let width = document.documentElement.clientWidth;
      let el = spec_chart_)delim" << div_id << R"delim(.chart;
      el.style.width = 0.8*width + "px";
      el.style.height = Math.min(500,Math.max(250, Math.min(0.4*width,height-175))) + "px";
      el.style.marginLeft = 0.05*width + "px";
      el.style.marginRight = 0.05*width + "px";
      
    )delim"
    << "  spec_chart_" << div_id << R"delim(.handleResize();
    };
    
    window.addEventListener('resize', resizeChart)delim" << i << R"delim();
    )delim" << endl;

    write_set_options_for_chart( output, div_id, options );

    output << "spec_chart_" << div_id << ".setReferenceLines( reference_lines_" << div_id << " );" << endl;

    output << "resizeChart" << i << "();" << endl;
    output << "</script>" << endl;



    output << "<table class=\"TopLinesTable\" style=\"\">" << endl;
    output << "<tr><th>Fit Energy (keV)</th><th>Fit Area</th><th>FitAreaUncert</th><th>ROI Lower</th><th>ROI Upper</th></tr>" << endl;
    output << "<caption>Top gamma lines in spectrum</caption>" << endl;
    for( size_t i = 0; i < fit_peaks.size(); ++i )
    {
      const PeakDef &p = fit_peaks[i];

      output << "<tr>"
      << "<td>" << p.mean() << "</td>"
      << "<td>" << p.amplitude() << "</td>"
      << "<td>" << p.amplitudeUncert() << "</td>"
      << "<td>" << p.lowerX() << "</td>"
      << "<td>" << p.upperX() << "</td>"
      << "</tr>"
      << endl;
    }

    output << "</table>" << endl;

    output << "</fieldset>" << endl;
  }//for( size_t i = 0; i < src_info.size(); ++i )

  output << "</body>" << endl;
  output << "</html>" << endl;

  cout << "Done" << endl;
}//void create_n42_peak_fits_for_dir()


void create_n42_peak_fits( const vector<DetectorInjectSet> &inject_sets, const vector<DataSrcInfo> &input_srcs )
{
  const bool isHPGe = true;

  //Generation [48], Best=-19.7537, Average=-19.7535, Best genes:
  //Using all live time and HPGe detectors, and cities - 20250912
  const FindCandidateSettings hpge_candidate_settings = ([](){
    FindCandidateSettings settings;
    settings.num_smooth_side_channels = 9;
    settings.smooth_polynomial_order = 2;
    settings.threshold_FOM = 0.758621;
    settings.more_scrutiny_FOM_threshold = 1.598265;
    settings.pos_sum_threshold_sf = 0.119178;
    settings.num_chan_fluctuate = 1;
    settings.more_scrutiny_coarser_FOM = 3.001943;
    settings.more_scrutiny_min_dev_from_line = 6.816465;
    settings.amp_to_apply_line_test_below = 6.000000;
    return settings;
  })();

  const InitialPeakFindSettings hpge_initial_fit_settings = ([](){;
    //Generation [138], Best=-36152.8, Average=-36151.4, - from optimization in Sep. 2025
    InitialPeakFindSettings settings;
    settings.initial_stat_threshold = 1.951264;
    settings.initial_hypothesis_threshold = 0.673169;
    settings.initial_min_nsigma_roi = 2.246770;
    settings.initial_max_nsigma_roi = 6.378162;
    settings.fwhm_fcn_form = InitialPeakFindSettings::FwhmFcnForm::SqrtPolynomialTwoCoefs;
    settings.search_roi_nsigma_deficit = 4.241748;
    //settings.search_stat_threshold = 8.051485;
    settings.search_hypothesis_threshold= 3.342207;
    settings.search_stat_significance = 3.0;
    settings.ROI_add_nsigma_required = 3.526017;
    settings.ROI_add_chi2dof_improve = 0.45;
    return settings;
  })();

  // Not final settings
  const FinalPeakFitSettings hpge_final_fit_settings = ([](){
    FinalPeakFitSettings settings;

    settings.require_combine_num_fwhm_near = 2.0;
    settings.not_allow_combine_num_fwhm_near = 4.0;
    settings.cont_type_peak_nsigma_threshold = 40.450658;
    settings.cont_type_left_right_nsigma = 8.423588;
    settings.cont_poly_order_increase_chi2dof_required = 0.671733;
    settings.cont_step_type_increase_chi2dof_required = 0.6;
    settings.skew_nsigma = 8.359958;
    settings.left_residual_sum_min_to_try_skew = 2.878255;
    settings.right_residual_sum_min_to_try_skew = 0.158640;
    settings.skew_improve_chi2_dof_threshold = 1.016935;
    settings.roi_extent_low_num_fwhm_base_highstat = 1.5;
    settings.roi_extent_high_num_fwhm_base_highstat = 1.5;
    settings.roi_extent_low_num_fwhm_base_lowstat = 1.5;
    settings.roi_extent_high_num_fwhm_base_lowstat = 1.5;
    settings.high_stat_threshold = 15.0;
    settings.roi_extent_low_num_fwhm_extra = 1.5;
    settings.roi_extent_high_num_fwhm_extra = 1.5;
    settings.roi_end_second_deriv_thresh = 3.5;
    settings.break_multi_roi_up_continuum_away_sigma = 5;
    settings.break_multi_roi_up_required_chi2dof_improve = 0.5;

    return settings;
  })();

  FinalFit_GA::FinalFitScore sum_score;


  auto do_fit_peaks = [&hpge_candidate_settings, &hpge_initial_fit_settings, &hpge_final_fit_settings, isHPGe, &sum_score]
                        ( const shared_ptr<const SpecUtils::Measurement> &spectrum, const DataSrcInfo &info )
                      -> vector<PeakDef>
  {
    const bool multithread = true;
    size_t num_add_candidates_fit_for = 0, num_add_candidates_accepted = 0; //Only for eval purposes
    const vector<PeakDef> initial_peaks = InitialFit_GA::initial_peak_find_and_fit( hpge_initial_fit_settings,
                                                                                   hpge_candidate_settings,
                                                                                   spectrum,
                                                                                   multithread,
                                                                                   num_add_candidates_fit_for,
                                                                                   num_add_candidates_accepted);
    //return initial_peaks;

    cout << "Input into final_peak_fit:" << endl;
    for( const PeakDef &p : initial_peaks )
      cout << "    {" << p.mean() << ", " << p.fwhm() << ", " << p.lowerX() << "-" << p.upperX() << "}" << endl;
    cout << endl;

    const vector<PeakDef> fit_peaks = FinalFit_GA::final_peak_fit( initial_peaks, hpge_final_fit_settings, isHPGe,
                                                                  spectrum, multithread );

    
    const FinalFit_GA::FinalFitScore score = FinalFit_GA::eval_final_peak_fit( hpge_final_fit_settings, info, initial_peaks, false );


    sum_score.area_score += score.area_score;
    sum_score.width_score += score.width_score;
    sum_score.position_score += score.position_score;
    sum_score.ignored_unexpected_peaks += score.ignored_unexpected_peaks;
    sum_score.unexpected_peaks_sum_significance += score.unexpected_peaks_sum_significance;
    sum_score.total_weight += score.total_weight;
    sum_score.num_peaks_used += score.num_peaks_used;

    //cout << info.location_name << "/"<< info.detector_name << "/" << info.live_time_name << "/"
    //<< info.src_info.src_name << "\n-> " << score.print( "score" ) << "--------" << endl << endl;


    return fit_peaks;
  };//do_fit_peaks


  ofstream output( "peak_fits.html" );

  // Lu177m_Unsh is probably the best test of settings...

  D3SpectrumExport::write_html_page_header( output, "Peak Fits", "InterSpec_resources" );

  output << "<body>" << endl;

  //Hack putting <style></style> block in HTML body, but wahtever, seems ot work
  output <<"<style>"
  << ".TopLinesTable{ margin-top: 25px; margin-left: auto; margin-right: auto; border-collapse: collapse; border: 1px solid black; }" << endl
  << "table, th, td{ border: 1px solid black; }" << endl
  << "fieldset{width: 90vw; margin-left: auto; margin-right: auto; margin-top: 20px;}" << endl
  << "</style>" << endl;


  for( size_t i = 0; i < input_srcs.size(); ++i )
  {
    const DataSrcInfo &info = input_srcs[i];
    const InjectSourceInfo &src_info = info.src_info;

    const string &detector_name = info.detector_name;
    const string &location_name = info.location_name;
    const string &live_time_name = info.live_time_name;

    //if( !SpecUtils::icontains( info.src_info.src_name, "Lu177m_Unsh") )
    //  continue;


    const vector<PeakTruthInfo> &source_lines = src_info.source_lines;
    const vector<PeakTruthInfo> &background_lines = src_info.background_lines;
    //src_info.spec_file;
    const vector<shared_ptr<const SpecUtils::Measurement>> &src_spectra = src_info.src_spectra;
    assert( !src_spectra.empty() );
    const shared_ptr<const SpecUtils::Measurement> &spectrum = src_spectra.front(); //We'll just plot the first spectrum only
    const shared_ptr<const SpecUtils::Measurement> &short_background = src_info.short_background;
    const shared_ptr<const SpecUtils::Measurement> &long_background = src_info.long_background;
    //src_info.src_no_poisson;
    //src_info.background_no_poisson;

    const vector<ExpectedPhotopeakInfo> &expected_photopeaks = info.expected_photopeaks;

    const vector<PeakDef> fit_peaks = do_fit_peaks( spectrum, info );
    vector<shared_ptr<const PeakDef>> fit_peaks_ptrs;
    for( const PeakDef &p : fit_peaks )
      fit_peaks_ptrs.push_back( make_shared<PeakDef>(p) );

    ReferenceLineInfo ref_line_info;
    ref_line_info.m_validity = ReferenceLineInfo::InputValidity::Valid;
    ref_line_info.m_source_type = ReferenceLineInfo::SourceType::OneOffSrcLines;
    ref_line_info.m_has_coincidences = false;

    double max_peak_area = 0.0;
    for( const ExpectedPhotopeakInfo &line : expected_photopeaks )
    {
      if( line.nsigma_over_background < 0.5 ) //Dont show 1-sigma peaks
        continue;

      ReferenceLineInfo::RefLine ref_line;
      ref_line.m_energy = line.effective_energy;
      ref_line.m_normalized_intensity = line.peak_area;
      max_peak_area = std::max( max_peak_area, line.peak_area );
      ref_line.m_drf_factor = 1.0;
      ref_line.m_shield_atten = 1.0;
      ref_line.m_particle_sf_applied = 1.0;
      ref_line.m_color = Wt::WColor( Wt::GlobalColor::darkBlue );
      ref_line.m_decay_intensity = line.peak_area;
      ref_line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
      ref_line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Normal;
      ref_line.m_attenuation_applies = false;

      ref_line_info.m_ref_lines.push_back( ref_line );
    }//for( const ExpectedPhotopeakInfo &line : expected_photopeaks )

    if( max_peak_area > 1.0 )
    {
      for( auto &p : ref_line_info.m_ref_lines )
        p.m_normalized_intensity /= max_peak_area;
    }

    string ref_line_json;
    ref_line_info.toJson( ref_line_json );

    std::string title = "";
    std::string dataTitle = "";
    bool useLogYAxis = true, showVerticalGridLines = false, showHorizontalGridLines = false;
    bool legendEnabled = true, compactXAxis = true;
    bool showPeakUserLabels = false, showPeakEnergyLabels = false, showPeakNuclideLabels = false, showPeakNuclideEnergyLabels = false;
    bool showEscapePeakMarker = false, showComptonPeakMarker = false, showComptonEdgeMarker = false, showSumPeakMarker = false;
    bool backgroundSubtract = false;
    float xMin = 0, xMax = 3000;
    std::map<std::string,std::string> refernce_lines_json;
    refernce_lines_json["TruthPeaks"] = ref_line_json;

    D3SpectrumExport::D3SpectrumChartOptions options( title, "Energy (keV)", "Counts/Channel",
                                                     dataTitle, useLogYAxis,
                                                     showVerticalGridLines, showHorizontalGridLines,
                                                     legendEnabled, compactXAxis,
                                                     showPeakUserLabels, showPeakEnergyLabels, showPeakNuclideLabels,
                                                     showPeakNuclideEnergyLabels, showEscapePeakMarker, showComptonPeakMarker,
                                                     showComptonEdgeMarker, showSumPeakMarker, backgroundSubtract,
                                                     xMin, xMax, refernce_lines_json );

    D3SpectrumExport::D3SpectrumOptions foreground_opts, background_opts;
    foreground_opts.line_color = "black";
    background_opts.line_color = "steelblue";
    foreground_opts.title = src_info.src_name;
    background_opts.title = "Background";
    foreground_opts.display_scale_factor = 1.0;
    background_opts.display_scale_factor = spectrum->live_time() / long_background->live_time();
    foreground_opts.spectrum_type = SpecUtils::SpectrumType::Foreground;
    background_opts.spectrum_type = SpecUtils::SpectrumType::Background;
    foreground_opts.peaks_json = PeakDef::peak_json( fit_peaks_ptrs, spectrum, Wt::WColor(), false );

    const string div_id = "chart_" + std::to_string(i);


    output << "<fieldset style=\"\">" << endl
    << "<legend>" << info.location_name << "/" << info.detector_name << "/" << info.live_time_name << "/" << src_info.src_name << "</legend>" << endl;

    output << "<div id=\"" << div_id << "\" class=\"chart\" oncontextmenu=\"return false;\"></div>" << endl;  // Adding the main chart div


    output << "<script>" << endl;

    D3SpectrumExport::write_js_for_chart( output, div_id, options.m_dataTitle, options.m_xAxisTitle, options.m_yAxisTitle );

    std::vector< std::pair<const SpecUtils::Measurement *,D3SpectrumExport::D3SpectrumOptions> > measurements;
    measurements.emplace_back( spectrum.get(), foreground_opts );
    measurements.emplace_back( long_background.get(), background_opts );

    write_and_set_data_for_chart( output, div_id, measurements );


    output << R"delim(
    const resizeChart)delim" << i << R"delim( = function(){
      let height = window.innerHeight;
      let width = document.documentElement.clientWidth;
      let el = spec_chart_)delim" << div_id << R"delim(.chart;
      el.style.width = 0.8*width + "px";
      el.style.height = Math.min(500,Math.max(250, Math.min(0.4*width,height-175))) + "px";
      el.style.marginLeft = 0.05*width + "px";
      el.style.marginRight = 0.05*width + "px";
      
    )delim"
    << "  spec_chart_" << div_id << R"delim(.handleResize();
    };
    
    window.addEventListener('resize', resizeChart)delim" << i << R"delim();
    )delim" << endl;

    write_set_options_for_chart( output, div_id, options );

    output << "spec_chart_" << div_id << ".setReferenceLines( reference_lines_" << div_id << " );" << endl;

    output << "resizeChart" << i << "();" << endl;
    output << "</script>" << endl;

    vector<ExpectedPhotopeakInfo> photopeaks = expected_photopeaks;
    std::sort( begin(photopeaks), end(photopeaks), [](const ExpectedPhotopeakInfo &lhs, const ExpectedPhotopeakInfo &rhs ){
      return lhs.peak_area > rhs.peak_area;
    });

    output << "<table class=\"TopLinesTable\" style=\"\">" << endl;
    output << "<tr><th>Energy (keV)</th><th>Fit Mean</th><th>Truth Area</th><th>Fit Area</th><th>Truth CPS</th><th>Fit CPS</th><th>FitAreaUncert</th><th># sigma off</th><th>ROI Lower</th><th>ROI Upper</th><th>Continuum Area</th></tr>" << endl;
    output << "<caption>Top gamma lines in spectrum</caption>" << endl;
    for( size_t i = 0; (i < photopeaks.size()) && (i < 20); ++i )
    {
      const ExpectedPhotopeakInfo &p = photopeaks[i];


      double nearest_energy = 1.0E6;
      const PeakDef *nearest_fit = nullptr;
      for( const PeakDef &peak : fit_peaks )
      {
        const double de = fabs(p.effective_energy - peak.mean());
        if( (de < p.effective_fwhm) && (de < nearest_energy) )
        {
          nearest_energy = de;
          nearest_fit = &peak;
        }
      }


      output << "<tr>"
      << "<td>" << p.effective_energy << "</td>"
      << "<td>" << (nearest_fit ? nearest_fit->mean() : 0.0) << "</td>"
      << "<td>" << p.peak_area << "</td>"
      << "<td>" << (nearest_fit ? nearest_fit->peakArea() : -999.9) << "</td>"
      << "<td>" << p.peak_area / spectrum->live_time() << "</td>"
      << "<td>" << (nearest_fit ? (nearest_fit->peakArea() / spectrum->live_time()) : -999.9) << "</td>"
      << "<td>" << (nearest_fit ? nearest_fit->peakAreaUncert() : -999.9) << "</td>"
      << "<td>" << (nearest_fit ? (fabs(nearest_fit->peakArea() - p.peak_area) / nearest_fit->peakAreaUncert()) : -999.9) << "</td>"
      << "<td>" << p.roi_lower << "</td>"
      << "<td>" << p.roi_upper << "</td>"
      << "<td>" << p.continuum_area << "</td>"
      << "</tr>"
      << endl;
    }

    if( photopeaks.size() > 20 )
      output << "<tr><td colspan=\"6\">Plus " << (photopeaks.size() - 20) << " more peaks</td></tr>" << endl;

    output << "</table>" << endl;

    output << "</fieldset>" << endl;
  }//for( size_t i = 0; i < src_info.size(); ++i )

  output << "</body>" << endl;
  output << "</html>" << endl;


  cout << sum_score.print( "sum_score" ) << "--------" << endl << endl;
}//void create_n42_peak_fits()




int main( int argc, char **argv )
{
  const double start_wall = SpecUtils::get_wall_time();
  const double start_cpu = SpecUtils::get_cpu_time();
  
  // Command line argument parsing
  namespace po = boost::program_options;
  
  string data_base_dir = "/Users/wcjohns/rad_ana/peak_area_optimization/peak_fit_accuracy_inject/";
  string static_data_dir;
  PeakFitImprove::sm_num_optimization_threads = std::max( 8u, std::thread::hardware_concurrency() > 2 ? std::thread::hardware_concurrency() - 2 : 1 );
  size_t number_threads_per_individual = 1;
  string action_str = "FinalFit"; //"CodeDev"; // "FinalFit"; //"InitialFit"; //"Candidate";
  bool debug_printout_arg = false;
  size_t ga_population = 1500;
  size_t ga_generation_max = 250;
  size_t ga_best_stall_max = 15;
  size_t ga_elite_count = 10;
  double ga_crossover_fraction = 0.7;
  double ga_mutation_rate = 0.4;
  double ga_mutate_threshold = 0.15;
  double ga_crossover_threshold = -1.0; // Will be set based on action if not specified
  
  po::options_description desc( "Allowed options" );
  desc.add_options()
    ("help,h", "produce help message")
    ("data-base-dir", po::value<string>( &data_base_dir )->default_value( data_base_dir ), "base directory for input data")
    ("number-threads", po::value<size_t>( &PeakFitImprove::sm_num_optimization_threads )->default_value( PeakFitImprove::sm_num_optimization_threads ), "number of individuals evaluated at the same time")
    ("num-threads-per-individual", po::value<size_t>( &PeakFitImprove::sm_num_threads_per_individual )->default_value( PeakFitImprove::sm_num_threads_per_individual ), "number of threads each individual can use")
    ("action", po::value<string>( &action_str )->default_value( action_str ), "optimization action: Candidate, InitialFit, FinalFit, CodeDev, AccuracyFromCsvsStudy")
    ("debug-printout", po::bool_switch( &debug_printout_arg ), "enable debug printout")
    ("static-data-dir", po::value<string>( &static_data_dir ), "static data directory (optional)")
    ("ga-population", po::value<size_t>( &ga_population )->default_value( ga_population ), "genetic algorithm population size")
    ("ga-generation-max", po::value<size_t>( &ga_generation_max )->default_value( ga_generation_max ), "genetic algorithm maximum generations")
    ("ga-best-stall-max", po::value<size_t>( &ga_best_stall_max )->default_value( ga_best_stall_max ), "genetic algorithm best stall maximum")
    ("ga-elite-count", po::value<size_t>( &ga_elite_count )->default_value( ga_elite_count ), "genetic algorithm elite count")
    ("ga-crossover-fraction", po::value<double>( &ga_crossover_fraction )->default_value( ga_crossover_fraction ), "Fraction of individuals that will have crossover applied to them every generation")
    ("ga-mutation-rate", po::value<double>( &ga_mutation_rate )->default_value( ga_mutation_rate ), "Fraction of individuals that will have mutation applied to them every generation")
    ("ga-mutate-threshold", po::value<double>( &ga_mutate_threshold )->default_value( ga_mutate_threshold ), "The fraction of genes that will be mutated.")
    ("ga-crossover-threshold", po::value<double>( &ga_crossover_threshold ), "The fraction of genes that will be crossed over (default depends on action: Candidate/InitialFit=0.25, FinalFit=0.15)");
  
  po::variables_map vm;
  
  try
  {
    po::store( po::parse_command_line( argc, argv, desc ), vm );
    po::notify( vm );
  }catch( const po::error &e )
  {
    cerr << "Error parsing command line: " << e.what() << endl;
    return -1;
  }
  
  if( vm.count( "help" ) )
  {
    cout << desc << endl;
    return 0;
  }
  
  // Validate and set up directories
  if( !SpecUtils::is_directory( data_base_dir ) )
  {
    cerr << "Error: data-base-dir '" << data_base_dir << "' does not exist or is not a directory." << endl;
    return -2;
  }
  
  if( !static_data_dir.empty() && !SpecUtils::is_directory( static_data_dir ) )
  {
    cerr << "Error: static-data-dir '" << static_data_dir << "' does not exist or is not a directory." << endl;
    return -3;
  }
  
  // Set crossover_threshold default based on action if not specified by user
  if( ga_crossover_threshold < 0.0 )
  {
    if( action_str == "FinalFit" )
      ga_crossover_threshold = 0.15;
    else  // Candidate, InitialFit, CodeDev, AccuracyFromCsvsStudy
      ga_crossover_threshold = 0.25;
  }
  
  // Set the optimization configuration
  PeakFitImprove::sm_ga_population = ga_population;
  PeakFitImprove::sm_ga_generation_max = ga_generation_max;
  PeakFitImprove::sm_ga_best_stall_max = ga_best_stall_max;
  PeakFitImprove::sm_ga_elite_count = ga_elite_count;
  PeakFitImprove::sm_ga_crossover_fraction = ga_crossover_fraction;
  PeakFitImprove::sm_ga_mutation_rate = ga_mutation_rate;
  PeakFitImprove::sm_ga_mutate_threshold = ga_mutate_threshold;
  PeakFitImprove::sm_ga_crossover_threshold = ga_crossover_threshold;
  
  // Parse action enum
  enum class OptimizationAction : int
  {
    Candidate,
    InitialFit,
    FinalFit,
    CodeDev,
    AccuracyFromCsvsStudy
  };//enum class OptimizationAction : int
  
  OptimizationAction action;
  if( action_str == "Candidate" )
    action = OptimizationAction::Candidate;
  else if( action_str == "InitialFit" )
    action = OptimizationAction::InitialFit;
  else if( action_str == "FinalFit" )
    action = OptimizationAction::FinalFit;
  else if( action_str == "CodeDev" )
    action = OptimizationAction::CodeDev;
  else if( action_str == "AccuracyFromCsvsStudy" )
    action = OptimizationAction::AccuracyFromCsvsStudy;
  else
  {
    cerr << "Error: invalid action '" << action_str << "'. Valid actions are: Candidate, InitialFit, FinalFit, CodeDev, AccuracyFromCsvsStudy" << endl;
    return -4;
  }
  
  cout << "Using " << PeakFitImprove::sm_num_optimization_threads
  << " threads for individual evaluation, with each individual using up to " << PeakFitImprove::sm_num_threads_per_individual
  << " threads for optimization."
  << endl;
  cout << "Data base directory: " << data_base_dir << endl;
  if( !static_data_dir.empty() )
    cout << "Static data directory: " << static_data_dir << endl;
  cout << "Debug printout: " << ( debug_printout_arg ? "enabled" : "disabled" ) << endl;
  cout << "Action: " << action_str << endl;
  cout << "GA configuration:" << endl;
  cout << "  Population: " << ga_population << endl;
  cout << "  Max generations: " << ga_generation_max << endl;
  cout << "  Best stall max: " << ga_best_stall_max << endl;
  cout << "  Elite count: " << ga_elite_count << endl;
  cout << "  Crossover fraction: " << ga_crossover_fraction << endl;
  cout << "  Mutation rate: " << ga_mutation_rate << endl;
  cout << "  Mutate threshold: " << ga_mutate_threshold << endl;
  cout << "  Crossover threshold: " << ga_crossover_threshold << endl;
  
  string datadir = static_data_dir;
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


#if( defined(SpecUtils_USE_WT_THREADPOOL) )
  // We need to start the WServer so we can populat the thread pool.
  const string docroot = ".";
  const string wt_config = "data/config/wt_config_localweb.xml";
  const int rval = InterSpecServer::start_server( argv[0], "user_data",
                                                   docroot.c_str(),
                                                   wt_config.c_str(),
                                                   static_cast<short int>(0) );

  cout << "ThreadPool will have " << Wt::WServer::instance()->ioService().threadCount() << " threads." << endl;
#endif //#if( defined(SpecUtils_USE_WT_THREADPOOL) )


  /*
   Known issues with "truth" CSV files
   - The truth lines are before random-summing, so the observed peaks will be smaller than the CSV says
   - The annihilation 511 keV will not be included (sources like Na22 will have 511 keV, but "truth" will claim be slightly smaller than actual)
   - Reaction lines not be accounted for, e.g. 488 keV Alpha-Li reaction
   - Single and double escape peaks not present in the "truth" lines.
   - I dont think sum lines are accounted for
   */

  const string base_dir = data_base_dir;

  vector<string> hpges {
    //"Detective-X",
    //"Detective-EX",
    //"Detective-X_noskew",
    //"Falcon 5000",
    //"Fulcrum40h",
    //"LANL_X",
    "HPGe_Planar_50%"
  };

  if( debug_printout_arg )
    hpges = vector<string>{ "Detective-X" };

#if( WRITE_ALL_SPEC_TO_HTML )
  hpges = vector<string>{ "Detective-X", "IdentiFINDER-R500-NaI" };
#endif


  vector<string> live_times{
    //"30_seconds",
    //"300_seconds",
    "1800_seconds"
  };

#if( WRITE_ALL_SPEC_TO_HTML )
  live_times = {"300_seconds"};
#endif

  vector<string> wanted_cities{
    "Livermore"
    //, "Baltimore",
    //"Denver"
  };


  const std::tuple<std::vector<DetectorInjectSet>,std::vector<DataSrcInfo>> &loaded_data
      = PeakFitImproveData::load_inject_data_with_truth_info( base_dir, hpges, live_times, wanted_cities );

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
      /*
      //   Apply best settings found; Generation 491, population best score: -18.1076, population average score: -14.9656
      best_settings.num_smooth_side_channels = 9;
      best_settings.smooth_polynomial_order = 2;
      best_settings.threshold_FOM = 0.758621;
      best_settings.more_scrutiny_FOM_threshold = 1.598265;
      best_settings.pos_sum_threshold_sf = 0.119178;
      best_settings.num_chan_fluctuate = 1;
      best_settings.more_scrutiny_coarser_FOM = best_settings.threshold_FOM + 1.451548;
      best_settings.more_scrutiny_min_dev_from_line = 5.866464;
      best_settings.amp_to_apply_line_test_below = 6;
       */

      //Generation [48], Best=-19.7537, Average=-19.7535, Best genes:
      //Using all live time and HPGe detectors, and cities - 20250912
      best_settings.num_smooth_side_channels = 9;
      best_settings.smooth_polynomial_order = 2;
      best_settings.threshold_FOM = 0.758621;
      best_settings.more_scrutiny_FOM_threshold = 1.598265;
      best_settings.pos_sum_threshold_sf = 0.119178;
      best_settings.num_chan_fluctuate = 1;
      best_settings.more_scrutiny_coarser_FOM = 3.001943;
      best_settings.more_scrutiny_min_dev_from_line = 6.816465;
      best_settings.amp_to_apply_line_test_below = 6.000000;

      std::function<double( const InitialPeakFindSettings &)> ga_eval_fcn 
            = [best_settings, &input_srcs]( const InitialPeakFindSettings &settings ) -> double {
        double score_sum = 0.0;
        for( const DataSrcInfo &info : input_srcs )
        {
          const bool multithread = (PeakFitImprove::sm_num_threads_per_individual > 1);
          const double score = InitialFit_GA::eval_initial_peak_find_and_fit( settings, best_settings, info, false ).find_weight;
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
      candidate_settings.threshold_FOM = 0.758621;
      candidate_settings.more_scrutiny_FOM_threshold = 1.598265;
      candidate_settings.pos_sum_threshold_sf = 0.119178;
      candidate_settings.num_chan_fluctuate = 1;
      candidate_settings.more_scrutiny_coarser_FOM = 3.001943;
      candidate_settings.more_scrutiny_min_dev_from_line = 6.816465;
      candidate_settings.amp_to_apply_line_test_below = 6.000000;
      best_settings = candidate_settings;


      InitialPeakFindSettings initial_fit_settings;
      //Generation [138], Best=-36152.8, Average=-36151.4, - from optimization in Sep. 2025
      initial_fit_settings.initial_stat_threshold = 1.951264;
      initial_fit_settings.initial_hypothesis_threshold = 0.673169;
      initial_fit_settings.initial_min_nsigma_roi = 2.246770;
      initial_fit_settings.initial_max_nsigma_roi = 6.378162;
      initial_fit_settings.fwhm_fcn_form = InitialPeakFindSettings::FwhmFcnForm::SqrtPolynomialTwoCoefs;
      initial_fit_settings.search_roi_nsigma_deficit = 4.241748;
      //initial_fit_settings.search_stat_threshold = 8.051485;
      initial_fit_settings.search_hypothesis_threshold= 3.342207;
      initial_fit_settings.search_stat_significance = 2.025582;
      initial_fit_settings.ROI_add_nsigma_required = 3.526017;
      initial_fit_settings.ROI_add_chi2dof_improve = 0.516229;

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
      cout << "InterSpec 1.0.dev, " << interspec_num_found << ", " << (total_num_expected - interspec_num_found)<< ", " << interspec_num_extra << endl;
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
      //create_n42_peak_fits_for_dir( "/Users/wcjohns/Downloads/spec 2" );
      create_n42_peak_fits( inject_sets, input_srcs );
      break;

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
            CandidatePeak_GA::find_candidate_peaks( info.src_info.src_spectra[0], 0, 0, best_settings );
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

          vector<PeakDef> candidate_peaks = CandidatePeak_GA::find_candidate_peaks( data, 0, 0, best_settings );

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
#if( USE_LM_PEAK_FIT )
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
      //fit_settings.search_stat_threshold = 3;//Reasonable search range of 0 and 8.
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

        auto worker = [&](){
          if( debug_printout_arg )
          {
            std::lock_guard<std::mutex> lock( score_mutex );
            cout << "Evaluating " << info.location_name << "/" << info.live_time_name << "/" << info.detector_name << "/" << info.src_info.src_name << endl;
          }

          const InitialFit_GA::PeakFindAndFitWeights weight = InitialFit_GA::eval_initial_peak_find_and_fit( fit_settings, best_settings, info, false );

          std::lock_guard<std::mutex> lock( score_mutex );
          sum_find_weight += weight.find_weight;
          sum_def_wanted_area_weight += weight.def_wanted_area_weight;
          sum_maybe_wanted_area_weight += weight.maybe_wanted_area_weight;
          sum_not_wanted_area_weight += weight.not_wanted_area_weight;
          sum_def_wanted_area_median_weight += weight.def_wanted_area_median_weight;
          sum_maybe_wanted_area_median_weight += weight.maybe_wanted_area_median_weight;
          sum_not_wanted_area_median_weight += weight.not_wanted_area_median_weight;

          if( debug_printout_arg )
            cout << "For " << info.location_name << "/" << info.live_time_name << "/" << info.detector_name << "/" << info.src_info.src_name
            << ", got weight=" << sum_find_weight << endl;
        };

        if( debug_printout_arg )
        {
          worker();
        }else
        {
          pool.post( std::move(worker) );
        }

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

#if( defined(SpecUtils_USE_WT_THREADPOOL) )
  InterSpecServer::killServer();
#endif

  return EXIT_SUCCESS;
}//int main( int argc, const char * argv[] )


