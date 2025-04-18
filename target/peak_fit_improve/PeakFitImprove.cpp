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

#include <boost/filesystem.hpp>
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
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakFit.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakDists.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/InterSpecServer.h"

#include "openGA.hpp"

using namespace std;

// While optimizing `FindCandidateSettings`, its better to not return a bunch of PeakDefs,
//  so while we're transitioning phases here we will return the vector of tuple results, as
//  well as PeakDefs
#define RETURN_PeakDef_Candidates 1

const bool debug_printout = false;


/** The weights and limits to apply to the various optimizations.
 
 These are all chosen using 'expert judgment', which means they may not be that well chosen.
 */
namespace JudgmentFactors
{
  // Lets not worry about super-small peaks, even where there is little to no background
  const double min_truth_peak_area = 5;
  
  //Photopeak clusters below this next number of sigma will be discarded, since we really
  //  shouldnt find these peaks.
  //  i.e., The threshold at which we will start punishing if a peak is not expected at this
  //  threshold; these source photopeaks have already been removed.
  const double min_truth_nsigma = 1.0;
  
  const double def_want_nsigma = 4;   // i.e., above 4 sigma, lets weight all peaks the same
  const double min_def_wanted_counts = 15; //i.e., if expected peak area is below 15 counts, we wont punish for not finding
  const double lower_want_nsigma = 2; // The number of sigma above which we will positively reward finding a peak
  // Between `def_want_nsigma` and `lower_want_nsigma` we will linearly weight for not finding a peak
  
  const double found_extra_punishment = 0.25; // 1/this-value gives the trade-off of finding extra peaks, verses not finding peaks
  
  // When a peak between lower_want_nsigma and def_want_nsigma is found, the minimum value we should assign
  const double min_initial_fit_maybe_want_score = 0.25;
  
  // Cost for fitting an extra peak after the initial proper fit
  const double initial_fit_extra_peak_punishment = 0.75;
  
  // Multiple for the fraction of additional candidates tried, that didnt stick around
  //  E.g., with a value of 1.0, if 100% of tried peaks failed, this would be equivalent of not
  //    finding one peaks
  const double extra_add_fits_punishment = 1.0;
}//namespace JudgmentFactors



struct FindCandidateSettings
{
  int num_smooth_side_channels = 4; // low res more
  int smooth_polynomial_order = 3;  // highres 3, lowres 2
  double threshold_FOM = 1.3;  // accept peaks higher than this FOM
  double more_scrutiny_FOM_threshold = 3.5; // Peaks bellow this get extra scrutiny
  float pos_sum_threshold_sf = -0.01f;
  
  /** For second-derivative, how many channels are required to be above threshold, in-order to signal a transition */
  size_t num_chan_fluctuate = 2;
  
  //float min_counts_per_channel = 1.0f;
  float more_scrutiny_coarser_FOM = 5.0f;
  
  /** The minimum Chi2 required, of at least one channel in ROI, to be above a straight
   line predicted by the channels on either side of the ROI.
   */
  float more_scrutiny_min_dev_from_line = 4.0;
  
  float amp_to_apply_line_test_below = 40;
  
  std::string print( const string &var_name ) const
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
  
  std::string to_json() const
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
};//struct FindCandidateSettings


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
//      "cont_type_step_improve" (refit peak with different steps, and this is the chi2/dof needed to accept solution)
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
struct InitialPeakFindSettings
{
  /** Filter significance, after initial fit, using `chi2_significance_test(...)` function
   
   These two parameters should be optimized separate from everything else
   
   `initial_stat_threshold`: This is how incompatible with background/continuum the data
   must be, before a peak is allowed to exist. A reasonable search range of values is maybe between 0 and 8.
   
   `initial_hypothesis_threshold`:  this specifies how well the peak must match in shape to a gaussian
   in order to keep the peak.  The higher this number, the more like a gaussian it must be. It is actually the ratio of
   the null hypothesis chi2 to the test hypothesis chi2.  A reasonable search range of values is maybe between -0.1 and 10.
   */
  double initial_stat_threshold, initial_hypothesis_threshold;
  
  double initial_min_nsigma_roi, initial_max_nsigma_roi;
  
  enum class FwhmFcnForm : int {
    Gadras,
    
    SqrtPolynomialTwoCoefs,  //FWHM = sqrt( Sum_i{A_i*pow(x/1000,i)} );
    SqrtPolynomialThreeCoefs,
    
    SqrtEnergyPlusInverse, // FWHM = `sqrt(A0 + A1*E + A2/E)`
    
    NumFwhmFcnForm
  };
  
  FwhmFcnForm fwhm_fcn_form = FwhmFcnForm::SqrtPolynomialTwoCoefs;
  
  
  /** Adds peaks based on now kinda known FWHM functional form.
   Lets just do something stupid, and slide along non-ROI areas, and have an ROI ~5 sigma wide,
   and draw a line using surrounding channels to fit a flat continuum, and if the deficit is
   some reasonable value (2.5 sigma? - need to do a time-tradeoff comparison to get something
   reasonable), then call `fit_amp_and_offset(...)` and `chi2_significance_test(...)`
   
   `search_roi_nsigma_deficit`: Reasonable search range of 1 to 10.
   `search_stat_threshold`: Reasonable search range of 0 and 8.
   `search_hypothesis_threshold`: Reasonable search range of 0 and 8. -0.1 and 10.
   `search_stat_significance`: Reasonable search range of 1 and 6.
   */
  double search_roi_nsigma_deficit, search_stat_threshold, search_hypothesis_threshold, search_stat_significance;
  
  /** Add peaks to ROIs.  WIP.
   `ROI_add_nsigma_required`: required previous and new peaks to be better than. Reasonable search range: 1 to 8
   `ROI_add_chi2dof_improve`: Reasonable search range: 0 to 8
   */
  double ROI_add_nsigma_required, ROI_add_chi2dof_improve;
  
  std::string print( const string &var_name ) const
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
};//struct InitialPeakFindSettings
  


struct FinalPeakFitSettings
{
  /** Combine ROIs, based on how near means are, and how much initial ROIs are overlapping.
   
   `combine_nsigma_near`: Reasonable search range of 1 to 10
   `combine_ROI_overlap_frac`: Reasonable search range of -1 to 1
   */
  double combine_nsigma_near, combine_ROI_overlap_frac;
  
   
  /** How many nsigma the peak needs to be to consider modifying the continuum type.
   
   Reasonable search range: 5 to 25. ???
   */
  double cont_type_peak_nsigma_threshold;
  
  /** How many stat higher the left is than the right, where stat is based of side edge area, to consider a stepped continuum.
   
   Reasonable search range: 1 to 20. ???
   */
  double cont_type_left_right_nsigma;
  
  /** To see if flat, linear, or bilinear stepped continuum
   //Cont P1: -1.6985 on left, -0.3472 on right
   */
  //double cont_type_sum_slopes;
  
  /** Refit peak with different steps, and this is the chi2/dof needed to accept solution
   
   Reasonable search range: 0 to 4. ???
   */
  double cont_type_step_improve;
  
  /** How much of an improvement quadratic or cubic continuum type needs to be to actually use.
   
   Reasonable search range: 0 to 4. ???
   */
  double cont_poly_order_increase_chi2dof_required;
  
  /** Check if we should add skew.
      Right now maybe just just ratio of area between continuum and data for -1 to -4 sigma from mean,
      compared to peak area, or maybe the ratio of sqrt of the areas.
   
   Reasonable search range: 0 to 10. ???
   */
  double skew_nsigma;
  
  /** How much adding skew needs to improve the ROI threshold in order to keep it.
   
   Reasonable search range: 0 to 4. ???
   */
  double skew_improve_chi2_dof_threshold;
  
  
  /** Determine ROI left and right base-widths for single peak ROIs.
   Width will be modified based on stat uncert of initial fit.
   
   Reasonable search range: 0.5 to 7. ???
  */
  double roi_extent_low_num_fwhm_base, roi_extent_high_num_fwhm_base;
  
  /** How to add/subtract based the multiple based of statistical significance.
   */
  enum class RoiExtentMultType : int
  {
    Linear, Sqrt
  };
  
  RoiExtentMultType roi_extent_mult_type;
  
  /** The multiple to add/subtract from ROI width
   
   The multiple can be peak-area within +-1 FWHM, and continuum-area in +- 1 FWHM
   
   Reasonable search range: 0 to 1. ???
   */
  double roi_extent_lower_side_stat_multiple, roi_extent_upper_side_stat_multiple;
  
  //- Add in upper-bounds on stat uncert for peak (e.g., no differenace above 60 sigma or something)
  //- Add in a lower-bound for min number of sigma before allowing step continuum
  //- Add in seperate extent for stepped peaks
  //- Add in a left vs right height, that is an "or" to just trying to fit a step
  //- Or maybe it should be height mult or divided by stat significance
  
  
  /** The multiple to add/subtract from multiple-peak ROI widths
   
   Use single peak peak value for starting, then have an additional add/subtract on each side
   
   Reasonable search range: -1 to 2. ???
   */
  double multi_roi_extent_lower_side_fwhm_mult, multi_roi_extent_upper_side_fwhm_mult;
  
  /** A value to see how many times each ROI should be refit for the final fit.
   
   Reasonable search range: 1 to 3.
   */
  int num_refit_final;
};//struct FinalPeakFitSettings


namespace CandidatePeak_GA
{
  std::function<double(const FindCandidateSettings &)> ns_ga_eval_fcn;
  
  struct CandidatePeakSolution
  {
    int num_smooth_side_channels;
    int smooth_polynomial_order;
    double threshold_FOM;
    double more_scrutiny_FOM_threshold_delta;
    double pos_sum_threshold_sf;
    int num_chan_fluctuate;
    double more_scrutiny_coarser_FOM_delta;
    double more_scrutiny_min_dev_from_line;
    int amp_to_apply_line_test_below;

    string to_string( const string &separator ) const
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
  };

  struct CandidatePeakCost
  {
    // This is where the results of simulation
    // is stored but not yet finalized.
    double objective1;
  };

  typedef EA::Genetic<CandidatePeakSolution,CandidatePeakCost> GA_Type;
  typedef EA::GenerationType<CandidatePeakSolution,CandidatePeakCost> Generation_Type;

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
      
      X_new.num_smooth_side_channels += shrink_scale*(rnd01()-rnd01()); //not multiplying by `mu`, because we can get stuck in a single int
      //in_range=in_range&&(X_new.num_smooth_side_channels>=2 && X_new.num_smooth_side_channels<15);
      X_new.num_smooth_side_channels = std::max( X_new.num_smooth_side_channels, 6 );
      X_new.num_smooth_side_channels = std::min( X_new.num_smooth_side_channels, 13 );
      in_range=in_range&&(X_new.num_smooth_side_channels>=6 && X_new.num_smooth_side_channels<=13);
      
      //X_new.smooth_polynomial_order+=mu*(rnd01()-rnd01()); //This is an int we could get stuck in...
      //in_range=in_range&&(X_new.smooth_polynomial_order>=2 && X_new.smooth_polynomial_order<3);
      
      X_new.threshold_FOM += mu * (rnd01() - rnd01());
      //in_range=in_range&&(X_new.threshold_FOM>=0.75 && X_new.threshold_FOM<3.5);
      in_range=in_range&&(X_new.threshold_FOM>=0.75 && X_new.threshold_FOM<=2);
      
      X_new.more_scrutiny_FOM_threshold_delta += mu*(rnd01() - rnd01());
      //in_range=in_range&&(X_new.more_scrutiny_FOM_threshold_delta>=-0.5 && X_new.more_scrutiny_FOM_threshold_delta<3.5);
      in_range=in_range&&(X_new.more_scrutiny_FOM_threshold_delta>=-0.1 && X_new.more_scrutiny_FOM_threshold_delta<=2);
      
      X_new.pos_sum_threshold_sf += 0.05 * mu * (rnd01() - rnd01()); //mu*(rnd01() - rnd01());
      //in_range=in_range&&(X_new.pos_sum_threshold_sf>=-0.1 && X_new.pos_sum_threshold_sf<0.1);
      in_range=in_range&&(X_new.pos_sum_threshold_sf>=-0.15 && X_new.pos_sum_threshold_sf<=0.15);
      
      //X_new.num_chan_fluctuate+=mu*(rnd01()-rnd01());  //THis is an int we could get stuck in...
      //in_range=in_range&&(X_new.num_chan_fluctuate>=1 && X_new.num_chan_fluctuate<4);
      
      X_new.more_scrutiny_coarser_FOM_delta += mu*(rnd01() - rnd01());
      //in_range=in_range&&(X_new.more_scrutiny_coarser_FOM_delta>=-0.1 && X_new.more_scrutiny_coarser_FOM_delta<5);
      in_range=in_range&&(X_new.more_scrutiny_coarser_FOM_delta>=-0.1 && X_new.more_scrutiny_coarser_FOM_delta<4);
      
      X_new.more_scrutiny_min_dev_from_line += 2 * mu*(rnd01()-rnd01());
      //in_range=in_range&&(X_new.more_scrutiny_min_dev_from_line>=0.0 && X_new.more_scrutiny_min_dev_from_line<10.0);
      in_range=in_range&&(X_new.more_scrutiny_min_dev_from_line>=0.0 && X_new.more_scrutiny_min_dev_from_line<7.0);
      
      X_new.amp_to_apply_line_test_below += 5 * mu*(rnd01()-rnd01());
      X_new.amp_to_apply_line_test_below = std::max( X_new.amp_to_apply_line_test_below, 0 );
      X_new.amp_to_apply_line_test_below = std::min( X_new.amp_to_apply_line_test_below, 100 );
      in_range=in_range&&(X_new.amp_to_apply_line_test_below>=0.0 && X_new.amp_to_apply_line_test_below<=100);
      //in_range=in_range&&(X_new.amp_to_apply_line_test_below>=0.0 && X_new.amp_to_apply_line_test_below<=75);
      
    } while(!in_range);
    return X_new;
  }

  CandidatePeakSolution crossover(
    const CandidatePeakSolution& X1,
    const CandidatePeakSolution& X2,
    const std::function<double(void)> &rnd01)
  {
    CandidatePeakSolution X_new;
    
    double r;
    r=rnd01();
    X_new.num_smooth_side_channels=r*X1.num_smooth_side_channels+(1.0-r)*X2.num_smooth_side_channels;
    //X_new.num_smooth_side_channels = X1.num_smooth_side_channels;
    
    r=rnd01();
    //X_new.smooth_polynomial_order=r*X1.smooth_polynomial_order+(1.0-r)*X2.smooth_polynomial_order;
    X_new.smooth_polynomial_order = X1.smooth_polynomial_order;
    
    r=rnd01();
    X_new.threshold_FOM=r*X1.threshold_FOM+(1.0-r)*X2.threshold_FOM;
    
    r=rnd01();
    X_new.more_scrutiny_FOM_threshold_delta=r*X1.more_scrutiny_FOM_threshold_delta+(1.0-r)*X2.more_scrutiny_FOM_threshold_delta;
    
    r=rnd01();
    X_new.pos_sum_threshold_sf=r*X1.pos_sum_threshold_sf+(1.0-r)*X2.pos_sum_threshold_sf;
    
    r=rnd01();
    //X_new.num_chan_fluctuate=r*X1.num_chan_fluctuate+(1.0-r)*X2.num_chan_fluctuate;
    X_new.num_chan_fluctuate = X1.num_chan_fluctuate;
    
    r=rnd01();
    X_new.more_scrutiny_coarser_FOM_delta=r*X1.more_scrutiny_coarser_FOM_delta+(1.0-r)*X2.more_scrutiny_coarser_FOM_delta;
    
    r=rnd01();
    X_new.more_scrutiny_min_dev_from_line=r*X1.more_scrutiny_min_dev_from_line+(1.0-r)*X2.more_scrutiny_min_dev_from_line;
    
    r=rnd01();
    X_new.amp_to_apply_line_test_below=r*X1.amp_to_apply_line_test_below+(1.0-r)*X2.amp_to_apply_line_test_below;
    
    return X_new;
  }

  double calculate_SO_total_fitness(const GA_Type::thisChromosomeType &X)
  {
    // finalize the cost
    double final_cost=0.0;
    final_cost+=X.middle_costs.objective1;
    return final_cost;
  }

  std::ofstream output_file;
  
  bool m_set_best_genes = false;
  CandidatePeakSolution m_best_genes;
  double m_best_total_cost = 1.0E99;
  
  void SO_report_generation(
    int generation_number,
    const EA::GenerationType<CandidatePeakSolution,CandidatePeakCost> &last_generation,
    const CandidatePeakSolution& best_genes)
  {
    bool best_yet = false;
    if( !m_set_best_genes || (last_generation.best_total_cost < m_best_total_cost) )
    {
      best_yet = true;
      m_set_best_genes = true;
      m_best_genes = best_genes;
      m_best_total_cost = last_generation.best_total_cost;
    }
    
    cout <<"Generation ["<<generation_number<<"], "
    <<"Best="<<last_generation.best_total_cost << ", "
    <<"Average="<<last_generation.average_cost << ", "
    <<"Best genes: {\n\t" <<best_genes.to_string("\n\t")  << "\n}\n"
    <<"Exe_time="<<last_generation.exe_time
    << endl << endl;

    output_file
      <<generation_number<<"\t"
      <<last_generation.average_cost<<"\t"
      <<last_generation.best_total_cost<<"\t"
      << "{" << best_genes.to_string(", ") << "}\n\n";
  }

  FindCandidateSettings do_ga_eval( std::function<double(const FindCandidateSettings &)> ga_eval_fcn )
  {
    assert( !RETURN_PeakDef_Candidates ); //dont want to do a static_assert - since I am still compiling all this code both ways...
    if( RETURN_PeakDef_Candidates )
    {
      cerr << "Please change 'RETURN_PeakDef_Candidates' to false and recompile (for efficiency)." << endl;
      exit(1);
    }
    
    
    assert( !!ga_eval_fcn );
    if( !ga_eval_fcn )
      throw runtime_error( "Invalid eval function passed in." );
    
    ns_ga_eval_fcn = ga_eval_fcn;
    
    
    output_file.open("results.txt");
    output_file<<"step"<<"\t"<<"cost_avg"<<"\t"<<"cost_best"<<"\t"<<"solution_best"<<"\n";

    EA::Chronometer timer;
    timer.tic();

    GA_Type ga_obj;
    ga_obj.problem_mode=EA::GA_MODE::SOGA; //Single objective genetic algorithm
    ga_obj.multi_threading=true;
    ga_obj.idle_delay_us=1; // switch between threads quickly
    ga_obj.dynamic_threading = true; //If false,  thread responsibilities are divided at the beginning
    ga_obj.verbose=false;
    ga_obj.population=1000;
    ga_obj.generation_max=1000;
    ga_obj.calculate_SO_total_fitness=calculate_SO_total_fitness;
    ga_obj.init_genes=init_genes;
    ga_obj.eval_solution=eval_solution;
    ga_obj.mutate=mutate;
    ga_obj.crossover=crossover;
    ga_obj.SO_report_generation=SO_report_generation;
    ga_obj.crossover_fraction=0.7;
    ga_obj.mutation_rate=0.2;
    ga_obj.best_stall_max=10;
    ga_obj.elite_count=10;
    ga_obj.N_threads = 8; //Keep some free cores on my M1 max so I can still use the computer
    EA::StopReason stop_reason = ga_obj.solve();
    
    cout << "Stop reason was: " << ga_obj.stop_reason_to_string( stop_reason) << endl;
    cout << "The problem is optimized in " << timer.toc() << " seconds." << endl;
    
    output_file.close();
    
    return genes_to_settings( m_best_genes );
  }
}//namespace CandidatePeak_GA



namespace InitialFit_GA
{
  atomic<size_t> ns_num_evals_this_generation( 0 );
  
  std::function<double( const InitialPeakFindSettings &)> ns_ga_eval_fcn;
  
  
  
  // I think we could actually just use `InitialPeakFindSettings` instead of defining this struck - but I'll wait on that until we get it up and going a bit
  struct InitialFitSolution
  {
    double initial_stat_threshold;
    double initial_hypothesis_threshold;
    double initial_min_nsigma_roi;
    double initial_max_nsigma_roi;
    int fwhm_fcn_form;
    double search_roi_nsigma_deficit;
    double search_stat_threshold;
    double search_hypothesis_threshold;
    double search_stat_significance;
    double ROI_add_nsigma_required;
    double ROI_add_chi2dof_improve;


    string to_string( const string &separator ) const
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
    
  };//struct InitialFitSolution

  static InitialPeakFindSettings genes_to_settings( const InitialFitSolution &solution )
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
  
  
  struct InitialFitCost
  {
    // This is where the results of simulation
    // is stored but not yet finalized.
    double objective1;
  };

  typedef EA::Genetic<InitialFitSolution,InitialFitCost> GA_Type;
  typedef EA::GenerationType<InitialFitSolution,InitialFitCost> Generation_Type;

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

  InitialFitSolution mutate(
    const InitialFitSolution& X_base,
    const std::function<double(void)> &rnd01,
    double shrink_scale)
  {
    InitialFitSolution X_new;
    const double mu = 0.2*shrink_scale; // mutation radius (adjustable)
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
      X_new.initial_stat_threshold += mu*(rnd01()-rnd01());
      in_range = in_range&&(X_new.initial_stat_threshold>=0.0 && X_new.initial_stat_threshold<8);
      X_new.initial_hypothesis_threshold+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.initial_hypothesis_threshold>=0.25 && X_new.initial_hypothesis_threshold<3.0);
      X_new.initial_min_nsigma_roi+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.initial_min_nsigma_roi>=2 && X_new.initial_min_nsigma_roi<8);
      X_new.initial_max_nsigma_roi+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.initial_max_nsigma_roi>=4 && X_new.initial_max_nsigma_roi<12);
      
      X_new.fwhm_fcn_form = X_base.fwhm_fcn_form; //unnecessary, but to be explicit
      //X_new.fwhm_fcn_form+=mu*(rnd01()-rnd01());
      //in_range=in_range&&(X_new.fwhm_fcn_form>=0 && X_new.fwhm_fcn_form<3);
      
      X_new.search_roi_nsigma_deficit+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.search_roi_nsigma_deficit>=3 && X_new.search_roi_nsigma_deficit < 10.0);
      X_new.search_stat_threshold+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.search_stat_threshold>=1 && X_new.search_stat_threshold<9);
      X_new.search_hypothesis_threshold+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.search_hypothesis_threshold>=-0.1 && X_new.search_hypothesis_threshold<10.0);
      X_new.search_stat_significance+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.search_stat_significance>=1 && X_new.search_stat_significance<6);
      X_new.ROI_add_nsigma_required+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.ROI_add_nsigma_required>=1 && X_new.ROI_add_nsigma_required<8);
      X_new.ROI_add_chi2dof_improve+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.ROI_add_chi2dof_improve>=0.0 && X_new.ROI_add_chi2dof_improve<8);
    } while(!in_range);
    return X_new;
  }

  InitialFitSolution crossover(
    const InitialFitSolution& X1,
    const InitialFitSolution& X2,
    const std::function<double(void)> &rnd01)
  {
    InitialFitSolution X_new;
    double r;
    r=rnd01();
    X_new.initial_stat_threshold=r*X1.initial_stat_threshold+(1.0-r)*X2.initial_stat_threshold;
    r=rnd01();
    X_new.initial_hypothesis_threshold=r*X1.initial_hypothesis_threshold+(1.0-r)*X2.initial_hypothesis_threshold;
    r=rnd01();
    X_new.initial_min_nsigma_roi=r*X1.initial_min_nsigma_roi+(1.0-r)*X2.initial_min_nsigma_roi;
    r=rnd01();
    X_new.initial_max_nsigma_roi=r*X1.initial_max_nsigma_roi+(1.0-r)*X2.initial_max_nsigma_roi;
    r=rnd01();
    X_new.fwhm_fcn_form=r*X1.fwhm_fcn_form+(1.0-r)*X2.fwhm_fcn_form;
    r=rnd01();
    X_new.search_roi_nsigma_deficit=r*X1.search_roi_nsigma_deficit+(1.0-r)*X2.search_roi_nsigma_deficit;
    r=rnd01();
    X_new.search_stat_threshold=r*X1.search_stat_threshold+(1.0-r)*X2.search_stat_threshold;
    r=rnd01();
    X_new.search_hypothesis_threshold=r*X1.search_hypothesis_threshold+(1.0-r)*X2.search_hypothesis_threshold;
    r=rnd01();
    X_new.search_stat_significance=r*X1.search_stat_significance+(1.0-r)*X2.search_stat_significance;
    r=rnd01();
    X_new.ROI_add_nsigma_required=r*X1.ROI_add_nsigma_required+(1.0-r)*X2.ROI_add_nsigma_required;
    r=rnd01();
    X_new.ROI_add_chi2dof_improve=r*X1.ROI_add_chi2dof_improve+(1.0-r)*X2.ROI_add_chi2dof_improve;
    return X_new;
  }

  double calculate_SO_total_fitness(const GA_Type::thisChromosomeType &X)
  {
    // finalize the cost
    double final_cost=0.0;
    final_cost+=X.middle_costs.objective1;
    return final_cost;
  }

  std::ofstream output_file;

  
  bool m_set_best_genes = false;
  InitialFitSolution m_best_genes;
  double m_best_total_cost = 1.0E99;
  
  
  void SO_report_generation(
    int generation_number,
    const EA::GenerationType<InitialFitSolution,InitialFitCost> &last_generation,
    const InitialFitSolution& best_genes)
  {
    bool best_yet = false;
    if( !m_set_best_genes || (last_generation.best_total_cost < m_best_total_cost) )
    {
      best_yet = true;
      m_set_best_genes = true;
      m_best_genes = best_genes;
      m_best_total_cost = last_generation.best_total_cost;
    }
    
    output_file
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
  
    output_file.open( "initial_fit_results.txt" );
    output_file << "step" << "\t" << "cost_avg" << "\t" << "cost_best" << "\t" << "solution_best" << "\n";

    EA::Chronometer timer;
    timer.tic();

    GA_Type ga_obj;
    ga_obj.problem_mode=EA::GA_MODE::SOGA;
    ga_obj.multi_threading = true;
    ga_obj.idle_delay_us = 1; // switch between threads quickly
    ga_obj.dynamic_threading = true;
    ga_obj.verbose = false;
    ga_obj.population = 100;
    ga_obj.generation_max = 450;
    ga_obj.calculate_SO_total_fitness=calculate_SO_total_fitness;
    ga_obj.init_genes=init_genes;
    ga_obj.eval_solution=eval_solution;
    ga_obj.mutate=mutate;
    ga_obj.crossover=crossover;
    ga_obj.SO_report_generation=SO_report_generation;
    ga_obj.crossover_fraction=0.7;
    ga_obj.mutation_rate=0.2;
    ga_obj.best_stall_max=10;
    ga_obj.elite_count=10;
    ga_obj.N_threads = 8; //Keep some free cores on my M1 max so I can still use the computer
    EA::StopReason stop_reason = ga_obj.solve();

    cout<<"The problem is optimized in "<<timer.toc()<<" seconds."<<endl;

    output_file.close();
    
    cout << "Stop reason was: " << ga_obj.stop_reason_to_string( stop_reason) << endl;
    cout << "The problem is optimized in " << timer.toc() << " seconds." << endl;
    
    output_file.close();
    
    return genes_to_settings( m_best_genes );
  }
}//namespace InitialFit_GA


namespace G2k
{
  struct G2kPeak
  {
    char multiplet;
    int PeakNum;
    float Energy, PeakSig, FWHM, ContinuumCounts, NetPeakArea, NetAreaError, NetCountRate;
  };
  
  vector<G2kPeak> g2k_peaks_from_file( ifstream &strm )
  {
    bool found_peaks_start = false;
    string line;
    while( !found_peaks_start && SpecUtils::safe_get_line(strm, line) )
    {
      SpecUtils::trim( line );
      if( !found_peaks_start
         && SpecUtils::icontains( line, "*                              Peak Search Report                               *") )
      {
        while( !found_peaks_start && SpecUtils::safe_get_line(strm, line) )
        {
          SpecUtils::trim( line );
          if( SpecUtils::istarts_with( line, "-----" ) )
          {
            found_peaks_start = true;
            break;
          }
        }
      }
    }//while( !found_peaks_start && SpecUtils::safe_get_line(strm, line) )
    
    if( !found_peaks_start )
      throw runtime_error( "Failed to find peaks part of report" );
    
    vector<G2kPeak> result;
    while( SpecUtils::safe_get_line(strm, line) )
    {
      SpecUtils::trim( line );
      if( line.empty() )
        continue;
      
      if( SpecUtils::istarts_with(line, "Errors" ) )
        break;
      
      vector<string> fields;
      SpecUtils::split( fields, line, " \t");
      
      assert( (fields.size() == 8) || (fields.size() == 9) );
      if( (fields.size() != 8) && (fields.size() != 9) )
        throw runtime_error( "Got line w/o 8 or 9 fields: '" + line + "'" );
      
      assert( (fields.size() == 8) || (fields[0] == "M") || (fields[0] == "m") );
      if( (fields.size() == 9) && ((fields[0] != "M") && (fields[0] != "m")) )
        throw runtime_error( "First field is not M or m: '" + line + "'" );
      
      const char mutltiplet = (fields.size() == 9) ? fields[0].at(0) : ' ';
      if( fields.size() == 9 )
      {
        line = line.substr( 1 );
        SpecUtils::trim( line );
      }
      
      vector<float> values;
      SpecUtils::split_to_floats( line.c_str(), values, " \t", false );
      assert( values.size() == 8 );
      
      G2kPeak peak;
      peak.multiplet = mutltiplet;
      peak.PeakNum = static_cast<int>( values[0] );
      peak.Energy = values[1];
      peak.PeakSig = values[2];
      peak.FWHM = values[3];
      peak.ContinuumCounts = values[4];
      peak.NetPeakArea = values[5];
      peak.NetAreaError = 0.5*values[6]; //the RPT files give uncert at 2-sigma
      peak.NetCountRate = values[7];
      
      assert( result.empty() || ((result.back().PeakNum + 1) == peak.PeakNum) );
      result.push_back( peak );
    }//while( SpecUtils::safe_get_line(strm, line) )
    
    return result;
  };//auto g2k_peaks_from_file
}//namespace G2k


std::vector<PeakDef> find_candidate_peaks( const std::shared_ptr<const SpecUtils::Measurement> data,
                          size_t start_channel,
                          size_t end_channel,
                          std::vector< std::tuple<float,float,float> > &results,
                          const FindCandidateSettings &settings )
{
  // Plan:
  //  - Should we split the spectrum up into N regions?  So we can refine have fewer smoothing bins lower in energy?
  //    - For HPGe it looks like we can avoid this (but out of time to investigate)
  std::vector<PeakDef> result_peaks;
  
  results.clear();
  
  const size_t nchannel = data->num_gamma_channels();
  
  //We will let one bin fluctate negative to avoid fluctations near threshold_FOM
  //Untested for HPGe data.
  //  Currently (20141209) for low res data, this should be kept the same as
  //  in find_roi_for_2nd_deriv_candidate(...) {although all this is under
  //  development}.
  const size_t nFluxuate = settings.num_chan_fluctuate;

  if( nchannel < (2*nFluxuate + 2) || (nchannel < 4) )
    return;

  if( start_channel >= nchannel )
    start_channel = 0;
  if( end_channel <= start_channel || end_channel >= (nchannel-1) )
    end_channel = nchannel - nFluxuate - 1;
  
  
  const double threshold_FOM = settings.threshold_FOM;
  const float pos_sum_threshold_sf = settings.pos_sum_threshold_sf;
  
    
  const int order = settings.smooth_polynomial_order; //highres ? 3 : 2;
  const size_t side_bins = settings.num_smooth_side_channels; //highres ? 4 : std::max( size_t(5), static_cast<size_t>( 0.022f*midenergy/midbinwidth + 0.5f ) );
  
  const vector<float> &spectrum = *data->gamma_counts();
  
  vector<float> second_deriv, smoothed_data, coarser_data, rougher_second_deriv;
  
  // Create a smoothed second order derivative; its negative-most values will correspond
  //  to peak centers
  smoothSpectrum( spectrum, static_cast<int>(side_bins), order, 2, second_deriv );
  
  // We will smooth the data a bit, so we can sum it to help determine stat significance
  smoothSpectrum( spectrum, static_cast<int>(side_bins), 1, 0, smoothed_data );
  
  // We will use `rougher_second_deriv` to evaluate sigma, and ROI extents, because if
  // `side_bins`, than the predicted sigma and ROI extents will be a lot wider than
  // expected.
  smoothSpectrum( spectrum, 3, 3, 2, rougher_second_deriv );
  
  
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
      
      //cout << "mean=" << mean << ", pos2ndDerivAreaSum=" << pos2ndDerivAreaSum
      //<< ", positive2ndDataSum=" << positive2ndDataSum << ", secondsum=" << secondsum
      //<< " - seconsumratios=" << -secondsum/pos2ndDerivAreaSum << endl;
      
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
      
      if( rougher_second_deriv[rough_index] < 0.0 )
      {
        //size_t nFluxuate = 1; //
        size_t lower_pos_index = rough_index, upper_pos_index = rough_index;
        while( lower_pos_index > nFluxuate )
        {
          bool all_pos = true;
          for( size_t ind = lower_pos_index; all_pos && (ind > (lower_pos_index - nFluxuate)); --ind )
            all_pos = (rougher_second_deriv[ind] >= -FLT_EPSILON);
          if( all_pos )
            break;
          
          lower_pos_index -= 1;
        }
        
        
        while( upper_pos_index < (end_channel - nFluxuate) )
        {
          bool all_pos = true;
          for( size_t ind = upper_pos_index; all_pos && (ind < (upper_pos_index + nFluxuate)); ++ind )
            all_pos = (rougher_second_deriv[ind] >= -FLT_EPSILON);
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
        const double rough_est_sigma = sqrt( std::max(rough_data_area,1.0) );
        
        
        rougher_FOM = 0.68*rougher_amplitude / rough_est_sigma;
        
        //cout << "mean=" << mean << " -> sigma=" << sigma
        //<< ", rougher_sigma=" << rougher_sigma
        //<< ", rougher_amp=" << rougher_amplitude
        //<< ", rougher_FOM=" << rougher_FOM << endl;
        
        
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
        }
        //assert( tester == (lower_pos_index - roi_begin_index) );
        //lower_data_counts /= (lower_pos_index - roi_begin_index);
          
        //tester = 0;
        for( ; ((roi_end_index+1) < end_channel) && (rougher_second_deriv[roi_end_index] > 0.0); ++roi_end_index )
        {
          //++tester;
          //upper_data_counts += spectrum[roi_end_index];
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
      }else
      {
        //cout << "mean=" << mean << " had a positive rougher 2nd deriv" << endl;
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
      

      double est_sigma = sqrt( std::max(data_area,1.0) );
      
      //In principle we would want to use the (true) continuums area to derive
      //  the est_sigma from, but for practical purposes I think this can give
      //  us false positives fairly often
      
      const double figure_of_merit = 0.68*amplitude / est_sigma;
      
      bool passed_higher_scrutiny = true;
      if( figure_of_merit < settings.more_scrutiny_FOM_threshold
         && (settings.more_scrutiny_FOM_threshold >= threshold_FOM) )
      {
        passed_higher_scrutiny = ((rougher_second_deriv[rough_index] < 0.0)
                                  && (rougher_FOM >= settings.more_scrutiny_coarser_FOM)
                                  );
        
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
        
        if( debug_printout )
          cout << "Straight-line Chi2 = " << chi2_dof << " and max_chi2=" << max_chi2 << " for energy=" << mean << endl;
        
        //vector<double> cont_equation{0.0, 0.0};
        //PeakContinuum::eqn_from_offsets( roi_begin_index, roi_end_index, mean, data,
        //                                cont_est_channels, cont_est_channels,
        //                                cont_equation[1], cont_equation[0] );
        
        passed_higher_scrutiny = (max_chi2 > settings.more_scrutiny_min_dev_from_line);
      }//if( check Chi2 of region )
      
      
      
      if( (figure_of_merit < 5) && (rougher_second_deriv[rough_index] >= 0.0) )
        passed_higher_scrutiny = false;
      
      //const double min_required_data = settings.min_counts_per_channel*2*deriv_sigma;
      
      if( (figure_of_merit > threshold_FOM) && passed_higher_scrutiny )
      {
        if( debug_printout )
          cout << "Accepted: energy=" << mean << ", FOM=" << figure_of_merit
          << ", amp=" << amplitude << ", FWHM=" << sigma*2.35482f
          << ", data_area=" << data_area << ", rougher_FOM=" << rougher_FOM
          //<< ", min_required=" << min_required_data
          << endl;
        
        results.push_back( std::tuple<float,float,float>{mean, sigma, amplitude} );
        
#if( RETURN_PeakDef_Candidates )
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
#endif
      }else
      {
       if( debug_printout )
        cout << "Rejected: energy=" << mean << ", FOM=" << figure_of_merit
        << ", amp=" << amplitude << ", FWHM=" << sigma*2.35482f
        << ", data_area=" << data_area  << ", rougher_FOM=" << rougher_FOM
        //<< ", min_required=" << min_required_data
        << endl;
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






struct PeakTruthInfo
{
  float energy;
  float cps;
  float area;
  float fwhm;
  
  /**
   Where:
   const double SC = 2.0 + 0.02*(lowSkew + highSkew)*std::pow(energy/661.0, resolutionPower );
   FullWidth = SC * GetFWHM(energy, resolutionOffset, resolution661, resolutionPower );
   */
  float full_width;
  
  /**
   "S.E.", "D.E.", "EscapeXRay", "Peak"
   */
  std::string label;
  
  PeakTruthInfo()
  : energy( 0.0 ), cps( 0.0 ), area( 0.0 ), fwhm( 0.0 ), full_width( 0.0 ), label()
  {
  }
  
  explicit PeakTruthInfo( const std::string &line )
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
           || (label == "Peak+S.E.") );
  }//PeakTruthInfo(...)
};//struct PeakTruthInfo


struct InjectSourceInfo
{
  string file_base_path;
  string src_name;
  
  vector<PeakTruthInfo> source_lines;
  vector<PeakTruthInfo> background_lines;
  
  /** The PCF file for this source. */
  shared_ptr<SpecUtils::SpecFile> spec_file;
  
  /** The source + background spectra, that are Poisson varied. */
  vector<shared_ptr<const SpecUtils::Measurement>> src_spectra;
  
  /** Background that is Poisson varied, and same duration as src spectra. */
  shared_ptr<const SpecUtils::Measurement> short_background;
  
  /** A longer background, with Poisson variation. */
  shared_ptr<const SpecUtils::Measurement> long_background;
  
  /** Source + background without Poisson variation. */
  shared_ptr<const SpecUtils::Measurement> src_no_poisson;
  
  /** Background, with no Poisson variation. */
  shared_ptr<const SpecUtils::Measurement> background_no_poisson;
};//struct InjectSourceInfo


struct DetectorInjectSet
{
  string detector_name;
  string location_name;
  string live_time_name;
  
  vector<InjectSourceInfo> source_infos;
};//struct DetectorInjectSet



struct ExpectedPhotopeakInfo
{
  /** An approximate detection sigma, i.e., 2.33 times this value is kinds 95% DL */
  double nsigma_over_background;
  
  double roi_lower;
  double roi_upper;
  double effective_energy;
  double peak_area;
  double continuum_area;
  
  vector<PeakTruthInfo> gamma_lines;
};//struct ExpectedPhotopeakInfo


struct DataSrcInfo
{
  string detector_name;
  string location_name;
  string live_time_name;
  
  InjectSourceInfo src_info;
  
  vector<ExpectedPhotopeakInfo> expected_photopeaks;
};//struct DataSrcInfo

ExpectedPhotopeakInfo create_expected_photopeak( const InjectSourceInfo &info, const vector<PeakTruthInfo> &lines )
{
  assert( !lines.empty() );
  
  ExpectedPhotopeakInfo peak;
  peak.gamma_lines = lines;
  
  peak.peak_area = 0.0;
  
  // We'll weight the mean and ROI lower/upper by each lines relative area
  //  - not really sure if this is the right thing to do...
  peak.roi_lower = 0.0;
  peak.roi_upper = 0.0;
  peak.effective_energy = 0.0;
  
  for( const auto &l : lines )
  {
    peak.peak_area += l.area;
    peak.roi_lower += l.area * (l.energy - 0.66*l.full_width); // Using 0.66 instead of 0.5 because skew, or whatever, I assume, so none of the dataset ever has more predicted than actually in the data
    peak.roi_upper += l.area * (l.energy + 0.66*l.full_width);
    peak.effective_energy += l.area * l.energy;
  }//for( const auto &l : lines )
  
  peak.roi_lower /= peak.peak_area;
  peak.roi_upper /= peak.peak_area;
  peak.effective_energy /= peak.peak_area;
  
  // Use the non-Poisson varied spectrum to calculate total area
  const double total_area = info.src_no_poisson->gamma_integral( static_cast<float>(peak.roi_lower),
                                                              static_cast<float>(peak.roi_upper) );

  if( !( ((total_area - peak.peak_area) > -8.0*sqrt(peak.peak_area))
         || ((total_area - peak.peak_area) > -5.0)
         || (peak.roi_upper > info.src_no_poisson->gamma_energy_max()) ) )
  {
    // This looks to trigger for very large are peak (like 500k counts), where there is a good
    //   amount of peak-skew.
    cerr << "Warning: create_expected_photopeak: For peak at " << peak.effective_energy << " keV,"
    << " got peak.peak_area=" << peak.peak_area << ", while histogram total_area=" << total_area
    << " between " << peak.roi_lower << " and " << peak.roi_upper << " keV, for \n\t"
    << info.file_base_path << "." << endl;
  }
  
  assert( (total_area > 1.0E5) //large skew with such areas
         || ((total_area - peak.peak_area) > -8.0*sqrt(peak.peak_area))
         || ((total_area - peak.peak_area) > -5.0)
         || (peak.roi_upper > info.src_no_poisson->gamma_energy_max()) 
         || (info.src_no_poisson->live_time() < 0.85*info.src_no_poisson->real_time()) );
  
  peak.continuum_area = std::max( 0.0, total_area - peak.peak_area );
  peak.nsigma_over_background = peak.peak_area / sqrt( std::max(1.0, peak.continuum_area) );
  
  
  //if( peak.effective_energy > 1407.5 && peak.effective_energy < 1408.5 )
  //{
  //  cout << "This one" << endl;
  //}
  
  assert( peak.roi_lower < peak.roi_upper );
  assert( peak.effective_energy < peak.roi_upper );
  assert( peak.effective_energy > peak.roi_lower );
  
  return peak;
}//ExpectedPhotopeakInfo create_expected_photopeak(...)




InjectSourceInfo parse_inject_source_files( const string &base_name )
{
  /*
  CSV file is ordered as:
   DRF  Detector/CZT/1.5cm-2cm-2cm
   Source  Br76_Sh
   Distance  50
   MeasurementTime  1800
   BackgroundLocation  Livermore
   InjectDesc  Source + Background
              
             
   Source Lines:
   Energy (keV)  CPS  Area  FWHM  FullWidth  Label
   100.3  8.39348e-09  1.51083e-05  4.51382  12.7386  D.E.
   178.7  1.89767e-07  0.00034158  5.49774  15.7744  EscapeXRay
   182.23  4.80212e-07  0.000864382  5.53458  15.8892  EscapeXRay
   209.7  0.000884833  1.5927  5.80637  16.7385  Peak
   
   ...
   
   Background Lines:
   Energy (keV)  CPS  Area  FWHM  FullWidth  Label
   0.46  1.65619e-08  2.98114e-05  2.98932  7.44028  EscapeXRay
   0.6  8.78397e-11  1.58111e-07  2.98932  7.47825  EscapeXRay
   3  9.11796e-11  1.64123e-07  2.98932  7.73021  EscapeXRay
   
   */
  
  const string pcf_filename = base_name + ".pcf";
  const string csv_filename = base_name + "_gamma_lines.csv";
  
  auto spec_file = make_shared<SpecUtils::SpecFile>();
  const bool loaded_pcf = spec_file->load_pcf_file( pcf_filename );
  assert( loaded_pcf );
  if( !loaded_pcf )
    throw runtime_error( "failed to load " + pcf_filename );
  
  const vector<shared_ptr<const SpecUtils::Measurement>> meass = spec_file->measurements();
  assert( meass.size() == 16 );
  if( meass.size() != 16 )
    throw runtime_error( "Unexpected number of measurements in " + pcf_filename );
  
  InjectSourceInfo info;
  info.file_base_path = base_name;
  info.src_name = SpecUtils::filename( base_name );
  info.spec_file = spec_file;
  
  if( meass[0]->title() != info.src_name )
    cerr << "Mismatch of src: '" << meass[0]->title() << "', vs '" << info.src_name << "'" << endl;
  
  assert( meass[0]->title() == info.src_name );
  assert( meass[1]->title() == "Background" );
  assert( meass[2]->title() == "Background" );
  assert( meass[3]->title() == info.src_name );
  assert( meass[4]->title() == "Background (non-Poisson)" );
  
  info.src_spectra.push_back( meass[0] );
  info.short_background = meass[1];
  info.long_background = meass[2];
  info.src_no_poisson = meass[3];
  info.background_no_poisson = meass[4];
  
  for( size_t i = 5; i < meass.size(); ++i )
  {
    assert( meass[i]->title() == info.src_name );
    info.src_spectra.push_back( meass[i] );
  }
  
  
  ifstream csv_strm( csv_filename.c_str(), (ios::binary | ios::in) );
  assert( csv_strm.good() );
  if( !csv_strm.good() )
    throw runtime_error( "Failed to open CSV: " + csv_filename );
  
  bool found_src_lines = false, found_background_lines = false;
  string line;
  while( !found_src_lines && SpecUtils::safe_get_line(csv_strm, line) )
    found_src_lines = SpecUtils::istarts_with(line, "Source Lines:");
  
  assert( found_src_lines );
  if( !found_src_lines )
    throw runtime_error( "Failed to find Source Lines in " + csv_filename );
  
  SpecUtils::safe_get_line(csv_strm, line);
  assert( line == "Energy (keV),CPS,Area,FWHM,FullWidth,Label" );
  if( line != "Energy (keV),CPS,Area,FWHM,FullWidth,Label" )
    throw runtime_error( "Unexpected foreground CSV header: " + line );
  
  while( SpecUtils::safe_get_line(csv_strm, line) )
  {
    found_background_lines = SpecUtils::istarts_with(line, "Background Lines:");
    if( found_background_lines )
      break;
    
    SpecUtils::trim( line );
    if( line.empty() )
      continue;
    
    info.source_lines.emplace_back( line );
  }//while( loop over foreground lines )
  
  assert( found_background_lines );
  SpecUtils::safe_get_line(csv_strm, line);
  assert( line == "Energy (keV),CPS,Area,FWHM,FullWidth,Label" );
  if( line != "Energy (keV),CPS,Area,FWHM,FullWidth,Label" )
    throw runtime_error( "Unexpected background CSV header: " + line );
  
  
  while( SpecUtils::safe_get_line(csv_strm, line) )
  {
    SpecUtils::trim( line );
    if( line.empty() )
      continue;
    
    info.background_lines.emplace_back( line );
  }//while( loop over background lines )
  
  return info;
}//InjectSourceInfo parse_inject_source_files( const string &base_name )


/**
 
 */
vector<PeakDef> initial_peak_find_and_fit( const InitialPeakFindSettings &fit_settings,
                              const FindCandidateSettings &candidate_settings,
                              const shared_ptr<const SpecUtils::Measurement> &data,
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
  const vector<PeakDef> candidate_peaks = find_candidate_peaks( data, 0, 0, dummy, candidate_settings );
  
  
#warning "initial_peak_find_and_fit: always assuming HPGe right now - for dev"
  //bool isHPGe = (data->num_gamma_channels() > 5000); // This will be wrong a lot...
  const bool isHPGe = true;
  // TODO: improve selection of HPGe here
  
  bool amplitudeOnly = false;
  vector<PeakDef> zeroth_fit_results, initial_fit_results;
  const std::vector<PeakDef> dummy_fixedpeaks;
  
  if( false )
  {//Begin optimization block - delte when done optimizing `fitPeaks(...)`
    const double start_wall = SpecUtils::get_wall_time();
    const double start_cpu = SpecUtils::get_cpu_time();
    
    for( size_t i = 0; i < 20; ++i )
    //while( true )
    {
      vector<PeakDef> zeroth_fit_results, initial_fit_results;
      fitPeaks( candidate_peaks,
               fit_settings.initial_stat_threshold, fit_settings.initial_hypothesis_threshold,
               data, zeroth_fit_results, dummy_fixedpeaks, amplitudeOnly, isHPGe );
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
  
  auto fit_peaks_per_roi = [&fit_settings, &data, isHPGe]( const vector<PeakDef> &input_peaks ) -> vector<PeakDef> {
    
    const bool amplitudeOnly = false;
    const std::vector<PeakDef> dummy_fixedpeaks;
    
    map<const PeakContinuum *,vector<PeakDef>> roi_to_peaks_map;
    for( const PeakDef &p : input_peaks )
      roi_to_peaks_map[p.continuum().get()].push_back( p );
    
    vector<vector<PeakDef>> fit_rois( roi_to_peaks_map.size() );
    
    SpecUtilsAsync::ThreadPool threadpool;
    
    size_t fit_rois_index = 0;
    for( auto &roi_peaks : roi_to_peaks_map )
    {
      const vector<PeakDef> &peaks = roi_peaks.second;
      
      vector<PeakDef> &results = fit_rois[fit_rois_index];

      //20250409: check per

      threadpool.post( boost::bind( &fitPeaks,
                                   boost::cref(peaks),
                                   fit_settings.initial_stat_threshold,
                                   fit_settings.initial_hypothesis_threshold,
                                   data,
                                   boost::ref( results ),
                                   dummy_fixedpeaks,
                                   amplitudeOnly,
                                   isHPGe ) );
      
      fit_rois_index += 1;
    }//for( size_t peakn = 0; peakn < seperated_peaks.size(); ++peakn )
    
    threadpool.join();
    
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
    if( debug_printout )
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
      
      if( debug_printout )
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
      }//if( debug_printout )
      
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



double eval_initial_peak_find_and_fit( const InitialPeakFindSettings &fit_settings,
                                      const FindCandidateSettings &candidate_settings,
                                      const DataSrcInfo &src_info )
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
                        = initial_peak_find_and_fit( fit_settings, candidate_settings, data,
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
  
  if( debug_printout )
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
  
  // blah - blah blah - so now we need to test this function by walking through a few spectra.
  if( debug_printout )
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
  
  
  return score;
}//double eval_initial_peak_find_and_fit(...)


vector<PeakDef> final_peak_fit( const vector<PeakDef> &pre_fit_peaks,
                               const FinalPeakFitSettings &final_fit_settings,
                               const DataSrcInfo &src_info )
{
  // blah blah blah - fit peaks
  assert( 0 );
  throw runtime_error( "final_peak_fit not implemented yet" );
  
  return {};
}//vector<PeakDef> final_peak_fit


double eval_final_peak_fit( const FindCandidateSettings &candidate_settings,
                           const InitialPeakFindSettings &initial_fit_settings,
                           const FinalPeakFitSettings &final_fit_settings,
                           const DataSrcInfo &src_info )
{
  // blah blah blah - fit peaks
  assert( 0 );
  throw runtime_error( "eval_final_peak_fit not implemented yet" );
  
  return 0.0;
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
        *result = initial_peak_find_and_fit( initial_fit_settings, candidate_settings, data, dummy1, dummy2 );
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

/*
  std::function<double( const InitialPeakFindSettings &)> ga_eval_fcn
        = [best_settings, &input_srcs]( const InitialPeakFindSettings &settings ) -> double {
    double score_sum = 0.0;
    for( const DataSrcInfo &info : input_srcs )
    {
      const double score = eval_initial_peak_find_and_fit( settings, best_settings, info );
      score_sum += score;
    }

    return -score_sum;
  };// set InitialFit_GA::ns_ga_eval_fcn

  const InitialPeakFindSettings best_initial_fit_settings = InitialFit_GA::do_ga_eval( ga_eval_fcn );
*/






  //TODO:
  // - Should time using LinearProblemSubSolveChi2Fcn - may be faster/better?
  // - Also, try using Ceres instead of Minuit.  Look for `USE_LM_PEAK_FIT`
  // -Compare using PeakFitLM::fit_peak_for_user_click_LM(...) vs
  //  `refitPeaksThatShareROI(...)`
  // - See what works better for peaks - `LinearProblemSubSolveChi2Fcn` or `PeakFitChi2Fcn` - e.g. more accurate answers for some default fit scenarios.
  
}//void do_final_peak_fit_ga_optimization(...)



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
  
  vector<DetectorInjectSet> inject_sets;
  
  for( boost::filesystem::directory_iterator detector_itr(base_dir);
      detector_itr != boost::filesystem::directory_iterator(); ++detector_itr )
  {
    //detector_itr->path().filename()
    
    if( !boost::filesystem::is_directory(detector_itr->status()) )
      continue;
    
    const boost::filesystem::path detector_path = detector_itr->path();
    
    vector<string> hpges {
      "Detective-X",
      "Detective-EX",
      "Detective-X_noskew",
      "Falcon 5000",
      "Fulcrum40h",
      "LANL_X",
      "HPGe_Planar_50%"
    };
    
    if( debug_printout )
      hpges = vector<string>{ "Detective-X" };
    
    bool is_wanted_det = false;
    for( const string &det : hpges )
      is_wanted_det |= (detector_path.filename() == det);
    
    if( !is_wanted_det )
    {
      cerr << "Skipping detector " << detector_path.filename() << endl;
      continue;
    }
    
    cout << "In detector directory: " << detector_path.filename() << endl;
    
    for( boost::filesystem::directory_iterator city_itr(detector_path);
        city_itr != boost::filesystem::directory_iterator(); ++city_itr )
    {
      if( !boost::filesystem::is_directory(city_itr->status()) )
        continue;
      
      const boost::filesystem::path city_path = city_itr->path();
      cout << "In city directory: " << city_path.filename() << endl;
      
      if( city_path.filename() != "Livermore" )
      {
        cout << "Skipping city " << city_itr->path().filename() << endl;
        continue;
      }
      
      // Now loop over {"30_seconds", "300_seconds", "1800_seconds"}
      for( boost::filesystem::directory_iterator livetime_itr(city_path);
          livetime_itr != boost::filesystem::directory_iterator(); ++livetime_itr )
      {
        if( !boost::filesystem::is_directory(livetime_itr->status()) )
          continue;
        
        const boost::filesystem::path livetime_path = livetime_itr->path();
        
        
        if( debug_printout )
        {
          vector<string> live_times{
            "30_seconds",
            "300_seconds",
            "1800_seconds"
          };
          
          bool is_wanted_lt = false;
          for( const string &lt : live_times )
            is_wanted_lt |= (livetime_path.filename() == lt);
          
          if( !is_wanted_lt )
          {
            cerr << "Skipping live-time " << livetime_path.filename() << endl;
            continue;
          }
        }//if( debug_printout )
        
        //cout << "In livetime: " << livetime_path << endl;
        
        // Now loop over PCF files
        vector<string> files_to_load_basenames;
        for( boost::filesystem::directory_iterator file_itr(livetime_path);
            file_itr != boost::filesystem::directory_iterator(); ++file_itr )
        {
          if( !boost::filesystem::is_regular_file(file_itr->status()) )
            continue;
          
          const string pcf_name = file_itr->path().string();
          if( !SpecUtils::iends_with(pcf_name, ".pcf") )
            continue;
          
          const string base_name = pcf_name.substr(0, pcf_name.size() - 4 );
          
          //base_name + "_raw_gamma_lines.txt"
          const string csv_name = base_name + "_gamma_lines.csv";
          
          if( !boost::filesystem::is_regular_file(csv_name) )
          {
            assert( 0 );
            cerr << "No gamma lines CSV for " << file_itr->path() << " - skipping" << endl;
            continue;
          }
          
          if( debug_printout )
          {
            const vector<string> wanted_sources{
              //"Ac225_Unsh",
              "Eu152_Sh",
              //"Fe59_Phant",
              //"Am241_Unsh"
            };
            bool wanted = wanted_sources.empty();
            for( const string &src : wanted_sources )
            {
              const bool want = (base_name.find(src) != string::npos);
              wanted |= want;
              if( want )
                cout << "Setting info for '" << src << "'" << endl;
            }
            if( !wanted )
              continue;
          }//if( debug_printout )
          
          files_to_load_basenames.push_back( base_name );
        }//for( loop over files in directory )
        
        const size_t num_srcs = files_to_load_basenames.size();
        
        inject_sets.push_back( DetectorInjectSet{} );
        
        DetectorInjectSet &injects = inject_sets.back();
        injects.detector_name = detector_path.filename().string();
        injects.location_name = city_path.filename().string();
        injects.live_time_name = livetime_path.filename().string();
        injects.source_infos.resize( num_srcs );
        
        SpecUtilsAsync::ThreadPool pool;
        
        for( size_t file_index = 0; file_index < num_srcs; ++file_index )
        {
          InjectSourceInfo *info = &(injects.source_infos[file_index]);
          const string base_name = files_to_load_basenames[file_index];
          
          pool.post( [info,base_name](){
            *info = parse_inject_source_files( base_name );
          } );
        }//for( loop over files to parse )
        
        pool.join();
        
        // A smoke check to make sure nothing got messed up
        assert( injects.source_infos.size() == num_srcs );
        for( size_t index = 0; index < num_srcs; ++index )
        {
          assert( !injects.source_infos[index].source_lines.empty() );
          assert( !injects.source_infos[index].background_lines.empty() );
          assert( injects.source_infos[index].file_base_path == files_to_load_basenames[index] );
          assert( !injects.source_infos[index].src_spectra.empty() );
          assert( injects.source_infos[index].src_spectra[0]->title() == SpecUtils::filename(files_to_load_basenames[index]) );
        }//for( size_t index = 0; index < injects.source_infos.size(); ++index )
        
        cout << "Parsed " << injects.source_infos.size() << " sources" << endl;
      }//for( loop over live-time directories, livetime_itr )
    }//for( loop over cities, city_itr )
  }//for( loop over detector types )
  
  
  
  
  vector<DataSrcInfo> input_srcs;
  size_t num_inputs = 0, num_accepted_inputs = 0;
  for( const DetectorInjectSet &inject_set : inject_sets )
  {
    for( const InjectSourceInfo &info : inject_set.source_infos )
    {
      // For the moment, we'll look for visible lines on the first Poisson varied source
      //cout << "For " << inject_set.detector_name << "/"
      //<< inject_set.live_time_name << "/" << inject_set.location_name
      //<< " source " << info.src_name << ": ";
      //cout << "For '" << info.file_base_path << "':" << endl;
      
      num_inputs += 1;
      
      // The "truth" lines are all before random-summing, so the observed peaks will be smaller
      //  in area than the truth CSV says, but if we keep dead-time low, this wont be noticeable.
      //  Right now have chosen 2%, fairly arbitrarily
      //  For EX-100
      //   1%: Lose 14 out of 223 files
      //   2%: Lose  8 out of 223 files
      //   5%: Lose  4 out of 223 files
      //  10%: Lose  4 out of 223 files
      
//#if( !RETURN_PeakDef_Candidates )
#warning "Only using <2% dead-time files"
      if( info.src_no_poisson->live_time() < 0.98*info.src_no_poisson->real_time() )
      {
        continue;
      }
//#else
//#warning "Using all Live-Times for peak search"
//#endif
      
      const shared_ptr<const SpecUtils::Measurement> meas = info.src_spectra[0];
      vector<PeakTruthInfo> lines = info.source_lines;
      lines.insert( end(lines), begin(info.background_lines), end(info.background_lines) );
      
      //for( const PeakTruthInfo &line : lines )
      //  cout << "Source line: " << line.energy << " keV, area=" << line.area << endl;
      
      /**
       - Sort lines from largest area to smallest.
       - cluster, using `~0.75*FWHM`
       - See what is hit, and whats not
       */
      
      std::sort( begin(lines), end(lines), []( const PeakTruthInfo &lhs, const PeakTruthInfo &rhs ){
        return (lhs.area < rhs.area);
      } );
      
      
      const double cluster_fwhm_multiple = 0.5;
      vector<vector<PeakTruthInfo>> clustered_lines;
      while( !lines.empty() )
      {
        const PeakTruthInfo main_line = std::move( lines.back() );
        lines.resize( lines.size() - 1 ); //remove the line we just grabbed
        
        vector<PeakTruthInfo> cluster;
        cluster.push_back( main_line );
        
        // Look all through `lines` for lines within `cluster_fwhm_multiple` of main_line
        deque<size_t> index_to_remove;
        for( size_t i = 0; i < lines.size(); ++i )
        {
          if( fabs(lines[i].energy - main_line.energy) < cluster_fwhm_multiple*main_line.fwhm )
          {
            index_to_remove.push_back( i );
            cluster.push_back( lines[i] );
          }
        }//for( loop over lines to cluster )
        
        for( auto iter = std::rbegin(index_to_remove); iter != std::rend(index_to_remove); ++iter )
        {
//          cout << "Removing " << *iter << ", which has energy " << lines[*iter].energy
//          << ", main line=" << main_line.energy << endl;
          lines.erase( begin(lines) + (*iter) );
        }
        
        clustered_lines.push_back( cluster );
      }//while( !lines.empty() )
      
      
      
      vector<ExpectedPhotopeakInfo> detectable_clusters;
      for( const vector<PeakTruthInfo> &cluster : clustered_lines )
      {
        
        assert( !cluster.empty() );
        if( cluster.empty() )
          continue;
        
        const PeakTruthInfo &main_line = cluster.front();
        if( (main_line.area < 5.0)  // Let be realistic, and require something
           || (fabs(main_line.energy - 478.0) < 1.2)  // Avoid Alpha-Li reaction
           || (fabs(main_line.energy - 511.0) < 1.0 ) // Avoid D.B. 511
              // Shouldn't there be some other broadened reaction lines here???
           || ((main_line.energy + 0.5*main_line.full_width) > info.src_no_poisson->gamma_energy_max()) // Make sure not off upper-end
           || ((main_line.energy - 0.5*main_line.full_width) < info.src_no_poisson->gamma_energy_min()) // Make sure not below lower-end
           || ((main_line.energy - 0.5*main_line.full_width) < 50.0) //eh, kinda arbitrary, maybe shouldnt limit?
           )
        {
          continue;
        }
        
        ExpectedPhotopeakInfo roi_info = create_expected_photopeak( info, cluster );
        
        if( (roi_info.peak_area > JudgmentFactors::min_truth_peak_area)
           && (roi_info.nsigma_over_background > JudgmentFactors::min_truth_nsigma) )
          detectable_clusters.push_back( std::move(roi_info) );
      }//for( const vector<PeakTruthInfo> &cluster : clustered_lines )
      
      std::sort( begin(detectable_clusters), end(detectable_clusters),
                []( const ExpectedPhotopeakInfo &lhs, const ExpectedPhotopeakInfo &rhs ) -> bool {
        return lhs.effective_energy < rhs.effective_energy;
      } );
      
      if( debug_printout )
      {
        for( const ExpectedPhotopeakInfo &roi_info : detectable_clusters )
        {
          cout << "Expected ROI: {energy: " << roi_info.effective_energy
          << ", fwhm: " << roi_info.gamma_lines.front().fwhm
          << ", PeakArea: " << roi_info.peak_area
          << ", ContinuumArea: " << roi_info.continuum_area
          << ", NSigma: " << roi_info.nsigma_over_background
          << ", ROI: [" << roi_info.roi_lower << ", " << roi_info.roi_upper << "]"
          << "}"
          << endl;
          
          //cout << "cluster: {";
          //for( const auto &c : cluster )
          //  cout << "{" << c.energy << "," << c.area << "}, ";
          //cout << "}" << endl;
        }//for( const ExpectedPhotopeakInfo &roi_info : detectable_clusters )
      }//if( debug_printout )
      
                
      if( !detectable_clusters.empty() )
      {
        num_accepted_inputs += 1;
        DataSrcInfo src_info;
        src_info.detector_name = inject_set.detector_name;
        src_info.location_name = inject_set.location_name;
        src_info.live_time_name = inject_set.live_time_name;
        
        src_info.src_info = info;
        src_info.expected_photopeaks = detectable_clusters;
        
        input_srcs.push_back( src_info );
      }//if( !detectable_clusters.empty() )
    }//for( const InjectSourceInfo &src : source_infos )
  }//for( const DetectorInjectSet &inject_set : inject_sets )
  
  cout << "Used " << num_accepted_inputs << " of total " << num_inputs << " input files." << endl;
  
   
  auto eval_candidate_settings = [&input_srcs]( const FindCandidateSettings settings,
                                                const bool write_n42 )
      -> tuple<double,size_t,size_t,size_t,size_t,size_t> //<score, num_peaks_found, num_possibly_accepted_peaks_not_found, num_extra_peaks>
  {
    double sum_score = 0.0;
    size_t num_possibly_accepted_peaks_not_found = 0, num_def_wanted_not_found = 0;
    size_t num_extra_peaks = 0, num_peaks_found = 0;
    size_t num_def_wanted_peaks_found = 0;
    
    for( const DataSrcInfo &info : input_srcs )
    {
      vector<tuple<float,float,float>> detected_peaks; //{mean, sigma, amplitude}
      shared_ptr<const SpecUtils::Measurement> src_spectrum = info.src_info.src_spectra.front();
      find_candidate_peaks( src_spectrum, 0, 0, detected_peaks, settings );
      
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
            if( debug_printout )
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
  auto eval_candidate_settings_fcn = [&]( const FindCandidateSettings settings, const bool write_n42 ){
    
    const tuple<double,size_t,size_t,size_t,size_t,size_t> result = eval_candidate_settings( settings, write_n42 );
  
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
  //eval_candidate_settings_fcn( best_settings, false );
  
  
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
      const auto ga_eval = [&eval_candidate_settings](const FindCandidateSettings &settings) -> double {
        const tuple<double,size_t,size_t,size_t,size_t,size_t> score = eval_candidate_settings( settings, false );
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
                      
                      pool.post( [eval_candidate_settings_fcn,settings,&score_mutex](){
                        try
                        {
                          eval_candidate_settings_fcn( settings, false );
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
      
      eval_candidate_settings_fcn( best_settings, true );
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
          const double score = eval_initial_peak_find_and_fit( settings, best_settings, info );
          score_sum += score;
        }
        
        return -score_sum;
      };// set InitialFit_GA::ns_ga_eval_fcn
      
      const InitialPeakFindSettings best_initial_fit_settings = InitialFit_GA::do_ga_eval( ga_eval_fcn );
      
      
      //eval_candidate_settings_fcn( best_settings, best_initial_fit_settings, true );
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
      
      
      do_final_peak_fit_ga_optimization( candidate_settings, initial_fit_settings, input_srcs );
      
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
        const vector<G2k::G2kPeak> g2k_peaks_initial = G2k::g2k_peaks_from_file( g2k_strm );
        vector<PeakDef> g2k_peaks;
        for( const G2k::G2kPeak &p : g2k_peaks_initial )
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
        eval_candidate_settings( candidate_settings, true );
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
            find_candidate_peaks( info.src_info.src_spectra[0], 0, 0, results, best_settings );
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
          vector<PeakDef> candidate_peaks = find_candidate_peaks( data, 0, 0, dummy, best_settings );
          
          for( PeakDef &p : candidate_peaks )
          {
            p.continuum()->setType( PeakContinuum::OffsetType::FlatStep );
            //p.continuum()->setType( PeakContinuum::OffsetType::Linear );
          }

          const bool isHPGe = true, amplitudeOnly = false;
          vector<PeakDef> zeroth_fit_results, initial_fit_results;
          const std::vector<PeakDef> dummy_fixedpeaks;
          
          
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
              fitPeaks( candidate_peaks, 0.0, 0.0, data, peaks, dummy_fixedpeaks, amplitudeOnly, isHPGe );
              
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
      
      double sum_weight = 0.0;
      const double start_wall = SpecUtils::get_wall_time();
      const double start_cpu = SpecUtils::get_cpu_time();
      for( const DataSrcInfo &info : input_srcs )
      {
        if( debug_printout )
          cout << "Evaluating " << info.location_name << "/" << info.live_time_name << "/" << info.detector_name << "/" << info.src_info.src_name << endl;
        
        const double weight = eval_initial_peak_find_and_fit( fit_settings, best_settings, info );
        sum_weight += weight;
        if( debug_printout )
          cout << "Got weight=" << weight << endl;
      }//for( const DataSrcInfo &info : input_srcs )
      const double end_wall = SpecUtils::get_wall_time();
      const double end_cpu = SpecUtils::get_cpu_time();

      cout << "Average weight of " << input_srcs.size() << " files is " << (sum_weight/input_srcs.size()) << endl;
      cout << "\t\teval took " << (end_wall - start_wall) << " s, wall and " << (end_cpu - start_cpu) << " s cpu" << endl;
      break;
    }//case OptimizationAction::CodeDev:
  }//switch( action )
  
  
  
  
  return EXIT_SUCCESS;
}//int main( int argc, const char * argv[] )


