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

bool debug_printout = false;


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

    string to_string() const
    {
      return
        string("{")
        +  "num_smooth_side_channels:"+std::to_string(num_smooth_side_channels)
        +", smooth_polynomial_order:"+std::to_string(smooth_polynomial_order)
        +", threshold_FOM:"+std::to_string(threshold_FOM)
        +", more_scrutiny_FOM_threshold_delta:"+std::to_string(more_scrutiny_FOM_threshold_delta)
        +", pos_sum_threshold_sf:"+std::to_string(pos_sum_threshold_sf)
        +", num_chan_fluctuate:"+std::to_string(num_chan_fluctuate)
        +", more_scrutiny_coarser_FOM_delta:"+std::to_string(more_scrutiny_coarser_FOM_delta)
        +", more_scrutiny_min_dev_from_line:"+std::to_string(more_scrutiny_min_dev_from_line)
        +", amp_to_apply_line_test_below:"+std::to_string(amp_to_apply_line_test_below)
        +"}";
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
    p.num_smooth_side_channels=2+13*rnd01();
    p.smooth_polynomial_order=2+1*rnd01();
    p.threshold_FOM=0.8+2.7*rnd01();
    p.more_scrutiny_FOM_threshold_delta=-0.5+4*rnd01();
    p.pos_sum_threshold_sf=-0.1+0.2*rnd01();
    p.num_chan_fluctuate=1+3*rnd01();
    p.more_scrutiny_coarser_FOM_delta=-0.1+5.1*rnd01();
    p.more_scrutiny_min_dev_from_line=0.0+10*rnd01();
    p.amp_to_apply_line_test_below=0.0+120*rnd01();
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
    do{
      in_range=true;
      X_new=X_base;
      X_new.num_smooth_side_channels+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.num_smooth_side_channels>=2 && X_new.num_smooth_side_channels<15);
      X_new.smooth_polynomial_order+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.smooth_polynomial_order>=2 && X_new.smooth_polynomial_order<3);
      X_new.threshold_FOM+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.threshold_FOM>=0.8 && X_new.threshold_FOM<3.5);
      X_new.more_scrutiny_FOM_threshold_delta+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.more_scrutiny_FOM_threshold_delta>=-0.5 && X_new.more_scrutiny_FOM_threshold_delta<3.5);
      X_new.pos_sum_threshold_sf+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.pos_sum_threshold_sf>=-0.1 && X_new.pos_sum_threshold_sf<0.1);
      X_new.num_chan_fluctuate+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.num_chan_fluctuate>=1 && X_new.num_chan_fluctuate<4);
      X_new.more_scrutiny_coarser_FOM_delta+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.more_scrutiny_coarser_FOM_delta>=-0.1 && X_new.more_scrutiny_coarser_FOM_delta<5);
      X_new.more_scrutiny_min_dev_from_line+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.more_scrutiny_min_dev_from_line>=0.0 && X_new.more_scrutiny_min_dev_from_line<10.0);
      X_new.amp_to_apply_line_test_below+=mu*(rnd01()-rnd01());
      in_range=in_range&&(X_new.amp_to_apply_line_test_below>=0.0 && X_new.amp_to_apply_line_test_below<120);
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
    r=rnd01();
    X_new.smooth_polynomial_order=r*X1.smooth_polynomial_order+(1.0-r)*X2.smooth_polynomial_order;
    r=rnd01();
    X_new.threshold_FOM=r*X1.threshold_FOM+(1.0-r)*X2.threshold_FOM;
    r=rnd01();
    X_new.more_scrutiny_FOM_threshold_delta=r*X1.more_scrutiny_FOM_threshold_delta+(1.0-r)*X2.more_scrutiny_FOM_threshold_delta;
    r=rnd01();
    X_new.pos_sum_threshold_sf=r*X1.pos_sum_threshold_sf+(1.0-r)*X2.pos_sum_threshold_sf;
    r=rnd01();
    X_new.num_chan_fluctuate=r*X1.num_chan_fluctuate+(1.0-r)*X2.num_chan_fluctuate;
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
  
  void SO_report_generation(
    int generation_number,
    const EA::GenerationType<CandidatePeakSolution,CandidatePeakCost> &last_generation,
    const CandidatePeakSolution& best_genes)
  {
    m_set_best_genes = true;
    m_best_genes = best_genes;
    
    cout
      <<"Generation ["<<generation_number<<"], "
      <<"Best="<<last_generation.best_total_cost<<", "
      <<"Average="<<last_generation.average_cost<<", "
      <<"Best genes=("<<best_genes.to_string()<<")"<<", "
      <<"Exe_time="<<last_generation.exe_time
      <<endl;

    output_file
      <<generation_number<<"\t"
      <<last_generation.average_cost<<"\t"
      <<last_generation.best_total_cost<<"\t"
      <<best_genes.to_string()<<"\n";
  }

  void do_ga_eval( std::function<double(const FindCandidateSettings &)> ga_eval_fcn )
  {
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
    ga_obj.dynamic_threading=false; //If false,  thread responsibilities are divided at the beginning
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
    EA::StopReason stop_reason = ga_obj.solve();
    
    cout << "Stop reason was: " << ga_obj.stop_reason_to_string( stop_reason) << endl;
    cout << "The problem is optimized in " << timer.toc() << " seconds." << endl;
    
    output_file.close();
    
    return genes_to_settings( m_best_genes );
  }
}//namespace CandidatePeak_GA














// While optimizing `FindCandidateSettings`, its better to not return a bunch of PeakDefs,
//  so while we're transitioning phases here we will return the vector of tuple results, as
//  well as PeakDefs
#define RETURN_PeakDef_Candidates 0

#if( RETURN_PeakDef_Candidates )
std::vector<PeakDef>
#else
void
#endif
find_candidate_peaks( const std::shared_ptr<const SpecUtils::Measurement> data,
                          size_t start_channel,
                          size_t end_channel,
                          std::vector< std::tuple<float,float,float> > &results,
                          const FindCandidateSettings &settings )
{
  // Plan:
  //  - Should we split the spectrum up into N regions?  So we can refine have fewer smoothing bins lower in energy?
  //    - For HPGe it looks like we can avoid this (but out of time to investigate)
  
#if( RETURN_PeakDef_Candidates )
  std::vector<PeakDef> result_peaks;
#endif
  
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
        
        secondzero = ((trial_second_zero > (minbin + 1)) && second_deriv[trial_second_zero-1] < 0.0) 
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
      
      assert( second_deriv[firstzero+1] <= 0.0 );
      assert( second_deriv[secondzero-1] <= 0.0 );
      
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
        
        assert( rougher_second_deriv[lower_pos_index] >= 0 );
        assert( rougher_second_deriv[lower_pos_index+1] <= 0 );
        assert( rougher_second_deriv[upper_pos_index] >= 0 );
        assert( rougher_second_deriv[upper_pos_index-1] <= 0 );
        
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
        
        assert( rougher_second_deriv[lower_pos_index] >= 0 );
        assert( rougher_second_deriv[upper_pos_index] >= 0 );
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
        
        assert( second_deriv[firstzero] >= 0.0 );
        assert( second_deriv[secondzero] >= 0.0 );
        
        // Make sure second_deriv[roi_begin_index] and second_deriv[roi_end_index] are positive,
        //  although they should already be
        for( ; (roi_begin_index > 0) && (second_deriv[roi_begin_index] < 0.0); --roi_begin_index )
        {
        }
        
        for( ; ((roi_end_index+1) < end_channel) && (second_deriv[roi_end_index] < 0.0); ++roi_end_index )
        {
        }
        
        assert( second_deriv[firstzero] >= 0.0 );
        assert( second_deriv[secondzero] >= 0.0 );
        
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
                                      : (end_channel - cont_est_channels);
        
        
        size_t lower_nchannel = 0, upper_nchannel = 0;
        double lower_cnts_per_chnl = 0.0, upper_cnts_per_chnl = 0.0;
        for( size_t i = lower_edge_start; i < roi_begin_index; ++i, ++lower_nchannel  )
          lower_cnts_per_chnl += spectrum[i];
        for( size_t i = roi_end_index + 1; i <= upper_edge_end; ++i, ++upper_nchannel )
          upper_cnts_per_chnl += spectrum[i];
        
        assert( lower_nchannel == (roi_begin_index - lower_edge_start) );
        assert( upper_nchannel == (upper_edge_end - roi_end_index) );
        
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
  
#if( RETURN_PeakDef_Candidates )
  return result_peaks;
#endif
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
    
  assert( ((total_area - peak.peak_area) > -0.000*peak.peak_area)
         || ((total_area - peak.peak_area) > -5.0)
         || (peak.roi_upper > info.src_no_poisson->gamma_energy_max()) );
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
  
  const string base_dir = "/Users/wcjohns/rad_ana/peak_area_optimization/peak_fit_accuracy_inject/";
  
  vector<DetectorInjectSet> inject_sets;
  
  for( boost::filesystem::directory_iterator detector_itr(base_dir);
      detector_itr != boost::filesystem::directory_iterator(); ++detector_itr )
  {
    //detector_itr->path().filename()
    
    if( !boost::filesystem::is_directory(detector_itr->status()) )
      continue;
    
    const boost::filesystem::path detector_path = detector_itr->path();
    
    const vector<string> hpges {
      "Detective-X",
      "Detective-EX",
      "Detective-X_noskew",
      "Falcon 5000",
      "Fulcrum40h",
      "LANL_X",
      "HPGe_Planar_50%"
    };
    
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
        
        
        const vector<string> live_times {
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
          
          /*
          debug_printout = true;
          const vector<string> wanted_sources{
            "Ac225_Unsh",
            //"Eu152_Sh",
            //"Fe59_Phant"
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
        */
          
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
  
  
  
  // Lets not worry about super-small peaks, even where there is little to no background
  const double min_peak_area = 5;
  
  //Photopeak clusters below this next number of sigma will be discarded, since we really
  //  shouldnt find these peaks
  const double not_expected_thresh_nsigma = 1.0;
  
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
      if( info.src_no_poisson->live_time() < 0.98*info.src_no_poisson->real_time() )
      {
        continue;
      }
      
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
        
        if( (roi_info.peak_area > min_peak_area)
           && (roi_info.nsigma_over_background > not_expected_thresh_nsigma) )
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
  
  const double def_want_nsigma = 4;   // i.e., above 4 sigma, lets weight all peaks the same
  const double min_def_wanted_counts = 10; //i.e., if expected peak area is below 10 counts, we wont punish for not finding
  const double lower_want_nsigma = 2; // The number of sigma above which we will positively reward finding a peak
  // Between `def_want_nsigma` and `lower_want_nsigma` we will linearly weight for not finding a peak
  //not_expected_thresh_nsigma - the threshold at which we will start punishing if a peak is not
  //                             expected at this threshold; these source photopeaks have already
  //                             been removed.
  const double found_extra_punishment = 0.25; // 1/this-value gives the trade-off of finding extra peaks, verses not finding peaks
  
  auto eval_settings = [&input_srcs, def_want_nsigma, min_def_wanted_counts, lower_want_nsigma, found_extra_punishment](
                                                            const FindCandidateSettings settings,
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
          
          if( (expected.nsigma_over_background > def_want_nsigma) && (expected.peak_area > min_def_wanted_counts) )
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
           if( expected.nsigma_over_background < lower_want_nsigma )
           score -= 0; //No punishment
           else if( expected.nsigma_over_background < def_want_nsigma )
           score -= ( (def_want_nsigma - expected.nsigma_over_background) / (expected.nsigma_over_background - lower_want_nsigma) );
           else
           score -= 1;
           */
        }else
        {
          // we got a match
          num_detected_expected += 1;
          expected_and_was_detected.push_back( expected );
          
          if( (expected.nsigma_over_background > def_want_nsigma) && (expected.peak_area > min_def_wanted_counts) )
            num_def_wanted_detected += 1;
          
          if( expected.nsigma_over_background > def_want_nsigma )
            score += 1.0;
          else if( expected.nsigma_over_background > lower_want_nsigma )
            score += ((def_want_nsigma - expected.nsigma_over_background) / (expected.nsigma_over_background - lower_want_nsigma));
          else
            score += 0.0; //No punishment, but no reward.
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
      
      score -= found_extra_punishment * num_detected_not_expected;
      
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
        
        string outdir = "/Users/wcjohns/rad_ana/InterSpec/target/peak_fit_improve/build_xcode/output_n42";
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
        
        const string out_n42 = SpecUtils::append_path( outdir, info.src_info.src_name ) + ".n42";
        
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
  };//eval_settings lambda
  
  
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
  auto eval_settings_fcn = [&]( const FindCandidateSettings settings, const bool write_n42 ){
    
    const tuple<double,size_t,size_t,size_t,size_t,size_t> result = eval_settings( settings, write_n42 );
  
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
  //eval_settings_fcn( best_settings, false );
  
  if( debug_printout )
  {
    cerr << "Setting best_settings instead of actually finding them!" << endl;
    best_settings.num_smooth_side_channels = 9;
    best_settings.smooth_polynomial_order = 2;
    best_settings.threshold_FOM = 1.300000;
    best_settings.more_scrutiny_FOM_threshold = 0.000000;
    best_settings.pos_sum_threshold_sf = 0.000000;
    best_settings.num_chan_fluctuate = 2;
    best_settings.more_scrutiny_coarser_FOM = 1.300000;
    best_settings.more_scrutiny_min_dev_from_line = 0.500000;
    best_settings.amp_to_apply_line_test_below = 5.000000;
  }else
  {
    //best_settings = CandidatePeak_GA::do_ga_eval( std::function<double(const FindCandidateSettings &)> ns_ga_eval_fcn )
    
    
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
                    
                    pool.post( [eval_settings_fcn,settings,&score_mutex](){
                      try
                      {
                        eval_settings_fcn( settings, false );
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
     
  }//if( debug_printout ) / else

  eval_settings_fcn( best_settings, true );
  cout << "Wrote N42s with best settings." << endl;
  
  const double end_wall = SpecUtils::get_wall_time();
  const double end_cpu = SpecUtils::get_cpu_time();
  auto now = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());
  
  stringstream outmsg;
  outmsg << "Finish Time: " << SpecUtils::to_iso_string(now) << endl
  << endl
  << best_settings.print("\tbest_settings") << endl 
  << endl
  << "found_extra_punishment = " << found_extra_punishment << endl
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
    ofstream best_settings_file( "best_settings.txt" );
    best_settings_file << outmsg.str();
  }
  
  cout << outmsg.str()<< endl;
  
  return EXIT_SUCCESS;
}//int main( int argc, const char * argv[] )


