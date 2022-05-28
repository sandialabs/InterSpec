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
#include <deque>
#include <limits>
#include <functional>

#include "ceres/ceres.h"

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakFit.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;

namespace
{

struct DoWorkOnDestruct
{
  std::function<void()> m_worker;
  DoWorkOnDestruct( std::function<void()> &&worker ) : m_worker( std::move(worker) ){}
  ~DoWorkOnDestruct(){ if(m_worker) m_worker(); }
};//struct DoWorkOnDestruct


struct RoiRangeChannels : public RelActCalcAuto::RoiRange
{
  // TODO: currently if the energy calibration is adjusted - we have to keep the total number of channels (i.e. the number of residuals constant) we have to move the first/last channel the same number of channels
  size_t first_channel, last_channel, num_channels;
  
  RoiRangeChannels( std::shared_ptr<const SpecUtils::Measurement> foreground, const RelActCalcAuto::RoiRange &roi_range )
  : RelActCalcAuto::RoiRange( roi_range ),
    first_channel( foreground->find_gamma_channel( roi_range.lower_energy ) ),
    last_channel( foreground->find_gamma_channel( roi_range.upper_energy ) ),
    num_channels( 1 + last_channel - first_channel )
  {
    if( roi_range.upper_energy <= roi_range.lower_energy )
      throw runtime_error( "RoiRangeChannels: " );
  }
};//struct RoiRangeChannels


struct NucInputGamma : public RelActCalcAuto::NucInputInfo
{
  vector<SandiaDecay::EnergyRatePair> nominal_gammas;
  
  static size_t remove_gamma( const double energy, vector<SandiaDecay::EnergyRatePair> &gammas )
  {
    const size_t ninitial = gammas.size();
    
    gammas.erase( std::remove_if(begin(gammas), end(gammas),
      [energy](const SandiaDecay::EnergyRatePair &e ){
        return fabs(e.energy - energy) < 0.001;
    }), end(gammas) );
    
    return gammas.size() - ninitial;
  }
  
  NucInputGamma( const RelActCalcAuto::NucInputInfo &info )
   : RelActCalcAuto::NucInputInfo( info )
  {
    if( !nuclide )
      throw runtime_error( "NucInputGamma: null Nuclide." );
    
    if( age < 0.0 )
      throw runtime_error( "NucInputGamma: age may not be negative (" + nuclide->symbol + ": "
                           + PhysicalUnits::printToBestTimeUnits(age)  + ")" );
    
    double nominal_age = info.age;
    
    if( !fit_age )
      nominal_age = PeakDef::defaultDecayTime( nuclide, nullptr );
    
    SandiaDecay::NuclideMixture mix;
    mix.addAgedNuclideByActivity( nuclide, SandiaDecay::Bq, nominal_age );
    nominal_gammas = mix.gammas( 0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy, true );
    
    for( double energy : gammas_to_exclude )
    {
      const size_t did_remove = remove_gamma( energy, nominal_gammas );
      assert( did_remove );
    }
  }//NucInputGamma constructor
};//struct NucInputGamma



struct RelActAutoCostFcn /* : ROOT::Minuit2::FCNBase() */
{
  RelActCalcAuto::Options m_options;
  std::vector<NucInputGamma> m_nuclides;
  std::vector<RoiRangeChannels> m_energy_ranges;
  std::vector<RelActCalcAuto::FloatingPeak> m_extra_peaks;
  
  std::shared_ptr<const SpecUtils::Measurement> m_foreground;
  std::shared_ptr<const SpecUtils::Measurement> m_background;
  std::shared_ptr<SpecUtils::Measurement> m_spectrum;
  
  float m_live_time;
  vector<float> m_channel_counts; //background subtracted channel counts
  vector<float> m_channel_count_uncerts; //e.g. sqrt( m_channel_counts[i] )
  std::shared_ptr<const SpecUtils::EnergyCalibration> m_energy_cal;
  
  /** We will punish */
  double m_rel_eff_punishment;
  
  /** Will either be null, or have FWHM info. */
  std::shared_ptr<const DetectorPeakResponse> m_drf;
  
  
  RelActAutoCostFcn( RelActCalcAuto::Options options,
                     vector<RelActCalcAuto::RoiRange> energy_ranges,
                     vector<RelActCalcAuto::NucInputInfo> nuclides,
                     vector<RelActCalcAuto::FloatingPeak> extra_peaks,
                     shared_ptr<const SpecUtils::Measurement> foreground,
                     shared_ptr<const SpecUtils::Measurement> background,
                     std::shared_ptr<const DetectorPeakResponse> drf
                           )
  : m_options( options ),
  m_nuclides{},
  m_energy_ranges{},
  m_extra_peaks( extra_peaks ),
  m_foreground( foreground ),
  m_background( background ),
  m_live_time( foreground ? foreground->live_time() : 0.0f ),
  m_channel_counts{},
  m_channel_count_uncerts{},
  m_energy_cal{},
  m_rel_eff_punishment( 1000.0 ),
  m_drf( nullptr )
  {
    if( !foreground
       || (foreground->num_gamma_channels() < 128)
       || !foreground->energy_calibration()
       || !foreground->energy_calibration()->valid() )
      throw runtime_error( "RelActAutoCostFcn: invalid foreground spectrum." );
    
    if( background && (background->num_gamma_channels() != foreground->num_gamma_channels()) )
      throw runtime_error( "RelActAutoCostFcn: Diff num background/foreground channels." );
    
    m_energy_cal = foreground->energy_calibration();
    
    assert( foreground->gamma_counts() );
    m_channel_counts = *foreground->gamma_counts();
    m_channel_count_uncerts.resize( m_channel_counts.size(), 0.0 );
    
    if( !background )
    {
      for( size_t i = 0; i < m_channel_counts.size(); ++i )
      {
        const double counts = m_channel_counts[i];
        m_channel_count_uncerts[i] = (counts <= 1.0) ? 1.0 : sqrt( counts );
      }
      
      m_spectrum = make_shared<SpecUtils::Measurement>( *foreground );
    }else
    {
      if( !background->energy_calibration() || !background->energy_calibration()->valid() )
        throw runtime_error( "RelActAutoCostFcn: invalid background spectrum." );
      
      if( (background->live_time() <= 0.0) || (foreground->live_time() <= 0.0) )
        throw runtime_error( "RelActAutoCostFcn: live-time missing from spectrum." );
      
      const double lt_sf = foreground->live_time() / background->live_time();
      const vector<float> &orig_back_counts = *background->gamma_counts();
      const vector<float> &back_energies = *background->energy_calibration()->channel_energies();
      const vector<float> &fore_energies = *foreground->energy_calibration()->channel_energies();
      
      vector<float> background_counts;
      SpecUtils::rebin_by_lower_edge( back_energies, orig_back_counts, fore_energies, background_counts );
      
      assert( background_counts.size() == m_channel_counts.size() );
      for( size_t i = 0; i < m_channel_counts.size(); ++i )
      {
        const double fore_counts = (m_channel_counts[i] < 0.0f) ? 0.0 : m_channel_counts[i];
        const double back_counts = (background_counts[i] < 0.0f) ? 0.0f : background_counts[i];
        const double uncert_2 = fore_counts*fore_counts + lt_sf*lt_sf*back_counts*back_counts;
        const double sub_val = fore_counts - lt_sf*back_counts;
        
        m_channel_counts[i] = static_cast<float>( std::max(sub_val,0.0) );
        m_channel_count_uncerts[i] = static_cast<float>( std::max( 1.0, sqrt(uncert_2) ) );
      }//for( loop over and set channel counts and uncertainties )
      
      m_spectrum = make_shared<SpecUtils::Measurement>( *foreground );
      m_spectrum->set_gamma_counts( make_shared<vector<float>>(m_channel_counts), foreground->live_time(), foreground->real_time() );
    }//if( !background ) / else
    
    
    
    
    
    //Need to initialize m_nuclides
    if( nuclides.empty() )
      throw runtime_error( "RelActAutoCostFcn: no nuclides specified." );
    
    for( const auto &n : nuclides )
    {
      if( !n.nuclide )
        throw runtime_error( "RelActAutoCostFcn: null Nuclide." );
      
      for( const auto &pn : m_nuclides )
      {
        if( n.nuclide == pn.nuclide )
          throw runtime_error( "RelActAutoCostFcn: duplicate nuclide (" + n.nuclide->symbol + ")." );
      }
      
      m_nuclides.emplace_back( n );
    }//for( const auto &n : nuclides )
    
    
    if( !drf || !drf->hasResolutionInfo() )
    {
      std::shared_ptr<DetectorPeakResponse> new_drf;
      if( drf && drf->isValid() )
        new_drf = make_shared<DetectorPeakResponse>( *drf );
      else
        new_drf = make_shared<DetectorPeakResponse>( "FLAT", "FLAT" );
      
      try
      {
        const std::vector<float> drf_coefs{ 0.0f, 0.0f, 0.0f, 0.0f }, uncerts;
        new_drf->fromExpOfLogPowerSeriesAbsEff( drf_coefs, uncerts,
                                               25*PhysicalUnits::cm,
                                               2*PhysicalUnits::cm,
                                               PhysicalUnits::keV,
                                               m_foreground->gamma_energy_min(),
                                               m_foreground->gamma_energy_max() );
        
        vector<shared_ptr<const PeakDef> > peaks = ExperimentalAutomatedPeakSearch::search_for_peaks( foreground, nullptr, {}, false );
          
        auto all_peaks = make_shared<deque<shared_ptr<const PeakDef>>>( begin(peaks), end(peaks) );
        
        new_drf->fitResolution( all_peaks, foreground, DetectorPeakResponse::kGadrasResolutionFcn );
        
        drf = new_drf;
      }catch( std::exception &e )
      {
        cerr << "RelActAutoCostFcn: error fitting FWHM for DRF: " << e.what() << endl;
        
        drf.reset();
      }//try / catch (setup drf FWHM info )
        
    }//if( we dont have fwhm info )
    
    m_drf = drf;
    
    
    //Need to initialize m_energy_ranges
    if( energy_ranges.empty() )
      throw runtime_error( "RelActAutoCostFcn: no energy ranges defined." );
    
    const bool highres = PeakFitUtils::is_high_res( foreground );
    
    for( size_t roi_index = 0; roi_index < energy_ranges.size(); ++roi_index )
    {
      const RelActCalcAuto::RoiRange &roi_range = energy_ranges[roi_index];
      
      if( (roi_range.lower_energy >= roi_range.upper_energy)
         || (roi_range.lower_energy < 0.0) )
        throw runtime_error( "RelActAutoCostFcn: mal-formed energy range" );
      
      // We'll check the ranges dont overlap (but they can touch)
      for( size_t other_roi = roi_index + 1; other_roi < energy_ranges.size(); ++other_roi )
      {
        const RelActCalcAuto::RoiRange &other_roi_range = energy_ranges[other_roi];
        
        if( (roi_range.lower_energy < other_roi_range.upper_energy)
           && (other_roi_range.lower_energy < roi_range.upper_energy) )
        {
          throw runtime_error( "RelActAutoCostFcn: input energy ranges are overlapping ["
                              + std::to_string(roi_range.lower_energy) + ", "
                              + std::to_string(roi_range.upper_energy) + "] and ["
                              + std::to_string(other_roi_range.lower_energy) + ", "
                              + std::to_string(other_roi_range.upper_energy) + "]."
                              );
        }
      }//for( loop over ROIs that come after roi_range )

      if( roi_range.force_full_range )
      {
        assert( !roi_range.allow_expand_for_peak_width );
        if( roi_range.allow_expand_for_peak_width )
          throw runtime_error( "RelActAutoCostFcn: RoiRange::force_full_range and RoiRange::allow_expand_for_peak_width can not both be true." );
        
        m_energy_ranges.emplace_back( foreground, roi_range );
        continue;
      }//if( roi_range.force_full_range )
      
      //We'll try to limit/break-up energy ranges
      const double min_br = numeric_limits<double>::epsilon();  //arbitrary
      const double num_fwhm_roi = 2.5; // arbitrary...
      
      vector<pair<double,double>> gammas_in_range;
      for( const auto &n : m_nuclides )
      {
        for( const auto &g : n.nominal_gammas )
        {
          if( (g.numPerSecond > min_br)
             && (g.energy >= roi_range.lower_energy)
             && (g.energy <= roi_range.upper_energy) )
          {
            double energy_sigma;
            float min_sigma, max_sigma;
            expected_peak_width_limits( g.energy, highres, min_sigma, max_sigma );
            
            if( m_drf && m_drf->hasResolutionInfo() )
            {
              energy_sigma = m_drf->peakResolutionSigma(g.energy);
              
              // A sanity check... maybe we dont want this?
              if( energy_sigma < min_sigma )
                energy_sigma = min_sigma;
              if( energy_sigma > max_sigma )
                energy_sigma = max_sigma;
            }else
            {
              energy_sigma = max_sigma;
            }
            
            double gamma_row_lower = g.energy - num_fwhm_roi*energy_sigma;
            double gamma_row_upper = g.energy + num_fwhm_roi*energy_sigma;
            
            if( !roi_range.allow_expand_for_peak_width )
            {
              gamma_row_lower = std::max( gamma_row_lower, roi_range.lower_energy );
              gamma_row_upper = std::min( gamma_row_upper, roi_range.upper_energy );
            }
            
            gammas_in_range.push_back( {gamma_row_lower,gamma_row_upper} );
          }
        }//for( const auto &g : n.nominal_gammas )
      }//for( const auto &n : m_nuclides )
      
      std::sort( begin(gammas_in_range), end(gammas_in_range) );
      
      for( size_t index = 0; index < gammas_in_range.size();  )
      {
        const size_t start_index = index;
        for( index += 1; index < gammas_in_range.size(); ++index )
        {
          if( gammas_in_range[index - 1].second < gammas_in_range[index].first )
            break;
        }
        
        RelActCalcAuto::RoiRange this_range = roi_range;
        this_range.lower_energy = gammas_in_range[start_index].first;
        this_range.upper_energy = gammas_in_range[index - 1].second;
        
        // TODO: its possible the the channels of the ranges could overlap - right now the overlapping channel will be double counted; should fix.
        
        m_energy_ranges.emplace_back( foreground, this_range );
      }//for( loop over gammas_in_range )
    }//for( const RelActCalcAuto::RoiRange &input : energy_ranges )
    
    // TODO: set m_rel_eff_punishment to the largest peak area divided m_live_time, or something like that
  }//RelActAutoCostFcn constructor.
  
  
  
  size_t number_parameters() const
  {
    size_t num_pars = 0;
    
    // Energy calibration; we will always have these, even if fixed values
    num_pars += 2;
    
    // The FWHM equation
    switch( m_options.fwhm_form )
    {
      case RelActCalcAuto::FwhmForm::Gadras:       num_pars += 3; break;
      case RelActCalcAuto::FwhmForm::Polynomial_2: num_pars += 2; break;
      case RelActCalcAuto::FwhmForm::Polynomial_3: num_pars += 3; break;
      case RelActCalcAuto::FwhmForm::Polynomial_4: num_pars += 4; break;
      case RelActCalcAuto::FwhmForm::Polynomial_5: num_pars += 5; break;
      case RelActCalcAuto::FwhmForm::Polynomial_6: num_pars += 6; break;
    }//switch( m_options.fwhm_form )
    
    // The Relative Eff coefficients
    num_pars += (m_options.rel_eff_eqn_order + 1);
    
    // The Activities; one parameter for activity, one for age (which may be fixed)
    num_pars += 2*m_nuclides.size();
    
    // Floating peaks; one parameter for amplitude, one for FWHM (which will usually be unused)
    num_pars += 2*m_extra_peaks.size();
    
    // Anything else?
    
    return num_pars;
  }//number_parameters()
  
  size_t number_residuals() const
  {
    // Number of gamma channels, plus one to anchor relative eff
    size_t num_resids = 1;
    for( const auto &r : m_energy_ranges )
      num_resids += r.num_channels;
    
    return num_resids;
  }//size_t number_residuals() const
  
  /** Solve the problem, using the Ceres optimizer. */
  static RelActCalcAuto::RelActAutoSolution solve_ceres( RelActCalcAuto::Options options,
                                                        std::vector<RelActCalcAuto::RoiRange> energy_ranges,
                                                        std::vector<RelActCalcAuto::NucInputInfo> nuclides,
                                                        std::vector<RelActCalcAuto::FloatingPeak> extra_peaks,
                                                        RelActCalcAuto::FwhmForm fwhm_form,
                                                        RelActCalc::RelEffEqnForm rel_eff_form,
                                                        size_t rel_eff_order,
                                                        std::shared_ptr<const SpecUtils::Measurement> foreground,
                                                        std::shared_ptr<const SpecUtils::Measurement> background,
                                                        const std::shared_ptr<const DetectorPeakResponse> drf
                                                        )
  {
    const auto start_time = std::chrono::high_resolution_clock::now();
    
    RelActCalcAuto::RelActAutoSolution solution;
    
    DoWorkOnDestruct setFinalTime( [&solution,start_time](){
      const auto end_time = std::chrono::high_resolution_clock::now();
      solution.m_num_microseconds_eval = std::chrono::duration<double, std::micro>(end_time - start_time).count();
    });
    
    solution.m_foreground = foreground;
    solution.m_background = background;
    solution.m_rel_eff_form = rel_eff_form;
    
    //std::vector<double> solution.m_rel_eff_coefficients;
    //std::vector<std::vector<double>> solution.m_rel_eff_covariance;
    
    //std::vector<NuclideRelAct> solution.m_rel_activities;
    //std::vector<std::vector<double>> solution.m_rel_act_covariance;
    
    solution.m_fwhm_form = fwhm_form;
    //std::vector<double> solution.m_fwhm_coefficients;
    //std::vector<std::vector<double>> solution.m_fwhm_covariance;
    
    //solution.m_floating_peaks = extra_peaks;
    
    //std::vector<PeakDef> solution.m_fit_peaks;
    solution.m_input_roi_ranges = energy_ranges;
    solution.m_fit_energy_cal_adjustments = options.fit_energy_cal;
    if( drf && drf->isValid() && drf->hasResolutionInfo() )
      solution.m_drf = drf;
    
    RelActAutoCostFcn *cost_functor = nullptr;
    try
    {
      cost_functor = new RelActAutoCostFcn( options, energy_ranges, nuclides,
                                       extra_peaks, foreground, background, drf );
    }catch( std::exception &e )
    {
      solution.m_status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
      solution.m_error_message = "Error initializing problem: " + string(e.what());
      
      return solution;
    }//try / catch
  
    auto cost_function = new ceres::DynamicNumericDiffCostFunction<RelActAutoCostFcn>( cost_functor );
    cost_function->SetNumResiduals( cost_functor->number_residuals() );
    
    
    solution.m_status = RelActCalcAuto::RelActAutoSolution::Status::FailToSolveProblem;
    solution.m_drf = cost_functor->m_drf;
    solution.m_spectrum = cost_functor->m_spectrum;
    
    const size_t num_pars = cost_functor->number_parameters();
    
    for( size_t i = 0; i < num_pars; ++i )
      cost_function->AddParameterBlock( 1 );
    
    vector<double> parameters( num_pars, 0.0 );
    double *pars = &parameters[0];
    
    vector<double *> parameter_blocks( num_pars );
    for( size_t i = 0; i < num_pars; ++i )
      parameter_blocks[i] = pars + i;
    
    
    ceres::Problem problem;
    
    // TODO: investigate using a LossFunction - probably really need it
    ceres::LossFunction *lossfcn = nullptr;
    problem.AddResidualBlock( cost_function, lossfcn, parameter_blocks );
    
    
    if( options.fit_energy_cal )
    {
      assert( !cost_functor->m_energy_ranges.empty() );
      
      const double lowest_energy = cost_functor->m_energy_ranges.front().lower_energy;
      const double highest_energy = cost_functor->m_energy_ranges.back().upper_energy;
      
      // We'll allow changing the gain by 1% (limit chosen fairly arbitrarily)
      problem.SetParameterLowerBound(pars + 1, 0, -0.01 );
      problem.SetParameterUpperBound(pars + 1, 0, +0.01 );
      
      if( (lowest_energy < 200) && (highest_energy > 600) )
      {
        //We'll allow changing the offset by 5 keV (limit chosen fairly arbitrarily)
        problem.SetParameterLowerBound(pars + 0, 0, -5.0 );
        problem.SetParameterUpperBound(pars + 0, 0, +5.0 );
      }else
      {
        //We'll only fit gain
        problem.SetParameterBlockConstant( pars + 0 );
      }
    }else
    {
      problem.SetParameterBlockConstant( pars + 0 );
      problem.SetParameterBlockConstant( pars + 1 );
    }
    
    /*
    // The FWHM equation
    switch( m_options.fwhm_form )
    {
      case RelActCalcAuto::FwhmForm::Gadras:       num_pars += 3; break;
      case RelActCalcAuto::FwhmForm::Polynomial_2: num_pars += 2; break;
      case RelActCalcAuto::FwhmForm::Polynomial_3: num_pars += 3; break;
      case RelActCalcAuto::FwhmForm::Polynomial_4: num_pars += 4; break;
      case RelActCalcAuto::FwhmForm::Polynomial_5: num_pars += 5; break;
      case RelActCalcAuto::FwhmForm::Polynomial_6: num_pars += 6; break;
    }//switch( m_options.fwhm_form )
    */
  }//void solve_ceres( ceres::Problem &problem )
  
  
  void eval( const std::vector<double> &x, double *residuals ) const
  {
    assert( x.size() == number_parameters() );
    
    blah blah blah - do work here
  }//void eval( const std::vector<double> &x, double *residuals ) const
  
  
  virtual double operator()( const std::vector<double> &x ) const
  {
    vector<double> residuals( number_residuals(), 0.0 );
    try
    {
      eval( x, residuals.data() );
    }catch( std::exception &e )
    {
      cerr << "RelActAutoCostFcn::operator() caught: " << e.what() << endl;
      return std::numeric_limits<double>::max();
    }
    
    double chi2 = 0.0;
    for( const double &d : residuals )
      chi2 += d*d;
    
    return chi2;
  }//operator() - for minuit
  
  
  // For Minuit2
  virtual double Up() const
  {
    return 1.0;
  }

  
  // The return value indicates whether the computation of the
  // residuals and/or jacobians was successful or not.
  bool operator()( double const *const *parameters, double *residuals ) const
  {
    try
    {
      vector<double> pars( number_parameters(), 0.0 );
      
      blah blah blah - do work getting parameters into order
        
      eval( pars, residuals );
    }catch( std::exception &e )
    {
      cerr << "RelActAutoCostFcn::operator() caught: " << e.what() << endl;
      return false;
    }
    
    return true;
  };//bool operator() - for Ceres
};//struct RelActAutoCostFcn

}


namespace RelActCalcAuto
{

Options::Options()
: fit_energy_cal( false ),
  rel_eff_eqn_type( RelActCalc::RelEffEqnForm::LnX ),
  rel_eff_eqn_order( 3 ),
  fwhm_form( FwhmForm::Polynomial_2 )
{
}

RelActAutoSolution::RelActAutoSolution()
: m_status( RelActAutoSolution::Status::NotInitiated ),
  m_error_message( "" ),
  m_foreground{ nullptr },
  m_background{ nullptr },
  m_spectrum{ nullptr },
  m_rel_eff_form( RelActCalc::RelEffEqnForm::LnX ),
  m_rel_eff_coefficients{},
  m_rel_eff_covariance{},
  m_rel_activities{},
  m_rel_act_covariance{},
  m_fwhm_form( FwhmForm::Gadras ),
  m_fwhm_coefficients{},
  m_fwhm_covariance{},
  m_floating_peaks{},
  m_fit_peaks{},
  m_input_roi_ranges{},
  m_energy_cal_adjustments{ 0.0, 0.0 },
  m_drf{ nullptr },
  m_fit_energy_cal_adjustments( false ),
  m_chi2( -1.0 ),
  m_dof( 0 ),
  m_num_function_eval( 0 ),
  m_num_microseconds_eval( 0 )
{
  
}



RelActAutoSolution solve( Options options,
                         std::vector<RoiRange> energy_ranges,
                         std::vector<NucInputInfo> nuclides,
                         std::vector<FloatingPeak> extra_peaks,
                         FwhmForm fwhm_form,
                         RelActCalc::RelEffEqnForm rel_eff_form,
                         size_t rel_eff_order,
                         std::shared_ptr<const SpecUtils::Measurement> foreground,
                         std::shared_ptr<const SpecUtils::Measurement> background
                         )
{
  return RelActAutoCostFcn::solve_ceres(
                     options,
                     energy_ranges,
                     nuclides,
                     extra_peaks,
                     fwhm_form,
                     rel_eff_form,
                     rel_eff_order,
                     foreground,
                     background );
}//RelActAutoSolution


}//namespace RelActCalcAuto
