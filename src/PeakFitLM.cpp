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


#ifdef _MSC_VER
#undef isinf
#undef isnan
#undef isfinite
#undef isnormal
#endif

#include "ceres/ceres.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "SpecUtils/SpecFile.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/PeakFitUtils.h"
#include "SpecUtils/EnergyCalibration.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;

namespace PeakFitLM
{
// TODO: can we make this a "CostFunctor" where operator() is templated to allow better differentiation?
// Takes 1 parameter blocks, of size number_parameters() - e.g., we are putting all parameters together:
// Emits number_residuals() residuals
//
// This is the somewhat equivalent of the MultiPeakFitChi2Fcn class
//
// TODO: currently will uses the entire first/last bin or ROI energy range, no matter how-little the ROI extends into those channels - should probably do some rounding or something.
struct PeakFitDiffCostFunction
{
  PeakFitDiffCostFunction(const std::shared_ptr<const SpecUtils::Measurement> data,
                          const std::vector< std::shared_ptr<const PeakDef> > &starting_peaks,
                          const double roi_lower_energy,
                          const double roi_upper_energy,
                          const PeakContinuum::OffsetType offset_type,
                          const double continuum_ref_energy )
  : m_data( data ),
    m_lower_channel( data->find_gamma_channel( roi_lower_energy ) ),
    m_upper_channel( data->find_gamma_channel( roi_upper_energy ) ),
    m_roi_lower_energy( roi_lower_energy ),
    m_roi_upper_energy( roi_upper_energy ),
    m_num_peaks( starting_peaks.size() ),
    // Make sure m_starting_peaks is sorted - use lambda trick so all member variables are const
    m_starting_peaks( ([starting_peaks]() -> vector<shared_ptr<const PeakDef>>{
      vector<shared_ptr<const PeakDef>> sorted_peaks = starting_peaks;
      std::sort( begin(sorted_peaks), end(sorted_peaks), &PeakDef::lessThanByMeanShrdPtr );
      return sorted_peaks;
    } )() ),
    // The number of parameters depend on how many sigmas are being fit; again, use lambda trick
    m_num_fit_sigmas( ([starting_peaks]() -> size_t{
      size_t n_fit_sigma = 0;
      for( const auto &p : starting_peaks )
        n_fit_sigma += p->fitFor(PeakDef::Sigma);
      return n_fit_sigma;
    })() ),
    m_offset_type( offset_type ),
    m_ref_energy( continuum_ref_energy ),
    m_ncalls( 0 )
  {
    assert( m_data );
    assert( m_data->gamma_counts() && !m_data->gamma_counts()->empty() );
    assert( m_upper_channel >= m_lower_channel );
    assert( m_upper_channel < m_data->num_gamma_channels() );
    assert( m_roi_lower_energy < m_roi_upper_energy );
    
    assert( m_num_peaks > 0 );
    assert( m_num_peaks < (m_upper_channel - m_lower_channel) );
    
    if( !m_data || !m_data->gamma_counts() || m_data->gamma_counts()->empty() )
      throw runtime_error( "PeakFitDiffCostFunction: !m_data or !gamma_counts()" );
    
    if( m_roi_lower_energy >= m_roi_upper_energy )
      throw runtime_error( "PeakFitDiffCostFunction: m_roi_lower_energy >= m_roi_upper_energy" );
    
    if( m_upper_channel < m_lower_channel )
      throw runtime_error( "PeakFitDiffCostFunction: m_upper_channel < m_lower_channel" );
      
    const std::vector<float> &counts = *m_data->gamma_counts();
    if( m_upper_channel >= counts.size() )
      throw runtime_error( "PeakFitDiffCostFunction: m_upper_channel >= counts.size()" );
    
    if( !m_num_peaks )
      throw runtime_error( "PeakFitDiffCostFunction: no peaks to fit" );
    
    for( size_t i = 1; i < starting_peaks.size(); ++i )
    {
      if( starting_peaks[i]->continuum() != starting_peaks[0]->continuum() )
        throw runtime_error( "PeakFitDiffCostFunction: input peaks must all share a continuum" );
    }
    
    auto continuum = m_starting_peaks[0]->continuum();
    // TODO: get rid of m_roi_lower_energy, m_roi_upper_energy, m_offset_type, and m_ref_energy - since they are all duplicated
    assert( continuum->lowerEnergy() == m_roi_lower_energy );
    assert( continuum->upperEnergy() == m_roi_upper_energy );
    assert( continuum->type() == m_offset_type );
    assert( continuum->referenceEnergy() == m_ref_energy );
    
    
    
    bool valid_offset_type = false;
    switch( offset_type )
    {
      case PeakContinuum::NoOffset:     case PeakContinuum::Constant:
      case PeakContinuum::Linear:       case PeakContinuum::Quadratic:
      case PeakContinuum::Cubic:        case PeakContinuum::FlatStep:
      case PeakContinuum::LinearStep:   case PeakContinuum::BiLinearStep:
      case PeakContinuum::External:
        valid_offset_type = true;
        break;
    }//switch( offset_type )
    
    assert( valid_offset_type );
    
    if( !valid_offset_type )
      throw runtime_error( "PeakFitDiffCostFunction: !valid_offset_type" );
    
    // It wouldnt be that hard to support external continuum, but we'll leave that for later, at
    //  the moment.
    assert( m_offset_type != PeakContinuum::OffsetType::External );
    if( m_offset_type == PeakContinuum::OffsetType::External )
      throw runtime_error( "PeakFitDiffCostFunction: PeakContinuum::OffsetType::External not supported" );
    
    if( dof() <= 0.0 )
      throw runtime_error( "PeakFitDiffCostFunction: not enough data channels to fit all the parameters" );
  }//PeakFitDiffCostFunction constructor
  
  double dof() const
  {
    const double num_channels = 1.0 + m_upper_channel - m_lower_channel;
    const double num_continuum_pars = PeakContinuum::num_parameters( m_offset_type );
    const double num_sigmas_fit = number_sigma_parameters();
    
    // fixed sigmas are handled by just not having them included in the parameters, so we wont
    //  count fixed sigmas in this next loop.
    size_t num_fixed_pars = 0;
    for( const auto &p : m_starting_peaks )
    {
      num_fixed_pars += !p->fitFor(PeakDef::GaussAmplitude);
      num_fixed_pars += !p->fitFor(PeakDef::Mean);
    }
    
    // Get the number of fixed continuum parameters.
    for( const bool do_fit_par : m_starting_peaks[0]->continuum()->fitForParameter() )
      num_fixed_pars += !do_fit_par;
    
    // TODO: MultiPeakFitChi2Fcn::dof() adds 1 to the degrees of freedom - not quite sure why that is, at the moment, or if we should do this here.
    return num_channels - 2.0*m_num_peaks + num_fixed_pars - num_sigmas_fit - num_continuum_pars;
  }//double dof() const
  
  size_t number_parameters() const
  {
    const size_t num_continuum_pars = PeakContinuum::num_parameters( m_offset_type );
    return num_continuum_pars + number_sigma_parameters() + 2*m_num_peaks;
  }
  
  size_t number_sigma_parameters() const
  {
    return std::min( m_num_fit_sigmas, size_t(2) );
  }
  
  size_t number_residuals() const
  {
    return m_num_peaks + m_upper_channel - m_lower_channel;
  }//size_t number_residuals() const
  
  
  vector<PeakDef> parametersToPeaks( const double * const params, const double * const uncertainties, double *residuals ) const
  {
    const size_t num_continuum_pars = PeakContinuum::num_parameters( m_offset_type );
    const size_t last_data_residual = m_upper_channel - m_lower_channel;
    
    auto continuum = std::make_shared<PeakContinuum>( *m_starting_peaks[0]->continuum() );
    assert( continuum->lowerEnergy() == m_roi_lower_energy );
    assert( continuum->upperEnergy() == m_roi_upper_energy );
    assert( continuum->type() == m_offset_type );
    assert( continuum->referenceEnergy() == m_ref_energy );
    
    continuum->setRange( m_roi_lower_energy, m_roi_upper_energy );
    continuum->setType( m_offset_type );
    continuum->setParameters( m_ref_energy, params, uncertainties );
    
    std::unique_ptr<vector<double>> local_residuals;
    if( !residuals )
    {
      local_residuals.reset( new vector<double>( number_residuals(), 0.0 ) );
      residuals = local_residuals->data();
    }
      
    const size_t num_sigmas_fit = number_sigma_parameters();
    
    vector<PeakDef> peaks( m_num_peaks );
    
    for( size_t i = 0; i < m_num_peaks; ++i )
    {
      const size_t start_par = num_continuum_pars + num_sigmas_fit + 2*i;
      assert( (start_par + 1) < number_parameters() );
      
      const double mean = params[start_par + 0];
      const double amp = params[start_par + 1];
      
      double sigma, sigma_uncert;
      
      if( !m_starting_peaks[i]->fitFor(PeakDef::Sigma) || (num_sigmas_fit == 0) )
      {
        sigma = m_starting_peaks[i]->sigma();
        sigma_uncert = m_starting_peaks[i]->sigmaUncert();
      }else
      {
        sigma = params[num_continuum_pars];
        sigma_uncert = uncertainties ? uncertainties[num_continuum_pars] : 0.0;
        
        if( num_sigmas_fit > 1 )
        {
          assert( (num_continuum_pars + 1) < number_parameters() );
          
          const double frac = (mean - m_roi_lower_energy) / (m_roi_upper_energy - m_roi_lower_energy);
          sigma = params[num_continuum_pars] + (frac * params[num_continuum_pars + 1]);
          
          if( uncertainties )
          {
            //TODO: I'm not entirely conviced uncert is correct (note: assuming 100% correlation)
            //TODO: use the covariance between these two parameters to better assign sigma uncert
            sigma_uncert += frac * uncertainties[num_continuum_pars + 1];
          }
        }//if( num_sigmas_fit > 1 )
      }//if( num_sigmas_fit == 0 ) / else
      
      peaks[i].setMean( mean );
      peaks[i].setSigma( sigma );
      peaks[i].setAmplitude( amp );
      
      if( uncertainties )
      {
        double mean_uncert = uncertainties[start_par + 0];
        double amp_uncert = uncertainties[start_par + 1];
        
        if( !m_starting_peaks[i]->fitFor(PeakDef::Mean) )
          mean_uncert = m_starting_peaks[i]->meanUncert();
        
        if( !m_starting_peaks[i]->fitFor(PeakDef::GaussAmplitude) )
          amp_uncert = m_starting_peaks[i]->amplitudeUncert();
        
        if( mean_uncert > 0.0 )
          peaks[i].setMeanUncert( mean_uncert );
          
        if( amp_uncert > 0.0 )
          peaks[i].setAmplitudeUncert( amp_uncert );
        
        if( sigma_uncert > 0.0 )
          peaks[i].setSigmaUncert( sigma_uncert );
      }//if( uncertainties )
      
      peaks[i].setContinuum( continuum );
    }//for( size_t i = 0; i < m_num_peaks; ++i )
    
    // Sort the peaks to make calculating the punishment for peaks being to close easier
    std::sort( begin(peaks), end(peaks), &PeakDef::lessThanByMean );
    
    double chi2 = 0.0;
    
    for( size_t channel = m_lower_channel; channel <= m_upper_channel; ++channel )
    {
      const float ndata = m_data->gamma_channel_content(channel);
      const float channel_lower_energy = m_data->gamma_channel_lower( channel );
      const float channel_upper_energy = m_data->gamma_channel_upper( channel );
      
      const double ncontinuum = peaks[0].continuum()->offset_integral(channel_lower_energy, channel_upper_energy, m_data);
      
      double peak_gaus_contribution = 0.0;
      for( const PeakDef &peak : peaks )
        peak_gaus_contribution += peak.gauss_integral( channel_lower_energy, channel_upper_energy );
      
      const size_t residual_num = channel - m_lower_channel;
      assert( residual_num < number_residuals() );
      
      if( ndata > 0.000001f )  // The cutoff is arbitrary, but needs to be greater than zero
        residuals[residual_num] = (ndata - ncontinuum - peak_gaus_contribution) / sqrt(ndata);
      else
        residuals[residual_num] = peak_gaus_contribution + ncontinuum; // TODO: This is ad-hoc, and just following PeakFitChi2Fcn implementation for the moment
      
      chi2 += residuals[residual_num] * residuals[residual_num];
    }//for( size_t channel = m_lower_channel; channel <= m_upper_channel; ++channel )
    
    const double chi2_dof = chi2 / dof();
    for( size_t i = 0; i < m_num_peaks; ++i )
      peaks[i].set_coefficient( chi2_dof, PeakDef::Chi2DOF );
    
    // Punish if the peaks are too close.
    for( size_t i = 1; i < peaks.size(); ++i )
    {
      const double sigma_i = peaks[i].gausPeak() ? peaks[i].sigma() : 0.25*peaks[i].roiWidth();
      const double sigma_j = peaks[i-1].gausPeak() ? peaks[i-1].sigma() : 0.25*peaks[i-1].roiWidth();
      const double sigma = 0.5*(sigma_j + sigma_i);
      
      const double dist = fabs(peaks[i-1].mean() - peaks[i].mean());
      double reldist = dist / sigma;
      if( reldist < 0.01 || IsInf(reldist) || IsNan(reldist) )
        reldist = 0.01;
      
      const double punishment_factor = 2 * (m_upper_channel - m_lower_channel) / m_num_peaks;
      
      assert( (last_data_residual + i) < number_residuals() );
      
      if( reldist < 1.0 )
      {
        residuals[last_data_residual + i] = (punishment_factor / reldist);
      }else if( reldist < 1.5 )
      {
        // At reldist==1, this will be equal to the above case,
        //  and by reldist==1.5: exp( -pow(5*0.5,3.0) )== 1.64E-7
        residuals[last_data_residual + i] = punishment_factor * std::exp( -std::pow(5*(reldist - 1),3.0) );
      }else
      {
        residuals[last_data_residual + i] = 0.0;
      }
    }//for( size_t i = 1; i < peaks.size(); ++i )
    
    // TODO: PeakFitChi2Fcn also punishes for peak being statistically insignificant
    
    
    // Now inherit all the assigned nuclides, colors, whatever from the starting peaks.
    //  However, ordering of peaks could swap, etc, so we need to match the original
    //  peaks up to the closest current peak.
    assert( m_starting_peaks.size() == peaks.size() );
    map<size_t,size_t> old_to_new_indexs;
    
    // If the mean of a peak is fixed, then we "know" it wont swap positions
    for( size_t index = 0; index < m_starting_peaks.size(); ++index )
    {
      if( !m_starting_peaks[index]->fitFor(PeakDef::Mean) )
        old_to_new_indexs[index] = index;
    }
    
    for( size_t new_index = 0; new_index < peaks.size(); ++new_index )
    {
      const double new_mean = peaks[new_index].mean();
      
      // Make sure we havent already assigned this new peak to an old peak (because mean was fixed)
      bool already_assigned = false;
      for( const auto &t : old_to_new_indexs )
        already_assigned = (already_assigned || (t.second == new_index));
      
      if( already_assigned )
        continue;
      
      int closest_orig = -1;
      for( size_t orig_index = 0; orig_index < m_starting_peaks.size(); ++orig_index )
      {
        // Check if we've already assigned this index
        if( old_to_new_indexs.count(orig_index) )
          continue;
        
        if( closest_orig < 0 )
        {
          closest_orig = static_cast<int>( orig_index );
        }else
        {
          const double prev_mean = m_starting_peaks[closest_orig]->mean();
          const double orig_mean = m_starting_peaks[orig_index]->mean();
          if( fabs(new_mean - orig_mean) < fabs(new_mean - prev_mean) )
            closest_orig = static_cast<int>( orig_index );
        }
      }//for( loop over original peaks )
      
      assert( closest_orig >= 0 );
      
      old_to_new_indexs[static_cast<size_t>(closest_orig)] = new_index;
    }//for( loop over new peaks )
    
    assert( old_to_new_indexs.size() == peaks.size() );
    
    for( const auto &old_to_new : old_to_new_indexs )
    {
      assert( old_to_new.first < peaks.size() );
      assert( old_to_new.second < peaks.size() );
      const bool inheritNonFitForValues = true; //I dont think this really matters
      peaks[old_to_new.second].inheritUserSelectedOptions( *m_starting_peaks[old_to_new.first],  inheritNonFitForValues );
    }//for( const auto &old_to_new : old_to_new_indexs )
    
    return peaks;
  }//parametersToPeaks(...)
  
  
  bool operator()( double const *const *parameters, double *residuals ) const
  {
    m_ncalls += 1;
    
    const size_t npar = number_parameters();
    vector<double> pars( npar );
    for( size_t i = 0; i < npar; ++i )
      pars[i] = parameters[i][0];
    
    parametersToPeaks( pars.data(), nullptr, residuals );
    
    return true;
  }//operator()
  
  
  
public:
  const std::shared_ptr<const SpecUtils::Measurement> m_data;
  const size_t m_lower_channel;
  const size_t m_upper_channel;
  const double m_roi_lower_energy;
  const double m_roi_upper_energy;
  
  const size_t m_num_peaks;
  const std::vector< std::shared_ptr<const PeakDef> > m_starting_peaks;
  const size_t m_num_fit_sigmas;
  const PeakContinuum::OffsetType m_offset_type;
  const double m_ref_energy;
  
  mutable std::atomic<unsigned int> m_ncalls;
};//struct PeakFitDiffCostFunction



void fit_peak_for_user_click_LM( PeakShrdVec &results,
                             double &chi2Dof,
                             const std::shared_ptr<const SpecUtils::Measurement> &dataH,
                             PeakShrdVec coFitPeaks,
                             const double mean0, const double sigma0,
                             const double area0,
                             const float roiLowerEnergy,
                             const float roiUpperEnergy )
{
  /** For this first go, we will have Ceres solve everything very stupidly.
   
   This is only tested to the "seems to work, for simple cases" level.
   
   We will:
     - Constrains the widths of the peaks to have a linearly related width
   
   In the future:
     - we should do like LinearProblemSubSolveChi2Fcn and just use matrices to solve for peak amplitudes and continuum parameters - see Options below for how I think this should be setup
      - Should separate the linear and non-linear parameter blocks now at least
     - Add punishment for stat insignificant peaks
     - See if we can template cost function to take advantage of the auto diff
       - It looks like to use ceres::AutoDiffCostFunction, we need to know how many parameters we will use at compile-time (so we could just have special cases for 1 and 2 peaks, and just use DynamicNumericDiffCostFunction otherwise); we also need to inline all the peak and continuum calculations to be in this file so we can evaluate ceres::Jets (we could also potentially template these function)
     - Need to deal with errors during solving, and during Covariance computation
     - Deal with setting number of threads reasonably.
     - Try out using a ceres::LossFunction
     - ...
   
   */
  
  typedef std::shared_ptr<const PeakDef> PeakDefShrdPtr;
  
  // Lets make sure all coFitPeaks share a continuum.
  for( size_t i = 1; i < coFitPeaks.size(); ++i )
  {
    if( coFitPeaks[i]->continuum() != coFitPeaks[0]->continuum() )
      throw runtime_error( "fit_peak_for_user_click_LM: input peaks all must share a single continuum" );
  }//for( size_t i = 1; i < coFitPeaks.size(); ++i )
  
  chi2Dof = DBL_MAX;
  results.clear();
  
  try
  {
    // TODO: check that repeaded calls to this function wont cause adding/removing channels due to find_gamma_channel(...) rounding type things - e.g., do we need to subtract off a tiny bit from roiUpperEnergy
    const size_t lower_channel = dataH->find_gamma_channel(roiLowerEnergy);
    const size_t upper_channel = dataH->find_gamma_channel(roiUpperEnergy);
    const double roi_lower = dataH->gamma_channel_lower( lower_channel );
    const double roi_upper = dataH->gamma_channel_upper( upper_channel );
    
    assert( roi_lower <= roiLowerEnergy );
    assert( roi_upper >= roiUpperEnergy );
    
    const size_t nchannels = dataH->num_gamma_channels();
    const size_t midbin = dataH->find_gamma_channel( mean0 );
    const float binwidth = dataH->gamma_channel_width( midbin );
    const bool highres = PeakFitUtils::is_high_res( dataH );
    const size_t num_fit_peaks = coFitPeaks.size() + 1;
    
    double reference_energy = roiLowerEnergy;
    
    //The below should probably go off the number of bins in the ROI
    PeakContinuum::OffsetType offset_type;
    if( highres )
      offset_type = (num_fit_peaks < 3) ? PeakContinuum::Linear : PeakContinuum::Quadratic;
    else
      offset_type = (num_fit_peaks < 2) ? PeakContinuum::Linear : PeakContinuum::Quadratic;
    
    //for( size_t i = 0; i < coFitPeaks.size(); ++i )
    if( coFitPeaks.size() )
    {
      //const auto continuum = coFitPeaks[i]->continuum();
      const auto continuum = coFitPeaks[0]->continuum();
      offset_type = std::max( offset_type, continuum->type() );
    }
    
    const size_t num_continuum_pars = PeakContinuum::num_parameters( offset_type );
    
    bool any_continuum_par_fixed = false;
    // TODO: we can get rid of continuum_parameter_fit... as this will be contained in `initial_continuum` as well.
    vector<bool> continuum_parameter_fit( num_continuum_pars, true );
    //for( size_t i = 0; i < coFitPeaks.size(); ++i )
    if( coFitPeaks.size() )
    {
      //const auto continuum = coFitPeaks[i]->continuum();
      const auto continuum = coFitPeaks[0]->continuum();
      const std::vector<bool> parFitFors = continuum->fitForParameter();
      assert( parFitFors.size() <= num_continuum_pars );
      
      for( size_t fit_for_index = 0; fit_for_index < parFitFors.size(); ++fit_for_index )
      {
        if( !parFitFors[fit_for_index] )
        {
          reference_energy = continuum->referenceEnergy();
          any_continuum_par_fixed = true;
          continuum_parameter_fit[fit_for_index] = false;
        }
      }//
    }//for( size_t i = 0; i < coFitPeaks.size(); ++i )
    
    
    std::shared_ptr<PeakDef> candidatepeak = std::make_shared<PeakDef>(mean0, sigma0, area0);
    
    // Make sure all initial peaks share a continuum
    //  TODO: this common continuum shares a good amount of information we also pass to PeakFitDiffCostFunction ... should consolidate
    shared_ptr<PeakContinuum> initial_continuum;
    if( coFitPeaks.size() )
    {
      initial_continuum = make_shared<PeakContinuum>( *coFitPeaks[0]->continuum() );
      
      candidatepeak->setContinuum( initial_continuum );
      for( size_t i = 0; i < coFitPeaks.size(); ++i )
      {
        auto newpeak = make_shared<PeakDef>( *coFitPeaks[i] );
        newpeak->setContinuum( initial_continuum );
        coFitPeaks[i] = newpeak;
      }
    }else
    {
      initial_continuum = candidatepeak->continuum();
    }//if( coFitPeaks.size() )
    
    assert( initial_continuum );
    
    initial_continuum->setRange( roiLowerEnergy, roiUpperEnergy );
    if( !any_continuum_par_fixed )
      initial_continuum->calc_linear_continuum_eqn( dataH, reference_energy, roiLowerEnergy, roiUpperEnergy, 2, 2 );
    initial_continuum->setType( offset_type );
    
    coFitPeaks.push_back( candidatepeak );
    std::sort( coFitPeaks.begin(), coFitPeaks.end(), &PeakDef::lessThanByMeanShrdPtr );

    
    assert( roiLowerEnergy < roiUpperEnergy );
    if( roiLowerEnergy >= roiUpperEnergy )
      throw runtime_error( "Invalid energy range (" + std::to_string(roiLowerEnergy) + ", " + std::to_string(roiUpperEnergy) + ")" );
      
    const double range = (roi_upper - roi_lower);
      
    double minsigma = binwidth;
    double maxsigma = 0.5*range;
      
    ceres::Problem problem;
    
    auto cost_functor = new PeakFitDiffCostFunction( dataH, coFitPeaks, roiLowerEnergy, roiUpperEnergy, offset_type, reference_energy );
    
    auto cost_function = new ceres::DynamicNumericDiffCostFunction<PeakFitDiffCostFunction>( cost_functor );
      
    const size_t num_fit_pars = cost_functor->number_parameters();
    const size_t num_sigmas_fit = cost_functor->number_sigma_parameters();
    
    for( size_t i = 0; i < num_fit_pars; ++i )
      cost_function->AddParameterBlock( 1 );
    cost_function->SetNumResiduals( static_cast<int>(cost_functor->number_residuals()) );
    
    vector<double> parameters( num_fit_pars, 0.0 );
    double *pars = &parameters[0];
    
    vector<double*> parameter_blocks( num_fit_pars );
    for( size_t i = 0; i < num_fit_pars; ++i )
      parameter_blocks[i] = pars + i;
      
    // TODO: investigate using a LossFunction; from a brief investigation it maybe really affects the amplitude uncertainty if not chosen carefully
    ceres::LossFunction *lossfcn = nullptr;
    //ceres::LossFunction *lossfcn = new ceres::CauchyLoss(0.5) or new HuberLoss(...), or ...
    
    problem.AddResidualBlock( cost_function, lossfcn, parameter_blocks );
    //problem.SetParameterBlockConstant( nullptr );
      
    assert( num_continuum_pars == initial_continuum->parameters().size() );
    
    for( size_t i = 0; i < num_continuum_pars; ++i )
    {
      parameters[i] = initial_continuum->parameters()[i];
      if( !continuum_parameter_fit[i] )
      {
        assert( !initial_continuum->fitForParameter()[i] );
        problem.SetParameterBlockConstant( pars + i );
      }
    }
    
    
    for( size_t i = 0; i < coFitPeaks.size(); ++i )
    {
      const double mean = coFitPeaks[i]->mean();
      const double sigma = coFitPeaks[i]->sigma();
      const double amp = coFitPeaks[i]->amplitude();
      
      const size_t start_par = num_continuum_pars + num_sigmas_fit + 2*i;
      parameters[start_par + 0] = mean;
      parameters[start_par + 1] = amp;
      
      if( !coFitPeaks[i]->fitFor(PeakDef::Mean) )
      {
        problem.SetParameterBlockConstant( pars + start_par + 0 );
      }else
      {
        if( coFitPeaks[i] == candidatepeak )
        {
          problem.SetParameterLowerBound( pars + start_par + 0, 0, mean - sigma );
          problem.SetParameterUpperBound( pars + start_par + 0, 0, mean + sigma );
        }else
        {
          problem.SetParameterLowerBound( pars + start_par + 0, 0, mean - 0.25*sigma );
          problem.SetParameterUpperBound( pars + start_par + 0, 0, mean + 0.25*sigma );
        }
      }//
      
      // Either set the amplitude as fixed, or the minimum amplitude to zero
      if( !coFitPeaks[i]->fitFor(PeakDef::GaussAmplitude) )
        problem.SetParameterBlockConstant( pars + start_par + 1 );
      else
        problem.SetParameterLowerBound( pars + start_par + 1, 0, 0.0 );
      
      if( !coFitPeaks[i]->fitFor(PeakDef::Sigma) )
      {
        minsigma = std::min( minsigma, coFitPeaks[i]->sigma() );
        maxsigma = std::max( maxsigma, coFitPeaks[i]->sigma() );
      }else
      {
        float lowersigma, uppersigma;
        expected_peak_width_limits( mean, highres, dataH, lowersigma, uppersigma );
        if( !i )
          minsigma = lowersigma;
        if( i == (coFitPeaks.size()-1) )
          maxsigma = 1.33*uppersigma;
      }
    }//for( size_t i = 0; i < coFitPeaks.size(); ++i )
    
    if( num_sigmas_fit > 0 )
    {
      problem.SetParameterLowerBound( pars + num_continuum_pars, 0, 0.75*minsigma );
      problem.SetParameterUpperBound( pars + num_continuum_pars, 0, 1.33*maxsigma );
      if( num_sigmas_fit > 1 )
      {
        problem.SetParameterLowerBound( pars + num_continuum_pars + 1, 0, -0.15 );
        problem.SetParameterUpperBound( pars + num_continuum_pars + 1, 0, +0.15 );
      }
    }//if( num_sigmas_fit > 0 )
    
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = false; //true;
    options.num_threads = 4;
    
    // TODO: separate the peak amplitude and continuum parameters into a sub-problem to allow faster/better iteration - not really sure how to do totdally
    //options.use_inner_iterations = true;
    //options.inner_iteration_ordering = make_shared<ceres::ParameterBlockOrdering>();
    //options.inner_iteration_ordering->AddElementToGroup
    //options.inner_iteration_tolerance = ...;
    
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    //std::cout << summary.BriefReport() << "\n";
    std::cout << summary.FullReport() << "\n";
    cout << "Took " << cost_functor->m_ncalls.load() << " calls to solve." << endl;
    cost_functor->m_ncalls = 0;
    
    switch( summary.termination_type )
    {
      case ceres::CONVERGENCE:
      case ceres::USER_SUCCESS:
        break;
        
      case ceres::NO_CONVERGENCE:
      case ceres::FAILURE:
      case ceres::USER_FAILURE:
        throw runtime_error( "The L-M ceres::Solver solving failed." );
        break;
    }//switch( summary.termination_type )
    
    
    //ceres::Solve(options, &problem, nullptr );
    
    
    ceres::Covariance::Options cov_options;
    cov_options.algorithm_type = ceres::CovarianceAlgorithmType::SPARSE_QR; //
    //cov_options.num_threads = 1;
    
    vector<double> uncertainties( num_fit_pars, 0.0 );
    double *uncertainties_ptr = uncertainties.data();
    
    
    ceres::Covariance covariance(cov_options);
    vector<pair<const double*, const double*> > covariance_blocks;
    for( size_t i = 0; i < num_fit_pars; ++i )
      covariance_blocks.push_back( make_pair( pars + i, pars + i ) );
    
    if( !covariance.Compute(covariance_blocks, &problem) )
    {
      cerr << "Failed to compute covariance!" << endl;
      uncertainties_ptr = nullptr;
    }else
    {
      for( size_t i = 0; i < num_fit_pars; ++i )
      {
        double variance = 0.0;
        covariance.GetCovarianceBlock( pars + i, pars + i, &variance );
        
        // TODO: the uncertainties are a bit larger than with Minuit method - should maybe check out.
        if( variance > 0.0 )
          uncertainties[i] = sqrt( variance );
      }
      
      /*
       // This is how we would do it if we had all parameters in a single block
      vector<double> covariance_xx(num_fit_pars * num_fit_pars, 0.0 );
      covariance.GetCovarianceBlock( pars, pars, covariance_xx.data() );
      
      for( size_t i = 0; i < num_fit_pars; ++i )
      {
        assert( (i*(num_fit_pars + 1)) < (num_fit_pars*num_fit_pars) );
        
        const double var = covariance_xx[i*(num_fit_pars + 1)];
        //cout << "Par " << i << " = " << parameters[i] << " +- " << sqrt(var) << endl;
        if( var > 0.0 )
          uncertainties[i] = sqrt( var );
      }//for( size_t i = 0; i < num_fit_pars; ++i )
       */
    }//if( we failed to computer covariance ) / else
    
    
    vector<double> residuals( cost_functor->number_residuals(), 0.0 );
    auto final_peaks = cost_functor->parametersToPeaks( &parameters[0], uncertainties_ptr, residuals.data() );
    
    results.resize( final_peaks.size() );
    for( size_t i = 0; i < final_peaks.size(); ++i )
    {
      chi2Dof = final_peaks[i].chi2dof();
      results[i] = make_shared<PeakDef>( final_peaks[i] );
    }
    
    // One example (two peaks) with 46 iterations had 2529 to solve, and 17 calls to get covariance.
    cout << "Took " << cost_functor->m_ncalls.load() << " calls to get covariance." << endl;
  }catch( std::exception &e )
  {
    cout << "fit_peak_for_user_click_LM caught: " << e.what() << endl;
    //assert( 0 );
  }//try / catch
}//void fit_peak_for_user_click_LM(...)


}//namespace PeakFitLM
